#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "pdecol/pdecol_solver.h"
#include "vmcw_pars.h"
#include "vmcw_mesh.h"
#include "pnm_writer.h"

// I use many global variables here. Maybe it would be better
// to pack everything inside a class

/*****************************/
// Solver parameters.

/* Number of equations. Can not be changed. */
const int npde = 7;

/* How many derivatives to calculate. Can not be changed. */
const int nder = 3;

/* Number of points. It can be changed at any time, but real change
happenes when the solver is (re)started. In the code it can be used only
to initialize the solver. In other places solver->get_npts() should be
used (it corresponds to the actual number of points in the running solver). */
int npts=257;

/* Solver accuracy. Can be changed during calculation. There is in idea that
if it is a power of 2 solver runs faster */
double acc = pow(2,-20);

/* Min time step (recommended 1e-10). It can be changes at any time,
but real change happenes when the solver is (re)started. */
double mindt = 1e-10;

/* Current time [s]. Changes only before doing solver->step(). */
double tcurr = 0;

/* Time step [s]. Can be changed during calculation. */
double tstep = 5e-3;

/* End of current sweep/wait command [s]. Calculation is done until this time
without reading new commands. Should be equal to "time" in the begining. */
double tend = tcurr;

/*****************************/
// NMR frequency, magnetic fields.

/* Gyromagnetic ratio. */
const double gyro = 20378.0;

/* NMR frequency [Hz]. The calculation is done at fixed frequency. Magnetic
fiels should be changed to observe the resonance. Default value 1MHz
corresponds to default magnetic field. */
double f0 = 1e6;

// Magnetic field.
// Total field is: H = 2pi*f0/gyro + H0 + HT*t + HG*x + HQ*x^2. */
//  H0 - uniform field measured from 2pi*f0/gyro [G],
//  HT - sweep rate, dH/dt [G/s],
//  HG - gradient, dH/dz [G/cm],
//  HQ - quadratic term, d^2H/dz^2 [G/cm^2].
double H0=0.0,  HT = 0.0,  HG = 0.1,  HQ = 0.0;

// Radio-frequncy field.
// Total field is: HR = (HR0 + HRT*t) * (1 + x/xRG + (x/xRQ)^2).
// Gradient and quadratic terms change proportianally with the field.
//  HR0 - Initial value [G],
//  HRT - sweep rate, dHr/dt [G/s],
//  HRGP - gradient profile [1/cm],
//  HRQP - quadratic profile [1/cm^2].
double HR0=1e-3, HRT=0.0, HRGP=0.0, HRQP=0.0;

// Type of initial conditions: 0 - plain, 1 - n-soliton, 2 - theta-soliton.
int icond_type = 0;

// Type of boundary conditions: 1 - open cell, 2 - no spin currents.
int bcond_type = 2;

/*****************************/
// He3 properties

// t1 relaxation time, initial value [s] and sweep rate [s/s]
double T10=1.0, T1T=0.0;

// Leggett-Takagi relaxation time tau_f, initial value [s] and sweep rate [s/s]
double TF0=1e-5, TFT=0.0;

// Spin-wave velocity, c_par, initial value [cm/s] and sweep rate [(cm/s)/s]
double CP0=500.0, CPT=0.0;

// Leggett frequency Omega_B, initial value [Hz] and sweep rate [Hz/s]
double LF0=1e5, LFT=0.0;

// Spin diffusion, initial value [cm^2/s] and sweep rate [(cm^2/s)/s]
double DF0=0.1, DFT=0.0;

/*****************************/
// files

// prefix for data files, command file name without extension
std::string pref;

// mesh, prof file counters
int cnt_mesh=0, cnt_prof=0;

/*****************************/

struct pars_t pars; // parameter structure

/* PDECOL solver pointer */
pdecol_solver *solver = NULL;

/* Container for PNM writers. Key is the file name. */
pnm_writers_t pnm_writers;

/******************************************************************/
// Set parameters for Leggett equations using main parameter structure.
// This function is called from F and BNDRY
extern "C" {
  void set_bulk_pars_(double *t, double *x,
                 double *Wr, double *Wz, double *W0,
                 double *WB, double *Cpar, double *dCpar,
                 double *Diff, double *Tf, double *T1){

    *W0 = 2*M_PI*f0;
    *Wz = *W0 + gyro*(H0 + HG*(*x) + HQ*(*x)*(*x) + HT*(*t));
    *Wr = gyro*(HR0 + HRT*(*t)) * (1.0 + (*x)*HRGP + (*x)*(*x)*HRQP);

    *WB = ( LF0 + LFT*(*t))*2*M_PI;
    *Cpar = CP0 + CPT*(*t); *dCpar = 0;
    *Diff = DF0 + DFT*(*t);
    *Tf   = TF0 + TFT*(*t);
    *T1   = T10 + T1T*(*t);

    // spatial modulation
    if (pars.AER){
      *Cpar *= 1.0 - 0.5*aer_step(&pars, *x,0);
      *dCpar =(*Cpar) * 0.5*aer_step(&pars, *x,1);
      *Diff *= 1.0 - 0.835 * aer_step(&pars, *x,0);
      *Tf   *= 1.0 - 0.5 * aer_step(&pars, *x,0);
    }
  }

  void set_bndry_pars_(double *t, double *x, double *W0,
                 double *Cpar, double *Diff, int *IBN){

    *W0 = 2*M_PI*f0;
    *Cpar = CP0 + CPT*(*t);
    *Diff = DF0 + DFT*(*t);

    // spatial modulation
    if (pars.AER){
      *Cpar *= 1.0 - 0.5*aer_step(&pars, *x,0);
      *Diff *= 1.0 - 0.835 * aer_step(&pars, *x,0);
    }
    // type of boundary condition
    *IBN = bcond_type;
  }
}

/******************************************************************/
/// UINIT function called from pdecol
extern "C" {
  void uinit_(double *x, double u[], int *n){
    // default
    u[0] = 0.0; // Mx
    u[1] = 0.0; // My
    u[2] = 1.0; // Mz
    u[3] = 0.0; // nx
    u[4] = 0.0; // ny
    u[5] = 1.0; // nz
    u[6] = acos(-0.25); // theta

    switch (icond_type) {

    case 1: // n-soliton at z=0
      if (*x<0.0)  u[5]=-1.0;
      if (*x==0.0) u[4]=1.0, u[6]=0.0;
      break;

    case 2: // theta-soliton at z=0
      if (*x<0.0)  u[6] = -acos(-0.25);
      if (*x==0.0) u[6] = 0.0;
      break;

    case 3: // n-down
      u[5]=-1.0;
      break;

    case 4: // t = -104
      u[6] = -acos(-0.25);
      break;

    case 10: // hpd (not really)
      u[6] = acos(-0.25);
      u[0] = sin(u[6]); // Mx
      u[2] = cos(u[6]); // Mz
      u[4] = 1.0;         // Ny
      u[5] = 0.0;         // Nz
      break;

    case 11: // inversed hpd
      u[6]=acos(-0.25);
      u[0] =-sin(u[6]); // Mx
      u[2] = cos(u[6]); // Mz
      u[4] =-1.0;         // Ny
      u[5] = 0.0;         // Nz
      break;

    case 12: // hpd 2pi-soliton (again, not really)
      double w = 0.1; // width
      double p = *x/w; // -1..1
      if (p<-1.0) p=-1.0;
      if (p>+1.0) p=+1.0;
      u[6]=acos(-0.28);
      u[0]=-cos(p*M_PI) * sin(u[6]); // Mx
      u[1]= sin(p*M_PI) * sin(u[6]); // Mx
      u[2]= cos(u[6]);  // Mz
      u[3]=-sin(p*M_PI); // nx
      u[4]=-cos(p*M_PI); // ny
      u[5]= 0.0;         // nz
      break;

    }
  }
}

/******************************************************************/

#ifdef HE3LIB
extern "C" {
#include <he3.h>
}
// set parameters using temperature, pressure and he3lib
void set_he3tp(double ttc, double p){
  double nu_b = he3_nu_b_(&ttc,&p);
  CP0 = he3_cpar_(&ttc,&p) - CPT*tcurr;
  LF0 = nu_b  - LFT*tcurr;
  DF0 = he3_diff_perp_zz_(&ttc,&p,&f0) - DFT*tcurr;
  double tr  = 1.2e-7/sqrt(1.0-ttc);
  TF0   = 1.0/ (4.0*M_PI*M_PI * nu_b*nu_b * tr) - TFT*tcurr;
}
#endif

/******************************************************************/

// Error class for exceptions
class Err {
  std::ostringstream s;
  public:
    Err(){}
    Err(const Err & o) { s << o.s.str(); }
    template <typename T>
    Err & operator<<(const T & o){ s << o; return *this; }
    std::string str()  const { return s.str(); }
};

// check number of arguments for a command, throw Err if it is wrong
void
check_nargs(int nargs, int n1, int n2=-1){
  if (n2<0 && nargs!=n1)
    throw Err() << "command requires "  << n1 << " arguments";
  if (n2>=0 && (nargs < n1 || nargs > n2))
    throw Err() << "command requires "  << n1 << " to " << n2 << " arguments";
}

// Read an argument of arbitrary type from string. Throw Err in case of
// error.
template <typename T> T
get_arg(const std::string &str){
  std::istringstream ss(str);
  T v; ss >> v;
  if (!ss.eof() || ss.fail())
    throw Err() << "can't parse argument \"" << str << "\"";
  return v;
}

// check that there is one argument in the list and read it.
// (for set, step commands)
template <typename T> T
get_one_arg(const std::vector<std::string> & args){
  check_nargs(args.size(), 1);
  return get_arg<T>(args[0]);
}

/******************************************************************/
void
write_profile(pdecol_solver *solver, const std::string & fname) {

  std::vector<double> xsol(npts);
  set_mesh(&pars, xsol);
  std::vector<double> usol = solver->values(xsol, nder);

  std::ofstream ss(fname);
  // print legend: # coord  U(0) U(1) ... U(0)' U(1)' ...
  ss << "# coord.     ";
  for (int d = 0; d<nder; d++){
    ss << "  ";
    for (int n = 0; n<npde; n++){
      ss << "U(" << n << ")";
      for (int i=0; i<d; i++) ss << "'";
      ss << "      ";
    }
  }
  ss << "\n" << std::scientific << std::setprecision(6);

  //print values
  for (int i=0; i< xsol.size(); i++){
    ss << xsol[i];
    for (int d = 0; d<nder; d++){
      ss << "  ";
      for (int n = 0; n<npde; n++){
        ss << " " << solver->get_value(usol, npts, i, n, d);
      }
    }
    ss << "\n";
  }
}

/******************************************************************/

// Process SWEEP command.
// Modify P0, PT, change global variable tend.
template <typename T>
void
cmd_sweep(const std::vector<std::string> & args, T *P0, T *PT, T factor=1){
  check_nargs(args.size(), 2);
  double VD = get_arg<double>(args[0]); // destination
  double R = get_arg<double>(args[1]);  // rate
  double VO = *P0/factor;                // old value
  int steps = abs(rint((VD-VO)/R/tstep));
  if (steps==0) throw Err() << "zero steps for sweep";
  *PT = (VD-VO)/(steps*tstep) * factor;
  *P0 -= tcurr*(*PT);
  tend  = tcurr + steps*tstep;
}


/// Read one or more commends from a stream.
/// Return 1 if file is finished, 0 if more calculations are needed.
int
read_cmd(std::istream &in_c, std::ostream & out_c){
  // reset sweeps
  H0  = H0  + tcurr*HT;  HT=0.0;
  HR0 = HR0 + tcurr*HRT; HRT=0.0;
  TF0 = TF0 + tcurr*TFT; TFT=0.0;
  T10 = T10 + tcurr*T1T; T1T=0.0;
  LF0 = LF0 + tcurr*LFT; LFT=0.0;
  CP0 = CP0 + tcurr*CPT; CPT=0.0;
  DF0 = DF0 + tcurr*DFT; DFT=0.0;

  // Read input string line by line
  // Stop reading is some command increase tend, then
  // calculation is needed.
  std::string line;
  while (tend <= tcurr){

    // return 1 at the end of command file
    if (!getline(in_c, line)){
      out_c << "End of command file, stop calculations\n";
      return 1;
    }

    // remove comments
    size_t c = line.find('#');
    if (c != std::string::npos){
      line = line.substr(0, c);
    }

    // logging to out_c
    out_c << "Command: " << line << "\n";

    // split the line into command and arguments
    std::istringstream in(line);
    std::string cmd;
    in >> cmd;
    if (!in) continue; // empty command

    std::vector<std::string> args;
    while (1) {
      std::string a;
      in >> a;
      if (!in) break;
      args.push_back(a);
    }
    int narg = args.size();

    // Commands can throw Err exceptions. In this case error message
    // from the exception should be printed and next command started.
    try {


      /*******************************************************/
      // solver

      // (Re)start the solver.
      if (cmd == "start") {
        check_nargs(narg, 0);

        // stop solver if it already exists
        if (solver) delete solver;
        // destroy all writers
        pnm_writers.clear();

        // initialize new solver
        tend=tcurr;
        std::vector<double> xbrpt(npts);
        set_mesh(&pars, xbrpt);
        solver = new pdecol_solver(tcurr, mindt, acc, xbrpt, npde);
        if (!solver) throw Err() << "can't create solver";

        // save the mesh
        std::ostringstream ss;
        ss << pref << ".mesh" << cnt_mesh << ".dat";
        save_mesh(&pars, xbrpt, ss.str());
        cnt_mesh++;
        continue;
      }

      // stop the solver
      if (cmd == "stop") {
        check_nargs(narg, 0);
        if (!solver) throw Err() << "solver is not running";
        else delete solver;
        pnm_writers.clear();
        solver=NULL;
        continue;
      }

      // exit the program
      if (cmd == "exit") {
        check_nargs(narg, 0);
        if (solver) delete solver;
        pnm_writers.clear();
        solver=NULL;
        return 1;
      }

      // write function profiles to a file
      if (cmd == "write_profile") {
        check_nargs(narg, 0,1);
        if (!solver) throw Err() << "solver is not running";

        std::string name;
        if (narg>0) {
          name = args[0];
        }
        else{
          std::ostringstream ss;
          ss << pref << ".prof" << cnt_prof << ".dat";
          name = ss.str();
        }
        write_profile(solver, name);
        cnt_prof++;
        continue;
      }

      // Do calculations for some time.
      if (cmd == "wait") {
        double dt = get_one_arg<double>(args);
        if (!solver) throw Err() << "solver is not running";
        tend = tcurr + dt*1e-3;
        continue;
      }

      // Change solver accuracy. If solver is not running,
      // the value will be used after start.
      if (cmd == "acc") {
        acc = get_one_arg<double>(args);
        if (solver) solver->ch_eps(acc);
        continue;
      }

      // Change accuracy (power of two). If solver is not running,
      // the value will be used after start.
      if (cmd == "acc2") {
        double v = get_one_arg<double>(args);
        acc=pow(2, -v);
        if (solver) solver->ch_eps(acc);
        continue;
      }

      // Change min. time step (recommended 1e-10). It can be changes at any time,
      // but real change happenes when the solver is (re)started.
      if (cmd == "mindt") { mindt = get_one_arg<double>(args); continue; }

      // Change number of points. It can be changed at any time,
      // but real change happenes when the solver is (re)started.
      if (cmd == "npts") { npts = get_one_arg<int>(args); continue; }

      // Change time step. Can be changed during calculation.
      if (cmd == "tstep") { tstep = get_one_arg<double>(args); continue; }

      /*******************************************************/
      // pnm writer

      // Initialize a pnm_writer. Solver should be startded.
      if (cmd == "pnm_start") {
        check_nargs(narg, 0,1);
        std::string name = narg>0 ? args[0]: pref + ".pic.pnm";
        if (!solver)
          throw Err() << "can't start pnm_writer if solver is not running";
        if (!pnm_writers.add(name, solver))
          throw Err() << "can't open create pnm_writer";
        continue;
      }

      // pnm_legend: start draw a legend
      if (cmd == "pnm_legend") {
        check_nargs(narg, 0,1);
        std::string name = narg>0 ? args[0]: pref + ".pic.pnm";
        if (!pnm_writers.legend(name, 50, 100))
          throw Err() << "no such writer";
        continue;
      }

      // pnm_hline: draw a horizontal line
      if (cmd == "pnm_hline") {
        check_nargs(narg, 0,1);
        std::string name = narg>0 ? args[0]: pref + ".pic.pnm";
        if (!pnm_writers.hline(name))
          throw Err() << "no such writer";
        continue;
      }

      // stop a pnm_writer
      if (cmd == "pnm_stop") {
        check_nargs(narg, 0,1);
        std::string name = narg>0 ? args[0]: pref + ".pic.pnm";
        if (!pnm_writers.del(name))
          throw Err() << "no such writer";
        continue;
      }


      /*******************************************************/
      // set NMR frequency
      if (cmd == "set_freq") {
        f0 = get_one_arg<double>(args); continue; }

      // Set uniform field [G, from larmor].
      if (cmd == "set_field") {
        H0 = get_one_arg<double>(args); continue; }

      // Uniform field step [G].
      if (cmd == "step_field") {
        H0 += get_one_arg<double>(args); continue; }

      // Set field gradient [G/cm].
      if (cmd == "set_field_grad") {
        HG = get_one_arg<double>(args); continue; }

      // Set field quadratic term [G/cm^2].
      if (cmd == "set_field_quad") {
        HQ = get_one_arg<double>(args); continue; }

      // Sweep uniform field: destination [G], rate [G/s].
      if (cmd == "sweep_field") {
        cmd_sweep(args, &H0, &HT); continue; }

      // Set uniform field in frequency shift units [Hz from NMR freq].
      if (cmd == "set_field_hz") {
        H0 = get_one_arg<double>(args) * 2*M_PI/gyro; continue; }

      // Uniform field step [Hz].
      if (cmd == "step_field_hz") {
        H0 += get_one_arg<double>(args) * 2*M_PI/gyro; continue; }

      // Sweep uniform field: destination [Hz], rate [Hz/s].
      if (cmd == "sweep_field_hz") {
        cmd_sweep(args, &H0, &HT, 2*M_PI/gyro); continue; }

      // Set uniform field in Larmor position units [cm].
      // Gradient term is used to convert field to cm. Quadratic term is not used.
      // Lower Larmor positon means higher field.
      if (cmd == "set_field_cm") {
        if (HG == 0.0) throw Err() << "can't set Larmor position if "
                                      "field gradient is zero";
        H0 = -get_one_arg<double>(args)*HG; continue; }

      // Uniform field step [cm].
      if (cmd == "step_field_cm") {
        if (HG == 0.0) throw Err() << "can't set Larmor position if "
                                      "field gradient is zero";
        H0 -= get_one_arg<double>(args)*HG; continue; }

      // Sweep uniform field: destination [cm], rate [cm/s].
      if (cmd == "sweep_field_cm") {
        if (HG == 0.0) throw Err() << "can't set Larmor position if "
                                      "field gradient is zero";
        cmd_sweep(args, &H0, &HT, -HG); continue; }

      /*******************************************************/

      // RF-field profile, gradient term [1/cm], quadratic term [1/cm^2].
      if (cmd == "set_rf_prof") {
        check_nargs(narg, 2);
        HRGP = get_arg<double>(args[0]);
        HRQP = get_arg<double>(args[1]);
        continue;
      }

      // Set/step/sweep RF field [G].
      if (cmd == "set_rf_field") {
        HR0 = get_one_arg<double>(args); continue; }
      if (cmd == "step_rf_field") {
        HR0 += get_one_arg<double>(args); continue; }
      if (cmd == "sweep_rf_field") {
        cmd_sweep(args, &HR0, &HRT); continue; }

      // Set/sweep relaxation time t_1 [s]
      if (cmd == "set_t1") {
        T10 = get_one_arg<double>(args); continue; }
      if (cmd == "sweep_t1") {
        cmd_sweep(args, &T10, &T1T); continue; }

      // Set/sweep Leggett-Takagi relaxation time tau_f [s]
      if (cmd == "set_tf") {
        TF0 = get_one_arg<double>(args); continue; }
      if (cmd == "sweep_tf") {
        cmd_sweep(args, &TF0, &TFT); continue; }

      // Set/sweep spin diffusion [cm^2/s]
      if (cmd == "set_diff") {
        DF0 = get_one_arg<double>(args); continue; }
      if (cmd == "sweep_diff") {
        cmd_sweep(args, &DF0, &DFT); continue; }

      // Set and sweep spin-wave velocity c_parallel [cm/s]
      if (cmd == "set_cpar") {
        CP0 = get_one_arg<double>(args); continue; }
      if (cmd == "sweep_cpar") {
        cmd_sweep(args, &CP0, &CPT); continue; }

      // Set and sweep Leggett frequency [Hz]
      if (cmd == "set_leggett_freq") {
        LF0 = get_one_arg<double>(args); continue; }
      if (cmd == "sweep_leggett_freq") {
        cmd_sweep(args, &LF0, &LFT); continue; }

      /*******************************************************/

      // commands
      if (cmd == "bcond_type"){ bcond_type = get_one_arg<int>(args); continue;}
      if (cmd == "icond_type"){ icond_type = get_one_arg<int>(args); continue;}

      if (cmd == "CELL_LEN")  { pars.CELL_LEN  = get_one_arg<double>(args); continue;}
      if (cmd == "XMESH_K")   { pars.XMESH_K   = get_one_arg<double>(args); continue;}
      if (cmd == "XMESH_ACC") { pars.XMESH_ACC = get_one_arg<double>(args); continue;}

      if (cmd == "AER")       { pars.AER     = get_one_arg<int>(args); continue;}
      if (cmd == "AER_LEN")   { pars.AER_LEN = get_one_arg<double>(args); continue;}
      if (cmd == "AER_CNT")   { pars.AER_CNT = get_one_arg<double>(args); continue;}
      if (cmd == "AER_TRW")   { pars.AER_TRW = get_one_arg<double>(args); continue;}

      /*******************************************************/

#ifdef HE3LIB
      if (cmd == "set_ttc_press") {
        check_nargs(narg, 2);
        double T = get_arg<double>(args[0]);
        double P = get_arg<double>(args[1]);
        set_he3tp(T, P);
        continue;
      }
#endif

      /*******************************************************/
      throw Err() << "skipping unknown command.\n";
    }
    catch (Err e){
      std::cerr << line << ": " << e.str() << "\n";
    }
  }
  return 0;
}
/********************************************************************/

// write integral magnetization
void
write_magn(std::ostream & s,
           const std::vector<double> zsol,
           const std::vector<double> usol){

  double smx=0, smy=0, smz=0, sz = 0;
  for (int i=0; i<zsol.size()-1; i++){
    double mx1 = usol[i*npde+0];
    double my1 = usol[i*npde+1];
    double mz1 = usol[i*npde+2];

    double mx2 = usol[(i+1)*npde+0];
    double my2 = usol[(i+1)*npde+1];
    double mz2 = usol[(i+1)*npde+2];
    double dz = zsol[i+1]-zsol[i];
    smx += dz*(mx1+mx2)/2;
    smy += dz*(my1+my2)/2;
    smz += dz*(mz1+mz2)/2;
    sz += dz;
  }
  smx/=sz;
  smy/=sz;
  smz/=sz;

  s << std::scientific << std::setprecision(6)
    << tcurr << "  " << H0 + tcurr*HT << "  "
    << smx << " " << smy << " " << smz << "\n";
}

void
write_pars(std::ostream & s){
  s << " T=" << tcurr*1000 << " ms, "
    << "H0=" << H0+HT*tcurr << " G, "
    << "HR=" << 1e3*(HR0+HRT*tcurr) << " mOe, "
    << "LF=" << 1e-3*(LF0+LFT*tcurr) << " kHz, "
    << "CP=" <<  (CP0+CPT*tcurr) << " cm/s, "
    << "DF=" << (DF0+DFT*tcurr) << " cm^2/s, "
    << "TF=" << 1e6*(TF0+TFT*tcurr) << " mks, "
    << "T1=" <<  (T10+T1T*tcurr) << " cm/s, "
    << "\n";
}



/********************************************************************/
int
main(int argc, char *argv[]){
try{

  // program should be run with one argument - name of command file
  if (argc!=2){
    std::cerr << "Usage: " << argv[0] << " <command file>\n";
    return 1;
  }

  // open command file
  std::ifstream in_c(argv[1]);
  if (!in_c.good()){
    std::cerr << "Can't open command file: " << argv[1] << "\n";
    return 1;
  }

  // find prefix (command file name without extension)
  const char * pos = rindex(argv[1], '.');
  pref = pos ? std::string(argv[1], pos-argv[1]) : argv[1];

  // set default parameters
  set_def_pars(&pars);


  std::ofstream out_m(pref + ".magn.dat"); // log total magnetization
  out_m << "# Integral magnetization log: T, LP, Mx, Mx, Mz\n";
  std::ofstream out_l(pref + ".run.log"); // log commands
  out_l << "# Commands and main parameters\n";

  write_pars(out_l);


  // main cycle
  int n=0;
  while (1) {

    // If we reach final time, read new cmd.
    // If it returns 1, finish the program
    if (tcurr >= tend && read_cmd(in_c, out_l)) break;

    // do the next step
    if (solver) {
      tcurr += tstep;
      solver->step(tcurr);

      std::vector<double> xsol(npts);
      set_mesh(&pars, xsol);

      std::vector<double> usol = solver->values(xsol, nder);
      // write results
      write_magn(out_m, xsol, usol);
      pnm_writers.write(xsol, usol);
      write_pars(out_l);
    }

    // flush files
    out_l.flush();
    out_m.flush();
  }

} catch (pdecol_solver::Err e){
  std::cerr << "Error: " << e.str() << "\n";
}
}
