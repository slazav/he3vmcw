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
double icond_type = 0;

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

struct pars_t pars; // parameter structure

/* PDECOL solver pointer */
pdecol_solver *solver = NULL;

/* Container for PNM writers. Key is the file name. */
pnm_writers_t pnm_writers;

/******************************************************************/
// Set parameters for Leggett equations using main parameter structure.
// This function is called from F and BNDRY
extern "C" {
  void set_pars_(double *t, double *x,
                 double *Wr, double *Wz, double *W0,
                 double *WB, double *Cpar, double *dCpar,
                 double *Diff, double *Tf, double *T1, int *IBN){

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
      int d0=0, d1=1;
      *Cpar *= 1.0 - 0.5*aer_step(&pars, *x,0);
      *dCpar =(*Cpar) * 0.5*aer_step(&pars, *x,1);
      *Diff *= 1.0 - 0.835 * aer_step(&pars, *x,0);
      *Tf   *= 1.0 - 0.5 * aer_step(&pars, *x,0);
    }
    // type of boundary condition
    *IBN = pars.IBN;
  }

  void set_icond_(int *IIN){
    *IIN = icond_type;
  }
}

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
read_arg(const std::string &str){
  std::istringstream ss(str);
  T v; ss >> v;
  if (!ss.eof() || ss.fail())
    throw Err() << "can't parse argument \"" << str << "\"";
  return v;
}

/******************************************************************/

// process SET command
// return 1 if command have been found (correct or wrong)
template <typename T>
bool
cmd_set(const std::vector<std::string> & args, T *ref){
  check_nargs(args.size(), 1);
  std::istringstream ss(args[0]);
  T v;
  ss >> v;
  if (!ss.eof() || ss.fail())
    throw Err() << "unreadable argument\n";
  *ref=v;
  return true;
}

/// Read one or more commends from a stream.
/// Return 1 if file is finished, 0 otherwise.
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

  // read input string line by line
  std::string line;
  while (getline(in_c, line)){

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
        save_mesh(&pars, xbrpt, "mesh.txt");
        solver = new pdecol_solver(tcurr, mindt, acc, xbrpt, npde);
        if (!solver) throw Err() << "can't create solver";
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
/*
      // write function profiles to a file
      if (cmd == "profile") {
        check_nargs(narg, 1);
        std::ofstream ss(args[0].c_str());
        if (!solver) throw Err() << "solver is not running";
        else solver->write_profile(ss);
        continue;
      }
*/
/*
      // save mesh into file
      if (cmd == "save_mesh") {
        check_nargs(narg, 1);
        if (!solver) throw Err() << "solver is not running";
        else save_mesh(&pars, solver->get_crd_vec(), args[0]);
        continue;
      }
*/
      // Do calculations for some time.
      if (cmd == "wait") {
        check_nargs(narg, 1);
        double dt = read_arg<double>(args[0]);
        if (!solver) throw Err() << "solver is not running";
        tend = tcurr + dt*1e-3;
        return 0;
      }

      // Change solver accuracy. If solver is not running,
      // the value will be used after start.
      if (cmd == "acc") {
        check_nargs(narg, 1);
        acc = read_arg<double>(args[0]);
        if (solver) solver->ch_eps(acc);
        continue;
      }

      // Change accuracy (power of two). If solver is not running,
      // the value will be used after start.
      if (cmd == "acc2") {
        check_nargs(narg, 1);
        double v = read_arg<double>(args[0]);
        acc=pow(2, -v);
        if (solver) solver->ch_eps(acc);
        continue;
      }

      // Change min. time step (recommended 1e-10). It can be changes at any time,
      // but real change happenes when the solver is (re)started.
      if (cmd == "mindt") {
        check_nargs(narg, 1);
        mindt = read_arg<double>(args[0]);
        continue;
      }

      // Change number of points. It can be changed at any time,
      // but real change happenes when the solver is (re)started.
      if (cmd == "npts") {
        check_nargs(narg, 1);
        npts = read_arg<int>(args[0]);
        continue;
      }

      // Change time step. Can be changed during calculation.
      if (cmd == "tstep") {
        check_nargs(narg, 1);
        tstep = read_arg<double>(args[0]);
        continue;
      }


      /*******************************************************/
      // pnm writer

      // Initialize a pnm_writer. Solver should be startded.
      if (cmd == "pnm_start") {
        check_nargs(narg, 1);
        if (!solver)
          throw Err() << "can't start pnm_writer if solver is not running";
        if (!pnm_writers.add(args[0], solver))
          throw Err() << "can't open create pnm_writer";
        continue;
      }

      // pnm_legend: start draw a legend
      if (cmd == "pnm_legend") {
        check_nargs(narg, 1);
        if (!pnm_writers.legend(args[0], 50, 100))
          throw Err() << "no such writer";
        continue;
      }

      // pnm_hline: draw a horizontal line
      if (cmd == "pnm_hline") {
        check_nargs(narg, 1);
        if (!pnm_writers.hline(args[0]))
          throw Err() << "no such writer";
        continue;
      }

      // stop a pnm_writer
      if (cmd == "pnm_stop") {
        check_nargs(narg, 1);
        if (!pnm_writers.del(args[0]))
          throw Err() << "no such writer";
        continue;
      }


      /*******************************************************/
      // set NMR frequency
      if (cmd == "set_freq") {
        check_nargs(narg, 1);
        f0 = read_arg<double>(args[0]);
        continue;
      }

      // Set uniform field [G, from larmor].
      if (cmd == "set_field") {
        check_nargs(narg, 1);
        H0 = read_arg<double>(args[0]);
        continue;
      }

      // Uniform field step [G].
      if (cmd == "step_field") {
        check_nargs(narg, 1);
        H0 += read_arg<double>(args[0]);
        continue;
      }

      // Set field gradient [G/cm].
      if (cmd == "set_field_grad") {
        check_nargs(narg, 1);
        HG = read_arg<double>(args[0]);
        continue;
      }

      // Set field quadratic term [G/cm^2].
      if (cmd == "set_field_quad") {
        check_nargs(narg, 1);
        HQ = read_arg<double>(args[0]);
        continue;
      }

      // Sweep uniform field: destination [G], rate [G/s].
      if (cmd == "sweep_field") {
        check_nargs(narg, 2);
        double v = read_arg<double>(args[0]); // destination, G
        double r = read_arg<double>(args[1]); // rate, G/s
        int steps = abs(rint((v-H0)/r/tstep));
        if (steps==0) throw Err() << "zero steps";
        HT = (v-H0)/(steps*tstep);
        H0 -= tcurr*HT;
        tend  = tcurr + steps*tstep;
        return 0;
      }

      // Set uniform field in frequency shift units [Hz from NMR freq].
      if (cmd == "set_field_hz") {
        check_nargs(narg, 1);
        H0 = read_arg<double>(args[0]) * 2*M_PI/gyro;
        continue;
      }

      // Uniform field step [Hz].
      if (cmd == "step_field_hz") {
        check_nargs(narg, 1);
        H0 += read_arg<double>(args[0]) * 2*M_PI/gyro;
        continue;
      }

      // Sweep uniform field: destination [Hz], rate [Hz/s].
      if (cmd == "sweep_field_hz") {
        check_nargs(narg, 2);
        double v = read_arg<double>(args[0]);
        double r = read_arg<double>(args[1]);
        double o = H0/2.0/M_PI*gyro;
        int steps = abs(rint((v-o)/r/tstep));
        if (steps==0) throw Err() << "zero steps";
        HT = (v-o)/(steps*tstep) * 2*M_PI/gyro;
        H0 -= tcurr*HT;
        tend  = tcurr + steps*tstep;
        return 0;
      }

      // Set uniform field in Larmor position units [cm].
      // Gradient term is used to convert field to cm. Quadratic term is not used.
      // Lower Larmor positon means higher field.
      if (cmd == "set_field_cm") {
        check_nargs(narg, 1);
        if (HG == 0.0) throw Err() << "can't set Larmor position if "
                                      "field gradient is zero";
        H0 = - read_arg<double>(args[0])*HG;
        continue;
      }

      // Uniform field step [cm].
      if (cmd == "step_field_cm") {
        check_nargs(narg, 1);
        if (HG == 0.0) throw Err() << "can't set Larmor position if "
                                      "field gradient is zero";
        H0 -= read_arg<double>(args[0])*HG;
        continue;
      }

      // Sweep uniform field: destination [cm], rate [cm/s].
      if (cmd == "sweep_field_cm") {
        check_nargs(narg, 2);
        if (HG == 0.0) throw Err() << "can't set Larmor position if "
                                      "field gradient is zero";
        double v = read_arg<double>(args[0]); // destination, cm
        double r = read_arg<double>(args[1]); // rate cm/s
        double o = -H0/HG;    // old value, cm
        int steps = abs(rint((v-o)/r/tstep));
        if (steps==0) throw Err() << "zero steps";
        HT = -(v-o)/(steps*tstep)*HG;
        H0 -= tcurr*HT;
        tend  = tcurr + steps*tstep;
        return 0;
      }
      /*******************************************************/

      // RF-field profile, gradient term [1/cm], quadratic term [1/cm^2].
      if (cmd == "set_rf_prof") {
        check_nargs(narg, 1,2);
        HRGP = read_arg<double>(args[0]);
        if (narg>1) HRQP = read_arg<double>(args[1]);
        continue;
      }

      // Set RF field [G].
      if (cmd == "set_rf_field") {
        check_nargs(narg, 1);
        HR0 = read_arg<double>(args[0]);
        continue;
      }

      // Do RF-field step [G].
      if (cmd == "step_rf_field") {
        check_nargs(narg, 1);
        HR0 += read_arg<double>(args[0]);
        continue;
      }

      // Sweep RF field: destination [G], rate [G/s].
      if (cmd == "sweep_rf_field") {
        check_nargs(narg, 2);
        double v = read_arg<double>(args[0]); // destination, G
        double r = read_arg<double>(args[1]); // rate, G/s
        int steps = abs(rint((v-HR0)/r/tstep));
        if (steps==0) throw Err() << "zero steps";
        HRT = (v-HR0)/(steps*tstep);
        HR0 -= tcurr*HRT;
        tend  = tcurr + steps*tstep;
        return 0;
      }


      // Set and sweep relaxation time t_1 [s]
      if (cmd == "set_t1") {
        check_nargs(narg, 1);
        T10 = read_arg<double>(args[0]);
        continue;
      }
      if (cmd == "sweep_t1") {
        check_nargs(narg, 2);
        double v = read_arg<double>(args[0]); // destination, G
        double r = read_arg<double>(args[1]); // rate, G/s
        int steps = abs(rint((v-T10)/r/tstep));
        if (steps==0) throw Err() << "zero steps";
        T1T = (v-T10)/(steps*tstep);
        T10 -= tcurr*T10;
        tend  = tcurr + steps*tstep;
        return 0;
      }

      // Setand sweep Leggett-Takagi relaxation time tau_f [s]
      if (cmd == "set_tf") {
        check_nargs(narg, 1);
        TF0 = read_arg<double>(args[0]);
        continue;
      }
      if (cmd == "sweep_tf") {
        check_nargs(narg, 2);
        double v = read_arg<double>(args[0]); // destination, G
        double r = read_arg<double>(args[1]); // rate, G/s
        int steps = abs(rint((v-TF0)/r/tstep));
        if (steps==0) throw Err() << "zero steps";
        TFT = (v-TF0)/(steps*tstep);
        TF0 -= tcurr*TF0;
        tend  = tcurr + steps*tstep;
        return 0;
      }

      // Set and sweep spin diffusion [cm^2/s]
      if (cmd == "set_diff") {
        check_nargs(narg, 1);
        DF0 = read_arg<double>(args[0]);
        continue;
      }
      if (cmd == "sweep_diff") {
        check_nargs(narg, 2);
        double v = read_arg<double>(args[0]); // destination, G
        double r = read_arg<double>(args[1]); // rate, G/s
        int steps = abs(rint((v-DF0)/r/tstep));
        if (steps==0) throw Err() << "zero steps";
        DFT = (v-DF0)/(steps*tstep);
        DF0 -= tcurr*DF0;
        tend  = tcurr + steps*tstep;
        return 0;
      }

      // Set and sweep spin-wave velocity c_parallel [cm/s]
      if (cmd == "set_cpar") {
        check_nargs(narg, 1);
        CP0 = read_arg<double>(args[0]);
        continue;
      }
      if (cmd == "sweep_cpar") {
        check_nargs(narg, 2);
        double v = read_arg<double>(args[0]); // destination, G
        double r = read_arg<double>(args[1]); // rate, G/s
        int steps = abs(rint((v-CP0)/r/tstep));
        if (steps==0) throw Err() << "zero steps";
        CPT = (v-CP0)/(steps*tstep);
        CP0 -= tcurr*CP0;
        tend  = tcurr + steps*tstep;
        return 0;
      }

      // Set and sweep Leggett frequency [Hz]
      if (cmd == "set_leggett_freq") {
        check_nargs(narg, 1);
        LF0 = read_arg<double>(args[0]);
        continue;
      }
      if (cmd == "sweep_leggett_freq") {
        check_nargs(narg, 2);
        double v = read_arg<double>(args[0]); // destination, G
        double r = read_arg<double>(args[1]); // rate, G/s
        int steps = abs(rint((v-LF0)/r/tstep));
        if (steps==0) throw Err() << "zero steps";
        LFT = (v-LF0)/(steps*tstep);
        LF0 -= tcurr*LF0;
        tend  = tcurr + steps*tstep;
        return 0;
      }

      /*******************************************************/

      // commands
      if (cmd == "IBN")        { cmd_set(args, &pars.IBN       ); continue;}
      if (cmd == "icond_type") { cmd_set(args, &icond_type     ); continue;}
      if (cmd == "CELL_LEN")  { cmd_set(args, &pars.CELL_LEN  ); continue;}
      if (cmd == "XMESH_K")   { cmd_set(args, &pars.XMESH_K   ); continue;}
      if (cmd == "XMESH_ACC") { cmd_set(args, &pars.XMESH_ACC ); continue;}
      if (cmd == "AER")       { cmd_set(args, &pars.AER       ); continue;}
      if (cmd == "AER_LEN")   { cmd_set(args, &pars.AER_LEN   ); continue;}
      if (cmd == "AER_CNT")   { cmd_set(args, &pars.AER_CNT   ); continue;}
      if (cmd == "AER_TRW")   { cmd_set(args, &pars.AER_TRW   ); continue;}

      /*******************************************************/

#ifdef HE3LIB
      if (cmd == "set_ttc_press") {
        check_nargs(narg, 2);
        double T = read_arg<double>(args[0]);
        double P = read_arg<double>(args[1]);
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
  // end of file
  out_c << "End of command file, stop calculations\n";
  return 1;
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

  if (argc!=2){
    std::cerr << "Usage: " << argv[0] << " <command file>\n";
    return 1;
  }

  std::ifstream in_c(argv[1]);   // read commands
  if (!in_c.good()){
    std::cerr << "Can't open command file: " << argv[1] << "\n";
    return 1;
  }

  // set default parameters
  set_def_pars(&pars);


  std::ofstream out_m("magn.dat"); // log total magnetization
  out_m << "# Integral magnetization log: T, LP, Mx, Mx, Mz\n";
  std::ofstream out_c("cmd_log.dat"); // log commands and main parameters
  out_c << "# Commands and main parameters\n";

  write_pars(out_c);


  // main cycle
  int n=0;
  while (1) {

    // If we reach final time, read new cmd.
    // If it returns 1, finish the program
    if (tcurr >= tend && read_cmd(in_c, out_c)) break;

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
      write_pars(out_c);
    }

    // flush files
    out_c.flush();
    out_m.flush();
  }

} catch (pdecol_solver::Err e){
  std::cerr << "Error: " << e.str() << "\n";
}
}
