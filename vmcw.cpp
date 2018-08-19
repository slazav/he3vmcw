#include <iostream>
#include <fstream>
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



struct pars_t pars; // parameter structure

/* PDECOL solver pointer */
pdecol_solver *solver = NULL;

/* Container for PNM writers. Kay is the file name. */
pnm_writers_t pnm_writers;

/******************************************************************/
// Set parameters for Leggett equations using main parameter structure.
// This function is called from F and BNDRY
extern "C" {
  void set_pars_(double *t, double *x,
                 double *Wr, double *Wz, double *W0,
                 double *WB, double *Cpar, double *dCpar,
                 double *Diff, double *Tf, double *T1, int *IBN){
    double gyro = 20378.0; // Gyromagnetic ratio. Use from he3lib?
    *Wr = gyro*(pars.HR0 + pars.HRG*(*x) + pars.HRQ*(*x)*(*x) + pars.HRT*(*t));
    *Wz = gyro*(pars.H + pars.grad*(*x));
    *W0 = gyro*(pars.H + pars.grad*(pars.LP0 + pars.LP_SWR*(*t)));
    *WB = (pars.LF0 + pars.LF_SWR*(*t))*2*M_PI;
    *Cpar = pars.CPAR0 + pars.CPAR_SWR*(*t);
    *dCpar = 0;
    *Diff = pars.DF0 + pars.DF_SWR*(*t);
    *Tf   = pars.TF0 + pars.TF_SWR*(*t);
    *T1   = 1/pars.t11;

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
}

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

/******************************************************************/

// process SET command
// return 1 if command have been found (correct or wrong)
template <typename T>
bool
cmd_set(const std::vector<std::string> & args, // splitted command line
             const char *name, // parameter name
             T *ref            // parameter reference
             ){
  if (args.size() < 1 || strcasecmp(args[0].c_str(),name)!=0) return false;

  check_nargs(args.size(), 2);
  std::istringstream ss(args[1]);
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
  pars.HR0=pars.HR0+tcurr*pars.HRT;          pars.HRT=0.0;
  pars.LP0=pars.LP0+tcurr*pars.LP_SWR;       pars.LP_SWR=0.0;
  pars.DF0=pars.DF0+tcurr*pars.DF_SWR;       pars.DF_SWR=0.0;
  pars.TF0=pars.TF0+tcurr*pars.TF_SWR;       pars.TF_SWR=0.0;
  pars.LF0=pars.LF0+tcurr*pars.LF_SWR;       pars.LF_SWR=0.0;
  pars.CPAR0=pars.CPAR0+tcurr*pars.CPAR_SWR; pars.CPAR_SWR = 0.0;

  // read input string line by line
  std::string line;
  while (getline(in_c, line)){

    // remove comments
    size_t c = line.find('#');
    if (c != std::string::npos){
      line = line.substr(0, c);
    }

    // split the line into words
    std::istringstream in(line);
    std::vector<std::string> args;
    while (1) {
      std::string a;
      in >> a;
      if (!in) break;
      args.push_back(a);
    }

    // logging to out_c
    out_c << "Command: " << line << "\n";

    // empty string
    if (args.size()<1) continue;

    // Commands can throw Err exceptions. In this case error message
    // from the exception should be printed and next command started.
    try {

      // commands
      if (cmd_set(args, "beta",      &pars.BETA      )) continue;
      if (cmd_set(args, "IBN",       &pars.IBN       )) continue;
      if (cmd_set(args, "CELL_LEN",  &pars.CELL_LEN  )) continue;
      if (cmd_set(args, "XMESH_K",   &pars.XMESH_K   )) continue;
      if (cmd_set(args, "XMESH_ACC", &pars.XMESH_ACC )) continue;
      if (cmd_set(args, "AER",       &pars.AER       )) continue;
      if (cmd_set(args, "AER_LEN",   &pars.AER_LEN   )) continue;
      if (cmd_set(args, "AER_CNT",   &pars.AER_CNT   )) continue;
      if (cmd_set(args, "AER_TRW",   &pars.AER_TRW   )) continue;

      if (cmd_set(args, "t1c",    &pars.T1C    )) continue;
      if (cmd_set(args, "H",      &pars.H      )) continue;
      if (cmd_set(args, "grad",   &pars.grad   )) continue;
      if (cmd_set(args, "Hr",     &pars.HR0    )) continue;
      if (cmd_set(args, "Hrg",    &pars.HRG    )) continue;
      if (cmd_set(args, "Hrq",    &pars.HRQ    )) continue;
      if (cmd_set(args, "DF0",    &pars.DF0    )) continue;
      if (cmd_set(args, "LF0",    &pars.LF0    )) continue;
      if (cmd_set(args, "CPAR",    &pars.CPAR0 )) continue;
      if (cmd_set(args, "tstep",  &tstep       )) continue;

//    if (args[0] == "temp_press") {
//      check_nargs(args.size(), 3);
//      pars.TTC   = atof(args[1].c_str());
//      pars.PRESS = atof(args[2].c_str());
//      set_he3pt_();
//      continue;
//    }

      /*******************************************************/
      // solver

      // (re)start the solver
      if (args[0] == "start") {
        check_nargs(args.size(), 1);

        // stop solver if it already exists
        if (solver) delete solver;

        // initialize new solver
        tend=tcurr;
        solver = new pdecol_solver(tcurr, mindt, acc, npts, npde);
        if (!solver) throw pdecol_solver::Err() << "can't create solver";
        // set up the mesh and save it into a file
        set_mesh(&pars, solver->get_crd_vec());
        continue;
      }

      // stop the solver
      if (args[0] == "stop") {
        check_nargs(args.size(), 1);
        if (solver) delete solver;
        solver=NULL;
        continue;
      }

      // exit the program
      if (args[0] == "exit") {
        check_nargs(args.size(), 1);
        if (solver) delete solver;
        solver=NULL;
        return 1;
      }

      // write function profiles to a file
      if (args[0] == "profile") {
        check_nargs(args.size(), 2);
        std::ofstream ss(args[1].c_str());
        if (solver) solver->write_profile(ss);
        continue;
      }

      // save mesh into file
      if (args[0] == "save_mesh") {
        check_nargs(args.size(), 2);
        if (solver) save_mesh(&pars, solver->get_crd_vec(), args[1]);
        continue;
      }

      // do calculations for some time
      if (args[0] == "wait") {
        check_nargs(args.size(), 2);
        double dt = atof(args[1].c_str());
        tend = tcurr + dt*1e-3;
        return 0;
      }

      // Change solver accuracy.
      if (args[0] == "acc") {
        check_nargs(args.size(), 2);
        acc = atof(args[1].c_str());
        if (solver) solver->ch_eps(acc);
        continue;
      }

      // change accuracy (power of two)
      if (args[0] == "acc2") {
        check_nargs(args.size(), 2);
        double v = atof(args[1].c_str());
        acc=pow(2, -v);
        if (solver) solver->ch_eps(acc);
        continue;
      }

      // Change min. time step (recommended 1e-10). It can be changes at any time,
      // but real change happenes when the solver is (re)started.
      if (args[0] == "mindt") {
        check_nargs(args.size(), 2);
        mindt = atof(args[1].c_str());
        continue;
      }

      // Change number of points. It can be changed at any time,
      // but real change happenes when the solver is (re)started.
      if (args[0] == "npts") {
        check_nargs(args.size(), 2);
        npts = atoi(args[1].c_str());
        continue;
      }

      /*******************************************************/
      // pnm writer

      // initialize a pnm_writer
      if (args[0] == "pnm_start") {
        check_nargs(args.size(), 2);
        if (!pnm_writers.add(args[1]))
          throw Err() << "can't open create pnm_writer";
        continue;
      }

      // pnm_legend: start draw a legend
      if (args[0] == "pnm_legend") {
        check_nargs(args.size(), 2);
        if (!pnm_writers.legend(args[1], 50, 100))
          throw Err() << "no such writer";
        continue;
      }

      // pnm_hline: draw a horizontal line
      if (args[0] == "pnm_hline") {
        check_nargs(args.size(), 2);
        if (!pnm_writers.hline(args[1]))
          throw Err() << "no such writer";
        continue;
      }

      // stop a pnm_writer
      if (args[0] == "pnm_stop") {
        check_nargs(args.size(), 2);
        if (!pnm_writers.del(args[1]))
          throw Err() << "no such writer";
        continue;
      }

      /*******************************************************/

      if (args[0] == "LP") {
        check_nargs(args.size(), 2);
        double v = atof(args[1].c_str());
        out_c << "larmor position: " << v << " cm\n";
        pars.LP0 = v - tcurr*pars.LP_SWR;
        continue;
      }

      if (args[0] == "LP_ADD") {
        check_nargs(args.size(), 2);
        double v = atof(args[1].c_str());
        out_c << "larmor position step: " << v << " cm\n";
        pars.LP0 += v;
        continue;
      }

      if (args[0] == "LP_SWEEP_TO") {
        check_nargs(args.size(), 3);
        double v = atof(args[1].c_str());
        double r = fabs(atof(args[2].c_str()));
        out_c << "sweep larmor position to " << v << " cm at " << r << " cm/s\n";
        int steps = abs(rint((v-pars.LP0)/r/tstep));

        if (steps==0){
          out_c << "Warning: zero steps, skip the command.\n";
          continue;
        }
        pars.LP0 = pars.LP0+tcurr*pars.LP_SWR;
        pars.LP_SWR = (v-pars.LP0)/(steps*tstep);
        pars.LP0   -= tcurr*pars.LP_SWR;
        tend  = tcurr + steps*tstep;
        out_c << "  real rate: " << pars.LP_SWR << " cm/s\n";
        return 0;
      }
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

  s << tcurr << "  " << pars.LP0 + tcurr*pars.LP_SWR << "  "
    << smx << " " << smy << " " << smz << "\n";
}

void
write_pars(std::ostream & s){
  s << " T=" << tcurr*1000 << " ms, "
    << "LP=" << pars.LP0+pars.LP_SWR*tcurr << " cm, "
    << "HR=" << 1e3*(pars.HR0+pars.HRT*tcurr) << " mOe, "
    << "TF=" << 1e6*(pars.TF0+pars.TF_SWR*tcurr) << " mks, "
    << "LF=" << 1e-3*(pars.LF0+pars.LF_SWR*tcurr) << " kHz, "
    << "DF=" << (pars.DF0+pars.DF_SWR*tcurr) << " cm^2/s, "
    << "CPAR=" <<  (pars.CPAR0+pars.CPAR_SWR*tcurr) << " cm/s, "
    << "\n";
}



/********************************************************************/
int
main(){
try{


  // set default parameters
  set_def_pars(&pars);

  std::ifstream in_c("cmd.txt");   // read commands
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

//    if (fabs(pars.TTC_ST) >  1e-5) {
//      pars.TTC += pars.TTC_ST;
//      set_he3pt_();
//    }

    // do the next step
    if (solver) {
      tcurr += tstep;
      solver->step(tcurr);

      // write results
      write_magn(out_m, solver->get_crd_vec(), solver->get_sol_vec());
      pnm_writers.write(solver->get_crd_vec(), solver->get_sol_vec(), npde);
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
