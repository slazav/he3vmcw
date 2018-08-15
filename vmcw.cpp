#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "vmcw_pdecol.h"
#include "vmcw_pars.h"
#include "vmcw_mesh.h"
#include "pnm_writer.h"

#define NPDE 7
#define NDER 3

struct pars_t pars;

/// set parameters for Leggett equations using main parameter structure
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
    *IBN = pars.IBN;
  }
}


bool
check_nargs(const std::string &line, int n1, int n2){
  if (n1==n2) return true;
  std::cerr << "  Warning: skip bad command (wrong number of arguments)\n";
  return false;
}

bool
check_init(const std::string &line, int stage){
  if (stage<1) return true;
  std::cerr << "  Warning: this command works only before starting calculations:\n";
  return false;
}

/// calculation stages
#define STAGE_INIT 1
#define STAGE_RUN   2

/// command types
#define CMD_INIT  1 // can appear before calculations, in STAGE_INIT
#define CMD_RUN   2 // can appear during calculations, in STAGE_RUN
#define CMD_SWEEP 4 // can be used for sweeps

// process SET command
// return 1 if command have been found (correct or wrong)
template <typename T>
bool
cmd_set(const std::vector<std::string> & args, // splitted command line
             const char *name, // parameter name
             T *ref,           // parameter reference
             const int stage,  // calculation stage
             const int cmdtype // type of command
             ){
  if (args.size() < 1 || strcasecmp(args[0].c_str(),name)!=0) return false;

  if (stage == STAGE_INIT && (cmdtype & CMD_INIT == 0)){
    std::cerr << "  Warning: skip bad command (can appear only before calculations)\n";
    return true;
  }
  if (stage == STAGE_RUN && (cmdtype & CMD_RUN == 0)){
    std::cerr << "  Warning: skip bad command (can appear only during calculations)\n";
    return true;
  }
  if (args.size() != 2){
    std::cerr << "  Warning: skip bad command (wrong number of arguments)\n";
    return true;
  }

  std::istringstream ss(args[1]);
  T v;
  ss >> v;
  if (!ss.eof() || ss.fail()){
    std::cerr << "  Warning: skip bad command (unreadable value)\n";
    return true;
  }
  *ref=v;
  return true;
}

/// Read one or more commends from a stream.
/// stage = 0: pre-configure (before solver is started)
/// stage = 1: configuration during solving
/// Return 1 if file is finished, 0 otherwise.
int
read_cmd(std::istream &in_c, std::ostream & out_c, int stage, pdecol_solver *solver){
  // reset sweeps
  pars.HR0=pars.HR0+pars.time*pars.HRT;          pars.HRT=0.0;
  pars.LP0=pars.LP0+pars.time*pars.LP_SWR;       pars.LP_SWR=0.0;
  pars.DF0=pars.DF0+pars.time*pars.DF_SWR;       pars.DF_SWR=0.0;
  pars.TF0=pars.TF0+pars.time*pars.TF_SWR;       pars.TF_SWR=0.0;
  pars.LF0=pars.LF0+pars.time*pars.LF_SWR;       pars.LF_SWR=0.0;
  pars.CPAR0=pars.CPAR0+pars.time*pars.CPAR_SWR; pars.CPAR_SWR = 0.0;

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

    // commands
    if (cmd_set(args, "beta",      &pars.BETA,      stage, CMD_INIT)) continue;
    if (cmd_set(args, "IBN",       &pars.IBN,       stage, CMD_INIT)) continue;
    if (cmd_set(args, "CELL_LEN",  &pars.CELL_LEN,  stage, CMD_INIT)) continue;
    if (cmd_set(args, "XMESH_K",   &pars.XMESH_K,   stage, CMD_INIT)) continue;
    if (cmd_set(args, "XMESH_ACC", &pars.XMESH_ACC, stage, CMD_INIT)) continue;
    if (cmd_set(args, "AER",       &pars.AER,       stage, CMD_INIT)) continue;
    if (cmd_set(args, "AER_LEN",   &pars.AER_LEN,   stage, CMD_INIT)) continue;
    if (cmd_set(args, "AER_CNT",   &pars.AER_CNT,   stage, CMD_INIT)) continue;
    if (cmd_set(args, "AER_TRW",   &pars.AER_TRW,   stage, CMD_INIT)) continue;

    if (cmd_set(args, "t1c",    &pars.T1C,    stage, CMD_INIT | CMD_RUN)) continue;
    if (cmd_set(args, "H",      &pars.H,      stage, CMD_INIT | CMD_RUN)) continue;
    if (cmd_set(args, "grad",   &pars.grad,   stage, CMD_INIT | CMD_RUN)) continue;
    if (cmd_set(args, "Hr",     &pars.HR0,    stage, CMD_INIT | CMD_RUN | CMD_SWEEP)) continue;
    if (cmd_set(args, "Hrg",    &pars.HRG,    stage, CMD_INIT | CMD_RUN | CMD_SWEEP)) continue;
    if (cmd_set(args, "Hrq",    &pars.HRQ,    stage, CMD_INIT | CMD_RUN | CMD_SWEEP)) continue;
    if (cmd_set(args, "DF0",    &pars.DF0,    stage, CMD_INIT | CMD_RUN | CMD_SWEEP)) continue;
    if (cmd_set(args, "LF0",    &pars.LF0,    stage, CMD_INIT | CMD_RUN | CMD_SWEEP)) continue;
    if (cmd_set(args, "CPAR",    &pars.CPAR0,    stage, CMD_INIT | CMD_RUN | CMD_SWEEP)) continue;
    if (cmd_set(args, "tstep",  &pars.tstep,  stage, CMD_INIT | CMD_RUN)) continue;

//    if (args[0] == "temp_press") {
//      if (!check_nargs(line, args.size(), 3)) continue;
//      pars.TTC   = atof(args[1].c_str());
//      pars.PRESS = atof(args[2].c_str());
//      set_he3pt_();
//      continue;
//    }


    if (args[0] == "start") {
      if (!check_nargs(line, args.size(), 1)) continue;
      if (!check_init(line, stage)) continue;
      out_c << "Start calculation\n";
      return 0;
    }
    if (args[0] == "stop") {
      if (!check_nargs(line, args.size(), 1)) continue;
      out_c << "Stop calculations\n";
      return 1;
    }

    if (args[0] == "profile") {
      if (!check_nargs(line, args.size(), 2)) continue;
      out_c << "Write function profiles\n";
      std::ofstream ss(args[1]);
      if (solver) solver->write_profile(ss);
      return 0;
    }


    if (args[0] == "wait") {
      if (!check_nargs(line, args.size(), 2)) continue;
      double dt = atof(args[1].c_str());
      pars.tend = pars.time + dt*1e-3;
      out_c << "Wait " << dt << " ms\n";
      return 0;
    }

    if (args[0] == "LP") {
      if (!check_nargs(line, args.size(), 2)) continue;
      double v = atof(args[1].c_str());
      out_c << "larmor position: " << v << " cm\n";
      pars.LP0 = v - pars.time*pars.LP_SWR;
      continue;
    }

    if (args[0] == "LP_ADD") {
      if (!check_nargs(line, args.size(), 2)) continue;
      double v = atof(args[1].c_str());
      out_c << "larmor position step: " << v << " cm\n";
      pars.LP0 += v;
      continue;
    }

    if (args[0] == "LP_SWEEP_TO") {
      if (!check_nargs(line, args.size(), 3)) continue;
      double v = atof(args[1].c_str());
      double r = fabs(atof(args[2].c_str()));
      out_c << "sweep larmor position to " << v << " cm at " << r << " cm/s\n";
      int steps = abs(rint((v-pars.LP0)/r/pars.tstep));

      if (steps==0){
        out_c << "Warning: zero steps, skip the command.\n";
        continue;
      }
      pars.LP_SWR = (v-pars.LP0)/(steps*pars.tstep);
      pars.LP0   -= pars.time*pars.LP_SWR;
      pars.tend  = pars.time + steps*pars.tstep;
      out_c << "  real rate: " << pars.LP_SWR << " cm/s\n";
      return 0;
    }

    out_c << "  Warning: skipping unknown command.\n";
    continue;
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
    double mx1 = usol[i*NPDE+0];
    double my1 = usol[i*NPDE+1];
    double mz1 = usol[i*NPDE+2];

    double mx2 = usol[(i+1)*NPDE+0];
    double my2 = usol[(i+1)*NPDE+1];
    double mz2 = usol[(i+1)*NPDE+2];
    double dz = zsol[i+1]-zsol[i];
    smx += dz*(mx1+mx2)/2;
    smy += dz*(my1+my2)/2;
    smz += dz*(mz1+mz2)/2;
    sz += dz;
  }
  smx/=sz;
  smy/=sz;
  smz/=sz;

  s << pars.time << "  " << pars.LP0 + pars.time*pars.LP_SWR << "  "
    << smx << " " << smy << " " << smz << "\n";
}

void
write_pars(std::ostream & s){
  s << " T=" << pars.time*1000 << " ms, "
    << "LP=" << pars.LP0+pars.LP_SWR*pars.time << " cm, "
    << "HR=" << 1e3*(pars.HR0+pars.HRT*pars.time) << " mOe, "
    << "TF=" << 1e6*(pars.TF0+pars.TF_SWR*pars.time) << " mks, "
    << "LF=" << 1e-3*(pars.LF0+pars.LF_SWR*pars.time) << " kHz, "
    << "DF=" << (pars.DF0+pars.DF_SWR*pars.time) << " cm^2/s, "
    << "CPAR=" <<  (pars.CPAR0+pars.CPAR_SWR*pars.time) << " cm/s, "
    << "\n";
}



/********************************************************************/
int
main(){
try{

  int npts=257;

  // set default parameters
  set_def_pars(&pars);

  // allocate memory
  std::vector<double> usol(NDER*npts*NPDE, 0.0);
  std::vector<double> xsol(npts, 0.0);

  std::ifstream in_c("cmd.txt");   // read commands
  std::ofstream out_m("magn.dat"); // log total magnetization
  out_m << "# Integral magnetization log: T, LP, Mx, Mx, Mz\n";
  std::ofstream out_c("cmd_log.dat"); // log commands and main parameters
  out_c << "# Commands and main parameters\n";
  std::ofstream out_p("cmd_log.pnm"); // log commands and main parameters
  pnm_writer pnmw(out_p);

  // read all commands until start or eof
  if (read_cmd(in_c, out_c, STAGE_INIT, NULL)) return 0;

  // set up the mesh and save it into a file
  set_mesh(&pars, xsol);
  save_mesh(&pars, xsol, "mesh.txt");

  // initialize the solver
  pdecol_solver solver(xsol, usol, pars.time, 1e-10, pow(2,-20), NPDE);

    write_pars(out_c);

  int n=0;
  while (1) {
//    if (fabs(pars.TTC_ST) >  1e-5) {
//      pars.TTC += pars.TTC_ST;
//      set_he3pt_();
//    }
    // If we reach final time, read new cmd.
    // If it returns 1, finish the program
    if (pars.time >= pars.tend &&
        read_cmd(in_c, out_c, STAGE_RUN, &solver)) break;

    // do the next step
    pars.time += pars.tstep;
    solver.step(pars.time);

    // write results
    write_magn(out_m, xsol, usol);
    pnmw.write(xsol, usol, NPDE);
    write_pars(out_c);

    // flush files
    out_c.flush();
    out_m.flush();
  }

} catch (pdecol_solver::Err e){
  std::cerr << "Error: " << e.str() << "\n";
}
}
