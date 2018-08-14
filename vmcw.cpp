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




#define NPDE 7
#define NDER 3

// fortran functions
extern "C"{
  void writemj_open_(double *usol, double *xsol);
//  void cmd_open_();
//  void cmd_read_();
  void monitor_(double *usol, double *xsol);
  void set_he3pt_();
//  int vmcw_f_(double *usol, double *xsol);
}

extern "C"{
  extern struct pars_t pars_;
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
  pars_.HR0=pars_.HR0+pars_.time*pars_.HR_SWR;       pars_.HR_SWR=0.0;
  pars_.LP0=pars_.LP0+pars_.time*pars_.LP_SWR;       pars_.LP_SWR=0.0;
  pars_.DF0=pars_.DF0+pars_.time*pars_.DF_SWR;       pars_.DF_SWR=0.0;
  pars_.TF0=pars_.TF0+pars_.time*pars_.TF_SWR;       pars_.TF_SWR=0.0;
  pars_.LF0=pars_.LF0+pars_.time*pars_.LF_SWR;       pars_.LF_SWR=0.0;
  pars_.CPAR0=pars_.CPAR0+pars_.time*pars_.CPAR_SWR; pars_.CPAR_SWR = 0.0;

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
    if (cmd_set(args, "beta",      &pars_.BETA,      stage, CMD_INIT)) continue;
    if (cmd_set(args, "IBN",       &pars_.IBN,       stage, CMD_INIT)) continue;
    if (cmd_set(args, "CELL_LEN",  &pars_.CELL_LEN,  stage, CMD_INIT)) continue;
    if (cmd_set(args, "XMESH_K",   &pars_.XMESH_K,   stage, CMD_INIT)) continue;
    if (cmd_set(args, "XMESH_ACC", &pars_.XMESH_ACC, stage, CMD_INIT)) continue;
    if (cmd_set(args, "AER",       &pars_.AER,       stage, CMD_INIT)) continue;
    if (cmd_set(args, "AER_LEN",   &pars_.AER_LEN,   stage, CMD_INIT)) continue;
    if (cmd_set(args, "AER_CNT",   &pars_.AER_CNT,   stage, CMD_INIT)) continue;
    if (cmd_set(args, "AER_TRW",   &pars_.AER_TRW,   stage, CMD_INIT)) continue;

    if (cmd_set(args, "t1c",    &pars_.T1C,    stage, CMD_INIT | CMD_RUN)) continue;
    if (cmd_set(args, "H",      &pars_.H,      stage, CMD_INIT | CMD_RUN)) continue;
    if (cmd_set(args, "grad",   &pars_.grad,   stage, CMD_INIT | CMD_RUN)) continue;
    if (cmd_set(args, "Hr",     &pars_.HR0,    stage, CMD_INIT | CMD_RUN | CMD_SWEEP)) continue;
    if (cmd_set(args, "tstep",  &pars_.tstep,  stage, CMD_INIT | CMD_RUN)) continue;

//    if (args[0] == "temp_press") {
//      if (!check_nargs(line, args.size(), 3)) continue;
//      pars_.TTC   = atof(args[1].c_str());
//      pars_.PRESS = atof(args[2].c_str());
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
      pars_.tend = pars_.time + dt*1e-3;
      out_c << "Wait " << dt << " ms\n";
      return 0;
    }

    if (args[0] == "LP") {
      if (!check_nargs(line, args.size(), 2)) continue;
      double v = atof(args[1].c_str());
      out_c << "larmor position: " << v << " cm\n";
      pars_.LP0 = v - pars_.time*pars_.LP_SWR;
      continue;
    }

    if (args[0] == "LP_ADD") {
      if (!check_nargs(line, args.size(), 2)) continue;
      double v = atof(args[1].c_str());
      out_c << "larmor position step: " << v << " cm\n";
      pars_.LP0 += v;
      continue;
    }

    if (args[0] == "LP_SWEEP_TO") {
      if (!check_nargs(line, args.size(), 3)) continue;
      double v = atof(args[1].c_str());
      double r = fabs(atof(args[2].c_str()));
      out_c << "sweep larmor position to " << v << " cm at " << r << " cm/s\n";
      int steps = abs(rint((v-pars_.LP0)/r/pars_.tstep));

      if (steps==0){
        out_c << "Warning: zero steps, skip the command.\n";
        continue;
      }
      pars_.LP_SWR = (v-pars_.LP0)/(steps*pars_.tstep);
      pars_.LP0   -= pars_.time*pars_.LP_SWR;
      pars_.tend  = pars_.time + steps*pars_.tstep;
      out_c << "  real rate: " << pars_.LP_SWR << " cm/s\n";
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

  s << pars_.time << "  " << pars_.LP0 + pars_.time*pars_.LP_SWR << "  "
    << smx << " " << smy << " " << smz << "\n";
}

void
write_pars(std::ostream & s){
  s << " T=" << pars_.time*1000 << " ms, "
    << "LP=" << pars_.LP0+pars_.LP_SWR*pars_.time << " cm, "
    << "HR=" << 1e3*(pars_.HR0+pars_.HR_SWR*pars_.time) << " mOe, "
    << "TF=" << 1e6*(pars_.TF0+pars_.TF_SWR*pars_.time) << " mks, "
    << "LF=" << 1e-3*(pars_.LF0+pars_.LF_SWR*pars_.time) << " kHz, "
    << "DF=" << (pars_.DF0+pars_.DF_SWR*pars_.time) << " cm^2/s, "
    << "CPAR=" <<  (pars_.CPAR0+pars_.CPAR_SWR*pars_.time) << " cm/s, "
    << "\n";
}



/********************************************************************/
int
main(){
try{

  int npts=257;

  // set default parameters
  set_def_pars(&pars_);

  // allocate memory
  std::vector<double> usol(NDER*npts*NPDE, 0.0);
  std::vector<double> xsol(npts, 0.0);

  std::ifstream in_c("cmd.txt");   // read commands
  std::ofstream out_m("magn.dat"); // log total magnetization
  out_m << "# Integral magnetization log: T, LP, Mx, Mx, Mz\n";
  std::ofstream out_c("cmd_log.dat"); // log commands and main parameters
  out_c << "# Commands and main parameters\n";

  // read all commands until start or eof
  if (read_cmd(in_c, out_c, STAGE_INIT, NULL)) return 0;

  // set up the mesh and save it into a file
  set_mesh(&pars_, xsol);
  save_mesh(&pars_, xsol, "mesh.txt");

  // initialize the solver
  pdecol_solver solver(xsol, usol, pars_.time, 1e-10, pow(2,-20), NPDE);

    write_pars(out_c);

  while (1) {
    if (fabs(pars_.TTC_ST) >  1e-5) {
      pars_.TTC += pars_.TTC_ST;
      set_he3pt_();
    }
    // If we reach final time, read new cmd.
    // If it returns 1, finish the program
    if (pars_.time >= pars_.tend &&
        read_cmd(in_c, out_c, STAGE_RUN, &solver)) return 0;

    // do the next step
    pars_.time += pars_.tstep;
    solver.step(pars_.time);

    // write results
    write_magn(out_m, xsol, usol);
    write_pars(out_c);

    // flush files
    out_c.flush();
    out_m.flush();
  }


} catch (pdecol_solver::Err e){
  std::cerr << "Error: " << e.str() << "\n";
}
}
