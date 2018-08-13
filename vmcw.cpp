#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "vmcw_pdecol.h"
#include "vmcw_pars.h"
#include "vmcw_mesh.h"

#define STAGE1 1
#define STAGE2 2

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
  std::cerr << "Warning: skip bad command (wrong number of arguments):\n"
            << "> " << line << "\n";
  return false;
}

bool
check_init(const std::string &line, int stage){
  if (stage<1) return true;
  std::cerr << "Warning: this command works only before starting calculations:\n"
            << "> " << line << "\n";
  return false;
}

/// Read one or more commends from a stream.
/// stage = 0: pre-configure (before solver is started)
/// stage = 1: configuration during solving
/// Return 1 if file is finished, 0 otherwise.
int
read_cmd(std::istream &s, int stage){
  // read input string line by line
  std::string line;
  while (getline(s, line)){

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

    // empty string
    if (args.size()<1) continue;

    // commands
    if (args[0] == "beta") {
      if (!check_nargs(line, args.size(), 2)) continue;
      if (!check_init(line, stage)) continue;
      pars_.BETA = atof(args[1].c_str());
      continue;
    }

    if (args[0] == "t1c") {
      if (!check_nargs(line, args.size(), 2)) continue;
      pars_.T1C = atof(args[1].c_str());
      continue;
    }

    if (args[0] == "temp_press") {
      if (!check_nargs(line, args.size(), 3)) continue;
      pars_.TTC   = atof(args[1].c_str());
      pars_.PRESS = atof(args[2].c_str());
      set_he3pt_();
      continue;
    }

    if (args[0] == "H") {
      if (!check_nargs(line, args.size(), 2)) continue;
      pars_.H = atof(args[1].c_str());
      continue;
    }

    if (args[0] == "grad") {
      if (!check_nargs(line, args.size(), 2)) continue;
      pars_.grad = atof(args[1].c_str());
      continue;
    }

    if (args[0] == "Hr") {
      if (!check_nargs(line, args.size(), 2)) continue;
      pars_.HR0 = atof(args[1].c_str());
      continue;
    }

    if (args[0] == "IBN") {
      if (!check_nargs(line, args.size(), 2)) continue;
      if (!check_init(line, stage)) continue;
      pars_.IBN = atoi(args[1].c_str());
      continue;
    }

    if (args[0] == "CELL_LEN") {
      if (!check_nargs(line, args.size(), 2)) continue;
      if (!check_init(line, stage)) continue;
      pars_.CELL_LEN = atof(args[1].c_str());
      std::cerr << "set cell length: " << pars_.CELL_LEN << " cm\n";
      continue;
    }

    if (args[0] == "XMESH_K") {
      if (!check_nargs(line, args.size(), 2)) continue;
      if (!check_init(line, stage)) continue;
      pars_.XMESH_K = atof(args[1].c_str());
      continue;
    }

    if (args[0] == "XMESH_ACC") {
      if (!check_nargs(line, args.size(), 2)) continue;
      if (!check_init(line, stage)) continue;
      pars_.XMESH_ACC = atof(args[1].c_str());
      continue;
    }

    if (args[0] == "AER") {
      if (!check_nargs(line, args.size(), 2)) continue;
      if (!check_init(line, stage)) continue;
      pars_.AER = atoi(args[1].c_str());
      continue;
    }

    if (args[0] == "AER_LEN") {
      if (!check_nargs(line, args.size(), 2)) continue;
      if (!check_init(line, stage)) continue;
      pars_.AER_LEN = atof(args[1].c_str());
      continue;
    }

    if (args[0] == "AER_CNT") {
      if (!check_nargs(line, args.size(), 2)) continue;
      if (!check_init(line, stage)) continue;
      pars_.AER_CNT = atof(args[1].c_str());
      continue;
    }

    if (args[0] == "AER_TRW") {
      if (!check_nargs(line, args.size(), 2)) continue;
      if (!check_init(line, stage)) continue;
      pars_.AER_TRW = atof(args[1].c_str());
      continue;
    }



    if (args[0] == "start") {
      if (!check_nargs(line, args.size(), 1)) continue;
      if (!check_init(line, stage)) continue;
      std::cerr << "Start calculation\n";
      return 0;
    }

    if (args[0] == "tstep") {
      if (!check_nargs(line, args.size(), 2)) continue;
      double dt = atof(args[1].c_str());
      std::cerr << "time step: " << dt << " ms\n";
      pars_.tstep = dt*1e-3;
      continue;
    }

    if (args[0] == "stop") {
      if (!check_nargs(line, args.size(), 1)) continue;
      std::cerr << "Stop calculations\n";
      return 1;
    }

    if (args[0] == "wait") {
      if (!check_nargs(line, args.size(), 2)) continue;
      double dt = atof(args[1].c_str());
      pars_.tend = pars_.time + dt*1e-3;
      std::cerr << "Wait " << dt << " ms\n";
      return 0;
    }

    if (args[0] == "LP") {
      if (!check_nargs(line, args.size(), 2)) continue;
      double v = atof(args[1].c_str());
      std::cerr << "larmor position: " << v << " cm\n";
      pars_.LP0 = v - pars_.time*pars_.LP_SWR;
      continue;
    }

    if (args[0] == "LP_ADD") {
      if (!check_nargs(line, args.size(), 2)) continue;
      double v = atof(args[1].c_str());
      std::cerr << "larmor position step: " << v << " cm\n";
      pars_.LP0 += v;
      continue;
    }

    if (args[0] == "LP_SWEEP_TO") {
      if (!check_nargs(line, args.size(), 3)) continue;
      double v = atof(args[1].c_str());
      double r = fabs(atof(args[2].c_str()));
      std::cerr << "sweep larmor position to " << v << " cm at " << r << " cm/s\n";
      int steps = abs(rint((v-pars_.LP0)/r/pars_.tstep));

      if (steps==0){
        std::cerr << "Warning: zero steps, skip the command\n";
        continue;
      }
      pars_.LP_SWR = (v-pars_.LP0)/(steps*pars_.tstep);
      pars_.LP0   -= pars_.time*pars_.LP_SWR;
      pars_.tend  = pars_.time + steps*pars_.tstep;
      std::cerr << "  real rate: " << pars_.LP_SWR << " cm/s\n";
      return 0;
    }

    std::cerr << "Warning: skip unknown command:\n"
              << "> " << line << "\n";
    continue;
  }
  // end of file
  std::cerr << "End of command file, stop calculations\n";
  return 1;
}


int
main(){
try{

  int npde=7;
  int npts=257;
  int nderv=3;

  // set default parameters
  set_def_pars(&pars_);

  // allocate memory
  std::vector<double> usol(nderv*npts*npde, 0.0);
  std::vector<double> xsol(npts, 0.0);

  std::ifstream cmdf("cmd.txt");

  // read all commands until start or eof
  if (read_cmd(cmdf, 0)) return 0;

  // set up the mesh and save it into a file
  set_mesh(&pars_, xsol);
  save_mesh(&pars_, xsol, "mesh.txt");

  // initialize the solver
  pdecol_solver solver(xsol, usol, pars_.time, 1e-10, pow(2,-20), npde);

  writemj_open_(usol.data(), xsol.data());
  while (1) {
    if (fabs(pars_.TTC_ST) >  1e-5) {
      pars_.TTC += pars_.TTC_ST;
      set_he3pt_();
    }
    // If we reach final time, read new cmd.
    // If it returns 1, finish the program
    if (pars_.time >= pars_.tend && read_cmd(cmdf, 1)) return 0;

    // do the next step
    pars_.time += pars_.tstep;
    solver.step(pars_.time);
    monitor_(usol.data(), xsol.data());
  }


} catch (pdecol_solver::Err e){
  std::cerr << "Error: " << e.str() << "\n";
}
}
