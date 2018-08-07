#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "vmcw_pdecol.h"
#include "vmcw_pars.h"
#include "vmcw_mesh.h"

// fortran functions
extern "C"{
  void writemj_open_(double *usol, double *xsol);
  void cmd_open_();
  void cmd_read_();
  void monitor_(double *usol, double *xsol);
  void set_he3pt_();
  int vmcw_f_(double *usol, double *xsol);
}

extern "C"{
  extern struct pars_t pars_;

  extern struct {
   double T, TSTEP, TEND;
  } timep_;
}

/// Read one or more commends from a stream.
/// stage = 0: pre-configure (before solver is started)
/// stage = 1: configuration during solving
/// Return 1 if stage is finished (solver started or stopped), 0 otherwise.
int
read_cmd(std::istream &s, int stage){
  while (!s.eof()){
    std::string l = getline(s);
  }
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

  ifstream cmdf("cmd.txt");
  read_cmd(cmdf, 0);

//  read_cfg_("vmcw.cfg");

  vmcw_f_(usol.data(), xsol.data());

  timep_.T = 0.0;
  timep_.TSTEP = 5e-3;
  timep_.TEND = timep_.T;
  pars_.LP0=0.0;
  pars_.HR0=1e-3;
  pars_.LP_SWR=0.0;
  pars_.HR_SWR=0.0;

  set_mesh(&pars_, xsol);
  save_mesh(&pars_, xsol, "mesh.txt");
  writemj_open_(usol.data(), xsol.data());
  cmd_open_();
  set_he3pt_();

  pdecol_solver solver(xsol, usol, timep_.T, 1e-10, pow(2,-20), npde);
  while (1) {
    if (fabs(pars_.TTC_ST) >  1e-5) {
      pars_.TTC += pars_.TTC_ST;
      set_he3pt_();
    }
    if (timep_.T >= timep_.TEND) cmd_read_();
    timep_.T += timep_.TSTEP;
    solver.step(timep_.T);
    monitor_(usol.data(), xsol.data());
  }


} catch (pdecol_solver::Err e){
  std::cerr << "Error: " << e.str() << "\n";
}
}