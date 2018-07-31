#include <iostream>
#include <cmath>
#include "vmcw_pdecol.h"

extern "C"{

  void pdecol_init_(double *T); // set PDECOL parameters
  void set_mesh_(double *xsol, int *npts);
  void save_mesh_(double *xsol, int *npts);
  void writemj_open_(double *usol, double *xsol);
  void cmd_open_();
  void cmd_read_();
  void pdecol_run_(double *t, double *usol, double *xsol);
  void monitor_(double *usol, double *xsol);
  void set_he3pt_();
  int vmcw_f_(double *usol, double *xsol);
}

extern "C"{
  int npde_;
  int npts_;
  int nderv_;

  extern struct {
    double *usol;
    double *xsol;
  } arrays_;

  extern struct {
    double t11, grad, H;
    double HR0,HR_SWR;
    double LP0,LP_SWR;
    double DF0,DF_SWR;
    double TF0,TF_SWR;
    double LF0,LF_SWR;
    double CPAR0, CPAR_SWR_;
    double AER, AER_LEN, AER_CNT, AER_TRW;
    double CELL_LEN;
    double XMESH_K, XMESH_ACC;
    double BETA, IBN;
    double PRESS, TTC, TTC_ST, T1C;
  } pars_;

  extern struct {
   double T, TSTEP, TEND;
  } timep_;
}


/// Read one or more commends from a stream.
/// stage = 0: pre-configure (before solver is started)
/// stage = 1: configuration during solving
/// Return 1 if stage is finished, 0 otherwise.
int
read_cmd(std::istream &s, int stage){
}

int
main(){
try{

  npde_=7;
  npts_=257;
  nderv_=3;

  std::cout << "npts: " << npts_ << "\n";

  std::vector<double> usol(nderv_*npts_*npde_, 0.0);
  std::vector<double> xsol(npts_, 0.0);
  arrays_.usol = usol.data();
  arrays_.xsol = xsol.data();

//  read_cfg_("vmcw.cfg");

  vmcw_f_(arrays_.usol, arrays_.xsol);

  timep_.T = 0.0;
  timep_.TSTEP = 5e-3;
  timep_.TEND = timep_.T;
  pars_.LP0=0.0;
  pars_.HR0=1e-3;
  pars_.LP_SWR=0.0;
  pars_.HR_SWR=0.0;

  set_mesh_(arrays_.xsol, &npts_);
  save_mesh_(arrays_.xsol, &npts_);
  writemj_open_(arrays_.usol, arrays_.xsol);
  cmd_open_();
  set_he3pt_();

  pdecol_solver solver(xsol, usol, timep_.T, 1e-10, pow(2,-20), npde_);
  while (1) {
    if (abs(pars_.TTC_ST) >  1e-5) {
      pars_.TTC += pars_.TTC_ST;
      set_he3pt_();
    }
    if (timep_.T >= timep_.TEND) cmd_read_();
    timep_.T += timep_.TSTEP;
    solver.step(timep_.T);
    monitor_(arrays_.usol,arrays_.xsol);
  }


} catch (pdecol_solver::Err e){
  std::cerr << "Error: " << e.str() << "\n";
}
}