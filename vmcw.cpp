#include <iostream>

extern "C"{ int vmcw_f_(); }

extern "C"{
  int npde_;
  int npts_;
  int nderv_;
  extern struct {
    double *usol_,*xsol_;
  } arrays_;
}

int
main(){

  npde_=7;
  npts_=254;
  nderv_=3;
  std::cout << "npts: " << npts_ << "\n";

  arrays_.usol_ = new double [nderv_*npts_*npde_];
  arrays_.xsol_ = new double [npts_];

  vmcw_f_();
  delete[] arrays_.usol_;
  delete[] arrays_.xsol_;

}