#include <vector>
#include <sstream>
#include <iostream>
#include <sys/time.h>
#include "vmcw_pdecol.h"

// C++ wrapper for PDECOL solver
// V.Zavjalov, 07.2018

/********************************************************************/
// PDECOL common blocks
extern "C"{
  extern struct {
    int LOUT; // file descriptor for error messages (3 by default)
  } iounit_;

  extern struct {
    double DTUSED; // STEP SIZE LAST USED (SUCCESSFULLY)
    int NQUSED;    // ORDER LAST USED (SUCCESSFULLY)
    int NSTEP;     // NUMBER OF STEPS TAKEN SO FAR
    int NFE;       // NUMBER OF RESIDUAL EVALUATIONS (RES CALLS) SO FAR
    int NJE;       // NUMBER OF MATRIX EVALUATIONS (PSETIB CALLS) SO FAR
  } gear0_;

  extern struct {
    int NOGAUS; // set to 0 by default
    int MAXDER; // set to 5 by default
  } option_;
}

/********************************************************************/
// PDECOL functions
extern "C"{
  void pdecol_(double *t0, double *t, double *dt, double *XSOL,
        double *EPS, int *NINT, int *KORD, int *NCC, int *NPDE, int *MF,
        int *INDEX, double *WORK, int *IWORK);

  void values_(double *XSOL, double *USOL,
               double *SCTCH, int *NDIM1, int *NDIM2, int *NPTS, int *NDERV, double *WORK);
}
/********************************************************************/
// wrapper class methods

pdecol_solver::pdecol_solver(
  std::vector<double> &XSOL_, std::vector<double> &USOL_,
  double t0_, double dt_, double EPS_,
  int NPDE_, int NDERV_, int KORD_, int NCC_, int MF_, int verbose_
): XSOL(XSOL_), USOL(USOL_), t0(t0_), dt(dt_), EPS(EPS_),
  NPDE(NPDE_), NDERV(NDERV_), KORD(KORD_), NCC(NCC_), MF(MF_), verbose(verbose_){

  NPTS = XSOL.size();
  NINT = NPTS-1;
  INDEX=1;  // type of call (first call)

  // some tests
  if (USOL.size() != NPTS*NPDE*(NDERV+1)) throw Err() <<
    "pdecol_solver: USOL.size() should be equal to NPTS*NPDE*NDERV";
  if (NINT<1 ) throw Err() << "pdecol_solver: NINT should be >= 1";
  if (KORD<3 || KORD>20) throw Err() << "pdecol_solver: KORD should be 3..20";
  if (NCC<2 || NCC>=KORD) throw Err() << "pdecol_solver: NCC should be 2..KORD-1";
  if (NPDE<1) throw Err() << "pdecol_solver: NPDE should be > 0";

  int ML = NPDE*(KORD-1)-1;
  int NCPTS = KORD*NINT-NCC*(NINT-1);
  int MAXDER = option_.MAXDER; // set in /OPTION/

  // prepare working arrays
  int IDIMWORK = KORD + 4*NPDE + 9*NPDE*NPDE + NCPTS*(3*KORD+2)+
                 NPDE*NCPTS*(3*ML + MAXDER + 7);
  WORK = std::vector<double>(IDIMWORK, 0.0);

  int IDIMIWORK = NCPTS*(NPDE+1);
  IWORK = std::vector<int>(IDIMIWORK, 0);
  IWORK[0] = IDIMWORK;
  IWORK[1] = IDIMIWORK;

  SCTCH = std::vector<double>(KORD*(NDERV+1), 0.0);

  iounit_.LOUT = 0; // send messages to stderr by default
}

// It is possible to change accuracy and method during calcultion.
// In this case INDEX is set to 4;
void
pdecol_solver::ch_eps(double new_eps){
  if (new_eps == EPS) return;
  EPS = new_eps;
  if (INDEX==0) INDEX = 4;
}

void
pdecol_solver::ch_mf(int new_mf){
  if (new_mf == MF) return;
  MF = new_mf;
  if (INDEX==0) INDEX = 4;
}

double
pdecol_solver::get_dtused() const { return gear0_.DTUSED;}

int
pdecol_solver::get_nqused() const { return gear0_.NQUSED;}

int
pdecol_solver::get_nstep() const { return gear0_.NSTEP;}

int
pdecol_solver::get_nfe() const { return gear0_.NFE;}

int
pdecol_solver::get_nje() const { return gear0_.NJE;}

int
pdecol_solver::step(double t) {

  struct timeval tt;
  gettimeofday(&tt, NULL);

  if (verbose){
    if (INDEX==1){
      std::cerr << "PDECOL first-call parameters (index==1):\n";
      std::cerr << "  t0:     " << t0 << " -- the inital value of T\n";
      std::cerr << "  dt:     " << t  << " -- the initial step size in T\n";
      std::cerr << "  xleft:  " << *(XSOL.begin())  << " -- left X value\n";
      std::cerr << "  xright: " << *(XSOL.rbegin())  << " -- right X value\n";
      std::cerr << "  eps:    " << EPS  << " -- the relative time error bound\n";
      std::cerr << "  nint:   " << NINT << " -- the number of subintervals\n";
      std::cerr << "  kord:   " << KORD << " -- the order of the piecewise polinomial space to be used\n";
      std::cerr << "  ncc:    " << NCC  << " -- the number of continuity conditins\n";
      std::cerr << "  npde:   " << NPDE << " -- the number of partial differential equations in the system\n";
      std::cerr << "  mf:     " << MF   << " -- the method flag\n";
      std::cerr << "  work size:  " << IWORK[0]  << " -- length of WORK array\n";
      std::cerr << "  iwork size: " << IWORK[1]  << " -- length of IWORK array\n";
    }
    if (INDEX==3){
      std::cerr << "PDECOL one-step call (index==3) parameters:\n";
      std::cerr << "  dt:     " << t  << " -- the maximum step size allowed\n";
    }
    if (INDEX==4){
      std::cerr << "PDECOL call with EPS or MF reset (index==4) parameters:\n";
      std::cerr << "  eps:    " << EPS  << " -- the maximum step size allowed\n";
      std::cerr << "  mf:     " << MF   << " -- the method flag\n";
    }
    std::cerr << tt.tv_sec << "." << tt.tv_usec << ": ";
    std::cerr << "PDECOL: INDEX: " << INDEX <<
                 " TOUT:" << t  << " ";
  }

  pdecol_(&t0,&t,&dt,XSOL.data(),
    &EPS,&NINT,&KORD,&NCC,&NPDE,&MF,&INDEX,
    WORK.data(), IWORK.data() );
  switch (INDEX){
  case  0: break;
  case -1: throw Err() <<
    "PDECOL: THE INTEGRATION WAS HALTED AFTER FAILING TO PASS THE "
    "ERROR TEST EVEN AFTER REDUCING DT BY A FACTOR OF "
    "1.E10 FROM ITS INITIAL VALUE";
  case -2: throw Err() <<
    "PDECOL: AFTER SOME INITIAL SUCCESS, THE INTEGRATION WAS "
    "HALTED EITHER BY REPEATED ERROR TEST FAILURES OR BY "
    "A TEST ON EPS.  TOO MUCH ACCURACY HAS BEEN REQUESTED.";
  case -3: throw Err() <<
    "PDECOL: THE INTEGRATION WAS HALTED AFTER FAILING TO ACHIEVE "
    "CORRECTOR CONVERGENCE EVEN AFTER REDUCING DT BY A "
    "FACTOR OF 1.E10 FROM ITS INITIAL VALUE.";
  case -4: throw Err() <<
    "PDECOL: SINGULAR MATRIX ENCOUNTERED.  PROBABLY DUE TO STORAGE "
    "OVERWRITES.";
  case -5: throw Err() <<
    "PDECOL: INDEX WAS 4 ON INPUT, BUT THE DESIRED CHANGES OF "
    "PARAMETERS WERE NOT IMPLEMENTED BECAUSE TOUT "
    "WAS NOT BEYOND T.  INTERPOLATION TO T = TOUT WAS "
    "PERFORMED AS ON A NORMAL RETURN.  TO TRY AGAIN, "
    "SIMPLY CALL AGAIN WITH INDEX = 4 AND A NEW TOUT.";
  case -6: throw Err() << "PDECOL: ILLEGAL INDEX VALUE.";
  case -7: throw Err() << "PDECOL: ILLEGAL EPS VALUE.";
  case -8: throw Err() <<
    "PDECOL: AN ATTEMPT TO INTEGRATE IN THE WRONG DIRECTION. THE "
    "SIGN OF DT IS WRONG RELATIVE TO T0 AND TOUT.";
  case -9: throw Err() << "PDECOL: DT = 0";
  case -10: throw Err() << "PDECOL: ILLEGAL NINT VALUE.";
  case -11: throw Err() << "PDECOL: ILLEGAL KORD VALUE.";
  case -12: throw Err() << "PDECOL: ILLEGAL NCC VALUE.";
  case -13: throw Err() << "PDECOL: ILLEGAL NPDE VALUE.";
  case -14: throw Err() << "PDECOL: ILLEGAL MF VALUE.";
  case -15: throw Err() << "PDECOL: ILLEGAL BREAKPOINTS - NOT STRICTLY INCREASING.";
  case -16: throw Err() << "PDECOL: INSUFFICIENT STORAGE FOR WORK OR IWORK.";
  default: throw Err() << "unknown error";
  }
  if (verbose) {
    std::cerr << " DT: " << get_dtused() << " NQ: " << get_nqused() << "\n";
  }
  values_(XSOL.data(),USOL.data(),
          SCTCH.data(),&NPDE,&NPTS,&NPTS,&NDERV,WORK.data());
}

