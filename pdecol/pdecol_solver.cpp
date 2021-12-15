#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include "pdecol_solver.h"

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

  // other blocks are used only to save/restore state
  extern struct {
    double T,DTC,DTMN,DTMX,EPSC,UROUND;
    int N,MFC,KFLAG,JSTART;
  } gear1_;

  extern struct {
    double EPSJ,R0;
    int ML,MU,MW,NM1,N0ML,N0W;
  } gear9_;

  extern struct {
    int v[18];
  } istart_;

  extern struct {
    int NINT,KORD,NCC,NPDE,NCPTS,NEQN,IQUAD;
  } sizes_;

  // internal data (originally use SAVE statements)
  extern struct {
    int ILEFT, MFLAG;
  } values_save_;

  extern struct {
    int ILEFT;
  } colpnt_save_;

  extern struct {
    double DELTAM, DELTAP;
    int J;
  } bsplvn_save_;

  extern struct {
    int ILO;
  } interv_save_;

  extern struct { // stifib function - core integrator
     double BND,CON,CRATE,D,D1,E,EDN,ENQ1,ENQ2,ENQ3,EPSOLD;
     double EUP,FN,HOLD,OLDL0,PR1,PR2,PR3,R1,RC,RH,RMAX,TOLD;
     int    I,IDOUB,IER,IREDO,IRET,IWEVAL,J,J1,J2,L,LMAX,M,MEO,METH;
     int    MFOLD,MIO,MITER,NEWQ,NOLD,NQ,NSTEPJ;
  } stifib_save_;

}

/********************************************************************/
// PDECOL functions
extern "C"{
  void pdecol_(double *t0, double *t, double *dt, double *XBKPT,
        double *EPS, int *NINT, int *KORD, int *NCC, int *NPDE, int *MF,
        int *INDEX, double *WORK, int *IWORK);

  void values_(double *XSOL, double *USOL,
               double *SCTCH, int *NDIM1, int *NDIM2, int *NPTS, int *NDERV, double *WORK);

  void findeq_(double *T, double *X, double *UVAL, int *NPTS, int *MSG_LVL, int *err);
}
/********************************************************************/
// wrapper class methods

// Constructor. Do the first call of PCECOL
pdecol_solver::pdecol_solver(
  double t0_, double mindt_, double EPS_, std::vector<double> & XBKPT_,
  int NPDE_, int KORD_, int NCC_, int MF_, int verbose_
): t0(t0_), mindt(mindt_), EPS(EPS_), XBKPT(XBKPT_), NPDE(NPDE_),
   KORD(KORD_), NCC(NCC_), MF(MF_), verbose(verbose_){

  // number of intervals
  NINT = XBKPT.size()-1;

  // some tests
  if (NINT<1 ) throw Err() << "pdecol_solver: NINT should be >= 1";
  if (NPDE<1) throw Err() << "pdecol_solver: NPDE should be > 0";
  if (KORD<3 || KORD>20) throw Err() << "pdecol_solver: KORD should be 3..20";
  if (NCC<2 || NCC>=KORD) throw Err() << "pdecol_solver: NCC should be 2..KORD-1";

  // prepare working arrays
  int ML = NPDE*(KORD-1)-1;
  int NCPTS = KORD*NINT-NCC*(NINT-1);
  int MAXDER = option_.MAXDER; // set in /OPTION/

#ifdef SOLVER_EPDE_DP
  int IDIMWORK = KORD + 4*NPDE + NCPTS*(3*KORD + 2) +
                 NPDE*NCPTS*(MAXDER + 6) + NPDE*NPDE*(13+KORD*(KORD-NCC)*NINT);
#endif
#ifdef SOLVER_PDE_DP
  int IDIMWORK = KORD + 4*NPDE + 9*NPDE*NPDE + NCPTS*(3*KORD+2)+
                 NPDE*NCPTS*(3*ML + MAXDER + 7);
#endif

  WORK = std::vector<double>(IDIMWORK, 0.0);

  int IDIMIWORK = NCPTS*(NPDE+1);
  IWORK = std::vector<int>(IDIMIWORK, 0);
  IWORK[0] = IDIMWORK;
  IWORK[1] = IDIMIWORK;

  // send messages to stderr by default
  iounit_.LOUT = 0;

  INDEX=1;  // type of call (first call)

  // Do first step with minimum length and reset time to zero
  // This is needed to set all solver parameters and avoid some
  // strange state when solver exists, but some functions
  // (like values()) are not available.
  double t = t0+mindt;
  pdecol_(&t0,&t,&mindt,XBKPT.data(),
    &EPS,&NINT,&KORD,&NCC,&NPDE,&MF,&INDEX,
    WORK.data(), IWORK.data() );
}

// restart running solver
void
pdecol_solver::restart() {
  t0 = t;
  IWORK[0] = WORK.size();
  IWORK[1] = IWORK.size();
  INDEX=1;
}

void
pdecol_solver::reset_time() {t0=t=gear1_.T = 0;}

// step
void
pdecol_solver::step(double t_, bool exact) {
  if (exact) INDEX=2;
  t=t_;
  // Here INDEX should be 1, 0 or 4. For call types 2, 3
  // additional methods should be done (?).
  if (verbose) print_index_info();

  // Run pdecol. Some parameters are used only in the first call INDEX=1.
  double dt = mindt; // don't want to change original value (for restarting)
  pdecol_(&t0,&t,&dt,XBKPT.data(),
    &EPS,&NINT,&KORD,&NCC,&NPDE,&MF,&INDEX,
    WORK.data(), IWORK.data() );

  // check INDEX after run (0 or error)
  check_error();

  if (verbose) {
    std::cerr << " DT: " << get_dtused() << " NQ: " << get_nqused() << "\n";
  }
}

std::vector<double>
pdecol_solver::find_eq(std::vector<double> & X) {
  int msg_lvl = 1;
  int n = X.size(), err;
  auto U = values(X,0);
  auto U1 = U;
  findeq_(&t0, X.data(), U.data(), &n, &msg_lvl, &err);
  switch (err) {
    case -1: throw Err() << "find_eq/tn: error in input parameters";
    case +2: throw Err() << "find_eq/tn: more then MAXFUN evaluations";
    case +3: std::cerr << "find_eq/tn: line search failed to find lower point\n"; break;
    case 101: throw Err() << "find_eq: too few points";
  }

  for (size_t i1 = 0; i1<X.size(); ++i1){
  std::cout << i1;
    for (size_t i2 = 0; i2 < NPDE; ++i2)
      std::cout << "   " << U1[i2 + i1*NPDE] << " " << U[i2 + i1*NPDE];
    std::cout << "\n";
  }

  return U;
}

// values
void
pdecol_solver::values(std::vector<double> & xsol,
                      std::vector<double> & usol, int NDERV){

  int NPTS = xsol.size();

  // array for values
  if (usol.size()<NPTS*NPDE*(NDERV+1)) usol.resize(NPTS*NPDE*(NDERV+1));

  // Working storage array
  std::vector<double> SCTCH(KORD*(NDERV+1), 0.0);

  values_(xsol.data(),usol.data(),
          SCTCH.data(),&NPDE,&NPTS,&NPTS,&NDERV,WORK.data());
}


std::vector<double>
pdecol_solver::values(std::vector<double> & xsol, int NDERV){
  std::vector<double> usol;
  values(xsol, usol, NDERV);
  return usol;
}

/********************************************************************/
// save/read state

void
pdecol_solver::save_state(const std::string & fname){
#define WRITE_DATA(S) ff.write((char *)&S, sizeof(S));
  int version = 2;
  std::ofstream ff (fname.c_str(),std::ofstream::binary);
  WRITE_DATA(version);
  WRITE_DATA(verbose);
  WRITE_DATA(INDEX);
  WRITE_DATA(NINT);
  WRITE_DATA(NPDE);
  WRITE_DATA(KORD);
  WRITE_DATA(NCC);
  WRITE_DATA(MF);
  WRITE_DATA(EPS);
  WRITE_DATA(t);
  WRITE_DATA(t0);
  WRITE_DATA(mindt);

  int s;
  s = XBKPT.size();
  WRITE_DATA(s);
  if (s>0) ff.write((char *)XBKPT.data(), s*sizeof(XBKPT[0]));

  s = WORK.size();
  WRITE_DATA(s);
  if (s>0) ff.write((char *)WORK.data(), s*sizeof(WORK[0]));

  s = IWORK.size();
  WRITE_DATA(s);
  if (s>0) ff.write((char *)IWORK.data(), s*sizeof(IWORK[0]));

  // pdecol structures
  WRITE_DATA(iounit_);
  WRITE_DATA(istart_);
  WRITE_DATA(sizes_);
  WRITE_DATA(gear0_);
  WRITE_DATA(gear1_);
  WRITE_DATA(gear9_);
  WRITE_DATA(option_);
  // only in version 2
  WRITE_DATA(values_save_);
  WRITE_DATA(colpnt_save_);
  WRITE_DATA(bsplvn_save_);
  WRITE_DATA(interv_save_);
  WRITE_DATA(stifib_save_);
}

void
pdecol_solver::load_state(const std::string & fname){
#define READ_DATA(S) ff.read((char *)&S, sizeof(S));
  std::ifstream ff (fname.c_str(),std::ifstream::binary);
  if (ff.fail()) throw Err() << "can't read file: " << fname;
  int version;
  ff.read((char *)&version, sizeof(version));

  int NPDE0 = NPDE;
  switch (version){
    case 1:
    case 2:
      READ_DATA(verbose);
      READ_DATA(INDEX);
      READ_DATA(NINT);
      READ_DATA(NPDE);
      READ_DATA(KORD);
      READ_DATA(NCC);
      READ_DATA(MF);
      READ_DATA(EPS);
      READ_DATA(t);
      READ_DATA(t0);
      READ_DATA(mindt);

      if (NPDE!=NPDE0) throw Err() << "can't load state with NPDE=" << NPDE;

      int s;
      READ_DATA(s);
      XBKPT.resize(s);
      if (s>0) ff.read((char *)XBKPT.data(), s*sizeof(XBKPT[0]));

      READ_DATA(s);
      WORK.resize(s);
      if (s>0) ff.read((char *)WORK.data(), s*sizeof(WORK[0]));

      READ_DATA(s);
      IWORK.resize(s);
      if (s>0) ff.read((char *)IWORK.data(), s*sizeof(IWORK[0]));

      // pdecol structures
      READ_DATA(iounit_);
      READ_DATA(istart_);
      READ_DATA(sizes_);
      READ_DATA(gear0_);
      READ_DATA(gear1_);
      READ_DATA(gear9_);
      READ_DATA(option_);
      if (version==1) break;
      // only in version 2
      READ_DATA(values_save_);
      READ_DATA(colpnt_save_);
      READ_DATA(bsplvn_save_);
      READ_DATA(interv_save_);
      READ_DATA(stifib_save_);
      break;
    default:
      throw Err() << "unsupported state file version: " << version;
  }
}


/********************************************************************/
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

/********************************************************************/
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


/********************************************************************/
void
pdecol_solver::check_error(){
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
  case -5:
    std::cerr <<
    "PDECOL: INDEX WAS 4 ON INPUT, BUT THE DESIRED CHANGES OF "
    "PARAMETERS WERE NOT IMPLEMENTED BECAUSE TOUT "
    "WAS NOT BEYOND T.  INTERPOLATION TO T = TOUT WAS "
    "PERFORMED AS ON A NORMAL RETURN.  TO TRY AGAIN, "
    "SIMPLY CALL AGAIN WITH INDEX = 4 AND A NEW TOUT.";
    INDEX=4; // try again with INDEX=4
    break;
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
}

void
pdecol_solver::print_index_info() const{
  struct timeval tt;
  gettimeofday(&tt, NULL);

  switch (INDEX){
  case 0:
  case 2:
    std::cerr << tt.tv_sec << "." << std::setw(6) << std::setfill('0') << tt.tv_usec << ": ";
    std::cerr << "PDECOL: INDEX: " << INDEX <<
                 " TOUT:" << std::scientific << std::setprecision(4) << t  << " ";
    break;

  case 1:
    std::cerr << "PDECOL first-call parameters (index==1):\n";
    std::cerr << "  t0:     " << t0 << " -- the inital value of T\n";
    std::cerr << "  mindt:  " << mindt << " -- minimum time step\n";
    std::cerr << "  xleft:  " << *(XBKPT.begin())  << " -- left X value\n";
    std::cerr << "  xright: " << *(XBKPT.rbegin())  << " -- right X value\n";
    std::cerr << "  eps:    " << EPS  << " -- the relative time error bound\n";
    std::cerr << "  nint:   " << NINT << " -- the number of subintervals\n";
    std::cerr << "  kord:   " << KORD << " -- the order of the piecewise polinomial space to be used\n";
    std::cerr << "  ncc:    " << NCC  << " -- the number of continuity conditins\n";
    std::cerr << "  npde:   " << NPDE << " -- the number of partial differential equations in the system\n";
    std::cerr << "  mf:     " << MF   << " -- the method flag\n";
    std::cerr << "  work size:  " << IWORK[0]  << " -- length of WORK array\n";
    std::cerr << "  iwork size: " << IWORK[1]  << " -- length of IWORK array\n";
    break;

  case 3:
    std::cerr << "PDECOL one-step call (index==3) parameters:\n";
    std::cerr << "  dt:     " << t  << " -- the maximum step size allowed\n";
    break;

  case 4:
    std::cerr << "PDECOL call with EPS or MF reset (index==4) parameters:\n";
    std::cerr << "  eps:    " << EPS  << " -- the relative time error bound\n";
    std::cerr << "  mf:     " << MF   << " -- the method flag\n";
    break;

  default:
    throw Err() << "bad or unknown INDEX setting in pdecol_solver::step: " << INDEX;
  }
}
