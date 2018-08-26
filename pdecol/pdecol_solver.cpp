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

  // other bloks are used only to save/restore state
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
}

/********************************************************************/
// PDECOL functions
extern "C"{
  void pdecol_(double *t0, double *t, double *dt, double *XBKPT,
        double *EPS, int *NINT, int *KORD, int *NCC, int *NPDE, int *MF,
        int *INDEX, double *WORK, int *IWORK);

  void values_(double *XSOL, double *USOL,
               double *SCTCH, int *NDIM1, int *NDIM2, int *NPTS, int *NDERV, double *WORK);
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
}

// restart running solver
void
pdecol_solver::restart() {
  t0 = t;
  IWORK[0] = WORK.size();
  IWORK[1] = IWORK.size();
  INDEX=1;
}

// step
int
pdecol_solver::step(double t_) {
  t=t_;
  // Here INDEX should be 1, 0 or 4. For call types 2, 3
  // additional methods should be done (?).
  if (verbose) print_index_info(INDEX);

  // Run pdecol. Some parameters are used only in the first call INDEX=1.
  double dt = mindt; // don't want to change original value (for restarting)
  pdecol_(&t0,&t,&dt,XBKPT.data(),
    &EPS,&NINT,&KORD,&NCC,&NPDE,&MF,&INDEX,
    WORK.data(), IWORK.data() );

  // check INDEX after run (0 or error)
  check_error(INDEX);

  if (verbose) {
    std::cerr << " DT: " << get_dtused() << " NQ: " << get_nqused() << "\n";
  }
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
  int version = 1;
  std::ofstream ff (fname.c_str(),std::ofstream::binary);
  ff.write((char *)&version, sizeof(version));
  ff.write((char *)&verbose, sizeof(verbose));
  ff.write((char *)&INDEX,   sizeof(INDEX));
  ff.write((char *)&NINT,    sizeof(NINT));
  ff.write((char *)&NPDE,    sizeof(NPDE));
  ff.write((char *)&KORD,    sizeof(KORD));
  ff.write((char *)&NCC,     sizeof(NCC));
  ff.write((char *)&MF,      sizeof(MF));
  ff.write((char *)&EPS,     sizeof(EPS));
  ff.write((char *)&t,       sizeof(t));
  ff.write((char *)&t0,      sizeof(t0));
  ff.write((char *)&mindt,   sizeof(mindt));

  int s;
  s = XBKPT.size();
  ff.write((char *)&s,   sizeof(s));
  if (s>0) ff.write((char *)XBKPT.data(), s*sizeof(XBKPT[0]));

  s = WORK.size();
  ff.write((char *)&s,   sizeof(s));
  if (s>0) ff.write((char *)WORK.data(), s*sizeof(WORK[0]));

  s = IWORK.size();
  ff.write((char *)&s,   sizeof(s));
  if (s>0) ff.write((char *)IWORK.data(), s*sizeof(IWORK[0]));

  // pdecol structures
  ff.write((char *)&iounit_, sizeof(iounit_));
  ff.write((char *)&istart_, sizeof(istart_));
  ff.write((char *)&sizes_,  sizeof(sizes_));
  ff.write((char *)&gear0_,  sizeof(gear0_));
  ff.write((char *)&gear1_,  sizeof(gear1_));
  ff.write((char *)&gear9_,  sizeof(gear9_));
  ff.write((char *)&option_, sizeof(option_));

}

void
pdecol_solver::load_state(const std::string & fname){
  std::ifstream ff (fname.c_str(),std::ifstream::binary);
  int version;
  ff.read((char *)&version, sizeof(version));

  switch (version){
    case 1:
      ff.read((char *)&verbose, sizeof(verbose));
      ff.read((char *)&INDEX,   sizeof(INDEX));
      ff.read((char *)&NINT,    sizeof(NINT));
      ff.read((char *)&NPDE,    sizeof(NPDE));
      ff.read((char *)&KORD,    sizeof(KORD));
      ff.read((char *)&NCC,     sizeof(NCC));
      ff.read((char *)&MF,      sizeof(MF));
      ff.read((char *)&EPS,     sizeof(EPS));
      ff.read((char *)&t,       sizeof(t));
      ff.read((char *)&t0,      sizeof(t0));
      ff.read((char *)&mindt,   sizeof(mindt));

      int s;
      ff.read((char *)&s,   sizeof(s));
      XBKPT.resize(s);
      if (s>0) ff.read((char *)XBKPT.data(), s*sizeof(XBKPT[0]));

      ff.read((char *)&s,   sizeof(s));
      WORK.resize(s);
      if (s>0) ff.read((char *)WORK.data(), s*sizeof(WORK[0]));

      ff.read((char *)&s,   sizeof(s));
      IWORK.resize(s);
      if (s>0) ff.read((char *)IWORK.data(), s*sizeof(IWORK[0]));

      // pdecol structures
      ff.read((char *)&iounit_, sizeof(iounit_));
      ff.read((char *)&istart_, sizeof(istart_));
      ff.read((char *)&sizes_,  sizeof(sizes_));
      ff.read((char *)&gear0_,  sizeof(gear0_));
      ff.read((char *)&gear1_,  sizeof(gear1_));
      ff.read((char *)&gear9_,  sizeof(gear9_));
      ff.read((char *)&option_, sizeof(option_));

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
pdecol_solver::check_error(const int index) const {
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
}

void
pdecol_solver::print_index_info(const int index) const{
  struct timeval tt;
  gettimeofday(&tt, NULL);

  switch (index){
  case 0:
  case 2:
    std::cerr << tt.tv_sec << "." << std::setw(6) << std::setfill('0') << tt.tv_usec << ": ";
    std::cerr << "PDECOL: INDEX: " << index <<
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
    throw Err() << "bad or unknown INDEX setting in pdecol_solver::step";
  }
}
