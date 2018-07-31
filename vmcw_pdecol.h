#ifndef VMCW_PDECOL
#define VMCW_PDECOL
#include <vector>
#include <sstream>

/** C++ wrapper for PDECOL solver.
             V.Zavjalov, 07.2018
    See further information in pde_dp.f
**/


/********************************************************************/
// wrapper class
class pdecol_solver {

  /************************************/
  public:

  /// Constructor. Allocate memory, initialize PDECOL
  pdecol_solver(
    std::vector<double> &XSOL_, std::vector<double> &USOL_,
    double t0_, double dt_, double EPS_,
    int NPDE_, int NDERV_=2, int KORD_=4, int NCC_=2, int MF_=22, int verbose=1
  );

  // It is possible to change accuracy and method during calcultion.
  // In this case INDEX has to be set to 4;

  /// change EPS
  void ch_eps(double new_eps);

  /// change MF
  void ch_mf(int new_mf);

  /// get step size last used (sucsessfully)
  double get_dtused() const;

  /// get order last used (sucsessfully)
  int get_nqused() const;

  /// get number of steps fo far
  int get_nstep() const;

  // get number of residual evaluations fo far (RES calls)
  int get_nfe() const;

  // get number of matrix evaluations fo far (PSETIB calls)
  int get_nje() const;

  /// Do calculation until time t.
  // TODO: some more exotic calculations can be done (INDEX=2,3)
  // TODO: change verbosity using iounit_.LOUT
  int step(double t);

  /************************************/
  // Error class for exceptions
  public:
    class Err {
      std::ostringstream s;
      public:
        Err(){}
        Err(const Err & o) { s << o.s.str(); }
        template <typename T>
        Err & operator<<(const T & o){ s << o; return *this; }
        std::string str()  const { return s.str(); }
    };

  /************************************/
  private:

  double t0, dt; // initial time (used only on first call), min time step
  std::vector<double> &XSOL; // reference to x-vector
  std::vector<double> &USOL; // reference to solution vector

  int INDEX; // type of call -- result
  double EPS;                // Accuracy. Can be changed during calculations
  int NPTS, NINT, NPDE, KORD, NCC; // num.points, num.intervals, num.eq,
                                   // polynom.order (recommended 4), number of cont.cond (recommended 2)
                                   // Used only on first call.
  int MF;                    // The method flag. Can be changed during calculations
                             // 11,12,21 or 22
  int NDERV;                 // How many derivatives return into USOL

  std::vector<double> WORK;  // FLOATING POINT WORKING ARRAY FOR PDECOL.
  std::vector<int>   IWORK;  // INTEGER WORKING ARRAY FOR PDECOL.
  std::vector<double> SCTCH; // WORKING STORAGE ARRAY FOR VALUES.

  int verbose;
};

#endif
