#ifndef VMCW_PDECOL
#define VMCW_PDECOL
#include <vector>
#include <sstream>

/** C++ wrapper for PDECOL solver.
             V.Zavjalov, 07.2018
    See further information in pde_dp.f

    Additional functions should be provided:
      subroutine F(T,X,U,UX,UXX,FV,NPDE)
      subroutine BNDRY(T,X,U,UX,DBDU,DBDUX,DZDT,NPDE)
      subroutine UINIT(XI,UI,NPDEI)
      subroutine DERIVF(T,X,U,UX,UXX,DFDU,DFDUX,DFDUXX,NPDE)
**/

/// TODO: save/restore state

/********************************************************************/
// wrapper class
class pdecol_solver {

  /************************************/
  public:

  /// Constructor. Allocate memory, initialize PDECOL
  /// arguments:
  ///   XSOL
  ///   USOL
  ///   t0 -- Initial time.
  ///   dt -- Min time step.
  ///   eps -- Accuracy. Can be changed during calculations (see ch_eps() below)
  ///   NPDE  -- number of differential equations
  ///   NDERV -- How many derivatives return into USOL.
  ///   KORD  -- polynom.order (recommended 4)
  ///   NCC   -- number of cont.cond (recommended 2)
  ///   MF    -- The method flag (11,12,21 or 22).
  ///            Can be changed during calculations (see ch_mf() below)
  ///   verbose -- verbosity level
  pdecol_solver(
    double t0, double dt, double eps,
    int NPTS, int NPDE, int NDERV=2, int KORD=4, int NCC=2, int MF=22, int verbose=1
  );

  /// change EPS during calculation
  void ch_eps(double new_eps);

  /// change MF during calculation
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

  // get reference to the coordinate vector
  std::vector<double> & get_crd_vec() {return XSOL;}

  // get reference to the solution vector
  std::vector<double> & get_sol_vec() {return USOL;}

  /// Do calculation until time t.
  // TODO: some more exotic calculations can be done (INDEX=2,3)
  int step(double t);

  /// Write functions and derivatives vs. the coordinate.
  void write_profile(std::ostream &ss);
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
  std::vector<double> XSOL; // reference to x-vector
  std::vector<double> USOL; // reference to solution vector

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
