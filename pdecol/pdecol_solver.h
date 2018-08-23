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
  ///   t0 -- Initial time.
  ///   dt -- Min time step (recommended 1e-10).
  ///   eps -- Accuracy. Can be changed during calculations (see ch_eps() below)
  ///   XBKPT -- strictly increasing array of piecewise polynomial breakpoints
  ///   NPDE  -- number of differential equations
  ///   KORD  -- polynom.order (recommended 4)
  ///   NCC   -- number of cont.cond (recommended 2)
  ///   MF    -- The method flag (11,12,21 or 22).
  ///            Can be changed during calculations (see ch_mf() below)
  ///   verbose -- verbosity level
  pdecol_solver(
    double t0, double dt, double eps, std::vector<double> & XBKPT,
    int NPDE, int KORD=4, int NCC=2, int MF=22, int verbose=1
  );

  /// Do calculation until time t.
  // TODO: some more exotic calculations can be done (INDEX=2,3)
  int step(double t);

  /// Get function values
  /// Arguments:
  ///   xsol  -- Cordinate values.
  ///   usol  -- storage to be filled.
  ///            Will be resized to xsol.size()*npde*(nderv+1) if needed.
  ///   nderv -- Number of derivatives to get.
  void values(std::vector<double> & xsol, std::vector<double> & usol, int nderv);

  /// Create the storage and get function values.
  std::vector<double> values(std::vector<double> & xsol, int nderv);

  /// Get value from the array returned by values().
  /// No range checking.
  /// Arguments:
  ///   NPTS -- size of X vector used in values()
  ///   ix   -- coordinate index, 0..NPTS-1
  ///   ie   -- equation index, 0..NPDE-1
  ///   id   -- derivative index, 0..nderv
  double get_value(const std::vector<double> & usol,
                   const int NPTS, const int ix,
                   const int ie, const int id) const{
    return usol[id*(NPDE*NPTS) + ix*NPDE + ie];}


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

  // get number of equations
  int get_npde() const {return NPDE;}


  // get coordinate span
  int get_xlen() const {return xmax-xmin;}

  // get min coordinate
  int get_xmin() const {return xmin;}

  // get min coordinate
  int get_xmax() const {return xmax;}

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

  void check_error(const int index); // throw Err if index!=0;

  int INDEX;  // type of call -- result
  double EPS; // Accuracy. Can be changed during calculations
  int NPDE;   // number of equations
  int KORD;   // polynom.order (used in the first call and in values())
  int MF;     // The method flag. Can be changed during calculations
              // 11,12,21 or 22
  double xmin, xmax; // min/max value of the coordinate

  std::vector<double> WORK;  // FLOATING POINT WORKING ARRAY FOR PDECOL.
  std::vector<int>   IWORK;  // INTEGER WORKING ARRAY FOR PDECOL.

  int verbose;
};

#endif
