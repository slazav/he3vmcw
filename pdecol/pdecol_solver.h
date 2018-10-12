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

/********************************************************************/
// wrapper class
class pdecol_solver {

  /************************************/
  public:

  /// Constructor. Allocate memory, initialize PDECOL
  /// arguments:
  ///   t0 -- Initial time.
  ///   mindt -- Min time step (recommended 1e-10).
  ///   eps   -- Accuracy. Can be changed during calculations (see ch_eps() below)
  ///   XBKPT -- strictly increasing array of piecewise polynomial breakpoints
  ///   NPDE  -- number of differential equations
  ///   KORD  -- polynom.order (recommended 4)
  ///   NCC   -- number of cont.cond (recommended 2)
  ///   MF    -- The method flag (11,12,21 or 22).
  ///            Can be changed during calculations (see ch_mf() below)
  ///   verbose -- verbosity level
  pdecol_solver(
    double t0, double mindt, double eps, std::vector<double> & XBKPT,
    int NPDE, int KORD=4, int NCC=2, int MF=22, int verbose=1
  );

  /// Do calculation until time t.
  // TODO: some more exotic calculations can be done (INDEX=2,3)
  // Set exect=true to go to the time t exactly (should be done before
  // changing something in functions).
  int step(double t, bool exact=false);

  void restart();

  void reset_time();

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

  /// save state to a file
  void save_state(const std::string & fname);

  /// restore state from a file
  void load_state(const std::string & fname);

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
  double get_xlen() const {return *XBKPT.rbegin()-*XBKPT.begin();}

  // get min coordinate
  double get_xmin() const {return *XBKPT.begin();}

  // get min coordinate
  double get_xmax() const {return *XBKPT.rbegin();}

  // get current time
  double get_t() const {return t;}

  // get array of breakpoints used in the first call
  std::vector<double> & get_xmesh() {return XBKPT;}

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

  // throw Err if index!=0 && index!=-5 If index==-5, reset INDEX to 4;
  void check_error();
  void print_index_info() const; // print call info

  int INDEX;  // type of call -- result
  double EPS; // Accuracy. Can be changed during calculations
  int NINT;   // number of intervals -- for the first call
  int NPDE;   // number of equations
  int KORD;   // polynom.order (used in the first call and in values())
  int NCC;    // number of cont.cond (recommended 2) -- for the first call
  int MF;     // The method flag. Can be changed during calculations
              // 11,12,21 or 22

  double t; // current time
  double t0, mindt;  // initial time and min.time step -- for the first call
  std::vector<double> XBKPT;

  std::vector<double> WORK;  // FLOATING POINT WORKING ARRAY FOR PDECOL.
  std::vector<int>   IWORK;  // INTEGER WORKING ARRAY FOR PDECOL.

  int verbose;
};

#endif
