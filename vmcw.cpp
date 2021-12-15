#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "pdecol/pdecol_solver.h"
#include "pnm_writer.h"

/* Number of equations. 6 or 7. */
const int npde = 6;

/* How many derivatives to calculate. Can not be changed. */
const int nder = 3;


/******************************************************************/
// Error class for exceptions
class Err {
  std::ostringstream s;
  public:
    Err(){}
    Err(const Err & o) { s << o.s.str(); }
    template <typename T>
    Err & operator<<(const T & o){ s << o; return *this; }
    std::string str()  const { return s.str(); }
};

// A global variable with the program parameters.
// It should be global because it contains initial conditions
// and He3 parameters which should be accessed from
// UINIT_, F_, BNDRY_ functions called by the solver.
struct pars_t {

  /* Number of points. It can be changed at any time, but real change
  happenes when the solver is (re)started. In the code it can be used only
  to initialize the solver. In other places solver->get_npts() should be
  used (it corresponds to the actual number of points in the running solver). */
  int npts=257;

  /* Solver accuracy. Can be changed during calculation. There is in idea that
  if it is a power of 2 solver runs faster */
  double acc = pow(2,-20);

  /* Min time step (recommended 1e-10). It can be changes at any time,
  but real change happenes when the solver is (re)started. */
  double mindt = 1e-10;

  /* Current time [s]. Changes only before doing solver->step(). */
  double tcurr = 0;

  /* Time step [s]. Can be changed during calculation. */
  double tstep = 5e-3;

  /* End of current sweep/wait command [s]. Calculation is done until this time
  without reading new commands. Should be equal to "time" in the begining. */
  double tend = tcurr;

  /*****************************/
  // NMR frequency, magnetic fields.

  /* Gyromagnetic ratio. */
  const double gyro = 20378.0;

  /* NMR frequency [Hz]. The calculation is done at fixed frequency. Magnetic
  fiels should be changed to observe the resonance. Default value 1MHz
  corresponds to default magnetic field. */
  double f0 = 1e6;

  // Magnetic field.
  // Total field is: H = 2pi*f0/gyro + H0 + HT*t + HG*x + HQ*x^2. */
  //  H0 - uniform field measured from 2pi*f0/gyro [G],
  //  HT - sweep rate, dH/dt [G/s],
  //  HG - gradient, dH/dz [G/cm],
  //  HQ - quadratic term, d^2H/dz^2 [G/cm^2].
  double H0=0.0,  HT = 0.0,  HG = 0.1,  HQ = 0.0;
  double get_Wz(const double x, const double t) const {
    return 2*M_PI*f0 + gyro*(H0 + HG*x + HQ*x*x + HT*t); }
  double get_Wz(const double x) const { return get_Wz(x, tcurr); }

  // Radio-frequncy field.
  // Total field is: HR = (HR0 + HRT*t) * (1 + x/xRG + (x/xRQ)^2).
  // Gradient and quadratic terms change proportianally with the field.
  //  HR0 - Initial value [G],
  //  HRT - sweep rate, dHr/dt [G/s],
  //  HRGP - gradient profile [1/cm],
  //  HRQP - quadratic profile [1/cm^2].
  double HR0=1e-3, HRT=0.0, HRGP=0.0, HRQP=0.0;
  double get_Wr(const double x, const double t) const {
    return gyro*(HR0 + HRT*t) * (1.0 + x*HRGP + x*x*HRQP); }
  double get_Wr(const double x) const { return get_Wr(x, tcurr);}

  /*****************************/

  // Type of boundary conditions on the left (x<0) and right (x>0) sides.
  // 1 - open cell, 2 - no spin currents, 3 - constant functions (NPD wall).
  int bctype_l = 2;
  int bctype_r = 2;

  // Data for initial conditions.
  // Array with n*(npde+1) values. n is arbitrarary number,
  // values contain z coordinate and npde function components.
  // If n==1 then uniform i.c. with values from init_data are used;
  // If n==0 then default uniform i.c. (nz=1, mz=1, theta=acos(-1/4)) are used.
  std::vector<double> init_data;

  /*****************************/
  // He3 properties

  // t1 relaxation time, initial value [s] and sweep rate [s/s]
  double T10=1.0, T1T=0.0;
  double get_T1(const double t) const {return T10 + T1T*t;}
  double get_T1() const {return get_T1(tcurr);}

  // Leggett-Takagi relaxation time tau_f, initial value [s] and sweep rate [s/s]
  double TF0=1e-5, TFT=0.0;
  double get_TF(const double t) const {return TF0 + TFT*t;}
  double get_TF() const {return get_TF(tcurr);}

  // Spin-wave velocity, c_par, initial value [cm/s] and sweep rate [(cm/s)/s]
  double CP0=500.0, CPT=0.0;
  double get_CP(const double t) const {return CP0 + CPT*tcurr;}
  double get_CP() const {return get_CP(tcurr);}

  // Leggett frequency Omega_B, initial value [Hz] and sweep rate [Hz/s]
  double LF0=1e5, LFT=0.0;
  double get_LF(const double t) const {return LF0 + LFT*tcurr;}
  double get_LF() const {return get_LF(tcurr);}

  // Spin diffusion, initial value [cm^2/s] and sweep rate [(cm^2/s)/s]
  double DF0=0.1, DFT=0.0;
  double get_DF(const double t) const {return DF0 + DFT*tcurr;}
  double get_DF() const {return get_DF(tcurr);}

  // Reset sweeps (done before executing every new command)
  void reset_sweeps() {
    H0  = H0  + tcurr*HT;  HT=0.0;
    HR0 = HR0 + tcurr*HRT; HRT=0.0;
    TF0 = get_TF(); TFT=0.0;
    T10 = get_T1(); T1T=0.0;
    LF0 = get_LF(); LFT=0.0;
    CP0 = get_CP(); CPT=0.0;
    DF0 = get_DF(); DFT=0.0;
  }

  /*****************************/
  // sweep parameter offset, rate, name
  double sweep_par_o=0, sweep_par_r=0;
  std::string sweep_par_n;

  /*****************************/
  // cell length, cm
  double cell_len = 0.4;


  /*****************************/
  // files

  // prefix for data files, command file name without extension
  std::string pref;

  // mesh, prof, magn file counters (for default filenames)
  int cnt_mesh=0, cnt_prof=0, cnt_magn=0;

  /*****************************/

  /* PDECOL solver pointer */
  pdecol_solver *solver = NULL;

  /* Container for PNM writers. Key is the file name. */
  pnm_writers_t pnm_writers;

  /* Stream for magnetization writer */
  std::ostream *out_m = NULL;


} pp;



/******************************************************************/
// Set parameters for Leggett equations using main parameter structure.
// This functions are called from F and BNDRY
extern "C" {
  void set_bulk_pars_(double *t, double *x,
                 double *Wr, double *Wz, double *W0,
                 double *WB, double *Cpar, double *Diff,
                 double *Tf, double *T1){

    *W0 = 2*M_PI*pp.f0;
    *Wz = pp.get_Wz(*x,*t);
    *Wr = pp.get_Wr(*x,*t);

    *WB   = pp.get_LF(*t)*2*M_PI;
    *Cpar = pp.get_CP(*t);
    *Diff = pp.get_DF(*t);
    *Tf   = pp.get_TF(*t);
    *T1   = pp.get_T1(*t);
  }

  void set_bndry_pars_(double *t, double *x, double *Wz,
                 double *Cpar, double *Diff, int *IBN){

    *Wz = pp.get_Wz(*x,*t);
    *Cpar = pp.get_CP(*t);
    *Diff = pp.get_DF(*t);

    // type of boundary condition
    *IBN = (*x<0)? pp.bctype_l:pp.bctype_r;
  }
}

/******************************************************************/
/// UINIT function called from pdecol
extern "C" {
  void uinit_(double *x, double u[], int *n){

    // default
    u[0] = 0.0; // Mx
    u[1] = 0.0; // My
    u[2] = 1.0; // Mz
    if (npde == 7) {
      u[3] = 0.0; // nx
      u[4] = 0.0; // ny
      u[5] = 1.0; // nz
      u[6] = acos(-0.25); // th
    }
    else {
      u[3] = 0.0; // nx*th
      u[4] = 0.0; // ny*th
      u[5] = acos(-0.25); // nz*th
    }


    double w; // width
    double p; // -1..1
    int nn;

    nn = pp.init_data.size()/(npde+1);
    if (nn<1) return;
    if (nn==1 || *x<= pp.init_data[0]*pp.cell_len){
      for (int iu = 0; iu < npde; iu++)
          u[iu] = pp.init_data[iu+1];
      return;
    }
    if (*x>= pp.init_data[(nn-1)*(npde+1)]*pp.cell_len){
      for (int iu = 0; iu < npde; iu++)
        u[iu] = pp.init_data[(nn-1)*(npde+1)+iu+1];
      return;
    }
    for (int i=0; i<nn-1; i++){
      double x1 = pp.init_data[i*(npde+1)]*pp.cell_len;
      double x2 = pp.init_data[(i+1)*(npde+1)]*pp.cell_len;
      if (x1 < *x && *x <= x2){
        for (int iu = 0; iu < npde; iu++){
          double u1 = pp.init_data[i*(npde+1)+iu+1];
          double u2 = pp.init_data[(i+1)*(npde+1)+iu+1];
          u[iu] = u1 + (u2-u1)*(*x-x1)/(x2-x1);
        }
        return;
      }
    }
  }
}

/******************************************************************/

#ifdef HE3LIB
extern "C" {
#include <he3.h>
}
// set parameters using temperature, pressure and he3lib
void set_he3tp(double ttc, double p){
  double nu_b = he3_nu_b_(&ttc,&p);
  pp.CP0 = he3_cpar_(&ttc,&p) - pp.CPT*pp.tcurr;
  pp.LF0 = nu_b  - pp.LFT*pp.tcurr;
  pp.DF0 = he3_diff_perp_zz_(&ttc,&p,&pp.f0) - pp.DFT*pp.tcurr;
  double tr  = 1.2e-7/sqrt(1.0-ttc);
  pp.TF0   = 1.0/ (4.0*M_PI*M_PI * nu_b*nu_b * tr) - pp.TFT*pp.tcurr;
}
#endif


  // Create uniform mesh
  std::vector<double>
  make_uniform_mesh(const int N, const double cell_len){
    if (N < 2) throw Err() << "Too few points for making mesh: " << N;
    std::vector<double> x(N);
    double x0 = -cell_len/2.0;
    double dx = cell_len/(N-1);
    for (int i=0; i<N; i++) x[i] = x0 + i*dx;
    return x;
  }

  // Create adaptive mesh (running solver should exist)
  // xmesh_k -- non-uniform mesh step: cell_len/(NPTS-1) / (1+xmesh_k*aer_step(x)')
  //            0 means that mesh is uniform
  //            1 means that mesh is twice mere dense if RMS of function derivatives is 1
  // xmesh_acc -- non-uniform mesh accuracy
  std::vector<double>
  make_adaptive_mesh(pdecol_solver* solver, const int N, const double cell_len,
                     const double xmesh_k, const double xmesh_acc = 1e-10){
    if (!solver)
      throw Err() << "Running solver is needed for making adaptive mesh";

    if (N < 2)
      throw Err() << "Too few points for making mesh: " << N;

    //  start with a uniform mesh
    std::vector<double> x(N);
    double x0 = -cell_len/2.0;
    double dx = cell_len/(N-1);
    for (int i=0; i<N; i++) x[i] = x0 + i*dx;

    int maxk=100;
    for (int k=0; k<maxk; k++){
      // get function derivatives int the  mesh points
      std::vector<double> usol = solver->values(x, nder);
      // Calculate point weights: w(i) = 1./(1 + k <U[i]'>)
      // and total weight. Use only 1..N-1 points
      std::vector<double> w(N);
      double sum = 0;
      for (int i=1; i<N; i++){
        w[i]=0;
        for (int j=0; j<npde; j++){
          w[i] += pow( solver->get_value(usol, N, i, j, 1), 2);
        }
        w[i] = 1./(1 + xmesh_k*sqrt(w[i]));
        sum+=w[i];
      }

      // Modify the mesh by changing dx -> dx*sum/w
      // Find max point shift
      x[0] = x0;
      double sh=0;
      for (int i=1; i<N; i++){
        double x1 = x[i-1] + dx * w[i-1]/sum * (N-1);
        if (sh < abs(x[i] - x1)) sh = abs(x[i] - x1);
        x[i] = x1;
      }
      if (sh<xmesh_acc) break;
    }
    return x;
  }

/******************************************************************/
// Writing data

// Write profile to <fname>.
// Use solver mesh if N==0, or build a uniform one with N points (N>2)
// Now it saves everything in the same way for npde 6 or 7.
void
write_profile(pdecol_solver *solver, const std::string & fname, const int N) {

  if (!pp.solver)
     throw Err() << "Running solver is needed for writing profile";

  // use mesh from the solver (if N=0) or build a uniform mesh with N points:
  std::vector<double> xsol = (N==0)?
      pp.solver->get_xmesh() : make_uniform_mesh(N, pp.cell_len);
  std::vector<double> usol = pp.solver->values(xsol, nder);

  std::ofstream ss(fname.c_str());
  // print legend: # coord  U(0) U(1) ... U(0)' U(1)' ...
  ss << "# coord.     ";
  for (int d = 0; d<nder; d++){
    ss << "  ";
    for (int n = 0; n<npde; n++){
      ss << "U(" << n << ")";
      for (int i=0; i<d; i++) ss << "'";
      ss << "      ";
    }
  }
  ss << "\n" << std::scientific << std::setprecision(6);

  //print values
  for (int i=0; i< xsol.size(); i++){
    ss << xsol[i];
    for (int d = 0; d<nder; d++){
      ss << "  ";
      for (int n = 0; n<npde; n++)
        ss << " " << pp.solver->get_value(usol, xsol.size(), i, n, d);
    }
    ss << "\n";
  }
}

// write integral magnetization
void
write_magn(std::ostream & s){

  if (!pp.solver)
     throw Err() << "Running solver is needed for writing magnetization";

  // Use mesh from the solver
  std::vector<double> xsol = pp.solver->get_xmesh();
  std::vector<double> usol = pp.solver->values(xsol, nder);

  double smx=0, smy=0, smz=0, sz = 0;
  for (int i=0; i<xsol.size()-1; i++){
    double mx1 = usol[i*npde+0];
    double my1 = usol[i*npde+1];
    double mz1 = usol[i*npde+2];

    double mx2 = usol[(i+1)*npde+0];
    double my2 = usol[(i+1)*npde+1];
    double mz2 = usol[(i+1)*npde+2];
    double dz = xsol[i+1]-xsol[i];
    smx += dz*(mx1+mx2)/2;
    smy += dz*(my1+my2)/2;
    smz += dz*(mz1+mz2)/2;
    sz += dz;
  }
  smx/=sz;
  smy/=sz;
  smz/=sz;

  s << std::scientific << std::setprecision(6)
    << pp.tcurr << "  " << pp.H0 + pp.tcurr*pp.HT << "  "
    << smx << " " << smy << " " << smz << "\n";
}

void
write_pars(std::ostream & s){
  s << " T=" << pp.tcurr*1000 << " ms, "
    << "H0=" << pp.H0+pp.HT*pp.tcurr << " G, "
    << "HR=" << 1e3*(pp.HR0+pp.HRT*pp.tcurr) << " mOe, "
    << "LF=" << 1e-3*pp.get_LF() << " kHz, "
    << "CP=" << pp.get_CP() << " cm/s, "
    << "DF=" << pp.get_DF() << " cm^2/s, "
    << "TF=" << 1e6*pp.get_TF() << " mks, "
    << "T1=" << pp.get_T1() << " s, "
    << "\n";
}


/******************************************************************/
// Some init_data-related functions.

// make uniform initial conditions in init_data
void
init_data_uniform(const double mx, const double my, const double mz,
                  const double nx, const double ny, const double nz,
                  const double th = acos(-0.25)) {
  pp.init_data.resize(npde+1);
  pp.init_data[0]=0;
  pp.init_data[1]=mx;
  pp.init_data[2]=my;
  pp.init_data[3]=mz;
  if (npde == 7) {
    pp.init_data[4]=nx;
    pp.init_data[5]=ny;
    pp.init_data[6]=nz;
    pp.init_data[7]=th;
  }
  else {
    pp.init_data[4]=nx*th;
    pp.init_data[5]=ny*th;
    pp.init_data[6]=nz*th;
  }
}

// make initial conditions with a simple soliton in init_data
// width >0 does not work yet (vectors do not rotate!)
void
init_data_soliton(double w, // soliton width
                  double mx1, double mx2, // mx on the left and write side
                  double my1, double my2, // my on the left and write side
                  double mz1, double mz2, // mz on the left and write side
                  double nx1, double nx2, // nx on the left and write side
                  double ny1, double ny2, // ny on the left and write side
                  double nz1, double nz2, // nz on the left and write side
                  double th1, double th2  // th on the left and write side
                 ) {
  if (w<0){
    w=-w;
    std::swap(mx1,mx2); std::swap(my1,my2); std::swap(mz1,mz2);
    std::swap(nx1,nx2); std::swap(ny1,ny2); std::swap(nz1,nz2); std::swap(th1,th2);
  }
  pp.init_data.resize(2*(npde+1));
  pp.init_data[0]=-w/2; pp.init_data[npde+1+0]=+w/2;
  pp.init_data[1]=mx1;  pp.init_data[npde+1+1]=mx2;
  pp.init_data[2]=my1;  pp.init_data[npde+1+2]=my2;
  pp.init_data[3]=mz1;  pp.init_data[npde+1+3]=mz2;
  if (npde==7){
    pp.init_data[4]=nx1; pp.init_data[npde+1+4]=nx2;
    pp.init_data[5]=ny1; pp.init_data[npde+1+5]=ny2;
    pp.init_data[6]=nz1; pp.init_data[npde+1+6]=nz2;
    pp.init_data[7]=th1; pp.init_data[npde+1+7]=th2;
  }
  else {
    pp.init_data[4]=nx1*th1; pp.init_data[npde+1+4]=nx2*th2;
    pp.init_data[5]=ny1*th1; pp.init_data[npde+1+5]=ny2*th2;
    pp.init_data[6]=nz1*th1; pp.init_data[npde+1+6]=nz2*th2;
  }
}

// set HPD initial condition (RF-field, freq shift, Leggett-Takagi relaxation is used).
// sn=+1/-1 and st=+1/-1 are ny and theta sign.
void
init_data_hpd(int sn=1, int st=1){

  double Wr, Wz, W0, WB, Cpar, Diff, Tf, T1;

  sn = (sn>0)? +1:-1;
  st = (st>0)? +1:-1;

  pp.init_data.resize(pp.npts*(npde+1));
  for (int i=0; i<pp.npts; i++){
    double x  = i/(pp.npts+1.0) - 0.5;
    double xl = x*pp.cell_len;
    // get local parameters
    set_bulk_pars_(&pp.tcurr, &xl, &Wr, &Wz, &W0,&WB, &Cpar, &Diff, &Tf, &T1);

    double h = Wr/W0;
    double d = -(Wz-W0)/W0;
    double b = (WB/W0)*(WB/W0);
    double wt = Tf * W0;
//    if (d<0) throw Err() << "init_data_hpd: d<0";

    double th = (st>0? 1:-1) * acos(-0.25 - 15.0/16 * d/b - sn*h/sqrt(15.0));

    double nx = -st*sqrt(15.0)/4.0 * d*d/h/b/wt;
    double nz = -st*sqrt(15.0/16.0) * d/b/wt;
    double ny = (sn>0? 1:-1) * sqrt(1-nx*nx-nz*nz);

    double wz = cos(th);
    double wx = ny*sin(th);
    double wy = b/h*nz*nz*wt;

    pp.init_data[(npde+1)*i+0] = x;
    pp.init_data[(npde+1)*i+1] = wx + h;
    pp.init_data[(npde+1)*i+2] = wy;
    pp.init_data[(npde+1)*i+3] = wz - d;
    if (npde == 7){
      pp.init_data[(npde+1)*i+4] = nx;
      pp.init_data[(npde+1)*i+5] = ny;
      pp.init_data[(npde+1)*i+6] = nz;
      pp.init_data[(npde+1)*i+7] = th;
    }
    else {
      pp.init_data[(npde+1)*i+4] = nx*th;
      pp.init_data[(npde+1)*i+5] = ny*th;
      pp.init_data[(npde+1)*i+6] = nz*th;
    }
  }
}

// set NPD initial condition (RF-field, freq shift).
// sn=+1/-1 and st=+1/-1 are ny and theta sign.
// Spin diff is not used, could be problems near resonance.
void
init_data_npd(int sn=1, int st=1){

  double Wr, Wz, W0, WB, Cpar, Diff, Tf, T1;

  sn = (sn>0)? +1:-1;
  st = (st>0)? +1:-1;

  pp.init_data.resize(pp.npts*(npde+1));
  for (int i=0; i<pp.npts; i++){
    double x  = i/(pp.npts+1.0) - 0.5;
    double xl = x*pp.cell_len;
    // get local parameters
    set_bulk_pars_(&pp.tcurr, &xl, &Wr, &Wz, &W0,&WB, &Cpar, &Diff, &Tf, &T1);

    double th = (st>0? 1:-1) * acos(-0.25);
    double wx = -Wr/sqrt(pow(Wr,2) + pow(Wz-W0,2));
    double wz = -(Wz-W0)/sqrt(pow(Wr,2) + pow(Wz-W0,2));
    if (wz<0) {wz=-wz; wx=-wx;}

    double nz = (sn>0? 1:-1) * sqrt((wz+0.25)/1.25);
    double nx = nz*wx /(wz + 1.0);
    double ny = sqrt(3.0/5.0)*wx /(wz + 1.0);

    // safety:
    double nn = sqrt(nx*nx+ny*ny+nz*nz);
    nx = nx/nn; ny=ny/nn; nz=nz/nn;

    std::cerr << "W> " << wx << " " << wz <<"\n";
    std::cerr << "N> " << nx << " " << nz <<"\n";

    // s/W0 = w - z + Wz/W0 z + Wr/W0 x
    pp.init_data[(npde+1)*i+0] = x;
    pp.init_data[(npde+1)*i+1] = wx + Wr/W0;
    pp.init_data[(npde+1)*i+2] = 0;
    pp.init_data[(npde+1)*i+3] = wz + (Wz-W0)/W0;
    if (npde == 7){
      pp.init_data[(npde+1)*i+4] = nx;
      pp.init_data[(npde+1)*i+5] = ny;
      pp.init_data[(npde+1)*i+6] = nz;
      pp.init_data[(npde+1)*i+7] = th;
    }
    else {
      pp.init_data[(npde+1)*i+4] = nx*th;
      pp.init_data[(npde+1)*i+5] = ny*th;
      pp.init_data[(npde+1)*i+6] = nz*th;
    }
  }
}

// Save current profile to init_data.
void
init_data_save(pdecol_solver *solver) {

  if (!solver)
     throw Err() << "Running solver is needed for writing magnetization";

  // get values using mesh from the running solver
  std::vector<double> xsol = pp.solver->get_xmesh();
  std::vector<double> usol = solver->values(xsol, 0); // no derivatives!

  pp.init_data = std::vector<double>(xsol.size()*(npde+1));
  // fill init_data array:
  for (int i=0; i< xsol.size(); i++){
    pp.init_data[i*(npde+1)] = xsol[i]/solver->get_xlen();
    for (int n = 0; n<npde; n++)
      pp.init_data[i*(npde+1)+n+1] =
        solver->get_value(usol, xsol.size(), i, n, 0);
  }
}


/******************************************************************/
// Helpers for command parsing


// check number of arguments for a command, throw Err if it is wrong
void
check_nargs(int nargs, int n1, int n2=-1){
  if (n2<0 && nargs!=n1)
    throw Err() << "command requires "  << n1 << " arguments";
  if (n2>=0 && (nargs < n1 || nargs > n2))
    throw Err() << "command requires "  << n1 << " to " << n2 << " arguments";
}

// Read an argument of arbitrary type from string. Throw Err in case of
// error.
template <typename T> T
get_arg(const std::string &str){
  std::istringstream ss(str);
  T v; ss >> v;
  if (!ss.eof() || ss.fail())
    throw Err() << "can't parse argument \"" << str << "\"";
  return v;
}

// check that there is one argument in the list and read it.
// (for set, step commands)
template <typename T> T
get_one_arg(const std::vector<std::string> & args){
  check_nargs(args.size(), 1);
  return get_arg<T>(args[0]);
}

// Process SWEEP command.
// Modify P0, PT, change global variable tend.
template <typename T>
void
cmd_sweep(const char *name, const std::vector<std::string> & args, T *P0, T *PT, T factor=1){
  check_nargs(args.size(), 2);
  double VD = get_arg<double>(args[0]); // destination
  double R = get_arg<double>(args[1]);  // rate
  double VO = *P0/factor;                // old value

  pp.sweep_par_o = VO;
  pp.sweep_par_r = R;
  pp.sweep_par_n = name;

  int steps = abs(rint((VD-VO)/R/pp.tstep));
  if (steps==0) throw Err() << "zero steps for sweep";
  *PT = (VD-VO)/(steps*pp.tstep) * factor;
  *P0 -= pp.tcurr*(*PT);
  pp.tend  = pp.tcurr + steps*pp.tstep;
}

/******************************************************************/
/// Read one or more commends from a stream.
/// Return 1 if file is finished, 0 if more calculations are needed.
int
read_cmd(std::istream &in_c, std::ostream & out_c){
  pp.reset_sweeps();

  // Read input string line by line
  // Stop reading is some command increase tend, then
  // calculation is needed.
  std::string line;
  while (pp.tend <= pp.tcurr){

    // return 1 at the end of command file
    if (!getline(in_c, line)){
      out_c << "End of command file, stop calculations\n";
      return 1;
    }

    // remove comments
    size_t c = line.find('#');
    if (c != std::string::npos){
      line = line.substr(0, c);
    }

    // split the line into command and arguments
    std::istringstream in(line);
    std::string cmd;
    in >> cmd;
    if (!in) continue; // empty command

    std::vector<std::string> args;
    while (1) {
      std::string a;
      in >> a;
      if (!in) break;
      args.push_back(a);
    }
    int narg = args.size();

    // logging to out_c
    out_c << "Command: " << line << "\n";
    std::cerr << line << "\n";

    /*******************************************************/
    // solver

    // (Re)start the solver.
    if (cmd == "start") {
      check_nargs(narg, 0);
      // stop solver if it already exists
      if (pp.solver) delete pp.solver;

      // destroy all writers
      pp.pnm_writers.clear();

      // initialize new solver
      pp.tend=pp.tcurr=0;
      std::vector<double> xbrpt = make_uniform_mesh(pp.npts, pp.cell_len);
      pp.solver = new pdecol_solver(pp.tcurr, pp.mindt, pp.acc, xbrpt, npde);
      if (!pp.solver) throw Err() << "can't start solver";
      continue;
    }

    // stop the solver
    if (cmd == "stop") {
      check_nargs(narg, 0);
      if (!pp.solver) throw Err() << "solver is not running";
      if (pp.solver) delete pp.solver;
      pp.solver = NULL;
      continue;
    }

    // exit the program
    if (cmd == "exit") {
      check_nargs(narg, 0);
      if (pp.solver) delete pp.solver;
      pp.pnm_writers.clear();
      pp.solver=NULL;
      return 1;
    }

    // save state to a file
    if (cmd == "save_state") {
      check_nargs(narg, 1);
      if (!pp.solver) throw Err() << "solver is not running";
      pp.solver->save_state(args[0]);
      continue;
    }

    // read state from a file
    if (cmd == "load_state") {
      check_nargs(narg, 1);
      if (!pp.solver){
        // create some solver (parameters are not important)
        std::vector<double> xbrpt(pp.npts,0.0);
        pp.solver = new pdecol_solver(pp.tcurr, pp.mindt, pp.acc, xbrpt, npde);
      }
      pp.solver->load_state(args[0]);
      pp.tcurr = pp.tend = pp.solver->get_t();
      continue;
    }

    // find equilibrium
    if (cmd == "find_eq") {
      check_nargs(narg, 0);
      if (!pp.solver) throw Err() << "solver is not running";
      std::vector<double> xsol = make_uniform_mesh(pp.npts,pp.cell_len);
      auto usol = pp.solver->find_eq(xsol);
      pp.init_data = std::vector<double>(xsol.size()*(npde+1));
      for (int i=0; i< xsol.size(); i++){
        pp.init_data[i*(npde+1)] = xsol[i]/pp.solver->get_xlen();
        for (int n = 0; n<npde; n++)
          pp.init_data[i*(npde+1)+n+1] = usol[i*npde + n];
      }
      pp.solver->restart();
      continue;
    }


    if (cmd == "reset_time") {
      check_nargs(narg, 0);
      if (!pp.solver) throw Err() << "solver is not running";
      pp.solver->reset_time();
      pp.tcurr = pp.tend = 0;
      continue;
    }

    // Do calculations for some time.
    if (cmd == "wait") {
      if (!pp.solver) throw Err() << "solver is not running";
      pp.sweep_par_o = 0.0;
      pp.sweep_par_r = 1.0;
      pp.sweep_par_n = 'T';
      pp.tend = pp.tcurr + get_one_arg<double>(args);
      continue;
    }

    // Change solver accuracy. If solver is not running,
    // the value will be used after start.
    if (cmd == "acc") {
      pp.acc = get_one_arg<double>(args);
      if (pp.solver) pp.solver->ch_eps(pp.acc);
      continue;
    }

    // Change accuracy (power of two). If solver is not running,
    // the value will be used after start.
    if (cmd == "acc2") {
      pp.acc=pow(2, -get_one_arg<double>(args));
      if (pp.solver) pp.solver->ch_eps(pp.acc);
      continue;
    }

    // Change min. time step (recommended 1e-10). It can be changes at any time,
    // but real change happenes when the solver is (re)started.
    if (cmd == "mindt") { pp.mindt = get_one_arg<double>(args); continue; }

    // Change number of points.
    if (cmd == "npts") { pp.npts = get_one_arg<int>(args); continue; }

    // Change time step. Can be changed during calculation.
    if (cmd == "tstep") { pp.tstep = get_one_arg<double>(args); continue; }

    // Change cell length. It can be changes at any time,
    // but real change happenes when the solver is (re)started.
    if (cmd == "cell_len")  { pp.cell_len = get_one_arg<double>(args); continue;}

    //if (cmd == "aer_len")   { pp.aer_len = get_one_arg<double>(args); continue;}
    //if (cmd == "aer_cnt")   { pp.aer_cnt = get_one_arg<double>(args); continue;}
    //if (cmd == "aer_trw")   { pp.aer_trw = get_one_arg<double>(args); continue;}

    /*******************************************************/
    // boubdary and initial conditions

    if (cmd == "bcond_type"){
      pp.bctype_l = pp.bctype_r = get_one_arg<int>(args); continue;}
    if (cmd == "bcond_type_l"){
      pp.bctype_l = get_one_arg<int>(args); continue;}
    if (cmd == "bcond_type_r"){
      pp.bctype_r = get_one_arg<int>(args); continue;}


    // Set "hpd" i.c. with ny=-1 or ny=+1 (default).
    // He3 parameters, cell size and field profile should be set before.
    // At each point equilibrium value is set (assuming system is uniform locally).
    // Spin currents are ignored, field profile, Leggett frequency and
    // Leggett-takagi relaxation are used.
    if (cmd == "set_icond_hpd") {
      check_nargs(narg, 0, 1);
      int ny = narg>0 ? get_arg<int>(args[0]) : 1;
      ny = ny>=0? 1:-1;
      init_data_hpd(ny);
      continue;
    }

    // Set "npd" i.c. with ny=-1 or ny=+1 (default).
    // He3 parameters, cell size and field profile should be set before.
    // At each point equilibrium value is set (assuming system is uniform locally).
    // Spin currents are ignored, field profile and Leggett frequency are used.
    if (cmd == "set_icond_npd") {
      check_nargs(narg, 0, 1);
      int ny = narg>0 ? get_arg<int>(args[0]) : 1;
      ny = ny>=0? 1:-1;
      init_data_npd(ny);
      continue;
    }

    // Simple HPD initial condition
    if (cmd == "set_icond_hpd_simple") {
      check_nargs(narg, 0, 1);
      int ny = narg>0 ? get_arg<int>(args[0]) : 1;
      ny = ny>=0? 1:-1;
      double th = acos(-0.25);
      init_data_uniform(ny*sin(th),0,cos(th), 0,ny,0, th);
      continue;
    }

    // set uniform i.c. with nz=-1 or nz=+1 (default)
    if (cmd == "set_icond_uniform" ||
        cmd == "set_icond_npd_simple") {
      check_nargs(narg, 0, 1);
      int nz = narg>0 ? get_arg<int>(args[0]) : 1;
      nz = nz>=0? 1:-1;
      init_data_uniform(0,0,1, 0,0,nz, acos(-0.25));
      continue;
    }

    // set i.c witn a simple n-soliton <width>
    // width >0 does not work yet
    if (cmd == "set_icond_nsol") {
      check_nargs(narg, 1, 1);
      double w = get_one_arg<double>(args);
      double th = acos(-0.25);
      init_data_soliton(w, 0,0, 0,0, 1,1,
                           0,0, 0,0, 1,1, -th,th);
      continue;
    }

    // set i.c witn a simple t-soliton <width>
    if (cmd == "set_icond_tsol") {
      check_nargs(narg, 1, 1);
      double w = get_one_arg<double>(args);
      double th = acos(-0.25);
      init_data_soliton(w, 0,0, 0,0, 1,1,
                           0,0, 0,0, 1,1, -th,th);
      continue;
    }

    // set the initial condtions
    if (cmd == "init") {
      double th0 = acos(-0.25);
      check_nargs(narg, 1, 2);
      std::string type = args[0];
      if (type == "NPD"){
        init_data_uniform(0,0,1, 0,0,1, th0);
        break;
      }
      if (type == "NPD-"){
        init_data_uniform(0,0,1, 0,0,1, -th0);
        break;
      }
      if (type == "th_soliton"){
        double w = (narg<2)? 0.01 : get_arg<double>(args[1]);
        init_data_soliton(w, 0,0, 0,0, 1,1,
                             0,0, 0,0, 1,1, th0, 2*M_PI-th0);
        break;
      }

    }


    // deform the current solution and restart the solver
    if (cmd == "deform") {
      if (!pp.solver) throw Err() << "solver is not running";
      check_nargs(narg, 1, 4);
      std::string type = args[0];
      init_data_save(pp.solver); // save current profile to the init data

      int N = pp.init_data.size()/(npde+1);
      // Be careful with first and last point!
      // Initial conditions should be compatable with boundary conditions.
      // It is better not to modify these points
      for (int i=0; i<N; i++){
        double  x = pp.init_data[i*(npde+1) + 0];
        double mx = pp.init_data[i*(npde+1) + 1];
        double my = pp.init_data[i*(npde+1) + 2];
        double mz = pp.init_data[i*(npde+1) + 3];
        double mm = sqrt(mx*mx + my*my + mz*mz);
        double am = atan2(my,mx);
        double bm = acos(mz/mm);

        double nx,ny,nz,th;
        nx = pp.init_data[i*(npde+1) + 4];
        ny = pp.init_data[i*(npde+1) + 5];
        nz = pp.init_data[i*(npde+1) + 6];
        th = sqrt(nx*nx + ny*ny + nz*nz);
        nx/=th; ny/=th; nz/=th;
        double an = atan2(ny,nx);
        double bn = acos(nz);
        if (npde==7) th = pp.init_data[i*(npde+1) + 7];

        // trivial
        if (type == "0") break;

        // rotate by PI
        if (type == "half_turn") {
          an+=M_PI; am+=M_PI;}

        // constant rotation aroung z axis (number of periods - parameter n)
        if (type == "rotation") {
          check_nargs(narg, 1, 2);
          int n = (narg<2)? 1 : get_arg<int>(args[1]);
          an -= 2*n*M_PI*(x+0.5);
          am -= 2*n*M_PI*(x+0.5);
        }

        // 2pi soliton with width w (orientation depends on w sign)
        if (type == "2pi_soliton") {
          check_nargs(narg, 1, 2);
          double w = (narg<2)? 0.1 : get_arg<double>(args[1]);
          an += 4*atan(exp(x/w));
          am += 4*atan(exp(x/w));
        }

        // inverse direction of n)
        if (type == "inverse_n") {
          an = -an;
          bn = M_PI - bn;
        }

        // theta soliton: HPD -> NPD- -> NPD+ -> HPD
        if (type == "th_soliton") {
          check_nargs(narg, 1, 2);
          double w = (narg<2)? 0.01 : get_arg<double>(args[1]);
          if (x/w>=-1.5 && x/w<-0.5){
            double k = x/w+1.5; // 0..1
            bn = bn*(1-k) + M_PI*k;
            bm = bm*(1-k);
          }
          if (x/w>=-0.5 && x/w<+0.5){
            double k = x/w+0.5; // 0..1
            bn = M_PI;
            th = th*(1-k) + (2*M_PI-th)*k;
          }
          if (x/w>=+0.5 && x/w<+1.5){
            double k = x/w-0.5; // 0..1
            th = 2*M_PI-th;
            bn = M_PI*(1-k)+k*bn;
            bm = bm*k;
            an=-an;
          }
          if (x/w>=+1.5){
            th = 2*M_PI-th;
            an=-an;
          }
          if (x/w>=-0.5 && x/w<0.5) bm = 0;
        }


        // theta soliton: HPD -> NPD- -> NPD+ -> HPD
        // Same, but with sharp theta step
        if (type == "th_soliton1") {
          check_nargs(narg, 1, 4);
          double w  = (narg<2)? 0.01 : get_arg<double>(args[1]);
          int   sn = (narg<3)? -1 : get_arg<int>(args[2]);
          int   t  = (narg<4)? 1  : get_arg<int>(args[3]);

          double wB   = pp.get_LF()*2*M_PI;
          double Cpar = pp.get_CP();
          // xiD = 13/24 K1/gD = 65/64 c_par^2/wB^2
          // here I use approximation K1 = K/4
          double xiD = sqrt(65.0/64.0) * Cpar/wB;
          double X = x*pp.cell_len; // cm

          // exact solution of theta-soliton:
          if (t == 1)
            // thL -> 2pi-thL
            th = M_PI + 2.0*atan(
              (1.0+cos(th))/sin(th) * tanh(sqrt(65.0/64.0)*X/2.0/xiD));
          else
            // thL -> -thL
            th = 2.0*atan(
              sin(th)/(1.0+cos(th)) * tanh(-sqrt(65.0/64.0)*X/2.0/xiD));


          // n-vector - just a linear changes to nz=0
          if (x/w>=-1 && x/w<0){
            double k = x/w+1; // 0..1
            bn = bn*(1-k) + (sn<0 ?M_PI*k:0);
            bm = bm*(1-k);
          }
          if (x/w>=0 && x/w<1){
            double k = (x/w); // 0..1
            an += M_PI;
            bn = (sn<0? M_PI*(1-k):0) + bn*k;
            bm = bm*k;
          }
          if (x/w>=1){
            an += M_PI;
          }

        }

        // theta soliton: HPD -> NPD- -> 0 -> NPD- -> HPD
        // Sharp theta step.
        if (type == "th_soliton2") {
          check_nargs(narg, 1, 2);
          double w = (narg<2)? 0.01 : get_arg<double>(args[1]);
          if (x/w>=-1 && x/w<0){
            double k = x/w+1; // 0..1
            bn = bn*(1-k) + M_PI*k;
            bm = bm*(1-k);
          }
          if (x/w>=0 && x/w<1){
            double k = (x/w); // 0..1
            bn = M_PI*(1-k) + bn*k;
            bm = bm*k;
            th = 2*M_PI-th;
          }
          if (x/w>=1){
            th = 2*M_PI-th;
          }
        }

        nx = sin(bn)*cos(an);
        ny = sin(bn)*sin(an);
        nz = cos(bn);

        mx = mm*sin(bm)*cos(am);
        my = mm*sin(bm)*sin(am);
        mz = mm*cos(bm);

        pp.init_data[i*(npde+1) + 0] =  x;
        pp.init_data[i*(npde+1) + 1] = mx;
        pp.init_data[i*(npde+1) + 2] = my;
        pp.init_data[i*(npde+1) + 3] = mz;
        if (npde==7) {
          pp.init_data[i*(npde+1) + 4] = nx;
          pp.init_data[i*(npde+1) + 5] = ny;
          pp.init_data[i*(npde+1) + 6] = nz;
          pp.init_data[i*(npde+1) + 7] = th;
        }
        else {
          pp.init_data[i*(npde+1) + 4] = nx*th;
          pp.init_data[i*(npde+1) + 5] = ny*th;
          pp.init_data[i*(npde+1) + 6] = nz*th;
        }
     }
      pp.solver->restart();
      continue;
    }

    if (cmd == "adaptive_mesh") {
      if (!pp.solver) throw Err() << "solver is not running";
      check_nargs(narg, 0,2);
      if (narg>=1)  pp.npts = get_arg<int>(args[0]);
      double xmesh_k = (narg<2)? 1 : get_arg<double>(args[1]);

      // destroy all writers
      pp.pnm_writers.clear();

      // initialize new mesh and fill init_data
      std::vector<double> xbrpt = make_adaptive_mesh(pp.solver, pp.npts, pp.cell_len, xmesh_k);

      std::vector<double> usol = pp.solver->values(xbrpt, 0); // no derivatives!
      pp.init_data = std::vector<double>(xbrpt.size()*(npde+1));
      //print values
      for (int i=0; i< xbrpt.size(); i++){
        pp.init_data[i*(npde+1)] = xbrpt[i]/pp.solver->get_xlen();
        for (int n = 0; n<npde; n++)
          pp.init_data[i*(npde+1)+n+1] = pp.solver->get_value(usol, xbrpt.size(), i, n, 0);
      }

      pp.solver->restart();
      continue;

    }

    /*******************************************************/
    // write profile

    // write function profiles to a file
    if (cmd == "write_profile") {
      check_nargs(narg, 0,2);
      if (!pp.solver) throw Err() << "solver is not running";

      std::string name;
      if (narg>0 && args[0] != "-") {
        name = args[0];
      }
      else{
        std::ostringstream ss;
        ss << pp.pref << ".prof" << pp.cnt_prof << ".dat";
        name = ss.str();
        pp.cnt_prof++;
      }
      int N = (narg>1) ? atoi(args[1].c_str()) : 0;
      write_profile(pp.solver, name, N);
      continue;
    }

    // write mesh
    if (cmd == "write_mesh") {
      check_nargs(narg, 0,1);
      if (!pp.solver) throw Err() << "solver is not running";

      std::string name = narg>0 ? args[0]:"";
      if (name =="") {
        std::ostringstream ss;
        ss << pp.pref << ".mesh" << pp.cnt_mesh << ".dat";
        name = ss.str();
        pp.cnt_mesh++;
      }
      auto xbrpt = pp.solver->get_xmesh();
      std::ofstream ff(name);
      ff << "# N X\n" << std::scientific << std::setprecision(6);
      for (int i=0; i<xbrpt.size(); i++)
        ff << i << " " << xbrpt[i] << "\n";

      continue;
    }

    /*******************************************************/
    // pnm writer

    // Initialize a pnm_writer. Solver should be startded.
    if (cmd == "pnm_start") {
      check_nargs(narg, 0,1);
      std::string name = narg>0 ? args[0]: pp.pref + ".pic.pnm";
      if (!pp.solver)
        throw Err() << "can't start pnm_writer if solver is not running";
      if (!pp.pnm_writers.add(name, pp.solver))
        throw Err() << "can't open create pnm_writer";
      continue;
    }

    // pnm_legend: start draw a legend
    if (cmd == "pnm_legend") {
      check_nargs(narg, 0,1);
      std::string name = narg>0 ? args[0]: pp.pref + ".pic.pnm";
      if (!pp.pnm_writers.legend(name, 50, 100))
        throw Err() << "no such writer";
      continue;
    }

    // pnm_hline: draw a horizontal line
    if (cmd == "pnm_hline") {
      check_nargs(narg, 0,1);
      std::string name = narg>0 ? args[0]: pp.pref + ".pic.pnm";
      if (!pp.pnm_writers.hline(name))
        throw Err() << "no such writer";
      continue;
    }

    // stop a pnm_writer
    if (cmd == "pnm_stop") {
      check_nargs(narg, 0,1);
      std::string name = narg>0 ? args[0]: pp.pref + ".pic.pnm";
      if (!pp.pnm_writers.del(name))
        throw Err() << "no such writer";
      continue;
    }


    // Start recording total magnetization to a file.
    if (cmd == "magn_start") {
      check_nargs(narg, 0,1);
      std::string name = narg>0 ? args[0]: "";
      if (name=="") {
        std::ostringstream ss;
        ss << pp.pref << ".magn" << pp.cnt_magn << ".dat";
        pp.cnt_magn++;
        name = ss.str();
      }

      if (pp.out_m) delete pp.out_m;
      pp.out_m = new std::ofstream(name);
      if (!pp.out_m || !pp.out_m->good())
        throw Err() << "Can't open file: " << name;
      *pp.out_m << "# Integral magnetization log: T, LP, Mx, Mx, Mz\n";
      continue;
    }

    // Start recording total magnetization.
    if (cmd == "pnm_stop") {
      if (pp.out_m) {
        delete pp.out_m;
        pp.out_m = NULL;
      }
      continue;
    }

    /*******************************************************/
    // set NMR frequency
    if (cmd == "set_freq") {
      pp.f0 = get_one_arg<double>(args); continue; }

    // Set uniform field [G, from larmor].
    if (cmd == "set_field") {
      pp.H0 = get_one_arg<double>(args); continue; }

    // Uniform field step [G].
    if (cmd == "step_field") {
      pp.H0 += get_one_arg<double>(args); continue; }

    // Set field gradient [G/cm].
    if (cmd == "set_field_grad") {
      pp.HG = get_one_arg<double>(args); continue; }

    // Set field quadratic term [G/cm^2].
    if (cmd == "set_field_quad") {
      pp.HQ = get_one_arg<double>(args); continue; }

    // Sweep uniform field: destination [G], rate [G/s].
    if (cmd == "sweep_field") {
      cmd_sweep("B(G)", args, &pp.H0, &pp.HT); continue; }

    // Set uniform field in frequency shift units [Hz from NMR freq].
    if (cmd == "set_field_hz") {
      pp.H0 = get_one_arg<double>(args) * 2*M_PI/pp.gyro; continue; }

    // Uniform field step [Hz].
    if (cmd == "step_field_hz") {
      pp.H0 += get_one_arg<double>(args) * 2*M_PI/pp.gyro; continue; }

    // Sweep uniform field: destination [Hz], rate [Hz/s].
    if (cmd == "sweep_field_hz") {
      cmd_sweep("B(Hz)", args, &pp.H0, &pp.HT, 2*M_PI/pp.gyro); continue; }

    // Set uniform field in Larmor position units [cm].
    // Gradient term is used to convert field to cm. Quadratic term is not used.
    // Lower Larmor positon means higher field.
    if (cmd == "set_field_cm") {
      if (pp.HG == 0.0) throw Err() << "can't set Larmor position if "
                                    "field gradient is zero";
      pp.H0 = -get_one_arg<double>(args)*pp.HG; continue; }

    // Uniform field step [cm].
    if (cmd == "step_field_cm") {
      if (pp.HG == 0.0) throw Err() << "can't set Larmor position if "
                                    "field gradient is zero";
      pp.H0 -= get_one_arg<double>(args)*pp.HG; continue; }

    // Sweep uniform field: destination [cm], rate [cm/s].
    if (cmd == "sweep_field_cm") {
      if (pp.HG == 0.0) throw Err() << "can't set Larmor position if "
                                    "field gradient is zero";
      cmd_sweep("B(cm)", args, &pp.H0, &pp.HT, -pp.HG); continue; }

    // Set/step/sweep RF field [G].
    if (cmd == "set_rf_field") {
      pp.HR0 = get_one_arg<double>(args); continue; }
    if (cmd == "step_rf_field") {
      pp.HR0 += get_one_arg<double>(args); continue; }
    if (cmd == "sweep_rf_field") {
      cmd_sweep("Brf(G)", args, &pp.HR0, &pp.HRT); continue; }

    // RF-field profile, gradient term [1/cm], quadratic term [1/cm^2].
    if (cmd == "set_rf_prof") {
      check_nargs(narg, 2);
      pp.HRGP = get_arg<double>(args[0]);
      pp.HRQP = get_arg<double>(args[1]);
      continue;
    }


    /*******************************************************/


    // Set/sweep relaxation time t_1 [s]
    if (cmd == "set_t1") {
      pp.T10 = get_one_arg<double>(args); continue; }
    if (cmd == "sweep_t1") {
      cmd_sweep("T1(s)", args, &pp.T10, &pp.T1T); continue; }

    // Set/sweep Leggett-Takagi relaxation time tau_f [s]
    if (cmd == "set_tf") {
      pp.TF0 = get_one_arg<double>(args); continue; }
    if (cmd == "sweep_tf") {
      cmd_sweep("TF(s)", args, &pp.TF0, &pp.TFT); continue; }

    // Set/sweep spin diffusion [cm^2/s]
    if (cmd == "set_diff") {
      pp.DF0 = get_one_arg<double>(args); continue; }
    if (cmd == "sweep_diff") {
      cmd_sweep("DIFF(cm2/s)", args, &pp.DF0, &pp.DFT); continue; }

    // Set and sweep spin-wave velocity c_parallel [cm/s]
    if (cmd == "set_cpar") {
      pp.CP0 = get_one_arg<double>(args); continue; }
    if (cmd == "sweep_cpar") {
      cmd_sweep("Cpar, cm/s", args, &pp.CP0, &pp.CPT); continue; }

    // Set and sweep Leggett frequency [Hz]
    if (cmd == "set_leggett_freq") {
      pp.LF0 = get_one_arg<double>(args); continue; }
    if (cmd == "sweep_leggett_freq") {
      cmd_sweep("fB(Hz)", args, &pp.LF0, &pp.LFT); continue; }

    /*******************************************************/


#ifdef HE3LIB
    if (cmd == "set_ttc_press") {
      check_nargs(narg, 2);
      double T = get_arg<double>(args[0]);
      double P = get_arg<double>(args[1]);
      set_he3tp(T, P);
      continue;
    }
#endif

    /*******************************************************/
    throw Err() << "skipping unknown command: " << cmd;
  }
  return 0;
}



/********************************************************************/
int
main(int argc, char *argv[]){
try{

  // program should be run with one argument - name of command file
  if (argc!=2){
    std::cerr << "Usage: " << argv[0] << " <command file>\n";
    return 1;
  }

  // open command file
  std::ifstream in_c(argv[1]);
  if (!in_c.good()){
    std::cerr << "Can't open command file: " << argv[1] << "\n";
    return 1;
  }

  // find prefix (command file name without extension)
  const char * pos1 = rindex(argv[1], '/');
  const char * pos2 = rindex(argv[1], '.');
  if (!pos1) pos1 = argv[1];
  pp.pref = pos2 && pos2>pos1+1 ? std::string(argv[1], pos2-argv[1]) : argv[1];

  std::ofstream out_l((pp.pref + ".run.log").c_str()); // log commands
  out_l << "# Commands and main parameters\n";

  write_pars(out_l);


  // main cycle
  int n=0;
  while (1) {

    // If we reach final time, read new cmd.
    // If it returns 1, finish the program
    if (pp.tcurr >= pp.tend && read_cmd(in_c, out_l)) break;
    // do the next step
    if (pp.solver) {
      pp.tcurr += pp.tstep;
//      pp.solver->step(pp.tcurr, (pp.tcurr>=pp.tend));
      std::cerr << pp.sweep_par_n << ":" << pp.sweep_par_o + pp.sweep_par_r*pp.tcurr << " ";
      pp.solver->step(pp.tcurr, false);

      // write magnetization (using mesh from the solver)
      if (pp.out_m) write_magn(*pp.out_m);

      // write pnm (using uniform mesh)
      std::vector<double> xsol = make_uniform_mesh(pp.npts, pp.cell_len);
      std::vector<double> usol = pp.solver->values(xsol, nder);
      pp.pnm_writers.write(xsol, usol);
       write_pars(out_l);
    }

    // flush files
    out_l.flush();
    if (pp.out_m) pp.out_m->flush();
  }

}
catch (const Err & e){
  std::cerr << "Error: " << e.str() << "\n";
}
catch (const pdecol_solver::Err & e){
  std::cerr << "Solver error: " << e.str() << "\n";
}
}
