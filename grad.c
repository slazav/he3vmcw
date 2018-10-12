/************************************************************
1D Gradient energy terms for Leggett equations.

See also http://github.com/slazav/he3en for 3D formulas etc.
All functions are inside <extern "C"> and have names with "_"
suffix for using in Fortran code.

*************************************************************/

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define DIM 3
#define FOR(i) for(i=0; i<DIM; i++)

// Fill a zero vector.
inline void fill_vec_zero_(double v[DIM]){
  memset(v, 0, DIM*sizeof(double));
}

// Fill a zero matrix.
inline void fill_matr_zero_(double v[DIM][DIM]){
  memset(v, 0, DIM*DIM*sizeof(double));
}

// Fill a zero matrix DIMxDIMxDIM.
inline void fill_matr3_zero_(double v[DIM][DIM][DIM]){
  memset(v, 0, DIM*DIM*DIM*sizeof(double));
}

// Print 3x3 matrix
inline void print_matr_(double v[DIM][DIM]){
  int i,j;
  FOR(i){ FOR(j){ printf(" %e", v[i][j]); } printf("\n"); }
}


// Fill e_{ijk} matrix.
inline void fill_ee_(double ee[DIM][DIM][DIM]){
  fill_matr3_zero_(ee);
  ee[0][1][2] = ee[1][2][0] = ee[2][0][1] = +1.0;
  ee[2][1][0] = ee[1][0][2] = ee[0][2][1] = -1.0;
}

// Fill delta_{ij} matrix.
inline void fill_dd_(double dd[DIM][DIM]){
  int i;
  fill_matr_zero_(dd);
  FOR(i) dd[i][i] = 1.0;
}

// Fill e_{ijk}*n_k matrix.
inline void fill_en_(double en[DIM][DIM], const double n[DIM]){
  en[0][1] = +n[2]; en[1][0] = -n[2];
  en[1][2] = +n[0]; en[2][1] = -n[0];
  en[2][0] = +n[1]; en[0][2] = -n[1];
  en[0][0] = en[1][1] = en[2][2] = 0.0;
}

// Fill a unit vector using polar and azimuthal abgles a,b
inline void fill_vec_ab_(double v[DIM], const double a, const double b){
  v[0]=sin(b)*cos(a);
  v[1]=sin(b)*sin(a);
  v[2]=cos(b);
}

// Fill a random vector.
void fill_vec_rnd_(double v[DIM], const double min, const double max){
  int i;
  FOR(i) v[i] = min+drand48()*(max-min);
}

// Fill a random unit vector
void fill_vec_urnd_(double v[DIM]){
  double a = 2.0*M_PI*drand48(); // azimuthal angle
  double b = M_PI*drand48(); // polar angle
  fill_vec_ab_(v,a,b);
}

/***********************************************************/
// Fill a rotation matrix using rotation axis (n) and
// rotation angle (theta)
void fill_R_nt_(double R[DIM][DIM], const double n[DIM], const double th){
  double ct=cos(th), st=sin(th);
  double dd[DIM][DIM];
  double en[DIM][DIM];
  int i,j;
  fill_dd_(dd);
  fill_en_(en, n);
  FOR(i) FOR(j) R[i][j] = ct*dd[i][j] + (1.-ct)*n[i]*n[j] - st*en[i][j];
}

// calculate rotation matrix gradient in n and theta coordinates
void fill_gR_nt_(double gR[DIM][DIM],
                     const double n[DIM], const double t,
                     const double gn[DIM], const double gt) {
   double dd[DIM][DIM];
   double en[DIM][DIM];
   double eg[DIM][DIM];
   int a,j;
   fill_dd_(dd);
   fill_en_(en, n);
   fill_en_(eg, gn);
   double ct=cos(t), ctm=(1.0-ct), st=sin(t);

   FOR(a) FOR(j) gR[a][j] =
       st*(n[a]*n[j] - dd[a][j])*gt
    + ctm*(n[j]*gn[a] + n[a]*gn[j])
    - ct*en[a][j]*gt - st*eg[a][j];
}

// calculate rotation matrix double gradient in n and theta coordinates
void fill_ggR_nt_(double ggR[DIM][DIM],
                     const double n[DIM], const double t,
                     const double gn[DIM], const double gt,
                     const double ggn[DIM], const double ggt) {
   double dd[DIM][DIM];
   double en[DIM][DIM];
   double eg[DIM][DIM];
   double egg[DIM][DIM];
   int a,j;
   fill_dd_(dd);
   fill_en_(en, n);
   fill_en_(eg, gn);
   fill_en_(egg, ggn);
   double ct=cos(t), ctm=(1.0-ct), st=sin(t);

   FOR(a) FOR(j) ggR[a][j] =
       + (n[a]*n[j]-dd[a][j])*(ct*gt*gt + st*ggt)
       + (n[a]*gn[j]+n[j]*gn[a])*2.0*st*gt
       + ctm*(gn[j]*gn[a] + gn[a]*gn[j])
       + ctm*(n[j]*ggn[a] + n[a]*ggn[j])
       + en[a][j]*(st*gt*gt - ct*ggt) - 2.0*ct*eg[a][j]*gt - st*egg[a][j];
}

/***********************************************************/
// Calculate gradient energies Ea, Eb (just by definition)
void fill_EG0_nt_(double Ea, double Eb, double Ec,
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt) {
  int a,j;
  double gR[DIM][DIM];
  fill_gR_nt_(gR,n,t,gn,gt);
  Ea=Eb=0.0;

  FOR(a) FOR(j){
    Ea += gR[a][j]*gR[a][j]; // K1
    Eb += gR[a][3]*gR[a][3]; // (K2+K3)
  }
}

/***********************************************************/
// Calculate spin currents Ja, Jb (just by definition)
void fill_JG0_nt_(double Ja[DIM], double Jb[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt) {
  int a,b,c,j;
  double R[DIM][DIM], gR[DIM][DIM];
  double ee[DIM][DIM][DIM];
  fill_ee_(ee);
  fill_vec_zero_(Ja);
  fill_vec_zero_(Jb);
  fill_R_nt_(R,n,t);
  fill_gR_nt_(gR,n,t,gn, gt);
  FOR(a) FOR(b) FOR(c) FOR(j) Ja[a] += ee[a][b][c]*R[c][j]*gR[b][j];
  FOR(a) FOR(b) FOR(c)        Jb[a] += ee[a][b][c]*R[c][3]*gR[b][3];
}

/***********************************************************/
// Calculate gradient torques Ta, Tb (just by definition)
void fill_TG0_nt_(double Ta[DIM], double Tb[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt,
                 const double ggn[DIM], const double ggt) {
  int a,b,c,j;
  double R[DIM][DIM], ggR[DIM][DIM];
  double ee[DIM][DIM][DIM];
  fill_ee_(ee);
  fill_vec_zero_(Ta);
  fill_vec_zero_(Tb);
  fill_R_nt_(R,n,t);
  fill_ggR_nt_(ggR, n,t, gn, gt, ggn, ggt);
  FOR(a) FOR(b) FOR(c) FOR(j) Ta[a] += ee[a][b][c]*R[c][j]*ggR[b][j];
  FOR(a) FOR(b) FOR(c)        Tb[a] += ee[a][b][c]*R[c][3]*ggR[b][3];
}

/***********************************************************/
// Calculate gradient energies Ea, Eb (just by definition)
void fill_EG1_nt_(double Ea, double Eb, double Ec,
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt) {
  int a,j;
  double gR[DIM][DIM];
  fill_gR_nt_(gR,n,t,gn,gt);
  Ea=Eb=0.0;

%  for a=1:3; for j=1:3; for k=1:3; for l=1:3; for m=1:3;
%    e1 = e1 +...
%      (((1-ct)*(dd(a,l)*n(j) + dd(j,l)*n(a)) - st*ee(a,j,l)) * gn(l,k) + ...
%      (st*(n(a)*n(j) - dd(a,j))*dd(j,l) - ct*ee(a,j,l)*n(l)) * gt(k)) * ...
%      (((1-ct)*(dd(a,m)*n(j) + dd(j,m)*n(a)) - st*ee(a,j,m)) * gn(m,k) + ...
%      (st*(n(a)*n(j) - dd(a,j))*dd(j,m) - ct*ee(a,j,m)*n(m)) * gt(k));
%    e2 = e2 +...
%      (((1-ct)*(dd(a,l)*n(j) + dd(j,l)*n(a)) - st*ee(a,j,l)) * gn(l,k) + ...
%      (st*(n(a)*n(j) - dd(a,j))*dd(j,l) - ct*ee(a,j,l)*n(l)) * gt(k)) * ...
%      (((1-ct)*(dd(a,m)*n(k) + dd(k,m)*n(a)) - st*ee(a,k,m)) * gn(m,j) + ...
%      (st*(n(a)*n(k) - dd(a,k))*dd(k,m) - ct*ee(a,k,m)*n(m)) *gt(j));
%    e3 = e3 +...
%      (((1-ct)*(dd(a,l)*n(j) + dd(j,l)*n(a)) - st*ee(a,j,l)) * gn(l,j) + ...
%      (st*(n(a)*n(j) - dd(a,j))*dd(j,l) - ct*ee(a,j,l)*n(l)) * gt(j)) * ...
%      (((1-ct)*(dd(a,m)*n(k) + dd(k,m)*n(a)) - st*ee(a,k,m)) * gn(m,k) + ...
%      (st*(n(a)*n(k) - dd(a,k))*dd(k,m) - ct*ee(a,k,m)*n(m)) *gt(k));
%  end; end; end; end; end;

}

/***********************************************************/
main(){
  int i,j;
  double a0, b0, t0, n0[DIM], R0[DIM][DIM];
  double am, bm, tm, nm[DIM], Rm[DIM][DIM];
  double ap, bp, tp, np[DIM], Rp[DIM][DIM];

  double gn[DIM], gt, gR[DIM][DIM], ga, gb;
  double ggn[DIM], ggt, ggR[DIM][DIM];
  double D = 1e-6; // grid size

  // make a random rotation matrix: R0
  // make left and right matricis Rm, Rp,
  // assuming gradients ~1 and grid size D
  t0 = M_PI*(2*drand48()-1);
  a0 = 2*M_PI*drand48();
  b0 = M_PI*drand48();

  am = a0 + drand48()*D;
  bm = b0 + drand48()*D;
  tm = t0 + drand48()*D;

  ap = a0 + drand48()*D;
  bp = b0 + drand48()*D;
  tp = t0 + drand48()*D;

  fill_vec_ab_(n0, a0,b0);
  fill_vec_ab_(nm, am,bm);
  fill_vec_ab_(np, ap,bp);

  fill_R_nt_(R0, n0,t0);
  fill_R_nt_(Rm, nm,tm);
  fill_R_nt_(Rp, np,tp);

  // Make first and second-order derivatives of n, th, R
  // For small D this should be accurate:
  gt  = (tp-tm)/D;
  ggt = (tp+tm-2.*t0)/D/D;
  FOR(i) {
    gn[i]  = (np[i]-nm[i])/D;
    ggn[i] = (np[i]+nm[i]-2.*n0[i])/D/D;
  }
  FOR(i) FOR(j) {
    gR[i][j]  = (Rp[i][j]-Rm[i][j])/D;
    ggR[i][j] = (Rp[i][j]+Rm[i][j]-2.*R0[i][j])/D/D;
  }

  // test 1: check formula for gR
  {
     int a,j;
     double r0,r1;
     double gR1[DIM][DIM];
     fill_gR_nt_(gR1, n0,t0, gn, gt);

     // calculate difference with gR:
     FOR(a) FOR(j){
       r1 += pow(gR1[a][j]-gR[a][j], 2);
       r0 += pow(gR[a][j], 2);
     }
     r1 = sqrt(r1)/sqrt(r0);
     printf("test gR:  %e\n", r1);
  }

  // test 2: check formula for ggR
  {
     int a,j;
     double r0,r1;
     double ggR1[DIM][DIM];
     fill_ggR_nt_(ggR1, n0,t0, gn, gt, ggn, ggt);

     // calculate difference with gR:
     FOR(a) FOR(j){
       r1 += pow(ggR1[a][j]-ggR[a][j], 2);
       r0 += pow(ggR[a][j], 2);
     }
     r1 = sqrt(r1)/sqrt(r0);
     printf("test ggR: %e\n", r1);
  }

}

}