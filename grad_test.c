extern "C" {
#include "grad.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

// check if absolute value is less then E
void checkz(double V1, double E, bool verb=false){
  assert(fabs(V1) < E);
}

// check if relative difference between two (non-zero) values is less then E
void check(double V1, double V2, double E, bool verb=false){
//  if (verb) printf("%e - %e\n", V1, V2);
  assert(fabs(V1-V2)/(fabs(V1)+fabs(V2)) < E);
}

void check3(double V1[DIM], double V2[DIM], double E, bool verb=false){
  int i;
  FOR(i){
//    if (verb) printf("[%d]: ", i);
    check(V1[i], V2[i], E, verb);
  }
}



main(){
  int i,j;
  double a0, b0, t0, n0[DIM], R0[DIM][DIM];
  double am, bm, tm, nm[DIM], Rm[DIM][DIM];
  double ap, bp, tp, np[DIM], Rp[DIM][DIM];

  double gn[DIM], gt, gR[DIM][DIM], ga, gb;
  double ggn[DIM], ggt, ggR[DIM][DIM];
  double D = 1e-6;  // grid size
  double E = 1e-6; // for comparison of exect formulas

  // initialize random number generator
  srand48(random()); 

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
     checkz(r1, D);
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
     checkz(r1, D);
     printf("test ggR: %e\n", r1);
  }

  // test G0 vs G1
  {
    double E0a, E0b, E1a, E1b;
    fill_EG0_nt_(&E0a, &E0b, n0, t0, gn, gt);
    fill_EG1_nt_(&E1a, &E1b, n0, t0, gn, gt);
    check(E0a, E1a, E);
    check(E0b, E1b, E);
    printf("test E0-E1: OK\n");
  }

  // test J0 vs J1
  {
    int i;
    double J0a[3], J0b[3], J1a[3], J1b[3];
    fill_JG0_nt_(J0a, J0b, n0, t0, gn, gt);
    fill_JG1_nt_(J1a, J1b, n0, t0, gn, gt);
    check3(J0a, J1a, E);
    check3(J0b, J1b, E);
    printf("test J0-J1: OK\n");
  }

  // test J0 vs J2
  {
    int i;
    double J0a[3], J0b[3], J2a[3], J2b[3];
    fill_JG0_nt_(J0a, J0b, n0, t0, gn, gt);
    fill_JG2_nt_(J2a, J2b, n0, t0, gn, gt);
    check3(J0a, J2a, E);
    check3(J0b, J2b, E);
    printf("test J0-J2: OK\n");
  }

  // test J0 vs JD
  {
    int i;
    double J0a[3], J0b[3], JD[3];
    fill_JG0_nt_(J0a, J0b, n0, t0, gn, gt);
    fill_JGD_nt_(JD, n0, t0, gn, gt);
    FOR(i) check(J0a[i]/2 + J0b[i], JD[i], E);
    printf("test J0-JD: OK\n");
  }

  // test T0 vs T1
  {
    int i;
    double T0a[3], T0b[3], T1a[3], T1b[3];
    fill_TG0_nt_(T0a, T0b, n0, t0, gn, gt, ggn, ggt);
    fill_TG1_nt_(T1a, T1b, n0, t0, gn, gt, ggn, ggt);
    check3(T0a, T1a, E);
    check3(T0b, T1b, E);
    printf("test T0-T1: OK\n");
  }

  // test T0 vs T2
  {
    int i;
    double T0a[3], T0b[3], T2a[3], T2b[3];
    fill_TG0_nt_(T0a, T0b, n0, t0, gn, gt, ggn, ggt);
    fill_TG2_nt_(T2a, T2b, n0, t0, gn, gt, ggn, ggt);
    check3(T0a, T2a, E);
    check3(T0b, T2b, E);
    printf("test T0-T2: OK\n");
  }

  // test T0 vs TD
  {
    int i;
    double T0a[3], T0b[3], TD[3];
    fill_TG0_nt_(T0a, T0b, n0, t0, gn, gt, ggn, ggt);
    fill_TGD_nt_(TD, n0, t0, gn, gt, ggn, ggt);
    FOR(i) check(T0a[i]/2 + T0b[i], TD[i], E, true);
    printf("test T0-TD: OK\n");
  }

  // test spin current derivatives
  // apply small change to n, th, gn, gt and check
  // that J(n+dn) = J(n) + dJ/dn dn
  {
    double dn[DIM], dgn[DIM], dt, dgt;
    double n1[DIM], gn1[DIM];
    double J0a[3], J0b[3], J1a[3], J1b[3];
    double DJa[3][3], DJb[3][3];
    double v;
    fill_vec_rnd_(dn, -D, D);
    fill_vec_rnd_(gn, -D, D);
    dt  = (2*drand48()-1)*D;
    dgt = (2*drand48()-1)*D;
    // dn should be perpendicular to n. Find projection:
    // dn = dn - n (n*dn)
    v = n0[0]*dn[0]+n0[1]*dn[1]+n0[2]*dn[2];
    dn[0] -= n0[0]*v;
    dn[1] -= n0[1]*v;
    dn[2] -= n0[2]*v;

    // test DJ/Dgn
    gn1[0]=gn[0]+dgn[0];
    gn1[1]=gn[1]+dgn[1];
    gn1[2]=gn[2]+dgn[2];
    fill_JG0_nt_(J0a, J0b, n0, t0, gn, gt);
    fill_JG0_nt_(J1a, J1b, n0, t0, gn1, gt);

    fill_DJGgn_nt_(DJa, DJb, n0, t0, gn, gt);

    FOR(i){
      J0a[i] += dgn[0]*DJa[i][0] + dgn[1]*DJa[i][1] + dgn[2]*DJa[i][2];
      J0b[i] += dgn[0]*DJb[i][0] + dgn[1]*DJb[i][1] + dgn[2]*DJb[i][2];
    }
    check3(J0a, J1a, D);
    check3(J0b, J1b, D);
    printf("test DJgn: OK\n");
  }
}
}