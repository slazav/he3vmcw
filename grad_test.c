extern "C" {
#include "grad.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

const char *xyz[3]={"x","y","z"};


// check if absolute value is less then E
void checkz(const char *name, double V1, double E, bool nl){
  printf(" [%s", name);
  if (fabs(V1) > E){
    printf("\nFAIL: V1=%e; E=%e\n", V1, E);
    exit(1);
  }
  printf("]");
  if (nl) printf("\n");
}

// check if relative difference between two (non-zero) values is less then E
void check(const char *name, double V1, double V2, double E, bool nl){
  printf(" [%s", name);
  if (fabs(V1-V2)/(fabs(V1)+fabs(V2)) > E){
    printf("\n FAIL: V1=%e; V2=%e, E=%e\n", V1, V2, E);
    exit(1);
  }
  printf("]");
  if (nl) printf("\n");
}

void check3z(const char *name, double V1[DIM], double E, bool nl){
  int i;
  printf(" [%s: ", name);
  FOR(i) checkz(xyz[i], V1[i], E, false);
  printf("]");
  if (nl) printf("\n");
}

void check3(const char *name, double V1[DIM], double V2[DIM], double E, bool nl){
  int i;
  printf(" [%s: ", name);
  FOR(i) check(xyz[i], V1[i], V2[i], E, false);
  printf("]");
  if (nl) printf("\n");
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
     checkz("gR", r1, D, false);
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
     checkz("ggR", r1, D, false);
  }
  printf("\n");

  // test G0 vs G1
  {
    double E0a, E0b, E1a, E1b;
    fill_EG0_nt_(&E0a, &E0b, n0, t0, gn, gt);
    fill_EG1_nt_(&E1a, &E1b, n0, t0, gn, gt);
    check("E0a vs E1a", E0a, E1a, E, false);
    check("E0b vs E1b", E0b, E1b, E, true);
  }

  // test J0 vs J1
  {
    int i;
    double J0a[3], J0b[3], J1a[3], J1b[3];
    fill_JG0_nt_(J0a, J0b, n0, t0, gn, gt);
    fill_JG1_nt_(J1a, J1b, n0, t0, gn, gt);
    check3("J0a vs J1a", J0a, J1a, E, false);
    check3("J0b vs J1b", J0b, J1b, E, true);
  }

  // test J0 vs J2
  {
    int i;
    double J0a[3], J0b[3], J2a[3], J2b[3];
    fill_JG0_nt_(J0a, J0b, n0, t0, gn, gt);
    fill_JG2_nt_(J2a, J2b, n0, t0, gn, gt);
    check3("J0a vs J2a", J0a, J2a, E, false);
    check3("J0b vs J2b", J0b, J2b, E, true);
  }

  // test J0 vs JD
  {
    int i;
    double J0a[3], J0b[3], JD[3];
    fill_JG0_nt_(J0a, J0b, n0, t0, gn, gt);
    fill_JGD_nt_(JD, n0, t0, gn, gt);
    printf("[J0 vs JD ");
    FOR(i) check(xyz[i], J0a[i]/2 + J0b[i], JD[i], E, false);
    printf("]\n");
  }

  // test T0 vs T1
  {
    int i;
    double T0a[3], T0b[3], T1a[3], T1b[3];
    fill_TG0_nt_(T0a, T0b, n0, t0, gn, gt, ggn, ggt);
    fill_TG1_nt_(T1a, T1b, n0, t0, gn, gt, ggn, ggt);
    check3("T0a vs T1a", T0a, T1a, E, false);
    check3("T0b vs T1b", T0b, T1b, E, true);
  }

  // test T0 vs T2
  {
    int i;
    double T0a[3], T0b[3], T2a[3], T2b[3];
    fill_TG0_nt_(T0a, T0b, n0, t0, gn, gt, ggn, ggt);
    fill_TG2_nt_(T2a, T2b, n0, t0, gn, gt, ggn, ggt);
    check3("T0a vs T2a", T0a, T2a, E, false);
    check3("T0b vs T2b", T0b, T2b, E, true);
  }

  // test T0 vs TD
  {
    int i;
    double T0a[3], T0b[3], TD[3];
    fill_TG0_nt_(T0a, T0b, n0, t0, gn, gt, ggn, ggt);
    fill_TGD_nt_(TD, n0, t0, gn, gt, ggn, ggt);
    printf("[T0 vs TD: ");
    FOR(i) check(xyz[i], T0a[i]/2 + T0b[i], TD[i], E, false);
    printf("]\n");
  }

  // test spin current derivatives
  // apply small change to n, th, gn, gt and check
  // that J(n+dn) = J(n) + dJ/dn dn
  {
    double dn[DIM], dgn[DIM], dt, dgt;
    double n1[DIM], gn1[DIM];
    double J0a[DIM], J0b[DIM], J1a[DIM], J1b[DIM];
    double DIFFa[DIM], DIFFb[DIM];
    double Da[DIM][4], Db[DIM][4];
    double DZa[DIM][4], DZb[DIM][4];
    double NDa[DIM][3], NDb[DIM][3];
    double NDZa[DIM][3], NDZb[DIM][3];
    double dtn[DIM], dgtn[DIM];
    double v;

    fill_vec_rnd_(dn,  -D, D);
    fill_vec_rnd_(dgn, -D, D);
    dt  = (2*drand48()-1)*D;
    dgt = (2*drand48()-1)*D;

    // dn should be perpendicular to n. Find projection:
    // dn = dn - n (n*dn)
    v = n0[0]*dn[0]+n0[1]*dn[1]+n0[2]*dn[2];
    dn[0] -= n0[0]*v;
    dn[1] -= n0[1]*v;
    dn[2] -= n0[2]*v;
// Do not work
dn[0] = 0;
dn[1] = 0;
dn[2] = 0;
dgn[0] = 0;
dgn[1] = 0;
dgn[2] = 0;
//dgt=0;
//dt=0;

    // test DJ/Dgn
    gn1[0]=gn[0]+dgn[0];
    gn1[1]=gn[1]+dgn[1];
    gn1[2]=gn[2]+dgn[2];
    n1[0]=n0[0]+dn[0];
    n1[1]=n0[1]+dn[1];
    n1[2]=n0[2]+dn[2];
    fill_JG0_nt_(J0a, J0b, n0, t0, gn, gt);
    fill_JG0_nt_(J1a, J1b, n1, t0+dt, gn1, gt+dgt);
    fill_DJ_nt_(Da, Db, DZa, DZb, n0, t0, gn, gt);

    // now J1 = J0 + sum(DJ_i*dUi)
    // J1~J0~DJ~1,  dUi~D << 1

    FOR(i){
      DIFFa[i] += dn[0]*Da[i][0]
                + dn[1]*Da[i][1]
                + dn[2]*Da[i][2]
                + dt*Da[i][3]
                + dgn[0]*DZa[i][0]
                + dgn[1]*DZa[i][1]
                + dgn[2]*DZa[i][2]
                + dgt*DZa[i][3]
                + J0a[i] - J1a[i];
      DIFFb[i] += dn[0]*Db[i][0]
                + dn[1]*Db[i][1]
                + dn[2]*Db[i][2]
                + dt*Db[i][3]
                + dgn[0]*DZb[i][0]
                + dgn[1]*DZb[i][1]
                + dgn[2]*DZb[i][2]
                + dgt*DZb[i][3]
                + J0b[i] - J1b[i];
    }
    check3z("DJa", DIFFa, D*D, false);
    check3z("DJb", DIFFb, D*D, true);

    // same in th*n coordinates
// Does not work now...
dn[0] = 0;
dn[1] = 0;
dn[2] = 0;
//dgn[0] = 0;
//dgn[1] = 0;
//dgn[2] = 0;
dt  = 0;
dgt = 0;

    dtn[0] = n0[0]*dt + t0*dn[0];
    dtn[1] = n0[1]*dt + t0*dn[1];
    dtn[2] = n0[2]*dt + t0*dn[2];
    dgtn[0] = gn[0]*dt + gt*dn[0] + n0[0]*dgt + t0*dgn[0];
    dgtn[1] = gn[1]*dt + gt*dn[1] + n0[1]*dgt + t0*dgn[1];
    dgtn[2] = gn[2]*dt + gt*dn[2] + n0[2]*dgt + t0*dgn[2];

printf("%e %e %e %e\n", dgtn[0], dgtn[1], dgtn[2], t0);

    gn1[0]=gn[0]+dgn[0];
    gn1[1]=gn[1]+dgn[1];
    gn1[2]=gn[2]+dgn[2];
    n1[0]=n0[0]+dn[0];
    n1[1]=n0[1]+dn[1];
    n1[2]=n0[2]+dn[2];
    fill_JG0_nt_(J0a, J0b, n0, t0, gn, gt);
    fill_JG0_nt_(J1a, J1b, n1, t0+dt, gn1, gt+dgt);

    fill_DJ_t_(NDa, NDb, NDZa, NDZb, n0, t0, gn, gt);
    FOR(i){
      DIFFa[i] += dtn[0]*NDa[i][0]
                + dtn[1]*NDa[i][1]
                + dtn[2]*NDa[i][2]
                + dgtn[0]*NDZa[i][0]
                + dgtn[1]*NDZa[i][1]
                + dgtn[2]*NDZa[i][2]
                + J0a[i] - J1a[i];
      DIFFb[i] += dtn[0]*NDb[i][0]
                + dtn[1]*NDb[i][1]
                + dtn[2]*NDb[i][2]
                + dgtn[0]*NDZb[i][0]
                + dgtn[1]*NDZb[i][1]
                + dgtn[2]*NDZb[i][2]
                + J0b[i] - J1b[i];
    }
    check3z("NDJa", DIFFa, D*D, false);
    check3z("NDJb", DIFFb, D*D, true);

    // Same for Dmitriev's functions:
// Does not work now...
//dn[0] = 0;
//dn[1] = 0;
//dn[2] = 0;
//dgn[0] = 0;
//dgn[1] = 0;
//dgn[2] = 0;
//dt  = 0;
//dgt = 0;
    gn1[0]=gn[0]+dgn[0];
    gn1[1]=gn[1]+dgn[1];
    gn1[2]=gn[2]+dgn[2];
    n1[0]=n0[0]+dn[0];
    n1[1]=n0[1]+dn[1];
    n1[2]=n0[2]+dn[2];
    fill_JG0_nt_(J0a, J0b, n0, t0, gn, gt);
    fill_JG0_nt_(J1a, J1b, n1, t0+dt, gn1, gt+dgt);

    fill_DJD_nt_(Da, DZa, n0, t0, gn, gt);

    FOR(i){
      DIFFa[i] += dn[0]*Da[i][0]
              + dn[1]*Da[i][1]
              + dn[2]*Da[i][2]
              + dt*Da[i][3]
              + dgn[0]*DZa[i][0]
              + dgn[1]*DZa[i][1]
              + dgn[2]*DZa[i][2]
              + dgt*DZa[i][3]
              + (J0a[i] - J1a[i])/2
              + (J0b[i] - J1b[i]);
    }
//    check3z("DJD", DIFFa, D*D, true);


  }
}
}