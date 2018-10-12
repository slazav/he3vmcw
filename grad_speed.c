extern "C" {
#include "grad.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

long get_time(){
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return ((long)tv.tv_sec)*1000000 + (long)tv.tv_usec;
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


  // test G
  {
    int i;
    double Ea, Eb;
    int N = 100000;
    long t1;

    t1 = get_time();
    for (i=0; i<N; i++) fill_EG0_nt_(&Ea, &Eb, n0, t0, gn, gt);
    t1 = get_time() - t1;
    printf("E0: %e us\n", (double)t1/(double)N);

    t1 = get_time();
    for (i=0; i<N; i++) fill_EG1_nt_(&Ea, &Eb, n0, t0, gn, gt);
    t1 = get_time() - t1;
    printf("E1: %e us\n", (double)t1/(double)N);
  }


  // test J
  {
    int i;
    int N = 100000;
    long t1;
    double Ja[3], Jb[3];

    t1 = get_time();
    for (i=0; i<N; i++) fill_JG0_nt_(Ja, Jb, n0, t0, gn, gt);
    t1 = get_time() - t1;
    printf("J0: %e us\n", (double)t1/(double)N);

    t1 = get_time();
    for (i=0; i<N; i++) fill_JG1_nt_(Ja, Jb, n0, t0, gn, gt);
    t1 = get_time() - t1;
    printf("J1: %e us\n", (double)t1/(double)N);

    t1 = get_time();
    for (i=0; i<N; i++) fill_JG2_nt_(Ja, Jb, n0, t0, gn, gt);
    t1 = get_time() - t1;
    printf("J2: %e us\n", (double)t1/(double)N);

    t1 = get_time();
    for (i=0; i<N; i++) fill_JGD_nt_(Ja, n0, t0, gn, gt);
    t1 = get_time() - t1;
    printf("JD: %e us\n", (double)t1/(double)N);

  }


  // test T
  {
    int i;
    int N = 100000;
    long t1;
    double Ta[3], Tb[3];

    t1 = get_time();
    for (i=0; i<N; i++) fill_TG0_nt_(Ta, Tb, n0, t0, gn, gt, ggn, ggt);
    t1 = get_time() - t1;
    printf("T0: %e us\n", (double)t1/(double)N);

    t1 = get_time();
    for (i=0; i<N; i++) fill_TG1_nt_(Ta, Tb, n0, t0, gn, gt, ggn, ggt);
    t1 = get_time() - t1;
    printf("T1: %e us\n", (double)t1/(double)N);

    t1 = get_time();
    for (i=0; i<N; i++) fill_TGD_nt_(Ta, n0, t0, gn, gt, ggn, ggt);
    t1 = get_time() - t1;
    printf("TD: %e us\n", (double)t1/(double)N);
  }


/*
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

  // test T0 vs TD
  {
    int i;
    double T0a[3], T0b[3], TD[3];
    fill_TG0_nt_(T0a, T0b, n0, t0, gn, gt, ggn, ggt);
    fill_TGD_nt_(TD, n0, t0, gn, gt, ggn, ggt);
    FOR(i) check(T0a[i]/2 + T0b[i], TD[i], E, true);
    printf("test T0-TD: OK\n");
  }
*/
}

}