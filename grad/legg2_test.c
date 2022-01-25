extern "C" {
#include "grad.h"
#include "check.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

// test second and third Leggett equations:
// (\dot n,\dot\theta) = F(w,n,\theta)
// w = G(n,\theta, \dot n, \dot\theta)

int
main(){

  int i,j;
  double th, st, ct, C;
  double n[DIM], w[DIM], wm[DIM], wp[DIM], w1[DIM];
  double zz[DIM], nxz[DIM], nxdn[DIM];
  double v1[DIM], dn[DIM], dth;

  double E = 1e-6; // for comparison

  // initialize random number generator
  srand48(random());

  // random unit vectors n,w
  fill_vec_urnd_(n);
  fill_vec_urnd_(w);

  // angle theta
  th = M_PI*(2*drand48()-1);
  st=sin(th);
  ct=cos(th);
  C = st/(1.0-ct);

  // w-z, w+z
  FOR(i) { wm[i] = w[i]; wp[i]=w[i];}
  wm[2] -=1;
  wp[2] +=1;

  // Leggett equations:
  // \dot n  = -0.5* n\times[ w+z + C n\times(w-z)]
  // \dot\th = n(w-z)
  fill_vxv_(v1, n,wm);
  FOR(i) {v1[i] = wp[i] + C*v1[i];}
  fill_vxv_(dn, n, v1);
  FOR(i) {dn[i] = -0.5*dn[i];}

  dth=0;
  FOR(i) {dth += n[i]*wm[i];}

  printf("n,th:   %f %f %f  %f\n", n[0],n[1],n[2],th);
  printf("dn,dth: %f %f %f  %f\n", dn[0],dn[1],dn[2],dth);
  printf("w:      %f %f %f\n", w[0],w[1],w[2]);


  // inverse:
  // w = n*dth + (1-ct)*nxdn + st*dn
  //   + ct*z + (1-ct)*n*nz + st*nxz;
  zz[0]=zz[1]=0; zz[2]=1;
  fill_vxv_(nxz,  n,zz);
  fill_vxv_(nxdn, n,dn);

  FOR(i) w1[i] = n[i]*dth + (1.0-ct)*nxdn[i] + st*dn[i]
               + ct*zz[i] + (1.0-ct)*n[i]*n[2] + st*nxz[i];

  printf("w1:      %f %f %f\n", w1[0],w1[1],w1[2]);
  check3("w vs w1", w, w1, E,true);

  // For soliton mass calculation:
  // S = w - z + H
  // S = S0 + dS (equilibrium + time-dependent)
  // dS * (S0 - H)  = -2.0*(1.0-ct)*nxdn[z]
  {
    double S0[DIM], dS[DIM];
    double r = 0, r1 = 0;
    FOR(i) S0[i] = ct*zz[i] + (1.0-ct)*n[i]*n[2] + st*nxz[i] - zz[i];
    FOR(i) dS[i] = n[i]*dth + (1.0-ct)*nxdn[i] + st*dn[i];
    FOR(i) r += S0[i]*dS[i];
    r1 = -2.0*(1.0-ct)*nxdn[2];
    check("r", r, r1, E, true);
  }

}
}