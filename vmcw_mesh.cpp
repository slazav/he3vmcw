#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include "vmcw_pars.h"

// function aer_step_ is used from fortran function F.
// It use the global parameter structure
extern "C"{ extern struct pars_t pars_;}

// Aerogel profile.
// Returns 1 in the central
// part of the cell with fermi steps to 0 on edges.
//   x is coordinate
//   d is derivative (0 or 1)
// Use parameters:
//   CELL_LEN - cell length, cm
//   AER      -- if >0 then do step
//   AER_LEN  -- aerogel length / cell length
//   AER_CNT  -- center of aerogel area / cell length
//   AER_TRW  -- transition width / cell length
extern "C" {
double
aer_step_(double *x, int *d){

  if (!pars_.AER) return 0.0;
  double L = pars_.CELL_LEN;
  double l = pars_.AER_LEN;
  double c = pars_.AER_CNT;
  double w = pars_.AER_TRW;

  double xx = *x/L - c;
  double a = (fabs(xx) - l/2.0)/w;
  double dadx = (xx>0? 1:-1)/w/L;

  if(a < 82.0){
    if (*d==0) return 1.0/(1.0+exp(a));
    else       return dadx * exp(a)/(1.0+exp(a))/(1.0+exp(a));
  }
  return 0.0;
}
}

// Set the mesh
void
set_mesh(struct pars_t *p, std::vector<double> & x){
  if (x.size() < 2) return;

  // start with homogenious mesh with dx intervals
  double L = p->CELL_LEN;
  int N = x.size();
  double dx = L/(N-1);
  x[0] = -L/2.0;
  for (int k=0; k<100; k++){
    for (int i=0; i<N-1; i++){
      int der = 1;
      x[i+1] = x[i] + dx/(1.0+p->XMESH_K*fabs(aer_step_(&(x[i]),&der)));
    }
    // scale the whole mesh to fit CELL_LEN
    double d = L - (x[N-1]-x[0]);
    dx+=d/(N+1);
    if (fabs(d)<p->XMESH_ACC) break;
  }
}

// Save the mesh to file
void
save_mesh(struct pars_t *p, const std::vector<double> &x, const char *fname){
  std::ofstream ff(fname);
  for (int i=0; i<x.size(); i++){
    int d0=0, d1=1;
    double xx = x[i];
    double a0 = aer_step_(&xx,&d0);
    double a1 = aer_step_(&xx,&d1);
    ff << i << " " << xx << " " << a0 << " " << a1 << "\n";
  }
}