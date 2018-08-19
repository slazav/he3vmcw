#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include "vmcw_pars.h"

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
double
aer_step(struct pars_t *p, double x, int d){

  if (!p->AER) return 0.0;
  double L = p->CELL_LEN;
  double l = p->AER_LEN;
  double c = p->AER_CNT;
  double w = p->AER_TRW;

  double xx = x/L - c;
  double a = (fabs(xx) - l/2.0)/w;
  double dadx = (xx>0? 1:-1)/w/L;

  if(a < 82.0){
    if (d==0) return 1.0/(1.0+exp(a));
    else      return dadx * exp(a)/(1.0+exp(a))/(1.0+exp(a));
  }
  return 0.0;
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
//std::cerr << ">>> " << x[0] << "\n";
  for (int k=0; k<100; k++){
    for (int i=0; i<N-1; i++){
      x[i+1] = x[i] + dx/(1.0+p->XMESH_K*fabs(aer_step(p,x[i],1)));
    }
    // scale the whole mesh to fit CELL_LEN
    double d = L - (x[N-1]-x[0]);
    dx+=d/(N+1);
    if (fabs(d)<p->XMESH_ACC) break;
  }
}

// Save the mesh to file
void
save_mesh(struct pars_t *p,
          const std::vector<double> &x,
          const std::string & fname){
  std::ofstream ff(fname);
  ff << "# N X AER AER'\n"
     << std::scientific << std::setprecision(6);
  for (int i=0; i<x.size(); i++){
    double xx = x[i];
    double a0 = aer_step(p,xx,0);
    double a1 = aer_step(p,xx,1);
    ff << i << " " << xx << " " << a0 << " " << a1 << "\n";
  }
}
