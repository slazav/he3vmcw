#ifndef VMCW_MESH_H
#define VMCW_MESH_H
#include "vmcw_pars.h"

double aer_step(double x, int d);

// Set the mesh
void set_mesh(struct pars_t *p, std::vector<double> & x);

// Save the mesh to file
void save_mesh(struct pars_t *p, const std::vector<double> &x, const char *fname);

#endif
