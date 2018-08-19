#ifndef VMCW_MESH_H
#define VMCW_MESH_H
#include "vmcw_pars.h"
#include <string>

double aer_step(pars_t *p, double x, int d);

// Set the mesh
void set_mesh(struct pars_t *p, std::vector<double> & x);

// Save the mesh to file
void save_mesh(struct pars_t *p, const std::vector<double> &x, const std::string &fname);

#endif
