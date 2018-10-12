#ifndef GRAD_H
#define GRAD_H
/************************************************************
1D Gradient energy terms for Leggett equations.

- Gradiend energy:
  E = K1 \nabla_j Rak \nabla_j Rak
    + K2 \nabla_j Raj \nabla_k Rak
    + K3 \nabla_j Rak \nabla_k Raj

- Spin current
  J_ak = e_{abc} Rcj [K1 \nabla_k Rbj + (K2+K3) \nabla_j Rbk]

- Gradient torque:
  T = -\nabla_m J_am =
    = -e_{abc} Rcj [K1 \nabla_k\nabla_k Rbj + (K2+K3) \nabla_k\nabla_j Rbk]

See also http://github.com/slazav/he3en for 3D formulas.
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

/// Fill a zero vector.
inline void fill_vec_zero_(double v[DIM]){
  memset(v, 0, DIM*sizeof(double)); }

/// Fill a zero matrix.
inline void fill_matr_zero_(double v[DIM][DIM]){
  memset(v, 0, DIM*DIM*sizeof(double)); }

/// Fill a zero matrix DIMxDIMxDIM.
inline void fill_matr3_zero_(double v[DIM][DIM][DIM]){
  memset(v, 0, DIM*DIM*DIM*sizeof(double)); }

/// Print 3x3 matrix
inline void print_matr_(double v[DIM][DIM]){
  int i,j;
  FOR(i){ FOR(j){ printf(" %e", v[i][j]); } printf("\n"); }
}


/// Fill e_{ijk} matrix.
inline void fill_ee_(double ee[DIM][DIM][DIM]){
  fill_matr3_zero_(ee);
  ee[0][1][2] = ee[1][2][0] = ee[2][0][1] = +1.0;
  ee[2][1][0] = ee[1][0][2] = ee[0][2][1] = -1.0;
}

/// Fill delta_{ij} matrix.
inline void fill_dd_(double dd[DIM][DIM]){
  int i;
  fill_matr_zero_(dd);
  FOR(i) dd[i][i] = 1.0;
}

/// Fill e_{ijk}*n_k matrix.
inline void fill_en_(double en[DIM][DIM], const double n[DIM]){
  en[0][1] = +n[2]; en[1][0] = -n[2];
  en[1][2] = +n[0]; en[2][1] = -n[0];
  en[2][0] = +n[1]; en[0][2] = -n[1];
  en[0][0] = en[1][1] = en[2][2] = 0.0;
}

/// Vector product
inline void fill_vxv_(double r[DIM], const double v1[DIM], const double v2[DIM]){
  r[0] = v1[1]*v2[2] - v1[2]*v2[1];
  r[1] = v1[2]*v2[0] - v1[0]*v2[2];
  r[2] = v1[0]*v2[1] - v1[1]*v2[0];
}


/// Fill a unit vector using polar and azimuthal abgles a,b
inline void fill_vec_ab_(double v[DIM], const double a, const double b){
  v[0]=sin(b)*cos(a);
  v[1]=sin(b)*sin(a);
  v[2]=cos(b);
}

/// Fill a random vector.
inline void fill_vec_rnd_(double v[DIM], const double min, const double max){
  int i;
  FOR(i) v[i] = min+drand48()*(max-min);
}

/// Fill a random unit vector
inline void fill_vec_urnd_(double v[DIM]){
  double a = 2.0*M_PI*drand48(); // azimuthal angle
  double b = M_PI*drand48(); // polar angle
  fill_vec_ab_(v,a,b);
}

/***********************************************************/
/// Fill a rotation matrix using rotation axis (n) and rotation angle (theta)
void fill_R_nt_(double R[DIM][DIM], const double n[DIM], const double th);

/// calculate rotation matrix gradient in n and theta coordinates
void fill_gR_nt_(double gR[DIM][DIM],
                     const double n[DIM], const double t,
                     const double gn[DIM], const double gt);

// calculate rotation matrix double gradient in n and theta coordinates
void fill_ggR_nt_(double ggR[DIM][DIM],
                     const double n[DIM], const double t,
                     const double gn[DIM], const double gt,
                     const double ggn[DIM], const double ggt);

/***********************************************************/
/// Calculate gradient energies Ea, Eb (just by definition)
void fill_EG0_nt_(double *Ea, double *Eb,
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt);

/// Calculate gradient energies Ea, Eb (v1, for n,th)
void fill_EG1_nt_(double *Ea, double *Eb,
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt);

/***********************************************************/
/// Calculate spin currents Ja, Jb (just by definition)
void fill_JG0_nt_(double Ja[DIM], double Jb[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt);

/// Calculate spin currents Ja, Jb (v1 - for n,th)
void fill_JG1_nt_(double Ja[DIM], double Jb[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt);

/// Calculate spin currents Ja, Jb (v2 - in coordinates)
void fill_JG2_nt_(double Ja[DIM], double Jb[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt);

/// Calculate spin current J = Ja/2 + Jb (as in Dmitriev's program)
void fill_JGD_nt_(double J[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt);

/***********************************************************/
/// Calculate gradient torques Ta, Tb (just by definition)
/// T = -\nabla J
void fill_TG0_nt_(double Ta[DIM], double Tb[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt,
                 const double ggn[DIM], const double ggt);

// Calculate gradient torques Ta, Tb (v1, in n and th)
void fill_TG1_nt_(double Ta[DIM], double Tb[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt,
                 const double ggn[DIM], const double ggt);

// Calculate gradient torques Ta, Tb (v2, in conponents)
void fill_TG2_nt_(double Ta[DIM], double Tb[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt,
                 const double ggn[DIM], const double ggt);

// Calculate gradient torque TD = Ta/2+Tb (as in Dmitriev's program)
void fill_TGD_nt_(double T[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt,
                 const double ggn[DIM], const double ggt);


}// extern
#endif
