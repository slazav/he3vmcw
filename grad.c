extern "C" {
#include "grad.h"

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
void fill_EG0_nt_(double *Ea, double *Eb,
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt) {
  int a,j;
  double gR[DIM][DIM];
  fill_gR_nt_(gR,n,t,gn,gt);
  *Ea=*Eb=0.0;

  FOR(a) FOR(j) *Ea += pow(gR[a][j],2); // K1
  FOR(a)        *Eb += pow(gR[a][2],2); // (K2+K3)
}

// v1, in n-th coordinates
void fill_EG1_nt_(double *Ea, double *Eb,
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt) {
  int a,j;
  double ct=cos(t), ctm=(1.0-ct), st=sin(t);
  double en[DIM][DIM];
  double eg[DIM][DIM];
  double dd[DIM][DIM];
  fill_en_(en, n);
  fill_en_(eg, gn);
  fill_dd_(dd);
  *Ea=*Eb=0.0;

  FOR(a) FOR(j) *Ea +=
      pow((ctm*(n[j]*gn[a] + n[a]*gn[j]) - st*eg[a][j])  +
      (st*(n[a]*n[j] - dd[a][j]) - ct*en[a][j]) * gt, 2);

  FOR(a) *Eb +=
      pow((ctm*(n[2]*gn[a] + n[a]*gn[2]) - st*eg[a][2]) +
      (st*(n[a]*n[2] - dd[a][2]) - ct*en[a][2]) * gt, 2);

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
  FOR(a) FOR(b) FOR(c)        Jb[a] += ee[a][b][c]*R[c][2]*gR[b][2];
}

// Calculate spin currents Ja, Jb (v1)
void fill_JG1_nt_(double Ja[DIM], double Jb[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt) {
  int a,b,c,j;
  double ct=cos(t), ctm=(1.0-ct), st=sin(t);
  double ee[DIM][DIM][DIM];
  double nxng[DIM];
  double en[DIM][DIM];
  double eg[DIM][DIM];
  fill_ee_(ee);
  fill_en_(en,n);
  fill_en_(eg,gn);
  fill_vxv_(nxng,n,gn);

  fill_vec_zero_(Ja);
  fill_vec_zero_(Jb);

  FOR(a) Ja[a] += -2*(n[a]*gt + st*gn[a] + ctm*nxng[a]);

  FOR(a) Jb[a] +=
     - n[a]*gt  *(1-ctm*n[2]*n[2])
     + st*ctm   *2*n[a]*n[2]*gn[2]
     + st       *en[2][a]*n[2]*gt
     + ct*ctm   *(en[2][a]*gn[2] + eg[2][a]*n[2])
     - ctm*ctm  *nxng[a]*n[2]*n[2]
     - st*gn[a] *(ct + ctm*n[2]*n[2]);
    ;
  Jb[2] +=  ct*n[2]*gt - st*(1-2*ct)*gn[2] - st*st*nxng[2];
}

/// Calculate spin currents Ja, Jb (v2 - in coordinates)
void fill_JG2_nt_(double Ja[DIM], double Jb[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt) {
  double ct=cos(t), ctm=(1.0-ct), st=sin(t);
  double nz2 = n[2]*n[2];

  Ja[0] = -2*(n[0]*gt + st*gn[0] + ctm *(n[1]*gn[2] - n[2]*gn[1]));
  Ja[1] = -2*(n[1]*gt + st*gn[1] + ctm *(n[2]*gn[0] - n[0]*gn[2]));
  Ja[2] = -2*(n[2]*gt + st*gn[2] + ctm *(n[0]*gn[1] - n[1]*gn[0]));

  Jb[0] =
   - (1-ctm*nz2)*n[0]*gt
   + st*n[1]*n[2]*gt
   -  st*(ct + ctm*nz2) *gn[0]
   + ctm*(ct + ctm*nz2) *n[2]*gn[1]
   + ctm*(ct - ctm*nz2) *n[1]*gn[2]
   + 2*st*ctm*n[0]*n[2]*gn[2];

  Jb[1] =
   - (1-ctm*nz2)*n[1]*gt
   - st*n[0]*n[2]*gt
   - ctm*(ct+ctm*nz2) *n[2]*gn[0]
   -  st*(ct+ctm*nz2) *gn[1]
   - ctm*(ct-ctm*nz2) *n[0]*gn[2]
   + 2*st*ctm*n[1]*n[2]*gn[2];

  Jb[2] =
   - ctm*(1 - nz2) * n[2]*gt
   + (st*st + ctm*ctm *nz2) *n[1]*gn[0]
   - (st*st + ctm*ctm *nz2) *n[0]*gn[1]
   - st*ctm*(1-nz2)*gn[2];
}

/// Calculate spin current J = Ja/2 + Jb (as in Dmitriev's program)
void fill_JGD_nt_(double J[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt) {
  double ct=cos(t), ctm=(1.0-ct), st=sin(t);
  double FTN=ctm*(n[0]*gn[1]-n[1]*gn[0]) - st*gn[2] - gt*n[2];
  J[0] = -2.0*(gt*n[0]+st*gn[0]+ctm*(n[1]*gn[2]-gn[1]*n[2]))
       - (ctm*n[0]*n[2]+n[1]*st)*FTN;
  J[1] = -2.0*(gt*n[1]+st*gn[1]-ctm*(n[0]*gn[2]-gn[0]*n[2]))
       - (ctm*n[1]*n[2]-n[0]*st)*FTN;
  J[2] = -2.0*(gt*n[2]+st*gn[2]+ctm*(n[0]*gn[1]-gn[0]*n[1]))
       - (ctm*n[2]*n[2]+ct)*FTN;
}

/// Spin current derivative, D(J)/D(U), D(J)/D(U') where U = nx,ny,nz,th
void fill_DJ_nt_(double DJDUa[DIM][4],  double DJDUb[DIM][4],
                 double DJDUXa[DIM][4], double DJDUXb[DIM][4],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt){

  double ct=cos(t), ctm=(1.0-ct), st=sin(t);
  double nz2 = n[2]*n[2];

  DJDUa[0][0] = -2*gt;
  DJDUa[0][1] = -2*ctm*gn[2];
  DJDUa[0][2] = +2*ctm*gn[1];

  DJDUa[1][0] = +2*ctm*gn[2];
  DJDUa[1][1] = -2*gt;
  DJDUa[1][2] = -2*ctm*gn[0];

  DJDUa[2][0] = -2*ctm*gn[1];
  DJDUa[2][1] = +2*ctm*gn[0];
  DJDUa[2][2] = -2*gt;

  // d(J)/d(n)
  DJDUa[0][3] = -2*(n[0]*gt + st*gn[0] + ctm *(n[1]*gn[2] - n[2]*gn[1]));
  DJDUa[1][3] = -2*(n[1]*gt + st*gn[1] + ctm *(n[2]*gn[0] - n[0]*gn[2]));
  DJDUa[2][3] = -2*(n[2]*gt + st*gn[2] + ctm *(n[0]*gn[1] - n[1]*gn[0]));

  DJDUb[0][0] = - (1-ctm*nz2)*gt + 2*st*ctm*n[2]*gn[2];
  DJDUb[0][1] = + st*n[2]*gt + ctm*(ct-ctm*nz2)*gn[2];
  DJDUb[0][2] =
   + 2.0*ctm*n[2]*n[0]*gt + st*n[1]*gt
   - 2.0*st*ctm*n[2]*gn[0]
   + ctm*(ct + 3.0*ctm*n[2]*n[2])*gn[1]
   - 2.0*ctm*ctm*n[2]*n[1]*gn[2]
   + 2.0*st*ctm*n[0]*gn[2];

  DJDUb[1][0] = - st*n[2]*gt - ctm*(ct-ctm*nz2)*gn[2];
  DJDUb[1][1] = - (1-ctm*nz2)*gt + 2*st*ctm*n[2]*gn[2];
  DJDUb[1][2] =
   + 2.0*ctm*n[2]*n[1]*gt - st*n[0]*gt
   - ctm*(ct + 3.0*ctm*n[2]*n[2])*gn[0]
   - 2.0*st*ctm*n[2] *gn[1]
   + 2.0*ctm*ctm*n[2]*n[0]*gn[2]
   + 2.0*st*ctm*n[1]*gn[2];

  DJDUb[2][0] = - (st*st + ctm*ctm *nz2)*gn[1];
  DJDUb[2][1] = + (st*st + ctm*ctm *nz2)*gn[0];
  DJDUb[2][2] =
   - ctm*(1 - 3*n[2]*n[2])*gt
   + 2.0*ctm*ctm*n[2]*n[1]*gn[0]
   - 2.0*ctm*ctm*n[2]*n[0]*gn[1]
   + 2.0*st*ctm*n[2]*gn[2];

  // d(J)/d(t) -- todo
  DJDUb[0][3] =
     st*nz2*n[0]*gt
   + ct*n[1]*n[2]*gt
   - ct*(ct + ctm*nz2) *gn[0]
   + st*(ct + ctm*nz2) *n[2]*gn[1]
   + st*(ct - ctm*nz2) *n[1]*gn[2]
   + 2*ct*ctm*n[0]*n[2]*gn[2]
   +  st*st*(1-nz2) *gn[0]
   - ctm*st*(1-nz2) *n[2]*gn[1]
   - ctm*st*(1+nz2) *n[1]*gn[2]
   + 2*st*st*n[0]*n[2]*gn[2];
  DJDUb[1][3] =
     st*nz2*n[1]*gt
   - ct*n[0]*n[2]*gt
   - st*(ct+ctm*nz2) *n[2]*gn[0]
   - ct*(ct+ctm*nz2) *gn[1]
   - st*(ct-ctm*nz2) *n[0]*gn[2]
   + 2.0*ct*ctm*n[1]*n[2]*gn[2]
   + ctm*st*(1-nz2) *n[2]*gn[0]
   + st*st*(1-nz2) *gn[1]
   + ctm*st*(1+nz2) *n[0]*gn[2]
   + 2.0*st*st*n[1]*n[2]*gn[2];
  DJDUb[2][3] =
   - st*(1 - nz2) * n[2]*gt
   + 2.0*st*(ct+ctm*nz2) *n[1]*gn[0]
   - 2.0*st*(ct+ctm*nz2) *n[0]*gn[1]
   - ct*ctm*(1-nz2)*gn[2]
   - st*st*(1-nz2)*gn[2];

  // d(J)/d(gn)
  DJDUXa[0][0] = -2*st;
  DJDUXa[0][1] = +2*ctm*n[2];
  DJDUXa[0][2] = -2*ctm*n[1];

  DJDUXa[1][0] = -2*ctm*n[2];
  DJDUXa[1][1] = -2*st;
  DJDUXa[1][2] = +2*ctm*n[0];

  DJDUXa[2][0] = +2*ctm*n[1];
  DJDUXa[2][1] = -2*ctm*n[0];
  DJDUXa[2][2] = -2*st;

  DJDUXb[0][0] = -st*(ct+ctm*nz2);
  DJDUXb[0][1] = +ctm*(ct+ctm*nz2)*n[2];
  DJDUXb[0][2] = +ctm*(ct-ctm*nz2)*n[1] + 2*st*ctm*n[0]*n[2];

  DJDUXb[1][0] = -ctm*(ct+ctm*nz2)*n[2];
  DJDUXb[1][1] = -st*(ct+ctm*nz2);
  DJDUXb[1][2] = -ctm*(ct-ctm*nz2)*n[0] + 2*st*ctm*n[1]*n[2];

  DJDUXb[2][0] = +(st*st+ctm*ctm*nz2)*n[1];
  DJDUXb[2][1] = -(st*st+ctm*ctm*nz2)*n[0];
  DJDUXb[2][2] = -st*ctm*(1-nz2);

  // d(J)/d(gt)
  DJDUXa[0][3] = -2*n[0];
  DJDUXa[1][3] = -2*n[1];
  DJDUXa[2][3] = -2*n[2];

  DJDUXb[0][3] = -(1-ctm*nz2)*n[0] + st*n[1]*n[2];
  DJDUXb[1][3] = -(1-ctm*nz2)*n[1] - st*n[0]*n[2];
  DJDUXb[2][3] = - ctm*(1 - nz2) * n[2];
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
  FOR(a) FOR(b) FOR(c) FOR(j) Ta[a] += -ee[a][b][c]*R[c][j]*ggR[b][j];
  FOR(a) FOR(b) FOR(c)        Tb[a] += -ee[a][b][c]*R[c][2]*ggR[b][2];
}

/***********************************************************/
// Calculate gradient torques Ta, Tb (v1, in n and th)
void fill_TG1_nt_(double Ta[DIM], double Tb[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt,
                 const double ggn[DIM], const double ggt) {
  double ct=cos(t), ctm=(1.0-ct), st=sin(t);
  int a,b,j,k,m;
  double ee[DIM][DIM][DIM];
  double dd[DIM][DIM];

  double en[DIM][DIM];
  fill_dd_(dd);
  fill_ee_(ee);
  fill_en_(en, n);
  fill_vec_zero_(Ta);
  fill_vec_zero_(Tb);


  FOR(a) FOR(j){
    Ta[a] +=
       + dd[j][0] *2*ctm   *gn[a]*gt
       + dd[j][0] *2*st    *ggn[a]
       + dd[j][0] *2*n[a]  *ggt
       -          2*st     *en[a][j]*gn[j]*gt
       +          2*ctm    *en[j][a]*ggn[j]
       -          2*st*ctm *n[a]*n[j]*ggn[j]
       -          2*ctm*st *n[a]*gn[j]*gn[j];
    Tb[a] +=
       + dd[j][0] *dd[a][2] *n[2]*st*gt*gt
       - dd[j][0] *dd[a][2] *n[2]*ct*ggt
       - dd[j][0] *dd[a][2] *2*(ct*ct-st*st)*gn[2]*gt
       - dd[j][0] *dd[a][2] *(2*ct-1)*st *ggn[2]
       - dd[j][0] *st*n[a]*n[2]*n[2] *gt*gt
       - dd[j][0] *2*st*ctm*n[a]*gn[2]*gn[2]
       + dd[j][0] *2*ct*ct *gn[a]*gt
       + dd[j][0] *2*ctm*ct* n[2]*n[2] *gn[a]*gt
       - dd[j][0] *4*st*st*  n[a]*n[2] *gn[2]*gt
       - dd[j][0] *2*ct*st*  en[2][a]  *gn[2]*gt
       + dd[j][0] *n[a]*ggt
       - dd[j][0] *ctm*n[a]*n[2]*n[2] * ggt
       - dd[j][0] *st*en[2][a]*n[2]*ggt
       + dd[j][0] *ct * st     *ggn[a]
       + dd[j][0] *ctm*st *n[2]*n[2] * ggn[a]
       - dd[j][0] *2*ctm*st *n[a]*n[2] * ggn[2]
       - dd[j][0] *ct*ctm *en[2][a] * ggn[2]
       - dd[a][2] *st*st *en[2][j]* ggn[j]
       - 2*ctm*st *en[a][j]*n[2]*n[2]*gn[j]*gt
       -  ct*ctm *ee[a][j][2] *n[2] * ggn[j]
       - ctm*ctm* en[a][j]*n[2]*n[2]*ggn[j];
     }
}

/***********************************************************/
// Calculate gradient torques Ta, Tb (v2, in conponents)
void fill_TG2_nt_(double Ta[DIM], double Tb[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt,
                 const double ggn[DIM], const double ggt) {
  double ct=cos(t), ctm=(1.0-ct), st=sin(t);
  int a,b,j,k,m;
  double nz2 = n[2]*n[2];
  double cp=ct+ctm*nz2, cm=ct-ctm*nz2;

  double xx=st*ctm*(n[0]*ggn[0]+gn[0]*gn[0] +
                    n[1]*ggn[1]+gn[1]*gn[1] +
                    n[2]*ggn[2]+gn[2]*gn[2]);
  double x0=st*gn[0]*gt+ctm*ggn[0];
  double x1=st*gn[1]*gt+ctm*ggn[1];
  double x2=st*gn[2]*gt+ctm*ggn[2];

  Ta[0] =  2*(ctm*gn[0]*gt + st*ggn[0] + n[0]*ggt
            - x1*n[2] + x2*n[1] - xx*n[0]);
  Ta[1] =  2*(ctm*gn[1]*gt + st*ggn[1] + n[1]*ggt
            + x0*n[2] - x2*n[0] - xx*n[1]);
  Ta[2] =  2*(ctm*gn[2]*gt + st*ggn[2] + n[2]*ggt
            - x0*n[1] + x1*n[0] - xx*n[2]);

  Tb[0] =
     + ((1.0- ctm*nz2)*n[0] - st*n[1]*n[2])*ggt
     - st*n[0]*nz2*gt*gt
     + 2*ct*cp*gn[0]*gt
     - 2*ctm*st*nz2*n[2]*gn[1]*gt
     - 2*st*cm*n[1]*gn[2]*gt
     - 4*st*st*n[0]*n[2]*gn[2]*gt
     - 2*st*ctm*n[0]*gn[2]*gn[2]
     +  st*cp*ggn[0]
     - ctm*cp*n[2]*ggn[1]
     - ctm*cm*n[1]*ggn[2]
     - 2*ctm*st*n[0]*n[2]*ggn[2];
  Tb[1] =
     + ((1.0- ctm*nz2)*n[1] + st*n[0]*n[2])*ggt
     - st*n[1]*nz2*gt*gt
     + 2*ctm*st*nz2*n[2]*gn[0]*gt
     + 2*ct*cp*gn[1]*gt
     + 2*st*cm*n[0]*gn[2]*gt
     - 4*st*st*  n[1]*n[2]*gn[2]*gt
     - 2*st*ctm*n[1]*gn[2]*gn[2]
     + ctm*cp*n[2]*ggn[0]
     + st*cp*ggn[1]
     + ctm*cm*n[0]*ggn[2]
     - 2*ctm*st*n[1]*n[2]*ggn[2];
  Tb[2] =
     + (1.0-nz2)*n[2]*(st*gt*gt + ctm*ggt)
     - 2*ctm*st*nz2*(n[1]*gn[0]-n[0]*gn[1])*gt
     - 2*(st*st - (ctm*ct+2*st*st)*nz2)*gn[2]*gt
     - 2*st*ctm*n[2]*gn[2]*gn[2]
     - (st*st+ctm*ctm*nz2) *(n[1]*ggn[0]-n[0]*ggn[1])
     + st*ctm*(1.0-nz2)*ggn[2];
}



// Calculate gradient torque TD = Ta/2+Tb (as in Dmitriev's program)
void fill_TGD_nt_(double T[DIM],
                 const double n[DIM], const double t,
                 const double gn[DIM], const double gt,
                 const double ggn[DIM], const double ggt) {
  double ct=cos(t), ctm=1.0-ct, ctp=1.0+ct, st=sin(t);

  double DD45=n[0]*gn[1]-n[1]*gn[0];
  double FTN=ctm*DD45 - st*gn[2] - gt*n[2];
  double DFTN=ctm*(n[0]*ggn[1]-ggn[0]*n[1])-st*ggn[2]-ggt*n[2]-
              ctp*gt*gn[2]+st*gt*DD45;

  T[0] = 2.0*(ggt*n[0]+ctp*gt*gn[0]+st*ggn[0]+
        st*gt*(n[1]*gn[2]-gn[1]*n[2])+ctm*(n[1]*ggn[2]-ggn[1]*n[2]))+
        (ctm*n[0]*n[2]+n[1]*st)*DFTN+(st*gt*n[0]*n[2]+
        ctm*(gn[0]*n[2]+n[0]*gn[2])+gn[1]*st+n[1]*ct*gt)*FTN;
  T[1] = 2.0*(ggt*n[1]+ctp*gt*gn[1]+st*ggn[1]-
        st*gt*(n[0]*gn[2]-gn[0]*n[2])-ctm*(n[0]*ggn[2]-ggn[0]*n[2]))+
        (ctm*n[1]*n[2]-n[0]*st)*DFTN+(st*gt*n[1]*n[2]+
        ctm*(gn[1]*n[2]+n[1]*gn[2])-gn[0]*st-n[0]*ct*gt)*FTN;
  T[2] = 2.0*(ggt*n[2]+ctp*gt*gn[2]+st*ggn[2]+
        st*gt*DD45+ctm*(n[0]*ggn[1]-ggn[0]*n[1]))+
        (ctm*n[2]*n[2]+ct)*DFTN+(st*gt*n[2]*n[2]+
        ctm*2.0*n[2]*gn[2]-st*gt)*FTN;
}

} // extern
