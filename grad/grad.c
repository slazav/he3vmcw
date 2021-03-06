extern "C" {
#include "grad.h"

/***********************************************************/
// Fill a rotation matrix using rotation axis (n) and
// rotation angle (theta)
void fill_r_nt_(double R[DIM][DIM], const double n[DIM], const double *th){
  double ct=cos(*th), st=sin(*th);
  double dd[DIM][DIM];
  double en[DIM][DIM];
  int i,j;
  fill_dd_(dd);
  fill_en_(en, n);
  FOR(i) FOR(j) R[i][j] = ct*dd[i][j] + (1.-ct)*n[i]*n[j] - st*en[i][j];
}

// calculate rotation matrix gradient in n and theta coordinates
void fill_gr_nt_(double gR[DIM][DIM],
                     const double n[DIM], const double *t,
                     const double gn[DIM], const double *gt) {
   double dd[DIM][DIM];
   double en[DIM][DIM];
   double eg[DIM][DIM];
   int a,j;
   fill_dd_(dd);
   fill_en_(en, n);
   fill_en_(eg, gn);
   double ct=cos(*t), ctm=(1.0-ct), st=sin(*t);

   FOR(a) FOR(j) gR[a][j] =
       st*(n[a]*n[j] - dd[a][j])*(*gt)
    + ctm*(n[j]*gn[a] + n[a]*gn[j])
    - ct*en[a][j]*(*gt) - st*eg[a][j];
}

// calculate rotation matrix double gradient in n and theta coordinates
void fill_ggr_nt_(double ggR[DIM][DIM],
                     const double n[DIM], const double *t,
                     const double gn[DIM], const double *gt,
                     const double ggn[DIM], const double *ggt) {
   double dd[DIM][DIM];
   double en[DIM][DIM];
   double eg[DIM][DIM];
   double egg[DIM][DIM];
   int a,j;
   fill_dd_(dd);
   fill_en_(en, n);
   fill_en_(eg, gn);
   fill_en_(egg, ggn);
   double ct=cos(*t), ctm=(1.0-ct), st=sin(*t);

   FOR(a) FOR(j) ggR[a][j] =
       + (n[a]*n[j]-dd[a][j])*(ct*(*gt)*(*gt) + st*(*ggt))
       + (n[a]*gn[j]+n[j]*gn[a])*2.0*st*(*gt)
       + ctm*(gn[j]*gn[a] + gn[a]*gn[j])
       + ctm*(n[j]*ggn[a] + n[a]*ggn[j])
       + en[a][j]*(st*(*gt)*(*gt) - ct*(*ggt))
       - 2.0*ct*eg[a][j]*(*gt) - st*egg[a][j];
}

/***********************************************************/
// Calculate gradient energies Ea, Eb (just by definition)
void fill_eg0_nt_(double *Ea, double *Eb,
                 const double n[DIM], const double *t,
                 const double gn[DIM], const double *gt) {
  int a,j;
  double gR[DIM][DIM];
  fill_gr_nt_(gR,n,t,gn,gt);
  *Ea=*Eb=0.0;

  FOR(a) FOR(j) *Ea += pow(gR[a][j],2); // K1
  FOR(a)        *Eb += pow(gR[a][2],2); // (K2+K3)
}

// v1, in n-th coordinates
void fill_eg1_nt_(double *Ea, double *Eb,
                 const double n[DIM], const double *t,
                 const double gn[DIM], const double *gt) {
  int a,j;
  double ct=cos(*t), ctm=(1.0-ct), st=sin(*t);
  double en[DIM][DIM];
  double eg[DIM][DIM];
  double dd[DIM][DIM];
  fill_en_(en, n);
  fill_en_(eg, gn);
  fill_dd_(dd);
  *Ea=*Eb=0.0;

  FOR(a) FOR(j) *Ea +=
      pow((ctm*(n[j]*gn[a] + n[a]*gn[j]) - st*eg[a][j])  +
      (st*(n[a]*n[j] - dd[a][j]) - ct*en[a][j]) * (*gt), 2);

  FOR(a) *Eb +=
      pow((ctm*(n[2]*gn[a] + n[a]*gn[2]) - st*eg[a][2]) +
      (st*(n[a]*n[2] - dd[a][2]) - ct*en[a][2]) * (*gt), 2);

}

/***********************************************************/
// Calculate spin currents Ja, Jb (just by definition)
void fill_jg0_nt_(double Ja[DIM], double Jb[DIM],
                 const double n[DIM], const double *t,
                 const double gn[DIM], const double *gt) {
  int a,b,c,j;
  double R[DIM][DIM], gR[DIM][DIM];
  double ee[DIM][DIM][DIM];
  fill_ee_(ee);
  fill_vec_zero_(Ja);
  fill_vec_zero_(Jb);
  fill_r_nt_(R,n,t);
  fill_gr_nt_(gR,n,t,gn, gt);
  FOR(a) FOR(b) FOR(c) FOR(j) Ja[a] += ee[a][b][c]*R[c][j]*gR[b][j];
  FOR(a) FOR(b) FOR(c)        Jb[a] += ee[a][b][c]*R[c][2]*gR[b][2];
}

// Calculate spin currents Ja, Jb (v1)
void fill_jg1_nt_(double Ja[DIM], double Jb[DIM],
                 const double n[DIM], const double *t,
                 const double gn[DIM], const double *gt) {
  int a,b,c,j;
  double ct=cos(*t), ctm=(1.0-ct), st=sin(*t);
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

  FOR(a) Ja[a] += -2*(n[a]*(*gt) + st*gn[a] + ctm*nxng[a]);

  FOR(a) Jb[a] +=
     - n[a]*(*gt)  *(1-ctm*n[2]*n[2])
     + st*ctm   *2*n[a]*n[2]*gn[2]
     + st       *en[2][a]*n[2]*(*gt)
     + ct*ctm   *(en[2][a]*gn[2] + eg[2][a]*n[2])
     - ctm*ctm  *nxng[a]*n[2]*n[2]
     - st*gn[a] *(ct + ctm*n[2]*n[2]);
    ;
  Jb[2] +=  ct*n[2]*(*gt) - st*(1-2*ct)*gn[2] - st*st*nxng[2];
}

/// Calculate spin currents Ja, Jb (v2 - in coordinates)
void fill_jg2_nt_(double Ja[DIM], double Jb[DIM],
                 const double n[DIM], const double *t,
                 const double gn[DIM], const double *gt) {
  double ct=cos(*t), ctm=(1.0-ct), st=sin(*t);
  double nz2 = n[2]*n[2];

  Ja[0] = -2*(n[0]*(*gt) + st*gn[0] + ctm *(n[1]*gn[2] - n[2]*gn[1]));
  Ja[1] = -2*(n[1]*(*gt) + st*gn[1] + ctm *(n[2]*gn[0] - n[0]*gn[2]));
  Ja[2] = -2*(n[2]*(*gt) + st*gn[2] + ctm *(n[0]*gn[1] - n[1]*gn[0]));

  Jb[0] =
   - (1-ctm*nz2)*n[0]*(*gt)
   + st*n[1]*n[2]*(*gt)
   -  st*(ct + ctm*nz2) *gn[0]
   + ctm*(ct + ctm*nz2) *n[2]*gn[1]
   + ctm*(ct - ctm*nz2) *n[1]*gn[2]
   + 2*st*ctm*n[0]*n[2]*gn[2];

  Jb[1] =
   - (1-ctm*nz2)*n[1]*(*gt)
   - st*n[0]*n[2]*(*gt)
   - ctm*(ct+ctm*nz2) *n[2]*gn[0]
   -  st*(ct+ctm*nz2) *gn[1]
   - ctm*(ct-ctm*nz2) *n[0]*gn[2]
   + 2*st*ctm*n[1]*n[2]*gn[2];

  Jb[2] =
   - ctm*(1 - nz2) * n[2]*(*gt)
   + (st*st + ctm*ctm *nz2) *n[1]*gn[0]
   - (st*st + ctm*ctm *nz2) *n[0]*gn[1]
   - st*ctm*(1-nz2)*gn[2];
}

/// Calculate spin current J = Ja/2 + Jb (from Dmitriev's program)
void fill_jgd_nt_(double J[DIM],
                 const double n[DIM], const double *t,
                 const double gn[DIM], const double *gt) {
  double ct=cos(*t), ctm=(1.0-ct), st=sin(*t);
  double FTN=ctm*(n[0]*gn[1]-n[1]*gn[0]) - st*gn[2] - (*gt)*n[2];
  J[0] = -2.0*((*gt)*n[0]+st*gn[0]+ctm*(n[1]*gn[2]-gn[1]*n[2]))
       - (ctm*n[0]*n[2]+n[1]*st)*FTN;
  J[1] = -2.0*((*gt)*n[1]+st*gn[1]-ctm*(n[0]*gn[2]-gn[0]*n[2]))
       - (ctm*n[1]*n[2]-n[0]*st)*FTN;
  J[2] = -2.0*((*gt)*n[2]+st*gn[2]+ctm*(n[0]*gn[1]-gn[0]*n[1]))
       - (ctm*n[2]*n[2]+ct)*FTN;
}

/// Spin current derivative, D(J)/D(U), D(J)/D(U') where U = nx,ny,nz,th
void fill_dj_nt_(double DJDUa[DIM][4],  double DJDUb[DIM][4],
                 double DJDUZa[DIM][4], double DJDUZb[DIM][4],
                 const double n[DIM], const double *t,
                 const double gn[DIM], const double *gt){

  double ct=cos(*t), ctm=(1.0-ct), st=sin(*t);
  double nz2 = n[2]*n[2];

  // d(J)/d(n)
  DJDUa[0][0] = -2*(*gt);
  DJDUa[0][1] = -2*ctm*gn[2];
  DJDUa[0][2] = +2*ctm*gn[1];

  DJDUa[1][0] = +2*ctm*gn[2];
  DJDUa[1][1] = -2*(*gt);
  DJDUa[1][2] = -2*ctm*gn[0];

  DJDUa[2][0] = -2*ctm*gn[1];
  DJDUa[2][1] = +2*ctm*gn[0];
  DJDUa[2][2] = -2*(*gt);

  DJDUb[0][0] = - (1-ctm*nz2)*(*gt) + 2*st*ctm*n[2]*gn[2];
  DJDUb[0][1] = + st*n[2]*(*gt) + ctm*(ct-ctm*nz2)*gn[2];
  DJDUb[0][2] =
   + 2.0*ctm*n[2]*n[0]*(*gt) + st*n[1]*(*gt)
   - 2.0*st*ctm*n[2]*gn[0]
   + ctm*(ct + 3.0*ctm*n[2]*n[2])*gn[1]
   - 2.0*ctm*ctm*n[2]*n[1]*gn[2]
   + 2.0*st*ctm*n[0]*gn[2];

  DJDUb[1][0] = - st*n[2]*(*gt) - ctm*(ct-ctm*nz2)*gn[2];
  DJDUb[1][1] = - (1-ctm*nz2)*(*gt) + 2*st*ctm*n[2]*gn[2];
  DJDUb[1][2] =
   + 2.0*ctm*n[2]*n[1]*(*gt) - st*n[0]*(*gt)
   - ctm*(ct + 3.0*ctm*n[2]*n[2])*gn[0]
   - 2.0*st*ctm*n[2] *gn[1]
   + 2.0*ctm*ctm*n[2]*n[0]*gn[2]
   + 2.0*st*ctm*n[1]*gn[2];

  DJDUb[2][0] = - (st*st + ctm*ctm *nz2)*gn[1];
  DJDUb[2][1] = + (st*st + ctm*ctm *nz2)*gn[0];
  DJDUb[2][2] = ctm*(
   - (1 - 3*n[2]*n[2])*(*gt)
   + 2.0*ctm*n[1]*n[2]*gn[0]
   - 2.0*ctm*n[0]*n[2]*gn[1]
   + 2.0*st*n[2]*gn[2]
  );

  // d(J)/d(t)
  DJDUa[0][3] = -2*(ct*gn[0] + st*(n[1]*gn[2] - n[2]*gn[1]));
  DJDUa[1][3] = -2*(ct*gn[1] + st*(n[2]*gn[0] - n[0]*gn[2]));
  DJDUa[2][3] = -2*(ct*gn[2] + st*(n[0]*gn[1] - n[1]*gn[0]));

  DJDUb[0][3] =
     (st*nz2*n[0] + ct*n[1]*n[2])*(*gt)
   - (ct*(ct+ctm*nz2)-st*st*(1-nz2))  *gn[0]
   + (st*(ct+ctm*nz2)-ctm*st*(1-nz2)) *n[2]*gn[1]
   + (st*(ct-ctm*nz2)-ctm*st*(1+nz2)) *n[1]*gn[2]
   + 2.0*(ct*ctm + st*st)  *n[0]*n[2]*gn[2];
  DJDUb[1][3] =
     (st*nz2*n[1] - ct*n[0]*n[2])*(*gt)
   - (st*(ct+ctm*nz2)-ctm*st*(1-nz2)) *n[2]*gn[0]
   - (ct*(ct+ctm*nz2)-st*st*(1-nz2))  *gn[1]
   - (st*(ct-ctm*nz2)-ctm*st*(1+nz2)) *n[0]*gn[2]
   + 2.0*(ct*ctm+st*st) *n[1]*n[2]*gn[2];
  DJDUb[2][3] =
   - st*(1 - nz2) * n[2]*(*gt)
   + 2.0*st*(ct+ctm*nz2) *n[1]*gn[0]
   - 2.0*st*(ct+ctm*nz2) *n[0]*gn[1]
   - (ct*ctm+st*st)*(1-nz2)*gn[2];

  // d(J)/d(gn)
  DJDUZa[0][0] = -2*st;
  DJDUZa[0][1] = +2*ctm*n[2];
  DJDUZa[0][2] = -2*ctm*n[1];

  DJDUZa[1][0] = -2*ctm*n[2];
  DJDUZa[1][1] = -2*st;
  DJDUZa[1][2] = +2*ctm*n[0];

  DJDUZa[2][0] = +2*ctm*n[1];
  DJDUZa[2][1] = -2*ctm*n[0];
  DJDUZa[2][2] = -2*st;

  DJDUZb[0][0] = -st*(ct+ctm*nz2);
  DJDUZb[0][1] = +ctm*(ct+ctm*nz2)*n[2];
  DJDUZb[0][2] = +ctm*(ct-ctm*nz2)*n[1] + 2*st*ctm*n[0]*n[2];

  DJDUZb[1][0] = -ctm*(ct+ctm*nz2)*n[2];
  DJDUZb[1][1] = -st*(ct+ctm*nz2);
  DJDUZb[1][2] = -ctm*(ct-ctm*nz2)*n[0] + 2*st*ctm*n[1]*n[2];

  DJDUZb[2][0] = +(st*st+ctm*ctm*nz2)*n[1];
  DJDUZb[2][1] = -(st*st+ctm*ctm*nz2)*n[0];
  DJDUZb[2][2] = -st*ctm*(1-nz2);

  // d(J)/d(gt)
  DJDUZa[0][3] = -2*n[0];
  DJDUZa[1][3] = -2*n[1];
  DJDUZa[2][3] = -2*n[2];

  DJDUZb[0][3] = -(1-ctm*nz2)*n[0] + st*n[1]*n[2];
  DJDUZb[1][3] = -(1-ctm*nz2)*n[1] - st*n[0]*n[2];
  DJDUZb[2][3] = - ctm*(1 - nz2) * n[2];
}

/// Spin current derivative, D(J)/D(U), D(J)/D(U') where U = nx,ny,nz,th
/// From Dmitriev's program
void fill_djd_nt_(double DJDU[DIM][4], double DJDUZ[DIM][4],
                 const double n[DIM], const double *t,
                 const double gn[DIM], const double *gt){
  int i,j;
  double ct=cos(*t), ctm=(1.0-ct), st=sin(*t);
  double nz2 = n[2]*n[2];

  double DD45=n[0]*gn[1]-gn[0]*n[1];
  double FTN=ctm*DD45-st*gn[2]-(*gt)*n[2];
  double CTF=ctm*FTN;
  double STF=st*FTN;
  double FTN4=ctm*gn[1];
  double FTN5=-ctm*gn[0];
  double FTN7=st*DD45-ct*gn[2];
  double FTNX4=-ctm*n[1];
  double FTNX5=ctm*n[0];
  double C46=ctm*n[0]*n[2]+n[1]*st;
  double C56=ctm*n[1]*n[2]-n[0]*st;
  double C66=ctm*n[2]*n[2]+ct;
  double C266=2.0-C66;

  DJDU[0][0] =  2.0*(*gt)+CTF*n[2]+C46*FTN4;
  DJDU[0][1] =  2.0*ctm*gn[2]+STF+C46*FTN5;
  DJDU[0][2] = -2.0*ctm*gn[1]+CTF*n[0]-C46*(*gt);
  DJDU[0][3] =  2.0*(ct*gn[0]+st*(n[1]*gn[2]-gn[1]*n[2]))+STF*n[0]*n[2]+n[1]*ct*FTN+C46*FTN7;
  DJDUZ[0][0] = 2.0*st+C46*FTNX4;
  DJDUZ[0][1] =-2.0*ctm*n[2]+C46*FTNX5;
  DJDUZ[0][2] = 2.0*ctm*n[1]-C46*st;
  DJDUZ[0][3] = 2.0*n[0]-C46*n[2];

  DJDU[1][0] =  -2.0*ctm*gn[2]-STF+C56*FTN4;
  DJDU[1][1] =   2.0*(*gt)+CTF*n[2]+C56*FTN5;
  DJDU[1][2] =   2.0*ctm*gn[0]+CTF*n[1]-C56*(*gt);
  DJDU[1][3] =   2.0*(ct*gn[1]-st*(n[0]*gn[2]-gn[0]*n[2]))+STF*n[1]*n[2]-n[0]*ct*FTN+C56*FTN7;
  DJDUZ[1][0] =  2.0*ctm*n[2]+C56*FTNX4;
  DJDUZ[1][1] =  2.0*st+C56*FTNX5;
  DJDUZ[1][2] = -2.0*ctm*n[0]-C56*st;
  DJDUZ[1][3] =  2.0*n[1]-C56*n[2];

  DJDU[2][0] =   2.0*ctm*gn[1]+C66*FTN4;
  DJDU[2][1] =  -2.0*ctm*gn[0]+C66*FTN5;
  DJDU[2][2] =   2.0*n[2]*CTF+C266*(*gt);
  DJDU[2][3] =   2.0*(ct*gn[2]+st*DD45)+ STF*(nz2-1.0)+C66*FTN7;
  DJDUZ[2][0] = -2.0*ctm*n[1]+C66*FTNX4;
  DJDUZ[2][1] =  2.0*ctm*n[0]+C66*FTNX5;
  DJDUZ[2][2] =  C266*st;
  DJDUZ[2][3] =  C266*n[2];

  // invert sign
  FOR(i) for (j=0;j<4; j++){
    DJDU[i][j] = -DJDU[i][j];
    DJDUZ[i][j] = -DJDUZ[i][j];
  }

}

/***********************************************************/
// Calculate gradient torques Ta, Tb (just by definition)
void fill_tg0_nt_(double Ta[DIM], double Tb[DIM],
                 const double n[DIM], const double *t,
                 const double gn[DIM], const double *gt,
                 const double ggn[DIM], const double *ggt) {
  int a,b,c,j;
  double R[DIM][DIM], ggR[DIM][DIM];
  double ee[DIM][DIM][DIM];
  fill_ee_(ee);
  fill_vec_zero_(Ta);
  fill_vec_zero_(Tb);
  fill_r_nt_(R,n,t);
  fill_ggr_nt_(ggR, n,t, gn, gt, ggn, ggt);
  FOR(a) FOR(b) FOR(c) FOR(j) Ta[a] += -ee[a][b][c]*R[c][j]*ggR[b][j];
  FOR(a) FOR(b) FOR(c)        Tb[a] += -ee[a][b][c]*R[c][2]*ggR[b][2];
}

/***********************************************************/
// Calculate gradient torques Ta, Tb (v1, in n and th)
void fill_tg1_nt_(double Ta[DIM], double Tb[DIM],
                 const double n[DIM], const double *t,
                 const double gn[DIM], const double *gt,
                 const double ggn[DIM], const double *ggt) {
  double ct=cos(*t), ctm=(1.0-ct), st=sin(*t);
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
       + dd[j][0] *2*ctm   *gn[a]*(*gt)
       + dd[j][0] *2*st    *ggn[a]
       + dd[j][0] *2*n[a]  *(*ggt)
       -          2*st     *en[a][j]*gn[j]*(*gt)
       +          2*ctm    *en[j][a]*ggn[j]
       -          2*st*ctm *n[a]*n[j]*ggn[j]
       -          2*ctm*st *n[a]*gn[j]*gn[j];
    Tb[a] +=
       + dd[j][0] *dd[a][2] *n[2]*st*(*gt)*(*gt)
       - dd[j][0] *dd[a][2] *n[2]*ct*(*ggt)
       - dd[j][0] *dd[a][2] *2*(ct*ct-st*st)*gn[2]*(*gt)
       - dd[j][0] *dd[a][2] *(2*ct-1)*st *ggn[2]
       - dd[j][0] *st*n[a]*n[2]*n[2] *(*gt)*(*gt)
       - dd[j][0] *2*st*ctm*n[a]*gn[2]*gn[2]
       + dd[j][0] *2*ct*ct *gn[a]*(*gt)
       + dd[j][0] *2*ctm*ct* n[2]*n[2] *gn[a]*(*gt)
       - dd[j][0] *4*st*st*  n[a]*n[2] *gn[2]*(*gt)
       - dd[j][0] *2*ct*st*  en[2][a]  *gn[2]*(*gt)
       + dd[j][0] *n[a]*(*ggt)
       - dd[j][0] *ctm*n[a]*n[2]*n[2] * (*ggt)
       - dd[j][0] *st*en[2][a]*n[2]*(*ggt)
       + dd[j][0] *ct * st     *ggn[a]
       + dd[j][0] *ctm*st *n[2]*n[2] * ggn[a]
       - dd[j][0] *2*ctm*st *n[a]*n[2] * ggn[2]
       - dd[j][0] *ct*ctm *en[2][a] * ggn[2]
       - dd[a][2] *st*st *en[2][j]* ggn[j]
       - 2*ctm*st *en[a][j]*n[2]*n[2]*gn[j]*(*gt)
       -  ct*ctm *ee[a][j][2] *n[2] * ggn[j]
       - ctm*ctm* en[a][j]*n[2]*n[2]*ggn[j];
     }
}

/// Spin current derivative, D(J)/D(U), D(J)/D(U') where U = tnx,tny,tnz
void fill_dj_t_(double DJDUa[DIM][DIM],  double DJDUb[DIM][DIM],
                double DJDUZa[DIM][DIM], double DJDUZb[DIM][DIM],
                const double n[DIM], const double *t,
                const double gn[DIM], const double *gt){
  int i,j,k;
  double ct=cos(*t), ctm=(1.0-ct), st=sin(*t);
  double nz2 = n[2]*n[2];
  double Da[DIM][4], Db[DIM][4], DZa[DIM][4], DZb[DIM][4];
  double dd[DIM][DIM];
  fill_dd_(dd);

  // calculate derivative in n,t coordinates
  fill_dj_nt_(Da,Db,DZa,DZb, n,t,gn,gt);
  FOR(i){
    // same conversion for 'a' and 'b' parts
    for (k=0;k<2; k++){
      // new derivatives (pointers)
      double *NDx  = ((k==0)? DJDUa[i]:DJDUb[i]);
      double *NDgx = ((k==0)? DJDUZa[i]:DJDUZb[i]);
      // old derivatives
      double Dnx  = (k==0)? Da[i][0]:Db[i][0];
      double Dny  = (k==0)? Da[i][1]:Db[i][1];
      double Dnz  = (k==0)? Da[i][2]:Db[i][2];
      double Dt   = (k==0)? Da[i][3]:Db[i][3];
      double Dgnx = (k==0)? DZa[i][0]:DZb[i][0];
      double Dgny = (k==0)? DZa[i][1]:DZb[i][1];
      double Dgnz = (k==0)? DZa[i][2]:DZb[i][2];
      double Dgt  = (k==0)? DZa[i][3]:DZb[i][3];

      NDx[0] = Dt*n[0] + Dgt*gn[0]
      + Dnx/(*t) - n[0]*(Dnx*n[0]+Dny*n[1]+Dnz*n[2])/(*t)
      + Dgnx*((*gt)/(*t)*(n[0]*n[0]-1.0)-2.0*n[0]*gn[0])/(*t)
      + Dgny*((*gt)/(*t)*(n[0]*n[1])-n[0]*gn[1]-n[1]*gn[0])/(*t)
      + Dgnz*((*gt)/(*t)*(n[0]*n[2])-n[0]*gn[2]-n[2]*gn[0])/(*t);
      NDx[1] = Dt*n[1] + Dgt*gn[1]
      + Dny/(*t) - n[1]*(Dnx*n[0]+Dny*n[1]+Dnz*n[2])/(*t)
      + Dgnx*((*gt)/(*t)*(n[1]*n[0])-n[1]*gn[0]-n[0]*gn[1])/(*t)
      + Dgny*((*gt)/(*t)*(n[1]*n[1]-1.0)-2.0*n[1]*gn[1])/(*t)
      + Dgnz*((*gt)/(*t)*(n[1]*n[2])-n[1]*gn[2]-n[2]*gn[1])/(*t);
      NDx[2] = Dt*n[2] + Dgt*gn[2]
      + Dnz/(*t) - n[2]*(Dnx*n[0]+Dny*n[1]+Dnz*n[2])/(*t)
      + Dgnx*((*gt)/(*t)*(n[2]*n[0])-n[2]*gn[0]-n[0]*gn[2])/(*t)
      + Dgny*((*gt)/(*t)*(n[2]*n[1])-n[2]*gn[1]-n[1]*gn[2])/(*t)
      + Dgnz*((*gt)/(*t)*(n[2]*n[2]-1.0)-2.0*n[2]*gn[2])/(*t);

      NDgx[0] = Dgt*n[0] + Dgnx/(*t) - n[0]*(Dgnx*n[0]+Dgny*n[1]+Dgnz*n[2])/(*t);
      NDgx[1] = Dgt*n[1] + Dgny/(*t) - n[1]*(Dgnx*n[0]+Dgny*n[1]+Dgnz*n[2])/(*t);
      NDgx[2] = Dgt*n[2] + Dgnz/(*t) - n[2]*(Dgnx*n[0]+Dgny*n[1]+Dgnz*n[2])/(*t);
    }
  }
}

/***********************************************************/
// Calculate gradient torques Ta, Tb (v2, in conponents)
void fill_tg2_nt_(double Ta[DIM], double Tb[DIM],
                 const double n[DIM], const double *t,
                 const double gn[DIM], const double *gt,
                 const double ggn[DIM], const double *ggt) {
  double ct=cos(*t), ctm=(1.0-ct), st=sin(*t);
  int a,b,j,k,m;
  double nz2 = n[2]*n[2];
  double cp=ct+ctm*nz2, cm=ct-ctm*nz2;

  double xx=st*ctm*(n[0]*ggn[0]+gn[0]*gn[0] +
                    n[1]*ggn[1]+gn[1]*gn[1] +
                    n[2]*ggn[2]+gn[2]*gn[2]);
  double x0=st*gn[0]*(*gt)+ctm*ggn[0];
  double x1=st*gn[1]*(*gt)+ctm*ggn[1];
  double x2=st*gn[2]*(*gt)+ctm*ggn[2];

  Ta[0] =  2*(ctm*gn[0]*(*gt) + st*ggn[0] + n[0]*(*ggt)
            - x1*n[2] + x2*n[1] - xx*n[0]);
  Ta[1] =  2*(ctm*gn[1]*(*gt) + st*ggn[1] + n[1]*(*ggt)
            + x0*n[2] - x2*n[0] - xx*n[1]);
  Ta[2] =  2*(ctm*gn[2]*(*gt) + st*ggn[2] + n[2]*(*ggt)
            - x0*n[1] + x1*n[0] - xx*n[2]);

  Tb[0] =
     + ((1.0- ctm*nz2)*n[0] - st*n[1]*n[2])*(*ggt)
     - st*n[0]*nz2*(*gt)*(*gt)
     + 2*ct*cp*gn[0]*(*gt)
     - 2*ctm*st*nz2*n[2]*gn[1]*(*gt)
     - 2*st*cm*n[1]*gn[2]*(*gt)
     - 4*st*st*n[0]*n[2]*gn[2]*(*gt)
     - 2*st*ctm*n[0]*gn[2]*gn[2]
     +  st*cp*ggn[0]
     - ctm*cp*n[2]*ggn[1]
     - ctm*cm*n[1]*ggn[2]
     - 2*ctm*st*n[0]*n[2]*ggn[2];
  Tb[1] =
     + ((1.0- ctm*nz2)*n[1] + st*n[0]*n[2])*(*ggt)
     - st*n[1]*nz2*(*gt)*(*gt)
     + 2*ctm*st*nz2*n[2]*gn[0]*(*gt)
     + 2*ct*cp*gn[1]*(*gt)
     + 2*st*cm*n[0]*gn[2]*(*gt)
     - 4*st*st*  n[1]*n[2]*gn[2]*(*gt)
     - 2*st*ctm*n[1]*gn[2]*gn[2]
     + ctm*cp*n[2]*ggn[0]
     + st*cp*ggn[1]
     + ctm*cm*n[0]*ggn[2]
     - 2*ctm*st*n[1]*n[2]*ggn[2];
  Tb[2] =
     + (1.0-nz2)*n[2]*(st*(*gt)*(*gt) + ctm*(*ggt))
     - 2*ctm*st*nz2*(n[1]*gn[0]-n[0]*gn[1])*(*gt)
     - 2*(st*st - (ctm*ct+2*st*st)*nz2)*gn[2]*(*gt)
     - 2*st*ctm*n[2]*gn[2]*gn[2]
     - (st*st+ctm*ctm*nz2) *(n[1]*ggn[0]-n[0]*ggn[1])
     + st*ctm*(1.0-nz2)*ggn[2];
}



// Calculate gradient torque TD = Ta/2+Tb (from Dmitriev's program)
void fill_tgd_nt_(double T[DIM],
                 const double n[DIM], const double *t,
                 const double gn[DIM], const double *gt,
                 const double ggn[DIM], const double *ggt) {
  double ct=cos(*t), ctm=1.0-ct, ctp=1.0+ct, st=sin(*t);

  double DD45=n[0]*gn[1]-n[1]*gn[0];
  double FTN=ctm*DD45 - st*gn[2] - (*gt)*n[2];
  double DFTN=ctm*(n[0]*ggn[1]-ggn[0]*n[1])-st*ggn[2]-(*ggt)*n[2]-
              ctp*(*gt)*gn[2]+st*(*gt)*DD45;

  T[0] = 2.0*((*ggt)*n[0]+ctp*(*gt)*gn[0]+st*ggn[0]+
        st*(*gt)*(n[1]*gn[2]-gn[1]*n[2])+ctm*(n[1]*ggn[2]-ggn[1]*n[2]))+
        (ctm*n[0]*n[2]+n[1]*st)*DFTN+(st*(*gt)*n[0]*n[2]+
        ctm*(gn[0]*n[2]+n[0]*gn[2])+gn[1]*st+n[1]*ct*(*gt))*FTN;
  T[1] = 2.0*((*ggt)*n[1]+ctp*(*gt)*gn[1]+st*ggn[1]-
        st*(*gt)*(n[0]*gn[2]-gn[0]*n[2])-ctm*(n[0]*ggn[2]-ggn[0]*n[2]))+
        (ctm*n[1]*n[2]-n[0]*st)*DFTN+(st*(*gt)*n[1]*n[2]+
        ctm*(gn[1]*n[2]+n[1]*gn[2])-gn[0]*st-n[0]*ct*(*gt))*FTN;
  T[2] = 2.0*((*ggt)*n[2]+ctp*(*gt)*gn[2]+st*ggn[2]+
        st*(*gt)*DD45+ctm*(n[0]*ggn[1]-ggn[0]*n[1]))+
        (ctm*n[2]*n[2]+ct)*DFTN+(st*(*gt)*n[2]*n[2]+
        ctm*2.0*n[2]*gn[2]-st*(*gt))*FTN;
}


/***********************************************************/
/// Calculate components of L=S*R (just by definition)
void fill_l0_nt_(double L[DIM], const double S[DIM], const double n[DIM], const double *t){
  int a,j;
  double R[DIM][DIM];
  fill_r_nt_(R,n,t);
  fill_vec_zero_(L);

  FOR(a) FOR(j) L[j] += R[a][j]*S[a];

}

/// Calculate components of L=S*R
void fill_l1_nt_(double L[DIM], const double S[DIM], const double n[DIM], const double *t){
  int a,j;
  double nxS[DIM];
  double ct=cos(*t), ctm=(1.0-ct), st=sin(*t);
  double nS = n[0]*S[0]+n[1]*S[1]+n[2]*S[2];
  fill_vxv_(nxS, n, S);

  FOR(j) L[j] = ct*S[j] + ctm*n[j]*nS - st*nxS[j];
}


/// Calculate derivatives dL/dS,dL/dN,dL/dTh
void fill_dl_nt_(double DL[DIM][7], const double S[DIM], const double n[DIM], const double *t){
  int a,j;
  double nxS[DIM];
  double ct=cos(*t), ctm=(1.0-ct), st=sin(*t);
  double nS = n[0]*S[0]+n[1]*S[1]+n[2]*S[2];
  double dd[DIM][DIM];
  double eS[DIM][DIM];
  double en[DIM][DIM];
  fill_dd_(dd);
  fill_en_(eS, S);
  fill_en_(en, n);
  fill_vxv_(nxS, n, S);

  // dL_j/dS_a = R[a][j]
  FOR(a) FOR(j) DL[j][a] = ct*dd[a][j] + (1.-ct)*n[a]*n[j] - st*en[a][j];

  // dL_j/dn_a
  FOR(a) FOR(j) DL[j][a+3] = ctm*(S[a]*n[j] + nS*dd[j][a]) - st*eS[j][a];

  // dL_j/dth
  FOR(j) DL[j][6] = st*(n[j]*nS - S[j]) - ct*nxS[j];
}

/// Calculate derivatives dL/dS,dL/dNTh
void fill_dl_t_(double DL[DIM][6], const double S[DIM], const double n[DIM], const double *t){
  int a,j;
  double nxS[DIM];
  double ct=cos(*t), ctm=(1.0-ct), st=sin(*t);
  double nS = n[0]*S[0]+n[1]*S[1]+n[2]*S[2];
  double dd[DIM][DIM];
  double eS[DIM][DIM];
  double en[DIM][DIM];
  fill_dd_(dd);
  fill_en_(eS, S);
  fill_en_(en, n);
  fill_vxv_(nxS, n, S);

  // dL_j/dS_a = R[a][j]
  FOR(a) FOR(j) DL[j][a] = ct*dd[a][j] + (1.-ct)*n[a]*n[j] - st*en[a][j];

  // dL_j/dnth_a = dL_j/dn_b * dn_b/dnth_a + dL_j/dth * dth/dnth_a
  //             = dL_j/dn_b * (d[a][b] - n[a]*n[b])/th + dL_j/dth * n[a]
  FOR(a) FOR(j) DL[j][a+3] =
      (ctm*(S[a]*n[j] + nS*dd[j][a] - 2.0*nS*n[j]*n[a]) - st*eS[j][a])/(*t)
     + st*(n[j]*nS - S[j])*n[a] + st*nxS[j]*n[a]/(*t);


  FOR(a) FOR(j) DL[j][a+3] =
      ctm/(*t)*(S[a]*n[j] + nS*dd[j][a] - 2.0*nS*n[j]*n[a])
     - st/(*t)*(eS[j][a] - nxS[j]*n[a])
     + (st*(n[j]*nS - S[j]) - ct*nxS[j])*n[a];
}


} // extern
