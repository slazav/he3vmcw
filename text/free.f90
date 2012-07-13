MODULE energies

  USE general

  USE text

  USE glob

  IMPLICIT NONE

  REAL (KIND=dp), SAVE :: dx
  REAL (KIND=dp), SAVE :: nub,nu0
  REAL (KIND=dp), DIMENSION(0:maxnpt), SAVE :: apsi

  REAL (KIND=dp) :: lsg=3._dp ! see Fig. 1 in Erkki's paper

  CONTAINS


    SUBROUTINE sfun(n,x,f,g)
      IMPLICIT NONE
      INTEGER :: i,n
      REAL (KIND=dp) :: f
      REAL (KIND=dp), DIMENSION(n) :: x,g
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta,ga,gb
      DO i=1,nmax
         alpha(i)=x(i+1)
         beta(i)=x(i+nmax+1)
      END DO
      alpha(0)=x(1)
      beta(0)=0._dp
      f=energy(alpha,beta)
      CALL egrad(alpha,beta,ga,gb)
      DO i=1,nmax
         g(i+1)=ga(i)
         g(i+nmax+1)=gb(i)
      END DO      
      g(1)=ga(0)
    END SUBROUTINE sfun


    FUNCTION energy(alpha,beta) RESULT(e)
      ! Calculates the textural free energy
      IMPLICIT NONE
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta      
      REAL (KIND=dp) :: e
      e=emagn(beta)
      e=e+spinorbit(beta)
      e=e+esurf(alpha(nmax),beta(nmax))
      e=e+ebend(alpha,beta)
      e=e+eflow(alpha,beta)
      e=e+evortex(alpha,beta)
    END FUNCTION energy


    FUNCTION emagn(beta) RESULT(e)
      ! Calculates the magnetic free energy
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: beta      
      REAL (KIND=dp) :: sq,rp,rm,bp,bm,e
      sq=SQRT(3._dp)
      DO i=0,nmax-1
         rp=(i+(3+sq)/6)*dx
         rm=(i+(3-sq)/6)*dx
         bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
         bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6
         e=e+0.5*dx*(rp*SIN(bp)**2+rm*SIN(bm)**2)
      END DO
    END FUNCTION emagn


    FUNCTION spinorbit(beta) RESULT(e)
      !calculates the spin-orbit free energy
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: beta      
      REAL (KIND=dp) :: sq,rp,rm,bp,bm,e,chia,blp,blm,apsip,apsim
      sq=SQRT(3._dp)
      chia=fchia(t,p)
      e=0.0_dp
      DO i=0,nmax-1
         rp=(i+(3+sq)/6)*dx
         rm=(i+(3-sq)/6)*dx
         bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
         bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6
        ! blp=ACOS(-1/4+5/4*(COS(bp))**2)
        ! blm=ACOS(-1/4+5/4*(COS(bm))**2)
         apsip=(3+sq)*apsi(i+1)/6+(3-sq)*apsi(i)/6
         apsim=(3-sq)*apsi(i+1)/6+(3+sq)*apsi(i)/6
         e=e+chia*dx*((nub/nu0)**2)*(0.5*(apsip**2)*rp*SIN(bp)**2+0.5*(apsim**2)*rm*SIN(bm)**2)
      END DO
    END FUNCTION spinorbit


    FUNCTION esurf(alpha,beta) RESULT(e)
      ! Calculates the surface free energy
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp) :: alpha,beta,dar,nr,nf,nz,e
      dar=fdar(t,p,r)
      nr=-SIN(beta)*COS(alpha)
      nf=SIN(beta)*SIN(alpha)
      nz=COS(beta)
      e=-5*dar*(SQRT(5.)*nz*nr-SQRT(3.)*nf)**2/16
    END FUNCTION esurf    

    
    FUNCTION eflow(alpha,beta) RESULT(e)
      ! Calculates the flow free energy
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta      
      REAL (KIND=dp) :: s,c,sq,ap,am,rp,rm,bp,bm,e
      REAL (KIND=dp) :: vd,nr,nf,nz,rzr,rzf,rzz
      REAL (KIND=dp) :: vrp,vfp,vzp,vrm,vfm,vzm
      sq=SQRT(3._dp)
      c=-0.25_dp
      s=SQRT(15.)/4.0_dp
      e=0.0_dp
      vd=fvd(t,p)
      DO i=0,nmax-1
         vrp=evrp(i)
         vfp=evfp(i)
         vzp=evzp(i)
         rp=(i+(3+sq)/6)*dx
         ap=(3+sq)*alpha(i+1)/6+(3-sq)*alpha(i)/6
         bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
         nr=-SIN(bp)*COS(ap)
         nf=SIN(bp)*SIN(ap)
         nz=COS(bp)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         e=e-dx*rp*(rzr*vrp+rzf*vfp+rzz*vzp)**2/(5*vd**2)
         vrm=evrm(i)
         vfm=evfm(i)
         vzm=evzm(i)
         rm=(i+(3-sq)/6)*dx
         am=(3-sq)*alpha(i+1)/6+(3+sq)*alpha(i)/6
         bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6
         nr=-SIN(bm)*COS(am)
         nf=SIN(bm)*SIN(am)
         nz=COS(bm)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         e=e-dx*rm*(rzr*vrm+rzf*vfm+rzz*vzm)**2/(5*vd**2)
      END DO
    END FUNCTION eflow


    FUNCTION evortex(alpha,beta) RESULT(e)
      ! Calculates the vortex free energy
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta      
      REAL (KIND=dp) :: s,c,sq,ap,am,rp,rm,bp,bm,e
      REAL (KIND=dp) :: wm,wp,nr,nf,nz,rzr,rzf,rzz
      REAL (KIND=dp) :: lrp,lfp,lzp,lrm,lfm,lzm
      sq=SQRT(3._dp)
      c=-0.25_dp
      s=SQRT(15.)/4.0_dp
      e=0.0_dp
      DO i=0,nmax-1
         lrp=elrp(i)
         lfp=elfp(i)
         lzp=elzp(i)
         wp=ewp(i)
         rp=(i+(3+sq)/6)*dx
         ap=(3+sq)*alpha(i+1)/6+(3-sq)*alpha(i)/6
         bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
         nr=-SIN(bp)*COS(ap)
         nf=SIN(bp)*SIN(ap)
         nz=COS(bp)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         e=e+dx*lo*rp*wp*(rzr*lrp+rzf*lfp+rzz*lzp)**2/10
         lrm=elrm(i)
         lfm=elfm(i)
         lzm=elzm(i)
         wm=ewm(i)
         rm=(i+(3-sq)/6)*dx
         am=(3-sq)*alpha(i+1)/6+(3+sq)*alpha(i)/6
         bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6
         nr=-SIN(bm)*COS(am)
         nf=SIN(bm)*SIN(am)
         nz=COS(bm)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         e=e+dx*lo*rm*wm*(rzr*lrm+rzf*lfm+rzz*lzm)**2/10
      END DO
    END FUNCTION evortex


    FUNCTION ebend(alpha,beta) RESULT(e)
      ! Calculates the bending free energy
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta      
      REAL (KIND=dp) :: da,db,e,con,help,xir,de
      REAL (KIND=dp) :: ap,am,rp,rm,bp,bm,sq
      sq=SQRT(3._dp)      
      xir=fxih(t,p,h)/r
      de=fdelta(t,p)
      e=0.0_dp
      con=4*(4+de)*xir**2/13
      DO i=0,nmax-1
         da=alpha(i+1)-alpha(i)
         db=beta(i+1)-beta(i)
         e=e+con*(i+0.5)*db**2
         rp=(i+(3+sq)/6)
         rm=(i+(3-sq)/6)
         bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
         bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6
         e=e+0.5*con*da**2*(rp*SIN(bp)**2+rm*SIN(bm)**2)
         e=e+0.5*con*(SIN(bp)**2/rp+SIN(bm)**2/rm)
      END DO
      con=-(2+de)*xir**2/26
      DO i=0,nmax-1
         da=alpha(i+1)-alpha(i)
         db=beta(i+1)-beta(i)
         rp=(i+(3+sq)/6)
         rm=(i+(3-sq)/6)
         ap=(3+sq)*alpha(i+1)/6+(3-sq)*alpha(i)/6
         am=(3-sq)*alpha(i+1)/6+(3+sq)*alpha(i)/6
         bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
         bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6
         help=(SQRT(5.)*SIN(ap)-sq*COS(bp)*COS(ap))*db
         help=help+SIN(bp)*(SQRT(5.)*COS(bp)*COS(ap)+sq*SIN(ap))*da
         help=help+SIN(bp)*(SQRT(5.)*COS(bp)*SIN(ap)-sq*COS(ap))/rp
         e=e+0.5*con*rp*help**2
         help=(SQRT(5.)*SIN(am)-sq*COS(bm)*COS(am))*db
         help=help+SIN(bm)*(SQRT(5.)*COS(bm)*COS(am)+sq*SIN(am))*da
         help=help+SIN(bm)*(SQRT(5.)*COS(bm)*SIN(am)-sq*COS(am))/rm
         e=e+0.5*con*rm*help**2
      END DO      
!
      e=e+4*(2+de)*xir**2*SIN(beta(nmax))**2/13
!
      e=e-2*lsg*xir**2*SIN(beta(nmax))**2/13
    END FUNCTION ebend


    SUBROUTINE egrad(alpha,beta,ga,gb)
      ! Calculates the first-order derivatives
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp), DIMENSION(0:nmax) :: alpha,beta,ga,gb     
      REAL (KIND=dp) :: nr,nf,nz,help,bn,an,con,bi,bip,bim,rp,rm
      REAL (KIND=dp) :: dap,dam,dbp,dbm,dar,xir,de,sq,bp,bm,chia
      REAL (KIND=dp) :: db,da,aim,ai,aip,ap,am
      REAL (KIND=dp) :: vd,rzr,rzf,rzz,s,c
      REAL (KIND=dp) :: vrp,vfp,vzp,vrm,vfm,vzm
      REAL (KIND=dp) :: lrp,lfp,lzp,lrm,lfm,lzm,wm,wp
      REAL (KIND=dp) :: blp,blm,apsip,apsim
      sq=SQRT(3._dp)
      dar=fdar(t,p,r)
      xir=fxih(t,p,h)/r
      de=fdelta(t,p)
      DO i=0,nmax
         ga(i)=0.0_dp
         gb(i)=0.0_dp
      END DO
      bn=beta(nmax)
      an=alpha(nmax)

     ! rp=(i+(3+sq)/6)*dx
     ! rm=(i+(3-sq)/6)*dx
     ! bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
     ! bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6
     ! e=e+0.5*dx*(rp*SIN(bp)**2+rm*SIN(bm)**2)

      DO i=1,nmax-1
         rp=(i+(3+sq)/6)*dx
         rm=(i+(3-sq)/6)*dx
         bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
         bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6
         gb(i)=gb(i)+0.5*dx*rp*SIN(2*bp)*(3-sq)/6
         gb(i)=gb(i)+0.5*dx*rm*SIN(2*bm)*(3+sq)/6
         rp=(i-1+(3+sq)/6)*dx
         rm=(i-1+(3-sq)/6)*dx
         bp=(3+sq)*beta(i)/6+(3-sq)*beta(i-1)/6
         bm=(3-sq)*beta(i)/6+(3+sq)*beta(i-1)/6
         gb(i)=gb(i)+0.5*dx*rp*SIN(2*bp)*(3+sq)/6
         gb(i)=gb(i)+0.5*dx*rm*SIN(2*bm)*(3-sq)/6
      END DO
      i=nmax
      rp=(i-1+(3+sq)/6)*dx
      rm=(i-1+(3-sq)/6)*dx
      bp=(3+sq)*beta(i)/6+(3-sq)*beta(i-1)/6
      bm=(3-sq)*beta(i)/6+(3+sq)*beta(i-1)/6
      gb(i)=gb(i)+0.5*dx*rp*SIN(2*bp)*(3+sq)/6
      gb(i)=gb(i)+0.5*dx*rm*SIN(2*bm)*(3-sq)/6      
!
      nr=-SIN(bn)*COS(an)
      nf=SIN(bn)*SIN(an)
      nz=COS(bn)
      help=SQRT(5.)*COS(2*bn)*COS(an)+SQRT(3.)*COS(bn)*SIN(an)
      gb(nmax)=gb(nmax)+5*dar*(SQRT(5.)*nz*nr-SQRT(3.)*nf)*help/8
      help=SQRT(5.)*nz*nf+SQRT(3.)*nr
      ga(nmax)=ga(nmax)-5*dar*(SQRT(5.)*nz*nr-SQRT(3.)*nf)*help/8
!
      gb(nmax)=gb(nmax)+4*(2+de)*xir**2*SIN(2*bn)/13
      con=4*(4+de)*xir**2/13
      DO i=1,nmax-1
         bim=beta(i-1)
         bi=beta(i)
         bip=beta(i+1)
         dap=alpha(i+1)-alpha(i)
         dam=alpha(i)-alpha(i-1)
         gb(i)=gb(i)-2*(i+0.5)*con*(bip-bi)
         gb(i)=gb(i)+2*(i-0.5)*con*(bi-bim)
!
         rp=(i+(3+sq)/6)
         rm=(i+(3-sq)/6)
         bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
         bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6
         gb(i)=gb(i)+0.5*con*dap**2*rp*SIN(2*bp)*(3-sq)/6
         gb(i)=gb(i)+0.5*con*dap**2*rm*SIN(2*bm)*(3+sq)/6
         ga(i)=ga(i)-con*dap*(rp*SIN(bp)**2+rm*SIN(bm)**2)
         rp=(i-1+(3+sq)/6)
         rm=(i-1+(3-sq)/6)
         bp=(3+sq)*beta(i)/6+(3-sq)*beta(i-1)/6
         bm=(3-sq)*beta(i)/6+(3+sq)*beta(i-1)/6
         gb(i)=gb(i)+0.5*con*dam**2*rp*SIN(2*bp)*(3+sq)/6
         gb(i)=gb(i)+0.5*con*dam**2*rm*SIN(2*bm)*(3-sq)/6
         ga(i)=ga(i)+con*dam*(rp*SIN(bp)**2+rm*SIN(bm)**2)
!
         rp=(i+(3+sq)/6)
         rm=(i+(3-sq)/6)
         bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
         bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6
         gb(i)=gb(i)+0.5*con*SIN(2*bp)*(3-sq)/(6*rp)
         gb(i)=gb(i)+0.5*con*SIN(2*bm)*(3+sq)/(6*rm)
         rp=(i-1+(3+sq)/6)
         rm=(i-1+(3-sq)/6)
         bp=(3+sq)*beta(i)/6+(3-sq)*beta(i-1)/6
         bm=(3-sq)*beta(i)/6+(3+sq)*beta(i-1)/6
         gb(i)=gb(i)+0.5*con*SIN(2*bp)*(3+sq)/(6*rp)
         gb(i)=gb(i)+0.5*con*SIN(2*bm)*(3-sq)/(6*rm)
      END DO
      i=0
      dap=alpha(i+1)-alpha(i)
      rp=(i+(3+sq)/6)
      rm=(i+(3-sq)/6)
      bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
      bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6      
      ga(i)=ga(i)-con*dap*(rp*SIN(bp)**2+rm*SIN(bm)**2)
      i=nmax
      bim=beta(i-1)
      bi=beta(i)
      gb(i)=gb(i)+2*(i-0.5)*con*(bi-bim)
      dam=alpha(i)-alpha(i-1)
      rp=(i-1+(3+sq)/6)
      rm=(i-1+(3-sq)/6)
      bp=(3+sq)*beta(i)/6+(3-sq)*beta(i-1)/6
      bm=(3-sq)*beta(i)/6+(3+sq)*beta(i-1)/6
      gb(i)=gb(i)+0.5*con*dam**2*rp*SIN(2*bp)*(3+sq)/6
      gb(i)=gb(i)+0.5*con*dam**2*rm*SIN(2*bm)*(3-sq)/6      
      ga(i)=ga(i)+con*dam*(rp*SIN(bp)**2+rm*SIN(bm)**2)
!
      gb(i)=gb(i)+0.5*con*SIN(2*bp)*(3+sq)/(6*rp)
      gb(i)=gb(i)+0.5*con*SIN(2*bm)*(3-sq)/(6*rm)      
!
      con=-(2+de)*xir**2/26
      DO i=0,nmax-1
         rp=(i+(3+sq)/6)
         rm=(i+(3-sq)/6)
         bi=beta(i)
         bip=beta(i+1)         
         ai=alpha(i)
         aip=alpha(i+1)         
         gb(i)=gb(i)+con*rp*((-bi + bip)*(-(sq* &
              COS(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.)) + &
              SQRT(5.)*SIN(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)) + &
              (-ai + aip)*(SQRT(5.)* &
              COS(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.) + &
              sq*SIN(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.))* &
              SIN(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.) + &
              ((-(sq*COS(((3 - sq)*ai)/6. & 
              + ((3 + sq)*aip)/6.)) + &
              SQRT(5.)*COS(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.)* &
              SIN(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.))* &
              SIN(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.))/rp)* &
              (sq*COS(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.) - &
              SQRT(5.)*SIN(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.) + &
              ((3 - sq)*(-ai + aip)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.)* &
              (SQRT(5.)*COS(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.) + &
              sq*SIN(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)))/ &
              6. + ((3 - sq)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.)* &
              (-(sq*COS(((3 - sq)*ai)/6. + &
              ((3 + sq)*aip)/6.)) + &
              SQRT(5.)*COS(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.)* &
              SIN(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)))/(6.*rp) &
              - ((-3 + sq)*(-bi + bip)* &
              COS(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)* &
              SIN(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.))/ &
              (2.*sq) + (SQRT(5.)*(-3 + sq)*(-ai + aip)* &
              COS(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)* &
              SIN(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.)**2)/6. + &
              (SQRT(5.)*(-3 + sq)* &
              SIN(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)* &
              SIN(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.)**2)/(6.*rp))
         ga(i)=ga(i)+con*rp*((-bi + bip)*(-(sq* &
                COS(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)* &
                COS(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.)) + &
             SQRT(5.)*SIN(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)) + &
          (-ai + aip)*(SQRT(5.)* &
              COS(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.) + &
             sq*SIN(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.))* &
           SIN(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.) + &
          ((-(sq*COS(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)) + &
               SQRT(5.)*COS(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.)* &
                SIN(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.))* &
             SIN(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.))/rp)* &
        ((-bi + bip)*((SQRT(5.)*(3 - sq)* &
                COS(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.))/6. - &
             ((-3 + sq)*COS(((3 - sq)*bi)/6. + &
                  ((3 + sq)*bip)/6.)* &
                SIN(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.))/ &
              (2.*sq)) - (SQRT(5.)* &
              COS(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.) + &
             sq*SIN(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.))* &
           SIN(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.) + &
          (((SQRT(5.)*(3 - sq)* &
                  COS(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.)* &
                  COS(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.))/6. - &
               ((-3 + sq)* &
                  SIN(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.))/ &
                (2.*sq))*SIN(((3 - sq)*bi)/6. + &
               ((3 + sq)*bip)/6.))/rp + &
          (-ai + aip)*(((3 - sq)* &
                COS(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.))/ &
              (2.*sq) + (SQRT(5.)*(-3 + sq)* &
                COS(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.)* &
                SIN(((3 - sq)*ai)/6. + ((3 + sq)*aip)/6.))/6.)* &
           SIN(((3 - sq)*bi)/6. + ((3 + sq)*bip)/6.))
         gb(i)=gb(i)+con*rm*((-bi + bip)*(-(sq* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.)) + &
              SQRT(5.)*SIN(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)) + &
              (-ai + aip)*(SQRT(5.)* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.) + &
              sq*SIN(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.))* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.) + &
              ((-(sq*COS(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)) + &
              SQRT(5.)*COS(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.)* &
              SIN(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.))* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.))/rm)* &
              (sq*COS(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.) - &
              SQRT(5.)*SIN(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.) + &
              ((3 + sq)*(-ai + aip)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.)* &
              (SQRT(5.)*COS(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.) + &
              sq*SIN(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)))/ &
              6. + ((3 + sq)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.)* &
              (-(sq*COS(((3 + sq)*ai)/6. + &
              ((3 - sq)*aip)/6.)) + &
              SQRT(5.)*COS(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.)* &
              SIN(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)))/(6.*rm) &
              - ((-3 - sq)*(-bi + bip)* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.))/ &
              (2.*sq) + (SQRT(5.)*(-3 - sq)*(-ai + aip)*  &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.)**2)/6. + &
              (SQRT(5.)*(-3 - sq)* &
              SIN(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.)**2)/(6.*rm)) 
         ga(i)=ga(i)+con*rm*((-bi + bip)*(-(sq* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.)) + &
              SQRT(5.)*SIN(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)) + &
              (-ai + aip)*(SQRT(5.)* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.) + &
              sq*SIN(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.))* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.) + &
              ((-(sq*COS(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)) + &
              SQRT(5.)*COS(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.)* &
              SIN(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.))* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.))/rm)* &
              ((-bi + bip)*((SQRT(5.)*(3 + sq)* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.))/6. - &
              ((-3 - sq)*COS(((3 + sq)*bi)/6. + &
              ((3 - sq)*bip)/6.)* &
              SIN(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.))/ &
              (2.*sq)) - (SQRT(5.)* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.) + &
              sq*SIN(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.))* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.) + &
              (((SQRT(5.)*(3 + sq)* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.))/6. - &
              ((-3 - sq)* &
              SIN(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.))/ &
              (2.*sq))*SIN(((3 + sq)*bi)/6. + &
              ((3 - sq)*bip)/6.))/rm + &
              (-ai + aip)*(((3 + sq)* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.))/ &
              (2.*sq) + (SQRT(5.)*(-3 - sq)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.)* &
              SIN(((3 + sq)*ai)/6. + ((3 - sq)*aip)/6.))/6.)* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bip)/6.))
      END DO
      DO i=1,nmax
         rp=(i-1+(3+sq)/6)
         rm=(i-1+(3-sq)/6)
         bi=beta(i)
         bim=beta(i-1)         
         ai=alpha(i)
         aim=alpha(i-1)
         gb(i)=gb(i)+con*rp*((bi - bim)*(-(sq* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.)) + &
              SQRT(5.)*SIN(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)) + &
              (ai - aim)*(SQRT(5.)*COS(((3 + sq)*ai)/6. + &
              ((3 - sq)*aim)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.) + &
              sq*SIN(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.))* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.) + &
              ((-(sq*COS(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)) + &
              SQRT(5.)*COS(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.)* &
              SIN(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.))* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.))/rp)* &
              (-(sq*COS(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.)) + &
              SQRT(5.)*SIN(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.) + &
              ((3 + sq)*(ai - aim)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.)* &
              (SQRT(5.)*COS(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.) + &
              sq*SIN(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)))/ &
              6. + ((3 + sq)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.)* &
              (-(sq*COS(((3 + sq)*ai)/6. + &
              ((3 - sq)*aim)/6.)) + &
              SQRT(5.)*COS(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.)* &
              SIN(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)))/(6.*rp) &
              - ((-3 - sq)*(bi - bim)* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.))/ &
              (2.*sq) + (SQRT(5.)*(-3 - sq)*(ai - aim)* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.)**2)/6. + &
              (SQRT(5.)*(-3 - sq)* &
              SIN(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.)**2)/(6.*rp))
         ga(i)=ga(i)+con*rp*((bi - bim)*(-(sq* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.)) + &
              SQRT(5.)*SIN(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)) + &
              (ai - aim)*(SQRT(5.)*COS(((3 + sq)*ai)/6. + &
              ((3 - sq)*aim)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.) + &
              sq*SIN(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.))* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.) + &
              ((-(sq*COS(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)) + &
              SQRT(5.)*COS(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.)* &
              SIN(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.))* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.))/rp)* &
              ((bi - bim)*((SQRT(5.)*(3 + sq)* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.))/6. - &
              ((-3 - sq)*COS(((3 + sq)*bi)/6. + &
              ((3 - sq)*bim)/6.)* &
              SIN(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.))/ &
              (2.*sq)) + (SQRT(5.)* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.) + &
              sq*SIN(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.))* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.) + &
              (((SQRT(5.)*(3 + sq)* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.))/6. - &
              ((-3 - sq)* &
              SIN(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.))/ &
              (2.*sq))*SIN(((3 + sq)*bi)/6. + &
              ((3 - sq)*bim)/6.))/rp + &
              (ai - aim)*(((3 + sq)* &
              COS(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.))/ &
              (2.*sq) + (SQRT(5.)*(-3 - sq)* &
              COS(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.)* &
              SIN(((3 + sq)*ai)/6. + ((3 - sq)*aim)/6.))/6.)* &
              SIN(((3 + sq)*bi)/6. + ((3 - sq)*bim)/6.))
         gb(i)=gb(i)+con*rm*((bi - bim)*(-(sq* &
              COS(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.)) + &
              SQRT(5.)*SIN(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)) + &
              (ai - aim)*(SQRT(5.)*COS(((3 - sq)*ai)/6. + &
              ((3 + sq)*aim)/6.)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.) + &
              sq*SIN(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.))* &
              SIN(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.) + &
              ((-(sq*COS(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)) + &
              SQRT(5.)*COS(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.)* &
              SIN(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.))* &
              SIN(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.))/rm)* &
              (-(sq*COS(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.)) + &
              SQRT(5.)*SIN(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.) + &
              ((3 - sq)*(ai - aim)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.)* &
              (SQRT(5.)*COS(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.) + &
              sq*SIN(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)))/ &
              6. + ((3 - sq)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.)* &
              (-(sq*COS(((3 - sq)*ai)/6. + &
              ((3 + sq)*aim)/6.)) + &
              SQRT(5.)*COS(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.)* &
              SIN(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)))/(6.*rm) &
              - ((-3 + sq)*(bi - bim)* &
              COS(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)* &
              SIN(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.))/ &
              (2.*sq) + (SQRT(5.)*(-3 + sq)*(ai - aim)* &
              COS(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)* &
              SIN(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.)**2)/6. + &
              (SQRT(5.)*(-3 + sq)* &
              SIN(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)* &
              SIN(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.)**2)/(6.*rm)) 
         ga(i)=ga(i)+con*rm*((bi - bim)*(-(sq* & 
              COS(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.)) + &
              SQRT(5.)*SIN(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)) + &
              (ai - aim)*(SQRT(5.)*COS(((3 - sq)*ai)/6. + &
              ((3 + sq)*aim)/6.)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.) + &
              sq*SIN(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.))* &
              SIN(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.) + &
              ((-(sq*COS(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)) + &
              SQRT(5.)*COS(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.)* &
              SIN(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.))* &
              SIN(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.))/rm)* &
              ((bi - bim)*((SQRT(5.)*(3 - sq)* &
              COS(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.))/6. - &
              ((-3 + sq)*COS(((3 - sq)*bi)/6. + &
              ((3 + sq)*bim)/6.)* &
              SIN(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.))/ &
              (2.*sq)) + (SQRT(5.)* &
              COS(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.) + &
              sq*SIN(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.))* &
              SIN(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.) + &
              (((SQRT(5.)*(3 - sq)* &
              COS(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.))/6. - &
              ((-3 + sq)* &
              SIN(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.))/ &
              (2.*sq))*SIN(((3 - sq)*bi)/6. + &
              ((3 + sq)*bim)/6.))/rm + &
              (ai - aim)*(((3 - sq)* &
              COS(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.))/ &
              (2.*sq) + (SQRT(5.)*(-3 + sq)* &
              COS(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.)* &
              SIN(((3 - sq)*ai)/6. + ((3 + sq)*aim)/6.))/6.)* &
              SIN(((3 - sq)*bi)/6. + ((3 + sq)*bim)/6.))     
      END DO
!
      gb(nmax)=gb(nmax)-2*lsg*xir**2*SIN(2*bn)/13
!
      c=-0.25_dp
      s=SQRT(15.)/4.0_dp
      vd=fvd(t,p)
      DO i=0,nmax-1
         vrp=evrp(i)
         vfp=evfp(i)
         vzp=evzp(i)
         rp=(i+(3+sq)/6)*dx
         ap=(3+sq)*alpha(i+1)/6+(3-sq)*alpha(i)/6
         bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
         nr=-SIN(bp)*COS(ap)
         nf=SIN(bp)*SIN(ap)
         nz=COS(bp)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=vrp*(-(1-c)*COS(2*bp)*COS(ap)-s*COS(bp)*SIN(ap))
         help=help+vfp*((1-c)*COS(2*bp)*SIN(ap)-s*COS(bp)*COS(ap))
         help=help+vzp*(-(1-c)*SIN(2*bp))
         gb(i)=gb(i)-((3-sq)/3)*dx*rp*(rzr*vrp+rzf*vfp+rzz*vzp)*help/(5*vd**2)
         help=vrp*((1-c)*SIN(bp)*COS(bp)*SIN(ap)-s*SIN(bp)*COS(ap))
         help=help+vfp*((1-c)*SIN(bp)*COS(bp)*COS(ap)+s*SIN(bp)*SIN(ap))
         ga(i)=ga(i)-((3-sq)/3)*dx*rp*(rzr*vrp+rzf*vfp+rzz*vzp)*help/(5*vd**2)
         vrm=evrm(i)
         vfm=evfm(i)
         vzm=evzm(i)
         rm=(i+(3-sq)/6)*dx
         am=(3-sq)*alpha(i+1)/6+(3+sq)*alpha(i)/6
         bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6
         nr=-SIN(bm)*COS(am)
         nf=SIN(bm)*SIN(am)
         nz=COS(bm)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=vrm*(-(1-c)*COS(2*bm)*COS(am)-s*COS(bm)*SIN(am))
         help=help+vfm*((1-c)*COS(2*bm)*SIN(am)-s*COS(bm)*COS(am))
         help=help+vzm*(-(1-c)*SIN(2*bm))
         gb(i)=gb(i)-((3+sq)/3)*dx*rm*(rzr*vrm+rzf*vfm+rzz*vzm)*help/(5*vd**2)
         help=vrm*((1-c)*SIN(bm)*COS(bm)*SIN(am)-s*SIN(bm)*COS(am))
         help=help+vfm*((1-c)*SIN(bm)*COS(bm)*COS(am)+s*SIN(bm)*SIN(am))
         ga(i)=ga(i)-((3+sq)/3)*dx*rm*(rzr*vrm+rzf*vfm+rzz*vzm)*help/(5*vd**2)
      END DO
      DO i=1,nmax
         vrp=evrp(i-1)
         vfp=evfp(i-1)
         vzp=evzp(i-1)
         rp=(i-1+(3+sq)/6)*dx
         ap=(3+sq)*alpha(i)/6+(3-sq)*alpha(i-1)/6
         bp=(3+sq)*beta(i)/6+(3-sq)*beta(i-1)/6
         nr=-SIN(bp)*COS(ap)
         nf=SIN(bp)*SIN(ap)
         nz=COS(bp)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=vrp*(-(1-c)*COS(2*bp)*COS(ap)-s*COS(bp)*SIN(ap))
         help=help+vfp*((1-c)*COS(2*bp)*SIN(ap)-s*COS(bp)*COS(ap))
         help=help+vzp*(-(1-c)*SIN(2*bp))
         gb(i)=gb(i)-((3+sq)/3)*dx*rp*(rzr*vrp+rzf*vfp+rzz*vzp)*help/(5*vd**2)
         help=vrp*((1-c)*SIN(bp)*COS(bp)*SIN(ap)-s*SIN(bp)*COS(ap))
         help=help+vfp*((1-c)*SIN(bp)*COS(bp)*COS(ap)+s*SIN(bp)*SIN(ap))
         ga(i)=ga(i)-((3+sq)/3)*dx*rp*(rzr*vrp+rzf*vfp+rzz*vzp)*help/(5*vd**2)
         vrm=evrm(i-1)
         vfm=evfm(i-1)
         vzm=evzm(i-1)
         rm=(i-1+(3-sq)/6)*dx
         am=(3-sq)*alpha(i)/6+(3+sq)*alpha(i-1)/6
         bm=(3-sq)*beta(i)/6+(3+sq)*beta(i-1)/6
         nr=-SIN(bm)*COS(am)
         nf=SIN(bm)*SIN(am)
         nz=COS(bm)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=vrm*(-(1-c)*COS(2*bm)*COS(am)-s*COS(bm)*SIN(am))
         help=help+vfm*((1-c)*COS(2*bm)*SIN(am)-s*COS(bm)*COS(am))
         help=help+vzm*(-(1-c)*SIN(2*bm))
         gb(i)=gb(i)-((3-sq)/3)*dx*rm*(rzr*vrm+rzf*vfm+rzz*vzm)*help/(5*vd**2)
         help=vrm*((1-c)*SIN(bm)*COS(bm)*SIN(am)-s*SIN(bm)*COS(am))
         help=help+vfm*((1-c)*SIN(bm)*COS(bm)*COS(am)+s*SIN(bm)*SIN(am))
         ga(i)=ga(i)-((3-sq)/3)*dx*rm*(rzr*vrm+rzf*vfm+rzz*vzm)*help/(5*vd**2)
      END DO
!
!
      DO i=0,nmax-1
         lrp=elrp(i)
         lfp=elfp(i)
         lzp=elzp(i)
         wp=ewp(i)
         rp=(i+(3+sq)/6)*dx
         ap=(3+sq)*alpha(i+1)/6+(3-sq)*alpha(i)/6
         bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
         nr=-SIN(bp)*COS(ap)
         nf=SIN(bp)*SIN(ap)
         nz=COS(bp)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=lrp*(-(1-c)*COS(2*bp)*COS(ap)-s*COS(bp)*SIN(ap))
         help=help+lfp*((1-c)*COS(2*bp)*SIN(ap)-s*COS(bp)*COS(ap))
         help=help+lzp*(-(1-c)*SIN(2*bp))
         gb(i)=gb(i)+lo*wp*((3-sq)/3)*dx*rp*(rzr*lrp+rzf*lfp+rzz*lzp)*help/10
         help=lrp*((1-c)*SIN(bp)*COS(bp)*SIN(ap)-s*SIN(bp)*COS(ap))
         help=help+lfp*((1-c)*SIN(bp)*COS(bp)*COS(ap)+s*SIN(bp)*SIN(ap))
         ga(i)=ga(i)+lo*wp*((3-sq)/3)*dx*rp*(rzr*lrp+rzf*lfp+rzz*lzp)*help/10
         lrm=elrm(i)
         lfm=elfm(i)
         lzm=elzm(i)
         wm=ewm(i)
         rm=(i+(3-sq)/6)*dx
         am=(3-sq)*alpha(i+1)/6+(3+sq)*alpha(i)/6
         bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6
         nr=-SIN(bm)*COS(am)
         nf=SIN(bm)*SIN(am)
         nz=COS(bm)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=lrm*(-(1-c)*COS(2*bm)*COS(am)-s*COS(bm)*SIN(am))
         help=help+lfm*((1-c)*COS(2*bm)*SIN(am)-s*COS(bm)*COS(am))
         help=help+lzm*(-(1-c)*SIN(2*bm))
         gb(i)=gb(i)+lo*wm*((3+sq)/3)*dx*rm*(rzr*lrm+rzf*lfm+rzz*lzm)*help/10
         help=lrm*((1-c)*SIN(bm)*COS(bm)*SIN(am)-s*SIN(bm)*COS(am))
         help=help+lfm*((1-c)*SIN(bm)*COS(bm)*COS(am)+s*SIN(bm)*SIN(am))
         ga(i)=ga(i)+lo*wm*((3+sq)/3)*dx*rm*(rzr*lrm+rzf*lfm+rzz*lzm)*help/10
      END DO
      DO i=1,nmax
         lrp=elrp(i-1)
         lfp=elfp(i-1)
         lzp=elzp(i-1)
         wp=ewp(i-1)
         rp=(i-1+(3+sq)/6)*dx
         ap=(3+sq)*alpha(i)/6+(3-sq)*alpha(i-1)/6
         bp=(3+sq)*beta(i)/6+(3-sq)*beta(i-1)/6
         nr=-SIN(bp)*COS(ap)
         nf=SIN(bp)*SIN(ap)
         nz=COS(bp)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=lrp*(-(1-c)*COS(2*bp)*COS(ap)-s*COS(bp)*SIN(ap))
         help=help+lfp*((1-c)*COS(2*bp)*SIN(ap)-s*COS(bp)*COS(ap))
         help=help+lzp*(-(1-c)*SIN(2*bp))
         gb(i)=gb(i)+lo*wp*((3+sq)/3)*dx*rp*(rzr*lrp+rzf*lfp+rzz*lzp)*help/10
         help=lrp*((1-c)*SIN(bp)*COS(bp)*SIN(ap)-s*SIN(bp)*COS(ap))
         help=help+lfp*((1-c)*SIN(bp)*COS(bp)*COS(ap)+s*SIN(bp)*SIN(ap))
         ga(i)=ga(i)+lo*wp*((3+sq)/3)*dx*rp*(rzr*lrp+rzf*lfp+rzz*lzp)*help/10
         lrm=elrm(i-1)
         lfm=elfm(i-1)
         lzm=elzm(i-1)
         wm=ewm(i-1)
         rm=(i-1+(3-sq)/6)*dx
         am=(3-sq)*alpha(i)/6+(3+sq)*alpha(i-1)/6
         bm=(3-sq)*beta(i)/6+(3+sq)*beta(i-1)/6
         nr=-SIN(bm)*COS(am)
         nf=SIN(bm)*SIN(am)
         nz=COS(bm)
         rzr=(1-c)*nz*nr-s*nf
         rzf=(1-c)*nz*nf+s*nr
         rzz=c+(1-c)*nz**2
         help=lrm*(-(1-c)*COS(2*bm)*COS(am)-s*COS(bm)*SIN(am))
         help=help+lfm*((1-c)*COS(2*bm)*SIN(am)-s*COS(bm)*COS(am))
         help=help+lzm*(-(1-c)*SIN(2*bm))
         gb(i)=gb(i)+lo*wm*((3-sq)/3)*dx*rm*(rzr*lrm+rzf*lfm+rzz*lzm)*help/10
         help=lrm*((1-c)*SIN(bm)*COS(bm)*SIN(am)-s*SIN(bm)*COS(am))
         help=help+lfm*((1-c)*SIN(bm)*COS(bm)*COS(am)+s*SIN(bm)*SIN(am))
         ga(i)=ga(i)+lo*wm*((3-sq)/3)*dx*rm*(rzr*lrm+rzf*lfm+rzz*lzm)*help/10
      END DO

      chia=fchia(t,p)
      DO i=1,nmax-1              
         rp=(i+(3+sq)/6)*dx
         rm=(i+(3-sq)/6)*dx
         bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
         bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6
         !blp=ACOS(-1/4+5/4*COS(bp)**2)
         !blm=ACOS(-1/4+5/4*COS(bm)**2)
         apsip=(3+sq)*apsi(i+1)/6+(3-sq)*apsi(i)/6
         apsim=(3-sq)*apsi(i+1)/6+(3+sq)*apsi(i)/6
         gb(i)=gb(i)+chia*dx*((nub/nu0)**2)*rp*(SIN(bp*2)*0.5*(apsip**2)*(3-sq)/6)
         gb(i)=gb(i)+chia*dx*((nub/nu0)**2)*rm*(SIN(bm*2)*0.5*(apsim**2)*(3+sq)/6)

         rp=(i-1+(3+sq)/6)*dx
         rm=(i-1+(3-sq)/6)*dx
         bp=(3+sq)*beta(i)/6+(3-sq)*beta(i-1)/6
         bm=(3-sq)*beta(i)/6+(3+sq)*beta(i-1)/6
         !blp=ACOS(-1/4+5/4*COS(bp)**2)
         !blm=ACOS(-1/4+5/4*COS(bm)**2)
         apsip=(3+sq)*apsi(i)/6+(3-sq)*apsi(i-1)/6
         apsim=(3-sq)*apsi(i)/6+(3+sq)*apsi(i-1)/6
         gb(i)=gb(i)+chia*dx*((nub/nu0)**2)*rp*(SIN(bp*2)*0.5*(apsip**2)*(3+sq)/6)
         gb(i)=gb(i)+chia*dx*((nub/nu0)**2)*rm*(SIN(bm*2)*0.5*(apsim**2)*(3-sq)/6)   
      END DO

      gb(0)=0._dp
    END SUBROUTINE egrad

    ! rp=(i+(3+sq)/6)*dx
    ! rm=(i+(3-sq)/6)*dx
    ! bp=(3+sq)*beta(i+1)/6+(3-sq)*beta(i)/6
    ! bm=(3-sq)*beta(i+1)/6+(3+sq)*beta(i)/6
    ! gb(i)=gb(i)+0.5*dx*rp*SIN(2*bp)*(3-sq)/6
    ! gb(i)=gb(i)+0.5*dx*rm*SIN(2*bm)*(3+sq)/6
    ! rp=(i-1+(3+sq)/6)*dx
    ! rm=(i-1+(3-sq)/6)*dx
    ! bp=(3+sq)*beta(i)/6+(3-sq)*beta(i-1)/6
    ! bm=(3-sq)*beta(i)/6+(3+sq)*beta(i-1)/6
    ! gb(i)=gb(i)+0.5*dx*rp*SIN(2*bp)*(3+sq)/6
    ! gb(i)=gb(i)+0.5*dx*rm*SIN(2*bm)*(3-sq)/6

END MODULE energies

