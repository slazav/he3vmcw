MODULE subs

  USE general

  IMPLICIT NONE

  INTEGER :: maxi=100

  REAL (KIND=dp), DIMENSION(16) :: plist = & 
       (/ 0._dp,1._dp,2._dp,3._dp,6._dp,9._dp,12._dp,15._dp, &
       18._dp,21._dp,24._dp,27._dp,30._dp,33._dp,34.39_dp,1.E5_dp /) 

  REAL (KIND=dp) :: mbare=5.01 !    10^(-24) g   
  REAL (KIND=dp) :: kb=1.38 !   10^(-16) erg/K  
  REAL (KIND=dp) :: navo=6.02 !         10^23   
  REAL (KIND=dp) :: hbar=1.05 !    10^(-27) erg s
  REAL (KIND=dp) :: gyro=-2.04 !  10^4 1/(Gauss s)
 

  CONTAINS



    FUNCTION bcsgap(t) RESULT(gap)
      ! BCS gap / (kB Tc) for pure 3He-B, t = T / Tc
      ! Newton iteration based on a note by EVT & RH
      IMPLICIT NONE
      INTEGER :: n
      INTEGER, PARAMETER :: m=30
      REAL (KIND=dp) :: t,root,y,dy,ynew,g,dg,gap
      dy=1.0
      ynew=1.7638*SQRT(1-t)/(2*pi)
      DO WHILE (ABS(dy) > 1.0E-8)
         y=ynew
         root=SQRT((m*t)**2+y**2)
         g=LOG((m*t+root)/(2*m))-(1.0_dp/m**2-m*(t/root)**3)/24
         dg=y/(root*(m*t+root))-m*t**3*y/(8*root**5)
         DO n=1,m
            root=SQRT((t*(n-0.5))**2+y**2)
            g=g+1.0_dp/(n-0.5)-t/root
            dg=dg+t*y/root**3
         END DO
         dy=g/dg
         ynew=ynew-dy
      END DO
      gap=2*pi*ynew
      IF (t >= 1.0) gap=0._dp
    END FUNCTION bcsgap

 

    FUNCTION trivial(t,p) RESULT(gap)
      ! Trivial strong-coupling correction to the
      ! BCS energy gap
      IMPLICIT NONE
      INTEGER :: it,ic
      REAL (KIND=dp) :: t,p,gap,dcpcn,wt1,wt2,wc1,wc2,corr
      REAL (KIND=dp), DIMENSION(5) :: c = (/ &
           1.43_dp,1.6_dp,1.8_dp,2._dp,2.2_dp /)
      REAL (KIND=dp), DIMENSION(0:10,5) :: x
      dcpcn=41.9/volf(p)+0.322
      x(10,1:5)=(/ 1._dp,1.056_dp,1.115_dp,1.171_dp,1.221_dp /)
      x(9,1:5)=(/ 1._dp,1.048_dp,1.097_dp,1.141_dp,1.18_dp /)
      x(8,1:5)=(/ 1._dp,1.041_dp,1.083_dp,1.119_dp,1.15_dp /)
      x(7,1:5)=(/ 1._dp,1.036_dp,1.072_dp,1.102_dp,1.128_dp /)
      x(6,1:5)=(/ 1._dp,1.032_dp,1.063_dp,1.089_dp,1.112_dp /)
      x(5,1:5)=(/ 1._dp,1.028_dp,1.056_dp,1.079_dp,1.099_dp /)
      x(4,1:5)=(/ 1._dp,1.026_dp,1.051_dp,1.073_dp,1.091_dp /)
      x(3,1:5)=(/ 1._dp,1.024_dp,1.049_dp,1.069_dp,1.086_dp /)
      x(2,1:5)=(/ 1._dp,1.024_dp,1.048_dp,1.068_dp,1.085_dp /)
      x(1,1:5)=(/ 1._dp,1.024_dp,1.048_dp,1.068_dp,1.085_dp /)
      x(0,1:5)=(/ 1._dp,1.024_dp,1.048_dp,1.068_dp,1.085_dp /)
      it=INT(t*10-1E-5)
      wt1=(0.1*(it+1)-t)/0.1
      wt2=1-wt1
      ic=1
      DO
         IF (dcpcn < c(ic+1)) EXIT
         ic=ic+1
         IF (ic == 4) EXIT
      END DO
      wc1=(c(ic+1)-dcpcn)/(c(ic+1)-c(ic))
      wc2=1-wc1
      corr=wt1*(wc1*x(it,ic)+wc2*x(it,ic+1))
      corr=corr+wt2*(wc1*x(it+1,ic)+wc2*x(it+1,ic+1))
      gap=bcsgap(t)*corr
    END FUNCTION trivial



    FUNCTION f0af(p) RESULT(f0a)
      ! F0a as a function of pressure p
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp) :: p,weight,f0a
      REAL (KIND=dp), DIMENSION(16) :: flist = & 
           (/ -0.698_dp,-0.71_dp,-0.718_dp,-0.724_dp,-0.734_dp, &
           -0.741_dp,-0.746_dp,-0.751_dp,-0.754_dp,-0.756_dp, &
           -0.757_dp,-0.757_dp,-0.754_dp,-0.755_dp,-0.753_dp, &
           -0.753_dp /)
      i=16
      DO WHILE (plist(i) > p) 
         i=i-1
      END DO
      weight=(plist(i+1)-p)/(plist(i+1)-plist(i))
      f0a=weight*flist(i)+(1-weight)*flist(i+1)
    END FUNCTION f0af
    


    FUNCTION f1sf(p) RESULT(f1s)
      ! F1s as a function of pressure p
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp) :: p,weight,f1s
      REAL (KIND=dp), DIMENSION(16) :: flist = & 
           (/ 5.39_dp,5.78_dp,6.14_dp,6.49_dp,7.45_dp, &
           8.31_dp,9.09_dp,9.85_dp,10.6_dp,11.34_dp, &
           12.07_dp,12.79_dp,13.5_dp,14.21_dp,14.56_dp, &
           14.56_dp /) 
      i=16
      DO WHILE (plist(i) > p) 
         i=i-1
      END DO
      weight=(plist(i+1)-p)/(plist(i+1)-plist(i))
      f1s=weight*flist(i)+(1-weight)*flist(i+1)
    END FUNCTION f1sf



    FUNCTION f1af(p) RESULT(f1a)
      ! F1a as a function of pressure p
      IMPLICIT NONE
      REAL (KIND=dp) :: p,f1a
      f1a=-0.5678+p*(-0.04753+p*(1.791E-3-p*2.273E-5))
    END FUNCTION f1af



    FUNCTION z3f(t,gap) RESULT(z3)
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp) :: t,gap,z3,y,help,corr1,corr2,sum
      y=gap/(2*pi)
      sum=0._dp
      DO i=1,maxi
         sum=sum+t/((t*(i-0.5))**2+y**2)**1.5
      END DO
      help=SQRT((maxi*t)**2+y**2)
      corr1=1/(help*(maxi*t+help))
      corr2=maxi*t**3/(8*help**5)
      z3=y**2*(sum+corr1-corr2)
    END FUNCTION z3f


    
    FUNCTION z5f(t,gap) RESULT(z5)
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp) :: t,gap,z5,y,help,corr1,corr2,sum
      y=gap/(2*pi)
      sum=0._dp
      DO i=1,maxi
         sum=sum+t/((t*(i-0.5))**2+y**2)**2.5
      END DO
      help=SQRT((maxi*t)**2+y**2)
      corr1=(maxi*t+2*help)/(3*help**3*(maxi*t+help)**2)
      corr2=5*maxi*t**3/(24*help**7)
      z5=y**4*(sum+corr1-corr2)
    END FUNCTION z5f    



    FUNCTION z7f(t,gap) RESULT(z7)
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp) :: t,gap,z7,y,help,mt,corr1,corr2,sum
      y=gap/(2*pi)
      sum=0._dp
      mt=maxi*t
      DO i=1,maxi
         sum=sum+t/((t*(i-0.5))**2+y**2)**3.5
      END DO
      help=SQRT(mt**2+y**2)
      corr1=(11*mt*mt+9*mt*help+8*y*y)/(15*help**5*(mt+help)**3)
      corr2=7*mt**3/(24*help**9)
      z7=y**6*(sum+corr1-corr2)
    END FUNCTION z7f



    FUNCTION yoshida(t,p) RESULT(y)
      ! Trivial-corrected Yoshida function
      IMPLICIT NONE
      REAL (KIND=dp) :: t,p,y
      y=1-z3f(t,trivial(t,p))
    END FUNCTION yoshida



    FUNCTION gdf(p) RESULT(gd)
      ! Dipole coefficient in units of 1e32 1/(erg cm^3)
      IMPLICIT NONE
      REAL (KIND=dp) :: p,gd
      gd=0.27733+p*(5.8087E-4+2.515E-4*p)
    END FUNCTION gdf



    FUNCTION tcf(p) RESULT(tc)
      ! Tc in mK as a function of pressure p
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp) :: p,weight,tc
      REAL (KIND=dp), DIMENSION(16) :: flist = & 
           (/ 0.929_dp,1.061_dp,1.181_dp,1.29_dp,1.56_dp, & 
           1.769_dp,1.934_dp,2.067_dp,2.177_dp,2.267_dp, &
           2.339_dp,2.395_dp,2.438_dp,2.474_dp,2.491_dp, &
           2.491_dp /) 
      i=16
      DO WHILE (plist(i) > p) 
         i=i-1
      END DO
      weight=(plist(i+1)-p)/(plist(i+1)-plist(i))
      tc=weight*flist(i)+(1-weight)*flist(i+1)
    END FUNCTION tcf



    FUNCTION volf(p) RESULT(vol)
      ! Specific volume in cm^3 as a function of pressure p
      IMPLICIT NONE
      INTEGER :: i
      REAL (KIND=dp) :: p,weight,vol
      REAL (KIND=dp), DIMENSION(16) :: flist = & 
           (/ 36.84_dp,35.74_dp,34.78_dp,33.95_dp,32.03_dp, &
           30.71_dp,29.71_dp,28.89_dp,28.18_dp,27.55_dp, &
           27.01_dp,26.56_dp,26.17_dp,25.75_dp,25.5_dp, &
           25.5_dp /) 
      i=16
      DO WHILE (plist(i) > p) 
         i=i-1
      END DO
      weight=(plist(i+1)-p)/(plist(i+1)-plist(i))
      vol=weight*flist(i)+(1-weight)*flist(i+1)
    END FUNCTION volf



    FUNCTION mefff(p) RESULT(meff)
      ! Effective mass in units of 10^(-24) g
      IMPLICIT NONE
      REAL (KIND=dp) :: p,meff
      meff=mbare*(1+f1sf(p)/3)
    END FUNCTION mefff



    FUNCTION rhof(p) RESULT(rho)
      ! Total density in units of g/cm^3
      IMPLICIT NONE
      REAL (KIND=dp) :: p,rho
      rho=0.1*mbare*navo/volf(p)
    END FUNCTION rhof



    FUNCTION n0f(p) RESULT(n0)
      ! One-spin DOS at the Fermi energy in 10^38 1/(erg cm^3)
      IMPLICIT NONE
      REAL (KIND=dp) :: p,n0
      n0=mefff(p)*(3*rhof(p)/(pi*mbare))**(1/3._dp)/(2*pi*hbar**2)
    END FUNCTION n0f



    FUNCTION vff(p) RESULT(vf)
      ! Fermi velocity in units of 10^3 cm/s
      IMPLICIT NONE
      REAL (KIND=dp) :: p,vf
      vf=100*hbar*(3*pi*pi*rhof(p)/mbare)**(1/3._dp)/mefff(p)
    END FUNCTION vff



    FUNCTION xiglf(t,p) RESULT(xigl)
      ! Extrapolated GL coherence length in 10^(-5) cm
      IMPLICIT NONE
      REAL (KIND=dp) :: t,p,xigl
      xigl=hbar*vff(p)/(SQRT(10.)*kb*tcf(p)*bcsgap(t))
    END FUNCTION xiglf


    FUNCTION vneq(t,p,omega,r) RESULT(nv)
      ! Equilibrium vortex number
      IMPLICIT NONE
      REAL (KIND=dp) :: t,p,omega,r,nv
      REAL (KIND=dp) :: kappa,help,rv,rc
      kappa=6.65E-4_dp
      rc=xiglf(t,p)*1.0E-5
      rv=SQRT(kappa/(2*pi*omega))
      help=1.0_dp-SQRT(kappa*LOG(rv/rc)/(4*pi*omega*r*r))
!      help=1.0_dp
      nv=2*pi*omega*r**2/kappa*help**2
    END FUNCTION vneq


END MODULE subs


