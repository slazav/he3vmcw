MODULE text

  USE general

  USE subs

  IMPLICIT NONE

  REAL (KIND=dp), SAVE :: flhvfix,chi


  CONTAINS


    FUNCTION fdelta(t,p) RESULT(de)
      ! Textural parameter delta
      IMPLICIT NONE
      REAL (KIND=dp) :: t,p,f1a,y,de
      y=yoshida(t,p)
      f1a=f1af(p)
      de=f1a*(1-y)/(3+f1a*y)
    END FUNCTION fdelta


    FUNCTION fdar(t,p,r) RESULT(dar)
      ! Textural parameter d/aR (R in cm)
      IMPLICIT NONE
      REAL (KIND=dp) :: t,p,r,dar
      dar=10*fd(t,p)/(fa(t,p)*r)
    END FUNCTION fdar

    FUNCTION fchia(t,p) RESULT(chia)
      !Textural parameter chi/a, chi=susceptibility
      !chi to program units: multiply by 10^14
      IMPLICIT NONE
      REAL (KIND=dp) :: t,p,chia, dummy
      dummy=chi/(fa(t,p))
      chia=dummy*1E14
      END FUNCTION fchia

    FUNCTION fvd(t,p) RESULT(vd)
      ! Textural parameter vd in cm/s
      IMPLICIT NONE
      REAL (KIND=dp) :: t,p,vd
      vd=0.001*SQRT(0.4*fa(t,p)/flhv(t,p))
    END FUNCTION fvd


    FUNCTION fxih(t,p,h) RESULT(xih)
      ! Magnetic coherence length in cm (h in Gauss)
      IMPLICIT NONE
      REAL (KIND=dp) :: t,p,h,xih
      xih=SQRT(10000*65*flg2(t,p)/(8*fa(t,p)*h*h))
    END FUNCTION fxih


    FUNCTION fa(t,p) RESULT(a)
      ! Textural parameter a in 10^(-14) erg/(Gauss^2 cm^3)
      IMPLICIT NONE
      REAL (KIND=dp) :: t,p,a,y,z3,z5,f0a,gd,gap
      gd=gdf(p)
      gap=trivial(t,p)
      z3=z3f(t,gap)
      z5=z5f(t,gap)
      f0a=f0af(p)
      y=1-z3
      IF (ABS(1-t)<1.0e-8_dp) THEN
         a=2.5*gd*5*(0.5*hbar*gyro/(1+f0a*(2+y)/3))**2
      ELSE
         a=2.5*gd*(5-3*z5/z3)*(0.5*hbar*gyro/(1+f0a*(2+y)/3))**2
      END IF
    END FUNCTION fa


    FUNCTION fd(t,p) RESULT(d)
      ! Textural parameter d in 10^(-13) erg/(cm^2 Gauss^2)
      IMPLICIT NONE
      REAL (KIND=dp) :: t,p,d,xi,y,f0a,n0,d0
      y=yoshida(t,p)
      f0a=f0af(p)
      n0=n0f(p)
      xi=0.0_dp
      d0=2.2+p*0.5/34.39
      IF (ABS(1-t)>1.0e-8_dp) xi=xiglf(t,p) 
      d=n0*(hbar*gyro)**2*xi*d0*(1-y)/(4*(1+f0a)*(3+f0a*(2+y)))
    END FUNCTION fd    


    FUNCTION flg2(t,p) RESULT(lg2)
      ! Textural parameter lg2 in 10^(-10) erg/cm
      IMPLICIT NONE
      REAL (KIND=dp) :: t,p,lg2,factor,y,f1a
      y=yoshida(t,p)
      f1a=f1af(p)
      factor=(1+f1a/3)*(1-y)/(1+f1a*(2+3*y)/15)
      lg2=1000*hbar*hbar*navo*factor/(40*mefff(p)*volf(p))
    END FUNCTION flg2 

    FUNCTION flhv(t,p) RESULT(lhv)
      ! Textural parameter lhv in 10^(-8) g/(cm^3 Gauss^2)
      IMPLICIT NONE
      REAL (KIND=dp) :: t,p,lhv
      lhv=flhvfix
    END FUNCTION flhv
   


    FUNCTION flhvtheor(t,p) RESULT(lhv)
      ! Textural parameter lhv in 10^(-8) g/(cm^3 Gauss^2)
      IMPLICIT NONE
      REAL (KIND=dp) :: t,p,lhv,y,f1s,f0a,gap,z3,z5,z7
      REAL (KIND=dp) :: help1,help2,help3
      y=yoshida(t,p)
      f1s=f1sf(p)
      f0a=f0af(p)
      gap=trivial(t,p)
      z3=z3f(t,gap)
      z5=z5f(t,gap)
      z7=z7f(t,gap)
      help1=z3-0.9*z5+0.9*z5**2/z3-1.5*z7
      help2=(1+f1s/3)/(1+f1s*y/3)**2
      help3=(0.5*hbar*gyro/(1+f0a*(2+y)/3))**2
      lhv=rhof(p)*help1*help2*help3/(gap*kb*tcf(p))**2
    END FUNCTION flhvtheor


END MODULE text


