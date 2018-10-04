! Leggett equations in rotating frame. 1D calculations along rotation axis.
! Spin currents, spin diffusion, leggett-takagi, extra Mz relaxation.
! input:
!   T - time
!   X - x-coord
!   U(7)   - Mx My Mz Nx Ny Nz Theta
!   Ux(7)  - dU/dx
!   Uxx(7) - d2U/dx2
!   NPDE=7
! output
!   Fv(7)  - dU/dt
! Parameters are set using set_pars function.

C-- F ---------- EVALUATION OF F ------------------------------------
      subroutine F(T,X,U,UX,UXX,FV,NPDE)
        implicit none
        integer NPDE
        real*8 T,X,U,UX,UXX,FV
        dimension U(NPDE),UX(NPDE),UXX(NPDE),FV(NPDE)
!       RF-field in (rad/s), constant field (rad/s), pumping freq (rad/s),
        real*8 Wr, Wz, W0
!       Leggett freq (rad/s), Cpar (cm/s) and it's z derivative,
!       D (cm^2/s), Tf (s), T1 (s)
        real*8 WB, Cpar, dCpar, Diff, Tf, T1

        real*8 UMx,UMy,UMz, UNx,UNy,UNz, Uth, UN
        real*8 UMxm,UMym,UMzm
        real*8 GNx,GNy,GNz, GMx,GMy,GMz, Gth
        real*8 GGNx,GGNy,GGNz, GGMx,GGMy,GGMz, GGth
        real*8 DD45,ST,CT,CTM,CT1,CTG,UT,AUT,AF,DAF,FTN,DFTN,B
        real*8 UJX,UJY,UJZ,DJX,DJY,DJZ ! spin current J and dJ/dz

!       set parameters:
        call set_bulk_pars(T,X, Wr,Wz,W0,WB,
     *                     Cpar, dCpar, Diff, Tf, T1)

C       Normilize n vector length
        UN=dsqrt(U(4)**2+U(5)**2+U(6)**2)
        UNx = U(4)/UN
        UNy = U(5)/UN
        UNz = U(6)/UN

        UMx = U(1)
        UMy = U(2)
        UMz = U(3)
        Uth = U(7)

        UMxm = UMx-Wr/Wz
        UMym = UMy
        UMzm = UMz-1.0D0

!        GMx = UX(1)
!        GMy = UX(2)
!        GMz = UX(3)
        GNx = UX(4)
        GNy = UX(5)
        GNz = UX(6)
        Gth = UX(7)

        GGMx = UXX(1)
        GGMy = UXX(2)
        GGMz = UXX(3)
        GGNx = UXX(4)
        GGNy = UXX(5)
        GGNz = UXX(6)
        GGth = UXX(7)

        ST=dsin(Uth)
        CT=dcos(Uth)

        CTM=1.0D0-CT
        CT1=1.0D0+CT
        CTG=ST/CTM      ! ctg(T/2) = sin(T)/(1-cos(T))
        UT=ST*(1.0D0+4.0D0*CT)*0.2666666D0 ! 4/15 sin(t)*(1+4*cos(t))

        AUT=UT*WB**2/Wz
        AF=-Cpar**2/Wz
        DAF = -2D0*Cpar*dCpar/Wz

C       nx ny' - ny nx'
        DD45=UNx*GNy-GNx*UNy
C       (1-ct)*(nx ny' - ny nx') - st*nz' - nz t'
        FTN=CTM*DD45-ST*GNz-Gth*UNz
C       FTN'
        DFTN=CTM*(UNx*GGNy-GGNx*UNy)-ST*GGNz-GGth*UNz-
     *   CT1*Gth*GNz+ST*Gth*DD45

C       components of spin current, Jiz (without AF factor!)
        UJX=2.0D0*(Gth*UNx+ST*GNx+CTM*(UNy*GNz-GNy*UNz))
     *     + (CTM*UNx*UNz+UNy*ST)*FTN
        UJY=2.0D0*(Gth*UNy+ST*GNy-CTM*(UNx*GNz-GNx*UNz))
     *     + (CTM*UNy*UNz-UNx*ST)*FTN
        UJZ=2.0D0*(Gth*UNz+ST*GNz+CTM*(UNx*GNy-GNx*UNy))
     *     + (CTM*UNz**2+CT)*FTN

C       derivative of the spin flow: d/dz(Jiz - Diff*Mi')
C       Note: here non-uniform spin-wave velocity (DAF!=0) can be used,
C             but spin diffusion should be uniform
        DJX = AF*(2.0D0*(GGth*UNx+CT1*Gth*GNx+ST*GGNx+
     *   ST*Gth*(UNy*GNz-GNy*UNz)+CTM*(UNy*GGNz-GGNy*UNz))+
     *   (CTM*UNx*UNz+UNy*ST)*DFTN+(ST*Gth*UNx*UNz+
     *   CTM*(GNx*UNz+UNx*GNz)+GNy*ST+UNy*CT*Gth)*FTN)
     *   + DAF*UJX - DIFF*GGMx
        DJY=AF*(2.0D0*(GGth*UNy+CT1*Gth*GNy+ST*GGNy-
     *   ST*Gth*(UNx*GNz-GNx*UNz)-CTM*(UNx*GGNz-GGNx*UNz))+
     *   (CTM*UNy*UNz-UNx*ST)*DFTN+(ST*Gth*UNy*UNz+
     *   CTM*(GNy*UNz+UNy*GNz)-GNx*ST-UNx*CT*Gth)*FTN)
     *   + DAF*UJY - DIFF*GGMy
        DJZ=AF*(2.0D0*(GGth*UNz+CT1*Gth*GNz+ST*GGNz+
     *   ST*Gth*DD45+CTM*(UNx*GGNy-GGNx*UNy))+
     *   (CTM*UNz**2+CT)*DFTN+(ST*Gth*UNz**2+
     *   CTM*2.0D0*UNz*GNz-ST*Gth)*FTN)
     *   + DAF*UJZ - DIFF*GGMz

        B = UNx*UMxm + UMym*UNy + UMzm*UNz

C       Leggett equations
        FV(1)=   (Wz-W0)*UMy          + AUT*UNx - DJX
        FV(2)= - (Wz-W0)*UMx + Wr*UMz + AUT*UNy - DJY
        FV(3)=               - Wr*UMy + AUT*UNz - DJZ - UMzm/T1
        FV(4)=-W0*UNy-0.5D0*Wz*(UMzm*UNy-UMym*UNz+CTG*(B*UNx-UMxm))
        FV(5)= W0*UNx-0.5D0*Wz*(UMxm*UNz-UMzm*UNx+CTG*(B*UNy-UMym))
        FV(6)=       -0.5D0*Wz*(UMym*UNx-UMxm*UNy+CTG*(B*UNz-UMzm))
        FV(7)= Wz*B + UT/Tf
        return
      end

C-- BNDRY ------ BOUNDARY CONDITIONS -- B(U,UX)=Z(T) ------------

      subroutine BNDRY(T,X,U,UX,DBDU,DBDUX,DZDT,NPDE)
        implicit none
        integer NPDE
        dimension U(NPDE),UX(NPDE),DZDT(NPDE),
     *   DBDU(NPDE,NPDE),DBDUX(NPDE,NPDE)
        real*8 T,X,U,UX,DZDT,DBDU,DBDUX
        real*8 UN,UNx,UNy,UNz, Uth
        real*8 GNx,GNy,GNz, GMx,GMy,GMz, Gth
        real*8 GGNx,GGNy,GGNz, GGMx,GGMy,GGMz, GGth
        real*8 ST,ST2,CT,CTM,CTM2,DD45,FTN,CTF,STF
        real*8 FTN4,FTN5,FTN7,FTNX4,FTNX5,C46,C56,C66,C266
        real*8 W,AF,DA
        integer I,J

!       paramaters, same as in F (only Cpar, Diff used)
        real*8 Cpar, Diff, W0
        integer IBN

!       set parameters:
        call set_bndry_pars(T,X, W0, Cpar, Diff, IBN)


        do I=1,NPDE
          DZDT(I)=0.0D0
          do J=1,NPDE
            DBDU(I,J)=0.0D0
            DBDUX(I,J)=0.0D0
          enddo
        enddo

!       Symmetric th-soliton, N and theta are constants (and should be
!       started from NPD initial conditions), dM/dz = 0 (symmetric waves only)
        if (IBN.EQ.3) THEN
          do J=1,3
            DBDUX(J,J)=1.0D0
          enddo
          do J=4,NPDE
            DBDU(J,J)=1.0D0
          enddo
          return
        endif

!       Closed cell: no spin flow through walls
!       Jiz - Diff Mi' = 0
        if(IBN.EQ.2)THEN       ! CLOSED CELL
C         fix n vector length
          UN=dsqrt(U(4)**2+U(5)**2+U(6)**2)
          UNx=U(4)/UN
          UNy=U(5)/UN
          UNz=U(6)/UN
          Uth=U(7)

!          GMx = UX(1)
!          GMy = UX(2)
!          GMz = UX(3)
          GNx = UX(4)
          GNy = UX(5)
          GNz = UX(6)
          Gth = UX(7)

          ST=dsin(Uth)
          ST2=2.0D0*ST
          CT=dcos(Uth)
          CTM=1.0D0-CT
          CTM2=2.0D0*CTM

          DD45=UNx*GNy-GNx*UNy
          FTN=CTM*DD45-ST*GNz-Gth*UNz

          CTF=CTM*FTN
          STF=ST*FTN
          FTN4=CTM*GNy
          FTN5=-CTM*GNx
          FTN7=ST*DD45-CT*GNz
          FTNX4=-CTM*UNy
          FTNX5=CTM*UNx
          C46=CTM*UNx*UNz+UNy*ST
          C56=CTM*UNy*UNz-UNx*ST             !!!!!!!!!!!
          C66=CTM*UNz**2+CT
          C266=2.0D0-C66

          AF=-Cpar**2/W0
          DA=-Diff/AF

          DBDUX(4,1)=DA
          DBDUX(5,2)=DA
          DBDUX(6,3)=DA

          DBDU(4,4)=2.0D0*Gth+CTF*UNz+C46*FTN4
          DBDU(4,5)=CTM2*GNz+STF+C46*FTN5
          DBDU(4,6)=-CTM2*GNy+CTF*UNx-C46*Gth
          DBDU(4,7)=2.0D0*(CT*GNx+ST*(UNy*GNz-GNy*UNz))+
     *     STF*UNx*UNz+UNy*CT*FTN+C46*FTN7

          DBDU(5,4)=-CTM2*GNz-STF+C56*FTN4
          DBDU(5,5)=2.0D0*Gth+CTF*UNz+C56*FTN5
          DBDU(5,6)=CTM2*GNx+CTF*UNy-C56*Gth
          DBDU(5,7)=2.0D0*(CT*GNy-ST*(UNx*GNz-GNx*UNz))+
     *     STF*UNy*UNz-UNx*CT*FTN+C56*FTN7

          DBDU(6,4)=CTM2*GNy+C66*FTN4
          DBDU(6,5)=-CTM2*GNx+C66*FTN5
          DBDU(6,6)=2.0D0*UNz*CTF+C266*Gth
          DBDU(6,7)=2.0D0*(CT*GNz+ST*DD45)+STF*(UNz**2-1.0D0)+C66*FTN7

          DBDUX(4,4)=ST2+C46*FTNX4
          DBDUX(4,5)=-CTM2*UNz+C46*FTNX5
          DBDUX(4,6)=CTM2*UNy-C46*ST
          DBDUX(4,7)=2.0D0*UNx-C46*UNz

          DBDUX(5,4)=CTM2*UNz+C56*FTNX4
          DBDUX(5,5)=ST2+C56*FTNX5
          DBDUX(5,6)=-CTM2*UNx-C56*ST
          DBDUX(5,7)=2.0D0*UNy-C56*UNz

          DBDUX(6,4)=-CTM2*UNy+C66*FTNX4
          DBDUX(6,5)=CTM2*UNx+C66*FTNX5
          DBDUX(6,6)=C266*ST
          DBDUX(6,7)=C266*UNz
C          DBDU(7,4)=GNx         !!
C          DBDU(7,5)=GNy         !!
C          DBDU(7,6)=GNz         !!
C          DBDUX(7,4)=UNx         !!
C          DBDUX(7,5)=UNy         !!
C          DBDUX(7,6)=UNz         !!
          return
        endif
        return
      end

C-- DERIVF ----- SET UP DERIVATIVES ---------------------------------
CC ATTENTION :   DERIVF IS WRONG
      subroutine DERIVF(T,X,U,UX,UXX,DFDU,DFDUX,DFDUXX,NPDE)
        implicit none
        dimension U(NPDE),UX(NPDE),UXX(NPDE),
     *       DFDU(NPDE,NPDE),DFDUX(NPDE,NPDE),DFDUXX(NPDE,NPDE)
        integer I,J,NPDE
        real*8  T,X,U,UX,UXX,DFDU,DFDUX,DFDUXX

        do I=1,NPDE
          do J=1,NPDE
            DFDU(I,J)=0.0D0
            DFDUX(I,J)=0.0D0
            DFDUXX(I,J)=0.0D0
          enddo
        enddo
        return
      end
