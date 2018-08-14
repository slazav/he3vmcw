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
! use PARS common block (see vmcw_pars.fh) to set parameters:



C-- F ---------- EVALUATION OF F ------------------------------------
      subroutine F(T,X,U,UX,UXX,FV,NPDE)
        include 'vmcw_pars.fh'
        include 'he3_const.fh'
        integer NPDE
        real*8 T,X,U,UX,UXX,FV
        dimension U(NPDE),UX(NPDE),UXX(NPDE),FV(NPDE)
        real*8 WY, WZ, WL2, DW
        real*8 UN,UNx,UNy,UNz,UMxm,UMym,UMzm
        real*8 DD45,ST,CT,CTM,CT1,CTG,UT,AUT,AF,DAF,FTN,DFTN,B
        real*8 UJX,UJY,UJZ,DJX,DJY,DJZ ! spinn current J and dJ/dz

        real*8 DIFF,TF1,W0,WL
        real*8 AER_STEP

C       Arguments:
C         T - time
C         X - x-coord
C         U   - Mx My Mz Nx Ny Nz T 
C         UX  - dU/dx
C         UXX - d2U/dx2
C         FV  - result
C       Parameters (from vmcw_pars.fh and he3_const.fh)
C         GAM -- gyromagnetic ratio
C         PI
C         HR0, HR_SWR -- RF-field
C         H, GRAD     -- constant field, gradient
C         LP0, LP_SWR -- Larmor position
C         LF0, LF_SWR -- Leggett freq.
C         CPAR0, CPAR_SWR -- spin-wave vel
C         TF0, TF_SWR
C         DF0, DF_SWR
C         T11

C       calculate freq
        WY = GAM*(HR0+HR_SWR*T)
!     *   *(1D0-(2.0*X/CELL_LEN-1.0)*(2.0*X/CELL_LEN-1.0)*0.5)
!     *   *(1D0+(2.0*X/CELL_LEN-1.0)*0.6)

        WL = GAM*(H + GRAD*X)
        W0= GAM*(H + GRAD*(LP0+LP_SWR*T))
        DW = WL-W0

C       fix n vector length
        UN=dsqrt(U(4)**2+U(5)**2+U(6)**2)
        UNx = U(4)/UN
        UNy = U(5)/UN
        UNz = U(6)/UN

        UMxm = U(1)-WY/WL
        UMym = U(2)
        UMzm = U(3)-1.0D0

        WL2 = 0.5D0*WL
        DD45=UNx*UX(5)-UX(4)*UNy       ! Nx Ny` - Nx` Ny
        ST=dsin(U(7))
        CT=dcos(U(7))
        CTM=1.0D0-CT
        CT1=1.0D0+CT
        CTG=ST/CTM      ! ctg(T/2) = sin(T)/(1-cos(T))
        UT=ST*(1.0D0+4.0D0*CT)*0.2666666D0

        AUT=UT*(LF0+LF_SWR*T)**2/WL*4.0D0*PI**2
        AF=-(CPAR0+CPAR_SWR*T)**2/WL
        TF1=TF0+TF_SWR*T
        DIFF=DF0+DF_SWR*T

C       spatial modulation
        AF=AF*(1D0 - 0.5D0 * AER_STEP(X,0))
        DAF=-AF*0.5D0 * AER_STEP(X,1)
        AUT=AUT*(1D0 - 0.835D0 * AER_STEP(X,0))
        TF1=TF1*(1D0 - 0.5D0 * AER_STEP(X,0))

C       something..
        FTN=CTM*DD45-ST*UX(6)-UX(7)*UNz
        DFTN=CTM*(UNx*UXX(5)-UXX(4)*UNy)-ST*UXX(6)-UXX(7)*UNz-
     *   CT1*UX(7)*UX(6)+ST*UX(7)*DD45   !!! dFTN/dz

C       components of spin current, Ji
        UJX=2.0D0*(UX(7)*UNx+ST*UX(4)+CTM*(UNy*UX(6)-UX(5)*UNz))+
     *   (CTM*UNx*UNz+UNy*ST)*FTN
        UJY=2.0D0*(UX(7)*UNy+ST*UX(5)-CTM*(UNx*UX(6)-UX(4)*UNz))+
     *   (CTM*UNy*UNz-UNx*ST)*FTN
        UJZ=2.0D0*(UX(7)*UNz+ST*UX(6)+CTM*(UNx*UX(5)-UX(4)*UNy))+
     *   (CTM*UNz**2+CT)*FTN

C       dJi/dz
        DJX=AF*(2.0D0*(UXX(7)*UNx+CT1*UX(7)*UX(4)+ST*UXX(4)+
     *   ST*UX(7)*(UNy*UX(6)-UX(5)*UNz)+CTM*(UNy*UXX(6)-UXX(5)*UNz))+
     *   (CTM*UNx*UNz+UNy*ST)*DFTN+(ST*UX(7)*UNx*UNz+
     *   CTM*(UX(4)*UNz+UNx*UX(6))+UX(5)*ST+UNy*CT*UX(7))*FTN)-
     *   DIFF*UXX(1)+DAF*UJX
        DJY=AF*(2.0D0*(UXX(7)*UNy+CT1*UX(7)*UX(5)+ST*UXX(5)-
     *   ST*UX(7)*(UNx*UX(6)-UX(4)*UNz)-CTM*(UNx*UXX(6)-UXX(4)*UNz))+
     *   (CTM*UNy*UNz-UNx*ST)*DFTN+(ST*UX(7)*UNy*UNz+
     *   CTM*(UX(5)*UNz+UNy*UX(6))-UX(4)*ST-UNx*CT*UX(7))*FTN)-   !!!!!!!!!
     *   DIFF*UXX(2)+DAF*UJY
        DJZ=AF*(2.0D0*(UXX(7)*UNz+CT1*UX(7)*UX(6)+ST*UXX(6)+
     *   ST*UX(7)*DD45+CTM*(UNx*UXX(5)-UXX(4)*UNy))+
     *   (CTM*UNz**2+CT)*DFTN+(ST*UX(7)*UNz**2+
     *   CTM*2.0D0*UNz*UX(6)-ST*UX(7))*FTN)-DIFF*UXX(3)+DAF*UJZ

        B = UNx*UMxm + UMym*UNy + UMzm*UNz

C       Leggett equations
        FV(1)=   DW*U(2)           + AUT*UNx - DJX
        FV(2)= - DW*U(1) + WY*U(3) + AUT*UNy - DJY
        FV(3)=           - WY*U(2) + AUT*UNz - DJZ - UMzm*T11
        FV(4)= - W0*UNy - WL2*(UMzm*UNy-UMym*UNz+CTG*(B*UNx-UMxm))
        FV(5)=   W0*UNx - WL2*(UMxm*UNz-UMzm*UNx+CTG*(B*UNy-UMym))
        FV(6)=           - WL2*(UMym*UNx-UMxm*UNy+CTG*(B*UNz-UMzm))
        FV(7)= WL*B + UT/TF1
        return
      end

C-- BNDRY ------ BOUNDARY CONDITIONS -- B(U,UX)=Z(T) ------------
      subroutine BNDRY(T,X,U,UX,DBDU,DBDUX,DZDT,NPDE)
        include 'vmcw_pars.fh'
        include 'he3_const.fh'

        integer NPDE
        dimension U(NPDE),UX(NPDE),DZDT(NPDE),
     *   DBDU(NPDE,NPDE),DBDUX(NPDE,NPDE)
        real*8 T,X,U,UX,DZDT,DBDU,DBDUX
        real*8 UN,UNx,UNy,UNz
        real*8 ST,ST2,CT,CTM,CTM2,DD45,FTN,CTF,STF
        real*8 FTN4,FTN5,FTN7,FTNX4,FTNX5,C46,C56,C66,C266
        real*8 W,AF,DA
        real*8 AER_STEP
        integer I,J

        do I=1,NPDE
          DZDT(I)=0.0D0
          do J=1,NPDE
            DBDU(I,J)=0.0D0
            DBDUX(I,J)=0.0D0
          enddo
        enddo

        if(IBN.EQ.2)THEN       ! CLOSED CELL

C         fix n vector length
          UN=dsqrt(U(4)**2+U(5)**2+U(6)**2)
          UNx=U(4)/UN
          UNy=U(5)/UN
          UNz=U(6)/UN

          ST=dsin(U(7))
          ST2=2.0D0*ST
          CT=dcos(U(7))
          CTM=1.0D0-CT
          CTM2=2.0D0*CTM
          DD45=UNx*UX(5)-UX(4)*UNy
          FTN=CTM*DD45-ST*UX(6)-UX(7)*UNz
          CTF=CTM*FTN
          STF=ST*FTN
          FTN4=CTM*UX(5)
          FTN5=-CTM*UX(4)
          FTN7=ST*DD45-CT*UX(6)
          FTNX4=-CTM*UNy
          FTNX5=CTM*UNx
          C46=CTM*UNx*UNz+UNy*ST
          C56=CTM*UNy*UNz-UNx*ST             !!!!!!!!!!!
          C66=CTM*UNz**2+CT
          C266=2.0D0-C66

          W=GAM*H
          AF=-(CPAR0+CPAR_SWR*T)**2/W

          AF=AF*(1D0 - 0.5D0 * AER_STEP(X,0))

          DA=-(DF0+DF_SWR*T)/AF

          DBDUX(4,1)=DA
          DBDUX(5,2)=DA
          DBDUX(6,3)=DA

          DBDU(4,4)=2.0D0*UX(7)+CTF*UNz+C46*FTN4
          DBDU(4,5)=CTM2*UX(6)+STF+C46*FTN5
          DBDU(4,6)=-CTM2*UX(5)+CTF*UNx-C46*UX(7)
          DBDU(4,7)=2.0D0*(CT*UX(4)+ST*(UNy*UX(6)-UX(5)*UNz))+
     *     STF*UNx*UNz+UNy*CT*FTN+C46*FTN7

          DBDU(5,4)=-CTM2*UX(6)-STF+C56*FTN4
          DBDU(5,5)=2.0D0*UX(7)+CTF*UNz+C56*FTN5
          DBDU(5,6)=CTM2*UX(4)+CTF*UNy-C56*UX(7)
          DBDU(5,7)=2.0D0*(CT*UX(5)-ST*(UNx*UX(6)-UX(4)*UNz))+
     *     STF*UNy*UNz-UNx*CT*FTN+C56*FTN7

          DBDU(6,4)=CTM2*UX(5)+C66*FTN4
          DBDU(6,5)=-CTM2*UX(4)+C66*FTN5
          DBDU(6,6)=2.0D0*UNz*CTF+C266*UX(7)
          DBDU(6,7)=2.0D0*(CT*UX(6)+ST*DD45)+STF*(UNz**2-1.0D0)+C66*FTN7

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
C          DBDU(7,4)=UX(4)         !!
C          DBDU(7,5)=UX(5)         !!
C          DBDU(7,6)=UX(6)         !!
C          DBDUX(7,4)=UNx         !!
C          DBDUX(7,5)=UNy         !!
C          DBDUX(7,6)=UNz         !!
        endif
        return
      end

C-- UINIT ------ INITIAL CONDITIONS ---------------------------------
      subroutine UINIT(XI,UI,NPDEI)
        include 'vmcw_pars.fh'
        include 'he3_const.fh'
        real*8 USOL, XSOL
        real*8 XI, UI
        integer NPDEI
        dimension UI(NPDEI)

        real*8 BET, DELTA, DELTAX, DELTAY,
     .         UCTG, UNX,UNY,UNZ, UMX,UMY,UMZ,
     .         UNZ2, USM
        integer I,J,K

        if (NPDEI.ne.7) then
          write(*,*) 'Error in UINIT: NPDEI!=7'
          stop
        endif

        BET=BETA*PI/180.0D0
        UMZ=dcos(BET)
        UMX=dsin(BET)*dsqrt(0.5D0)
        UMY=UMX
        if (UMZ.GE.-0.25D0)THEN
          UNZ2=0.8D0*(0.25D0+dcos(BET))
          UNZ=dsqrt(UNZ2)
          DELTA=(25.0D0*UNZ2+15.0D0)/16.0D0
          DELTAX=dsin(BET)*dsqrt(0.5D0)*
     *     (UNZ*1.25D0-dsqrt(15.0D0)*0.25D0)
          DELTAY=dsin(BET)*dsqrt(0.5D0)*
     *     (UNZ*1.25D0+dsqrt(15.0D0)*0.25D0)
          UNX=DELTAX/DELTA
          UNY=DELTAY/DELTA
          UCTG=dacos(-0.25D0)
        else
          UNZ=0.0D0
          UNX=-dsqrt(0.5D0)
          UNY=dsqrt(0.5D0)
          UCTG=BET
        endif

        USM=dsqrt(UNX**2+UNY**2+UNZ**2)
        UI(1)=UMX             ! Mx
        UI(2)=UMY             ! My
        UI(3)=UMZ             ! Mz
        UI(4)=UNX/USM         ! Nx   !!
        UI(5)=UNY/USM         ! Ny   !!
        UI(6)=UNZ/USM         ! Nz
        UI(7)=UCTG            ! TETA         !!!!!+/-
        return
      end
C-- DERIVF ----- SET UP DERIVATIVES ---------------------------------
      subroutine DERIVF(T,X,U,UX,UXX,DFDU,DFDUX,DFDUXX,NPDE)
        include 'vmcw_pars.fh'
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
CCC ATTENTION :   DERIVF IS WRONG
        return
      end
