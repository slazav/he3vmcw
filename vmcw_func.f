! Leggett equations in rotating frame. 1D calculations along rotation axis.
! Spin currents, spin diffusion, leggett-takagi, extra Mz relaxation.
! input:
!   T - time
!   X - x-coord
!   U(6)   - Mx My Mz Nx*Th Ny*Th Nz*Th
!   Ux(6)  - dU/dx
!   Uxx(6) - d2U/dx2
!   NPDE=6
! output
!   Fv(6)  - dU/dt
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

        real*8 UMx,UMy,UMz, UNx,UNy,UNz, Uth, UDU
        real*8 UMxm,UMym,UMzm
        real*8 GNx,GNy,GNz, GMx,GMy,GMz, Gth
        real*8 GGNx,GGNy,GGNz, GGMx,GGMy,GGMz, GGth
        real*8 dNx, dNy, dNz, dTh
        real*8 DD45,ST,CT,CTM,CT1,CTG,UT,AUT,AF,DAF,FTN,DFTN,B
        real*8 UJX,UJY,UJZ,DJX,DJY,DJZ ! spin current J and dJ/dz

!       set parameters:
        call set_bulk_pars(T,X, Wr,Wz,W0,WB,
     *                     Cpar, dCpar, Diff, Tf, T1)

        UMx = U(1)
        UMy = U(2)
        UMz = U(3)
        UMxm = UMx-Wr/Wz
        UMym = UMy
        UMzm = UMz-1.0D0
!        GMx = UX(1)
!        GMy = UX(2)
!        GMz = UX(3)
        GGMx = UXX(1)
        GGMy = UXX(2)
        GGMz = UXX(3)


C       U(4,5,6) = n*th
C       th = |U|
C       n = U/|U|
        Uth=dsqrt(U(4)**2+U(5)**2+U(6)**2)
        UNx = U(4)/Uth
        UNy = U(5)/Uth
        UNz = U(6)/Uth

C       dTh = d|U| = (U*dU)/|U|
C       dn  = dU/|U| - U d|U|/|U|^2
        Gth = (U(4)*UX(4)+U(5)*UX(5)+U(6)*UX(6))/Uth
        GNx = UX(4)/Uth - U(4)*Gth/Uth**2
        GNy = UX(5)/Uth - U(5)*Gth/Uth**2
        GNz = UX(6)/Uth - U(6)*Gth/Uth**2

C       d2|U| = (dU*dU)/|U| + (U*d2U)/|U| - (U*dU)^2/|U|^3
        GGth = (UX(4)**2+UX(5)**2+U(6)**2
     *         +U(4)*UXX(4)+U(5)*UXX(5)+U(6)*UXX(6)
     *         - Gth**2)/Uth

C       d2n = d2U/|U| - 2 dU d|U| / |U|^2
C           - U d2|U|/|U|^2 + 2 U d|U|^2/|U|^3
        GGNx = UXX(4)/Uth - 2D0*UX(4)*Gth/Uth**2
     *         - U(4)*GGth/Uth**2 + 2D0*U(4)*Gth**2/Uth**3
        GGNy = UXX(5)/Uth - 2D0*UX(5)*Gth/Uth**2
     *         - U(5)*GGth/Uth**2 + 2D0*U(5)*Gth**2/Uth**3
        GGNz = UXX(6)/Uth - 2D0*UX(6)*Gth/Uth**2
     *         - U(6)*GGth/Uth**2 + 2D0*U(6)*Gth**2/Uth**3

        ST=dsin(Uth)
        CT=dcos(Uth)

        if (dabs(ST).lt.1D-6)
     *      write (*,*) "Warning: ST: ", ST

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
        dNx = -W0*UNy-0.5D0*Wz*(UMzm*UNy-UMym*UNz+CTG*(B*UNx-UMxm))
        dNy =  W0*UNx-0.5D0*Wz*(UMxm*UNz-UMzm*UNx+CTG*(B*UNy-UMym))
        dNz =        -0.5D0*Wz*(UMym*UNx-UMxm*UNy+CTG*(B*UNz-UMzm))
        dTh =  Wz*B + UT/Tf

        FV(4) = dNx*Uth + dTh*UNx
        FV(5) = dNy*Uth + dTh*UNy
        FV(6) = dNz*Uth + dTh*UNz

        return
      end

C-- BNDRY ------ BOUNDARY CONDITIONS -- B(U,UX)=Z(T) ------------

      subroutine BNDRY(T,X,U,UX,DBDU,DBDUX,DZDT,NPDE)
        implicit none
        integer NPDE
        dimension U(NPDE),UX(NPDE),DZDT(NPDE),
     *   DBDU(NPDE,NPDE),DBDUX(NPDE,NPDE)
        real*8 T,X,U,UX,DZDT,DBDU,DBDUX
        real*8 UNx,UNy,UNz, Uth, UGU
        real*8 GNx,GNy,GNz, GMx,GMy,GMz, Gth
        real*8 GGNx,GGNy,GGNz, GGMx,GGMy,GGMz, GGth
        real*8 ST,ST2,CT,CTM,CTM2,DD45,FTN,CTF,STF
        real*8 FTN4,FTN5,FTN7,FTNX4,FTNX5,C46,C56,C66,C266
        real*8 dBdnx,dBdny,dBdnz,dBdth
        real*8 dBdgnx,dBdgny,dBdgnz,dBdgth
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
        if (IBN.EQ.1) THEN
          do J=1,NPDE
            DBDUX(J,J)=1.0D0
          enddo
          return
        endif

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

          Uth=dsqrt(U(4)**2+U(5)**2+U(6)**2)
          UNx = U(4)/Uth
          UNy = U(5)/Uth
          UNz = U(6)/Uth

C         d|U| = (U*dU)/|U|
C         dn = dU/|U| - U d|U|/|U|^2
          Gth = (U(4)*UX(4)+U(5)*UX(5)+U(6)*UX(6))/Uth
          GNx = UX(4)/UTh - U(4)*Gth/Uth**2
          GNy = UX(5)/UTh - U(5)*Gth/Uth**2
          GNz = UX(6)/UTh - U(6)*Gth/Uth**2

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

!          U = n th
!          GU = GN Uth + UN GTh

!         B(4)
!          dBdnx = Gth + (1-CTM*UNz**2)*Gth + 2*ST*CTM*UNz*GNz
!          dBdny = CTM*GNz - ST*UNz*Gth - CTM*(CT-CTM*UNz**2)*GNz
!          dBdnz = -CTM*GNy - 2D0*CTM*UNz*UNx*Gth - ST*UNy*Gth
!     *          + 2D0*ST*CTM*(UNz*GNx-UNx*GNz)
!     *          + 2D0*CTM**2*(UNy*GNz-UNz*GNy)*UNz
!     *          - CTM*(CT+CTM*UNz**2)*GNy
!          dBdth = CT*GNx + ST *(UNy*GNz - UNz*GNy)
!     *          - ST*UNz**2*UNx*Gth - CT*UNy*UNz*Gth
!     *          + CT*(CT+CTM*UNz**2)*GNx
!     *          - ST*(CT+CTM*UNz**2)*UNz*GNy
!     *          - ST*(CT-CTM*UNz**2)*UNy*GNz
!     *          - 2D0*CT*CTM*UNx*UNz*GNz
!     *          -  ST*ST*(1D0-UNz**2)*GNx
!     *          + CTM*ST*(1D0-UNz**2)*UNz*GNy
!     *          + CTM*ST*(1D0+UNz**2)*UNy*GNz
!     *          - 2D0*ST**2*UNx*UNz*GNz
!          dBdgnx =   ST*(1D0 + CT + CTM*UNz**2)
!          dBdgny = -CTM*(1D0 + CT + CTM*UNz**2)*UNz
!          dBdgnz =  CTM*(1D0 - CT + CTM*UNz**2)*UNy
!     *           - 2D0*ST*CTM*UNx*UNz*GNz
!          dBdgth = UNx + (1D0-CTM*UNz**2)*UNx - ST*UNy*UNz

          dBdnx =  2.0D0*Gth+CTF*UNz+C46*FTN4
          dBdny =  CTM2*GNz+STF+C46*FTN5
          dBdnz = -CTM2*GNy+CTF*UNx-C46*Gth
          dBdth =  2.0D0*(CT*GNx+ST*(UNy*GNz-GNy*UNz))+
     *             STF*UNx*UNz+UNy*CT*FTN+C46*FTN7
          dBdgnx =  ST2+C46*FTNX4
          dBdgny = -CTM2*UNz+C46*FTNX5
          dBdgnz =  CTM2*UNy-C46*ST
          dBdgth =  2.0D0*UNx-C46*UNz

!  qi = ni*th
!  th = |q|
!  ni = qi/|q|
! d(th)/dqi  = ni
! d(nj)/dqi  = (dij - ni nj)/th
!
!  Gth = (Gq*q)/|q|
!  Gni = Gqi/|q| - qi(Gq*q)/|q|^3
! d(Gth)/dqi  = Gni
! d(Gni)/dqj  = (ni nj - dij)Gth/th^2 -(nj Gni + ni Gnj)/th + 
! d(Gth)/dGqi = ni
! d(Gni)/dGqj = (dij - ni nj)/th
          DBDU(4,4) = dBdth*UNx + dBdnx*(1D0-Unx*Unx)/Uth
     *              - dBdny*Unx*Uny/Uth - dBdnz*Unx*Unz/Uth
     *              + dBdgth*GNx
     *              + dBdgnx*(Gth*(UNx**2-1D0)/Uth-2D0*UNx*GNx)/Uth
     *              + dBdgny*(Gth*(UNx*UNy)/Uth-UNx*GNy-UNy*GNx)/Uth
     *              + dBdgnz*(Gth*(UNx*UNz)/Uth-UNx*GNz-UNz*GNx)/Uth
          DBDU(4,5) = dBdth*UNy + dBdny*(1D0-Uny*Uny)/Uth
     *              - dBdnx*Uny*Unx/Uth - dBdnz*Uny*Unz/Uth
     *              + dBdgth*GNy
     *              + dBdgnx*(Gth*(UNy*UNx)/Uth-UNy*GNx-UNx*GNy)/Uth
     *              + dBdgny*(Gth*(UNy**2-1D0)/Uth-2D0*UNy*GNy)/Uth
     *              + dBdgnz*(Gth*(UNy*UNz)/Uth-UNy*GNz-UNz*GNy)/Uth
          DBDU(4,6) = dBdth*UNz + dBdnz*(1D0-Unz*Unz)/Uth
     *              - dBdnx*Unz*Unx/Uth - dBdny*Unz*Uny/Uth
     *              + dBdgth*GNz
     *              + dBdgnx*(Gth*(UNz*UNx)/Uth-UNz*GNx-UNx*GNz)/Uth
     *              + dBdgny*(Gth*(UNz*UNy)/Uth-UNz*GNy-UNy*GNz)/Uth
     *              + dBdgnz*(Gth*(UNz**2-1D0)/Uth-2D0*UNz*GNz)/Uth
          DBDUX(4,4) = dBdgth*UNx + dBdgnx*(1D0-Unx*Unx)/Uth
     *               - dBdgny*Unx*Uny/Uth - dBdgnz*Unx*Unz/Uth
          DBDUX(4,5) = dBdgth*UNy + dBdgny*(1D0-Uny*Uny)/Uth
     *               - dBdgnx*Uny*Unx/Uth - dBdgnz*Uny*Unz/Uth
          DBDUX(4,6) = dBdgth*UNz + dBdgnz*(1D0-Unz*Unz)/Uth
     *               - dBdgnx*Unz*Unx/Uth - dBdgny*Unz*Uny/Uth
!         write (*,*) DBDU(4,4), DBDU(4,5), DBDU(4,6)
!         write (*,*) DBDUX(4,4), DBDUX(4,5), DBDUX(4,6)

!         B(5)
!          dBdnx = -CTM*GNz + ST*UNz*Gth + CTM*(CT-CTM*UNz**2)*GNz
!          dBdny = Gth + (1-CTM*UNz**2)*Gth - 2*ST*CTM*UNz*GNz
!          dBdnz = CTM*GNx - 2D0*CTM*UNz*UNy*Gth + ST*UNx*Gth
!     *          + 2D0*ST*CTM*(UNz*GNy-UNy*GNz)
!     *          + 2D0*CTM**2*(UNz*GNx-UNx*GNz)*UNz
!     *          + CTM*(CT+CTM*UNz**2)*GNx
!          dBdth = CT*GNy + ST*(UNz*GNx - UNx*GNz)
!     *          - ST*UNz**2*UNy*Gth + CT*UNx*UNz*Gth
!     *          + CT*(CT+CTM*UNz**2)*GNy
!     *          + ST*(CT+CTM*UNz**2)*UNz*GNx
!     *          + ST*(CT-CTM*UNz**2)*UNx*GNz
!     *          - 2D0*CT*CTM*UNy*UNz*GNz
!     *          -  ST*ST*(1D0-UNz**2)*GNy
!     *          - CTM*ST*(1D0-UNz**2)*UNz*GNx
!     *          - CTM*ST*(1D0+UNz**2)*UNx*GNz
!     *          - 2D0*ST**2*UNy*UNz*GNz
!          dBdgnx = CTM*(1D0+CT+CTM*UNz**2)*UNz
!          dBdgny =  ST*(1D0+CT+CTM*UNz**2)
!          dBdgnz =-CTM*(1D0-CT+CTM*UNz**2)*UNx
!     *           - 2D0*ST*CTM*UNy*UNz
!          dBdgth = UNy + (1D0-CTM*UNz**2)*UNy + ST*UNx*UNz

          dBdnx = -CTM2*GNz-STF+C56*FTN4
          dBdny =  2.0D0*Gth+CTF*UNz+C56*FTN5
          dBdnz =  CTM2*GNx+CTF*UNy-C56*Gth
          dBdth =  2.0D0*(CT*GNy-ST*(UNx*GNz-GNx*UNz))+
     *             STF*UNy*UNz-UNx*CT*FTN+C56*FTN7
          dBdgnx = CTM2*UNz+C56*FTNX4
          dBdgny = ST2+C56*FTNX5
          dBdgnz = -CTM2*UNx-C56*ST
          dBdgth = 2.0D0*UNy-C56*UNz

          DBDU(5,4) = dBdth*UNx + dBdnx*(1D0-Unx*Unx)/Uth
     *              - dBdny*Unx*Uny/Uth - dBdnz*Unx*Unz/Uth
     *              + dBdgth*GNx
     *              + dBdgnx*(Gth*(UNx**2-1D0)/Uth-2D0*UNx*GNx)/Uth
     *              + dBdgny*(Gth*(UNx*UNy)/Uth-UNx*GNy-UNy*GNx)/Uth
     *              + dBdgnz*(Gth*(UNx*UNz)/Uth-UNx*GNz-UNz*GNx)/Uth
          DBDU(5,5) = dBdth*UNy + dBdny*(1D0-Uny*Uny)/Uth
     *              - dBdnx*Uny*Unx/Uth - dBdnz*Uny*Unz/Uth
     *              + dBdgth*GNy
     *              + dBdgnx*(Gth*(UNy*UNx)/Uth-UNy*GNx-UNx*GNy)/Uth
     *              + dBdgny*(Gth*(UNy**2-1D0)/Uth-2D0*UNy*GNy)/Uth
     *              + dBdgnz*(Gth*(UNy*UNz)/Uth-UNy*GNz-UNz*GNy)/Uth
          DBDU(5,6) = dBdth*UNz + dBdnz*(1D0-Unz*Unz)/Uth
     *              - dBdnx*Unz*Unx/Uth - dBdny*Unz*Uny/Uth
     *              + dBdgth*GNz
     *              + dBdgnx*(Gth*(UNz*UNx)/Uth-UNz*GNx-UNx*GNz)/Uth
     *              + dBdgny*(Gth*(UNz*UNy)/Uth-UNz*GNy-UNy*GNz)/Uth
     *              + dBdgnz*(Gth*(UNz**2-1D0)/Uth-2D0*UNz*GNz)/Uth
          DBDUX(5,4) = dBdgth*UNx + dBdgnx*(1D0-Unx*Unx)/Uth
     *               - dBdgny*Unx*Uny/Uth - dBdgnz*Unx*Unz/Uth
          DBDUX(5,5) = dBdgth*UNy + dBdgny*(1D0-Uny*Uny)/Uth
     *               - dBdgnx*Uny*Unx/Uth - dBdgnz*Uny*Unz/Uth
          DBDUX(5,6) = dBdgth*UNz + dBdgnz*(1D0-Unz*Unz)/Uth
     *               - dBdgnx*Unz*Unx/Uth - dBdgny*Unz*Uny/Uth

!          write (*,*) DBDU(5,4), DBDU(5,5), DBDU(5,6)
!          write (*,*) DBDUX(5,4), DBDUX(5,5), DBDUX(5,6)

!         B(6)
!          dBdnx = CTM*GNy + (ST**2+CTM**2*UNz**2)*GNy
!          dBdny = -CTM*GNx - (ST*ST + CTM**2 *UNz**2)*GNx
!          dBdnz = Gth - 2D0*CTM*UNz**2*Gth
!     *          - 2D0*ST*CTM*UNz*GNz
!     *          + CTM*(1D0 - UNz**2)*Gth
!     *          + 2D0*CTM**2*(UNx*GNy-UNy*GNx)*UNz
!          dBdth = CT*GNz + ST*(UNx*GNy-UNy*GNx)
!     *          + ST*(1D0-UNz**2)*UNz*Gth
!     *          - (2D0*ST*CT + 2D0*CTM*ST*UNz**2)*(UNy*GNx-UNx*GNy)
!     *          + (CT*CTM+ST*ST) * (1-UNz**2) *GNz
!          dBdgnx = - (CTM + ST*ST + CTM**2*UNz**2) *UNy
!          dBdgny =   (CTM + ST*ST + CTM**2*UNz**2) *UNx
!          dBdgnz = ST + ST*CTM*(1-UNz**2)
!          dBdgth = UNz + CTM*(1D0-UNz**2)*UNz

          dBdnx =  CTM2*GNy+C66*FTN4
          dBdny = -CTM2*GNx+C66*FTN5
          dBdnz =  2.0D0*UNz*CTF+C266*Gth
          dBdth =  2.0D0*(CT*GNz+ST*DD45)+
     *             STF*(UNz**2-1.0D0)+C66*FTN7
          dBdgnx= -CTM2*UNy+C66*FTNX4
          dBdgny=  CTM2*UNx+C66*FTNX5
          dBdgnz=  C266*ST
          dBdgth=  C266*UNz

          DBDU(6,4) = dBdth*UNx + dBdnx*(1D0-Unx*Unx)/Uth
     *              - dBdny*Unx*Uny/Uth - dBdnz*Unx*Unz/Uth
     *              + dBdgth*GNx
     *              + dBdgnx*(Gth*(UNx**2-1D0)/Uth-2D0*UNx*GNx)/Uth
     *              + dBdgny*(Gth*(UNx*UNy)/Uth-UNx*GNy-UNy*GNx)/Uth
     *              + dBdgnz*(Gth*(UNx*UNz)/Uth-UNx*GNz-UNz*GNx)/Uth
          DBDU(6,5) = dBdth*UNy + dBdny*(1D0-Uny*Uny)/Uth
     *              - dBdnx*Uny*Unx/Uth - dBdnz*Uny*Unz/Uth
     *              + dBdgth*GNy
     *              + dBdgnx*(Gth*(UNy*UNx)/Uth-UNy*GNx-UNx*GNy)/Uth
     *              + dBdgny*(Gth*(UNy**2-1D0)/Uth-2D0*UNy*GNy)/Uth
     *              + dBdgnz*(Gth*(UNy*UNz)/Uth-UNy*GNz-UNz*GNy)/Uth
          DBDU(6,6) = dBdth*UNz + dBdnz*(1D0-Unz*Unz)/Uth
     *              - dBdnx*Unz*Unx/Uth - dBdny*Unz*Uny/Uth
     *              + dBdgth*GNz
     *              + dBdgnx*(Gth*(UNz*UNx)/Uth-UNz*GNx-UNx*GNz)/Uth
     *              + dBdgny*(Gth*(UNz*UNy)/Uth-UNz*GNy-UNy*GNz)/Uth
     *              + dBdgnz*(Gth*(UNz**2-1D0)/Uth-2D0*UNz*GNz)/Uth
          DBDUX(6,4) = dBdgth*UNx + dBdgnx*(1D0-Unx*Unx)/Uth
     *               - dBdgny*Unx*Uny/Uth - dBdgnz*Unx*Unz/Uth
          DBDUX(6,5) = dBdgth*UNy + dBdgny*(1D0-Uny*Uny)/Uth
     *               - dBdgnx*Uny*Unx/Uth - dBdgnz*Uny*Unz/Uth
          DBDUX(6,6) = dBdgth*UNz + dBdgnz*(1D0-Unz*Unz)/Uth
     *               - dBdgnx*Unz*Unx/Uth - dBdgny*Unz*Uny/Uth

!          write (*,*) DBDU(6,4), DBDU(6,5), DBDU(6,6)
!          write (*,*) DBDUX(6,4), DBDUX(6,5), DBDUX(6,6)
!          write(*,*) dBdgth, UNx, dBdgth, UNy

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
