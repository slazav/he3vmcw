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

        real*8 UMx,UMy,UMz, UN(3), Uth, UDU
        real*8 UMxm,UMym,UMzm
        real*8 GN(3), GMx,GMy,GMz, Gth
        real*8 GGN(3), GGMx,GGMy,GGMz, GGth
        real*8 dNx, dNy, dNz, dTh
        real*8 ST,CT,CTGT,UT,AUT,AF,DAF,B
        real*8 UJ(3),DJ(3) ! spin current J and dJ/dz

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
        UN(1) = U(4)/Uth
        UN(2) = U(5)/Uth
        UN(3) = U(6)/Uth

C       dTh = d|U| = (U*dU)/|U|
C       dn  = dU/|U| - U d|U|/|U|^2
        Gth = (U(4)*UX(4)+U(5)*UX(5)+U(6)*UX(6))/Uth
        GN(1) = UX(4)/Uth - U(4)*Gth/Uth**2
        GN(2) = UX(5)/Uth - U(5)*Gth/Uth**2
        GN(3) = UX(6)/Uth - U(6)*Gth/Uth**2

C       d2|U| = (dU*dU)/|U| + (U*d2U)/|U| - (U*dU)^2/|U|^3
        GGth = (UX(4)**2+UX(5)**2+U(6)**2
     *         +U(4)*UXX(4)+U(5)*UXX(5)+U(6)*UXX(6)
     *         - Gth**2)/Uth

C       d2n = d2U/|U| - 2 dU d|U| / |U|^2
C           - U d2|U|/|U|^2 + 2 U d|U|^2/|U|^3
        GGN(1) = UXX(4)/Uth - 2D0*UX(4)*Gth/Uth**2
     *         - U(4)*GGth/Uth**2 + 2D0*U(4)*Gth**2/Uth**3
        GGN(2) = UXX(5)/Uth - 2D0*UX(5)*Gth/Uth**2
     *         - U(5)*GGth/Uth**2 + 2D0*U(5)*Gth**2/Uth**3
        GGN(3) = UXX(6)/Uth - 2D0*UX(6)*Gth/Uth**2
     *         - U(6)*GGth/Uth**2 + 2D0*U(6)*Gth**2/Uth**3

        ST=dsin(Uth)
        CT=dcos(Uth)

        if (dabs(ST).lt.1D-6)
     *      write (*,*) "Warning: ST: ", ST

        if (dabs(Uth).gt.1D-3) then
          CTGT=ST/(1.0D0-CT)*Uth      ! T*ctg(T/2) = sin(T)/(1-cos(T))
        else
          CTGT=2D0-Uth**2/6D0
        endif
        UT=ST*(1.0D0+4.0D0*CT)*0.2666666D0 ! 4/15 sin(t)*(1+4*cos(t))

        AUT=UT*WB**2/Wz
        AF=-Cpar**2/Wz
        DAF = -2D0*Cpar*dCpar/Wz

C       Spin currents
        call fill_jgd_nt(UJ, UN, Uth, GN, Gth)
C       Gradient torque
        call fill_tgd_nt(DJ, UN, Uth, GN, Gth, GGN, GGth)


C       derivative of the spin flow: d/dz(Jiz - Diff*Mi')
C       Note: here non-uniform spin-wave velocity (DAF!=0) can be used,
C             but spin diffusion should be uniform
        DJ(1) = AF*DJ(1) + DAF*UJ(3) - DIFF*GGMx
        DJ(2) = AF*DJ(2) + DAF*UJ(3) - DIFF*GGMy
        DJ(3) = AF*DJ(3) + DAF*UJ(3) - DIFF*GGMz

        B = UN(1)*UMxm + UMym*UN(2) + UMzm*UN(3)

C       Leggett equations
        FV(1)=   (Wz-W0)*UMy          + AUT*UN(1) - DJ(1)
        FV(2)= - (Wz-W0)*UMx + Wr*UMz + AUT*UN(2) - DJ(2)
        FV(3)=               - Wr*UMy + AUT*UN(3) - DJ(3) - UMzm/T1
!        dNx = -W0*UN(2)-0.5D0*Wz*(UMzm*UN(2)-UMym*UN(3)+CTG*(B*UN(1)-UMxm))
!        dNy =  W0*UN(1)-0.5D0*Wz*(UMxm*UN(3)-UMzm*UN(1)+CTG*(B*UN(2)-UMym))
!        dNz =        -0.5D0*Wz*(UMym*UN(1)-UMxm*UN(2)+CTG*(B*UN(3)-UMzm))
!        dTh =  Wz*B + UT/Tf

        FV(4) = -W0*UN(2)*Uth - 0.5D0*Wz*(
     *              (UMzm*UN(2)-UMym*UN(3))*Uth
     *             +CTGT*(B*UN(1)-UMxm)) + (Wz*B + UT/Tf)*UN(1)
        FV(5) = W0*UN(1)*Uth - 0.5D0*Wz*(
     *              (UMxm*UN(3)-UMzm*UN(1))*Uth
     *             +CTGT*(B*UN(2)-UMym)) + (Wz*B + UT/Tf)*UN(2)
        FV(6) = -0.5D0*Wz*(
     *         UMym*UN(1)*Uth
     *        -UMxm*UN(2)*Uth
     *        +CTGT*(B*UN(3)-UMzm)) + (Wz*B + UT/Tf)*UN(3)

        return
      end

C-- BNDRY ------ BOUNDARY CONDITIONS -- B(U,UX)=Z(T) ------------

      subroutine BNDRY(T,X,U,UX,DBDU,DBDUX,DZDT,NPDE)
        implicit none
        integer NPDE
        dimension U(NPDE),UX(NPDE),DZDT(NPDE),
     *   DBDU(NPDE,NPDE),DBDUX(NPDE,NPDE)
        real*8 T,X,U,UX,DZDT,DBDU,DBDUX
        real*8 UJx,UJy,UJz
        real*8 UN(3), Uth, UGU
        real*8 GN(3), GMx,GMy,GMz, Gth
        real*8 GGN(3), GGMx,GGMy,GGMz, GGth
        real*8 ST,ST2,CT,CTM,CTM2,DD45,FTN,CTF,STF
        real*8 FTN4,FTN5,FTN7,FTNX4,FTNX5,C46,C56,C66,C266
        real*8 dBdnx,dBdny,dBdnz,dBdth
        real*8  dBdGNx, dBdGNy, dBdGNz,dBdgth
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

!       Zero derivatives
        if (IBN.EQ.1) THEN
          do J=1,NPDE
            DBDUX(J,J)=1.0D0
          enddo
          return
        endif

!       Zero derivatives of M, fixed n
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
          UN(1) = U(4)/Uth
          UN(2) = U(5)/Uth
          UN(3) = U(6)/Uth

C         d|U| = (U*dU)/|U|
C         dn = dU/|U| - U d|U|/|U|^2
          Gth = (U(4)*UX(4)+U(5)*UX(5)+U(6)*UX(6))/Uth
          GN(1) = UX(4)/UTh - U(4)*Gth/Uth**2
          GN(2) = UX(5)/UTh - U(5)*Gth/Uth**2
          GN(3) = UX(6)/UTh - U(6)*Gth/Uth**2

          ST=dsin(Uth)
          CT=dcos(Uth)
          CTM=1.0D0-CT

          AF=-Cpar**2/W0
          DA=-Diff/AF

          DBDUX(1,1)=DA
          DBDUX(2,2)=DA
          DBDUX(3,3)=DA

!          U = n th
!          GU = GN Uth + UN GTh

!         B(4) == JX
          dBdnx = Gth + (1D0-CTM*UN(3)**2)*Gth
     *          - 2D0*ST*CTM*UN(3)*GN(3)
          dBdny = CTM*GN(3) - ST*UN(3)*Gth
     *          - CTM*(CT-CTM*UN(3)**2)*GN(3)
          dBdnz = -CTM*GN(2) - 2D0*CTM*UN(3)*UN(1)*Gth - ST*UN(2)*Gth
     *          + 2D0*ST*CTM*(UN(3)*GN(1)-UN(1)*GN(3))
     *          + 2D0*CTM**2*(UN(2)*GN(3)-UN(3)*GN(2))*UN(3)
     *          - CTM*(CT+CTM*UN(3)**2)*GN(2)
          dBdth = CT*GN(1) + ST *(UN(2)*GN(3) - UN(3)*GN(2))
     *          - ST*UN(3)**2*UN(1)*Gth - CT*UN(2)*UN(3)*Gth
     *          + CT*(CT+CTM*UN(3)**2)*GN(1)
     *          - ST*(CT+CTM*UN(3)**2)*UN(3)*GN(2)
     *          - ST*(CT-CTM*UN(3)**2)*UN(2)*GN(3)
     *          - 2D0*CT*CTM*UN(1)*UN(3)*GN(3)
     *          -  ST*ST*(1D0-UN(3)**2)*GN(1)
     *          + CTM*ST*(1D0-UN(3)**2)*UN(3)*GN(2)
     *          + CTM*ST*(1D0+UN(3)**2)*UN(2)*GN(3)
     *          - 2D0*ST**2*UN(1)*UN(3)*GN(3)
           dBdGNx =   ST*(1D0 + CT + CTM*UN(3)**2)
           dBdGNy = -CTM*(1D0 + CT + CTM*UN(3)**2)*UN(3)
           dBdGNz =  CTM*(1D0 - CT + CTM*UN(3)**2)*UN(2)
     *           - 2D0*ST*CTM*UN(1)*UN(3)
          dBdgth = UN(1) + (1D0-CTM*UN(3)**2)*UN(1) - ST*UN(2)*UN(3)

          ! see He3B_04_HPD_Q2.tex
          DBDU(1,4) = dBdth*UN(1) + dBdgth*GN(1)
     *  + dBdnx/Uth - UN(1)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(1)**2-1D0)-2D0*UN(1)*GN(1))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(1)*UN(2))-UN(1)*GN(2)-UN(2)*GN(1))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(1)*UN(3))-UN(1)*GN(3)-UN(3)*GN(1))/Uth
          DBDU(1,5) = dBdth*UN(2) + dBdgth*GN(2)
     *  + dBdny/Uth - UN(2)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(2)*UN(1))-UN(2)*GN(1)-UN(1)*GN(2))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(2)**2-1D0)-2D0*UN(2)*GN(2))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(2)*UN(3))-UN(2)*GN(3)-UN(3)*GN(2))/Uth
          DBDU(1,6) = dBdth*UN(3) + dBdgth*GN(3)
     *  + dBdnz/Uth - UN(3)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(3)*UN(1))-UN(3)*GN(1)-UN(1)*GN(3))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(3)*UN(2))-UN(3)*GN(2)-UN(2)*GN(3))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(3)**2-1D0)-2D0*UN(3)*GN(3))/Uth
          DBDUX(1,4) = dBdgth*UN(1) + dBdGNx/Uth
     *    - UN(1)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth
          DBDUX(1,5) = dBdgth*UN(2) +  dBdGNy/Uth
     *    - UN(2)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth
          DBDUX(1,6) = dBdgth*UN(3) +  dBdGNz/Uth
     *    - UN(3)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth
!         write (*,*) DBDU(1,4), DBDU(1,5), DBDU(1,6)
!         write (*,*) DBDUX(1,4), DBDUX(1,5), DBDUX(1,6)

!         B(2) == JY
          dBdnx = -CTM*GN(3) + ST*UN(3)*Gth
     *            +CTM*(CT-CTM*UN(3)**2)*GN(3)
          dBdny = Gth + (1D0-CTM*UN(3)**2)*Gth
     *          - 2D0*ST*CTM*UN(3)*GN(3)
          dBdnz = CTM*GN(1) - 2D0*CTM*UN(3)*UN(2)*Gth + ST*UN(1)*Gth
     *          + 2D0*ST*CTM*(UN(3)*GN(2)-UN(2)*GN(3))
     *          + 2D0*CTM**2*(UN(3)*GN(1)-UN(1)*GN(3))*UN(3)
     *          + CTM*(CT+CTM*UN(3)**2)*GN(1)
          dBdth = CT*GN(2) + ST*(UN(3)*GN(1) - UN(1)*GN(3))
     *          - ST*UN(3)**2*UN(2)*Gth + CT*UN(1)*UN(3)*Gth
     *          + CT*(CT+CTM*UN(3)**2)*GN(2)
     *          + ST*(CT+CTM*UN(3)**2)*UN(3)*GN(1)
     *          + ST*(CT-CTM*UN(3)**2)*UN(1)*GN(3)
     *          - 2D0*CT*CTM*UN(2)*UN(3)*GN(3)
     *          -  ST*ST*(1D0-UN(3)**2)*GN(2)
     *          - CTM*ST*(1D0-UN(3)**2)*UN(3)*GN(1)
     *          - CTM*ST*(1D0+UN(3)**2)*UN(1)*GN(3)
     *          - 2D0*ST**2*UN(2)*UN(3)*GN(3)
           dBdGNx = CTM*(1D0+CT+CTM*UN(3)**2)*UN(3)
           dBdGNy =  ST*(1D0+CT+CTM*UN(3)**2)
           dBdGNz =-CTM*(1D0-CT+CTM*UN(3)**2)*UN(1)
     *           - 2D0*ST*CTM*UN(2)*UN(3)
          dBdgth = UN(2) + (1D0-CTM*UN(3)**2)*UN(2) + ST*UN(1)*UN(3)

          DBDU(2,4) = dBdth*UN(1) + dBdgth*GN(1)
     *  + dBdnx/Uth - UN(1)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(1)**2-1D0)-2D0*UN(1)*GN(1))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(1)*UN(2))-UN(1)*GN(2)-UN(2)*GN(1))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(1)*UN(3))-UN(1)*GN(3)-UN(3)*GN(1))/Uth
          DBDU(2,5) = dBdth*UN(2) + dBdgth*GN(2)
     *  + dBdny/Uth - UN(2)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(2)*UN(1))-UN(2)*GN(1)-UN(1)*GN(2))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(2)**2-1D0)-2D0*UN(2)*GN(2))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(2)*UN(3))-UN(2)*GN(3)-UN(3)*GN(2))/Uth
          DBDU(2,6) = dBdth*UN(3) + dBdgth*GN(3)
     *  + dBdnz/Uth - UN(3)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(3)*UN(1))-UN(3)*GN(1)-UN(1)*GN(3))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(3)*UN(2))-UN(3)*GN(2)-UN(2)*GN(3))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(3)**2-1D0)-2D0*UN(3)*GN(3))/Uth
          DBDUX(2,4) = dBdgth*UN(1) + dBdGNx/Uth
     *    - UN(1)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth
          DBDUX(2,5) = dBdgth*UN(2) + dBdGNy/Uth
     *    - UN(2)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth
          DBDUX(2,6) = dBdgth*UN(3) + dBdGNz/Uth
     *    - UN(3)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth

!          write (*,*) DBDU(2,4), DBDU(2,5), DBDU(2,6)
!          write (*,*) DBDUX(2,4), DBDUX(2,5), DBDUX(2,6)

!         B(3) == JZ
          dBdnx = CTM*GN(2) + (ST**2+CTM**2*UN(3)**2)*GN(2)
          dBdny = -CTM*GN(1) - (ST**2 + CTM**2*UN(3)**2)*GN(1)
          dBdnz = Gth - 2D0*CTM*UN(3)**2*Gth
     *        - 2D0*ST*CTM*UN(3)*GN(3)
     *        + CTM*(1D0 - UN(3)**2)*Gth
     *        + 2D0*CTM**2*(UN(1)*GN(2)-UN(2)*GN(1))*UN(3)
          dBdth = CT*GN(3) + ST*(UN(1)*GN(2)-UN(2)*GN(1))
     *        + ST*(1D0-UN(3)**2)*UN(3)*Gth
     *        - 2D0*ST*(CT + CTM*UN(3)**2)*(UN(2)*GN(1)-UN(1)*GN(2))
     *        + (CT*CTM+ST*ST) * (1D0-UN(3)**2) *GN(3)
           dBdGNx = - (CTM + ST*ST + CTM**2*UN(3)**2) *UN(2)
           dBdGNy =   (CTM + ST*ST + CTM**2*UN(3)**2) *UN(1)
           dBdGNz = ST + ST*CTM*(1D0-UN(3)**2)
          dBdgth = UN(3) + CTM*(1D0-UN(3)**2)*UN(3)

          DBDU(3,4) = dBdth*UN(1) + dBdgth*GN(1)
     *  + dBdnx/Uth - UN(1)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(1)**2-1D0)-2D0*UN(1)*GN(1))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(1)*UN(2))-UN(1)*GN(2)-UN(2)*GN(1))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(1)*UN(3))-UN(1)*GN(3)-UN(3)*GN(1))/Uth
          DBDU(3,5) = dBdth*UN(2) + dBdgth*GN(2)
     *  + dBdny/Uth - UN(2)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(2)*UN(1))-UN(2)*GN(1)-UN(1)*GN(2))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(2)**2-1D0)-2D0*UN(2)*GN(2))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(2)*UN(3))-UN(2)*GN(3)-UN(3)*GN(2))/Uth
          DBDU(3,6) = dBdth*UN(3) + dBdgth*GN(3)
     *  + dBdnz/Uth - UN(3)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(3)*UN(1))-UN(3)*GN(1)-UN(1)*GN(3))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(3)*UN(2))-UN(3)*GN(2)-UN(2)*GN(3))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(3)**2-1D0)-2D0*UN(3)*GN(3))/Uth
          DBDUX(3,4) = dBdgth*UN(1) + dBdGNx/Uth
     *    - UN(1)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth
          DBDUX(3,5) = dBdgth*UN(2) + dBdGNy/Uth
     *    - UN(2)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth
          DBDUX(3,6) = dBdgth*UN(3) + dBdGNz/Uth
     *    - UN(3)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth

!         Test spin currents
!          DD45=UN(1)*GN(2)-GN(1)*UN(2)
!          FTN=CTM*DD45-ST*GN(3)-Gth*UN(3)
!          UJX=2.0D0*(Gth*UN(1)+ST*GN(1)+CTM*(UN(2)*GN(3)-GN(2)*UN(3)))
!     *     + (CTM*UN(1)*UN(3)+UN(2)*ST)*FTN + DA*UX(1)
!          UJY=2.0D0*(Gth*UN(2)+ST*GN(2)-CTM*(UN(1)*GN(3)-GN(1)*UN(3)))
!     *     + (CTM*UN(2)*UN(3)-UN(1)*ST)*FTN + DA*UX(2)
!          UJZ=2.0D0*(Gth*UN(3)+ST*GN(3)+CTM*(UN(1)*GN(2)-GN(1)*UN(2)))
!     *     + (CTM*UN(3)**2+CT)*FTN + DA*UX(3)
!          if (dabs(UJX).gt.1D-6) write(*,*) "UJX: ", UJX
!          if (dabs(UJY).gt.1D-6) write(*,*) "UJY: ", UJY
!          if (dabs(UJZ).gt.1D-6) write(*,*) "UJZ: ", UJZ

!          write (*,*) DBDU(3,4), DBDU(3,5), DBDU(3,6)
!          write (*,*) DBDUX(3,4), DBDUX(3,5), DBDUX(3,6)
!          write(*,*) dBdgth, UN(1), dBdgth, UN(2)
          return
        endif


!       Closed cell: no spin flow through walls
!       Dmitriev's version
!       Jiz - Diff Mi' = 0
        if(IBN.EQ.4)THEN       ! CLOSED CELL

          Uth=dsqrt(U(4)**2+U(5)**2+U(6)**2)
          UN(1) = U(4)/Uth
          UN(2) = U(5)/Uth
          UN(3) = U(6)/Uth

C         d|U| = (U*dU)/|U|
C         dn = dU/|U| - U d|U|/|U|^2
          Gth = (U(4)*UX(4)+U(5)*UX(5)+U(6)*UX(6))/Uth
          GN(1) = UX(4)/UTh - U(4)*Gth/Uth**2
          GN(2) = UX(5)/UTh - U(5)*Gth/Uth**2
          GN(3) = UX(6)/UTh - U(6)*Gth/Uth**2

          ST=dsin(Uth)
          ST2=2.0D0*ST
          CT=dcos(Uth)
          CTM=1.0D0-CT
          CTM2=2.0D0*CTM

          DD45=UN(1)*GN(2)-GN(1)*UN(2)
          FTN=CTM*DD45-ST*GN(3)-Gth*UN(3)

          CTF=CTM*FTN
          STF=ST*FTN
          FTN4=CTM*GN(2)
          FTN5=-CTM*GN(1)
          FTN7=ST*DD45-CT*GN(3)
          FTNX4=-CTM*UN(2)
          FTNX5=CTM*UN(1)
          C46=CTM*UN(1)*UN(3)+UN(2)*ST
          C56=CTM*UN(2)*UN(3)-UN(1)*ST             !!!!!!!!!!!
          C66=CTM*UN(3)**2+CT
          C266=2.0D0-C66

          AF=-Cpar**2/W0
          DA=-Diff/AF

          DBDUX(1,1)=DA
          DBDUX(2,2)=DA
          DBDUX(3,3)=DA

!          U = n th
!          GU = GN Uth + UN GTh

!         B(4) == JX
          dBdnx =  2.0D0*Gth+CTF*UN(3)+C46*FTN4
          dBdny =  CTM2*GN(3)+STF+C46*FTN5
          dBdnz = -CTM2*GN(2)+CTF*UN(1)-C46*Gth
          dBdth =  2.0D0*(CT*GN(1)+ST*(UN(2)*GN(3)-GN(2)*UN(3)))+
     *             STF*UN(1)*UN(3)+UN(2)*CT*FTN+C46*FTN7
           dBdGNx =  ST2+C46*FTNX4
           dBdGNy = -CTM2*UN(3)+C46*FTNX5
           dBdGNz =  CTM2*UN(2)-C46*ST
          dBdgth =  2.0D0*UN(1)-C46*UN(3)

          ! see He3B_04_HPD_Q2.tex
          DBDU(1,4) = dBdth*UN(1) + dBdgth*GN(1)
     *  + dBdnx/Uth - UN(1)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(1)**2-1D0)-2D0*UN(1)*GN(1))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(1)*UN(2))-UN(1)*GN(2)-UN(2)*GN(1))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(1)*UN(3))-UN(1)*GN(3)-UN(3)*GN(1))/Uth
          DBDU(1,5) = dBdth*UN(2) + dBdgth*GN(2)
     *  + dBdny/Uth - UN(2)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(2)*UN(1))-UN(2)*GN(1)-UN(1)*GN(2))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(2)**2-1D0)-2D0*UN(2)*GN(2))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(2)*UN(3))-UN(2)*GN(3)-UN(3)*GN(2))/Uth
          DBDU(1,6) = dBdth*UN(3) + dBdgth*GN(3)
     *  + dBdnz/Uth - UN(3)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(3)*UN(1))-UN(3)*GN(1)-UN(1)*GN(3))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(3)*UN(2))-UN(3)*GN(2)-UN(2)*GN(3))/Uth
     *    +  dBdGNz*(Gth/Uth*(UN(3)**2-1D0)-2D0*UN(3)*GN(3))/Uth
          DBDUX(1,4) = dBdgth*UN(1) +  dBdGNx/Uth
     *    - UN(1)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth
          DBDUX(1,5) = dBdgth*UN(2) +  dBdGNy/Uth
     *    - UN(2)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth
          DBDUX(1,6) = dBdgth*UN(3) +  dBdGNz/Uth
     *    - UN(3)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth
!         write (*,*) DBDU(1,4), DBDU(1,5), DBDU(1,6)
!         write (*,*) DBDUX(1,4), DBDUX(1,5), DBDUX(1,6)

!         B(2) == JY
          dBdnx = -CTM2*GN(3)-STF+C56*FTN4
          dBdny =  2.0D0*Gth+CTF*UN(3)+C56*FTN5
          dBdnz =  CTM2*GN(1)+CTF*UN(2)-C56*Gth
          dBdth =  2.0D0*(CT*GN(2)-ST*(UN(1)*GN(3)-GN(1)*UN(3)))+
     *             STF*UN(2)*UN(3)-UN(1)*CT*FTN+C56*FTN7
           dBdGNx = CTM2*UN(3)+C56*FTNX4
           dBdGNy = ST2+C56*FTNX5
           dBdGNz = -CTM2*UN(1)-C56*ST
          dBdgth = 2.0D0*UN(2)-C56*UN(3)

          DBDU(2,4) = dBdth*UN(1) + dBdgth*GN(1)
     *  + dBdnx/Uth - UN(1)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(1)**2-1D0)-2D0*UN(1)*GN(1))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(1)*UN(2))-UN(1)*GN(2)-UN(2)*GN(1))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(1)*UN(3))-UN(1)*GN(3)-UN(3)*GN(1))/Uth
          DBDU(2,5) = dBdth*UN(2) + dBdgth*GN(2)
     *  + dBdny/Uth - UN(2)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(2)*UN(1))-UN(2)*GN(1)-UN(1)*GN(2))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(2)**2-1D0)-2D0*UN(2)*GN(2))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(2)*UN(3))-UN(2)*GN(3)-UN(3)*GN(2))/Uth
          DBDU(2,6) = dBdth*UN(3) + dBdgth*GN(3)
     *  + dBdnz/Uth - UN(3)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(3)*UN(1))-UN(3)*GN(1)-UN(1)*GN(3))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(3)*UN(2))-UN(3)*GN(2)-UN(2)*GN(3))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(3)**2-1D0)-2D0*UN(3)*GN(3))/Uth
          DBDUX(2,4) = dBdgth*UN(1) +  dBdGNx/Uth
     *    - UN(1)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth
          DBDUX(2,5) = dBdgth*UN(2) +  dBdGNy/Uth
     *    - UN(2)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth
          DBDUX(2,6) = dBdgth*UN(3) +  dBdGNz/Uth
     *    - UN(3)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth

!          write (*,*) DBDU(2,4), DBDU(2,5), DBDU(2,6)
!          write (*,*) DBDUX(2,4), DBDUX(2,5), DBDUX(2,6)

!         B(3) == JZ
          dBdnx =  CTM2*GN(2)+C66*FTN4
          dBdny = -CTM2*GN(1)+C66*FTN5
          dBdnz =  2.0D0*UN(3)*CTF+C266*Gth
          dBdth =  2.0D0*(CT*GN(3)+ST*DD45)+
     *             STF*(UN(3)**2-1.0D0)+C66*FTN7
           dBdGNx= -CTM2*UN(2)+C66*FTNX4
           dBdGNy=  CTM2*UN(1)+C66*FTNX5
           dBdGNz=  C266*ST
          dBdgth=  C266*UN(3)

          DBDU(3,4) = dBdth*UN(1) + dBdgth*GN(1)
     *  + dBdnx/Uth - UN(1)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(1)**2-1D0)-2D0*UN(1)*GN(1))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(1)*UN(2))-UN(1)*GN(2)-UN(2)*GN(1))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(1)*UN(3))-UN(1)*GN(3)-UN(3)*GN(1))/Uth
          DBDU(3,5) = dBdth*UN(2) + dBdgth*GN(2)
     *  + dBdny/Uth - UN(2)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(2)*UN(1))-UN(2)*GN(1)-UN(1)*GN(2))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(2)**2-1D0)-2D0*UN(2)*GN(2))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(2)*UN(3))-UN(2)*GN(3)-UN(3)*GN(2))/Uth
          DBDU(3,6) = dBdth*UN(3) + dBdgth*GN(3)
     *  + dBdnz/Uth - UN(3)*(dBdnx*UN(1)+dBdny*UN(2)+dBdnz*UN(3))/Uth
     *  +  dBdGNx*(Gth/Uth*(UN(3)*UN(1))-UN(3)*GN(1)-UN(1)*GN(3))/Uth
     *  +  dBdGNy*(Gth/Uth*(UN(3)*UN(2))-UN(3)*GN(2)-UN(2)*GN(3))/Uth
     *  +  dBdGNz*(Gth/Uth*(UN(3)**2-1D0)-2D0*UN(3)*GN(3))/Uth
          DBDUX(3,4) = dBdgth*UN(1) +  dBdGNx/Uth
     *    - UN(1)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth
          DBDUX(3,5) = dBdgth*UN(2) +  dBdGNy/Uth
     *    - UN(2)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth
          DBDUX(3,6) = dBdgth*UN(3) +  dBdGNz/Uth
     *    - UN(3)*( dBdGNx*UN(1)+ dBdGNy*UN(2)+ dBdGNz*UN(3))/Uth

!         Test spin currents
!          UJX=2.0D0*(Gth*UN(1)+ST*GN(1)+CTM*(UN(2)*GN(3)-GN(2)*UN(3)))
!     *     + (CTM*UN(1)*UN(3)+UN(2)*ST)*FTN + DA*UX(1)
!          UJY=2.0D0*(Gth*UN(2)+ST*GN(2)-CTM*(UN(1)*GN(3)-GN(1)*UN(3)))
!     *     + (CTM*UN(2)*UN(3)-UN(1)*ST)*FTN + DA*UX(2)
!          UJZ=2.0D0*(Gth*UN(3)+ST*GN(3)+CTM*(UN(1)*GN(2)-GN(1)*UN(2)))
!     *     + (CTM*UN(3)**2+CT)*FTN + DA*UX(3)
!          if (dabs(UJX).gt.1D-6) write(*,*) "UJX: ", UJX
!          if (dabs(UJY).gt.1D-6) write(*,*) "UJY: ", UJY
!          if (dabs(UJZ).gt.1D-6) write(*,*) "UJZ: ", UJZ

!          write (*,*) DBDU(3,4), DBDU(3,5), DBDU(3,6)
!          write (*,*) DBDUX(3,4), DBDUX(3,5), DBDUX(3,6)
!          write(*,*) dBdgth, UN(1), dBdgth, UN(2)
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
