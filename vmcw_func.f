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
        real*8 dTh
        real*8 ST,CT,CTG,CTGT
        real*8 UT,AUT,AF,DAF,B
        real*8 UJ(3),DJ(3) ! spin current J and dJ/dz
        integer i

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

C       For NPDE==6 it is a value of theta, for NPDE==7 it is
C       just a correction of n length (and theta will come from U(7))
        Uth=dsqrt(U(4)**2+U(5)**2+U(6)**2)
        UN(1) = U(4)/Uth
        UN(2) = U(5)/Uth
        UN(3) = U(6)/Uth

        if (NPDE.EQ.7) then
          Uth = U(7)
          GN(1) = UX(4)
          GN(2) = UX(5)
          GN(3) = UX(6)
          Gth   = UX(7)
          GGN(1) = UXX(4)
          GGN(2) = UXX(5)
          GGN(3) = UXX(6)
          GGth   = UXX(7)

C         GN and GGN can be fixed by the following way:
C           gn = gn - n*(n*gn);
C           ggn = ggn - n*(ggn*n + gn*gn)
C         (it makes them correct and do not change correct values).
C         NOTE: this works only with BC=1!
!          if (1.eq.1) then
!            UDU = UN(1)*GN(1) + UN(2)*GN(2) + UN(3)*GN(3)
!            DO i=1,3
!              GN(i) = GN(i) - UN(i)*UDU
!            ENDDO
!
!            UDU = UN(1)*GGN(1) + UN(2)*GGN(2) + UN(3)*GGN(3)
!     *          + GN(1)*GN(1) + GN(2)*GN(2) + GN(3)*GN(3)
!            DO i=1,3
!               GGN(i) = GGN(i) - UN(i)*UDU
!            ENDDO
!          ENDIF

        else
C         TODO: move to grad.c
C         U(4,5,6) = n*th
C         th = |U|
C         n = U/|U|
C         dTh = d|U| = (U*dU)/|U|
C         dn  = dU/|U| - U d|U|/|U|^2
          Gth = (U(4)*UX(4)+U(5)*UX(5)+U(6)*UX(6))/Uth
          GN(1) = UX(4)/Uth - U(4)*Gth/Uth**2
          GN(2) = UX(5)/Uth - U(5)*Gth/Uth**2
          GN(3) = UX(6)/Uth - U(6)*Gth/Uth**2

C         d2|U| = (dU*dU)/|U| + (U*d2U)/|U| - (U*dU)^2/|U|^3
          GGth = (UX(4)**2+UX(5)**2+U(6)**2
     *           +U(4)*UXX(4)+U(5)*UXX(5)+U(6)*UXX(6)
     *         - Gth**2)/Uth

C         d2n = d2U/|U| - 2 dU d|U| / |U|^2
C             - U d2|U|/|U|^2 + 2 U d|U|^2/|U|^3
          GGN(1) = UXX(4)/Uth - 2D0*UX(4)*Gth/Uth**2
     *           - U(4)*GGth/Uth**2 + 2D0*U(4)*Gth**2/Uth**3
          GGN(2) = UXX(5)/Uth - 2D0*UX(5)*Gth/Uth**2
     *           - U(5)*GGth/Uth**2 + 2D0*U(5)*Gth**2/Uth**3
          GGN(3) = UXX(6)/Uth - 2D0*UX(6)*Gth/Uth**2
     *           - U(6)*GGth/Uth**2 + 2D0*U(6)*Gth**2/Uth**3
        endif


        ST=dsin(Uth)
        CT=dcos(Uth)
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
        DJ(1) = AF*DJ(1) + DAF*UJ(1) - DIFF*GGMx
        DJ(2) = AF*DJ(2) + DAF*UJ(2) - DIFF*GGMy
        DJ(3) = AF*DJ(3) + DAF*UJ(3) - DIFF*GGMz

        B = UN(1)*UMxm + UMym*UN(2) + UMzm*UN(3)

C       Leggett equations
        FV(1)=   (Wz-W0)*UMy          + AUT*UN(1) - DJ(1)
        FV(2)= - (Wz-W0)*UMx + Wr*UMz + AUT*UN(2) - DJ(2)
        FV(3)=               - Wr*UMy + AUT*UN(3) - DJ(3) - UMzm/T1

        if (NPDE.EQ.7) then
C         if cos(th) is near 1 print a warning
          if (dabs(1D0-CT).lt.1D-6) write (*,*) "Warning: CT: ", CT
          CTG = ST/(1D0-CT)
          FV(4) = -W0*UN(2)
     *      -0.5D0*Wz*(UMzm*UN(2)-UMym*UN(3)+CTG*(B*UN(1)-UMxm))
          FV(5) =  W0*UN(1)
     *      -0.5D0*Wz*(UMxm*UN(3)-UMzm*UN(1)+CTG*(B*UN(2)-UMym))
          FV(6) =
     *      -0.5D0*Wz*(UMym*UN(1)-UMxm*UN(2)+CTG*(B*UN(3)-UMzm))
          FV(7) =  Wz*B + UT/Tf

!          UDU = UN(1)*FV(4) + UN(2)*FV(5) + UN(3)*FV(6)
!          DO i=1,3
!            FV(i+3) = FV(i+3) - UN(i)*UDU
!          ENDDO
        else

C         avoid divergency in cos(th)=1
          if (dabs(1D0-CT).gt.1D-3) then
            CTGT=ST/(1.0D0-CT)*Uth      ! T*ctg(T/2) = sin(T)/(1-cos(T))
          else
            CTGT=2D0-Uth**2/6D0
          endif

C         d(n*th) = th*dn + n*dth
          dTh = Wz*B + UT/Tf
          FV(4) = -W0*UN(2)*Uth - 0.5D0*Wz*
     *       ((UMzm*UN(2)-UMym*UN(3))*Uth + CTGT*(B*UN(1)-UMxm))
     *       + dTh*UN(1)
          FV(5) = W0*UN(1)*Uth - 0.5D0*Wz*
     *       ((UMxm*UN(3)-UMzm*UN(1))*Uth + CTGT*(B*UN(2)-UMym))
     *       + dTh*UN(2)
          FV(6) = -0.5D0*Wz*
     *       ((UMym*UN(1)-UMxm*UN(2))*Uth + CTGT*(B*UN(3)-UMzm))
     *        + dTh*UN(3)
        endif

        return
      end


C-- BNDRY ------ BOUNDARY CONDITIONS -- B(U,UX)=Z(T) ------------

      subroutine BNDRY(T,X,U,UX,DBDU,DBDUX,DZDT,NPDE)
        implicit none
        integer NPDE
        dimension U(NPDE),UX(NPDE),DZDT(NPDE),
     *   DBDU(NPDE,NPDE),DBDUX(NPDE,NPDE)
        real*8 T,X,U,UX,DZDT,DBDU,DBDUX
        real*8 UN(3), Uth
        real*8 GN(3), Gth, UDU
        real*8 DJDU(4,3), DJDUZ(4,3)
        real*8 DJDUa(3,3), DJDUb(3,3), DJDUZa(3,3), DJDUZb(3,3)
        real*8 DA
!      real*8 U4,U5,U6,STF,ST2,ST,FTNX5,FTNX4,FTN7,FTN5,FTN4,FTN
!      real*8 DD45,CTF,CTM,CTM2,C266,C46,C56,C66,CT
        integer I,J

!       paramaters, same as in F (only Cpar, Diff used)
        real*8 Cpar, Diff, Wz
        integer IBN

!       set parameters:
        call set_bndry_pars(T,X, Wz, Cpar, Diff, IBN)

!       set everything to zero:
        do I=1,NPDE
          DZDT(I)=0.0D0
          do J=1,NPDE
            DBDU(I,J)=0.0D0
            DBDUX(I,J)=0.0D0
          enddo
        enddo

!       BC1: Zero derivatives
        if (IBN.EQ.1) THEN
          do J=1,NPDE
            DBDUX(J,J)=1.0D0
          enddo
          return
        endif

!       BC3: Zero derivatives of M, fixed n and theta
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

!         see similar formulas in FUNC
          Uth=dsqrt(U(4)**2+U(5)**2+U(6)**2)
          UN(1) = U(4)/Uth
          UN(2) = U(5)/Uth
          UN(3) = U(6)/Uth
          if (NPDE.EQ.7) then
            Uth = U(7)
            GN(1) = UX(4)
            GN(2) = UX(5)
            GN(3) = UX(6)
            Gth   = UX(7)

!!!         this compensation does not work!
!              UDU = UN(1)*GN(1) + UN(2)*GN(2) + UN(3)*GN(3)
!              DO i=1,3
!                GN(i) = GN(i) - UN(i)*UDU
!              ENDDO
          else
            Gth = (U(4)*UX(4)+U(5)*UX(5)+U(6)*UX(6))/Uth
            GN(1) = UX(4)/Uth - U(4)*Gth/Uth**2
            GN(2) = UX(5)/Uth - U(5)*Gth/Uth**2
            GN(3) = UX(6)/Uth - U(6)*Gth/Uth**2
          endif

          DA=Wz*Diff/Cpar**2

          if (NPDE.EQ.7) then
            call fill_djd_nt(DJDU, DJDUZ, UN,Uth,GN,Gth)
            do I=1,3
              DBDUX(I+3,I) = DA
              do J=1,4
                DBDU(I+3,J+3)  = -DJDU(J,I)
                DBDUX(I+3,J+3) = -DJDUZ(J,I)
              enddo
            enddo
          else
            call fill_dj_t(DJDUa, DJDUb, DJDUZa, DJDUZb, UN,Uth,GN,Gth)
            do I=1,3
              DBDUX(I+0,I) = DA
              do J=1,3
                DBDU(I+0,J+3)  = -(DJDUa(J,I)/2D0 + DJDUb(J,I))
                DBDUX(I+0,J+3) = -(DJDUZa(J,I)/2D0 + DJDUZb(J,I))
              enddo
            enddo
          endif
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
