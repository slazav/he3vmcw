C===============================================
C Finding equilibrium:
C Use values on the user-providid mesh, skipping 2 points
C from each side (which should be kept in accordance with
C boundary conditions) as free parameters.
C Minimize value of SUM(F()^2)
C Parameters:
C  T             - time (for F function)
C  X(NPTS)       - X breakpoints
C  U(NPDE, NPTS) - function values (input, output)
C  NPTS          - size of X and U arrays
C  MSG_LVL       - TN message level

      SUBROUTINE FINDEQ (T, X, UVAL, NPTS, MSG_LVL, ERROR)
        IMPLICIT NONE

C       Parameters of PDECOL solver (we need NPDE)
        COMMON /SIZES/ NINT,KORD,NCC,NPDE,NCPTS,NEQN,IQUAD
        INTEGER NINT,KORD,NCC,NPDE,NCPTS,NEQN,IQUAD
        INTEGER NPTS

C       Use pointers in a comon block to pass
C       parameters to min function
        REAL*8,TARGET :: UVAL(NPDE,NPTS), X(NPTS), T
        REAL*8, POINTER :: UVALP(:,:), XP(:), TP
        COMMON /FINDEQ_CMN/ TP, XP, UVALP, NPDE1, NPTS1
        INTEGER NPDE1, NPTS1

        INTEGER ERROR,MSG_LVL,N,I,J,IX
        REAL*8 WORK_TN(NPDE*(NPTS-4)*14)

        REAL*8 FF, GG(NPDE*(NPTS-4)), UU(NPDE*(NPTS-4))

        DOUBLE PRECISION ETA, ACCRCY, XTOL, STEPMX, DSQRT, MCHPR1
        INTEGER MAXIT, MAXFUN
        EXTERNAL FINDEQ_MINFUNC

        TP => T
        XP => X
        UVALP => UVAL
        NPDE1 = NPDE
        NPTS1 = NPTS


        if (NPTS.LT.5) then
          ERROR=101
          return
        endif
        N=NPDE*(NPTS-4)

C       Transfer values from UVAL to UU array
        do I=3,NPTS-2
          do J=1,NPDE
            IX = J + NPDE*(I-3)
            UU(IX) = UVAL(J,I)
            UU(IX) = 1D-1
          enddo
        enddo

        FF = 0D0

C SET UP PARAMETERS FOR THE OPTIMIZATION ROUTINE
C ETA    - SEVERITY OF THE LINESEARCH
C MAXFUN - MAXIMUM ALLOWABLE NUMBER OF FUNCTION EVALUATIONS
C XTOL   - DESIRED ACCURACY FOR THE SOLUTION X*
C STEPMX - MAXIMUM ALLOWABLE STEP IN THE LINESEARCH
C ACCRCY - ACCURACY OF COMPUTED FUNCTION VALUES
C MSGLVL - DETERMINES QUANTITY OF PRINTED OUTPUT
C          0 = NONE, 1 = ONE LINE PER MAJOR ITERATION.
C MAXIT  - MAXIMUM NUMBER OF INNER ITERATIONS PER STEP

        MAXIT = 10*N/2
        IF (MAXIT .GT. 500) MAXIT = 10*500
        IF (MAXIT .LE. 0) MAXIT = 1
        MAXFUN = 10*150*N
        ETA = .25D0
        STEPMX = 1.D-2
        ACCRCY = 1.D2*MCHPR1()
        XTOL = DSQRT(ACCRCY)


!        CALL FINDEQ_MINFUNC_TEST(N,UU)

        CALL LMQN (ERROR, N, UU, FF, GG, WORK_TN, 14*N,
     *     FINDEQ_MINFUNC, MSG_LVL,
     *     MAXIT, MAXFUN, ETA, STEPMX, ACCRCY, XTOL)

C       Transfer values from UU to UVAL
        do I=3,NPTS-2
          do J=1,NPDE
            IX = J + NPDE*(I-3)
            UVAL(J,I) = UU(IX)
          enddo
!          write (*,*) I, UVAL(1,I), GG(1+NPDE*(I-3))
        enddo
      END

!====================================================================
C  Test FINDEQ_MINFUNC: modify a values and cmpare gradients
      SUBROUTINE FINDEQ_MINFUNC_TEST(N,UU)
        INTEGER N,I
        REAL*8 UU(N),FF1,GG1(N), FF2,GG2(N)
        REAL*8 U1, R, GF

        CALL FINDEQ_MINFUNC(N,UU,FF1,GG1)
        R = 1D-5
        do I=1,N
          U1 = UU(I)
          UU(I) = UU(I) + R
          CALL FINDEQ_MINFUNC(N,UU,FF2,GG2)
          UU(I) = U1
          GF = (FF2-FF1)/R
          write(*,*) I, GF, GG1(I)

        enddo
      end


!====================================================================
C  Function for minimization (SUM(UDOT)) and its derivatives
C  I = 2..NPTS-1
C  RETURN SUM(F^2)
      SUBROUTINE FINDEQ_MINFUNC(N,UU,FF,GG)
        IMPLICIT NONE
        REAL*8 RET
        INTEGER I,J,K,IX
        INTEGER N
        REAL*8 UU(N), FF, GG(N)


        REAL*8, POINTER :: UVAL(:,:), X(:), T
        COMMON /FINDEQ_CMN/ T, X, UVAL, NPDE, NPTS
        INTEGER NPDE, NPTS

        REAL*8 R,U0
        REAL*8 F1, F2
        REAL*8 FINDEQ_MINFUNC_F

C       Check that N = (NPTS-4)*NPDE
        if (N.ne.(NPTS-4)*NPDE) then
           write(*,*) 'FINDEQ_MINFUNC E0'
           return
        endif

C       Transfer values from UU to UVAL array
        do I=3,NPTS-2
          do J=1,NPDE
            IX = J + NPDE*(I-3)
            UVAL(J,I) = UU(IX)
          enddo
        enddo

        FF = FINDEQ_MINFUNC_F(0)

        R = 1D-4
        do I=3,NPTS-2
          do J=1,NPDE
            U0 = UVAL(J,I)
            UVAL(J,I) = U0 - R
            F1 = FINDEQ_MINFUNC_F(I)
            UVAL(J,I) = U0 + R
            F2 = FINDEQ_MINFUNC_F(I)
            UVAL(J,I) = U0

            IX = J + NPDE*(I-3)
            GG(IX) = (F2-F1)/2D0/R
          enddo ! J=1,NPDE
        enddo ! I=3,NPTS-2

!        write(*,*) '> ', FF, UU(1), GG(1)

        RETURN
      END



!====================================================================
C  Function for minimization (SUM(UDOT^2)) on the given grid
C  without gradients.
C  IF K=0 then calculate the whole sum, if K = 3..NPTS-2,
C  Calculate K-1..K+1 part
C
      FUNCTION FINDEQ_MINFUNC_F(K)
        IMPLICIT NONE
        INTEGER I,J,K, IMIN, IMAX
        REAL*8 FINDEQ_MINFUNC_F

        REAL*8, POINTER :: UVAL(:,:), X(:), T
        COMMON /FINDEQ_CMN/ T, X, UVAL, NPDE, NPTS
        INTEGER NPDE, NPTS

        REAL*8 U(NPDE),DU(NPDE),DDU(NPDE)
        REAL*8 A,B,DET
        REAL*8 FI(NPDE)

        if (K.eq.0) then
          IMIN=3
          IMAX=NPTS-2
        else
          IMIN = MAX(3,K-1)
          IMAX = MIN(NPTS-2,K+1)
        endif

        FINDEQ_MINFUNC_F = 0D0
        do I=IMIN, IMAX
          do J=1,NPDE

C         We have function U(i) on the grid X(i) (i=1..NPTS)
C         Near point i we approximate U with a quadratic function
C         U = A*x^2 + B*x + C
C
C         To calculate function derivatives in x_i we solve system
C              A*X(i+1)^2 + B*X(i+1) + C = U(i+1)
C              A*X(i)^2   + B*X(i)   + C = U(i)
C              A*X(i-1)^2 + B*X(i-1) + C = U(i-1)
C
C         Then derivatives dU/dX and d2U/dxdx at i:
C           DU(i) = 2*A*X(I) + B,  DDU(i) = 2*A
C
            DET = X(I+1)**2*(X(  I)-X(I-1))
     *          + X(  I)**2*(X(I-1)-X(I+1))
     *          + X(I-1)**2*(X(I+1)-X(  I))

            A = UVAL(J,I+1)*(X(  I)-X(I-1))
     *        + UVAL(J,  I)*(X(I-1)-X(I+1))
     *        + UVAL(J,I-1)*(X(I+1)-X(  I))
            B = X(I+1)**2*(UVAL(J,  I)-UVAL(J,I-1))
     *        + X(  I)**2*(UVAL(J,I-1)-UVAL(J,I+1))
     *        + X(I-1)**2*(UVAL(J,I+1)-UVAL(J,  I))
            A = A/DET
            B = B/DET
            U(J)   = UVAL(J,I)
            DU(J)  = 2D0*A*X(I) + B
            DDU(J) = 2D0*A

          enddo ! J=1,NPDE

          CALL F(T,X(I),U,DU,DDU,FI,NPDE)

          do J=1,NPDE
            FINDEQ_MINFUNC_F = FINDEQ_MINFUNC_F + FI(J)**2
          enddo ! J=1,NPDE

        enddo ! I=IMIN,IMAX
        FINDEQ_MINFUNC_F = FINDEQ_MINFUNC_F/(NPTS-4)/NPDE * 1D-12

        RETURN
      END

