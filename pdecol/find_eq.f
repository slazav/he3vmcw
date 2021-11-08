!====================================================================
C      AN EXTENSION TO PDECOL ALGORITHM (540R) TO FIND EQUILIBRIUM:
C      F( T, X, U, UX, UXX ) = 0
C      V.ZAVJALOV, 2021

      SUBROUTINE FINDEQ (T, WORK, IWORK, MSG_LVL)
        IMPLICIT NONE

C       Use pointers in a comon block to pass
C       T, WORK, IWORK into FINDEQ_MINFUNC
        REAL*8,TARGET :: WORK(:), T
        REAL*8, POINTER :: WORKP(:), TP
        INTEGER,TARGET :: IWORK(:)
        INTEGER,POINTER :: IWORKP(:)
        COMMON /FINDEQ_CMN/ TP, WORKP, IWORKP

C       Parameters of PDECOL solver
        COMMON /SIZES/ NINT,KORD,NCC,NPDE,NCPTS,NEQN,IQUAD
        INTEGER NINT,KORD,NCC,NPDE,NCPTS,NEQN,IQUAD

        REAL*8  C(NEQN),F,G(NEQN)
        INTEGER ERROR,MSG_LVL
        REAL*8 WORK_TN(14*NEQN)

        EXTERNAL FINDEQ_MINFUNC

        TP => T
        WORKP => WORK
        IWORKP => IWORK

        CALL FINDEQ_GETC(C, WORK)

        CALL TN(ERROR, NEQN, C, F, G, WORK_TN, 14*NEQN,
     *          FINDEQ_MINFUNC, MSG_LVL)

        CALL FINDEQ_SETC(C, WORK)
      END

!====================================================================
C  GET VALUES OF C FOR PDECOL SOLVER
      SUBROUTINE FINDEQ_GETC(C,WORK)
        IMPLICIT NONE
        REAL*8 C(NEQN)

C       Parameters of PDECOL solver
        COMMON /SIZES/ NINT,KORD,NCC,NPDE,NCPTS,NEQN,IQUAD
        INTEGER NINT,KORD,NCC,NPDE,NCPTS,NEQN,IQUAD
        COMMON /OPTION/ NOGAUS,MAXDER
        INTEGER NOGAUS,MAXDER

C       Work arrays and indices
        REAL*8  WORK(KORD+NPDE*(4+9*NPDE)+(KORD+(NINT-1)*(KORD-NCC))*
     *             (3*KORD+2+NPDE*(3*(KORD-1)*NPDE+MAXDER+4)))
        COMMON /ISTART/ IW1,IW2,IW3,IW4,IW5,IW6,IW7,IW8,IW9,IW10,
     *                  IW11,IW12,IW13,IW14,IW15,IW16,IW17,IW18
        INTEGER IW1,IW2,IW3,IW4,IW5,IW6,IW7,IW8,IW9,IW10,
     *          IW11,IW12,IW13,IW14,IW15,IW16,IW17,IW18
        INTEGER J

        DO J=1,NEQN
          C(J) = WORK(IW10+J-1)
        ENDDO
      END

!====================================================================
C  SET VALUES OF C FOR PDECOL SOLVER, SET C TIME DERIVATIVES TO 0
      SUBROUTINE FINDEQ_SETC(C,WORK)
        IMPLICIT NONE
        REAL*8 C(NEQN)

C       Parameters of PDECOL solver
        COMMON /SIZES/ NINT,KORD,NCC,NPDE,NCPTS,NEQN,IQUAD
        INTEGER NINT,KORD,NCC,NPDE,NCPTS,NEQN,IQUAD
        COMMON /OPTION/ NOGAUS,MAXDER
        INTEGER NOGAUS,MAXDER

C       Work arrays and indices
        REAL*8  WORK(KORD+NPDE*(4+9*NPDE)+(KORD+(NINT-1)*(KORD-NCC))*
     *             (3*KORD+2+NPDE*(3*(KORD-1)*NPDE+MAXDER+4)))
        COMMON /ISTART/ IW1,IW2,IW3,IW4,IW5,IW6,IW7,IW8,IW9,IW10,
     *                  IW11,IW12,IW13,IW14,IW15,IW16,IW17,IW18
        INTEGER IW1,IW2,IW3,IW4,IW5,IW6,IW7,IW8,IW9,IW10,
     *          IW11,IW12,IW13,IW14,IW15,IW16,IW17,IW18
        INTEGER I,J

        DO J=1,NEQN
          WORK(IW10+J-1) = C(J)
        ENDDO

        DO I=1,MAXDER
          DO J=1,NEQN
            WORK(IW10 + NEQN*I + J - 1) = 0D0
          ENDDO
        ENDDO

      END

!====================================================================
C  CALCULATE F=sum(GFUN(Ci)^2) and Gi = dF/dCi.
C  FOR PDECOL SOLVER
      SUBROUTINE FINDEQ_MINFUNC(C,F,G)
        IMPLICIT NONE
        REAL*8 C(NEQN), F, G(NEQN)

C       Parameters of PDECOL solver
        COMMON /SIZES/ NINT,KORD,NCC,NPDE,NCPTS,NEQN,IQUAD
        INTEGER NINT,KORD,NCC,NPDE,NCPTS,NEQN,IQUAD
        COMMON /OPTION/ NOGAUS,MAXDER
        INTEGER NOGAUS,MAXDER

        COMMON /GEAR9/ EPSJ,R0,ML,MU,MW,NM1,N0ML,N0W
        REAL*8 EPSJ,R0
        INTEGER ML,MU,MW,NM1,N0ML,N0W

C       Work arrays and indices
        INTEGER,POINTER :: IWORK(:)
        REAL*8,POINTER :: WORK(:)
        COMMON /ISTART/ IW1,IW2,IW3,IW4,IW5,IW6,IW7,IW8,IW9,IW10,
     *                  IW11,IW12,IW13,IW14,IW15,IW16,IW17,IW18
        INTEGER IW1,IW2,IW3,IW4,IW5,IW6,IW7,IW8,IW9,IW10,
     *          IW11,IW12,IW13,IW14,IW15,IW16,IW17,IW18
        REAL*8, POINTER :: T
        COMMON /FINDEQ_CMN/ T, WORK, IWORK

        REAL*8 R,RINV,CJ
        REAL*8 FF1(NPDE), FF2(NPDE)
        INTEGER N,I,J

        CALL GFUN (T, C, FF1, NPDE, NCPTS, WORK(IW1), WORK,
     *             WORK(IW14), WORK(IW15), WORK(IW16),
     *             WORK(IW3), WORK(IW9), IWORK)
        F = 0D0
        DO I=1,NPDE
          F = F + FF1(I)*FF1(I)
        ENDDO
        DO J=1,NEQN
          R = EPSJ * WORK(IW4+J-1) ! Step, same as in PDECOD/DIFFF
          R = DMAX1(R,R0)
          CJ = C(J)
          C(J) = C(J) + R
          RINV = 1D0 / R
          CALL GFUN (T, C, FF2, NPDE, NCPTS, WORK(IW1), WORK,
     *               WORK(IW14), WORK(IW15), WORK(IW16),
     *               WORK(IW3), WORK(IW9), IWORK)
          G(J) = 0D0
          DO I=1,NPDE
            G(J) = G(J) + 2D0*FF1(I)*(FF2(I)-FF1(I)) * RINV
          ENDDO
          C(J) = CJ
        ENDDO
        RETURN
      END
