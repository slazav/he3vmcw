C--------------- FFT ---------------------------------------
      PARAMETER (NNX=120000,ND=15000)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL FF(NNX),F1(ND),F2(ND),F3(ND),TT(NNX)
      CHARACTER*100 STR2
      CHARACTER*100 I_FILE
      CHARACTER*100 O_FILE
      PI2=8.0*DATAN(1.0D0)             ! 2*PI

      OPEN(89,FILE='FFT.PAR')
      READ(89,*)I_FILE
      READ(89,*)O_FILE
      READ(89,*)FC
      READ(89,*)DF
      READ(89,*)NF
      READ(89,*)T1
      READ(89,*)T2
      CLOSE(89)

      II=0
      open(88,FILE=I_FILE)
      read(88,*, END=11) STR2
      do NP=0,NNX
        if (NP.GE.NNX) then
          write (*,*)'Too many data points..'
          close(88)
          stop
        endif
        read(88,*, END=11) TTA,FFA
        if (TTA.GT.T1.AND.TTA.LE.T2) then
          TT(II)=TTA-T1
          FF(II)=FFA
          II=II+1
        endif
      enddo
   11 close(88)
      WRITE(*,*)'Number of points: ', NP+1, ' used: ', II+1

      DT=TT(2)-TT(1)
      if (2*NF+1.GT.ND) then
        write (*,*)'Too many points to calculate..'
        stop
      endif

      AA=0.0
      do I=0,II-1
        AA=AA+FF(I)
      enddo
      AA=AA/II
      do I=0,II-1
        FF(I)=FF(I)-AA
      enddo

      FI=FC-DF*NF
      do K=1,2*NF+1
        F1(K)=0.0
        F2(K)=0.0
      enddo
      FMAX=0
      AMAX=0
      open(87,FILE=O_FILE)
      do K=1,2*NF+1
        FCUR=FI+DF*(K-1)
        do I=0,II
          TI=DT*(I-1)
          F1(K)=F1(K)+FF(I)*DSIN(PI2*FCUR*TI)
          F2(K)=F2(K)+FF(I)*DCOS(PI2*FCUR*TI)
        enddo
        F11=F1(K)**2
        F22=F2(K)**2
        F3(K)=DSQRT(F11+F22)
        if (AMAX.LT.F3(K)) then
          AMAX=F3(K)
          FMAX=FCUR
        endif
        write(87,101)FCUR,F3(K)
      enddo
      write (*,'("MAX: "3(1PE15.6))') (T1+T2)/2, AMAX, FMAX
  101 format(4(1PE25.16))
      end
