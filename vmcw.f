C---------------- CB=0.0 !!!!!!!!!!
        implicit REAL*8(A-H,O-Z)
        include 'par.fh'
        common /IOUNIT/ LOUT
        common /TIMEP/ T
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),U(NPDE,NPTS),X(NPTS)
        common /CH_PAR/ CLE,SLP,SWR,EPS,MF,IC,MJW,IBN,BETA
        common /BLK_UMU/ T11,GW,W,W0,TOLD,AA,TF,AF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,XS,PI,DTW,FF,CCC,CK,DTW1
        dimension SCTCH(KORD*(NDERV+1)),WORK(IDIMWORK),IWORK(IDIMIWORK)
C--------------- INITIALIZATION -------------------------------------
        GAMMA=2.0378D+4    ! GAMMA
        open(54,FILE='cfg')
        read(54,*)BETA   ! TIPPING ANGLE (DEGREES)
        read(54,*)IBN    ! TYPE OF BOUND. COND.: 1-OPEN CELL 2-CLOSED CELL
        read(54,*)TEMP   ! TEMPERATURE (K)          (0.9299)
        read(54,*)H      ! FIELD (OE)               (110)
        read(54,*)GRAD   ! GRADIENT H (OE/CM)       (-0.2)
        read(54,*)CLE    ! CELL LENGTH (CM)         (0.4)
        read(54,*)HY     ! RF FIELD (OE)            (0.06)
        read(54,*)SLP    ! STARTING LARMOR POSITION (CM) (-0.1)
        read(54,*)SWR    ! SWEEP RATE (CM/SEC)      (0.03)
        read(54,*)DC     ! DIFFUSION * T^2          (1.38E-6)
        read(54,*)T1C    ! T1*T                     (1.0E-3)
        read(54,*)DTW1   ! STEP IN TIME FOR TIME DEPENDENCIES (1.0E-5)
        read(54,*)NEPS   ! RELATIVE TIME ERROR BOUND(1.0E-3)
        read(54,*)MF     ! METHOD OF PDE SOLUTION    (11)
                         ! 11 - ADAMS METHOD; USE DERIVF FOR JACOBIAN
                         ! 12 - ADAMS METHOD; CALCULATE JACOBIAN
                         ! 21 - BWD METHOD; USE DERIVF FOR JACOBIAN
                         ! 22 - BWD METHOD; CALCULATE JACOBIAN
        read(54,*)MJW    ! WRITE TO MJ FILE:  0-DON'T  1-WRITE
        read(54,*)TSW    ! TIME TO STOP SWEEP
        read(54,*)TW     ! TIME TO WRITE MJ
        read(54,*)TS     ! TIME OF X OSC
        read(54,*)XS     ! AMPL OF X-OSC
        read(54,*)FF     ! FF
        read(54,*)CCC    ! CCC (AERO KRUTIZNA  1000)
        read(54,*)CK     ! CK (KRUTIZNI KUCHNOST' 0.01)
        read(54,*)TOUT   ! TIME TO STOP COMPUTATION (SEC) (0.0001)
        close(54)

        FLP=SLP+SWR*TSW
        SH1=GRAD*FLP/H-XS/GAMMA/H
        SH2=GRAD*(FLP-CLE)/H-XS/GAMMA/H
        HMAL=HY/H
        open(76,FILE='shift')
        write(76,68)TSW,SH1,SH2,HMAL,HMAL/DSQRT(15.0D0)
        close(76)

        DTW=DTW1
        EPS=2.0D0**(-NEPS)
        open(44,FILE='vmcw.t')
        open(21,FILE='vmcw.mj')
        open(24,FILE='vmcwx.mj')
        open(47,FILE='pulse.t')

        open(10,FILE='nul')
        LOUT=10                     ! DIAGNOSTIC'S OUTPUT
        T=0.0                        ! STARTING TIME
        IC=0                       ! TIME DEPENDENCIES COUNTER
        call SET_MESH()
        call SET_ICOND()
        IWORK(1)=IDIMWORK
        IWORK(2)=IDIMIWORK

        open(49,FILE='aero.x')
        do J=1,NPTS
          CC=FER(X(J),CLE*0.5D0,CCC)
          write(49,*)X(J),CC
        enddo
        close(49)

        DELT=20.0D0/CCC
        TOLD=0.0D0
        DT=1.D-10                    ! INITIAL STEP SIZE IN T
        T0=T                         ! STARTING TIME
        INDEX=1                      ! TYPE OF CALL (FIRST CALL)
        do I=1,IDIMWORK
          WORK(I)=0.0D0
        enddo
C--------------- COMPUTE PARAMETERS ----------------------------
        W=GAMMA*H                    ! OMEGA (RAD/SEC)
        WY=GAMMA*HY                  ! rf-OMEGA (RAD/SEC)
        GW=GAMMA*GRAD                ! GRADIENT OMEGA (RAD/SEC/CM)
        T1=T1C*TEMP                  ! RELAXATION TIME
        T11=1.0D0/T1
        PI=4.0*DATAN(1.0D0)
        W0=GW*SLP
        TTC=TEMP/0.00093D0             ! T/TC  0 BAR
        TETC=DSQRT(1.0D0-TTC)
        FLEG=330460.0D0*TETC           ! LEGGETT FREQ.(HZ) 0 BAR
        CPAR=1300.0D0*TETC             ! SPIN WAVES VELOCITY  (2000)
        TF=0.0000005D0/TETC      !!!!!!!!!!         ! TAU EFFECTIVE (L-T) SECONDS
        AA=FLEG*FLEG/W*4.0D0*PI**2
        AF0=-CPAR**2/W
        DIFF=DC/TEMP/TEMP               ! DIFFUSION
        DW=GW*SWR
C----------------MAIN LOOP -------------------------------------------
   2    CONTINUE
        if(T.GE.TOUT) THEN
          write(*,*) 'FINAL TIME REACHED..'
          stop
        endif
        DTW2=DTW1*0.02D0
        if(T.GT.TS-DTW1)DTW=DTW2
        T=T+DTW
        call PDECOL(T0,T,DT,X,EPS,NINT,KORD,NCC,NPDE,MF,
     +   INDEX,WORK,IWORK)
        write(LOUT,*) 'TIME=',T,'INDEX=',INDEX
        if(INDEX.NE.0) THEN
          open(46,FILE='vmcw.err')
          write(46,*) 'INTEGRATION FAILED; INDEX=', INDEX
          close(46)
          close(44)
          close(21)
          close(24)
          close(47)
          close(48)
          stop
        endif
        call VALUES(X,USOL,SCTCH,NPDE,NPTS,NPTS,2,WORK)
        do I=1,NPTS
          do J=1,NPDE
            U(J,I)=USOL(J,I,1)
          enddo
        enddo
        call MONITOR()
        goto 2
   68   format(F8.4, 4(1PE18.9))
      end
C-- F ---------- EVALUATION OF F ------------------------------------
      subroutine F(T,X,U,UX,UXX,FV,NPDE)
        implicit REAL*8(A-H,O-Z)
        dimension U(NPDE),UX(NPDE),UXX(NPDE),FV(NPDE)
        common /CH_PAR/ CLE,SLP,SWR,EPS,MF,IC,MJW,IBN,BETA
        common /BLK_UMU/ T11,GW,W,W0,TOLD,AA,TF,AF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,XS,PI,DTW,FF,CCC,CK,DTW1
        if(T.GE.TSW)THEN
          WZ=W0+DW*TSW
        else
          WZ=W0+DW*T
        endif
        WR=W+X*GW
        WZR=WZ+W
        XZA=X*GW-WZ
        if(T.GE.TS)THEN
          XZ=XZA+XS
        else
          XZ=XZA
        endif
        UM=DSQRT(U(4)**2+U(5)**2+U(6)**2)
        U4=U(4)/UM
        U5=U(5)/UM
        U6=U(6)/UM
        WR2=0.5D0*WR
        WYT1=WY/WR-U(1)
        U31=U(3)-1.0D0
        USN=U(2)*U5+U31*U6-U4*WYT1
        DD45=U4*UX(5)-UX(4)*U5
        ST=DSIN(U(7))
        CT=DCOS(U(7))
        CTM=1.0D0-CT
        CT1=1.0D0+CT
        CTG=ST/CTM
        UT=ST*(1.0D0+4.0D0*CT)*0.2666666D0
        AUT0=AA*UT
        AF=AF0-AF0*0.5D0*FER(X,CLE*0.5D0,CCC)
        DAF=-AF0*0.5D0*DFER(X,CLE*0.5D0,CCC)
        AUT=AUT0-AUT0*0.835D0*FER(X,CLE*0.5D0,CCC)
        TF0=TF-TF*0.5D0*FER(X,CLE*0.5D0,CCC)
        FTN=CTM*DD45-ST*UX(6)-UX(7)*U6
        DFTN=CTM*(U4*UXX(5)-UXX(4)*U5)-ST*UXX(6)-UXX(7)*U6-
     *   CT1*UX(7)*UX(6)+ST*UX(7)*DD45

        UJX=2.0D0*(UX(7)*U(4)+ST*UX(4)+CTM*(U(5)*UX(6)-UX(5)*U(6)))+
     *   (CTM*U(4)*U(6)+U(5)*ST)*FTN
        UJY=2.0D0*(UX(7)*U(5)+ST*UX(5)-CTM*(U(4)*UX(6)-UX(4)*U(6)))+
     *   (CTM*U(5)*U(6)-U(4)*ST)*FTN
        UJZ=2.0D0*(UX(7)*U(6)+ST*UX(6)+CTM*(U(4)*UX(5)-UX(4)*U(5)))+
     *   (CTM*U(6)**2+CT)*FTN

        DJX=AF*(2.0D0*(UXX(7)*U4+CT1*UX(7)*UX(4)+ST*UXX(4)+
     *   ST*UX(7)*(U5*UX(6)-UX(5)*U6)+CTM*(U5*UXX(6)-UXX(5)*U6))+
     *   (CTM*U4*U6+U5*ST)*DFTN+(ST*UX(7)*U4*U6+
     *   CTM*(UX(4)*U6+U4*UX(6))+UX(5)*ST+U5*CT*UX(7))*FTN)-
     *   DIFF*UXX(1)+DAF*UJX
        DJY=AF*(2.0D0*(UXX(7)*U5+CT1*UX(7)*UX(5)+ST*UXX(5)-
     *   ST*UX(7)*(U4*UX(6)-UX(4)*U6)-CTM*(U4*UXX(6)-UXX(4)*U6))+
     *   (CTM*U5*U6-U4*ST)*DFTN+(ST*UX(7)*U5*U6+
     *   CTM*(UX(5)*U6+U5*UX(6))-UX(4)*ST-U4*CT*UX(7))*FTN)-   !!!!!!!!!
     *   DIFF*UXX(2)+DAF*UJY
        DJZ=AF*(2.0D0*(UXX(7)*U6+CT1*UX(7)*UX(6)+ST*UXX(6)+
     *   ST*UX(7)*DD45+CTM*(U4*UXX(5)-UXX(4)*U5))+
     *   (CTM*U6**2+CT)*DFTN+(ST*UX(7)*U6**2+
     *   CTM*2.0D0*U6*UX(6)-ST*UX(7))*FTN)-DIFF*UXX(3)+DAF*UJZ
        FV(1)=U(2)*XZ+AUT*U4-DJX
        FV(2)=-U(1)*XZ+AUT*U5+WY*U(3)-DJY
        FV(3)=AUT*U6-U31*T11-WY*U(2)-DJZ
        FV(4)=-WZR*U5-WR2*(U31*U5-U(2)*U6+CTG*
     *   ((U5*U(2)+U6*U31)*U4+WYT1*(U5**2+U6**2)))
        FV(5)=WZR*U4-WR2*(-WYT1*U6-U31*U4+CTG*
     *   ((U6*U31-U4*WYT1)*U5-U(2)*(U4**2+U6**2)))
        FV(6)=-WR2*(U(2)*U4+U5*WYT1+CTG*
     *   ((U5*U(2)-U4*WYT1)*U6-U31*(U4**2+U5**2)))
        FV(7)=WR*USN+UT/TF0                           !U(7)=TETA
        return
      end
C-- BNDRY ------ BOUNDARY CONDITIONS -- B(U,UX)=Z(T) ------------
      subroutine BNDRY(T,X,U,UX,DBDU,DBDUX,DZDT,NPDE)
        implicit REAL*8(A-H,O-Z)
        dimension U(NPDE),UX(NPDE),DZDT(NPDE),
     *   DBDU(NPDE,NPDE),DBDUX(NPDE,NPDE)
        common /CH_PAR/ CLE,SLP,SWR,EPS,MF,IC,MJW,IBN,BETA
        common /BLK_UMU/ T11,GW,W,W0,TOLD,AA,TF,AF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,XS,PI,DTW,FF,CCC,CK,DTW1
        do I=1,NPDE
          DZDT(I)=0.0D0
          do J=1,NPDE
            DBDU(I,J)=0.0D0
            DBDUX(I,J)=0.0D0
          enddo
        enddo
        UM=DSQRT(U(4)**2+U(5)**2+U(6)**2)
        U4=U(4)/UM
        U5=U(5)/UM
        U6=U(6)/UM
        ST=DSIN(U(7))
        ST2=2.0D0*ST
        CT=DCOS(U(7))
        CTM=1.0D0-CT
        CTM2=2.0D0*CTM
        DD45=U4*UX(5)-UX(4)*U5
        FTN=CTM*DD45-ST*UX(6)-UX(7)*U6
        CTF=CTM*FTN
        STF=ST*FTN
        FTN4=CTM*UX(5)
        FTN5=-CTM*UX(4)
        FTN7=ST*DD45-CT*UX(6)
        FTNX4=-CTM*U5
        FTNX5=CTM*U4
        C46=CTM*U4*U6+U5*ST
        C56=CTM*U5*U6-U4*ST             !!!!!!!!!!!
        C66=CTM*U6**2+CT
        C266=2.0D0-C66
        AF=AF0-AF0*0.5D0*FER(X,CLE*0.5D0,CCC)
        DA=-DIFF/AF
        if(IBN.EQ.2)THEN       ! CLOSED CELL
          DBDUX(4,1)=DA
          DBDUX(5,2)=DA
          DBDUX(6,3)=DA
          DBDU(4,4)=2.0D0*UX(7)+CTF*U6+C46*FTN4
          DBDU(4,5)=CTM2*UX(6)+STF+C46*FTN5
          DBDU(4,6)=-CTM2*UX(5)+CTF*U4-C46*UX(7)
          DBDU(4,7)=2.0D0*(CT*UX(4)+ST*(U5*UX(6)-UX(5)*U6))+
     *     STF*U4*U6+U5*CT*FTN+C46*FTN7
          DBDU(5,4)=-CTM2*UX(6)-STF+C56*FTN4
          DBDU(5,5)=2.0D0*UX(7)+CTF*U6+C56*FTN5
          DBDU(5,6)=CTM2*UX(4)+CTF*U5-C56*UX(7)
          DBDU(5,7)=2.0D0*(CT*UX(5)-ST*(U4*UX(6)-UX(4)*U6))+
     *     STF*U5*U6-U4*CT*FTN+C56*FTN7
          DBDU(6,4)=CTM2*UX(5)+C66*FTN4
          DBDU(6,5)=-CTM2*UX(4)+C66*FTN5
          DBDU(6,6)=2.0D0*U6*CTF+C266*UX(7)
          DBDU(6,7)=2.0D0*(CT*UX(6)+ST*DD45)+STF*(U6**2-1.0D0)+C66*FTN7
          DBDUX(4,4)=ST2+C46*FTNX4
          DBDUX(4,5)=-CTM2*U6+C46*FTNX5
          DBDUX(4,6)=CTM2*U5-C46*ST
          DBDUX(4,7)=2.0D0*U4-C46*U6
          DBDUX(5,4)=CTM2*U6+C56*FTNX4
          DBDUX(5,5)=ST2+C56*FTNX5
          DBDUX(5,6)=-CTM2*U4-C56*ST
          DBDUX(5,7)=2.0D0*U5-C56*U6
          DBDUX(6,4)=-CTM2*U5+C66*FTNX4
          DBDUX(6,5)=CTM2*U4+C66*FTNX5
          DBDUX(6,6)=C266*ST
          DBDUX(6,7)=C266*U6
C          DBDU(7,4)=UX(4)         !!
C          DBDU(7,5)=UX(5)         !!
C          DBDU(7,6)=UX(6)         !!
C          DBDUX(7,4)=U4         !!
C          DBDUX(7,5)=U5         !!
C          DBDUX(7,6)=U6         !!
        endif
        return
      end
C-- SET_MESH --- SET UP THE MESH ------------------------------------
      subroutine SET_MESH()
        implicit REAL*8(A-H,O-Z)
        include 'par.fh'
        dimension DXX(NPTS),XX(NPTS)
        common /SIG/ XM(NPTS)
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),U(NPDE,NPTS),X(NPTS)
        common /CH_PAR/ CLE,SLP,SWR,EPS,MF,IC,MJW,IBN,BETA
        common /BLK_UMU/ T11,GW,W,W0,TOLD,AA,TF,AF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,XS,PI,DTW,FF,CCC,CK,DTW1
        open(69,FILE='ress')
        open(57,FILE='mesh')
        CLE2=CLE*0.5D0
        DX=CLE/(NPTS-1)
        do J=1,NPTS
          XX(J)=DX*(J-1)
          CC=FER(XX(J),CLE2,CCC)
          CCD=ABS(DFER(XX(J),CLE2,CCC))
          write(69,*)XX(J),CC,CCD
        enddo
        close(69)
        KK=1
        EP=1.0D-20
   44   CONTINUE
        X(1)=0.0D0
        XT=0.0D0
        DXX(1)=DX/(1.0D0+CK*ABS(DFER(DX*0.5D0,CLE2,CCC)))
        do J=2,(NPTS-1)/2+1
          XT=XT+DXX(J-1)
          DXX(J)=DX/(1.0D0+CK*ABS(DFER(XT+DX*0.5D0,CLE2,CCC)))
        enddo
        KK=KK+1
        if(DABS(CLE2-XT).GT.EP.AND.KK.LT.2000)THEN
          DK=XT-CLE2
          DD=DK/(NPTS-1)*2
          write(*,*)'DK=',DK,'  K=',KK-1
          DX=DX-DD
          goto 44
        else
          do J=1,(NPTS-1)/2
            X(J+1)=X(J)+DXX(J)
          enddo
          LL=(NPTS-1)/2+1

          UDXX=1.0D0-(X(LL)-CLE2)/CLE2
C         UFF=(X(LL)-CLE2)/(NPTS-1)*2.0D0
          do J=1,(NPTS-1)/2
            DXX(J)=DXX(J)*UDXX
C           X(J+1)=X(J)+DXX(J)-UFF
            X(J+1)=X(J)+DXX(J)
C           write(57,55)J,X(J),DXX(J),UDXX
          enddo
        endif
        do J=NPTS,(NPTS-1)/2+1,-1
          X(J)=CLE-X(NPTS-J+1)
        enddo
        X(1)=0.0D0
        X(NPTS)=CLE
        do J=1,NPTS
          write(57,55)J,X(J)
        enddo
        do J=2,NPTS-1
          XM(J)=(1.0D5*X(J+1)-1.0D5*X(J-1))*0.5D0
        enddo
        XM(1)=(1.0D5*X(2)-1.0D5*X(1))*0.5D0
        XM(NPTS)=(1.0D5*X(NPTS)-1.0D5*X(NPTS-1))*0.5D0
        do J=1,NPTS
cc        write(57,55)J,X(J),XM(J)
        enddo
        close(57)

   55   format(I5, 7(1PE25.16))
      end
C-- SET_ICOND -- INITIAL CONDITIONS ---------------------------------
      subroutine SET_ICOND()
        implicit REAL*8(A-H,O-Z)
        include 'par.fh'
        common /CH_PAR/ CLE,SLP,SWR,EPS,MF,IC,MJW,IBN,BETA
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),U(NPDE,NPTS),X(NPTS)
        PI=4.0D0*DATAN(1.0D0)
        BET=BETA*PI/180.0D0
        UMZ=DCOS(BET)
        UMX=DSIN(BET)*DSQRT(0.5D0)
        UMY=UMX
        if(UMZ.GE.-0.25D0)THEN
          UNZ2=0.8D0*(0.25D0+DCOS(BET))
          UNZ=DSQRT(UNZ2)
          DELTA=(25.0D0*UNZ2+15.0D0)/16.0D0
          DELTAX=DSIN(BET)*DSQRT(0.5D0)*
     *     (UNZ*1.25D0-DSQRT(15.0D0)*0.25D0)
          DELTAY=DSIN(BET)*DSQRT(0.5D0)*
     *     (UNZ*1.25D0+DSQRT(15.0D0)*0.25D0)
          UNX=DELTAX/DELTA
          UNY=DELTAY/DELTA
          UCTG=DACOS(-0.25D0)
        else
          UNZ=0.0D0
          UNX=-DSQRT(0.5D0)
          UNY=DSQRT(0.5D0)
          UCTG=BET
        endif
        do I=1,NPTS
          USOL(1,I,1)=UMX             ! Mx
          USOL(2,I,1)=UMY             ! My
          USOL(3,I,1)=UMZ             ! Mz
          USOL(4,I,1)=UNX             ! Nx   !!
          USOL(5,I,1)=UNY             ! Ny   !!
          USOL(6,I,1)=UNZ             ! Nz
          USOL(7,I,1)=UCTG            ! TETA         !!!!!+/-
        enddo
        do I=1,NPTS
          do J=1,NPDE
            do K=2,3
              USOL(J,I,K)=0.0
            enddo
          enddo
        enddo
        do I=1,NPTS
          do J=1,NPDE
            U(J,I)=USOL(J,I,1)
          enddo
        enddo
        return
      end
C-- USP(X) ----- CSI OF SOLUTION ------------------------------------
      double precision function USP(XI,I)
        implicit REAL*8(A-H,O-Z)
        include 'par.fh'
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),U(NPDE,NPTS),X(NPTS)
        common /CH_PAR/ CLE,SLP,SWR,EPS,MF,IC,MJW,IBN,BETA
        do K=1,NPTS
          USM=DSQRT(USOL(5,K,1)**2+USOL(6,K,1)**2+USOL(4,K,1)**2)
          USOL(4,K,1)=USOL(4,K,1)/USM
          USOL(5,K,1)=USOL(5,K,1)/USM
          USOL(6,K,1)=USOL(6,K,1)/USM
        enddo
C      KL=1
C      KH=NPTS
C   19    IF(KH-KL.GT.1)THEN
C        K=(KH+KL)/2
C          IF(X(K).GT.XI)THEN
C          KH=K
C          ELSE
C          KL=K
C          ENDIF
C        GOTO 19
C        ENDIF
C      H=X(KH)-X(KL)
C      A=(X(KH)-XI)/H
C      B=(XI-X(KL))/H
C      YAL=USOL(I,KL,1)
C      YAH=USOL(I,KH,1)
C      Y2L=USOL(I,KL,3)
C      Y2H=USOL(I,KH,3)
C      USP=A*YAL+B*YAH+((A**3-A)*Y2L+(B**3-B)*Y2H)*(H**2)/6.0D0
C
CC      J=INT(XI/CLE*NINT+1.1)
CC      USP=USOL(I,J,1)
        AA=1.0D20
        IK=1
        do II=1,NPTS
          BB=(X(II)-XI)**2
          if(BB.LE.AA)THEN
            AA=BB
            IK=II
          endif
        enddo
        J=IK
        USP=USOL(I,J,1)
        return
      end
C-- UINIT ------ INITIAL CONDITIONS ---------------------------------
      subroutine UINIT(XI,UI,NPDEI)
        implicit REAL*8(A-H,O-Z)
        dimension UI(NPDEI)
        do I=1,NPDEI
          UI(I)=USP(XI,I)
        enddo
        return
      end
C-- DERIVF ----- SET UP DERIVATIVES ---------------------------------
      subroutine DERIVF(T,X,U,UX,UXX,DFDU,DFDUX,DFDUXX,NPDE)
        implicit REAL*8(A-H,O-Z)
        dimension U(NPDE),UX(NPDE),UXX(NPDE),
     *       DFDU(NPDE,NPDE),DFDUX(NPDE,NPDE),DFDUXX(NPDE,NPDE)
        common /CH_PAR/ CLE,SLP,SWR,EPS,MF,IC,MJW,IBN,BETA
        common /BLK_UMU/ T11,GW,W,W0,TOLD,AA,TF,AF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,XS,PI,DTW,FF,CCC,CK,DTW1
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
C-- FERMI_ST --- FERMI STEP -----------------------------------------
      double precision function FER(X1,X0,C)
        implicit REAL*8(A-H,O-Z)
        common /CH_PAR/ CLE,SLP,SWR,EPS,MF,IC,MJW,IBN,BETA
        common /BLK_UMU/ T11,GW,W,W0,TOLD,AA,TF,AF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,XS,PI,DTW,FF,CCC,CK,DTW1
        CLE2=CLE*0.5D0
        ARG=C*(ABS(X1-X0)-CLE2*FF)
        if(ARG.LT.82.0D0)THEN
          FER=1.0D0/(1.0D0+DEXP(ARG))
        else
          FER=0.0D0
        endif
        return
      end
C----- DFERMI STEP -----------------------------------------
      double precision function DFER(X1,X0,C)
        implicit REAL*8(A-H,O-Z)
        common /CH_PAR/ CLE,SLP,SWR,EPS,MF,IC,MJW,IBN,BETA
        common /BLK_UMU/ T11,GW,W,W0,TOLD,AA,TF,AF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,XS,PI,DTW,FF,CCC,CK,DTW1
        CLE2=CLE*0.5D0
        ARG=C*(ABS(X1-X0)-CLE2*FF)
        UU=DSIGN(1.0D0,X0-X1)
        if(ARG.LT.82.0D0)THEN
          DFER=C/((1.0D0+DEXP(ARG))**2)*DEXP(ARG)*UU
        else
          DFER=0.0D0
        endif
        return
      end
C-- MONITOR ---- MONITORING THE SOLUTION ----------------------------
      subroutine MONITOR()
        implicit REAL*8(A-H,O-Z)
        include 'par.fh'
        common /SIG/ XM(NPTS)
        common /TIMEP/ T
        common /GEAR0/ DTUSED,NQUSED,NSTEP,NFE,NJE
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),U(NPDE,NPTS),X(NPTS)
        common /TM_ARR/ TMMS(ITP),TMLP(ITP),TMAB(ITP),TMDS(ITP),
     *   TMPC(ITP),TMZ(ITP)
        common /CH_PAR/ CLE,SLP,SWR,EPS,MF,IC,MJW,IBN,BETA
        common /BLK_UMU/ T11,GW,W,W0,TOLD,AA,TF,AF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,XS,PI,DTW,FF,CCC,CK,DTW1
        IC = IC+1
C--------------- COMPUTE TIME DEPENDENCIES --------------------------
        TMMS(IC)=T*1000.0D0
        if(T.GE.TSW)THEN
          TMLP(IC)=SLP+SWR*TSW
        else
          TMLP(IC)=SLP+SWR*T
        endif
        TMAB(IC)=0.0D0
        TMDS(IC)=0.0D0
        TMZ(IC)=0.0D0
        do I=1,NPTS
          TMAB(IC)=TMAB(IC)+USP(X(I),1)*XM(I)
          TMDS(IC)=TMDS(IC)+USP(X(I),2)*XM(I)
          TMZ(IC)=TMZ(IC)+USP(X(I),3)*XM(I)
        enddo
        TMPC(IC)=(DSQRT(TMAB(IC)**2+TMDS(IC)**2))
        if(T.LE.TS)THEN
          write(44,61)
     +     TMMS(IC),TMLP(IC),TMAB(IC),TMDS(IC),TMPC(IC),TMZ(IC)
        endif
        if(T.GE.TS-DTW1)THEN
C          write(47,63)TMMS(IC),TMAB(IC),TMDS(IC),TMZ(IC)
          write(47,69)TMMS(IC),TMDS(IC)
        endif
C--------------- SHOW INFORMATION -----------------------------------
        write(*,12) ' TIME=',T*1000.,' (MSEC)  COUNT=',IC
        if(MJW.EQ.1.OR.T.GE.TW) CALL WRITE_MJ()    !
c        write(*,*) '  X    ','   MX    ','  MY     ','   MZ   ',
c     *'   NX    ','  NY   ','   NZ   '
c        DO I=1,NPTS,(NPTS-1)/21
c          write(*,72)X(I),U(1,I),U(2,I),U(3,I),U(4,I),U(5,I),U(6,I)
c        ENDDO
c        write(*,*)
        if(IC.EQ.1) RETURN
  12    format(A,F10.3,A,I5)
  72    format(F6.3, 6(1PE11.3))
   61   format(F7.1, 6(1PE14.6))
   63   format(F8.3, 3(1PE25.16))
   69   format(F8.3, 1PE25.16)
  181   format(A)
      end
C-- WRITE_MJ --- WRITE SPINS & CURRENTS TO VMCW ------------------
      subroutine WRITE_MJ()
        implicit REAL*8(A-H,O-Z)
        include 'par.fh'
        common /BLK_UMU/ T11,GW,W,W0,TOLD,AA,TF,AF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,XS,PI,DTW,FF,CCC,CK,DTW1
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),U(NPDE,NPTS),X(NPTS)
        do I=1, NPTS
          CT=DCOS(USOL(7,I,1))
          ST=DSIN(USOL(7,I,1))
          CTM=1.0D0-CT
          DD45=USOL(4,I,1)*USOL(5,I,2)-USOL(4,I,2)*USOL(5,I,1)
          FTN=CTM*DD45-ST*USOL(6,I,2)-USOL(7,I,2)*USOL(6,I,1)
          U4=USOL(4,I,1)
          U5=USOL(5,I,1)
          U6=USOL(6,I,1)
          UX1=USOL(1,I,2)
          UX2=USOL(2,I,2)
          UX3=USOL(3,I,2)
          UX4=USOL(4,I,2)
          UX5=USOL(5,I,2)
          UX6=USOL(6,I,2)
          UX7=USOL(7,I,2)

          UJX=2.0D0*(UX7*U4+ST*UX4+CTM*(U5*UX6-UX5*U6))+
     *     (CTM*U4*U6+U5*ST)*FTN
          UJY=2.0D0*(UX7*U5+ST*UX5-CTM*(U4*UX6-UX4*U6))+
     *     (CTM*U5*U6-U4*ST)*FTN
          UJZ=2.0D0*(UX7*U6+ST*UX6+CTM*(U4*UX5-UX4*U5))+
     *     (CTM*U6**2+CT)*FTN
          UFX=UJX*AF-DIFF*UX1
          UFY=UJY*AF-DIFF*UX2
          UFZ=UJZ*AF-DIFF*UX3

C          write(21,101)X(I),U(1,I),U(2,I),U(3,I),
C     *      U(4,I),U(5,I),U(6,I),U(7,I)
          write(24,102)X(I),USOL(1,I,2),USOL(2,I,2),USOL(3,I,2),
     *      USOL(4,I,2),USOL(5,I,2),USOL(6,I,2),USOL(7,I,2),UFX,UFY,UFZ
          write(21,101)X(I),USOL(1,I,1),USOL(2,I,1),USOL(3,I,1),
     *      USOL(4,I,1),USOL(5,I,1),USOL(6,I,1),USOL(7,I,1)

        enddo
  100   format(A)
  101   format(F10.6, 7(1PE25.16))
  102   format(F10.6, 10(1PE25.16))
      end
