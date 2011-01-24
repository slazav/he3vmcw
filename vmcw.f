C---------------- CB=0.0 !!!!!!!!!!
        implicit REAL*8(A-H,O-Z)
        include 'par.fh'
        common /TIMEP/ T
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /CH_PAR/ SLP,SWR,BETA,MJW,IBN

        common /CFG_AER/  AER, AER_LEN, AER_CNT, AER_TRW
        common /CFG_CELL/ CELL_LEN
        common /CFG_MESH/ XMESH_K,XMESH_ACC


        common /BLK_UMU/ T11,GW,W,W0,AA,TF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,TIME_STEP,XS,PI
        dimension SCTCH(KORD*(NDERV+1)),WORK(IDIMWORK),IWORK(IDIMIWORK)
        character*64 CFG_KEY

C--------------- INITIALIZATION -------------------------------------
        GAMMA=2.0378D+4    ! GAMMA

        open(54,FILE='vmcw.cfg')
   11   read(54,*,END=12) CFG_KEY, CFG_VAL
        if (CFG_KEY.EQ.'BETA') then
          BETA=CFG_VAL ! INITIAL TIPPING ANGLE (DEGREES)
        elseif (CFG_KEY.EQ.'SLP') then
          SLP=CFG_VAL  ! STARTING LARMOR POSITION (CM) (-0.1)
        elseif (CFG_KEY.EQ.'SWR') then
          SWR=CFG_VAL  ! SWEEP RATE (CM/SEC)      (0.03)

        elseif (CFG_KEY.EQ.'DC') then
          DC=CFG_VAL   ! DIFFUSION * T^2          (1.38E-6)
        elseif (CFG_KEY.EQ.'T1C') then
          T1C=CFG_VAL  ! T1*T                     (1.0E-3)
        elseif (CFG_KEY.EQ.'PDECOL_ACC_LOG2') then
          PDECOL_ACC_LOG2=CFG_VAL ! RELATIVE TIME ERROR BOUND

        elseif (CFG_KEY.EQ.'TEMP') then
          TEMP=CFG_VAL ! TEMPERATURE (K)          (0.9299)
        elseif (CFG_KEY.EQ.'H') then
          H=CFG_VAL    ! FIELD (OE)               (110)
        elseif (CFG_KEY.EQ.'GRAD') then
          GRAD=CFG_VAL ! GRADIENT H (OE/CM)       (-0.2)
        elseif (CFG_KEY.EQ.'HY') then
          HY=CFG_VAL   ! RF FIELD (OE)            (0.06)

        elseif (CFG_KEY.EQ.'IBN') then
          IBN=CFG_VAL  ! TYPE OF BOUND. COND.: 1-OPEN CELL 2-CLOSED CELL
        elseif (CFG_KEY.EQ.'MJW') then
          MJW=CFG_VAL  ! WRITE TO MJ FILE:  0-DON'T  1-WRITE
        elseif (CFG_KEY.EQ.'TSW') then
          TSW=CFG_VAL  ! TIME TO STOP SWEEP
        elseif (CFG_KEY.EQ.'TW') then
          TW=CFG_VAL   ! TIME TO WRITE MJ
        elseif (CFG_KEY.EQ.'TS') then
          TS=CFG_VAL   ! TIME OF X OSC
        elseif (CFG_KEY.EQ.'XS') then
          XS=CFG_VAL   ! AMPL OF X-OSC
        elseif (CFG_KEY.EQ.'TIME') then
          TIME=CFG_VAL ! TIME TO STOP COMPUTATION (SEC) (0.0001)
        elseif (CFG_KEY.EQ.'TIME_STEP') then
          TIME_STEP=CFG_VAL ! STEP IN TIME FOR TIME DEPENDENCIES

C       CFG_CELL parameter group:
        elseif (CFG_KEY.EQ.'CELL_LEN') then
          CELL_LEN=CFG_VAL
C       CFG_MESH parameter group:
        elseif (CFG_KEY.EQ.'XMESH_K') then
          XMESH_K=CFG_VAL
        elseif (CFG_KEY.EQ.'XMESH_ACC') then
          XMESH_ACC=CFG_VAL
C       CFG_AER parameter group:
        elseif (CFG_KEY.EQ.'AER') then
          AER=CFG_VAL
        elseif (CFG_KEY.EQ.'AER_LEN') then
          AER_LEN=CFG_VAL
        elseif (CFG_KEY.EQ.'AER_CNT') then
          AER_CNT=CFG_VAL
        elseif (CFG_KEY.EQ.'AER_TRW') then
          AER_TRW=CFG_VAL

        else
          write(*,*) 'warning: unknown parameter in cfg-file: ', CFG_KEY
        endif
        goto 11
   12   close(54)

        FLP=SLP+SWR*TSW
        SH1=GRAD*FLP/H-XS/GAMMA/H
        SH2=GRAD*(FLP-CELL_LEN)/H-XS/GAMMA/H
        HMAL=HY/H
        open(76,FILE='shift')
        write(76,68)TSW,SH1,SH2,HMAL,HMAL/DSQRT(15.0D0)
        close(76)

        DTW=TIME_STEP
        open(44,FILE='vmcw.t')
        open(21,FILE='vmcw.mj')
        open(24,FILE='vmcwx.mj')
        open(47,FILE='pulse.t')


C       PDECOL parameters:
        INDEX=1                      ! TYPE OF CALL (FIRST CALL)
        MF=22
        IWORK(1)=IDIMWORK
        IWORK(2)=IDIMIWORK
        do I=1,IDIMWORK
          WORK(I)=0.0D0
        enddo
        PDECOL_ACC = 2.0D0**(-PDECOL_ACC_LOG2)
        T=0.0                        ! TIME
        T0=T                         ! STARTING TIME FOR PDECOL
        DT=1.D-10                    ! INITIAL STEP SIZE IN T

        call SET_MESH()
        call SET_ICOND()

        open(54,FILE='aer_step.dat')
        write(54,*), '#  I    X(I) STEP(X) STEP''(X)'
        do J=1,NPTS
          write(54,'(I4," ",F7.5," ",F7.5," ",e12.5e2)')
     +     J, X(J), AER_STEP(X(J),0), AER_STEP(X(J),1)
        enddo
        close(54)

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
        TF=0.0000005D0/TETC            ! TAU EFFECTIVE (L-T) SECONDS
        AA=FLEG*FLEG/W*4.0D0*PI**2
        AF0=-CPAR**2/W
        DIFF=DC/TEMP/TEMP               ! DIFFUSION
        DW=GW*SWR
C----------------MAIN LOOP -------------------------------------------
   2    CONTINUE

          if(T.GE.TIME) THEN
            write(*,*) 'FINAL TIME REACHED..'
            stop
          endif

          if(T.GT.TS-TIME_STEP) DTW = TIME_STEP*0.02D0

          T=T+DTW

          call PDECOL(T0,T,DT,X,PDECOL_ACC,NINT,KORD,NCC,NPDE,MF,
     +                INDEX,WORK,IWORK)
          if(INDEX.NE.0) THEN
            write(*,*) 'INTEGRATION FAILED; INDEX=', INDEX
            stop
          endif

          call VALUES(X,USOL,SCTCH,NPDE,NPTS,NPTS,2,WORK)
          call MONITOR()
        goto 2
   68   format(F8.4, 4(1PE18.9))
      end
C-- F ---------- EVALUATION OF F ------------------------------------
      subroutine F(T,X,U,UX,UXX,FV,NPDE)
        implicit REAL*8(A-H,O-Z)
        dimension U(NPDE),UX(NPDE),UXX(NPDE),FV(NPDE)
        common /CH_PAR/ SLP,SWR,BETA,MJW,IBN
        common /BLK_UMU/ T11,GW,W,W0,AA,TF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,TIME_STEP,XS,PI
C	T - time
C	X - x-coord
C	U   - Mx My Mz Nx Ny Nz T 
C	UX  - dU/dx
C	UXX - d2U/dx2
C	FV  - result

C	W = GAMMA H
C	WY - RF-field

C       calculate freq
        if(T.GE.TSW)THEN
          WZ=W0+DW*TSW
        else
          WZ=W0+DW*T
        endif
        WR=W+X*GW
        WZR=WZ+W
        XZ=X*GW-WZ

C	x-field step an TS
        if(T.GE.TS)THEN
          XZ=XZ+XS
        endif

C	fix n vector length
        UN=DSQRT(U(4)**2+U(5)**2+U(6)**2)
        UNx = U(4)/UN
        UNy = U(5)/UN
        UNz = U(6)/UN

        UMxm = U(1)-WY/WR
        UMym = U(2)
        UMzm = U(3)-1.0D0

        WR2 = 0.5D0*WR
        DD45=UNx*UX(5)-UX(4)*UNy       ! Nx Ny` - Nx` Ny
        ST=DSIN(U(7))
        CT=DCOS(U(7))
        CTM=1.0D0-CT
        CT1=1.0D0+CT
        CTG=ST/CTM      ! ctg(T/2) = sin(T)/(1-cos(T))
        UT=ST*(1.0D0+4.0D0*CT)*0.2666666D0
        AUT0=AA*UT

        AF=AF0 - AF0*0.5D0 * AER_STEP(X,0)
        DAF=-AF0*0.5D0 * AER_STEP(X,1)
        AUT=AUT0 - AUT0*0.835D0 * AER_STEP(X,0)
    	TF0=TF-TF*0.5D0 * AER_STEP(X,0)

        FTN=CTM*DD45-ST*UX(6)-UX(7)*UNz
        DFTN=CTM*(UNx*UXX(5)-UXX(4)*UNy)-ST*UXX(6)-UXX(7)*UNz-
     *   CT1*UX(7)*UX(6)+ST*UX(7)*DD45

        UJX=2.0D0*(UX(7)*UNx+ST*UX(4)+CTM*(UNy*UX(6)-UX(5)*UNz))+
     *   (CTM*UNx*UNz+UNy*ST)*FTN
        UJY=2.0D0*(UX(7)*UNy+ST*UX(5)-CTM*(UNx*UX(6)-UX(4)*UNz))+
     *   (CTM*UNy*UNz-UNx*ST)*FTN
        UJZ=2.0D0*(UX(7)*UNz+ST*UX(6)+CTM*(UNx*UX(5)-UX(4)*UNy))+
     *   (CTM*UNz**2+CT)*FTN

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

        FV(1)=   XZ*U(2)           + AUT*UNx - DJX
        FV(2)= - XZ*U(1) + WY*U(3) + AUT*UNy - DJY
        FV(3)=           - WY*U(2) + AUT*UNz - DJZ - UMzm*T11

        FV(4)=-WZR*UNy - WR2*(UMzm*UNy-UMym*UNz+CTG*(B*UNx-UMxm))
        FV(5)= WZR*UNx - WR2*(UMxm*UNz-UMzm*UNx+CTG*(B*UNy-UMym))
        FV(6)=         - WR2*(UMym*UNx-UMxm*UNy+CTG*(B*UNz-UMzm))
        FV(7)=WR*B + UT/TF0
        return
      end

C-- BNDRY ------ BOUNDARY CONDITIONS -- B(U,UX)=Z(T) ------------
      subroutine BNDRY(T,X,U,UX,DBDU,DBDUX,DZDT,NPDE)
        implicit REAL*8(A-H,O-Z)
        dimension U(NPDE),UX(NPDE),DZDT(NPDE),
     *   DBDU(NPDE,NPDE),DBDUX(NPDE,NPDE)
        common /CH_PAR/ SLP,SWR,BETA,MJW,IBN
        common /BLK_UMU/ T11,GW,W,W0,AA,TF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,TIME_STEP,XS,PI
        do I=1,NPDE
          DZDT(I)=0.0D0
          do J=1,NPDE
            DBDU(I,J)=0.0D0
            DBDUX(I,J)=0.0D0
          enddo
        enddo

C	fix n vector length
        UN=DSQRT(U(4)**2+U(5)**2+U(6)**2)
        UNx=U(4)/UN
        UNy=U(5)/UN
        UNz=U(6)/UN

        ST=DSIN(U(7))
        ST2=2.0D0*ST
        CT=DCOS(U(7))
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
        AF=AF0-AF0*0.5D0 * AER_STEP(X,0)
        DA=-DIFF/AF
        if(IBN.EQ.2)THEN       ! CLOSED CELL
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
C-- SET_ICOND -- INITIAL CONDITIONS ---------------------------------
      subroutine SET_ICOND()
        implicit REAL*8(A-H,O-Z)
        include 'par.fh'
        common /CH_PAR/ SLP,SWR,BETA,MJW,IBN
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
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
        return
      end
C-- USP(X) ----- CSI OF SOLUTION ------------------------------------
      double precision function USP(XI,I)
        implicit REAL*8(A-H,O-Z)
        include 'par.fh'
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /CH_PAR/ SLP,SWR,BETA,MJW,IBN
        do K=1,NPTS
          USM=DSQRT(USOL(5,K,1)**2+USOL(6,K,1)**2+USOL(4,K,1)**2)
          USOL(4,K,1)=USOL(4,K,1)/USM
          USOL(5,K,1)=USOL(5,K,1)/USM
          USOL(6,K,1)=USOL(6,K,1)/USM
        enddo
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
        common /CH_PAR/ SLP,SWR,BETA,MJW,IBN
        common /BLK_UMU/ T11,GW,W,W0,AA,TF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,TIME_STEP,XS,PI
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

C------- AEROGEL STEP -----------------------------------------
C       Aerogel density function. Returns 1 in the central
C       part of the cell with fermi steps to 0 on edges.
C         X        -- coord, cm
C         D        -- derivative order (0|1)
C         CELL_LEN -- cell length, cm
C         AER      -- if >0 then do step
C         AER_LEN  -- aerogel length / cell length
C         AER_CNT  -- center of aerogel area / cell length
C         AER_TRW  -- transition width / cell length
      double precision function AER_STEP(X,D)
        implicit REAL*8(A-H,O-Z)
        integer D
        common /CFG_AER/  AER, AER_LEN, AER_CNT, AER_TRW
        common /CFG_CELL/ CELL_LEN
        if (AER.LE.0) then
          AER_STEP=0.0D0
          return
        endif
        ARG=(ABS(X/CELL_LEN - AER_CNT) - AER_LEN*0.5D0) / AER_TRW
        if(ARG.LT.82.0D0)THEN
          if (D.eq.0) then
            AER_STEP=1.0D0/(1.0D0+DEXP(ARG))
          elseif (D.eq.1) then
            STEP_SIGN=DSIGN(1.0D0,X/CELL_LEN-AER_CNT)
            AER_STEP=STEP_SIGN*DEXP(ARG)/AER_TRW/CELL_LEN/
     +                (1.0D0+DEXP(ARG))**2
          else
            AER_STEP=0.0D0
          endif
        else
          AER_STEP=0.0D0
        endif
        return
      end
C-- SET_MESH --- SET UP THE MESH ------------------------------------
      subroutine SET_MESH()
C	Set mesh according with AER_STEP function
        implicit REAL*8(A-H,O-Z)
        include 'par.fh'
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /CFG_CELL/ CELL_LEN
        common /CFG_MESH/ XMESH_K,XMESH_ACC
        X(1)=0
C	start with homogenious mesh with DX intervals
        DX=CELL_LEN/(NPTS-1)
        do I=1,100
C         build mesh with scaled intervals
          do J=1,NPTS-1
            X(J+1)=X(J) +
     +        DX/(1.0D0+XMESH_K*ABS(AER_STEP(X(J),1)))
          enddo
C	  scale the whole mesh to fit CELL_LEN
          DELTA=CELL_LEN - X(NPTS)
          DX = DX + DELTA/(NPTS+1)
          if (ABS(DELTA).LT.XMESH_ACC) return
        enddo
        write(*,*) 'warning: low mesh accuracy: ', ABS(DELTA)
      end

C-- MONITOR ---- MONITORING THE SOLUTION ----------------------------
      subroutine MONITOR()
        implicit REAL*8(A-H,O-Z)
        include 'par.fh'
        common /TIMEP/ T
        common /GEAR0/ DTUSED,NQUSED,NSTEP,NFE,NJE
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /CH_PAR/ SLP,SWR,BETA,MJW,IBN
        common /BLK_UMU/ T11,GW,W,W0,AA,TF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,TIME_STEP,XS,PI
C--------------- COMPUTE TIME DEPENDENCIES --------------------------
        TMMS=T*1000.0D0
        if(T.GE.TSW)THEN
          TMLP=SLP+SWR*TSW
        else
          TMLP=SLP+SWR*T
        endif
        TMAB=0.0D0
        TMDS=0.0D0
        TMZ=0.0D0
        do I=1,NPTS-1
          TMAB=TMAB + USOL(1,I,1) * (X(I+1)-X(I)) *1E5
          TMDS=TMDS + USOL(2,I,1) * (X(I+1)-X(I)) *1E5
          TMZ=TMZ   + USOL(3,I,1) * (X(I+1)-X(I)) *1E5
        enddo
        TMPC=(DSQRT(TMAB**2+TMDS**2))
        if(T.LE.TS)THEN
          write(44,61)
     +     TMMS,TMLP,TMAB,TMDS,TMPC,TMZ
        endif
        if(T.GE.TS-TIME_STEP) THEN
          write(47,69)TMMS,TMDS
        endif
C--------------- SHOW INFORMATION -----------------------------------
        write(*,'(A,F6.1,A)') ' TIME=',T*1000.,' ms'
        if(MJW.EQ.1.OR.T.GE.TW) CALL WRITE_MJ()    !
   61   format(F7.1, 6(1PE14.6))
   69   format(F8.3, 1PE25.16)
      end
C-- WRITE_MJ --- WRITE SPINS & CURRENTS TO VMCW ------------------
      subroutine WRITE_MJ()
        implicit REAL*8(A-H,O-Z)
        include 'par.fh'
        common /BLK_UMU/ T11,GW,W,W0,AA,TF,DIFF,WY,DW,TSW,TW,
     +   AF0,TS,TIME_STEP,XS,PI
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /TIMEP/ T
        do I=1, NPTS, 128
          CT=DCOS(USOL(7,I,1))
          ST=DSIN(USOL(7,I,1))
          CTM=1.0D0-CT
          UNx=USOL(4,I,1)
          UNy=USOL(5,I,1)
          UNz=USOL(6,I,1)

          UX1=USOL(1,I,2)
          UX2=USOL(2,I,2)
          UX3=USOL(3,I,2)
          UX4=USOL(4,I,2)
          UX5=USOL(5,I,2)
          UX6=USOL(6,I,2)
          UX7=USOL(7,I,2)

          DD45 = UNx*UX5 - UX4*UNy
          FTN = CTM*DD45 - ST*UX6 - UX7*UNz

          UJX=2.0D0*(UX7*UNx+ST*UX4+CTM*(UNy*UX6-UX5*UNz))+
     *     (CTM*UNx*UNz+UNy*ST)*FTN
          UJY=2.0D0*(UX7*UNy+ST*UX5-CTM*(UNx*UX6-UX4*UNz))+
     *     (CTM*UNy*UNz-UNx*ST)*FTN
          UJZ=2.0D0*(UX7*UNz+ST*UX6+CTM*(UNx*UX5-UX4*UNy))+
     *     (CTM*UNz**2+CT)*FTN
          AF=AF0-AF0*0.5D0 * AER_STEP(X,0)
          UFX=UJX*AF-DIFF*UX1
          UFY=UJY*AF-DIFF*UX2
          UFZ=UJZ*AF-DIFF*UX3

          write(21,101) T*1000.,X(I),
     *      USOL(1,I,1),USOL(2,I,1),USOL(3,I,1),
     *      USOL(4,I,1),USOL(5,I,1),USOL(6,I,1),USOL(7,I,1)
          write(24,102) T*1000., X(I),
     *      USOL(1,I,2),USOL(2,I,2),USOL(3,I,2),
     *      USOL(4,I,2),USOL(5,I,2),USOL(6,I,2),USOL(7,I,2),
     *      UFX,UFY,UFZ

        enddo
        write(24,*)''
        write(21,*)''
  101   format(F7.1 F10.6, 7(1PE15.6))
  102   format(F7.1 F10.6, 10(1PE15.6))
      end
