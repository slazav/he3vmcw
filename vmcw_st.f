
C---------------- CB=0.0 !!!!!!!!!!
        implicit real*8(A-H,O-Z)
        include 'par.fh'
        include 'he3_const.fh'
        common /TIMEP/ T, TEND
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /CH_PAR/ SLP,SWR,BETA,MJW,IBN

        common /CFG_AER/  AER, AER_LEN, AER_CNT, AER_TRW
        common /CFG_CELL/ CELL_LEN
        common /CFG_MESH/ XMESH_K,XMESH_ACC

        common /BLK_UMU/ T11,GW,W,AA,TF,DIFF,WY,TW,
     +   AF0,TS,TIME_STEP,XS
        dimension SCTCH(KORD*(NDERV+1)),WORK(IDIMWORK),IWORK(IDIMIWORK)
        character*64 CFG_KEY

        integer FILES_MJ(NPTS) ! files for writing MJ along the cell
        integer FILE_MMJ       ! file for writing mean values
        common /FILES/ FILES_MJ, FILE_MMJ

        character CMD_FILE_NAME*20  ! file for reading commands
        integer   CMD_FILE          ! file descriptor
        data CMD_FILE_NAME/'vmcw.cmd'/, CMD_FILE/200/
        common /CMD/ CMD_FILE, CMD_FILE_NAME

        real*8 WRITEMJ_XSTEP
        common /CFG_WRITE/ WRITEMJ_XSTEP

C--------------- INITIALIZATION -------------------------------------

        WRITEMJ_XSTEP=0.1D0

        open(54,FILE='vmcw.cfg')
   11   read(54,*,END=12) CFG_KEY, CFG_VAL
        if (CFG_KEY.EQ.'BETA') then
          BETA=CFG_VAL ! INITIAL TIPPING ANGLE (DEGREES)
        elseif (CFG_KEY.EQ.'SLP') then
          SLP=CFG_VAL  ! STARTING LARMOR POSITION (CM) (-0.1)

        elseif (CFG_KEY.EQ.'T1C') then
          T1C=CFG_VAL  ! T1*T                     (1.0E-3)
        elseif (CFG_KEY.EQ.'PDECOL_ACC_LOG2') then
          PDECOL_ACC_LOG2=CFG_VAL ! RELATIVE TIME ERROR BOUND

        elseif (CFG_KEY.EQ.'TEMP') then
          TEMP=CFG_VAL ! TEMPERATURE (mK)
        elseif (CFG_KEY.EQ.'PRESS') then
          PRESS=CFG_VAL ! PRESS (bar)

        elseif (CFG_KEY.EQ.'H') then
          H=CFG_VAL    ! FIELD (OE)               (110)
        elseif (CFG_KEY.EQ.'GRAD') then
          GRAD=CFG_VAL ! GRADIENT H (OE/CM)       (-0.2)
        elseif (CFG_KEY.EQ.'HY') then
          HY=CFG_VAL   ! RF FIELD (OE)            (0.06)

        elseif (CFG_KEY.EQ.'IBN') then
          IBN=INT(CFG_VAL)  ! TYPE OF BOUND. COND.: 1-OPEN CELL 2-CLOSED CELL
        elseif (CFG_KEY.EQ.'MJW') then
          MJW=INT(CFG_VAL)  ! WRITE TO MJ FILE:  0-DON'T  1-WRITE
        elseif (CFG_KEY.EQ.'TW') then
          TW=CFG_VAL   ! TIME TO WRITE MJ
        elseif (CFG_KEY.EQ.'TS') then
          TS=CFG_VAL   ! TIME OF X OSC
        elseif (CFG_KEY.EQ.'XS') then
          XS=CFG_VAL   ! AMPL OF X-OSC
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
C       CFG_WRITE parameter group:
        elseif (CFG_KEY.EQ.'WRITEMJ_XSTEP') then
          WRITEMJ_XSTEP=CFG_VAL
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


        call SET_MESH()
        TIME=0D0
        SWR=0D0
        call SET_ICOND()

        T0=T                         ! STARTING TIME FOR PDECOL
        DT=1.D-10                    ! INITIAL STEP SIZE IN T


        call WRITEMJ_OPEN()

C--------------- COMPUTE PARAMETERS ----------------------------
        PI=4.0D0*DATAN(1.0D0)
        W=GAM*H                    ! OMEGA (RAD/SEC)
        WY=GAM*HY                  ! rf-OMEGA (RAD/SEC)
        GW=GAM*GRAD                ! GRADIENT OMEGA (RAD/SEC/CM)
        T1=T1C*TEMP/1000D0         ! RELAXATION TIME
        T11=1.0D0/T1

        TTC=TEMP/TCF(PRESS)             ! T/TC  19.5 BAR
        TETC=DSQRT(1.0D0-TTC)
C        FLEG=330460.0D0*TETC           ! LEGGETT FREQ.(HZ) 0 BAR
        CPAR=CPARF(PRESS,TEMP)          ! SPIN WAVES VELOCITY
        FLEG=DSQRT(LF2F(PRESS,TTC))     ! LEGGETT FREQ
        DIFF=DF(PRESS,TEMP)             ! SPIN DIFFUSION

        TF=1.2D-7/TETC                  ! TAU EFFECTIVE (L-T) SECONDS WV pic.10.5 20bar

        AA=FLEG*FLEG/W*4.0D0*PI**2
        AF0=-CPAR**2/W

        write(*,'(" P: ",F5.2," bar (Tc = ",F5.3," mK, Tab = ",F5.3,
     *            " mK), T: ",F5.3," mK = ",F5.3," Tc")'),
     *    PRESS, TCF(PRESS), TABF(PRESS), TEMP, TTC
        write(*,'(" F_legg: ", F9.3, " kHz,  ",
     .            " D: ", E9.2, " cm^2/s, ",
     .            " T_lt: ", E8.2, " s, ",
     .            " C_par: ", F6.1, " cm/s ")'),
     .    FLEG/1D3, DIFF, TF, CPAR

        call CMD_OPEN()
C----------------MAIN LOOP -------------------------------------------
   2    CONTINUE

          if(T.GE.TEND) call CMD_READ()

          if(T.GT.TS-TIME_STEP) then
            T = T + TIME_STEP * 0.02D0
          else
            T = T + TIME_STEP
          endif
          call PDECOL(T0,T,DT,X,PDECOL_ACC,NINT,KORD,NCC,NPDE,MF,
     +                INDEX,WORK,IWORK)
          if(INDEX.NE.0) THEN
            write(*,*) 'INTEGRATION FAILED; INDEX=', INDEX
            stop
          endif

          call VALUES(X,USOL,SCTCH,NPDE,NPTS,NPTS,2,WORK)
          call MONITOR()
C          call STATE_DUMP()
        goto 2
   68   format(F8.4, 4(1PE18.9))
      end
C-- F ---------- EVALUATION OF F ------------------------------------
      subroutine F(T,X,U,UX,UXX,FV,NPDE)
        implicit real*8(A-H,O-Z)
        dimension U(NPDE),UX(NPDE),UXX(NPDE),FV(NPDE)
        common /CH_PAR/ SLP,SWR,BETA,MJW,IBN
        common /BLK_UMU/ T11,GW,W,AA,TF,DIFF,WY,TW,
     +   AF0,TS,TIME_STEP,XS
C       T - time
C       X - x-coord
C       U   - Mx My Mz Nx Ny Nz T 
C       UX  - dU/dx
C       UXX - d2U/dx2
C       FV  - result

C       W = GAM H
C       WY - RF-field

C       calculate freq
        WZ=GW*(SLP+SWR*T)
        WR=W+X*GW
        WZR=WZ+W
        XZ=X*GW-WZ

C       x-field step an TS
        if(T.GE.TS)THEN
          if (X.GE.0.09D0) XZ=XZ+XS
        endif

C       fix n vector length
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

        FV(4)= - WZR*UNy - WR2*(UMzm*UNy-UMym*UNz+CTG*(B*UNx-UMxm))
        FV(5)=   WZR*UNx - WR2*(UMxm*UNz-UMzm*UNx+CTG*(B*UNy-UMym))
        FV(6)=           - WR2*(UMym*UNx-UMxm*UNy+CTG*(B*UNz-UMzm))
        FV(7)= WR*B + UT/TF0
        return
      end

C-- BNDRY ------ BOUNDARY CONDITIONS -- B(U,UX)=Z(T) ------------
      subroutine BNDRY(T,X,U,UX,DBDU,DBDUX,DZDT,NPDE)
        implicit real*8(A-H,O-Z)
        dimension U(NPDE),UX(NPDE),DZDT(NPDE),
     *   DBDU(NPDE,NPDE),DBDUX(NPDE,NPDE)
        common /CH_PAR/ SLP,SWR,BETA,MJW,IBN
        common /BLK_UMU/ T11,GW,W,AA,TF,DIFF,WY,TW,
     +   AF0,TS,TIME_STEP,XS
        do I=1,NPDE
          DZDT(I)=0.0D0
          do J=1,NPDE
            DBDU(I,J)=0.0D0
            DBDUX(I,J)=0.0D0
          enddo
        enddo

        if(IBN.EQ.2)THEN       ! CLOSED CELL

C         fix n vector length
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
        implicit real*8(A-H,O-Z)
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
              USOL(J,I,K)=0D0
            enddo
          enddo
        enddo
        return
      end
C-- USP(X) ----- CSI OF SOLUTION ------------------------------------
      double precision function USP(XI,I)
        implicit real*8(A-H,O-Z)
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
        implicit real*8(A-H,O-Z)
        dimension UI(NPDEI)
        do I=1,NPDEI
          UI(I)=USP(XI,I)
        enddo
        return
      end
C-- DERIVF ----- SET UP DERIVATIVES ---------------------------------
      subroutine DERIVF(T,X,U,UX,UXX,DFDU,DFDUX,DFDUXX,NPDE)
        implicit real*8(A-H,O-Z)
        dimension U(NPDE),UX(NPDE),UXX(NPDE),
     *       DFDU(NPDE,NPDE),DFDUX(NPDE,NPDE),DFDUXX(NPDE,NPDE)
        common /CH_PAR/ SLP,SWR,BETA,MJW,IBN
        common /BLK_UMU/ T11,GW,W,AA,TF,DIFF,WY,TW,
     +   AF0,TS,TIME_STEP,XS
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
        implicit real*8(A-H,O-Z)
        integer D
        common /CFG_AER/  AER, AER_LEN, AER_CNT, AER_TRW
        common /CFG_CELL/ CELL_LEN
        if (AER.LE.0D0) then
          AER_STEP=0D0
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
C       Set mesh according with AER_STEP function
        implicit real*8(A-H,O-Z)
        include 'par.fh'
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /CFG_CELL/ CELL_LEN
        common /CFG_MESH/ XMESH_K,XMESH_ACC
        X(1)=0D0
C       start with homogenious mesh with DX intervals
        DX=CELL_LEN/DFLOAT(NPTS-1)
        do I=1,100
C         build mesh with scaled intervals
          do J=1,NPTS-1
            X(J+1)=X(J) +
     +        DX/(1.0D0+XMESH_K*ABS(AER_STEP(X(J),1)))
          enddo
C         scale the whole mesh to fit CELL_LEN
          DELTA=CELL_LEN - X(NPTS)
          DX = DX + DELTA/DFLOAT(NPTS+1)
          if (ABS(DELTA).LT.XMESH_ACC) return
        enddo
        write(*,*) 'warning: low mesh accuracy: ', ABS(DELTA)
      end

      subroutine SAVE_MESH(FNAME)
        include 'par.fh'
        common /CFG_AER/  AER, AER_LEN, AER_CNT, AER_TRW
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        character FNAME*128
        integer FILE_AER
        FILE_AER=54
        open(FILE_AER,FILE=FNAME)
        write(FILE_AER,*), '#  I    X(I) STEP(X) STEP''(X)'
        do J=1,NPTS
          write(FILE_AER,'(I4," ",F7.5," ",F7.5," ",e12.5e2)')
     +     J, X(J), AER_STEP(X(J),0), AER_STEP(X(J),1)
        enddo
        close(FILE_AER)
      end

C-- MONITOR ---- MONITORING THE SOLUTION ----------------------------
      subroutine MONITOR()
        implicit real*8(A-H,O-Z)
        include 'par.fh'
        common /TIMEP/ T, TEND
        common /GEAR0/ DTUSED,NQUSED,NSTEP,NFE,NJE
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /CH_PAR/ SLP,SWR,BETA,MJW,IBN
        common /BLK_UMU/ T11,GW,W,AA,TF,DIFF,WY,TW,
     +   AF0,TS,TIME_STEP,XS
C--------------- COMPUTE TIME DEPENDENCIES --------------------------
        TMMS=T*1000.0D0

        TMLP=GW*(SLP+SWR*T)
        TMAB=0.0D0
        TMDS=0.0D0
        TMZ=0.0D0
        do I=1,NPTS-1
          TMAB=TMAB + USOL(1,I,1) * (X(I+1)-X(I)) *1D5
          TMDS=TMDS + USOL(2,I,1) * (X(I+1)-X(I)) *1D5
          TMZ=TMZ   + USOL(3,I,1) * (X(I+1)-X(I)) *1D5
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
        write(*,'(A,F6.1,A)') ' TIME=',T*1000D0,' ms'
        if(MJW.EQ.1.OR.T.GE.TW) CALL WRITE_MJ()    !
        call WRITEMJ_DO()
   61   format(F7.1, 6(1PE14.6))
   69   format(F9.2, 1PE25.16)
      end
C-- WRITE_MJ --- WRITE SPINS & CURRENTS TO VMCW ------------------
      subroutine WRITE_MJ()
        implicit real*8(A-H,O-Z)
        include 'par.fh'
        common /BLK_UMU/ T11,GW,W,AA,TF,DIFF,WY,TW,
     +   AF0,TS,TIME_STEP,XS
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /TIMEP/ T, TEND
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

          write(21,101) T*1000D0,X(I),
     *      USOL(1,I,1),USOL(2,I,1),USOL(3,I,1),
     *      USOL(4,I,1),USOL(5,I,1),USOL(6,I,1),USOL(7,I,1)
          write(24,102) T*1000D0, X(I),
     *      USOL(1,I,2),USOL(2,I,2),USOL(3,I,2),
     *      USOL(4,I,2),USOL(5,I,2),USOL(6,I,2),USOL(7,I,2),
     *      UFX,UFY,UFZ

        enddo
        write(24,*)''
        write(21,*)''
  101   format(F7.1 F10.6, 7(1PE15.6))
  102   format(F7.1 F10.6, 10(1PE15.6))
      end

      subroutine WRITEMJ_OPEN()
        include 'par.fh'
        implicit real*8(A-H,O-Z)
        integer FILES_MJ(NPTS), FILE_MMJ
        common /FILES/ FILES_MJ, FILE_MMJ
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /CFG_CELL/ CELL_LEN
        common /CFG_WRITE/ WRITEMJ_XSTEP
        common /TIMEP/ T, TEND
        real*8 X0
        character*64 FNAME
        integer I
        FILE_MMJ=1000
        if (T.GT.1D-10)  then
          open(FILE_MMJ, FILE='mj_mean.dat', ACCESS='APPEND')
        else
          open(FILE_MMJ, FILE='mj_mean.dat')
        endif
        X0=0D0
        do I=1,NPTS
          if (X(I).GE.X0) then
            FILES_MJ(I)=1000+I
            write(FNAME,'(A,F6.4,A)') 'mj',X(I),'.dat'
            if (T.GT.1D-10)  then
              open(FILES_MJ(I),FILE=FNAME, ACCESS='APPEND')
            else
              open(FILES_MJ(I),FILE=FNAME)
              write(FILES_MJ(I),'(A,F6.3AI4)') '# X= ', X(I), ' I = ', I
            endif
            X0 = X0 + WRITEMJ_XSTEP*CELL_LEN
          else
            FILES_MJ(I)=0
          endif
        enddo
      end

      subroutine WRITEMJ_DO()
        implicit real*8(A-H,O-Z)
        include 'par.fh'
        common /BLK_UMU/ T11,GW,W,AA,TF,DIFF,WY,TW,
     +   AF0,TS,TIME_STEP,XS
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /TIMEP/ T, TEND
        integer FILES_MJ(NPTS), FILE_MMJ
        common /FILES/ FILES_MJ, FILE_MMJ

        do I=1, NPTS
          if (FILES_MJ(I).NE.0) then

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
     *       (CTM*U4*U6+U5*ST)*FTN
            UJY=2.0D0*(UX7*U5+ST*UX5-CTM*(U4*UX6-UX4*U6))+
     *       (CTM*U5*U6-U4*ST)*FTN
            UJZ=2.0D0*(UX7*U6+ST*UX6+CTM*(U4*UX5-UX4*U5))+
     *       (CTM*U6**2+CT)*FTN
            UFX=UJX*AF-DIFF*UX1
            UFY=UJY*AF-DIFF*UX2
            UFZ=UJZ*AF-DIFF*UX3

            write(FILES_MJ(I),101) T*1000D0,
     *        USOL(1,I,1),USOL(2,I,1),USOL(3,I,1),
     *        USOL(4,I,1),USOL(5,I,1),USOL(6,I,1),USOL(7,I,1)
C            write(24,102) T*1000., X(I),
C     *        USOL(1,I,2),USOL(2,I,2),USOL(3,I,2),
C     *        USOL(4,I,2),USOL(5,I,2),USOL(6,I,2),USOL(7,I,2),
C     *        UFX,UFY,UFZ
            flush(FILES_MJ(I))
          endif
        enddo
C       write(24,*)''
  101   format(F7.1 F10.6, 7(1PE15.6))
  102   format(F7.1 F10.6, 10(1PE15.6))
      end

CCC   STATE DUMP/RESTORE

      subroutine STATE_DUMP(FNAME)
        implicit real*8(A-H,O-Z)
        include 'par.fh'
        common /BLK_UMU/ T11,GW,W,AA,TF,DIFF,WY,TW,
     +   AF0,TS,TIME_STEP,XS
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /TIMEP/ T, TEND
        character FNAME*128
        NFILE=2001

        open (NFILE, FILE=FNAME)
        write (NFILE, '("# T,s  LP,cm")')
        write (NFILE, '(F8.5," ",F8.5)') T, SLP+T*SWR
        write (NFILE, '("# USOL(i,j,k):")')
        do I=1, NPTS
          do ID=1, NDERV
            do IP=1, NPDE
              write (NFILE, '(E18.10E3) ',advance='no')
     *                USOL(IP,I,ID)
            enddo
          enddo
          write (NFILE,'()')
        enddo
        close (NFILE)
      end

      subroutine STATE_RESTORE(FNAME)
        implicit real*8(A-H,O-Z)
        include 'par.fh'
        common /BLK_UMU/ T11,GW,W,AA,TF,DIFF,WY,TW,
     +   AF0,TS,TIME_STEP,XS
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /TIMEP/ T, TEND
        character FNAME*128

        NFILE=2001

        open (NFILE, FILE=FNAME)
        read (NFILE, '(A)', END=201)
        read (NFILE, *, END=202) T, SLP
        SWR=0D0
        read (NFILE, '(A)', END=202)
        do I=1, NPTS
          do ID=1, NDERV
            do IP=1, NPDE
              read (NFILE, '(E18.10E3)',END=202,advance='no')
     *                USOL(IP,I,ID)
            enddo
          enddo
          read (NFILE, '()',END=202)
        enddo
        close (NFILE)
        write(*,*) 'Resore state: T=', T
        return

 202    write (*, '("Error: incomplete file vmcw_usol.dat")')
 201    close (NFILE)
        T=0D0
        TIME=0D0
        SWR=0D0
        call SET_ICOND()
      end

CCC   CMD PROCESSING

      subroutine CMD_OPEN()
        common /CMD/ CMD_FILE, CMD_FILE_NAME
        integer CMD_FILE
        character CMD_FILE_NAME*20
        open (CMD_FILE, FILE=CMD_FILE_NAME)
      end

      subroutine CMD_READ()
        common /CMD/ CMD_FILE, CMD_FILE_NAME
        integer CMD_FILE
        character CMD_LINE*128, CMD*20, FNAME*128
        common /TIMEP/ T, TEND
        common /CH_PAR/ SLP,SWR,BETA,MJW,IBN
        real*8 T,TEND,TIME,ARG,SWR,SLP


 303    read (CMD_FILE,'(A)', END=302) CMD_LINE

C       No args
        if (CMD.EQ.'STOP') then
          write(*,*) 'Command: STOP'
          close (CMD_FILE)
          stop
        endif

C       1 str arg
        read (CMD_LINE,*, ERR=302) CMD, FNAME

        if (CMD.eq.'SAVE'.or.CMD.eq.'DUMP') then
          write(*,'(A, A, A, A)') 'Command: SAVE, ',
     *       ' FILE = ', FNAME, 's'
          call STATE_DUMP(FNAME)
          return
        endif

        if (CMD.eq.'RESTORE') then
          write(*,'(A, A, A ,A)') 'Command: RESTORE, ',
     *       ' FILE = ', FNAME, 's'
          call STATE_RESTORE(FNAME)
          return
        endif

C       1 arg
        read (CMD_LINE,*, ERR=302) CMD, TIME

        if (CMD.EQ.'WAIT') then
          write(*,'(A, A, F5.3,A)') 'Command: WAIT, ',
     *       ' T = ', TIME, 's'
          TEND=T+TIME
          SLP=SLP+T*SWR
          SWR=0D0
          return
        endif

C       2 arg
        read (CMD_LINE,*, ERR=302) CMD, TIME, ARG

        if (CMD.EQ.'SWEEP_LP') then
          write(*,'(A, F5.3, A, F5.3, A)') 'Command: SWEEP_LP, ',
     *       ARG, ' cm/s, T = ', TIME, 's'
          TEND=T+TIME
          SLP=SLP+T*(SWR-ARG)
          SWR=ARG
          return
        endif

        write(*,*) 'Unknown command: ', CMD_LINE
        goto 303

 302    write(*,*) 'All commands processed.'
        close (CMD_FILE)
        stop
 301    write(*,*) 'Error in CMD file.'
        close (CMD_FILE)
        stop
      end

      subroutine CMD_CLOSE()
        common /CMD/ CMD_FILE, CMD_FILE_NAME
        integer CMD_FILE
        close (CMD_FILE)
      end
