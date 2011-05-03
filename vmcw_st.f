C---------------- CB=0.0 !!!!!!!!!!
        implicit real*8(A-H,O-Z)
        include 'par.fh'
        include 'he3_const.fh'
        common /TIMEP/ T, TSTEP, TEND
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /CH_PAR/ LP0,LP_SWR,BETA,IBN, TTC, TTC_ST
        real*8 LP0,LP_SWR,BETA

        common /CFG_AER/  AER, AER_LEN, AER_CNT, AER_TRW
        common /CFG_CELL/ CELL_LEN
        common /CFG_MESH/ XMESH_K,XMESH_ACC

        common /BLK_UMU/ T11,GRAD,H,AA,TF,DIFF,HR0,HR_SWR,AF0
        character*64 CFG_KEY

        integer FILES_MJ(NPTS) ! files for writing MJ along the cell
        common /FILES/ FILES_MJ

        character CMD_FILE_NAME*20  ! file for reading commands
        integer   CMD_FILE          ! file descriptor
        common /CMD_FILE/ CMD_FILE, INTERACTIVE, CMD_FILE_NAME
        data CMD_FILE_NAME/'vmcw.cmd'/, CMD_FILE/200/,INTERACTIVE/0/

        integer   M_FILE ! file for writing Mx,My,Mz
        common /M_FILE/ M_FILE
        data M_FILE /201/

        real*8 WRITEMJ_XSTEP
        common /CFG_WRITE/ WRITEMJ_XSTEP

C       PDECOL parameters
        common /PDECOL_DATA/ INDEX,MF,SCTCH, WORK,IWORK,
     *    PDECOL_ACC, PDECOL_ACC_LOG2,T0,DT
        integer INDEX,IWORK(IDIMIWORK)
        real*8 SCTCH(KORD*(NDERV+1)),WORK(IDIMWORK)
        real*8 PDECOL_ACC, PDECOL_ACC_LOG2,T0,DT


C--------------- INITIALIZATION -------------------------------------

        WRITEMJ_XSTEP=0.1D0

        open(54,FILE='vmcw.cfg')
   11   read(54,*,END=12) CFG_KEY, CFG_VAL
        if (CFG_KEY.EQ.'BETA') then
          BETA=CFG_VAL ! INITIAL TIPPING ANGLE (DEGREES)

        elseif (CFG_KEY.EQ.'T1C') then
          T1C=CFG_VAL  ! T1*T                     (1.0E-3)
        elseif (CFG_KEY.EQ.'PDECOL_ACC_LOG2') then
          PDECOL_ACC_LOG2=CFG_VAL ! RELATIVE TIME ERROR BOUND

        elseif (CFG_KEY.EQ.'TEMP') then
          TTC=CFG_VAL ! TEMPERATURE (T/TC)
        elseif (CFG_KEY.EQ.'PRESS') then
          PRESS=CFG_VAL ! PRESS (bar)

        elseif (CFG_KEY.EQ.'H') then
          H=CFG_VAL    ! FIELD (OE)               (110)
        elseif (CFG_KEY.EQ.'GRAD') then
          GRAD=CFG_VAL ! GRADIENT H (OE/CM)       (-0.2)
        elseif (CFG_KEY.EQ.'HR') then
          HR0=CFG_VAL   ! RF FIELD (OE)            (0.06)

        elseif (CFG_KEY.EQ.'IBN') then
          IBN=INT(CFG_VAL)  ! TYPE OF BOUND. COND.: 1-OPEN CELL 2-CLOSED CELL

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
        elseif (CFG_KEY.EQ.'INTERACTIVE') then
          INTERACTIVE=int(CFG_VAL)


        else
          write(*,'(A,A20)')
     *     'warning: unknown parameter in cfg-file: ', CFG_KEY
        endif
        goto 11
   12   close(54)

        T=0D0
        TSTEP=5D-3
        TEND=0D0

        LP0=0D0
        HR0=1D-3
        LP_SWR=0D0
        HR_SWR=0D0

        call PDECOL_INIT(T) ! set PDECOL parameters

        call SET_MESH()
        call SAVE_MESH('mesh.dat')
        call SET_ICOND()

        call WRITEMJ_OPEN()

C--------------- COMPUTE PARAMETERS ----------------------------
        PI=4.0D0*DATAN(1.0D0)
        call CMD_OPEN()
        call SET_HE3PT(PRESS,TTC,T1C)
C----------------MAIN LOOP -------------------------------------------
   2    CONTINUE

          if (dabs(TTC_ST).ge.1D-5) then
            TTC=TTC+TTC_ST
            call SET_HE3PT1(PRESS,TTC,T1C)
          endif

          if(T.GE.TEND) call CMD_READ()
          T=T+TSTEP

          call PDECOL(T0,T,DT,X,PDECOL_ACC,NINT,KORD,NCC,NPDE,MF,
     +                INDEX,WORK,IWORK)
          if(INDEX.NE.0) THEN
            write(*,*) 'INTEGRATION FAILED; INDEX=', INDEX
            stop
          endif

          call VALUES(X,USOL,SCTCH,NPDE,NPTS,NPTS,2,WORK)
          call MONITOR()
        goto 2
      end
C-- F ---------- EVALUATION OF F ------------------------------------
      subroutine F(T,X,U,UX,UXX,FV,NPDE)
        implicit real*8(A-H,O-Z)
        dimension U(NPDE),UX(NPDE),UXX(NPDE),FV(NPDE)
        common /CH_PAR/ LP0,LP_SWR,BETA,IBN,TTC,TTC_ST
        real*8 LP0,LP_SWR,BETA
        common /BLK_UMU/ T11,GRAD,H,AA,TF,DIFF,HR0,HR_SWR,AF0
        common /HE3DATA/ PCP,TPCP,PA,ANA,PI1,HC,R,AKB,GAM,AM3
C       T - time
C       X - x-coord
C       U   - Mx My Mz Nx Ny Nz T 
C       UX  - dU/dx
C       UXX - d2U/dx2
C       FV  - result
C       calculate freq

        WY = GAM*(HR0+HR_SWR*T)
        WZ = GAM*GRAD*(LP0+LP_SWR*T)
        WR = GAM*(H+GRAD*X)
        WZR= GAM*H + WZ
        XZ=  GAM*GRAD*(X-LP0-LP_SWR*T)

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

C        AF=AF0 - AF0*0.5D0 * AER_STEP(X,0)
C        DAF=-AF0*0.5D0 * AER_STEP(X,1)
C        AUT=AUT0 - AUT0*0.835D0 * AER_STEP(X,0)
C        TF0=TF-TF*0.5D0 * AER_STEP(X,0)

        AUT=AUT0 - AUT0*0.0D0 * AER_STEP(X,0)
        AF=AF0
        DAF=0D0
        TF0=TF

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
        common /CH_PAR/ LP0,LP_SWR,BETA,IBN,TTC,TTC_ST
        real*8 LP0,LP_SWR,BETA
        common /BLK_UMU/ T11,GRAD,H,AA,TF,DIFF,HR0,HR_SWR,AF0
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
        common /CH_PAR/ LP0,LP_SWR,BETA,IBN,TTC,TTC_ST
        real*8 LP0,LP_SWR,BETA
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
        common /CH_PAR/ LP0,LP_SWR,BETA,IBN,TTC,TTC_ST
        real*8 LP0,LP_SWR,BETA
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
        USP=USOL(I,IK,1)
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
        common /CH_PAR/ LP0,LP_SWR,BETA,IBN,TTC,TTC_ST
        real*8 LP0,LP_SWR,BETA
        common /BLK_UMU/ T11,GRAD,H,AA,TF,DIFF,HR0,HR_SWR,AF0
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
        implicit real*8(A-H,O-Z)
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
        common /TIMEP/ T, TSTEP, TEND
        common /GEAR0/ DTUSED,NQUSED,NSTEP,NFE,NJE
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /CH_PAR/ LP0,LP_SWR,BETA,IBN,TTC,TTC_ST
        real*8 LP0,LP_SWR,BETA
        common /BLK_UMU/ T11,GRAD,H,AA,TF,DIFF,HR0,HR_SWR,AF0
        common /CFG_CELL/ CELL_LEN

        integer   M_FILE
        common /M_FILE/ M_FILE

C--------------- COMPUTE TIME DEPENDENCIES --------------------------
        TMMS=T*1000.0D0

        TMLP=(LP0+LP_SWR*T)/CELL_LEN
        TMAB=0.0D0
        TMDS=0.0D0
        TMZ=0.0D0
        do I=1,NPTS-1
          TMAB=TMAB + USOL(1,I,1) * (X(I+1)-X(I))
          TMDS=TMDS + USOL(2,I,1) * (X(I+1)-X(I))
          TMZ=TMZ   + USOL(3,I,1) * (X(I+1)-X(I))
        enddo
        TMAB=TMAB/CELL_LEN
        TMDS=TMDS/CELL_LEN
        TMZ=TMZ/CELL_LEN
        TMPC=(DSQRT(TMAB**2+TMDS**2))
        write(M_FILE, 61)
     +   TMMS,TMLP,TMAB,TMDS,TMPC,TMZ
C--------------- SHOW INFORMATION -----------------------------------
        write(*,'(A,F6.1,A,A,F6.3,A,F6.3,A,A,F6.3,A,A,F6.3)')
     *    ' T=',TMMS,' ms, ',
     *    'LP=', LP0+LP_SWR*T , ' cm = ', TMLP , ' cell, ',
     *    'HR=',  1D3*(HR0+HR_SWR*T) , ' mOe, ',
     *    'TTC=',  TTC
        call WRITEMJ_DO()
   61   format(F7.1, 6(F12.8))
   69   format(F9.2, 1PE25.16)
      end
C-- WRITE_MJ --- WRITE SPINS & CURRENTS TO VMCW ------------------

      subroutine WRITEMJ_OPEN()
        include 'par.fh'
        implicit real*8(A-H,O-Z)
        integer FILES_MJ(NPTS)
        common /FILES/ FILES_MJ
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /CFG_CELL/ CELL_LEN
        common /CFG_WRITE/ WRITEMJ_XSTEP
        common /TIMEP/ T, TSTEP, TEND
        real*8 X0
        character*64 FNAME
        integer I
        X0=0D0
        do I=1,NPTS
          if (X(I).GE.X0) then
            FILES_MJ(I)=1000+I
            write(FNAME,'(A,F6.4,A)') 'mj',X(I),'.dat'
            open(FILES_MJ(I),FILE=FNAME)
            write(FILES_MJ(I),'(A,F6.3AI4)') '# X= ', X(I), ' I = ', I
            X0 = X0 + WRITEMJ_XSTEP*CELL_LEN
          else
            FILES_MJ(I)=0
          endif
        enddo
      end

      subroutine WRITEMJ_DO()
        implicit real*8(A-H,O-Z)
        include 'par.fh'
        common /BLK_UMU/ T11,GRAD,H,AA,TF,DIFF,HR0,HR_SWR,AF0
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /TIMEP/ T, TSTEP, TEND
        integer FILES_MJ(NPTS)
        common /FILES/ FILES_MJ
        common /CH_PAR/ LP0,LP_SWR,BETA,IBN,TTC,TTC_ST
        real*8 LP0,LP_SWR,BETA
        common /CFG_CELL/ CELL_LEN

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

            write(FILES_MJ(I),101) T*1000D0, (LP0+LP_SWR*T)/CELL_LEN,
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
        common /BLK_UMU/ T11,GRAD,H,AA,TF,DIFF,HR0,HR_SWR,AF0
        common /CH_PAR/ LP0,LP_SWR,BETA,IBN,TTC,TTC_ST
        real*8 LP0,LP_SWR,BETA
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /TIMEP/ T, TSTEP, TEND
        character FNAME*128
        NFILE=2001

        open (NFILE, FILE=FNAME)
        write (NFILE, '("# T,s  LP,cm")')
        write (NFILE, '(3(F10.5))') T,
     *   LP0+T*LP_SWR, HR0+T*HR_SWR
        write (NFILE, '("# USOL(i,j,k):")')
        do I=1, NPTS
          do ID=1, NDERV
            do IP=1, NPDE
              write (NFILE, '(E19.10E3)', advance='no')
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
        common /BLK_UMU/ T11,GRAD,H,AA,TF,DIFF,HR0,HR_SWR,AF0
        common /CH_PAR/ LP0,LP_SWR,BETA,IBN,TTC,TTC_ST
        real*8 LP0,LP_SWR,BETA
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /TIMEP/ T, TSTEP, TEND
        character FNAME*128

        NFILE=2001

        open (NFILE, FILE=FNAME)
        read (NFILE, '(A)', END=201)
        read (NFILE, *, END=202) T, LP0, HR0
        read (NFILE, '(A)', END=202)
        LP_SWR=0D0
        HR_SWR=0D0
        do I=1, NPTS
          do ID=1, NDERV
            do IP=1, NPDE
              read (NFILE, '(E19.10E3)',
     *                ERR=202, END=202, advance='no')
     *                USOL(IP,I,ID)
            enddo
          enddo
          read (NFILE, '()',ERR=202, END=202)
        enddo
        close (NFILE)
        write(*,'("Restore state:", 3(A,F7.4,A))')
     *   ' T = ', T, ' s',
     *   ' LP = ', LP0, ' cm',
     *   ' HR = ', HR0*1D3, ' mOe'
        call PDECOL_INIT(T) ! reset PDECOL parameters
        return

 202    write (*, '("Error: incomplete file: ",A)'), FNAME
        close (NFILE)
        T=0D0
        LP_SWR=0D0
        HR_SWR=0D0
        call SET_ICOND()
        return
 201    write (*, '("Error: skip file: ",A)'), FNAME
        close (NFILE)
      end

CCC   CMD PROCESSING

      subroutine CMD_OPEN()
        implicit real*8(A-H,O-Z)
        common /CMD_FILE/ CMD_FILE,INTERACTIVE,CMD_FILE_NAME
        integer CMD_FILE
        character CMD_FILE_NAME*20
        if (INTERACTIVE.eq.0) open (CMD_FILE, FILE=CMD_FILE_NAME)
      end

      subroutine CMD_CLOSE()
        common /CMD_FILE/ CMD_FILE,INTERACTIVE, CMD_FILE_NAME
        integer CMD_FILE
        if (INTERACTIVE.eq.0) close (CMD_FILE)
      end

      subroutine CMD_READ()
        implicit real*8(A-H,O-Z)

        integer CMD_FILE
        common /CMD_FILE/ CMD_FILE, INTERACTIVE, CMD_FILE_NAME

        integer   M_FILE
        common /M_FILE/ M_FILE

        character CMD_LINE*128, CMD*128, FNAME*128
        common /TIMEP/ T, TSTEP, TEND
        common /CH_PAR/ LP0,LP_SWR,BETA,IBN,TTC,TTC_ST
        real*8 LP0,LP_SWR,BETA
        real*8 T,TSTEP,TEND,ARG1,ARG2
        common /BLK_UMU/ T11,GRAD,H,AA,TF,DIFF,HR0,HR_SWR,AF0

        real*8 LP,HR,STEPS

 303    continue
        if (INTERACTIVE.eq.0) then
          read (CMD_FILE,'(A)', ERR=302,END=302) CMD_LINE
        else
          write(*,'("> ")',advance='no')
          read (*,'(A)', ERR=302,END=302) CMD_LINE
        endif

        if (CMD_LINE(1:1).eq.'#') goto 303
        if (CMD_LINE.eq.'') goto 303

C       No args
        CMD=CMD_LINE

CC command STOP -- stop program
        if (CMD.EQ.'STOP') then
          write(*,*) '> STOP'
          close (CMD_FILE)
          stop
        endif

CC command WRITE_M_STOP -- stop writing Mx,My,Mz to file
        if (CMD.eq.'WRITE_M_STOP') then
          write(*,'(A)') '> WRITE_M_STOP'
          close(M_FILE)
          goto 303 ! next command
        endif

C       1 str arg
        read (CMD_LINE,*, ERR=301, END=301) CMD, FNAME

CC command SAVE/DUMP <filename> -- save state to file
        if (CMD.eq.'SAVE'.or.CMD.eq.'DUMP') then
          write(*,'(A, A, A30)') '> SAVE, ',
     *       ' FILE = ', FNAME
          call STATE_DUMP(FNAME)
          goto 303 ! next command
        endif

CC command RESTORE <filename> -- restore state from file
        if (CMD.eq.'RESTORE') then
          write(*,'(A, A, A30)') '> RESTORE, ',
     *       ' FILE = ', FNAME
          call STATE_RESTORE(FNAME)
          goto 303 ! next command
        endif

CC command WRITE_M <filename> -- write Mx,My,Mz to file
        if (CMD.eq.'WRITE_M') then
          write(*,'(A, A, A30)') '> WRITE_M, ',
     *       ' FILE = ', FNAME
          PI=3.1415926D0
          W=20378D0*H
          open(M_FILE,FILE=FNAME,ERR=310)
          write(M_FILE,*), '# AA:     ', AA
          write(M_FILE,*), '# FLEGG:  ', dsqrt(W*AA/4.0D0/PI**2)
          write(M_FILE,*), '# AF0:    ', AF0
          write(M_FILE,*), '# CPAR:   ', dsqrt(-W*AF0)
          write(M_FILE,*), '# Diff:   ', DIFF
          write(M_FILE,*), '# TF:     ', TF
          write(M_FILE,*), '# H:      ', H
          write(M_FILE,*), '# GRAD:   ', GRAD
          write(M_FILE,*), '# HR:     ', HR0+T*HR_SWR
          write(M_FILE,*), '# LP:     ', LP0+T*LP_SWR

          write(M_FILE,*), '# TIME  LP  TMAB  TMDS  TMPC  TMZ'
          goto 303 ! next command
 310      write (*,*) 'Error: can''t open file: ', FNAME
          goto 303 ! next command
        endif

C       1 real arg
        read (CMD_LINE,*, ERR=301, END=301) CMD, ARG1

CC command WAIT <time, ms> -- wait
        if (CMD.EQ.'WAIT') then
          write(*,'("> ", A12, F8.2, A)') CMD, ARG1, ' ms'
          LP0=LP0+T*LP_SWR
          LP_SWR=0D0
          HR0=HR0+T*HR_SWR
          HR_SWR=0D0
          TEND=T+ARG1*1D-3
          return
        endif

CC command TSTEP <step, ms> -- set time step
        if (CMD.EQ.'TSTEP') then
          write(*,'("> ", A12, F8.5, A)') CMD, ARG1, ' ms'
          TSTEP=ARG1*1D-3
          goto 303 ! next command
        endif

CC command LP <position, cm> -- set larmor position
        if (CMD.EQ.'LP') then
          write(*,'("> ", A12, F8.5, A)') CMD, ARG1, ' cm'
          LP0=ARG1-T*LP_SWR
          goto 303 ! next command
        endif

CC command LP_ADD <distance, cm> -- add value to larmor position
        if (CMD.EQ.'LP_ADD') then
          write(*,'("> ", A12, F8.5, A)') CMD, ARG1, ' cm'
          LP0=LP0+ARG1
          goto 303 ! next command
        endif

CC command LP_RATE <rate, cm/s> -- set larmor position sweep rate
C        if (CMD.EQ.'LP_RATE') then
C          write(*,'("> ", A12, F5.3, A)') CMD, ARG1, ' cm/s'
C          LP0=LP0+T*(LP_SWR-ARG1)
C          LP_SWR=ARG1
C          goto 303 ! next command
C        endif

CC command HR <mOe> -- set rf field
        if (CMD.EQ.'HR') then
          write(*,'("> ", A12, F5.3, A)') CMD, ARG1, ' mOe'
          HR0=ARG1*1D-3 - T*HR_SWR
          goto 303 ! next command
        endif

CC command TTC_CHANGE <Tstep/Tc> -- set TTC step
        if (CMD.EQ.'TTC_CHANGE') then
          write(*,'("> ", A12, F5.3, A)') CMD, ARG1, ' Tc'
          TTC_ST=ARG1
          goto 303 ! next command
        endif

CC command HR_ADD <mOe> -- add value to RF-field
        if (CMD.EQ.'HR_ADD') then
          write(*,'("> ", A12, F5.3, A)') CMD, ARG1, ' mOe'
          HR0=HR0+ARG1*1D-3
          goto 303 ! next command
        endif

CC command HR_RATE <rate, mOe/s> -- set RF-field sweep rate
C        if (CMD.EQ.'HR_RATE') then
C          write(*,'("> ", A12, F5.3, A)') CMD, ARG1, ' mOe/s'
C          ARG1=ARG1*1D-3
C          HR0=HR0+T*(HR_SWR-ARG1)
C          HR_SWR=ARG1
C          goto 303 ! next command
C        endif

C       2 real args
        read (CMD_LINE,*, ERR=301, END=301) CMD, ARG1, ARG2


CC command LP_SWEEP_TO <value, cm> <rate, cm/s> -- sweep larmor position
        if (CMD.EQ.'LP_SWEEP_TO') then
          write(*,'("> ", A12, F5.3, A, F5.3, A)') CMD,
     *       ARG1, ' cm with rate ', ARG2, 'cm/s'
          LP = LP0 + T*LP_SWR
          ARG2=dabs(ARG2)
          STEPS = dabs(dfloat(int((ARG1-LP)/ARG2/TSTEP)))
          write(*,*) ' do ', int(STEPS), ' time steps'
          if (STEPS.eq.0D0) goto 303
          LP_SWR = (ARG1-LP)/(STEPS*TSTEP)
          write(*,'(A,F8.5,A)') ' real rate: ', LP_SWR,  ' cm/s'
          LP0  = LP - T*LP_SWR
          HR0=HR0+T*HR_SWR
          HR_SWR=0D0
          TEND = T + STEPS*TSTEP
          return
        endif

CC command HR_SWEEP_TO <value, mOe> <rate, mOe/s> -- sweep RF-field
        if (CMD.EQ.'HR_SWEEP_TO') then
          write(*,'("> ", A12, F5.3, A, F5.3, A)') CMD,
     *       ARG1, ' mOe with rate ', ARG2, ' mOe/s'
          ARG1=ARG1*1D-3
          ARG2=ARG2*1D-3
          HR = HR0 + T*HR_SWR
          ARG2=dabs(ARG2)
          STEPS = dabs(dfloat(int((ARG1-HR)/ARG2/TSTEP)))
          write(*,*) ' do ', int(STEPS), ' time steps'
          if (STEPS.eq.0D0) goto 303
          HR_SWR = (ARG1-HR)/(STEPS*TSTEP)
          write(*,'(A,F8.5,A)') ' real rate: ', HR_SWR*1D3, ' mOe/s'
          HR0  = HR - T*HR_SWR
          LP0=LP0+T*LP_SWR
          LP_SWR=0D0
          TEND = T + STEPS*TSTEP
          return
        endif


 301    write(*,*) 'Unknown command: ', CMD_LINE
        goto 303

 302    write(*,*) 'All commands processed.'
        close (CMD_FILE)
        stop
      end


      subroutine PDECOL_INIT(T)
        implicit real*8(A-H,O-Z)
        real*8 T
        include 'par.fh'
        common /PDECOL_DATA/ INDEX,MF,SCTCH,WORK,IWORK,
     *    PDECOL_ACC,PDECOL_ACC_LOG2,T0,DT
        integer INDEX,IWORK(IDIMIWORK)
        real*8 SCTCH(KORD*(NDERV+1)),WORK(IDIMWORK)
        real*8 PDECOL_ACC, PDECOL_ACC_LOG2, T0, DT

        INDEX=1  ! TYPE OF CALL (FIRST CALL)
        MF=22
        IWORK(1)=IDIMWORK
        IWORK(2)=IDIMIWORK
        do I=1,IDIMWORK
          WORK(I)=0.0D0
        enddo
        PDECOL_ACC = 2.0D0**(-PDECOL_ACC_LOG2)
        T0=T       ! STARTING TIME FOR PDECOL
        DT=1.D-10  ! INITIAL STEP SIZE IN T
      end


      subroutine SET_HE3PT(PRESS, TTC, T1C)
        implicit real*8(A-H,O-Z)
        include 'he3_const.fh'
        common /BLK_UMU/ T11,GRAD,H,AA,TF,DIFF,HR0,HR_SWR,AF0

        W=GAM*H
        TEMP=TTC*TCF(PRESS)
        T1=T1C*TEMP/1000D0         ! RELAXATION TIME
        T11=1.0D0/T1

        TETC=DSQRT(1.0D0-TTC)
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
     *            " D: ", E9.2, " cm^2/s, ",
     *            " T_lt: ", E8.2, " s, ",
     *            " C_par: ", F6.1, " cm/s ")'),
     *              FLEG/1D3, DIFF, TF, CPAR

      end

      subroutine SET_HE3PT1(PRESS, TTC, T1C)
        implicit real*8(A-H,O-Z)
        include 'he3_const.fh'
        common /BLK_UMU/ T11,GRAD,H,AA,TF,DIFF,HR0,HR_SWR,AF0

        W=GAM*H
        TEMP=TTC*TCF(PRESS)
        T1=T1C*TEMP/1000D0         ! RELAXATION TIME
        T11=1.0D0/T1
        TETC=DSQRT(1.0D0-TTC)

C        CPAR=CPARF(PRESS,TEMP)          ! SPIN WAVES VELOCITY
C        AF0=-CPAR**2/W

        FLEG=DSQRT(LF2F(PRESS,TTC))     ! LEGGETT FREQ
        AA=FLEG*FLEG/W*4.0D0*PI**2

C        DIFF=DF(PRESS,TEMP)             ! SPIN DIFFUSION

C        TF=1.2D-7/TETC                  ! TAU EFFECTIVE (L-T) SECONDS WV pic.10.5 20bar


        write(*,'(" P: ",F5.2," bar (Tc = ",F5.3," mK, Tab = ",F5.3,
     *            " mK), T: ",F5.3," mK = ",F5.3," Tc")'),
     *    PRESS, TCF(PRESS), TABF(PRESS), TEMP, TTC
        write(*,'(" F_legg: ", F9.3, " kHz,  ",
     *            " D: ", E9.2, " cm^2/s, ",
     *            " T_lt: ", E8.2, " s, ",
     *            " C_par: ", F6.1, " cm/s ")'),
     *              FLEG/1D3, DIFF, TF, CPAR

      end
