C---------------- CB=0.0 !!!!!!!!!!
        include 'vmcw.fh'
        include 'par.fh'
        include 'he3_const.fh'
        common /TIMEP/ T, TSTEP, TEND
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),XSOL(NPTS)

        character*64 CFG_KEY

        integer FILES_MJ(NPTS), FILES_MJ0
        common /FILES/ FILES_MJ, FILES_MJ0

        character CMD_FILE_NAME*8  ! file for reading commands
        integer   CMD_FILE          ! file descriptor
        common /CMD_FILE/ CMD_FILE, INTERACTIVE, CMD_FILE_NAME
        data CMD_FILE_NAME/'vmcw.cmd'/, CMD_FILE/200/,INTERACTIVE/0/

        integer   M_FILE ! file for writing Mx,My,Mz
        common /M_FILE/ M_FILE
        data M_FILE /201/

        real*8 WRITEMJ_XSTEP
        common /CFG_WRITE/ WRITEMJ_XSTEP

        common /PDECOL_DATA/ INDEX,MF,SCTCH,WORK,IWORK,
     *    PDECOL_ACC,PDECOL_ACC_LOG2,T0,DT
        integer INDEX,IWORK(IDIMIWORK)
        real*8 SCTCH(KORD*(NDERV+1)),WORK(IDIMWORK)
        real*8 PDECOL_ACC, PDECOL_ACC_LOG2, T0, DT

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

        LP0=0D0
        HR0=1D-3
        LP_SWR=0D0
        HR_SWR=0D0

        call pdecol_init(T) ! set PDECOL parameters

        call SET_MESH(XSOL, NPTS)
        call SAVE_MESH(XSOL, NPTS, 'mesh.dat')
        call SET_ICOND()

        call WRITEMJ_OPEN()

C--------------- COMPUTE PARAMETERS ----------------------------
        call CMD_OPEN()
        call SET_HE3PT(PRESS,TTC,T1C)
C----------------MAIN LOOP -------------------------------------------
        TEND=T
   2    CONTINUE
          if (dabs(TTC_ST).ge.1D-5) then
            TTC=TTC+TTC_ST
            call SET_HE3PT(PRESS,TTC,T1C)
          endif
          if (T.ge.TEND) call CMD_READ()
          T=T+TSTEP
          call pdecol_run(T)
          call MONITOR()
        goto 2
      end

      subroutine pdecol_init(T)
        include 'vmcw.fh'
        include 'par.fh'
        real*8 T

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

      subroutine pdecol_run(T)
        include 'vmcw.fh'
        include 'par.fh'
        real*8 T

        common /PDECOL_DATA/ INDEX,MF,SCTCH,WORK,IWORK,
     *    PDECOL_ACC,PDECOL_ACC_LOG2,T0,DT
        integer INDEX,IWORK(IDIMIWORK)
        real*8 SCTCH(KORD*(NDERV+1)),WORK(IDIMWORK)
        real*8 PDECOL_ACC, PDECOL_ACC_LOG2, T0, DT

        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),XSOL(NPTS)

        call PDECOL(T0,T,DT,XSOL,PDECOL_ACC,NINT,KORD,NCC,NPDE,MF,
     +              INDEX,WORK,IWORK)
        if(INDEX.NE.0) THEN
          write(*,*) 'INTEGRATION FAILED; INDEX=', INDEX
          stop
        endif

        call VALUES(XSOL,USOL,SCTCH,NPDE,NPTS,NPTS,2,WORK)
      end

C-- F ---------- EVALUATION OF F ------------------------------------
      subroutine F(T,X,U,UX,UXX,FV,NPDE)
        include 'vmcw.fh'
        include 'he3_const.fh'

        dimension U(NPDE),UX(NPDE),UXX(NPDE),FV(NPDE)
C       T - time
C       X - x-coord
C       U   - Mx My Mz Nx Ny Nz T 
C       UX  - dU/dx
C       UXX - d2U/dx2
C       FV  - result

C       calculate freq
        WY = GAM*(HR0+HR_SWR*T)

        WL = GAM*(H + GRAD*X)
        W0= GAM*(H + GRAD*(LP0+LP_SWR*T))
        DW = WL-W0

C       fix n vector length
        UN=DSQRT(U(4)**2+U(5)**2+U(6)**2)
        UNx = U(4)/UN
        UNy = U(5)/UN
        UNz = U(6)/UN

        UMxm = U(1)-WY/WL
        UMym = U(2)
        UMzm = U(3)-1.0D0

        WL2 = 0.5D0*WL
        DD45=UNx*UX(5)-UX(4)*UNy       ! Nx Ny` - Nx` Ny
        ST=DSIN(U(7))
        CT=DCOS(U(7))
        CTM=1.0D0-CT
        CT1=1.0D0+CT
        CTG=ST/CTM      ! ctg(T/2) = sin(T)/(1-cos(T))
        UT=ST*(1.0D0+4.0D0*CT)*0.2666666D0

        AUT=UT*(LF0+LF_SWR*T)**2/WL*4.0D0*PI**2
        AF=-(CPAR0+CPAR_SWR*T)**2/WL
        TF1=TF0+TF_SWR*T
        DIFF=DF0+DF_SWR*T


        AF=AF*(1D0 - 0.5D0 * AER_STEP(X,0))
        DAF=-AF*0.5D0 * AER_STEP(X,1)
        AUT=AUT*(1D0 - 0.835D0 * AER_STEP(X,0))
        TF1=TF1*(1D0 - 0.5D0 * AER_STEP(X,0))

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

        FV(1)=   DW*U(2)           + AUT*UNx - DJX
        FV(2)= - DW*U(1) + WY*U(3) + AUT*UNy - DJY
        FV(3)=           - WY*U(2) + AUT*UNz - DJZ - UMzm*T11

        FV(4)= - W0*UNy - WL2*(UMzm*UNy-UMym*UNz+CTG*(B*UNx-UMxm))
        FV(5)=   W0*UNx - WL2*(UMxm*UNz-UMzm*UNx+CTG*(B*UNy-UMym))
        FV(6)=           - WL2*(UMym*UNx-UMxm*UNy+CTG*(B*UNz-UMzm))
        FV(7)= WL*B + UT/TF1
        return
      end

C-- BNDRY ------ BOUNDARY CONDITIONS -- B(U,UX)=Z(T) ------------
      subroutine BNDRY(T,X,U,UX,DBDU,DBDUX,DZDT,NPDE)
        include 'vmcw.fh'
        include 'he3_const.fh'
        dimension U(NPDE),UX(NPDE),DZDT(NPDE),
     *   DBDU(NPDE,NPDE),DBDUX(NPDE,NPDE)
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

          W=GAM*H
          AF=-(CPAR0+CPAR_SWR*T)**2/W

          AF=AF*(1D0 - 0.5D0 * AER_STEP(X,0))

          DA=-(DF0+DF_SWR*T)/AF

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
        include 'vmcw.fh'
        include 'par.fh'
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
        include 'vmcw.fh'
        include 'par.fh'
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),XSOL(NPTS)
        do K=1,NPTS
          USM=DSQRT(USOL(5,K,1)**2+USOL(6,K,1)**2+USOL(4,K,1)**2)
          USOL(4,K,1)=USOL(4,K,1)/USM
          USOL(5,K,1)=USOL(5,K,1)/USM
          USOL(6,K,1)=USOL(6,K,1)/USM
        enddo
        AA=1.0D20
        IK=1
        do II=1,NPTS
          BB=(XSOL(II)-XI)**2
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
        include 'vmcw.fh'
        dimension UI(NPDEI)
        do I=1,NPDEI
          UI(I)=USP(XI,I)
        enddo
        return
      end
C-- DERIVF ----- SET UP DERIVATIVES ---------------------------------
      subroutine DERIVF(T,X,U,UX,UXX,DFDU,DFDUX,DFDUXX,NPDE)
        include 'vmcw.fh'
        dimension U(NPDE),UX(NPDE),UXX(NPDE),
     *       DFDU(NPDE,NPDE),DFDUX(NPDE,NPDE),DFDUXX(NPDE,NPDE)
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

CCC   CMD PROCESSING

      subroutine CMD_OPEN()
        include 'vmcw.fh'
        common /CMD_FILE/ CMD_FILE,INTERACTIVE,CMD_FILE_NAME
        integer CMD_FILE
        character CMD_FILE_NAME*8
        if (INTERACTIVE.eq.0) open (CMD_FILE, FILE=CMD_FILE_NAME)
      end

      subroutine CMD_CLOSE()
        include 'vmcw.fh'
        common /CMD_FILE/ CMD_FILE,INTERACTIVE, CMD_FILE_NAME
        integer CMD_FILE
        if (INTERACTIVE.eq.0) close (CMD_FILE)
      end

      subroutine STOP_SWEEP()
        include 'vmcw.fh'
        common /TIMEP/ T, TSTEP, TEND
        HR0=HR0+T*HR_SWR
        HR_SWR=0D0
        LP0=LP0+T*LP_SWR
        LP_SWR=0D0
        DF0=DF0+T*DF_SWR
        DF_SWR=0D0
        TF0=TF0+T*TF_SWR
        TF_SWR=0D0
        LF0=LF0+T*LF_SWR
        LF_SWR=0D0
        CPAR0=CPAR0+T*CPAR_SWR
        CPAR_SWR=0D0
      end

      subroutine CMD_READ()
        include 'vmcw.fh'

        integer CMD_FILE
        common /CMD_FILE/ CMD_FILE, INTERACTIVE, CMD_FILE_NAME

        integer   M_FILE
        common /M_FILE/ M_FILE

        character CMD_LINE*128, CMD*128, FNAME*128
        common /TIMEP/ T, TSTEP, TEND
        real*8 T,TSTEP,TEND,ARG1,ARG2

        real*8 LP,HR,STEPS

        call STOP_SWEEP !!!

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
          open(M_FILE,FILE=FNAME,ERR=310)
          write(M_FILE,*), '# FLEGG:  ', LF0 + LF_SWR
          write(M_FILE,*), '# CPAR:   ', CPAR0 + CPAR_SWR*T
          write(M_FILE,*), '# Diff:   ', DF0 + DF_SWR*T
          write(M_FILE,*), '# TF:     ', TF0 + TF_SWR*T
          write(M_FILE,*), '# H:      ', H
          write(M_FILE,*), '# GRAD:   ', GRAD
          write(M_FILE,*), '# HR:     ', HR0 + T*HR_SWR
          write(M_FILE,*), '# LP:     ', LP0 + T*LP_SWR

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
          if (ARG1.eq.0D0) goto 303 ! next command
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
CC command DF <cm^2/s> -- set diffusion
        if (CMD.EQ.'DF') then
          write(*,'("> ", A12, F5.3, A)') CMD, ARG1, ' cm^2/s'
          DF0=ARG1 - T*DF_SWR
          goto 303 ! next command
        endif
CC command TF <s> -- set T_eff
        if (CMD.EQ.'TF') then
          write(*,'("> ", A12, F5.3, A)') CMD, ARG1, ' s'
          TF0=ARG1 - T*TF_SWR
          goto 303 ! next command
        endif
CC command LF <s> -- set leggett freq
        if (CMD.EQ.'LF') then
          write(*,'("> ", A12, F5.3, A)') CMD, ARG1, ' Hz'
          LF0=ARG1 - T*LF_SWR
          goto 303 ! next command
        endif
CC command CPAR <s> -- set C_par
        if (CMD.EQ.'CPAR') then
          write(*,'("> ", A12, F5.3, A)') CMD, ARG1, ' cm/s'
          CPAR0=ARG1 - T*CPAR_SWR
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
          ARG2=dabs(ARG2)
          STEPS = dabs(dfloat(int((ARG1-LP0)/ARG2/TSTEP)))
          write(*,*) ' do ', int(STEPS), ' time steps'
          if (STEPS.eq.0D0) goto 303
          LP_SWR = (ARG1-LP0)/(STEPS*TSTEP)
          write(*,'(A,F8.5,A)') ' real rate: ', LP_SWR,  ' cm/s'
          LP0  = LP0 - T*LP_SWR
          TEND = T + (STEPS - 1D-4)*TSTEP
          return
        endif

CC command HR_SWEEP_TO <value, mOe> <rate, mOe/s> -- sweep RF-field
        if (CMD.EQ.'HR_SWEEP_TO') then
          write(*,'("> ", A12, F5.3, A, F5.3, A)') CMD,
     *       ARG1, ' mOe with rate ', ARG2, ' mOe/s'
          ARG1=ARG1*1D-3
          ARG2=ARG2*1D-3
          ARG2=dabs(ARG2)
          STEPS = dabs(dfloat(int((ARG1-HR0)/ARG2/TSTEP)))
          write(*,*) ' do ', int(STEPS), ' time steps'
          if (STEPS.eq.0D0) goto 303
          HR_SWR = (ARG1-HR0)/(STEPS*TSTEP)
          write(*,'(A,F8.5,A)') ' real rate: ', HR_SWR*1D3, ' mOe/s'
          HR0  = HR0 - T*HR_SWR
          TEND = T + (STEPS - 1D-4)*TSTEP
          return
        endif

CC command DF_SWEEP_TO <diff value, cm^2/s> <rate, cm^2/s^2> -- sweep Diff
        if (CMD.EQ.'DF_SWEEP_TO') then
          write(*,'("> ", A12, F5.3, A, F5.3, A)') CMD,
     *       ARG1, ' cm^2/s with rate ', ARG2, ' cm^2/s^2'
          ARG2=dabs(ARG2)
          STEPS = dabs(dfloat(int((ARG1-DF0)/ARG2/TSTEP)))
          write(*,*) ' do ', int(STEPS), ' time steps'
          if (STEPS.eq.0D0) goto 303
          DF_SWR = (ARG1-DF0)/(STEPS*TSTEP)
          write(*,'(A,F8.5,A)') ' real rate: ', DF_SWR, ' cm^2/s^2'
          DF0  = DF0 - T*DF_SWR
          TEND = T + (STEPS - 1D-4)*TSTEP
          return
        endif

CC command TF_SWEEP_TO <t_eff value, s> <rate> -- sweep t_eff
        if (CMD.EQ.'TF_SWEEP_TO') then
          write(*,'("> ", A12, F5.3, A, F5.3, A)') CMD,
     *       ARG1, ' s with rate ', ARG2, ' s/s'
          ARG2=dabs(ARG2)
          STEPS = dabs(dfloat(int((ARG1-TF0)/ARG2/TSTEP)))
          write(*,*) ' do ', int(STEPS), ' time steps'
          if (STEPS.eq.0D0) goto 303
          TF_SWR = (ARG1-TF0)/(STEPS*TSTEP)
          write(*,'(A,F8.5,A)') ' real rate: ', TF_SWR, ' s/s'
          TF0  = TF0 - T*TF_SWR
          TEND = T + (STEPS - 1D-4)*TSTEP
          return
        endif

CC command LF_SWEEP_TO <legg value, Hz> <rate, Hz/s> -- sweep leggett freq
        if (CMD.EQ.'LF_SWEEP_TO') then
          write(*,'("> ", A12, F5.3, A, F5.3, A)') CMD,
     *       ARG1, ' Hz with rate ', ARG2, ' Hz/s'
          ARG2=dabs(ARG2)
          STEPS = dabs(dfloat(int((ARG1-LF0)/ARG2/TSTEP)))
          write(*,*) ' do ', int(STEPS), ' time steps'
          if (STEPS.eq.0D0) goto 303
          LF_SWR = (ARG1-LF0)/(STEPS*TSTEP)
          write(*,'(A,F8.5,A)') ' real rate: ', LF_SWR, ' s/s'
          LF0  = LF0 - T*LF_SWR
          TEND = T + (STEPS - 1D-4)*TSTEP
          return
        endif

CC command CPAR_SWEEP_TO <legg value, Hz> <rate, Hz/s> -- sweep leggett freq
        if (CMD.EQ.'CPAR_SWEEP_TO') then
          write(*,'("> ", A12, F5.3, A, F5.3, A)') CMD,
     *       ARG1, ' cm/s with rate ', ARG2, ' cm/s^2'
          ARG2=dabs(ARG2)
          STEPS = dabs(dfloat(int((ARG1-CPAR0)/ARG2/TSTEP)))
          write(*,*) ' do ', int(STEPS), ' time steps'
          if (STEPS.eq.0D0) goto 303
          CPAR_SWR = (ARG1-CPAR0)/(STEPS*TSTEP)
          write(*,'(A,F8.5,A)') ' real rate: ', CPAR_SWR, ' cm/s^2'
          CPAR0  = CPAR0 - T*CPAR_SWR
          TEND = T + (STEPS - 1D-4)*TSTEP
          return
        endif

CC command T_P <t> <p> -- set T/P
        if (CMD.EQ.'P_T') then
          write(*,'("> ", A12, F5.3, A, F5.3, A)')
     +      CMD, ARG1, ' bar, ' , ARG2, ' Tc'
          call SET_HE3PT(ARG1, ARG2, 6.0D+14)
          goto 303 ! next command
        endif


 301    write(*,*) 'Unknown command: ', CMD_LINE
        goto 303

 302    write(*,*) 'All commands processed.'
        close (CMD_FILE)
        stop
      end




      subroutine SET_HE3PT(PRESS, TTC, T1C)
        include 'vmcw.fh'
        include 'he3_const.fh'

        call STOP_SWEEP

        TEMP=TTC*TCF(PRESS)
        T1=T1C*TEMP/1000D0         ! RELAXATION TIME
        T11=1.0D0/T1

        CPAR0=CPARF(PRESS,TEMP)          ! SPIN WAVES VELOCITY
        LF0  =DSQRT(LF2F(PRESS,TTC))     ! LEGGETT FREQ
        DF0  =DF(PRESS,TEMP)             ! SPIN DIFFUSION

        TR=1.2D-7/DSQRT(1.0D0-TTC)                  !
C        TR=1.2D-7/DSQRT(1.0D0-0.94D0)!!!
        TF0=1D0/ (4D0*PI**2 *LF2F(PRESS,TTC)*TR)   ! TAU EFFECTIVE (L-T) SECONDS WV pic.10.5 20bar



        write(*,'(" P: ",F5.2," bar (Tc = ",F5.3," mK, Tab = ",F5.3,
     *            " mK), T: ",F5.3," mK = ",F5.3," Tc")'),
     *    PRESS, TCF(PRESS), TABF(PRESS), TEMP, TTC
        write(*,'(" F_legg: ", F9.3, " kHz,  ",
     *            " D: ", E9.2, " cm^2/s, ",
     *            " T_lt: ", E12.6, " s, ",
     *            " C_par: ", F6.1, " cm/s ")'),
     *              LF0/1D3, DF0, TF0, CPAR0

      end

