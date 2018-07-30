C---------------- CB=0.0 !!!!!!!!!!
      function VMCW_F(USOL, XSOL)
        include 'vmcw.fh'
        include 'par.fh'
        include 'he3_const.fh'

        include 'vmcw_cmd.fh'
        integer VMCW_F

        data CMD_FILE_NAME/'vmcw.cmd'/, CMD_FILE/200/,INTERACTIVE/0/

        real*8 T, TSTEP, TEND
        common /TIMEP/ T, TSTEP, TEND

        real*8 USOL(NPDE,NPTS,NDERV),XSOL(NPTS), TEST(10)

        character*64 CFG_KEY
        real*8 CFG_VAL

        integer FILES_MJ(NPTS), FILES_MJ0
        common /FILES/ FILES_MJ, FILES_MJ0

        integer   M_FILE ! file for writing Mx,My,Mz
        common /M_FILE/ M_FILE
        data M_FILE /201/

        real*8 WRITEMJ_XSTEP
        common /CFG_WRITE/ WRITEMJ_XSTEP

C--------------- INITIALIZATION -------------------------------------

        VMCW_F=0
        WRITEMJ_XSTEP=0.1D0

        open(54,FILE='vmcw.cfg')
   11   read(54,*,END=12) CFG_KEY, CFG_VAL
        if (CFG_KEY.EQ.'BETA') then
          BETA=CFG_VAL ! INITIAL TIPPING ANGLE (DEGREES)

        elseif (CFG_KEY.EQ.'T1C') then
          T1C=CFG_VAL  ! T1*T                     (1.0E-3)

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

c        T=0D0
c        TSTEP=5D-3

c        LP0=0D0
c        HR0=1D-3
c        LP_SWR=0D0
c        HR_SWR=0D0

c          write(*,'(A,A20)') 'init'
c        call pdecol_init(T) ! set PDECOL parameters

c        call SET_MESH(XSOL, NPTS)
c        call SAVE_MESH(XSOL, NPTS)
c        call WRITEMJ_OPEN(USOL,XSOL)

C--------------- COMPUTE PARAMETERS ----------------------------
c        call CMD_OPEN()
c        call SET_HE3PT()
C----------------MAIN LOOP -------------------------------------------
c        TEND=T
c   2    CONTINUE
c          if (dabs(TTC_ST).ge.1D-5) then
c            TTC=TTC+TTC_ST
c            call SET_HE3PT()
c          endif
c          if (T.ge.TEND) call CMD_READ()
c          T=T+TSTEP
c          call pdecol_run(T,USOL, XSOL)
c          call MONITOR(USOL,XSOL)
c        goto 2
      end



      subroutine SET_HE3PT()
        include 'vmcw.fh'
        include 'he3_const.fh'
        real*8 TEMP, T1, TR

        call STOP_SWEEP

        TEMP=TTC*TCF(PRESS)
        T1=T1C*TEMP/1000D0         ! RELAXATION TIME
        T11=1.0D0/T1

        CPAR0=CPARF(PRESS,TEMP)          ! SPIN WAVES VELOCITY
        LF0  =dsqrt(LF2F(PRESS,TTC))     ! LEGGETT FREQ
        DF0  =DF(PRESS,TEMP)             ! SPIN DIFFUSION

        TR=1.2D-7/dsqrt(1.0D0-TTC)                  !
C        TR=1.2D-7/dsqrt(1.0D0-0.94D0)!!!
        TF0=1D0/ (4D0*PI**2 *LF2F(PRESS,TTC)*TR)   ! TAU EFFECTIVE (L-T) SECONDS WV pic.10.5 20bar



        write(*,'(" P: ",F5.2," bar (Tc = ",F5.3," mK, Tab = ",F5.3,
     *            " mK), T: ",F5.3," mK = ",F5.3," Tc")')
     *    PRESS, TCF(PRESS), TABF(PRESS), TEMP, TTC
        write(*,'(" F_legg: ", F9.3, " kHz,  ",
     *            " D: ", E9.2, " cm^2/s, ",
     *            " T_lt: ", E12.6, " s, ",
     *            " C_par: ", F6.1, " cm/s ")')
     *              LF0/1D3, DF0, TF0, CPAR0

      end

