CCC   CMD PROCESSING

      subroutine CMD_OPEN()
        include 'vmcw.fh'
        include 'vmcw_cmd.fh'
        if (INTERACTIVE.eq.0) open (CMD_FILE, FILE=CMD_FILE_NAME)
      end

      subroutine CMD_CLOSE()
        include 'vmcw.fh'
        include 'vmcw_cmd.fh'
        if (INTERACTIVE.eq.0) close (CMD_FILE)
      end

      subroutine STOP_SWEEP()
        include 'vmcw.fh'
        real*8 T, TSTEP, TEND
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
        include 'vmcw_cmd.fh'

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
          call STATE_DUMP(FNAME, USOL, XSOL)
          goto 303 ! next command
        endif

CC command RESTORE <filename> -- restore state from file
        if (CMD.eq.'RESTORE') then
          write(*,'(A, A, A30)') '> RESTORE, ',
     *       ' FILE = ', FNAME
          call STATE_RESTORE(FNAME, USOL, XSOL)
          goto 303 ! next command
        endif

CC command WRITE_M <filename> -- write Mx,My,Mz to file
        if (CMD.eq.'WRITE_M') then
          write(*,'(A, A, A30)') '> WRITE_M, ',
     *       ' FILE = ', FNAME
          open(M_FILE,FILE=FNAME,ERR=310)
          write(M_FILE,*) '# FLEGG:  ', LF0 + LF_SWR
          write(M_FILE,*) '# CPAR:   ', CPAR0 + CPAR_SWR*T
          write(M_FILE,*) '# Diff:   ', DF0 + DF_SWR*T
          write(M_FILE,*) '# TF:     ', TF0 + TF_SWR*T
          write(M_FILE,*) '# H:      ', H
          write(M_FILE,*) '# GRAD:   ', GRAD
          write(M_FILE,*) '# HR:     ', HR0 + T*HR_SWR
          write(M_FILE,*) '# LP:     ', LP0 + T*LP_SWR

          write(M_FILE,*) '# TIME  LP  TMAB  TMDS  TMPC  TMZ'
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
        call WRITEMJ_CLOSE()
        stop 'DONE'
      end
