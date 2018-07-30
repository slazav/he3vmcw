CCC   STATE DUMP/RESTORE (does not work?)

      subroutine STATE_DUMP(FNAME, USOL, XSOL)
        include 'vmcw.fh'
        include 'par.fh'
        dimension USOL(NPDE,NPTS,NDERV),X(NPTS)
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

      subroutine STATE_RESTORE(FNAME, USOL, XSOL)
        include 'vmcw.fh'
        include 'par.fh'
        dimension USOL(NPDE,NPTS,NDERV),X(NPTS)
        common /TIMEP/ T, TSTEP, TEND
        character FNAME*128

        NFILE=2001

        open (NFILE, FILE=FNAME)
        read (NFILE, '(A)', END=201)
        read (NFILE, *, END=202) T, LP0, HR0
        T=0D0
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

 202    write (*, '("Error: incomplete file: ",A)') FNAME
        close (NFILE)
        T=0D0
        LP_SWR=0D0
        HR_SWR=0D0
        return
 201    write (*, '("Error: skip file: ",A)') FNAME
        close (NFILE)
      end
