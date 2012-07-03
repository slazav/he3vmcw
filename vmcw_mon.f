C-- MONITOR ---- MONITORING THE SOLUTION ----------------------------
      subroutine MONITOR()
        include 'vmcw.fh'
        include 'par.fh'
        common /TIMEP/ T, TSTEP, TEND
        common /GEAR0/ DTUSED,NQUSED,NSTEP,NFE,NJE
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)

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
        write(*,'(7(A,F8.3,A))')
     *    ' T=',TMMS,' ms, ',
     *    'LP=', LP0+LP_SWR*T , ' cm, ',
     *    'HR=',  1D3*(HR0+HR_SWR*T) , ' mOe, ',
     *    'TF=',  1D6*(TF0+TF_SWR*T) , ' mks, ',
     *    'LF=',  1D-3*(LF0+LF_SWR*T) , ' kHz, ',
     *    'DF=',  (DF0+DF_SWR*T) , ' cm^2/s, ',
     *    'CPAR=',  (CPAR0+CPAR_SWR*T) , ' cm/s, '
        call WRITEMJ_DO()
   61   format(F7.1, 6(F12.8))
   69   format(F9.2, 1PE25.16)
      end
C-- WRITE_MJ --- WRITE SPINS & CURRENTS TO VMCW ------------------

      subroutine WRITEMJ_OPEN()
        include 'vmcw.fh'
        include 'par.fh'
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)
        integer FILES_MJ(NPTS), FILES_MJ0
        common /FILES/ FILES_MJ, FILES_MJ0
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
        FILES_MJ0=1000
        open(FILES_MJ0,FILE='mj_all.dat')
      end

      subroutine WRITEMJ_DO()
        include 'vmcw.fh'
        include 'par.fh'
        common /TIMEP/ T, TSTEP, TEND
        integer FILES_MJ(NPTS), FILES_MJ0
        common /FILES/ FILES_MJ, FILES_MJ0
        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),X(NPTS)

        DIFF=DF0+DF_SWR*T

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

            write(FILES_MJ0,103) X(I), T*1000D0,
     *        (LP0+LP_SWR*T)/CELL_LEN,
     *        USOL(1,I,1),USOL(2,I,1),USOL(3,I,1),
     *        USOL(4,I,1),USOL(5,I,1),USOL(6,I,1),USOL(7,I,1)

C            write(24,102) T*1000., X(I),
C     *        USOL(1,I,2),USOL(2,I,2),USOL(3,I,2),
C     *        USOL(4,I,2),USOL(5,I,2),USOL(6,I,2),USOL(7,I,2),
C     *        UFX,UFY,UFZ
          endif
        enddo
        write(FILES_MJ0,*)
        flush(FILES_MJ0)
        flush(FILES_MJ(I))
C       write(24,*)''
  101   format(F7.1 F10.6, 7(1PE15.6))
  102   format(F7.1 F10.6, 10(1PE15.6))
  103   format(10(1PE15.6))
      end

