      subroutine pdecol_init(T)
        include 'vmcw_pdecol.fh'
        real*8 T

        INDEX=1  ! TYPE OF CALL (FIRST CALL)
        MF=22
        IWORK(1)=IDIMWORK
        IWORK(2)=IDIMIWORK
        do I=1,IDIMWORK
          WORK(I)=0.0D0
        enddo
        T0=T       ! STARTING TIME FOR PDECOL
        DT=1.D-10  ! INITIAL STEP SIZE IN T
      end

      subroutine pdecol_run(T,USOL,XSOL)
        include 'vmcw_pdecol.fh'
        real*8 T,USOL(NPDE,NPTS,NDERV),XSOL(NPTS)

        call PDECOL(T0,T,DT,XSOL,ACC,NINT,KORD,NCC,NPDE,MF,
     +              INDEX,WORK,IWORK)
        if(INDEX.ne.0) THEN
          write(*,*) 'INTEGRATION FAILED; INDEX=', INDEX
          stop
        endif

        call VALUES(XSOL,USOL,SCTCH,NPDE,NPTS,NPTS,2,WORK)
      end
