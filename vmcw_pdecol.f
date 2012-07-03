      subroutine pdecol_init(T)
        include 'par.fh'
        real*8 T

        common /PDECOL_DATA/ INDEX,MF,SCTCH,WORK,IWORK,T0,DT
        integer INDEX, IWORK(IDIMIWORK)
        real*8 SCTCH(KORD*(NDERV+1)),WORK(IDIMWORK),T0,DT

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

      subroutine pdecol_run(T)
        include 'vmcw.fh'
        include 'par.fh'
        real*8 T

        common /PDECOL_DATA/ INDEX,MF,SCTCH,WORK,IWORK,T0,DT
        integer INDEX,IWORK(IDIMIWORK)
        real*8 SCTCH(KORD*(NDERV+1)),WORK(IDIMWORK),T0,DT

        common /ARRAYS/ USOL(NPDE,NPTS,NDERV),XSOL(NPTS)
        real*8 USOL,XSOL

        call PDECOL(T0,T,DT,XSOL,ACC,NINT,KORD,NCC,NPDE,MF,
     +              INDEX,WORK,IWORK)
        if(INDEX.ne.0) THEN
          write(*,*) 'INTEGRATION FAILED; INDEX=', INDEX
          stop
        endif

        call VALUES(XSOL,USOL,SCTCH,NPDE,NPTS,NPTS,2,WORK)
      end
