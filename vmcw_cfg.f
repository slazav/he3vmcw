      subroutine read_cfg(fname)
        include 'vmcw.fh'
        character*(*) fname

c       local vars
        integer fid/54/
        real*8    CFG_VAL
        character CFG_KEY(64)

        open(fid,FILE=fname)
   11   read(fid,*,END=12) CFG_KEY, CFG_VAL
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
   12   close(fid)

        LP0=0D0
        HR0=1D-3
        LP_SWR=0D0
        HR_SWR=0D0
        WRITEMJ_XSTEP=0.1D0
      end
