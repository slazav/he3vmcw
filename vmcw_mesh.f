C------- AEROGEL STEP -----------------------------------------
C       Aerogel density function. Returns 1 in the central
C       part of the cell with fermi steps to 0 on edges.
C         X        -- coord, cm
C         D        -- derivative order (0|1)
C 
C         CELL_LEN -- cell length, cm
C         AER      -- if >0 then do step
C         AER_LEN  -- aerogel length / cell length
C         AER_CNT  -- center of aerogel area / cell length
C         AER_TRW  -- transition width / cell length
      real*8 function AER_STEP(X,D)
        include 'vmcw_pars.fh'
        integer D
        real*8  X
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
      subroutine SET_MESH(X, N)
        include 'vmcw_pars.fh'
C       Set mesh according with AER_STEP function
        real*8 DELTA,DX,X, AER_STEP
        dimension X(N)
        X(1)=-CELL_LEN/2D0
C       start with homogenious mesh with DX intervals
        DX=CELL_LEN/DFLOAT(N-1)
        do I=1,100
C         build mesh with scaled intervals
          do J=1,N-1
            X(J+1)=X(J) +
     +        DX/(1.0D0+XMESH_K*ABS(AER_STEP(X(J),1)))
          enddo
C         scale the whole mesh to fit CELL_LEN
          DELTA=CELL_LEN - (X(N)-X(1))
          DX = DX + DELTA/DFLOAT(N+1)
          if (ABS(DELTA).LT.XMESH_ACC) return
        enddo
        write(*,*) 'warning: low mesh accuracy: ', ABS(DELTA)
      end

      subroutine SAVE_MESH(X, N)
        integer fid/54/, N
        real*8 X, AER_STEP
        dimension X(N)
        open(fid,FILE='mesh.dat')
        write(fid,*) '#  I    X(I) STEP(X) STEP''(X)'
        do J=1,N
          write(fid,'(I4," ",F7.5," ",F7.5," ",e12.5e2)')
     +     J, X(J), AER_STEP(X(J),0), AER_STEP(X(J),1)
        enddo
        close(fid)
      end
