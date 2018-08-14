      subroutine SET_HE3PT()
        include 'he3_const.fh'
        real*8 TEMP, T1, TR

!        call STOP_SWEEP

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

