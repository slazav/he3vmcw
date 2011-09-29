CC    T_c [mK] vs P [bars] -- TCF(P)
C     Greywall. Phys. Rev.B v.33 #11 p.7520
      function TCF(P)
        real*8 TCF,MCF,P
        TCF=   .92938375D0
     .       + .13867188D0*P
     .       - .69302185D-2*P**2
     .       + .25685169D-3*P**3
     .       - .57248644D-5*P**4
     .       + .53010918D-7*P**5
        if (P.GT.MCF(TCF)) then
          TCF=0D0
        endif
        return
      end

CC    T_ab [mK] vs P [bars] -- TABF(P)
C     Greywall. Phys. Rev.B v.33 #11 p.7520
      function TABF(P)
        real*8 TABF,MCF,P,PR,TCF
        common /HE3DATA/ PCP,TPCP,PA,ANA,PI,HC,R,AKB,GAM,AM3
        real*8 PCP,TPCP,PA,ANA,PI,HC,R,AKB,GAM,AM3
        PR=P-PCP
        if (PR.LT.0D0) then
          TABF=TCF(P)
        else
          TABF= TPCP
     .       - .10322623D-1*PR
     .       - .53633181D-2*PR**2
     .       + .83437032D-3*PR**3
     .       - .61709783D-4*PR**4
     .       + .17038992D-5*PR**5
        endif
        if (P.GT.MCF(TABF)) then
          TABF=0D0
        endif
        return
      end

CC    MOLAR VOLUME V[cm**3/mole] vs P [bars] -- MVF(P)
C     Greywall. Phys. Rev.B v.33 #11 p.7520 ref. 27
C     That is Wheatley Rev.Mod.Phys. 47,415(1975).
      function MVF(P)
        real*8 MVF,P
        MVF=   36.837231D0
     .       - 1.1803474D0*P
     .       + 0.0834214D0*P**2
     .       - 0.3885962D-2*P**3
     .       + 0.94759780D-4*P**4
     .       - 0.91253577D-6*P**5
        return
      end

CC    MELTING PRESSURE P [bars] vs T [mK] -- MCF(T)
C     Greywall. Phys. Rev.B v.33 #11 p.7520
      function MCF(T)
        real*8 MCF,T
        common /HE3DATA/ PCP,TPCP,PA,ANA,PI,HC,R,AKB,GAM,AM3
        real*8 PCP,TPCP,PA,ANA,PI,HC,R,AKB,GAM,AM3
        MCF=PA-.19652970D-1*T**(-3)-.61880268D-1*T**(-2)-
     .      .78803055D-1*T**(-1)+.13050600D0-.43519381D-1*T+
     .      .13752791D-3*T**2-.17180436D-6*T**3-
     .      .22093906D-9*T**4+.85450245D-12*T**5
        return
      end

CC    R-Gas constant GAMMA=C/RT [1/(K*mol)] vs P [bar] -- GAMMAF(P)
C     Greywall. Phys. Rev.B v.33 #11 p.7520
      function GAMMAF(P)
        real*8 GAMMAF,P
        GAMMAF =  .27840464D+1
     .          + .69575243D-1*P
     .          - .14738303D-2*P**2
     .          + .46153498D-4*P**3
     .          - .53785385D-6*P**4
        return
      end

CC    LEGGETT FREQ^2, Lf^2 [Hz^2] vs P [bar], T/Tc -- LF2F(P,T)
C     Ahonen (18.7,21.1,25.4,29.0,32 bars), Osheroff(MP).
C     Ahonen et.al.  JLTP. v.25 p.421(1976)
C     Osheroff       PRL v.34. p.190
C     From 0 to 18.7 bar data are not reliable.
C     See file YOM3.
C     Least squares 2-D fitting. From fil:: YOM3
C     F= vs. T= & P=
C     Polinmms of the orders : 4,  1
      function LF2F(P,T)
        real*8 LF2F,P,T,TABF,TCF
        real*8 WORK(5),A(10)
        real*8 XMIN,XMAX,YMIN,YMAX,F
        DATA XMIN/0.266D0/,XMAX/1.000D0/
        DATA YMIN/0.000D0/,YMAX/34.36D0/
        DATA K/4/,L/1/,NA/10/,NWORK/5/
        DATA A/  1.3169660D+11,  5.0310062D+10,
     .          -6.6371420D+10, -2.0536020D+10,
     .          -5.2457041D+09, -5.1179683D+09,
     .           5.8596372D+09,  3.1320310D+08,
     .          -6.9788020D+07,  2.0270360D+08/
        if (T.LT.TABF(P)/TCF(P)) THEN
          IFAIL=1
          call E02CBE(1,1,K,L,T,XMIN,XMAX,P,YMIN,YMAX,F,
     ,                A,NA,WORK,NWORK,IFAIL)
          LF2F=F
          if (IFAIL.EQ.2) THEN
            print *,'P out of range:',YMIN,YMAX
            LF2F=0D0
          elseif (IFAIL.EQ.3) THEN
            print *,'T out of range:',XMIN,XMAX
            LF2F=0D0
          elseif (IFAIL.NE.0) THEN
            print *,'Error:',IFAIL
            LF2F=0D0
          endif
          if (P.LT.18.7D0) print *,'Data not reliable.'
        else
          print *,'No data for A-phase.'
          LF2F=0D0
        endif
        return
      end

C--   Attention density of state must be multiplyied by ANA later.
      function DNDEF(P)
        real*8 DNDEF,P,GAMMAF
        common /HE3DATA/ PCP,TPCP,PA,ANA,PI,HC,R,AKB,GAM,AM3
        real*8 PCP,TPCP,PA,ANA,PI,HC,R,AKB,GAM,AM3
        DNDEF=3D0*GAMMAF(P)/AKB/PI**2
        return
      end

CC    EFFECTIVE MASS [g] vs P [bar] -- MAF(P)
      function MAF(P)
        real*8 MAF,P,MVF,DNDEF,PF,PFF
        PF=PFF(P)
C       MAF=GAMMAF(P)*R*HC*MVF(P)*(HC/PF)*3
C       print *,GAMMAF(P),GAMMAF(P)*R*HC,R
        MAF=DNDEF(P)*PF/3D0*PF
C       MAF=DNDEF(P)*PF/3./MVF(P)*PF
        return
      end

CC    FERMI MOMENTUM [sgs] vs P [bar] -- PFF(P)
      function PFF(P)
        common /HE3DATA/ PCP,TPCP,PA,ANA,PI,HC,R,AKB,GAM,AM3
        real*8 PCP,TPCP,PA,ANA,PI,HC,R,AKB,GAM,AM3
        real*8 PFF,P,MVF
        PFF=HC*(3D0*PI**2*ANA/MVF(P))**.3333333D0
        return
      end

CC    FERMI VELOSITY [cm/s] vs P [bar] -- VFF(P)
      function VFF(P)
        real*8 VFF,P,PFF,MAF
        VFF=PFF(P)/MAF(P)
        return
      end

CC    Osheroff's spin wave velocity. S [cm/s] vs P [bar], T [mK] -- SF(P,T)
C     Osheroff's spin wave velocity if recalculated to arbitrary pressure
C     by taking value of Osheroff for melting pressure (1100 cm/sek)
C     and assuming S is proportional to Fermi velocity.
      function SF(P,T)
        real*8 SF,P,T,VFF,TCF
        SF=1100D0/VFF(34.3D0)*VFF(P)*SQRT(1D0-T/TCF(P))
      end

CC    Parallel Fomin spin wave velocity. Cpar [cm/c] vs P [bar], T [mK] -- CPARF(P,T)
      function CPARF(P,T)
        real*8 CPARF,P,T,SF
        CPARF=SF(P,T)*SQRT(2.0D0)
      end
CC    Perp Fomin spin wave velocity. Cper [cm/c] vs P [bar], T [mK] -- CPERF(P,T)
      function CPERF(P,T)
        real*8 CPERF,P,T,SF
        CPERF=SF(P,T)*SQRT(1.5D0)
      end

CC    SPIN DIFFUSION COEFF, SUPERFLUID D [cm**2/sec] vs P [bar], T[mk] -- DSUPF(P,T)
C--   Our data. Was measured at pressures .6,11,14.8,20,29 bar
C     at T .45-.67,  .45-.8  .45-.81  .45-.94  .45-.67 Tc respectivly.
C     Least squares 2-D fitting. From fil:: MUHLR1
C     D/C= vs. T= & P=
C     Polinoms of the orders : 4,  4
      function DSUPF(P,T)
        real*8 DSUPF,P,T,CPARF,TCF,TR
        real*8 WORK(5),A(25)
        real*8 XMIN,XMAX,YMIN,YMAX,F
        DATA XMIN/0.440D0/,XMAX/0.9400001D0/
        DATA YMIN/0.000D0/,YMAX/29.00000D0/
        DATA K/4/, L/4/, NA/25/, NWORK/5/
        DATA A/  6.6012963D-03, -2.2342042D-03,  4.2030680D-04
     ,  ,  1.6044091D-04, -8.2734930D-04,  1.7867290D-03, -6.6435340D-04
     ,  , -1.8627312D-04,  4.3953900D-04, -6.3450170D-04, -4.1724560D-04
     ,  ,  3.9544643D-04, -3.2521144D-04,  1.7599921D-04, -1.9650970D-04
     ,  ,  1.2244173D-04,  3.2786340D-05,  6.5047061D-08, -3.2710520D-05
     ,  , -6.1170401D-05,  3.7204154D-05,  9.9630890D-06,  1.2296370D-08
     ,  , -9.9427400D-06, -1.8589210D-05/
        TR=T/TCF(P)
        IFAIL=1
        call E02CBE(1,1,K,L,TR,XMIN,XMAX,P,YMIN,YMAX,F,
     ,              A,NA,WORK,NWORK,IFAIL)
        if (IFAIL.EQ.2) THEN
          print *,'Y out of range.'
        else if (IFAIL.EQ.3) THEN
          print *,'X out of range.'
        else if (IFAIL.NE.0) THEN
          print *,'Error:',IFAIL
        end if
        DSUPF=F*CPARF(P,T)**(.6666666666666D0)
      end
C--

CC    SPIN DIFFUSION COEFF, NORMAL D [cm**2/sec] vs P [bar], T[mk] -- DNORF(P,T)
C--   Brewer data for D*T**2
C--   JLTP v.56 #5/6 p.617
C     Least squares fitting. From file: YDIFD
C     DT2=                  vs.  P=
C     Polinom of the order : 4
C     Residual: 0.000
      function DNORF(P,T)
        real*8 DNORF,P,T,DT2F, XCAP,F
        DATA M1/4/
        real*8 XMAX, XMIN
        DATA XMIN/0.0D0/,XMAX/27.75999D0/
        real*8 A(4)
        DATA A/1.0603673D-06, -5.4074212D-07,
     .         2.2418920D-07, -6.4375400D-08/
        IFAIL=1
        XCAP=((P-XMIN)-(XMAX-P))/(XMAX-XMIN)
        call E02AEE(M1,A,XCAP,F,IFAIL)
        if (IFAIL.NE.0) print *,'Error in E02AEE :',IFAIL
C       .89 accounts fo Grewall scale.
        DT2F=(F*0.89D0**2)
        DNORF=DT2F/(T*1D-3)**2
        return
      end

CC    SPIN DIFFUSION COEFF D [cm**2/sec] vs P [bar], T[mk] -- DF(P,T)
      function DF(P,T)
        real*8 DF, P,T, DNORF,DSUPF,TCF
        if (T.GT.TCF(P)) THEN
          DF=DNORF(P,T)
        else
          DF=DSUPF(P,T)
C          print *,'Superflow region. Check if data out range.'
        end if
        return
      end

CC    SUSEPTIBILITY [sgs] vs P [bar], T [mK] -- HIF(P,T)
      function HIF(P,T)
        real*8 HIF,P,T,MVF,Y,YOSHIDF,TTC,TCF,DNDEF,Z0F,Z0
        common /HE3DATA/ PCP,TPCP,PA,ANA,PI,HC,R,AKB,GAM,AM3
        real*8 PCP,TPCP,PA,ANA,PI,HC,R,AKB,GAM,AM3
        HIF=0.25D0*GAM**2*HC*DNDEF(P)*HC*ANA/(1D0+Z0F(P)/4D0)/MVF(P)
        TTC=T/TCF(P)
        if (TTC.LT.1D0) then
          Z0=Z0F(P)
          Y=YOSHIDF(TTC)
C         print *,Z,Y,((1.+Z0/4.)*(2.+Y)/3.)/(1.+Z0/4.*(2.+Y)/3.)
          HIF=HIF*((1D0+Z0/4D0)*(2D0+Y)/3D0)/(1D0+Z0/4D0*(2D0+Y)/3D0)
        end if
      end


CC    YOSIDA FUNCTION vs T/Tc -- YOSHIDF(TTC)
C     Sourse B-phase notebook.
C     Least squares fitting. From file: YOSHID
C     Polinom of the order : 5
C     Residual: 0.000
      function YOSHIDF(TTC)
        real*8 YOSHIDF,TTC
        real*8 A(5)
        DATA M1/5/
        real*8 XMIN,XMAX,XCAP
        DATA XMIN/9.9999994D-02/,XMAX/1.000000D0/
        DATA A/.7349111D0, .5123515D0, .1371038D0,
     ,         -1.4855450D-02, -4.5979050D-03/
        if (TTC.GE.0.1D0) then
          IFAIL=1
          XCAP=((TTC-XMIN)-(XMAX-TTC))/(XMAX-XMIN)
          call E02AEE(M1,A,XCAP,YOSHIDF,IFAIL)
          if (IFAIL.NE.0)print *,'Error in E02AEE :',IFAIL
        else
          YOSHIDF=0D0
        end if
        return
      end

CC    Z0 vs P [bar] -- Z0F(P)
C--   Sourse B-phase notebook.
C     Least squares fitting. From file: YZ0
C     Polinom of the order : 5
C     Residual: 0.000
      function Z0F(P)
        real*8 Z0F,P
        real*8 A(5)
        DATA M1/ 5/
        real*8 XMIN,XMAX,XCAP
        DATA XMIN/0.000000D0/,XMAX/34.36000D0/
        DATA A/-5.762472D0, -0.1136529D0, 5.5511940D-02,
     *       -1.7914600D-02, 4.0055060D-03/
        IFAIL=1
        XCAP=((P-XMIN)-(XMAX-P))/(XMAX-XMIN)
        call E02AEE(M1,A,XCAP,Z0F,IFAIL)
        if (IFAIL.NE.0) print *,'Error in E02AEE :',IFAIL
      end

CC    F0A vs P [bar] -- F0AF(P)
C-    Sourse - Greywall, Phys.Rev.  v.27  5  (1983)
C     Least squares fitting. From file: YF0A
C     Polinom of the order : 7
C     Residual: 0.000
      function F0AF(P)
        real*8 F0AF,P
        real*8 A(7)
        DATA M1/7/
        real*8 XMIN,XMAX,XCAP
        DATA XMIN/0.000000D0/,XMAX/34.36000D0/
        DATA A/-1.489332D0, -2.3159460D-02,  1.3571171D-02,
     .         -4.2908173D-03,  1.4413130D-03, -1.1601811D-03,
     .          9.9658221D-04/
        IFAIL=1
        XCAP=((P-XMIN)-(XMAX-P))/(XMAX-XMIN)
        call E02AEE(M1,A,XCAP,F0AF,IFAIL)
        if (IFAIL.NE.0)print *,'Error in E02AEE :',IFAIL
      end


      include 'libs/E02AEE.FOR'
      include 'libs/E02CBE.FOR'
      include 'libs/M01AGE.FOR'
      include 'libs/P01AAE.FOR'
      include 'libs/X02AAE.FOR'
      include 'libs/X04AAE.FOR'
