      FUNCTION FMIINT0(L,EIJ,A,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       FFFFFFFF MM       MM IIII IIII NN    NN TTTTTTTT 000000        C
C       FF       MMM     MMM  II   II  NNN   NN    TT   00   000       C
C       FF       MMMM   MMMM  II   II  NNNN  NN    TT   00  0000       C
C       FFFFFF   MM MM MM MM  II   II  NN NN NN    TT   00 00 00       C
C       FF       MM  MMM  MM  II   II  NN  NNNN    TT   0000  00       C
C       FF       MM   M   MM  II   II  NN   NNN    TT   000   00       C
C       FF       MM       MM IIII IIII NN    NN    TT    000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  FMIINT0 CALCULATES AN AUXILLIARY FERMI POTENTIAL INTEGRAL OVER      C
C  A PAIR OF ATOMIC BASIS FUNCTIONS, WHERE L IS A POSITIVE INTEGER.    C
C                                                                      C
C    FMIINT0(L,λ,ζ) = ∫^{∞}_{0} r^2L+1 exp(-λ r^2) r*V_fermi(r) dr.    C
C                                                                      C
C -------------------------------------------------------------------- C
C  THIS IS MY OWN SERIES EXPANSION METHOD -- NO EFFECT AFTER L=2.      C
C**********************************************************************C
      INCLUDE 'param.h'
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In FMIINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In FMIINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
C     FACTORS NEEDED FOR ALL PARAMETERS N
      U   = A/C
      U3  = U*U*U
      X   = EIJ*C*C
      X12 = DSQRT(X)
      C2L = C*C
      C2L = C2L**(L+1)
C
C     FERMI NORMALISATION CONSTANT
      FNRM = 1.0D0 + PI*PI*U*U - 6.0D0*U3*POLYLOG(3,-1.0D0/U)
      FNRM = 1.0D0/FNRM
C
C     INTEGRAL TYPE: INCOMPLETE GAMMA FUNCTION
      X1A = GAMLWR(2*L+3,X)*(3.0D0+PI*PI*U*U)/(2.0D0*X12)
      X1B =-GAMLWR(2*L+5,X)/(2.0D0*X12*X12*X12)
      X1C = GAMUPR(2*L+2,X)*(1.0D0+PI*PI*U*U)
C
C     INTEGRAL TYPE: LOGARITHMIC INTEGRAL
      X2A =-6.0D0*GAMHLF(2*L)*U3*POLYLOG(3,-1.0D0/U)
C
C     TRANSFER DATA TO FMIINT0
      X1 = X1A+X1B+X1C
      X2 = X2A
      FMIINT0 = 0.5D0*FNRM*(X1+X2)/(EIJ**(L+1))
C
C     ALGORITHMS AVAILABLE FOR THE REMAINING SET OF INTEGRALS
      U02 = U*U
      U04 = U*U*U02
      U06 = U*U*U04
      U08 = U*U*U06
      U10 = U*U*U08
      U12 = U*U*U10
      U14 = U*U*U12
      U16 = U*U*U14
      U18 = U*U*U16

      X0  = 1.0D0
      X1  = X*X0
      X2  = X*X1
      X3  = X*X2
      X4  = X*X3
      X5  = X*X4
      X6  = X*X5
      X7  = X*X6

      P04 = POLYLOG( 4,0.0D0)
      P06 = POLYLOG( 6,0.0D0)
      P08 = POLYLOG( 8,0.0D0)
      P10 = POLYLOG(10,0.0D0)
      P12 = POLYLOG(12,0.0D0)
      P14 = POLYLOG(14,0.0D0)
      P16 = POLYLOG(16,0.0D0)
      P18 = POLYLOG(18,0.0D0)
      
      S05 = POLYLOG( 5,-1.0D0/U)
      S07 = POLYLOG( 7,-1.0D0/U)
      S09 = POLYLOG( 9,-1.0D0/U)
      S11 = POLYLOG(11,-1.0D0/U)
      S13 = POLYLOG(13,-1.0D0/U)
      S15 = POLYLOG(15,-1.0D0/U)
      S17 = POLYLOG(17,-1.0D0/U)
      S19 = POLYLOG(19,-1.0D0/U)
      S21 = POLYLOG(21,-1.0D0/U)
      S23 = POLYLOG(23,-1.0D0/U)
      
      Y04 = 0.0D0
      Y06 = 0.0D0
      Y08 = 0.0D0
      Y10 = 0.0D0
      Y12 = 0.0D0
      Y14 = 0.0D0
      Y16 = 0.0D0
      Y18 = 0.0D0

      E05 = 0.0D0
      E07 = 0.0D0
      E09 = 0.0D0
      E11 = 0.0D0
      E13 = 0.0D0
      E15 = 0.0D0
      E17 = 0.0D0
      E19 = 0.0D0
      E21 = 0.0D0
      E23 = 0.0D0

      Y34 = 0.0D0
      IF(L.EQ.0.AND.U*DSQRT(X).LT.0.1D0) THEN
C
        Y04 = 8.0D0 - 12.0D0*X1 + 8.0D0*X2 - 10.0D0*X3/3.0D0
     &      + X4 - 7.0D0*X5/30.0D0 + 2.0D0*X6/45.0D0 - X7/140.0D0
        Y06 =-72.0D0 + 160.0D0*X1 - 140.0D0*X2 + 72.0D0*X3 
     &      - 77.0D0*X4/3.0D0 + 104.0D0*X5/15.0D0 - 3.0D0*X6/2.0D0
        Y08 = 960.0D0 - 2800.0D0*X1 + 3024.0D0*X2 - 1848.0D0*X3 
     &      + 2288.0D0*X4/3.0D0 - 234.0D0*X5
        Y10 =-16800.0D0 + 60480.0D0*X1 - 77616.0D0*X2 + 54912.0D0*X3 
     &      - 25740.0D0*X4
        Y12 = 362880.0D0 - 1552320.0D0*X1 + 2306304.0D0*X2 
     &      - 1853280.0D0*X3
        Y14 =-9313920.0D0 + 46126080.0D0*X1 - 77837760.0D0*X2
        Y16 = 276756480.0D0 - 1556755200.0D0*X1
        Y18 =-9340531200.0D0
        
        YP   = U02*X0*Y04*P04 + U04*X1*Y06*P06 + U06*X2*Y08*P08
     &       + U08*X3*Y10*P10 + U10*X4*Y12*P12 + U12*X5*Y14*P14
     &       + U14*X6*Y16*P16 + U16*X7*Y18*P18

        E05 = 1.0D0
        E07 =-9.0D0*X1*U02
        E09 = 120.0D0*X2*U04
        E11 =-2100.0D0*X3*U06
        E13 = 45360.0D0*X4*U08
        E15 =-1164240.0D0*X5*U10
        E17 = 34594560.0D0*X6*U12
        E19 =-1167566400.0D0*X7*U14
        
        YE  = E05*S05 + E07*S07 + E09*S09 + E11*S11
     &      + E13*S13 + E15*S15 + E17*S17 + E19*S19

        Y34  = 3.0D0*A*A*(YP + 4.0D0*U02*U*YE)
C
      ELSEIF(L.EQ.1.AND.U*DSQRT(X).LT.0.1D0) THEN
C
        Y04 = 12.0D0 - 16.0D0*X1 + 10.0D0*X2 - 4.0D0*X3
     &      + 7.0D0*X4/6.0D0 - 4.0D0*X5/15.0D0 + X6/20.0D0 - X7/126.0D0
        Y06 = 72.0D0 - 320.0D0*X1 + 420.0D0*X2 - 288.0D0*X3 
     &      + 385.0D0*X4/3.0D0 - 208.0D0*X5/5.0D0 + 21.0D0*X6/2.0D0
        Y08 =-1920.0D0 + 8400.0D0*X1 - 12096.0D0*X2 + 9240.0D0*X3 
     &      - 4576.0D0*X4 + 1638.0D0*X5
        Y10 = 50400.0D0 - 241920.0D0*X1 + 388080.0D0*X2 - 329472.0D0*X3 
     &      + 180180.0D0*X4
        Y12 =-1451520.0D0 + 7761600.0D0*X1 - 13837824.0D0*X2 
     &      + 12972960.0D0*X3
        Y14 = 46569600.0D0 - 276756480.0D0*X1 + 544864320.0D0*X2
        Y16 =-1660538880.0D0 + 10897286400.0D0*X1
        Y18 = 65383718400.0D0
        
        YP   = U02*X0*Y04*P04 + U04*X0*Y06*P06 + U06*X1*Y08*P08
     &       + U08*X2*Y10*P10 + U10*X3*Y12*P12 + U12*X4*Y14*P14
     &       + U14*X5*Y16*P16 + U16*X6*Y18*P18

        E07 = 3.0D0
        E09 =-80.0D0*X1*U02
        E11 = 2100.0D0*X2*U04
        E13 =-60480.0D0*X3*U06
        E15 = 1940400.0D0*X4*U08
        E17 =-69189120.0D0*X5*U10
        E19 = 2724321600.0D0*X6*U12
        E21 =-117621505000.0D0*X7*U14
        
        YE  = E07*S07 + E09*S09 + E11*S11 + E13*S13
     &      + E15*S15 + E17*S17 + E19*S19 + E21*S21

        Y34  = 3.0D0*U*U*C2L*(YP + 12.0D0*U04*U*YE)
C
      ELSEIF(L.EQ.2.AND.U*DSQRT(X).LT.0.1D0) THEN
C
        Y04 = 16.0D0 - 20.0D0*X1 + 12.0D0*X2 - 14.0D0*X3/3.0D0
     &      + 3.0D0*X4/10.0D0 - 3.0D0*X5/10.0D0 + X6/18.0D0
     &      - 11.0D0*X7/1260.0D0
        Y06 = 320.0D0 - 840.0D0*X1 + 864.0D0*X2 - 1540.0D0*X3/3.0D0
     &      + 208.0D0*X4 - 63.0D0*X5 + 136.0D0*X6/9.0D0
        Y08 = 1920.0D0 - 16800.0D0*X1 + 36288.0D0*X2 - 36960.0D0*X3 
     &      + 22880.0D0*X4 - 9828.0D0*X5
        Y10 =-100800.0D0 + 825760.0D0*X1 - 1552320.0D0*X2
     &      + 1647360.0D0*X3 - 1081080.0D0*X4
        Y12 = 4354560.0D0 - 31046400.0D0*X1 + 69189120.0D0*X2 
     &      - 77837760.0D0*X3
        Y14 =-186278400.0D0 + 1383782400.0D0*X1 - 3269185920.0D0*X2
        Y16 = 8302694400.0D0 - 65383708400.0D0*X1
        Y18 =-392302310400.0D0
        
        YP   = U02*X0*Y04*P04 + U04*X0*Y06*P06 + U06*X0*Y08*P08
     &       + U08*X1*Y10*P10 + U10*X2*Y12*P12 + U12*X3*Y14*P14
     &       + U14*X4*Y16*P16 + U16*X5*Y18*P18

        E09 = 2.0D0
        E11 =-105.0D0*X1*U02
        E13 = 4536.0D0*X2*U04
        E15 =-194040.0D0*X3*U06
        E17 = 8648640.0D0*X4*U08
        E19 =-408648240.0D0*X5*U10
        E21 = 20583763200.0D0*X6*U12
        E23 =-1106230245120.0D0*X7*U14
        
        YE  = E09*S09 + E11*S11 + E13*S13 + E15*S15
     &      + E17*S17 + E19*S19 + E21*S21 + E23*S23

        Y34  = 3.0D0*U*U*C2L*C2L*(YP + 480.0D0*U06*U*YE)
        y34  = 0.0d0
C
      ELSE
        Y34 = 0.0D0
      ENDIF

C      
      FMIINT0 = FMIINT0 + FNRM*Y34
C
      RETURN
      END
