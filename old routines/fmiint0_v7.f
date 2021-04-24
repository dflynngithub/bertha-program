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
C                        ∞                                             C
C      FMIINT0(L,λ,ζ) = ∫ r^2L+1 exp(-λ r^2) r*V_fermi(r) dr.          C
C                        0                                             C
C  THIS ROUTINE USES MY OWN RECIPE FOR THE ANNOYING INTEGRAL CLASSES,  C
C  FOR THE FIRST FEW L CASES AND ON THE CONDITION THAT A√λ < 0.05D0.   C
C -------------------------------------------------------------------- C
C  DFNOTE: THIS VERSION IS FINISHED (EXPANSION IN N/2A√λ) BUT TRUNCATED.C
C          IT'S NOT SITTING IN BERTHA BECAUSE I FOUND A NICER METHOD.  C
C**********************************************************************C
      INCLUDE 'param.h'
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12
C
      DIMENSION U(0:5),X(0:6),P(8)
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In FMIINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In FMIINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
C     FACTORS NEEDED FOR ALL PARAMETERS N
      U(0) = 1.0D0
      DO K=1,5
        U(K) = U(K-1)*A/C
      ENDDO
      
      X(0) = 1.0D0
      DO K=1,6
        X(K) = X(K-1)*EIJ*C*C
      ENDDO
C
      P(4) =-(PI**4)*7.0D0/720.0D0
      P(6) =-(PI**6)*31.0D0/30240.0D0
      P(8) =-(PI**8)*127.0D0/12096.0D2
C
C     FERMI NORMALISATION CONSTANT
      FNRM = 1.0D0 + PI*PI*U(2) - 6.0D0*U(3)*POLYLOG(3,-1.0D0/U(1))
      FNRM = 1.0D0/FNRM
C
C     INTEGRAL TYPE: INCOMPLETE GAMMA FUNCTION
      X1A = GAMLWR(2*L+3,X(1))*(3.0D0+PI*PI*U(2))/(2.0D0*DSQRT(X(1)))
      X1B =-GAMLWR(2*L+5,X(1))/(2.0D0*X(1)*DSQRT(X(1)))
      X1C = GAMUPR(2*L+2,X(1))*(1.0D0+PI*PI*U(2))
C
C     INTEGRAL TYPE: LOGARITHMIC INTEGRAL
      X2A =-6.0D0*GAMHLF(2*L)*U(3)*POLYLOG(3,-1.0D0/U(1))
C
C     TRANSFER DATA TO FMIINT0
      X1 = X1A+X1B+X1C
      X2 = X2A
      FMIINT0 = 0.5D0*FNRM*(X1+X2)/(EIJ**(L+1))
C
      TOL = 0.05D0
C
C     ALGORITHMS AVAILABLE FOR THE REMAINING SET OF INTEGRALS...
C
C     EXPAND CORRECT EXPRESSION FOR V3+V4 GIVEN LARGE X = N/(2*A*√λ)
      IF(A*DSQRT(EIJ).LT.TOL) THEN
C
        IF(L.EQ.0) THEN
C
          Y04 = 8.0D0 - 12.0D0*X(1) + 8.0D0*X(2) - 1.0D1*X(3)/3.0D0
     &        + X(4) - 7.0D0*X(5)/3.0D1 + 2.0D0*X(6)/45.0D0 
          Y06 =-72.0D0*X(1) + 16.0D1*X(2) - 14.0D1*X(3) + 72.0D0*X(4) 
     &        - 77.0D0*X(5)/3.0D0 + 104.0D0*X(6)/15.0D0
          Y08 = 96.0D1*X(2) - 28.0D2*X(3) + 3024.0D0*X(4)
     &        - 1848.0D0*X(5) + 2288.0D0*X(6)/3.0D0
C
        ELSEIF(L.EQ.1) THEN
C
          Y04 = 12.0D0 - 16.0D0*X(1) + 1.0D1*X(2) - 4.0D0*X(3)
     &        + 7.0D0*X(4)/6.0D0 - 4.0D0*X(5)/15.0D0
          Y06 = 72.0D0 - 32.0D1*X(1) + 42.0D1*X(2) - 288.0D0*X(3) 
     &        + 385.0D0*X(4)/3.0D0 - 41.6D0*X(5)
          Y08 =-192.0D1*X(1) + 84.0D2*X(2) - 12096.0D0*X(3)
     &        + 924.0D1*X(4) - 4576.0D0*X(5)
C
        ELSEIF(L.EQ.2) THEN
C
          Y04 = 16.0D0 - 2.0D1*X(1) + 12.0D0*X(2) - 14.0D0*X(3)/3.0D0
     &        - 2.03125D0*X(4)
          Y06 = 32.0D1 - 84.0D1*X(1) + 864.0D0*X(2)
     &        - 154.0D1*X(3)/3.0D0 - 316.875D0*X(4)
          Y08 = 192.0D1 - 168.0D2*X(1) + 36288.0D0*X(2) - 3696.0D1*X(3)
     &        - 34856.25D0*X(4)
C
        ELSEIF(L.EQ.3) THEN
C
          Y04 = 2.0D1 - 24.0D0*X(1) - 117.21875D0*X(2)
     &        - 2639.0D0*X(3)/192.0D0
          Y06 = 84.0D1 - 1728.0D0*X(1) - 12894.0625D0*X(2)
     &        - 2144.1875D0*X(3)
          Y08 = 168.0D2 - 72576.0D0*X(1) - 928372.5D0*X(2)
     &        - 235860.625D0*X(3)
C
        ELSEIF(L.EQ.4) THEN
C
          Y04 =-3584.515625D0 - 1471.40625D0*X(1) - 100.078125D0*X(2)
          Y06 =-258085.125D0 - 161854.6875D0*X(1) - 15612.1875D0*X(2)
          Y08 =-10839575.25D0 - 11653537.5D0*X(1) - 1717340.625D0*X(2)
C
        ELSEIF(L.EQ.5) THEN
C
          Y04 =-15816.6640625D0 - 7954.33203125D0*X(1)
          Y06 =-1739833.046875D0 - 1240875.7968750*X(1)
          Y08 =-125267979.375D0 - 136496337.65625D0*X(1)
C
        ELSEIF(L.EQ.6) THEN
C
          Y04 =-68526.642578125D0
          Y06 =-10690156.2421875D0
          Y08 =-1175917186.640625D0
C
C       POLYNOMIAL EXPANSION VANISHES FOR ANY L > 6.
        ELSE
C
          Y04 = 0.0D0
          Y06 = 0.0D0
          Y08 = 0.0D0
C
        ENDIF
C
C       COMBINE EXPRESSIONS
        YP   = Y04*P(4)*U(1) + Y06*P(6)*U(3) + Y08*P(8)*U(5)
C     
C       INCLUDE MULTIPLICATIVE FACTOR
        Y34 = 3.0D0*U(3)*(C**(2*L+2))*YP
C
C     NUMERICAL QUADRATURE FOR V3+V4 (WORKS FOR ANY SET OF PARAMETERS)
      ELSE
C
C       INTEGRATION GRID DETAILS
        NLIN = 100
        NEXP = 1000
        RMAX = 10.0D0
        ELIN = 1.0D0/DFLOAT(NLIN)
        EEXP = DLOG(RMAX/C)/DFLOAT(NEXP)
C
C       NUMERICALLY INTEGRATE REMAINING TERMS FROM 0 TO C
        X0C = 0.0D0
        DO N=0,NLIN
          ZN  = ELIN*DFLOAT(N)
          Y1  = ZN**(2*L+1)
          Y2  = DEXP(-X(1)*ZN*ZN)
          Y3  =-ZN*POLYLOG(2,(ZN-1.0D0)/U(1))
          Y4  = 2.0D0*U(1)*POLYLOG(3,(ZN-1.0D0)/U(1))
          X0C = X0C + EXTINT11(Y1*Y2*(Y3+Y4),N,NLIN)
        ENDDO
        X0C = 15.0D0*U(2)*(C**(2*L+2))*ELIN*X0C/299376.0D0
C
C       NUMERICALLY INTEGRATE REMAINING TERMS FROM C TO INFINITY (EXP.)
        XCI = 0.0D0
        DO N=0,NEXP
          TN  = EEXP*DFLOAT(N)
          ZN  = DEXP(TN)
          Y1  = ZN**(2*L+2)
          Y2  = DEXP(-X(1)*ZN*ZN)
          Y3  = ZN*POLYLOG(2,(1.0D0-ZN)/U(1))
          Y4  = 2.0D0*U(1)*POLYLOG(3,(1.0D0-ZN)/U(1))
          XCI = XCI + EXTINT11(Y1*Y2*(Y3+Y4),N,NEXP)
        ENDDO
        XCI = 15.0D0*U(2)*(C**(2*L+2))*EEXP*XCI/299376.0D0
        
        Y34 = X0C+XCI
C
      ENDIF
C      
      FMIINT0 = FMIINT0 + FNRM*Y34
C
      RETURN
      END

