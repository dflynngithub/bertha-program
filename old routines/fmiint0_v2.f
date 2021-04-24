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
C DFNOTE: THE ANNOYING TERMS WERE CALCULATED BY GAUSSIAN QUADRATURE    C
C         USING AN 11-POINT NEWTON-COTES SCHEME, WITH UNIFORM GRID     C
C         BETWEEN R=0 AND R=C, AND EXPONENTIAL BETWEEN R=C AND R=INF.  C
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
      write(*,*) a
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
      IF(X.LT.1.0D0) THEN
C     EXPAND THE GAUSSIAN AS A TAYLOR SERIES
C
C       TRUNCATION DETAILS
        NGSS = NTR
C
C       INTEGRATE REMAINING TERMS FROM 0 TO C
        X0C = 0.0D0
        PHK = 1.0D0
        XUP = 1.0D0
        A2L = (A*A)**(L+1)
        DO K=0,NGSS
C
          Y1  = GAMLOG(2*(2*L+2*K+4))-GAMLOG(2*(K+1))
          Y2  = POLYLOG(2*L+2*K+5,-1.0D0/U)
C
          X0C = X0C + PHK*XUP*DFLOAT(2*L+2*K+1)*DEXP(Y1)*Y2
C
C         UPDATE PHASES AND THINGS
          PHK =-PHK
          XUP = XUP*X*U*U
C
        ENDDO
        X0C = 3.0D0*U3*A2L*X0C
C
C       INTEGRATE REMAINING TERMS FROM C TO INFINITY
        XCI = 0.0D0
        PHK = 1.0D0
        XK  = 1.0D0
        DO K=0,NGSS
C
          U2I = 1.0D0
          TKL = 0.0D0
          DO I=0,K+L
C
C           PREPARATION
            V1 = GAMLOG(2*(2*L+2*K    +3)) - GAMLOG(2*(K+1))
     &         - GAMLOG(2*(2*L+2*K-2*I+2))
            V2 = 2.0D0*U
            V3 = POLYLOG(2*I+4,0.0D0)
            V4 = DFLOAT(2*L+2*K+3)/DFLOAT(2*L+2*K+2-2*I)
            V5 = POLYLOG(2*I+3,0.0D0)
C
C           CONTRIBUTION TO SUM
            TKL = TKL + U2I*DEXP(V1)*(V2*V3-V4*V5)
C
C           UPDATE POWERS OF U
            U2I = U2I*U*U
C
          ENDDO
C
C         CONTRIBUTION TO SUM
          XCI = XCI + PHK*XK*TKL
C
C         EXTRA TERM
          Y1 = GAMLOG(2*(2*L+2*K+4)) - GAMLOG(2*(K+1))
          Y2 = POLYLOG(2*L+2*K+5,-1.0D0/U)
          Y3 = (U*U)**(L+K+1)
C
          XCI = XCI - PHK*XK*DEXP(Y1)*Y2*Y3
C
C         UPDATE PHASES AND THINGS
          PHK =-PHK
          XK  = XK*X
C
        ENDDO
        XCI = 6.0D0*U3*C2L*XCI
        XCI = 0.0D0
C
      ELSEIF(A.GT.1.0D0) THEN
C     EXPAND THE EXPONENTIAL AS A TAYLOR SERIES
C
C       TRUNCATION DETAILS
        NEXP = NTR
C
C       INTEGRATE REMAINING TERMS FROM 0 TO C
        X0C = 0.0D0
        XUP = 1.0D0
        DO K=0,NEXP
C
          Y1 =-GAMLOG(2*(K+1))
          Y2 = 2.0D0*X12*U*POLYLOG(3-K,-1.0D0/U)
          Y3 = GAMLWR(2*L+K+2,X)
          Y4 = POLYLOG(2-K,-1.0D0/U)
          Y5 = GAMLWR(2*L+K+3,X)
C
          X0C = X0C + DFLOAT(L+K+2)*DEXP(Y1)*(Y2*Y3-Y4*Y5)/XUP
C
C         UPDATE SOME THINGS
          XUP = XUP*X12*U
C
        ENDDO
        X0C = 1.5D0*U*U*X0C/(X12*(EIJ**(L+1)))
C
C       INTEGRATE REMAINING TERMS FROM C TO INFINITY
        XCI = 0.0D0
        UK  = 1.0D0
        DO K=0,NEXP
C
          XRI = 1.0D0
          TKL = 0.0D0
          PHI = 1.0D0
          DO I=0,K
C
C           PREPARATION
            V1 = GAMLOG(2*(K+1)) - GAMLOG(2*(I+1)) - GAMLOG(2*(K-I+1))
            V2 = 2.0D0*X12*U*POLYLOG(3-K,0.0D0)
            V3 = GAMUPR(2*L+I+2,X)
            V4 = POLYLOG(2-K,0.0D0)
            V5 = GAMUPR(2*L+I+3,X)
C
C           CONTRIBUTION TO SUM
            TKL = TKL + PHI*DEXP(V1)*(V2*V3+V4*V5)/XRI
C
C           UPDATE POWERS OF U AND STUFF
            XRI = XRI*X12
            PHI =-PHI
C
          ENDDO
C
C         CONTRIBUTION TO SUM
          XCI = XCI + TKL/UK
C
C         UPDATE PHASES AND THINGS
          UK  = UK*U
C
        ENDDO
        XCI = 1.5D0*U*U*XCI/(X12*(EIJ**(L+1)))
C
      ELSE
C     NO OTHER CHOICE -- QUADRATURE!
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
          Y2  = DEXP(-X*ZN*ZN)
          Y3  =-ZN*POLYLOG(2,(ZN-1.0D0)/U)
          Y4  = 2.0D0*U*POLYLOG(3,(ZN-1.0D0)/U)
          X0C = X0C + EXTINT11(Y1*Y2*(Y3+Y4),N,NLIN)
        ENDDO
        X0C = 15.0D0*U*U*C2L*ELIN*X0C/299376.0D0
C
C       NUMERICALLY INTEGRATE REMAINING TERMS FROM C TO INFINITY (EXP.)
        XCI = 0.0D0
        DO N=0,NEXP
          TN  = EEXP*DFLOAT(N)
          ZN  = DEXP(TN)
          Y1  = ZN**(2*L+2)
          Y2  = DEXP(-X*ZN*ZN)
          Y3  = ZN*POLYLOG(2,(1.0D0-ZN)/U)
          Y4  = 2.0D0*U*POLYLOG(3,(1.0D0-ZN)/U)
          XCI = XCI + EXTINT11(Y1*Y2*(Y3+Y4),N,NEXP)
        ENDDO
        XCI = 15.0D0*U*U*C2L*EEXP*XCI/299376.0D0
C      
      ENDIF
C      
      FMIINT0 = FMIINT0 + FNRM*(X0C+XCI)
C
      RETURN
      END

