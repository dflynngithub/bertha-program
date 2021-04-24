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
C**********************************************************************C
      INCLUDE 'param.h'
C
      COMMON/FCTS/RFACT(0:20),SFACT(0:20)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In FMIINT0: order L must be positive. L = ',L
        WRITE(7, *) 'In FMIINT0: order L must be positive. L = ',L
        STOP
      ENDIF
C
C     FACTORS NEEDED FOR ALL PARAMETERS N
      U   = A/C
      X   = EIJ*C*C
      X12 = DSQRT(X)
C
C     FERMI NORMALISATION CONSTANT
      FNRM = 1.0D0 + PI*PI*U*U - 6.0D0*U*U*U*POLYLOG(-1.0D0/U,3)
      FNRM = 1.0D0/FNRM
C
C     INTEGRAL TYPE: INCOMPLETE GAMMA FUNCTION
      X1A = GAMLWR(2*L+3,X)*(3.0D0+PI*PI*U*U)/(2.0D0*X12)
      X1B =-GAMLWR(2*L+5,X)/(2.0D0*X12*X12*X12)
      X1C = GAMUPR(2*L+2,X)*(1.0D0+PI*PI*U*U)
C
C     INTEGRAL TYPE: LOGARITHMIC INTEGRAL
      X2A =-6.0D0*GAMHLF(2*L)*U*U*U*POLYLOG(-1.0D0/U,3)
C
C     INTEGRAL TYPE: INVOLVES S2(R)
      PHN = 1.0D0
      GSS = DEXP(-X)
      UX3 = 3.0D0*U*U/X12
C
      X3A = 0.0D0
      X3B = 0.0D0
      X3C = 0.0D0
      
      IF(X.GT.0.1D0) THEN
      DO N=1,NTR
C
C       PHASE FACTORS AND EXPONENTIALS
        PHN =-PHN
        RN2 = DFLOAT(N*N)

        ENU = DEXP(-DFLOAT(N)/U)
        XUN = DFLOAT(N*N)/(4.0D0*X*U*U)
        BGR = 1.0D0
        BLS = 1.0D0
        WGR = (DSQRT(XUN) + X12)**2
        WLS = (DSQRT(XUN) - X12)**2
        PMD = (DSQRT(XUN)      )**(2*L+3)
        PGR = (DSQRT(XUN) + X12)**(2*L+3)
        PLS = (DSQRT(XUN) - X12)**(2*L+3)
C
C       SERIES EXPANSION FOR 2L+2 TERM
        PHI = 1.0D0
        SM3 = 0.0D0
        SL3 = 0.0D0
        SG3 = 0.0D0
        DO I=0,2*L+2
C
C         BINOMIAL COEFFICIENT
          BLI = RFACT(2*L+2)/(RFACT(I)*RFACT(2*L+2-I))
C
C         INDEX FOR CONFLUENT HYPERGEOMETRIC FUNCTION
          AH = 0.5D0*DFLOAT(2*L+3-I)
C
C         UPDATE MIDDLE TERM
          SM3 = SM3 + BLI    *PMD    *GCF(AH,XUN)
C
C         UPDATE LOWER AND HIGHER TERMS
          SG3 = SG3 + BLI*PHI*PGR*BGR*GCF(AH,WGR)
          SL3 = SL3 + BLI*PHI*PLS*BLS*GCF(AH,WLS)
C
C         UPDATE BGR AND BLS
          BGR = BGR/(1.0D0 + 2.0D0*X*U/DFLOAT(N))
          BLS = BLS/(1.0D0 - 2.0D0*X*U/DFLOAT(N))
          PHI =-PHI
C
        ENDDO
C
C       UPDATE COUNTERS
        X3A = X3A - UX3*ENU*PHN*SM3/RN2
        X3B = X3B + UX3*GSS*PHN*SG3/RN2
        X3C = X3C - UX3*GSS*PHN*SL3/RN2
C
      ENDDO
      ENDIF
C
C     INTEGRAL TYPE: INVOLVES S3(R)
      PHN = 1.0D0
      GSS = DEXP(-X)
      U36 = 6.0D0*U*U*U
C
      X4A = 0.0D0
      X4B = 0.0D0
      X4C = 0.0D0
      
      IF(X.GT.0.1D0) THEN
      DO N=1,NTR
C
C       PHASE FACTORS AND EXPONENTIALS
        PHN =-PHN
        RN3 = DFLOAT(N*N*N)

        ENU = DEXP(-DFLOAT(N)/U)
        XUN = DFLOAT(N*N)/(4.0D0*X*U*U)
        BGR = 1.0D0
        BLS = 1.0D0
        WGR = (DSQRT(XUN) + X12)**2
        WLS = (DSQRT(XUN) - X12)**2
        PMD = (DSQRT(XUN)      )**(2*L+2)
        PGR = (DSQRT(XUN) + X12)**(2*L+2)
        PLS = (DSQRT(XUN) - X12)**(2*L+2)
C
C       SERIES EXPANSION FOR 2L+2 TERM
        PHI = 1.0D0
        SM4 = 0.0D0
        SL4 = 0.0D0
        SG4 = 0.0D0
        DO I=0,2*L+1
C
C         BINOMIAL COEFFICIENT
          BLI = RFACT(2*L+1)/(RFACT(I)*RFACT(2*L+1-I))
C
C         INDEX FOR CONFLUENT HYPERGEOMETRIC FUNCTION
          AH = 0.5D0*DFLOAT(2*L+2-I)
C
C         UPDATE MIDDLE TERM
          SM4 = SM4 + BLI    *PMD    *GCF(AH,XUN)
C
C         UPDATE LOWER AND HIGHER TERMS
          SG4 = SG4 + BLI*PHI*PGR*BGR*GCF(AH,WGR)
          SL4 = SL4 + BLI*PHI*PLS*BLS*GCF(AH,WLS)
C
C         UPDATE BGR AND BLS
          BGR = BGR/(1.0D0 + 2.0D0*X*U/DFLOAT(N))
          BLS = BLS/(1.0D0 - 2.0D0*X*U/DFLOAT(N))
          PHI =-PHI
C
        ENDDO
C
C       UPDATE COUNTERS
        X4A = X4A + U36*ENU*PHN*SM4/RN3
        X4B = X4B + U36*GSS*PHN*SG4/RN3
        X4C = X4C - U36*GSS*PHN*SL4/RN3
C
      ENDDO
      ENDIF
C
C     TRANSFER DATA TO FMIINT0
      X1 = X1A+X1B+X1C
      X2 = X2A
      X3 = X3A+X3B+X3C
      X4 = X4A+X4B+X4C
C
      FMIINT0 = 0.5D0*FNRM*(X1+X2+X3+X4)/(EIJ**(L+1))
C
      RETURN
      END
