      SUBROUTINE FMINC0(VMAT,EXL,IZ,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         FFFFFFFF MM       MM IIII NN    NN  CCCCCC   000000          C
C         FF       MMM     MMM  II  NNN   NN CC    CC 00   000         C
C         FF       MMMM   MMMM  II  NNNN  NN CC       00  0000         C
C         FFFFFF   MM MM MM MM  II  NN NN NN CC       00 00 00         C
C         FF       MM  MMM  MM  II  NN  NNNN CC       0000  00         C
C         FF       MM   M   MM  II  NN   NNN CC    CC 000   00         C
C         FF       MM       MM IIII NN    NN  CCCCCC   000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  FMINC0 CALCULATES ATOM-CENTRED NUCLEAR ATTRACTION MATRIX ELEMENTS   C
C  FOR A FERMI NUCLEAR CHARGE MODEL.                                   C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RN(MBS*MBS,4),VMAT(MBD,MBD),EXL(MBS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,MFT),XNUC(MCT,MFT),NNUC(MCT),NMDL(MCT)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DETERMINE THE LQN
      IF(KQN.LT.0) THEN
        LQN =-KQN-1
      ELSE
        LQN = KQN
      ENDIF
      RL = DFLOAT(LQN)
      TL = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     GAUSSIAN WIDTH PARAMETER
      APRM = AFMI(IZ)
      CPRM = 5.0D0*RNUC(IZ)*RNUC(IZ)/3.0D0
     &     - 7.0D0*PI*PI*AFMI(IZ)*AFMI(IZ)/3.0D0
      CPRM = DSQRT(CPRM)
C
C     WARN USER IF FERMI MODEL IS NOT PHYSICAL
      IF(CPRM.LT.APRM) THEN
        WRITE(6, *) 'In FRMINC0: nucleus is too small for this model.'
        WRITE(7, *) 'In FRMINC0: nucleus is too small for this model.'
      ENDIF
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
          T52 = RL+2.5D0
          E52 = EIJ**T52
C
C         MATRIX ELEMENT DEPENDS ON HAMILTONIAN
          IF(HMLT.EQ.'NORL') THEN
C
C           NON-RELATIVISTIC CASE
            ULL = FMIINT0(LQN  ,EIJ,APRM,CPRM)
C
C           TRANSFER INTO NUCLEAR ATTRACTION MATRIX
            VMAT(IBAS,JBAS) = -ZNUC(IZ)*RN(M,1)*ULL
C
          ELSE
C
C           RELATIVISTIC CASE
            ULL = FMIINT0(LQN  ,EIJ,APRM,CPRM)
C
            SSS = 2.0D0*GAMHLF(2*LQN+5)*EPR/E52
C
            VSA = 4.0D0*EPR*FMIINT0(LQN+1,EIJ,APRM,CPRM)
            IF(KQN.LT.0) THEN
              VSN = 0.0D0
            ELSE
              VSN =-2.0D0*EIJ*TL*FMIINT0(LQN  ,EIJ,APRM,CPRM)
     &                    +TL*TL*FMIINT0(LQN-1,EIJ,APRM,CPRM)
            ENDIF
            USS = VSA+VSN
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NBAS
            LBAS = JBAS+NBAS
C
C           TRANSFER INTO NUCLEAR ATTRACTION MATRIX
            VMAT(IBAS,JBAS) = -ZNUC(IZ)*RN(M,1)*ULL
            VMAT(KBAS,LBAS) = -ZNUC(IZ)*RN(M,4)*USS
     &                        -2.0D0*EMSS*CV*CV*RN(M,4)*SSS
C
          ENDIF
C
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
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
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/FCTS/RFACT(0:80),SFACT(0:80)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
      DIMENSION U(0:60),P(24)
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In FMIINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In FMIINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
C     MULTIPLES OF THE FACTOR U
      U(0) = 1.0D0
      DO K=1,60
        U(K) = U(K-1)*A/C
      ENDDO
C
C     SPECIAL VALUES OF THE RIEMANN-ZETA FUNCTION
      P( 4) =-(PI**4)*7.0D0/720.0D0
      P( 6) =-(PI**6)*31.0D0/30240.0D0
      P( 8) =-(PI**8)*127.0D0/12096.0D2
      P(10) =-(PI**10)*73.0D0/684288.0D1
      P(12) =-(PI**12)*1414477.0D0/1307674368.0D3
      P(14) =-(PI**14)*8191.0D0/747242496.0D2
      P(16) =-(PI**16)*16931177.0D0/152437569184.0D4
      P(18) =-(PI**18)*5749691557.0D0/5109094217170944.0D3
      P(20) =-(PI**20)*91546277357.0D0/8028576626982912.0D5
      P(22) =-(PI**22)*3324754717.0D0/28777755182432256.0D4
      P(24) =-(PI**24)*1982765468311237.0D0/16938241367317436694528.0D5
C
C     CORRECT DIMENSIONS ARISE FROM THESE TERMS
      EL1 = 1.0D0/(EIJ**(L+1))
      C2L = C**(2*L+2)
C
C     DIMENSIONLESS GAUSSIAN PARAMETER
      SIGMA = EIJ*C*C
C
C     FERMI NORMALISATION CONSTANT
      FNRM = 1.0D0 + PI*PI*U(2) - 6.0D0*U(3)*POLYLOG(3,-1.0D0/U(1))
      FNRM = 1.0D0/FNRM
C
C     INTEGRAL TYPE: INCOMPLETE GAMMA FUNCTIONS
      FA  = 1.5D0 + PI*PI*U(2)*0.5D0
      FB  =-0.5D0
      FC  = 1.0D0 + PI*PI*U(2)
      Y1A = 0.5D0*GAMLWR(2*L+3,SIGMA)/(DSQRT(SIGMA))
      Y1B = 0.5D0*GAMLWR(2*L+5,SIGMA)/(SIGMA*DSQRT(SIGMA))
      Y1C = 0.5D0*GAMUPR(2*L+2,SIGMA)
      Y1  = FA*Y1A + FB*Y1B + FC*Y1C
C
C     INTEGRAL TYPE: SEMI-INFINITE GAMMA FUNCTIONS
      FD  =-6.0D0*U(3)*POLYLOG(3,-1.0D0/U(1))
      Y2D = 0.5D0*RFACT(L)
      
      Y2  = FD*Y2D
C
C     INTEGRAL TYPE: POLYLOG FUNCTIONS WITH GAUSSIANS AND POLYNOMIALS
C                    (EVALUATION METHOD DEPENDS ON SIGMA PARAMETER)
      Y34 = 0.0D0
C
C     EXPAND THE GAUSSIAN FACTOR IN TAYLOR SERIES AND ASSUME LARGE N/U
      IF(SIGMA.LT.1.0D0) THEN
C
C       MAXIMUM EXPANSION ORDER FOR GAUSSIAN
        KMAX = 20
C
C       INITIALISE BINS FOR GAUSSIAN EXPANSION SERIES     
        PHK = 1.0D0
        SGK = 1.0D0
        FTK = 1.0D0
C
C       TAYLOR EXPANSION OVER THE GAUSSIAN FOR SMALL SIGMA
        DO K=0,KMAX
C
          BIN = U(1)*U(2*L+2*K)*POLYLOG(2*L+2*K+5,-1.0D0/U(1))
          DO J=0,MIN(L+K,10)
C
            PR = 1.0D0/RFACT(2*L+2*K-2*J+1)
            BIN = BIN + 2.0D0*U(2*J)*P(2*J+4)*PR
C
          ENDDO
C             
          Y34 = Y34 + PHK*SGK*DFLOAT(L+K+2)*RFACT(2*L+2*K+1)*BIN/FTK
C
C         UPDATE PHASES, FACTORIALS AND POLYNOMIALS IN K
          PHK =-PHK
          SGK = SGK*SIGMA
          FTK = FTK*DFLOAT(K+1)
C          
        ENDDO
C     
C       INCLUDE MULTIPLICATIVE FACTOR
        Y34 = 6.0D0*U(4)*Y34
C
C     NUMERICAL QUADRATURE (WORKS FOR ANY CASE BUT COMPUTATIONALLY SLOW)
      ELSE
C
C       INTEGRATION GRID DETAILS
        NLIN = 100
        NEXP = 1000
        RMAX = 15.0D0
        ELIN = 1.0D0/DFLOAT(NLIN)
        EEXP = DLOG(RMAX/C)/DFLOAT(NEXP)
C
C       NUMERICALLY INTEGRATE REMAINING TERMS FROM 0 TO C (LINEAR)
        X0C = 0.0D0
        DO N=0,NLIN
          ZN  = ELIN*DFLOAT(N)
          Z1  = ZN**(2*L+1)
          Z2  = DEXP(-SIGMA*ZN*ZN)
          Z3  =-ZN*POLYLOG(2,(ZN-1.0D0)/U(1))
          Z4  = 2.0D0*U(1)*POLYLOG(3,(ZN-1.0D0)/U(1))
          X0C = X0C + EXTINT11(Z1*Z2*(Z3+Z4),N,NLIN)
        ENDDO
        X0C = 15.0D0*U(2)*ELIN*X0C/299376.0D0
C
C       NUMERICALLY INTEGRATE REMAINING TERMS FROM C TO INFINITY (EXP.)
        XCI = 0.0D0
        DO N=0,NEXP
          TN  = EEXP*DFLOAT(N)
          ZN  = DEXP(TN)
          Z1  = ZN**(2*L+2)
          Z2  = DEXP(-SIGMA*ZN*ZN)
          Z3  = ZN*POLYLOG(2,(1.0D0-ZN)/U(1))
          Z4  = 2.0D0*U(1)*POLYLOG(3,(1.0D0-ZN)/U(1))
          XCI = XCI + EXTINT11(Z1*Z2*(Z3+Z4),N,NEXP)
        ENDDO
        XCI = 15.0D0*U(2)*EEXP*XCI/299376.0D0
C        
        Y34 = X0C+XCI
C
      ENDIF
C
C     COMBINE ALL TERMS AND APPLY NORMALISATION FACTOR
      FMIINT0 = FNRM*(EL1*Y1 + EL1*Y2 + C2L*Y34)
C
      RETURN
      END

