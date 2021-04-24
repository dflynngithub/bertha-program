      SUBROUTINE UNINC0(HMAT,EXL,ZCRG,ZRMS,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          UU    UU NN    NN IIII NN    NN  CCCCCC   000000            C
C          UU    UU NNN   NN  II  NNN   NN CC    CC 00   000           C
C          UU    UU NNNN  NN  II  NNNN  NN CC       00  0000           C
C          UU    UU NN NN NN  II  NN NN NN CC       00 00 00           C
C          UU    UU NN  NNNN  II  NN  NNNN CC       0000  00           C
C          UU    UU NN   NNN  II  NN   NNN CC    CC 000   00           C
C           UUUUUU  NN    NN IIII NN    NN  CCCCCC   000000            C
C                                                                      C
C -------------------------------------------------------------------- C
C  UNINC0 CALCULATES ATOM-CENTRED NUCLEAR ATTRACTION MATRIX ELEMENTS   C
C  FOR A UNIFORMLY-CHARGED SPHERICAL NUCLEAR MODEL.                    C
C**********************************************************************C
      INCLUDE 'param.h'
C
      CHARACTER*4 HMLT
C
      DIMENSION RN(MBS*MBS,4),HMAT(MBD,MBD),EXL(MBS)
C
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,TW12
      COMMON/PHYS/CV,PRTN,GFREE,GFRMI,WEIN,CHZ,CFM,CDB
      COMMON/PRMS/HMLT,IOPT,IMOL,INEW,ILEV,ISML,ISWZ,IEQS,IERC,IPAR,ICOR
C
C     DETERMINE THE LQN
      IF(KQN.LT.0) THEN
        LQN =-KQN-1
      ELSE
        LQN = KQN
      ENDIF
      TL = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
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
          E02 = EIJ**(LQN)
          E22 = E02*EIJ
          E42 = E22*EIJ
          E52 = E42*DSQRT(EIJ)
C
          X   = 5.0D0*EIJ*ZRMS*ZRMS/3.0D0
C
C         MATRIX ELEMENT DEPENDS ON HAMILTONIAN
          IF(HMLT.EQ.'NORL') THEN
C
C           NON-RELATIVISTIC CASE
            ULL = GAMUPR(LQN+1,X) + PHIX(LQN+1,X)
            ULL = 0.5D0*ULL/E22
C
C           TRANSFER INTO DIRAC MATRIX
            HMAT(IBAS,JBAS) = HMAT(IBAS,JBAS) - ZCRG*RN(M,1)*ULL
C
          ELSE
C
C           RELATIVISTIC CASE
            ULL = GAMUPR(LQN+1,X ) + PHIX(LQN+1,X)
            ULL = 0.5D0*ULL/E22
C
            SSS = 2.0D0*GAMHLF(2*LQN+5)*EPR/E52
            VS2 = GAMUPR(LQN+2,X ) + PHIX(LQN+2,X)
            VS2 = 2.0D0*EPR*VS2/E42
C
            IF(KQN.LT.0) THEN
              VS0 = 0.0D0
              VS1 = 0.0D0
            ELSE
C
              VS0 = GAMUPR(LQN  ,X ) + PHIX(LQN  ,X)
              VS1 = GAMUPR(LQN+1,X ) + PHIX(LQN+1,X)
C
              VS0 = 0.5D0*TL*TL*VS0/E02
              VS1 =-TL*VS1/E02
C
            ENDIF
            USS = VS0+VS1+VS2
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NBAS
            LBAS = JBAS+NBAS
C
C           TRANSFER INTO ARRAY
            HMAT(IBAS,JBAS) = HMAT(IBAS,JBAS) - ZCRG*RN(M,1)*ULL
            HMAT(KBAS,LBAS) = HMAT(KBAS,LBAS) - ZCRG*RN(M,3)*USS
     &                                  -2.0D0*CV*CV*RN(M,3)*SSS
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
      FUNCTION PHIX(L,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                   PPPPPPP  HH    HH IIII XX     XX                   C
C                   PP    PP HH    HH  II   XX   XX                    C
C                   PP    PP HH    HH  II    XX XX                     C
C                   PP    PP HHHHHHHH  II     XXX                      C
C                   PPPPPPP  HH    HH  II    XX XX                     C
C                   PP       HH    HH  II   XX   XX                    C
C                   PP       HH    HH IIII XX     XX                   C
C                                                                      C
C -------------------------------------------------------------------- C
C  PHIX IS AN AUXILLIARY FUNCTION WHICH EVALUATES AS A POWER SERIES    C
C  OF TRUNCTATED LENGTH N:0->20 THE COMBINATION OF LOWER INCOMPLETE    C
C  GAMMA FUNCTIONS NEEDED FOR ATOMIC UNIFORM NUCLEUS MATRIX ELEMENTS:  C
C                                                                      C
C  Φ(ℓ,x) = [3*x*γ(ℓ+1/2,x)-γ(ℓ+3/2,x)]/8*π*x^(3/2) with ℓ>0 and x>0.  C
C                                                                      C
C**********************************************************************C
      INCLUDE 'param.h'
C
      COMMON/FCTS/RFACT(0:20),SFACT(0:20)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,TW12
C
C     ROUTINE ONLY ALLOWS PARAMETERS L>0
      IF(L.LT.1) THEN
        WRITE(6, *) 'In PHIX: order L must be at least 1. L = ',L
        WRITE(7, *) 'In PHIX: order L must be at least 1. L = ',L
      ENDIF
C
C     PRE-MULTIPLYING FACTOR
      FC1 = X**(L)/PI
C
C     INITIALISE POLYNOMIAL COUNTERS
      PLY = 0.0D0
      PHS = 1.0D0
      XPW = 1.0D0
C
C     LOOP OVER REQUIRED POLYNOMIAL DEGREES (CHOOSE 20)
      DO N=0,20
        RAT = DFLOAT(L+N+2)/DFLOAT((2*L+2*N+1)*(2*L+2*N+3))
        PLY = PLY + PHS*XPW*RAT/RFACT(N)
        XPW = XPW*X
        PHS =-PHS
      ENDDO
C
C     VALUE OF INTEGER-ORDER UPPER INCOMPLETE GAMMA FUNCTION
      PHIX = FC1*PLY
C
      RETURN
      END
C
C
      FUNCTION GAMUPR(L,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        GGGGGG     AA    MM       MM UU    UU PPPPPPP  RRRRRRR        C
C       GG    GG   AAAA   MMM     MMM UU    UU PP    PP RR    RR       C
C       GG        AA  AA  MMMM   MMMM UU    UU PP    PP RR    RR       C
C       GG       AA    AA MM MM MM MM UU    UU PP    PP RR    RR       C
C       GG   GGG AAAAAAAA MM  MMM  MM UU    UU PPPPPPP  RRRRRRR        C
C       GG    GG AA    AA MM   M   MM UU    UU PP       RR    RR       C
C        GGGGGG  AA    AA MM       MM  UUUUUU  PP       RR    RR       C
C                                                                      C
C -------------------------------------------------------------------- C
C  GAMUPR RETURNS THE UPPER INCOMPLETE GAMMA FUNCTION FOR AN INTEGER   C
C  ARGUMENT L ONLY. NO SERIES OR CONTINUED FRACTION APPROXIMATIONS ARE C
C  NECESSARY -- AN ALGEBRAIC SOLUTION HAS BEEN DEDUCED.                C
C                                                                      C
C                GAMUPR(ℓ,x) = Γ(ℓ,x) with ℓ>0 and x>0.                C
C                                                                      C
C**********************************************************************C
      INCLUDE 'param.h'
C
      COMMON/FCTS/RFACT(0:20),SFACT(0:20)
C
C     ROUTINE ONLY ALLOWS PARAMETERS L>0
      IF(N.LT.1) THEN
        WRITE(6, *) 'In GAMUPR: order L must be at least 1. L = ',L
        WRITE(7, *) 'In GAMUPR: order L must be at least 1. L = ',L
      ENDIF
C
C     EXPONENTIAL FACTOR
      FC1 = DEXP(-X)
C
C     INITIALISE POLYNOMIAL COUNTERS
      POLY = 0.0D0
      XPOW = 1.0D0
C
C     LOOP OVER REQUIRED POLYNOMIAL DEGREES
      DO I=0,L-1
        POLY = POLY + XPOW*RFACT(L-1)/RFACT(I)
        XPOW = XPOW*X
      ENDDO
C
C     VALUE OF INTEGER-ORDER UPPER INCOMPLETE GAMMA FUNCTION
      GAMUPR = FC1*POLY
C
      RETURN
      END

