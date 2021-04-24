      SUBROUTINE OVRLP0(OVAP,EXL,KQN,NFUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          OOOOOO  VV    VV RRRRRRR  LL      PPPPPPP   000000          C
C         OO    OO VV    VV RR    RR LL      PP    PP 00    00         C
C         OO    OO VV    VV RR    RR LL      PP    PP 00    00         C
C         OO    OO VV    VV RR    RR LL      PP    PP 00    00         C
C         OO    OO  VV  VV  RRRRRRR  LL      PPPPPPP  00    00         C
C         OO    OO   VVVV   RR    RR LL      PP       00    00         C
C          OOOOOO     VV    RR    RR LLLLLLL PP        000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  ONEELRE0 CALCULATES THE DIRAC AND OVERLAP MATRICES FOR SYMMETRY     C
C  TYPE KQN, USING EVEN-TEMPERED SGTFS.                                C
C**********************************************************************C
      PARAMETER(MDM=2500,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26)
C
      CHARACTER*4 HMLTN
C
      DIMENSION RN(MBS*MBS,4),OVAP(2*MBS,2*MBS),EXL(MBS)
C
      COMMON/GMFN/GAMMAL(100),GAMMAF(100)
      COMMON/PRMS/CV,HMLTN,INEW,ITREE,IMOL,ILEV,IEQS
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     DETERMINE THE LQN
      IF(KQN.LT.0) THEN
        LQN =-(KQN+1)
      ELSE
        LQN =  KQN
      ENDIF
      RL = DFLOAT(LQN)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NFUN,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NFUN
        EI = EXL(IBAS)
        DO JBAS=1,NFUN
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
          T32 = RL+1.5D0
          T52 = RL+2.5D0
          E32 = EIJ**T32
          E52 = EIJ**T52
          SLL = 2.0D0*RN(M,1)*GAMMAF(2*LQN+3)/E32
          SSS = 8.0D0*RN(M,3)*GAMMAF(2*LQN+5)*EPR/E52
C
C         OVERLAP MATRIX ELEMENTS
          OVAP(IBAS     ,JBAS     ) = SLL
          IF(HMLTN.EQ.'NORL') GOTO 50
          OVAP(IBAS+NFUN,JBAS     ) = 0.0D0
          OVAP(JBAS     ,IBAS+NFUN) = 0.0D0
          OVAP(IBAS+NFUN,JBAS+NFUN) = SSS
50        CONTINUE
C         
        ENDDO
      ENDDO
C
      RETURN
      END
