      SUBROUTINE BREIT0(B11,B21,B12,B22,D1,D2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           BBBBBBB  RRRRRRR  EEEEEEEE IIII TTTTTTTT 000000            C
C           BB    BB RR    RR EE        II     TT   00    00           C
C           BB    BB RR    RR EE        II     TT   00    00           C
C           BBBBBBB  RR    RR EEEEEE    II     TT   00    00           C
C           BB    BB RRRRRRR  EE        II     TT   00    00           C
C           BB    BB RR    RR EE        II     TT   00    00           C
C           BBBBBBB  RR    RR EEEEEEEE IIII    TT    000000            C
C                                                                      C
C -------------------------------------------------------------------- C
C  BREIT0 CONSTRUCTS THE ATOMIC BREIT MATRIX FROM RADIAL DIRECT AND    C
C  EXCHANGE INTEGRALS AND A MEAN-FIELD CHARGE DENSITY.                 C
C -------------------------------------------------------------------- C
C  DFNOTE: NOT FINISHED YET. EXCHANGE TERMS INCORRECT FOR S-TYPE GTFS. C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MAXAB=20)
C
      DIMENSION B11(2*MBS,2*MBS),B21(2*MBS,2*MBS),
     &          B12(2*MBS,2*MBS),B22(2*MBS,2*MBS)
      DIMENSION D1(MB2,3),D2(MB2,3),RN(MB2,4)
C
      COMMON/BLOC/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/BSIS/EXLA(MBS),EXLB(MBS)
      COMMON/BTRE/RKLLSS(MB2,4),RKSLSL(MB2,4),RKSSLL(MB2,4),
     &            RMSLSL(MB2,4)
      COMMON/IJ/EIJ(MAXAB),RNIJ(3),EI,EJ
      COMMON/KL/EK(MB2),EKL0(MB2),EIK(MB2,MAXAB),EJL(MB2,MAXAB),
     &          EL(MB2),EKL1(MB2),EKL(MB2,MAXAB),RNKL(MB2,4),
     &          B1(MB2,MAXAB,MAXAB),B2(MB2,MAXAB,MAXAB),
     &          B3(MB2,MAXAB,MAXAB),B4(MB2,MAXAB,MAXAB),
     &          B5(MB2,MAXAB,MAXAB),B6(MB2,MAXAB,MAXAB)
C
C     GENERATE 'KL' EXPONENT POWERS FOR LATER B-INTEGRAL GENERATION
      CALL KLSET
C
C     GENERATE A BATCH OF NORMALISATION CONSTANTS
      CALL RNORMA(RN,EXLA,NFUNA,LQNA)
C
C     INITIALISE BREIT MATRIX
      DO IBAS=1,2*NFUNA
        DO JBAS=1,2*NFUNA
          B11(IBAS,JBAS) = 0.0D0
          B21(IBAS,JBAS) = 0.0D0
          B12(IBAS,JBAS) = 0.0D0
          B22(IBAS,JBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     ITERATE OVER ALL MATRIX ELEMENTS
      IJ = 0
      DO IBAS=1,NFUNA
        EI = EXLA(IBAS)
        DO JBAS=1,NFUNA
C         BASIS EXPONENT COMBINATIONS FOR LATER B-INTEGRAL GENERATION
          IJ = IJ+1
C
          KBAS = IBAS + NFUNA
          LBAS = JBAS + NFUNA
C
          EJ = EXLA(JBAS)
          EIJ0 = EI + EJ
          EIJR = DSQRT(EIJ0)
          EIJA = EIJ0**(-LQNA)
          DO N=1,6
            EIJ(N) = EIJA
            EIJA   = EIJA/EIJR
          ENDDO
          RNIJ(1) = RN(IJ,1)
          RNIJ(2) = RN(IJ,2)
          RNIJ(3) = RN(IJ,3)
C
C         GENERATE A BATCH OF BETA INTEGRALS
          CALL BINTGEN
C
C         GENERATE BATCHES OF RADIAL INTEGRALS (J AND K MATRICES)
          CALL BII0
C
C    (11) KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 10
C
C         INITIALISE COUNTERS
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0
C
C         SUM OVER MEAN FIELD CONTRIBUTIONS FOR THIS BASIS PAIR
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,1)*D1(M,3)
            BSL = BSL + RKSLSL(M,1)*D1(M,2) + RMSLSL(M,1)*D1(M,2)
            BSS = BSS + RKSSLL(M,1)*D1(M,1)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B11(IBAS,JBAS) = BLL
          B11(KBAS,JBAS) = BSL
          B11(JBAS,KBAS) = BSL
          B11(KBAS,LBAS) = BSS
C
10        CONTINUE
C
C    (21) KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNB.EQ.0) GOTO 11
C
C         INITIALISE COUNTERS
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0
C
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,2)*D1(M,3)
            BSL = BSL + RKSLSL(M,2)*D1(M,2) + RMSLSL(M,2)*D1(M,2)
            BSS = BSS + RKSSLL(M,2)*D1(M,1)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B21(IBAS,JBAS) = BLL
          B21(KBAS,JBAS) = BSL
          B21(JBAS,KBAS) = BSL
          B21(KBAS,LBAS) = BSS
C
11        CONTINUE

C
C    (12) KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0) GOTO 12
C
C         INITIALISE COUNTERS          
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0
C
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,3)*D2(M,3)
            BSL = BSL + RKSLSL(M,3)*D2(M,2) + RMSLSL(M,3)*D2(M,2)
            BSS = BSS + RKSSLL(M,3)*D2(M,1)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B12(IBAS,JBAS) = BLL
          B12(KBAS,JBAS) = BSL
          B12(JBAS,KBAS) = BSL
          B12(KBAS,LBAS) = BSS
C
12        CONTINUE
C
C    (22) KQNA < 0 AND KQNB < 0  CONTRIBUTIONS (CANNOT SKIP)
C
C         INITIALISE COUNTERS          
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0
C
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,4)*D2(M,3)
            BSL = BSL + RKSLSL(M,4)*D2(M,2) + RMSLSL(M,4)*D2(M,2)
            BSS = BSS + RKSSLL(M,4)*D2(M,1)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B22(IBAS,JBAS) = BLL
          B22(KBAS,JBAS) = BSL
          B22(JBAS,KBAS) = BSL
          B22(KBAS,LBAS) = BSS
C
        ENDDO
      ENDDO
C
      RETURN
      END

