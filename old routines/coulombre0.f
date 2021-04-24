      SUBROUTINE COULOMBRE0(G11,G21,G12,G22,D1,D2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    CCCCCC  LL      MM       MM BBBBBBB  RRRRRRR  EEEEEEEE  000000    C
C   CC    CC LL      MMM     MMM BB    BB RR    RR EE       00    00   C
C   CC       LL      MMMM   MMMM BB    BB RR    RR EE       00    00   C
C   CC       LL      MM MM MM MM BBBBBBB  RR    RR EEEEEE   00    00   C
C   CC       LL      MM  MMM  MM BB    BB RRRRRRR  EE       00    00   C
C   CC    CC LL      MM   M   MM BB    BB RR    RR EE       00    00   C
C    CCCCCC  LLLLLLL MM       MM BBBBBBB  RR    RR EEEEEEEE  000000    C
C                                                                      C
C -------------------------------------------------------------------- C
C  COULOMBRE0 CONSTRUCTS THE ATOMIC COULOMB MATRIX FROM RADIAL DIRECT  C
C  AND EXCHANGE INTEGRALS AND A MEAN-FIELD CHARGE DENSITY.             C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MNU=MKP+1,MAB=2*MNU+5)
C
      DIMENSION G11(2*MBS,2*MBS),G21(2*MBS,2*MBS),
     &          G12(2*MBS,2*MBS),G22(2*MBS,2*MBS)
      DIMENSION D1(MB2,3),D2(MB2,3),RN(MB2,4)
C
      COMMON/BLOC/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/BSIS/EXLA(MBS),EXLB(MBS)
      COMMON/CLRE/RKLLLL(MB2,4),RKSSSS(MB2,4),RKSLSL(MB2,4),
     &            RJLLLL(MB2,4),RJSSSS(MB2,4),RJLLSS(MB2,4),
     &            RJSSLL(MB2,4)
      COMMON/IJ/EIJ(MAB),RNIJ(3),EI,EJ
      COMMON/KL/EK(MB2),EKL0(MB2),EIK(MB2,MAB),EJL(MB2,MAB),
     &          EL(MB2),EKL1(MB2),EKL(MB2,MAB),RNKL(MB2,4),
     &          B1(MB2,MAB,MAB),B2(MB2,MAB,MAB),B3(MB2,MAB,MAB),
     &          B4(MB2,MAB,MAB),B5(MB2,MAB,MAB),B6(MB2,MAB,MAB)
C
C     GENERATE 'KL' EXPONENT POWERS FOR LATER B-INTEGRAL GENERATION
      CALL KLSET
C
C     GENERATE A BATCH OF NORMALISATION CONSTANTS
      CALL RNORM0(RN,EXLA,NFUNA,LQNA)
C
C     INITIALISE COULOMB MATRIX
      DO IBAS=1,2*NFUNA
        DO JBAS=1,2*NFUNA
          G11(IBAS,JBAS) = 0.0D0
          G21(IBAS,JBAS) = 0.0D0
          G12(IBAS,JBAS) = 0.0D0
          G22(IBAS,JBAS) = 0.0D0
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
          CALL BINT
C
C         GENERATE BATCH OF RADIAL INTEGRALS (J AND K MATRICES)
          CALL ERIRE0
C
C    (11) KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 10
C
C         INITIALISE COUNTERS          
          GLL = 0.0D0
          GSL = 0.0D0
          GSS = 0.0D0
C
C         SUM OVER MEAN FIELD CONTRIBUTIONS FOR THIS BASIS PAIR
          DO M=1,MAXM
            GLL = GLL + RJLLLL(M,1)*D1(M,1) - RKLLLL(M,1)*D1(M,1)
     &                + RJLLSS(M,1)*D1(M,3)
            GSL = GSL                       - RKSLSL(M,1)*D1(M,2)
            GSS = GSS + RJSSSS(M,1)*D1(M,3) - RKSSSS(M,1)*D1(M,3)
     &                + RJSSLL(M,1)*D1(M,1)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO COULOMB MATRIX
          G11(IBAS,JBAS) = GLL
          G11(KBAS,JBAS) = GSL
          G11(JBAS,KBAS) = GSL
          G11(KBAS,LBAS) = GSS
C
10        CONTINUE
C
C    (21) KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNB.EQ.0) GOTO 11
C
C         INITIALISE COUNTERS          
          GLL = 0.0D0
          GSL = 0.0D0
          GSS = 0.0D0
C
          DO M=1,MAXM
            GLL = GLL + RJLLLL(M,2)*D1(M,1) - RKLLLL(M,2)*D1(M,1)
     &                + RJLLSS(M,2)*D1(M,3)
            GSL = GSL                       - RKSLSL(M,2)*D1(M,2)
            GSS = GSS + RJSSSS(M,2)*D1(M,3) - RKSSSS(M,2)*D1(M,3)
     &                + RJSSLL(M,2)*D1(M,1)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO COULOMB MATRIX
          G21(IBAS,JBAS) = GLL
          G21(KBAS,JBAS) = GSL
          G21(JBAS,KBAS) = GSL
          G21(KBAS,LBAS) = GSS
C
11        CONTINUE

C
C    (12) KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0) GOTO 12
C
C         INITIALISE COUNTERS          
          GLL = 0.0D0
          GSL = 0.0D0
          GSS = 0.0D0
C
          DO M=1,MAXM
            GLL = GLL + RJLLLL(M,3)*D2(M,1) - RKLLLL(M,3)*D2(M,1)
     &                + RJLLSS(M,3)*D2(M,3)
            GSL = GSL                       - RKSLSL(M,3)*D2(M,2)
            GSS = GSS + RJSSSS(M,3)*D2(M,3) - RKSSSS(M,3)*D2(M,3)
     &                + RJSSLL(M,3)*D2(M,1)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO COULOMB MATRIX
          G12(IBAS,JBAS) = GLL
          G12(KBAS,JBAS) = GSL
          G12(JBAS,KBAS) = GSL
          G12(KBAS,LBAS) = GSS
C
12        CONTINUE
C
C    (22) KQNA < 0 AND KQNB < 0  CONTRIBUTIONS (CANNOT SKIP)
C
C         INITIALISE COUNTERS          
          GLL = 0.0D0
          GSL = 0.0D0
          GSS = 0.0D0
C
          DO M=1,MAXM
            GLL = GLL + RJLLLL(M,4)*D2(M,1) - RKLLLL(M,4)*D2(M,1)
     &                + RJLLSS(M,4)*D2(M,3)
            GSL = GSL                       - RKSLSL(M,4)*D2(M,2)
            GSS = GSS + RJSSSS(M,4)*D2(M,3) - RKSSSS(M,4)*D2(M,3)
     &                + RJSSLL(M,4)*D2(M,1)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO COULOMB MATRIX
          G22(IBAS,JBAS) = GLL
          G22(KBAS,JBAS) = GSL
          G22(JBAS,KBAS) = GSL
          G22(KBAS,LBAS) = GSS
C
        ENDDO
      ENDDO
C
      RETURN
      END

