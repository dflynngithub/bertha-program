      SUBROUTINE BREIT0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           BBBBBBB  RRRRRRR  EEEEEEEE IIII TTTTTTTT 000000            C
C           BB    BB RR    RR EE        II     TT   00   000           C
C           BB    BB RR    RR EE        II     TT   00  0000           C
C           BBBBBBB  RR    RR EEEEEE    II     TT   00 00 00           C
C           BB    BB RRRRRRR  EE        II     TT   0000  00           C
C           BB    BB RR    RR EE        II     TT   000   00           C
C           BBBBBBB  RR    RR EEEEEEEE IIII    TT    000000            C
C                                                                      C
C -------------------------------------------------------------------- C
C  BREIT0 CONSTRUCTS THE ATOMIC BREIT MATRIX FROM RADIAL DIRECT AND    C
C  EXCHANGE INTEGRALS AND A MEAN-FIELD CHARGE DENSITY.                 C
C**********************************************************************C
      PARAMETER(MBS=44,MBD=2*MBS,MB2=MBS*MBS,MKP=13,MNU=MKP+1,
     &                                                      MAB=2*MNU+6)
C
      DIMENSION RKLLSS(MB2,4),RKSLLS(MB2,4),RKSSLL(MB2,4),RMSLLS(MB2,4)
C
      COMMON/ATMB/B11(MBD,MBD),B21(MBD,MBD),B12(MBD,MBD),B22(MBD,MBD)
      COMMON/ATMD/DLL1(MB2),DSL1(MB2),DSS1(MB2),DLS1(MB2),
     &            DLL2(MB2),DSL2(MB2),DSS2(MB2),DLS2(MB2)
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
      COMMON/XABE/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/XBIJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ
      COMMON/XCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
C
C     INITIALISE BREIT MATRIX
      DO IBAS=1,MBD
        DO JBAS=1,MBD
          B11(IBAS,JBAS) = 0.0D0
          B21(IBAS,JBAS) = 0.0D0
          B12(IBAS,JBAS) = 0.0D0
          B22(IBAS,JBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     COEFFICIENTS FOR INTEGRAL ASSEMBLY
      C3 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+3)
      C5 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+5)
      C7 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+7)
C
      V1 = 1.0D0
      V2 = 2.0D0
      V4 = 4.0D0
C
      TI = DFLOAT(2*LQNA+1)
      TJ = DFLOAT(2*LQNA+1)
      TK = DFLOAT(2*LQNB+1)
      TL = DFLOAT(2*LQNB+1)
C
      T0000 = 1.0D0
      T1000 = TI
      T0100 = TJ
      T0010 = TK
      T0001 = TL
      T1100 = TI*TJ
      T1001 = TI*TL
      T0011 = TK*TL
C
C     EVALUATE CLOSED-SHELL BREIT INTERACTION ANGULAR INTEGRALS
      CALL ANGBRT0
C
C     BASIS SET INTERMEDIATES FOR THE KL-PAIRS
      CALL KLSET
C
C     ITERATE OVER ALL MATRIX ELEMENTS
      DO IBAS=1,NBASA
        DO JBAS=1,NBASA
C
C         GAUSSIAN EXPONENTS FOR THIS PAIR
          EI = EXLA(IBAS)
          EJ = EXLA(JBAS)
C
C         BASIS SET INTERMEDIATES FOR THE IJ-PAIRS
          CALL IJSET
C
C         GENERATE BATCH OF RADIAL INTEGRALS (J AND K MATRICES)
          CALL BII0(RKLLSS,RKSLLS,RKSSLL,RMSLLS)
C
C         SMALL-COMPONENT MATRIX ADDRESSES
          KBAS = IBAS + NBASA
          LBAS = JBAS + NBASA
C
C    (22) KQNA < 0 AND KQNB < 0 CONTRIBUTIONS (CANNOT SKIP)
C
C         INITIALISE COUNTERS
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0
C
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,1)*DSS2(M)
            BSL = BSL + RKSLLS(M,1)*DLS2(M) + RMSLLS(M,1)*DLS2(M)
            BSS = BSS + RKSSLL(M,1)*DLL2(M)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B22(IBAS,JBAS) = BLL
          B22(KBAS,JBAS) = BSL
          B22(JBAS,KBAS) = BSL
          B22(KBAS,LBAS) = BSS
C
C    (21) KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNB.EQ.0) GOTO 200
C
C         INITIALISE COUNTERS
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0
C
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,2)*DSS1(M)
            BSL = BSL + RKSLLS(M,2)*DLS1(M) + RMSLLS(M,2)*DLS1(M)
            BSS = BSS + RKSSLL(M,2)*DLL1(M)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B21(IBAS,JBAS) = BLL
          B21(KBAS,JBAS) = BSL
          B21(JBAS,KBAS) = BSL
          B21(KBAS,LBAS) = BSS
C
200       CONTINUE
C
C    (12) KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0) GOTO 300
C
C         INITIALISE COUNTERS
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0
C
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,3)*DSS2(M)
            BSL = BSL + RKSLLS(M,3)*DLS2(M) + RMSLLS(M,3)*DLS2(M)
            BSS = BSS + RKSSLL(M,3)*DLL2(M)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B12(IBAS,JBAS) = BLL
          B12(KBAS,JBAS) = BSL
          B12(JBAS,KBAS) = BSL
          B12(KBAS,LBAS) = BSS
C
300       CONTINUE
C
C    (11) KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 400
C
C         INITIALISE COUNTERS
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0
C
C         SUM OVER MEAN FIELD CONTRIBUTIONS FOR THIS BASIS PAIR
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,4)*DSS1(M)
            BSL = BSL + RKSLLS(M,4)*DLS1(M) + RMSLLS(M,4)*DLS1(M)
            BSS = BSS + RKSSLL(M,4)*DLL1(M)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B11(IBAS,JBAS) = BLL
          B11(KBAS,JBAS) = BSL
          B11(JBAS,KBAS) = BSL
          B11(KBAS,LBAS) = BSS
C
400       CONTINUE
C
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE BII0(RKLLSS,RKSLLS,RKSSLL,RMSLLS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                      BBBBBBB IIII IIII 000000                        C
C                      BB    BB II   II 00   000                       C
C                      BB    BB II   II 00  0000                       C
C                      BBBBBBB  II   II 00 00 00                       C
C                      BB    BB II   II 0000  00                       C
C                      BB    BB II   II 000   00                       C
C                      BBBBBBB IIII IIII 000000                        C
C                                                                      C
C -------------------------------------------------------------------- C
C  BII0 EVALUATES A DIRECT AND EXCHANGE BATCH OF BREIT INTERACTION     C
C  INTEGRALS OF ALL COMPONENT LABEL COMBINATIONS L AND S IN THE ATOMIC C
C  (RELATIVISTIC) SCF PROCEDURE.                                       C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C    RKLLSS(M,N) - EXCHANGE BII OVERLAP {LL,SS}                        C
C    RKSLLS(M,N) - EXCHANGE BII OVERLAP {SL,SL}                        C
C    RKSSLL(M,N) - EXCHANGE BII OVERLAP {SS,LL}                        C
C    RMSLLS(M,N) - SEMI-RANGE BII OVERLAP {SL,SL}                      C
C -------------------------------------------------------------------- C
C    N=1 - KQN(A)<0, KQN(B)<0 (TYPICAL LABEL 22)                       C
C    N=2 - KQN(A)<0, KQN(B)>0 (TYPICAL LABEL 12)                       C
C    N=3 - KQN(A)>0, KQN(B)<0 (TYPICAL LABEL 21)                       C
C    N=4 - KQN(A)>0, KQN(B)>0 (TYPICAL LABEL 11)                       C
C**********************************************************************C
      PARAMETER(MBS=44,MB2=MBS*MBS,MKP=13,MNU=MKP+1,MAB=2*MNU+6)
C
      DIMENSION XK(MB2,2),IAA(2),IBB(2)
      DIMENSION RTIK0(MBS),RTJL0(MBS),PTIK0(MBS),PTJL0(MBS)
      DIMENSION BETA(MB2),BTA1(MB2),XROOT(MB2),BIN(MB2),TRM(MB2)
      DIMENSION BXU(MB2,-MAB:MAB,-MAB:MAB),BXL(MB2,-MAB:MAB,-MAB:MAB)
      DIMENSION RKLLSS(MB2,4),RKSLLS(MB2,4),RKSSLL(MB2,4),RMSLLS(MB2,4)
C
      COMMON/ANGL/BK(MNU,4),ELL(MNU,4),ESS(MNU,4),ESL(MNU,4),GSL(MNU,4)
      COMMON/XABE/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/XBIJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ
      COMMON/XBIK/EIK(MB2,-MAB:MAB),IKIND(MB2)
      COMMON/XBJL/EJL(MB2,-MAB:MAB),JLIND(MB2)
      COMMON/XBKL/EKL(MB2,-MAB:MAB),RNKL(MB2,4),EK(MB2),EL(MB2)
      COMMON/XCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
C
C**********************************************************************C
C     PREPARATION OF EXPONENT POWERS AND BETA FUNCTION ARGUMENTS.      C
C**********************************************************************C
C
C     TENSOR ORDER LIMITS
      NUI = NUS(1)
      NUF = NUS(NUNUM)
C
C     GAUSSIAN EXPONENT FOR (IBAS,JBAS) OVERLAP
      TIJ0 = EI+EJ
C
C     INITIALISE THE ARRAY XK(M,2) FOR INCOMPLETE BETA FUNCTION ARGS.
      DO M=1,MAXM
        TKL0  = EK(M)+EL(M)
        TIJKL = TIJ0+TKL0
        XK(M,1) = (EI+EK(M))/TIJKL
        XK(M,2) = (EJ+EL(M))/TIJKL
      ENDDO
C
C     LOWEST EXPONENT POWER
      IPOWER = LQNA+LQNB-NUF
C
C     A BLOCK OF BASIS EXPONENT PRODUCTS
      DO KBAS=1,NBASB
        RTIK0(KBAS) = DSQRT(EI+EXLB(KBAS))
        RTJL0(KBAS) = DSQRT(EJ+EXLB(KBAS))
        PTIK0(KBAS) = RTIK0(KBAS)**(-IPOWER)
        PTJL0(KBAS) = RTJL0(KBAS)**(-IPOWER)
      ENDDO
C
C     CALCULATE A FULL SET OF EXPONENT OVERLAPS FOR EXCHANGE
      DO M=1,MAXM
        RTIK = RTIK0(IKIND(M))
        RTJL = RTJL0(JLIND(M))
        EIK(M,-NUF) = PTIK0(IKIND(M))
        EJL(M,-NUF) = PTJL0(JLIND(M))
        DO IPOW=-NUF+1,NUF+4
          EIK(M,IPOW) = EIK(M,IPOW-1)/RTIK
          EJL(M,IPOW) = EJL(M,IPOW-1)/RTJL
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INLINE BETA FUNCTION CODE FOR EXCHANGE TERMS                     C
C**********************************************************************C
C
C     NUMBER OF LEVELS NEEDED TO ACCOUNT FOR ALL EXCHANGE INTEGRALS
      NVALS = (NUF-NUI)/2+2
C
C     LOOP OVER ORDER FOR FIRST INDEX
      DO NX=1,NVALS
        IAA(1) = LQNA+LQNB+NUI+2*NX
        IAA(2) = LQNA+LQNB+NUI+2*NX
C
C       LOOP OVER ORDER FOR SECOND INDEX
        DO NY=1,NVALS
          IBB(1) = LQNA+LQNB-NUF+2*NY-1
          IBB(2) = LQNA+LQNB-NUF+2*NY-1
C
C         LOOP OVER BETA INTEGRAL TYPE
          DO IBETA=1,2
            IA = (IAA(IBETA)-1)/2
            IB =  IBB(IBETA)   /2
C
C           CASE 1: IB > 1
            IF(IB.GT.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BTA1(M) = (DSQRT(XK(M,IBETA))**IX)/X
              ENDDO
              RA = X
              RB = DFLOAT(1-IB)
              RC = 1.0D0+X
              RD = 1.0D0
              RCD = RC*RD
              FCT = RA*RB/RCD
              DO M=1,MAXM
                TRM(M) = FCT*XK(M,IBETA)
                BIN(M) = 1.0D0+TRM(M)
              ENDDO
              RA = RA+1.0D0
              RB = RB+1.0D0
              RC = RC+1.0D0
              RD = RD+1.0D0
              DO J=2,IB-1
                RCD = RC*RD
                FCT = RA*RB/RCD
                DO M=1,MAXM
                  TRM(M) = FCT*TRM(M)*XK(M,IBETA)
                  BIN(M) = BIN(M)+TRM(M)
                ENDDO
                RA = RA+1.0D0
                RB = RB+1.0D0
                RC = RC+1.0D0
                RD = RD+1.0D0
              ENDDO
              DO M=1,MAXM
                BETA(M) = BTA1(M)*BIN(M)
              ENDDO
C
C           CASE 2: IB = 1
            ELSEIF(IB.EQ.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BETA(M) = (DSQRT(XK(M,IBETA))**IX)/X
              ENDDO
C
C           CASE 3: IB = 0
            ELSEIF(IB.EQ.0) THEN
              DO M=1,MAXM
                XROOT(M) = DSQRT(XK(M,IBETA))
                DEN      = (1.0D0-XROOT(M))
                RAT      = (1.0D0+XROOT(M))/DEN
                BTA1(M)  = DLOG(RAT)
                BIN(M)   = 1.0D0
                TRM(M)   = XK(M,IBETA)
              ENDDO
C
C             CASE A: IA > 1
              IF(IA.GT.1) THEN
                DO K=2,IA
                  KK = 2*K-1
                  RK = DFLOAT(KK)
                  X  = 1.0D0/RK
                  DO M=1,MAXM
                    BIN(M) = BIN(M)+X*TRM(M)
                    TRM(M) = TRM(M)*XK(M,IBETA)
                  ENDDO
                ENDDO
                DO M=1,MAXM
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)*BIN(M)
                ENDDO
C             CASE B: IA = 1
              ELSEIF(IA.EQ.1) THEN
                DO M=1,MAXM
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)
                ENDDO
C             CASE C: IA = 0
              ELSE
                DO M=1,MAXM
                  BETA(M) = BTA1(M)
                ENDDO
              ENDIF
C
            ENDIF
C
C           SORT THE BETA INTEGRAL INTO THE UPPER/LOWER ARRAYS
            IF(IBETA.EQ.1) THEN
              DO M=1,MAXM
                BXU(M, NUI+2*NX  ,-NUF+2*NY-1) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXM
                BXL(M, NUI+2*NX  ,-NUF+2*NY-1) = BETA(M)
              ENDDO
            ENDIF
C
C         END LOOPS OVER IBETA AND INDICES NX, NY
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     RADIAL INTEGRALS OVER SPINORS (SEPARATED BY TENSOR ORDER)        C
C**********************************************************************C
C
C     EMPTY COUNTER ARRAYS FOR DIRECT AND EXCHANGE INTEGRALS
      DO M=1,MAXM
        DO N=1,4
          RKLLSS(M,N) = 0.0D0
          RKSLLS(M,N) = 0.0D0
          RKSSLL(M,N) = 0.0D0
          RMSLLS(M,N) = 0.0D0
        ENDDO
      ENDDO
C
      E0000 = 1.0D0
      E1000 = EI
      E0100 = EJ
      E1100 = EI*EJ
C
C     INITIATE LOOP OVER K,L BASIS FUNCTIONS
      DO M=1,MAXM
C
C       MORE VALUE PREPARATION
        E0010 = EK(M)
        E0001 = EL(M)
        E1001 = EI*EL(M)
        E0011 = EK(M)*EL(M)
C
C**********************************************************************C
C       EXCHANGE INTEGRAL MATRICES: RKLLSS, RKSLLS, RKSSLL             C
C**********************************************************************C
C
C       LOOP OVER THE TENSOR ORDERS OF THE BREIT INTERACTION
        DO LTEN=1,NUNUM
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(LTEN)
C
C         TEMPORARY STORAGE OF RAW RK(LTEN,M)
          B21 = EIK(M,-NU+1)*EJL(M, NU+2)*BXL(M, NU+2,-NU+1)
     &        + EIK(M, NU+2)*EJL(M,-NU+1)*BXU(M, NU+2,-NU+1)
          B23 = EIK(M,-NU+3)*EJL(M, NU+2)*BXL(M, NU+2,-NU+3)
     &        + EIK(M, NU+4)*EJL(M,-NU+1)*BXU(M, NU+4,-NU+1)
          B41 = EIK(M,-NU+1)*EJL(M, NU+4)*BXL(M, NU+4,-NU+1)
     &        + EIK(M, NU+2)*EJL(M,-NU+3)*BXU(M, NU+2,-NU+3)
          B43 = EIK(M,-NU+3)*EJL(M, NU+4)*BXL(M, NU+4,-NU+3)
     &        + EIK(M, NU+4)*EJL(M,-NU+3)*BXU(M, NU+4,-NU+3)
C
C         FILL RK ARRAYS FOR THIS LQNA,LQNB BLOCK
C
C         KQNA < 0 AND KQNB < 0 CONTRIBUTIONS (CANNOT SKIP)
          RKLL = V4*T0000*E0011*C7*B43
          RKSL = V4*T0000*E1001*C7*B43
          RKSS = V4*T0000*E1100*C7*B43
C
          RKLLSS(M,1) = RKLLSS(M,1) + ELL(LTEN,1)*RKLL
          RKSLLS(M,1) = RKSLLS(M,1) + ESL(LTEN,1)*RKSL
          RKSSLL(M,1) = RKSSLL(M,1) + ESS(LTEN,1)*RKSS
C
C         KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNB.EQ.0) GOTO 202
          RKLL = V4*T0000*E0011*C7*B43 - V2*T0001*E0010*C5*B23
     &         - V2*T0010*E0001*C5*B41 + V1*T0011*E0000*C3*B21
          RKSL = V4*T0000*E1001*C7*B43 - V2*T0001*E1000*C5*B23
          RKSS = V4*T0000*E1100*C7*B43
C
          RKLLSS(M,2) = RKLLSS(M,2) + ELL(LTEN,2)*RKLL
          RKSLLS(M,2) = RKSLLS(M,2) + ESL(LTEN,2)*RKSL
          RKSSLL(M,2) = RKSSLL(M,2) + ESS(LTEN,2)*RKSS
202       CONTINUE
C
C         KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0) GOTO 203
          RKLL = V4*T0000*E0011*C7*B43
          RKSL = V4*T0000*E1001*C7*B43 - V2*T1000*E0001*C5*B41
          RKSS = V4*T0000*E1100*C7*B43 - V2*T0100*E1000*C5*B23
     &         - V2*T1000*E0100*C5*B41 + V1*T1100*E0000*C3*B21
C
          RKLLSS(M,3) = RKLLSS(M,3) + ELL(LTEN,3)*RKLL
          RKSLLS(M,3) = RKSLLS(M,3) + ESL(LTEN,3)*RKSL
          RKSSLL(M,3) = RKSSLL(M,3) + ESS(LTEN,3)*RKSS
203       CONTINUE
C
C         KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 204
          RKLL = V4*T0000*E0011*C7*B43 - V2*T0001*E0010*C5*B23
     &         - V2*T0010*E0001*C5*B41 + V1*T0011*E0000*C3*B21
          RKSL = V4*T0000*E1001*C7*B43 - V2*T1000*E0001*C5*B41
     &         - V2*T0001*E1000*C5*B23 + V1*T1001*E0000*C3*B21
          RKSS = V4*T0000*E1100*C7*B43 - V2*T0100*E1000*C5*B23
     &         - V2*T1000*E0100*C5*B41 + V1*T1100*E0000*C3*B21
C
          RKLLSS(M,4) = RKLLSS(M,4) + ELL(LTEN,4)*RKLL
          RKSLLS(M,4) = RKSLLS(M,4) + ESL(LTEN,4)*RKSL
          RKSSLL(M,4) = RKSSLL(M,4) + ESS(LTEN,4)*RKSS
204       CONTINUE
C
        ENDDO
C
C**********************************************************************C
C       HALF-RANGE EXCHANGE INTEGRAL MATRICES: RMSLLS                  C
C**********************************************************************C
C
C       LOOP OVER THE TENSOR ORDERS OF THE COULOMB INTERACTION
        DO LTEN=1,NUNUM
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(LTEN)
C
C         TEMPORARY STORAGE OF RAW RM(LTEN,M) FOR RMSLLS
          B21 = EIK(M,-NU+1)*EJL(M, NU+2)*BXL(M, NU+2,-NU+1)
     &        - EIK(M, NU+2)*EJL(M,-NU+1)*BXU(M, NU+2,-NU+1)
          B23 = EIK(M,-NU+3)*EJL(M, NU+2)*BXL(M, NU+2,-NU+3)
     &        - EIK(M, NU+4)*EJL(M,-NU+1)*BXU(M, NU+4,-NU+1)
          B41 = EIK(M,-NU+1)*EJL(M, NU+4)*BXL(M, NU+4,-NU+1)
     &        - EIK(M, NU+2)*EJL(M,-NU+3)*BXU(M, NU+2,-NU+3)
          B43 = EIK(M,-NU+3)*EJL(M, NU+4)*BXL(M, NU+4,-NU+3)
     &        - EIK(M, NU+4)*EJL(M,-NU+3)*BXU(M, NU+4,-NU+3)
C
C         FILL RM ARRAYS FOR THIS LQNA,LQNB BLOCK
C
C         KQNA < 0 AND KQNB < 0 CONTRIBUTIONS (CANNOT SKIP)
          RMSL = V4*T0000*E1001*C7*B43

          RMSLLS(M,1) = RMSLLS(M,1) + GSL(LTEN,1)*RMSL
C
C         KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNB.EQ.0) GOTO 302
          RMSL = V4*T0000*E1001*C7*B43 - V2*T0001*E1000*C5*B23
C
          RMSLLS(M,2) = RMSLLS(M,2) + GSL(LTEN,2)*RMSL
302       CONTINUE
C
C         KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0) GOTO 303
          RMSL = V4*T0000*E1001*C7*B43 - V2*T1000*E0001*C5*B41
C
          RMSLLS(M,3) = RMSLLS(M,3) + GSL(LTEN,3)*RMSL
303       CONTINUE
C
C         KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 304
          RMSL = V4*T0000*E1001*C7*B43 - V2*T1000*E0001*C5*B41
     &         - V2*T0001*E1000*C5*B23 + V1*T1001*E0000*C3*B21
C
          RMSLLS(M,4) = RMSLLS(M,4) + GSL(LTEN,4)*RMSL
304       CONTINUE
C
        ENDDO
C
      ENDDO
C
C**********************************************************************C
C     APPLY NORMALISATION FACTORS                                      C
C**********************************************************************C
C
      DO M=1,MAXM
        RNLLSS = RNIJ(1)*RNKL(M,3)
        RNSSLL = RNIJ(3)*RNKL(M,1)
        T0SLLS = RNIJ(2)*RNKL(M,4)
        DO N=1,4
          RKLLSS(M,N) = RNLLSS*RKLLSS(M,N)
          RKSSLL(M,N) = RNSSLL*RKSSLL(M,N)
          RKSLLS(M,N) = T0SLLS*RKSLLS(M,N)
          RMSLLS(M,N) = T0SLLS*RMSLLS(M,N)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ANGBRT0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        AA    NN    NN  GGGGGG  BBBBBBB  RRRRRRR TTTTTTTT 000000      C
C       AAAA   NNN   NN GG    GG BB    BB RR    RR   TT   00   000     C
C      AA  AA  NNNN  NN GG       BB    BB RR    RR   TT   00  0000     C
C     AA    AA NN NN NN GG       BBBBBBB  RR    RR   TT   00 00 00     C
C     AAAAAAAA NN  NNNN GG   GGG BB    BB RRRRRRR    TT   0000  00     C
C     AA    AA NN   NNN GG    GG BB    BB RR    RR   TT   000   00     C
C     AA    AA NN    NN  GGGGGG  BBBBBBB  RR    RR   TT    000000      C
C                                                                      C
C -------------------------------------------------------------------- C
C  ANGBRT0 EVALUATES ANGULAR COEFFICIENTS OF THE ATOMIC CLOSED SHELL   C
C  BREIT INTERACTION FOR ALL (K1,K2) VALUES IN THE MANIFOLD (L1,L2).   C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C    NUNUM - NUMBER OF NU VALUES THAT SATISFY PARITY RESTRICTION RULE. C
C    NUI - MINIMUM NU VALUE IN THIS MANIFOLD.                          C
C    NUF - MAXIMUM NU VALUE IN THIS MANIFOLD.                          C
C    ELL(MNU,4) - ANGULAR TERMS FOR CLOSED BREIT EK(JA,JB;LL) TERMS    C
C    ESS(MNU,4) - ANGULAR TERMS FOR CLOSED BREIT EK(JA,JB;SS) TERMS    C
C    ESL(MNU,4) - ANGULAR TERMS FOR CLOSED BREIT EK(JA,JB;SL) TERMS    C
C    GSL(MNU,4) - ANGULAR TERMS FOR CLOSED BREIT EK(JA,JB;SL) TERMS    C
C**********************************************************************C
      PARAMETER(MBS=44,MKP=13,MNU=MKP+1)
C
      DIMENSION SCOEF(4,2)
C
      COMMON/ANGL/BK(MNU,4),ELL(MNU,4),ESS(MNU,4),ESL(MNU,4),GSL(MNU,4)
      COMMON/XABE/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
C
C     INITIALISE COEFFICIENT ARRAYS
      DO LTEN=1,MNU
        DO N=1,4
          ELL(LTEN,N) = 0.0D0
          ESS(LTEN,N) = 0.0D0
          ESL(LTEN,N) = 0.0D0
          GSL(LTEN,N) = 0.0D0
        ENDDO
        NUS(LTEN) = 0
      ENDDO
      NUNUM = 0
C
      NUI = IABS(LQNA-LQNB  )
      NUF = IABS(LQNA+LQNB+1)
C
C     SPECIFY ALLOWED NU VALUES
      NUNUM = 0
      DO NU=NUI-1,NUF
        IF(NU.GE.0.AND.MOD(LQNA+LQNB+NU,2).EQ.1) THEN
          NUNUM = NUNUM+1
          NUS(NUNUM) = NU
        ENDIF
      ENDDO
C
      NUI = NUS(1)
      NUF = NUS(NUNUM)
C
C**********************************************************************C
C     (1) KQNA < 0 AND KQNB < 0   (CANNOT SKIP)                        C
C**********************************************************************C
C
C     KQNA AND KQNB
      KQNA =-LQNA-1
      KQNB =-LQNB-1
C
C     JQNA AND JQNB
      JQNA = 2*IABS(KQNA)-1
      JQNB = 2*IABS(KQNB)-1
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      LTEN = 1
      DO NU=NUI,NUF
C
C       SKIP IF THIS VALUE OF NU EXCEEDS TRIANGLE RULE FOR KQN BLOCK
        IF(NU.GT.(JQNA+JQNB)/2) GOTO 101
C
C       RAW SQUARED 3J-SYMBOL
        RAW = SYM3JSQ(JQNA,JQNB,NU)
C
C       DETERMINE PARITY OF COMBINATION LQNA,LQNB,NU
        IF(MOD(LQNA+LQNB+NU,2).EQ.0) THEN
          IPAR = 1
        ELSE
          IPAR = 0
        ENDIF
C
        IF(IPAR.EQ.0.AND.NU.NE.0) THEN
C       CASE 1: ANGULAR COEFFICIENTS OF ODD PARITY
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          RNU  = DFLOAT(NU*(NU+1))
          COEF = DFLOAT((KQNA+KQNB)*(KQNA+KQNB))/RNU
C
C         SAVE THIS TENSOR ORDER NU TO THE LTEN ELEMENT OF 'NUS'
          NUS(LTEN) = NU
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,1) = ELL(LTEN,1) + COEF*RAW
          ESS(LTEN,1) = ESS(LTEN,1) + COEF*RAW
          ESL(LTEN,1) = ESL(LTEN,1) + COEF*RAW
C
        ELSEIF(IPAR.EQ.1) THEN
C       CASE 2: ANGULAR COEFFICIENTS OF EVEN-PARITY
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          CALL BRCOEF0(SCOEF,KQNA,KQNB,NU)
C
C         SAVE THIS TENSOR ORDER NU TO THE LTEN ELEMENT OF 'NUS'
          NUS(LTEN) = NU-1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,1) = ELL(LTEN,1) - SCOEF(1,1)*RAW
          ESL(LTEN,1) = ESL(LTEN,1) - SCOEF(2,1)*RAW
          ESS(LTEN,1) = ESS(LTEN,1) - SCOEF(3,1)*RAW
          GSL(LTEN,1) = GSL(LTEN,1) - SCOEF(4,1)*RAW
C
C         STEP PARAMETER FOR ADDITIONS TO COEFFICIENTS
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SAVE THIS TENSOR ORDER NU TO THE LTEN ELEMENT OF 'NUS'
          NUS(LTEN) = NU+1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,1) = ELL(LTEN,1) - SCOEF(1,2)*RAW
          ESL(LTEN,1) = ESL(LTEN,1) - SCOEF(2,2)*RAW
          ESS(LTEN,1) = ESS(LTEN,1) - SCOEF(3,2)*RAW
          GSL(LTEN,1) = GSL(LTEN,1) - SCOEF(4,2)*RAW
C
        ENDIF
C
101     CONTINUE
C
      ENDDO
C
C     ADJUST NUNUM IF NECESSARY
      NUI   = MIN(NUI,NUS(   1))
      NUF   = MAX(NUF,NUS(LTEN))
      NUNUM = MAX(LTEN,NUNUM)
C
C     APPROPRIATE NUI GIVEN ODD SELECTION RULE
      NUI = MOD(LQNA+LQNB+1,2)
C
C**********************************************************************C
C     (2) KQNA < 0 AND KQNB > 0   (SKIP IF POSSIBLE)                   C
C**********************************************************************C
C
      IF(LQNB.EQ.0) GOTO 200
C
C     KQNA AND KQNB
      KQNA =-LQNA-1
      KQNB = LQNB
C
C     JQNA AND JQNB
      JQNA = 2*IABS(KQNA)-1
      JQNB = 2*IABS(KQNB)-1
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      LTEN = 1
      DO NU=NUI,NUF
C
C       SKIP IF THIS VALUE OF NU EXCEEDS TRIANGLE RULE FOR KQN BLOCK
        IF(NU.GT.(JQNA+JQNB)/2) GOTO 201
C
C       RAW SQUARED 3J-SYMBOL
        RAW = SYM3JSQ(JQNA,JQNB,NU)
C
C       TEST WHETHER PARITY OF 'LQNA+LQNB+NU' IS ODD OR EVEN
        IF(MOD(LQNA+LQNB+NU,2).EQ.0) THEN
          IPAR = 1
        ELSE
          IPAR = 0
        ENDIF
C
        IF(IPAR.EQ.0.AND.NU.NE.0) THEN
C       CASE 1: ANGULAR COEFFICIENTS OF ODD PARITY
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          RNU  = DFLOAT(NU*(NU+1))
          COEF = DFLOAT((KQNA+KQNB)*(KQNA+KQNB))/RNU
C
C         SAVE THIS TENSOR ORDER NU TO THE LTEN ELEMENT OF 'NUS'
          NUS(LTEN) = NU
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,2) = ELL(LTEN,2) + COEF*RAW
          ESS(LTEN,2) = ESS(LTEN,2) + COEF*RAW
          ESL(LTEN,2) = ESL(LTEN,2) + COEF*RAW
C
        ELSEIF(IPAR.EQ.1) THEN
C       CASE 2: ANGULAR COEFFICIENTS OF EVEN-PARITY
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          CALL BRCOEF0(SCOEF,KQNA,KQNB,NU)
C
C         SAVE THIS TENSOR ORDER NU TO THE LTEN ELEMENT OF 'NUS'
          NUS(LTEN) = NU-1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,2) = ELL(LTEN,2) - SCOEF(1,1)*RAW
          ESL(LTEN,2) = ESL(LTEN,2) - SCOEF(2,1)*RAW
          ESS(LTEN,2) = ESS(LTEN,2) - SCOEF(3,1)*RAW
          GSL(LTEN,2) = GSL(LTEN,2) - SCOEF(4,1)*RAW
C
C         STEP PARAMETER FOR ADDITIONS TO COEFFICIENTS
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SAVE THIS TENSOR ORDER NU TO THE LTEN ELEMENT OF 'NUS'
          NUS(LTEN) = NU+1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,2) = ELL(LTEN,2) - SCOEF(1,2)*RAW
          ESL(LTEN,2) = ESL(LTEN,2) - SCOEF(2,2)*RAW
          ESS(LTEN,2) = ESS(LTEN,2) - SCOEF(3,2)*RAW
          GSL(LTEN,2) = GSL(LTEN,2) - SCOEF(4,2)*RAW
C
        ENDIF
C
201     CONTINUE
C
      ENDDO
C
C     ADJUST NUNUM IF NECESSARY
      NUI   = MIN(NUI,NUS(   1))
      NUF   = MAX(NUF,NUS(LTEN))
      NUNUM = MAX(LTEN,NUNUM)
C
200   CONTINUE
C
C**********************************************************************C
C     (3) KQNA > 0 AND KQNB < 0   (SKIP IF POSSIBLE)                   C
C**********************************************************************C
C
      IF(LQNA.EQ.0) GOTO 300
C
C     KQNA AND KQNB
      KQNA = LQNA
      KQNB =-LQNB-1
C
C     JQNA AND JQNB
      JQNA = 2*IABS(KQNA)-1
      JQNB = 2*IABS(KQNB)-1
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      LTEN = 1
      DO NU=NUI,NUF
C
C       SKIP IF THIS VALUE OF NU EXCEEDS TRIANGLE RULE FOR KQN BLOCK
        IF(NU.GT.(JQNA+JQNB)/2) GOTO 301
C
C       RAW SQUARED 3J-SYMBOL
        RAW = SYM3JSQ(JQNA,JQNB,NU)
C
C       TEST WHETHER PARITY OF 'LQNA+LQNB+NU' IS ODD OR EVEN
        IF(MOD(LQNA+LQNB+NU,2).EQ.0) THEN
          IPAR = 1
        ELSE
          IPAR = 0
        ENDIF
C
        IF(IPAR.EQ.0.AND.NU.NE.0) THEN
C       CASE 1: ANGULAR COEFFICIENTS OF ODD PARITY
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          RNU  = DFLOAT(NU*(NU+1))
          COEF = DFLOAT((KQNA+KQNB)*(KQNA+KQNB))/RNU
C
C         SAVE THIS TENSOR ORDER NU TO THE LTEN ELEMENT OF 'NUS'
          NUS(LTEN) = NU
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,3) = ELL(LTEN,3) + COEF*RAW
          ESS(LTEN,3) = ESS(LTEN,3) + COEF*RAW
          ESL(LTEN,3) = ESL(LTEN,3) + COEF*RAW
C
        ELSEIF(IPAR.EQ.1) THEN
C       CASE 2: ANGULAR COEFFICIENTS OF EVEN-PARITY
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          CALL BRCOEF0(SCOEF,KQNA,KQNB,NU)
C
C         SAVE THIS TENSOR ORDER NU TO THE LTEN ELEMENT OF 'NUS'
          NUS(LTEN) = NU-1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,3) = ELL(LTEN,3) - SCOEF(1,1)*RAW
          ESL(LTEN,3) = ESL(LTEN,3) - SCOEF(2,1)*RAW
          ESS(LTEN,3) = ESS(LTEN,3) - SCOEF(3,1)*RAW
          GSL(LTEN,3) = GSL(LTEN,3) - SCOEF(4,1)*RAW
C
C         STEP PARAMETER FOR ADDITIONS TO COEFFICIENTS
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SAVE THIS TENSOR ORDER NU TO THE LTEN ELEMENT OF 'NUS'
          NUS(LTEN) = NU+1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,3) = ELL(LTEN,3) - SCOEF(1,2)*RAW
          ESL(LTEN,3) = ESL(LTEN,3) - SCOEF(2,2)*RAW
          ESS(LTEN,3) = ESS(LTEN,3) - SCOEF(3,2)*RAW
          GSL(LTEN,3) = GSL(LTEN,3) - SCOEF(4,2)*RAW
C
        ENDIF
C
301     CONTINUE
C
      ENDDO
C
C     ADJUST NUNUM IF NECESSARY
      NUI   = MIN(NUI,NUS(   1))
      NUF   = MAX(NUF,NUS(LTEN))
      NUNUM = MAX(LTEN,NUNUM)
C
300   CONTINUE
C
C**********************************************************************C
C     (4) KQNA > 0 AND KQNB > 0   (SKIP IF POSSIBLE)                   C
C**********************************************************************C
C
      IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 400
C
C     KQNA AND KQNB
      KQNA = LQNA
      KQNB = LQNB
C
C     JQNA AND JQNB
      JQNA = 2*IABS(KQNA)-1
      JQNB = 2*IABS(KQNB)-1
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      LTEN = 1
      DO NU=NUI,NUF
C
C       SKIP IF THIS VALUE OF NU EXCEEDS TRIANGLE RULE FOR KQN BLOCK
        IF(NU.GT.(JQNA+JQNB)/2) GOTO 401
C
C       RAW SQUARED 3J-SYMBOL
        RAW = SYM3JSQ(JQNA,JQNB,NU)
C
C       TEST WHETHER PARITY OF 'LQNA+LQNB+NU' IS ODD OR EVEN
        IF(MOD(LQNA+LQNB+NU,2).EQ.0) THEN
          IPAR = 1
        ELSE
          IPAR = 0
        ENDIF
C
        IF(IPAR.EQ.0.AND.NU.NE.0) THEN
C       CASE 1: ANGULAR COEFFICIENTS OF ODD PARITY
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          RNU  = DFLOAT(NU*(NU+1))
          COEF = DFLOAT((KQNA+KQNB)*(KQNA+KQNB))/RNU
C
C         SAVE THIS TENSOR ORDER NU TO THE LTEN ELEMENT OF 'NUS'
          NUS(LTEN) = NU
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,4) = ELL(LTEN,4) + COEF*RAW
          ESS(LTEN,4) = ESS(LTEN,4) + COEF*RAW
          ESL(LTEN,4) = ESL(LTEN,4) + COEF*RAW
C
        ELSEIF(IPAR.EQ.1) THEN
C       CASE 2: ANGULAR COEFFICIENTS OF EVEN-PARITY
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          CALL BRCOEF0(SCOEF,KQNA,KQNB,NU)
C
C         SAVE THIS TENSOR ORDER NU TO THE LTEN ELEMENT OF 'NUS'
          NUS(LTEN) = NU-1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,4) = ELL(LTEN,4) - SCOEF(1,1)*RAW
          ESL(LTEN,4) = ESL(LTEN,4) - SCOEF(2,1)*RAW
          ESS(LTEN,4) = ESS(LTEN,4) - SCOEF(3,1)*RAW
          GSL(LTEN,4) = GSL(LTEN,4) - SCOEF(4,1)*RAW
C
C         STEP PARAMETER FOR ADDITIONS TO COEFFICIENTS
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SAVE THIS TENSOR ORDER NU TO THE LTEN ELEMENT OF 'NUS'
          NUS(LTEN) = NU+1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,4) = ELL(LTEN,4) - SCOEF(1,2)*RAW
          ESL(LTEN,4) = ESL(LTEN,4) - SCOEF(2,2)*RAW
          ESS(LTEN,4) = ESS(LTEN,4) - SCOEF(3,2)*RAW
          GSL(LTEN,4) = GSL(LTEN,4) - SCOEF(4,2)*RAW
C
        ENDIF
C
401     CONTINUE
C
      ENDDO
C
C     ADJUST NUNUM IF NECESSARY
      NUI   = MIN(NUI,NUS(   1))
      NUF   = MAX(NUF,NUS(LTEN))
      NUNUM = MAX(LTEN,NUNUM)
C
400   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE BRCOEF0(SCOEFF,KA,KB,NU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    BBBBBBB  RRRRRRR   CCCCCC   OOOOOO  EEEEEEEE FFFFFFFF 000000      C
C    BB    BB RR    RR CC    CC OO    OO EE       FF      00   000     C
C    BB    BB RR    RR CC       OO    OO EE       FF      00  0000     C
C    BBBBBBB  RR    RR CC       OO    OO EEEEEE   FFFFFF  00 00 00     C
C    BB    BB RRRRRRR  CC       OO    OO EE       FF      0000  00     C
C    BB    BB RR    RR CC    CC OO    OO EE       FF      000   00     C
C    BBBBBBB  RR    RR  CCCCCC   OOOOOO  EEEEEEEE FF       000000      C
C                                                                      C
C -------------------------------------------------------------------- C
C  BRCOEF0 EVALUATES THE INTERMEDIATE COEFFICIENTS OF THE BREIT        C
C  INTERACTION FOR CLOSED SHELLS (TABLE 3 OF GRANT AND PYPER 1976).    C
C**********************************************************************C
      DIMENSION SCOEFF(4,2)
C
      RU  = DFLOAT(NU)
      RK  = DFLOAT(KB-KA)
C
      IF(NU.GT.0) THEN
        RM = RU-1.0D0
        B1 = (RM+2.0D0)/(2.0D0*   (2.0D0*RM+1.0D0))
        C1 =-(RM-1.0D0)/(2.0D0*RU*(2.0D0*RM+1.0D0))
        SCOEFF(1,1) =-(RK+RU)*(C1*RK+B1)
        SCOEFF(2,1) =  B1*RU - C1*RK*RK
        SCOEFF(3,1) =-(RK-RU)*(C1*RK-B1)
        SCOEFF(4,1) =- RK    *(C1*RU-B1)
      ELSE
        SCOEFF(1,1) = 0.0D0
        SCOEFF(2,1) = 0.0D0
        SCOEFF(3,1) = 0.0D0
        SCOEFF(4,1) = 0.0D0
      ENDIF
      IF(NU+1.GT.1) THEN
        RP = RU+1.0D0
        B2 = (RP-1.0D0)/(2.0D0*   (2.0D0*RP+1.0D0))
        C2 = (RP+2.0D0)/(2.0D0*RP*(2.0D0*RP+1.0D0))
        SCOEFF(1,2) =-(RK-RP)*(C2*RK+B2)
        SCOEFF(2,2) =- B2*RP - C2*RK*RK 
        SCOEFF(3,2) =-(RK+RP)*(C2*RK-B2)
        SCOEFF(4,2) =- RK    *(C2*RP+B2)
      ELSE
        SCOEFF(1,2) = 0.0D0
        SCOEFF(2,2) = 0.0D0
        SCOEFF(3,2) = 0.0D0
        SCOEFF(4,2) = 0.0D0
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE IJSET
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C               IIII     JJJJ SSSSSS  EEEEEEEE TTTTTTTT                C
C                II       JJ SS    SS EE          TT                   C
C                II       JJ SS       EE          TT                   C
C                II       JJ  SSSSSS  EEEEEE      TT                   C
C                II       JJ       SS EE          TT                   C
C                II JJ    JJ SS    SS EE          TT                   C
C               IIII JJJJJJ   SSSSSS  EEEEEEEE    TT                   C
C                                                                      C
C -------------------------------------------------------------------- C
C  IJSET GENERATES BASIS SET INTERMEDIATES FOR IJ-PAIRS TO BE USED     C
C  IN THE CONSTRUCTION OF ATOMIC TWO-ELECTRON INTEGRALS.               C
C**********************************************************************C
      PARAMETER(MBS=44,MKP=13,MNU=MKP+1,MAB=2*MNU+6)
C
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
      COMMON/XABE/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/XBIJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ
C
      DATA TWOLOG/6.93147180559945309D-1/
C
C     NORMALISATION CONSTANTS FOR EXPONENTS EI AND EJ
      RL = DFLOAT(LQNA)
      G1 = TWOLOG-GAMLOG(2*LQNA+3)
      G2 = TWOLOG-GAMLOG(2*LQNA+5)
      R1 = RL+1.5D0
      R2 = RL+0.5D0
C
      ELOG = DLOG(2.0D0*EI)
      RNLI = DEXP(0.5D0*(G1+R1*ELOG))
      RNSI = DEXP(0.5D0*(G2+R2*ELOG))
C
      ELOG = DLOG(2.0D0*EJ)
      RNLJ = DEXP(0.5D0*(G1+R1*ELOG))
      RNSJ = DEXP(0.5D0*(G2+R2*ELOG))
C
C     NORMALISATION PAIRS
      RNIJ(1) = RNLI*RNLJ
      RNIJ(2) = RNSI*RNLJ
      RNIJ(3) = RNSI*RNSJ
      RNIJ(4) = RNLI*RNSJ
C
C     POWERS OF THE EXPONENT SUM
      EIJ0 = EI+EJ
      EIJR = DSQRT(EIJ0)
      EIJA = EIJ0**(-LQNA)
      DO N=0,5
        EIJ(N) = EIJA
        EIJA   = EIJA/EIJR
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE KLSET
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             KK    KK LL       SSSSSS  EEEEEEEE TTTTTTTT              C
C             KK   KK  LL      SS    SS EE          TT                 C
C             KK  KK   LL      SS       EE          TT                 C
C             KKKKK    LL       SSSSSS  EEEEEE      TT                 C
C             KK  KK   LL            SS EE          TT                 C
C             KK   KK  LL      SS    SS EE          TT                 C
C             KK    KK LLLLLLLL SSSSSS  EEEEEEEE    TT                 C
C                                                                      C
C -------------------------------------------------------------------- C
C  IJSET GENERATES BASIS SET INTERMEDIATES FOR KL-PAIRS TO BE USED     C
C  IN THE CONSTRUCTION OF ATOMIC TWO-ELECTRON INTEGRALS.               C
C**********************************************************************C
      PARAMETER(MBS=44,MB2=MBS*MBS,MKP=13,MNU=MKP+1,MAB=2*MNU+6)
C
      COMMON/XABE/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/XBIK/EIK(MB2,-MAB:MAB),IKIND(MB2)
      COMMON/XBJL/EJL(MB2,-MAB:MAB),JLIND(MB2)
      COMMON/XBKL/EKL(MB2,-MAB:MAB),RNKL(MB2,4),EK(MB2),EL(MB2)
C
C     GENERATE INDICES AND EXPONENT COMBINATIONS
      M = 0
      DO KBAS=1,NBASB
        EK0 = EXLB(KBAS)
        DO LBAS=1, NBASB
          M   = M+1
          EL0 = EXLB(LBAS)
          IKIND(M) = KBAS
          JLIND(M) = LBAS
          EK(M)    = EK0
          EL(M)    = EL0
          EKL0     = EK0+EL0
          EKLR     = DSQRT(EKL0)
          EKPW     = EKL0**LQNB
          EKLA     = 1.0D0/EKPW
          DO N=0,5
            EKL(M,N) = EKLA
            EKLA     = EKLA/EKLR
          ENDDO
        ENDDO
      ENDDO
C
C     NORMALISATION CONSTANTS
      CALL RNORM0(RNKL,EXLB,NBASB,LQNB)
C
      RETURN
      END
C
C
      SUBROUTINE RNORM0(RN,EXL,NBAS,LQN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       RRRRRRR  NN    NN  OOOOOO  RRRRRRR  MM       MM  000000        C
C       RR    RR NNN   NN OO    OO RR    RR MMM     MMM 00   000       C
C       RR    RR NNNN  NN OO    OO RR    RR MMMM   MMMM 00  0000       C
C       RR    RR NN NN NN OO    OO RR    RR MM MM MM MM 00 00 00       C
C       RRRRRRR  NN  NNNN OO    OO RRRRRRR  MM  MMM  MM 0000  00       C
C       RR    RR NN   NNN OO    OO RR    RR MM   M   MM 000   00       C
C       RR    RR NN    NN  OOOOOO  RR    RR MM       MM  000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  RNORM0 EVALUATES NORMALISATION CONSTANTS OF ALL VARIETIES.          C
C**********************************************************************C
      PARAMETER(MBS=44,MB2=MBS*MBS)
C
      DIMENSION RN(MB2,4),EXL(MBS),RNL(MBS),RNS(MBS)
C
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
C
      DATA TWOLOG/6.93147180559945309D-1/
C
      RL = DFLOAT(LQN)
      G1 = TWOLOG-GAMLOG(2*LQN+3)
      G2 = TWOLOG-GAMLOG(2*LQN+5)
      R1 = RL+1.5D0
      R2 = RL+0.5D0
      DO IBAS=1,NBAS
        ELOG      = DLOG(2.0D0*EXL(IBAS))
        RNL(IBAS) = DEXP(0.5D0*(G1+R1*ELOG))
        RNS(IBAS) = DEXP(0.5D0*(G2+R2*ELOG))
      ENDDO
C
C     RN(M,1) ARE THE LL NORMALISATION CONSTANTS
C     RN(M,2) ARE THE SL NORMALISATION CONSTANTS
C     RN(M,3) ARE THE SS NORMALISATION CONSTANTS
C     RN(M,4) ARE THE LS NORMALISATION CONSTANTS
C
      M = 0
      DO IBAS=1,NBAS
        DO JBAS=1,NBAS
          M = M+1
          RN(M,1) = RNL(IBAS)*RNL(JBAS)
          RN(M,2) = RNS(IBAS)*RNL(JBAS)
          RN(M,3) = RNS(IBAS)*RNS(JBAS)
          RN(M,4) = RNL(IBAS)*RNS(JBAS)
        ENDDO
      ENDDO
C
      RETURN
      END

