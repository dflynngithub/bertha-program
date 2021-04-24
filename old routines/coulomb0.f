      SUBROUTINE COULOMB0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    CCCCCC   OOOOOO  UU    UU LL      MM       MM BBBBBBB   000000    C
C   CC    CC OO    OO UU    UU LL      MMM     MMM BB    BB 00   000   C
C   CC       OO    OO UU    UU LL      MMMM   MMMM BB    BB 00  0000   C
C   CC       OO    OO UU    UU LL      MM MM MM MM BBBBBBB  00 00 00   C
C   CC       OO    OO UU    UU LL      MM  MMM  MM BB    BB 0000  00   C
C   CC    CC OO    OO UU    UU LL      MM   M   MM BB    BB 000   00   C
C    CCCCCC   OOOOOO   UUUUUU  LLLLLLL MM       MM BBBBBBB   000000    C
C                                                                      C
C -------------------------------------------------------------------- C
C  COULOMB0 CONSTRUCTS THE ATOMIC COULOMB MATRIX FROM RADIAL DIRECT    C
C  AND EXCHANGE INTEGRALS AND A MEAN-FIELD CHARGE DENSITY.             C
C**********************************************************************C
      PARAMETER(MBS=44,MBD=2*MBS,MB2=MBS*MBS,MKP=13,MNU=MKP+1,
     &                                                      MAB=2*MNU+6)
C
      CHARACTER*4 HMLTN
C
      DIMENSION RJLLLL(MB2,4),RJSSSS(MB2,4),RJLLSS(MB2,4),RJSSLL(MB2,4),
     &          RKLLLL(MB2,4),RKSSSS(MB2,4),RKSLSL(MB2,4)
C
      COMMON/ATMC/G11(MBD,MBD),G21(MBD,MBD),G12(MBD,MBD),G22(MBD,MBD)
      COMMON/ATMD/DLL1(MB2),DSL1(MB2),DSS1(MB2),DLS1(MB2),
     &            DLL2(MB2),DSL2(MB2),DSS2(MB2),DLS2(MB2)
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/XABE/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/XBIJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ
      COMMON/XCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
C
C     INITIALISE COULOMB MATRIX
      DO IBAS=1,MBD
        DO JBAS=1,MBD
          G11(IBAS,JBAS) = 0.0D0
          G21(IBAS,JBAS) = 0.0D0
          G12(IBAS,JBAS) = 0.0D0
          G22(IBAS,JBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     COEFFICIENTS FOR INTEGRAL ASSEMBLY
      C1 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+1)
      C3 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+3)
      C5 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+5)
      C7 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+7)
      C9 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+9)
C
      V1 = 1.0D0
      V2 = 2.0D0
      V4 = 4.0D0
      V8 = 8.0D0
      VS = 1.6D1
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
      T1010 = TI*TK
      T1001 = TI*TL
      T0110 = TJ*TK
      T0101 = TJ*TL
      T0011 = TK*TL
      T1110 = TI*TJ*TK
      T1101 = TI*TJ*TL
      T1011 = TI*TK*TL
      T0111 = TJ*TK*TL
      T1111 = TI*TJ*TK*TL
C
C     EVALUATE CLOSED-SHELL ELECTRON REPULSION ANGULAR INTEGRALS
      CALL ANGCLM0
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
          CALL ERI0(RJLLLL,RJSSSS,RJLLSS,RJSSLL,RKLLLL,RKSSSS,RKSLSL)
C
          IF(HMLTN.EQ.'NORL') THEN
C         NON-RELATIVISTIC HAMILTONIAN
C
C           INITIALISE COUNTER
            GLL = 0.0D0
C
C           BUILD THE FOCK MATRIX
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,1)*DLL2(M) - RKLLLL(M,1)*DLL2(M)
            ENDDO
C
C           TRANSFER COUNTER VALUE TO COULOMB MATRIX
            G22(IBAS,JBAS) = GLL
C
          ELSE
C         RELATIVISTIC HAMILTONIAN
C
C           SMALL-COMPONENT MATRIX ADDRESSES
            KBAS = IBAS + NBASA
            LBAS = JBAS + NBASA
C
C    (22)   KQNA < 0 AND KQNB < 0  CONTRIBUTIONS (CANNOT SKIP)
C
C           INITIALISE COUNTERS
            GLL = 0.0D0
            GSL = 0.0D0
            GSS = 0.0D0
C
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,1)*DLL2(M) - RKLLLL(M,1)*DLL2(M)
     &                  + RJLLSS(M,1)*DSS2(M)
              GSL = GSL                       - RKSLSL(M,1)*DSL2(M)
              GSS = GSS + RJSSSS(M,1)*DSS2(M) - RKSSSS(M,1)*DSS2(M)
     &                  + RJSSLL(M,1)*DLL2(M)
            ENDDO
C
C           TRANSFER COUNTER VALUES TO COULOMB MATRIX
            G22(IBAS,JBAS) = GLL
            G22(KBAS,JBAS) = GSL
            G22(JBAS,KBAS) = GSL
            G22(KBAS,LBAS) = GSS
C
C    (21)   KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNB.EQ.0) GOTO 200
C
C           INITIALISE COUNTERS
            GLL = 0.0D0
            GSL = 0.0D0
            GSS = 0.0D0
C
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,2)*DLL1(M) - RKLLLL(M,2)*DLL1(M)
     &                  + RJLLSS(M,2)*DSS1(M)
              GSL = GSL                       - RKSLSL(M,2)*DSL1(M)
              GSS = GSS + RJSSSS(M,2)*DSS1(M) - RKSSSS(M,2)*DSS1(M)
     &                  + RJSSLL(M,2)*DLL1(M)
            ENDDO
C
C           TRANSFER COUNTER VALUES TO COULOMB MATRIX
            G21(IBAS,JBAS) = GLL
            G21(KBAS,JBAS) = GSL
            G21(JBAS,KBAS) = GSL
            G21(KBAS,LBAS) = GSS
C
200         CONTINUE
C
C    (12)   KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNA.EQ.0) GOTO 300
C
C           INITIALISE COUNTERS
            GLL = 0.0D0
            GSL = 0.0D0
            GSS = 0.0D0
C
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,3)*DLL2(M) - RKLLLL(M,3)*DLL2(M)
     &                  + RJLLSS(M,3)*DSS2(M)
              GSL = GSL                       - RKSLSL(M,3)*DSL2(M)
              GSS = GSS + RJSSSS(M,3)*DSS2(M) - RKSSSS(M,3)*DSS2(M)
     &                  + RJSSLL(M,3)*DLL2(M)
            ENDDO
C
C           TRANSFER COUNTER VALUES TO COULOMB MATRIX
            G12(IBAS,JBAS) = GLL
            G12(KBAS,JBAS) = GSL
            G12(JBAS,KBAS) = GSL
            G12(KBAS,LBAS) = GSS
C
300         CONTINUE
C
C    (11)   KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 400
C
C           INITIALISE COUNTERS
            GLL = 0.0D0
            GSL = 0.0D0
            GSS = 0.0D0
C
C           SUM OVER MEAN FIELD CONTRIBUTIONS FOR THIS BASIS PAIR
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,4)*DLL1(M) - RKLLLL(M,4)*DLL1(M)
     &                  + RJLLSS(M,4)*DSS1(M)
              GSL = GSL                       - RKSLSL(M,4)*DSL1(M)
              GSS = GSS + RJSSSS(M,4)*DSS1(M) - RKSSSS(M,4)*DSS1(M)
     &                  + RJSSLL(M,4)*DLL1(M)
            ENDDO
C
C           TRANSFER COUNTER VALUES TO COULOMB MATRIX
            G11(IBAS,JBAS) = GLL
            G11(KBAS,JBAS) = GSL
            G11(JBAS,KBAS) = GSL
            G11(KBAS,LBAS) = GSS
C
400         CONTINUE
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
      SUBROUTINE ERI0(RJLLLL,RJSSSS,RJLLSS,RJSSLL,RKLLLL,RKSSSS,RKSLSL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                    EEEEEEEE RRRRRRR  IIII 000000                     C
C                    EE       RR    RR  II 00   000                    C
C                    EE       RR    RR  II 00  0000                    C
C                    EEEEEE   RR    RR  II 00 00 00                    C
C                    EE       RRRRRRR   II 0000  00                    C
C                    EE       RR    RR  II 000   00                    C
C                    EEEEEEEE RR    RR IIII 000000                     C
C                                                                      C
C -------------------------------------------------------------------- C
C  ERI0 EVALUATES A DIRECT AND EXCHANGE BATCH OF ELECTRON REPULSION    C
C  INTEGRALS OF ALL COMPONENT LABEL COMBINATIONS L AND S IN THE ATOMIC C
C  SCF PROCEDURE FOR A USER-SPECIFIED HAMILTONIAN.                     C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C    RJLLLL(M,N) - DIRECT ERI OVERLAP {LL,LL}                          C
C    RJSSSS(M,N) - DIRECT ERI OVERLAP {SS,SS}                          C
C    RJLLSS(M,N) - DIRECT ERI OVERLAP {LL,SS}                          C
C    RJSSLL(M,N) - DIRECT ERI OVERLAP {SS,LL}                          C
C    RKLLLL(M,N) - EXCHANGE ERI OVERLAP {LL,LL}                        C
C    RKSSSS(M,N) - EXCHANGE ERI OVERLAP {SS,SS}                        C
C    RKSLSL(M,N) - EXCHANGE ERI OVERLAP {SL,SL}                        C
C -------------------------------------------------------------------- C
C    N=1 - KQN(A)<0, KQN(B)<0 (TYPICAL LABEL 22)                       C
C    N=2 - KQN(A)<0, KQN(B)>0 (TYPICAL LABEL 12)                       C
C    N=3 - KQN(A)>0, KQN(B)<0 (TYPICAL LABEL 21)                       C
C    N=4 - KQN(A)>0, KQN(B)>0 (TYPICAL LABEL 11)                       C
C**********************************************************************C
      PARAMETER(MBS=44,MB2=MBS*MBS,MKP=13,MNU=MKP+1,MAB=2*MNU+6)
C
      CHARACTER*4 HMLTN
C
      DIMENSION XJ(MB2,2),XK(MB2,2),IAA(2),IBB(2)
      DIMENSION BETA(MB2),BTA1(MB2),XROOT(MB2),BIN(MB2),TRM(MB2)
      DIMENSION BDU(MB2,-MAB:MAB,-MAB:MAB),BDL(MB2,-MAB:MAB,-MAB:MAB)
      DIMENSION BXU(MB2,-MAB:MAB,-MAB:MAB),BXL(MB2,-MAB:MAB,-MAB:MAB)
C
      DIMENSION RTIK0(MBS),RTJL0(MBS),PTIK0(MBS),PTJL0(MBS)
      DIMENSION RJLLLL(MB2,4),RJSSSS(MB2,4),RJLLSS(MB2,4),RJSSLL(MB2,4),
     &          RKLLLL(MB2,4),RKSSSS(MB2,4),RKSLSL(MB2,4)
C
      COMMON/ANGL/BK(MNU,4),ELL(MNU,4),ESS(MNU,4),ESL(MNU,4),GSL(MNU,4)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
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
C     INITIALISE THE ARRAY XJ(M,2) FOR INCOMPLETE BETA FUNCTION ARGS.
      DO M=1,MAXM
        TKL0    = EK(M)+EL(M)
        TIJKL   = TIJ0+TKL0
        XJ(M,1) = TIJ0/TIJKL
        XJ(M,2) = TKL0/TIJKL
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
        RTIK  = RTIK0(IKIND(M))
        RTJL  = RTJL0(JLIND(M))
        EIK(M,-NUF) = PTIK0(IKIND(M))
        EJL(M,-NUF) = PTJL0(JLIND(M))
        DO IPOW=-NUF+1,NUF+5
          EIK(M,IPOW) = EIK(M,IPOW-1)/RTIK
          EJL(M,IPOW) = EJL(M,IPOW-1)/RTJL
        ENDDO
      ENDDO
C
C     CALCULATE LIST OF BETA FUNCTION ARGUMENTS, XK(MB2,2) (Z AND Z')
      DO M=1,MAXM
        TKL0    = EK(M)+EL(M)
        TIJKL   = EI+EJ+TKL0
        XK(M,1) = (EI+EK(M))/TIJKL
        XK(M,2) = (EJ+EL(M))/TIJKL
      ENDDO
C
C**********************************************************************C
C     INLINE BETA FUNCTION CODE FOR DIRECT TERMS                       C
C**********************************************************************C
C
C     NUMBER OF LEVELS NEEDED TO ACCOUNT FOR ALL DIRECT INTEGRALS
      NVALS = 3
C
C     LOOP OVER ORDER FOR FIRST INDEX
      DO NX=1,NVALS
C
C       TWICE THE ACTUAL FAMILY VALUE
        IAA(1) = 2*LQNA+2*NX-1
        IAA(2) = 2*LQNB+2*NX-1
C
C       LOOP OVER ORDER FOR SECOND INDEX
        DO NY=1,NVALS
          IBB(1) = 2*LQNB+2*NY-2
          IBB(2) = 2*LQNA+2*NY-2
C
C         LOOP OVER BETA INTEGRAL TYPE
          DO IBETA=1,2
            IA =(IAA(IBETA)-1)/2
            IB = IBB(IBETA)   /2
C
C           BEGIN CONDITIONAL STATEMENT OVER IB VALUES
C           CASE 1: IB > 1
            IF(IB.GT.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BTA1(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
              RA  = X
              RB  = DFLOAT(1-IB)
              RC  = X+1.0D0
              RD  = 1.0D0
              RCD = RC*RD
              FCT = RA*RB/RCD
              DO M=1,MAXM
                TRM(M) = FCT*XJ(M,IBETA)
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
                  TRM(M) = FCT*TRM(M)*XJ(M,IBETA)
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
                BETA(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
C
C           CASE 3: IB = 0
            ELSEIF(IB.EQ.0) THEN
              DO M=1,MAXM
                XROOT(M) = DSQRT(XJ(M,IBETA))
                DEN      =  1.0D0-XROOT(M)
                RAT      = (1.0D0+XROOT(M))/DEN
                BTA1(M)  = DLOG(RAT)
                BIN(M)   = 1.0D0
                TRM(M)   = XJ(M,IBETA)
              ENDDO
              IF(IA.GT.1) THEN
                DO K=2,IA
                  KK = 2*K-1
                  RK = DFLOAT(KK)
                  X  = 1.0D0/RK
                  DO M=1,MAXM
                    BIN(M) = BIN(M)+X*TRM(M)
                    TRM(M) = TRM(M)*XJ(M,IBETA)
                  ENDDO
                ENDDO
                DO M=1,MAXM
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)*BIN(M)
                ENDDO
              ELSEIF(IA.EQ.1) THEN
                DO M=1,MAXM
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)
                ENDDO
              ELSE
                DO M=1,MAXM
                  BETA(M) = BTA1(M)
                ENDDO
              ENDIF
C
            ENDIF
C
C           SORT THE BETA INTEGRAL INTO THE CORRECT ARRAY
            IF(IBETA.EQ.1) THEN
              DO M=1,MAXM
                BDU(M,2*NX-1,2*NY-2) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXM
                BDL(M,2*NX-1,2*NY-2) = BETA(M)
              ENDDO
            ENDIF
C
C         END LOOPS OVER IBETA AND INDICES NX, NY
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INLINE BETA FUNCTION CODE FOR EXCHANGE TERMS                     C
C**********************************************************************C
C
C     NUMBER OF LEVELS NEEDED TO ACCOUNT FOR ALL EXCHANGE INTEGRALS
      NVALS = (NUF-NUI)/2+3
C
C     LOOP OVER ORDER FOR FIRST INDEX
      DO NX=1,NVALS
        IAA(1) = LQNA+LQNB+NUI+2*NX-1
        IAA(2) = LQNA+LQNB+NUI+2*NX-1
C
C       LOOP OVER ORDER FOR SECOND INDEX
        DO NY=1,NVALS
          IBB(1) = LQNA+LQNB-NUF+2*NY-2
          IBB(2) = LQNA+LQNB-NUF+2*NY-2
C
C         LOOP OVER BETA INTEGRAL TYPE
          DO IBETA=1,2
            IA = (IAA(IBETA)-1)/2
            IB =  IBB(IBETA)   /2
C
C           BEGIN CONDITIONAL STATEMENT OVER IB VALUES
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
              ELSEIF(IA.EQ.1) THEN
                DO M=1,MAXM
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)
                ENDDO
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
                BXU(M, NUI+2*NX-1,-NUF+2*NY-2) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXM
                BXL(M, NUI+2*NX-1,-NUF+2*NY-2) = BETA(M)
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
          RJLLLL(M,N) = 0.0D0
          RJLLSS(M,N) = 0.0D0
          RJSSLL(M,N) = 0.0D0
          RJSSSS(M,N) = 0.0D0
          RKLLLL(M,N) = 0.0D0
          RKSLSL(M,N) = 0.0D0
          RKSSSS(M,N) = 0.0D0
        ENDDO
      ENDDO
C
C     VALUE PREPARATION
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
        E1010 = EI*EK(M)
        E1001 = EI*EL(M)
        E0110 = EJ*EK(M)
        E0101 = EJ*EL(M)
        E0011 = EK(M)*EL(M)
        E1110 = EI*EJ*EK(M)
        E1101 = EI*EJ*EL(M)
        E1011 = EI*EK(M)*EL(M)
        E0111 = EJ*EK(M)*EL(M)
        E1111 = EI*EJ*EK(M)*EL(M)
C
C**********************************************************************C
C       DIRECT INTEGRAL MATRICES: RJSSSS, RJLLSS, RJSSLL, RJLLLL       C
C**********************************************************************C
C
        IF(HMLTN.EQ.'NORL') THEN
C       NON-RELATIVISTIC HAMILTONIAN
C
          B32 = EIJ(2)*EKL(M,3)*BDL(M,3,2) + EIJ(3)*EKL(M,2)*BDU(M,3,2)
C
          RJLLLL(M,1) = V1*T0000*E0000*C5*B32
C
        ELSE
C       RELATIVISTIC HAMILTONIAN
C
C         TEMPORARY STORAGE OF RAW RJ(1,M) (LTEN=1 BECAUSE NU=0 ONLY)
          B10 = EIJ(0)*EKL(M,1)*BDL(M,1,0) + EIJ(1)*EKL(M,0)*BDU(M,1,0)
          B12 = EIJ(2)*EKL(M,1)*BDL(M,1,2) + EIJ(3)*EKL(M,0)*BDU(M,3,0)
          B14 = EIJ(4)*EKL(M,1)*BDL(M,1,4) + EIJ(5)*EKL(M,0)*BDU(M,5,0)
          B30 = EIJ(0)*EKL(M,3)*BDL(M,3,0) + EIJ(1)*EKL(M,2)*BDU(M,1,2)
          B32 = EIJ(2)*EKL(M,3)*BDL(M,3,2) + EIJ(3)*EKL(M,2)*BDU(M,3,2)
          B34 = EIJ(4)*EKL(M,3)*BDL(M,3,4) + EIJ(5)*EKL(M,2)*BDU(M,5,2)
          B50 = EIJ(0)*EKL(M,5)*BDL(M,5,0) + EIJ(1)*EKL(M,4)*BDU(M,1,4)
          B52 = EIJ(2)*EKL(M,5)*BDL(M,5,2) + EIJ(3)*EKL(M,4)*BDU(M,3,4)
          B54 = EIJ(4)*EKL(M,5)*BDL(M,5,4) + EIJ(5)*EKL(M,4)*BDU(M,5,4)
C
C         FILL RJ ARRAYS FOR THIS LQNA,LQNB BLOCK
C
C         KQNA < 0 AND KQNB < 0 CONTRIBUTIONS (CANNOT SKIP)
          RJLLLL(M,1) = V1*T0000*E0000*C5*B32
          RJLLSS(M,1) = V4*T0000*E0011*C7*B52
          RJSSLL(M,1) = V4*T0000*E1100*C7*B34
          RJSSSS(M,1) = VS*T0000*E1111*C9*B54
C
C         KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNB.EQ.0) GOTO 102
          RJLLLL(M,2) = V1*T0000*E0000*C5*B32
          RJLLSS(M,2) = V4*T0000*E0011*C7*B52
     &                - V2*T0001*E0010*C5*B32 - V2*T0010*E0001*C5*B32
     &                + V1*T0011*E0000*C3*B12
          RJSSLL(M,2) = V4*T0000*E1100*C7*B34
          RJSSSS(M,2) = VS*T0000*E1111*C9*B54
     &                - V8*T0001*E1110*C7*B34 - V8*T0010*E1101*C7*B34
     &                + V4*T0011*E1100*C5*B14
102       CONTINUE
C
C         KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0) GOTO 103
          RJLLLL(M,3) = V1*T0000*E0000*C5*B32
          RJLLSS(M,3) = V4*T0000*E0011*C7*B52
          RJSSLL(M,3) = V4*T0000*E1100*C7*B34
     &                - V2*T0100*E1000*C5*B32 - V2*T1000*E0100*C5*B32
     &                + V1*T1100*E0000*C3*B30
          RJSSSS(M,3) = VS*T0000*E1111*C9*B54
     &                - V8*T0100*E1011*C7*B52 - V8*T1000*E0111*C7*B52
     &                + V4*T1100*E0011*C5*B50
103       CONTINUE
C
C         KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 104
          RJLLLL(M,4) = V1*T0000*E0000*C5*B32
          RJLLSS(M,4) = V4*T0000*E0011*C7*B52
     &                - V2*T0001*E0010*C5*B32 - V2*T0010*E0001*C5*B32
     &                + V1*T0011*E0000*C3*B12
          RJSSLL(M,4) = V4*T0000*E1100*C7*B34
     &                - V2*T0100*E1000*C5*B32 - V2*T1000*E0100*C5*B32
     &                + V1*T1100*E0000*C3*B30
          RJSSSS(M,4) = VS*T0000*E1111*C9*B54
     &                - V8*T0001*E1110*C7*B34 - V8*T0010*E1101*C7*B34
     &                - V8*T0100*E1011*C7*B52 - V8*T1000*E0111*C7*B52
     &                + V4*T1100*E0011*C5*B50 + V4*T0011*E1100*C5*B14
     &                + V4*T0110*E1001*C5*B32 + V4*T1001*E0110*C5*B32
     &                + V4*T1010*E0101*C5*B32 + V4*T0101*E1010*C5*B32
     &                - V2*T0111*E1000*C3*B12 - V2*T1011*E0100*C3*B12
     &                - V2*T1101*E0010*C3*B30 - V2*T1110*E0001*C3*B30
     &                + V1*T1111*E0000*C1*B10
104       CONTINUE
C
        ENDIF
C
C**********************************************************************C
C       EXCHANGE INTEGRAL MATRICES: RKSSSS, RKSLSL, RKLLLL             C
C**********************************************************************C
C
C       LOOP OVER THE TENSOR ORDERS OF THE COULOMB INTERACTION
        DO LTEN=1,NUNUM
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(LTEN)
C
          IF(HMLTN.EQ.'NORL') THEN
C         NON-RELATIVISTIC HAMILTONIAN
C
            B32 = EIK(M,-NU+2)*EJL(M, NU+3)*BXL(M, NU+3,-NU+2)
     &          + EIK(M, NU+3)*EJL(M,-NU+2)*BXU(M, NU+3,-NU+2)
C
            RKLL = V1*T0000*E0000*C5*B32
C
            RKLLLL(M,1) = RKLLLL(M,1) + BK(LTEN,1)*RKLL
C
          ELSE
C         RELATIVISTIC HAMILTONIAN
C
C           TEMPORARY STORAGE OF RAW RK(LTEN,M)
            B10 = EIK(M,-NU  )*EJL(M, NU+1)*BXL(M, NU+1,-NU  )
     &          + EIK(M, NU+1)*EJL(M,-NU  )*BXU(M, NU+1,-NU  )
            B12 = EIK(M,-NU+2)*EJL(M, NU+1)*BXL(M, NU+1,-NU+2)
     &          + EIK(M, NU+3)*EJL(M,-NU  )*BXU(M, NU+3,-NU  )
            B14 = EIK(M,-NU+4)*EJL(M, NU+1)*BXL(M, NU+1,-NU+4)
     &          + EIK(M, NU+5)*EJL(M,-NU  )*BXU(M, NU+5,-NU  )
            B30 = EIK(M,-NU  )*EJL(M, NU+3)*BXL(M, NU+3,-NU  )
     &          + EIK(M, NU+1)*EJL(M,-NU+2)*BXU(M, NU+1,-NU+2)
            B32 = EIK(M,-NU+2)*EJL(M, NU+3)*BXL(M, NU+3,-NU+2)
     &          + EIK(M, NU+3)*EJL(M,-NU+2)*BXU(M, NU+3,-NU+2)
            B34 = EIK(M,-NU+4)*EJL(M, NU+3)*BXL(M, NU+3,-NU+4)
     &          + EIK(M, NU+5)*EJL(M,-NU+2)*BXU(M, NU+5,-NU+2)
            B50 = EIK(M,-NU  )*EJL(M, NU+5)*BXL(M, NU+5,-NU  )
     &          + EIK(M, NU+1)*EJL(M,-NU+4)*BXU(M, NU+1,-NU+4)
            B52 = EIK(M,-NU+2)*EJL(M, NU+5)*BXL(M, NU+5,-NU+2)
     &          + EIK(M, NU+3)*EJL(M,-NU+4)*BXU(M, NU+3,-NU+4)
            B54 = EIK(M,-NU+4)*EJL(M, NU+5)*BXL(M, NU+5,-NU+4)
     &          + EIK(M, NU+5)*EJL(M,-NU+4)*BXU(M, NU+5,-NU+4)
C
C           KQNA < 0 AND KQNB < 0 CONTRIBUTIONS (CANNOT SKIP)
            RKLL = V1*T0000*E0000*C5*B32
            RKSL = V4*T0000*E1010*C7*B34
            RKSS = VS*T0000*E1111*C9*B54
C
            RKLLLL(M,1) = RKLLLL(M,1) + BK(LTEN,1)*RKLL
            RKSLSL(M,1) = RKSLSL(M,1) + BK(LTEN,1)*RKSL
            RKSSSS(M,1) = RKSSSS(M,1) + BK(LTEN,1)*RKSS
C
C           KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNB.EQ.0) GOTO 202
            RKLL = V1*T0000*E0000*C5*B32
            RKSL = V4*T0000*E1010*C7*B34 - V2*T0010*E1000*C5*B32
            RKSS = VS*T0000*E1111*C9*B54 - V8*T0001*E1110*C7*B34
     &           - V8*T0010*E1101*C7*B52 + V4*T0011*E1100*C5*B32
C
            RKLLLL(M,2) = RKLLLL(M,2) + BK(LTEN,2)*RKLL
            RKSLSL(M,2) = RKSLSL(M,2) + BK(LTEN,2)*RKSL
            RKSSSS(M,2) = RKSSSS(M,2) + BK(LTEN,2)*RKSS
202         CONTINUE
C
C           KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNA.EQ.0) GOTO 203
            RKLL = V1*T0000*E0000*C5*B32
            RKSL = V4*T0000*E1010*C7*B34 - V2*T0100*E0010*C5*B32
            RKSS = VS*T0000*E1111*C9*B54 - V8*T0100*E1011*C7*B34
     &           - V8*T1000*E0111*C7*B52 + V4*T1100*E0011*C5*B32
C
            RKLLLL(M,3) = RKLLLL(M,3) + BK(LTEN,3)*RKLL
            RKSLSL(M,3) = RKSLSL(M,3) + BK(LTEN,3)*RKSL
            RKSSSS(M,3) = RKSSSS(M,3) + BK(LTEN,3)*RKSS
203         CONTINUE
C
C           KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 204
            RKLL = V1*T0000*E0000*C5*B32
            RKSL = V4*T0000*E1010*C7*B34 - V2*T1000*E0010*C5*B32
     &           - V2*T0010*E1000*C5*B32 + V1*T1010*E0000*C3*B30
            RKSS = VS*T0000*E1111*C9*B54
     &           - V8*T0001*E1110*C7*B34 - V8*T0010*E1101*C7*B52
     &           - V8*T0100*E1011*C7*B34 - V8*T1000*E0111*C7*B52
     &           + V4*T1100*E0011*C5*B32 + V4*T0011*E1100*C5*B32
     &           + V4*T1001*E0110*C5*B32 + V4*T0110*E1001*C5*B32
     &           + V4*T0101*E1010*C5*B14 + V4*T1010*E0101*C5*B50
     &           - V2*T1101*E0010*C3*B12 - V2*T0111*E1000*C3*B12
     &           - V2*T1110*E0001*C3*B30 - V2*T1011*E0100*C3*B30
     &           + V1*T1111*E0000*C1*B10
C
            RKLLLL(M,4) = RKLLLL(M,4) + BK(LTEN,4)*RKLL
            RKSLSL(M,4) = RKSLSL(M,4) + BK(LTEN,4)*RKSL
            RKSSSS(M,4) = RKSSSS(M,4) + BK(LTEN,4)*RKSS
204         CONTINUE
C
          ENDIF
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     APPLY NORMALISATION FACTORS                                      C
C**********************************************************************C
C
      IF(HMLTN.EQ.'NORL') THEN
C     NON-RELATIVISTIC HAMILTONIAN
C
        DO M=1,MAXM
          RNLLLL = RNIJ(1)*RNKL(M,1)
          RJLLLL(M,1) = RJLLLL(M,1)*RNLLLL
          RKLLLL(M,1) = RKLLLL(M,1)*RNLLLL
        ENDDO
C
      ELSE
C     RELATIVISTIC HAMILTONIAN
C
        DO M=1,MAXM
          RNLLLL = RNIJ(1)*RNKL(M,1)
          RNLLSS = RNIJ(1)*RNKL(M,3)
          RNSLSL = RNIJ(2)*RNKL(M,2)
          RNSSLL = RNIJ(3)*RNKL(M,1)
          RNSSSS = RNIJ(3)*RNKL(M,3)
          DO N=1,4
            RJLLLL(M,N) = RJLLLL(M,N)*RNLLLL
            RJLLSS(M,N) = RJLLSS(M,N)*RNLLSS
            RJSSLL(M,N) = RJSSLL(M,N)*RNSSLL
            RJSSSS(M,N) = RJSSSS(M,N)*RNSSSS
            RKLLLL(M,N) = RKLLLL(M,N)*RNLLLL
            RKSLSL(M,N) = RKSLSL(M,N)*RNSLSL
            RKSSSS(M,N) = RKSSSS(M,N)*RNSSSS
          ENDDO
        ENDDO
C
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE ANGCLM0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     AA    NN    NN  GGGGGG   CCCCCC  LL       MM       MM  000000    C
C    AAAA   NNN   NN GG    GG CC    CC LL       MMM     MMM 00   000   C
C   AA  AA  NNNN  NN GG       CC       LL       MMMM   MMMM 00  0000   C
C  AA    AA NN NN NN GG       CC       LL       MM MM MM MM 00 00 00   C
C  AAAAAAAA NN  NNNN GG   GGG CC       LL       MM  MMM  MM 0000  00   C
C  AA    AA NN   NNN GG    GG CC    CC LL       MM   M   MM 000   00   C
C  AA    AA NN    NN  GGGGGG   CCCCCC  LLLLLLLL MM       MM  000000    C
C                                                                      C
C -------------------------------------------------------------------- C
C  ANGCLM0 EVALUATES THE ANGULAR COEFFICIENTS OF THE COULOMB           C
C  INTERACTIONS FOR CLOSED SHELLS IN THE (L1,L2) MANIFOLD.             C
C**********************************************************************C
      PARAMETER(MBS=44,MKP=13,MNU=MKP+1)
C
      CHARACTER*4 HMLTN
C
      COMMON/ANGL/BK(MNU,4),ELL(MNU,4),ESS(MNU,4),ESL(MNU,4),GSL(MNU,4)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/XABE/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
C
C     CALCULATE KQNA AND 2*JQNA VALUES FROM LQNA
      KLA =-LQNA-1
      KRA = LQNA
      JLA = 2*IABS(KLA)-1
      JRA = 2*IABS(KRA)-1
C
C     CALCULATE KQNB AND 2*JQNB VALUES FROM LQNB
      KLB =-LQNB-1
      KRB = LQNB
      JLB = 2*IABS(KLB)-1
      JRB = 2*IABS(KRB)-1
C
C     START AND END PARAMETERS FROM TRIANGLE RULE
      NUI = IABS(LQNA-LQNB)
      NUF = LQNA+LQNB+1
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      LTEN = 0
      DO NU=NUI,NUF
C
C       TEST WHETHER 'LQNA+LQNB+NU' ODD OR EVEN
        IF(MOD(LQNA+LQNB+NU,2).EQ.0) THEN
          IPARAB = 1
        ELSE
          IPARAB = 0
        ENDIF
C
        IF(IPARAB.EQ.1) THEN
C       ONLY ANGULAR COEFFICIENTS OF EVEN PARITY ARE NON-ZERO
C
C         SAVE THIS TENSOR ORDER
          LTEN      = LTEN+1
          NUS(LTEN) = NU
C
          IF(HMLTN.EQ.'NORL') THEN
            BK(LTEN,1) = 0.5D0*ABC000(LQNA,LQNB,NU)
          ELSE
            BK(LTEN,1) = SYM3JSQ(JLA,JLB,NU)
            BK(LTEN,2) = SYM3JSQ(JLA,JRB,NU)
            BK(LTEN,3) = SYM3JSQ(JRA,JLB,NU)
            BK(LTEN,4) = SYM3JSQ(JRA,JRB,NU)
          ENDIF
        ENDIF
C
      ENDDO
C
C     NUMBER OF SURVIVING TENSOR ORDERS
      NUNUM = LTEN
C
      RETURN
      END
C
C
      FUNCTION ABC000(L1,L2,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           AA    BBBBBBB   CCCCCC   000000   000000   000000          C
C          AAAA   BB    BB CC    CC 00   000 00   000 00   000         C
C         AA  AA  BB    BB CC       00  0000 00  0000 00  0000         C
C        AA    AA BBBBBBB  CC       00 00 00 00 00 00 00 00 00         C
C        AAAAAAAA BB    BB CC       0000  00 0000  00 0000  00         C
C        AA    AA BB    BB CC    CC 000   00 000   00 000   00         C
C        AA    AA BBBBBBB   CCCCCC   000000   000000   000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  ABC000 EVALUATES THE NON-RELATIVISTIC 3-J SYMBOL FOR ATOMIC COULOMB C
C  ANGULAR COEFFICIENT ROUTINES, TAKEN FROM BRINK AND SATCHLER.        C
C  L1,L2, AND K MUST BE EQUAL TO THE ACTUAL (INTEGER) ANGULAR MOMENTA  C
C  OF THE ELECTRON AND PHOTON.                                         C
C**********************************************************************C
      COMMON/FCTS/RFACT(0:20),SFACT(0:20)
C
C     TRIANGLE INEQUALITY RESTRICTIONS
      IF(K.LT.IABS(L1-L2).OR.K.GT.(L1+L2)) THEN
        ABC000 = 0.0D0
        RETURN
      ENDIF
      LLK = L1+L2+K
C
C     PARITY SELECTION RULE
      IF((LLK/2)*2.NE.LLK) THEN
        ABC000 = 0.0D0
        RETURN
      ENDIF
C
      RF1 = RFACT(  L1+L2-K   )
      RF2 = RFACT(- L1+L2+K   )
      RF3 = RFACT(  L1-L2+K   )
      RF4 = RFACT(  L1+L2+K +1)
      RF5 = RFACT(( L1+L2+K)/2)
      RF6 = RFACT(( L1+L2-K)/2)
      RF7 = RFACT(( L1-L2+K)/2)
      RF8 = RFACT((-L1+L2+K)/2)
C
      T1 = RF1*RF2*RF3
      T2 = T1/RF4
      T3 = RF6*RF7*RF8
      T4 = RF5/T3
C
      ABC000 = T2*T4*T4
C
      RETURN
      END
C
C
      FUNCTION SYM3JSQ(J1,J2,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   SSSSSS  YY    YY MM       MM  333333       JJJJ SSSSSS   QQQQQQ    C
C  SS    SS YY    YY MMM     MMM 33    33       JJ SS    SS QQ    QQ   C
C  SS       YY    YY MMMM   MMMM       33       JJ SS       QQ    QQ   C
C   SSSSSS   YY  YY  MM MM MM MM    3333        JJ  SSSSSS  QQ    QQ   C
C        SS   YYYY   MM  MMM  MM       33       JJ       SS QQ   QQQ   C
C  SS    SS    YY    MM   M   MM 33    33 JJ    JJ SS    SS QQ    QQ   C
C   SSSSSS     YY    MM       MM  333333   JJJJJJ   SSSSSS   QQQQQQ Q  C
C                                                                      C
C -------------------------------------------------------------------- C
C  SYM3JSQ EVALUATES THE SQUARE OF A 3-J SYMBOL,   /  j   K   j' \^2   C
C  WHERE j = J1/2 AND j' = J2/2, FOR THE           \-1/2  0  1/2 /     C
C  COULOMB/BREIT ANGULAR COEFFICIENT ROUTINES.                         C
C**********************************************************************C
      COMMON/FCTS/RFACT(0:20),SFACT(0:20)
C
C     TRIANGLE RULE RESTRICTIONS
      IF(K.LT.IABS((J1-J2)/2).OR.K.GT.(J1+J2)/2) THEN
        SYM3JSQ = 0.0D0
        RETURN
      ELSEIF(J1.LE.0.OR.J2.LE.0) THEN
        SYM3JSQ = 0.0D0
        RETURN
      ENDIF
C
C     VARIABLE WHICH DEPENDS ON PARITY OF ARGUMENTS
      JJK = (J1+J2)/2 + K
      IF((JJK/2)*2.EQ.JJK) THEN
        M = K
      ELSE
        M = K+1
      ENDIF
C
      RN1 = RFACT(( J1+J2)/2 - K)
      RN2 = RFACT((-J1+J2)/2 + K)
      RN3 = RFACT(( J1-J2)/2 + K)
      RN4 = SFACT(( J1+J2)/2 + M)
      RD1 = DFLOAT(J1+1)
      RD2 = DFLOAT(J2+1)
      RD3 = RFACT(( J1+J2)/2 + K + 1)
      RD4 = SFACT(( J1+J2)/2 - M    )
      RD5 = SFACT(( J1-J2)/2 + M - 1)
      RD6 = SFACT((-J1+J2)/2 + M - 1)
      PHS = (-1.0D0)**((J2-(3*J1))/2+M)
C
      RNUM  = RN1*RN2*RN3*(RN4)**2
      RDEN  = RD1*RD2*RD3*(RD4*RD5*RD6)**2
C
      SYM3JSQ = PHS*RNUM/RDEN
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

