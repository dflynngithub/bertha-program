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
C  ▶ NUNUM - NUMBER OF NU VALUES THAT SATISFY PARITY RESTRICTION RULE. C
C  ▶ NUI - MINIMUM NU VALUE IN THIS MANIFOLD.                          C
C  ▶ NUF - MAXIMUM NU VALUE IN THIS MANIFOLD.                          C
C  ▶ ELL(MNU,4) - ANGULAR TERMS FOR CLOSED BREIT EK(JA,JB;LL) TERMS    C
C  ▶ ESS(MNU,4) - ANGULAR TERMS FOR CLOSED BREIT EK(JA,JB;SS) TERMS    C
C  ▶ ESL(MNU,4) - ANGULAR TERMS FOR CLOSED BREIT EK(JA,JB;SL) TERMS    C
C  ▶ GSL(MNU,4) - ANGULAR TERMS FOR CLOSED BREIT EK(JA,JB;SL) TERMS    C
C**********************************************************************C
      INCLUDE 'param.h'
C
      DIMENSION SCOEFF(4,2)
C
      COMMON/B0QN/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
      COMMON/XTNS/BK(MNU,4),ELL(MNU,4),ESS(MNU,4),ESL(MNU,4),GSL(MNU,4)
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
C     SPECIFY ALLOWED NU VALUES (PARITY LQNA+LQNB+NU MUST BE ODD)
      NUNUM = 0
      DO NU=0,IABS(LQNA+LQNB+1)
        IF(NU.GE.0.AND.MOD(LQNA+LQNB+NU,2).EQ.1) THEN
          NUNUM = NUNUM+1
          NUS(NUNUM) = NU
        ENDIF
      ENDDO
C
C     BREIT TENSOR ORDER LIMITS BASED ON PARITY CHECK
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
C       ANGULAR FACTOR FOR THIS (KQNA,KQNB,NU)
        AJJN = ANGSQJJ(JQNA,JQNB,NU)
C
C       PARITY OF 'LQNA+LQNB+NU' (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQNA+LQNB+NU,2)
C
C       CASE 1: ODD-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.1.AND.NU.NE.0) THEN
C
C         INTERMEDIATE (KQN,NU) COEFFICIENTS
          RCOEFF = DFLOAT((KQNA+KQNB)*(KQNA+KQNB))/DFLOAT(NU*(NU+1))
C
C         THIS CONTRIBUTION BELONGS IN "NU" BIN
          NUS(LTEN) = NU
C
C         CONTRIBUTIONS TO TENSOR ORDER "NU"
          ELL(LTEN,1) = ELL(LTEN,1) + AJJN*RCOEFF
          ESS(LTEN,1) = ESS(LTEN,1) + AJJN*RCOEFF
          ESL(LTEN,1) = ESL(LTEN,1) + AJJN*RCOEFF
C
        ENDIF
C
C       CASE 1: EVEN-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.0) THEN
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          CALL BRCOEF0(SCOEFF,KQNA,KQNB,NU)
C
C         FIRST CONTRIBUTION BELONGS IN "NU-1" BIN
          NUS(LTEN) = NU-1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,1) = ELL(LTEN,1) + AJJN*SCOEFF(1,1)
          ESL(LTEN,1) = ESL(LTEN,1) + AJJN*SCOEFF(2,1)
          ESS(LTEN,1) = ESS(LTEN,1) + AJJN*SCOEFF(3,1)
          GSL(LTEN,1) = GSL(LTEN,1) + AJJN*SCOEFF(4,1)
C
C         INCREASE TENSOR LIST LENGTH IF ENTRIES ARE NON-ZERO
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SECOND CONTRIBUTION BELONGS IN "NU+1" BIN
          NUS(LTEN) = NU+1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,1) = ELL(LTEN,1) + AJJN*SCOEFF(1,2)
          ESL(LTEN,1) = ESL(LTEN,1) + AJJN*SCOEFF(2,2)
          ESS(LTEN,1) = ESS(LTEN,1) + AJJN*SCOEFF(3,2)
          GSL(LTEN,1) = GSL(LTEN,1) + AJJN*SCOEFF(4,2)
C
        ENDIF
C
101     CONTINUE
C
      ENDDO
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
C       ANGULAR FACTOR FOR THIS (KQNA,KQNB,NU)
        AJJN = ANGSQJJ(JQNA,JQNB,NU)
C
C       PARITY OF 'LQNA+LQNB+NU' (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQNA+LQNB+NU,2)
C
C       CASE 1: ODD-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.1.AND.NU.NE.0) THEN
C
C         INTERMEDIATE (KQN,NU) COEFFICIENTS
          RCOEFF = DFLOAT((KQNA+KQNB)*(KQNA+KQNB))/DFLOAT(NU*(NU+1))
C
C         THIS CONTRIBUTION BELONGS IN "NU" BIN
          NUS(LTEN) = NU
C
C         CONTRIBUTIONS TO TENSOR ORDER "NU"
          ELL(LTEN,2) = ELL(LTEN,2) + AJJN*RCOEFF
          ESS(LTEN,2) = ESS(LTEN,2) + AJJN*RCOEFF
          ESL(LTEN,2) = ESL(LTEN,2) + AJJN*RCOEFF
C
        ENDIF
C
C       CASE 1: EVEN-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.0) THEN
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          CALL BRCOEF0(SCOEFF,KQNA,KQNB,NU)
C
C         FIRST CONTRIBUTION BELONGS IN "NU-1" BIN
          NUS(LTEN) = NU-1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,2) = ELL(LTEN,2) + AJJN*SCOEFF(1,1)
          ESL(LTEN,2) = ESL(LTEN,2) + AJJN*SCOEFF(2,1)
          ESS(LTEN,2) = ESS(LTEN,2) + AJJN*SCOEFF(3,1)
          GSL(LTEN,2) = GSL(LTEN,2) + AJJN*SCOEFF(4,1)
C
C         INCREASE TENSOR LIST LENGTH IF ENTRIES ARE NON-ZERO
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SECOND CONTRIBUTION BELONGS IN "NU+1" BIN
          NUS(LTEN) = NU+1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,2) = ELL(LTEN,2) + AJJN*SCOEFF(1,2)
          ESL(LTEN,2) = ESL(LTEN,2) + AJJN*SCOEFF(2,2)
          ESS(LTEN,2) = ESS(LTEN,2) + AJJN*SCOEFF(3,2)
          GSL(LTEN,2) = GSL(LTEN,2) + AJJN*SCOEFF(4,2)
C
        ENDIF
C
201     CONTINUE
C
      ENDDO
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
C       ANGULAR FACTOR FOR THIS (KQNA,KQNB,NU)
        AJJN = ANGSQJJ(JQNA,JQNB,NU)
C
C       PARITY OF 'LQNA+LQNB+NU' (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQNA+LQNB+NU,2)
C
C       CASE 1: ODD-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.1.AND.NU.NE.0) THEN
C
C         INTERMEDIATE (KQN,NU) COEFFICIENTS
          RCOEFF = DFLOAT((KQNA+KQNB)*(KQNA+KQNB))/DFLOAT(NU*(NU+1))
C
C         THIS CONTRIBUTION BELONGS IN "NU" BIN
          NUS(LTEN) = NU
C
C         CONTRIBUTIONS TO TENSOR ORDER "NU"
          ELL(LTEN,3) = ELL(LTEN,3) + AJJN*RCOEFF
          ESS(LTEN,3) = ESS(LTEN,3) + AJJN*RCOEFF
          ESL(LTEN,3) = ESL(LTEN,3) + AJJN*RCOEFF
C
        ENDIF
C
C       CASE 1: EVEN-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.0) THEN
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          CALL BRCOEF0(SCOEFF,KQNA,KQNB,NU)
C
C         FIRST CONTRIBUTION BELONGS IN "NU-1" BIN
          NUS(LTEN) = NU-1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,3) = ELL(LTEN,3) + AJJN*SCOEFF(1,1)
          ESL(LTEN,3) = ESL(LTEN,3) + AJJN*SCOEFF(2,1)
          ESS(LTEN,3) = ESS(LTEN,3) + AJJN*SCOEFF(3,1)
          GSL(LTEN,3) = GSL(LTEN,3) + AJJN*SCOEFF(4,1)
C
C         INCREASE TENSOR LIST LENGTH IF ENTRIES ARE NON-ZERO
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SECOND CONTRIBUTION BELONGS IN "NU+1" BIN
          NUS(LTEN) = NU+1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,3) = ELL(LTEN,3) + AJJN*SCOEFF(1,2)
          ESL(LTEN,3) = ESL(LTEN,3) + AJJN*SCOEFF(2,2)
          ESS(LTEN,3) = ESS(LTEN,3) + AJJN*SCOEFF(3,2)
          GSL(LTEN,3) = GSL(LTEN,3) + AJJN*SCOEFF(4,2)
C
        ENDIF
C
301     CONTINUE
C
      ENDDO
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
C       ANGULAR FACTOR FOR THIS (KQNA,KQNB,NU)
        AJJN = ANGSQJJ(JQNA,JQNB,NU)
C
C       PARITY OF 'LQNA+LQNB+NU' (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQNA+LQNB+NU,2)
C
C       CASE 1: ODD-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.1.AND.NU.NE.0) THEN
C
C         INTERMEDIATE (KQN,NU) COEFFICIENTS
          RCOEFF = DFLOAT((KQNA+KQNB)*(KQNA+KQNB))/DFLOAT(NU*(NU+1))
C
C         THIS CONTRIBUTION BELONGS IN "NU" BIN
          NUS(LTEN) = NU
C
C         CONTRIBUTIONS TO TENSOR ORDER "NU"
          ELL(LTEN,4) = ELL(LTEN,4) + AJJN*RCOEFF
          ESS(LTEN,4) = ESS(LTEN,4) + AJJN*RCOEFF
          ESL(LTEN,4) = ESL(LTEN,4) + AJJN*RCOEFF
C
        ENDIF
C
C       CASE 1: EVEN-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.0) THEN
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          CALL BRCOEF0(SCOEFF,KQNA,KQNB,NU)
C
C         FIRST CONTRIBUTION BELONGS IN "NU-1" BIN
          NUS(LTEN) = NU-1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,4) = ELL(LTEN,4) + AJJN*SCOEFF(1,1)
          ESL(LTEN,4) = ESL(LTEN,4) + AJJN*SCOEFF(2,1)
          ESS(LTEN,4) = ESS(LTEN,4) + AJJN*SCOEFF(3,1)
          GSL(LTEN,4) = GSL(LTEN,4) + AJJN*SCOEFF(4,1)
C
C         INCREASE TENSOR LIST LENGTH IF ENTRIES ARE NON-ZERO
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SECOND CONTRIBUTION BELONGS IN "NU+1" BIN
          NUS(LTEN) = NU+1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,4) = ELL(LTEN,4) + AJJN*SCOEFF(1,2)
          ESL(LTEN,4) = ESL(LTEN,4) + AJJN*SCOEFF(2,2)
          ESS(LTEN,4) = ESS(LTEN,4) + AJJN*SCOEFF(3,2)
          GSL(LTEN,4) = GSL(LTEN,4) + AJJN*SCOEFF(4,2)
C
        ENDIF
C
401     CONTINUE
C
      ENDDO
C
400   CONTINUE
C
      RETURN
      END