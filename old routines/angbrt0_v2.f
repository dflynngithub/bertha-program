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
      PARAMETER(MKP=13,MNU=MKP+1)
C
      DIMENSION SCOEF(4,2)
C
      COMMON/ANGL/BK(MNU,4),ELL(MNU,4),ESS(MNU,4),ESL(MNU,4),GSL(MNU,4)
      COMMON/BSQN/NBASA,NBASB,LQNA,LQNB,MAXM
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
      DO LTEN=1,NUNUM
C
C       VALUE OF NU
        NU = NUS(LTEN)
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
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,1) = ELL(LTEN,1) - SCOEF(1,1)*RAW
          ESL(LTEN,1) = ESL(LTEN,1) - SCOEF(2,1)*RAW
          ESS(LTEN,1) = ESS(LTEN,1) - SCOEF(3,1)*RAW
          GSL(LTEN,1) = GSL(LTEN,1) - SCOEF(4,1)*RAW
C
C         STEP PARAMETER FOR ADDITIONS TO COEFFICIENTS
          IF(NU.GT.0) THEN
            MTEN = LTEN+1
          ENDIF
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(MTEN,1) = ELL(MTEN,1) - SCOEF(1,2)*RAW
          ESL(MTEN,1) = ESL(MTEN,1) - SCOEF(2,2)*RAW
          ESS(MTEN,1) = ESS(MTEN,1) - SCOEF(3,2)*RAW
          GSL(MTEN,1) = GSL(MTEN,1) - SCOEF(4,2)*RAW
C
        ENDIF
C
101     CONTINUE
C
      ENDDO
C
      WRITE(*,*) 'KQNA < 0, KQNB < 0:'
      WRITE(*,*) 'NUI = ',NUI,'NUF = ',NUF,'NUNUM = ',NUNUM
      DO L=1,NUNUM
        WRITE(*,*) NUS(L),ELL(L,1),ESS(L,1),ESL(L,1),GSL(L,1)
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
      DO LTEN=1,NUNUM
C
C       VALUE OF NU
        NU = NUS(LTEN)
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
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,2) = ELL(LTEN,2) - SCOEF(1,1)*RAW
          ESL(LTEN,2) = ESL(LTEN,2) - SCOEF(2,1)*RAW
          ESS(LTEN,2) = ESS(LTEN,2) - SCOEF(3,1)*RAW
          GSL(LTEN,2) = GSL(LTEN,2) - SCOEF(4,1)*RAW
C
C         STEP PARAMETER FOR ADDITIONS TO COEFFICIENTS
          IF(NU.GT.0) THEN
            MTEN = LTEN+1
          ENDIF
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(MTEN,2) = ELL(MTEN,2) - SCOEF(1,2)*RAW
          ESL(MTEN,2) = ESL(MTEN,2) - SCOEF(2,2)*RAW
          ESS(MTEN,2) = ESS(MTEN,2) - SCOEF(3,2)*RAW
          GSL(MTEN,2) = GSL(MTEN,2) - SCOEF(4,2)*RAW
C
        ENDIF
C
201     CONTINUE
C
      ENDDO
C
      WRITE(*,*) 'KQNA < 0, KQNB > 0:'
      WRITE(*,*) 'NUI = ',NUI,'NUF = ',NUF,'NUNUM = ',NUNUM
      DO L=1,NUNUM
        WRITE(*,*) NUS(L),ELL(L,2),ESS(L,2),ESL(L,2),GSL(L,2)
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
      DO LTEN=1,NUNUM
C
C       VALUE OF NU
        NU = NUS(LTEN)
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
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,3) = ELL(LTEN,3) - SCOEF(1,1)*RAW
          ESL(LTEN,3) = ESL(LTEN,3) - SCOEF(2,1)*RAW
          ESS(LTEN,3) = ESS(LTEN,3) - SCOEF(3,1)*RAW
          GSL(LTEN,3) = GSL(LTEN,3) - SCOEF(4,1)*RAW
C
C         STEP PARAMETER FOR ADDITIONS TO COEFFICIENTS
          IF(NU.GT.0) THEN
            MTEN = LTEN+1
          ENDIF
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(MTEN,3) = ELL(MTEN,3) - SCOEF(1,2)*RAW
          ESL(MTEN,3) = ESL(MTEN,3) - SCOEF(2,2)*RAW
          ESS(MTEN,3) = ESS(MTEN,3) - SCOEF(3,2)*RAW
          GSL(MTEN,3) = GSL(MTEN,3) - SCOEF(4,2)*RAW
C
        ENDIF
C
301     CONTINUE
C
      ENDDO
C
      WRITE(*,*) 'KQNA > 0, KQNB < 0:'
      WRITE(*,*) 'NUI = ',NUI,'NUF = ',NUF,'NUNUM = ',NUNUM
      DO L=1,NUNUM
        WRITE(*,*) NUS(L),ELL(L,3),ESS(L,3),ESL(L,3),GSL(L,3)
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
      DO LTEN=1,NUNUM
C
C       VALUE OF NU
        NU = NUS(LTEN)
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
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,4) = ELL(LTEN,4) - SCOEF(1,1)*RAW
          ESL(LTEN,4) = ESL(LTEN,4) - SCOEF(2,1)*RAW
          ESS(LTEN,4) = ESS(LTEN,4) - SCOEF(3,1)*RAW
          GSL(LTEN,4) = GSL(LTEN,4) - SCOEF(4,1)*RAW
C
C         STEP PARAMETER FOR ADDITIONS TO COEFFICIENTS
          IF(NU.GT.0) THEN
            MTEN = LTEN+1
          ENDIF
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(MTEN,4) = ELL(MTEN,4) - SCOEF(1,2)*RAW
          ESL(MTEN,4) = ESL(MTEN,4) - SCOEF(2,2)*RAW
          ESS(MTEN,4) = ESS(MTEN,4) - SCOEF(3,2)*RAW
          GSL(MTEN,4) = GSL(MTEN,4) - SCOEF(4,2)*RAW
C
        ENDIF
C
401     CONTINUE
C
      ENDDO
C
      WRITE(*,*) 'KQNA > 0, KQNB > 0:'
      WRITE(*,*) 'NUI = ',NUI,'NUF = ',NUF,'NUNUM = ',NUNUM
      DO L=1,NUNUM
        WRITE(*,*) NUS(L),ELL(L,4),ESS(L,4),ESL(L,4),GSL(L,4)
      ENDDO
C
400   CONTINUE
C
      RETURN
      END

