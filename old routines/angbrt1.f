      SUBROUTINE ANGBRT1(DKAB,DKCD,EMAT,KQN,LQN,ISEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          AA    NN    NN  GGGGGG  BBBBBBB  RRRRRRR TTTTTTTT 11        C
C         AAAA   NNN   NN GG    GG BB    BB RR    RR   TT   111        C
C        AA  AA  NNNN  NN GG       BB    BB RR    RR   TT    11        C
C       AA    AA NN NN NN GG       BBBBBBB  RR    RR   TT    11        C
C       AAAAAAAA NN  NNNN GG   GGG BB    BB RRRRRRR    TT    11        C
C       AA    AA NN   NNN GG    GG BB    BB RR    RR   TT    11        C
C       AA    AA NN    NN  GGGGGG  BBBBBBB  RR    RR   TT   1111       C
C                                                                      C
C -------------------------------------------------------------------- C
C  ANGBRT1 PERFORMS THE ANGULAR ANALYSIS FOR THE EVALUATION OF TWO     C
C  ELECTRON INTEGRALS USING RACAH ALGEBRA TECHNIQUES. (OPEN-SHELL.)    C
C -------------------------------------------------------------------- C
C  DFNOTE: NOT WORKING YET...                                          C
C**********************************************************************C
      INCLUDE 'param.h'
C
      DIMENSION KQN(4),LQN(4),JQN(4)
      DIMENSION EMAT(MNU,8),SCOEFF(8,2)
      DIMENSION DKAB(MNU,MKP+1,MKP+1),DKCD(MNU,MKP+1,MKP+1)
C
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
C
C     INITIALISE BREIT COEFFICIENTS
      DO INU=1,MNU
        DO IMU=1,8
          EMAT(INU,IMU) = 0.0D0
        ENDDO
        NUS(INU) = 0
      ENDDO
C
C**********************************************************************C
C     PARITY ANALYSIS: CHECK UNDERLYING LQN COMBINATIONS AND EXIT      C
C     IF THERE IS NO MULTIPOLE EXPANSION OF THE INTERACTION.           C
C**********************************************************************C
C
C     A AND B: PARITY OF 'LQN(1)+LQN(2)' (0 IF EVEN, 1 IF ODD)
      IPARAB = MOD(LQN(1)+LQN(2),2)
C
C     C AND D: PARITY OF 'LQN(3)+LQN(4)' (0 IF EVEN, 1 IF ODD)
      IPARCD = MOD(LQN(3)+LQN(4),2)
C
C     LQN SELECTION RULE: BOTH LQN PAIRS MUST BE OF SAME SYMMETRY
      IF(IPARAB.EQ.IPARCD) THEN
        ISEL = 1
      ELSE
        ISEL = 0
        RETURN
      ENDIF
C
C**********************************************************************C
C     MULTIPOLE EXPANSION UPPER/LOWER LIMITS SATISFY TRIANGLE RULE.    C
C**********************************************************************C
C
C     ASSIGN JQN VALUES
      DO N=1,4
        JQN(N) = 2*IABS(KQN(N))-1
      ENDDO
C
C     TENSOR LIMITS BY TRIANGLE RULE
      NUI = MAX0(IABS(JQN(1)-JQN(2))/2,IABS(JQN(3)-JQN(4))/2)
      NUF = MIN0(    (JQN(1)+JQN(2))/2,    (JQN(3)+JQN(4))/2)
C
C     JQN SELECTION RULE: TRIANGLE RULE MUST PROVIDE VALID NU RANGE
      IF(NUI.GT.NUF) THEN
        ISEL = 0
        RETURN
      ELSE
        ISEL = 1
      ENDIF
C
C**********************************************************************C
C     LQN SELECTION RULE PARITY ANALYSIS: PAIRS BOTH EVEN OR BOTH ODD. C
C     ALSO CALCULATE BREIT FACTORS IN THIS SECTION.                    C
C**********************************************************************C
C
      ISEL = 0
      LTEN = 1
      DO NU=NUI,NUF
C
C       A AND B: PARITY OF 'LQN(1)+LQN(2)+NU' (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQN(1)+LQN(2)+NU,2)
C
C       C AND D: PARITY OF 'LQN(3)+LQN(4)+NU' (0 IF EVEN, 1 IF ODD)
        IPARCD = MOD(LQN(3)+LQN(4)+NU,2)
C
C       CASE 1: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH ODD (UNLESS NU=0)
        IF(IPARAB.EQ.1.AND.IPARCD.EQ.1.AND.NU.NE.0) THEN
C
C         THERE IS A NON-TRIVIAL TENSOR EXPANSION
          ISEL = 1
C
C         BREIT COEFFICIENT IS SIMPLE ENOUGH TO CALCULATE HERE
          RNU  = DFLOAT(NU*(NU+1))
          COEF =-DFLOAT((KQN(1)+KQN(2))*(KQN(3)+KQN(4)))/RNU
C
C         SAVE THIS PARTICULAR NU VALUE TO A LIST
          NUS(LTEN) = NU
C
C         COPY ACROSS TO ARRAY
          DO IMU=1,8
            EMAT(LTEN,IMU) = EMAT(LTEN,IMU) + COEF
          ENDDO
C
        ENDIF
C
C       CASE 2: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH EVEN
        IF(IPARAB.EQ.0.AND.IPARCD.EQ.0) THEN
C
C         THERE IS A NON-TRIVIAL TENSOR EXPANSION
          ISEL = 1
C
C         CALCULATE BREIT COEFFICIENTS WITH CALL TO BRCOEF1
          CALL BRCOEF1(SCOEFF,KQN,NU)
C
C         SAVE THIS PARTICULAR NU VALUE TO A LIST
          NUS(LTEN) = NU-1
C
C         COPY ACROSS TO ARRAY
          DO IMU=1,8
            EMAT(LTEN,IMU) = EMAT(LTEN,IMU) + SCOEFF(IMU,1)
          ENDDO
C
C         INCREASE TENSOR ORDER ONLY IF NU NONZERO
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SAVE THIS PARTICULAR NU VALUE TO A LIST
          NUS(LTEN) = NU+1
C          
C         COPY ACROSS TO ARRAY
          DO IMU=1,8
            EMAT(LTEN,IMU) = EMAT(LTEN,IMU) + SCOEFF(IMU,2)
          ENDDO
C
        ENDIF
C
      ENDDO
C
C     RESET THE LOWER AND UPPER LIMITS OF THE EXPANSION
      NUI = NUS(1)
      NUF = NUS(LTEN)
      NUNUM = LTEN
C
C     IF NO NU VALUES WERE ALLOWED DURING THIS RUN, EXIT PROCEDURE
      IF(ISEL.EQ.0) THEN
        RETURN
      ENDIF
C
C**********************************************************************C
C     ANGULAR FACTORS: EVALUATE AN ANGULAR FACTOR FOR EVERY |MQN|      C
C     COMBINATION IN THE MULTIPOLE EXPANSION OVER ALLOWED NU VALUES.   C
C**********************************************************************C
C
C     LOOP OVER MQN(A) AND MQN(B) VALUES
      DO MA=1,IABS(KQN(1))
        MJA = 2*MA-1
        DO MB=1,IABS(KQN(2))
          MJB = 2*MB-1
C
C         LOOP OVER ALL SURVIVING TENSOR ORDERS
          DO LTEN=1,NUNUM
C
C           READ THE ACTUAL ORDER NU
            NU = NUS(LTEN)
C
C           GENERATE AN ANGULAR FACTOR FOR ALL |MQN| SIGNS
            DKAB(LTEN,MJA  ,MJB  ) = DK(JQN(1),-MJA,JQN(2),-MJB,NU)
            DKAB(LTEN,MJA+1,MJB  ) = DK(JQN(1), MJA,JQN(2),-MJB,NU)
            DKAB(LTEN,MJA  ,MJB+1) = DK(JQN(1),-MJA,JQN(2), MJB,NU)
            DKAB(LTEN,MJA+1,MJB+1) = DK(JQN(1), MJA,JQN(2), MJB,NU)
C
          ENDDO
C
        ENDDO
      ENDDO
C
C     THE ARGUMENTS OF DK COEFFICIENTS ARE REVERSED IN THE CASE OF THE
C     CD PAIRS IN ORDER TO ACCOMMODATE THE RELATION:
C               DK(J,M,J',M',L) = ((-1)^Q)*DK(J',M',J,M,L).
C
C     LOOP OVER MQN(C) AND MQN(D) VALUES
      DO MC=1,IABS(KQN(3))
        MJC = 2*MC-1
        DO MD=1,IABS(KQN(4))
          MJD = 2*MD-1
C
C         LOOP OVER ALL SURVIVING TENSOR ORDERS
          DO LTEN=1,NUNUM
C
C           READ THE ACTUAL ORDER NU
            NU = NUS(LTEN)
C
C           GENERATE AN ANGULAR FACTOR FOR ALL |MQN| SIGNS
            DKCD(LTEN,MJC  ,MJD  ) = DK(JQN(4),-MJD,JQN(3),-MJC,NU)
            DKCD(LTEN,MJC+1,MJD  ) = DK(JQN(4),-MJD,JQN(3), MJC,NU)
            DKCD(LTEN,MJC  ,MJD+1) = DK(JQN(4), MJD,JQN(3),-MJC,NU)
            DKCD(LTEN,MJC+1,MJD+1) = DK(JQN(4), MJD,JQN(3), MJC,NU)
C
          ENDDO
C
        ENDDO
      ENDDO
C
      RETURN
      END

