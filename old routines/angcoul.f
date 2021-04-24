      SUBROUTINE ANGCOUL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       AA    NN    NN  GGGGGG   CCCCCC   OOOOOO  UU    UU LL          C
C      AAAA   NNN   NN GG    GG CC    CC OO    OO UU    UU LL          C
C     AA  AA  NNNN  NN GG       CC       OO    OO UU    UU LL          C
C    AA    AA NN NN NN GG       CC       OO    OO UU    UU LL          C
C    AAAAAAAA NN  NNNN GG   GGG CC       OO    OO UU    UU LL          C
C    AA    AA NN   NNN GG    GG CC    CC OO    OO UU    UU LL          C
C    AA    AA NN    NN  GGGGGG   CCCCCC   OOOOOO   UUUUUU  LLLLLLLL    C
C                                                                      C
C -------------------------------------------------------------------- C
C  ANGCOUL EVALUATES THE ANGULAR COEFFICIENTS OF THE COULOMB           C
C  INTERACTIONS FOR CLOSED SHELLS IN THE (L1,L2) MANIFOLD.             C
C**********************************************************************C
      PARAMETER(MKP=9,MNU=MKP+1)
C
      COMMON/ANGL/BK(MNU,4),DK(MNU,4),HK(MNU,4),FK(MNU,4),GM(MNU,4),
     &            NUS(MNU),NUNUM,NUMIN,NUMAX
      COMMON/BLOC/NFUNA,NFUNB,LQNA,LQNB,MAXM
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
C     GENERATE LIST OF FACTORIALS
      CALL DFACT
C
C     START AND END PARAMETERS FROM TRIANGLE RULE
      NUMIN = IABS(LQNA-LQNB)
      NUMAX = LQNA+LQNB+1
      NUNUM = 0
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      DO NU=NUMIN,NUMAX
C
C       TEST WHETHER 'LQNA+LQNB+NU' ODD OR EVEN
        LTEST = LQNA+LQNB+NU
        LEVEN = 2*(LTEST/2)
C
C       ANGULAR COEFFICIENTS OF ODD PARITY ARE ZERO
        IF(LTEST.NE.LEVEN) GOTO 14
C
C       ANGULAR COEFFICIENTS OF EVEN PARITY ARE NON-ZERO
        NUNUM       = NUNUM+1
        NUS(NUNUM)  = NU
        BK(NUNUM,1) = SYM3JSQ(JRA,JRB,NU)
        BK(NUNUM,2) = SYM3JSQ(JLA,JRB,NU)
        BK(NUNUM,3) = SYM3JSQ(JRA,JLB,NU)
        BK(NUNUM,4) = SYM3JSQ(JLA,JLB,NU)
C
14      CONTINUE
      ENDDO
C
      RETURN
      END

