      SUBROUTINE ANGNR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                AA    NN    NN  GGGGGG  NN    NN RRRRRRR              C
C               AAAA   NNN   NN GG    GG NNN   NN RR    RR             C
C              AA  AA  NNNN  NN GG       NNNN  NN RR    RR             C
C             AA    AA NN NN NN GG       NN NN NN RR    RR             C
C             AAAAAAAA NN  NNNN GG   GGG NN  NNNN RRRRRRR              C
C             AA    AA NN   NNN GG    GG NN   NNN RR    RR             C
C             AA    AA NN    NN  GGGGGG  NN    NN RR    RR             C
C                                                                      C
C -------------------------------------------------------------------- C
C  ANGNR EVALUATES THE NON-RELATIVISTIC ANGULAR COEFFICIENTS OF THE    C
C  COULOMB INTERACTIONS FOR CLOSED SHELLS (K1,K2), WITH LQNA AND LQNB. C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C    BCFS(I) = B-COEFFICIENTS                                          C
C    NUS(I)  = CORRESPONDING NU-VALUES                                 C
C**********************************************************************C
      PARAMETER(MKP=9,MNU=MKP+1)
C
      COMMON/ANGL/BK(MNU,4),DK(MNU,4),HK(MNU,4),FK(MNU,4),GM(MNU,4),
     &            NUS(MNU),NUNUM,NUMIN,NUMAX
      COMMON/BLOC/NFUNA,NFUNB,LQNA,LQNB,MAXM
C
C     GENERATE A LIST OF FACTORIALS
      CALL DFACT
C
C**********************************************************************C
C     OVERWRITE THE VECTOR NUS WITH THE NUNUM VALUES OF THE TENSOR     C
C     ORDER WHICH ARE COMMON TO ALL FOUR CASES.                        C
C**********************************************************************C
C
      NUMIN = IABS(LQNA-LQNB)
      NUMAX = LQNA+LQNB
      NUNUM = 0
      DO IV=NUMIN+1,NUMAX+1
        NU    = IV-1
        LTEST = LQNA+LQNB+NU
        LEVEN = 2*(LTEST/2)
        IF(LTEST.EQ.LEVEN) THEN
          NUNUM       = NUNUM+1
          NUS(NUNUM)  = NU
          BK(NUNUM,4) = 0.5D0*ABC000(LQNA,LQNB,NU)
        ENDIF
      ENDDO
C
      RETURN
      END

