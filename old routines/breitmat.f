      SUBROUTINE BREITMAT(K1,K2,ELL,ESS,ESL,GSL,NUS,NUNUM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                          SWIRLES MODULE 23:                          C
C                          ------------------                          C
C BBBBBBB  RRRRRRR  EEEEEEEE IIII TTTTTTTT MM       MM    AA  TTTTTTTT C
C BB    BB RR    RR EE        II     TT    MMM     MMM   AAAA    TT    C
C BB    BB RR    RR EE        II     TT    MMMM   MMMM  AA  AA   TT    C
C BBBBBBB  RR    RR EEEEEE    II     TT    MM MM MM MM AA    AA  TT    C
C BB    BB RRRRRRR  EE        II     TT    MM  MMM  MM AAAAAAAA  TT    C
C BB    BB RR    RR EE        II     TT    MM   M   MM AA    AA  TT    C
C BBBBBBB  RR    RR EEEEEEEE IIII    TT    MM       MM AA    AA  TT    C
C                                                                      C
C -------------------------------------------------------------------- C
C  BREITMAT EVALUATES THE MATRIX COEFFICIENTS OF THE LOW-FREQUENCY     C
C  INTERACTION FOR RELATIVISTIC CLOSED SHELLS.                         C
C**********************************************************************C
      DIMENSION S(4,2),NUS(1),ELL(1),ESS(1),ESL(1),GSL(1)
C
C     CALCULATE LQN1 AND 2*JQN1 FROM KQN1
      IF(K1.LT.0) THEN
        L1 = IABS(K1)-1
      ELSE
        L1 = K1
      ENDIF
      J1 = 2*IABS(K1)-1
C
C     CALCULATE LQN2 AND 2*JQN2 FROM KQN2
      IF(K2.LT.0) THEN
        L2 = IABS(K2)-1
      ELSE
        L2 = K2
      ENDIF
      J2 = 2*IABS(K2)-1
C
C     INITIALISE COEFFICIENT ARRAYS
      DO I=1,10
        ELL(I) = 0.0D0
        ESL(I) = 0.0D0
        ESS(I) = 0.0D0
        GSL(I) = 0.0D0
      ENDDO
C
C     START AND END PARAMETERS FROM TRIANGLE RULE
      NUI = IABS(J1-J2)/2
      NUF = (J1+J2)/2
      NUNUM = 1
C
C**********************************************************************C      
C     NOTE THAT THESE COEFFICIENTS DIFFER BY A FACTOR                  C
C     OF (-1) COMPARED WITH THE FORMULAE OF GRANT AND PYPER FOR        C
C     CONSISTENCY WITH THE USUAL FORM OF THE DHFB EQUATIONS            C
C**********************************************************************C
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      DO NU=NUI,NUF
C
C       GENERATE SQUARE OF 3J-SYMBOL
        CALL SYM3JSQ(J1,J2,NU,VALUE)
C
C       DETERMINE PARITY OF COMBINATION L1,L2,NU
        LTEST = L1+L2+NU
        LEVEN = 2*(LTEST/2)
C
C       STEP PARAMETER FOR ADDITIONS TO COEFFICIENTS
        NSTEP = 1
        IF(NU.LE.0) THEN
          NSTEP = 0
        ENDIF
C        
        IF(LTEST.NE.LEVEN.AND.NU.NE.0) THEN
C       CASE 1: L1+L2+NU IS ODD
          SYMFAC = DFLOAT((K1+K2)*(K1+K2))/DFLOAT(NU*(NU+1))
          SYMFAC = SYMFAC*VALUE
          ELL(NUNUM) = ELL(NUNUM) + SYMFAC
          ESS(NUNUM) = ESS(NUNUM) + SYMFAC
          ESL(NUNUM) = ESL(NUNUM) + SYMFAC
          NUS(NUNUM) = NU
        ELSEIF(LTEST.EQ.LEVEN) THEN
C       CASE 2: L1+L2+NU IS EVEN
          CALL EXCHNG(S,K1,K2,NU)
          NUS(NUNUM) = NU-1
          ELL(NUNUM) = ELL(NUNUM) - S(1,1)*VALUE
          ESL(NUNUM) = ESL(NUNUM) - S(2,1)*VALUE
          ESS(NUNUM) = ESS(NUNUM) - S(3,1)*VALUE
          GSL(NUNUM) = GSL(NUNUM) - S(4,1)*VALUE
C
          NUNUM = NUNUM + NSTEP
          NUS(NUNUM) = NU+1
          ELL(NUNUM) = ELL(NUNUM) - S(1,2)*VALUE
          ESL(NUNUM) = ESL(NUNUM) - S(2,2)*VALUE
          ESS(NUNUM) = ESS(NUNUM) - S(3,2)*VALUE
          GSL(NUNUM) = GSL(NUNUM) - S(4,2)*VALUE
        ENDIF
      ENDDO
C
C     WRITE BREIT COEFFICIENTS TO TERMINAL
      WRITE(6,69) K1,K2
69    FORMAT(//2X,'MATRIX BREIT   COEFFICIENTS:',2X,'KAPPA(A) = ',I4,
     &                                12X,'KAPPA(B) = ',I4/2X,62('=')/)
C
      WRITE(6,67)
67    FORMAT(/2X,'NU',6X,'ELL',12X,'ESS',12X,'ESL',12X,'GSL'/
     &                                                     12X,62('-'))
C
      DO I=1,NUNUM
        WRITE(6,68) NUS(I),ELL(I),ESS(I),ESL(I),GSL(I)
      ENDDO
68    FORMAT(2X,I2,4(2X,G13.6))
C
      RETURN
      END

