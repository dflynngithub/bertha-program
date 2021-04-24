      SUBROUTINE BRCOEF0(SCOEFF,KQNA,KQNB,NU)
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
      RK  = DFLOAT(KQNB-KQNA)
C
      IF(NU.GT.0) THEN
        RM = RU-1.0D0
        B1 = (RM+2.0D0)/(2.0D0*   (2.0D0*RM+1.0D0))
        C1 =-(RM-1.0D0)/(2.0D0*RU*(2.0D0*RM+1.0D0))
        SCOEFF(1,1) = (RK+RU)*(C1*RK+B1)
        SCOEFF(2,1) = -B1*RU + C1*RK*RK
        SCOEFF(3,1) = (RK-RU)*(C1*RK-B1)
        SCOEFF(4,1) =  RK    *(C1*RU-B1)
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
        SCOEFF(1,2) = (RK-RP)*(C2*RK+B2)
        SCOEFF(2,2) =  B2*RP + C2*RK*RK 
        SCOEFF(3,2) = (RK+RP)*(C2*RK-B2)
        SCOEFF(4,2) =  RK    *(C2*RP+B2)
      ELSE
        SCOEFF(1,2) = 0.0D0
        SCOEFF(2,2) = 0.0D0
        SCOEFF(3,2) = 0.0D0
        SCOEFF(4,2) = 0.0D0
      ENDIF
C
      RETURN
      END

