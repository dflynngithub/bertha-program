      SUBROUTINE BRCOEF1(SCOEFF,KQN,NU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      BBBBBBB  RRRRRRR   CCCCCC   OOOOOO  EEEEEEEE FFFFFFFF 11        C
C      BB    BB RR    RR CC    CC OO    OO EE       FF      111        C
C      BB    BB RR    RR CC       OO    OO EE       FF       11        C
C      BBBBBBB  RR    RR CC       OO    OO EEEEEE   FFFFFF   11        C
C      BB    BB RRRRRRR  CC       OO    OO EE       FF       11        C
C      BB    BB RR    RR CC    CC OO    OO EE       FF       11        C
C      BBBBBBB  RR    RR  CCCCCC   OOOOOO  EEEEEEEE FF      1111       C
C                                                                      C
C -------------------------------------------------------------------- C
C  BRCOEF1 EVALUATES THE INTERMEDIATE COEFFICIENTS OF THE BREIT        C
C  INTERACTION IN THE GENERAL CASE (TABLE 2 OF GRANT AND PYPER 1976).  C
C**********************************************************************C
      DIMENSION SCOEFF(8,2),KQN(4)
C
      RNU = DFLOAT(NU)
      RK1 = DFLOAT(KQN(2)-KQN(1))
      RK2 = DFLOAT(KQN(4)-KQN(3))
C
      IF(NU-1.GE.0) THEN
        B1 = DFLOAT(NU+1)/DFLOAT(2   *(2*NU+1))
        C1 =-DFLOAT(NU-2)/DFLOAT(2*NU*(2*NU+1))
        SCOEFF(1,1) = (RNU+RK1)*(B1+C1*RK2)
        SCOEFF(2,1) = (RNU+RK2)*(B1+C1*RK1)
        SCOEFF(3,1) = (RNU-RK1)*(B1-C1*RK2)
        SCOEFF(4,1) = (RNU-RK2)*(B1-C1*RK1)
        SCOEFF(5,1) =-(RNU+RK1)*(B1-C1*RK2)
        SCOEFF(6,1) =-(RNU-RK2)*(B1+C1*RK1)
        SCOEFF(7,1) =-(RNU-RK1)*(B1+C1*RK2)
        SCOEFF(8,1) =-(RNU+RK2)*(B1-C1*RK1)
      ELSE
        SCOEFF(1,1) = 0.0D0
        SCOEFF(2,1) = 0.0D0
        SCOEFF(3,1) = 0.0D0
        SCOEFF(4,1) = 0.0D0
        SCOEFF(5,1) = 0.0D0
        SCOEFF(6,1) = 0.0D0
        SCOEFF(7,1) = 0.0D0
        SCOEFF(8,1) = 0.0D0
      ENDIF
C
      IF(NU+1.GE.1) THEN
        B2 = DFLOAT(NU  )/DFLOAT(2       *(2*NU+3))
        C2 = DFLOAT(NU+3)/DFLOAT(2*(NU+1)*(2*NU+3))
        SCOEFF(1,2) = ( RK2-RNU-1.0D0)*(B2+C2*RK1)
        SCOEFF(2,2) = ( RK1-RNU-1.0D0)*(B2+C2*RK2)
        SCOEFF(3,2) = (-RK2-RNU-1.0D0)*(B2-C2*RK1)
        SCOEFF(4,2) = (-RK1-RNU-1.0D0)*(B2-C2*RK2)
        SCOEFF(5,2) =-(-RK2-RNU-1.0D0)*(B2+C2*RK1)
        SCOEFF(6,2) =-( RK1-RNU-1.0D0)*(B2-C2*RK2)
        SCOEFF(7,2) =-( RK2-RNU-1.0D0)*(B2-C2*RK1)
        SCOEFF(8,2) =-(-RK1-RNU-1.0D0)*(B2+C2*RK2)
      ELSE
        SCOEFF(1,2) = 0.0D0
        SCOEFF(2,2) = 0.0D0
        SCOEFF(3,2) = 0.0D0
        SCOEFF(4,2) = 0.0D0
        SCOEFF(5,2) = 0.0D0
        SCOEFF(6,2) = 0.0D0
        SCOEFF(7,2) = 0.0D0
        SCOEFF(8,2) = 0.0D0
      ENDIF
C
      RETURN
      END

