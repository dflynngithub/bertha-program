      SUBROUTINE SCOEFF(SCOEFF,KQN,NU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                          SWIRLES MODULE 28:                          C
C                          ------------------                          C
C         SSSSSS   CCCCCC   OOOOOO  EEEEEEEE FFFFFFFF FFFFFFFF         C
C        SS    SS CC    CC OO    OO EE       FF       FF               C
C        SS       CC       OO    OO EE       FF       FF               C
C         SSSSSS  CC       OO    OO EEEEEE   FFFFFF   FFFFFF           C
C              SS CC       OO    OO EE       FF       FF               C
C        SS    SS CC    CC OO    OO EE       FF       FF               C
C         SSSSSS   CCCCCC   OOOOOO  EEEEEEEE FF       FF               C
C                                                                      C
C -------------------------------------------------------------------- C
C  SCOEFF EVALUATES THE INTERMEDIATE COEFFICIENTS OF THE LOW-FREQUENCY C
C  BREIT INTERACTION (TABLE 2 OF GRANT AND PYPER 1976).                C
C**********************************************************************C
      DIMENSION SCOEFF(8,2),KQN(4)
C
      RNU  = DFLOAT(NU)
      DKCA = DFLOAT(KQN(3)-KQN(1))
      DKDB = DFLOAT(KQN(4)-KQN(2))
      IF(NU.GT.0) THEN
        B1 = DFLOAT(NU-1)/DFLOAT(2   *(2*NU-1))
        C1 =-DFLOAT(NU-2)/DFLOAT(2*NU*(2*NU-1))
        SCOEFF(1,1) = (RNU+DKCA)*(B1+C1*DKDB)
        SCOEFF(2,1) = (RNU+DKDB)*(B1+C1*DKCA)
        SCOEFF(3,1) = (RNU-DKCA)*(B1-C1*DKDB)
        SCOEFF(4,1) = (RNU-DKDB)*(B1-C1*DKCA)
        SCOEFF(5,1) =-(RNU+DKCA)*(B1-C1*DKDB)
        SCOEFF(6,1) =-(RNU-DKDB)*(B1+C1*DKCA)
        SCOEFF(7,1) =-(RNU-DKCA)*(B1+C1*DKDB)
        SCOEFF(8,1) =-(RNU+DKDB)*(B1-C1*DKCA)
      ELSE
        DO I=1,8
          SCOEFF(I,1) = 0.0D0
        ENDDO
      ENDIF
      IF(NU+1.GE.1) THEN
        B2 = DFLOAT(NU  )/DFLOAT( 2      *(2*NU+3))
        C2 = DFLOAT(NU+3)/DFLOAT((2*NU+2)*(2*NU+3))
        SCOEFF(1,2) = ( DKDB-RNU-1.0D0)*(B2+C2*DKCA)
        SCOEFF(2,2) = ( DKCA-RNU-1.0D0)*(B2+C2*DKDB)
        SCOEFF(3,2) = (-DKDB-RNU-1.0D0)*(B2-C2*DKCA)
        SCOEFF(4,2) = (-DKCA-RNU-1.0D0)*(B2-C2*DKDB)
        SCOEFF(5,2) =-(-DKDB-RNU-1.0D0)*(B2+C2*DKCA)
        SCOEFF(6,2) =-( DKCA-RNU-1.0D0)*(B2-C2*DKDB)
        SCOEFF(7,2) =-( DKDB-RNU-1.0D0)*(B2-C2*DKCA)
        SCOEFF(8,2) =-(-DKCA-RNU-1.0D0)*(B2+C2*DKDB)
      ELSE
        DO I=1,8
          SCOEFF(I,2) = 0.0D0
        ENDDO
      ENDIF
C
      RETURN
      END
