THESE ROUTINES WERE USED TO EXTERNALLY CALCULATE THE MAGNITUDE OF A LIST
OF E-COEFFS AND R-INTS WITHIN ERI AND BII. DONE INSIDE THE ROUTINES NOW.

C   [I] RDASUM: SUM OF MAGNITUDES OF A REAL VECTOR.                    C
C   [J] CDASUM: SUM OF MAGNITUDES OF A COMPLEX VECTOR.                 C

      FUNCTION RDASUM(N,DX,INCX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       RRRRRRR  DDDDDDD     AA     SSSSSS  UU    UU MM       MM       C
C       RR    RR DD    DD   AAAA   SS    SS UU    UU MMM     MMM       C
C       RR    RR DD    DD  AA  AA  SS       UU    UU MMMM   MMMM       C
C       RR    RR DD    DD AA    AA  SSSSSS  UU    UU MM MM MM MM       C
C       RRRRRRR  DD    DD AAAAAAAA       SS UU    UU MM  MMM  MM       C
C       RR    RR DD    DD AA    AA SS    SS UU    UU MM   M   MM       C
C       RR    RR DDDDDDD  AA    AA  SSSSSS   UUUUUU  MM       MM       C
C                                                                      C
C -------------------------------------------------------------------- C
C  RDASUM RETURNS THE SUM OF MAGNITUDES OF A VECTOR OF DOUBLES DX(N).  C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    N     - NUMBER OF ELEMENTS IN INPUT VECTOR(S).                    C
C    DX    - DOUBLE PRECISION VECTOR WITH N ELEMENTS.                  C
C    INCX  - STORAGE SPACING BETWEEN ELEMENTS OF DX.                   C
C  OUTPUT:                                                             C
C    RDASUM - DOUBLE PRECISION RESULT (ZERO IF N < 0).                 C
C**********************************************************************C
C
      DIMENSION DX(N)
C
C     INITIALISE COUNTER
      RDASUM = 0.D0
C
C     IF VECTOR LENGTH IS ZERO, RESULT IS ZERO
      IF(N.LE.0) THEN
        RETURN
      ENDIF
C
C     DECISION TREE FOR INCREMENT STEP SIZE
      IF(INCX.EQ.1) THEN
        GOTO 20
      ENDIF

      NS = N*INCX
      DO I=1,NS,INCX
        RDASUM = RDASUM + DABS(DX(I))
      ENDDO
      RETURN
C
20    CONTINUE
C
C     CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
      M = MOD(N,6)
C
      IF(M.EQ.0) GOTO 40
C
      DO I=1,M
        RDASUM = RDASUM + DABS(DX(I))
      ENDDO
C      
      IF(N.LT.6) RETURN
C
   40 CONTINUE
C   
      MP1 = M+1
      DO I = MP1,N,6
        RDASUM = RDASUM + DABS(DX(I  )) + DABS(DX(I+1)) + DABS(DX(I+2))
     &                  + DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION CDASUM(N,DX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        CCCCCC  DDDDDDD     AA     SSSSSS  UU    UU MM       MM       C
C       CC    CC DD    DD   AAAA   SS    SS UU    UU MMM     MMM       C
C       CC       DD    DD  AA  AA  SS       UU    UU MMMM   MMMM       C
C       CC       DD    DD AA    AA  SSSSSS  UU    UU MM MM MM MM       C
C       CC       DD    DD AAAAAAAA       SS UU    UU MM  MMM  MM       C
C       CC    CC DD    DD AA    AA SS    SS UU    UU MM   M   MM       C
C        CCCCCC  DDDDDDD  AA    AA  SSSSSS   UUUUUU  MM       MM       C
C                                                                      C
C -------------------------------------------------------------------- C
C  CDASUM RETURNS THE SUM OF MAGNITUDES OF A COMPLEX VECTOR DX(N).     C
C**********************************************************************C
C
      COMPLEX*16 DX(N)
C
C     INITIALISE COUNTER
      CDASUM = 0.0D0
C
C     IF VECTOR LENGTH IS ZERO, RESULT IS ZERO
      IF(N.LE.0) RETURN
C
C     CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
      M = MOD(N,6)
C
      IF(M.EQ.0) GOTO 40
C
      DO I=1,M
        CDASUM = CDASUM + ABS(DX(I))
      ENDDO
C
      IF(N.LT.6) RETURN
C
40    CONTINUE
C
      MP1 = M+1
      DO I=MP1,N,6
        CDASUM = CDASUM + ABS(DX(I  )) + ABS(DX(I+1)) + ABS(DX(I+2))
     &                  + ABS(DX(I+3)) + ABS(DX(I+4)) + ABS(DX(I+5))
      ENDDO
C
      RETURN
      END

