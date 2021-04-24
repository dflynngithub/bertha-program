      SUBROUTINE DENSTY0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    DDDDDDD  EEEEEEEE NN    NN  SSSSSS TTTTTTTT YY    YY  000000      C
C    DD    DD EE       NNN   NN SS    SS   TT    YY    YY 00    00     C
C    DD    DD EE       NNNN  NN SS         TT     YY  YY  00    00     C
C    DD    DD EEEEEE   NN NN NN  SSSSSS    TT      YYYY   00    00     C
C    DD    DD EE       NN  NNNN       SS   TT       YY    00    00     C
C    DD    DD EE       NN   NNN SS    SS   TT       YY    00    00     C
C    DDDDDDD  EEEEEEEE NN    NN  SSSSSS    TT       YY     000000      C
C                                                                      C
C -------------------------------------------------------------------- C
C  DENSTY0 IS A STARTING DENSITY ROUTINE FOR USE ONLY WHEN ITREE=1.    C
C**********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26)
C
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM),C(MDM,MDM)
      COMPLEX*16 SUM
C
      COMMON/COEF/C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/SHLL/ALPHA,BETA,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR,IOCCM0
C
C     CONSTRUCT THE CLOSED-SHELL AND TOTAL DENSITY ARRAYS
      DO I=1,NDIM
        DO J=1,NDIM
          SUM = DCMPLX(0.0D0,0.0D0)
          DO IOCC=1,IOCCM0
            SUM = SUM + DCONJG(C(I,IOCC+NSHIFT))*C(J,IOCC+NSHIFT)
          ENDDO
          DENC(I,J) = SUM
          DENO(I,J) = 0.0D0
          DENT(I,J) = DENC(I,J) + DENO(I,J)
        ENDDO
      ENDDO
C     
C     UPDATE OCCUPATION COUNTER
      IOCCM0 = NOCC
C
      RETURN
      END
