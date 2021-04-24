      SUBROUTINE DENSTY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         DDDDDDD  EEEEEEEE NN    NN  SSSSSS TTTTTTTT YY    YY         C
C         DD    DD EE       NNN   NN SS    SS   TT    YY    YY         C
C         DD    DD EE       NNNN  NN SS         TT     YY  YY          C
C         DD    DD EEEEEE   NN NN NN  SSSSSS    TT      YYYY           C
C         DD    DD EE       NN  NNNN       SS   TT       YY            C
C         DD    DD EE       NN   NNN SS    SS   TT       YY            C
C         DDDDDDD  EEEEEEEE NN    NN  SSSSSS    TT       YY            C
C                                                                      C
C -------------------------------------------------------------------- C
C  DENSTY GENERATES DENSITY MATRICES FROM THE EXPANSION COEFFS C.      C
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
C     MAKE CLOSED-SHELL DENSITY AND EMPTY OPEN-SHELL DENSITY (RSCF 81)
      DO I=1,NDIM
        DO J=1,NDIM
          SUM = DCMPLX(0.0D0,0.0D0)
          DO IOCC=1,NCLS
            ICL = ICLS(IOCC)
            SUM = SUM + DCONJG(C(I,ICL+NSHIFT))*C(J,ICL+NSHIFT)
          ENDDO
          DENC(I,J) = SUM
          DENO(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     MAKE THE OPEN-SHELL DENSITY (RSCF 82)
      IF(NOPN.NE.0) THEN
        DO I=1,NDIM
          DO J=1,NDIM
            SUM = DCMPLX(0.0D0,0.0D0)
            DO IOCC=1,NOPN
              IOP = IOPN(IOCC)
              SUM = SUM + FOPEN*DCONJG(C(I,IOP+NSHIFT))*C(J,IOP+NSHIFT)
            ENDDO
            DENO(I,J) = SUM
          ENDDO
        ENDDO
      ENDIF
C
C     MAKE THE TOTAL DENSITY MATRIX (RSCF 83)
      DO I=1,NDIM
        DO J=1,NDIM
          DENT(I,J) = DENC(I,J) + DENO(I,J)
        ENDDO
      ENDDO
C     
      RETURN
      END
