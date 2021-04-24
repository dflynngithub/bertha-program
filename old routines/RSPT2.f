      SUBROUTINE RSPT2(V1,V2,E3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C              RRRRRRR   SSSSSS  PPPPPPP TTTTTTTT 222222               C
C              RR    RR SS    SS PP    PP   TT   22    22              C
C              RR    RR SS       PP    PP   TT         22              C
C              RR    RR  SSSSSS  PP    PP   TT       22                C
C              RRRRRRR        SS PPPPPPP    TT     22                  C
C              RR    RR SS    SS PP         TT   22                    C
C              RR    RR  SSSSSS  PP         TT   22222222              C
C                                                                      C
C -------------------------------------------------------------------- C
C  RSPT2 APPLIES 2ND ORDER RAYLEIGH-SCHRODINGER PERTURBATION THEORY TO C
C  A SET OF 1ST ORDER MATRIX ELEMENTS IN V1, AND OUTPUTS RESULTS TO V2.C
C -------------------------------------------------------------------- C
C  DUE TO THE 2-FOLD MQN SYMMETRY, EITHER APPLY DEGENERATE RSPT OR     C
C  SIMPLY IGNORE THE RELEVANT PAIR ORBITAL.                            C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 V1(MDM,MDM),V2(MDM,MDM),E3(MDM)
C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
C
C     TOLERANCE VALUE FOR VANISHING MATRIX ELEMENTS
      TOL = 1.0D-10
C
C     LOOP OVER ALL ORBITAL COMBINATIONS (NEGATIVE AND POSITIVE ENERGY)
      DO IOCC=1,NDIM
        DO JOCC=1,NDIM
C
C         APPLY 2ND ORDER RSPT
          V1(IOCC,JOCC) = DCMPLX(0.0D0,0.0D0)
C
C         VANISHING MATRIX ELEMENTS
          TMP1 = DREAL(V2(IOCC,JOCC))
          TMP2 = DIMAG(V2(IOCC,JOCC))
          IF(DABS(TMP1).LT.TOL) THEN
            TMP1 = 0.0D0
          ENDIF
          IF(DABS(TMP2).LT.TOL) THEN
            TMP2 = 0.0D0
          ENDIF
          V2(IOCC,JOCC) = DCMPLX(TMP1,TMP2)
C
        ENDDO
      ENDDO
C
C     THIRD ORDER ENERGY CORRECTION
      DO IOCC=1,NDIM
        E3(IOCC) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
      RETURN
      END
