      FUNCTION RINT(N,ZETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                   RRRRRRR  IIII NN    NN TTTTTTTT                    C
C                   RR    RR  II  NNN   NN    TT                       C
C                   RR    RR  II  NNNN  NN    TT                       C
C                   RR    RR  II  NN NN NN    TT                       C
C                   RRRRRRR   II  NN  NNNN    TT                       C
C                   RR    RR  II  NN   NNN    TT                       C
C                   RR    RR IIII NN    NN    TT                       C
C                                                                      C
C -------------------------------------------------------------------- C
C  RINT CALCULATES THE INTEGRALS REQUIRED FOR ATOMIC BASIS FUNCTION    C
C  OVERLAPS WITH A LOCAL NUCLEUS WHOSE CHARGE DISTRIBUTION IS GAUSSIAN.C
C -------------------------------------------------------------------- C
C  RINT(N,ZETA) = INT{R^N*EXP(-ZETA*R^2)*ERF(SQRT(PNUC)*R)}.           C
C**********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26)
C
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR,IOCCM0
C
      IF(N.EQ.1) THEN
        RINT = 0.5D0*DSQRT(PNUC)/ZETA/DSQRT(PNUC+ZETA)
      ELSEIF(N.EQ.3) THEN
        T1   = DSQRT(PNUC)
        T2   = T1*T1
        T4   = ZETA*ZETA
        T8   = DSQRT(PNUC+ZETA)
        T9   = T8*T8
        T11  = 1.0D0/T9/T8
        RINT = 0.50D0/T4*T2*T1*T11 + 0.75D0*T1/ZETA*T11
      ELSEIF(N.EQ.5) THEN
        T1   = DSQRT(PNUC)
        T2   = T1*T1
        T3   = T2*T2
        T5   = ZETA*ZETA
        T10  = DSQRT(PNUC+ZETA)
        T11  = T10*T10
        T12  = T11*T11
        T14  = 1.0D0/T12/T10
        RINT = T3*T1/T5/ZETA*T14 + 2.5D0*T2*T1/T5*T14 
     &       + 1.875D0*T1/ZETA*T14
      ELSEIF(N.EQ.7) THEN
        T1   = DSQRT(PNUC)
        T2   = T1*T1
        T3   = T2*T1
        T4   = T2*T2
        T6   = ZETA*ZETA
        T7   = T6*T6
        T11  = DSQRT(PNUC+ZETA)
        T12  = T11*T11
        T14  = T12*T12
        T16  = 1.0D0/T14/T12/T11
        RINT = 3.0D0*T4*T3/T7*T16 + 1.05D1*T4*T1/T6/ZETA*T16
     &       + 1.3125D1*T3/T6*T16 + 6.5625D0*T1/ZETA*T16
      ELSEIF(N.EQ.9) THEN
        T1   = DSQRT(PNUC)
        T2   = T1*T1
        T3   = T2*T2
        T4   = T3*T3
        T6   = ZETA*ZETA
        T7   = T6*T6
        T12  = DSQRT(PNUC+ZETA)
        T13  = T12*T12
        T14  = T13*T13
        T15  = T14*T14
        T17  = 1.0D0/T15/T12
        T19  = T2*T1
        RINT = 1.2D1*T4*T1/T7/ZETA*T17  + 5.4D1*T3*T19/T7*T17
     &       + 9.45D1*T3*T1/T6/ZETA*T17 + 7.875D3*T19/T6*T17
     &       + 2.953125D1*T1/ZETA*T17
      ELSEIF(N.EQ.11) THEN
        T1   = DSQRT(PNUC)
        T2   = T1*T1
        T3   = T2*T1
        T4   = T2*T2
        T5   = T4*T4
        T7   = ZETA*ZETA
        T8   = T7*T7
        T13  = DSQRT(PNUC+ZETA)
        T14  = T13*T13
        T16  = T14*T14
        T17  = T16*T16
        T19  = 1.0D0/T17/T14/T13
        RINT = 6.0D1*T5*T3/T8/T7*T19 + 3.3D2*T5*T1/T8/ZETA*T19
     &       + 7.425D2*T4*T3/T8*T19  + 8.6625D2*T4*T1/T7/ZETA*T19
     &       + 5.4140625D2*T3/T7*T19 + 1.62421875D2*T1/ZETA*T19
      ELSE 
        WRITE(6, *) 'In RINT: invalad order N. N = ',N
        WRITE(7, *) 'In RINT: invalad order N. N = ',N
      ENDIF
C
      RETURN
      END

