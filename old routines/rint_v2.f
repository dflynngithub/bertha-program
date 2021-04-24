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
C     ROUTINE ONLY ALLOWS ODD PARAMETERS N
      IF(MOD(N,2).NE.1) THEN
        WRITE(6, *) 'In RINT: parameter N must be odd. N = ',N
        WRITE(7, *) 'In RINT: parameter N must be odd. N = ',N
      ENDIF
C
C     FACTORS NEEDED FOR ALL PARAMETERS N
      X   = ZETA/PNUC
      X5  = X*X*X*X*X
      T0  = PNUC+ZETA
      RAT = PNUC/T0
      TRM = 0.5D0*DSQRT(PNUC)/ZETA/DSQRT(T0)
      DO I=1,(N-1)/2
        TRM = 0.5D0*TRM*RAT/ZETA
      ENDDO
C
      IF(N.EQ.1) THEN
        RINT = TRM
      ELSEIF(N.EQ.3) THEN
        POLY = 2.0D0 + 3.0D0*X
        RINT = TRM*POLY
      ELSEIF(N.EQ.5) THEN
        POLY = 8.0D0 + 20.0D0*X + 15.0D0*X*X
        RINT = TRM*POLY
      ELSEIF(N.EQ.7) THEN
        POLY = 16.0D0 + 56.0D0*X + 70.0D0*X*X + 35.0D0*X*X*X
        RINT = 3.0D0*TRM*POLY
      ELSEIF(N.EQ.9) THEN
        POLY = 128.0D0 + 576.0D0*X + 1008.0D0*X*X + 840.0D0*X*X*X 
     &       + 315.0D0*X*X*X*X
        RINT = 3.0D0*TRM*POLY
      ELSEIF(N.EQ.11) THEN
        POLY = 256.0D0 + 1408.0D0*X + 3168.0D0*X*X + 3696.0D0*X*X*X 
     &       + 2310.0D0*X*X*X*X + 693.0D0*X*X*X*X*X
        RINT = 15.0D0*TRM*POLY
      ELSEIF(N.EQ.13) THEN
        POLY = 1024.0D0 + 6656.0D0*X + 18304.0D0*X*X 
     &       + 27456.0D0*X*X*X + 24024.0D0*X*X*X*X 
     &       + 12012.0D0*X5 + 3003.0D0*X5*X
        RINT = 45.0D0*TRM*POLY
      ELSEIF(N.EQ.15) THEN
        POLY = 2048.0D0 + 15360.0D0*X + 49920.0D0*X*X 
     &       + 91520.0D0*X*X*X+ 102960.0D0*X*X*X*X + 72072.0D0*X5 
     &       + 30030.0D0*X5*X + 6435.0D0*X5*X*X
        RINT = 315.0D0*TRM*POLY
      ELSEIF(N.EQ.17) THEN
        POLY = 32768.0D0 + 278528.0D0*X + 1044480.0D0*X*X 
     &       + 2263040.0D0*X*X*X + 3111680.0D0*X*X*X*X 
     &       + 2800512.0D0*X5 + 1633632.0D0*X5*X + 583440.0D0*X5*X*X
     &       + 109395.0D0*X5*X*X*X
        RINT = 315.0D0*TRM*POLY
      ELSEIF(N.EQ.19) THEN
        POLY = 65536.0D0 + 6222592.0D0*X + 2646016.0D0*X*X 
     &       + 6615040.0D0*X*X*X + 10749440.0D0*X*X*X*X 
     &       + 11824384.0D0*X5 + 8868288.0D0*X5*X + 4434144.0D0*X5*X*X 
     &       + 1385670.0D0*X5*X*X*X + 230945.0D0*X5*X*X*X*X
        RINT = 2835.0D0*TRM*POLY
      ELSEIF(N.EQ.21) THEN
        POLY = 262144.0D0 + 2752512.0D0*X + 13074432.0D0*X*X 
     &       + 37044224.0D0*X*X*X + 69457920.0D0*X*X*X*X 
     &       + 90295296.0D0*X5 + 82770688.0D0*X5*X 
     &       + 53209728.0D0*X5*X*X + 23279256.0D0*X5*X*X*X 
     &       + 6466460.0D0*X5*X*X*X*X + 969969.0D0*X5*X5
        RINT = 14175.0D0*TRM*POLY
      ELSE 
        WRITE(6, *) 'In RINT: parameter N too large. N = ',N
        WRITE(7, *) 'In RINT: parameter N too large. N = ',N
      ENDIF
C
      RETURN
      END
