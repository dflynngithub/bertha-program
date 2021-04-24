      FUNCTION UEHINT0(N,ZETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       UU    UU EEEEEEEE HH    HH IIII NN    NN TTTTTTTT 000000       C
C       UU    UU EE       HH    HH  II  NNN   NN    TT   00    00      C
C       UU    UU EE       HH    HH  II  NNNN  NN    TT   00    00      C
C       UU    UU EEEEEE   HHHHHHHH  II  NN NN NN    TT   00    00      C
C       UU    UU EE       HH    HH  II  NN  NNNN    TT   00    00      C
C       UU    UU EE       HH    HH  II  NN   NNN    TT   00    00      C
C        UUUUUU  EEEEEEEE HH    HH IIII NN    NN    TT    000000       C
C                                                                      C
C -------------------------------------------------------------------- C
C  UEHINT0 CALCULATES AN ATOMIC UEHLING POTENTIAL OVERLAP INTEGRAL     C
C  FOR A LOCAL GAUSSIAN NUCLEAR CHARGE DISTRIBUTION.                   C
C -------------------------------------------------------------------- C
C  DFNOTE: AT THE MOMENT THIS IS JUST A COPY OF BNUCINT0.              C
C**********************************************************************C
      PARAMETER(MDM=2500,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26)
C
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     ROUTINE ONLY ALLOWS ODD PARAMETERS N
      IF(MOD(N,2).NE.1) THEN
        WRITE(6, *) 'In UEHINT0: order N must be odd. N = ',N
        WRITE(7, *) 'In UEHINT0: order N must be odd. N = ',N
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
        TRM = TRM
      ELSEIF(N.EQ.3) THEN
        VA  = 2.0D0 + 3.0D0*X
        TRM = TRM*VA
      ELSEIF(N.EQ.5) THEN
        VA  = 8.0D0 + 20.0D0*X + 15.0D0*X*X
        TRM = TRM*VA
      ELSEIF(N.EQ.7) THEN
        VA  = 16.0D0 + 56.0D0*X + 70.0D0*X*X + 35.0D0*X*X*X
        TRM = 3.0D0*TRM*VA
      ELSEIF(N.EQ.9) THEN
        VA  = 128.0D0 + 576.0D0*X + 1008.0D0*X*X + 840.0D0*X*X*X 
        VB  = 315.0D0*X*X*X*X
        TRM = 3.0D0*TRM*(VA+VB)
      ELSEIF(N.EQ.11) THEN
        VA  = 256.0D0 + 1408.0D0*X + 3168.0D0*X*X + 3696.0D0*X*X*X 
        VB  = 2310.0D0*X*X*X*X + 693.0D0*X*X*X*X*X
        TRM = 15.0D0*TRM*(VA+VB)
      ELSEIF(N.EQ.13) THEN
        VA  = 1024.0D0 + 6656.0D0*X + 18304.0D0*X*X 
        VB  = 27456.0D0*X*X*X + 24024.0D0*X*X*X*X 
        VC  = 12012.0D0*X5 + 3003.0D0*X5*X
        TRM = 45.0D0*TRM*(VA+VB+VC)
      ELSEIF(N.EQ.15) THEN
        VA  = 2048.0D0 + 15360.0D0*X + 49920.0D0*X*X 
        VB  = 91520.0D0*X*X*X+ 102960.0D0*X*X*X*X + 72072.0D0*X5 
        VC  = 30030.0D0*X5*X + 6435.0D0*X5*X*X
        TRM = 315.0D0*TRM*(VA+VB+VC)
      ELSEIF(N.EQ.17) THEN
        VA  = 32768.0D0 + 278528.0D0*X + 1044480.0D0*X*X 
        VB  = 2263040.0D0*X*X*X + 3111680.0D0*X*X*X*X 
        VC  = 2800512.0D0*X5 + 1633632.0D0*X5*X + 583440.0D0*X5*X*X
        VD  = 109395.0D0*X5*X*X*X
        TRM = 315.0D0*TRM*(VA+VB+VC+VD)
      ELSEIF(N.EQ.19) THEN
        VA  = 65536.0D0 + 6222592.0D0*X + 2646016.0D0*X*X 
        VB  = 6615040.0D0*X*X*X + 10749440.0D0*X*X*X*X 
        VC  = 11824384.0D0*X5 + 8868288.0D0*X5*X + 4434144.0D0*X5*X*X 
        VD  = 1385670.0D0*X5*X*X*X + 230945.0D0*X5*X*X*X*X
        TRM = 2835.0D0*TRM*(VA+VB+VC+VD)
      ELSEIF(N.EQ.21) THEN
        VA  = 262144.0D0 + 2752512.0D0*X + 13074432.0D0*X*X 
        VB  = 37044224.0D0*X*X*X + 69457920.0D0*X*X*X*X 
        VC  = 90295296.0D0*X5 + 82770688.0D0*X5*X 
        VD  = 53209728.0D0*X5*X*X + 23279256.0D0*X5*X*X*X 
        VE  = 6466460.0D0*X5*X*X*X*X + 969969.0D0*X5*X5
        TRM = 14175.0D0*TRM*(VA+VB+VC+VD+VE)
      ELSE 
        WRITE(6, *) 'In UEHINT0: order N too large. N = ',N
        WRITE(7, *) 'In UEHINT0: order N too large. N = ',N
      ENDIF
C
C     TRANSFER DATA TO UEHINT0
      UEHINT0 = TRM
C
      RETURN
      END
