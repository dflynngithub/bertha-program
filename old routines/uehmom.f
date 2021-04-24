      FUNCTION UEHMOM(IZ,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     UU    UU EEEEEEEE HH    HH MM       MM  OOOOOO  MM       MM      C
C     UU    UU EE       HH    HH MMM     MMM OO    OO MMM     MMM      C
C     UU    UU EE       HH    HH MMMM   MMMM OO    OO MMMM   MMMM      C
C     UU    UU EEEEEE   HHHHHHHH MM MM MM MM OO    OO MM MM MM MM      C
C     UU    UU EE       HH    HH MM  MMM  MM OO    OO MM  MMM  MM      C
C     UU    UU EE       HH    HH MM   M   MM OO    OO MM   M   MM      C
C      UUUUUU  EEEEEEEE HH    HH MM       MM  OOOOOO  MM       MM      C
C                                                                      C
C -------------------------------------------------------------------- C
C  UEHMOM EVALUATES THE KTH MOMENT OF THE POLARISED NUCLEAR CHARGE     C
C  DENSITY DUE TO THE UEHLING INTERACTION.                             C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5  NMDL
C
      DIMENSION C0(7)
      DIMENSION RHO(0:NSRC)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,MFT),XNUC(MCT,MFT),NNUC(MCT),NMDL(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     MULTIPLIER FOR COMPTON WAVELENGTH (USE FOR MUON OR TAUON FIELD)
      CMPF = 1.0D0*CMPW
      CMP2 = CMPF*CMPF
C
C     SQUARE OF RMS RADIUS
      RN2 = RNUC(IZ)*RNUC(IZ)
C
C     CHECK THAT SPECIFIED MOMENT IS SUPPORTED BY THIS FUNCTION
      IF(K.LT.-1.OR.K.GT.6) THEN
        WRITE(6, *) 'In UEHMOM: cannot support this moment! K=',K
        WRITE(7, *) 'In UEHMOM: cannot support this moment! K=',K
        UEHMOM = 0.0D0
        RETURN
      ENDIF
C
C     CASE FOR EVEN MOMENTS
      IF(MOD(K,2).EQ.0) THEN
        IF(K.EQ.0) THEN
          UEHMOM = 0.0D0
        ELSEIF(K.EQ.2) THEN
          UEHMOM = 2.0D0*CMP2/(5.0D0*PI*CV)
        ELSEIF(K.EQ.4) THEN
          UEHMOM = 4.0D0*CMP2*RN2/(3.0D0*PI*CV)
     &           + 6.0D0*CMP2*CMP2/(7.0D0*PI*CV)
        ELSEIF(K.EQ.6) THEN
          UEHMOM = 14.0D0*CMP2*RNUCMOM(IZ,4)/(5.0D0*PI*CV)
     &           +  6.0D0*CMP2*CMP2*RN2/(PI*CV)
     &           + 16.0D0*CMP2*CMP2*CMP2/(3.0D0*PI*CV)
        ENDIF
        RETURN
      ENDIF
C
C     CASE FOR ODD MOMENTS
C
C     CHI_N(0) VALUES
      C0(2) = 9.0D0*PI/32.0D0
      C0(3) = 0.4D0
      C0(4) = 5.0D0*PI/64.0D0
      C0(5) = 6.0D0/35.0D0
      C0(6) = 21.0D0*PI/512.0D0
      C0(7) = 32.0D0/315.0D0
C
C     SEARCH FOR CHARGE RADIUS SMAX FOR WHICH RHO(SMAX) < 1.0D-16
      SMAX = 0.0D0
60    SMAX = SMAX + 0.1D0/CFM
      PMAX = SMAX*RHONUC(NMDL(IZ),IZ,SMAX)
      IF(DABS(PMAX).GT.1.0D-16) GOTO 60
C
C     SOURCE CHARGE STEP SIZE
      HS = SMAX/DFLOAT(NSRC)
C
C     EVALUATE CHARGE DENSITY ON UNIFORMLY-SPACED GRID
      DO M=0,NSRC
        S = HS*DFLOAT(M)
        RHO(M) = RHONUC(NMDL(IZ),IZ,S)
      ENDDO
C
      IF(K.EQ.-1) THEN
C                 ~
C       EVALUATE <R^-1>
        R1R = 0.0D0
        DO M=1,NSRC
C
C         SOURCE CHARGE RADIUS
          S  = HS*DFLOAT(M)
C
C         INTEGRAND VALUE
          X0 = 2.0D0*S/CMPF
          R1R = R1R + EXTINT11(S*RHO(M)*CHIFNC(X0,1,0),M,NSRC)
C
        ENDDO
C
C       MULTIPLICATIVE FACTORS
        R1R = 5.0D0*HS*R1R/2.99376D+5
        R1R = 4.0D0*PI*R1R
C
C       MULTIPLICATIVE FACTORS FOR VUEH(0)
        UEHMOM =-2.0D0*R1R/(3.0D0*PI*CV)
C
      ELSEIF(K.EQ.1) THEN
C                 ~
C       EVALUATE <R^-1>
        R3R = 0.0D0
        DO M=1,NSRC
C
C         SOURCE CHARGE RADIUS
          S = HS*DFLOAT(M)
C
C         INTEGRAND VALUE
          X0  = 2.0D0*S/CMPF
          BR3 = C0(3)-CHIFNCQUAD(X0,3)
          R3R = R3R + EXTINT11(S*RHO(M)*BR3,M,NSRC)
C
        ENDDO
C
C       MULTIPLICATIVE FACTORS
        R3R = 5.0D0*HS*R3R/2.99376D+5
C
C       MULTIPLICATIVE FACTORS FOR VUEH(0)
        UEHMOM = 4.0D0*CMP2*R3R/(3.0D0*CV)
C
      ELSEIF(K.EQ.3) THEN
C                 ~          ~
C       EVALUATE <R^-1> AND <R^+1>
        R5R = 0.0D0
        R1F = 0.0D0
        DO M=1,NSRC
C
C         SOURCE CHARGE RADIUS
          S = HS*DFLOAT(M)
C
C         INTEGRAND VALUE
          X0  = 2.0D0*S/CMPF
          BR5 = 1.0D0 - CHIFNC(X0,5,0)/C0(5)
          R5R = R5R + EXTINT11(S*RHO(M)*BR5,M,NSRC)
          R1F = R1F + EXTINT11(S*RHO(M)*CHIFNC(X0,1,1),M,NSRC)
C
        ENDDO
C
C       MULTIPLICATIVE FACTORS
        R5R = 5.0D0*HS*R5R/2.99376D+5
        R5R = 4.0D0*PI*R5R
        R1F = 5.0D0*HS*R1F/2.99376D+5
        R1F = 4.0D0*PI*R1F
C
C       MULTIPLICATIVE FACTORS FOR VUEH(0)
        UEHMOM = 9.0D0*CMP2*CMP2*R5R/(35.0D0*PI*CV)
     &         + 4.0D0*CMP2*R1F/(5.0D0*PI*CV)
C
      ELSEIF(K.EQ.5) THEN
C                 ~       ~          ~
C       EVALUATE <R^-1>, <R^+1> AND <R^+3>
        R7R = 0.0D0
        R1F = 0.0D0
        R3F = 0.0D0
        DO M=1,NSRC
C
C         SOURCE CHARGE RADIUS
          S = HS*DFLOAT(M)
C
C         INTEGRAND VALUE
          X0  = 2.0D0*S/CMPF
          BR7 = 1.0D0 - CHIFNC(X0,7,0)/C0(7)
          R7R = R7R + EXTINT11(S*RHO(M)*BR7,M,NSRC)
          R1F = R1F + EXTINT11(S*RHO(M)*CHIFNC(X0,1,1),M,NSRC)
          R3F = R3F + EXTINT11(S*RHO(M)*CHIFNC(X0,3,1),M,NSRC)
C
        ENDDO
C
C       MULTIPLICATIVE FACTORS
        R7R = 5.0D0*HS*R7R/2.99376D+5
        R7R = 4.0D0*PI*R7R
        R1F = 5.0D0*HS*R1F/2.99376D+5
        R1F = 4.0D0*PI*R1F
        R3F = 5.0D0*HS*R3F/2.99376D+5
        R3F = 4.0D0*PI*R3F
C
C       MULTIPLICATIVE FACTORS FOR VUEH(0)
        UEHMOM = 16.0D0*CMP2*CMP2*CMP2*R7R/(21.0D0*PI*CV)
     &         + 18.0D0*CMP2*CMP2*R1F/(7.0D0*PI*CV)
     &         +  2.0D0*CMP2*R3F/(PI*CV)
C
      ENDIF
C
      RETURN
      END
