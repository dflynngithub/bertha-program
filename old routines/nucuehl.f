C
C
      SUBROUTINE NUCUEHL(IZ,NFT,NZR,IWRT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    NN    NN UU    UU  CCCCCC  UU    UU EEEEEEEE HH    HH LL          C
C    NNN   NN UU    UU CC    CC UU    UU EE       HH    HH LL          C
C    NNNN  NN UU    UU CC       UU    UU EE       HH    HH LL          C
C    NN NN NN UU    UU CC       UU    UU EEEEEE   HHHHHHHH LL          C
C    NN  NNNN UU    UU CC       UU    UU EE       HH    HH LL          C
C    NN   NNN UU    UU CC    CC UU    UU EE       HH    HH LL          C
C    NN    NN  UUUUUU   CCCCCC   UUUUUU  EEEEEEEE HH    HH LLLLLLLL    C
C                                                                      C
C -------------------------------------------------------------------- C
C  NUCUEHL GENERATES A BEST-FIT GAUSSIAN SET FOR THE UEHLING POTENTIAL C
C  ARISING FROM NUCLEUS IZ, USING NFT TOTAL GAUSSIANS.                 C
C -------------------------------------------------------------------- C
C  MATCHING CRITERIA AVAILABLE FOR GAUSSIAN AMPLITUDES:                C
C   > RADIAL ARGUMENTS SCALED BY NUCLEAR RMS RADIUS RNUC(IZ).          C
C   > UEHLING POTENTIAL FIRST SCALED TO R*R*V(R).                      C
C   > POTENTIAL MATCHED TO NFT POINTS, FROM 3.0D0*RNUC TO 200.0D0*RNUC.C
C   > MATCHING POINTS ARE EXPONENTIALLY SPACED.                        C
C   > GAUSSIAN EXPONENTS FROM MODIFIED GEOMETRIC SERIES WITH THE SAME  C
C     BETA AND GAMMA, BUT ALPHA IS SCALED BY RNUC*RNUC.                C
C   > WHEN THIS IS COMPLETE, POTENTIAL V(R) NEAR R=0.0D0 IS AUGMENTED  C
C     WITH SOME MORE GAUSSIANS WITH LARGER EXPONENTS IF V(R) IS NEEDED.C
C**********************************************************************C
      INCLUDE 'param.h'
      PARAMETER(NPTS=100,NTRY=1000)
C
      CHARACTER*4  HMLT
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(2)
C
      DIMENSION VF(0:NPTS),VG(0:NPTS)
C
      DIMENSION XS(0:NLW),Y1S(0:NLW),D1S(0:NLW),Y2S(0:NLW),D2S(0:NLW),
     &          XB(0:NUP),Y1B(0:NUP),D1B(0:NUP),Y2B(0:NUP),D2B(0:NUP)
C
      DIMENSION X(MFT,MFT),Z(MFT),Y(MFT),IPIV(MFT)
      DIMENSION XUEH(MCT,MFT),FUEH(MCT,MFT)
      DIMENSION XORI(MCT,MFT),FORI(MCT,MFT)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),CFMI(MCT),
     &            FNUC(MCT,MFT),XNUC(MCT,MFT),RNUC(MCT),NNUC(MCT)
      COMMON/CONV/CHZ,CFM,CDB
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/PRMS/HMLT,IOPT,IMOL,INEW,ILEV,ISML,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/UEHL/RAD(0:NRAD),VVAC(MCT,0:NRAD),RORI,RMID,RMAX,NLIN,NEXP
C
C     MULTIPLIER FOR COMPTON WAVELENGTH (USE FOR MUON OR TAUON FIELD)
      CMPF = CMPW/1.0D0
C
C**********************************************************************C
C     GENERATE AND IMPORT THE INTERPOLATED CHI FUNCTIONS               C
C**********************************************************************C
C
C     GENERATE UEHLING ANCILLIARY FUNCTIONS AND WRITE TO FILE
C     CALL CHIGEN
C
C     IMPORT SPLINE DATA FOR CHI_1(X) GIVEN 0 <= X <= XSPL
      OPEN(UNIT=50,FILE='chi1_small.dat',STATUS='UNKNOWN')
      REWIND(UNIT=50)
      READ(50,*) XSPL
      DO N=0,NLW
        READ(50,*) XS(N),Y1S(N),D1S(N)
      ENDDO
      CLOSE(UNIT=50)
C
C     IMPORT SPLINE DATA FOR CHI_1(X) GIVEN X >= XSPL
      OPEN(UNIT=51,FILE='chi1_big.dat',STATUS='UNKNOWN')
      REWIND(UNIT=51)
      READ(51,*) XSPL
      DO N=0,NUP
        READ(51,*) XB(N),Y1B(N),D1B(N)
      ENDDO
      CLOSE(UNIT=51)
C
C     IMPORT SPLINE DATA FOR CHI_2(X) GIVEN 0 <= X <= XSPL
      OPEN(UNIT=52,FILE='chi2_small.dat',STATUS='UNKNOWN')
      REWIND(UNIT=52)
      READ(52,*) XSPL
      DO N=0,NLW
        READ(52,*) XS(N),Y2S(N),D2S(N)
      ENDDO
      CLOSE(UNIT=52)
C
C     IMPORT SPLINE DATA FOR CHI_2(X) GIVEN X >= XSPL
      OPEN(UNIT=53,FILE='chi2_big.dat',STATUS='UNKNOWN')
      REWIND(UNIT=53)
      READ(53,*) XSPL
      DO N=0,NUP
        READ(53,*) XB(N),Y2B(N),D2B(N)
      ENDDO
      CLOSE(UNIT=53)
C
C**********************************************************************C
C     INTEGRATION PARAMETERS AND UEHLING POTENTIAL AT THE ORIGIN       C
C**********************************************************************C
C
C     SEARCH FOR CHARGE RADIUS SMAX FOR WHICH RHO(SMAX) < 1.0D-16
      SMAX = 0.0D0
31    SMAX = SMAX + 0.1D0/CFM
      PMAX = RHONUC(IZ,SMAX)
      IF(DABS(PMAX).GT.1.0D-16) GOTO 31
C
C     SOURCE CHARGE STEP SIZE
      HS = SMAX/DFLOAT(NSRC)
C
C     UEHLING POTENTIAL ORIGIN VALUE
      V0 = 0.0D0
      DO M=1,NSRC
C
C       SOURCE CHARGE RADIUS
        S  = HS*DFLOAT(M)
C
C       CHI_1(X) VALUE (ARGUMENT WILL ALWAYS BE LESS THAN 0.6D0)
        X0 = 2.0D0*S/CMPF
        CALL SPLNINT(XS,Y1S,D1S,X0,C0,NLW)
        C0 = C0*DEXP(-X0)
C
C       CONTRIBUTION TO INTEGRAND
        V0 = V0 + EXTINT11(S*RHONUC(IZ,S)*C0,M,NSRC)
C
      ENDDO
C
C     INTEGRATION WEIGHTING FACTORS
      V0 = 5.0D0*HS*V0/2.99376D+5
      V0 =-8.0D0*V0/(3.0D0*CV)
C
C**********************************************************************C
C     WEIGHTED FERMI POTENTIAL R*R*V(R) -- MATCHING AND GRID VALUES.   C
C**********************************************************************C
C
C     LARGEST RADIUS TO COMPARE AGAINST
      RLIM = 250.0D0
C
C     MATCHING POINT DETAILS
      R0 = 3.0D0
      RN = 200.0D0
      HR = DLOG(RN/R0)/DFLOAT(NFT-1)
C
C     STORE R*R*V(R) ON THE UNIFORM GRID IN VF(0:NPTS)
      VF(0) = 0.0D0
      DO IPTS=1,NPTS
C
C       RADIUS ON UNIFORM SCALE
        R = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C       INITIALISE POTENTIAL VALUE
        VF(IPTS) = 0.0D0
C
C       INTEGRATE OVER CHARGE SOURCE
        DO M=0,NSRC
C
C         SOURCE RADIUS
          S  = HS*DFLOAT(M)
C
C         CHI FUNCTION ARGUMENTS
          XM = 2.0D0*DABS(R-S)/CMPF
          XP = 2.0D0*DABS(R+S)/CMPF
C
C         COMPONENTS OF INTEGRAND
          IF(0.5D0*(XM+XP).LT.XSPL) THEN
            CALL SPLNINT(XS,Y2S,D2S,XM,CM,NLW)
            CALL SPLNINT(XS,Y2S,D2S,XP,CP,NLW)    
          ELSE
            CALL SPLNINT(XB,Y2B,D2B,XM,CM,NUP)
            CALL SPLNINT(XB,Y2B,D2B,XP,CP,NUP)
          ENDIF
C
C         PERFORM THE MAPPING CHI(X) = SPLINE(X)*EXP(-X)
          CM = CM*DEXP(-XM)
          CP = CP*DEXP(-XP)
C
C         CONTRIBUTION TO INTEGRAND
          VF(IPTS) = VF(IPTS)+EXTINT11(S*RHONUC(IZ,S)*(CM-CP),M,NSRC)
C
        ENDDO
C
C       INTEGRATION WEIGHTING FACTORS
        VF(IPTS) = 5.0D0*HS*VF(IPTS)/2.99376D+5
C
C       OTHER FACTORS AND WEIGHTED POTENTIAL R*R*V(R)
        VF(IPTS) =-2.0D0*CMPF*R*VF(IPTS)/(3.0D0*CV)
C
      ENDDO
C
C     STORE R*R*V(R) MATCHING VALUES IN Y(NFT) -- MATRIX EQN SOLUTIONS
      DO IFT=1,NFT
C
C       RADIUS ON AN EXPONENTIAL SCALE
        R = R0*RNUC(IZ)*DEXP(DFLOAT(IFT-1)*HR)
C
C       INITIALISE VALUE
        Y(IFT) = 0.0D0
C
C       INTEGRATE OVER CHARGE SOURCE
        DO M=0,NSRC
C
C         SET SOURCE RADIUS
          S  = HS*DFLOAT(M)
C
C         CHI FUNCTION ARGUMENTS
          XM = 2.0D0*DABS(R-S)/CMPF
          XP = 2.0D0*DABS(R+S)/CMPF
C
C         COMPONENTS OF INTEGRAND
          IF(0.5D0*(XM+XP).LT.XSPL) THEN
            CALL SPLNINT(XS,Y2S,D2S,XM,CM,NLW)
            CALL SPLNINT(XS,Y2S,D2S,XP,CP,NLW)    
          ELSE
            CALL SPLNINT(XB,Y2B,D2B,XM,CM,NUP)
            CALL SPLNINT(XB,Y2B,D2B,XP,CP,NUP)
          ENDIF
C
C         PERFORM THE MAPPING CHI(X) = SPLINE(X)*EXP(-X)
          CM = CM*DEXP(-XM)
          CP = CP*DEXP(-XP)
C
C         CONTRIBUTION TO INTEGRAND
          Y(IFT) = Y(IFT) + EXTINT11(S*RHONUC(IZ,S)*(CM-CP),M,NSRC)
C
        ENDDO
C
C       INTEGRATION WEIGHTING FACTORS
        Y(IFT) = 5.0D0*HS*Y(IFT)/2.99376D+5
C
C       OTHER FACTORS AND WEIGHTED POTENTIAL R*R*V(R)
        Y(IFT) =-2.0D0*CMPF*R*Y(IFT)/(3.0D0*CV)
C
      ENDDO
C
C**********************************************************************C
C     GIVEN MODIFIED GEOMETRIC SET OF NFT PARAMETERS, FIND BEST ALPHA. C
C**********************************************************************C
C
C     GEOMETRIC PARAMETER INFORMATION
      ALPH = 4.00D-5/(RNUC(IZ)**2)
      BETA = 1.60D+0
      GAMA = 0.15D+0
C
C     INITIALISE BEST-FIT R-SQUARED AND INDEX FOR ALPHA
      RSQBIG = 0.0D0
      ITRYBG = 0
C
C     ITERATE OVER ALPHA SEARCH
      DO ITRY=0,NTRY
C
C       STORE GEOMETRIC SET OF PARAMETERS IN XUEH
        XI = ALPH*(1.0D0 + DFLOAT(ITRY)/DFLOAT(NTRY))
        DO IFT=1,NFT
          XUEH(IZ,IFT) = XI
          RT = DFLOAT(IFT)/DFLOAT(NFT+1)
          XI = BETA*XI*(1.0D0+GAMA*RT*RT)
        ENDDO
C
C       TRANSFER MATCHING VALUES R*R*V(R) INTO Z VECTOR
        DO IFT=1,NFT
          Z(IFT) = Y(IFT)
        ENDDO
C
C
C       SET UP MATRIX EQUATIONS FOR SET OF FITTING GAUSSIANS
        DO IFT=1,NFT
C
C         RADIUS ON AN EXPONENTIAL SCALE
          R = R0*RNUC(IZ)*DEXP(DFLOAT(IFT-1)*HR)
C
C         GAUSSIAN POTENTIALS AT THIS RADIUS
          DO JFT=1,NFT
            X(IFT,JFT) = R*R*DEXP(-XUEH(IZ,JFT)*R*R)
          ENDDO
C
        ENDDO
C
C       SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
        CALL DGESV(NFT,1,X,NFT,IPIV,Z,NFT,INFO)
C
C       TRANSFER AMPLITUDES VALUES INTO FUEH ARRAY
        DO IFT=1,NFT
          FUEH(IZ,IFT) = Z(IFT)
        ENDDO
C
C       BEST-FIT GAUSSIAN POTENTIAL VALUES ACROSS UNIFORM GRID
        DO IPTS=0,NPTS
C
C         RADIUS ON UNIFORM SCALE
          R = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C         INITIALISE BEST-FIT POTENTIAL AT THIS RADIUS
          V = 0.0D0
          DO IFT=1,NFT
            V = V + FUEH(IZ,IFT)*R*R*DEXP(-XUEH(IZ,IFT)*R*R)
          ENDDO
          VG(IPTS) = V
C
        ENDDO
C
C       AVERAGE POTENTIAL VALUE
        YB = 0.0D0
        DO IPTS=0,NPTS
          YB = YB + VG(IPTS)
        ENDDO
        YB = YB/DFLOAT(NPTS+1)
C
C       PREPARATION FOR R-SQUARED VALLUE
        SRES = 0.0D0
        STOT = 0.0D0
        DO IPTS=0,NPTS
          SRES = SRES + (VG(IPTS)-VF(IPTS))**2
          STOT = STOT + (VG(IPTS)-YB      )**2
        ENDDO
C
C       R-SQUARED VALUE AND BEST-FIT CHECK
        RSQ = 1.0D0 - SRES/STOT
        IF(RSQ.GT.RSQBIG) THEN
          ITRYBG = ITRY
          RSQBIG = RSQ
        ENDIF
C
      ENDDO
C
C     IF NO SOLUTION CAN BE FOUND, ALERT THE USER
      IF(RSQBIG.LT.0.0D0) THEN
        WRITE(6, *) 'In NUCUEHL: best-fit search failed. IZ =',IZ
        WRITE(7, *) 'In NUCUEHL: best-fit search failed. IZ =',IZ
        STOP
      ENDIF
C
C**********************************************************************C
C     BEST-FIT ALPHA HAS NOW BEEN FOUND -- RESTORE THESE RESULTS.      C
C**********************************************************************C
C
C     STORE GEOMETRIC SET OF PARAMETERS IN XUEH
      XI = ALPH*(1.0D0 + DFLOAT(ITRYBG)/DFLOAT(NTRY))
      DO IFT=1,NFT
        XUEH(IZ,IFT) = XI
        RT = DFLOAT(IFT)/DFLOAT(NFT+1)
        XI = BETA*XI*(1.0D0+GAMA*RT*RT)
      ENDDO
C
C     TRANSFER MATCHING VALUES R*R*V(R) INTO Z VECTOR
      DO IFT=1,NFT
        Z(IFT) = Y(IFT)
      ENDDO
C
C     SET UP MATRIX EQUATIONS FOR SET OF FITTING GAUSSIANS
      DO IFT=1,NFT
C
C       RADIUS ON AN EXPONENTIAL SCALE
        R = R0*RNUC(IZ)*DEXP(DFLOAT(IFT-1)*HR)
C
C       GAUSSIAN POTENTIALS AT THIS RADIUS
        DO JFT=1,NFT
          X(IFT,JFT) = R*R*DEXP(-XUEH(IZ,JFT)*R*R)
        ENDDO
C
      ENDDO
C
C     SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
      CALL DGESV(NFT,1,X,NFT,IPIV,Z,NFT,INFO)
C
C     TRANSFER AMPLITUDES VALUES INTO FUEH ARRAY
      DO IFT=1,NFT
        FUEH(IZ,IFT) = Z(IFT)
      ENDDO
C
C     BEST-FIT GAUSSIAN POTENTIAL VALUES ACROSS UNIFORM GRID
      DO IPTS=0,NPTS
C
C       RADIUS ON UNIFORM SCALE
        R = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C       INITIALISE BEST-FIT POTENTIAL AT THIS RADIUS
        V = 0.0D0
        DO IFT=1,NFT
          V = V + FUEH(IZ,IFT)*R*R*DEXP(-XUEH(IZ,IFT)*R*R)
        ENDDO
        VG(IPTS) = V
C
      ENDDO
C
      IF(NZR.EQ.0) GOTO 25
C
C**********************************************************************C
C     FERMI POTENTIAL V(R) VALUES -- MATCHING AND GRID VALUES.         C
C**********************************************************************C
C
C     LARGEST RADIUS TO COMPARE AGAINST
      RLIM = 3.0D0
C
C     MATCHING POINT DETAILS (DON'T INCLUDE R0 FROM LAST SECTION)
      R0 = 0.0D0
      RN = 3.0D0
      HR = (RN-R0)/DFLOAT(NZR)
C
C     STORE V(R) ON THE UNIFORM GRID IN VF(0:NPTS)
      VF(0) = V0
C
C     TAKE DIFFERENCE BETWEEN THIS AND EXISTING SOLUTION
      DO IFT=1,NFT
        VF(0) = VF(0) - FUEH(IZ,IFT)
      ENDDO

      DO IPTS=1,NPTS
C
C       RADIUS ON UNIFORM SCALE
        
        R = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C       INITIALISE POTENTIAL VALUE
        VF(IPTS) = 0.0D0
C
C       INTEGRATE OVER CHARGE SOURCE
        DO M=0,NSRC
C
C         SOURCE RADIUS
          S  = HS*DFLOAT(M)
C
C         CHI FUNCTION ARGUMENTS
          XM = 2.0D0*DABS(R-S)/CMPF
          XP = 2.0D0*DABS(R+S)/CMPF
C
C         COMPONENTS OF INTEGRAND
          IF(0.5D0*(XM+XP).LT.XSPL) THEN
            CALL SPLNINT(XS,Y2S,D2S,XM,CM,NLW)
            CALL SPLNINT(XS,Y2S,D2S,XP,CP,NLW)    
          ELSE
            CALL SPLNINT(XB,Y2B,D2B,XM,CM,NUP)
            CALL SPLNINT(XB,Y2B,D2B,XP,CP,NUP)
          ENDIF
C
C         PERFORM THE MAPPING CHI(X) = SPLINE(X)*EXP(-X)
          CM = CM*DEXP(-XM)
          CP = CP*DEXP(-XP)
C
C         CONTRIBUTION TO INTEGRAND
          VF(IPTS) = VF(IPTS)+EXTINT11(S*RHONUC(IZ,S)*(CM-CP),M,NSRC)
C
        ENDDO
C
C       INTEGRATION WEIGHTING FACTORS
        VF(IPTS) = 5.0D0*HS*VF(IPTS)/2.99376D+5
C
C       OTHER FACTORS AND POTENTIAL V(R)
        VF(IPTS) =-2.0D0*CMPF*VF(IPTS)/(3.0D0*CV*R)
C
C       TAKE DIFFERENCE BETWEEN THIS AND EXISTING SOLUTION
        DO IFT=1,NFT
          VF(IPTS) = VF(IPTS) - FUEH(IZ,IFT)*DEXP(-XUEH(IZ,IFT)*R*R)
        ENDDO
C
      ENDDO
C
C     STORE V(R) MATCHING VALUES IN Y(NZR) -- MATRIX EQN SOLUTIONS
      Y(1) = V0
C
C     TAKE DIFFERENCE BETWEEN THIS AND EXISTING SOLUTION
      DO IFT=1,NFT
        Y(1) = Y(1) - FUEH(IZ,IFT)
      ENDDO

      DO IZR=2,NZR
C
C       RADIUS ON A UNIFORM SCALE
        R = RNUC(IZ)*(R0 + DFLOAT(IZR-1)*HR)
C
C       INITIALISE VALUE
        Y(IZR) = 0.0D0
C
C       INTEGRATE OVER CHARGE SOURCE
        DO M=0,NSRC
C
C         SET SOURCE RADIUS
          S  = HS*DFLOAT(M)
C
C         CHI FUNCTION ARGUMENTS
          XM = 2.0D0*DABS(R-S)/CMPF
          XP = 2.0D0*DABS(R+S)/CMPF
C
C         COMPONENTS OF INTEGRAND
          IF(0.5D0*(XM+XP).LT.XSPL) THEN
            CALL SPLNINT(XS,Y2S,D2S,XM,CM,NLW)
            CALL SPLNINT(XS,Y2S,D2S,XP,CP,NLW)    
          ELSE
            CALL SPLNINT(XB,Y2B,D2B,XM,CM,NUP)
            CALL SPLNINT(XB,Y2B,D2B,XP,CP,NUP)
          ENDIF
C
C         PERFORM THE MAPPING CHI(X) = SPLINE(X)*EXP(-X)
          CM = CM*DEXP(-XM)
          CP = CP*DEXP(-XP)
C
C         CONTRIBUTION TO INTEGRAND
          Y(IZR) = Y(IZR) + EXTINT11(S*RHONUC(IZ,S)*(CM-CP),M,NSRC)
C
        ENDDO
C
C       INTEGRATION WEIGHTING FACTORS
        Y(IZR) = 5.0D0*HS*Y(IZR)/2.99376D+5
C
C       OTHER FACTORS AND POTENTIAL V(R)
        Y(IZR) =-2.0D0*CMPF*Y(IZR)/(3.0D0*CV*R)
C
C       TAKE DIFFERENCE BETWEEN THIS AND EXISTING SOLUTION
        DO IFT=1,NFT
          Y(IZR) = Y(IZR) - FUEH(IZ,IFT)*DEXP(-XUEH(IZ,IFT)*R*R)
        ENDDO
C
      ENDDO
      
      DO IPTS=0,NPTS
        VF(IPTS) =-VF(IPTS)
      ENDDO
      DO IZR=1,NZR
        Y(IZR) =-Y(IZR)
      ENDDO
C
C**********************************************************************C
C     GIVEN GEOMETRIC SET OF NZR PARAMETERS, FIND BEST ALPHA.          C
C**********************************************************************C
C
C     GEOMETRIC PARAMETER INFORMATION
      ALPH = TWLG/(0.2D0*RLIM*RNUC(IZ))**2
      BETA = 1.60D+0
C
C     INITIALISE BEST-FIT R-SQUARED AND INDEX FOR ALPHA
      RSQBIG = 0.0D0
      ITRYBG = 0
C
C     ITERATE OVER ALPHA SEARCH
      DO ITRY=0,NTRY
C
C       STORE GEOMETRIC SET OF PARAMETERS IN XORI
        XI = ALPH*(1.0D0 + DFLOAT(ITRY)/DFLOAT(NTRY))
        DO IZR=1,NZR
          XORI(IZ,IZR) = XI
          XI = BETA*XI
        ENDDO
C
C       TRANSFER MATCHING VALUES V(R) INTO Z VECTOR
        DO IZR=1,NZR
          Z(IZR) = Y(IZR)
        ENDDO
C
C       SET UP MATRIX EQUATIONS FOR SET OF FITTING GAUSSIANS
        DO IZR=1,NZR
C
C         RADIUS ON A UNIFORM SCALE
          R = RNUC(IZ)*(R0 + DFLOAT(IZR-1)*HR)
C
C         GAUSSIAN POTENTIALS AT THIS RADIUS
          DO JFT=1,NZR
            X(IZR,JFT) = DEXP(-XORI(IZ,JFT)*R*R)
          ENDDO
C
        ENDDO
C
C       SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
        CALL DGESV(NZR,1,X,NZR,IPIV,Z,NZR,INFO)
C
C       TRANSFER AMPLITUDES VALUES INTO FORI ARRAY
        DO IZR=1,NZR
          FORI(IZ,IZR) = Z(IZR)
        ENDDO
C
C       BEST-FIT GAUSSIAN POTENTIAL VALUES ACROSS UNIFORM GRID
        DO IPTS=0,NPTS
C
C         RADIUS ON UNIFORM SCALE
          R = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C         INITIALISE BEST-FIT POTENTIAL AT THIS RADIUS
          V = 0.0D0
          DO IZR=1,NZR
            V = V + FORI(IZ,IZR)*DEXP(-XORI(IZ,IZR)*R*R)
          ENDDO
          VG(IPTS) = V
C
        ENDDO
C
C       AVERAGE POTENTIAL VALUE
        YB = 0.0D0
        DO IPTS=0,NPTS
          YB = YB + VG(IPTS)
        ENDDO
        YB = YB/DFLOAT(NPTS+1)
C
C       PREPARATION FOR R-SQUARED VALLUE
        SRES = 0.0D0
        STOT = 0.0D0
        DO IPTS=0,NPTS
          SRES = SRES + (VG(IPTS)-VF(IPTS))**2
          STOT = STOT + (VG(IPTS)-YB      )**2
        ENDDO
C
C       R-SQUARED VALUE AND BEST-FIT CHECK
        RSQ = 1.0D0 - SRES/STOT
        IF(RSQ.GT.RSQBIG) THEN
          ITRYBG = ITRY
          RSQBIG = RSQ
        ENDIF
      
        WRITE(*,*) ITRY,XORI(IZ,1)*RNUC(IZ)*RNUC(IZ),RSQ
C
      ENDDO
C
C     IF NO SOLUTION CAN BE FOUND, ALERT THE USER
      IF(RSQBIG.LT.0.0D0) THEN
        WRITE(6, *) 'In NUCUEHL: best-fit search failed. IZ =',IZ
        WRITE(7, *) 'In NUCUEHL: best-fit search failed. IZ =',IZ
        STOP
      ENDIF
C
C**********************************************************************C
C     ***BEST-FIT ALPHA HAS NOW BEEN FOUND -- RESTORE THESE RESULTS.      C
C**********************************************************************C
C
C     STORE GEOMETRIC SET OF PARAMETERS IN XORI
      XI = ALPH*(1.0D0 + DFLOAT(ITRYBG)/DFLOAT(NTRY))
      DO IZR=1,NZR
        XORI(IZ,IZR) = XI
        RT = DFLOAT(IZR)/DFLOAT(NZR+1)
        XI = BETA*XI*(1.0D0+GAMA*RT*RT)
      ENDDO
C
C     TRANSFER MATCHING VALUES V(R) INTO Z VECTOR
      DO IZR=1,NZR
        Z(IZR) = Y(IZR)
      ENDDO
C
C     SET UP MATRIX EQUATIONS FOR SET OF FITTING GAUSSIANS
      DO IZR=1,NZR
C
C       RADIUS ON A UNIFORM SCALE
        R = RNUC(IZ)*(R0 + DFLOAT(IZR-1)*HR)
C
C       GAUSSIAN POTENTIALS AT THIS RADIUS
        DO JFT=1,NZR
          X(IZR,JFT) = DEXP(-XORI(IZ,JFT)*R*R)
        ENDDO
C
      ENDDO
C
C     SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
      CALL DGESV(NZR,1,X,NZR,IPIV,Z,NZR,INFO)
C
C     TRANSFER AMPLITUDES VALUES INTO FORI ARRAY
      DO IZR=1,NZR
        FORI(IZ,IZR) = Z(IZR)
      ENDDO
C
C     BEST-FIT GAUSSIAN POTENTIAL VALUES ACROSS UNIFORM GRID
      DO IPTS=0,NPTS
C
C       RADIUS ON UNIFORM SCALE
        R = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C       INITIALISE BEST-FIT POTENTIAL AT THIS RADIUS
        V = 0.0D0
        DO IZR=1,NZR
          V = V + FORI(IZ,IZR)*DEXP(-XORI(IZ,IZR)*R*R)
        ENDDO
        VG(IPTS) = V
C
      ENDDO
C
25    CONTINUE
C
C**********************************************************************C
C     POINT-WISE EVALUATION OF R*R*V(R) ON RADIAL GRID.                C
C**********************************************************************C
C
C     NUMBER OF DATA POINTS IN UNIFORMLY-SPACED AND EXPONENTIAL REGION
      NLIN = NRAD/10 - MOD(NRAD/10,10)
      NEXP = NRAD-NLIN
C
C     GENERATE RADIAL GRID (LINEAR FROM FEMTOMETERS, EXPONENTIAL IN AU)
      RORI =  0.0D0/CFM
      RMID = 10.0D0/CFM
      RMAX =  0.5D0
C
      HL = (RMID-RORI)/DFLOAT(NLIN)
      HE = DLOG(RMAX/RMID)/DFLOAT(NEXP)
C
      DO N=0,NLIN
        RAD(N) = RORI + HL*DFLOAT(N)
      ENDDO
C
      DO N=0,NEXP
        RAD(N+NLIN) = RMID*DEXP(HE*DFLOAT(N))
      ENDDO
C
C     ORIGIN VALUE IS A SPECIAL CASE
      VVAC(IZ,0) = 0.0D0
C
C     VALUE OF R*R*V(R) AT EACH OF THESE RADII
      DO N=1,NRAD
C
C       RADIUS R AND INITIALISE COUNTER FOR POTENTIAL
        R = RAD(N)
C
C       INTIALISE POTENTIAL COUNTER
        V = 0.0D0
C
C       INTEGRATE OVER CHARGE SOURCE
        DO M=0,NSRC
C
C         SET SOURCE RADIUS
          S  = HS*DFLOAT(M)
C
C         CHI FUNCTION ARGUMENTS
          XM = 2.0D0*DABS(R-S)/CMPF
          XP = 2.0D0*DABS(R+S)/CMPF
C
C         COMPONENTS OF INTEGRAND
          IF(0.5D0*(XM+XP).LT.XSPL) THEN
            CALL SPLNINT(XS,Y2S,D2S,XM,CM,NLW)
            CALL SPLNINT(XS,Y2S,D2S,XP,CP,NLW)    
          ELSE
            CALL SPLNINT(XB,Y2B,D2B,XM,CM,NUP)
            CALL SPLNINT(XB,Y2B,D2B,XP,CP,NUP)
          ENDIF
C
C         PERFORM THE MAPPING CHI(X) = SPLINE(X)*EXP(-X)
          CM = CM*DEXP(-XM)
          CP = CP*DEXP(-XP)
C
C         CONTRIBUTION TO INTEGRAND
          V = V + EXTINT11(S*RHONUC(IZ,S)*(CM-CP),M,NSRC)
C
        ENDDO
C
C       INTEGRATION WEIGHTING FACTORS
        V = 5.0D0*HS*V/2.99376D+5
C
C       OTHER FACTORS AND FINAL VALUE R*R*V(R)
        VVAC(IZ,N) =-2.0D0*CMPF*R*V/(3.0D0*CV)
C
      ENDDO
C
C     SKIP OUTPUT SECTION UNLESS PROMPTED
      IF(IWRT.EQ.0) GOTO 20
C
C**********************************************************************C
C     PRINT RESULTS IF PROMPTED                                        C
C**********************************************************************C
C
C     SOLUTION SET
35    FORMAT(1X,A,I2)
36    FORMAT(1X,A,F10.8,1X,A)
37    FORMAT(1X,A,I2,A,F10.5,A,F8.5,A)
      WRITE(6, *) ''
      WRITE(7, *) ''
      WRITE(6,35) 'Uehling gaussian basis set for IZ = ',IZ
      WRITE(7,35) 'Uehling gaussian basis set for IZ = ',IZ
      WRITE(6,36) 'Here ξ0 = 1/R_n^2 and R_n  = ',RNUC(IZ),'a0'
      WRITE(7,36) 'Here ξ0 = 1/R_n^2 and R_n  = ',RNUC(IZ),'a0'
      WRITE(6,36) 'Least-squares best fit: R^2 = ',RSQBIG
      WRITE(7,36) 'Least-squares best fit: R^2 = ',RSQBIG
      WRITE(6, *) REPEAT('-',45)
      WRITE(7, *) REPEAT('-',45)
      DO IFT=1,NFT
        FR = FUEH(IZ,IFT)
        XI = RNUC(IZ)*RNUC(IZ)*XUEH(IZ,IFT)
        WRITE(6,37) 'Gaussian ',IFT,': ',FR,'*exp(-',XI,' ξ0 r^2)'
        WRITE(7,37) 'Gaussian ',IFT,': ',FR,'*exp(-',XI,' ξ0 r^2)'
      ENDDO
      DO IZR=1,NZR
        FR = FORI(IZ,IZR)
        XI = RNUC(IZ)*RNUC(IZ)*XORI(IZ,IZR)
        WRITE(6,37) 'Gaussian ',IZR,': ',FR,'*exp(-',XI,' ξ0 r^2)'
        WRITE(7,37) 'Gaussian ',IZR,': ',FR,'*exp(-',XI,' ξ0 r^2)'
      ENDDO
C
      WRITE(6, *) REPEAT('-',45)
      WRITE(7, *) REPEAT('-',45)
C
C**********************************************************************C
C     PLOTTING SECTION                                                 C
C**********************************************************************C
C
C     DETAILS COMMON TO ALL PLOTS
      XAXIS  = 'r/RNUC'
      KEY(1) = 'Exact Uehling'
      KEY(2) = 'Best fit'
C
C     UNIFORM GRID VALUES R*R*V(R)
      XOUT   = 'Uehling-uniform'
      TITLE  = 'Uehling r^{2}*V(r) on uniform grid'
      YAXIS  = 'r^{2}*V(r)/Z'
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
      DO IPTS=0,NPTS
        R = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
        WRITE(8, *) R/RNUC(IZ),VF(IPTS),VG(IPTS)
      ENDDO
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C     NUMERICAL EVALUATION FOR V(R)
      XOUT   = 'Uehling-pointwise'
      TITLE  = 'Uehling V(r) on piecewise grid'
      YAXIS  = 'V(r)/Z'
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
      DO N=0,NRAD
        R = RAD(N)
        IF(R/RNUC(IZ).GT.5.0D0) GOTO 123
C       TRUE UEHLING POTENTIAL
        IF(N.EQ.0) THEN
          VTR = V0
        ELSE
          VTR = VVAC(IZ,N)/(R*R)
        ENDIF
C       APPROXIMATE VALUES
        V = 0.0D0
        DO IFT=1,NFT
          V = V + FUEH(IZ,IFT)*DEXP(-XUEH(IZ,IFT)*R*R)
        ENDDO
        DO IZR=1,NZR
          V = V + FORI(IZ,IZR)*DEXP(-XORI(IZ,IZR)*R*R)
        ENDDO
        WRITE(8, *) R/RNUC(IZ),VTR,V
      ENDDO
123   CONTINUE
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C     UEHLING DIFFERENCES
      XOUT   = 'Uehling-difference'
      TITLE  = 'Difference between V(r) on piecewise grid'
      YAXIS  = 'Delta V(r)/Z'
      KEY(1) = 'Difference'
      KEY(2) = ' '
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
      DO N=0,NRAD
        R = RAD(N)
        IF(R/RNUC(IZ).GT.5.0D0) GOTO 124
C       TRUE UEHLING POTENTIAL
        IF(N.EQ.0) THEN
          VTR = V0
        ELSE
          VTR = VVAC(IZ,N)/(R*R)
        ENDIF
C       APPROXIMATE VALUES
        V = 0.0D0
        DO IFT=1,NFT
          V = V + FUEH(IZ,IFT)*DEXP(-XUEH(IZ,IFT)*R*R)
        ENDDO
        DO IZR=1,NZR
          V = V + FORI(IZ,IFT)*DEXP(-XORI(IZ,IFT)*R*R)
        ENDDO
        WRITE(8, *) R/RNUC(IZ),VTR-V,0.0D0
c       WRITE(*, *) R/RNUC(IZ),VTR-V
      ENDDO
124   CONTINUE
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C     OVERWRITE UEHLING POTENTIAL ON THE GRID
C
C     VALUE OF R*R*V(R) AT EACH OF THESE RADII
      DO N=0,NRAD
C
C       RADIUS R
        R = RAD(N)
C
C       BEST-FIT UEHLING POTENTIAL
        VVAC(IZ,N) = 0.0D0
        DO IFT=1,NFT
          XI = XUEH(IZ,IFT)
          VVAC(IZ,N) = VVAC(IZ,N) + R*R*FUEH(IZ,IFT)*DEXP(-XI*R*R)
        ENDDO
        DO IZR=1,NZR
          XI = XORI(IZ,IFT)
          VVAC(IZ,N) = VVAC(IZ,N) + R*R*FORI(IZ,IFT)*DEXP(-XI*R*R)
        ENDDO
C      
      ENDDO
C
20    CONTINUE
C
      RETURN
      END

