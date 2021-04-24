      SUBROUTINE NUCUEHL(IZ,NFT,IWRT)
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
      PARAMETER(NPTS=6000)
C
      CHARACTER*3  ATRM
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(2)
C
      DIMENSION XS(0:NLW),Y2S(0:NLW),D2S(0:NLW),Y1S(0:NLW),D1S(0:NLW),
     &          XB(0:NUP),Y2B(0:NUP),D2B(0:NUP)
      DIMENSION VF(0:NPTS),VG(0:NPTS),RU(0:NPTS)
      DIMENSION X(NFT,NFT),Z(NFT),Y(NFT),RM(NFT),IPIV(NFT)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),CFMI(MCT),
     &            RNUC(MCT),FNUC(MCT,MFT),XNUC(MCT,MFT),NNUC(MCT)
      COMMON/BQED/RUEH(MCT),FUEH(MCT,MFT),XUEH(MCT,MFT),NUEH(MCT)
      COMMON/CONV/CHZ,CFM,CDB
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/UEHL/RAD(0:NRAD),VVAC(MCT,0:NRAD),RORI,RMID,RMAX,NLIN,NEXP
C
      WRITE(6, *) 'Generating nuclear Uehling potentials...'
      WRITE(7, *) 'Generating nuclear Uehling potentials...'
C
C     MULTIPLIER FOR COMPTON WAVELENGTH (USE FOR MUON OR TAUON FIELD)
      CMPF = CMPW/1.0D0
C
C     CHECK THAT THERE ARE ENOUGH FITTING FUNCTIONS
      IF(NFT.LT.20) THEN
        WRITE(6, *) 'In NUCUEHL: need more fitting functions. NFT =',NFT
        WRITE(7, *) 'In NUCUEHL: need more fitting functions. NFT =',NFT
        STOP
      ENDIF
C
C     NUMBER OF FITTING FUNCTIONS
      NUEH(IZ) = NFT
C
C     RMS MOMENT OF UEHLING CHARGE DENSITY (MIGHT COME IN HANDY)
      RUEH(IZ) = 8.0D0*CMPF*CMPF/(5.0D0*CV)
C
C**********************************************************************C
C     GENERATE AND IMPORT THE INTERPOLATED CHI FUNCTIONS               C
C**********************************************************************C
C
C     GENERATE UEHLING ANCILLIARY FUNCTIONS AND WRITE TO FILE
C     CALL CHIGEN
C
C     IMPORT SPLINE DATA FOR X*CHI_1(X) GIVEN 0 <= X <= XSPL
      OPEN(UNIT=50,FILE='chi1x.dat',STATUS='UNKNOWN')
      REWIND(UNIT=50)
      READ(50,*) XSPL
      DO N=0,NLW
        READ(50,*) XS(N),Y1S(N),D1S(N)
      ENDDO
      CLOSE(UNIT=50)
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
      PMAX = SMAX*RHONUC(IZ,SMAX)
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
C       X*CHI_1(X) VALUE (ARGUMENT WILL ALWAYS BE LESS THAN 0.6D0)
        X0 = 2.0D0*S/CMPF
        CALL SPLNINT(XS,Y1S,D1S,X0,C0,NLW)
C
C       CONTRIBUTION TO INTEGRAND
        V0 = V0 + EXTINT11(RHONUC(IZ,S)*C0,M,NSRC)
C
      ENDDO
C
C     INTEGRATION WEIGHTING FACTORS
      V0 = 5.0D0*HS*V0/2.99376D+5
      V0 =-4.0D0*CMPF*V0/(3.0D0*CV)
C
C**********************************************************************C
C     UEHLING POTENTIAL VALUES ON THE TESTING GRID.                    C
C**********************************************************************C
C
C     UEHLING POTENTIAL SAMPLING VALUES (FOR R-SQUARED CALCULATION)
      RLIM = 250.0D0
C
C     STORE R*R*V(R) ON THE UNIFORM GRID IN VF(0:NPTS)
      RU(0) = 0.0D0
      VF(0) = 0.0D0
C
      DO IPTS=1,NPTS
C
C       RADIUS ON UNIFORM SCALE
        R = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C       STORE FOR LATER
        RU(IPTS) = R
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
C       OTHER FACTORS NEEDED FOR V(R)
        VF(IPTS) =-2.0D0*CMPF*VF(IPTS)/(3.0D0*CV*R)
C
      ENDDO
C
C**********************************************************************C
C     UEHLING POTENTIAL VALUES ON THE MATCHING GRID.                   C
C**********************************************************************C
C
C     DIVIDE NFT POINTS INTO DIRECT AND WEIGHTED REGIONS
      NDCT = NFT/6
      NWGT = NFT-NDCT
C
C     UEHLING POTENTIAL MATCHING VALUES
      R0 = 0.0D0
      RC = 3.0D0
      RN = 200.0D0
      HD = (RC-R0)/DFLOAT(NDCT)
      HW = DLOG(RN/RC)/DFLOAT(NWGT-1)
C
C     STORE R*R*V(R) MATCHING VALUES IN Y(NFT) -- MATRIX EQN SOLUTIONS
      RM(1) = 0.0D0
      Y(1)  = V0
C
      DO IFT=2,NFT
C
C       RADIUS DEPENDS ON WHICH REGION WE'RE CONSIDERING
        IF(IFT.LE.NDCT) THEN
          R = RNUC(IZ)*(R0 + DFLOAT(IFT-1)*HD)
        ELSE
          R = RNUC(IZ)*RC*DEXP(DFLOAT(IFT-NDCT-1)*HW)
        ENDIF
C
C       STORE RADIUS VALUE FOR LATER
        RM(IFT) = R
C
C       INITIALISE VALUE FOR THE POTENTIAL HERE
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
C       OTHER FACTORS AND UEHLING POTENTIAL V(R)
        Y(IFT) =-2.0D0*CMPF*Y(IFT)/(3.0D0*CV*R)
C
C       DATA POINTS IN THE WEIGHTED REGION MUST BE MULTIPLIED BY R*R
        IF(IFT.GT.NDCT) THEN
          Y(IFT) = R*R*Y(IFT)
        ENDIF
C
      ENDDO
C
C**********************************************************************C
C     GIVEN MODIFIED GEOMETRIC SET OF NFT PARAMETERS, FIND BEST ALPHA. C
C**********************************************************************C
C
C     GEOMETRIC PARAMETER SEARCH DETAILS
      ALPH = 2.00D-5
      BETA = 1.60D+0
      GAMA = 0.21D+0
C
C     NUMBER OF INCREMENTS IN SEARCHES
      NAF = 1200
C
C     INITIALISE BEST-FIT R-SQUARED AND INDEX FOR ALPHA
      RSBG = 0.0D0
      IABG = 0
C
C     ITERATE OVER ALPHA SEARCH
      DO IA=0,NAF
C
C       TRIAL ALPHA PARAMETER AND FROZEN BETA PARAMETER
        AF = ALPH*(1.0D0 + 3.0D0*DFLOAT(IA)/DFLOAT(NAF))/(RNUC(IZ)**2)
        BF = BETA
C
C       STORE GEOMETRIC SET OF PARAMETERS IN XUEH
        XI = AF
        DO IFT=1,NFT
          XUEH(IZ,IFT) = XI
          RT = DFLOAT(IFT)/DFLOAT(NFT+1)
          XI = BF*XI*(1.0D0 + GAMA*RT*RT)
        ENDDO
C
C       TRANSFER MATCHING VALUES V(R) AND R*R*V(R) INTO Z VECTOR
        DO IFT=1,NFT
          Z(IFT) = Y(IFT)
        ENDDO
C
C       SET UP MATRIX EQUATIONS FOR SET OF FITTING GAUSSIANS
        DO IFT=1,NFT
C
C         MATCHING RADIUS
          R = RM(IFT)
C
C         GAUSSIAN POTENTIALS AT THIS RADIUS
          DO JFT=1,NFT
            IF(IFT.LE.NDCT) THEN
              X(IFT,JFT) =     DEXP(-XUEH(IZ,JFT)*R*R)
            ELSE
              X(IFT,JFT) = R*R*DEXP(-XUEH(IZ,JFT)*R*R)
            ENDIF
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
          R = RU(IPTS)
C
C         INITIALISE BEST-FIT POTENTIAL AT THIS RADIUS
          VG(IPTS) = 0.0D0
          DO IFT=1,NFT
            VG(IPTS) = VG(IPTS) + FUEH(IZ,IFT)*DEXP(-XUEH(IZ,IFT)*R*R)
          ENDDO
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
C       PREPARATION FOR R-SQUARED VALUE
        SRES = 0.0D0
        STOT = 0.0D0
        DO IPTS=0,NPTS
          SRES = SRES + RU(IPTS)*RU(IPTS)*(VG(IPTS)-VF(IPTS))**2
          STOT = STOT + RU(IPTS)*RU(IPTS)*(VG(IPTS)-YB      )**2
        ENDDO
C
C       R-SQUARED VALUE AND BEST-FIT CHECK
        RSQ = 1.0D0 - SRES/STOT
        IF(RSQ.GT.RSBG) THEN
          IABG = IA
          RSBG = RSQ
        ENDIF
C
      ENDDO
C
C     IF NO SOLUTION CAN BE FOUND, ALERT THE USER
      IF(RSBG.LT.0.0D0) THEN
        WRITE(6, *) 'In NUCUEHL: best-fit search failed. IZ =',IZ
        WRITE(7, *) 'In NUCUEHL: best-fit search failed. IZ =',IZ
        STOP
      ENDIF
C
C**********************************************************************C
C     BEST-FIT ALPHA HAS NOW BEEN FOUND -- RESTORE THESE RESULTS.      C
C**********************************************************************C
C
C     BEST-FIT VALUES FOR ALPHA, BETA AND GAMMA
      AF = ALPH*(1.0D0 + 3.0D0*DFLOAT(IABG)/DFLOAT(NAF))
      BF = BETA
      GF = GAMA
C
C     STORE GEOMETRIC SET OF PARAMETERS IN XUEH
      XI = AF/(RNUC(IZ)**2)
      DO IFT=1,NFT
        XUEH(IZ,IFT) = XI
        RT = DFLOAT(IFT)/DFLOAT(NFT+1)
        XI = BF*XI*(1.0D0 + GF*RT*RT)
      ENDDO
C
C     TRANSFER MATCHING VALUES V(R) AND R*R*V(R) INTO Z VECTOR
      DO IFT=1,NFT
        Z(IFT) = Y(IFT)
      ENDDO
C
C     SET UP MATRIX EQUATIONS FOR SET OF FITTING GAUSSIANS
      DO IFT=1,NFT
C
C       MATCHING RADIUS
        R = RM(IFT)
C
C       GAUSSIAN POTENTIALS AT THIS RADIUS
        DO JFT=1,NFT
          IF(IFT.LE.NDCT) THEN
            X(IFT,JFT) =     DEXP(-XUEH(IZ,JFT)*R*R)
          ELSE
            X(IFT,JFT) = R*R*DEXP(-XUEH(IZ,JFT)*R*R)
          ENDIF
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
37    FORMAT(1X,A,F10.8,A,F6.4,A,F6.4)
38    FORMAT(1X,A,I2,A,F10.5,A,F8.5,A)
      WRITE(6, *) ''
      WRITE(7, *) ''
      WRITE(6,35) 'Uehling gaussian basis set for IZ = ',IZ
      WRITE(7,35) 'Uehling gaussian basis set for IZ = ',IZ
      WRITE(6, *) REPEAT('-',45)
      WRITE(7, *) REPEAT('-',45)
      WRITE(6,36) 'Here ξ0 = 1/R_n^2 and R_n  = ',RNUC(IZ),'a0'
      WRITE(7,36) 'Here ξ0 = 1/R_n^2 and R_n  = ',RNUC(IZ),'a0'
      WRITE(6,37) 'α = ',AF,' ξ0,  β = ',BF,',  γ = ',GF
      WRITE(7,37) 'α = ',AF,' ξ0,  β = ',BF,',  γ = ',GF
      WRITE(6, *) REPEAT('-',45)
      WRITE(7, *) REPEAT('-',45)
      WRITE(6,36) 'Least-squares best fit: R^2 = ',RSBG
      WRITE(7,36) 'Least-squares best fit: R^2 = ',RSBG
      WRITE(6, *) REPEAT('-',45)
      WRITE(7, *) REPEAT('-',45)
      DO IFT=1,NFT
        FR = FUEH(IZ,IFT)
        XI = RNUC(IZ)*RNUC(IZ)*XUEH(IZ,IFT)
        WRITE(6,38) 'Gaussian ',IFT,': ',FR,'*exp(-',XI,' ξ0 r^2)'
        WRITE(7,38) 'Gaussian ',IFT,': ',FR,'*exp(-',XI,' ξ0 r^2)'
      ENDDO
C
      WRITE(6, *) REPEAT('-',45)
      WRITE(7, *) REPEAT('-',45)
      RETURN
C
C**********************************************************************C
C     BEST-FIT UEHLING POTENTIAL ON THE TESTING GRID.                  C
C**********************************************************************C
C
C     BEST-FIT GAUSSIAN POTENTIAL VALUES ACROSS UNIFORM GRID
      DO IPTS=0,NPTS
C
C       RADIUS ON UNIFORM SCALE
        R = RU(IPTS)
C
C       INITIALISE BEST-FIT POTENTIAL AT THIS RADIUS
        VG(IPTS) = 0.0D0
        DO IFT=1,NFT
          VG(IPTS) = VG(IPTS) + FUEH(IZ,IFT)*DEXP(-XUEH(IZ,IFT)*R*R)
        ENDDO
C
      ENDDO
C
C**********************************************************************C
C     EXACT UEHLING POTENTIAL (WEIGHTED BY R*R) ON COMPOSITE GRID.     C
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
C**********************************************************************C
C     PLOTTING SECTION                                                 C
C**********************************************************************C
C
C     NAME TAG FOR PDF DOCUMENTS
      IF(INT(ANUC(IZ)).LT.10) THEN
        WRITE(ATRM,'(A,I1)') '00',INT(ANUC(IZ))
      ELSEIF(INT(ANUC(IZ)).LT.100) THEN
        WRITE(ATRM,'(A,I2)') '0',INT(ANUC(IZ))
      ELSE
        WRITE(ATRM,'(I3)') INT(ANUC(IZ))
      ENDIF
C
C     DETAILS COMMON TO ALL PLOTS
      XAXIS  = 'r/RNUC'
      KEY(1) = 'Exact Uehling'
      KEY(2) = 'Best fit'
C
C     UNIFORM GRID VALUES R*R*V(R)
      XOUT   = 'Uehling-uniform'//TRIM(ATRM)
      TITLE  = 'Uehling r^{2}*V(r) on uniform grid'
      YAXIS  = 'r^{2}*V(r)/Z'
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
      DO IPTS=0,NPTS
        R = RU(IPTS)
        WRITE(8, *) R/RNUC(IZ),R*R*VF(IPTS),R*R*VG(IPTS)
      ENDDO
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C     NUMERICAL EVALUATION FOR V(R)
      XOUT   = 'Uehling-pointwise'//TRIM(ATRM)
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
        WRITE(8, *) R/RNUC(IZ),VTR,V
      ENDDO
123   CONTINUE
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
20    CONTINUE
C
      RETURN
      END

