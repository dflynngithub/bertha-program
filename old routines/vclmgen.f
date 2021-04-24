      SUBROUTINE VCLMGEN(IZ,NFT,IWRT,RSQBIG,RMAX,ASTART,BSTART,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  VV    VV  CCCCCC  LL       MM       MM  GGGGGG  EEEEEEEE NN    NN   C
C  VV    VV CC    CC LL       MMM     MMM GG    GG EE       NNN   NN   C
C  VV    VV CC       LL       MMMM   MMMM GG       EE       NNNN  NN   C
C  VV    VV CC       LL       MM MM MM MM GG       EEEEEE   NN NN NN   C
C   VV  VV  CC       LL       MM  MMM  MM GG   GGG EE       NN  NNNN   C
C    VVVV   CC    CC LL       MM   M   MM GG    GG EE       NN   NNN   C
C     VV     CCCCCC  LLLLLLLL MM       MM  GGGGGG  EEEEEEEE NN    NN   C
C                                                                      C
C -------------------------------------------------------------------- C
C  VCLMGEN EVALUATES THE DETAILS OF NFT GAUSSIAN NUCLEAR FUNCTIONS     C
C  FOR EACH CENTRE IN THE MOLECULE, SO AS TO BEST FIT THE FERMI MODEL. C
C  EXPONENTS XNUC ARE SELECTED MANUALLY, AND FRACTIONS FNUC CALCULATED.C
C -------------------------------------------------------------------- C
C  MATCHING CRITERIA AVAILABLE FOR GAUSSIAN AMPLITUDES:                C
C   ▶ TOTAL NORMALISED CHARGE OF UNITY <R^0>.                          C
C   ▶ SECOND MOMENT <R^2> SAME AS EMPIRICAL FORMULA OR AVAILABLE DATA. C
C   ▶ ANY OTHER INTEGER MOMENT CAN BE MATCHED WITH FERMI MODEL.        C
C   ▶ ANY OTHER REAL VALUED MOMENT (K > -3) CAN BE MATCHED WITH DATA.  C
C   ▶ FORCED AGREEMENT WITH FERMI POTENTIAL VFERMI(R) AT ANY RADIUS R. C
C**********************************************************************C
      INCLUDE 'parameters.h'
      PARAMETER(NPTS=1000,NTRY=1000)
C
      CHARACTER*5  NMDL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(4)
C
      DIMENSION CFMI(MCT)
      DIMENSION RADS(0:NPTS),POTENTIALS(4,0:NPTS),DENSITIES(4,0:NPTS)
      DIMENSION PAC(11),RAC(11),PRC(11),PLG(11)
      DIMENSION RF(-2:8),RG(-2:8),R0(-2:8),RU(-2:8)
      DIMENSION VU(0:NPTS),VF(0:NPTS),VG(0:NPTS),V0(0:NPTS)
C
      DIMENSION X(NFT,NFT),Z(NFT),Y(NFT)
      DIMENSION IPIV(NFT)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,MFT),XNUC(MCT,MFT),NNUC(MCT),NMDL(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     NUMBER OF FITTING FUNCTIONS
      NNUC(IZ) = NFT
C
C     NUCLEAR RADIUS
      IF(INT(ZNUC(IZ)).EQ.1) THEN
C       PROTON HAS A VERY SMALL RADIUS (IN FM)
        RFM = PRAD
      ELSEIF(ZNUC(IZ).GT.1.0D0.AND.ZNUC(IZ).LT.190.0D0) THEN
C       EMPIRICAL FORMULA FOR NUCLEAR RMS RADIUS IF Z < 90 (IN FM)
        RFM = 0.836D0*(ANUC(IZ)**(1.0D0/3.0D0)) + 0.57D0
C       RFM = 1.20D0*(ANUC(IZ)**(1.0D0/3.0D0))
      ELSE
C       EMPIRICAL FORMULA FOR NUCLEAR RMS RADIUS IF Z > 90 (IN FM)
        RFM = 0.770D0*(ANUC(IZ)**(1.0D0/3.0D0)) + 0.98D0
      ENDIF
C
C     CONVERT RESULT TO ATOMIC UNITS
      RNUC(IZ) = RFM/CFM
C
C     POINT-NUCLEUS MODEL
      IF(NMDL(IZ).EQ.'POINT') THEN
        RNUC(IZ) = 0.0D0
        RSQBIG   = 1.0D0
        GOTO 20
      ENDIF
C
C     GAUSSIAN MODEL
      IF(NMDL(IZ).EQ.'GAUSS') THEN
        FNUC(IZ,1) = 1.0D0
        XNUC(IZ,1) = 1.5D0/(RNUC(IZ)*RNUC(IZ))
        RSQBIG     = 1.0D0
        GOTO 20
      ENDIF
C
C**********************************************************************C
C     UNIFORM NUCLEAR MODEL DETAILS.                                   C
C**********************************************************************C
C
C     HAVEN'T DONE THIS YET. SKIP AHEAD TO FERMI.
C
C**********************************************************************C
C     FERMI NUCLEAR MODEL DETAILS.                                     C
C**********************************************************************C
C
C     FERMI SKIN THICKNESS AND HALF-DENSITY PARAMETERS (IN A0)
      TFMI(IZ) = 2.30D0/CFM
      AFMI(IZ) = 0.25D0*TFMI(IZ)/THLG
      CFMI(IZ) = 5.0D0*RNUC(IZ)*RNUC(IZ)/3.0D0
     &         - 7.0D0*PI*PI*AFMI(IZ)*AFMI(IZ)/3.0D0
      CFMI(IZ) = DSQRT(CFMI(IZ))
      AFMI(IZ) = 1.03900D-05
      CFMI(IZ) = 1.39190D-04
C
C     CHECK WHETHER FERMI MODEL IS POSSIBLE
      EPS = DSQRT(1.4D0)*PI*0.25D0*TFMI(IZ)/THLG
      IF(EPS.GT.RNUC(IZ)) THEN
        NNUC(IZ)   = 1
        FNUC(IZ,1) = 1.0D0
        XNUC(IZ,1) = 1.5D0/(RNUC(IZ)*RNUC(IZ))
        RSQBIG     = 1.0D0
        NMDL(IZ)   = 'GAUSS'
        GOTO 20
      ENDIF
C
C     THESE PARAMETERS MAKE FOR NEATER EXPRESSIONS
      PST = 1.0D0
      RST = 1.0D0
      DO I=1,11
        PST = PST*PI*AFMI(IZ)/CFMI(IZ)
        RST = RST*AFMI(IZ)/CFMI(IZ)
        PAC(I) = PST
        RAC(I) = RST
        PLG(I) = POLYLOG(I,-1.0D0/RAC(1))
      ENDDO
C
C     FERMI MODEL NORMALISATION CONSTANT RHO0
      RHO0 = 1.0D0 + PAC(2) - 6.0D0*RAC(3)*PLG(3)
      RHO0 = 3.0D0/(RHO0*4.0D0*PI*(CFMI(IZ)**3))
C
C     FERMI NORMALISATION CHARGE DENSITY
      DO I=1,11
        PRC(I) = 4.0D0*PI*RHO0*(CFMI(IZ)**(I))/DFLOAT(I)
      ENDDO
C
C     FERMI MODEL RADIAL MOMENTS
      RF(-2) = PRC(1)*(1.0D0 - RAC(1)*PLG(1))
      RF(-1) = PRC(2)*(1.0D0 + PAC(2)/3.0D0 + 2.0D0*RAC(2)*PLG(2))
      RF( 0) = PRC(3)*(1.0D0 + PAC(2) - 6.0D0*RAC(3)*PLG(3))
      RF( 1) = PRC(4)*(1.0D0 + 2.0D0*PAC(2) + 7.0D0*PAC(4)/15.0D0
     &       + 24.0D0*RAC(4)*PLG(4))
      RF( 2) = PRC(5)*(1.0D0 + 10.0D0*PAC(2)/3.0D0
     &       + 7.0D0*PAC(4)/3.0D0 - 120.0D0*RAC(5)*PLG(5))
      RF( 3) = PRC(6)*(1.0D0 + 5.0D0*PAC(2) + 7.0D0*PAC(4)
     &       + 31.0D0*PAC(6)/21.0D0 + 720.0D0*RAC(6)*PLG(6))
      RF( 4) = PRC(7)*(1.0D0 + 7.0D0*PAC(2) + 49.0D0*PAC(4)/3.0D0
     &       + 31.0D0*PAC(6)/3.0D0 - 5040.0D0*RAC(7)*PLG(7))
      RF( 5) = PRC(8)*(1.0D0 + 28.0D0*PAC(2) + 98.0D0*PAC(4)/3.0D0
     &       + 124.0D0*PAC(6)/3.0D0 + 127.0D0*PAC(8)/15.0D0 
     &       + 40320.0D0*RAC(8)*PLG(8))
      RF( 6) = PRC(9)*(1.0D0 + 12.0D0*PAC(2) + 294.0D0*PAC(4)/5.0D0
     &       + 124.0D0*PAC(6)/3.0D0 + 381.0D0*PAC(8)
     &       - 362880.0D0*RAC(9)*PLG(9))
      RF( 7) = PRC(10)*(1.0D0 + 150.0D0*PAC(2) + 98.0D0*PAC(4)
     &       + 310.0D0*PAC(6) + 381.0D0*PAC(8)  + 2555.0D0*PAC(10)
     &       + 3628800.0D0*RAC(10)*PLG(10))
      RF( 8) = PRC(11)*(1.0D0 + 55.0D0*PAC(2) + 154.0D0*PAC(4)
     &       + 682.0D0*PAC(6) + 1397.0D0*PAC(8) + 2555.0D0*PAC(10)
     &       - 39916800.0D0*RAC(11)*PLG(11))
C
C     SAVE THE FERMI POTENTIAL VALUES
      DO IPTS=0,NPTS
        R = 5.0D0*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
        VF(IPTS) = VNFERMI(R,CFMI(IZ),TFMI(IZ))
      ENDDO
C
C     SET UP MATRIX EQUATIONS FOR FERMI DISTRIBUTION (TO BE MATCHED)
      Y(1) = 1.0D0
      Y(2) = RF(2)
      ITG = 3
      DO IFT=3,NNUC(IZ)
C       R = 2.0D0*DFLOAT(IFT-ITG)*RNUC(IZ)/DFLOAT(NNUC(IZ)-ITG)
        R = RMAX*DFLOAT(IFT-ITG)*RNUC(IZ)/DFLOAT(NNUC(IZ)-ITG)
        Y(IFT) = VNFERMI(R,CFMI(IZ),TFMI(IZ))
      ENDDO
C
C**********************************************************************C
C     GIVEN GEOMETRIC BETA AND NNUC(IZ)=NFT, FIND OPTIMAL ALPHA VALUE. C
C**********************************************************************C
C
      RSQBIG = 0.0D0
      ITRYBG = 0

C      ASTART = 0.10D0
C      BSTART = 0.25D0
C      BETA   = 1.50D0
C
      DO ITRY=0,NTRY
C
C       SET OF GAUSSIAN EXPONENTS (SEARCH WITHIN A BETA WINDOW)
        XI0 = ASTART + BSTART*DFLOAT(ITRY)/DFLOAT(NTRY)
        XI  = XI0/(RNUC(IZ)**2)
        DO IFT=1,NNUC(IZ)
          XNUC(IZ,IFT) = XI
          XI = BETA*XI
        ENDDO
C
C       RESET FERMI Z MATRIX
        DO IFT=1,NNUC(IZ)
          Z(IFT) = Y(IFT)
        ENDDO
C
C       SET UP MATRIX EQUATIONS FOR SET OF FITTING GAUSSIANS
C
C       FIRST EQUATION ASSERTS NORMALISATION
        DO JFT=1,NNUC(IZ)
          X(1,JFT) = 1.0D0
        ENDDO
C
C       SECOND EQUATION ASSERTS R^2 VALUE
        DO JFT=1,NNUC(IZ)
          X(2,JFT) = 1.5D0/XNUC(IZ,JFT)
        ENDDO
C
C       ALL OTHER EQUATIONS MATCH TO THE FERMI NUCLEUS
        DO IFT=3,NNUC(IZ)
C
C         MATCHING RADIUS: LINEAR GRID BETWEEN 0.0D0 AND 2.0D0*RNUC(IZ)
          R = 2.0D0*DFLOAT(IFT-ITG)*RNUC(IZ)/DFLOAT(NNUC(IZ)-ITG)
C
C         ERF(R)/R POTENTIALS AT THIS RADIUS
          DO JFT=1,NNUC(IZ)
            X(IFT,JFT) = VNGAUSS(R,XNUC(IZ,JFT),1.0D0)
          ENDDO
C
        ENDDO
C
C       SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
        CALL DGESV(NNUC(IZ),1,X,NNUC(IZ),IPIV,Z,NNUC(IZ),INFO)
C
C       TRANSFER THE X VALUES TO FRACTIONAL ARRAY
        DO IFT=1,NNUC(IZ)
          FNUC(IZ,IFT) = Z(IFT)
        ENDDO
C
C       BEST-FIT GAUSSIAN POTENTIAL
        DO IPTS=0,NPTS
          R = 5.0D0*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
          VG(IPTS) = 0.0D0
          DO IFT=1,NNUC(IZ)
            VG(IPTS) = VG(IPTS) + VNGAUSS(R,XNUC(IZ,IFT),FNUC(IZ,IFT))
          ENDDO
        ENDDO
C
C       AVERAGE EVALUATED POTENTIAL OVER ALL RADII
        YB = 0.0D0
        DO IPTS=0,NPTS
          YB = YB + VG(IPTS)
        ENDDO
        YB = YB/DFLOAT(NPTS+1)
C
C       BITS AND PIECES FOR THE R-SQUARED VALUE
        SRES = 0.0D0
        STOT = 0.0D0
        DO IPTS=0,NPTS
          SRES = SRES + (VG(IPTS)-VF(IPTS))**2
          STOT = STOT + (VG(IPTS)-YB      )**2
        ENDDO
C
C       DECIDE WHETHER THIS IS THE BEST FIT SO FAR
        RSQ = 1.0D0 - SRES/STOT
        IF(RSQ.GT.RSQBIG) THEN
          ITRYBG = ITRY
          RSQBIG = RSQ
        ENDIF
C
      ENDDO
C
C**********************************************************************C
C     END OF BEST-FIT DETERMINATION -- NOW IMPLEMENT THE BEST ALPHA.   C
C**********************************************************************C
C
C     SET OF GAUSSIAN EXPONENTS
      XI0 = ASTART + BSTART*DFLOAT(ITRYBG)/DFLOAT(NTRY)
      XI  = XI0/(RNUC(IZ)**2)
      DO IFT=1,NNUC(IZ)
        XNUC(IZ,IFT) = XI
        XI = BETA*XI
      ENDDO
C
C     RESET FERMI Z MATRIX
      DO IFT=1,NNUC(IZ)
        Z(IFT) = Y(IFT)
      ENDDO
C
C     SET UP MATRIX EQUATIONS FOR SET OF FITTING GAUSSIANS
C
C     FIRST EQUATION ASSERTS NORMALISATION
      DO JFT=1,NNUC(IZ)
        X(1,JFT) = 1.0D0
      ENDDO
C
C     SECOND EQUATION ASSERTS R^2 VALUE
      DO JFT=1,NNUC(IZ)
        X(2,JFT) = 1.5D0/XNUC(IZ,JFT)
      ENDDO
C
C     ALL OTHER EQUATIONS MATCH TO THE FERMI NUCLEUS
      DO IFT=3,NNUC(IZ)
C
C       MATCHING RADIUS: LINEAR GRID BETWEEN 0.0D0 AND 2.0D0*RNUC(IZ)
        R = 2.0D0*DFLOAT(IFT-ITG)*RNUC(IZ)/DFLOAT(NNUC(IZ)-ITG)
C
C       ERF(R)/R POTENTIALS AT THIS RADIUS
        DO JFT=1,NNUC(IZ)
          X(IFT,JFT) = VNGAUSS(R,XNUC(IZ,JFT),1.0D0)
        ENDDO
C
      ENDDO
C
C     SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
      CALL DGESV(NNUC(IZ),1,X,NNUC(IZ),IPIV,Z,NNUC(IZ),INFO)
C
C     TRANSFER THE X VALUES TO FRACTIONAL ARRAY
      DO IFT=1,NNUC(IZ)
        FNUC(IZ,IFT) = Z(IFT)
      ENDDO
C
C     QUIT HERE IF NO WRITTEN RESULTS ARE REQUIRED
      IF(IWRT.EQ.0) GOTO 20
C
C     WRITE THE BEST-FIT SOLUTION
35    FORMAT(1X,A,I2)
36    FORMAT(1X,A,F10.8,1X,A)
37    FORMAT(1X,A,I2,A,F10.5,A,F8.5,A)
      WRITE(6, *) ''
      WRITE(7, *) ''
      WRITE(6,35) 'Nuclear gaussian basis set for IZ = ',IZ
      WRITE(7,35) 'Nuclear gaussian basis set for IZ = ',IZ
      WRITE(6,36) 'Here ξ0 = 3/2R_n^2 and R_n  = ',RNUC(IZ),'a0'
      WRITE(7,36) 'Here ξ0 = 3/2R_n^2 and R_n  = ',RNUC(IZ),'a0'
      WRITE(6,36) 'Least-squares best fit: R^2 = ',RSQBIG
      WRITE(7,36) 'Least-squares best fit: R^2 = ',RSQBIG
      WRITE(6, *) REPEAT('-',45)
      WRITE(7, *) REPEAT('-',45)
      DO IFT=1,NNUC(IZ)
        XMULT = 2.0D0*RNUC(IZ)*RNUC(IZ)*XNUC(IZ,IFT)/3.0D0
        WRITE(6,37) 'Gaussian ',IFT,': ',FNUC(IZ,IFT),'*exp(-',XMULT,
     &                                                       ' ξ0 r^2)'
        WRITE(7,37) 'Gaussian ',IFT,': ',FNUC(IZ,IFT),'*exp(-',XMULT,
     &                                                       ' ξ0 r^2)'
      ENDDO
      WRITE(6, *) REPEAT('-',45)
      WRITE(7, *) REPEAT('-',45)
C
C**********************************************************************C
C     RADIAL MOMENTS FOR THIS AND OTHER MODELS.                        C
C**********************************************************************C
C
C     SINGLE GAUSSIAN NUCLEAR CHARGE
      DO I=-2,8
        FC1 = 4.0D0/PI12
        FC2 = (2.0D0/3.0D0)**(0.5D0*I)
        FC3 = RNUC(IZ)**I
        FC4 = GAMHLF(5+I)
        FC5 = 3.0D0+I
        R0(I) = FC1*FC2*FC3*FC4/FC5
      ENDDO
C
C     UNIFORMLY-CHARGED NUCLEUS
      DO I=-2,8
        FC1 = 3.0D0
        FC2 = (5.0D0/3.0D0)**(0.5D0*I)
        FC3 = RNUC(IZ)**I
        FC5 = 3.0D0+I
        RU(I) = FC1*FC2*FC3/FC5
      ENDDO
C
C     FINITE SET OF GAUSSIANS
      DO I=-2,8
        RG(I) = 0.0D0
        DO IFT=1,NNUC(IZ)
          XI  = XNUC(IZ,IFT)
          FC  = FNUC(IZ,IFT)
          FC1 = FC*2.0D0/PI12
          FC2 = GAMHLF(3+I)
          FC3 = XI**(0.5D0*I)
          RG(I) = RG(I) + FC1*FC2/FC3
        ENDDO
      ENDDO
C
C     CREATE TABLE FOR THE MOMENTS
49    FORMAT(2X,A,7X,A,12X,A,12X,A,12X,A)
52    FORMAT(2X,'<r^',I2,'>^1/',I2,4(2X,ES16.8))
53    FORMAT(2X,'<r^',I2,'>',5X,4(2X,ES16.8))
      WRITE(6, *) ''
      WRITE(7, *) ''
      WRITE(6,49) 'Moment (fm or e)','Gauss0','Unifrm',' Fermi','Bestft'
      WRITE(7,49) 'Moment (fm or e)','Gauss0','Unifrm',' Fermi','Bestft'
      WRITE(6, *) REPEAT('-',84)
      WRITE(7, *) REPEAT('-',84)
      CF = 1.0D0/(CFM**2)
      DO I=-2,8
        IF(I.NE.0) THEN
          RU(I) = (RU(I)*CF)**(1.0D0/DFLOAT(I))
          RF(I) = (RF(I)*CF)**(1.0D0/DFLOAT(I))
          RG(I) = (RG(I)*CF)**(1.0D0/DFLOAT(I))
          R0(I) = (R0(I)*CF)**(1.0D0/DFLOAT(I))
          WRITE(6,52) I,I,R0(I),RU(I),RF(I),RG(I)
          WRITE(7,52) I,I,R0(I),RU(I),RF(I),RG(I)
        ELSE
          RU(I) = RU(I)*CF
          RF(I) = RF(I)*CF
          RG(I) = RG(I)*CF
          R0(I) = R0(I)*CF
          WRITE(6,53) I,  R0(I),RU(I),RF(I),RG(I)
          WRITE(7,53) I,  R0(I),RU(I),RF(I),RG(I)
        ENDIF
        CF = CF*CFM
      ENDDO
      WRITE(6, *) REPEAT('-',84)
      WRITE(7, *) REPEAT('-',84)
      WRITE(6, *) ''
      WRITE(7, *) ''
C
C**********************************************************************C
C     SCALAR POTENTIALS FOR THIS AND OTHER MODELS.                     C
C**********************************************************************C
C
C     POTENTIALS FOR ALL EXACT NUCLEAR MODELS
      DO IPTS=0,NPTS
C
C       RADIAL VALUE
        R = 5.0D0*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C       FERMI POTENTIAL
        VF(IPTS) = VNFERMI(R,CFMI(IZ),TFMI(IZ))
C
C       SINGLE GAUSSIAN POTENTIAL
        V0(IPTS) = VNGAUSS(R,1.5D0/(RNUC(IZ)*RNUC(IZ)),1.0D0)
C
C       UNIFORM POTENTIAL
        IF(R.LT.DSQRT(5.0D0/3.0D0)*RNUC(IZ)) THEN
          CU = DSQRT(5.0D0/3.0D0)*RNUC(IZ)
          VU(IPTS) = 0.5D0*(3.0D0 - (R/CU)**2)/CU
        ELSE
          VU(IPTS) = 1.0D0/R
        ENDIF
C
C       FINITE SET OF GAUSSIANS
        VG(IPTS) = 0.0D0
        DO IFT=1,NNUC(IZ)
          VG(IPTS) = VG(IPTS) + VNGAUSS(R,XNUC(IZ,IFT),FNUC(IZ,IFT))
        ENDDO
C
C       WRITE ALL THIS INTO A BIG ARRAY
        RADS(IPTS) = R/RNUC(IZ)
        POTENTIALS(1,IPTS) = VU(IPTS)
        POTENTIALS(2,IPTS) = VF(IPTS)
        POTENTIALS(3,IPTS) = V0(IPTS)
        POTENTIALS(4,IPTS) = VG(IPTS)
C
      ENDDO
C
C     CREATE A TABLE FOR THE POTENTIALS
      GOTO 23
50    FORMAT(1X,F8.4,'R',3X,F20.10,2X,F20.10,2X,F20.10)
51    FORMAT(2X,A,13X,A,11X,A,11X,A)
      WRITE(6,51) 'Radius r','V_Fermi(r)','V_GaussN(r)','V_Gauss0(r)'
      WRITE(7,51) 'Radius r','V_Fermi(r)','V_GaussN(r)','V_Gauss0(r)'
      WRITE(6, *) REPEAT('-',76)
      WRITE(7, *) REPEAT('-',76)
      DO IPTS=0,NPTS
        R = 5.0D0*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
        WRITE(6,50) R/RNUC(IZ),VF(IPTS),VG(IPTS),V0(IPTS)
        WRITE(7,50) R/RNUC(IZ),VF(IPTS),VG(IPTS),V0(IPTS)
      ENDDO
      WRITE(6, *) REPEAT('-',76)
      WRITE(7, *) REPEAT('-',76)
23    CONTINUE
C
C     FILE NAME AND TITLES
      XOUT   = 'Potentials'
      TITLE  = 'Potentials'
      XAXIS  = 'r/RNUC'
      YAXIS  = '-V(r)/Z'
      KEY(1) = 'Uniform'
      KEY(2) = 'Fermi'
      KEY(3) = 'Gauss(1)'
      KEY(4) = 'Gauss(N)'
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C       BEGIN LOOP OVER ALL POINTS
        DO IPTS=0,NPTS
C
C         PRINT DISTANCE ALONG PATH AND FOUR-CURRENT
          WRITE(8, *) RADS(IPTS),(ABS(POTENTIALS(ITP,IPTS)),ITP=1,4)
C
        ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,4,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C**********************************************************************C
C     RADIAL CHARGE DENSITIES FOR THIS AND OTHER MODELS.               C
C**********************************************************************C
C
C     CHARGE DENSITIES FOR ALL EXACT NUCLEAR MODELS
      DO IPTS=0,NPTS
C
C       RADIAL VALUE
        R = 5.0D0*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C       CHARGE DENSITY AT THIS RADIUS
        RADS(IPTS) = R/RNUC(IZ)
        DENSITIES(1,IPTS) = RHONUC('UNIFM',IZ,R)
        DENSITIES(2,IPTS) = RHONUC('FERMI',IZ,R)
        DENSITIES(3,IPTS) = RHONUC('GAUSS',IZ,R)
        DENSITIES(4,IPTS) = 0.0D0
        DO IFT=1,NNUC(IZ)
          RHOI = FNUC(IZ,IFT)*XNUC(IZ,IFT)*DSQRT(XNUC(IZ,IFT))/PI32
          RHOI = RHOI*DEXP(-XNUC(IZ,IFT)*R*R)
          DENSITIES(4,IPTS) = DENSITIES(4,IPTS) + RHOI
        ENDDO
C
      ENDDO
C
C     FILE NAME AND TITLES
      XOUT   = 'Densities'
      TITLE  = 'Densities'
      XAXIS  = 'r/RNUC'
      YAXIS  = '-rho(r)/Z'
      KEY(1) = 'Uniform'
      KEY(2) = 'Fermi'
      KEY(3) = 'Gauss(1)'
      KEY(4) = 'Gauss(N)'
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C       BEGIN LOOP OVER ALL POINTS
        DO IPTS=0,NPTS
C
C         PRINT DISTANCE ALONG PATH AND FOUR-CURRENT
          WRITE(8, *) RADS(IPTS),(ABS(DENSITIES(ITP,IPTS)),ITP=1,4)
C
        ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,4,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
20    CONTINUE
C
      RETURN
      END
