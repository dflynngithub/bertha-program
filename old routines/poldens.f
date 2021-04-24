      SUBROUTINE POLDENS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    PPPPPPP   OOOOOO  LL       DDDDDDD  EEEEEEEE NN    NN  SSSSSS     C
C    PP    PP OO    OO LL       DD    DD EE       NNN   NN SS    SS    C
C    PP    PP OO    OO LL       DD    DD EE       NNNN  NN SS          C
C    PP    PP OO    OO LL       DD    DD EEEEEE   NN NN NN  SSSSSS     C
C    PPPPPPP  OO    OO LL       DD    DD EE       NN  NNNN       SS    C
C    PP       OO    OO LL       DD    DD EE       NN   NNN SS    SS    C
C    PP        OOOOOO  LLLLLLLL DDDDDDD  EEEEEEEE NN    NN  SSSSSS     C
C                                                                      C
C -------------------------------------------------------------------- C
C     POLDENS MAKES A BUNCH OF POLARISED DENSITY PLOTS.                C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      LOGICAL RADPLOT
C
      CHARACTER*5  NMDL
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(1)
C
      DIMENSION POLC(0:NRAD)
      DIMENSION ETMP(0:NRAD),RK(7)
      DIMENSION C0(7),CNS(7),CND(7)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,MFT),XNUC(MCT,MFT),NNUC(MCT),NMDL(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/VPOL/RAD(0:NRAD),VUEH(0:NRAD),EUEH(0:NRAD),PUEH(0:NRAD)
C
C     MULTIPLIER FOR COMPTON WAVELENGTH (USE FOR MUON OR TAUON FIELD)
      CMPF = 1.0D0*CMPW
      
      iz = 1
C
C     RADIAL PLOTTING INFORMATION
      R0 = 0.0D0/CFM
      RN = 100.0D0/CFM
      HR = (RN-R0)/DFLOAT(NRAD)
C
      S0 = 0.0D0
      SM = 0.0D0
31    SN = SN + 0.001D0/CFM
      PM = SN*rhonuc(NMDL(IZ),IZ,SN)
      IF(DABS(PM).GT.1.0D-10) GOTO 31
c      SM = DSQRT(5.0D0/3.0D0)*RNUC(IZ)
      HS = (SN-S0)/DFLOAT(NSRC)
C
CC     CHARGE SOURCE EXTENT AND UNITLESS PARAMETER
C      SN = DSQRT(5.0D0/3.0D0)*RNUC(IZ)
C      X  = 2.0D0*SM/CMPF
C
      RADPLOT = .TRUE.
      IF(.NOT.RADPLOT) GOTO 50
C
C     GENERATE RADIAL GRID AND CHARGE DENSITY
      DO N=0,NRAD
C
C       SET RADIUS
        R = R0 + HR*DFLOAT(N)
        RAD(N) = R
C
C       INITIALISE COUNTER FOR RHOPOL(R)
        RHOPOL = 0.0D0
C
C       MULTIPLICATIVE FACTOR
        RMLT = 8.0D0*R/(3.0D0*CV*CMPF)
C
C       INITIALISE COUNTERS
        RHO1 = 0.0D0
        RHO2 = 0.0D0
C
C       TWO INTEGRALS BUT SEPARATE OUT THE FIRST ONE
        DO M=0,NSRC
C
C         SET SOURCE CHARGE RADIUS
          S = S0 + HS*DFLOAT(M)
C
C         VARIOUS ARGUMENTS
          X1 = 2.0D0*DABS(R-S)/CMPF
          X2 = 2.0D0*DABS(R+S)/CMPF
C
C         ADD TO COUNTERS
          IF(X1.GT.0.0D0) THEN
            QINT = S*RHONUC(NMDL(IZ),IZ,S)*CHIFNCquad(X1,0)
            RHO1 = RHO1 + EXTINT11(QINT,M,NSRC)
          ENDIF
          
          IF(X2.GT.0.0D0) THEN
            QINT = S*RHONUC(NMDL(IZ),IZ,S)*CHIFNCquad(X2,0)
            RHO2 = RHO2 + EXTINT11(QINT,M,NSRC)
          ENDIF
C
        ENDDO
C
        RHOPOL = RHOPOL + 5.0D0*HS*(RHO1-RHO2)/2.99376D+5
c
        POLC(N) = RMLT*RHOPOL
c       WRITE(6, *) rad(n)/cfm,rmlt,polc(n)
        WRITE(6, *) rad(n)*cfm,rho1,rho2,polc(n)
C        
      ENDDO
C
C     FILE NAME AND TITLES
      XOUT   = 'Pueh_uni_exct'
      TITLE  = 'Weighted Uehling radial charge density rho(r)'
      XAXIS  = 'r [fm]'
      YAXIS  = '4*pi*r^2rho(r) [C/fm]'
      KEY(1) = 'Gaussian'
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C       BEGIN LOOP OVER ALL POINTS
        DO N=0,300
C
C         PRINT DISTANCE ALONG PATH AND FOUR-CURRENT
          WRITE(8, *) RAD(N),POLC(N)/CFM
C
        ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,1,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
      RETURN
50    CONTINUE
C
C     CHARGE SOURCE EXTENT AND UNITLESS PARAMETER
      SM = DSQRT(5.0D0/3.0D0)*RNUC(IZ)
      X  = 2.0D0*SM/CMPF
C
C     CHI_N(0) VALUES
      C0(2) = 9.0D0*PI/32.0D0
      C0(3) = 0.4D0
      C0(4) = 5.0D0*PI/64.0D0
      C0(5) = 6.0D0/35.0D0
      C0(6) = 21.0D0*PI/512.0D0
      C0(7) = 32.0D0/315.0D0
C
C     CHI_N(X) AND CHI_N(2X) VALUES
      DO N=2,7
        CNS(N) = CHIFNC(      X,N,0)
        CND(N) = CHIFNC(2.0D0*X,N,0)
      ENDDO
C
C     SEARCH FOR START OF POSITIVE CHARGE
      RZR = DSQRT(5.0D0/3.0D0)*RNUC(IZ)*CFM
C
C     TOTAL CHARGE IN EACH REGION (MAGNITUDE)
c     ***
      QP =-C0(4) + CND(4) + 2.0D0*X*CND(3) + X*X*C0(2) + X*X*CND(2)
      QP = QP/(CV*PI*X*X*X)
C
C     CALCULATION OF MOMENTS
C
C     <R^-2> UNDEFINED FOR NOW
      RK(2) = 0.0D0
C
C     <R^-1>
      RK(3) =-C0(3) + CNS(3) + X*CNS(2)
      RK(3) = 4.0D0*RK(3)/(PI*X*X*X*CFM)
C
C     <R^0>
      RK(4) = 0.0D0
C
C     <R^1>
      RK(5) =-C0(5) + CNS(5) + X*CNS(4) + 0.5D0*X*X*C0(3)
      RK(5) = 2.0D0*RK(5)*CMPF*CFM/(CV*PI*X*X*X)
C
C     <R^2>
      RK(6) = 2.0D0*((CMPF*CFM)**2)/(5.0D0*CV*PI)
C
C     <R^3>
c     ***
      RK(7) =-C0(7) + CNS(7) + X*CNS(6) + 0.5D0*X*X*C0(5)
     &      + 0.125D0*X*X*X*X*C0(3)

c      WRITE(6, *) RNUC(IZ)*CFM,-C0(7),CNS(7),X*CNS(6),0.5D0*X*X*C0(5),
c     &           0.125D0*X*X*X*X*C0(3),rk(7)
c      IF(RNUC(IZ)*CFM.EQ.2.0D0) STOP

      RK(7) = 6.0D0*RK(7)/(CV*PI)
      RK(7) = RK(7)*((CMPF*CFM/X)**3)
C
C     VUEH(0) DEFINED AS IDENTICAL TO <R^-1>
      RK(1) = RK(3)
C
      AVAL = (RNUC(IZ)*CFM - 0.570D0)/0.836D0
      AVAL = AVAL**3
      
      WRITE(6, *) RNUC(IZ)*CFM,AVAL,(RK(I),I=1,7),QP,RZR
C
      RETURN
      END
