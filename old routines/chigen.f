      SUBROUTINE CHIGEN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            CCCCCC  HH    HH IIII GGGGGG  EEEEEEEE NN    NN           C
C           CC    CC HH    HH  II GG    GG EE       NNN   NN           C
C           CC       HH    HH  II GG       EE       NNNN  NN           C
C           CC       HHHHHHHH  II GG       EEEEEE   NN NN NN           C
C           CC       HH    HH  II GG   GGG EE       NN  NNNN           C
C           CC    CC HH    HH  II GG    GG EE       NN   NNN           C
C            CCCCCC  HH    HH IIII GGGGGG  EEEEEEEE NN    NN           C
C                                                                      C
C -------------------------------------------------------------------- C
C  CHIGEN EXPLICITLY GENERATES THE FUNCTION CHI_2(X) AND SAVES IT TO   C
C  FILE, FOR LATER SPLINE INTERPOLATION IN UEHLING CALCULATIONS.       C
C  IT USES THE EXPONENTIAL INTEGRAL EXPANSION METHOD.                  C
C  THE SPLINE FUNCTION IS ACTUALLY χ_2(X)*EXP(X), LEAVING A POLYNOMIAL.C
C**********************************************************************C
      INCLUDE 'param.h'
C
      DIMENSION XS(0:NLW),Y1S(0:NLW),D1S(0:NLW),Y2S(0:NLW),D2S(0:NLW),
     &          XB(0:NUP),Y1B(0:NUP),D1B(0:NUP),Y2B(0:NUP),D2B(0:NUP),
     &                    YOS(0:NLW),DOS(0:NLW)
C
C     JOINING POINT FOR BOTH SPLINE FUNCTIONS
      XSPL = 0.60D0
C
C**********************************************************************C
C     INTERPOLATING FUNCTIONS IN REGION 0.0D0 < X <~ 0.6D0             C
C**********************************************************************C
C
C     REGION GRID LENGTHS AND ENDPOINTS
      NCOR = 100
      XBEG = 0.00D0
      XCOR = 1.00D-4
      XEND = 1.25D0*XSPL
C
C     VERY SMALL X: UNIFORMLY-SPACED GRID
      DO N=0,NCOR
        XS(N)  = XBEG + (XCOR-XBEG)*DFLOAT(N)/DFLOAT(NCOR)
        IF(N.EQ.0) THEN
          Y1S(N) = 1.0D30
        ELSE
          Y1S(N) = CHIFNC(XS(N),1,0)*DEXP(XS(N))
        ENDIF
        YOS(N) = CHI1SF(XS(N))
        Y2S(N) = CHIFNC(XS(N),2,0)*DEXP(XS(N))
      ENDDO
C
C     LARGER X: EXPONENTIALLY-SPACED GRID
      SPC = DLOG(XEND/XCOR)/DFLOAT(NLW-NCOR)
      DO N=0,NLW-NCOR
        XS(N+NCOR)  = XCOR*DEXP(N*SPC)
        YOS(N+NCOR) = CHI1SF(XS(N+NCOR))
        Y1S(N+NCOR) = CHIFNC(XS(N+NCOR),1,0)*DEXP(XS(N+NCOR))
        Y2S(N+NCOR) = CHIFNC(XS(N+NCOR),2,0)*DEXP(XS(N+NCOR))
      ENDDO
C
C     DERIVATIVES AT SPLINE ENDS
      DO0 = 1.0D30
      DON = CHI1SD(XEND)
      D10 =-1.0D30
      D1N = (CHIDRV(XEND,1,0)+CHIFNC(XEND,1,0))*DEXP(XEND)
      D20 =-1.0D30
      D2N = (CHIDRV(XEND,2,0)+CHIFNC(XEND,2,0))*DEXP(XEND)
C
C     GENERATE SPLINES
      CALL SPLNGEN(XS,YOS,DOS,DO0,DON,NLW)
      CALL SPLNGEN(XS,Y1S,D1S,D10,D1N,NLW)
      CALL SPLNGEN(XS,Y2S,D2S,D20,D2N,NLW)
C
C     WRITE SPLINE DATA TO FILE
      OPEN(UNIT=50,FILE='chi1s.dat',STATUS='UNKNOWN')
      REWIND(UNIT=50)
        WRITE(50,*) XSPL
        DO N=0,NLW
          WRITE(50,*) XS(N),YOS(N),DOS(N)
        ENDDO
      CLOSE(UNIT=50)
C
C     WRITE SPLINE DATA TO FILE
      OPEN(UNIT=50,FILE='chi1_small.dat',STATUS='UNKNOWN')
      REWIND(UNIT=50)
        WRITE(50,*) XSPL
        DO N=0,NLW
          WRITE(50,*) XS(N),Y1S(N),D1S(N)
        ENDDO
      CLOSE(UNIT=50)
C
C     WRITE SPLINE DATA TO FILE
      OPEN(UNIT=52,FILE='chi2_small.dat',STATUS='UNKNOWN')
      REWIND(UNIT=52)
        WRITE(52,*) XSPL
        DO N=0,NLW
          WRITE(52,*) XS(N),Y2S(N),D2S(N)
        ENDDO
      CLOSE(UNIT=52)
C
C**********************************************************************C
C     INTERPOLATING FUNCTIONS IN REGION X >~ 0.6D0                     C
C**********************************************************************C
C
C     LARGE X ENDPOINTS
      XBEG =   0.90D0*XSPL
      XEND = 250.00D0
C
C     USE EXPONENTIALLY-SPACED GRID
      SPC = DLOG(XEND/XBEG)/DFLOAT(NUP)
      DO N=0,NUP
        XB(N)  = XBEG*DEXP(N*SPC)
        Y1B(N) = CHIFNC(XB(N),1,1)*DEXP(XB(N))
        Y2B(N) = CHIFNC(XB(N),2,1)*DEXP(XB(N))
      ENDDO
C
C     DERIVATIVES AT SPLINE ENDS
      D10 = (CHIDRV(XBEG,1,1)+CHIFNC(XBEG,1,1))*DEXP(XBEG)
      D1N = (CHIDRV(XEND,1,1)+CHIFNC(XEND,1,1))*DEXP(XEND)
      D20 = (CHIDRV(XBEG,2,1)+CHIFNC(XBEG,2,1))*DEXP(XBEG)
      D2N = (CHIDRV(XEND,2,1)+CHIFNC(XEND,2,1))*DEXP(XEND)
C
C     GENERATE SPLINE
      CALL SPLNGEN(XB,Y1B,D1B,D10,D1N,NUP)
      CALL SPLNGEN(XB,Y2B,D2B,D20,D2N,NUP)
C
C     WRITE SPLINE DATA TO FILE
      OPEN(UNIT=51,FILE='chi1_big.dat',STATUS='UNKNOWN')
      REWIND(UNIT=51)
        WRITE(51,*) XSPL
        DO N=0,NUP
          WRITE(51,*) XB(N),Y1B(N),D1B(N)
        ENDDO
      CLOSE(UNIT=51)
C
C     WRITE SPLINE DATA TO FILE
      OPEN(UNIT=53,FILE='chi2_big.dat',STATUS='UNKNOWN')
      REWIND(UNIT=53)
        WRITE(53,*) XSPL
        DO N=0,NUP
          WRITE(53,*) XB(N),Y2B(N),D2B(N)
        ENDDO
      CLOSE(UNIT=53)
C
      RETURN
      END
C
C
      FUNCTION CHIFNC(X,N,MODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           CCCCCC  HH    HH IIII FFFFFFFF NN    NN  CCCCCC            C
C          CC    CC HH    HH  II  FF       NNN   NN CC    CC           C          
C          CC       HH    HH  II  FF       NNNN  NN CC                 C
C          CC       HHHHHHHH  II  FFFFFF   NN NN NN CC                 C
C          CC       HH    HH  II  FF       NN  NNNN CC                 C
C          CC    CC HH    HH  II  FF       NN   NNN CC    CC           C
C           CCCCCC  HH    HH IIII FF       NN    NN  CCCCCC            C
C                                                                      C
C -------------------------------------------------------------------- C      
C     EVALUATE THE FUNCTION CHI_N (X) FOR PARTICULAR VALUE OF X,       C
C     USING A TRUNCATED SERIES WITH NTR TERMS.                         C
C**********************************************************************C
      INCLUDE 'param.h'
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     CHECK THAT N IS WITHIN ALLOWED VALUES
      IF(N.LE.1.AND.X.EQ.0.0D0) THEN
        WRITE(6, *) 'In CHIFNC: divergent behaviour when N < 1. N = ',N
        WRITE(7, *) 'In CHIFNC: divergent behaviour when N < 1. N = ',N
        CHIFNC = 1.0D30
      ENDIF
C
C     INTEGRATION PARAMETERS (UNIFORM GRID)
      TBEG = 1.0D0
      TEND = 1.0D4
      HT   = (TEND/TBEG)/DFLOAT(KTR)
C
      CHIFNC = 0.0D0
      DO K=0,KTR
        T = TBEG + HT*DFLOAT(K)
        TR = 1.0D0/T
        Y1 = TR**N
        Y2 = DSQRT(1.0D0 - TR*TR)
        Y3 = 1.0D0 + 0.5D0*TR*TR
        Y4 = DEXP(-X*T)
        Z  = Y1*Y2*Y3*Y4
        CHIFNC = CHIFNC + 5.0D0*HT*EXTINT11(Z,K,KTR)/2.99376D+5
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION CHIFNCSER(X,N,MODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           CCCCCC  HH    HH IIII FFFFFFFF NN    NN  CCCCCC            C
C          CC    CC HH    HH  II  FF       NNN   NN CC    CC           C          
C          CC       HH    HH  II  FF       NNNN  NN CC                 C
C          CC       HHHHHHHH  II  FFFFFF   NN NN NN CC                 C
C          CC       HH    HH  II  FF       NN  NNNN CC                 C
C          CC    CC HH    HH  II  FF       NN   NNN CC    CC           C
C           CCCCCC  HH    HH IIII FF       NN    NN  CCCCCC            C
C                                                                      C
C -------------------------------------------------------------------- C      
C     EVALUATE THE FUNCTION CHI_N (X) FOR PARTICULAR VALUE OF X,       C
C     USING A TRUNCATED SERIES WITH NTR TERMS.                         C
C                   ∞      _______                                     C
C         χ_N(X) = ∫ t^-N √1-1/t^2 (1 + 1/2t^2) exp(-x t) dt.          C
C                   1                                                  C
C                                                                      C
C -------------------------------------------------------------------- C
C     MODE = 0: USE EXACT AND APPROX. CHI_N (0.0D0) TO MIN. ERROR      C
C     MODE = 1: DO NOT USE EXACT VALUE (LARGE X).                      C
C**********************************************************************C
      INCLUDE 'param.h'
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     CHECK THAT N IS WITHIN ALLOWED VALUES
      IF(N.LE.1.AND.X.EQ.0.0D0) THEN
        WRITE(6, *) 'In CHIFNC: divergent behaviour when N < 1. N = ',N
        WRITE(7, *) 'In CHIFNC: divergent behaviour when N < 1. N = ',N
        STOP
      ENDIF
C
C     SPECIAL CASE: X = 0.0D0
      IF(X.EQ.0.0D0) THEN
        CHIFNC = 0.75D0*PI12*GAMHLF(N+3)/(DFLOAT(N-1)*GAMHLF(N+4))
        RETURN
      ENDIF
C
C     χ_N(X) WITH SUM OF EXPONENTIAL FUNCTIONS
C
C     INITIALISE THE CHIFNC COUNTER
      CHIFNC = 0.0D0
C
C     LOOP FROM K=0 TO TRUNCATED INFINITY
      PHS =-1.0D0
      DO K=0,KTR
C
        RK = DFLOAT(K)
        IF(K.EQ.0) THEN
          BNM = 1.0D0
        ELSE
          BNM = (1.5D0-RK)*BNM/RK
        ENDIF
C
        PHS =-PHS
        ENK = EXPINTE(2*K+N,X) + 0.5D0*EXPINTE(2*K+N+2,X)
C
        CHIFNC = CHIFNC + PHS*BNM*ENK
C
      ENDDO
C
C     χ_N(X) NEAR ZERO, TAKING A CLEVER DIFFERENCE FROM χ_N(0)
      IF(MODE.EQ.0.AND.N.GT.1) THEN
C
C       EXACT VALUE AT X = 0.0D0
        TRU = 0.75D0*PI12*GAMHLF(N+3)/(DFLOAT(N-1)*GAMHLF(N+4))
C
C       INITIALISE THE COUNTER FOR RESIDUAL
        RSD = 0.0D0
C
        PHS =-1.0D0
        DO K=0,KTR
C
          RK = DFLOAT(K)
          IF(K.EQ.0) THEN
            BNM = 1.0D0
          ELSE
            BNM = (1.5D0-RK)*BNM/RK
          ENDIF
C
          PHS =-PHS
          RAT = DFLOAT(6*K+3*N+1)/DFLOAT(2*(2*K+N-1)*(2*K+N+1))
C
          RSD = RSD + PHS*BNM*RAT
C
        ENDDO
C
        CHIFNC = CHIFNC+TRU-RSD
C
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION CHIDRV(X,N,MODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           CCCCCC  HH    HH IIII DDDDDDD  RRRRRRR  VV    VV           C
C          CC    CC HH    HH  II  DD    DD RR    RR VV    VV           C
C          CC       HH    HH  II  DD    DD RR    RR VV    VV           C
C          CC       HHHHHHHH  II  DD    DD RR    RR VV    VV           C
C          CC       HH    HH  II  DD    DD RRRRRRR   VV  VV            C
C          CC    CC HH    HH  II  DD    DD RR    RR   VVVV             C
C           CCCCCC  HH    HH IIII DDDDDDD  RR    RR    VV              C
C                                                                      C
C -------------------------------------------------------------------- C
C     CHI2DERIV EVALUATES THE INSTANTANEOUS DERIVATIVE OF CHI_N(X)     C
C     USING A SERIES OF EXPONENTIAL INTEGRALS, OF LENGTH NTR.          C
C**********************************************************************C
      INCLUDE 'param.h'
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     CHECK THAT N IS WITHIN ALLOWED VALUES
      IF(N.LE.2.AND.X.EQ.0.0D0) THEN
        WRITE(6, *) 'In CHIDRV: divergent behaviour when N < 2. N = ',N
        WRITE(7, *) 'In CHIDRV: divergent behaviour when N < 2. N = ',N
        STOP
      ENDIF
C
C     SPECIAL CASE: X = 0.0D0
      IF(X.EQ.0.0D0) THEN
        CHIDRV =-0.75D0*PI12*GAMHLF(N+2)/(DFLOAT(N-2)*GAMHLF(N+3))
        RETURN
      ENDIF
C
C     χ'_N(X) WITH SUM OF EXPONENTIAL FUNCTIONS
C
C     INITIALISE THE CHIDRV COUNTER
      CHIDRV = 0.0D0
C
C     LOOP FROM K=0 TO TRUNCATED INFINITY
      PHS =-1.0D0
      DO K=0,KTR
C
        RK = DFLOAT(K)
        IF(K.EQ.0) THEN
          BNM = 1.0D0
        ELSE
          BNM = (1.5D0-RK)*BNM/RK
        ENDIF
C
        PHS =-PHS
        ENK = EXPINTE(2*K+N-1,X) + 0.5D0*EXPINTE(2*K+N+1,X)
C
        CHIDRV = CHIDRV - PHS*BNM*ENK
C
      ENDDO
C
C     χ'_N(X) NEAR ZERO, TAKING A CLEVER DIFFERENCE FROM χ'_N(0)
      IF(MODE.EQ.0.AND.N.GT.2) THEN
C
C       EXACT VALUE AT X = 0.0D0
        TRU =-0.75D0*PI12*GAMHLF(N+2)/(DFLOAT(N-2)*GAMHLF(N+3))
C
C       INITIALISE THE COUNTER FOR RESIDUAL
        RSD = 0.0D0
C
        PHS =-1.0D0
        DO K=0,KTR
C
          RK = DFLOAT(K)
          IF(K.EQ.0) THEN
            BNM = 1.0D0
          ELSE
            BNM = (1.5D0-RK)*BNM/RK
          ENDIF
C
          PHS =-PHS
          RAT = DFLOAT(6*K+3*N-2)/DFLOAT(2*(2*K+N-2)*(2*K+N))
C
          RSD = RSD - PHS*BNM*RAT
C
        ENDDO
C
        CHIDRV = CHIDRV+TRU-RSD
C
      ENDIF
C
      RETURN
      END
