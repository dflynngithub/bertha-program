      SUBROUTINE EILS3(ELS,EXL,XYZ,KQN,MQN,NBS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C               EEEEEEEE IIII LL       SSSSSS   333333                 C
C               EE        II  LL      SS    SS 33    33                C
C               EE        II  LL      SS             33                C
C               EEEEEE    II  LL       SSSSSS    33333                 C
C               EE        II  LL            SS       33                C
C               EE        II  LL      SS    SS 33    33                C
C               EEEEEEEE IIII LLLLLLLL SSSSSS   333333                 C
C                                                                      C
C -------------------------------------------------------------------- C
C  EILS3 GENERATES A VECTOR BLOCK OF RAW EILS-COEFFICIENTS FOR A       C
C  PARTICULAR COMBINATION OF MQNS (WITH SIGN), TO BE USED BY EILSB3.   C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL - PAIR OF BASIS SET PARAMETER LISTS.                          C
C  ▶ XYZ - PAIR OF BASIS FUNCTION CARTESIAN CENTRES.                   C
C  ▶ KQN - PAIR OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.               C
C  ▶ MQN - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.            C
C  ▶ NBS - PAIR OF BASIS FUNCTION BLOCK LENGTHS.                       C
C  ▶ IQ  - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.             C
C  OUTPUT:                                                             C
C  ▶ ELS - RAW EQ-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).              C
C**********************************************************************C
      INCLUDE 'param.h'
C
      DIMENSION LQLAB(2),MQLAB(2)
      DIMENSION EXL(MBS,2),XYZ(3,2),KQN(2),LQN(2),JQN(2),MQN(2),NBS(2)
      DIMENSION RN(MBS,2)
      DIMENSION T2(MB2),T0(MB2)
C
      COMPLEX*16 CONE
      COMPLEX*16 ELS(MB2,MEQ,3),ESG(MB2,MEQ),ENSG(MB2,MEQ)
C
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     CLEBSCH-GORDAN SENSITIVITY PARAMETER
      DATA SENS/1.0D-14/
C
C     INITIALISE THE COEFFICIENTS TO ZERO
      DO M=1,MB2
        DO ITUV=1,MEQ
          DO IQ=1,3
            ELS(M,ITUV,IQ) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C     CALCULATE LQN VALUES
      IF(KQN(1).LT.0) THEN
        LQN(1) =-KQN(1)-1
      ELSE
        LQN(1) = KQN(1)
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        LQN(2) =-KQN(2)-1
      ELSE
        LQN(2) = KQN(2)
      ENDIF
C
C     CALCULATE JQN VALUES
      JQN(1) = 2*IABS(KQN(1))-1
      JQN(2) = 2*IABS(KQN(2))-1
C
C     CALCULATE THE APPROPRIATE CAB GORDAN FACTORS
      IF(KQN(1).LT.0) THEN
        CAU = DSQRT(DFLOAT(JQN(1)+MQN(1)  )/DFLOAT(2*JQN(1)  ))
        CAL = DSQRT(DFLOAT(JQN(1)-MQN(1)  )/DFLOAT(2*JQN(1)  ))
      ELSE
        CAU =-DSQRT(DFLOAT(JQN(1)-MQN(1)+2)/DFLOAT(2*JQN(1)+4))
        CAL = DSQRT(DFLOAT(JQN(1)+MQN(1)+2)/DFLOAT(2*JQN(1)+4))
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        CBU =-DSQRT(DFLOAT(JQN(2)-MQN(2)+2)/DFLOAT(2*JQN(2)+4))
        CBL = DSQRT(DFLOAT(JQN(2)+MQN(2)+2)/DFLOAT(2*JQN(2)+4))
      ELSE
        CBU = DSQRT(DFLOAT(JQN(2)+MQN(2)  )/DFLOAT(2*JQN(2)  ))
        CBL = DSQRT(DFLOAT(JQN(2)-MQN(2)  )/DFLOAT(2*JQN(2)  ))
      ENDIF
C
C     DETERMINE THE NUMBER OF PAIRS IN THIS BLOCK
      MAXM = NBS(1)*NBS(2)
C
C     KINETIC PRE-FACTORS FOR THIS BLOCK
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
          T2(M) =-2.0D0*EXL(JBAS,2)
          T0(M) = DFLOAT(2*LQN(2)+1)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INITIALISE COMMON GEOMETRIC INFORMATION                          C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTRES A AND B
      DX  = (XYZ(1,1)-XYZ(1,2))**2
      DY  = (XYZ(2,1)-XYZ(2,2))**2
      DZ  = (XYZ(3,1)-XYZ(3,2))**2
      AB2 = DX+DY+DZ
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
C
          P(M)    = EXL(IBAS,1)+EXL(JBAS,2)
          P2(M)   = 2.0D0*P(M)
          P22(M)  = P2(M)*P2(M)
          PX      = (EXL(IBAS,1)*XYZ(1,1)+EXL(JBAS,2)*XYZ(1,2))/P(M)
          PY      = (EXL(IBAS,1)*XYZ(2,1)+EXL(JBAS,2)*XYZ(2,2))/P(M)
          PZ      = (EXL(IBAS,1)*XYZ(3,1)+EXL(JBAS,2)*XYZ(3,2))/P(M)
          PAX(M)  = PX-XYZ(1,1)
          PAY(M)  = PY-XYZ(2,1)
          PAZ(M)  = PZ-XYZ(3,1)
          PBX(M)  = PX-XYZ(1,2)
          PBY(M)  = PY-XYZ(2,2)
          PBZ(M)  = PZ-XYZ(3,2)
          PA2(M)  = PAX(M)*PAX(M) + PAY(M)*PAY(M) + PAZ(M)*PAZ(M)
          PB2(M)  = PBX(M)*PBX(M) + PBY(M)*PBY(M) + PBZ(M)*PBZ(M)
          RKAB(M) = DEXP(-EXL(IBAS,1)*EXL(JBAS,2)*AB2/P(M))
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CASE 1: KQN(2).LT.0                                              C
C**********************************************************************C
C
      IF(KQN(2).GT.0) GOTO 100
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)
      LQLAB(2) = LQN(2)+1
C
C >>  TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MQLAB(1) = (MQN(1)-1)/2
      MQLAB(2) = (MQN(2)-1)/2
C
C     CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
      CAB = CAU*CBU
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM0  = LQLAB(1)+LQLAB(2)
        NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C       CONTRIBUTION TO ELSZ: UPPER-UPPER
        DO M=1,MAXM
          TK = CAB*T2(M)
          DO ITUV=1,NTUV0
            ELS(M,ITUV,3) = ELS(M,ITUV,3) + TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
C >>  TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MQLAB(1) = (MQN(1)+1)/2
      MQLAB(2) = (MQN(2)+1)/2
C
C     CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
      CAB = CAL*CBL
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM0  = LQLAB(1)+LQLAB(2)
        NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C       CONTRIBUTION TO ELSZ: LOWER-LOWER
        DO M=1,MAXM
          TK = CAB*T2(M)
          DO ITUV=1,NTUV0
            ELS(M,ITUV,3) = ELS(M,ITUV,3) - TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ELSX AND ELSY
C
C >>  TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MQLAB(1) = (MQN(1)-1)/2
      MQLAB(2) = (MQN(2)+1)/2
C
C     CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
      CAB = CAU*CBL
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM0  = LQLAB(1)+LQLAB(2)
        NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C       CONTRIBUTION TO ELSX AND ELSY: UPPER-LOWER
        DO M=1,MAXM
          TK = CAB*T2(M)
          DO ITUV=1,NTUV0
            ELS(M,ITUV,1) = ELS(M,ITUV,1) +      TK*ESG(M,ITUV)
            ELS(M,ITUV,2) = ELS(M,ITUV,2) - CONE*TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
C >>  TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MQLAB(1) = (MQN(1)+1)/2
      MQLAB(2) = (MQN(2)-1)/2
C
C     CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
      CAB = CAL*CBU
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM0  = LQLAB(1)+LQLAB(2)
        NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C       CONTRIBUTION TO ELSX AND ELSY: LOWER-UPPER
        DO M=1,MAXM
          TK = CAB*T2(M)
          DO ITUV=1,NTUV0
            ELS(M,ITUV,1) = ELS(M,ITUV,1) +      TK*ESG(M,ITUV)
            ELS(M,ITUV,2) = ELS(M,ITUV,2) + CONE*TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
100   CONTINUE
C
C**********************************************************************C
C     CASE 2: KQN(2).GT.0                                              C
C**********************************************************************C
C
      IF(KQN(2).LT.0) GOTO 200
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)
      LQLAB(2) = LQN(2)-1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ELS0 AND ELSZ
C
C >>  TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MQLAB(1) = (MQN(1)-1)/2
      MQLAB(2) = (MQN(2)-1)/2
C
C     CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
      CAB = CAU*CBU
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM0  = LQLAB(1)+LQLAB(2)
        NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C       CONTRIBUTION TO ELSZ: UPPER-UPPER FOR [N=0,N'=0]
        DO M=1,MAXM
          TK = CAB*T0(M)
          DO ITUV=1,NTUV0
            ELS(M,ITUV,3) = ELS(M,ITUV,3) + TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
C       INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
        CALL STEPN(ESG,ENSG,LAM0,MAXM,2)
C
C       INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
        LAM2  = LQLAB(1)+LQLAB(2)+2
        NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C       CONTRIBUTION TO ELSZ: UPPER-UPPER FOR [N=0,N'=1]
        DO M=1,MAXM
          TK = CAB*T2(M)
          DO ITUV=1,NTUV2
            ELS(M,ITUV,3) = ELS(M,ITUV,3) + TK*ENSG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
C >>  TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MQLAB(1) = (MQN(1)+1)/2
      MQLAB(2) = (MQN(2)+1)/2
C
C     CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
      CAB = CAL*CBL
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM0  = LQLAB(1)+LQLAB(2)
        NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C       CONTRIBUTION TO ELSZ: LOWER-LOWER FOR [N=0,N'=0]
        DO M=1,MAXM
          TK = CAB*T0(M)
          DO ITUV=1,NTUV0
            ELS(M,ITUV,3) = ELS(M,ITUV,3) - TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
C       INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
        CALL STEPN(ESG,ENSG,LAM0,MAXM,2)
C
C       INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
        LAM2  = LQLAB(1)+LQLAB(2)+2
        NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C       CONTRIBUTION TO ELSZ: LOWER-LOWER FOR [N=0,N'=1]
        DO M=1,MAXM
          TK = CAB*T2(M)
          DO ITUV=1,NTUV2
            ELS(M,ITUV,3) = ELS(M,ITUV,3) - TK*ENSG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ELSX AND ELSY
C
C >>  TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MQLAB(1) = (MQN(1)-1)/2
      MQLAB(2) = (MQN(2)+1)/2
C
C     CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
      CAB = CAU*CBL
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM0  = LQLAB(1)+LQLAB(2)
        NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C       CONTRIBUTION TO ELSX AND ELSY: UPPER-LOWER FOR [N=0,N'=0]
        DO M=1,MAXM
          TK = CAB*T0(M)
          DO ITUV=1,NTUV0
            ELS(M,ITUV,1) = ELS(M,ITUV,1) +      TK*ESG(M,ITUV)
            ELS(M,ITUV,2) = ELS(M,ITUV,2) - CONE*TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
C       INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
        CALL STEPN(ESG,ENSG,LAM0,MAXM,2)
C
C       INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
        LAM2  = LQLAB(1)+LQLAB(2)+2
        NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C       CONTRIBUTION TO ELSX AND ELSY: UPPER-LOWER FOR [N=0,N'=1]
        DO M=1,MAXM
          TK = CAB*T2(M)
          DO ITUV=1,NTUV2
            ELS(M,ITUV,1) = ELS(M,ITUV,1) +      TK*ENSG(M,ITUV)
            ELS(M,ITUV,2) = ELS(M,ITUV,2) - CONE*TK*ENSG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
C >>  TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MQLAB(1) = (MQN(1)+1)/2
      MQLAB(2) = (MQN(2)-1)/2
C
C     CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
      CAB = CAL*CBU
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM0  = LQLAB(1)+LQLAB(2)
        NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C       CONTRIBUTION TO ELSX AND ELSY: LOWER-UPPER FOR [N=0,N'=0]
        DO M=1,MAXM
          TK = CAB*T0(M)
          DO ITUV=1,NTUV0
            ELS(M,ITUV,1) = ELS(M,ITUV,1) +      TK*ESG(M,ITUV)
            ELS(M,ITUV,2) = ELS(M,ITUV,2) + CONE*TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
C       INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
        CALL STEPN(ESG,ENSG,LAM0,MAXM,2)
C
C       INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
        LAM2  = LQLAB(1)+LQLAB(2)+2
        NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C       CONTRIBUTION TO ELSX AND ELSY: LOWER-UPPER FOR [N=0,N'=1]
        DO M=1,MAXM
          TK = CAB*T2(M)
          DO ITUV=1,NTUV2
            ELS(M,ITUV,1) = ELS(M,ITUV,1) +      TK*ENSG(M,ITUV)
            ELS(M,ITUV,2) = ELS(M,ITUV,2) + CONE*TK*ENSG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
200   CONTINUE
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
C     MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
      LAM2  = LQN(1)+LQN(2)+1
      NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C     GENERATE RNLS NORMALISATION CONSTANTS
      CALL RNLS(RN,EXL,LQN,NBS)
C
C     NORMALISE THE ELLQ COEFFICIENT BLOCK
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M   = M+1
          RLS = RN(IBAS,1)*RN(JBAS,2)
          DO ITUV=1,NTUV2
            DO IQ=1,3
              ELS(M,ITUV,IQ) = RLS*ELS(M,ITUV,IQ)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     BRING THESE COEFFICIENTS INTO THE ELSQ VALUES AND ALSO FACTOR i
      DO IQ=1,3
        DO M=1,MAXM
          DO ITUV=1,NTUV2
            ELS(M,ITUV,IQ) = CONE*ELS(M,ITUV,IQ)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END

