      SUBROUTINE RMAKE(RC,QP,APH,MAXM,LAMBDA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C           RRRRRRR  MM       MM    AA    KK    KK EEEEEEEE           C
C           RR    RR MMM     MMM   AAAA   KK   KK  EE                 C
C           RR    RR MMMM   MMMM  AA  AA  KK  KK   EE                 C
C           RR    RR MM MM MM MM AA    AA KKKKK    EEEEEE             C
C           RRRRRRR  MM  MMM  MM AAAAAAAA KK  KK   EE                 C
C           RR    RR MM   M   MM AA    AA KK   KK  EE                 C
C           RR    RR MM       MM AA    AA KK    KK EEEEEEEE           C
C                                                                     C
C ------------------------------------------------------------------- C
C     RMAKE GENERATES A COMPLETE SET OF R-INTEGRALS REQUIRED IN THE   C
C     FINITE SUM REPRESENTATION OF A MULTI-CENTRE GAUSSIAN OVERLAP.   C
C*********************************************************************C
      PARAMETER(MKP=9,MBS=26,MB2=MBS*MBS,IL4=2*(MKP-1),
     &          MLL=(MKP+8)*(MKP+9)*(MKP+10)/6,MLM=30,MRC=969)
C     
      DIMENSION FS(MB2,MLM),APH(MB2),QP(MB2,3),RC(MB2,MRC),RC2(MB2,MRC)
      DIMENSION F0(MB2,MLM),F1(MB2,MLM),F2(MB2,MLM),F3(MB2,MLM)     
      DIMENSION X0(MB2),X1(MB2),X2(MB2),X3(MB2),
     &          I0(MB2),I1(MB2),I2(MB2),I3(MB2)
C
      COMMON/ACSS/INABCD(0:IL4,0:IL4,0:IL4),
     &            IVEC(MLL),JVEC(MLL),KVEC(MLL),LAMVEC(MLL)
C
C*********************************************************************C
C     THE FIRST STEP OF THIS SUBROUTINE IS TO EVALUATE THE REQUIRED   C
C     BOYS INTEGRALS. THIS IS DIVIDED INTO CASES, DEPENDING ON THE    C
C     MAGNITUDE OF THE ARGUMENT.                                      C
C ------------------------------------------------------------------- C
C        FS_M (X) = INT_{0}^{1} T^{2M} EXP(-X*T^{2}) DT               C
C ------------------------------------------------------------------- C
C   EVALUATED FOR ALL VALUES OF M IN THE RANGE 0 < M < LAMBDA.        C
C             FOR ALL VALUES OF X IN THE RANGE X > 0.                 C
C*********************************************************************C
C
      N0 = 0
      N1 = 0
      N2 = 0
      N3 = 0
C
C     FOR EACH PAIR OF BASIS FUNCTIONS (EXPONENTS EI AND EJ IN 'M'),
C     DETERMINE THE BEST WAY TO EVALUATE THE BOYS FUNCTION
      DO M=1,MAXM
C       X = EIJ*(P-C)^2
        X = APH(M)*(QP(M,1)*QP(M,1)+QP(M,2)*QP(M,2)+QP(M,3)*QP(M,3))
C       CASE 0: IF X IS ALMOST ZERO (SO WHEN Q=P OR EIJ<<1)
        IF(X.LE.1.00D-11) THEN
          N0     = N0 + 1
          X0(N0) = X
          I0(N0) = M
C       CASE 1: IF X IS SMALLER THAN X=18
        ELSEIF(X.GT.1.00D-11.AND.X.LE.1.80D+01) THEN
          N1     = N1 + 1
          X1(N1) = X
          I1(N1) = M
C       CASE 2: IF X IS SMALLER THAN ABOUT 30
        ELSEIF(X.GT.1.80D+01.AND.X.LE.3.00D+01) THEN
          N2     = N2 + 1
          X2(N2) = X
          I2(N2) = M
C       CASE 3: IF X IS LARGER THAN ABOUT 30
        ELSE
          N3     = N3 + 1
          X3(N3) = X
          I3(N3) = M
        ENDIF
      ENDDO
C
C*********************************************************************C
C     CASE 0: ARGUMENT OF THE BOYS FUNCTION IS X = 0.                 C
C             THE VALUE OF THIS FUNCTION IS 2N+1 (DONE IN FUNFX).     C
C*********************************************************************C 
C
      IF(N0.NE.0) THEN
        CALL FUNFX(F0,X0,N0,LAMBDA,0)     
        DO JJ=1,LAMBDA+1
          DO M=1,N0
            FS(I0(M),JJ) = F0(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C*********************************************************************C
C     CASE 1: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=18.     C
C             EVALUATE WITH LOCAL POLYNOMIAL EXPANSION OF ORDER 5,    C
C             AND RECURRENCE IN DIRECTION OF DECREASING M.            C
C*********************************************************************C
C
      IF(N1.NE.0) THEN
        CALL FUNFX(F1,X1,N1,LAMBDA,1)
        DO JJ=1,LAMBDA+1
          DO M=1,N1
            FS(I1(M),JJ) = F1(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C*********************************************************************C
C     CASE 2: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=30.     C
C             EVALUATE USING ASYMPTOTIC FORMULA WITH EXPONENTIAL.     C
C*********************************************************************C
C
      IF(N2.NE.0) THEN
        CALL FUNFX(F2,X2,N2,LAMBDA,2)
        DO JJ=1,LAMBDA+1
          DO M=1,N2
            FS(I2(M),JJ) = F2(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C    
C*********************************************************************C
C     CASE 3: ARGUMENT OF THE BOYS FUNCTION IS LARGER THAN X=30.      C
C             EVALUATE USING ASYMPTOTIC FORMULA WITHOUT EXPONENTIAL.  C
C*********************************************************************C
C
      IF(N3.NE.0) THEN
        CALL FUNFX(F3,X3,N3,LAMBDA,3)
        DO JJ=1,LAMBDA+1
          DO M=1,N3
            FS(I3(M),JJ) = F3(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C*********************************************************************C
C     THE SECOND STEP OF THIS SUBROUTINE IS TO EVALUATE THE ACTUAL    C
C     R-COEFFICIENTS, BASED ON THE CORRESPONDING BOYS FUNCTIONS.      C
C*********************************************************************C
   
C*********************************************************************C
C     CONSTRUCT TOP-LEVEL                                             C
C*********************************************************************C
      DO M=1,MAXM
        RC(M,1) = ((-2.0D0*APH(M))**(LAMBDA))*FS(M,LAMBDA+1)
      ENDDO
C
      IF(MOD(LAMBDA,2).EQ.0) THEN
        ITUVMIN = 1
      ELSE
        ITUVMIN = 2
      ENDIF
C
      ITUV = -1
      DO ILEVEL=LAMBDA-1,ITUVMIN,-2
        ITUV = ITUV+1
        DO IT=0,ITUV
          RIT = DFLOAT(IT)
        DO IU=0,ITUV-IT
          RIU = DFLOAT(IU)
        DO IV=0,ITUV-IT-IU
          RIV = DFLOAT(IV)
C
          N1 = INABCD(IT+1,IU  ,IV  )
          N2 = INABCD(IT  ,IU+1,IV  )
          N3 = INABCD(IT  ,IU  ,IV+1)
C
          IF(IT.NE.0) THEN
            IF(IU.NE.0) THEN
              IF(IV.NE.0) THEN      
C               CASE (1 1 1)
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
              ELSE
C               CASE (1 1 0)
                K1 = INABCD(IT  ,IU  ,IV)
                M1 = INABCD(IT-1,IU  ,IV)
                M2 = INABCD(IT  ,IU-1,IV)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
                ENDDO
              ENDIF
            ELSE
              IF(IV.NE.0) THEN
C               CASE (1 0 1)
                K1 = INABCD(IT  ,IU,IV  )
                M1 = INABCD(IT-1,IU,IV  )
                M3 = INABCD(IT  ,IU,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
              ELSE
C               CASE (1 0 0)
                K1 = INABCD(IT  ,IU,IV)
                M1 = INABCD(IT-1,IU,IV)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
                ENDDO
              ENDIF
            ENDIF
          ELSE
            IF(IU.NE.0) THEN
              IF(IV.NE.0) THEN
C               CASE (0 1 1) 
                K1 = INABCD(IT,IU  ,IV  )
                M2 = INABCD(IT,IU-1,IV  )
                M3 = INABCD(IT,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
              ELSE
C               CASE (0 1 0)
                K1 = INABCD(IT,IU  ,IV)
                M2 = INABCD(IT,IU-1,IV)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
                ENDDO
              ENDIF
            ELSE
              IF(IV.NE.0) THEN
C               CASE (0 0 1)
                K1 = INABCD(IT,IU,IV  )
                M3 = INABCD(IT,IU,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
              ELSE
C               CASE (0 0 0)
                K1 = INABCD(IT,IU,IV)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
                ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        ENDDO
        ENDDO
C       ADD IN (TI=0,IU=0,IV=0) CASE
        DO M=1,MAXM
          RC2(M,1) = ((-2.0D0*APH(M))**(ILEVEL))*FS(M,ILEVEL+1)
        ENDDO
C
        ITUV = ITUV+1
        DO IT=0,ITUV
          RIT = DFLOAT(IT)
        DO IU=0,ITUV-IT
          RIU = DFLOAT(IU)
        DO IV=0,ITUV-IT-IU
          RIV = DFLOAT(IV)
C
          N1 = INABCD(IT+1,IU  ,IV  )
          N2 = INABCD(IT  ,IU+1,IV  )
          N3 = INABCD(IT  ,IU  ,IV+1)
C
          IF(IT.NE.0) THEN
            IF(IU.NE.0) THEN
              IF(IV.NE.0) THEN      
C               CASE (1 1 1)
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
              ELSE
C               CASE (1 1 0)
                K1 = INABCD(IT  ,IU  ,IV)
                M1 = INABCD(IT-1,IU  ,IV)
                M2 = INABCD(IT  ,IU-1,IV)
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1)
                ENDDO
              ENDIF
            ELSE
              IF(IV.NE.0) THEN
C               CASE (1 0 1)
                K1 = INABCD(IT  ,IU,IV  )
                M1 = INABCD(IT-1,IU,IV  )
                M3 = INABCD(IT  ,IU,IV-1)
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
              ELSE
C               CASE (1 0 0)
                K1 = INABCD(IT  ,IU,IV)
                M1 = INABCD(IT-1,IU,IV)
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1)
                ENDDO
              ENDIF
            ENDIF
          ELSE
            IF(IU.NE.0) THEN
              IF(IV.NE.0) THEN
C               CASE (0 1 1)
                K1=INABCD(IT,IU  ,IV  )
                M2=INABCD(IT,IU-1,IV  )
                M3=INABCD(IT,IU  ,IV-1)
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
              ELSE
C               CASE (0 1 0)
                K1 = INABCD(IT,IU  ,IV)
                M2 = INABCD(IT,IU-1,IV)
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1)
                ENDDO
              ENDIF
            ELSE
              IF(IV.NE.0) THEN
C               CASE (0 0 1)
                K1 = INABCD(IT,IU,IV  )
                M3 = INABCD(IT,IU,IV-1)
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
              ELSE
C               CASE (0 0 0)
                K1 = INABCD(IT,IU,IV)
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1)
                ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        ENDDO
        ENDDO
C       ADD IN (TI=0,IU=0,IV=0) CASE
        DO M=1,MAXM
          RC(M,1) = ((-2.0D0*APH(M))**(ILEVEL-1))*FS(M,ILEVEL)
        ENDDO
      ENDDO
C
      IF(MOD(LAMBDA,2).EQ.1) THEN
C
        ITUV = ITUV+1
        DO IT=0,ITUV
          RIT = DFLOAT(IT)
        DO IU=0,ITUV-IT
          RIU = DFLOAT(IU)
        DO IV=0,ITUV-IT-IU
          RIV = DFLOAT(IV)
C
          N1 = INABCD(IT+1,IU,IV)
          N2 = INABCD(IT,IU+1,IV)
          N3 = INABCD(IT,IU,IV+1)
C
          IF(IT.NE.0) THEN
            IF(IU.NE.0) THEN
              IF(IV.NE.0) THEN      
C               CASE (1 1 1)
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
              ELSE
C               CASE (1 1 0)
                K1 = INABCD(IT  ,IU  ,IV)
                M1 = INABCD(IT-1,IU  ,IV)
                M2 = INABCD(IT  ,IU-1,IV)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
                ENDDO
              ENDIF
            ELSE
              IF(IV.NE.0) THEN
C               CASE (1 0 1)
                K1 = INABCD(IT,IU,IV)
                M1 = INABCD(IT-1,IU,IV)
                M3 = INABCD(IT,IU,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
              ELSE
C               CASE (1 0 0)
                K1 = INABCD(IT  ,IU,IV)
                M1 = INABCD(IT-1,IU,IV)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
                ENDDO
              ENDIF
            ENDIF
          ELSE
            IF(IU.NE.0) THEN
              IF(IV.NE.0) THEN
C               CASE (0 1 1)
                K1 = INABCD(IT,IU  ,IV  )
                M2 = INABCD(IT,IU-1,IV  )
                M3 = INABCD(IT,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
              ELSE
C               CASE (0 1 0)
                K1 = INABCD(IT,IU  ,IV)
                M2 = INABCD(IT,IU-1,IV)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
                ENDDO
              ENDIF
            ELSE
              IF(IV.NE.0) THEN
C               CASE (0 0 1)
                K1 = INABCD(IT,IU,IV  )
                M3 = INABCD(IT,IU,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
              ELSE
C               CASE (0 0 0)
                K1=INABCD(IT,IU,IV)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
                ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        ENDDO
        ENDDO
C       ADD IN (TI=0,IU=0,IV=0) CASE
        DO M=1,MAXM
C         RC2(M,1)=((-2.0D0*APH(M))**(ILEVEL))*FS(M,ILEVEL+1)
          RC2(M,1) = FS(M,1)
        ENDDO
C       WRITE ARRAY RC2 INTO RC
        ITMAX = ((LAMBDA+1)*(LAMBDA+2)*(LAMBDA+3))/6
        DO IT=1,ITMAX
        DO M=1,MAXM
          RC(M,IT) = RC2(M,IT)
        ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE FUNFX(FX,X,N,LAMBDA,ITYPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C           FFFFFFFF UU    UU NN    NN FFFFFFFF MM       MM           C
C           FF       UU    UU NNN   NN FF       MMM     MMM           C
C           FF       UU    UU NNNN  NN FF       MMMM   MMMM           C
C           FFFFFF   UU    UU NN NN NN FFFFFF   MM MM MM MM           C
C           FF       UU    UU NN  NNNN FF       MM  MMM  MM           C
C           FF       UU    UU NN   NNN FF       MM   M   MM           C
C           FF        UUUUUU  NN    NN FF       MM       MM           C
C                                                                     C
C ------------------------------------------------------------------- C
C     FUNFX EVALUATES INTEGRAL [INT_{0}^{1} U^{2M} EXP(-TU^{2}) DU]   C
C     FOR VARIABLE X > 0 FOR ALL ORDERS 0 < M < LAMBDA.               C
C ------------------------------------------------------------------- C
C     ITYPE = 0  SPECIAL CASE X = 0.0D0                               C
C     ITYPE = 1  POWER SERIES AND REVERSE RECURRENCE                  C
C                (ONLY MSER TERMS WILL BE USED, SO USE MUST SUPPLY    C
C                A VALUE APPROPRIATE TO THE MAX VALUE OF X IN BATCH). C
C     ITYPE = 2  ASYMPTOTIC EXPANSION AND FORWARD RECURRENCE.         C        
C     ITYPE = 3  ASYMPTOTIC EXPANSION AND FORWARD RECURRENCE.         C
C                ALL TERMS DEPENDING ON EXP(-X) ARE OMITTED TO AVOID  C
C                NUMERICAL UNDERFLOW PROBLEMS. MSER NOT REQUIRED.     C
C*********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MLM=30,MSER=60,NINT=1801,ISTEP=100)

      DIMENSION FX(MB2,MLM),X(MB2),XLAMBDA(MB2),XX2(MB2),
     &          XEXP(MB2),XROOT(MB2),HSTEP(MB2),XSTEP(MB2),JM(MB2)
C
      COMMON/FNFM/FM(NINT,MLM),T(NINT)
C     
      DATA PIROOT,A0,B0/8.862269254527580D-1,4.994501191201870D-1,
     &                  4.551838436668326D-1/
C
C*********************************************************************C
C     ITYPE = 0: SPECIAL CASE FOR T = 0.0D0                           C
C*********************************************************************C
C
      IF(ITYPE.EQ.0) THEN
        DO JJ=1,LAMBDA+1
          VALUE = 1.0D0/DFLOAT(2*JJ-1)
          DO M=1,N
            FX(M,JJ) = VALUE
          ENDDO
        ENDDO
        RETURN
C
C*********************************************************************C
C     ITYPE = 1: POWER SERIES EVALUATION (INITIALIZE AT M = LAMBDA)   C
C*********************************************************************C
C
      ELSEIF(ITYPE.EQ.1) THEN
C        DO M=1,N
C          TEXP(M)      = DEXP(-T(M))
C          TT2(M)       = 2.0D0*T(M)
C          TLAMBDA(M)   = 1.0D0
C          FM(M,LAMBDA+1) = 1.0D0
C        ENDDO
CC
CC       LOOP OVER TERMS IN THE POWER SERIES
CC       (CONVERGENCE IS ACHIEVED WHEN EXP(-T)*TERM(JJ)< 1.0D-14 WHERE
CC       TERM(JJ) IS THE JJTH TERM IN THE POWER SERIES)
C        DO JJ=1,MSER
C          DLAMBDA = DFLOAT(2*(LAMBDA+JJ)+1)
C          DO M=1,N
C            TLAMBDA(M)     = TLAMBDA(M)*(TT2(M)/DLAMBDA)
C            FM(M,LAMBDA+1) = FM(M,LAMBDA+1)+TLAMBDA(M)
C          ENDDO
C        ENDDO
        DO M=1,N
          IM             = (DFLOAT(ISTEP)*X(M)) + 0.5D0
          IM             = IM + 1
          JM(M)          = IM
          FX(M,LAMBDA+1) = FM(IM,LAMBDA+1)
          HSTEP(M)       = X(M) - T(IM)
          XSTEP(M)       =-HSTEP(M)
          XEXP(M)        = DEXP(-X(M))
          XX2(M)         = 2.0D0*X(M)
        ENDDO

        DO ITERM=1,3
          DO M=1,N
            FX(M,LAMBDA+1) = FX(M,LAMBDA+1) 
     &                       + XSTEP(M)*FM(JM(M),LAMBDA+ITERM+1)
            XSTEP(M) = -XSTEP(M)*HSTEP(M)/DFLOAT(ITERM+1)
          ENDDO
        ENDDO
C
C       NOW COMPLETE TABLE BY BACKWARDS RECURRENCE
        DO I=1,LAMBDA
          MIND  = LAMBDA-I+1
          COEFF = DFLOAT(2*MIND-1)
          DO M=1,N
            FX(M,MIND) = (XX2(M)*FX(M,MIND+1) + XEXP(M))/COEFF
          ENDDO
        ENDDO
        RETURN
C
C*********************************************************************C
C     ITYPE = 2: ASYMPTOTIC EXPANSION WITH VERY LARGE ARGUMENT        C
C*********************************************************************C
C
      ELSEIF(ITYPE.EQ.2) THEN
C
C       INITIALIZE THE ASYMPTOTIC EXPANSION
        DO M=1,N
          XEXP(M)  = DEXP(-X(M))
          XX2(M)   = 2.0D0*X(M)
          XROOT(M) = DSQRT(X(M)) 
        ENDDO
        DO M=1,N
          FX(M,1) = A0/(B0+X(M))
        ENDDO
C
C       RESCALE BY THE PREFACTOR
        DO M=1,N
          FX(M,1) = (PIROOT/XROOT(M))-(XEXP(M)*FX(M,1))
        ENDDO
C
C       NOW COMPLETE TABLE BY FORWARD RECURRENCE
        DO MIND=1,LAMBDA
          COEFF = DFLOAT(2*MIND-1)
          DO M=1,N
            FX(M,MIND+1) = (COEFF*FX(M,MIND)-XEXP(M))/XX2(M)
          ENDDO
        ENDDO
        RETURN
C
C*********************************************************************C
C     ITYPE = 3: ASYMPTOTIC EXPANSION WITH VERY LARGE ARGUMENT        C
C*********************************************************************C
C
      ELSEIF(ITYPE.EQ.3) THEN
C
C       INITIALIZE THE ASYMPTOTIC EXPANSION
        DO M=1,N
          XX2(M)  = 2.0D0*X(M)
          FX(M,1) = PIROOT/DSQRT(X(M))
        ENDDO
C
C       NOW COMPLETE TABLE BY FORWARD RECURRENCE
        DO MIND=1,LAMBDA
          COEFF = DFLOAT(2*MIND-1)
          DO M=1,N
            FX(M,MIND+1) = (COEFF*FX(M,MIND))/XX2(M)
          ENDDO
        ENDDO
C
C*********************************************************************C
C     ITYPE OUT OF RANGE: INVALID INPUT TO FUNFX                      C
C*********************************************************************C
C
      ELSE
91      FORMAT(2X,'In FUNFX: invalid type (must be 0-3)',I4)
        WRITE(6,91) ITYPE
        WRITE(7,91) ITYPE
        STOP
      ENDIF
C
      RETURN
      END

