      SUBROUTINE RMAKE(RC,QP,APH,MAXM,LAM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           RRRRRRR  MM       MM    AA    KK    KK EEEEEEEE            C
C           RR    RR MMM     MMM   AAAA   KK   KK  EE                  C
C           RR    RR MMMM   MMMM  AA  AA  KK  KK   EE                  C
C           RR    RR MM MM MM MM AA    AA KKKKK    EEEEEE              C
C           RRRRRRR  MM  MMM  MM AAAAAAAA KK  KK   EE                  C
C           RR    RR MM   M   MM AA    AA KK   KK  EE                  C
C           RR    RR MM       MM AA    AA KK    KK EEEEEEEE            C
C                                                                      C
C -------------------------------------------------------------------- C
C  RMAKE GENERATES A COMPLETE SET OF R-INTEGRALS REQUIRED IN THE       C
C  FINITE SUM REPRESENTATION OF A MULTI-CENTER GAUSSIAN OVERLAP.       C
C**********************************************************************C
      PARAMETER(MDM=2500,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26,MB2=MBS*MBS,
     &          IL4=2*(MKP-1),MRC=(IL4+1)*(IL4+2)*(IL4+3)/6,MLM=30)
C     
      DIMENSION FS(MB2,MLM),APH(MB2),QP(MB2,3),RC(MB2,MRC),RC2(MB2,MRC)
      DIMENSION F1(MB2,MLM),F2(MB2,MLM),F3(MB2,MLM),F4(MB2,MLM)     
      DIMENSION X1(MB2),X2(MB2),X3(MB2),X4(MB2),
     &          I1(MB2),I2(MB2),I3(MB2),I4(MB2)
C
      COMMON/ACSS/INABCD(0:IL4,0:IL4,0:IL4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
C
C**********************************************************************C
C     THE FIRST STEP OF THIS ROUTINE IS TO EVALUATE THE REQUIRED       C
C     BOYS INTEGRALS. THIS IS DIVIDED INTO CASES, DEPENDING ON THE     C
C     MAGNITUDE OF THE ARGUMENT.                                       C
C -------------------------------------------------------------------- C
C        FS_M (X) = INT_{0}^{1} T^{2M} EXP(-X*T^{2}) DT                C
C -------------------------------------------------------------------- C
C   EVALUATED FOR ALL VALUES OF M IN THE RANGE 0 < M < LAM.            C
C             FOR ALL VALUES OF X IN THE RANGE X > 0.                  C
C**********************************************************************C
C
      N1 = 0
      N2 = 0
      N3 = 0
      N4 = 0
C
C     FOR EACH PAIR OF BASIS FUNCTIONS (EXPONENTS EI AND EJ IN 'M'),
C     DETERMINE THE BEST WAY TO EVALUATE THE BOYS FUNCTION
      DO M=1,MAXM
        X = APH(M)*(QP(M,1)**2 + QP(M,2)**2 + QP(M,3)**2)
        IF(X.LE.1.0D-11) THEN
C       CASE 1: IF X IS ALMOST ZERO (SO WHEN Q=P OR EIJ<<1)
          N1     = N1 + 1
          X1(N1) = X
          I1(N1) = M
        ELSEIF(X.GT.1.00D-11.AND.X.LE.1.70D+01) THEN
C       CASE 2: IF X IS SMALLER THAN X=17.0D0
          N2     = N2 + 1
          X2(N2) = X
          I2(N2) = M
        ELSEIF(X.GT.1.70D+01.AND.X.LE.3.00D+01) THEN
C       CASE 3: IF X IS SMALLER THAN ABOUT 30.0D0
          N3     = N3 + 1
          X3(N3) = X
          I3(N3) = M
        ELSE
C       CASE 4: IF X IS LARGER THAN ABOUT 30.0D0
          N4     = N4 + 1
          X4(N4) = X
          I4(N4) = M
        ENDIF
      ENDDO
C
C     CASE 1: ARGUMENT OF THE BOYS FUNCTION IS X=0.0D0.
C     THE VALUE OF THIS FUNCTION IS 2N+1 (DONE IN FUNFM).
      IF(N1.NE.0) THEN
        CALL FUNFM(F1,X1,N1,LAM,1)     
        DO JJ=1,LAM+1
          DO M=1,N1
            FS(I1(M),JJ) = F1(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 2: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=17.0D0.
C     EVALUATE WITH LOCAL POLYNOMIAL EXPANSION OF ORDER 5,
C     AND RECURRENCE IN DIRECTION OF DECREASING M.
      IF(N2.NE.0) THEN
        CALL FUNFM(F2,X2,N2,LAM,2)
        DO JJ=1,LAM+1
          DO M=1,N2
            FS(I2(M),JJ) = F2(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 3: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=30.0D0.
C     EVALUATE USING ASYMPTOTIC FORMULA WITH EXPONENTIAL.
      IF(N3.NE.0) THEN
        CALL FUNFM(F3,X3,N3,LAM,3)
        DO JJ=1,LAM+1
          DO M=1,N3
            FS(I3(M),JJ) = F3(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 4: ARGUMENT OF THE BOYS FUNCTION IS LARGER THAN X=30.0D0.
C     EVALUATE USING ASYMPTOTIC FORMULA WITHOUT EXPONENTIAL.
      IF(N4.NE.0) THEN
        CALL FUNFM(F4,X4,N4,LAM,4)
        DO JJ=1,LAM+1
          DO M=1,N4
            FS(I4(M),JJ) = F4(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C**********************************************************************C
C     THE SECOND STEP OF THIS ROUTINE IS TO EVALUATE THE ACTUAL        C
C     R-COEFFICIENTS, BASED ON THE CORRESPONDING BOYS FUNCTIONS.       C
C**********************************************************************C
C
C     CONSTRUCT TOP-LEVEL
      DO M=1,MAXM
        RC(M,1) = ((-2.0D0*APH(M))**LAM)*FS(M,LAM+1)
      ENDDO
C
      IF(MOD(LAM,2).EQ.0) THEN
        ITUVMIN = 1
      ELSE
        ITUVMIN = 2
      ENDIF
C
      ITUV = -1
      DO ILEVEL=LAM-1,ITUVMIN,-2
        ITUV = ITUV+1
        DO IT=0,ITUV
          RIT = DFLOAT(IT)
          DO IU=0,ITUV-IT
            RIU = DFLOAT(IU)
            DO IV=0,ITUV-IT-IU
              RIV = DFLOAT(IV)
C
              N2 = INABCD(IT+1,IU  ,IV  )
              N3 = INABCD(IT  ,IU+1,IV  )
              N4 = INABCD(IT  ,IU  ,IV+1)
C
              IF(IT.NE.0) THEN
                IF(IU.NE.0) THEN
                  IF(IV.NE.0) THEN      
C                   CASE (1 1 1)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M1 = INABCD(IT-1,IU  ,IV  )
                    M2 = INABCD(IT  ,IU-1,IV  )
                    M3 = INABCD(IT  ,IU  ,IV-1)
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                    ENDDO
                  ELSE
C                   CASE (1 1 0)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M1 = INABCD(IT-1,IU  ,IV  )
                    M2 = INABCD(IT  ,IU-1,IV  )
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1)
                    ENDDO
                  ENDIF
                ELSE
                  IF(IV.NE.0) THEN
C                   CASE (1 0 1)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M1 = INABCD(IT-1,IU  ,IV  )
                    M3 = INABCD(IT  ,IU  ,IV-1)
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                    ENDDO
                  ELSE
C                   CASE (1 0 0)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M1 = INABCD(IT-1,IU  ,IV  )
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1)
                    ENDDO
                  ENDIF
                ENDIF
              ELSE
                IF(IU.NE.0) THEN
                  IF(IV.NE.0) THEN
C                   CASE (0 1 1) 
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M2 = INABCD(IT  ,IU-1,IV  )
                    M3 = INABCD(IT  ,IU  ,IV-1)
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                    ENDDO
                  ELSE
C                   CASE (0 1 0)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M2 = INABCD(IT  ,IU-1,IV  )
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1)
                    ENDDO
                  ENDIF
                ELSE
                  IF(IV.NE.0) THEN
C                   CASE (0 0 1)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M3 = INABCD(IT  ,IU  ,IV-1)
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                    ENDDO
                  ELSE
C                   CASE (0 0 0)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1)
                    ENDDO
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C       ADD IN (IT=0,IU=0,IV=0) CASE
        DO M=1,MAXM
          RC2(M,1) = ((-2.0D0*APH(M))**ILEVEL)*FS(M,ILEVEL+1)
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
              N2 = INABCD(IT+1,IU  ,IV  )
              N3 = INABCD(IT  ,IU+1,IV  )
              N4 = INABCD(IT  ,IU  ,IV+1)
C
              IF(IT.NE.0) THEN
                IF(IU.NE.0) THEN
                  IF(IV.NE.0) THEN      
C                   CASE (1 1 1)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M1 = INABCD(IT-1,IU  ,IV  )
                    M2 = INABCD(IT  ,IU-1,IV  )
                    M3 = INABCD(IT  ,IU  ,IV-1)
                    DO M=1,MAXM
                      RC(M,N2) = -QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                      RC(M,N3) = -QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                      RC(M,N4) = -QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                    ENDDO
                  ELSE
C                   CASE (1 1 0)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M1 = INABCD(IT-1,IU  ,IV  )
                    M2 = INABCD(IT  ,IU-1,IV  )
                    DO M=1,MAXM
                      RC(M,N2) = -QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                      RC(M,N3) = -QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                      RC(M,N4) = -QP(M,3)*RC2(M,K1)
                    ENDDO
                  ENDIF
                ELSE
                  IF(IV.NE.0) THEN
C                   CASE (1 0 1)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M1 = INABCD(IT-1,IU  ,IV  )
                    M3 = INABCD(IT  ,IU  ,IV-1)
                    DO M=1,MAXM
                      RC(M,N2) = -QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                      RC(M,N3) = -QP(M,2)*RC2(M,K1)
                      RC(M,N4) = -QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                    ENDDO
                  ELSE
C                   CASE (1 0 0)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M1 = INABCD(IT-1,IU  ,IV  )
                    DO M=1,MAXM
                      RC(M,N2) = -QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                      RC(M,N3) = -QP(M,2)*RC2(M,K1)
                      RC(M,N4) = -QP(M,3)*RC2(M,K1)
                    ENDDO
                  ENDIF
                ENDIF
              ELSE
                IF(IU.NE.0) THEN
                  IF(IV.NE.0) THEN
C                   CASE (0 1 1)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M2 = INABCD(IT  ,IU-1,IV  )
                    M3 = INABCD(IT  ,IU  ,IV-1)
                    DO M=1,MAXM
                      RC(M,N2) = -QP(M,1)*RC2(M,K1)
                      RC(M,N3) = -QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                      RC(M,N4) = -QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                    ENDDO
                  ELSE
C                   CASE (0 1 0)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M2 = INABCD(IT  ,IU-1,IV  )
                    DO M=1,MAXM
                      RC(M,N2) = -QP(M,1)*RC2(M,K1)
                      RC(M,N3) = -QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                      RC(M,N4) = -QP(M,3)*RC2(M,K1)
                    ENDDO
                  ENDIF
                ELSE
                  IF(IV.NE.0) THEN
C                   CASE (0 0 1)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M3 = INABCD(IT  ,IU  ,IV-1)
                    DO M=1,MAXM
                      RC(M,N2) = -QP(M,1)*RC2(M,K1)
                      RC(M,N3) = -QP(M,2)*RC2(M,K1)
                      RC(M,N4) = -QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                    ENDDO
                  ELSE
C                   CASE (0 0 0)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    DO M=1,MAXM
                      RC(M,N2) = -QP(M,1)*RC2(M,K1)
                      RC(M,N3) = -QP(M,2)*RC2(M,K1)
                      RC(M,N4) = -QP(M,3)*RC2(M,K1)
                    ENDDO
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
C       ADD IN THE (IT=0,IU=0,IV=0) CASE
        DO M=1,MAXM
          RC(M,1) = ((-2.0D0*APH(M))**(ILEVEL-1))*FS(M,ILEVEL)
        ENDDO
      ENDDO
C
      IF(MOD(LAM,2).EQ.1) THEN
C
        ITUV = ITUV+1
        DO IT=0,ITUV
          RIT = DFLOAT(IT)
          DO IU=0,ITUV-IT
            RIU = DFLOAT(IU)
            DO IV=0,ITUV-IT-IU
              RIV = DFLOAT(IV)
C
              N2 = INABCD(IT+1,IU  ,IV  )
              N3 = INABCD(IT  ,IU+1,IV  )
              N4 = INABCD(IT  ,IU  ,IV+1)
C
              IF(IT.NE.0) THEN
                IF(IU.NE.0) THEN
                  IF(IV.NE.0) THEN      
C                   CASE (1 1 1)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M1 = INABCD(IT-1,IU  ,IV  )
                    M2 = INABCD(IT  ,IU-1,IV  )
                    M3 = INABCD(IT  ,IU  ,IV-1)
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                    ENDDO
                  ELSE
C                   CASE (1 1 0)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M1 = INABCD(IT-1,IU  ,IV  )
                    M2 = INABCD(IT  ,IU-1,IV  )
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1)
                    ENDDO
                  ENDIF
                ELSE
                  IF(IV.NE.0) THEN
C                   CASE (1 0 1)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M1 = INABCD(IT-1,IU  ,IV  )
                    M3 = INABCD(IT  ,IU  ,IV-1)
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                    ENDDO
                  ELSE
C                   CASE (1 0 0)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M1 = INABCD(IT-1,IU  ,IV  )
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1)
                    ENDDO
                  ENDIF
                ENDIF
              ELSE
                IF(IU.NE.0) THEN
                  IF(IV.NE.0) THEN
C                   CASE (0 1 1)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M2 = INABCD(IT  ,IU-1,IV  )
                    M3 = INABCD(IT  ,IU  ,IV-1)
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                    ENDDO
                  ELSE
C                   CASE (0 1 0)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M2 = INABCD(IT  ,IU-1,IV  )
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1)
                    ENDDO
                  ENDIF
                ELSE
                  IF(IV.NE.0) THEN
C                   CASE (0 0 1)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    M3 = INABCD(IT  ,IU  ,IV-1)
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                    ENDDO
                  ELSE
C                   CASE (0 0 0)
                    K1 = INABCD(IT  ,IU  ,IV  )
                    DO M=1,MAXM
                      RC2(M,N2) = -QP(M,1)*RC(M,K1)
                      RC2(M,N3) = -QP(M,2)*RC(M,K1)
                      RC2(M,N4) = -QP(M,3)*RC(M,K1)
                    ENDDO
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
C       ADD IN (IT=0,IU=0,IV=0) CASE
        DO M=1,MAXM
          RC2(M,1)=((-2.0D0*APH(M))**(ILEVEL))*FS(M,ILEVEL+1)
        ENDDO
C
C       WRITE ARRAY RC2 INTO RC
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
        DO ITUV=1,NTUV
          DO M=1,MAXM
            RC(M,ITUV) = RC2(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE WMAKE(WC,RP,APH,MAXM,LAM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         WW         WW MM       MM    AA    KK    KK EEEEEEEE         C
C         WW         WW MMM     MMM   AAAA   KK   KK  EE               C
C         WW         WW MMMM   MMMM  AA  AA  KK  KK   EE               C
C         WW    W    WW MM MM MM MM AA    AA KKKKK    EEEEEE           C
C          WW  WWW  WW  MM  MMM  MM AAAAAAAA KK  KK   EE               C
C           WWWW WWWW   MM   M   MM AA    AA KK   KK  EE               C
C            WW   WW    MM       MM AA    AA KK    KK EEEEEEEE         C
C                                                                      C
C -------------------------------------------------------------------- C
C  FOR A PARTICULAR BLOCK OF BASIS FUNCTION OVERLAPS, MAKE THE W       C
C  FUNCTIONS W(APH,RP;A,B,C) AS A FUNCTION OF (X,Y,Z) -- THIS          C
C  COMES TO THE POTENTIAL FROM A HGTF OVERLAP CHARGE SOURCE.           C
C**********************************************************************C
      PARAMETER(MDM=2500,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26,MB2=MBS*MBS,
     &          IL4=2*(MKP-1),MRC=(IL4+1)*(IL4+2)*(IL4+3)/6,MLM=30)
C
      DIMENSION FS(MB2,MLM),APH(MB2),RP(MB2,3),WC(MB2,MRC)
      DIMENSION X1(MB2),X2(MB2),X3(MB2),X4(MB2),
     &          I1(MB2),I2(MB2),I3(MB2),I4(MB2)
      DIMENSION F1(MB2,MLM),F2(MB2,MLM),F3(MB2,MLM),F4(MB2,MLM)     
      DIMENSION H(0:MKP-1,0:MMV-1,3)
C
      COMMON/ACCSS/INABCD(0:IL4,0:IL4,0:IL4),
     &             IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
C
      DATA PI/3.141592653589793D0/      
C
C**********************************************************************C
C     THE FIRST STEP OF THIS ROUTINE IS TO EVALUATE THE REQUIRED       C
C     BOYS INTEGRALS. THIS IS DIVIDED INTO CASES, DEPENDING ON THE     C
C     MAGNITUDE OF THE ARGUMENT.                                       C
C -------------------------------------------------------------------- C
C        FS_M (X) = INT_{0}^{1} T^{2M} EXP(-X*T^{2}) DT                C
C -------------------------------------------------------------------- C
C   EVALUATED FOR ALL VALUES OF M IN THE RANGE 0 < M < LAM.            C
C             FOR ALL VALUES OF X IN THE RANGE X > 0.                  C
C**********************************************************************C
C
C     CALCULATE GAMMA FUNCTION VALUES FOR LATER USE
      CALL GAMMAS
C
      N1 = 0
      N2 = 0
      N3 = 0
      N4 = 0
C
C     FOR EACH PAIR OF BASIS FUNCTIONS (EXPONENTS EI AND EJ IN 'M'),
C     DETERMINE THE BEST WAY TO EVALUATE THE BOYS FUNCTION
      DO M=1,MAXM
C       X = EIJ*(P-C)^2
        X = APH(M)*(RP(M,1)*RP(M,1)+RP(M,2)*RP(M,2)+RP(M,3)*RP(M,3))
C       CASE 1: IF X IS ALMOST ZERO (SO WHEN Q=P OR EIJ<<1)
        IF(X.LE.1.00D-11) THEN
          N1     = N1 + 1
          X1(N1) = X
          I1(N1) = M1
C       CASE 2: IF X IS SMALLER THAN X=17.0D0
        ELSEIF(X.GT.1.00D-11.AND.X.LE.1.70D+01) THEN
          N2     = N2 + 1
          X2(N2) = X
          I2(N2) = M
C       CASE 3: IF X IS SMALLER THAN ABOUT 30.0D0
        ELSEIF(X.GT.1.70D+01.AND.X.LE.3.00D+01) THEN
          N3     = N3 + 1
          X3(N3) = X
          I3(N3) = M
C       CASE 4: IF X IS LARGER THAN ABOUT 30.0D0
        ELSE
          N4     = N4 + 1
          X4(N4) = X
          I4(N4) = M
        ENDIF
      ENDDO
C
C     CASE 1: ARGUMENT OF THE BOYS FUNCTION IS X=0.0D0.
C     THE VALUE OF THIS FUNCTION IS 2N+1 (DONE IN FUNFM).
      IF(N1.NE.0) THEN
        CALL FUNFM(F1,X1,N1,LAM,1)     
        DO JJ=1,LAM+1
          DO M=1,N1
            FS(I1(M),JJ) = F1(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 2: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=17.0D0.
C     EVALUATE WITH LOCAL POLYNOMIAL EXPANSION OF ORDER 5,
C     AND RECURRENCE IN DIRECTION OF DECREASING M.
      IF(N2.NE.0) THEN
        CALL FUNFM(F2,X2,N2,LAM,2)
        DO JJ=1,LAM+1
          DO M=1,N2
            FS(I2(M),JJ) = F2(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 3: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=30.0D0.
C     EVALUATE USING ASYMPTOTIC FORMULA WITH EXPONENTIAL.
      IF(N3.NE.0) THEN
        CALL FUNFM(F3,X3,N3,LAM,3)
        DO JJ=1,LAM+1
          DO M=1,N3
            FS(I3(M),JJ) = F3(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 4: ARGUMENT OF THE BOYS FUNCTION IS LARGER THAN X=30.0D0.
C     EVALUATE USING ASYMPTOTIC FORMULA WITHOUT EXPONENTIAL.
      IF(N4.NE.0) THEN
        CALL FUNFM(F4,X4,N4,LAM,4)
        DO JJ=1,LAM+1
          DO M=1,N4
            FS(I4(M),JJ) = F4(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C**********************************************************************C
C     THE SECOND STEP OF THIS SUBROUTINE IS TO EVALUATE THE ACTUAL     C
C     W-FUNCTIONS, BASED ON THE CORRESPONDING BOYS INTEGRALS.          C
C**********************************************************************C
C
C     INITIALISE THE WC ARRAY
      DO M=1,MAXM
        DO ITUV=1,LAM
          WC(M,ITUV) = 0.0D0
        ENDDO
      ENDDO
C
C     NUMBER OF ADDRESSES IADD FOR LAM VALUE
      NTUV = ((LAM+1)*(LAM+2)*(LAM+3))/6
C
C     LOOP OVER BASIS FUNCTION PAIRS
      DO M=1,MAXM
C       
C       GENERATE ALL NECESSARY HERMITE POLYNOMIALS FOR LATER REFERENCE
        DO ILAM=0,LAM
C         TOTAL NUMBER OF TERMS IN HERMITE EXPANSION
          NTERMS = (ILAM-MOD(ILAM,2))/2
          DO N=0,NTERMS
            H(ILAM,N,1) = HCOEFF(APH(M),RP(M,1),ILAM,N)
            H(ILAM,N,2) = HCOEFF(APH(M),RP(M,2),ILAM,N)
            H(ILAM,N,3) = HCOEFF(APH(M),RP(M,3),ILAM,N)
          ENDDO
        ENDDO
C
C       START A LOOP THAT RUNS OVER ALL COMBINATIONS ITUV={T,U,V}
        DO ITUV=1,NTUV
C
C         ESTABLISH PARAMETERS T,U,V AND LAM FOR THIS VALUE ITUV
          IT   = IVEC(ITUV)
          IU   = JVEC(ITUV)
          IV   = KVEC(ITUV)
          ILAM = LAMVEC(ITUV)
C
C         FOR GIVEN ITUV, START THREE NESTED LOOPS FOR X,Y,Z DECOMP.
          STORE = 0.0D0
          DO IA=0,(IT-MOD(IT,2))/2
            IPOW = (IT + MOD(IT,2))/2 + IA
            DO IB=0,(IU-MOD(IU,2))/2
              JPOW = (IU + MOD(IU,2))/2 + IB
              DO IC=0,(IV-MOD(IV,2))/2
                KPOW = (IV + MOD(IV,2))/2 + IC
C               
C               PRODUCT OF HERMITE EXPANSION COEFFICIENTS
                HTMS = H(IT,IA,1)*H(IU,IB,2)*H(IV,IC,3)
C
C               BOYS FUNCTION OF RELEVANT INDEX NBOYS
                NBOYS = IPOW + JPOW + KPOW
                FBOYS = FS(M,NBOYS+1)
C               ADD TO W FUNCTION BIN THE RELEVANT PRODUCT
                STORE = STORE + HTMS*FBOYS
C
              ENDDO
            ENDDO
          ENDDO
C
C         MULTIPLY RESULT BY 2*PI/LAM
          WC(M,ITUV) = 2.0D0*PI*STORE/APH(M)
C
        ENDDO
C     
      ENDDO
C      
      RETURN
      END
C
C
      SUBROUTINE FUNFM(FX,X,N,LAM,ITP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           FFFFFFFF UU    UU NN    NN FFFFFFFF MM       MM            C
C           FF       UU    UU NNN   NN FF       MMM     MMM            C
C           FF       UU    UU NNNN  NN FF       MMMM   MMMM            C
C           FFFFFF   UU    UU NN NN NN FFFFFF   MM MM MM MM            C
C           FF       UU    UU NN  NNNN FF       MM  MMM  MM            C
C           FF       UU    UU NN   NNN FF       MM   M   MM            C
C           FF        UUUUUU  NN    NN FF       MM       MM            C
C                                                                      C
C -------------------------------------------------------------------- C
C  FUNFM EVALUATES INTEGRAL [INT_{0}^{1} U^{2M} EXP(-T*U^{2}) dU]      C
C  FOR VARIABLE X > 0 FOR ALL ORDERS 0 < M < LAM.                      C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    ITP = 0 - SPECIAL CASE X = 0.0D0.                                 C
C    ITP = 1 - POWER SERIES AND REVERSE RECURRENCE.                    C
C              (ONLY MSER TERMS WILL BE USED, SO USE MUST SUPPLY A     C
C               VALUE APPROPRIATE TO THE MAX VALUE OF X IN BATCH).     C
C    ITP = 2 - ASYMPTOTIC EXPANSION AND FORWARD RECURRENCE.            C
C    ITP = 3 - ASYMPTOTIC EXPANSION AND FORWARD RECURRENCE.            C
C              ALL TERMS DEPENDING ON EXP(-X) ARE OMITTED TO AVOID     C
C              NUMERICAL UNDERFLOW PROBLEMS. MSER NOT REQUIRED.        C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MLM=30,MSER=60,NINT=1801,ISTEP=100)

      DIMENSION FX(MB2,MLM),X(MB2),XLAM(MB2),XX2(MB2),
     &          XEXP(MB2),XROOT(MB2),HSTEP(MB2),XSTEP(MB2),JM(MB2)
C
      COMMON/FNFM/FM(NINT,MLM),T(NINT)
C     
      DATA PIROOT,A0,B0/8.862269254527580D-1, 4.994501191201870D-1,
     &                                        4.551838436668326D-1/
C
C     ITP = 1: SPECIAL CASE FOR T = 0.0D0
      IF(ITP.EQ.1) THEN
        DO JJ=1,LAM+1
          DEN   = DFLOAT(2*JJ-1)
          VALUE = 1.0D0/DEN
          DO M=1,N
            FX(M,JJ) = VALUE
          ENDDO
        ENDDO
        RETURN
C
C     ITP = 2: POWER SERIES EVALUATION (INITIALIZE AT M = LAM)
      ELSEIF(ITP.EQ.2) THEN
        DO M=1,N
          IM          = (DFLOAT(ISTEP)*X(M)) + 0.5D0
          IM          = IM + 1
          JM(M)       = IM
          FX(M,LAM+1) = FM(IM,LAM+1)
          HSTEP(M)    = X(M) - T(IM)
          XSTEP(M)    =-HSTEP(M)
          XEXP(M)     = DEXP(-X(M))
          XX2(M)      = 2.0D0*X(M)
        ENDDO

        DO ITERM=1,3
          RI1 = DFLOAT(ITERM+1)
          DO M=1,N
            FX(M,LAM+1) = FX(M,LAM+1) + XSTEP(M)*FM(JM(M),LAM+ITERM+1)
            XSTEP(M)    =-XSTEP(M)*HSTEP(M)/RI1
          ENDDO
        ENDDO
C
C       NOW COMPLETE TABLE BY BACKWARDS RECURRENCE
        DO I=1,LAM
          MIND  = LAM-I+1
          COEFF = DFLOAT(2*MIND-1)
          DO M=1,N
            FX(M,MIND) = (XX2(M)*FX(M,MIND+1) + XEXP(M))/COEFF
          ENDDO
        ENDDO
        RETURN
C
C     ITP = 3: ASYMPTOTIC EXPANSION WITH VERY LARGE ARGUMENT
      ELSEIF(ITP.EQ.3) THEN
C
C       INITIALIZE THE ASYMPTOTIC EXPANSION
        DO M=1,N
          XEXP(M)  = DEXP(-X(M))
          XX2(M)   = 2.0D0*X(M)
          XROOT(M) = DSQRT(X(M)) 
        ENDDO
        DO M=1,N
          B0X     = B0+X(M)
          FX(M,1) = A0/B0X
        ENDDO
C
C       RESCALE BY THE PREFACTOR
        DO M=1,N
          FX(M,1) = (PIROOT/XROOT(M))-(XEXP(M)*FX(M,1))
        ENDDO
C
C       NOW COMPLETE TABLE BY FORWARD RECURRENCE
        DO MIND=1,LAM
          COEFF = DFLOAT(2*MIND-1)
          DO M=1,N
            FX(M,MIND+1) = (COEFF*FX(M,MIND)-XEXP(M))/XX2(M)
          ENDDO
        ENDDO
        RETURN
C
C     ITP = 4: ASYMPTOTIC EXPANSION WITH VERY LARGE ARGUMENT
      ELSEIF(ITP.EQ.4) THEN
C
C       INITIALIZE THE ASYMPTOTIC EXPANSION
        DO M=1,N
          XX2(M)  = 2.0D0*X(M)
          FX(M,1) = PIROOT/DSQRT(X(M))
        ENDDO
C
C       NOW COMPLETE TABLE BY FORWARD RECURRENCE
        DO MIND=1,LAM
          COEFF = DFLOAT(2*MIND-1)
          DO M=1,N
            FX(M,MIND+1) = (COEFF*FX(M,MIND))/XX2(M)
          ENDDO
        ENDDO
C
C     ITP OUT OF RANGE: INVALID INPUT TO FUNFM
      ELSE
20      FORMAT(2X,'In FUNFM: invalid type (must be 1-4)',I4)
        WRITE(6,20) ITP
        WRITE(7,20) ITP
        STOP
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE BOYSGEN(LAMBDA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    BBBBBBB   OOOOOO  YY    YY  SSSSSS   GGGGGG  EEEEEEEE NN    NN    C
C    BB    BB OO    OO YY    YY SS    SS GG    GG EE       NNN   NN    C
C    BB    BB OO    OO  YY  YY  SS       GG       EE       NNNN  NN    C
C    BBBBBBB  OO    OO   YYYY    SSSSSS  GG       EEEEEE   NN NN NN    C
C    BB    BB OO    OO    YY          SS GG   GGG EE       NN  NNNN    C
C    BB    BB OO    OO    YY    SS    SS GG    GG EE       NN   NNN    C
C    BBBBBBB   OOOOOO     YY     SSSSSS   GGGGGG  EEEEEEEE NN    NN    C
C                                                                      C
C -------------------------------------------------------------------- C
C  BOYSGEN PRODUCES A DATA FILE WHICH CONTAINS A FAMILY OF BOYS        C
C  FUNCTIONS OVER A SPECIFIED REGION, GIVEN A MAXIMUM FAMILY           C
C  PARAMETER DETERMINED BY LAMBDA.                                     C
C**********************************************************************C
      PARAMETER (MBS=26,MB2=MBS*MBS,MLM=30)
C
      DIMENSION X1(MB2),X2(MB2),X3(MB2),X4(MB2)
      DIMENSION F1(MB2,MLM),F2(MB2,MLM),F3(MB2,MLM),
     &          F4(MB2,MLM),FS(MB2,MLM)
C
C     EVALUATION PARAMETERS
      NTOT = 400
      XMIN = 0.00D0
      XMAX = 4.00D1
      HSTP = (XMAX-XMIN)/NTOT
C
C**********************************************************************C
C     THE FIRST STEP OF THIS SUBROUTINE IS TO EVALUATE THE REQUIRED    C
C     BOYS INTEGRALS. THIS IS DIVIDED INTO CASES, DEPENDING ON THE     C
C     MAGNITUDE OF THE ARGUMENT.                                       C
C -------------------------------------------------------------------- C
C        FS_M (X) = INT_{0}^{1} T^{2M} EXP(-X*T^{2}) DT                C
C -------------------------------------------------------------------- C
C   EVALUATED FOR ALL VALUES OF M IN THE RANGE 0 < M < LAMBDA.         C
C             FOR ALL VALUES OF X IN THE RANGE X > 0.                  C
C**********************************************************************C
C
      OPEN(UNIT=8,FILE='plots/boysfunction.dat',STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO NX=0,NTOT
C
      X = XMIN + HSTP*NX
C      
C     DETERMINE THE BEST WAY TO EVALUATE THE BOYS FUNCTION
      IF(X.LE.1.00D-11) THEN
C     CASE 0: ARGUMENT OF THE BOYS FUNCTION IS X = 0.
C             THE VALUE OF THIS FUNCTION IS (2N+1).
        N1     = 1
        X1(N1) = X
        CALL FUNFM(F1,X1,N1,LAMBDA,1)
        DO JJ=1,LAMBDA+1
          FS(N1,JJ) = F1(N1,JJ)
        ENDDO
      ELSEIF(X.GT.1.00D-11.AND.X.LE.1.70D+01) THEN
C     CASE 1: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=17.
C             EVALUATE WITH LOCAL POLYNOMIAL EXPANSION OF ORDER 5,
C             AND RECURRENCE IN DIRECTION OF DECREASING M.
        N2     = 1
        X2(N2) = X
        CALL FUNFM(F2,X2,N2,LAMBDA,2)
        DO JJ=1,LAMBDA+1
          FS(N2,JJ) = F2(N2,JJ)
        ENDDO
      ELSEIF(X.GT.1.70D+01.AND.X.LE.3.00D+01) THEN
C     CASE 2: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=30.
C             EVALUATE USING ASYMPTOTIC FORMULA WITH EXPONENTIAL.
        N3     = 1
        X3(N3) = X
        CALL FUNFM(F3,X3,N3,LAMBDA,3)
        DO JJ=1,LAMBDA+1
          FS(N3,JJ) = F3(N3,JJ)
        ENDDO
      ELSE
C     CASE 3: ARGUMENT OF THE BOYS FUNCTION IS LARGER THAN X=30.
C             EVALUATE USING ASYMPTOTIC FORMULA WITHOUT EXPONENTIAL.
        N4     = 1
        X4(N4) = X
        CALL FUNFM(F4,X4,N4,LAMBDA,4)
        DO JJ=1,LAMBDA+1
          FS(N4,JJ) = F4(N4,JJ)
        ENDDO
      ENDIF
C      
      WRITE(8,*) X,(FS(1,L),L=1,LAMBDA+1)
C      
      ENDDO
      CLOSE(UNIT=8)
C      
      RETURN
      END
C
C
      SUBROUTINE GFINIT(MMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             GGGGGG  FFFFFFFF IIII NN    NN IIII TTTTTTTT             C
C            GG    GG FF        II  NNN   NN  II     TT                C
C            GG       FF        II  NNNN  NN  II     TT                C
C            GG       FFFFFF    II  NN NN NN  II     TT                C
C            GG   GGG FF        II  NN  NNNN  II     TT                C
C            GG    GG FF        II  NN   NNN  II     TT                C
C             GGGGGG  FF       IIII NN    NN IIII    TT                C
C                                                                      C
C -------------------------------------------------------------------- C
C  GFINIT EVALUATES THE BOYS INCOMPLETE GAMMA FUNCTION INTEGRAL        C
C          F_M(T) = INT_{0}^{1} U^{2M} EXP(-TU^{2}) dU.                C
C                                                                      C
C  FOR 0 < M < MMAX AND 0 < T < 18.0 in STEPS OF 0.01                  C
C  USING A POWER SERIES REPRESENTATION. VALUES ARE STORED FOR USE      C
C  IN SUBSEQUENT EVALUATIONS USING THE TAYLOR SERIES                   C
C  F_N(X+H) = F_N(X) - F_{N+1}(X)H + (1/2!)F_{N+2}(X)H^2 - ...         C
C -------------------------------------------------------------------- C
C  THE VALUES X ARE CHOSEN FOR A STEP 0 < H < 0.01.                    C
C**********************************************************************C
      PARAMETER(NINT=1801,ISTEP=100,MLM=30,MSER=60)
C
      DIMENSION TMMAX(NINT),TT2(NINT),TEXP(NINT)
C
      COMMON/FNFM/FM(NINT,MLM),T(NINT)
C
C     SPECIAL CASE FOR T = 0.0D0
      T(1) = 0.0D0
      DO K=1,MMAX+1
        MVAL    = K-1
        R2MV    = DFLOAT(2*MVAL+1)
        VALUE   = 1.0D0/R2MV
        FM(1,K) = VALUE
      ENDDO
C
C     POWER SERIES EVALUATION: INITIALIZE THE POWER SERIES FOR M = MMAX
      RISTEP = DFLOAT(ISTEP)
      DO M=2,NINT
        T(M)         = DFLOAT(M-1)/RISTEP
        TEXP(M)      = DEXP(-T(M))
        TT2(M)       = 2.0D0*T(M)
        TMMAX(M)     = 1.0D0
        FM(M,MMAX+1) = 1.0D0
      ENDDO
C
C     LOOP OVER TERMS IN THE POWER SERIES.
C     CONVERGENCE ACHIEVED WHEN EXP(-T)*TERM(K)< 1.0D-14 WHERE TERM(K) 
C     IS THE KTH TERM IN THE POWER SERIES. NOTE THAT THE TERMS ARE 
C     POSITIVE SO THAT THERE IS NO NEED TO TEST FOR ABSOLUTE VALUE.
      DO K=1,MSER
        DMMAX = DFLOAT(2*(MMAX+K)+1)
        DO M=2,NINT
          TMMAX(M)     = TMMAX(M)*(TT2(M)/DMMAX)
          FM(M,MMAX+1) = FM(M,MMAX+1)+TMMAX(M)
        ENDDO
      ENDDO
C
C     RESCALE BY THE PREFACTOR
      DEN = DFLOAT((2*MMAX)+1)
      DO M=2,NINT
        FM(M,MMAX+1) = FM(M,MMAX+1)*TEXP(M)/DEN
      ENDDO
C
C     NOW COMPLETE TABLE BY BACKWARDS RECURRENCE
      DO I=1,MMAX
        MIND  = MMAX-I+1
        MVAL  = MIND-1
        COEFF = DFLOAT(MVAL+MVAL+1)
        DO M=2,NINT
          FM(M,MIND) = (TT2(M)*FM(M,MIND+1)+TEXP(M))/COEFF
        ENDDO
      ENDDO
C
      RETURN
      END
      
      
      SUBROUTINE ESETLL(LAMLL,NLL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        EEEEEEEE  SSSSSS  EEEEEEEE TTTTTTTT LL       LL               C
C        EE       SS    SS EE          TT    LL       LL               C
C        EE       SS       EE          TT    LL       LL               C
C        EEEEEE    SSSSSS  EEEEEE      TT    LL       LL               C
C        EE             SS EE          TT    LL       LL               C
C        EE       SS    SS EE          TT    LL       LL               C
C        EEEEEEEE  SSSSSS  EEEEEEEE    TT    LLLLLLLL LLLLLLLL         C
C                                                                      C
C -------------------------------------------------------------------- C
C  ESETLL CONSTRUCTS ALL ELL0 COEFFICIENTS FOR A SYSTEM WITH CALLS TO  C
C  EMAKELL AND SAVES THEM TO AN EXTERNAL DATA FILE.                    C
C**********************************************************************C
      PARAMETER(MDM=2500,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26,MB2=MBS*MBS,
     &                                        MLL=MKP*(MKP+1)*(MKP+2)/6)
C
      DIMENSION NLL(0:MKP),NSS(0:MKP),NLS(0:MKP)
C
      COMMON/NQTT/NADE0LL,NADE0SS,NADEILS
C     COMMON/E0LL/E0LLFILE(N0LLSAV,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
C     COMMON/E0SS/E0SSFILE(N0SSSAV,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
C     COMMON/EILS/EILSFILE(NILSSAV,8),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)

      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C**********************************************************************C
C     LOOP OVER FOCK BLOCK AND COUNT ALL REQUIRED EQ-WORDS.            C
C**********************************************************************C
C
C     INITIALISE MAXIMUM LAMBDA
      LAMLL = 0
C
C     INITIALISE TOTAL COEFFICIENT COUNTERS
      NADE0LL =  NADE0LL + NTUVLL*MAXAB
      NADE0SS =  NADE0SS + NTUVSS*MAXAB
      NADEILS =  NADEILS + NTUVLS*MAXAB*3
C
C     INITIALISE LAMBDA COEFFICIENT COUNTERS
      DO ILAM=0,MKP
        NLL(ILAM) = 0
        NSS(ILAM) = 0
        NLS(ILAM) = 0
      ENDDO
C
C     LOOP OVER CENTERS A AND B
      DO ICNTA=1,NCNT
        DO ICNTB=1,NCNT
C
C         LOOP OVER KQNA VALUES
          DO KA=1,NKAP(ICNTA)
C
C           QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
            IF(KVALS(KA,ICNTA).GT.0) THEN
             LQNA = KVALS(KA,ICNTA)
            ELSE
             LQNA =-KVALS(KA,ICNTA)-1
            ENDIF
            NFUNA = NFUNCT(LQNA+1,ICNTA)
C
C           LOOP OVER KQNB VALUES
            DO KB=1,NKAP(ICNTB)
C
C             QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
              IF(KVALS(KB,ICNTB).GT.0) THEN
                LQNB = KVALS(KB,ICNTB)
              ELSE
                LQNB =-KVALS(KB,ICNTB)-1
              ENDIF
              NFUNB = NFUNCT(LQNB+1,ICNTB)
C
C             LOOP OVER |MQNA| VALUES
              DO MA=1,IABS(KVALS(KA,ICNTA))
                MJA = 2*MA-1
C
C               LOOP OVER |MQNB| VALUES
                DO MB=1,IABS(KVALS(KB,ICNTB))
                  MJB = 2*MB-1
C
C                 NUMBER OF BASIS FUNCTION OVERLAPS
                  MAXAB = NFUNA*NFUNB
C
C                 CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
                  LAMAB  = LQNA+LQNB
                  NTUVLL = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
                  NTUVSS = (LAMAB+3)*(LAMAB+4)*(LAMAB+5)/6
                  NTUVLS = (LAMAB+2)*(LAMAB+3)*(LAMAB+4)/6
C
C                 UPDATE LARGEST LAMBDA VALUE
                  IF(LAMAB.GT.LAMLL) THEN
                    LAMLL = LAMAB
                  ENDIF
C
C                 INCREASE NUMBER OF WORDS FOR THIS LAMBDA VALUE
                  NLL(LAMAB  ) = NLL(LAMAB  ) + NTUVLL*MAXAB
                  NSS(LAMAB+2) = NSS(LAMAB+2) + NTUVSS*MAXAB
                  NLS(LAMAB+1) = NLS(LAMAB+1) + NTUVLS*MAXAB*3
C
C                 INCREASE TOTAL NUMBER OF WORDS
                  NADE0LL =  NADE0LL + NTUVLL*MAXAB
                  NADE0SS =  NADE0SS + NTUVSS*MAXAB
                  NADEILS =  NADEILS + NTUVLS*MAXAB*3
C
C               END LOOPS OVER |MQNA| AND |MQNB|
                ENDDO
              ENDDO
C
C           END LOOPS OVER KQNA AND KQNB
            ENDDO
          ENDDO
C
C       END LOOPS OVER ICNTA AND ICNTB
        ENDDO
      ENDDO
C
C     DOUBLE LOOP OVER FOCK BLOCK COMPLETE
C
C**********************************************************************C
C     SUMMARY OF WORD ANALYSIS                                         C
C**********************************************************************C
C
C     SECTION TITLE
20    FORMAT(1X,A,6X,A,7X,A,14X,A,16X,A)
21    FORMAT(1X,A,10X,I2,7X,I5,9X,I10,15X,F10.3)
22    FORMAT(1X,A,31X,I10,15X,F10.3)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',22),'E-coefficient word analysis'
      WRITE(7, *) REPEAT(' ',22),'E-coefficient word analysis'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,20) 'Type','Lambda','Depth','Words','Size (MB)'
      WRITE(7,20) 'Type','Lambda','Depth','Words','Size (MB)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     INITIALISE OVERALL WORD COUNTER AND SIZE COUNTER
      NADETOT = 0
      SPCETOT = 0.0D0
C
C     E0LL ANALYSIS
      DO ILAM=0,LAMLL
        IF(NLL(ILAM).EQ.0) GOTO 200
        NTUVLL = (ILAM+1)*(ILAM+2)*(ILAM+3)/6
        SPCELL = NLL(ILAM)*1.6D-5
        WRITE(6,21) 'E0LL',ILAM,NTUVLL,NLL(ILAM),SPCELL
        WRITE(7,21) 'E0LL',ILAM,NTUVLL,NLL(ILAM),SPCELL
        NADETOT = NADETOT + NLL(ILAM)
        SPCETOT = SPCETOT + SPCELL
200     CONTINUE
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
      IF(HMLTN.EQ.'NORL') GOTO 100
C
C     E0SS ANALYSIS
      DO ILAM=2,LAMLL+2
        IF(NSS(ILAM).EQ.0) GOTO 210
        NTUVSS = (ILAM+1)*(ILAM+2)*(ILAM+3)/6
        SPCESS = NSS(ILAM)*1.6D-5
        WRITE(6,21) 'E0SS',ILAM,NTUVSS,NSS(ILAM),SPCESS
        WRITE(7,21) 'E0SS',ILAM,NTUVSS,NSS(ILAM),SPCESS
        NADETOT = NADETOT + NSS(ILAM)
        SPCETOT = SPCETOT + SPCESS
210     CONTINUE
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
      IF(HMLTN.EQ.'BARE'.OR.HMLTN.EQ.'DHFR') GOTO 100
C
C     EILS ANALYSIS
      DO ILAM=1,LAMLL+1
        IF(NLS(ILAM).EQ.0) GOTO 220
        NTUVLS = (ILAM+1)*(ILAM+2)*(ILAM+3)/6
        SPCELS = NSS(ILAM)*1.6D-5
        WRITE(6,21) 'EILS',ILAM,NTUVLS,NLS(ILAM),SPCELS
        WRITE(7,21) 'EILS',ILAM,NTUVLS,NLS(ILAM),SPCELS
        NADETOT = NADETOT + NLS(ILAM)
        SPCETOT = SPCETOT + SPCELS
220     CONTINUE
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
100   CONTINUE
C
C     SUMMARY OF TOTALS
      WRITE(6,22) 'Total:',NADETOT,SPCETOT
      WRITE(7,22) 'Total:',NADETOT,SPCETOT
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C**********************************************************************C
C     GENERATE COMPLETE BATCH OF EQ-COEFFS
C**********************************************************************C
C
C     SECTION TITLE
      WRITE(6, *) REPEAT(' ',18),'Generating E-coefficient data files'
      WRITE(7, *) REPEAT(' ',18),'Generating E-coefficient data files'
C
C     E0LL COEFFICIENTS
      CALL CPU_TIME(TDM1)
      CALL ESETLL
      CALL CPU_TIME(TDM2)
      TELL = TELL+TDM2-TDM1
C
      IF(HMLTN.EQ.'NORL') GOTO 300
C
C     E0SS COEFFICIENTS
      CALL ESETSS
      CALL CPU_TIME(TDM3)
      TESS = TELL+TDM3-TDM2
C
      IF(HMLTN.EQ.'BARE'.OR.HMLTN.EQ.'DHFR') GOTO 300
C
C     EILS COEFFICIENTS
      CALL ESETLS
      CALL CPU_TIME(TDM4)
      TELS = TELS+TDM4-TDM3
C
300   CONTINUE
C
C     END OF SECTION
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
      RETURN
      END

