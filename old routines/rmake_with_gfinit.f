c
c     the GFINIT routine should be called some time during the CARDIN
c     routine, perhaps after biggest LQN has been determined:
C     GENERATE MINIMAL LIST OF BOYS INTEGRALS FOR USE IN FUNFM/RMAKE
      CALL GFINIT(4*LBIG+5)
c
C   [E] GFINIT: GAMMA FUNCTION INTEGRALS FOR USE IN FUNFM.             C
c
c
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
      PARAMETER(MLM=30,MSER=60,NINT=1801,ISTEP=100)
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
      DEN = DFLOAT(2*MMAX+1)
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

c
c
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
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MLM=30,
     &                      ML4=2*(MKP+1),MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C     
      DIMENSION FS(MB2,MLM),APH(MB2),QP(MB2,3),RC(MB2,MRC),RC2(MB2,MRC)
      DIMENSION F1(MB2,MLM),F2(MB2,MLM),F3(MB2,MLM),F4(MB2,MLM)     
      DIMENSION X1(MB2),X2(MB2),X3(MB2),X4(MB2),
     &          I1(MB2),I2(MB2),I3(MB2),I4(MB2)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
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
c
c
      SUBROUTINE FUNFM(FX,X,N,LAM,ITYPE)
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
C    ITYPE = 1 - SPECIAL CASE X = 0.0D0.                               C
C    ITYPE = 2 - POWER SERIES AND REVERSE RECURRENCE.                  C
C                (ONLY MSER TERMS WILL BE USED, SO USE MUST SUPPLY A   C
C                 VALUE APPROPRIATE TO THE MAX VALUE OF X IN BATCH).   C
C    ITYPE = 3 - ASYMPTOTIC EXPANSION AND FORWARD RECURRENCE.          C
C    ITYPE = 4 - ASYMPTOTIC EXPANSION AND FORWARD RECURRENCE.          C
C                ALL TERMS DEPENDING ON EXP(-X) ARE OMITTED TO AVOID   C
C                NUMERICAL UNDERFLOW PROBLEMS. MSER NOT REQUIRED.      C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MLM=30,MSER=60,NINT=1801,ISTEP=100)
C
      DIMENSION FX(MB2,MLM),X(MB2),XLAM(MB2),XX2(MB2),
     &          XEXP(MB2),XROOT(MB2),HSTEP(MB2),XSTEP(MB2),JM(MB2)
C
      COMMON/FNFM/FM(NINT,MLM),T(NINT)
C     
      DATA PIROOT,A0,B0/8.862269254527580D-1,4.994501191201870D-1,
     &                                       4.551838436668326D-1/
C
C**********************************************************************C
C     ITYPE = 1: SPECIAL CASE FOR T = 0.0D0                            C
C**********************************************************************C
C
      IF(ITYPE.EQ.1) THEN
        DO K=1,LAM+1
          MVAL  = K-1
          VALUE = 1.0D0/DFLOAT(MVAL+MVAL+1)
          DO M=1,N
            FX(M,K) = VALUE
          ENDDO
        ENDDO
        RETURN
C
C**********************************************************************C
C     ITYPE = 2: POWER SERIES EVALUATION                               C
C**********************************************************************C
C
      ELSEIF(ITYPE.EQ.2) THEN
C
C       INITIALIZE THE POWER SERIES FOR M = LAM
        DO M=1,N
          IM          = (DFLOAT(ISTEP)*X(M)) + 0.5D0
          IM          = IM+1
          JM(M)       = IM
          FX(M,LAM+1) = FM(IM,LAM+1)
          HSTEP(M)    = X(M) - T(IM)
          XSTEP(M)    =-HSTEP(M)
          XEXP(M)     = DEXP(-X(M))
          XX2(M)      = 2.0D0*X(M)
        ENDDO

        DO ITERM=1,3
          DO M=1,N
            FX(M,LAM+1) = FX(M,LAM+1) + XSTEP(M)*FM(JM(M),LAM+ITERM+1)
            XSTEP(M)    =-XSTEP(M)*HSTEP(M)/DFLOAT(ITERM+1)
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
C**********************************************************************C
C     ITYPE = 3: ASYMPTOTIC EXPANSION WITH VERY LARGE ARGUMENT.        C
C**********************************************************************C
C
      ELSEIF(ITYPE.EQ.3) THEN
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
        DO MIND=1,LAM
          COEFF = DFLOAT(2*MIND-1)
          DO M=1,N
            FX(M,MIND+1) = (COEFF*FX(M,MIND)-XEXP(M))/XX2(M)
          ENDDO
        ENDDO
        RETURN
C
C**********************************************************************C
C     ITYPE = 4: ASYMPTOTIC EXPANSION WITH VERY LARGE ARGUMENT         C
C**********************************************************************C
C
      ELSEIF(ITYPE.EQ.4) THEN
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
C**********************************************************************C
C     ITYPE OUT OF RANGE: INVALID INPUT TO FUNFM                       C
C**********************************************************************C
C
      ELSE
91      FORMAT(2X,'In FUNFM: invalid type (must be 1-4)',I4)
        WRITE(6,91) ITYPE
        WRITE(7,91) ITYPE
        STOP
      ENDIF
C
      RETURN
      END

