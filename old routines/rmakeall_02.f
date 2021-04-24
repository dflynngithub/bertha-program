      SUBROUTINE FUNFM(FM,T,N,LAM,ITYPE)
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
C  FOR VARIABLE T > 0 FOR ALL ORDERS 0 < M < LAM.                      C
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
      PARAMETER(MBS=26,MB2=MBS*MBS,MLM=30,MSER=60)
C
      DIMENSION FM(MB2,MLM),T(MB2),TLAM(MB2),TT2(MB2),
     &          TXP(MB2),TRT(MB2)
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
            FM(M,K) = VALUE
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
C       INITIALISE THE POWER SERIES FOR M = LAM
        DO M=1,N
          TXP(M)      = DEXP(-T(M))
          TT2(M)      = 2.0D0*T(M)
          TLAM(M)     = 1.0D0
          FM(M,LAM+1) = 1.0D0
        ENDDO
C
C       LOOP OVER TERMS IN THE POWER SERIES
        DO K=1,MSER
          DLAM = DFLOAT(2*(LAM+K)+1)
          DO M=1,N
            TLAM(M)     = TLAM(M)*(TT2(M)/DLAM)
            FM(M,LAM+1) = FM(M,LAM+1) + TLAM(M)
          ENDDO
        ENDDO
C
C       RESCALE BY THE PREFACTOR
        DEN = DFLOAT((2*LAM)+1)
        DO M=1,N
          FM(M,LAM+1) = FM(M,LAM+1)*TXP(M)/DEN
        ENDDO
C
C       NOW COMPLETE TABLE BY BACKWARDS RECURRENCE
        DO I=1,LAM
          MIND  = LAM-I+1
          MVAL  = MIND-1
          COEFF = DFLOAT(MVAL+MVAL+1)
          DO M=1,N
            FM(M,MIND) = (TT2(M)*FM(M,MIND+1) + TXP(M))/COEFF
          ENDDO
        ENDDO
C
C**********************************************************************C
C     ITYPE = 3: ASYMPTOTIC EXPANSION WITH VERY LARGE ARGUMENT.        C
C**********************************************************************C
C
      ELSEIF(ITYPE.EQ.3) THEN
C
C       INITIALISE THE ASYMPTOTIC EXPANSION
        DO M=1,N
          TXP(M) = DEXP(-T(M))
          TT2(M) = 2.0D0*T(M)
          TRT(M) = DSQRT(T(M))
        ENDDO
C
C       SEED VALUES
        DO M=1,N
          FM(M,1) = A0/(B0+T(M))
        ENDDO
C
C       RESCALE BY THE PREFACTOR
        DO M=1,N
          FM(M,1) = (PIROOT/TRT(M)) - (TXP(M)*FM(M,1))
        ENDDO
C
C       NOW COMPLETE TABLE BY FORWARD RECURRENCE
        DO MIND=1,LAM
          MVAL  = MIND-1
          COEFF = DFLOAT(MVAL+MVAL+1)
          DO M=1,N
            FM(M,MIND+1) = (COEFF*FM(M,MIND) - TXP(M))/TT2(M)
          ENDDO
        ENDDO
C
C**********************************************************************C
C     ITYPE = 4: ASYMPTOTIC EXPANSION WITH VERY LARGE ARGUMENT         C
C**********************************************************************C
C
      ELSEIF(ITYPE.EQ.4) THEN
C
C     INITIALISE THE ASYMPTOTIC EXPANSION
      DO M=1,N
        TT2(M)  = 2.0D0*T(M)
        FM(M,1) = PIROOT/DSQRT(T(M))
      ENDDO
C
C     NOW COMPLETE TABLE BY FORWARD RECURRENCE
      DO MIND=1,LAM
        MVAL  = MIND-1
        COEFF = DFLOAT(MVAL+MVAL+1)
        DO M=1,N
          FM(M,MIND+1) = (COEFF*FM(M,MIND))/TT2(M)
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
C
C
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
C  FINITE SUM REPRESENTATION OF A MULTI-CENTRE GAUSSIAN OVERLAP.       C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MLM=30,MKP=9,
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
C
        X = APH(M)*(QP(M,1)*QP(M,1)+QP(M,2)*QP(M,2)+QP(M,3)*QP(M,3))
C
C       CASE 1: IF X IS ALMOST ZERO (SO WHEN Q=P OR EIJ<<1)
        IF(X.LE.1.0D-11) THEN
          N1     = N1+1
          X1(N1) = X
          I1(N1) = M
C
C       CASE 2: IF X IS SMALLER THAN 17.0D0
        ELSEIF(X.GT.1.0D-11.AND.X.LE.17.0D0) THEN
          N2     = N2+1
          X2(N2) = X
          I2(N2) = M
C
C       CASE 3: IF X IS SMALLER THAN 30.0D0
        ELSEIF(X.GT.17.0D0.AND.X.LE.30.0D0) THEN
          N3     = N3+1
          X3(N3) = X
          I3(N3) = M
C
C       CASE 4: IF X IS LARGER THAN 30.0D0
        ELSE
          N4     = N4+1
          X4(N4) = X
          I4(N4) = M
        ENDIF
C
      ENDDO
C
C     EVALUATE THE BOYS INTEGRALS -- A BATCH FOR EACH ITYPE
C
C     CASE 1: ARGUMENT OF THE BOYS FUNCTION IS X=0.0D0.
C     THE VALUE OF THIS FUNCTION IS 2N+1 (DONE IN FUNFM).
      IF(N1.GT.0) THEN
        CALL FUNFM(F1,X1,N1,LAM,1)
        DO K=1,LAM+1
          DO M=1,N1
            FS(I1(M),K) = F1(M,K)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 2: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=17.0D0.
C     EVALUATE WITH LOCAL POLYNOMIAL EXPANSION OF ORDER 5,
C     AND RECURRENCE IN DIRECTION OF DECREASING M.
      IF(N2.GT.0) THEN
        CALL FUNFM(F2,X2,N2,LAM,2)
        DO K=1,LAM+1
          DO M=1,N2
            FS(I2(M),K) = F2(M,K)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 3: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=30.0D0.
C     EVALUATE USING ASYMPTOTIC FORMULA WITH EXPONENTIAL.
      IF(N3.GT.0) THEN
        CALL FUNFM(F3,X3,N3,LAM,3)
        DO K=1,LAM+1
          DO M=1,N3
            FS(I3(M),K) = F3(M,K)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 4: ARGUMENT OF THE BOYS FUNCTION IS LARGER THAN X=30.0D0.
C     EVALUATE USING ASYMPTOTIC FORMULA WITHOUT EXPONENTIAL.
      IF(N4.GT.0) THEN
        CALL FUNFM(F4,X4,N4,LAM,4)
        DO K=1,LAM+1
          DO M=1,N4
            FS(I4(M),K) = F4(M,K)
          ENDDO
        ENDDO
      ENDIF
C
C**********************************************************************C
C     WITH THE FULL SET OF BOYS' INTEGRALS WE NOW APPLY RECURRENCE     C
C     RELATIONS TO F_N (X) AND EVALUATE THE R-INTEGRALS.               C
C**********************************************************************C
C
C     CONSTRUCT TOP LEVEL (FOR MAXIMUM LAM VALUE)
      DO M=1,MAXM
        RC(M,1)=((-2.0D0*APH(M))**(LAM))*FS(M,LAM+1)
      ENDDO
C
C     MINIMUM LEVEL ILEV BASED ON LAM VALUE
      IF(MOD(LAM,2).EQ.0) THEN
       ITUVMIN = 1
      ELSE
       ITUVMIN = 2
      ENDIF
C
C     INITIALISE ITUV COUNTER (RELATES TO # CARTESIAN INDICES FOR lam)
      ITUV=-1
C
C     MAIN LOOP: LEVEL 'ILEV' STARTING AT LAM-1 AND WORKING BACKWARDS
      DO ILEV=LAM-1,ITUVMIN,-2
C
C       UPDATE ITUV COUNTER
        ITUV = ITUV+1
C
C       LOOP OVER ALL UNIQUE (IT,IU,IV)
        DO IT=0,ITUV
          RIT = DFLOAT(IT)
C
          DO IU=0,ITUV-IT
            RIU=DFLOAT(IU)
C
            DO IV=0,ITUV-IT-IU
              RIV=DFLOAT(IV)
C
C             R-INTEGRAL CARTESIAN DESTINATION ADDRESSES
              N1 = INABCD(IT+1,IU  ,IV  )
              N2 = INABCD(IT  ,IU+1,IV  )
              N3 = INABCD(IT  ,IU  ,IV+1)
C
C             RECURRENCE RELATIONS DIFFERENT IF ANY (IT,IU,IV) ARE ZERO

C             CASE (IT,IU,IV)
              IF(IT.NE.0.AND.IU.NE.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE (IT,IU, 0)
              ELSEIF(IT.NE.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
               ENDDO
C
C             CASE (IT, 0,IV)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE (IT, 0, 0)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE ( 0,IU,IV)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE ( 0,IU, 0)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE ( 0, 0,IV)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE ( 0, 0, 0)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
                ENDDO
C
C             ALL POSSIBILITIES ACCOUNTED FOR -- END LOOP
              ENDIF
C
C           END LOOPS OVER (IT,IU,IV) ADDRESSES FOR THIS ILEV
            ENDDO
          ENDDO
        ENDDO
C
C       ADD IN (IT=0,IU=0,IV=0) CASE
        DO M=1,MAXM
          RC2(M,1) = ((-2.0D0*APH(M))**(ILEV))*FS(M,ILEV+1)
        ENDDO
C
C
C       UPDATE ITUV COUNTER
        ITUV = ITUV+1
C
C       LOOP OVER ALL UNIQUE (IT,IU,IV)
        DO IT=0,ITUV
          RIT=DFLOAT(IT)
          DO IU=0,ITUV-IT
            RIU=DFLOAT(IU)
            DO IV=0,ITUV-IT-IU
              RIV=DFLOAT(IV)
C
C             R-INTEGRAL CARTESIAN DESTINATION ADDRESSES
              N1 = INABCD(IT+1,IU  ,IV  )
              N2 = INABCD(IT  ,IU+1,IV  )
              N3 = INABCD(IT  ,IU  ,IV+1)
C
C             RECURRENCE RELATIONS DIFFERENT IF ANY (IT,IU,IV) ARE ZERO

C             CASE (IT,IU,IV)
              IF(IT.NE.0.AND.IU.NE.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC(M,N1) =-QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                  RC(M,N2) =-QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) =-QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
C
C             CASE (IT,IU, 0)
              ELSEIF(IT.NE.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC(M,N1) =-QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                  RC(M,N2) =-QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) =-QP(M,3)*RC2(M,K1)
                ENDDO
C
C             CASE (IT, 0,IV)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC(M,N1) =-QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                  RC(M,N2) =-QP(M,2)*RC2(M,K1)
                  RC(M,N3) =-QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
C
C             CASE (IT, 0, 0)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                DO M=1,MAXM
                 RC(M,N1) =-QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                 RC(M,N2) =-QP(M,2)*RC2(M,K1)
                 RC(M,N3) =-QP(M,3)*RC2(M,K1)
               ENDDO
C
C             CASE ( 0,IU,IV)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC(M,N1) =-QP(M,1)*RC2(M,K1)
                  RC(M,N2) =-QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) =-QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
C
C             CASE ( 0,IU, 0)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC(M,N1) =-QP(M,1)*RC2(M,K1)
                  RC(M,N2) =-QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) =-QP(M,3)*RC2(M,K1)
                ENDDO
C
C             CASE ( 0, 0,IV)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC(M,N1) =-QP(M,1)*RC2(M,K1)
                  RC(M,N2) =-QP(M,2)*RC2(M,K1)
                  RC(M,N3) =-QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
C
C             CASE ( 0, 0, 0)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                DO M=1,MAXM
                  RC(M,N1) =-QP(M,1)*RC2(M,K1)
                  RC(M,N2) =-QP(M,2)*RC2(M,K1)
                  RC(M,N3) =-QP(M,3)*RC2(M,K1)
                ENDDO
C
C             ALL POSSIBILITIES ACCOUNTED FOR -- END LOOP
              ENDIF
C
C           END LOOPS OVER (IT,IU,IV) ADDRESSES FOR THIS ILEV
            ENDDO
          ENDDO
        ENDDO
C
C       ADD IN (IT=0,IU=0,IV=0) CASE
C
        DO M=1,MAXM
          RC(M,1) = ((-2.0D0*APH(M))**(ILEV-1))*FS(M,ILEV)
        ENDDO
C
      ENDDO
C
C
C     AN ADDITIONAL LOOP OVER ADDRESSES (WHEN LAM IS ODD)
      IF(MOD(LAM,2).EQ.1) THEN
C
C       UPDATE ITUV COUNTER
        ITUV = ITUV+1
C
C       LOOP OVER ALL UNIQUE (IT,IU,IV)
        DO IT=0,ITUV
          RIT=DFLOAT(IT)
          DO IU=0,ITUV-IT
            RIU=DFLOAT(IU)
            DO IV=0,ITUV-IT-IU
              RIV=DFLOAT(IV)
C
C             R-INTEGRAL CARTESIAN DESTINATION ADDRESSES
              N1 = INABCD(IT+1,IU  ,IV  )
              N2 = INABCD(IT  ,IU+1,IV  )
              N3 = INABCD(IT  ,IU  ,IV+1)
C
C             RECURRENCE RELATIONS DIFFERENT IF ANY (IT,IU,IV) ARE ZERO

C             CASE (IT,IU,IV)
              IF(IT.NE.0.AND.IU.NE.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE (IT,IU, 0)
              ELSEIF(IT.NE.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE (IT, 0,IV)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE (IT, 0, 0)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE ( 0,IU,IV)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE ( 0,IU, 0)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE ( 0, 0,IV)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE ( 0, 0, 0)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
               ENDDO
C
C             ALL POSSIBILITIES ACCOUNTED FOR -- END LOOP
              ENDIF
C
C           END LOOPS OVER (IT,IU,IV) ADDRESSES FOR THIS ILEV
            ENDDO
          ENDDO
        ENDDO
C
C       ADD IN (IT=0,IU=0,IV=0) CASE
C
        DO M=1,MAXM
          RC2(M,1) = FS(M,1)
        ENDDO
C
C       MOVE THE RC2 ARRAY INTO RC
C
        ITMAX = (LAM+1)*(LAM+2)*(LAM+3)/6
        DO IT=1,ITMAX
          DO M=1,MAXM
            RC(M,IT) = RC2(M,IT)
          ENDDO
        ENDDO
C
C     END IF STATEMENT FOR THE ODD LAM CASE
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
C  WMAKE CALCULATES A BATCH OF HERMITE COULOMB INTEGRALS (THE INTEGRAL C
C  OVER A PRODUCT OF A POINT-NUCLEUS POTENTIAL AND ALL HGTFS).         C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MLM=30,
     &                      ML4=2*(MKP+1),MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      DIMENSION FS(MB2,MLM),APH(MB2),RP(MB2,3),WC(MB2,MRC)
      DIMENSION X1(MB2),X2(MB2),X3(MB2),X4(MB2),
     &          I1(MB2),I2(MB2),I3(MB2),I4(MB2)
      DIMENSION F1(MB2,MLM),F2(MB2,MLM),F3(MB2,MLM),F4(MB2,MLM)
      DIMENSION H(0:MKP-1,0:(MKP-1)/2,3)
C
      COMMON/ACCSS/INABCD(0:ML4,0:ML4,0:ML4),
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
      CALL GAMGEN
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
          I1(N1) = M
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

