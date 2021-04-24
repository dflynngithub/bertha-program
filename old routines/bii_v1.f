      SUBROUTINE BII(RR,XYZ,KQN,MQN,EXPT,NBAS,IBAS,JBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                          BBBBBBB  IIII IIII                          C
C                          BB    BB  II   II                           C
C                          BB    BB  II   II                           C
C                          BBBBBBB   II   II                           C
C                          BB    BB  II   II                           C
C                          BB    BB  II   II                           C
C                          BBBBBBB  IIII IIII                          C
C                                                                      C
C -------------------------------------------------------------------- C
C  ERI GENERATES A BATCH OF MOLECULAR ELECTRON REPULSION INTEGRALS BY  C
C  MEANS OF THE MCMURCHIE-DAVIDSION ALGORITHM (DOUBLE FINITE SUM OVER  C
C  EQ-COEFFICIENTS AND INTEGRALS OVER A PAIR OF HGTFS.)                C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    XYZ(3,4)  - COORDINATES OF THE 4 NUCLEAR CENTERS.                 C
C    KQN(4)    - KQN QUANTUM NUMBERS OF THE CENTERS.                   C
C    MQN(4)    - |MQN| QUANTUM NUMBERS OF THE CENTERS.                 C
C    EXPT(I,4) - EXPONENTS IN BASIS FUNCTION BLOCK FOR A,B,C,D         C
C    NBAS(J)   - NUMBER OF FUNCTIONS ON CENTER J.                      C
C    IBAS,JBAS - COMPONENT LABEL INDEX FOR BASIS FUNCTION PAIR ON AB . C
C    IEAB,IECD - 0 DON'T RECALCULATE E(AB/CD)-COEFFICIENTS             C
C                1 DO    RECALCULATE E(AB/CD)-COEFFICIENTS             C
C  OUTPUT:                                                             C
C    RR(MB2,16) - BII'S FOR BLOCK AB, ALL 16 MQN PROJECTION COMBOS.    C
C -------------------------------------------------------------------- C
C  DFNOTE: THERE ARE NO SHORTCUTS FOR INUCAB OR INUCCD CONDITIONS.     C
C          IBSCR DATA HASN'T BEEN UTILISED YET.                        C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MCT=15,MKP=9,MFL=7000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      CHARACTER*4 HMLTN
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
      DIMENSION PQ(MB2,3),APH(MB2),PRE(MB2)
      DIMENSION IAB11(MEQ,3),IAB21(MEQ,3),ICD11(MEQ,3),ICD21(MEQ,3)
      DIMENSION RC(MB2,MRC),IRC(MRC)
C
      COMPLEX*16 CONE
      COMPLEX*16 RR(MB2,16),Q1(MB2),Q2(MB2)
      COMPLEX*16 EAB11(MB2,MEQ,3),EAB21(MB2,MEQ,3),BAB11(MB2,MEQ),
     &           ECD11(MB2,MEQ,3),ECD21(MB2,MEQ,3),BAB21(MB2,MEQ)
C
      SAVE EAB11,EAB21,IAB11,IAB21
      SAVE ECD11,ECD21,ICD11,ICD21
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/EILS/EILSFL(MFL,24),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/IBSC/IBSCR(MB2),IBMAP(MB2)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS
      COMMON/MAKE/IEAB,IECD,IRIJ(MBS,MBS)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,IEQS,IERC,IPAR,ICOR,ILEV
      COMMON/TMMD/TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR
      COMMON/TSCF/TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRR,TCC1,TCC2,TCMC,
     &            TB1B,TB1R,TB1F,TB1M,TBEC,TBRM,TBRR,TBC1,TBC2,TBMC,
     &            THMX,TC1T,TC2T,TB1T,TB2T,TEIG,TSCR,TTOT,
     &            TC1S,TC2S,TB1S,TB2S
C
      DATA ROOTPI5,SENS/17.4934183276248628D0,1.0D-10/
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     EVALUATE LQNS FOR BASIS FUNCTIONS A, B, C, D
      DO N=1,4
        IF(KQN(N).LT.0) THEN
         LQN(N) =-KQN(N)-1
        ELSE
         LQN(N) = KQN(N)
        ENDIF
      ENDDO
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NBAS(1)*NBAS(2)
      MAXCD = NBAS(3)*NBAS(4)
C
C     PHASE FACTOR FOR AB AND CD PAIR OVERLAPS
      IPHSAB = 1
      IPHSCD =-1
C
C     VRS MAXIMUM LAMBDA LEVEL FOR EQ-COEFFICIENT ADDRESSES
      LAMAB = LQN(1)+LQN(2)+1
      LAMCD = LQN(3)+LQN(4)+1
C
C     VRS TOTAL LENGTH OF EQ-COEFFICIENT LISTS
      NTUVAB = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
      NTUVCD = (LAMCD+1)*(LAMCD+2)*(LAMCD+3)/6
C
C     VRS MAXIMUM LAMBDA LEVEL FOR R-INTEGRAL ADDRESSES
      LAMABCD = LAMAB+LAMCD
C
C     VRS TOTAL LENGTH OF R-INTEGRAL LISTS
      NTUVABCD = (LAMABCD+1)*(LAMABCD+2)*(LAMABCD+3)/6
C
C**********************************************************************C
C     IF ASKED TO RECALCULATE E(AB) COEFFICIENTS, DO THIS FIRST        C
C**********************************************************************C
C
      IF(IEAB.EQ.1) THEN
C
C       GENERATE ELS(AB) COEFFICIENTS
        CALL CPU_TIME(TDM1)
        IF(IEQS.EQ.0) THEN
          CALL EMAKEB3(EAB11,EAB21,EXPT,XYZ,KQN,MQN,NBAS,IPHSAB,1,2)
        ELSEIF(IEQS.EQ.1) THEN
          DO ITUV=1,NTUVAB
            IAD = IABLS + (ITUV-1)*MAXAB
            DO M=1,MAXAB
              EAB11(M,ITUV,1)=DCMPLX(EILSFL(IAD+M, 1),EILSFL(IAD+M, 2))
              EAB21(M,ITUV,1)=DCMPLX(EILSFL(IAD+M, 3),EILSFL(IAD+M, 4))
              EAB11(M,ITUV,2)=DCMPLX(EILSFL(IAD+M, 5),EILSFL(IAD+M, 6))
              EAB21(M,ITUV,2)=DCMPLX(EILSFL(IAD+M, 7),EILSFL(IAD+M, 8))
              EAB11(M,ITUV,3)=DCMPLX(EILSFL(IAD+M, 9),EILSFL(IAD+M,10))
              EAB21(M,ITUV,3)=DCMPLX(EILSFL(IAD+M,11),EILSFL(IAD+M,12))
            ENDDO
          ENDDO
        ENDIF
        CALL CPU_TIME(TDM2)
        TELS = TELS + TDM2 - TDM1
C
C       SCREENING PROCEDURE: NORM SUM OF E-COEFFICIENT LIST FOR EACH IAB
        DO IQ=1,3
          DO IAB=1,NTUVAB
C
C           11 OVERLAP (AB PAIRS)
            SUM = 0.0D0
            DO M=1,MAXAB
              SUM = SUM + ABS(EAB11(M,IAB,IQ))
            ENDDO
C
            IF(SUM.LE.SENS) THEN
              IAB11(IAB,IQ) = 0
            ELSE
              IAB11(IAB,IQ) = 1
            ENDIF
C
C           21 OVERLAP (AB PAIRS)
            SUM = 0.0D0
            DO M=1,MAXAB
              SUM = SUM + ABS(EAB21(M,IAB,IQ))
            ENDDO
C
            IF(SUM.LE.SENS) THEN
              IAB21(IAB,IQ) = 0
            ELSE
              IAB21(IAB,IQ) = 1
            ENDIF
C
          ENDDO
        ENDDO
C
C       DO NOT CALCULATE AGAIN UNTIL ASKED EXTERNALLY
        IEAB = 0
C
      ENDIF
C
C**********************************************************************C
C     IF ASKED TO RECALCULATE E(CD) COEFFICIENTS, DO THIS NEXT         C
C**********************************************************************C
C
      IF(IECD.EQ.1) THEN
C
C       GENERATE ELS(CD) COEFFICIENTS
        CALL CPU_TIME(TDM1)
        IF(IEQS.EQ.0) THEN
          CALL EMAKEB3(ECD11,ECD21,EXPT,XYZ,KQN,MQN,NBAS,IPHSCD,3,4)
        ELSEIF(IEQS.EQ.1) THEN
          DO ITUV=1,NTUVCD
            IAD = ICDLS + (ITUV-1)*MAXCD
            DO M=1,MAXCD
              ECD11(M,ITUV,1)=DCMPLX(EILSFL(IAD+M,13),EILSFL(IAD+M,14))
              ECD21(M,ITUV,1)=DCMPLX(EILSFL(IAD+M,15),EILSFL(IAD+M,16))
              ECD11(M,ITUV,2)=DCMPLX(EILSFL(IAD+M,17),EILSFL(IAD+M,18))
              ECD21(M,ITUV,2)=DCMPLX(EILSFL(IAD+M,19),EILSFL(IAD+M,20))
              ECD11(M,ITUV,3)=DCMPLX(EILSFL(IAD+M,21),EILSFL(IAD+M,22))
              ECD21(M,ITUV,3)=DCMPLX(EILSFL(IAD+M,23),EILSFL(IAD+M,24))
            ENDDO
          ENDDO
        ENDIF
        CALL CPU_TIME(TDM2)
        TELS = TELS + TDM2 - TDM1
C
C       SCREENING PROCEDURE: NORM SUM OF E-COEFFICIENT LIST FOR EACH ICD
        DO IQ=1,3
          DO ICD=1,NTUVCD
C
C           11 OVERLAP (CD PAIRS)
            SUM = 0.0D0
            DO M=1,MAXCD
              SUM = SUM + ABS(ECD11(M,ICD,IQ))
            ENDDO
C
            IF(SUM.LE.SENS) THEN
              ICD11(ICD,IQ) = 0
            ELSE
              ICD11(ICD,IQ) = 1
            ENDIF
C
C           21 OVERLAP (CD PAIRS)
            SUM = 0.0D0
            DO M=1,MAXCD
              SUM = SUM + ABS(ECD21(M,ICD,IQ))
            ENDDO
C
            IF(SUM.LE.SENS) THEN
              ICD21(ICD,IQ) = 0
            ELSE
              ICD21(ICD,IQ) = 1
            ENDIF
C
          ENDDO
        ENDDO
C
C       DO NOT CALCULATE AGAIN UNTIL ASKED EXTERNALLY
        IECD = 0
C
      ENDIF
C
C**********************************************************************C
C     R-INTEGRAL EVALUATION                                            C
C**********************************************************************C
C
C     GAUSSIAN OVERLAP VALUES
      EIJ = EXPT(IBAS,1) + EXPT(JBAS,2)
      PX  = (XYZ(1,1)*EXPT(IBAS,1) + XYZ(1,2)*EXPT(JBAS,2))/EIJ
      PY  = (XYZ(2,1)*EXPT(IBAS,1) + XYZ(2,2)*EXPT(JBAS,2))/EIJ
      PZ  = (XYZ(3,1)*EXPT(IBAS,1) + XYZ(3,2)*EXPT(JBAS,2))/EIJ
C
C     AUXILLIARY DATA FOR RMAKE ROUTINE
      M = 0
      N = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
          IF(IBSCR(M).EQ.1) THEN
            N   = N+1
            EKL = EXPT(KBAS,3) + EXPT(LBAS,4)
            QX  = (XYZ(1,3)*EXPT(KBAS,3) + XYZ(1,4)*EXPT(LBAS,4))/EKL
            QY  = (XYZ(2,3)*EXPT(KBAS,3) + XYZ(2,4)*EXPT(LBAS,4))/EKL
            QZ  = (XYZ(3,3)*EXPT(KBAS,3) + XYZ(3,4)*EXPT(LBAS,4))/EKL
            EIJKL   = EIJ+EKL
            APH(N)  = EIJ*EKL/EIJKL
            PQ(N,1) = QX-PX
            PQ(N,2) = QY-PY
            PQ(N,3) = QZ-PZ
            EMX     = DSQRT(EIJ+EKL)*EIJ*EKL
            PRE(N)  = 2.0D0*ROOTPI5/EMX
          ENDIF
        ENDDO
      ENDDO
C
C     M WAS THE LONGER COUNTER, AND N CAN BE SHORTENED BY USING IBSCR
C     MAXN COUNTS # OVERLAPS THAT SURVIVE SCREENING
      MAXM = M
      MAXN = N
C
C     GENERATE R-INTEGRALS AND RECORD TIME
      CALL CPU_TIME(TDM1)
      CALL RMAKE(RC,PQ,APH,MAXN,LAMABCD+2)
      CALL CPU_TIME(TDM2)
      TRBR = TRBR + TDM2 - TDM1
C
C     SCREENING PROCEDURE: NORM SUM OF R-INTEGRAL LIST FOR EACH IABCD
      DO IABCD=1,NTUVABCD
C
C       VECTOR MAGNITUDE OF FULL R-INTEGRAL LIST
        SUM = 0.0D0
        DO MR=1,MAXN
          SUM = SUM + DABS(RC(MR,IABCD))
        ENDDO
C
        IF(SUM.LE.SENS) THEN
          IRC(IABCD) = 0
        ELSE
          IRC(IABCD) = 1
        ENDIF
C
      ENDDO
C
C     INITIALISE RR ARRAY
      DO ITG=1,16
        DO M=1,MAXCD
          RR(M,ITG) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CONSTRUCT INTERMEDIATE MATRICES FOR MCMURCHIE-DAVIDSON           C
C**********************************************************************C
C
C     BEGIN FIRST LOOP: CARTESIAN INDEX ICMP FOR CENTER AB (USE 6000)
      DO 6000 ICMP=1,3
C
C     CARTESIAN INDICES FOR {IA,IB,IC}
      IF(ICMP.EQ.1) THEN
        IDX = 1
        IDY = 0
        IDZ = 0
      ELSEIF(ICMP.EQ.2) THEN
        IDX = 0
        IDY = 1
        IDZ = 0
      ELSEIF(ICMP.EQ.3) THEN
        IDX = 0
        IDY = 0
        IDZ = 1
      ENDIF
C
C     INITIATE LOOP OVER ALL INDICES {IA,IB,IC} FOR AB
      DO IAB=1,NTUVAB
C
        IAB1 = 0
        IAB2 = 0
C
C       SCREENING MARKERS
        IF((IAB11(IAB,1)+IAB11(IAB,2)+IAB11(IAB,3)).NE.0) THEN
          IAB1 = 1
        ENDIF
C
        IF((IAB21(IAB,1)+IAB21(IAB,2)+IAB21(IAB,3)).NE.0) THEN
          IAB2 = 1
        ENDIF
C
C       INITIALISE IF ANY E-COEFF. FOR THIS {IA,IB,IC} PASSES TEST
        IF(IAB1+IAB2.GT.0) THEN
          DO N=1,MAXN
            BAB11(N,IAB) = DCMPLX(0.0D0,0.0D0)
            BAB21(N,IAB) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDIF
C
C       INITIATE LOOP OVER ALL INDICES {IA',IB',IC'} FOR CD
        DO JCMP=1,3
C
C         CARTESIAN INDICES FOR {IA',IB',IC'}
          IF(JCMP.EQ.1) THEN
            JDX = 1
            JDY = 0
            JDZ = 0
          ELSEIF(JCMP.EQ.2) THEN
            JDX = 0
            JDY = 1
            JDZ = 0
          ELSEIF(JCMP.EQ.3) THEN
            JDX = 0
            JDY = 0
            JDZ = 1
          ENDIF
C
C ***     SPECIAL CASE: CARTESIAN INDICES ARE EQUAL {BXX, BYY, BZZ}
          IF(ICMP.EQ.JCMP) THEN
C
C           BEGIN LOOP OVER ALL (CD) ADDRESSES FOR A GIVEN (AB) ADDRESS
            DO ICD=1,NTUVCD
C
C             OVERALL ADDRESS OF THIS {IA,IB,IC},{IA',IB',IC'}
              IRABCD = INABCD(IVEC(IAB)+IVEC(ICD),JVEC(IAB)+JVEC(ICD),
     &                                            KVEC(IAB)+KVEC(ICD))
C
C             SCREEN CALCULATION IF THE BIGGEST R-INTEGRAL IS TOO SMALL
              IF(IRC(IRABCD).EQ.0) GOTO 797
C
C             IF E11(CD) IS BIG ENOUGH, CONTINUE TO ADD ALL OVERLAPS
              IF(ICD11(ICD,JCMP).EQ.1) THEN
                DO N=1,MAXN
                  BAB11(N,IAB) = BAB11(N,IAB)
     &                         - ECD11(IBMAP(N),ICD,JCMP)*RC(N,IRABCD)
                ENDDO
              ENDIF
C
C             IF E21(CD) IS BIG ENOUGH, CONTINUE TO ADD ALL OVERLAPS
              IF(ICD21(ICD,JCMP).EQ.1) THEN
                DO N=1,MAXN
                  BAB21(N,IAB) = BAB21(N,IAB)
     &                         - ECD21(IBMAP(N),ICD,JCMP)*RC(N,IRABCD)
                ENDDO
              ENDIF
C
797         CONTINUE
            ENDDO
C
C ***     END CONDITIONAL OVER LIKE INDICES
          ENDIF
C
C         LOOP OVER FINITE SUM INDICES {IA',IB',IC'} FOR CD
          DO ICD=1,NTUVCD
C
            IF(JCMP.EQ.1) THEN
              RTP = DFLOAT(IVEC(IAB) + IVEC(ICD))
            ELSEIF(JCMP.EQ.2) THEN
              RTP = DFLOAT(JVEC(IAB) + JVEC(ICD))
            ELSE
              RTP = DFLOAT(KVEC(IAB) + KVEC(ICD))
            ENDIF
C
C           STARTING ADDRESSES FOR CARTESIAN COMPONENTS
            I1 = IVEC(IAB) + IVEC(ICD) + IDX + JDX
            J1 = JVEC(IAB) + JVEC(ICD) + IDY + JDY
            K1 = KVEC(IAB) + KVEC(ICD) + IDZ + JDZ
            IADR1 = INABCD(I1,J1,K1)
C
            I2 = IVEC(IAB) + IVEC(ICD) + IDX
            J2 = JVEC(IAB) + JVEC(ICD) + IDY
            K2 = KVEC(IAB) + KVEC(ICD) + IDZ
            IADR2 = INABCD(I2,J2,K2)
C
            I3 = IVEC(IAB) + IVEC(ICD) + IDX - JDX
            J3 = JVEC(IAB) + JVEC(ICD) + IDY - JDY
            K3 = KVEC(IAB) + KVEC(ICD) + IDZ - JDZ
            IF((I3.GE.0).AND.(J3.GE.0).AND.(K3.GE.0)) THEN
              IADR3 = INABCD(I3,J3,K3)
            ELSE
              IADR3 = 1
              RTP   = 0.0D0
            ENDIF
C
C           IF E11(CD) IS BIG ENOUGH, CONTINUE TO ADD ALL OVERLAPS
            IF(ICD11(ICD,JCMP).EQ.1) THEN
              DO N=1,MAXN
                T1 = 0.5D0*RC(N,IADR1)/APH(N)
                T2 = RC(N,IADR2)*PQ(N,JCMP)
                T3 = RC(N,IADR3)*RTP
                BAB11(N,IAB) = BAB11(N,IAB)
     &                         + ECD11(IBMAP(N),ICD,JCMP)*(T1-T2+T3)
              ENDDO
            ENDIF
C
C           IF E21(CD) IS BIG ENOUGH, CONTINUE TO ADD ALL OVERLAPS
            IF(ICD21(ICD,JCMP).EQ.1) THEN
              DO N=1,MAXN
                T1 = RC(N,IADR1)*0.5D0/APH(N)
                T2 = RC(N,IADR2)*PQ(N,JCMP)
                T3 = RC(N,IADR3)*RTP
                BAB21(N,IAB) = BAB21(N,IAB)
     &                       + ECD21(IBMAP(N),ICD,JCMP)*(T1-T2+T3)
              ENDDO
            ENDIF
C
C         END LOOP OVER INDICES {IA',IB',IC'} FOR CD
          ENDDO
C
C       END LOOP OVER CD INDEX {JX,JY,JZ}
        ENDDO
C
C     END LOOP OVER INDICES {IA,IB,IC} FOR AB
      ENDDO
C
C**********************************************************************C
C     GENERATE ALL POSSIBLE TWO-ELECTRON INTEGRALS FROM THE            C
C     EAB COEFFICIENTS AND THE G-ARRAYS                                C
C**********************************************************************C
C
C     CALCULATE PHASE FACTORS FOR MQN AND KQN COMBINATIONS (P1,P2,P3)
      P1 = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
      P2 = DFLOAT((-1)**((MQN(3)-MQN(4))/2))
      P1 =-(P1*DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2))))
      P2 =-(P2*DFLOAT((KQN(3)*KQN(4))/IABS(KQN(3)*KQN(4))))
      P3 = P1*P2
C
      IJ = (IBAS-1)*NBAS(2) + JBAS
C
C**********************************************************************C
C     INTEGRAL BATCH 1: ( - - || - - )                                 C
C**********************************************************************C
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY
      DO N=1,MAXN
        Q1(N) = DCMPLX(0.0D0,0.0D0)
        Q2(N) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {IA,IB,IC}
      DO IAB=1,NTUVAB
        IF(IAB11(IAB,ICMP).EQ.1) THEN
          DO N=1,MAXN
            Q1(N) = Q1(N) +      EAB11(IJ,IAB,ICMP)*DREAL(BAB11(N,IAB))
            Q2(N) = Q2(N) + CONE*EAB11(IJ,IAB,ICMP)*DIMAG(BAB11(N,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO N=1,MAXN
        RR(N,1 ) = RR(N,1 ) +    (Q1(N)+Q2(N))*PRE(N)
        RR(N,4 ) = RR(N,4 ) + P2*(Q1(N)-Q2(N))*PRE(N)
        RR(N,13) = P3*DCONJG(RR(N,4))
        RR(N,16) = P3*DCONJG(RR(N,1))
      ENDDO
C
C**********************************************************************C
C     INTEGRAL BATCH 2: ( - - || + - )                                 C
C**********************************************************************C
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY
      DO N=1,MAXN
        Q1(N) = DCMPLX(0.0D0,0.0D0)
        Q2(N) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {IA,IB,IC}
      DO IAB=1,NTUVAB
        IF(IAB11(IAB,ICMP).EQ.1) THEN
          DO N=1,MAXN
            Q1(N) = Q1(N) +      EAB11(IJ,IAB,ICMP)*DREAL(BAB21(N,IAB))
            Q2(N) = Q2(N) + CONE*EAB11(IJ,IAB,ICMP)*DIMAG(BAB21(N,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO N=1,MAXN
        RR(N,3 ) = RR(N,3 ) +    (Q1(N)+Q2(N))*PRE(N)
        RR(N,2 ) = RR(N,2 ) - P2*(Q1(N)-Q2(N))*PRE(N)
        RR(N,15) =-P3*DCONJG(RR(N,2))
        RR(N,14) =-P3*DCONJG(RR(N,3))
      ENDDO
C
C**********************************************************************C
C     INTEGRAL BATCH 3: ( + - || - - )                                 C
C**********************************************************************C
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY
      DO N=1,MAXN
        Q1(N) = DCMPLX(0.0D0,0.0D0)
        Q2(N) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {IA,IB,IC}
      DO IAB=1,NTUVAB
        IF(IAB21(IAB,ICMP).EQ.1) THEN
          DO N=1,MAXN
            Q1(N) = Q1(N) +      EAB21(IJ,IAB,ICMP)*DREAL(BAB11(N,IAB))
            Q2(N) = Q2(N) + CONE*EAB21(IJ,IAB,ICMP)*DIMAG(BAB11(N,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO N=1,MAXN
        RR(N,9 ) = RR(N,9 ) +    (Q1(N)+Q2(N))*PRE(N)
        RR(N,12) = RR(N,12) + P2*(Q1(N)-Q2(N))*PRE(N)
        RR(N,5 ) =-P3*DCONJG(RR(N,12))
        RR(N,8 ) =-P3*DCONJG(RR(N,9 ))
      ENDDO
C
C**********************************************************************C
C     INTEGRAL BATCH 4: ( + - || + - )                                 C
C**********************************************************************C
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY
      DO N=1,MAXN
        Q1(N) = DCMPLX(0.0D0,0.0D0)
        Q2(N) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {IA,IB,IC}
      DO IAB=1,NTUVAB
        IF(IAB21(IAB,ICMP).EQ.1) THEN
          DO N=1,MAXN
            Q1(N) = Q1(N) +      EAB21(IJ,IAB,ICMP)*DREAL(BAB21(N,IAB))
            Q2(N) = Q2(N) + CONE*EAB21(IJ,IAB,ICMP)*DIMAG(BAB21(N,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO N=1,MAXN
        RR(N,11) = RR(N,11) +    (Q1(N)+Q2(N))*PRE(N)
        RR(N,10) = RR(N,10) - P2*(Q1(N)-Q2(N))*PRE(N)
        RR(N,7 ) = P3*DCONJG(RR(N,10))
        RR(N,6 ) = P3*DCONJG(RR(N,11))
      ENDDO
C
C     END LOOP OVER INDICES {IX,IY,IZ}
6000  CONTINUE
C
C**********************************************************************C
C     BREIT OVERLAP ARRAY NOW FULLY CONSTRUCTED                        C
C**********************************************************************C
C
C     INCLUDE THE OUTSIDE FACTOR OF (1/2)
C     DFNOTE: GOT RID OF THIS WHEN I CHANGED ICNTA TERMINAL
      DO ITG=1,16
        DO N=1,MAXN
          RR(IBMAP(N),ITG) = 0.5D0*RR(IBMAP(N),ITG)
        ENDDO
      ENDDO
C
      RETURN
      END

