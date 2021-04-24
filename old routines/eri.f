      SUBROUTINE ERI(RR,XYZ,KQN,MQN,EXPT,NFUNS,ITQN,IBAS,JBAS,IEAB,IECD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                        EEEEEEEE RRRRRRR  IIII                        C
C                        EE       RR    RR  II                         C
C                        EE       RR    RR  II                         C
C                        EEEEEE   RR    RR  II                         C
C                        EE       RRRRRRR   II                         C
C                        EE       RR    RR  II                         C
C                        EEEEEEEE RR    RR IIII                        C
C                                                                      C
C -------------------------------------------------------------------- C
C  ERI GENERATES BLOCKS OF ELECTRON REPULSION INTEGRALS OVER           C
C  KINETICALLY BALANCED G-SPINOR BASIS FUNCTIONS. THE DENSITIES ARE    C
C  EXPANDED IN A BASIS OF HERMITE GAUSSIANS AND THE INTEGRALS ARE      C
C  GENERATED USING THE MCMURCHIE-DAVIDSON ALGORITHM.                   C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    XYZ(3,4)  - COORDINATES OF THE 4 NUCLEAR CENTERS.                 C
C    KQN(4)    - KQN QUANTUM NUMBERS OF THE CENTERS.                   C
C    MQN(4)    - |MQN| QUANTUM NUMBERS OF THE CENTERS.                 C
C    EXPT(I,4) - EXPONENTS IN BASIS FUNCTION BLOCK FOR A,B,C,D         C
C    NFUNS(J)  - NUMBER OF FUNCTIONS ON CENTER J.                      C
C    ITQN(2)   - SPINOR LABEL PAIR FOR AB (I=1) AND CD (I=2).          C
C                ITQN(I) = {LL,LS,SL,SS}.                              C
C    IBAS,JBAS - COMPONENT LABEL INDEX FOR BASIS FUNCTION PAIR ON AB . C
C    IEAB,IECD - 0 DON'T RECALCULATE E(AB/CD)-COEFFICIENTS             C
C                1 DO    RECALCULATE E(AB/CD)-COEFFICIENTS             C
C  OUTPUT:                                                             C
C    RR(MB2,16) - ERI'S FOR BLOCK AB, ALL 16 MQN PROJECTION COMBOS.    C
C -------------------------------------------------------------------- C
C  DFNOTE: THERE ARE NO SHORTCUTS FOR INUCAB OR INUCCD CONDITIONS.     C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 CONE
      COMPLEX*16 RR(MB2,16),Q1(MB2),Q2(MB2)
      COMPLEX*16 EAB11(MB2,MEQ),ECD11(MB2,MEQ),GAB11(MB2,MEQ),
     &           EAB21(MB2,MEQ),ECD21(MB2,MEQ),GAB21(MB2,MEQ)
C
      DIMENSION KQN(4),LQN(4),MQN(4),ITQN(2),NFUNS(4),EXPT(MBS,4)
      DIMENSION PQ(MB2,3),APH(MB2),PRE(MB2),RC(MB2,MRC),XYZ(3,4)
      DIMENSION IAB11(MEQ),IAB21(MEQ),ICD11(MEQ),ICD21(MEQ),IRC(MRC)
C
      SAVE EAB11,EAB21,ECD11,ECD21
      SAVE IAB11,IAB21,ICD11,ICD21
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &             IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,IEQS,IPAR,ICOR,ILEV
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
      DATA ROOTPI,SENS/1.7724538509055160D0,1.0D-10/
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     ILLEGAL COMPONENT OVERLAP CHECKER
      DO IT=1,2
        IF(ITQN(IT).NE.1.AND.ITQN(IT).NE.4) THEN
          WRITE(6, *) 'In ERI: incorrect call. ITQN(',IT,') = ',ITQN(IT)
          WRITE(7, *) 'In ERI: incorrect call. ITQN(',IT,') = ',ITQN(IT)
          STOP
        ENDIF
      ENDDO
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
      MAXAB = NFUNS(1)*NFUNS(2)
      MAXCD = NFUNS(3)*NFUNS(4)
C
C     PHASE FACTOR FOR AB AND CD PAIR OVERLAPS
      IPHSAB = 1
      IPHSCD =-1
C
C     VRS MAXIMUM LAMBDA LEVEL FOR EQ-COEFFICIENT ADDRESSES
      IF(ITQN(1).EQ.1) THEN
        LAMAB = LQN(1)+LQN(2)
      ELSEIF(ITQN(1).EQ.4) THEN
        LAMAB = LQN(1)+LQN(2)+2
      ENDIF
C
      IF(ITQN(2).EQ.1) THEN
        LAMCD = LQN(3)+LQN(4)
      ELSEIF(ITQN(2).EQ.4) THEN
        LAMCD = LQN(3)+LQN(4)+2
      ENDIF
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
      IF(IEAB.EQ.0) GOTO 100
C
C     GENERATE THE COEFFICIENT LIST
      CALL CPU_TIME(TDM1)
      IF(IEQS.EQ.0) THEN
C       CALCULATE FROM SCRATCH
        IF(ITQN(1).EQ.1) THEN
          CALL EMAKELL(EAB11,EAB21,EXPT,XYZ,KQN,MQN,NFUNS,IPHSAB,1,2,0)         
        ELSEIF(ITQN(1).EQ.4) THEN
          CALL EMAKESS(EAB11,EAB21,EXPT,XYZ,KQN,MQN,NFUNS,IPHSAB,1,2,0)         
        ENDIF
      ELSEIF(IEQS.EQ.1) THEN
C       READ FROM LOCAL EQ-COEFFICIENT FILE
        IF(ITQN(1).EQ.1) THEN
          DO ITUV=1,NTUVAB
            IAD = IABLL + (ITUV-1)*MAXAB
            DO M=1,MAXAB
              EAB11(M,ITUV) = DCMPLX(E0LLFL(IAD+M,1),E0LLFL(IAD+M,2))
              EAB21(M,ITUV) = DCMPLX(E0LLFL(IAD+M,3),E0LLFL(IAD+M,4))
            ENDDO       
          ENDDO       
        ELSEIF(ITQN(1).EQ.4) THEN
          DO ITUV=1,NTUVAB
            IAD = IABSS + (ITUV-1)*MAXAB
            DO M=1,MAXAB
              EAB11(M,ITUV) = DCMPLX(E0SSFL(IAD+M,1),E0SSFL(IAD+M,2))
              EAB21(M,ITUV) = DCMPLX(E0SSFL(IAD+M,3),E0SSFL(IAD+M,4))
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C
C     RECORD THE TIME TAKEN TO GENERATE/READ THE E(AB) COEFFICIENTS
      CALL CPU_TIME(TDM2)
      IF(ITQN(1).EQ.1) THEN
        TELL = TELL + TDM2 - TDM1       
      ELSEIF(ITQN(1).EQ.4) THEN
        TESS = TESS + TDM2 - TDM1
      ENDIF
C
C     SCREENING: CHECK E(AB) LIST MAGNITUDE BY IAB ENTRY IN FINITE SUM
      DO IAB=1,NTUVAB
C
C       11 OVERLAP (AB PAIRS)
        SUM = 0.0D0
        DO M=1,MAXAB
          SUM = SUM + ABS(EAB11(M,IAB))
        ENDDO
C
        IF(SUM.LE.SENS) THEN
          IAB11(IAB) = 0
        ELSE
          IAB11(IAB) = 1
        ENDIF
C
C       21 OVERLAP (AB PAIRS)
        SUM = 0.0D0
        DO M=1,MAXAB
          SUM = SUM + ABS(EAB21(M,IAB))
        ENDDO
C
        IF(SUM.LE.SENS) THEN
          IAB21(IAB) = 0
        ELSE
          IAB21(IAB) = 1
        ENDIF
        
      ENDDO
C
C     DO NOT CALCULATE AGAIN UNTIL ASKED EXTERNALLY
      IEAB = 0
C
100   CONTINUE
C
C**********************************************************************C
C     IF ASKED TO RECALCULATE E(CD) COEFFICIENTS, DO THIS NEXT         C
C**********************************************************************C
C
      IF(IECD.EQ.0) GOTO 200
C
C     GENERATE THE COEFFICIENT LIST
      CALL CPU_TIME(TDM1)
      IF(IEQS.EQ.0) THEN
C       CALCULATE FROM SCRATCH
        IF(ITQN(2).EQ.1) THEN
          CALL EMAKELL(ECD11,ECD21,EXPT,XYZ,KQN,MQN,NFUNS,IPHSCD,3,4,0)
        ELSEIF(ITQN(2).EQ.4) THEN
          CALL EMAKESS(ECD11,ECD21,EXPT,XYZ,KQN,MQN,NFUNS,IPHSCD,3,4,0)
        ENDIF
      ELSEIF(IEQS.EQ.1) THEN
C       READ FROM LOCAL EQ-COEFFICIENT FILE
        IF(ITQN(2).EQ.1) THEN
          DO ITUV=1,NTUVCD
            IAD = ICDLL + (ITUV-1)*MAXCD
            DO M=1,MAXCD
              ECD11(M,ITUV) = DCMPLX(E0LLFL(IAD+M,5),E0LLFL(IAD+M,6))
              ECD21(M,ITUV) = DCMPLX(E0LLFL(IAD+M,7),E0LLFL(IAD+M,8))
            ENDDO       
          ENDDO       
        ELSEIF(ITQN(2).EQ.4) THEN
          DO ITUV=1,NTUVCD
            IAD = ICDSS + (ITUV-1)*MAXCD
            DO M=1,MAXCD
              ECD11(M,ITUV) = DCMPLX(E0SSFL(IAD+M,5),E0SSFL(IAD+M,6))
              ECD21(M,ITUV) = DCMPLX(E0SSFL(IAD+M,7),E0SSFL(IAD+M,8))
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C
C     RECORD THE TIME TAKEN TO GENERATE/READ THE E(CD) COEFFICIENTS
      CALL CPU_TIME(TDM2)
      IF(ITQN(2).EQ.1) THEN
        TELL = TELL + TDM2 - TDM1       
      ELSEIF(ITQN(2).EQ.4) THEN
        TESS = TESS + TDM2 - TDM1
      ENDIF

C
C     SCREENING: CHECK E(CD) LIST MAGNITUDE BY ICD ENTRY IN FINITE SUM
      DO ICD=1,NTUVCD
C
C       11 OVERLAP (CD PAIRS)
        SUM = 0.0D0
        DO M=1,MAXCD
          SUM = SUM + ABS(ECD11(M,ICD))
        ENDDO
C
        IF(SUM.LE.SENS) THEN
          ICD11(ICD) = 0
        ELSE
          ICD11(ICD) = 1
        ENDIF
C
C       21 OVERLAP (CD PAIRS)
        SUM = 0.0D0
        DO M=1,MAXCD
          SUM = SUM + ABS(ECD21(M,ICD))
        ENDDO
C
        IF(SUM.LE.SENS) THEN
          ICD21(ICD) = 0
        ELSE
          ICD21(ICD) = 1
        ENDIF
C
      ENDDO
C
C     DO NOT CALCULATE AGAIN UNTIL ASKED EXTERNALLY
      IECD = 0
C
200   CONTINUE
C
C**********************************************************************C
C     HERMITE TWO-BODY INTEGRAL PREPARATION AND EVALUATION - R(ABCD)   C
C**********************************************************************C
C
C     GAUSSIAN OVERLAP CENTER
      EIJ = EXPT(IBAS,1) + EXPT(JBAS,2)
      PX  = (XYZ(1,1)*EXPT(IBAS,1) + XYZ(1,2)*EXPT(JBAS,2))/EIJ
      PY  = (XYZ(2,1)*EXPT(IBAS,1) + XYZ(2,2)*EXPT(JBAS,2))/EIJ
      PZ  = (XYZ(3,1)*EXPT(IBAS,1) + XYZ(3,2)*EXPT(JBAS,2))/EIJ
C
C     AUXILLIARY DATA FOR RMAKE ROUTINE
      M = 0
      DO KBAS=1,NFUNS(3)
        DO LBAS=1,NFUNS(4)
          M = M+1         
          EKL = EXPT(KBAS,3) + EXPT(LBAS,4)
          QX  = (XYZ(1,3)*EXPT(KBAS,3) + XYZ(1,4)*EXPT(LBAS,4))/EKL
          QY  = (XYZ(2,3)*EXPT(KBAS,3) + XYZ(2,4)*EXPT(LBAS,4))/EKL
          QZ  = (XYZ(3,3)*EXPT(KBAS,3) + XYZ(3,4)*EXPT(LBAS,4))/EKL
          EIJKL   = EIJ+EKL
          APH(M)  = EIJ*EKL/EIJKL
          PQ(M,1) = QX-PX
          PQ(M,2) = QY-PY
          PQ(M,3) = QZ-PZ
          EMX     = DSQRT(EIJ+EKL)*EIJ*EKL
          PRE(M)  = 2.0D0*(ROOTPI**5)/EMX
        ENDDO
      ENDDO
C
C     GENERATE R-INTEGRALS
      CALL CPU_TIME (TDM1)
      CALL RMAKE(RC,PQ,APH,MAXCD,LAMABCD)
C
C     RECORD THE TIME TAKEN TO GENERATE THE R(ABCD) INTEGRALS
      CALL CPU_TIME(TDM2)
      IF(ITQN(1).EQ.1.AND.ITQN(2).EQ.1) THEN
        TRLL = TRLL + TDM2 - TDM1
      ELSEIF(ITQN(1).EQ.4.AND.ITQN(2).EQ.4) THEN
        TRSS = TRSS + TDM2 - TDM1
      ELSE
        TRLS = TRLS + TDM2 - TDM1
      ENDIF
C
C     SCREENING PROCEDURE: R(ABCD) LIST MAGNITUDE BY IABCD INDEX
      DO IABCD=1,NTUVABCD
C
C       VECTOR MAGNITUDE OF FULL R-INTEGRAL LIST
        SUM = 0.0D0
        DO M=1,MAXCD
          SUM = SUM + DABS(RC(M,IABCD))
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
C**********************************************************************C
C     CONSTRUCT INTERMEDIATE MATRICES FOR MCMURCHIE-DAVIDSON           C
C**********************************************************************C
C
C     LOOP OVER ALL ADDRESSES FOR E(AB) FINITE EXPANSION
      DO IAB=1,NTUVAB
C
C       INITIALISE THE INTERMEDIATE ERI COMPONENTS G(AB)
        DO M=1,MAXCD
          GAB11(M,IAB) = DCMPLX(0.0D0,0.0D0)
          GAB21(M,IAB) = DCMPLX(0.0D0,0.0D0)
        ENDDO
C
C       SKIP ENTIRE PROCESS IF E(AB) HAS BEEN SCREENED
        IF(IAB11(IAB)+IAB21(IAB).EQ.0) GOTO 799
C
C       LOOP OVER ALL FINITE EXPANSION ADDRESSES FOR E(CD)
        DO ICD=1,NTUVCD
C
C         CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
          IRABCD = INABCD(IVEC(IAB)+IVEC(ICD),JVEC(IAB)+JVEC(ICD),
     &                                        KVEC(IAB)+KVEC(ICD))
C
C         SKIP THIS STEP IF THE R-INTEGRAL LIST FOR THIS AB/CD SCREENED
          IF(IRC(IRABCD).EQ.0) GOTO 798

          IABC = IVEC(ICD) + JVEC(ICD) + KVEC(ICD)
          IF(MOD(IABC,2).EQ.2) THEN
            PR =-1.0D0
          ELSE
            PR = 1.0D0
          ENDIF
C
C         CONTRIBUTIONS TO G(AB|--) FROM EACH E(CD|--) ADDRESS
          IF(ICD11(ICD).EQ.1) THEN
            DO M=1,MAXCD
              GAB11(M,IAB) = GAB11(M,IAB) + ECD11(M,ICD)*RC(M,IRABCD)
            ENDDO
          ENDIF
C
C         CONTRIBUTIONS TO G(AB|+-) FROM EACH E(CD|+-) ADDRESS
          IF(ICD21(ICD).EQ.1) THEN
            DO M=1,MAXCD
              GAB21(M,IAB) = GAB21(M,IAB) + ECD21(M,ICD)*RC(M,IRABCD)
            ENDDO
          ENDIF
C
C         SKIP POINT FOR R(ABCD) SCREENING
798       CONTINUE
C
C       END LOOP OVER E(CD) FINITE EXPANSION ADDRESSES
        ENDDO
C
C       SKIP POINT FOR E(AB) SCREENING
799     CONTINUE
C
C     END LOOP OVER E(AB) FINITE EXPANSION ADDRESSES
      ENDDO
C
C**********************************************************************C
C     GENERATE ALL POSSIBLE TWO-ELECTRON INTEGRALS FROM THE            C
C     EAB COEFFICIENTS AND THE G-ARRAYS                                C
C**********************************************************************C
C
C     CALCULATE PHASES FOR BASIS FUNCTION OVERLAP COMBINATIONS
      PAB = ISIGN(1,KQN(1)*KQN(2))*(-1)**((MQN(1)-MQN(2))/2)
      PCD = ISIGN(1,KQN(3)*KQN(4))*(-1)**((MQN(3)-MQN(4))/2)
      PABCD = PAB*PCD
C
      IJ = (IBAS-1)*NFUNS(2) + JBAS
C
C**********************************************************************C
C     INTEGRAL BATCH 1: ( - - || - - )                                 C
C**********************************************************************C
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY
      DO M=1,MAXCD
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     CONTRIBUTIONS TO RR(--||--) FROM EACH E(AB|--) ADDRESS
      DO IAB=1,NTUVAB
        IF(IAB11(IAB).EQ.1) THEN
          DO M=1,MAXCD
            Q1(M) = Q1(M) +      EAB11(IJ,IAB)*DREAL(GAB11(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB11(IJ,IAB)*DIMAG(GAB11(M,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO M=1,MAXCD
        RR(M,1 ) =     (Q1(M)+Q2(M))*PRE(M)
        RR(M,4 ) = PCD*(Q1(M)-Q2(M))*PRE(M)
        RR(M,13) = PABCD*DCONJG(RR(M,4))
        RR(M,16) = PABCD*DCONJG(RR(M,1))
      ENDDO
C
C**********************************************************************C
C     INTEGRAL BATCH 2: ( - - || + - )                                 C
C**********************************************************************C
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY
      DO M=1,MAXCD
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     CONTRIBUTIONS TO RR(--||+-) FROM EACH E(AB|--) ADDRESS
      DO IAB=1,NTUVAB
        IF(IAB11(IAB).EQ.1) THEN
          DO M=1,MAXCD
            Q1(M) = Q1(M) +      EAB11(IJ,IAB)*DREAL(GAB21(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB11(IJ,IAB)*DIMAG(GAB21(M,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO M=1,MAXCD
        RR(M,3 ) =     (Q1(M)+Q2(M))*PRE(M)
        RR(M,2 ) =-PCD*(Q1(M)-Q2(M))*PRE(M)
        RR(M,15) =-PABCD*DCONJG(RR(M,2))
        RR(M,14) =-PABCD*DCONJG(RR(M,3))
      ENDDO
C
C**********************************************************************C
C     INTEGRAL BATCH 3: ( + - || - - )                                 C
C**********************************************************************C
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY
      DO M=1,MAXCD
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     CONTRIBUTIONS TO RR(+-||--) FROM EACH E(AB|+-) ADDRESS
      DO IAB=1,NTUVAB
       IF(IAB21(IAB).EQ.1) THEN
          DO M=1,MAXCD
            Q1(M) = Q1(M) +      EAB21(IJ,IAB)*DREAL(GAB11(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB21(IJ,IAB)*DIMAG(GAB11(M,IAB))
          ENDDO
       ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO M=1,MAXCD
        RR(M,9 ) =     (Q1(M)+Q2(M))*PRE(M)
        RR(M,12) = PCD*(Q1(M)-Q2(M))*PRE(M)
        RR(M,5 ) =-PABCD*DCONJG(RR(M,12))
        RR(M,8 ) =-PABCD*DCONJG(RR(M,9 ))
      ENDDO
C
C**********************************************************************C
C     INTEGRAL BATCH 4: ( + - || + - )                                 C
C**********************************************************************C
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY
      DO M=1,MAXCD
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     CONTRIBUTIONS TO RR(+-||+-) FROM EACH E(AB|+-) ADDRESS
      DO IAB=1,NTUVAB
        IF(IAB21(IAB).EQ.1) THEN
          DO M=1,MAXCD
            Q1(M) = Q1(M) +      EAB21(IJ,IAB)*DREAL(GAB21(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB21(IJ,IAB)*DIMAG(GAB21(M,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO M=1,MAXCD
        RR(M,11) =     (Q1(M)+Q2(M))*PRE(M)
        RR(M,10) =-PCD*(Q1(M)-Q2(M))*PRE(M)
        RR(M,7 ) = PABCD*DCONJG(RR(M,10))
        RR(M,6 ) = PABCD*DCONJG(RR(M,11))
      ENDDO
C
      RETURN
      END

