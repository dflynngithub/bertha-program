      SUBROUTINE BII(RR,XYZ,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN,BGAUNT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                          BBBBBBB IIII IIII                           C
C                          BB    BB II   II                            C
C                          BB    BB II   II                            C
C                          BBBBBBB  II   II                            C
C                          BB    BB II   II                            C
C                          BB    BB II   II                            C
C                          BBBBBBB IIII IIII                           C
C                                                                      C
C -------------------------------------------------------------------- C
C  ERI GENERATES A BATCH OF MOLECULAR ELECTRON REPULSION INTEGRALS BY  C
C  MEANS OF THE MCMURCHIE-DAVIDSION ALGORITHM (DOUBLE FINITE SUM OVER  C
C  EQ-COEFFICIENTS AND INTEGRALS OVER A PAIR OF HGTFS.)                C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ XYZ  - FULL SET OF CARTESIAN BASIS CENTERS.                       C
C  ▶ KQN  - FULL SET OF RELATIVISTIC LABELS.                           C
C  ▶ MQN  - FULL SET OF MAGNETIC QUANTUM NUMBERS (MAGNITUDE).          C
C  ▶ NBAS - FULL SET OF EXPONENT LIST LENGTHS.                         C
C  ▶ EXL  - FULL LISTS OF EXPONENTS IN THE BLOCK.                      C
C  ▶ IBAS - 1ST BASIS FUNCTION (HELD CONSTANT DURING ROUTINE).         C
C  ▶ JBAS - 2ND BASIS FUNCTION (HELD CONSTANT DURING ROUTINE).         C
C  ▶ ITN  - COMPONENT OVERLAP COMBINATION.                             C
C  ▶ MGNT - GAUNT INTERACTION OVERRIDE OPTION.                         C
C  OUTPUT:                                                             C
C  ▶ RR   - BII'S FOR BLOCK AB, ALL 16 MQN SIGN COMBINATIONS.          C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL BGAUNT
C
      DIMENSION ITN(2)
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
      DIMENSION PQ(MB2,3),APH(MB2),PRE(MB2),T(MB2),RC(MB2,MRC)
      DIMENSION IABR11(MEQ,3),IABI11(MEQ,3),IABR21(MEQ,3),IABI21(MEQ,3)
      DIMENSION ICDR11(MEQ,3),ICDI11(MEQ,3),ICDR21(MEQ,3),ICDI21(MEQ,3)
      DIMENSION IRC(MRC)
      DIMENSION RCTTFL(20*MFL),IRCTTFL(MFL)
      DIMENSION GABR11(MB2,MEQ),GABI11(MB2,MEQ),
     &          GABR21(MB2,MEQ),GABI21(MB2,MEQ)
      DIMENSION QR1(MB2),QI1(MB2),QR2(MB2),QI2(MB2)
      DIMENSION IDX(3),JDX(3)
C
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 EAB11(MB2,MEQ,3),EAB21(MB2,MEQ,3),
     &           ECD11(MB2,MEQ,3),ECD21(MB2,MEQ,3)
C
      SAVE EAB11,EAB21,IABR11,IABI11,IABR21,IABI21
      SAVE ECD11,ECD21,ICDR11,ICDI11,ICD221,ICDI21
      SAVE RCTTFL,IRCTTFL
C
      COMMON/EILS/EILSFL(MFL,12),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EISL/EISLFL(MFL,12),IADISL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS,IABSL,ICDSL
      COMMON/IRCM/IEAB,IECD,IRIJ(MBS,MBS)
      COMMON/ISCR/IMTX(MB2,11),ISCR(MB2),IMAP(MB2),IBCH,ITOG,MAXN
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/TMMD/TELL,TESS,TELS,TESL,TRLL,TRSS,TRLS,TRSL,TRBR
      COMMON/TSCF/TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRR,TCC1,TCC2,TCMC,
     &            TB1B,TB1R,TB1F,TB1M,TBEC,TBRM,TBRR,TBC1,TBC2,TBMC,
     &            TQMX,THMX,TC1T,TC2T,TCVT,TB1T,TB2T,TEIG,TSCR,TTOT,
     &            TC1S,TC2S,TB1S,TB2S
C
C     EQ-COEFFICIENT SENSITIVITY PARAMETER
      DATA SENS/1.0D-10/
C
C     ILLEGAL COMPONENT OVERLAP CHECKER
      DO IT=1,2
        IF(ITN(IT).NE.2.AND.ITN(IT).NE.3) THEN
          WRITE(6, *) 'In BII: illegal component overlaps in ITN.'
          WRITE(7, *) 'In BII: illegal component overlaps in ITN.'
          STOP
        ENDIF
      ENDDO
C
C     EVALUATE LQNS FOR BASIS FUNCTIONS (A,B,C,D)
      DO N=1,4
        IF(KQN(N).LT.0) THEN
         LQN(N) =-KQN(N)-1
        ELSE
         LQN(N) = KQN(N)
        ENDIF
      ENDDO
C
C     INTEGRAL SKIPPING ON MOLECULAR GROUP SYMMETRY CLASS BASIS
      IF(SHAPE.EQ.'ATOMIC'.OR.SHAPE.EQ.'DIATOMIC'.OR.
     &                        SHAPE.EQ.'LINEAR') THEN
        ISYM = 1
      ELSE
        ISYM = 0
      ENDIF
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
C     VRS MAXIMUM LAMBDA LEVEL FOR CONTRACTED R-INTEGRAL BATCH
      IF(BGAUNT) THEN
C       GAUNT INTEGRALS ONLY
        LAMABCD = LAMAB+LAMCD
      ELSE
C       FULL BREIT INTERACTION
        LAMABCD = LAMAB+LAMCD+2
      ENDIF
C
C     VRS TOTAL LENGTH OF EQ-COEFFICIENT LISTS AND R-INTEGRAL BATCH
      NTUVAB   = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
      NTUVCD   = (LAMCD+1)*(LAMCD+2)*(LAMCD+3)/6
      NTUVABCD = (LAMABCD+1)*(LAMABCD+2)*(LAMABCD+3)/6
C
C     LIST ADDRESS FOR (AB|  ) AND GAUSSIAN EXPONENT FOR AB OVERLAP
      IJ  = (IBAS-1)*NBAS(2)+JBAS
      EIJ = EXL(IBAS,1)+EXL(JBAS,2)
C
C     INITIALISE RR ARRAY
      DO M=1,MAXCD
        DO ITG=1,16
          RR(M,ITG) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     GENERATE NEW BATCH OF E(AB|  ) COEFFICIENTS IF PROMPTED          C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
      IF(IEAB.EQ.0) GOTO 100
C
C     START TIME
      CALL CPU_TIME(TDM1)
C
C     OPTION 1: READ FROM LOCAL EQ-COEFFICIENT FILE
      IF(EQFILE) THEN
        IF(ITN(1).EQ.2) THEN
          DO ITUV=1,NTUVAB
            IAD = IABLS + (ITUV-1)*MAXAB
            DO M=1,MAXAB
            EAB11(M,ITUV,1) = DCMPLX(EILSFL(IAD+M, 1),EILSFL(IAD+M, 2))
            EAB21(M,ITUV,1) = DCMPLX(EILSFL(IAD+M, 3),EILSFL(IAD+M, 4))
            EAB11(M,ITUV,2) = DCMPLX(EILSFL(IAD+M, 5),EILSFL(IAD+M, 6))
            EAB21(M,ITUV,2) = DCMPLX(EILSFL(IAD+M, 7),EILSFL(IAD+M, 8))
            EAB11(M,ITUV,3) = DCMPLX(EILSFL(IAD+M, 9),EILSFL(IAD+M,10))
            EAB21(M,ITUV,3) = DCMPLX(EILSFL(IAD+M,11),EILSFL(IAD+M,12))
            ENDDO
          ENDDO
        ELSEIF(ITN(1).EQ.3) THEN
          DO ITUV=1,NTUVAB
            IAD = IABSL + (ITUV-1)*MAXAB
            DO M=1,MAXAB
            EAB11(M,ITUV,1) = DCMPLX(EISLFL(IAD+M, 1),EISLFL(IAD+M, 2))
            EAB21(M,ITUV,1) = DCMPLX(EISLFL(IAD+M, 3),EISLFL(IAD+M, 4))
            EAB11(M,ITUV,2) = DCMPLX(EISLFL(IAD+M, 5),EISLFL(IAD+M, 6))
            EAB21(M,ITUV,2) = DCMPLX(EISLFL(IAD+M, 7),EISLFL(IAD+M, 8))
            EAB11(M,ITUV,3) = DCMPLX(EISLFL(IAD+M, 9),EISLFL(IAD+M,10))
            EAB21(M,ITUV,3) = DCMPLX(EISLFL(IAD+M,11),EISLFL(IAD+M,12))
            ENDDO
          ENDDO
        ENDIF
C
C     OPTION 2: CALCULATE FROM SCRATCH
      ELSE
        IF(ITN(1).EQ.2) THEN
          CALL EILSB3(EAB11,EAB21,EXL,XYZ,KQN,MQN,NBAS,IPHSAB,1,2)
        ELSEIF(ITN(1).EQ.3) THEN
          CALL EISLB3(EAB11,EAB21,EXL,XYZ,KQN,MQN,NBAS,IPHSAB,1,2)
        ENDIF
      ENDIF
C
C     RECORD THE TIME TAKEN TO GENERATE/READ THE E(AB|  ) COEFFICIENTS
      CALL CPU_TIME(TDM2)
      IF(ITN(1).EQ.2) THEN
        TELS = TELS+TDM2-TDM1
      ELSEIF(ITN(1).EQ.3) THEN
        TESL = TESL+TDM2-TDM1
      ENDIF
C
C     SCREENING PROCEDURE: NORM SUM OF EQ-COEFFICIENT LIST FOR EACH IAB
      DO ICMP=1,3
        DO IAB=1,NTUVAB
C
C         Re{E(AB|--)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXAB
            ER = DREAL(EAB11(M,IAB,ICMP))
            SUM = SUM + DABS(ER)
            IF(SUM.GT.SENS) THEN
              IABR11(IAB,ICMP) = 1
              GOTO 101
            ENDIF
          ENDDO
          IABR11(IAB,ICMP) = 0
101       CONTINUE
C
C         Im{E(AB|--)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXAB
            EI = DIMAG(EAB11(M,IAB,ICMP))
            SUM = SUM + DABS(EI)
            IF(SUM.GT.SENS) THEN
              IABI11(IAB,ICMP) = 1
              GOTO 102
            ENDIF
          ENDDO
          IABI11(IAB,ICMP) = 0
102       CONTINUE
C
C         Re{E(AB|+-)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXAB
            ER = DREAL(EAB21(M,IAB,ICMP))
            SUM = SUM + DABS(ER)
            IF(SUM.GT.SENS) THEN
              IABR21(IAB,ICMP) = 1
              GOTO 103
            ENDIF
          ENDDO
          IABR21(IAB,ICMP) = 0
103       CONTINUE
C
C         Im{E(AB|+-)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXAB
            EI = DIMAG(EAB21(M,IAB,ICMP))
            SUM = SUM + DABS(EI)
            IF(SUM.GT.SENS) THEN
              IABI21(IAB,ICMP) = 1
              GOTO 104
            ENDIF
          ENDDO
          IABI21(IAB,ICMP) = 0
104       CONTINUE
C
        ENDDO
      ENDDO
C
C     DO NOT CALCULATE AGAIN UNTIL PROMPTED EXTERNALLY
      IEAB = 0
C
100   CONTINUE
C
C**********************************************************************C
C     GENERATE NEW BATCH OF E(CD| -) COEFFICIENTS IF PROMPTED          C
C**********************************************************************C
C
      IF(IECD.EQ.0) GOTO 200
C
C     START TIME
      CALL CPU_TIME(TDM1)
C
C     OPTION 1: READ FROM LOCAL EQ-COEFFICIENT FILE
      IF(EQFILE) THEN
        IF(ITN(2).EQ.2) THEN
          DO ITUV=1,NTUVCD
            IAD = ICDLS + (ITUV-1)*MAXCD
            Z   = DFLOAT((-1)**(ILAM(ITUV)))
            DO M=1,MAXCD
            ECD11(M,ITUV,1)=Z*DCMPLX(EILSFL(IAD+M, 1),EILSFL(IAD+M, 2))
            ECD21(M,ITUV,1)=Z*DCMPLX(EILSFL(IAD+M, 3),EILSFL(IAD+M, 4))
            ECD11(M,ITUV,2)=Z*DCMPLX(EILSFL(IAD+M, 5),EILSFL(IAD+M, 6))
            ECD21(M,ITUV,2)=Z*DCMPLX(EILSFL(IAD+M, 7),EILSFL(IAD+M, 8))
            ECD11(M,ITUV,3)=Z*DCMPLX(EILSFL(IAD+M, 9),EILSFL(IAD+M,10))
            ECD21(M,ITUV,3)=Z*DCMPLX(EILSFL(IAD+M,11),EILSFL(IAD+M,12))
            ENDDO
          ENDDO
        ELSEIF(ITN(2).EQ.3) THEN
          DO ITUV=1,NTUVCD
            IAD = ICDSL + (ITUV-1)*MAXCD
            Z   = DFLOAT((-1)**(ILAM(ITUV)))
            DO M=1,MAXCD
            ECD11(M,ITUV,1)=Z*DCMPLX(EISLFL(IAD+M, 1),EISLFL(IAD+M, 2))
            ECD21(M,ITUV,1)=Z*DCMPLX(EISLFL(IAD+M, 3),EISLFL(IAD+M, 4))
            ECD11(M,ITUV,2)=Z*DCMPLX(EISLFL(IAD+M, 5),EISLFL(IAD+M, 6))
            ECD21(M,ITUV,2)=Z*DCMPLX(EISLFL(IAD+M, 7),EISLFL(IAD+M, 8))
            ECD11(M,ITUV,3)=Z*DCMPLX(EISLFL(IAD+M, 9),EISLFL(IAD+M,10))
            ECD21(M,ITUV,3)=Z*DCMPLX(EISLFL(IAD+M,11),EISLFL(IAD+M,12))
            ENDDO
          ENDDO
        ENDIF
C
C     OPTION 2: CALCULATE FROM SCRATCH
      ELSE
        IF(ITN(2).EQ.2) THEN
          CALL EILSB3(ECD11,ECD21,EXL,XYZ,KQN,MQN,NBAS,IPHSCD,3,4)
        ELSEIF(ITN(2).EQ.3) THEN
          CALL EISLB3(ECD11,ECD21,EXL,XYZ,KQN,MQN,NBAS,IPHSCD,3,4)
        ENDIF
      ENDIF
C
C     RECORD THE TIME TAKEN TO GENERATE/READ THE E(CD|  ) COEFFICIENTS
      CALL CPU_TIME(TDM2)
      IF(ITN(2).EQ.2) THEN
        TELS = TELS+TDM2-TDM1
      ELSEIF(ITN(2).EQ.3) THEN
        TESL = TESL+TDM2-TDM1
      ENDIF
C
C     SCREENING PROCEDURE: NORM SUM OF EQ-COEFFICIENT LIST FOR EACH ICD
      DO JCMP=1,3
        DO ICD=1,NTUVCD
C
C         Re{E(CD|--)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXCD
            ER = DREAL(ECD11(M,ICD,JCMP))
            SUM = SUM + DABS(ER)
            IF(SUM.GT.SENS) THEN
              ICDR11(ICD,JCMP) = 1
              GOTO 201
            ENDIF
          ENDDO
          ICDR11(ICD,JCMP) = 0
201       CONTINUE
C
C         Im{E(CD|--)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXCD
            EI = DIMAG(ECD11(M,ICD,JCMP))
            SUM = SUM + DABS(EI)
            IF(SUM.GT.SENS) THEN
              ICDI11(ICD,JCMP) = 1
              GOTO 202
            ENDIF
          ENDDO
          ICDI11(ICD,JCMP) = 0
202       CONTINUE
C
C         Re{E(CD|+-)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXCD
            ER = DREAL(ECD21(M,ICD,JCMP))
            SUM = SUM + DABS(ER)
            IF(SUM.GT.SENS) THEN
              ICDR21(ICD,JCMP) = 1
              GOTO 203
            ENDIF
          ENDDO
          ICDR21(ICD,JCMP) = 0
203       CONTINUE
C
C         Im{E(CD|+-)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXCD
            EI = DIMAG(ECD21(M,ICD,JCMP))
            SUM = SUM + DABS(EI)
            IF(SUM.GT.SENS) THEN
              ICDI21(ICD,JCMP) = 1
              GOTO 204
            ENDIF
          ENDDO
          ICDI21(ICD,JCMP) = 0
204       CONTINUE
C
        ENDDO
      ENDDO
C
C     DO NOT CALCULATE AGAIN UNTIL ASKED EXTERNALLY
      IECD = 0
C
200   CONTINUE
      CALL CPU_TIME(T2)
      TBEC = TBEC + T2 - T1
C
C**********************************************************************C
C     GENERATE NEW BATCH OF RC(AB|CD) FROM SCRATCH OR READ-IN          C
C**********************************************************************C
C
C     FACTORS NEEDED IN BOTH CASES
C
C     GAUSSIAN OVERLAP CENTRE
      PX = (XYZ(1,1)*EXL(IBAS,1)+XYZ(1,2)*EXL(JBAS,2))/EIJ
      PY = (XYZ(2,1)*EXL(IBAS,1)+XYZ(2,2)*EXL(JBAS,2))/EIJ
      PZ = (XYZ(3,1)*EXL(IBAS,1)+XYZ(3,2)*EXL(JBAS,2))/EIJ
C
C     AUXILLIARY DATA FOR RMAKE ROUTINE
      M = 0
      N = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
          IF((RCFILE.AND.IRIJ(IBAS,JBAS).EQ.0.AND.ISCR(M).EQ.1)
     &                    .OR.(RCFILE.AND.IRIJ(IBAS,JBAS).EQ.1)
     &                       .OR.(.NOT.RCFILE.AND.ISCR(M).EQ.1)) THEN
            N   = N+1
            EKL = EXL(KBAS,3)+EXL(LBAS,4)
            QX  = (XYZ(1,3)*EXL(KBAS,3)+XYZ(1,4)*EXL(LBAS,4))/EKL
            QY  = (XYZ(2,3)*EXL(KBAS,3)+XYZ(2,4)*EXL(LBAS,4))/EKL
            QZ  = (XYZ(3,3)*EXL(KBAS,3)+XYZ(3,4)*EXL(LBAS,4))/EKL
            APH(N)  = EIJ*EKL/(EIJ+EKL)
            PQ(N,1) = QX-PX
            PQ(N,2) = QY-PY
            PQ(N,3) = QZ-PZ
            EMX     = DSQRT(EIJ+EKL)*EIJ*EKL
            PRE(N)  = 2.0D0*PI52/EMX
          ENDIF
        ENDDO
      ENDDO
C
C     START TIME
      CALL CPU_TIME(TDM1)
C
C     READ FROM LOCAL RC(AB|CD) FILE
      IF(RCFILE.AND.IRIJ(IBAS,JBAS).EQ.0) THEN
C
        CALL CPU_TIME(T1)
C
C       STARTING ADDRESS FOR SAVED R(AB|CD) INTEGRALS
        IADRTT = (IJ-1)*NBAS(3)*NBAS(4)*NTUVABCD
C
C       READ RC(AB|CD) INTEGRALS FROM THIS STARTING POINT
        N = 0
        DO N=1,MAXN
          DO IABCD=1,NTUVABCD
            IAD = IADRTT + (IMAP(N)-1)*NTUVABCD + IABCD
            RC(N,IABCD) = RCTTFL(IAD)
          ENDDO
        ENDDO
C
C       STARTING ADDRESS FOR SCREENING FLAGS
        IADSCR = (IJ-1)*NTUVABCD
C
C       READ SCREENING FLAGS FROM THIS STARTING POINT
        DO IABCD=1,NTUVABCD
          IAD = IADSCR + IABCD
          IRC(IABCD) = IRCTTFL(IAD)
        ENDDO
C
C       RECORD TIME SPENT READING R-INTEGRALS
        CALL CPU_TIME(T2)
        TBRR = TBRR+T2-T1
C
      ENDIF
C
C     CALCULATE FROM SCRATCH
      IF(.NOT.RCFILE.OR.IRIJ(IBAS,JBAS).EQ.1) THEN
C
        CALL CPU_TIME(T1)
C
C       GENERATE R-INTEGRALS
        CALL RMAKE(RC,PQ,APH,MAXN,LAMABCD)
C
C       SCREENING: TEST RC(AB|CD) COLUMNS WITH INDEX (T+T',U+U',V+V')
        DO IABCD=1,NTUVABCD
C
C         SUM OF RC(AB|CD) MAGNITUDES
          SUM = 0.0D0
          DO N=1,MAXN
            SUM = SUM + DABS(RC(N,IABCD))
            IF(SUM.GT.SENS) THEN
              IRC(IABCD) = 1
              GOTO 301
            ENDIF
          ENDDO
          IRC(IABCD) = 0
301       CONTINUE
C
        ENDDO
C
C       SAVE THIS SET TO APPROPRIATE CLASS ADDRESS
        IF(RCFILE) THEN
C
C         TEST WHETHER FINAL ADDRESS IS STILL INSIDE ARRAY BOUNDS
          ILIM = IJ*NBAS(3)*NBAS(4)*NTUVABCD
C
          IF(ILIM.GT.20*MFL) THEN
C           OUT OF BOUNDS: PRINT WARNING BUT KEEP GOING
            WRITE(6, *) 'In BII: RCTT words exceed allocated limit.'
            WRITE(7, *) 'In BII: RCTT words exceed allocated limit.'
            GOTO 300
          ELSE
C           DO NOT CALCULATE AGAIN UNTIL PROMPTED EXTERNALLY
            IRIJ(IBAS,JBAS) = 0
          ENDIF
C
C         STARTING ADDRESS FOR SAVED R(AB|CD) INTEGRALS
          IADRTT = (IJ-1)*NBAS(3)*NBAS(4)*NTUVABCD

C         COPY THIS BATCH OF INTEGRALS TO A SAVED LIST
          DO N=1,MAXN
            DO IABCD=1,NTUVABCD
              IAD = IADRTT + (N-1)*NTUVABCD + IABCD
              RCTTFL(IAD) = RC(N,IABCD)
            ENDDO
          ENDDO
C
C         STARTING ADDRESS FOR SCREENING FLAGS
          IADSCR = (IJ-1)*NTUVABCD
C
C         COPY SCREENING MARKERS TO A SAVED LIST
          DO IABCD=1,NTUVABCD
            IAD = IADSCR + IABCD
            IRCTTFL(IAD) = IRC(IABCD)
          ENDDO
C
C         SHORTEN THE CURRENT RC LIST WITH IMAP FROM SCREENING
          M = 0
          N = 0
          DO KBAS=1,NBAS(3)
            DO LBAS=1,NBAS(4)
              M = M+1
              IF(ISCR(M).EQ.1) THEN
                N       = N+1
                APH(N)  = APH(M)
                PQ(N,1) = PQ(M,1)
                PQ(N,2) = PQ(M,2)
                PQ(N,3) = PQ(M,3)
                PRE(N)  = PRE(M)
                DO IABCD=1,NTUVABCD
                  RC(N,IABCD) = RC(M,IABCD)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
C
C         SHORTEN MAXN
          MAXN = N
C
        ENDIF
300     CONTINUE
C
        CALL CPU_TIME(T2)
        TBRM = TBRM+T2-T1
C
      ENDIF
C
C     RECORD THE TIME TAKEN TO GENERATE THE RC(AB|CD) BATCH
      CALL CPU_TIME(TDM2)
      TRBR = TRBR+TDM2-TDM1
C
C**********************************************************************C
C     PERFORM FIRST CONTRACTION: G(AB| -) = E(CD| -)*RC(AB|CD).        C
C     THIS YIELDS ALL MQN SIGN POSSIBILITIES FOR C AND D.              C
C**********************************************************************C
C
C     LOOP OVER CARTESIAN INDEX ICMP FOR CENTRE AB (USE INDEX 6000)
      DO 6000 ICMP=1,3
C
C     TIME AT START OF FIRST CONTRACTION FOR THIS ICMP INDEX
      CALL CPU_TIME(T1I)
C
C     CARTESIAN INDEX ICMP AS A VECTOR, IDX
      CALL NCART(IDX,ICMP)
C
C     LOOP OVER ALL ADDRESSES FOR E(AB| -) FINITE EXPANSION
      DO IAB=1,NTUVAB
C
C       RESET CONTRACTION STORAGE ARRAYS G(AB| -)
        DO N=1,MAXN
          GABR11(N,IAB) = 0.0D0
          GABI11(N,IAB) = 0.0D0
          GABR21(N,IAB) = 0.0D0
          GABI21(N,IAB) = 0.0D0
        ENDDO
C
C       SKIP ENTIRE PROCESS IF E(AB| -) PASSES SCREENING CONDITION
        IF(IABR11(IAB,ICMP)+IABI11(IAB,ICMP)
     &                   +IABR21(IAB,ICMP)+IABI21(IAB,ICMP).EQ.0) THEN
          GOTO 401
        ENDIF
C
C >>>>> GAUNT INTERACTION
C
C       LOOP OVER ALL FINITE EXPANSION ADDRESSES FOR E(CD| -)
        DO ICD=1,NTUVCD
C
C         CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
          IRABCD = IABC(IA(IAB)+IA(ICD),IB(IAB)+IB(ICD),IC(IAB)+IC(ICD))
C
C         SKIP THIS STEP IF THE RC(AB|CD) PASSES SCREENING CONDITION
          IF(IRC(IRABCD).EQ.0) GOTO 402
C
C         SKIP THIS STEP IF THE E(CD) PASSES SCREENING CONDITION
          IF(ICDR11(ICD,ICMP)+ICDI11(ICD,ICMP)
     &                   +ICDR21(ICD,ICMP)+ICDI21(ICD,ICMP).EQ.0) THEN
            GOTO 402
          ENDIF
C
C         CONTRIBUTIONS TO Re{G(AB|--)} FROM EACH Re{E(CD|--)} ADDRESS
          IF(ISYM.EQ.1.AND.IABR11(IAB,ICMP).EQ.0) GOTO 411
          IF(ICDR11(ICD,ICMP).EQ.1) THEN
            DO N=1,MAXN
              GABR11(N,IAB) = GABR11(N,IAB)
     &                    + DREAL(ECD11(IMAP(N),ICD,ICMP))*RC(N,IRABCD)
            ENDDO
          ENDIF
411       CONTINUE
C
C         CONTRIBUTIONS TO Im{G(AB|--)} FROM EACH Im{E(CD|--)} ADDRESS
          IF(ISYM.EQ.1.AND.IABI11(IAB,ICMP).EQ.0) GOTO 412
          IF(ICDI11(ICD,ICMP).EQ.1) THEN
            DO N=1,MAXN
              GABI11(N,IAB) = GABI11(N,IAB)
     &                    + DIMAG(ECD11(IMAP(N),ICD,ICMP))*RC(N,IRABCD)
            ENDDO
          ENDIF
412       CONTINUE
C
C         CONTRIBUTIONS TO Re{G(AB|+-)} FROM EACH Re{E(CD|+-)} ADDRESS
          IF(ISYM.EQ.1.AND.IABR21(IAB,ICMP).EQ.0) GOTO 413
          IF(ICDR21(ICD,ICMP).EQ.1) THEN
            DO N=1,MAXN
              GABR21(N,IAB) = GABR21(N,IAB)
     &                    + DREAL(ECD21(IMAP(N),ICD,ICMP))*RC(N,IRABCD)
            ENDDO
          ENDIF
413       CONTINUE
C
C         CONTRIBUTIONS TO Im{G(AB|+-)} FROM EACH Im{E(CD|+-)} ADDRESS
          IF(ISYM.EQ.1.AND.IABI21(IAB,ICMP).EQ.0) GOTO 414
          IF(ICDI21(ICD,ICMP).EQ.1) THEN
            DO N=1,MAXN
              GABI21(N,IAB) = GABI21(N,IAB)
     &                    + DIMAG(ECD21(IMAP(N),ICD,ICMP))*RC(N,IRABCD)
            ENDDO
          ENDIF
414       CONTINUE
C
C         SKIP POINT FOR RC(AB|CD) AND E(CD) SCREENING
402       CONTINUE
C
C >>>>>   GAUGE TERM (REQUIRES ADDITIONAL CARTESIAN SUM Q')
          IF(BGAUNT) GOTO 450
C
C         LOOP OVER CARTESIAN INDEX JCMP FOR CENTRE CD
          DO JCMP=1,3
C
C           SKIP THIS STEP IF THE E(CD) PASSES SCREENING CONDITION
            IF(ICDR11(ICD,JCMP)+ICDI11(ICD,JCMP)
     &                     +ICDR21(ICD,JCMP)+ICDI21(ICD,JCMP).EQ.0) THEN
              GOTO 403
            ENDIF
C
C           CARTESIAN INDEX JCMP AS A VECTOR, JDX
            CALL NCART(JDX,JCMP)
C
C           NEW ADDRESS DEPENDING ON JCMP CARTESIAN INDEX
            IF(JCMP.EQ.1) THEN
              RTP = DFLOAT(IA(IAB)+IA(ICD))
            ELSEIF(JCMP.EQ.2) THEN
              RTP = DFLOAT(IB(IAB)+IB(ICD))
            ELSEIF(JCMP.EQ.3) THEN
              RTP = DFLOAT(IC(IAB)+IC(ICD))
            ENDIF
C
C           FIRST CONTRIBUTION ADDRESS
            I1 = IA(IAB)+IA(ICD)+IDX(1)+JDX(1)
            J1 = IB(IAB)+IB(ICD)+IDX(2)+JDX(2)
            K1 = IC(IAB)+IC(ICD)+IDX(3)+JDX(3)
C
C           CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
            IADR1 = IABC(I1,J1,K1)
C
C           SECOND CONTRIBUTION ADDRESS
            I2 = IA(IAB)+IA(ICD)+IDX(1)
            J2 = IB(IAB)+IB(ICD)+IDX(2)
            K2 = IC(IAB)+IC(ICD)+IDX(3)
C
C           CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
            IADR2 = IABC(I2,J2,K2)
C
C           THIRD CONTRIBUTION ADDRESS
            I3 = IA(IAB)+IA(ICD)+IDX(1)-JDX(1)
            J3 = IB(IAB)+IB(ICD)+IDX(2)-JDX(2)
            K3 = IC(IAB)+IC(ICD)+IDX(3)-JDX(3)
C
C           CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
            IF(I3.GE.0.AND.J3.GE.0.AND.K3.GE.0) THEN
              IADR3 = IABC(I3,J3,K3)
            ELSE
              IADR3 = 0
            ENDIF
C
C           SKIP THIS STEP IF THE RC(AB|CD) PASSES SCREENING CONDITION
            IF(IADR3.NE.0) THEN
              IF(IRC(IADR1)+IRC(IADR2)+IRC(IADR3).EQ.0) GOTO 403
            ELSE
              IF(IRC(IADR1)+IRC(IADR2).EQ.0) GOTO 403
            ENDIF
C
C           PRE-FACTORS FOR THE UPCOMING CONTRACTION
            DO N=1,MAXN
              T1 = RC(N,IADR1)*0.5D0/APH(N)
              T2 = RC(N,IADR2)*PQ(N,JCMP)
              IF(I3.GE.0.AND.J3.GE.0.AND.K3.GE.0) THEN
                T3 = RC(N,IADR3)*RTP
              ELSE
                T3 = 0.0D0
              ENDIF
              T(N) = T1-T2+T3
            ENDDO
C
C           CONTRIBUTIONS TO Re{G(AB|--)} FROM EACH Re{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABR11(IAB,ICMP).EQ.0) GOTO 415
            IF(ICDR11(ICD,JCMP).EQ.1) THEN
              DO N=1,MAXN
                GABR11(N,IAB) = GABR11(N,IAB)
     &                            - DREAL(ECD11(IMAP(N),ICD,JCMP))*T(N)
              ENDDO
            ENDIF
415         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|--)} FROM EACH Im{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABI11(IAB,ICMP).EQ.0) GOTO 416
            IF(ICDI11(ICD,JCMP).EQ.1) THEN
              DO N=1,MAXN
                GABI11(N,IAB) = GABI11(N,IAB)
     &                            - DIMAG(ECD11(IMAP(N),ICD,JCMP))*T(N)
              ENDDO
            ENDIF
416         CONTINUE
C
C           CONTRIBUTIONS TO Re{G(AB|--)} FROM EACH Re{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABR21(IAB,ICMP).EQ.0) GOTO 417
            IF(ICDR21(ICD,JCMP).EQ.1) THEN
              DO N=1,MAXN
                GABR21(N,IAB) = GABR21(N,IAB)
     &                            - DREAL(ECD21(IMAP(N),ICD,JCMP))*T(N)
              ENDDO
            ENDIF
417         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|--)} FROM EACH Im{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABI21(IAB,ICMP).EQ.0) GOTO 418
            IF(ICDI21(ICD,JCMP).EQ.1) THEN
              DO N=1,MAXN
                GABI21(N,IAB) = GABI21(N,IAB)
     &                            - DIMAG(ECD21(IMAP(N),ICD,JCMP))*T(N)
              ENDDO
            ENDIF
418         CONTINUE
C
C         SKIP POINT FOR E(CD) SCREENING
403       CONTINUE
C
C         END LOOP OVER CARTESIAN INDEX JCMP FOR CENTRE CD
          ENDDO
C
C         SKIP POINT FOR GAUNT INTERACTION ONLY
450       CONTINUE
C
C       END LOOP OVER E(CD|  ) FINITE EXPANSION ADDRESSES
        ENDDO
C
C       SKIP POINT FOR E(AB|  ) SCREENING
401     CONTINUE
C
C     END LOOP OVER E(AB|  ) FINITE EXPANSION ADDRESSES
      ENDDO
C
C     TIME AT END OF FIRST CONTRACTION FOR THIS ICMP INDEX
      CALL CPU_TIME(T1F)
      TBC1 = TBC1+T1F-T1I
C
C**********************************************************************C
C     PERFORM SECOND CONTRACTION: ( -| -) = E(AB| -)*G(AB| -).         C
C     THIS YIELDS A FULL BATCH OF TWO-ELECTRON INTEGRALS (16 PERM'NS). C
C**********************************************************************C
C
C     CALCULATE PHASES FOR BASIS FUNCTION OVERLAP COMBINATIONS
      PAB =-ISIGN(1,KQN(1)*KQN(2))*(-1)**((MQN(1)-MQN(2))/2)
      PCD =-ISIGN(1,KQN(3)*KQN(4))*(-1)**((MQN(3)-MQN(4))/2)
C
      PABCD = PAB*PCD
C
C     1ST SET: ( 1) = (--|--)   ( 4) = (--|++)
C              (16) = (++|++)   (13) = (++|--)
C
C     RESET CONTRACTION STORAGE LISTS
      DO N=1,MAXN
        QR1(N) = 0.0D0
        QI1(N) = 0.0D0
        QR2(N) = 0.0D0
        QI2(N) = 0.0D0
      ENDDO
C
C     RAW CONTRACTION (--|--) = E(AB|--)*(Re{G(AB|--)} + i*Im{G(AB|--)})
      DO IAB=1,NTUVAB
        IF(IABR11(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB11(IJ,IAB,ICMP))*GABR11(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB11(IJ,IAB,ICMP))*GABI11(N,IAB)
          ENDDO
        ENDIF
        IF(IABI11(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB11(IJ,IAB,ICMP))*GABR11(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB11(IJ,IAB,ICMP))*GABI11(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     ADD THIS ICMP TERM TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N,1 ) = RR(N,1 ) +     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N,4 ) = RR(N,4 ) + PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
      ENDDO
C
      IF(ISYM.EQ.1) GOTO 501
C
C     2ND SET: ( 3) = (--|+-)   ( 2) = (--|-+)
C              (14) = (++|-+)   (15) = (++|+-)
C
C     RESET CONTRACTION STORAGE LISTS
      DO N=1,MAXN
        QR1(N) = 0.0D0
        QI1(N) = 0.0D0
        QR2(N) = 0.0D0
        QI2(N) = 0.0D0
      ENDDO
C
C     RAW CONTRACTION (--|+-) = E(AB|--)*(Re{G(AB|+-)} + i*Im{G(AB|+-)})
      DO IAB=1,NTUVAB
        IF(IABR11(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB11(IJ,IAB,ICMP))*GABR21(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB11(IJ,IAB,ICMP))*GABI21(N,IAB)
          ENDDO
        ENDIF
        IF(IABI11(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB11(IJ,IAB,ICMP))*GABR21(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB11(IJ,IAB,ICMP))*GABI21(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     ADD THIS ICMP TERM TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N,3 ) = RR(N,3 ) +     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N,2 ) = RR(N,2 ) - PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
      ENDDO
C
C     3RD SET: ( 9) = (+-|--)   (12) = (+-|++)
C              ( 8) = (-+|++)   ( 5) = (-+|--)
C
C     RESET CONTRACTION STORAGE LISTS
      DO N=1,MAXN
        QR1(N) = 0.0D0
        QI1(N) = 0.0D0
        QR2(N) = 0.0D0
        QI2(N) = 0.0D0
      ENDDO
C
C     RAW CONTRACTION (+-|--) = E(AB|+-)*(Re{G(AB|--)} + i*Im{G(AB|--)})
      DO IAB=1,NTUVAB
        IF(IABR21(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB21(IJ,IAB,ICMP))*GABR11(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB21(IJ,IAB,ICMP))*GABI11(N,IAB)
          ENDDO
        ENDIF
        IF(IABI21(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB21(IJ,IAB,ICMP))*GABR11(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB21(IJ,IAB,ICMP))*GABI11(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     ADD THIS ICMP TERM TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N,9 ) = RR(N,9 ) +     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N,12) = RR(N,12) + PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
      ENDDO
C
501   CONTINUE
C
C     4TH SET: (11) = (+-|+-)   (10) = (+-|-+)
C              ( 6) = (-+|-+)   ( 7) = (-+|+-)
C
C     RESET CONTRACTION STORAGE LISTS
      DO N=1,MAXN
        QR1(N) = 0.0D0
        QI1(N) = 0.0D0
        QR2(N) = 0.0D0
        QI2(N) = 0.0D0
      ENDDO
C
C     RAW CONTRACTION (+-|+-) = E(AB|+-)*(Re{G(AB|+-)} + i*Im{G(AB|+-)})
      DO IAB=1,NTUVAB
        IF(IABR21(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB21(IJ,IAB,ICMP))*GABR21(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB21(IJ,IAB,ICMP))*GABI21(N,IAB)
          ENDDO
        ENDIF
        IF(IABI21(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB21(IJ,IAB,ICMP))*GABR21(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB21(IJ,IAB,ICMP))*GABI21(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     ADD THIS ICMP TERM TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N,11) = RR(N,11) +     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N,10) = RR(N,10) - PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
      ENDDO
C
C     TIME AT END OF SECOND CONTRACTION FOR THIS ICMP INDEX
      CALL CPU_TIME(T2F)
      TBC2 = TBC2+T2F-T1F
C
C     END LOOP OVER CARTESIAN INDICES {IX,IY,IZ}
6000  CONTINUE
C
C     HALF OF THE RR ARRAY CAN BE GENERATED WITH PHASE RELATIONS
      DO N=1,MAXN
        RR(N,16) = PABCD*DCONJG(RR(N,1 ))
        RR(N,13) = PABCD*DCONJG(RR(N,4 ))
        RR(N,14) =-PABCD*DCONJG(RR(N,3 ))
        RR(N,15) =-PABCD*DCONJG(RR(N,2 ))
        RR(N,8 ) =-PABCD*DCONJG(RR(N,9 ))
        RR(N,5 ) =-PABCD*DCONJG(RR(N,12))
        RR(N,6 ) = PABCD*DCONJG(RR(N,11))
        RR(N,7 ) = PABCD*DCONJG(RR(N,10))
      ENDDO     
C
C**********************************************************************C
C     BREIT INTEGRAL BATCH NOW FULLY CONSTRUCTED                       C
C**********************************************************************C
C
C     INCLUDE THE OUTSIDE FACTOR OF (-1/2) AND MOVE TO FULL ARRAY
      DO N=1,MAXN
        DO ITG=1,16
          RR(N,ITG) =-0.5D0*PRE(N)*RR(N,ITG)
        ENDDO
      ENDDO
C
      RETURN
      END

