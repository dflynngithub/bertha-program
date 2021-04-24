      SUBROUTINE ERI(RR,XYZ,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN)
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
C  ▶ ITN  - COMPONENT OVERLAP (T,T') FOR CD. ITN(I) = {LL,LS,SL,SS}.   C
C  OUTPUT:                                                             C
C  ▶ RR   - ERI'S FOR BLOCK AB, ALL 16 MQN SIGN COMBINATIONS.          C
C**********************************************************************C
      INCLUDE 'param.h'
C
      CHARACTER*4 HMLT
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4),LQN(4),ITN(2)
      DIMENSION PQ(MB2,3),APH(MB2),PRE(MB2),RC(MB2,MRC)
      DIMENSION IABR11(MEQ),IABR21(MEQ),IABI11(MEQ),IABI21(MEQ)
      DIMENSION ICDR11(MEQ),ICDR21(MEQ),ICDI11(MEQ),ICDI21(MEQ)
      DIMENSION IRC(MRC)
      DIMENSION RCTTFL(20*MFL),IRCTTFL(MFL)
      DIMENSION GABR11(MB2,MEQ),GABI11(MB2,MEQ),
     &          GABR21(MB2,MEQ),GABI21(MB2,MEQ)
      DIMENSION QR1(MB2),QI1(MB2),QR2(MB2),QI2(MB2)
C
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 EAB11(MB2,MEQ),ECD11(MB2,MEQ),
     &           EAB21(MB2,MEQ),ECD21(MB2,MEQ)
C
      SAVE EAB11,EAB21,IABR11,IABI11,IABR21,IABI21
      SAVE ECD11,ECD21,ICDR11,ICDI11,ICDR21,ICDI21
      SAVE RCTTFL,IRCTTFL
C
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS,IABSL,ICDSL
      COMMON/IRCM/IEAB,IECD,IRIJ(MBS,MBS)
      COMMON/ISCR/IMTX(MB2,11),ISCR(MB2),IMAP(MB2),IBCH,ITOG,MAXN
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PRMS/HMLT,ICLC,INEW,IATM,ILIN
      COMMON/TGGL/ILEV,ICB1,ISML,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/TMMD/TELL,TESS,TELS,TESL,TRLL,TRSS,TRLS,TRSL,TRBR
      COMMON/TSCF/TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRR,TCC1,TCC2,TCMC,
     &            TB1B,TB1R,TB1F,TB1M,TBEC,TBRM,TBRR,TBC1,TBC2,TBMC,
     &            TQMX,THMX,TC1T,TC2T,TB1T,TB2T,TEIG,TSCR,TTOT,
     &            TC1S,TC2S,TB1S,TB2S
C
C     EQ-COEFFICIENT SENSITIVITY PARAMETER
      DATA SENS/1.0D-10/
C
C     ILLEGAL COMPONENT OVERLAP CHECKER
      DO IT=1,2
        IF(ITN(IT).NE.1.AND.ITN(IT).NE.4) THEN
          WRITE(6, *) 'In ERI: illegal component overlaps in ITN.'
          WRITE(7, *) 'In ERI: illegal component overlaps in ITN.'
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
C     INTEGRAL SKIPPING ON MOLECULAR GROUP SYMMETRY CLASS BASIS
      IF(IATM.EQ.1.OR.ILIN.EQ.1) THEN
        ISYM = 1
      ELSE
        ISYM = 0
      ENDIF
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NBAS(1)*NBAS(2)
      MAXCD = NBAS(3)*NBAS(4)
C
C     PHASE FACTORS FOR AB AND CD PAIR OVERLAPS
      IPHSAB = 1
      IPHSCD =-1
C
C     VRS MAXIMUM LAMBDA LEVEL FOR EQ(AB)-COEFFICIENTS
      IF(ITN(1).EQ.1) THEN
        LAMAB = LQN(1)+LQN(2)
      ELSEIF(ITN(1).EQ.4) THEN
        LAMAB = LQN(1)+LQN(2)+2
      ENDIF
C
C     VRS MAXIMUM LAMBDA LEVEL FOR EQ(CD)-COEFFICIENTS
      IF(ITN(2).EQ.1) THEN
        LAMCD = LQN(3)+LQN(4)
      ELSEIF(ITN(2).EQ.4) THEN
        LAMCD = LQN(3)+LQN(4)+2
      ENDIF
C
C     VRS MAXIMUM LAMBDA LEVEL FOR CONTRACTED R-INTEGRAL BATCH
      LAMABCD = LAMAB+LAMCD
C
C     VRS MAXIMUM LAMBDA LEVEL FOR RCTTFL SAVED LIST (ONLY IF MCNT.GT.2)
      IF(IERC.EQ.0.OR.HMLT.EQ.'NORL') THEN
        LAMABCDFL = LAMABCD
      ELSEIF(IERC.EQ.1.AND.HMLT.NE.'NORL') THEN
        IF(ILEV.EQ.1) THEN
          LAMABCDFL = LQN(1)+LQN(2)+LQN(3)+LQN(4)
        ELSEIF(ILEV.EQ.2) THEN
          LAMABCDFL = LQN(1)+LQN(2)+LQN(3)+LQN(4)+2
        ELSEIF(ILEV.EQ.3) THEN
          LAMABCDFL = LQN(1)+LQN(2)+LQN(3)+LQN(4)+4
        ENDIF
      ENDIF
C
C     VRS TOTAL LENGTH OF EQ-COEFFICIENT LISTS AND R-INTEGRAL BATCH
      NTUVAB     = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
      NTUVCD     = (LAMCD+1)*(LAMCD+2)*(LAMCD+3)/6
      NTUVABCD   = (LAMABCD+1)*(LAMABCD+2)*(LAMABCD+3)/6
      NTUVABCDFL = (LAMABCDFL+1)*(LAMABCDFL+2)*(LAMABCDFL+3)/6
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
      IF(IEQS.EQ.1) THEN
        IF(ITN(1).EQ.1) THEN
          DO ITUV=1,NTUVAB
            IAD = IABLL + (ITUV-1)*MAXAB
            DO M=1,MAXAB
              EAB11(M,ITUV) = DCMPLX(E0LLFL(IAD+M,1),E0LLFL(IAD+M,2))
              EAB21(M,ITUV) = DCMPLX(E0LLFL(IAD+M,3),E0LLFL(IAD+M,4))
            ENDDO
          ENDDO
        ELSEIF(ITN(1).EQ.4) THEN
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
C     OPTION 2: CALCULATE FROM SCRATCH
      IF(IEQS.EQ.0) THEN
        IF(ITN(1).EQ.1) THEN
          CALL EQLLMK(EAB11,EAB21,EXL,XYZ,KQN,MQN,NBAS,IPHSAB,1,2,0)
        ELSEIF(ITN(1).EQ.4) THEN
          CALL EQSSMK(EAB11,EAB21,EXL,XYZ,KQN,MQN,NBAS,IPHSAB,1,2,0)
        ENDIF
      ENDIF
C
C     RECORD THE TIME TAKEN TO GENERATE/READ THE E(AB|  ) COEFFICIENTS
      CALL CPU_TIME(TDM2)
      IF(ITN(1).EQ.1) THEN
        TELL = TELL+TDM2-TDM1
      ELSEIF(ITN(1).EQ.4) THEN
        TESS = TESS+TDM2-TDM1
      ENDIF
C
C     SCREENING: TEST E(AB| -) COLUMNS OF CARTESIAN INDEX (T ,U ,V )
      DO IAB=1,NTUVAB
C
C       Re{E(AB|--)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXAB
          ERAB11 = DREAL(EAB11(M,IAB))
          SUM = SUM + DABS(ERAB11)
          IF(SUM.GT.SENS) THEN
            IABR11(IAB) = 1
            GOTO 101
          ENDIF
        ENDDO
        IABR11(IAB) = 0
101     CONTINUE
C
C       Im{E(AB|--)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXAB
          EIAB11 = DIMAG(EAB11(M,IAB))
          SUM = SUM + DABS(EIAB11)
          IF(SUM.GT.SENS) THEN
            IABI11(IAB) = 1
            GOTO 102
          ENDIF
        ENDDO
        IABI11(IAB) = 0
102     CONTINUE
C
C       Re{E(AB|+-)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXAB
          ERAB21 = DREAL(EAB21(M,IAB))
          SUM = SUM + DABS(ERAB21)
          IF(SUM.GT.SENS) THEN
            IABR21(IAB) = 1
            GOTO 103
          ENDIF
        ENDDO
        IABR21(IAB) = 0
103     CONTINUE
C
C       Im{E(AB|+-)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXAB
          EIAB21 = DIMAG(EAB21(M,IAB))
          SUM = SUM + DABS(EIAB21)
          IF(SUM.GT.SENS) THEN
            IABI21(IAB) = 1
            GOTO 104
          ENDIF
        ENDDO
        IABI21(IAB) = 0
104     CONTINUE
C
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
      IF(IEQS.EQ.1) THEN
        IF(ITN(2).EQ.1) THEN
          DO ITUV=1,NTUVCD
            IAD = ICDLL + (ITUV-1)*MAXCD
            DO M=1,MAXCD
              ECD11(M,ITUV) = DCMPLX(E0LLFL(IAD+M,5),E0LLFL(IAD+M,6))
              ECD21(M,ITUV) = DCMPLX(E0LLFL(IAD+M,7),E0LLFL(IAD+M,8))
            ENDDO
          ENDDO
        ELSEIF(ITN(2).EQ.4) THEN
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
C     OPTION 2: CALCULATE FROM SCRATCH
      IF(IEQS.EQ.0) THEN
        IF(ITN(2).EQ.1) THEN
          CALL EQLLMK(ECD11,ECD21,EXL,XYZ,KQN,MQN,NBAS,IPHSCD,3,4,0)
        ELSEIF(ITN(2).EQ.4) THEN
          CALL EQSSMK(ECD11,ECD21,EXL,XYZ,KQN,MQN,NBAS,IPHSCD,3,4,0)
        ENDIF
      ENDIF
C
C     RECORD THE TIME TAKEN TO GENERATE/READ THE E(CD|  ) COEFFICIENTS
      CALL CPU_TIME(TDM2)
      IF(ITN(2).EQ.1) THEN
        TELL = TELL+TDM2-TDM1
      ELSEIF(ITN(2).EQ.4) THEN
        TESS = TESS+TDM2-TDM1
      ENDIF
C
C     SCREENING: TEST E(CD| -) COLUMNS OF CARTESIAN INDEX (T',U',V')
      DO ICD=1,NTUVCD
C
C       Re{E(CD|--)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXCD
          ERCD11 = DREAL(ECD11(M,ICD))
          SUM = SUM + DABS(ERCD11)
          IF(SUM.GT.SENS) THEN
            ICDR11(ICD) = 1
            GOTO 201
          ENDIF
        ENDDO
        ICDR11(ICD) = 0
201     CONTINUE
C
C       Im{E(CD|--)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXCD
          EICD11 = DIMAG(ECD11(M,ICD))
          SUM = SUM + DABS(EICD11)
          IF(SUM.GT.SENS) THEN
            ICDI11(ICD) = 1
            GOTO 202
          ENDIF
        ENDDO
        ICDI11(ICD) = 0
202     CONTINUE
C
C       Re{E(CD|+-)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXCD
          ERCD21 = DREAL(ECD21(M,ICD))
          SUM = SUM + DABS(ERCD21)
          IF(SUM.GT.SENS) THEN
            ICDR21(ICD) = 1
            GOTO 203
          ENDIF
        ENDDO
        ICDR21(ICD) = 0
203     CONTINUE
C
C       Im{E(CD|+-)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXCD
          EICD21 = DIMAG(ECD21(M,ICD))
          SUM = SUM + DABS(EICD21)
          IF(SUM.GT.SENS) THEN
            ICDI21(ICD) = 1
            GOTO 204
          ENDIF
        ENDDO
        ICDI21(ICD) = 0
204     CONTINUE
C
      ENDDO
C
C     DO NOT CALCULATE AGAIN UNTIL PROMPTED EXTERNALLY
      IECD = 0
C
200   CONTINUE
      CALL CPU_TIME(T2)
      TCEC = TCEC + T2 - T1
C
C**********************************************************************C
C     GENERATE NEW BATCH OF RC(AB|CD) FROM SCRATCH OR READ-IN          C
C**********************************************************************C
C
C     START TIME
      CALL CPU_TIME(TDM1)
C
C     READ FROM LOCAL RC(AB|CD) FILE
      IF(IERC.EQ.1.AND.IRIJ(IBAS,JBAS).EQ.0) THEN
C
        CALL CPU_TIME(T1)
C
C       STARTING ADDRESS FOR SAVED R(AB|CD) INTEGRALS
        IADRTT = (IJ-1)*NBAS(3)*NBAS(4)*NTUVABCDFL
C
C       READ RC(AB|CD) INTEGRALS FROM THIS STARTING POINT
        DO N=1,MAXN
          DO IABCD=1,NTUVABCD
            IAD = IADRTT + (IMAP(N)-1)*NTUVABCDFL + IABCD
            RC(N,IABCD) = RCTTFL(IAD)
          ENDDO
        ENDDO
C
C       STARTING ADDRESS FOR SCREENING FLAGS
        IADSCR = (IJ-1)*NTUVABCDFL
C
C       READ SCREENING FLAGS FROM THIS STARTING POINT
        DO IABCD=1,NTUVABCD
          IAD = IADSCR + IABCD
          IRC(IABCD) = IRCTTFL(IAD)
        ENDDO
C
C       RECORD TIME SPENT READING R-INTEGRALS
        CALL CPU_TIME(T2)
        TCRR = TCRR+T2-T1
C
C       NORMALISATION FACTORS FOR THIS BATCH
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.1) THEN
              N   = N+1
              EKL = EXL(KBAS,3)+EXL(LBAS,4)
              EMX = DSQRT(EIJ+EKL)*EIJ*EKL
              PRE(N) = 2.0D0*PI52/EMX
            ENDIF
          ENDDO
        ENDDO
C
      ENDIF
C
C     CALCULATE FROM SCRATCH
      IF(IERC.EQ.0.OR.IRIJ(IBAS,JBAS).EQ.1) THEN
C
        CALL CPU_TIME(T1)
C
C       GAUSSIAN OVERLAP CENTRE
        PX = (XYZ(1,1)*EXL(IBAS,1)+XYZ(1,2)*EXL(JBAS,2))/EIJ
        PY = (XYZ(2,1)*EXL(IBAS,1)+XYZ(2,2)*EXL(JBAS,2))/EIJ
        PZ = (XYZ(3,1)*EXL(IBAS,1)+XYZ(3,2)*EXL(JBAS,2))/EIJ
C
C       AUXILLIARY DATA FOR RMAKE ROUTINE
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF((IERC.EQ.1.AND.IRIJ(IBAS,JBAS).EQ.1)
     &                           .OR.(IERC.EQ.0.AND.ISCR(M).EQ.1)) THEN
              N   = N+1
              EKL = EXL(KBAS,3)+EXL(LBAS,4)
              EMX = DSQRT(EIJ+EKL)*EIJ*EKL
              APH(N) = EIJ*EKL/(EIJ+EKL)
              PRE(N) = 2.0D0*PI52/EMX
              QX = (XYZ(1,3)*EXL(KBAS,3) + XYZ(1,4)*EXL(LBAS,4))/EKL
              QY = (XYZ(2,3)*EXL(KBAS,3) + XYZ(2,4)*EXL(LBAS,4))/EKL
              QZ = (XYZ(3,3)*EXL(KBAS,3) + XYZ(3,4)*EXL(LBAS,4))/EKL
              PQ(N,1) = QX-PX
              PQ(N,2) = QY-PY
              PQ(N,3) = QZ-PZ
            ENDIF
          ENDDO
        ENDDO
C
C       EXTEND MAXN IF GENERATING A FULL SET OF INTEGRALS
        MAXN = MAXCD
C
C       GENERATE R-INTEGRALS
        CALL RMAKE(RC,PQ,APH,MAXN,LAMABCDFL)
C
C       SCREENING: TEST RC(AB|CD) COLUMNS WITH INDEX (T+T',U+U',V+V')
        DO IABCDFL=1,NTUVABCDFL
C
C         SUM OF RC(AB|CD) MAGNITUDES
          SUM = 0.0D0
          DO N=1,MAXN
            SUM = SUM + DABS(RC(N,IABCDFL))
            IF(SUM.GT.SENS) THEN
              IRC(IABCDFL) = 1
              GOTO 301
            ENDIF
          ENDDO
          IRC(IABCDFL) = 0
301       CONTINUE
C
        ENDDO
C
C       SAVE THIS SET TO APPROPRIATE CLASS ADDRESS
        IF(IERC.EQ.1) THEN
C
C         TEST WHETHER FINAL ADDRESS IS STILL INSIDE ARRAY BOUNDS
          ILIM = IJ*NBAS(3)*NBAS(4)*NTUVABCDFL
C
          IF(ILIM.GT.20*MFL) THEN
C           OUT OF BOUNDS: PRINT WARNING BUT KEEP GOING
            WRITE(6, *) 'In ERI: RCTT words exceed allocated limit.'
            WRITE(7, *) 'In ERI: RCTT words exceed allocated limit.'
            GOTO 300
          ELSE
C           DO NOT CALCULATE AGAIN UNTIL PROMPTED EXTERNALLY
            IRIJ(IBAS,JBAS) = 0
          ENDIF
C
C         STARTING ADDRESS FOR SAVED R(AB|CD) INTEGRALS
          IADRTT = (IJ-1)*NBAS(3)*NBAS(4)*NTUVABCDFL

C         COPY THIS BATCH OF INTEGRALS TO A SAVED LIST
          DO N=1,MAXN
            DO IABCDFL=1,NTUVABCDFL
              IAD = IADRTT + (N-1)*NTUVABCDFL + IABCDFL
              RCTTFL(IAD)  = RC(N,IABCDFL)
            ENDDO
          ENDDO
C
C         STARTING ADDRESS FOR SCREENING FLAGS
          IADSCR = (IJ-1)*NTUVABCDFL
C
C         COPY SCREENING MARKERS TO A SAVED LIST
          DO IABCDFL=1,NTUVABCDFL
            IAD = IADSCR + IABCDFL
            IRCTTFL(IAD) = IRC(IABCDFL)
          ENDDO
C
C         SHORTEN THE CURRENT RC LIST WITH IMAP FROM SCREENING
          M = 0
          N = 0
          DO KBAS=1,NBAS(3)
            DO LBAS=1,NBAS(4)
              M = M+1
              IF(ISCR(M).EQ.1) THEN
                N      = N+1
                PRE(N) = PRE(M)
                DO IABCD=1,NTUVABCD
                  RC(N,IABCD) = RC(M,IABCD)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
C
        ENDIF
300     CONTINUE
C
        CALL CPU_TIME(T2)
        TCRM = TCRM+T2-T1
C
      ENDIF
C
C     RECORD THE TIME TAKEN TO GENERATE THE RC(AB|CD) BATCH
      CALL CPU_TIME(TDM2)
      IF(ITN(1).EQ.1.AND.ITN(2).EQ.1) THEN
        TRLL = TRLL+TDM2-TDM1
      ELSEIF(ITN(1).EQ.4.AND.ITN(2).EQ.4) THEN
        TRSS = TRSS+TDM2-TDM1
      ELSE
        TRLS = TRLS+TDM2-TDM1
      ENDIF
C
C**********************************************************************C
C     PERFORM FIRST CONTRACTION: G(AB| -) = E(CD| -)*RC(AB|CD).        C
C     THIS YIELDS ALL MQN SIGN POSSIBILITIES FOR C AND D.              C
C**********************************************************************C
C
C     TIME AT START OF FIRST CONTRACTION
      CALL CPU_TIME(T1)
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
        IF(IABR11(IAB)+IABI11(IAB)+IABR21(IAB)+IABI21(IAB).EQ.0) THEN
          GOTO 401
        ENDIF
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
          IF(ICDR11(ICD)+ICDI11(ICD)+ICDR21(ICD)+ICDI21(ICD).EQ.0) THEN
            GOTO 402
          ENDIF
C
C         CONTRIBUTIONS TO Re{G(AB|--)} FROM EACH Re{E(CD|--)} ADDRESS
          IF(ISYM.EQ.1.AND.IABR11(IAB).EQ.0) GOTO 411
          IF(ICDR11(ICD).EQ.1) THEN
            DO N=1,MAXN
              GABR11(N,IAB) = GABR11(N,IAB)
     &                       + DREAL(ECD11(IMAP(N),ICD))*RC(N,IRABCD)
            ENDDO
          ENDIF
411       CONTINUE
C
C         CONTRIBUTIONS TO Im{G(AB|--)} FROM EACH Im{E(CD|--)} ADDRESS
          IF(ISYM.EQ.1.AND.IABI11(IAB).EQ.0) GOTO 412
          IF(ICDI11(ICD).EQ.1) THEN
            DO N=1,MAXN
              GABI11(N,IAB) = GABI11(N,IAB)
     &                         + DIMAG(ECD11(IMAP(N),ICD))*RC(N,IRABCD)
            ENDDO
          ENDIF
412       CONTINUE
C
C         CONTRIBUTIONS TO Re{G(AB|+-)} FROM EACH Re{E(CD|+-)} ADDRESS
          IF(ISYM.EQ.1.AND.IABR21(IAB).EQ.0) GOTO 413
          IF(ICDR21(ICD).EQ.1) THEN
            DO N=1,MAXN
              GABR21(N,IAB) = GABR21(N,IAB)
     &                         + DREAL(ECD21(IMAP(N),ICD))*RC(N,IRABCD)
            ENDDO
          ENDIF
413       CONTINUE
C
C         CONTRIBUTIONS TO Im{G(AB|+-)} FROM EACH Im{E(CD|+-)} ADDRESS
          IF(ISYM.EQ.1.AND.IABI21(IAB).EQ.0) GOTO 414
          IF(ICDI21(ICD).EQ.1) THEN
            DO N=1,MAXN
              GABI21(N,IAB) = GABI21(N,IAB)
     &                         + DIMAG(ECD21(IMAP(N),ICD))*RC(N,IRABCD)
            ENDDO
          ENDIF
414       CONTINUE
C
C         SKIP POINT FOR RC(AB|CD) SCREENING
402       CONTINUE
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
C     TIME AT END OF FIRST CONTRACTION
      CALL CPU_TIME(T2)
      TCC1 = TCC1+T2-T1
C
C**********************************************************************C
C     PERFORM SECOND CONTRACTION: ( -| -) = E(AB| -)*G(AB| -).         C
C     THIS YIELDS A FULL BATCH OF TWO-ELECTRON INTEGRALS (16 PERM'NS). C
C**********************************************************************C
C
C     CALCULATE PHASES FOR BASIS FUNCTION OVERLAP COMBINATIONS
      PAB = ISIGN(1,KQN(1)*KQN(2))*(-1)**((MQN(1)-MQN(2))/2)
      PCD = ISIGN(1,KQN(3)*KQN(4))*(-1)**((MQN(3)-MQN(4))/2)
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
        IF(IABR11(IAB).EQ.1) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB11(IJ,IAB))*GABR11(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB11(IJ,IAB))*GABI11(N,IAB)
          ENDDO
        ENDIF
        IF(IABI11(IAB).EQ.1) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB11(IJ,IAB))*GABR11(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB11(IJ,IAB))*GABI11(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     APPLY PHASE RELATIONS AND NORMALISATION FACTORS TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N, 1) =     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N, 4) = PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
        RR(N,16) = PABCD*DCONJG(RR(N, 1))
        RR(N,13) = PABCD*DCONJG(RR(N, 4))
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
        IF(IABR11(IAB).EQ.1) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB11(IJ,IAB))*GABR21(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB11(IJ,IAB))*GABI21(N,IAB)
          ENDDO
        ENDIF
        IF(IABI11(IAB).EQ.1) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB11(IJ,IAB))*GABR21(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB11(IJ,IAB))*GABI21(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     APPLY PHASE RELATIONS AND NORMALISATION FACTORS TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N, 3) =     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N, 2) =-PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
        RR(N,14) =-PABCD*DCONJG(RR(N, 3))
        RR(N,15) =-PABCD*DCONJG(RR(N, 2))
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
        IF(IABR21(IAB).EQ.1) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB21(IJ,IAB))*GABR11(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB21(IJ,IAB))*GABI11(N,IAB)
          ENDDO
        ENDIF
        IF(IABI21(IAB).EQ.1) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB21(IJ,IAB))*GABR11(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB21(IJ,IAB))*GABI11(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     APPLY PHASE RELATIONS AND NORMALISATION FACTORS TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N, 9) =     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N,12) = PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
        RR(N, 8) =-PABCD*DCONJG(RR(N, 9))
        RR(N, 5) =-PABCD*DCONJG(RR(N,12))
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
        IF(IABR21(IAB).EQ.1) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB21(IJ,IAB))*GABR21(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB21(IJ,IAB))*GABI21(N,IAB)
          ENDDO
        ENDIF
        IF(IABI21(IAB).EQ.1) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB21(IJ,IAB))*GABR21(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB21(IJ,IAB))*GABI21(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     APPLY PHASE RELATIONS AND NORMALISATION FACTORS TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N,11) =     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N,10) =-PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
        RR(N, 6) = PABCD*DCONJG(RR(N,11))
        RR(N, 7) = PABCD*DCONJG(RR(N,10))
      ENDDO
C
C     TIME AT END OF SECOND CONTRACTION
      CALL CPU_TIME(T3)
      TCC2 = TCC2+T3-T2
C
C**********************************************************************C
C     COULOMB INTEGRAL BATCH NOW FULLY CONSTRUCTED                     C
C**********************************************************************C
C
C     INCLUDE THE R-INTEGRAL NORMALISATION FACTOR
      DO N=1,MAXN
        DO ITG=1,16
          RR(N,ITG) = PRE(N)*RR(N,ITG)
        ENDDO
      ENDDO
C
      RETURN
      END

