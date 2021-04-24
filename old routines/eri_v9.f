      SUBROUTINE ERI(RR,XYZ,ICNT,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN)
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
C  ▶ XYZ  - FULL SET OF CARTESIAN BASIS CENTRES.                       C
C  ▶ ICNT - FULL SET OF BASIS CENTRE ORIGINS (NUCLEAR LABELS).         C
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
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL L0CASE
C
      DIMENSION EXL(MBS,4),XYZ(3,4)
      DIMENSION ICNT(4),KQN(4),MQN(4),NBAS(4),LQN(4),ITN(2)
      DIMENSION PQ(MB2,3),APH(MB2),PRE(MB2),RC(MB2,MRC)
      DIMENSION IABR11(MB2,MEQ),IABI11(MB2,MEQ),
     &          IABR21(MB2,MEQ),IABI21(MB2,MEQ)
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
      COMMON/E0LL/E0LLFL(MFL,4),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,4),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS,IABSL,ICDSL
      COMMON/IRCM/IEAB,IECD,NCD,IGAB,IRIJ(MBS,MBS)
      COMMON/ISCR/IMTX(MB2,11),ISCR(MB2),IMAP(MB2),IBCH,ITOG,MAXN
      COMMON/LSHF/SHLEV(4),SHLV,ILEV
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/TMMD/TELL,TESS,TELS,TESL,TRLL,TRSS,TRLS,TRSL,TRBR
      COMMON/TSCF/TC1A,TC1I,TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRW,TCC1,
     &            TCC2,TCMC,TB1A,TB1I,TB1B,TB1R,TB1F,TB1M,TB1T,TBEC,
     &            TBRM,TBRW,TBC1,TBC2,TBMC,TSMX,TUMX,THMX,TAMX,TC1T,
     &            TC2T,TCVT,TB2T,TACC,TEIG,TSCR,TTOT,TC2S,TB2S
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
C     EVALUATE LQNS FOR BASIS FUNCTIONS (A,B,C,D)
      DO N=1,4
        LQN(N) = LVAL(KQN(N))
      ENDDO
C
C     SPECIAL CASE FOR S-TYPE OVERLAPS (ONLY EVER NEEDED ONCE)
      IF(LQN(1)+LQN(2)+LQN(3)+LQN(4).EQ.0) THEN
        L0CASE = .TRUE.
      ELSE
        L0CASE = .FALSE.
      ENDIF
C
C     INTEGRAL SKIPPING ON MOLECULAR GROUP SYMMETRY CLASS BASIS
      IF(SHAPE.EQ.'ATOMIC') THEN
        ISYM = 2
      ELSEIF(SHAPE.EQ.'DIATOM'.OR.SHAPE.EQ.'LINEAR') THEN
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
C     MCMURCHIE-DAVIDSON MAXINUM ORDER FOR EQ(AB)-COEFFICIENTS
      IF(ITN(1).EQ.1) THEN
        LAMAB = LQN(1)+LQN(2)
      ELSEIF(ITN(1).EQ.4) THEN
        LAMAB = LQN(1)+LQN(2)+2
      ENDIF
C
C     MCMURCHIE-DAVIDSON MAXIMUM ORDER FOR EQ(CD)-COEFFICIENTS
      IF(ITN(2).EQ.1) THEN
        LAMCD = LQN(3)+LQN(4)
      ELSEIF(ITN(2).EQ.4) THEN
        LAMCD = LQN(3)+LQN(4)+2
      ENDIF
C
C     MCMURCHIE-DAVIDSON MAXIMUM ORDER FOR CONTRACTED R-INTEGRAL BATCH
      LAMABCD = LAMAB+LAMCD
C
C     MAXIMUM LAMBDA ORDER FOR RCTTFL SAVED LIST (ONLY IF MCNT.GT.2)
      IF(.NOT.RCFILE.OR.HMLT.EQ.'NORL') THEN
        LAMABCDFL = LAMABCD
      ELSEIF(RCFILE.AND.HMLT.NE.'NORL') THEN
        IF(ILEV.EQ.1) THEN
          LAMABCDFL = LQN(1)+LQN(2)+LQN(3)+LQN(4)
        ELSEIF(ILEV.EQ.2) THEN
          LAMABCDFL = LQN(1)+LQN(2)+LQN(3)+LQN(4)+2
        ELSEIF(ILEV.EQ.3) THEN
          LAMABCDFL = LQN(1)+LQN(2)+LQN(3)+LQN(4)+4
        ENDIF
      ENDIF
C
C     UTILISE EXPANSION SUBSET (WHEN RCFILE.EQ.FALSE)
      IF((ITN(1)+ITN(2))/3.EQ.ILEV) THEN
        ISKIP = 0
      ELSE
        ISKIP = 1
      ENDIF
C
C     PAIR AND GROUP EQ-COEFFICIENT LIST AND R-INTEGRAL BATCH LENGTHS
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
      CALL SYSTEM_CLOCK(ICL1,RATE)
      IF(IEAB.EQ.0) GOTO 100
C
C     START TIME
      CALL SYSTEM_CLOCK(ICL3,RATE)
C
      IF(EQFILE) THEN
C
C       OPTION 1: READ FROM LOCAL EQ-COEFFICIENT FILE
        IF(ITN(1).EQ.1) THEN
          DO IAB=1,NTUVAB
            IAD = IABLL + (IAB-1)*MAXAB
            DO M=1,MAXAB
              EAB11(M,IAB) = DCMPLX(E0LLFL(IAD+M,1),E0LLFL(IAD+M,2))
              EAB21(M,IAB) = DCMPLX(E0LLFL(IAD+M,3),E0LLFL(IAD+M,4))
            ENDDO
          ENDDO
        ELSEIF(ITN(1).EQ.4) THEN
          DO IAB=1,NTUVAB
            IAD = IABSS + (IAB-1)*MAXAB
            DO M=1,MAXAB
              EAB11(M,IAB) = DCMPLX(E0SSFL(IAD+M,1),E0SSFL(IAD+M,2))
              EAB21(M,IAB) = DCMPLX(E0SSFL(IAD+M,3),E0SSFL(IAD+M,4))
            ENDDO
          ENDDO
        ENDIF
C
      ELSE
C
C       OPTION 2: CALCULATE FROM SCRATCH
        IF(ITN(1).EQ.1) THEN
          CALL EQLLMK(EAB11,EAB21,EXL,XYZ,KQN,MQN,NBAS,IPHSAB,1,2,0)
        ELSEIF(ITN(1).EQ.4) THEN
          CALL EQSSMK(EAB11,EAB21,EXL,XYZ,KQN,MQN,NBAS,IPHSAB,1,2,0)
        ENDIF
C
      ENDIF
C
C     RECORD THE TIME TAKEN TO GENERATE/READ THE E(AB|  ) COEFFICIENTS
      CALL SYSTEM_CLOCK(ICL4)
      IF(ITN(1).EQ.1) THEN
        TELL = TELL + DFLOAT(ICL4-ICL3)/RATE
      ELSEIF(ITN(1).EQ.4) THEN
        TESS = TESS + DFLOAT(ICL4-ICL3)/RATE
      ENDIF
C
C     SCREENING: TEST E(AB| -) COLUMNS OF CARTESIAN INDEX (T ,U ,V )
      DO IAB=1,NTUVAB
C
C       Re{E(AB|--)} COEFFICIENTS
        DO M=1,MAXAB
          ERAB11 = DREAL(EAB11(M,IAB))
          IF(DABS(ERAB11).GT.SENS) THEN
            IABR11(M,IAB) = 1
          ELSE
            IABR11(M,IAB) = 0
          ENDIF
        ENDDO
C
C       Im{E(AB|--)} COEFFICIENTS
        DO M=1,MAXAB
          EIAB11 = DIMAG(EAB11(M,IAB))
          IF(DABS(EIAB11).GT.SENS) THEN
            IABI11(M,IAB) = 1
          ELSE
            IABI11(M,IAB) = 0
          ENDIF
        ENDDO
C
C       Re{E(AB|+-)} COEFFICIENTS
        DO M=1,MAXAB
          ERAB21 = DREAL(EAB21(M,IAB))
          IF(DABS(ERAB21).GT.SENS) THEN
            IABR21(M,IAB) = 1
          ELSE
            IABR21(M,IAB) = 0
          ENDIF
        ENDDO
C
C       Im{E(AB|+-)} COEFFICIENTS
        DO M=1,MAXAB
          EIAB21 = DIMAG(EAB21(M,IAB))
          IF(DABS(EIAB21).GT.SENS) THEN
            IABI21(M,IAB) = 1
          ELSE
            IABI21(M,IAB) = 0
          ENDIF
        ENDDO
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
      CALL SYSTEM_CLOCK(ICL3,RATE)
C
C     OPTION 1: READ FROM LOCAL EQ-COEFFICIENT FILE
      IF(EQFILE) THEN
        IF(ITN(2).EQ.1) THEN
          DO ICD=1,NTUVCD
            IAD = ICDLL + (ICD-1)*MAXCD
            Z   = DFLOAT((-1)**(ILAM(ICD)))
            DO M=1,MAXCD
              ECD11(M,ICD) = Z*DCMPLX(E0LLFL(IAD+M,1),E0LLFL(IAD+M,2))
              ECD21(M,ICD) = Z*DCMPLX(E0LLFL(IAD+M,3),E0LLFL(IAD+M,4))
            ENDDO
          ENDDO
        ELSEIF(ITN(2).EQ.4) THEN
          DO ICD=1,NTUVCD
            IAD = ICDSS + (ICD-1)*MAXCD
            Z   = DFLOAT((-1)**(ILAM(ICD)))
            DO M=1,MAXCD
              ECD11(M,ICD) = Z*DCMPLX(E0SSFL(IAD+M,1),E0SSFL(IAD+M,2))
              ECD21(M,ICD) = Z*DCMPLX(E0SSFL(IAD+M,3),E0SSFL(IAD+M,4))
            ENDDO
          ENDDO
        ENDIF
C
C     OPTION 2: CALCULATE FROM SCRATCH
      ELSE
        IF(ITN(2).EQ.1) THEN
          CALL EQLLMK(ECD11,ECD21,EXL,XY Z,KQN,MQN,NBAS,IPHSCD,3,4,0)
        ELSEIF(ITN(2).EQ.4) THEN
          CALL EQSSMK(ECD11,ECD21,EXL,XYZ,KQN,MQN,NBAS,IPHSCD,3,4,0)
        ENDIF
      ENDIF
C
C     RECORD THE TIME TAKEN TO GENERATE/READ THE E(CD|  ) COEFFICIENTS
      CALL SYSTEM_CLOCK(ICL4)
      IF(ITN(2).EQ.1) THEN
        TELL = TELL + DFLOAT(ICL4-ICL3)/RATE
      ELSEIF(ITN(2).EQ.4) THEN
        TESS = TESS + DFLOAT(ICL4-ICL3)/RATE
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
      CALL SYSTEM_CLOCK(ICL2)
      TCEC = TCEC + DFLOAT(ICL2-ICL1)/RATE
C
C**********************************************************************C
C     GENERATE NEW BATCH OF RC(AB|CD) INTEGRALS IF PROMPTED            C
C**********************************************************************C
C
C     START TIME
      CALL SYSTEM_CLOCK(ICL1,RATE)
C
C     SKIP IF A HIGHER LEVEL OF THIS BATCH EXISTS IN RC
      IF(.NOT.RCFILE.AND.ISKIP.EQ.0) GOTO 300
C
C     SKIP IF A HIGHER LEVEL OF THIS BATCH EXISTS IN RC
      IF(L0CASE.AND.ISKIP.EQ.0) GOTO 300
C
C     SKIP IF INTEGRAL BATCH EXISTS IN FILE
      IF(RCFILE.AND.IRIJ(IBAS,JBAS).EQ.0) GOTO 300
C
      CALL SYSTEM_CLOCK(ICL3,RATE)
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
          IF(RCFILE.AND.L0CASE.AND.ISCR(M).EQ.0) GOTO 301
          IF(.NOT.RCFILE.AND.ISCR(M).EQ.0) GOTO 301
          N   = N+1
          EKL = EXL(KBAS,3)+EXL(LBAS,4)
          QX  = (XYZ(1,3)*EXL(KBAS,3)+XYZ(1,4)*EXL(LBAS,4))/EKL
          QY  = (XYZ(2,3)*EXL(KBAS,3)+XYZ(2,4)*EXL(LBAS,4))/EKL
          QZ  = (XYZ(3,3)*EXL(KBAS,3)+XYZ(3,4)*EXL(LBAS,4))/EKL
          PQ(N,1) = QX-PX
          PQ(N,2) = QY-PY
          PQ(N,3) = QZ-PZ
          APH(N)  = EIJ*EKL/(EIJ+EKL)
301       CONTINUE
        ENDDO
      ENDDO
C
C     BATCH SIZE AND EXPANSION LENGTH DEPENDS ON MODE
      IF(RCFILE.AND..NOT.L0CASE) THEN
        MBCH = MAXCD
        MLAM = LAMABCDFL
        MNTV = NTUVABCDFL
      ELSE
        MBCH = MAXN
        MLAM = LAMABCD
        MNTV = NTUVABCD
      ENDIF
C
C     GENERATE R-INTEGRALS
      CALL RMAKE(RC,PQ,APH,MBCH,MLAM)
C
C     SCREENING: TEST RC(AB|CD) COLUMNS WITH INDEX (T+T',U+U',V+V')
      DO INTV=1,MNTV
C
C       SUM OF RC(AB|CD) MAGNITUDES
        SUM = 0.0D0
        DO N=1,MBCH
          SUM = SUM + DABS(RC(N,INTV))
          IF(SUM.GT.SENS) THEN
            IRC(INTV) = 1
            GOTO 302
          ENDIF
        ENDDO
        IRC(INTV) = 0
302     CONTINUE
C
      ENDDO
      CALL SYSTEM_CLOCK(ICL4)
      TCRM = TCRM + DFLOAT(ICL4-ICL3)/RATE
C
C     CONTINUE ONLY IF INTEGRALS ARE TO BE SAVED TO LARGE FILE
      IF(L0CASE) GOTO 300
      IF(.NOT.RCFILE) GOTO 300
C
      CALL SYSTEM_CLOCK(ICL3,RATE)
C
C     TEST WHETHER FINAL ADDRESS IS STILL INSIDE ARRAY BOUNDS
      IF(20*MFL.LT.IJ*MAXCD*NTUVABCDFL) THEN
C       OUT OF BOUNDS: PRINT WARNING BUT KEEP GOING
        WRITE(6, *) 'In ERI: RCTT words exceed allocated limit.'
        WRITE(7, *) 'In ERI: RCTT words exceed allocated limit.'
        GOTO 300
      ELSE
C       DO NOT CALCULATE AGAIN UNTIL PROMPTED EXTERNALLY
        IRIJ(IBAS,JBAS) = 0
      ENDIF
C
C     STARTING ADDRESS FOR THIS BATCH OF SAVED R(AB|CD) INTEGRALS
      IADRTT = (IJ-1)*MAXCD*NTUVABCDFL
C
C     COPY THIS BATCH OF INTEGRALS TO A SAVED LIST
      DO IABCDFL=1,NTUVABCDFL
        IAD = IADRTT + MAXCD*(IABCDFL-1)
        DO M=1,MAXCD
          RCTTFL(IAD+M) = RC(M,IABCDFL)
        ENDDO
      ENDDO
C
C     STARTING ADDRESS FOR THIS BATCH OF SCREENING FLAGS
      IADSCR = (IJ-1)*NTUVABCDFL
C
C     COPY SCREENING MARKERS TO A SAVED LIST
      DO IABCDFL=1,NTUVABCDFL
        IRCTTFL(IADSCR+IABCDFL) = IRC(IABCDFL)
      ENDDO
C
      CALL SYSTEM_CLOCK(ICL4)
      TCRW = TCRW + DFLOAT(ICL4-ICL3)/RATE
C
300   CONTINUE
C
C     RECORD THE TIME TAKEN TO GENERATE THE RC(AB|CD) BATCH
      CALL SYSTEM_CLOCK(ICL2)
      IF(ITN(1).EQ.1.AND.ITN(2).EQ.1) THEN
        TRLL = TRLL + DFLOAT(ICL2-ICL1)/RATE
      ELSEIF(ITN(1).EQ.4.AND.ITN(2).EQ.4) THEN
        TRSS = TRSS + DFLOAT(ICL2-ICL1)/RATE
      ELSE
        TRLS = TRLS + DFLOAT(ICL2-ICL1)/RATE
      ENDIF
C
C**********************************************************************C
C     PERFORM FIRST CONTRACTION: G(AB| -) = E(CD| -)*RC(AB|CD).        C
C     THIS YIELDS ALL MQN SIGN POSSIBILITIES FOR C AND D.              C
C**********************************************************************C
C
C     TIME AT START OF FIRST CONTRACTION
      CALL SYSTEM_CLOCK(ICL1,RATE)
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
C       SKIP ENTIRE PROCESS IF E(AB| -) FAILS SCREENING CONDITION
        IABALL = IABR11(IJ,IAB) + IABI11(IJ,IAB)
     &         + IABR21(IJ,IAB) + IABI21(IJ,IAB)
        IF(IABALL.EQ.0) GOTO 401
C
C       CALCULATE DIRECTLY FROM SMALL RC(AB|CD) LOCAL ARRAY
        IF(.NOT.RCFILE.OR.L0CASE) THEN
C
C         LOOP OVER ALL FINITE EXPANSION ADDRESSES FOR E(CD| -)
          DO ICD=1,NTUVCD
C
C           SKIP THIS STEP IF THE E(CD) FAILS SCREENING CONDITION
            ICDALL = ICDR11(ICD)+ICDI11(ICD)+ICDR21(ICD)+ICDI21(ICD)
            IF(ICDALL.EQ.0) GOTO 402
C
C           CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
            IRABCD = IABC(IA(IAB)+IA(ICD),IB(IAB)+IB(ICD),
     &                                    IC(IAB)+IC(ICD))
C
C           SKIP THIS STEP IF THE RC(AB|CD) FAILS SCREENING CONDITION
            IF(IRC(IRABCD).EQ.0) GOTO 402
C
C           CONTRIBUTIONS TO Re{G(AB|--)} FROM EACH Re{E(CD|--)} ADDRESS
            IF(IABR11(IJ,IAB).EQ.0.AND.ISYM.GT.1) GOTO 411
            IF(ICDR11(ICD).EQ.0) GOTO 411
            DO N=1,MAXN
              GABR11(N,IAB) = GABR11(N,IAB)
     &                         + DREAL(ECD11(IMAP(N),ICD))*RC(N,IRABCD)
            ENDDO
411         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|--)} FROM EACH Im{E(CD|--)} ADDRESS
            IF(IABI11(IJ,IAB).EQ.0.AND.ISYM.GT.1) GOTO 412
            IF(ICDI11(ICD).EQ.0) GOTO 412
            DO N=1,MAXN
              GABI11(N,IAB) = GABI11(N,IAB)
     &                         + DIMAG(ECD11(IMAP(N),ICD))*RC(N,IRABCD)
            ENDDO
412         CONTINUE
C
C           CONTRIBUTIONS TO Re{G(AB|+-)} FROM EACH Re{E(CD|+-)} ADDRESS
            IF(IABR21(IJ,IAB).EQ.0.AND.ISYM.GT.1) GOTO 413
            IF(ICDR21(ICD).EQ.0) GOTO 413
            DO N=1,MAXN
              GABR21(N,IAB) = GABR21(N,IAB)
     &                         + DREAL(ECD21(IMAP(N),ICD))*RC(N,IRABCD)
            ENDDO
413         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|+-)} FROM EACH Im{E(CD|+-)} ADDRESS
            IF(IABI21(IJ,IAB).EQ.0.AND.ISYM.GT.1) GOTO 414
            IF(ICDI21(ICD).EQ.0) GOTO 414
            DO N=1,MAXN
              GABI21(N,IAB) = GABI21(N,IAB)
     &                         + DIMAG(ECD21(IMAP(N),ICD))*RC(N,IRABCD)
            ENDDO
414         CONTINUE
C
C           SKIP POINT FOR RC(AB|CD) SCREENING
402         CONTINUE
C
C         END LOOP OVER E(CD|  ) FINITE EXPANSION ADDRESSES
          ENDDO
C
C       CALCULATE DIRECTLY FROM LARGE RC(AB|CD) SCRATCH ARRAY
        ELSE
C
C         LOOP OVER ALL FINITE EXPANSION ADDRESSES FOR E(CD| -)
          DO ICD=1,NTUVCD
C
C           SKIP THIS STEP IF THE E(CD) FAILS SCREENING CONDITION
            ICDALL = ICDR11(ICD)+ICDI11(ICD)+ICDR21(ICD)+ICDI21(ICD)
            IF(ICDALL.EQ.0) GOTO 432
C
C           CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
            IRABCD = IABC(IA(IAB)+IA(ICD),IB(IAB)+IB(ICD),
     &                                    IC(IAB)+IC(ICD))
C
C           STARTING ADDRESS FOR THIS BATCH OF SAVED R(AB|CD) INTEGRALS
            IAD = (IJ-1)*MAXCD*NTUVABCDFL + MAXCD*(IRABCD-1)
C
C           SKIP THIS STEP IF THE RC(AB|CD) FAILS SCREENING CONDITION
            IF(IRCTTFL((IJ-1)*NTUVABCDFL+IRABCD).EQ.0) GOTO 432
C
C           CONTRIBUTIONS TO Re{G(AB|--)} FROM EACH Re{E(CD|--)} ADDRESS
            IF(IABR11(IJ,IAB).EQ.0.AND.ISYM.GT.1) GOTO 421
            IF(ICDR11(ICD).EQ.0) GOTO 421
            DO N=1,MAXN
              GABR11(N,IAB) = GABR11(N,IAB)
     &                  + DREAL(ECD11(IMAP(N),ICD))*RCTTFL(IAD+IMAP(N))
            ENDDO
421         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|--)} FROM EACH Im{E(CD|--)} ADDRESS
            IF(IABI11(IJ,IAB).EQ.0.AND.ISYM.GT.1) GOTO 422
            IF(ICDI11(ICD).EQ.0) GOTO 422
            DO N=1,MAXN
              GABI11(N,IAB) = GABI11(N,IAB)
     &                  + DIMAG(ECD11(IMAP(N),ICD))*RCTTFL(IAD+IMAP(N))
            ENDDO
422         CONTINUE
C
C           CONTRIBUTIONS TO Re{G(AB|+-)} FROM EACH Re{E(CD|+-)} ADDRESS
            IF(IABR21(IJ,IAB).EQ.0.AND.ISYM.GT.1) GOTO 423
            IF(ICDR21(ICD).EQ.0) GOTO 423
            DO N=1,MAXN
              GABR21(N,IAB) = GABR21(N,IAB)
     &                  + DREAL(ECD21(IMAP(N),ICD))*RCTTFL(IAD+IMAP(N))
            ENDDO
423         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|+-)} FROM EACH Im{E(CD|+-)} ADDRESS
            IF(IABI21(IJ,IAB).EQ.0.AND.ISYM.GT.1) GOTO 424
            IF(ICDI21(ICD).EQ.0) GOTO 424
            DO N=1,MAXN
              GABI21(N,IAB) = GABI21(N,IAB)
     &                  + DIMAG(ECD21(IMAP(N),ICD))*RCTTFL(IAD+IMAP(N))
            ENDDO
424         CONTINUE
C
C           SKIP POINT FOR RC(AB|CD) SCREENING
432         CONTINUE
C
C         END LOOP OVER E(CD|  ) FINITE EXPANSION ADDRESSES
          ENDDO
C
        ENDIF
C
C       SKIP POINT FOR E(AB|  ) SCREENING
401     CONTINUE
C
C     END LOOP OVER E(AB|  ) FINITE EXPANSION ADDRESSES
      ENDDO
C
C     TIME AT END OF FIRST CONTRACTION
      CALL SYSTEM_CLOCK(ICL2)
      TCC1 = TCC1 + DFLOAT(ICL2-ICL1)/RATE
C
C**********************************************************************C
C     PERFORM SECOND CONTRACTION: ( -| -) = E(AB| -)*G(AB| -).         C
C     THIS YIELDS A FULL BATCH OF TWO-ELECTRON INTEGRALS (16 PERM'NS). C
C**********************************************************************C
C
      CALL SYSTEM_CLOCK(ICL1)
C
C     CALCULATE PHASES FOR BASIS FUNCTION OVERLAP COMBINATIONS
      PAB = ISIGN(1,KQN(1)*KQN(2))*(-1)**((MQN(1)-MQN(2))/2)
      PCD = ISIGN(1,KQN(3)*KQN(4))*(-1)**((MQN(3)-MQN(4))/2)
C
      PABCD = PAB*PCD
C
C     SPECIAL CASE: LINEAR MOLECULE OR ATOM
      IF(ISYM.NE.1.AND.ISYM.NE.2) GOTO 501
C
C     1ST SET: ( 1) = (--|--)   ( 4) = (--|++)
C              (16) = (++|++)   (13) = (++|--)
C
C     RESET CONTRACTION STORAGE LISTS
      DO N=1,MAXN
        QR1(N) = 0.0D0
        QR2(N) = 0.0D0
      ENDDO
C
C     RAW CONTRACTION (--|--) = E(AB|--)*(Re{G(AB|--)} + i*Im{G(AB|--)})
      DO IAB=1,NTUVAB
        IF(IABR11(IJ,IAB).EQ.1) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB11(IJ,IAB))*GABR11(N,IAB)
          ENDDO
        ENDIF
        IF(IABI11(IJ,IAB).EQ.1) THEN
          DO N=1,MAXN
            QR2(N) = QR2(N) - DIMAG(EAB11(IJ,IAB))*GABI11(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     APPLY PHASE RELATIONS AND NORMALISATION FACTORS TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N, 1) =     DCMPLX(QR1(N)+QR2(N),0.0D0)
        RR(N, 4) = PCD*DCMPLX(QR1(N)-QR2(N),0.0D0)
        RR(N,16) = PABCD*DCONJG(RR(N, 1))
        RR(N,13) = PABCD*DCONJG(RR(N, 4))
      ENDDO
C
C     4TH SET: (11) = (+-|+-)   (10) = (+-|-+)
C              ( 6) = (-+|-+)   ( 7) = (-+|+-)
C
C     RESET CONTRACTION STORAGE LISTS
      DO N=1,MAXN
        QR1(N) = 0.0D0
        QR2(N) = 0.0D0
      ENDDO
C
C     RAW CONTRACTION (+-|+-) = E(AB|+-)*(Re{G(AB|+-)} + i*Im{G(AB|+-)})
      DO IAB=1,NTUVAB
        IF(IABR21(IJ,IAB).EQ.1) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB21(IJ,IAB))*GABR21(N,IAB)
          ENDDO
        ENDIF
        IF(IABI21(IJ,IAB).EQ.1) THEN
          DO N=1,MAXN
            QR2(N) = QR2(N) - DIMAG(EAB21(IJ,IAB))*GABI21(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     APPLY PHASE RELATIONS AND NORMALISATION FACTORS TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N,11) =     DCMPLX(QR1(N)+QR2(N),0.0D0)
        RR(N,10) =-PCD*DCMPLX(QR1(N)-QR2(N),0.0D0)
        RR(N, 6) = PABCD*DCONJG(RR(N,11))
        RR(N, 7) = PABCD*DCONJG(RR(N,10))
      ENDDO
C
      GOTO 502
C
C     GENERAL CASE
501   CONTINUE
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
        IF(IABR11(IJ,IAB).EQ.1) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB11(IJ,IAB))*GABR11(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB11(IJ,IAB))*GABI11(N,IAB)
          ENDDO
        ENDIF
        IF(IABI11(IJ,IAB).EQ.1) THEN
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
        IF(IABR11(IJ,IAB).EQ.1) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB11(IJ,IAB))*GABR21(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB11(IJ,IAB))*GABI21(N,IAB)
          ENDDO
        ENDIF
        IF(IABI11(IJ,IAB).EQ.1) THEN
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
        IF(IABR21(IJ,IAB).EQ.1) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB21(IJ,IAB))*GABR11(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB21(IJ,IAB))*GABI11(N,IAB)
          ENDDO
        ENDIF
        IF(IABI21(IJ,IAB).EQ.1) THEN
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
        IF(IABR21(IJ,IAB).EQ.1) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB21(IJ,IAB))*GABR21(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB21(IJ,IAB))*GABI21(N,IAB)
          ENDDO
        ENDIF
        IF(IABI21(IJ,IAB).EQ.1) THEN
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
502   CONTINUE
C
C     TIME AT END OF SECOND CONTRACTION
      CALL SYSTEM_CLOCK(ICL2)
      TCC2 = TCC2 + DFLOAT(ICL2-ICL1)/RATE
C
C**********************************************************************C
C     COULOMB INTEGRAL BATCH NOW FULLY CONSTRUCTED                     C
C**********************************************************************C
C
C     CALCULATE THE R-INTEGRAL NORMALISATION FACTOR
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
          IF(ISCR(M).EQ.1) THEN
            EKL = EXL(KBAS,3)+EXL(LBAS,4)
            EMX = DSQRT(EIJ+EKL)*EIJ*EKL
            PRE(M) = 2.0D0*PI52/EMX
          ENDIF
        ENDDO
      ENDDO
C
C     INCLUDE THE R-INTEGRAL NORMALISATION FACTOR
      DO N=1,MAXN
        DO ITG=1,16
          RR(N,ITG) = PRE(IMAP(N))*RR(N,ITG)
        ENDDO
      ENDDO
C
      RETURN
      END

