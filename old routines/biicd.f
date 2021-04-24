      SUBROUTINE BIICD(RR,XYZ,ICNT,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN,GAUNT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                  BBBBBBB IIII IIII CCCCCC  DDDDDDD                   C
C                  BB    BB II   II CC    CC DD    DD                  C
C                  BB    BB II   II CC       DD    DD                  C
C                  BBBBBBB  II   II CC       DD    DD                  C
C                  BB    BB II   II CC       DD    DD                  C
C                  BB    BB II   II CC    CC DD    DD                  C
C                  BBBBBBB IIII IIII CCCCCC  DDDDDDD                   C
C                                                                      C
C -------------------------------------------------------------------- C
C  BIICD GENERATES A BATCH OF MOLECULAR BREIT INTERACTION INTEGRALS BY C
C  MEANS OF THE MCMURCHIE-DAVIDSION ALGORITHM (DOUBLE FINITE SUM OVER  C
C  EQ-COEFFICIENTS AND INTEGRALS OVER A PAIR OF HGTFS.)                C
C -------------------------------------------------------------------- C
C  THIS IS A SPECIAL VERSION OF BII, WHICH CONTRACTS OVER THE (AB)     C
C  PAIRS BEFORE (CD), AND RECYCLES THE FIRST CONTRACTION OVER (CD).    C
C  THE ROUTINE EXPLICITLY MAKES USE OF EILS COEFFICIENTS ONLY, AND     C
C  CALLS THE COEFFICIENTS DIRECTLY FROM A GLOBAL FILE. IT CANNOT BE    C
C  INVOKED IF THE SET OF EILS OR R(AB|CD) BATCHES EXCEEDS MAX SPACE.   C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ XYZ    - FULL SET OF CARTESIAN BASIS CENTRES.                     C
C  ▶ ICNT   - FULL SET OF BASIS CENTRE ORIGINS (NUCLEAR LABELS).       C
C  ▶ KQN    - FULL SET OF RELATIVISTIC LABELS.                         C
C  ▶ MQN    - FULL SET OF MAGNETIC QUANTUM NUMBERS (MAGNITUDE).        C
C  ▶ NBAS   - FULL SET OF EXPONENT LIST LENGTHS.                       C
C  ▶ EXL    - FULL LISTS OF EXPONENTS IN THE BLOCK.                    C
C  ▶ IBAS   - 1ST BASIS FUNCTION (HELD CONSTANT DURING ROUTINE).       C
C  ▶ JBAS   - 2ND BASIS FUNCTION (HELD CONSTANT DURING ROUTINE).       C
C  ▶ ITN    - COMPONENT OVERLAP COMBINATION.                           C
C  ▶ GAUNT  - GAUNT INTERACTION OVERRIDE OPTION.                       C
C  OUTPUT:                                                             C
C  ▶ RR     - BII'S FOR BLOCK AB, ALL 16 MQN SIGN COMBINATIONS.        C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL GAUNT
C
      DIMENSION EXL(MBS,4),XYZ(3,4)
      DIMENSION ICNT(4),KQN(4),LQN(4),MQN(4),NBAS(4),ITN(2)
      DIMENSION PQ(MB2,3),APH(MB2),PRE(MB2),T(MB2),RC(MB2,MRC)
      DIMENSION IABR11(MB2,MEQ,3),IABI11(MB2,MEQ,3),
     &          IABR21(MB2,MEQ,3),IABI21(MB2,MEQ,3)
      DIMENSION ICDR11(MEQ,3,4),ICDI11(MEQ,3,4),
     &          ICDR21(MEQ,3,4),ICDI21(MEQ,3,4)
      DIMENSION ICDR11KM(4,MKP,MKP,MEQ,3,4),ICDI11KM(4,MKP,MKP,MEQ,3,4),
     &          ICDR21KM(4,MKP,MKP,MEQ,3,4),ICDI21KM(4,MKP,MKP,MEQ,3,4)
      DIMENSION IRC(MRC)
      DIMENSION RCTTFL(20*MFL),IRCTTFL(MFL)
      DIMENSION GCDR11(MB2,MEQ,3),GCDI11(MB2,MEQ,3),
     &          GCDR21(MB2,MEQ,3),GCDI21(MB2,MEQ,3)
      DIMENSION QR1(MB2),QI1(MB2),QR2(MB2),QI2(MB2)
C
      COMPLEX*16 RR(MB2,16)
C
      COMMON/EILS/EILSFL(MFL,12),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EISL/EISLFL(MFL,12),IADISL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/ICDB/ICDR11KM,ICDR11,ICDI11KM,ICDI11,
     &            ICDR21KM,ICDR21,ICDI21KM,ICDI21
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS,IABSL,ICDSL
      COMMON/IRCM/IEAB,IECD,NCD,IGAB,IRIJ(MBS,MBS)
      COMMON/ISCR/IMTX(MB2,11),ISCR(MB2),IMAP(MB2),IBCH,ITOG,MAXN
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
        IF(ITN(IT).NE.2) THEN
          WRITE(6, *) 'In BIICD: illegal component overlaps in ITN.'
          WRITE(7, *) 'In BIICD: illegal component overlaps in ITN.'
          STOP
        ENDIF
      ENDDO
C
C     EVALUATE LQNS FOR BASIS FUNCTIONS (A,B,C,D)
      DO N=1,4
        LQN(N) = LVAL(KQN(N))
      ENDDO
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
C     PHASE FACTOR FOR AB AND CD PAIR OVERLAPS
      IPHSAB = 1
      IPHSCD =-1
C
C     MCMURCHIE-DAVIDSON MAXIMUM ORDER FOR EQ-COEFFICIENTS
      LAMAB = LQN(1)+LQN(2)+1
      LAMCD = LQN(3)+LQN(4)+1
C
C     MCMURCHIE-DAVIDSON MAXIMUM ORDER FOR CONTRACTED R-INTEGRAL BATCH
      IF(GAUNT) THEN
C       GAUNT INTEGRALS ONLY
        LAMABCD = LAMAB+LAMCD
      ELSE
C       FULL BREIT INTERACTION
        LAMABCD = LAMAB+LAMCD+2
      ENDIF
C
C     MCMURCHIE-DAVIDSON EQ-COEFFICIENT AND R-INTEGRAL LIST LENGTHS
      NTUVAB   = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
      NTUVCD   = (LAMCD+1)*(LAMCD+2)*(LAMCD+3)/6
      NTUVABCD = (LAMABCD+1)*(LAMABCD+2)*(LAMABCD+3)/6
C
C     LIST ADDRESS FOR (AB|  ) AND GAUSSIAN EXPONENT FOR AB OVERLAP
      IJ  = (IBAS-1)*NBAS(2)+JBAS
      EIJ = EXL(IBAS,1)+EXL(JBAS,2)
C
C**********************************************************************C
C     GENERATE NEW BATCH OF E(AB|  ) SCREENING VALUES                  C
C**********************************************************************C
C
      CALL SYSTEM_CLOCK(ICL1,RATE)
      IF(IEAB.EQ.0) GOTO 100
C
C     SCREENING PROCEDURE: NORM SUM OF EQ-COEFFICIENT LIST FOR EACH IAB
      DO IX=1,3
        DO IAB=1,NTUVAB
          MAB = IABLS + (IAB-1)*MAXAB
C
C         Re{E(AB|--)} COEFFICIENTS
          DO M=1,MAXAB
            IF(DABS(EILSFL(MAB+M,4*(IX-1)+1)).GT.SENS) THEN
              IABR11(M,IAB,IX) = 1
            ELSE
              IABR11(M,IAB,IX) = 0
            ENDIF
          ENDDO
C
C         Im{E(AB|--)} COEFFICIENTS
          DO M=1,MAXAB
            IF(DABS(EILSFL(MAB+M,4*(IX-1)+2)).GT.SENS) THEN
              IABI11(M,IAB,IX) = 1
            ELSE
              IABI11(M,IAB,IX) = 0
            ENDIF
          ENDDO
C
C         Re{E(AB|+-)} COEFFICIENTS
          DO M=1,MAXAB
            IF(DABS(EILSFL(MAB+M,4*(IX-1)+3)).GT.SENS) THEN
              IABR21(M,IAB,IX) = 1
            ELSE
              IABR21(M,IAB,IX) = 0
            ENDIF
          ENDDO
C
C         Im{E(AB|+-)} COEFFICIENTS
          DO M=1,MAXAB
            IF(DABS(EILSFL(MAB+M,4*(IX-1)+4)).GT.SENS) THEN
              IABI21(M,IAB,IX) = 1
            ELSE
              IABI21(M,IAB,IX) = 0
            ENDIF
          ENDDO
C
        ENDDO
      ENDDO
C
C     DO NOT CALCULATE AGAIN UNTIL PROMPTED EXTERNALLY
      IEAB = 0
C
100   CONTINUE
C
      CALL SYSTEM_CLOCK(ICL2)
      TBEC = TBEC + DFLOAT(ICL2-ICL1)/RATE
C
C**********************************************************************C
C     BATCH OF E(CD|  ) SCREENING VALUES IS GENERATED ELSEWHERE        C
C**********************************************************************C
C
C**********************************************************************C
C     GENERATE NEW BATCH OF RC(AB|CD) INTEGRALS IF PROMPTED            C
C**********************************************************************C
C
C     START TIME
      CALL SYSTEM_CLOCK(ICL1,RATE)
C
C     GAUSSIAN OVERLAP CENTRE
      PX = (XYZ(1,1)*EXL(IBAS,1)+XYZ(1,2)*EXL(JBAS,2))/EIJ
      PY = (XYZ(2,1)*EXL(IBAS,1)+XYZ(2,2)*EXL(JBAS,2))/EIJ
      PZ = (XYZ(3,1)*EXL(IBAS,1)+XYZ(3,2)*EXL(JBAS,2))/EIJ
C
C     AUXILLIARY DATA FOR RMAKE ROUTINE
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M   = M+1
          EKL = EXL(KBAS,3)+EXL(LBAS,4)
          QX  = (XYZ(1,3)*EXL(KBAS,3)+XYZ(1,4)*EXL(LBAS,4))/EKL
          QY  = (XYZ(2,3)*EXL(KBAS,3)+XYZ(2,4)*EXL(LBAS,4))/EKL
          QZ  = (XYZ(3,3)*EXL(KBAS,3)+XYZ(3,4)*EXL(LBAS,4))/EKL
          PQ(M,1) = QX-PX
          PQ(M,2) = QY-PY
          PQ(M,3) = QZ-PZ
          APH(M)  = EIJ*EKL/(EIJ+EKL)
        ENDDO
      ENDDO
C
C     SKIP IF INTEGRAL BATCH EXISTS IN FILE
      IF(RCFILE.AND.IRIJ(IBAS,JBAS).EQ.0) GOTO 300
C
      CALL SYSTEM_CLOCK(ICL3)
C
C     BATCH SIZE AND EXPANSION LENGTH
      MBCH = MAXCD
C
C     GENERATE R-INTEGRALS
      CALL RMAKE(RC,PQ,APH,MBCH,LAMABCD)
C
C     SCREENING: TEST RC(AB|CD) COLUMNS WITH INDEX (T+T',U+U',V+V')
      DO IABCD=1,NTUVABCD
C
C       SUM OF RC(AB|CD) MAGNITUDES
        SUM = 0.0D0
        DO N=1,MBCH
          SUM = SUM + DABS(RC(N,IABCD))
          IF(SUM.GT.SENS) THEN
            IRC(IABCD) = 1
            GOTO 301
          ENDIF
        ENDDO
        IRC(IABCD) = 0
301     CONTINUE
C
      ENDDO
C
      CALL SYSTEM_CLOCK(ICL4)
      TBRM = TBRM + DFLOAT(ICL4-ICL3)/RATE
C
      CALL SYSTEM_CLOCK(ICL3,RATE)
C
C     TEST WHETHER FINAL ADDRESS IS STILL INSIDE ARRAY BOUNDS
      IF(20*MFL.LT.IJ*MAXCD*NTUVABCD) THEN
C       OUT OF BOUNDS: PRINT WARNING BUT KEEP GOING
        WRITE(6, *) 'In BIICD: RCTT words exceed allocated limit.'
        WRITE(7, *) 'In BIICD: RCTT words exceed allocated limit.'
        GOTO 300
      ELSE
C       DO NOT CALCULATE AGAIN UNTIL PROMPTED EXTERNALLY
        IRIJ(IBAS,JBAS) = 0
      ENDIF
C
C     STARTING ADDRESS FOR THIS BATCH OF SAVED R(AB|CD) INTEGRALS
      IADRTT = (IJ-1)*MAXCD*NTUVABCD
C
C     COPY THIS BATCH OF INTEGRALS TO A SAVED LIST
      DO IABCD=1,NTUVABCD
        IAD = IADRTT + MAXCD*(IABCD-1)
        DO M=1,MAXCD
          RCTTFL(IAD+M) = RC(M,IABCD)
        ENDDO
      ENDDO
C
C     STARTING ADDRESS FOR SCREENING FLAGS
      IADSCR = (IJ-1)*NTUVABCD
C
C     COPY SCREENING MARKERS TO A SAVED LIST
      DO IABCD=1,NTUVABCD
        IRCTTFL(IADSCR+IABCD) = IRC(IABCD)
      ENDDO
C
      CALL SYSTEM_CLOCK(ICL4)
      TBRW = TBRW + DFLOAT(ICL4-ICL3)/RATE
C
300   CONTINUE
C
C     RECORD THE TIME TAKEN TO GENERATE THE RC(AB|CD) BATCH
      CALL SYSTEM_CLOCK(ICL2)
      TRBR = TRBR + DFLOAT(ICL2-ICL1)/RATE
C
C**********************************************************************C
C     PERFORM FIRST CONTRACTION: G(CD| -) = E(AB| -)*RC(AB|CD).        C
C     THIS YIELDS ALL MQN SIGN POSSIBILITIES FOR A AND B.              C
C**********************************************************************C
C
      IF(IGAB.EQ.0) GOTO 400
C
C     LOOP OVER CARTESIAN INDEX IX FOR CENTRE CD
      DO IX=1,3
C
C       TIME AT START OF FIRST CONTRACTION FOR THIS IX INDEX
        CALL SYSTEM_CLOCK(ICL1,RATE)
C
C       LOOP OVER ALL ADDRESSES FOR E(CD| -) FINITE EXPANSION
        DO ICD=1,NTUVCD
C
C         RESET CONTRACTION STORAGE ARRAYS G(CD| -)
          DO N=1,MAXN
            GCDR11(N,ICD,IX) = 0.0D0
            GCDI11(N,ICD,IX) = 0.0D0
            GCDR21(N,ICD,IX) = 0.0D0
            GCDI21(N,ICD,IX) = 0.0D0
          ENDDO
C
C         SKIP ENTIRE PROCESS IF E(CD| -) FAILS SCREENING CONDITION
          ICDALL = ICDR11(ICD,IX,2) + ICDI11(ICD,IX,2)
     &           + ICDR21(ICD,IX,2) + ICDI21(ICD,IX,2)
          IF(ICDALL.EQ.0) GOTO 401
C
C >>>>>   GAUNT INTERACTION
C
C         LOOP OVER ALL FINITE EXPANSION ADDRESSES FOR E(AB| -)
          DO IAB=1,NTUVAB
C
C           RELATIVE EQ(AB) LIST STARTING ADDRESS
            MAB = IABLS + (IAB-1)*MAXAB
C
C           CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
            IRABCD = IABC(IA(IAB)+IA(ICD),IB(IAB)+IB(ICD),
     &                                    IC(IAB)+IC(ICD))
C
C           ARE ANY EQ-COEFFICIENTS AT THIS ADDRESS NON-ZERO?
            IABALL = IABR11(IJ,IAB,IX) + IABI11(IJ,IAB,IX)
     &             + IABR21(IJ,IAB,IX) + IABI21(IJ,IAB,IX)
C
C           SKIP THIS STEP IF THE E(AB) FAILS SCREENING CONDITION
            IF(IABALL.EQ.0) GOTO 402
C
C           RELATIVE RC(AB|CD) LIST STARTING ADDRESS
            MABCD = (IJ-1)*MAXCD*NTUVABCD + MAXCD*(IRABCD-1)
C
C           SKIP THIS STEP IF THE RC(AB|CD) FAILS SCREENING CONDITION
            IF(IRCTTFL((IJ-1)*NTUVABCD+IRABCD).EQ.0) GOTO 450
C
C           CONTRIBUTIONS TO Re{G(CD|--)} FROM EACH Re{E(AB|--)}
            IF(ICDR11(ICD,IX,2).EQ.0.AND.ISYM.GE.1) GOTO 451
            IF(IABR11(IJ,IAB,IX).EQ.0) GOTO 451
            DO N=1,MAXN
              GCDR11(N,ICD,IX) = GCDR11(N,ICD,IX)
     &                + EILSFL(MAB+IJ,4*(IX-1)+1)*RCTTFL(MABCD+IMAP(N))
            ENDDO
451         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(CD|--)} FROM EACH Im{E(AB|--)}
            IF(ICDI11(ICD,IX,2).EQ.0.AND.ISYM.GE.1) GOTO 452
            IF(IABI11(IJ,IAB,IX).EQ.0) GOTO 452
            DO N=1,MAXN
              GCDI11(N,ICD,IX) = GCDI11(N,ICD,IX)
     &                + EILSFL(MAB+IJ,4*(IX-1)+2)*RCTTFL(MABCD+IMAP(N))
            ENDDO
452         CONTINUE
C
C           CONTRIBUTIONS TO Re{G(CD|+-)} FROM EACH Re{E(AB|+-)}
            IF(ICDR21(ICD,IX,2).EQ.0.AND.ISYM.GE.1) GOTO 453
            IF(IABR21(IJ,IAB,IX).EQ.0) GOTO 453
            DO N=1,MAXN
              GCDR21(N,ICD,IX) = GCDR21(N,ICD,IX)
     &                + EILSFL(MAB+IJ,4*(IX-1)+3)*RCTTFL(MABCD+IMAP(N))
            ENDDO
453         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(CD|+-)} FROM EACH Im{E(AB|+-)}
            IF(ICDI21(ICD,IX,2).EQ.0.AND.ISYM.GE.1) GOTO 454
            IF(IABI21(IJ,IAB,IX).EQ.0) GOTO 454
            DO N=1,MAXN
              GCDI21(N,ICD,IX) = GCDI21(N,ICD,IX)
     &                + EILSFL(MAB+IJ,4*(IX-1)+4)*RCTTFL(MABCD+IMAP(N))
            ENDDO
454         CONTINUE
C
C           SKIP POINT FOR RC(AB|CD) AND E(AB) SCREENING
450         CONTINUE
C
C           SKIP POINT FOR ALL EQ(AB) FAILING SCREENING CONDITION
402         CONTINUE
C
C >>>>>     GAUGE TERM (REQUIRES ADDITIONAL CARTESIAN SUM Q')
            IF(GAUNT) GOTO 480
C
C           LOOP OVER CARTESIAN INDEX JX FOR CENTRE AB
            DO JX=1,3
C
C             SKIP THIS STEP IF THE E(AB) FAILS SCREENING CONDITION
              IF(IABR11(IJ,IAB,JX)+IABI11(IJ,IAB,JX)
     &          +IABR21(IJ,IAB,JX)+IABI21(IJ,IAB,JX).EQ.0) THEN
                GOTO 403
              ENDIF
C
C             NEW ADDRESS DEPENDING ON JX CARTESIAN INDEX
              IF(JX.EQ.1) THEN
                RTP = DFLOAT(IA(IAB)+IA(ICD))
              ELSEIF(JX.EQ.2) THEN
                RTP = DFLOAT(IB(IAB)+IB(ICD))
              ELSEIF(JX.EQ.3) THEN
                RTP = DFLOAT(IC(IAB)+IC(ICD))
              ENDIF
C
C             FIRST CONTRIBUTION ADDRESS
              I1 = IA(IAB)+IA(ICD)+KRONECK(IX,1)+KRONECK(JX,1)
              J1 = IB(IAB)+IB(ICD)+KRONECK(IX,2)+KRONECK(JX,2)
              K1 = IC(IAB)+IC(ICD)+KRONECK(IX,3)+KRONECK(JX,3)
C
C             CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
              IADR1 = IABC(I1,J1,K1)
C
C             SECOND CONTRIBUTION ADDRESS
              I2 = IA(IAB)+IA(ICD)+KRONECK(IX,1)
              J2 = IB(IAB)+IB(ICD)+KRONECK(IX,2)
              K2 = IC(IAB)+IC(ICD)+KRONECK(IX,3)
C
C             CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
              IADR2 = IABC(I2,J2,K2)
C
C             THIRD CONTRIBUTION ADDRESS
              I3 = IA(IAB)+IA(ICD)+KRONECK(IX,1)-KRONECK(JX,1)
              J3 = IB(IAB)+IB(ICD)+KRONECK(IX,2)-KRONECK(JX,2)
              K3 = IC(IAB)+IC(ICD)+KRONECK(IX,3)-KRONECK(JX,3)
C
C             CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
              IF(I3.GE.0.AND.J3.GE.0.AND.K3.GE.0) THEN
                IADR3 = IABC(I3,J3,K3)
              ELSE
                IADR3 = 0
              ENDIF
C
C             SKIP THIS STEP IF RC(AB|CD) FAILS SCREENING CONDITION
              IA1 = (IJ-1)*MAXCD*NTUVABCD + MAXCD*(IADR1-1)
              IA2 = (IJ-1)*MAXCD*NTUVABCD + MAXCD*(IADR2-1)
              IA3 = (IJ-1)*MAXCD*NTUVABCD + MAXCD*(IADR3-1)
              IF(IADR3.NE.0) THEN
                IF(IA1+IA2+IA3.EQ.0) GOTO 403
              ELSE
                IF(IA1+IA2.EQ.0) GOTO 403
              ENDIF
C
C             PRE-FACTORS FOR THE UPCOMING CONTRACTION
              DO N=1,MAXN
                T1 = RCTTFL(IA1+IMAP(N))*0.5D0/APH(IMAP(N))
                T2 = RCTTFL(IA2+IMAP(N))*PQ(IMAP(N),JX)
                IF(I3.GE.0.AND.J3.GE.0.AND.K3.GE.0) THEN
                  T3 = RCTTFL(IA3+IMAP(N))*RTP
                ELSE
                  T3 = 0.0D0
                ENDIF
                T(N) = T1-T2+T3
              ENDDO
C
C             CONTRIBUTIONS TO Re{G(CD|--)} FROM EACH Re{E(AB|--)}
              IF(ICDR11(ICD,IX,2).EQ.0.AND.ISYM.GE.1) GOTO 415
              IF(IABR11(IJ,IAB,JX).EQ.0) GOTO 415
              DO N=1,MAXN
                GCDR11(N,ICD,IX) = GCDR11(N,ICD,IX)
     &                                 - EILSFL(MAB+IJ,4*(JX-1)+1)*T(N)
              ENDDO
415           CONTINUE
C
C             CONTRIBUTIONS TO Im{G(CD|--)} FROM EACH Im{E(AB|--)}
              IF(ICDI11(ICD,IX,2).EQ.0.AND.ISYM.GE.1) GOTO 416
              IF(IABI11(IJ,IAB,JX).EQ.0) GOTO 416
              DO N=1,MAXN
                GCDI11(N,ICD,IX) = GCDI11(N,ICD,IX)
     &                                 - EILSFL(MAB+IJ,4*(JX-1)+2)*T(N)
              ENDDO
416           CONTINUE
C
C             CONTRIBUTIONS TO Re{G(CD|--)} FROM EACH Re{E(AB|--)}
              IF(ICDR21(ICD,IX,2).EQ.0.AND.ISYM.GE.1) GOTO 417
              IF(IABR21(IJ,IAB,JX).EQ.0) GOTO 417
              DO N=1,MAXN
                GCDR21(N,ICD,IX) = GCDR21(N,ICD,IX)
     &                                 - EILSFL(MAB+IJ,4*(JX-1)+3)*T(N)
              ENDDO
417           CONTINUE
C
C             CONTRIBUTIONS TO Im{G(CD|--)} FROM EACH Im{E(AB|--)}
              IF(ICDI21(ICD,IX,2).EQ.0.AND.ISYM.GE.1) GOTO 418
              IF(IABI21(IJ,IAB,JX).EQ.0) GOTO 418
              DO N=1,MAXN
                GCDI21(N,ICD,IX) = GCDI21(N,ICD,IX)
     &                                 - EILSFL(MAB+IJ,4*(JX-1)+4)*T(N)
              ENDDO
418           CONTINUE
C
C             SKIP POINT FOR E(AB) SCREENING
403           CONTINUE
C
C             END LOOP OVER CARTESIAN INDEX JX FOR CENTRE AB
              ENDDO
C
C             SKIP POINT FOR GAUNT INTERACTION ONLY
480           CONTINUE
C
C           END LOOP OVER E(AB|  ) FINITE EXPANSION ADDRESSES
            ENDDO
C
C         SKIP POINT FOR E(CD|  ) SCREENING
401       CONTINUE
C
C       END LOOP OVER E(CD|  ) FINITE EXPANSION ADDRESSES
        ENDDO
C
C       TIME AT END OF FIRST CONTRACTION FOR THIS IX INDEX
        CALL SYSTEM_CLOCK(ICL2)
        TBC1 = TBC1 + DFLOAT(ICL2-ICL1)/RATE
C
C     END LOOP OVER CARTESIAN INDEX IX FOR CENTRE CD
      ENDDO
C
C     FIRST CONTRACTION DOES NOT NEED TO BE RECALCULATED
      IGAB = 0
400   CONTINUE
C
C**********************************************************************C
C     PERFORM SECOND CONTRACTION: ( -| -) = E(CD| -)*G(CD| -).         C
C     THIS YIELDS A FULL BATCH OF TWO-ELECTRON INTEGRALS (16 PERM'NS). C
C**********************************************************************C
C
C     INITIALISE RR ARRAY
      DO M=1,MAXCD
        DO ITG=1,16
          RR(M,ITG) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     CALCULATE PHASES FOR BASIS FUNCTION OVERLAP COMBINATIONS
      PAB =-ISIGN(1,KQN(1)*KQN(2))*(-1)**((MQN(1)-MQN(2))/2)
      PCD =-ISIGN(1,KQN(3)*KQN(4))*(-1)**((MQN(3)-MQN(4))/2)
C
      PABCD = PAB*PCD
C
C     SHORTHAND FOR MQN VALUES IN (CD) BLOCK
      MC = (MQN(3)+1)/2
      MD = (MQN(4)+1)/2
C
C     LOOP OVER CARTESIAN INDEX IX FOR CENTRE CD
      DO IX=1,3
C
        CALL SYSTEM_CLOCK(ICL1,RATE)
C
C       SPECIAL CASE: LINEAR MOLECULE OR ATOM
        IF(ISYM.NE.1.AND.ISYM.NE.2) GOTO 501
C
C       1ST SET: ( 1) = (--|--)   ( 4) = (--|++)
C                (16) = (++|++)   (13) = (++|--)
C
C       RESET CONTRACTION STORAGE LISTS
        DO N=1,MAXN
          QR1(N) = 0.0D0
          QR2(N) = 0.0D0
        ENDDO
C
C       RAW CONTRACTION (--|--) = E(CD|--)*(Re{G(CD|--)}+i*Im{G(CD|--)})
        DO ICD=1,NTUVCD
          MCD = ICDLS + (ICD-1)*MAXCD
          Z   = DFLOAT((-1)**(ILAM(ICD)))
          IF(ICDR11KM(NCD,MC,MD,ICD,IX,2).NE.0) THEN
            DO N=1,MAXN
              QR1(N) = QR1(N)
     &              + Z*EILSFL(MCD+IMAP(N),4*(IX-1)+1)*GCDR11(N,ICD,IX)
            ENDDO
          ENDIF
          IF(ICDI11KM(NCD,MC,MD,ICD,IX,2).NE.0) THEN
            DO N=1,MAXN
              QR2(N) = QR2(N)
     &              - Z*EILSFL(MCD+IMAP(N),4*(IX-1)+2)*GCDI11(N,ICD,IX)
            ENDDO
          ENDIF
        ENDDO
C
C       ADD THIS IX TERM TO RAW CONTRACTION
        DO N=1,MAXN
          RR(N,1 ) = RR(N,1 ) +     DCMPLX(QR1(N)+QR2(N),0.0D0)
          RR(N,13) = RR(N,13) + PAB*DCMPLX(QR1(N)-QR2(N),0.0D0)
        ENDDO
C
C       4TH SET: (11) = (+-|+-)   (10) = (+-|-+)
C                ( 6) = (-+|-+)   ( 7) = (-+|+-)
C
C       RESET CONTRACTION STORAGE LISTS
        DO N=1,MAXN
          QR1(N) = 0.0D0
          QR2(N) = 0.0D0
        ENDDO
C
C       RAW CONTRACTION (+-|+-) = E(CD|+-)*(Re{G(CD|+-)}+i*Im{G(CD|+-)})
        DO ICD=1,NTUVCD
          MCD = ICDLS + (ICD-1)*MAXCD
          Z   = DFLOAT((-1)**(ILAM(ICD)))
          IF(ICDR21KM(NCD,MC,MD,ICD,IX,2).NE.0) THEN
            DO N=1,MAXN
              QR1(N) = QR1(N)
     &              + Z*EILSFL(MCD+IMAP(N),4*(IX-1)+3)*GCDR21(N,ICD,IX)
            ENDDO
          ENDIF
          IF(ICDI21KM(NCD,MC,MD,ICD,IX,2).NE.0) THEN
            DO N=1,MAXN
              QR2(N) = QR2(N)
     &              - Z*EILSFL(MCD+IMAP(N),4*(IX-1)+4)*GCDI21(N,ICD,IX)
            ENDDO
          ENDIF
        ENDDO
C
C       ADD THIS IX TERM TO RAW CONTRACTION
        DO N=1,MAXN
          RR(N,11) = RR(N,11) +     DCMPLX(QR1(N)+QR2(N),0.0D0)
          RR(N,7 ) = RR(N,7 ) - PAB*DCMPLX(QR1(N)-QR2(N),0.0D0)
        ENDDO
C
        GOTO 502
C
C       GENERAL CASE
501     CONTINUE
C
C       1ST SET: ( 1) = (--|--)   ( 4) = (--|++)
C                (16) = (++|++)   (13) = (++|--)
C
C       RESET CONTRACTION STORAGE LISTS
        DO N=1,MAXN
          QR1(N) = 0.0D0
          QI1(N) = 0.0D0
          QR2(N) = 0.0D0
          QI2(N) = 0.0D0
        ENDDO
C
C       RAW CONTRACTION (--|--) = E(CD|--)*(Re{G(CD|--)}+i*Im{G(CD|--)})
        DO ICD=1,NTUVCD
          MCD = ICDLS + (ICD-1)*MAXCD
          Z   = DFLOAT((-1)**(ILAM(ICD)))
          IF(ICDR11KM(NCD,MC,MD,ICD,IX,2).NE.0) THEN
            DO N=1,MAXN
              QR1(N) = QR1(N)
     &              + Z*EILSFL(MCD+IMAP(N),4*(IX-1)+1)*GCDR11(N,ICD,IX)
              QI2(N) = QI2(N)
     &              + Z*EILSFL(MCD+IMAP(N),4*(IX-1)+1)*GCDI11(N,ICD,IX)
            ENDDO
          ENDIF
          IF(ICDI11KM(NCD,MC,MD,ICD,IX,2).NE.0) THEN
            DO N=1,MAXN
              QI1(N) = QI1(N)
     &              + Z*EILSFL(MCD+IMAP(N),4*(IX-1)+2)*GCDR11(N,ICD,IX)
              QR2(N) = QR2(N)
     &              - Z*EILSFL(MCD+IMAP(N),4*(IX-1)+2)*GCDI11(N,ICD,IX)
            ENDDO
          ENDIF
        ENDDO
C
C       ADD THIS IX TERM TO RAW CONTRACTION
        DO N=1,MAXN
          RR(N,1 ) = RR(N,1 ) +     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
          RR(N,13) = RR(N,13) + PAB*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
        ENDDO
C
C       2ND SET: ( 3) = (--|+-)   ( 2) = (--|-+)
C                (14) = (++|-+)   (15) = (++|+-)
C
C       RESET CONTRACTION STORAGE LISTS
        DO N=1,MAXN
          QR1(N) = 0.0D0
          QI1(N) = 0.0D0
          QR2(N) = 0.0D0
          QI2(N) = 0.0D0
        ENDDO
C
C       RAW CONTRACTION (--|+-) = E(CD|--)*(Re{G(CD|+-)}+i*Im{G(CD|+-)})
        DO ICD=1,NTUVCD
          MCD = ICDLS + (ICD-1)*MAXCD
          Z   = DFLOAT((-1)**(ILAM(ICD)))
          IF(ICDR21KM(NCD,MC,MD,ICD,IX,2).NE.0) THEN
            DO N=1,MAXN
              QR1(N) = QR1(N)
     &              + Z*EILSFL(MCD+IMAP(N),4*(IX-1)+3)*GCDR11(N,ICD,IX)
              QI2(N) = QI2(N)
     &              + Z*EILSFL(MCD+IMAP(N),4*(IX-1)+3)*GCDI11(N,ICD,IX)
            ENDDO
          ENDIF
          IF(ICDI21KM(NCD,MC,MD,ICD,IX,2).NE.0) THEN
            DO N=1,MAXN
              QI1(N) = QI1(N)
     &              + Z*EILSFL(MCD+IMAP(N),4*(IX-1)+4)*GCDR11(N,ICD,IX)
              QR2(N) = QR2(N)
     &              - Z*EILSFL(MCD+IMAP(N),4*(IX-1)+4)*GCDI11(N,ICD,IX)
            ENDDO
          ENDIF
        ENDDO
C
C       ADD THIS IX TERM TO RAW CONTRACTION
        DO N=1,MAXN
          RR(N,3 ) = RR(N,3 ) +     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
          RR(N,15) = RR(N,15) + PAB*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
        ENDDO
C
C       3RD SET: ( 9) = (+-|--)   (12) = (+-|++)
C                ( 8) = (-+|++)   ( 5) = (-+|--)
C
C       RESET CONTRACTION STORAGE LISTS
        DO N=1,MAXN
          QR1(N) = 0.0D0
          QI1(N) = 0.0D0
          QR2(N) = 0.0D0
          QI2(N) = 0.0D0
        ENDDO
C
C       RAW CONTRACTION (+-|--) = E(CD|+-)*(Re{G(CD|--)}+i*Im{G(CD|--)})
        DO ICD=1,NTUVCD
          MCD = ICDLS + (ICD-1)*MAXCD
          Z   = DFLOAT((-1)**(ILAM(ICD)))
          IF(ICDR11KM(NCD,MC,MD,ICD,IX,2).NE.0) THEN
            DO N=1,MAXN
              QR1(N) = QR1(N)
     &              + Z*EILSFL(MCD+IMAP(N),4*(IX-1)+1)*GCDR21(N,ICD,IX)
              QI2(N) = QI2(N)
     &              + Z*EILSFL(MCD+IMAP(N),4*(IX-1)+1)*GCDI21(N,ICD,IX)
            ENDDO
          ENDIF
          IF(ICDI11KM(NCD,MC,MD,ICD,IX,2).NE.0) THEN
            DO N=1,MAXN
              QI1(N) = QI1(N)
     &              + Z*EILSFL(MCD+IMAP(N),4*(IX-1)+2)*GCDR21(N,ICD,IX)
              QR2(N) = QR2(N)
     &              - Z*EILSFL(MCD+IMAP(N),4*(IX-1)+2)*GCDI21(N,ICD,IX)
            ENDDO
          ENDIF
        ENDDO
C
C       ADD THIS IX TERM TO RAW CONTRACTION
        DO N=1,MAXN
          RR(N,9 ) = RR(N,9 ) +     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
          RR(N,5 ) = RR(N,5 ) - PAB*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
        ENDDO
C
C       4TH SET: (11) = (+-|+-)   (10) = (+-|-+)
C                ( 6) = (-+|-+)   ( 7) = (-+|+-)
C
C       RESET CONTRACTION STORAGE LISTS
        DO N=1,MAXN
          QR1(N) = 0.0D0
          QI1(N) = 0.0D0
          QR2(N) = 0.0D0
          QI2(N) = 0.0D0
        ENDDO
C
C       RAW CONTRACTION (+-|+-) = E(CD|+-)*(Re{G(CD|+-)}+i*Im{G(CD|+-)})
        DO ICD=1,NTUVCD
          MCD = ICDLS + (ICD-1)*MAXCD
          Z   = DFLOAT((-1)**(ILAM(ICD)))
          IF(ICDR21KM(NCD,MC,MD,ICD,IX,2).NE.0) THEN
            DO N=1,MAXN
              QR1(N) = QR1(N)
     &              + Z*EILSFL(MCD+IMAP(N),4*(IX-1)+3)*GCDR21(N,ICD,IX)
              QI2(N) = QI2(N)
     &              + Z*EILSFL(MCD+IMAP(N),4*(IX-1)+3)*GCDI21(N,ICD,IX)
            ENDDO
          ENDIF
          IF(ICDI21KM(NCD,MC,MD,ICD,IX,2).NE.0) THEN
            DO N=1,MAXN
              QI1(N) = QI1(N)
     &              + Z*EILSFL(MCD+IMAP(N),4*(IX-1)+4)*GCDR21(N,ICD,IX)
              QR2(N) = QR2(N)
     &              - Z*EILSFL(MCD+IMAP(N),4*(IX-1)+4)*GCDI21(N,ICD,IX)
            ENDDO
          ENDIF
        ENDDO
C
C       ADD THIS IX TERM TO RAW CONTRACTION
        DO N=1,MAXN
          RR(N,11) = RR(N,11) +     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
          RR(N,7 ) = RR(N,7 ) - PAB*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
        ENDDO
C
502     CONTINUE
C
C       TIME AT END OF SECOND CONTRACTION FOR THIS IX INDEX
        CALL SYSTEM_CLOCK(ICL2)
        TBC2 = TBC2 + DFLOAT(ICL2-ICL1)/RATE
C
C     END LOOP OVER CARTESIAN INDEX IX FOR CENTRE CD
      ENDDO
C
C     HALF OF THE RR ARRAY CAN BE GENERATED WITH PHASE RELATIONS
      DO N=1,MAXN
        RR(N,16) = PABCD*DCONJG(RR(N,1 ))
        RR(N,4 ) = PABCD*DCONJG(RR(N,13))
        RR(N,14) =-PABCD*DCONJG(RR(N,3 ))
        RR(N,2 ) =-PABCD*DCONJG(RR(N,15))
        RR(N,8 ) =-PABCD*DCONJG(RR(N,9 ))
        RR(N,12) =-PABCD*DCONJG(RR(N,5 ))
        RR(N,6 ) = PABCD*DCONJG(RR(N,11))
        RR(N,10) = PABCD*DCONJG(RR(N,7 ))
      ENDDO 
C
C**********************************************************************C
C     BREIT INTEGRAL BATCH NOW FULLY CONSTRUCTED                       C
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
C     INCLUDE THE OUTSIDE FACTOR OF (-1/2) AND MOVE TO FULL ARRAY
      IF(GAUNT) THEN
        BFAC =-1.0D0
      ELSE
        BFAC =-0.5D0
      ENDIF
C
C     INCLUDE THE OUTSIDE FACTOR OF (-1/2) AND MOVE TO FULL ARRAY
      DO N=1,MAXN
        DO ITG=1,16
          RR(N,ITG) = BFAC*PRE(IMAP(N))*RR(N,ITG)
        ENDDO
      ENDDO
C
      RETURN
      END

