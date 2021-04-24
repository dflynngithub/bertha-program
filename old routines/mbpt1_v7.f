      SUBROUTINE MBPT1(MINO,NUMO,G2INT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C              MM       MM BBBBBBB  PPPPPPP TTTTTTTT 11                C
C              MMM     MMM BB    BB PP    PP   TT   111                C
C              MMMM   MMMM BB    BB PP    PP   TT    11                C
C              MM MM MM MM BBBBBBB  PP    PP   TT    11                C
C              MM  MMM  MM BB    BB PPPPPPP    TT    11                C
C              MM   M   MM BB    BB PP         TT    11                C
C              MM       MM BBBBBBB  PP         TT   1111               C
C                                                                      C
C -------------------------------------------------------------------- C
C  MBPT1 EVALUATES ZERO- AND FIRST-ORDER ENERGIES FOR ALL OCCUPIED     C
C  SOLUTIONS TO A CONVERGED HARTREE-FOCK PROBLEM.                      C
C -------------------------------------------------------------------- C
C INPUT:                                                               C
C  ▶ MINO  - LOWEST OCCUPIED STATE TO ACCOUNT FOR. (FULL: 1)           C
C  ▶ NUMO  - NUMBER OF OCCUPIED STATES TO ACCOUNT FOR. (FULL: NOCC)    C
C  ▶ G2INT - NAME OF TWO-BODY OPERATOR ('COULM' OR 'BREIT').           C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL BGAUNT
C
      CHARACTER*5 G2INT
      CHARACTER*16 HMS
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 TITLE
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4),LQN(4),JQN(4)
      DIMENSION INDEX(MCT,-(MEL+1):MEL,MKP),ITN(2)
      DIMENSION EAB1(NUMO,NUMO,6),EA1(NUMO,6)
C
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 B1(MBS*NUMO,8),B2(MBS*NUMO,8)
      COMPLEX*16 DB(MB2,NUMO*NUMO,4)
      COMPLEX*16 ADB1(MBS,NUMO*NUMO*NUMO,2),ADB2(MBS,NUMO*NUMO*NUMO,2)
      COMPLEX*16 CADB((NUMO+1)*NUMO*NUMO*NUMO/2)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VUEH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),WDIR(MDM,MDM),
     &           WXCH(MDM,MDM),CPLE(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/E0LL/E0LLFL(MFL,4),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,4),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EILS/EILSFL(MFL,12),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EISL/EISLFL(MFL,12),IADISL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS,IABSL,ICDSL
      COMMON/IRCM/IEAB,IECD,IRIJ(MBS,MBS)
      COMMON/ISCR/IMTX(MB2,11),ISCR(MB2),IMAP(MB2),IBCH,ITOG,MAXN
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VUEH,QDIR,
     &            QXCH,WDIR,WXCH,CPLE
C
C     TURN OFF RC(AB|CD) LOCAL FILE PROCESS
      RCFILE = .FALSE.
C
C     WARNINGS BASED ON INVALID HMLT VS. G2INT COMBINATIONS
      IF(G2INT.EQ.'COULM') THEN
        IF(HMLT.EQ.'BARE') THEN
          WRITE(6, *) 'In MBPT1: HMLT = BARE but G2INT = COULM.'
          WRITE(7, *) 'In MBPT1: HMLT = BARE but G2INT = COULM.'
        ENDIF
      ELSEIF(G2INT.EQ.'BREIT') THEN
        IF(HMLT.EQ.'NORL') THEN
          WRITE(6, *) 'In MBPT1: HMLT = NORL but G2INT = BREIT.'
          WRITE(7, *) 'In MBPT1: HMLT = NORL but G2INT = BREIT.'
          RETURN
        ELSEIF(HMLT.EQ.'BARE') THEN
          WRITE(6, *) 'In MBPT1: HMLT = BARE but G2INT = BREIT.'
          WRITE(7, *) 'In MBPT1: HMLT = BARE but G2INT = BREIT.'
        ELSEIF(HMLT.EQ.'DHFR') THEN
          WRITE(6, *) 'In MBPT1: HMLT = DHFR but G2INT = BREIT.'
          WRITE(7, *) 'In MBPT1: HMLT = DHFR but G2INT = BREIT.'
        ENDIF
      ENDIF
C
C     COMPONENT OVERLAP LABELS TO LOOP OVER
      IF(HMLT.EQ.'BARE') THEN
        RETURN
      ELSEIF(HMLT.EQ.'NORL') THEN
        IF(G2INT.EQ.'COULM') THEN
          ITSTRT = 1
          ITSTOP = 1
          ITSKIP = 1
        ENDIF
      ELSE
        IF(G2INT.EQ.'COULM') THEN
          ITSTRT = 4
          ITSTOP = 1
          ITSKIP =-3
        ELSEIF(G2INT.EQ.'BREIT') THEN
          ITSTRT = 2
          ITSTOP = 3
          ITSKIP = 1        
        ENDIF
      ENDIF
C
C     INITIALISE TIME COUNTERS
      T1EL = 0.0D0
      TERI = 0.0D0
      TCN1 = 0.0D0
      TCN2 = 0.0D0
      TCN3 = 0.0D0
      TCN4 = 0.0D0
      TSUM = 0.0D0
C
      CALL CPU_TIME(TBEG)
C
C     CLEAR ENERGY COUNTERS
      DO N=1,6
        DO IOCCB=1,NUMO
          EA1(IOCCB,N) = 0.0D0
          DO IOCCA=1,NUMO
            EAB1(IOCCA,IOCCB,N) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     CLEAR THE ARRAY FOR (AR|BS) VALUES
      M = 0
      DO IOCCC=1,NUMO
        DO IOCCA=1,NUMO
          DO IOCCB=1,IOCCA
            DO IOCCD=1,NUMO
              M = M+1
              CADB(M) = DCMPLX(0.0D0,0.0D0)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     ORDERED INDEX OF (ICNT,KQN,MQN) COMBINATIONS                     C
C**********************************************************************C
C
      ICOUNT = 0
C
C     LOOP OVER NUCLEAR CENTRES
      DO ICT=1,NCNT
C
C       LOOP OVER KAPPA VALUES FOR THIS NUCLEAR CENTRE
        DO KN=1,NKAP(ICT)
C
C         IMPORT KAPPA, MAXIMUM MQN
          KAPPA = KAPA(KN,ICT)
          MJMAX = 2*IABS(KAPPA)-1
C
C         LOOP OVER MQN VALUES AND RECORD INDEX
          DO MJ=1,MJMAX,2
            ICOUNT              = ICOUNT+1
            INDEX(ICT,KAPPA,MJ) = ICOUNT
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     ONE-BODY ENERGIES (INSIGNIFICANT COMPUTATIONAL COST)             C
C**********************************************************************C
C
C     RECORD TIME AT THE START OF THIS PROCESS
      CALL CPU_TIME(T1)
C
C     CALCULATE ONE-BODY MATRIX REPS
      CALL ONEEL
C
C     LOOP OVER ALL OCCUPIED ORBITAL PAIRS AND CONTRACT
      E1H = 0.0D0
      DO IOCCA=1,NUMO
        DO IOCCB=1,IOCCA
C
C         FOCK MATRIX ADDRESS FOR IOCCA AND IOCCB
          IA = MINO-1+IOCCA+NSKP
          IB = MINO-1+IOCCB+NSKP
C
C         ONE-BODY ENERGY
          EN = 0.0D0
          ET = 0.0D0
          IF(IOCCA.EQ.IOCCB) THEN
            DO J=1,NDIM
              DO I=1,NDIM
                EN = EN + DREAL(HNUC(I,J)*DCONJG(COEF(I,IA))*COEF(J,IA))
                ET = ET + DREAL(HKIN(I,J)*DCONJG(COEF(I,IA))*COEF(J,IA))
              ENDDO
            ENDDO
          ENDIF
          EAB1(IOCCA,IOCCB,1) = EN
          EAB1(IOCCA,IOCCB,2) = ET
          EAB1(IOCCA,IOCCB,3) = EN+ET
          E1H = E1H + EAB1(IOCCA,IOCCB,3)
C
        ENDDO
      ENDDO
C
C     RECORD TIME AT THE END OF THIS PROCESS
      CALL CPU_TIME(T2)
      T1EL = T1EL+T2-T1
C
C**********************************************************************C
C     LOOP OVER ATOMIC CENTRES A AND B (USE INDEX 1000)                C
C**********************************************************************C
C
C     LOOP OVER CENTRE A
      DO 1000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 1000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C**********************************************************************C
C     LOOP OVER KQN SYMMETRY TYPES A AND B (USE INDEX 2000)            C
C**********************************************************************C
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,ICNTA)
        IF(KQN(1).LT.0) THEN
          LQN(1) =-KQN(1)-1
        ELSE
          LQN(1) = KQN(1)
        ENDIF
        JQN(1) = 2*IABS(KQN(1))-1
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),ICNTA)
        DO IBAS=1,NBAS(1)
          EXL(IBAS,1) = BEXL(IBAS,LQN(1),ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KAPA(KB,ICNTB)
        IF(KQN(2).LT.0) THEN
          LQN(2) =-KQN(2)-1
        ELSE
          LQN(2) = KQN(2)
        ENDIF
        JQN(2) = 2*IABS(KQN(2))-1
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1, NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
        ENDDO
C
C**********************************************************************C
C     LOOP OVER |MQN| PROJECTIONS A AND B (INDEX 3000)                 C
C**********************************************************************C
C
C     LOOP OVER |MQN(A)| VALUES
      DO 3000 MA=1,IABS(KQN(1))
        MQN(1) = 2*MA-1
C
C     LOOP OVER |MQN(B)| VALUES
      DO 3000 MB=1,IABS(KQN(2))
        MQN(2) = 2*MB-1
C
C     INDEX ASSIGNMENT
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
C
C     SKIP CONTRIBUTIONS THAT ARISE BY PERMUTATION OF INTEGRALS
      IF(G2INT.EQ.'COULM') THEN
        IF(IQ1.LT.IQ2) GOTO 3001
      ENDIF
C
C     EQ-COEFFICIENT STARTING ADDRESSES FOR (AB) PAIR
      IABLL = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB)
      IABLS = IADILS(ICNTA,ICNTB,KA,KB,MA,MB)
      IABSL = IADISL(ICNTA,ICNTB,KA,KB,MA,MB)
      IABSS = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB)
C
C     EQ-COEFFICIENT PHASE FACTORS FOR PERMUTATION OF R-INTEGRALS
      IF(G2INT.EQ.'COULM') THEN
        PAB1 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
        PAB2 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      ELSEIF(G2INT.EQ.'BREIT') THEN
        PAB1 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
        PAB2 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      ENDIF
C
C**********************************************************************C
C     LOOP OVER COMPONENT OVERLAP LABELS A AND B (INDEX 4000)          C
C**********************************************************************C
C
C     COMPONENT LABEL FOR A AND B: TT = LL(1) or SS(4) <- COULOMB
C                                  TT = LS(2) or SL(3) <- BREIT
      DO 4000 IT1=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITN(1) = IT1
C
C       CALCULATE STARTING ADDRESS
        IF(ITN(1).EQ.1) THEN
          NADDA = 0
          NADDB = 0
        ELSEIF(ITN(1).EQ.2) THEN
          NADDA = 0
          NADDB = NSKP
        ELSEIF(ITN(1).EQ.3) THEN
          NADDA = NSKP
          NADDB = 0
        ELSEIF(ITN(1).EQ.4) THEN
          NADDA = NSKP
          NADDB = NSKP
        ENDIF
C
C       FLAG READ-IN OF E0(AB) COEFFICIENTS FOR THIS COMPONENT LABEL
        IEAB = 1
C
C**********************************************************************C
C     FOCK MATRIX STARTING ADDRESSES                                   C
C**********************************************************************C
C
C     FOCK ADDRESS FOR EACH BASIS FUNCTION (WITH SPIN PROJECTION)
      NA1 = LRGE(ICNTA,KA,2*MA-1) + NADDA
      NA2 = LRGE(ICNTA,KA,2*MA  ) + NADDA
C
      NB1 = LRGE(ICNTB,KB,2*MB-1) + NADDB
      NB2 = LRGE(ICNTB,KB,2*MB  ) + NADDB
C
C     CLEAR ARRAY FOR THE COMPLETED CONTRACTION OVER BLOCKS C AND D
      DO MDB=1,NUMO*NUMO
        DO MIJ=1,NBAS(1)*NBAS(2)
          DO IJSPIN=1,4
            DB(MIJ,MDB,IJSPIN) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     LOOP OVER ATOMIC CENTRES C AND D (USE INDEX 5000)                C
C**********************************************************************C
C
C     LOOP OVER CENTRE C
      DO 5000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE C
        XYZ(1,3) = BXYZ(1,ICNTC)
        XYZ(2,3) = BXYZ(2,ICNTC)
        XYZ(3,3) = BXYZ(3,ICNTC)
C
C     LOOP OVER CENTRE D
      DO 5000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE D
        XYZ(1,4) = BXYZ(1,ICNTD)
        XYZ(2,4) = BXYZ(2,ICNTD)
        XYZ(3,4) = BXYZ(3,ICNTD)
C
C     NUMBER OF NUCLEAR CENTRES INVOLVED IN THIS OVERLAP
      MCNT = NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
C
C**********************************************************************C
C     LOOP OVER KQN SYMMETRY TYPES C AND D (USE INDEX 6000)            C
C**********************************************************************C
C
C     LOOP OVER KQN(C) VALUES
      DO 6000 KC=1,NKAP(ICNTC)
C
C       QUANTUM NUMBERS FOR BLOCK C
        KQN(3) = KAPA(KC,ICNTC)
        IF(KQN(3).LT.0) THEN
          LQN(3) =-KQN(3)-1
        ELSE
          LQN(3) = KQN(3)
        ENDIF
        JQN(3) = 2*IABS(KQN(3))-1
C
C       BASIS EXPONENTS FOR BLOCK C
        NBAS(3) = NFNC(LQN(3),ICNTC)
        DO KBAS=1,NBAS(3)
          EXL(KBAS,3) = BEXL(KBAS,LQN(3),ICNTC)
        ENDDO
C
C     LOOP OVER KQN(D) VALUES
      DO 6000 KD=1,NKAP(ICNTD)
C
C       QUANTUM NUMBERS FOR BLOCK D
        KQN(4) = KAPA(KD,ICNTD)
        IF(KQN(4).LT.0) THEN
          LQN(4) =-KQN(4)-1
        ELSE
          LQN(4) = KQN(4)
        ENDIF
        JQN(4) = 2*IABS(KQN(4))-1
C
C       BASIS EXPONENTS FOR BLOCK D
        NBAS(4) = NFNC(LQN(4),ICNTD)
        DO LBAS=1,NBAS(4)
          EXL(LBAS,4) = BEXL(LBAS,LQN(4),ICNTD)
        ENDDO
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON KQN                           C
C**********************************************************************C
C
C     ATOM-CENTRED SELECTION RULES
      IF(MCNT.EQ.1) THEN
C
C       LQN PAIR PARITY (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQN(1)+LQN(2),2)
        IPARCD = MOD(LQN(3)+LQN(4),2)
C
C       LQN PAIR PARITY SELECTION RULE
        IF(IPARAB.NE.IPARCD) THEN
          GOTO 6001
        ENDIF
C
C       JQN TRIANGLE RULE CHECK FOR MULTIPOLE EXPANSION (ATOM-CENTRED)
        NUI = MAX0(IABS(JQN(1)-JQN(2))/2,IABS(JQN(3)-JQN(4))/2)
        NUF = MIN0(    (JQN(1)+JQN(2))/2,    (JQN(3)+JQN(4))/2)
        IF(NUI.GT.NUF) THEN
          GOTO 6001
        ENDIF
C
C       ADDITIONAL LQN SELECTION RULE PARITY ANALYSIS
        ISELK = 0
        DO NU=NUI,NUF
C
C         A AND B: LQN(1)+LQN(2)+NU EVEN OR ODD (0 IF EVEN, 1 IF ODD)
          IPARAB = MOD(LQN(1)+LQN(2)+NU,2)
C
C         C AND D: LQN(3)+LQN(4)+NU EVEN OR ODD (0 IF EVEN, 1 IF ODD)
          IPARCD = MOD(LQN(3)+LQN(4)+NU,2)
C
C         CASE 1: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH ODD (UNLESS NU=0)
          IF(IPARAB.EQ.1.AND.IPARCD.EQ.1.AND.NU.NE.0) ISELK = 1
C
C         CASE 2: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH EVEN
          IF(IPARAB.EQ.0.AND.IPARCD.EQ.0) ISELK = 1
C
        ENDDO
        IF(ISELK.EQ.0) GOTO 6001
C
      ENDIF
C
C     ADDITIONAL SELECTION RULES FOR CLOSED-SHELL ATOMS
C      IF(MCNT.EQ.1) THEN
C        ISELK = 0
C        IF(KA.EQ.KB.AND.KC.EQ.KD) ISELK = 1
C        IF(KA.EQ.KC.AND.KB.EQ.KD) ISELK = 1
C        IF(KA.EQ.KD.AND.KB.EQ.KC) ISELK = 1
C        IF(ISELK.EQ.0) GOTO 6001
C      ENDIF
C
C**********************************************************************C
C     LOOP OVER |MQN| PROJECTIONS C AND D (INDEX 7000)                 C
C**********************************************************************C
C
C     LOOP OVER |MQN(C)| VALUES
      DO 7000 MC=1,IABS(KQN(3))
        MQN(3) = 2*MC-1
C
C     LOOP OVER |MQN(D)| VALUES
      DO 7000 MD=1,IABS(KQN(4))
        MQN(4) = 2*MD-1
C
C     INDEX ASSIGNMENT
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
C     SKIP CONTRIBUTIONS THAT ARISE BY PERMUTATION OF INTEGRALS
      IF(G2INT.EQ.'COULM') THEN
        IF(IQ3.LT.IQ4) GOTO 7001
      ENDIF
C
C     EQ-COEFFICIENT STARTING ADDRESSES FOR (CD) PAIR
      ICDLL = IAD0LL(ICNTC,ICNTD,KC,KD,MC,MD)
      ICDLS = IADILS(ICNTC,ICNTD,KC,KD,MC,MD)
      ICDSL = IADISL(ICNTC,ICNTD,KC,KD,MC,MD)
      ICDSS = IAD0SS(ICNTC,ICNTD,KC,KD,MC,MD)
C
C     EQ-COEFFICIENT PHASE FACTORS FOR PERMUTATION OF R-INTEGRALS
      IF(G2INT.EQ.'COULM') THEN
        PCD1 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
        PCD2 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
      ELSEIF(G2INT.EQ.'BREIT') THEN
        PCD1 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
        PCD2 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
      ENDIF
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON MQN                           C
C**********************************************************************C
C
C     SPIN PROJECTION CONSERVED ALONG Z-AXIS FOR LINEAR MOLECULES
      IF(ISYM.EQ.1) THEN
        ISELM = 0
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) ISELM = 1
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) ISELM = 1
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) ISELM = 1
        IF(ISELM.EQ.0) GOTO 7001
      ENDIF
C
C     ATOM-CENTRED SELECTION RULES (ONLY APPLIES IF RACAH1 SWITCHED OFF)
      IF(MCNT.EQ.1) THEN
        ISELM = 0
        DO ISGN1=1,2
          DO ISGN2=1,2
            DO ISGN3=1,2
              DO ISGN4=1,2
                MMJA = MQN(1)*((-1)**ISGN1)
                MMJB = MQN(2)*((-1)**ISGN2)
                MMJC = MQN(3)*((-1)**ISGN3)
                MMJD = MQN(4)*((-1)**ISGN4)
                IF(MMJA-MMJB.EQ.MMJD-MMJC) ISELM = 1
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IF(ISELM.EQ.0) GOTO 7001
      ENDIF
C
C**********************************************************************C
C     LOOP OVER COMPONENT OVERLAP LABELS C AND D (INDEX 8000)          C
C**********************************************************************C
C
C     COMPONENT LABEL FOR C AND D: TT = LL(1) or SS(4) <- COULOMB
C                                  TT = LS(2) or SL(3) <- BREIT
      DO 8000 IT2=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITN(2) = IT2
C
C       CALCULATE STARTING ADDRESS
        IF(ITN(2).EQ.1) THEN
          NADDC = 0
          NADDD = 0
        ELSEIF(ITN(2).EQ.2) THEN
          NADDC = 0
          NADDD = NSKP
        ELSEIF(ITN(2).EQ.3) THEN
          NADDC = NSKP
          NADDD = 0
        ELSEIF(ITN(2).EQ.4) THEN
          NADDC = NSKP
          NADDD = NSKP
        ENDIF
C
C       FLAG READ-IN OF E0(CD) COEFFICIENTS FOR THIS COMPONENT LABEL
        IECD = 1
C
C**********************************************************************C
C     FOCK MATRIX STARTING ADDRESSES                                   C
C**********************************************************************C
C
C     FOCK ADDRESS FOR EACH BASIS FUNCTION (WITH SPIN PROJECTION)
      NC1 = LRGE(ICNTC,KC,2*MC-1) + NADDC
      NC2 = LRGE(ICNTC,KC,2*MC  ) + NADDC
C
      ND1 = LRGE(ICNTD,KD,2*MD-1) + NADDD
      ND2 = LRGE(ICNTD,KD,2*MD  ) + NADDD
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS IN BLOCKS A AND B (INDEX 9000)         C
C**********************************************************************C
C
      DO 9000 IBAS=1,NBAS(1)
      DO 9000 JBAS=1,NBAS(2)
C
C     LIST ADDRESS FOR THIS IBAS,JBAS COMBINATION
      MIJ = (IBAS-1)*NBAS(2)+JBAS
C
C     RESET SCREENING COUNTERS
      DO M=1,NBAS(3)*NBAS(4)
        IMAP(M) = M
        ISCR(M) = 1
      ENDDO
      MAXN = NBAS(3)*NBAS(4)
C
C     OVERRIDE GAUNT REPLACEMENT TOGGLE
      BGAUNT = .FALSE.
C
C     BATCH OF TWO-BODY INTEGRALS (IJ|KL) FOR FIXED (IJ)
      CALL CPU_TIME(T1)
      IF(G2INT.EQ.'COULM') THEN
        CALL ERI(RR,XYZ,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN)
      ELSEIF(G2INT.EQ.'BREIT') THEN
        CALL BII(RR,XYZ,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN,BGAUNT)
      ENDIF
      CALL CPU_TIME(T2)
      TERI = TERI+T2-T1
C
C     CLEAR ARRAY FOR FIRST CONTRACTION (DIRECT)
      DO MKB=1,NBAS(3)*NUMO
        DO IJKSPIN=1,8
          B1(MKB,IJKSPIN) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     CLEAR ARRAY FOR FIRST CONTRACTION (SWAP)
      DO MLB=1,NBAS(4)*NUMO
        DO IJLSPIN=1,8
          B2(MLB,IJLSPIN) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     FIRST CONTRACTION:                                               C
C     (IJ;T|KL;T') -> (IJ;T|KB;T')  AND  (IJ;T|LK;T') -> (IJ;T|LB;T')  C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     FIRST CONTRACTION (NORMAL): (IJ;T|KL;T') -> (IJ;T|KB;T')
C
C     LOOP OVER BASIS FUNCTIONS IN BLOCK C AND OCCUPIED STATES IOCCB
      DO KBAS=1,NBAS(3)
        DO IOCCB=1,NUMO
C
C         FOCK MATRIX ADDRESS FOR IOCCB
          IB = MINO-1+IOCCB+NSKP
C
C         LIST ADDRESS FOR THIS KBAS,IOCCB COMBINATION
          MKB = (KBAS-1)*NUMO + IOCCB
C
C         LOOP OVER BASIS FUNCTIONS IN BLOCK D AND CONTRACT OVER ERI
          DO LBAS=1,NBAS(4)
C
C           LIST ADDRESS FOR THIS KBAS,LBAS COMBINATION
            M = (KBAS-1)*NBAS(4) + LBAS
C
C           (--|-B) = (--|--) + (--|-+)
            B1(MKB,1) = B1(MKB,1) +      RR(M, 1)*COEF(ND1+LBAS,IB)
     &                            +      RR(M, 2)*COEF(ND2+LBAS,IB)
C           (+-|-B) = (+-|--) + (+-|-+)
            B1(MKB,2) = B1(MKB,2) +      RR(M, 9)*COEF(ND1+LBAS,IB)
     &                            +      RR(M,10)*COEF(ND2+LBAS,IB)
C           (-+|-B) = (-+|--) + (-+|-+)
            B1(MKB,3) = B1(MKB,3) +      RR(M, 5)*COEF(ND1+LBAS,IB)
     &                            +      RR(M, 6)*COEF(ND2+LBAS,IB)
C           (++|-B) = (++|--) + (++|-+)
            B1(MKB,4) = B1(MKB,4) +      RR(M,13)*COEF(ND1+LBAS,IB)
     &                            +      RR(M,14)*COEF(ND2+LBAS,IB)
C           (--|+B) = (--|+-) + (--|++)
            B1(MKB,5) = B1(MKB,5) +      RR(M, 3)*COEF(ND1+LBAS,IB)
     &                            +      RR(M, 4)*COEF(ND2+LBAS,IB)
C           (+-|+B) = (+-|+-) + (+-|++)
            B1(MKB,6) = B1(MKB,6) +      RR(M,11)*COEF(ND1+LBAS,IB)
     &                            +      RR(M,12)*COEF(ND2+LBAS,IB)
C           (-+|+B) = (-+|+-) + (-+|++)
            B1(MKB,7) = B1(MKB,7) +      RR(M, 7)*COEF(ND1+LBAS,IB)
     &                            +      RR(M, 8)*COEF(ND2+LBAS,IB)
C           (++|+B) = (++|+-) + (++|++)
            B1(MKB,8) = B1(MKB,8) +      RR(M,15)*COEF(ND1+LBAS,IB)
     &                            +      RR(M,16)*COEF(ND2+LBAS,IB)
C
          ENDDO
        ENDDO
      ENDDO
C
C     FIRST CONTRACTION (SWAP): (IJ;T|LK;T') -> (IJ;T|LB;T')
      IF(G2INT.EQ.'COULM'.AND.IQ3.EQ.IQ4) GOTO 9001
      IF(G2INT.EQ.'BREIT') GOTO 9001
C
C     LOOP OVER BASIS FUNCTIONS IN BLOCK D AND OCCUPIED STATES IOCCB
      DO LBAS=1,NBAS(4)
        DO IOCCB=1,NUMO
C
C         FOCK MATRIX ADDRESS FOR IOCCB
          IB = MINO-1+IOCCB+NSKP
C
C         LIST ADDRESS FOR THIS LBAS,IOCCB COMBINATION
          MLB = (LBAS-1)*NUMO + IOCCB
C
C         LOOP OVER BASIS FUNCTIONS IN BLOCK C AND CONTRACT OVER ERI
          DO KBAS=1,NBAS(3)
C
C           LIST ADDRESS FOR THIS KBAS,LBAS COMBINATION
            M = (KBAS-1)*NBAS(4) + LBAS
C
C           (--|+B) = PAB*{(--|B+)} = PAB*{(--|++) + ((--|-+))}
            B2(MLB,1) = B2(MLB,1) + PCD1*RR(M, 4)*COEF(NC1+KBAS,IB)
     &                            + PCD2*RR(M, 2)*COEF(NC2+KBAS,IB)
C
C           (+-|+B) = PAB*{(+-|B+)} = PAB*{(+-|++) + ((+-|-+))}
            B2(MLB,2) = B2(MLB,2) + PCD1*RR(M,12)*COEF(NC1+KBAS,IB)
     &                            + PCD2*RR(M,10)*COEF(NC2+KBAS,IB)
C
C           (-+|+B) = PAB*{(-+|B+)} = PAB*{(-+|++) + ((-+|-+))}
            B2(MLB,3) = B2(MLB,3) + PCD1*RR(M, 8)*COEF(NC1+KBAS,IB)
     &                            + PCD2*RR(M, 6)*COEF(NC2+KBAS,IB)
C
C           (++|+B) = PAB*{(++|B+)} = PAB*{(++|++) + ((++|-+))}
            B2(MLB,4) = B2(MLB,4) + PCD1*RR(M,16)*COEF(NC1+KBAS,IB)
     &                            + PCD2*RR(M,14)*COEF(NC2+KBAS,IB)
C
C           (--|-B) = PAB*{(--|B-)} = PAB*{(--|+-) + ((--|--))}
            B2(MLB,5) = B2(MLB,5) + PCD2*RR(M, 3)*COEF(NC1+KBAS,IB)
     &                            + PCD1*RR(M, 1)*COEF(NC2+KBAS,IB)
C
C           (+-|-B) = PAB*{(+-|B-)} = PAB*{(+-|+-) + ((+-|--))}
            B2(MLB,6) = B2(MLB,6) + PCD2*RR(M,11)*COEF(NC1+KBAS,IB)
     &                            + PCD1*RR(M, 9)*COEF(NC2+KBAS,IB)
C
C           (-+|-B) = PAB*{(-+|B-)} = PAB*{(-+|+-) + ((-+|--))}
            B2(MLB,7) = B2(MLB,7) + PCD2*RR(M, 7)*COEF(NC1+KBAS,IB)
     &                            + PCD1*RR(M, 5)*COEF(NC2+KBAS,IB)
C
C           (++|-B) = PAB*{(++|B-)} = PAB*{(++|+-) + ((++|--))}
            B2(MLB,8) = B2(MLB,8) + PCD2*RR(M,15)*COEF(NC1+KBAS,IB)
     &                            + PCD1*RR(M,13)*COEF(NC2+KBAS,IB)
C
          ENDDO
        ENDDO
      ENDDO
C
C     SKIP POINT FOR IQ3=IQ4
9001  CONTINUE
C
      CALL CPU_TIME(T2)
      TCN1 = TCN1+T2-T1
C
C**********************************************************************C
C     SECOND CONTRACTION:                                   ~          C
C     (IJ;T|KB;T') -> (IJ;T|DB)  AND  (IJ;T|LB;T') -> (IJ;T|DB)        C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     SECOND CONTRACTION (NORMAL): (IJ;T|KB) -> (IJ;T|DB)
C
C     LOOP OVER OCCUPIED STATES IOCCB AND IOCCD
      DO IOCCB=1,NUMO
        DO IOCCD=1,NUMO
C
C         FOCK MATRIX ADDRESS FOR IOCCD
          ID = MINO-1+IOCCD+NSKP
C
C         LIST ADDRESS FOR THIS IOCCD,IOCCB COMBINATION IN B1
          MDB = (IOCCB-1)*NUMO + IOCCD
C
C         LOOP OVER BASIS FUNCTIONS IN BLOCK C AND CONTRACT OVER B1
          DO KBAS=1,NBAS(3)
C
C           LIST ADDRESS FOR THIS KBAS,IOCCB COMBINATION
            MKB = (KBAS-1)*NUMO + IOCCB
C
C           (--|DB) = (--|-B) + (--|+B)
            DB(MIJ,MDB,1) = DB(MIJ,MDB,1)
     &                    + B1(MKB,1)*DCONJG(COEF(NC1+KBAS,ID))
     &                    + B1(MKB,5)*DCONJG(COEF(NC2+KBAS,ID))
C
C           (-+|DB) = (-+|-B) + (-+|+B)
            DB(MIJ,MDB,2) = DB(MIJ,MDB,2)
     &                    + B1(MKB,3)*DCONJG(COEF(NC1+KBAS,ID))
     &                    + B1(MKB,7)*DCONJG(COEF(NC2+KBAS,ID))
C
C           (+-|DB) = (+-|-B) + (+-|+B)
            DB(MIJ,MDB,3) = DB(MIJ,MDB,3)
     &                    + B1(MKB,2)*DCONJG(COEF(NC1+KBAS,ID))
     &                    + B1(MKB,6)*DCONJG(COEF(NC2+KBAS,ID))
C
C           (++|DB) = (++|-B) + (++|+B)
            DB(MIJ,MDB,4) = DB(MIJ,MDB,4)
     &                    + B1(MKB,4)*DCONJG(COEF(NC1+KBAS,ID))
     &                    + B1(MKB,8)*DCONJG(COEF(NC2+KBAS,ID))
C
          ENDDO
        ENDDO
      ENDDO
C
C     SECOND CONTRACTION (SWAP): (IJ;T|LB) -> (IJ;T|DB)
      IF(G2INT.EQ.'COULM'.AND.IQ3.EQ.IQ4) GOTO 9002
      IF(G2INT.EQ.'BREIT') GOTO 9002
C
C     LOOP OVER OCCUPIED STATES IOCCB AND IOCCD
      DO IOCCB=1,NUMO
        DO IOCCD=1,NUMO
C
C         FOCK MATRIX ADDRESS FOR IOCCD
          ID = MINO-1+IOCCD+NSKP
C
C         LIST ADDRESS FOR THIS IOCCD,IOCCB COMBINATION IN B1
          MDB = (IOCCB-1)*NUMO + IOCCD
C
C         LOOP OVER BASIS FUNCTIONS IN BLOCK D AND CONTRACT OVER B2
          DO LBAS=1,NBAS(4)
C
C           LIST ADDRESS FOR THIS LBAS,IOCCB COMBINATION
            MLB = (LBAS-1)*NUMO + IOCCB
C
C           (--|DB) = (--|+B) + (--|-B)
            DB(MIJ,MDB,1) = DB(MIJ,MDB,1)
     &                    + B2(MLB,1)*DCONJG(COEF(ND1+LBAS,ID))
     &                    + B2(MLB,5)*DCONJG(COEF(ND2+LBAS,ID))
C
C           (-+|DB) = (-+|+B) + (-+|-B)
            DB(MIJ,MDB,2) = DB(MIJ,MDB,2)
     &                    + B2(MLB,3)*DCONJG(COEF(ND1+LBAS,ID))
     &                    + B2(MLB,7)*DCONJG(COEF(ND2+LBAS,ID))
C
C           (+-|DB) = (+-|+B) + (+-|-B)
            DB(MIJ,MDB,3) = DB(MIJ,MDB,3)
     &                    + B2(MLB,2)*DCONJG(COEF(ND1+LBAS,ID))
     &                    + B2(MLB,6)*DCONJG(COEF(ND2+LBAS,ID))
C
C           (++|DB) = (++|+B) + (++|-B)
            DB(MIJ,MDB,4) = DB(MIJ,MDB,4)
     &                    + B2(MLB,4)*DCONJG(COEF(ND1+LBAS,ID))
     &                    + B2(MLB,8)*DCONJG(COEF(ND2+LBAS,ID))
C
          ENDDO
        ENDDO
      ENDDO
C
C     SKIP POINT FOR IQ3=IQ4
9002  CONTINUE
C
      CALL CPU_TIME(T2)
      TCN2 = TCN2+T2-T1
C
C     END LOOP OVER BASIS PAIR (IBAS,JBAS)
9000  CONTINUE
C
C     END LOOP OVER COMPONENT OVERLAP T'T'
8000  CONTINUE
C
C     SKIP POINT FOR IQ3.LT.IQ4
7001  CONTINUE
C
C     ALL CONTRIBUTIONS FROM BLOCK (C,D) NOW ACCOUNTED FOR
7000  CONTINUE
C     SKIP POINT FOR KQN SELECTION RULES
6001  CONTINUE
6000  CONTINUE
5001  CONTINUE
5000  CONTINUE
C
C**********************************************************************C
C     THIRD CONTRACTION:                                               C
C     (IJ;T|DB) -> (IA;T|DB)  AND  (JI;T|DB) -> (JA;T|DB)              C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     THIRD CONTRACTION (NORMAL): (IJ;T|DB) -> (IA;T|DB)
C
C     CLEAR ARRAY FOR THIRD CONTRACTION (NORMAL)
      DO MADB=1,NUMO*NUMO*NUMO
        DO IBAS=1,NBAS(1)
          DO ISPIN=1,2
            ADB1(IBAS,MADB,ISPIN) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C     LOOP OVER OCCUPIED STATES IOCCA
      DO IOCCA=1,NUMO
C
C       FOCK MATRIX ADDRESS FOR IOCCA
        IA = MINO-1+IOCCA+NSKP
C
C       LOOP OVER OCCUPIED STATES IOCCB AND IOCCD
        DO IOCCB=1,IOCCA
          DO IOCCD=1,NUMO
C
C           LIST ADDRESS FOR THIS IOCCD,IOCCB COMBINATION
            MDB = (IOCCB-1)*NUMO+IOCCD
C
C           LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MDB
            MADB = (IOCCA-1)*NUMO*NUMO + MDB
C
C           LOOP OVER BASIS FUNCTIONS IN A AND B, CONTRACT OVER DB
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
C
C               LIST ADDRESS FOR THIS IBAS,JBAS COMBINATION
                MIJ = (IBAS-1)*NBAS(2)+JBAS
C
C               (-A|DB) = (--|DB) + (-+|DB)
                ADB1(IBAS,MADB,1) = ADB1(IBAS,MADB,1)
     &                      +      DB(MIJ,MDB,1)*COEF(NB1+JBAS,IA)
     &                      +      DB(MIJ,MDB,2)*COEF(NB2+JBAS,IA)
C               (+A|DB) = (+-|DB) + (++|DB)
                ADB1(IBAS,MADB,2) = ADB1(IBAS,MADB,2)
     &                      +      DB(MIJ,MDB,3)*COEF(NB1+JBAS,IA)
     &                      +      DB(MIJ,MDB,4)*COEF(NB2+JBAS,IA)
C
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     THIRD CONTRACTION (SWAP): (JI;T|DB) -> (JA;T|DB)
      IF(G2INT.EQ.'COULM'.AND.IQ1.EQ.IQ2) GOTO 4001
      IF(G2INT.EQ.'BREIT') GOTO 4001
C
C     CLEAR ARRAY FOR THIRD CONTRACTION (SWAP)
      DO MADB=1,NUMO*NUMO*NUMO
        DO JBAS=1,NBAS(2)
          DO JSPIN=1,2
            ADB2(JBAS,MADB,JSPIN) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C     LOOP OVER OCCUPIED STATES IOCCA
      DO IOCCA=1,NUMO
C
C       FOCK MATRIX ADDRESS FOR IOCCA
        IA = MINO-1+IOCCA+NSKP
C
C       LOOP OVER OCCUPIED STATES IOCCB AND IOCCD
        DO IOCCB=1,IOCCA
          DO IOCCD=1,NUMO
C
C           LIST ADDRESS FOR THIS IOCCD,IOCCB COMBINATION
            MDB = (IOCCB-1)*NUMO+IOCCD
C
C           LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MDB
            MADB = (IOCCA-1)*NUMO*NUMO + MDB
C
C           LOOP OVER BASIS FUNCTIONS IN A AND B, CONTRACT OVER DB
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
C
C               LIST ADDRESS FOR THIS IBAS,JBAS COMBINATION
                MIJ = (IBAS-1)*NBAS(2)+JBAS
C
C               (+A|DB) = PCD*{(+A|DB)} = PCD*{(++|DB) + (-+|DB)}
                ADB2(JBAS,MADB,1) = ADB2(JBAS,MADB,1)
     &                      + PAB1*DB(MIJ,MDB,4)*COEF(NA1+IBAS,IA)
     &                      + PAB2*DB(MIJ,MDB,2)*COEF(NA2+IBAS,IA)
C               (-A|DB) = PCD*{(-A|DB)} = PCD*{(+-|DB) + (--|DB)}
                ADB2(JBAS,MADB,2) = ADB2(JBAS,MADB,2)
     &                      + PAB2*DB(MIJ,MDB,3)*COEF(NA1+IBAS,IA)
     &                      + PAB1*DB(MIJ,MDB,1)*COEF(NA2+IBAS,IA)
C
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     SKIP POINT FOR IQ1=IQ2
4001  CONTINUE
C
      CALL CPU_TIME(T2)
      TCN3 = TCN3+T2-T1
C
C**********************************************************************C
C     FOURTH CONTRACTION:                      ~                       C
C     (IA;T|DB) -> (CA|DB)  AND  (JA;T|DB) -> (CA|DB)                  C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     FOURTH CONTRACTION (NORMAL): (IA;T|DB) -> (CA|DB)
C
C     LOOP OVER OCCUPIED STATES IOCCA AND IOCCC
      DO IOCCA=1,NUMO
        DO IOCCC=1,NUMO
C
C         FOCK MATRIX ADDRESS FOR IOCCC
          IC = MINO-1+IOCCC+NSKP
C
C         LOOP OVER OCCUPIED STATES IOCCB AND IOCCD
          DO IOCCB=1,IOCCA
            DO IOCCD=1,NUMO
C
C             LIST ADDRESS FOR THIS IOCCD,IOCCB COMBINATION
              MDB = (IOCCB-1)*NUMO + IOCCD
C
C             LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MDB
              MADB = (IOCCA-1)*NUMO*NUMO + MDB
C
C             LIST ADDRESS FOR THIS IOCCC,IOCCA AND THE ABOVE MDB
              MCADB = (IOCCA-1)*NUMO*IOCCA*NUMO/2
     &              + (IOCCC-1)*NUMO*IOCCA + MDB
C
C             LOOP OVER BASIS FUNCTIONS IN BLOCK A, CONTRACT OVER ADB1
              DO IBAS=1,NBAS(1)
C
C               (CA|DB) = (-A|DB) + (+A|DB)
                CADB(MCADB) = CADB(MCADB)
     &                     + ADB1(IBAS,MADB,1)*DCONJG(COEF(NA1+IBAS,IC))
     &                     + ADB1(IBAS,MADB,2)*DCONJG(COEF(NA2+IBAS,IC))
C
              ENDDO
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C                                              ~
C     FOURTH CONTRACTION (SWAP): (JA;T|DB) -> (CA|DB)
      IF(G2INT.EQ.'COULM'.AND.IQ1.EQ.IQ2) GOTO 4002
      IF(G2INT.EQ.'BREIT') GOTO 4002
C
C     LOOP OVER OCCUPIED STATES IOCCA AND IOCCC
      DO IOCCA=1,NUMO
        DO IOCCC=1,NUMO
C
C         FOCK MATRIX ADDRESS FOR IOCCC
          IC = MINO-1+IOCCC+NSKP
C
C         LOOP OVER OCCUPIED STATES IOCCB AND IOCCD
          DO IOCCB=1,IOCCA
            DO IOCCD=1,NUMO
C
C             LIST ADDRESS FOR THIS IOCCD,IOCCB COMBINATION
              MDB = (IOCCB-1)*NUMO + IOCCD
C
C             LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MDB
              MADB = (IOCCA-1)*NUMO*NUMO + MDB
C
C             LIST ADDRESS FOR THIS IOCCC,IOCCA AND THE ABOVE MDB
              MCADB = (IOCCA-1)*NUMO*IOCCA*NUMO/2
     &              + (IOCCC-1)*NUMO*IOCCA + MDB
C
C             LOOP OVER BASIS FUNCTIONS IN BLOCK B, CONTRACT OVER ACB2
              DO JBAS=1,NBAS(2)
C
C               (CA|DB) = (+A|DB) + (-A|DB)
                CADB(MCADB) = CADB(MCADB)
     &                     + ADB2(JBAS,MADB,1)*DCONJG(COEF(NB1+JBAS,IC))
     &                     + ADB2(JBAS,MADB,2)*DCONJG(COEF(NB2+JBAS,IC))
C
              ENDDO
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     SKIP POINT FOR IQ1=IQ2
4002  CONTINUE
C
      CALL CPU_TIME(T2)
      TCN4 = TCN4+T2-T1
C
C     END LOOP OVER COMPONENT OVERLAP TT
4000  CONTINUE
C
C     SKIP POINT FOR IQ1.LT.IQ2
3001  CONTINUE
C
C     ALL CONTRIBUTIONS FROM BLOCK (A,B) NOW ACCOUNTED FOR
3000  CONTINUE
2000  CONTINUE
1000  CONTINUE
C
C**********************************************************************C
C     CALCULATE SECOND-ORDER PAIR CORRELATION ENERGY FROM (CA|DB)      C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     FOR EACH IOCCA,IOCCB PAIR, SUM OVER IVRTR AND IVRTS CONTRIBUTIONS
C
C     LOOP OVER OCCUPIED STATES IOCCA AND VIRTUAL STATES IVRTR
      DO IOCCC=1,NUMO
        DO IOCCA=1,NUMO
C
C         LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
          DO IOCCD=1,NUMO
            DO IOCCB=1,IOCCA
C
C             MAIN LIST ADDRESS FOR THIS IOCCA,IVRTR,IOCCB,IVRTS
              MCADB = (IOCCA-1)*NUMO*IOCCA*NUMO/2
     &              + (IOCCC-1)*NUMO*IOCCA + (IOCCB-1)*NUMO + IOCCD
C
              IF(IOCCA.EQ.IOCCC.AND.IOCCB.EQ.IOCCD) THEN
                EAB1(IOCCA,IOCCB,4) = DREAL(CADB(MCADB))
              ENDIF
C
              IF(IOCCB.EQ.IOCCC.AND.IOCCA.EQ.IOCCD) THEN
                EAB1(IOCCA,IOCCB,5) =-DREAL(CADB(MCADB))
              ENDIF
C
C             ADD TO DIRECT AND EXCHANGE BINS
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     FILL IN THE OTHER HALF OF THE ARRAY AND CALCULATE TOTALS
      E1D = 0.0D0
      E1X = 0.0D0
      E1S = 0.0D0
      DO IOCCA=1,NUMO
        DO IOCCB=1,IOCCA
C
C         INTERMEDIATE VALUES
          EAB1DIR = EAB1(IOCCA,IOCCB,4)
          EAB1XCH = EAB1(IOCCA,IOCCB,5)
          EAB1SUM = EAB1DIR + EAB1XCH
C
C         PUT THESE INTO EAB1 AND ADD CONTRIBUTION TO E1
          EAB1(IOCCA,IOCCB,6) = EAB1SUM
          IF(IOCCA.NE.IOCCB) THEN
            EAB1(IOCCB,IOCCA,4) = EAB1DIR
            EAB1(IOCCB,IOCCA,5) = EAB1XCH
            EAB1(IOCCB,IOCCA,6) = EAB1SUM
            E1D = E1D +       EAB1DIR
            E1X = E1X +       EAB1XCH
            E1S = E1S +       EAB1SUM
          ELSE
            E1D = E1D + 0.5D0*EAB1DIR
            E1X = E1X + 0.5D0*EAB1XCH
            E1S = E1S + 0.5D0*EAB1SUM
          ENDIF
        ENDDO
      ENDDO
C
C     WRITE RESULTS OF EAB ENERGIES TO AN EXTERNAL FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_MBPT1.dat',STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO IOCCA=1,NUMO
        DO IOCCB=1,NUMO
          WRITE(8, *) (EAB1(IOCCA,IOCCB,N),N=1,6)
        ENDDO
      ENDDO
      CLOSE(UNIT=8)
C
C**********************************************************************C
C     CALCULATE SECOND-ORDER SINGLE ORBITAL ENERGY                     C
C**********************************************************************C
C
C     FOR EACH IOCCA, SUM OVER THE IOCCB CONTRIBUTIONS
      DO IOCCA=1,NUMO
        DO N=1,6
          EA1(IOCCA,N) = 0.0D0
          DO IOCCB=1,NUMO
            EA1(IOCCA,N) = EA1(IOCCA,N) + EAB1(IOCCA,IOCCB,N)
          ENDDO
        ENDDO
      ENDDO
C
      CALL CPU_TIME(T2)
      TSUM = TSUM+T2-T1
C
C**********************************************************************C
C     TERMINAL OUTPUT SUMMARY                                          C
C**********************************************************************C
C
C     MBPT1 PAIRWISE SUMMARY
20    FORMAT(1X,A,9X,A,9X,A,9X,A,10X,A)
21    FORMAT(' (',I2,',',I2,')',3X,F13.7,5X,F11.7,5X,F11.7,4X,F13.7)
221   FORMAT(' (',I2,',',I2,')',3(F16.10,3X))
121   FORMAT(' (',A,',',A,')',3X,F16.10,5X,F16.10,5X,F16.10)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',25),'MBPT1 pairwise summary'
      WRITE(7, *) REPEAT(' ',25),'MBPT1 pairwise summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) '( a, b)','E1H(ab)','E1J(ab)','E1K(ab)','E1G(ab)'
      WRITE(7,20) '( a, b)','E1H(ab)','E1J(ab)','E1K(ab)','E1G(ab)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO IOCCA=1,NUMO
        IANUM = IOCCA+MINO-1
        DO IOCCB=1,IOCCA
          IBNUM = IOCCB+MINO-1
          WRITE(6,21) IANUM,IBNUM,(EAB1(IOCCA,IOCCB,N),N=3,6)
          WRITE(7,21) IANUM,IBNUM,(EAB1(IOCCA,IOCCB,N),N=3,6)
        ENDDO
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     MBPT1 SINGLE-PARTICLE SUMMARY
30    FORMAT(1X,A,10X,A,10X,A,10X,A,10X,A)
31    FORMAT('  ',I2,'    ',3X,F13.7,5X,F11.7,5X,F11.7,4X,F13.7)

      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',21),'MBPT1 single particle summary'
      WRITE(7, *) REPEAT(' ',21),'MBPT1 single particle summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,30) '  a    ','E1H(a)','E1J(a)','E1K(a)',' E1G(a)'
      WRITE(7,30) '  a    ','E1H(a)','E1J(a)','E1K(a)',' E1G(a)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO IOCCA=1,NUMO
        IANUM = IOCCA+MINO-1
        WRITE(6,31) IANUM,(EA1(IOCCA,N),N=3,6)
        WRITE(7,31) IANUM,(EA1(IOCCA,N),N=3,6)
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     MBPT1 TOTAL FIRST-ORDER INTERACTION
32    FORMAT(' total  ',3X,F13.7,5X,F11.7,5X,F11.7,4X,F13.7)

      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',19),'MBPT1 first order molecular energy'
      WRITE(7, *) REPEAT(' ',19),'MBPT1 first order molecular energy'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,32) E1H,E1D,E1X,E1S
      WRITE(7,32) E1H,E1D,E1X,E1S
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     MBPT1 LABOUR ANALYSIS
      CALL CPU_TIME(TFIN)
      TTOT = TFIN - TBEG
      TOTH = TTOT - (T1EL + TERI + TCN1 + TCN2 + TSUM)

40    FORMAT(1X,A,15X,A)
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',26),'MBPT1 labour analysis'
      WRITE(7, *) REPEAT(' ',26),'MBPT1 labour analysis'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
      WRITE(6,40) 'One-body terms - EH(A,B)                 ',HMS(T1EL)
      WRITE(7,40) 'One-body terms - EH(A,B)                 ',HMS(T1EL)
      WRITE(6,40) 'ERI construction - (IJ|KL)               ',HMS(TERI)
      WRITE(7,40) 'ERI construction - (IJ|KL)               ',HMS(TERI)
      WRITE(6,40) '1st contraction  - (IJ|KB)               ',HMS(TCN1)
      WRITE(7,40) '1st contraction  - (IJ|KB)               ',HMS(TCN1)
      WRITE(6,40) '2nd contraction  - (IJ|BB) and (IJ|AB)   ',HMS(TCN2)
      WRITE(7,40) '2nd contraction  - (IJ|BB) and (IJ|AB)   ',HMS(TCN2)
      WRITE(6,40) '3rd contraction  - (IA|BB) and (IA|AB)   ',HMS(TCN3)
      WRITE(7,40) '3rd contraction  - (IA|BB) and (IA|AB)   ',HMS(TCN3)
      WRITE(6,40) '4th contraction  - (AA|BB) and (BA|AB)   ',HMS(TCN4)
      WRITE(7,40) '4th contraction  - (AA|BB) and (BA|AB)   ',HMS(TCN4)
      WRITE(6,40) 'Other                                    ',HMS(TOTH)
      WRITE(7,40) 'Other                                    ',HMS(TOTH)
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,40) 'Total MBPT1 time                         ',HMS(TTOT)
      WRITE(7,40) 'Total MBPT1 time                         ',HMS(TTOT)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
      RETURN
      END

