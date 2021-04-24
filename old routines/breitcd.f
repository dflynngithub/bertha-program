      SUBROUTINE BREITCD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      BBBBBBB  RRRRRRR  EEEEEEEE IIII TTTTTTTT CCCCCC  DDDDDDD        C
C      BB    BB RR    RR EE        II     TT   CC    CC DD    DD       C
C      BB    BB RR    RR EE        II     TT   CC       DD    DD       C
C      BBBBBBB  RR    RR EEEEEE    II     TT   CC       DD    DD       C
C      BB    BB RRRRRRR  EE        II     TT   CC       DD    DD       C
C      BB    BB RR    RR EE        II     TT   CC    CC DD    DD       C
C      BBBBBBB  RR    RR EEEEEEEE IIII    TT    CCCCCC  DDDDDDD        C
C                                                                      C
C -------------------------------------------------------------------- C
C  BREITCD GENERATES ALL MANY-CENTRE BREIT INTERACTION INTEGRALS IN    C
C  BATCHES AND ADDS THEM TO THE SCF CLOSED/OPEN-SHELL BREIT MATRIX.    C
C  CALCULATIONS ARE MADE WITH A RELATIVISTIC MCMURCHIE-DAVIDSON SCHEME.C
C -------------------------------------------------------------------- C
C  THIS IS A SPECIAL VERSION OF THE USUAL BREIT ROUTINE, DESIGNED TO   C
C  RECYCLE THE INTERMEDIATE INTEGRALS G(  |CD) BY REARRANGING THE LOOP C
C  ORDER WITHIN A GIVEN SET OF LQN'S. LOOP OVER (AB) BASIS FUNCTIONS   C
C  FIRST, AND THEN USE THE SAME INTERMEDIATE INTEGRALS FOR ALL (CD)    C
C  BASIS FUNCTIONS. THE SWAP OF CONTRACTION ORDER IS DONE IN BIICD.    C
C -------------------------------------------------------------------- C
C  TODO: THIS ROUTINE COULD BENEFIT FROM PARALLELISATION -- OPENMP.    C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL GAUNT
C
      CHARACTER*80 TITLE
C
      DIMENSION EXL(MBS,4),XYZ(3,4)
      DIMENSION ICNT(4),LQN(4),KQN(4),JQN(4),MQN(4),NBAS(4)
      DIMENSION INDEX(MCT,-(MEL+1):MEL,MKP),ITN(2)
      DIMENSION MAPTTTT(4,4)
C
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VANM(MDM,MDM),
     &           VSLF(MDM,MDM),VUEH(MDM,MDM),VWKR(MDM,MDM),
     &           VKSB(MDM,MDM),QDIR(MDM,MDM),QXCH(MDM,MDM),
     &           WDIR(MDM,MDM),WXCH(MDM,MDM),CPLE(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EILS/EILSFL(MFL,12),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EISL/EISLFL(MFL,12),IADISL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/I2EL/PAB1,PAB2,PCD1,PCD2,NA1,NB1,NC1,ND1,NA2,NB2,NC2,ND2,
     &            IBAS,JBAS,MCNT,NADDAB,NADDCD,NBAS,IQL,IQR
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS,IABSL,ICDSL
      COMMON/IRCM/IEAB,IECD,NCD,IGAB,IRIJ(MBS,MBS)
      COMMON/ISCR/IMTX(MB2,11),ISCR(MB2),IMAP(MB2),IBCH,ITOG,MAXN
      COMMON/LSHF/SHLEV(4),SHLV,ILEV
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VANM,VSLF,
     &            VUEH,VWKR,VKSB,QDIR,QXCH,WDIR,WXCH,CPLE
      COMMON/SHLL/ACFF,BCFF,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
      COMMON/SWRZ/GDSC(MDM,MDM),BDSC(MDM,MDM)
      COMMON/T2EL/F2ES(5,9),T2ES(5,9),N2EB(5,9),N2EI(5,9),N2ES(5,9)
      COMMON/TSCF/TC1A,TC1I,TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRW,TCC1,
     &            TCC2,TCMC,TB1A,TB1I,TB1B,TB1R,TB1F,TB1M,TB1T,TBEC,
     &            TBRM,TBRW,TBC1,TBC2,TBMC,TSMX,TUMX,THMX,TAMX,TC1T,
     &            TC2T,TCVT,TB2T,TACC,TEIG,TSCR,TTOT,TC2S,TB2S
C
C     TWO-ELECTRON COMPONENT OVERLAP ADDRESSES
      DATA MAPTTTT/1,0,0,2,0,5,6,0,0,7,8,0,3,0,0,4/
C
C     INTEGRAL SCREENING SENSITIVITY PARAMETER
      DATA SENS/1.0D-12/
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
C     COMPONENT OVERLAP LABELS TO LOOP OVER (INTEGRAL SYMMETRY MATTERS)
      IF(INTSYM) THEN
        ITSTRT = 2
        ITSTOP = 2
        ITSKIP = 1
      ELSE
        ITSTRT = 2
        ITSTOP = 3
        ITSKIP = 1
      ENDIF
C
C     INITIALISE STORAGE MATRICES
      DO I=1,NDIM
        DO J=1,NDIM
          BDIR(I,J) = DCMPLX(0.0D0,0.0D0)
          BXCH(I,J) = DCMPLX(0.0D0,0.0D0)
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
            LENIQ               = ICOUNT
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     LOOP OVER ALL NUCLEAR CENTRES (USE INDEX 1000)                   C
C**********************************************************************C
C
C     LOOP OVER CENTRE A
      DO 1000 ICNTA=1,NCNT
        ICNT(1) = ICNTA
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 1000 ICNTB=1,NCNT
        ICNT(2) = ICNTB
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C     LOOP OVER CENTRE C
      DO 1000 ICNTC=1,NCNT
        ICNT(3) = ICNTC
C
C       CARTESIAN COORDINATES OF CENTRE C
        XYZ(1,3) = BXYZ(1,ICNTC)
        XYZ(2,3) = BXYZ(2,ICNTC)
        XYZ(3,3) = BXYZ(3,ICNTC)
C
C     LOOP OVER CENTRE D
      DO 1000 ICNTD=1,NCNT
        ICNT(4) = ICNTD
C
C       CARTESIAN COORDINATES OF CENTRE D
        XYZ(1,4) = BXYZ(1,ICNTD)
        XYZ(2,4) = BXYZ(2,ICNTD)
        XYZ(3,4) = BXYZ(3,ICNTD)
C
C     NUMBER OF NUCLEAR CENTRES INVOLVED IN THIS OVERLAP
      MCNT = NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
C
C     SKIP ONE-CENTRE CONTRIBUTIONS (DEFER TO RACAH ALGEBRA ROUTINE)
      IF(MCNT.EQ.1.AND.RACAH1) THEN
        GOTO 1001
      ENDIF
C
C     SKIP MULTI-CENTRE CONTRIBUTIONS IN STAGE 1
      IF(MCNT.NE.1.AND.ILEV.LT.3) GOTO 1001
C
C     SKIP INTEGRALS UNLESS THEY COME FROM TWO 1-CENTRE DENSITIES
      IF(VSDNSB) THEN
        IF(ICNTA.NE.ICNTB.OR.ICNTC.NE.ICNTD) GOTO 1001
      ENDIF
C
C     REPLACE BREIT INTERACTION WITH GAUNT ONLY
      IF(GAUNT2) THEN
        GAUNT = .TRUE.
      ELSE
        GAUNT = .FALSE.
      ENDIF
C
C**********************************************************************C
C     LOOP OVER ALL LQN ORBITAL TYPES (USE INDEX 2000)                 C
C**********************************************************************C
C
C     LOOP OVER LQN(A) VALUES
      DO 2000 LA=0,(NKAP(ICNTA)-1)/2
C
C       QUANTUM NUMBERS FOR BLOCK A
        LQN(1) = LA
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),ICNTA)
        DO IBAS=1,NBAS(1)
          EXL(IBAS,1) = BEXL(IBAS,LQN(1),ICNTA)
        ENDDO
C
C     LOOP OVER LQN(B) VALUES
      DO 2000 LB=0,(NKAP(ICNTB)-1)/2
C
C       QUANTUM NUMBERS FOR BLOCK B
        LQN(2) = LB
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
        ENDDO
C
C     LOOP OVER LQN(C) VALUES
      DO 2000 LC=0,(NKAP(ICNTC)-1)/2
C
C       QUANTUM NUMBERS FOR BLOCK C
        LQN(3) = LC
C
C       BASIS EXPONENTS FOR BLOCK C
        NBAS(3) = NFNC(LQN(3),ICNTC)
        DO KBAS=1,NBAS(3)
          EXL(KBAS,3) = BEXL(KBAS,LQN(3),ICNTC)
        ENDDO
C
C     LOOP OVER LQN(D) VALUES
      DO 2000 LD=0,(NKAP(ICNTD)-1)/2
C
C       QUANTUM NUMBERS FOR BLOCK D
        LQN(4) = LD
C
C       BASIS EXPONENTS FOR BLOCK D
        NBAS(4) = NFNC(LQN(4),ICNTD)
        DO LBAS=1,NBAS(4)
          EXL(LBAS,4) = BEXL(LBAS,LQN(4),ICNTD)
        ENDDO
C
C     THIS UNIQUELY DEFINES A FULL SET OF RC(AB|CD) INTEGRALS -- RESET
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          IRIJ(IBAS,JBAS) = 1
        ENDDO
      ENDDO
C
C     EVALUATE SPARSITY LISTS FOR UPCOMING EQ(CD) BATCH
      CALL EILSCD0(ICNT,LQN,NBAS)
      IF(.NOT.INTSYM) CALL EISLCD0(ICNT,LQN,NBAS)
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON KQN                           C
C     NOTE THAT IF 'RACAH' IS SWITCHED ON, THIS WON'T BE CALLED.       C
C**********************************************************************C
C
C     ATOM-CENTRED SELECTION RULES (ONLY APPLIES IF RACAH1 SWITCHED OFF)
      IF(MCNT.EQ.1) THEN
C
C       LQN PAIR PARITY (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQN(1)+LQN(2),2)
        IPARCD = MOD(LQN(3)+LQN(4),2)
C
C       LQN PAIR PARITY SELECTION RULE
        IF(IPARAB.NE.IPARCD) THEN
          GOTO 2001
        ENDIF
C
      ENDIF
C
C**********************************************************************C
C     LOOP OVER ALL KQN(AB) SYMMETRY TYPES (USE INDEX 3000)            C
C**********************************************************************C
C
C     LOOP OVER KQN(A) VALUES
      DO 3000 NA=KRONECK(LA,0),1
        KA = 2*LA+NA
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,ICNTA)
        JQN(1) = 2*IABS(KQN(1))-1
C
C     LOOP OVER KQN(B) VALUES
      DO 3000 NB=KRONECK(LB,0),1
        KB = 2*LB+NB
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KAPA(KB,ICNTB)
        JQN(2) = 2*IABS(KQN(2))-1
C
C**********************************************************************C
C     LOOP OVER ALL |MQN|(AB) SYMMETRY TYPES (USE INDEX 4000)          C
C**********************************************************************C
C
C     LOOP OVER |MQN(A)| VALUES
      DO 4000 MA=1,IABS(KQN(1))
        MQN(1) = 2*MA-1
C
C     LOOP OVER |MQN(B)| VALUES
      DO 4000 MB=1,IABS(KQN(2))
        MQN(2) = 2*MB-1
C
C     EQ-COEFFICIENT STARTING ADDRESSES FOR (AB) PAIR
      IEAB  = 1
      IABLS = IADILS(ICNTA,ICNTB,KA,KB,MA,MB)
      IABSL = IADISL(ICNTA,ICNTB,KA,KB,MA,MB)
C
C**********************************************************************C
C     LOOP OVER ALL (AB) COMPONENT OVERLAP OPTIONS T1 (USE INDEX 5000) C
C**********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR A AND B: TT = LS (2) or SL (3)
      DO 5000 IT1=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR BII ROUTINE LATER
        ITN(1) = IT1
C
C       CALCULATE STARTING ADDRESS
        IF(IT1.EQ.2) THEN
          NADDA = 0
          NADDB = NSKP
        ELSEIF(IT1.EQ.3) THEN
          NADDA = NSKP
          NADDB = 0
        ENDIF
C
C       FOCK ADDRESS FOR EACH BASIS FUNCTION (WITH SPIN PROJECTION)
        NA1 = LRGE(ICNTA,KA,2*MA-1) + NADDA
        NA2 = LRGE(ICNTA,KA,2*MA  ) + NADDA
        NB1 = LRGE(ICNTB,KB,2*MB-1) + NADDB
        NB2 = LRGE(ICNTB,KB,2*MB  ) + NADDB
C
C**********************************************************************C
C     LOOP OVER ALL (AB) BASIS FUNCTIONS (IBAS,JBAS) (USE INDEX 6000)  C
C**********************************************************************C
C
C     RECORD TIME AT START OF BATCH
      CALL SYSTEM_CLOCK(ICL1,RATE)
C
      DO 6000 IBAS=1,NBAS(1)
      DO 6000 JBAS=1,NBAS(2)
C
C     PERFORM FIRST CONTRACTION AGAIN
      IGAB = 1
C
C**********************************************************************C
C     LOOP OVER ALL KQN(CD) SYMMETRY TYPES (USE INDEX 7000)            C
C**********************************************************************C
C
C     LOOP OVER KQN(C) VALUES
      DO 7000 NC=KRONECK(LC,0),1
        KC = 2*LC+NC
C
C       QUANTUM NUMBERS FOR BLOCK C
        KQN(3) = KAPA(KC,ICNTC)
        JQN(3) = 2*IABS(KQN(3))-1
C
C     LOOP OVER KQN(D) VALUES
      DO 7000 ND=KRONECK(LD,0),1
        KD = 2*LD+ND
C
C       QUANTUM NUMBERS FOR BLOCK D
        KQN(4) = KAPA(KD,ICNTD)
        JQN(4) = 2*IABS(KQN(4))-1
C
C     UNIQUE INDEX FOR THIS PAIR OF KAPPAS
      NCD = 2*NC+ND+1
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON KQN                           C
C     NOTE THAT IF 'RACAH' IS SWITCHED ON, THIS WON'T BE CALLED.       C
C**********************************************************************C
C
C     ATOM-CENTRED SELECTION RULES (ONLY APPLIES IF RACAH1 SWITCHED OFF)
      IF(MCNT.EQ.1) THEN
C
C       JQN TRIANGLE RULE CHECK FOR MULTIPOLE EXPANSION (ATOM-CENTRED)
        NUI = MAX0(IABS(JQN(1)-JQN(2))/2,IABS(JQN(3)-JQN(4))/2)
        NUF = MIN0(    (JQN(1)+JQN(2))/2,    (JQN(3)+JQN(4))/2)
        IF(NUI.GT.NUF) THEN
          GOTO 7001
        ENDIF
C
C       ADDITIONAL LQN SELECTION RULE PARITY ANALYSIS
        ISELK = 0
        DO NU=NUI,NUF
C
C         A AND B: LQN(1)+LQN(2)+NU EVEN OR ODD (0 IF EVEN, 1 IF ODD)
          IPARAB = MOD(LQN(1)+LQN(2)+NU,2)
          IPARCD = MOD(LQN(3)+LQN(4)+NU,2)
C
C         CASE 1: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH ODD (UNLESS NU=0)
          IF(IPARAB.EQ.1.AND.IPARCD.EQ.1.AND.NU.NE.0) ISELK = 1
C
C         CASE 2: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH EVEN
          IF(IPARAB.EQ.0.AND.IPARCD.EQ.0) ISELK = 1
C
        ENDDO
        IF(ISELK.EQ.0) GOTO 7001
C
      ENDIF
C
C**********************************************************************C
C     LOOP OVER ALL |MQN|(CD) SYMMETRY TYPES (USE INDEX 8000)          C
C**********************************************************************C
C
C     LOOP OVER |MQN(C)| VALUES
      DO 8000 MC=1,IABS(KQN(3))
        MQN(3) = 2*MC-1
C
C     LOOP OVER |MQN(D)| VALUES
      DO 8000 MD=1,IABS(KQN(4))
        MQN(4) = 2*MD-1
C
C     EQ-COEFFICIENT STARTING ADDRESSES FOR (CD) PAIR
      IECD  = 1
      ICDLS = IADILS(ICNTC,ICNTD,KC,KD,MC,MD)
      ICDSL = IADISL(ICNTC,ICNTD,KC,KD,MC,MD)
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON MQN                           C
C**********************************************************************C
C
C     SPIN PROJECTION CONSERVED ALONG Z-AXIS FOR LINEAR MOLECULES
      IF(ISYM.EQ.1.OR.ISYM.EQ.2) THEN
        ISELM = 0
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4).AND.NOPN.EQ.0) ISELM=1
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) ISELM=1
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) ISELM=1
        IF(ISELM.EQ.0) GOTO 8001
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
        IF(ISELM.EQ.0) GOTO 8001
      ENDIF
C
C**********************************************************************C
C     IDENTIFICATION OF BII SYMMETRIES AVAILABLE TO THIS BLOCK         C
C**********************************************************************C
C
C     STARTING INDEX VALUES
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
C     COMBINED BLOCK INDEX IN A TWO-FUNCTION LIST
      IQL = (IQ1*(IQ1-1))/2 + IQ2
      IQR = (IQ3*(IQ3-1))/2 + IQ4
C     IQL = (IQ1-1)*LENIQ + IQ2
C     IQR = (IQ3-1)*LENIQ + IQ4
C
C     SKIP CONTRIBUTIONS THAT ARISE BY PERMUTATION OF INTEGRALS
      IF(INTSYM.AND.IQL.LT.IQR) GOTO 8001
C
C     EQ-COEFFICIENT PHASE FACTORS FOR PERMUTATION OF R-INTEGRALS
C     NB! OPPOSITE PHASE AS IN THE LL/SS CASE SEEN IN COULOMB
      PAB1 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PAB2 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      PCD1 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PCD2 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
C**********************************************************************C
C     LOOP OVER ALL (CD) COMPONENT OVERLAP OPTIONS T2 (USE INDEX 9000) C        C
C**********************************************************************C
C                                              _
C     LOOP OVER COMPONENT LABEL FOR C AND D: T'T' = LS (2) or SL (3)
      DO 9000 IT2=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITN(2) = IT2
C
C       CALCULATE STARTING ADDRESS
        IF(IT2.EQ.2) THEN
          NADDC = 0
          NADDD = NSKP
        ELSEIF(IT2.EQ.3) THEN
          NADDC = NSKP
          NADDD = 0
        ENDIF
C
C       FOCK ADDRESS FOR EACH BASIS FUNCTION (WITH SPIN PROJECTION)
        NC1 = LRGE(ICNTC,KC,2*MC-1) + NADDC
        NC2 = LRGE(ICNTC,KC,2*MC  ) + NADDC
        ND1 = LRGE(ICNTD,KD,2*MD-1) + NADDD
        ND2 = LRGE(ICNTD,KD,2*MD  ) + NADDD
C
C       FLAG READ-IN OF E0(CD) COEFFICIENTS FOR THIS COMPONENT LABEL
        IECD = 1
C
C       COMPONENT OVERLAP INDEX {(LS|LS)=5,(LS|SL)=6,(SL|LS)=7,(SL|SL)=8}
        ITT = MAPTTTT(IT1,IT2)
C
C       UPDATE COUNTER FOR NUMBER OF BLOCKS
        N2EB(MCNT,ITT) = N2EB(MCNT,ITT)+1
C
C**********************************************************************C
C     CONSTRUCTION OF WDIR/WMAT CONTRIBUTIONS FOR THIS BATCH           C
C**********************************************************************C
C
C     SCHWARZ SCREENING (ECONOMIC ONLY WHEN SCREENING FRACTION BIG)
      IF(SCHWRZ) THEN
        ITOG = 1
      ELSE
        ITOG = 0
      ENDIF
C
C     CALL THE SCREENING ROUTINE
      CALL SCHWARZ(BDSC,SENS,TB2S)
C
C     UPDATE COUNTER FOR NUMBER OF INTEGRALS AND SCREENED INTEGRALS
      N2EI(MCNT,ITT) = N2EI(MCNT,ITT) + NBAS(3)*NBAS(4)
      N2ES(MCNT,ITT) = N2ES(MCNT,ITT) + NBAS(3)*NBAS(4)-MAXN
C
C     CONDITIONAL TO SKIP THIS BATCH
      IF(IBCH.EQ.1) THEN
C
C       GENERATE A BATCH OF BREIT INTERACTION INTEGRALS
C        IF(RCFILE) THEN
          CALL BIICD(RR,XYZ,ICNT,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN,GAUNT)
C        ELSE
C          CALL BIIF(RR,XYZ,ICNT,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN,GAUNT)
C        ENDIF
C
C       MULTIPLY BY DENSITY ELEMENTS AND ADD TO WDIR/WMAT
        CALL BRTMAT(RR,TBMC)
C
      ENDIF
C
C     RECORD TIME AT END OF BATCH
      CALL SYSTEM_CLOCK(ICL2)
      T2ES(MCNT,ITT) = T2ES(MCNT,ITT) + DFLOAT(ICL2-ICL1)/RATE
C
9000  CONTINUE
8001  CONTINUE
8000  CONTINUE
7001  CONTINUE
7000  CONTINUE
6000  CONTINUE
5000  CONTINUE
4000  CONTINUE
3000  CONTINUE
2001  CONTINUE
2000  CONTINUE
1001  CONTINUE
1000  CONTINUE
C
C**********************************************************************C
C     COMPLETE CONSTRUCTION OF BREIT MATRIX BY MATRIX CONJUGATION.     C
C**********************************************************************C
C
      DO I=1,NDIM-NSKP
        DO J=NSKP+1,NDIM
          BDIR(J,I) = DCONJG(BDIR(I,J))
          BXCH(J,I) = DCONJG(BXCH(I,J))
          WDIR(J,I) = DCONJG(WDIR(I,J))
          WXCH(J,I) = DCONJG(WXCH(I,J))
        ENDDO
      ENDDO
C
C     MULTIPLY OPEN MATRIX BY ANGULAR COEFFICIENTS (LIFTED FROM COULOMB)
      IF(NOPN.GT.0) THEN
        DO J=1,NDIM
          DO I=1,NDIM
            WDIR(I,J) = ACFF*WDIR(I,J)
            WXCH(I,J) = BCFF*WXCH(I,J)
          ENDDO
        ENDDO
      ENDIF
C
c      TITLE = 'BREIT-BACKWARD'
c      CALL ZGNUMAP(BXCH,TITLE,NDIM)
C
      RETURN
      END

