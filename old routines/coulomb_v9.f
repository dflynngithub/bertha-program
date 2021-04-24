      SUBROUTINE COULOMB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    CCCCCC   OOOOOO  UU    UU LL      OOOOOO  MM       MM BBBBBBB     C
C   CC    CC OO    OO UU    UU LL     OO    OO MMM     MMM BB    BB    C
C   CC       OO    OO UU    UU LL     OO    OO MMMM   MMMM BB    BB    C
C   CC       OO    OO UU    UU LL     OO    OO MM MM MM MM BBBBBBB     C
C   CC       OO    OO UU    UU LL     OO    OO MM  MMM  MM BB    BB    C
C   CC    CC OO    OO UU    UU LL     OO    OO MM   M   MM BB    BB    C
C    CCCCCC   OOOOOO   UUUUUU  LLLLLLL OOOOOO  MM       MM BBBBBBB     C
C                                                                      C
C -------------------------------------------------------------------- C
C  COULOMB GENERATES ALL MANY-CENTER ELECTRON REPULSION INTEGRALS IN   C
C  BATCHES AND ADDS THEM TO THE SCF CLOSED/OPEN-SHELL COULOMB MATRIX.  C
C  CALCULATIONS ARE MADE WITH A RELATIVISTIC MCMURCHIE-DAVIDSON SCHEME.C
C -------------------------------------------------------------------- C
C  DFNOTE: THIS ROUTINE COULD BENEFIT FROM PARALLELISATION -- OPENMP.  C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=15,MKP=9,MFL=7000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      CHARACTER*4 HMLTN
      CHARACTER*8 SHAPE
C
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 RR(MB2,16)
C
      DIMENSION GDSC(MDM,MDM),BDSC(MDM,MDM)
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4),LQN(4),ITN(2)
      DIMENSION ISCF(11,6),IFLG(11)
      DIMENSION INDEX(MCT,-MKP:MKP,2*(MKP+1)*MKP),ICNT(4)
C
      COMMON/BLOC/PAB1,PAB2,PCD1,PCD2,NA1,NB1,NC1,ND1,NA2,NB2,NC2,ND2,
     &            IBAS,JBAS,ILIN,NADDAB,NADDCD
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/GEOM/SHAPE
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS
      COMMON/MAKE/IEAB,IECD,IRIJ(MBS,MBS)
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,IEQS,IERC,IPAR,ICOR,ILEV
      COMMON/QNMS/LABICN(MDM),LABKQN(MDM),LABMQN(MDM)
      COMMON/SCRN/F2ES(5,7),T2ES(5,7),N2EB(5,7),N2ET(5,7),N2ES(5,7)
      COMMON/SHLL/ACFF,BCFF,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN,NOELEC
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/SWRZ/GDSC,BDSC
      COMMON/TSCF/TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRR,TCC1,TCC2,TCMC,
     &            TB1B,TB1R,TB1F,TB1M,TBEC,TBRM,TBRR,TBC1,TBC2,TBMC,
     &            THMX,TC1T,TC2T,TB1T,TB2T,TEIG,TSCR,TTOT,
     &            TC1S,TC2S,TB1S,TB2S
C
C     ISCF TELLS WHICH INTEGRALS TO INCLUDE BASED ON OVERLAP COMBINATION
      DATA ISCF/1,1,1,1,1,1,1,1,0,0,0,
     &          1,1,0,0,1,1,1,0,0,0,0,
     &          1,0,1,1,1,0,1,0,0,0,0,
     &          1,1,1,0,1,1,0,0,0,0,0,
     &          1,0,1,0,1,0,0,0,0,0,0,
     &          1,0,0,0,1,0,0,0,0,0,0/
C
C     LINK THE OPENMP ROUTINE LIBRARY
      INCLUDE 'omp_lib.h'
C
C     SAVED BATCHES OF R(AB|CD) INTEGRALS IN ERI
      IF(IEQS.EQ.0) THEN
        IERC = 0
      ELSE
        IERC = 1
      ENDIF
C
C     COMPONENT OVERLAP LABELS TO LOOP OVER
      IF(HMLTN.EQ.'BARE') THEN
        RETURN
      ELSEIF(HMLTN.EQ.'NORL') THEN
        ITSTRT = 1
        ITSTOP = 1
        ITSKIP = 1
      ELSE
        ITSTRT = 1
        ITSTOP = 4
        ITSKIP = 3
      ENDIF
C
C     LINEAR MOLECULE SHORTCUT OPTION
      IF(SHAPE.EQ.'DIATOMIC'.OR.SHAPE.EQ.'LINEAR') THEN
        ILIN = 1
      ELSE
        ILIN = 0
      ENDIF
C
C     INITIALISE STORAGE MATRICES
      DO I=1,NDIM
        DO J=1,NDIM
          GDIR(I,J) = DCMPLX(0.0D0,0.0D0)
          GXCH(I,J) = DCMPLX(0.0D0,0.0D0)
          QDIR(I,J) = DCMPLX(0.0D0,0.0D0)
          QXCH(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     ORDERED INDEX OF (ICNT,KQN,MQN) COMBINATIONS                     C
C**********************************************************************C
C
      ICOUNT = 0
C
C     LOOP OVER NUCLEAR CENTERS
      DO ICT=1,NCNT
C
C       LOOP OVER KAPPA VALUES FOR THIS NUCLEAR CENTER
        DO KN=1,NKAP(ICT)
C
C         IMPORT KAPPA, MAXIMUM MQN
          KAPPA = KVALS(KN,ICT)
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
C     LOOP OVER ALL ATOMIC CENTERS (USE INDEX 1000)                    C
C**********************************************************************C
C
C     LOOP OVER CENTER A
      DO 1000 ICNTA=1,NCNT
        ICNT(1) = ICNTA
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 1000 ICNTB=1,ICNTA
        ICNT(2) = ICNTB
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C     LOOP OVER CENTER C
      DO 1000 ICNTC=1,NCNT
        ICNT(3) = ICNTC
C
C       CARTESIAN COORDINATES OF CENTER C
        XYZ(1,3) = COORD(1,ICNTC)
        XYZ(2,3) = COORD(2,ICNTC)
        XYZ(3,3) = COORD(3,ICNTC)
C
C     LOOP OVER CENTER D
      DO 1000 ICNTD=1,NCNT
        ICNT(4) = ICNTD
C
C       CARTESIAN COORDINATES OF CENTER D
        XYZ(1,4) = COORD(1,ICNTD)
        XYZ(2,4) = COORD(2,ICNTD)
        XYZ(3,4) = COORD(3,ICNTD)
C
C     NUMBER OF NUCLEAR CENTERS INVOLVED IN THIS OVERLAP
      MCNT = NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
C
C     SKIP ONE-CENTER CONTRIBUTIONS (DEFER TO COULOMB1)
      IF(MCNT.EQ.1) THEN
        GOTO 1001
      ENDIF
C
C**********************************************************************C
C     LOOP OVER ALL KQN SYMMETRY TYPES (USE INDEX 2000)                C
C**********************************************************************C
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
        NBAS(1) = NFUNCT(LQN(1)+1,ICNTA)
        DO IBAS=1,NBAS(1)
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NBAS(2) = NFUNCT(LQN(2)+1,ICNTB)
        DO JBAS=1,NBAS(2)
          EXPT(JBAS,2) = EXPSET(JBAS,LQN(2)+1,ICNTB)
        ENDDO
C
C     LOOP OVER KQN(C) VALUES
      DO 2000 KC=1,NKAP(ICNTC)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR C
        KQN(3) = KVALS(KC,ICNTC)
        IF(KQN(3).GT.0) THEN
          LQN(3) = KQN(3)
        ELSE
          LQN(3) =-KQN(3)-1
        ENDIF
C
        NBAS(3) = NFUNCT(LQN(3)+1,ICNTC)
        DO KBAS=1,NBAS(3)
          EXPT(KBAS,3) = EXPSET(KBAS,LQN(3)+1,ICNTC)
        ENDDO
C
C     LOOP OVER KQN(D) VALUES
      DO 2000 KD=1,NKAP(ICNTD)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR D
        KQN(4) = KVALS(KD,ICNTD)
        IF(KQN(4).GT.0) THEN
          LQN(4) = KQN(4)
        ELSE
          LQN(4) =-KQN(4)-1
        ENDIF
C
        NBAS(4) = NFUNCT(LQN(4)+1,ICNTD)
        DO LBAS=1,NBAS(4)
          EXPT(LBAS,4) = EXPSET(LBAS,LQN(4)+1,ICNTD)
        ENDDO
C
C     THIS UNIQUELY DEFINES A FULL SET OF RC(AB|CD) INTEGRALS -- RESET
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          IRIJ(IBAS,JBAS) = 1
        ENDDO
      ENDDO
C
C**********************************************************************C
C     LOOP OVER ALL |MQN| PROJECTIONS (INDEX 3000)                     C
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
C     EQ-COEFFICIENT STARTING ADDRESSES FOR (AB) PAIR
      IABLL = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB)
      IABSS = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB)
C
C     LOOP OVER |MQN(C)| VALUES
      DO 3000 MC=1,IABS(KQN(3))
        MQN(3) = 2*MC-1
C
C     LOOP OVER |MQN(D)| VALUES
      DO 3000 MD=1,IABS(KQN(4))
        MQN(4) = 2*MD-1
C
C     EQ-COEFFICIENT STARTING ADDRESSES FOR (CD) PAIR
      ICDLL = IAD0LL(ICNTC,ICNTD,KC,KD,MC,MD)
      ICDSS = IAD0SS(ICNTC,ICNTD,KC,KD,MC,MD)
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON MQN                           C
C**********************************************************************C
C
C     SPIN PROJECTION CONSERVED ALONG Z-AXIS FOR LINEAR MOLECULES
      IF(ILIN.EQ.1) THEN
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) GOTO 3003
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) GOTO 3003
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) GOTO 3003
        GOTO 3001
      ENDIF
3003  CONTINUE
C
C**********************************************************************C
C     IDENTIFICATION OF ERI SYMMETRIES AVAILABLE TO THIS BLOCK         C
C**********************************************************************C
C
C     CALCULATE BLOCK INDICES FOR {ABCD} COMBINATIONS
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
C     COMBINED BLOCK INDEX IN A TWO-FUNCTION LIST
      IQL = (IQ1*(IQ1-1))/2 + IQ2
      IQR = (IQ3*(IQ3-1))/2 + IQ4
C
C     SKIP CONTRIBUTIONS THAT ARISE BY PERMUTATION OF INTEGRALS
      IF(IQ1.LT.IQ2) GOTO 3001
      IF(IQ3.LT.IQ4) GOTO 3001
      IF(IQL.LT.IQR) GOTO 3001
C
      IF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQL.GT.IQR) THEN
C       IQ1 > IQ2, IQ3 > IQ4, IQL > IQR
        ITSCF = 1
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQL.EQ.IQR) THEN
C       IQ1 > IQ2, IQ3 > IQ4, IQL = IQR
        ITSCF = 2
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.EQ.IQ4.AND.IQL.GT.IQR) THEN
C       IQ1 > IQ2, IQ3 = IQ4, IQL > IQR
        ITSCF = 3
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.GT.IQ4.AND.IQL.GT.IQR) THEN
C       IQ1 = IQ2, IQ3 > IQ4, IQL > IQR
        ITSCF = 4
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQL.GT.IQR) THEN
C       IQ1 = IQ2, IQ3 = IQ4, IQL > IQR
        ITSCF = 5
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQL.EQ.IQR) THEN
C       IQ1 = IQ2, IQ3 = IQ4, IQL = IQR
        ITSCF = 6
      ELSE
C       ALL OTHER CASES GENERATED BY THE ABOVE, SO SKIP
        GOTO 3001
      ENDIF
C
C     READ IN FLAG VALUES FROM ISCF DATA BLOCK
      DO ISYM=1,11
        IFLG(ISYM) = ISCF(ISYM,ITSCF)
      ENDDO
C
C     INCLUDE SPECIAL CASES FOR MATCHING BLOCKS
      IF(IQ2.EQ.IQ3) THEN
C       IQ2 = IQ3
        IF(ITSCF.EQ.1.OR.ITSCF.EQ.3.OR.ITSCF.EQ.4) THEN
          IFLG( 9) = 1
        ENDIF
      ENDIF
C
      IF(ITSCF.EQ.1) THEN
        IF(IQ1.EQ.IQ3) THEN
C         IQ1 = IQ3
          IFLG(10) = 1
        ELSEIF(IQ2.EQ.IQ4) THEN
C         IQ2 = IQ4
          IFLG(11) = 1
        ENDIF
      ENDIF
C
C     EQ-COEFFICIENT PHASE FACTORS FOR PERMUTATION OF R-INTEGRALS
      PAB1 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PAB2 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      PCD1 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PCD2 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
C**********************************************************************C
C     LOOP OVER COMPONENT OVERLAP OPTIONS (INDEX 4000)                 C
C**********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR A AND B: TT = LL (1) or SS (4)
      DO 4000 IT1=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITN(1) = IT1
C
C       CALCULATE STARTING ADDRESS
        IF(IT1.EQ.1) THEN
          NADDAB = 0
        ELSE
          NADDAB = NSHIFT
        ENDIF
C
C       FOCK ADDRESS FOR EACH BASIS FUNCTION (WITH SPIN PROJECTION)
        NA1 = LARGE(ICNTA,KA,2*MA-1) + NADDAB
        NA2 = LARGE(ICNTA,KA,2*MA  ) + NADDAB
        NB1 = LARGE(ICNTB,KB,2*MB-1) + NADDAB
        NB2 = LARGE(ICNTB,KB,2*MB  ) + NADDAB
C
C       FLAG READ-IN OF E0(AB) COEFFICIENTS FOR THIS COMPONENT LABEL
        IEAB = 1
C
C     LOOP OVER COMPONENT LABEL FOR C AND D: T'T' = LL (1) or SS (4)
      DO 4000 IT2=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITN(2) = IT2
C
C       CALCULATE STARTING ADDRESS
        IF(IT2.EQ.1) THEN
          NADDCD = 0
        ELSE
          NADDCD = NSHIFT
        ENDIF
C
C       FOCK ADDRESS FOR EACH BASIS FUNCTION (WITH SPIN PROJECTION)
        NC1 = LARGE(ICNTC,KC,2*MC-1) + NADDCD
        NC2 = LARGE(ICNTC,KC,2*MC  ) + NADDCD
        ND1 = LARGE(ICNTD,KD,2*MD-1) + NADDCD
        ND2 = LARGE(ICNTD,KD,2*MD  ) + NADDCD
C
C       FLAG READ-IN OF E0(CD) COEFFICIENTS FOR THIS COMPONENT LABEL
        IECD = 1
C
C     COMPONENT OVERLAP INDEX {(LL|LL)=1,(LL|SS)=2,(SS|LL)=3,(SS|SS)=4}
      ITT = (2*IT1+IT2)/3
C
C     STAGE 1: INCLUDE ONLY (LL|LL) REPULSION INTEGRALS
      IF(ILEV.EQ.1.AND.ITT.GT.1) THEN
        GOTO 4001
      ENDIF
C
C     STAGE 2: INCLUDE ONLY (LL|SS) AND (SS|LL) REPULSION INTEGRALS
      IF(ILEV.EQ.2.AND.ITT.GT.3) THEN
        GOTO 4001
      ENDIF
C
C     UPDATE COUNTER FOR NUMBER OF CLASSES
      N2EB(MCNT,ITT) = N2EB(MCNT,ITT)+1
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS (IBAS,JBAS) TO CONSTRUCT GMAT/QMAT     C
C**********************************************************************C
C
C     START OF PARALLEL REGION
C!$OMP PARALLEL DO COLLAPSE(2)
C!$OMP&  PRIVATE(RR,ISKP,IFLG)
C!$OMP&  SHARED(XYZ,KQN,MQN,EXPT,NBAS,ITN)
C
C     LOOP OVER ELEMENTS OF FOCK MATRIX BLOCK (IBAS,JBAS)
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
C
C         RECORD TIME AT START OF (IBAS,JBAS) BATCH
          CALL CPU_TIME(TI)
C
C         UPDATE COUNTER FOR NUMBER OF BATCHES
          N2ET(MCNT,ITT) = N2ET(MCNT,ITT) + 1
C
C         RESET SCREENING COUNTERS
          ISKP = 0
C
C         SCHWARZ SCREENING (ECONOMIC ONLY WHEN SCREENING FRACTION BIG)
          CALL CPU_TIME(TBCH1)
          IF(ITT.EQ.4) THEN
            CALL SCHWARZ(GDSC,ICNT,KQN,MQN,NBAS,ITN,IFLG,ISKP)
          ENDIF
          CALL CPU_TIME(T2)
          TC2S = TC2S + T2 - TBCH1
C
C         CONDITIONAL TO SKIP THIS BATCH
          IF(ISKP.EQ.0) THEN
C
C           GENERATE BATCH OF ELECTRON REPULSION INTEGRALS
            CALL ERI(RR,XYZ,KQN,MQN,EXPT,NBAS,ITN,IBAS,JBAS)
C
C           MULTIPLY BY DENSITY ELEMENTS AND ADD TO GMAT/QMAT
            CALL COULMAT(RR,IFLG,NBAS,TCMC)
C
          ELSE
C
C           ADD TO NUMBER OF SKIPPED BATCHES
            N2ES(MCNT,ITT) = N2ES(MCNT,ITT)+1
C
          ENDIF
C
C         RECORD TIME AT END OF (IBAS,JBAS) BATCH
          CALL CPU_TIME(TF)
          T2ES(MCNT,ITT) = T2ES(MCNT,ITT) + TF - TI

C       END LOOP OVER (IBAS,JBAS)
        ENDDO
      ENDDO
C
C     END OF PARALLEL REGION
C!$OMP END PARALLEL DO
C
4001  CONTINUE
4000  CONTINUE
3001  CONTINUE
3000  CONTINUE
2000  CONTINUE
1001  CONTINUE
1000  CONTINUE
C
C**********************************************************************C
C     COMPLETE CONSTRUCTION OF ALL MATRICES BY CONJUGATION.            C
C**********************************************************************C
C
C     LOOP OVER LOWER TRIANGLE OF EACH TT' BLOCK
      DO J=1,NDIM-NSHIFT
        DO I=1,J
C
C         SMALL-COMPONENT ADDRESSES
          K = I + NSHIFT
          L = J + NSHIFT
C
C         SKIP DIAGONAL PARTS OF EACH SUB-BLOCK
          IF(LABICN(I).NE.LABICN(J)) GOTO 400
          IF(LABKQN(I).NE.LABKQN(J)) GOTO 400
          IF(IABS(LABMQN(I)).NE.IABS(LABMQN(J))) GOTO 400
          GOTO 401
400       CONTINUE
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF LL BLOCK
          GDIR(I,J) = GDIR(I,J) + DCONJG(GDIR(J,I))
          GDIR(J,I) =             DCONJG(GDIR(I,J))
          GXCH(I,J) = GXCH(I,J) + DCONJG(GXCH(J,I))
          GXCH(J,I) =             DCONJG(GXCH(I,J))
          QDIR(I,J) = QDIR(I,J) + DCONJG(QDIR(J,I))
          QDIR(J,I) =             DCONJG(QDIR(I,J))
          QXCH(I,J) = QXCH(I,J) + DCONJG(QXCH(J,I))
          QXCH(J,I) =             DCONJG(QXCH(I,J))
C
C         IF HMLTN = 'NORL' SKIP THE NEXT FEW CALCULATIONS
          IF(HMLTN.EQ.'NORL') GOTO 401
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF SS BLOCK
          GDIR(K,L) = GDIR(K,L) + DCONJG(GDIR(L,K))
          GDIR(L,K) =             DCONJG(GDIR(K,L))
          GXCH(K,L) = GXCH(K,L) + DCONJG(GXCH(L,K))
          GXCH(L,K) =             DCONJG(GXCH(K,L))
          QDIR(K,L) = QDIR(K,L) + DCONJG(QDIR(L,K))
          QDIR(L,K) =             DCONJG(QDIR(K,L))
          QXCH(K,L) = QXCH(K,L) + DCONJG(QXCH(L,K))
          QXCH(L,K) =             DCONJG(QXCH(K,L))
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF LS BLOCK
          GXCH(I,L) = GXCH(I,L) + DCONJG(GXCH(L,I))
          GXCH(L,I) =             DCONJG(GXCH(I,L))
          QXCH(I,L) = QXCH(I,L) + DCONJG(QXCH(L,I))
          QXCH(L,I) =             DCONJG(QXCH(I,L))
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF SL BLOCK
          GXCH(K,J) = GXCH(K,J) + DCONJG(GXCH(J,K))
          GXCH(J,K) =             DCONJG(GXCH(K,J))
          QXCH(K,J) = QXCH(K,J) + DCONJG(QXCH(J,K))
          QXCH(J,K) =             DCONJG(QXCH(K,J))
C
401       CONTINUE
        ENDDO
      ENDDO
C
C     MULTIPLY OPEN MATRIX BY ANGULAR COEFFICIENTS
      DO J=1,NDIM
        DO I=1,NDIM
          QDIR(I,J) = ACFF*QDIR(I,J)
          QXCH(I,J) = BCFF*QXCH(I,J)
        ENDDO
      ENDDO
C
      RETURN
      END

