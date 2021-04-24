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
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=7000000)
C
      CHARACTER*4 HMLTN
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4),LQN(4),ITN(2)
      DIMENSION ICNT(4)
      DIMENSION INDEX(MCT,-(MKP+1)/2:(MKP+1)/2,MKP)
      DIMENSION ISCF(11,6),IFLG(11),ISCR(11)
C
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
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
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
      COMMON/TSCF/TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRR,TCC1,TCC2,TCMC,
     &            TB1B,TB1R,TB1F,TB1M,TBEC,TBRM,TBRR,TBC1,TBC2,TBMC,
     &            THMX,TC1T,TC2T,TB1T,TB2T,TEIG,TSCR,TTOT,
     &            TC1S,TC2S,TB1S,TB2S
C
      DATA SENS/1.0D-12/
      DATA ISCF/1,1,1,1,1,1,1,1,0,0,0,
     &          1,1,0,0,1,1,1,0,0,0,0,
     &          1,0,1,1,1,0,1,0,0,0,0,
     &          1,1,1,0,1,1,0,0,0,0,0,
     &          1,0,1,0,1,0,0,0,0,0,0,
     &          1,0,0,0,1,0,0,0,0,0,0/
C
C     SAVED BATCHES OF R(AB|CD) INTEGRALS IN ERI
      IF(IEQS.EQ.0) THEN
        IERC = 0
      ELSE
        IERC = 1
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
C     CONSTRUCT ORDERED INDEX SYSTEM FOR ALL POSSIBLE {XYZ,KQN,MQN}
      ICOUNT = 0
C
C     LOOP OVER ATOMIC CENTERS
      DO IC=1,NCNT
C
C       LOOP OVER KAPPA VALUES FOR THIS ATOMIC CENTER
        DO KN=1,NKAP(IC)
C
C         IMPORT KAPPA, MAXIMUM MQN
          KAPPA = KVALS(KN,IC)
          MJMAX = 2*IABS(KAPPA)-1
C
C         LOOP OVER MQN VALUES AND RECORD INDEX
          DO MJ=1,MJMAX,2
            ICOUNT = ICOUNT+1
            INDEX(IC,KAPPA,MJ) = ICOUNT
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     LOOP OVER COMPONENT OVERLAP OPTIONS (INDEX 1000)                 C
C**********************************************************************C
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
C     LOOP OVER COMPONENT LABEL FOR A AND B: TT = LL (1) or SS (4)
      DO 1000 IT1=ITSTRT,ITSTOP,ITSKIP
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
C     LOOP OVER COMPONENT LABEL FOR C AND D: T'T' = LL (1) or SS (4)
      DO 1000 IT2=ITSTRT,ITSTOP,ITSKIP
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
C     COMPONENT OVERLAP INDEX {(LL|LL)=1,(LL|SS)=2,(SS|LL)=3,(SS|SS)=4}
      ITT = (2*IT1+IT2)/3
C
C**********************************************************************C
C     LOOP OVER ATOMIC CENTERS WITH KQN A AND B (USE INDEX 2000)       C
C**********************************************************************C
C
C     LOOP OVER CENTER A
      DO 2000 ICNTA=1,NCNT
        ICNT(1) = ICNTA
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 2000 ICNTB=1,ICNTA
        ICNT(2) = ICNTB
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 0
        ENDIF
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
        DO JBAS=1, NBAS(2)
          EXPT(JBAS,2) = EXPSET(JBAS,LQN(2)+1,ICNTB)
        ENDDO
C
C**********************************************************************C
C     LOOP OVER ATOMIC CENTERS WITH KQN C AND D (USE INDEX 3000)       C
C**********************************************************************C
C
C     LOOP OVER CENTER C
      DO 3000 ICNTC=1,NCNT
        ICNT(3) = ICNTC
C
C       CARTESIAN COORDINATES OF CENTER C
        XYZ(1,3) = COORD(1,ICNTC)
        XYZ(2,3) = COORD(2,ICNTC)
        XYZ(3,3) = COORD(3,ICNTC)
C
C     LOOP OVER CENTER D
      DO 3000 ICNTD=1,NCNT
        ICNT(4) = ICNTD
C
C       CARTESIAN COORDINATES OF CENTER D
        XYZ(1,4) = COORD(1,ICNTD)
        XYZ(2,4) = COORD(2,ICNTD)
        XYZ(3,4) = COORD(3,ICNTD)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER C AND D
        IF(ICNTC.EQ.ICNTD) THEN
          INUCCD = 1
        ELSE
          INUCCD = 0
        ENDIF
C
C     NUMBER OF NUCLEAR CENTERS INVOLVED IN THIS OVERLAP
      MCNT = NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
C      
C     LOOP OVER KQN(C) VALUES
      DO 3000 KC=1,NKAP(ICNTC)
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
      DO 3000 KD=1,NKAP(ICNTD)
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
C *** STAGES: DECISION TREE FOR SKIPPING MULTI-CENTER CONTRIBUTIONS
      IF(MCNT.EQ.1) THEN
C
C >>    ATOM-CENTERED INTEGRALS: DEFER CALCULATION TO COULOMB1 ROUTINE
        GOTO 3001
C
      ELSEIF(MCNT.NE.1) THEN
C
C >>    STAGE 1: INCLUDE ONLY (LL|LL) REPULSION INTEGRALS
        IF(ILEV.EQ.1.AND.IT1+IT2.GT.2) THEN
          GOTO 3001
        ENDIF
C
C >>    STAGE 2: INCLUDE ONLY (LL|SS) AND (SS|LL) REPULSION INTEGRALS
        IF(ILEV.EQ.2.AND.IT1+IT2.GT.5) THEN
          GOTO 3001
        ENDIF
C
C *** END OF ILEV DECISION TREE
      ENDIF
C
C     RESET RC(AB|CD) CLASS CALCULATION INDICATORS
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          IRIJ(IBAS,JBAS) = 1
        ENDDO
      ENDDO
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON XYZ, KQN AND LQN              C
C**********************************************************************C
C
C     LINEAR MOLECULE OPTION
      IF(ZPROJ(XYZ).LT.SENS) THEN
        ILIN = 1
      ELSE
        ILIN = 0
      ENDIF
C
C**********************************************************************C
C     LOOP OVER ALL |MQN| VALUES (INDEX 4000)                          C
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
C     CALCULATE NEW BLOCK OF E(AB|  ) COEFFS AT NEXT OPPORTUNITY
      IEAB  = 1
      IABLL = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB)
      IABSS = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB)
C
C     LOOP OVER |MQN(C)| VALUES
      DO 4000 MC=1,IABS(KQN(3))
        MQN(3) = 2*MC-1
C
C     LOOP OVER |MQN(D)| VALUES
      DO 4000 MD=1,IABS(KQN(4))
        MQN(4) = 2*MD-1
C
C     CALCULATE NEW BLOCK OF E(CD|  ) COEFFS AT NEXT OPPORTUNITY
      IECD  = 1
      ICDLL = IAD0LL(ICNTC,ICNTD,KC,KD,MC,MD)
      ICDSS = IAD0SS(ICNTC,ICNTD,KC,KD,MC,MD)
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON MQN                           C
C**********************************************************************C
C
C     SPIN PROJECTION CONSERVED ALONG Z-AXIS
      IF(ILIN.EQ.1) THEN
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) GOTO 5503
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) GOTO 5503
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) GOTO 5503
        GOTO 4001
      ENDIF
5503  CONTINUE
C
C**********************************************************************C
C     INDEX ASSIGNMENT AND IDENTIFICATION OF ERI SYMMETRY CLASSES      C
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
      IF(IQ1.LT.IQ2) GOTO 4001
      IF(IQ3.LT.IQ4) GOTO 4001
      IF(IQL.LT.IQR) GOTO 4001
C
      IF(IQ1.GT.IQ2) THEN
C       IQ1 > IQ2
        IF(IQ3.GT.IQ4) THEN
C         IQ3 > IQ4
          IF(IQL.GT.IQR) THEN
C           IQL > IQR
            ITSCF = 1
          ELSEIF(IQL.EQ.IQR) THEN
C           IQL = IQR
            ITSCF = 2
          ENDIF
        ELSEIF(IQ3.EQ.IQ4) THEN
C         IQ3 = IQ4
          IF(IQL.GT.IQR) THEN
C           IQL > IQR
            ITSCF = 3
          ELSE
C           IQL = IQR
            GOTO 4001
          ENDIF
        ENDIF
      ELSEIF(IQ1.EQ.IQ2) THEN
C       IQ1 = IQ2
        IF(IQ3.GT.IQ4) THEN
C         IQ3 > IQ4
          IF(IQL.GT.IQR) THEN
C           IQL > IQR          
            ITSCF = 4
          ELSE
C           IQL = IQR
            GOTO 4001
          ENDIF
        ELSEIF(IQ3.EQ.IQ4) THEN
C         IQ3 = IQ4
          IF(IQL.GT.IQR) THEN
C           IQL > IQR
            ITSCF = 5
          ELSEIF(IQL.EQ.IQR) THEN
C           IQL = IQR
            ITSCF = 6
          ENDIF
        ENDIF
      ENDIF
C
C     READ IN FLAG VALUES FROM ISCF DATA BLOCK
      DO N=1,11
        IFLG(N) = ISCF(N,ITSCF)
      ENDDO
C
C     INCLUDE SPECIAL CASES FOR MATCHING BLOCKS
      IF(ITSCF.EQ.1) THEN
        IF(IQ1.EQ.IQ3) THEN
C         IQ1 = IQ3
          IFLG(10) = 1
        ENDIF
        IF(IQ2.EQ.IQ4) THEN
C         IQ2 = IQ3
          IFLG(11) = 1
        ENDIF
        IF(IQ2.EQ.IQ3) THEN
C         IQ2 = IQ3
          IFLG( 9) = 1
        ENDIF
      ELSEIF(ITSCF.EQ.3) THEN
        IF(IQ2.EQ.IQ3) THEN
C         IQ2 = IQ3
          IFLG( 9) = 1
        ENDIF
      ELSEIF(ITSCF.EQ.4) THEN
        IF(IQ2.EQ.IQ3) THEN
C         IQ2 = IQ3
          IFLG( 9) = 1
        ENDIF
      ENDIF
C
C     UPDATE COUNTER FOR NUMBER OF CLASSES
      N2EB(MCNT,ITT) = N2EB(MCNT,ITT) + 1
C
C**********************************************************************C
C     PHASE FACTORS APPROPRIATE TO INTEGRAL SYMMETRY CLASSES           C
C**********************************************************************C
C
C     EQ-COEFFICIENT PHASE FACTORS FOR PERMUTATION OF R-INTEGRALS
      PAB1 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PAB2 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      PCD1 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PCD2 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
C     FOCK ADDRESS FOR EACH BASIS FUNCTION (WITH SPIN PROJECTION)
      NA1 = LARGE(ICNTA,KA,2*MA-1) + NADDAB
      NB1 = LARGE(ICNTB,KB,2*MB-1) + NADDAB
      NC1 = LARGE(ICNTC,KC,2*MC-1) + NADDCD
      ND1 = LARGE(ICNTD,KD,2*MD-1) + NADDCD
C
      NA2 = LARGE(ICNTA,KA,2*MA  ) + NADDAB
      NB2 = LARGE(ICNTB,KB,2*MB  ) + NADDAB
      NC2 = LARGE(ICNTC,KC,2*MC  ) + NADDCD
      ND2 = LARGE(ICNTD,KD,2*MD  ) + NADDCD
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS IN BLOCKS A AND B (INDEX 5000)         C
C**********************************************************************C
C
C     LOOP OVER ELEMENTS OF FOCK MATRIX BLOCK
      DO 5000 IBAS=1,NBAS(1)
      DO 5000 JBAS=1,NBAS(2)
C
C     UPDATE COUNTER FOR NUMBER OF BATCHES
      N2ET(MCNT,ITT) = N2ET(MCNT,ITT) + 1
C
C     RESET SCREENING COUNTERS
      DO ISYM=1,11
        ISCR(ISYM) = 1
      ENDDO
C
C     SCHWARZ SCREENING: SWITCH OFF IFLG VALUES OR SKIP ERI ENTIRELY
      CALL CPU_TIME(TBCH1)
c     CALL SCHWARZ(ICNT,KQN,MQN,NBAS,ITN,IBAS,JBAS,IFLG,ISCR,ISKP)
      CALL CPU_TIME(TX2)
      TC2S = TC2S + TX2 - TBCH1
C
C     CONDITIONAL TO SKIP THIS BATCH
      IF(ISKP.EQ.1) THEN
        N2ES(MCNT,ITT) = N2ES(MCNT,ITT)+1
        GOTO 5001
      ENDIF
C
C**********************************************************************C
C     INCLUDE A BATCH OF ELECTRON REPULSION INTEGRALS INTO GMAT/QMAT.  C
C     DEPENDING ON CURRENT COMBINATION OF MQN VALUES, TAKE ADVANTAGE   C
C     OF INTEGRAL PERMUTATION SYMMETRY (HENCE MINIMISE CALLS TO ERI):  C
C              ( MA, MB|-MD,-MC) =     PCD*( MA, MB| MC, MD)           C
C              (-MB,-MA| MC, MD) = PAB*    ( MA, MB| MC, MD)           C
C              ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)           C
C**********************************************************************C
C
C     GENERATE BATCH OF ELECTRON REPULSION INTEGRALS
      CALL ERI(RR,XYZ,KQN,MQN,EXPT,NBAS,ITN,IBAS,JBAS)
C
C     START OF MATRIX ARITHMETIC CODE
      CALL CPU_TIME(T1)
C
C     IF THE COMBINATION OF ATOMIC CENTERS RUNS ALONG Z-AXIS
C
C     1ST CASE (DIRECT):  ( MA, MB| MC, MD)                         
      IF(IFLG(1).EQ.1.AND.ISCR(1).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 1)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 4)*DENC(NC2+KBAS,ND2+LBAS)
C
            GDIR(NA2+IBAS,NB2+JBAS) = GDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(M,13)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,16)*DENC(NC2+KBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     2ND CASE (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(2).EQ.1.AND.ISCR(2).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENC(ND2+LBAS,NC2+KBAS)
C
            GDIR(NA2+IBAS,NB2+JBAS) = GDIR(NA2+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(M,16)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENC(ND2+LBAS,NC2+KBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     3RD CASE (DIRECT):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(3).EQ.1.AND.ISCR(3).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NC1+KBAS,ND1+LBAS) = GDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 1)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M,13)*DENC(NA2+IBAS,NB2+JBAS)
C
            GDIR(NC2+KBAS,ND2+LBAS) = GDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(M, 4)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M,16)*DENC(NA2+IBAS,NB2+JBAS)
C 
          ENDDO
        ENDDO
      ENDIF
C
C     4TH CASE (DIRECT):  ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(4).EQ.1.AND.ISCR(4).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NC1+KBAS,ND1+LBAS) = GDIR(NC1+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,13)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENC(NB2+JBAS,NA2+IBAS)
C
            GDIR(NC2+KBAS,ND2+LBAS) = GDIR(NC2+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,16)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENC(NB2+JBAS,NA2+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     5TH CASE (EXCHNG):  ( MA, MB| MC, MD)
      IF(IFLG(5).EQ.1.AND.ISCR(5).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C  
            GXCH(NA1+IBAS,ND1+LBAS) = GXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 1)*DENC(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 7)*DENC(NC2+KBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,ND2+LBAS) = GXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(M,10)*DENC(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M,16)*DENC(NC2+KBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     6TH CASE (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(6).EQ.1.AND.ISCR(6).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NA1+IBAS,NC1+KBAS) = GXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 4)*DENC(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 7)*DENC(ND2+LBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,NC2+KBAS) = GXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,10)*DENC(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,13)*DENC(ND2+LBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     7TH CASE (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(7).EQ.1.AND.ISCR(7).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NB1+JBAS,ND1+LBAS) = GXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,13)*DENC(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENC(NC2+KBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,ND2+LBAS) = GXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M,10)*DENC(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENC(NC2+KBAS,NA2+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     8TH CASE (EXCHNG):  ( MB, MA| MD, MC) = PAB*    (-MA,-MB| MC, MD)
C                                           = PAB*PCD*(-MA,-MB|-MC,-MD)
      IF(IFLG(8).EQ.1.AND.ISCR(8).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NB1+JBAS,NC1+KBAS) = GXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB1*PCD1*RR(M,16)*DENC(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(M, 7)*DENC(ND2+LBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,NC2+KBAS) = GXCH(NB2+JBAS,NC2+KBAS)
     &                     + PAB1*PCD1*RR(M, 1)*DENC(ND2+LBAS,NA2+IBAS)
     &                     + PAB2*PCD2*RR(M,10)*DENC(ND1+LBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     9TH CASE (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(9).EQ.1.AND.ISCR(9).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NC1+KBAS,NB1+JBAS) = GXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 1)*DENC(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M,10)*DENC(NA2+IBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NB2+JBAS) = GXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(M, 7)*DENC(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M,16)*DENC(NA2+IBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     10TH CASE (EXCHNG): ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(10).EQ.1.AND.ISCR(10).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NC1+KBAS,NA1+IBAS) = GXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,13)*DENC(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M,10)*DENC(NB2+JBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NA2+IBAS) = GXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 7)*DENC(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M, 4)*DENC(NB2+JBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     11TH CASE (EXCHNG): ( MD, MC| MA, MB) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(11).EQ.1.AND.ISCR(11).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(ND1+LBAS,NB1+JBAS) = GXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENC(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENC(NA2+IBAS,NC2+KBAS)
C
            GXCH(ND2+LBAS,NB2+JBAS) = GXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M, 7)*DENC(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENC(NA2+IBAS,NC2+KBAS)
C
          ENDDO
        ENDDO
      ENDIF      
C
C     ADD MORE CONTRIBUTIONS FOR GENERAL MOLECULAR BATCH
      IF(ILIN.EQ.1) GOTO 5002
C
C     1ST CASE (DIRECT):  ( MA, MB| MC, MD)
      IF(IFLG(1).EQ.1.AND.ISCR(1).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 2)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 3)*DENC(NC2+KBAS,ND1+LBAS)
C
            GDIR(NA1+IBAS,NB2+JBAS) = GDIR(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 5)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 6)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 7)*DENC(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 8)*DENC(NC2+KBAS,ND2+LBAS)
C
            GDIR(NA2+IBAS,NB1+JBAS) = GDIR(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M, 9)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,10)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,11)*DENC(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M,12)*DENC(NC2+KBAS,ND2+LBAS)
C
            GDIR(NA2+IBAS,NB2+JBAS) = GDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(M,14)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,15)*DENC(NC2+KBAS,ND1+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     2ND CASE (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(2).EQ.1.AND.ISCR(2).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 2)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 3)*DENC(ND2+LBAS,NC1+KBAS)
C
            GDIR(NA1+IBAS,NB2+JBAS) = GDIR(NA1+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 8)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 7)*DENC(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 5)*DENC(ND2+LBAS,NC2+KBAS)
C
            GDIR(NA2+IBAS,NB1+JBAS) = GDIR(NA2+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,12)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENC(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENC(ND2+LBAS,NC2+KBAS)
C
            GDIR(NA2+IBAS,NB2+JBAS) = GDIR(NA2+IBAS,NB2+JBAS)
     &                     +      PCD2*RR(M,14)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENC(ND2+LBAS,NC1+KBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     3RD CASE (DIRECT):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(3).EQ.1.AND.ISCR(3).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NC1+KBAS,ND1+LBAS) = GDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 5)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 9)*DENC(NA2+IBAS,NB1+JBAS)
C
            GDIR(NC1+KBAS,ND2+LBAS) = GDIR(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 2)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 6)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,10)*DENC(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,14)*DENC(NA2+IBAS,NB2+JBAS)
C
            GDIR(NC2+KBAS,ND1+LBAS) = GDIR(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 3)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 7)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,11)*DENC(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,15)*DENC(NA2+IBAS,NB2+JBAS)
C
            GDIR(NC2+KBAS,ND2+LBAS) = GDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(M, 8)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,12)*DENC(NA2+IBAS,NB1+JBAS)
C 
          ENDDO
        ENDDO
      ENDIF
C
C     4TH CASE (DIRECT):  ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(4).EQ.1.AND.ISCR(4).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NC1+KBAS,ND1+LBAS) = GDIR(NC1+KBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 5)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 9)*DENC(NB2+JBAS,NA1+IBAS)
C
            GDIR(NC1+KBAS,ND2+LBAS) = GDIR(NC1+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,14)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,10)*DENC(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 2)*DENC(NB2+JBAS,NA2+IBAS)
C
            GDIR(NC2+KBAS,ND1+LBAS) = GDIR(NC2+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,15)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENC(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENC(NB2+JBAS,NA2+IBAS)
C
            GDIR(NC2+KBAS,ND2+LBAS) = GDIR(NC2+KBAS,ND2+LBAS)
     &                     + PAB2*     RR(M, 8)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENC(NB2+JBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     5TH CASE (EXCHNG):  ( MA, MB| MC, MD)
      IF(IFLG(5).EQ.1.AND.ISCR(5).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C  
            GXCH(NA1+IBAS,ND1+LBAS) = GXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 5)*DENC(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 3)*DENC(NC2+KBAS,NB1+JBAS)
C
            GXCH(NA1+IBAS,ND2+LBAS) = GXCH(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M, 2)*DENC(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 6)*DENC(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 4)*DENC(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M, 8)*DENC(NC2+KBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,ND1+LBAS) = GXCH(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M, 9)*DENC(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M,13)*DENC(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M,11)*DENC(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M,15)*DENC(NC2+KBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,ND2+LBAS) = GXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(M,14)*DENC(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M,12)*DENC(NC2+KBAS,NB1+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     6TH CASE (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(6).EQ.1.AND.ISCR(6).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NA1+IBAS,NC1+KBAS) = GXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 8)*DENC(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M, 3)*DENC(ND2+LBAS,NB1+JBAS)
C
            GXCH(NA1+IBAS,NC2+KBAS) = GXCH(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 2)*DENC(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 6)*DENC(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 1)*DENC(ND2+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 5)*DENC(ND2+LBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,NC1+KBAS) = GXCH(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,12)*DENC(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,16)*DENC(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M,11)*DENC(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M,15)*DENC(ND2+LBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,NC2+KBAS) = GXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,14)*DENC(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 9)*DENC(ND2+LBAS,NB1+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     7TH CASE (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(7).EQ.1.AND.ISCR(7).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NB1+JBAS,ND1+LBAS) = GXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 5)*DENC(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(M,15)*DENC(NC2+KBAS,NA1+IBAS)
C
            GXCH(NB1+JBAS,ND2+LBAS) = GXCH(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,14)*DENC(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENC(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(M,16)*DENC(NC2+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 8)*DENC(NC2+KBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,ND1+LBAS) = GXCH(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 9)*DENC(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENC(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENC(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENC(NC2+KBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,ND2+LBAS) = GXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M, 2)*DENC(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENC(NC2+KBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     8TH CASE (EXCHNG):  ( MB, MA| MD, MC) = PAB*    (-MA,-MB| MC, MD)
C                                           = PAB*PCD*(-MA,-MB|-MC,-MD)
      IF(IFLG(8).EQ.1.AND.ISCR(8).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NB1+JBAS,NC1+KBAS) = GXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB2*PCD1*RR(M, 8)*DENC(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD2*RR(M,15)*DENC(ND2+LBAS,NA1+IBAS)
C
            GXCH(NB1+JBAS,NC2+KBAS) = GXCH(NB1+JBAS,NC2+KBAS)
     &                     + PAB1*PCD2*RR(M,14)*DENC(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(M, 6)*DENC(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD1*RR(M,13)*DENC(ND2+LBAS,NA1+IBAS)
     &                     + PAB2*PCD1*RR(M, 5)*DENC(ND2+LBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,NC1+KBAS) = GXCH(NB2+JBAS,NC1+KBAS)
     &                     + PAB2*PCD1*RR(M,12)*DENC(ND1+LBAS,NA1+IBAS)
     &                     + PAB1*PCD1*RR(M, 4)*DENC(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD2*RR(M,11)*DENC(ND2+LBAS,NA1+IBAS)
     &                     + PAB1*PCD2*RR(M, 3)*DENC(ND2+LBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,NC2+KBAS) = GXCH(NB2+JBAS,NC2+KBAS)
     &                     + PAB1*PCD2*RR(M, 2)*DENC(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD1*RR(M, 9)*DENC(ND2+LBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     9TH CASE (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(9).EQ.1.AND.ISCR(9).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NC1+KBAS,NB1+JBAS) = GXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 2)*DENC(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M, 9)*DENC(NA2+IBAS,ND1+LBAS)
C
            GXCH(NC1+KBAS,NB2+JBAS) = GXCH(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 5)*DENC(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 6)*DENC(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,13)*DENC(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M,14)*DENC(NA2+IBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NB1+JBAS) = GXCH(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M, 3)*DENC(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 4)*DENC(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,11)*DENC(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M,12)*DENC(NA2+IBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NB2+JBAS) = GXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(M, 8)*DENC(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,15)*DENC(NA2+IBAS,ND1+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     10TH CASE (EXCHNG): ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(10).EQ.1.AND.ISCR(10).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NC1+KBAS,NA1+IBAS) = GXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,14)*DENC(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M, 9)*DENC(NB2+JBAS,ND1+LBAS)
C
            GXCH(NC1+KBAS,NA2+IBAS) = GXCH(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 5)*DENC(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 6)*DENC(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M, 1)*DENC(NB2+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M, 2)*DENC(NB2+JBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NA1+IBAS) = GXCH(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,15)*DENC(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,16)*DENC(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M,11)*DENC(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M,12)*DENC(NB2+JBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NA2+IBAS) = GXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 8)*DENC(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M, 3)*DENC(NB2+JBAS,ND1+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     11TH CASE (EXCHNG): ( MD, MC| MA, MB) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(11).EQ.1.AND.ISCR(11).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(ND1+LBAS,NB1+JBAS) = GXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 2)*DENC(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(M,12)*DENC(NA2+IBAS,NC1+KBAS)
C
            GXCH(ND1+LBAS,NB2+JBAS) = GXCH(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 8)*DENC(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENC(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(M,16)*DENC(NA2+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,14)*DENC(NA2+IBAS,NC2+KBAS)
C
            GXCH(ND2+LBAS,NB1+JBAS) = GXCH(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 3)*DENC(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENC(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENC(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENC(NA2+IBAS,NC2+KBAS)
C
            GXCH(ND2+LBAS,NB2+JBAS) = GXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 5)*DENC(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENC(NA2+IBAS,NC1+KBAS)
C
          ENDDO
        ENDDO
      ENDIF      
C
C     SKIP POINT FOR LINEAR BATCH
5002  CONTINUE
C
C**********************************************************************C
C     OPEN-SHELL CALCULATIONS ONLY                                     C
C**********************************************************************C
C
      IF(NOPN.EQ.0) GOTO 5100
C
C     IF THE COMBINATION OF ATOMIC CENTERS RUNS ALONG Z-AXIS
C
C     1ST CASE (DIRECT):  ( MA, MB| MC, MD)
      IF(IFLG(1).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NA1+IBAS,NB1+JBAS) = QDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 1)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 4)*DENO(NC2+KBAS,ND2+LBAS)
C
            QDIR(NA2+IBAS,NB2+JBAS) = QDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(M,13)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,16)*DENO(NC2+KBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     2ND CASE (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(2).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NA1+IBAS,NB1+JBAS) = QDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENO(ND2+LBAS,NC2+KBAS)
C
            QDIR(NA2+IBAS,NB2+JBAS) = QDIR(NA2+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(M,16)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENO(ND2+LBAS,NC2+KBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     3RD CASE (DIRECT):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(3).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NC1+KBAS,ND1+LBAS) = QDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 1)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M,13)*DENO(NA2+IBAS,NB2+JBAS)
C
            QDIR(NC2+KBAS,ND2+LBAS) = QDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(M, 4)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M,16)*DENO(NA2+IBAS,NB2+JBAS)
C 
          ENDDO
        ENDDO
      ENDIF
C
C     4TH CASE (DIRECT):  ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(4).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NC1+KBAS,ND1+LBAS) = QDIR(NC1+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,13)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENO(NB2+JBAS,NA2+IBAS)
C
            QDIR(NC2+KBAS,ND2+LBAS) = QDIR(NC2+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,16)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENO(NB2+JBAS,NA2+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     5TH CASE (EXCHNG):  ( MA, MB| MC, MD)
      IF(IFLG(5).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C  
            QXCH(NA1+IBAS,ND1+LBAS) = QXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 1)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 7)*DENO(NC2+KBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,ND2+LBAS) = QXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(M,10)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M,16)*DENO(NC2+KBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     6TH CASE (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(6).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NA1+IBAS,NC1+KBAS) = QXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 4)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 7)*DENO(ND2+LBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,NC2+KBAS) = QXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,10)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,13)*DENO(ND2+LBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     7TH CASE (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(7).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NB1+JBAS,ND1+LBAS) = QXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,13)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENO(NC2+KBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,ND2+LBAS) = QXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M,10)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENO(NC2+KBAS,NA2+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     8TH CASE (EXCHNG):  ( MB, MA| MD, MC) = PAB*    (-MA,-MB| MC, MD)
C                                           = PAB*PCD*(-MA,-MB|-MC,-MD)
      IF(IFLG(8).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NB1+JBAS,NC1+KBAS) = QXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB1*PCD1*RR(M,16)*DENO(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(M, 7)*DENO(ND2+LBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,NC2+KBAS) = QXCH(NB2+JBAS,NC2+KBAS)
     &                     + PAB1*PCD1*RR(M, 1)*DENO(ND2+LBAS,NA2+IBAS)
     &                     + PAB2*PCD2*RR(M,10)*DENO(ND1+LBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     9TH CASE (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(9).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NC1+KBAS,NB1+JBAS) = QXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 1)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M,10)*DENO(NA2+IBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NB2+JBAS) = QXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(M, 7)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M,16)*DENO(NA2+IBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     10TH CASE (EXCHNG): ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(10).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NC1+KBAS,NA1+IBAS) = QXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,13)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M,10)*DENO(NB2+JBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NA2+IBAS) = QXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 7)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M, 4)*DENO(NB2+JBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     11TH CASE (EXCHNG): ( MD, MC| MA, MB) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(11).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(ND1+LBAS,NB1+JBAS) = QXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENO(NA2+IBAS,NC2+KBAS)
C
            QXCH(ND2+LBAS,NB2+JBAS) = QXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M, 7)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENO(NA2+IBAS,NC2+KBAS)
C
          ENDDO
        ENDDO
      ENDIF      
C
C     ADD MORE CONTRIBUTIONS FOR GENERAL MOLECULAR BATCH
      IF(ILIN.EQ.1) GOTO 5003
C
C     1ST CASE (DIRECT):  ( MA, MB| MC, MD)
      IF(IFLG(1).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NA1+IBAS,NB1+JBAS) = QDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 2)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 3)*DENO(NC2+KBAS,ND1+LBAS)
C
            QDIR(NA1+IBAS,NB2+JBAS) = QDIR(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 5)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 6)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 7)*DENO(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 8)*DENO(NC2+KBAS,ND2+LBAS)
C
            QDIR(NA2+IBAS,NB1+JBAS) = QDIR(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M, 9)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,10)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,11)*DENO(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M,12)*DENO(NC2+KBAS,ND2+LBAS)
C
            QDIR(NA2+IBAS,NB2+JBAS) = QDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(M,14)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,15)*DENO(NC2+KBAS,ND1+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     2ND CASE (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(2).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NA1+IBAS,NB1+JBAS) = QDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 2)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 3)*DENO(ND2+LBAS,NC1+KBAS)
C
            QDIR(NA1+IBAS,NB2+JBAS) = QDIR(NA1+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 8)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 7)*DENO(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 5)*DENO(ND2+LBAS,NC2+KBAS)
C
            QDIR(NA2+IBAS,NB1+JBAS) = QDIR(NA2+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,12)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENO(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENO(ND2+LBAS,NC2+KBAS)
C
            QDIR(NA2+IBAS,NB2+JBAS) = QDIR(NA2+IBAS,NB2+JBAS)
     &                     +      PCD2*RR(M,14)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENO(ND2+LBAS,NC1+KBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     3RD CASE (DIRECT):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(3).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NC1+KBAS,ND1+LBAS) = QDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 5)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 9)*DENO(NA2+IBAS,NB1+JBAS)
C
            QDIR(NC1+KBAS,ND2+LBAS) = QDIR(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 2)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 6)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,10)*DENO(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,14)*DENO(NA2+IBAS,NB2+JBAS)
C
            QDIR(NC2+KBAS,ND1+LBAS) = QDIR(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 3)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 7)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,11)*DENO(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,15)*DENO(NA2+IBAS,NB2+JBAS)
C
            QDIR(NC2+KBAS,ND2+LBAS) = QDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(M, 8)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,12)*DENO(NA2+IBAS,NB1+JBAS)
C 
          ENDDO
        ENDDO
      ENDIF
C
C     4TH CASE (DIRECT):  ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(4).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NC1+KBAS,ND1+LBAS) = QDIR(NC1+KBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 5)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 9)*DENO(NB2+JBAS,NA1+IBAS)
C
            QDIR(NC1+KBAS,ND2+LBAS) = QDIR(NC1+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,14)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,10)*DENO(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 2)*DENO(NB2+JBAS,NA2+IBAS)
C
            QDIR(NC2+KBAS,ND1+LBAS) = QDIR(NC2+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,15)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENO(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENO(NB2+JBAS,NA2+IBAS)
C
            QDIR(NC2+KBAS,ND2+LBAS) = QDIR(NC2+KBAS,ND2+LBAS)
     &                     + PAB2*     RR(M, 8)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENO(NB2+JBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     5TH CASE (EXCHNG):  ( MA, MB| MC, MD)
      IF(IFLG(5).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C  
            QXCH(NA1+IBAS,ND1+LBAS) = QXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 5)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 3)*DENO(NC2+KBAS,NB1+JBAS)
C
            QXCH(NA1+IBAS,ND2+LBAS) = QXCH(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M, 2)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 6)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 4)*DENO(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M, 8)*DENO(NC2+KBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,ND1+LBAS) = QXCH(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M, 9)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M,13)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M,11)*DENO(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M,15)*DENO(NC2+KBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,ND2+LBAS) = QXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(M,14)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M,12)*DENO(NC2+KBAS,NB1+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     6TH CASE (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(6).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NA1+IBAS,NC1+KBAS) = QXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 8)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M, 3)*DENO(ND2+LBAS,NB1+JBAS)
C
            QXCH(NA1+IBAS,NC2+KBAS) = QXCH(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 2)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 6)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 1)*DENO(ND2+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 5)*DENO(ND2+LBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,NC1+KBAS) = QXCH(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,12)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,16)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M,11)*DENO(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M,15)*DENO(ND2+LBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,NC2+KBAS) = QXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,14)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 9)*DENO(ND2+LBAS,NB1+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     7TH CASE (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(7).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NB1+JBAS,ND1+LBAS) = QXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 5)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(M,15)*DENO(NC2+KBAS,NA1+IBAS)
C
            QXCH(NB1+JBAS,ND2+LBAS) = QXCH(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,14)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(M,16)*DENO(NC2+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 8)*DENO(NC2+KBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,ND1+LBAS) = QXCH(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 9)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENO(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENO(NC2+KBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,ND2+LBAS) = QXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M, 2)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENO(NC2+KBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     8TH CASE (EXCHNG):  ( MB, MA| MD, MC) = PAB*    (-MA,-MB| MC, MD)
C                                           = PAB*PCD*(-MA,-MB|-MC,-MD)
      IF(IFLG(8).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NB1+JBAS,NC1+KBAS) = QXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB2*PCD1*RR(M, 8)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD2*RR(M,15)*DENO(ND2+LBAS,NA1+IBAS)
C
            QXCH(NB1+JBAS,NC2+KBAS) = QXCH(NB1+JBAS,NC2+KBAS)
     &                     + PAB1*PCD2*RR(M,14)*DENO(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(M, 6)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD1*RR(M,13)*DENO(ND2+LBAS,NA1+IBAS)
     &                     + PAB2*PCD1*RR(M, 5)*DENO(ND2+LBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,NC1+KBAS) = QXCH(NB2+JBAS,NC1+KBAS)
     &                     + PAB2*PCD1*RR(M,12)*DENO(ND1+LBAS,NA1+IBAS)
     &                     + PAB1*PCD1*RR(M, 4)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD2*RR(M,11)*DENO(ND2+LBAS,NA1+IBAS)
     &                     + PAB1*PCD2*RR(M, 3)*DENO(ND2+LBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,NC2+KBAS) = QXCH(NB2+JBAS,NC2+KBAS)
     &                     + PAB1*PCD2*RR(M, 2)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD1*RR(M, 9)*DENO(ND2+LBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     9TH CASE (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(9).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NC1+KBAS,NB1+JBAS) = QXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 2)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M, 9)*DENO(NA2+IBAS,ND1+LBAS)
C
            QXCH(NC1+KBAS,NB2+JBAS) = QXCH(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 5)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 6)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,13)*DENO(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M,14)*DENO(NA2+IBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NB1+JBAS) = QXCH(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M, 3)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 4)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,11)*DENO(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M,12)*DENO(NA2+IBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NB2+JBAS) = QXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(M, 8)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,15)*DENO(NA2+IBAS,ND1+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     10TH CASE (EXCHNG): ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(10).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NC1+KBAS,NA1+IBAS) = QXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,14)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M, 9)*DENO(NB2+JBAS,ND1+LBAS)
C
            QXCH(NC1+KBAS,NA2+IBAS) = QXCH(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 5)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 6)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M, 1)*DENO(NB2+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M, 2)*DENO(NB2+JBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NA1+IBAS) = QXCH(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,15)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,16)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M,11)*DENO(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M,12)*DENO(NB2+JBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NA2+IBAS) = QXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 8)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M, 3)*DENO(NB2+JBAS,ND1+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     11TH CASE (EXCHNG): ( MD, MC| MA, MB) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(11).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(ND1+LBAS,NB1+JBAS) = QXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 2)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(M,12)*DENO(NA2+IBAS,NC1+KBAS)
C
            QXCH(ND1+LBAS,NB2+JBAS) = QXCH(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 8)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(M,16)*DENO(NA2+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,14)*DENO(NA2+IBAS,NC2+KBAS)
C
            QXCH(ND2+LBAS,NB1+JBAS) = QXCH(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 3)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENO(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENO(NA2+IBAS,NC2+KBAS)
C
            QXCH(ND2+LBAS,NB2+JBAS) = QXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 5)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENO(NA2+IBAS,NC1+KBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     SKIP POINT FOR LINEAR BATCH
5003  CONTINUE
C     SKIPPING POINT FOR CLOSED SYSTEMS
5100  CONTINUE
C     RECORD MATRIX ARITHMETIC TIME
      CALL CPU_TIME(T2)
      TCMC = TCMC+T2-T1
C     SKIPPING POINT FOR INTEGRAL SCREENING
5001  CONTINUE
C     RECORD CPU TIME AT END OF BATCH AND ADD TO TIME COUNTER
      CALL CPU_TIME(TBCH2)
      T2ES(MCNT,ITT) = T2ES(MCNT,ITT) + TBCH2 - TBCH1
C     END LOOP OVER IBAS AND JBAS
5000  CONTINUE
C     SKIPPING POINT FOR MQN SELECTION RULES AND INTEGRAL SYMMETRY
4001  CONTINUE
C     END LOOP OVER |MQN| VALUES A AND B
4000  CONTINUE
C     SKIPPING POINT FOR INCLUSION LEVELS
3001  CONTINUE
C     END LOOP OVER CENTERS AND KQNS C AND D
3000  CONTINUE
C     END LOOP OVER CENTERS AND KQNS A AND B
2000  CONTINUE
C     END LOOP OVER COMPONENT OVERLAPS
1000  CONTINUE
C
C**********************************************************************C
C     COMPLETE CONSTRUCTION OF ALL MULTI-CENTER CONTRIBUTIONS.         C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C     LOOP OVER LOWER TRIANGLE OF EACH TT' BLOCK
      DO J=1,NDIM-NSHIFT
        DO I=1,J
C
C         SMALL-COMPONENT ADDRESSES
          K = I + NSHIFT
          L = J + NSHIFT
C
C         SKIP DIAGONAL PARTS OF EACH SUB-BLOCK
          IF(LABICN(I).NE.LABICN(J)) GOTO 7000
          IF(LABKQN(I).NE.LABKQN(J)) GOTO 7000
          IF(IABS(LABMQN(I)).NE.IABS(LABMQN(J))) GOTO 7000
          GOTO 7001
7000      CONTINUE
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
          IF(HMLTN.EQ.'NORL') GOTO 7001
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
7001      CONTINUE
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
      CALL CPU_TIME(T2)
      TCMC = TCMC+T2-T1
C
      RETURN
      END
