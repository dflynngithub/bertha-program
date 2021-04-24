      SUBROUTINE BREIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C               BBBBBBB  RRRRRRR  EEEEEEEE IIII TTTTTTTT               C
C               BB    BB RR    RR EE        II     TT                  C
C               BB    BB RR    RR EE        II     TT                  C
C               BBBBBBB  RR    RR EEEEEE    II     TT                  C
C               BB    BB RRRRRRR  EE        II     TT                  C
C               BB    BB RR    RR EE        II     TT                  C
C               BBBBBBB  RR    RR EEEEEEEE IIII    TT                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  BREIT GENERATES ELECTRON INTERACTION INTEGRALS IN BATCHES AND       C
C  CALCULATES THE SCF BREIT MATRIX (B). INTEGRAL SYMMETRY IS PARTIALLY C
C  PARTIALLY EXPLOITED, BUT WITH ROOM FOR IMPROVEMENT.                 C
C  (GEOMETRIC SYMM, R-INT SYMM, E-COEFF SYMM).                         C
C -------------------------------------------------------------------- C
C  DFNOTE: THIS ROUTINE COULD BENEFIT FROM PARALLELISATION -- OPENMP.  C
C          LOOK INTO OPEN SHELL EXTENSIONS AND CASES WITH ZERO DIRECT. C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=15,MKP=9,MFL=7000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      CHARACTER*4 HMLTN
      CHARACTER*8 SHAPE
C
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
C
      DIMENSION GDSC(MDM,MDM),BDSC(MDM,MDM)
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
      DIMENSION ISCF(11,6),IFLG(11)
      DIMENSION INDEX(MCT,-MKP:MKP,2*(MKP+1)*MKP),ICNT(4)
C
      COMMON/BLOC/PAB1,PAB2,PCD1,PCD2,NA1,NB1,NC1,ND1,NA2,NB2,NC2,ND2,
     &            IBAS,JBAS,ILIN,NADDAB,NADDCD
      COMMON/EILS/EILSFL(MFL,24),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/GEOM/SHAPE
      COMMON/IBSC/IBSCR(MB2),IBMAP(MB2)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS
      COMMON/MAKE/IEAB,IECD,IRIJ(MBS,MBS)
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,IEQS,IERC,IPAR,ICOR,ILEV
      COMMON/QNMS/LABICN(MDM),LABKQN(MDM),LABMQN(MDM)
      COMMON/SCRN/F2ES(5,7),T2ES(5,7),N2EB(5,7),N2ET(5,7),N2ES(5,7)
      COMMON/SHLL/ALPHA,BETA,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN,NOELEC
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
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
C     OVERRIDE ISCF UNTIL INTEGRAL SYMMETRIES HAVE BEEN IDENTIFIED
      DO ITSCF=1,2
        DO ISYM=1,11
          ISCF(ISYM,ITSCF) = 1
        ENDDO
      ENDDO
C
C     SKIP CASES WHEN IQL=IQR
      ISCF( 3,2) = 0
      ISCF( 4,2) = 0
      ISCF( 9,2) = 0
      ISCF(10,2) = 0
      ISCF(11,2) = 0
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
          BDIR(I,J) = DCMPLX(0.0D0,0.0D0)
          BXCH(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     TWO-ELECTRON COMPONENT OVERLAP INDEX
      ITT = 5
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
          IF(KAPPA.GT.0) THEN
            LQNN = KAPPA
          ELSE
            LQNN =-KAPPA-1
          ENDIF
          NFUN = NFUNCT(LQNN+1,ICT)
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
C     FIRST LAYER OF LOOPS, OVER CENTERS A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     LOOP OVER CENTER A
      DO 1000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 1000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C     LOOP OVER CENTER C
      DO 1000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER C
        XYZ(1,3) = COORD(1,ICNTC)
        XYZ(2,3) = COORD(2,ICNTC)
        XYZ(3,3) = COORD(3,ICNTC)
C
C     LOOP OVER CENTER D
      DO 1000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER D
        XYZ(1,4) = COORD(1,ICNTD)
        XYZ(2,4) = COORD(2,ICNTD)
        XYZ(3,4) = COORD(3,ICNTD)
C
C     NUMBER OF NUCLEAR CENTERS INVOLVED IN THIS OVERLAP
      MCNT = NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
C
C     SKIP ONE-CENTER CONTRIBUTIONS (DEFER TO BREIT1)
C      IF(MCNT.EQ.1) THEN
C        GOTO 1001
C      ENDIF
C
C     SKIP MULTI-CENTER CONTRIBUTIONS IN STAGE 1
      IF(MCNT.NE.1.AND.ILEV.EQ.1) GOTO 1001
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
      IEAB  = 1
      IABLS = IADILS(ICNTA,ICNTB,KA,KB,MA,MB)
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
      IECD  = 1
      ICDLS = IADILS(ICNTC,ICNTD,KC,KD,MC,MD)
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
C     STARTING INDEX VALUES
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
C     IF(IQ1.LT.IQ2) GOTO 3001
C     IF(IQ3.LT.IQ4) GOTO 3001
C
      IF(IQL.GT.IQR) THEN
C       IQL > IQR
        ITSCF = 1
      ELSEIF(IQL.EQ.IQR) THEN
C       IQL = IQR
        ITSCF = 2
      ELSEIF(IQL.LT.IQR) THEN
C       IQL < IQR
        GOTO 3001
      ENDIF
C
C     READ IN FLAG VALUES FROM ISCF DATA BLOCK
      DO ISYM=1,11
        IFLG(ISYM) = ISCF(ISYM,ITSCF)
      ENDDO
C
C     EQ-COEFFICIENT PHASE FACTORS FOR PERMUTATION OF R-INTEGRALS
C     NB! OPPOSITE PHASE AS IN THE LL/SS CASE SEEN IN SCF
      PAB1 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PAB2 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      PCD1 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PCD2 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
C     FURTHER DEFINE STARTING ADDRESSES FOR {ABCD} BASIS OVERLAPS
      NA1 = LARGE(ICNTA,KA,2*MA-1)
      NA2 = LARGE(ICNTA,KA,2*MA  )
      NB1 = LARGE(ICNTB,KB,2*MB-1) + NSHIFT
      NB2 = LARGE(ICNTB,KB,2*MB  ) + NSHIFT
      NC1 = LARGE(ICNTC,KC,2*MC-1)
      NC2 = LARGE(ICNTC,KC,2*MC  )
      ND1 = LARGE(ICNTD,KD,2*MD-1) + NSHIFT
      ND2 = LARGE(ICNTD,KD,2*MD  ) + NSHIFT
C
C     UPDATE COUNTER FOR NUMBER OF BLOCKS
      N2EB(MCNT,ITT) = N2EB(MCNT,ITT) + 1
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS (IBAS,JBAS) TO CONSTRUCT BMAT/WMAT     C
C**********************************************************************C
C
C     START OF PARALLEL REGION
C!$OMP PARALLEL DO COLLAPSE(2)
C!$OMP&  PRIVATE(RR,ISKP,IFLG)
C!$OMP&  SHARED(XYZ,KQN,MQN,EXPT,NBAS,ITN)
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
          DO M=1,NBAS(3)*NBAS(4)
            IBMAP(M) = M
            IBSCR(M) = 1
          ENDDO
C
C         SCHWARZ SCREENING (ECONOMIC ONLY WHEN SCREENING FRACTION BIG)
          CALL CPU_TIME(T1)
          CALL SCHWARZ(BDSC,NBAS,ISKP,IBMAP,IBSCR)
          CALL CPU_TIME(T2)
          TB2S = TB2S + T2 - T1
C
C         CONDITIONAL TO SKIP THIS BATCH
          IF(ISKP.EQ.0) THEN
C
C           GENERATE A BATCH OF BREIT INTERACTION INTEGRALS
            CALL BII(RR,XYZ,KQN,MQN,EXPT,NBAS,IBAS,JBAS)
C
C           MULTIPLY BY DENSITY ELEMENTS AND ADD TO BMAT/WMAT
            CALL BRTMAT(RR,IFLG,NBAS,TBMC)
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
C
        ENDDO
      ENDDO
C     END OF PARALLEL REGION
C!$OMP END PARALLEL DO
C
3001  CONTINUE
3000  CONTINUE
2000  CONTINUE
1001  CONTINUE
1000  CONTINUE
C
C**********************************************************************C
C     COMPLETE CONSTRUCTION OF BREIT MATRIX BY MATRIX CONJUGATION.     C
C**********************************************************************C
C
C     NOTE: THIS WILL NEED TO BE UPDATED FOR (LL),(SS) BLOCKS
      DO I=1,NSHIFT
        DO J=NSHIFT+1,NDIM
          BDIR(J,I) = DCONJG(BDIR(I,J))
          BXCH(J,I) = DCONJG(BXCH(I,J))
        ENDDO
      ENDDO
C
C     MULTIPLY OPEN MATRIX BY ANGULAR COEFFICIENTS (LIFTED FROM COULOMB)
C      DO J=1,NDIM
C        DO I=1,NDIM
C          WDIR(I,J) = ACFF*WDIR(I,J)
C          WXCH(I,J) = BCFF*WXCH(I,J)
C        ENDDO
C      ENDDO
C
      RETURN
      END
