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
      DIMENSION INDEX(MCT,-MKP:MKP,2*(MKP+1)*MKP),ICNT(4)
C
      COMMON/BLOC/PAB1,PAB2,PCD1,PCD2,NA1,NB1,NC1,ND1,NA2,NB2,NC2,ND2,
     &            IBAS,JBAS,ILIN,NADDAB,NADDCD
      COMMON/DENS/DENC,DENO,DENT
      COMMON/EILS/EILSFL(MFL,24),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/IBSC/IBSCR(MB2),IBMAP(MB2)
      COMMON/GEOM/SHAPE
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
      LENIQ = ICOUNT
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
C     LENIQ IS THE NUMBER OF BLOCKS TO BE CALCULATED
      IQL = (IQ1-1)*LENIQ + IQ2
      IQR = (IQ3-1)*LENIQ + IQ4
C

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
C
      NC1 = LARGE(ICNTC,KC,2*MC-1)
      NC2 = LARGE(ICNTC,KC,2*MC  )
      ND1 = LARGE(ICNTD,KD,2*MD-1) + NSHIFT
      ND2 = LARGE(ICNTD,KD,2*MD  ) + NSHIFT
C
C**********************************************************************C
C     SKIP BATCHES WHICH CONFORM TO ERI PERMUTATION SYMMETRIES         C
C**********************************************************************C
C
      IF(IQL.LT.IQR) GOTO 3001
C
C     UPDATE COUNTER FOR NUMBER OF BLOCKS
      N2EB(MCNT,ITT) = N2EB(MCNT,ITT) + 1
C
C**********************************************************************C
C     THIRD LAYER OF LOOPS, OVER BASIS FUNCTIONS A AND B (USE 4000)    C
C**********************************************************************C
C
      DO 5000 IBAS=1,NBAS(1)
      DO 5000 JBAS=1,NBAS(2)
C
C**********************************************************************C
C     SKIP BATCHES THAT PASS THE SCHWARZ SCREENING PROCEDURE           C
C**********************************************************************C
C
C       RECORD CPU TIME
        CALL CPU_TIME(TBCH1)
C
C       UPDATE COUNTER FOR NUMBER OF BATCHES
        N2ET(MCNT,ITT) = N2ET(MCNT,ITT) + 1
C
C       UPDATE TIME SPENT SCREENING
        CALL CPU_TIME(TX1)
        CALL CPU_TIME(TX2)
        TSCR = TSCR + TX2 - TX1
C
C       OVERRIDE SCREENING TESTS FOR THE SAKE OF CERTAINTY
        DO M=1,NBAS(3)*NBAS(4)
          IBMAP(M) = M
          IBSCR(M) = 1
        ENDDO
C
C       GENERATE A BATCH OF BREIT INTERACTION INTEGRALS
        CALL BII(RR,XYZ,KQN,MQN,EXPT,NBAS,IBAS,JBAS)

C       1ST CASE (DIRECT):  ( MA, MB| MC, MD)
        IF(NOPN.EQ.0) GOTO 5100
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BDIR(NA1+IBAS,NB1+JBAS) = BDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 1)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 2)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 3)*DENC(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 4)*DENC(NC2+KBAS,ND2+LBAS)
C
            BDIR(NA1+IBAS,NB2+JBAS) = BDIR(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 5)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 6)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 7)*DENC(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 8)*DENC(NC2+KBAS,ND2+LBAS)
C
            BDIR(NA2+IBAS,NB1+JBAS) = BDIR(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M, 9)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,10)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,11)*DENC(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M,12)*DENC(NC2+KBAS,ND2+LBAS)
C
            BDIR(NA2+IBAS,NB2+JBAS) = BDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(M,13)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,14)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,15)*DENC(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M,16)*DENC(NC2+KBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
C
C       2ND CASE (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BDIR(NA1+IBAS,NB1+JBAS) = BDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 2)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 3)*DENC(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENC(ND2+LBAS,NC2+KBAS)
C
            BDIR(NA1+IBAS,NB2+JBAS) = BDIR(NA1+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 8)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 7)*DENC(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 5)*DENC(ND2+LBAS,NC2+KBAS)
C
            BDIR(NA2+IBAS,NB1+JBAS) = BDIR(NA2+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,12)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENC(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENC(ND2+LBAS,NC2+KBAS)
C
            BDIR(NA2+IBAS,NB2+JBAS) = BDIR(NA2+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(M,16)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,14)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENC(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENC(ND2+LBAS,NC2+KBAS)
C
          ENDDO
        ENDDO
C
C ***   THERE IS AN EXTRA TERM FOR ELEMENTS SATISFYING IQL = IQR
        IF(IQL.EQ.IQR) GOTO 5200

C       3RD CASE (DIRECT):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BDIR(NC1+KBAS,ND1+LBAS) = BDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 1)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 5)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 9)*DENC(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,13)*DENC(NA2+IBAS,NB2+JBAS)
C
            BDIR(NC1+KBAS,ND2+LBAS) = BDIR(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 2)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 6)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,10)*DENC(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,14)*DENC(NA2+IBAS,NB2+JBAS)
C
            BDIR(NC2+KBAS,ND1+LBAS) = BDIR(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 3)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 7)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,11)*DENC(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,15)*DENC(NA2+IBAS,NB2+JBAS)
C
            BDIR(NC2+KBAS,ND2+LBAS) = BDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(M, 4)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 8)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,12)*DENC(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,16)*DENC(NA2+IBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
C
C       4TH CASE (DIRECT):  ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                             = PAB*    (-MA,-MB| MC, MD)
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BDIR(NC1+KBAS,ND1+LBAS) = BDIR(NC1+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,13)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 5)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 9)*DENC(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENC(NB2+JBAS,NA2+IBAS)
C
            BDIR(NC1+KBAS,ND2+LBAS) = BDIR(NC1+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,14)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,10)*DENC(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 2)*DENC(NB2+JBAS,NA2+IBAS)
C
            BDIR(NC2+KBAS,ND1+LBAS) = BDIR(NC2+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,15)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENC(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENC(NB2+JBAS,NA2+IBAS)
C
            BDIR(NC2+KBAS,ND2+LBAS) = BDIR(NC2+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,16)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 8)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENC(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENC(NB2+JBAS,NA2+IBAS)
C
          ENDDO
        ENDDO
C
5200    CONTINUE
C
5100    CONTINUE
C
C       5TH CASE (EXCHNG):  ( MA, MB| MC, MD)
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BXCH(NA1+IBAS,ND1+LBAS) = BXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 1)*DENC(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 5)*DENC(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 3)*DENC(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M, 7)*DENC(NC2+KBAS,NB2+JBAS)
C
            BXCH(NA1+IBAS,ND2+LBAS) = BXCH(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M, 2)*DENC(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 6)*DENC(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 4)*DENC(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M, 8)*DENC(NC2+KBAS,NB2+JBAS)
C
            BXCH(NA2+IBAS,ND1+LBAS) = BXCH(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M, 9)*DENC(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M,13)*DENC(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M,11)*DENC(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M,15)*DENC(NC2+KBAS,NB2+JBAS)
C
            BXCH(NA2+IBAS,ND2+LBAS) = BXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(M,10)*DENC(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M,14)*DENC(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M,12)*DENC(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M,16)*DENC(NC2+KBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
C
C       6TH CASE (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BXCH(NA1+IBAS,NC1+KBAS) = BXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 4)*DENC(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 8)*DENC(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M, 3)*DENC(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 7)*DENC(ND2+LBAS,NB2+JBAS)
C
            BXCH(NA1+IBAS,NC2+KBAS) = BXCH(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 2)*DENC(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 6)*DENC(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 1)*DENC(ND2+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 5)*DENC(ND2+LBAS,NB2+JBAS)
C
            BXCH(NA2+IBAS,NC1+KBAS) = BXCH(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,12)*DENC(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,16)*DENC(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M,11)*DENC(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M,15)*DENC(ND2+LBAS,NB2+JBAS)
C
            BXCH(NA2+IBAS,NC2+KBAS) = BXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,10)*DENC(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M,14)*DENC(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 9)*DENC(ND2+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,13)*DENC(ND2+LBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
C
C       7TH CASE (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BXCH(NB1+JBAS,ND1+LBAS) = BXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,13)*DENC(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 5)*DENC(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(M,15)*DENC(NC2+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENC(NC2+KBAS,NA2+IBAS)
C
            BXCH(NB1+JBAS,ND2+LBAS) = BXCH(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,14)*DENC(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENC(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(M,16)*DENC(NC2+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 8)*DENC(NC2+KBAS,NA2+IBAS)
C
            BXCH(NB2+JBAS,ND1+LBAS) = BXCH(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 9)*DENC(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENC(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENC(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENC(NC2+KBAS,NA2+IBAS)
C
            BXCH(NB2+JBAS,ND2+LBAS) = BXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M,10)*DENC(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 2)*DENC(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENC(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENC(NC2+KBAS,NA2+IBAS)
C
          ENDDO
        ENDDO
C
C     8TH CASE (EXCHNG):  ( MB, MA| MD, MC) = PAB*    (-MA,-MB| MC, MD)
C                                           = PAB*PCD*(-MA,-MB|-MC,-MD)
C     *** BLANK ***
C
C       THERE IS AN EXTRA TERM FOR ELEMENTS SATISFYING IQL = IQR
        IF(IQL.EQ.IQR) GOTO 5300
C
C       9TH CASE (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BXCH(NC1+KBAS,NB1+JBAS) = BXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 1)*DENC(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 2)*DENC(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M, 9)*DENC(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M,10)*DENC(NA2+IBAS,ND2+LBAS)
C
            BXCH(NC1+KBAS,NB2+JBAS) = BXCH(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 5)*DENC(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 6)*DENC(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,13)*DENC(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M,14)*DENC(NA2+IBAS,ND2+LBAS)
C
            BXCH(NC2+KBAS,NB1+JBAS) = BXCH(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M, 3)*DENC(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 4)*DENC(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,11)*DENC(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M,12)*DENC(NA2+IBAS,ND2+LBAS)
C
            BXCH(NC2+KBAS,NB2+JBAS) = BXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(M, 7)*DENC(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 8)*DENC(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,15)*DENC(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M,16)*DENC(NA2+IBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
C
C       10TH CASE (EXCHNG): ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                             = PAB*    (-MA,-MB| MC, MD)
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BXCH(NC1+KBAS,NA1+IBAS) = BXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,13)*DENC(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,14)*DENC(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M, 9)*DENC(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M,10)*DENC(NB2+JBAS,ND2+LBAS)
C
            BXCH(NC1+KBAS,NA2+IBAS) = BXCH(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 5)*DENC(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 6)*DENC(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M, 1)*DENC(NB2+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M, 2)*DENC(NB2+JBAS,ND2+LBAS)
C
            BXCH(NC2+KBAS,NA1+IBAS) = BXCH(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,15)*DENC(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,16)*DENC(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M,11)*DENC(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M,12)*DENC(NB2+JBAS,ND2+LBAS)
C
            BXCH(NC2+KBAS,NA2+IBAS) = BXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 7)*DENC(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 8)*DENC(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M, 3)*DENC(NB2+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M, 4)*DENC(NB2+JBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
C
C       11TH CASE (EXCHNG): ( MD, MC| MA, MB) =         ( MA, MB| MD, MC)
C                                             =     PCD*( MA, MB|-MC,-MD)
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BXCH(ND1+LBAS,NB1+JBAS) = BXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENC(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 2)*DENC(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(M,12)*DENC(NA2+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENC(NA2+IBAS,NC2+KBAS)
C
            BXCH(ND1+LBAS,NB2+JBAS) = BXCH(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 8)*DENC(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENC(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(M,16)*DENC(NA2+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,14)*DENC(NA2+IBAS,NC2+KBAS)
C
            BXCH(ND2+LBAS,NB1+JBAS) = BXCH(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 3)*DENC(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENC(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENC(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENC(NA2+IBAS,NC2+KBAS)
C
            BXCH(ND2+LBAS,NB2+JBAS) = BXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M, 7)*DENC(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 5)*DENC(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENC(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENC(NA2+IBAS,NC2+KBAS)
C
          ENDDO
        ENDDO
C
5300   CONTINUE
C
C**********************************************************************C
C     END OF BREIT MATRIX CONSTRUCTION                                 C
C**********************************************************************C
C
C     RECORD CPU TIME AT END OF BATCH AND ADD TO TIME COUNTER
      CALL CPU_TIME(TBCH2)
      T2ES(MCNT,ITT) = T2ES(MCNT,ITT) + TBCH2 - TBCH1
C
5000  CONTINUE
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
      DO I=1,NSHIFT
        DO J=NSHIFT+1,NDIM
          BDIR(J,I) = DCONJG(BDIR(I,J))
          BXCH(J,I) = DCONJG(BXCH(I,J))
        ENDDO
      ENDDO
C
      RETURN
      END
