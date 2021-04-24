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
      DIMENSION IBINDX(MCT,MKP)
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
C**********************************************************************C
C     ORDERED INDEX OF (ICNT,KQN,MQN) COMBINATIONS                     C
C**********************************************************************C
C
      ICOUNT = 0
      IKOUNT = 0
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
C         RECORD CUMULATIVE INDEX COUNTER FOR NUMBER OF BASIS FUNCTIONS
          IBINDX(ICT,KN) = IKOUNT
          IKOUNT         = IKOUNT+NFUN
C
C         LOOP OVER MQN VALUES AND RECORD INDEX
          DO MJ=1,MJMAX,2
            ICOUNT              = ICOUNT+1
            INDEX(ICT,KAPPA,MJ) = ICOUNT
            LENIQ               = ICOUNT
          ENDDO
c
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
C     LOOP OVER KQN(A) VALUES
      DO 1000 KA=1,NKAP(ICNTA)
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
      DO 1000 KB=1,NKAP(ICNTB)
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
C     LOOP OVER |MQN(A)| VALUES
      DO 1000 MA=1,IABS(KQN(1))
        MJA    = 2*MA-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 1000 MB=1,IABS(KQN(2))
        MJB    = 2*MB-1
        MQN(2) = MJB
C
C     CALCULATE NEW BLOCK OF E(AB) COEFFS AT NEXT OPPORTUNITY
      IEAB  = 1
      IABLS = IADILS(ICNTA,ICNTB,KA,KB,MA,MB)
C
C     TWO-ELECTRON COMPONENT OVERLAP INDEX
      ITT = 5
C
C**********************************************************************C
C     SECOND LAYER OF LOOPS, OVER CENTERS C AND D (USE INDEX 3000)     C
C**********************************************************************C
C
C     LOOP OVER CENTER C
      DO 2000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER C
        XYZ(1,3) = COORD(1,ICNTC)
        XYZ(2,3) = COORD(2,ICNTC)
        XYZ(3,3) = COORD(3,ICNTC)
C
C     LOOP OVER CENTER D
      DO 2000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER D
        XYZ(1,4) = COORD(1,ICNTD)
        XYZ(2,4) = COORD(2,ICNTD)
        XYZ(3,4) = COORD(3,ICNTD)
C
C       NUMBER OF NUCLEAR CENTERS INVOLVED IN THIS OVERLAP
        MCNT = NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
C
C       SKIP MULTI-CENTER CONTRIBUTIONS IN STAGE 1
        IF(MCNT.NE.1.AND.ILEV.EQ.1) GOTO 2200
C
C     LOOP OVER KQN(C) VALUES
      DO 2500 KC=1,NKAP(ICNTC)
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
      DO 2500 KD=1,NKAP(ICNTD)
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
C     LOOP OVER |MQN(C)| VALUES
      DO 2500 MC=1,IABS(KQN(3))
        MJC    = 2*MC-1
        MQN(3) = MJC
C
C     LOOP OVER |MQN(D)| VALUES
      DO 2500 MD=1,IABS(KQN(4))
        MJD    = 2*MD-1
        MQN(4) = MJD
C
C     CALCULATE NEW BLOCK OF E(CD) COEFFS AT NEXT OPPORTUNITY
      IECD  = 1
      ICDLS = IADILS(ICNTC,ICNTD,KC,KD,MC,MD)
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
        GOTO 2999
      ENDIF
5503  CONTINUE
C
C**********************************************************************C
C     FOR THIS CHOICE OF A,B,C AND D, COMPUTE INDICES                  C
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
C     FURTHER DEFINE STARTING ADDRESSES FOR {ABCD} BASIS OVERLAPS
      IA1 = LARGE(ICNTA,KA,MJA  )
      IB1 = LARGE(ICNTB,KB,MJB  ) + NSHIFT
      IC1 = LARGE(ICNTC,KC,MJC  )
      ID1 = LARGE(ICNTD,KD,MJD  ) + NSHIFT
C
      IA2 = LARGE(ICNTA,KA,MJA+1)
      IB2 = LARGE(ICNTB,KB,MJB+1) + NSHIFT
      IC2 = LARGE(ICNTC,KC,MJC+1)
      ID2 = LARGE(ICNTD,KD,MJD+1) + NSHIFT
C
      JA1 = LARGE(ICNTA,KA,MJA  )
      JB1 = LARGE(ICNTB,KB,MJB  ) + NSHIFT
      JC1 = LARGE(ICNTC,KC,MJC  )
      JD1 = LARGE(ICNTD,KD,MJD  ) + NSHIFT
C
      JA2 = LARGE(ICNTA,KA,MJA+1)
      JB2 = LARGE(ICNTB,KB,MJB+1) + NSHIFT
      JC2 = LARGE(ICNTC,KC,MJC+1)
      JD2 = LARGE(ICNTD,KD,MJD+1) + NSHIFT
C
C     NOT SURE WHAT THESE ARE FOR
      IAA = IBINDX(ICNTA,KA)
      JBB = IBINDX(ICNTB,KB) + NSHIFT
C
      ICC = IBINDX(ICNTC,KC)
      JDD = IBINDX(ICNTD,KD) + NSHIFT
C
C**********************************************************************C
C     SKIP BATCHES WHICH CONFORM TO ERI PERMUTATION SYMMETRIES         C
C**********************************************************************C
C
      IF(IQL.LT.IQR) GOTO 2999
C
C
C     UPDATE COUNTER FOR NUMBER OF BLOCKS
      N2EB(MCNT,ITT) = N2EB(MCNT,ITT) + 1
C
C
C**********************************************************************C
C     FOR THIS CHOICE OF A,B,C AND D, COMPUTE PHASES                   C
C**********************************************************************C
C
C     CALCULATE KQN PHASE FACTORS FOR PERMUTING INTEGRALS.
C     NB! OPPOSITE PHASE AS IN THE LL/SS CASE SEEN IN SCF
      IF((KQN(1)*KQN(2)).GT.0) THEN
        PKAB =-1.0D0
      ELSE
        PKAB = 1.0D0
      ENDIF
C
      IF((KQN(3)*KQN(4)).GT.0) THEN
        PKCD =-1.0D0
      ELSE
        PKCD = 1.0D0
      ENDIF
C
C     CALCULATE MQN PHASE FACTORS FOR PERMUTING INTEGRALS
      PMAB1 = DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PMAB2 = DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      PMCD1 = DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PMCD2 = DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
      F1 = PKAB*PMAB1
      F2 = PKAB*PMAB2
      G1 = PKCD*PMCD1
      G2 = PKCD*PMCD2
C
C**********************************************************************C
C     THIRD LAYER OF LOOPS, OVER BASIS FUNCTIONS A AND B (USE 4000)    C
C**********************************************************************C
C
      DO 3000 IBAS=1,NBAS(1)
      DO 3000 JBAS=1,NBAS(2)
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
C
C       FIRST BATCH
        IF(NOPN.EQ.0) GOTO 301
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BDIR(IA1+IBAS,JB1+JBAS) = BDIR(IA1+IBAS,JB1+JBAS)
     &                            +    RR(M, 1)*DENT(IC1+KBAS,JD1+LBAS)
     &                            +    RR(M, 2)*DENT(IC1+KBAS,JD2+LBAS)
     &                            +    RR(M, 3)*DENT(IC2+KBAS,JD1+LBAS)
     &                            +    RR(M, 4)*DENT(IC2+KBAS,JD2+LBAS)
     &                            + G1*RR(M, 4)*DENT(JD1+LBAS,IC1+KBAS)
     &                            + G2*RR(M, 2)*DENT(JD1+LBAS,IC2+KBAS)
     &                            + G2*RR(M, 3)*DENT(JD2+LBAS,IC1+KBAS)
     &                            + G1*RR(M, 1)*DENT(JD2+LBAS,IC2+KBAS)
C
            BDIR(IA1+IBAS,JB2+JBAS) = BDIR(IA1+IBAS,JB2+JBAS)
     &                            +    RR(M, 5)*DENT(IC1+KBAS,JD1+LBAS)
     &                            +    RR(M, 6)*DENT(IC1+KBAS,JD2+LBAS)
     &                            +    RR(M, 7)*DENT(IC2+KBAS,JD1+LBAS)
     &                            +    RR(M, 8)*DENT(IC2+KBAS,JD2+LBAS)
     &                            + G1*RR(M, 8)*DENT(JD1+LBAS,IC1+KBAS)
     &                            + G2*RR(M, 6)*DENT(JD1+LBAS,IC2+KBAS)
     &                            + G2*RR(M, 7)*DENT(JD2+LBAS,IC1+KBAS)
     &                            + G1*RR(M, 5)*DENT(JD2+LBAS,IC2+KBAS)
C
            BDIR(IA2+IBAS,JB1+JBAS) = BDIR(IA2+IBAS,JB1+JBAS)
     &                            +    RR(M, 9)*DENT(IC1+KBAS,JD1+LBAS)
     &                            +    RR(M,10)*DENT(IC1+KBAS,JD2+LBAS)
     &                            +    RR(M,11)*DENT(IC2+KBAS,JD1+LBAS)
     &                            +    RR(M,12)*DENT(IC2+KBAS,JD2+LBAS)
     &                            + G1*RR(M,12)*DENT(JD1+LBAS,IC1+KBAS)
     &                            + G2*RR(M,10)*DENT(JD1+LBAS,IC2+KBAS)
     &                            + G2*RR(M,11)*DENT(JD2+LBAS,IC1+KBAS)
     &                            + G1*RR(M, 9)*DENT(JD2+LBAS,IC2+KBAS)
C
            BDIR(IA2+IBAS,JB2+JBAS) = BDIR(IA2+IBAS,JB2+JBAS)
     &                            +    RR(M,13)*DENT(IC1+KBAS,JD1+LBAS)
     &                            +    RR(M,14)*DENT(IC1+KBAS,JD2+LBAS)
     &                            +    RR(M,15)*DENT(IC2+KBAS,JD1+LBAS)
     &                            +    RR(M,16)*DENT(IC2+KBAS,JD2+LBAS)
     &                            + G1*RR(M,16)*DENT(JD1+LBAS,IC1+KBAS)
     &                            + G2*RR(M,14)*DENT(JD1+LBAS,IC2+KBAS)
     &                            + G2*RR(M,15)*DENT(JD2+LBAS,IC1+KBAS)
     &                            + G1*RR(M,13)*DENT(JD2+LBAS,IC2+KBAS)
C
          ENDDO
        ENDDO
301     CONTINUE
C
C       SECOND BATCH
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BXCH(IA1+IBAS,JD1+LBAS) = BXCH(IA1+IBAS,JD1+LBAS)
     &                            +    RR(M, 1)*DENT(IC1+KBAS,JB1+JBAS)
     &                            +    RR(M, 5)*DENT(IC1+KBAS,JB2+JBAS)
     &                            +    RR(M, 3)*DENT(IC2+KBAS,JB1+JBAS)
     &                            +    RR(M, 7)*DENT(IC2+KBAS,JB2+JBAS)
C
            BXCH(IA1+IBAS,JD2+LBAS) = BXCH(IA1+IBAS,JD2+LBAS)
     &                            +    RR(M, 2)*DENT(IC1+KBAS,JB1+JBAS)
     &                            +    RR(M, 6)*DENT(IC1+KBAS,JB2+JBAS)
     &                            +    RR(M, 4)*DENT(IC2+KBAS,JB1+JBAS)
     &                            +    RR(M, 8)*DENT(IC2+KBAS,JB2+JBAS)
C
            BXCH(IA2+IBAS,JD1+LBAS) = BXCH(IA2+IBAS,JD1+LBAS)
     &                            +    RR(M, 9)*DENT(IC1+KBAS,JB1+JBAS)
     &                            +    RR(M,13)*DENT(IC1+KBAS,JB2+JBAS)
     &                            +    RR(M,11)*DENT(IC2+KBAS,JB1+JBAS)
     &                            +    RR(M,15)*DENT(IC2+KBAS,JB2+JBAS)
C
            BXCH(IA2+IBAS,JD2+LBAS) = BXCH(IA2+IBAS,JD2+LBAS)
     &                            +    RR(M,10)*DENT(IC1+KBAS,JB1+JBAS)
     &                            +    RR(M,14)*DENT(IC1+KBAS,JB2+JBAS)
     &                            +    RR(M,12)*DENT(IC2+KBAS,JB1+JBAS)
     &                            +    RR(M,16)*DENT(IC2+KBAS,JB2+JBAS)
C
          ENDDO
        ENDDO
C
C       THIRD BATCH
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BXCH(IA1+IBAS,JC1+KBAS) = BXCH(IA1+IBAS,JC1+KBAS)
     &                            + G1*RR(M, 4)*DENT(ID1+LBAS,JB1+JBAS)
     &                            + G1*RR(M, 8)*DENT(ID1+LBAS,JB2+JBAS)
     &                            + G2*RR(M, 3)*DENT(ID2+LBAS,JB1+JBAS)
     &                            + G2*RR(M, 7)*DENT(ID2+LBAS,JB2+JBAS)
C
            BXCH(IA1+IBAS,JC2+KBAS) = BXCH(IA1+IBAS,JC2+KBAS)
     &                            + G2*RR(M, 2)*DENT(ID1+LBAS,JB1+JBAS)
     &                            + G2*RR(M, 6)*DENT(ID1+LBAS,JB2+JBAS)
     &                            + G1*RR(M, 1)*DENT(ID2+LBAS,JB1+JBAS)
     &                            + G1*RR(M, 5)*DENT(ID2+LBAS,JB2+JBAS)
C
            BXCH(IA2+IBAS,JC1+KBAS) = BXCH(IA2+IBAS,JC1+KBAS)
     &                            + G1*RR(M,12)*DENT(ID1+LBAS,JB1+JBAS)
     &                            + G1*RR(M,16)*DENT(ID1+LBAS,JB2+JBAS)
     &                            + G2*RR(M,11)*DENT(ID2+LBAS,JB1+JBAS)
     &                            + G2*RR(M,15)*DENT(ID2+LBAS,JB2+JBAS)
C
            BXCH(IA2+IBAS,JC2+KBAS) = BXCH(IA2+IBAS,JC2+KBAS)
     &                            + G2*RR(M,10)*DENT(ID1+LBAS,JB1+JBAS)
     &                            + G2*RR(M,14)*DENT(ID1+LBAS,JB2+JBAS)
     &                            + G1*RR(M, 9)*DENT(ID2+LBAS,JB1+JBAS)
     &                            + G1*RR(M,13)*DENT(ID2+LBAS,JB2+JBAS)
C
          ENDDO
        ENDDO
C
C       FOURTH BATCH
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BXCH(IB1+JBAS,JD1+LBAS) = BXCH(IB1+JBAS,JD1+LBAS)
     &                            + F1*RR(M,13)*DENT(IC1+KBAS,JA1+IBAS)
     &                            + F2*RR(M, 5)*DENT(IC1+KBAS,JA2+IBAS)
     &                            + F1*RR(M,15)*DENT(IC2+KBAS,JA1+IBAS)
     &                            + F2*RR(M, 7)*DENT(IC2+KBAS,JA2+IBAS)
C
            BXCH(IB1+JBAS,JD2+LBAS) = BXCH(IB1+JBAS,JD2+LBAS)
     &                            + F1*RR(M,14)*DENT(IC1+KBAS,JA1+IBAS)
     &                            + F2*RR(M, 6)*DENT(IC1+KBAS,JA2+IBAS)
     &                            + F1*RR(M,16)*DENT(IC2+KBAS,JA1+IBAS)
     &                            + F2*RR(M, 8)*DENT(IC2+KBAS,JA2+IBAS)
C
            BXCH(IB2+JBAS,JD1+LBAS) = BXCH(IB2+JBAS,JD1+LBAS)
     &                            + F2*RR(M, 9)*DENT(IC1+KBAS,JA1+IBAS)
     &                            + F1*RR(M, 1)*DENT(IC1+KBAS,JA2+IBAS)
     &                            + F2*RR(M,11)*DENT(IC2+KBAS,JA1+IBAS)
     &                            + F1*RR(M, 3)*DENT(IC2+KBAS,JA2+IBAS)
C
            BXCH(IB2+JBAS,JD2+LBAS) = BXCH(IB2+JBAS,JD2+LBAS)
     &                            + F2*RR(M,10)*DENT(IC1+KBAS,JA1+IBAS)
     &                            + F1*RR(M, 2)*DENT(IC1+KBAS,JA2+IBAS)
     &                            + F2*RR(M,12)*DENT(IC2+KBAS,JA1+IBAS)
     &                            + F1*RR(M, 4)*DENT(IC2+KBAS,JA2+IBAS)
C
          ENDDO
        ENDDO
C
C ***   THERE IS AN EXTRA TERM FOR ELEMENTS SATISFYING IQL = IQR
        IF(IQL.EQ.IQR) GOTO 100
C
C       FIFTH BATCH
        IF(NOPN.EQ.0) GOTO 302
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BDIR(IC1+KBAS,JD1+LBAS) = BDIR(IC1+KBAS,JD1+LBAS)
     &                            +    RR(M, 1)*DENT(IA1+IBAS,JB1+JBAS)
     &                            +    RR(M, 5)*DENT(IA1+IBAS,JB2+JBAS)
     &                            +    RR(M, 9)*DENT(IA2+IBAS,JB1+JBAS)
     &                            +    RR(M,13)*DENT(IA2+IBAS,JB2+JBAS)
     &                            + F1*RR(M,13)*DENT(JB1+JBAS,IA1+IBAS)
     &                            + F2*RR(M, 5)*DENT(JB1+JBAS,IA2+IBAS)
     &                            + F2*RR(M, 9)*DENT(JB2+JBAS,IA1+IBAS)
     &                            + F1*RR(M, 1)*DENT(JB2+JBAS,IA2+IBAS)
C
            BDIR(IC1+KBAS,JD2+LBAS) = BDIR(IC1+KBAS,JD2+LBAS)
     &                            +    RR(M, 2)*DENT(IA1+IBAS,JB1+JBAS)
     &                            +    RR(M, 6)*DENT(IA1+IBAS,JB2+JBAS)
     &                            +    RR(M,10)*DENT(IA2+IBAS,JB1+JBAS)
     &                            +    RR(M,14)*DENT(IA2+IBAS,JB2+JBAS)
     &                            + F1*RR(M,14)*DENT(JB1+JBAS,IA1+IBAS)
     &                            + F2*RR(M, 6)*DENT(JB1+JBAS,IA2+IBAS)
     &                            + F2*RR(M,10)*DENT(JB2+JBAS,IA1+IBAS)
     &                            + F1*RR(M, 2)*DENT(JB2+JBAS,IA2+IBAS)
C
            BDIR(IC2+KBAS,JD1+LBAS) = BDIR(IC2+KBAS,JD1+LBAS)
     &                            +    RR(M, 3)*DENT(IA1+IBAS,JB1+JBAS)
     &                            +    RR(M, 7)*DENT(IA1+IBAS,JB2+JBAS)
     &                            +    RR(M,11)*DENT(IA2+IBAS,JB1+JBAS)
     &                            +    RR(M,15)*DENT(IA2+IBAS,JB2+JBAS)
     &                            + F1*RR(M,15)*DENT(JB1+JBAS,IA1+IBAS)
     &                            + F2*RR(M, 7)*DENT(JB1+JBAS,IA2+IBAS)
     &                            + F2*RR(M,11)*DENT(JB2+JBAS,IA1+IBAS)
     &                            + F1*RR(M, 3)*DENT(JB2+JBAS,IA2+IBAS)
C
            BDIR(IC2+KBAS,JD2+LBAS) = BDIR(IC2+KBAS,JD2+LBAS)
     &                            +    RR(M, 4)*DENT(IA1+IBAS,JB1+JBAS)
     &                            +    RR(M, 8)*DENT(IA1+IBAS,JB2+JBAS)
     &                            +    RR(M,12)*DENT(IA2+IBAS,JB1+JBAS)
     &                            +    RR(M,16)*DENT(IA2+IBAS,JB2+JBAS)
     &                            + F1*RR(M,16)*DENT(JB1+JBAS,IA1+IBAS)
     &                            + F2*RR(M, 8)*DENT(JB1+JBAS,IA2+IBAS)
     &                            + F2*RR(M,12)*DENT(JB2+JBAS,IA1+IBAS)
     &                            + F1*RR(M, 4)*DENT(JB2+JBAS,IA2+IBAS)
C
          ENDDO
        ENDDO
302     CONTINUE
C
C       SIXTH BATCH
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BXCH(IC1+KBAS,JA1+IBAS) = BXCH(IC1+KBAS,JA1+IBAS)
     &                            + F1*RR(M,13)*DENT(JB1+JBAS,ID1+LBAS)
     &                            + F1*RR(M,14)*DENT(JB1+JBAS,ID2+LBAS)
     &                            + F2*RR(M, 9)*DENT(JB2+JBAS,ID1+LBAS)
     &                            + F2*RR(M,10)*DENT(JB2+JBAS,ID2+LBAS)
C
            BXCH(IC1+KBAS,JA2+IBAS) = BXCH(IC1+KBAS,JA2+IBAS)
     &                            + F2*RR(M, 5)*DENT(JB1+JBAS,ID1+LBAS)
     &                            + F2*RR(M, 6)*DENT(JB1+JBAS,ID2+LBAS)
     &                            + F1*RR(M, 1)*DENT(JB2+JBAS,ID1+LBAS)
     &                            + F1*RR(M, 2)*DENT(JB2+JBAS,ID2+LBAS)
C
            BXCH(IC2+KBAS,JA1+IBAS) = BXCH(IC2+KBAS,JA1+IBAS)
     &                            + F1*RR(M,15)*DENT(JB1+JBAS,ID1+LBAS)
     &                            + F1*RR(M,16)*DENT(JB1+JBAS,ID2+LBAS)
     &                            + F2*RR(M,11)*DENT(JB2+JBAS,ID1+LBAS)
     &                            + F2*RR(M,12)*DENT(JB2+JBAS,ID2+LBAS)
C
            BXCH(IC2+KBAS,JA2+IBAS) = BXCH(IC2+KBAS,JA2+IBAS)
     &                            + F2*RR(M, 7)*DENT(JB1+JBAS,ID1+LBAS)
     &                            + F2*RR(M, 8)*DENT(JB1+JBAS,ID2+LBAS)
     &                            + F1*RR(M, 3)*DENT(JB2+JBAS,ID1+LBAS)
     &                            + F1*RR(M, 4)*DENT(JB2+JBAS,ID2+LBAS)
C
          ENDDO
        ENDDO
C
C       SEVENTH BATCH
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BXCH(ID1+LBAS,JB1+JBAS) = BXCH(ID1+LBAS,JB1+JBAS)
     &                            + G1*RR(M, 4)*DENT(JA1+IBAS,IC1+KBAS)
     &                            + G2*RR(M, 2)*DENT(JA1+IBAS,IC2+KBAS)
     &                            + G1*RR(M,12)*DENT(JA2+IBAS,IC1+KBAS)
     &                            + G2*RR(M,10)*DENT(JA2+IBAS,IC2+KBAS)
C
            BXCH(ID1+LBAS,JB2+JBAS) = BXCH(ID1+LBAS,JB2+JBAS)
     &                            + G1*RR(M, 8)*DENT(JA1+IBAS,IC1+KBAS)
     &                            + G2*RR(M, 6)*DENT(JA1+IBAS,IC2+KBAS)
     &                            + G1*RR(M,16)*DENT(JA2+IBAS,IC1+KBAS)
     &                            + G2*RR(M,14)*DENT(JA2+IBAS,IC2+KBAS)
C
            BXCH(ID2+LBAS,JB1+JBAS) = BXCH(ID2+LBAS,JB1+JBAS)
     &                            + G2*RR(M, 3)*DENT(JA1+IBAS,IC1+KBAS)
     &                            + G1*RR(M, 1)*DENT(JA1+IBAS,IC2+KBAS)
     &                            + G2*RR(M,11)*DENT(JA2+IBAS,IC1+KBAS)
     &                            + G1*RR(M, 9)*DENT(JA2+IBAS,IC2+KBAS)
C
            BXCH(ID2+LBAS,JB2+JBAS) = BXCH(ID2+LBAS,JB2+JBAS)
     &                            + G2*RR(M, 7)*DENT(JA1+IBAS,IC1+KBAS)
     &                            + G1*RR(M, 5)*DENT(JA1+IBAS,IC2+KBAS)
     &                            + G2*RR(M,15)*DENT(JA2+IBAS,IC1+KBAS)
     &                            + G1*RR(M,13)*DENT(JA2+IBAS,IC2+KBAS)
C
          ENDDO
        ENDDO
C
C       EIGHTH BATCH
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            BXCH(IC1+KBAS,JB1+JBAS) = BXCH(IC1+KBAS,JB1+JBAS)
     &                            +    RR(M, 1)*DENT(JA1+IBAS,ID1+LBAS)
     &                            +    RR(M, 2)*DENT(JA1+IBAS,ID2+LBAS)
     &                            +    RR(M, 9)*DENT(JA2+IBAS,ID1+LBAS)
     &                            +    RR(M,10)*DENT(JA2+IBAS,ID2+LBAS)
C
            BXCH(IC1+KBAS,JB2+JBAS) = BXCH(IC1+KBAS,JB2+JBAS)
     &                            +    RR(M, 5)*DENT(JA1+IBAS,ID1+LBAS)
     &                            +    RR(M, 6)*DENT(JA1+IBAS,ID2+LBAS)
     &                            +    RR(M,13)*DENT(JA2+IBAS,ID1+LBAS)
     &                            +    RR(M,14)*DENT(JA2+IBAS,ID2+LBAS)
C
            BXCH(IC2+KBAS,JB1+JBAS) = BXCH(IC2+KBAS,JB1+JBAS)
     &                            +    RR(M, 3)*DENT(JA1+IBAS,ID1+LBAS)
     &                            +    RR(M, 4)*DENT(JA1+IBAS,ID2+LBAS)
     &                            +    RR(M,11)*DENT(JA2+IBAS,ID1+LBAS)
     &                            +    RR(M,12)*DENT(JA2+IBAS,ID2+LBAS)
C
            BXCH(IC2+KBAS,JB2+JBAS) = BXCH(IC2+KBAS,JB2+JBAS)
     &                            +    RR(M, 7)*DENT(JA1+IBAS,ID1+LBAS)
     &                            +    RR(M, 8)*DENT(JA1+IBAS,ID2+LBAS)
     &                            +    RR(M,15)*DENT(JA2+IBAS,ID1+LBAS)
     &                            +    RR(M,16)*DENT(JA2+IBAS,ID2+LBAS)
C
          ENDDO
        ENDDO
C
100     CONTINUE
C
C**********************************************************************C
C     END OF BREIT MATRIX CONSTRUCTION                                 C
C**********************************************************************C
C
C     RECORD CPU TIME AT END OF BATCH AND ADD TO TIME COUNTER
      CALL CPU_TIME(TBCH2)
      T2ES(MCNT,ITT) = T2ES(MCNT,ITT) + TBCH2 - TBCH1
C
3000  CONTINUE
2999  CONTINUE
2500  CONTINUE
2200  CONTINUE
2000  CONTINUE
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

