      SUBROUTINE COULOMBOLD
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
C  COULOMB GENERATES ELECTRON REPULSION INTEGRALS IN BATCHES AND       C
C  CALCULATES THE SCF COULOMB MATRIX (G), APPLYING IT DIRECTLY TO THE  C
C  FOCK MATRIX. IN THE CASE OF OPEN SUBSHELLS, THE Q MATRIX IS ALSO    C
C  INCLUDED. INTEGRAL SYMMETRY PARTIALLY EXPLOITED, BUT WITH ROOM FOR  C
C  IMPROVEMENT (GEOMETRIC SYMM, R-INT SYMM, E-COEFF SYMM).             C
C -------------------------------------------------------------------- C
C  DFNOTE: THIS ROUTINE COULD BENEFIT FROM PARALLELISATION -- OPENMP.  C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=4,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 T1(MDM),T2(MDM),T3(MDM)
      COMPLEX*16 RR(MB2,16)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NFUNS(4),LQN(4)
      DIMENSION ITQN(2),IFLG(11),ISCF(11,6)
      DIMENSION INDEX(MCT,-MKP:MKP,2*(MKP+1)*MKP)
      DIMENSION GDSC(MDM,MDM),GXSC(MDM,MDM),BDSC(MDM,MDM),BXSC(MDM,MDM)
C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS
      COMMON/MTRX/OVAP,HNUC,HKIN,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,IEQS,IPAR,ICOR,ILEV
      COMMON/QNMS/ICNBAS(MDM),KQNBAS(MDM),MQNBAS(MDM)
      COMMON/SCRN/N2ET(5,7),N2ES(5,7),F2ES(5,7),T2ES(5,7),TSCR
      COMMON/SHLL/ACFF,BCFF,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN,NOELEC
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/SWRZ/GDSC,GXSC,BDSC,BXSC
C
C     ISCF TELLS WHICH INTEGRALS TO INCLUDE BASED ON OVERLAP COMBINATION
      DATA ISCF/1,0,1,0,1,0,1,0,1,1,0,
     &          1,0,0,0,1,0,0,0,1,1,0,
     &          0,1,1,0,0,0,0,0,1,1,0,
     &          1,0,0,1,1,0,0,0,1,0,0,
     &          0,1,0,1,0,0,0,0,1,0,0,
     &          0,1,0,0,0,0,0,0,1,0,0/
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
C     BARE NUCLEUS APPROXIMATION: NO MEAN-FIELD COULOMB MATRIX
      IF(HMLTN.EQ.'BARE') RETURN
C
C**********************************************************************C
C     INDEXING ROUTINE: SO WE CAN SET BASIS FUNCTION LABELS BY BLOCK   C
C**********************************************************************C
      ICOUNT = 0
C
C     LOOP OVER ATOMIC CENTERS
      DO ICNT=1,NCNT
C
C       LOOP OVER KAPPA VALUES FOR THIS NUCLEAR CENTER
        DO KN=1,NKAP(ICNT)
C
C         IMPORT KAPPA, MAXIMUM MQN
          KAPPA = KVALS(KN,ICNT)
          MJMAX  = 2*IABS(KAPPA)-1
C
C         LOOP OVER MQN VALUES AND RECORD INDEX
          DO MJ=1,MJMAX,2
            ICOUNT = ICOUNT+1
            INDEX(ICNT,KAPPA,MJ) = ICOUNT
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     SET THE LIMITS FOR RUNNING OVER DENSITY COMBINATIONS             C
C**********************************************************************C
C
      IF(HMLTN.EQ.'NORL') THEN
        ITSTRT = 1
        ITSTOP = 1
        ITSKIP = 1
      ELSE
        ITSTRT = 1
        ITSTOP = 4
        ITSKIP = 3
      ENDIF
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTERS A AND B (USE INDEX 1000)      C
C**********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR A AND B: TT = LL (1) or SS (4)
      DO 1000 IT1=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITQN(1) = IT1
C
C       CALCULATE STARTING ADDRESS
        IF(IT1.EQ.1) THEN
          NADDAB = 0
        ELSE
          NADDAB = NSHIFT
        ENDIF
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
      DO 1000 ICNTB=1,ICNTA
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
        NFUNS(1) = NFUNCT(LQN(1)+1,ICNTA)
        DO IBAS=1,NFUNS(1)
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
        NFUNS(2) = NFUNCT(LQN(2)+1,ICNTB)
        DO JBAS=1, NFUNS(2)
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
      IEAB = 1
      IABLL = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB)
      IABSS = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB)
C
C**********************************************************************C
C      SECOND LAYER OF LOOPS, OVER CENTERS C AND D (INDEX 2000/2500)   C
C**********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR C AND D: T'T' = LL (1) or SS (4)
      DO 2000 IT2=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITQN(2) = IT2
C
C       CALCULATE STARTING ADDRESS
        IF(IT2.EQ.1) THEN
          NADDCD = 0
        ELSE
          NADDCD = NSHIFT
        ENDIF
C
C       TWO-ELECTRON COMPONENT OVERLAP INDEX
        IF(IT1.EQ.1.AND.IT2.EQ.1) THEN
C         (LL|LL)
          ITT = 1
        ELSEIF(IT1.EQ.1.AND.IT2.EQ.4) THEN
C         (LL|SS)
          ITT = 2
        ELSEIF(IT1.EQ.4.AND.IT2.EQ.1) THEN
C         (SS|LL)
          ITT = 3
        ELSEIF(IT1.EQ.4.AND.IT2.EQ.4) THEN
C         (SS|SS)
          ITT = 4
        ENDIF
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
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER C AND D
        IF(ICNTC.EQ.ICNTD) THEN
          INUCCD = 1
        ELSE
          INUCCD = 0
        ENDIF
C
C       NUMBER OF NUCLEAR CENTERS INVOLVED IN THIS OVERLAP
        MCNT = NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
C
C ***   STAGES: DECISION TREE FOR SKIPPING MULTI-CENTER CONTRIBUTIONS
        IF(MCNT.NE.1) THEN
C
C >>      STAGE 1: INCLUDE ONLY (LL|LL) REPULSION INTEGRALS
          IF(ILEV.EQ.1.AND.IT1+IT2.GT.2) THEN
            GOTO 2200
          ENDIF
C
C >>      STAGE 2: INCLUDE ONLY (LL|SS) AND (SS|LL) REPULSION INTEGRALS
          IF(ILEV.EQ.2.AND.IT1+IT2.GT.5) THEN
            GOTO 2200
          ENDIF
C
C ***   END OF ILEV DECISION TREE
        ENDIF
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
        NFUNS(3) = NFUNCT(LQN(3)+1,ICNTC)
        DO KBAS=1,NFUNS(3)
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
        NFUNS(4) = NFUNCT(LQN(4)+1,ICNTD)
        DO LBAS=1,NFUNS(4)
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
      IECD = 1
      ICDLL = IAD0LL(ICNTC,ICNTD,KC,KD,MC,MD)
      ICDSS = IAD0SS(ICNTC,ICNTD,KC,KD,MC,MD)
C
C**********************************************************************C
C     FOR THIS CHOICE OF A,B,C AND D, COMPUTE ADDRESSES AND PHASES     C
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
C     FURTHER DEFINE STARTING ADDRESSES FOR {ABCD} BASIS LABELS
      IA1 = LARGE(ICNTA,KA,MJA  ) + NADDAB
      IB1 = LARGE(ICNTB,KB,MJB  ) + NADDAB
      IC1 = LARGE(ICNTC,KC,MJC  ) + NADDCD
      ID1 = LARGE(ICNTD,KD,MJD  ) + NADDCD
C
      IA2 = LARGE(ICNTA,KA,MJA+1) + NADDAB
      IB2 = LARGE(ICNTB,KB,MJB+1) + NADDAB
      IC2 = LARGE(ICNTC,KC,MJC+1) + NADDCD
      ID2 = LARGE(ICNTD,KD,MJD+1) + NADDCD
C
      JA1 = LARGE(ICNTA,KA,MJA  ) + NADDAB
      JB1 = LARGE(ICNTB,KB,MJB  ) + NADDAB
      JC1 = LARGE(ICNTC,KC,MJC  ) + NADDCD
      JD1 = LARGE(ICNTD,KD,MJD  ) + NADDCD
C
      JA2 = LARGE(ICNTA,KA,MJA+1) + NADDAB
      JB2 = LARGE(ICNTB,KB,MJB+1) + NADDAB
      JC2 = LARGE(ICNTC,KC,MJC+1) + NADDCD
      JD2 = LARGE(ICNTD,KD,MJD+1) + NADDCD
C
C     CALCULATE KQN PHASE FACTORS FOR PERMUTING INTEGRALS
      IF((KQN(1)*KQN(2)).GT.0) THEN 
        PKAB = 1.0D0
      ELSE
        PKAB =-1.0D0
      ENDIF
C        
      IF((KQN(3)*KQN(4)).GT.0) THEN 
        PKCD = 1.0D0
      ELSE
        PKCD =-1.0D0
      ENDIF
C
C     CALCULATE MQN PHASE FACTORS FOR PERMUTING INTEGRALS
      PMAB1 = DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PMAB2 = DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      PMCD1 = DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PMCD2 = DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
C     COMBINATIONS OF PHASE FACTORS TO BE USED FOR COULOMB
      F1 = PKAB*PMAB1
      F2 = PKAB*PMAB2
      G1 = PKCD*PMCD1
      G2 = PKCD*PMCD2
C
C**********************************************************************C
C     SKIP BATCHES WHICH TRANSGRESS AGAINST QED SELECTION RULES        C
C**********************************************************************C
C
CC     LQN SELECTION RULES FOR ATOMIC ELEMENTS
C      IF(MCNT.EQ.1) THEN
C        IF(LQN(1).EQ.LQN(2).AND.LQN(3).EQ.LQN(4)) GOTO 2996
C        IF(LQN(1).EQ.LQN(3).AND.LQN(2).EQ.LQN(4)) GOTO 2996
C        IF(LQN(1).EQ.LQN(4).AND.LQN(2).EQ.LQN(3)) GOTO 2996
C        GOTO 2999
C      ENDIF
C2996  CONTINUE
CC
CC    MQN SELECTION RULES FOR ATOMIC ELEMENTS (CAN DO BETTER THAN THIS)
C      IF(MCNT.EQ.1) THEN
C        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) GOTO 2997
C        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) GOTO 2997
C        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) GOTO 2997
C        GOTO 2999
C      ENDIF    
C2997  CONTINUE
C
C     DIATOMIC MOLECULES CARRY STRICT SELECTION RULES ON MQNS
      IF(NCNT.LE.2) THEN
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) GOTO 2998
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) GOTO 2998
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) GOTO 2998
        GOTO 2999
      ENDIF
2998  CONTINUE
C
C**********************************************************************C
C     SKIP BATCHES WHICH CONFORM TO ERI PERMUTATION SYMMETRIES         C
C**********************************************************************C
C
C     DECISION TREE FOR SKIPPING CONTRIBUTIONS DUE TO INTEGRAL SYMMETRY
      IF(IQ1.LT.IQ2) GOTO 2999
      IF(IQ3.LT.IQ4) GOTO 2999
      IF(IQL.LT.IQR) GOTO 2999
C
C     INDICATE BLOCKS TO BE INCLUDED AHEAD GIVEN A,B,C,D BASIS QNMS...
C     A =/= B AND C =/= D WITH AB LIST VALUE =/= CD LIST VALUE
      IF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQL.GT.IQR) THEN
        ITSCF = 1
C     A=/=B AND C=/=D WITH IND(AB) = IND(CD)
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQL.EQ.IQR) THEN
        ITSCF = 2
C     A=/=B AND C = D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.EQ.IQ4.AND.IQL.GT.IQR) THEN
        ITSCF = 3
C     A = B AND C=/=D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.GT.IQ4.AND.IQL.GT.IQR) THEN
        ITSCF = 4
C     A = B AND C = D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQL.GT.IQR) THEN
        ITSCF = 5
C     A = B AND C = D WITH IND(AB) = IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQL.EQ.IQR) THEN
        ITSCF = 6
C     COMBINATION OF A,B,C,D NOT TO BE INCLUDED -- USE MATRIX CONJ LATER
      ELSE
        GO TO 2999
      ENDIF
C
C     READ IN FLAG VALUES FROM ISCF DATA BLOCK
      DO N=1,11
        IFLG(N) = ISCF(N,ITSCF)
      ENDDO
C
C     INCLUDE SPECIAL CASES...
C     A=/=B AND C=/=D WITH IND(AB)=/=IND(CD)
      IF(ITSCF.EQ.1) THEN
C       A = C
        IF(IQ1.EQ.IQ3) IFLG( 6) = 1
C       B = D
        IF(IQ2.EQ.IQ4) IFLG(11) = 1
C       B = C
        IF(IQ2.EQ.IQ3) IFLG( 8) = 1
C     A=/=B AND C = D WITH IND(AB)=/=IND(CD)
      ELSEIF(ITSCF.EQ.3) THEN
C       B = C
        IF(IQ2.EQ.IQ3) IFLG( 8) = 1
C     A = B AND C=/=D WITH IND(AB)=/=IND(CD)
      ELSEIF(ITSCF.EQ.4) THEN
C       B = C
        IF(IQ2.EQ.IQ3) IFLG( 8) = 1
      ENDIF
C
C**********************************************************************C
C     THIRD LAYER OF LOOPS, OVER BASIS FUNCTIONS A AND B (3000)        C
C -------------------------------------------------------------------- C
C     THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH         C
C     GENERATE CLOSED AND OPEN COULOMB MATRICES FROM SPINOR INTEGRALS. C
C     THESE INCLUDE PHASE FACTORS FOR THE PERMUTATION OF               C
C     KQN(1) <-> KQN(2) AND MQN(1) <-> MQN(2)                          C
C -------------------------------------------------------------------- C
C     (RSCF 86, 87)                                                    C
C**********************************************************************C
C
C     LOOP OVER ELEMENTS OF FOCK MATRIX BLOCK
      DO 3000 IBAS=1,NFUNS(1)
      DO 3000 JBAS=1,NFUNS(2)
C
C**********************************************************************C
C     SKIP BATCHES THAT PASS THE SCHWARZ SCREENING PROCEDURE           C
C**********************************************************************C
C
C       RECORD CPU TIME
        CALL CPU_TIME(TBCH1)
C
C       UPDATE COUNTER FOR NUMBER OF BLOCKS
        N2ET(MCNT,ITT) = N2ET(MCNT,ITT) + 1
C
        GOTO 1111
C       MATRIX ELEMENT UPPER BOUND TOLERANCE VALUE
        SENS = 1.0D-12
C
        CALL CPU_TIME(TX1)
C
C       TESTING: ONLY CHECK THE DIRECT TERMS
        COULAB = MAX(GDSC(IA1+IBAS,JB1+JBAS),
     &               GDSC(IA1+IBAS,JB2+JBAS),
     &               GDSC(IA2+IBAS,JB1+JBAS),
     &               GDSC(IA2+IBAS,JB2+JBAS))
C
C       DENSITY CONTRIBUTIONS
        DENTAB = MAX(ABS(DENT(IA1+IBAS,JB1+JBAS)),
     &               ABS(DENT(IA1+IBAS,JB2+JBAS)),
     &               ABS(DENT(IA2+IBAS,JB1+JBAS)),
     &               ABS(DENT(IA2+IBAS,JB2+JBAS)))
        DENTBA = MAX(ABS(DENT(JB1+JBAS,IA1+IBAS)),
     &               ABS(DENT(JB1+JBAS,IA2+IBAS)),
     &               ABS(DENT(JB2+JBAS,IA1+IBAS)),
     &               ABS(DENT(JB2+JBAS,IA2+IBAS)))
C
C       LOOP OVER SECOND PAIR OF BASIS FUNCTIONS
        M = 0
        DO KBAS=1,NFUNS(3)
          DO LBAS=1,NFUNS(4)
            M = M+1
C
C           DIRECT CONTRIBUTIONS
            COULCD = MAX(GDSC(IC1+KBAS,JD1+LBAS),
     &                   GDSC(IC1+KBAS,JD2+LBAS),
     &                   GDSC(IC2+KBAS,JD1+LBAS),
     &                   GDSC(IC2+KBAS,JD2+LBAS))
C
C           DENSITY CONTRIBUTIONS
            DENTCD = MAX(ABS(DENT(IC1+KBAS,JD1+LBAS)),
     &                   ABS(DENT(IC1+KBAS,JD2+LBAS)),
     &                   ABS(DENT(IC2+KBAS,JD1+LBAS)),
     &                   ABS(DENT(IC2+KBAS,JD2+LBAS)))
            DENTDC = MAX(ABS(DENT(JD1+LBAS,IC1+KBAS)),
     &                   ABS(DENT(JD1+LBAS,IC2+KBAS)),
     &                   ABS(DENT(JD2+LBAS,IC1+KBAS)),
     &                   ABS(DENT(JD2+LBAS,IC2+KBAS)))
C
C           MATRIX ELEMENT UPPER BOUND
            EPS1 = COULAB*COULCD*DENTAB
            EPS2 = COULAB*COULCD*DENTBA
            EPS3 = COULAB*COULCD*DENTCD
            EPS4 = COULAB*COULCD*DENTDC

            EPSTOT = MAX(EPS1,EPS2,EPS3,EPS4)
C
C           TEST WHETHER THE UPPER BOUNDS ESTIMATES ARE TOO BIG
            IF(EPSTOT.GT.SENS) THEN
C             MUST CALCULATE THIS BLOCK
              GOTO 3997
            ENDIF
C
          ENDDO
        ENDDO
C
C       ALL MATRIX ELEMENT ESTIMATES BELOW THE TOLERANCE VALUE
C
C       UPDATE COUNTER FOR NUMBER OF BLOCKS
        N2ES(MCNT,ITT) = N2ES(MCNT,ITT) + 1
C
C       SKIP THIS BLOCK
        GOTO 3999
3997    CONTINUE
C
C       UPDATE TIME SPENT SCREENING
        CALL CPU_TIME(TX2)
        TSCR = TSCR + TX2 - TX1
1111    CONTINUE
C
C**********************************************************************C
C     ALL SCHWARZ SCREENING AND INTEGRAL SYMMETRIES FINISHED           C
C**********************************************************************C
C
C       GENERATE BATCH OF ELECTRON REPULSION INTEGRALS
        CALL ERI(RR,XYZ,KQN,MQN,EXPT,NFUNS,ITQN,IBAS,JBAS,IEAB,IECD)
C
C**********************************************************************C
C     CONSTRUCT CLOSED-SHELL COULOMB INTERACTION MATRICES, GDIR/GXCH.  C
C**********************************************************************C
C
C       THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH
C       GENERATE THE GDIR AND GXCH MATRIX FROM THE SPINOR INTEGRALS.
C
C       FIRST IFLG BATCH (DIRECT)
        IF(IFLG(1).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GDIR(IA1+IBAS,JB1+JBAS) = GDIR(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 1)*DENT(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M, 2)*DENT(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M, 3)*DENT(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M, 4)*DENT(IC2+KBAS,JD2+LBAS)
     &                         +    G1*RR(M, 4)*DENT(JD1+LBAS,IC1+KBAS)
     &                         +    G2*RR(M, 2)*DENT(JD1+LBAS,IC2+KBAS)
     &                         +    G2*RR(M, 3)*DENT(JD2+LBAS,IC1+KBAS)
     &                         +    G1*RR(M, 1)*DENT(JD2+LBAS,IC2+KBAS)
C
              GDIR(IA1+IBAS,JB2+JBAS) = GDIR(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M, 5)*DENT(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M, 6)*DENT(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M, 7)*DENT(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M, 8)*DENT(IC2+KBAS,JD2+LBAS)
     &                         +    G1*RR(M, 8)*DENT(JD1+LBAS,IC1+KBAS)
     &                         +    G2*RR(M, 6)*DENT(JD1+LBAS,IC2+KBAS)
     &                         +    G2*RR(M, 7)*DENT(JD2+LBAS,IC1+KBAS)
     &                         +    G1*RR(M, 5)*DENT(JD2+LBAS,IC2+KBAS)

C
              GDIR(IA2+IBAS,JB1+JBAS) = GDIR(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M, 9)*DENT(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M,10)*DENT(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M,11)*DENT(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M,12)*DENT(IC2+KBAS,JD2+LBAS)
     &                         +    G1*RR(M,12)*DENT(JD1+LBAS,IC1+KBAS)
     &                         +    G2*RR(M,10)*DENT(JD1+LBAS,IC2+KBAS)
     &                         +    G2*RR(M,11)*DENT(JD2+LBAS,IC1+KBAS)
     &                         +    G1*RR(M, 9)*DENT(JD2+LBAS,IC2+KBAS)
C
              GDIR(IA2+IBAS,JB2+JBAS) = GDIR(IA2+IBAS,JB2+JBAS)
     &                         +       RR(M,13)*DENT(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M,14)*DENT(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M,15)*DENT(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M,16)*DENT(IC2+KBAS,JD2+LBAS)
     &                         +    G1*RR(M,16)*DENT(JD1+LBAS,IC1+KBAS)
     &                         +    G2*RR(M,14)*DENT(JD1+LBAS,IC2+KBAS)
     &                         +    G2*RR(M,15)*DENT(JD2+LBAS,IC1+KBAS)
     &                         +    G1*RR(M,13)*DENT(JD2+LBAS,IC2+KBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SECOND IFLG BATCH (DIRECT)
        IF(IFLG(2).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GDIR(IA1+IBAS,JB1+JBAS) = GDIR(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 1)*DENT(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M, 2)*DENT(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M, 3)*DENT(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M, 4)*DENT(IC2+KBAS,JD2+LBAS)
C
              GDIR(IA1+IBAS,JB2+JBAS) = GDIR(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M, 5)*DENT(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M, 6)*DENT(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M, 7)*DENT(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M, 8)*DENT(IC2+KBAS,JD2+LBAS)
C
              GDIR(IA2+IBAS,JB1+JBAS) = GDIR(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M, 9)*DENT(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M,10)*DENT(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M,11)*DENT(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M,12)*DENT(IC2+KBAS,JD2+LBAS)
C
              GDIR(IA2+IBAS,JB2+JBAS) = GDIR(IA2+IBAS,JB2+JBAS)
     &                         +       RR(M,13)*DENT(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M,14)*DENT(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M,15)*DENT(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M,16)*DENT(IC2+KBAS,JD2+LBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       THIRD IFLG BATCH (DIRECT)
        IF(IFLG(3).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GDIR(IC1+KBAS,JD1+LBAS) = GDIR(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M, 1)*DENT(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 5)*DENT(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M, 9)*DENT(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,13)*DENT(IA2+IBAS,JB2+JBAS)
     &                         +    F1*RR(M,13)*DENT(JB1+JBAS,IA1+IBAS)
     &                         +    F2*RR(M, 5)*DENT(JB1+JBAS,IA2+IBAS)
     &                         +    F2*RR(M, 9)*DENT(JB2+JBAS,IA1+IBAS)
     &                         +    F1*RR(M, 1)*DENT(JB2+JBAS,IA2+IBAS)
C
              GDIR(IC1+KBAS,JD2+LBAS) = GDIR(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M, 2)*DENT(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 6)*DENT(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M,10)*DENT(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,14)*DENT(IA2+IBAS,JB2+JBAS)
     &                         +    F1*RR(M,14)*DENT(JB1+JBAS,IA1+IBAS)
     &                         +    F2*RR(M, 6)*DENT(JB1+JBAS,IA2+IBAS)
     &                         +    F2*RR(M,10)*DENT(JB2+JBAS,IA1+IBAS)
     &                         +    F1*RR(M, 2)*DENT(JB2+JBAS,IA2+IBAS)
C
              GDIR(IC2+KBAS,JD1+LBAS) = GDIR(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M, 3)*DENT(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 7)*DENT(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M,11)*DENT(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,15)*DENT(IA2+IBAS,JB2+JBAS)
     &                         +    F1*RR(M,15)*DENT(JB1+JBAS,IA1+IBAS)
     &                         +    F2*RR(M, 7)*DENT(JB1+JBAS,IA2+IBAS)
     &                         +    F2*RR(M,11)*DENT(JB2+JBAS,IA1+IBAS)
     &                         +    F1*RR(M, 3)*DENT(JB2+JBAS,IA2+IBAS)
C
              GDIR(IC2+KBAS,JD2+LBAS) = GDIR(IC2+KBAS,JD2+LBAS)
     &                         +       RR(M, 4)*DENT(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 8)*DENT(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M,12)*DENT(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,16)*DENT(IA2+IBAS,JB2+JBAS)
     &                         +    F1*RR(M,16)*DENT(JB1+JBAS,IA1+IBAS)
     &                         +    F2*RR(M, 8)*DENT(JB1+JBAS,IA2+IBAS)
     &                         +    F2*RR(M,12)*DENT(JB2+JBAS,IA1+IBAS)
     &                         +    F1*RR(M, 4)*DENT(JB2+JBAS,IA2+IBAS)
C 
            ENDDO
          ENDDO
        ENDIF
C
C       FOURTH IFLG BATCH (DIRECT)
        IF(IFLG(4).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GDIR(IC1+KBAS,JD1+LBAS) = GDIR(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M, 1)*DENT(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 5)*DENT(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M, 9)*DENT(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,13)*DENT(IA2+IBAS,JB2+JBAS)
C
              GDIR(IC1+KBAS,JD2+LBAS) = GDIR(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M, 2)*DENT(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 6)*DENT(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M,10)*DENT(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,14)*DENT(IA2+IBAS,JB2+JBAS)
C
              GDIR(IC2+KBAS,JD1+LBAS) = GDIR(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M, 3)*DENT(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 7)*DENT(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M,11)*DENT(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,15)*DENT(IA2+IBAS,JB2+JBAS)
C
              GDIR(IC2+KBAS,JD2+LBAS) = GDIR(IC2+KBAS,JD2+LBAS)
     &                         +       RR(M, 4)*DENT(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 8)*DENT(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M,12)*DENT(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,16)*DENT(IA2+IBAS,JB2+JBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       FIFTH IFLG BATCH (EXCHANGE)
        IF(IFLG(5).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GXCH(IA1+IBAS,JC1+KBAS) = GXCH(IA1+IBAS,JC1+KBAS)
     &                         +    G1*RR(M, 4)*DENT(ID1+LBAS,JB1+JBAS)
     &                         +    G1*RR(M, 8)*DENT(ID1+LBAS,JB2+JBAS)
     &                         +    G2*RR(M, 3)*DENT(ID2+LBAS,JB1+JBAS)
     &                         +    G2*RR(M, 7)*DENT(ID2+LBAS,JB2+JBAS)
C
              GXCH(IA1+IBAS,JC2+KBAS) = GXCH(IA1+IBAS,JC2+KBAS)
     &                         +    G2*RR(M, 2)*DENT(ID1+LBAS,JB1+JBAS)
     &                         +    G2*RR(M, 6)*DENT(ID1+LBAS,JB2+JBAS)
     &                         +    G1*RR(M, 1)*DENT(ID2+LBAS,JB1+JBAS)
     &                         +    G1*RR(M, 5)*DENT(ID2+LBAS,JB2+JBAS)
C
              GXCH(IA2+IBAS,JC1+KBAS) = GXCH(IA2+IBAS,JC1+KBAS)
     &                         +    G1*RR(M,12)*DENT(ID1+LBAS,JB1+JBAS)
     &                         +    G1*RR(M,16)*DENT(ID1+LBAS,JB2+JBAS)
     &                         +    G2*RR(M,11)*DENT(ID2+LBAS,JB1+JBAS)
     &                         +    G2*RR(M,15)*DENT(ID2+LBAS,JB2+JBAS)
C
              GXCH(IA2+IBAS,JC2+KBAS) = GXCH(IA2+IBAS,JC2+KBAS)
     &                         +    G2*RR(M,10)*DENT(ID1+LBAS,JB1+JBAS)
     &                         +    G2*RR(M,14)*DENT(ID1+LBAS,JB2+JBAS)
     &                         +    G1*RR(M, 9)*DENT(ID2+LBAS,JB1+JBAS)
     &                         +    G1*RR(M,13)*DENT(ID2+LBAS,JB2+JBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SIXTH IFLG BATCH (EXCHANGE)
        IF(IFLG(6).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GXCH(IC1+KBAS,JA1+IBAS) = GXCH(IC1+KBAS,JA1+IBAS)
     &                         +    F1*RR(M,13)*DENT(JB1+JBAS,ID1+LBAS)
     &                         +    F1*RR(M,14)*DENT(JB1+JBAS,ID2+LBAS)
     &                         +    F2*RR(M, 9)*DENT(JB2+JBAS,ID1+LBAS)
     &                         +    F2*RR(M,10)*DENT(JB2+JBAS,ID2+LBAS)
C
              GXCH(IC1+KBAS,JA2+IBAS) = GXCH(IC1+KBAS,JA2+IBAS)
     &                         +    F2*RR(M, 5)*DENT(JB1+JBAS,ID1+LBAS)
     &                         +    F2*RR(M, 6)*DENT(JB1+JBAS,ID2+LBAS)
     &                         +    F1*RR(M, 1)*DENT(JB2+JBAS,ID1+LBAS)
     &                         +    F1*RR(M, 2)*DENT(JB2+JBAS,ID2+LBAS)
C
              GXCH(IC2+KBAS,JA1+IBAS) = GXCH(IC2+KBAS,JA1+IBAS)
     &                         +    F1*RR(M,15)*DENT(JB1+JBAS,ID1+LBAS)
     &                         +    F1*RR(M,16)*DENT(JB1+JBAS,ID2+LBAS)
     &                         +    F2*RR(M,11)*DENT(JB2+JBAS,ID1+LBAS)
     &                         +    F2*RR(M,12)*DENT(JB2+JBAS,ID2+LBAS)
C
              GXCH(IC2+KBAS,JA2+IBAS) = GXCH(IC2+KBAS,JA2+IBAS)
     &                         +    F2*RR(M, 7)*DENT(JB1+JBAS,ID1+LBAS)
     &                         +    F2*RR(M, 8)*DENT(JB1+JBAS,ID2+LBAS)
     &                         +    F1*RR(M, 3)*DENT(JB2+JBAS,ID1+LBAS)
     &                         +    F1*RR(M, 4)*DENT(JB2+JBAS,ID2+LBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SEVENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(7).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GXCH(IB1+JBAS,JC1+KBAS) = GXCH(IB1+JBAS,JC1+KBAS)
     &                         + F1*G1*RR(M,16)*DENT(ID1+LBAS,JA1+IBAS)
     &                         + F2*G1*RR(M, 8)*DENT(ID1+LBAS,JA2+IBAS)
     &                         + F1*G2*RR(M,15)*DENT(ID2+LBAS,JA1+IBAS)
     &                         + F2*G2*RR(M, 7)*DENT(ID2+LBAS,JA2+IBAS)
C
              GXCH(IB1+JBAS,JC2+KBAS) = GXCH(IB1+JBAS,JC2+KBAS)
     &                         + F1*G2*RR(M,14)*DENT(ID1+LBAS,JA1+IBAS)
     &                         + F2*G2*RR(M, 6)*DENT(ID1+LBAS,JA2+IBAS)
     &                         + F1*G1*RR(M,13)*DENT(ID2+LBAS,JA1+IBAS)
     &                         + F2*G1*RR(M, 5)*DENT(ID2+LBAS,JA2+IBAS)
C
              GXCH(IB2+JBAS,JC1+KBAS) = GXCH(IB2+JBAS,JC1+KBAS)
     &                         + F2*G1*RR(M,12)*DENT(ID1+LBAS,JA1+IBAS)
     &                         + F1*G1*RR(M, 4)*DENT(ID1+LBAS,JA2+IBAS)
     &                         + F2*G2*RR(M,11)*DENT(ID2+LBAS,JA1+IBAS)
     &                         + F1*G2*RR(M, 3)*DENT(ID2+LBAS,JA2+IBAS)
C
              GXCH(IB2+JBAS,JC2+KBAS) = GXCH(IB2+JBAS,JC2+KBAS)
     &                         + F2*G2*RR(M,10)*DENT(ID1+LBAS,JA1+IBAS)
     &                         + F1*G2*RR(M, 2)*DENT(ID1+LBAS,JA2+IBAS)
     &                         + F2*G1*RR(M, 9)*DENT(ID2+LBAS,JA1+IBAS)
     &                         + F1*G1*RR(M, 1)*DENT(ID2+LBAS,JA2+IBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       EIGHTH IFLG BATCH (EXCHANGE)
        IF(IFLG(8).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GXCH(IC1+KBAS,JB1+JBAS) = GXCH(IC1+KBAS,JB1+JBAS)
     &                         +       RR(M, 1)*DENT(JA1+IBAS,ID1+LBAS)
     &                         +       RR(M, 2)*DENT(JA1+IBAS,ID2+LBAS)
     &                         +       RR(M, 9)*DENT(JA2+IBAS,ID1+LBAS)
     &                         +       RR(M,10)*DENT(JA2+IBAS,ID2+LBAS)
C
              GXCH(IC1+KBAS,JB2+JBAS) = GXCH(IC1+KBAS,JB2+JBAS)
     &                         +       RR(M, 5)*DENT(JA1+IBAS,ID1+LBAS)
     &                         +       RR(M, 6)*DENT(JA1+IBAS,ID2+LBAS)
     &                         +       RR(M,13)*DENT(JA2+IBAS,ID1+LBAS)
     &                         +       RR(M,14)*DENT(JA2+IBAS,ID2+LBAS)
C
              GXCH(IC2+KBAS,JB1+JBAS) = GXCH(IC2+KBAS,JB1+JBAS)
     &                         +       RR(M, 3)*DENT(JA1+IBAS,ID1+LBAS)
     &                         +       RR(M, 4)*DENT(JA1+IBAS,ID2+LBAS)
     &                         +       RR(M,11)*DENT(JA2+IBAS,ID1+LBAS)
     &                         +       RR(M,12)*DENT(JA2+IBAS,ID2+LBAS)
C
              GXCH(IC2+KBAS,JB2+JBAS) = GXCH(IC2+KBAS,JB2+JBAS)
     &                         +       RR(M, 7)*DENT(JA1+IBAS,ID1+LBAS)
     &                         +       RR(M, 8)*DENT(JA1+IBAS,ID2+LBAS)
     &                         +       RR(M,15)*DENT(JA2+IBAS,ID1+LBAS)
     &                         +       RR(M,16)*DENT(JA2+IBAS,ID2+LBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       NINTH IFLG BATCH (EXCHANGE)
        IF(IFLG(9).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C  
              GXCH(IA1+IBAS,JD1+LBAS) = GXCH(IA1+IBAS,JD1+LBAS)
     &                         +       RR(M, 1)*DENT(IC1+KBAS,JB1+JBAS)
     &                         +       RR(M, 5)*DENT(IC1+KBAS,JB2+JBAS)
     &                         +       RR(M, 3)*DENT(IC2+KBAS,JB1+JBAS)
     &                         +       RR(M, 7)*DENT(IC2+KBAS,JB2+JBAS)
C
              GXCH(IA1+IBAS,JD2+LBAS) = GXCH(IA1+IBAS,JD2+LBAS)
     &                         +       RR(M, 2)*DENT(IC1+KBAS,JB1+JBAS)
     &                         +       RR(M, 6)*DENT(IC1+KBAS,JB2+JBAS)
     &                         +       RR(M, 4)*DENT(IC2+KBAS,JB1+JBAS)
     &                         +       RR(M, 8)*DENT(IC2+KBAS,JB2+JBAS)
C
              GXCH(IA2+IBAS,JD1+LBAS) = GXCH(IA2+IBAS,JD1+LBAS)
     &                         +       RR(M, 9)*DENT(IC1+KBAS,JB1+JBAS)
     &                         +       RR(M,13)*DENT(IC1+KBAS,JB2+JBAS)
     &                         +       RR(M,11)*DENT(IC2+KBAS,JB1+JBAS)
     &                         +       RR(M,15)*DENT(IC2+KBAS,JB2+JBAS)
C
              GXCH(IA2+IBAS,JD2+LBAS) = GXCH(IA2+IBAS,JD2+LBAS)
     &                         +       RR(M,10)*DENT(IC1+KBAS,JB1+JBAS)
     &                         +       RR(M,14)*DENT(IC1+KBAS,JB2+JBAS)
     &                         +       RR(M,12)*DENT(IC2+KBAS,JB1+JBAS)
     &                         +       RR(M,16)*DENT(IC2+KBAS,JB2+JBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       TENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(10).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C
              GXCH(IB1+JBAS,JD1+LBAS) = GXCH(IB1+JBAS,JD1+LBAS)
     &                         +    F1*RR(M,13)*DENT(IC1+KBAS,JA1+IBAS)
     &                         +    F2*RR(M, 5)*DENT(IC1+KBAS,JA2+IBAS)
     &                         +    F1*RR(M,15)*DENT(IC2+KBAS,JA1+IBAS)
     &                         +    F2*RR(M, 7)*DENT(IC2+KBAS,JA2+IBAS)
C
              GXCH(IB1+JBAS,JD2+LBAS) = GXCH(IB1+JBAS,JD2+LBAS)
     &                         +    F1*RR(M,14)*DENT(IC1+KBAS,JA1+IBAS)
     &                         +    F2*RR(M, 6)*DENT(IC1+KBAS,JA2+IBAS)
     &                         +    F1*RR(M,16)*DENT(IC2+KBAS,JA1+IBAS)
     &                         +    F2*RR(M, 8)*DENT(IC2+KBAS,JA2+IBAS)
C
              GXCH(IB2+JBAS,JD1+LBAS) = GXCH(IB2+JBAS,JD1+LBAS)
     &                         +    F2*RR(M, 9)*DENT(IC1+KBAS,JA1+IBAS)
     &                         +    F1*RR(M, 1)*DENT(IC1+KBAS,JA2+IBAS)
     &                         +    F2*RR(M,11)*DENT(IC2+KBAS,JA1+IBAS)
     &                         +    F1*RR(M, 3)*DENT(IC2+KBAS,JA2+IBAS)
C
              GXCH(IB2+JBAS,JD2+LBAS) = GXCH(IB2+JBAS,JD2+LBAS)
     &                         +    F2*RR(M,10)*DENT(IC1+KBAS,JA1+IBAS)
     &                         +    F1*RR(M, 2)*DENT(IC1+KBAS,JA2+IBAS)
     &                         +    F2*RR(M,12)*DENT(IC2+KBAS,JA1+IBAS)
     &                         +    F1*RR(M, 4)*DENT(IC2+KBAS,JA2+IBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       ELEVENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(11).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C
              GXCH(ID1+LBAS,JB1+JBAS) = GXCH(ID1+LBAS,JB1+JBAS)
     &                         +    G1*RR(M, 4)*DENT(JA1+IBAS,IC1+KBAS)
     &                         +    G2*RR(M, 2)*DENT(JA1+IBAS,IC2+KBAS)
     &                         +    G1*RR(M,12)*DENT(JA2+IBAS,IC1+KBAS)
     &                         +    G2*RR(M,10)*DENT(JA2+IBAS,IC2+KBAS)
C
              GXCH(ID1+LBAS,JB2+JBAS) = GXCH(ID1+LBAS,JB2+JBAS)
     &                         +    G1*RR(M, 8)*DENT(JA1+IBAS,IC1+KBAS)
     &                         +    G2*RR(M, 6)*DENT(JA1+IBAS,IC2+KBAS)
     &                         +    G1*RR(M,16)*DENT(JA2+IBAS,IC1+KBAS)
     &                         +    G2*RR(M,14)*DENT(JA2+IBAS,IC2+KBAS)
C
              GXCH(ID2+LBAS,JB1+JBAS) = GXCH(ID2+LBAS,JB1+JBAS)
     &                         +    G2*RR(M, 3)*DENT(JA1+IBAS,IC1+KBAS)
     &                         +    G1*RR(M, 1)*DENT(JA1+IBAS,IC2+KBAS)
     &                         +    G2*RR(M,11)*DENT(JA2+IBAS,IC1+KBAS)
     &                         +    G1*RR(M, 9)*DENT(JA2+IBAS,IC2+KBAS)
C
              GXCH(ID2+LBAS,JB2+JBAS) = GXCH(ID2+LBAS,JB2+JBAS)
     &                         +    G2*RR(M, 7)*DENT(JA1+IBAS,IC1+KBAS)
     &                         +    G1*RR(M, 5)*DENT(JA1+IBAS,IC2+KBAS)
     &                         +    G2*RR(M,15)*DENT(JA2+IBAS,IC1+KBAS)
     &                         +    G1*RR(M,13)*DENT(JA2+IBAS,IC2+KBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C**********************************************************************C
C     CONSTRUCT OPEN-SHELL COULOMB INTERACTION MATRICES, QDIR/QXCH.    C
C**********************************************************************C
C
        IF(NOPN.EQ.0) GOTO 399
C
C       THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH
C       GENERATE THE GDIR AND GXCH MATRIX FROM THE SPINOR INTEGRALS.
C
C       FIRST IFLG BATCH (DIRECT)
        IF(IFLG(1).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QDIR(IA1+IBAS,JB1+JBAS) = QDIR(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 1)*DENO(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M, 2)*DENO(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M, 3)*DENO(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M, 4)*DENO(IC2+KBAS,JD2+LBAS)
     &                         +    G1*RR(M, 4)*DENO(JD1+LBAS,IC1+KBAS)
     &                         +    G2*RR(M, 2)*DENO(JD1+LBAS,IC2+KBAS)
     &                         +    G2*RR(M, 3)*DENO(JD2+LBAS,IC1+KBAS)
     &                         +    G1*RR(M, 1)*DENO(JD2+LBAS,IC2+KBAS)
C
              QDIR(IA1+IBAS,JB2+JBAS) = QDIR(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M, 5)*DENO(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M, 6)*DENO(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M, 7)*DENO(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M, 8)*DENO(IC2+KBAS,JD2+LBAS)
     &                         +    G1*RR(M, 8)*DENO(JD1+LBAS,IC1+KBAS)
     &                         +    G2*RR(M, 6)*DENO(JD1+LBAS,IC2+KBAS)
     &                         +    G2*RR(M, 7)*DENO(JD2+LBAS,IC1+KBAS)
     &                         +    G1*RR(M, 5)*DENO(JD2+LBAS,IC2+KBAS)
C
              QDIR(IA2+IBAS,JB1+JBAS) = QDIR(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M, 9)*DENO(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M,10)*DENO(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M,11)*DENO(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M,12)*DENO(IC2+KBAS,JD2+LBAS)
     &                         +    G1*RR(M,12)*DENO(JD1+LBAS,IC1+KBAS)
     &                         +    G2*RR(M,10)*DENO(JD1+LBAS,IC2+KBAS)
     &                         +    G2*RR(M,11)*DENO(JD2+LBAS,IC1+KBAS)
     &                         +    G1*RR(M, 9)*DENO(JD2+LBAS,IC2+KBAS)
C
              QDIR(IA2+IBAS,JB2+JBAS) = QDIR(IA2+IBAS,JB2+JBAS)
     &                         +       RR(M,13)*DENO(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M,14)*DENO(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M,15)*DENO(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M,16)*DENO(IC2+KBAS,JD2+LBAS)
     &                         +    G1*RR(M,16)*DENO(JD1+LBAS,IC1+KBAS)
     &                         +    G2*RR(M,14)*DENO(JD1+LBAS,IC2+KBAS)
     &                         +    G2*RR(M,15)*DENO(JD2+LBAS,IC1+KBAS)
     &                         +    G1*RR(M,13)*DENO(JD2+LBAS,IC2+KBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SECOND IFLG BATCH (DIRECT)
        IF(IFLG(2).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QDIR(IA1+IBAS,JB1+JBAS) = QDIR(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 1)*DENO(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M, 2)*DENO(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M, 3)*DENO(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M, 4)*DENO(IC2+KBAS,JD2+LBAS)
C
              QDIR(IA1+IBAS,JB2+JBAS) = QDIR(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M, 5)*DENO(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M, 6)*DENO(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M, 7)*DENO(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M, 8)*DENO(IC2+KBAS,JD2+LBAS)
C
              QDIR(IA2+IBAS,JB1+JBAS) = QDIR(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M, 9)*DENO(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M,10)*DENO(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M,11)*DENO(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M,12)*DENO(IC2+KBAS,JD2+LBAS)
C
              QDIR(IA2+IBAS,JB2+JBAS) = QDIR(IA2+IBAS,JB2+JBAS)
     &                         +       RR(M,13)*DENO(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M,14)*DENO(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M,15)*DENO(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M,16)*DENO(IC2+KBAS,JD2+LBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       THIRD IFLG BATCH (DIRECT)
        IF(IFLG(3).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QDIR(IC1+KBAS,JD1+LBAS) = QDIR(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M, 1)*DENO(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 5)*DENO(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M, 9)*DENO(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,13)*DENO(IA2+IBAS,JB2+JBAS)
     &                         +    F1*RR(M,13)*DENO(JB1+JBAS,IA1+IBAS)
     &                         +    F2*RR(M, 5)*DENO(JB1+JBAS,IA2+IBAS)
     &                         +    F2*RR(M, 9)*DENO(JB2+JBAS,IA1+IBAS)
     &                         +    F1*RR(M, 1)*DENO(JB2+JBAS,IA2+IBAS)
C
              QDIR(IC1+KBAS,JD2+LBAS) = QDIR(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M, 2)*DENO(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 6)*DENO(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M,10)*DENO(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,14)*DENO(IA2+IBAS,JB2+JBAS)
     &                         +    F1*RR(M,14)*DENO(JB1+JBAS,IA1+IBAS)
     &                         +    F2*RR(M, 6)*DENO(JB1+JBAS,IA2+IBAS)
     &                         +    F2*RR(M,10)*DENO(JB2+JBAS,IA1+IBAS)
     &                         +    F1*RR(M, 2)*DENO(JB2+JBAS,IA2+IBAS)
C
              QDIR(IC2+KBAS,JD1+LBAS) = QDIR(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M, 3)*DENO(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 7)*DENO(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M,11)*DENO(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,15)*DENO(IA2+IBAS,JB2+JBAS)
     &                         +    F1*RR(M,15)*DENO(JB1+JBAS,IA1+IBAS)
     &                         +    F2*RR(M, 7)*DENO(JB1+JBAS,IA2+IBAS)
     &                         +    F2*RR(M,11)*DENO(JB2+JBAS,IA1+IBAS)
     &                         +    F1*RR(M, 3)*DENO(JB2+JBAS,IA2+IBAS)
C
              QDIR(IC2+KBAS,JD2+LBAS) = QDIR(IC2+KBAS,JD2+LBAS)
     &                         +       RR(M, 4)*DENO(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 8)*DENO(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M,12)*DENO(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,16)*DENO(IA2+IBAS,JB2+JBAS)
     &                         +    F1*RR(M,16)*DENO(JB1+JBAS,IA1+IBAS)
     &                         +    F2*RR(M, 8)*DENO(JB1+JBAS,IA2+IBAS)
     &                         +    F2*RR(M,12)*DENO(JB2+JBAS,IA1+IBAS)
     &                         +    F1*RR(M, 4)*DENO(JB2+JBAS,IA2+IBAS)
C 
            ENDDO
          ENDDO
        ENDIF
C
C       FOURTH IFLG BATCH (DIRECT)
        IF(IFLG(4).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QDIR(IC1+KBAS,JD1+LBAS) = QDIR(IC1+KBAS,JD1+LBAS)
     &                         +       RR(M, 1)*DENO(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 5)*DENO(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M, 9)*DENO(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,13)*DENO(IA2+IBAS,JB2+JBAS)
C
              QDIR(IC1+KBAS,JD2+LBAS) = QDIR(IC1+KBAS,JD2+LBAS)
     &                         +       RR(M, 2)*DENO(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 6)*DENO(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M,10)*DENO(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,14)*DENO(IA2+IBAS,JB2+JBAS)
C
              QDIR(IC2+KBAS,JD1+LBAS) = QDIR(IC2+KBAS,JD1+LBAS)
     &                         +       RR(M, 3)*DENO(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 7)*DENO(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M,11)*DENO(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,15)*DENO(IA2+IBAS,JB2+JBAS)
C
              QDIR(IC2+KBAS,JD2+LBAS) = QDIR(IC2+KBAS,JD2+LBAS)
     &                         +       RR(M, 4)*DENO(IA1+IBAS,JB1+JBAS)
     &                         +       RR(M, 8)*DENO(IA1+IBAS,JB2+JBAS)
     &                         +       RR(M,12)*DENO(IA2+IBAS,JB1+JBAS)
     &                         +       RR(M,16)*DENO(IA2+IBAS,JB2+JBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       FIFTH IFLG BATCH (EXCHANGE)
        IF(IFLG(5).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QXCH(IA1+IBAS,JC1+KBAS) = QXCH(IA1+IBAS,JC1+KBAS)
     &                         +    G1*RR(M, 4)*DENO(ID1+LBAS,JB1+JBAS)
     &                         +    G1*RR(M, 8)*DENO(ID1+LBAS,JB2+JBAS)
     &                         +    G2*RR(M, 3)*DENO(ID2+LBAS,JB1+JBAS)
     &                         +    G2*RR(M, 7)*DENO(ID2+LBAS,JB2+JBAS)
C
              QXCH(IA1+IBAS,JC2+KBAS) = QXCH(IA1+IBAS,JC2+KBAS)
     &                         +    G2*RR(M, 2)*DENO(ID1+LBAS,JB1+JBAS)
     &                         +    G2*RR(M, 6)*DENO(ID1+LBAS,JB2+JBAS)
     &                         +    G1*RR(M, 1)*DENO(ID2+LBAS,JB1+JBAS)
     &                         +    G1*RR(M, 5)*DENO(ID2+LBAS,JB2+JBAS)
C
              QXCH(IA2+IBAS,JC1+KBAS) = QXCH(IA2+IBAS,JC1+KBAS)
     &                         +    G1*RR(M,12)*DENO(ID1+LBAS,JB1+JBAS)
     &                         +    G1*RR(M,16)*DENO(ID1+LBAS,JB2+JBAS)
     &                         +    G2*RR(M,11)*DENO(ID2+LBAS,JB1+JBAS)
     &                         +    G2*RR(M,15)*DENO(ID2+LBAS,JB2+JBAS)
C
              QXCH(IA2+IBAS,JC2+KBAS) = QXCH(IA2+IBAS,JC2+KBAS)
     &                         +    G2*RR(M,10)*DENO(ID1+LBAS,JB1+JBAS)
     &                         +    G2*RR(M,14)*DENO(ID1+LBAS,JB2+JBAS)
     &                         +    G1*RR(M, 9)*DENO(ID2+LBAS,JB1+JBAS)
     &                         +    G1*RR(M,13)*DENO(ID2+LBAS,JB2+JBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SIXTH IFLG BATCH (EXCHANGE)
        IF(IFLG(6).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QXCH(IC1+KBAS,JA1+IBAS) = QXCH(IC1+KBAS,JA1+IBAS)
     &                         +    F1*RR(M,13)*DENO(JB1+JBAS,ID1+LBAS)
     &                         +    F1*RR(M,14)*DENO(JB1+JBAS,ID2+LBAS)
     &                         +    F2*RR(M, 9)*DENO(JB2+JBAS,ID1+LBAS)
     &                         +    F2*RR(M,10)*DENO(JB2+JBAS,ID2+LBAS)
C
              QXCH(IC1+KBAS,JA2+IBAS) = QXCH(IC1+KBAS,JA2+IBAS)
     &                         +    F2*RR(M, 5)*DENO(JB1+JBAS,ID1+LBAS)
     &                         +    F2*RR(M, 6)*DENO(JB1+JBAS,ID2+LBAS)
     &                         +    F1*RR(M, 1)*DENO(JB2+JBAS,ID1+LBAS)
     &                         +    F1*RR(M, 2)*DENO(JB2+JBAS,ID2+LBAS)
C
              QXCH(IC2+KBAS,JA1+IBAS) = QXCH(IC2+KBAS,JA1+IBAS)
     &                         +    F1*RR(M,15)*DENO(JB1+JBAS,ID1+LBAS)
     &                         +    F1*RR(M,16)*DENO(JB1+JBAS,ID2+LBAS)
     &                         +    F2*RR(M,11)*DENO(JB2+JBAS,ID1+LBAS)
     &                         +    F2*RR(M,12)*DENO(JB2+JBAS,ID2+LBAS)
C
              QXCH(IC2+KBAS,JA2+IBAS) = QXCH(IC2+KBAS,JA2+IBAS)
     &                         +    F2*RR(M, 7)*DENO(JB1+JBAS,ID1+LBAS)
     &                         +    F2*RR(M, 8)*DENO(JB1+JBAS,ID2+LBAS)
     &                         +    F1*RR(M, 3)*DENO(JB2+JBAS,ID1+LBAS)
     &                         +    F1*RR(M, 4)*DENO(JB2+JBAS,ID2+LBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SEVENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(7).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QXCH(IB1+JBAS,JC1+KBAS) = QXCH(IB1+JBAS,JC1+KBAS)
     &                         + F1*G1*RR(M,16)*DENO(ID1+LBAS,JA1+IBAS)
     &                         + F2*G1*RR(M, 8)*DENO(ID1+LBAS,JA2+IBAS)
     &                         + F1*G2*RR(M,15)*DENO(ID2+LBAS,JA1+IBAS)
     &                         + F2*G2*RR(M, 7)*DENO(ID2+LBAS,JA2+IBAS)
C
              QXCH(IB1+JBAS,JC2+KBAS) = QXCH(IB1+JBAS,JC2+KBAS)
     &                         + F1*G2*RR(M,14)*DENO(ID1+LBAS,JA1+IBAS)
     &                         + F2*G2*RR(M, 6)*DENO(ID1+LBAS,JA2+IBAS)
     &                         + F1*G1*RR(M,13)*DENO(ID2+LBAS,JA1+IBAS)
     &                         + F2*G1*RR(M, 5)*DENO(ID2+LBAS,JA2+IBAS)
C
              QXCH(IB2+JBAS,JC1+KBAS) = QXCH(IB2+JBAS,JC1+KBAS)
     &                         + F2*G1*RR(M,12)*DENO(ID1+LBAS,JA1+IBAS)
     &                         + F1*G1*RR(M, 4)*DENO(ID1+LBAS,JA2+IBAS)
     &                         + F2*G2*RR(M,11)*DENO(ID2+LBAS,JA1+IBAS)
     &                         + F1*G2*RR(M, 3)*DENO(ID2+LBAS,JA2+IBAS)

C
              QXCH(IB2+JBAS,JC2+KBAS) = QXCH(IB2+JBAS,JC2+KBAS)
     &                         + F2*G2*RR(M,10)*DENO(ID1+LBAS,JA1+IBAS)
     &                         + F1*G2*RR(M, 2)*DENO(ID1+LBAS,JA2+IBAS)
     &                         + F2*G1*RR(M, 9)*DENO(ID2+LBAS,JA1+IBAS)
     &                         + F1*G1*RR(M, 1)*DENO(ID2+LBAS,JA2+IBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       EIGHTH IFLG BATCH (EXCHANGE)
        IF(IFLG(8).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QXCH(IC1+KBAS,JB1+JBAS) = QXCH(IC1+KBAS,JB1+JBAS)
     &                         +       RR(M, 1)*DENO(JA1+IBAS,ID1+LBAS)
     &                         +       RR(M, 2)*DENO(JA1+IBAS,ID2+LBAS)
     &                         +       RR(M, 9)*DENO(JA2+IBAS,ID1+LBAS)
     &                         +       RR(M,10)*DENO(JA2+IBAS,ID2+LBAS)
C
              QXCH(IC1+KBAS,JB2+JBAS) = QXCH(IC1+KBAS,JB2+JBAS)
     &                         +       RR(M, 5)*DENO(JA1+IBAS,ID1+LBAS)
     &                         +       RR(M, 6)*DENO(JA1+IBAS,ID2+LBAS)
     &                         +       RR(M,13)*DENO(JA2+IBAS,ID1+LBAS)
     &                         +       RR(M,14)*DENO(JA2+IBAS,ID2+LBAS)
C
              QXCH(IC2+KBAS,JB1+JBAS) = QXCH(IC2+KBAS,JB1+JBAS)
     &                         +       RR(M, 3)*DENO(JA1+IBAS,ID1+LBAS)
     &                         +       RR(M, 4)*DENO(JA1+IBAS,ID2+LBAS)
     &                         +       RR(M,11)*DENO(JA2+IBAS,ID1+LBAS)
     &                         +       RR(M,12)*DENO(JA2+IBAS,ID2+LBAS)
C
              QXCH(IC2+KBAS,JB2+JBAS) = QXCH(IC2+KBAS,JB2+JBAS)
     &                         +       RR(M, 7)*DENO(JA1+IBAS,ID1+LBAS)
     &                         +       RR(M, 8)*DENO(JA1+IBAS,ID2+LBAS)
     &                         +       RR(M,15)*DENO(JA2+IBAS,ID1+LBAS)
     &                         +       RR(M,16)*DENO(JA2+IBAS,ID2+LBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       NINTH IFLG BATCH (EXCHANGE)
        IF(IFLG(9).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C  
              QXCH(IA1+IBAS,JD1+LBAS) = QXCH(IA1+IBAS,JD1+LBAS)
     &                         +       RR(M, 1)*DENO(IC1+KBAS,JB1+JBAS)
     &                         +       RR(M, 5)*DENO(IC1+KBAS,JB2+JBAS)
     &                         +       RR(M, 3)*DENO(IC2+KBAS,JB1+JBAS)
     &                         +       RR(M, 7)*DENO(IC2+KBAS,JB2+JBAS)
C
              QXCH(IA1+IBAS,JD2+LBAS) = QXCH(IA1+IBAS,JD2+LBAS)
     &                         +       RR(M, 2)*DENO(IC1+KBAS,JB1+JBAS)
     &                         +       RR(M, 6)*DENO(IC1+KBAS,JB2+JBAS)
     &                         +       RR(M, 4)*DENO(IC2+KBAS,JB1+JBAS)
     &                         +       RR(M, 8)*DENO(IC2+KBAS,JB2+JBAS)
C
              QXCH(IA2+IBAS,JD1+LBAS) = QXCH(IA2+IBAS,JD1+LBAS)
     &                         +       RR(M, 9)*DENO(IC1+KBAS,JB1+JBAS)
     &                         +       RR(M,13)*DENO(IC1+KBAS,JB2+JBAS)
     &                         +       RR(M,11)*DENO(IC2+KBAS,JB1+JBAS)
     &                         +       RR(M,15)*DENO(IC2+KBAS,JB2+JBAS)
C
              QXCH(IA2+IBAS,JD2+LBAS) = QXCH(IA2+IBAS,JD2+LBAS)
     &                         +       RR(M,10)*DENO(IC1+KBAS,JB1+JBAS)
     &                         +       RR(M,14)*DENO(IC1+KBAS,JB2+JBAS)
     &                         +       RR(M,12)*DENO(IC2+KBAS,JB1+JBAS)
     &                         +       RR(M,16)*DENO(IC2+KBAS,JB2+JBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       TENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(10).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C
              QXCH(IB1+JBAS,JD1+LBAS) = QXCH(IB1+JBAS,JD1+LBAS)
     &                         +    F1*RR(M,13)*DENO(IC1+KBAS,JA1+IBAS)
     &                         +    F2*RR(M, 5)*DENO(IC1+KBAS,JA2+IBAS)
     &                         +    F1*RR(M,15)*DENO(IC2+KBAS,JA1+IBAS)
     &                         +    F2*RR(M, 7)*DENO(IC2+KBAS,JA2+IBAS)
C
              QXCH(IB1+JBAS,JD2+LBAS) = QXCH(IB1+JBAS,JD2+LBAS)
     &                         +    F1*RR(M,14)*DENO(IC1+KBAS,JA1+IBAS)
     &                         +    F2*RR(M, 6)*DENO(IC1+KBAS,JA2+IBAS)
     &                         +    F1*RR(M,16)*DENO(IC2+KBAS,JA1+IBAS)
     &                         +    F2*RR(M, 8)*DENO(IC2+KBAS,JA2+IBAS)
C
              QXCH(IB2+JBAS,JD1+LBAS) = QXCH(IB2+JBAS,JD1+LBAS)
     &                         +    F2*RR(M, 9)*DENO(IC1+KBAS,JA1+IBAS)
     &                         +    F1*RR(M, 1)*DENO(IC1+KBAS,JA2+IBAS)
     &                         +    F2*RR(M,11)*DENO(IC2+KBAS,JA1+IBAS)
     &                         +    F1*RR(M, 3)*DENO(IC2+KBAS,JA2+IBAS)
C
              QXCH(IB2+JBAS,JD2+LBAS) = QXCH(IB2+JBAS,JD2+LBAS)
     &                         +    F2*RR(M,10)*DENO(IC1+KBAS,JA1+IBAS)
     &                         +    F1*RR(M, 2)*DENO(IC1+KBAS,JA2+IBAS)
     &                         +    F2*RR(M,12)*DENO(IC2+KBAS,JA1+IBAS)
     &                         +    F1*RR(M, 4)*DENO(IC2+KBAS,JA2+IBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       ELEVENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(11).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C
              QXCH(ID1+LBAS,JB1+JBAS) = QXCH(ID1+LBAS,JB1+JBAS)
     &                         +    G1*RR(M, 4)*DENO(JA1+IBAS,IC1+KBAS)
     &                         +    G2*RR(M, 2)*DENO(JA1+IBAS,IC2+KBAS)
     &                         +    G1*RR(M,12)*DENO(JA2+IBAS,IC1+KBAS)
     &                         +    G2*RR(M,10)*DENO(JA2+IBAS,IC2+KBAS)
C
              QXCH(ID1+LBAS,JB2+JBAS) = QXCH(ID1+LBAS,JB2+JBAS)
     &                         +    G1*RR(M, 8)*DENO(JA1+IBAS,IC1+KBAS)
     &                         +    G2*RR(M, 6)*DENO(JA1+IBAS,IC2+KBAS)
     &                         +    G1*RR(M,16)*DENO(JA2+IBAS,IC1+KBAS)
     &                         +    G2*RR(M,14)*DENO(JA2+IBAS,IC2+KBAS)
C
              QXCH(ID2+LBAS,JB1+JBAS) = QXCH(ID2+LBAS,JB1+JBAS)
     &                         +    G2*RR(M, 3)*DENO(JA1+IBAS,IC1+KBAS)
     &                         +    G1*RR(M, 1)*DENO(JA1+IBAS,IC2+KBAS)
     &                         +    G2*RR(M,11)*DENO(JA2+IBAS,IC1+KBAS)
     &                         +    G1*RR(M, 9)*DENO(JA2+IBAS,IC2+KBAS)
C
              QXCH(ID2+LBAS,JB2+JBAS) = QXCH(ID2+LBAS,JB2+JBAS)
     &                         +    G2*RR(M, 7)*DENO(JA1+IBAS,IC1+KBAS)
     &                         +    G1*RR(M, 5)*DENO(JA1+IBAS,IC2+KBAS)
     &                         +    G2*RR(M,15)*DENO(JA2+IBAS,IC1+KBAS)
     &                         +    G1*RR(M,13)*DENO(JA2+IBAS,IC2+KBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C     SKIPPING POINT FOR CLOSED SYSTEMS
399   CONTINUE
C
C     CLOSE ALL LOOPS OVER BASIS FUNCTIONS A,B,C,D
3999  CONTINUE
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
          IF(ICNBAS(I).NE.ICNBAS(J)) GOTO 400
          IF(KQNBAS(I).NE.KQNBAS(J)) GOTO 400
          IF(IABS(MQNBAS(I)).NE.IABS(MQNBAS(J))) GOTO 400
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
C**********************************************************************C
C     SPECIFY OPEN-SHELL COUPLING MATRIX R AND ADD TO COULOMB MATRIX   C
C -------------------------------------------------------------------- C
C               R = [S.D(O).Q + Q.D(O).S]  (RSCF 89)                   C
C**********************************************************************C
C         
      DO J=1,NDIM
        DO K=1,NDIM
          T1(K) = DCMPLX(0.0D0,0.0D0)
          T2(K) = DCMPLX(0.0D0,0.0D0)
          DO L=1,NDIM
            T1(K) = T1(K) + DENT(K,L)*QDIR(L,J) - DENT(K,L)*QDIR(L,J)
            T2(K) = T2(K) + DENT(K,L)*OVAP(L,J)
          ENDDO
        ENDDO
C
        DO I=1,NDIM
          T3(I) = DCMPLX(0.0D0,0.0D0)
          DO K=1,NDIM
            T3(I) = T3(I) + T1(K)*OVAP(I,K) 
     &                    + T2(K)*QDIR(I,K) - T2(K)*QXCH(I,K)
          ENDDO
        ENDDO
C
C       ADD THE PROJECTOR AND Q-MATRIX TO THE COULOMB MATRIX  (RSCF 91)
C       DO I=1,NDIM 
C         FOCK(I,J) = FOCK(I,J) - QMAT(I,J) + T3(I)
C       ENDDO
C       DFNOTE: OKAY... SO I HAVE TO FIX ENERGY TERM
C
        DO I=1,NDIM
          QDIR(I,J) = QDIR(I,J) - T3(I)
          QXCH(I,J) = QXCH(I,J) + T3(I)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INTEGRAL BLOCK AND TIMING COUNTERS                               C
C**********************************************************************C
C
C     ADD ONE-, TWO- AND MANY-CENTER BINS TO (TT|TT) TOTALS
      DO MCNT=1,4
        DO ITT=1,4
          N2ET(5,ITT) = N2ET(5,ITT) + N2ET(MCNT,ITT)
          N2ES(5,ITT) = N2ES(5,ITT) + N2ES(MCNT,ITT)
          T2ES(5,ITT) = T2ES(5,ITT) + T2ES(MCNT,ITT)
        ENDDO
      ENDDO
C
C     ADD ALL (TT|TT) RESULTS TO TOTALS
      DO MCNT=1,5
        DO ITT=1,4
          N2ET(MCNT,7) = N2ET(MCNT,7) + N2ET(MCNT,ITT)
          N2ES(MCNT,7) = N2ES(MCNT,7) + N2ES(MCNT,ITT)
          T2ES(MCNT,7) = T2ES(MCNT,7) + T2ES(MCNT,ITT)       
        ENDDO
      ENDDO
C
      RETURN
      END

