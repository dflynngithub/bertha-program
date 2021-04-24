      SUBROUTINE MP2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                    MM       MM PPPPPPP   222222                      C
C                    MMM     MMM PP    PP 22    22                     C
C                    MMMM   MMMM PP    PP       22                     C
C                    MM MM MM MM PP    PP     22                       C
C                    MM  MMM  MM PPPPPPP    22                         C
C                    MM   M   MM PP       22                           C
C                    MM       MM PP       22222222                     C
C                                                                      C
C -------------------------------------------------------------------- C
C  MP2 EVALUATES SECOND ORDER COULOMB DIRECT/EXCHANGE ENERGIES         C
C  FOR ALL OCCUPIED SOLUTIONS TO A CONVERGED HARTREE-FOCK PROBLEM.     C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=15,MKP=9,MFL=7000000)
C
      CHARACTER*4  HMLTN
      CHARACTER*8  SHAPE
      CHARACTER*16 HMS
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4),LQN(4),ITN(2)
      DIMENSION ISCF(11,6),IFLG(11),ISCR(11)
      DIMENSION INDEX(MCT,-(MKP+1)/2:(MKP+1)/2,MKP),ICNT(4)
      DIMENSION EAB2(NOCC,NOCC,4),EA2(NOCC,4)
C
      COMPLEX*16 ETMP1,ETMP2,ETMP3,ETMP4,RNUMV,RNUMT,RNUMJ,RNUMK
      COMPLEX*16 ABRS(NVRT,NVRT),BARS(NVRT,NVRT),AVR(NVRT),ATR(NVRT)
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 C(MDM,MDM)
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
C
      COMMON/COEF/C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EMDR,EMXC,EUEH
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/GEOM/SHAPE
      COMMON/MAKE/IEAB,IECD,IRIJ(MBS,MBS)
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,IEQS,IERC,IPAR,ICOR,ILEV
      COMMON/QNMS/LABICN(MDM),LABKQN(MDM),LABMQN(MDM)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
C
C     ISCF TELLS WHICH INTEGRALS TO INCLUDE BASED ON OVERLAP COMBINATION
      DATA ISCF/1,1,1,1,1,1,1,1,0,0,0,
     &          1,1,0,0,1,1,1,0,0,0,0,
     &          1,0,1,1,1,0,1,0,0,0,0,
     &          1,1,1,0,1,1,0,0,0,0,0,
     &          1,0,1,0,1,0,0,0,0,0,0,
     &          1,0,0,0,1,0,0,0,0,0,0/
C
C     LINEAR MOLECULE SHORTCUT OPTION
      IF(SHAPE.EQ.'DIATOMIC'.OR.SHAPE.EQ.'LINEAR') THEN
        ILIN = 1
      ELSE
        ILIN = 0
      ENDIF
C
C     IMPORT MBPT1 PAIR RESULTS
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_MP1.dat',STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO IOCCA=1,NOCC
        DO IOCCB=1,NOCC
          READ(8, *) Q1,Q2,Q3,Q4,Q5,EAB2(IOCCA,IOCCB,1)
        ENDDO
      ENDDO
      CLOSE(UNIT=8)
C
C     INITIALISE MP2 ENERGY COUNTERS
      DO N=2,4
        DO IOCCA=1,NOCC
          EA2(IOCCA,N) = 0.0D0
          DO IOCCB=1,NOCC
            EAB2(IOCCA,IOCCB,N) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     INDEXING ROUTINE: SO WE CAN SET BASIS FUNCTION LABELS BY BLOCK
      ICOUNT = 0
C
C     LOOP OVER ATOMIC CENTERS
      DO ICT=1,NCNT
C
C       LOOP OVER KAPPA VALUES FOR THIS NUCLEAR CENTER
        DO KN=1,NKAP(ICT)
C
C         IMPORT KAPPA, MAXIMUM MQN
          KAPPA = KVALS(KN,ICT)
          MJMAX  = 2*IABS(KAPPA)-1
C
C         LOOP OVER MQN VALUES AND RECORD INDEX
          DO MJ=1,MJMAX,2
            ICOUNT = ICOUNT+1
            INDEX(ICT,KAPPA,MJ) = ICOUNT
          ENDDO
        ENDDO
      ENDDO
C
C     SET THE LIMITS FOR RUNNING OVER DENSITY COMBINATIONS
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
20    FORMAT(1X,A,11X,A,11X,A,11X,A,11X,A)
21    FORMAT(' (',I2,',',I2,')',3X,F13.7,5X,F11.7,5X,F11.7,4X,F13.7)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',24),'MP2 pair-wise summary'
      WRITE(7, *) REPEAT(' ',24),'MP2 pair-wise summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) '( a, b)','E2(H)','E2(J)','E2(K)','E2(ab)'
      WRITE(7,20) '( a, b)','E2(H)','E2(J)','E2(K)','E2(ab)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C**********************************************************************C
C     LOOP OVER UNIQUE OCCUPIED ORBITALS (A,B)       (INDEX 1000)      C
C**********************************************************************C
C
      CALL CPU_TIME(TDM1)
      DO 1000 IOCCA=1,NOCC
      DO 1000 IOCCB=1,IOCCA
C
      IA = IOCCA + NSHIFT
      IB = IOCCB + NSHIFT
C
C     INITIALISE STORAGE ARRAYS FOR (AR|BS) AND (AR|SB)
      DO IR=1,NVRT
        DO IS=1,NVRT
          ABRS(IR,IS) = DCMPLX(0.0D0,0.0D0)
          BARS(IR,IS) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     LOOP OVER UNIQUE VIRTUAL ORBITALS (R,S)        (INDEX 2000)      C
C**********************************************************************C
C
C     GOTO 2001
      DO 2000 IR=1,NVRT
      DO 2000 IS=1,NVRT
C
        IOCCR = NSHIFT + NOCC + IR
        IOCCS = NSHIFT + NOCC + IS
C
C       INITIALISE COULOMB DIRECT AND EXCHANGE ARRAYS
        DO J=1,NDIM
          DO I=1,NDIM
            GDIR(I,J) = DCMPLX(0.0D0,0.0D0)
            GXCH(I,J) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
C
C**********************************************************************C
C     LOOP OVER SPINORS A, B BY BLOCK                (INDEX 3000)      C
C**********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR A AND B: TT = LL (1) or SS (4)
      DO 3000 IT1=ITSTRT,ITSTOP,ITSKIP
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
C     LOOP OVER CENTER A
      DO 3000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 3000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
      DO 3000 KA=1,NKAP(ICNTA)
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
      DO 3000 KB=1,NKAP(ICNTB)
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
      DO 3000 MA=1,IABS(KQN(1))
        MJA    = 2*MA-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 3000 MB=1,IABS(KQN(2))
        MJB    = 2*MB-1
        MQN(2) = MJB
C
C     CALCULATE NEW BLOCK OF E(AB) COEFFS AT NEXT OPPORTUNITY
      IEAB = 1
C
C**********************************************************************C
C     LOOP OVER SPINORS C, D BY BLOCK                (INDEX 4000)      C
C**********************************************************************C
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
C     LOOP OVER CENTER C
      DO 4000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER C
        XYZ(1,3) = COORD(1,ICNTC)
        XYZ(2,3) = COORD(2,ICNTC)
        XYZ(3,3) = COORD(3,ICNTC)
C
C     LOOP OVER CENTER D
      DO 4000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER D
        XYZ(1,4) = COORD(1,ICNTD)
        XYZ(2,4) = COORD(2,ICNTD)
        XYZ(3,4) = COORD(3,ICNTD)
C
C     LOOP OVER KQN(C) VALUES
      DO 4000 KC=1,NKAP(ICNTC)
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
      DO 4000 KD=1,NKAP(ICNTD)
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
      DO 4000 MC=1,IABS(KQN(3))
        MJC    = 2*MC-1
        MQN(3) = MJC
C
C     LOOP OVER |MQN(D)| VALUES
      DO 4000 MD=1,IABS(KQN(4))
        MJD    = 2*MD-1
        MQN(4) = MJD
C
C     CALCULATE NEW BLOCK OF E(CD) COEFFS AT NEXT OPPORTUNITY
      IECD = 1
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
        GOTO 4000
      ENDIF
5503  CONTINUE
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
      F1 = PKAB*PMAB1
      F2 = PKAB*PMAB2
      G1 = PKCD*PMCD1
      G2 = PKCD*PMCD2
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
        GO TO 4000
      ENDIF
C
C     READ IN FLAG VALUES FROM ISCF DATA BLOCK
      DO M=1,11
        IFLG(M) = ISCF(M,ITSCF)
      ENDDO
C
C     INCLUDE SPECIAL CASES
      IF(ITSCF.EQ.1.AND.IQ1.EQ.IQ3) IFLG(10) = 1
      IF(ITSCF.EQ.1.AND.IQ2.EQ.IQ3) IFLG( 9) = 1
      IF(ITSCF.EQ.3.AND.IQ2.EQ.IQ3) IFLG( 9) = 1
      IF(ITSCF.EQ.4.AND.IQ2.EQ.IQ3) IFLG( 9) = 1
      IF(ITSCF.EQ.1.AND.IQ2.EQ.IQ4) IFLG(11) = 1
C
C**********************************************************************C
C     FOURTH LAYER OF LOOPS, OVER BASIS FUNCTIONS A AND B (5000)       C
C -------------------------------------------------------------------- C
C     DEPENDING ON CURRENT COMBINATION OF MQN VALUES, TAKE ADVANTAGE OF
C     QUANTUM NUMBER PERMUTATION SYMMETRY (HENCE MINIMISE CALLS TO ERI):
C              ( MA, MB|-MD,-MC) =     PCD*( MA, MB| MC, MD)
C              (-MB,-MA| MC, MD) = PAB*    ( MA, MB| MC, MD)
C              (-MA,-MB|-MC,-MD) = PAB*PCD*( MA, MB| MC, MD)*
C              ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
C -------------------------------------------------------------------- C
C     (RSCF 86, 87)                                                    C
C**********************************************************************C
C
C     LOOP OVER ELEMENTS OF FOCK MATRIX BLOCK
      DO 5000 IBAS=1,NBAS(1)
      DO 5000 JBAS=1,NBAS(2)
C
C       GENERATE BATCH OF ELECTRON REPULSION INTEGRALS
        CALL ERI(RR,XYZ,KQN,MQN,EXPT,NBAS,ITN,IBAS,JBAS)
C
C     DEPENDING ON CURRENT COMBINATION OF MQN VALUES, TAKE ADVANTAGE OF
C     QUANTUM NUMBER PERMUTATION SYMMETRY (HENCE MINIMISE CALLS TO ERI):
C              ( MA, MB|-MD,-MC) =     PCD*( MA, MB| MC, MD)
C              (-MB,-MA| MC, MD) = PAB*    ( MA, MB| MC, MD)
C              (-MA,-MB|-MC,-MD) = PAB*PCD*( MA, MB| MC, MD)*
C              ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
C
C       1ST CASE (DIRECT)
        IF(IFLG(1).EQ.1) THEN
          M = 0
          DO KBAS=1,NBAS(3)
            DO LBAS=1,NBAS(4)
              M = M+1
C
              GDIR(IA1+IBAS,JB1+JBAS) = GDIR(IA1+IBAS,JB1+JBAS)
     &          +       RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IR)
     &          +       RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IR)
C
              GDIR(IA1+IBAS,JB2+JBAS) = GDIR(IA1+IBAS,JB2+JBAS)
     &          +       RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IR)
     &          +       RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IR)
C
              GDIR(IA2+IBAS,JB1+JBAS) = GDIR(IA2+IBAS,JB1+JBAS)
     &          +       RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IR)
     &          +       RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IR)
C
              GDIR(IA2+IBAS,JB2+JBAS) = GDIR(IA2+IBAS,JB2+JBAS)
     &          +       RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IR)
     &          +       RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       2ND CASE (DIRECT)
        IF(IFLG(2).EQ.1) THEN
          M = 0
          DO KBAS=1,NBAS(3)
            DO LBAS=1,NBAS(4)
              M = M+1
C
              GDIR(IA1+IBAS,JB1+JBAS) = GDIR(IA1+IBAS,JB1+JBAS)
     &          +    G1*RR(M, 1)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M, 2)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M, 3)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IR)
     &          +    G1*RR(M, 4)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IR)
C
              GDIR(IA1+IBAS,JB2+JBAS) = GDIR(IA1+IBAS,JB2+JBAS)
     &          +    G1*RR(M, 5)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M, 6)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M, 7)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IR)
     &          +    G1*RR(M, 8)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IR)
C
              GDIR(IA2+IBAS,JB1+JBAS) = GDIR(IA2+IBAS,JB1+JBAS)
     &          +    G1*RR(M, 9)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M,10)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M,11)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IR)
     &          +    G1*RR(M,12)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IR)
C
              GDIR(IA2+IBAS,JB2+JBAS) = GDIR(IA2+IBAS,JB2+JBAS)
     &          +    G1*RR(M,13)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M,14)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M,15)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IR)
     &          +    G1*RR(M,16)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       3RD CASE (DIRECT)
        IF(IFLG(3).EQ.1) THEN
          M = 0
          DO KBAS=1,NBAS(3)
            DO LBAS=1,NBAS(4)
              M = M+1
C
              GDIR(IC1+KBAS,JD1+LBAS) = GDIR(IC1+KBAS,JD1+LBAS)
     &          +       RR(M, 1)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 5)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M, 9)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,13)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IR)
C
              GDIR(IC1+KBAS,JD2+LBAS) = GDIR(IC1+KBAS,JD2+LBAS)
     &          +       RR(M, 2)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 6)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M,10)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,14)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IR)
C
              GDIR(IC2+KBAS,JD1+LBAS) = GDIR(IC2+KBAS,JD1+LBAS)
     &          +       RR(M, 3)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 7)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M,11)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,15)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IR)
C
              GDIR(IC2+KBAS,JD2+LBAS) = GDIR(IC2+KBAS,JD2+LBAS)
     &          +       RR(M, 4)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 8)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M,12)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,16)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       4TH CASE (DIRECT)
        IF(IFLG(4).EQ.1) THEN
          M = 0
          DO KBAS=1,NBAS(3)
            DO LBAS=1,NBAS(4)
              M = M+1
C
              GDIR(IC1+KBAS,JD1+LBAS) = GDIR(IC1+KBAS,JD1+LBAS)
     &          +    F1*RR(M, 1)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M, 5)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M, 9)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IR)
     &          +    F1*RR(M,13)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IR)
C
              GDIR(IC1+KBAS,JD2+LBAS) = GDIR(IC1+KBAS,JD2+LBAS)
     &          +    F1*RR(M, 2)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M, 6)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M,10)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IR)
     &          +    F1*RR(M,14)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IR)
C
              GDIR(IC2+KBAS,JD1+LBAS) = GDIR(IC2+KBAS,JD1+LBAS)
     &          +    F1*RR(M, 3)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M, 7)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M,11)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IR)
     &          +    F1*RR(M,15)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IR)
C
              GDIR(IC2+KBAS,JD2+LBAS) = GDIR(IC2+KBAS,JD2+LBAS)
     &          +    F1*RR(M, 4)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M, 8)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M,12)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IR)
     &          +    F1*RR(M,16)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       5TH CASE (EXCHANGE)
        IF(IFLG(5).EQ.1) THEN
          DO LBAS=1,NBAS(4)
            DO KBAS=1,NBAS(3)
              M = (KBAS-1)*NBAS(4)+LBAS
C
              GXCH(IA1+IBAS,JD1+LBAS) = GXCH(IA1+IBAS,JD1+LBAS)
     &          +       RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IR)
C
              GXCH(IA1+IBAS,JD2+LBAS) = GXCH(IA1+IBAS,JD2+LBAS)
     &          +       RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IR)
C
              GXCH(IA2+IBAS,JD1+LBAS) = GXCH(IA2+IBAS,JD1+LBAS)
     &          +       RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IR)
C
              GXCH(IA2+IBAS,JD2+LBAS) = GXCH(IA2+IBAS,JD2+LBAS)
     &          +       RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       6TH CASE (EXCHANGE)
        IF(IFLG(6).EQ.1) THEN
          M = 0
          DO KBAS=1,NBAS(3)
            DO LBAS=1,NBAS(4)
              M = M+1
C
              GXCH(IA1+IBAS,JC1+KBAS) = GXCH(IA1+IBAS,JC1+KBAS)
     &          +    G2*RR(M, 3)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G1*RR(M, 4)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G2*RR(M, 7)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IR)
     &          +    G1*RR(M, 8)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IR)
C
              GXCH(IA1+IBAS,JC2+KBAS) = GXCH(IA1+IBAS,JC2+KBAS)
     &          +    G1*RR(M, 1)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G2*RR(M, 2)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G1*RR(M, 5)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IR)
     &          +    G2*RR(M, 6)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IR)
C
              GXCH(IA2+IBAS,JC1+KBAS) = GXCH(IA2+IBAS,JC1+KBAS)
     &          +    G2*RR(M,11)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G1*RR(M,12)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G2*RR(M,15)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IR)
     &          +    G1*RR(M,16)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IR)
C
              GXCH(IA2+IBAS,JC2+KBAS) = GXCH(IA2+IBAS,JC2+KBAS)
     &          +    G1*RR(M, 9)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G2*RR(M,10)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G1*RR(M,13)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IR)
     &          +    G2*RR(M,14)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       7TH CASE (EXCHANGE)
        IF(IFLG(7).EQ.1) THEN
          DO LBAS=1,NBAS(4)
            DO KBAS=1,NBAS(3)
              M = (KBAS-1)*NBAS(4)+LBAS
C
              GXCH(IB1+JBAS,JD1+LBAS) = GXCH(IB1+JBAS,JD1+LBAS)
     &          +    F2*RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F2*RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F1*RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IR)
     &          +    F1*RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IR)
C
              GXCH(IB1+JBAS,JD2+LBAS) = GXCH(IB1+JBAS,JD2+LBAS)
     &          +    F2*RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F2*RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F1*RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IR)
     &          +    F1*RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IR)
C
              GXCH(IB2+JBAS,JD1+LBAS) = GXCH(IB2+JBAS,JD1+LBAS)
     &          +    F1*RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F1*RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F2*RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IR)
     &          +    F2*RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IR)
C
              GXCH(IB2+JBAS,JD2+LBAS) = GXCH(IB2+JBAS,JD2+LBAS)
     &          +    F1*RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F1*RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F2*RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IR)
     &          +    F2*RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       8TH CASE (EXCHANGE)
        IF(IFLG(8).EQ.1) THEN
          M = 0
          DO KBAS=1,NBAS(3)
            DO LBAS=1,NBAS(4)
              M = M+1
C
              GXCH(IB1+JBAS,JC1+KBAS) = GXCH(IB1+JBAS,JC1+KBAS)
     &          + F2*G2*RR(M, 7)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F2*G1*RR(M, 8)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F1*G2*RR(M,15)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IR)
     &          + F1*G1*RR(M,16)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IR)
C
              GXCH(IB1+JBAS,JC2+KBAS) = GXCH(IB1+JBAS,JC2+KBAS)
     &          + F2*G1*RR(M, 5)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F2*G2*RR(M, 6)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F1*G1*RR(M,13)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IR)
     &          + F1*G2*RR(M,14)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IR)
C
              GXCH(IB2+JBAS,JC1+KBAS) = GXCH(IB2+JBAS,JC1+KBAS)
     &          + F1*G2*RR(M, 3)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F1*G1*RR(M, 4)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F2*G2*RR(M,11)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IR)
     &          + F2*G1*RR(M,12)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IR)
C
              GXCH(IB2+JBAS,JC2+KBAS) = GXCH(IB2+JBAS,JC2+KBAS)
     &          + F1*G1*RR(M, 1)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F1*G2*RR(M, 2)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F2*G1*RR(M, 9)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IR)
     &          + F2*G2*RR(M,10)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       9TH CASE (EXCHANGE)
        IF(IFLG(9).EQ.1) THEN
          M = 0
          DO KBAS=1,NBAS(3)
            DO LBAS=1,NBAS(4)
              M = M+1
C
              GXCH(IC1+KBAS,JB1+JBAS) = GXCH(IC1+KBAS,JB1+JBAS)
     &          +       RR(M, 1)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M, 2)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IR)
     &          +       RR(M, 9)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M,10)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IR)
C
              GXCH(IC1+KBAS,JB2+JBAS) = GXCH(IC1+KBAS,JB2+JBAS)
     &          +       RR(M, 5)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M, 6)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IR)
     &          +       RR(M,13)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M,14)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IR)
C
              GXCH(IC2+KBAS,JB1+JBAS) = GXCH(IC2+KBAS,JB1+JBAS)
     &          +       RR(M, 3)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M, 4)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IR)
     &          +       RR(M,11)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M,12)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IR)
C
              GXCH(IC2+KBAS,JB2+JBAS) = GXCH(IC2+KBAS,JB2+JBAS)
     &          +       RR(M, 7)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M, 8)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IR)
     &          +       RR(M,15)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M,16)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       10TH CASE (EXCHANGE)
        IF(IFLG(10).EQ.1) THEN
          M = 0
          DO KBAS=1,NBAS(3)
            DO LBAS=1,NBAS(4)
              M = M+1
C
              GXCH(IC1+KBAS,JA1+IBAS) = GXCH(IC1+KBAS,JA1+IBAS)
     &          +    F2*RR(M, 9)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F2*RR(M,10)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IR)
     &          +    F1*RR(M,13)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F1*RR(M,14)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IR)
C
              GXCH(IC1+KBAS,JA2+IBAS) = GXCH(IC1+KBAS,JA2+IBAS)
     &          +    F1*RR(M, 1)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F1*RR(M, 2)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IR)
     &          +    F2*RR(M, 5)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F2*RR(M, 6)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IR)
C
              GXCH(IC2+KBAS,JA1+IBAS) = GXCH(IC2+KBAS,JA1+IBAS)
     &          +    F2*RR(M,11)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F2*RR(M,12)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IR)
     &          +    F1*RR(M,15)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F1*RR(M,16)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IR)
C
              GXCH(IC2+KBAS,JA2+IBAS) = GXCH(IC2+KBAS,JA2+IBAS)
     &          +    F1*RR(M, 3)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F1*RR(M, 4)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IR)
     &          +    F2*RR(M, 7)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F2*RR(M, 8)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       11TH CASE (EXCHANGE)
        IF(IFLG(11).EQ.1) THEN
          DO LBAS=1,NBAS(4)
            DO KBAS=1,NBAS(3)
              M = (KBAS-1)*NBAS(4)+LBAS
C
              GXCH(ID1+LBAS,JB1+JBAS) = GXCH(ID1+LBAS,JB1+JBAS)
     &          +    G2*RR(M, 2)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G1*RR(M, 4)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IR)
     &          +    G2*RR(M,10)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G1*RR(M,12)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IR)
C
              GXCH(ID1+LBAS,JB2+JBAS) = GXCH(ID1+LBAS,JB2+JBAS)
     &          +    G2*RR(M, 6)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G1*RR(M, 8)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IR)
     &          +    G2*RR(M,14)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G1*RR(M,16)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IR)
C
              GXCH(ID2+LBAS,JB1+JBAS) = GXCH(ID2+LBAS,JB1+JBAS)
     &          +    G1*RR(M, 1)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M, 3)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IR)
     &          +    G1*RR(M, 9)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M,11)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IR)
C
              GXCH(ID2+LBAS,JB2+JBAS) = GXCH(ID2+LBAS,JB2+JBAS)
     &          +    G1*RR(M, 5)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M, 7)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IR)
     &          +    G1*RR(M,13)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M,15)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C     END LOOP OVER IAOCC,IBOCC BLOCK ADDRESSES
5000  CONTINUE
C     END LOOPS OVER A,B,C,D OVERLAP BLOCKS
4000  CONTINUE
3000  CONTINUE
C
C     COMPLETE CONSTRUCTION OF GDIR AND GXCH BY MATRIX CONJUGATION
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
C
C         IF HMLTN = 'NORL' SKIP THE NEXT FEW CALCULATIONS
          IF(HMLTN.EQ.'NORL') GOTO 401
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF LS BLOCK
          GDIR(I,L) = GDIR(I,L) + DCONJG(GDIR(L,I))
          GDIR(L,I) =             DCONJG(GDIR(I,L))
          GXCH(I,L) = GXCH(I,L) + DCONJG(GXCH(L,I))
          GXCH(L,I) =             DCONJG(GXCH(I,L))
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF SL BLOCK
          GDIR(K,J) = GDIR(K,J) + DCONJG(GDIR(J,K))
          GDIR(J,K) =             DCONJG(GDIR(K,J))
          GXCH(K,J) = GXCH(K,J) + DCONJG(GXCH(J,K))
          GXCH(J,K) =             DCONJG(GXCH(K,J))
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF SS BLOCK
          GDIR(K,L) = GDIR(K,L) + DCONJG(GDIR(L,K))
          GDIR(L,K) =             DCONJG(GDIR(K,L))
          GXCH(K,L) = GXCH(K,L) + DCONJG(GXCH(L,K))
          GXCH(L,K) =             DCONJG(GXCH(K,L))
C
401       CONTINUE
        ENDDO
      ENDDO
C
C     END LOOP OVER UNIQUE VIRTUAL ORBITAL COMBINATIONS (R,S)
2000  CONTINUE
C2001  CONTINUE
C
C**********************************************************************C
C     CALCULATE INTERACTION ENERGY FROM ORBITAL PAIR (A,B)             C
C**********************************************************************C
C
C     ALL MATRIX ELEMENTS (AR|BS) AND (AR|SB) ARE STORED FOR SUMMATION
      ETMP3 = DCMPLX(0.0D0,0.0D0)
      ETMP4 = DCMPLX(0.0D0,0.0D0)
      DO IOCCR=1,NVRT
        DO IOCCS=1,NVRT
C
          IR = NSHIFT + NOCC + IOCCR
          IS = NSHIFT + NOCC + IOCCS
C
          RNUMJ = DCONJG(ABRS(IOCCR,IOCCS))*ABRS(IOCCR,IOCCS)
          RNUMK = DCONJG(BARS(IOCCR,IOCCS))*ABRS(IOCCR,IOCCS)
          RDEN  = EIGEN(IA)+EIGEN(IB)-EIGEN(IR)-EIGEN(IS)

          ETMP3 = ETMP3 + RNUMJ/RDEN
          ETMP4 = ETMP4 + RNUMK/RDEN
        ENDDO
      ENDDO
      EAB2(IOCCA,IOCCB,2) = DREAL(ETMP3)
      EAB2(IOCCA,IOCCB,3) = DREAL(ETMP4)
      EAB2(IOCCA,IOCCB,4) = DREAL(ETMP3 - ETMP4)
C
C     OUTPUT ENERGIES TO TERMINAL
      WRITE(6,21) IOCCA,IOCCB,(EAB2(IOCCA,IOCCB,N),N=1,4)
      WRITE(7,21) IOCCA,IOCCB,(EAB2(IOCCA,IOCCB,N),N=1,4)
C
      IF(IOCCA.NE.IOCCB) THEN
        DO N=1,4
          EAB2(IOCCB,IOCCA,N) = EAB2(IOCCA,IOCCB,N)
        ENDDO
      ENDIF
C
C     END LOOP OVER UNIQUE OCCUPIED ORBITAL COMBINATIONS (A,B)
1000  CONTINUE
C
C**********************************************************************C
C     END OF CALCULATION OVER OCCUPIED ORBITALS (A,B)                  C
C**********************************************************************C
C
C     RECORD TIME AT THE END OF MAIN CALCULATION
      CALL CPU_TIME(TDM2)
C
C     WRITE RESULTS OF MP2 ENERGIES TO AN EXTERNAL FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_MP2.dat',STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO IOCCA=1,NOCC
        DO IOCCB=1,NOCC
          WRITE(8, *) (EAB2(IOCCA,IOCCB,N),N=1,4)
        ENDDO
      ENDDO
      CLOSE(UNIT=8)
C
C     CALCULATE MBPT1 SINGLE-PARTICLE ENERGIES AND MOLECULAR TOTALS
      EDIR2 = 0.0D0
      EXCH2 = 0.0D0
      DO IOCCA=1,NOCC
        DO N=1,3
          EA2(IOCCA,N) = 0.0D0
          DO IOCCB=1,NOCC
            EA2(IOCCA,N) = EA2(IOCCA,N) + EAB2(IOCCA,IOCCB,N)
          ENDDO
        ENDDO
        EA2(IOCCA,4) = EA2(IOCCA,2) - EA2(IOCCA,3)
        EDIR2 = EDIR2 + EA2(IOCCA,2)*0.5D0
        EXCH2 = EXCH2 + EA2(IOCCA,3)*0.5D0
      ENDDO
      ETOT2 = EDIR2 - EXCH2
C
C     ORBITAL SUMMARIES
30    FORMAT(1X,A,11X,A,11X,A,11X,A,11X,A)
31    FORMAT('  ',I2,'    ',3X,F13.7,5X,F11.7,5X,F11.7,4X,F13.7)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',21),'MP2 single particle summary'
      WRITE(7, *) REPEAT(' ',21),'MP2 single particle summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,30) '  a    ','E1(a)','E2(J)','E2(K)',' E2(a)'
      WRITE(7,30) '  a    ','E1(a)','E2(J)','E2(K)',' E2(a)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO IOCCA=1,NOCC
        WRITE(6,31) IOCCA,(EA2(IOCCA,N),N=1,4)
        WRITE(7,31) IOCCA,(EA2(IOCCA,N),N=1,4)
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     TOTAL ENERGIES
32    FORMAT(1X,A,4X,A,2X,'=',10X,F18.8,' au')
33    FORMAT(1X,A,5X,'=',23X,A)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',20),'MP2 molecular summary'
      WRITE(7, *) REPEAT(' ',20),'MP2 molecular summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,32) 'Correlation Coulomb direct   ','E2(J)',EDIR2
      WRITE(7,32) 'Correlation Coulomb direct   ','E2(J)',EDIR2
      WRITE(6,32) 'Correlation Coulomb exchange ','E2(K)',EXCH2
      WRITE(7,32) 'Correlation Coulomb exchange ','E2(K)',EXCH2
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,32) 'Correlation Coulomb total    ','E2(G)',EDIR2-EXCH2
      WRITE(7,32) 'Correlation Coulomb total    ','E2(G)',EDIR2-EXCH2
      WRITE(6,32) 'Hartree-Fock molecular energy','E2   ',ETOT+ETOT2
      WRITE(7,32) 'Hartree-Fock molecular energy','E2   ',ETOT+ETOT2
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,33) 'MP2 time                   ',HMS(TDM2-TDM1)
      WRITE(7,33) 'MP2 time                   ',HMS(TDM2-TDM1)
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
      RETURN
      END

