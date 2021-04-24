      SUBROUTINE MBPT1
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
C    MBPT1 EVALUATES ZERO AND FIRST ORDER ENERGIES FOR ALL OCCUPIED    C
C    SOLUTIONS TO A CONVERGED HARTREE-FOCK PROBLEM.                    C
C**********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=MKP,MBS=26,MB2=MBS*MBS,
     &                                        MLL=MKP*(MKP+1)*(MKP+2)/6)
C
      CHARACTER*4 HMLTN
      CHARACTER*15 TIMEHMS
      CHARACTER*40 FILNAM,STRING      
C
      COMPLEX*16 C(MDM,MDM),RR(MB2,16),ETMP1,ETMP2,ETMP3,ETMP4
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           GDIR(MDM,MDM),GXCH(MDM,MDM),QDIR(MDM,MDM),
     &           QXCH(MDM,MDM),BDIR(MDM,MDM),BXCH(MDM,MDM),FOCK(MDM,MDM)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NFUNS(4),LQN(4)
      DIMENSION ITQN(2),IFLG(11),ISCF(11,6)
      DIMENSION INDEX(MCT,-MKP:MKP,2*(MKP+1)*MKP)
      DIMENSION T1(MDM),T2(MDM),T3(MDM)
C
      COMMON/COEF/C
      COMMON/ILLM/ILLAD(MCT,MCT,MKP,MKP,MKP,MKP),IABLL,ICDLL
      COMMON/ISSM/ISSAD(MCT,MCT,MKP,MKP,MKP,MKP),IABSS,ICDSS
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EMDR,EMXC
      COMMON/FLNM/STRING,LN
      COMMON/MTRX/OVAP,HNUC,HKIN,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITER,IALL,IRUN,IEQS
      COMMON/QNMS/ICNLAB(MDM),KQNLAB(MDM),MQNLAB(MDM)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR,IOCCM0
C
      DIMENSION EAB1(NOCC,NOCC,6),EA1(NOCC,6)
C
C     ISCF TELLS WHICH INTEGRALS TO INCLUDE BASED ON OVERLAP COMBINATION
      DATA ISCF/1,0,1,0,1,0,1,0,1,1,0,
     &          1,0,0,0,1,0,0,0,1,1,0,
     &          0,1,1,0,0,0,0,0,1,1,0,
     &          1,0,0,1,1,0,0,0,1,0,0,
     &          0,1,0,1,0,0,0,0,1,0,0,
     &          0,1,0,0,0,0,0,0,1,0,0/
C
C
C     TABLE HEADINGS FOR DISPLAY OF RESULTS
10    FORMAT(1X,A,9X,A,8X,A,8X,A,9X,A)
11    FORMAT(' (',I2,',',I2,')',1X,F13.7,2X,F11.7,2X,F11.7,2X,F13.7)
12    FORMAT('  ',I2,'    ',1X,F13.7,2X,F11.7,2X,F11.7,2X,F13.7)
84    FORMAT(1X,A,4X,A,2X,'=',F18.8,' au')
87    FORMAT(1X,A,5X,'=',12X,A)
C
      WRITE(6, *) 
      WRITE(7, *) 
      WRITE(6, *) REPEAT(' ',19),'MBPT1 pair-wise summary'
      WRITE(7, *) REPEAT(' ',19),'MBPT1 pair-wise summary'
      WRITE(6, *) REPEAT('-',62)
      WRITE(7, *) REPEAT('-',62)
      WRITE(6,10) '( a, b)','E1(H)','E1(J)','E1(K)','E1(ab)'
      WRITE(7,10) '( a, b)','E1(H)','E1(J)','E1(K)','E1(ab)'
      WRITE(6, *) REPEAT('-',62)
      WRITE(7, *) REPEAT('-',62)
C
      CALL CPU_TIME(TBGN)
C
C     INITIALISE ENERGY COUNTERS
      DO N=1,6
        DO IOCCB=1,NOCC
          EA1(IOCCB,N) = 0.0D0
          DO IOCCA=1,NOCC
            EAB1(IOCCA,IOCCB,N) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     CALCULATE ONE-BODY MATRIX REPS
      CALL ONEEL
C
C     INDEXING ROUTINE: SO WE CAN SET BASIS FUNCTION LABELS BY BLOCK
      ICOUNT = 0
C
C     LOOP OVER NUCLEAR CENTRES
      DO ICNT=1,NCNT
C
C       LOOP OVER KAPPA VALUES FOR THIS NUCLEAR CENTRE
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
C     SET THE LIMITS FOR RUNNING OVER DENSITY COMBINATIONS
      IF(HMLTN.EQ.'NORL') THEN
        ITSTRT = 1
        ITSTOP = 1
        ITSKIP = 1
      ELSEIF(HMLTN.EQ.'DHFR'.OR.HMLTN.EQ.'DHFP'.OR.HMLTN.EQ.'DHFB') THEN
        ITSTRT = 1
        ITSTOP = 4
        ITSKIP = 3
      ENDIF
C
C**********************************************************************C
C     ONE-BODY ENERGIES                                                C
C**********************************************************************C
C
      DO 1000 IOCCA=1,NOCC
      DO 1000 IOCCB=1,IOCCA
C
        IA = IOCCA + NSHIFT
        IB = IOCCB + NSHIFT
C
C       ONE-BODY ENERGY
        ETMP1 = DCMPLX(0.0D0,0.0D0)
        ETMP2 = DCMPLX(0.0D0,0.0D0)
        IF(IOCCA.NE.IOCCB) GOTO 90
        DO J=1,NDIM
          DO I=1,NDIM
            ETMP1 = ETMP1 + DCONJG(C(I,IA))*C(J,IA)*HNUC(I,J)
            ETMP2 = ETMP2 + DCONJG(C(I,IA))*C(J,IA)*HKIN(I,J)
          ENDDO
        ENDDO
90      CONTINUE
        EAB1(IOCCA,IOCCB,1) = DREAL(ETMP1)
        EAB1(IOCCA,IOCCB,2) = DREAL(ETMP2)
        EAB1(IOCCA,IOCCB,3) = DREAL(ETMP1 + ETMP2)
C
C       INITIALISE COULOMB DIRECT AND EXCHANGE ARRAYS
        DO J=1,NDIM
          DO I=1,NDIM
            GDIR(I,J) = DCMPLX(0.0D0,0.0D0)
            GXCH(I,J) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
cC
c100   CONTINUE
C
C**********************************************************************C
C     LOOP OVER UNIQUE OCCUPIED ORBITALS (A,B)       (INDEX 1000)      C
C**********************************************************************C
C
c      DO 1000 IOCCA=1,NOCC
c      DO 1000 IOCCB=1,IOCCA
cC
c        IA = IOCCA + NSHIFT
c        IB = IOCCB + NSHIFT
C
C**********************************************************************C
C     LOOP OVER SPINORS A, B BY BLOCK                (INDEX 3000)      C
C**********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR A AND B: T = {L} OR {L,S}
      DO 3000 IT1=ITSTRT,ITSTOP,ITSKIP
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
C     LOOP OVER CENTRE A
      DO 3000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 3000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTRE B
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
        NFUNS(1) = NFUNCT(LQN(1)+1,ICNTA)
        DO IBAS=1,NFUNS(1)
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
        NFUNS(2) = NFUNCT(LQN(2)+1,ICNTB)
        DO JBAS=1,NFUNS(2)
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
      IABLL = ILLAD(ICNTA,ICNTB,KA,KB,MA,MB)
      IABSS = ISSAD(ICNTA,ICNTB,KA,KB,MA,MB)
C
C**********************************************************************C
C     FOR THIS CHOICE OF A AND B, COMPUTE ADDRESSES AND PHASES         C
C**********************************************************************C
C
C     CALCULATE BLOCK INDICES FOR {AB} COMBINATIONS
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
C
      IQ12 = (IQ1*(IQ1-1))/2 + IQ2
C
C     FURTHER DEFINE STARTING ADDRESSES FOR {ABCD} BASIS LABELS
      IA1 = LARGE(ICNTA,KA,MJA  ) + NADDAB
      IB1 = LARGE(ICNTB,KB,MJB  ) + NADDAB
C
      IA2 = LARGE(ICNTA,KA,MJA+1) + NADDAB
      IB2 = LARGE(ICNTB,KB,MJB+1) + NADDAB
C
      JA1 = LARGE(ICNTA,KA,MJA  ) + NADDAB
      JB1 = LARGE(ICNTB,KB,MJB  ) + NADDAB
C
      JA2 = LARGE(ICNTA,KA,MJA+1) + NADDAB
      JB2 = LARGE(ICNTB,KB,MJB+1) + NADDAB
C
C     CALCULATE KQN PHASE FACTORS FOR PERMUTING INTEGRALS
      IF((KQN(1)*KQN(2)).GT.0) THEN 
        PKAB = 1.0D0
      ELSE
        PKAB =-1.0D0
      ENDIF
C
C     CALCULATE MQN PHASE FACTORS FOR PERMUTING INTEGRALS
      PMAB1 = DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PMAB2 = DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
C
      F1 = PKAB*PMAB1
      F2 = PKAB*PMAB2
C
C**********************************************************************C
C     LOOP OVER SPINORS C, D BY BLOCK                (INDEX 4000)      C
C**********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR C AND D: T = {L} OR {L,S}
      DO 4000 IT2=ITSTRT,ITSTOP,ITSKIP
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
C     LOOP OVER CENTRE C
      DO 4000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE C
        XYZ(1,3) = COORD(1,ICNTC)
        XYZ(2,3) = COORD(2,ICNTC)
        XYZ(3,3) = COORD(3,ICNTC)
C
C     LOOP OVER CENTRE D
      DO 4000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE D
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
        NFUNS(3) = NFUNCT(LQN(3)+1,ICNTC)
        DO KBAS=1,NFUNS(3)
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
        NFUNS(4) = NFUNCT(LQN(4)+1,ICNTD)
        DO LBAS=1,NFUNS(4)
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
      ICDLL = ILLAD(ICNTC,ICNTD,KC,KD,MC,MD)
      ICDSS = ISSAD(ICNTC,ICNTD,KC,KD,MC,MD)
C
C**********************************************************************C
C     FOR THIS CHOICE OF C AND D, COMPUTE ADDRESSES AND PHASES         C
C**********************************************************************C
C
C     CALCULATE BLOCK INDICES FOR {CD} COMBINATIONS
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
      IQ34 = (IQ3*(IQ3-1))/2 + IQ4
C
C     FURTHER DEFINE STARTING ADDRESSES FOR {CD} BASIS LABELS
      IC1 = LARGE(ICNTC,KC,MJC  ) + NADDCD
      ID1 = LARGE(ICNTD,KD,MJD  ) + NADDCD
C
      IC2 = LARGE(ICNTC,KC,MJC+1) + NADDCD
      ID2 = LARGE(ICNTD,KD,MJD+1) + NADDCD
C
      JC1 = LARGE(ICNTC,KC,MJC  ) + NADDCD
      JD1 = LARGE(ICNTD,KD,MJD  ) + NADDCD
C
      JC2 = LARGE(ICNTC,KC,MJC+1) + NADDCD
      JD2 = LARGE(ICNTD,KD,MJD+1) + NADDCD
C        
      IF((KQN(3)*KQN(4)).GT.0) THEN 
        PKCD = 1.0D0
      ELSE
        PKCD =-1.0D0
      ENDIF
C
C     CALCULATE MQN PHASE FACTORS FOR PERMUTING INTEGRALS
      PMCD1 = DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PMCD2 = DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
      G1 = PKCD*PMCD1
      G2 = PKCD*PMCD2
C
C**********************************************************************C
C     SKIP BATCHES WHICH CONFORM TO INTEGRAL SYMMETRIES                C
C**********************************************************************C
C
C     DIATOMIC MOLECULES CARRY STRICT SELECTION RULES ON MQNS
      IF(NCNT.LE.2) THEN
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) GOTO 3999
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) GOTO 3999
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) GOTO 3999
        GOTO 4000
      ENDIF
3999  CONTINUE
C
C     DECISION TREE FOR SKIPPING CONTRIBUTIONS DUE TO INTEGRAL SYMMETRY
      IF(IQ1.LT.IQ2)   GOTO 4000
      IF(IQ3.LT.IQ4)   GOTO 4000
      IF(IQ12.LT.IQ34) GOTO 4000
C
C     INDICATE BLOCKS TO BE INCLUDED AHEAD GIVEN A,B,C,D BASIS QNMS...
C     A =/= B AND C =/= D WITH AB LIST VALUE =/= CD LIST VALUE
      IF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 1
C     A=/=B AND C=/=D WITH IND(AB) = IND(CD)
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.EQ.IQ34) THEN
        ITSCF = 2
C     A=/=B AND C = D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 3
C     A = B AND C=/=D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 4
C     A = B AND C = D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 5
C     A = B AND C = D WITH IND(AB) = IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.EQ.IQ34) THEN
        ITSCF = 6
C     COMBINATION OF A,B,C,D NOT TO BE INCLUDED -- USE MATRIX CONJ LATER
      ELSE
        GO TO 4000
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
C
C**********************************************************************C
C     FOURTH LAYER OF LOOPS, OVER BASIS FUNCTIONS A AND B (5000)       C
C -------------------------------------------------------------------- C
C     THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH         C
C     GENERATE THE GDIR AND GXCH MATRICES FROM THE SPINOR INTEGRALS.   C
C     THESE INCLUDE IMPLICIT PHASE FACTORS FOR THE PERMUTATION OF      C
C     KQN(1) <-> KQN(2) AND MQN(1) <-> MQN(2)                          C
C -------------------------------------------------------------------- C
C     (RSCF 86, 87)                                                    C
C**********************************************************************C
C
C     LOOP OVER ELEMENTS OF FOCK MATRIX BLOCK
      DO 5000 IBAS=1,NFUNS(1)
      DO 5000 JBAS=1,NFUNS(2)
C
C       GENERATE BATCH OF ELECTRON REPULSION INTEGRALS
        CALL ERI(RR,XYZ,KQN,MQN,EXPT,NFUNS,ITQN,IBAS,JBAS,IEAB,IECD)
C
C       THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH
C       GENERATE THE GDIR/GXCH MATRIX FROM THE SPINOR INTEGRALS.
C
C       FIRST IFLG BATCH (DIRECT)
        IF(IFLG(1).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GDIR(IA1+IBAS,JB1+JBAS) = GDIR(IA1+IBAS,JB1+JBAS)
     &           +       RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
     &           +    G1*RR(M, 4)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M, 2)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IB)
     &           +    G2*RR(M, 3)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M, 1)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IB)
C
              GDIR(IA1+IBAS,JB2+JBAS) = GDIR(IA1+IBAS,JB2+JBAS)
     &           +       RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
     &           +    G1*RR(M, 8)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M, 6)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IB)
     &           +    G2*RR(M, 7)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M, 5)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IB)
C
              GDIR(IA2+IBAS,JB1+JBAS) = GDIR(IA2+IBAS,JB1+JBAS) 
     &           +       RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
     &           +    G1*RR(M,12)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M,10)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IB)
     &           +    G2*RR(M,11)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M, 9)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IB)
C
              GDIR(IA2+IBAS,JB2+JBAS) = GDIR(IA2+IBAS,JB2+JBAS) 
     &           +       RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
     &           +    G1*RR(M,16)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M,14)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IB)
     &           +    G2*RR(M,15)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M,13)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IB)
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
     &           +       RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
C
              GDIR(IA1+IBAS,JB2+JBAS) = GDIR(IA1+IBAS,JB2+JBAS) 
     &           +       RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
C
              GDIR(IA2+IBAS,JB1+JBAS) = GDIR(IA2+IBAS,JB1+JBAS) 
     &           +       RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
C
              GDIR(IA2+IBAS,JB2+JBAS) = GDIR(IA2+IBAS,JB2+JBAS) 
     &           +       RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
C
            ENDDO
          ENDDO
        ENDIF

C       THIRD IFLG BATCH (DIRECT)
        IF(IFLG(3).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GDIR(IC1+KBAS,JD1+LBAS) = GDIR(IC1+KBAS,JD1+LBAS) 
     &           +       RR(M, 1)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 5)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M, 9)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,13)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
     &           +    F1*RR(M,13)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F2*RR(M, 5)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IB)
     &           +    F2*RR(M, 9)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F1*RR(M, 1)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IB)
C
              GDIR(IC1+KBAS,JD2+LBAS) = GDIR(IC1+KBAS,JD2+LBAS)
     &           +       RR(M, 2)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 6)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,10)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,14)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
     &           +    F1*RR(M,14)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F2*RR(M, 6)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IB)
     &           +    F2*RR(M,10)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F1*RR(M, 2)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IB)
C
              GDIR(IC2+KBAS,JD1+LBAS) = GDIR(IC2+KBAS,JD1+LBAS) 
     &           +       RR(M, 3)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 7)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,11)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,15)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
     &           +    F1*RR(M,15)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F2*RR(M, 7)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IB)
     &           +    F2*RR(M,11)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F1*RR(M, 3)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IB)
C
              GDIR(IC2+KBAS,JD2+LBAS) = GDIR(IC2+KBAS,JD2+LBAS) 
     &           +       RR(M, 4)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 8)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,12)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,16)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
     &           +    F1*RR(M,16)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F2*RR(M, 8)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IB)
     &           +    F2*RR(M,12)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F1*RR(M, 4)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IB)
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
     &           +       RR(M, 1)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 5)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M, 9)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,13)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
C
              GDIR(IC1+KBAS,JD2+LBAS) = GDIR(IC1+KBAS,JD2+LBAS)
     &           +       RR(M, 2)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 6)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,10)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,14)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
C
              GDIR(IC2+KBAS,JD1+LBAS) = GDIR(IC2+KBAS,JD1+LBAS) 
     &           +       RR(M, 3)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 7)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,11)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,15)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
C
              GDIR(IC2+KBAS,JD2+LBAS) = GDIR(IC2+KBAS,JD2+LBAS) 
     &           +       RR(M, 4)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 8)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,12)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,16)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
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
     &           +    G1*RR(M, 4)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G1*RR(M, 8)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IB)
     &           +    G2*RR(M, 3)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G2*RR(M, 7)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IB)
C
              GXCH(IA1+IBAS,JC2+KBAS) = GXCH(IA1+IBAS,JC2+KBAS)
     &           +    G2*RR(M, 2)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G2*RR(M, 6)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IB)
     &           +    G1*RR(M, 1)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G1*RR(M, 5)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IB)
C
              GXCH(IA2+IBAS,JC1+KBAS) = GXCH(IA2+IBAS,JC1+KBAS)
     &           +    G1*RR(M,12)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G1*RR(M,16)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IB)
     &           +    G2*RR(M,11)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G2*RR(M,15)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IB)
C
              GXCH(IA2+IBAS,JC2+KBAS) = GXCH(IA2+IBAS,JC2+KBAS)
     &           +    G2*RR(M,10)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G2*RR(M,14)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IB)
     &           +    G1*RR(M, 9)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G1*RR(M,13)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IB)
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
     &           +    F1*RR(M,13)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F1*RR(M,14)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IB)
     &           +    F2*RR(M, 9)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F2*RR(M,10)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IB)
C
              GXCH(IC1+KBAS,JA2+IBAS) = GXCH(IC1+KBAS,JA2+IBAS)
     &           +    F2*RR(M, 5)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F2*RR(M, 6)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IB)
     &           +    F1*RR(M, 1)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F1*RR(M, 2)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IB)
C
              GXCH(IC2+KBAS,JA1+IBAS) = GXCH(IC2+KBAS,JA1+IBAS)
     &           +    F1*RR(M,15)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F1*RR(M,16)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IB)
     &           +    F2*RR(M,11)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F2*RR(M,12)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IB)
C
              GXCH(IC2+KBAS,JA2+IBAS) = GXCH(IC2+KBAS,JA2+IBAS)
     &           +    F2*RR(M, 7)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F2*RR(M, 8)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IB)
     &           +    F1*RR(M, 3)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F1*RR(M, 4)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IB)
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
     &           + F1*G1*RR(M,16)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F2*G1*RR(M, 8)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IB)
     &           + F1*G2*RR(M,15)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F2*G2*RR(M, 7)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IB)
C
              GXCH(IB1+JBAS,JC2+KBAS) = GXCH(IB1+JBAS,JC2+KBAS)
     &           + F1*G2*RR(M,14)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F2*G2*RR(M, 6)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IB)
     &           + F1*G1*RR(M,13)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F2*G1*RR(M, 5)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IB)
C
              GXCH(IB2+JBAS,JC1+KBAS) = GXCH(IB2+JBAS,JC1+KBAS)
     &           + F2*G1*RR(M,12)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F1*G1*RR(M, 4)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IB)
     &           + F2*G2*RR(M,11)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F1*G2*RR(M, 3)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IB)
C
              GXCH(IB2+JBAS,JC2+KBAS) = GXCH(IB2+JBAS,JC2+KBAS)
     &           + F2*G2*RR(M,10)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F1*G2*RR(M, 2)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IB)
     &           + F2*G1*RR(M, 9)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F1*G1*RR(M, 1)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IB)
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
     &           +       RR(M, 1)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M, 2)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IB)
     &           +       RR(M, 9)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M,10)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IB)
C
              GXCH(IC1+KBAS,JB2+JBAS) = GXCH(IC1+KBAS,JB2+JBAS)
     &           +       RR(M, 5)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M, 6)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IB)
     &           +       RR(M,13)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M,14)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IB)
C
              GXCH(IC2+KBAS,JB1+JBAS) = GXCH(IC2+KBAS,JB1+JBAS)
     &           +       RR(M, 3)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M, 4)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IB)
     &           +       RR(M,11)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M,12)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IB)
C
              GXCH(IC2+KBAS,JB2+JBAS) = GXCH(IC2+KBAS,JB2+JBAS)
     &           +       RR(M, 7)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M, 8)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IB)
     &           +       RR(M,15)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M,16)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IB)
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
     &           +       RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IB)
C
              GXCH(IA1+IBAS,JD2+LBAS) = GXCH(IA1+IBAS,JD2+LBAS)
     &           +       RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IB)
C
              GXCH(IA2+IBAS,JD1+LBAS) = GXCH(IA2+IBAS,JD1+LBAS)
     &           +       RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IB)
C
              GXCH(IA2+IBAS,JD2+LBAS) = GXCH(IA2+IBAS,JD2+LBAS)
     &           +       RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IB)
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
     &           +    F1*RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F2*RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IB)
     &           +    F1*RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F2*RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IB)
C
              GXCH(IB1+JBAS,JD2+LBAS) = GXCH(IB1+JBAS,JD2+LBAS)
     &           +    F1*RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F2*RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IB)
     &           +    F1*RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F2*RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IB)
C
              GXCH(IB2+JBAS,JD1+LBAS) = GXCH(IB2+JBAS,JD1+LBAS)
     &           +    F2*RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F1*RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IB)
     &           +    F2*RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F1*RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IB)
C
              GXCH(IB2+JBAS,JD2+LBAS) = GXCH(IB2+JBAS,JD2+LBAS)
     &           +    F2*RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F1*RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IB)
     &           +    F2*RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F1*RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IB)
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
     &           +    G1*RR(M, 4)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M, 2)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IB)
     &           +    G1*RR(M,12)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M,10)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IB)
C
              GXCH(ID1+LBAS,JB2+JBAS) = GXCH(ID1+LBAS,JB2+JBAS)
     &           +    G1*RR(M, 8)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M, 6)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IB)
     &           +    G1*RR(M,16)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M,14)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IB)
C
              GXCH(ID2+LBAS,JB1+JBAS) = GXCH(ID2+LBAS,JB1+JBAS)
     &           +    G2*RR(M, 3)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M, 1)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IB)
     &           +    G2*RR(M,11)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M, 9)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IB)
C
              GXCH(ID2+LBAS,JB2+JBAS) = GXCH(ID2+LBAS,JB2+JBAS)
     &           +    G2*RR(M, 7)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M, 5)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IB)
     &           +    G2*RR(M,15)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M,13)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IB)
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
          IF(ICNLAB(I).NE.ICNLAB(J)) GOTO 400
          IF(KQNLAB(I).NE.KQNLAB(J)) GOTO 400
          IF(MQNLAB(I).NE.MQNLAB(J)) GOTO 400
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
C
C**********************************************************************C
C     CALCULATE INTERACTION ENERGY FROM ORBITAL PAIR (A,B)             C
C**********************************************************************C

C
      ETMP3 = DCMPLX(0.0D0,0.0D0)
      ETMP4 = DCMPLX(0.0D0,0.0D0)
      DO J=1,NDIM
        DO I=1,NDIM
          ETMP3 = ETMP3 + DCONJG(C(I,IA))*C(J,IA)*GDIR(I,J)
          ETMP4 = ETMP4 + DCONJG(C(I,IA))*C(J,IA)*GXCH(I,J)
        ENDDO
      ENDDO
      EAB1(IOCCA,IOCCB,4) = DREAL(ETMP3)
      EAB1(IOCCA,IOCCB,5) = DREAL(ETMP4)
      EAB1(IOCCA,IOCCB,6) = DREAL(ETMP1 + ETMP2 + ETMP3 - ETMP4)
C
C     OUTPUT ENERGIES TO TERMINAL
      WRITE(6,11) IOCCA,IOCCB,(EAB1(IOCCA,IOCCB,N),N=3,6)
      WRITE(7,11) IOCCA,IOCCB,(EAB1(IOCCA,IOCCB,N),N=3,6)
C
      IF(IOCCA.NE.IOCCB) THEN
        DO N=1,6
          EAB1(IOCCB,IOCCA,N) = EAB1(IOCCA,IOCCB,N)
        ENDDO
      ENDIF
C
C     END LOOP OVER UNIQUE OCCUPIED ORBITAL COMBINATIONS (A,B)
1000  CONTINUE
C
      CALL CPU_TIME(TFIN)
C
C     WRITE RESULTS OF EAB ENERGIES TO AN EXTERNAL FILE
      OPEN(UNIT=10,FILE=STRING(:LN)//'_MBPT1.dat',STATUS='UNKNOWN')
      REWIND(UNIT=10)
      DO IOCCA=1,NOCC
        DO IOCCB=1,NOCC
          WRITE(10,*) (EAB1(IOCCA,IOCCB,N),N=1,6)
        ENDDO
      ENDDO
      CLOSE(UNIT=10)
C
      WRITE(6, *) REPEAT('-',62)
      WRITE(7, *) REPEAT('-',62)
C
C     CALCULATE MBPT1 SINGLE-PARTICLE ENERGIES AND MOLECULAR TOTALS
      EENC1 = 0.0D0
      EKIN1 = 0.0D0
      EBAR1 = 0.0D0
      EDIR1 = 0.0D0
      EXCH1 = 0.0D0
      DO IOCCA=1,NOCC
        DO N=1,5
          EA1(IOCCA,N) = 0.0D0
          DO IOCCB=1,NOCC
            EA1(IOCCA,N) = EA1(IOCCA,N) + EAB1(IOCCA,IOCCB,N)
          ENDDO      
        ENDDO
        EA1(IOCCA,6) = EA1(IOCCA,3) + EA1(IOCCA,4) - EA1(IOCCA,5)
        EENC1 = EENC1 + EA1(IOCCA,1)
        EKIN1 = EKIN1 + EA1(IOCCA,2)
        EBAR1 = EBAR1 + EA1(IOCCA,3)
        EDIR1 = EDIR1 + EA1(IOCCA,4)*0.5D0
        EXCH1 = EXCH1 + EA1(IOCCA,5)*0.5D0
      ENDDO
      ETOT1 = EBAR1 + EDIR1 - EXCH1 + ENUC
C
C     ORBITAL SUMMARIES
      WRITE(6, *) 
      WRITE(7, *) 
      WRITE(6, *) REPEAT(' ',16),'MBPT1 single particle summary'
      WRITE(7, *) REPEAT(' ',16),'MBPT1 single particle summary'
      WRITE(6, *) REPEAT('-',62)
      WRITE(7, *) REPEAT('-',62)
      WRITE(6,10) '  a    ','E1(H)','E1(J)','E1(K)',' E1(a)'
      WRITE(7,10) '  a    ','E1(H)','E1(J)','E1(K)',' E1(a)'
      WRITE(6, *) REPEAT('-',62)
      WRITE(7, *) REPEAT('-',62)  
      DO IOCCA=1,NOCC
        WRITE(6,12) IOCCA,(EA1(IOCCA,N),N=3,6)
        WRITE(7,12) IOCCA,(EA1(IOCCA,N),N=3,6)
      ENDDO
      WRITE(6, *) REPEAT('-',62)
      WRITE(7, *) REPEAT('-',62)
C
C     TOTAL ENERGIES
      WRITE(6, *) 
      WRITE(7, *) 
      WRITE(6, *) REPEAT(' ',20),'MBPT1 molecular summary'
      WRITE(7, *) REPEAT(' ',20),'MBPT1 molecular summary'
      WRITE(6, *) REPEAT('-',62)
      WRITE(7, *) REPEAT('-',62)
      WRITE(6,84) 'Hartree-Fock electron-nucleus','E1(V)',EENC1
      WRITE(7,84) 'Hartree-Fock electron-nucleus','E1(V)',EENC1
      WRITE(6,84) 'Hartree-Fock electron kinetic','E1(T)',EKIN1
      WRITE(7,84) 'Hartree-Fock electron kinetic','E1(T)',EKIN1
      WRITE(6,84) 'Hartree-Fock Coulomb direct  ','E1(J)',EDIR1
      WRITE(7,84) 'Hartree-Fock Coulomb direct  ','E1(J)',EDIR1
      WRITE(6,84) 'Hartree-Fock Coulomb exchange','E1(K)',EXCH1
      WRITE(7,84) 'Hartree-Fock Coulomb exchange','E1(K)',EXCH1
      WRITE(6, *) REPEAT('-',62)
      WRITE(7, *) REPEAT('-',62)
      WRITE(6,84) 'Nuclear repulsion            ','E0(N)',ENUC
      WRITE(7,84) 'Nuclear repulsion            ','E0(N)',ENUC
      WRITE(6,84) 'Hartree-Fock one-electron    ','E1(H)',EBAR1
      WRITE(7,84) 'Hartree-Fock one-electron    ','E1(H)',EBAR1
      WRITE(6,84) 'Hartree-Fock Coulomb total   ','E1(G)',EDIR1-EXCH1
      WRITE(7,84) 'Hartree-Fock Coulomb total   ','E1(G)',EDIR1-EXCH1
      WRITE(6,84) 'Hartree-Fock molecular energy','E1   ',ETOT1
      WRITE(7,84) 'Hartree-Fock molecular energy','E1   ',ETOT1
      WRITE(6, *) REPEAT('-',62)
      WRITE(7, *) REPEAT('-',62)
      WRITE(6,87) 'MBPT1 time                   ',TIMEHMS(TFIN-TBGN)
      WRITE(7,87) 'MBPT1 time                   ',TIMEHMS(TFIN-TBGN)
      WRITE(6, *) REPEAT('-',62)
      WRITE(7, *) REPEAT('-',62)
C
      RETURN
      END

