      SUBROUTINE COULOMB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C    CCCCCC   OOOOOO  UU    UU LL      OOOOOO  MM       MM BBBBBBB    C
C   CC    CC OO    OO UU    UU LL     OO    OO MMM     MMM BB    BB   C
C   CC       OO    OO UU    UU LL     OO    OO MMMM   MMMM BB    BB   C
C   CC       OO    OO UU    UU LL     OO    OO MM MM MM MM BBBBBBB    C
C   CC       OO    OO UU    UU LL     OO    OO MM  MMM  MM BB    BB   C
C   CC    CC OO    OO UU    UU LL     OO    OO MM   M   MM BB    BB   C
C    CCCCCC   OOOOOO   UUUUUU  LLLLLLL OOOOOO  MM       MM BBBBBBB    C
C                                                                     C
C ------------------------------------------------------------------- C
C    COULOMB GENERATES ELECTRON REPULSION INTEGRALS IN BATCHES AND    C
C    CALCULATES THE SCF COULOMB MATRIX (G), APPLYING IT DIRECTLY TO   C
C    THE FOCK MATRIX. IN THE CASE OF OPEN SUBSHELLS, THE Q MATRIX IS  C
C    ALSO INCLUDED. INTEGRAL SYMMETRY PARTIALLY EXPLOITED, BUT WITH   C
C    ROOM FOR IMPROVEMENT (GEOMETRIC SYMM, R-INT SYMM, E-COEFF SYMM). C
C ------------------------------------------------------------------- C
C    NOTE: THIS ROUTINE COULD BENEFIT FROM PARALLELISATION -- OPENMP. C
C*********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=MKP,MBS=26,MB2=MBS*MBS,
     &          IL4=2*(MKP-1),MLL=(MKP+8)*(MKP+9)*(MKP+10)/6,MRC=969)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           GDIR(MDM,MDM),GXCH(MDM,MDM),QDIR(MDM,MDM),
     &           QXCH(MDM,MDM),BDIR(MDM,MDM),BXCH(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 RR(MB2,16)
C
      DIMENSION EXPT(MBS,4),KQN(4),MQN(4),NFUNS(4),LQN(4)
      DIMENSION ITQN(2),IFLG(11),ISCF(11,6)
      DIMENSION INDEX(MCT,-MKP:MKP,2*(MKP+1)*MKP)
      DIMENSION T1(MDM),T2(MDM),T3(MDM)
C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/LBL2/LDIAG(100),NDIG
      COMMON/MTRX/OVAP,HNUC,HKIN,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITER,IALL,IRUN
      COMMON/QNMS/ICNLAB(MDM),KQNLAB(MDM),MQNLAB(MDM)
      COMMON/SHLL/ACFF,BCFF,FOPEN,ICLOSE(100),IOPEN(6),NCLOSE,NOPEN
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR,IOCCM0
      COMMON/XYZ4/XYZ(3,4)
C
C     ISCF TELLS WHICH INTEGRALS TO INCLUDE BASED ON OVERLAP COMBINATION
      DATA ISCF/1,1,1,1,0,1,0,1,0,0,0,
     &          1,0,1,1,0,1,0,0,0,0,0,
     &          0,1,1,0,0,1,0,0,0,1,0,
     &          1,0,1,1,0,0,0,0,0,0,1,
     &          0,0,1,0,0,0,0,0,0,1,1,
     &          0,0,1,0,0,0,0,0,0,1,0/  
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
C*********************************************************************C
C     INDEXING ROUTINE: SO WE CAN SET BASIS FUNCTION LABELS BY BLOCK  C
C*********************************************************************C
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
C*********************************************************************C
C     SET THE LIMITS FOR RUNNING OVER DENSITY COMBINATIONS            C
C*********************************************************************C
C
      IF(HMLTN.EQ.'NORL') THEN
        ITSTRT = 1
        ITSTOP = 1
        ITSKIP = 1
      ELSEIF(HMLTN.EQ.'DHFR'.OR.HMLTN.EQ.'DHFB') THEN
        ITSTRT = 1
        ITSTOP = 4
        ITSKIP = 3
      ENDIF
C
C*********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 1000)     C
C*********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR A AND B: T = {L} OR {L,S}
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
C     LOOP OVER CENTRE A
      DO 1000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 1000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTRE/MULTI-CENTRE OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 2
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
        NFUNA    = NFUNCT(LQN(1)+1,ICNTA)
        NFUNS(1) = NFUNA
C
        DO IBAS=1,NFUNA
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
        NFUNB    = NFUNCT(LQN(2)+1,ICNTB)
        NFUNS(2) = NFUNB
C
        DO JBAS=1,NFUNB
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
C
C*********************************************************************C
C      SECOND LAYER OF LOOPS, OVER CENTRES C AND D (INDEX 2000/2500)  C
C*********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR C AND D: T' = {L} OR {L,S}
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

C     LOOP OVER CENTRE C
      DO 2000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE C
        XYZ(1,3) = COORD(1,ICNTC)
        XYZ(2,3) = COORD(2,ICNTC)
        XYZ(3,3) = COORD(3,ICNTC)
C
C     LOOP OVER CENTRE D
      DO 2000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE D
        XYZ(1,4) = COORD(1,ICNTD)
        XYZ(2,4) = COORD(2,ICNTD)
        XYZ(3,4) = COORD(3,ICNTD)
C
C       PARAMETER FOR SINGLE-CENTRE/MULTI-CENTRE OVERLAP OVER C AND D
        IF(ICNTC.EQ.ICNTD) THEN
          INUCCD = 1
        ELSE
          INUCCD = 2
        ENDIF
C
C       PARAMETER FOR ATOMIC OR MULTICNTRE INTEGRAL
        IF(INUCAB*INUCCD.EQ.1.AND.ICNTA.EQ.ICNTC) THEN
          IATOM = 1
        ELSE
          IATOM = 0
        ENDIF
C
C ***   STAGES: DECISION TREE FOR SKIPPING MULTI-CENTRE CONTRIBUTIONS
        IF(IATOM.EQ.0) THEN
C
C >>      STAGE 0: INCLUDE ONLY (LL|LL) REPULSION INTEGRALS
          IF(IALL.EQ.0.AND.IT1+IT2.GT.2) THEN
            GOTO 2200
          ENDIF
C
C >>      STAGE 1: INCLUDE ONLY (LL|SS) AND (SS|LL) REPULSION INTEGRALS
          IF(IALL.EQ.1.AND.IT1+IT2.GT.5) THEN
            GOTO 2200
          ENDIF
C
C ***   END OF IALL DECISION TREE
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
        NFUNC    = NFUNCT(LQN(3)+1,ICNTC)
        NFUNS(3) = NFUNC
C
        DO KBAS=1,NFUNC
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
        NFUND    = NFUNCT(LQN(4)+1,ICNTD)
        NFUNS(4) = NFUND
C
        DO LBAS=1,NFUND
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
C
C*********************************************************************C
C     FOR THIS CHOICE OF A,B,C AND D, COMPUTE ADDRESSES AND PHASES    C
C*********************************************************************C
C
C     CALCULATE BLOCK INDICES FOR {ABCD} COMBINATIONS
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
      IQ12 = (IQ1*(IQ1-1))/2 + IQ2
      IQ34 = (IQ3*(IQ3-1))/2 + IQ4
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
C*********************************************************************C
C     SKIP BATCHES WHICH CONFORM TO INTEGRAL SYMMETRIES               C
C*********************************************************************C
C
      IF(NCNT.LE.2) THEN
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) GOTO 2998
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) GOTO 2998
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) GOTO 2998
        GOTO 2999
      ENDIF
2998  CONTINUE

C     DECISION TREE FOR SKIPPING CONTRIBUTIONS DUE TO INTEGRAL SYMMETRY
      IF(IQ1.LT.IQ2)   GOTO 2999
      IF(IQ3.LT.IQ4)   GOTO 2999
      IF(IQ12.LT.IQ34) GOTO 2999
C
C     LABEL VARIOUS TYPES OF OVERLAP COMBINATIONS USING ITSCF OR SKIP
      IF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 1
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.EQ.IQ34) THEN
        ITSCF = 2
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 3
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 4
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 5
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.EQ.IQ34) THEN
        ITSCF = 6
      ELSE
        GO TO 2999
      ENDIF
C
C     READ IN FLAG VALUES FROM ISCF DATA BLOCK
      DO M=1,11
        IFLG(M) = ISCF(M,ITSCF)
      ENDDO
C
C     INCLUDE SPECIAL CASES
      IF(ITSCF.EQ.1.AND.IQ1.EQ.IQ3) IFLG(5)=1
      IF(ITSCF.EQ.1.AND.IQ2.EQ.IQ4) IFLG(7)=1
      IF(ITSCF.EQ.1.AND.IQ2.EQ.IQ3) IFLG(9)=1
      IF(ITSCF.EQ.3.AND.IQ2.EQ.IQ3) IFLG(9)=1
      IF(ITSCF.EQ.4.AND.IQ2.EQ.IQ3) IFLG(9)=1
C
C*********************************************************************C
C     THIRD LAYER OF LOOPS, OVER BASIS FUNCTIONS A AND B (3000)       C
C ------------------------------------------------------------------- C
C     THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH        C
C     GENERATE CLOSED AND OPEN COULOMB MATRICES FROM SPINOR INTEGRALS.C
C     THESE INCLUDE IMPLICIT PHASE FACTORS FOR THE PERMUTATION OF     C
C     KQN(1) <-> KQN(2) AND MQN(1) <-> MQN(2)                         C
C ------------------------------------------------------------------- C
C     (RSCF 86, 87)                                                   C
C*********************************************************************C
C
C     LOOP OVER ELEMENTS OF FOCK MATRIX BLOCK
      DO 3000 I=1,NFUNA
      DO 3000 J=1,NFUNB
C
C       GENERATE BATCH OF ELECTRON REPULSION INTEGRALS
        CALL ERI(RR,XYZ,KQN,MQN,EXPT,NFUNS,ITQN,I,J,IEAB,IECD)
C
C       THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH
C       GENERATE THE GDIR AND GXCH MATRIX FROM THE SPINOR INTEGRALS.
C
C*********************************************************************C
C     CONSTRUCT CLOSED-SHELL COULOMB INTERACTION MATRICES, QDIR/QXCH. C
C*********************************************************************C
C
C       ATOMIC/DIATOMIC MATRIX ELEMENTS...
C
C       FIRST IFLG BATCH
        IF(IFLG(1).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GDIR(IA1+I,JB1+J) = GDIR(IA1+I,JB1+J)
     &                          +            RR(M, 1)*DENC(IC1+K,JD1+L)
     &                          +            RR(M, 4)*DENC(IC2+K,JD2+L)
     &                          +         G1*RR(M, 4)*DENC(JD1+L,IC1+K)
     &                          +         G1*RR(M, 1)*DENC(JD2+L,IC2+K)
C
              GDIR(IA2+I,JB2+J) = GDIR(IA2+I,JB2+J)
     &                          +            RR(M,13)*DENC(IC1+K,JD1+L)
     &                          +            RR(M,16)*DENC(IC2+K,JD2+L)
     &                          +         G1*RR(M,16)*DENC(JD1+L,IC1+K)
     &                          +         G1*RR(M,13)*DENC(JD2+L,IC2+K)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SECOND IFLG BATCH
        IF(IFLG(2).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GDIR(IC1+K,JD1+L) = GDIR(IC1+K,JD1+L)
     &                          +            RR(M, 1)*DENC(IA1+I,JB1+J)
     &                          +            RR(M,13)*DENC(IA2+I,JB2+J)
     &                          +         F1*RR(M,13)*DENC(JB1+J,IA1+I)
     &                          +         F1*RR(M, 1)*DENC(JB2+J,IA2+I)
C
              GDIR(IC2+K,JD2+L) = GDIR(IC2+K,JD2+L)
     &                          +            RR(M, 4)*DENC(IA1+I,JB1+J)
     &                          +            RR(M,16)*DENC(IA2+I,JB2+J)
     &                          +         F1*RR(M,16)*DENC(JB1+J,IA1+I)
     &                          +         F1*RR(M, 4)*DENC(JB2+J,IA2+I)
C
            ENDDO
          ENDDO
        ENDIF
C
C       THIRD IFLG BATCH          
        IF(IFLG(3).EQ.1) THEN
          DO L=1,NFUND
            DO K=1,NFUNC
              M = (K-1)*NFUND + L
C  
              GXCH(IA1+I,JD1+L) = GXCH(IA1+I,JD1+L)
     &                          +            RR(M, 1)*DENC(IC1+K,JB1+J) 
     &                          +            RR(M, 7)*DENC(IC2+K,JB2+J)
C
              GXCH(IA2+I,JD2+L) = GXCH(IA2+I,JD2+L)
     &                          +            RR(M,10)*DENC(IC1+K,JB1+J) 
     &                          +            RR(M,16)*DENC(IC2+K,JB2+J)
C
            ENDDO
          ENDDO
        ENDIF
C
C       FOURTH IFLG BATCH          
        IF(IFLG(4).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GXCH(IA1+I,JC1+K) = GXCH(IA1+I,JC1+K)
     &                          +         G1*RR(M, 4)*DENC(ID1+L,JB1+J) 
     &                          +         G2*RR(M, 7)*DENC(ID2+L,JB2+J)
C
              GXCH(IA2+I,JC2+K) = GXCH(IA2+I,JC2+K)
     &                          +         G2*RR(M,10)*DENC(ID1+L,JB1+J) 
     &                          +         G1*RR(M,13)*DENC(ID2+L,JB2+J)
C
            ENDDO
          ENDDO
        ENDIF
C
C       FIFTH IFLG BATCH
        IF(IFLG(5).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GXCH(IC1+K,JA1+I) = GXCH(IC1+K,JA1+I)
     &                          +         F1*RR(M,13)*DENC(JB1+J,ID1+L)
     &                          +         F2*RR(M,10)*DENC(JB2+J,ID2+L)
C     
              GXCH(IC2+K,JA2+I) = GXCH(IC2+K,JA2+I)
     &                          +         F2*RR(M, 7)*DENC(JB1+J,ID1+L)
     &                          +         F1*RR(M, 4)*DENC(JB2+J,ID2+L)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SIXTH IFLG BATCH
        IF(IFLG(6).EQ.1) THEN
          DO L=1,NFUND
            DO K=1,NFUNC
              M = (K-1)*NFUND + L
C
              GXCH(IB1+J,JD1+L) = GXCH(IB1+J,JD1+L)
     &                          +         F1*RR(M,13)*DENC(IC1+K,JA1+I) 
     &                          +         F2*RR(M, 7)*DENC(IC2+K,JA2+I)
C
              GXCH(IB2+J,JD2+L) = GXCH(IB2+J,JD2+L)
     &                          +         F2*RR(M,10)*DENC(IC1+K,JA1+I) 
     &                          +         F1*RR(M, 4)*DENC(IC2+K,JA2+I)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SEVENTH IFLG BATCH
        IF(IFLG(7).EQ.1) THEN
          DO L=1,NFUND
            DO K=1,NFUNC
              M = (K-1)*NFUND + L
C
              GXCH(ID1+L,JB1+J) = GXCH(ID1+L,JB1+J)
     &                          +         G1*RR(M, 4)*DENC(JA1+I,IC1+K)
     &                          +         G2*RR(M,10)*DENC(JA2+I,IC2+K)
C
              GXCH(ID2+L,JB2+J) = GXCH(ID2+L,JB2+J)
     &                          +         G2*RR(M, 7)*DENC(JA1+I,IC1+K)
     &                          +         G1*RR(M,13)*DENC(JA2+I,IC2+K)
C
            ENDDO
          ENDDO
        ENDIF
C
C       EIGHTH IFLG BATCH
        IF(IFLG(8).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GXCH(IB1+J,JC1+K) = GXCH(IB1+J,JC1+K)
     &                          +      F1*G1*RR(M,16)*DENC(ID1+L,JA1+I) 
     &                          +      F2*G2*RR(M, 7)*DENC(ID2+L,JA2+I)
C
              GXCH(IB2+J,JC2+K) = GXCH(IB2+J,JC2+K)
     &                          +      F2*G2*RR(M,10)*DENC(ID1+L,JA1+I) 
     &                          +      F1*G1*RR(M, 1)*DENC(ID2+L,JA2+I)
C
            ENDDO
          ENDDO
        ENDIF
C
C       NINTH IFLG BATCH
        IF(IFLG(9).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GXCH(IC1+K,JB1+J) = GXCH(IC1+K,JB1+J)
     &                          +            RR(M, 1)*DENC(JA1+I,ID1+L)
     &                          +            RR(M,10)*DENC(JA2+I,ID2+L)
C
              GXCH(IC2+K,JB2+J) = GXCH(IC2+K,JB2+J)
     &                          +            RR(M, 7)*DENC(JA1+I,ID1+L)
     &                          +            RR(M,16)*DENC(JA2+I,ID2+L)
C
            ENDDO
          ENDDO           
        ENDIF
C
C       TENTH IFLG BATCH
        IF(IFLG(10).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GDIR(IA1+I,JB1+J) = GDIR(IA1+I,JB1+J)
     &                          +            RR(M, 1)*DENC(IC1+K,JD1+L) 
     &                          +            RR(M, 4)*DENC(IC2+K,JD2+L)
C
              GDIR(IA2+I,JB2+J) = GDIR(IA2+I,JB2+J)
     &                          +            RR(M,13)*DENC(IC1+K,JD1+L) 
     &                          +            RR(M,16)*DENC(IC2+K,JD2+L)
C
            ENDDO
          ENDDO
        ENDIF
C
C       ELEVENTH IFLG BATCH
        IF(IFLG(11).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GDIR(IC1+K,JD1+L) = GDIR(IC1+K,JD1+L)
     &                          +            RR(M, 1)*DENC(IA1+I,JB1+J) 
     &                          +            RR(M,13)*DENC(IA2+I,JB2+J)
C
              GDIR(IC2+K,JD2+L) = GDIR(IC2+K,JD2+L)
     &                          +            RR(M, 4)*DENC(IA1+I,JB1+J) 
     &                          +            RR(M,16)*DENC(IA2+I,JB2+J)
C
            ENDDO
          ENDDO
        ENDIF
C
C       ADDITIONAL MATRIX ELEMENTS (ASSUMES NO GEOMETRIC SYMMETRY)...
        IF(NCNT.LE.2) GOTO 390
C
C       FIRST IFLG BATCH
        IF(IFLG(1).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GDIR(IA1+I,JB1+J) = GDIR(IA1+I,JB1+J)
     &                          +            RR(M, 2)*DENC(IC1+K,JD2+L)
     &                          +            RR(M, 3)*DENC(IC2+K,JD1+L) 
     &                          +         G2*RR(M, 2)*DENC(JD1+L,IC2+K)
     &                          +         G2*RR(M, 3)*DENC(JD2+L,IC1+K) 
C
              GDIR(IA1+I,JB2+J) = GDIR(IA1+I,JB2+J)
     &                          +            RR(M, 5)*DENC(IC1+K,JD1+L) 
     &                          +            RR(M, 6)*DENC(IC1+K,JD2+L)
     &                          +            RR(M, 7)*DENC(IC2+K,JD1+L) 
     &                          +            RR(M, 8)*DENC(IC2+K,JD2+L)
     &                          +         G1*RR(M, 8)*DENC(JD1+L,IC1+K) 
     &                          +         G2*RR(M, 6)*DENC(JD1+L,IC2+K)
     &                          +         G2*RR(M, 7)*DENC(JD2+L,IC1+K) 
     &                          +         G1*RR(M, 5)*DENC(JD2+L,IC2+K)

C
              GDIR(IA2+I,JB1+J) = GDIR(IA2+I,JB1+J)
     &                          +            RR(M, 9)*DENC(IC1+K,JD1+L) 
     &                          +            RR(M,10)*DENC(IC1+K,JD2+L)
     &                          +            RR(M,11)*DENC(IC2+K,JD1+L) 
     &                          +            RR(M,12)*DENC(IC2+K,JD2+L)
     &                          +         G1*RR(M,12)*DENC(JD1+L,IC1+K) 
     &                          +         G2*RR(M,10)*DENC(JD1+L,IC2+K)
     &                          +         G2*RR(M,11)*DENC(JD2+L,IC1+K) 
     &                          +         G1*RR(M, 9)*DENC(JD2+L,IC2+K)
C
              GDIR(IA2+I,JB2+J) = GDIR(IA2+I,JB2+J)
     &                          +            RR(M,14)*DENC(IC1+K,JD2+L)
     &                          +            RR(M,15)*DENC(IC2+K,JD1+L) 
     &                          +         G2*RR(M,14)*DENC(JD1+L,IC2+K)
     &                          +         G2*RR(M,15)*DENC(JD2+L,IC1+K) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       SECOND IFLG BATCH
        IF(IFLG(2).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GDIR(IC1+K,JD1+L) = GDIR(IC1+K,JD1+L)
     &                          +            RR(M, 5)*DENC(IA1+I,JB2+J)
     &                          +            RR(M, 9)*DENC(IA2+I,JB1+J) 
     &                          +         F2*RR(M, 5)*DENC(JB1+J,IA2+I)
     &                          +         F2*RR(M, 9)*DENC(JB2+J,IA1+I) 
C
              GDIR(IC1+K,JD2+L) = GDIR(IC1+K,JD2+L)
     &                          +            RR(M, 2)*DENC(IA1+I,JB1+J) 
     &                          +            RR(M, 6)*DENC(IA1+I,JB2+J)
     &                          +            RR(M,10)*DENC(IA2+I,JB1+J) 
     &                          +            RR(M,14)*DENC(IA2+I,JB2+J)
     &                          +         F1*RR(M,14)*DENC(JB1+J,IA1+I) 
     &                          +         F2*RR(M, 6)*DENC(JB1+J,IA2+I)
     &                          +         F2*RR(M,10)*DENC(JB2+J,IA1+I) 
     &                          +         F1*RR(M, 2)*DENC(JB2+J,IA2+I)
C
              GDIR(IC2+K,JD1+L) = GDIR(IC2+K,JD1+L)
     &                          +            RR(M, 3)*DENC(IA1+I,JB1+J) 
     &                          +            RR(M, 7)*DENC(IA1+I,JB2+J)
     &                          +            RR(M,11)*DENC(IA2+I,JB1+J) 
     &                          +            RR(M,15)*DENC(IA2+I,JB2+J)
     &                          +         F1*RR(M,15)*DENC(JB1+J,IA1+I) 
     &                          +         F2*RR(M, 7)*DENC(JB1+J,IA2+I)
     &                          +         F2*RR(M,11)*DENC(JB2+J,IA1+I) 
     &                          +         F1*RR(M, 3)*DENC(JB2+J,IA2+I)
C
              GDIR(IC2+K,JD2+L) = GDIR(IC2+K,JD2+L)
     &                          +            RR(M, 8)*DENC(IA1+I,JB2+J)
     &                          +            RR(M,12)*DENC(IA2+I,JB1+J) 
     &                          +         F2*RR(M, 8)*DENC(JB1+J,IA2+I)
     &                          +         F2*RR(M,12)*DENC(JB2+J,IA1+I) 
C 
            ENDDO
          ENDDO
        ENDIF
C
C       THIRD IFLG BATCH          
        IF(IFLG(3).EQ.1) THEN
          DO L=1,NFUND
            DO K=1,NFUNC
              M = (K-1)*NFUND + L
C  
              GXCH(IA1+I,JD1+L) = GXCH(IA1+I,JD1+L)
     &                          +            RR(M, 5)*DENC(IC1+K,JB2+J)
     &                          +            RR(M, 3)*DENC(IC2+K,JB1+J) 
C
              GXCH(IA1+I,JD2+L) = GXCH(IA1+I,JD2+L)
     &                          +            RR(M, 2)*DENC(IC1+K,JB1+J) 
     &                          +            RR(M, 6)*DENC(IC1+K,JB2+J)
     &                          +            RR(M, 4)*DENC(IC2+K,JB1+J) 
     &                          +            RR(M, 8)*DENC(IC2+K,JB2+J)
C
              GXCH(IA2+I,JD1+L) = GXCH(IA2+I,JD1+L)
     &                          +            RR(M, 9)*DENC(IC1+K,JB1+J) 
     &                          +            RR(M,13)*DENC(IC1+K,JB2+J)
     &                          +            RR(M,11)*DENC(IC2+K,JB1+J) 
     &                          +            RR(M,15)*DENC(IC2+K,JB2+J)
C
              GXCH(IA2+I,JD2+L) = GXCH(IA2+I,JD2+L)
     &                          +            RR(M,14)*DENC(IC1+K,JB2+J)
     &                          +            RR(M,12)*DENC(IC2+K,JB1+J) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       FOURTH IFLG BATCH          
        IF(IFLG(4).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GXCH(IA1+I,JC1+K) = GXCH(IA1+I,JC1+K)
     &                          +         G1*RR(M, 8)*DENC(ID1+L,JB2+J)
     &                          +         G2*RR(M, 3)*DENC(ID2+L,JB1+J)
C
              GXCH(IA1+I,JC2+K) = GXCH(IA1+I,JC2+K)
     &                          +         G2*RR(M, 2)*DENC(ID1+L,JB1+J) 
     &                          +         G2*RR(M, 6)*DENC(ID1+L,JB2+J)
     &                          +         G1*RR(M, 1)*DENC(ID2+L,JB1+J) 
     &                          +         G1*RR(M, 5)*DENC(ID2+L,JB2+J)
C
              GXCH(IA2+I,JC1+K) = GXCH(IA2+I,JC1+K)
     &                          +         G1*RR(M,12)*DENC(ID1+L,JB1+J) 
     &                          +         G1*RR(M,16)*DENC(ID1+L,JB2+J)
     &                          +         G2*RR(M,11)*DENC(ID2+L,JB1+J) 
     &                          +         G2*RR(M,15)*DENC(ID2+L,JB2+J)
C
              GXCH(IA2+I,JC2+K) = GXCH(IA2+I,JC2+K)
     &                          +         G2*RR(M,14)*DENC(ID1+L,JB2+J)
     &                          +         G1*RR(M, 9)*DENC(ID2+L,JB1+J) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       FIFTH IFLG BATCH
        IF(IFLG(5).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GXCH(IC1+K,JA1+I) = GXCH(IC1+K,JA1+I)
     &                          +         F1*RR(M,14)*DENC(JB1+J,ID2+L)
     &                          +         F2*RR(M, 9)*DENC(JB2+J,ID1+L) 
C
              GXCH(IC1+K,JA2+I) = GXCH(IC1+K,JA2+I)
     &                          +         F2*RR(M, 5)*DENC(JB1+J,ID1+L) 
     &                          +         F2*RR(M, 6)*DENC(JB1+J,ID2+L)
     &                          +         F1*RR(M, 1)*DENC(JB2+J,ID1+L) 
     &                          +         F1*RR(M, 2)*DENC(JB2+J,ID2+L)
C
              GXCH(IC2+K,JA1+I) = GXCH(IC2+K,JA1+I)
     &                          +         F1*RR(M,15)*DENC(JB1+J,ID1+L) 
     &                          +         F1*RR(M,16)*DENC(JB1+J,ID2+L)
     &                          +         F2*RR(M,11)*DENC(JB2+J,ID1+L) 
     &                          +         F2*RR(M,12)*DENC(JB2+J,ID2+L)
C
              GXCH(IC2+K,JA2+I) = GXCH(IC2+K,JA2+I)
     &                          +         F2*RR(M, 8)*DENC(JB1+J,ID2+L)
     &                          +         F1*RR(M, 3)*DENC(JB2+J,ID1+L) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       SIXTH IFLG BATCH
        IF(IFLG(6).EQ.1) THEN
          DO L=1,NFUND
            DO K=1,NFUNC
              M = (K-1)*NFUND + L
C
              GXCH(IB1+J,JD1+L) = GXCH(IB1+J,JD1+L)
     &                          +         F2*RR(M, 5)*DENC(IC1+K,JA2+I)
     &                          +         F1*RR(M,15)*DENC(IC2+K,JA1+I) 
C
              GXCH(IB1+J,JD2+L) = GXCH(IB1+J,JD2+L)
     &                          +         F1*RR(M,14)*DENC(IC1+K,JA1+I) 
     &                          +         F2*RR(M, 6)*DENC(IC1+K,JA2+I)
     &                          +         F1*RR(M,16)*DENC(IC2+K,JA1+I) 
     &                          +         F2*RR(M, 8)*DENC(IC2+K,JA2+I)
C
              GXCH(IB2+J,JD1+L) = GXCH(IB2+J,JD1+L)
     &                          +         F2*RR(M, 9)*DENC(IC1+K,JA1+I) 
     &                          +         F1*RR(M, 1)*DENC(IC1+K,JA2+I)
     &                          +         F2*RR(M,11)*DENC(IC2+K,JA1+I) 
     &                          +         F1*RR(M, 3)*DENC(IC2+K,JA2+I)
C
              GXCH(IB2+J,JD2+L) = GXCH(IB2+J,JD2+L)
     &                          +         F1*RR(M, 2)*DENC(IC1+K,JA2+I)
     &                          +         F2*RR(M,12)*DENC(IC2+K,JA1+I) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       SEVENTH IFLG BATCH
        IF(IFLG(7).EQ.1) THEN
          DO L=1,NFUND
            DO K=1,NFUNC
              M = (K-1)*NFUND + L
C
              GXCH(ID1+L,JB1+J) = GXCH(ID1+L,JB1+J)
     &                          +         G2*RR(M, 2)*DENC(JA1+I,IC2+K)
     &                          +         G1*RR(M,12)*DENC(JA2+I,IC1+K) 
C
              GXCH(ID1+L,JB2+J) = GXCH(ID1+L,JB2+J)
     &                          +         G1*RR(M, 8)*DENC(JA1+I,IC1+K) 
     &                          +         G2*RR(M, 6)*DENC(JA1+I,IC2+K)
     &                          +         G1*RR(M,16)*DENC(JA2+I,IC1+K) 
     &                          +         G2*RR(M,14)*DENC(JA2+I,IC2+K)
C
              GXCH(ID2+L,JB1+J) = GXCH(ID2+L,JB1+J)
     &                          +         G2*RR(M, 3)*DENC(JA1+I,IC1+K) 
     &                          +         G1*RR(M, 1)*DENC(JA1+I,IC2+K)
     &                          +         G2*RR(M,11)*DENC(JA2+I,IC1+K) 
     &                          +         G1*RR(M, 9)*DENC(JA2+I,IC2+K)
C
              GXCH(ID2+L,JB2+J) = GXCH(ID2+L,JB2+J)
     &                          +         G1*RR(M, 5)*DENC(JA1+I,IC2+K)
     &                          +         G2*RR(M,15)*DENC(JA2+I,IC1+K) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       EIGHTH IFLG BATCH
        IF(IFLG(8).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GXCH(IB1+J,JC1+K) = GXCH(IB1+J,JC1+K)
     &                          +      F2*G1*RR(M, 8)*DENC(ID1+L,JA2+I)
     &                          +      F1*G2*RR(M,15)*DENC(ID2+L,JA1+I) 
C
              GXCH(IB1+J,JC2+K) = GXCH(IB1+J,JC2+K)
     &                          +      F1*G2*RR(M,14)*DENC(ID1+L,JA1+I) 
     &                          +      F2*G2*RR(M, 6)*DENC(ID1+L,JA2+I)
     &                          +      F1*G1*RR(M,13)*DENC(ID2+L,JA1+I) 
     &                          +      F2*G1*RR(M, 5)*DENC(ID2+L,JA2+I)
C
              GXCH(IB2+J,JC1+K) = GXCH(IB2+J,JC1+K)
     &                          +      F2*G1*RR(M,12)*DENC(ID1+L,JA1+I) 
     &                          +      F1*G1*RR(M, 4)*DENC(ID1+L,JA2+I)
     &                          +      F2*G2*RR(M,11)*DENC(ID2+L,JA1+I) 
     &                          +      F1*G2*RR(M, 3)*DENC(ID2+L,JA2+I)

C
              GXCH(IB2+J,JC2+K) = GXCH(IB2+J,JC2+K)
     &                          +      F1*G2*RR(M, 2)*DENC(ID1+L,JA2+I)
     &                          +      F2*G1*RR(M, 9)*DENC(ID2+L,JA1+I) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       NINTH IFLG BATCH
        IF(IFLG(9).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GXCH(IC1+K,JB1+J) = GXCH(IC1+K,JB1+J)
     &                          +            RR(M, 2)*DENC(JA1+I,ID2+L)
     &                          +            RR(M, 9)*DENC(JA2+I,ID1+L) 
C
              GXCH(IC1+K,JB2+J) = GXCH(IC1+K,JB2+J)
     &                          +            RR(M, 5)*DENC(JA1+I,ID1+L) 
     &                          +            RR(M, 6)*DENC(JA1+I,ID2+L)
     &                          +            RR(M,13)*DENC(JA2+I,ID1+L) 
     &                          +            RR(M,14)*DENC(JA2+I,ID2+L)
C
              GXCH(IC2+K,JB1+J) = GXCH(IC2+K,JB1+J)
     &                          +            RR(M, 3)*DENC(JA1+I,ID1+L) 
     &                          +            RR(M, 4)*DENC(JA1+I,ID2+L)
     &                          +            RR(M,11)*DENC(JA2+I,ID1+L) 
     &                          +            RR(M,12)*DENC(JA2+I,ID2+L)
C
              GXCH(IC2+K,JB2+J) = GXCH(IC2+K,JB2+J)
     &                          +            RR(M, 8)*DENC(JA1+I,ID2+L)
     &                          +            RR(M,15)*DENC(JA2+I,ID1+L) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       TENTH IFLG BATCH
        IF(IFLG(10).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GDIR(IA1+I,JB1+J) = GDIR(IA1+I,JB1+J)
     &                          +            RR(M, 2)*DENC(IC1+K,JD2+L)
     &                          +            RR(M, 3)*DENC(IC2+K,JD1+L) 
C
              GDIR(IA1+I,JB2+J) = GDIR(IA1+I,JB2+J)
     &                          +            RR(M, 5)*DENC(IC1+K,JD1+L) 
     &                          +            RR(M, 6)*DENC(IC1+K,JD2+L)
     &                          +            RR(M, 7)*DENC(IC2+K,JD1+L) 
     &                          +            RR(M, 8)*DENC(IC2+K,JD2+L)
C
              GDIR(IA2+I,JB1+J) = GDIR(IA2+I,JB1+J)
     &                          +            RR(M, 9)*DENC(IC1+K,JD1+L) 
     &                          +            RR(M,10)*DENC(IC1+K,JD2+L)
     &                          +            RR(M,11)*DENC(IC2+K,JD1+L) 
     &                          +            RR(M,12)*DENC(IC2+K,JD2+L)
C
              GDIR(IA2+I,JB2+J) = GDIR(IA2+I,JB2+J)
     &                          +            RR(M,14)*DENC(IC1+K,JD2+L)
     &                          +            RR(M,15)*DENC(IC2+K,JD1+L) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       ELEVENTH IFLG BATCH
        IF(IFLG(11).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              GDIR(IC1+K,JD1+L) = GDIR(IC1+K,JD1+L)
     &                          +            RR(M, 5)*DENC(IA1+I,JB2+J)
     &                          +            RR(M, 9)*DENC(IA2+I,JB1+J) 
C
              GDIR(IC1+K,JD2+L) = GDIR(IC1+K,JD2+L)
     &                          +            RR(M, 2)*DENC(IA1+I,JB1+J) 
     &                          +            RR(M, 6)*DENC(IA1+I,JB2+J)
     &                          +            RR(M,10)*DENC(IA2+I,JB1+J) 
     &                          +            RR(M,14)*DENC(IA2+I,JB2+J)
C
              GDIR(IC2+K,JD1+L) = GDIR(IC2+K,JD1+L)
     &                          +            RR(M, 3)*DENC(IA1+I,JB1+J) 
     &                          +            RR(M, 7)*DENC(IA1+I,JB2+J)
     &                          +            RR(M,11)*DENC(IA2+I,JB1+J) 
     &                          +            RR(M,15)*DENC(IA2+I,JB2+J)
C
              GDIR(IC2+K,JD2+L) = GDIR(IC2+K,JD2+L)
     &                          +            RR(M, 8)*DENC(IA1+I,JB2+J)
     &                          +            RR(M,12)*DENC(IA2+I,JB1+J) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       SKIP POINT FOR ATOMICS/DIATOMICS
390     CONTINUE
C
C*********************************************************************C
C     CONSTRUCT OPEN-SHELL COULOMB INTERACTION MATRICES, QDIR/QXCH.   C
C*********************************************************************C
C
        IF(NOPEN.EQ.0) GOTO 399
C
C       ATOMIC/DIATOMIC MATRIX ELEMENTS...
C
C       FIRST IFLG BATCH
        IF(IFLG(1).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QDIR(IA1+I,JB1+J) = QDIR(IA1+I,JB1+J)
     &                          +       ACFF*RR(M, 1)*DENO(IC1+K,JD1+L)
     &                          +       ACFF*RR(M, 4)*DENO(IC2+K,JD2+L)
     &                          +    G1*ACFF*RR(M, 4)*DENO(JD1+L,IC1+K)
     &                          +    G1*ACFF*RR(M, 1)*DENO(JD2+L,IC2+K)
C
              QDIR(IA2+I,JB2+J) = QDIR(IA2+I,JB2+J)
     &                          +       ACFF*RR(M,13)*DENO(IC1+K,JD1+L)
     &                          +       ACFF*RR(M,16)*DENO(IC2+K,JD2+L)
     &                          +    G1*ACFF*RR(M,16)*DENO(JD1+L,IC1+K)
     &                          +    G1*ACFF*RR(M,13)*DENO(JD2+L,IC2+K)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SECOND IFLG BATCH
        IF(IFLG(2).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QDIR(IC1+K,JD1+L) = QDIR(IC1+K,JD1+L)
     &                          +       ACFF*RR(M, 1)*DENO(IA1+I,JB1+J)
     &                          +       ACFF*RR(M,13)*DENO(IA2+I,JB2+J)
     &                          +    F1*ACFF*RR(M,13)*DENO(JB1+J,IA1+I)
     &                          +    F1*ACFF*RR(M, 1)*DENO(JB2+J,IA2+I)
C
              QDIR(IC2+K,JD2+L) = QDIR(IC2+K,JD2+L)
     &                          +       ACFF*RR(M, 4)*DENO(IA1+I,JB1+J)
     &                          +       ACFF*RR(M,16)*DENO(IA2+I,JB2+J)
     &                          +    F1*ACFF*RR(M,16)*DENO(JB1+J,IA1+I)
     &                          +    F1*ACFF*RR(M, 4)*DENO(JB2+J,IA2+I)
C
            ENDDO
          ENDDO
        ENDIF
C
C       THIRD IFLG BATCH          
        IF(IFLG(3).EQ.1) THEN
          DO L=1,NFUND
            DO K=1,NFUNC
              M = (K-1)*NFUND + L
C  
              QXCH(IA1+I,JD1+L) = QXCH(IA1+I,JD1+L)
     &                          +       BCFF*RR(M, 1)*DENO(IC1+K,JB1+J) 
     &                          +       BCFF*RR(M, 7)*DENO(IC2+K,JB2+J)
C
              QXCH(IA2+I,JD2+L) = QXCH(IA2+I,JD2+L)
     &                          +       BCFF*RR(M,10)*DENO(IC1+K,JB1+J) 
     &                          +       BCFF*RR(M,16)*DENO(IC2+K,JB2+J)
C
            ENDDO
          ENDDO
        ENDIF
C
C       FOURTH IFLG BATCH          
        IF(IFLG(4).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QXCH(IA1+I,JC1+K) = QXCH(IA1+I,JC1+K)
     &                          +    G1*BCFF*RR(M, 4)*DENO(ID1+L,JB1+J) 
     &                          +    G2*BCFF*RR(M, 7)*DENO(ID2+L,JB2+J)
C
              QXCH(IA2+I,JC2+K) = QXCH(IA2+I,JC2+K)
     &                          +    G2*BCFF*RR(M,10)*DENO(ID1+L,JB1+J) 
     &                          +    G1*BCFF*RR(M,13)*DENO(ID2+L,JB2+J)
C
            ENDDO
          ENDDO
        ENDIF
C
C       FIFTH IFLG BATCH
        IF(IFLG(5).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QXCH(IC1+K,JA1+I) = QXCH(IC1+K,JA1+I)
     &                          +    F1*BCFF*RR(M,13)*DENO(JB1+J,ID1+L)
     &                          +    F2*BCFF*RR(M,10)*DENO(JB2+J,ID2+L)
C     
              QXCH(IC2+K,JA2+I) = QXCH(IC2+K,JA2+I)
     &                          +    F2*BCFF*RR(M, 7)*DENO(JB1+J,ID1+L)
     &                          +    F1*BCFF*RR(M, 4)*DENO(JB2+J,ID2+L)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SIXTH IFLG BATCH
        IF(IFLG(6).EQ.1) THEN
          DO L=1,NFUND
            DO K=1,NFUNC
              M = (K-1)*NFUND + L
C
              QXCH(IB1+J,JD1+L) = QXCH(IB1+J,JD1+L)
     &                          +    F1*BCFF*RR(M,13)*DENO(IC1+K,JA1+I) 
     &                          +    F2*BCFF*RR(M, 7)*DENO(IC2+K,JA2+I)
C
              QXCH(IB2+J,JD2+L) = QXCH(IB2+J,JD2+L)
     &                          +    F2*BCFF*RR(M,10)*DENO(IC1+K,JA1+I) 
     &                          +    F1*BCFF*RR(M, 4)*DENO(IC2+K,JA2+I)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SEVENTH IFLG BATCH
        IF(IFLG(7).EQ.1) THEN
          DO L=1,NFUND
            DO K=1,NFUNC
              M = (K-1)*NFUND + L
C
              QXCH(ID1+L,JB1+J) = QXCH(ID1+L,JB1+J)
     &                          +    G1*BCFF*RR(M, 4)*DENO(JA1+I,IC1+K)
     &                          +    G2*BCFF*RR(M,10)*DENO(JA2+I,IC2+K)
C
              QXCH(ID2+L,JB2+J) = QXCH(ID2+L,JB2+J)
     &                          +    G2*BCFF*RR(M, 7)*DENO(JA1+I,IC1+K)
     &                          +    G1*BCFF*RR(M,13)*DENO(JA2+I,IC2+K)
C
            ENDDO
          ENDDO
        ENDIF
C
C       EIGHTH IFLG BATCH
        IF(IFLG(8).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QXCH(IB1+J,JC1+K) = QXCH(IB1+J,JC1+K)
     &                          + F1*G1*BCFF*RR(M,16)*DENO(ID1+L,JA1+I) 
     &                          + F2*G2*BCFF*RR(M, 7)*DENO(ID2+L,JA2+I)
C
              QXCH(IB2+J,JC2+K) = QXCH(IB2+J,JC2+K)
     &                          + F2*G2*BCFF*RR(M,10)*DENO(ID1+L,JA1+I) 
     &                          + F1*G1*BCFF*RR(M, 1)*DENO(ID2+L,JA2+I)
C
            ENDDO
          ENDDO
        ENDIF
C
C       NINTH IFLG BATCH
        IF(IFLG(9).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QXCH(IC1+K,JB1+J) = QXCH(IC1+K,JB1+J)
     &                          +       BCFF*RR(M, 1)*DENO(JA1+I,ID1+L)
     &                          +       BCFF*RR(M,10)*DENO(JA2+I,ID2+L)
C
              QXCH(IC2+K,JB2+J) = QXCH(IC2+K,JB2+J)
     &                          +       BCFF*RR(M, 7)*DENO(JA1+I,ID1+L)
     &                          +       BCFF*RR(M,16)*DENO(JA2+I,ID2+L)
C
            ENDDO
          ENDDO           
        ENDIF
C
C       TENTH IFLG BATCH
        IF(IFLG(10).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QDIR(IA1+I,JB1+J) = QDIR(IA1+I,JB1+J)
     &                          +       ACFF*RR(M, 1)*DENO(IC1+K,JD1+L) 
     &                          +       ACFF*RR(M, 4)*DENO(IC2+K,JD2+L)
C
              QDIR(IA2+I,JB2+J) = QDIR(IA2+I,JB2+J)
     &                          +       ACFF*RR(M,13)*DENO(IC1+K,JD1+L) 
     &                          +       ACFF*RR(M,16)*DENO(IC2+K,JD2+L)
C
            ENDDO
          ENDDO
        ENDIF
C
C       ELEVENTH IFLG BATCH
        IF(IFLG(11).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QDIR(IC1+K,JD1+L) = QDIR(IC1+K,JD1+L)
     &                          +       ACFF*RR(M, 1)*DENO(IA1+I,JB1+J) 
     &                          +       ACFF*RR(M,13)*DENO(IA2+I,JB2+J)
C
              QDIR(IC2+K,JD2+L) = QDIR(IC2+K,JD2+L)
     &                          +       ACFF*RR(M, 4)*DENO(IA1+I,JB1+J) 
     &                          +       ACFF*RR(M,16)*DENO(IA2+I,JB2+J)
C
            ENDDO
          ENDDO
        ENDIF
C
C       ADDITIONAL MATRIX ELEMENTS (ASSUMES NO GEOMETRIC SYMMETRY)...
c        IF(NCNT.LE.2) GOTO 395
C
C       FIRST IFLG BATCH
        IF(IFLG(1).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QDIR(IA1+I,JB1+J) = QDIR(IA1+I,JB1+J)
     &                          +       ACFF*RR(M, 2)*DENO(IC1+K,JD2+L)
     &                          +       ACFF*RR(M, 3)*DENO(IC2+K,JD1+L) 
     &                          +    G2*ACFF*RR(M, 2)*DENO(JD1+L,IC2+K)
     &                          +    G2*ACFF*RR(M, 3)*DENO(JD2+L,IC1+K) 
C
              QDIR(IA1+I,JB2+J) = QDIR(IA1+I,JB2+J)
     &                          +       ACFF*RR(M, 5)*DENO(IC1+K,JD1+L) 
     &                          +       ACFF*RR(M, 6)*DENO(IC1+K,JD2+L)
     &                          +       ACFF*RR(M, 7)*DENO(IC2+K,JD1+L) 
     &                          +       ACFF*RR(M, 8)*DENO(IC2+K,JD2+L)
     &                          +    G1*ACFF*RR(M, 8)*DENO(JD1+L,IC1+K) 
     &                          +    G2*ACFF*RR(M, 6)*DENO(JD1+L,IC2+K)
     &                          +    G2*ACFF*RR(M, 7)*DENO(JD2+L,IC1+K) 
     &                          +    G1*ACFF*RR(M, 5)*DENO(JD2+L,IC2+K)

C
              QDIR(IA2+I,JB1+J) = QDIR(IA2+I,JB1+J)
     &                          +       ACFF*RR(M, 9)*DENO(IC1+K,JD1+L) 
     &                          +       ACFF*RR(M,10)*DENO(IC1+K,JD2+L)
     &                          +       ACFF*RR(M,11)*DENO(IC2+K,JD1+L) 
     &                          +       ACFF*RR(M,12)*DENO(IC2+K,JD2+L)
     &                          +    G1*ACFF*RR(M,12)*DENO(JD1+L,IC1+K) 
     &                          +    G2*ACFF*RR(M,10)*DENO(JD1+L,IC2+K)
     &                          +    G2*ACFF*RR(M,11)*DENO(JD2+L,IC1+K) 
     &                          +    G1*ACFF*RR(M, 9)*DENO(JD2+L,IC2+K)
C
              QDIR(IA2+I,JB2+J) = QDIR(IA2+I,JB2+J)
     &                          +       ACFF*RR(M,14)*DENO(IC1+K,JD2+L)
     &                          +       ACFF*RR(M,15)*DENO(IC2+K,JD1+L) 
     &                          +    G2*ACFF*RR(M,14)*DENO(JD1+L,IC2+K)
     &                          +    G2*ACFF*RR(M,15)*DENO(JD2+L,IC1+K) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       SECOND IFLG BATCH
        IF(IFLG(2).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QDIR(IC1+K,JD1+L) = QDIR(IC1+K,JD1+L)
     &                          +       ACFF*RR(M, 5)*DENO(IA1+I,JB2+J)
     &                          +       ACFF*RR(M, 9)*DENO(IA2+I,JB1+J) 
     &                          +    F2*ACFF*RR(M, 5)*DENO(JB1+J,IA2+I)
     &                          +    F2*ACFF*RR(M, 9)*DENO(JB2+J,IA1+I) 
C
              QDIR(IC1+K,JD2+L) = QDIR(IC1+K,JD2+L)
     &                          +       ACFF*RR(M, 2)*DENO(IA1+I,JB1+J) 
     &                          +       ACFF*RR(M, 6)*DENO(IA1+I,JB2+J)
     &                          +       ACFF*RR(M,10)*DENO(IA2+I,JB1+J) 
     &                          +       ACFF*RR(M,14)*DENO(IA2+I,JB2+J)
     &                          +    F1*ACFF*RR(M,14)*DENO(JB1+J,IA1+I) 
     &                          +    F2*ACFF*RR(M, 6)*DENO(JB1+J,IA2+I)
     &                          +    F2*ACFF*RR(M,10)*DENO(JB2+J,IA1+I) 
     &                          +    F1*ACFF*RR(M, 2)*DENO(JB2+J,IA2+I)
C
              QDIR(IC2+K,JD1+L) = QDIR(IC2+K,JD1+L)
     &                          +       ACFF*RR(M, 3)*DENO(IA1+I,JB1+J) 
     &                          +       ACFF*RR(M, 7)*DENO(IA1+I,JB2+J)
     &                          +       ACFF*RR(M,11)*DENO(IA2+I,JB1+J) 
     &                          +       ACFF*RR(M,15)*DENO(IA2+I,JB2+J)
     &                          +    F1*ACFF*RR(M,15)*DENO(JB1+J,IA1+I) 
     &                          +    F2*ACFF*RR(M, 7)*DENO(JB1+J,IA2+I)
     &                          +    F2*ACFF*RR(M,11)*DENO(JB2+J,IA1+I) 
     &                          +    F1*ACFF*RR(M, 3)*DENO(JB2+J,IA2+I)
C
              QDIR(IC2+K,JD2+L) = QDIR(IC2+K,JD2+L)
     &                          +       ACFF*RR(M, 8)*DENO(IA1+I,JB2+J)
     &                          +       ACFF*RR(M,12)*DENO(IA2+I,JB1+J) 
     &                          +    F2*ACFF*RR(M, 8)*DENO(JB1+J,IA2+I)
     &                          +    F2*ACFF*RR(M,12)*DENO(JB2+J,IA1+I) 
C 
            ENDDO
          ENDDO
        ENDIF
C
C       THIRD IFLG BATCH          
        IF(IFLG(3).EQ.1) THEN
          DO L=1,NFUND
            DO K=1,NFUNC
              M = (K-1)*NFUND + L
C  
              QXCH(IA1+I,JD1+L) = QXCH(IA1+I,JD1+L)
     &                          +       BCFF*RR(M, 5)*DENO(IC1+K,JB2+J)
     &                          +       BCFF*RR(M, 3)*DENO(IC2+K,JB1+J) 
C
              QXCH(IA1+I,JD2+L) = QXCH(IA1+I,JD2+L)
     &                          +       BCFF*RR(M, 2)*DENO(IC1+K,JB1+J) 
     &                          +       BCFF*RR(M, 6)*DENO(IC1+K,JB2+J)
     &                          +       BCFF*RR(M, 4)*DENO(IC2+K,JB1+J) 
     &                          +       BCFF*RR(M, 8)*DENO(IC2+K,JB2+J)
C
              QXCH(IA2+I,JD1+L) = QXCH(IA2+I,JD1+L)
     &                          +       BCFF*RR(M, 9)*DENO(IC1+K,JB1+J) 
     &                          +       BCFF*RR(M,13)*DENO(IC1+K,JB2+J)
     &                          +       BCFF*RR(M,11)*DENO(IC2+K,JB1+J) 
     &                          +       BCFF*RR(M,15)*DENO(IC2+K,JB2+J)
C
              QXCH(IA2+I,JD2+L) = QXCH(IA2+I,JD2+L)
     &                          +       BCFF*RR(M,14)*DENO(IC1+K,JB2+J)
     &                          +       BCFF*RR(M,12)*DENO(IC2+K,JB1+J) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       FOURTH IFLG BATCH          
        IF(IFLG(4).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QXCH(IA1+I,JC1+K) = QXCH(IA1+I,JC1+K)
     &                          +    G1*BCFF*RR(M, 8)*DENO(ID1+L,JB2+J)
     &                          +    G2*BCFF*RR(M, 3)*DENO(ID2+L,JB1+J)
C
              QXCH(IA1+I,JC2+K) = QXCH(IA1+I,JC2+K)
     &                          +    G2*BCFF*RR(M, 2)*DENO(ID1+L,JB1+J) 
     &                          +    G2*BCFF*RR(M, 6)*DENO(ID1+L,JB2+J)
     &                          +    G1*BCFF*RR(M, 1)*DENO(ID2+L,JB1+J) 
     &                          +    G1*BCFF*RR(M, 5)*DENO(ID2+L,JB2+J)
C
              QXCH(IA2+I,JC1+K) = QXCH(IA2+I,JC1+K)
     &                          +    G1*BCFF*RR(M,12)*DENO(ID1+L,JB1+J) 
     &                          +    G1*BCFF*RR(M,16)*DENO(ID1+L,JB2+J)
     &                          +    G2*BCFF*RR(M,11)*DENO(ID2+L,JB1+J) 
     &                          +    G2*BCFF*RR(M,15)*DENO(ID2+L,JB2+J)
C
              QXCH(IA2+I,JC2+K) = QXCH(IA2+I,JC2+K)
     &                          +    G2*BCFF*RR(M,14)*DENO(ID1+L,JB2+J)
     &                          +    G1*BCFF*RR(M, 9)*DENO(ID2+L,JB1+J) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       FIFTH IFLG BATCH
        IF(IFLG(5).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QXCH(IC1+K,JA1+I) = QXCH(IC1+K,JA1+I)
     &                          +    F1*BCFF*RR(M,14)*DENO(JB1+J,ID2+L)
     &                          +    F2*BCFF*RR(M, 9)*DENO(JB2+J,ID1+L) 
C
              QXCH(IC1+K,JA2+I) = QXCH(IC1+K,JA2+I)
     &                          +    F2*BCFF*RR(M, 5)*DENO(JB1+J,ID1+L) 
     &                          +    F2*BCFF*RR(M, 6)*DENO(JB1+J,ID2+L)
     &                          +    F1*BCFF*RR(M, 1)*DENO(JB2+J,ID1+L) 
     &                          +    F1*BCFF*RR(M, 2)*DENO(JB2+J,ID2+L)
C
              QXCH(IC2+K,JA1+I) = QXCH(IC2+K,JA1+I)
     &                          +    F1*BCFF*RR(M,15)*DENO(JB1+J,ID1+L) 
     &                          +    F1*BCFF*RR(M,16)*DENO(JB1+J,ID2+L)
     &                          +    F2*BCFF*RR(M,11)*DENO(JB2+J,ID1+L) 
     &                          +    F2*BCFF*RR(M,12)*DENO(JB2+J,ID2+L)
C
              QXCH(IC2+K,JA2+I) = QXCH(IC2+K,JA2+I)
     &                          +    F2*BCFF*RR(M, 8)*DENO(JB1+J,ID2+L)
     &                          +    F1*BCFF*RR(M, 3)*DENO(JB2+J,ID1+L) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       SIXTH IFLG BATCH
        IF(IFLG(6).EQ.1) THEN
          DO L=1,NFUND
            DO K=1,NFUNC
              M = (K-1)*NFUND + L
C
              QXCH(IB1+J,JD1+L) = QXCH(IB1+J,JD1+L)
     &                          +    F2*BCFF*RR(M, 5)*DENO(IC1+K,JA2+I)
     &                          +    F1*BCFF*RR(M,15)*DENO(IC2+K,JA1+I) 
C
              QXCH(IB1+J,JD2+L) = QXCH(IB1+J,JD2+L)
     &                          +    F1*BCFF*RR(M,14)*DENO(IC1+K,JA1+I) 
     &                          +    F2*BCFF*RR(M, 6)*DENO(IC1+K,JA2+I)
     &                          +    F1*BCFF*RR(M,16)*DENO(IC2+K,JA1+I) 
     &                          +    F2*BCFF*RR(M, 8)*DENO(IC2+K,JA2+I)
C
              QXCH(IB2+J,JD1+L) = QXCH(IB2+J,JD1+L)
     &                          +    F2*BCFF*RR(M, 9)*DENO(IC1+K,JA1+I) 
     &                          +    F1*BCFF*RR(M, 1)*DENO(IC1+K,JA2+I)
     &                          +    F2*BCFF*RR(M,11)*DENO(IC2+K,JA1+I) 
     &                          +    F1*BCFF*RR(M, 3)*DENO(IC2+K,JA2+I)
C
              QXCH(IB2+J,JD2+L) = QXCH(IB2+J,JD2+L)
     &                          +    F1*BCFF*RR(M, 2)*DENO(IC1+K,JA2+I)
     &                          +    F2*BCFF*RR(M,12)*DENO(IC2+K,JA1+I) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       SEVENTH IFLG BATCH
        IF(IFLG(7).EQ.1) THEN
          DO L=1,NFUND
            DO K=1,NFUNC
              M = (K-1)*NFUND + L
C
              QXCH(ID1+L,JB1+J) = QXCH(ID1+L,JB1+J)
     &                          +    G2*BCFF*RR(M, 2)*DENO(JA1+I,IC2+K)
     &                          +    G1*BCFF*RR(M,12)*DENO(JA2+I,IC1+K) 
C
              QXCH(ID1+L,JB2+J) = QXCH(ID1+L,JB2+J)
     &                          +    G1*BCFF*RR(M, 8)*DENO(JA1+I,IC1+K) 
     &                          +    G2*BCFF*RR(M, 6)*DENO(JA1+I,IC2+K)
     &                          +    G1*BCFF*RR(M,16)*DENO(JA2+I,IC1+K) 
     &                          +    G2*BCFF*RR(M,14)*DENO(JA2+I,IC2+K)
C
              QXCH(ID2+L,JB1+J) = QXCH(ID2+L,JB1+J)
     &                          +    G2*BCFF*RR(M, 3)*DENO(JA1+I,IC1+K) 
     &                          +    G1*BCFF*RR(M, 1)*DENO(JA1+I,IC2+K)
     &                          +    G2*BCFF*RR(M,11)*DENO(JA2+I,IC1+K) 
     &                          +    G1*BCFF*RR(M, 9)*DENO(JA2+I,IC2+K)
C
              QXCH(ID2+L,JB2+J) = QXCH(ID2+L,JB2+J)
     &                          +    G1*BCFF*RR(M, 5)*DENO(JA1+I,IC2+K)
     &                          +    G2*BCFF*RR(M,15)*DENO(JA2+I,IC1+K) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       EIGHTH IFLG BATCH
        IF(IFLG(8).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QXCH(IB1+J,JC1+K) = QXCH(IB1+J,JC1+K)
     &                          + F2*G1*BCFF*RR(M, 8)*DENO(ID1+L,JA2+I)
     &                          + F1*G2*BCFF*RR(M,15)*DENO(ID2+L,JA1+I) 
C
              QXCH(IB1+J,JC2+K) = QXCH(IB1+J,JC2+K)
     &                          + F1*G2*BCFF*RR(M,14)*DENO(ID1+L,JA1+I) 
     &                          + F2*G2*BCFF*RR(M, 6)*DENO(ID1+L,JA2+I)
     &                          + F1*G1*BCFF*RR(M,13)*DENO(ID2+L,JA1+I) 
     &                          + F2*G1*BCFF*RR(M, 5)*DENO(ID2+L,JA2+I)
C
              QXCH(IB2+J,JC1+K) = QXCH(IB2+J,JC1+K)
     &                          + F2*G1*BCFF*RR(M,12)*DENO(ID1+L,JA1+I) 
     &                          + F1*G1*BCFF*RR(M, 4)*DENO(ID1+L,JA2+I)
     &                          + F2*G2*BCFF*RR(M,11)*DENO(ID2+L,JA1+I) 
     &                          + F1*G2*BCFF*RR(M, 3)*DENO(ID2+L,JA2+I)

C
              QXCH(IB2+J,JC2+K) = QXCH(IB2+J,JC2+K)
     &                          + F1*G2*BCFF*RR(M, 2)*DENO(ID1+L,JA2+I)
     &                          + F2*G1*BCFF*RR(M, 9)*DENO(ID2+L,JA1+I) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       NINTH IFLG BATCH
        IF(IFLG(9).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QXCH(IC1+K,JB1+J) = QXCH(IC1+K,JB1+J)
     &                          +       BCFF*RR(M, 2)*DENO(JA1+I,ID2+L)
     &                          +       BCFF*RR(M, 9)*DENO(JA2+I,ID1+L) 
C
              QXCH(IC1+K,JB2+J) = QXCH(IC1+K,JB2+J)
     &                          +       BCFF*RR(M, 5)*DENO(JA1+I,ID1+L) 
     &                          +       BCFF*RR(M, 6)*DENO(JA1+I,ID2+L)
     &                          +       BCFF*RR(M,13)*DENO(JA2+I,ID1+L) 
     &                          +       BCFF*RR(M,14)*DENO(JA2+I,ID2+L)
C
              QXCH(IC2+K,JB1+J) = QXCH(IC2+K,JB1+J)
     &                          +       BCFF*RR(M, 3)*DENO(JA1+I,ID1+L) 
     &                          +       BCFF*RR(M, 4)*DENO(JA1+I,ID2+L)
     &                          +       BCFF*RR(M,11)*DENO(JA2+I,ID1+L) 
     &                          +       BCFF*RR(M,12)*DENO(JA2+I,ID2+L)
C
              QXCH(IC2+K,JB2+J) = QXCH(IC2+K,JB2+J)
     &                          +       BCFF*RR(M, 8)*DENO(JA1+I,ID2+L)
     &                          +       BCFF*RR(M,15)*DENO(JA2+I,ID1+L) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       TENTH IFLG BATCH
        IF(IFLG(10).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QDIR(IA1+I,JB1+J) = QDIR(IA1+I,JB1+J)
     &                          +       ACFF*RR(M, 2)*DENO(IC1+K,JD2+L)
     &                          +       ACFF*RR(M, 3)*DENO(IC2+K,JD1+L) 
C
              QDIR(IA1+I,JB2+J) = QDIR(IA1+I,JB2+J)
     &                          +       ACFF*RR(M, 5)*DENO(IC1+K,JD1+L) 
     &                          +       ACFF*RR(M, 6)*DENO(IC1+K,JD2+L)
     &                          +       ACFF*RR(M, 7)*DENO(IC2+K,JD1+L) 
     &                          +       ACFF*RR(M, 8)*DENO(IC2+K,JD2+L)
C
              QDIR(IA2+I,JB1+J) = QDIR(IA2+I,JB1+J)
     &                          +       ACFF*RR(M, 9)*DENO(IC1+K,JD1+L) 
     &                          +       ACFF*RR(M,10)*DENO(IC1+K,JD2+L)
     &                          +       ACFF*RR(M,11)*DENO(IC2+K,JD1+L) 
     &                          +       ACFF*RR(M,12)*DENO(IC2+K,JD2+L)
C
              QDIR(IA2+I,JB2+J) = QDIR(IA2+I,JB2+J)
     &                          +       ACFF*RR(M,14)*DENO(IC1+K,JD2+L)
     &                          +       ACFF*RR(M,15)*DENO(IC2+K,JD1+L) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       ELEVENTH IFLG BATCH
        IF(IFLG(11).EQ.1) THEN
          M = 0
          DO K=1,NFUNC
            DO L=1,NFUND
              M = M+1
C
              QDIR(IC1+K,JD1+L) = QDIR(IC1+K,JD1+L)
     &                          +       ACFF*RR(M, 5)*DENO(IA1+I,JB2+J)
     &                          +       ACFF*RR(M, 9)*DENO(IA2+I,JB1+J) 
C
              QDIR(IC1+K,JD2+L) = QDIR(IC1+K,JD2+L)
     &                          +       ACFF*RR(M, 2)*DENO(IA1+I,JB1+J) 
     &                          +       ACFF*RR(M, 6)*DENO(IA1+I,JB2+J)
     &                          +       ACFF*RR(M,10)*DENO(IA2+I,JB1+J) 
     &                          +       ACFF*RR(M,14)*DENO(IA2+I,JB2+J)
C
              QDIR(IC2+K,JD1+L) = QDIR(IC2+K,JD1+L)
     &                          +       ACFF*RR(M, 3)*DENO(IA1+I,JB1+J) 
     &                          +       ACFF*RR(M, 7)*DENO(IA1+I,JB2+J)
     &                          +       ACFF*RR(M,11)*DENO(IA2+I,JB1+J) 
     &                          +       ACFF*RR(M,15)*DENO(IA2+I,JB2+J)
C
              QDIR(IC2+K,JD2+L) = QDIR(IC2+K,JD2+L)
     &                          +       ACFF*RR(M, 8)*DENO(IA1+I,JB2+J)
     &                          +       ACFF*RR(M,12)*DENO(IA2+I,JB1+J) 
C
            ENDDO
          ENDDO
        ENDIF
C
C       SKIP POINT FOR ATOMICS/DIATOMICS
395     CONTINUE
C
399   CONTINUE
C
C     CLOSE ALL LOOPS OVER BASIS FUNCTIONS A,B,C,D
3000  CONTINUE
2999  CONTINUE
2500  CONTINUE
2200  CONTINUE
2000  CONTINUE
1000  CONTINUE
C
C*********************************************************************C
C     COMPLETE CONSTRUCTION OF ALL MATRICES BY CONJUGATION.           C
C*********************************************************************C
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
C*********************************************************************C
C     SPECIFY OPEN-SHELL COUPLING MATRIX R AND ADD TO COULOMB MATRIX  C
C ------------------------------------------------------------------- C
C               R = [S.D(O).Q + Q.D(O).S]  (RSCF 89)                  C
C*********************************************************************C
C         
      DO J=1,NDIM
        DO K=1,NDIM
          T1(K) = 0.0D0
          T2(K) = 0.0D0
          DO L=1,NDIM
            T1(K) = T1(K) + DENT(K,L)*QDIR(L,J) - DENT(K,L)*QDIR(L,J)
            T2(K) = T2(K) + DENT(K,L)*OVAP(L,J)
          ENDDO
        ENDDO
C
        DO I=1,NDIM
          T3(I) = 0.0D0
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
      RETURN
      END

