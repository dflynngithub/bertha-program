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
      PARAMETER (MXDM=1600,NCENTM=6,MAXB=26,MAXB2=MAXB*MAXB,NKAPM=9,
     &           LMAX=(NKAPM-1)/2,LAMMX=2*LMAX,IL4=2*LAMMX,MAXMV=NKAPM,
     &           MAXMV2=2*MAXMV,MLL=((NKAPM+8)*(NKAPM+9)*(NKAPM+10)/6))
C
      CHARACTER*4 TREE
C
      COMPLEX*16 C(MXDM,MXDM)
      COMPLEX*16 DENC(MXDM,MXDM),DENO(MXDM,MXDM),DENT(MXDM,MXDM)
      COMPLEX*16 FOCK(MXDM,MXDM),OVAP(MXDM,MXDM),QMAT(MXDM,MXDM)
      COMPLEX*16 RR(MAXB2,16)
C
      COMMON/COEFF/C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/ENERGY/ETOT,ENUC,EONE,ECLM,EBRT
      COMMON/FMAT/FOCK,OVAP,QMAT
      COMMON/IDECID/TREE
      COMMON/LABELS/LARGE(NCENTM,NKAPM,MAXMV2)
      COMMON/LBL2/LDIAG(100),NDIG
      COMMON/LOCGEO/XYZ(3,4)
      COMMON/SHELL/ACFF,BCFF,FOPEN,ICLOSE(100),IOPEN(6),NCLOSE,NOPEN
      COMMON/SPEC/COORD(3,NCENTM),ZNUC(NCENTM),CNUC(NCENTM),
     &           EXPSET(MAXB,NKAPM,NCENTM),IZNUC(NCENTM),ICRGE(NCENTM),
     &           KVALS(NKAPM,NCENTM),NKAPPA(NCENTM),LMAXX(NCENTM),
     &           NFUNCT(NKAPM,NCENTM),NCENT,NDIM,NSHIFT,
     &           NOCC,ITER,IALL,IRUN,IOCCM0
      COMMON/SYMARR/KAPLAB(MXDM),ICNLAB(MXDM),IMLAB(MXDM)
C
      DIMENSION EXPT(MAXB,4),KQN(4),MQN(4),NFUNS(4),LQN(4)
      DIMENSION ITQN(2),IFLG(11),ISCF(11,6)
      DIMENSION INDEX(NCENTM,-NKAPM:NKAPM,2*(NKAPM+1)*NKAPM)
      DIMENSION T1(MXDM),T2(MXDM),T3(MXDM)
C
C     ISCF TELLS WHICH INTEGRALS TO INCLUDE BASED ON OVERLAP COMBINATION
      DATA ISCF/1,1,1,1,0,1,0,1,0,0,0,
     &          1,0,1,1,0,1,0,0,0,0,0,
     &          0,1,1,0,0,1,0,0,0,1,0,
     &          1,0,1,1,0,0,0,0,0,0,1,
     &          0,0,1,0,0,0,0,0,0,1,1,
     &          0,0,1,0,0,0,0,0,0,1,0/
C
C     INITIALISE TWO-ELECTRON COULOMB ENERGY COUNTER
      ECLM  = 0.0D0
C
C*********************************************************************C
C     INDEXING ROUTINE: SO WE CAN SET BASIS FUNCTION LABELS BY BLOCK  C
C*********************************************************************C
      ICOUNT = 0
      LAMMAX = 0  
C
C     LOOP OVER NUCLEAR CENTRES
      DO ICENT=1,NCENT
C
        IF(LMAXX(ICENT).GT.LAMMAX) LAMMAX = LMAXX(ICENT)
C
C       LOOP OVER KAPPA VALUES FOR THIS NUCLEAR CENTRE
        DO KN=1,NKAPPA(ICENT)
C
C         IMPORT KAPPA, MAXIMUM MQN
          KAPPA = KVALS(KN,ICENT)
          MJMAX  = 2*IABS(KAPPA)-1
C
C         LOOP OVER MQN VALUES AND RECORD INDEX
          DO MJ=1,MJMAX,2
            ICOUNT = ICOUNT+1
            INDEX(ICENT,KAPPA,MJ) = ICOUNT
          ENDDO
        ENDDO
      ENDDO
C
C*********************************************************************C
C     GENERATE SOME BOYS INTEGRALS FOR UPCOMING RMAKE/FUNFX ROUTINES  C
C*********************************************************************C
C
       CALL GFINIT(4*LAMMAX+10)       
C
C*********************************************************************C
C     SET THE LIMITS FOR RUNNING OVER DENSITY COMBINATIONS            C
C*********************************************************************C
C
      IF(TREE.EQ.'NORL') THEN
        ITSTRT = 1
        ITSTOP = 1
        ITSKIP = 1
      ELSEIF(TREE.EQ.'DHFR'.OR.TREE.EQ.'DHFB') THEN
        ITSTRT = 1
        ITSTOP = 4
        ITSKIP = 3
      ENDIF
C
C*********************************************************************C
C     FIRST LAYER OF LOOPS FOR COMPONENTS IN TREE (USE INDEX 1000)    C
C*********************************************************************C
C
C     SELECT COMPONENT TYPES: NORL = (LL|LL) 
C                             DHF  = (LL|LL),(LL|SS),(SS|LL),(SS|SS)
      DO 1000 IT1=ITSTRT,ITSTOP,ITSKIP
      DO 1000 IT2=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITQN(1) = IT1
        ITQN(2) = IT2
C
C       GENERATE STARTING ADDRESSES FOR THIS IT1, IT2 COMBINATION
        IF(IT1.EQ.1) THEN
          NADDAB = 0
        ELSE
          NADDAB = NSHIFT
        ENDIF
C
        IF(IT2.EQ.1) THEN
          NADDCD = 0
        ELSE
          NADDCD = NSHIFT
        ENDIF
C
C*********************************************************************C
C     SECOND LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)    C
C*********************************************************************C
C
C     LOOP OVER CENTRE A
      DO 2000 ICENTA=1,NCENT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = COORD(1,ICENTA)
        XYZ(2,1) = COORD(2,ICENTA)
        XYZ(3,1) = COORD(3,ICENTA)
C
C     LOOP OVER CENTRE B
      DO 2000 ICENTB=1,ICENTA
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = COORD(1,ICENTB)
        XYZ(2,2) = COORD(2,ICENTB)
        XYZ(3,2) = COORD(3,ICENTB)
C
C       PARAMETER FOR SINGLE-CENTRE/MULTI-CENTRE OVERLAP OVER A AND B
        IF(ICENTA.EQ.ICENTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 2
        ENDIF
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAPPA(ICENTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICENTA)       
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNA    = NFUNCT(LQN(1)+1,ICENTA)
        NFUNS(1) = NFUNA
C
        DO IBAS=1,NFUNA
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICENTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAPPA(ICENTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2) = KVALS(KB,ICENTB)       
        IF(KQN(2).GT.0) THEN
         LQN(2) = KQN(2)
        ELSE
         LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNB    = NFUNCT(LQN(2)+1,ICENTB)
        NFUNS(2) = NFUNB
C
        DO IBAS=1,NFUNB
          EXPT(IBAS,2) = EXPSET(IBAS,LQN(2)+1,ICENTB)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 2000 MA=1,IABS(KQN(1))
        MJA    = 2*MA-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 2000 MB=1,IABS(KQN(2))
        MJB    = 2*MB-1
        MQN(2) = MJB
C
C     CALCULATE NEW BLOCK OF E(AB) COEFFS AT NEXT OPPORTUNITY
      IEAB = 1
C
C*********************************************************************C
C      THIRD LAYER OF LOOPS, OVER CENTRES C AND D (USE INDEX 3000)    C
C*********************************************************************C
C     LOOP OVER CENTRE C
      DO 3000 ICENTC=1,NCENT
C
C       CARTESIAN COORDINATES OF CENTRE C
        XYZ(1,3) = COORD(1,ICENTC)
        XYZ(2,3) = COORD(2,ICENTC)
        XYZ(3,3) = COORD(3,ICENTC)
C
C     LOOP OVER CENTRE D
      DO 3000 ICENTD=1,NCENT
C
C       CARTESIAN COORDINATES OF CENTRE D
        XYZ(1,4) = COORD(1,ICENTD)
        XYZ(2,4) = COORD(2,ICENTD)
        XYZ(3,4) = COORD(3,ICENTD)
C
C       PARAMETER FOR SINGLE-CENTRE/MULTI-CENTRE OVERLAP OVER C AND D
        IF(ICENTC.EQ.ICENTD) THEN
          INUCCD = 1
        ELSE
          INUCCD = 2
        ENDIF
C
C       PARAMETER FOR ATOMIC OR MULTICENTRE INTEGRAL
        IF(INUCAB*INUCCD.EQ.1.AND.ICENTA.EQ.ICENTC) THEN
          IATOM = 1
        ELSE
          IATOM = 0
        ENDIF
C
C     LOOP OVER KQN(C) VALUES
      DO 3000 KC=1,NKAPPA(ICENTC)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR C
        KQN(3) = KVALS(KC,ICENTC)
        IF(KQN(3).GT.0) THEN
          LQN(3) = KQN(3)
        ELSE
          LQN(3) =-KQN(3)-1
        ENDIF
C         
        NFUNC    = NFUNCT(LQN(3)+1,ICENTC)
        NFUNS(3) = NFUNC
C
        DO IBAS=1,NFUNC
          EXPT(IBAS,3) = EXPSET(IBAS,LQN(3)+1,ICENTC)
        ENDDO
C
C     LOOP OVER KQN(D) VALUES
      DO 3000 KD=1,NKAPPA(ICENTD)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR D
        KQN(4) = KVALS(KD,ICENTD)
        IF(KQN(4).GT.0) THEN
          LQN(4) = KQN(4)
        ELSE
          LQN(4) =-KQN(4)-1
        ENDIF
C
        NFUND    = NFUNCT(LQN(4)+1,ICENTD)
        NFUNS(4) = NFUND
C
        DO IBAS=1,NFUND
          EXPT(IBAS,4) = EXPSET(IBAS,LQN(4)+1,ICENTD)
        ENDDO
C
C     LOOP OVER |MQN(C)| VALUES
      DO 3000 MC=1,IABS(KQN(3))
        MJC    = 2*MC-1
        MQN(3) = MJC
C
C     LOOP OVER |MQN(D)| VALUES
      DO 3000 MD=1,IABS(KQN(4))
        MJD    = 2*MD-1
        MQN(4) = MJD
C
C
C       STAGES: DECISION TREE FOR SKIPPING MULTI-CENTRE CONTRIBUTIONS
        IF(IATOM.EQ.0) THEN
C
C >>      STAGE 0: INCLUDE ONLY (LL|LL) REPULSION INTEGRALS
          IF(IALL.EQ.0.AND.IT1+IT2.GT.2) THEN
            GOTO 2999
          ENDIF
C
C >>      STAGE 1: INCLUDE ONLY (LL|SS) AND (SS|LL) REPULSION INTEGRALS
          IF(IALL.EQ.1.AND.IT1+IT2.GT.5) THEN
            GOTO 2999
          ENDIF
C
        ENDIF
C
C     CALCULATE NEW BLOCK OF E(AB) COEFFS AT NEXT OPPORTUNITY
      IECD = 1
C
C*********************************************************************C
C     FOR THIS CHOICE OF A,B,C AND D, COMPUTE ADDRESSES AND PHASES    C
C*********************************************************************C
C
C     CALCULATE BLOCK INDICES FOR {ABCD} COMBINATIONS
      IQ1 = INDEX(ICENTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICENTB,KQN(2),MQN(2))
      IQ3 = INDEX(ICENTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICENTD,KQN(4),MQN(4))
C
      IQ12 = (IQ1*(IQ1-1))/2 + IQ2
      IQ34 = (IQ3*(IQ3-1))/2 + IQ4
C
C     FURTHER DEFINE STARTING ADDRESSES FOR {ABCD} BASIS LABELS
      IA1 = LARGE(ICENTA,KA,MJA)   + NADDAB
      IB1 = LARGE(ICENTB,KB,MJB)   + NADDAB
      IC1 = LARGE(ICENTC,KC,MJC)   + NADDCD
      ID1 = LARGE(ICENTD,KD,MJD)   + NADDCD
C
      IA2 = LARGE(ICENTA,KA,MJA+1) + NADDAB
      IB2 = LARGE(ICENTB,KB,MJB+1) + NADDAB
      IC2 = LARGE(ICENTC,KC,MJC+1) + NADDCD
      ID2 = LARGE(ICENTD,KD,MJD+1) + NADDCD
C
      JA1 = LARGE(ICENTA,KA,MJA)   + NADDAB
      JB1 = LARGE(ICENTB,KB,MJB)   + NADDAB
      JC1 = LARGE(ICENTC,KC,MJC)   + NADDCD
      JD1 = LARGE(ICENTD,KD,MJD)   + NADDCD
C
      JA2 = LARGE(ICENTA,KA,MJA+1) + NADDAB
      JB2 = LARGE(ICENTB,KB,MJB+1) + NADDAB
      JC2 = LARGE(ICENTC,KC,MJC+1) + NADDCD
      JD2 = LARGE(ICENTD,KD,MJD+1) + NADDCD
C
C     CALCULATE KQN PHASE FACTORS FOR PERMUTING INTEGRALS
      IF((KQN(1)*KQN(2)).GT.0) THEN 
        RKFAC1 = 1.0D0
      ELSE
        RKFAC1 =-1.0D0
      ENDIF
C        
      IF((KQN(3)*KQN(4)).GT.0) THEN 
        RKFAC2 = 1.0D0
      ELSE
        RKFAC2 =-1.0D0
      ENDIF
C
C     CALCULATE MQN PHASE FACTORS FOR PERMUTING INTEGRALS
      FACAB1 = DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      FACAB2 = DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      FACCD1 = DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      FACCD2 = DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
C*********************************************************************C
C     SKIP BATCHES WHICH CONFORM TO INTEGRAL SYMMETRIES               C
C*********************************************************************C
C
      IF(NCENT.LT.3) THEN
        IF (MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) GOTO 3998
        IF (MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) GOTO 3998
        IF (MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) GOTO 3998
        GOTO 3999
      ENDIF
3998  CONTINUE

C     DECISION TREE FOR SKIPPING CONTRIBUTIONS DUE TO INTEGRAL SYMMETRY
      IF (IQ1.LT.IQ2)   GOTO 3999
      IF (IQ3.LT.IQ4)   GOTO 3999
      IF (IQ12.LT.IQ34) GOTO 3999
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
        GOTO 3999
      ENDIF
      
      WRITE(*,*) 'IT,IQ',IT1,IT2,IQ1,IQ2,IQ3,IQ4
C
C     READ IN FLAG VALUES FROM ISCF DATA BLOCK
      DO M=1,11
        IFLG(M) = ISCF(M,ITSCF)
      ENDDO
C
C     INCLUDE SPECIAL CASES
      IF(ITSCF.EQ.1.AND.IQ1.EQ.IQ3) IFLG(5) = 1
      IF(ITSCF.EQ.1.AND.IQ2.EQ.IQ4) IFLG(7) = 1
      IF(ITSCF.EQ.1.AND.IQ2.EQ.IQ3) IFLG(9) = 1
      IF(ITSCF.EQ.3.AND.IQ2.EQ.IQ3) IFLG(9) = 1
      IF(ITSCF.EQ.4.AND.IQ2.EQ.IQ3) IFLG(9) = 1
C
C
C*********************************************************************C
C     FOURTH LAYER OF LOOPS, OVER BASIS FUNCTIONS A AND B (4000)      C
C ------------------------------------------------------------------- C
C     THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH        C
C     GENERATE THE FOCK MATRIX FROM THE SPINOR INTEGRALS.             C
C     THESE INCLUDE IMPLICIT PHASE FACTORS FOR THE PERMUTATION OF     C
C     KQN(1) <-> KQN(2) AND MQN(1) <-> MQN(2)                         C
C ------------------------------------------------------------------- C
C     (RSCF 86, 87)                                                   C
C*********************************************************************C  
C
C     NORMALLY WE WOULD USE IBAS AND JBAS, BUT SHORTEN TO I AND J FOR SPACE
      DO 4000 I=1,NFUNA
      DO 4000 J=1,NFUNB
C
C       GENERATE BATCH OF ELECTRON REPULSION INTEGRALS
        CALL ERI(RR,XYZ,KQN,MQN,EXPT,NFUNS,ITQN,I,J,IEAB,IECD)
C
C       THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH
C       GENERATE THE FOCK MATRIX FROM THE SPINOR INTEGRALS.
C
        IF(NCENT.LE.2) THEN
C       DIATOMIC CASE:
C         FIRST IFLG BATCH
          IF(IFLG(1).EQ.1) THEN
            F1 = RKFAC2*FACCD1
            F2 = RKFAC2*FACCD2
            M  = 0
            DO K=1,NFUNC
              DO L=1,NFUND
                M = M+1
C
                FOCK(IA1+I,JB1+J) = FOCK(IA1+I,JB1+J)
     &            + RR(M,1 )*(DENC(IC1+K,JD1+L) + F1*DENC(IC2+K,JD2+L))
     &            + RR(M,4 )*(DENC(IC1+K,JD1+L) + F1*DENC(IC2+K,JD2+L))
C
                FOCK(IA2+I,JB2+J) = FOCK(IA2+I,JB2+J)
     &            + RR(M,13)*(DENC(IC1+K,JD1+L) + F1*DENC(IC2+K,JD2+L))
     &            + RR(M,16)*(DENC(IC1+K,JD1+L) + F1*DENC(IC2+K,JD2+L))
              ENDDO
            ENDDO
          ENDIF
C
C         SECOND IFLG BATCH
          IF(IFLG(2).EQ.1) THEN
            F1 = RKFAC1*FACAB1
            F2 = RKFAC1*FACAB2
            M  = 0
            DO K=1,NFUNC
              DO L=1,NFUND
                M = M+1
C
                FOCK(IC1+K,JD1+L) = FOCK(IC1+K,JD1+L)
     &            + RR(M,1 )*(DENC(IA1+I,JB1+J) + F1*DENC(IA2+I,JB2+J))
     &            + RR(M,13)*(DENC(IA1+I,JB1+J) + F1*DENC(IA2+I,JB2+J))
C
                FOCK(IC2+K,JD2+L) = FOCK(IC2+K,JD2+L)
     &            + RR(M,4 )*(DENC(IA1+I,JB1+J) + F1*DENC(IA2+I,JB2+J))
     &            + RR(M,16)*(DENC(IA1+I,JB1+J) + F1*DENC(IA2+I,JB2+J))
              ENDDO
            ENDDO
          ENDIF
C
C         THIRD IFLG BATCH          
          IF(IFLG(3).EQ.1) THEN
            DO L=1,NFUND
              DO K=1,NFUNC
                M1 = (K-1)*NFUND + L
C  
                FOCK(IA1+I,JD1+L) = FOCK(IA1+I,JD1+L)
     &                                    - RR(M1,1 )*DENC(IC1+K,JB1+J) 
     &                                    - RR(M1,7 )*DENC(IC2+K,JB2+J)
C
                FOCK(IA2+I,JD2+L) = FOCK(IA2+I,JD2+L)
     &                                    - RR(M1,10)*DENC(IC1+K,JB1+J) 
     &                                    - RR(M1,16)*DENC(IC2+K,JB2+J)
              ENDDO
            ENDDO
          ENDIF
C
C         FOURTH IFLG BATCH          
          IF(IFLG(4).EQ.1) THEN
            F1 = RKFAC2*FACCD1
            F2 = RKFAC2*FACCD2
            M  = 0
            DO K=1,NFUNC
              DO L=1,NFUND
                M = M+1
C
                FOCK(IA1+I,JC1+K) = FOCK(IA1+I,JC1+K)
     &                                  - F1*RR(M,4 )*DENC(ID1+L,JB1+J)
     &                                  - F2*RR(M,7 )*DENC(ID2+L,JB2+J)
C
                FOCK(IA2+I,JC2+K) = FOCK(IA2+I,JC2+K)
     &                                  - F2*RR(M,10)*DENC(ID1+L,JB1+J)
     &                                  - F1*RR(M,13)*DENC(ID2+L,JB2+J)
              ENDDO
            ENDDO
          ENDIF
C
C         FIFTH IFLG BATCH
          IF(IFLG(5).EQ.1) THEN
            F1 = RKFAC1*FACAB1
            F2 = RKFAC1*FACAB2
            M  = 0
            DO K=1,NFUNC
              DO L=1,NFUND
                M = M+1
                
                FOCK(IC1+K,JA1+I) = FOCK(IC1+K,JA1+I)
     &                                  - F1*RR(M,13)*DENC(ID1+L,JB1+J)
     &                                  - F2*RR(M,10)*DENC(ID2+L,JB2+J)
C     
                FOCK(IC2+K,JA2+I) = FOCK(IC2+K,JA2+I)
     &                                  - F2*RR(M,7 )*DENC(ID1+L,JB1+J)
     &                                  - F1*RR(M,4 )*DENC(ID2+L,JB2+J)
              ENDDO
            ENDDO
          ENDIF
C
C         SIXTH IFLG BATCH
          IF(IFLG(6).EQ.1) THEN
            F1 = RKFAC1*FACAB1
            F2 = RKFAC1*FACAB2
            DO L=1,NFUND
              DO K=1,NFUNC
                M = (K-1)*NFUND + L
C
                FOCK(IB1+J,JD1+L) = FOCK(IB1+J,JD1+L)
     &                                  - F1*RR(M,13)*DENC(IC1+K,JA1+I)
     &                                  - F2*RR(M,7 )*DENC(IC2+K,JA2+I)
C
                FOCK(IB2+J,JD2+L) = FOCK(IB2+J,JD2+L)
     &                                  - F2*RR(M,10)*DENC(IC1+K,JA1+I)
     &                                  - F1*RR(M,4 )*DENC(IC2+K,JA2+I)
              ENDDO
            ENDDO
          ENDIF
C
C         SEVENTH IFLG BATCH
          IF(IFLG(7).EQ.1) THEN
            F1 = RKFAC2*FACCD1
            F2 = RKFAC2*FACCD2
            DO L=1,NFUND
              DO K=1,NFUNC
                M = (K-1)*NFUND + L
C
                FOCK(ID1+L,JB1+J) = FOCK(ID1+L,JB1+J)
     &                                  - F1*RR(M,4 )*DENC(IC1+K,JA1+I)
     &                                  - F2*RR(M,10)*DENC(IC2+K,JA2+I)
C
                FOCK(ID2+L,JB2+J) = FOCK(ID2+L,JB2+J)
     &                                  - F2*RR(M,7 )*DENC(IC1+K,JA1+I)
     &                                  - F1*RR(M,13)*DENC(IC2+K,JA2+I)
              ENDDO
            ENDDO
          ENDIF
C
C         EIGHTH IFLG BATCH
          IF(IFLG(8).EQ.1) THEN
            F11 = RKFAC1*RKFAC2*FACAB1*FACCD1
            F12 = RKFAC1*RKFAC2*FACAB1*FACCD2
            F21 = RKFAC1*RKFAC2*FACAB2*FACCD1
            F22 = RKFAC1*RKFAC2*FACAB2*FACCD2
            M   = 0
            DO K=1,NFUNC
              DO L=1,NFUND
                M = M+1
                
                FOCK(IB1+J,JC1+K) = FOCK(IB1+J,JC1+K)
     &                                 - F11*RR(M,16)*DENC(ID1+L,JA1+I)
     &                                 - F22*RR(M,7 )*DENC(ID2+L,JA2+I)
C
                  FOCK(IB2+J,JC2+K) = FOCK(IB2+J,JC2+K)
     &                                 - F22*RR(M,10)*DENC(ID1+L,JA1+I)
     &                                 - F11*RR(M,1 )*DENC(ID2+L,JA2+I)
              ENDDO
            ENDDO
          ENDIF
C
C         NINTH IFLG BATCH
          IF(IFLG(9).EQ.1) THEN
            M = 0
            DO K=1,NFUNC
              DO L=1,NFUND
                M = M+1
C
                FOCK(IC1+K,JB1+J) = FOCK(IC1+K,JB1+J)
     &                                     - RR(M,1 )*DENC(ID1+L,JA1+I)
     &                                     - RR(M,10)*DENC(ID2+L,JA2+I)
C
                FOCK(IC2+K,JB2+J) = FOCK(IC2+K,JB2+J)
     &                                     - RR(M,7 )*DENC(ID1+L,JA1+I)
     &                                     - RR(M,16)*DENC(ID2+L,JA2+I)
              ENDDO
            ENDDO
          ENDIF
C
C         TENTH IFLG BATCH
          IF(IFLG(10).EQ.1) THEN
            M=0
            DO K=1,NFUNC
              DO L=1,NFUND
                M=M+1
C
                FOCK(IA1+I,JB1+J) = FOCK(IA1+I,JB1+J)
     &                                     + RR(M,1 )*DENC(IC1+K,JD1+L)
     &                                     + RR(M,4 )*DENC(IC2+K,JD2+L)
C
                FOCK(IA2+I,JB2+J) = FOCK(IA2+I,JB2+J)
     &                                     + RR(M,13)*DENC(IC1+K,JD1+L)
     &                                     + RR(M,16)*DENC(IC2+K,JD2+L)
              ENDDO
            ENDDO
          ENDIF
C
C         ELEVENTH IFLG BATCH
          IF(IFLG(11).EQ.1) THEN
            M = 0
            DO K=1,NFUNC
              DO L=1,NFUND
                M = M+1
C
                FOCK(IC1+K,JD1+L) = FOCK(IC1+K,JD1+L)
     &                                     + RR(M,1 )*DENC(IA1+I,JB1+J)
     &                                     + RR(M,13)*DENC(IA2+I,JB2+J)
C
                FOCK(IC2+K,JD2+L) = FOCK(IC2+K,JD2+L)
     &                                     + RR(M,4 )*DENC(IA1+I,JB1+J)
     &                                     + RR(M,16)*DENC(IA2+I,JB2+J)
              ENDDO
            ENDDO
          ENDIF
C          
        ELSE
C       GENERAL CASE (NO GEOMETRIC SYMMETRY):
C         FIRST IFLG BATCH
          IF(IFLG(1).EQ.1) THEN
            F1 = RKFAC2*FACCD1
            F2 = RKFAC2*FACCD2
            M  = 0
            DO K=1,NFUNC
              DO L=1,NFUND
                M = M+1
C
                FOCK(IA1+I,JB1+J) = FOCK(IA1+I,JB1+J)
     &            + RR(M,1 )*(DENC(IC1+K,JD1+L) + F1*DENC(JD2+L,IC2+K))
     &            + RR(M,2 )*(DENC(IC1+K,JD2+L) + F2*DENC(JD1+L,IC2+K))
     &            + RR(M,3 )*(DENC(IC2+K,JD1+L) + F2*DENC(JD2+L,IC1+K))
     &            + RR(M,4 )*(DENC(IC2+K,JD2+L) + F1*DENC(JD1+L,IC1+K))
C
                FOCK(IA1+I,JB2+J) = FOCK(IA1+I,JB2+J)
     &            + RR(M,5 )*(DENC(IC1+K,JD1+L) + F1*DENC(JD2+L,IC2+K))
     &            + RR(M,6 )*(DENC(IC1+K,JD2+L) + F2*DENC(JD1+L,IC2+K))
     &            + RR(M,7 )*(DENC(IC2+K,JD1+L) + F2*DENC(JD2+L,IC1+K))
     &            + RR(M,8 )*(DENC(IC2+K,JD2+L) + F1*DENC(JD1+L,IC1+K))
C
                FOCK(IA2+I,JB1+J) = FOCK(IA2+I,JB1+J)
     &            + RR(M,9 )*(DENC(IC1+K,JD1+L) + F1*DENC(JD2+L,IC2+K))
     &            + RR(M,10)*(DENC(IC1+K,JD2+L) + F2*DENC(JD1+L,IC2+K))
     &            + RR(M,11)*(DENC(IC2+K,JD1+L) + F2*DENC(JD2+L,IC1+K))
     &            + RR(M,12)*(DENC(IC2+K,JD2+L) + F1*DENC(JD1+L,IC1+K))
C
                FOCK(IA2+I,JB2+J) = FOCK(IA2+I,JB2+J)
     &            + RR(M,13)*(DENC(IC1+K,JD1+L) + F1*DENC(JD2+L,IC2+K))
     &            + RR(M,14)*(DENC(IC1+K,JD2+L) + F2*DENC(JD1+L,IC2+K))
     &            + RR(M,15)*(DENC(IC2+K,JD1+L) + F2*DENC(JD2+L,IC1+K))
     &            + RR(M,16)*(DENC(IC2+K,JD2+L) + F1*DENC(JD1+L,IC1+K))

c      IF(IA1+I+JB1+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',1,FOCK(IA1+I,JB1+J)
c      ENDIF
c      IF(IA1+I+JB2+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',1,FOCK(IA1+I,JB2+J)
c      ENDIF
c      IF(IA2+I+JB1+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',1,FOCK(IA2+I,JB1+J)
c      ENDIF
c      IF(IA2+I+JB2+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',1,FOCK(IA2+I,JB2+J)
c      ENDIF  
              ENDDO
            ENDDO 
          ENDIF
C
C         SECOND IFLG BATCH
          IF(IFLG(2).EQ.1) THEN
            F1 = RKFAC1*FACAB1
            F2 = RKFAC1*FACAB2
            M  = 0
            DO K=1,NFUNC
              DO L=1,NFUND
                M = M+1
C
                FOCK(IC1+K,JD1+L) = FOCK(IC1+K,JD1+L)
     &            + RR(M,1 )*(DENC(IA1+I,JB1+J) + F1*DENC(JB2+J,IA2+I))
     &            + RR(M,5 )*(DENC(IA1+I,JB2+J) + F2*DENC(JB1+J,IA2+I))
     &            + RR(M,9 )*(DENC(IA2+I,JB1+J) + F2*DENC(JB2+J,IA1+I))
     &            + RR(M,13)*(DENC(IA2+I,JB2+J) + F1*DENC(JB1+J,IA1+I))
C
                FOCK(IC1+K,JD2+L) = FOCK(IC1+K,JD2+L)
     &            + RR(M,2 )*(DENC(IA1+I,JB1+J) + F1*DENC(JB2+J,IA2+I))
     &            + RR(M,6 )*(DENC(IA1+I,JB2+J) + F2*DENC(JB1+J,IA2+I))
     &            + RR(M,10)*(DENC(IA2+I,JB1+J) + F2*DENC(JB2+J,IA1+I))
     &            + RR(M,14)*(DENC(IA2+I,JB2+J) + F1*DENC(JB1+J,IA1+I)) 
C
                FOCK(IC2+K,JD1+L) = FOCK(IC2+K,JD1+L)
     &            + RR(M,3 )*(DENC(IA1+I,JB1+J) + F1*DENC(JB2+J,IA2+I))
     &            + RR(M,7 )*(DENC(IA1+I,JB2+J) + F2*DENC(JB1+J,IA2+I))
     &            + RR(M,11)*(DENC(IA2+I,JB1+J) + F2*DENC(JB2+J,IA1+I))
     &            + RR(M,15)*(DENC(IA2+I,JB2+J) + F1*DENC(JB1+J,IA1+I))
C
                FOCK(IC2+K,JD2+L) = FOCK(IC2+K,JD2+L)
     &            + RR(M,4 )*(DENC(IA1+I,JB1+J) + F1*DENC(JB2+J,IA2+I))
     &            + RR(M,8 )*(DENC(IA1+I,JB2+J) + F2*DENC(JB1+J,IA2+I))
     &            + RR(M,12)*(DENC(IA2+I,JB1+J) + F2*DENC(JB2+J,IA1+I))
     &            + RR(M,16)*(DENC(IA2+I,JB2+J) + F1*DENC(JB1+J,IA1+I))

c      IF(IC1+K+JD1+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',2,FOCK(IC1+K,JD1+L)
c      ENDIF
c      IF(IC1+K+JD2+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',2,FOCK(IC1+K,JD2+L)
c      ENDIF
c      IF(IC2+K+JD1+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',2,FOCK(IC2+K,JD1+L)
c      ENDIF
c      IF(IC2+K+JD2+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',2,FOCK(IC2+K,JD2+L)
c      ENDIF  
              ENDDO
            ENDDO

          ENDIF
C
C         THIRD IFLG BATCH          
          IF(IFLG(3).EQ.1) THEN
            DO L=1,NFUND
              DO K=1,NFUNC
                M1 = (K-1)*NFUND + L
C  
                FOCK(IA1+I,JD1+L) = FOCK(IA1+I,JD1+L)
     &                                    - RR(M1,1 )*DENC(IC1+K,JB1+J) 
     &                                    - RR(M1,3 )*DENC(IC2+K,JB1+J)
     &                                    - RR(M1,5 )*DENC(IC1+K,JB2+J) 
     &                                    - RR(M1,7 )*DENC(IC2+K,JB2+J)
C
                FOCK(IA1+I,JD2+L) = FOCK(IA1+I,JD2+L)
     &                                    - RR(M1,2 )*DENC(IC1+K,JB1+J) 
     &                                    - RR(M1,4 )*DENC(IC2+K,JB1+J)
     &                                    - RR(M1,6 )*DENC(IC1+K,JB2+J) 
     &                                    - RR(M1,8 )*DENC(IC2+K,JB2+J)
C
                FOCK(IA2+I,JD1+L) = FOCK(IA2+I,JD1+L)
     &                                    - RR(M1,9 )*DENC(IC1+K,JB1+J) 
     &                                    - RR(M1,11)*DENC(IC2+K,JB1+J)
     &                                    - RR(M1,13)*DENC(IC1+K,JB2+J) 
     &                                    - RR(M1,15)*DENC(IC2+K,JB2+J)
C
                FOCK(IA2+I,JD2+L) = FOCK(IA2+I,JD2+L)
     &                                    - RR(M1,10)*DENC(IC1+K,JB1+J) 
     &                                    - RR(M1,12)*DENC(IC2+K,JB1+J)
     &                                    - RR(M1,14)*DENC(IC1+K,JB2+J) 
     &                                    - RR(M1,16)*DENC(IC2+K,JB2+J)

c      IF(IA1+I+JD1+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',3,FOCK(IA1+I,JD1+L)
c      ENDIF
c      IF(IA1+I+JD2+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',3,FOCK(IA1+I,JD2+L)
c      ENDIF
c      IF(IA2+I+JD1+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',3,FOCK(IA2+I,JD1+L)
c      ENDIF
c      IF(IA2+I+JD2+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',3,FOCK(IA2+I,JD2+L)
c      ENDIF  

              ENDDO
            ENDDO
          ENDIF
C
C         FOURTH IFLG BATCH          
          IF(IFLG(4).EQ.1) THEN
            F1 = RKFAC2*FACCD1
            F2 = RKFAC2*FACCD2
            M  = 0
            DO K=1,NFUNC
              DO L=1,NFUND
                M = M+1
C
                FOCK(IA1+I,JC1+K) = FOCK(IA1+I,JC1+K)
     &                                  - F1*RR(M,4 )*DENC(ID1+L,JB1+J)
     &                                  - F1*RR(M,8 )*DENC(ID1+L,JB2+J)
     &                                  - F2*RR(M,3 )*DENC(ID2+L,JB1+J)
     &                                  - F2*RR(M,7 )*DENC(ID2+L,JB2+J)
     
                FOCK(IA1+I,JC2+K) = FOCK(IA1+I,JC2+K)
     &                                  - F2*RR(M,2 )*DENC(ID1+L,JB1+J)
     &                                  - F2*RR(M,6 )*DENC(ID1+L,JB2+J)
     &                                  - F1*RR(M,1 )*DENC(ID2+L,JB1+J)
     &                                  - F1*RR(M,5 )*DENC(ID2+L,JB2+J)
     
                FOCK(IA2+I,JC1+K) = FOCK(IA2+I,JC1+K)
     &                                  - F1*RR(M,12)*DENC(ID1+L,JB1+J)
     &                                  - F1*RR(M,16)*DENC(ID1+L,JB2+J)
     &                                  - F2*RR(M,11)*DENC(ID2+L,JB1+J)
     &                                  - F2*RR(M,15)*DENC(ID2+L,JB2+J)
     
                FOCK(IA2+I,JC2+K) = FOCK(IA2+I,JC2+K)
     &                                  - F2*RR(M,10)*DENC(ID1+L,JB1+J)
     &                                  - F2*RR(M,14)*DENC(ID1+L,JB2+J)
     &                                  - F1*RR(M,9 )*DENC(ID2+L,JB1+J)
     &                                  - F1*RR(M,13)*DENC(ID2+L,JB2+J)

c      IF(IA1+I+JC1+K.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',4,FOCK(IA1+I,JC1+K)
c      ENDIF
c      IF(IA1+I+JC2+K.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',4,FOCK(IA1+I,JC2+K)
c      ENDIF
c      IF(IA2+I+JC1+K.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',4,FOCK(IA2+I,JC1+K)
c      ENDIF
c      IF(IA2+I+JC2+K.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',4,FOCK(IA2+I,JC2+K)
c      ENDIF 

              ENDDO
            ENDDO

          ENDIF
C
C         FIFTH IFLG BATCH
          IF(IFLG(5).EQ.1) THEN
            F1 = RKFAC1*FACAB1
            F2 = RKFAC1*FACAB2
            M  = 0
            DO K=1,NFUNC
              DO L=1,NFUND
                M = M+1
                
                FOCK(IC1+K,JA1+I) = FOCK(IC1+K,JA1+I)
     &                                  - F1*RR(M,13)*DENC(JB1+J,ID1+L)
     &                                  - F1*RR(M,14)*DENC(JB1+J,ID2+L)
     &                                  - F2*RR(M,9 )*DENC(JB2+J,ID1+L)
     &                                  - F2*RR(M,10)*DENC(JB2+J,ID2+L)
     
                FOCK(IC1+K,JA2+I) = FOCK(IC1+K,JA2+I)
     &                                  - F2*RR(M,5 )*DENC(JB1+J,ID1+L)
     &                                  - F2*RR(M,6 )*DENC(JB1+J,ID2+L)
     &                                  - F1*RR(M,1 )*DENC(JB2+J,ID1+L)
     &                                  - F1*RR(M,2 )*DENC(JB2+J,ID2+L)
     
                FOCK(IC2+K,JA1+I) = FOCK(IC2+K,JA1+I)
     &                                  - F1*RR(M,15)*DENC(JB1+J,ID1+L)
     &                                  - F1*RR(M,16)*DENC(JB1+J,ID2+L)
     &                                  - F2*RR(M,11)*DENC(JB2+J,ID1+L)
     &                                  - F2*RR(M,12)*DENC(JB2+J,ID2+L)
     
                FOCK(IC2+K,JA2+I) = FOCK(IC2+K,JA2+I)
     &                                  - F2*RR(M,7 )*DENC(JB1+J,ID1+L)
     &                                  - F2*RR(M,8 )*DENC(JB1+J,ID2+L)
     &                                  - F1*RR(M,3 )*DENC(JB2+J,ID1+L)
     &                                  - F1*RR(M,4 )*DENC(JB2+J,ID2+L)

c      IF(IC1+K+JA1+I.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',5,FOCK(IC1+K,JA1+I)
c      ENDIF
c      IF(IC1+K+JA2+I.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',5,FOCK(IC1+K,JA2+I)
c      ENDIF
c      IF(IC2+K+JA1+I.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',5,FOCK(IC2+K,JA1+I)
c      ENDIF
c      IF(IC2+K+JA2+I.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',5,FOCK(IC2+K,JA2+I)
c      ENDIF  

              ENDDO
            ENDDO

          ENDIF
C
C         SIXTH IFLG BATCH
          IF(IFLG(6).EQ.1) THEN
            F1 = RKFAC1*FACAB1
            F2 = RKFAC1*FACAB2
            DO L=1,NFUND
              DO K=1,NFUNC
                M = (K-1)*NFUND + L
C
                FOCK(IB1+J,JD1+L) = FOCK(IB1+J,JD1+L)
     &                                  - F1*RR(M,13)*DENC(IC1+K,JA1+I)
     &                                  - F2*RR(M,5 )*DENC(IC1+K,JA2+I)
     &                                  - F1*RR(M,15)*DENC(IC2+K,JA1+I)
     &                                  - F2*RR(M,7 )*DENC(IC2+K,JA2+I)
C
                FOCK(IB1+J,JD2+L) = FOCK(IB1+J,JD2+L)
     &                                  - F1*RR(M,14)*DENC(IC1+K,JA1+I)
     &                                  - F2*RR(M,6 )*DENC(IC1+K,JA2+I)
     &                                  - F1*RR(M,16)*DENC(IC2+K,JA1+I)
     &                                  - F2*RR(M,8 )*DENC(IC2+K,JA2+I)
C
                FOCK(IB2+J,JD1+L) = FOCK(IB2+J,JD1+L)
     &                                  - F2*RR(M,9 )*DENC(IC1+K,JA1+I)
     &                                  - F1*RR(M,1 )*DENC(IC1+K,JA2+I)
     &                                  - F2*RR(M,11)*DENC(IC2+K,JA1+I)
     &                                  - F1*RR(M,3 )*DENC(IC2+K,JA2+I)
C
                FOCK(IB2+J,JD2+L) = FOCK(IB2+J,JD2+L)
     &                                  - F2*RR(M,10)*DENC(IC1+K,JA1+I)
     &                                  - F1*RR(M,2 )*DENC(IC1+K,JA2+I)
     &                                  - F2*RR(M,12)*DENC(IC2+K,JA1+I)
     &                                  - F1*RR(M,4 )*DENC(IC2+K,JA2+I)

c      IF(IB1+J+JD1+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',6,FOCK(IB1+J,JD1+L)
c      ENDIF
c      IF(IB1+J+JD2+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',6,FOCK(IB1+J,JD2+L)
c      ENDIF
c      IF(IB2+J+JD1+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',6,FOCK(IB2+J,JD1+L)
c      ENDIF
c      IF(IB2+J+JD2+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',6,FOCK(IB2+J,JD2+L)
c      ENDIF 

              ENDDO
            ENDDO

          ENDIF
C
C         SEVENTH IFLG BATCH
          IF(IFLG(7).EQ.1) THEN
            F1 = RKFAC2*FACCD1
            F2 = RKFAC2*FACCD2
            DO L=1,NFUND
              DO K=1,NFUNC
                M = (K-1)*NFUND + L
C
                FOCK(ID1+L,JB1+J) = FOCK(ID1+L,JB1+J)
     &                                  - F1*RR(M,4 )*DENC(JA1+I,IC1+K)
     &                                  - F2*RR(M,2 )*DENC(JA1+I,IC2+K)
     &                                  - F1*RR(M,12)*DENC(JA2+I,IC1+K)
     &                                  - F2*RR(M,10)*DENC(JA2+I,IC2+K)
C
                FOCK(ID1+L,JB2+J) = FOCK(ID1+L,JB2+J)
     &                                  - F1*RR(M,8 )*DENC(JA1+I,IC1+K)
     &                                  - F2*RR(M,6 )*DENC(JA1+I,IC2+K)
     &                                  - F1*RR(M,16)*DENC(JA2+I,IC1+K)
     &                                  - F2*RR(M,14)*DENC(JA2+I,IC2+K)
C
                FOCK(ID2+L,JB1+J) = FOCK(ID2+L,JB1+J)
     &                                  - F2*RR(M,3 )*DENC(JA1+I,IC1+K)
     &                                  - F1*RR(M,1 )*DENC(JA1+I,IC2+K)
     &                                  - F2*RR(M,11)*DENC(JA2+I,IC1+K)
     &                                  - F1*RR(M,9 )*DENC(JA2+I,IC2+K)
C
                FOCK(ID2+L,JB2+J) = FOCK(ID2+L,JB2+J)
     &                                  - F2*RR(M,7 )*DENC(JA1+I,IC1+K)
     &                                  - F1*RR(M,5 )*DENC(JA1+I,IC2+K)
     &                                  - F2*RR(M,15)*DENC(JA2+I,IC1+K)
     &                                  - F1*RR(M,13)*DENC(JA2+I,IC2+K)

c      IF(ID1+L+JB1+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',7,FOCK(ID1+L,JB1+J)
c      ENDIF
c      IF(ID1+L+JB2+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',7,FOCK(ID1+L,JB2+J)
c      ENDIF
c      IF(ID2+L+JB1+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',7,FOCK(ID2+L,JB1+J)
c      ENDIF
c      IF(ID2+L+JB2+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',7,FOCK(ID2+L,JB2+J)
c      ENDIF 

              ENDDO
            ENDDO

          ENDIF
C
C         EIGHTH IFLG BATCH
          IF(IFLG(8).EQ.1) THEN
            F11 = RKFAC1*RKFAC2*FACAB1*FACCD1
            F12 = RKFAC1*RKFAC2*FACAB1*FACCD2
            F21 = RKFAC1*RKFAC2*FACAB2*FACCD1
            F22 = RKFAC1*RKFAC2*FACAB2*FACCD2
            M   = 0
            DO K=1,NFUNC
              DO L=1,NFUND
                M = M+1
                
                FOCK(IB1+J,JC1+K) = FOCK(IB1+J,JC1+K)
     &                                 - F11*RR(M,16)*DENC(ID1+L,JA1+I)
     &                                 - F21*RR(M,8 )*DENC(ID1+L,JA2+I)
     &                                 - F12*RR(M,15)*DENC(ID2+L,JA1+I)
     &                                 - F22*RR(M,7 )*DENC(ID2+L,JA2+I)
C
                FOCK(IB1+J,JC2+K) = FOCK(IB1+J,JC2+K)
     &                                 - F12*RR(M,14)*DENC(ID1+L,JA1+I)
     &                                 - F22*RR(M,6 )*DENC(ID1+L,JA2+I)
     &                                 - F11*RR(M,13)*DENC(ID2+L,JA1+I)
     &                                 - F21*RR(M,5 )*DENC(ID2+L,JA2+I)
C
                FOCK(IB2+J,JC1+K) = FOCK(IB2+J,JC1+K)
     &                                 - F21*RR(M,12)*DENC(ID1+L,JA1+I)
     &                                 - F11*RR(M,4 )*DENC(ID1+L,JA2+I)
     &                                 - F22*RR(M,11)*DENC(ID2+L,JA1+I)
     &                                 - F12*RR(M,3 )*DENC(ID2+L,JA2+I)
C
                FOCK(IB2+J,JC2+K) = FOCK(IB2+J,JC2+K)
     &                                 - F22*RR(M,10)*DENC(ID1+L,JA1+I)
     &                                 - F12*RR(M,2 )*DENC(ID1+L,JA2+I)
     &                                 - F21*RR(M,9 )*DENC(ID2+L,JA1+I)
     &                                 - F11*RR(M,1 )*DENC(ID2+L,JA2+I)

c      IF(IB1+J+JC1+K.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',8,FOCK(IB1+J,JC1+K)
c      ENDIF
c      IF(IB1+J+JC2+K.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',8,FOCK(IB1+J,JC2+K)
c      ENDIF
c      IF(IB2+J+JC1+K.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',8,FOCK(IB2+J,JC1+K)
c      ENDIF
c      IF(IB2+J+JC2+K.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',8,FOCK(IB2+J,JC2+K)
c      ENDIF
      
              ENDDO
            ENDDO

          ENDIF
C
C         NINTH IFLG BATCH
          IF(IFLG(9).EQ.1) THEN
            M = 0
            DO K=1,NFUNC
              DO L=1,NFUND
                M = M+1
C
                FOCK(IC1+K,JB1+J) = FOCK(IC1+K,JB1+J)
     &                                     - RR(M,1 )*DENC(JA1+I,ID1+L)
     &                                     - RR(M,2 )*DENC(JA1+I,ID2+L)
     &                                     - RR(M,9 )*DENC(JA2+I,ID1+L)
     &                                     - RR(M,10)*DENC(JA2+I,ID2+L)
C
                FOCK(IC1+K,JB2+J) = FOCK(IC1+K,JB2+J)
     &                                     - RR(M,5 )*DENC(JA1+I,ID1+L)
     &                                     - RR(M,6 )*DENC(JA1+I,ID2+L)
     &                                     - RR(M,13)*DENC(JA2+I,ID1+L)
     &                                     - RR(M,14)*DENC(JA2+I,ID2+L)
C
                FOCK(IC2+K,JB1+J) = FOCK(IC2+K,JB1+J)
     &                                     - RR(M,3 )*DENC(JA1+I,ID1+L)
     &                                     - RR(M,4 )*DENC(JA1+I,ID2+L)
     &                                     - RR(M,11)*DENC(JA2+I,ID1+L)
     &                                     - RR(M,12)*DENC(JA2+I,ID2+L)
C
                FOCK(IC2+K,JB2+J) = FOCK(IC2+K,JB2+J)
     &                                     - RR(M,7 )*DENC(JA1+I,ID1+L)
     &                                     - RR(M,8 )*DENC(JA1+I,ID2+L)
     &                                     - RR(M,15)*DENC(JA2+I,ID1+L)
     &                                     - RR(M,16)*DENC(JA2+I,ID2+L)

c      IF(IC1+K+JB1+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',9,FOCK(IC1+K,JB1+J)
c      ENDIF
c      IF(IC1+K+JB2+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',9,FOCK(IC1+K,JB2+J)
c      ENDIF
c      IF(IC2+K+JB1+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',9,FOCK(IC2+K,JB1+J)
c      ENDIF
c      IF(IC2+K+JB2+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',9,FOCK(IC2+K,JB2+J)
c      ENDIF

              ENDDO
            ENDDO

          ENDIF
C
C         TENTH IFLG BATCH
          IF(IFLG(10).EQ.1) THEN
            M=0
            DO K=1,NFUNC
              DO L=1,NFUND
                M=M+1
C
                FOCK(IA1+I,JB1+J) = FOCK(IA1+I,JB1+J)
     &                                     + RR(M,1 )*DENC(IC1+K,JD1+L)
     &                                     + RR(M,2 )*DENC(IC1+K,JD2+L)
     &                                     + RR(M,3 )*DENC(IC2+K,JD1+L)
     &                                     + RR(M,4 )*DENC(IC2+K,JD2+L)
C
                FOCK(IA1+I,JB2+J) = FOCK(IA1+I,JB2+J)
     &                                     + RR(M,5 )*DENC(IC1+K,JD1+L)
     &                                     + RR(M,6 )*DENC(IC1+K,JD2+L)
     &                                     + RR(M,7 )*DENC(IC2+K,JD1+L)
     &                                     + RR(M,8 )*DENC(IC2+K,JD2+L)
C
                FOCK(IA2+I,JB1+J) = FOCK(IA2+I,JB1+J)
     &                                     + RR(M,9 )*DENC(IC1+K,JD1+L)
     &                                     + RR(M,10)*DENC(IC1+K,JD2+L)
     &                                     + RR(M,11)*DENC(IC2+K,JD1+L)
     &                                     + RR(M,12)*DENC(IC2+K,JD2+L)
C
                FOCK(IA2+I,JB2+J) = FOCK(IA2+I,JB2+J)
     &                                     + RR(M,13)*DENC(IC1+K,JD1+L)
     &                                     + RR(M,14)*DENC(IC1+K,JD2+L)
     &                                     + RR(M,15)*DENC(IC2+K,JD1+L)
     &                                     + RR(M,16)*DENC(IC2+K,JD2+L)

c      IF(IA1+I+JB1+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',10,FOCK(IA1+I,JB1+J)
c      ENDIF
c      IF(IA1+I+JB2+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',10,FOCK(IA1+I,JB2+J)
c      ENDIF
c      IF(IA2+I+JB1+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',10,FOCK(IA2+I,JB1+J)
c      ENDIF
c      IF(IA2+I+JB2+J.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',10,FOCK(IA2+I,JB2+J)
c      ENDIF

              ENDDO
            ENDDO


          ENDIF
C
C         ELEVENTH IFLG BATCH
          IF(IFLG(11).EQ.1) THEN
            M = 0
            DO K=1,NFUNC
              DO L=1,NFUND
                M = M+1
C
                FOCK(IC1+K,JD1+L) = FOCK(IC1+K,JD1+L)
     &                                     + RR(M,1 )*DENC(IA1+I,JB1+J)
     &                                     + RR(M,5 )*DENC(IA1+I,JB2+J)
     &                                     + RR(M,9 )*DENC(IA2+I,JB1+J)
     &                                     + RR(M,13)*DENC(IA2+I,JB2+J)
C
                FOCK(IC1+K,JD2+L) = FOCK(IC1+K,JD2+L)
     &                                     + RR(M,2 )*DENC(IA1+I,JB1+J)
     &                                     + RR(M,6 )*DENC(IA1+I,JB2+J)
     &                                     + RR(M,10)*DENC(IA2+I,JB1+J)
     &                                     + RR(M,14)*DENC(IA2+I,JB2+J)

                FOCK(IC2+K,JD1+L) = FOCK(IC2+K,JD1+L)
     &                                     + RR(M,3 )*DENC(IA1+I,JB1+J)
     &                                     + RR(M,7 )*DENC(IA1+I,JB2+J)
     &                                     + RR(M,11)*DENC(IA2+I,JB1+J)
     &                                     + RR(M,15)*DENC(IA2+I,JB2+J)
C
                FOCK(IC2+K,JD2+L) = FOCK(IC2+K,JD2+L)
     &                                     + RR(M,4 )*DENC(IA1+I,JB1+J)
     &                                     + RR(M,8 )*DENC(IA1+I,JB2+J)
     &                                     + RR(M,12)*DENC(IA2+I,JB1+J)
     &                                     + RR(M,16)*DENC(IA2+I,JB2+J)

c      IF(IC1+K+JD1+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',11,FOCK(IC1+K,JD1+L)
c      ENDIF
c      IF(IC1+K+JD2+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',11,FOCK(IC1+K,JD2+L)
c      ENDIF
c      IF(IC2+K+JD1+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',11,FOCK(IC2+K,JD1+L)
c      ENDIF
c      IF(IC2+K+JD2+L.EQ.2) THEN
c        WRITE(*,*) 'fockupdate',11,FOCK(IC2+K,JD2+L)
c      ENDIF

              ENDDO
            ENDDO
          ENDIF
        ENDIF
C
C*********************************************************************C
C     CONSTRUCTION OF THE OPEN-SHELL PART OF THE QMAT MATRIX          C
C*********************************************************************C
        IF(NOPEN.NE.0) THEN
C
          IF(NCENT.LE.2) THEN
C         DIATOMIC CASE:
C
C           FIRST IFLG BATCH
            IF(IFLG(1).EQ.1) THEN
              F1 = RKFAC2*FACCD1
              F2 = RKFAC2*FACCD2
              M  = 0
              DO K=1,NFUNC
                DO L=1,NFUND
                  M = M+1
C
                  QMAT(IA1+I,JB1+J) = QMAT(IA1+I,JB1+J)
     &                + ACFF*(RR(M,1 ) + F1*RR(M,4 ))*DENO(IC1+K,JD1+L)
     &                + ACFF*(RR(M,4 ) + F1*RR(M,1 ))*DENO(IC2+K,JD2+L)
C
                  QMAT(IA2+I,JB2+J) = QMAT(IA2+I,JB2+J)
     &                + ACFF*(RR(M,13) + F1*RR(M,16))*DENO(IC1+K,JD1+L)
     &                + ACFF*(RR(M,16) + F1*RR(M,13))*DENO(IC2+K,JD2+L)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND IFLG BATCH
            IF(IFLG(2).EQ.1) THEN
              F1 = RKFAC1*FACAB1
              F2 = RKFAC1*FACAB2
              M  = 0
              DO K=1,NFUNC
                DO L=1,NFUND
                  M = M+1
C
                  QMAT(IC1+K,JD1+L) = QMAT(IC1+K,JD1+L)
     &                + ACFF*(RR(M,1 ) + F1*RR(M,13))*DENO(IA1+I,JB1+J)
     &                + ACFF*(RR(M,13) + F1*RR(M,1 ))*DENO(IA2+I,JB2+J)
C
                  QMAT(IC2+K,JD2+L) = QMAT(IC2+K,JD2+L)
     &                + ACFF*(RR(M,4 ) + F1*RR(M,16))*DENO(IA1+I,JB1+J)
     &                + ACFF*(RR(M,16) + F1*RR(M,4 ))*DENO(IA2+I,JB2+J)
                ENDDO
              ENDDO
            ENDIF
C
C           THIRD IFLG BATCH
            IF(IFLG(3).EQ.1) THEN
              DO L=1,NFUND
                DO K=1,NFUNC
                  M1 = (K-1)*NFUND + L
C
                  QMAT(IA1+I,JD1+L) = QMAT(IA1+I,JD1+L)
     &                               - BCFF*RR(M1,1 )*DENO(IC1+K,JB1+J)
     &                               - BCFF*RR(M1,7 )*DENO(IC2+K,JB2+J)
C
                  QMAT(IA2+I,JD2+L) = QMAT(IA2+I,JD2+L)
     &                               - BCFF*RR(M1,10)*DENO(IC1+K,JB1+J)
     &                               - BCFF*RR(M1,16)*DENO(IC2+K,JB2+J)
                ENDDO
              ENDDO
            ENDIF
C
C           FOURTH IFLG BATCH
            IF(IFLG(4).EQ.1) THEN
              F1 = RKFAC2*FACCD1
              F2 = RKFAC2*FACCD2
              M  = 0
              DO K=1,NFUNC
                DO L=1,NFUND
                  M = M+1
C
                  QMAT(IA1+I,JC1+K) = QMAT(IA1+I,JC1+K)
     &                             - F1*BCFF*RR(M,4 )*DENO(ID1+L,JB1+J)
     &                             - F2*BCFF*RR(M,7 )*DENO(ID2+L,JB2+J)
C
                  QMAT(IA2+I,JC2+K) = QMAT(IA2+I,JC2+K)
     &                             - F2*BCFF*RR(M,10)*DENO(ID1+L,JB1+J)
     &                             - F1*BCFF*RR(M,13)*DENO(ID2+L,JB2+J)
                ENDDO
              ENDDO
            ENDIF
C
C           FIFTH IFLG BATCH
            IF(IFLG(5).EQ.1) THEN
              F1 = RKFAC1*FACAB1
              F2 = RKFAC1*FACAB2
              M  = 0
              DO K=1,NFUNC
                DO L=1,NFUND
                  M = M+1
C
                  QMAT(IC1+K,JA1+I) = QMAT(IC1+K,JA1+I)
     &                             - F1*BCFF*RR(M,13)*DENO(ID1+L,JB1+J)
     &                             - F2*BCFF*RR(M,10)*DENO(ID2+L,JB2+J)
C
                 QMAT(IC2+K,JA2+I) = QMAT(IC2+K,JA2+I)
     &                             - F2*BCFF*RR(M,7 )*DENO(ID1+L,JB1+J)
     &                             - F1*BCFF*RR(M,4 )*DENO(ID2+L,JB2+J)
                ENDDO
              ENDDO
            ENDIF
C
C           SIXTH IFLG BATCH
            IF(IFLG(6).EQ.1) THEN
              F1 = RKFAC1*FACAB1
              F2 = RKFAC1*FACAB2
              DO L=1,NFUND
                DO K=1,NFUNC
                  M = (K-1)*NFUND+L
C
                  QMAT(IB1+J,JD1+L) = QMAT(IB1+J,JD1+L)
     &                             - F1*BCFF*RR(M,13)*DENO(IC1+K,JA1+I)
     &                             - F2*BCFF*RR(M,7 )*DENO(IC2+K,JA2+I)
C
                  QMAT(IB2+J,JD2+L) = QMAT(IB2+J,JD2+L)
     &                             - F2*BCFF*RR(M,10)*DENO(IC1+K,JA1+I)
     &                             - F1*BCFF*RR(M,4 )*DENO(IC2+K,JA2+I)
                ENDDO
              ENDDO
            ENDIF
C
C           SEVENTH IFLG BATCH
            IF(IFLG(7).EQ.1) THEN
              F1 = RKFAC2*FACCD1
              F2 = RKFAC2*FACCD2
              DO L=1,NFUND
                DO K=1,NFUNC
                  M = (K-1)*NFUND+L
C
                  QMAT(ID1+L,JB1+J) = QMAT(ID1+L,JB1+J)
     &                             - F1*BCFF*RR(M,4 )*DENO(IC1+K,JA1+I)
     &                             - F2*BCFF*RR(M,10)*DENO(IC2+K,JA2+I)
C
                  QMAT(ID2+L,JB2+J) = QMAT(ID2+L,JB2+J)
     &                             - F2*BCFF*RR(M,7 )*DENO(IC1+K,JA1+I)
     &                             - F1*BCFF*RR(M,13)*DENO(IC2+K,JA2+I)
                ENDDO
              ENDDO
            ENDIF
C
C           EIGHTH IFLG BATCH
            IF(IFLG(8).EQ.1) THEN
              F11 = RKFAC1*RKFAC2*FACAB1*FACCD1
              F22 = RKFAC1*RKFAC2*FACAB2*FACCD2
              M   = 0
              DO K=1,NFUNC
                DO L=1,NFUND
                  M = M+1
C
                  QMAT(IB1+J,JC1+K) = QMAT(IB1+J,JC1+K)
     &                            - F11*BCFF*RR(M,16)*DENO(ID1+L,JA1+I)
     &                            - F22*BCFF*RR(M,7 )*DENO(ID2+L,JA2+I)
C
                  QMAT(IB2+J,JC2+K) = QMAT(IB2+J,JC2+K)
     &                            - F22*BCFF*RR(M,10)*DENO(ID1+L,JA1+I)
     &                            - F11*BCFF*RR(M,1 )*DENO(ID2+L,JA2+I)
                  ENDDO
                ENDDO
              ENDIF
C
C           NINTH IFLG BATCH
            IF(IFLG(9).EQ.1) THEN
              M = 0
              DO K=1,NFUNC
                DO L=1,NFUND
                  M = M+1
C
                  QMAT(IC1+K,JB1+J) = QMAT(IC1+K,JB1+J)
     &                                - BCFF*RR(M,1 )*DENO(ID1+L,JA1+I)
     &                                - BCFF*RR(M,10)*DENO(ID2+L,JA2+I)
C
                  QMAT(IC2+K,JB2+J) = QMAT(IC2+K,JB2+J)
     &                                - BCFF*RR(M,7 )*DENO(ID1+L,JA1+I)
     &                                - BCFF*RR(M,16)*DENO(ID2+L,JA2+I)
                ENDDO
              ENDDO
            ENDIF
C
C           TENTH IFLG BATCH
            IF(IFLG(10).EQ.1) THEN
              M = 0
              DO K=1,NFUNC
                DO L=1,NFUND
                  M = M+1
C
                  QMAT(IA1+I,JB1+J) = QMAT(IA1+I,JB1+J)
     &                                + ACFF*RR(M,1 )*DENO(IC1+K,JD1+L)
     &                                + ACFF*RR(M,4 )*DENO(IC2+K,JD2+L)
C
                  QMAT(IA2+I,JB2+J) = QMAT(IA2+I,JB2+J)
     &                                + ACFF*RR(M,13)*DENO(IC1+K,JD1+L)
     &                                + ACFF*RR(M,16)*DENO(IC2+K,JD2+L)
                ENDDO
              ENDDO
            ENDIF
C
C           ELEVENTH IFLG BATCH
            IF(IFLG(11).EQ.1) THEN
              M = 0
              DO K=1,NFUNC
                DO L=1,NFUND
                  M = M+1
C
                  QMAT(IC1+K,JD1+L) = QMAT(IC1+K,JD1+L)
     &                                + ACFF*RR(M,1 )*DENO(IA1+I,JB1+J)
     &                                + ACFF*RR(M,13)*DENO(IA2+I,JB2+J)
C
                  QMAT(IC2+K,JD2+L) = QMAT(IC2+K,JD2+L)
     &                                + ACFF*RR(M,4 )*DENO(IA1+I,JB1+J)
     &                                + ACFF*RR(M,16)*DENO(IA2+I,JB2+J)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ELSE
C         WRITE GENERAL CODE...
        ENDIF
C*********************************************************************C
C     END OF Q-MATRIX CONSTRUCTION                                    C
C*********************************************************************C
C
C     CLOSE ALL LOOPS
4000  CONTINUE
2999  CONTINUE
3999  CONTINUE
3000  CONTINUE
2000  CONTINUE
1000  CONTINUE
C
C*********************************************************************C
C     COMPLETE CONSTRUCTION OF FOCK MATRIX BY MATRIX CONJUGATION.     C
C     (THIS IS ONLY COMPLETE FOR THE DIATOMIC CASE.)                  C
C*********************************************************************C
C
C     DON'T DO THIS IF IT IS THE FIRST ITERATION OR ATOMIC
      IF((ITER.LE.1.OR.IRUN.NE.0).AND.NCENT.EQ.1) GOTO 5000
C
C     NON-RELATIVISTIC CALCULATIONS, DIMENSION NDIM
      IF(TREE.EQ.'NORL') THEN
        DO IDENT=2,NDIG
          DO JL=LDIAG(IDENT)+1,NDIM
            DO IL=LDIAG(IDENT-1)+1,LDIAG(IDENT)
              FOCK(IL,JL) = FOCK(IL,JL) + FOCK(JL,IL)
              FOCK(JL,IL) = FOCK(IL,JL)
            ENDDO
          ENDDO
        ENDDO
C     RELATIVISTIC CALCULATIONS, DIMENSION 2*NDIM
      ELSEIF(TREE.EQ.'DHFR'.OR.TREE.EQ.'DHFB') THEN
        DO IDENT=2,NDIG
          DO JL=LDIAG(IDENT)+1,NSHIFT
            DO IL=LDIAG(IDENT-1)+1,LDIAG(IDENT)
CC        THIS CODE COMES FROM BERTHA
C         DO JL=1,NSHIFT
C           DO IL=1,JL
C             IF((ICNLAB(IL).NE.ICNLAB(JL)).OR.(KAPLAB(IL).NE.KAPLAB(JL)).OR.(IMLAB(IL).NE.IMLAB(JL))) THEN
              IS = IL + NSHIFT
              JS = JL + NSHIFT
C
              FOCK(IL,JL) = FOCK(IL,JL) + DCONJG(FOCK(JL,IL))
              FOCK(JL,IL) =               DCONJG(FOCK(IL,JL))
              FOCK(IL,JS) = FOCK(IL,JS) + DCONJG(FOCK(JS,IL))
              FOCK(JS,IL) =               DCONJG(FOCK(IL,JS))
              FOCK(IS,JL) = FOCK(IS,JL) + DCONJG(FOCK(JL,IS))
              FOCK(JL,IS) =               DCONJG(FOCK(IS,JL))
              FOCK(IS,JS) = FOCK(IS,JS) + DCONJG(FOCK(JS,IS))
              FOCK(JS,IS) =               DCONJG(FOCK(IS,JS))
C             ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
5000  CONTINUE
C
C*********************************************************************C
C     COMPLETE CONSTRUCTION OF Q-MATRIX AND RELATED OBJECTS.          C
C     THIS CODE APPLIES ONLY TO OPEN-SHELL SYSTEMS                    C
C*********************************************************************C
C
      IF(NOPEN.EQ.0) GOTO 5001
C
C     COMPLETE CONSTRUCTION OF THE Q-MATRIX AND CLOSED-OPEN
C     COUPLING MATRICES FOR OPEN-SHELL SYSTEMS
C
C     NON-RELATIVISTIC CALCULATIONS, DIMENSION NDIM
      IF(TREE.EQ.'NORL') THEN
        DO IDENT=2,NDIG
          DO JL=LDIAG(IDENT)+1,NDIM
            DO IL=LDIAG(IDENT-1)+1,LDIAG(IDENT)
              QMAT(IL,JL) = QMAT(IL,JL) + QMAT(JL,IL)
              QMAT(JL,IL) = QMAT(IL,JL)
            ENDDO
          ENDDO
        ENDDO
C     RELATIVISTIC CALCULATIONS, DIMENSION 2*NDIM
      ELSEIF(TREE.EQ.'DHFR'.OR.TREE.EQ.'DHFB') THEN
        DO IDENT=2,NDIG
          DO J=LDIAG(IDENT)+1,NSHIFT
            DO I=LDIAG(IDENT-1)+1,LDIAG(IDENT)
              IS = IL + NSHIFT
              JS = JL + NSHIFT
C
              QMAT(IL,JL) = QMAT(IL,JL) + QMAT(JL,IL)
              QMAT(JL,IL) = QMAT(IL,JL)
              QMAT(IL,JS) = QMAT(IL,JS) + QMAT(JS,IL)
              QMAT(JS,IL) = QMAT(IL,JS)
              QMAT(IS,JL) = QMAT(IS,JL) + QMAT(JL,IS)
              QMAT(JL,IS) = QMAT(IS,JL)
              QMAT(IS,JS) = QMAT(IS,JS) + QMAT(JS,IS)
              QMAT(JS,IS) = QMAT(IS,JS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
C     SPECIFY THE OPEN-SHELL COUPLING MATRIX R AND ADD TO FOCK MATRIX
C     R = [S.D(O).Q + Q.D(O).S]  (RSCF 89)
      DO JL=1,NDIM
        DO KL=1,NDIM
          T1(KL) = 0.0D0
          T2(KL) = 0.0D0
          DO LL=1,NDIM
            T1(KL) = T1(KL) + DENT(KL,LL)*QMAT(LL,JL)
            T2(KL) = T2(KL) + DENT(KL,LL)*OVAP(LL,JL)
          ENDDO
        ENDDO
C
        DO IL=1,NDIM
          T3(IL) = 0.0D0
          DO KL=1,NDIM
            T3(IL) = T3(IL) + OVAP(IL,KL)*T1(KL) + QMAT(IL,KL)*T2(KL)
          ENDDO
        ENDDO
C
C       ADD THE PROJECTOR AND Q-MATRIX TO THE FOCK MATRIX  (RSCF 91)
        DO IL=1,NDIM 
          FOCK(IL,JL) = FOCK(IL,JL) - QMAT(IL,JL) + T3(IL)
        ENDDO
      ENDDO
C
5001  CONTINUE
C
C*********************************************************************C
C     CALCULATE THE TOTAL ENERGY AND PART THAT COMES FROM QMAT        C
C*********************************************************************C
C
C     INITIALISE THE TOTAL ENERGY BIN
      ETMP = 0.0D0
C
C     GENERATE TOTAL ENERGY FROM FOCK AND QMAT MATRICES
      IF(NOPEN.EQ.0) THEN
        DO IL=1,NDIM
          DO JL=1,NDIM
            ETMP = ETMP + 0.5D0*DENT(IL,JL)*FOCK(IL,JL)
          ENDDO
        ENDDO
C
C     CONTRIBUTIONS FROM Q MATRIX AND PROJECTION MATRIX
      ELSEIF(NOPEN.NE.0) THEN
        DO IL=1,NDIM
          DO JL=1,NDIM
            ETMP = ETMP - 0.5D0*QMAT(IL,JL)*DENO(IL,JL)
     &                  + 0.5D0*QMAT(IL,JL)*DENT(IL,JL)*(FOPEN-1.0D0)
          ENDDO
        ENDDO
      ENDIF
C
C     UPDATE TWO-ELECTRON MEAN-FIELD COULOMB ENERGY
      ECLM = ETMP
C
      RETURN
      END
C
C      
      SUBROUTINE ERI(RR,XYZ,KQN,MQN,EXPT,NFUNS,ITQN,I,J,IEAB,IECD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C                       EEEEEEEE RRRRRRR  IIIII                       C
C                       EE       RR    RR  II                         C
C                       EE       RR    RR  II                         C
C                       EEEEEE   RR    RR  II                         C
C                       EE       RRRRRRR   II                         C
C                       EE       RR    RR  II                         C
C                       EEEEEEEE RR    RR IIIII                       C
C                                                                     C
C ------------------------------------------------------------------- C
C     ERI GENERATES BLOCKS OF ELECTRON REPULSION INTEGRALS            C
C     OVER KINETICALLY BALANCED G-SPINOR BASIS FUNCTIONS              C
C                                                                     C
C     THE DENSITIES ARE EXPANDED IN A BASIS OF HERMITE GAUSSIANS      C
C     AND THE INTEGRALS ARE GENERATED USING THE MCMURCHIE-            C
C     DAVIDSON ALGORITHM                                              C
C ------------------------------------------------------------------- C
C                      INPUT PARAMETERS                               C
C                                                                     C
C     XYZ(3,4)    COORDINATES OF THE 4 NUCLEAR CENTRES                C
C     KQN(4)      KQN QUANTUM NUMBERS OF THE CENTRES                  C
C     MQN(4)      |MQN| QUANTUM NUMBERS OF THE CENTRES                C
C     EXPT(I,J)   EXPONENTS ON CENTRE J                               C
C     NFUNS(J)    NUMBER OF FUNCTIONS ON CENTRE J                     C
C     ITQN(2)     COMPONENT PAIRS: ITQN(I)=1 - > LL                   C
C                                  ITQN(I)=2 - > LS                   C
C                                  ITQN(I)=3 - > SL                   C
C                                  ITQN(I)=4 - > SS                   C
C                                  I=1       - > AB                   C
C                                  I=2       - > CD                   C
C     I,J         INDEX FOR BASIS FUNCTION PAIR ON AB                 C
C     IEAB,IECD   0 DON'T RECALCULATE E(AB/CD)-COEFFICIENTS           C
C                 1 DO    RECALCULATE E(AB/CD)-COEFFICIENTS           C
C ------------------------------------------------------------------- C
C     NOTE: THERE ARE NO SHORTCUTS FOR INUCAB OR INUCCD CONDITIONS.   C
C           IF IATOM = 1 (4 BASIS FUNCTIONS ON SAME CENTRE), SHOULD   C
C           EVALUATE WITH RACAH ALGEBRA INSTEAD!                      C
C*********************************************************************C
      PARAMETER (MXDM=1600,NCENTM=6,MAXB=26,MAXB2=MAXB*MAXB,NKAPM=9,
     &           LMAX=(NKAPM-1)/2,LAMMX=2*LMAX,IL4=2*LAMMX,MAXMV=NKAPM,
     &           MAXMV2=2*MAXMV,MLL=((NKAPM+8)*(NKAPM+9)*(NKAPM+10))/6,
     &           MAXR=969)
C
      COMPLEX*16 CONE
      COMPLEX*16 RR(MAXB2,16),Q1(MAXB2),Q2(MAXB2)
      COMPLEX*16 EAB11(MAXB2,MLL),ECD11(MAXB2,MLL),GAB11(MAXB2,MLL),
     &           EAB21(MAXB2,MLL),ECD21(MAXB2,MLL),GAB21(MAXB2,MLL)
C
      COMMON/ACCSS/INABCD(0:IL4,0:IL4,0:IL4),
     &             IVEC(MLL),JVEC(MLL),KVEC(MLL),LAMVEC(MLL)
      COMMON/ESTOR/EAB11,EAB21
C
      DIMENSION XYZ(3,4),RC(MAXB2,MAXR),
     &          PQ(MAXB2,3),EXPT(MAXB,4),ALPHA(MAXB2),PREFAC(MAXB2)
      DIMENSION KQN(4),LQN(4),MQN(4),ITQN(2),NFUNS(4)
      DIMENSION IAB11(MLL),IAB21(MLL),ICD11(MLL),ICD21(MLL),IRC(MAXR)
C
      DATA ROOTPI,SENS/1.7724538509055160D0,1.0D-14/
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     EVALUATE LQNS FOR BASIS FUNCTIONS A, B, C, D 
      DO L=1,4
        IF(KQN(L).LT.0) THEN
          LQN(L) =-KQN(L)-1
        ELSE
          LQN(L) = KQN(L)
        ENDIF
      ENDDO
C
C     COUNTER FOR BRA OR KET OVERLAPS
      IALTAB = 1
      IALTCD =-1
C
C     CALCULATE FINITE SUM TRUNCATION VALUES DEPENDING ON OVERLAP TYPE
      IF(ITQN(1).EQ.1.AND.ITQN(2).EQ.1) THEN
        LAMAB = LQN(1) + LQN(2)
        LAMCD = LQN(3) + LQN(4)
      ELSEIF(ITQN(1).EQ.1.AND.ITQN(2).EQ.4) THEN
        LAMAB = LQN(1) + LQN(2)
        LAMCD = LQN(3) + LQN(4) + 2
      ELSEIF(ITQN(1).EQ.4.AND.ITQN(2).EQ.1) THEN
        LAMAB = LQN(1) + LQN(2) + 2
        LAMCD = LQN(3) + LQN(4)
      ELSEIF(ITQN(1).EQ.4.AND.ITQN(2).EQ.4) THEN
        LAMAB = LQN(1) + LQN(2) + 2
        LAMCD = LQN(3) + LQN(4) + 2
      ELSE
        WRITE(6,*) 'In ERI: incorrect call.'
        STOP
      ENDIF
C
C     NUMBER OF FINITE EXPANSION INDICES FOR EACH CENTRE
      NTUVAB = ((LAMAB+1)*(LAMAB+2)*(LAMAB+3))/6
      NTUVCD = ((LAMCD+1)*(LAMCD+2)*(LAMCD+3))/6
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB  = NFUNS(1)*NFUNS(2)
      MAXCD  = NFUNS(3)*NFUNS(4)     
C
C*********************************************************************C
C     IF ASKED TO RECALCULATE E(AB) COEFFICIENTS, DO THIS FIRST       C
C*********************************************************************C
C
      IF(IEAB.EQ.1) THEN
C
C       AB PAIRS      
        CALL RNORMF(EXPT,LQN,NFUNS,1,2)
C        
        IF (ITQN(1).EQ.1) THEN
          CALL EMAKELL(EAB11,EAB21,EXPT,KQN,MQN,NFUNS,IALTAB,1,2,0)
        ELSEIF(ITQN(1).EQ.4) THEN                         
          CALL EMAKESS(EAB11,EAB21,EXPT,KQN,MQN,NFUNS,IALTAB,1,2,0)
        ENDIF
C
C       SCREENING PROCEDURE: TEST MAGNITUDES OF E-COEFFICIENT LISTS
        DO IAB=1,NTUVAB
C
C         11 OVERLAP (AB PAIRS)
          TEST = CDASUM(MAXAB,EAB11(1,IAB))
          IF(TEST.LE.SENS) THEN
            IAB11(IAB) = 0
          ELSE
            IAB11(IAB) = 1
          ENDIF
C
C         21 OVERLAP (AB PAIRS)
          TEST = CDASUM(MAXAB,EAB21(1,IAB))
          IF(TEST.LE.SENS) THEN
            IAB21(IAB) = 0
          ELSE
            IAB21(IAB) = 1
          ENDIF
        
        ENDDO
C
C       DO NOT CALCULATE AGAIN UNTIL ASKED EXTERNALLY
        IEAB = 0
C
      ENDIF
C
C*********************************************************************C
C     IF ASKED TO RECALCULATE E(CD) COEFFICIENTS, DO THIS NEXT        C
C*********************************************************************C
C
      IF(IECD.EQ.1) THEN
C
C       CD PAIRS
        CALL RNORMF(EXPT,LQN,NFUNS,3,4)
C        
        IF(ITQN(2).EQ.1) THEN
          CALL EMAKELL(ECD11,ECD21,EXPT,KQN,MQN,NFUNS,IALTCD,3,4,0)
        ELSEIF(ITQN(2).EQ.4) THEN
          CALL EMAKESS(ECD11,ECD21,EXPT,KQN,MQN,NFUNS,IALTCD,3,4,0)
        ENDIF
C
C       SCREENING PROCEDURE: TEST MAGNITUDES OF E-COEFFICIENT LISTS       
        DO ICD=1,NTUVCD
C
C         11 OVERLAP (CD PAIRS)
          TEST = CDASUM(MAXCD,ECD11(1,ICD))
          IF(TEST.LE.SENS) THEN
            ICD11(ICD) = 0
          ELSE
            ICD11(ICD) = 1
          ENDIF
C
C         21 OVERLAP (CD PAIRS)
          TEST = CDASUM(MAXCD,ECD21(1,ICD))
          IF(TEST.LE.SENS) THEN
            ICD21(ICD) = 0
          ELSE
            ICD21(ICD) = 1
          ENDIF
C
        ENDDO
C
C       DO NOT CALCULATE AGAIN UNTIL ASKED EXTERNALLY
        IECD = 0
C
      ENDIF
C
C*********************************************************************C
C     R-INTEGRAL EVALUATION                                           C
C*********************************************************************C
C
C     GAUSSIAN OVERLAP VALUES  
      EIJ = EXPT(I,1) + EXPT(J,2)
      PX  = (XYZ(1,1)*EXPT(I,1)+XYZ(1,2)*EXPT(J,2))/EIJ
      PY  = (XYZ(2,1)*EXPT(I,1)+XYZ(2,2)*EXPT(J,2))/EIJ
      PZ  = (XYZ(3,1)*EXPT(I,1)+XYZ(3,2)*EXPT(J,2))/EIJ
C
      M = 0
      DO K=1,NFUNS(3)
        DO L=1,NFUNS(4)
          M = M+1
          
          EKL = EXPT(K,3)+EXPT(L,4)
          QX  = (XYZ(1,3)*EXPT(K,3)+XYZ(1,4)*EXPT(L,4))/EKL
          QY  = (XYZ(2,3)*EXPT(K,3)+XYZ(2,4)*EXPT(L,4))/EKL
          QZ  = (XYZ(3,3)*EXPT(K,3)+XYZ(3,4)*EXPT(L,4))/EKL
          
          ALPHA(M)  = (EIJ*EKL)/(EIJ+EKL)
          PQ(M,1)   = QX-PX
          PQ(M,2)   = QY-PY
          PQ(M,3)   = QZ-PZ
          PREFAC(M) = 2.0D0*(ROOTPI**5)/(DSQRT(EIJ+EKL)*EIJ*EKL)
        ENDDO
      ENDDO
C      
      MAXCD = NFUNS(3)*NFUNS(4)
C
      CALL RMAKE(RC,PQ,ALPHA,MAXCD,LAMAB+LAMCD)
C
C     INITIALIZE ARRAY TO IMPLEMENT SPARSENESS IN R-VECTOR
      LAMABCD  = LAMAB+LAMCD
      NTUVABCD = ((LAMABCD+1)*(LAMABCD+2)*(LAMABCD+3))/6
C      
      DO MRC=1,NTUVABCD
        TEST = DASUM(MAXCD,RC(1,MRC),1)
        IF(TEST.LE.SENS) THEN
          IRC(MRC) = 0
        ELSE
          IRC(MRC) = 1
        ENDIF
      ENDDO     
C
C*********************************************************************C
C     CONSTRUCT INTERMEDIATE MATRICES FOR MCMURCHIE-DAVIDSON          C
C*********************************************************************C
C
      DO IAB=1,NTUVAB
C
        IAB1 = 0
        IAB2 = 0
C
C       SCREENING MARKERS
        IF(IAB11(IAB).EQ.1) THEN
          IAB1 = 1
        ENDIF
C
        IF(IAB21(IAB).EQ.1) THEN
          IAB2 = 1
        ENDIF
C
        IF(IAB1+IAB2.GT.0) THEN
          DO M=1,MAXCD
            GAB11(M,IAB) = DCMPLX(0.0D0,0.0D0)
            GAB21(M,IAB) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDIF
C
C       CALCULATE OVERALL ADDRESS AND ADD TO THE G-ARRAY BINS
        DO ICD=1,NTUVCD
          IRABCD = INABCD(IVEC(IAB)+IVEC(ICD),
     &                    JVEC(IAB)+JVEC(ICD),
     &                    KVEC(IAB)+KVEC(ICD))
C
C         SKIP THIS STEP IF THE R-INTEGRAL IS NOT BIG ENOUGH
          IF(IRC(IRABCD).EQ.0) GOTO 798
C
          IF(ICD11(ICD).EQ.1) THEN
            DO M=1,MAXCD
              GAB11(M,IAB) = GAB11(M,IAB) + ECD11(M,ICD)*RC(M,IRABCD)
            ENDDO
          ENDIF

          IF(ICD21(ICD).EQ.1) THEN
            DO M=1,MAXCD
              GAB21(M,IAB) = GAB21(M,IAB) + ECD21(M,ICD)*RC(M,IRABCD)
            ENDDO
          ENDIF
C 
798       CONTINUE
        ENDDO
      ENDDO
C
C*********************************************************************C
C     GENERATE ALL POSSIBLE TWO-ELECTRON INTEGRALS FROM THE           C
C     EAB COEFFICIENTS AND THE G-ARRAYS                               C
C*********************************************************************C
C
C     CALCULATE PHASES FOR BASIS FUNCTION OVERLAP COMBINATIONS
      P1 = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
      P2 = DFLOAT((-1)**((MQN(3)-MQN(4))/2))
      P1 = P1*DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
      P2 = P2*DFLOAT((KQN(3)*KQN(4))/IABS(KQN(3)*KQN(4)))
      P3 = P1*P2
C
C*********************************************************************C
C     INTEGRAL BATCH 1: ( - - || - - )                                C
C*********************************************************************C
C
      NTUVAB = ((LAMAB+1)*(LAMAB+2)*(LAMAB+3))/6
      IJ     = (I-1)*NFUNS(2) + J
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSENESS 
      DO M=1,MAXCD
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {ALPH,BETA,GAMA}
      DO IAB=1,NTUVAB
        IF(IAB11(IAB).EQ.1) THEN
          DO M=1,MAXCD
            Q1(M) = Q1(M) +      EAB11(IJ,IAB)*DREAL(GAB11(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB11(IJ,IAB)*DIMAG(GAB11(M,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO M=1,MAXCD
        RR(M,1 ) =    (Q1(M)+Q2(M))*PREFAC(M)
        RR(M,4 ) = P2*(Q1(M)-Q2(M))*PREFAC(M)
        RR(M,13) = P3*DCONJG(RR(M,4))
        RR(M,16) = P3*DCONJG(RR(M,1))
      ENDDO
C
C*********************************************************************C
C     INTEGRAL BATCH 2: ( - - || + - )                                C
C*********************************************************************C
C
      NTUVAB = ((LAMAB+1)*(LAMAB+2)*(LAMAB+3))/6
      IJ     = (I-1)*NFUNS(2) + J
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSENESS 
      DO M=1,MAXCD
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {ALPH,BETA,GAMA}
      DO IAB=1,NTUVAB
        IF(IAB11(IAB).EQ.1) THEN
          DO M=1,MAXCD
            Q1(M) = Q1(M) +      EAB11(IJ,IAB)*DREAL(GAB21(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB11(IJ,IAB)*DIMAG(GAB21(M,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO M=1,MAXCD
        RR(M,3 ) =    (Q1(M)+Q2(M))*PREFAC(M)
        RR(M,2 ) =-P2*(Q1(M)-Q2(M))*PREFAC(M)
        RR(M,15) =-P3*DCONJG(RR(M,2))
        RR(M,14) =-P3*DCONJG(RR(M,3))
      ENDDO
C
C*********************************************************************C
C     INTEGRAL BATCH 3: ( + - || - - )                                C
C*********************************************************************C
C
      NTUVAB = ((LAMAB+1)*(LAMAB+2)*(LAMAB+3))/6
      IJ     = (I-1)*NFUNS(2) + J
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSENESS 
      DO M=1,MAXCD
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {ALPH,BETA,GAMA}
      DO IAB=1,NTUVAB
        IF(IAB21(IAB).EQ.1) THEN
          DO M=1,MAXCD
            Q1(M) = Q1(M) +      EAB21(IJ,IAB)*DREAL(GAB11(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB21(IJ,IAB)*DIMAG(GAB11(M,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO M=1,MAXCD
        RR(M,9 ) =    (Q1(M)+Q2(M))*PREFAC(M)
        RR(M,12) = P2*(Q1(M)-Q2(M))*PREFAC(M)
        RR(M,5 ) =-P3*DCONJG(RR(M,12))
        RR(M,8 ) =-P3*DCONJG(RR(M,9 ))      
      ENDDO
C
C*********************************************************************C
C     INTEGRAL BATCH 4: ( + - || + - )                                C
C*********************************************************************C
C
      NTUVAB = ((LAMAB+1)*(LAMAB+2)*(LAMAB+3))/6
      IJ  = (I-1)*NFUNS(2) + J
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSENESS 
      DO M=1,MAXCD
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {ALPH,BETA,GAMA}
      DO IAB=1,NTUVAB
        IF(IAB21(IAB).EQ.1) THEN
          DO M=1,MAXCD
            Q1(M) = Q1(M) +      EAB21(IJ,IAB)*DREAL(GAB21(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB21(IJ,IAB)*DIMAG(GAB21(M,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO M=1,MAXCD
        RR(M,11) =    (Q1(M)+Q2(M))*PREFAC(M)
        RR(M,10) =-P2*(Q1(M)-Q2(M))*PREFAC(M)
        RR(M,7 ) = P3*DCONJG(RR(M,10))
        RR(M,6 ) = P3*DCONJG(RR(M,11))      
      ENDDO
C
      RETURN
      END

