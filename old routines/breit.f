      SUBROUTINE BREIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C              BBBBBBB  RRRRRRR  EEEEEEEE IIII TTTTTTTT               C
C              BB    BB RR    RR EE        II     TT                  C
C              BB    BB RR    RR EE        II     TT                  C
C              BBBBBBB  RRRRRRR  EEEEEE    II     TT                  C
C              BB    BB RR   RR  EE        II     TT                  C
C              BB    BB RR    RR EE        II     TT                  C
C              BBBBBBB  RR    RR EEEEEEEE IIII    TT                  C
C                                                                     C
C ------------------------------------------------------------------- C
C    BREIT GENERATES ELECTRON INTERACTION INTEGRALS IN BATCHES AND    C
C    CALCULATES THE SCF BREIT MATRIX (B). INTEGRAL SYMMETRY IS        C
C    PARTIALLY EXPLOITED, BUT WITH ROOM FOR IMPROVEMENT (GEOMETRIC    C
C    SYMM, R-INT SYMM, E-COEFF SYMM).                                 C
C ------------------------------------------------------------------- C
C    NOTE: THIS ROUTINE COULD BENEFIT FROM PARALLELISATION -- OPENMP. C
C          LOOK INTO OPEN SHELL EXTENSIONS AND CASES WITH ZERO DIRECT.C
C*********************************************************************C
      PARAMETER (MXDM=1600,NCENTM=6,MAXB=26,MAXB2=MAXB*MAXB,NKAPM=9,
     &           LMAX=(NKAPM-1)/2,LAMMX=2*LMAX,IL4=2*LAMMX,MAXMV=NKAPM,
     &           MAXMV2=2*MAXMV,MAXR=969)
C
      COMPLEX*16 DENC(MXDM,MXDM),DENO(MXDM,MXDM),DENT(MXDM,MXDM)
      COMPLEX*16 FOCK(MXDM,MXDM),OVAP(MXDM,MXDM),QMAT(MXDM,MXDM)
      COMPLEX*16 BMAT(MXDM,MXDM)
      COMPLEX*16 RR(MAXB2,16)
C
      COMMON/BRTVAL/BMAT
      COMMON/BSCR/BRX(MXDM,MXDM),IBINDX(NCENTM,NKAPM)
      COMMON/DENS/DENC,DENO,DENT
      COMMON/ENERGY/ETOT,ENUC,EONE,ECLM,EBRT
      COMMON/FMAT/FOCK,OVAP,QMAT
      COMMON/IBSC/IBSCR(MAXB2),IBMAP(MAXB2)
      COMMON/LABELS/LARGE(NCENTM,NKAPM,MAXMV2)
      COMMON/LOCGEO/XYZ(3,4)
      COMMON/SHELL/ALPHA,BETA,FOPEN,ICLOSE(100),IOPEN(6),NCLOSE,NOPEN
      COMMON/SPEC/COORD(3,NCENTM),ZNUC(NCENTM),CNUC(NCENTM),
     &            EXPSET(MAXB,NKAPM,NCENTM),IZNUC(NCENTM),ICRGE(NCENTM),
     &            KVALS(NKAPM,NCENTM),NKAPPA(NCENTM),LMAXX(NCENTM),
     &            NFUNCT(NKAPM,NCENTM),NCENT,NDIM,NSHIFT,
     &            NOCC,ITER,IALL,IRUN,IOCCM0
      COMMON/SYMARR/KAPLAB(MXDM),ICNLAB(MXDM),IMLAB(MXDM)
C
      DIMENSION EXPT(MAXB,4),KQN(4),LQN(4),MQN(4),NFUNS(4)
      DIMENSION INDEX(NCENTM,-NKAPM:NKAPM,2*(NKAPM+1)*NKAPM)
      DIMENSION RMX(MAXB2)
C
C     INITIALISE TWO-ELECTRON BREIT ENERGY COUNTER
      EBRT = 0.0D0
C
C     INITIALISE THE BREIT MATRIX
      DO J=1,NDIM
        DO I=1,NDIM
          BMAT(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C*********************************************************************C
C     INDEXING ROUTINE: SO WE CAN SET BASIS FUNCTION LABELS BY BLOCK  C
C*********************************************************************C
      ICOUNT = 0
      IKOUNT = 0
C
C     LOOP OVER NUCLEAR CENTRES
      DO ICENT=1,NCENT
C
C       LOOP OVER KAPPA VALUES FOR THIS NUCLEAR CENTRE
        DO KN=1,NKAPPA(ICENT)
C
C         IMPORT KAPPA, LQN, MAXIMUM MQN AND NUMBER OF BASIS FUNCTIONS
          KAPPA = KVALS(KN,ICENT)
          MJMAX = 2*IABS(KAPPA)-1
          IF(KAPPA.GT.0) THEN
            LQNN = KAPPA
          ELSE
            LQNN =-KAPPA-1
          ENDIF
          NFUN = NFUNCT(LQNN+1,ICENT)
C
C         RECORD CUMULATIVE INDEX COUNTER FOR NUMBER OF BASIS FUNCTIONS
          IBINDX(ICENT,KN) = IKOUNT
          IKOUNT           = IKOUNT + NFUN
C
C         LOOP OVER MQN VALUES AND RECORD INDEX (WILL BE MORE THAN IKOUNT)
          DO MJ=1,MJMAX,2
            ICOUNT                = ICOUNT+1
            INDEX(ICENT,KAPPA,MJ) = ICOUNT
            LENIQ                 = ICOUNT
          ENDDO
        ENDDO
      ENDDO
C
C*********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)     C
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
      DO 2000 ICENTB=1,NCENT
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
        KQN(2)=KVALS(KB,ICENTB)
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
        MJA    = (2*MA)-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 2000 MB=1,IABS(KQN(2))
        MJB    = (2*MB)-1
        MQN(2) = MJB
C
C*********************************************************************C
C     SECOND LAYER OF LOOPS, OVER CENTRES C AND D (USE INDEX 1000)    C
C*********************************************************************C
C
C     LOOP OVER CENTRE C
      DO 1000 ICENTC=1,NCENT
C
C       CARTESIAN COORDINATES OF CENTRE C
        XYZ(1,3) = COORD(1,ICENTC)
        XYZ(2,3) = COORD(2,ICENTC)
        XYZ(3,3) = COORD(3,ICENTC)
C
C     LOOP OVER CENTRE D
      DO 1000 ICENTD=1,NCENT
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
      DO 1000 KC=1,NKAPPA(ICENTC)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR C
        KQN(3)=KVALS(KC,ICENTC)
        IF(KQN(3).GT.0) THEN
          LQN(3)= KQN(3)
        ELSE
          LQN(3)=-KQN(3)-1
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
      DO 1000 KD=1,NKAPPA(ICENTD)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR D
        KQN(4)=KVALS(KD,ICENTD)
        IF(KQN(4).GT.0) THEN
         LQN(4)= KQN(4)
        ELSE
         LQN(4)=-KQN(4)-1
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
      DO 1000 MC=1,IABS(KQN(3))
        MJC    = (2*MC)-1
        MQN(3) = MJC
C
C     LOOP OVER |MQN(D)| VALUES
      DO 1000 MD=1,IABS(KQN(4))
        MJD    = (2*MD)-1
        MQN(4) = MJD
C
C*********************************************************************C
C     FOR THIS CHOICE OF A,B,C AND D, COMPUTE INDICES AND PHASES      C
C*********************************************************************C
C
C     STARTING INDEX VALUES
      IQ1 = INDEX(ICENTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICENTB,KQN(2),MQN(2))
      IQ3 = INDEX(ICENTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICENTD,KQN(4),MQN(4))
C
C     LENIQ IS THE NUMBER OF BLOCKS TO BE CALCULATED
      IQ12 = (IQ1-1)*LENIQ + IQ2
      IQ34 = (IQ3-1)*LENIQ + IQ4
C
      IF (IQ12.LT.IQ34) GOTO 1999
C
C     IN SCF THIS IS WHERE WE WOULD SET THE COMBINATION LIMITS BASED
C     ON THE TREE, BUT BREIT DISAPPEARS IN NON-REL LIMIT SO NO NEED.
C
C*********************************************************************C
C     SHOULD BE ABLE TO APPLY SELECTION RULES BASED ON INUCAB,        C
C     INUCCD, IATOM AND MOLECULAR GEOMETRY HERE (LOOK AT DIATOM_02.F) C
C*********************************************************************C
C
C     CALCULATE KQN PHASE FACTORS FOR PERMUTING INTEGRALS.
C     NB! OPPOSITE PHASE AS IN THE LL/SS CASE SEEN IN SCF
      IF((KQN(1)*KQN(2)).GT.0) THEN 
        RKFAC1 =-1.0D0
      ELSE
        RKFAC1 = 1.0D0
      ENDIF
C
      IF((KQN(3)*KQN(4)).GT.0) THEN 
        RKFAC2 =-1.0D0
      ELSE
        RKFAC2 = 1.0D0
      ENDIF
C
C     CALCULATE MQN PHASE FACTORS FOR PERMUTING INTEGRALS
      FACAB1 = DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      FACAB2 = DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      FACCD1 = DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      FACCD2 = DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
C     STARTING ADDRESSES FOR BASIS OVERLAPS
      IA1 = LARGE(ICENTA,KA,MJA  )
      IB1 = LARGE(ICENTB,KB,MJB  ) + NSHIFT
      IC1 = LARGE(ICENTC,KC,MJC  )
      ID1 = LARGE(ICENTD,KD,MJD  ) + NSHIFT
C
      IA2 = LARGE(ICENTA,KA,MJA+1)
      IB2 = LARGE(ICENTB,KB,MJB+1) + NSHIFT
      IC2 = LARGE(ICENTC,KC,MJC+1)
      ID2 = LARGE(ICENTD,KD,MJD+1) + NSHIFT
C
      JA1 = LARGE(ICENTA,KA,MJA  )
      JB1 = LARGE(ICENTB,KB,MJB  ) + NSHIFT
      JC1 = LARGE(ICENTC,KC,MJC  )
      JD1 = LARGE(ICENTD,KD,MJD  ) + NSHIFT
C
      JA2 = LARGE(ICENTA,KA,MJA+1)
      JB2 = LARGE(ICENTB,KB,MJB+1) + NSHIFT
      JC2 = LARGE(ICENTC,KC,MJC+1)
      JD2 = LARGE(ICENTD,KD,MJD+1) + NSHIFT
C
C     NOT SURE WHAT THESE ARE FOR, OR WHETHER WE REALLY NEED COMMON/BSCR
      IAA = IBINDX(ICENTA,KA)
      JBB = IBINDX(ICENTB,KB) + NSHIFT
C
      ICC = IBINDX(ICENTC,KC)
      JDD = IBINDX(ICENTD,KD) + NSHIFT
C
C     TOTAL NUMBER OF BASIS FUNCTION OVERLAPS     
      MAXMAB  = NFUNA*NFUNB
      MAXMCD  = NFUNC*NFUND
C
C     THESE COUNTERS TELL THE BRTINT PROCEDURE TO GENERATE NEW E-COEFFS
      IEMAKE = 1
C
C*********************************************************************C
C     THIRD LAYER OF LOOPS, OVER BASIS FUNCTIONS A AND B (USE 3000)   C
C*********************************************************************C 
C
      DO 3000 I=1,NFUNA
      DO 3000 J=1,NFUNB
C
C
C*********************************************************************C
C     OVERRIDE SCREENING TESTS FOR THE SAKE OF CERTAINTY              C
C*********************************************************************C
C
      DO M=1,MAXMCD
        IBMAP(M) = M
        IBSCR(M) = 1
      ENDDO
C
C     GENERATE A BATCH OF BREIT INTERACTION INTEGRALS
      CALL BII(RR,XYZ,KQN,MQN,EXPT,NFUNS,I,J,IEMAKE,INUCAB,INUCCD)
C
C     E-COEFFICIENTS HAVE NOW BEEN GENERATED SO TURN OFF SWITCH
      IEMAKE = 0
C
C     FIRST BATCH
      F1 = RKFAC2*FACCD1
      F2 = RKFAC2*FACCD2
      M  = 0
      DO K=1,NFUNC
        DO L=1,NFUND
          M = M+1
C
          BMAT(IA1+I,JB1+J) = BMAT(IA1+I,JB1+J)
     &         + RR(M,1 )*(DENT(IC1+K,JD1+L) + F1*DENT(JD2+L,IC2+K))
     &         + RR(M,2 )*(DENT(IC1+K,JD2+L) + F2*DENT(JD1+L,IC2+K))
     &         + RR(M,3 )*(DENT(IC2+K,JD1+L) + F2*DENT(JD2+L,IC1+K))
     &         + RR(M,4 )*(DENT(IC2+K,JD2+L) + F1*DENT(JD1+L,IC1+K))
C
          BMAT(IA1+I,JB2+J) = BMAT(IA1+I,JB2+J)
     &         + RR(M,5 )*(DENT(IC1+K,JD1+L) + F1*DENT(JD2+L,IC2+K))
     &         + RR(M,6 )*(DENT(IC1+K,JD2+L) + F2*DENT(JD1+L,IC2+K))
     &         + RR(M,7 )*(DENT(IC2+K,JD1+L) + F2*DENT(JD2+L,IC1+K))
     &         + RR(M,8 )*(DENT(IC2+K,JD2+L) + F1*DENT(JD1+L,IC1+K))
C
          BMAT(IA2+I,JB1+J) = BMAT(IA2+I,JB1+J)
     &         + RR(M,9 )*(DENT(IC1+K,JD1+L) + F1*DENT(JD2+L,IC2+K))
     &         + RR(M,10)*(DENT(IC1+K,JD2+L) + F2*DENT(JD1+L,IC2+K))
     &         + RR(M,11)*(DENT(IC2+K,JD1+L) + F2*DENT(JD2+L,IC1+K))
     &         + RR(M,12)*(DENT(IC2+K,JD2+L) + F1*DENT(JD1+L,IC1+K))
C
          BMAT(IA2+I,JB2+J) = BMAT(IA2+I,JB2+J)
     &         + RR(M,13)*(DENT(IC1+K,JD1+L) + F1*DENT(JD2+L,IC2+K))
     &         + RR(M,14)*(DENT(IC1+K,JD2+L) + F2*DENT(JD1+L,IC2+K))
     &         + RR(M,15)*(DENT(IC2+K,JD1+L) + F2*DENT(JD2+L,IC1+K))
     &         + RR(M,16)*(DENT(IC2+K,JD2+L) + F1*DENT(JD1+L,IC1+K))
        ENDDO
      ENDDO
C
C     SECOND BATCH
      DO L=1,NFUND
        DO K=1,NFUNC
          M1 = (K-1)*NFUND + L
C
          BMAT(IA1+I,JD1+L) = BMAT(IA1+I,JD1+L)
     &                                - RR(M1,1 )*DENT(IC1+K,JB1+J)
     &                                - RR(M1,3 )*DENT(IC2+K,JB1+J)
     &                                - RR(M1,5 )*DENT(IC1+K,JB2+J)
     &                                - RR(M1,7 )*DENT(IC2+K,JB2+J)
C
          BMAT(IA1+I,JD2+L) = BMAT(IA1+I,JD2+L)
     &                                - RR(M1,2 )*DENT(IC1+K,JB1+J)
     &                                - RR(M1,4 )*DENT(IC2+K,JB1+J)
     &                                - RR(M1,6 )*DENT(IC1+K,JB2+J)
     &                                - RR(M1,8 )*DENT(IC2+K,JB2+J)
C
          BMAT(IA2+I,JD1+L) = BMAT(IA2+I,JD1+L)
     &                                - RR(M1,9 )*DENT(IC1+K,JB1+J)
     &                                - RR(M1,11)*DENT(IC2+K,JB1+J)
     &                                - RR(M1,13)*DENT(IC1+K,JB2+J)
     &                                - RR(M1,15)*DENT(IC2+K,JB2+J)
C
          BMAT(IA2+I,JD2+L) = BMAT(IA2+I,JD2+L)
     &                                - RR(M1,10)*DENT(IC1+K,JB1+J)
     &                                - RR(M1,12)*DENT(IC2+K,JB1+J)
     &                                - RR(M1,14)*DENT(IC1+K,JB2+J)
     &                                - RR(M1,16)*DENT(IC2+K,JB2+J)
        ENDDO
      ENDDO
C
C     THIRD BATCH
      F1 = RKFAC2*FACCD1
      F2 = RKFAC2*FACCD2
      M = 0
      DO K=1,NFUNC
        DO L=1,NFUND
          M = M+1
          BMAT(IA1+I,JC1+K) = BMAT(IA1+I,JC1+K)
     &                              - F1*RR(M,4 )*DENT(ID1+L,JB1+J)
     &                              - F1*RR(M,8 )*DENT(ID1+L,JB2+J)
     &                              - F2*RR(M,3 )*DENT(ID2+L,JB1+J)
     &                              - F2*RR(M,7 )*DENT(ID2+L,JB2+J)
          BMAT(IA1+I,JC2+K) = BMAT(IA1+I,JC2+K)
     &                              - F2*RR(M,2 )*DENT(ID1+L,JB1+J)
     &                              - F2*RR(M,6 )*DENT(ID1+L,JB2+J)
     &                              - F1*RR(M,1 )*DENT(ID2+L,JB1+J)
     &                              - F1*RR(M,5 )*DENT(ID2+L,JB2+J)
          BMAT(IA2+I,JC1+K) = BMAT(IA2+I,JC1+K)
     &                              - F1*RR(M,12)*DENT(ID1+L,JB1+J)
     &                              - F1*RR(M,16)*DENT(ID1+L,JB2+J)
     &                              - F2*RR(M,11)*DENT(ID2+L,JB1+J)
     &                              - F2*RR(M,15)*DENT(ID2+L,JB2+J)
          BMAT(IA2+I,JC2+K) = BMAT(IA2+I,JC2+K)
     &                              - F2*RR(M,10)*DENT(ID1+L,JB1+J)
     &                              - F2*RR(M,14)*DENT(ID1+L,JB2+J)
     &                              - F1*RR(M,9 )*DENT(ID2+L,JB1+J)
     &                              - F1*RR(M,13)*DENT(ID2+L,JB2+J)
        ENDDO
      ENDDO
C
C     FOURTH BATCH
      F1 = RKFAC1*FACAB1
      F2 = RKFAC1*FACAB2
      DO L=1,NFUND
        DO K=1,NFUNC
        M = (K-1)*NFUND + L
        BMAT(IB1+J,JD1+L) = BMAT(IB1+J,JD1+L)
     &                              - F1*RR(M,13)*DENT(IC1+K,JA1+I)
     &                              - F2*RR(M,5 )*DENT(IC1+K,JA2+I)
     &                              - F1*RR(M,15)*DENT(IC2+K,JA1+I)
     &                              - F2*RR(M,7 )*DENT(IC2+K,JA2+I)
        BMAT(IB1+J,JD2+L) = BMAT(IB1+J,JD2+L)
     &                              - F1*RR(M,14)*DENT(IC1+K,JA1+I)
     &                              - F2*RR(M,6 )*DENT(IC1+K,JA2+I)
     &                              - F1*RR(M,16)*DENT(IC2+K,JA1+I)
     &                              - F2*RR(M,8 )*DENT(IC2+K,JA2+I)
        BMAT(IB2+J,JD1+L) = BMAT(IB2+J,JD1+L)
     &                              - F2*RR(M,9 )*DENT(IC1+K,JA1+I)
     &                              - F1*RR(M,1 )*DENT(IC1+K,JA2+I)
     &                              - F2*RR(M,11)*DENT(IC2+K,JA1+I)
     &                              - F1*RR(M,3 )*DENT(IC2+K,JA2+I)
        BMAT(IB2+J,JD2+L) = BMAT(IB2+J,JD2+L)
     &                              - F2*RR(M,10)*DENT(IC1+K,JA1+I)
     &                              - F1*RR(M,2 )*DENT(IC1+K,JA2+I)
     &                              - F2*RR(M,12)*DENT(IC2+K,JA1+I)
     &                              - F1*RR(M,4 )*DENT(IC2+K,JA2+I)
        ENDDO
      ENDDO
C
C *** BEGIN CONDITIONAL FOR BLOCK COMBINATIONS
      IF (IQ12.NE.IQ34) THEN
C
C       FIRST BATCH
        F1 = RKFAC1*FACAB1
        F2 = RKFAC1*FACAB2
        M  = 0
        DO K=1,NFUNC
          DO L=1,NFUND
            M = M+1
            BMAT(IC1+K,JD1+L) = BMAT(IC1+K,JD1+L)
     &           + RR(M,1 )*(DENT(IA1+I,JB1+J) + F1*DENT(JB2+J,IA2+I))
     &           + RR(M,13)*(DENT(IA2+I,JB2+J) + F1*DENT(JB1+J,IA1+I))       
     &           + RR(M,5 )*(DENT(IA1+I,JB2+J) + F2*DENT(JB1+J,IA2+I))
     &           + RR(M,9 )*(DENT(IA2+I,JB1+J) + F2*DENT(JB2+J,IA1+I))
C
            BMAT(IC1+K,JD2+L) = BMAT(IC1+K,JD2+L)
     &           + RR(M,2 )*(DENT(IA1+I,JB1+J) + F1*DENT(JB2+J,IA2+I))
     &           + RR(M,14)*(DENT(IA2+I,JB2+J) + F1*DENT(JB1+J,IA1+I))
     &           + RR(M,6 )*(DENT(IA1+I,JB2+J) + F2*DENT(JB1+J,IA2+I))
     &           + RR(M,10)*(DENT(IA2+I,JB1+J) + F2*DENT(JB2+J,IA1+I))
C
            BMAT(IC2+K,JD1+L) = BMAT(IC2+K,JD1+L)
     &           + RR(M,3 )*(DENT(IA1+I,JB1+J) + F1*DENT(JB2+J,IA2+I))
     &           + RR(M,15)*(DENT(IA2+I,JB2+J) + F1*DENT(JB1+J,IA1+I))
     &           + RR(M,7 )*(DENT(IA1+I,JB2+J) + F2*DENT(JB1+J,IA2+I))
     &           + RR(M,11)*(DENT(IA2+I,JB1+J) + F2*DENT(JB2+J,IA1+I))
C
            BMAT(IC2+K,JD2+L) = BMAT(IC2+K,JD2+L)
     &           + RR(M,4 )*(DENT(IA1+I,JB1+J) + F1*DENT(JB2+J,IA2+I))
     &           + RR(M,16)*(DENT(IA2+I,JB2+J) + F1*DENT(JB1+J,IA1+I))
     &           + RR(M,8 )*(DENT(IA1+I,JB2+J) + F2*DENT(JB1+J,IA2+I))
     &           + RR(M,12)*(DENT(IA2+I,JB1+J) + F2*DENT(JB2+J,IA1+I))
          ENDDO
        ENDDO
C
C       SECOND BATCH
        F1 = RKFAC1*FACAB1
        F2 = RKFAC1*FACAB2
        M  = 0
        DO K=1,NFUNC
          DO L=1,NFUND
            M = M+1
C
            BMAT(IC1+K,JA1+I) = BMAT(IC1+K,JA1+I)
     &                                - F1*RR(M,13)*DENT(JB1+J,ID1+L)
     &                                - F1*RR(M,14)*DENT(JB1+J,ID2+L)
     &                                - F2*RR(M,9 )*DENT(JB2+J,ID1+L)
     &                                - F2*RR(M,10)*DENT(JB2+J,ID2+L)
C
            BMAT(IC1+K,JA2+I) = BMAT(IC1+K,JA2+I)
     &                                - F2*RR(M,5 )*DENT(JB1+J,ID1+L)
     &                                - F2*RR(M,6 )*DENT(JB1+J,ID2+L)
     &                                - F1*RR(M,1 )*DENT(JB2+J,ID1+L)
     &                                - F1*RR(M,2 )*DENT(JB2+J,ID2+L)
C
            BMAT(IC2+K,JA1+I) = BMAT(IC2+K,JA1+I)
     &                                - F1*RR(M,15)*DENT(JB1+J,ID1+L)
     &                                - F1*RR(M,16)*DENT(JB1+J,ID2+L)
     &                                - F2*RR(M,11)*DENT(JB2+J,ID1+L)
     &                                - F2*RR(M,12)*DENT(JB2+J,ID2+L)
C
            BMAT(IC2+K,JA2+I) = BMAT(IC2+K,JA2+I)
     &                                - F2*RR(M,7 )*DENT(JB1+J,ID1+L)
     &                                - F2*RR(M,8 )*DENT(JB1+J,ID2+L)
     &                                - F1*RR(M,3 )*DENT(JB2+J,ID1+L)
     &                                - F1*RR(M,4 )*DENT(JB2+J,ID2+L)
        ENDDO
      ENDDO
C
C       THIRD BATCH
        F1 = RKFAC2*FACCD1
        F2 = RKFAC2*FACCD2
        DO L=1,NFUND
          DO K=1,NFUNC
            M = (K-1)*NFUND + L
C
            BMAT(ID1+L,JB1+J) = BMAT(ID1+L,JB1+J)
     &                                - F1*RR(M,4 )*DENT(JA1+I,IC1+K)
     &                                - F2*RR(M,2 )*DENT(JA1+I,IC2+K)
     &                                - F1*RR(M,12)*DENT(JA2+I,IC1+K)
     &                                - F2*RR(M,10)*DENT(JA2+I,IC2+K)
C
            BMAT(ID1+L,JB2+J) = BMAT(ID1+L,JB2+J)
     &                                - F1*RR(M,8 )*DENT(JA1+I,IC1+K)
     &                                - F2*RR(M,6 )*DENT(JA1+I,IC2+K)
     &                                - F1*RR(M,16)*DENT(JA2+I,IC1+K)
     &                                - F2*RR(M,14)*DENT(JA2+I,IC2+K)
C
            BMAT(ID2+L,JB1+J) = BMAT(ID2+L,JB1+J)
     &                                - F2*RR(M,3 )*DENT(JA1+I,IC1+K)
     &                                - F1*RR(M,1 )*DENT(JA1+I,IC2+K)
     &                                - F2*RR(M,11)*DENT(JA2+I,IC1+K)
     &                                - F1*RR(M,9 )*DENT(JA2+I,IC2+K)
C
            BMAT(ID2+L,JB2+J) = BMAT(ID2+L,JB2+J)
     &                                - F2*RR(M,7 )*DENT(JA1+I,IC1+K)
     &                                - F1*RR(M,5 )*DENT(JA1+I,IC2+K)
     &                                - F2*RR(M,15)*DENT(JA2+I,IC1+K)
     &                                - F1*RR(M,13)*DENT(JA2+I,IC2+K)
          ENDDO
        ENDDO
C
C       FOURTH BATCH
        M = 0
        DO K=1,NFUNC
          DO L=1,NFUND
            M = M+1
C
            BMAT(IC1+K,JB1+J) = BMAT(IC1+K,JB1+J)
     &                                   - RR(M,1 )*DENT(JA1+I,ID1+L)
     &                                   - RR(M,2 )*DENT(JA1+I,ID2+L)
     &                                   - RR(M,9 )*DENT(JA2+I,ID1+L)
     &                                   - RR(M,10)*DENT(JA2+I,ID2+L)
C
            BMAT(IC1+K,JB2+J) = BMAT(IC1+K,JB2+J)
     &                                   - RR(M,5 )*DENT(JA1+I,ID1+L)
     &                                   - RR(M,6 )*DENT(JA1+I,ID2+L)
     &                                   - RR(M,13)*DENT(JA2+I,ID1+L)
     &                                   - RR(M,14)*DENT(JA2+I,ID2+L)
C
            BMAT(IC2+K,JB1+J) = BMAT(IC2+K,JB1+J)
     &                                   - RR(M,3 )*DENT(JA1+I,ID1+L)
     &                                   - RR(M,4 )*DENT(JA1+I,ID2+L)
     &                                   - RR(M,11)*DENT(JA2+I,ID1+L)
     &                                   - RR(M,12)*DENT(JA2+I,ID2+L)
C
            BMAT(IC2+K,JB2+J) = BMAT(IC2+K,JB2+J)
     &                                   - RR(M,7 )*DENT(JA1+I,ID1+L)
     &                                   - RR(M,8 )*DENT(JA1+I,ID2+L)
     &                                   - RR(M,15)*DENT(JA2+I,ID1+L)
     &                                   - RR(M,16)*DENT(JA2+I,ID2+L)
          ENDDO
        ENDDO
C
C *** END CONDITIONAL OVER BLOCK COMBINATIONS
      ENDIF
C
C*********************************************************************C
C     END OF BREIT MATRIX CONSTRUCTION                                C
C*********************************************************************C      
C
3000  CONTINUE
1999  CONTINUE
1000  CONTINUE
2000  CONTINUE
C
C*********************************************************************C
C     COMPLETE CONSTRUCTION OF BREIT MATRIX BY MATRIX CONJUGATION.    C
C*********************************************************************C
C
      DO I=1,NSHIFT
        DO J=NSHIFT+1,NDIM
          BMAT(J,I) = DCONJG(BMAT(I,J))
        ENDDO
      ENDDO
C
C     INITIALISE TEMPORARY COUNTER FOR BREIT ENERGY
      ETMP = 0.0D0
C    
C     USE BMAT AND DENSITY MATRIX TO CALCULATE BREIT ENERGY
      DO I=1,NDIM
        DO J=1,NDIM
C         DFNOTE: WHY ONLY THE REAL COMPONENT OF THE DENSITY MATRIX?        
          ETMP = ETMP + 0.5D0*DREAL(DENT(I,J)*BMAT(I,J))
        ENDDO
      ENDDO
      EBRT = ETMP
C
      RETURN
      END
C
C
      SUBROUTINE BII(BINTG,XYZ,KQN,MQN,EXPT,NFUNS,IBAS,JBAS,
     &                                             IEMAKE,INUCAB,INUCCD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C 
C                                                                     C
C                         BBBBBBB  IIII IIII                          C
C                         BB    BB  II   II                           C
C                         BB    BB  II   II                           C
C                         BBBBBBB   II   II                           C
C                         BB    BB  II   II                           C
C                         BB    BB  II   II                           C
C                         BBBBBBB  IIII IIII                          C
C                                                                     C
C ------------------------------------------------------------------- C
C     BII CONSTRUCTS INTERMEDIATE MATRICES FOR USE IN THE             C
C     MCMURCHIE-DAVIDSON EMPLOYMENT OF THE BREIT INTERACTION          C
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
C     I,J         INDEX FOR BASIS FUNCTION PAIR ON AB                 C
C     IEMAKE      IEMAKE = 0 DON'T RECALCULATE E-COEFFICIENTS         C
C                 IEMAKE = 1 DO    RECALCULATE E-COEFFICIENTS         C
C ------------------------------------------------------------------- C
C     TRY TO USE INUC COMBINATIONS TO MAKE ATOMIC CALCULATIONS EASIER C
C ------------------------------------------------------------------- C
C     DFNOTE: TESTING PROCEDURE COULD BE MADE MORE EFFICIENT.         C
C             APPLY TEST FOR REAL AND IMAG. COMPONENTS LIKE OLD VRSN  C
C*********************************************************************C
      PARAMETER (MXDM=1600,NCENTM=6,MAXB=26,MAXB2=MAXB*MAXB,NKAPM=9,
     &           LMAX=(NKAPM-1)/2,LAMMX=2*LMAX,IL4=2*LAMMX,MAXMV=NKAPM,
     &           MAXMV2=2*MAXMV,MLL=((NKAPM+8)*(NKAPM+9)*(NKAPM+10))/6,
     &           MAXR=969)
C
      COMPLEX*16 CONE
      COMPLEX*16 BINTG(MAXB2,16),Q1(MAXB2),Q2(MAXB2)
      COMPLEX*16 E11X(MAXB2,MLL),E11Y(MAXB2,MLL),E11Z(MAXB2,MLL),
     &           E21X(MAXB2,MLL),E21Y(MAXB2,MLL),E21Z(MAXB2,MLL)
      COMPLEX*16 EAB11(MAXB2,MLL,3),EAB21(MAXB2,MLL,3),
     &           ECD11(MAXB2,MLL,3),ECD21(MAXB2,MLL,3)
      COMPLEX*16 GAB11(MAXB2,MLL),GAB21(MAXB2,MLL)
C
      COMMON/IBSC/IBSCR(MAXB2),IBMAP(MAXB2)
      COMMON/ACCSS/INABCD(0:IL4,0:IL4,0:IL4),
     &             IVEC(MLL),JVEC(MLL),KVEC(MLL),LAMVEC(MLL)
C
      DIMENSION KQN(4),LQN(4),MQN(4),NFUNS(4)
      DIMENSION IAB11(MLL,3),IAB21(MLL,3),ICD11(MLL,3),ICD21(MLL,3)
      DIMENSION XYZ(3,4),RC(MAXB2,MAXR),PQ(MAXB2,3),EXPT(MAXB,4),
     &          APH(MAXB2),PREFAC(MAXB2),IRC(MAXR),ICART(3),JCART(3)
    
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
C     CALCULATE LAMBDA VALUES FOR EACH OVERLAP AND NUMBER OF INDICES
      LAMAB  = LQN(1) + LQN(2) + 1
      LAMCD  = LQN(3) + LQN(4) + 1
C
      NTUVAB = ((LAMAB+1)*(LAMAB+2)*(LAMAB+3))/6
      NTUVCD = ((LAMCD+1)*(LAMCD+2)*(LAMCD+3))/6
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB  = NFUNS(1)*NFUNS(2)
      MAXCD  = NFUNS(3)*NFUNS(4)
C
C*********************************************************************C
C     IF ASKED TO RECALCULATE E COEFFICIENTS, DO THIS FIRST           C
C*********************************************************************C
C
C *** BEGIN CONDITIONAL OVER IEMAKE INDEX
      IF(IEMAKE.EQ.1) THEN
C
C       START WITH AB PAIRS
C
C       CALCULATE NORMALISATION CONSTANTS
        CALL RNORMF(EXPT,LQN,NFUNS,1,2)
C
C       GENERATE ELS(AB) COEFFICIENTS
        CALL EMAKELS(E11X,E21X,EXPT,KQN,MQN,NFUNS,IALTAB,1,2,1)
        CALL EMAKELS(E11Y,E21Y,EXPT,KQN,MQN,NFUNS,IALTAB,1,2,2)
        CALL EMAKELS(E11Z,E21Z,EXPT,KQN,MQN,NFUNS,IALTAB,1,2,3)
C
C
C       COPY TEMPORARY VALUES ACROSS TO A TOTAL ARRAY
        DO M=1,MAXAB
          DO IAB=1,NTUVAB
            EAB11(M,IAB,1) = E11X(M,IAB)
            EAB21(M,IAB,1) = E21X(M,IAB)
            EAB11(M,IAB,2) = E11Y(M,IAB)
            EAB21(M,IAB,2) = E21Y(M,IAB)
            EAB11(M,IAB,3) = E11Z(M,IAB)
            EAB21(M,IAB,3) = E21Z(M,IAB)
          ENDDO
        ENDDO
C
C       SCREENING PROCEDURE: TEST MAGNITUDES OF E-COEFFICIENT LISTS
        DO ICP=1,3
          DO IAB=1,NTUVAB
C
C           11 OVERLAP (AB PAIRS)
            TEST = CDASUM(MAXAB,EAB11(1,IAB,ICP))
            IF(TEST.LE.SENS) THEN
              IAB11(IAB,ICP) = 0
            ELSE
              IAB11(IAB,ICP) = 1
            ENDIF
C
C           21 OVERLAP (AB PAIRS)
            TEST = CDASUM(MAXAB,EAB21(1,IAB,ICP))
            IF(TEST.LE.SENS) THEN
              IAB21(IAB,ICP) = 0
            ELSE
              IAB21(IAB,ICP) = 1
            ENDIF
C
          ENDDO
        ENDDO
C
C       MOVE ON TO CD PAIRS
C
C       CALCULATE NORMALISATION CONSTANTS
        CALL RNORMF(EXPT,LQN,NFUNS,3,4)
C
C       GENERATE ELS(CD) COEFFICIENTS
        CALL EMAKELS(E11X,E21X,EXPT,KQN,MQN,NFUNS,IALTCD,3,4,1)
        CALL EMAKELS(E11Y,E21Y,EXPT,KQN,MQN,NFUNS,IALTCD,3,4,2)
        CALL EMAKELS(E11Z,E21Z,EXPT,KQN,MQN,NFUNS,IALTCD,3,4,3)
C
C       COPY TEMPORARY VALUES ACROSS TO A TOTAL ARRAY
        DO M=1,MAXCD
          DO ICD=1,NTUVCD
            ECD11(M,ICD,1) = E11X(M,ICD)
            ECD21(M,ICD,1) = E21X(M,ICD)
            ECD11(M,ICD,2) = E11Y(M,ICD)
            ECD21(M,ICD,2) = E21Y(M,ICD)
            ECD11(M,ICD,3) = E11Z(M,ICD)
            ECD21(M,ICD,3) = E21Z(M,ICD)
          ENDDO
        ENDDO
C
C       SCREENING PROCEDURE: TEST MAGNITUDES OF E-COEFFICIENT LISTS
        DO ICP=1,3
          DO ICD=1,NTUVCD
C
C           11 OVERLAP (CD PAIRS)
            TEST = CDASUM(MAXCD,ECD11(1,ICD,ICP))
            IF(TEST.LE.SENS) THEN
              ICD11(ICD,ICP) = 0
            ELSE
              ICD11(ICD,ICP) = 1
            ENDIF
C
C           21 OVERLAP (CD PAIRS)
            TEST = CDASUM(MAXCD,ECD21(1,ICD,ICP))
            IF(TEST.LE.SENS) THEN
              ICD21(ICD,ICP) = 0
            ELSE
              ICD21(ICD,ICP) = 1
            ENDIF
C
        ENDDO
      ENDDO
C
C *** END CONDITIONAL OVER IEMAKE INDEX
      ENDIF
C
C*********************************************************************C
C     R-INTEGRAL EVALUATION                                           C
C*********************************************************************C
C
C     GAUSSIAN OVERLAP VALUES
      EIJ = EXPT(IBAS,1) + EXPT(JBAS,2)
      PX  = (XYZ(1,1)*EXPT(IBAS,1) + XYZ(1,2)*EXPT(JBAS,2))/EIJ
      PY  = (XYZ(2,1)*EXPT(IBAS,1) + XYZ(2,2)*EXPT(JBAS,2))/EIJ
      PZ  = (XYZ(3,1)*EXPT(IBAS,1) + XYZ(3,2)*EXPT(JBAS,2))/EIJ
C
      M = 0
      MPRIME = 0
      DO K=1,NFUNS(3)
        DO L=1,NFUNS(4)
          M = M+1
          IF(IBSCR(M).EQ.1) THEN
            MPRIME = MPRIME+1
            EKL = EXPT(K,3) + EXPT(L,4)
            QX  = (XYZ(1,3)*EXPT(K,3)+XYZ(1,4)*EXPT(L,4))/EKL
            QY  = (XYZ(2,3)*EXPT(K,3)+XYZ(2,4)*EXPT(L,4))/EKL
            QZ  = (XYZ(3,3)*EXPT(K,3)+XYZ(3,4)*EXPT(L,4))/EKL
            APH(MPRIME)    = (EIJ*EKL)/(EIJ+EKL)
            PQ(MPRIME,1)   = QX-PX
            PQ(MPRIME,2)   = QY-PY
            PQ(MPRIME,3)   = QZ-PZ
            PREFAC(MPRIME) = 2.0D0*(ROOTPI**5)/(DSQRT(EIJ+EKL)*EIJ*EKL)
          ENDIF
        ENDDO
      ENDDO
C
C     MAXM WAS ONCE CALLED MMXNEW -- IT COUNTS # BASIS FUNCTION
C     OVERLAPS GIVEN THE SCREENING TEST IN BREIT SUBROUTINE
      MMXNEW = MPRIME
      MAXM = M
C
      CALL RMAKE(RC,PQ,APH,MMXNEW,LAMAB+LAMCD+2)
C
C     INITIALIZE ARRAY TO IMPLEMENT SPARSENESS IN R-VECTOR
      LAMABCD  = LAMAB + LAMCD
      NTUVABCD = ((LAMABCD+1)*(LAMABCD+2)*(LAMABCD+3))/6
C
      DO MRC=1,NTUVABCD
        TEST = DASUM(MMXNEW,RC(1,MRC),1)
        IF(TEST.LE.SENS) THEN
          IRC(MRC) = 0
        ELSE
          IRC(MRC) = 1
        ENDIF
      ENDDO
C
C     INITIALISE BINTG ARRAY
      DO ITG=1,16
        DO M=1,MAXCD
          BINTG(M,ITG) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C*********************************************************************C
C     BEGIN FIRST LOOP: CARTESIAN INDEX {IX,IY,IZ} FOR CENTRE AB      C
C*********************************************************************C
C
      DO 6000 ICMP=1,3
C
C     INITIATE LOOP OVER ALL FINITE SUM INDICES {ALPH,BETA,GAMA} FOR AB
      DO IAB=1,NTUVAB
C
        IAB1 = 0
        IAB2 = 0
        IAB3 = 0
        IAB4 = 0
C
C       SCREENING MARKERS
        IF((IAB11(IAB,1)+IAB11(IAB,2)+IAB11(IAB,3)).NE.0) THEN
          IAB1 = 1
        ENDIF
C
        IF((IAB21(IAB,1)+IAB21(IAB,2)+IAB21(IAB,3)).NE.0) THEN
          IAB2 = 1
        ENDIF
C
C       INITIALISE IF ANY E-COEFF. FOR THIS {ALPH,BETA,GAMA}.GT.0.0D0
        IF((IAB1+IAB2).GT.0) THEN
          DO M=1,MMXNEW
            GAB11(M,IAB) = DCMPLX(0.0D0,0.0D0)
            GAB21(M,IAB) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDIF
C
C*********************************************************************C
C       FOR A GIVEN AB OVERLAP M, SUM UP ALL COMBINATIONS OVER CD     C
C       THAT CONTRIBUTE TO THE G-ARRAYS FOR INDEX {ALPH,BETA,GAMA}.   C
C*********************************************************************C
C
C       BEGIN LOOP OVER CD INDEX {JX,JY,JZ} FOR CENTRE CD
        DO JCMP=1,3
C
C         IF THE CARTESIAN INDICES ARE THE SAME {BXX, BYY, BZZ}
          IF(ICMP.EQ.JCMP) THEN
C
C           BEGIN LOOP OVER ALL (CD) ADDRESSES FOR A GIVEN (AB) ADDRESS
            DO ICD=1,NTUVCD
C
C             OVERALL ADDRESS OF THIS {ALPH,BETA,GAMA},{ALPH',BETA',GAMA'}
              IRABCD = INABCD(IVEC(IAB)+IVEC(ICD),
     &                        JVEC(IAB)+JVEC(ICD),
     &                        KVEC(IAB)+KVEC(ICD))
C
C             DON'T BOTHER CALCULATING IF THE BIGGEST R-INTEGRAL IS SMALL
              IF(IRC(IRABCD).EQ.0) GOTO 797
C
C             IF E11(CD) IS BIG ENOUGH, CONTINUE TO ADD ALL OVERLAPS
              IF(ICD11(ICD,JCMP).EQ.1) THEN
                DO M=1,MMXNEW
                  GAB11(M,IAB) = GAB11(M,IAB)
     &                       - ECD11(IBMAP(M),ICD,JCMP)*RC(M,IRABCD)
                ENDDO
              ENDIF
C
C             IF E21(CD) IS BIG ENOUGH, CONTINUE TO ADD ALL OVERLAPS
              IF(ICD21(ICD,JCMP).EQ.1) THEN
                DO M=1,MMXNEW
                  GAB21(M,IAB) = GAB21(M,IAB)
     &                       - ECD21(IBMAP(M),ICD,JCMP)*RC(M,IRABCD)
                ENDDO
              ENDIF
C
797         CONTINUE
            ENDDO
C
C         END CONDITIONAL OVER LIKE INDICES
          ENDIF
C
C         CARTESIAN VECTOR FOR IDX,JDX VALUES
          CALL NCART(ICART,ICMP)
          CALL NCART(JCART,JCMP)
C
C         LOOP OVER FINITE SUM INDICES {ALPH',BETA',GAMA'} FOR CD
          DO ICD=1,NTUVCD
C
            IF(JCMP.EQ.1) THEN
              RTP = DFLOAT(IVEC(IAB)+IVEC(ICD))
            ELSEIF(JCMP.EQ.2) THEN
              RTP = DFLOAT(JVEC(IAB)+JVEC(ICD))
            ELSE
              RTP = DFLOAT(KVEC(IAB)+KVEC(ICD))
            ENDIF
C
C           STARTING ADDRESSES FOR REFERENCE IN A MOMENT
            I1 = IVEC(IAB) + IVEC(ICD) + ICART(1) + JCART(1)
            J1 = JVEC(IAB) + JVEC(ICD) + ICART(2) + JCART(2)
            K1 = KVEC(IAB) + KVEC(ICD) + ICART(3) + JCART(3)
            N1 = INABCD(I1,J1,K1)
C
            I1 = IVEC(IAB) + IVEC(ICD) + ICART(1)
            J1 = JVEC(IAB) + JVEC(ICD) + ICART(2)
            K1 = KVEC(IAB) + KVEC(ICD) + ICART(3)
            N2 = INABCD(I1,J1,K1)
C
            I1 = IVEC(IAB) + IVEC(ICD) + ICART(1) - JCART(1)
            J1 = JVEC(IAB) + JVEC(ICD) + ICART(2) - JCART(2)
            K1 = KVEC(IAB) + KVEC(ICD) + ICART(3) - JCART(3)
C
            IF((I1.GE.0).AND.(J1.GE.0).AND.(K1.GE.0)) THEN
              N3 = INABCD(I1,J1,K1)
            ELSE
              N3  = 1
              RTP = 0.0D0
            ENDIF
C
C           IF E11(CD) IS BIG ENOUGH, CONTINUE TO ADD ALL OVERLAPS
            IF(ICD11(ICD,JCMP).EQ.1) THEN
              DO M=1,MMXNEW
                T1 = RC(M,N1)*0.5D0/APH(M)
                T2 = RC(M,N2)*PQ(M,JCMP)
                T3 = RC(M,N3)*RTP
                GAB11(M,IAB) = GAB11(M,IAB)
     &                         + ECD11(IBMAP(M),ICD,JCMP)*(T1-T2+T3)
              ENDDO
            ENDIF
C
C           IF E21(CD) IS BIG ENOUGH, CONTINUE TO ADD ALL OVERLAPS
            IF(ICD21(ICD,JCMP).EQ.1) THEN
              DO M=1,MMXNEW
                T1 = RC(M,N1)*0.5D0/APH(M)
                T2 = RC(M,N2)*PQ(M,JCMP)
                T3 = RC(M,N3)*RTP
                GAB21(M,IAB) = GAB21(M,IAB)
     &                         + ECD21(IBMAP(M),ICD,JCMP)*(T1-T2+T3)
              ENDDO
            ENDIF
C
C         END LOOP OVER {ALPH',BETA',GAMA'}
          ENDDO
C
C       END LOOP OVER CD INDEX {JX,JY,JZ}
        ENDDO
C
C     END LOOP OVER FINITE SUM INDICES {ALPH,BETA,GAMA} FOR AB
      ENDDO
C
C*********************************************************************C
C     GENERATE ALL POSSIBLE TWO-ELECTRON INTEGRALS FROM THE           C
C     EAB COEFFICIENTS AND THE G-ARRAYS                               C
C*********************************************************************C
C
C     CALCULATE PHASE FACTORS FOR MQN AND KQN COMBINATIONS (P1,P2,P3)
      P1 = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
      P2 = DFLOAT((-1)**((MQN(3)-MQN(4))/2))
C
      P1 =-(P1*DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2))))
      P2 =-(P2*DFLOAT((KQN(3)*KQN(4))/IABS(KQN(3)*KQN(4))))
C
      P3 = P1*P2
C
C*********************************************************************C
C     INTEGRAL BATCH 1: ( - - || - - )                                C
C*********************************************************************C
C
      NTUVAB = ((LAMAB+1)*(LAMAB+2)*(LAMAB+3))/6
      IJ     = (IBAS-1)*NFUNS(2) + JBAS
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSENESS 
      DO M=1,MMXNEW
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {ALPH,BETA,GAMA}
      DO IAB=1,NTUVAB
        IF(IAB11(IAB,ICMP).EQ.1) THEN
          DO M=1,MMXNEW
            Q1(M) = Q1(M) +      EAB11(IJ,IAB,ICMP)*DREAL(GAB11(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB11(IJ,IAB,ICMP)*DIMAG(GAB11(M,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE BINTG ARRAY
      DO M=1,MMXNEW
        BINTG(M,1 ) = BINTG(M,1 )     +    (Q1(M)+Q2(M))*PREFAC(M)
        BINTG(M,4 ) = BINTG(M,4 )     + P2*(Q1(M)-Q2(M))*PREFAC(M)
        BINTG(M,13) = P3*DCONJG(BINTG(M,4))
        BINTG(M,16) = P3*DCONJG(BINTG(M,1))
      ENDDO
C
C*********************************************************************C
C     INTEGRAL BATCH 2: ( - - || + - )                                C
C*********************************************************************C
C
      NTUVAB = ((LAMAB+1)*(LAMAB+2)*(LAMAB+3))/6
      IJ     = (IBAS-1)*NFUNS(2)+JBAS
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSENESS 
      DO M=1,MMXNEW
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {ALPH,BETA,GAMA}
      DO IAB=1,NTUVAB
        IF(IAB11(IAB,ICMP).EQ.1) THEN
          DO M=1,MMXNEW
            Q1(M) = Q1(M) +      EAB11(IJ,IAB,ICMP)*DREAL(GAB21(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB11(IJ,IAB,ICMP)*DIMAG(GAB21(M,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE BINTG ARRAY
      DO M=1,MMXNEW
        BINTG(M,3 ) = BINTG(M,3 )     +    (Q1(M)+Q2(M))*PREFAC(M)
        BINTG(M,2 ) = BINTG(M,2 )     - P2*(Q1(M)-Q2(M))*PREFAC(M)
        BINTG(M,15) =-P3*DCONJG(BINTG(M,2))
        BINTG(M,14) =-P3*DCONJG(BINTG(M,3))
      ENDDO
C
C*********************************************************************C
C     INTEGRAL BATCH 3: ( + - || - - )                                C
C*********************************************************************C
C
      NTUVAB = ((LAMAB+1)*(LAMAB+2)*(LAMAB+3))/6
      IJ     = (IBAS-1)*NFUNS(2) + JBAS
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSENESS 
      DO M=1,MMXNEW
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {ALPH,BETA,GAMA}
      DO IAB=1,NTUVAB
        IF(IAB21(IAB,ICMP).EQ.1) THEN
          DO M=1,MMXNEW
            Q1(M) = Q1(M) +      EAB21(IJ,IAB,ICMP)*DREAL(GAB11(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB21(IJ,IAB,ICMP)*DIMAG(GAB11(M,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C
C     FILL THIS BATCH OF THE BINTG ARRAY
      DO M=1,MMXNEW
        BINTG(M,9 ) = BINTG(M,9 )     +    (Q1(M)+Q2(M))*PREFAC(M)
        BINTG(M,12) = BINTG(M,12)     + P2*(Q1(M)-Q2(M))*PREFAC(M)
        BINTG(M,5 ) =-P3*DCONJG(BINTG(M,12))
        BINTG(M,8 ) =-P3*DCONJG(BINTG(M,9 ))      
      ENDDO
C
C*********************************************************************C
C     INTEGRAL BATCH 4: ( + - || + - )                                C
C*********************************************************************C
C
      NTUVAB = ((LAMAB+1)*(LAMAB+2)*(LAMAB+3))/6
      IJ     = (IBAS-1)*NFUNS(2)+JBAS
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSENESS 
      DO M=1,MMXNEW
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {ALPH,BETA,GAMA}
      DO IAB=1,NTUVAB
        IF(IAB21(IAB,ICMP).EQ.1) THEN
          DO M=1,MMXNEW
            Q1(M) = Q1(M) +      EAB21(IJ,IAB,ICMP)*DREAL(GAB21(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB21(IJ,IAB,ICMP)*DIMAG(GAB21(M,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE BINTG ARRAY
      DO M=1,MMXNEW
        BINTG(M,11) = BINTG(M,11)      +    (Q1(M)+Q2(M))*PREFAC(M)
        BINTG(M,10) = BINTG(M,10)      - P2*(Q1(M)-Q2(M))*PREFAC(M)
        BINTG(M,7 ) =  P3*DCONJG(BINTG(M,10))
        BINTG(M,6 ) = P3*DCONJG(BINTG(M,11))      
      ENDDO
C
C     END LOOP OVER INDICES {IX,IY,IZ}
6000  CONTINUE
C
C*********************************************************************C
C     BINTG ARRAY NOW FULLY CONSTRUCTED                               C
C*********************************************************************C
C
C     INCLUDE THE OUTSIDE FACTOR OF (1/2)
      DO ITG=1,16
        DO MPRIME=1,MMXNEW
          BINTG(IBMAP(MPRIME),ITG) = 0.5D0*BINTG(IBMAP(MPRIME),ITG)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE NCART(IVECT,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C            NN    NN  CCCCCCC     AA    RRRRRRR  TTTTTTTT            C
C            NNN   NN CC     CC   AAAA   RR    RR    TT               C
C            NNNN  NN CC         AA  AA  RR    RR    TT               C
C            NN NN NN CC        AA    AA RR    RR    TT               C
C            NN  NNNN CC        AAAAAAAA RRRRRRR     TT               C
C            NN   NNN CC     CC AA    AA RR    RR    TT               C
C            NN    NN  CCCCCCC  AA    AA RR    RR    TT               C
C                                                                     C
C ------------------------------------------------------------------- C
C     NCART RETURNS THE CARTESIAN INDEX FROM THE INDEX VALUE IND.     C
C*********************************************************************C
C
      DIMENSION IVECT(3)
C
      IF(IND.EQ.1) THEN
        IVECT(1) = 1
        IVECT(2) = 0
        IVECT(3) = 0
      ELSEIF(IND.EQ.2) THEN
        IVECT(1) = 0
        IVECT(2) = 1
        IVECT(3) = 0
      ELSEIF(IND.EQ.3) THEN
        IVECT(1) = 0
        IVECT(2) = 0
        IVECT(3) = 1
      ELSE
        WRITE(*,*) 'In NCART: supplied index not valid',IND
        RETURN
      ENDIF
C
      RETURN
      END

