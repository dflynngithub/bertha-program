C     THIS IS A VERSION I KEPT BEFORE CHANGING ARRAY STYLES
C
      SUBROUTINE MBPT2(MINO,NUMO,MINV,NUMV,G2INT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            MM       MM BBBBBBB  PPPPPPP TTTTTTTT 222222              C
C            MMM     MMM BB    BB PP    PP   TT   22    22             C
C            MMMM   MMMM BB    BB PP    PP   TT         22             C
C            MM MM MM MM BBBBBBB  PP    PP   TT       22               C
C            MM  MMM  MM BB    BB PPPPPPP    TT     22                 C
C            MM   M   MM BB    BB PP         TT   22                   C
C            MM       MM BBBBBBB  PP         TT   22222222             C
C                                                                      C
C -------------------------------------------------------------------- C
C  MBPT2 CALCULATES SECOND-ORDER PAIR CORRELATION CORRECTIONS OVER A   C
C  USER-SPECIFIED TWO-BODY OPERATOR. IT USES A RELATIVISTIC ADAPTION   C
C  OF THE DIRECT MP2 ALGORITHM OF HEAD-GORDON, POPLE AND FRISCH (1988).C
C  THIS IS NOT MOLLER-PLESSET PERTURBATION THEORY, SO CALL IT MBPT2.   C
C  ELECTRON REPULATION INTEGRALS ARE GENERATED ONCE, AT THE EXPENSE OF C
C  ADDITIONAL MEMORY STORAGE -- THEY ARE CONTRACTED IN FOUR STEPS.     C
C  THIS ALGORITHM EXPLOITS SYMMETRIES BETWEEN (A,B) AND (C,D) IN THE   C
C  INTEGRALS (AB|CD), BUT NOT THE SWAP (AB|CD)<->(CD|AB). STRUCTURE    C
C  MIMICS THAT OF 'COULOMB', BUT DOES NOT IMPLEMENT SELECTION RULES.   C
C  CONTRACTION ADDRESSES ARE CALCULATED EXPLICITLY RATHER THAN JUST    C
C  WITH UPDATING COUNTERS -- THIS IS TO GUIDE THE USER IN TRACKING.    C
C -------------------------------------------------------------------- C
C INPUT:                                                               C
C ▶ MINO  - LOWEST OCCUPIED STATE TO ACCOUNT FOR. (FULL: 1)            C
C ▶ NUMO  - NUMBER OF OCCUPIED STATES TO ACCOUNT FOR. (FULL: NOCC)     C
C ▶ MINV  - LOWEST VIRTUAL STATE TO ACCOUNT FOR. (FULL: NOCC+1)        C
C ▶ NUMV  - NUMBER OF VIRTUAL STATES TO ACCOUNT FOR. (FULL: NVRT)      C
C ▶ G2INT - NAME OF TWO-BODY OPERATOR ('COULM' OR 'BREIT').            C
C -------------------------------------------------------------------- C
C ▶ LOOP STRUCTURE NOT LIKE COULOMB -- MUST CONTRACT C AND D ASAP.     C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL BGAUNT
C
      CHARACTER*5 G2INT
      CHARACTER*16 HMS
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4),LQN(4),JQN(4)
      DIMENSION INDEX(MCT,-(MEL+1):MEL,MKP),ITN(2)
      DIMENSION EAB2(NUMO,NUMO,6),EA2(NUMO,6)
C
      INTEGER*8 INEST,NNEST
C
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 B(MBS*NUMO,8),SB(MB2,NUMO*NUMV,4)
      COMPLEX*16 ASB1(MBS,NUMO*NUMO*NUMV,2),ASB2(MBS,NUMO*NUMO*NUMV,2)
      COMPLEX*16 RASB((NUMO+1)*NUMO*NUMV*NUMV/2)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/E0LL/E0LLFL(MFL,4),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,4),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EILS/EILSFL(MFL,12),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EISL/EISLFL(MFL,12),IADISL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS,IABSL,ICDSL
      COMMON/IRCM/IEAB,IECD,IRIJ(MBS,MBS)
      COMMON/ISCR/IMTX(MB2,11),ISCR(MB2),IMAP(MB2),IBCH,ITOG,MAXN
C
C     TURN OFF RC(AB|CD) LOCAL FILE PROCESS
      RCFILE = .FALSE.
C
C     INTEGRAL SKIPPING ON MOLECULAR GROUP SYMMETRY CLASS BASIS
      IF(SHAPE.EQ.'ATOMIC'.OR.SHAPE.EQ.'DIATOMIC'.OR.
     &                        SHAPE.EQ.'LINEAR') THEN
        ISYM = 1
      ELSE
        ISYM = 0
      ENDIF
C
C     WARNINGS BASED ON INVALID HMLT VS. G2INT COMBINATIONS
      IF(G2INT.EQ.'COULM') THEN
        IF(HMLT.EQ.'BARE') THEN
          WRITE(6, *) 'In MBPT2: HMLT = BARE but G2INT = COULM.'
          WRITE(7, *) 'In MBPT2: HMLT = BARE but G2INT = COULM.'
        ENDIF
      ELSEIF(G2INT.EQ.'BREIT') THEN
        IF(HMLT.EQ.'NORL') THEN
          WRITE(6, *) 'In MBPT2: HMLT = NORL but G2INT = BREIT.'
          WRITE(7, *) 'In MBPT2: HMLT = NORL but G2INT = BREIT.'
          RETURN
        ELSEIF(HMLT.EQ.'BARE') THEN
          WRITE(6, *) 'In MBPT2: HMLT = BARE but G2INT = BREIT.'
          WRITE(7, *) 'In MBPT2: HMLT = BARE but G2INT = BREIT.'
        ELSEIF(HMLT.EQ.'DHFR') THEN
          WRITE(6, *) 'In MBPT2: HMLT = DHFR but G2INT = BREIT.'
          WRITE(7, *) 'In MBPT2: HMLT = DHFR but G2INT = BREIT.'
        ELSEIF(HMLT.EQ.'DHFP') THEN
          WRITE(6, *) 'In MBPT2: HMLT = DHFP but G2INT = BREIT.'
          WRITE(7, *) 'In MBPT2: HMLT = DHFP but G2INT = BREIT.'
        ENDIF
      ENDIF
C
C     COMPONENT OVERLAP LABELS TO LOOP OVER
      IF(HMLT.EQ.'BARE') THEN
        RETURN
      ELSEIF(HMLT.EQ.'NORL') THEN
        IF(G2INT.EQ.'COULM') THEN
          ITSTRT = 1
          ITSTOP = 1
          ITSKIP = 1
        ENDIF
      ELSE
        IF(G2INT.EQ.'COULM') THEN
          ITSTRT = 4
          ITSTOP = 1
          ITSKIP =-3
        ELSEIF(G2INT.EQ.'BREIT') THEN
          ITSTRT = 2
          ITSTOP = 3
          ITSKIP = 1        
        ENDIF
      ENDIF
C
C     INITIALISE TIME COUNTERS
      TERI = 0.0D0
      TCN1 = 0.0D0
      TCN2 = 0.0D0
      TCN3 = 0.0D0
      TCN4 = 0.0D0
      TSUM = 0.0D0
C
      CALL CPU_TIME(TBEG)
C
C     IMPORT MBPT1 PAIR RESULTS FOR E1(ab)
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_MBPT1.dat',STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO IOCCA=1,NUMO
        DO IOCCB=1,NUMO
          READ(8, *) Q1,Q2,Q3,(EAB2(IOCCA,IOCCB,N),N=1,3)
        ENDDO
      ENDDO
      CLOSE(UNIT=8)
C
C     CLEAR THE REMAINDER OF THE CORRELATION PAIR ENERGY VALUES
      DO IOCCA=1,NUMO
        DO IOCCB=1,IOCCA
          DO N=4,6
            EAB2(IOCCA,IOCCB,N) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     CLEAR THE ARRAY FOR (AR|BS) VALUES
C      M = 0
C      DO IVRTR=1,NUMV
C        DO IOCCA=1,NUMO
C          DO IOCCB=1,IOCCA
C            DO IVRTS=1,NUMV
C              M = M+1
C              RASB(M) = DCMPLX(0.0D0,0.0D0)
C            ENDDO
C          ENDDO
C        ENDDO
C      ENDDO
       NRASB = (NUMO+1)*NUMO/2
       NRASB = NRASB*NUMV*NUMV
       WRITE(*,*) NUMO,NUMV,(NUMO+1)*NUMO/2,NUMV*NUMV,NRASB
       DO M=1,NRASB
         RASB(M) = DCMPLX(0.0D0,0.0D0)
       ENDDO
C
C**********************************************************************C
C     ORDERED INDEX OF (ICNT,KQN,MQN) COMBINATIONS                     C
C**********************************************************************C
C
      ICOUNT = 0
C
C     LOOP OVER NUCLEAR CENTRES
      DO ICT=1,NCNT
C
C       LOOP OVER KAPPA VALUES FOR THIS NUCLEAR CENTRE
        DO KN=1,NKAP(ICT)
C
C         IMPORT KAPPA, MAXIMUM MQN
          KAPPA = KAPA(KN,ICT)
          MJMAX = 2*IABS(KAPPA)-1
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
C     ORDERED COUNTER FOR NESTED LOOP PROGRESS                         C
C**********************************************************************C
C
      NNEST = 0
C
C     COMPONENT IT1, BLOCKS A AND B
      DO IT1=ITSTRT,ITSTOP,ITSKIP
      DO ICNTA=1,NCNT
      DO KA=1,NKAP(ICNTA)
      DO MA=1,IABS(KAPA(KA,ICNTA))
      DO ICNTB=1,ICNTA
      DO KB=1,NKAP(ICNTB)
      DO MB=1,IABS(KAPA(KB,ICNTB))
C
        NBAS(1) = NFNC(LVAL(KAPA(KA,ICNTA)),ICNTA)
        NBAS(2) = NFNC(LVAL(KAPA(KB,ICNTB)),ICNTB)
C
C       BASIS FUNCTIONS IN A AND B BLOCKS
        DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
C
C         COMPONENT IT2, BLOCK C
          DO IT2=ITSTRT,ITSTOP,ITSKIP
          DO ICNTC=1,NCNT
          DO KC=1,NKAP(ICNTC)
          DO MC=1,IABS(KAPA(KC,ICNTC))
C
C           BLOCK D
            DO ICNTD=1,NCNT
            DO KD=1,NKAP(ICNTD)
            DO MD=1,IABS(KAPA(KD,ICNTD))
C
              NBAS(3) = NFNC(LVAL(KAPA(KC,ICNTC)),ICNTC)
              NBAS(4) = NFNC(LVAL(KAPA(KD,ICNTD)),ICNTD)
C
C             FIRST CONTRACTION
              NNEST = NNEST + NBAS(4)*NUMO*NBAS(3)
C
            ENDDO
            ENDDO
            ENDDO
C
C           SECOND CONTRACTION
            NNEST = NNEST + NUMO*NUMV*NBAS(3)
C
          ENDDO
          ENDDO
          ENDDO
          ENDDO
C
        ENDDO
        ENDDO
C
C       THIRD CONTRACTION
        NNEST = NNEST + NUMO*(NUMO+1)*NUMV*NBAS(1)*NBAS(2)/2 
C
C       FOURTH CONTRACTION
        NNEST = NNEST + NUMO*NUMV*(NUMO+1)*NUMV*NBAS(1)/2 
C
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      
      WRITE(*,*) NNEST
C
C     START CONTRACTION COUNTER
      INEST = 0
      
      WRITE(*,*) 'START'
C
C**********************************************************************C
C     LOOP OVER COMPONENT OVERLAP LABELS TT(AB) (USE INDEX 1000)       C
C**********************************************************************C
C
C     COMPONENT LABEL FOR A AND B: TT = LL(1) or SS(4) <- COULOMB
C                                  TT = LS(2) or SL(3) <- BREIT
      DO 1000 IT1=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITN(1) = IT1
C
C       CALCULATE STARTING ADDRESS
        IF(ITN(1).EQ.1) THEN
          NADDA = 0
          NADDB = 0
        ELSEIF(ITN(1).EQ.2) THEN
          NADDA = 0
          NADDB = NSKP
        ELSEIF(ITN(1).EQ.3) THEN
          NADDA = NSKP
          NADDB = 0
        ELSEIF(ITN(1).EQ.4) THEN
          NADDA = NSKP
          NADDB = NSKP
        ENDIF
C
C**********************************************************************C
C     LOOP OVER BLOCKS OF A (USE INDEX 2000)                           C
C**********************************************************************C
C
C     LOOP OVER CENTRE A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,ICNTA)
        IF(KQN(1).LT.0) THEN
          LQN(1) =-KQN(1)-1
        ELSE
          LQN(1) = KQN(1)
        ENDIF
        JQN(1) = 2*IABS(KQN(1))-1
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),ICNTA)
        DO IBAS=1,NBAS(1)
          EXL(IBAS,1) = BEXL(IBAS,LQN(1),ICNTA)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 2000 MA=1,IABS(KQN(1))
        MQN(1) = 2*MA-1
C
C     INDEX ASSIGNMENT
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
C
C     STARTING FOCK ADDRESS FOR BASIS FUNCTIONS ON A
      NA1 = LRGE(ICNTA,KA,2*MA-1) + NADDA
      NA2 = LRGE(ICNTA,KA,2*MA  ) + NADDA
C
C**********************************************************************C
C     LOOP OVER BLOCKS OF B (USE INDEX 3000)                           C
C**********************************************************************C
C
C     LOOP OVER CENTRE B
      DO 3000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C     LOOP OVER KQN(B) VALUES
      DO 3000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KAPA(KB,ICNTB)
        IF(KQN(2).LT.0) THEN
          LQN(2) =-KQN(2)-1
        ELSE
          LQN(2) = KQN(2)
        ENDIF
        JQN(2) = 2*IABS(KQN(2))-1
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
        ENDDO
C
C     LOOP OVER |MQN(B)| VALUES
      DO 3000 MB=1,IABS(KQN(2))
        MQN(2) = 2*MB-1
C
C     INDEX ASSIGNMENT
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
C
C     STARTING FOCK ADDRESS FOR BASIS FUNCTIONS ON B
      NB1 = LRGE(ICNTB,KB,2*MB-1) + NADDB
      NB2 = LRGE(ICNTB,KB,2*MB  ) + NADDB
C
      WRITE(*,*) ' >',IT1,ICNTA,KA,MA,ICNTB,KB,MB
C
C     EQ-COEFFICIENT STARTING ADDRESSES FOR (AB) PAIR
      IF(IT1.EQ.1) THEN
        IABLL = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB)
      ELSEIF(IT1.EQ.2) THEN
        IABLS = IADILS(ICNTA,ICNTB,KA,KB,MA,MB)
      ELSEIF(IT1.EQ.3) THEN
        IABSL = IADISL(ICNTA,ICNTB,KA,KB,MA,MB)
      ELSEIF(IT1.EQ.4) THEN
        IABSS = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB)
      ENDIF
C
C     FLAG READ-IN OF E0(AB) COEFFICIENTS FOR THIS COMPONENT LABEL
      IEAB = 1
C
C     EQ-COEFFICIENT PHASE FACTORS FOR PERMUTATION OF R-INTEGRALS
      IF(G2INT.EQ.'COULM') THEN
        PAB1 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
        PAB2 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      ELSEIF(G2INT.EQ.'BREIT') THEN
        PAB1 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
        PAB2 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      ENDIF
C
C     SKIP CONTRIBUTIONS THAT ARISE BY PERMUTATION OF INTEGRALS
      IF(G2INT.EQ.'COULM') THEN
        WRITE(*,*) 'GOTCHA'
        IF(IQ1.LT.IQ2) GOTO 3001
      ENDIF
      WRITE(*,*) 'DONT GOTCHA'
C
C     CLEAR ARRAY FOR SECOND CONTRACTION (KL->SB)
      DO MIJ=1,NBAS(1)*NBAS(2)
        DO MSB=1,NUMO*NUMV
          DO IJSPIN=1,4
            SB(MIJ,MSB,IJSPIN) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
        WRITE(*,*) MIJ,NBAS(1)*NBAS(2)
      ENDDO
      
      WRITE(*,*) '>>',IT1,ICNTA,KA,MA,ICNTB,KB,MB
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS (IBAS,JBAS) (INDEX 4000)               C
C**********************************************************************C
C
      DO 4000 IBAS=1,NBAS(1)
      DO 4000 JBAS=1,NBAS(2)
C
C     LIST ADDRESS FOR THIS IBAS,JBAS COMBINATION
      MIJ = (IBAS-1)*NBAS(2)+JBAS
C
C**********************************************************************C
C     LOOP OVER COMPONENT OVERLAP LABELS TT(CD) (USE INDEX 5000)       C
C**********************************************************************C
C
C     COMPONENT LABEL FOR C AND D: TT = LL(1) or SS(4) <- COULOMB
C                                  TT = LS(2) or SL(3) <- BREIT
      DO 5000 IT2=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITN(2) = IT2
C
C       CALCULATE STARTING ADDRESS
        IF(ITN(2).EQ.1) THEN
          NADDC = 0
          NADDD = 0
        ELSEIF(ITN(2).EQ.2) THEN
          NADDC = 0
          NADDD = NSKP
        ELSEIF(ITN(2).EQ.3) THEN
          NADDC = NSKP
          NADDD = 0
        ELSEIF(ITN(2).EQ.4) THEN
          NADDC = NSKP
          NADDD = NSKP
        ENDIF
C
C**********************************************************************C
C     LOOP OVER BLOCKS OF C (USE INDEX 6000)                           C
C**********************************************************************C
C
C     LOOP OVER CENTRE C
      DO 6000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE C
        XYZ(1,3) = BXYZ(1,ICNTC)
        XYZ(2,3) = BXYZ(2,ICNTC)
        XYZ(3,3) = BXYZ(3,ICNTC)
C
C     LOOP OVER KQN(C) VALUES
      DO 6000 KC=1,NKAP(ICNTC)
C
C       QUANTUM NUMBERS FOR BLOCK C
        KQN(3) = KAPA(KC,ICNTC)
        IF(KQN(3).LT.0) THEN
          LQN(3) =-KQN(3)-1
        ELSE
          LQN(3) = KQN(3)
        ENDIF
        JQN(3) = 2*IABS(KQN(3))-1
C
C       BASIS EXPONENTS FOR BLOCK C
        NBAS(3) = NFNC(LQN(3),ICNTC)
        DO KBAS=1,NBAS(3)
          EXL(KBAS,3) = BEXL(KBAS,LQN(3),ICNTC)
        ENDDO
C
C     LOOP OVER |MQN(C)| VALUES
      DO 6000 MC=1,IABS(KQN(3))
        MQN(3) = 2*MC-1
C
C     INDEX ASSIGNMENT
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
C
C     STARTING FOCK ADDRESS FOR BASIS FUNCTIONS ON C
      NC1 = LRGE(ICNTC,KC,2*MC-1) + NADDC
      NC2 = LRGE(ICNTC,KC,2*MC  ) + NADDC
C
C     CLEAR ARRAY FOR FIRST CONTRACTION (L->B)
      DO MKB=1,NBAS(3)*NUMO
        DO IJKSPIN=1,8
          B(MKB,IJKSPIN) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     LOOP OVER BLOCKS OF D (USE INDEX 7000)                           C
C**********************************************************************C
C
C     LOOP OVER CENTRE D
      DO 7000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE D
        XYZ(1,4) = BXYZ(1,ICNTD)
        XYZ(2,4) = BXYZ(2,ICNTD)
        XYZ(3,4) = BXYZ(3,ICNTD)
C
C     NUMBER OF NUCLEAR CENTRES INVOLVED IN THIS OVERLAP
      MCNT = NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
C
C     LOOP OVER KQN(D) VALUES
      DO 7000 KD=1,NKAP(ICNTD)
C
C       QUANTUM NUMBERS FOR BLOCK D
        KQN(4) = KAPA(KD,ICNTD)
        IF(KQN(4).LT.0) THEN
          LQN(4) =-KQN(4)-1
        ELSE
          LQN(4) = KQN(4)
        ENDIF
        JQN(4) = 2*IABS(KQN(4))-1
C
C       BASIS EXPONENTS FOR BLOCK D
        NBAS(4) = NFNC(LQN(4),ICNTD)
        DO LBAS=1,NBAS(4)
          EXL(LBAS,4) = BEXL(LBAS,LQN(4),ICNTD)
        ENDDO
C
C     LOOP OVER |MQN(D)| VALUES
      DO 7000 MD=1,IABS(KQN(4))
        MQN(4) = 2*MD-1
C
C     INDEX ASSIGNMENT
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
C     STARTING FOCK ADDRESS FOR BASIS FUNCTIONS ON D
      ND1 = LRGE(ICNTD,KD,2*MD-1) + NADDD
      ND2 = LRGE(ICNTD,KD,2*MD  ) + NADDD
C
C     EQ-COEFFICIENT STARTING ADDRESSES FOR (CD) PAIR
      IF(IT2.EQ.1) THEN
        ICDLL = IAD0LL(ICNTC,ICNTD,KC,KD,MC,MD)
      ELSEIF(IT2.EQ.2) THEN
        ICDLS = IADILS(ICNTC,ICNTD,KC,KD,MC,MD)
      ELSEIF(IT2.EQ.3) THEN
        ICDSL = IADISL(ICNTC,ICNTD,KC,KD,MC,MD)
      ELSEIF(IT2.EQ.4) THEN
        ICDSS = IAD0SS(ICNTC,ICNTD,KC,KD,MC,MD)
      ENDIF
C
C     FLAG READ-IN OF E0(CD) COEFFICIENTS FOR THIS COMPONENT LABEL
      IECD = 1
C
C     EQ-COEFFICIENT PHASE FACTORS FOR PERMUTATION OF R-INTEGRALS
      IF(G2INT.EQ.'COULM') THEN
        PCD1 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
        PCD2 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
      ELSEIF(G2INT.EQ.'BREIT') THEN
        PCD1 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
        PCD2 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
      ENDIF
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON KQN                           C
C**********************************************************************C
C
      GOTO 1111
C     ATOM-CENTRED SELECTION RULES
      IF(MCNT.EQ.1) THEN
C
C       LQN PAIR PARITY (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQN(1)+LQN(2),2)
        IPARCD = MOD(LQN(3)+LQN(4),2)
C
C       LQN PAIR PARITY SELECTION RULE
        IF(IPARAB.NE.IPARCD) THEN
          GOTO 7001
        ENDIF
C
C       JQN TRIANGLE RULE CHECK FOR MULTIPOLE EXPANSION (ATOM-CENTRED)
        NUI = MAX0(IABS(JQN(1)-JQN(2))/2,IABS(JQN(3)-JQN(4))/2)
        NUF = MIN0(    (JQN(1)+JQN(2))/2,    (JQN(3)+JQN(4))/2)
        IF(NUI.GT.NUF) THEN
          GOTO 7001
        ENDIF
C
C       ADDITIONAL LQN SELECTION RULE PARITY ANALYSIS
        ISELK = 0
        DO NU=NUI,NUF
C
C         A AND B: LQN(1)+LQN(2)+NU EVEN OR ODD (0 IF EVEN, 1 IF ODD)
          IPARAB = MOD(LQN(1)+LQN(2)+NU,2)
C
C         C AND D: LQN(3)+LQN(4)+NU EVEN OR ODD (0 IF EVEN, 1 IF ODD)
          IPARCD = MOD(LQN(3)+LQN(4)+NU,2)
C
C         CASE 1: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH ODD (UNLESS NU=0)
          IF(IPARAB.EQ.1.AND.IPARCD.EQ.1.AND.NU.NE.0) ISELK = 1
C
C         CASE 2: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH EVEN
          IF(IPARAB.EQ.0.AND.IPARCD.EQ.0) ISELK = 1
C
        ENDDO
        IF(ISELK.EQ.0) GOTO 7001
C
      ENDIF
1111  CONTINUE
C
C     ADDITIONAL SELECTION RULES FOR CLOSED-SHELL ATOMS
C      IF(MCNT.EQ.1) THEN
C        ISELK = 0
C        IF(KA.EQ.KB.AND.KC.EQ.KD) ISELK = 1
C        IF(KA.EQ.KC.AND.KB.EQ.KD) ISELK = 1
C        IF(KA.EQ.KD.AND.KB.EQ.KC) ISELK = 1
C        IF(ISELK.EQ.0) GOTO 7001
C      ENDIF
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON MQN                           C
C**********************************************************************C
C
      GOTO 1112
C     SPIN PROJECTION CONSERVED ALONG Z-AXIS FOR LINEAR MOLECULES
      IF(ISYM.EQ.1) THEN
        ISELM = 0
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) ISELM = 1
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) ISELM = 1
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) ISELM = 1
        IF(ISELM.EQ.0) GOTO 7001
      ENDIF
C
C     ATOM-CENTRED SELECTION RULES (ONLY APPLIES IF RACAH1 SWITCHED OFF)
      IF(MCNT.EQ.1) THEN
        ISELM = 0
        DO ISGN1=1,2
          DO ISGN2=1,2
            DO ISGN3=1,2
              DO ISGN4=1,2
                MMJA = MQN(1)*((-1)**ISGN1)
                MMJB = MQN(2)*((-1)**ISGN2)
                MMJC = MQN(3)*((-1)**ISGN3)
                MMJD = MQN(4)*((-1)**ISGN4)
                IF(MMJA-MMJB.EQ.MMJD-MMJC) ISELM = 1
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IF(ISELM.EQ.0) GOTO 7001
      ENDIF
1112  CONTINUE
C
C**********************************************************************C
C     ALL LOOPS NOW COMPLETE -- GENERATE BATCH OF ERIs AND CONTRACT    C
C**********************************************************************C
C
C     RESET SCREENING COUNTERS
      DO M=1,NBAS(3)*NBAS(4)
        IMAP(M) = M
        ISCR(M) = 1
      ENDDO
      MAXN = NBAS(3)*NBAS(4)
C
C     OVERRIDE GAUNT REPLACEMENT TOGGLE
      BGAUNT = .FALSE.
C
C     BATCH OF ELECTRON INTERACTION INTEGRALS (IJ|KL) FOR FIXED (IJ)
      WRITE(*,*) 'IN'
      CALL CPU_TIME(T1)
      IF(G2INT.EQ.'COULM') THEN
        CALL ERI(RR,XYZ,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN)
      ELSEIF(G2INT.EQ.'BREIT') THEN
        CALL BII(RR,XYZ,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN,BGAUNT)
      ENDIF
      CALL CPU_TIME(T2)
      TERI = TERI+T2-T1
      WRITE(*,*) 'OUT'
C
C**********************************************************************C
C     FIRST CONTRACTION:                                               C
C     (IJ;T|KL;T') -> (IJ;T|KB;T')                                     C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     FIRST CONTRACTION (NORMAL): (IJ;T|KL;T') -> (IJ;T|KB;T')
C
C     LOOP OVER BASIS FUNCTIONS IN BLOCK C AND OCCUPIED STATES IOCCB
      DO KBAS=1,NBAS(3)
        DO IOCCB=1,NUMO
C
C         FOCK MATRIX ADDRESS FOR IOCCB
          IB = MINO-1+IOCCB+NSKP
C
C         LIST ADDRESS FOR THIS KBAS,IOCCB COMBINATION
          MKB = (KBAS-1)*NUMO + IOCCB
C
C         CONTRACT OVER ALL LBAS IN BLOCK D
          DO LBAS=1,NBAS(4)
C
C           LIST ADDRESS FOR THIS KBAS,LBAS COMBINATION
            M = (KBAS-1)*NBAS(4) + LBAS
C
C           (--|-B) = (--|--) + (--|-+)
            B(MKB,1) = B(MKB,1) + RR(M, 1)*COEF(ND1+LBAS,IB)
     &                          + RR(M, 2)*COEF(ND2+LBAS,IB)
C
C           (+-|-B) = (+-|--) + (+-|-+)
            B(MKB,2) = B(MKB,2) + RR(M, 9)*COEF(ND1+LBAS,IB)
     &                          + RR(M,10)*COEF(ND2+LBAS,IB)
C
C           (-+|-B) = (-+|--) + (-+|-+)
            B(MKB,3) = B(MKB,3) + RR(M, 5)*COEF(ND1+LBAS,IB)
     &                          + RR(M, 6)*COEF(ND2+LBAS,IB)
C
C           (++|-B) = (++|--) + (++|-+)
            B(MKB,4) = B(MKB,4) + RR(M,13)*COEF(ND1+LBAS,IB)
     &                          + RR(M,14)*COEF(ND2+LBAS,IB)
C
C           (--|+B) = (--|+-) + (--|++)
            B(MKB,5) = B(MKB,5) + RR(M, 3)*COEF(ND1+LBAS,IB)
     &                          + RR(M, 4)*COEF(ND2+LBAS,IB)
C
C           (+-|+B) = (+-|+-) + (+-|++)
            B(MKB,6) = B(MKB,6) + RR(M,11)*COEF(ND1+LBAS,IB)
     &                          + RR(M,12)*COEF(ND2+LBAS,IB)
C
C           (-+|+B) = (-+|+-) + (-+|++)
            B(MKB,7) = B(MKB,7) + RR(M, 7)*COEF(ND1+LBAS,IB)
     &                          + RR(M, 8)*COEF(ND2+LBAS,IB)
C
C           (++|+B) = (++|+-) + (++|++)
            B(MKB,8) = B(MKB,8) + RR(M,15)*COEF(ND1+LBAS,IB)
     &                          + RR(M,16)*COEF(ND2+LBAS,IB)
C
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NBAS(4)*NUMO*NBAS(3)
C
      CALL CPU_TIME(T2)
      TCN1 = TCN1+T2-T1
C
C     SKIP POINT FOR KQN AND MQN SELECTION RULES
7001  CONTINUE
C
C     FIRST CONTRACTION COMPLETE FOR THIS D BLOCK
7000  CONTINUE
C
C**********************************************************************C
C     SECOND CONTRACTION:                                              C
C     (IJ;T|KB;T') -> (IJ;T|SB)                                        C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     SECOND CONTRACTION (NORMAL): (IJ;T|KB;T') -> (IJ;T|SB;T')
C
C     LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
      DO IOCCB=1,NUMO
        DO IVRTS=1,NUMV
C
C         FOCK MATRIX ADDRESS FOR IVRTS
          IS = MINV-1+IVRTS+NSKP
C
C         LIST ADDRESS FOR THIS IVRTS,IOCCB COMBINATION IN B1
          MSB = (IOCCB-1)*NUMV + IVRTS
C
C         CONTRACT OVER ALL KBAS IN BLOCK C
          DO KBAS=1,NBAS(3)
C
C           LIST ADDRESS FOR THIS KBAS,IOCCB COMBINATION
            MKB = (KBAS-1)*NUMO + IOCCB
C
C           (--|SB) = (--|-B) + (--|+B)
            SB(MIJ,MSB,1) = SB(MIJ,MSB,1)
     &                             + B(MKB,1)*DCONJG(COEF(NC1+KBAS,IS))
     &                             + B(MKB,5)*DCONJG(COEF(NC2+KBAS,IS))
C
C           (-+|SB) = (-+|-B) + (-+|+B)
            SB(MIJ,MSB,2) = SB(MIJ,MSB,2)
     &                             + B(MKB,3)*DCONJG(COEF(NC1+KBAS,IS))
     &                             + B(MKB,7)*DCONJG(COEF(NC2+KBAS,IS))
C
C           (+-|SB) = (+-|-B) + (+-|+B)
            SB(MIJ,MSB,3) = SB(MIJ,MSB,3)
     &                             + B(MKB,2)*DCONJG(COEF(NC1+KBAS,IS))
     &                             + B(MKB,6)*DCONJG(COEF(NC2+KBAS,IS))
C
C           (++|SB) = (++|-B) + (++|+B)
            SB(MIJ,MSB,4) = SB(MIJ,MSB,4)
     &                             + B(MKB,4)*DCONJG(COEF(NC1+KBAS,IS))
     &                             + B(MKB,8)*DCONJG(COEF(NC2+KBAS,IS))
C
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NUMO*NUMV*NBAS(3)
C
      CALL CPU_TIME(T2)
      TCN2 = TCN2+T2-T1
C
C     SECOND CONTRACTION COMPLETE FOR THIS C BLOCK AND T'T' OVERLAP
6000  CONTINUE
5000  CONTINUE
C
C     SECOND CONTRACTION COMPLETE FOR ALL (IBAS,JBAS)
4000  CONTINUE
C
C     WRITE PERCENTAGE COMPLETION
      WRITE(*,*) 'GET TO HERE'
      WRITE(6,*) 'Completion: ',DFLOAT(INEST)/DFLOAT(NNEST),' %'
      WRITE(7,*) 'Completion: ',DFLOAT(INEST)/DFLOAT(NNEST),' %'
C
C**********************************************************************C
C     THIRD CONTRACTION:                                               C
C     (IJ;T|SB) -> (IA;T|SB)  AND  (JI;T|SB) -> (JA;T|SB)              C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     THIRD CONTRACTION (DIRECT): (IJ;T|SB) -> (IA;T|SB)
C
C     CLEAR ARRAY FOR THIRD CONTRACTION (DIRECT)
      DO MASB=1,NUMV*NUMO*NUMO
        DO IBAS=1,NBAS(1)
          DO ISPIN=1,2
            ASB1(IBAS,MASB,ISPIN) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C     LOOP OVER OCCUPIED STATES IOCCA
      DO IOCCA=1,NUMO
C
C       FOCK MATRIX ADDRESS FOR IOCCA
        IA = MINO-1+IOCCA+NSKP
C
C       LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
        DO IOCCB=1,IOCCA
          DO IVRTS=1,NUMV
C
C           LIST ADDRESS FOR THIS IVRTS,IOCCB COMBINATION
            MSB = (IOCCB-1)*NUMV+IVRTS
C
C           LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MSB
            MASB = (IOCCA-1)*NUMV*NUMO + MSB
C
C           CONTRACT OVER ALL (IBAS,JBAS) IN BLOCKS A AND B
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
C
C               LIST ADDRESS FOR THIS IBAS,JBAS COMBINATION
                MIJ = (IBAS-1)*NBAS(2)+JBAS
C
C               (-A|SB) = (--|SB) + (-+|SB)
                ASB1(IBAS,MASB,1) = ASB1(IBAS,MASB,1)
     &                            +     SB(MIJ,MSB,1)*COEF(NB1+JBAS,IA)
     &                            +     SB(MIJ,MSB,2)*COEF(NB2+JBAS,IA)
C
C               (+A|SB) = (+-|SB) + (++|SB)
                ASB1(IBAS,MASB,2) = ASB1(IBAS,MASB,2)
     &                            +     SB(MIJ,MSB,3)*COEF(NB1+JBAS,IA)
     &                            +     SB(MIJ,MSB,4)*COEF(NB2+JBAS,IA)
C
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NUMO*(NUMO+1)*NUMV*NBAS(1)*NBAS(2)/2 
C
C     THIRD CONTRACTION (SWAP): (JI;T|SB) -> (JA;T|SB)
      IF(G2INT.EQ.'COULM'.AND.IQ1.EQ.IQ2) GOTO 3002
      IF(G2INT.EQ.'BREIT') GOTO 3002
C
C     CLEAR ARRAY FOR THIRD CONTRACTION (SWAP)
      DO MASB=1,NUMV*NUMO*NUMO
        DO JBAS=1,NBAS(2)
          DO JSPIN=1,2
            ASB2(JBAS,MASB,JSPIN) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C     LOOP OVER OCCUPIED STATES IOCCA
      DO IOCCA=1,NUMO
C
C       FOCK MATRIX ADDRESS FOR IOCCA
        IA = MINO-1+IOCCA+NSKP
C
C       LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
        DO IOCCB=1,IOCCA
          DO IVRTS=1,NUMV
C
C           LIST ADDRESS FOR THIS IVRTS,IOCCB COMBINATION
            MSB = (IOCCB-1)*NUMV+IVRTS
C
C           LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MSB
            MASB = (IOCCA-1)*NUMV*NUMO + MSB
C
C           CONTRACT OVER ALL (IBAS,JBAS) IN BLOCKS A AND B
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
C
C               LIST ADDRESS FOR THIS IBAS,JBAS COMBINATION
                MIJ = (IBAS-1)*NBAS(2)+JBAS
C
C               (+A|SB) = PCD*{(+A|SB)} = PCD*{(++|SB) + (-+|SB)}
                ASB2(JBAS,MASB,1) = ASB2(JBAS,MASB,1)
     &                           + PAB1*SB(MIJ,MSB,4)*COEF(NA1+IBAS,IA)
     &                           + PAB2*SB(MIJ,MSB,2)*COEF(NA2+IBAS,IA)
C
C               (-A|SB) = PCD*{(-A|SB)} = PCD*{(+-|SB) + (--|SB)}
                ASB2(JBAS,MASB,2) = ASB2(JBAS,MASB,2)
     &                           + PAB2*SB(MIJ,MSB,3)*COEF(NA1+IBAS,IA)
     &                           + PAB1*SB(MIJ,MSB,1)*COEF(NA2+IBAS,IA)
C
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     SKIP POINT FOR IQ1=IQ2
3002  CONTINUE
C
      CALL CPU_TIME(T2)
      TCN3 = TCN3+T2-T1
C
C**********************************************************************C
C     FOURTH CONTRACTION:                      ~                       C
C     (IA;T|SB) -> (RA|SB)  AND  (JA;T|SB) -> (RA|SB)                  C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     FOURTH CONTRACTION (DIRECT): (IA;T|SB) -> (RA|SB)
C
C     LOOP OVER OCCUPIED STATES IOCCA AND VIRTUAL STATES IVRTR
      DO IOCCA=1,NUMO
        DO IVRTR=1,NUMV
C
C         FOCK MATRIX ADDRESS FOR IVRTR
          IR = MINV-1+IVRTR+NSKP
C
C         LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
          DO IOCCB=1,IOCCA
            DO IVRTS=1,NUMV
C
C             LIST ADDRESS FOR THIS IVRTS,IOCCB COMBINATION
              MSB = (IOCCB-1)*NUMV + IVRTS
C
C             LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MSB
              MASB = (IOCCA-1)*NUMV*NUMO + MSB
C
C             LIST ADDRESS FOR THIS IVRTR,IOCCA AND THE ABOVE MSB
              MRASB = (IOCCA-1)*NUMV*IOCCA*NUMV/2
     &              + (IVRTR-1)*NUMV*IOCCA + MSB
C
C             CONTRACT OVER ALL IBAS IN BLOCK A
              DO IBAS=1,NBAS(1)
C
C               (RA|SB) = (-A|SB) + (+A|SB)
                RASB(MRASB) = RASB(MRASB)
     &                    + ASB1(IBAS,MASB,1)*DCONJG(COEF(NA1+IBAS,IR))
     &                    + ASB1(IBAS,MASB,2)*DCONJG(COEF(NA2+IBAS,IR))
C
              ENDDO
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NUMO*NUMV*(NUMO+1)*NUMV*NBAS(1)/2 
C                                              ~
C     FOURTH CONTRACTION (SWAP): (JA;T|SB) -> (RA|SB)
      IF(G2INT.EQ.'COULM'.AND.IQ1.EQ.IQ2) GOTO 3003
      IF(G2INT.EQ.'BREIT') GOTO 3003
C
C     LOOP OVER OCCUPIED STATES IOCCA AND VIRTUAL STATES IVRTR
      DO IOCCA=1,NUMO
        DO IVRTR=1,NUMV
C
C         FOCK MATRIX ADDRESS FOR IVRTR
          IR = MINV-1+IVRTR+NSKP
C
C         LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
          DO IOCCB=1,IOCCA
            DO IVRTS=1,NUMV
C
C             LIST ADDRESS FOR THIS IVRTS,IOCCB COMBINATION
              MSB = (IOCCB-1)*NUMV+IVRTS
C
C             LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MSB
              MASB = (IOCCA-1)*NUMV*NUMO + MSB
C
C             LIST ADDRESS FOR THIS IVRTR,IOCCA AND THE ABOVE MSB
              MRASB = (IOCCA-1)*NUMV*IOCCA*NUMV/2
     &              + (IVRTR-1)*NUMV*IOCCA + MSB
C
C             CONTRACT OVER ALL JBAS IN BLOCK B
              DO JBAS=1,NBAS(2)
C
C               (PA|SB) = (+A|SB) + (-A|SB)
                RASB(MRASB) = RASB(MRASB)
     &                    + ASB2(JBAS,MASB,1)*DCONJG(COEF(NB1+JBAS,IR))
     &                    + ASB2(JBAS,MASB,2)*DCONJG(COEF(NB2+JBAS,IR))
C
              ENDDO
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     SKIP POINT FOR IQ1=IQ2
3003  CONTINUE
C
      CALL CPU_TIME(T2)
      TCN4 = TCN4+T2-T1
C
C     SKIP POINT FOR IQ1.LT.IQ2
3001  CONTINUE

      WRITE(*,*) 'END UP HERE'
C
C     ALL CONTRIBUTIONS FROM THIS CLASS (A,B) BLOCK NOW ACCOUNTED FOR
3000  CONTINUE
2000  CONTINUE
1000  CONTINUE
C
C**********************************************************************C
C     CALCULATE SECOND-ORDER PAIR CORRELATION ENERGY                   C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     FOR EACH IOCCA,IOCCB PAIR, SUM OVER IVRTR AND IVRTS CONTRIBUTIONS
C
C     LOOP OVER OCCUPIED STATES IOCCA AND VIRTUAL STATES IVRTR
      DO IOCCA=1,NUMO
        DO IVRTR=1,NUMV
C
C         FOCK MATRIX ADDRESSES
          IA = MINO-1+IOCCA+NSKP
          IR = MINV-1+IVRTR+NSKP
C
C         LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
          DO IOCCB=1,IOCCA
            DO IVRTS=1,NUMV
C
C             FOCK MATRIX ADDRESSES
              IB = MINO-1+IOCCB+NSKP
              IS = MINV-1+IVRTS+NSKP
C
C             MAIN LIST ADDRESS FOR THIS IOCCA,IVRTR,IOCCB,IVRTS
              MRASB = (IOCCA-1)*NUMV*IOCCA*NUMV/2
     &              + (IVRTR-1)*NUMV*IOCCA + (IOCCB-1)*NUMV + IVRTS
C
C             MAIN LIST ADDRESS FOR THIS IOCCA,IVRTS,IOCCB,IVRTR
              MSARB = (IOCCA-1)*NUMV*IOCCA*NUMV/2
     &              + (IVRTS-1)*NUMV*IOCCA + (IOCCB-1)*NUMV + IVRTR
C
C             MAKE SURE THIS ADDRESS DOESN'T EXCEED BOUNDS
              IF(MRASB.GT.NRASB) THEN
                WRITE(6,*) 'In MBPT2: address exceeds dimension.',MRASB
                WRITE(7,*) 'In MBPT2: address exceeds dimension.',MRASB
C              ELSEIF(MSARB.GT.NRASB) THEN
C                WRITE(6,*) 'In MBPT2: address exceeds dimension.',MSARB
C                WRITE(7,*) 'In MBPT2: address exceeds dimension.',MSARB
              ENDIF
C
C             NUMERATOR FOR DIRECT AND EXCHANGE CONTRIBUTIONS
              RNUMD = DREAL(RASB(MRASB)*DCONJG(RASB(MRASB)))
              RNUMX =-DREAL(RASB(MRASB)*DCONJG(RASB(MSARB)))
C
C             DENOMINATOR FOR BOTH CONTRIBUTIONS
              EABRS = EIGN(IA) + EIGN(IB) - EIGN(IR) - EIGN(IS)
C
C             ADD TO DIRECT AND EXCHANGE BINS
              EAB2(IOCCA,IOCCB,4) = EAB2(IOCCA,IOCCB,4) + RNUMD/EABRS
              EAB2(IOCCA,IOCCB,5) = EAB2(IOCCA,IOCCB,5) + RNUMX/EABRS
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     FILL IN THE OTHER HALF OF THE ARRAY AND CALCULATE TOTALS
      E1S = 0.0D0
      E2D = 0.0D0
      E2X = 0.0D0
      E2S = 0.0D0
      DO IOCCA=1,NUMO
        DO IOCCB=1,IOCCA
C
C         INTERMEDIATE VALUES
          EAB1TOT = EAB2(IOCCA,IOCCB,3)
          EAB2DIR = EAB2(IOCCA,IOCCB,4)
          EAB2XCH = EAB2(IOCCA,IOCCB,5)
          EAB2SUM = EAB2DIR + EAB2XCH
C
C         PUT THESE INTO EAB2 AND ADD CONTRIBUTION TO E2
          EAB2(IOCCA,IOCCB,6) = EAB2SUM
          IF(IOCCA.NE.IOCCB) THEN
            EAB2(IOCCB,IOCCA,3) = EAB1TOT
            EAB2(IOCCB,IOCCA,4) = EAB2DIR
            EAB2(IOCCB,IOCCA,5) = EAB2XCH
            EAB2(IOCCB,IOCCA,6) = EAB2SUM
            E1S = E1S +       EAB1TOT
            E2D = E2D +       EAB2DIR
            E2X = E2X +       EAB2XCH
            E2S = E2S +       EAB2SUM
          ELSE
            E1S = E1S + 0.5D0*EAB1TOT
            E2D = E2D + 0.5D0*EAB2DIR
            E2X = E2X + 0.5D0*EAB2XCH
            E2S = E2S + 0.5D0*EAB2SUM
          ENDIF
        ENDDO
      ENDDO
C
C     WRITE RESULTS OF EAB ENERGIES TO AN EXTERNAL FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_MBPT2.dat',STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO IOCCA=1,NUMO
        DO IOCCB=1,NUMO
          WRITE(8, *) (EAB2(IOCCA,IOCCB,N),N=1,6)
        ENDDO
      ENDDO
      CLOSE(UNIT=8)
C
C**********************************************************************C
C     CALCULATE SECOND-ORDER SINGLE ORBITAL ENERGY                     C
C**********************************************************************C
C
C     FOR EACH IOCCA, SUM OVER THE IOCCB CONTRIBUTIONS
      DO IOCCA=1,NUMO
        DO N=1,6
          EA2(IOCCA,N) = 0.0D0
          DO IOCCB=1,NUMO
            EA2(IOCCA,N) = EA2(IOCCA,N) + EAB2(IOCCA,IOCCB,N)
          ENDDO
        ENDDO
      ENDDO
C
      CALL CPU_TIME(T2)
      TSUM = TSUM+T2-T1
C
C**********************************************************************C
C     TERMINAL OUTPUT SUMMARY                                          C
C**********************************************************************C
C
C     MBPT2 PAIRWISE SUMMARY
20    FORMAT(1X,A,9X,A,9X,A,9X,A,9X,A)
21    FORMAT(' (',I2,',',I2,')',3X,F13.7,5X,F11.7,5X,F11.7,4X,F13.7)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',25),'MBPT2 pairwise summary'
      WRITE(7, *) REPEAT(' ',25),'MBPT2 pairwise summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) '( a, b)','E1G(ab)','E2J(ab)','E2K(ab)','E2G(ab)'
      WRITE(7,20) '( a, b)','E1G(ab)','E2J(ab)','E2K(ab)','E2G(ab)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO IOCCA=1,NUMO
        IANUM = IOCCA+MINO-1
        DO IOCCB=1,IOCCA
          IBNUM = IOCCB+MINO-1
          WRITE(6,21) IANUM,IBNUM,(EAB2(IOCCA,IOCCB,N),N=3,6)
          WRITE(7,21) IANUM,IBNUM,(EAB2(IOCCA,IOCCB,N),N=3,6)
        ENDDO
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     MBPT2 SINGLE-PARTICLE SUMMARY
30    FORMAT(1X,A,10X,A,10X,A,10X,A,10X,A)
31    FORMAT('  ',I2,'    ',3X,F13.7,5X,F11.7,5X,F11.7,4X,F13.7)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',21),'MBPT2 single particle summary'
      WRITE(7, *) REPEAT(' ',21),'MBPT2 single particle summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,30) '  a    ','E1G(a)','E2J(a)','E2K(a)',' E2G(a)'
      WRITE(7,30) '  a    ','E1G(a)','E2J(a)','E2K(a)',' E2G(a)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO IOCCA=1,NUMO
        IANUM = IOCCA+MINO-1
        WRITE(6,31) IANUM,(EA2(IOCCA,N),N=3,6)
        WRITE(7,31) IANUM,(EA2(IOCCA,N),N=3,6)
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     MBPT2 TOTAL SECOND ORDER CORRELATION ENERGY
32    FORMAT(' total  ',3X,F13.7,5X,F11.7,5X,F11.7,4X,F13.7)

      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',17),'MBPT2 second order correlation energy'
      WRITE(7, *) REPEAT(' ',17),'MBPT2 second order correlation energy'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,32) E1S,E2D,E2X,E2S
      WRITE(7,32) E1S,E2D,E2X,E2S
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     MBPT2 LABOUR ANALYSIS
      CALL CPU_TIME(TFIN)
      TTOT = TFIN - TBEG
      TOTH = TTOT - (TERI + TCN1 + TCN2 + TCN3 + TCN4 + TSUM)

40    FORMAT(1X,A,24X,A)
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',26),'MBPT2 labour analysis'
      WRITE(7, *) REPEAT(' ',26),'MBPT2 labour analysis'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
      WRITE(6,40) 'ERI construction - (IJ|KL)      ',HMS(TERI)
      WRITE(7,40) 'ERI construction - (IJ|KL)      ',HMS(TERI)
      WRITE(6,40) '1st contraction  - (IJ|KB)      ',HMS(TCN1)
      WRITE(7,40) '1st contraction  - (IJ|KB)      ',HMS(TCN1)
      WRITE(6,40) '2nd contraction  - (IJ|SB)      ',HMS(TCN2)
      WRITE(7,40) '2nd contraction  - (IJ|SB)      ',HMS(TCN2)
      WRITE(6,40) '3rd contraction  - (IA|SB)      ',HMS(TCN3)
      WRITE(7,40) '3rd contraction  - (IA|SB)      ',HMS(TCN3)
      WRITE(6,40) '4th contraction  - (RA|SB)      ',HMS(TCN4)
      WRITE(7,40) '4th contraction  - (RA|SB)      ',HMS(TCN4)
      WRITE(6,40) 'Virtual orbital sum - E2(A,B)   ',HMS(TSUM)
      WRITE(7,40) 'Virtual orbital sum - E2(A,B)   ',HMS(TSUM)
      WRITE(6,40) 'Other                           ',HMS(TOTH)
      WRITE(7,40) 'Other                           ',HMS(TOTH)
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,40) 'Total MBPT2 time                ',HMS(TTOT)
      WRITE(7,40) 'Total MBPT2 time                ',HMS(TTOT)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
      RETURN
      END

