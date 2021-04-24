      SUBROUTINE MBPT2
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
C  DFNOTE: THIS VERSION WAS ADAPTED FROM SOME OLD BERTHA CODE.         C
C -------------------------------------------------------------------- C
C  MBPT2 CALCULATES THE SECOND-ORDER PAIR CORRELATION ENERGIES OF A    C
C  CLOSED-SHELL SYSTEM USING A RELATIVISTIC ADAPTION OF THE DIRECT MP2 C
C  ALGORITHM OF HEAD-GORDON, POPLE AND FRISCH. SPINOR INTEGRALS ARE    C
C  GENERATED IN ATOMIC SYMMETRY BLOCKS AND TRANSFORMED INTO MATRIX     C
C  ELEMENTS BY SEQUENTIAL ONE-INDEX TRANSFORMATIONS. MBPT2 EXPLOITS    C
C  THE SYMMETRY (AB|CD)=(AB|DC) BUT NOT THE SYMMETRY (AB|CD)=(CD|AB)   C
C  (AND PERMUTATIONS). THE LOOP STRUCTURE IS BASED ON THE FOCK MATRIX  C
C  ROUTINE SCFMAT, WHICH LOOPS OVER THE UNIQUE LIST OF BLOCK LABELS.   C  
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=7000000)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 C(MDM,MDM),RR(MB2,16)
C
      COMMON/COEF/C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS
      COMMON/MAKE/IEAB,IECD,IERC,IRIJ(MBS,MBS)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,IEQS,IPAR,ICOR,ILEV
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND,
     &            TSCR
C
      COMPLEX*16 RASB(NOCC*NOCC*NVRT*NVRT)
      COMPLEX*16 ASB1(MBS,NOCC*NOCC*NVRT,2),ASB2(MBS,NOCC*NOCC*NVRT,2)
      COMPLEX*16 SB(MB2,NOCC*NVRT,4)
      COMPLEX*16 B1(NOCC*MBS,8),B2(NOCC*MBS,8)
C
      DIMENSION EAB2(NOCC,NOCC,6)
      DIMENSION KQN(4),LQN(4),MQN(4),NBAS(4),ITQN(2)
      DIMENSION EXPT(MBS,4),XYZ(3,4)
      DIMENSION INDEX(MCT,-(MKP+1)/2:(MKP+1)/2,MKP)
C
C     CONSTRUCT ORDERED INDEX SYSTEM FOR ALL POSSIBLE {XYZ,KQN,MQN}
      ICOUNT = 0
C
C     LOOP OVER ATOMIC CENTERS
      DO ICNT=1,NCNT
C
C       LOOP OVER KAPPA VALUES FOR THIS ATOMIC CENTER
        DO KN=1,NKAP(ICNT)
C
C         IMPORT KAPPA, MAXIMUM MQN
          KAPPA = KVALS(KN,ICNT)
          MJMAX = 2*IABS(KAPPA)-1
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
C     LOOP OVER COMPONENT OVERLAP OPTIONS (INDEX 1000)                 C
C**********************************************************************C
C
C     COMPONENT OVERLAP LABELS TO LOOP OVER
      IF(HMLTN.EQ.'BARE') THEN
        RETURN
      ELSEIF(HMLTN.EQ.'NORL') THEN
        ITSTRT = 1
        ITSTOP = 1
        ITSKIP = 1
      ELSE
        ITSTRT = 1
        ITSTOP = 4
        ITSKIP = 3
      ENDIF
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
C     LOOP OVER COMPONENT LABEL FOR C AND D: T'T' = LL (1) or SS (4)
      DO 1000 IT2=ITSTRT,ITSTOP,ITSKIP
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
C     COMPONENT OVERLAP INDEX {(LL|LL)=1,(LL|SS)=2,(SS|LL)=3,(SS|SS)=4}
      ITT = (2*IT1+IT2)/3
C
C**********************************************************************C
C     LOOP OVER ATOMIC CENTERS WITH KQN A AND B (USE INDEX 2000)       C
C**********************************************************************C
C
C     LOOP OVER CENTER A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 2000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
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
        DO JBAS=1, NBAS(2)
          EXPT(JBAS,2) = EXPSET(JBAS,LQN(2)+1,ICNTB)
        ENDDO
C
C**********************************************************************C
C     LOOP OVER ATOMIC CENTERS WITH KQN C AND D (USE INDEX 3000)       C
C**********************************************************************C
C
C     LOOP OVER CENTER C
      DO 3000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER C
        XYZ(1,3) = COORD(1,ICNTC)
        XYZ(2,3) = COORD(2,ICNTC)
        XYZ(3,3) = COORD(3,ICNTC)
C
C     LOOP OVER CENTER D
      DO 3000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER D
        XYZ(1,4) = COORD(1,ICNTD)
        XYZ(2,4) = COORD(2,ICNTD)
        XYZ(3,4) = COORD(3,ICNTD)
C      
C     LOOP OVER KQN(C) VALUES
      DO 3000 KC=1,NKAP(ICNTC)
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
      DO 3000 KD=1,NKAP(ICNTD)
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
C     RESET RC(AB|CD) CLASS CALCULATION INDICATORS
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          IRIJ(IBAS,JBAS) = 1
        ENDDO
      ENDDO
C
C     BUILD NEW RC(AB|CD) CLASS AT NEXT OPPORTUNITY
      IF(IEQS.EQ.0) THEN
C       E-COEFFICIENTS NOT SAVED TO LOCAL FILE -> DON'T SAVE RC(AB|CD)
        IERC = 0
      ELSEIF(LQN(1)+LQN(2)+LQN(3)+LQN(4).EQ.0) THEN
C       TURN OFF RC(AB|CD) LOCAL FILE PROCESS (ONE MQN COMBINATION)
        IERC = 0
      ELSE
C       ALLOW RC(AB|CD) LOCAL FILE PROCESS
        IERC = 1
      ENDIF
C
C**********************************************************************C
C     LOOP OVER ALL |MQN| VALUES (INDEX 4000)                          C
C**********************************************************************C
C
C     LOOP OVER |MQN(A)| VALUES
      DO 4000 MA=1,IABS(KQN(1))
        MQN(1) = 2*MA-1
C
C     LOOP OVER |MQN(B)| VALUES
      DO 4000 MB=1,IABS(KQN(2))
        MQN(2) = 2*MB-1
C
C     CALCULATE NEW BLOCK OF E(AB|  ) COEFFS AT NEXT OPPORTUNITY
      IEAB  = 1
      IABLL = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB)
      IABSS = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB)
C
C     LOOP OVER |MQN(C)| VALUES
      DO 4000 MC=1,IABS(KQN(3))
        MQN(3) = 2*MC-1
C
C     LOOP OVER |MQN(D)| VALUES
      DO 4000 MD=1,IABS(KQN(4))
        MQN(4) = 2*MD-1
C
C     CALCULATE NEW BLOCK OF E(CD|  ) COEFFS AT NEXT OPPORTUNITY
      IECD  = 1
      ICDLL = IAD0LL(ICNTC,ICNTD,KC,KD,MC,MD)
      ICDSS = IAD0SS(ICNTC,ICNTD,KC,KD,MC,MD)
C
C**********************************************************************C
C     CANNOT USE SELECTION RULES IN MBPT (DOUBLE EXCITATIONS!)         C
C**********************************************************************C
C
C**********************************************************************C
C     INDEX ASSIGNMENT AND IDENTIFICATION OF ERI SYMMETRY CLASSES      C
C**********************************************************************C
C
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
C     SKIP CONTRIBUTIONS THAT ARISE BY PERMUTATION OF INTEGRALS
      IF(IQ1.LT.IQ2) GOTO 4001
      IF(IQ3.LT.IQ4) GOTO 4001
C
C**********************************************************************C
C     PHASE FACTORS APPROPRIATE TO INTEGRAL SYMMETRY CLASSES           C
C**********************************************************************C
C
C     EQ-COEFFICIENT PHASE FACTORS FOR PERMUTATION OF R-INTEGRALS
      PAB1 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PAB2 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      PCD1 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PCD2 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
C     FOCK ADDRESS FOR EACH BASIS FUNCTION (WITH SPIN PROJECTION)
      NA1 = LARGE(ICNTA,KA,2*MA-1) + NADDAB
      NB1 = LARGE(ICNTB,KB,2*MB-1) + NADDAB
      NC1 = LARGE(ICNTC,KC,2*MC-1) + NADDCD
      ND1 = LARGE(ICNTD,KD,2*MD-1) + NADDCD
C
      NA2 = LARGE(ICNTA,KA,2*MA  ) + NADDAB
      NB2 = LARGE(ICNTB,KB,2*MB  ) + NADDAB
      NC2 = LARGE(ICNTC,KC,2*MC  ) + NADDCD
      ND2 = LARGE(ICNTD,KD,2*MD  ) + NADDCD
C
C     RESET AN INTERMEDIATE ARRAY FOR (IBAS,JBAS) CONTRIBUTIONS
      DO MSB=1,NOCC*NVRT
        DO M=1,NBAS(1)*NBAS(2) 
          DO KL=1,4
            SB(M,MSB,KL) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS IN BLOCKS A AND B (INDEX 5000)         C
C**********************************************************************C
C
      DO 5000 IBAS=1,NBAS(1)
      DO 5000 JBAS=1,NBAS(2)
C
C     LIST INDEX FOR THIS OVERLAP
      MIJ = (IBAS-1)*NBAS(2) + JBAS
C
C**********************************************************************C
C     CALCULATE A BATCH OF INTEGRALS (IJ|KL) FOR FIXED (IJ)            C
C**********************************************************************C
C
C     GENERATE BATCH OF ELECTRON REPULSION INTEGRALS
      CALL ERI(RR,XYZ,KQN,MQN,EXPT,NBAS,ITQN,IBAS,JBAS)
C
      DO MKB=1,NBAS(3)*NOCC
        DO KL=1,8
          B1(MKB,KL) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO

      DO MLB=1,NBAS(4)*NOCC
        DO KL=1,8
          B2(MLB,KL) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     PERFORM ONE-INDEX CONTRACTION:    (IJ;T|KL;T') -> (IJ;T|KB;T')   C
C                                       (IJ;T|LK;T') -> (IJ;T|LB;T')   C
C**********************************************************************C
C
C     DIRECT EVALUATION: (IJ;T|KL;T') -> (IJ;T|KB;T')
      MKB = 0
      DO KBAS=1,NBAS(3)
        DO IOCCB=1,NOCC
          IB = IOCCB+NSHIFT
C
C         LIST INDEX FOR THIS KBAS,IOCCB PAIR
          MKB = MKB+1
C
C         LOOP OVER BASIS FUNCTIONS IN BLOCK D AND CONTRACT ERI BLOCK
          DO LBAS=1,NBAS(4)
C
C           LIST INDEX FOR THIS KBAS,LBAS PAIR
            M = (KBAS-1)*NBAS(4)  +LBAS
C
C           (--|-B) = (--|--) + (--|-+)
            B1(MKB,1) = B1(MKB,1) +      RR(M, 1)*C(ND1+LBAS,IB)
     &                            +      RR(M, 2)*C(ND2+LBAS,IB)
C           (+-|-B) = (+-|--) + (+-|-+)
            B1(MKB,2) = B1(MKB,2) +      RR(M, 9)*C(ND1+LBAS,IB)
     &                            +      RR(M,10)*C(ND2+LBAS,IB)
C           (-+|-B) = (-+|--) + (-+|-+)
            B1(MKB,3) = B1(MKB,3) +      RR(M, 5)*C(ND1+LBAS,IB)
     &                            +      RR(M, 6)*C(ND2+LBAS,IB)
C           (++|-B) = (++|--) + (++|-+)
            B1(MKB,4) = B1(MKB,4) +      RR(M,13)*C(ND1+LBAS,IB)
     &                            +      RR(M,14)*C(ND2+LBAS,IB)
C           (--|+B) = (--|+-) + (--|++)
            B1(MKB,5) = B1(MKB,5) +      RR(M, 3)*C(ND1+LBAS,IB)
     &                            +      RR(M, 4)*C(ND2+LBAS,IB)
C           (+-|+B) = (+-|+-) + (+-|++)
            B1(MKB,6) = B1(MKB,6) +      RR(M,11)*C(ND1+LBAS,IB)
     &                            +      RR(M,12)*C(ND2+LBAS,IB)
C           (-+|+B) = (-+|+-) + (-+|++)
            B1(MKB,7) = B1(MKB,7) +      RR(M, 7)*C(ND1+LBAS,IB)
     &                            +      RR(M, 8)*C(ND2+LBAS,IB)
C           (++|+B) = (++|+-) + (++|++)
            B1(MKB,8) = B1(MKB,8) +      RR(M,15)*C(ND1+LBAS,IB)
     &                            +      RR(M,16)*C(ND2+LBAS,IB)
          ENDDO
C
        ENDDO
      ENDDO
C
C     INTEGRAL SYMMETRY CASE: (IJ;T|LK;T') -> (IJ;T|LB;T')
      IF(IQ3.EQ.IQ4) GOTO 5100
C
      MLB = 0
      DO LBAS=1,NBAS(4)
        DO IOCCB=1,NOCC
          IB = IOCCB+NSHIFT
C
C         LIST INDEX FOR THIS LBAS,IOCCB PAIR
          MLB = MLB+1
C
C         LOOP OVER BASIS FUNCTIONS IN BLOCK C AND CONTRACT ERI BLOCK
          DO KBAS=1,NBAS(3)
C
C           LIST INDEX FOR THIS KBAS,LBAS PAIR
            M = (KBAS-1)*NBAS(4) + LBAS
C
C           (--|+B) = PAB*{(--|B+)} = PAB*{(--|++) + ((--|-+))}
            B2(MLB,1) = B2(MLB,1) + PAB1*RR(M, 4)*C(NC1+KBAS,IB)
     &                            + PAB2*RR(M, 2)*C(NC2+KBAS,IB)
C           (+-|+B) = PAB*{(+-|B+)} = PAB*{(+-|++) + ((+-|-+))}
            B2(MLB,2) = B2(MLB,2) + PAB1*RR(M,12)*C(NC1+KBAS,IB)
     &                            + PAB2*RR(M,10)*C(NC2+KBAS,IB)
C           (-+|+B) = PAB*{(-+|B+)} = PAB*{(-+|++) + ((-+|-+))}
            B2(MLB,3) = B2(MLB,3) + PAB1*RR(M, 8)*C(NC1+KBAS,IB)
     &                            + PAB2*RR(M, 6)*C(NC2+KBAS,IB)
C           (++|+B) = PAB*{(++|B+)} = PAB*{(++|++) + ((++|-+))}
            B2(MLB,4) = B2(MLB,4) + PAB1*RR(M,16)*C(NC1+KBAS,IB)
     &                            + PAB2*RR(M,14)*C(NC2+KBAS,IB)
C           (--|-B) = PAB*{(--|B-)} = PAB*{(--|+-) + ((--|--))}
            B2(MLB,5) = B2(MLB,5) + PAB2*RR(M, 3)*C(NC1+KBAS,IB)
     &                            + PAB1*RR(M, 1)*C(NC2+KBAS,IB)
C           (+-|-B) = PAB*{(+-|B-)} = PAB*{(+-|+-) + ((+-|--))}
            B2(MLB,6) = B2(MLB,6) + PAB2*RR(M,11)*C(NC1+KBAS,IB)
     &                            + PAB1*RR(M, 9)*C(NC2+KBAS,IB)
C           (-+|-B) = PAB*{(-+|B-)} = PAB*{(-+|+-) + ((-+|--))}
            B2(MLB,7) = B2(MLB,7) + PAB2*RR(M, 7)*C(NC1+KBAS,IB)
     &                            + PAB1*RR(M, 5)*C(NC2+KBAS,IB)
C           (++|-B) = PAB*{(++|B-)} = PAB*{(++|+-) + ((++|--))}
            B2(MLB,8) = B2(MLB,8) + PAB2*RR(M,15)*C(NC1+KBAS,IB)
     &                            + PAB1*RR(M,13)*C(NC2+KBAS,IB)
          ENDDO
C
        ENDDO
      ENDDO
C
C     SKIP POINT FOR IQ3 = IQ4
5100  CONTINUE
C
C**********************************************************************C
C     PERFORM ONE-INDEX CONTRACTION:    (IJ;T|KB;T') -> (IJ;T|SB;T')   C
C                                       (IJ;T|LB;T') -> (IJ;T|$B;T')   C
C**********************************************************************C
C
C     DIRECT EVALUATION: (IJ;T|KB;T') -> (IJ;T|SB;T')
      DO IOCCB=1,NOCC
        DO IVRTS=1,NVRT
          IS  = IVRTS+NOCC+NSHIFT
C
C         LIST INDEX FOR THIS IOCCB,IVRTS PAIR
          MSB = (IOCCB-1)*NVRT + IVRTS
C
          DO KBAS=1,NBAS(3)
C
C           LOCATION OF RELEVANT INDEX IN B1
            MKB = (KBAS-1)*NOCC + IOCCB
C
C           (--|SB) = (--|-B) + (--|+B)
            SB(MIJ,MSB,1) = SB(MIJ,MSB,1)
     &                             + B1(MKB,1)*DCONJG(C(NC1+KBAS,IS))
     &                             + B1(MKB,5)*DCONJG(C(NC2+KBAS,IS))
C           (-+|SB) = (-+|-B) + (-+|+B)
            SB(MIJ,MSB,2) = SB(MIJ,MSB,2)
     &                             + B1(MKB,3)*DCONJG(C(NC1+KBAS,IS))
     &                             + B1(MKB,7)*DCONJG(C(NC2+KBAS,IS))
C           (+-|SB) = (+-|-B) + (+-|+B)
            SB(MIJ,MSB,3) = SB(MIJ,MSB,3)
     &                             + B1(MKB,2)*DCONJG(C(NC1+KBAS,IS))
     &                             + B1(MKB,6)*DCONJG(C(NC2+KBAS,IS))
C           (++|SB) = (++|-B) + (++|+B)
            SB(MIJ,MSB,4) = SB(MIJ,MSB,4)
     &                             + B1(MKB,4)*DCONJG(C(NC1+KBAS,IS))
     &                             + B1(MKB,8)*DCONJG(C(NC2+KBAS,IS))
          ENDDO
        ENDDO
      ENDDO

C
C     INTEGRAL SYMMETRY CASE: (IJ;T|LB;T') -> (IJ;T|$B;T')
      IF(IQ3.EQ.IQ4) GOTO 5200
C
      DO IOCCB=1,NOCC
        DO IVRTS=1,NVRT
          IS  = IVRTS+NOCC+NSHIFT
C
C         LIST INDEX FOR THIS IOCCB,IVRTS PAIR
          MSB = (IOCCB-1)*NVRT + IVRTS
          DO LBAS=1,NBAS(4)
C
C           LOCATION OF RELEVANT INDEX IN B1
            MLB = (LBAS-1)*NOCC + IOCCB
C
C           (--|$B) = (--|+B) + (--|-B)
            SB(MIJ,MSB,1) = SB(MIJ,MSB,1)
     &                             + B2(MLB,1)*DCONJG(C(ND1+LBAS,IS))
     &                             + B2(MLB,5)*DCONJG(C(ND2+LBAS,IS))
C           (-+|$B) = (-+|+B) + (-+|-B)
            SB(MIJ,MSB,2) = SB(MIJ,MSB,2)
     &                             + B2(MLB,3)*DCONJG(C(ND1+LBAS,IS))
     &                             + B2(MLB,7)*DCONJG(C(ND2+LBAS,IS))
C           (+-|$B) = (+-|+B) + (+-|-B)
            SB(MIJ,MSB,3) = SB(MIJ,MSB,3)
     &                             + B2(MLB,2)*DCONJG(C(ND1+LBAS,IS))
     &                             + B2(MLB,6)*DCONJG(C(ND2+LBAS,IS))
C           (++|$B) = (++|+B) + (++|-B)
            SB(MIJ,MSB,4) = SB(MIJ,MSB,4)
     &                             + B2(MLB,4)*DCONJG(C(ND1+LBAS,IS))
     &                             + B2(MLB,8)*DCONJG(C(ND2+LBAS,IS))
          ENDDO
        ENDDO
      ENDDO
C
C     SKIP POINT FOR IQ3 = IQ4
5200  CONTINUE
C
C     SKIP POINT FOR INTEGRAL SCREENING
5001  CONTINUE
C     END LOOP OVER IBAS AND JBAS
5000  CONTINUE
C
C**********************************************************************C
C     SB ARRAY NOW CONTAINS THE PARTIALLY-TRANSFORMED LIST: (IJ,SB).   C
C     IT IS A FULL LIST SO THERE IS NO MORE NEED FOR THE '$' NOTATION. C
C**********************************************************************C
C
C**********************************************************************C
C     PERFORM ONE-INDEX CONTRACTION:    (IJ;T|SB;T') -> (IA;T|SB;T')   C
C                                       (JI;T|SB;T') -> (JA;T|SB;T')   C
C**********************************************************************C
C
C     CLEAR INTERMEDIATE ARRAYS
      DO IBAS=1,NBAS(1)
        DO MASB=1,NVRT*NOCC*NOCC
          ASB1(IBAS,MASB,1) = DCMPLX(0.0D0,0.0D0)
          ASB1(IBAS,MASB,2) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     DIRECT EVALUATION: (IJ;T|SB;T') -> (IA;T|SB;T')
      DO IOCCA=1,NOCC
        IA = IOCCA+NSHIFT
        DO IOCCB=1,IOCCA
          DO IVRTS=1,NVRT
C
C           LIST INDEX FOR THIS IVRTS,IOCCB PAIR
            MSB  = (IOCCB-1)*NVRT + IVRTS
C
C           LIST INDEX FOR THIS IOCCA,IVRTS,IOCCB PAIR
            MASB = (IOCCA-1)*NVRT*NOCC + MSB
C
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                MIJ = (IBAS-1)*NBAS(2) + JBAS
C
C               (-A|SB) = (--|SB) + (-+|SB)
                ASB1(IBAS,MASB,1) = ASB1(IBAS,MASB,1)
     &                              +      SB(MIJ,MSB,1)*C(NB1+JBAS,IA)
     &                              +      SB(MIJ,MSB,2)*C(NB2+JBAS,IA)
C               (+A|SB) = (+-|SB) + (++|SB)
                ASB1(IBAS,MASB,2) = ASB1(IBAS,MASB,2)
     &                              +      SB(MIJ,MSB,3)*C(NB1+JBAS,IA)
     &                              +      SB(MIJ,MSB,4)*C(NB2+JBAS,IA)
C
              ENDDO
            ENDDO
          ENDDO
        ENDDO     
      ENDDO
C
C     INTEGRAL SYMMETRY CASE: (JI;T|SB;T') -> (JA;T|SB;T')
      IF(IQ1.EQ.IQ2) GOTO 4100
C
C     CLEAR INTERMEDIATE ARRAYS
      DO JBAS=1,NBAS(2)
        DO MASB=1,NVRT*NOCC*NOCC
          ASB2(JBAS,MASB,1) = DCMPLX(0.0D0,0.0D0)
          ASB2(JBAS,MASB,2) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
      DO IOCCA=1,NOCC
        IA = IOCCA+NSHIFT
        DO IOCCB=1,IOCCA
          DO IVRTS=1,NVRT
C
C           LIST INDEX FOR THIS IVRTS,IOCCB PAIR
            MSB  = (IOCCB-1)*NVRT + IVRTS
C
C           LIST INDEX FOR THIS IOCCA,IVRTS,IOCCB PAIR
            MASB = (IOCCA-1)*NVRT*NOCC + MSB
C
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                MIJ = (IBAS-1)*NBAS(2) + JBAS
C
C               (+A|SB) = PCD*{(+A|SB)} = PCD*{(++|SB) + (-+|SB)}
                ASB2(JBAS,MASB,1) = ASB2(JBAS,MASB,1)
     &                              + PCD1*SB(MIJ,MSB,4)*C(NA1+IBAS,IA)
     &                              + PCD2*SB(MIJ,MSB,2)*C(NA2+IBAS,IA)
C               (-A|SB) = PCD*{(-A|SB)} = PCD*{(+-|SB) + (--|SB)}
                ASB2(JBAS,MASB,2) = ASB2(JBAS,MASB,2)
     &                              + PCD2*SB(MIJ,MSB,3)*C(NA1+IBAS,IA)
     &                              + PCD1*SB(MIJ,MSB,1)*C(NA2+IBAS,IA)
C
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDDO
C
C     SKIP POINT FOR IQ1 = IQ2
4100  CONTINUE
C
C**********************************************************************C
C     PERFORM ONE-INDEX CONTRACTION:    (IA;T|SB;T') -> (RA;T|SB;T')   C
C                                       (JA;T|SB;T') -> (PA;T|SB;T')   C
C**********************************************************************C
C
C     DIRECT EVALUATION: (IA;T|SB;T') -> (RA;T|SB;T')
      M = 0
      DO IOCCA=1,NOCC
        DO IVRTR=1,NVRT
          IR = IVRTR+NOCC+NSHIFT
          DO IOCCB=1,IOCCA
            DO IVRTS=1,NVRT
C
C             LIST INDEX FOR THIS IVRTS,IOCCB PAIR
              MSB  = (IOCCB-1)*NVRT + IVRTS
C
C             LIST INDEX FOR THIS IVRTR,IOCCA,IVRTS,IOCCB COMBINATION
              MRASB = (IOCCA-1)*NVRT*NOCC + MSB
C
              M = M+1
              DO IBAS=1,NBAS(1)
C               (RA|SB) = (-A|SB) + (+A|SB)
                RASB(M) = RASB(M)
     &                    + ASB1(IBAS,MRASB,1)*DCONJG(C(NA1+IBAS,IR))
     &                    + ASB1(IBAS,MRASB,2)*DCONJG(C(NA2+IBAS,IR))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO    
C
C     INTEGRAL SYMMETRY CASE: (JA;T|SB;T') -> (PA;T|SB;T')
      IF(IQ1.EQ.IQ2) GOTO 4200
C
      M = 0
      DO IOCCA=1,NOCC
        DO IVRTR=1,NVRT
          DO IOCCB=1,IOCCA
            DO IVRTS=1,NVRT
C
C             LIST INDEX FOR THIS IVRTS,IOCCB PAIR
              MSB  = (IOCCB-1)*NVRT + IVRTS
C
C             LIST INDEX FOR THIS IVRTR,IOCCA,IVRTS,IOCCB COMBINATION
              MRASB = (IOCCA-1)*NVRT*NOCC + MSB
C
              M = M+1
              DO JBAS=1,NBAS(2)
C               (PA|SB) = (+A|SB) + (-A|SB)
                RASB(M) = RASB(M)
     &                    + ASB2(JBAS,MRASB,1)*DCONJG(C(NB1+JBAS,IR))
     &                    + ASB2(JBAS,MRASB,2)*DCONJG(C(NB2+JBAS,IR))

              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     SKIP POINT FOR IQ1 = IQ2
4200  CONTINUE
C
C**********************************************************************C
C     ALL CONTRIBUTIONS FROM THIS CLASS (A,B,C,D) NOW ACCOUNTED FOR.   C
C     RASB ARRAY CONTAINS FULL LIST, SO DROP THE 'P' NOTATION.         C
C**********************************************************************C
C
4001  CONTINUE
C     END LOOP OVER ALL |MQN| VALUES A, B, C, D
4000  CONTINUE
C     SKIPPING POINT FOR INCLUSION LEVELS
3001  CONTINUE
      WRITE(*,*) 'ICENT',ICNTA,ICNTB,ICNTC,ICNTD
      WRITE(*,*) 'KQNS ',KQN(1),KQN(2),KQN(3),KQN(4)
      WRITE(*,*) '---------------------------------------'
C     END LOOP OVER CENTERS AND KQNS C AND D
3000  CONTINUE
C     END LOOP OVER CENTERS AND KQNS A AND B
2000  CONTINUE
C     END LOOP OVER COMPONENT OVERLAPS
1000  CONTINUE
C
C
C*********************************************************************C
C     CALCULATE SECOND-ORDER PAIR CORRELATION ENERGY                  C
C*********************************************************************C
C
      DO IOCCA=1,NOCC
        DO IOCCB=1,IOCCA
C         EMPTY THE COUNTERS
          DO N=1,6
            EAB2(IOCCA,IOCCB,N) = 0.0D0
          ENDDO
          DO IVRTR=1,NVRT
            DO IVRTS=1,NVRT
C
C             LIST ADDRESS FOR THIS (R,A,S,B) COMBINATION
              MABRS = NVRT*NVRT*IOCCA*(IOCCA-1)/2
     &              + (IVRTR-1)*NVRT*IOCCA + (IOCCB-1)*NVRT+IVRTS
C
C             LIST ADDRESS FOR THIS (R,B,S,A) COMBINATION
              MBARS = NVRT*NVRT*IOCCA*(IOCCA-1)/2
     &              + (IVRTS-1)*NVRT*IOCCA + (IOCCB-1)*NVRT+IVRTR
C
C             NUMERATOR FOR DIRECT AND EXCHANGE CONTRIBUTIONS
              RNUM1 = DREAL(RASB(MABRS)*DCONJG(RASB(MABRS)))
              RNUM2 =-DREAL(RASB(MABRS)*DCONJG(RASB(MBARS)))
C
C             FOCK ADDRESSES
              IA = IOCCA+NSHIFT
              IB = IOCCB+NSHIFT
              IR = IVRTR+NOCC+NSHIFT
              IS = IVRTS+NOCC+NSHIFT
C
C             DENOMINATOR FOR BOTH CONTRIBUTIONS
              DABRS = EIGEN(IA) + EIGEN(IB) - EIGEN(IR) - EIGEN(IS)
C
C             ADD TO DIRECT AND EXCHANGE BINS
              EAB2(IOCCA,IOCCB,4) = EAB2(IOCCA,IOCCB,4) + RNUM1/DABRS
              EAB2(IOCCA,IOCCB,5) = EAB2(IOCCA,IOCCB,5) + RNUM2/DABRS
C
            ENDDO
          ENDDO
          EAB2(IOCCA,IOCCB,6) = EAB2(IOCCA,IOCCB,4)+EAB2(IOCCA,IOCCB,5)
          IF(IOCCA.NE.IOCCB) THEN
            EAB2(IOCCB,IOCCA,4) = EAB2(IOCCA,IOCCB,4)
            EAB2(IOCCB,IOCCA,5) = EAB2(IOCCA,IOCCB,5)
            EAB2(IOCCB,IOCCA,6) = EAB2(IOCCA,IOCCB,6)
          ENDIF
        ENDDO
      ENDDO
C
C**********************************************************************C
C     OUTPUT SECOND-ORDER PAIR CORRELATION ENERGIES AND CALCULATE      C
C     THE TOTAL CORRELATION ENERGY                                     C
C**********************************************************************C
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
C     WRITE PAIR-WISE TERMS
      DO IOCCA=1,NOCC
        DO IOCCB=1,IOCCA
          WRITE(6,21) IOCCA,IOCCB,0.0D0,(EAB2(IOCCA,IOCCB,N),N=4,6)
        ENDDO
      ENDDO
C
C     TOTAL CORRELATION ENERGY
      E2 = 0.0D0
      DO IOCCA=1,NOCC
        DO IOCCB=1,NOCC
          E2 = E2 + 0.5D0*EAB2(IOCCA,IOCCB,6)
        ENDDO
      ENDDO
C    
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
22    FORMAT('Second order correlation energy = ',F13.7)
      WRITE(6,22) E2
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
      RETURN
      END

