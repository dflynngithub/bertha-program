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
      INCLUDE 'param.h'
C
      character*2 klab
      CHARACTER*4 HMLT
      CHARACTER*80 TITLE
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),JQN(4),MQN(4),NBAS(4)
      DIMENSION INDEX(MCT,-(MEL+1):MEL,MKP),ICNT(4)
      
      DIMENSION EB(MKP,MKP,5)
C
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VUEH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),WDIR(MDM,MDM),
     &           WXCH(MDM,MDM),CPLE(MDM,MDM)
C
      COMMON/CMON/ICNTA,ICNTB,ICNTC,ICNTD
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/DENS/DENC,DENO,DENT
      COMMON/EILS/EILSFL(MFL,24),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/I2EL/PAB1,PAB2,PCD1,PCD2,NA1,NB1,NC1,ND1,NA2,NB2,NC2,ND2,
     &            IBAS,JBAS,MCNT,NADDAB,NADDCD,NBAS,IQL,IQR
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS
      COMMON/IRCM/IEAB,IECD,IRIJ(MBS,MBS)
      COMMON/ISCR/IMTX(MB2,11),ISCR(MB2),IMAP(MB2),IBCH,ITOG,MAXN
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VUEH,QDIR,
     &            QXCH,WDIR,WXCH,CPLE
      COMMON/PRMS/HMLT,ICLC,INEW,IATM,ILIN
      COMMON/SHLL/ALPH,BETA,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
      COMMON/SWRZ/GDSC(MDM,MDM),BDSC(MDM,MDM)
      COMMON/T2EL/F2ES(5,6),T2ES(5,6),N2EB(5,6),N2EI(5,6),N2ES(5,6)
      COMMON/TGGL/ILEV,ICB1,ISML,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/TSCF/TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRR,TCC1,TCC2,TCMC,
     &            TB1B,TB1R,TB1F,TB1M,TBEC,TBRM,TBRR,TBC1,TBC2,TBMC,
     &            TQMX,THMX,TC1T,TC2T,TB1T,TB2T,TEIG,TSCR,TTOT,
     &            TC1S,TC2S,TB1S,TB2S

      COMMON/QNMS/LABICN(MDM),LABKQN(MDM),LABMQN(MDM)
C
C     INTEGRAL SCREENING SENSITIVITY PARAMETER
      DATA SENS/1.0D-12/
C
C     DFWED
      DO ITT=1,5
        DO KA=1,MKP
          DO KB=1,MKP
            EB(KA,KB,ITT) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     INTEGRAL SKIPPING ON MOLECULAR GROUP SYMMETRY CLASS BASIS
      IF(IATM.EQ.1.OR.ILIN.EQ.1) THEN
        ISYM = 1
      ELSE
        ISYM = 0
      ENDIF
CC
C      KABKEEP = 1
C      KCDKEEP = 1
Cc
C      do 9999 kab=1,3
CC
C      KQNA = KAPA(kab,1)
CC
C      DO 8888 KCD=1,3
CC
C      KQNB = KAPA(KCD,1)
C
C     INITIALISE STORAGE MATRICES
      DO I=1,NDIM
        DO J=1,NDIM
          BDIR(I,J) = DCMPLX(0.0D0,0.0D0)
          BXCH(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
c
c      OPEN(UNIT=8,FILE='breit/Rn_2019.dat',STATUS='UNKNOWN')
c      REWIND(UNIT=8)
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
            LENIQ               = ICOUNT
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     LOOP OVER ALL NUCLEAR CENTRES (USE INDEX 1000)                   C
C**********************************************************************C
C
C     LOOP OVER CENTRE A
      DO 1000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 1000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C     LOOP OVER CENTRE C
      DO 1000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE C
        XYZ(1,3) = BXYZ(1,ICNTC)
        XYZ(2,3) = BXYZ(2,ICNTC)
        XYZ(3,3) = BXYZ(3,ICNTC)
C
C     LOOP OVER CENTRE D
      DO 1000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE D
        XYZ(1,4) = BXYZ(1,ICNTD)
        XYZ(2,4) = BXYZ(2,ICNTD)
        XYZ(3,4) = BXYZ(3,ICNTD)
C
C     NUMBER OF NUCLEAR CENTRES INVOLVED IN THIS OVERLAP
      MCNT = NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
C
C     SKIP ONE-CENTRE CONTRIBUTIONS (DEFER TO RACAH ALGEBRA ROUTINE)
      IF(MCNT.EQ.1.AND.ICB1.EQ.1) THEN
        GOTO 1001
      ENDIF
C
C     SKIP MULTI-CENTRE CONTRIBUTIONS IN STAGE 1
      IF(MCNT.NE.1.AND.ILEV.EQ.1) GOTO 1001
C
C**********************************************************************C
C     LOOP OVER ALL KQN SYMMETRY TYPES (USE INDEX 2000)                C
C**********************************************************************C
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
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
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
C     LOOP OVER KQN(C) VALUES
      DO 2000 KC=1,NKAP(ICNTC)
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
C     LOOP OVER KQN(D) VALUES
      DO 2000 KD=1,NKAP(ICNTD)
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
C     THIS UNIQUELY DEFINES A FULL SET OF RC(AB|CD) INTEGRALS -- RESET
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          IRIJ(IBAS,JBAS) = 1
        ENDDO
      ENDDO
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON KQN                           C
C**********************************************************************C
C
C     ATOM-CENTRED SELECTION RULES (ONLY APPLIES IF ICB1 SWITCHED OFF)
      IF(MCNT.EQ.1) THEN
C
C       LQN PAIR PARITY
        IPARAB = MOD(LQN(1)+LQN(2)+1,2)
        IPARCD = MOD(LQN(3)+LQN(4)+1,2)
C
C       LQN PAIR PARITY SELECTION RULE
        IF(IPARAB.NE.IPARCD) THEN
          GOTO 2001
        ENDIF
C
C       JQN TRIANGLE RULE CHECK FOR MULTIPOLE EXPANSION (ATOM-CENTRED)
        NUI = MAX0(IABS(JQN(1)-JQN(2))/2,IABS(JQN(3)-JQN(4))/2)
        NUF = MIN0(    (JQN(1)+JQN(2))/2,    (JQN(3)+JQN(4))/2)
        IF(NUI.GT.NUF) THEN
          GOTO 2001
        ENDIF
C
C       ADDITIONAL LQN SELECTION RULE PARITY ANALYSIS
        ISELK = 0
        DO NU=NUI,NUF
C
C         A AND B: LQN(1)+LQN(2)+NU EVEN OR ODD
          IPARAB = MOD(LQN(1)+LQN(2)+1+NU,2)
          IPARCD = MOD(LQN(3)+LQN(4)+1+NU,2)
C
C         CASE 1: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH ODD (UNLESS NU=0)
          IF(IPARAB.EQ.0.AND.IPARCD.EQ.0.AND.NU.NE.0) ISELK = 1
C
C         CASE 2: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH EVEN
          IF(IPARAB.EQ.1.AND.IPARCD.EQ.1) ISELK = 1
C
        ENDDO
        IF(ISELK.EQ.0) GOTO 2001
C
      ENDIF
CC
CC     ADDITIONAL SELECTION RULES FOR CLOSED-SHELL ATOMS
C      IF(MCNT.EQ.1) THEN
C        iselk = 0
C        if(ka.eq.kb.and.kc.eq.kd) iselk = 1
C        if(ka.eq.kc.and.kb.eq.kd) iselk = 1
C        if(ka.eq.kd.and.kb.eq.kc) iselk = 1
C        if(iselk.eq.0) goto 2001
C      ENDIF
CC      
C      ISELK=0
CC     IF(KA.EQ.KAB.AND.KB.EQ.KAB.AND.KC.EQ.KCD.AND.KD.EQ.KCD) ISELK=1
C      IF(KA.EQ.KAB.AND.KC.EQ.KAB.AND.KB.EQ.KCD.AND.KD.EQ.KCD) ISELK=1
C      IF(KA.EQ.KAB.AND.KD.EQ.KAB.AND.KC.EQ.KCD.AND.KB.EQ.KCD) ISELK=1
CC      IF(KB.EQ.KAB.AND.KC.EQ.KAB) ISELK=1
CC      IF(KA.EQ.KAB.AND.KD.EQ.KAB) ISELK=1
CC      IF(KB.EQ.KAB.AND.KD.EQ.KAB) ISELK=1
C      IF(ISELK.EQ.0) GOTO 2001
CC
C15    FORMAT(1X,'( ',A,1X,A,' || ',A,1X,A,' )')
c      WRITE(*,15) (KLAB(KQN(N)),N=1,4)
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
      IF(ISYM.EQ.1) THEN
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) GOTO 3003
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) GOTO 3003
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) GOTO 3003
        GOTO 3001
      ENDIF
3003  CONTINUE
C
C     ATOM-CENTRED SELECTION RULES (ONLY APPLIES IF ICB1 SWITCHED OFF)
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
        IF(ISELM.EQ.0) GOTO 3001
      ENDIF
CC
CC     ADDITIONAL SELECTION RULES FOR CLOSED-SHELL ATOMS
C      IF(MCNT.EQ.1) THEN
C        iselm = 0
Cc       if(ma.eq.mb.and.mc.eq.md) iselm = 1
C        if(ma.eq.mc.and.mb.eq.md) iselm = 1
C        if(ma.eq.md.and.mb.eq.mc) iselm = 1
C        if(iselm.eq.0) goto 3001
C      ENDIF
CC      
C      if(ma+mb+mc+md.gt.4) goto 3001
C
C      WRITE(*,*) KLAB(KQN(1)),KLAB(KQN(2)),KLAB(KQN(3)),KLAB(KQN(4)),
C     &            MQN(1),MQN(2),MQN(3),MQN(4)
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
C     COMBINED BLOCK INDEX IN A TWO-FUNCTION LIST
      IQL = (IQ1-1)*LENIQ + IQ2
      IQR = (IQ3-1)*LENIQ + IQ4
C
C     SKIP CONTRIBUTIONS THAT ARISE BY PERMUTATION OF INTEGRALS
      IF(IQL.LT.IQR) GOTO 3001
C
C     DFWED: TRY TO FIND A RELATION SO THIS CAN BE USED
C     SKIP CONTRIBUTIONS THAT ARISE BY PERMUTATION OF INTEGRALS
C      IF(IQ1.LT.IQ2) GOTO 3001
C      IF(IQ3.LT.IQ4) GOTO 3001
C
C     EQ-COEFFICIENT PHASE FACTORS FOR PERMUTATION OF R-INTEGRALS
C     NB! OPPOSITE PHASE AS IN THE LL/SS CASE SEEN IN SCF
      PAB1 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PAB2 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      PCD1 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PCD2 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
C     FURTHER DEFINE STARTING ADDRESSES FOR {ABCD} BASIS OVERLAPS
      NA1 = LRGE(ICNTA,KA,2*MA-1)
      NA2 = LRGE(ICNTA,KA,2*MA  )
      NB1 = LRGE(ICNTB,KB,2*MB-1)+NSKP
      NB2 = LRGE(ICNTB,KB,2*MB  )+NSKP
      NC1 = LRGE(ICNTC,KC,2*MC-1)
      NC2 = LRGE(ICNTC,KC,2*MC  )
      ND1 = LRGE(ICNTD,KD,2*MD-1)+NSKP
      ND2 = LRGE(ICNTD,KD,2*MD  )+NSKP
C
C     UPDATE COUNTER FOR NUMBER OF BLOCKS
      N2EB(MCNT,ITT) = N2EB(MCNT,ITT)+1
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS (IBAS,JBAS) TO CONSTRUCT WDIR/WMAT     C
C**********************************************************************C
C
C     RECORD TIME AT START OF BATCH
      CALL CPU_TIME(TI)
C
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
C
C         SCHWARZ SCREENING (ECONOMIC ONLY WHEN SCREENING FRACTION BIG)
          IF(ISWZ.EQ.0) THEN
            ITOG = 0
          ENDIF
C
C         CALL THE SCREENING ROUTINE
          CALL SCHWARZ(BDSC,SENS,TB2S)
C
C         UPDATE COUNTER FOR NUMBER OF INTEGRALS AND SCREENED INTEGRALS
          N2EI(MCNT,ITT) = N2EI(MCNT,ITT) + NBAS(3)*NBAS(4)
          N2ES(MCNT,ITT) = N2ES(MCNT,ITT) + NBAS(3)*NBAS(4)-MAXN
C
C         CONDITIONAL TO SKIP THIS BATCH
          IF(IBCH.EQ.1) THEN
C
C           GENERATE A BATCH OF BREIT INTERACTION INTEGRALS
            CALL BII(RR,XYZ,KQN,MQN,NBAS,EXL,IBAS,JBAS)
C
C           MULTIPLY BY DENSITY ELEMENTS AND ADD TO WDIR/WMAT
            CALL BRTMAT(RR,TBMC)
C
          ENDIF
C
        ENDDO
      ENDDO
C
C     RECORD TIME AT END OF BATCH
      CALL CPU_TIME(TF)
      T2ES(MCNT,ITT) = T2ES(MCNT,ITT) + TF - TI
C
3001  CONTINUE
3000  CONTINUE
2001  CONTINUE
2000  CONTINUE
1001  CONTINUE
1000  CONTINUE
C
C**********************************************************************C
C     COMPLETE CONSTRUCTION OF BREIT MATRIX BY MATRIX CONJUGATION.     C
C**********************************************************************C
C
      DO I=1,NDIM-NSKP
        DO J=NSKP+1,NDIM
          BDIR(J,I) = DCONJG(BDIR(I,J))
          BXCH(J,I) = DCONJG(BXCH(I,J))
          WDIR(J,I) = DCONJG(WDIR(I,J))
          WXCH(J,I) = DCONJG(WXCH(I,J))
        ENDDO
      ENDDO
CC
CC     CALCULATION OF BREIT ENERGY
C      ELLSS = 0.0D0
C      ESSLL = 0.0D0
C      ESLLS = 0.0D0
C      ELSSL = 0.0D0
C      DO IL=1,NSKP
C        DO JL=1,NSKP
CC
CC         SMALL-COMPONENT LABELS
C          IS = IL+NSKP
C          JS = JL+NSKP
CC
C          ELLSS = ELLSS - 0.5D0*BXCH(IL,JL)*DENT(IL,JL)
C          ESSLL = ESSLL - 0.5D0*BXCH(IS,JS)*DENT(IS,JS)
C          ESLLS = ESLLS - 0.5D0*BXCH(IS,JL)*DENT(IS,JL)
C          ELSSL = ELSSL - 0.5D0*BXCH(IL,JS)*DENT(IL,JS)
CC
C        ENDDO
C      ENDDO
C      EBTOT = ELLSS+ESSLL+ESLLS+ELSSL

      GOTO 111
C
C     CALCULATION OF BREIT ENERGY
      DO IL=1,NSKP
        DO JL=1,NSKP
C
C         SMALL-COMPONENT LABELS
          IS = IL+NSKP
          JS = JL+NSKP
C
          EB(KAB,KCD,1) = EB(KAB,KCD,1)
     &                  - 0.5D0*DREAL(BXCH(IL,JL)*DENT(IL,JL))
          EB(KAB,KCD,2) = EB(KAB,KCD,2)
     &                  - 0.5D0*DREAL(BXCH(IL,JS)*DENT(IL,JS))
          EB(KAB,KCD,3) = EB(KAB,KCD,3)
     &                  - 0.5D0*DREAL(BXCH(IS,JL)*DENT(IS,JL))
          EB(KAB,KCD,4) = EB(KAB,KCD,4)
     &                  - 0.5D0*DREAL(BXCH(IS,JS)*DENT(IS,JS))
C
        ENDDO
      ENDDO
C
C11    FORMAT(A,1X,F18.12)
C      WRITE(*,11) 'E(LL|SS) = ',ELLSS
C      WRITE(*,11) 'E(LS|LS) = ',ELSSL
C      WRITE(*,11) 'E(SL|LS) = ',ESLLS
C      WRITE(*,11) 'E(SS|LL) = ',ESSLL
C      WRITE(*,11) 'TOTAL    = ',ELLSS+ESSLL+ESLLS+ELSSL
C      WRITE(*,*) REPEAT('-',72)
C
8888  CONTINUE
9999  continue
c      CLOSE(UNIT=8)
12    format(1X,'(',i2,',',i2,')',2x,5(f16.10))
      ELLTOT = 0.0D0
      ELSTOT = 0.0D0
      ESLTOT = 0.0D0
      ESSTOT = 0.0D0
      DO KAB=1,3
        DO KCD=1,3
          EB(KAB,KCD,5) = EB(KAB,KCD,1) + EB(KAB,KCD,2) 
     &                  + EB(KAB,KCD,3) + EB(KAB,KCD,4)
          write(*,12) KAPA(KAB,1),KAPA(KCD,1),(EB(KAB,KCD,N),N=1,5)
          ELLTOT = ELLTOT + EB(KAB,KCD,1)
          ELSTOT = ELSTOT + EB(KAB,KCD,2)
          ESLTOT = ESSTOT + EB(KAB,KCD,3)
          ESSTOT = ESSTOT + EB(KAB,KCD,4)
        ENDDO
      ENDDO
      ETTTOT = ELLTOT + ELSTOT + ESLTOT + ESSTOT
C
13    format(A,2x,5(f16.10))
      WRITE(*,*) REPEAT('-',73)
      write(*,13) ' Atomic ',ELLTOT,ELSTOT,ESLTOT,ESSTOT,ETTTOT
      WRITE(*,*) REPEAT('-',73)
111   CONTINUE

c      STOP
C
C     MULTIPLY OPEN MATRIX BY ANGULAR COEFFICIENTS (LIFTED FROM COULOMB)
C      IF(NOPN.GT.0) THEN
C        DO J=1,NDIM
C          DO I=1,NDIM
C            WDIR(I,J) = ACFF*WDIR(I,J)
C            WXCH(I,J) = BCFF*WXCH(I,J)
C          ENDDO
C        ENDDO
C      ENDIF

c      stop
C
      RETURN
      END

