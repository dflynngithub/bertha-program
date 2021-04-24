      SUBROUTINE BREIT1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             BBBBBBB  RRRRRRR  EEEEEEEE IIII TTTTTTTT 11              C
C             BB    BB RR    RR EE        II     TT   111              C
C             BB    BB RR    RR EE        II     TT    11              C
C             BBBBBBB  RR    RR EEEEEE    II     TT    11              C
C             BB    BB RRRRRRR  EE        II     TT    11              C
C             BB    BB RR    RR EE        II     TT    11              C
C             BBBBBBB  RR    RR EEEEEEEE IIII    TT   1111             C
C                                                                      C
C -------------------------------------------------------------------- C
C  BREIT1 CONSTRUCTS ALL ONE-CENTER CONTRIBUTIONS TO THE MOLECULAR     C
C  MEAN-FIELD BREIT MATRIX WITH RACAH ALGEBRA AND BETA INTEGRALS.      C
C  THIS ROUTINE IS SIMILAR TO THE IN-LINE CONSTRUCTION OF MEAN-FIELD   C
C  ATOMIC BREIT MATRIX IN HFSCF0, BUT WITH MQN STRUCTURE AS WELL.      C
C -------------------------------------------------------------------- C
C  DFNOTE: NOT WORKING YET...                                          C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MNU=MKP+1,
     &                                                      MAB=2*MNU+6)
C
      CHARACTER*4 HMLTN
C
      DIMENSION EXPT(MBS,4),KQN(4),MQN(4),NBAS(4),LQN(4)
      DIMENSION DKAB(MNU,MKP+1,MKP+1),DKCD(MNU,MKP+1,MKP+1)
      DIMENSION EMAT(MNU,8),ANGFAC(8)
      DIMENSION XLSLS(MB2),XSLLS(MB2),XLSSL(MB2),XSLSL(MB2)
C
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           GDIR(MDM,MDM),GXCH(MDM,MDM),QDIR(MDM,MDM),
     &           QXCH(MDM,MDM),BDIR(MDM,MDM),BXCH(MDM,MDM),
     &           VUEH(MDM,MDM),FOCK(MDM,MDM)
C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SHLL/ACFF,BCFF,FOPN,ICLS(500),IOPN(6),NCLS,NOPN,NOELEC
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
      COMMON/T2EL/F2ES(5,6),T2ES(5,6),N2EB(5,6),N2EI(5,6),N2ES(5,6)
      COMMON/TSCF/TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRR,TCC1,TCC2,TCMC,
     &            TB1B,TB1R,TB1F,TB1M,TBEC,TBRM,TBRR,TBC1,TBC2,TBMC,
     &            THMX,TC1T,TC2T,TB1T,TB2T,TEIG,TSCR,TTOT,
     &            TC1S,TC2S,TB1S,TB2S
      COMMON/XBIJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ
      COMMON/XCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
      COMMON/XRKB/RKLSLS(MB2,MNU,2),RKSLSL(MB2,MNU,2),
     &            RKSLLS(MB2,MNU,2),RKLSSL(MB2,MNU,2)
C
      DATA SENS/1.0D-12/
C
      CALL CPU_TIME(TBCH1)
C
C**********************************************************************C
C     LOOP OVER ATOMIC CENTERS (USE INDEX 1000)                        C
C**********************************************************************C
C
C     LOOP OVER ALL CENTERS
      DO 1000 ICNT=1,NCNT
C
C**********************************************************************C
C     LOOP OVER KQN SYMMETRY TYPE FOR A,B,C,D BLOCKS (USE INDEX 2000)  C
C**********************************************************************C
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNT)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KVALS(KA,ICNT)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFUNCT(LQN(1)+1,ICNT)
        DO IBAS=1,NBAS(1)
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNT)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNT)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KVALS(KB,ICNT)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFUNCT(LQN(2)+1,ICNT)
        DO JBAS=1,NBAS(2)
          EXPT(JBAS,2) = EXPSET(JBAS,LQN(2)+1,ICNT)
        ENDDO
C
C     LOOP OVER KQN(C) VALUES
      DO 2000 KC=1,NKAP(ICNT)
C
C       QUANTUM NUMBERS FOR BLOCK C
        KQN(3) = KVALS(KC,ICNT)
        IF(KQN(3).GT.0) THEN
          LQN(3) = KQN(3)
        ELSE
          LQN(3) =-KQN(3)-1
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK C
        NBAS(3) = NFUNCT(LQN(3)+1,ICNT)
        DO KBAS=1,NBAS(3)
          EXPT(KBAS,3) = EXPSET(KBAS,LQN(3)+1,ICNT)
        ENDDO
C
C     LOOP OVER KQN(D) VALUES
      DO 2000 KD=1,NKAP(ICNT)
C
C       QUANTUM NUMBERS FOR BLOCK D
        KQN(4) = KVALS(KD,ICNT)
        IF(KQN(4).GT.0) THEN
          LQN(4) = KQN(4)
        ELSE
          LQN(4) =-KQN(4)-1
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK D
        NBAS(4) = NFUNCT(LQN(4)+1,ICNT)
        DO LBAS=1,NBAS(4)
          EXPT(LBAS,4) = EXPSET(LBAS,LQN(4)+1,ICNT)
        ENDDO
C
C**********************************************************************C
C     PREPARE INTERMEDIATE DATA FOR USE IN RKBRT1                      C
C**********************************************************************C
C
C     COEFFICIENTS FOR INTEGRAL ASSEMBLY
      C3 = 0.25D0*GAMHLF(LQN(1)+LQN(2)+LQN(3)+LQN(4)+3)
      C5 = 0.25D0*GAMHLF(LQN(1)+LQN(2)+LQN(3)+LQN(4)+5)
      C7 = 0.25D0*GAMHLF(LQN(1)+LQN(2)+LQN(3)+LQN(4)+7)
C
      V1 = 1.0D0
      V2 = 2.0D0
      V4 = 4.0D0
C
      TI = DFLOAT(LQN(1)+KQN(1)+1)
      TJ = DFLOAT(LQN(2)+KQN(2)+1)
      TK = DFLOAT(LQN(3)+KQN(3)+1)
      TL = DFLOAT(LQN(4)+KQN(4)+1)
C
      T0000 = 1.0D0
      T1000 = TI
      T0100 = TJ
      T0010 = TK
      T0001 = TL
      T1010 = TI*TK
      T1001 = TI*TL
      T0110 = TJ*TK
      T0101 = TJ*TL
C
C     ANGULAR COEFFICIENTS
      CALL CPU_TIME(T1)
      CALL ANGBRT1(DKAB,DKCD,EMAT,KQN,LQN,NUS,NUF,NUI,NUNUM,ISEL)
      CALL CPU_TIME(T2)
      TB1B = TB1B+T2-T1
C
C     EXIT THIS COMBINATION IF IT VIOLATES A SELECTION RULE
      IF(ISEL.EQ.0) GOTO 2001
C
C     BASIS SET INTERMEDIATES FOR THE KL-PAIRS
      CALL CPU_TIME(T1)
      CALL PREPKL(LQN,KC,NBAS(3),KD,NBAS(4),ICNT)
      CALL CPU_TIME(T2)
      TB1B = TB1B+T2-T1
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS IN BLOCKS A AND B (INDEX 3000)         C
C**********************************************************************C
C
C     UPDATE COUNTER FOR NUMBER OF CLASSES
      N2EB(1,5) = N2EB(1,5) + 1
C
      DO 3000 IBAS=1,NBAS(1)
      DO 3000 JBAS=1,NBAS(2)
C
C       UPDATE COUNTER FOR NUMBER OF INTEGRALS
        N2EI(1,5) = N2EI(1,5)+NBAS(3)*NBAS(4)
C
C       GAUSSIAN EXPONENTS FOR THIS PAIR
        EI = EXPT(IBAS,1)
        EJ = EXPT(JBAS,2)
C
C       BASIS SET INTERMEDIATES FOR THE IJ-PAIRS
        CALL CPU_TIME(T1)
        CALL PREPIJ(LQN)
        CALL CPU_TIME(T2)
        TB1B = TB1B+T2-T1
C
C       BATCH OF RADIAL INTEGRALS (EFFECTIVE INTERACTION STRENGTHS)
        CALL CPU_TIME(T1)
        CALL RKBRT1(KQN,LQN,NBAS(3)*NBAS(4))
        CALL CPU_TIME(T2)
        TB1R = TB1R+T2-T1
C
C**********************************************************************C
C     LOOP OVER |KQN| MAGNITUDES FOR A,B,C,D BLOCKS (USE INDEX 4000)   C
C**********************************************************************C
C
      DO 4000 MA=1,IABS(KQN(1))
        MQN(1) = 2*MA-1
C
      DO 4000 MB=1,IABS(KQN(2))
        MQN(2) = 2*MB-1
C
      DO 4000 MC=1,IABS(KQN(3))
        MQN(3) = 2*MC-1
C
      DO 4000 MD=1,IABS(KQN(4))
        MQN(4) = 2*MD-1
C
C     TRANSFORM THE X(L) INTO G-SPINOR MATRIX ELEMENTS USING
C     THE TENSOR EXPANSION IN {L,Q}
C
C     LOOP OVER THE SIGNS OF |MQN| AND DETERMINE FOCK ADDRESSES
      DO 5000 ISGN1=1,2
        MMJA = MQN(1)*((-1)**ISGN1)
        IMJA = MQN(1)+ISGN1-1
C
      DO 5000 ISGN2=1,2
        MMJB = MQN(2)*((-1)**ISGN2)
        IMJB = MQN(2)+ISGN2-1
C
      DO 5000 ISGN3=1,2
        MMJC = MQN(3)*((-1)**ISGN3)
        IMJC = MQN(3)+ISGN3-1
C
      DO 5000 ISGN4=1,2
        MMJD = MQN(4)*((-1)**ISGN4)
        IMJD = MQN(4)+ISGN4-1
C
C     STARTING FOCK ADDRESS FOR EACH BASIS LIST
      NAL = LARGE(ICNT,KA,IMJA)
      NBL = LARGE(ICNT,KB,IMJB)
      NCL = LARGE(ICNT,KC,IMJC)
      NDL = LARGE(ICNT,KD,IMJD)
C
      NAS = LARGE(ICNT,KA,IMJA) + NSHIFT
      NBS = LARGE(ICNT,KB,IMJB) + NSHIFT
      NCS = LARGE(ICNT,KC,IMJC) + NSHIFT
      NDS = LARGE(ICNT,KD,IMJD) + NSHIFT
C
C     APPLY ANGULAR MQN SELECTION RULE
      IF(MMJA-MMJB.NE.MMJD-MMJC) GOTO 5001
C
C     RESET CONTRACTED RADIAL ARRAYS
      CALL CPU_TIME(T1)
      DO M=1,NBAS(3)*NBAS(4)
        XLSLS(M) = 0.0D0
        XLSSL(M) = 0.0D0
        XSLLS(M) = 0.0D0
        XSLSL(M) = 0.0D0
      ENDDO
C
C     ASSEMBLE INTEGRAL BY SUMMING OVER EXPANSION POWER L
      DO LTEN=1,NUNUM
C
C       ANGULAR COEFFICIENTS
        DKABCD = DKAB(LTEN,IMJA,IMJB)*DKCD(LTEN,IMJC,IMJD)
        DO IMU=1,8
          ANGFAC(IMU) = DKABCD*EMAT(LTEN,IMU)
        ENDDO
C
C       SCREENING OF INTEGRAL BASED ON ANGULAR COEFFICIENT
        ANGSUM = 0.0D0
        DO IMU=1,8
          ANGSUM = ANGSUM + ANGFAC(IMU)
        ENDDO
        IF(DABS(ANGSUM).LE.SENS) GOTO 5003
C
        DO M=1,NBAS(3)*NBAS(4)
          XLSLS(M) = XLSLS(M) + ANGFAC(1)*RKLSLS(M,LTEN,1)
     &                        + ANGFAC(2)*RKLSLS(M,LTEN,2)
          XSLSL(M) = XSLSL(M) + ANGFAC(3)*RKSLSL(M,LTEN,1)
     &                        + ANGFAC(4)*RKSLSL(M,LTEN,2)
          XLSSL(M) = XLSSL(M) + ANGFAC(5)*RKLSSL(M,LTEN,1)
     &                        + ANGFAC(6)*RKLSSL(M,LTEN,2)
          XSLLS(M) = XSLLS(M) + ANGFAC(7)*RKSLLS(M,LTEN,1)
     &                        + ANGFAC(8)*RKSLLS(M,LTEN,2)
        ENDDO
C
C       SKIP POINT FOR ANGULAR SCREENING
5003    CONTINUE
C
      ENDDO
C
C     FULL-INTEGRAL CONSTRUCTION COMPLETE
      CALL CPU_TIME(T2)
      TB1F = TB1F+T2-T1
C
C**********************************************************************C
C     ADD THIS BATCH OF R-INTEGRALS TO CLOSED-SHELL COULOMB MATRIX     C
C -------------------------------------------------------------------- C
C  DFNOTE: NO BDIR CONTRIBUTIONS FOR CLOSED-SHELL SYSTEMS.             C
C**********************************************************************C
C
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
C
C         DIRECT CONTRIBUTIONS
          BDIR(NAL+IBAS,NBS+JBAS) = BDIR(NAL+IBAS,NBS+JBAS)
     &                          +      XLSLS(M)*DENC(NCL+KBAS,NDS+LBAS)
     &                          +      XLSSL(M)*DENC(NCS+KBAS,NDL+LBAS)
C
          BDIR(NAS+IBAS,NBL+JBAS) = BDIR(NAS+IBAS,NBL+JBAS)
     &                          +      XSLLS(M)*DENC(NCL+KBAS,NDS+LBAS)
     &                          +      XSLSL(M)*DENC(NCS+KBAS,NDL+LBAS)
C
C
C         EXCHANGE CONTRIBUTIONS
          BXCH(NAL+IBAS,NDL+LBAS) = BXCH(NAL+IBAS,NDL+LBAS)
     &                          +      XLSLS(M)*DENC(NCS+KBAS,NBS+JBAS)
C
          BXCH(NAL+IBAS,NDS+LBAS) = BXCH(NAL+IBAS,NDS+LBAS)
     &                          +      XLSSL(M)*DENC(NCS+KBAS,NBL+JBAS)
C
          BXCH(NAS+IBAS,NDL+LBAS) = BXCH(NAS+IBAS,NDL+LBAS)
     &                          +      XSLLS(M)*DENC(NCL+KBAS,NBS+JBAS)

          BXCH(NAS+IBAS,NDS+LBAS) = BXCH(NAS+IBAS,NDS+LBAS)
     &                          +      XSLSL(M)*DENC(NCL+KBAS,NBL+JBAS)
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     ADD THIS BATCH OF R-INTEGRALS TO OPEN-SHELL SCF MATRIX           C
C**********************************************************************C
C
      IF(NOPN.EQ.0) GOTO 6100
C
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
C
C         DIRECT CONTRIBUTIONS
          QDIR(NAL+IBAS,NBS+JBAS) = QDIR(NAL+IBAS,NBS+JBAS)
     &                          + ACFF*XLSLS(M)*DENO(NCL+KBAS,NDS+LBAS)
     &                          + ACFF*XLSSL(M)*DENO(NCS+KBAS,NDL+LBAS)
C
          QDIR(NAS+IBAS,NBL+JBAS) = QDIR(NAS+IBAS,NBL+JBAS)
     &                          + ACFF*XSLLS(M)*DENO(NCL+KBAS,NDS+LBAS)
     &                          + ACFF*XSLSL(M)*DENO(NCS+KBAS,NDL+LBAS)
C
C         EXCHANGE CONTRIBUTIONS
          QXCH(NAL+IBAS,NDL+LBAS) = QXCH(NAL+IBAS,NDL+LBAS)
     &                          + BCFF*XLSLS(M)*DENO(NCS+KBAS,NBS+JBAS)
C
          QXCH(NAL+IBAS,NDS+LBAS) = QXCH(NAL+IBAS,NDS+LBAS)
     &                          + BCFF*XLSSL(M)*DENO(NCS+KBAS,NBL+JBAS)
C
          QXCH(NAS+IBAS,NDL+LBAS) = QXCH(NAS+IBAS,NDL+LBAS)
     &                          + BCFF*XSLLS(M)*DENO(NCL+KBAS,NBS+JBAS)

          QXCH(NAS+IBAS,NDS+LBAS) = QXCH(NAS+IBAS,NDS+LBAS)
     &                          + BCFF*XSLSL(M)*DENO(NCL+KBAS,NBL+JBAS)
C
        ENDDO
      ENDDO
C
6100  CONTINUE
C
C     MATRIX MULTIPLICATION STEP COMPLETE
      CALL CPU_TIME(T3)
      TB1M = TB1M+T3-T2
C
C**********************************************************************C
C     COMPLETE CONSTRUCTION OF ALL ONE-CENTER CONTRIBUTIONS.           C
C**********************************************************************C
C
C     EARLY EXIT FOR MQN SELECTION RULES
5001  CONTINUE
C     END LOOP OVER ALL |MQN| SIGNS
5000  CONTINUE
C     END LOOP OVER ALL |MQN| MAGNITUDES
4000  CONTINUE
C     END LOOP OVER IBAS AND JBAS
3000  CONTINUE
C     EARLY EXIT FOR KQN SELECTION RULES
2001  CONTINUE
C     END LOOP OVER ALL KQNS
2000  CONTINUE
C     END LOOP OVER ATOMIC CENTERS
1000  CONTINUE
C
C     CALCULATION OF BREIT ENERGY
      ELLSS = 0.0D0
      ESSLL = 0.0D0
      ESLLS = 0.0D0
      ELSSL = 0.0D0
      DO IL=1,NSHIFT
        DO JL=1,NSHIFT
C
C         SMALL-COMPONENT LABELS
          IS = IL+NSHIFT
          JS = JL+NSHIFT
C
          ELLSS = ELLSS - 0.5D0*BXCH(IL,JL)*DENC(IL,JL)
          ESSLL = ESSLL - 0.5D0*BXCH(IS,JS)*DENC(IS,JS)
          ESLLS = ESLLS - 0.5D0*BXCH(IS,JL)*DENC(IS,JL)
          ELSSL = ELSSL - 0.5D0*BXCH(IL,JS)*DENC(IL,JS)
C
        ENDDO
      ENDDO
C
      WRITE(*,*) 'E(LL|SS) = ',ELLSS
      WRITE(*,*) 'E(SS|LL) = ',ESSLL
      WRITE(*,*) 'E(SL|LS) = ',ESLLS
      WRITE(*,*) 'E(LS|SL) = ',ELSSL
      WRITE(*,*) ' '
      WRITE(*,*) 'E_total  = ',ELLSS+ESSLL+ESLLS+ELSSL
      WRITE(*,*) 'E_correct= ',6.38264D-5
      WRITE(*,*) '|Delta|  = ',DABS(6.38264D-5-ELLSS-ESSLL-ESLLS-ELSSL)
      WRITE(*,*) REPEAT('-',72)
C
C      DO I=1,NDIM
C        DO J=1,NDIM
C          ETMP( 6) = ETMP( 6) + 0.5D0*DENT(I,J)*BDIR(I,J)
C          ETMP( 7) = ETMP( 7) + 0.5D0*DENT(I,J)*BXCH(I,J)
C        ENDDO
C      ENDDO

C
C     RECORD CPU TIME AT END OF BATCH AND ADD TO APPROPRIATE COUNTER
      CALL CPU_TIME(TBCH2)
      T2ES(1,5) = T2ES(1,5)+TBCH2-TBCH1
C
      RETURN
      END
C
C
      SUBROUTINE RKBRT1(KQN,LQN,MAXM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            RRRRRRR  KK    KK BBBBBBB  RRRRRRR TTTTTTTT 11            C
C            RR    RR KK   KK  BB    BB RR    RR   TT   111            C
C            RR    RR KK  KK   BB    BB RR    RR   TT    11            C
C            RR    RR KKKKK    BBBBBBB  RR    RR   TT    11            C
C            RRRRRRR  KK  KK   BB    BB RRRRRRR    TT    11            C
C            RR    RR KK   KK  BB    BB RR    RR   TT    11            C
C            RR    RR KK    KK BBBBBBB  RR    RR   TT   1111           C
C                                                                      C
C -------------------------------------------------------------------- C
C  RKBRT1 EVALUATES A BATCH OF GENERAL (OPEN-SHELL) ONE-CENTER RADIAL  C
C  INTEGRALS OVER THE BREIT INTERACTION, SEPARATING RESULTS BY THE     C
C  ALLOWED TENSOR ORDERS RK(ABCD).                                     C
C -------------------------------------------------------------------- C
C  DFNOTE: NOT WORKING YET...                                          C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MNU=MKP+1,MAB=2*MNU+6)
C
      DIMENSION KQN(4),LQN(4)
C
      DIMENSION XJ(MB2,2),IAA(2),IBB(2)
      DIMENSION BETA(MB2),BTA1(MB2),XROOT(MB2),BIN(MB2),TRM(MB2)
      DIMENSION BDU(MB2,-MAB:MAB,-MAB:MAB),BDL(MB2,-MAB:MAB,-MAB:MAB)
C
      COMMON/XBIJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ
      COMMON/XBKL/EKL(MB2,-MAB:MAB),RNKL(MB2,4),EK(MB2),EL(MB2)
      COMMON/XCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
      COMMON/XRKB/RKLSLS(MB2,MNU,2),RKSLSL(MB2,MNU,2),
     &            RKSLLS(MB2,MNU,2),RKLSSL(MB2,MNU,2)
C
C**********************************************************************C
C     PREPARATION OF EXPONENT POWERS AND BETA FUNCTION ARGUMENTS.      C
C**********************************************************************C
C
C     TENSOR ORDER LIMITS
      NUI = NUS(1)
      NUF = NUS(NUNUM)
C
C     GAUSSIAN EXPONENT FOR (IBAS,JBAS) OVERLAP
      TIJ0 = EI+EJ
C
C     INITIALISE THE ARRAY XJ(M,2) FOR INCOMPLETE BETA FUNCTION ARGS.
      DO M=1,MAXM
        TKL0    = EK(M)+EL(M)
        TIJKL   = TIJ0+TKL0
        XJ(M,1) = TIJ0/TIJKL
        XJ(M,2) = TKL0/TIJKL
      ENDDO
C
C**********************************************************************C
C     INLINE BETA FUNCTION CODE FOR DIRECT TERMS                       C
C**********************************************************************C
C
C     NUMBER OF LEVELS NEEDED TO ACCOUNT FOR ALL EXCHANGE INTEGRALS
      NVALS = (NUF-NUI)/2+2
C
C     LOOP OVER EXPANSION TERMINALS FOR FIRST PAIR
      DO NX=1,NVALS
        IAA(1) = LQN(1)+LQN(2)+NUI+2*NX
        IAA(2) = LQN(3)+LQN(4)+NUI+2*NX
C
C       LOOP OVER EXPANSION TERMINALS FOR SECOND PAIR
        DO NY=1,NVALS
          IBB(1) = LQN(3)+LQN(4)-NUF+2*NY-1
          IBB(2) = LQN(1)+LQN(2)-NUF+2*NY-1
C
C         LOOP OVER BETA INTEGRAL COMBINATIONS
          DO IBETA=1,2
            IA =(IAA(IBETA)-1)/2
            IB = IBB(IBETA)   /2
C
C           BEGIN CONDITIONAL STATEMENT OVER IB VALUES
C           CASE 1: IB > 1
            IF(IB.GT.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BTA1(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
              RA  = X
              RB  = DFLOAT(1-IB)
              RC  = X+1.0D0
              RD  = 1.0D0
              RCD = RC*RD
              FCT = RA*RB/RCD
              DO M=1,MAXM
                TRM(M) = FCT*XJ(M,IBETA)
                BIN(M) = 1.0D0+TRM(M)
              ENDDO
              RA = RA+1.0D0
              RB = RB+1.0D0
              RC = RC+1.0D0
              RD = RD+1.0D0
              DO IT=2,IB-1
                RCD = RC*RD
                FCT = RA*RB/RCD
                DO M=1,MAXM
                  TRM(M) = FCT*TRM(M)*XJ(M,IBETA)
                  BIN(M) = BIN(M)+TRM(M)
                ENDDO
                RA = RA+1.0D0
                RB = RB+1.0D0
                RC = RC+1.0D0
                RD = RD+1.0D0
              ENDDO
              DO M=1,MAXM
                BETA(M) = BTA1(M)*BIN(M)
              ENDDO
C
C           CASE 2: IB = 1
            ELSEIF(IB.EQ.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BETA(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
C
C           CASE 3: IB = 0
            ELSEIF(IB.EQ.0) THEN
              DO M=1,MAXM
                XROOT(M) = DSQRT(XJ(M,IBETA))
                DEN      = (1.0D0-XROOT(M))
                RAT      = (1.0D0+XROOT(M))/DEN
                BTA1(M)  = DLOG(RAT)
                BIN(M)   = 1.0D0
                TRM(M)   = XJ(M,IBETA)
              ENDDO
              IF(IA.GT.1) THEN
                DO K=2,IA
                  KK = 2*K-1
                  RK = DFLOAT(KK)
                  X  = 1.0D0/RK
                  DO M=1,MAXM
                    BIN(M) = BIN(M)+X*TRM(M)
                    TRM(M) = TRM(M)*XJ(M,IBETA)
                  ENDDO
                ENDDO
                DO M=1,MAXM
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)*BIN(M)
                ENDDO
              ELSEIF(IA.EQ.1) THEN
                DO M=1,MAXM
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)
                ENDDO
              ELSE
                DO M=1,MAXM
                  BETA(M) = BTA1(M)
                ENDDO
              ENDIF
C
C           END CONDITIONAL STATEMENT OVER IB VALUES
            ENDIF
C
            IF(IBETA.EQ.1) THEN
              DO M=1,MAXM
                BDU(M, NUI+2*NX  ,-NUF+2*NY-1) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXM
                BDL(M, NUI+2*NX  ,-NUF+2*NY-1) = BETA(M)
              ENDDO
            ENDIF

          ENDDO

        ENDDO
      ENDDO
C
C**********************************************************************C
C     RADIAL INTEGRALS OVER SPINORS (SEPARATED BY TENSOR ORDER)        C
C**********************************************************************C
C
C     EMPTY COUNTER ARRAYS FOR DIRECT AND EXCHANGE INTEGRALS
      DO M=1,MAXM
        DO LTEN=1,NUNUM
          DO ILU=1,2
            RKLSLS(M,LTEN,ILU) = 0.0D0
            RKSLSL(M,LTEN,ILU) = 0.0D0
            RKSLLS(M,LTEN,ILU) = 0.0D0
            RKLSSL(M,LTEN,ILU) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     VALUE PREPARATION
      E0000 = 1.0D0
      E1000 = EI
      E0100 = EJ
C
C     INITIATE LOOP OVER K,L BASIS FUNCTIONS
      DO M=1,MAXM
C
C       MORE VALUE PREPARATION
        E0010 = EK(M)
        E0001 = EL(M)
        E1010 = EI*EK(M)
        E1001 = EI*EL(M)
        E0110 = EJ*EK(M)
        E0101 = EJ*EL(M)
C
C       LOOP OVER THE TENSOR ORDERS OF THE COULOMB INTERACTION
        DO LTEN=1,NUNUM
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(LTEN)
C
C         TEMPORARY STORAGE OF RAW RJ(LTEN,M)
          B21L = EIJ(-NU+1)*EKL(M, NU+2)*BDL(M, NU+2,-NU+1)
          B21U = EIJ( NU+2)*EKL(M,-NU+1)*BDU(M, NU+2,-NU+1)
          B23L = EIJ(-NU+3)*EKL(M, NU+2)*BDL(M, NU+2,-NU+3)
          B23U = EIJ( NU+4)*EKL(M,-NU+1)*BDU(M, NU+4,-NU+1)
          B41L = EIJ(-NU+1)*EKL(M, NU+4)*BDL(M, NU+4,-NU+1)
          B41U = EIJ( NU+2)*EKL(M,-NU+3)*BDU(M, NU+2,-NU+3)
          B43L = EIJ(-NU+3)*EKL(M, NU+4)*BDL(M, NU+4,-NU+3)
          B43U = EIJ( NU+4)*EKL(M,-NU+3)*BDU(M, NU+4,-NU+3)
C
C         EFFECTIVE INTERACTION STRENGTH RKLSLS(M,LTEN)
          RKLSLS(M,LTEN,1) 
     &         = V4*T0000*E0101*C7*B43L - V2*T0001*E0100*C5*B23L
     &         - V2*T0100*E0001*C5*B41L + V1*T0101*E0000*C3*B21L
C
          RKLSLS(M,LTEN,2) 
     &         = V4*T0000*E0101*C7*B43U - V2*T0001*E0100*C5*B23U
     &         - V2*T0100*E0001*C5*B41U + V1*T0101*E0000*C3*B21U
C
C         EFFECTIVE INTERACTION STRENGTH RKSLSL(M,LTEN)
          RKSLSL(M,LTEN,1)
     &         = V4*T0000*E1010*C7*B43L - V2*T0010*E1000*C5*B23L
     &         - V2*T1000*E0010*C5*B41L + V1*T1010*E0000*C3*B21L
C
          RKSLSL(M,LTEN,2)
     &         = V4*T0000*E1010*C7*B43U - V2*T0010*E1000*C5*B23U
     &         - V2*T1000*E0010*C5*B41U + V1*T1010*E0000*C3*B21U
C
C         EFFECTIVE INTERACTION STRENGTH RKSLLS(M,LTEN)
          RKSLLS(M,LTEN,1) 
     &         = V4*T0000*E1001*C7*B43L - V2*T1000*E0001*C5*B41L
     &         - V2*T0001*E1000*C5*B23L + V1*T1001*E0000*C3*B21L
C
          RKSLLS(M,LTEN,2) 
     &         = V4*T0000*E1001*C7*B43U - V2*T1000*E0001*C5*B41U
     &         - V2*T0001*E1000*C5*B23U + V1*T1001*E0000*C3*B21U
C
C         EFFECTIVE INTERACTION STRENGTH RKLSSL(M,LTEN)
          RKLSSL(M,LTEN,1) 
     &         = V4*T0000*E0110*C7*B43L - V2*T0100*E0010*C5*B23L
     &         - V2*T0010*E0100*C5*B41L + V1*T0110*E0000*C3*B21L
C
          RKLSSL(M,LTEN,2) 
     &         = V4*T0000*E0110*C7*B43U - V2*T0100*E0010*C5*B23U
     &         - V2*T0010*E0100*C5*B41U + V1*T0110*E0000*C3*B21U
C
C         SKIP POINT FOR NON-RELATIVISTIC HAMILTONIANS
999       CONTINUE
C
C       END LOOP OVER TENSOR ORDERS
        ENDDO
C
C     END LOOP OVER K,L BASIS FUNCTIONS
      ENDDO
C
C**********************************************************************C
C     APPLY NORMALISATION FACTORS                                      C
C**********************************************************************C
C
      DO M=1,MAXM
        RNLSLS = RNIJ(4)*RNKL(M,4)
        RNSLSL = RNIJ(2)*RNKL(M,2)
        RNSLLS = RNIJ(2)*RNKL(M,4)
        RNLSSL = RNIJ(4)*RNKL(M,2)
        DO LTEN=1,NUNUM
          DO ILU=1,2
            RKLSLS(M,LTEN,ILU) = RNLSLS*RKLSLS(M,LTEN,ILU)
            RKSLSL(M,LTEN,ILU) = RNSLSL*RKSLSL(M,LTEN,ILU)
            RKSLLS(M,LTEN,ILU) = RNSLLS*RKSLLS(M,LTEN,ILU)
            RKLSSL(M,LTEN,ILU) = RNLSSL*RKLSSL(M,LTEN,ILU)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ANGBRT1(DKAB,DKCD,EMAT,KQN,LQN,NUS,NUF,NUI,NUNUM,ISEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          AA    NN    NN  GGGGGG  BBBBBBB  RRRRRRR TTTTTTTT 11        C
C         AAAA   NNN   NN GG    GG BB    BB RR    RR   TT   111        C
C        AA  AA  NNNN  NN GG       BB    BB RR    RR   TT    11        C
C       AA    AA NN NN NN GG       BBBBBBB  RR    RR   TT    11        C
C       AAAAAAAA NN  NNNN GG   GGG BB    BB RRRRRRR    TT    11        C
C       AA    AA NN   NNN GG    GG BB    BB RR    RR   TT    11        C
C       AA    AA NN    NN  GGGGGG  BBBBBBB  RR    RR   TT   1111       C
C                                                                      C
C -------------------------------------------------------------------- C
C  ANGBRT1 PERFORMS THE ANGULAR ANALYSIS FOR THE EVALUATION OF TWO     C
C  ELECTRON INTEGRALS USING RACAH ALGEBRA TECHNIQUES. (OPEN-SHELL.)    C
C -------------------------------------------------------------------- C
C  DFNOTE: NOT WORKING YET...                                          C
C**********************************************************************C
      PARAMETER(MKP=9,MNU=MKP+1)
C
      DIMENSION KQN(4),LQN(4),JQN(4),NUS(MNU)
      DIMENSION DKAB(MNU,MKP+1,MKP+1),DKCD(MNU,MKP+1,MKP+1)
      DIMENSION EMAT(MNU,8),SCOEFF(8,2)
C
C     INITIALISE BREIT COEFFICIENTS
      DO INU=1,MNU
        DO IMU=1,8
          EMAT(INU,IMU) = 0.0D0
        ENDDO
      ENDDO
C
C**********************************************************************C
C     PARITY ANALYSIS: CHECK UNDERLYING LQN COMBINATIONS AND EXIT      C
C     IF THERE IS NO MULTIPOLE EXPANSION OF THE INTERACTION.           C
C**********************************************************************C
C
C     A AND B: LQN(1)+LQN(2) EVEN OR ODD
      IF(MOD(LQN(1)+LQN(2),2).EQ.0) THEN
        IPARAB = 1
      ELSE
        IPARAB = 0
      ENDIF
C
C     C AND D: LQN(3)+LQN(4) EVEN OR ODD
      IF(MOD(LQN(3)+LQN(4),2).EQ.0) THEN
        IPARCD = 1
      ELSE
        IPARCD = 0
      ENDIF
C
C     LQN SELECTION RULE: BOTH LQN PAIRS MUST BE OF SAME SYMMETRY
      IF(IPARAB.EQ.IPARCD) THEN
        ISEL = 1
      ELSE
        ISEL = 0
        RETURN
      ENDIF
C
C**********************************************************************C
C     MULTIPOLE EXPANSION UPPER/LOWER LIMITS SATISFY TRIANGLE RULE.    C
C**********************************************************************C
C
C     ASSIGN JQN VALUES
      DO N=1,4
        JQN(N) = 2*IABS(KQN(N))-1
      ENDDO
C
      NUI = MAX0(IABS(JQN(1)-JQN(2))/2,IABS(JQN(3)-JQN(4))/2)
      NUF = MIN0(    (JQN(1)+JQN(2))/2,    (JQN(3)+JQN(4))/2)
C
C     JQN SELECTION RULE: TRIANGLE RULE MUST PROVIDE VALID NU RANGE
      IF(NUI.GT.NUF) THEN
        ISEL = 0
        RETURN
      ELSE
        ISEL = 1
      ENDIF
C
C**********************************************************************C
C     LQN SELECTION RULE PARITY ANALYSIS: PAIRS BOTH EVEN OR BOTH ODD. C
C     ALSO CALCULATE BREIT FACTORS IN THIS SECTION.                    C
C**********************************************************************C
C
      NUNUM = 0
      DO NU=NUI,NUF
C
C       A AND B: LQN(1)+LQN(2)+NU EVEN OR ODD
        IF(MOD(LQN(1)+LQN(2)+NU,2).EQ.0) THEN
          IPARAB = 1
        ELSE
          IPARAB = 0
        ENDIF
C
C       C AND D: LQN(3)+LQN(4)+NU EVEN OR ODD
        IF(MOD(LQN(3)+LQN(4)+NU,2).EQ.0) THEN
          IPARCD = 1
        ELSE
          IPARCD = 0
        ENDIF
C
C       CASE 1: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH ODD (UNLESS NU=0)
        IF(IPARAB.EQ.0.AND.IPARCD.EQ.0.AND.NU.NE.0) THEN
C
C         BREIT COEFFICIENT IS SIMPLE ENOUGH TO CALCULATE HERE
          RNU  = DFLOAT(NU*(NU+1))
          COEF =-DFLOAT((KQN(1)+KQN(3))*(KQN(2)+KQN(4)))/RNU
C
C         INCREASE TOTAL NUMBER OF NU VALUES TO BE FACILITATED
          NUNUM = NUNUM+1
C
C         SAVE THIS PARTICULAR NU VALUE TO A LIST
          NUS(NUNUM) = NU
C
C         COPY ACROSS TO ARRAY
          DO IMU=1,8
            EMAT(NUNUM,IMU) = EMAT(NUNUM,IMU) + COEF
          ENDDO
C
        ENDIF
C
C       CASE 2: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH EVEN
        IF(IPARAB.EQ.1.AND.IPARCD.EQ.1) THEN
C
C         CALCULATE BREIT COEFFICIENTS WITH CALL TO BRCOEF1
          CALL BRCOEF1(SCOEFF,KQN,NU)
C
C         INCREASE TOTAL NUMBER OF NU VALUES TO BE FACILITATED
          NUNUM = NUNUM+1
C
C         SAVE THIS PARTICULAR NU VALUE TO A LIST
          NUS(NUNUM) = NU-1
C
C         COPY ACROSS TO ARRAY
          DO IMU=1,8
            EMAT(NUNUM,IMU) = EMAT(NUNUM,IMU) + SCOEFF(IMU,1)
          ENDDO
C
C         INCREASE TENSOR ORDER ONLY IF NU NONZERO
          IF(NU.GT.0) THEN
            NUNUM = NUNUM+1
          ENDIF
C
C         SAVE THIS PARTICULAR NU VALUE TO A LIST
          NUS(NUNUM) = NU+1
C          
C         COPY ACROSS TO ARRAY
          DO IMU=1,8
            EMAT(NUNUM,IMU) = EMAT(NUNUM,IMU) + SCOEFF(IMU,2)
          ENDDO
C
        ENDIF
C
      ENDDO
C
C     IF NO NU VALUES WERE ALLOWED DURING THIS RUN, EXIT PROCEDURE
      IF(NUNUM.EQ.0) THEN
        ISEL = 0
        RETURN
      ELSE
        ISEL = 1
      ENDIF
C
C     RE-SET THE LOWER AND UPPER LIMITS OF THE EXPANSION
      NUI = NUS(1)
      NUF = NUS(NUNUM)
C
C**********************************************************************C
C     ANGULAR FACTORS: EVALUATE AN ANGULAR FACTOR FOR EVERY |MQN|      C
C     COMBINATION IN THE MULTIPOLE EXPANSION OVER ALLOWED NU VALUES.   C
C**********************************************************************C
C
C     LOOP OVER MQN(A) AND MQN(B) VALUES
      DO MA=1,IABS(KQN(1))
        MJA = 2*MA-1
        DO MB=1,IABS(KQN(2))
          MJB = 2*MB-1
C
C         LOOP OVER ALL SURVIVING TENSOR ORDERS
          DO LTEN=1,NUNUM
C
C           READ THE ACTUAL ORDER NU
            NU = NUS(LTEN)
C
C           GENERATE AN ANGULAR FACTOR FOR ALL |MQN| SIGNS
            DKAB(LTEN,MJA  ,MJB  ) = DK(JQN(1),-MJA,JQN(2),-MJB,NU)
            DKAB(LTEN,MJA+1,MJB  ) = DK(JQN(1),+MJA,JQN(2),-MJB,NU)
            DKAB(LTEN,MJA  ,MJB+1) = DK(JQN(1),-MJA,JQN(2),+MJB,NU)
            DKAB(LTEN,MJA+1,MJB+1) = DK(JQN(1),+MJA,JQN(2),+MJB,NU)
C
          ENDDO
C
        ENDDO
      ENDDO
C
C     THE ARGUMENTS OF DK COEFFICIENTS ARE REVERSED IN THE CASE OF THE
C     CD PAIRS IN ORDER TO ACCOMMODATE THE RELATION:
C               DK(J,M,J',M',L) = ((-1)^Q)*DK(J'M',J,M,L).
C
C     LOOP OVER MQN(C) AND MQN(D) VALUES
      DO MC=1,IABS(KQN(3))
        MJC = 2*MC-1
        DO MD=1,IABS(KQN(4))
          MJD = 2*MD-1
C
C         LOOP OVER ALL SURVIVING TENSOR ORDERS
          DO LTEN=1,NUNUM
C
C           READ THE ACTUAL ORDER NU
            NU = NUS(LTEN)
C
C           GENERATE AN ANGULAR FACTOR FOR ALL |MQN| SIGNS
            DKCD(LTEN,MJC  ,MJD  ) = DK(JQN(4),-MJD,JQN(3),-MJC,NU)
            DKCD(LTEN,MJC+1,MJD  ) = DK(JQN(4),-MJD,JQN(3),+MJC,NU)
            DKCD(LTEN,MJC  ,MJD+1) = DK(JQN(4),+MJD,JQN(3),-MJC,NU)
            DKCD(LTEN,MJC+1,MJD+1) = DK(JQN(4),+MJD,JQN(3),+MJC,NU)
C
          ENDDO
C
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE BRCOEF1(SCOEFF,KQN,NU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      BBBBBBB  RRRRRRR   CCCCCC   OOOOOO  EEEEEEEE FFFFFFFF 11        C
C      BB    BB RR    RR CC    CC OO    OO EE       FF      111        C
C      BB    BB RR    RR CC       OO    OO EE       FF       11        C
C      BBBBBBB  RR    RR CC       OO    OO EEEEEE   FFFFFF   11        C
C      BB    BB RRRRRRR  CC       OO    OO EE       FF       11        C
C      BB    BB RR    RR CC    CC OO    OO EE       FF       11        C
C      BBBBBBB  RR    RR  CCCCCC   OOOOOO  EEEEEEEE FF      1111       C
C                                                                      C
C -------------------------------------------------------------------- C
C  BRCOEF1 EVALUATES THE INTERMEDIATE COEFFICIENTS OF THE BREIT        C
C  INTERACTION IN THE GENERAL CASE (TABLE 2 OF GRANT AND PYPER 1976).  C
C**********************************************************************C
      DIMENSION SCOEFF(8,2),KQN(4)
C
      RNU = DFLOAT(NU)
      RK1 = DFLOAT(KQN(3)-KQN(1))
      RK2 = DFLOAT(KQN(4)-KQN(2))
C
      IF(NU-1.GE.0) THEN
        B1 = DFLOAT(NU+1)/DFLOAT(2   *(2*NU+1))
        C1 =-DFLOAT(NU-2)/DFLOAT(2*NU*(2*NU+1))
        SCOEFF(1,1) = (RNU+RK1)*(B1+C1*RK2)
        SCOEFF(2,1) = (RNU+RK2)*(B1+C1*RK1)
        SCOEFF(3,1) = (RNU-RK1)*(B1-C1*RK2)
        SCOEFF(4,1) = (RNU-RK2)*(B1-C1*RK1)
        SCOEFF(5,1) =-(RNU+RK1)*(B1-C1*RK2)
        SCOEFF(6,1) =-(RNU-RK2)*(B1+C1*RK1)
        SCOEFF(7,1) =-(RNU-RK1)*(B1+C1*RK2)
        SCOEFF(8,1) =-(RNU+RK2)*(B1-C1*RK1)
      ELSE
        SCOEFF(1,1) = 0.0D0
        SCOEFF(2,1) = 0.0D0
        SCOEFF(3,1) = 0.0D0
        SCOEFF(4,1) = 0.0D0
        SCOEFF(5,1) = 0.0D0
        SCOEFF(6,1) = 0.0D0
        SCOEFF(7,1) = 0.0D0
        SCOEFF(8,1) = 0.0D0
      ENDIF
C
      IF(NU+1.GE.1) THEN
        B2 = DFLOAT(NU  )/DFLOAT(2       *(2*NU+3))
        C2 = DFLOAT(NU+3)/DFLOAT(2*(NU+1)*(2*NU+3))
        SCOEFF(1,2) = ( RK2-RNU-1.0D0)*(B2+C2*RK1)
        SCOEFF(2,2) = ( RK1-RNU-1.0D0)*(B2+C2*RK2)
        SCOEFF(3,2) = (-RK2-RNU-1.0D0)*(B2-C2*RK1)
        SCOEFF(4,2) = (-RK1-RNU-1.0D0)*(B2-C2*RK2)
        SCOEFF(5,2) =-(-RK2-RNU-1.0D0)*(B2+C2*RK1)
        SCOEFF(6,2) =-( RK1-RNU-1.0D0)*(B2-C2*RK2)
        SCOEFF(7,2) =-( RK2-RNU-1.0D0)*(B2-C2*RK1)
        SCOEFF(8,2) =-(-RK1-RNU-1.0D0)*(B2+C2*RK2)
      ELSE
        SCOEFF(1,2) = 0.0D0
        SCOEFF(2,2) = 0.0D0
        SCOEFF(3,2) = 0.0D0
        SCOEFF(4,2) = 0.0D0
        SCOEFF(5,2) = 0.0D0
        SCOEFF(6,2) = 0.0D0
        SCOEFF(7,2) = 0.0D0
        SCOEFF(8,2) = 0.0D0
      ENDIF
C
      RETURN
      END

