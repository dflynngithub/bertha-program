      SUBROUTINE BREIT1old(IZ)
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
C  BREIT1 CONSTRUCTS A ONE-CENTRE CONTRIBUTION TO THE MOLECULAR        C
C  MEAN-FIELD BREIT MATRIX WITH RACAH ALGEBRA AND BETA INTEGRALS.      C
C  THIS ROUTINE IS SIMILAR TO THE IN-LINE CONSTRUCTION OF MEAN-FIELD   C
C  ATOMIC BREIT MATRIX IN HFSCF0, BUT WITH MQN STRUCTURE AS WELL.      C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      DIMENSION ANGFAC(8)
      DIMENSION BCOF(MNU,8,MKP+1,MKP+1,MKP+1,MKP+1)
      DIMENSION RJLSLS(MB2,MNU,2),RJSLSL(MB2,MNU,2),
     &          RJLSSL(MB2,MNU,2),RJSLLS(MB2,MNU,2)
      DIMENSION XLSLS(MB2),XSLSL(MB2),XLSSL(MB2),XSLLS(MB2)
C
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VANM(MDM,MDM),
     &           VSLF(MDM,MDM),VUEH(MDM,MDM),VWKR(MDM,MDM),
     &           VKSB(MDM,MDM),QDIR(MDM,MDM),QXCH(MDM,MDM),
     &           WDIR(MDM,MDM),WXCH(MDM,MDM),CPLE(MDM,MDM)
C
      COMMON/B0IJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ,IBAS,JBAS
      COMMON/B1QN/EXL(MBS,4),MQN(4),KQN(4),LQN(4),NBAS(4),MAXM
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/DENS/DENC,DENO,DENT
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VANM,VSLF,
     &            VUEH,VWKR,VKSB,QDIR,QXCH,WDIR,WXCH,CPLE
      COMMON/RCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/SHLL/ACFF,BCFF,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
      COMMON/T2EL/F2ES(5,9),T2ES(5,9),N2EB(5,9),N2EI(5,9),N2ES(5,9)
      COMMON/TSCF/TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRW,TCC1,TCC2,TCMC,
     &            TB1B,TB1R,TB1F,TB1M,TBEC,TBRM,TBRW,TBC1,TBC2,TBMC,
     &            TSMX,TUMX,THMX,TAMX,TC1T,TC2T,TCVT,TB1T,TB2T,TACC,
     &            TEIG,TSCR,TTOT,TC1S,TC2S,TB1S,TB2S
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
C
C     ANGULAR FACTOR SENSITIVITY PARAMETER
      DATA SENS/1.0D-10/
C
      CALL CPU_TIME(TBCH1)
C
C     INTEGRAL SKIPPING ON MOLECULAR GROUP SYMMETRY CLASS BASIS
      IF(SHAPE.EQ.'ATOMIC') THEN
        ISYM = 2
      ELSEIF(SHAPE.EQ.'DIATOM'.OR.SHAPE.EQ.'LINEAR') THEN
        ISYM = 1
      ELSE
        ISYM = 0
      ENDIF
C
C**********************************************************************C
C     LOOP OVER KQN SYMMETRY TYPE FOR A,B,C,D BLOCKS (USE INDEX 1000)  C
C**********************************************************************C
C
C     LOOP OVER KQN(A) VALUES
      DO 1000 KA=1,NKAP(IZ)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,IZ)
        LQN(1) = LVAL(KQN(1))
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),IZ)
        DO IBAS=1,NBAS(1)
          EXL(IBAS,1) = BEXL(IBAS,LQN(1),IZ)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 1000 KB=1,NKAP(IZ)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KAPA(KB,IZ)
        LQN(2) = LVAL(KQN(2))
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),IZ)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),IZ)
        ENDDO
C
C     LOOP OVER KQN(C) VALUES
      DO 1000 KC=1,NKAP(IZ)
C
C       QUANTUM NUMBERS FOR BLOCK C
        KQN(3) = KAPA(KC,IZ)
        LQN(3) = LVAL(KQN(3))
C
C       BASIS EXPONENTS FOR BLOCK C
        NBAS(3) = NFNC(LQN(3),IZ)
        DO KBAS=1,NBAS(3)
          EXL(KBAS,3) = BEXL(KBAS,LQN(3),IZ)
        ENDDO
C
C     LOOP OVER KQN(D) VALUES
      DO 1000 KD=1,NKAP(IZ)
C
C       QUANTUM NUMBERS FOR BLOCK D
        KQN(4) = KAPA(KD,IZ)
        LQN(4) = LVAL(KQN(4))
C
C       BASIS EXPONENTS FOR BLOCK D
        NBAS(4) = NFNC(LQN(4),IZ)
        DO LBAS=1,NBAS(4)
          EXL(LBAS,4) = BEXL(LBAS,LQN(4),IZ)
        ENDDO
C
C       NUMBER OF BASIS FUNCTIONS IN (CD) BLOCK
        MAXM = NBAS(3)*NBAS(4)
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
      T1100 = TI*TJ
      T1010 = TI*TK
      T1001 = TI*TL
      T0110 = TJ*TK
      T0101 = TJ*TL
      T0011 = TK*TL
      T1110 = TI*TJ*TK
      T1101 = TI*TJ*TL
      T1011 = TI*TK*TL
      T0111 = TJ*TK*TL
      T1111 = TI*TJ*TK*TL
C
C     ANGULAR COEFFICIENTS
      CALL CPU_TIME(T1)
      CALL ANGBRT1old(BCOF,KQN,LQN,ISEL)
      CALL CPU_TIME(T2)
      TB1B = TB1B+T2-T1
C
C     EXIT THIS COMBINATION IF IT VIOLATES A SELECTION RULE
      IF(ISEL.EQ.0) GOTO 1001
C
C     BASIS SET INTERMEDIATES FOR THE KL-PAIRS
      CALL CPU_TIME(T1)
      CALL KLSET1(IZ)
      CALL CPU_TIME(T2)
      TB1B = TB1B+T2-T1
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS IN BLOCKS A AND B (INDEX 2000)         C
C**********************************************************************C
C
C     UPDATE COUNTER FOR NUMBER OF CLASSES
      N2EB(1,5) = N2EB(1,5)+1
C
      DO 2000 IBAS=1,NBAS(1)
      DO 2000 JBAS=1,NBAS(2)
C
C       UPDATE COUNTER FOR NUMBER OF INTEGRALS
        N2EI(1,5) = N2EI(1,5) + NBAS(3)*NBAS(4)
C
C       BASIS SET INTERMEDIATES FOR THE IJ-PAIRS
        CALL CPU_TIME(T1)
        CALL IJSET1
        CALL CPU_TIME(T2)
        TB1B = TB1B+T2-T1
C
C       BATCH OF RADIAL INTEGRALS (EFFECTIVE INTERACTION STRENGTHS)
        CALL CPU_TIME(T1)
        CALL RKBRT1old(RJLSLS,RJSLSL,RJLSSL,RJSLLS)
        CALL CPU_TIME(T2)
        TB1R = TB1R+T2-T1
C
C       CAN NOW CONTRACT THESE RADIAL INTEGRALS OVER ANGULAR COMPONENTS
C       OF G-SPINOR BASIS FUNCTIONS USING A TENSOR EXPANSION IN {L,Q}
C
C**********************************************************************C
C     LOOP OVER |KQN| MAGNITUDES FOR A,B,C,D BLOCKS (USE INDEX 3000)   C
C**********************************************************************C
C
      DO 3000 MA=1,IABS(KQN(1))
        MQN(1) = 2*MA-1
C
      DO 3000 MB=1,IABS(KQN(2))
        MQN(2) = 2*MB-1
C
      DO 3000 MC=1,IABS(KQN(3))
        MQN(3) = 2*MC-1
C
      DO 3000 MD=1,IABS(KQN(4))
        MQN(4) = 2*MD-1
C
C**********************************************************************C
C     LOOP OVER THE SIGNS OF |MQN| FOR A,B,C,D BLOCKS (USE INDEX 4000) C
C**********************************************************************C
C
      DO 4000 ISGN1=1,2
        MMJA = MQN(1)*((-1)**ISGN1)
        IMJA = MQN(1)+ISGN1-1
C
      DO 4000 ISGN2=1,2
        MMJB = MQN(2)*((-1)**ISGN2)
        IMJB = MQN(2)+ISGN2-1
C
      DO 4000 ISGN3=1,2
        MMJC = MQN(3)*((-1)**ISGN3)
        IMJC = MQN(3)+ISGN3-1
C
      DO 4000 ISGN4=1,2
        MMJD = MQN(4)*((-1)**ISGN4)
        IMJD = MQN(4)+ISGN4-1
C
C     STARTING FOCK ADDRESS FOR EACH BASIS LIST
      NAL = LRGE(IZ,KA,IMJA)
      NBL = LRGE(IZ,KB,IMJB)
      NCL = LRGE(IZ,KC,IMJC)
      NDL = LRGE(IZ,KD,IMJD)
C
      NAS = LRGE(IZ,KA,IMJA) + NSKP
      NBS = LRGE(IZ,KB,IMJB) + NSKP
      NCS = LRGE(IZ,KC,IMJC) + NSKP
      NDS = LRGE(IZ,KD,IMJD) + NSKP
C
C     APPLY ANGULAR MQN SELECTION RULE
      IF(MMJA-MMJB.NE.MMJD-MMJC) GOTO 4001
      IF(ISYM.EQ.1.OR.ISYM.EQ.2) THEN
        IF(ISGN1.EQ.ISGN2.AND.ISGN3.EQ.ISGN4) GOTO 4002
        IF(ISGN1.EQ.ISGN4.AND.ISGN3.EQ.ISGN2) GOTO 4002     
        GOTO 4001
      ENDIF
4002  CONTINUE
C
C     RESET CONTRACTED RADIAL ARRAYS
      CALL CPU_TIME(T1)
      DO M=1,NBAS(3)*NBAS(4)
        XLSLS(M) = 0.0D0
        XSLSL(M) = 0.0D0
        XLSSL(M) = 0.0D0
        XSLLS(M) = 0.0D0
      ENDDO
C
C     ASSEMBLE INTEGRAL BY SUMMING OVER EXPANSION POWER L
      DO LTEN=1,NUNUM
C
C       ANGULAR COEFFICIENTS
        DO IMU=1,8
          ANGFAC(IMU) = BCOF(LTEN,IMU,IMJA,IMJB,IMJC,IMJD)
        ENDDO
C
C       SCREENING OF INTEGRAL BASED ON ANGULAR COEFFICIENT
        ANGSUM = 0.0D0
        DO IMU=1,8
          ANGSUM = ANGSUM + DABS(ANGFAC(IMU))
        ENDDO
        IF(DABS(ANGSUM).LE.SENS) GOTO 4003
C
        DO M=1,NBAS(3)*NBAS(4)
C
C         FACTORS LEADING TO BXCH(LS)
          XLSLS(M) = XLSLS(M) + ANGFAC(1)*RJLSLS(M,LTEN,1)
     &                        + ANGFAC(2)*RJLSLS(M,LTEN,2)
C
C         FACTORS LEADING TO BXCH(SL)
          XSLSL(M) = XSLSL(M) + ANGFAC(3)*RJSLSL(M,LTEN,1)
     &                        + ANGFAC(4)*RJSLSL(M,LTEN,2)
C
C         FACTORS LEADING TO BXCH(LL)
          XLSSL(M) = XLSSL(M) + ANGFAC(5)*RJLSSL(M,LTEN,1)
     &                        + ANGFAC(6)*RJLSSL(M,LTEN,2)
C
C         FACTORS LEADING TO BXCH(SS)
          XSLLS(M) = XSLLS(M) + ANGFAC(7)*RJSLLS(M,LTEN,1)
     &                        + ANGFAC(8)*RJSLLS(M,LTEN,2)
C
        ENDDO
C
C       SKIP POINT FOR ANGULAR SCREENING
4003    CONTINUE
C
      ENDDO
C
C     FULL-INTEGRAL CONSTRUCTION COMPLETE
      CALL CPU_TIME(T2)
      TB1F = TB1F+T2-T1
C
C**********************************************************************C
C     ADD THIS BATCH OF R-INTEGRALS TO CLOSED-SHELL BREIT MATRIX       C
C -------------------------------------------------------------------- C
C          (NO BDIR CONTRIBUTIONS FOR CLOSED-SHELL SYSTEMS.)           C
C**********************************************************************C
C
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
C
C         EXCHANGE CONTRIBUTIONS
          BXCH(NAL+IBAS,NDS+LBAS) = BXCH(NAL+IBAS,NDS+LBAS)
     &                  +      XLSLS(M)*DCONJG(DENT(NBS+JBAS,NCL+KBAS))
C
          BXCH(NAS+IBAS,NDL+LBAS) = BXCH(NAS+IBAS,NDL+LBAS)
     &                  +      XSLSL(M)*DCONJG(DENT(NBL+JBAS,NCS+KBAS))
C
          BXCH(NAL+IBAS,NDL+LBAS) = BXCH(NAL+IBAS,NDL+LBAS)
     &                  +      XLSSL(M)*DCONJG(DENT(NBS+JBAS,NCS+KBAS))
C
          BXCH(NAS+IBAS,NDS+LBAS) = BXCH(NAS+IBAS,NDS+LBAS)
     &                  +      XSLLS(M)*DCONJG(DENT(NBL+JBAS,NCL+KBAS))
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     ADD THIS BATCH OF R-INTEGRALS TO OPEN-SHELL BREIT MATRIX.        C
C     THIS ALSO REQUIRES THE CLOSED-SHELL DIRECT MATRIX ELEMENTS.      C
C**********************************************************************C
C
C
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
C
C         DIRECT CONTRIBUTIONS TO CLOSED-SHELL MATRIX
          BDIR(NAL+IBAS,NBS+JBAS) = BDIR(NAL+IBAS,NBS+JBAS)
     &                  +      XLSLS(M)*DCONJG(DENT(NCL+KBAS,NDS+LBAS))
     &                  +      XLSSL(M)*DCONJG(DENT(NCS+KBAS,NDL+LBAS))
C
          BDIR(NAS+IBAS,NBL+JBAS) = BDIR(NAS+IBAS,NBL+JBAS)
     &                  +      XSLSL(M)*DCONJG(DENT(NCS+KBAS,NDL+LBAS))
     &                  +      XSLLS(M)*DCONJG(DENT(NCL+KBAS,NDS+LBAS))
C
        ENDDO
      ENDDO

      IF(NOPN.EQ.0) GOTO 5100
C
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
C
C         DIRECT CONTRIBUTIONS TO OPEN-SHELL MATRIX
          WDIR(NAL+IBAS,NBS+JBAS) = WDIR(NAL+IBAS,NBS+JBAS)
     &                  + ACFF*XLSLS(M)*DCONJG(DENO(NCL+KBAS,NDS+LBAS))
     &                  + ACFF*XLSSL(M)*DCONJG(DENO(NCS+KBAS,NDL+LBAS))
C
          WDIR(NAS+IBAS,NBL+JBAS) = WDIR(NAS+IBAS,NBL+JBAS)
     &                  + ACFF*XSLSL(M)*DCONJG(DENO(NCS+KBAS,NDL+LBAS))
     &                  + ACFF*XSLLS(M)*DCONJG(DENO(NCL+KBAS,NDS+LBAS))
C
C         EXCHANGE CONTRIBUTIONS TO OPEN-SHELL MATRIX
          WXCH(NAL+IBAS,NDS+LBAS) = WXCH(NAL+IBAS,NDS+LBAS)
     &                  + BCFF*XLSLS(M)*DCONJG(DENO(NBS+JBAS,NCL+KBAS))
C
          WXCH(NAS+IBAS,NDL+LBAS) = WXCH(NAS+IBAS,NDL+LBAS)
     &                  + BCFF*XSLSL(M)*DCONJG(DENO(NBL+JBAS,NCS+KBAS))
C
          WXCH(NAL+IBAS,NDL+LBAS) = WXCH(NAL+IBAS,NDL+LBAS)
     &                  + BCFF*XLSSL(M)*DCONJG(DENO(NBS+JBAS,NCS+KBAS))
C
          WXCH(NAS+IBAS,NDS+LBAS) = WXCH(NAS+IBAS,NDS+LBAS)
     &                  + BCFF*XSLLS(M)*DCONJG(DENO(NBL+JBAS,NCL+KBAS))
C
        ENDDO
      ENDDO
C
5100  CONTINUE
C
C     MATRIX MULTIPLICATION STEP COMPLETE
      CALL CPU_TIME(T3)
      TB1M = TB1M+T3-T2
C
C**********************************************************************C
C     COMPLETE CONSTRUCTION OF ALL ONE-CENTRE CONTRIBUTIONS.           C
C**********************************************************************C
C
C     EARLY EXIT FOR MQN SELECTION RULES
4001  CONTINUE
C     END LOOP OVER ALL |MQN| SIGNS
4000  CONTINUE
C     END LOOP OVER ALL |MQN| MAGNITUDES
3000  CONTINUE
C     END LOOP OVER IBAS AND JBAS
2000  CONTINUE
C     EARLY EXIT FOR KQN SELECTION RULES
1001  CONTINUE
C     END LOOP OVER ALL KQNS
1000  CONTINUE
C
C     RECORD CPU TIME AT END OF BATCH AND ADD TO APPROPRIATE COUNTER
      CALL CPU_TIME(TBCH2)
      IF(INTSYM) THEN
        T2ES(1,5) = T2ES(1,5)+TBCH2-TBCH1
      ELSE
        T2ES(1,5) = T2ES(1,5) + 0.25D0*(TBCH2-TBCH1)
        T2ES(1,6) = T2ES(1,6) + 0.25D0*(TBCH2-TBCH1)
        T2ES(1,7) = T2ES(1,7) + 0.25D0*(TBCH2-TBCH1)
        T2ES(1,8) = T2ES(1,8) + 0.25D0*(TBCH2-TBCH1)
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE RKBRT1old(RJLSLS,RJSLSL,RJLSSL,RJSLLS)
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
C  RKBRT1 EVALUATES A BATCH OF GENERAL (OPEN-SHELL) ONE-CENTRE RADIAL  C
C  INTEGRALS OVER THE BREIT INTERACTION, SEPARATING RESULTS BY THE     C
C  ALLOWED TENSOR ORDERS RK(ABCD).                                     C
C -------------------------------------------------------------------- C
C  RJLSLS(MB2,MNU,IUL) - JLSLS INTEGRAL LIST OF TENSOR TYPE MNU        C
C  RJSLSL(MB2,MNU,IUL) - JSLSL INTEGRAL LIST OF TENSOR TYPE MNU        C
C  RJLSSL(MB2,MNU,IUL) - JLSSL INTEGRAL LIST OF TENSOR TYPE MNU        C
C  RJSLLS(MB2,MNU,IUL) - JSLLS INTEGRAL LIST OF TENSOR TYPE MNU        C
C -------------------------------------------------------------------- C
C  INDEX 1 FOR UPPER INTEGRALS, INDEX 2 FOR LOWER INTEGRALS.           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION XJ(MB2,2),IAA(2),IBB(2)
      DIMENSION BETA(MB2),BTA1(MB2),XROOT(MB2),BIN(MB2),TRM(MB2)
      DIMENSION BDU(MB2,-MAB:MAB,-MAB:MAB),BDL(MB2,-MAB:MAB,-MAB:MAB)
C
      DIMENSION RJLSLS(MB2,MNU,2),RJSLSL(MB2,MNU,2),
     &          RJLSSL(MB2,MNU,2),RJSLLS(MB2,MNU,2)
C
      COMMON/B0IJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ,IBAS,JBAS
      COMMON/B0KL/EKL(MB2,-MAB:MAB),RNKL(MB2,4),EK(MB2),EL(MB2)
      COMMON/B1QN/EXL(MBS,4),MQN(4),KQN(4),LQN(4),NBAS(4),MAXM
      COMMON/RCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
C
C**********************************************************************C
C     PREPARATION OF EXPONENT POWERS AND BETA FUNCTION ARGUMENTS.      C
C**********************************************************************C
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
C     NUMBER OF LEVELS NEEDED TO ACCOUNT FOR ALL BETA INTEGRALS
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
          DO IUL=1,2
            RJLSLS(M,LTEN,IUL) = 0.0D0
            RJSLSL(M,LTEN,IUL) = 0.0D0
            RJLSSL(M,LTEN,IUL) = 0.0D0
            RJSLLS(M,LTEN,IUL) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     VALUE PREPARATION
      E0000 = 1.0D0
      E1000 = EI
      E0100 = EJ
      E1100 = EI*EJ
C
C     LOOP OVER K,L BASIS FUNCTIONS
      DO M=1,MAXM
C
C       MORE VALUE PREPARATION
        E0010 = EK(M)
        E0001 = EL(M)
        E1010 = EI*EK(M)
        E1001 = EI*EL(M)
        E0110 = EJ*EK(M)
        E0101 = EJ*EL(M)
        E0011 = EK(M)*EL(M)
C
C       LOOP OVER THE TENSOR ORDERS OF THE BREIT INTERACTION
        DO LTEN=1,NUNUM
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(LTEN)
C
C         TEMPORARY STORAGE OF RAW RJ(LTEN,M)
          B21L = EIJ(-NU+1)*EKL(M, NU+2)*BDL(M, NU+2,-NU+1)
          B21U = EIJ( NU+2)*EKL(M,-NU+1)*BDU(M, NU+2,-NU+1)
          B23L = EIJ(-NU+3)*EKL(M, NU+2)*BDL(M, NU+2,-NU+3)
          B41U = EIJ( NU+4)*EKL(M,-NU+1)*BDU(M, NU+4,-NU+1)
          B41L = EIJ(-NU+1)*EKL(M, NU+4)*BDL(M, NU+4,-NU+1)
          B23U = EIJ( NU+2)*EKL(M,-NU+3)*BDU(M, NU+2,-NU+3)
          B43L = EIJ(-NU+3)*EKL(M, NU+4)*BDL(M, NU+4,-NU+3)
          B43U = EIJ( NU+4)*EKL(M,-NU+3)*BDU(M, NU+4,-NU+3)
C
C         EFFECTIVE INTERACTION STRENGTH RJLSLS(M,LTEN)
          RJLSLS(M,LTEN,1) 
     &         = V4*T0000*E0101*C7*B43U - V2*T0001*E0100*C5*B41U
     &         - V2*T0100*E0001*C5*B23U + V1*T0101*E0000*C3*B21U
C
          RJLSLS(M,LTEN,2) 
     &         = V4*T0000*E0101*C7*B43L - V2*T0001*E0100*C5*B23L
     &         - V2*T0100*E0001*C5*B41L + V1*T0101*E0000*C3*B21L
C
C         EFFECTIVE INTERACTION STRENGTH RJSLSL(M,LTEN)
          RJSLSL(M,LTEN,1)
     &         = V4*T0000*E1010*C7*B43U - V2*T0010*E1000*C5*B41U
     &         - V2*T1000*E0010*C5*B23U + V1*T1010*E0000*C3*B21U
C
          RJSLSL(M,LTEN,2)
     &         = V4*T0000*E1010*C7*B43L - V2*T0010*E1000*C5*B23L
     &         - V2*T1000*E0010*C5*B41L + V1*T1010*E0000*C3*B21L
C
C         EFFECTIVE INTERACTION STRENGTH RJLSSL(M,LTEN)
          RJLSSL(M,LTEN,1) 
     &         = V4*T0000*E0110*C7*B43U - V2*T0010*E0100*C5*B41U
     &         - V2*T0100*E0010*C5*B23U + V1*T0110*E0000*C3*B21U
C
          RJLSSL(M,LTEN,2)
     &         = V4*T0000*E0110*C7*B43L - V2*T0010*E0100*C5*B23L
     &         - V2*T0100*E0010*C5*B41L + V1*T0110*E0000*C3*B21L
C
C         EFFECTIVE INTERACTION STRENGTH RJSLLS(M,LTEN)
          RJSLLS(M,LTEN,1) 
     &         = V4*T0000*E1001*C7*B43U - V2*T0001*E1000*C5*B41U
     &         - V2*T1000*E0001*C5*B23U + V1*T1001*E0000*C3*B21U
C
          RJSLLS(M,LTEN,2) 
     &         = V4*T0000*E1001*C7*B43L - V2*T0001*E1000*C5*B23L
     &         - V2*T1000*E0001*C5*B41L + V1*T1001*E0000*C3*B21L
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
        RNLSLS = RNIJ(2)*RNKL(M,2)
        RNSLSL = RNIJ(3)*RNKL(M,3)
        RNSLLS = RNIJ(3)*RNKL(M,2)
        RNLSSL = RNIJ(2)*RNKL(M,3)
        DO LTEN=1,NUNUM
          DO IUL=1,2
            RJLSLS(M,LTEN,IUL) = RNLSLS*RJLSLS(M,LTEN,IUL)
            RJSLSL(M,LTEN,IUL) = RNSLSL*RJSLSL(M,LTEN,IUL)
            RJLSSL(M,LTEN,IUL) = RNLSSL*RJLSSL(M,LTEN,IUL)
            RJSLLS(M,LTEN,IUL) = RNSLLS*RJSLLS(M,LTEN,IUL)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ANGBRT1old(BCOF,KQN,LQN,ISEL)
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
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION KQN(4),LQN(4),JQN(4)
      DIMENSION SCOEFF(8,2)
      DIMENSION DKAB(0:MNU,MKP+1,MKP+1),DKCD(0:MNU,MKP+1,MKP+1)
      DIMENSION AJJN(MKP+1,MKP+1,MKP+1,MKP+1)
      DIMENSION BCOF(MNU,8,MKP+1,MKP+1,MKP+1,MKP+1)
C
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
C
C     INITIALISE BREIT COEFFICIENTS
      DO LTEN=1,MNU
        DO IMU=1,8
          DO IMJA=1,2*IABS(KQN(1))
            DO IMJB=1,2*IABS(KQN(2))
              DO IMJC=1,2*IABS(KQN(3))
                DO IMJD=1,2*IABS(KQN(4))
                  BCOF(LTEN,IMU,IMJA,IMJB,IMJC,IMJD) = 0.0D0
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        NUS(LTEN) = 0
      ENDDO
      NUNUM = 0
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
C     TENSOR LIMITS BY TRIANGLE RULE
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
C     LQN PARITY ANALYSIS: CHECK PARITY OF LQN COMBINATIONS.           C
C     EXIT IF THERE IS NO MULTIPOLE EXPANSION FOR THIS CASE.           C
C**********************************************************************C
C
C     A AND B: PARITY OF 'LQN(1)+LQN(2)' (0 IF EVEN, 1 IF ODD)
      IPARAB = MOD(LQN(1)+LQN(2),2)
C
C     C AND D: PARITY OF 'LQN(3)+LQN(4)' (0 IF EVEN, 1 IF ODD)
      IPARCD = MOD(LQN(3)+LQN(4),2)
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
C     WIGNER-ECKART: EVALUATE AN ANGULAR FACTOR FOR EVERY |MQN|        C
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
          DO NU=NUI,NUF
C
C           GENERATE AN ANGULAR FACTOR FOR ALL |MQN| SIGNS
            DKAB(NU,MJA  ,MJB  ) = DK(JQN(1),-MJA,JQN(2),-MJB,NU)
            DKAB(NU,MJA+1,MJB  ) = DK(JQN(1), MJA,JQN(2),-MJB,NU)
            DKAB(NU,MJA  ,MJB+1) = DK(JQN(1),-MJA,JQN(2), MJB,NU)
            DKAB(NU,MJA+1,MJB+1) = DK(JQN(1), MJA,JQN(2), MJB,NU)
C
          ENDDO
C
        ENDDO
      ENDDO
C
C     THE ARGUMENTS OF DK COEFFICIENTS ARE REVERSED IN THE CASE OF THE
C     CD PAIRS IN ORDER TO ACCOMMODATE THE RELATION:
C               DK(J,M,J',M',L) = ((-1)^Q)*DK(J',M',J,M,L).
C
C     LOOP OVER MQN(C) AND MQN(D) VALUES
      DO MC=1,IABS(KQN(3))
        MJC = 2*MC-1
        DO MD=1,IABS(KQN(4))
          MJD = 2*MD-1
C
C         LOOP OVER ALL SURVIVING TENSOR ORDERS
          DO NU=NUI,NUF
C
C           GENERATE AN ANGULAR FACTOR FOR ALL |MQN| SIGNS
            DKCD(NU,MJC  ,MJD  ) = DK(JQN(4),-MJD,JQN(3),-MJC,NU)
            DKCD(NU,MJC+1,MJD  ) = DK(JQN(4),-MJD,JQN(3), MJC,NU)
            DKCD(NU,MJC  ,MJD+1) = DK(JQN(4), MJD,JQN(3),-MJC,NU)
            DKCD(NU,MJC+1,MJD+1) = DK(JQN(4), MJD,JQN(3), MJC,NU)
C
          ENDDO
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CALCULATE COMPONENT-TYPE OVERLAP COEFFICIENTS FOR TENSOR ORDERS  C
C     NU (AND NU±1 WHERE APPLICABLE), COUPLING IN ANGULAR FACTORS.     C
C**********************************************************************C
C
      ISEL = 0
      LTEN = 1
      DO NU=NUI,NUF
C
C       A AND B: PARITY OF 'LA+LB+NU' (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQN(1)+LQN(2)+NU,2)
C
C       C AND D: PARITY OF 'LC+LD+NU' (0 IF EVEN, 1 IF ODD)
        IPARCD = MOD(LQN(3)+LQN(4)+NU,2)
C
C       JQN AND NU ANGULAR FACTORS FOR THIS COMBINATION
        DO IMJA=1,2*IABS(KQN(1))
          DO IMJB=1,2*IABS(KQN(2))
            DO IMJC=1,2*IABS(KQN(3))
              DO IMJD=1,2*IABS(KQN(4))
                AJJN(IMJA,IMJB,IMJC,IMJD)
     &           = DKAB(NU,IMJA,IMJB)*DKCD(NU,IMJC,IMJD)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C
C       CASE 1: LA+LB+NU AND LC+LD+NU ARE BOTH ODD (AND NU=/=0)
        IF(IPARAB.EQ.1.AND.IPARCD.EQ.1.AND.NU.NE.0) THEN
C
C         THERE IS A NON-TRIVIAL TENSOR EXPANSION
          ISEL = 1
C
C         INTERMEDIATE (KQN,NU) COEFFICIENTS
          RNU    = DFLOAT(NU*(NU+1))
          RCOEFF =-DFLOAT((KQN(1)+KQN(2))*(KQN(3)+KQN(4)))/RNU
C
C         THIS CONTRIBUTION BELONGS IN "NU" BIN
          NUS(LTEN) = NU
C
C         CONTRIBUTIONS TO TENSOR ORDER "NU"
          DO IMU=1,8
C
            DO IMJA=1,2*IABS(KQN(1))
              DO IMJB=1,2*IABS(KQN(2))
                DO IMJC=1,2*IABS(KQN(3))
                  DO IMJD=1,2*IABS(KQN(4))
                    BCOF(LTEN,IMU,IMJA,IMJB,IMJC,IMJD)
     &               = BCOF(LTEN,IMU,IMJA,IMJB,IMJC,IMJD)
     &               +          AJJN(IMJA,IMJB,IMJC,IMJD)*RCOEFF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
C
          ENDDO
C
        ENDIF
C
C       CASE 2: LA+LB+NU AND LC+LD+NU ARE BOTH EVEN
        IF(IPARAB.EQ.0.AND.IPARCD.EQ.0) THEN
C
C         THERE IS A NON-TRIVIAL TENSOR EXPANSION
          ISEL = 1
C
C         CALCULATE BREIT COEFFICIENTS WITH CALL TO BRCOEF1
          CALL BRCOEF1(SCOEFF,KQN,NU)
C
C         FIRST CONTRIBUTION BELONGS IN "NU-1" BIN
          NUS(LTEN) = NU-1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          DO IMU=1,8
C
            DO IMJA=1,2*IABS(KQN(1))
              DO IMJB=1,2*IABS(KQN(2))
                DO IMJC=1,2*IABS(KQN(3))
                  DO IMJD=1,2*IABS(KQN(4))
                    BCOF(LTEN,IMU,IMJA,IMJB,IMJC,IMJD)
     &               = BCOF(LTEN,IMU,IMJA,IMJB,IMJC,IMJD)
     &               +          AJJN(IMJA,IMJB,IMJC,IMJD)*SCOEFF(IMU,1)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
C
          ENDDO
C
C         INCREASE TENSOR LIST LENGTH IF ENTRIES ARE NON-ZERO
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SECOND CONTRIBUTION BELONGS IN "NU+1" BIN
          NUS(LTEN) = NU+1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          DO IMU=1,8
C
            DO IMJA=1,2*IABS(KQN(1))
              DO IMJB=1,2*IABS(KQN(2))
                DO IMJC=1,2*IABS(KQN(3))
                  DO IMJD=1,2*IABS(KQN(4))
                    BCOF(LTEN,IMU,IMJA,IMJB,IMJC,IMJD)
     &               = BCOF(LTEN,IMU,IMJA,IMJB,IMJC,IMJD)
     &         +                AJJN(IMJA,IMJB,IMJC,IMJD)*SCOEFF(IMU,2)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
C
          ENDDO
C
        ENDIF
C
      ENDDO
C
C     RESET THE LOWER AND UPPER LIMITS OF THE EXPANSION
      NUI = NUS(1)
      NUF = NUS(LTEN)
      NUNUM = LTEN
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
C  COEFFICIENTS WE CONFIRMED BY REFERENCING (CHEN AND JOHNSON 1993).   C
C -------------------------------------------------------------------- C
C  MAPPING OF SCOEFF ADDRESSES TO BREIT INTEGRALS (NU±1)               C
C   1: LSLS (UPPER)    2: LSLS (LOWER)                                 C
C   3: SLSL (UPPER)    4: SLSL (LOWER)                                 C
C   5: LSSL (UPPER)    6: LSSL (LOWER)                                 C
C   7: SLLS (UPPER)    8: SLLS (LOWER)                                 C
C**********************************************************************C
      DIMENSION SCOEFF(8,2),KQN(4)
C
      KCA = KQN(2)-KQN(1)
      KDB = KQN(4)-KQN(3)
C
      IF(NU.GE.1) THEN
        B1 = DFLOAT(NU*(NU+1))/DFLOAT(2*NU*(2*NU-1))
        C1 = DFLOAT(   (NU-2))/DFLOAT(2*NU*(2*NU-1))
        SCOEFF(1,1) = DFLOAT(NU  +KCA)*(B1-C1*KDB)
        SCOEFF(2,1) = DFLOAT(NU  +KDB)*(B1-C1*KCA)
        SCOEFF(3,1) = DFLOAT(NU  -KCA)*(B1+C1*KDB)
        SCOEFF(4,1) = DFLOAT(NU  -KDB)*(B1+C1*KCA)
        SCOEFF(5,1) =-DFLOAT(NU  +KCA)*(B1+C1*KDB)
        SCOEFF(6,1) =-DFLOAT(NU  -KDB)*(B1-C1*KCA)
        SCOEFF(7,1) =-DFLOAT(NU  -KCA)*(B1-C1*KDB)
        SCOEFF(8,1) =-DFLOAT(NU  +KDB)*(B1+C1*KCA)
      ELSE
        DO I=1,8
          SCOEFF(I,1) = 0.0D0
        ENDDO
      ENDIF
C
      IF(NU.GE.0) THEN
        B2 = DFLOAT(NU*(NU+1))/DFLOAT((2*NU+2)*(2*NU+3))
        C2 = DFLOAT(    NU+3 )/DFLOAT((2*NU+2)*(2*NU+3))
        SCOEFF(1,2) =-DFLOAT(NU+1-KDB)*(B2+C2*KCA)
        SCOEFF(2,2) =-DFLOAT(NU+1-KCA)*(B2+C2*KDB)
        SCOEFF(3,2) =-DFLOAT(NU+1+KDB)*(B2-C2*KCA)
        SCOEFF(4,2) =-DFLOAT(NU+1+KCA)*(B2-C2*KDB)
        SCOEFF(5,2) = DFLOAT(NU+1+KDB)*(B2+C2*KCA)
        SCOEFF(6,2) = DFLOAT(NU+1-KCA)*(B2-C2*KDB)
        SCOEFF(7,2) = DFLOAT(NU+1-KDB)*(B2-C2*KCA)
        SCOEFF(8,2) = DFLOAT(NU+1+KCA)*(B2+C2*KDB)
      ELSE
        DO I=1,8
          SCOEFF(I,2) = 0.0D0
        ENDDO
      ENDIF
C
      RETURN
      END






      SUBROUTINE BREIT1(IZ)
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
C  BREIT1 CONSTRUCTS A ONE-CENTRE CONTRIBUTION TO THE MOLECULAR        C
C  MEAN-FIELD BREIT MATRIX WITH RACAH ALGEBRA AND BETA INTEGRALS.      C
C  THIS ROUTINE IS SIMILAR TO THE IN-LINE CONSTRUCTION OF MEAN-FIELD   C
C  ATOMIC BREIT MATRIX IN HFSCF0, BUT WITH MQN STRUCTURE AS WELL.      C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*16 HMS
C
      DIMENSION ANGFAC(8)
      DIMENSION BCOF(MNU,8,MKP+1,MKP+1,MKP+1,MKP+1)
      DIMENSION RJLSLS(MB2,MNU,2),RJSLSL(MB2,MNU,2),
     &          RJLSSL(MB2,MNU,2),RJSLLS(MB2,MNU,2)
      DIMENSION XLSLS(MB2),XSLSL(MB2),XLSSL(MB2),XSLLS(MB2)
C
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VANM(MDM,MDM),
     &           VSLF(MDM,MDM),VUEH(MDM,MDM),VWKR(MDM,MDM),
     &           VKSB(MDM,MDM),QDIR(MDM,MDM),QXCH(MDM,MDM),
     &           WDIR(MDM,MDM),WXCH(MDM,MDM),CPLE(MDM,MDM)
C
      COMMON/B0IJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ,IBAS,JBAS
      COMMON/B1QN/EXL(MBS,4),MQN(4),KQN(4),LQN(4),NBAS(4),MAXAB,MAXCD
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/DENS/DENC,DENO,DENT
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VANM,VSLF,
     &            VUEH,VWKR,VKSB,QDIR,QXCH,WDIR,WXCH,CPLE
      COMMON/RCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/SHLL/ACFF,BCFF,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
      COMMON/T2EL/F2ES(5,9),T2ES(5,9),N2EB(5,9),N2EI(5,9),N2ES(5,9)
      COMMON/TSCF/TC1I,TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRW,TCC1,TCC2,
     &            TCMC,TB1I,TB1B,TB1R,TB1F,TB1M,TBEC,TBRM,TBRW,TBC1,
     &            TBC2,TBMC,TSMX,TUMX,THMX,TAMX,TC1T,TC2T,TCVT,TB1T,
     &            TB2T,TACC,TEIG,TSCR,TTOT,TC1S,TC2S,TB1S,TB2S
      COMMON/XNUS/INU(MNU,16),NUS(MNU),NUI,NUF,NUNUM,K4AD
C
C     ANGULAR FACTOR SENSITIVITY PARAMETER
      DATA SENS/1.0D-10/
C
      CALL CPU_TIME(TBCH1)
C
C     INTEGRAL SKIPPING ON MOLECULAR GROUP SYMMETRY CLASS BASIS
      IF(SHAPE.EQ.'ATOMIC') THEN
        ISYM = 2
      ELSEIF(SHAPE.EQ.'DIATOM'.OR.SHAPE.EQ.'LINEAR') THEN
        ISYM = 1
      ELSE
        ISYM = 0
      ENDIF
C
C**********************************************************************C
C     LOOP OVER ALL LQN ORBITAL TYPES (USE INDEX 1000)                 C
C**********************************************************************C
C
C     LOOP OVER LQN(A) VALUES
      DO 1000 LA=0,(NKAP(IZ)-1)/2
C
C       QUANTUM NUMBERS FOR BLOCK A
        LQN(1) = LA
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),IZ)
        DO IBAS=1,NBAS(1)
          EXL(IBAS,1) = BEXL(IBAS,LQN(1),IZ)
        ENDDO
C
C     LOOP OVER LQN(B) VALUES
      DO 1000 LB=0,(NKAP(IZ)-1)/2
C
C       QUANTUM NUMBERS FOR BLOCK B
        LQN(2) = LB
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),IZ)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),IZ)
        ENDDO
C
C       NUMBER OF BASIS FUNCTIONS IN (AB) BLOCK
        MAXAB = NBAS(1)*NBAS(2)
C
C     LOOP OVER LQN(C) VALUES
      DO 1000 LC=0,(NKAP(IZ)-1)/2
C
C       QUANTUM NUMBERS FOR BLOCK C
        LQN(3) = LC
C
C       BASIS EXPONENTS FOR BLOCK C
        NBAS(3) = NFNC(LQN(3),IZ)
        DO KBAS=1,NBAS(3)
          EXL(KBAS,3) = BEXL(KBAS,LQN(3),IZ)
        ENDDO
C
C     LOOP OVER LQN(D) VALUES
      DO 1000 LD=0,(NKAP(IZ)-1)/2
C
C       QUANTUM NUMBERS FOR BLOCK D
        LQN(4) = LD
C
C       BASIS EXPONENTS FOR BLOCK D
        NBAS(4) = NFNC(LQN(4),IZ)
        DO LBAS=1,NBAS(4)
          EXL(LBAS,4) = BEXL(LBAS,LQN(4),IZ)
        ENDDO
C
C       NUMBER OF BASIS FUNCTIONS IN (CD) BLOCK
        MAXCD = NBAS(3)*NBAS(4)
C
C       DETERMINE THE TENSOR ORDERS REQUIRED FOR THIS LQN BLOCK
        CALL TNSBRT1(LQN,ISEL)
C
        IF(ISEL.EQ.0) GOTO 1001
C
C       BASIS SET INTERMEDIATES FOR THE KL-PAIRS
        CALL CPU_TIME(T1)
        CALL KLSET1(IZ)
        CALL CPU_TIME(T2)
        TB1I = TB1I+T2-T1
C
C**********************************************************************C
C     LOOP OVER ALL KQN SYMMETRY TYPES FOR THESE LQNS (USE INDEX 2000) C
C**********************************************************************C
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 NA=1,KRONECK(LA,0),-1
        KA = 2*LA+NA
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,IZ)
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 NB=1,KRONECK(LB,0),-1
        KB = 2*LB+NB
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KAPA(KB,IZ)
C
C     LOOP OVER KQN(C) VALUES
      DO 2000 NC=1,KRONECK(LC,0),-1
        KC = 2*LC+NC
C
C       QUANTUM NUMBERS FOR BLOCK C
        KQN(3) = KAPA(KC,IZ)
C
C     LOOP OVER KQN(D) VALUES
      DO 2000 ND=1,KRONECK(LD,0),-1
        KD = 2*LD+ND
C
C       QUANTUM NUMBERS FOR BLOCK D
        KQN(4) = KAPA(KD,IZ)
C
C     UNIQUE ADDRESS FOR THIS KQN COMBINATION WITHIN THE LQN BLOCK
      K4AD = 16 - (8*NA + 4*NB + 2*NC + ND)
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
      T1100 = TI*TJ
      T1010 = TI*TK
      T1001 = TI*TL
      T0110 = TJ*TK
      T0101 = TJ*TL
      T0011 = TK*TL
      T1110 = TI*TJ*TK
      T1101 = TI*TJ*TL
      T1011 = TI*TK*TL
      T0111 = TJ*TK*TL
      T1111 = TI*TJ*TK*TL
C
C     ANGULAR COEFFICIENTS
      CALL CPU_TIME(T1)
      CALL ANGBRT1(BCOF,KQN,LQN,ISEL)
      CALL CPU_TIME(T2)
      TB1I = TB1I+T2-T1
C
C     EXIT THIS COMBINATION IF IT VIOLATES A SELECTION RULE
      IF(ISEL.EQ.0) GOTO 2001
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS IN BLOCKS A AND B (INDEX 3000)         C
C**********************************************************************C
C
C     UPDATE COUNTER FOR NUMBER OF CLASSES
      N2EB(1,5) = N2EB(1,5)+1
C
      DO 3000 IBAS=1,NBAS(1)
      DO 3000 JBAS=1,NBAS(2)
C
C       UPDATE COUNTER FOR NUMBER OF INTEGRALS
        N2EI(1,5) = N2EI(1,5) + NBAS(3)*NBAS(4)
C
C       BASIS SET INTERMEDIATES FOR THE IJ-PAIRS
        CALL CPU_TIME(T1)
        CALL IJSET1
        CALL CPU_TIME(T2)
        TB1I = TB1I+T2-T1
C
        IF(K4AD.EQ.1) THEN
          CALL CPU_TIME(T1)
          CALL BTABRT1
          CALL CPU_TIME(T2)
          TB1B = TB1B+T2-T1
        ENDIF
C
C       BATCH OF RADIAL INTEGRALS (EFFECTIVE INTERACTION STRENGTHS)
        CALL CPU_TIME(T1)
        CALL RKBRT1(RJLSLS,RJSLSL,RJLSSL,RJSLLS)
        CALL CPU_TIME(T2)
        TB1R = TB1R+T2-T1
C
C       CAN NOW CONTRACT THESE RADIAL INTEGRALS OVER ANGULAR COMPONENTS
C       OF G-SPINOR BASIS FUNCTIONS USING A TENSOR EXPANSION IN {L,Q}
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
C**********************************************************************C
C     LOOP OVER THE SIGNS OF |MQN| FOR A,B,C,D BLOCKS (USE INDEX 5000) C
C**********************************************************************C
C
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
      NAL = LRGE(IZ,KA,IMJA)
      NBL = LRGE(IZ,KB,IMJB)
      NCL = LRGE(IZ,KC,IMJC)
      NDL = LRGE(IZ,KD,IMJD)
C
      NAS = LRGE(IZ,KA,IMJA) + NSKP
      NBS = LRGE(IZ,KB,IMJB) + NSKP
      NCS = LRGE(IZ,KC,IMJC) + NSKP
      NDS = LRGE(IZ,KD,IMJD) + NSKP
C
C     APPLY ANGULAR MQN SELECTION RULE
      IF(MMJA-MMJB.NE.MMJD-MMJC) GOTO 5001
      IF(ISYM.EQ.1.OR.ISYM.EQ.2) THEN
        IF(ISGN1.EQ.ISGN2.AND.ISGN3.EQ.ISGN4.AND.NOPN.NE.0) GOTO 5002
        IF(ISGN1.EQ.ISGN4.AND.ISGN3.EQ.ISGN2) GOTO 5002
        GOTO 5001
      ENDIF
5002  CONTINUE
C
C     RESET CONTRACTED RADIAL ARRAYS
      CALL CPU_TIME(T1)
      DO M=1,NBAS(3)*NBAS(4)
        XLSLS(M) = 0.0D0
C       XSLSL(M) = 0.0D0
        XLSSL(M) = 0.0D0
        XSLLS(M) = 0.0D0
      ENDDO
C
C     ASSEMBLE INTEGRAL BY SUMMING OVER EXPANSION POWER L
      DO LTEN=1,NUNUM
C
C       SKIP IF THIS LQN BLOCK TENSOR ORDER ISN'T IN THE KQN BLOCK
        IF(INU(LTEN,K4AD).EQ.0) GOTO 5003
C
C       ANGULAR COEFFICIENTS
        DO IMU=1,8
          ANGFAC(IMU) = BCOF(LTEN,IMU,IMJA,IMJB,IMJC,IMJD)
        ENDDO
C
C       SCREENING OF INTEGRAL BASED ON ANGULAR COEFFICIENT
        ANGSUM = 0.0D0
        DO IMU=1,8
          ANGSUM = ANGSUM + DABS(ANGFAC(IMU))
        ENDDO
        IF(DABS(ANGSUM).LE.SENS) GOTO 5003
C
C       FACTORS LEADING TO BXCH(LS)
        DO M=1,NBAS(3)*NBAS(4)
          XLSLS(M) = XLSLS(M) + ANGFAC(1)*RJLSLS(M,LTEN,1)
     &                        + ANGFAC(2)*RJLSLS(M,LTEN,2)
        ENDDO
C
C       FACTORS LEADING TO BXCH(SL)
C        IF(NOPN.NE.0) THEN
C          DO M=1,NBAS(3)*NBAS(4)
C            XSLSL(M) = XSLSL(M) + ANGFAC(3)*RJSLSL(M,LTEN,1)
C     &                          + ANGFAC(4)*RJSLSL(M,LTEN,2)
C          ENDDO
C        ENDIF
C
C       FACTORS LEADING TO BXCH(LL)
        IF(NAL.LE.NDL) THEN
          DO M=1,NBAS(3)*NBAS(4)
            XLSSL(M) = XLSSL(M) + ANGFAC(5)*RJLSSL(M,LTEN,1)
     &                          + ANGFAC(6)*RJLSSL(M,LTEN,2)
          ENDDO
        ENDIF
C
C       FACTORS LEADING TO BXCH(SS)
        IF(NAS.LE.NDS) THEN
          DO M=1,NBAS(3)*NBAS(4)
            XSLLS(M) = XSLLS(M) + ANGFAC(7)*RJSLLS(M,LTEN,1)
     &                          + ANGFAC(8)*RJSLLS(M,LTEN,2)
          ENDDO
        ENDIF
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
C     ADD THIS BATCH OF R-INTEGRALS TO CLOSED-SHELL BREIT MATRIX       C
C -------------------------------------------------------------------- C
C          (NO BDIR CONTRIBUTIONS FOR CLOSED-SHELL SYSTEMS.)           C
C**********************************************************************C
C
C     CLOSED-SHELL EXCHANGE MATRIX BXCH
      IF(ISYM.EQ.1.OR.ISYM.EQ.2) THEN
        IF(ISGN1.NE.ISGN4.OR.ISGN3.NE.ISGN2) GOTO 6001
      ENDIF
C
C     CLOSED-SHELL EXCHANGE MATRIX BLOCK BXCH(LL)
      IF(NAL.LE.NDL) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            BXCH(NAL+IBAS,NDL+LBAS) = BXCH(NAL+IBAS,NDL+LBAS)
     &                  +      XLSSL(M)*DCONJG(DENT(NBS+JBAS,NCS+KBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     CLOSED-SHELL EXCHANGE MATRIX BLOCK BXCH(LS)
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
          BXCH(NAL+IBAS,NDS+LBAS) = BXCH(NAL+IBAS,NDS+LBAS)
     &                  +      XLSLS(M)*DCONJG(DENT(NBS+JBAS,NCL+KBAS))
        ENDDO
      ENDDO
C
CC     CLOSED-SHELL EXCHANGE MATRIX BLOCK BXCH(SL)
C      M = 0
C      DO KBAS=1,NBAS(3)
C        DO LBAS=1,NBAS(4)
C          M = M+1
C          BXCH(NAS+IBAS,NDL+LBAS) = BXCH(NAS+IBAS,NDL+LBAS)
C     &                  +      XSLSL(M)*DCONJG(DENT(NBL+JBAS,NCS+KBAS))
C        ENDDO
C      ENDDO
C
C     CLOSED-SHELL EXCHANGE MATRIX BLOCK BXCH(SS)
      IF(NAS.LE.NDS) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            BXCH(NAS+IBAS,NDS+LBAS) = BXCH(NAS+IBAS,NDS+LBAS)
     &                  +      XSLLS(M)*DCONJG(DENT(NBL+JBAS,NCL+KBAS))
          ENDDO
        ENDDO
      ENDIF
6001  CONTINUE
C
C**********************************************************************C
C     ADD THIS BATCH OF R-INTEGRALS TO OPEN-SHELL BREIT MATRIX.        C
C     THIS ALSO REQUIRES THE CLOSED-SHELL DIRECT MATRIX ELEMENTS.      C
C**********************************************************************C
C
      IF(NOPN.EQ.0) GOTO 6100
C
C     CLOSED-SHELL DIRECT MATRIX BDIR
      IF(ISYM.EQ.1.OR.ISYM.EQ.2) THEN
        IF(ISGN1.NE.ISGN2.OR.ISGN3.NE.ISGN4) GOTO 6002
      ENDIF
C
C     CLOSED-SHELL DIRECT MATRIX BLOCK BDIR(LS)
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
          BDIR(NAL+IBAS,NBS+JBAS) = BDIR(NAL+IBAS,NBS+JBAS)
     &                  +      XLSLS(M)*DCONJG(DENT(NCL+KBAS,NDS+LBAS))
     &                  +      XLSSL(M)*DCONJG(DENT(NCS+KBAS,NDL+LBAS))
C
        ENDDO
      ENDDO
CC
CC     CLOSED-SHELL DIRECT MATRIX BLOCK BDIR(SL)
C      M = 0
C      DO KBAS=1,NBAS(3)
C        DO LBAS=1,NBAS(4)
C          M = M+1
C          BDIR(NAS+IBAS,NBL+JBAS) = BDIR(NAS+IBAS,NBL+JBAS)
C     &                  +      XSLSL(M)*DCONJG(DENT(NCS+KBAS,NDL+LBAS))
C     &                  +      XSLLS(M)*DCONJG(DENT(NCL+KBAS,NDS+LBAS))
C        ENDDO
C      ENDDO
C
C     OPEN-SHELL DIRECT MATRIX WDIR
C
C     OPEN-SHELL DIRECT MATRIX BLOCK WDIR(LS)
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
          WDIR(NAL+IBAS,NBS+JBAS) = WDIR(NAL+IBAS,NBS+JBAS)
     &                  + ACFF*XLSLS(M)*DCONJG(DENO(NCL+KBAS,NDS+LBAS))
     &                  + ACFF*XLSSL(M)*DCONJG(DENO(NCS+KBAS,NDL+LBAS))
        ENDDO
      ENDDO
CC
CC     OPEN-SHELL DIRECT MATRIX BLOCK WDIR(SL)
C      M = 0
C      DO KBAS=1,NBAS(3)
C        DO LBAS=1,NBAS(4)
C          M = M+1
C          WDIR(NAS+IBAS,NBL+JBAS) = WDIR(NAS+IBAS,NBL+JBAS)
C     &                  + ACFF*XSLSL(M)*DCONJG(DENO(NCS+KBAS,NDL+LBAS))
C     &                  + ACFF*XSLLS(M)*DCONJG(DENO(NCL+KBAS,NDS+LBAS))
C        ENDDO
C      ENDDO
C
6002  CONTINUE
C
C     OPEN-SHELL DIRECT MATRIX WXCH
      IF(ISYM.EQ.1.OR.ISYM.EQ.2) THEN
        IF(ISGN1.NE.ISGN4.OR.ISGN3.NE.ISGN2) GOTO 6003
      ENDIF
C
C     OPEN-SHELL EXCHANGE MATRIX BLOCK WXCH(LL)
      IF(NAL.LE.NDL) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            WXCH(NAL+IBAS,NDL+LBAS) = WXCH(NAL+IBAS,NDL+LBAS)
     &                  + BCFF*XLSSL(M)*DCONJG(DENO(NBS+JBAS,NCS+KBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     OPEN-SHELL EXCHANGE MATRIX BLOCK WXCH(LS)
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
          WXCH(NAL+IBAS,NDS+LBAS) = WXCH(NAL+IBAS,NDS+LBAS)
     &                  + BCFF*XLSLS(M)*DCONJG(DENO(NBS+JBAS,NCL+KBAS))
        ENDDO
      ENDDO
CC
CC     OPEN-SHELL EXCHANGE MATRIX BLOCK WXCH(SL)
C      M = 0
C      DO KBAS=1,NBAS(3)
C        DO LBAS=1,NBAS(4)
C          M = M+1
C          WXCH(NAS+IBAS,NDL+LBAS) = WXCH(NAS+IBAS,NDL+LBAS)
C     &                  + BCFF*XSLSL(M)*DCONJG(DENO(NBL+JBAS,NCS+KBAS))
C        ENDDO
C      ENDDO
C
C     OPEN-SHELL EXCHANGE MATRIX BLOCK WXCH(SS)
      IF(NAS.LE.NDS) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            WXCH(NAS+IBAS,NDS+LBAS) = WXCH(NAS+IBAS,NDS+LBAS)
     &                  + BCFF*XSLLS(M)*DCONJG(DENO(NBL+JBAS,NCL+KBAS))
          ENDDO
        ENDDO
      ENDIF
6003  CONTINUE
C
6100  CONTINUE
C
C     MATRIX MULTIPLICATION STEP COMPLETE
      CALL CPU_TIME(T3)
      TB1M = TB1M+T3-T2
C
C**********************************************************************C
C     COMPLETE CONSTRUCTION OF ALL ONE-CENTRE CONTRIBUTIONS.           C
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
C     EARLY EXIT FOR LQN SELECTION RULES
1001  CONTINUE
C     END LOOP OVER ALL LQNS
1000  CONTINUE
C
C     RECORD CPU TIME AT END OF BATCH AND ADD TO APPROPRIATE COUNTER
      CALL CPU_TIME(TBCH2)
      IF(INTSYM) THEN
        T2ES(1,5) = T2ES(1,5)+TBCH2-TBCH1
      ELSE
        T2ES(1,5) = T2ES(1,5) + 0.25D0*(TBCH2-TBCH1)
        T2ES(1,6) = T2ES(1,6) + 0.25D0*(TBCH2-TBCH1)
        T2ES(1,7) = T2ES(1,7) + 0.25D0*(TBCH2-TBCH1)
        T2ES(1,8) = T2ES(1,8) + 0.25D0*(TBCH2-TBCH1)
      ENDIF
C
C     MATRIX HERMITICITY (CLOSED-SHELL)
      DO I=1,NDIM-NSKP
        DO J=NSKP+1,NDIM
          BXCH(J,I) = DCONJG(BXCH(I,J))
        ENDDO
      ENDDO
C
C     MATRIX HERMITICITY (OPEN-SHELL)
      IF(NOPN.NE.0) THEN
        DO I=1,NDIM-NSKP
          DO J=NSKP+1,NDIM
            BDIR(J,I) = DCONJG(BDIR(I,J))
            WDIR(J,I) = DCONJG(WDIR(I,J))
            WXCH(J,I) = DCONJG(WXCH(I,J))
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE BTABRT1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        BBBBBBB TTTTTTTT   AA    BBBBBBB  RRRRRRR TTTTTTTT 11         C
C        BB    BB   TT     AAAA   BB    BB RR    RR   TT   111         C
C        BB    BB   TT    AA  AA  BB    BB RR    RR   TT    11         C
C        BBBBBBB    TT   AA    AA BBBBBBB  RR    RR   TT    11         C
C        BB    BB   TT   AAAAAAAA BB    BB RRRRRRR    TT    11         C
C        BB    BB   TT   AA    AA BB    BB RR    RR   TT    11         C
C        BBBBBBB    TT   AA    AA BBBBBBB  RR    RR   TT   1111        C
C                                                                      C
C -------------------------------------------------------------------- C
C  BTABRT1 EVALUATES A BATCH OF BETA INTEGRALS REQUIRED FOR THE        C
C  ONE-CENTRE RADIAL INTEGRALS OVER THE BREIT INTERACTION.             C
C -------------------------------------------------------------------- C
C  INDEX 1 FOR UPPER INTEGRALS, INDEX 2 FOR LOWER INTEGRALS.           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION XJ(MB2,2),IAA(2),IBB(2)
      DIMENSION BETA(MB2),BTA1(MB2),XROOT(MB2),BIN(MB2),TRM(MB2)
C
      DIMENSION RJLSLS(MB2,MNU,2),RJSLSL(MB2,MNU,2),
     &          RJLSSL(MB2,MNU,2),RJSLLS(MB2,MNU,2)
C
      COMMON/B0IJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ,IBAS,JBAS
      COMMON/B0KL/EKL(MB2,-MAB:MAB),RNKL(MB2,4),EK(MB2),EL(MB2)
      COMMON/B1QN/EXL(MBS,4),MQN(4),KQN(4),LQN(4),NBAS(4),MAXAB,MAXCD
      COMMON/BETA/BDU(MB2*MB2,MNU,MNU),BDL(MB2*MB2,MNU,MNU)
      COMMON/RCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/XNUS/INU(MNU,16),NUS(MNU),NUI,NUF,NUNUM,K4AD
C
C**********************************************************************C
C     PREPARATION OF EXPONENT POWERS AND BETA FUNCTION ARGUMENTS.      C
C**********************************************************************C
C
C     UNIQUE FIXED ADDRESS FOR THE IBAS,JBAS COMBINATION
      MIJ = (IBAS-1)*NBAS(2)+JBAS-1
C
C     GAUSSIAN EXPONENT FOR (IBAS,JBAS) OVERLAP
      TIJ0 = EI+EJ
C
C     INITIALISE THE ARRAY XJ(M,2) FOR INCOMPLETE BETA FUNCTION ARGS.
      DO M=1,MAXCD
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
C     NUMBER OF LEVELS NEEDED TO ACCOUNT FOR ALL BETA INTEGRALS
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
              DO M=1,MAXCD
                BTA1(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
              RA  = X
              RB  = DFLOAT(1-IB)
              RC  = X+1.0D0
              RD  = 1.0D0
              RCD = RC*RD
              FCT = RA*RB/RCD
              DO M=1,MAXCD
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
                DO M=1,MAXCD
                  TRM(M) = FCT*TRM(M)*XJ(M,IBETA)
                  BIN(M) = BIN(M)+TRM(M)
                ENDDO
                RA = RA+1.0D0
                RB = RB+1.0D0
                RC = RC+1.0D0
                RD = RD+1.0D0
              ENDDO
              DO M=1,MAXCD
                BETA(M) = BTA1(M)*BIN(M)
              ENDDO
C
C           CASE 2: IB = 1
            ELSEIF(IB.EQ.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXCD
                BETA(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
C
C           CASE 3: IB = 0
            ELSEIF(IB.EQ.0) THEN
              DO M=1,MAXCD
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
                  DO M=1,MAXCD
                    BIN(M) = BIN(M)+X*TRM(M)
                    TRM(M) = TRM(M)*XJ(M,IBETA)
                  ENDDO
                ENDDO
                DO M=1,MAXCD
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)*BIN(M)
                ENDDO
              ELSEIF(IA.EQ.1) THEN
                DO M=1,MAXCD
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)
                ENDDO
              ELSE
                DO M=1,MAXCD
                  BETA(M) = BTA1(M)
                ENDDO
              ENDIF
C
C           END CONDITIONAL STATEMENT OVER IB VALUES
            ENDIF
C
            IF(IBETA.EQ.1) THEN
              DO M=1,MAXCD
                BDU(MIJ*MAXCD+M,NX,NY) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXCD
                BDL(MIJ*MAXCD+M,NX,NY) = BETA(M)
              ENDDO
            ENDIF

          ENDDO

        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE RKBRT1(RJLSLS,RJSLSL,RJLSSL,RJSLLS)
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
C  RKBRT1 EVALUATES A BATCH OF GENERAL (OPEN-SHELL) ONE-CENTRE RADIAL  C
C  INTEGRALS OVER THE BREIT INTERACTION, SEPARATING RESULTS BY THE     C
C  ALLOWED TENSOR ORDERS RK(ABCD).                                     C
C -------------------------------------------------------------------- C
C  RJLSLS(MB2,MNU,IUL) - JLSLS INTEGRAL LIST OF TENSOR TYPE MNU        C
C  RJSLSL(MB2,MNU,IUL) - JSLSL INTEGRAL LIST OF TENSOR TYPE MNU        C
C  RJLSSL(MB2,MNU,IUL) - JLSSL INTEGRAL LIST OF TENSOR TYPE MNU        C
C  RJSLLS(MB2,MNU,IUL) - JSLLS INTEGRAL LIST OF TENSOR TYPE MNU        C
C -------------------------------------------------------------------- C
C  INDEX 1 FOR UPPER INTEGRALS, INDEX 2 FOR LOWER INTEGRALS.           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RJLSLS(MB2,MNU,2),RJSLSL(MB2,MNU,2),
     &          RJLSSL(MB2,MNU,2),RJSLLS(MB2,MNU,2)
C
      COMMON/B0IJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ,IBAS,JBAS
      COMMON/B0KL/EKL(MB2,-MAB:MAB),RNKL(MB2,4),EK(MB2),EL(MB2)
      COMMON/B1QN/EXL(MBS,4),MQN(4),KQN(4),LQN(4),NBAS(4),MAXAB,MAXCD
      COMMON/BETA/BDU(MB2*MB2,MNU,MNU),BDL(MB2*MB2,MNU,MNU)
      COMMON/RCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/XNUS/INU(MNU,16),NUS(MNU),NUI,NUF,NUNUM,K4AD
C
C**********************************************************************C
C     RADIAL INTEGRALS OVER SPINORS (SEPARATED BY TENSOR ORDER)        C
C**********************************************************************C
C
C     EMPTY COUNTER ARRAYS FOR DIRECT AND EXCHANGE INTEGRALS
      DO M=1,MAXCD
        DO LTEN=1,NUNUM
          DO IUL=1,2
            RJLSLS(M,LTEN,IUL) = 0.0D0
C           RJSLSL(M,LTEN,IUL) = 0.0D0
            RJLSSL(M,LTEN,IUL) = 0.0D0
            RJSLLS(M,LTEN,IUL) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     VALUE PREPARATION
      E0000 = 1.0D0
      E1000 = EI
      E0100 = EJ
      E1100 = EI*EJ
C
C     LOOP OVER K,L BASIS FUNCTIONS
      DO M=1,MAXCD
C
C       MORE VALUE PREPARATION
        E0010 = EK(M)
        E0001 = EL(M)
        E1010 = EI*EK(M)
        E1001 = EI*EL(M)
        E0110 = EJ*EK(M)
        E0101 = EJ*EL(M)
        E0011 = EK(M)*EL(M)
C
C       UNIQUE FIXED ADDRESS FOR THE IBAS,JBAS COMBINATION
        MIJKL = ((IBAS-1)*NBAS(2)+JBAS-1)*MAXCD + M
C
C       LOOP OVER THE TENSOR ORDERS OF THE BREIT INTERACTION
        DO LTEN=1,NUNUM
C
C         SCREEN LQN BLOCK TENSOR ORDERS THAT ARE NOT IN THIS KQN BLOCK
          IF(INU(LTEN,K4AD).EQ.0) GOTO 50
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(LTEN)
C
C         TEMPORARY STORAGE OF RAW RJ(LTEN,M)
          B21L = EIJ(-NU+1)*EKL(M, NU+2)*BDL(MIJKL,LTEN  ,NUNUM-LTEN+1)
          B21U = EIJ( NU+2)*EKL(M,-NU+1)*BDU(MIJKL,LTEN  ,NUNUM-LTEN+1)
          B23L = EIJ(-NU+3)*EKL(M, NU+2)*BDL(MIJKL,LTEN  ,NUNUM-LTEN+2)
          B41U = EIJ( NU+4)*EKL(M,-NU+1)*BDU(MIJKL,LTEN+1,NUNUM-LTEN+1)
          B41L = EIJ(-NU+1)*EKL(M, NU+4)*BDL(MIJKL,LTEN+1,NUNUM-LTEN+1)
          B23U = EIJ( NU+2)*EKL(M,-NU+3)*BDU(MIJKL,LTEN  ,NUNUM-LTEN+2)
          B43L = EIJ(-NU+3)*EKL(M, NU+4)*BDL(MIJKL,LTEN+1,NUNUM-LTEN+2)
          B43U = EIJ( NU+4)*EKL(M,-NU+3)*BDU(MIJKL,LTEN+1,NUNUM-LTEN+2)
C
C         EFFECTIVE INTERACTION STRENGTH RJLSLS(M,LTEN)
          RJLSLS(M,LTEN,1) 
     &         = V4*T0000*E0101*C7*B43U - V2*T0001*E0100*C5*B41U
     &         - V2*T0100*E0001*C5*B23U + V1*T0101*E0000*C3*B21U
C
          RJLSLS(M,LTEN,2) 
     &         = V4*T0000*E0101*C7*B43L - V2*T0001*E0100*C5*B23L
     &         - V2*T0100*E0001*C5*B41L + V1*T0101*E0000*C3*B21L
CC
CC         EFFECTIVE INTERACTION STRENGTH RJSLSL(M,LTEN)
C          RJSLSL(M,LTEN,1)
C     &         = V4*T0000*E1010*C7*B43U - V2*T0010*E1000*C5*B41U
C     &         - V2*T1000*E0010*C5*B23U + V1*T1010*E0000*C3*B21U
CC
C          RJSLSL(M,LTEN,2)
C     &         = V4*T0000*E1010*C7*B43L - V2*T0010*E1000*C5*B23L
C     &         - V2*T1000*E0010*C5*B41L + V1*T1010*E0000*C3*B21L
C
C         EFFECTIVE INTERACTION STRENGTH RJLSSL(M,LTEN)
          RJLSSL(M,LTEN,1) 
     &         = V4*T0000*E0110*C7*B43U - V2*T0010*E0100*C5*B41U
     &         - V2*T0100*E0010*C5*B23U + V1*T0110*E0000*C3*B21U
C
          RJLSSL(M,LTEN,2)
     &         = V4*T0000*E0110*C7*B43L - V2*T0010*E0100*C5*B23L
     &         - V2*T0100*E0010*C5*B41L + V1*T0110*E0000*C3*B21L
C
C         EFFECTIVE INTERACTION STRENGTH RJSLLS(M,LTEN)
          RJSLLS(M,LTEN,1) 
     &         = V4*T0000*E1001*C7*B43U - V2*T0001*E1000*C5*B41U
     &         - V2*T1000*E0001*C5*B23U + V1*T1001*E0000*C3*B21U
C
          RJSLLS(M,LTEN,2) 
     &         = V4*T0000*E1001*C7*B43L - V2*T0001*E1000*C5*B23L
     &         - V2*T1000*E0001*C5*B41L + V1*T1001*E0000*C3*B21L
C
50        CONTINUE
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
      DO M=1,MAXCD
        RNLSLS = RNIJ(2)*RNKL(M,2)
C       RNSLSL = RNIJ(3)*RNKL(M,3)
        RNSLLS = RNIJ(3)*RNKL(M,2)
        RNLSSL = RNIJ(2)*RNKL(M,3)
        DO LTEN=1,NUNUM
          IF(INU(LTEN,K4AD).NE.0) THEN
            DO IUL=1,2
              RJLSLS(M,LTEN,IUL) = RNLSLS*RJLSLS(M,LTEN,IUL)
C             RJSLSL(M,LTEN,IUL) = RNSLSL*RJSLSL(M,LTEN,IUL)
              RJLSSL(M,LTEN,IUL) = RNLSSL*RJLSSL(M,LTEN,IUL)
              RJSLLS(M,LTEN,IUL) = RNSLLS*RJSLLS(M,LTEN,IUL)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE TNSBRT1(LQN,ISEL,IWRIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       TTTTTTTT NN    NN  SSSSSS  BBBBBBB  RRRRRRR TTTTTTTT 11        C
C          TT    NNN   NN SS    SS BB    BB RR    RR   TT   111        C
C          TT    NNNN  NN SS       BB    BB RR    RR   TT    11        C
C          TT    NN NN NN  SSSSSS  BBBBBBB  RR    RR   TT    11        C
C          TT    NN  NNNN       SS BB    BB RRRRRRR    TT    11        C
C          TT    NN   NNN SS    SS BB    BB RR    RR   TT    11        C
C          TT    NN    NN  SSSSSS  BBBBBBB  RR    RR   TT   1111       C
C                                                                      C
C -------------------------------------------------------------------- C
C  TNSBRT1 PERFORMS THE TENSOR EXPANSION ANALYSIS FOR THE ONE-CENTRE   C
C  BREIT INTERACTION, FOR A GIVEN LQN COMBINATION BLOCK.               C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION LQN(4)
C
      COMMON/XNUS/INU(MNU,16),NUS(MNU),NUI,NUF,NUNUM,K4AD
C
C     INITIALISE TENSOR EXPANSION ARRAY
      DO LTEN=1,MNU
        NUS(LTEN) = 0
        DO I=1,16
          INU(LTEN,I) = 0
        ENDDO
      ENDDO
      NUNUM = 0
C
C**********************************************************************C
C     LQN PARITY ANALYSIS: CHECK PARITY OF LQN COMBINATIONS.           C
C     EXIT IF THERE IS NO MULTIPOLE EXPANSION FOR THIS CASE.           C
C**********************************************************************C
C
C     A AND B: PARITY OF 'LQN(1)+LQN(2)' (0 IF EVEN, 1 IF ODD)
      IPARAB = MOD(LQN(1)+LQN(2),2)
C
C     C AND D: PARITY OF 'LQN(3)+LQN(4)' (0 IF EVEN, 1 IF ODD)
      IPARCD = MOD(LQN(3)+LQN(4),2)
C
C     LQN SELECTION RULE: BOTH LQN PAIRS MUST BE OF SAME SYMMETRY
      IF(IPARAB.NE.IPARCD) THEN
        ISEL = 0
        RETURN
      ELSE
        ISEL = 1
      ENDIF
C
C**********************************************************************C
C     MULTIPOLE EXPANSION UPPER/LOWER LIMITS SATISFY TRIANGLE RULE.    C
C**********************************************************************C
C
C     TENSOR LIMITS BY TRIANGLE RULE
      NUI = MAX0(IABS(LQN(1)-LQN(2)),IABS(LQN(3)-LQN(4)))
      NUF = MIN0(    (LQN(1)+LQN(2)),    (LQN(3)+LQN(4)))
C
C     TRIANGLE RULE MUST PROVIDE VALID NU RANGE
      IF(NUI.GT.NUF) THEN
        ISEL = 0
        RETURN
      ELSE
        ISEL = 1
      ENDIF
C
C**********************************************************************C
C     CALCULATE COMPONENT-TYPE OVERLAP COEFFICIENTS FOR TENSOR ORDERS  C
C     NU (AND NU±1 WHERE APPLICABLE), COUPLING IN ANGULAR FACTORS.     C
C**********************************************************************C
C
      ISEL = 0
      LTEN = 1
      DO NU=NUI,NUF
C
C       A AND B: PARITY OF 'LA+LB+NU' (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQN(1)+LQN(2)+NU,2)
C
C       C AND D: PARITY OF 'LC+LD+NU' (0 IF EVEN, 1 IF ODD)
        IPARCD = MOD(LQN(3)+LQN(4)+NU,2)
C
C       CASE 1: LA+LB+NU AND LC+LD+NU ARE BOTH ODD (AND NU=/=0)
        IF(IPARAB.EQ.1.AND.IPARCD.EQ.1.AND.NU.NE.0) THEN
C
C         THERE IS A NON-TRIVIAL TENSOR EXPANSION
          ISEL = 1
C
C         THIS CONTRIBUTION BELONGS IN "NU" BIN
          NUS(LTEN) = NU
C
        ENDIF
C
C       CASE 2: LA+LB+NU AND LC+LD+NU ARE BOTH EVEN
        IF(IPARAB.EQ.0.AND.IPARCD.EQ.0) THEN
C
C         THERE IS A NON-TRIVIAL TENSOR EXPANSION
          ISEL = 1
C
C         FIRST CONTRIBUTION BELONGS IN "NU-1" BIN
          NUS(LTEN) = NU-1
C
C         INCREASE TENSOR LIST LENGTH IF ENTRIES ARE NON-ZERO
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SECOND CONTRIBUTION BELONGS IN "NU+1" BIN
          NUS(LTEN) = NU+1
C
        ENDIF
C
      ENDDO
C
C     RESET THE LOWER AND UPPER LIMITS OF THE EXPANSION
      NUI = NUS(1)
      NUF = NUS(LTEN)
      NUNUM = LTEN
      
      IF(IWRIT.NE.0) THEN
        WRITE(*,*) REPEAT('=',53)
        WRITE(*,*) 'LQN: ',(LQN(N),N=1,4)
        WRITE(*,*) REPEAT('-',53)
        WRITE(*,*) 'ZZ:  ',(NUS(NU),NU=1,LTEN)
        WRITE(*,*) REPEAT('-',53)
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE ANGBRT1(BCOF,KQN,LQN,ISEL,IWRIT)
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
C  WE USE TEMPORARY COUNTERS FOR TENSOR ORDERS AND TRANSFER AT THE END.C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION KQN(4),LQN(4),JQN(4)
      DIMENSION SCOEFF(8,2)
      DIMENSION DKAB(0:MNU,MKP+1,MKP+1),DKCD(0:MNU,MKP+1,MKP+1)
      DIMENSION AJJN(MKP+1,MKP+1,MKP+1,MKP+1)
      DIMENSION BCOF(MNU,8,MKP+1,MKP+1,MKP+1,MKP+1)
C
      COMMON/XNUS/INU(MNU,16),NUS(MNU),NUI,NUF,NUNUM,K4AD
C
C     INITIALISE BREIT COEFFICIENTS
      DO K=1,MNU
        DO IMU=1,8
          DO IMJA=1,2*IABS(KQN(1))
            DO IMJB=1,2*IABS(KQN(2))
              DO IMJC=1,2*IABS(KQN(3))
                DO IMJD=1,2*IABS(KQN(4))
                  BCOF(K,IMU,IMJA,IMJB,IMJC,IMJD) = 0.0D0
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
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
C     TENSOR LIMITS BY TRIANGLE RULE
      NUIT = MAX0(IABS(JQN(1)-JQN(2))/2,IABS(JQN(3)-JQN(4))/2)
      NUFT = MIN0(    (JQN(1)+JQN(2))/2,    (JQN(3)+JQN(4))/2)
C
C     JQN SELECTION RULE: TRIANGLE RULE MUST PROVIDE VALID NU RANGE
      IF(NUIT.GT.NUFT) THEN
        ISEL = 0
        RETURN
      ELSE
        ISEL = 1
      ENDIF
C
C**********************************************************************C
C     LQN PARITY ANALYSIS: CHECK PARITY OF LQN COMBINATIONS.           C
C     EXIT IF THERE IS NO MULTIPOLE EXPANSION FOR THIS CASE.           C
C**********************************************************************C
C
C     A AND B: PARITY OF 'LQN(1)+LQN(2)' (0 IF EVEN, 1 IF ODD)
      IPARAB = MOD(LQN(1)+LQN(2),2)
C
C     C AND D: PARITY OF 'LQN(3)+LQN(4)' (0 IF EVEN, 1 IF ODD)
      IPARCD = MOD(LQN(3)+LQN(4),2)
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
C     WIGNER-ECKART: EVALUATE AN ANGULAR FACTOR FOR EVERY |MQN|        C
C     COMBINATION IN THE MULTIPOLE EXPANSION OVER ALLOWED NU VALUES.   C
C**********************************************************************C
C
C     LOOP OVER MQN(A) AND MQN(B) VALUES
      DO MA=1,IABS(KQN(1))
        MJA = 2*MA-1
        DO MB=1,IABS(KQN(2))
          MJB = 2*MB-1
C
C         LOOP OVER ALL SURVIVING TENSOR ORDERS FOR THIS KQN BLOCK
          DO NU=NUIT,NUFT
C
C           GENERATE AN ANGULAR FACTOR FOR ALL |MQN| SIGNS
            DKAB(NU,MJA  ,MJB  ) = DK(JQN(1),-MJA,JQN(2),-MJB,NU)
            DKAB(NU,MJA+1,MJB  ) = DK(JQN(1), MJA,JQN(2),-MJB,NU)
            DKAB(NU,MJA  ,MJB+1) = DK(JQN(1),-MJA,JQN(2), MJB,NU)
            DKAB(NU,MJA+1,MJB+1) = DK(JQN(1), MJA,JQN(2), MJB,NU)
C
          ENDDO
C
        ENDDO
      ENDDO
C
C     THE ARGUMENTS OF DK COEFFICIENTS ARE REVERSED IN THE CASE OF THE
C     CD PAIRS IN ORDER TO ACCOMMODATE THE RELATION:
C               DK(J,M,J',M',L) = ((-1)^Q)*DK(J',M',J,M,L).
C
C     LOOP OVER MQN(C) AND MQN(D) VALUES
      DO MC=1,IABS(KQN(3))
        MJC = 2*MC-1
        DO MD=1,IABS(KQN(4))
          MJD = 2*MD-1
C
C         LOOP OVER ALL SURVIVING TENSOR ORDERS FOR THIS KQN BLOCK
          DO NU=NUIT,NUFT
C
C           GENERATE AN ANGULAR FACTOR FOR ALL |MQN| SIGNS
            DKCD(NU,MJC  ,MJD  ) = DK(JQN(4),-MJD,JQN(3),-MJC,NU)
            DKCD(NU,MJC+1,MJD  ) = DK(JQN(4),-MJD,JQN(3), MJC,NU)
            DKCD(NU,MJC  ,MJD+1) = DK(JQN(4), MJD,JQN(3),-MJC,NU)
            DKCD(NU,MJC+1,MJD+1) = DK(JQN(4), MJD,JQN(3), MJC,NU)
C
          ENDDO
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CALCULATE COMPONENT-TYPE OVERLAP COEFFICIENTS FOR TENSOR ORDERS  C
C     NU (AND NU±1 WHERE APPLICABLE), COUPLING IN ANGULAR FACTORS.     C
C**********************************************************************C
C
      ISEL = 0
      DO NU=NUIT,NUFT
C
C       A AND B: PARITY OF 'LA+LB+NU' (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQN(1)+LQN(2)+NU,2)
C
C       C AND D: PARITY OF 'LC+LD+NU' (0 IF EVEN, 1 IF ODD)
        IPARCD = MOD(LQN(3)+LQN(4)+NU,2)
C
C       JQN AND NU ANGULAR FACTORS FOR THIS COMBINATION
        DO IMJA=1,2*IABS(KQN(1))
          DO IMJB=1,2*IABS(KQN(2))
            DO IMJC=1,2*IABS(KQN(3))
              DO IMJD=1,2*IABS(KQN(4))
                AJJN(IMJA,IMJB,IMJC,IMJD)
     &           = DKAB(NU,IMJA,IMJB)*DKCD(NU,IMJC,IMJD)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C
C       CASE 1: LA+LB+NU AND LC+LD+NU ARE BOTH ODD (AND NU=/=0)
        IF(IPARAB.EQ.1.AND.IPARCD.EQ.1.AND.NU.NE.0) THEN
C
C         THERE IS A NON-TRIVIAL TENSOR EXPANSION
          ISEL = 1
C
C         FIND THE APPROPRIATE INDEX K FOR THIS NU
          K = KFIND(NU)
          INU(K,K4AD) = 1
C
C         INTERMEDIATE (KQN,NU) COEFFICIENTS
          RNU    = DFLOAT(NU*(NU+1))
          RCOEFF =-DFLOAT((KQN(1)+KQN(2))*(KQN(3)+KQN(4)))/RNU
C
C         CONTRIBUTIONS TO TENSOR ORDER "NU"
          DO IMU=1,8
C
            DO IMJA=1,2*IABS(KQN(1))
              DO IMJB=1,2*IABS(KQN(2))
                DO IMJC=1,2*IABS(KQN(3))
                  DO IMJD=1,2*IABS(KQN(4))
                    BCOF(K,IMU,IMJA,IMJB,IMJC,IMJD)
     &               = BCOF(K,IMU,IMJA,IMJB,IMJC,IMJD)
     &               +          AJJN(IMJA,IMJB,IMJC,IMJD)*RCOEFF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
C
          ENDDO
C
        ENDIF
C
C       CASE 2: LA+LB+NU AND LC+LD+NU ARE BOTH EVEN
        IF(IPARAB.EQ.0.AND.IPARCD.EQ.0) THEN
C
C         THERE IS A NON-TRIVIAL TENSOR EXPANSION
          ISEL = 1
C
C         CALCULATE BREIT COEFFICIENTS WITH CALL TO BRCOEF1
          CALL BRCOEF1(SCOEFF,KQN,NU)
C
C         FIRST CONTRIBUTION BELONGS IN "NU-1" BIN
          IF(NU.GT.0) THEN
            K = KFIND(NU-1)
            INU(K,K4AD) = 1
          ENDIF
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          DO IMU=1,8
C
            DO IMJA=1,2*IABS(KQN(1))
              DO IMJB=1,2*IABS(KQN(2))
                DO IMJC=1,2*IABS(KQN(3))
                  DO IMJD=1,2*IABS(KQN(4))
                    BCOF(K,IMU,IMJA,IMJB,IMJC,IMJD)
     &               = BCOF(K,IMU,IMJA,IMJB,IMJC,IMJD)
     &               +          AJJN(IMJA,IMJB,IMJC,IMJD)*SCOEFF(IMU,1)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
C
          ENDDO
C
C         SECOND CONTRIBUTION BELONGS IN "NU+1" BIN
          K = KFIND(NU+1)
          INU(K,K4AD) = 1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          DO IMU=1,8
C
            DO IMJA=1,2*IABS(KQN(1))
              DO IMJB=1,2*IABS(KQN(2))
                DO IMJC=1,2*IABS(KQN(3))
                  DO IMJD=1,2*IABS(KQN(4))
                    BCOF(K,IMU,IMJA,IMJB,IMJC,IMJD)
     &               = BCOF(K,IMU,IMJA,IMJB,IMJC,IMJD)
     &         +                AJJN(IMJA,IMJB,IMJC,IMJD)*SCOEFF(IMU,2)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
C
          ENDDO
C
        ENDIF
C
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION KFIND(NU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C               KK    KK FFFFFFFF IIII NN    NN DDDDDDD                C
C               KK   KK  FF        II  NNN   NN DD    DD               C
C               KK  KK   FF        II  NNNN  NN DD    DD               C
C               KKKKK    FFFFFF    II  NN NN NN DD    DD               C
C               KK  KK   FF        II  NN  NNNN DD    DD               C
C               KK   KK  FF        II  NN   NNN DD    DD               C
C               KK    KK FF       IIII NN    NN DDDDDDD                C
C                                                                      C
C -------------------------------------------------------------------- C
C  KFIND DETERMINES THE TENSOR INDEX FOR A GIVEN VALUE N (LQN BLOCK).  C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5 STRING
      COMMON/XNUS/INU(MNU,16),NUS(MNU),NUI,NUF,NUNUM,K4AD
C
      KFIND = 0
      DO K=1,NUNUM
        IF(NU.EQ.NUS(K)) THEN
          KFIND = K
          GOTO 50
        ENDIF
      ENDDO
      WRITE(6, *) 'In KFIND: could not find tensor index.',NU
      WRITE(7, *) 'In KFIND: could not find tensor index.',NU
50    CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE BTACLM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      BBBBBBB TTTTTTTT   AA     CCCCCC  LL       MM       MM  11      C
C      BB    BB   TT     AAAA   CC    CC LL       MMM     MMM 111      C
C      BB    BB   TT    AA  AA  CC       LL       MMMM   MMMM  11      C
C      BBBBBBB    TT   AA    AA CC       LL       MM MM MM MM  11      C
C      BB    BB   TT   AAAAAAAA CC       LL       MM  MMM  MM  11      C
C      BB    BB   TT   AA    AA CC    CC LL       MM   M   MM  11      C
C      BBBBBBB    TT   AA    AA  CCCCCC  LLLLLLLL MM       MM 1111     C
C                                                                      C
C -------------------------------------------------------------------- C
C  BTACLM1 EVALUATES A BATCH OF BETA INTEGRALS REQUIRED FOR THE        C
C  ONE-CENTRE RADIAL INTEGRALS OVER THE COULOMB INTERACTION.           C
C -------------------------------------------------------------------- C
C  INDEX 1 FOR UPPER INTEGRALS, INDEX 2 FOR LOWER INTEGRALS.           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION XJ(MB2,2),IAA(2),IBB(2)
      DIMENSION BETA(MB2),BTA1(MB2),XROOT(MB2),BIN(MB2),TRM(MB2)
C
      DIMENSION RJLSLS(MB2,MNU,2),RJSLSL(MB2,MNU,2),
     &          RJLSSL(MB2,MNU,2),RJSLLS(MB2,MNU,2)
C
      COMMON/B0IJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ,IBAS,JBAS
      COMMON/B0KL/EKL(MB2,-MAB:MAB),RNKL(MB2,4),EK(MB2),EL(MB2)
      COMMON/B1QN/EXL(MBS,4),MQN(4),KQN(4),LQN(4),NBAS(4),MAXAB,MAXCD
      COMMON/BETA/BDU(MB2*MB2,MNU,MNU),BDL(MB2*MB2,MNU,MNU)
      COMMON/RCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/XNUS/INU(MNU,16),NUS(MNU),NUI,NUF,NUNUM,K4AD
C
C**********************************************************************C
C     PREPARATION OF EXPONENT POWERS AND BETA FUNCTION ARGUMENTS.      C
C**********************************************************************C
C
C     UNIQUE FIXED ADDRESS FOR THE IBAS,JBAS COMBINATION
      MIJ = (IBAS-1)*NBAS(2)+JBAS-1
C
C     GAUSSIAN EXPONENT FOR (IBAS,JBAS) OVERLAP
      TIJ0 = EI+EJ
C
C     INITIALISE THE ARRAY XJ(M,2) FOR INCOMPLETE BETA FUNCTION ARGS.
      DO M=1,MAXCD
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
      NVALS = (NUF-NUI)/2+3
C
C     LOOP OVER EXPANSION TERMINALS FOR FIRST PAIR
      DO NX=1,NVALS
        IAA(1) = LQN(1)+LQN(2)+NUI+2*NX-1
        IAA(2) = LQN(3)+LQN(4)+NUI+2*NX-1
C
C       LOOP OVER EXPANSION TERMINALS FOR SECOND PAIR
        DO NY=1,NVALS
          IBB(1) = LQN(3)+LQN(4)-NUF+2*NY-2
          IBB(2) = LQN(1)+LQN(2)-NUF+2*NY-2
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
              DO M=1,MAXCD
                BTA1(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
              RA  = X
              RB  = DFLOAT(1-IB)
              RC  = X+1.0D0
              RD  = 1.0D0
              RCD = RC*RD
              FCT = RA*RB/RCD
              DO M=1,MAXCD
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
                DO M=1,MAXCD
                  TRM(M) = FCT*TRM(M)*XJ(M,IBETA)
                  BIN(M) = BIN(M)+TRM(M)
                ENDDO
                RA = RA+1.0D0
                RB = RB+1.0D0
                RC = RC+1.0D0
                RD = RD+1.0D0
              ENDDO
              DO M=1,MAXCD
                BETA(M) = BTA1(M)*BIN(M)
              ENDDO
C
C           CASE 2: IB = 1
            ELSEIF(IB.EQ.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXCD
                BETA(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
C
C           CASE 3: IB = 0
            ELSEIF(IB.EQ.0) THEN
              DO M=1,MAXCD
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
                  DO M=1,MAXCD
                    BIN(M) = BIN(M)+X*TRM(M)
                    TRM(M) = TRM(M)*XJ(M,IBETA)
                  ENDDO
                ENDDO
                DO M=1,MAXCD
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)*BIN(M)
                ENDDO
              ELSEIF(IA.EQ.1) THEN
                DO M=1,MAXCD
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)
                ENDDO
              ELSE
                DO M=1,MAXCD
                  BETA(M) = BTA1(M)
                ENDDO
              ENDIF
C
C           END CONDITIONAL STATEMENT OVER IB VALUES
            ENDIF
C
            IF(IBETA.EQ.1) THEN
              DO M=1,MAXCD
                BDU(M, NUI+2*NX-1,-NUF+2*NY-2) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXCD
                BDL(M, NUI+2*NX-1,-NUF+2*NY-2) = BETA(M)
              ENDDO
            ENDIF

          ENDDO

        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ANGCLM1(DKAB,DKCD,KQN,LQN,ISEL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       AA    NN    NN  GGGGGG   CCCCCC  LL       MM       MM  11      C
C      AAAA   NNN   NN GG    GG CC    CC LL       MMM     MMM 111      C
C     AA  AA  NNNN  NN GG       CC       LL       MMMM   MMMM  11      C
C    AA    AA NN NN NN GG       CC       LL       MM MM MM MM  11      C
C    AAAAAAAA NN  NNNN GG   GGG CC       LL       MM  MMM  MM  11      C
C    AA    AA NN   NNN GG    GG CC    CC LL       MM   M   MM  11      C
C    AA    AA NN    NN  GGGGGG   CCCCCC  LLLLLLLL MM       MM 1111     C
C                                                                      C
C -------------------------------------------------------------------- C
C  ANGCLM1 PERFORMS THE ANGULAR ANALYSIS FOR THE EVALUATION OF TWO     C
C  ELECTRON INTEGRALS USING RACAH ALGEBRA TECHNIQUES. (OPEN-SHELL.)    C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION KQN(4),LQN(4),JQN(4)
      DIMENSION DKAB(MNU,MKP+1,MKP+1),DKCD(MNU,MKP+1,MKP+1)
C
      COMMON/XNUS/INU(MNU,16),NUS(MNU),NUI,NUF,NUNUM,K4AD
C
C**********************************************************************C
C     PARITY ANALYSIS: CHECK UNDERLYING LQN COMBINATIONS AND EXIT      C
C     IF THERE IS NO MULTIPOLE EXPANSION OF THE INTERACTION.           C
C**********************************************************************C
C
C     A AND B: PARITY OF 'LQN(1)+LQN(2)' (0 IF EVEN, 1 IF ODD)
      IPARAB = MOD(LQN(1)+LQN(2),2)
C
C     C AND D: PARITY OF 'LQN(3)+LQN(4)' (0 IF EVEN, 1 IF ODD)
      IPARCD = MOD(LQN(3)+LQN(4),2)
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
C     TENSOR LIMITS BY TRIANGLE RULE
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
C     FURTHER PARITY ANALYSIS: RUN OVER ALLOWED TENSOR ORDERS NU AND   C
C     SCREEN THOSE WHICH ARE NOT ALLOWED ON PARITY GROUNDS.            C
C**********************************************************************C
C
      NUNUM = 0
      DO NU=NUI,NUF
C
C       A AND B: PARITY OF 'LQN(1)+LQN(2)+NU' (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQN(1)+LQN(2)+NU,2)
C
C       C AND D: PARITY OF 'LQN(3)+LQN(4)+NU' (0 IF EVEN, 1 IF ODD)
        IPARCD = MOD(LQN(3)+LQN(4)+NU,2)
C
C       LQN SELECTION RULE: WHEN SUM OF A AND B, C AND D EVEN
        IF(IPARAB.EQ.0.AND.IPARCD.EQ.0) THEN
C
C         INCREASE TOTAL NUMBER OF NU VALUES TO BE FACILITATED
          NUNUM = NUNUM+1
C
C         SAVE THIS PARTICULAR NU VALUE TO A LIST
          NUS(NUNUM) = NU
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
C     WIGNER-ECKART: EVALUATE AN ANGULAR FACTOR FOR EVERY |MQN|        C
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
            DKAB(LTEN,MJA+1,MJB  ) = DK(JQN(1), MJA,JQN(2),-MJB,NU)
            DKAB(LTEN,MJA  ,MJB+1) = DK(JQN(1),-MJA,JQN(2), MJB,NU)
            DKAB(LTEN,MJA+1,MJB+1) = DK(JQN(1), MJA,JQN(2), MJB,NU)
C
          ENDDO
C
        ENDDO
      ENDDO
C
C     THE ARGUMENTS OF DK COEFFICIENTS ARE REVERSED IN THE CASE OF THE
C     CD PAIRS IN ORDER TO ACCOMMODATE THE RELATION:
C               DK(J,M,J',M',L) = ((-1)^Q)*DK(J',M',J,M,L).
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
            DKCD(LTEN,MJC+1,MJD  ) = DK(JQN(4),-MJD,JQN(3), MJC,NU)
            DKCD(LTEN,MJC  ,MJD+1) = DK(JQN(4), MJD,JQN(3),-MJC,NU)
            DKCD(LTEN,MJC+1,MJD+1) = DK(JQN(4), MJD,JQN(3), MJC,NU)
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
      SUBROUTINE KLSET1(IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            KK    KK LL       SSSSSS  EEEEEEEE TTTTTTTT 11            C
C            KK   KK  LL      SS    SS EE          TT   111            C
C            KK  KK   LL      SS       EE          TT    11            C
C            KKKKK    LL       SSSSSS  EEEEEE      TT    11            C
C            KK  KK   LL            SS EE          TT    11            C
C            KK   KK  LL      SS    SS EE          TT    11            C
C            KK    KK LLLLLLLL SSSSSS  EEEEEEEE    TT   1111           C
C                                                                      C
C -------------------------------------------------------------------- C
C  KLSET1 PREPARES BASIS SET INTERMEDIATES FOR ALL IJ-PAIRS OF BASIS   C
C  FUNCTIONS, FOR USE IN ROUTINES RKCLM1/RKBRT1.                       C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RNLC(MBS),RNSC(MBS),RNLD(MBS),RNSD(MBS)
      DIMENSION EKL0(MB2)
C
C     COMMON/B0IJ/EIJ(MB2,-MAB:MAB),RNIJ(MB2,4),EI(MB2),EJ(MB2),MAXAB
      COMMON/B0KL/EKL(MB2,-MAB:MAB),RNKL(MB2,4),EK(MB2),EL(MB2),MAXCD
      COMMON/B1QN/EXL(MBS,4),MQN(4),KQN(4),LQN(4),NBAS(4),IBAS,JBAS
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/XNUS/INU(MNU,16),NUS(MNU),NUI,NUF,NUNUM,K4AD
C
C     NORMALISATION CONSTANTS FOR EXPONENTS IN LIST C
      RL = DFLOAT(LQN(3))
      G1 = TWLG-GAMLOG(2*LQN(3)+3)
      G2 = TWLG-GAMLOG(2*LQN(3)+5)
      R1 = RL+1.5D0
      R2 = RL+0.5D0
      DO KBAS=1,NBAS(3)
        EKV        = EXL(KBAS,3)
        ELOG       = DLOG(2.0D0*EKV)
        RNLC(KBAS) = DEXP(0.5D0*(G1+R1*ELOG))
        RNSC(KBAS) = DEXP(0.5D0*(G2+R2*ELOG))
      ENDDO
C
C     NORMALISATION CONSTANTS FOR EXPONENTS IN LIST D
      RL = DFLOAT(LQN(4))
      G1 = TWLG-GAMLOG(2*LQN(4)+3)
      G2 = TWLG-GAMLOG(2*LQN(4)+5)
      R1 = RL+1.5D0
      R2 = RL+0.5D0
      DO LBAS=1,NBAS(4)
        ELV        = EXL(LBAS,4)
        ELOG       = DLOG(2.0D0*ELV)
        RNLD(LBAS) = DEXP(0.5D0*(G1+R1*ELOG))
        RNSD(LBAS) = DEXP(0.5D0*(G2+R2*ELOG))
      ENDDO
C
C     LIST OF EXPONENTS AND NORMALISATION COEFFICIENTS IN THE BLOCK
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
          EK(M)     = EXL(KBAS,3)
          EL(M)     = EXL(LBAS,4)
          EKL0(M)   = EK(M)+EL(M)
          RNKL(M,1) = RNLC(KBAS)*RNLD(LBAS)
          RNKL(M,2) = RNLC(KBAS)*RNSD(LBAS)
          RNKL(M,3) = RNSC(KBAS)*RNLD(LBAS)
          RNKL(M,4) = RNSC(KBAS)*RNSD(LBAS)
        ENDDO
      ENDDO
C
C     TENSOR ORDER LIMITS
      NUI = NUS(1)
      NUF = NUS(NUNUM)
C
C     LOWEST EXPONENT POWER
      IPOWER = LQN(3)+LQN(4)-NUF
C
C     GAUSSIAN OVERLAPS TO ALL REQUIRED POWERS
      DO M=1,MAXCD
C
C       SEED GAUSSIAN EXPONENT VALUE
        EKLR = DSQRT(EKL0(M))
        EKL(M,-NUF) = EKLR**(-IPOWER)
C
C       LOOP OVER ALL TENSOR ORDERS AND DIVIDE
        DO IPOW=-NUF+1,NUF+5
          EKL(M,IPOW) = EKL(M,IPOW-1)/EKLR
        ENDDO
C
      ENDDO
C
      RETURN
      END

