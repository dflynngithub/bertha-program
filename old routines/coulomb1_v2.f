      SUBROUTINE COULOMB1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  CCCCCC   OOOOOO  UU    UU LL      OOOOOO  MM       MM BBBBBBB   11  C
C CC    CC OO    OO UU    UU LL     OO    OO MMM     MMM BB    BB 111  C
C CC       OO    OO UU    UU LL     OO    OO MMMM   MMMM BB    BB  11  C
C CC       OO    OO UU    UU LL     OO    OO MM MM MM MM BBBBBBB   11  C
C CC       OO    OO UU    UU LL     OO    OO MM  MMM  MM BB    BB  11  C
C CC    CC OO    OO UU    UU LL     OO    OO MM   M   MM BB    BB  11  C
C  CCCCCC   OOOOOO   UUUUUU  LLLLLLL OOOOOO  MM       MM BBBBBBB  1111 C
C                                                                      C
C -------------------------------------------------------------------- C
C  COULOMB1 CONSTRUCTS ALL ONE-CENTRE CONTRIBUTIONS TO THE MOLECULAR   C
C  MEAN-FIELD COULOMB MATRIX BY USING A RACAH ALGEBRA DECOMPOSITION    C
C  OF TWO-ELECTRON MATRIX ELEMENTS.                                    C
C -------------------------------------------------------------------- C
C  (THIS ROUTINE IS SIMILAR TO THE IN-LINE CONSTRUCTION OF MEAN-FIELD  C
C   ATOMIC COULOMB MATRIX IN HFSCF0, BUT WITH MQN STRUCTURE AS WELL.)  C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=36,MB2=MBS*MBS,MCT=15,MKP=13,MNU=MKP+1,
     &                                                      MAB=2*MNU+6)
C
      CHARACTER*4 HMLTN
C
      DIMENSION EXPT(MBS,4),KQN(4),MQN(4),NBAS(4),LQN(4)
      DIMENSION DKAB(MNU,MKP+1,MKP+1),DKCD(MNU,MKP+1,MKP+1)
      DIMENSION RLL(MB2),RSL(MB2),RLS(MB2),RSS(MB2)
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
      COMMON/XBIJ/EIJ1(MAB),EIJ2(MAB),RNIJ(4),EI,EJ
      COMMON/XCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/XMAT/XLLLL(MB2,MNU),XLLSS(MB2,MNU),
     &            XSSLL(MB2,MNU),XSSSS(MB2,MNU)
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
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
C     PREPARE INTERMEDIATE DATA FOR USE IN RABCD                       C
C**********************************************************************C
C
C     COEFFICIENTS FOR INTEGRAL ASSEMBLY
      C1 = 0.25D0*GAMHLF(LQN(1)+LQN(2)+LQN(3)+LQN(4)+1)
      C3 = 0.25D0*GAMHLF(LQN(1)+LQN(2)+LQN(3)+LQN(4)+3)
      C5 = 0.25D0*GAMHLF(LQN(1)+LQN(2)+LQN(3)+LQN(4)+5)
      C7 = 0.25D0*GAMHLF(LQN(1)+LQN(2)+LQN(3)+LQN(4)+7)
      C9 = 0.25D0*GAMHLF(LQN(1)+LQN(2)+LQN(3)+LQN(4)+9)
C
      V1 = 1.0D0
      V2 = 2.0D0
      V4 = 4.0D0
      V8 = 8.0D0
      VS = 1.6D1
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
      CALL ANGCLM1(DKAB,DKCD,KQN,LQN,NUS,NUF,NUI,NUNUM,ISEL)
      CALL CPU_TIME(T2)
      TC1B = TC1B+T2-T1
C
C     EXIT THIS COMBINATION IF NO CONTRIBUTING MULTIPOLES IN ANGLES
      IF(ISEL.EQ.0) GOTO 2001
C
C     BASIS SET INTERMEDIATES FOR THE KL-PAIRS
      CALL CPU_TIME(T1)
      CALL PREPKL(LQN,KC,NBAS(3),KD,NBAS(4),ICNT)
      CALL CPU_TIME(T2)
      TC1B = TC1B+T2-T1
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS IN BLOCKS A AND B (INDEX 3000)         C
C**********************************************************************C
C
C     UPDATE COUNTER FOR NUMBER OF CLASSES
      IF(HMLTN.EQ.'NORL') THEN
        N2EB(1,1) = N2EB(1,1) + 1
      ELSE
        N2EB(1,1) = N2EB(1,1) + 1
        N2EB(1,2) = N2EB(1,2) + 1
        N2EB(1,3) = N2EB(1,3) + 1
        N2EB(1,4) = N2EB(1,4) + 1
      ENDIF
C
      DO 3000 IBAS=1,NBAS(1)
      DO 3000 JBAS=1,NBAS(2)
C
C       UPDATE COUNTER FOR NUMBER OF INTEGRALS
        IF(HMLTN.EQ.'NORL') THEN
          N2EI(1,1) = N2EI(1,1)+NBAS(3)*NBAS(4)
        ELSE
          N2EI(1,1) = N2EI(1,1)+NBAS(3)*NBAS(4)
          N2EI(1,2) = N2EI(1,2)+NBAS(3)*NBAS(4)
          N2EI(1,3) = N2EI(1,3)+NBAS(3)*NBAS(4)
          N2EI(1,4) = N2EI(1,4)+NBAS(3)*NBAS(4)
        ENDIF
C
C       GAUSSIAN EXPONENTS FOR THIS PAIR
        EI = EXPT(IBAS,1)
        EJ = EXPT(JBAS,2)
C
C       BASIS SET INTERMEDIATES FOR THE IJ-PAIRS
        CALL CPU_TIME(T1)
        CALL PREPIJ(LQN)
        CALL CPU_TIME(T2)
        TC1B = TC1B+T2-T1
C
C       BATCH OF RADIAL INTEGRALS (EFFECTIVE INTERACTION STRENGTHS)
        CALL CPU_TIME(T1)
        CALL RABCD(KQN,LQN,NBAS(3)*NBAS(4))
        CALL CPU_TIME(T2)
        TC1R = TC1R+T2-T1
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
        RLL(M) = 0.0D0
        RLS(M) = 0.0D0
        RSL(M) = 0.0D0
        RSS(M) = 0.0D0
      ENDDO
C
C     ASSEMBLE INTEGRAL BY SUMMING OVER EXPANSION POWER L
      DO LTEN=1,NUNUM
C
C       ANGULAR COEFFICIENT
        ANGFAC = DKAB(LTEN,IMJA,IMJB)*DKCD(LTEN,IMJC,IMJD)
C
C       SCREENING OF INTEGRAL BASED ON ANGULAR COEFFICIENT
        IF(DABS(ANGFAC).LE.SENS) GOTO 5003
C
        IF(HMLTN.EQ.'NORL') THEN
C       NON-RELATIVISTIC HAMILTONIAN
C
          DO M=1,NBAS(3)*NBAS(4)
            RLL(M) = RLL(M) + ANGFAC*XLLLL(M,LTEN)
          ENDDO
C
        ELSE
C       RELATIVISTIC HAMILTONIAN
C
          DO M=1,NBAS(3)*NBAS(4)
            RLL(M) = RLL(M) + ANGFAC*XLLLL(M,LTEN)
            RLS(M) = RLS(M) + ANGFAC*XLLSS(M,LTEN)
            RSL(M) = RSL(M) + ANGFAC*XSSLL(M,LTEN)
            RSS(M) = RSS(M) + ANGFAC*XSSSS(M,LTEN)
          ENDDO
C
        ENDIF
C
C       SKIP POINT FOR ANGULAR SCREENING
5003    CONTINUE
C
      ENDDO
C
C     FULL-INTEGRAL CONSTRUCTION COMPLETE
      CALL CPU_TIME(T2)
      TC1F = TC1F+T2-T1
C
C**********************************************************************C
C     ADD THIS BATCH OF R-INTEGRALS TO CLOSED-SHELL COULOMB MATRIX     C
C**********************************************************************C
C
      IF(HMLTN.EQ.'NORL') THEN
C     NON-RELATIVISTIC HAMILTONIAN
C
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
C           DIRECT CONTRIBUTIONS
            GDIR(NAL+IBAS,NBL+JBAS) = GDIR(NAL+IBAS,NBL+JBAS)
     &                    +              RLL(M)*DENC(NCL+KBAS,NDL+LBAS)
C
C           EXCHANGE CONTRIBUTIONS
            GXCH(NAL+IBAS,NDL+LBAS) = GXCH(NAL+IBAS,NDL+LBAS)
     &                    +              RLL(M)*DENC(NCL+KBAS,NBL+JBAS)
C
          ENDDO
        ENDDO
C
      ELSE
C     RELATIVISTIC HAMILTONIAN
C
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
C           DIRECT CONTRIBUTIONS
            GDIR(NAL+IBAS,NBL+JBAS) = GDIR(NAL+IBAS,NBL+JBAS)
     &                    +              RLL(M)*DENC(NCL+KBAS,NDL+LBAS)
     &                    +              RLS(M)*DENC(NCS+KBAS,NDS+LBAS)
C
            GDIR(NAS+IBAS,NBS+JBAS) = GDIR(NAS+IBAS,NBS+JBAS)
     &                    +              RSL(M)*DENC(NCL+KBAS,NDL+LBAS)
     &                    +              RSS(M)*DENC(NCS+KBAS,NDS+LBAS)
C
C           EXCHANGE CONTRIBUTIONS
            GXCH(NAL+IBAS,NDL+LBAS) = GXCH(NAL+IBAS,NDL+LBAS)
     &                    +              RLL(M)*DENC(NCL+KBAS,NBL+JBAS)
C
            GXCH(NAL+IBAS,NDS+LBAS) = GXCH(NAL+IBAS,NDS+LBAS)
     &                    +              RLS(M)*DENC(NCS+KBAS,NBL+JBAS)
C
            GXCH(NAS+IBAS,NDL+LBAS) = GXCH(NAS+IBAS,NDL+LBAS)
     &                    +              RSL(M)*DENC(NCL+KBAS,NBS+JBAS)

            GXCH(NAS+IBAS,NDS+LBAS) = GXCH(NAS+IBAS,NDS+LBAS)
     &                    +              RSS(M)*DENC(NCS+KBAS,NBS+JBAS)
C
          ENDDO
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C     ADD THIS BATCH OF R-INTEGRALS TO OPEN-SHELL COULOMB MATRIX       C
C**********************************************************************C
C
      IF(NOPN.EQ.0) GOTO 6100
C
      IF(HMLTN.EQ.'NORL') THEN
C     NON-RELATIVISTIC HAMILTONIAN
C
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
C           DIRECT CONTRIBUTIONS
            QDIR(NAL+IBAS,NBL+JBAS) = QDIR(NAL+IBAS,NBL+JBAS)
     &                    +         ACFF*RLL(M)*DENO(NCL+KBAS,NDL+LBAS)
C
C           EXCHANGE CONTRIBUTIONS
            QXCH(NAL+IBAS,NDL+LBAS) = QXCH(NAL+IBAS,NDL+LBAS)
     &                    +         BCFF*RLL(M)*DENO(NCL+KBAS,NBL+JBAS)
C
          ENDDO
        ENDDO
C
      ELSE
C     RELATIVISTIC HAMILTONIAN
C
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
C           DIRECT CONTRIBUTIONS
            QDIR(NAL+IBAS,NBL+JBAS) = QDIR(NAL+IBAS,NBL+JBAS)
     &                    +         ACFF*RLL(M)*DENO(NCL+KBAS,NDL+LBAS)
     &                    +         ACFF*RLS(M)*DENO(NCS+KBAS,NDS+LBAS)
C
            QDIR(NAS+IBAS,NBS+JBAS) = QDIR(NAS+IBAS,NBS+JBAS)
     &                    +         ACFF*RSL(M)*DENO(NCL+KBAS,NDL+LBAS)
     &                    +         ACFF*RSS(M)*DENO(NCS+KBAS,NDS+LBAS)
C
C           EXCHANGE CONTRIBUTIONS
            QXCH(NAL+IBAS,NDL+LBAS) = QXCH(NAL+IBAS,NDL+LBAS)
     &                    +         BCFF*RLL(M)*DENO(NCL+KBAS,NBL+JBAS)
C
            QXCH(NAL+IBAS,NDS+LBAS) = QXCH(NAL+IBAS,NDS+LBAS)
     &                    +         BCFF*RLS(M)*DENO(NCS+KBAS,NBL+JBAS)
C
            QXCH(NAS+IBAS,NDL+LBAS) = QXCH(NAS+IBAS,NDL+LBAS)
     &                    +         BCFF*RSL(M)*DENO(NCL+KBAS,NBS+JBAS)

            QXCH(NAS+IBAS,NDS+LBAS) = QXCH(NAS+IBAS,NDS+LBAS)
     &                    +         BCFF*RSS(M)*DENO(NCS+KBAS,NBS+JBAS)
C
          ENDDO
        ENDDO
C
      ENDIF
C
6100  CONTINUE
C
C     MATRIX MULTIPLICATION STEP COMPLETE
      CALL CPU_TIME(T3)
      TC1M = TC1M+T3-T2
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
C     RECORD CPU TIME AT END OF BATCH AND ADD TO APPROPRIATE COUNTER
      CALL CPU_TIME(TBCH2)
      IF(HMLTN.EQ.'NORL') THEN
        T2ES(1,1) = T2ES(1,1) + TBCH2 - TBCH1
      ELSE
        T2ES(1,1) = T2ES(1,1) + 0.25D0*(TBCH2-TBCH1)
        T2ES(1,2) = T2ES(1,2) + 0.25D0*(TBCH2-TBCH1)
        T2ES(1,3) = T2ES(1,3) + 0.25D0*(TBCH2-TBCH1)
        T2ES(1,4) = T2ES(1,4) + 0.25D0*(TBCH2-TBCH1)
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE RABCD(KQN,LQN,MAXM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             RRRRRRR     AA    BBBBBBB   CCCCCC  DDDDDDD              C
C             RR    RR   AAAA   BB    BB CC    CC DD    DD             C
C             RR    RR  AA  AA  BB    BB CC       DD    DD             C
C             RR    RR AA    AA BBBBBBB  CC       DD    DD             C
C             RRRRRRR  AAAAAAAA BB    BB CC       DD    DD             C
C             RR    RR AA    AA BB    BB CC    CC DD    DD             C
C             RR    RR AA    AA BBBBBBB   CCCCCC  DDDDDDD              C
C                                                                      C
C -------------------------------------------------------------------- C
C  RABCD EVALUATES THE EFFECTIVE INTERACTION STRENGTHS OF THE COULOMB  C
C  AND BREIT INTERACTIONS IN A G-SPINOR BASIS SET. THE ALGORITHM IS    C
C  SEGMENTED TO TAKE INTO ACCOUNT THE SIGNS OF KQN IN THE T=S CASE.    C
C**********************************************************************C
      PARAMETER(MBS=36,MB2=MBS*MBS,MKP=13,MNU=MKP+1,MAB=2*MNU+6)
C
      CHARACTER*4 HMLTN
C
      DIMENSION XJ(MB2,2),IAA(2),IBB(2)
      DIMENSION BETA(MB2),BTA1(MB2),XROOT(MB2),BIN(MB2),TRM(MB2)
      DIMENSION BDU(MB2,MAB,MAB),BDL(MB2,MAB,MAB)
C
      DIMENSION KQN(4),LQN(4)
C
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/XBIJ/EIJ1(MAB),EIJ2(MAB),RNIJ(4),EI,EJ
      COMMON/XBKL/EK(MB2),EKL0(MB2),
     &            EL(MB2),EKL1(MB2,MAB),EKL2(MB2,MAB),RNKL(MB2,4)
      COMMON/XCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/XMAT/XLLLL(MB2,MNU),XLLSS(MB2,MNU),
     &            XSSLL(MB2,MNU),XSSSS(MB2,MNU)
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
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
                BDU(M,NX,NY) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXM
                BDL(M,NX,NY) = BETA(M)
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
          XLLLL(M,LTEN) = 0.0D0
          XLLSS(M,LTEN) = 0.0D0
          XSSLL(M,LTEN) = 0.0D0
          XSSSS(M,LTEN) = 0.0D0
        ENDDO
      ENDDO
C
C     VALUE PREPARATION
      E0000 = 1.0D0
      E1000 = EI
      E0100 = EJ
      E1100 = EI*EJ
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
        E0011 = EK(M)*EL(M)
        E1110 = EI*EJ*EK(M)
        E1101 = EI*EJ*EL(M)
        E1011 = EI*EK(M)*EL(M)
        E0111 = EJ*EK(M)*EL(M)
        E1111 = EI*EJ*EK(M)*EL(M)
C
C       LOOP OVER THE TENSOR ORDERS OF THE COULOMB INTERACTION
        DO LTEN=1,NUNUM
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(LTEN)
C
C         INDEX OFFSETS FOR TENSOR ORDER
          IQ1 = (-NUI+NU)/2
          IQ2 = ( NUF-NU)/2
C
C**********************************************************************C
C     SPECIAL CASE: NON-RELATIVISTIC HAMILTONIAN                       C
C**********************************************************************C
C
          IF(HMLTN.EQ.'NORL') THEN
C         NON-RELATIVISTIC HAMILTONIAN
C
C           BETA INTEGRALS WITH POWERS OF EXPONENTIALS
            B22 = EIJ2(IQ2+2)*EKL2(M,IQ1+2)*BDL(M,IQ1+2,IQ2+2)
     &          + EIJ1(IQ1+2)*EKL1(M,IQ2+2)*BDU(M,IQ1+2,IQ2+2)
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XLLLL(M,LTEN) = C5*B22
C
C           SKIP PAST THE RELATIVISTIC STAGES
            GOTO 999
C
          ENDIF
C
C**********************************************************************C
C     LARGER CASE: RELATIVISTIC HAMILTONIAN (T TERMS GIVE ZEROES)      C
C**********************************************************************C
C
C         TEMPORARY STORAGE OF RAW RJ(LTEN,M)
          B11 = EIJ2(IQ2+1)*EKL2(M,IQ1+1)*BDL(M,IQ1+1,IQ2+1)
     &        + EIJ1(IQ1+1)*EKL1(M,IQ2+1)*BDU(M,IQ1+1,IQ2+1)
          B12 = EIJ2(IQ2+2)*EKL2(M,IQ1+1)*BDL(M,IQ1+1,IQ2+2)
     &        + EIJ1(IQ1+2)*EKL1(M,IQ2+1)*BDU(M,IQ1+2,IQ2+1)
          B13 = EIJ2(IQ2+3)*EKL2(M,IQ1+1)*BDL(M,IQ1+1,IQ2+3)
     &        + EIJ1(IQ1+3)*EKL1(M,IQ2+1)*BDU(M,IQ1+3,IQ2+1)
          B21 = EIJ2(IQ2+1)*EKL2(M,IQ1+2)*BDL(M,IQ1+2,IQ2+1)
     &        + EIJ1(IQ1+1)*EKL1(M,IQ2+2)*BDU(M,IQ1+1,IQ2+2)
          B22 = EIJ2(IQ2+2)*EKL2(M,IQ1+2)*BDL(M,IQ1+2,IQ2+2)
     &        + EIJ1(IQ1+2)*EKL1(M,IQ2+2)*BDU(M,IQ1+2,IQ2+2)
          B23 = EIJ2(IQ2+3)*EKL2(M,IQ1+2)*BDL(M,IQ1+2,IQ2+3)
     &        + EIJ1(IQ1+3)*EKL1(M,IQ2+2)*BDU(M,IQ1+3,IQ2+2)
          B31 = EIJ2(IQ2+1)*EKL2(M,IQ1+3)*BDL(M,IQ1+3,IQ2+1)
     &        + EIJ1(IQ1+1)*EKL1(M,IQ2+3)*BDU(M,IQ1+1,IQ2+3)
          B32 = EIJ2(IQ2+2)*EKL2(M,IQ1+3)*BDL(M,IQ1+3,IQ2+2)
     &        + EIJ1(IQ1+2)*EKL1(M,IQ2+3)*BDU(M,IQ1+2,IQ2+3)
          B33 = EIJ2(IQ2+3)*EKL2(M,IQ1+3)*BDL(M,IQ1+3,IQ2+3)
     &        + EIJ1(IQ1+3)*EKL1(M,IQ2+3)*BDU(M,IQ1+3,IQ2+3)
C
C         EFFECTIVE INTERACTION STRENGTH XLLLL(M,LTEN)
          XLLLL(M,LTEN) = V1*T0000*E0000*C5*B22
C
C         EFFECTIVE INTERACTION STRENGTH XLLSS(M,LTEN)
          XLLSS(M,LTEN) = V4*T0000*E0011*C7*B32 - V2*T0001*E0010*C5*B22
     &                  - V2*T0010*E0001*C5*B22 + V1*T0011*E0000*C3*B12
C
C         EFFECTIVE INTERACTION STRENGTH XSSLL(M,LTEN)
          XSSLL(M,LTEN) = V4*T0000*E1100*C7*B23 - V2*T0100*E1000*C5*B22
     &                  - V2*T1000*E0100*C5*B22 + V1*T1100*E0000*C3*B21
C
C         EFFECTIVE INTERACTION STRENGTH XSSSS(M,LTEN)
          XSSSS(M,LTEN) = VS*T0000*E1111*C9*B33
     &                  - V8*T0001*E1110*C7*B23 - V8*T0010*E1101*C7*B23
     &                  - V8*T0100*E1011*C7*B32 - V8*T1000*E0111*C7*B32
     &                  + V4*T1100*E0011*C5*B31 + V4*T0011*E1100*C5*B13
     &                  + V4*T1001*E0110*C5*B22 + V4*T0110*E1001*C5*B22
     &                  + V4*T0101*E1010*C5*B22 + V4*T1010*E0101*C5*B22
     &                  - V2*T1101*E0010*C3*B21 - V2*T0111*E1000*C3*B12
     &                  - V2*T1110*E0001*C3*B21 - V2*T1011*E0100*C3*B12
     &                  + V1*T1111*E0000*C1*B11
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
      IF(HMLTN.EQ.'NORL') THEN
C     NON-RELATIVISTIC HAMILTONIAN
C
        DO M=1,MAXM
          RNLLLL = RNIJ(1)*RNKL(M,1)
          DO LTEN=1,NUNUM
            XLLLL(M,LTEN) = RNLLLL*XLLLL(M,LTEN)
          ENDDO
        ENDDO
C
      ELSE
C     RELATIVISTIC HAMILTONIAN
C
        DO M=1,MAXM
          RNLLLL = RNIJ(1)*RNKL(M,1)
          RNLLSS = RNIJ(1)*RNKL(M,3)
          RNSSLL = RNIJ(3)*RNKL(M,1)
          RNSSSS = RNIJ(3)*RNKL(M,3)
          DO LTEN=1,NUNUM
            XLLLL(M,LTEN) = RNLLLL*XLLLL(M,LTEN)
            XLLSS(M,LTEN) = RNLLSS*XLLSS(M,LTEN)
            XSSLL(M,LTEN) = RNSSLL*XSSLL(M,LTEN)
            XSSSS(M,LTEN) = RNSSSS*XSSSS(M,LTEN)
          ENDDO
        ENDDO
C
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE PREPIJ(LQN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          PPPPPPP  RRRRRRR  EEEEEEEE PPPPPPP  IIII    JJJJJ           C
C          PP    PP RR    RR EE       PP    PP  II       JJ            C
C          PP    PP RR    RR EE       PP    PP  II       JJ            C
C          PP    PP RR    RR EEEEEE   PP    PP  II       JJ            C
C          PPPPPPP  RRRRRRR  EE       PPPPPPP   II       JJ            C
C          PP       RR    RR EE       PP        II JJ    JJ            C
C          PP       RR    RR EEEEEEEE PP       IIII JJJJJJ             C
C                                                                      C
C -------------------------------------------------------------------- C
C  PREPIJ INITIALISES VARIABLES REQUIRED BY SUBROUTINE RABCD.          C
C  STRUCTURE SIMILAR TO RNORM0 ROUTINE, EXCEPT FOR ONLY ONE PAIR.      C
C**********************************************************************C
      PARAMETER(MBS=36,MB2=MBS*MBS,MKP=13,MNU=MKP+1,MAB=2*MNU+6)
C
      DIMENSION LQN(4)
C
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
      COMMON/XBIJ/EIJ1(MAB),EIJ2(MAB),RNIJ(4),EI,EJ
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
C
      DATA TWOLOG/6.93147180559945309D-1/
C
C     NORMALISATION CONSTANTS FOR EXPONENT EI
      RL = DFLOAT(LQN(1))
      G1 = TWOLOG-GAMLOG(2*LQN(1)+3)
      G2 = TWOLOG-GAMLOG(2*LQN(1)+5)
      R1 = RL+1.5D0
      R2 = RL+0.5D0
C
      ELOG = DLOG(2.0D0*EI)
      RNLI = DEXP(0.5D0*(G1+R1*ELOG))
      RNSI = DEXP(0.5D0*(G2+R2*ELOG))
C
C     NORMALISATION CONSTANTS FOR EXPONENT EJ
      RL = DFLOAT(LQN(2))
      G1 = TWOLOG-GAMLOG(2*LQN(2)+3)
      G2 = TWOLOG-GAMLOG(2*LQN(2)+5)
      R1 = RL+1.5D0
      R2 = RL+0.5D0
C
      ELOG = DLOG(2.0D0*EJ)
      RNLJ = DEXP(0.5D0*(G1+R1*ELOG))
      RNSJ = DEXP(0.5D0*(G2+R2*ELOG))
C
C     COMPOSITE NORMALISATION CONSTANTS
      RNIJ(1) = RNLI*RNLJ
      RNIJ(2) = RNSI*RNLJ
      RNIJ(3) = RNSI*RNSJ
      RNIJ(4) = RNLI*RNSJ
C
C     NUMBER OF EIJ POWERS REQUIRED FOR COULOMB SUM
      NOVALS = (NUS(NUNUM)-NUS(1))/2 + 3
C
      IPOWER1 = LQN(1)+LQN(2)+NUS(1)+1
      IPOWER2 = LQN(1)+LQN(2)-NUS(NUNUM)
C
C     SEED EXPONENT POWER
      EIJ0 = EI+EJ
      EIJR = DSQRT(EIJ0)
C
C     STARTING POINT FOR EIJ1 AND EIJ2
      EIJ1A = EIJR**(-IPOWER1)
      EIJ2A = EIJR**(-IPOWER2)
C
C     SAVE VALUES TO ARRAYS UP TO A FINAL VLAUE
      DO IPOW=1,NOVALS
        EIJ1(IPOW) = EIJ1A
        EIJ2(IPOW) = EIJ2A
        EIJ1A = EIJ1A/EIJ0
        EIJ2A = EIJ2A/EIJ0
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE PREPKL(LQN,KC,NBASC,KD,NBASD,ICNT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        PPPPPPP  RRRRRRR  EEEEEEEE PPPPPPP  KK    KK LL               C
C        PP    PP RR    RR EE       PP    PP KK   KK  LL               C
C        PP    PP RR    RR EE       PP    PP KK  KK   LL               C
C        PP    PP RR    RR EEEEEE   PP    PP KKKKK    LL               C
C        PPPPPPP  RRRRRRR  EE       PPPPPPP  KK  KK   LL               C
C        PP       RR    RR EE       PP       KK   KK  LL               C
C        PP       RR    RR EEEEEEEE PP       KK    KK LLLLLLLL         C
C                                                                      C
C -------------------------------------------------------------------- C
C  PREPKL INITIALISES VARIABLES REQUIRED BY SUBROUTINE RABCD.          C
C  STRUCTURE SIMILAR TO RNORM0 ROUTINE, EXCEPT OVERLAPS ARE DIFFERENT. C
C**********************************************************************C
      PARAMETER(MCT=15,MBS=36,MB2=MBS*MBS,MKP=13,MNU=MKP+1,MAB=2*MNU+6)
C
      DIMENSION RNLC(MBS),RNSC(MBS),RNLD(MBS),RNSD(MBS)
      DIMENSION LQN(4)
C
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
      COMMON/XBKL/EK(MB2),EKL0(MB2),
     &            EL(MB2),EKL1(MB2,MAB),EKL2(MB2,MAB),RNKL(MB2,4)
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
C
      DATA TWOLOG/6.93147180559945309D-1/
C
C     BLOCK LIST LENGTH
      MAXM = NBASC*NBASD
C
C     NORMALISATION CONSTANTS FOR EXPONENTS IN LIST C
      RLQN = DFLOAT(LQN(3))
      G1   = TWOLOG - GAMLOG(2*LQN(3)+3)
      G2   = TWOLOG - GAMLOG(2*LQN(3)+5)
      R1   = RLQN + 1.5D0
      R2   = RLQN + 0.5D0
      DO KBAS=1,NBASC
        EKV        = EXPSET(KBAS,LQN(3)+1,ICNT)
        ELOG       = DLOG(2.0D0*EKV)
        RNLC(KBAS) = DEXP(0.5D0*(G1 + R1*ELOG))
        RNSC(KBAS) = DEXP(0.5D0*(G2 + R2*ELOG))
      ENDDO
C
C     NORMALISATION CONSTANTS FOR EXPONENTS IN LIST D
      RLQN = DFLOAT(LQN(4))
      G1   = TWOLOG - GAMLOG(2*LQN(4)+3)
      G2   = TWOLOG - GAMLOG(2*LQN(4)+5)
      R1   = RLQN + 1.5D0
      R2   = RLQN + 0.5D0
      DO LBAS=1,NBASD
        ELV        = EXPSET(LBAS,LQN(4)+1,ICNT)
        ELOG       = DLOG(2.0D0*ELV)
        RNLD(LBAS) = DEXP(0.5D0*(G1 + R1*ELOG))
        RNSD(LBAS) = DEXP(0.5D0*(G2 + R2*ELOG))
      ENDDO
C
C     LIST OF EXPONENTS AND NORMALISATION COEFFICIENTS IN THE BLOCK
      M = 0
      DO KBAS=1,NBASC
        DO LBAS=1,NBASD
          M = M+1
          EK(M)     = EXPSET(KBAS,LQN(3)+1,ICNT)
          EL(M)     = EXPSET(LBAS,LQN(4)+1,ICNT)
          EKL0(M)   = EK(M)+EL(M)
          RNKL(M,1) = RNLC(KBAS)*RNLD(LBAS)
          RNKL(M,2) = RNSC(KBAS)*RNLD(LBAS)
          RNKL(M,3) = RNSC(KBAS)*RNSD(LBAS)
          RNKL(M,4) = RNLC(KBAS)*RNSD(LBAS)
        ENDDO
      ENDDO
C
C     NUMBER OF EKL POWERS REQUIRED FOR COULOMB SUM
      NOVALS  = ((NUS(NUNUM)-NUS(1))/2) + 3
C
      IPOWER1 = LQN(3)+LQN(4)-NUS(NUNUM)
      IPOWER2 = LQN(3)+LQN(4)+NUS(    1)+1
C
C     GAUSSIAN OVERLAPS TO ALL REQUIRED POWERS
      DO M=1,MAXM
C
C       SEED GAUSSIAN EXPONENT VALUE
        EKLR = DSQRT(EKL0(M))
C
C       STARTING POINT FOR EKL1 AND EKL2
        EKL1A = EKLR**(-IPOWER1)
        EKL2A = EKLR**(-IPOWER2)
        DO IPOWER=1,NOVALS
          EKL1(M,IPOWER) = EKL1A
          EKL2(M,IPOWER) = EKL2A
          EKL1A = EKL1A/EKL0(M)
          EKL2A = EKL2A/EKL0(M)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ANGCLM1(DKAB,DKCD,KQN,LQN,NUS,NUF,NUI,NUNUM,ISEL)
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
      PARAMETER(MKP=13,MNU=MKP+1)
C
      DIMENSION KQN(4),LQN(4),JQN(4),NUS(MNU)
      DIMENSION DKAB(MNU,MKP+1,MKP+1),DKCD(MNU,MKP+1,MKP+1)
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
C     FURTHER PARITY ANALYSIS: RUN OVER ALLOWED TENSOR ORDERS NU AND   C
C     SCREEN THOSE WHICH ARE NOT ALLOWED ON PARITY GROUNDS.            C
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
C       LQN SELECTION RULE: WHEN SUM OF A AND B, C AND D EVEN
        IF(IPARAB.EQ.1.AND.IPARCD.EQ.1) THEN
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
      NUI = NUS(    1)
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
C           READ THE ACTUAL ODER NU
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
C           READ THE ACTUAL ODER NU
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
