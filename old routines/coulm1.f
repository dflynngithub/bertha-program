      SUBROUTINE COULM1(IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          CCCCCC   OOOOOO  UU    UU LL       MM       MM  11          C
C         CC    CC OO    OO UU    UU LL       MMM     MMM 111          C
C         CC       OO    OO UU    UU LL       MMMM   MMMM  11          C
C         CC       OO    OO UU    UU LL       MM MM MM MM  11          C
C         CC       OO    OO UU    UU LL       MM  MMM  MM  11          C
C         CC    CC OO    OO UU    UU LL       MM   M   MM  11          C
C          CCCCCC   OOOOOO   UUUUUU  LLLLLLLL MM       MM 1111         C
C                                                                      C
C -------------------------------------------------------------------- C
C  COULM1 CONSTRUCT A ONE-CENTRE CONTRIBUTION TO THE MOLECULAR         C
C  MEAN-FIELD COULOMB MATRIX WITH RACAH ALGEBRA AND BETA INTEGRALS.    C
C  THIS ROUTINE IS SIMILAR TO THE IN-LINE CONSTRUCTION OF MEAN-FIELD   C
C  ATOMIC COULOMB MATRIX IN HFSCF0, BUT WITH MQN STRUCTURE AS WELL.    C
C**********************************************************************C
      INCLUDE 'param.h'
C
      CHARACTER*4 HMLT
C
      DIMENSION DKAB(MNU,MKP+1,MKP+1),DKCD(MNU,MKP+1,MKP+1)
      DIMENSION RJLLLL(MB2,MNU),RJLLSS(MB2,MNU),
     &          RJSSLL(MB2,MNU),RJSSSS(MB2,MNU)
      DIMENSION XLLLL(MB2),XSSLL(MB2),XLLSS(MB2),XSSSS(MB2)
C
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VUEH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),WDIR(MDM,MDM),
     &           WXCH(MDM,MDM),CPLE(MDM,MDM)
C
      COMMON/B0IJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ,IBAS,JBAS
      COMMON/B1QN/EXL(MBS,4),MQN(4),KQN(4),LQN(4),NBAS(4),MAXM
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/DENS/DENC,DENO,DENT
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VUEH,QDIR,
     &            QXCH,WDIR,WXCH,CPLE
      COMMON/PRMS/HMLT,ICLC,INEW,IATM,ILIN
      COMMON/RCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/SHLL/ACFF,BCFF,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
      COMMON/T2EL/F2ES(5,6),T2ES(5,6),N2EB(5,6),N2EI(5,6),N2ES(5,6)
      COMMON/TSCF/TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRR,TCC1,TCC2,TCMC,
     &            TB1B,TB1R,TB1F,TB1M,TBEC,TBRM,TBRR,TBC1,TBC2,TBMC,
     &            TQMX,THMX,TC1T,TC2T,TB1T,TB2T,TEIG,TSCR,TTOT,
     &            TC1S,TC2S,TB1S,TB2S
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
C
C     ANGULAR FACTOR SENSITIVITY PARAMETER
      DATA SENS/1.0D-10/
C
      CALL CPU_TIME(TBCH1)
C
C     INTEGRAL SKIPPING ON MOLECULAR GROUP SYMMETRY CLASS BASIS
      IF(IATM.EQ.1.OR.ILIN.EQ.1) THEN
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
        IF(KQN(1).LT.0) THEN
          LQN(1) =-KQN(1)-1
        ELSE
          LQN(1) = KQN(1)
        ENDIF
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
        IF(KQN(2).LT.0) THEN
          LQN(2) =-KQN(2)-1
        ELSE
          LQN(2) = KQN(2)
        ENDIF
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
        IF(KQN(3).LT.0) THEN
          LQN(3) =-KQN(3)-1
        ELSE
          LQN(3) = KQN(3)
        ENDIF
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
        IF(KQN(4).LT.0) THEN
          LQN(4) =-KQN(4)-1
        ELSE
          LQN(4) = KQN(4)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK D
        NBAS(4) = NFNC(LQN(4),IZ)
        DO LBAS=1,NBAS(4)
          EXL(LBAS,4) = BEXL(LBAS,LQN(4),IZ)
        ENDDO
C
C     NUMBER OF BASIS FUNCTION COMBINATIONS IN (CD) BLOCK
      MAXM = NBAS(3)*NBAS(4)
C
C**********************************************************************C
C     PREPARE INTERMEDIATE DATA FOR USE IN RKCLM1                      C
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
      CALL ANGCLM1(DKAB,DKCD,KQN,LQN,ISEL)
      CALL CPU_TIME(T2)
      TC1B = TC1B+T2-T1
C
C     EXIT THIS COMBINATION IF IT VIOLATES A SELECTION RULE
      IF(ISEL.EQ.0) GOTO 1001
C
C     BASIS SET INTERMEDIATES FOR THE KL-PAIRS
      CALL CPU_TIME(T1)
      CALL KLSET1(IZ)
      CALL CPU_TIME(T2)
      TC1B = TC1B+T2-T1
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS IN BLOCKS A AND B (INDEX 2000)         C
C**********************************************************************C
C
C     UPDATE COUNTER FOR NUMBER OF CLASSES
      IF(HMLT.EQ.'NORL') THEN
        N2EB(1,1) = N2EB(1,1) + 1
      ELSE
        N2EB(1,1) = N2EB(1,1) + 1
        N2EB(1,2) = N2EB(1,2) + 1
        N2EB(1,3) = N2EB(1,3) + 1
        N2EB(1,4) = N2EB(1,4) + 1
      ENDIF
C
      DO 2000 IBAS=1,NBAS(1)
      DO 2000 JBAS=1,NBAS(2)
C
C       UPDATE COUNTER FOR NUMBER OF INTEGRALS
        IF(HMLT.EQ.'NORL') THEN
          N2EI(1,1) = N2EI(1,1) + NBAS(3)*NBAS(4)
        ELSE
          N2EI(1,1) = N2EI(1,1) + NBAS(3)*NBAS(4)
          N2EI(1,2) = N2EI(1,2) + NBAS(3)*NBAS(4)
          N2EI(1,3) = N2EI(1,3) + NBAS(3)*NBAS(4)
          N2EI(1,4) = N2EI(1,4) + NBAS(3)*NBAS(4)
        ENDIF
C
C       BASIS SET INTERMEDIATES FOR THE IJ-PAIRS
        CALL CPU_TIME(T1)
        CALL IJSET1
        CALL CPU_TIME(T2)
        TC1B = TC1B+T2-T1
C
C       BATCH OF RADIAL INTEGRALS (EFFECTIVE INTERACTION STRENGTHS)
        CALL CPU_TIME(T1)
        CALL RKCLM1(RJLLLL,RJLLSS,RJSSLL,RJSSSS)
        CALL CPU_TIME(T2)
        TC1R = TC1R+T2-T1
C
C       CAN NOW CONTRACT THESE RADIAL INTEGRALS OVER ANGULAR COMPONENTS
C       OF G-SPINOR BASIS FUNCTIONS USING A TENSOR EXPANSION IN {L,Q}
C
C**********************************************************************C
C     LOOP OVER |MQN| MAGNITUDES FOR A,B,C,D BLOCKS (USE INDEX 3000)   C
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
      NAS = LRGE(IZ,KA,IMJA)+NSKP
      NBS = LRGE(IZ,KB,IMJB)+NSKP
      NCS = LRGE(IZ,KC,IMJC)+NSKP
      NDS = LRGE(IZ,KD,IMJD)+NSKP
C
C     APPLY ANGULAR MQN SELECTION RULE
      IF(MMJA-MMJB.NE.MMJD-MMJC) GOTO 4001
      IF(ISYM.EQ.1) THEN
        IF(ISGN1.EQ.ISGN2.AND.ISGN3.EQ.ISGN4) GOTO 4002
        IF(ISGN1.EQ.ISGN4.AND.ISGN3.EQ.ISGN2) GOTO 4002     
        GOTO 4001
      ENDIF
4002  CONTINUE
C
C     INITIALISE EFFECTIVE INTERACTION STRENGTH VALUES
      CALL CPU_TIME(T1)
      DO M=1,NBAS(3)*NBAS(4)
        XLLLL(M) = 0.0D0
        XLLSS(M) = 0.0D0
        XSSLL(M) = 0.0D0
        XSSSS(M) = 0.0D0
      ENDDO
C
C     ASSEMBLE INTEGRAL BY SUMMING OVER EXPANSION POWER L
      DO LTEN=1,NUNUM
C
C       ANGULAR COEFFICIENT
C       DFNOTE: ANGFAC CAN BE EXPRESSED AS THE SQUARE ROOT OF A 
C               RATIONAL NUMBER -- SMALL ERRORS AT PRESENT.
        ANGFAC = DKAB(LTEN,IMJA,IMJB)*DKCD(LTEN,IMJC,IMJD)
C
C       SCREENING OF INTEGRAL BASED ON ANGULAR COEFFICIENT
        IF(DABS(ANGFAC).LE.SENS) GOTO 4003
C
        IF(HMLT.EQ.'NORL') THEN
C       NON-RELATIVISTIC HAMILTONIAN
C
          DO M=1,NBAS(3)*NBAS(4)
            XLLLL(M) = XLLLL(M) + ANGFAC*RJLLLL(M,LTEN)
          ENDDO
C
        ELSE
C       RELATIVISTIC HAMILTONIAN
C
          DO M=1,NBAS(3)*NBAS(4)
            XLLLL(M) = XLLLL(M) + ANGFAC*RJLLLL(M,LTEN)
            XLLSS(M) = XLLSS(M) + ANGFAC*RJLLSS(M,LTEN)
            XSSLL(M) = XSSLL(M) + ANGFAC*RJSSLL(M,LTEN)
            XSSSS(M) = XSSSS(M) + ANGFAC*RJSSSS(M,LTEN)
          ENDDO
C
        ENDIF
C
C       SKIP POINT FOR ANGULAR SCREENING
4003    CONTINUE
C
      ENDDO
C
C     FULL-INTEGRAL CONSTRUCTION COMPLETE
      CALL CPU_TIME(T2)
      TC1F = TC1F+T2-T1
C
C**********************************************************************C
C     ADD THIS BATCH OF R-INTEGRALS TO CLOSED-SHELL COULOMB MATRIX     C
C -------------------------------------------------------------------- C
C     DFNOTE: BE AWARE THAT THERE ARE ISSUES WITH NUMERICAL ACCURACY   C
C             IN THIS METHOD, STEMMING FROM ACCUMULATING ERRORS IN THE C
C             FLOATING-POINT REPRESENTATION OF INTEGERS AND RATIONALS. C
C             (EG. ANGFAC^2 IS A RATIONAL NUMBER).                     C
C**********************************************************************C
C
C     NON-RELATIVISTIC HAMILTONIAN
      IF(HMLT.EQ.'NORL') THEN
C
C       DIRECT CONTRIBUTIONS
        IF(ISYM.EQ.1) THEN
          IF(ISGN1.NE.ISGN2.OR.ISGN3.NE.ISGN4) GOTO 5001
        ENDIF
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NAL+IBAS,NBL+JBAS) = GDIR(NAL+IBAS,NBL+JBAS)
     &                          +      XLLLL(M)*DENT(NCL+KBAS,NDL+LBAS)
C
          ENDDO
        ENDDO
5001    CONTINUE
C
C       EXCHANGE CONTRIBUTIONS
        IF(ISYM.EQ.1) THEN
          IF(ISGN1.NE.ISGN4.OR.ISGN3.NE.ISGN2) GOTO 5002
        ENDIF
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NAL+IBAS,NDL+LBAS) = GXCH(NAL+IBAS,NDL+LBAS) 
     &                          +      XLLLL(M)*DENT(NCL+KBAS,NBL+JBAS)
C
          ENDDO
        ENDDO
5002    CONTINUE
C
C     RELATIVISTIC HAMILTONIAN
      ELSE
C
C       DIRECT CONTRIBUTIONS
        IF(ISYM.EQ.1) THEN
          IF(ISGN1.NE.ISGN2.OR.ISGN3.NE.ISGN4) GOTO 5003
        ENDIF
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NAL+IBAS,NBL+JBAS) = GDIR(NAL+IBAS,NBL+JBAS)
     &                          +      XLLLL(M)*DENT(NCL+KBAS,NDL+LBAS)
     &                          +      XLLSS(M)*DENT(NCS+KBAS,NDS+LBAS)
C
            GDIR(NAS+IBAS,NBS+JBAS) = GDIR(NAS+IBAS,NBS+JBAS)
     &                          +      XSSLL(M)*DENT(NCL+KBAS,NDL+LBAS)
     &                          +      XSSSS(M)*DENT(NCS+KBAS,NDS+LBAS)
C
          ENDDO
        ENDDO
5003    CONTINUE
C
C       EXCHANGE CONTRIBUTIONS
        IF(ISYM.EQ.1) THEN
          IF(ISGN1.NE.ISGN4.OR.ISGN3.NE.ISGN2) GOTO 5004
        ENDIF
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NAL+IBAS,NDL+LBAS) = GXCH(NAL+IBAS,NDL+LBAS)
     &                          +      XLLLL(M)*DENT(NCL+KBAS,NBL+JBAS)
C
            GXCH(NAL+IBAS,NDS+LBAS) = GXCH(NAL+IBAS,NDS+LBAS)
     &                          +      XLLSS(M)*DENT(NCS+KBAS,NBL+JBAS)
C
            GXCH(NAS+IBAS,NDL+LBAS) = GXCH(NAS+IBAS,NDL+LBAS)
     &                          +      XSSLL(M)*DENT(NCL+KBAS,NBS+JBAS)
C
            GXCH(NAS+IBAS,NDS+LBAS) = GXCH(NAS+IBAS,NDS+LBAS)
     &                          +      XSSSS(M)*DENT(NCS+KBAS,NBS+JBAS)
C
          ENDDO
        ENDDO
5004    CONTINUE
C
      ENDIF
C
C**********************************************************************C
C     ADD THIS BATCH OF R-INTEGRALS TO OPEN-SHELL COULOMB MATRIX       C
C**********************************************************************C
C
      IF(NOPN.EQ.0) GOTO 5100
C
C     NON-RELATIVISTIC HAMILTONIAN
      IF(HMLT.EQ.'NORL') THEN
C
C       DIRECT CONTRIBUTIONS
        IF(ISYM.EQ.1) THEN
          IF(ISGN1.NE.ISGN2.OR.ISGN3.NE.ISGN4) GOTO 5005
        ENDIF
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NAL+IBAS,NBL+JBAS) = QDIR(NAL+IBAS,NBL+JBAS)
     &                          + ACFF*XLLLL(M)*DENO(NCL+KBAS,NDL+LBAS)
C
          ENDDO
        ENDDO
5005    CONTINUE
C
C       EXCHANGE CONTRIBUTIONS
        IF(ISYM.EQ.1) THEN
          IF(ISGN1.NE.ISGN4.OR.ISGN3.NE.ISGN2) GOTO 5006
        ENDIF
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C           
            QXCH(NAL+IBAS,NDL+LBAS) = QXCH(NAL+IBAS,NDL+LBAS)
     &                          + BCFF*XLLLL(M)*DENO(NCL+KBAS,NBL+JBAS)
C
          ENDDO
        ENDDO
5006    CONTINUE
C
C     RELATIVISTIC HAMILTONIAN
      ELSE
C
C       DIRECT CONTRIBUTIONS
        IF(ISYM.EQ.1) THEN
          IF(ISGN1.NE.ISGN2.OR.ISGN3.NE.ISGN4) GOTO 5007
        ENDIF
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NAL+IBAS,NBL+JBAS) = QDIR(NAL+IBAS,NBL+JBAS)
     &                          + ACFF*XLLLL(M)*DENO(NCL+KBAS,NDL+LBAS)
     &                          + ACFF*XLLSS(M)*DENO(NCS+KBAS,NDS+LBAS)
C
            QDIR(NAS+IBAS,NBS+JBAS) = QDIR(NAS+IBAS,NBS+JBAS)
     &                          + ACFF*XSSLL(M)*DENO(NCL+KBAS,NDL+LBAS)
     &                          + ACFF*XSSSS(M)*DENO(NCS+KBAS,NDS+LBAS)
C
          ENDDO
        ENDDO
5007    CONTINUE
C
C       EXCHANGE CONTRIBUTIONS
        IF(ISYM.EQ.1) THEN
          IF(ISGN1.NE.ISGN4.OR.ISGN3.NE.ISGN2) GOTO 5008
        ENDIF
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NAL+IBAS,NDL+LBAS) = QXCH(NAL+IBAS,NDL+LBAS)
     &                          + BCFF*XLLLL(M)*DENO(NCL+KBAS,NBL+JBAS)
C
            QXCH(NAL+IBAS,NDS+LBAS) = QXCH(NAL+IBAS,NDS+LBAS)
     &                          + BCFF*XLLSS(M)*DENO(NCS+KBAS,NBL+JBAS)
C
            QXCH(NAS+IBAS,NDL+LBAS) = QXCH(NAS+IBAS,NDL+LBAS)
     &                          + BCFF*XSSLL(M)*DENO(NCL+KBAS,NBS+JBAS)
            
            QXCH(NAS+IBAS,NDS+LBAS) = QXCH(NAS+IBAS,NDS+LBAS)
     &                          + BCFF*XSSSS(M)*DENO(NCS+KBAS,NBS+JBAS)
C
          ENDDO
        ENDDO
5008    CONTINUE
C
      ENDIF
C
5100  CONTINUE
C
C     MATRIX MULTIPLICATION STEP COMPLETE
      CALL CPU_TIME(T3)
      TC1M = TC1M+T3-T2
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
      IF(HMLT.EQ.'NORL') THEN
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
