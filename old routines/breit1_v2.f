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
C -------------------------------------------------------------------- C
C  DFNOTE: NOT WORKING YET...                                          C
C**********************************************************************C
      INCLUDE 'param.h'
C
      DIMENSION DKAB(MNU,MKP+1,MKP+1),DKCD(MNU,MKP+1,MKP+1)
      DIMENSION EMAT(MNU,8),ANGFAC(8)
      DIMENSION RJLSSL(MB2,MNU,2),RJSLLS(MB2,MNU,2),
     &          RJSLSL(MB2,MNU,2),RJLSLS(MB2,MNU,2)
      DIMENSION XLSSL(MB2),XSLLS(MB2),XSLSL(MB2),XLSLS(MB2)
      
      dimension jqn(4)
      dimension zssll(10,10),zllss(10,10),
     &          zlssl(10,10),zslls(10,10),ztotal(10,10)
C
      CHARACTER*2 KLAB
C
      complex*16 sum
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VUEH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),WDIR(MDM,MDM),
     &           WXCH(MDM,MDM),CPLE(MDM,MDM)
      
      complex*16 bsav(mdm,mdm,10)
C
      COMMON/B0IJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ,IBAS,JBAS
      COMMON/B1QN/EXL(MBS,4),MQN(4),KQN(4),LQN(4),NBAS(4),MAXM
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/DENS/DENC,DENO,DENT
      COMMON/EIGC/COEF
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VUEH,QDIR,
     &            QXCH,WDIR,WXCH,CPLE
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
c     empty bsav
      do i=1,ndim
        do j=1,ndim
          do iocc=1,10
            bsav(i,j,iocc) = dcmplx(0.0d0,0.0d0)
          enddo
        enddo
      enddo      
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
        JQN(1) = 2*IABS(KQN(1))-1
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
        JQN(2) = 2*IABS(KQN(2))-1
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
        JQN(3) = 2*IABS(KQN(3))-1
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
        JQN(4) = 2*IABS(KQN(4))-1
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
      
C      if(kqn(1).ne.kqn(4)) goto 1001
C      if(kqn(2).ne.kqn(3)) goto 1001
C
C     ANGULAR COEFFICIENTS
      CALL CPU_TIME(T1)
      CALL ANGBRT1(DKAB,DKCD,EMAT,KQN,LQN,ISEL)
      CALL CPU_TIME(T2)
      TB1B = TB1B+T2-T1
C
C      if(isel.ne.0) then
C        write(*,*) kqn(1),kqn(2),kqn(3),kqn(4),'|',nunum,'|',
C     &                                        (nus(lten),lten=1,nunum)
C        do lten=1,nunum
C          write(*,*) nus(lten),emat(lten,1),emat(lten,3),
C     &                         emat(lten,5),emat(lten,7)
C          write(*,*) '           ',emat(lten,2),emat(lten,4),
C     &                             emat(lten,5),emat(lten,8)
C        enddo
C      endif
C      if(kqn(1).eq.-1.and.kqn(2).eq.1.and.kqn(3).eq.1.and.kqn(4).eq.-1)
C     & then
C       nunum=2
C       nui=0
C       nuf=2
C       nus(1)=0
C       nus(2)=2
C      endif
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
        CALL RKBRT1(RJLSSL,RJSLLS,RJSLSL,RJLSLS)
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
C     RESET CONTRACTED RADIAL ARRAYS
      CALL CPU_TIME(T1)
      DO M=1,NBAS(3)*NBAS(4)
        XLSSL(M) = 0.0D0
        XSLLS(M) = 0.0D0
        XSLSL(M) = 0.0D0
        XLSLS(M) = 0.0D0
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
        IF(DABS(ANGSUM).LE.SENS) GOTO 4003
C
        DO M=1,NBAS(3)*NBAS(4)
C
C         FACTORS LEADING TO BXCH(LL)
          XLSSL(M) = XLSSL(M) + ANGFAC(8)*RJLSSL(M,LTEN,1)
     &                        + ANGFAC(7)*RJLSSL(M,LTEN,2)
C
C         FACTORS LEADING TO BXCH(SS)
          XSLLS(M) = XSLLS(M) + ANGFAC(6)*RJSLLS(M,LTEN,1)
     &                        + ANGFAC(5)*RJSLLS(M,LTEN,2)
C
C         FACTORS LEADING TO BXCH(SL)
          XSLSL(M) = XSLSL(M) + ANGFAC(4)*RJSLSL(M,LTEN,1)
     &                        + ANGFAC(3)*RJSLSL(M,LTEN,2)
C
C         FACTORS LEADING TO BXCH(LS)
          XLSLS(M) = XLSLS(M) + ANGFAC(2)*RJLSLS(M,LTEN,1)
     &                        + ANGFAC(1)*RJLSLS(M,LTEN,2)
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
          BXCH(NAL+IBAS,NDL+LBAS) = BXCH(NAL+IBAS,NDL+LBAS)
     &                          +      XLSSL(M)*DENT(NBS+JBAS,NCS+KBAS)
C
          BXCH(NAL+IBAS,NDS+LBAS) = BXCH(NAL+IBAS,NDS+LBAS)
     &                          +      XLSLS(M)*DENT(NBS+JBAS,NCL+KBAS)
C
          BXCH(NAS+IBAS,NDL+LBAS) = BXCH(NAS+IBAS,NDL+LBAS)
     &                          +      XSLSL(M)*DENT(NBL+JBAS,NCS+KBAS)

          BXCH(NAS+IBAS,NDS+LBAS) = BXCH(NAS+IBAS,NDS+LBAS)
     &                          +      XSLLS(M)*DENT(NBL+JBAS,NCL+KBAS)
C
        ENDDO
      ENDDO

C
C     FOR SYMMETRY TYPE REDUCTION
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
          do iocc=1,10
          
          iad = iocc+NSKP
           
          bsav(NAL+IBAS,NDL+LBAS,iocc) = bsav(NAL+IBAS,NDL+LBAS,iocc)
     &     + XLSSL(M)*DCONJG(COEF(NBS+JBAS,iad))*COEF(NCS+KBAS,iad)
C
          bsav(NAL+IBAS,NDS+LBAS,iocc) = bsav(NAL+IBAS,NDS+LBAS,iocc)
     &     + XLSLS(M)*DCONJG(COEF(NBS+JBAS,iad))*COEF(NCL+KBAS,iad)
C
          bsav(NAS+IBAS,NDL+LBAS,iocc) = bsav(NAS+IBAS,NDL+LBAS,iocc)
     &     + XSLSL(M)*DCONJG(COEF(NBL+JBAS,iad))*COEF(NCS+KBAS,iad)

          bsav(NAS+IBAS,NDS+LBAS,iocc) = bsav(NAS+IBAS,NDS+LBAS,iocc)
     &     + XSLLS(M)*DCONJG(COEF(NBL+JBAS,iad))*COEF(NCL+KBAS,iad)

          enddo
        ENDDO
      ENDDO

C
C**********************************************************************C
C     ADD THIS BATCH OF R-INTEGRALS TO OPEN-SHELL BREIT MATRIX.        C
C     THIS ALSO REQUIRES THE CLOSED-SHELL DIRECT MATRIX ELEMENTS.      C
C**********************************************************************C
C
      IF(NOPN.EQ.0) GOTO 5100
C
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
C
C         DIRECT CONTRIBUTIONS TO CLOSED-SHELL MATRIX
          BDIR(NAL+IBAS,NBS+JBAS) = BDIR(NAL+IBAS,NBS+JBAS)
     &                          +      XLSLS(M)*DENT(NCL+KBAS,NDS+LBAS)
     &                          +      XLSSL(M)*DENT(NCS+KBAS,NDL+LBAS)
C
          BDIR(NAS+IBAS,NBL+JBAS) = BDIR(NAS+IBAS,NBL+JBAS)
     &                          +      XSLLS(M)*DENT(NCL+KBAS,NDS+LBAS)
     &                          +      XSLSL(M)*DENT(NCS+KBAS,NDL+LBAS)
C
C         DIRECT CONTRIBUTIONS TO OPEN-SHELL MATRIX
          WDIR(NAL+IBAS,NBS+JBAS) = WDIR(NAL+IBAS,NBS+JBAS)
     &                          + ACFF*XLSLS(M)*DENO(NCL+KBAS,NDS+LBAS)
     &                          + ACFF*XLSSL(M)*DENO(NCS+KBAS,NDL+LBAS)
C
          WDIR(NAS+IBAS,NBL+JBAS) = WDIR(NAS+IBAS,NBL+JBAS)
     &                          + ACFF*XSLLS(M)*DENO(NCL+KBAS,NDS+LBAS)
     &                          + ACFF*XSLSL(M)*DENO(NCS+KBAS,NDL+LBAS)
C
C         EXCHANGE CONTRIBUTIONS TO OPEN-SHELL MATRIX
          WXCH(NAL+IBAS,NDL+LBAS) = WXCH(NAL+IBAS,NDL+LBAS)
     &                          + BCFF*XLSSL(M)*DENO(NBS+JBAS,NCS+KBAS)
C
          WXCH(NAL+IBAS,NDS+LBAS) = WXCH(NAL+IBAS,NDS+LBAS)
     &                          + BCFF*XLSLS(M)*DENO(NBS+JBAS,NCL+KBAS)
C
          WXCH(NAS+IBAS,NDL+LBAS) = WXCH(NAS+IBAS,NDL+LBAS)
     &                          + BCFF*XSLSL(M)*DENO(NBL+JBAS,NCS+KBAS)

          WXCH(NAS+IBAS,NDS+LBAS) = WXCH(NAS+IBAS,NDS+LBAS)
     &                          + BCFF*XSLLS(M)*DENO(NBL+JBAS,NCL+KBAS)
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
C     CALCULATION OF BREIT ENERGY
      ELLSS = 0.0D0
      ESSLL = 0.0D0
      ESLLS = 0.0D0
      ELSSL = 0.0D0
      DO IL=1,NSKP
        DO JL=1,NSKP
C
C         SMALL-COMPONENT LABELS
          IS = IL+NSKP
          JS = JL+NSKP
C
          ELLSS = ELLSS - 0.5D0*BXCH(IL,JL)*DENT(IL,JL)
          ESSLL = ESSLL - 0.5D0*BXCH(IS,JS)*DENT(IS,JS)
          ESLLS = ESLLS - 0.5D0*BXCH(IS,JL)*DENT(IS,JL)
          ELSSL = ELSSL - 0.5D0*BXCH(IL,JS)*DENT(IL,JS)
C
        ENDDO
      ENDDO
C
      do iocc=1,10
        do jocc=1,10
          jad = jocc+NSKP
          zllss(iocc,jocc) = 0.0d0
          zssll(iocc,jocc) = 0.0d0
          zslls(iocc,jocc) = 0.0d0
          zlssl(iocc,jocc) = 0.0d0
          do il=1,NSKP
            do jl=1,NSKP
              is = il+NSKP
              js = jl+NSKP
              zllss(iocc,jocc) = zllss(iocc,jocc) 
     &            - 0.5D0*bsav(il,jl,iocc)*COEF(il,jad)*COEF(jl,jad)
              zssll(iocc,jocc) = zssll(iocc,jocc) 
     &            - 0.5D0*bsav(is,js,iocc)*COEF(is,jad)*COEF(js,jad)
              zslls(iocc,jocc) = zslls(iocc,jocc) 
     &            - 0.5D0*bsav(is,jl,iocc)*COEF(is,jad)*COEF(jl,jad)
              zlssl(iocc,jocc) = zlssl(iocc,jocc) 
     &            - 0.5D0*bsav(il,js,iocc)*COEF(il,jad)*COEF(js,jad)
            enddo
          enddo
          ztotal(iocc,jocc) = zllss(iocc,jocc) + zssll(iocc,jocc)
     &                      + zslls(iocc,jocc) + zlssl(iocc,jocc)
c          write(*,*) iocc,jocc,zllss(iocc,jocc),zssll(iocc,jocc),
c     &               zslls(iocc,jocc),zlssl(iocc,jocc),ztotal(iocc,jocc)
        enddo
      enddo
c
601   FORMAT(1X,'(KA,KB)',9X,
     &              'E(LL|SS)',8X,'E(LS|SL)',8X,'E(SL|LS)',8X,
     &              'E(SS|LL)',10X,'E(TOT)')
602   FORMAT(1X,'(',A,',',A,')',3X,F14.10,2X,F14.10,2X,F14.10,2X,
     &                               F14.10,2X,F14.10,2X,F14.10,2X)
603   FORMAT(1X,'total',5X,F14.10,2X,F14.10,2X,F14.10,2X,F14.10,2X,
     &                                         F14.10,2X,F14.10,2X)

C
C     PRINT A HEADER
      WRITE(*,*) REPEAT('=',88)
      WRITE(*,601)
      WRITE(*,*) REPEAT('-',88)

      ellssnet  = 0.0d0
      essllnet  = 0.0d0
      esllsnet  = 0.0d0
      elsslnet  = 0.0d0
      etotalnet = 0.0d0


c     (s_1/2|s_1/2)
      ellss  = 0.0d0
      essll  = 0.0d0
      eslls  = 0.0d0
      elssl  = 0.0d0
      etotal = 0.0d0     
      do iocc=1,4
        do jocc=1,4
          ellss  = ellss  + zllss(iocc,jocc)
          essll  = essll  + zssll(iocc,jocc)
          eslls  = eslls  + zslls(iocc,jocc)
          elssl  = elssl  + zlssl(iocc,jocc)
          etotal = etotal + ztotal(iocc,jocc)
        enddo
      enddo
      write(*,602) KLAB(-1),KLAB(-1),ellss,elssl,eslls,essll,etotal
c
      ellssnet  = ellssnet + ellss
      essllnet  = essllnet + essll
      esllsnet  = esllsnet + eslls
      elsslnet  = elsslnet + elssl
      etotalnet = etotalnet + etotal
C
c     (s_1/2|p_1/2)
      ellss  = 0.0d0
      essll  = 0.0d0
      eslls  = 0.0d0
      elssl  = 0.0d0
      etotal = 0.0d0     
      do iocc=1,4
        do jocc=5,6
          ellss  = ellss  + zllss(iocc,jocc)
          essll  = essll  + zssll(iocc,jocc)
          eslls  = eslls  + zslls(iocc,jocc)
          elssl  = elssl  + zlssl(iocc,jocc)
          etotal = etotal + ztotal(iocc,jocc)
        enddo
      enddo
      write(*,602) KLAB(-1),KLAB(+1),ellss,elssl,eslls,essll,etotal
c
      ellssnet  = ellssnet + ellss
      essllnet  = essllnet + essll
      esllsnet  = esllsnet + eslls
      elsslnet  = elsslnet + elssl
      etotalnet = etotalnet + etotal
C
c     (s_1/2|p_3/2)
      ellss  = 0.0d0
      essll  = 0.0d0
      eslls  = 0.0d0
      elssl  = 0.0d0
      etotal = 0.0d0     
      do iocc=1,4
        do jocc=7,10
          ellss  = ellss  + zllss(iocc,jocc)
          essll  = essll  + zssll(iocc,jocc)
          eslls  = eslls  + zslls(iocc,jocc)
          elssl  = elssl  + zlssl(iocc,jocc)
          etotal = etotal + ztotal(iocc,jocc)
        enddo
      enddo
      write(*,602) KLAB(-1),KLAB(-2),ellss,elssl,eslls,essll,etotal
c
      ellssnet  = ellssnet + ellss
      essllnet  = essllnet + essll
      esllsnet  = esllsnet + eslls
      elsslnet  = elsslnet + elssl
      etotalnet = etotalnet + etotal
C
c     (p_1/2|s_1/2)
      ellss  = 0.0d0
      essll  = 0.0d0
      eslls  = 0.0d0
      elssl  = 0.0d0
      etotal = 0.0d0     
      do iocc=5,6
        do jocc=1,4
          ellss  = ellss  + zllss(iocc,jocc)
          essll  = essll  + zssll(iocc,jocc)
          eslls  = eslls  + zslls(iocc,jocc)
          elssl  = elssl  + zlssl(iocc,jocc)
          etotal = etotal + ztotal(iocc,jocc)
        enddo
      enddo
      write(*,602) KLAB(+1),KLAB(-1),ellss,elssl,eslls,essll,etotal
c
      ellssnet  = ellssnet + ellss
      essllnet  = essllnet + essll
      esllsnet  = esllsnet + eslls
      elsslnet  = elsslnet + elssl
      etotalnet = etotalnet + etotal
C
c     (p_1/2|p_1/2)
      ellss  = 0.0d0
      essll  = 0.0d0
      eslls  = 0.0d0
      elssl  = 0.0d0
      etotal = 0.0d0     
      do iocc=5,6
        do jocc=5,6
          ellss  = ellss  + zllss(iocc,jocc)
          essll  = essll  + zssll(iocc,jocc)
          eslls  = eslls  + zslls(iocc,jocc)
          elssl  = elssl  + zlssl(iocc,jocc)
          etotal = etotal + ztotal(iocc,jocc)
        enddo
      enddo
      write(*,602) KLAB(+1),KLAB(+1),ellss,elssl,eslls,essll,etotal
c
      ellssnet  = ellssnet + ellss
      essllnet  = essllnet + essll
      esllsnet  = esllsnet + eslls
      elsslnet  = elsslnet + elssl
      etotalnet = etotalnet + etotal
C
c     (p_1/2|p_3/2)
      ellss  = 0.0d0
      essll  = 0.0d0
      eslls  = 0.0d0
      elssl  = 0.0d0
      etotal = 0.0d0     
      do iocc=5,6
        do jocc=7,10
          ellss  = ellss  + zllss(iocc,jocc)
          essll  = essll  + zssll(iocc,jocc)
          eslls  = eslls  + zslls(iocc,jocc)
          elssl  = elssl  + zlssl(iocc,jocc)
          etotal = etotal + ztotal(iocc,jocc)
        enddo
      enddo
      write(*,602) KLAB(+1),KLAB(-2),ellss,elssl,eslls,essll,etotal
c
      ellssnet  = ellssnet + ellss
      essllnet  = essllnet + essll
      esllsnet  = esllsnet + eslls
      elsslnet  = elsslnet + elssl
      etotalnet = etotalnet + etotal
C
c     (p_3/2|s_1/2)
      ellss  = 0.0d0
      essll  = 0.0d0
      eslls  = 0.0d0
      elssl  = 0.0d0
      etotal = 0.0d0     
      do iocc=7,10
        do jocc=1,4
          ellss  = ellss  + zllss(iocc,jocc)
          essll  = essll  + zssll(iocc,jocc)
          eslls  = eslls  + zslls(iocc,jocc)
          elssl  = elssl  + zlssl(iocc,jocc)
          etotal = etotal + ztotal(iocc,jocc)
        enddo
      enddo
      write(*,602) KLAB(-2),KLAB(-1),ellss,elssl,eslls,essll,etotal
c
      ellssnet  = ellssnet + ellss
      essllnet  = essllnet + essll
      esllsnet  = esllsnet + eslls
      elsslnet  = elsslnet + elssl
      etotalnet = etotalnet + etotal
C
c     (p_3/2|p_1/2)
      ellss  = 0.0d0
      essll  = 0.0d0
      eslls  = 0.0d0
      elssl  = 0.0d0
      etotal = 0.0d0     
      do iocc=7,10
        do jocc=5,6
          ellss  = ellss  + zllss(iocc,jocc)
          essll  = essll  + zssll(iocc,jocc)
          eslls  = eslls  + zslls(iocc,jocc)
          elssl  = elssl  + zlssl(iocc,jocc)
          etotal = etotal + ztotal(iocc,jocc)
        enddo
      enddo
      write(*,602) KLAB(-2),KLAB(+1),ellss,elssl,eslls,essll,etotal
c
      ellssnet  = ellssnet + ellss
      essllnet  = essllnet + essll
      esllsnet  = esllsnet + eslls
      elsslnet  = elsslnet + elssl
      etotalnet = etotalnet + etotal
C
c     (p_3/2|p_3/2)
      ellss  = 0.0d0
      essll  = 0.0d0
      eslls  = 0.0d0
      elssl  = 0.0d0
      etotal = 0.0d0     
      do iocc=7,10
        do jocc=7,10
          ellss  = ellss  + zllss(iocc,jocc)
          essll  = essll  + zssll(iocc,jocc)
          eslls  = eslls  + zslls(iocc,jocc)
          elssl  = elssl  + zlssl(iocc,jocc)
          etotal = etotal + ztotal(iocc,jocc)
        enddo
      enddo
      write(*,602) KLAB(-2),KLAB(-2),ellss,elssl,eslls,essll,etotal
c
      ellssnet  = ellssnet + ellss
      essllnet  = essllnet + essll
      esllsnet  = esllsnet + eslls
      elsslnet  = elsslnet + elssl
      etotalnet = etotalnet + etotal

      WRITE(*,*) REPEAT('-',88)
      write(*,603) ellssnet,elsslnet,esllsnet,essllnet,etotalnet
      WRITE(*,*) REPEAT('=',88)
C
C     RECORD CPU TIME AT END OF BATCH AND ADD TO APPROPRIATE COUNTER
      CALL CPU_TIME(TBCH2)
      T2ES(1,5) = T2ES(1,5)+TBCH2-TBCH1
      
      STOP
C
      RETURN
      END
