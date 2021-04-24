      SUBROUTINE BREIT2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           BBBBBBB  RRRRRRR  EEEEEEEE IIII TTTTTTTT 222222            C
C           BB    BB RR    RR EE        II     TT   22    22           C
C           BB    BB RR    RR EE        II     TT         22           C
C           BBBBBBB  RR    RR EEEEEE    II     TT        22            C
C           BB    BB RRRRRRR  EE        II     TT      22              C
C           BB    BB RR    RR EE        II     TT    22                C
C           BBBBBBB  RR    RR EEEEEEEE IIII    TT   22222222           C
C                                                                      C
C -------------------------------------------------------------------- C
C  BREIT2 GENERATES ALL TWO-CENTRE BREIT INTERACTION MATRIX ELEMENTS   C
C  WHICH ARISE FROM A CLASSICAL PAIR OF ONE-CENTRE CURRENT DENSITIES,  C
C  BY FIRST REDUCING SOURCES TO IDEAL MAGNETIC DIPOLES ON EACH CENTRE. C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*80 TITLE
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
      DIMENSION RAB(MCT,MCT)
C
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VANM(MDM,MDM),
     &           VSLF(MDM,MDM),VUEH(MDM,MDM),VWKR(MDM,MDM),
     &           VKSB(MDM,MDM),QDIR(MDM,MDM),QXCH(MDM,MDM),
     &           WDIR(MDM,MDM),WXCH(MDM,MDM),CPLE(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 SMAT(MDM,MDM,3)
      COMPLEX*16 UMAT(MDM,MDM,3)
      COMPLEX*16 ELS11(MB2,MEQ),ELS21(MB2,MEQ)
      COMPLEX*16 ESL11(MB2,MEQ),ESL21(MB2,MEQ)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/DENS/DENC,DENO,DENT
      COMMON/EILS/EILSFL(MFL,12),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EISL/EISLFL(MFL,12),IADISL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/LSHF/SHLEV(4),SHLV,ILEV
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VANM,VSLF,
     &            VUEH,VWKR,VKSB,QDIR,QXCH,WDIR,WXCH,CPLE
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/QNMS/LABICN(MDM),LABKQN(MDM),LABMQN(MDM)
C
C     SKIP MULTI-CENTRE CONTRIBUTIONS IN STAGES 1 AND 2
      IF(ILEV.LT.3) RETURN
C
C     DISTANCE BETWEEN FUNCTION ORIGINS
      DO IA=1,NCNT
        XA = BXYZ(1,IA)
        YA = BXYZ(2,IA)
        ZA = BXYZ(3,IA)
        DO IB=1,NCNT
          XB = BXYZ(1,IB)
          YB = BXYZ(2,IB)
          ZB = BXYZ(3,IB)
          RAB(IA,IB) = DSQRT((XB-XA)**2 + (YB-YA)**2 + (ZB-ZA)**2)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     BASIS PAIR CURRENT DENSITY (MATRIX S)                            C
C**********************************************************************C
C
C     INITIALISE THE OVERLAP ARRAY
      DO I=1,NDIM
        DO J=1,NDIM
          DO IQ=1,3
            SMAT(I,J,IQ) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
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
C     IF(ICNTA.NE.ICNTB) GOTO 1001

      WRITE(*,*) ICNTA,ICNTB
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,ICNTA)
        LQN(1) = LVAL(KQN(1))
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
        LQN(2) = LVAL(KQN(2))
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
        ENDDO
C
C     IF(KA.NE.3.OR.KB.NE.3) GOTO 2001
C
C     LOOP OVER |MQN(A)| VALUES
      DO 3000 MA=1,IABS(KQN(1))
        MJA    = 2*MA-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 3000 MB=1,IABS(KQN(2))
        MJB    = 2*MB-1
        MQN(2) = MJB
C
C     CALCULATE COMPONENT OFFSETS
      NA1L = LRGE(ICNTA,KA,2*MA-1)
      NA2L = LRGE(ICNTA,KA,2*MA  )
      NB1L = LRGE(ICNTB,KB,2*MB-1)
      NB2L = LRGE(ICNTB,KB,2*MB  )
C
      NA1S = NA1L+NSKP
      NA2S = NA2L+NSKP
      NB1S = NB1L+NSKP
      NB2S = NB2L+NSKP
C
C     THIS PHASE RELATES EQ22 AND EQ12 COEFFS TO EQ11 AND EQ21
      PHS = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
     &                    *DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
C
C     LOOP OVER CARTESIAN INDICES
      DO 4000 IQ=1,3
C
C       GENERATE ELSQ COEFFICIENTS
        CALL CPU_TIME(TDM1)
        CALL EQLSMK(ELS11,ELS21,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,IQ)
        CALL CPU_TIME(TDM2)
        TELS = TELS+TDM2-TDM1
C
C       GENERATE ESLQ COEFFICIENTS
        CALL CPU_TIME(TDM1)
        CALL EQSLMK(ESL11,ESL21,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,IQ)
        CALL CPU_TIME(TDM2)
        TESL = TESL+TDM2-TDM1
C
C       OVERLAP MATRIX ELEMENTS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M = M+1
C
C           EXPONENT COMBINATIONS
            EIJ = EXL(IBAS,1)+EXL(JBAS,2)
            ERT = DSQRT(PI/EIJ)**3
C
C           MATRIX ELEMENTS
            SMAT(NA1L+IBAS,NB1S+JBAS,IQ) = ERT*ELS11(M,1)
            SMAT(NA2L+IBAS,NB1S+JBAS,IQ) = ERT*ELS21(M,1)
            SMAT(NA1L+IBAS,NB2S+JBAS,IQ) =-PHS*ERT*DCONJG(ELS21(M,1))
            SMAT(NA2L+IBAS,NB2S+JBAS,IQ) = PHS*ERT*DCONJG(ELS11(M,1))
C
            SMAT(NA1S+IBAS,NB1L+JBAS,IQ) = ERT*ESL11(M,1)
            SMAT(NA2S+IBAS,NB1L+JBAS,IQ) = ERT*ESL21(M,1)
            SMAT(NA1S+IBAS,NB2L+JBAS,IQ) =-PHS*ERT*DCONJG(ESL21(M,1))
            SMAT(NA2S+IBAS,NB2L+JBAS,IQ) = PHS*ERT*DCONJG(ESL11(M,1))
C
          ENDDO
        ENDDO
C
C     END LOOP OVER CARTESIAN INDICES
4000  CONTINUE
C
C     END LOOP OVER CENTRES A AND B
3000  CONTINUE
2001  CONTINUE
2000  CONTINUE
1001  CONTINUE
1000  CONTINUE
C
C**********************************************************************C
C     PERFORM ONE-INDEX CONTRACTION (K) OVER SMAT AND DENS             C
C**********************************************************************C
C
      DO IQ=1,3
        DO JL=1,NDIM-NSKP
          JS = JL+NSKP
          IB = LABICN(JL)
C         J INDEX (NU) SHOULD BE ON CENTRE 1
          IF(IB.NE.1) GOTO 20
          DO LL=1,NDIM-NSKP
            LS = LL+NSKP
            ID = LABICN(LL)
C           L INDEX (TAU) SHOULD BE ON CENTRE 2
            IF(ID.NE.2) GOTO 30
            DO KL=1,NDIM-NSKP
              KS = KL+NSKP
              IC = LABICN(KL)
C             K INDEX (SIGMA) SHOULD BE ON CENTRE 2
              IF(IC.NE.2) GOTO 40
              UMAT(JL,LL,IQ) = UMAT(JL,LL,IQ)+SMAT(KS,JL,IQ)*DENT(KS,LL)
              UMAT(JL,LS,IQ) = UMAT(JL,LS,IQ)+SMAT(KS,JL,IQ)*DENT(KS,LS)
              UMAT(JS,LL,IQ) = UMAT(JS,LL,IQ)+SMAT(KL,JS,IQ)*DENT(KL,LL)
              UMAT(JS,LS,IQ) = UMAT(JS,LS,IQ)+SMAT(KL,JS,IQ)*DENT(KL,LS)
40            CONTINUE
            ENDDO
30          CONTINUE
          ENDDO
20        CONTINUE
        ENDDO
      ENDDO
C
C**********************************************************************C
C     PERFORM ONE-INDEX CONTRACTION (L) OVER SMAT AND UMAT             C
C**********************************************************************C
C
      DO IL=1,NDIM-NSKP
        IS = IL+NSKP
        IA = LABICN(IL)
C       I INDEX (MU) SHOULD BE ON CENTRE 1
        IF(IA.NE.1) GOTO 50
        DO JL=1,NDIM-NSKP
          JS  = JL+NSKP
          IB = LABICN(JL)
C         J INDEX (NU) SHOULD BE ON CENTRE 1
          IF(IB.NE.1) GOTO 60
          DO LL=1,NDIM-NSKP
            LS = LL+NSKP
            ID = LABICN(LL)
C           L INDEX (TAU) SHOULD BE ON CENTRE 2
            IF(ID.NE.2) GOTO 70
            DO IQ=1,3
              BXCH(IL,JL) = BXCH(IL,JL)
     &                   + 0.5D0*SMAT(IL,LS,IQ)*UMAT(JL,LS,IQ)/RAB(1,2)
              BXCH(IL,JS) = BXCH(IL,JS)
     &                   + 0.5D0*SMAT(IL,LS,IQ)*UMAT(JS,LS,IQ)/RAB(1,2)
              BXCH(IS,JL) = BXCH(IS,JL)
     &                   + 0.5D0*SMAT(IS,LL,IQ)*UMAT(JL,LL,IQ)/RAB(1,2)
              BXCH(IS,JS) = BXCH(IS,JS)
     &                   + 0.5D0*SMAT(IS,LL,IQ)*UMAT(JS,LL,IQ)/RAB(1,2)
            ENDDO
70          CONTINUE
          ENDDO
60        CONTINUE
        ENDDO
50      CONTINUE
      ENDDO
C
      RETURN
      END

