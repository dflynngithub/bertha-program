      SUBROUTINE SLFINT(ITER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           SSSSSS  LL       FFFFFFFF IIII NN    NN TTTTTTTT           C
C          SS    SS LL       FF        II  NNN   NN    TT              C
C          SS       LL       FF        II  NNNN  NN    TT              C
C           SSSSSS  LL       FFFFFF    II  NN NN NN    TT              C
C                SS LL       FF        II  NN  NNNN    TT              C
C          SS    SS LL       FF        II  NN   NNN    TT              C
C           SSSSSS  LLLLLLLL FF       IIII NN    NN    TT              C
C                                                                      C
C -------------------------------------------------------------------- C
C  SLFINT CONSTRUCTS THE MULTI-CENTRE ELECTRON SELF-INTERACTION        C
C  MATRIX ELEMENTS BASED ON THE LOW-ENERGY QUINEY MODEL 2019.          C
C -------------------------------------------------------------------- C
C ▶ THE `LL' COMPONENT OF THIS CAN BE USED IN THE 'NORL' TREE OPTION,  C
C   BUT THE INPUT DECK ASSUMES THAT 'DHFQ' IMPLIES RELATIVISTIC.       C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RKNT(MDM)
C
      COMPLEX*16 DOT
      COMPLEX*16 VIJ(MDM,MDM)
      COMPLEX*16 OLSX(MDM,MDM),OLSY(MDM,MDM),OLSZ(MDM,MDM),
     &           OSLX(MDM,MDM),OSLY(MDM,MDM),OSLZ(MDM,MDM)
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VANM(MDM,MDM),
     &           VSLF(MDM,MDM),VUEH(MDM,MDM),QDIR(MDM,MDM),
     &           QXCH(MDM,MDM),WDIR(MDM,MDM),WXCH(MDM,MDM),CPLE(MDM,MDM)
      COMPLEX*16 CNTC(MDM,MDM),BTHE(MDM,MDM)
      COMPLEX*16 OPRJ(MDM,MDM),ELWM(MDM)
      COMPLEX*16 PINX(MDM,MDM),PINY(MDM,MDM),PINZ(MDM,MDM)
      COMPLEX*16 PMNX(MDM,MDM),PMNY(MDM,MDM),PMNZ(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,MFT),XNUC(MCT,MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/LSHF/SHLEV(4),SHLV,ILEV
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VANM,VSLF,
     &            VUEH,QDIR,QXCH,WDIR,WXCH,CPLE
      COMMON/NTRX/CNTC,BTHE
      COMMON/OTTI/OLSX,OLSY,OLSZ,OSLX,OSLY,OSLZ
C
      DATA DEPS/0.01D0/
C
C     CUTOFF PARAMETER (Zα)^2*mc < CUTK < mc
C
C     MAXIMUM NUCLEAR CHARGE IN MOLECULE (IMPORTANT FOR OTHER CUTK VALS)
      ZMAX = 0.0D0
      DO IZ=1,NCNT
        IF(ZNUC(IZ).GT.ZMAX) THEN
          ZMAX = ZNUC(IZ)
        ENDIF
      ENDDO
      CUTK = ZMAX*ZMAX/CV
C
C     DFNOTE: HARRY SAYS JUST MAKE THE CUTOFF ENERGY M*CV*CV FOR NOW
      CUTK = CV
C
C     AMPLITUDE FOR FREE-WAVE SELF-INTERACTION TERM (B&S 19.3)
      AG2 = (DLOG(EMSS*CV/CUTK)-TWLG+11.0D0/24.0D0)/(3.0D0*CV*PI)
      AG2 =-AG2/(EMSS*EMSS*CV*CV)
C
C     AMPLITUDE FOR LOW-ENERGY CONTRIBUTION
      ALW =-2.0D0/(3.0D0*PI*CV)
C
C**********************************************************************C
C     HIGH-ENERGY (NUCLEAR CONTACT) CONTRIBUTION                       C
C**********************************************************************C
C
      IF(ITER.NE.1.AND.ILEV.LT.4) GOTO 201
C
C     INITIALISE HIGH-ENERGY SELF-INTERACTION MATRIX
      DO I=1,NDIM
        DO J=1,NDIM
          CNTC(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     LOOP OVER ALL NUCLEAR CENTRES
      DO IZ=1,NCNT
C
C       LARGE-LARGE CONTRIBUTION (FACTOR -4π FROM POISSON'S EQUATION)
        CALL VNCOLAP(VIJ,IZ,1,0,1,2)
        DO I=1,NDIM-NSKP
          DO J=1,NDIM-NSKP
            CNTC(I,J) = CNTC(I,J) - AG2*4.0D0*PI*ZNUC(IZ)*VIJ(I,J)
          ENDDO
        ENDDO
C
C       SMALL-SMALL CONTRIBUTION (FACTOR -4π FROM POISSON'S EQUATION)
        CALL VNCOLAP(VIJ,IZ,4,0,1,2)
        DO I=1,NDIM-NSKP
          DO J=1,NDIM-NSKP
            K = I+NSKP
            L = J+NSKP
            CNTC(K,L) = CNTC(K,L) - AG2*4.0D0*PI*ZNUC(IZ)*VIJ(I,J)
          ENDDO
        ENDDO
C
C     END LOOP OVER NUCLEAR CENTRES
      ENDDO
C
201   CONTINUE
C
C**********************************************************************C
C     LOW-ENERGY BETHE FORMULA E(<) = <M|Γ|M> FOR OCCUPIED STATES M.   C
C**********************************************************************C
C
C     GOTO 301
C
C     GENERATE ALL NECESSARY OVERLAP INTEGRALS
C     IF(ITER.EQ.1.OR.ILEV.EQ.4) THEN
        CALL VMOMNT0(OLSX,2,1,1,2)
        CALL VMOMNT0(OLSY,2,2,1,2)
        CALL VMOMNT0(OLSZ,2,3,1,2)
        CALL VMOMNT0(OSLX,3,1,1,2)
        CALL VMOMNT0(OSLY,3,2,1,2)
        CALL VMOMNT0(OSLZ,3,3,1,2)
C     ENDIF
C
C     INITIALISE INTERMEDIATE STORAGE ARRAYS
      DO I=1,NDIM
        DO N=1,NDIM
          PINX(I,N) = DCMPLX(0.0D0,0.0D0)
          PINY(I,N) = DCMPLX(0.0D0,0.0D0)
          PINZ(I,N) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     CONTRACT ON ALL POSITIVE-ENERGY STATES -> PINX,PINY,PINZ
      CALL CPU_TIME(TI)
      DO N=1,NDIM-NSKP
C
C       ADDRESS FOR POSITIVE-ENERGY STATE N
        NP = N+NSKP
        DO I=1,NDIM-NSKP
C
C         CONTRACT OVER THIS INDEX FOR ALL POSITIVE-ENERGY STATES N
          DO J=1,NDIM-NSKP
C
C           COUPLE OLS ELEMENTS WITH SMALL-COMPONENT COEFFICIENTS
            PINX(I     ,N) = PINX(I     ,N) + COEF(J+NSKP,NP)*OLSX(I,J)
            PINY(I     ,N) = PINY(I     ,N) + COEF(J+NSKP,NP)*OLSY(I,J)
            PINZ(I     ,N) = PINZ(I     ,N) + COEF(J+NSKP,NP)*OLSZ(I,J)
C
C           COUPLE OSL ELEMENTS WITH LARGE-COMPONENT COEFFICIENTS
            PINX(I+NSKP,N) = PINX(I+NSKP,N) + COEF(J     ,NP)*OSLX(I,J)
            PINY(I+NSKP,N) = PINY(I+NSKP,N) + COEF(J     ,NP)*OSLY(I,J)
            PINZ(I+NSKP,N) = PINZ(I+NSKP,N) + COEF(J     ,NP)*OSLZ(I,J)
C
          ENDDO
C
        ENDDO
      ENDDO
C
C     INITIALISE INTERMEDIATE STORAGE ARRAYS
      DO M=1,NDIM
        DO N=1,NDIM
          PMNX(M,N) = DCMPLX(0.0D0,0.0D0)
          PMNY(M,N) = DCMPLX(0.0D0,0.0D0)
          PMNZ(M,N) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     CONTRACT ON ALL OCCUPIED-ENERGY STATES -> PMNX,PMNY,PMNZ
      CALL CPU_TIME(TI)
      DO M=1,NOCC
C
C       ADDRESS FOR POSITIVE-ENERGY STATE M
        MP = M+NSKP
        DO N=1,NDIM-NSKP
C
C         CONTRACT OVER THIS INDEX FOR ALL POSITIVE-ENERGY STATES N
          DO I=1,NDIM
            PMNX(M,N) = PMNX(M,N) + DCONJG(COEF(I,MP))*PINX(I,N)
            PMNY(M,N) = PMNY(M,N) + DCONJG(COEF(I,MP))*PINY(I,N)
            PMNZ(M,N) = PMNZ(M,N) + DCONJG(COEF(I,MP))*PINZ(I,N)
          ENDDO
C
        ENDDO
      ENDDO
C
      DO M=1,NDIM
        ELWM(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     CONSTRUCT FINISHED MATRIX ELMENTS FOR EACH OCCUPIED M
      DO M=1,NOCC
C
C       INITIALISE THE FINAL COUNTER
        ELWM(M) = DCMPLX(0.0D0,0.0D0)
C
C       STATE ENERGY FOR M
        EM = EIGN(M+NSKP)
C
C       LOOP OVER POSITIVE-ENERGY STATES N
        DO N=1,NDIM-NSKP
C
C         STATE ENERGY FOR N
          EN = EIGN(N+NSKP)
C
C         ENERGY DIFFERENCE
          DEMN = EM-EN
C
C         MANUALLY SKIP CASES OF DEGENERATE ENERGY LEVELS
          IF(DABS(DEMN).GT.DEPS) THEN
C
C           LOGARITHMIC ENERGY TERM
            DTRM = (CV*CUTK-DEMN)/DABS(DEMN)
            DTRM = DEMN*DLOG(DTRM)
C
C           VECTOR DOT PRODUCT OF TRANSITION CURRENTS
            DOT = PMNX(M,N)*DCONJG(PMNX(M,N))
     &          + PMNY(M,N)*DCONJG(PMNY(M,N))
     &          + PMNZ(M,N)*DCONJG(PMNZ(M,N))
C
C           MATRIX ELEMENT
            ELWM(M) = ELWM(M) + ALW*DOT*DTRM
C
          ENDIF
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     APPLY PROJECTION OPERATOR TO GENERATE VSLF (MEAN-FIELD)          C
C**********************************************************************C
C
C     NEW SET OF PROJECTED OVERLAP MATRICES
      DO I=1,NDIM
        DO M=1,NOCC
          OPRJ(I,M) = 0.0D0
          DO J=1,NDIM
            OPRJ(I,M) = OPRJ(I,M) + COEF(J,M+NSKP)*OVLP(I,J)
          ENDDO
        ENDDO
      ENDDO
C
C     COUPLE THE LOW-ENERGY <M||M>
      DO I=1,NDIM
        DO J=1,NDIM
          BTHE(I,J) = DCMPLX(0.0D0,0.0D0)
          DO M=1,NOCC
            BTHE(I,J) = BTHE(I,J) + DCONJG(OPRJ(I,M))*ELWM(M)*OPRJ(J,M)
          ENDDO
        ENDDO
      ENDDO
C
301   CONTINUE
C
C**********************************************************************C
C     TOTAL MATRIX VSLF                                                C
C**********************************************************************C
C
      DO I=1,NDIM
        DO J=1,NDIM
          VSLF(I,J) = CNTC(I,J) + BTHE(I,J)
        ENDDO
      ENDDO
C
      RETURN
      END
