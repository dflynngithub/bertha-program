      SUBROUTINE ONEEL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C               OOOOOO  NN    NN EEEEEEE EEEEEEE LL                    C
C              OO    OO NNN   NN EE      EE      LL                    C
C              OO    OO NNNN  NN EE      EE      LL                    C
C              OO    OO NN NN NN EEEEEE  EEEEEE  LL                    C
C              OO    OO NN  NNNN EE      EE      LL                    C
C              OO    OO NN   NNN EE      EE      LL                    C
C               OOOOOO  NN    NN EEEEEEE EEEEEEE LLLLLLL               C
C                                                                      C
C -------------------------------------------------------------------- C
C  ONEEL CONSTRUCTS A FULL SET OF MULTI-CENTRE OVERLAP, KINETIC AND    C
C  NUCLEAR ATTRACTION BASIS FUNCTION MATRIX ELEMENTS.                  C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
      DIMENSION RC(MB2,MRC),CP(MB2,3),APH(MB2),PNC(MB2)
C
      COMPLEX*16 CTMP1,CTMP2,CTMP3,CTMP4
      COMPLEX*16 E11A,E11B,E11C,TRM11,E21A,E21B,E21C,TRM21
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 SLL(MBS,MBS,4),SSS(MBS,MBS,4),
     &           VLL(MBS,MBS,4),VSS(MBS,MBS,4),
     &           TLL(MBS,MBS,4),TLS(MBS,MBS,4),TSL(MBS,MBS,4)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VUEH(MDM,MDM),
     &           VSLF(MDM,MDM),QDIR(MDM,MDM),QXCH(MDM,MDM),
     &           WDIR(MDM,MDM),WXCH(MDM,MDM),CPLE(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,MFT),XNUC(MCT,MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/E0LL/E0LLFL(MFL,4),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,4),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VUEH,VSLF,
     &            QDIR,QXCH,WDIR,WXCH,CPLE
      COMMON/TMMD/TELL,TESS,TELS,TESL,TRLL,TRSS,TRLS,TRSL,TRBR
C
C     INITIALISE STORAGE MATRICES
      DO I=1,NDIM
        DO J=1,NDIM
          HNUC(I,J) = DCMPLX(0.0D0,0.0D0)
          HKIN(I,J) = DCMPLX(0.0D0,0.0D0)
          OVLP(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)      C
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
C     LOOP OVER CENTRE B
      DO 2000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
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
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 2000 MA=1,IABS(KQN(1))
        MJA    = 2*MA-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 2000 MB=1,IABS(KQN(2))
        MJB    = 2*MB-1
        MQN(2) = MJB
C
C**********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS (MA,MB).  C
C     FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE ORDERED:             C
C -------------------------------------------------------------------- C
C  11: = (-|MQN(A)|,-|MQN(B)|) -> 1  12: = (-|MQN(A)|,+|MQN(B)|) -> 2  C
C  21: = (+|MQN(A)|,-|MQN(B)|) -> 3  22: = (+|MQN(A)|,+|MQN(B)|) -> 4  C
C**********************************************************************C
C
C     CONSTRUCTION OF ONE-ELECTRON MATRICES BY TT' BLOCKS...
C
C     THIS PHASE RELATES EQ22 AND EQ12 COEFFS TO EQ11 AND EQ21
      PHS = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
     &                    *DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXM = NBAS(1)*NBAS(2)
C
C     INITIALISE STORAGE ARRAYS
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          DO IM=1,4
            SLL(IBAS,JBAS,IM) = DCMPLX(0.0D0,0.0D0)
            SSS(IBAS,JBAS,IM) = DCMPLX(0.0D0,0.0D0)
            TLL(IBAS,JBAS,IM) = DCMPLX(0.0D0,0.0D0)
            TLS(IBAS,JBAS,IM) = DCMPLX(0.0D0,0.0D0)
            TSL(IBAS,JBAS,IM) = DCMPLX(0.0D0,0.0D0)
            VLL(IBAS,JBAS,IM) = DCMPLX(0.0D0,0.0D0)
            VSS(IBAS,JBAS,IM) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     PART 1: THE LL MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMLL  = LQN(1)+LQN(2)
      NTUVLL = (LAMLL+1)*(LAMLL+2)*(LAMLL+3)/6
C
C     GENERATE ELL0 COEFFICIENTS (IPHS = +1)
      CALL CPU_TIME(TDM1)
      IF(EQFILE) THEN
        DO ITUV=1,NTUVLL
          IAD = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXM
          DO M=1,MAXM
            E11(M,ITUV) = DCMPLX(E0LLFL(IAD+M,1),E0LLFL(IAD+M,2))
            E21(M,ITUV) = DCMPLX(E0LLFL(IAD+M,3),E0LLFL(IAD+M,4))
          ENDDO
        ENDDO
      ELSE
        CALL EQLLMK(E11,E21,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,0)
      ENDIF
      CALL CPU_TIME(TDM2)
      TELL = TELL+TDM2-TDM1
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
          EIJ   = EXL(IBAS,1)+EXL(JBAS,2)
          EROOT = DSQRT(PI/EIJ)**3
          SLL(IBAS,JBAS,1) = EROOT*E11(M,1)
          SLL(IBAS,JBAS,3) = EROOT*E21(M,1)
          SLL(IBAS,JBAS,2) =-PHS*DCONJG(SLL(IBAS,JBAS,3))
          SLL(IBAS,JBAS,4) = PHS*DCONJG(SLL(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     NUCLEAR ATTRACTION MATRIX ELEMENTS
      DO IZ=1,NCNT
C
C       NUCLEAR COORDINATES
        CX = BXYZ(1,IZ)
        CY = BXYZ(2,IZ)
        CZ = BXYZ(3,IZ)
C
C       GAUSSIAN PRODUCT THEOREM OVER BASIS FUNCTIONS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M   = M+1
            EIJ = EXL(IBAS,1)+EXL(JBAS,2)
            PX  = (XYZ(1,1)*EXL(IBAS,1) + XYZ(1,2)*EXL(JBAS,2))/EIJ
            PY  = (XYZ(2,1)*EXL(IBAS,1) + XYZ(2,2)*EXL(JBAS,2))/EIJ
            PZ  = (XYZ(3,1)*EXL(IBAS,1) + XYZ(3,2)*EXL(JBAS,2))/EIJ
            CP(M,1) = CX-PX
            CP(M,2) = CY-PY
            CP(M,3) = CZ-PZ
          ENDDO
        ENDDO
C
        IF(NNUC(IZ).EQ.0) THEN
C       POINT-NUCLEUS APPROXIMATION
C
C         PREPARE ELEMENTS FOR RMAKE
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
C
C             POINT-NUCLEUS EFFECTIVE PARAMETERS
              EIJ    = EXL(IBAS,1)+EXL(JBAS,2)
              APH(M) = EIJ
              PNC(M) = 2.0D0*ZNUC(IZ)*PI/EIJ
C
            ENDDO
          ENDDO
C
C         GENERATE A BATCH OF R-INTEGRALS
          CALL CPU_TIME(TDM1)
          CALL RMAKE(RC,CP,APH,MAXM,LAMLL)
          CALL CPU_TIME(TDM2)
          TRLL = TRLL + TDM2 - TDM1
C
C         NUCLEAR ATTRACTION INTEGRALS AS A FINITE SUM OF ELL0 AND RC
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              DO ITUV=1,NTUVLL
                VLL(IBAS,JBAS,1) = VLL(IBAS,JBAS,1)
     &                           - PNC(M)*E11(M,ITUV)*RC(M,ITUV)
                VLL(IBAS,JBAS,3) = VLL(IBAS,JBAS,3)
     &                           - PNC(M)*E21(M,ITUV)*RC(M,ITUV)
              ENDDO
              VLL(IBAS,JBAS,2) =-PHS*DCONJG(VLL(IBAS,JBAS,3))
              VLL(IBAS,JBAS,4) = PHS*DCONJG(VLL(IBAS,JBAS,1))
            ENDDO
          ENDDO
C
        ELSE
C       BEST-FIT EXPANSION
C
C         LOOP OVER GAUSSIANS FOR THIS NUCLEAR CENTRE
          DO IFT=1,NNUC(IZ)
C
C           GAUSSIAN EXPONENT AND FRACTIONAL CONTRIBUTION
            XI = XNUC(IZ,IFT)
            FC = FNUC(IZ,IFT)
C
C           PREPARE ELEMENTS FOR RMAKE
            M = 0
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                M = M+1
C
C               FINITE-NUCLEUS BOYS EXPONENT AND MULTIPLIER
                EIJ    = EXL(IBAS,1)+EXL(JBAS,2)
                ESM    = EIJ+XI
                APH(M) = EIJ*XI/ESM
                PNC(M) = 2.0D0*PI*ZNUC(IZ)*FC*DSQRT(XI/ESM)/EIJ
C
              ENDDO
            ENDDO
C
C           GENERATE A BATCH OF R-INTEGRALS
            CALL CPU_TIME(TDM1)
            CALL RMAKE(RC,CP,APH,MAXM,LAMLL)
            CALL CPU_TIME(TDM2)
            TRLL = TRLL + TDM2 - TDM1
C
C           NUCLEAR ATTRACTION INTEGRALS AS A FINITE SUM OF ELL0 AND RC
            M = 0
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                M = M+1
                DO ITUV=1,NTUVLL
                  VLL(IBAS,JBAS,1) = VLL(IBAS,JBAS,1)
     &                             - PNC(M)*E11(M,ITUV)*RC(M,ITUV)
                  VLL(IBAS,JBAS,3) = VLL(IBAS,JBAS,3)
     &                             - PNC(M)*E21(M,ITUV)*RC(M,ITUV)
                ENDDO
                VLL(IBAS,JBAS,2) =-PHS*DCONJG(VLL(IBAS,JBAS,3))
                VLL(IBAS,JBAS,4) = PHS*DCONJG(VLL(IBAS,JBAS,1))
              ENDDO
            ENDDO
C
C         END LOOP OVER NUCLEAR BASIS FOR IZ
          ENDDO
C
        ENDIF
C
C     END LOOP OVER CENTRES IZ
      ENDDO
C
C     CONSTRUCT NON-RELATIVISTIC KINETIC ENERGY INTEGRALS (IOS 91)
      IF(HMLT.EQ.'NORL') THEN
        RL2 = DFLOAT(2*LQN(2)+3)
        M   = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M   = M+1
C
C           GAUSSIAN PRODUCT THEOREM DETAILS
            EJ  = EXL(JBAS,2)
            EIJ = EXL(IBAS,1) + EXL(JBAS,2)
            ERT = DSQRT(PI/EIJ)**3
            PX = (XYZ(1,1)*EXL(IBAS,1) + XYZ(1,2)*EXL(JBAS,2))/EIJ
            PY = (XYZ(2,1)*EXL(IBAS,1) + XYZ(2,2)*EXL(JBAS,2))/EIJ
            PZ = (XYZ(3,1)*EXL(IBAS,1) + XYZ(3,2)*EXL(JBAS,2))/EIJ
C
C           COORDINATES RELATIVE TO NUCLEAR CENTRE
            PBX = PX-XYZ(1,2)
            PBY = PY-XYZ(2,2)
            PBZ = PZ-XYZ(3,2)
            PB2 = PBX*PBX + PBY*PBY + PBZ*PBZ
C
            E0FC = EJ*RL2 - 2.0D0*EJ*EJ*PB2 - 3.0D0*EJ*EJ/EIJ
            E1FC = 4.0D0*EJ*EJ
C
C           TRUNCATE EXPRESSION DEPENDING ON LAMLL VALUE
C           ALL COMBINATIONS ALLOW FOR THE LAMLL = 0 MANIFOLD
            TRM11 = E0FC*E11(M,IABC(0,0,0))
            TRM21 = E0FC*E21(M,IABC(0,0,0))
C           IF LAMLL > 0 PROVIDE SECOND BUNCH OF TERMS
            IF(LAMLL.GE.1) THEN
              E11A = E11(M,IABC(1,0,0))
              E21A = E21(M,IABC(1,0,0))
              E11B = E11(M,IABC(0,1,0))
              E21B = E21(M,IABC(0,1,0))
              E11C = E11(M,IABC(0,0,1))
              E21C = E21(M,IABC(0,0,1))
              TRM11 = TRM11 - E1FC*(PBX*E11A + PBY*E11B + PBZ*E11C)
              TRM21 = TRM21 - E1FC*(PBX*E21A + PBY*E21B + PBZ*E21C)
            ENDIF
C           IF LAMLL > 1 PROVIDE FINAL BUNCH OF TERMS
            IF(LAMLL.GE.2) THEN
              E11A = E11(M,IABC(2,0,0))
              E21A = E21(M,IABC(2,0,0))
              E11B = E11(M,IABC(0,2,0))
              E21B = E21(M,IABC(0,2,0))
              E11C = E11(M,IABC(0,0,2))
              E21C = E21(M,IABC(0,0,2))
              TRM11 = TRM11 - E1FC*(E11A + E11B + E11C)
              TRM21 = TRM21 - E1FC*(E21A + E21B + E21C)
            ENDIF
            TLL(IBAS,JBAS,1) = ERT*TRM11
            TLL(IBAS,JBAS,3) = ERT*TRM21
            TLL(IBAS,JBAS,2) =-PHS*DCONJG(TLL(IBAS,JBAS,3))
            TLL(IBAS,JBAS,4) = PHS*DCONJG(TLL(IBAS,JBAS,1))
C            
          ENDDO
        ENDDO
C       NON-RELATIVISTIC HAMILTONIAN MATRICES COMPLETE
        GOTO 500
      ENDIF
C
C**********************************************************************C
C     PART 2: THE SS MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMSS  = LQN(1)+LQN(2)+2
      NTUVSS = (LAMSS+1)*(LAMSS+2)*(LAMSS+3)/6
C
C     GENERATE ESS0 COEFFICIENTS (IPHS = +1)
      CALL CPU_TIME(TDM1)
      IF(EQFILE) THEN
        DO ITUV=1,NTUVSS
          IAD = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXM
          DO M=1,MAXM
            E11(M,ITUV) = DCMPLX(E0SSFL(IAD+M,1),E0SSFL(IAD+M,2))
            E21(M,ITUV) = DCMPLX(E0SSFL(IAD+M,3),E0SSFL(IAD+M,4))
          ENDDO
        ENDDO
      ELSE
        CALL EQSSMK(E11,E21,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,0)
      ENDIF
      CALL CPU_TIME(TDM2)
      TESS = TESS+TDM2-TDM1
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
          EIJ   = EXL(IBAS,1)+EXL(JBAS,2)
          EROOT = DSQRT(PI/EIJ)**3
          SSS(IBAS,JBAS,1) = EROOT*E11(M,1)
          SSS(IBAS,JBAS,3) = EROOT*E21(M,1)
          SSS(IBAS,JBAS,2) =-PHS*DCONJG(SSS(IBAS,JBAS,3))
          SSS(IBAS,JBAS,4) = PHS*DCONJG(SSS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     NUCLEAR ATTRACTION MATRIX ELEMENTS
      DO IZ=1,NCNT
C
C       NUCLEAR COORDINATES
        CX = BXYZ(1,IZ)
        CY = BXYZ(2,IZ)
        CZ = BXYZ(3,IZ)
C
C       GAUSSIAN PRODUCT THEOREM OVER BASIS FUNCTIONS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M   = M+1
            EIJ = EXL(IBAS,1)+EXL(JBAS,2)
            PX  = (XYZ(1,1)*EXL(IBAS,1) + XYZ(1,2)*EXL(JBAS,2))/EIJ
            PY  = (XYZ(2,1)*EXL(IBAS,1) + XYZ(2,2)*EXL(JBAS,2))/EIJ
            PZ  = (XYZ(3,1)*EXL(IBAS,1) + XYZ(3,2)*EXL(JBAS,2))/EIJ
            CP(M,1) = CX-PX
            CP(M,2) = CY-PY
            CP(M,3) = CZ-PZ
          ENDDO
        ENDDO
C
        IF(NNUC(IZ).EQ.0) THEN
C       POINT-NUCLEUS APPROXIMATION
C
C         PREPARE ELEMENTS FOR RMAKE
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
C
C             POINT-NUCLEUS EFFECTIVE PARAMETERS
              EIJ = EXL(IBAS,1)+EXL(JBAS,2)
              APH(M) = EIJ
              PNC(M) = 2.0D0*ZNUC(IZ)*PI/EIJ
C
            ENDDO
          ENDDO
C
C         GENERATE A BATCH OF R-INTEGRALS
          CALL CPU_TIME(TDM1)
          CALL RMAKE(RC,CP,APH,MAXM,LAMSS)
          CALL CPU_TIME(TDM2)
          TRSS = TRSS + TDM2 - TDM1
C
C         NUCLEAR ATTRACTION INTEGRALS AS A FINITE SUM OF ESS0 AND RC
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              DO ITUV=1,NTUVSS
                VSS(IBAS,JBAS,1) = VSS(IBAS,JBAS,1)
     &                           + PNC(M)*E11(M,ITUV)*RC(M,ITUV)
                VSS(IBAS,JBAS,3) = VSS(IBAS,JBAS,3)
     &                           + PNC(M)*E21(M,ITUV)*RC(M,ITUV)
              ENDDO
              VSS(IBAS,JBAS,2) =-PHS*DCONJG(VSS(IBAS,JBAS,3))
              VSS(IBAS,JBAS,4) = PHS*DCONJG(VSS(IBAS,JBAS,1))
            ENDDO
          ENDDO
C
        ELSE
C       BEST-FIT EXPANSION
C
C         LOOP OVER GAUSSIANS FOR THIS NUCLEAR CENTRE
          DO IFT=1,NNUC(IZ)
C
C           GAUSSIAN EXPONENT AND FRACTIONAL CONTRIBUTION
            XI = XNUC(IZ,IFT)
            FC = FNUC(IZ,IFT)
C
C           PREPARE ELEMENTS FOR RMAKE
            M = 0
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                M = M+1
C
C               FINITE-NUCLEUS BOYS EXPONENT AND MULTIPLIER
                EIJ    = EXL(IBAS,1)+EXL(JBAS,2)
                ESM    = EIJ+XI
                APH(M) = EIJ*XI/ESM
                PNC(M) = 2.0D0*PI*ZNUC(IZ)*FC*DSQRT(XI/ESM)/EIJ
C
              ENDDO
            ENDDO
C
C           GENERATE A BATCH OF R-INTEGRALS
            CALL CPU_TIME(TDM1)
            CALL RMAKE(RC,CP,APH,MAXM,LAMSS)
            CALL CPU_TIME(TDM2)
            TRSS = TRSS + TDM2 - TDM1
C
C           NUCLEAR ATTRACTION INTEGRALS AS A FINITE SUM OF ESS0 AND RC
            M = 0
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                M = M+1
                DO ITUV=1,NTUVSS
                  VSS(IBAS,JBAS,1) = VSS(IBAS,JBAS,1)
     &                             + PNC(M)*E11(M,ITUV)*RC(M,ITUV)
                  VSS(IBAS,JBAS,3) = VSS(IBAS,JBAS,3)
     &                             + PNC(M)*E21(M,ITUV)*RC(M,ITUV)
                ENDDO
                VSS(IBAS,JBAS,2) =-PHS*DCONJG(VSS(IBAS,JBAS,3))
                VSS(IBAS,JBAS,4) = PHS*DCONJG(VSS(IBAS,JBAS,1))
              ENDDO
            ENDDO
C
C         END LOOP OVER NUCLEAR BASIS FOR IZ
          ENDDO
C
        ENDIF
C
C     END LOOP OVER CENTRES IZ
      ENDDO
C
C     SUBTRACT THE SS OVERLAP MATRIX AND FINISH CONSTRUCTION
      CV2 = 2.0D0*EMSS*CV*CV
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          VSS(IBAS,JBAS,1) =-VSS(IBAS,JBAS,1) - CV2*SSS(IBAS,JBAS,1)
          VSS(IBAS,JBAS,3) =-VSS(IBAS,JBAS,3) - CV2*SSS(IBAS,JBAS,3)
          VSS(IBAS,JBAS,2) =-PHS*DCONJG(VSS(IBAS,JBAS,3))
          VSS(IBAS,JBAS,4) = PHS*DCONJG(VSS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C**********************************************************************C
C     PART 3: THE SL MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMSS  = LQN(1)+LQN(2)+2
      NTUVSS = (LAMSS+1)*(LAMSS+2)*(LAMSS+3)/6
C
C     KINETIC MATRIX ELEMENTS
      FACT = CV*DSQRT(DFLOAT(2*LQN(2)+3))
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
          EIJ    = EXL(IBAS,1) + EXL(JBAS,2)
          EJRT   = FACT*DSQRT(EXL(JBAS,2))
          EROOT  = DSQRT(PI/EIJ)**3
          TSL(IBAS,JBAS,1) = EJRT*EROOT*E11(M,1)
          TSL(IBAS,JBAS,3) = EJRT*EROOT*E21(M,1)
          TSL(IBAS,JBAS,2) =-PHS*DCONJG(TSL(IBAS,JBAS,3))
          TSL(IBAS,JBAS,4) = PHS*DCONJG(TSL(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C
C**********************************************************************C
C     PART 4: THE LS MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMSS  = LQN(1)+LQN(2)+2
      NTUVSS = (LAMSS+1)*(LAMSS+2)*(LAMSS+3)/6
C
C     GENERATE ESS0 COEFFICIENTS (IPHS = +1)
C
      CALL CPU_TIME(TDM1)
C
C     TODO: THE INDEX SWAP '2 1' MAKES FILE IMPORT ANNOYING. I'M LAZY.
C             JUST GENERATE THESE EQ'S AS A BATCH (DOESN'T TAKE LONG).
C
C     IF(EQFILE) THEN
C       DO ITUV=1,NTUVSS
C         IAD = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXM
C         M = 0
C         DO IBAS=1,NBAS(1)
C           DO JBAS=1,NBAS(2)
C             M = M+1
C             N = (JBAS-1)*NBAS(2) + IBAS
C             PI^{SL}_{IJ} = PI^{LS}_{JI}* (THE NEXT LINES ARE WRONG)
C             E11(M,ITUV) = DCMPLX(E0SSFL(IAD+N,1),-E0SSFL(IAD+N,2))
C             E21(M,ITUV) = DCMPLX(E0SSFL(IAD+N,3),-E0SSFL(IAD+N,4))
C           ENDDO
C         ENDDO
C       ENDDO
C     ELSE
        CALL EQSSMK(E11,E21,EXL,XYZ,KQN,MQN,NBAS,+1,2,1,0)
C     ENDIF
      CALL CPU_TIME(TDM2)
      TESS = TESS+TDM2-TDM1
C
C     KINETIC MATRIX ELEMENTS
      FACT = CV*DSQRT(DFLOAT(2*LQN(1)+3))
      M = 0
      DO JBAS=1,NBAS(2)
        DO IBAS=1,NBAS(1)
          M = M+1
          EIJ    = EXL(JBAS,2) + EXL(IBAS,1)
          EIRT   = FACT*DSQRT(EXL(IBAS,1))
          EROOT  = DSQRT(PI/EIJ)**3
          TLS(IBAS,JBAS,1) = EIRT*EROOT*E11(M,1)
          TLS(IBAS,JBAS,3) = EIRT*EROOT*E21(M,1)
          TLS(IBAS,JBAS,2) =-PHS*DCONJG(TLS(IBAS,JBAS,3))
          TLS(IBAS,JBAS,4) = PHS*DCONJG(TLS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     GENERATE LS MATRICES FROM THE ABOVE SL MATRICES
      M = 0
      DO JBAS=1,NBAS(2)
        DO IBAS=1,NBAS(1)
          M = M + 1
          CTMP1 = DCONJG(TLS(IBAS,JBAS,1))
          CTMP2 = DCONJG(TLS(IBAS,JBAS,2))
          CTMP3 = DCONJG(TLS(IBAS,JBAS,3))
          CTMP4 = DCONJG(TLS(IBAS,JBAS,4))
C
          TLS(IBAS,JBAS,1) = CTMP1
          TLS(IBAS,JBAS,2) = CTMP3
          TLS(IBAS,JBAS,3) = CTMP2
          TLS(IBAS,JBAS,4) = CTMP4
        ENDDO
      ENDDO
C
500   CONTINUE
C
C**********************************************************************C
C     WE NOW HAVE ALL PIECES OF HNUC AND HKIN FOR THIS BLOCK.          C
C**********************************************************************C
C
C     CALCULATE COMPONENT OFFSETS
      IL1 = LRGE(ICNTA,KA,MJA  )
      IL2 = LRGE(ICNTA,KA,MJA+1)
      JL1 = LRGE(ICNTB,KB,MJB  )
      JL2 = LRGE(ICNTB,KB,MJB+1)
C
      IS1 = IL1+NSKP
      IS2 = IL2+NSKP
      JS1 = JL1+NSKP
      JS2 = JL2+NSKP
C
C     LL OVERLAP BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=1,NBAS(1)
            OVLP(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            OVLP(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            OVLP(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            OVLP(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
C
            OVLP(JL1+JBAS,IL1+IBAS) = DCONJG(OVLP(IL1+IBAS,JL1+JBAS))
            OVLP(JL2+JBAS,IL1+IBAS) = DCONJG(OVLP(IL1+IBAS,JL2+JBAS))
            OVLP(JL1+JBAS,IL2+IBAS) = DCONJG(OVLP(IL2+IBAS,JL1+JBAS))
            OVLP(JL2+JBAS,IL2+IBAS) = DCONJG(OVLP(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=JBAS,NBAS(1)
            OVLP(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            OVLP(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            OVLP(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            OVLP(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
C
            OVLP(JL1+JBAS,IL1+IBAS) = DCONJG(OVLP(IL1+IBAS,JL1+JBAS))
            OVLP(JL2+JBAS,IL1+IBAS) = DCONJG(OVLP(IL1+IBAS,JL2+JBAS))
            OVLP(JL1+JBAS,IL2+IBAS) = DCONJG(OVLP(IL2+IBAS,JL1+JBAS))
            OVLP(JL2+JBAS,IL2+IBAS) = DCONJG(OVLP(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     LL NUCLEAR POTENTIAL BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=1,NBAS(1)
            HNUC(IL1+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,1)
            HNUC(IL1+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,2)
            HNUC(IL2+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,3)
            HNUC(IL2+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,4)
C
            HNUC(JL1+JBAS,IL1+IBAS) = DCONJG(HNUC(IL1+IBAS,JL1+JBAS))
            HNUC(JL2+JBAS,IL1+IBAS) = DCONJG(HNUC(IL1+IBAS,JL2+JBAS))
            HNUC(JL1+JBAS,IL2+IBAS) = DCONJG(HNUC(IL2+IBAS,JL1+JBAS))
            HNUC(JL2+JBAS,IL2+IBAS) = DCONJG(HNUC(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=JBAS,NBAS(1)
            HNUC(IL1+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,1)
            HNUC(IL1+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,2)
            HNUC(IL2+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,3)
            HNUC(IL2+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,4)
C
            HNUC(JL1+JBAS,IL1+IBAS) = DCONJG(HNUC(IL1+IBAS,JL1+JBAS))
            HNUC(JL2+JBAS,IL1+IBAS) = DCONJG(HNUC(IL1+IBAS,JL2+JBAS))
            HNUC(JL1+JBAS,IL2+IBAS) = DCONJG(HNUC(IL2+IBAS,JL1+JBAS))
            HNUC(JL2+JBAS,IL2+IBAS) = DCONJG(HNUC(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     NON-RELATIVISTIC HAMILTONIAN HAS A KINETIC MATRIX IN THE LL BLOCK
      IF(HMLT.EQ.'NORL') THEN
C
C       LL KINETIC BLOCKS
        IF(IL1.GT.JL1) THEN
          DO JBAS=1,NBAS(2)
            DO IBAS=1,NBAS(1)
              HKIN(IL1+IBAS,JL1+JBAS) = TLL(IBAS,JBAS,1)
              HKIN(IL1+IBAS,JL2+JBAS) = TLL(IBAS,JBAS,2)
              HKIN(IL2+IBAS,JL1+JBAS) = TLL(IBAS,JBAS,3)
              HKIN(IL2+IBAS,JL2+JBAS) = TLL(IBAS,JBAS,4)
C
              HKIN(JL1+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JL1+JBAS))
              HKIN(JL2+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JL2+JBAS))
              HKIN(JL1+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JL1+JBAS))
              HKIN(JL2+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JL2+JBAS))
            ENDDO
          ENDDO
        ENDIF
C
        IF(IL1.EQ.JL1) THEN
          DO JBAS=1,NBAS(2)
            DO IBAS=JBAS,NBAS(1)
              HKIN(IL1+IBAS,JL1+JBAS) = TLL(IBAS,JBAS,1)
              HKIN(IL1+IBAS,JL2+JBAS) = TLL(IBAS,JBAS,2)
              HKIN(IL2+IBAS,JL1+JBAS) = TLL(IBAS,JBAS,3)
              HKIN(IL2+IBAS,JL2+JBAS) = TLL(IBAS,JBAS,4)
C
              HKIN(JL1+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JL1+JBAS))
              HKIN(JL2+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JL2+JBAS))
              HKIN(JL1+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JL1+JBAS))
              HKIN(JL2+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JL2+JBAS))
            ENDDO
          ENDDO
        ENDIF
C
C       NON-RELATIVISTIC MATRIX CONSTRUCTION COMPLETE
        GOTO 600
C
      ENDIF
C
C     SS OVERLAP BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=1,NBAS(1)
            OVLP(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            OVLP(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            OVLP(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            OVLP(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
C
            OVLP(JS1+JBAS,IS1+IBAS) = DCONJG(OVLP(IS1+IBAS,JS1+JBAS))
            OVLP(JS2+JBAS,IS1+IBAS) = DCONJG(OVLP(IS1+IBAS,JS2+JBAS))
            OVLP(JS1+JBAS,IS2+IBAS) = DCONJG(OVLP(IS2+IBAS,JS1+JBAS))
            OVLP(JS2+JBAS,IS2+IBAS) = DCONJG(OVLP(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=JBAS,NBAS(1)
            OVLP(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            OVLP(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            OVLP(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            OVLP(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
C
            OVLP(JS1+JBAS,IS1+IBAS) = DCONJG(OVLP(IS1+IBAS,JS1+JBAS))
            OVLP(JS2+JBAS,IS1+IBAS) = DCONJG(OVLP(IS1+IBAS,JS2+JBAS))
            OVLP(JS1+JBAS,IS2+IBAS) = DCONJG(OVLP(IS2+IBAS,JS1+JBAS))
            OVLP(JS2+JBAS,IS2+IBAS) = DCONJG(OVLP(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     SS NUCLEAR POTENTIAL BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=1,NBAS(1)
            HNUC(IS1+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,1)
            HNUC(IS1+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,2)
            HNUC(IS2+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,3)
            HNUC(IS2+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,4)
C
            HNUC(JS1+JBAS,IS1+IBAS) = DCONJG(HNUC(IS1+IBAS,JS1+JBAS))
            HNUC(JS2+JBAS,IS1+IBAS) = DCONJG(HNUC(IS1+IBAS,JS2+JBAS))
            HNUC(JS1+JBAS,IS2+IBAS) = DCONJG(HNUC(IS2+IBAS,JS1+JBAS))
            HNUC(JS2+JBAS,IS2+IBAS) = DCONJG(HNUC(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=JBAS,NBAS(1)
            HNUC(IS1+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,1)
            HNUC(IS1+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,2)
            HNUC(IS2+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,3)
            HNUC(IS2+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,4)
C
            HNUC(JS1+JBAS,IS1+IBAS) = DCONJG(HNUC(IS1+IBAS,JS1+JBAS))
            HNUC(JS2+JBAS,IS1+IBAS) = DCONJG(HNUC(IS1+IBAS,JS2+JBAS))
            HNUC(JS1+JBAS,IS2+IBAS) = DCONJG(HNUC(IS2+IBAS,JS1+JBAS))
            HNUC(JS2+JBAS,IS2+IBAS) = DCONJG(HNUC(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     LS BLOCKS
      IF(IL1.GE.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=1,NBAS(1)
            HKIN(IL1+IBAS,JS1+JBAS) = TLS(IBAS,JBAS,1)
            HKIN(IL1+IBAS,JS2+JBAS) = TLS(IBAS,JBAS,2)
            HKIN(IL2+IBAS,JS1+JBAS) = TLS(IBAS,JBAS,3)
            HKIN(IL2+IBAS,JS2+JBAS) = TLS(IBAS,JBAS,4)
C
            HKIN(JS1+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JS1+JBAS))
            HKIN(JS2+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JS2+JBAS))
            HKIN(JS1+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JS1+JBAS))
            HKIN(JS2+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     SL BLOCKS
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=1,NBAS(1)
            HKIN(IS1+IBAS,JL1+JBAS) = TSL(IBAS,JBAS,1)
            HKIN(IS1+IBAS,JL2+JBAS) = TSL(IBAS,JBAS,2)
            HKIN(IS2+IBAS,JL1+JBAS) = TSL(IBAS,JBAS,3)
            HKIN(IS2+IBAS,JL2+JBAS) = TSL(IBAS,JBAS,4)
C
            HKIN(JL1+JBAS,IS1+IBAS) = DCONJG(HKIN(IS1+IBAS,JL1+JBAS))
            HKIN(JL2+JBAS,IS1+IBAS) = DCONJG(HKIN(IS1+IBAS,JL2+JBAS))
            HKIN(JL1+JBAS,IS2+IBAS) = DCONJG(HKIN(IS2+IBAS,JL1+JBAS))
            HKIN(JL2+JBAS,IS2+IBAS) = DCONJG(HKIN(IS2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
600   CONTINUE
C
C     END LOOPS OVER BASIS PAIRS A,B
2000  CONTINUE
C
      RETURN
      END

