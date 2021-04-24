      SUBROUTINE PVIOLTNTEMP(NEUT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       PPPPPPP  VV    VV IIII OOOOOO  LL      TTTTTTTT NN    NN       C
C       PP    PP VV    VV  II OO    OO LL         TT    NNN   NN       C
C       PP    PP VV    VV  II OO    OO LL         TT    NNNN  NN       C
C       PP    PP VV    VV  II OO    OO LL         TT    NN NN NN       C
C       PPPPPPP   VV  VV   II OO    OO LL         TT    NN  NNNN       C
C       PP         VVVV    II OO    OO LL         TT    NN   NNN       C
C       PP          VV    IIII OOOOOO  LLLLLLLL   TT    NN    NN       C
C                                                                      C
C -------------------------------------------------------------------- C
C  PVIOLTN PERFORMS A P-ODD EFFECTIVE OPERATOR ANALYSIS.               C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    NEUT: LIST OF NEUTRON NUMBERS FOR EACH OF THE NUCLEAR CENTRES.    C
C -------------------------------------------------------------------- C
C DFNOTE: THE IMAGINARY COMPONENTS OF EORB DO NOT SUM TO ZERO -- A     C
C         RELIC OF INTER-ATOMIC MATRIX ELEMENTS AND THE ASSIGNMENT OF  C
C         PHASE `CONE' TO THE SMALL COMPONENT TO ENSURE REAL ATOM-     C
C         CENTRED RESULTS. FOR DISCUSSIONS ON THIS, SEE:               C
C         Phys. Chem. Chem. Phys., 2011, 13, 864–876                   C
C**********************************************************************C
      INCLUDE 'param.h'
C
      CHARACTER*2 ELMT(120),ELA
      CHARACTER*4 HMLT
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 TITLE
C
      COMPLEX*16 EORB(MCT,MDM)
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 VNUCLS(MDM,MDM),VNUCSL(MDM,MDM)
      COMPLEX*16 VBLOCK(MDM,MDM)
C
      DIMENSION ECNT(MCT),EIOCC(MDM)
      DIMENSION NEUT(MCT),LF(MCT)
      DIMENSION EIRL(MCT,MDM),ERL(MCT,MDM)

      dimension vrl(mdm,mdm)
C
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/MATH/PI,PI12,PI52,PILG,TWLG,TW12
      COMMON/MDLV/ELMT
      COMMON/PHYS/CV,PRTN,GFREE,GFRMI,WEIN,CHZ,CFM
      COMMON/PRMS/HMLT,IOPT,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/QNMS/LABICN(MDM),LABKQN(MDM),LABMQN(MDM)
      COMMON/SPEC/BSET(MBS,MKP,MCT),BXYZ(3,MCT),ZNUC(MCT),AMSS(MCT),
     &            CNUC(MCT),PNUC,LRGE(MCT,MKP,MKP+1),NFNC(MKP,MCT),
     &            KAPA(MKP,MCT),IZNC(MCT),IQNC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSKP,NOCC,NVRT
C
C     THIS IS STRICTLY A RELATIVISTIC PHENOMENON
      IF(HMLT.EQ.'NORL') THEN
        WRITE(6, *) 'In PVIOLTN: not possible for this Hamiltonian!'
        WRITE(7, *) 'In PVIOLTN: not possible for this Hamiltonian!'
      ENDIF
C
      IZ=1
C
C     GENERATE NUCLEAR ATTRACTION OVERLAP INTEGRALS
      CALL VNCOLAPATOM(VNUCLS,IZ,2,0,1,2)
      CALL VNCOLAPATOM(VNUCSL,IZ,3,0,1,2)
C
C     LOOP OVER OCCUPIED ORBITALS
      DO IOCC=1,4
        DO JOCC=31,33
C
C         POSITIVE-SPECTRUM ADDRESS
          ISKP = IOCC+NSKP
          JSKP = JOCC+NSKP
C
C         RESET VBLOCK VALUE
          VBLOCK(IOCC,JOCC) = DCMPLX(0.0D0,0.0D0)
C
C         LOOP OVER BASIS SET
          DO I=1,NDIM-NSKP
            DO J=1,NDIM-NSKP
C
              VBLOCK(IOCC,JOCC) = VBLOCK(IOCC,JOCC)
     &        + DCONJG(COEF(I     ,ISKP))*COEF(J+NSKP,JSKP)*VNUCLS(I,J)
     &        + DCONJG(COEF(I+NSKP,ISKP))*COEF(J     ,JSKP)*VNUCSL(I,J)
C
            ENDDO
          ENDDO
C
        ENDDO
      ENDDO
C
      DO IOCC=1,4
        WRITE(*,*) (DIMAG(VBLOCK(IOCC,JOCC)),JOCC=31,33)
      ENDDO
C
      WRITE(6, *) REPEAT('=',LTB)
      WRITE(7, *) REPEAT('=',LTB)
C
      RETURN
      END
C
C
      SUBROUTINE VNCOLAPATOM(VIJ,IZ,ITT,IQ,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    VV    VV NN    NN  CCCCCC   OOOOOO  LL          AA    PPPPPPP     C
C    VV    VV NNN   NN CC    CC OO    OO LL         AAAA   PP    PP    C
C    VV    VV NNNN  NN CC       OO    OO LL        AA  AA  PP    PP    C
C    VV    VV NN NN NN CC       OO    OO LL       AA    AA PP    PP    C
C     VV  VV  NN  NNNN CC       OO    OO LL       AAAAAAAA PPPPPPP     C
C      VVVV   NN   NNN CC    CC OO    OO LL       AA    AA PP          C
C       VV    NN    NN  CCCCCC   OOOOOO  LLLLLLLL AA    AA PP          C
C                                                                      C
C -------------------------------------------------------------------- C
C  VNCOLAPATOM CONSTRUCTS A MATRIX OF (μ,T|rho_nuc|ν,T') OVERLAP       C
C  INTEGRALS OVER ALL BASIS FUNCTIONS, AND SAVES THE RESULT TO VIJ.    C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    IZ                -> NUCLEAR CENTRE OF INTEREST.                  C
C    ITT   = {1,2,3,4} -> {LL,LS,SL,SS} COMPONENT COMBINATION.         C
C    IQ    = {0,1,2,3} -> {0,X,Y,Z} THE COUPLING σ MATRIX.             C
C    I1,I2 = {1,2,3,4} -> BASIS FUNCTION OVERLAPS FOR EQ-COEFFS.       C
C**********************************************************************C
      INCLUDE 'param.h'
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
      DIMENSION CP(MB2,3),APH(MB2),PNC(MB2)
C
      COMPLEX*16 CONE,TA1,TA2
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 VIJ(MDM,MDM)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
      COMMON/MATH/PI,PI12,PI52,PILG,TWLG,TW12
      COMMON/SPEC/BSET(MBS,MKP,MCT),BXYZ(3,MCT),ZNUC(MCT),AMSS(MCT),
     &            CNUC(MCT),PNUC,LRGE(MCT,MKP,MKP+1),NFNC(MKP,MCT),
     &            KAPA(MKP,MCT),IZNC(MCT),IQNC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSKP,NOCC,NVRT
      COMMON/TMMD/TELL,TESS,TELS,TESL,TRLL,TRSS,TRLS,TRSL,TRBR
C
C     UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     INITIALISE THE OVERLAP ARRAY
      DO I=1,NDIM-NSKP
        DO J=1,NDIM-NSKP
          VIJ(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     LOOP OVER ATOMIC CENTRES A AND B (USE INDEX 1000)                C
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
C     MUST BE ATOM-CENTRED
      IF(ICNTA.NE.IZ.OR.ICNTB.NE.IZ) GOTO 1000
C
C**********************************************************************C
C     LOOP OVER KQNS FOR CENTRES A AND B (USE INDEX 2000)              C
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
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1)+1,ICNTA)
        DO IBAS=1,NBAS(1)
          EXPT(IBAS,1) = BSET(IBAS,LQN(1)+1,ICNTA)
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
        NBAS(2) = NFNC(LQN(2)+1,ICNTB)
        DO JBAS=1,NBAS(2)
          EXPT(JBAS,2) = BSET(JBAS,LQN(2)+1,ICNTB)
        ENDDO
C
C     MUST HAVE KAPPA OVERLAP (-1,1)
      IF(KQN(1).NE.-1.OR.KQN(2).NE.1) GOTO 2000
C
C     NUMBER OF UNIQUE ADDRESSES IN THIS EXPANSION
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C**********************************************************************C
C     LOOP OVER MQNS FOR CENTRES A AND B (USE INDEX 3000)              C
C**********************************************************************C
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
C**********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS (MA,MB).  C
C     FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE ORDERED:             C
C -------------------------------------------------------------------- C
C  11: = (-|MQN(A)|,-|MQN(B)|) -> 1  12: = (-|MQN(A)|,+|MQN(B)|) -> 2  C
C  21: = (+|MQN(A)|,-|MQN(B)|) -> 3  22: = (+|MQN(A)|,+|MQN(B)|) -> 4  C
C**********************************************************************C
C
C     CALCULATE COMPONENT OFFSETS
      NA1 = LRGE(ICNTA,KA,2*MA-1)
      NA2 = LRGE(ICNTA,KA,2*MA  )
      NB1 = LRGE(ICNTB,KB,2*MB-1)
      NB2 = LRGE(ICNTB,KB,2*MB  )
C
C**********************************************************************C
C     GENERATE OVERLAP MATRICES                                        C
C**********************************************************************C
C
C     CALCULATE OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
C
C         EXPONENT COMBINATIONS
          EI  = EXPT(IBAS,1)
          EJ  = EXPT(JBAS,2)
          EIJ = EI+EJ
          EPR = EI*EJ
          E12 = DSQRT(EPR)
          E14 = DSQRT(E12)
          E34 = E14*E14*E14
C
C          EI12 = DSQRT(EI)
C          EI14 = DSQRT(EI12)
C          EI34 = EI14*EI14*EI14
C          RNL  = TW12*EI34/DSQRT(GAMHLF(3))
C
C          EJ12 = DSQRT(EJ)
C          EJ14 = DSQRT(EJ12)
C          EJ34 = EJ14*EJ14*EJ14
C          RNS  = TW12*EJ34/DSQRT(GAMHLF(7))
C
C          EJ3 = 2.0D0*EJ
C          RNL = DSQRT(2.0D0*DSQRT(EI3)/GAMHLF(3))
C          RNS = DSQRT(2.0D0*DSQRT(EJ3)/GAMHLF(7))
C
C         INCLUDE NUCLEAR EXPONENT
          ESM = CNUC(IZ)+EIJ
          EBT = CNUC(IZ)+EI
          EBT = EBT/ESM
          EFC = CNUC(IZ)/ESM
          EMX = CNUC(IZ)/ESM
          EMX = DSQRT(EMX*EMX*EMX)
C
C         PRE-FACTOR
          PRE = 4.0D0*DSQRT(6.0D0/5.0D0)/(PI*PI12)
C
C         FIRST TERM
          IF(ITT.EQ.2) THEN
            TA1 = PRE*E34*EMX*EBT
          ENDIF
C          
C         NEXT TERM
          AX  = E14*EPR
          AY  = 1.0D0/ESM
          AZ  = CNUC(IZ)/ESM
          AZ  = DSQRT(AZ*AZ*AZ)
          PRE = 4.0D0*TW12/(PI*PI12)
          
          IF(ITT.EQ.3) THEN
            TA1 = PRE*AX*AY*AZ
          ENDIF
C
C         MATRIX ELEMENTS
          VIJ(NA1+IBAS,NB1+JBAS) = CONE*TA1
          VIJ(NA2+IBAS,NB2+JBAS) = CONE*TA1
C
        ENDDO
      ENDDO
C
C     END LOOP OVER CENTRES A AND B
3000  CONTINUE
2000  CONTINUE
1000  CONTINUE
C
      RETURN
      END
