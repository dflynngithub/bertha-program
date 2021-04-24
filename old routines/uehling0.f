      SUBROUTINE UEHLING0(UMAT,EXL,ZCRG,KQN,NFUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  UU    UU EEEEEEEE HH    HH LL      IIII NN    NN  GGGGGG   000000   C
C  UU    UU EE       HH    HH LL       II  NNN   NN GG    GG 00    00  C
C  UU    UU EE       HH    HH LL       II  NNNN  NN GG       00    00  C
C  UU    UU EEEEEE   HHHHHHHH LL       II  NN NN NN GG       00    00  C
C  UU    UU EE       HH    HH LL       II  NN  NNNN GG   GGG 00    00  C
C  UU    UU EE       HH    HH LL       II  NN   NNN GG    GG 00    00  C
C   UUUUUU  EEEEEEEE HH    HH LLLLLLL IIII NN    NN  GGGGGG   000000   C
C                                                                      C
C -------------------------------------------------------------------- C
C  UEHLING0 GENERATES ATOMIC UEHLING INTERACTION MATRIX FOR KQNA.      C
C**********************************************************************C
      PARAMETER(MDM=2500,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26)
C
      CHARACTER*4 HMLTN
C
      DIMENSION RN(MBS*MBS,4),UMAT(2*MBS,2*MBS),EXL(MBS)
C
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,IEQS,ILEV
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
      DATA ROOTPI/1.7724538509055160D0/
C
C     DETERMINE THE LQN
      IF(KQN.LT.0) THEN
        LQN =-(KQN+1)
      ELSE
        LQN =  KQN
      ENDIF
      RL = DFLOAT(LQN)
      G  = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NFUN,LQN)
C
C     PRE-FACTOR COMMON TO ALL MATRIX ELEMENTS
      VCF  = DSQRT(PNUC)/ROOTPI
      VCF  = (VCF)**3
      VCF  =-2.0D0*ZCRG*VCF/(3.0D0*CV*CV)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NFUN
        EI = EXL(IBAS)
        DO JBAS=1,NFUN
          M    = M+1
          EJ   = EXL(JBAS)
          EIJ  = EI+EJ
          EPR  = EI*EJ
          T52  = RL+2.5D0
          E52  = EIJ**T52
C
C         LL OVERLAP
          UMAT(IBAS,JBAS) = VCF*RN(M,1)*UEHINT0(LQN+1,EIJ)
C
C         SS OVERLAP
          IF(HMLTN.EQ.'NORL') GOTO 50
C
C         SMALL COMPONENT BLOCK ADDRESSES
          KBAS = IBAS+NFUN
          LBAS = JBAS+NFUN
C
          VSA = 4.0D0*EPR*UEHINT0(LQN+2,EIJ)
          IF(KQN.GT.0) THEN
            VSN =-2.0D0*EIJ*G*UEHINT0(LQN+1,EIJ) 
     &                  + G*G*UEHINT0(LQN  ,EIJ)
          ELSE
            VSN = 0.0D0
          ENDIF
          UMAT(KBAS,LBAS) = VCF*RN(M,3)*(VSA+VSN)
C
50        CONTINUE
C         
        ENDDO
      ENDDO
C      
      RETURN
      END
