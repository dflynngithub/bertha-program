      SUBROUTINE MNPOLE(VIJ,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       MM       MM NN    NN PPPPPPP   OOOOOO  LL       EEEEEEEE       C
C       MMM     MMM NNN   NN PP    PP OO    OO LL       EE             C
C       MMMM   MMMM NNNN  NN PP    PP OO    OO LL       EE             C
C       MM MM MM MM NN NN NN PP    PP OO    OO LL       EEEEEE         C
C       MM  MMM  MM NN  NNNN PPPPPPP  OO    OO LL       EE             C
C       MM   M   MM NN   NNN PP       OO    OO LL       EE             C
C       MM       MM NN    NN PP        OOOOOO  LLLLLLLL EEEEEEEE       C
C                                                                      C
C -------------------------------------------------------------------- C
C  MNPOLE CONSTRUCTS A MATRIX OF (u,T|SIG_Q|v,T') OVERLAP INTEGRALS    C
C  OVER ALL BASIS FUNCTION PAIRS, AND SAVES THE RESULT TO VIJ.         C
C**********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=MKP,MBS=26,MB2=MBS*MBS,
     &          MLL=MKP*(MKP+1)*(MKP+2)/6,IL4=2*(MKP-1),
     &          MRC=(IL4+1)*(IL4+2)*(IL4+3)/6)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 CONE
      COMPLEX*16 VIJ(MDM,MDM)
      COMPLEX*16 E11(MB2,MLL),E21(MB2,MLL)
      COMPLEX*16 SLL(MBS,MBS,4),SSS(MBS,MBS,4)
C
      DIMENSION RC(MB2,MRC),EXPT(MBS,4),XYZ(3,4),APH(MB2),CP(MB2,3),
     &          PNC(MB2),KQN(4),LQN(4),MQN(4),NFUNS(4)
      DIMENSION INDEX(MCT,-MKP:MKP,2*(MKP+1)*MKP)
C
      COMMON/ACSS/INABCD(0:IL4,0:IL4,0:IL4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/PRMS/CV,HMLTN,ITER,IIL,ITREE,IEQS
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR,IOCCM0
      COMMON/TIME/TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TSTRT,TATOM
C
      DATA PI/3.1415926535897932D0/
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     INITIALISE THE OVERLAP ARRAY
      DO I=1,NDIM
        DO J=1,NDIM
          VIJ(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTERS A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     LOOP OVER CENTER A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 0
        ENDIF
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
         LQN(1) = KQN(1)
        ELSE
         LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNA    = NFUNCT(LQN(1)+1,ICNTA)
        NFUNS(1) = NFUNA
C
        DO IBAS=1,NFUNA
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2)=KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNB    = NFUNCT(LQN(2)+1,ICNTB)
        NFUNS(2) = NFUNB
C
        DO IBAS=1,NFUNB
          EXPT(IBAS,2) = EXPSET(IBAS,LQN(2)+1,ICNTB)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 2000 MA=1,IABS(KQN(1))
        MJA    = (2*MA)-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 2000 MB=1,IABS(KQN(2))
        MJB    = (2*MB)-1
        MQN(2) = MJB
C
C**********************************************************************C
C     IMPLEMENT ANY AVAILABLE SELECTION RULES HERE                     C
C**********************************************************************C
C
C     SELECTION RULES TO BE MADE BASED ON GEOMETRIC SYMMETRY,
C     ATOMIC COORDINATES AND QUANTUM NUMBER PAIRS. THE IDEA IS TO
C     SKIP CALCULATIONS THAT SHOULD PRODUCE NO OR NEGLIGIBLE EFFECT.
      IF(NCNT.LE.2) THEN
        IF(MQN(1).NE.MQN(2)) GOTO 2000
      ENDIF
C
C**********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS           C
C     OF (MA,MB). FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE          C
C     ORDERED                                                          C
C     11: = (-|MQN(A)|,-|MQN(B)|)                                      C
C     12: = (-|MQN(A)|,+|MQN(B)|)                                      C
C     21: = (+|MQN(A)|,-|MQN(B)|)                                      C
C     22: = (+|MQN(A)|,+|MQN(B)|)                                      C
C**********************************************************************C
C
C     CALCULATE COMPONENT OFFSETS
      IL1 = LARGE(ICNTA,KA,MJA  )
      IL2 = LARGE(ICNTA,KA,MJA+1)
      JL1 = LARGE(ICNTB,KB,MJB  )
      JL2 = LARGE(ICNTB,KB,MJB+1)
C
      IS1 = IL1 + NSHIFT
      IS2 = IL2 + NSHIFT
      JS1 = JL1 + NSHIFT
      JS2 = JL2 + NSHIFT
C
C     CONSTRUCTION OF ONE-ELECTRON OVERLAP MATRIX BY TT' BLOCKS...
C
C     THE PHASE GENERATES E022 AND E012 COEFFS FROM E011 AND E021
      FASE = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
     &                    *DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXM = NFUNA*NFUNB
C
C**********************************************************************C
C     PART 1: THE LL MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)
      NTUVLL = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ELL0 COEFFICIENTS
      IALT = 1
      CALL CPU_TIME(TBEG)
      CALL EMAKELL(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,1,2,0)
      CALL CPU_TIME(TFIN)
      TELL = TELL + TFIN - TBEG
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          EIJ    = EXPT(IBAS,1)+EXPT(JBAS,2)
          EROOT  = DSQRT(PI/EIJ)**3
          SLL(IBAS,JBAS,1) = EROOT*E11(M,1)
          SLL(IBAS,JBAS,3) = EROOT*E21(M,1)
          SLL(IBAS,JBAS,2) =-FASE*DCONJG(SLL(IBAS,JBAS,3))
          SLL(IBAS,JBAS,4) = FASE*DCONJG(SLL(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     NON-REL OVERLAP CALCULATIONS COMPLETE
      IF(HMLTN.EQ.'NORL') GOTO 500
C
C**********************************************************************C
C     PART 2: THE SS MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)+2
      NTUVSS = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ESS0 COEFFICIENTS
      IALT = 1
      CALL CPU_TIME(TBEG)
      CALL EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,1,2,0)
      CALL CPU_TIME(TFIN)
      TESS = TESS + TFIN - TBEG
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          EIJ     = EXPT(IBAS,1)+EXPT(JBAS,2)
          EROOT   = DSQRT(PI/EIJ)**3
          SSS(IBAS,JBAS,1) = EROOT*E11(M,1)
          SSS(IBAS,JBAS,3) = EROOT*E21(M,1)
          SSS(IBAS,JBAS,2) =-FASE*DCONJG(SSS(IBAS,JBAS,3))
          SSS(IBAS,JBAS,4) = FASE*DCONJG(SSS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
500   CONTINUE
C
C**********************************************************************C
C     WE NOW HAVE ALL PIECES OF THE OVERLAP MATRIX FOR THIS BLOCK OF   C
C     BASIS FUNCTIONS -- NOW OVERLAY THE RESULTS INTO VIJ.             C
C**********************************************************************C
C
C     CALCULATE COMPONENT OFFSETS
      IL1 = LARGE(ICNTA,KA,MJA  )
      IL2 = LARGE(ICNTA,KA,MJA+1)
      JL1 = LARGE(ICNTB,KB,MJB  )
      JL2 = LARGE(ICNTB,KB,MJB+1)
C
      IS1 = IL1 + NSHIFT
      IS2 = IL2 + NSHIFT
      JS1 = JL1 + NSHIFT
      JS2 = JL2 + NSHIFT
C
C     LL BLOCKS
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            VIJ(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            VIJ(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            VIJ(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            VIJ(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
            VIJ(JL1+JBAS,IL1+IBAS) = DCONJG(VIJ(IL1+IBAS,JL1+JBAS))
            VIJ(JL2+JBAS,IL1+IBAS) = DCONJG(VIJ(IL1+IBAS,JL2+JBAS))
            VIJ(JL1+JBAS,IL2+IBAS) = DCONJG(VIJ(IL2+IBAS,JL1+JBAS))
            VIJ(JL2+JBAS,IL2+IBAS) = DCONJG(VIJ(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            VIJ(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            VIJ(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            VIJ(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            VIJ(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
            VIJ(JL1+JBAS,IL1+IBAS) = DCONJG(VIJ(IL1+IBAS,JL1+JBAS))
            VIJ(JL2+JBAS,IL1+IBAS) = DCONJG(VIJ(IL1+IBAS,JL2+JBAS))
            VIJ(JL1+JBAS,IL2+IBAS) = DCONJG(VIJ(IL2+IBAS,JL1+JBAS))
            VIJ(JL2+JBAS,IL2+IBAS) = DCONJG(VIJ(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     NON-REL OVERLAP MATRIX COMPLETE
      IF(HMLTN.EQ.'NORL') GOTO 600
C
C     SS BLOCKS
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            VIJ(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            VIJ(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            VIJ(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            VIJ(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
            VIJ(JS1+JBAS,IS1+IBAS) = DCONJG(VIJ(IS1+IBAS,JS1+JBAS))
            VIJ(JS2+JBAS,IS1+IBAS) = DCONJG(VIJ(IS1+IBAS,JS2+JBAS))
            VIJ(JS1+JBAS,IS2+IBAS) = DCONJG(VIJ(IS2+IBAS,JS1+JBAS))
            VIJ(JS2+JBAS,IS2+IBAS) = DCONJG(VIJ(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            VIJ(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            VIJ(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            VIJ(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            VIJ(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
            VIJ(JS1+JBAS,IS1+IBAS) = DCONJG(VIJ(IS1+IBAS,JS1+JBAS))
            VIJ(JS2+JBAS,IS1+IBAS) = DCONJG(VIJ(IS1+IBAS,JS2+JBAS))
            VIJ(JS1+JBAS,IS2+IBAS) = DCONJG(VIJ(IS2+IBAS,JS1+JBAS))
            VIJ(JS2+JBAS,IS2+IBAS) = DCONJG(VIJ(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
600   CONTINUE
C
2000  CONTINUE
C
      RETURN
      END
