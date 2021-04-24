      SUBROUTINE MBPT2TEST
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***********************************************************************C
C                                                                       C
C          MM       MM  BBBBBB   PPPPPP  TTTTTTTT         2222          C
C          MMMM   MMMM  BB   BB  PP   PP    TT           22  22         C
C          MM MM MM MM  BB   BB  PP   PP    TT               22         C
C          MM  MMM  MM  BBBBBB   PPPPPP     TT    #####     22          C
C          MM   M   MM  BB   BB  PP         TT             22           C
C          MM       MM  BB   BB  PP         TT            22            C
C          MM       MM  BBBBBB   PP         TT           2222222        C
C                                                                       C
C      MBPT2 CALCULATES THE SECOND-ORDER PAIR CORRELATION ENERGIES      C
C      OF A CLOSED SHELL SYSTEM USING A RELATIVISTIC ADAPTION OF        C
C      THE DIRECT MP2 ALGORITHM OF HEAD-GORDON, POPLE AND FRISCH        C
C      SPINOR INTEGRALS ARE GENERATED IN ATOMIC STMMETRY BLOCKS         C
C      AND TRANSFORMED INTO MATRIX ELEMENTS BY SEQUENTIAL ONE-INDEX     C
C      TRANSFORMATIONS. MBPT EXPLOITS THE SYMMETRY (AB|CD)=(AB|DC)      C
C      BUT NOT THE SYMMETRY (AB|CD)=(CD|AB) (AND PERMUTATIONS)          C
C                                                                       C
C      THE LOOP STRUCTURE IS BASED ON THE FOCK MATRIX ROUTINE SCFMAT,   C
C      WHICH PASSES THROUGH THE UNIQUE LIST OF SPINOR BLOCK LABELS      C  
C                                                                       C
C***********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=7000000)
C
      PARAMETER(MAXOCC=14,MAXVRT=(MDM/2)-MAXOCC)
C
      CHARACTER*4 HMLTN
      COMPLEX*16 C(MDM,MDM),RR(MB2,16)
C

      DIMENSION TEMP1(MAXOCC*MBS,8),
     &          TEMP2(MAXOCC*MBS,8),RIJ(MB2,4,MAXOCC*MAXVRT),
     &          RMAT(MAXOCC*MAXOCC*MAXVRT*MAXVRT),
     &          V1(MBS,MAXOCC*MAXOCC*MAXVRT,2),
     &          V2(MBS,MAXOCC*MAXOCC*MAXVRT,2),
     &          EPAIR(MAXOCC,MAXOCC)
C
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS
      COMMON/MAKE/IEAB,IECD,IERC,IRIJ(MBS,MBS)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,IEQS,IPAR,ICOR,ILEV
C
      COMMON/COEF/C
      COMMON/EIGN/EIGEN(MDM)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),
      DIMENSION KQN(4),LQN(4),MQN(4),NBAS(4),ITQN(2)
      DIMENSION INDEX(MCT,-(MKP+1)/2:(MKP+1)/2,MKP)
C
      NVIRT=NDIM-NSHIFT-NOCC
C
C***********************************************************************C
C     SET LIMIT ON THE LOOP OVER COMPONENTS                             C
C***********************************************************************C
      IF(HMLTN.EQ.'NORL') THEN
       ITSTOP=1
      ELSE
       ITSTOP=2
      ENDIF
C-----------------------------------------------------------------------C
C     INDEXING ROUTINE
C-----------------------------------------------------------------------C
      ICOUNT=0
      DO 50 ICNT=1,NCNT
      DO 50 KN=1,NKAP(ICNT)
      KQNN=KVALS(KN,ICNT)
      MJMAX=2*IABS(KQNN)-1
      DO 50 MJ=1,MJMAX,2
      ICOUNT=ICOUNT+1
      INDEX(ICNT,KQNN,MJ)=ICOUNT
50    CONTINUE
C-----------------------------------------------------------------------C
C
      DO 2000 ICNTA=1,NCNT
      XYZ(1,1)=COORD(1,ICNTA)
      XYZ(2,1)=COORD(2,ICNTA)
      XYZ(3,1)=COORD(3,ICNTA)
      DO 2000 ICNTB=1,ICNTA
      IF(ICNTA.EQ.ICNTB) THEN
       INUCAB=1
      ELSE
       INUCAB=2
      ENDIF
      XYZ(1,2)=COORD(1,ICNTB)
      XYZ(2,2)=COORD(2,ICNTB)
      XYZ(3,2)=COORD(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
C
      DO 2000 KA=1,NKAP(ICNTA)
      KQN(1)=KVALS(KA,ICNTA)
      IF(KQN(1).GT.0) THEN
       LQN(1)=KQN(1)
      ELSE
       LQN(1)=-KQN(1)-1
      ENDIF
C
      NFUNA=NFUNCT(LQN(1)+1,ICNTA)
      NBAS(1)=NFUNA
      DO 70 IBAS=1,NFUNA
      EXPT(IBAS,1)=EXPSET(IBAS,LQN(1)+1,ICNTA)
70    CONTINUE
C
      IBASB=0
C
C     LOOP OVER KQN(B) VALUES
C
      DO 2000 KB=1,NKAP(ICNTB)
      KQN(2)=KVALS(KB,ICNTB)
      IF(KQN(2).GT.0) THEN
       LQN(2)=KQN(2)
      ELSE
       LQN(2)=-KQN(2)-1
      ENDIF
      NFUNB=NFUNCT(LQN(2)+1,ICNTB)
      NBAS(2)=NFUNB
      DO 71 IBAS=1,NFUNB
      EXPT(IBAS,2)=EXPSET(IBAS,LQN(2)+1,ICNTB)
71    CONTINUE
C
C     LOOP OVER |MA| VALUES
C
      DO 2000 MA=1,IABS(KQN(1))
      MJA=(2*MA)-1
      MQN(1)=MJA
C
C     LOOP OVER |MB| VALUES
C
      DO 2000 MB=1,IABS(KQN(2))
      MJB=(2*MB)-1
      MQN(2)=MJB

C
C     CALCULATE NEW BLOCK OF E(AB|  ) COEFFS AT NEXT OPPORTUNITY
      IEAB  = 1
      IABLL = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB)
      IABSS = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB)

C
C     LOOP OVER CENTRES C AND D
C
      DO 1000 ICNTC=1,NCNT
      XYZ(1,3)=COORD(1,ICNTC)
      XYZ(2,3)=COORD(2,ICNTC)
      XYZ(3,3)=COORD(3,ICNTC)
C
      DO 1000 ICNTD=1,NCNT
      XYZ(1,4)=COORD(1,ICNTD)
      XYZ(2,4)=COORD(2,ICNTD)
      XYZ(3,4)=COORD(3,ICNTD)
C
C     TEST WHETHER THIS COMBINATION IS OF ATOM-IN-MOLECULE TYPE
C
      IF(ICNTC.EQ.ICNTD) THEN
       INUCCD=1
      ELSE
       INUCCD=2
      ENDIF
C
C     LOOP OVER KQN(C) VALUES
C
      DO 1000 KC=1,NKAP(ICNTC)
      KQN(3)=KVALS(KC,ICNTC)
      IF(KQN(3).GT.0) THEN
       LQN(3)=KQN(3)
      ELSE
       LQN(3)=-KQN(3)-1
      ENDIF
      NFUNC=NFUNCT(LQN(3)+1,ICNTC)
      NBAS(3)=NFUNC
      DO 72 IBAS=1,NFUNC
      EXPT(IBAS,3)=EXPSET(IBAS,LQN(3)+1,ICNTC)
72    CONTINUE
C
C     LOOP OVER KQN(D) VALUES
C
      DO 1000 KD=1,NKAP(ICNTD)
      KQN(4)=KVALS(KD,ICNTD)
      IF(KQN(4).GT.0) THEN
       LQN(4)=KQN(4)
      ELSE
       LQN(4)=-KQN(4)-1
      ENDIF
      NFUND=NFUNCT(LQN(4)+1,ICNTD)
      NBAS(4)=NFUND
      DO 73 IBAS=1,NFUND
      EXPT(IBAS,4)=EXPSET(IBAS,LQN(4)+1,ICNTD)
73    CONTINUE

C
C     RESET RC(AB|CD) CLASS CALCULATION INDICATORS
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          IRIJ(IBAS,JBAS) = 1
        ENDDO
      ENDDO
C
C     BUILD NEW RC(AB|CD) CLASS AT NEXT OPPORTUNITY
      IF(IEQS.EQ.0) THEN
C       E-COEFFICIENTS NOT SAVED TO LOCAL FILE -> DON'T SAVE RC(AB|CD)
        IERC = 0
      ELSEIF(LQN(1)+LQN(2)+LQN(3)+LQN(4).EQ.0) THEN
C       TURN OFF RC(AB|CD) LOCAL FILE PROCESS (ONE MQN COMBINATION)
        IERC = 0
      ELSE
C       ALLOW RC(AB|CD) LOCAL FILE PROCESS
        IERC = 1
      ENDIF


C
C     LOOP OVER |MC| AND |MD| VALUES
C
      DO 1000 MC=1,IABS(KQN(3))
      MJC=(2*MC)-1
      MQN(3)=MJC
C
      DO 1000 MD=1,IABS(KQN(4))
      MJD=(2*MD)-1
      MQN(4)=MJD
C
C
C     CALCULATE NEW BLOCK OF E(CD|  ) COEFFS AT NEXT OPPORTUNITY
      IECD  = 1
      ICDLL = IAD0LL(ICNTC,ICNTD,KC,KD,MC,MD)
      ICDSS = IAD0SS(ICNTC,ICNTD,KC,KD,MC,MD)

C
      MAXM=NFUNC*NFUND  
C
      IF(ICNTA.EQ.ICNTB) THEN
       INUCAB=1
      ELSE
       INUCAB=2
      ENDIF
      IF(ICNTC.EQ.ICNTD) THEN
       INUCCD=1
      ELSE
       INUCCD=2
      ENDIF
      IF(INUCAB*INUCCD.EQ.1.AND.ICNTA.EQ.ICNTC) THEN
       IATOM=1
      ELSE
       IATOM=0
      ENDIF

C**********************************************************************C
C     IMPLEMENT DIATOMIC SELECTION RULE                                C
C**********************************************************************C
      ITT=0
C
      DO 463 II=1,7,2
      IM=0
      IF (MQN(1).EQ.II) THEN
      IM=IM+1
      ENDIF
      IF (MQN(2).EQ.II) THEN
      IM=IM+1
      ENDIF
      IF (MQN(3).EQ.II) THEN
      IM=IM+1
      ENDIF
      IF (MQN(4).EQ.II) THEN
      IM=IM+1
      ENDIF
      ITT=ITT+IM*IM
463   CONTINUE
C
      IF (ITT.NE.8.AND.ITT.NE.16) GOTO 1999
C
C**********************************************************************C
C
      IQ1=INDEX(ICNTA,KQN(1),MQN(1))
      IQ2=INDEX(ICNTB,KQN(2),MQN(2))
      IQ3=INDEX(ICNTC,KQN(3),MQN(3))
      IQ4=INDEX(ICNTD,KQN(4),MQN(4))
      IF (IQ1.LT.IQ2) GOTO 1999
      IF (IQ3.LT.IQ4) GOTO 1999
C
C     LOOP OVER COMPONENTS.
C
      DO 4000 IT1=1,ITSTOP
      DO 4000 IT2=1,ITSTOP
      ITQN(1)=3*IT1-2
      ITQN(2)=3*IT2-2
      IT14=(IT1-1)*4
C
C     CALCULATE PHASE FACTORS FOR PERMUTING INTEGRALS.
C
      IF((KQN(1)*KQN(2)).GT.0) THEN 
       RKFAC1=1.0D0
      ELSE
       RKFAC1=-1.0D0
      ENDIF
      IF((KQN(3)*KQN(4)).GT.0) THEN 
       RKFAC2=1.0D0
      ELSE
       RKFAC2=-1.0D0
      ENDIF
C
      FACAB1=DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      FACAB2=DFLOAT((-1)**((MQN(1)+MQN(2))/2))
      FACCD1=DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      FACCD2=DFLOAT((-1)**((MQN(3)+MQN(4))/2))
C
C
      IF(ITQN(1).EQ.1) THEN
       NADDAB=0
      ELSE
       NADDAB=NSHIFT
      ENDIF
      IF(ITQN(2).EQ.1) THEN
       NADDCD=0
      ELSE
       NADDCD=NSHIFT
      ENDIF
C
      IA1=LARGE(ICNTA,KA,MJA)+NADDAB
      IB1=LARGE(ICNTB,KB,MJB)+NADDAB
      IC1=LARGE(ICNTC,KC,MJC)+NADDCD
      ID1=LARGE(ICNTD,KD,MJD)+NADDCD
C
      IA2=LARGE(ICNTA,KA,MJA+1)+NADDAB
      IB2=LARGE(ICNTB,KB,MJB+1)+NADDAB
      IC2=LARGE(ICNTC,KC,MJC+1)+NADDCD
      ID2=LARGE(ICNTD,KD,MJD+1)+NADDCD
C
C
      IEMAKE=1
      ITYPE=1
      MAXM1=NFUNA*NFUNB
      MAXM2=NFUNC*NFUND
      ICOUNT=0
C
      DO 4003 IOV=1,NOCC*NVIRT
      DO 4003 KL=1,4
      DO 4003 M=1,NFUNA*NFUNB
      RIJ(M,KL,IOV)=0.0D0
4003  CONTINUE
C
      DO 3000 IBAS=1,NFUNA
      DO 5000 JBAS=1,NFUNB
      IJBAS=(IBAS-1)*NFUNB+JBAS

C*********************************************************************C
C     CALCULATE A BATCH OF INTEGRALS (IJ|KL) FOR FIXED (IJ)           C
C*********************************************************************C
      CALL ERI(RR,XYZ,KQN,MQN,EXPT,NBAS,ITQN,IBAS,JBAS)
      IEMAKE=0
C
C     PERMUTATIONAL PHASE FACTORS
C
      F1=RKFAC2*FACCD1
      F2=RKFAC2*FACCD2
      G1=RKFAC1*FACAB1
      G2=RKFAC1*FACAB2
C
      DO 4001 KL=1,8
      DO 4001 KIV=1,NFUNC*NOCC
      TEMP1(KIV,KL)=0.0D0
4001  CONTINUE
      DO 4002 KL=1,8
      DO 4002 KIV=1,NFUND*NOCC
      TEMP2(KIV,KL)=0.0D0
4002  CONTINUE
C
C*********************************************************************C
C     PERFORM THE ONE-INDEX TRANSFORMATION OVER ALL POSSIBILITIES     C
C                                                                     C
C                      (IJ|KL) -> (IJ|KO)                             C
C                                                                     C
C*********************************************************************C
C
      NVSHFT=NSHIFT+NOCC
C
      DO 4010 K=1,NFUNC
      DO 4010 IOCC=1,NOCC
      KIO=(K-1)*NOCC+IOCC
      DO 4009 L=1,NFUND
      M=(K-1)*NFUND+L
      TEMP1(KIO,1)=TEMP1(KIO,1)+RR(M,1) *C(ID1+L,IOCC+NSHIFT)
     &                         +RR(M,2) *C(ID2+L,IOCC+NSHIFT)
      TEMP1(KIO,2)=TEMP1(KIO,2)+RR(M,9) *C(ID1+L,IOCC+NSHIFT)
     &                         +RR(M,10)*C(ID2+L,IOCC+NSHIFT)
      TEMP1(KIO,3)=TEMP1(KIO,3)+RR(M,5) *C(ID1+L,IOCC+NSHIFT)
     &                         +RR(M,6) *C(ID2+L,IOCC+NSHIFT)
      TEMP1(KIO,4)=TEMP1(KIO,4)+RR(M,13)*C(ID1+L,IOCC+NSHIFT)
     &                         +RR(M,14)*C(ID2+L,IOCC+NSHIFT)
      TEMP1(KIO,5)=TEMP1(KIO,5)+RR(M,3) *C(ID1+L,IOCC+NSHIFT)
     &                         +RR(M,4) *C(ID2+L,IOCC+NSHIFT)
      TEMP1(KIO,6)=TEMP1(KIO,6)+RR(M,11)*C(ID1+L,IOCC+NSHIFT)
     &                         +RR(M,12)*C(ID2+L,IOCC+NSHIFT)
      TEMP1(KIO,7)=TEMP1(KIO,7)+RR(M,7) *C(ID1+L,IOCC+NSHIFT)
     &                         +RR(M,8) *C(ID2+L,IOCC+NSHIFT)
      TEMP1(KIO,8)=TEMP1(KIO,8)+RR(M,15)*C(ID1+L,IOCC+NSHIFT)
     &                         +RR(M,16)*C(ID2+L,IOCC+NSHIFT)
4009  CONTINUE
4010  CONTINUE
      
C
C*********************************************************************C
C     PERFORM THE ONE-INDEX TRANSFORMATION OVER ALL POSSIBILITIES     C
C                                                                     C
C                      (IJ|LK) -> (IJ|LO)                             C
C                                                                     C
C                ONLY PERFORMED IF IQ3.NEQ.IQ4                        C
C*********************************************************************C
C
C
      IF(IQ3.NE.IQ4) THEN
      DO 4011 L=1,NFUND
      DO 4011 IOCC=1,NOCC
      DO 4008 K=1,NFUNC
      LIO=(L-1)*NOCC+IOCC
      M=(K-1)*NFUND+L
      TEMP2(LIO,1)=TEMP2(LIO,1)+F1*RR(M,4) *C(IC1+K,IOCC+NSHIFT)
     &                         +F2*RR(M,2) *C(IC2+K,IOCC+NSHIFT)
      TEMP2(LIO,2)=TEMP2(LIO,2)+F1*RR(M,12)*C(IC1+K,IOCC+NSHIFT)
     &                         +F2*RR(M,10)*C(IC2+K,IOCC+NSHIFT)
      TEMP2(LIO,3)=TEMP2(LIO,3)+F1*RR(M,8) *C(IC1+K,IOCC+NSHIFT)
     &                         +F2*RR(M,6) *C(IC2+K,IOCC+NSHIFT)
      TEMP2(LIO,4)=TEMP2(LIO,4)+F1*RR(M,16)*C(IC1+K,IOCC+NSHIFT)
     &                         +F2*RR(M,14)*C(IC2+K,IOCC+NSHIFT)
      TEMP2(LIO,5)=TEMP2(LIO,5)+F2*RR(M,3) *C(IC1+K,IOCC+NSHIFT)
     &                         +F1*RR(M,1) *C(IC2+K,IOCC+NSHIFT)
      TEMP2(LIO,6)=TEMP2(LIO,6)+F2*RR(M,11)*C(IC1+K,IOCC+NSHIFT)
     &                         +F1*RR(M,9) *C(IC2+K,IOCC+NSHIFT)
      TEMP2(LIO,7)=TEMP2(LIO,7)+F2*RR(M,7) *C(IC1+K,IOCC+NSHIFT)
     &                         +F1*RR(M,5) *C(IC2+K,IOCC+NSHIFT)
      TEMP2(LIO,8)=TEMP2(LIO,8)+F2*RR(M,15)*C(IC1+K,IOCC+NSHIFT)
     &                         +F1*RR(M,13)*C(IC2+K,IOCC+NSHIFT)
4008  CONTINUE
4011  CONTINUE
      ENDIF
C
C*********************************************************************C
C     PERFORM THE ONE-INDEX TRANSFORMATION OVER ALL POSSIBILITIES     C
C                                                                     C
C                      (IJ|KO) -> (IJ|VO)                             C
C                                                                     C
C*********************************************************************C
C
       DO 4020 IOCC=1,NOCC
       DO 4020 IV=1,NVIRT
       IOV=(IOCC-1)*NVIRT+IV
       DO 4019 K=1,NFUNC
       KIO=(K-1)*NOCC+IOCC
       RIJ(IJBAS,1,IOV)=RIJ(IJBAS,1,IOV)
     &                  +TEMP1(KIO,1)*C(IC1+K,IV+NVSHFT)
     &                  +TEMP1(KIO,5)*C(IC2+K,IV+NVSHFT)
       RIJ(IJBAS,2,IOV)=RIJ(IJBAS,2,IOV)
     &                  +TEMP1(KIO,3)*C(IC1+K,IV+NVSHFT)
     &                  +TEMP1(KIO,7)*C(IC2+K,IV+NVSHFT)
       RIJ(IJBAS,3,IOV)=RIJ(IJBAS,3,IOV)
     &                  +TEMP1(KIO,2)*C(IC1+K,IV+NVSHFT)
     &                  +TEMP1(KIO,6)*C(IC2+K,IV+NVSHFT)
       RIJ(IJBAS,4,IOV)=RIJ(IJBAS,4,IOV)
     &                  +TEMP1(KIO,4)*C(IC1+K,IV+NVSHFT)
     &                  +TEMP1(KIO,8)*C(IC2+K,IV+NVSHFT)
4019  CONTINUE
4020  CONTINUE
C
C
C*********************************************************************C
C     PERFORM THE ONE-INDEX TRANSFORMATION OVER ALL POSSIBILITIES     C
C                                                                     C
C                      (IJ|LO) -> (IJ|VO)                             C
C                                                                     C
C                ONLY PERFORMED IF IQ3.NEQ.IQ4                        C
C*********************************************************************C
C
       IF(IQ3.NE.IQ4) THEN
       DO 4025 IOCC=1,NOCC
       DO 4025 IV=1,NVIRT
       IVO=(IOCC-1)*NVIRT+IV
       DO 4024 L=1,NFUND
       LIO=(L-1)*NOCC+IOCC
       RIJ(IJBAS,1,IVO)=RIJ(IJBAS,1,IVO)
     &                   +TEMP2(LIO,1)*C(ID1+L,IV+NVSHFT)
     &                   +TEMP2(LIO,5)*C(ID2+L,IV+NVSHFT)
       RIJ(IJBAS,2,IVO)=RIJ(IJBAS,2,IVO)
     &                   +TEMP2(LIO,3)*C(ID1+L,IV+NVSHFT)
     &                   +TEMP2(LIO,7)*C(ID2+L,IV+NVSHFT)
       RIJ(IJBAS,3,IVO)=RIJ(IJBAS,3,IVO)
     &                   +TEMP2(LIO,2)*C(ID1+L,IV+NVSHFT)
     &                   +TEMP2(LIO,6)*C(ID2+L,IV+NVSHFT)
       RIJ(IJBAS,4,IVO)=RIJ(IJBAS,4,IVO)
     &                   +TEMP2(LIO,4)*C(ID1+L,IV+NVSHFT)
     &                   +TEMP2(LIO,8)*C(ID2+L,IV+NVSHFT)
4024  CONTINUE
4025  CONTINUE
      ENDIF
C
5000  CONTINUE
3000  CONTINUE
C
C*********************************************************************C
C     THE ARRAY RIJ CONTAINS THE PARTIALLY-TRANSFORMED LIST           C
C                            (IJ,OV)                                  C 
C*********************************************************************C
C
C
C*********************************************************************C
C     PERFORM THE ONE-INDEX TRANSFORMATION OVER ALL POSSIBILITIES     C
C                                                                     C
C                      (IJ|VO) -> (IO|VO)                             C
C                                                                     C
C*********************************************************************C
C
      DO 4004 KL=1,NVIRT*NOCC*NOCC
      DO 4004 M=1,NFUNA
      V1(M,KL,1)=0.0D0
      V1(M,KL,2)=0.0D0
4004  CONTINUE
C
      DO 4031 IOCC1=1,NOCC
      DO 4030 IOCC2=1,IOCC1
      DO 4030 IVIRT2=1,NVIRT
      IVO=(IOCC2-1)*NVIRT+IVIRT2
      IOVO=(IOCC1-1)*NVIRT*NOCC+IVO
      DO 4030 IBAS=1,NFUNA
      DO 4030 JBAS=1,NFUNB
      IJBAS=(IBAS-1)*NFUNB+JBAS
      V1(IBAS,IOVO,1)=V1(IBAS,IOVO,1)
     &                     +RIJ(IJBAS,1,IVO)*C(IB1+JBAS,IOCC1+NSHIFT)
     &                     +RIJ(IJBAS,2,IVO)*C(IB2+JBAS,IOCC1+NSHIFT)  
      V1(IBAS,IOVO,2)=V1(IBAS,IOVO,2)
     &                     +RIJ(IJBAS,3,IVO)*C(IB1+JBAS,IOCC1+NSHIFT)
     &                     +RIJ(IJBAS,4,IVO)*C(IB2+JBAS,IOCC1+NSHIFT) 
4030  CONTINUE
4031  CONTINUE
C
C*********************************************************************C
C     PERFORM THE ONE-INDEX TRANSFORMATION OVER ALL POSSIBILITIES     C
C                                                                     C
C                      (JI|VO) -> (JO|VO)                             C
C                                                                     C
C                ONLY PERFORMED IF IQ1.NEQ.IQ2                        C
C*********************************************************************C
C
      IF(IQ1.NE.IQ2) THEN
      DO 4005 KL=1,NVIRT*NOCC*NOCC
      DO 4005 M=1,NFUNB
      V2(M,KL,1)=0.0D0
      V2(M,KL,2)=0.0D0
4005  CONTINUE
C
      DO 4033 IOCC1=1,NOCC
      DO 4032 IOCC2=1,IOCC1
      DO 4032 IVIRT2=1,NVIRT
      IVO=(IOCC2-1)*NVIRT+IVIRT2
      IOVO=(IOCC1-1)*NVIRT*NOCC+IVO
      DO 4032 IBAS=1,NFUNA
      DO 4032 JBAS=1,NFUNB
      IJBAS=(IBAS-1)*NFUNB+JBAS
      V2(JBAS,IOVO,1)=V2(JBAS,IOVO,1)
     &       +RIJ(IJBAS,4,IVO)*C(IA1+IBAS,IOCC1+NSHIFT)*G1
     &       +RIJ(IJBAS,2,IVO)*C(IA2+IBAS,IOCC1+NSHIFT)*G2
      V2(JBAS,IOVO,2)=V2(JBAS,IOVO,2)
     &       +RIJ(IJBAS,3,IVO)*C(IA1+IBAS,IOCC1+NSHIFT)*G2
     &       +RIJ(IJBAS,1,IVO)*C(IA2+IBAS,IOCC1+NSHIFT)*G1
4032  CONTINUE
4033  CONTINUE
C
      ENDIF
C
C*********************************************************************C
C     PERFORM THE ONE-INDEX TRANSFORMATION OVER ALL POSSIBILITIES     C
C                                                                     C
C                      (IO|VO) -> (VO|VO)                             C
C                                                                     C
C*********************************************************************C
C 
      M=0
      DO 4040 IOCC1=1,NOCC
      DO 4040 IVIRT1=1,NVIRT
      DO 4040 IOCC2=1,IOCC1
      DO 4040 IVIRT2=1,NVIRT
      IVO=(IOCC2-1)*NVIRT+IVIRT2
      IOVO=(IOCC1-1)*NVIRT*NOCC+IVO
C
      M=M+1
      DO 4040 IBAS=1,NFUNA
      RMAT(M)=RMAT(M)+V1(IBAS,IOVO,1)*C(IA1+IBAS,IVIRT1+NVSHFT)
     &               +V1(IBAS,IOVO,2)*C(IA2+IBAS,IVIRT1+NVSHFT)
4040  CONTINUE
C
C*********************************************************************C
C     PERFORM THE ONE-INDEX TRANSFORMATION OVER ALL POSSIBILITIES     C
C                                                                     C
C                      (JO|VO) -> (VO|VO)                             C
C                                                                     C
C                ONLY PERFORMED IF IQ1.NEQ.IQ2                        C
C*********************************************************************C
C
      IF(IQ1.NE.IQ2) THEN
C
      M=0
      DO 4041 IOCC1=1,NOCC
      DO 4041 IVIRT1=1,NVIRT
      DO 4041 IOCC2=1,IOCC1
      DO 4041 IVIRT2=1,NVIRT
      IVO=(IOCC2-1)*NVIRT+IVIRT2
      IOVO=(IOCC1-1)*NVIRT*NOCC+IVO
C
      M=M+1
      DO 4041 JBAS=1,NFUNB
      RMAT(M)=RMAT(M)+V2(JBAS,IOVO,1)*C(IB1+JBAS,IVIRT1+NVSHFT)
     &               +V2(JBAS,IOVO,2)*C(IB2+JBAS,IVIRT1+NVSHFT)
4041  CONTINUE
C
      ENDIF
C*********************************************************************C
C
4000  CONTINUE
1999  CONTINUE
1000  CONTINUE
2000  CONTINUE
C
C*********************************************************************C
C     CALCULATE SECOND-ORDER PAIR CORRELATION ENERGY                  C
C*********************************************************************C
C
      DO 7000 IOCC1=1,NOCC
      DO 7000 IOCC2=1,IOCC1
      EPAIR(IOCC1,IOCC2)=0.0D0
      DO 7000 IVIRT1=1,NVIRT
      DO 7000 IVIRT2=1,NVIRT
C
      IABRS=NVIRT*NVIRT*IOCC1*(IOCC1-1)/2+(IVIRT1-1)*NVIRT*IOCC1+
     &      (IOCC2-1)*NVIRT+IVIRT2
C
      IBARS=NVIRT*NVIRT*IOCC1*(IOCC1-1)/2+(IVIRT2-1)*NVIRT*IOCC1+
     &      (IOCC2-1)*NVIRT+IVIRT1
C
      DABRS=EIGEN(IOCC1+NSHIFT)+EIGEN(IOCC2+NSHIFT)
     &     -EIGEN(IVIRT1+NVSHFT)-EIGEN(IVIRT2+NVSHFT)
C
      EPAIR(IOCC1,IOCC2)=EPAIR(IOCC1,IOCC2)
     &   + (RMAT(IABRS)*RMAT(IABRS)
     &     -RMAT(IABRS)*RMAT(IBARS))/DABRS
C
7000  CONTINUE
C
C********************************************************************C
C     OUTPUT SECOND-ORDER PAIR CORRELATION ENERGIES AND CALCULATE    C
C     THE TOTAL CORRELATION ENERGY                                   C
C********************************************************************C
C
      WRITE(6,7999)
7999  FORMAT(2X,'SECOND ORDER PAIR CORRELATION ENERGIES'//)
      E2TOT=0.0D0
      DO 8000 IOCC1=1,NOCC
      DO 8000 IOCC2=1,IOCC1
      IF(IOCC1.NE.IOCC2) THEN
       FACTOR=1.0D0
      ELSE
       FACTOR=5.0D-1
      ENDIF
      E2TOT=E2TOT+EPAIR(IOCC1,IOCC2)*FACTOR
      WRITE(6,*) IOCC1,IOCC2,EPAIR(IOCC1,IOCC2)
8000  CONTINUE
8001  FORMAT(2X,2(3X,I2),6X, G18.12)
      WRITE(6,8002) E2TOT
8002  FORMAT(//2X,'TOTAL SECOND ORDER ENERGY = ',G18.12)
C
C     
      RETURN
      END

