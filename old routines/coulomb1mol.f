C
C
      SUBROUTINE COULOMB1MOL
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
C  COULOMB GENERATES ALL MANY-CENTER ELECTRON REPULSION INTEGRALS IN   C
C  BATCHES AND ADDS THEM TO THE SCF CLOSED/OPEN-SHELL COULOMB MATRIX.  C
C  CALCULATIONS ARE MADE WITH A RELATIVISTIC MCMURCHIE-DAVIDSON SCHEME.C
C -------------------------------------------------------------------- C
C  THIS IS THE EXPERIMENTAL VERSION. IT CONSTRUCTS ALL THE ONE-CENTER  C
C  COULOMB INTERACTION CONTRIBUTIONS TO THE FOCK MATRIX, BUT BY USING  C
C  THE MOLECULAR INTEGRAL CODES (MCMURCHIE-DAVIDSON ALGORITHM).        C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      CHARACTER*4  HMLTN
      CHARACTER*16 HMS
C
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 RR(MB2,16)
C
      DIMENSION KQN(4),MQN(4),NFUNS(4),LQN(4),ITQN(2),ISEL(16)
      DIMENSION EXPT(MBS,4),XYZ(3,4)
      DIMENSION INDEX(MCT,-MKP:MKP,2*(MKP+1)*MKP),IFLG(11),ISCF(11,6)
      DATA ISCF/1,0,1,0,1,0,1,0,1,1,0,
     &          1,0,0,0,1,0,0,0,1,1,0,
     &          0,1,1,0,0,0,0,0,1,1,0,
     &          1,0,0,1,1,0,0,0,1,0,0,
     &          0,1,0,1,0,0,0,0,1,0,0,
     &          0,1,0,0,0,0,0,0,1,0,0/
C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS
      COMMON/MTRX/OVAP,HNUC,HKIN,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,IEQS,IPAR,ICOR,ILEV
      COMMON/SCRN/F2ES(5,7),T2ES(5,7),TSCR,N2EB(5,7),N2ET(5,7),N2ES(5,7)
      COMMON/SHLL/ACFF,BCFF,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN,NOELEC
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TOGL/IEAB,IECD,IERC,IRIJ(MBS,MBS)
C
      DATA SENS/1.0D-12/
C
      TIMER = 0.0D0
C
C**********************************************************************C
C     LOOP OVER COMPONENT OVERLAP OPTIONS (INDEX 1000)                 C
C**********************************************************************C
C
C     COMPONENT OVERLAP LABELS TO LOOP OVER
      IF(HMLTN.EQ.'BARE') THEN
        RETURN
      ELSEIF(HMLTN.EQ.'NORL') THEN
        ITSTRT = 1
        ITSTOP = 1
        ITSKIP = 1
      ELSE
        ITSTRT = 1
        ITSTOP = 4
        ITSKIP = 3
      ENDIF
C
C     LOOP OVER COMPONENT LABEL FOR A AND B: TT = LL (1) or SS (4)
      DO 1000 IT1=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITQN(1) = IT1
C
C       CALCULATE STARTING ADDRESS
        IF(IT1.EQ.1) THEN
          NADDAB = 0
        ELSE
          NADDAB = NSHIFT
        ENDIF
C
C     LOOP OVER COMPONENT LABEL FOR C AND D: T'T' = LL (1) or SS (4)
      DO 1000 IT2=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITQN(2) = IT2
C
C       CALCULATE STARTING ADDRESS
        IF(IT2.EQ.1) THEN
          NADDCD = 0
        ELSE
          NADDCD = NSHIFT
        ENDIF
C
C     COMPONENT OVERLAP INDEX {(LL|LL)=1,(LL|SS)=2,(SS|LL)=3,(SS|SS)=4}
      ITT = (2*IT1+IT2)/3
C
C**********************************************************************C
C     LOOP OVER ALL ATOMIC CENTERS (USE INDEX 2000)                    C
C**********************************************************************C
C
C     LOOP OVER CENTER A
      DO 2000 ICNT=1,NCNT
C
C       CARTESIAN COORDINATES OF THIS CENTER
        DO N=1,4
          XYZ(1,N) = COORD(1,ICNT)
          XYZ(2,N) = COORD(2,ICNT)
          XYZ(3,N) = COORD(3,ICNT)
        ENDDO
C
C**********************************************************************C
C     LOOP OVER ALL KQN VALUES (USE INDEX 3000)                        C
C**********************************************************************C
C
C     LOOP OVER KQN(A) VALUES
      DO 3000 KA=1,NKAP(ICNT)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNT)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNS(1) = NFUNCT(LQN(1)+1,ICNT)
        DO IBAS=1,NFUNS(1)
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNT)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 3000 KB=1,NKAP(ICNT)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2) = KVALS(KB,ICNT)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNS(2) = NFUNCT(LQN(2)+1,ICNT)
        DO JBAS=1, NFUNS(2)
          EXPT(JBAS,2) = EXPSET(JBAS,LQN(2)+1,ICNT)
        ENDDO
C
C     LOOP OVER KQN(C) VALUES
      DO 3000 KC=1,NKAP(ICNT)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR C
        KQN(3) = KVALS(KC,ICNT)
        IF(KQN(3).GT.0) THEN
          LQN(3) = KQN(3)
        ELSE
          LQN(3) =-KQN(3)-1
        ENDIF
C         
        NFUNS(3) = NFUNCT(LQN(3)+1,ICNT)
        DO KBAS=1,NFUNS(3)
          EXPT(KBAS,3) = EXPSET(KBAS,LQN(3)+1,ICNT)
        ENDDO      
C
C     LOOP OVER KQN(D) VALUES
      DO 3000 KD=1,NKAP(ICNT)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR D
        KQN(4) = KVALS(KD,ICNT)
        IF(KQN(4).GT.0) THEN
          LQN(4) = KQN(4)
        ELSE
          LQN(4) =-KQN(4)-1
        ENDIF
C
        NFUNS(4) = NFUNCT(LQN(4)+1,ICNT)
        DO LBAS=1,NFUNS(4)
          EXPT(LBAS,4) = EXPSET(LBAS,LQN(4)+1,ICNT)
        ENDDO
C
C     RESET RC(AB|CD) CLASS CALCULATION INDICATORS
      DO IBAS=1,NFUNS(1)
        DO JBAS=1,NFUNS(2)
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
C**********************************************************************C
C     LOOP OVER ALL |MQN| VALUES (INDEX 4000)                          C
C**********************************************************************C
C
C     LOOP OVER |MQN(A)| VALUES
      DO 4000 MA=1,IABS(KQN(1))
        MQN(1) = 2*MA-1
C
C     LOOP OVER |MQN(B)| VALUES
      DO 4000 MB=1,IABS(KQN(2))
        MQN(2) = 2*MB-1
C
C     CALCULATE NEW BLOCK OF E(AB|  ) COEFFS AT NEXT OPPORTUNITY
      IEAB  = 1
      IABLL = IAD0LL(ICNT,ICNT,KA,KB,MA,MB)
      IABSS = IAD0SS(ICNT,ICNT,KA,KB,MA,MB)
C
C     LOOP OVER |MQN(C)| VALUES
      DO 4000 MC=1,IABS(KQN(3))
        MQN(3) = 2*MC-1
C
C     LOOP OVER |MQN(D)| VALUES
      DO 4000 MD=1,IABS(KQN(4))
        MQN(4) = 2*MD-1
C
C     CALCULATE NEW BLOCK OF E(CD|  ) COEFFS AT NEXT OPPORTUNITY
      IECD  = 1
      ICDLL = IAD0LL(ICNT,ICNT,KC,KD,MC,MD)
      ICDSS = IAD0SS(ICNT,ICNT,KC,KD,MC,MD)
C
C
C**********************************************************************C
C     INDEX ASSIGNMENT AND IDENTIFICATION OF ERI SYMMETRY CLASSES      C
C**********************************************************************C
C
C     CALCULATE BLOCK INDICES FOR {ABCD} COMBINATIONS
      IQ1 = INDEX(ICNT,KQN(1),MQN(1))
      IQ2 = INDEX(ICNT,KQN(2),MQN(2))
      IQ3 = INDEX(ICNT,KQN(3),MQN(3))
      IQ4 = INDEX(ICNT,KQN(4),MQN(4))
C
C     COMBINED BLOCK INDEX IN A TWO-FUNCTION LIST
      IQL = (IQ1*(IQ1-1))/2 + IQ2
      IQR = (IQ3*(IQ3-1))/2 + IQ4
C
C     GMAT AND QMAT ARE HERMITIAN, SO SKIP UPPER TRI-DIAGONAL BITS
      IF(IQ1.LT.IQ2) GOTO 4001
      IF(IQ3.LT.IQ4) GOTO 4001
      IF(IQL.LT.IQR) GOTO 4001
C
      IF(IQ1.GT.IQ2) THEN
C       IQ1 > IQ2
        IF(IQ3.GT.IQ4) THEN
C         IQ3 > IQ4
          IF(IQL.GT.IQR) THEN
C           IQL > IQR
            ITSCF = 1
          ELSEIF(IQL.EQ.IQR) THEN
C           IQL = IQR
            ITSCF = 2
          ENDIF
        ELSEIF(IQ3.EQ.IQ4) THEN
C         IQ3 = IQ4
          IF(IQL.GT.IQR) THEN
C           IQL > IQR
            ITSCF = 3
          ELSE
C           IQL = IQR
            GOTO 4001
          ENDIF
        ENDIF
      ELSEIF(IQ1.EQ.IQ2) THEN
C       IQ1 = IQ2
        IF(IQ3.GT.IQ4) THEN
C         IQ3 > IQ4
          IF(IQL.GT.IQR) THEN
C           IQL > IQR          
            ITSCF = 4
          ELSE
C           IQL = IQR
            GOTO 4001
          ENDIF
        ELSEIF(IQ3.EQ.IQ4) THEN
C         IQ3 = IQ4
          IF(IQL.GT.IQR) THEN
C           IQL > IQR
            ITSCF = 5
          ELSEIF(IQL.EQ.IQR) THEN
C           IQL = IQR
            ITSCF = 6
          ENDIF
        ENDIF
      ENDIF
C
C     READ IN FLAG VALUES FROM ISCF DATA BLOCK
      DO N=1,11
        IFLG(N) = ISCF(N,ITSCF)
      ENDDO
C
C     INCLUDE SPECIAL CASES FOR MATCHING BLOCKS
      IF(ITSCF.EQ.1) THEN
        IF(IQ1.EQ.IQ3) THEN
C         IQ1 = IQ3
          IFLG( 6) = 1
        ENDIF
        IF(IQ2.EQ.IQ4) THEN
C         IQ2 = IQ3
          IFLG(11) = 1
        ENDIF
        IF(IQ2.EQ.IQ3) THEN
C         IQ2 = IQ3
          IFLG( 8) = 1
        ENDIF
      ELSEIF(ITSCF.EQ.3) THEN
        IF(IQ2.EQ.IQ3) THEN
C         IQ2 = IQ3
          IFLG( 8) = 1
        ENDIF
      ELSEIF(ITSCF.EQ.4) THEN
        IF(IQ2.EQ.IQ3) THEN
C         IQ2 = IQ3
          IFLG( 8) = 1
        ENDIF
      ENDIF
C
C     UPDATE COUNTER FOR NUMBER OF CLASSES
      N2EB(1,ITT) = N2EB(1,ITT) + 1
C
C
C**********************************************************************C
C     PHASE FACTORS APPROPRIATE TO INTEGRAL SYMMETRY CLASSES           C
C**********************************************************************C
C
C     EQ-COEFFICIENT PHASE FACTORS THAT YIELD NEW R-INTEGRALS
      PAB1 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PAB2 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      PCD1 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PCD2 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
C     FOCK ADDRESS FOR EACH BASIS FUNCTION (WITH SPIN PROJECTION)
      NA1 = LARGE(ICNT,KA,2*MA-1) + NADDAB
      NB1 = LARGE(ICNT,KB,2*MB-1) + NADDAB
      NC1 = LARGE(ICNT,KC,2*MC-1) + NADDCD
      ND1 = LARGE(ICNT,KD,2*MD-1) + NADDCD
C
      NA2 = LARGE(ICNT,KA,2*MA  ) + NADDAB
      NB2 = LARGE(ICNT,KB,2*MB  ) + NADDAB
      NC2 = LARGE(ICNT,KC,2*MC  ) + NADDCD
      ND2 = LARGE(ICNT,KD,2*MD  ) + NADDCD
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS IN BLOCKS A AND B (INDEX 5000)         C
C**********************************************************************C
C
C     LOOP OVER ELEMENTS OF FOCK MATRIX BLOCK
      DO 5000 IBAS=1,NFUNS(1)
      DO 5000 JBAS=1,NFUNS(2)
C
      CALL CPU_TIME(TBCH1)
C     UPDATE COUNTER FOR NUMBER OF BATCHES
      N2ET(1,ITT) = N2ET(1,ITT) + 1
C
C**********************************************************************C
C     GENERATE A BATCH OF ERIS OVER ALL |MQN| SIGNS                    C
C**********************************************************************C
C
C     GENERATE BATCH OF ELECTRON REPULSION INTEGRALS
      CALL CPU_TIME(T1)
      CALL ERI(RR,XYZ,KQN,MQN,EXPT,NFUNS,ITQN,IBAS,JBAS)
      CALL CPU_TIME(T2)
      TIMER = TIMER + T2 - T1
C
C     FIRST IFLG BATCH (DIRECT)
      IF(IFLG(1).EQ.1) THEN
        M = 0
        DO KBAS=1,NFUNS(3)
          DO LBAS=1,NFUNS(4)
            M = M+1
C
            GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 1)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 2)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 3)*DENC(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 4)*DENC(NC2+KBAS,ND2+LBAS)
     &                     +      PCD1*RR(M, 4)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 2)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 3)*DENC(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENC(ND2+LBAS,NC2+KBAS)
C
            GDIR(NA1+IBAS,NB2+JBAS) = GDIR(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 5)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 6)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 7)*DENC(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 8)*DENC(NC2+KBAS,ND2+LBAS)
     &                     +      PCD1*RR(M, 8)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 7)*DENC(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 5)*DENC(ND2+LBAS,NC2+KBAS)

C
            GDIR(NA2+IBAS,NB1+JBAS) = GDIR(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M, 9)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,10)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,11)*DENC(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M,12)*DENC(NC2+KBAS,ND2+LBAS)
     &                     +      PCD1*RR(M,12)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENC(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENC(ND2+LBAS,NC2+KBAS)
C
            GDIR(NA2+IBAS,NB2+JBAS) = GDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(M,13)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,14)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,15)*DENC(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M,16)*DENC(NC2+KBAS,ND2+LBAS)
     &                     +      PCD1*RR(M,16)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,14)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENC(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENC(ND2+LBAS,NC2+KBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     SECOND IFLG BATCH (DIRECT)
      IF(IFLG(2).EQ.1) THEN
        M = 0
        DO KBAS=1,NFUNS(3)
          DO LBAS=1,NFUNS(4)
            M = M+1
C
            GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 1)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 2)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 3)*DENC(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 4)*DENC(NC2+KBAS,ND2+LBAS)
C
            GDIR(NA1+IBAS,NB2+JBAS) = GDIR(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 5)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 6)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 7)*DENC(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 8)*DENC(NC2+KBAS,ND2+LBAS)
C
            GDIR(NA2+IBAS,NB1+JBAS) = GDIR(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M, 9)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,10)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,11)*DENC(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M,12)*DENC(NC2+KBAS,ND2+LBAS)
C
            GDIR(NA2+IBAS,NB2+JBAS) = GDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(M,13)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,14)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,15)*DENC(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M,16)*DENC(NC2+KBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     THIRD IFLG BATCH (DIRECT)
      IF(IFLG(3).EQ.1) THEN
        M = 0
        DO KBAS=1,NFUNS(3)
          DO LBAS=1,NFUNS(4)
            M = M+1
C
            GDIR(NC1+KBAS,ND1+LBAS) = GDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 1)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 5)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 9)*DENC(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,13)*DENC(NA2+IBAS,NB2+JBAS)
     &                     + PAB1*     RR(M,13)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 5)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 9)*DENC(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENC(NB2+JBAS,NA2+IBAS)
C
            GDIR(NC1+KBAS,ND2+LBAS) = GDIR(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 2)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 6)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,10)*DENC(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,14)*DENC(NA2+IBAS,NB2+JBAS)
     &                     + PAB1*     RR(M,14)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,10)*DENC(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 2)*DENC(NB2+JBAS,NA2+IBAS)
C
            GDIR(NC2+KBAS,ND1+LBAS) = GDIR(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 3)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 7)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,11)*DENC(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,15)*DENC(NA2+IBAS,NB2+JBAS)
     &                     + PAB1*     RR(M,15)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENC(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENC(NB2+JBAS,NA2+IBAS)
C
            GDIR(NC2+KBAS,ND2+LBAS) = GDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(M, 4)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 8)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,12)*DENC(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,16)*DENC(NA2+IBAS,NB2+JBAS)
     &                     + PAB1*     RR(M,16)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 8)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENC(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENC(NB2+JBAS,NA2+IBAS)
C 
          ENDDO
        ENDDO
      ENDIF
C
C     FOURTH IFLG BATCH (DIRECT)
      IF(IFLG(4).EQ.1) THEN
        M = 0
        DO KBAS=1,NFUNS(3)
          DO LBAS=1,NFUNS(4)
            M = M+1
C
            GDIR(NC1+KBAS,ND1+LBAS) = GDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 1)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 5)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 9)*DENC(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,13)*DENC(NA2+IBAS,NB2+JBAS)
C
            GDIR(NC1+KBAS,ND2+LBAS) = GDIR(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 2)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 6)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,10)*DENC(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,14)*DENC(NA2+IBAS,NB2+JBAS)
C
            GDIR(NC2+KBAS,ND1+LBAS) = GDIR(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 3)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 7)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,11)*DENC(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,15)*DENC(NA2+IBAS,NB2+JBAS)
C
            GDIR(NC2+KBAS,ND2+LBAS) = GDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(M, 4)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 8)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,12)*DENC(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,16)*DENC(NA2+IBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     FIFTH IFLG BATCH (EXCHANGE)
      IF(IFLG(5).EQ.1) THEN
        M = 0
        DO KBAS=1,NFUNS(3)
          DO LBAS=1,NFUNS(4)
            M = M+1
C
            GXCH(NA1+IBAS,NC1+KBAS) = GXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 4)*DENC(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 8)*DENC(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M, 3)*DENC(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 7)*DENC(ND2+LBAS,NB2+JBAS)
C
            GXCH(NA1+IBAS,NC2+KBAS) = GXCH(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 2)*DENC(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 6)*DENC(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 1)*DENC(ND2+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 5)*DENC(ND2+LBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,NC1+KBAS) = GXCH(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,12)*DENC(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,16)*DENC(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M,11)*DENC(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M,15)*DENC(ND2+LBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,NC2+KBAS) = GXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,10)*DENC(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M,14)*DENC(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 9)*DENC(ND2+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,13)*DENC(ND2+LBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     SIXTH IFLG BATCH (EXCHANGE)
      IF(IFLG(6).EQ.1) THEN
        M = 0
        DO KBAS=1,NFUNS(3)
          DO LBAS=1,NFUNS(4)
            M = M+1
C
            GXCH(NC1+KBAS,NA1+IBAS) = GXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,13)*DENC(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,14)*DENC(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M, 9)*DENC(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M,10)*DENC(NB2+JBAS,ND2+LBAS)
C
            GXCH(NC1+KBAS,NA2+IBAS) = GXCH(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 5)*DENC(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 6)*DENC(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M, 1)*DENC(NB2+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M, 2)*DENC(NB2+JBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NA1+IBAS) = GXCH(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,15)*DENC(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,16)*DENC(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M,11)*DENC(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M,12)*DENC(NB2+JBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NA2+IBAS) = GXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 7)*DENC(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 8)*DENC(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M, 3)*DENC(NB2+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M, 4)*DENC(NB2+JBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     SEVENTH IFLG BATCH (EXCHANGE)
      IF(IFLG(7).EQ.1) THEN
        M = 0
        DO KBAS=1,NFUNS(3)
          DO LBAS=1,NFUNS(4)
            M = M+1
C
            GXCH(NB1+JBAS,NC1+KBAS) = GXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB1*PCD1*RR(M,16)*DENC(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD1*RR(M, 8)*DENC(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD2*RR(M,15)*DENC(ND2+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(M, 7)*DENC(ND2+LBAS,NA2+IBAS)
C
            GXCH(NB1+JBAS,NC2+KBAS) = GXCH(NB1+JBAS,NC2+KBAS)
     &                     + PAB1*PCD2*RR(M,14)*DENC(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(M, 6)*DENC(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD1*RR(M,13)*DENC(ND2+LBAS,NA1+IBAS)
     &                     + PAB2*PCD1*RR(M, 5)*DENC(ND2+LBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,NC1+KBAS) = GXCH(NB2+JBAS,NC1+KBAS)
     &                     + PAB2*PCD1*RR(M,12)*DENC(ND1+LBAS,NA1+IBAS)
     &                     + PAB1*PCD1*RR(M, 4)*DENC(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD2*RR(M,11)*DENC(ND2+LBAS,NA1+IBAS)
     &                     + PAB1*PCD2*RR(M, 3)*DENC(ND2+LBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,NC2+KBAS) = GXCH(NB2+JBAS,NC2+KBAS)
     &                     + PAB1*PCD1*RR(M, 1)*DENC(ND2+LBAS,NA2+IBAS)
     &                     + PAB1*PCD2*RR(M, 2)*DENC(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD1*RR(M, 9)*DENC(ND2+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(M,10)*DENC(ND1+LBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     EIGHTH IFLG BATCH (EXCHANGE)
      IF(IFLG(8).EQ.1) THEN
        M = 0
        DO KBAS=1,NFUNS(3)
          DO LBAS=1,NFUNS(4)
            M = M+1
C
            GXCH(NC1+KBAS,NB1+JBAS) = GXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 1)*DENC(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 2)*DENC(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M, 9)*DENC(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M,10)*DENC(NA2+IBAS,ND2+LBAS)
C
            GXCH(NC1+KBAS,NB2+JBAS) = GXCH(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 5)*DENC(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 6)*DENC(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,13)*DENC(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M,14)*DENC(NA2+IBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NB1+JBAS) = GXCH(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M, 3)*DENC(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 4)*DENC(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,11)*DENC(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M,12)*DENC(NA2+IBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NB2+JBAS) = GXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(M, 7)*DENC(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 8)*DENC(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,15)*DENC(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M,16)*DENC(NA2+IBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     NINTH IFLG BATCH (EXCHANGE)
      IF(IFLG(9).EQ.1) THEN
        DO LBAS=1,NFUNS(4)
          DO KBAS=1,NFUNS(3)
            M = (KBAS-1)*NFUNS(4)+LBAS
C  
            GXCH(NA1+IBAS,ND1+LBAS) = GXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 1)*DENC(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 5)*DENC(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 3)*DENC(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M, 7)*DENC(NC2+KBAS,NB2+JBAS)
C
            GXCH(NA1+IBAS,ND2+LBAS) = GXCH(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M, 2)*DENC(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 6)*DENC(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 4)*DENC(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M, 8)*DENC(NC2+KBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,ND1+LBAS) = GXCH(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M, 9)*DENC(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M,13)*DENC(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M,11)*DENC(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M,15)*DENC(NC2+KBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,ND2+LBAS) = GXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(M,10)*DENC(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M,14)*DENC(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M,12)*DENC(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M,16)*DENC(NC2+KBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     TENTH IFLG BATCH (EXCHANGE)
      IF(IFLG(10).EQ.1) THEN
        DO LBAS=1,NFUNS(4)
          DO KBAS=1,NFUNS(3)
            M = (KBAS-1)*NFUNS(4)+LBAS
C
            GXCH(NB1+JBAS,ND1+LBAS) = GXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,13)*DENC(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 5)*DENC(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(M,15)*DENC(NC2+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENC(NC2+KBAS,NA2+IBAS)
C
            GXCH(NB1+JBAS,ND2+LBAS) = GXCH(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,14)*DENC(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENC(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(M,16)*DENC(NC2+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 8)*DENC(NC2+KBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,ND1+LBAS) = GXCH(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 9)*DENC(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENC(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENC(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENC(NC2+KBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,ND2+LBAS) = GXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M,10)*DENC(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 2)*DENC(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENC(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENC(NC2+KBAS,NA2+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     ELEVENTH IFLG BATCH (EXCHANGE)
      IF(IFLG(11).EQ.1) THEN
        DO LBAS=1,NFUNS(4)
          DO KBAS=1,NFUNS(3)
            M = (KBAS-1)*NFUNS(4)+LBAS
C
            GXCH(ND1+LBAS,NB1+JBAS) = GXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENC(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 2)*DENC(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(M,12)*DENC(NA2+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENC(NA2+IBAS,NC2+KBAS)
C
            GXCH(ND1+LBAS,NB2+JBAS) = GXCH(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 8)*DENC(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENC(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(M,16)*DENC(NA2+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,14)*DENC(NA2+IBAS,NC2+KBAS)
C
            GXCH(ND2+LBAS,NB1+JBAS) = GXCH(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 3)*DENC(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENC(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENC(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENC(NA2+IBAS,NC2+KBAS)
C
            GXCH(ND2+LBAS,NB2+JBAS) = GXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M, 7)*DENC(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 5)*DENC(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENC(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENC(NA2+IBAS,NC2+KBAS)
C
          ENDDO
        ENDDO
      ENDIF      

cC
cC     LOOP OVER THE SIGNS OF |MQN|
c      DO 6000 ISGN1=1,2
c        MMJA = MQN(1)*((-1)**ISGN1)
c        IMJA = MQN(1)+ISGN1-1
cC
c      DO 6000 ISGN2=1,2
c        MMJB = MQN(2)*((-1)**ISGN2)
c        IMJB = MQN(2)+ISGN2-1
cC
c      DO 6000 ISGN3=1,2
c        MMJC = MQN(3)*((-1)**ISGN3)
c        IMJC = MQN(3)+ISGN3-1
cC
c      DO 6000 ISGN4=1,2
c        MMJD = MQN(4)*((-1)**ISGN4)
c        IMJD = MQN(4)+ISGN4-1
cC
cC     APPLY ANGULAR MQN SELECTION RULE
c      IF(MMJA-MMJB.NE.MMJD-MMJC) GO TO 6001
cC
cC     STARTING FOCK ADDRESS FOR EACH BASIS LIST
c      NAT = LARGE(ICNT,KA,IMJA) + NADDAB
c      NBT = LARGE(ICNT,KB,IMJB) + NADDAB
c      NCT = LARGE(ICNT,KC,IMJC) + NADDCD
c      NDT = LARGE(ICNT,KD,IMJD) + NADDCD
cC
cC     RR ADDRESS
c      NR = 8*(ISGN1-1) + 4*(ISGN2-1) + 2*(ISGN3-1) + ISGN4
cC
cC     DIRECT CONTRIBUTIONS
c      M = 0
c      DO KBAS=1,NFUNS(3)
c        DO LBAS=1,NFUNS(4)
c          M = M+1
c          GDIR(NAT+IBAS,NBT+JBAS) = GDIR(NAT+IBAS,NBT+JBAS)
c     &                     +           RR(M,NR)*DENC(NCT+KBAS,NDT+LBAS)
c        ENDDO
cC
cC     EXCHANGE CONTRIBUTIONS
c      M = 0
c      DO KBAS=1,NFUNS(3)
c        DO LBAS=1,NFUNS(4)
c          M = M+1 
c          GXCH(NAT+IBAS,NDT+LBAS) = GXCH(NAT+IBAS,NDT+LBAS)
c     &                     +           RR(M,NR)*DENC(NCT+KBAS,NBT+JBAS)
c        ENDDO
c      ENDDO
c
C        M = 0
C        DO KBAS=1,NFUNS(3)
C          DO LBAS=1,NFUNS(4)
C            M = M+1
CC
C            GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
C     &                     +           RR(M, 1)*DENC(NC1+KBAS,ND1+LBAS)
C     &                     +           RR(M, 2)*DENC(NC1+KBAS,ND2+LBAS)
C     &                     +           RR(M, 3)*DENC(NC2+KBAS,ND1+LBAS)
C     &                     +           RR(M, 4)*DENC(NC2+KBAS,ND2+LBAS)
CC
C            GDIR(NA1+IBAS,NB2+JBAS) = GDIR(NA1+IBAS,NB2+JBAS)
C     &                     +           RR(M, 5)*DENC(NC1+KBAS,ND1+LBAS)
C     &                     +           RR(M, 6)*DENC(NC1+KBAS,ND2+LBAS)
C     &                     +           RR(M, 7)*DENC(NC2+KBAS,ND1+LBAS)
C     &                     +           RR(M, 8)*DENC(NC2+KBAS,ND2+LBAS)
CC
C            GDIR(NA2+IBAS,NB1+JBAS) = GDIR(NA2+IBAS,NB1+JBAS)
C     &                     +           RR(M, 9)*DENC(NC1+KBAS,ND1+LBAS)
C     &                     +           RR(M,10)*DENC(NC1+KBAS,ND2+LBAS)
C     &                     +           RR(M,11)*DENC(NC2+KBAS,ND1+LBAS)
C     &                     +           RR(M,12)*DENC(NC2+KBAS,ND2+LBAS)
CC
C            GDIR(NA2+IBAS,NB2+JBAS) = GDIR(NA2+IBAS,NB2+JBAS)
C     &                     +           RR(M,13)*DENC(NC1+KBAS,ND1+LBAS)
C     &                     +           RR(M,14)*DENC(NC1+KBAS,ND2+LBAS)
C     &                     +           RR(M,15)*DENC(NC2+KBAS,ND1+LBAS)
C     &                     +           RR(M,16)*DENC(NC2+KBAS,ND2+LBAS)
CC
C          ENDDO
C        ENDDO
C
C        DO LBAS=1,NFUNS(4)
C          DO KBAS=1,NFUNS(3)
C            M = (KBAS-1)*NFUNS(4)+LBAS
CC  
C            GXCH(NA1+IBAS,ND1+LBAS) = GXCH(NA1+IBAS,ND1+LBAS)
C     &                     +           RR(M, 1)*DENC(NC1+KBAS,NB1+JBAS)
C     &                     +           RR(M, 5)*DENC(NC1+KBAS,NB2+JBAS)
C     &                     +           RR(M, 3)*DENC(NC2+KBAS,NB1+JBAS)
C     &                     +           RR(M, 7)*DENC(NC2+KBAS,NB2+JBAS)
CC
C            GXCH(NA1+IBAS,ND2+LBAS) = GXCH(NA1+IBAS,ND2+LBAS)
C     &                     +           RR(M, 2)*DENC(NC1+KBAS,NB1+JBAS)
C     &                     +           RR(M, 6)*DENC(NC1+KBAS,NB2+JBAS)
C     &                     +           RR(M, 4)*DENC(NC2+KBAS,NB1+JBAS)
C     &                     +           RR(M, 8)*DENC(NC2+KBAS,NB2+JBAS)
CC
C            GXCH(NA2+IBAS,ND1+LBAS) = GXCH(NA2+IBAS,ND1+LBAS)
C     &                     +           RR(M, 9)*DENC(NC1+KBAS,NB1+JBAS)
C     &                     +           RR(M,13)*DENC(NC1+KBAS,NB2+JBAS)
C     &                     +           RR(M,11)*DENC(NC2+KBAS,NB1+JBAS)
C     &                     +           RR(M,15)*DENC(NC2+KBAS,NB2+JBAS)
CC
C            GXCH(NA2+IBAS,ND2+LBAS) = GXCH(NA2+IBAS,ND2+LBAS)
C     &                     +           RR(M,10)*DENC(NC1+KBAS,NB1+JBAS)
C     &                     +           RR(M,14)*DENC(NC1+KBAS,NB2+JBAS)
C     &                     +           RR(M,12)*DENC(NC2+KBAS,NB1+JBAS)
C     &                     +           RR(M,16)*DENC(NC2+KBAS,NB2+JBAS)
CC
C          ENDDO
C        ENDDO

C
C**********************************************************************C
C     OPEN-SHELL CALCULATIONS ONLY                                     C
C**********************************************************************C
C
c      IF(NOPN.EQ.0) GOTO 5100
cC
cC     DIRECT CONTRIBUTIONS
c      M = 0
c      DO KBAS=1,NFUNS(3)
c        DO LBAS=1,NFUNS(4)
c          M = M+1
c          QDIR(NAT+IBAS,NBT+JBAS) = QDIR(NAT+IBAS,NBT+JBAS)
c     &                     +           RR(M,NR)*DENO(NCT+KBAS,NDT+LBAS)
c        ENDDO
c      ENDDO
cC
cC     EXCHANGE CONTRIBUTIONS
c      M = 0
c      DO KBAS=1,NFUNS(3)
c        DO LBAS=1,NFUNS(4)
c          M = M+1
c          QXCH(NAT+IBAS,NDT+LBAS) = QXCH(NAT+IBAS,NDT+LBAS)
c     &                     +           RR(M,NR)*DENO(NCT+KBAS,NBT+JBAS)
c        ENDDO
c      ENDDO
C
cC     SKIPPING POINT FOR CLOSED SYSTEMS
c5100  CONTINUE
cC     SKIPPING POINT FOR MQN SELECTION RULE
c6001  CONTINUE
cC     END LOOP OVER |MQN| SIGNS
c6000  CONTINUE
C     RECORD CPU TIME AT END OF BATCH AND ADD TO TIME COUNTER
      CALL CPU_TIME(TBCH2)
      T2ES(1,ITT) = T2ES(1,ITT) + TBCH2 - TBCH1
C     END LOOP OVER IBAS AND JBAS
5000  CONTINUE
C     SKIPPING POINT FOR MQN INTEGRAL SYMMETRY
4001  CONTINUE
C     END LOOP OVER ALL |MQN| VALUES
4000  CONTINUE
C     END LOOP OVER KQN SYMMETRY TYPES
3000  CONTINUE
C     END LOOP OVER ATOMIC CENTER
2000  CONTINUE
C     END LOOP OVER COMPONENT OVERLAPS
1000  CONTINUE
C
      WRITE(*,*) 'TIME IN ONE-CENTER ERI',HMS(TIMER)
C
      RETURN
      END
C
C     
      SUBROUTINE ERI1(RR,KQN,MQN,EXPT,NFUNS,ITQN,IBAS,JBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                    EEEEEEEE RRRRRRR  IIII   11                       C
C                    EE       RR    RR  II   111                       C
C                    EE       RR    RR  II    11                       C
C                    EEEEEE   RR    RR  II    11                       C
C                    EE       RRRRRRR   II    11                       C
C                    EE       RR    RR  II    11                       C
C                    EEEEEEEE RR    RR IIII 111111                     C
C                                                                      C
C -------------------------------------------------------------------- C
C  ERI1 GENERATES A BATCH OF ATOMIC ELECTRON REPULSION INTEGRALS BY    C
C  MEANS OF INCOMPLETE BETA FUNCTIONS AND WIGNER-3J SYMBOLS.           C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    XYZ(3,4)    - COORDINATES OF THE NUCLEAR CENTERS IN THIS BLOCK.   C
C    KQN(4)      - KQN RELATIVISTIC LABELS OF THE CENTERS.             C
C    MQN(4)      - |MQN| QUANTUM NUMBERS OF THE CENTERS.               C
C    EXPT(MBS,4) - LIST OF EXPONENTS IN THE BLOCK.                     C
C    NFUNS(4)    - NUMBER OF FUNCTIONS ON CENTER J.                    C
C    ITQN(2)     - SPINOR LABEL PAIR FOR AB (I=1) AND CD (I=2).        C
C                  ITQN(I) = {LL,LS,SL,SS}.                            C
C    IBAS,JBAS   - COMPONENT LABEL INDEX FOR AB BASIS FUNCTIONS.       C
C    IEAB,IECD - 0 DON'T RECALCULATE E(AB/CD)-COEFFICIENTS             C
C                1 DO    RECALCULATE E(AB/CD)-COEFFICIENTS             C
C  OUTPUT:                                                             C
C    RR(MB2,16) - ERI'S FOR BLOCK AB, ALL 16 MQN PROJECTION COMBOS.    C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MNU=MKP+1,MAB=2*MNU+5)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 RR(MB2,16)
C
      DIMENSION KQN(4),LQN(4),MQN(4),ITQN(2),NFUNS(4),EXPT(MBS,4)
      DIMENSION PQ(MB2,3),APH(MB2),PRE(MB2),XYZ(3,4)
C
      COMMON/ANGL/BK(MNU,4),DK(MNU,4),HK(MNU,4),FK(MNU,4),GM(MNU,4),
     &            NUS(MNU),NUNUM,NUMIN,NUMAX
      COMMON/BSIJ/EIJ(MAB),RNIJ(3),EI,EJ
      COMMON/BSKL/EK(MB2),EKL0(MB2),EIK(MB2,MAB),EJL(MB2,MAB),
     &            EL(MB2),EKL1(MB2),EKL(MB2,MAB),RNKL(MB2,4),
     &            B1(MB2,MAB,MAB),B2(MB2,MAB,MAB),B3(MB2,MAB,MAB),
     &            B4(MB2,MAB,MAB),B5(MB2,MAB,MAB),B6(MB2,MAB,MAB)
      COMMON/BSQN/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/CLRE/RKLLLL(MB2,4),RKSSSS(MB2,4),RKSLSL(MB2,4),
     &            RJLLLL(MB2,4),RJSSSS(MB2,4),RJLLSS(MB2,4),
     &            RJSSLL(MB2,4)
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,IEQS,IPAR,ICOR,ILEV
C
C     EMPTY COUNTER ARRAYS FOR DIRECT AND EXCHANGE INTEGRALS
      DO M=1,MAXM
        DO N=1,4
          RJLLLL(M,N) = 0.0D0
          RJLLSS(M,N) = 0.0D0
          RJSSLL(M,N) = 0.0D0
          RJSSSS(M,N) = 0.0D0
          RKLLLL(M,N) = 0.0D0
          RKSLSL(M,N) = 0.0D0
          RKSSSS(M,N) = 0.0D0
        ENDDO
      ENDDO
C
C     PREPARE VALUES FOR UPCOMING CALCULATIONS
      C1 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+1)
      C3 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+3)
      C5 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+5)
      C7 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+7)
      C9 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+9)
C
      V1 = 1.0D0
      V2 = 2.0D0
      V4 = 4.0D0
      V8 = 8.0D0
      VS = 1.6D1
C
      F  = DFLOAT(2*LQNA+1)
      G  = DFLOAT(2*LQNB+1)
C
      F0G0 = 1.0D0
      F1G0 = F
      F0G1 = G
      F1G1 = F*G
      F2G0 = F*F
      F0G2 = G*G
      F2G1 = F*F*G
      F1G2 = F*G*G
      F2G2 = F*F*G*G
C
C**********************************************************************C
C     AN (LQNA,LQNB) COMBINATION HAS 1, 2 OR 4 (KQNA,KQNB) SUB-BLOCKS  C
C     SMALL-COMPONENT CONTRIBUTIONS DEPEND ON KQN SYMMETRY TYPE.       C
C**********************************************************************C
C
C     INITIATE LOOP OVER K,L BASIS FUNCTIONS
      DO M=1,MAXM
C
C       MORE VALUE PREPARATION
        E0000 = 1.0D0
        E1000 = EI
        E0100 = EJ
        E0010 = EK(M)
        E0001 = EL(M)
        E1100 = EI*EJ
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
C**********************************************************************C
C       DIRECT INTEGRAL MATRICES: RJSSSS, RJLLSS, RJSSLL, RJLLLL       C
C**********************************************************************C
C
C       TEMPORARY STORAGE OF VALUES
C       BXY MARKS 'B' BETA COMBINATION AND 'X' ONTO EFFECTIVE LQN STORE
        B00 = EIJ(2)*EKL(M,1)*B1(M,1,1) + EIJ(1)*EKL(M,2)*B2(M,1,1)
        B02 = EIJ(2)*EKL(M,3)*B1(M,1,2) + EIJ(1)*EKL(M,4)*B2(M,2,1)
        B04 = EIJ(2)*EKL(M,5)*B1(M,1,3) + EIJ(1)*EKL(M,6)*B2(M,3,1)
        B20 = EIJ(4)*EKL(M,1)*B1(M,2,1) + EIJ(3)*EKL(M,2)*B2(M,1,2)
        B22 = EIJ(4)*EKL(M,3)*B1(M,2,2) + EIJ(3)*EKL(M,4)*B2(M,2,2)
        B24 = EIJ(4)*EKL(M,5)*B1(M,2,3) + EIJ(3)*EKL(M,6)*B2(M,3,2)
        B40 = EIJ(6)*EKL(M,1)*B1(M,3,1) + EIJ(5)*EKL(M,2)*B2(M,1,3)
        B42 = EIJ(6)*EKL(M,3)*B1(M,3,2) + EIJ(5)*EKL(M,4)*B2(M,2,3)
        B44 = EIJ(6)*EKL(M,5)*B1(M,3,3) + EIJ(5)*EKL(M,6)*B2(M,3,3)
C
C       FILL RJ ARRAYS FOR THIS LQNA,LQNB BLOCK
C
C >>>   LQNA  =  0 AND LQNB  =  0 (REQUIRED FOR ALL BLOCKS)
        RJLLLL(M,4) = V1*F0G0*E0000*C5*B22
        IF(HMLTN.EQ.'NORL') GOTO 101
        RJLLSS(M,4) = V4*F0G0*E0011*C7*B24
        RJSSLL(M,4) = V4*F0G0*E1100*C7*B42
        RJSSSS(M,4) = VS*F0G0*E1111*C9*B44
C
C >>>   LQNA =/= 0                (NEED KQNA > 0 BLOCK)
        IF(LQNA.EQ.0) GOTO 103
        RJLLLL(M,3) = RJLLLL(M,4)
        RJLLSS(M,3) = RJLLSS(M,4)
        RJSSLL(M,3) = V4*F0G0*E1100*C7*B42
     &              - V2*F1G0*E1000*C5*B22 - V2*F1G0*E0100*C5*B22
     &              + V1*F2G0*E0000*C3*B02
        RJSSSS(M,3) = VS*F0G0*E1111*C9*B44
     &              - V8*F1G0*E1011*C7*B24 - V8*F1G0*E0111*C7*B24
     &              + V4*F2G0*E0011*C5*B04
103     CONTINUE
C
C >>>                  LQNB =/= 0 (NEED KQNB > 0 BLOCK)
        IF(LQNB.EQ.0) GOTO 102
        RJLLLL(M,2) = RJLLLL(M,4)
        RJLLSS(M,2) = V4*F0G0*E0011*C7*B24 
     &              - V2*F0G1*E0010*C5*B22 - V2*F0G1*E0001*C5*B22
     &              + V1*F0G2*E0000*C3*B20
        RJSSLL(M,2) = RJSSLL(M,4)
        RJSSSS(M,2) = VS*F0G0*E1111*C9*B44
     &              - V8*F0G1*E1110*C7*B42 - V8*F0G1*E1101*C7*B42
     &              + V4*F0G2*E1100*C5*B40
102     CONTINUE
C
C >>>   LQNA =/= 0 AND LQNB =/= 0 (NEED KQNA,KQNB > 0 BLOCK)
        IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 101
        RJLLLL(M,1) = RJLLLL(M,4)
        RJLLSS(M,1) = RJLLSS(M,2)
        RJSSLL(M,1) = RJSSLL(M,3)
        RJSSSS(M,1) = VS*F0G0*E1111*C9*B44
     &              - V8*F0G1*E1110*C7*B42 - V8*F0G1*E1101*C7*B42
     &              - V8*F1G0*E1011*C7*B24 - V8*F1G0*E0111*C7*B24
     &              + V4*F2G0*E0011*C5*B04 + V4*F0G2*E1100*C5*B40
     &              + V4*F1G1*E1001*C5*B22 + V4*F1G1*E0110*C5*B22
     &              + V4*F1G1*E0101*C5*B22 + V4*F1G1*E1010*C5*B22
     &              - V2*F1G2*E1000*C3*B20 - V2*F1G2*E0100*C3*B20
     &              - V2*F2G1*E0010*C3*B02 - V2*F2G1*E0001*C3*B02
     &              + V1*F2G2*E0000*C1*B00
101     CONTINUE
C
C
C**********************************************************************C
C       EXCHANGE INTEGRAL MATRICES: RKSSSS, RKSLSL, RKLLLL             C
C**********************************************************************C
C
C       LOOP OVER ORDERS FOR BETA INTEGRALS
        DO IV=1,NUNUM
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(IV)
C
C         INDEX OFFSETS FOR THE EXPONENT LISTS
          IA = NUMAX+NU+1
          IB = NUMAX-NU+1
C
C         INDEX OFFSETS FOR BETA INTEGRAL ARRAYS
          NX = (-NUMIN+NU)/2
          NY = ( NUMAX-NU)/2
C
C         TEMPORARY STORAGE OF VALUES
C         BXY MARKS 'B' BETA COMBINATION AND 'X' TO EFFECTIVE LQN STORE
C         ONTO WHICH NU IS ADDED OR SUBTRACTED.         
          B00 = EIK(M,IA+1)*EJL(M,IB  )*B3(M,NX+1,NY+1)
     &        + EIK(M,IB  )*EJL(M,IA+1)*B4(M,NX+1,NY+1)
          B02 = EIK(M,IA+1)*EJL(M,IB+2)*B3(M,NX+1,NY+2)
     &        + EIK(M,IB  )*EJL(M,IA+3)*B4(M,NX+2,NY+1)
          B04 = EIK(M,IA+1)*EJL(M,IB+4)*B3(M,NX+1,NY+3)
     &        + EIK(M,IB  )*EJL(M,IA+5)*B4(M,NX+3,NY+1)
          B20 = EIK(M,IA+3)*EJL(M,IB  )*B3(M,NX+2,NY+1)
     &        + EIK(M,IB+2)*EJL(M,IA+1)*B4(M,NX+1,NY+2)
          B22 = EIK(M,IA+3)*EJL(M,IB+2)*B3(M,NX+2,NY+2)
     &        + EIK(M,IB+2)*EJL(M,IA+3)*B4(M,NX+2,NY+2)
          B24 = EIK(M,IA+3)*EJL(M,IB+4)*B3(M,NX+2,NY+3)
     &        + EIK(M,IB+2)*EJL(M,IA+5)*B4(M,NX+3,NY+2)
          B40 = EIK(M,IA+5)*EJL(M,IB  )*B3(M,NX+3,NY+1)
     &        + EIK(M,IB+4)*EJL(M,IA+1)*B4(M,NX+1,NY+3)
          B42 = EIK(M,IA+5)*EJL(M,IB+2)*B3(M,NX+3,NY+2)
     &        + EIK(M,IB+4)*EJL(M,IA+3)*B4(M,NX+2,NY+3)
          B44 = EIK(M,IA+5)*EJL(M,IB+4)*B3(M,NX+3,NY+3)
     &        + EIK(M,IB+4)*EJL(M,IA+5)*B4(M,NX+3,NY+3)
C
C >>>     LQNA  =  0 AND LQNB  =  0 (REQUIRED FOR ALL BLOCKS)
          RKLLLL(M,4) = RKLLLL(M,4) + BK(IV,4)*V1*F0G0*E0000*C5*B22
          RKSLSL(M,4) = RKSLSL(M,4) + BK(IV,4)*V4*F0G0*E1010*C7*B42
          RKSSSS(M,4) = RKSSSS(M,4) + BK(IV,4)*VS*F0G0*E1111*C9*B44
C
C >>>     LQNA =/= 0                (NEED KQNA > 0 BLOCK)
          IF(LQNA.EQ.0) GOTO 203
          RKLL = V1*F0G0*E0000*C5*B22
          RKSL = V4*F0G0*E1010*C7*B42 - V2*F1G0*E0010*C5*B22
          RKSS = VS*F0G0*E1111*C9*B44 - V8*F1G0*E1011*C7*B42 
     &         - V8*F1G0*E0111*C7*B24 + V4*F2G0*E0011*C5*B22
          RKLLLL(M,3) = RKLLLL(M,3) + BK(IV,3)*RKLL
          RKSLSL(M,3) = RKSLSL(M,3) + BK(IV,3)*RKSL
          RKSSSS(M,3) = RKSSSS(M,3) + BK(IV,3)*RKSS
203       CONTINUE
C
C >>>                  LQNB =/= 0 (NEED KQNB > 0 BLOCK)
          IF(LQNB.EQ.0) GOTO 202
          RKLL = V1*F0G0*E0000*C5*B22
          RKSL = V4*F0G0*E1010*C7*B42 - V2*F0G1*E1000*C5*B22
          RKSS = VS*F0G0*E1111*C9*B44 - V8*F0G1*E1110*C7*B42 
     &         - V8*F0G1*E1101*C7*B24 + V4*F0G2*E1100*C5*B22
          RKLLLL(M,2) = RKLLLL(M,2) + BK(IV,2)*RKLL
          RKSLSL(M,2) = RKSLSL(M,2) + BK(IV,2)*RKSL
          RKSSSS(M,2) = RKSSSS(M,2) + BK(IV,2)*RKSS
202       CONTINUE
C
C >>>     LQNA =/= 0 AND LQNB =/= 0 (NEED KQNA,KQNB > 0 BLOCK)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 201
          RKLL = V1*F0G0*E0000*C5*B22
          RKSL = V4*F0G0*E1010*C7*B42 - V2*F1G0*E0010*C5*B22 
     &         - V2*F0G1*E1000*C5*B22 + V1*F1G1*E0000*C3*B02
          RKSS = VS*F0G0*E1111*C9*B44
     &         - V8*F0G1*E1110*C7*B42 - V8*F0G1*E1101*C7*B24 
     &         - V8*F1G0*E1011*C7*B42 - V8*F1G0*E0111*C7*B24
     &         + V4*F2G0*E0011*C5*B22 + V4*F0G2*E1100*C5*B22 
     &         + V4*F1G1*E0110*C5*B22 + V4*F1G1*E1001*C5*B22
     &         + V4*F1G1*E1010*C5*B40 + V4*F1G1*E0101*C5*B04
     &         - V2*F2G1*E0010*C3*B20 - V2*F1G2*E1000*C3*B20
     &         - V2*F2G1*E0001*C3*B02 - V2*F1G2*E0100*C3*B02
     &         + V1*F2G2*E0000*C1*B00
          RKLLLL(M,1) = RKLLLL(M,1) + BK(IV,1)*RKLL
          RKSLSL(M,1) = RKSLSL(M,1) + BK(IV,1)*RKSL
          RKSSSS(M,1) = RKSSSS(M,1) + BK(IV,1)*RKSS 
201       CONTINUE
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     NORMALISE ACCUMULATED INTEGRALS (DIRECT AND EXCHANGE)            C
C**********************************************************************C
C
      DO M=1,MAXM
        T0LLLL = RNIJ(1)*RNKL(M,1)
        T0LLSS = RNIJ(1)*RNKL(M,3)
        T0SLSL = RNIJ(2)*RNKL(M,2)
        T0SSLL = RNIJ(3)*RNKL(M,1)
        T0SSSS = RNIJ(3)*RNKL(M,3)
        DO N=1,4
          RJLLLL(M,N) = RJLLLL(M,N)*T0LLLL
          RJLLSS(M,N) = RJLLSS(M,N)*T0LLSS
          RJSSLL(M,N) = RJSSLL(M,N)*T0SSLL
          RJSSSS(M,N) = RJSSSS(M,N)*T0SSSS
          RKLLLL(M,N) = RKLLLL(M,N)*T0LLLL
          RKSLSL(M,N) = RKSLSL(M,N)*T0SLSL
          RKSSSS(M,N) = RKSSSS(M,N)*T0SSSS
        ENDDO
      ENDDO
C
      RETURN
      END

