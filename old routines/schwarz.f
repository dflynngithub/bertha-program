      SUBROUTINE SCHWARZ(ICNT,KQN,MQN,NBAS,ITN,IFLG,ISCR,ISKP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   SSSSSS   CCCCCC  HH    HH WW         WW    AA    RRRRRRR  ZZZZZZZZ C
C  SS    SS CC    CC HH    HH WW         WW   AAAA   RR    RR      ZZ  C
C  SS       CC       HH    HH WW         WW  AA  AA  RR    RR     ZZ   C
C   SSSSSS  CC       HHHHHHHH WW    W    WW AA    AA RR    RR    ZZ    C
C        SS CC       HH    HH  WW  WWW  WW  AAAAAAAA RRRRRRR    ZZ     C
C  SS    SS CC    CC HH    HH   WWWW WWWW   AA    AA RR    RR  ZZ      C
C   SSSSSS   CCCCCC  HH    HH    WW   WW    AA    AA RR    RR ZZZZZZZZ C
C                                                                      C
C -------------------------------------------------------------------- C
C  SCHWARZ APPROXIMATES THE UPPER BOUND TO A BLOCK OF TWO-ELECTRON     C
C  INTEGRALS AND ITS CONTRIBUTION TO THE GDIR/GXCH MATRIX.             C
C -------------------------------------------------------------------- C
C  DFNOTE: ELABORATE FLAG CONSIDERATIONS SEEM TO COST MORE TIME THAN   C
C  IT SAVES LATER. MIGHT BE BETTER TO SEARCH FOR HIGHEST CONTRIBUTION  C
C  FROM ALL DENC OVERLAPS (NA1 -> ND2 POSSIBILITIES) AND ALL GDSC      C
C  CHOICES, THEN TO LOOK AT THE PRODUCT. (NEED DENO CHOICES AS WELL.)  C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=15,MKP=9,MFL=7000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      CHARACTER*4 HMLTN
C
      DIMENSION ICNT(4),KQN(4),MQN(4),NBAS(4),ITN(2),KNDX(4)
      DIMENSION IFLG(11),ISCR(11)
      DIMENSION GDSC(MDM,MDM),BDSC(MDM,MDM)
      DIMENSION RR(MB2,16),X(16)
C
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
C
      COMMON/BLOC/PAB1,PAB2,PCD1,PCD2,NA1,NB1,NC1,ND1,NA2,NB2,NC2,ND2,
     &            IBAS,JBAS,ILIN
      COMMON/DENS/DENC,DENO,DENT
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,IEQS,IERC,IPAR,ICOR,ILEV
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
      COMMON/SWRZ/GDSC,BDSC
C
      SENS = 1.0D-09
C
C     DATA FOR ARRAY ADDRESSES
      IF(ITN(1).EQ.1) THEN
        NADDAB = 0
      ELSE
        NADDAB = NSHIFT
      ENDIF
C
      IF(ITN(2).EQ.1) THEN
        NADDCD = 0
      ELSE
        NADDCD = NSHIFT
      ENDIF
C
      DO N=1,4
        KNDX(N) = 2*IABS(KQN(N))
        IF(KQN(N).LT.0) THEN
          KNDX(N) = KNDX(N)-1
        ENDIF
      ENDDO
C
C     FOCK ADDRESS FOR EACH BASIS FUNCTION (WITH SPIN PROJECTION)
      NA1 = LARGE(ICNT(1),KNDX(1),MQN(1)  ) + NADDAB
      NB1 = LARGE(ICNT(2),KNDX(2),MQN(2)  ) + NADDAB
      NC1 = LARGE(ICNT(3),KNDX(3),MQN(3)  ) + NADDCD
      ND1 = LARGE(ICNT(4),KNDX(4),MQN(4)  ) + NADDCD
C
      NA2 = LARGE(ICNT(1),KNDX(1),MQN(1)+1) + NADDAB
      NB2 = LARGE(ICNT(2),KNDX(2),MQN(2)+1) + NADDAB
      NC2 = LARGE(ICNT(3),KNDX(3),MQN(3)+1) + NADDCD
      ND2 = LARGE(ICNT(4),KNDX(4),MQN(4)+1) + NADDCD
C
C     CONSTRUCT INTEGRAL UPPER BOUNDS
C      M = 0
C      DO KBAS=1,NBAS(3)
C        DO LBAS=1,NBAS(4)
C          M = M+1
CC
C          RR1(M, 1) = GDSC(NA1+IBAS,NB1+JBAS)*GDSC(NC1+KBAS,ND1+LBAS)
C          RR1(M, 2) = GDSC(NA1+IBAS,NB1+JBAS)*GDSC(NC1+KBAS,ND2+LBAS)
C          RR1(M, 3) = GDSC(NA1+IBAS,NB1+JBAS)*GDSC(NC2+KBAS,ND1+LBAS)
C          RR1(M, 4) = GDSC(NA1+IBAS,NB1+JBAS)*GDSC(NC2+KBAS,ND2+LBAS)
C          RR1(M, 5) = GDSC(NA1+IBAS,NB2+JBAS)*GDSC(NC1+KBAS,ND1+LBAS)
C          RR1(M, 6) = GDSC(NA1+IBAS,NB2+JBAS)*GDSC(NC1+KBAS,ND2+LBAS)
C          RR1(M, 7) = GDSC(NA1+IBAS,NB2+JBAS)*GDSC(NC2+KBAS,ND1+LBAS)
C          RR1(M, 8) = GDSC(NA1+IBAS,NB2+JBAS)*GDSC(NC2+KBAS,ND2+LBAS)
C          RR1(M, 9) = GDSC(NA2+IBAS,NB1+JBAS)*GDSC(NC1+KBAS,ND1+LBAS)
C          RR1(M,10) = GDSC(NA2+IBAS,NB1+JBAS)*GDSC(NC1+KBAS,ND2+LBAS)
C          RR1(M,11) = GDSC(NA2+IBAS,NB1+JBAS)*GDSC(NC2+KBAS,ND1+LBAS)
C          RR1(M,12) = GDSC(NA2+IBAS,NB1+JBAS)*GDSC(NC2+KBAS,ND2+LBAS)
C          RR1(M,13) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(NC1+KBAS,ND1+LBAS)
C          RR1(M,14) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(NC1+KBAS,ND2+LBAS)
C          RR1(M,15) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(NC2+KBAS,ND1+LBAS)
C          RR1(M,16) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(NC2+KBAS,ND2+LBAS)
CC
C          RR2(M, 1) = GDSC(NA1+IBAS,ND1+LBAS)*GDSC(NC1+KBAS,NB1+JBAS)
C          RR2(M, 2) = GDSC(NA1+IBAS,ND1+LBAS)*GDSC(NC1+KBAS,NB2+JBAS)
C          RR2(M, 3) = GDSC(NA1+IBAS,ND1+LBAS)*GDSC(NC2+KBAS,NB1+JBAS)
C          RR2(M, 4) = GDSC(NA1+IBAS,ND1+LBAS)*GDSC(NC2+KBAS,NB2+JBAS)
C          RR2(M, 5) = GDSC(NA1+IBAS,ND2+LBAS)*GDSC(NC1+KBAS,NB1+JBAS)
C          RR2(M, 6) = GDSC(NA1+IBAS,ND2+LBAS)*GDSC(NC1+KBAS,NB2+JBAS)
C          RR2(M, 7) = GDSC(NA1+IBAS,ND2+LBAS)*GDSC(NC2+KBAS,NB1+JBAS)
C          RR2(M, 8) = GDSC(NA1+IBAS,ND2+LBAS)*GDSC(NC2+KBAS,NB2+JBAS)
C          RR2(M, 9) = GDSC(NA2+IBAS,ND1+LBAS)*GDSC(NC1+KBAS,NB1+JBAS)
C          RR2(M,10) = GDSC(NA2+IBAS,ND1+LBAS)*GDSC(NC1+KBAS,NB2+JBAS)
C          RR2(M,11) = GDSC(NA2+IBAS,ND1+LBAS)*GDSC(NC2+KBAS,NB1+JBAS)
C          RR2(M,12) = GDSC(NA2+IBAS,ND1+LBAS)*GDSC(NC2+KBAS,NB2+JBAS)
C          RR2(M,13) = GDSC(NA2+IBAS,ND2+LBAS)*GDSC(NC1+KBAS,NB1+JBAS)
C          RR2(M,14) = GDSC(NA2+IBAS,ND2+LBAS)*GDSC(NC1+KBAS,NB2+JBAS)
C          RR2(M,15) = GDSC(NA2+IBAS,ND2+LBAS)*GDSC(NC2+KBAS,NB1+JBAS)
C          RR2(M,16) = GDSC(NA2+IBAS,ND2+LBAS)*GDSC(NC2+KBAS,NB2+JBAS)
CC
C          RR3(M, 1) = GDSC(NA1+IBAS,NC1+KBAS)*GDSC(ND1+LBAS,NB1+JBAS)
C          RR3(M, 2) = GDSC(NA1+IBAS,NC1+KBAS)*GDSC(ND1+LBAS,NB2+JBAS)
C          RR3(M, 3) = GDSC(NA1+IBAS,NC1+KBAS)*GDSC(ND2+LBAS,NB1+JBAS)
C          RR3(M, 4) = GDSC(NA1+IBAS,NC1+KBAS)*GDSC(ND2+LBAS,NB2+JBAS)
C          RR3(M, 5) = GDSC(NA1+IBAS,NC2+KBAS)*GDSC(ND1+LBAS,NB1+JBAS)
C          RR3(M, 6) = GDSC(NA1+IBAS,NC2+KBAS)*GDSC(ND1+LBAS,NB2+JBAS)
C          RR3(M, 7) = GDSC(NA1+IBAS,NC2+KBAS)*GDSC(ND2+LBAS,NB1+JBAS)
C          RR3(M, 8) = GDSC(NA1+IBAS,NC2+KBAS)*GDSC(ND2+LBAS,NB2+JBAS)
C          RR3(M, 9) = GDSC(NA2+IBAS,NC1+KBAS)*GDSC(ND1+LBAS,NB1+JBAS)
C          RR3(M,10) = GDSC(NA2+IBAS,NC1+KBAS)*GDSC(ND1+LBAS,NB2+JBAS)
C          RR3(M,11) = GDSC(NA2+IBAS,NC1+KBAS)*GDSC(ND2+LBAS,NB1+JBAS)
C          RR3(M,12) = GDSC(NA2+IBAS,NC1+KBAS)*GDSC(ND2+LBAS,NB2+JBAS)
C          RR3(M,13) = GDSC(NA2+IBAS,NC2+KBAS)*GDSC(ND1+LBAS,NB1+JBAS)
C          RR3(M,14) = GDSC(NA2+IBAS,NC2+KBAS)*GDSC(ND1+LBAS,NB2+JBAS)
C          RR3(M,15) = GDSC(NA2+IBAS,NC2+KBAS)*GDSC(ND2+LBAS,NB1+JBAS)
C          RR3(M,16) = GDSC(NA2+IBAS,NC2+KBAS)*GDSC(ND2+LBAS,NB2+JBAS)
CC
C        ENDDO
C      ENDDO
C
C**********************************************************************C
C     FOR EACH IFLG OPTION, FIND BIGGEST CONTRIBUTION TO FOCK MATRIX.  C
C**********************************************************************C
C
      ITAL = 0
C
C     1ST FLAG (DIRECT):  ( MA, MB| MC, MD)
      IF(IFLG(1).EQ.1) THEN
C
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
C            X( 1) = RR1(M, 1)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X( 2) = RR1(M, 2)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X( 3) = RR1(M, 3)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X( 4) = RR1(M, 4)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X( 5) = RR1(M, 5)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X( 6) = RR1(M, 6)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X( 7) = RR1(M, 7)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X( 8) = RR1(M, 8)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X( 9) = RR1(M, 9)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X(10) = RR1(M,10)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X(11) = RR1(M,11)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X(12) = RR1(M,12)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X(13) = RR1(M,13)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X(14) = RR1(M,14)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X(15) = RR1(M,15)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X(16) = RR1(M,16)*ABS(DENC(NC2+KBAS,ND2+LBAS))
C
            X( 1) = GDSC(NA1+IBAS,NB1+JBAS)*GDSC(NC1+KBAS,ND1+LBAS)
     &                                 *ABS(DENC(NC1+KBAS,ND1+LBAS))
            X( 2) = GDSC(NA1+IBAS,NB1+JBAS)*GDSC(NC1+KBAS,ND2+LBAS)
     &                                 *ABS(DENC(NC1+KBAS,ND2+LBAS))
            X( 3) = GDSC(NA1+IBAS,NB1+JBAS)*GDSC(NC2+KBAS,ND1+LBAS)
     &                                 *ABS(DENC(NC2+KBAS,ND1+LBAS))
            X( 4) = GDSC(NA1+IBAS,NB1+JBAS)*GDSC(NC2+KBAS,ND2+LBAS)
     &                                 *ABS(DENC(NC2+KBAS,ND2+LBAS))
C
            X( 5) = GDSC(NA1+IBAS,NB2+JBAS)*GDSC(NC1+KBAS,ND1+LBAS)
     &                                 *ABS(DENC(NC1+KBAS,ND1+LBAS))
            X( 6) = GDSC(NA1+IBAS,NB2+JBAS)*GDSC(NC1+KBAS,ND2+LBAS)
     &                                 *ABS(DENC(NC1+KBAS,ND2+LBAS))
            X( 7) = GDSC(NA1+IBAS,NB2+JBAS)*GDSC(NC2+KBAS,ND1+LBAS)
     &                                 *ABS(DENC(NC2+KBAS,ND1+LBAS))
            X( 8) = GDSC(NA1+IBAS,NB2+JBAS)*GDSC(NC2+KBAS,ND2+LBAS)
     &                                 *ABS(DENC(NC2+KBAS,ND2+LBAS))
C
            X( 9) = GDSC(NA2+IBAS,NB1+JBAS)*GDSC(NC1+KBAS,ND1+LBAS)
     &                                 *ABS(DENC(NC1+KBAS,ND1+LBAS))
            X(10) = GDSC(NA2+IBAS,NB1+JBAS)*GDSC(NC1+KBAS,ND2+LBAS)
     &                                 *ABS(DENC(NC1+KBAS,ND2+LBAS))
            X(11) = GDSC(NA2+IBAS,NB1+JBAS)*GDSC(NC2+KBAS,ND1+LBAS)
     &                                 *ABS(DENC(NC2+KBAS,ND1+LBAS))
            X(12) = GDSC(NA2+IBAS,NB1+JBAS)*GDSC(NC2+KBAS,ND2+LBAS)
     &                                 *ABS(DENC(NC2+KBAS,ND2+LBAS))
C
            X(13) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(NC1+KBAS,ND1+LBAS)
     &                                 *ABS(DENC(NC1+KBAS,ND1+LBAS))
            X(14) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(NC1+KBAS,ND2+LBAS)
     &                                 *ABS(DENC(NC1+KBAS,ND2+LBAS))
            X(15) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(NC2+KBAS,ND1+LBAS)
     &                                 *ABS(DENC(NC2+KBAS,ND1+LBAS))
            X(16) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(NC2+KBAS,ND2+LBAS)
     &                                 *ABS(DENC(NC2+KBAS,ND2+LBAS))
C
C           EARLY EXIT IF ANY OF THESE LARGER THAN SENS
            DO N=1,16
              IF(DABS(X(N)).GT.SENS) THEN
                ISCR(1) = 1
                ITAL    = ITAL+1
                GOTO 101
              ENDIF
            ENDDO
C
          ENDDO
        ENDDO
C
C       LOOP COMPLETED WITHOUT EARLY EXIT
        ISCR(1) = 0
C
C       EXIT IFLG(1) CHECK
101     CONTINUE
C
      ENDIF
C
C     2ND FLAG (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(2).EQ.1) THEN
C
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
CC
C            X( 1) = RR1(M, 4)*ABS(DENC(ND1+LBAS,NC1+KBAS))
C            X( 2) = RR1(M, 2)*ABS(DENC(ND1+LBAS,NC2+KBAS))
C            X( 3) = RR1(M, 3)*ABS(DENC(ND2+LBAS,NC1+KBAS))
C            X( 4) = RR1(M, 1)*ABS(DENC(ND2+LBAS,NC2+KBAS))
CC
C            X( 5) = RR1(M, 8)*ABS(DENC(ND1+LBAS,NC1+KBAS))
C            X( 6) = RR1(M, 6)*ABS(DENC(ND1+LBAS,NC2+KBAS))
C            X( 7) = RR1(M, 7)*ABS(DENC(ND2+LBAS,NC1+KBAS))
C            X( 8) = RR1(M, 5)*ABS(DENC(ND2+LBAS,NC2+KBAS))
CC
C            X( 9) = RR1(M,12)*ABS(DENC(ND1+LBAS,NC1+KBAS))
C            X(10) = RR1(M,10)*ABS(DENC(ND1+LBAS,NC2+KBAS))
C            X(11) = RR1(M,11)*ABS(DENC(ND2+LBAS,NC1+KBAS))
C            X(12) = RR1(M, 9)*ABS(DENC(ND1+LBAS,NC1+KBAS))
CC
C            X(13) = RR1(M,16)*ABS(DENC(ND1+LBAS,NC1+KBAS))
C            X(14) = RR1(M,14)*ABS(DENC(ND1+LBAS,NC2+KBAS))
C            X(15) = RR1(M,15)*ABS(DENC(ND2+LBAS,NC1+KBAS))
C            X(16) = RR1(M,13)*ABS(DENC(ND1+LBAS,NC1+KBAS))
CC
            X( 1) = GDSC(NA1+IBAS,NB1+JBAS)*GDSC(ND2+LBAS,NC2+KBAS)
     &                                 *ABS(DENC(ND1+LBAS,NC1+KBAS))
            X( 2) = GDSC(NA1+IBAS,NB1+JBAS)*GDSC(ND2+LBAS,NC1+KBAS)
     &                                 *ABS(DENC(ND1+LBAS,NC2+KBAS))
            X( 3) = GDSC(NA1+IBAS,NB1+JBAS)*GDSC(ND1+LBAS,NC2+KBAS)
     &                                 *ABS(DENC(ND2+LBAS,NC1+KBAS))
            X( 4) = GDSC(NA1+IBAS,NB1+JBAS)*GDSC(ND1+LBAS,NC1+KBAS)
     &                                 *ABS(DENC(ND2+LBAS,NC2+KBAS))
C
            X( 5) = GDSC(NA1+IBAS,NB2+JBAS)*GDSC(ND2+LBAS,NC2+KBAS)
     &                                 *ABS(DENC(ND1+LBAS,NC1+KBAS))
            X( 6) = GDSC(NA1+IBAS,NB2+JBAS)*GDSC(ND2+LBAS,NC1+KBAS)
     &                                 *ABS(DENC(ND1+LBAS,NC2+KBAS))
            X( 7) = GDSC(NA1+IBAS,NB2+JBAS)*GDSC(ND1+LBAS,NC2+KBAS)
     &                                 *ABS(DENC(ND2+LBAS,NC1+KBAS))
            X( 8) = GDSC(NA1+IBAS,NB2+JBAS)*GDSC(ND1+LBAS,NC1+KBAS)
     &                                 *ABS(DENC(ND2+LBAS,NC2+KBAS))
C
            X( 9) = GDSC(NA2+IBAS,NB1+JBAS)*GDSC(ND2+LBAS,NC2+KBAS)
     &                                 *ABS(DENC(ND1+LBAS,NC1+KBAS))
            X(10) = GDSC(NA2+IBAS,NB1+JBAS)*GDSC(ND2+LBAS,NC1+KBAS)
     &                                 *ABS(DENC(ND1+LBAS,NC2+KBAS))
            X(11) = GDSC(NA2+IBAS,NB1+JBAS)*GDSC(ND1+LBAS,NC2+KBAS)
     &                                 *ABS(DENC(ND2+LBAS,NC1+KBAS))
            X(12) = GDSC(NA2+IBAS,NB1+JBAS)*GDSC(ND2+LBAS,NC2+KBAS)
     &                                 *ABS(DENC(ND1+LBAS,NC1+KBAS))
C
            X(13) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(ND2+LBAS,NC2+KBAS)
     &                                 *ABS(DENC(ND1+LBAS,NC1+KBAS))
            X(14) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(ND2+LBAS,NC1+KBAS)
     &                                 *ABS(DENC(ND1+LBAS,NC2+KBAS))
            X(15) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(ND1+LBAS,NC2+KBAS)
     &                                 *ABS(DENC(ND2+LBAS,NC1+KBAS))
            X(16) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(ND2+LBAS,NC2+KBAS)
     &                                 *ABS(DENC(ND1+LBAS,NC1+KBAS))
C
C           EARLY EXIT IF ANY OF THESE LARGER THAN SENS
            DO N=1,16
              IF(DABS(X(N)).GT.SENS) THEN
                ISCR(2) = 1
                ITAL    = ITAL+1
                GOTO 102
              ENDIF
            ENDDO
C
          ENDDO
        ENDDO
C
C       LOOP COMPLETED WITHOUT EARLY EXIT
        ISCR(2) = 0
C
C       EXIT IFLG(2) CHECK
102     CONTINUE
C
      ENDIF
C
C     3RD FLAG (DIRECT):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(3).EQ.1) THEN
C
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
CC
C            X( 1) = RR1(M, 1)*ABS(DENC(NA1+IBAS,NB1+JBAS))
C            X( 2) = RR1(M, 5)*ABS(DENC(NA1+IBAS,NB2+JBAS))
C            X( 3) = RR1(M, 9)*ABS(DENC(NA2+IBAS,NB1+JBAS))
C            X( 4) = RR1(M,13)*ABS(DENC(NA2+IBAS,NB2+JBAS))
CC
C            X( 5) = RR1(M, 2)*ABS(DENC(NA1+IBAS,NB1+JBAS))
C            X( 6) = RR1(M, 6)*ABS(DENC(NA1+IBAS,NB2+JBAS))
C            X( 7) = RR1(M,10)*ABS(DENC(NA2+IBAS,NB1+JBAS))
C            X( 8) = RR1(M,14)*ABS(DENC(NA2+IBAS,NB2+JBAS))
CC
C            X( 9) = RR1(M, 3)*ABS(DENC(NA1+IBAS,NB1+JBAS))
C            X(10) = RR1(M, 7)*ABS(DENC(NA1+IBAS,NB2+JBAS))
C            X(11) = RR1(M,11)*ABS(DENC(NA2+IBAS,NB1+JBAS))
C            X(12) = RR1(M,15)*ABS(DENC(NA2+IBAS,NB2+JBAS))
CC
C            X(13) = RR1(M, 4)*ABS(DENC(NA1+IBAS,NB1+JBAS))
C            X(14) = RR1(M, 8)*ABS(DENC(NA1+IBAS,NB2+JBAS))
C            X(15) = RR1(M,12)*ABS(DENC(NA2+IBAS,NB1+JBAS))
C            X(16) = RR1(M,16)*ABS(DENC(NA2+IBAS,NB2+JBAS))
CC
            X( 1) = GDSC(NC1+KBAS,ND1+LBAS)*GDSC(NA1+IBAS,NB1+JBAS)
     &                                 *ABS(DENC(NA1+IBAS,NB1+JBAS))
            X( 2) = GDSC(NC1+KBAS,ND1+LBAS)*GDSC(NA1+IBAS,NB2+JBAS)
     &                                 *ABS(DENC(NA1+IBAS,NB2+JBAS))
            X( 3) = GDSC(NC1+KBAS,ND1+LBAS)*GDSC(NA2+IBAS,NB1+JBAS)
     &                                 *ABS(DENC(NA2+IBAS,NB1+JBAS))
            X( 4) = GDSC(NC1+KBAS,ND1+LBAS)*GDSC(NA2+IBAS,NB2+JBAS)
     &                                 *ABS(DENC(NA2+IBAS,NB2+JBAS))
C
            X( 5) = GDSC(NC1+KBAS,ND2+LBAS)*GDSC(NA1+IBAS,NB1+JBAS)
     &                                 *ABS(DENC(NA1+IBAS,NB1+JBAS))
            X( 6) = GDSC(NC1+KBAS,ND2+LBAS)*GDSC(NA1+IBAS,NB2+JBAS)
     &                                 *ABS(DENC(NA1+IBAS,NB2+JBAS))
            X( 7) = GDSC(NC1+KBAS,ND2+LBAS)*GDSC(NA2+IBAS,NB1+JBAS)
     &                                 *ABS(DENC(NA2+IBAS,NB1+JBAS))
            X( 8) = GDSC(NC1+KBAS,ND2+LBAS)*GDSC(NA2+IBAS,NB2+JBAS)
     &                                 *ABS(DENC(NA2+IBAS,NB2+JBAS))
C
            X( 9) = GDSC(NC2+KBAS,ND1+LBAS)*GDSC(NA1+IBAS,NB1+JBAS)
     &                                 *ABS(DENC(NA1+IBAS,NB1+JBAS))
            X(10) = GDSC(NC2+KBAS,ND1+LBAS)*GDSC(NA1+IBAS,NB2+JBAS)
     &                                 *ABS(DENC(NA1+IBAS,NB2+JBAS))
            X(11) = GDSC(NC2+KBAS,ND1+LBAS)*GDSC(NA2+IBAS,NB1+JBAS)
     &                                 *ABS(DENC(NA2+IBAS,NB1+JBAS))
            X(12) = GDSC(NC2+KBAS,ND1+LBAS)*GDSC(NA2+IBAS,NB2+JBAS)
     &                                 *ABS(DENC(NA2+IBAS,NB2+JBAS))
C
            X(13) = GDSC(NC2+KBAS,ND2+LBAS)*GDSC(NA1+IBAS,NB1+JBAS)
     &                                 *ABS(DENC(NA1+IBAS,NB1+JBAS))
            X(14) = GDSC(NC2+KBAS,ND2+LBAS)*GDSC(NA1+IBAS,NB2+JBAS)
     &                                 *ABS(DENC(NA1+IBAS,NB2+JBAS))
            X(15) = GDSC(NC2+KBAS,ND2+LBAS)*GDSC(NA2+IBAS,NB1+JBAS)
     &                                 *ABS(DENC(NA2+IBAS,NB1+JBAS))
            X(16) = GDSC(NC2+KBAS,ND2+LBAS)*GDSC(NA2+IBAS,NB2+JBAS)
     &                                 *ABS(DENC(NA2+IBAS,NB2+JBAS))
C
C           EARLY EXIT IF ANY OF THESE LARGER THAN SENS
            DO N=1,16
              IF(DABS(X(N)).GT.SENS) THEN
                ISCR(3) = 1
                ITAL    = ITAL+1
                GOTO 103
              ENDIF
            ENDDO
C
          ENDDO
        ENDDO
C
C       LOOP COMPLETED WITHOUT EARLY EXIT
        ISCR(3) = 0
C
C       EXIT IFLG(3) CHECK
103     CONTINUE
C
      ENDIF
C
C     4TH FLAG (DIRECT):  ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(4).EQ.1) THEN
C
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
CC
C            X( 1) = RR1(M,13)*ABS(DENC(NB1+JBAS,NA1+IBAS))
C            X( 2) = RR1(M, 5)*ABS(DENC(NB1+JBAS,NA2+IBAS))
C            X( 3) = RR1(M, 9)*ABS(DENC(NB2+JBAS,NA1+IBAS))
C            X( 4) = RR1(M, 1)*ABS(DENC(NB2+JBAS,NA2+IBAS))
CC
C            X( 5) = RR1(M,14)*ABS(DENC(NB1+JBAS,NA1+IBAS))
C            X( 6) = RR1(M, 6)*ABS(DENC(NB1+JBAS,NA2+IBAS))
C            X( 7) = RR1(M,10)*ABS(DENC(NB2+JBAS,NA1+IBAS))
C            X( 8) = RR1(M, 2)*ABS(DENC(NB2+JBAS,NA2+IBAS))
CC
C            X( 9) = RR1(M,15)*ABS(DENC(NB1+JBAS,NA1+IBAS))
C            X(10) = RR1(M, 7)*ABS(DENC(NB1+JBAS,NA2+IBAS))
C            X(11) = RR1(M,11)*ABS(DENC(NB2+JBAS,NA1+IBAS))
C            X(12) = RR1(M, 3)*ABS(DENC(NB2+JBAS,NA2+IBAS))
CC
C            X(13) = RR1(M,16)*ABS(DENC(NB1+JBAS,NA1+IBAS))
C            X(14) = RR1(M, 8)*ABS(DENC(NB1+JBAS,NA2+IBAS))
C            X(15) = RR1(M,12)*ABS(DENC(NB2+JBAS,NA1+IBAS))
C            X(16) = RR1(M, 4)*ABS(DENC(NB2+JBAS,NA2+IBAS))
CC
            X( 1) = GDSC(NC1+KBAS,ND1+LBAS)*GDSC(NB1+JBAS,NA1+IBAS)
     &                                 *ABS(DENC(NB1+JBAS,NA1+IBAS))
            X( 2) = GDSC(NC1+KBAS,ND1+LBAS)*GDSC(NB1+JBAS,NA2+IBAS)
     &                                 *ABS(DENC(NB1+JBAS,NA2+IBAS))
            X( 3) = GDSC(NC1+KBAS,ND1+LBAS)*GDSC(NB2+JBAS,NA1+IBAS)
     &                                 *ABS(DENC(NB2+JBAS,NA1+IBAS))
            X( 4) = GDSC(NC1+KBAS,ND1+LBAS)*GDSC(NB2+JBAS,NA2+IBAS)
     &                                 *ABS(DENC(NB2+JBAS,NA2+IBAS))
C
            X( 5) = GDSC(NC1+KBAS,ND2+LBAS)*GDSC(NB1+JBAS,NA1+IBAS)
     &                                 *ABS(DENC(NB1+JBAS,NA1+IBAS))
            X( 6) = GDSC(NC1+KBAS,ND2+LBAS)*GDSC(NB1+JBAS,NA2+IBAS)
     &                                 *ABS(DENC(NB1+JBAS,NA2+IBAS))
            X( 7) = GDSC(NC1+KBAS,ND2+LBAS)*GDSC(NB2+JBAS,NA1+IBAS)
     &                                 *ABS(DENC(NB2+JBAS,NA1+IBAS))
            X( 8) = GDSC(NC1+KBAS,ND2+LBAS)*GDSC(NB2+JBAS,NA2+IBAS)
     &                                 *ABS(DENC(NB2+JBAS,NA2+IBAS))
C
            X( 9) = GDSC(NC2+KBAS,ND1+LBAS)*GDSC(NB1+JBAS,NA1+IBAS)
     &                                 *ABS(DENC(NB1+JBAS,NA1+IBAS))
            X(10) = GDSC(NC2+KBAS,ND1+LBAS)*GDSC(NB1+JBAS,NA2+IBAS)
     &                                 *ABS(DENC(NB1+JBAS,NA2+IBAS))
            X(11) = GDSC(NC2+KBAS,ND1+LBAS)*GDSC(NB2+JBAS,NA1+IBAS)
     &                                 *ABS(DENC(NB2+JBAS,NA1+IBAS))
            X(12) = GDSC(NC2+KBAS,ND1+LBAS)*GDSC(NB2+JBAS,NA2+IBAS)
     &                                 *ABS(DENC(NB2+JBAS,NA2+IBAS))
c
            X(13) = GDSC(NC2+KBAS,ND2+LBAS)*GDSC(NB1+JBAS,NA1+IBAS)
     &                                 *ABS(DENC(NB1+JBAS,NA1+IBAS))
            X(14) = GDSC(NC2+KBAS,ND2+LBAS)*GDSC(NB1+JBAS,NA2+IBAS)
     &                                 *ABS(DENC(NB1+JBAS,NA2+IBAS))
            X(15) = GDSC(NC2+KBAS,ND2+LBAS)*GDSC(NB2+JBAS,NA1+IBAS)
     &                                 *ABS(DENC(NB2+JBAS,NA1+IBAS))
            X(16) = GDSC(NC2+KBAS,ND2+LBAS)*GDSC(NB2+JBAS,NA2+IBAS)
     &                                 *ABS(DENC(NB2+JBAS,NA2+IBAS))
C
C           EARLY EXIT IF ANY OF THESE LARGER THAN SENS
            DO N=1,16
              IF(DABS(X(N)).GT.SENS) THEN
                ISCR(4) = 1
                ITAL    = ITAL+1
                GOTO 104
              ENDIF
            ENDDO
C
          ENDDO
        ENDDO
C
C       LOOP COMPLETED WITHOUT EARLY EXIT
        ISCR(4) = 0
C
C       EXIT IFLG(4) CHECK
104     CONTINUE
C
      ENDIF
C
C     5TH FLAG (EXCHNG):  ( MA, MB| MC, MD)
      IF(IFLG(5).EQ.1) THEN
C
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
C
C            X( 1) = RR2(M, 1)*ABS(DENC(NC1+KBAS,NB1+JBAS))
C            X( 2) = RR2(M, 2)*ABS(DENC(NC1+KBAS,NB2+JBAS))
C            X( 3) = RR2(M, 3)*ABS(DENC(NC2+KBAS,NB1+JBAS))
C            X( 4) = RR2(M, 4)*ABS(DENC(NC2+KBAS,NB2+JBAS))
CC
C            X( 5) = RR2(M, 5)*ABS(DENC(NC1+KBAS,NB1+JBAS))
C            X( 6) = RR2(M, 6)*ABS(DENC(NC1+KBAS,NB2+JBAS))
C            X( 7) = RR2(M, 7)*ABS(DENC(NC2+KBAS,NB1+JBAS))
C            X( 8) = RR2(M, 8)*ABS(DENC(NC2+KBAS,NB2+JBAS))
CC
C            X( 9) = RR2(M, 9)*ABS(DENC(NC1+KBAS,NB1+JBAS))
C            X(10) = RR2(M,10)*ABS(DENC(NC1+KBAS,NB2+JBAS))
C            X(11) = RR2(M,11)*ABS(DENC(NC2+KBAS,NB1+JBAS))
C            X(12) = RR2(M,12)*ABS(DENC(NC2+KBAS,NB2+JBAS))
CC
C            X(13) = RR2(M,13)*ABS(DENC(NC1+KBAS,NB1+JBAS))
C            X(14) = RR2(M,14)*ABS(DENC(NC1+KBAS,NB2+JBAS))
C            X(15) = RR2(M,15)*ABS(DENC(NC2+KBAS,NB1+JBAS))
C            X(16) = RR2(M,16)*ABS(DENC(NC2+KBAS,NB2+JBAS))

C
            X( 1) = GDSC(NA1+IBAS,ND1+LBAS)*GDSC(NC1+KBAS,NB1+JBAS)
     &                                 *ABS(DENC(NC1+KBAS,NB1+JBAS))
            X( 2) = GDSC(NA1+IBAS,ND1+LBAS)*GDSC(NC1+KBAS,NB2+JBAS)
     &                                 *ABS(DENC(NC1+KBAS,NB2+JBAS))
            X( 3) = GDSC(NA1+IBAS,ND1+LBAS)*GDSC(NC2+KBAS,NB1+JBAS)
     &                                 *ABS(DENC(NC2+KBAS,NB1+JBAS))
            X( 4) = GDSC(NA1+IBAS,ND1+LBAS)*GDSC(NC2+KBAS,NB2+JBAS)
     &                                 *ABS(DENC(NC2+KBAS,NB2+JBAS))
C
            X( 5) = GDSC(NA1+IBAS,ND2+LBAS)*GDSC(NC1+KBAS,NB1+JBAS)
     &                                 *ABS(DENC(NC1+KBAS,NB1+JBAS))
            X( 6) = GDSC(NA1+IBAS,ND2+LBAS)*GDSC(NC1+KBAS,NB2+JBAS)
     &                                 *ABS(DENC(NC1+KBAS,NB2+JBAS))
            X( 7) = GDSC(NA1+IBAS,ND2+LBAS)*GDSC(NC2+KBAS,NB1+JBAS)
     &                                 *ABS(DENC(NC2+KBAS,NB1+JBAS))
            X( 8) = GDSC(NA1+IBAS,ND2+LBAS)*GDSC(NC2+KBAS,NB2+JBAS)
     &                                 *ABS(DENC(NC2+KBAS,NB2+JBAS))
C
            X( 9) = GDSC(NA2+IBAS,ND1+LBAS)*GDSC(NC1+KBAS,NB1+JBAS)
     &                                 *ABS(DENC(NC1+KBAS,NB1+JBAS))
            X(10) = GDSC(NA2+IBAS,ND1+LBAS)*GDSC(NC1+KBAS,NB2+JBAS)
     &                                 *ABS(DENC(NC1+KBAS,NB2+JBAS))
            X(11) = GDSC(NA2+IBAS,ND1+LBAS)*GDSC(NC2+KBAS,NB1+JBAS)
     &                                 *ABS(DENC(NC2+KBAS,NB1+JBAS))
            X(12) = GDSC(NA2+IBAS,ND1+LBAS)*GDSC(NC2+KBAS,NB2+JBAS)
     &                                 *ABS(DENC(NC2+KBAS,NB2+JBAS))
C
            X(13) = GDSC(NA2+IBAS,ND2+LBAS)*GDSC(NC1+KBAS,NB1+JBAS)
     &                                 *ABS(DENC(NC1+KBAS,NB1+JBAS))
            X(14) = GDSC(NA2+IBAS,ND2+LBAS)*GDSC(NC1+KBAS,NB2+JBAS)
     &                                 *ABS(DENC(NC1+KBAS,NB2+JBAS))
            X(15) = GDSC(NA2+IBAS,ND2+LBAS)*GDSC(NC2+KBAS,NB1+JBAS)
     &                                 *ABS(DENC(NC2+KBAS,NB1+JBAS))
            X(16) = GDSC(NA2+IBAS,ND2+LBAS)*GDSC(NC2+KBAS,NB2+JBAS)
     &                                 *ABS(DENC(NC2+KBAS,NB2+JBAS))
C
C           EARLY EXIT IF ANY OF THESE LARGER THAN SENS
            DO N=1,16
              IF(DABS(X(N)).GT.SENS) THEN
                ISCR(5) = 1
                ITAL    = ITAL+1
                GOTO 105
              ENDIF
            ENDDO
C
          ENDDO
        ENDDO
C
C       LOOP COMPLETED WITHOUT EARLY EXIT
        ISCR(5) = 0
C
C       EXIT IFLG(5) CHECK
105     CONTINUE
C
      ENDIF
C
C     6TH FLAG (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(6).EQ.1) THEN
C
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
C
C            X( 1) = RR3(M, 1)*ABS(DENC(ND1+LBAS,NB1+JBAS))
C            X( 2) = RR3(M, 2)*ABS(DENC(ND1+LBAS,NB2+JBAS))
C            X( 3) = RR3(M, 3)*ABS(DENC(ND2+LBAS,NB1+JBAS))
C            X( 4) = RR3(M, 4)*ABS(DENC(ND2+LBAS,NB2+JBAS))
CC
C            X( 5) = RR3(M, 5)*ABS(DENC(ND1+LBAS,NB1+JBAS))
C            X( 6) = RR3(M, 6)*ABS(DENC(ND1+LBAS,NB2+JBAS))
C            X( 7) = RR3(M, 7)*ABS(DENC(ND2+LBAS,NB1+JBAS))
C            X( 8) = RR3(M, 8)*ABS(DENC(ND2+LBAS,NB2+JBAS))
CC
C            X( 9) = RR3(M, 9)*ABS(DENC(ND1+LBAS,NB1+JBAS))
C            X(10) = RR3(M,10)*ABS(DENC(ND1+LBAS,NB2+JBAS))
C            X(11) = RR3(M,11)*ABS(DENC(ND2+LBAS,NB1+JBAS))
C            X(12) = RR3(M,12)*ABS(DENC(ND2+LBAS,NB2+JBAS))
CC
C            X(13) = RR3(M,13)*ABS(DENC(ND1+LBAS,NB1+JBAS))
C            X(14) = RR3(M,14)*ABS(DENC(ND1+LBAS,NB2+JBAS))
C            X(15) = RR3(M,15)*ABS(DENC(ND2+LBAS,NB1+JBAS))
C            X(16) = RR3(M,16)*ABS(DENC(ND2+LBAS,NB2+JBAS))

C
            X( 1) = GDSC(NA1+IBAS,NC1+KBAS)*GDSC(ND1+LBAS,NB1+JBAS)
     &                                 *ABS(DENC(ND1+LBAS,NB1+JBAS))
            X( 2) = GDSC(NA1+IBAS,NC1+KBAS)*GDSC(ND1+LBAS,NB2+JBAS)
     &                                 *ABS(DENC(ND1+LBAS,NB2+JBAS))
            X( 3) = GDSC(NA1+IBAS,NC1+KBAS)*GDSC(ND2+LBAS,NB1+JBAS)
     &                                 *ABS(DENC(ND2+LBAS,NB1+JBAS))
            X( 4) = GDSC(NA1+IBAS,NC1+KBAS)*GDSC(ND2+LBAS,NB2+JBAS)
     &                                 *ABS(DENC(ND2+LBAS,NB2+JBAS))
C
            X( 5) = GDSC(NA1+IBAS,NC2+KBAS)*GDSC(ND1+LBAS,NB1+JBAS)
     &                                 *ABS(DENC(ND1+LBAS,NB1+JBAS))
            X( 6) = GDSC(NA1+IBAS,NC2+KBAS)*GDSC(ND1+LBAS,NB2+JBAS)
     &                                 *ABS(DENC(ND1+LBAS,NB2+JBAS))
            X( 7) = GDSC(NA1+IBAS,NC2+KBAS)*GDSC(ND2+LBAS,NB1+JBAS)
     &                                 *ABS(DENC(ND2+LBAS,NB1+JBAS))
            X( 8) = GDSC(NA1+IBAS,NC2+KBAS)*GDSC(ND2+LBAS,NB2+JBAS)
     &                                 *ABS(DENC(ND2+LBAS,NB2+JBAS))
C
            X( 9) = GDSC(NA2+IBAS,NC1+KBAS)*GDSC(ND1+LBAS,NB1+JBAS)
     &                                 *ABS(DENC(ND1+LBAS,NB1+JBAS))
            X(10) = GDSC(NA2+IBAS,NC1+KBAS)*GDSC(ND1+LBAS,NB2+JBAS)
     &                                 *ABS(DENC(ND1+LBAS,NB2+JBAS))
            X(11) = GDSC(NA2+IBAS,NC1+KBAS)*GDSC(ND2+LBAS,NB1+JBAS)
     &                                 *ABS(DENC(ND2+LBAS,NB1+JBAS))
            X(12) = GDSC(NA2+IBAS,NC1+KBAS)*GDSC(ND2+LBAS,NB2+JBAS)
     &                                 *ABS(DENC(ND2+LBAS,NB2+JBAS))
C
            X(13) = GDSC(NA2+IBAS,NC2+KBAS)*GDSC(ND1+LBAS,NB1+JBAS)
     &                                 *ABS(DENC(ND1+LBAS,NB1+JBAS))
            X(14) = GDSC(NA2+IBAS,NC2+KBAS)*GDSC(ND1+LBAS,NB2+JBAS)
     &                                 *ABS(DENC(ND1+LBAS,NB2+JBAS))
            X(15) = GDSC(NA2+IBAS,NC2+KBAS)*GDSC(ND2+LBAS,NB1+JBAS)
     &                                 *ABS(DENC(ND2+LBAS,NB1+JBAS))
            X(16) = GDSC(NA2+IBAS,NC2+KBAS)*GDSC(ND2+LBAS,NB2+JBAS)
     &                                 *ABS(DENC(ND2+LBAS,NB2+JBAS))
C
C           EARLY EXIT IF ANY OF THESE LARGER THAN SENS
            DO N=1,16
              IF(DABS(X(N)).GT.SENS) THEN
                ISCR(6) = 1
                ITAL    = ITAL+1
                GOTO 106
              ENDIF
            ENDDO
C
          ENDDO
        ENDDO
C
C       LOOP COMPLETED WITHOUT EARLY EXIT
        ISCR(6) = 0
C
C       EXIT IFLG(6) CHECK
106     CONTINUE
C
      ENDIF
C
C     7TH FLAG (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(7).EQ.1) THEN
C
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
C
C            X( 1) = RR3(M, 1)*ABS(DENC(NC1+KBAS,NA1+IBAS))
C            X( 2) = RR3(M, 9)*ABS(DENC(NC1+KBAS,NA2+IBAS))
C            X( 3) = RR3(M, 5)*ABS(DENC(NC2+KBAS,NA1+IBAS))
C            X( 4) = RR3(M,13)*ABS(DENC(NC2+KBAS,NA2+IBAS))
CC
C            X( 5) = RR3(M, 3)*ABS(DENC(NC1+KBAS,NA1+IBAS))
C            X( 6) = RR3(M,11)*ABS(DENC(NC1+KBAS,NA2+IBAS))
C            X( 7) = RR3(M, 7)*ABS(DENC(NC2+KBAS,NA1+IBAS))
C            X( 8) = RR3(M,15)*ABS(DENC(NC2+KBAS,NA2+IBAS))
CC
C            X( 9) = RR3(M, 2)*ABS(DENC(NC1+KBAS,NA1+IBAS))
C            X(10) = RR3(M,10)*ABS(DENC(NC1+KBAS,NA2+IBAS))
C            X(11) = RR3(M, 6)*ABS(DENC(NC2+KBAS,NA1+IBAS))
C            X(12) = RR3(M,14)*ABS(DENC(NC2+KBAS,NA2+IBAS))
CC
C            X(13) = RR3(M, 4)*ABS(DENC(NC1+KBAS,NA1+IBAS))
C            X(14) = RR3(M,12)*ABS(DENC(NC1+KBAS,NA2+IBAS))
C            X(15) = RR3(M, 8)*ABS(DENC(NC2+KBAS,NA1+IBAS))
C            X(16) = RR3(M,16)*ABS(DENC(NC2+KBAS,NA2+IBAS))

C
            X( 1) = GDSC(NB1+JBAS,ND1+LBAS)*GDSC(NC1+KBAS,NA1+IBAS)
     &                                 *ABS(DENC(NC1+KBAS,NA1+IBAS))
            X( 2) = GDSC(NB1+JBAS,ND1+LBAS)*GDSC(NC1+KBAS,NA2+IBAS)
     &                                 *ABS(DENC(NC1+KBAS,NA2+IBAS))
            X( 3) = GDSC(NB1+JBAS,ND1+LBAS)*GDSC(NC2+KBAS,NA1+IBAS)
     &                                 *ABS(DENC(NC2+KBAS,NA1+IBAS))
            X( 4) = GDSC(NB1+JBAS,ND1+LBAS)*GDSC(NC2+KBAS,NA2+IBAS)
     &                                 *ABS(DENC(NC2+KBAS,NA2+IBAS))
C
            X( 5) = GDSC(NB1+JBAS,ND2+LBAS)*GDSC(NC1+KBAS,NA1+IBAS)
     &                                 *ABS(DENC(NC1+KBAS,NA1+IBAS))
            X( 6) = GDSC(NB1+JBAS,ND2+LBAS)*GDSC(NC1+KBAS,NA2+IBAS)
     &                                 *ABS(DENC(NC1+KBAS,NA2+IBAS))
            X( 7) = GDSC(NB1+JBAS,ND2+LBAS)*GDSC(NC2+KBAS,NA1+IBAS)
     &                                 *ABS(DENC(NC2+KBAS,NA1+IBAS))
            X( 8) = GDSC(NB1+JBAS,ND2+LBAS)*GDSC(NC2+KBAS,NA2+IBAS)
     &                                 *ABS(DENC(NC2+KBAS,NA2+IBAS))
C
            X( 9) = GDSC(NB2+JBAS,ND1+LBAS)*GDSC(NC1+KBAS,NA1+IBAS)
     &                                 *ABS(DENC(NC1+KBAS,NA1+IBAS))
            X(10) = GDSC(NB2+JBAS,ND1+LBAS)*GDSC(NC1+KBAS,NA2+IBAS)
     &                                 *ABS(DENC(NC1+KBAS,NA2+IBAS))
            X(11) = GDSC(NB2+JBAS,ND1+LBAS)*GDSC(NC2+KBAS,NA1+IBAS)
     &                                 *ABS(DENC(NC2+KBAS,NA1+IBAS))
            X(12) = GDSC(NB2+JBAS,ND1+LBAS)*GDSC(NC2+KBAS,NA2+IBAS)
     &                                 *ABS(DENC(NC2+KBAS,NA2+IBAS))
C
            X(13) = GDSC(NB2+JBAS,ND2+LBAS)*GDSC(NC1+KBAS,NA1+IBAS)
     &                                 *ABS(DENC(NC1+KBAS,NA1+IBAS))
            X(14) = GDSC(NB2+JBAS,ND2+LBAS)*GDSC(NC1+KBAS,NA2+IBAS)
     &                                 *ABS(DENC(NC1+KBAS,NA2+IBAS))
            X(15) = GDSC(NB2+JBAS,ND2+LBAS)*GDSC(NC2+KBAS,NA1+IBAS)
     &                                 *ABS(DENC(NC2+KBAS,NA1+IBAS))
            X(16) = GDSC(NB2+JBAS,ND2+LBAS)*GDSC(NC2+KBAS,NA2+IBAS)
     &                                 *ABS(DENC(NC2+KBAS,NA2+IBAS))
C
C           EARLY EXIT IF ANY OF THESE LARGER THAN SENS
            DO N=1,16
              IF(DABS(X(N)).GT.SENS) THEN
                ISCR(7) = 1
                ITAL    = ITAL+1
                GOTO 107
              ENDIF
            ENDDO
C
          ENDDO
        ENDDO
C
C       LOOP COMPLETED WITHOUT EARLY EXIT
        ISCR(7) = 0
C
C       EXIT IFLG(7) CHECK
107     CONTINUE
C
      ENDIF
C
C     8TH FLAG (EXCHNG):  ( MB, MA| MD, MC) = PAB*    (-MA,-MB| MC, MD)
C                                           = PAB*PCD*(-MA,-MB|-MC,-MD)
      IF(IFLG(8).EQ.1) THEN
C
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
C
C            X( 1) = RR2(M, 1)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X( 2) = RR2(M, 9)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X( 3) = RR2(M, 5)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X( 4) = RR2(M,13)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X( 5) = RR2(M, 3)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X( 6) = RR2(M,11)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X( 7) = RR2(M, 7)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X( 8) = RR2(M,15)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X( 9) = RR2(M, 2)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X(10) = RR2(M,10)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X(11) = RR2(M, 6)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X(12) = RR2(M,14)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X(13) = RR2(M, 4)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X(14) = RR2(M,12)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X(15) = RR2(M, 8)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X(16) = RR2(M,16)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
            X( 1) = GDSC(NB1+JBAS,NC1+KBAS)*GDSC(ND1+LBAS,NA1+IBAS)
     &                                 *ABS(DENC(ND1+LBAS,NA1+IBAS))
            X( 2) = GDSC(NB1+JBAS,NC1+KBAS)*GDSC(ND1+LBAS,NA2+IBAS)
     &                                 *ABS(DENC(ND1+LBAS,NA2+IBAS))
            X( 3) = GDSC(NB1+JBAS,NC1+KBAS)*GDSC(ND2+LBAS,NA1+IBAS)
     &                                 *ABS(DENC(ND2+LBAS,NA1+IBAS))
            X( 4) = GDSC(NB1+JBAS,NC1+KBAS)*GDSC(ND2+LBAS,NA2+IBAS)
     &                                 *ABS(DENC(ND2+LBAS,NA2+IBAS))
C
            X( 5) = GDSC(NB1+JBAS,NC2+KBAS)*GDSC(ND1+LBAS,NA1+IBAS)
     &                                 *ABS(DENC(ND1+LBAS,NA1+IBAS))
            X( 6) = GDSC(NB1+JBAS,NC2+KBAS)*GDSC(ND1+LBAS,NA2+IBAS)
     &                                 *ABS(DENC(ND1+LBAS,NA2+IBAS))
            X( 7) = GDSC(NB1+JBAS,NC2+KBAS)*GDSC(ND2+LBAS,NA1+IBAS)
     &                                 *ABS(DENC(ND2+LBAS,NA1+IBAS))
            X( 8) = GDSC(NB1+JBAS,NC2+KBAS)*GDSC(ND2+LBAS,NA2+IBAS)
     &                                 *ABS(DENC(ND2+LBAS,NA2+IBAS))
C
            X( 9) = GDSC(NB2+JBAS,NC1+KBAS)*GDSC(ND1+LBAS,NA1+IBAS)
     &                                 *ABS(DENC(ND1+LBAS,NA1+IBAS))
            X(10) = GDSC(NB2+JBAS,NC1+KBAS)*GDSC(ND1+LBAS,NA2+IBAS)
     &                                 *ABS(DENC(ND1+LBAS,NA2+IBAS))
            X(11) = GDSC(NB2+JBAS,NC1+KBAS)*GDSC(ND2+LBAS,NA1+IBAS)
     &                                 *ABS(DENC(ND2+LBAS,NA1+IBAS))
            X(12) = GDSC(NB2+JBAS,NC1+KBAS)*GDSC(ND2+LBAS,NA2+IBAS)
     &                                 *ABS(DENC(ND2+LBAS,NA2+IBAS))
C
            X(13) = GDSC(NB2+JBAS,NC2+KBAS)*GDSC(ND1+LBAS,NA1+IBAS)
     &                                 *ABS(DENC(ND1+LBAS,NA1+IBAS))
            X(14) = GDSC(NB2+JBAS,NC2+KBAS)*GDSC(ND1+LBAS,NA2+IBAS)
     &                                 *ABS(DENC(ND1+LBAS,NA2+IBAS))
            X(15) = GDSC(NB2+JBAS,NC2+KBAS)*GDSC(ND2+LBAS,NA1+IBAS)
     &                                 *ABS(DENC(ND2+LBAS,NA1+IBAS))
            X(16) = GDSC(NB2+JBAS,NC2+KBAS)*GDSC(ND2+LBAS,NA2+IBAS)
     &                                 *ABS(DENC(ND2+LBAS,NA2+IBAS))
C
C           EARLY EXIT IF ANY OF THESE LARGER THAN SENS
            DO N=1,16
              IF(DABS(X(N)).GT.SENS) THEN
                ISCR(8) = 1
                ITAL    = ITAL+1
                GOTO 108
              ENDIF
            ENDDO
C
          ENDDO
        ENDDO
C
C       LOOP COMPLETED WITHOUT EARLY EXIT
        ISCR(8) = 0
C
C       EXIT IFLG(8) CHECK
108     CONTINUE
C
      ENDIF
C
C     9TH FLAG (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(9).EQ.1) THEN
C
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
C
C            X( 1) = RR2(M, 1)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X( 2) = RR2(M, 5)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X( 3) = RR2(M, 9)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X( 4) = RR2(M,13)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X( 5) = RR2(M, 2)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X( 6) = RR2(M, 6)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X( 7) = RR2(M,10)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X( 8) = RR2(M,14)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X( 9) = RR2(M, 3)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X(10) = RR2(M, 7)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X(11) = RR2(M,11)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X(12) = RR2(M,15)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X(13) = RR2(M, 4)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X(14) = RR2(M, 8)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X(15) = RR2(M,12)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X(16) = RR2(M,16)*ABS(DENC(NC2+KBAS,ND2+LBAS))

C
            X( 1) = GDSC(NC1+KBAS,NB1+JBAS)*GDSC(NA1+IBAS,ND1+LBAS)
     &                                 *ABS(DENC(NA1+IBAS,ND1+LBAS))
            X( 2) = GDSC(NC1+KBAS,NB1+JBAS)*GDSC(NA1+IBAS,ND2+LBAS)
     &                                 *ABS(DENC(NA1+IBAS,ND2+LBAS))
            X( 3) = GDSC(NC1+KBAS,NB1+JBAS)*GDSC(NA2+IBAS,ND1+LBAS)
     &                                 *ABS(DENC(NA2+IBAS,ND1+LBAS))
            X( 4) = GDSC(NC1+KBAS,NB1+JBAS)*GDSC(NA2+IBAS,ND2+LBAS)
     &                                 *ABS(DENC(NA2+IBAS,ND2+LBAS))
C
            X( 5) = GDSC(NC1+KBAS,NB2+JBAS)*GDSC(NA1+IBAS,ND1+LBAS)
     &                                 *ABS(DENC(NA1+IBAS,ND1+LBAS))
            X( 6) = GDSC(NC1+KBAS,NB2+JBAS)*GDSC(NA1+IBAS,ND2+LBAS)
     &                                 *ABS(DENC(NA1+IBAS,ND2+LBAS))
            X( 7) = GDSC(NC1+KBAS,NB2+JBAS)*GDSC(NA2+IBAS,ND1+LBAS)
     &                                 *ABS(DENC(NA2+IBAS,ND1+LBAS))
            X( 8) = GDSC(NC1+KBAS,NB2+JBAS)*GDSC(NA2+IBAS,ND2+LBAS)
     &                                 *ABS(DENC(NA2+IBAS,ND2+LBAS))
C
            X( 9) = GDSC(NC2+KBAS,NB1+JBAS)*GDSC(NA1+IBAS,ND1+LBAS)
     &                                 *ABS(DENC(NA1+IBAS,ND1+LBAS))
            X(10) = GDSC(NC2+KBAS,NB1+JBAS)*GDSC(NA1+IBAS,ND2+LBAS)
     &                                 *ABS(DENC(NA1+IBAS,ND2+LBAS))
            X(11) = GDSC(NC2+KBAS,NB1+JBAS)*GDSC(NA2+IBAS,ND1+LBAS)
     &                                 *ABS(DENC(NA2+IBAS,ND1+LBAS))
            X(12) = GDSC(NC2+KBAS,NB1+JBAS)*GDSC(NA2+IBAS,ND2+LBAS)
     &                                 *ABS(DENC(NA2+IBAS,ND2+LBAS))
C
            X(13) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(NC1+KBAS,ND1+LBAS)
     &                                 *ABS(DENC(NC1+KBAS,ND1+LBAS))
            X(14) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(NC1+KBAS,ND2+LBAS)
     &                                 *ABS(DENC(NC1+KBAS,ND2+LBAS))
            X(15) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(NC2+KBAS,ND1+LBAS)
     &                                 *ABS(DENC(NC2+KBAS,ND1+LBAS))
            X(16) = GDSC(NA2+IBAS,NB2+JBAS)*GDSC(NC2+KBAS,ND2+LBAS)
     &                                 *ABS(DENC(NC2+KBAS,ND2+LBAS))
C
C           EARLY EXIT IF ANY OF THESE LARGER THAN SENS
            DO N=1,16
              IF(DABS(X(N)).GT.SENS) THEN
                ISCR(9) = 1
                ITAL    = ITAL+1
                GOTO 109
              ENDIF
            ENDDO
C
          ENDDO
        ENDDO
C
C       LOOP COMPLETED WITHOUT EARLY EXIT
        ISCR(9) = 0
C
C       EXIT IFLG(9) CHECK
109     CONTINUE
C
      ENDIF
C
C     10TH FLAG (EXCHNG): ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(10).EQ.1) THEN
C
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
C
C            X( 1) = RR3(M, 1)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X( 2) = RR3(M, 3)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X( 3) = RR3(M, 2)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X( 4) = RR3(M, 4)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X( 5) = RR3(M, 9)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X( 6) = RR3(M,11)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X( 7) = RR3(M,10)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X( 8) = RR3(M,12)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X( 9) = RR3(M, 5)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X(10) = RR3(M, 7)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X(11) = RR3(M, 6)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X(12) = RR3(M, 8)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X(13) = RR3(M,13)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X(14) = RR3(M,15)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X(15) = RR3(M,14)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X(16) = RR3(M,16)*ABS(DENC(NC2+KBAS,ND2+LBAS))

C
            X( 1) = GDSC(NC1+KBAS,NA1+IBAS)*GDSC(NB1+JBAS,ND1+LBAS)
     &                                 *ABS(DENC(NB1+JBAS,ND1+LBAS))
            X( 2) = GDSC(NC1+KBAS,NA1+IBAS)*GDSC(NB1+JBAS,ND2+LBAS)
     &                                 *ABS(DENC(NB1+JBAS,ND2+LBAS))
            X( 3) = GDSC(NC1+KBAS,NA1+IBAS)*GDSC(NB2+JBAS,ND1+LBAS)
     &                                 *ABS(DENC(NB2+JBAS,ND1+LBAS))
            X( 4) = GDSC(NC1+KBAS,NA1+IBAS)*GDSC(NB2+JBAS,ND2+LBAS)
     &                                 *ABS(DENC(NB2+JBAS,ND2+LBAS))
C
            X( 5) = GDSC(NC1+KBAS,NA2+IBAS)*GDSC(NB1+JBAS,ND1+LBAS)
     &                                 *ABS(DENC(NB1+JBAS,ND1+LBAS))
            X( 6) = GDSC(NC1+KBAS,NA2+IBAS)*GDSC(NB1+JBAS,ND2+LBAS)
     &                                 *ABS(DENC(NB1+JBAS,ND2+LBAS))
            X( 7) = GDSC(NC1+KBAS,NA2+IBAS)*GDSC(NB2+JBAS,ND1+LBAS)
     &                                 *ABS(DENC(NB2+JBAS,ND1+LBAS))
            X( 8) = GDSC(NC1+KBAS,NA2+IBAS)*GDSC(NB2+JBAS,ND2+LBAS)
     &                                 *ABS(DENC(NB2+JBAS,ND2+LBAS))
C
            X( 9) = GDSC(NC2+KBAS,NA1+IBAS)*GDSC(NB1+JBAS,ND1+LBAS)
     &                                 *ABS(DENC(NB1+JBAS,ND1+LBAS))
            X(10) = GDSC(NC2+KBAS,NA1+IBAS)*GDSC(NB1+JBAS,ND2+LBAS)
     &                                 *ABS(DENC(NB1+JBAS,ND2+LBAS))
            X(11) = GDSC(NC2+KBAS,NA1+IBAS)*GDSC(NB2+JBAS,ND1+LBAS)
     &                                 *ABS(DENC(NB2+JBAS,ND1+LBAS))
            X(12) = GDSC(NC2+KBAS,NA1+IBAS)*GDSC(NB2+JBAS,ND2+LBAS)
     &                                 *ABS(DENC(NB2+JBAS,ND2+LBAS))
C
            X(13) = GDSC(NC2+KBAS,NA2+IBAS)*GDSC(NB1+JBAS,ND1+LBAS)
     &                                 *ABS(DENC(NB1+JBAS,ND1+LBAS))
            X(14) = GDSC(NC2+KBAS,NA2+IBAS)*GDSC(NB1+JBAS,ND2+LBAS)
     &                                 *ABS(DENC(NB1+JBAS,ND2+LBAS))
            X(15) = GDSC(NC2+KBAS,NA2+IBAS)*GDSC(NB2+JBAS,ND1+LBAS)
     &                                 *ABS(DENC(NB2+JBAS,ND1+LBAS))
            X(16) = GDSC(NC2+KBAS,NA2+IBAS)*GDSC(NB2+JBAS,ND2+LBAS)
     &                                 *ABS(DENC(NB2+JBAS,ND2+LBAS))
C
C           EARLY EXIT IF ANY OF THESE LARGER THAN SENS
            DO N=1,16
              IF(DABS(X(N)).GT.SENS) THEN
                ISCR(10) = 1
                ITAL    = ITAL+1
                GOTO 110
              ENDIF
            ENDDO
C
          ENDDO
        ENDDO
C
C       LOOP COMPLETED WITHOUT EARLY EXIT
        ISCR(10) = 0
C
C       EXIT IFLG(10) CHECK
110     CONTINUE
C
      ENDIF
C
C     11TH FLAG (EXCHNG): ( MD, MC| MA, MB) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(11).EQ.1) THEN
C
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
C
C            X( 1) = RR3(M, 1)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X( 2) = RR3(M, 5)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X( 3) = RR3(M, 9)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X( 4) = RR3(M,13)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X( 5) = RR3(M, 2)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X( 6) = RR3(M, 6)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X( 7) = RR3(M,10)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X( 8) = RR3(M,14)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X( 9) = RR3(M, 3)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X(10) = RR3(M, 7)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X(11) = RR3(M,11)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X(12) = RR3(M,15)*ABS(DENC(NC2+KBAS,ND2+LBAS))
CC
C            X(13) = RR3(M, 4)*ABS(DENC(NC1+KBAS,ND1+LBAS))
C            X(14) = RR3(M, 8)*ABS(DENC(NC1+KBAS,ND2+LBAS))
C            X(15) = RR3(M,12)*ABS(DENC(NC2+KBAS,ND1+LBAS))
C            X(16) = RR3(M,16)*ABS(DENC(NC2+KBAS,ND2+LBAS))

C
            X( 1) = GDSC(ND1+LBAS,NB1+JBAS)*GDSC(NA1+IBAS,NC1+KBAS)
     &                                 *ABS(DENC(NA1+IBAS,NC1+KBAS))
            X( 2) = GDSC(ND1+LBAS,NB1+JBAS)*GDSC(NA1+IBAS,NC2+KBAS)
     &                                 *ABS(DENC(NA1+IBAS,NC2+KBAS))
            X( 3) = GDSC(ND1+LBAS,NB1+JBAS)*GDSC(NA2+IBAS,NC1+KBAS)
     &                                 *ABS(DENC(NA2+IBAS,NC1+KBAS))
            X( 4) = GDSC(ND1+LBAS,NB1+JBAS)*GDSC(NA2+IBAS,NC2+KBAS)
     &                                 *ABS(DENC(NA2+IBAS,NC2+KBAS))
C
            X( 5) = GDSC(ND1+LBAS,NB2+JBAS)*GDSC(NA1+IBAS,NC1+KBAS)
     &                                 *ABS(DENC(NA1+IBAS,NC1+KBAS))
            X( 6) = GDSC(ND1+LBAS,NB2+JBAS)*GDSC(NA1+IBAS,NC2+KBAS)
     &                                 *ABS(DENC(NA1+IBAS,NC2+KBAS))
            X( 7) = GDSC(ND1+LBAS,NB2+JBAS)*GDSC(NA2+IBAS,NC1+KBAS)
     &                                 *ABS(DENC(NA2+IBAS,NC1+KBAS))
            X( 8) = GDSC(ND1+LBAS,NB2+JBAS)*GDSC(NA2+IBAS,NC2+KBAS)
     &                                 *ABS(DENC(NA2+IBAS,NC2+KBAS))
C
            X( 9) = GDSC(ND2+LBAS,NB1+JBAS)*GDSC(NA1+IBAS,NC1+KBAS)
     &                                 *ABS(DENC(NA1+IBAS,NC1+KBAS))
            X(10) = GDSC(ND2+LBAS,NB1+JBAS)*GDSC(NA1+IBAS,NC2+KBAS)
     &                                 *ABS(DENC(NA1+IBAS,NC2+KBAS))
            X(11) = GDSC(ND2+LBAS,NB1+JBAS)*GDSC(NA2+IBAS,NC1+KBAS)
     &                                 *ABS(DENC(NA2+IBAS,NC1+KBAS))
            X(12) = GDSC(ND2+LBAS,NB1+JBAS)*GDSC(NA2+IBAS,NC2+KBAS)
     &                                 *ABS(DENC(NA2+IBAS,NC2+KBAS))
C
            X(13) = GDSC(ND2+LBAS,NB2+JBAS)*GDSC(NA1+IBAS,NC1+KBAS)
     &                                 *ABS(DENC(NA1+IBAS,NC1+KBAS))
            X(14) = GDSC(ND2+LBAS,NB2+JBAS)*GDSC(NA1+IBAS,NC2+KBAS)
     &                                 *ABS(DENC(NA1+IBAS,NC2+KBAS))
            X(15) = GDSC(ND2+LBAS,NB2+JBAS)*GDSC(NA2+IBAS,NC1+KBAS)
     &                                 *ABS(DENC(NA2+IBAS,NC1+KBAS))
            X(16) = GDSC(ND2+LBAS,NB2+JBAS)*GDSC(NA2+IBAS,NC2+KBAS)
     &                                 *ABS(DENC(NA2+IBAS,NC2+KBAS))
C
C           EARLY EXIT IF ANY OF THESE LARGER THAN SENS
            DO N=1,16
              IF(DABS(X(N)).GT.SENS) THEN
                ISCR(11) = 1
                ITAL    = ITAL+1
                GOTO 111
              ENDIF
            ENDDO
C
          ENDDO
        ENDDO
C
C       LOOP COMPLETED WITHOUT EARLY EXIT
        ISCR(11) = 0
C
C       EXIT IFLG(11) CHECK
111     CONTINUE
C
      ENDIF
C
200   CONTINUE
C
C     IF ALL FLAGS HAVE BEEN SWITCHED OFF, SKIP THE BATCH
      IF(ITAL.EQ.0) THEN
        ISKP = 1
      ELSE
        ISKP = 0
      ENDIF
C
      RETURN
      END

