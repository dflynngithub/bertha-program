      SUBROUTINE SPECSUM(NLIST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   SSSSSS  PPPPPPP  EEEEEEEE  CCCCCC   SSSSSS  UU    UU MM       MM   C
C  SS    SS PP    PP EE       CC    CC SS    SS UU    UU MMM     MMM   C
C  SS       PP    PP EE       CC       SS       UU    UU MMMM   MMMM   C
C   SSSSSS  PP    PP EEEEEE   CC        SSSSSS  UU    UU MM MM MM MM   C
C        SS PPPPPPP  EE       CC             SS UU    UU MM  MMM  MM   C
C  SS    SS PP       EE       CC    CC SS    SS UU    UU MM   M   MM   C
C   SSSSSS  PP       EEEEEEEE  CCCCCC   SSSSSS   UUUUUU  MM       MM   C
C                                                                      C
C -------------------------------------------------------------------- C
C  SPECSUM SUMMARISES THE ENERGY LEVELS OF AN SCF ITERATION, USING     C
C  ATOMIC TERM SYMBOLS, ATOMIC CENTERS AND FRACTIONAL POPULATIONS.     C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    NLIST - SPECIFIES THE LENGTH OF LISTING ALGORITHM LIMIT.          C
C -------------------------------------------------------------------- C
C  DFNOTE: LOOK AT 'operator.f' ROUTINE 'QASSIGN' FOR MORE OPTIONS.    C
C**********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26,MN=20)
C
      COMPLEX*16 C(MDM,MDM)
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           GDIR(MDM,MDM),GXCH(MDM,MDM),QDIR(MDM,MDM),
     &           QXCH(MDM,MDM),BDIR(MDM,MDM),BXCH(MDM,MDM),FOCK(MDM,MDM)
C
      CHARACTER*1 ELLTERM
      CHARACTER*2 ELMNT(120)
C
      COMMON/ATOM/ELMNT
      COMMON/COEF/C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/MTRX/OVAP,HNUC,HKIN,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/SHLL/ALPHA,BETA,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN
      COMMON/QBRK/IVLC(MDM),IVLM(MDM),IVLK(MDM),IVLL(MDM),
     &            IVLJ(MDM),IVLN(MDM),PRTY(MDM),
     &            NKM(MCT,MN,MKP,2*MMV),ISORT(MDM),
     &            NMAX(MCT),KMAX(MCT,MN)
      COMMON/QNMS/ICNLAB(MDM),KQNLAB(MDM),MQNLAB(MDM)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR,IOCCM0
C
C     INITIALISE ATOMIC TERM SYMBOL ARRAYS
      DO IOCC=1,NDIM-NSHIFT
        IVLC(IOCC) = 0
        IVLM(IOCC) = 0
        IVLK(IOCC) = 0
        IVLL(IOCC) = 0
        IVLJ(IOCC) = 0
        IVLN(IOCC) = 0
        PRTY(IOCC) = 0.0D0
      ENDDO
C
C     AND ALSO THE NSHELL, KAPPA, MJ ARRAY
      DO ICNT=1,MCT
        NMAX(ICNT) = 0
        DO INPR=1,MN
          KMAX(ICNT,INPR) = 0
          DO IKAP=1,MKP
            DO IMV=1,2*MMV
              NKM(ICNT,INPR,IKAP,IMV) = 0
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     PRINT TITLE TO TERMINAL/FILE
10    FORMAT(1X,'Orb.',1X,'Center',2X,'Term sym.',3X,'m_j',10X,
     &                                    'Energy (au)',3X,'Population')
11    FORMAT(1X,I3,2X,I2,'(',A,')',2X,I2,A,'_',I1,'/2',4X,I2,'/2',3X,
     &                                                  F18.12,3X,F10.8)
C
      WRITE(6,10)
      WRITE(7,10)
      WRITE(6, *) REPEAT('-',62)
      WRITE(7, *) REPEAT('-',62)
C
C     FOR ALL POSTITIVE-ENERGY ORBITALS (OCCUPIED, VIRTUAL OR CONTIVUM)
      DO IOCC=1,NLIST
C
C       IGNORE NEGATIVE SPECTRUM
        MOCC = IOCC + NSHIFT
C
C       COUNTERS NEEDED THROUGHOUT THESE LOOPS
        BIN  = 0.0D0
        TOT  = 0.0D0
C
C       LOOP OVER NUCLEAR CENTERS
        DO ICNT=1,NCNT
C
C         LOOP OVER |MQN|
          DO IMV=1,MMV
            MQN = 2*IMV-1
C
C           MQN<0 BLOCKS

C           ALL KQN VALUES AVAILABLE FOR THE NUCLEAR CENTER
            DO IKAP=1,NKAP(ICNT)
              KQN = KVALS(IKAP,ICNT)
              IF(KQN.GT.0) THEN
                LQN = KQN
              ELSE
                LQN =-KQN-1
              ENDIF
C             NUMBER OF BASIS FUNCTIONS OF THIS TYPE AND COEFF ADDRESS
              NFUN  = NFUNCT(LQN+1,ICNT)
              MQMAX = 2*IABS(KQN)-1
              IADD  = LARGE(ICNT,IKAP,MQN)
              IF(MQN.LE.MQMAX) THEN
C               SUM UP MAGNITUDES OF EXPANSION COEFFICIENTS OF THIS TYPE
                DO IFN=1,NFUN
                  BIN = BIN + ABS(C(IADD+IFN       ,MOCC))
                  BIN = BIN + ABS(C(IADD+IFN+NSHIFT,MOCC))
                  TOT = TOT + ABS(C(IADD+IFN       ,MOCC))
                  TOT = TOT + ABS(C(IADD+IFN+NSHIFT,MOCC))
                ENDDO
C               PRINT OUT THE SUM OF THIS BLOCK
C               IF THIS SUM IS THE BIGGEST SO FAR, ATTRIBUTE Q.N. VALUES
                IF(BIN.GT.PRTY(IOCC)) THEN
                  PRTY(IOCC) = BIN
                  IVLC(IOCC) = ICNT
                  IVLM(IOCC) =-MQN
                  IVLK(IOCC) = KQN
                  IVLL(IOCC) = LQN
                  IVLJ(IOCC) = 2*IABS(KQN)-1
C                 THE NQN MUST BE GREATER THAN LQN, AND OBEYS PAULI
                  NSHELL     =  LQN + 1
                  DO JORB=1,IOCC-1
                    IF(IVLC(IOCC).NE.IVLC(JORB)) GOTO 1
                    IF(IVLK(IOCC).NE.IVLK(JORB)) GOTO 1
                    IF(IVLM(IOCC).NE.IVLM(JORB)) GOTO 1
                    NSHELL = NSHELL + 1
1                   CONTINUE
                  ENDDO
                  IF(NSHELL.LE.LQN) NSHELL = NSHELL + LQN
                  IVLN(IOCC) = NSHELL
C                 DETERMINE UPPER KQN FOR THIS NQN AND NUCLEUS
                  IF(IKAP.GT.KMAX(ICNT,NSHELL)) THEN
                    KMAX(ICNT,NSHELL) = IKAP
                  ENDIF
C                 GIVE ADDRESS IN THE NKM MATRIX FOR EASE OF USE LATER
                  NKM(ICNT,NSHELL,KA,MQN) = IOCC
                ENDIF
                BIN = 0.0D0
              ENDIF
            ENDDO
C
C           MQN>0 BLOCKS
C
C           ALL KQN VALUES AVAILABLE FOR THE NUCLEAR CENTER
            DO IKAP=1,NKAP(ICNT)
              KQN = KVALS(IKAP,ICNT)
              IF(KQN.GT.0) THEN
                LQN = KQN
              ELSE
                LQN =-KQN-1
              ENDIF
C             NUMBER OF BASIS FUNCTIONS OF THIS TYPE AND COEFF ADDRESS
              NFUN  = NFUNCT(LQN+1,ICNT)
              MQMAX = 2*IABS(KQN)-1
              IADD  = LARGE(ICNT,IKAP,MQN+1)
              IF(MQN.LE.MQMAX) THEN
C               SUM UP MAGNITUDES OF EXPANSION COEFFICIENTS OF THIS TYPE
                DO IFN=1,NFUN
                  BIN = BIN + ABS(C(IADD+IFN       ,MOCC))
                  BIN = BIN + ABS(C(IADD+IFN+NSHIFT,MOCC))
                  TOT = TOT + ABS(C(IADD+IFN       ,MOCC))
                  TOT = TOT + ABS(C(IADD+IFN+NSHIFT,MOCC))
                ENDDO
C               PRINT OUT THE SUM OF THIS BLOCK
C               IF THIS SUM IS THE BIGGEST SO FAR, ATTRIBUTE Q.N. VALUES
                IF(BIN.GT.PRTY(IOCC)) THEN
                  PRTY(IOCC) =  BIN
                  IVLC(IOCC) =  ICNT
                  IVLM(IOCC) =  MQN
                  IVLK(IOCC) =  KQN
                  IVLL(IOCC) =  LQN
                  IVLJ(IOCC) =  2*IABS(KQN)-1
C                 THE NQN MUST BE GREATER THAN LQN, AND OBEYS PAULI
                  NSHELL     =  LQN + 1
                  DO JORB=1,IOCC-1
                    IF(IVLC(IOCC).NE.IVLC(JORB)) GOTO 2
                    IF(IVLK(IOCC).NE.IVLK(JORB)) GOTO 2
                    IF(IVLM(IOCC).NE.IVLM(JORB)) GOTO 2
                    NSHELL = NSHELL + 1
2                   CONTINUE
                  ENDDO
                  IVLN(IOCC) = NSHELL
C                 DETERMINE UPPER KQN FOR THIS NQN AND NUCLEUS
                  IF(IKAP.GT.KMAX(ICNT,NSHELL)) THEN
                    KMAX(ICNT,NSHELL) = IKAP
                  ENDIF
C                 GIVE ADDRESS IN THE NKM MATRIX FOR EASE OF USE LATER
                  NKM(ICNT,NSHELL,IKAP,MQN+1) = IOCC
                ENDIF
                BIN = 0.0D0
              ENDIF
            ENDDO
          ENDDO 
C       DETERMINE UPPER NQN FOR THIS NUCLEUS
        IF(NSHELL.GT.NMAX(ICNT)) NMAX(ICNT) = NSHELL
        NSHELL = 0
        ENDDO
C
C       DIVIDE BLOCK SUM COUNTER BY THE TOTAL SUM FOR THE WHOLE EVECTOR
        PRTY(IOCC) = PRTY(IOCC)/TOT
C
C       SUMMARY OF RESULTS
        WRITE(6,11) IOCC,IVLC(IOCC),ELMNT(IZNUC(IVLC(IOCC))),
     &              IVLN(IOCC),ELLTERM(IVLL(IOCC)),IVLJ(IOCC),
     &              IVLM(IOCC),EIGEN(IOCC+NSHIFT),PRTY(IOCC)
        WRITE(7,11) IOCC,IVLC(IOCC),ELMNT(IZNUC(IVLC(IOCC))),
     &              IVLN(IOCC),ELLTERM(IVLL(IOCC)),IVLJ(IOCC),
     &              IVLM(IOCC),EIGEN(IOCC+NSHIFT),PRTY(IOCC)
        IF(IOCC.EQ.NOCC) THEN
          WRITE(6, *) REPEAT('-',62)
          WRITE(7, *) REPEAT('-',62)
        ENDIF
      ENDDO
C
      WRITE(6, *) REPEAT('=',62)
      WRITE(7, *) REPEAT('=',62)
C
      RETURN
      END
