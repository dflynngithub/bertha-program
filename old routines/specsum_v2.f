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
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26)
C
      CHARACTER*1 ELLTERM
      CHARACTER*2 ELMNT(120)
C
      COMPLEX*16 C(MDM,MDM)
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           GDIR(MDM,MDM),GXCH(MDM,MDM),QDIR(MDM,MDM),
     &           QXCH(MDM,MDM),BDIR(MDM,MDM),BXCH(MDM,MDM),FOCK(MDM,MDM)
C
      DIMENSION FRAC(MDM,MCT,MKP,2*MMV),NPR(MCT,MKP,2*MMV)
C
      COMMON/ATOM/ELMNT
      COMMON/COEF/C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/MTRX/OVAP,HNUC,HKIN,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/QNMS/ICNLAB(MDM),KQNLAB(MDM),MQNLAB(MDM)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR,IOCCM0
C
C     INITIALISE COUNTER ARRAYS
      DO ICNT=1,NCNT
        DO IKAP=1,NKAP(ICNT)
          IKQN = KVALS(IKAP,ICNT)
          NMV  = 2*IABS(IKQN)
          DO IMV=1,NMV
            DO IOCC=1,NLIST
              FRAC(IOCC,ICNT,IKAP,IMV) = 0.0D0
            ENDDO
            NPR(ICNT,IKAP,IMV)
          ENDDO
        ENDDO
      ENDDO
C
C     LOOP OVER OCCUPIED ORBITALS
      DO IOCC=1,NLIST
C
C       IGNORE NEGATIVE SPECTRUM
        MOCC = IOCC + NSHIFT
C
C       LOOP OVER FOCK MATRIX ADDRESS BLOCKS
        DO I=1,NDIM-NSHIFT
          IS = I+NSHIFT
          ICNT = ICNLAB(I)
          IKQN = KQNLAB(I)
          IMQN = MQNLAB(I)
          IF(IKQN.LT.0) THEN
            IKAP =-2*IKQN-1
          ELSE
            IKAP = 2*IKQN
          ENDIF
          IMV = (IABS(IMQN)+1)/2
          DO J=1,NDIM-NSHIFT
            JS = J+NSHIFT
            JCNT = ICNLAB(J)
            JKQN = KQNLAB(J)
            JMQN = MQNLAB(J)
C
C           LARGE AND SMALL CONTRIBUTIONS
            TMP = DREAL(DCONJG(C(I ,MOCC))*C(J ,MOCC)*OVAP(I ,J ))
     &          + DREAL(DCONJG(C(IS,MOCC))*C(JS,MOCC)*OVAP(IS,JS))

C           DECIDE WHERE TO PUT CONTRIBUTION
            IF(ICNT.EQ.JCNT.AND.IKQN.EQ.JKQN.AND.IMQN.EQ.JMQN) THEN
              IF(IMQN.LT.0) THEN
                FRAC(IOCC,ICNT,IKAP,IMV  ) = FRAC(IOCC,ICNT,IKAP,IMV  ) 
     &                                     + TMP
              ELSE
                FRAC(IOCC,ICNT,IKAP,IMV+1) = FRAC(IOCC,ICNT,IKAP,IMV+1) 
     &                                     + TMP
              ENDIF
            ENDIF
C          
          ENDDO
C
        ENDDO
C
C     END LOOP OVER OCCUPIED ORBITALS
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
C     SUMMARISE RESULTS
      DO IOCC=1,NLIST
C
C       SEARCH FOR HIGHEST POPULATED ATOMIC STATE
        PLTN = 0.0D0
        DO ICNT=1,NCNT
          DO IKAP=1,NKAP(ICNT)
            IKQN = KVALS(IKAP,ICNT)
            NMV  = 2*IABS(IKQN)
            DO IMV=1,NMV
              IF(FRAC(IOCC,ICNT,IKAP,IMV).GT.PLTN) THEN
                PLTN = FRAC(IOCC,ICNT,IKAP,IMV)
                NPR(ICNT,IKAP,IMV) = NPR(ICNT,IKAP,IMV)+1
                KNQN = NPR(ICNT,IKAP,IMV)
                KCNT = ICNT
                KKQN = IKQN
                IF(KKQN.LT.0) THEN
                  KLQN =-KKQN-1
                ELSE
                  KLQN = KKQN
                ENDIF
                IF(MOD(IMV,2).EQ.1) THEN
                  KMQN = IMV
                ELSE
                  KMQN = IMV-1
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
C       OUTPUT TO TERMINAL
        WRITE(6,11) KKQN,KCNT,ELMNT(IZNUC(KCNT)),KNQN,ELLTERM(KLQN),
     &              IABS(KKQN),KMQN,EIGEN(IOCC+NSHIFT),PLTN
        WRITE(7,11) KKQN,KCNT,ELMNT(IZNUC(KCNT)),KNQN,ELLTERM(KLQN),
     &              IABS(KKQN),KMQN,EIGEN(IOCC+NSHIFT),PLTN
        IF(IOCC.EQ.NOCC) THEN
          WRITE(6, *) REPEAT('-',62)
          WRITE(7, *) REPEAT('-',62)
        ENDIF
C
      ENDDO
C
      WRITE(6, *) REPEAT('=',62)
      WRITE(7, *) REPEAT('=',62)
C
      RETURN
      END
