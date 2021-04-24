      SUBROUTINE SPECTRM0(IBND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    SSSSSS  PPPPPPP   CCCCCC TTTTTTTT RRRRRRR  MM       MM  000000    C
C   SS    SS PP    PP CC    CC   TT    RR    RR MMM     MMM 00   000   C
C   SS       PP    PP CC         TT    RR    RR MMMM   MMMM 00  0000   C
C    SSSSSS  PP    PP CC         TT    RR    RR MM MM MM MM 00 00 00   C
C         SS PPPPPPP  CC         TT    RRRRRRR  MM  MMM  MM 0000  00   C
C   SS    SS PP       CC    CC   TT    RR    RR MM   M   MM 000   00   C
C    SSSSSS  PP        CCCCCC    TT    RR    RR MM       MM  000000    C
C                                                                      C
C -------------------------------------------------------------------- C
C  SPECTRM SUMMARISES THE ENERGY LEVELS OF AN SCF ITERATION, USING     C
C  ATOMIC TERM SYMBOLS, ATOMIC CENTRES AND FRACTIONAL POPULATIONS.     C
C  IT ALSO APPLIES SYMMETRY LABELS TO EACH STATE AND SORTS MQN M'FOLDS.C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    IBND - NUMBER OF BOUND ORBITALS TO INCLUDE.                       C
C -------------------------------------------------------------------- C
C  DFNOTE: FIGURE OUT HOW TO SORT THROUGH MANIFOLD IN ATOMIC CASE.     C
C**********************************************************************C
      PARAMETER(MDM=1500,MBS=40,MCT=4,MKP=11)
C
      CHARACTER*1 CHL,CHS,ELLTERM
      CHARACTER*2 ELMT(120),ELA
      CHARACTER*4 HMLTN
C
      DIMENSION FRAC(MDM,MCT,MKP,MKP+1),NPR(MCT,MKP,MKP+1)
      DIMENSION ICNLST(MDM),KQNLST(MDM),MQNLST(MDM),NQNLST(MDM),
     &          POPLST(MDM)
C
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           GDIR(MDM,MDM),GXCH(MDM,MDM),QDIR(MDM,MDM),
     &           QXCH(MDM,MDM),BDIR(MDM,MDM),BXCH(MDM,MDM),
     &           VUEH(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 CTEMP,ROT1,ROT2
C
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGEN(MDM)
      COMMON/MDLV/ELMT
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/OCPD/RAVG(MDM,4),IOCPN(MDM),IOCCM0
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/QNMS/LABICN(MDM),LABKQN(MDM),LABMQN(MDM)
      COMMON/SPEC/BSET(MBS,MKP,MCT),BXYZ(3,MCT),ZNUC(MCT),AMSS(MCT),
     &            CNUC(MCT),PNUC,LRGE(MCT,MKP,MKP+1),NFNC(MKP,MCT),
     &            KAPA(MKP,MCT),IZNC(MCT),IQNC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSKP,NOCC,NVRT
C
      DATA PI/3.1415926535897932D0/
C
C     SENSITIVITY TOLERANCE PARAMETER
      EPS = 1.0D-14
C
C     INITIALISE COUNTER ARRAYS
      DO ICNT=1,NCNT
        DO IKAP=1,NKAP(ICNT)
          IKQN = KAPA(IKAP,ICNT)
          IF(IKQN.LT.0) THEN
            ILQN =-IKQN-1
          ELSE
            ILQN = IKQN
          ENDIF
          NMV = 2*IABS(IKQN)
          DO IMV=1,NMV
            DO IOCC=1,IBND
              FRAC(IOCC,ICNT,IKAP,IMV) = 0.0D0
            ENDDO
            NPR(ICNT,IKAP,IMV) = ILQN
          ENDDO
        ENDDO
      ENDDO
C
C     CALCULATE FRACTIONAL OCCUPATION FOR EACH SOLUTION VECTOR
C
C     LOOP OVER OCCUPIED ORBITALS
      DO IOCC=1,IBND
C
C       IGNORE NEGATIVE SPECTRUM
        MOCC = IOCC+NSKP
C
C       LOOP OVER FOCK MATRIX ADDRESS BLOCKS
        DO I=1,NDIM-NSKP
          K    = I+NSKP
          ICNT = LABICN(I)
          IKQN = LABKQN(I)
          IMQN = LABMQN(I)
          IF(IKQN.LT.0) THEN
            IKAP =-2*IKQN-1
          ELSE
            IKAP = 2*IKQN
          ENDIF
          IF(IMQN.LT.0) THEN
            IMV = IABS(IMQN)
          ELSE
            IMV = IMQN+1
          ENDIF
C
          DO J=1,NDIM-NSKP
            L    = J+NSKP
            JCNT = LABICN(J)
            JKQN = LABKQN(J)
            JMQN = LABMQN(J)
C
C           LARGE AND SMALL CONTRIBUTIONS
            TMP = DREAL(DCONJG(COEF(I,MOCC))*COEF(J,MOCC)*OVAP(I,J))
     &          + DREAL(DCONJG(COEF(K,MOCC))*COEF(L,MOCC)*OVAP(K,L))
C
C           DECIDE WHERE TO PUT CONTRIBUTION
            IF(ICNT.EQ.JCNT.AND.IKQN.EQ.JKQN.AND.IMQN.EQ.JMQN) THEN
              FRAC(IOCC,ICNT,IKAP,IMV) = FRAC(IOCC,ICNT,IKAP,IMV) + TMP
            ENDIF
C
          ENDDO
        ENDDO
      ENDDO
C
C     IDENTIFY CONVENTIONAL DIRAC QUANTUM NUMBERS BASED ON FRAC.
C
C     LOOP OVER POSITIVE ENERGY SPECTRUM IN PAIRS
      DO IOCC=1,NDIM-NSKP
C
C       SEARCH FOR HIGHEST POPULATED ATOMIC STATE
        PLTN = 0.0D0
        DO ICNT=1,NCNT
          DO IKAP=1,NKAP(ICNT)
            IKQN = KAPA(IKAP,ICNT)
            NMV  = 2*IABS(IKQN)
            DO IMV=1,NMV
              IF(FRAC(IOCC,ICNT,IKAP,IMV).GT.PLTN) THEN
                PLTN = FRAC(IOCC,ICNT,IKAP,IMV)
                KCNT = ICNT
                KKQN = IKQN
                KKAP = IKAP
                KMV  = IMV
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        NPR(KCNT,KKAP,KMV) = NPR(KCNT,KKAP,KMV)+1
C
C       APPLY DIRAC LABELS ACCORDING TO HIGHEST POPULATED STATE
        NQNLST(IOCC) = NPR(KCNT,KKAP,KMV)
        ICNLST(IOCC) = KCNT
        KQNLST(IOCC) = KKQN
        POPLST(IOCC) = FRAC(IOCC,KCNT,KKAP,KMV)
        IF(MOD(KMV,2).EQ.1) THEN
          MQNLST(IOCC) =-KMV
        ELSE
          MQNLST(IOCC) = KMV-1
        ENDIF
C
C     END LOOP OVER ORBITALS
      ENDDO
C
C     PRINT TITLE TO TERMINAL/FILE
20    FORMAT(1X,'Z',2X,'|',1X,'nl',6X,'Energy (au)',4X,'Pop |',4X,
     &                    '<1/r>',7X,'<r>',5X,'<r^2>',5X,'<r^3>')
21    FORMAT(1X,A,' |',I2,A,A,1X,F15.9,2X,F5.3,' |',F9.5,3(F10.5))
C
      WRITE(6,20)
      WRITE(7,20)
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     SUMMARISE RESULTS
      DO IOCC=1,IBND
C
        MOCC = IOCC+NSKP
C
        ICNT = ICNLST(IOCC)
        NQN  = NQNLST(IOCC)
        KQN  = KQNLST(IOCC)
        MQN  = MQNLST(IOCC)
        ELA  = ELMT(IZNC(ICNT))
        IF(KQN.LT.0) THEN
          LQN =-KQN-1
        ELSE
          LQN = KQN
        ENDIF
        CHL  = ELLTERM(LQN)

        IF(IOCC.LE.IBND) THEN
          PLTN = POPLST(IOCC)
        ELSE
          PLTN = 0.0D0
        ENDIF
C
        IF(HMLTN.EQ.'NORL') THEN
          PLTN = 0.5D0*PLTN
        ENDIF
C
        IF(KQN.EQ.-1) THEN
          CHS = ' '
        ELSEIF(KQN.LT.0) THEN
          CHS = '-'
        ELSEIF(KQN.GT.0) THEN
          CHS = '+'
        ENDIF
C
C       OUTPUT TO TERMINAL
        WRITE(6,21) ELA,NQN,CHL,CHS,EIGEN(MOCC),PLTN,
     &                                           (RAVG(MOCC,IK),IK=1,4)
        WRITE(7,21) ELA,NQN,CHL,CHS,EIGEN(MOCC),PLTN,
     &                                           (RAVG(MOCC,IK),IK=1,4)
C
      ENDDO
C
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
      RETURN
      END
