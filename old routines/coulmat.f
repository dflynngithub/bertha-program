      SUBROUTINE COULMAT(RR,IFLG,NBAS,TCMC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    CCCCCC   OOOOOO  UU    UU LL       MM       MM    AA   TTTTTTTT   C
C   CC    CC OO    OO UU    UU LL       MMM     MMM   AAAA     TT      C
C   CC       OO    OO UU    UU LL       MMMM   MMMM  AA  AA    TT      C
C   CC       OO    OO UU    UU LL       MM MM MM MM AA    AA   TT      C
C   CC       OO    OO UU    UU LL       MM  MMM  MM AAAAAAAA   TT      C
C   CC    CC OO    OO UU    UU LL       MM   M   MM AA    AA   TT      C
C    CCCCCC   OOOOOO   UUUUUU  LLLLLLLL MM       MM AA    AA   TT      C
C                                                                      C
C -------------------------------------------------------------------- C
C  COULMAT MULTIPLIES A MOLECULAR ERI BATCH BY DENSITY ELEMENTS AND    C
C  ADDS THE CONTRIBUTIONS TO THE OPEN/CLOSED SCF COULOMB MATRICES.     C
C  DEPENDING ON THE COMBINATION OF MQN VALUES, CAN TAKE ADVANTAGE OF   C
C  INTEGRAL PERMUTATION SYMMETRIES (MINIMISING CALLS TO ERI):          C
C              ( MA, MB|-MD,-MC) =     PCD*( MA, MB| MC, MD)           C
C              (-MB,-MA| MC, MD) = PAB*    ( MA, MB| MC, MD)           C
C              ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)           C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    RR(MB2,16) - ERI'S FOR BLOCK AB, ALL 16 MQN SIGN COMBINATIONS.    C
C    IFLG(11)   - INTEGRAL SYMMETRY FLAGS.                             C
C    NBAS(4)    - NUMBER OF BASIS FUNCTIONS IN BLOCK (ABCD).           C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=15,MKP=9)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
C
      DIMENSION IFLG(11),NBAS(4)
C
      COMMON/BLOC/PAB1,PAB2,PCD1,PCD2,NA1,NB1,NC1,ND1,NA2,NB2,NC2,ND2,
     &            IBAS,JBAS,ILIN,NADDAB,NADDCD
      COMMON/DENS/DENC,DENO,DENT
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,IEQS,IERC,IPAR,ICOR,ILEV
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
C
C     TIME AT START OF ROUTINE
      CALL CPU_TIME(T1)
C
C**********************************************************************C
C     CLOSED-SHELL CONTRIBUTIONS...                                    C
C**********************************************************************C
C
C     1ST CASE (DIRECT):  ( MA, MB| MC, MD)
      IF(IFLG(1).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 1)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 4)*DENC(NC2+KBAS,ND2+LBAS)
C
            GDIR(NA2+IBAS,NB2+JBAS) = GDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(M,13)*DENC(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,16)*DENC(NC2+KBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     2ND CASE (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(2).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
              GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENC(ND2+LBAS,NC2+KBAS)
C
              GDIR(NA2+IBAS,NB2+JBAS) = GDIR(NA2+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(M,16)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENC(ND2+LBAS,NC2+KBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     3RD CASE (DIRECT):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(3).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NC1+KBAS,ND1+LBAS) = GDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 1)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M,13)*DENC(NA2+IBAS,NB2+JBAS)
C
            GDIR(NC2+KBAS,ND2+LBAS) = GDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(M, 4)*DENC(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M,16)*DENC(NA2+IBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     4TH CASE (DIRECT):  ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(4).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NC1+KBAS,ND1+LBAS) = GDIR(NC1+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,13)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENC(NB2+JBAS,NA2+IBAS)
C
            GDIR(NC2+KBAS,ND2+LBAS) = GDIR(NC2+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,16)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENC(NB2+JBAS,NA2+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     5TH CASE (EXCHNG):  ( MA, MB| MC, MD)
      IF(IFLG(5).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NA1+IBAS,ND1+LBAS) = GXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 1)*DENC(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 7)*DENC(NC2+KBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,ND2+LBAS) = GXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(M,10)*DENC(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M,16)*DENC(NC2+KBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     6TH CASE (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(6).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NA1+IBAS,NC1+KBAS) = GXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 4)*DENC(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 7)*DENC(ND2+LBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,NC2+KBAS) = GXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,10)*DENC(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,13)*DENC(ND2+LBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     7TH CASE (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(7).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NB1+JBAS,ND1+LBAS) = GXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,13)*DENC(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENC(NC2+KBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,ND2+LBAS) = GXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M,10)*DENC(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENC(NC2+KBAS,NA2+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     8TH CASE (EXCHNG):  ( MB, MA| MD, MC) = PAB*    (-MA,-MB| MC, MD)
C                                           = PAB*PCD*(-MA,-MB|-MC,-MD)
      IF(IFLG(8).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NB1+JBAS,NC1+KBAS) = GXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB1*PCD1*RR(M,16)*DENC(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(M, 7)*DENC(ND2+LBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,NC2+KBAS) = GXCH(NB2+JBAS,NC2+KBAS)
     &                     + PAB1*PCD1*RR(M, 1)*DENC(ND2+LBAS,NA2+IBAS)
     &                     + PAB2*PCD2*RR(M,10)*DENC(ND1+LBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     9TH CASE (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(9).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NC1+KBAS,NB1+JBAS) = GXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 1)*DENC(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M,10)*DENC(NA2+IBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NB2+JBAS) = GXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(M, 7)*DENC(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M,16)*DENC(NA2+IBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     10TH CASE (EXCHNG): ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(10).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NC1+KBAS,NA1+IBAS) = GXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,13)*DENC(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M,10)*DENC(NB2+JBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NA2+IBAS) = GXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 7)*DENC(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M, 4)*DENC(NB2+JBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     11TH CASE (EXCHNG): ( MD, MC| MA, MB) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(11).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(ND1+LBAS,NB1+JBAS) = GXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENC(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENC(NA2+IBAS,NC2+KBAS)
C
            GXCH(ND2+LBAS,NB2+JBAS) = GXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M, 7)*DENC(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENC(NA2+IBAS,NC2+KBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     ADD MORE CONTRIBUTIONS FOR GENERAL MOLECULAR BATCH
      IF(ILIN.EQ.1) GOTO 5002
C
C     1ST CASE (DIRECT):  ( MA, MB| MC, MD)
      IF(IFLG(1).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 2)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 3)*DENC(NC2+KBAS,ND1+LBAS)
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
     &                     +           RR(M,14)*DENC(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,15)*DENC(NC2+KBAS,ND1+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     2ND CASE (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(2).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
              GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 2)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 3)*DENC(ND2+LBAS,NC1+KBAS)
C
              GDIR(NA1+IBAS,NB2+JBAS) = GDIR(NA1+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 8)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 7)*DENC(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 5)*DENC(ND2+LBAS,NC2+KBAS)
C
              GDIR(NA2+IBAS,NB1+JBAS) = GDIR(NA2+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,12)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENC(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENC(ND2+LBAS,NC2+KBAS)
C
              GDIR(NA2+IBAS,NB2+JBAS) = GDIR(NA2+IBAS,NB2+JBAS)
     &                     +      PCD2*RR(M,14)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENC(ND2+LBAS,NC1+KBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     3RD CASE (DIRECT):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(3).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NC1+KBAS,ND1+LBAS) = GDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 5)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 9)*DENC(NA2+IBAS,NB1+JBAS)
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
     &                     +           RR(M, 8)*DENC(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,12)*DENC(NA2+IBAS,NB1+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     4TH CASE (DIRECT):  ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(4).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NC1+KBAS,ND1+LBAS) = GDIR(NC1+KBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 5)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 9)*DENC(NB2+JBAS,NA1+IBAS)
C
            GDIR(NC1+KBAS,ND2+LBAS) = GDIR(NC1+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,14)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,10)*DENC(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 2)*DENC(NB2+JBAS,NA2+IBAS)
C
            GDIR(NC2+KBAS,ND1+LBAS) = GDIR(NC2+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,15)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENC(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENC(NB2+JBAS,NA2+IBAS)
C
            GDIR(NC2+KBAS,ND2+LBAS) = GDIR(NC2+KBAS,ND2+LBAS)
     &                     + PAB2*     RR(M, 8)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENC(NB2+JBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     5TH CASE (EXCHNG):  ( MA, MB| MC, MD)
      IF(IFLG(5).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NA1+IBAS,ND1+LBAS) = GXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 5)*DENC(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 3)*DENC(NC2+KBAS,NB1+JBAS)
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
     &                     +           RR(M,14)*DENC(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M,12)*DENC(NC2+KBAS,NB1+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     6TH CASE (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(6).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NA1+IBAS,NC1+KBAS) = GXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 8)*DENC(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M, 3)*DENC(ND2+LBAS,NB1+JBAS)
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
     &                     +      PCD2*RR(M,14)*DENC(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 9)*DENC(ND2+LBAS,NB1+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     7TH CASE (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(7).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NB1+JBAS,ND1+LBAS) = GXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 5)*DENC(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(M,15)*DENC(NC2+KBAS,NA1+IBAS)
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
     &                     + PAB1*     RR(M, 2)*DENC(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENC(NC2+KBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     8TH CASE (EXCHNG):  ( MB, MA| MD, MC) = PAB*    (-MA,-MB| MC, MD)
C                                           = PAB*PCD*(-MA,-MB|-MC,-MD)
      IF(IFLG(8).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NB1+JBAS,NC1+KBAS) = GXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB2*PCD1*RR(M, 8)*DENC(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD2*RR(M,15)*DENC(ND2+LBAS,NA1+IBAS)
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
     &                     + PAB1*PCD2*RR(M, 2)*DENC(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD1*RR(M, 9)*DENC(ND2+LBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     9TH CASE (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(9).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NC1+KBAS,NB1+JBAS) = GXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 2)*DENC(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M, 9)*DENC(NA2+IBAS,ND1+LBAS)
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
     &                     +           RR(M, 8)*DENC(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,15)*DENC(NA2+IBAS,ND1+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     10TH CASE (EXCHNG): ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(10).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(NC1+KBAS,NA1+IBAS) = GXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,14)*DENC(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M, 9)*DENC(NB2+JBAS,ND1+LBAS)
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
     &                     + PAB2*     RR(M, 8)*DENC(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M, 3)*DENC(NB2+JBAS,ND1+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     11TH CASE (EXCHNG): ( MD, MC| MA, MB) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(11).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GXCH(ND1+LBAS,NB1+JBAS) = GXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 2)*DENC(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(M,12)*DENC(NA2+IBAS,NC1+KBAS)
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
     &                     +      PCD1*RR(M, 5)*DENC(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENC(NA2+IBAS,NC1+KBAS)
C
          ENDDO
        ENDDO
      ENDIF

C
5002  CONTINUE
C
C**********************************************************************C
C     OPEN-SHELL CONTRIBUTIONS...                                      C
C**********************************************************************C
C
      IF(NOPN.EQ.0) GOTO 5100
C
C     1ST CASE (DIRECT):  ( MA, MB| MC, MD)
      IF(IFLG(1).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NA1+IBAS,NB1+JBAS) = QDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 1)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 4)*DENO(NC2+KBAS,ND2+LBAS)
C
            QDIR(NA2+IBAS,NB2+JBAS) = QDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(M,13)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,16)*DENO(NC2+KBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     2ND CASE (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(2).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NA1+IBAS,NB1+JBAS) = QDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENO(ND2+LBAS,NC2+KBAS)
C
            QDIR(NA2+IBAS,NB2+JBAS) = QDIR(NA2+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(M,16)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENO(ND2+LBAS,NC2+KBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     3RD CASE (DIRECT):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(3).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NC1+KBAS,ND1+LBAS) = QDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 1)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M,13)*DENO(NA2+IBAS,NB2+JBAS)
C
            QDIR(NC2+KBAS,ND2+LBAS) = QDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(M, 4)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M,16)*DENO(NA2+IBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     4TH CASE (DIRECT):  ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(4).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NC1+KBAS,ND1+LBAS) = QDIR(NC1+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,13)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENO(NB2+JBAS,NA2+IBAS)
C
            QDIR(NC2+KBAS,ND2+LBAS) = QDIR(NC2+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,16)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENO(NB2+JBAS,NA2+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     5TH CASE (EXCHNG):  ( MA, MB| MC, MD)
      IF(IFLG(5).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NA1+IBAS,ND1+LBAS) = QXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 1)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 7)*DENO(NC2+KBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,ND2+LBAS) = QXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(M,10)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M,16)*DENO(NC2+KBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     6TH CASE (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(6).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NA1+IBAS,NC1+KBAS) = QXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 4)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 7)*DENO(ND2+LBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,NC2+KBAS) = QXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,10)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,13)*DENO(ND2+LBAS,NB2+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     7TH CASE (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(7).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NB1+JBAS,ND1+LBAS) = QXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,13)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENO(NC2+KBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,ND2+LBAS) = QXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M,10)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENO(NC2+KBAS,NA2+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     8TH CASE (EXCHNG):  ( MB, MA| MD, MC) = PAB*    (-MA,-MB| MC, MD)
C                                           = PAB*PCD*(-MA,-MB|-MC,-MD)
      IF(IFLG(8).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NB1+JBAS,NC1+KBAS) = QXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB1*PCD1*RR(M,16)*DENO(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(M, 7)*DENO(ND2+LBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,NC2+KBAS) = QXCH(NB2+JBAS,NC2+KBAS)
     &                     + PAB2*PCD2*RR(M,10)*DENO(ND1+LBAS,NA1+IBAS)
     &                     + PAB1*PCD1*RR(M, 1)*DENO(ND2+LBAS,NA2+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     9TH CASE (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(9).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NC1+KBAS,NB1+JBAS) = QXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 1)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M,10)*DENO(NA2+IBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NB2+JBAS) = QXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(M, 7)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M,16)*DENO(NA2+IBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     10TH CASE (EXCHNG): ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(10).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NC1+KBAS,NA1+IBAS) = QXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,13)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M,10)*DENO(NB2+JBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NA2+IBAS) = QXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 7)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M, 4)*DENO(NB2+JBAS,ND2+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     11TH CASE (EXCHNG): ( MD, MC| MA, MB) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(11).EQ.1) THEN
        M = 0
        DO LBAS=1,NBAS(4)
          DO KBAS=1,NBAS(3)
            M = M+1
C
            QXCH(ND1+LBAS,NB1+JBAS) = QXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENO(NA2+IBAS,NC2+KBAS)
C
            QXCH(ND2+LBAS,NB2+JBAS) = QXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M, 7)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENO(NA2+IBAS,NC2+KBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     ADD MORE CONTRIBUTIONS FOR GENERAL MOLECULAR BATCH
      IF(ILIN.EQ.1) GOTO 5003
C
C     1ST CASE (DIRECT):  ( MA, MB| MC, MD)
      IF(IFLG(1).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NA1+IBAS,NB1+JBAS) = QDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 2)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 3)*DENO(NC2+KBAS,ND1+LBAS)
C
            QDIR(NA1+IBAS,NB2+JBAS) = QDIR(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 5)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 6)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 7)*DENO(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 8)*DENO(NC2+KBAS,ND2+LBAS)
C
            QDIR(NA2+IBAS,NB1+JBAS) = QDIR(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M, 9)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,10)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,11)*DENO(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M,12)*DENO(NC2+KBAS,ND2+LBAS)
C
            QDIR(NA2+IBAS,NB2+JBAS) = QDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(M,14)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,15)*DENO(NC2+KBAS,ND1+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     2ND CASE (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(2).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NA1+IBAS,NB1+JBAS) = QDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 2)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 3)*DENO(ND2+LBAS,NC1+KBAS)
C
            QDIR(NA1+IBAS,NB2+JBAS) = QDIR(NA1+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 8)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 7)*DENO(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 5)*DENO(ND2+LBAS,NC2+KBAS)
C
            QDIR(NA2+IBAS,NB1+JBAS) = QDIR(NA2+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,12)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENO(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENO(ND2+LBAS,NC2+KBAS)
C
            QDIR(NA2+IBAS,NB2+JBAS) = QDIR(NA2+IBAS,NB2+JBAS)
     &                     +      PCD2*RR(M,14)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENO(ND2+LBAS,NC1+KBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     3RD CASE (DIRECT):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(3).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NC1+KBAS,ND1+LBAS) = QDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 5)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 9)*DENO(NA2+IBAS,NB1+JBAS)
C
            QDIR(NC1+KBAS,ND2+LBAS) = QDIR(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 2)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 6)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,10)*DENO(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,14)*DENO(NA2+IBAS,NB2+JBAS)
C
            QDIR(NC2+KBAS,ND1+LBAS) = QDIR(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 3)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 7)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,11)*DENO(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,15)*DENO(NA2+IBAS,NB2+JBAS)
C
            QDIR(NC2+KBAS,ND2+LBAS) = QDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(M, 8)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,12)*DENO(NA2+IBAS,NB1+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     4TH CASE (DIRECT):  ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(4).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QDIR(NC1+KBAS,ND1+LBAS) = QDIR(NC1+KBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 5)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 9)*DENO(NB2+JBAS,NA1+IBAS)
C
            QDIR(NC1+KBAS,ND2+LBAS) = QDIR(NC1+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,14)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,10)*DENO(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 2)*DENO(NB2+JBAS,NA2+IBAS)
C
            QDIR(NC2+KBAS,ND1+LBAS) = QDIR(NC2+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,15)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENO(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENO(NB2+JBAS,NA2+IBAS)
C
            QDIR(NC2+KBAS,ND2+LBAS) = QDIR(NC2+KBAS,ND2+LBAS)
     &                     + PAB2*     RR(M, 8)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENO(NB2+JBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     5TH CASE (EXCHNG):  ( MA, MB| MC, MD)
      IF(IFLG(5).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NA1+IBAS,ND1+LBAS) = QXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 5)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 3)*DENO(NC2+KBAS,NB1+JBAS)
C
            QXCH(NA1+IBAS,ND2+LBAS) = QXCH(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M, 2)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 6)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 4)*DENO(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M, 8)*DENO(NC2+KBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,ND1+LBAS) = QXCH(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M, 9)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M,13)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M,11)*DENO(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M,15)*DENO(NC2+KBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,ND2+LBAS) = QXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(M,14)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M,12)*DENO(NC2+KBAS,NB1+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     6TH CASE (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(6).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NA1+IBAS,NC1+KBAS) = QXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 8)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M, 3)*DENO(ND2+LBAS,NB1+JBAS)
C
            QXCH(NA1+IBAS,NC2+KBAS) = QXCH(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 2)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 6)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 1)*DENO(ND2+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 5)*DENO(ND2+LBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,NC1+KBAS) = QXCH(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,12)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,16)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(M,11)*DENO(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M,15)*DENO(ND2+LBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,NC2+KBAS) = QXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,14)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 9)*DENO(ND2+LBAS,NB1+JBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     7TH CASE (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(7).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NB1+JBAS,ND1+LBAS) = QXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 5)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(M,15)*DENO(NC2+KBAS,NA1+IBAS)
C
            QXCH(NB1+JBAS,ND2+LBAS) = QXCH(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,14)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(M,16)*DENO(NC2+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 8)*DENO(NC2+KBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,ND1+LBAS) = QXCH(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 9)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENO(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENO(NC2+KBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,ND2+LBAS) = QXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M, 2)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENO(NC2+KBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     8TH CASE (EXCHNG):  ( MB, MA| MD, MC) = PAB*    (-MA,-MB| MC, MD)
C                                           = PAB*PCD*(-MA,-MB|-MC,-MD)
      IF(IFLG(8).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NB1+JBAS,NC1+KBAS) = QXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB2*PCD1*RR(M, 8)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD2*RR(M,15)*DENO(ND2+LBAS,NA1+IBAS)
C
            QXCH(NB1+JBAS,NC2+KBAS) = QXCH(NB1+JBAS,NC2+KBAS)
     &                     + PAB1*PCD2*RR(M,14)*DENO(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(M, 6)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD1*RR(M,13)*DENO(ND2+LBAS,NA1+IBAS)
     &                     + PAB2*PCD1*RR(M, 5)*DENO(ND2+LBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,NC1+KBAS) = QXCH(NB2+JBAS,NC1+KBAS)
     &                     + PAB2*PCD1*RR(M,12)*DENO(ND1+LBAS,NA1+IBAS)
     &                     + PAB1*PCD1*RR(M, 4)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD2*RR(M,11)*DENO(ND2+LBAS,NA1+IBAS)
     &                     + PAB1*PCD2*RR(M, 3)*DENO(ND2+LBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,NC2+KBAS) = QXCH(NB2+JBAS,NC2+KBAS)
     &                     + PAB1*PCD2*RR(M, 2)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD1*RR(M, 9)*DENO(ND2+LBAS,NA1+IBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     9TH CASE (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(9).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NC1+KBAS,NB1+JBAS) = QXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(M, 2)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M, 9)*DENO(NA2+IBAS,ND1+LBAS)
C
            QXCH(NC1+KBAS,NB2+JBAS) = QXCH(NC1+KBAS,NB2+JBAS)
     &                     +           RR(M, 5)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 6)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,13)*DENO(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M,14)*DENO(NA2+IBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NB1+JBAS) = QXCH(NC2+KBAS,NB1+JBAS)
     &                     +           RR(M, 3)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(M, 4)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,11)*DENO(NA2+IBAS,ND1+LBAS)
     &                     +           RR(M,12)*DENO(NA2+IBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NB2+JBAS) = QXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(M, 8)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(M,15)*DENO(NA2+IBAS,ND1+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     10TH CASE (EXCHNG): ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(10).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            QXCH(NC1+KBAS,NA1+IBAS) = QXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,14)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M, 9)*DENO(NB2+JBAS,ND1+LBAS)
C
            QXCH(NC1+KBAS,NA2+IBAS) = QXCH(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 5)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M, 6)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M, 1)*DENO(NB2+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M, 2)*DENO(NB2+JBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NA1+IBAS) = QXCH(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(M,15)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,16)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(M,11)*DENO(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(M,12)*DENO(NB2+JBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NA2+IBAS) = QXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 8)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(M, 3)*DENO(NB2+JBAS,ND1+LBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     11TH CASE (EXCHNG): ( MD, MC| MA, MB) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(11).EQ.1) THEN
        M = 0
        DO LBAS=1,NBAS(4)
          DO KBAS=1,NBAS(3)
            M = M+1
C
            QXCH(ND1+LBAS,NB1+JBAS) = QXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 2)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(M,12)*DENO(NA2+IBAS,NC1+KBAS)
C
            QXCH(ND1+LBAS,NB2+JBAS) = QXCH(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 8)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(M,16)*DENO(NA2+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,14)*DENO(NA2+IBAS,NC2+KBAS)
C
            QXCH(ND2+LBAS,NB1+JBAS) = QXCH(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(M, 3)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENO(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENO(NA2+IBAS,NC2+KBAS)
C
            QXCH(ND2+LBAS,NB2+JBAS) = QXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 5)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENO(NA2+IBAS,NC1+KBAS)
C
          ENDDO
        ENDDO
      ENDIF
C
C     SKIP POINT FOR LINEAR BATCH
5003  CONTINUE
C     SKIPPING POINT FOR CLOSED SYSTEMS
5100  CONTINUE
C
C     TIME AT END OF ROUTINE
      CALL CPU_TIME(T2)
C
C     ADD TO THE LINKED TIME INDEX
      TCMC = TCMC+T2-T1
C
      RETURN
      END

