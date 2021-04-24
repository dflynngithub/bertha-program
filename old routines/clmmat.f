      SUBROUTINE CLMMAT(RR,IFLG,TADD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       CCCCCC  LL       MM       MM MM       MM    AA   TTTTTTTT      C
C      CC    CC LL       MMM     MMM MMM     MMM   AAAA     TT         C
C      CC       LL       MMMM   MMMM MMMM   MMMM  AA  AA    TT         C
C      CC       LL       MM MM MM MM MM MM MM MM AA    AA   TT         C
C      CC       LL       MM  MMM  MM MM  MMM  MM AAAAAAAA   TT         C
C      CC    CC LL       MM   M   MM MM   M   MM AA    AA   TT         C
C       CCCCCC  LLLLLLLL MM       MM MM       MM AA    AA   TT         C
C                                                                      C
C -------------------------------------------------------------------- C
C  CLMMAT MULTIPLIES A MOLECULAR ERI BATCH BY DENSITY ELEMENTS AND     C
C  ADDS THE CONTRIBUTIONS TO THE OPEN/CLOSED SCF COULOMB MATRICES.     C
C  DEPENDING ON THE COMBINATION OF MQN VALUES, CAN TAKE ADVANTAGE OF   C
C  INTEGRAL PERMUTATION SYMMETRIES (MINIMISING CALLS TO ERI):          C
C              ( MA, MB|-MD,-MC) =     PCD*( MA, MB| MC, MD)           C
C              (-MB,-MA| MC, MD) = PAB*    ( MA, MB| MC, MD)           C
C              ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)           C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ RR(MB2,16) - ERI'S FOR BLOCK AB, ALL 16 MQN SIGN COMBINATIONS.    C
C  ▶ NBAS(4)    - NUMBER OF BASIS FUNCTIONS IN BLOCK (ABCD).           C
C  ▶ IFLG(11)   - INTEGRAL SYMMETRY FLAGS.                             C
C**********************************************************************C
      INCLUDE 'param.h'
C
      DIMENSION IFLG(11),NBAS(4)
C
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 FOCK(MDM,MDM),OVAP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VUEH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),WDIR(MDM,MDM),
     &           WXCH(MDM,MDM),CPLE(MDM,MDM)
C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/I2EL/PAB1,PAB2,PCD1,PCD2,NA1,NB1,NC1,ND1,NA2,NB2,NC2,ND2,
     &            IBAS,JBAS,MCNT,ILIN,NADDAB,NADDCD,NBAS,IQL,IQR
      COMMON/ISCR/IMTX(MB2,11),ISCR(MB2),IMAP(MB2),IBCH,ITOG,MAXN
      COMMON/MTRX/FOCK,OVAP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VUEH,QDIR,
     &            QXCH,WDIR,WXCH,CPLE
      COMMON/SHLL/ALPH,BETA,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
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
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 101
            N = N+1
            IF(IMTX(M, 1).EQ.0) GOTO 101
C
            GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N, 1)*DENT(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N, 4)*DENT(NC2+KBAS,ND2+LBAS)
C
            GDIR(NA2+IBAS,NB2+JBAS) = GDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(N,13)*DENT(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N,16)*DENT(NC2+KBAS,ND2+LBAS)
C
101         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     2ND CASE (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(2).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 102
            N = N+1
            IF(IMTX(M, 2).EQ.0) GOTO 102
C
              GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(N, 4)*DENT(ND1+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 1)*DENT(ND2+LBAS,NC2+KBAS)
C
              GDIR(NA2+IBAS,NB2+JBAS) = GDIR(NA2+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(N,16)*DENT(ND1+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(N,13)*DENT(ND2+LBAS,NC2+KBAS)
C
102         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     3RD CASE (DIRECT):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(3).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 103
            N = N+1
            IF(IMTX(M, 3).EQ.0) GOTO 103
C
            GDIR(NC1+KBAS,ND1+LBAS) = GDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N, 1)*DENT(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N,13)*DENT(NA2+IBAS,NB2+JBAS)
C
            GDIR(NC2+KBAS,ND2+LBAS) = GDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(N, 4)*DENT(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N,16)*DENT(NA2+IBAS,NB2+JBAS)
C
103         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     4TH CASE (DIRECT):  ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(4).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 104
            N = N+1
            IF(IMTX(M, 4).EQ.0) GOTO 104
C
            GDIR(NC1+KBAS,ND1+LBAS) = GDIR(NC1+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,13)*DENT(NB1+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 1)*DENT(NB2+JBAS,NA2+IBAS)
C
            GDIR(NC2+KBAS,ND2+LBAS) = GDIR(NC2+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(N,16)*DENT(NB1+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 4)*DENT(NB2+JBAS,NA2+IBAS)
C
104         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     5TH CASE (EXCHNG):  ( MA, MB| MC, MD)
      IF(IFLG(5).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 105
            N = N+1
            IF(IMTX(M, 5).EQ.0) GOTO 105
C
            GXCH(NA1+IBAS,ND1+LBAS) = GXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 1)*DENT(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N, 7)*DENT(NC2+KBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,ND2+LBAS) = GXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(N,10)*DENT(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N,16)*DENT(NC2+KBAS,NB2+JBAS)
C
105         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     6TH CASE (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(6).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 106
            N = N+1
            IF(IMTX(M, 6).EQ.0) GOTO 106
C
            GXCH(NA1+IBAS,NC1+KBAS) = GXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 4)*DENT(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 7)*DENT(ND2+LBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,NC2+KBAS) = GXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,10)*DENT(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N,13)*DENT(ND2+LBAS,NB2+JBAS)
C
106         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     7TH CASE (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(7).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 107
            N = N+1
            IF(IMTX(M, 7).EQ.0) GOTO 107
C
            GXCH(NB1+JBAS,ND1+LBAS) = GXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,13)*DENT(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 7)*DENT(NC2+KBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,ND2+LBAS) = GXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(N,10)*DENT(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 4)*DENT(NC2+KBAS,NA2+IBAS)
C
107         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     8TH CASE (EXCHNG):  ( MB, MA| MD, MC) = PAB*    (-MA,-MB| MC, MD)
C                                           = PAB*PCD*(-MA,-MB|-MC,-MD)
      IF(IFLG(8).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 108
            N = N+1
            IF(IMTX(M, 8).EQ.0) GOTO 108
C
            GXCH(NB1+JBAS,NC1+KBAS) = GXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB1*PCD1*RR(N,16)*DENT(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(N, 7)*DENT(ND2+LBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,NC2+KBAS) = GXCH(NB2+JBAS,NC2+KBAS)
     &                     + PAB1*PCD1*RR(N, 1)*DENT(ND2+LBAS,NA2+IBAS)
     &                     + PAB2*PCD2*RR(N,10)*DENT(ND1+LBAS,NA1+IBAS)
C
108         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     9TH CASE (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(9).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 109
            N = N+1
            IF(IMTX(M, 9).EQ.0) GOTO 109
C
            GXCH(NC1+KBAS,NB1+JBAS) = GXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N, 1)*DENT(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N,10)*DENT(NA2+IBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NB2+JBAS) = GXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(N, 7)*DENT(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N,16)*DENT(NA2+IBAS,ND2+LBAS)
C
109         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     10TH CASE (EXCHNG): ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(10).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 110
            N = N+1
            IF(IMTX(M,10).EQ.0) GOTO 110
C
            GXCH(NC1+KBAS,NA1+IBAS) = GXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N,13)*DENT(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N,10)*DENT(NB2+JBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NA2+IBAS) = GXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N, 7)*DENT(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N, 4)*DENT(NB2+JBAS,ND2+LBAS)
C
110         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     11TH CASE (EXCHNG): ( MD, MC| MA, MB) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(11).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 111
            N = N+1
            IF(IMTX(M,11).EQ.0) GOTO 111
C
            GXCH(ND1+LBAS,NB1+JBAS) = GXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N, 4)*DENT(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(N,10)*DENT(NA2+IBAS,NC2+KBAS)
C
            GXCH(ND2+LBAS,NB2+JBAS) = GXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(N, 7)*DENT(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N,13)*DENT(NA2+IBAS,NC2+KBAS)
111         CONTINUE
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
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 201
            N = N+1
            IF(IMTX(M, 1).EQ.0) GOTO 201
C
            GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N, 2)*DENT(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N, 3)*DENT(NC2+KBAS,ND1+LBAS)
C
            GDIR(NA1+IBAS,NB2+JBAS) = GDIR(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N, 5)*DENT(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N, 6)*DENT(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N, 7)*DENT(NC2+KBAS,ND1+LBAS)
     &                     +           RR(N, 8)*DENT(NC2+KBAS,ND2+LBAS)
C
            GDIR(NA2+IBAS,NB1+JBAS) = GDIR(NA2+IBAS,NB1+JBAS)
     &                     +           RR(N, 9)*DENT(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N,10)*DENT(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N,11)*DENT(NC2+KBAS,ND1+LBAS)
     &                     +           RR(N,12)*DENT(NC2+KBAS,ND2+LBAS)
C
            GDIR(NA2+IBAS,NB2+JBAS) = GDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(N,14)*DENT(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N,15)*DENT(NC2+KBAS,ND1+LBAS)
C
201         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     2ND CASE (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(2).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 202
            N = N+1
            IF(IMTX(M, 2).EQ.0) GOTO 202
C
              GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 2)*DENT(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(N, 3)*DENT(ND2+LBAS,NC1+KBAS)
C
              GDIR(NA1+IBAS,NB2+JBAS) = GDIR(NA1+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 8)*DENT(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(N, 6)*DENT(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(N, 7)*DENT(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 5)*DENT(ND2+LBAS,NC2+KBAS)
C
              GDIR(NA2+IBAS,NB1+JBAS) = GDIR(NA2+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(N,12)*DENT(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(N,10)*DENT(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,11)*DENT(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 9)*DENT(ND2+LBAS,NC2+KBAS)
C
              GDIR(NA2+IBAS,NB2+JBAS) = GDIR(NA2+IBAS,NB2+JBAS)
     &                     +      PCD2*RR(N,14)*DENT(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,15)*DENT(ND2+LBAS,NC1+KBAS)
C
202         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     3RD CASE (DIRECT):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(3).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 203
            N = N+1
            IF(IMTX(M, 3).EQ.0) GOTO 203
C
            GDIR(NC1+KBAS,ND1+LBAS) = GDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N, 5)*DENT(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N, 9)*DENT(NA2+IBAS,NB1+JBAS)
C
            GDIR(NC1+KBAS,ND2+LBAS) = GDIR(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N, 2)*DENT(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N, 6)*DENT(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N,10)*DENT(NA2+IBAS,NB1+JBAS)
     &                     +           RR(N,14)*DENT(NA2+IBAS,NB2+JBAS)
C
            GDIR(NC2+KBAS,ND1+LBAS) = GDIR(NC2+KBAS,ND1+LBAS)
     &                     +           RR(N, 3)*DENT(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N, 7)*DENT(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N,11)*DENT(NA2+IBAS,NB1+JBAS)
     &                     +           RR(N,15)*DENT(NA2+IBAS,NB2+JBAS)
C
            GDIR(NC2+KBAS,ND2+LBAS) = GDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(N, 8)*DENT(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N,12)*DENT(NA2+IBAS,NB1+JBAS)
C
203         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     4TH CASE (DIRECT):  ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(4).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 204
            N = N+1
            IF(IMTX(M, 4).EQ.0) GOTO 204
C
            GDIR(NC1+KBAS,ND1+LBAS) = GDIR(NC1+KBAS,ND1+LBAS)
     &                     + PAB2*     RR(N, 5)*DENT(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(N, 9)*DENT(NB2+JBAS,NA1+IBAS)
C
            GDIR(NC1+KBAS,ND2+LBAS) = GDIR(NC1+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(N,14)*DENT(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 6)*DENT(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,10)*DENT(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 2)*DENT(NB2+JBAS,NA2+IBAS)
C
            GDIR(NC2+KBAS,ND1+LBAS) = GDIR(NC2+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,15)*DENT(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 7)*DENT(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,11)*DENT(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 3)*DENT(NB2+JBAS,NA2+IBAS)
C
            GDIR(NC2+KBAS,ND2+LBAS) = GDIR(NC2+KBAS,ND2+LBAS)
     &                     + PAB2*     RR(N, 8)*DENT(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,12)*DENT(NB2+JBAS,NA1+IBAS)
C
204         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     5TH CASE (EXCHNG):  ( MA, MB| MC, MD)
      IF(IFLG(5).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 205
            N = N+1
            IF(IMTX(M, 5).EQ.0) GOTO 205
C
            GXCH(NA1+IBAS,ND1+LBAS) = GXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 5)*DENT(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N, 3)*DENT(NC2+KBAS,NB1+JBAS)
C
            GXCH(NA1+IBAS,ND2+LBAS) = GXCH(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N, 2)*DENT(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N, 6)*DENT(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N, 4)*DENT(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N, 8)*DENT(NC2+KBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,ND1+LBAS) = GXCH(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N, 9)*DENT(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N,13)*DENT(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N,11)*DENT(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N,15)*DENT(NC2+KBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,ND2+LBAS) = GXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(N,14)*DENT(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N,12)*DENT(NC2+KBAS,NB1+JBAS)
C
205         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     6TH CASE (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(6).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 206
            N = N+1
            IF(IMTX(M, 6).EQ.0) GOTO 206
C
            GXCH(NA1+IBAS,NC1+KBAS) = GXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 8)*DENT(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(N, 3)*DENT(ND2+LBAS,NB1+JBAS)
C
            GXCH(NA1+IBAS,NC2+KBAS) = GXCH(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N, 2)*DENT(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 6)*DENT(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 1)*DENT(ND2+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N, 5)*DENT(ND2+LBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,NC1+KBAS) = GXCH(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N,12)*DENT(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N,16)*DENT(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(N,11)*DENT(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N,15)*DENT(ND2+LBAS,NB2+JBAS)
C
            GXCH(NA2+IBAS,NC2+KBAS) = GXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,14)*DENT(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 9)*DENT(ND2+LBAS,NB1+JBAS)
C
206         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     7TH CASE (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(7).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 207
            N = N+1
            IF(IMTX(M, 7).EQ.0) GOTO 207
C
            GXCH(NB1+JBAS,ND1+LBAS) = GXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N, 5)*DENT(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(N,15)*DENT(NC2+KBAS,NA1+IBAS)
C
            GXCH(NB1+JBAS,ND2+LBAS) = GXCH(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(N,14)*DENT(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 6)*DENT(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(N,16)*DENT(NC2+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 8)*DENT(NC2+KBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,ND1+LBAS) = GXCH(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N, 9)*DENT(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 1)*DENT(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,11)*DENT(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 3)*DENT(NC2+KBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,ND2+LBAS) = GXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(N, 2)*DENT(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,12)*DENT(NC2+KBAS,NA1+IBAS)
C
207         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     8TH CASE (EXCHNG):  ( MB, MA| MD, MC) = PAB*    (-MA,-MB| MC, MD)
C                                           = PAB*PCD*(-MA,-MB|-MC,-MD)
      IF(IFLG(8).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 208
            N = N+1
            IF(IMTX(M, 8).EQ.0) GOTO 208
C
            GXCH(NB1+JBAS,NC1+KBAS) = GXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB2*PCD1*RR(N, 8)*DENT(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD2*RR(N,15)*DENT(ND2+LBAS,NA1+IBAS)
C
            GXCH(NB1+JBAS,NC2+KBAS) = GXCH(NB1+JBAS,NC2+KBAS)
     &                     + PAB1*PCD2*RR(N,14)*DENT(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(N, 6)*DENT(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD1*RR(N,13)*DENT(ND2+LBAS,NA1+IBAS)
     &                     + PAB2*PCD1*RR(N, 5)*DENT(ND2+LBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,NC1+KBAS) = GXCH(NB2+JBAS,NC1+KBAS)
     &                     + PAB2*PCD1*RR(N,12)*DENT(ND1+LBAS,NA1+IBAS)
     &                     + PAB1*PCD1*RR(N, 4)*DENT(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD2*RR(N,11)*DENT(ND2+LBAS,NA1+IBAS)
     &                     + PAB1*PCD2*RR(N, 3)*DENT(ND2+LBAS,NA2+IBAS)
C
            GXCH(NB2+JBAS,NC2+KBAS) = GXCH(NB2+JBAS,NC2+KBAS)
     &                     + PAB1*PCD2*RR(N, 2)*DENT(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD1*RR(N, 9)*DENT(ND2+LBAS,NA1+IBAS)
C
208         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     9TH CASE (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(9).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 209
            N = N+1
            IF(IMTX(M, 9).EQ.0) GOTO 209
C
            GXCH(NC1+KBAS,NB1+JBAS) = GXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N, 2)*DENT(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N, 9)*DENT(NA2+IBAS,ND1+LBAS)
C
            GXCH(NC1+KBAS,NB2+JBAS) = GXCH(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N, 5)*DENT(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 6)*DENT(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N,13)*DENT(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N,14)*DENT(NA2+IBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NB1+JBAS) = GXCH(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N, 3)*DENT(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 4)*DENT(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N,11)*DENT(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N,12)*DENT(NA2+IBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NB2+JBAS) = GXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(N, 8)*DENT(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N,15)*DENT(NA2+IBAS,ND1+LBAS)
C
209         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     10TH CASE (EXCHNG): ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(10).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 210
            N = N+1
            IF(IMTX(M,10).EQ.0) GOTO 210
C
            GXCH(NC1+KBAS,NA1+IBAS) = GXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N,14)*DENT(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(N, 9)*DENT(NB2+JBAS,ND1+LBAS)
C
            GXCH(NC1+KBAS,NA2+IBAS) = GXCH(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N, 5)*DENT(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N, 6)*DENT(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(N, 1)*DENT(NB2+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N, 2)*DENT(NB2+JBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NA1+IBAS) = GXCH(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N,15)*DENT(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,16)*DENT(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(N,11)*DENT(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N,12)*DENT(NB2+JBAS,ND2+LBAS)
C
            GXCH(NC2+KBAS,NA2+IBAS) = GXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N, 8)*DENT(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(N, 3)*DENT(NB2+JBAS,ND1+LBAS)
C
210         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     11TH CASE (EXCHNG): ( MD, MC| MA, MB) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(11).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 211
            N = N+1
            IF(IMTX(M,11).EQ.0) GOTO 211
C
            GXCH(ND1+LBAS,NB1+JBAS) = GXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 2)*DENT(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(N,12)*DENT(NA2+IBAS,NC1+KBAS)
C
            GXCH(ND1+LBAS,NB2+JBAS) = GXCH(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 8)*DENT(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(N, 6)*DENT(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(N,16)*DENT(NA2+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(N,14)*DENT(NA2+IBAS,NC2+KBAS)
C
            GXCH(ND2+LBAS,NB1+JBAS) = GXCH(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 3)*DENT(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 1)*DENT(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,11)*DENT(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 9)*DENT(NA2+IBAS,NC2+KBAS)
C
            GXCH(ND2+LBAS,NB2+JBAS) = GXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 5)*DENT(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,15)*DENT(NA2+IBAS,NC1+KBAS)
C
211         CONTINUE
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
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 301
            N = N+1
            IF(IMTX(M, 1).EQ.0) GOTO 301
C
            QDIR(NA1+IBAS,NB1+JBAS) = QDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N, 1)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N, 4)*DENO(NC2+KBAS,ND2+LBAS)
C
            QDIR(NA2+IBAS,NB2+JBAS) = QDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(N,13)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N,16)*DENO(NC2+KBAS,ND2+LBAS)
C
301         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     2ND CASE (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(2).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 302
            N = N+1
            IF(IMTX(M, 2).EQ.0) GOTO 302
C
            QDIR(NA1+IBAS,NB1+JBAS) = QDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(N, 4)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 1)*DENO(ND2+LBAS,NC2+KBAS)
C
            QDIR(NA2+IBAS,NB2+JBAS) = QDIR(NA2+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(N,16)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(N,13)*DENO(ND2+LBAS,NC2+KBAS)
C
302         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     3RD CASE (DIRECT):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(3).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 303
            N = N+1
            IF(IMTX(M, 3).EQ.0) GOTO 303
C
            QDIR(NC1+KBAS,ND1+LBAS) = QDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N, 1)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N,13)*DENO(NA2+IBAS,NB2+JBAS)
C
            QDIR(NC2+KBAS,ND2+LBAS) = QDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(N, 4)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N,16)*DENO(NA2+IBAS,NB2+JBAS)
C
303         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     4TH CASE (DIRECT):  ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(4).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 304
            N = N+1
            IF(IMTX(M, 4).EQ.0) GOTO 304
C
            QDIR(NC1+KBAS,ND1+LBAS) = QDIR(NC1+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,13)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 1)*DENO(NB2+JBAS,NA2+IBAS)
C
            QDIR(NC2+KBAS,ND2+LBAS) = QDIR(NC2+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(N,16)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 4)*DENO(NB2+JBAS,NA2+IBAS)
C
304         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     5TH CASE (EXCHNG):  ( MA, MB| MC, MD)
      IF(IFLG(5).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 305
            N = N+1
            IF(IMTX(M, 5).EQ.0) GOTO 305
C
            QXCH(NA1+IBAS,ND1+LBAS) = QXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 1)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N, 7)*DENO(NC2+KBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,ND2+LBAS) = QXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(N,10)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N,16)*DENO(NC2+KBAS,NB2+JBAS)
C
305         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     6TH CASE (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(6).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 306
            N = N+1
            IF(IMTX(M, 6).EQ.0) GOTO 306
C
            QXCH(NA1+IBAS,NC1+KBAS) = QXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 4)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 7)*DENO(ND2+LBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,NC2+KBAS) = QXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,10)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N,13)*DENO(ND2+LBAS,NB2+JBAS)
C
306         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     7TH CASE (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(7).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 307
            N = N+1
            IF(IMTX(M, 7).EQ.0) GOTO 307
C
            QXCH(NB1+JBAS,ND1+LBAS) = QXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,13)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 7)*DENO(NC2+KBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,ND2+LBAS) = QXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(N,10)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 4)*DENO(NC2+KBAS,NA2+IBAS)
C
307         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     8TH CASE (EXCHNG):  ( MB, MA| MD, MC) = PAB*    (-MA,-MB| MC, MD)
C                                           = PAB*PCD*(-MA,-MB|-MC,-MD)
      IF(IFLG(8).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 308
            N = N+1
            IF(IMTX(M, 8).EQ.0) GOTO 308
C
            QXCH(NB1+JBAS,NC1+KBAS) = QXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB1*PCD1*RR(N,16)*DENO(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(N, 7)*DENO(ND2+LBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,NC2+KBAS) = QXCH(NB2+JBAS,NC2+KBAS)
     &                     + PAB2*PCD2*RR(N,10)*DENO(ND1+LBAS,NA1+IBAS)
     &                     + PAB1*PCD1*RR(N, 1)*DENO(ND2+LBAS,NA2+IBAS)
C
308         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     9TH CASE (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(9).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 309
            N = N+1
            IF(IMTX(M, 9).EQ.0) GOTO 309
C
            QXCH(NC1+KBAS,NB1+JBAS) = QXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N, 1)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N,10)*DENO(NA2+IBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NB2+JBAS) = QXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(N, 7)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N,16)*DENO(NA2+IBAS,ND2+LBAS)
C
309         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     10TH CASE (EXCHNG): ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(10).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 310
            N = N+1
            IF(IMTX(M,10).EQ.0) GOTO 310
C
            QXCH(NC1+KBAS,NA1+IBAS) = QXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N,13)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N,10)*DENO(NB2+JBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NA2+IBAS) = QXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N, 7)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N, 4)*DENO(NB2+JBAS,ND2+LBAS)
C
310         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     11TH CASE (EXCHNG): ( MD, MC| MA, MB) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(11).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 311
            N = N+1
            IF(IMTX(M,11).EQ.0) GOTO 311
C
            QXCH(ND1+LBAS,NB1+JBAS) = QXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N, 4)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(N,10)*DENO(NA2+IBAS,NC2+KBAS)
C
            QXCH(ND2+LBAS,NB2+JBAS) = QXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(N, 7)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N,13)*DENO(NA2+IBAS,NC2+KBAS)
C
311         CONTINUE
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
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 401
            N = N+1
            IF(IMTX(M, 1).EQ.0) GOTO 401
C
            QDIR(NA1+IBAS,NB1+JBAS) = QDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N, 2)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N, 3)*DENO(NC2+KBAS,ND1+LBAS)
C
            QDIR(NA1+IBAS,NB2+JBAS) = QDIR(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N, 5)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N, 6)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N, 7)*DENO(NC2+KBAS,ND1+LBAS)
     &                     +           RR(N, 8)*DENO(NC2+KBAS,ND2+LBAS)
C
            QDIR(NA2+IBAS,NB1+JBAS) = QDIR(NA2+IBAS,NB1+JBAS)
     &                     +           RR(N, 9)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N,10)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N,11)*DENO(NC2+KBAS,ND1+LBAS)
     &                     +           RR(N,12)*DENO(NC2+KBAS,ND2+LBAS)
C
            QDIR(NA2+IBAS,NB2+JBAS) = QDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(N,14)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N,15)*DENO(NC2+KBAS,ND1+LBAS)
C
401        CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     2ND CASE (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(2).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 402
            N = N+1
            IF(IMTX(M, 2).EQ.0) GOTO 402
C
            QDIR(NA1+IBAS,NB1+JBAS) = QDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 2)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(N, 3)*DENO(ND2+LBAS,NC1+KBAS)
C
            QDIR(NA1+IBAS,NB2+JBAS) = QDIR(NA1+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 8)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(N, 6)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(N, 7)*DENO(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 5)*DENO(ND2+LBAS,NC2+KBAS)
C
            QDIR(NA2+IBAS,NB1+JBAS) = QDIR(NA2+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(N,12)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(N,10)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,11)*DENO(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 9)*DENO(ND2+LBAS,NC2+KBAS)
C
            QDIR(NA2+IBAS,NB2+JBAS) = QDIR(NA2+IBAS,NB2+JBAS)
     &                     +      PCD2*RR(N,14)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,15)*DENO(ND2+LBAS,NC1+KBAS)
C
402         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     3RD CASE (DIRECT):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(3).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 403
            N = N+1
            IF(IMTX(M, 3).EQ.0) GOTO 403
C
            QDIR(NC1+KBAS,ND1+LBAS) = QDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N, 5)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N, 9)*DENO(NA2+IBAS,NB1+JBAS)
C
            QDIR(NC1+KBAS,ND2+LBAS) = QDIR(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N, 2)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N, 6)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N,10)*DENO(NA2+IBAS,NB1+JBAS)
     &                     +           RR(N,14)*DENO(NA2+IBAS,NB2+JBAS)
C
            QDIR(NC2+KBAS,ND1+LBAS) = QDIR(NC2+KBAS,ND1+LBAS)
     &                     +           RR(N, 3)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N, 7)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N,11)*DENO(NA2+IBAS,NB1+JBAS)
     &                     +           RR(N,15)*DENO(NA2+IBAS,NB2+JBAS)
C
            QDIR(NC2+KBAS,ND2+LBAS) = QDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(N, 8)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N,12)*DENO(NA2+IBAS,NB1+JBAS)
C
403         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     4TH CASE (DIRECT):  ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(4).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 404
            N = N+1
            IF(IMTX(M, 4).EQ.0) GOTO 404
C
            QDIR(NC1+KBAS,ND1+LBAS) = QDIR(NC1+KBAS,ND1+LBAS)
     &                     + PAB2*     RR(N, 5)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(N, 9)*DENO(NB2+JBAS,NA1+IBAS)
C
            QDIR(NC1+KBAS,ND2+LBAS) = QDIR(NC1+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(N,14)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 6)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,10)*DENO(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 2)*DENO(NB2+JBAS,NA2+IBAS)
C
            QDIR(NC2+KBAS,ND1+LBAS) = QDIR(NC2+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,15)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 7)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,11)*DENO(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 3)*DENO(NB2+JBAS,NA2+IBAS)
C
            QDIR(NC2+KBAS,ND2+LBAS) = QDIR(NC2+KBAS,ND2+LBAS)
     &                     + PAB2*     RR(N, 8)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,12)*DENO(NB2+JBAS,NA1+IBAS)
C
404         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     5TH CASE (EXCHNG):  ( MA, MB| MC, MD)
      IF(IFLG(5).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 405
            N = N+1
            IF(IMTX(M, 5).EQ.0) GOTO 405
C
            QXCH(NA1+IBAS,ND1+LBAS) = QXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 5)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N, 3)*DENO(NC2+KBAS,NB1+JBAS)
C
            QXCH(NA1+IBAS,ND2+LBAS) = QXCH(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N, 2)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N, 6)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N, 4)*DENO(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N, 8)*DENO(NC2+KBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,ND1+LBAS) = QXCH(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N, 9)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N,13)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N,11)*DENO(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N,15)*DENO(NC2+KBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,ND2+LBAS) = QXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(N,14)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N,12)*DENO(NC2+KBAS,NB1+JBAS)
C
405         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     6TH CASE (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(6).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 406
            N = N+1
            IF(IMTX(M, 6).EQ.0) GOTO 406
C
            QXCH(NA1+IBAS,NC1+KBAS) = QXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 8)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(N, 3)*DENO(ND2+LBAS,NB1+JBAS)
C
            QXCH(NA1+IBAS,NC2+KBAS) = QXCH(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N, 2)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 6)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 1)*DENO(ND2+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N, 5)*DENO(ND2+LBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,NC1+KBAS) = QXCH(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N,12)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N,16)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(N,11)*DENO(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N,15)*DENO(ND2+LBAS,NB2+JBAS)
C
            QXCH(NA2+IBAS,NC2+KBAS) = QXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,14)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 9)*DENO(ND2+LBAS,NB1+JBAS)
C
406         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     7TH CASE (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(7).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 407
            N = N+1
            IF(IMTX(M, 7).EQ.0) GOTO 407
C
            QXCH(NB1+JBAS,ND1+LBAS) = QXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N, 5)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(N,15)*DENO(NC2+KBAS,NA1+IBAS)
C
            QXCH(NB1+JBAS,ND2+LBAS) = QXCH(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(N,14)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 6)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(N,16)*DENO(NC2+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 8)*DENO(NC2+KBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,ND1+LBAS) = QXCH(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N, 9)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 1)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,11)*DENO(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 3)*DENO(NC2+KBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,ND2+LBAS) = QXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(N, 2)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,12)*DENO(NC2+KBAS,NA1+IBAS)
C
407         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     8TH CASE (EXCHNG):  ( MB, MA| MD, MC) = PAB*    (-MA,-MB| MC, MD)
C                                           = PAB*PCD*(-MA,-MB|-MC,-MD)
      IF(IFLG(8).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 408
            N = N+1
            IF(IMTX(M, 8).EQ.0) GOTO 408
C
            QXCH(NB1+JBAS,NC1+KBAS) = QXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB2*PCD1*RR(N, 8)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD2*RR(N,15)*DENO(ND2+LBAS,NA1+IBAS)
C
            QXCH(NB1+JBAS,NC2+KBAS) = QXCH(NB1+JBAS,NC2+KBAS)
     &                     + PAB1*PCD2*RR(N,14)*DENO(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(N, 6)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD1*RR(N,13)*DENO(ND2+LBAS,NA1+IBAS)
     &                     + PAB2*PCD1*RR(N, 5)*DENO(ND2+LBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,NC1+KBAS) = QXCH(NB2+JBAS,NC1+KBAS)
     &                     + PAB2*PCD1*RR(N,12)*DENO(ND1+LBAS,NA1+IBAS)
     &                     + PAB1*PCD1*RR(N, 4)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD2*RR(N,11)*DENO(ND2+LBAS,NA1+IBAS)
     &                     + PAB1*PCD2*RR(N, 3)*DENO(ND2+LBAS,NA2+IBAS)
C
            QXCH(NB2+JBAS,NC2+KBAS) = QXCH(NB2+JBAS,NC2+KBAS)
     &                     + PAB1*PCD2*RR(N, 2)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD1*RR(N, 9)*DENO(ND2+LBAS,NA1+IBAS)
C
408         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     9TH CASE (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(9).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 409
            N = N+1
            IF(IMTX(M,9).EQ.0) GOTO 409
C
            QXCH(NC1+KBAS,NB1+JBAS) = QXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N, 2)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N, 9)*DENO(NA2+IBAS,ND1+LBAS)
C
            QXCH(NC1+KBAS,NB2+JBAS) = QXCH(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N, 5)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 6)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N,13)*DENO(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N,14)*DENO(NA2+IBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NB1+JBAS) = QXCH(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N, 3)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 4)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N,11)*DENO(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N,12)*DENO(NA2+IBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NB2+JBAS) = QXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(N, 8)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N,15)*DENO(NA2+IBAS,ND1+LBAS)
C
409         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     10TH CASE (EXCHNG): ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(10).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 410
            N = N+1
            IF(IMTX(M,10).EQ.0) GOTO 410
C
            QXCH(NC1+KBAS,NA1+IBAS) = QXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N,14)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(N, 9)*DENO(NB2+JBAS,ND1+LBAS)
C
            QXCH(NC1+KBAS,NA2+IBAS) = QXCH(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N, 5)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N, 6)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(N, 1)*DENO(NB2+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N, 2)*DENO(NB2+JBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NA1+IBAS) = QXCH(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N,15)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,16)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(N,11)*DENO(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N,12)*DENO(NB2+JBAS,ND2+LBAS)
C
            QXCH(NC2+KBAS,NA2+IBAS) = QXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N, 8)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(N, 3)*DENO(NB2+JBAS,ND1+LBAS)
C
410         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C     11TH CASE (EXCHNG): ( MD, MC| MA, MB) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(11).EQ.1) THEN
        M = 0
        N = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
            IF(ISCR(M).EQ.0) GOTO 411
            N = N+1
            IF(IMTX(M,11).EQ.0) GOTO 411
C
            QXCH(ND1+LBAS,NB1+JBAS) = QXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 2)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(N,12)*DENO(NA2+IBAS,NC1+KBAS)
C
            QXCH(ND1+LBAS,NB2+JBAS) = QXCH(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 8)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(N, 6)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(N,16)*DENO(NA2+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(N,14)*DENO(NA2+IBAS,NC2+KBAS)
C
            QXCH(ND2+LBAS,NB1+JBAS) = QXCH(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 3)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 1)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,11)*DENO(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 9)*DENO(NA2+IBAS,NC2+KBAS)
C
            QXCH(ND2+LBAS,NB2+JBAS) = QXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 5)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,15)*DENO(NA2+IBAS,NC1+KBAS)
C
411         CONTINUE
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
      TADD = TADD+T2-T1
C
      RETURN
      END

