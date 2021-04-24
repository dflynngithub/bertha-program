      SUBROUTINE BRTMAT(RR,IFLG,TADD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        BBBBBBB  RRRRRRR TTTTTTTT MM       MM    AA   TTTTTTTT        C
C        BB    BB RR    RR   TT    MMM     MMM   AAAA     TT           C
C        BB    BB RR    RR   TT    MMMM   MMMM  AA  AA    TT           C
C        BBBBBBB  RR    RR   TT    MM MM MM MM AA    AA   TT           C
C        BB    BB RRRRRRR    TT    MM  MMM  MM AAAAAAAA   TT           C
C        BB    BB RR    RR   TT    MM   M   MM AA    AA   TT           C
C        BBBBBBB  RR    RR   TT    MM       MM AA    AA   TT           C
C                                                                      C
C -------------------------------------------------------------------- C
C  BRTMAT MULTIPLIES A MOLECULAR ERI BATCH BY DENSITY ELEMENTS AND     C
C  ADDS THE CONTRIBUTIONS TO THE OPEN/CLOSED SCF BREIT MATRICES.       C
C  DEPENDING ON THE COMBINATION OF MQN VALUES, CAN TAKE ADVANTAGE OF   C
C  INTEGRAL PERMUTATION SYMMETRIES (MINIMISING CALLS TO BII):          C
C              ( MA, MB|-MD,-MC) =     PCD*( MA, MB| MC, MD)           C
C              (-MB,-MA| MC, MD) = PAB*    ( MA, MB| MC, MD)           C
C              ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)           C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ RR(MB2,16) - ERI'S FOR BLOCK AB, ALL 16 MQN SIGN COMBINATIONS.    C
C  ▶ NBAS(4)    - NUMBER OF BASIS FUNCTIONS IN BLOCK (ABCD).           C
C  ▶ IFLG(11)   - INTEGRAL SYMMETRY FLAGS.                             C
C -------------------------------------------------------------------- C
C  THIS VERSION HAS NO MOLECULAR GEOMETRY SHORTCUTS.                   C
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
     &            IBAS,JBAS,ILIN,NADDAB,NADDCD,NBAS
      COMMON/ISCR/IMTX(MB2,11),ISCR(MB2),IMAP(MB2),IBCH,ITOG,MAXN
      COMMON/MTRX/FOCK,OVAP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VUEH,QDIR,
     &            QXCH,WDIR,WXCH,CPLE
      COMMON/SHLL/ALPH,BETA,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
C
C     TIME AT START OF ROUTINE
      CALL CPU_TIME(T1)
C
C**********************************************************************C
C     CLOSED-SHELL CONTRIBUTIONS... (NO DIRECT TERMS)                  C
C**********************************************************************C
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
            BXCH(NA1+IBAS,ND1+LBAS) = BXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 1)*DENT(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N, 5)*DENT(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N, 3)*DENT(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N, 7)*DENT(NC2+KBAS,NB2+JBAS)
C
            BXCH(NA1+IBAS,ND2+LBAS) = BXCH(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N, 2)*DENT(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N, 6)*DENT(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N, 4)*DENT(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N, 8)*DENT(NC2+KBAS,NB2+JBAS)
C
            BXCH(NA2+IBAS,ND1+LBAS) = BXCH(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N, 9)*DENT(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N,13)*DENT(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N,11)*DENT(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N,15)*DENT(NC2+KBAS,NB2+JBAS)
C
            BXCH(NA2+IBAS,ND2+LBAS) = BXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(N,10)*DENT(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N,14)*DENT(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N,12)*DENT(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N,16)*DENT(NC2+KBAS,NB2+JBAS)
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
            BXCH(NA1+IBAS,NC1+KBAS) = BXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 4)*DENT(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N, 8)*DENT(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(N, 3)*DENT(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 7)*DENT(ND2+LBAS,NB2+JBAS)
C
            BXCH(NA1+IBAS,NC2+KBAS) = BXCH(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N, 2)*DENT(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 6)*DENT(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 1)*DENT(ND2+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N, 5)*DENT(ND2+LBAS,NB2+JBAS)
C
            BXCH(NA2+IBAS,NC1+KBAS) = BXCH(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N,12)*DENT(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N,16)*DENT(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(N,11)*DENT(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N,15)*DENT(ND2+LBAS,NB2+JBAS)
C
            BXCH(NA2+IBAS,NC2+KBAS) = BXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,10)*DENT(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N,14)*DENT(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 9)*DENT(ND2+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N,13)*DENT(ND2+LBAS,NB2+JBAS)
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
            BXCH(NB1+JBAS,ND1+LBAS) = BXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,13)*DENT(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 5)*DENT(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(N,15)*DENT(NC2+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 7)*DENT(NC2+KBAS,NA2+IBAS)
C
            BXCH(NB1+JBAS,ND2+LBAS) = BXCH(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(N,14)*DENT(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 6)*DENT(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(N,16)*DENT(NC2+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 8)*DENT(NC2+KBAS,NA2+IBAS)
C
            BXCH(NB2+JBAS,ND1+LBAS) = BXCH(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N, 9)*DENT(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 1)*DENT(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,11)*DENT(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 3)*DENT(NC2+KBAS,NA2+IBAS)
C
            BXCH(NB2+JBAS,ND2+LBAS) = BXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(N,10)*DENT(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 2)*DENT(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,12)*DENT(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 4)*DENT(NC2+KBAS,NA2+IBAS)
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
            BXCH(NB1+JBAS,NC1+KBAS) = BXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB1*PCD1*RR(N,16)*DENT(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD1*RR(N, 8)*DENT(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD2*RR(N,15)*DENT(ND2+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(N, 7)*DENT(ND2+LBAS,NA2+IBAS)
C
            BXCH(NB1+JBAS,NC2+KBAS) = BXCH(NB1+JBAS,NC2+KBAS)
     &                     + PAB1*PCD2*RR(N,14)*DENT(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(N, 6)*DENT(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD1*RR(N,13)*DENT(ND2+LBAS,NA1+IBAS)
     &                     + PAB2*PCD1*RR(N, 5)*DENT(ND2+LBAS,NA2+IBAS)
C
            BXCH(NB2+JBAS,NC1+KBAS) = BXCH(NB2+JBAS,NC1+KBAS)
     &                     + PAB2*PCD1*RR(N,12)*DENT(ND1+LBAS,NA1+IBAS)
     &                     + PAB1*PCD1*RR(N, 4)*DENT(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD2*RR(N,11)*DENT(ND2+LBAS,NA1+IBAS)
     &                     + PAB1*PCD2*RR(N, 3)*DENT(ND2+LBAS,NA2+IBAS)
C
            BXCH(NB2+JBAS,NC2+KBAS) = BXCH(NB2+JBAS,NC2+KBAS)
     &                     + PAB1*PCD1*RR(N, 1)*DENT(ND2+LBAS,NA2+IBAS)
     &                     + PAB1*PCD2*RR(N, 2)*DENT(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD1*RR(N, 9)*DENT(ND2+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(N,10)*DENT(ND1+LBAS,NA1+IBAS)
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
            BXCH(NC1+KBAS,NB1+JBAS) = BXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N, 1)*DENT(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 2)*DENT(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N, 9)*DENT(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N,10)*DENT(NA2+IBAS,ND2+LBAS)
C
            BXCH(NC1+KBAS,NB2+JBAS) = BXCH(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N, 5)*DENT(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 6)*DENT(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N,13)*DENT(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N,14)*DENT(NA2+IBAS,ND2+LBAS)
C
            BXCH(NC2+KBAS,NB1+JBAS) = BXCH(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N, 3)*DENT(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 4)*DENT(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N,11)*DENT(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N,12)*DENT(NA2+IBAS,ND2+LBAS)
C
            BXCH(NC2+KBAS,NB2+JBAS) = BXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(N, 7)*DENT(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 8)*DENT(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N,15)*DENT(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N,16)*DENT(NA2+IBAS,ND2+LBAS)
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
            BXCH(NC1+KBAS,NA1+IBAS) = BXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N,13)*DENT(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,14)*DENT(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(N, 9)*DENT(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N,10)*DENT(NB2+JBAS,ND2+LBAS)
C
            BXCH(NC1+KBAS,NA2+IBAS) = BXCH(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N, 5)*DENT(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N, 6)*DENT(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(N, 1)*DENT(NB2+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N, 2)*DENT(NB2+JBAS,ND2+LBAS)
C
            BXCH(NC2+KBAS,NA1+IBAS) = BXCH(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N,15)*DENT(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,16)*DENT(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(N,11)*DENT(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N,12)*DENT(NB2+JBAS,ND2+LBAS)
C
            BXCH(NC2+KBAS,NA2+IBAS) = BXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N, 7)*DENT(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N, 8)*DENT(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(N, 3)*DENT(NB2+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N, 4)*DENT(NB2+JBAS,ND2+LBAS)
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
            BXCH(ND1+LBAS,NB1+JBAS) = BXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N, 4)*DENT(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(N, 2)*DENT(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(N,12)*DENT(NA2+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(N,10)*DENT(NA2+IBAS,NC2+KBAS)
C
            BXCH(ND1+LBAS,NB2+JBAS) = BXCH(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 8)*DENT(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(N, 6)*DENT(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(N,16)*DENT(NA2+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(N,14)*DENT(NA2+IBAS,NC2+KBAS)
C
            BXCH(ND2+LBAS,NB1+JBAS) = BXCH(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 3)*DENT(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 1)*DENT(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,11)*DENT(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 9)*DENT(NA2+IBAS,NC2+KBAS)
C
            BXCH(ND2+LBAS,NB2+JBAS) = BXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(N, 7)*DENT(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 5)*DENT(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,15)*DENT(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N,13)*DENT(NA2+IBAS,NC2+KBAS)
C
211         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
C**********************************************************************C
C     OPEN-SHELL CONTRIBUTIONS... (NOTE: WILL NEED WDIR/WXCH INSTEAD)  C
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
            IF(ISCR(M).EQ.0) GOTO 401
            N = N+1
            IF(IMTX(M, 1).EQ.0) GOTO 401
C
            WDIR(NA1+IBAS,NB1+JBAS) = WDIR(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N, 1)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N, 2)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N, 3)*DENO(NC2+KBAS,ND1+LBAS)
     &                     +           RR(N, 4)*DENO(NC2+KBAS,ND2+LBAS)
C
            WDIR(NA1+IBAS,NB2+JBAS) = WDIR(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N, 5)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N, 6)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N, 7)*DENO(NC2+KBAS,ND1+LBAS)
     &                     +           RR(N, 8)*DENO(NC2+KBAS,ND2+LBAS)
C
            WDIR(NA2+IBAS,NB1+JBAS) = WDIR(NA2+IBAS,NB1+JBAS)
     &                     +           RR(N, 9)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N,10)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N,11)*DENO(NC2+KBAS,ND1+LBAS)
     &                     +           RR(N,12)*DENO(NC2+KBAS,ND2+LBAS)
C
            WDIR(NA2+IBAS,NB2+JBAS) = WDIR(NA2+IBAS,NB2+JBAS)
     &                     +           RR(N,13)*DENO(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N,14)*DENO(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N,15)*DENO(NC2+KBAS,ND1+LBAS)
     &                     +           RR(N,16)*DENO(NC2+KBAS,ND2+LBAS)
C
401         CONTINUE
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
            WDIR(NA1+IBAS,NB1+JBAS) = WDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(N, 4)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(N, 2)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(N, 3)*DENO(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 1)*DENO(ND2+LBAS,NC2+KBAS)
C
            WDIR(NA1+IBAS,NB2+JBAS) = WDIR(NA1+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 8)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(N, 6)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(N, 7)*DENO(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 5)*DENO(ND2+LBAS,NC2+KBAS)
C
            WDIR(NA2+IBAS,NB1+JBAS) = WDIR(NA2+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(N,12)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(N,10)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,11)*DENO(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 9)*DENO(ND2+LBAS,NC2+KBAS)
C
            WDIR(NA2+IBAS,NB2+JBAS) = WDIR(NA2+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(N,16)*DENO(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(N,14)*DENO(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,15)*DENO(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(N,13)*DENO(ND2+LBAS,NC2+KBAS)
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
            WDIR(NC1+KBAS,ND1+LBAS) = WDIR(NC1+KBAS,ND1+LBAS)
     &                     +           RR(N, 1)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N, 5)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N, 9)*DENO(NA2+IBAS,NB1+JBAS)
     &                     +           RR(N,13)*DENO(NA2+IBAS,NB2+JBAS)
C
            WDIR(NC1+KBAS,ND2+LBAS) = WDIR(NC1+KBAS,ND2+LBAS)
     &                     +           RR(N, 2)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N, 6)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N,10)*DENO(NA2+IBAS,NB1+JBAS)
     &                     +           RR(N,14)*DENO(NA2+IBAS,NB2+JBAS)
C
            WDIR(NC2+KBAS,ND1+LBAS) = WDIR(NC2+KBAS,ND1+LBAS)
     &                     +           RR(N, 3)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N, 7)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N,11)*DENO(NA2+IBAS,NB1+JBAS)
     &                     +           RR(N,15)*DENO(NA2+IBAS,NB2+JBAS)
C
            WDIR(NC2+KBAS,ND2+LBAS) = WDIR(NC2+KBAS,ND2+LBAS)
     &                     +           RR(N, 4)*DENO(NA1+IBAS,NB1+JBAS)
     &                     +           RR(N, 8)*DENO(NA1+IBAS,NB2+JBAS)
     &                     +           RR(N,12)*DENO(NA2+IBAS,NB1+JBAS)
     &                     +           RR(N,16)*DENO(NA2+IBAS,NB2+JBAS)
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
            WDIR(NC1+KBAS,ND1+LBAS) = WDIR(NC1+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,13)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 5)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(N, 9)*DENO(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 1)*DENO(NB2+JBAS,NA2+IBAS)
C
            WDIR(NC1+KBAS,ND2+LBAS) = WDIR(NC1+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(N,14)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 6)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,10)*DENO(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 2)*DENO(NB2+JBAS,NA2+IBAS)
C
            WDIR(NC2+KBAS,ND1+LBAS) = WDIR(NC2+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,15)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 7)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,11)*DENO(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 3)*DENO(NB2+JBAS,NA2+IBAS)
C
            WDIR(NC2+KBAS,ND2+LBAS) = WDIR(NC2+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(N,16)*DENO(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 8)*DENO(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,12)*DENO(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 4)*DENO(NB2+JBAS,NA2+IBAS)
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
            WXCH(NA1+IBAS,ND1+LBAS) = WXCH(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 1)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N, 5)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N, 3)*DENO(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N, 7)*DENO(NC2+KBAS,NB2+JBAS)
C
            WXCH(NA1+IBAS,ND2+LBAS) = WXCH(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N, 2)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N, 6)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N, 4)*DENO(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N, 8)*DENO(NC2+KBAS,NB2+JBAS)
C
            WXCH(NA2+IBAS,ND1+LBAS) = WXCH(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N, 9)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N,13)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N,11)*DENO(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N,15)*DENO(NC2+KBAS,NB2+JBAS)
C
            WXCH(NA2+IBAS,ND2+LBAS) = WXCH(NA2+IBAS,ND2+LBAS)
     &                     +           RR(N,10)*DENO(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N,14)*DENO(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N,12)*DENO(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N,16)*DENO(NC2+KBAS,NB2+JBAS)
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
            WXCH(NA1+IBAS,NC1+KBAS) = WXCH(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 4)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N, 8)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(N, 3)*DENO(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 7)*DENO(ND2+LBAS,NB2+JBAS)
C
            WXCH(NA1+IBAS,NC2+KBAS) = WXCH(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N, 2)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 6)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 1)*DENO(ND2+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N, 5)*DENO(ND2+LBAS,NB2+JBAS)
C
            WXCH(NA2+IBAS,NC1+KBAS) = WXCH(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N,12)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N,16)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(N,11)*DENO(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N,15)*DENO(ND2+LBAS,NB2+JBAS)
C
            WXCH(NA2+IBAS,NC2+KBAS) = WXCH(NA2+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,10)*DENO(ND1+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N,14)*DENO(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 9)*DENO(ND2+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N,13)*DENO(ND2+LBAS,NB2+JBAS)
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
            WXCH(NB1+JBAS,ND1+LBAS) = WXCH(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,13)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 5)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(N,15)*DENO(NC2+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 7)*DENO(NC2+KBAS,NA2+IBAS)
C
            WXCH(NB1+JBAS,ND2+LBAS) = WXCH(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(N,14)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 6)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB1*     RR(N,16)*DENO(NC2+KBAS,NA1+IBAS)
     &                     + PAB2*     RR(N, 8)*DENO(NC2+KBAS,NA2+IBAS)
C
            WXCH(NB2+JBAS,ND1+LBAS) = WXCH(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N, 9)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 1)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,11)*DENO(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 3)*DENO(NC2+KBAS,NA2+IBAS)
C
            WXCH(NB2+JBAS,ND2+LBAS) = WXCH(NB2+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(N,10)*DENO(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 2)*DENO(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N,12)*DENO(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N, 4)*DENO(NC2+KBAS,NA2+IBAS)
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
            WXCH(NB1+JBAS,NC1+KBAS) = WXCH(NB1+JBAS,NC1+KBAS)
     &                     + PAB1*PCD1*RR(N,16)*DENO(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD1*RR(N, 8)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD2*RR(N,15)*DENO(ND2+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(N, 7)*DENO(ND2+LBAS,NA2+IBAS)
C
            WXCH(NB1+JBAS,NC2+KBAS) = WXCH(NB1+JBAS,NC2+KBAS)
     &                     + PAB1*PCD2*RR(N,14)*DENO(ND1+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(N, 6)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB1*PCD1*RR(N,13)*DENO(ND2+LBAS,NA1+IBAS)
     &                     + PAB2*PCD1*RR(N, 5)*DENO(ND2+LBAS,NA2+IBAS)
C
            WXCH(NB2+JBAS,NC1+KBAS) = WXCH(NB2+JBAS,NC1+KBAS)
     &                     + PAB2*PCD1*RR(N,12)*DENO(ND1+LBAS,NA1+IBAS)
     &                     + PAB1*PCD1*RR(N, 4)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD2*RR(N,11)*DENO(ND2+LBAS,NA1+IBAS)
     &                     + PAB1*PCD2*RR(N, 3)*DENO(ND2+LBAS,NA2+IBAS)
C
            WXCH(NB2+JBAS,NC2+KBAS) = WXCH(NB2+JBAS,NC2+KBAS)
     &                     + PAB1*PCD1*RR(N, 1)*DENO(ND2+LBAS,NA2+IBAS)
     &                     + PAB1*PCD2*RR(N, 2)*DENO(ND1+LBAS,NA2+IBAS)
     &                     + PAB2*PCD1*RR(N, 9)*DENO(ND2+LBAS,NA1+IBAS)
     &                     + PAB2*PCD2*RR(N,10)*DENO(ND1+LBAS,NA1+IBAS)
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
            IF(IMTX(M, 9).EQ.0) GOTO 409
C
            WXCH(NC1+KBAS,NB1+JBAS) = WXCH(NC1+KBAS,NB1+JBAS)
     &                     +           RR(N, 1)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 2)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N, 9)*DENO(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N,10)*DENO(NA2+IBAS,ND2+LBAS)
C
            WXCH(NC1+KBAS,NB2+JBAS) = WXCH(NC1+KBAS,NB2+JBAS)
     &                     +           RR(N, 5)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 6)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N,13)*DENO(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N,14)*DENO(NA2+IBAS,ND2+LBAS)
C
            WXCH(NC2+KBAS,NB1+JBAS) = WXCH(NC2+KBAS,NB1+JBAS)
     &                     +           RR(N, 3)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 4)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N,11)*DENO(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N,12)*DENO(NA2+IBAS,ND2+LBAS)
C
            WXCH(NC2+KBAS,NB2+JBAS) = WXCH(NC2+KBAS,NB2+JBAS)
     &                     +           RR(N, 7)*DENO(NA1+IBAS,ND1+LBAS)
     &                     +           RR(N, 8)*DENO(NA1+IBAS,ND2+LBAS)
     &                     +           RR(N,15)*DENO(NA2+IBAS,ND1+LBAS)
     &                     +           RR(N,16)*DENO(NA2+IBAS,ND2+LBAS)
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
            WXCH(NC1+KBAS,NA1+IBAS) = WXCH(NC1+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N,13)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,14)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(N, 9)*DENO(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N,10)*DENO(NB2+JBAS,ND2+LBAS)
C
            WXCH(NC1+KBAS,NA2+IBAS) = WXCH(NC1+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N, 5)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N, 6)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(N, 1)*DENO(NB2+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N, 2)*DENO(NB2+JBAS,ND2+LBAS)
C
            WXCH(NC2+KBAS,NA1+IBAS) = WXCH(NC2+KBAS,NA1+IBAS)
     &                     + PAB1*     RR(N,15)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N,16)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB2*     RR(N,11)*DENO(NB2+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N,12)*DENO(NB2+JBAS,ND2+LBAS)
C
            WXCH(NC2+KBAS,NA2+IBAS) = WXCH(NC2+KBAS,NA2+IBAS)
     &                     + PAB2*     RR(N, 7)*DENO(NB1+JBAS,ND1+LBAS)
     &                     + PAB2*     RR(N, 8)*DENO(NB1+JBAS,ND2+LBAS)
     &                     + PAB1*     RR(N, 3)*DENO(NB2+JBAS,ND1+LBAS)
     &                     + PAB1*     RR(N, 4)*DENO(NB2+JBAS,ND2+LBAS)
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
            WXCH(ND1+LBAS,NB1+JBAS) = WXCH(ND1+LBAS,NB1+JBAS)
     &                     +      PCD1*RR(N, 4)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(N, 2)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(N,12)*DENO(NA2+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(N,10)*DENO(NA2+IBAS,NC2+KBAS)
C
            WXCH(ND1+LBAS,NB2+JBAS) = WXCH(ND1+LBAS,NB2+JBAS)
     &                     +      PCD1*RR(N, 8)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(N, 6)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD1*RR(N,16)*DENO(NA2+IBAS,NC1+KBAS)
     &                     +      PCD2*RR(N,14)*DENO(NA2+IBAS,NC2+KBAS)
C
            WXCH(ND2+LBAS,NB1+JBAS) = WXCH(ND2+LBAS,NB1+JBAS)
     &                     +      PCD2*RR(N, 3)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 1)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,11)*DENO(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 9)*DENO(NA2+IBAS,NC2+KBAS)
C
            WXCH(ND2+LBAS,NB2+JBAS) = WXCH(ND2+LBAS,NB2+JBAS)
     &                     +      PCD2*RR(N, 7)*DENO(NA1+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N, 5)*DENO(NA1+IBAS,NC2+KBAS)
     &                     +      PCD2*RR(N,15)*DENO(NA2+IBAS,NC1+KBAS)
     &                     +      PCD1*RR(N,13)*DENO(NA2+IBAS,NC2+KBAS)
C
411         CONTINUE
          ENDDO
        ENDDO
      ENDIF
C
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

