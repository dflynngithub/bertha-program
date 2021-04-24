C
C     GENERATE BATCH OF ELECTRON REPULSION INTEGRALS
      CALL ERI(RR,XYZ,KQN,MQN,EXPT,NBAS,ITN,IBAS,JBAS)
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
C     2ND CASE (DIRECT):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(2).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
C
            GDIR(NA1+IBAS,NB1+JBAS) = GDIR(NA1+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 2)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 3)*DENC(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENC(ND2+LBAS,NC2+KBAS)
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
     &                     +      PCD1*RR(M,16)*DENC(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,14)*DENC(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENC(ND2+LBAS,NC1+KBAS)
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
     &                     + PAB2*     RR(M, 5)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 9)*DENC(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENC(NB2+JBAS,NA2+IBAS)
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
     &                     + PAB1*     RR(M,16)*DENC(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 8)*DENC(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENC(NB2+JBAS,NA1+IBAS)
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
C     6TH CASE (EXCHNG):  ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      IF(IFLG(6).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
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
C     7TH CASE (EXCHNG):  ( MB, MA| MC, MD) = PAB*    (-MA,-MB| MC, MD)
      IF(IFLG(7).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
            M = M+1
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
C     9TH CASE (EXCHNG):  ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      IF(IFLG(9).EQ.1) THEN
        M = 0
        DO KBAS=1,NBAS(3)
          DO LBAS=1,NBAS(4)
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
