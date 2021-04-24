      SUBROUTINE BREIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***********************************************************************C
C                                                                       C
C            BBBBBBB   RRRRRRR   EEEEEEE  IIII  TTTTTTTT                C
C            BB    BB  RR    RR  EE        II      TT                   C
C            BB    BB  RR    RR  EE        II      TT                   C
C            BBBBBBB   RRRRRRR   EEEEEE    II      TT                   C
C            BB    BB  RR   RR   EE        II      TT                   C
C            BB    BB  RR    RR  EE        II      TT                   C
C            BBBBBBB   RR    RR  EEEEEEE  IIII     TT                   C
C                                                                       C
C                                                                       C
C      BREIT GENERATES ELECTRON REPULSION INTEGRALS IN BATCHES, AND     C
C      THEN USES THEM TO CONSTRUCT THE OPEN AND CLOSED-SHELL FOCK       C
C      MATRICES, EXPLOITING INTEGRAL SYMMETRY                           C
C***********************************************************************C
      INCLUDE 'param.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
      DIMENSION INDEX(MCT,-(MEL+1):MEL,MKP)
C
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 BMAT(MDM,MDM)
      COMPLEX*16 RR(MB2,16)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/DENS/DENC,DENO,DENT
      COMMON/MTRX/BMAT
C
C-----------------------------------------------------------------------C
C     INDEXING ROUTINE
C-----------------------------------------------------------------------C
      ICOUNT=0
      IKOUNT=0
C
      DO ICT=1,NCNT
        DO KN=1,NKAP(ICT)
          KAPPA=KAPA(KN,ICT)
          MJMAX=2*IABS(KAPPA)-1
          IF(KAPPA.GT.0) THEN
            LORB=KAPPA
          ELSE
            LORB=-KAPPA-1
          ENDIF
          NFUN=NFNC(LORB,ICT)
          IKOUNT=IKOUNT+NFUN
          DO MJ=1,MJMAX,2
            ICOUNT=ICOUNT+1
            INDEX(ICT,KAPPA,MJ)=ICOUNT
            LENIQ=ICOUNT
          ENDDO
C
        ENDDO
      ENDDO
C-----------------------------------------------------------------------C
C     CLEAR OUT BMAT MATRIX
C
      DO J=1,NDIM
        DO I=1,NDIM
          BMAT(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C-----------------------------------------------------------------------C
C
      DO 2000 ICNTA=1,NCNT
      XYZ(1,1)=BXYZ(1,ICNTA)
      XYZ(2,1)=BXYZ(2,ICNTA)
      XYZ(3,1)=BXYZ(3,ICNTA)

      DO 2000 ICNTB=1,NCNT
      XYZ(1,2)=BXYZ(1,ICNTB)
      XYZ(2,2)=BXYZ(2,ICNTB)
      XYZ(3,2)=BXYZ(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
C
      DO 2000 KA=1,NKAP(ICNTA)
      KQN(1)=KAPA(KA,ICNTA)
      IF(KQN(1).GT.0) THEN
        LQN(1)=KQN(1)
      ELSE
        LQN(1)=-KQN(1)-1
      ENDIF
C
      NBAS(1)=NFNC(LQN(1),ICNTA)
      DO IBAS=1,NBAS(1)
        EXL(IBAS,1)=BEXL(IBAS,LQN(1),ICNTA)
      ENDDO
C
C     LOOP OVER KQN(B) VALUES
C
      DO 2000 KB=1,NKAP(ICNTB)
      KQN(2)=KAPA(KB,ICNTB)
      IF(KQN(2).GT.0) THEN
        LQN(2)=KQN(2)
      ELSE
        LQN(2)=-KQN(2)-1
      ENDIF
C
      NBAS(2)=NFNC(LQN(2),ICNTB)
      DO JBAS=1,NBAS(2)
        EXL(JBAS,2)=BEXL(JBAS,LQN(2),ICNTB)
      ENDDO
C
C     LOOP OVER |MA| VALUES
C
      DO 2000 MA=1,IABS(KQN(1))
        MQN(1) = 2*MA-1
C
C     LOOP OVER |MB| VALUES
C
      DO 2000 MB=1,IABS(KQN(2))
        MQN(2) = 2*MB-1
C
C     LOOP OVER CENTRES C AND D
C
      DO 1000 ICNTC=1,NCNT
      XYZ(1,3)=BXYZ(1,ICNTC)
      XYZ(2,3)=BXYZ(2,ICNTC)
      XYZ(3,3)=BXYZ(3,ICNTC)
C
      DO 1000 ICNTD=1,NCNT
      XYZ(1,4)=BXYZ(1,ICNTD)
      XYZ(2,4)=BXYZ(2,ICNTD)
      XYZ(3,4)=BXYZ(3,ICNTD)
C
C     LOOP OVER KQN(C) VALUES
C
      DO 1000 KC=1,NKAP(ICNTC)
      KQN(3)=KAPA(KC,ICNTC)
      IF(KQN(3).GT.0) THEN
        LQN(3)=KQN(3)
      ELSE
        LQN(3)=-KQN(3)-1
      ENDIF
C
      NBAS(3)=NFNC(LQN(3),ICNTC)
      DO KBAS=1,NBAS(3)
        EXL(KBAS,3)=BEXL(KBAS,LQN(3),ICNTC)
      ENDDO
C
C     LOOP OVER KQN(D) VALUES
C
      DO 1000 KD=1,NKAP(ICNTD)
      KQN(4)=KAPA(KD,ICNTD)
      IF(KQN(4).GT.0) THEN
        LQN(4)=KQN(4)
      ELSE
        LQN(4)=-KQN(4)-1
      ENDIF
C
      NBAS(4)=NFNC(LQN(4),ICNTD)
      DO LBAS=1,NBAS(4)
        EXL(LBAS,4)=BEXL(LBAS,LQN(4),ICNTD)
      ENDDO
C
C     LOOP OVER |MC| AND |MD| VALUES
C
      DO 1000 MC=1,IABS(KQN(3))
        MQN(3) = 2*MC-1
C
      DO 1000 MD=1,IABS(KQN(4))
        MQN(4) = 2*MD-1
C
      MAXCD = NBAS(3)*NBAS(4)  
C
C**********************************************************************C
C
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
      IQ12 = (IQ1-1)*LENIQ+IQ2
      IQ34 = (IQ3-1)*LENIQ+IQ4
C
      IF(IQ12.LT.IQ34) GOTO 1999
C
C     CALCULATE PHASE FACTORS FOR PERMUTING INTEGRALS.
C     NB! OPPOSITE PHASE AS IN THE LL/SS CASE SEEN IN SCF
      PAB1 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PAB2 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      PCD1 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PCD2 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
      NA1 = LRGE(ICNTA,KA,2*MA-1)
      NA2 = LRGE(ICNTA,KA,2*MA  )
      NB1 = LRGE(ICNTB,KB,2*MB-1)+NSKP
      NB2 = LRGE(ICNTB,KB,2*MB  )+NSKP
      NC1 = LRGE(ICNTC,KC,2*MC-1)
      NC2 = LRGE(ICNTC,KC,2*MC  )
      ND1 = LRGE(ICNTD,KD,2*MD-1)+NSKP
      ND2 = LRGE(ICNTD,KD,2*MD  )+NSKP
C
      DO 3000 IBAS=1,NBAS(1)
      DO 3000 JBAS=1,NBAS(2)
C
C*********************************************************************C
C
      CALL BRTINT(RR,XYZ,KQN,MQN,EXL,NBAS,IBAS,JBAS)
C
C     BATCH 01: B^{LS,LS}
C     DIRECT INTEGRALS    ( MA, MB| MD, MC) =         ( MA, MB| MD, MC)
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
C
          BMAT(NA1+IBAS,NB1+JBAS) = BMAT(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 1)*DENT(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 2)*DENT(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 3)*DENT(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 4)*DENT(NC2+KBAS,ND2+LBAS)
C
          BMAT(NA1+IBAS,NB2+JBAS) = BMAT(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 5)*DENT(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 6)*DENT(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 7)*DENT(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 8)*DENT(NC2+KBAS,ND2+LBAS)
C
          BMAT(NA2+IBAS,NB1+JBAS) = BMAT(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M, 9)*DENT(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,10)*DENT(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,11)*DENT(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M,12)*DENT(NC2+KBAS,ND2+LBAS)
C
          BMAT(NA2+IBAS,NB2+JBAS) = BMAT(NA2+IBAS,NB2+JBAS)
     &                     +           RR(M,13)*DENT(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M,14)*DENT(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M,15)*DENT(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M,16)*DENT(NC2+KBAS,ND2+LBAS)
C
        ENDDO
      ENDDO
C
C     BATCH 02: B^{LS,SL}
C     DIRECT INTEGRALS    ( MA, MB| MD, MC) =     PCD*( MA, MB|-MC,-MD)
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
C
          BMAT(NA1+IBAS,NB1+JBAS) = BMAT(NA1+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENT(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 2)*DENT(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 3)*DENT(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENT(ND2+LBAS,NC2+KBAS)
C
          BMAT(NA1+IBAS,NB2+JBAS) = BMAT(NA1+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(M, 8)*DENT(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENT(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M, 7)*DENT(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 5)*DENT(ND2+LBAS,NC2+KBAS)
C
          BMAT(NA2+IBAS,NB1+JBAS) = BMAT(NA2+IBAS,NB1+JBAS)
     &                     +      PCD1*RR(M,12)*DENT(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENT(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENT(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENT(ND2+LBAS,NC2+KBAS)
C
          BMAT(NA2+IBAS,NB2+JBAS) = BMAT(NA2+IBAS,NB2+JBAS)
     &                     +      PCD1*RR(M,16)*DENT(ND1+LBAS,NC1+KBAS)
     &                     +      PCD2*RR(M,14)*DENT(ND1+LBAS,NC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENT(ND2+LBAS,NC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENT(ND2+LBAS,NC2+KBAS)
C
        ENDDO
      ENDDO
C
C     BATCH 03: B^{LS,SL}        ~       ~
C     EXCHANGE INTEGRALS  ( MA, MD| MC, MB) =         ( MA, MB| MC, MD)
      DO LBAS=1,NBAS(4)
        DO KBAS=1,NBAS(3)
          M=(K-1)*NBAS(4)+LBAS
C
          BMAT(NA1+IBAS,ND1+LBAS) = BMAT(NA1+IBAS,ND1+LBAS)
     &                     -           RR(M, 1)*DENT(NC1+KBAS,NB1+JBAS)
     &                     -           RR(M, 3)*DENT(NC2+KBAS,NB1+JBAS)
     &                     -           RR(M, 5)*DENT(NC1+KBAS,NB2+JBAS)
     &                     -           RR(M, 7)*DENT(NC2+KBAS,NB2+JBAS)
C
          BMAT(NA1+IBAS,ND2+LBAS) = BMAT(NA1+IBAS,ND2+LBAS)
     &                     -           RR(M, 2)*DENT(NC1+KBAS,NB1+JBAS)
     &                     -           RR(M, 4)*DENT(NC2+KBAS,NB1+JBAS)
     &                     -           RR(M, 6)*DENT(NC1+KBAS,NB2+JBAS)
     &                     -           RR(M, 8)*DENT(NC2+KBAS,NB2+JBAS)
C
          BMAT(NA2+IBAS,ND1+LBAS) = BMAT(NA2+IBAS,ND1+LBAS)
     &                     -           RR(M, 9)*DENT(NC1+KBAS,NB1+JBAS)
     &                     -           RR(M,11)*DENT(NC2+KBAS,NB1+JBAS)
     &                     -           RR(M,13)*DENT(NC1+KBAS,NB2+JBAS)
     &                     -           RR(M,15)*DENT(NC2+KBAS,NB2+JBAS)
C
          BMAT(NA2+IBAS,ND2+LBAS) = BMAT(NA2+IBAS,ND2+LBAS)
     &                     -           RR(M,10)*DENT(NC1+KBAS,NB1+JBAS)
     &                     -           RR(M,12)*DENT(NC2+KBAS,NB1+JBAS)
     &                     -           RR(M,14)*DENT(NC1+KBAS,NB2+JBAS)
     &                     -           RR(M,16)*DENT(NC2+KBAS,NB2+JBAS)
C
        ENDDO
      ENDDO
C
C     BATCH 04: B^{LL,SS}        ~       ~
C     EXCHANGE INTEGRALS  ( MA, MC| MD, MB) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
C
          BMAT(NA1+IBAS,NC1+K) = BMAT(NA1+IBAS,NC1+K)
     &                     -      PCD1*RR(M, 4)*DENT(ND1+LBAS,NB1+JBAS)
     &                     -      PCD1*RR(M, 8)*DENT(ND1+LBAS,NB2+JBAS)
     &                     -      PCD1*RR(M, 3)*DENT(ND2+LBAS,NB1+JBAS)
     &                     -      PCD1*RR(M, 7)*DENT(ND2+LBAS,NB2+JBAS)
C
          BMAT(NA1+IBAS,NC2+K) = BMAT(NA1+IBAS,NC2+K)
     &                     -      PCD1*RR(M, 2)*DENT(ND1+LBAS,NB1+JBAS)
     &                     -      PCD1*RR(M, 6)*DENT(ND1+LBAS,NB2+JBAS)
     &                     -      PCD1*RR(M, 1)*DENT(ND2+LBAS,NB1+JBAS)
     &                     -      PCD1*RR(M, 5)*DENT(ND2+LBAS,NB2+JBAS)
C
          BMAT(NA2+IBAS,NC1+K) = BMAT(NA2+IBAS,NC1+K)
     &                     -      PCD1*RR(M,12)*DENT(ND1+LBAS,NB1+JBAS)
     &                     -      PCD1*RR(M,16)*DENT(ND1+LBAS,NB2+JBAS)
     &                     -      PCD1*RR(M,11)*DENT(ND2+LBAS,NB1+JBAS)
     &                     -      PCD1*RR(M,15)*DENT(ND2+LBAS,NB2+JBAS)
C
          BMAT(NA2+IBAS,NC2+K) = BMAT(NA2+IBAS,NC2+K)
     &                     -      PCD1*RR(M,10)*DENT(ND1+LBAS,NB1+JBAS)
     &                     -      PCD1*RR(M,14)*DENT(ND1+LBAS,NB2+JBAS)
     &                     -      PCD1*RR(M, 9)*DENT(ND2+LBAS,NB1+JBAS)
     &                     -      PCD1*RR(M,13)*DENT(ND2+LBAS,NB2+JBAS)
C
        ENDDO
      ENDDO
C
C     BATCH 05: B^{SS,LL}        ~       ~
C     EXCHANGE INTEGRALS  ( MB, MD| MC, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      DO LBAS=1,NBAS(4)
        DO KBAS=1,NBAS(3)
          M=(K-1)*NBAS(4)+LBAS
C
          BMAT(NB1+JBAS,ND1+LBAS) = BMAT(NB1+JBAS,ND1+LBAS)
     &                     - PAB1*     RR(M,13)*DENT(NC1+KBAS,NA1+IBAS)
     &                     - PAB2*     RR(M, 5)*DENT(NC1+KBAS,NA2+IBAS)
     &                     - PAB1*     RR(M,15)*DENT(NC2+KBAS,NA1+IBAS)
     &                     - PAB2*     RR(M, 7)*DENT(NC2+KBAS,NA2+IBAS)
C
          BMAT(NB1+JBAS,ND2+LBAS) = BMAT(NB1+JBAS,ND2+LBAS)
     &                     - PAB1*     RR(M,14)*DENT(NC1+KBAS,NA1+IBAS)
     &                     - PAB2*     RR(M, 6)*DENT(NC1+KBAS,NA2+IBAS)
     &                     - PAB1*     RR(M,16)*DENT(NC2+KBAS,NA1+IBAS)
     &                     - PAB2*     RR(M, 8)*DENT(NC2+KBAS,NA2+IBAS)
C
          BMAT(NB2+JBAS,ND1+LBAS) = BMAT(NB2+JBAS,ND1+LBAS)
     &                     - PAB2*     RR(M, 9)*DENT(NC1+KBAS,NA1+IBAS)
     &                     - PAB1*     RR(M, 1)*DENT(NC1+KBAS,NA2+IBAS)
     &                     - PAB2*     RR(M,11)*DENT(NC2+KBAS,NA1+IBAS)
     &                     - PAB1*     RR(M, 3)*DENT(NC2+KBAS,NA2+IBAS)
C
          BMAT(NB2+JBAS,ND2+LBAS) = BMAT(NB2+JBAS,ND2+LBAS)
     &                     - PAB2*     RR(M,10)*DENT(NC1+KBAS,NA1+IBAS)
     &                     - PAB1*     RR(M, 2)*DENT(NC1+KBAS,NA2+IBAS)
     &                     - PAB2*     RR(M,12)*DENT(NC2+KBAS,NA1+IBAS)
     &                     - PAB1*     RR(M, 4)*DENT(NC2+KBAS,NA2+IBAS)
C
        ENDDO
      ENDDO
C
C-------------------------------------------------------------
C
      IF(IQ12.EQ.IQ34) GOTO 222
C
C     BATCH 06: B^{LS,LS}
C     DIRECT INTEGRALS    ( MC, MD| MA, MB) =         ( MA, MB| MC, MD)
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
C
          BMAT(NC1+KBAS,ND1+LBAS) = BMAT(NC1+KBAS,ND1+LBAS)
     &                     +           RR(M, 1)*DENT(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 5)*DENT(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M, 9)*DENT(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,13)*DENT(NA2+IBAS,NB2+JBAS)
C
          BMAT(NC1+KBAS,ND2+LBAS) = BMAT(NC1+KBAS,ND2+LBAS)
     &                     +           RR(M, 2)*DENT(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 6)*DENT(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,10)*DENT(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,14)*DENT(NA2+IBAS,NB2+JBAS)      
C
          BMAT(NC2+KBAS,ND1+LBAS) = BMAT(NC2+KBAS,ND1+LBAS)
     &                     +           RR(M, 3)*DENT(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 7)*DENT(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,11)*DENT(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,15)*DENT(NA2+IBAS,NB2+JBAS)
C
          BMAT(NC2+KBAS,ND2+LBAS) = BMAT(NC2+KBAS,ND2+LBAS)
     &                     +           RR(M, 4)*DENT(NA1+IBAS,NB1+JBAS)
     &                     +           RR(M, 8)*DENT(NA1+IBAS,NB2+JBAS)
     &                     +           RR(M,12)*DENT(NA2+IBAS,NB1+JBAS)
     &                     +           RR(M,16)*DENT(NA2+IBAS,NB2+JBAS)
C
        ENDDO
      ENDDO
C
C     BATCH 07: B^{LS,LS}
C     DIRECT INTEGRALS    ( MC, MD| MB, MA) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
C
          BMAT(NC1+KBAS,ND1+LBAS) = BMAT(NC1+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,13)*DENT(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 5)*DENT(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M, 9)*DENT(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENT(NB2+JBAS,NA2+IBAS)
C
          BMAT(NC1+KBAS,ND2+LBAS) = BMAT(NC1+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,14)*DENT(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENT(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,10)*DENT(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 2)*DENT(NB2+JBAS,NA2+IBAS)
C
          BMAT(NC2+KBAS,ND1+LBAS) = BMAT(NC2+KBAS,ND1+LBAS)
     &                     + PAB1*     RR(M,15)*DENT(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENT(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENT(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENT(NB2+JBAS,NA2+IBAS)
C
          BMAT(NC2+KBAS,ND2+LBAS) = BMAT(NC2+KBAS,ND2+LBAS)
     &                     + PAB1*     RR(M,16)*DENT(NB1+JBAS,NA1+IBAS)
     &                     + PAB2*     RR(M, 8)*DENT(NB1+JBAS,NA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENT(NB2+JBAS,NA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENT(NB2+JBAS,NA2+IBAS)
C
        ENDDO
      ENDDO
C
C     BATCH 08: B^{LL,SS}    ~       ~    
C     EXCHANGE INTEGRALS  ( MC, MA| MB, MD) =         ( MB, MA| MC, MD)
C                                           = PAB*    (-MA,-MB| MC, MD)
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
C
          BMAT(NC1+KBAS,NA1+IBAS) = BMAT(NC1+KBAS,NA1+IBAS)
     &                     - PAB1*     RR(M,13)*DENT(NB1+JBAS,ND1+LBAS)
     &                     - PAB1*     RR(M,14)*DENT(NB1+JBAS,ND2+LBAS)
     &                     - PAB2*     RR(M, 9)*DENT(NB2+JBAS,ND1+LBAS)
     &                     - PAB2*     RR(M,10)*DENT(NB2+JBAS,ND2+LBAS)
C
          BMAT(NC1+KBAS,NA2+IBAS) = BMAT(NC1+KBAS,NA2+IBAS)
     &                     - PAB2*     RR(M, 5)*DENT(NB1+JBAS,ND1+LBAS)
     &                     - PAB2*     RR(M, 6)*DENT(NB1+JBAS,ND2+LBAS)
     &                     - PAB1*     RR(M, 1)*DENT(NB2+JBAS,ND1+LBAS)
     &                     - PAB1*     RR(M, 2)*DENT(NB2+JBAS,ND2+LBAS)
C
          BMAT(NC2+KBAS,NA1+IBAS) = BMAT(NC2+KBAS,NA1+IBAS)
     &                     - PAB1*     RR(M,15)*DENT(NB1+JBAS,ND1+LBAS)
     &                     - PAB1*     RR(M,16)*DENT(NB1+JBAS,ND2+LBAS)
     &                     - PAB2*     RR(M,11)*DENT(NB2+JBAS,ND1+LBAS)
     &                     - PAB2*     RR(M,12)*DENT(NB2+JBAS,ND2+LBAS)
C
          BMAT(NC2+KBAS,NA2+IBAS) = BMAT(NC2+KBAS,NA2+IBAS)
     &                     - PAB2*     RR(M, 7)*DENT(NB1+JBAS,ND1+LBAS)
     &                     - PAB2*     RR(M, 8)*DENT(NB1+JBAS,ND2+LBAS)
     &                     - PAB1*     RR(M, 3)*DENT(NB2+JBAS,ND1+LBAS)
     &                     - PAB1*     RR(M, 4)*DENT(NB2+JBAS,ND2+LBAS)
C
        ENDDO
      ENDDO
C
C     BATCH 09: B^{SS,LL}    ~       ~    
C     EXCHANGE INTEGRALS  ( MD, MB| MA, MC) =         ( MA, MB| MD, MC)
C                                           =     PCD*( MA, MB|-MC,-MD)
      DO LBAS=1,NBAS(4)
        DO KBAS=1,NBAS(3)
          M=(K-1)*NBAS(4)+LBAS
C
          BMAT(ND1+LBAS,NB1+JBAS) = BMAT(ND1+LBAS,NB1+JBAS)
     &                     -      PCD1*RR(M, 4)*DENT(NA1+IBAS,NC1+KBAS)
     &                     -      PCD2*RR(M, 2)*DENT(NA1+IBAS,NC2+KBAS)
     &                     -      PCD1*RR(M,12)*DENT(NA2+IBAS,NC1+KBAS)
     &                     -      PCD2*RR(M,10)*DENT(NA2+IBAS,NC2+KBAS)
C
          BMAT(ND1+LBAS,NB2+JBAS) = BMAT(ND1+LBAS,NB2+JBAS)
     &                     -      PCD1*RR(M, 8)*DENT(NA1+IBAS,NC1+KBAS)
     &                     -      PCD2*RR(M, 6)*DENT(NA1+IBAS,NC2+KBAS)
     &                     -      PCD1*RR(M,16)*DENT(NA2+IBAS,NC1+KBAS)
     &                     -      PCD2*RR(M,14)*DENT(NA2+IBAS,NC2+KBAS)
C
          BMAT(ND2+LBAS,NB1+JBAS) = BMAT(ND2+LBAS,NB1+JBAS)
     &                     -      PCD2*RR(M, 3)*DENT(NA1+IBAS,NC1+KBAS)
     &                     -      PCD1*RR(M, 1)*DENT(NA1+IBAS,NC2+KBAS)
     &                     -      PCD2*RR(M,11)*DENT(NA2+IBAS,NC1+KBAS)
     &                     -      PCD1*RR(M, 9)*DENT(NA2+IBAS,NC2+KBAS)
C
          BMAT(ND2+LBAS,NB2+JBAS) = BMAT(ND2+LBAS,NB2+JBAS)
     &                     -      PCD2*RR(M, 7)*DENT(NA1+IBAS,NC1+KBAS)
     &                     -      PCD1*RR(M, 5)*DENT(NA1+IBAS,NC2+KBAS)
     &                     -      PCD2*RR(M,15)*DENT(NA2+IBAS,NC1+KBAS)
     &                     -      PCD1*RR(M,13)*DENT(NA2+IBAS,NC2+KBAS)
C
        ENDDO
      ENDDO
C
C     BATCH 10: B^{SS,LL}    ~       ~    
C     EXCHANGE INTEGRALS  ( MC, MB| MA, MD) =         ( MA, MB| MD, MC)
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
C
          BMAT(NC1+KBAS,NB1+JBAS) = BMAT(NC1+KBAS,NB1+JBAS)
     &                     -           RR(M, 1)*DENT(NA1+IBAS,ND1+LBAS)
     &                     -           RR(M, 2)*DENT(NA1+IBAS,ND2+LBAS)
     &                     -           RR(M, 9)*DENT(NA2+IBAS,ND1+LBAS)
     &                     -           RR(M,10)*DENT(NA2+IBAS,ND2+LBAS)
C
          BMAT(NC1+KBAS,NB2+JBAS) = BMAT(NC1+KBAS,NB2+JBAS)
     &                     -           RR(M, 5)*DENT(NA1+IBAS,ND1+LBAS)
     &                     -           RR(M, 6)*DENT(NA1+IBAS,ND2+LBAS)
     &                     -           RR(M,13)*DENT(NA2+IBAS,ND1+LBAS)
     &                     -           RR(M,14)*DENT(NA2+IBAS,ND2+LBAS)
C
          BMAT(NC2+KBAS,NB1+JBAS) = BMAT(NC2+KBAS,NB1+JBAS)
     &                     -           RR(M, 3)*DENT(NA1+IBAS,ND1+LBAS)
     &                     -           RR(M, 4)*DENT(NA1+IBAS,ND2+LBAS)
     &                     -           RR(M,11)*DENT(NA2+IBAS,ND1+LBAS)
     &                     -           RR(M,12)*DENT(NA2+IBAS,ND2+LBAS)
C
          BMAT(NC2+KBAS,NB2+JBAS) = BMAT(NC2+KBAS,NB2+JBAS)
     &                     -           RR(M, 7)*DENT(NA1+IBAS,ND1+LBAS)
     &                     -           RR(M, 8)*DENT(NA1+IBAS,ND2+LBAS)
     &                     -           RR(M,15)*DENT(NA2+IBAS,ND1+LBAS)
     &                     -           RR(M,16)*DENT(NA2+IBAS,ND2+LBAS)
C
        ENDDO
      ENDDO
C
222   CONTINUE
C
3000  CONTINUE
4001  CONTINUE
4000  CONTINUE
1999  CONTINUE
1000  CONTINUE
2000  CONTINUE
C
C     COMPLETE BREIT MATRIX BY HERMITICITY
      DO I=1,NSKP
        DO J=NSKP+1,NDIM
          BMAT(J,I) = DCONJG(BMAT(I,J))
        ENDDO
      ENDDO
C
C     BREIT INTERACTION ENERGY
      BENRGY=0.0D0
      DO I=1,NDIM
        DO J=1,NDIM
          BENRGY = BENRGY + 0.5D0*DREAL(DENT(I,J)*BMAT(I,J))
        ENDDO
      ENDDO
C
C      WRITE(6,*) 'BREIT ENERGY = ',BENRGY
C
      RETURN
      END
C
C
      SUBROUTINE BRTINT(RINTG,XYZ,KQN,MQN,EXL,NBAS,IBAS,JBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C****************************************************************C
C               ERI: ELECTRON REPULSION INTEGRALS                C
C                                                                C
C     ERI GENERATES BLOCKS OF ELECTRON REPULSION INTEGRALS       C
C     OVER KINETICALLY BALANCED G-SPINOR BASIS FUNCTIONS         C
C                                                                C
C     THE DENSITIES ARE EXPANDED IN A BASIS OF HERMITE GAUSSIANS C
C     AND THE INTEGRALS ARE GENERATED USING THE MCMURCHIE-       C
C     DAVIDSON ALGORITHM                                         C
C                                                                C
C                                                                C
C                      INPUT PARAMETERS                          C
C                                                                C
C     XYZ(3,4)    COORDINATES OF THE 4 NUCLEAR CENTRES           C
C     KQN(4)      KQN QUANTUM NUMBERS OF THE CENTRES             C
C     MQN(4)      |M|   QUANTUM NUMBERS OF THE CENTRES           C
C     EXL(I,J)    EXPONENTS ON CENTRE J                          C
C     NBAS(J)    NUMBER OF FUNCTIONS ON CENTRE J                 C
C     I,J         INDEX FOR BASIS FUNCTION PAIR ON AB            C
C     IEMAKE      IEMAKE=0 DON'T RECALCULATE E-COEFFICIENTS      C
C                 IEMAKE=1 DO    RECALCULATE E-COEFFICIENTS      C
C****************************************************************C
      INCLUDE 'param.h'
C
      COMPLEX*16 RINTG(MB2,16)
C
      DIMENSION RR(MB2,16),RI(MB2,16)
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
      DIMENSION PQ(MB2,3),APH(MB2),PRE(MB2),RC(MB2,MRC)
C
      DIMENSION GABR11(MB2,MEQ),GABI11(MB2,MEQ),
     &          GABR21(MB2,MEQ),GABI21(MB2,MEQ)
      DIMENSION QR1(MB2),QI1(MB2),QR2(MB2),QI2(MB2) 
      DIMENSION IDX(3),JDX(3)
C
      SAVE LAMAB,LAMCD
C
      COMMON/ESAVE/ERAB11(MB2,MEQ,3),EIAB11(MB2,MEQ,3),
     &             ERAB21(MB2,MEQ,3),EIAB21(MB2,MEQ,3),
     &             ERCD11(MB2,MEQ,3),EICD11(MB2,MEQ,3),
     &             ERCD21(MB2,MEQ,3),EICD21(MB2,MEQ,3)
      COMMON/BCOMP/IABR11(MEQ,3),IABR21(MEQ,3),
     &             IABI11(MEQ,3),IABI21(MEQ,3),
     &             ICDR11(MEQ,3),ICDR21(MEQ,3),
     &             ICDI11(MEQ,3),ICDI21(MEQ,3),IRC(MRC)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
      DO N=1,4
        IF(KQN(N).LT.0) THEN
          LQN(N)=-KQN(N)-1
        ELSE
          LQN(N)=KQN(N)
        ENDIF
      ENDDO
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NBAS(1)*NBAS(2)
      MAXCD = NBAS(3)*NBAS(4)
C
C     PHASE FACTOR FOR AB AND CD PAIR OVERLAPS
      IPHSAB = 1
      IPHSCD =-1
C
C     VRS MAXIMUM LAMBDA LEVEL FOR EQ-COEFFICIENT ADDRESSES
      LAMAB = LQN(1)+LQN(2)+1
      LAMCD = LQN(3)+LQN(4)+1
C
C     VRS MAXIMUM LAMBDA LEVEL FOR CONTRACTED R-INTEGRAL BATCH
      LAMABCD = LAMAB+LAMCD+2
C
C     VRS TOTAL LENGTH OF EQ-COEFFICIENT LISTS AND R-INTEGRAL BATCH
      NTUVAB   = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
      NTUVCD   = (LAMCD+1)*(LAMCD+2)*(LAMCD+3)/6
      NTUVABCD = (LAMABCD+1)*(LAMABCD+2)*(LAMABCD+3)/6
C
C     LIST ADDRESS FOR (AB|  ) AND GAUSSIAN EXPONENT FOR AB OVERLAP
      IJ  = (IBAS-1)*NBAS(2)+JBAS
      EIJ = EXL(IBAS,1)+EXL(JBAS,2)
C
C     INITIALISE RR ARRAY
      DO M=1,MAXCD
        DO ITG=1,16
          RR(M,ITG) = 0.0D0
          RI(M,ITG) = 0.0D0
        ENDDO
      ENDDO
C
C**********************************************************************C
C     GENERATE NEW BATCH OF E(AB|  ) COEFFICIENTS IF PROMPTED          C
C**********************************************************************C
C
      CALL RNORMF(EXL,LQN,NBAS,1,2)
      CALL EMAKELS(ERAB11,EIAB11,ERAB21,EIAB21,EXL,
     &        KQN,MQN,NBAS,IPHSAB,1,2)
C
      DO ICP=1,3
        DO IAB=1,NTUVAB
C
C         REAL(1,1)
          TEST=DASUM(MAXAB,ERAB11(1,IAB,ICP),1)
          IF(TEST.LE.1.0D-14) THEN
            IABR11(IAB,ICP)=0
          ELSE
            IABR11(IAB,ICP)=1
          ENDIF
C
C         IMAGINARY (1,1)
          TEST=DASUM(MAXAB,EIAB11(1,IAB,ICP),1)
          IF(TEST.LE.1.0D-14) THEN
            IABI11(IAB,ICP)=0
          ELSE
            IABI11(IAB,ICP)=1
          ENDIF
C
C         REAL(2,1)
          TEST=DASUM(MAXAB,ERAB21(1,IAB,ICP),1)
          IF(TEST.LE.1.0D-14) THEN
            IABR21(IAB,ICP)=0
          ELSE
            IABR21(IAB,ICP)=1
          ENDIF
C
C         IMAGINARY (2,1)
          TEST=DASUM(MAXAB,EIAB21(1,IAB,ICP),1)
          IF(TEST.LE.1.0D-14) THEN
            IABI21(IAB,ICP)=0
          ELSE
            IABI21(IAB,ICP)=1
          ENDIF
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     GENERATE NEW BATCH OF E(CD| -) COEFFICIENTS IF PROMPTED          C
C**********************************************************************C
C
      CALL RNORMF(EXL,LQN,NBAS,3,4)
      CALL EMAKELS(ERCD11,EICD11,ERCD21,EICD21,EXL,
     &        KQN,MQN,NBAS,IPHSCD,3,4)
C
      DO ICP=1,3
        DO ICD=1,NTUVCD
C
C         REAL(1,1)
          TEST=DASUM(MAXCD,ERCD11(1,ICD,ICP),1)
          IF(TEST.LE.1.0D-14) THEN
            ICDR11(ICD,ICP)=0
          ELSE
            ICDR11(ICD,ICP)=1
          ENDIF
C     
C         IMAGINARY (1,1)
          TEST=DASUM(MAXCD,EICD11(1,ICD,ICP),1)
          IF(TEST.LE.1.0D-14) THEN
            ICDI11(ICD,ICP)=0
          ELSE
            ICDI11(ICD,ICP)=1
          ENDIF
C
C         REAL(2,1)
          TEST=DASUM(MAXCD,ERCD21(1,ICD,ICP),1)
          IF(TEST.LE.1.0D-14) THEN
            ICDR21(ICD,ICP)=0
          ELSE
            ICDR21(ICD,ICP)=1
          ENDIF
C     
C         IMAGINARY (2,1)
          TEST=DASUM(MAXCD,EICD21(1,ICD,ICP),1)
          IF(TEST.LE.1.0D-14) THEN
            ICDI21(ICD,ICP)=0
          ELSE
            ICDI21(ICD,ICP)=1
          ENDIF
C
        ENDDO
      ENDDO
C
C****************************************************************C
C     EVALUATE GEOMETRIC PARAMETERS FOR THE R-INTEGRALS          C
C****************************************************************C
C
      PX = (XYZ(1,1)*EXL(IBAS,1)+XYZ(1,2)*EXL(JBAS,2))/EIJ
      PY = (XYZ(2,1)*EXL(IBAS,1)+XYZ(2,2)*EXL(JBAS,2))/EIJ
      PZ = (XYZ(3,1)*EXL(IBAS,1)+XYZ(3,2)*EXL(JBAS,2))/EIJ

      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
C
          EKL = EXL(K,3)+EXL(L,4)
          QX  = (XYZ(1,3)*EXL(K,3)+XYZ(1,4)*EXL(L,4))/EKL
          QY  = (XYZ(2,3)*EXL(K,3)+XYZ(2,4)*EXL(L,4))/EKL
          QZ  = (XYZ(3,3)*EXL(K,3)+XYZ(3,4)*EXL(L,4))/EKL
C
          APH(M)  = (EIJ*EKL)/(EIJ+EKL)
          PQ(M,1) = QX-PX
          PQ(M,2) = QY-PY
          PQ(M,3) = QZ-PZ
          EMX     = DSQRT(EIJ+EKL)*EIJ*EKL
          PRE(M)  = 2.0D0*PI52/EMX
C
        ENDDO
      ENDDO
C
      CALL RMAKE(RC,PQ,APH,MAXCD,LAMABCD)
C
C     INITIALIZE ARRAY TO IMPLEMENT SPARSENESS IN R-VECTOR
      DO MRC=1,NTUVABCD
        TEST=DASUM(MAXCD,RC(1,MRC),1)
        IF(TEST.LE.1.0D-14) THEN
          IRC(MRC)=0
        ELSE
          IRC(MRC)=1
        ENDIF
      ENDDO
C
C****************************************************************C
C     CONSTRUCT INTERMEDIATE MATRICES FOR MCMURCHIE-DAVIDSON     C
C****************************************************************C
C
      DO 6000 ICMP=1,3
C
      DO IAB=1,NTUVAB
C
C----------------------------------------------------------------C
C
        IAB1=0
        IAB2=0
        IAB3=0
        IAB4=0
C
        IF((IABR11(IAB,1)+IABR11(IAB,2)+IABR11(IAB,3)).NE.0) THEN
        IAB1=1
        ENDIF
C
        IF((IABI11(IAB,1)+IABI11(IAB,2)+IABI11(IAB,3)).NE.0) THEN
        IAB2=1
        ENDIF
C
        IF((IABR21(IAB,1)+IABR21(IAB,2)+IABR21(IAB,3)).NE.0) THEN
        IAB3=1
        ENDIF
C
        IF((IABI21(IAB,1)+IABI21(IAB,2)+IABI21(IAB,3)).NE.0) THEN
        IAB4=1
        ENDIF
C----------------------------------------------------------------C
C
C
        IF((IAB1+IAB2+IAB3+IAB4).GT.0) THEN
          DO M=1,MAXCD
            GABR11(M,IAB)=0.0D0
            GABI11(M,IAB)=0.0D0
            GABR21(M,IAB)=0.0D0
            GABI21(M,IAB)=0.0D0
          ENDDO
        ENDIF
C
C-------------------------------------------------
C
        DO 6100 JCMP=1,3
C
        IF(ICMP.EQ.JCMP) THEN
C
          DO ICD=1,NTUVCD
C
            IRABCD=INABCD(IVEC(IAB)+IVEC(ICD),
     &             JVEC(IAB)+JVEC(ICD),KVEC(IAB)+KVEC(ICD))
C
            IF(IRC(IRABCD).EQ.0) GOTO 797
C
            IF(ICDR11(ICD,JCMP).EQ.1) THEN
              DO M=1,MAXCD
                GABR11(M,IAB)=
     &           GABR11(M,IAB)-ERCD11(M,ICD,JCMP)*RC(M,IRABCD)
              ENDDO
            ENDIF
C
            IF(ICDI11(ICD,JCMP).EQ.1) THEN
              DO M=1,MAXCD
                GABI11(M,IAB)=
     &           GABI11(M,IAB)-EICD11(M,ICD,JCMP)*RC(M,IRABCD)
              ENDDO
            ENDIF
C
            IF(ICDR21(ICD,JCMP).EQ.1) THEN
              DO M=1,MAXCD
                GABR21(M,IAB)=
     &           GABR21(M,IAB)-ERCD21(M,ICD,JCMP)*RC(M,IRABCD)
              ENDDO
            ENDIF
C
            IF(ICDI21(ICD,JCMP).EQ.1) THEN
              DO M=1,MAXCD
                GABI21(M,IAB)=
     &           GABI21(M,IAB)-EICD21(M,ICD,JCMP)*RC(M,IRABCD)
              ENDDO
            ENDIF
C
797         CONTINUE
          ENDDO
C
        ENDIF
C
C-------------------------------------------------
C
        IF(ICMP.EQ.1) THEN
          IDX(1)=1
          IDX(2)=0
          IDX(3)=0
        ELSEIF(ICMP.EQ.2) THEN
          IDX(1)=0
          IDX(2)=1
          IDX(3)=0
        ELSEIF(ICMP.EQ.3) THEN
          IDX(1)=0
          IDX(2)=0
          IDX(3)=1
        ENDIF
C
        IF(JCMP.EQ.1) THEN
          JDX(1)=1
          JDX(2)=0
          JDX(3)=0
        ELSEIF(JCMP.EQ.2) THEN
          JDX(1)=0
          JDX(2)=1
          JDX(3)=0
        ELSEIF(JCMP.EQ.3) THEN
          JDX(1)=0
          JDX(2)=0
          JDX(3)=1
        ENDIF
C
C-------------------------------------------------
C
        DO ICD=1,NTUVCD
C
          IF(JCMP.EQ.1) THEN
            RTP=DFLOAT(IVEC(IAB)+IVEC(ICD))
          ELSEIF(JCMP.EQ.2) THEN
            RTP=DFLOAT(JVEC(IAB)+JVEC(ICD))
          ELSE
            RTP=DFLOAT(KVEC(IAB)+KVEC(ICD))
          ENDIF
C
          I1=IVEC(IAB)+IVEC(ICD)+IDX(1)+JDX(1)
          J1=JVEC(IAB)+JVEC(ICD)+IDX(2)+JDX(2)
          K1=KVEC(IAB)+KVEC(ICD)+IDX(3)+JDX(3)
          IADR1=INABCD(I1,J1,K1)
C
          I1=IVEC(IAB)+IVEC(ICD)+IDX(1)
          J1=JVEC(IAB)+JVEC(ICD)+IDX(2)
          K1=KVEC(IAB)+KVEC(ICD)+IDX(3)
          IADR2=INABCD(I1,J1,K1)
C
          I1=IVEC(IAB)+IVEC(ICD)+IDX(1)-JDX(1)
          J1=JVEC(IAB)+JVEC(ICD)+IDX(2)-JDX(2)
          K1=KVEC(IAB)+KVEC(ICD)+IDX(3)-JDX(3)
C
          IF((I1.GE.0).AND.(J1.GE.0).AND.(K1.GE.0)) THEN
            IADR3=INABCD(I1,J1,K1)
          ELSE
            IADR3=1
            RTP=0.0D0
          ENDIF
C
          IF(ICDR11(ICD,JCMP).EQ.1) THEN
            DO M=1,MAXCD
            GABR11(M,IAB)=GABR11(M,IAB)+ERCD11(M,ICD,JCMP)*
     &       (0.5D0*RC(M,IADR1)/APH(M)-
     &        PQ(M,JCMP)*RC(M,IADR2)+
     &        RTP*RC(M,IADR3))
            ENDDO
          ENDIF
C
          IF(ICDI11(ICD,JCMP).EQ.1) THEN
            DO M=1,MAXCD
            GABI11(M,IAB)=GABI11(M,IAB)+EICD11(M,ICD,JCMP)*
     &       (0.5D0*RC(M,IADR1)/APH(M)-
     &        PQ(M,JCMP)*RC(M,IADR2)+
     &        RTP*RC(M,IADR3))
            ENDDO
          ENDIF
C
          IF(ICDR21(ICD,JCMP).EQ.1) THEN
            DO M=1,MAXCD
            GABR21(M,IAB)=GABR21(M,IAB)+ERCD21(M,ICD,JCMP)*
     &       (0.5D0*RC(M,IADR1)/APH(M)-
     &        PQ(M,JCMP)*RC(M,IADR2)+
     &        RTP*RC(M,IADR3))
            ENDDO
          ENDIF
C
          IF(ICDI21(ICD,JCMP).EQ.1) THEN
            DO M=1,MAXCD
            GABI21(M,IAB)=GABI21(M,IAB)+EICD21(M,ICD,JCMP)*
     &       (0.5D0*RC(M,IADR1)/APH(M)-
     &        PQ(M,JCMP)*RC(M,IADR2)+
     &        RTP*RC(M,IADR3))
            ENDDO
          ENDIF
C
        ENDDO
C
6100    CONTINUE
C
      ENDDO
C
C****************************************************************C
C     GENERATE ALL POSSIBLE TWO-ELECTRON INTEGRALS FROM THE      C
C     EAB COEFFICIENTS AND THE G-ARRAYS                          C
C****************************************************************C
C
C     CALCULATE PHASES FOR BASIS FUNCTION OVERLAP COMBINATIONS
      PAB =-ISIGN(1,KQN(1)*KQN(2))*(-1)**((MQN(1)-MQN(2))/2)
      PCD =-ISIGN(1,KQN(3)*KQN(4))*(-1)**((MQN(3)-MQN(4))/2)
C
      PABCD = PAB*PCD
C
C****************************************************************C
C
C     Integral : ( - - || - - )
C
      DO M=1,MAXCD
        QR1(M)=0.0D0
        QI1(M)=0.0D0
        QR2(M)=0.0D0
        QI2(M)=0.0D0
      ENDDO
C
      DO IAB=1,NTUVAB
C
        IF(IABR11(IAB,ICMP).EQ.1) THEN
          DO M=1,MAXCD
            QR1(M)=QR1(M)+ERAB11(IJ,IAB,ICMP)*GABR11(M,IAB)
            QI2(M)=QI2(M)+ERAB11(IJ,IAB,ICMP)*GABI11(M,IAB)
          ENDDO
        ENDIF
C
        IF(IABI11(IAB,ICMP).EQ.1) THEN
          DO M=1,MAXCD
            QR2(M)=QR2(M)-EIAB11(IJ,IAB,ICMP)*GABI11(M,IAB)
            QI1(M)=QI1(M)+EIAB11(IJ,IAB,ICMP)*GABR11(M,IAB)
          ENDDO
        ENDIF
C
      ENDDO
C
C
      DO M=1,MAXCD
        RR(M, 1) = RR(M, 1) +     (QR1(M)+QR2(M))*PRE(M)
        RI(M, 1) = RI(M, 1) +     (QI1(M)+QI2(M))*PRE(M)
        RR(M, 4) = RR(M, 4) + PCD*(QR1(M)-QR2(M))*PRE(M)
        RI(M, 4) = RI(M, 4) + PCD*(QI1(M)-QI2(M))*PRE(M)
        RR(M,13) = PABCD*RR(M, 4)
        RI(M,13) =-PABCD*RI(M, 4)
        RR(M,16) = PABCD*RR(M, 1)
        RI(M,16) =-PABCD*RI(M, 1)
      ENDDO
C
C****************************************************************C
C
C     Integral : ( - - || + - )
C
      DO M=1,MAXCD
        QR1(M)=0.0D0
        QI1(M)=0.0D0
        QR2(M)=0.0D0
        QI2(M)=0.0D0
      ENDDO
C
      DO IAB=1,NTUVAB
C
        IF(IABR11(IAB,ICMP).EQ.1) THEN
          DO M=1,MAXCD
            QR1(M)=QR1(M)+ERAB11(IJ,IAB,ICMP)*GABR21(M,IAB)
            QI2(M)=QI2(M)+ERAB11(IJ,IAB,ICMP)*GABI21(M,IAB)
          ENDDO
        ENDIF
C
        IF(IABI11(IAB,ICMP).EQ.1) THEN
          DO M=1,MAXCD
            QR2(M)=QR2(M)-EIAB11(IJ,IAB,ICMP)*GABI21(M,IAB)
            QI1(M)=QI1(M)+EIAB11(IJ,IAB,ICMP)*GABR21(M,IAB)
          ENDDO
        ENDIF
C
      ENDDO
C
C
      DO M=1,MAXCD
        RR(M, 3) = RR(M, 3) +     (QR1(M)+QR2(M))*PRE(M)
        RI(M, 3) = RI(M, 3) +     (QI1(M)+QI2(M))*PRE(M)
        RR(M, 2) = RR(M, 2) - PCD*(QR1(M)-QR2(M))*PRE(M)
        RI(M, 2) = RI(M, 2) - PCD*(QI1(M)-QI2(M))*PRE(M)
        RR(M,15) =-PABCD*RR(M,2)
        RI(M,15) = PABCD*RI(M,2)
        RR(M,14) =-PABCD*RR(M,3)
        RI(M,14) = PABCD*RI(M,3)
      ENDDO
C
C****************************************************************C
C
C     Integral : ( + - || - - )
C
      DO M=1,MAXCD
        QR1(M)=0.0D0
        QI1(M)=0.0D0
        QR2(M)=0.0D0
        QI2(M)=0.0D0
      ENDDO
C
      DO IAB=1,NTUVAB
C
        IF(IABR21(IAB,ICMP).EQ.1) THEN
          DO M=1,MAXCD
            QR1(M)=QR1(M)+ERAB21(IJ,IAB,ICMP)*GABR11(M,IAB)
            QI2(M)=QI2(M)+ERAB21(IJ,IAB,ICMP)*GABI11(M,IAB)
          ENDDO
        ENDIF
C
        IF(IABI21(IAB,ICMP).EQ.1) THEN
          DO M=1,MAXCD
            QR2(M)=QR2(M)-EIAB21(IJ,IAB,ICMP)*GABI11(M,IAB)
            QI1(M)=QI1(M)+EIAB21(IJ,IAB,ICMP)*GABR11(M,IAB)
          ENDDO
        ENDIF
C
      ENDDO
C
C
      DO M=1,MAXCD
        RR(M, 9) = RR(M, 9) +     (QR1(M)+QR2(M))*PRE(M)
        RI(M, 9) = RI(M, 9) +     (QI1(M)+QI2(M))*PRE(M)
        RR(M,12) = RR(M,12) + PCD*(QR1(M)-QR2(M))*PRE(M)
        RI(M,12) = RI(M,12) + PCD*(QI1(M)-QI2(M))*PRE(M)
        RR(M, 5) =-PABCD*RR(M,12)
        RI(M, 5) = PABCD*RI(M,12)
        RR(M, 8) =-PABCD*RR(M, 9)
        RI(M, 8) = PABCD*RI(M, 9)
      ENDDO
C
C****************************************************************C
C
C     Integral : ( + - || + - )
C
      DO M=1,MAXCD
        QR1(M)=0.0D0
        QI1(M)=0.0D0
        QR2(M)=0.0D0
        QI2(M)=0.0D0
      ENDDO
C
      DO IAB=1,NTUVAB
C
        IF(IABR21(IAB,ICMP).EQ.1) THEN
          DO M=1,MAXCD
            QR1(M)=QR1(M)+ERAB21(IJ,IAB,ICMP)*GABR21(M,IAB)
            QI2(M)=QI2(M)+ERAB21(IJ,IAB,ICMP)*GABI21(M,IAB)
          ENDDO
        ENDIF
C
        IF(IABI21(IAB,ICMP).EQ.1) THEN
          DO M=1,MAXCD
          QR2(M)=QR2(M)-EIAB21(IJ,IAB,ICMP)*GABI21(M,IAB)
          QI1(M)=QI1(M)+EIAB21(IJ,IAB,ICMP)*GABR21(M,IAB)
          ENDDO
        ENDIF
C
      ENDDO
C
C
      DO M=1,MAXCD
        RR(M,11) = RR(M,11) +     (QR1(M)+QR2(M))*PRE(M)
        RI(M,11) = RI(M,11) +     (QI1(M)+QI2(M))*PRE(M)
        RR(M,10) = RR(M,10) - PCD*(QR1(M)-QR2(M))*PRE(M)
        RI(M,10) = RI(M,10) - PCD*(QI1(M)-QI2(M))*PRE(M)
        RR(M, 7) = PABCD*RR(M,10)
        RI(M, 7) =-PABCD*RI(M,10)
        RR(M, 6) = PABCD*RR(M,11)
        RI(M, 6) =-PABCD*RI(M,11)
      ENDDO
C
C
6000  CONTINUE
C
C     WRITE INTEGRALS INTO COMPLEX INTEGRAL ARRAY
C
      DO ITG=1,16
        DO M=1,MAXCD
          RINTG(M,ITG)=DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     *** INCLUDE THE OUTSIDE FACTOR OF (1/2) ***
C
      DO ITG=1,16
        DO M=1,MAXCD
          RINTG(M,ITG) = 0.5D0*DCMPLX(RR(M,ITG),RI(M,ITG))
        ENDDO
      ENDDO
C
      RETURN
      END
