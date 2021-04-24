      SUBROUTINE FOCKBAS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          FFFFFFF LL         AA    BBBBBBB  EEEEEEE LL                C
C          FF      LL        AAAA   BB    BB EE      LL                C
C          FF      LL       AA  AA  BB    BB EE      LL                C
C          FFFFF   LL      AA    AA BBBBBBB  EEEEE   LL                C
C          FF      LL      AAAAAAAA BB    BB EE      LL                C
C          FF      LL      AA    AA BB    BB EE      LL                C
C          FF      LLLLLLL AA    AA BBBBBBB  EEEEEEE LLLLLLL           C
C                                                                      C
C -------------------------------------------------------------------- C
C  FOCKBAS CALCUATES THE ADDRESSES OF THE FOCK MATRIX BLOCKS.          C
C**********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26)
C
      COMMON/LBL2/LDIAG(500),NDIG
      COMMON/QNMS/ICNBAS(MDM),KQNBAS(MDM),MQNBAS(MDM)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR,IOCCM0
C
C     QUANTUM NUMBER LABELS
      ILST = 0
      IDIG = 1
      DO ICNT=1,NCNT
        DO IM=1,MMV
          MQN = 2*IM-1
C
C         LABEL MQN<0 BLOCKS
          DO KA=1,NKAP(ICNT)
            KQN = KVALS(KA,ICNT)
            IF(KQN.GT.0) THEN
              LQN = KQN
            ELSE
              LQN =-KQN-1
            ENDIF
            NBAS  = NFUNCT(LQN+1,ICNT)
            MQMAX = 2*IABS(KQN)-1
            IF(MQMAX.GE.MQN) THEN
              LARGE(ICNT,KA,MQN) = ILST
              LDIAG(IDIG)        = ILST
              IDIG               = IDIG + 1
              DO IBAS=1,NBAS
                ICNBAS(ILST+IBAS) = ICNT
                KQNBAS(ILST+IBAS) = KQN
                MQNBAS(ILST+IBAS) =-MQN
              ENDDO
              ILST = ILST+NBAS
            ENDIF
          ENDDO
C
C         LABEL MQN>0 BLOCKS
          DO KA=1,NKAP(ICNT)
            KQN = KVALS(KA,ICNT)
            IF(KQN.GT.0) THEN
              LQN =  KQN
            ELSE
              LQN = -KQN-1
            ENDIF
            NBAS  = NFUNCT(LQN+1,ICNT)
            MQMAX = 2*IABS(KQN)-1
            IF(MQMAX.GE.MQN) THEN
              LARGE(ICNT,KA,MQN+1) = ILST
              LDIAG(IDIG)          = ILST
              IDIG                 = IDIG + 1
              DO IBAS=1,NBAS
                ICNBAS(ILST+IBAS) = ICNT
                KQNBAS(ILST+IBAS) = KQN
                MQNBAS(ILST+IBAS) = MQN
              ENDDO
              ILST = ILST + NBAS
            ENDIF
          ENDDO

        ENDDO
      ENDDO
      NDIG = IDIG-1
C
      RETURN
      END

