      SUBROUTINE MULLIKN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     MM       MM UU    UU LL       LL       IIII KK    KK NN    NN    C
C     MMM     MMM UU    UU LL       LL        II  KK   KK  NNN   NN    C
C     MMMM   MMMM UU    UU LL       LL        II  KK  KK   NNNN  NN    C
C     MM MM MM MM UU    UU LL       LL        II  KKKKK    NN NN NN    C
C     MM  MMM  MM UU    UU LL       LL        II  KK  KK   NN  NNNN    C
C     MM   M   MM UU    UU LL       LL        II  KK   KK  NN   NNN    C
C     MM       MM  UUUUUU  LLLLLLLL LLLLLLLL IIII KK    KK NN    NN    C
C                                                                      C
C -------------------------------------------------------------------- C
C  MULLIKN CALCULATES A MULLIKEN POPULATION ANALYSIS ON A DIATOMIC     C
C  SYSTEM, AS DESCRIBED IN:                                            C
C  (*) J.Chem.Phys., 23: 1833, 1841, 2338, 2343 (1955).                C
C  (*) J.Chem.Phys., 36: 3428 (1962).                                  C
C -------------------------------------------------------------------- C
C  CONTRIBUTIONS FOR EACH ORBITAL ARE GROUPED BY ATOMIC CENTER, KAPPA  C
C  SYMMETRY AND MQN. ANY CHARGE DENSITY OTHER THAN THIS IS IN 'BORD'.  C
C**********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26)
C
      CHARACTER*2 ELMNT(120),ELA
      CHARACTER*4 HMLTN
C
      COMPLEX*16 C(MDM,MDM),OLAPLL(MDM,MDM),OLAPSS(MDM,MDM)
      COMPLEX*16 TMP1,TMP2
C
      DIMENSION FRAC(MDM,MCT,MKP,MMV,3)
      DIMENSION BORD(MDM,3),TOTL(MDM,3),RHO(MCT,3)
C
      COMMON/ATOM/ELMNT
      COMMON/COEF/C
      COMMON/PRMS/CV,HMLTN,ITER,IIL,ITREE,IEQS
      COMMON/QNMS/ICNLAB(MDM),KQNLAB(MDM),MQNLAB(MDM)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR,IOCCM0
C
C     INITIALISE COUNTER MATRICES
      DO IT=1,3
        DO IOCC=1,NOCC
          DO ICNT=1,NCNT
            DO IKAP=1,NKAP(ICNT)
              IKQN = KVALS(IKAP,ICNT)
              NMV  = IABS(IKQN)
              DO IMV=1,NMV
                FRAC(IOCC,ICNT,IKAP,IMV,IT) = 0.0D0
              ENDDO
            ENDDO
          ENDDO
          BORD(IOCC,IT) = 0.0D0
          TOTL(IOCC,IT) = 0.0D0
        ENDDO
        DO ICNT=1,NCNT
          RHO(ICNT,IT) = 0.0D0
        ENDDO
      ENDDO
C
C     GENERATE DIRECT OVERLAP MATRICES
      CALL MNPOLE(OLAPLL,1,0,1,2)
      CALL MNPOLE(OLAPSS,2,0,1,2)
C
C     LOOP OVER OCCUPIED ORBITALS
      DO IOCC=1,NOCC
C
C       IGNORE NEGATIVE SPECTRUM
        MOCC = IOCC + NSHIFT
C
C       LOOP OVER FOCK MATRIX ADDRESS BLOCKS
        DO I=1,NDIM-NSHIFT
          ICNT = ICNLAB(I)
          IKQN = KQNLAB(I)
          IMQN = MQNLAB(I)
          IF(IKQN.LT.0) THEN
            IKAP =-2*IKQN-1
          ELSE
            IKAP = 2*IKQN
          ENDIF
          IMV = (IMQN+1)/2
          DO J=1,NDIM-NSHIFT
            JCNT = ICNLAB(J)
            JKQN = KQNLAB(J)
            JMQN = MQNLAB(J)
C
C           LARGE AND SMALL CONTRIBUTIONS
            TMP1 = DCONJG(C(I       ,MOCC))*C(J       ,MOCC)*OLAPLL(I,J)
            TMP2 = DCONJG(C(I+NSHIFT,MOCC))*C(J+NSHIFT,MOCC)*OLAPSS(I,J)
C
C           UPDATE CHARGE ON CENTER ICNT
            RHO(ICNT,1) = RHO(ICNT,1) + DREAL(TMP1)
            RHO(ICNT,2) = RHO(ICNT,2) + DREAL(TMP2)
            RHO(ICNT,3) = RHO(ICNT,3) + DREAL(TMP1) + DREAL(TMP2)

C           DECIDE WHERE TO PUT CONTRIBUTION
            IF(ICNT.EQ.JCNT.AND.IKQN.EQ.JKQN.AND.IMQN.EQ.JMQN) THEN
              FRAC(IOCC,ICNT,IKAP,IMV,1) = FRAC(IOCC,ICNT,IKAP,IMV,1) 
     &                                   + DREAL(TMP1)
              FRAC(IOCC,ICNT,IKAP,IMV,2) = FRAC(IOCC,ICNT,IKAP,IMV,2) 
     &                                   + DREAL(TMP2)
              FRAC(IOCC,ICNT,IKAP,IMV,3) = FRAC(IOCC,ICNT,IKAP,IMV,1) 
     &                                   + FRAC(IOCC,ICNT,IKAP,IMV,2)
            ELSE
              BORD(IOCC,1) = BORD(IOCC,1) + DREAL(TMP1)
              BORD(IOCC,2) = BORD(IOCC,2) + DREAL(TMP2)
              BORD(IOCC,3) = BORD(IOCC,3) + DREAL(TMP1) + DREAL(TMP2)
            ENDIF
C          
          ENDDO
C
        ENDDO
C
C       TOTAL OCCUPANCIES
        DO ICNT=1,NCNT
            DO IKAP=1,NKAP(ICNT)
              IKQN = KVALS(IKAP,ICNT)
              NMV  = IABS(IKQN)
              DO IMV=1,NMV
                TOTL(IOCC,1) = TOTL(IOCC,1) + FRAC(IOCC,ICNT,IKAP,IMV,1)
                TOTL(IOCC,2) = TOTL(IOCC,2) + FRAC(IOCC,ICNT,IKAP,IMV,2)
                TOTL(IOCC,3) = TOTL(IOCC,3) + FRAC(IOCC,ICNT,IKAP,IMV,3)
              ENDDO
            ENDDO
        ENDDO
        TOTL(IOCC,1) = TOTL(IOCC,1) + BORD(IOCC,1)
        TOTL(IOCC,2) = TOTL(IOCC,2) + BORD(IOCC,2)
        TOTL(IOCC,3) = TOTL(IOCC,3) + BORD(IOCC,3)
C
C     END LOOP OVER OCCUPIED ORBITALS
      ENDDO
C
C     PRINT RESULTS OF ANALYSIS
      WRITE(6, *) 'Mulliken population analysis:'
      WRITE(7, *) 'Mulliken population analysis:'
      WRITE(6,*) REPEAT('-',62)
      WRITE(7,*) REPEAT('-',62)

10    FORMAT(1X,'Total charge on center ',I2,' = ',F15.10)
11    FORMAT(1X,'Total charge on molecule  = ',F15.10)
      SUM = 0.0D0
      DO ICNT=1,NCNT
        WRITE(6,10) ICNT,ZNUC(ICNT)-RHO(ICNT,3)
        WRITE(7,10) ICNT,ZNUC(ICNT)-RHO(ICNT,3)
        SUM = SUM + ZNUC(ICNT)-RHO(ICNT,3)
      ENDDO
      WRITE(6,11) SUM
      WRITE(7,11) SUM
      WRITE(6, *) REPEAT('-',62)
      WRITE(7, *) REPEAT('-',62)
C
12    FORMAT(' Orb.',2X,'Cent.',3X,'KQN',3X,'|MQN|',7X,
     &                                     'Q(L)',9X,'Q(S)',7X,'Q(TOT)')
13    FORMAT(1X,I3,2X,I2,' (',A,')',3X,I2,3X,I2,'/2',1X,F11.8,2X,
     &                                                   F11.8,2X,F11.8)
14    FORMAT(1X,I3,3X,I2,' (',A,')',2X,I2,3X,I2,'/2',1X,F11.8,2X,
     &                                                   F11.8,2X,F11.8)
15    FORMAT(1X,I3,3X,A ,13X,F11.8,2X,F11.8,2X,F11.8)
C
      WRITE(6,12)
      WRITE(7,12)
      DO IOCC=1,NOCC
        WRITE(6, *) REPEAT('-',62)
        WRITE(7, *) REPEAT('-',62)
        DO ICNT=1,NCNT
          ELA = ELMNT(IZNUC(ICNT))
          DO IKAP=1,NKAP(ICNT)
            IKQN = KVALS(IKAP,ICNT)
            NMV  = IABS(IKQN)
            DO IMV=1,NMV
              IF(FRAC(IOCC,ICNT,IKAP,IMV,3).GT.1.0D-9) THEN
                IF(NCNT.LT.10) THEN
                  WRITE(6,13) IOCC,ICNT,ELA,KVALS(IKAP,ICNT),2*IMV-1,
     &                              (FRAC(IOCC,ICNT,IKAP,IMV,IT),IT=1,3)
                  WRITE(7,13) IOCC,ICNT,ELA,KVALS(IKAP,ICNT),2*IMV-1,
     &                              (FRAC(IOCC,ICNT,IKAP,IMV,IT),IT=1,3)
                ELSE
                  WRITE(6,14) IOCC,ICNT,ELA,KVALS(IKAP,ICNT),2*IMV-1,
     &                              (FRAC(IOCC,ICNT,IKAP,IMV,IT),IT=1,3)
                  WRITE(7,14) IOCC,ICNT,ELA,KVALS(IKAP,ICNT),2*IMV-1,
     &                              (FRAC(IOCC,ICNT,IKAP,IMV,IT),IT=1,3)         
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        WRITE(6,15) IOCC,'Other ',(BORD(IOCC,IT),IT=1,3)
        WRITE(7,15) IOCC,'Other ',(BORD(IOCC,IT),IT=1,3)
        WRITE(6,15) IOCC,'Total ',(TOTL(IOCC,IT),IT=1,3)
        WRITE(7,15) IOCC,'Total ',(TOTL(IOCC,IT),IT=1,3)        
      ENDDO
C
      RETURN
      END
