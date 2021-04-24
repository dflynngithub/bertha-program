      SUBROUTINE MULLIKNKEEP(IVIR)
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
C  CONTRIBUTIONS FOR EACH ORBITAL ARE GROUPED BY ATOMIC CENTRE, KQN    C
C  SYMMETRY AND MQN. ANY CHARGE DENSITY OTHER THAN THIS IS IN 'BRD'.   C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    IVIR - NUMBER OF VIRTUAL ORBITALS TO INCLUDE IN LIST.             C
C**********************************************************************C
      INCLUDE 'param.h'
C
      CHARACTER*2 ELMT(120),ELA,ELB
      CHARACTER*4 HMLT
C
      DIMENSION SPNET(MCT,MCT,MDM),SPGRS(MCT,MDM)
      DIMENSION SPCNT(MDM),SPOLP(MDM),SPTOT(MDM)
      DIMENSION ATNET(MCT,MCT),ATOLP(MCT),ATGRS(MCT)
C
      COMPLEX*16 VL,VS
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 OLAPLL(MDM,MDM),OLAPSS(MDM,MDM)
C
      COMMON/EIGC/COEF
      COMMON/MDLV/ELMT
      COMMON/PRMS/HMLT,IOPT,IMOL,INEW,ILEV,ISML,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/QNMS/LABICN(MDM),LABKQN(MDM),LABMQN(MDM)
      COMMON/SPEC/BSET(MBS,MKP,MCT),BXYZ(3,MCT),ZNUC(MCT),ANUC(MCT),
     &            RNUC(MCT),PNUC,LRGE(MCT,MKP,MKP+1),NFNC(MKP,MCT),
     &            KAPA(MKP,MCT),IZNC(MCT),IQNC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSKP,NOCC,NVRT
C
C     GENERATE DIRECT OVERLAP MATRICES (ZEROTH MOMENT)
      CALL VMOMNT0(OLAPLL,1,0,1,2)
      IF(HMLT.NE.'NORL') THEN
        CALL VMOMNT0(OLAPSS,4,0,1,2)
      ENDIF
C
C     INTIALISE NET AND GROSS SPINOR CHARGES
      DO IOCC=1,NOCC+IVIR
        SPTOT(IOCC) = 0.0D0
        SPCNT(IOCC) = 0.0D0
        SPOLP(IOCC) = 0.0D0
        DO ICNT=1,NCNT
          SPGRS(ICNT,IOCC) = 0.0D0
          DO JCNT=1,NCNT
            SPNET(ICNT,JCNT,IOCC) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     INTIALISE NET AND GROSS ATOMIC CHARGES
      DO ICNT=1,NCNT
        ATGRS(ICNT) = 0.0D0
        ATOLP(ICNT) = 0.0D0
        DO JCNT=1,NCNT
          ATNET(ICNT,JCNT) = 0.0D0
        ENDDO
      ENDDO
C
C     LOOP OVER OCCUPIED (AND VIRTUAL) SPINORS
      DO IOCC=1,NOCC+IVIR
C
C       FOCK ADDRESS FOR THIS CENTRE
        IAD = IOCC+NSKP
C
C       NET SPINOR CHARGES: LOOP OVER ALL BASIS FUNCTIONS
        DO I=1,NDIM-NSKP
          DO J=1,NDIM-NSKP
C
C           SMALL-COMPONENT ADDRESSES
            K = I+NSKP
            L = J+NSKP
C
C           IDENTIFY ATOMIC CENTRES
            ICNT = LABICN(I)
            JCNT = LABICN(J)
C
C           LARGE-COMPONENT ELECTRON DENSITY
            VL = DCONJG(COEF(I,IAD))*COEF(J,IAD)*OLAPLL(I,J)
C
C           SMALL-COMPONENT ELECTRON DENSITY
            IF(HMLT.NE.'NORL') THEN
              VS = DCONJG(COEF(K,IAD))*COEF(L,IAD)*OLAPSS(I,J)
            ELSE
              VS = DCMPLX(0.0D0,0.0D0)
            ENDIF
C
C           UPDATE ELECTRONIC DENSITY
            SPNET(ICNT,JCNT,IOCC) = SPNET(ICNT,JCNT,IOCC) + DREAL(VL+VS)
C
          ENDDO
        ENDDO
C
C       GROSS SPINOR CHANGES: LOOP OVER ATOMIC CENTRES
        DO ICNT=1,NCNT
          SPGRS(ICNT,IOCC) = SPNET(ICNT,ICNT,IOCC)
          DO JCNT=1,NCNT
            IF(JCNT.NE.ICNT) THEN
              SPGRS(ICNT,IOCC) = SPGRS(ICNT,IOCC)+SPNET(ICNT,JCNT,IOCC)
              SPOLP(IOCC)      = SPOLP(IOCC)     +SPNET(ICNT,JCNT,IOCC)
            ELSE
              SPCNT(IOCC) = SPCNT(IOCC) + SPNET(ICNT,JCNT,IOCC)
            ENDIF
          ENDDO
          SPTOT(IOCC) = SPTOT(IOCC) + SPGRS(ICNT,IOCC)
        ENDDO
C
      ENDDO
C
C     NET ATOMIC CHARGES: LOOP OVER ALL PAIRS OF CENTRES
      DO ICNT=1,NCNT
        DO JCNT=1,NCNT
C
C         ADD CONTRIBUTION FROM EACH SPINOR TO ATOMIC CENTRE
          DO IOCC=1,NOCC
            ATNET(ICNT,JCNT) = ATNET(ICNT,JCNT) + SPNET(ICNT,JCNT,IOCC)
          ENDDO
C
        ENDDO
      ENDDO
C
C     GROSS ATOMIC CHARGES: LOOP OVER ALL CENTRES
      DO ICNT=1,NCNT
C
C       ADD CONTRIBUTION FROM EACH SPINOR TO ATOMIC CENTRE
        DO IOCC=1,NOCC
          ATGRS(ICNT) = ATGRS(ICNT) + SPGRS(ICNT,IOCC)
          DO JCNT=1,NCNT
            IF(ICNT.NE.JCNT) THEN
              ATOLP(ICNT) = ATOLP(ICNT) + SPNET(ICNT,JCNT,IOCC)
            ENDIF
          ENDDO
        ENDDO
C
      ENDDO
C
C     CENTRED ATOMIC CHARGE AND CROSS-CENTRE ATOMIC CHARGE
      ATCNT = 0.0D0
      ATCRS = 0.0D0
      DO ICNT=1,NCNT
        DO JCNT=1,NCNT
          IF(ICNT.EQ.JCNT) THEN
            ATCNT = ATCNT + ATNET(ICNT,JCNT)
          ENDIF
        ENDDO
        ATCRS = ATCRS + ATOLP(ICNT)
      ENDDO
C
C     TOTAL ATOMIC CHARGE
      ATTOT = 0.0D0
      DO ICNT=1,NCNT
        ATTOT = ATTOT + ATGRS(ICNT)
      ENDDO
C
C     FINAL LIST:
50    FORMAT(1X,A,1X,'|',3X,A,13X,A,2X,'|',3X,A,9X,A,2X,'|',9X,A)
51    FORMAT(1X,I3,2X'|',2X,I2,' (',A,')',2X,F12.8,2X,'|',25X,'|',
     &                                                        2X,F12.8)
52    FORMAT(1X,I3,2X,'|',2X,I2,' (',A,')',2X,F12.8,2X,'|',
     &                       2X,I2,' (',A,')',2X,F12.8,2X,'|',2X,F12.8)
53    FORMAT(6X,'|',25X,'|',2X,I2,' (',A,')',2X,F12.8,2X,'|')
54    FORMAT(6X,'|',2X,I2,' (',A,')',2X,F12.8,2X,'|',25X,'|',2X,F12.8)
55    FORMAT(6X,'|',2X,I2,' (',A,')',2X,F12.8,2X,'|',2X,I2,' (',A,')',
     &                                        2X,F12.8,2X,'|',2X,F12.8)
56    FORMAT(1X,'Tot.',' |',11X,F12.8,2X,'|',11X,F12.8,2X,'|',2X,F12.8)
      WRITE(6,50) 'Orb.','C(1)','Net','C(2)','Overlap','Gross'
      WRITE(7,50) 'Orb.','C(1)','Net','C(2)','Overlap','Gross'
      WRITE(6,*) REPEAT('=',72)
      WRITE(7,*) REPEAT('=',72)
C
      DO IOCC=1,NOCC+IVIR
C
C       SPINOR DENSITY ON EACH CENTRE
        DO ICNT=1,NCNT
C
C         ELEMENT NAME FOR ICNT
          ELA = ELMT(IZNC(ICNT))
C
C         OFF-DIAGONAL CENTRE COUNTER
          ICT = 0
          DO JCNT=1,NCNT
C
C           SPECIAL CASE FOR ATOMIC CALCULATION
            IF(NCNT.EQ.1) THEN
              WRITE(6,51) IOCC,ICNT,ELA,SPNET(ICNT,ICNT,IOCC),
     &                                                 SPGRS(ICNT,IOCC)
              WRITE(7,51) IOCC,ICNT,ELA,SPNET(ICNT,ICNT,IOCC),
     &                                                 SPGRS(ICNT,IOCC)
            ENDIF
C
C           GENERAL MOLECULE
            IF(JCNT.NE.ICNT) THEN
C
C             ELEMENT NAME FOR JCNT
              ELB = ELMT(IZNC(JCNT))
C
C             UPDATE OFF-DIAGONAL CENTRE COUNTER AND WRITE RESULTS
              ICT = ICT+1
              IF(ICT.EQ.1) THEN
                WRITE(6,52) IOCC,ICNT,ELA,SPNET(ICNT,ICNT,IOCC),
     &                  JCNT,ELB,SPNET(ICNT,JCNT,IOCC),SPGRS(ICNT,IOCC)
                WRITE(7,52) IOCC,ICNT,ELA,SPNET(ICNT,ICNT,IOCC),
     &                  JCNT,ELB,SPNET(ICNT,JCNT,IOCC),SPGRS(ICNT,IOCC)
              ELSE
                WRITE(6,53) JCNT,ELB,SPNET(ICNT,JCNT,IOCC)
                WRITE(7,53) JCNT,ELB,SPNET(ICNT,JCNT,IOCC)
              ENDIF
            ENDIF
C
          ENDDO
C
        ENDDO
C
        WRITE(6, *) REPEAT('-',72)
        WRITE(7, *) REPEAT('-',72)
        WRITE(6,56) SPCNT(IOCC),SPOLP(IOCC),SPTOT(IOCC)
        WRITE(7,56) SPCNT(IOCC),SPOLP(IOCC),SPTOT(IOCC)
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
C
      ENDDO
C
      WRITE(6, *) ''
      WRITE(7, *) ''
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',29),'MOLECULE TOTAL'
      WRITE(7, *) REPEAT(' ',29),'MOLECULE TOTAL'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,50) '    ','C(1)','Net','C(2)','Overlap','Gross'
      WRITE(7,50) '    ','C(1)','Net','C(2)','Overlap','Gross'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     SPINOR DENSITY ON EACH CENTRE
      DO ICNT=1,NCNT
C
C       ELEMENT NAME FOR ICNT
        ELA = ELMT(IZNC(ICNT))
C
C       OFF-DIAGONAL CENTRE COUNTER
        ICT = 0
        DO JCNT=1,NCNT
C
C         SPECIAL CASE FOR ATOMIC CALCULATION
          IF(NCNT.EQ.1) THEN
            WRITE(6,54) ICNT,ELA,ATNET(ICNT,ICNT),ATGRS(ICNT)
            WRITE(7,54) ICNT,ELA,ATNET(ICNT,ICNT),ATGRS(ICNT)
          ENDIF
C
C         GENERAL MOLECULE
          IF(JCNT.NE.ICNT) THEN
C
C           ELEMENT NAME FOR JCNT
            ELB = ELMT(IZNC(JCNT))
C
C           UPDATE OFF-DIAGONAL CENTRE COUNTER AND WRITE RESULTS
            ICT = ICT+1
            IF(ICT.EQ.1) THEN
              WRITE(6,55) ICNT,ELA,ATNET(ICNT,ICNT),
     &                            JCNT,ELB,ATNET(ICNT,JCNT),ATGRS(ICNT)
              WRITE(7,55) ICNT,ELA,ATNET(ICNT,ICNT),
     &                            JCNT,ELB,ATNET(ICNT,JCNT),ATGRS(ICNT)
            ELSE
              WRITE(6,53) JCNT,ELB,ATNET(ICNT,JCNT)
              WRITE(7,53) JCNT,ELB,ATNET(ICNT,JCNT)
            ENDIF
          ENDIF
C
        ENDDO
      ENDDO

      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,56) ATCNT,ATTOT-ATCNT,ATTOT
      WRITE(7,56) ATCNT,ATTOT-ATCNT,ATTOT
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) ''
      WRITE(7, *) ''
C

57    FORMAT(1X,A,2X,A,2X,'|',7X,A,8X,A,5X,A,7X,A)
58    FORMAT(5X,I2,7X,A,2X,'|',F12.8,2X,F12.8,2X,F12.8,2X,F12.8)
59    FORMAT(1X,A,12X,'|',F12.8,2X,F12.8,2X,F12.8,2X,F12.8)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,57) 'Centre','Element','Z_nuc','Q_atom','Q_overlap',
     &                                                         'Q_total'
      WRITE(7,57) 'Centre','Element','Z_nuc','Q_atom','Q_overlap',
     &                                                         'Q_total'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      ZTOT = 0.0D0
      DO ICNT=1,NCNT
        ZTOT = ZTOT+ZNUC(ICNT)
        WRITE(6,58) ICNT,ELMT(IZNC(ICNT)),ZNUC(ICNT),ATNET(ICNT,ICNT),
     &                               ATOLP(ICNT),ZNUC(ICNT)-ATGRS(ICNT)
        WRITE(7,58) ICNT,ELMT(IZNC(ICNT)),ZNUC(ICNT),ATNET(ICNT,ICNT),
     &                               ATOLP(ICNT),ZNUC(ICNT)-ATGRS(ICNT)
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,59) 'Total',ZTOT,ATCNT,ATCRS,ZTOT-ATTOT
      WRITE(7,59) 'Total',ZTOT,ATCNT,ATCRS,ZTOT-ATTOT
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) ''
      WRITE(7, *) ''
C
      RETURN
      END

