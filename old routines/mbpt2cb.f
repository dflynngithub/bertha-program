C
C
      SUBROUTINE MBPT2CB(MINO,NUMO,MINV,NUMV,G2INTA,G2INTB,IOSAV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   MM       MM BBBBBBB  PPPPPPP TTTTTTTT 222222   CCCCCC  BBBBBBB     C
C   MMM     MMM BB    BB PP    PP   TT   22    22 CC    CC BB    BB    C
C   MMMM   MMMM BB    BB PP    PP   TT         22 CC       BB    BB    C
C   MM MM MM MM BBBBBBB  PP    PP   TT       22   CC       BBBBBBB     C
C   MM  MMM  MM BB    BB PPPPPPP    TT     22     CC       BB    BB    C
C   MM   M   MM BB    BB PP         TT   22       CC    CC BB    BB    C
C   MM       MM BBBBBBB  PP         TT   22222222  CCCCCC  BBBBBBB     C
C                                                                      C
C -------------------------------------------------------------------- C
C  MBPT2CB READS CALCULATED MBPT1 AND MBPT2 RESULTS UNDER THE COULOMB  C
C  AND BREIT/GAUNT INTERACTIONS (OR ANY OTHER TWO-BODY OPERATORS), AND C
C  THEN CALCULATES DIRECT/EXCHANGE CONTRIBUTIONS FOR THEIR SUM.        C
C -------------------------------------------------------------------- C
C INPUT:                                                               C
C ▶ MINO   - LOWEST OCCUPIED STATE TO ACCOUNT FOR. (FULL: 1)           C
C ▶ NUMO   - NUMBER OF OCCUPIED STATES TO ACCOUNT FOR. (FULL: NOCC)    C
C ▶ MINV   - LOWEST VIRTUAL STATE TO ACCOUNT FOR. (FULL: NOCC+1)       C
C ▶ NUMV   - NUMBER OF VIRTUAL STATES TO ACCOUNT FOR. (FULL: NVRT)     C
C ▶ G2INTA - NAME OF 1ST 2-BODY OPERATOR ('COULM', 'GAUNT' OR 'BREIT').C
C ▶ G2INTB - NAME OF 2ND 2-BODY OPERATOR ('COULM', 'GAUNT' OR 'BREIT').C
C ▶ IOSAV  - I/O OPTION (READ RESULTS, WRITE RESULTS, DO NOTHING).     C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL FILETHERE
C
      CHARACTER*5  G2INT1,G2INT2
      CHARACTER*6  IOSAV
      CHARACTER*16 HMS
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*60 RASBOUT,EAB1OUT,EAB2OUT
C
      DIMENSION EAB1(NUMO,NUMO,6),EAB2(NUMO,NUMO,6)
      DIMENSION EA1(NUMO,6),EA1(NUMO,6)
C
      INTEGER*8 INEST,NNEST,NSB
C
      COMPLEX*16 CADB1((NUMO+1)*NUMO*NUMO*NUMO/2),
     &           CADB2((NUMO+1)*NUMO*NUMO*NUMO/2)
      COMPLEX*16 RASB1((NUMO+1)*NUMO*NUMV*NUMV/2),
     &           RASB2((NUMO+1)*NUMO*NUMV*NUMV/2)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGE/EIGN(MDM)
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
C
C     WARNINGS BASED ON INVALID HMLT VS. G2INT COMBINATIONS
10    FORMAT(1X,'In MBPT2CB: HMLT = ',A,' but G2INT1 = ',A,
     &                                          ' and G2INT2 = ',A,'.')
      ELSEIF(HMLT.EQ.'NORL'.OR.HMLT.EQ.'DHFR'.OR.HMLT.EQ.'DHFP') THEN
        IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') THEN
          WRITE(6,10) HMLT,G2INTA,G2INTB
          WRITE(7,10) HMLT,G2INTA,G2INTB
        ENDIF
      ENDIF
C
C     INITIALISE TIME COUNTERS
      TERI = 0.0D0
      TCN1 = 0.0D0
      TCN2 = 0.0D0
      TCN3 = 0.0D0
      TCN4 = 0.0D0
      TSUM = 0.0D0
C
      CALL CPU_TIME(TBEG)
C
C     INITIALISE THE ARRAY OF ENERGY VALUES
      DO IOCCA=1,NUMO
        DO IOCCB=1,IOCCA
          DO N=1,6
            EAB2(IOCCA,IOCCB,N) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     IMPORT MBPT1 PAIR RESULTS FOR E1(ab) OR INITIALISE ARRAY
      EAB1OUT = TRIM(OUTFL)//'_EAB1_'//G2INT//'.dat'
      INQUIRE(FILE=TRIM(EAB1OUT),EXIST=FILETHERE)
      IF(FILETHERE) THEN
        OPEN(UNIT=8,FILE=TRIM(EAB1OUT),STATUS='OLD')
        REWIND(UNIT=8)
        DO IOCCA=1,NUMO
          DO IOCCB=1,NUMO
            READ(8, *) Q1,Q2,Q3,(EAB2(IOCCA,IOCCB,N),N=1,3)
          ENDDO
        ENDDO
      ENDIF
      CLOSE(UNIT=8)
C
C     INITIALISE THE ARRAY FOR (AR|BS) VALUES
      M = 0
      DO IVRTR=1,NUMV
        DO IOCCA=1,NUMO
          DO IOCCB=1,IOCCA
            DO IVRTS=1,NUMV
              M = M+1
              RASB(M) = DCMPLX(0.0D0,0.0D0)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     TRY TO READ CABD ARRAY FROM EXISTING FILE
      IF(IOSAV.EQ.'READIN') THEN
        RASBOUT = TRIM(OUTFL)//'_RASB_'//G2INT1//'.dat'
        INQUIRE(FILE=TRIM(RASBOUT),EXIST=FILETHERE)
        IF(FILETHERE) THEN
          OPEN(UNIT=8,FILE=TRIM(RASBOUT),STATUS='OLD')
          REWIND(UNIT=8)
C         READ MINO AND NUMO
          READ(8, *) MINOTMP,NUMOTMP,MINVTMP,NUMVTMP
C         IF DIFFERENT TO INPUT, WRITE TO LOG AND UPDATE
          IF(MINOTMP.NE.MINO) THEN
            WRITE(6, *) 'In MBPT2CB: input MINO = ',MINO
            WRITE(6, *) '             file MINO = ',MINOTMP
          ENDIF
          IF(NUMOTMP.NE.NUMO) THEN
            WRITE(6, *) 'In MBPT2CB: input NUMO = ',NUMO
            WRITE(6, *) '             file NUMO = ',NUMOTMP
          ENDIF
          IF(MINVTMP.NE.MINV) THEN
            WRITE(6, *) 'In MBPT2CB: input MINV = ',MINV
            WRITE(6, *) '             file MINV = ',MINVTMP
          ENDIF
          IF(NUMVTMP.NE.NUMV) THEN
            WRITE(6, *) 'In MBPT2CB: input NUMV = ',NUMV
            WRITE(6, *) '             file NUMV = ',NUMVTMP
          ENDIF
          MINO = MINOTMP
          NUMO = NUMOTMP
          MINV = MINVTMP
          NUMV = NUMVTMP
C         LOOP OVER OCCUPIED STATES IOCCA AND VIRTUAL STATES IVRTR
          DO IOCCA=1,NUMO
            DO IVRTR=1,NUMV
C             LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
              DO IOCCB=1,IOCCA
                DO IVRTS=1,NUMV
C                 MAIN LIST ADDRESS FOR THIS IOCCA,IVRTR,IOCCB,IVRTS
                  MRASB = (IOCCA-1)*NUMV*IOCCA*NUMV/2
     &                  + (IVRTR-1)*NUMV*IOCCA + (IOCCB-1)*NUMV + IVRTS
                  READ(8, *) RASB(MRASB)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CLOSE(UNIT=8)
C         SKIP CALCULATION PHASE
          GOTO 100
        ELSE
C       IF THE FILE DOESN'T YET EXIST, CALCULATE AND SAVE IT
          WRITE(*, 6) 'In MBPT2CB: tried to read RASB but failed.
          WRITE(*, 7) 'In MBPT2CB: tried to read RASB but failed.
          IOSAV = 'WRTOUT'
        ENDIF
      ENDIF
C
C**********************************************************************C
C     CALCULATE SECOND-ORDER PAIR CORRELATION ENERGY                   C
C**********************************************************************C
C
C     WRITE CADB ARRAY TO EXTERNAL FILE IF PROMPTED
      IF(IOSAV.EQ.'WRTOUT') THEN
        RASBOUT = TRIM(OUTFL)//'_RASB_'//G2INT//'.dat'
        OPEN(UNIT=8,FILE=TRIM(RASBOUT),STATUS='UNKNOWN')
        REWIND(UNIT=8)
C       WRITE MINO AND NUMO
        WRITE(8, *) MINO,NUMO,MINV,NUMV
C       LOOP OVER OCCUPIED STATES IOCCA AND VIRTUAL STATES IVRTR
        DO IOCCA=1,NUMO
          DO IVRTR=1,NUMV
C           LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
            DO IOCCB=1,IOCCA
              DO IVRTS=1,NUMV
C               MAIN LIST ADDRESS FOR THIS IOCCA,IVRTR,IOCCB,IVRTS
                MRASB = (IOCCA-1)*NUMV*IOCCA*NUMV/2
     &                + (IVRTR-1)*NUMV*IOCCA + (IOCCB-1)*NUMV + IVRTS
                WRITE(8, *) RASB(MRASB)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CLOSE(UNIT=8)
      ENDIF
C
C     SKIP POINT FOR READIN OPTION
100   CONTINUE
C
      CALL CPU_TIME(T1)
C
C     FOR EACH IOCCA,IOCCB PAIR, SUM OVER IVRTR AND IVRTS CONTRIBUTIONS
C
C     LOOP OVER OCCUPIED STATES IOCCA AND VIRTUAL STATES IVRTR
      DO IOCCA=1,NUMO
        DO IVRTR=1,NUMV
C
C         FOCK MATRIX ADDRESSES
          IA = MINO-1+IOCCA+NSKP
          IR = MINV-1+IVRTR+NSKP
C
C         LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
          DO IOCCB=1,IOCCA
            DO IVRTS=1,NUMV
C
C             FOCK MATRIX ADDRESSES
              IB = MINO-1+IOCCB+NSKP
              IS = MINV-1+IVRTS+NSKP
C
C             MAIN LIST ADDRESS FOR THIS IOCCA,IVRTR,IOCCB,IVRTS
              MRASB = (IOCCA-1)*NUMV*IOCCA*NUMV/2
     &              + (IVRTR-1)*NUMV*IOCCA + (IOCCB-1)*NUMV + IVRTS
C
C             MAIN LIST ADDRESS FOR THIS IOCCA,IVRTS,IOCCB,IVRTR
              MSARB = (IOCCA-1)*NUMV*IOCCA*NUMV/2
     &              + (IVRTS-1)*NUMV*IOCCA + (IOCCB-1)*NUMV + IVRTR
C
C             NUMERATOR FOR DIRECT AND EXCHANGE CONTRIBUTIONS
              RNUMD = DREAL(RASB(MRASB)*DCONJG(RASB(MRASB)))
              RNUMX =-DREAL(RASB(MRASB)*DCONJG(RASB(MSARB)))
C
C             DENOMINATOR FOR BOTH CONTRIBUTIONS
              EABRS = EIGN(IA) + EIGN(IB) - EIGN(IR) - EIGN(IS)
C
C             ADD TO DIRECT AND EXCHANGE BINS
              EAB2(IOCCA,IOCCB,4) = EAB2(IOCCA,IOCCB,4) + RNUMD/EABRS
              EAB2(IOCCA,IOCCB,5) = EAB2(IOCCA,IOCCB,5) + RNUMX/EABRS
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     FILL IN THE OTHER HALF OF THE ARRAY AND CALCULATE TOTALS
      E1S = 0.0D0
      E2D = 0.0D0
      E2X = 0.0D0
      E2S = 0.0D0
      DO IOCCA=1,NUMO
        DO IOCCB=1,IOCCA
C
C         INTERMEDIATE VALUES
          EAB1TOT = EAB2(IOCCA,IOCCB,3)
          EAB2DIR = EAB2(IOCCA,IOCCB,4)
          EAB2XCH = EAB2(IOCCA,IOCCB,5)
          EAB2SUM = EAB2DIR + EAB2XCH
C
C         PUT THESE INTO EAB2 AND ADD CONTRIBUTION TO E2
          EAB2(IOCCA,IOCCB,6) = EAB2SUM
          IF(IOCCA.NE.IOCCB) THEN
            EAB2(IOCCB,IOCCA,3) = EAB1TOT
            EAB2(IOCCB,IOCCA,4) = EAB2DIR
            EAB2(IOCCB,IOCCA,5) = EAB2XCH
            EAB2(IOCCB,IOCCA,6) = EAB2SUM
            E1S = E1S +       EAB1TOT
            E2D = E2D +       EAB2DIR
            E2X = E2X +       EAB2XCH
            E2S = E2S +       EAB2SUM
          ELSE
            E1S = E1S + 0.5D0*EAB1TOT
            E2D = E2D + 0.5D0*EAB2DIR
            E2X = E2X + 0.5D0*EAB2XCH
            E2S = E2S + 0.5D0*EAB2SUM
          ENDIF
        ENDDO
      ENDDO
C
C     IMPORT MBPT1 PAIR RESULTS FOR E1(ab) OR INITIALISE ARRAY
      EAB2OUT = TRIM(OUTFL)//'_EAB2_'//G2INT//'.dat'
      INQUIRE(FILE=TRIM(EAB2OUT),EXIST=FILETHERE)
      IF(FILETHERE) THEN
        OPEN(UNIT=8,FILE=TRIM(EAB2OUT),STATUS='OLD')
        REWIND(UNIT=8)
        DO IOCCA=1,NUMO
          DO IOCCB=1,NUMO
            WRITE(8, *) (EAB2(IOCCA,IOCCB,N),N=1,6)
          ENDDO
        ENDDO
      ENDIF
      CLOSE(UNIT=8)
C
C**********************************************************************C
C     CALCULATE SECOND-ORDER SINGLE ORBITAL ENERGY                     C
C**********************************************************************C
C
C     FOR EACH IOCCA, SUM OVER THE IOCCB CONTRIBUTIONS
      DO IOCCA=1,NUMO
        DO N=1,6
          EA2(IOCCA,N) = 0.0D0
          DO IOCCB=1,NUMO
            EA2(IOCCA,N) = EA2(IOCCA,N) + EAB2(IOCCA,IOCCB,N)
          ENDDO
        ENDDO
      ENDDO
C
      CALL CPU_TIME(T2)
      TSUM = TSUM+T2-T1
C
C**********************************************************************C
C     TERMINAL OUTPUT SUMMARY                                          C
C**********************************************************************C
C
C     MBPT2CB PAIRWISE SUMMARY
20    FORMAT(1X,A,9X,A,9X,A,9X,A,9X,A)
21    FORMAT(' (',I2,',',I2,')',3X,F13.7,5X,F11.7,5X,F11.7,4X,F13.7)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',24),'MBPT2CB pairwise summary'
      WRITE(7, *) REPEAT(' ',24),'MBPT2CB pairwise summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) '( a, b)','E1G(ab)','E2J(ab)','E2K(ab)','E2G(ab)'
      WRITE(7,20) '( a, b)','E1G(ab)','E2J(ab)','E2K(ab)','E2G(ab)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO IOCCA=1,NUMO
        IANUM = IOCCA+MINO-1
        DO IOCCB=1,IOCCA
          IBNUM = IOCCB+MINO-1
          WRITE(6,21) IANUM,IBNUM,(EAB2(IOCCA,IOCCB,N),N=3,6)
          WRITE(7,21) IANUM,IBNUM,(EAB2(IOCCA,IOCCB,N),N=3,6)
        ENDDO
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     MBPT2CB SINGLE-PARTICLE SUMMARY
30    FORMAT(1X,A,10X,A,10X,A,10X,A,10X,A)
31    FORMAT('  ',I2,'    ',3X,F13.7,5X,F11.7,5X,F11.7,4X,F13.7)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',20),'MBPT2CB single particle summary'
      WRITE(7, *) REPEAT(' ',20),'MBPT2CB single particle summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,30) '  a    ','E1G(a)','E2J(a)','E2K(a)',' E2G(a)'
      WRITE(7,30) '  a    ','E1G(a)','E2J(a)','E2K(a)',' E2G(a)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO IOCCA=1,NUMO
        IANUM = IOCCA+MINO-1
        WRITE(6,31) IANUM,(EA2(IOCCA,N),N=3,6)
        WRITE(7,31) IANUM,(EA2(IOCCA,N),N=3,6)
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     MBPT2CB TOTAL SECOND ORDER CORRELATION ENERGY
32    FORMAT(' total  ',3X,F13.7,5X,F11.7,5X,F11.7,4X,F13.7)

      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',17),'MBPT2CB 2nd order correlation energy'
      WRITE(7, *) REPEAT(' ',17),'MBPT2CB 2nd order correlation energy'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,32) E1S,E2D,E2X,E2S
      WRITE(7,32) E1S,E2D,E2X,E2S
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     MBPT2 LABOUR ANALYSIS
      CALL CPU_TIME(TFIN)
      TTOT = TFIN - TBEG
      TOTH = TTOT - (TERI + TCN1 + TCN2 + TCN3 + TCN4 + TSUM)
C
      RETURN
      END


