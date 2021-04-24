C
C  This was just after I wrote the analytic derivative routines (still
C  with a little bug in the latter). I decided to split them up into
C  separate trees.
C
C**********************************************************************C
C ==================================================================== C
C  [12] EQ-COEFFS: FINITE BASIS OVERLAP SPIN-STRUCTURE FACTORS.        C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] EQSAVE: MAIN ROUTINE FOR BUILDING A GLOBAL FILE OF EQ-COEFFS.  C
C   [B] E0LLGN: GENERATE FULL SET OF E0LL COEFFICIENTS AND SAVE.       C
C   [C] E0SSGN: GENERATE FULL SET OF E0SS COEFFICIENTS AND SAVE.       C
C   [D] EILSGN: GENERATE FULL SET OF EILS COEFFICIENTS AND SAVE.       C
C   [E] EISLGN: GENERATE FULL SET OF EISL COEFFICIENTS AND SAVE.       C
C -------------------------------------------------------------------- C
C   [A] EQLLMK: GENERATE A BATCH OF EQLL COEFFICIENTS: (--) AND (+-).  C
C   [B] EQLSMK: GENERATE A BATCH OF EQLS COEFFICIENTS: (--) AND (+-).  C
C   [C] EQSLMK: GENERATE A BATCH OF EQSL COEFFICIENTS: (--) AND (+-).  C
C   [D] EQSSMK: GENERATE A BATCH OF EQSS COEFFICIENTS: (--) AND (+-).  C
C   [E] GQLLMK: GENERATE A BATCH OF GQLL COEFFICIENTS: (--) AND (+-).  C
C   [F] GQLSMK: GENERATE A BATCH OF GQLS COEFFICIENTS: (--) AND (+-).  C
C   [G] GQSLMK: GENERATE A BATCH OF GQSL COEFFICIENTS: (--) AND (+-).  C
C   [H] GQSSMK: GENERATE A BATCH OF GQSS COEFFICIENTS: (--) AND (+-).  C
C   [I] EILSB3: GENERATE A VECTOR BATCH OF EILS COEFFICIENTS FOR BREIT.C
C   [J] EISLB3: GENERATE A VECTOR BATCH OF EISL COEFFICIENTS FOR BREIT.C
C   [K] EGQLL: A RAW BLOCK OF EQLL/GQLL COEFFICIENTS FOR EQLLMK/GQLLMK.C
C   [L] EGQLS: A RAW BLOCK OF EQLS/GQLS COEFFICIENTS FOR EQLSMK/GQLSMK.C
C   [M] EGQSL: A RAW BLOCK OF EQSL/GQSL COEFFICIENTS FOR EQSLMK/GQSLMK.C
C   [N] EGQSS: A RAW BLOCK OF EQSS/GQSS COEFFICIENTS FOR EQSSMK/GQSSMK.C
C   [O] ESGTF: SET OF ES-COEFFS OVER SPHERICAL HARMONICS AND HGTFS.    C
C   [P] GSGTF: SET OF GS-COEFFS OVER SPHERICAL HARMONICS AND HGTFS.    C
C   [Q] EVRS: EXPANSION COEFFS IN HGTF OVERLAPS, CALLED IN ESGTF.      C
C   [R] GVRS: EXPANSION COEFFS IN HGTF OVERLAPS, CALLED IN GSGTF.      C
C   [S] ESTEPLM: SIMULTANEOUS INCREASE IN (L,M) FOR USE IN EVRS.       C
C   [T] GSTEPLM: SIMULTANEOUS INCREASE IN (L,M) FOR USE IN GVRS.       C
C   [U] ESTEPL: INCREMENT IN L FOR USE IN EVRS.                        C
C   [V] GSTEPL: INCREMENT IN L FOR USE IN GVRS.                        C
C   [W] ESTEPN: INCREMENT IN N FOR USE IN EVRS.                        C
C   [X] GSTEPN: INCREMENT IN N FOR USE IN GVRS.                        C
C -------------------------------------------------------------------- C
C   [A] RNLL: A BLOCK OF LL NORMALISATION COEFFS.                      C
C   [B] RNSS: A BLOCK OF SS NORMALISATION COEFFS.                      C
C   [C] RNLS: A BLOCK OF LS NORMALISATION COEFFS.                      C
C   [D] DNORM: NORM FOR A REAL OR COMPLEX PART OF EQ-COEFF LIST.       C
C**********************************************************************C
C
C
      SUBROUTINE EQSAVE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        EEEEEEEE  QQQQQQ    SSSSSS     AA    VV    VV EEEEEEEE        C
C        EE       QQ    QQ  SS    SS   AAAA   VV    VV EE              C
C        EE      QQ      QQ SS        AA  AA  VV    VV EE              C
C        EEEEEE  QQ      QQ  SSSSSS  AA    AA VV    VV EEEEEE          C
C        EE      QQ      QQ       SS AAAAAAAA  VV  VV  EE              C
C        EE       QQ    QQ  SS    SS AA    AA   VVVV   EE              C
C        EEEEEEEE  QQQQQQ QQ SSSSSS  AA    AA    VV    EEEEEEEE        C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQSAVE CONSTRUCTS A SET OF COMMON ARRAYS FOR ALL REQUIRED EQTT'     C
C  COEFFICIENTS IN A CALCULATION THAT RESTS WITHIN QED.                C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      DIMENSION NTUV(0:MKP+1,4),NTRM(0:MKP+1,4),NWRD(0:MKP+1,4),
     &          SPCE(0:MKP+1,4)
      DIMENSION NTUVT(4),NTRMT(4),NWRDT(4),SPCET(4)
C
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/TMMD/TELL,TESS,TELS,TESL,TRLL,TRSS,TRLS,TRSL,TRBR
C
C     MAXIMUM REQUIRED COMPONENT TYPES
      IF(HMLT.EQ.'NORL') THEN
        ITTMAX = 1
      ELSEIF(HMLT.EQ.'BARE'.OR.HMLT.EQ.'DHFR') THEN
        ITTMAX = 2
      ELSEIF(TREE.NE.'MBPTN') THEN
        ITTMAX = 3
      ELSE
        ITTMAX = 4
      ENDIF
C
C**********************************************************************C
C     LOOP OVER FOCK BLOCK AND COUNT ALL REQUIRED EQ-WORDS.            C
C**********************************************************************C
C
C     INITIALISE MAXIMUM LAMBDA
      LAMMX = 0
C
C     INITIALISE TOTAL COEFFICIENT COUNTERS
      DO ITT=1,4
        DO LAM=0,MKP+1
          NTUV(LAM,ITT) = 0
          NTRM(LAM,ITT) = 0
          NWRD(LAM,ITT) = 0
          SPCE(LAM,ITT) = 0.0D0
        ENDDO
        NTUVT(ITT) = 0
        NTRMT(ITT) = 0
        NWRDT(ITT) = 0
        SPCET(ITT) = 0.0D0
      ENDDO
C
C     LOOP OVER CENTRES A AND B
      DO ICNTA=1,NCNT
        DO ICNTB=1,NCNT
C
C         LOOP OVER KQNA VALUES
          DO KA=1,NKAP(ICNTA)
C
C           QUANTUM NUMBERS FOR BLOCK A
            IF(KAPA(KA,ICNTA).LT.0) THEN
              LQNA =-KAPA(KA,ICNTA)-1
            ELSE
              LQNA = KAPA(KA,ICNTA)
            ENDIF
            NBASA = NFNC(LQNA,ICNTA)
C
C           LOOP OVER KQNB VALUES
            DO KB=1,NKAP(ICNTB)
C
C             QUANTUM NUMBERS FOR BLOCK B
              IF(KAPA(KB,ICNTB).LT.0) THEN
                LQNB =-KAPA(KB,ICNTB)-1
              ELSE
                LQNB = KAPA(KB,ICNTB)
              ENDIF
              NBASB = NFNC(LQNB,ICNTB)
C
C             LOOP OVER |MQNA| VALUES
              DO MA=1,IABS(KAPA(KA,ICNTA))
                MJA = 2*MA-1
C
C               LOOP OVER |MQNB| VALUES
                DO MB=1,IABS(KAPA(KB,ICNTB))
                  MJB = 2*MB-1
C
C                 NUMBER OF BASIS FUNCTION OVERLAPS
                  MAXAB = NBASA*NBASB
C
C                 CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
                  LAM  = LQNA+LQNB
C
C                 NUMBER OF TERMS IN EXPANSIONS OF THIS LENGTH
                  NTUV(LAM  ,1) = (LAM+1)*(LAM+2)*(LAM+3)/6
                  NTUV(LAM+2,2) = (LAM+3)*(LAM+4)*(LAM+5)/6
                  NTUV(LAM+1,3) = (LAM+2)*(LAM+3)*(LAM+4)/6
                  NTUV(LAM+1,4) = (LAM+2)*(LAM+3)*(LAM+4)/6
C
C                 INCREASE NUMBER OF WORDS FOR THIS LAMBDA VALUE
                  NTRM(LAM  ,1) = NTRM(LAM  ,1) + NTUV(LAM  ,1)*MAXAB
                  NTRM(LAM+2,2) = NTRM(LAM+2,2) + NTUV(LAM+2,2)*MAXAB
                  NTRM(LAM+1,3) = NTRM(LAM+1,3) + NTUV(LAM+1,3)*MAXAB
                  NTRM(LAM+1,4) = NTRM(LAM+1,4) + NTUV(LAM+1,4)*MAXAB
C
C                 UPDATE LARGEST LAMBDA VALUE
                  IF(LAM.GT.LAMMX) THEN
                    LAMMX = LAM
                  ENDIF
C
C               END LOOPS OVER |MQNA| AND |MQNB|
                ENDDO
              ENDDO
C
C           END LOOPS OVER KQNA AND KQNB
            ENDDO
          ENDDO
C
C       END LOOPS OVER ICNTA AND ICNTB
        ENDDO
      ENDDO
C
C     NUMBER OF WORDS IN SET AND SPACE REQUIRED
      DO LAM=0,LAMMX+2
        NWRD(LAM,1) =  4*NTRM(LAM,1)
        NWRD(LAM,2) =  4*NTRM(LAM,2)
        NWRD(LAM,3) = 12*NTRM(LAM,3)
        NWRD(LAM,4) = 12*NTRM(LAM,4)
        DO ITT=1,ITTMAX
          SPCE(LAM,ITT) = 8.0D-6*NWRD(LAM,ITT)
        ENDDO
      ENDDO
C
C     CALCULATE TOTALS
      DO ITT=1,ITTMAX
        DO LAM=0,LAMMX+2
          NTUVT(ITT) = NTUVT(ITT) + NTUV(LAM,ITT)
          NTRMT(ITT) = NTRMT(ITT) + NTRM(LAM,ITT)
          NWRDT(ITT) = NWRDT(ITT) + NWRD(LAM,ITT)
          SPCET(ITT) = SPCET(ITT) + SPCE(LAM,ITT)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     SUMMARY OF WORD ANALYSIS                                         C
C**********************************************************************C
C
C     SECTION TITLE
20    FORMAT(1X,A,9X,A,4X,A,6X,A,8X,A,8X,A,6X,A)
21    FORMAT(1X,A,8X,I2,4X,I5,3X,I9,7X,I2,3X,I10,5X,F10.3)
22    FORMAT(1X,A,13X,I10,7X,I2,3X,I10,5X,F10.3)
23    FORMAT(1X,A,9X,I2,3X,I5,3X,I9,7X,I2,3X,I10,5X,F10.3)
24    FORMAT(1X,A,5X,A,3X,A,6X,A,8X,A,8X,A,6X,A)
25    FORMAT(1X,A,40X,I10,5X,F10.3)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',22),'Eq-coefficient word analysis'
      WRITE(7, *) REPEAT(' ',22),'Eq-coefficient word analysis'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,20) 'Type','Λ','Terms','Length','#','Words','Size (MB)'
      WRITE(7,20) 'Type','Λ','Terms','Length','#','Words','Size (MB)'
C
C     E0LL ANALYSIS
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO LAM=0,LAMMX
        IF(LAM.EQ.0) THEN
          WRITE(6,21) 'E0LL',LAM,NTUV(LAM,1),NTRM(LAM,1), 4,NWRD(LAM,1),
     &                                                      SPCE(LAM,1)
          WRITE(7,21) 'E0LL',LAM,NTUV(LAM,1),NTRM(LAM,1), 4,NWRD(LAM,1),
     &                                                      SPCE(LAM,1)
        ELSE
          WRITE(6,21) '    ',LAM,NTUV(LAM,1),NTRM(LAM,1), 4,NWRD(LAM,1),
     &                                                      SPCE(LAM,1)
          WRITE(7,21) '    ',LAM,NTUV(LAM,1),NTRM(LAM,1), 4,NWRD(LAM,1),
     &                                                      SPCE(LAM,1)
        ENDIF
      ENDDO
C
      IF(HMLT.EQ.'NORL') GOTO 100
C
C     E0SS ANALYSIS
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO LAM=2,LAMMX+2
        IF(LAM.EQ.2) THEN
          WRITE(6,21) 'E0SS',LAM,NTUV(LAM,2),NTRM(LAM,2), 4,NWRD(LAM,2),
     &                                                      SPCE(LAM,2)
          WRITE(7,21) 'E0SS',LAM,NTUV(LAM,2),NTRM(LAM,2), 4,NWRD(LAM,2),
     &                                                      SPCE(LAM,2)
        ELSE
          WRITE(6,21) '    ',LAM,NTUV(LAM,2),NTRM(LAM,2), 4,NWRD(LAM,2),
     &                                                      SPCE(LAM,2)
          WRITE(7,21) '    ',LAM,NTUV(LAM,2),NTRM(LAM,2), 4,NWRD(LAM,2),
     &                                                      SPCE(LAM,2)
        ENDIF
      ENDDO
C
      IF(HMLT.EQ.'BARE'.OR.HMLT.EQ.'DHFR') GOTO 100
C
C     E0SS ANALYSIS
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO LAM=1,LAMMX+1
        IF(LAM.EQ.1) THEN
          WRITE(6,21) 'EILS',LAM,NTUV(LAM,3),NTRM(LAM,3),12,NWRD(LAM,3),
     &                                                      SPCE(LAM,3)
          WRITE(7,21) 'EILS',LAM,NTUV(LAM,3),NTRM(LAM,3),12,NWRD(LAM,3),
     &                                                      SPCE(LAM,3)
        ELSE
          WRITE(6,21) '    ',LAM,NTUV(LAM,3),NTRM(LAM,3),12,NWRD(LAM,3),
     &                                                      SPCE(LAM,3)
          WRITE(7,21) '    ',LAM,NTUV(LAM,3),NTRM(LAM,3),12,NWRD(LAM,3),
     &                                                      SPCE(LAM,3)
        ENDIF
      ENDDO
C
      IF(TREE.NE.'MBPTN') GOTO 100
C
C     EISL ANALYSIS
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO LAM=1,LAMMX+1
        IF(LAM.EQ.1) THEN
          WRITE(6,21) 'EISL',LAM,NTUV(LAM,4),NTRM(LAM,4),12,NWRD(LAM,4),
     &                                                      SPCE(LAM,4)
          WRITE(7,21) 'EISL',LAM,NTUV(LAM,4),NTRM(LAM,4),12,NWRD(LAM,4),
     &                                                      SPCE(LAM,4)
        ELSE
          WRITE(6,21) '    ',LAM,NTUV(LAM,4),NTRM(LAM,4),12,NWRD(LAM,4),
     &                                                      SPCE(LAM,4)
          WRITE(7,21) '    ',LAM,NTUV(LAM,4),NTRM(LAM,4),12,NWRD(LAM,4),
     &                                                      SPCE(LAM,4)
        ENDIF
      ENDDO
C
100   CONTINUE
C
C     SUMMARY OF TOTALS
      NWRDNET = 0
      SPCENET = 0.0D0
      DO ITT=1,ITTMAX
        NWRDNET = NWRDNET + NWRDT(ITT)
        SPCENET = SPCENET + SPCET(ITT)
      ENDDO
C
C     FIGURE OUT HOW MUCH SPACE IS ALLOWED BY PARAMETERS
      IF(HMLT.EQ.'NORL') THEN
        MWRD = 4
        MARR = 1
      ENDIF
      IF(HMLT.EQ.'DHFR') THEN
        MWRD = 8
        MARR = 2
      ENDIF
      IF(HMLT.EQ.'DHFP'.OR.HMLT.EQ.'DHFB'.OR.HMLT.EQ.'DHFQ') THEN
        MWRD = 20
        MARR = 3
      ENDIF
      IF(TREE.EQ.'MBPTN') THEN
        IF(HMLT.EQ.'DHFP'.OR.HMLT.EQ.'DHFB'.OR.HMLT.EQ.'DHFQ') THEN
          MWRD = 32
          MARR = 4
        ENDIF
      ENDIF
      SPCEMFL = 8.0D-6*MARR*MWRD*MFL
C
C     SUMMARISE TOTALS BY OVERLAP TYPE
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,24) 'Total','Λ_max','Terms','Length','#',
     &                                              'Words','Size (MB)'
      WRITE(7,24) 'Total','Λ_max','Terms','Length','#',
     &                                              'Words','Size (MB)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,23) 'E0LL',LAMMX  ,NTUVT(1),NTRMT(1), 4,NWRDT(1),SPCET(1)
      WRITE(7,23) 'E0LL',LAMMX  ,NTUVT(1),NTRMT(1), 4,NWRDT(1),SPCET(1)
      IF(HMLT.EQ.'NORL') GOTO 200
      WRITE(6,23) 'E0SS',LAMMX+2,NTUVT(2),NTRMT(2), 4,NWRDT(2),SPCET(2)
      WRITE(7,23) 'E0SS',LAMMX+2,NTUVT(2),NTRMT(2), 4,NWRDT(2),SPCET(2)
      IF(HMLT.EQ.'DHFR') GOTO 200
      WRITE(6,23) 'EILS',LAMMX+1,NTUVT(3),NTRMT(3),12,NWRDT(3),SPCET(3)
      WRITE(7,23) 'EILS',LAMMX+1,NTUVT(3),NTRMT(3),12,NWRDT(3),SPCET(3)
      IF(TREE.NE.'MBPTN') GOTO 200
      IF(HMLT.EQ.'NORL'.OR.HMLT.EQ.'DHFR') GOTO 200
      WRITE(6,23) 'EISL',LAMMX+1,NTUVT(4),NTRMT(4),12,NWRDT(4),SPCET(4)
      WRITE(7,23) 'EISL',LAMMX+1,NTUVT(4),NTRMT(4),12,NWRDT(4),SPCET(4)
200   CONTINUE
C      
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,25) '       ',NWRDNET,SPCENET
      WRITE(7,25) '       ',NWRDNET,SPCENET
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,22) 'parameters.h',MFL,MWRD,MWRD*MFL,SPCEMFL
      WRITE(7,22) 'parameters.h',MFL,MWRD,MWRD*MFL,SPCEMFL
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     OPTION WHEN NUMBER OF WORDS EXCEEDS ALLOCATED SIZE LIMIT
      IF(NTRMT(1).GT.MFL) THEN
        WRITE(6, *) 'In EQSAVE: E0LL words exceed allocated limit.'
        WRITE(7, *) 'In EQSAVE: E0LL words exceed allocated limit.'
        GOTO 150
      ENDIF
C
      IF(HMLT.NE.'NORL') THEN
        IF(NTRMT(2).GT.MFL) THEN
          WRITE(6, *) 'In EQSAVE: E0SS words exceed allocated limit.'
          WRITE(7, *) 'In EQSAVE: E0SS words exceed allocated limit.'
          GOTO 150
        ENDIF
      ENDIF
C     NO NEED TO ASK ABOUT EILS OR EISL! E0SS IS LARGER THAN EITHER.
C
C     SIZE LIMITS ARE ALL OK -- SKIP TO BATCH GENERATION
      GOTO 250
C
C     ONE OF THE CLASSES EXCEEDS WORD LIMIT
150   CONTINUE
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     HAVE TO GENERATE E COEFFICIENTS BY BATCH
      WRITE(6, *) 'In EQSAVE: Eq-coefficients to be generated by batch.'
      WRITE(7, *) 'In EQSAVE: Eq-coefficients to be generated by batch.'
C
C     FLIP THE EQ-GENERATION TOGGLE AND EXIT
      EQFILE = .FALSE.
      GOTO 300
C
250   CONTINUE
C
C**********************************************************************C
C     GENERATE COMPLETE BATCHES OF EQ-COEFFS                           C
C**********************************************************************C
C
C     SECTION TITLE
      WRITE(6, *) REPEAT(' ',18),'Generating Eq-coefficient data files'
      WRITE(7, *) REPEAT(' ',18),'Generating Eq-coefficient data files'
C
C     E0LL COEFFICIENTS
      CALL CPU_TIME(TDM1)
      CALL E0LLGN
      CALL CPU_TIME(TDM2)
      TELL = TELL+TDM2-TDM1
C
      IF(HMLT.EQ.'NORL') GOTO 300
C
C     E0SS COEFFICIENTS
      CALL E0SSGN
      CALL CPU_TIME(TDM3)
      TESS = TELL+TDM3-TDM2
C
      IF(HMLT.EQ.'BARE'.OR.HMLT.EQ.'DHFR') GOTO 300
C
C     EILS COEFFICIENTS
      CALL EILSGN
      CALL CPU_TIME(TDM4)
      TELS = TELS+TDM4-TDM3
C
      IF(TREE.NE.'MBPTN') GOTO 300
C
C     EISL COEFFICIENTS
      CALL EISLGN
      CALL CPU_TIME(TDM5)
      TESL = TESL+TDM5-TDM4
C
300   CONTINUE
C
C     END OF SECTION
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
      RETURN
      END
C
C
      SUBROUTINE E0LLGN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         EEEEEEEE 000000  LL       LL       GGGGGG  NN    NN          C
C         EE      00   000 LL       LL      GG    GG NNN   NN          C
C         EE      00  0000 LL       LL      GG       NNNN  NN          C
C         EEEEEE  00 00 00 LL       LL      GG       NN NN NN          C
C         EE      0000  00 LL       LL      GG   GGG NN  NNNN          C
C         EE      000   00 LL       LL      GG    GG NN   NNN          C
C         EEEEEEEE 000000  LLLLLLLL LLLLLLLL GGGGGG  NN    NN          C
C                                                                      C
C -------------------------------------------------------------------- C
C  E0LLGN GENERATES A FULL SET OF E0LL COEFFICIENTS AND SAVES TO A     C
C  COMMON ARRAY E0LL, INCLUDING AN ADDRESS INDEX.                      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/E0LL/E0LLFL(MFL,4),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     INITIALISE COUNT INDEX
      ICOUNT = 0
C
C     LOOP OVER CENTRE A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,ICNTA)
        IF(KQN(1).LT.0) THEN
          LQN(1) =-KQN(1)-1
        ELSE
          LQN(1) = KQN(1)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),ICNTA)
        DO IBAS=1,NBAS(1)
          EXL(IBAS,1) = BEXL(IBAS,LQN(1),ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KAPA(KB,ICNTB)
        IF(KQN(2).LT.0) THEN
          LQN(2) =-KQN(2)-1
        ELSE
          LQN(2) = KQN(2)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 2000 MA=1,IABS(KQN(1))
        MJA    = 2*MA-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 2000 MB=1,IABS(KQN(2))
        MJB    = 2*MB-1
        MQN(2) = MJB
C
C**********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS (MA,MB).  C
C     FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE ORDERED:             C
C -------------------------------------------------------------------- C
C  11: = (-|MQN(A)|,-|MQN(B)|) -> 1  12: = (-|MQN(A)|,+|MQN(B)|) -> 2  C
C  21: = (+|MQN(A)|,-|MQN(B)|) -> 3  22: = (+|MQN(A)|,+|MQN(B)|) -> 4  C
C**********************************************************************C
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NBAS(1)*NBAS(2)
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMAB  = LQN(1)+LQN(2)
      NTUVLL = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
C
C     INDEX TO TRACK START OF THIS BLOCK OF COEFFICIENTS
      IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB) = ICOUNT
C
C     GENERATE ELL0(AB) COEFFICIENTS
      CALL EQLLMK(E11,E21,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,0)
C
C     WRITE ELL0(AB) TO THE MAIN ARRAY
      DO ITUV=1,NTUVLL
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          E0LLFL(IAD+M,1) = DREAL(E11(M,ITUV))
          E0LLFL(IAD+M,2) = DIMAG(E11(M,ITUV))
          E0LLFL(IAD+M,3) = DREAL(E21(M,ITUV))
          E0LLFL(IAD+M,4) = DIMAG(E21(M,ITUV))
        ENDDO
      ENDDO
C
C     INCREASE COUNT INDEX
      ICOUNT = ICOUNT + NTUVLL*MAXAB
C
C     END LOOPS OVER FOCK BLOCK
2000  CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE E0SSGN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         EEEEEEEE 000000   SSSSSS   SSSSSS   GGGGGG  NN    NN         C
C         EE      00   000 SS    SS SS    SS GG    GG NNN   NN         C
C         EE      00  0000 SS       SS       GG       NNNN  NN         C
C         EEEEEE  00 00 00  SSSSSS   SSSSSS  GG       NN NN NN         C
C         EE      0000  00       SS       SS GG   GGG NN  NNNN         C
C         EE      000   00 SS    SS SS    SS GG    GG NN   NNN         C
C         EEEEEEEE 000000   SSSSSS   SSSSSS   GGGGGG  NN    NN         C
C                                                                      C
C -------------------------------------------------------------------- C
C  E0SSGN GENERATES A FULL SET OF E0SS COEFFICIENTS AND SAVES TO A     C
C  COMMON ARRAY E0SS, INCLUDING AN ADDRESS INDEX.                      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/E0SS/E0SSFL(MFL,4),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     INITIALISE COUNT INDEX
      ICOUNT = 0
C
C     LOOP OVER CENTRE A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,ICNTA)
        IF(KQN(1).LT.0) THEN
          LQN(1) =-KQN(1)-1
        ELSE
          LQN(1) = KQN(1)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),ICNTA)
        DO IBAS=1,NBAS(1)
          EXL(IBAS,1) = BEXL(IBAS,LQN(1),ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KAPA(KB,ICNTB)
        IF(KQN(2).LT.0) THEN
          LQN(2) =-KQN(2)-1
        ELSE
          LQN(2) = KQN(2)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 2000 MA=1,IABS(KQN(1))
        MJA    = 2*MA-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 2000 MB=1,IABS(KQN(2))
        MJB    = 2*MB-1
        MQN(2) = MJB
C
C**********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS (MA,MB).  C
C     FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE ORDERED:             C
C -------------------------------------------------------------------- C
C  11: = (-|MQN(A)|,-|MQN(B)|) -> 1  12: = (-|MQN(A)|,+|MQN(B)|) -> 2  C
C  21: = (+|MQN(A)|,-|MQN(B)|) -> 3  22: = (+|MQN(A)|,+|MQN(B)|) -> 4  C
C**********************************************************************C
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NBAS(1)*NBAS(2)
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMAB  = LQN(1)+LQN(2)+2
      NTUVSS = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
C
C     INDEX TO TRACK START OF THIS BLOCK OF COEFFICIENTS
      IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB) = ICOUNT
C
C     GENERATE ESS0(AB) COEFFICIENTS
      CALL EQSSMK(E11,E21,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,0)
C
C     WRITE ESS0(AB) TO THE MAIN ARRAY
      DO ITUV=1,NTUVSS
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          E0SSFL(IAD+M,1) = DREAL(E11(M,ITUV))
          E0SSFL(IAD+M,2) = DIMAG(E11(M,ITUV))
          E0SSFL(IAD+M,3) = DREAL(E21(M,ITUV))
          E0SSFL(IAD+M,4) = DIMAG(E21(M,ITUV))
        ENDDO
      ENDDO
C
C     INCREASE COUNT INDEX
      ICOUNT = ICOUNT + NTUVSS*MAXAB
C
C     END LOOPS OVER FOCK BLOCK
2000  CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE EILSGN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           EEEEEEEE IIII LL       SSSSSS   GGGGGG  NN    NN           C
C           EE        II  LL      SS    SS GG    GG NNN   NN           C
C           EE        II  LL      SS       GG       NNNN  NN           C
C           EEEEEE    II  LL       SSSSSS  GG       NN NN NN           C
C           EE        II  LL            SS GG   GGG NN  NNNN           C
C           EE        II  LL      SS    SS GG    GG NN   NNN           C
C           EEEEEEEE IIII LLLLLLLL SSSSSS   GGGGGG  NN    NN           C
C                                                                      C
C -------------------------------------------------------------------- C
C  EILSGN GENERATES A FULL SET OF EILS COEFFICIENTS AND SAVES TO A     C
C  COMMON ARRAY EILS, INCLUDING AN ADDRESS INDEX.                      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
C
      COMPLEX*16 E11X(MB2,MEQ),E11Y(MB2,MEQ),E11Z(MB2,MEQ),
     &           E21X(MB2,MEQ),E21Y(MB2,MEQ),E21Z(MB2,MEQ)
C
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EILS/EILSFL(MFL,12),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     INITIALISE COUNT INDEX
      ICOUNT = 0
C
C     LOOP OVER CENTRE A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,ICNTA)
        IF(KQN(1).LT.0) THEN
          LQN(1) =-KQN(1)-1
        ELSE
          LQN(1) = KQN(1)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),ICNTA)
        DO IBAS=1,NBAS(1)
          EXL(IBAS,1) = BEXL(IBAS,LQN(1),ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KAPA(KB,ICNTB)
        IF(KQN(2).LT.0) THEN
          LQN(2) =-KQN(2)-1
        ELSE
          LQN(2) = KQN(2)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 2000 MA=1,IABS(KQN(1))
        MJA    = 2*MA-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 2000 MB=1,IABS(KQN(2))
        MJB    = 2*MB-1
        MQN(2) = MJB
C
C**********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS (MA,MB).  C
C     FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE ORDERED:             C
C -------------------------------------------------------------------- C
C  11: = (-|MQN(A)|,-|MQN(B)|) -> 1  12: = (-|MQN(A)|,+|MQN(B)|) -> 2  C
C  21: = (+|MQN(A)|,-|MQN(B)|) -> 3  22: = (+|MQN(A)|,+|MQN(B)|) -> 4  C
C**********************************************************************C
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NBAS(1)*NBAS(2)
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMAB  = LQN(1)+LQN(2)+1
      NTUVLS = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
C
C     INDEX TO TRACK START OF THIS BLOCK OF COEFFICIENTS
      IADILS(ICNTA,ICNTB,KA,KB,MA,MB) = ICOUNT
C
C     GENERATE EILS(AB) COEFFICIENTS
      CALL EQLSMK(E11X,E21X,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,1)
      CALL EQLSMK(E11Y,E21Y,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,2)
      CALL EQLSMK(E11Z,E21Z,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,3)
C
C     WRITE EILS(AB) TO THE MAIN ARRAY
      DO ITUV=1,NTUVLS
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          EILSFL(IAD+M, 1) = DREAL(E11X(M,ITUV))
          EILSFL(IAD+M, 2) = DIMAG(E11X(M,ITUV))
          EILSFL(IAD+M, 3) = DREAL(E21X(M,ITUV))
          EILSFL(IAD+M, 4) = DIMAG(E21X(M,ITUV))
          EILSFL(IAD+M, 5) = DREAL(E11Y(M,ITUV))
          EILSFL(IAD+M, 6) = DIMAG(E11Y(M,ITUV))
          EILSFL(IAD+M, 7) = DREAL(E21Y(M,ITUV))
          EILSFL(IAD+M, 8) = DIMAG(E21Y(M,ITUV))
          EILSFL(IAD+M, 9) = DREAL(E11Z(M,ITUV))
          EILSFL(IAD+M,10) = DIMAG(E11Z(M,ITUV))
          EILSFL(IAD+M,11) = DREAL(E21Z(M,ITUV))
          EILSFL(IAD+M,12) = DIMAG(E21Z(M,ITUV))
        ENDDO
      ENDDO
C
C     INCREASE COUNT INDEX
      ICOUNT = ICOUNT + NTUVLS*MAXAB
C
C     END LOOPS OVER FOCK BLOCK
2000  CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE EISLGN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           EEEEEEEE IIII  SSSSSS  LL       GGGGGG  NN    NN           C
C           EE        II  SS    SS LL      GG    GG NNN   NN           C
C           EE        II  SS       LL      GG       NNNN  NN           C
C           EEEEEE    II   SSSSSS  LL      GG       NN NN NN           C
C           EE        II        SS LL      GG   GGG NN  NNNN           C
C           EE        II  SS    SS LL      GG    GG NN   NNN           C
C           EEEEEEEE IIII  SSSSSS  LLLLLLLL GGGGGG  NN    NN           C
C                                                                      C
C -------------------------------------------------------------------- C
C  EISLGN GENERATES A FULL SET OF EISL COEFFICIENTS AND SAVES TO A     C
C  COMMON ARRAY EISL, INCLUDING AN ADDRESS INDEX.                      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
C
      COMPLEX*16 E11X(MB2,MEQ),E11Y(MB2,MEQ),E11Z(MB2,MEQ),
     &           E21X(MB2,MEQ),E21Y(MB2,MEQ),E21Z(MB2,MEQ)
C
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EISL/EISLFL(MFL,12),IADISL(MCT,MCT,MKP,MKP,MKP,MKP)
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     INITIALISE COUNT INDEX
      ICOUNT = 0
C
C     LOOP OVER CENTRE A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,ICNTA)
        IF(KQN(1).LT.0) THEN
          LQN(1) =-KQN(1)-1
        ELSE
          LQN(1) = KQN(1)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),ICNTA)
        DO IBAS=1,NBAS(1)
          EXL(IBAS,1) = BEXL(IBAS,LQN(1),ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KAPA(KB,ICNTB)
        IF(KQN(2).LT.0) THEN
          LQN(2) =-KQN(2)-1
        ELSE
          LQN(2) = KQN(2)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 2000 MA=1,IABS(KQN(1))
        MJA    = 2*MA-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 2000 MB=1,IABS(KQN(2))
        MJB    = 2*MB-1
        MQN(2) = MJB
C
C**********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS (MA,MB).  C
C     FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE ORDERED:             C
C -------------------------------------------------------------------- C
C  11: = (-|MQN(A)|,-|MQN(B)|) -> 1  12: = (-|MQN(A)|,+|MQN(B)|) -> 2  C
C  21: = (+|MQN(A)|,-|MQN(B)|) -> 3  22: = (+|MQN(A)|,+|MQN(B)|) -> 4  C
C**********************************************************************C
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NBAS(1)*NBAS(2)
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMAB  = LQN(1)+LQN(2)+1
      NTUVSL = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
C
C     INDEX TO TRACK START OF THIS BLOCK OF COEFFICIENTS
      IADISL(ICNTA,ICNTB,KA,KB,MA,MB) = ICOUNT
C
C     GENERATE EISL(AB) COEFFICIENTS
      CALL EQSLMK(E11X,E21X,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,1)
      CALL EQSLMK(E11Y,E21Y,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,2)
      CALL EQSLMK(E11Z,E21Z,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,3)
C
C     WRITE EISL(AB) TO THE MAIN ARRAY
      DO ITUV=1,NTUVSL
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          EISLFL(IAD+M, 1) = DREAL(E11X(M,ITUV))
          EISLFL(IAD+M, 2) = DIMAG(E11X(M,ITUV))
          EISLFL(IAD+M, 3) = DREAL(E21X(M,ITUV))
          EISLFL(IAD+M, 4) = DIMAG(E21X(M,ITUV))
          EISLFL(IAD+M, 5) = DREAL(E11Y(M,ITUV))
          EISLFL(IAD+M, 6) = DIMAG(E11Y(M,ITUV))
          EISLFL(IAD+M, 7) = DREAL(E21Y(M,ITUV))
          EISLFL(IAD+M, 8) = DIMAG(E21Y(M,ITUV))
          EISLFL(IAD+M, 9) = DREAL(E11Z(M,ITUV))
          EISLFL(IAD+M,10) = DIMAG(E11Z(M,ITUV))
          EISLFL(IAD+M,11) = DREAL(E21Z(M,ITUV))
          EISLFL(IAD+M,12) = DIMAG(E21Z(M,ITUV))
        ENDDO
      ENDDO
C
C     INCREASE COUNT INDEX
      ICOUNT = ICOUNT + NTUVSL*MAXAB
C
C     END LOOPS OVER FOCK BLOCK
2000  CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE EQLLMK(E11,E21,EXL,XYZ,KQN,MQN,NBS,IPHS,IA1,IA2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      EEEEEEEE  QQQQQQ    LL       LL       MM       MM KK    KK      C
C      EE       QQ    QQ   LL       LL       MMM     MMM KK   KK       C
C      EE      QQ      QQ  LL       LL       MMMM   MMMM KK  KK        C
C      EEEEEE  QQ      QQ  LL       LL       MM MM MM MM KKKKK         C
C      EE      QQ      QQ  LL       LL       MM  MMM  MM KK  KK        C
C      EE       QQ    QQ   LL       LL       MM   M   MM KK   KK       C
C      EEEEEEEE  QQQQQQ QQ LLLLLLLL LLLLLLLL MM       MM KK    KK      C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQLLMK GENERATES A BATCH OF EQLL COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.        C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - BASIS INDICES TO CONSTRUCT EQ-BLOCK FROM.               C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 G11(MB2,MEQ),G21(MB2,MEQ)
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,4),XY2(3,4),KQ2(4),MQ2(4),NB2(4)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        XY2(IX,1) = XYZ(IX,IA1)
        XY2(IX,2) = XYZ(IX,IA2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KQ2(1) = KQN(IA1)
      KQ2(2) = KQN(IA2)
C
C     CALCULATE LQN VALUES
      IF(KQ2(1).LT.0) THEN
        LA =-KQ2(1)-1
      ELSE
        LA = KQ2(1)
      ENDIF
      IF(KQ2(2).LT.0) THEN
        LB =-KQ2(2)-1
      ELSE
        LB = KQ2(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NB2(1) = NBS(IA1)
      DO IBAS=1,NB2(1)
        EX2(IBAS,1) = EXL(IBAS,IA1)
      ENDDO
C
      NB2(2) = NBS(IA2)
      DO JBAS=1,NB2(2)
        EX2(JBAS,2) = EXL(JBAS,IA2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMLL  = LA+LB
      NTUVLL = (LAMLL+1)*(LAMLL+2)*(LAMLL+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR MQN PAIRS (-|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) =-MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQLL(E11,G11,EX2,XY2,KQ2,MQ2,NB2,IQ,0,0,LR)
C                                    ~
C     MULTIPLY EQLL BY PHASE TERM IF EQLL ARE NEEDED
      DO ITUV=1,NTUVLL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR MQN PAIRS (+|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) = MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQLL(E21,G21,EX2,XY2,KQ2,MQ2,NB2,IQ,0,0,LR)
C                                    ~
C     MULTIPLY EQLL BY PHASE TERM IF EQLL ARE NEEDED
      DO ITUV=1,NTUVLL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV) = PHS*E21(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT E22 AND E12 ARE RELATED TO THESE BY PHASE FACTORS.
C
      RETURN
      END
C
C
      SUBROUTINE EQLSMK(E11,E21,EXL,XYZ,KQN,MQN,NBS,IPHS,IA1,IA2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      EEEEEEEE  QQQQQQ    LL       SSSSSS  MM       MM KK    KK       C
C      EE       QQ    QQ   LL      SS    SS MMM     MMM KK   KK        C
C      EE      QQ      QQ  LL      SS       MMMM   MMMM KK  KK         C
C      EEEEEE  QQ      QQ  LL       SSSSSS  MM MM MM MM KKKKK          C
C      EE      QQ      QQ  LL            SS MM  MMM  MM KK  KK         C
C      EE       QQ    QQ   LL      SS    SS MM   M   MM KK   KK        C
C      EEEEEEEE  QQQQQQ QQ LLLLLLLL SSSSSS  MM       MM KK    KK       C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQLSMK GENERATES A BATCH OF EQLS COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB+1 -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.      C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - BASIS INDICES TO CONSTRUCT EQ-BLOCK FROM.               C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 G11(MB2,MEQ),G21(MB2,MEQ)
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        XY2(IX,1) = XYZ(IX,IA1)
        XY2(IX,2) = XYZ(IX,IA2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KQ2(1) = KQN(IA1)
      KQ2(2) = KQN(IA2)
C
C     CALCULATE LQN VALUES
      IF(KQ2(1).LT.0) THEN
        LA =-KQ2(1)-1
      ELSE
        LA = KQ2(1)
      ENDIF
      IF(KQ2(2).LT.0) THEN
        LB =-KQ2(2)-1
      ELSE
        LB = KQ2(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NB2(1) = NBS(IA1)
      DO IBAS=1,NB2(1)
        EX2(IBAS,1) = EXL(IBAS,IA1)
      ENDDO
C
      NB2(2) = NBS(IA2)
      DO JBAS=1,NB2(2)
        EX2(JBAS,2) = EXL(JBAS,IA2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMLS  = LA+LB+1
      NTUVLS = (LAMLS+1)*(LAMLS+2)*(LAMLS+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR MQN PAIRS (-|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) =-MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQLS(E11,G11,EX2,XY2,KQ2,MQ2,NB2,IQ,0,0,LR)
C                                    ~
C     MULTIPLY EQLS BY PHASE TERM IF EQLS ARE NEEDED
      DO ITUV=1,NTUVLS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR MQN PAIRS (+|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) = MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQLS(E21,G21,EX2,XY2,KQ2,MQ2,NB2,IQ,0,0,LR)
C                                    ~
C     MULTIPLY EQLS BY PHASE TERM IF EQLS ARE NEEDED
      DO ITUV=1,NTUVLS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV) = PHS*E21(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT E22 AND E12 ARE RELATED TO THESE BY PHASE FACTORS.
C     THERE IS ALSO A RELATION WHICH ALLOWS US TO OBTAIN ESLQ FROM ELSQ.
C
      RETURN
      END
C
C
      SUBROUTINE EQSLMK(E11,E21,EXL,XYZ,KQN,MQN,NBS,IPHS,IA1,IA2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      EEEEEEEE  QQQQQQ    SSSSSS  LL       MM       MM KK    KK       C
C      EE       QQ    QQ  SS    SS LL       MMM     MMM KK   KK        C
C      EE      QQ      QQ SS       LL       MMMM   MMMM KK  KK         C
C      EEEEEE  QQ      QQ  SSSSSS  LL       MM MM MM MM KKKKK          C
C      EE      QQ      QQ       SS LL       MM  MMM  MM KK  KK         C
C      EE       QQ    QQ  SS    SS LL       MM   M   MM KK   KK        C
C      EEEEEEEE  QQQQQQ QQ SSSSSS  LLLLLLLL MM       MM KK    KK       C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQSLMK GENERATES A BATCH OF EQSL COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB+1 -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.      C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - BASIS INDICES TO CONSTRUCT EQ-BLOCK FROM.               C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 G11(MB2,MEQ),G21(MB2,MEQ)
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        XY2(IX,1) = XYZ(IX,IA1)
        XY2(IX,2) = XYZ(IX,IA2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KQ2(1) = KQN(IA1)
      KQ2(2) = KQN(IA2)
C
C     CALCULATE LQN VALUES
      IF(KQ2(1).LT.0) THEN
        LA =-KQ2(1)-1
      ELSE
        LA = KQ2(1)
      ENDIF
      IF(KQ2(2).LT.0) THEN
        LB =-KQ2(2)-1
      ELSE
        LB = KQ2(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NB2(1) = NBS(IA1)
      DO IBAS=1,NB2(1)
        EX2(IBAS,1) = EXL(IBAS,IA1)
      ENDDO
C
      NB2(2) = NBS(IA2)
      DO JBAS=1,NB2(2)
        EX2(JBAS,2) = EXL(JBAS,IA2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMSL  = LA+LB+1
      NTUVSL = (LAMSL+1)*(LAMSL+2)*(LAMSL+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR MQN PAIRS (-|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) =-MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQSL(E11,G11,EX2,XY2,KQ2,MQ2,NB2,IQ,0,0,LR)
C                                    ~
C     MULTIPLY EQSL BY PHASE TERM IF EQSL ARE NEEDED
      DO ITUV=1,NTUVSL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR MQN PAIRS (+|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) = MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQSL(E21,G21,EX2,XY2,KQ2,MQ2,NB2,IQ,0,0,LR)
C                                    ~
C     MULTIPLY EQSL BY PHASE TERM IF EQSL ARE NEEDED
      DO ITUV=1,NTUVSL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV) = PHS*E21(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT E22 AND E12 ARE RELATED TO THESE BY PHASE FACTORS.
C
      RETURN
      END
C
C
      SUBROUTINE EQSSMK(E11,E21,EXL,XYZ,KQN,MQN,NBS,IPHS,IA1,IA2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      EEEEEEEE  QQQQQQ    SSSSSS   SSSSSS  MM       MM KK    KK       C
C      EE       QQ    QQ  SS    SS SS    SS MMM     MMM KK   KK        C
C      EE      QQ      QQ SS       SS       MMMM   MMMM KK  KK         C
C      EEEEEE  QQ      QQ  SSSSSS   SSSSSS  MM MM MM MM KKKKK          C
C      EE      QQ      QQ       SS       SS MM  MMM  MM KK  KK         C
C      EE       QQ    QQ  SS    SS SS    SS MM   M   MM KK   KK        C
C      EEEEEEEE  QQQQQQ QQ SSSSSS   SSSSSS  MM       MM KK    KK       C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQSSMK GENERATES A BATCH OF EQSS COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB+2 -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.      C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - BASIS INDICES TO CONSTRUCT EQ-BLOCK FROM.               C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 G11(MB2,MEQ),G21(MB2,MEQ)
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        XY2(IX,1) = XYZ(IX,IA1)
        XY2(IX,2) = XYZ(IX,IA2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KQ2(1) = KQN(IA1)
      KQ2(2) = KQN(IA2)
C
C     CALCULATE LQN VALUES
      IF(KQ2(1).LT.0) THEN
        LA =-KQ2(1)-1
      ELSE
        LA = KQ2(1)
      ENDIF
      IF(KQ2(2).LT.0) THEN
        LB =-KQ2(2)-1
      ELSE
        LB = KQ2(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NB2(1) = NBS(IA1)
      DO IBAS=1,NB2(1)
        EX2(IBAS,1) = EXL(IBAS,IA1)
      ENDDO
C
      NB2(2) = NBS(IA2)
      DO JBAS=1,NB2(2)
        EX2(JBAS,2) = EXL(JBAS,IA2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMSS  = LA+LB+2
      NTUVSS = (LAMSS+1)*(LAMSS+2)*(LAMSS+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR MQN PAIRS (-|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) =-MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQSS(E11,G11,EX2,XY2,KQ2,MQ2,NB2,IQ,0,0,LR)
C                                    ~
C     MULTIPLY EQSS BY PHASE TERM IF EQSS ARE NEEDED
      DO ITUV=1,NTUVSS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR MQN PAIRS (+|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) = MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQSS(E21,G21,EX2,XY2,KQ2,MQ2,NB2,IQ,0,0,LR)
C                                    ~
C     MULTIPLY EQSS BY PHASE TERM IF EQSS ARE NEEDED
      DO ITUV=1,NTUVSS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV) = PHS*E21(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT E22 AND E12 ARE RELATED TO THESE BY SIMPLE FACTORS.
C
      RETURN
      END
C
C
      SUBROUTINE EILSB3(E11,E21,EXL,XYZ,KQN,MQN,NBS,IPHS,IA1,IA2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           EEEEEEEE IIII LL       SSSSSS  BBBBBBB   333333            C
C           EE        II  LL      SS    SS BB    BB 33    33           C
C           EE        II  LL      SS       BB    BB       33           C
C           EEEEEE    II  LL       SSSSSS  BBBBBBB    33333            C
C           EE        II  LL            SS BB    BB       33           C
C           EE        II  LL      SS    SS BB    BB 33    33           C
C           EEEEEEEE IIII LLLLLLLL SSSSSS  BBBBBBB   333333            C
C                                                                      C
C -------------------------------------------------------------------- C
C  EILSB3 GENERATES A BATCH OF EILS COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB+1 -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.      C
C  THIS ACTUALLY MAKES A VECTOR LIST OF EQLS NEEDED FOR BREIT.         C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 E11(MB2,MEQ,3),E21(MB2,MEQ,3)
      COMPLEX*16 G11(MB2,MEQ),G21(MB2,MEQ)
      COMPLEX*16 E11X(MB2,MEQ),E11Y(MB2,MEQ),E11Z(MB2,MEQ),
     &           E21X(MB2,MEQ),E21Y(MB2,MEQ),E21Z(MB2,MEQ)
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        XY2(IX,1) = XYZ(IX,IA1)
        XY2(IX,2) = XYZ(IX,IA2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KQ2(1) = KQN(IA1)
      KQ2(2) = KQN(IA2)
C
C     CALCULATE LQN VALUES
      IF(KQ2(1).LT.0) THEN
        LA =-KQ2(1)-1
      ELSE
        LA = KQ2(1)
      ENDIF
      IF(KQ2(2).LT.0) THEN
        LB =-KQ2(2)-1
      ELSE
        LB = KQ2(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NB2(1) = NBS(IA1)
      DO IBAS=1,NB2(1)
        EX2(IBAS,1) = EXL(IBAS,IA1)
      ENDDO
C
      NB2(2) = NBS(IA2)
      DO JBAS=1,NB2(2)
        EX2(JBAS,2) = EXL(JBAS,IA2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMLS  = LA+LB+1
      NTUVLS = (LAMLS+1)*(LAMLS+2)*(LAMLS+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR MQN PAIRS (-|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) =-MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQLS(E11X,G11,EX2,XY2,KQ2,MQ2,NB2,1,0,0,LR)
      CALL EGQLS(E11Y,G11,EX2,XY2,KQ2,MQ2,NB2,2,0,0,LR)
      CALL EGQLS(E11Z,G11,EX2,XY2,KQ2,MQ2,NB2,3,0,0,LR)
C                                    ~
C     MULTIPLY EQLS BY PHASE TERM IF EQLS ARE NEEDED
      DO ITUV=1,NTUVLS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV,1) = PHS*E11X(M,ITUV)
          E11(M,ITUV,2) = PHS*E11Y(M,ITUV)
          E11(M,ITUV,3) = PHS*E11Z(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR MQN PAIRS (+|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) = MQN(IA1)
      MQ2(2) =-MQN(IA2)

C     GENERATE THE RAW COEFFICIENTS
      CALL EGQLS(E21X,G21,EX2,XY2,KQ2,MQ2,NB2,1,0,0,LR)
      CALL EGQLS(E21Y,G21,EX2,XY2,KQ2,MQ2,NB2,2,0,0,LR)
      CALL EGQLS(E21Z,G21,EX2,XY2,KQ2,MQ2,NB2,3,0,0,LR)
C                                    ~
C     MULTIPLY EQLS BY PHASE TERM IF EQLS ARE NEEDED
      DO ITUV=1,NTUVLS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV,1) = PHS*E21X(M,ITUV)
          E21(M,ITUV,2) = PHS*E21Y(M,ITUV)
          E21(M,ITUV,3) = PHS*E21Z(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT E22 AND E12 ARE RELATED TO THESE BY PHASE FACTORS.
C     THERE IS ALSO A RELATION WHICH ALLOWS US TO OBTAIN ESLQ FROM ELSQ.
C
      RETURN
      END
C
C
      SUBROUTINE EISLB3(E11,E21,EXL,XYZ,KQN,MQN,NBS,IPHS,IA1,IA2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          EEEEEEEE IIII  SSSSSS  LL       BBBBBBB   333333            C
C          EE        II  SS    SS LL       BB    BB 33    33           C
C          EE        II  SS       LL       BB    BB       33           C
C          EEEEEE    II   SSSSSS  LL       BBBBBBB    33333            C
C          EE        II        SS LL       BB    BB       33           C
C          EE        II  SS    SS LL       BB    BB 33    33           C
C          EEEEEEEE IIII  SSSSSS  LLLLLLLL BBBBBBB   333333            C
C                                                                      C
C -------------------------------------------------------------------- C
C  EISLB3 GENERATES A BATCH OF EISL COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB+1 -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.      C
C  THIS ACTUALLY MAKES A VECTOR LIST OF EQSL NEEDED FOR BREIT (MBPT).  C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 E11(MB2,MEQ,3),E21(MB2,MEQ,3)
      COMPLEX*16 G11(MB2,MEQ),G21(MB2,MEQ)
      COMPLEX*16 E11X(MB2,MEQ),E11Y(MB2,MEQ),E11Z(MB2,MEQ),
     &           E21X(MB2,MEQ),E21Y(MB2,MEQ),E21Z(MB2,MEQ)
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        XY2(IX,1) = XYZ(IX,IA1)
        XY2(IX,2) = XYZ(IX,IA2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KQ2(1) = KQN(IA1)
      KQ2(2) = KQN(IA2)
C
C     CALCULATE LQN VALUES
      IF(KQ2(1).LT.0) THEN
        LA =-KQ2(1)-1
      ELSE
        LA = KQ2(1)
      ENDIF
      IF(KQ2(2).LT.0) THEN
        LB =-KQ2(2)-1
      ELSE
        LB = KQ2(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NB2(1) = NBS(IA1)
      DO IBAS=1,NB2(1)
        EX2(IBAS,1) = EXL(IBAS,IA1)
      ENDDO
C
      NB2(2) = NBS(IA2)
      DO JBAS=1,NB2(2)
        EX2(JBAS,2) = EXL(JBAS,IA2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMSL  = LA+LB+1
      NTUVSL = (LAMSL+1)*(LAMSL+2)*(LAMSL+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR MQN PAIRS (-|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) =-MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQSL(E11X,G11,EX2,XY2,KQ2,MQ2,NB2,1,0,0,LR)
      CALL EGQSL(E11Y,G11,EX2,XY2,KQ2,MQ2,NB2,2,0,0,LR)
      CALL EGQSL(E11Z,G11,EX2,XY2,KQ2,MQ2,NB2,3,0,0,LR)
C                                    ~
C     MULTIPLY EQSL BY PHASE TERM IF EQSL ARE NEEDED
      DO ITUV=1,NTUVSL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV,1) = PHS*E11X(M,ITUV)
          E11(M,ITUV,2) = PHS*E11Y(M,ITUV)
          E11(M,ITUV,3) = PHS*E11Z(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR MQN PAIRS (+|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) = MQN(IA1)
      MQ2(2) =-MQN(IA2)

C     GENERATE THE RAW COEFFICIENTS
      CALL EGQSL(E21X,G21,EX2,XY2,KQ2,MQ2,NB2,1,0,0,LR)
      CALL EGQSL(E21Y,G21,EX2,XY2,KQ2,MQ2,NB2,2,0,0,LR)
      CALL EGQSL(E21Z,G21,EX2,XY2,KQ2,MQ2,NB2,3,0,0,LR)
C                                    ~
C     MULTIPLY EQSL BY PHASE TERM IF EQSL ARE NEEDED
      DO ITUV=1,NTUVSL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV,1) = PHS*E21X(M,ITUV)
          E21(M,ITUV,2) = PHS*E21Y(M,ITUV)
          E21(M,ITUV,3) = PHS*E21Z(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT E22 AND E12 ARE RELATED TO THESE BY PHASE FACTORS.
C     THERE IS ALSO A RELATION WHICH ALLOWS US TO OBTAIN ELSQ FROM ESLQ.
C
      RETURN
      END
C
C
      SUBROUTINE GQLLMK(E11,E21,G11,G21,EXL,XYZ,KQN,MQN,NBS,IPHS,
     &                                                IA1,IA2,IQ,NX,LR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      GGGGGG    QQQQQQ    LL       LL       MM       MM KK    KK      C
C     GG    GG  QQ    QQ   LL       LL       MMM     MMM KK   KK       C
C     GG       QQ      QQ  LL       LL       MMMM   MMMM KK  KK        C
C     GG       QQ      QQ  LL       LL       MM MM MM MM KKKKK         C
C     GG   GGG QQ      QQ  LL       LL       MM  MMM  MM KK  KK        C
C     GG    GG  QQ    QQ   LL       LL       MM   M   MM KK   KK       C
C      GGGGGG    QQQQQQ QQ LLLLLLLL LLLLLLLL MM       MM KK    KK      C
C                                                                      C
C -------------------------------------------------------------------- C
C  GQLLMK GENERATES A BATCH OF GQLL COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.        C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - BASIS INDICES TO CONSTRUCT EQ-BLOCK FROM.               C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  ▶ NX      - CARTESIAN PARTIAL DERIVATIVE COMPONENT.                 C
C  ▶ LR      - BASIS FUNCTION THAT PARTIAL DERIVATIVE APPLIES TO (L/R).C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C  ▶ G11     - UNIQUE GQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ G21     - UNIQUE GQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 G11(MB2,MEQ),G21(MB2,MEQ)
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,4),XY2(3,4),KQ2(4),MQ2(4),NB2(4)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        XY2(IX,1) = XYZ(IX,IA1)
        XY2(IX,2) = XYZ(IX,IA2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KQ2(1) = KQN(IA1)
      KQ2(2) = KQN(IA2)
C
C     CALCULATE LQN VALUES
      IF(KQ2(1).LT.0) THEN
        LA =-KQ2(1)-1
      ELSE
        LA = KQ2(1)
      ENDIF
      IF(KQ2(2).LT.0) THEN
        LB =-KQ2(2)-1
      ELSE
        LB = KQ2(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NB2(1) = NBS(IA1)
      DO IBAS=1,NB2(1)
        EX2(IBAS,1) = EXL(IBAS,IA1)
      ENDDO
C
      NB2(2) = NBS(IA2)
      DO JBAS=1,NB2(2)
        EX2(JBAS,2) = EXL(JBAS,IA2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMLL  = LA+LB+1
      NTUVLL = (LAMLL+1)*(LAMLL+2)*(LAMLL+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE G11, FOR MQN PAIRS (-|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) =-MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQLL(E11,G11,EX2,XY2,KQ2,MQ2,NB2,IQ,1,NX,LR)
C                                    ~
C     MULTIPLY GQLL BY PHASE TERM IF GQLL ARE NEEDED
      DO ITUV=1,NTUVLL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
          G11(M,ITUV) = PHS*G11(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE G21, FOR MQN PAIRS (+|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) = MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQLL(E21,G21,EX2,XY2,KQ2,MQ2,NB2,IQ,1,NX,LR)
C                                    ~
C     MULTIPLY GQLL BY PHASE TERM IF GQLL ARE NEEDED
      DO ITUV=1,NTUVLL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV) = PHS*E21(M,ITUV)
          G21(M,ITUV) = PHS*G21(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT G22 AND G12 ARE RELATED TO THESE BY PHASE FACTORS.
C
      RETURN
      END
C
C
      SUBROUTINE GQLSMK(E11,E21,G11,G21,EXL,XYZ,KQN,MQN,NBS,IPHS,
     &                                                IA1,IA2,IQ,NX,LR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       GGGGGG    QQQQQQ    LL       SSSSSS  MM       MM KK    KK      C
C      GG    GG  QQ    QQ   LL      SS    SS MMM     MMM KK   KK       C
C      GG       QQ      QQ  LL      SS       MMMM   MMMM KK  KK        C
C      GG       QQ      QQ  LL       SSSSSS  MM MM MM MM KKKKK         C
C      GG   GGG QQ      QQ  LL            SS MM  MMM  MM KK  KK        C
C      GG    GG  QQ    QQ   LL      SS    SS MM   M   MM KK   KK       C
C       GGGGGG    QQQQQQ QQ LLLLLLLL SSSSSS  MM       MM KK    KK      C
C                                                                      C
C -------------------------------------------------------------------- C
C  GQLSMK GENERATES A BATCH OF EQLS COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB+1 -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.      C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - BASIS INDICES TO CONSTRUCT EQ-BLOCK FROM.               C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  ▶ NX      - CARTESIAN PARTIAL DERIVATIVE COMPONENT.                 C
C  ▶ LR      - BASIS FUNCTION THAT PARTIAL DERIVATIVE APPLIES TO (L/R).C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C  ▶ G11     - UNIQUE GQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ G21     - UNIQUE GQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 G11(MB2,MEQ),G21(MB2,MEQ)
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        XY2(IX,1) = XYZ(IX,IA1)
        XY2(IX,2) = XYZ(IX,IA2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KQ2(1) = KQN(IA1)
      KQ2(2) = KQN(IA2)
C
C     CALCULATE LQN VALUES
      IF(KQ2(1).LT.0) THEN
        LA =-KQ2(1)-1
      ELSE
        LA = KQ2(1)
      ENDIF
      IF(KQ2(2).LT.0) THEN
        LB =-KQ2(2)-1
      ELSE
        LB = KQ2(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NB2(1) = NBS(IA1)
      DO IBAS=1,NB2(1)
        EX2(IBAS,1) = EXL(IBAS,IA1)
      ENDDO
C
      NB2(2) = NBS(IA2)
      DO JBAS=1,NB2(2)
        EX2(JBAS,2) = EXL(JBAS,IA2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMLS  = LA+LB+2
      NTUVLS = (LAMLS+1)*(LAMLS+2)*(LAMLS+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE G11, FOR MQN PAIRS (-|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) =-MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQLS(E11,G11,EX2,XY2,KQ2,MQ2,NB2,IQ,1,NX,LR)
C                                    ~
C     MULTIPLY GQLS BY PHASE TERM IF GQLS ARE NEEDED
      DO ITUV=1,NTUVLS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
          G11(M,ITUV) = PHS*G11(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE G21, FOR MQN PAIRS (+|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) = MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQLS(E21,G21,EX2,XY2,KQ2,MQ2,NB2,IQ,1,NX,LR)
C                                    ~
C     MULTIPLY GQLS BY PHASE TERM IF GQLS ARE NEEDED
      DO ITUV=1,NTUVLS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV) = PHS*E21(M,ITUV)
          G21(M,ITUV) = PHS*G21(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT G22 AND G12 ARE RELATED TO THESE BY PHASE FACTORS.
C     THERE IS ALSO A RELATION WHICH ALLOWS US TO OBTAIN ESLQ FROM ELSQ.
C
      RETURN
      END
C
C
      SUBROUTINE GQSLMK(E11,E21,G11,G21,EXL,XYZ,KQN,MQN,NBS,IPHS,
     &                                                IA1,IA2,IQ,NX,LR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       GGGGGG    QQQQQQ    SSSSSS  LL       MM       MM KK    KK      C
C      GG    GG  QQ    QQ  SS    SS LL       MMM     MMM KK   KK       C
C      GG       QQ      QQ SS       LL       MMMM   MMMM KK  KK        C
C      GG       QQ      QQ  SSSSSS  LL       MM MM MM MM KKKKK         C
C      GG   GGG QQ      QQ       SS LL       MM  MMM  MM KK  KK        C
C      GG    GG  QQ    QQ  SS    SS LL       MM   M   MM KK   KK       C
C       GGGGGG    QQQQQQ QQ SSSSSS  LLLLLLLL MM       MM KK    KK      C
C                                                                      C
C -------------------------------------------------------------------- C
C  GQSLMK GENERATES A BATCH OF EQSL COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB+1 -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.      C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - BASIS INDICES TO CONSTRUCT EQ-BLOCK FROM.               C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  ▶ NX      - CARTESIAN PARTIAL DERIVATIVE COMPONENT.                 C
C  ▶ LR      - BASIS FUNCTION THAT PARTIAL DERIVATIVE APPLIES TO (L/R).C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C  ▶ G11     - UNIQUE GQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ G21     - UNIQUE GQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 G11(MB2,MEQ),G21(MB2,MEQ)
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        XY2(IX,1) = XYZ(IX,IA1)
        XY2(IX,2) = XYZ(IX,IA2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KQ2(1) = KQN(IA1)
      KQ2(2) = KQN(IA2)
C
C     CALCULATE LQN VALUES
      IF(KQ2(1).LT.0) THEN
        LA =-KQ2(1)-1
      ELSE
        LA = KQ2(1)
      ENDIF
      IF(KQ2(2).LT.0) THEN
        LB =-KQ2(2)-1
      ELSE
        LB = KQ2(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NB2(1) = NBS(IA1)
      DO IBAS=1,NB2(1)
        EX2(IBAS,1) = EXL(IBAS,IA1)
      ENDDO
C
      NB2(2) = NBS(IA2)
      DO JBAS=1,NB2(2)
        EX2(JBAS,2) = EXL(JBAS,IA2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMSL  = LA+LB+2
      NTUVSL = (LAMSL+1)*(LAMSL+2)*(LAMSL+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE G11, FOR MQN PAIRS (-|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) =-MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQSL(E11,G11,EX2,XY2,KQ2,MQ2,NB2,IQ,1,NX,LR)
C                                    ~
C     MULTIPLY GQSL BY PHASE TERM IF GQSL ARE NEEDED
      DO ITUV=1,NTUVSL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
          G11(M,ITUV) = PHS*G11(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE G21, FOR MQN PAIRS (+|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) = MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQSL(E21,G21,EX2,XY2,KQ2,MQ2,NB2,IQ,1,NX,LR)
C                                    ~
C     MULTIPLY GQSL BY PHASE TERM IF GQSL ARE NEEDED
      DO ITUV=1,NTUVSL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV) = PHS*E21(M,ITUV)
          G21(M,ITUV) = PHS*G21(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT G22 AND G12 ARE RELATED TO THESE BY PHASE FACTORS.
C
      RETURN
      END
C
C
      SUBROUTINE GQSSMK(E11,E21,G11,G21,EXL,XYZ,KQN,MQN,NBS,IPHS,
     &                                                IA1,IA2,IQ,NX,LR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       GGGGGG    QQQQQQ    SSSSSS   SSSSSS  MM       MM KK    KK      C
C      GG    GG  QQ    QQ  SS    SS SS    SS MMM     MMM KK   KK       C
C      GG       QQ      QQ SS       SS       MMMM   MMMM KK  KK        C
C      GG       QQ      QQ  SSSSSS   SSSSSS  MM MM MM MM KKKKK         C
C      GG   GGG QQ      QQ       SS       SS MM  MMM  MM KK  KK        C
C      GG    GG  QQ    QQ  SS    SS SS    SS MM   M   MM KK   KK       C
C       GGGGGG    QQQQQQ QQ SSSSSS   SSSSSS  MM       MM KK    KK      C
C                                                                      C
C -------------------------------------------------------------------- C
C  GQSSMK GENERATES A BATCH OF EQSS COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB+2 -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.      C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - BASIS INDICES TO CONSTRUCT EQ-BLOCK FROM.               C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  ▶ NX      - CARTESIAN PARTIAL DERIVATIVE COMPONENT.                 C
C  ▶ LR      - BASIS FUNCTION THAT PARTIAL DERIVATIVE APPLIES TO (L/R).C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C  ▶ G11     - UNIQUE GQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ G21     - UNIQUE GQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 G11(MB2,MEQ),G21(MB2,MEQ)
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        XY2(IX,1) = XYZ(IX,IA1)
        XY2(IX,2) = XYZ(IX,IA2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KQ2(1) = KQN(IA1)
      KQ2(2) = KQN(IA2)
C
C     CALCULATE LQN VALUES
      IF(KQ2(1).LT.0) THEN
        LA =-KQ2(1)-1
      ELSE
        LA = KQ2(1)
      ENDIF
      IF(KQ2(2).LT.0) THEN
        LB =-KQ2(2)-1
      ELSE
        LB = KQ2(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NB2(1) = NBS(IA1)
      DO IBAS=1,NB2(1)
        EX2(IBAS,1) = EXL(IBAS,IA1)
      ENDDO
C
      NB2(2) = NBS(IA2)
      DO JBAS=1,NB2(2)
        EX2(JBAS,2) = EXL(JBAS,IA2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMSS  = LA+LB+3
      NTUVSS = (LAMSS+1)*(LAMSS+2)*(LAMSS+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE G11, FOR MQN PAIRS (-|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) =-MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQSS(E11,G11,EX2,XY2,KQ2,MQ2,NB2,IQ,1,NX,LR)
C                                    ~
C     MULTIPLY GQSS BY PHASE TERM IF GQSS ARE NEEDED
      DO ITUV=1,NTUVSS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
          G11(M,ITUV) = PHS*G11(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE G21, FOR MQN PAIRS (+|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) = MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EGQSS(E21,G21,EX2,XY2,KQ2,MQ2,NB2,IQ,1,NX,LR)
C                                    ~
C     MULTIPLY GQSS BY PHASE TERM IF GQSS ARE NEEDED
      DO ITUV=1,NTUVSS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV) = PHS*E21(M,ITUV)
          G21(M,ITUV) = PHS*G21(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT G22 AND G12 ARE RELATED TO THESE BY SIMPLE FACTORS.
C
      RETURN
      END
C
C
      SUBROUTINE EGQLL(ELL,GLL,EXL,XYZ,KQN,MQN,NBS,IQ,NDRV,NX,LR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            EEEEEEEE GGGGGG    QQQQQQ    LL       LL                  C
C            EE      GG    GG  QQ    QQ   LL       LL                  C
C            EE      GG       QQ      QQ  LL       LL                  C
C            EEEEEE  GG       QQ      QQ  LL       LL                  C
C            EE      GG   GGG QQ      QQ  LL       LL                  C
C            EE      GG    GG  QQ    QQ   LL       LL                  C
C            EEEEEEEE GGGGGG    QQQQQQ QQ LLLLLLLL LLLLLLLL            C
C                                                                      C
C -------------------------------------------------------------------- C
C  EGQLL GENERATES A BLOCK OF RAW EQLL/GQLL-COEFFICIENTS FOR A GIVEN   C
C  COMBINATION OF MQNS (WITH SIGN), TO BE USED BY EQLLMK/GQLLMK.       C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ XYZ  - PAIR OF BASIS FUNCTION CARTESIAN CENTRES.                  C
C  ▶ KQN  - PAIR OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.              C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ NBS  - PAIR OF BASIS FUNCTION BLOCK LENGTHS.                      C
C  ▶ IQ   - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.            C
C  ▶ NDRV - NUMBER OF DERIVATIVE OPERATORS IN OVERLAP.                 C
C  ▶ NX   - CARTESIAN PARTIAL DERIVATIVE COMPONENT.                    C
C  ▶ LR   - BASIS FUNCTION THAT PARTIAL DERIVATIVE APPLIES TO (L/R).   C
C  OUTPUT:                                                             C
C  ▶ ELL  - RAW EQLL-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C  ▶ GLL  - RAW GQLL-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 CONE
      COMPLEX*16 ELL(MB2,MEQ),GLL(MB2,MEQ)
      COMPLEX*16 ESG(MB2,MEQ),GSG(MB2,MEQ)
C
      DIMENSION LQLAB(2),MQLAB(2)
      DIMENSION EXL(MBS,2),XYZ(3,2),KQN(2),LQN(2),JQN(2),MQN(2),NBS(2)
      DIMENSION RN(MBS,2)
C
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     CLEBSCH-GORDAN SENSITIVITY PARAMETER
      DATA SENS/1.0D-14/
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     SIGN MULITPLIER FOR σ COUPLING MATRICES
      IF(IQ.EQ.0.OR.IQ.EQ.1) THEN
        SIG = 1.0D0
      ELSEIF(IQ.EQ.2.OR.IQ.EQ.3) THEN
        SIG =-1.0D0
      ENDIF
C
C     INITIALISE COEFFICIENT STORAGE ARRAY TO ZERO
      DO ITUV=1,MEQ
        DO M=1,MB2
          ELL(M,ITUV) = DCMPLX(0.0D0,0.0D0)
          GLL(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     CALCULATE LQN VALUES
      IF(KQN(1).LT.0) THEN
        LQN(1) =-KQN(1)-1
      ELSE
        LQN(1) = KQN(1)
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        LQN(2) =-KQN(2)-1
      ELSE
        LQN(2) = KQN(2)
      ENDIF
C
C     CALCULATE JQN VALUES
      JQN(1) = 2*IABS(KQN(1))-1
      JQN(2) = 2*IABS(KQN(2))-1
C
C     CALCULATE THE APPROPRIATE CLEBSCH-GORDAN FACTORS
      IF(KQN(1).LT.0) THEN
        CAU = DSQRT(DFLOAT(JQN(1)+MQN(1)  )/DFLOAT(2*JQN(1)  ))
        CAL = DSQRT(DFLOAT(JQN(1)-MQN(1)  )/DFLOAT(2*JQN(1)  ))
      ELSE
        CAU =-DSQRT(DFLOAT(JQN(1)-MQN(1)+2)/DFLOAT(2*JQN(1)+4))
        CAL = DSQRT(DFLOAT(JQN(1)+MQN(1)+2)/DFLOAT(2*JQN(1)+4))
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        CBU = DSQRT(DFLOAT(JQN(2)+MQN(2)  )/DFLOAT(2*JQN(2)  ))
        CBL = DSQRT(DFLOAT(JQN(2)-MQN(2)  )/DFLOAT(2*JQN(2)  ))
      ELSE
        CBU =-DSQRT(DFLOAT(JQN(2)-MQN(2)+2)/DFLOAT(2*JQN(2)+4))
        CBL = DSQRT(DFLOAT(JQN(2)+MQN(2)+2)/DFLOAT(2*JQN(2)+4))
      ENDIF
C
C     BASIS FUNCTION OVERLAP LIST LENGTH
      MAXM = NBS(1)*NBS(2)
C
C**********************************************************************C
C     INITIALISE COMMON GEOMETRIC INFORMATION                          C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTRES A AND B
      ABDST(1) = XYZ(1,1)-XYZ(1,2)
      ABDST(2) = XYZ(2,1)-XYZ(2,2)
      ABDST(3) = XYZ(3,1)-XYZ(3,2)
      AB2 = ABDST(1)*ABDST(1) + ABDST(2)*ABDST(2) + ABDST(3)*ABDST(3)
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
C
          EIBS(M) = EXL(IBAS,1)
          EJBS(M) = EXL(JBAS,2)
          P(M)    = EXL(IBAS,1)+EXL(JBAS,2)
          P2(M)   = 2.0D0*P(M)
          P22(M)  = P2(M)*P2(M)
          UV(M)   = EXL(IBAS,1)*EXL(JBAS,2)
          PX      = (EXL(IBAS,1)*XYZ(1,1)+EXL(JBAS,2)*XYZ(1,2))/P(M)
          PY      = (EXL(IBAS,1)*XYZ(2,1)+EXL(JBAS,2)*XYZ(2,2))/P(M)
          PZ      = (EXL(IBAS,1)*XYZ(3,1)+EXL(JBAS,2)*XYZ(3,2))/P(M)
          PAX(M)  = PX-XYZ(1,1)
          PAY(M)  = PY-XYZ(2,1)
          PAZ(M)  = PZ-XYZ(3,1)
          PBX(M)  = PX-XYZ(1,2)
          PBY(M)  = PY-XYZ(2,2)
          PBZ(M)  = PZ-XYZ(3,2)
          PA2(M)  = PAX(M)*PAX(M)+PAY(M)*PAY(M)+PAZ(M)*PAZ(M)
          PB2(M)  = PBX(M)*PBX(M)+PBY(M)*PBY(M)+PBZ(M)*PBZ(M)
          RKAB(M) = DEXP(-(EXL(IBAS,1)*EXL(JBAS,2)*AB2)/P(M))
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CASE 1: ALL KQN(1) AND KQN(2) TYPES.                             C
C**********************************************************************C
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)
      LQLAB(2) = LQN(2)
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR E/GLL0 AND E/GLLZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQLL (UPPER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLL
            DO M=1,MAXM
              DO ITUV=1,NTUV0
                ELL(M,ITUV) = ELL(M,ITUV) +     CAB*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQLL (UPPER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLL
            DO ITUV=1,NTUV0
              DO M=1,MAXM
                ELL(M,ITUV) = ELL(M,ITUV) +     CAB*ESG(M,ITUV)
                GLL(M,ITUV) = GLL(M,ITUV) +     CAB*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQLL (LOWER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLL
            DO M=1,MAXM
              DO ITUV=1,NTUV0
                ELL(M,ITUV) = ELL(M,ITUV) + SIG*CAB*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQLL (LOWER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLL
            DO ITUV=1,NTUV0
              DO M=1,MAXM
                ELL(M,ITUV) = ELL(M,ITUV) + SIG*CAB*ESG(M,ITUV)
                GLL(M,ITUV) = GLL(M,ITUV) + SIG*CAB*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR E/GLLX AND E/GLLY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQLL (UPPER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLL
            DO M=1,MAXM
              DO ITUV=1,NTUV0
                ELL(M,ITUV) = ELL(M,ITUV) + SIG*CAB*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQLL (UPPER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLL
            DO ITUV=1,NTUV0
              DO M=1,MAXM
                ELL(M,ITUV) = ELL(M,ITUV) + SIG*CAB*ESG(M,ITUV)
                GLL(M,ITUV) = GLL(M,ITUV) + SIG*CAB*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQLL (LOWER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLL
            DO M=1,MAXM
              DO ITUV=1,NTUV0
                ELL(M,ITUV) = ELL(M,ITUV) +     CAB*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQLL (LOWER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLL
            DO ITUV=1,NTUV0
              DO M=1,MAXM
                ELL(M,ITUV) = ELL(M,ITUV) +     CAB*ESG(M,ITUV)
                GLL(M,ITUV) = GLL(M,ITUV) +     CAB*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
C     INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
      LAM0  = LQLAB(1)+LQLAB(2)+NDRV
      NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C     GENERATE RNLL NORMALISATION CONSTANTS
      CALL RNLL(RN,EXL,LQN,NBS)
C
C     NORMALISE THE ELLQ/GLLQ COEFFICIENT BLOCK
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M   = M+1
          RLL = RN(IBAS,1)*RN(JBAS,2)
          DO ITUV=1,NTUV0
            ELL(M,ITUV) = RLL*ELL(M,ITUV)
            GLL(M,ITUV) = RLL*GLL(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     σ_Y SPECIAL CASE: MULTIPLY ELLY AND GLLY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV0
            ELL(M,ITUV) = CONE*ELL(M,ITUV)
            GLL(M,ITUV) = CONE*GLL(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE EGQLS(ELS,GLS,EXL,XYZ,KQN,MQN,NBS,IQ,NDRV,NX,LR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            EEEEEEEE GGGGGG    QQQQQQ    LL       SSSSSS              C
C            EE      GG    GG  QQ    QQ   LL      SS    SS             C
C            EE      GG       QQ      QQ  LL      SS                   C
C            EEEEEE  GG       QQ      QQ  LL       SSSSSS              C
C            EE      GG   GGG QQ      QQ  LL            SS             C
C            EE      GG    GG  QQ    QQ   LL      SS    SS             C
C            EEEEEEEE GGGGGG    QQQQQQ QQ LLLLLLLL SSSSSS              C
C                                                                      C
C -------------------------------------------------------------------- C
C  EGQLS GENERATES A BLOCK OF RAW EQLS/GQLS-COEFFICIENTS FOR A GIVEN   C
C  COMBINATION OF MQNS (WITH SIGN), TO BE USED BY EQLSMK/GQLSMK.       C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ XYZ  - PAIR OF BASIS FUNCTION CARTESIAN CENTRES.                  C
C  ▶ KQN  - PAIR OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.              C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ NBS  - PAIR OF BASIS FUNCTION BLOCK LENGTHS.                      C
C  ▶ IQ   - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.            C
C  ▶ NDRV - NUMBER OF DERIVATIVE OPERATORS IN OVERLAP.                 C
C  ▶ NX   - CARTESIAN PARTIAL DERIVATIVE COMPONENT.                    C
C  ▶ LR   - BASIS FUNCTION THAT PARTIAL DERIVATIVE APPLIES TO (L/R).   C
C  OUTPUT:                                                             C
C  ▶ ELS  - RAW EQLS-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C  ▶ GLS  - RAW GQLS-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 CONE
      COMPLEX*16 ELS(MB2,MEQ),GLS(MB2,MEQ)
      COMPLEX*16 ESG(MB2,MEQ),GSG(MB2,MEQ)
      COMPLEX*16 ENSG(MB2,MEQ),GNSG(MB2,MEQ)
C
      DIMENSION LQLAB(2),MQLAB(2)
      DIMENSION EXL(MBS,2),XYZ(3,2),KQN(2),LQN(2),JQN(2),MQN(2),NBS(2)
      DIMENSION RN(MBS,2)
      DIMENSION T2(MB2),T0(MB2)
C
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     CLEBSCH-GORDAN SENSITIVITY PARAMETER
      DATA SENS/1.0D-14/
C
C     SIGN MULITPLIER FOR σ COUPLING MATRICES
      IF(IQ.EQ.0.OR.IQ.EQ.1) THEN
        SIG = 1.0D0
      ELSEIF(IQ.EQ.2.OR.IQ.EQ.3) THEN
        SIG =-1.0D0
      ENDIF
C
C     INITIALISE COEFFICIENT STORAGE ARRAY TO ZERO
      DO M=1,MB2
        DO ITUV=1,MEQ
          ELS(M,ITUV) = DCMPLX(0.0D0,0.0D0)
          GLS(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     CALCULATE LQN VALUES
      IF(KQN(1).LT.0) THEN
        LQN(1) =-KQN(1)-1
      ELSE
        LQN(1) = KQN(1)
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        LQN(2) =-KQN(2)-1
      ELSE
        LQN(2) = KQN(2)
      ENDIF
C
C     CALCULATE JQN VALUES
      JQN(1) = 2*IABS(KQN(1))-1
      JQN(2) = 2*IABS(KQN(2))-1
C
C     CALCULATE THE APPROPRIATE CLEBSCH-GORDAN FACTORS
      IF(KQN(1).LT.0) THEN
        CAU = DSQRT(DFLOAT(JQN(1)+MQN(1)  )/DFLOAT(2*JQN(1)  ))
        CAL = DSQRT(DFLOAT(JQN(1)-MQN(1)  )/DFLOAT(2*JQN(1)  ))
      ELSE
        CAU =-DSQRT(DFLOAT(JQN(1)-MQN(1)+2)/DFLOAT(2*JQN(1)+4))
        CAL = DSQRT(DFLOAT(JQN(1)+MQN(1)+2)/DFLOAT(2*JQN(1)+4))
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        CBU =-DSQRT(DFLOAT(JQN(2)-MQN(2)+2)/DFLOAT(2*JQN(2)+4))
        CBL = DSQRT(DFLOAT(JQN(2)+MQN(2)+2)/DFLOAT(2*JQN(2)+4))
      ELSE
        CBU = DSQRT(DFLOAT(JQN(2)+MQN(2)  )/DFLOAT(2*JQN(2)  ))
        CBL = DSQRT(DFLOAT(JQN(2)-MQN(2)  )/DFLOAT(2*JQN(2)  ))
      ENDIF
C
C     DETERMINE THE NUMBER OF FUNCTIONS ON EACH CENTRE
      MAXM = NBS(1)*NBS(2)
C
C     KINETIC PRE-FACTORS FOR THIS BLOCK
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
          T2(M) =-2.0D0*EXL(JBAS,2)
          T0(M) = DFLOAT(2*LQN(2)+1)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INITIALISE COMMON GEOMETRIC INFORMATION                          C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTRES A AND B
      ABDST(1) = XYZ(1,1)-XYZ(1,2)
      ABDST(2) = XYZ(2,1)-XYZ(2,2)
      ABDST(3) = XYZ(3,1)-XYZ(3,2)
      AB2 = ABDST(1)*ABDST(1) + ABDST(2)*ABDST(2) + ABDST(3)*ABDST(3)
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
C
          EIBS(M) = EXL(IBAS,1)
          EJBS(M) = EXL(JBAS,2)
          P(M)    = EXL(IBAS,1)+EXL(JBAS,2)
          P2(M)   = 2.0D0*P(M)
          P22(M)  = P2(M)*P2(M)
          UV(M)   = EXL(IBAS,1)*EXL(JBAS,2)
          PX      = (EXL(IBAS,1)*XYZ(1,1)+EXL(JBAS,2)*XYZ(1,2))/P(M)
          PY      = (EXL(IBAS,1)*XYZ(2,1)+EXL(JBAS,2)*XYZ(2,2))/P(M)
          PZ      = (EXL(IBAS,1)*XYZ(3,1)+EXL(JBAS,2)*XYZ(3,2))/P(M)
          PAX(M)  = PX-XYZ(1,1)
          PAY(M)  = PY-XYZ(2,1)
          PAZ(M)  = PZ-XYZ(3,1)
          PBX(M)  = PX-XYZ(1,2)
          PBY(M)  = PY-XYZ(2,2)
          PBZ(M)  = PZ-XYZ(3,2)
          PA2(M)  = PAX(M)*PAX(M) + PAY(M)*PAY(M) + PAZ(M)*PAZ(M)
          PB2(M)  = PBX(M)*PBX(M) + PBY(M)*PBY(M) + PBZ(M)*PBZ(M)
          RKAB(M) = DEXP(-EXL(IBAS,1)*EXL(JBAS,2)*AB2/P(M))
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CASE 1: KQN(2).LT.0                                              C
C**********************************************************************C
C
      IF(KQN(2).GT.0) GOTO 100
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)
      LQLAB(2) = LQN(2)+1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR E/GLS0 AND E/GLSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQLS (UPPER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQLS (UPPER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLS
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
                GLS(M,ITUV) = GLS(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQLS (LOWER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQLS (LOWER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLS
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GLS(M,ITUV) = GLS(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR E/GLSX AND E/GLSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQLS (UPPER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQLS (UPPER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLS
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GLS(M,ITUV) = GLS(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQLS (LOWER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQLS (LOWER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLS
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
                GLS(M,ITUV) = GLS(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
      ENDIF
C
100   CONTINUE
C
C**********************************************************************C
C     CASE 2: KQN(2).GT.0                                              C
C**********************************************************************C
C
      IF(KQN(2).LT.0) GOTO 200
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)
      LQLAB(2) = LQN(2)-1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR E/GLS0 AND E/GLSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQLS (UPPER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ELS(M,ITUV) = ELS(M,ITUV) +     TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQLS (UPPER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLS
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
                GLS(M,ITUV) = GLS(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,2,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ELS(M,ITUV) = ELS(M,ITUV) +     TK*ENSG(M,ITUV)
                GLS(M,ITUV) = GLS(M,ITUV) +     TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQLS (LOWER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQLS (LOWER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLS
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GLS(M,ITUV) = GLS(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,2,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
                GLS(M,ITUV) = GLS(M,ITUV) + SIG*TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR E/GLSX AND E/GLSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQLS (UPPER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQLS (UPPER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLS
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GLS(M,ITUV) = GLS(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,2,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
                GLS(M,ITUV) = GLS(M,ITUV) + SIG*TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQLS (LOWER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ELS(M,ITUV) = ELS(M,ITUV) +     TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQLS (LOWER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLS
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
                GLS(M,ITUV) = GLS(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,2,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQLS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ELS(M,ITUV) = ELS(M,ITUV) +     TK*ENSG(M,ITUV)
                GLS(M,ITUV) = GLS(M,ITUV) +     TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF

200   CONTINUE
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
C     MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
      LAM2  = LQN(1)+LQN(2)+NDRV+2
      NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C     GENERATE RNLS NORMALISATION CONSTANTS
      CALL RNLS(RN,EXL,LQN,NBS)
C
C     NORMALISE THE ELSQ/GLSQ COEFFICIENT BLOCK AND MULTIPLY BY FACTOR i
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M   = M+1
          RLS = RN(IBAS,1)*RN(JBAS,2)
          DO ITUV=1,NTUV2
            ELS(M,ITUV) = CONE*RLS*ELS(M,ITUV)
            GLS(M,ITUV) = CONE*RLS*GLS(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     σ_Y SPECIAL CASE: MULTIPLY ELSY AND GLSY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV2
            ELS(M,ITUV) = CONE*ELS(M,ITUV)
            GLS(M,ITUV) = CONE*GLS(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE EGQSL(ESL,GSL,EXL,XYZ,KQN,MQN,NBS,IQ,NDRV,NX,LR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            EEEEEEEE GGGGGG    QQQQQQ    SSSSSS  LL                   C
C            EE      GG    GG  QQ    QQ  SS    SS LL                   C
C            EE      GG       QQ      QQ SS       LL                   C
C            EEEEEE  GG       QQ      QQ  SSSSSS  LL                   C
C            EE      GG   GGG QQ      QQ       SS LL                   C
C            EE      GG    GG  QQ    QQ  SS    SS LL                   C
C            EEEEEEEE GGGGGG    QQQQQQ QQ SSSSSS  LLLLLLLL             C
C                                                                      C
C -------------------------------------------------------------------- C
C  EGQSL GENERATES A BLOCK OF RAW EQSL/GQSL-COEFFICIENTS FOR A GIVEN   C
C  COMBINATION OF MQNS (WITH SIGN), TO BE USED BY EQSLMK/GQSLMK.       C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ XYZ  - PAIR OF BASIS FUNCTION CARTESIAN CENTRES.                  C
C  ▶ KQN  - PAIR OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.              C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ NBS  - PAIR OF BASIS FUNCTION BLOCK LENGTHS.                      C
C  ▶ IQ   - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.            C
C  ▶ NDRV - NUMBER OF DERIVATIVE OPERATORS IN OVERLAP.                 C
C  ▶ NX   - CARTESIAN PARTIAL DERIVATIVE COMPONENT.                    C
C  ▶ LR   - BASIS FUNCTION THAT PARTIAL DERIVATIVE APPLIES TO (L/R).   C
C  OUTPUT:                                                             C
C  ▶ ESL  - RAW EQSL-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C  ▶ GSL  - RAW GQSL-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 CONE
      COMPLEX*16 ESL(MB2,MEQ),GSL(MB2,MEQ)
      COMPLEX*16 ESG(MB2,MEQ),GSG(MB2,MEQ)
      COMPLEX*16 ENSG(MB2,MEQ),GNSG(MB2,MEQ)
C
      DIMENSION LQLAB(2),MQLAB(2)
      DIMENSION EXL(MBS,2),XYZ(3,2),KQN(2),LQN(2),JQN(2),MQN(2),NBS(2)
      DIMENSION RN(MBS,2)
      DIMENSION T2(MB2),T0(MB2)
C
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     CLEBSCH-GORDAN SENSITIVITY PARAMETER
      DATA SENS/1.0D-14/
C
C     SIGN MULITPLIER FOR σ COUPLING MATRICES
      IF(IQ.EQ.0.OR.IQ.EQ.1) THEN
        SIG = 1.0D0
      ELSEIF(IQ.EQ.2.OR.IQ.EQ.3) THEN
        SIG =-1.0D0
      ENDIF
C
C     INITIALISE COEFFICIENT STORAGE ARRAY TO ZERO
      DO M=1,MB2
        DO ITUV=1,MEQ
          ESL(M,ITUV) = DCMPLX(0.0D0,0.0D0)
          GSL(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     CALCULATE LQN VALUES
      IF(KQN(1).LT.0) THEN
        LQN(1) =-KQN(1)-1
      ELSE
        LQN(1) = KQN(1)
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        LQN(2) =-KQN(2)-1
      ELSE
        LQN(2) = KQN(2)
      ENDIF
C
C     CALCULATE JQN VALUES
      JQN(1) = 2*IABS(KQN(1))-1
      JQN(2) = 2*IABS(KQN(2))-1
C
C     CALCULATE THE APPROPRIATE CLEBSCH-GORDAN FACTORS
      IF(KQN(1).LT.0) THEN
        CAU =-DSQRT(DFLOAT(JQN(1)-MQN(1)+2)/DFLOAT(2*JQN(1)+4))
        CAL = DSQRT(DFLOAT(JQN(1)+MQN(1)+2)/DFLOAT(2*JQN(1)+4))
      ELSE
        CAU = DSQRT(DFLOAT(JQN(1)+MQN(1)  )/DFLOAT(2*JQN(1)  ))
        CAL = DSQRT(DFLOAT(JQN(1)-MQN(1)  )/DFLOAT(2*JQN(1)  ))
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        CBU = DSQRT(DFLOAT(JQN(2)+MQN(2)  )/DFLOAT(2*JQN(2)  ))
        CBL = DSQRT(DFLOAT(JQN(2)-MQN(2)  )/DFLOAT(2*JQN(2)  ))
      ELSE
        CBU =-DSQRT(DFLOAT(JQN(2)-MQN(2)+2)/DFLOAT(2*JQN(2)+4))
        CBL = DSQRT(DFLOAT(JQN(2)+MQN(2)+2)/DFLOAT(2*JQN(2)+4))
      ENDIF
C
C     DETERMINE THE NUMBER OF FUNCTIONS ON EACH CENTRE
      MAXM = NBS(1)*NBS(2)
C
C     KINETIC PRE-FACTORS FOR THIS BLOCK
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
          T2(M) =-2.0D0*EXL(IBAS,1)
          T0(M) = DFLOAT(2*LQN(1)+1)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INITIALISE COMMON GEOMETRIC INFORMATION                          C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTRES A AND B
      ABDST(1) = XYZ(1,1)-XYZ(1,2)
      ABDST(2) = XYZ(2,1)-XYZ(2,2)
      ABDST(3) = XYZ(3,1)-XYZ(3,2)
      AB2 = ABDST(1)*ABDST(1) + ABDST(2)*ABDST(2) + ABDST(3)*ABDST(3)
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
C
          EIBS(M) = EXL(IBAS,1)
          EJBS(M) = EXL(JBAS,2)
          P(M)    = EXL(IBAS,1)+EXL(JBAS,2)
          P2(M)   = 2.0D0*P(M)
          P22(M)  = P2(M)*P2(M)
          UV(M)   = EXL(IBAS,1)*EXL(JBAS,2)
          PX      = (EXL(IBAS,1)*XYZ(1,1)+EXL(JBAS,2)*XYZ(1,2))/P(M)
          PY      = (EXL(IBAS,1)*XYZ(2,1)+EXL(JBAS,2)*XYZ(2,2))/P(M)
          PZ      = (EXL(IBAS,1)*XYZ(3,1)+EXL(JBAS,2)*XYZ(3,2))/P(M)
          PAX(M)  = PX-XYZ(1,1)
          PAY(M)  = PY-XYZ(2,1)
          PAZ(M)  = PZ-XYZ(3,1)
          PBX(M)  = PX-XYZ(1,2)
          PBY(M)  = PY-XYZ(2,2)
          PBZ(M)  = PZ-XYZ(3,2)
          PA2(M)  = PAX(M)*PAX(M) + PAY(M)*PAY(M) + PAZ(M)*PAZ(M)
          PB2(M)  = PBX(M)*PBX(M) + PBY(M)*PBY(M) + PBZ(M)*PBZ(M)
          RKAB(M) = DEXP(-(EXL(IBAS,1)*EXL(JBAS,2)*AB2)/P(M))
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CASE 1: KQN(2).LT.0                                              C
C**********************************************************************C
C
      IF(KQN(1).GT.0) GOTO 100
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)+1
      LQLAB(2) = LQN(2)
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR E/GSL0 AND E/GSLZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQSL (UPPER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSL (UPPER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSL
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
                GSL(M,ITUV) = GSL(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQSL (LOWER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSL (LOWER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSL
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GSL(M,ITUV) = GSL(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR E/GSLX AND E/GSLY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQSL (UPPER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSL (UPPER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSL
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GSL(M,ITUV) = GSL(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQSL (LOWER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSL (LOWER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSL
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
                GSL(M,ITUV) = GSL(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF
C
100   CONTINUE
C
C**********************************************************************C
C     CASE 2: KQN(1).GT.0                                              C
C**********************************************************************C
C
      IF(KQN(1).LT.0) GOTO 200
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)-1
      LQLAB(2) = LQN(2)
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR E/GSL0 AND E/GSLZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQSL (UPPER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ESL(M,ITUV) = ESL(M,ITUV) +     TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSL (UPPER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSL
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
                GSL(M,ITUV) = GSL(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSL [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ESL(M,ITUV) = ESL(M,ITUV) +     TK*ENSG(M,ITUV)
                GSL(M,ITUV) = GSL(M,ITUV) +     TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQSL (LOWER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSL (LOWER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSL
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GSL(M,ITUV) = GSL(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSL [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ENSG(M,ITUV)
                GSL(M,ITUV) = GSL(M,ITUV) + SIG*TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR E/GSLX AND E/GSLY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQSL (UPPER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSL (UPPER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSL
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GSL(M,ITUV) = GSL(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSL [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ENSG(M,ITUV)
                GSL(M,ITUV) = GSL(M,ITUV) + SIG*TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQSL (LOWER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ESL(M,ITUV) = ESL(M,ITUV) +     TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSL (LOWER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSL
            DO M=1,MAXM
              TK = CAB*T0(M)
              DO ITUV=1,NTUV0
                ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
                GSL(M,ITUV) = GSL(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSL [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T2(M)
              DO ITUV=1,NTUV2
                ESL(M,ITUV) = ESL(M,ITUV) +     TK*ENSG(M,ITUV)
                GSL(M,ITUV) = GSL(M,ITUV) +     TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF

200   CONTINUE
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
C     MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
      LAM2  = LQN(1)+LQN(2)+NDRV+2
      NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C     GENERATE RNSL NORMALISATION CONSTANTS
      CALL RNSL(RN,EXL,LQN,NBS)
C
C     NORMALISE THE ESLQ/GSLQ COEFFICIENT BLOCK AND MULTIPLY BY FACTOR -i
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M   = M+1
          RSL = RN(IBAS,1)*RN(JBAS,2)
          DO ITUV=1,NTUV2
            ESL(M,ITUV) =-CONE*RSL*ESL(M,ITUV)
            GSL(M,ITUV) =-CONE*RSL*GSL(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     σ_Y SPECIAL CASE: MULTIPLY ESLY AND GSLY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV2
            ESL(M,ITUV) = CONE*ESL(M,ITUV)
            GSL(M,ITUV) = CONE*GSL(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE EGQSS(ESS,GSS,EXL,XYZ,KQN,MQN,NBS,IQ,NDRV,NX,LR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            EEEEEEEE GGGGGG    QQQQQQ    SSSSSS   SSSSSS              C
C            EE      GG    GG  QQ    QQ  SS    SS SS    SS             C
C            EE      GG       QQ      QQ SS       SS                   C
C            EEEEEE  GG       QQ      QQ  SSSSSS   SSSSSS              C
C            EE      GG   GGG QQ      QQ       SS       SS             C
C            EE      GG    GG  QQ    QQ  SS    SS SS    SS             C
C            EEEEEEEE GGGGGG    QQQQQQ QQ SSSSSS   SSSSSS              C
C                                                                      C
C -------------------------------------------------------------------- C
C  EGQSS GENERATES A BLOCK OF RAW EQSS/GQSS-COEFFICIENTS FOR A GIVEN   C
C  COMBINATION OF MQNS (WITH SIGN), TO BE USED BY EQSSMK/GQSSMK.       C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ XYZ  - PAIR OF BASIS FUNCTION CARTESIAN CENTRES.                  C
C  ▶ KQN  - PAIR OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.              C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ NBS  - PAIR OF BASIS FUNCTION BLOCK LENGTHS.                      C
C  ▶ IQ   - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.            C
C  ▶ NDRV - NUMBER OF DERIVATIVE OPERATORS IN OVERLAP.                 C
C  ▶ NX   - CARTESIAN PARTIAL DERIVATIVE COMPONENT.                    C
C  ▶ LR   - BASIS FUNCTION THAT PARTIAL DERIVATIVE APPLIES TO (L/R).   C
C  OUTPUT:                                                             C
C  ▶ ESS  - RAW EQSS-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C  ▶ GSS  - RAW GQSS-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 CONE
      COMPLEX*16 ESS(MB2,MEQ),GSS(MB2,MEQ)
      COMPLEX*16 ESG(MB2,MEQ),GSG(MB2,MEQ)
      COMPLEX*16 ENSG(MB2,MEQ),GNSG(MB2,MEQ)
C
      DIMENSION LQLAB(2),MQLAB(2)
      DIMENSION EXL(MBS,2),XYZ(3,2),KQN(2),LQN(2),JQN(2),MQN(2),NBS(2)
      DIMENSION RN(MBS,2)
      DIMENSION T22(MB2),T20(MB2),T02(MB2),T00(MB2)
C
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     CLEBSCH-GORDAN SENSITIVITY PARAMETER
      DATA SENS/1.0D-14/
C
C     SIGN MULITPLIER FOR σ COUPLING MATRICES
      IF(IQ.EQ.0.OR.IQ.EQ.1) THEN
        SIG = 1.0D0
      ELSEIF(IQ.EQ.2.OR.IQ.EQ.3) THEN
        SIG =-1.0D0
      ENDIF
C
C     INITIALISE COEFFICIENT STORAGE ARRAY TO ZERO
      DO M=1,MB2
        DO ITUV=1,MEQ
          ESS(M,ITUV) = DCMPLX(0.0D0,0.0D0)
          GSS(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     CALCULATE LQN VALUES
      IF(KQN(1).LT.0) THEN
        LQN(1) =-KQN(1)-1
      ELSE
        LQN(1) = KQN(1)
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        LQN(2) =-KQN(2)-1
      ELSE
        LQN(2) = KQN(2)
      ENDIF
C
C     CALCULATE JQN VALUES
      JQN(1) = 2*IABS(KQN(1))-1
      JQN(2) = 2*IABS(KQN(2))-1
C
C     CALCULATE THE APPROPRIATE CLEBSCH-GORDAN FACTORS
      IF(KQN(1).LT.0) THEN
        CAU =-DSQRT(DFLOAT(JQN(1)-MQN(1)+2)/DFLOAT(2*JQN(1)+4))
        CAL = DSQRT(DFLOAT(JQN(1)+MQN(1)+2)/DFLOAT(2*JQN(1)+4))
      ELSE
        CAU = DSQRT(DFLOAT(JQN(1)+MQN(1)  )/DFLOAT(2*JQN(1)  ))
        CAL = DSQRT(DFLOAT(JQN(1)-MQN(1)  )/DFLOAT(2*JQN(1)  ))
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        CBU =-DSQRT(DFLOAT(JQN(2)-MQN(2)+2)/DFLOAT(2*JQN(2)+4))
        CBL = DSQRT(DFLOAT(JQN(2)+MQN(2)+2)/DFLOAT(2*JQN(2)+4))
      ELSE
        CBU = DSQRT(DFLOAT(JQN(2)+MQN(2)  )/DFLOAT(2*JQN(2)  ))
        CBL = DSQRT(DFLOAT(JQN(2)-MQN(2)  )/DFLOAT(2*JQN(2)  ))
      ENDIF
C
C     DETERMINE THE NUMBER OF FUNCTIONS ON EACH CENTRE
      MAXM = NBS(1)*NBS(2)
C
C     KINETIC PRE-FACTORS FOR THIS BLOCK
      RL1 = DFLOAT(2*LQN(1)+1)
      RL2 = DFLOAT(2*LQN(2)+1)
C
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
          T22(M) = 4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
          T20(M) =-2.0D0*RL2*EXL(IBAS,1)
          T02(M) =-2.0D0*RL1*EXL(JBAS,2)
          T00(M) = RL1*RL2
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INITIALISE COMMON GEOMETRIC INFORMATION                          C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTRES A AND B
      ABDST(1) = XYZ(1,1)-XYZ(1,2)
      ABDST(2) = XYZ(2,1)-XYZ(2,2)
      ABDST(3) = XYZ(3,1)-XYZ(3,2)
      AB2 = ABDST(1)*ABDST(1) + ABDST(2)*ABDST(2) + ABDST(3)*ABDST(3)
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
C
          EIBS(M) = EXL(IBAS,1)
          EJBS(M) = EXL(JBAS,2)
          P(M)    = EXL(IBAS,1)+EXL(JBAS,2)
          P2(M)   = 2.0D0*P(M)
          P22(M)  = P2(M)*P2(M)
          UV(M)   = EXL(IBAS,1)*EXL(JBAS,2)
          PX      = (EXL(IBAS,1)*XYZ(1,1)+EXL(JBAS,2)*XYZ(1,2))/P(M)
          PY      = (EXL(IBAS,1)*XYZ(2,1)+EXL(JBAS,2)*XYZ(2,2))/P(M)
          PZ      = (EXL(IBAS,1)*XYZ(3,1)+EXL(JBAS,2)*XYZ(3,2))/P(M)
          PAX(M)  = PX-XYZ(1,1)
          PAY(M)  = PY-XYZ(2,1)
          PAZ(M)  = PZ-XYZ(3,1)
          PBX(M)  = PX-XYZ(1,2)
          PBY(M)  = PY-XYZ(2,2)
          PBZ(M)  = PZ-XYZ(3,2)
          PA2(M)  = PAX(M)*PAX(M)+PAY(M)*PAY(M)+PAZ(M)*PAZ(M)
          PB2(M)  = PBX(M)*PBX(M)+PBY(M)*PBY(M)+PBZ(M)*PBZ(M)
          RKAB(M) = DEXP(-(EXL(IBAS,1)*EXL(JBAS,2)*AB2)/P(M))
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CASE 1: KQN(1).LT.0 AND KQN(2).LT.0                              C
C**********************************************************************C
C
      IF(KQN(1).GT.0.OR.KQN(2).GT.0) GOTO 100
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)+1
      LQLAB(2) = LQN(2)+1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR E/GSS0 AND E/GSSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQSS (UPPER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (UPPER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQSS (LOWER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (LOWER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR E/GSSX AND E/GSSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQSS (UPPER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (UPPER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         EQSS (LOWER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (LOWER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF
C
100   CONTINUE
C
C**********************************************************************C
C     CASE 2: KQN(1).LT.0 AND KQN(2).GT.0                              C
C**********************************************************************C
C
      IF(KQN(1).GT.0.OR.KQN(2).LT.0) GOTO 200
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)+1
      LQLAB(2) = LQN(2)-1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR E/GSS0 AND E/GSSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQSS (UPPER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +    TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (UPPER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,2,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +    TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +    TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQSS (LOWER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (LOWER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,2,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR E/GSSX AND E/GSSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQSS (UPPER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (UPPER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,2,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQSS (LOWER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (LOWER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,2,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF
C
200   CONTINUE
C
C**********************************************************************C
C     CASE 3: KQN(1).GT.0 AND KQN(2).LT.0                              C
C**********************************************************************C
C
      IF(KQN(1).LT.0.OR.KQN(2).GT.0) GOTO 300
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)-1
      LQLAB(2) = LQN(2)+1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR E/GSS0 AND E/GSSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQSS (UPPER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (UPPER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQSS (LOWER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (LOWER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR E/GSSX AND E/GSSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQSS (UPPER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (UPPER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         EQSS (LOWER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (LOWER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF

300   CONTINUE      
C
C**********************************************************************C
C     CASE 4: KQN(1).GT.0 AND KQN(2).GT.0                              C
C**********************************************************************C
C
      IF(KQN(1).LT.0.OR.KQN(2).LT.0) GOTO 400
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)-1
      LQLAB(2) = LQN(2)-1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR E/GSS0 AND E/GSSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          LAM4  = LQLAB(1)+LQLAB(2)+NDRV+4
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
          NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C         EQSS (UPPER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T00(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ESG)
            CALL ESTEPN(ENSG,ESG,LAM2,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV4
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (UPPER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T00(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO GQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,2,NX,LR)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO GQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GSG)
            CALL GSTEPN(ENSG,ESG,GNSG,GSG,LAM2,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO GQSS [NA=1,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV4
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          LAM4  = LQLAB(1)+LQLAB(2)+NDRV+4
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
          NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C         EQSS (LOWER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T00(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ESG)
            CALL ESTEPN(ENSG,ESG,LAM2,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV4
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (LOWER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T00(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,2,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GSG)
            CALL GSTEPN(ENSG,ESG,GNSG,GSG,LAM2,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=1,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV4
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR E/GSSX AND E/GSSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          LAM4  = LQLAB(1)+LQLAB(2)+NDRV+4
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
          NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C         EQSS (UPPER-LOWER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T00(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ESG)
            CALL ESTEPN(ENSG,ESG,LAM2,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV4
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (UPPER-LOWER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T00(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,2,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GSG)
            CALL GSTEPN(ENSG,ESG,GNSG,GSG,LAM2,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=1,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV4
                ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) + SIG*TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)+NDRV
          LAM2  = LQLAB(1)+LQLAB(2)+NDRV+2
          LAM4  = LQLAB(1)+LQLAB(2)+NDRV+4
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
          NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C         EQSS (LOWER-UPPER) COEFFICIENTS
          IF(NDRV.EQ.0) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
            CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
            DO M=1,MAXM
              TK = CAB*T00(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
            CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ESG)
            CALL ESTEPN(ENSG,ESG,LAM2,MAXM,1)
C
C           ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV4
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
              ENDDO
            ENDDO
C
C         GQSS (LOWER-UPPER) COEFFICIENTS
          ELSEIF(NDRV.EQ.1) THEN
C
C           GENERATE RAW SPHERICAL GAUSSIAN (GS) COEFFICIENTS
            CALL GSGTF(ESG,GSG,LQLAB,MQLAB,MAXM,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS
            DO M=1,MAXM
              TK = CAB*T00(M)
              DO ITUV=1,NTUV0
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=1,NB=0]
            DO M=1,MAXM
              TK = CAB*T20(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO GNSG)
            CALL GSTEPN(ESG,ENSG,GSG,GNSG,LAM0,MAXM,2,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=0,NB=1]
            DO M=1,MAXM
              TK = CAB*T02(M)
              DO ITUV=1,NTUV2
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GNSG(M,ITUV)
              ENDDO
            ENDDO
C
C           INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO GSG)
            CALL GSTEPN(ENSG,ESG,GNSG,GSG,LAM2,MAXM,1,NX,LR)
C
C           ADD THIS ANGULAR-(GS) CONTRIBUTION TO GQSS [NA=1,NB=1]
            DO M=1,MAXM
              TK = CAB*T22(M)
              DO ITUV=1,NTUV4
                ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
                GSS(M,ITUV) = GSS(M,ITUV) +     TK*GSG(M,ITUV)
              ENDDO
            ENDDO
C
          ENDIF
C
        ENDIF
      ENDIF
C
400   CONTINUE
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
C     GENERATE RNSS NORMALISATION CONSTANTS
      CALL RNSS(RN,EXL,LQN,NBS)
C
C     MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
      LAM4  = LQN(1)+LQN(2)+NDRV+4
      NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C     NORMALISE THE ESSQ/GSSQ COEFFICIENT BLOCK
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M   = M+1
          RSS = RN(IBAS,1)*RN(JBAS,2)
          DO ITUV=1,NTUV4
            ESS(M,ITUV) = RSS*ESS(M,ITUV)
            GSS(M,ITUV) = RSS*GSS(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     σ_Y SPECIAL CASE: MULTIPLY ESSY AND GSSY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV4
            ESS(M,ITUV) = CONE*ESS(M,ITUV)
            GSS(M,ITUV) = CONE*GSS(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE ESGTF(ESG,LQN,MQN,MAXM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C              EEEEEEEE SSSSSS   GGGGGG TTTTTTTT FFFFFFFF              C
C              EE      SS    SS GG    GG   TT    FF                    C
C              EE      SS       GG         TT    FF                    C
C              EEEEEE   SSSSSS  GG         TT    FFFFFF                C
C              EE            SS GG   GGG   TT    FF                    C
C              EE      SS    SS GG    GG   TT    FF                    C
C              EEEEEEEE SSSSSS   GGGGGG    TT    FF                    C
C                                                                      C
C -------------------------------------------------------------------- C
C  ESGTF CONSTRUCTS THE EXPANSION COEFFICIENTS OF THE OVERLAP DENSITY  C
C  OF TWO SPHERICAL HARMONIC FUNCTIONS IN AN AUXILIARY HGTF BASIS.     C
C                                                                      C
C  THE OVERLAP DENSITY IS DEFINED BY Y*[L,M]Y[L',M'], WHERE Y[L,M] ARE C
C  SPHERICAL HARMONICS FOLLOWING THE CONDON-SHORTLEY PHASE CONVENTION. C
C                                                                      C
C  THE REQUIRED COEFFICIENTS ARE GENERATED BY A CALL TO EVRS, WHICH IS C
C  CONSTRUCTED ACCORDING TO THE RECURRENCE RELATIONS DEFINED BY        C
C  V.R.SAUNDERS. THE OUTPUT OF EVRS IS THEN ADJUSTED TO INCLUDE THE    C
C  ANGULAR NORMALISATION CONSTANTS, AS WELL AS A PHASE FACTOR TO       C
C  CONVERT FROM THE SCHIFF TO THE CONDON-SHORTLEY PHASE CONVENTION.    C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ LQN  - PAIR OF BASIS SET ORBITAL QUANTUM NUMBERS.                 C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ MAXM - SIZE OF THIS BLOCK.                                        C
C  OUTPUT:                                                             C
C  ▶ ESG  - EXPANSION COEFFICIENTS FOR EACH OVERLAP IN THE BLOCK.      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION LQN(2),MQN(2),MQNLAB(2)
C
      COMPLEX*16 ESG(MB2,MEQ)
C
      COMMON/FCTS/RFACT(0:80),SFACT(0:80)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
      LAM  = LQN(1)+LQN(2)
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     EVRS IS CALLED WITH THE SIGN OF MQN(1) REVERSED
C     TO AFFECT COMPLEX CONJUGATION, ALONG WITH THE REQUISITE
C     PHASE, WHICH IS CALCULATED LATER.
      MQNLAB(1) =-MQN(1)
      MQNLAB(2) = MQN(2)
C
C     GENERATE RAW COEFFICIENTS WITH EVRS
      CALL EVRS(ESG,LQN,MQNLAB,MAXM)
C
C     MQN PHASES
      PHS1 = (-1.0D0)**((MQN(1)+IABS(MQN(1)))/2)
      PHS2 = (-1.0D0)**((MQN(2)+IABS(MQN(2)))/2)
C
C     CG COEFFICIENTS
      PI4 = 0.25D0/PI
      DGL = DFLOAT((2*LQN(1)+1)*(2*LQN(2)+1))
      CG1 = RFACT(LQN(1)-IABS(MQN(1)))/RFACT(LQN(1)+IABS(MQN(1)))
      CG2 = RFACT(LQN(2)-IABS(MQN(2)))/RFACT(LQN(2)+IABS(MQN(2)))
C
C     ANGULAR NORMALISATION CONSTANT
      ANG = PHS1*PHS2*PI4*DSQRT(DGL*CG1*CG2)
C
C     APPLY ANGULAR FACTOR TO RAW COEFFICIENTS
      DO M=1,MAXM
        DO ITUV=1,NTUV
          ESG(M,ITUV) = ANG*ESG(M,ITUV)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE GSGTF(ESG,GSG,LQN,MQN,MAXM,NX,LR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C              GGGGGG   SSSSSS   GGGGGG TTTTTTTT FFFFFFFF              C
C             GG    GG SS    SS GG    GG   TT    FF                    C
C             GG       SS       GG         TT    FF                    C
C             GG        SSSSSS  GG         TT    FFFFFF                C
C             GG   GGG       SS GG   GGG   TT    FF                    C
C             GG    GG SS    SS GG    GG   TT    FF                    C
C              GGGGGG   SSSSSS   GGGGGG    TT    FF                    C
C                                                                      C
C -------------------------------------------------------------------- C
C  GSGTF CONSTRUCTS THE EXPANSION COEFFICIENTS OF THE OVERLAP DENSITY  C
C  OF TWO SPHERICAL HARMONIC FUNCTIONS IN AN AUXILIARY HGTF BASIS.     C
C                                                                      C
C  THE OVERLAP DENSITY IS DEFINED BY Y*[L,M]Y[L',M'], WHERE Y[L,M] ARE C
C  SPHERICAL HARMONICS FOLLOWING THE CONDON-SHORTLEY PHASE CONVENTION. C
C                                                                      C
C  THE REQUIRED COEFFICIENTS ARE GENERATED BY A CALL TO GVRS, WHICH IS C
C  CONSTRUCTED ACCORDING TO THE RECURRENCE RELATIONS DEFINED BY        C
C  V.R.SAUNDERS. THE OUTPUT OF GVRS IS THEN ADJUSTED TO INCLUDE THE    C
C  ANGULAR NORMALISATION CONSTANTS, AS WELL AS A PHASE FACTOR TO       C
C  CONVERT FROM THE SCHIFF TO THE CONDON-SHORTLEY PHASE CONVENTION.    C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ LQN  - PAIR OF BASIS SET ORBITAL QUANTUM NUMBERS.                 C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ MAXM - SIZE OF THIS BLOCK.                                        C
C  ▶ NX   - CARTESIAN DERIVATIVE DIRECTION.                            C
C  ▶ LR   - BASIS FUNCTION THAT PARTIAL DERIVATIVE APPLIES TO (L/R).   C
C  OUTPUT:                                                             C
C  ▶ ESG  - EXPANSION COEFFICIENTS FOR EACH OVERLAP IN THE BLOCK.      C
C  ▶ GSG  - DERIVATIVE COEFFICIENTS FOR EACH OVERLAP IN THE BLOCK.     C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      DIMENSION LQN(2),MQN(2),MQNLAB(2)
C
      COMPLEX*16 ESG(MB2,MEQ),GSG(MB2,MEQ)
C
      COMMON/FCTS/RFACT(0:80),SFACT(0:80)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
      NDRV = 1
      LAM  = LQN(1)+LQN(2)+NDRV
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GVRS IS CALLED WITH THE SIGN OF MQN(1) REVERSED
C     TO AFFECT COMPLEX CONJUGATION, ALONG WITH THE REQUISITE
C     PHASE, WHICH IS CALCULATED LATER.
      MQNLAB(1) =-MQN(1)
      MQNLAB(2) = MQN(2)
C
C     GENERATE RAW COEFFICIENTS WITH GVRS
      CALL GVRS(ESG,GSG,LQN,MQNLAB,MAXM,NX,LR)
C
C     MQN PHASES
      PHS1 = (-1.0D0)**(MQN(1)+IABS(MQN(1)))
      PHS2 = (-1.0D0)**(MQN(2)+IABS(MQN(2)))
C
C     CG COEFFICIENTS
      PI4 = 0.25D0/PI
      DGL = DFLOAT((2*LQN(1)+1)*(2*LQN(2)+1))
      CG1 = RFACT(LQN(1)-IABS(MQN(1)))/RFACT(LQN(1)+IABS(MQN(1)))
      CG2 = RFACT(LQN(2)-IABS(MQN(2)))/RFACT(LQN(2)+IABS(MQN(2)))
C
C     ANGULAR NORMALISATION CONSTANT      
      ANG = PHS1*PHS2*PI4*DSQRT(DGL*CG1*CG2)
C
C     APPLY MULTIPLICATIVE CONSTANT TO RAW ANGULAR COEFFICIENTS
      DO M=1,MAXM
        DO ITUV=1,NTUV
          ESG(M,ITUV) = ANG*ESG(M,ITUV)
          GSG(M,ITUV) = ANG*GSG(M,ITUV)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE EVRS(ESG,LQN,MQN,MAXM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 EEEEEEEE VV    VV RRRRRRR   SSSSSS                   C
C                 EE       VV    VV RR    RR SS    SS                  C
C                 EE       VV    VV RR    RR SS                        C
C                 EEEEEE   VV    VV RR    RR  SSSSSS                   C
C                 EE        VV  VV  RRRRRRR        SS                  C
C                 EE         VVVV   RR    RR SS    SS                  C
C                 EEEEEEEE    VV    RR    RR  SSSSSS                   C
C                                                                      C
C -------------------------------------------------------------------- C
C  EVRS EVALUATES THE EXPANSION COEFFICIENTS OF THE OVERLAP CHARGE     C
C  DENSITY OF SGTFS IN AN AUXILIARY HGTF. COEFFICIENTS ARE EVALUATED   C
C  USING THE RECURRENCE RELATIONS DEFINED BY VIC SAUNDERS IN:          C
C                                                                      C
C  V.R.SAUNDERS, "MOLECULAR INTEGRALS FOR GAUSSIAN-TYPE FUNCTIONS",    C
C  METHODS OF COMPUTATIONAL MOLECULAR PHYSICS, DIERCKSEN AND WILSON,   C
C  pp 1-26, REIDEL PUBLISHING, DORDRECHT (1983).                       C
C                                                                      C
C  THE EQ-COEFFS IN THIS PROCEDURE ARE FOR AN UN-NORMALISED SGTF.      C
C  THE COEFFICIENTS ARE DETERMINED ACCORDING TO THE SAME RULES AS      C
C  DEFINED IN THE ABOVE ARTICLE. CONSEQUENTLY, IT SHOULD BE NOTED THAT C
C  THE COEFFICIENTS ARE THOSE OF SPHERICAL HARMONIC FUNCTIONS THAT ARE C
C   ▶ UN-NORMALISED                                                    C
C   ▶ SATISFY THE SCHIFF PHASE CONVENTION.                             C
C                                                                      C
C  THE OUTLINE FOR THE GENERATION OF EQ-COEFFICIENTS IS TAKEN FROM p16 C
C  OF THE ABOVE ARTICLE. EQUATION NUMBERS ARE GIVEN IN COMMENTS.       C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ LQN  - PAIR OF BASIS SET ORBITAL QUANTUM NUMBERS.                 C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ MAXM - SIZE OF THIS BLOCK.                                        C
C  OUTPUT:                                                             C
C  ▶ ESG  - UN-NORMALISED EXPANSION COEFFICIENTS FOR THIS BLOCK.       C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMPLEX*16 ESG(MB2,MEQ),ETEMP(MB2*MRC)
C
      DIMENSION NBAS(2),LQN(2),MQN(2)
C
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     IMPORT LQN AND BASIS PAIR MQNS FOR LOCAL USE
      LQNA = LQN(1)
      LQNB = LQN(2)
      MQNA = MQN(1)
      MQNB = MQN(2)
C
C     INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
      LMAX = LQN(1)+LQN(2)
      NTUV = (LMAX+1)*(LMAX+2)*(LMAX+3)/6
C
C     INITIALISE ESG AND ARRAY
      DO M=1,MAXM
        DO ITUV=1,NTUV
          ESG(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     EXIT IF LQN,MQN ORDERS YIELD ZERO CG-COEFFICIENTS
      IF(IABS(MQN(1)).GT.LQN(1).OR.IABS(MQN(2)).GT.LQN(2)) RETURN
C
C     CHECK THAT LMAX IS WITHIN THE BOUNDS OF MKP
      IF(LMAX.GT.MKP+1) THEN
        WRITE(6, *) 'In EVRS: LMAX exceeds MKP+1 parameter.',LMAX,MKP+1
        WRITE(7, *) 'In EVRS: LMAX exceeds MKP+1 parameter.',LMAX,MKP+1
        STOP
      ENDIF
C
c     INITIALISE TEMPORARY ARRAY
      DO M=1,MB2*MRC
        ETEMP(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     SET INITIAL VALUES TO E[0,0;0,0;0,0,0,0] = RKAB
      DO M=1,MAXM
        ETEMP(M) = DCMPLX(RKAB(M),0.0D0)
      ENDDO
C
C     STEP 1:
C     GENERATE E[|MQNA|,MQNA;0,0] FROM E[0,0;0,0] USING
C     SIMULTANEOUS STEP OF LQN AND MQN ON CENTRE A
      ISTART = 0
      LAM    = 0
      CALL ESTEPLM(ETEMP,LAM,ISTART,MQNA,MAXM,1)
C
C     STEP 2:
C     GENERATE E[LQNA,MQNA;0,0] FROM E[|MQNA|,MQNA;0,0]
C     USING THE STEP OF LQN ONLY ON CENTRE A
      CALL ESTEPL(ETEMP,LAM,ISTART,LQNA,MQNA,MAXM,1)
C
C     STEP 3:
C     GENERATE E[LQNA,MQNA;|MQNB|,MQNB] FROM E[LQNA,MQNA;0,0]
C     USING SIMULTANEOUS STEP OF LQN AND MQN ON CENTRE B
      CALL ESTEPLM(ETEMP,LAM,ISTART,MQNB,MAXM,2)
C
C     STEP 4:
C     GENERATE E[LQNA,MQNA;LQNB,MQNB] FROM E[LQNA,MQNA;|MQNB|,MQNB]
C     USING THE STEP OF LQN ONLY ON CENTRE B
      CALL ESTEPL(ETEMP,LAM,ISTART,LQNB,MQNB,MAXM,2)
C
C     NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS LAM
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     STEP 5:
C     COPY FINAL BLOCK OF ENTRIES AS THE REQUIRED OUTPUT
      K = 0
      DO ITUV=1,NTUV
        DO M=1,MAXM
          K = K+1
          ESG(M,ITUV) = ETEMP(ISTART+K)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE GVRS(ESG,GSG,LQN,MQN,MAXM,NX,LR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                  GGGGGG  VV    VV RRRRRRR   SSSSSS                   C
C                 GG    GG VV    VV RR    RR SS    SS                  C
C                 GG       VV    VV RR    RR SS                        C
C                 GG       VV    VV RR    RR  SSSSSS                   C
C                 GG   GGG  VV  VV  RRRRRRR        SS                  C
C                 GG    GG   VVVV   RR    RR SS    SS                  C
C                  GGGGGG     VV    RR    RR  SSSSSS                   C
C                                                                      C
C -------------------------------------------------------------------- C
C  GVRS EVALUATES THE EXPANSION COEFFICIENTS OF THE OVERLAP CHARGE     C
C  DENSITY OF SGTFS IN AN AUXILIARY HGTF. COEFFICIENTS ARE EVALUATED   C
C  USING THE RECURRENCE RELATIONS DEFINED BY VIC SAUNDERS IN:          C
C                                                                      C
C  V.R.SAUNDERS, "MOLECULAR INTEGRALS FOR GAUSSIAN-TYPE FUNCTIONS",    C
C  METHODS OF COMPUTATIONAL MOLECULAR PHYSICS, DIERCKSEN AND WILSON,   C
C  pp 1-26, REIDEL PUBLISHING, DORDRECHT (1983).                       C
C                                                                      C
C  THE EQ-COEFFS IN THIS PROCEDURE ARE FOR AN UN-NORMALISED SGTF.      C
C  THE COEFFICIENTS ARE DETERMINED ACCORDING TO THE SAME RULES AS      C
C  DEFINED IN THE ABOVE ARTICLE. CONSEQUENTLY, IT SHOULD BE NOTED THAT C
C  THE COEFFICIENTS ARE THOSE OF SPHERICAL HARMONIC FUNCTIONS THAT ARE C
C   ▶ UN-NORMALISED                                                    C
C   ▶ SATISFY THE SCHIFF PHASE CONVENTION.                             C
C                                                                      C
C  THE OUTLINE FOR THE GENERATION OF EQ-COEFFICIENTS IS TAKEN FROM p16 C
C  OF THE ABOVE ARTICLE. EQUATION NUMBERS ARE GIVEN IN COMMENTS.       C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ LQN  - PAIR OF BASIS SET ORBITAL QUANTUM NUMBERS.                 C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ MAXM - SIZE OF THIS BLOCK.                                        C
C  ▶ NX   - CARTESIAN DERIVATIVE DIRECTION.                            C
C  ▶ LR   - BASIS FUNCTION THAT PARTIAL DERIVATIVE APPLIES TO (L/R).   C
C  OUTPUT:                                                             C
C  ▶ ESG  - UN-NORMALISED EXPANSION COEFFICIENTS FOR THIS BLOCK.       C
C  ▶ GSG  - UN-NORMALISED DERIVATIVE COEFFICIENTS FOR THIS BLOCK.      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 ESG(MB2,MEQ),ETEMP(MB2*MRC)
      COMPLEX*16 GSG(MB2,MEQ),GTEMP(MB2*MRC)
C
      DIMENSION NBAS(2),LQN(2),MQN(2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
      NDRV = 1
      LMAX = LQN(1)+LQN(2)+NDRV
      NTUV = (LMAX+1)*(LMAX+2)*(LMAX+3)/6
C
C     INITIALISE ESG AND GSG ARRAYS
      DO M=1,MAXM
        DO ITUV=1,NTUV
          ESG(M,ITUV) = DCMPLX(0.0D0,0.0D0)
          GSG(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     IMPORT LQN AND BASIS PAIR MQNS FOR LOCAL USE
      LQNA = LQN(1)
      LQNB = LQN(2)
      MQNA = MQN(1)
      MQNB = MQN(2)
C
C     EXIT IF LQN,MQN ORDERS YIELD ZERO CG-COEFFICIENTS
      IF(IABS(MQN(1)).GT.LQN(1).OR.IABS(MQN(2)).GT.LQN(2)) RETURN
C
C     CHECK THAT LMAX IS WITHIN THE BOUNDS OF MKP
      IF(LMAX.GT.MKP+1) THEN
        WRITE(6, *) 'In GVRS: LMAX exceeds MKP+1 parameter.',LMAX,MKP+1
        WRITE(7, *) 'In GVRS: LMAX exceeds MKP+1 parameter.',LMAX,MKP+1
        STOP
      ENDIF
C
C     INITIALISE TEMPORARY ARRAYS
      DO M=1,MRC*MB2
        ETEMP(M) = DCMPLX(0.0D0,0.0D0)
        GTEMP(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     SET STARTING POINT E[0,0;0,0;0,0,0] = RKAB
C     SET STARTING POINT G[0,0;0,0;0,0,0] = ∂/∂A_NX  E[0,0;0,0;0,0,0]
      DO M=1,MAXM
        ETEMP(M) = DCMPLX(RKAB(M),0.0D0)
        GTEMP(M) = DCMPLX(2.0D0*UV(M)*ABDST(NX)*RKAB(M)/P(M),0.0D0)
      ENDDO
C
C     RESULT OF CARTESIAN KRONECKER DELTAS
      MX = KRONECK(1,NX)
      MY = KRONECK(2,NX)
      MZ = KRONECK(3,NX)
C
C     SET 1ST RECURRENCE G[0,0;0,0;1,0,0] = EI/(EI+EJ) E[0,0;0,0;0,0,0]
      ISTART = (IABC(MX,MY,MZ)-1)*MAXM
      IF(LR.EQ.'L') THEN
        DO M=1,MAXM
          K = ISTART+M
          GTEMP(K) = DCMPLX(EIBS(M)*RKAB(M)/P(M),0.0D0)
        ENDDO
      ELSEIF(LR.EQ.'R') THEN
        DO M=1,MAXM
          K = ISTART+M
          GTEMP(K) = DCMPLX(EJBS(M)*RKAB(M)/P(M),0.0D0)
        ENDDO
      ENDIF
C
C     INITIALISE STARTING ADDRESS AND EXPANSION DEGREE
      ISTART = 0
      LAM    = 1
C
C     STEP 1:
C     GENERATE E/G[|MQNA|,MQNA;0,0] FROM E/G[0,0;0,0] USING
C     SIMULTANEOUS STEP OF LQN AND MQN ON CENTRE A
      CALL GSTEPLM(ETEMP,GTEMP,LAM,ISTART,MQNA,MAXM,1,NX,LR)
C
C     STEP 2:
C     GENERATE E/G[LQNA,MQNA;0,0] FROM E/G[|MQNA|,MQNA;0,0]
C     USING THE STEP OF LQN ONLY ON CENTRE A
      CALL GSTEPL(ETEMP,GTEMP,LAM,ISTART,LQNA,MQNA,MAXM,1,NX,LR)
C
C     STEP 3:
C     GENERATE E/G[LQNA,MQNA;|MQNB|,MQNB] FROM E/G[LQNA,MQNA;0,0]
C     SIMULTANEOUS STEP OF LQN AND MQN ON CENTRE B
      CALL GSTEPLM(ETEMP,GTEMP,LAM,ISTART,MQNB,MAXM,2,NX,LR)
C
C     STEP 4:
C     GENERATE E/G[LQNA,MQNA;LQNB,MQNB] FROM E/G[LQNA,MQNA;|MQNB|,MQNB]
C     USING THE STEP OF LQN ONLY ON CENTRE B
      CALL GSTEPL(ETEMP,GTEMP,LAM,ISTART,LQNB,MQNB,MAXM,2,NX,LR)
C
C     NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS DEGREE LAM
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     STEP 5:
C     COPY FINAL BLOCK OF ENTRIES AS THE REQUIRED OUTPUT
      K = 0
      DO ITUV=1,NTUV
        DO M=1,MAXM
          K = K+1
          ESG(M,ITUV) = ETEMP(ISTART+K)
          GSG(M,ITUV) = GTEMP(ISTART+K)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ESTEPLM(ETEMP,LAM,ISTART,MQN,MAXM,IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   EEEEEEEE SSSSSS TTTTTTTT EEEEEEEE PPPPPPP  LL       MM       MM    C
C   EE      SS    SS   TT    EE       PP    PP LL       MMM     MMM    C
C   EE      SS         TT    EE       PP    PP LL       MMMM   MMMM    C
C   EEEEEE   SSSSSS    TT    EEEEEE   PP    PP LL       MM MM MM MM    C
C   EE            SS   TT    EE       PPPPPPP  LL       MM  MMM  MM    C
C   EE      SS    SS   TT    EE       PP       LL       MM   M   MM    C
C   EEEEEEEE SSSSSS    TT    EEEEEEEE PP       LLLLLLLL MM       MM    C
C                                                                      C
C -------------------------------------------------------------------- C
C  SIMULTANEOUSLY INCREMENT THE QUANTUM NUMBERS LQN & MQN, STARTING    C
C  WITH E[0,0], USING THE RECURSION ALGORITHM OF V.R.SAUNDERS IN       C
C  `MOLECULAR INTEGRALS FOR GAUSSIAN TYPE FUNCTIONS' 1983.             C
C -------------------------------------------------------------------- C
C       E[0,0;IT,IU,IV] -> E[|MQN|+1,±(|MQN|+1);IT,IU,IV]  Eq.(64)     C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ LAM    - LENGTH OF THE INPUT HGTF EXPANSION.                      C
C  ▶ ISTART - COUNTER TO TRACK THE START OF THIS PORTION OF THE LIST.  C
C  ▶ MQN    - MAGNETIC QUANTUM NUMBER.                                 C
C  ▶ MAXM   - NUMBER OF EXPONENT/DENSITY PAIRS.                        C
C  ▶ IZ     - CENTRE TO STEP UP.                                       C
C  OUTPUT:                                                             C
C  ▶ ETEMP  - UN-NORMALISED EXPANSION COEFFICIENTS FOR THIS BLOCK.     C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION PX(MB2),PY(MB2)
C
      COMPLEX*16 CONE
      COMPLEX*16 ETEMP(MB2*MRC)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     IF |MQN|.EQ.0 THEN NO INCREMENT IN (LQN,MQN) IS REQUIRED
      IF(IABS(MQN).EQ.0) RETURN
C
C     IMPORT GEOMETRIC VALUES FOR CENTRE OF INTEREST
      DO M=1,MAXM
        IF(IZ.EQ.1) THEN
          PX(M) = PAX(M)
          PY(M) = PAY(M)
        ELSEIF(IZ.EQ.2) THEN
          PX(M) = PBX(M)
          PY(M) = PBY(M)
        ENDIF
      ENDDO
C
C     PHASE TERM FOR SIGN OF MQN
      PHM = DFLOAT(ISIGN(1,MQN))
C
C     LOOP OVER ALL MAGNETIC NUMBERS UP TO THIS MQN (USE LQN AS COUNTER)
      DO LQN=0,IABS(MQN)-1
C
C       DEGENERACY COUNTER 2*LQN+1
        R2L1 = DFLOAT(2*LQN+1)
C
C       NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS LAM
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C       ISTART LABELS THE PREVIOUS LQN VALUE
C       JSTART LABELS THE CURRENT  LQN VALUE
C       KSTART LABELS THE NEXT     LQN VALUE
        JSTART = ISTART
        KSTART = JSTART + NTUV*MAXM
C
C**********************************************************************C
C                          INDEX MAPPINGS:                             C
C -------------------------------------------------------------------- C
C       J0-> E[LQN  ,LQN  ;IT  ,IU  ,IV]                               C
C       K0-> E[LQN+1,LQN+1;IT  ,IU  ,IV]                               C
C       K1-> E[LQN+1,LQN+1;IT+1,IU  ,IV]                               C
C       K2-> E[LQN+1,LQN+1;IT  ,IU+1,IV]                               C
C       K3-> E[LQN+1,LQN+1;IT-1,IU  ,IV]                               C
C       K4-> E[LQN+1,LQN+1;IT  ,IU-1,IV]                               C
C**********************************************************************C
C
C       LOOP OVER THE HGTF INDICES OF THE SEED LAYER
        DO IOUTER=0,LAM
          DO IT=0,IOUTER
            DO IU=0,IOUTER-IT
              IV = IOUTER-IT-IU
C
C             STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (LQN)
              J0 = JSTART + (IABC(IT  ,IU  ,IV)-1)*MAXM
              K0 = KSTART + (IABC(IT  ,IU  ,IV)-1)*MAXM
              K1 = KSTART + (IABC(IT+1,IU  ,IV)-1)*MAXM
              K2 = KSTART + (IABC(IT  ,IU+1,IV)-1)*MAXM
              IF(IT.NE.0) K3 = KSTART + (IABC(IT-1,IU  ,IV  )-1)*MAXM
              IF(IU.NE.0) K4 = KSTART + (IABC(IT  ,IU-1,IV  )-1)*MAXM
C
C             INVOKE RECURRENCE RELATIONS ON LAYER (LQN+1)
C
C             TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
              DO M=1,MAXM
                T1 = R2L1/P2(M)
                TX = R2L1*PX(M)
                TY = R2L1*PY(M)
                ETEMP(K0+M) = ETEMP(K0+M) +          TX*ETEMP(J0+M)
     &                                    + PHM*CONE*TY*ETEMP(J0+M)
                ETEMP(K1+M) = ETEMP(K1+M) +          T1*ETEMP(J0+M)
                ETEMP(K2+M) = ETEMP(K2+M) + PHM*CONE*T1*ETEMP(J0+M)
              ENDDO
C
C             SPECIAL CASE EXCLUDES IT=0
              IF(IT.NE.0) THEN
                FC = R2L1*DFLOAT(IT)
                DO M=1,MAXM
                  ETEMP(K3+M) = ETEMP(K3+M) +          FC*ETEMP(J0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IU=0
              IF(IU.NE.0) THEN
                FC = R2L1*DFLOAT(IU)
                DO M=1,MAXM
                  ETEMP(K4+M) = ETEMP(K4+M) + PHM*CONE*FC*ETEMP(J0+M)
                ENDDO
              ENDIF
C
C           END OF LOOPS OVER HGTF INDICES
            ENDDO
          ENDDO
        ENDDO
C
C       UPDATE 'PREVIOUS' START VALUE
        ISTART = ISTART + NTUV*MAXM
C
C       UPDATE LAMBDA VALUE
        LAM = LAM+1
C
C     END OF LOOP OVER MQN COUNTER
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE GSTEPLM(ETEMP,GTEMP,LAM,ISTART,MQN,MAXM,IZ,NX,LR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    GGGGGG   SSSSSS TTTTTTTT EEEEEEEE PPPPPPP  LL       MM       MM   C
C   GG    GG SS    SS   TT    EE       PP    PP LL       MMM     MMM   C
C   GG       SS         TT    EE       PP    PP LL       MMMM   MMMM   C
C   GG        SSSSSS    TT    EEEEEE   PP    PP LL       MM MM MM MM   C
C   GG   GGG       SS   TT    EE       PPPPPPP  LL       MM  MMM  MM   C
C   GG    GG SS    SS   TT    EE       PP       LL       MM   M   MM   C
C    GGGGGG   SSSSSS    TT    EEEEEEEE PP       LLLLLLLL MM       MM   C
C                                                                      C
C -------------------------------------------------------------------- C
C  SIMULTANEOUSLY INCREMENT THE QUANTUM NUMBERS LQN & MQN, STARTING    C
C  WITH E[0,0], USING THE RECURSION ALGORITHM OF V.R.SAUNDERS IN       C
C  `MOLECULAR INTEGRALS FOR GAUSSIAN TYPE FUNCTIONS' 1983 (A) AND      C
C  `ANALYTICAL HF GRADIENTS FOR PERIODIC SYSTEMS', JOURNAL OF QUANTUM  C
C  CHEMISTRY, VOL. 82, 1-13, 2001 (B)                                  C
C -------------------------------------------------------------------- C
C       E[0,0;IT,IU,IV] -> E[|MQN|+1,±(|MQN|+1);IT,IU,IV]  Eq.(64 A)   C
C       G[0,0;IT,IU,IV] -> G[|MQN|+1,±(|MQN|+1);IT,IU,IV]  Eq.(64 A)   C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ LAM    - LENGTH OF THE INPUT HGTF EXPANSION.                      C
C  ▶ ISTART - COUNTER TO TRACK THE START OF THIS PORTION OF THE LIST.  C
C  ▶ MQN    - MAGNETIC QUANTUM NUMBER.                                 C
C  ▶ MAXM   - NUMBER OF EXPONENT/DENSITY PAIRS.                        C
C  ▶ IZ     - CENTRE TO STEP UP.                                       C
C  ▶ NX     - CARTESIAN DERIVATIVE DIRECTION.                          C
C  ▶ LR     - BASIS FUNCTION THAT PARTIAL DERIVATIVE APPLIES TO (L/R). C
C  OUTPUT:                                                             C
C  ▶  ETEMP - UN-NORMALISED EXPANSION COEFFICIENTS FOR THIS BLOCK.     C
C  ▶  GTEMP - UN-NORMALISED DERIVATIVE COEFFICIENTS FOR THIS BLOCK.    C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 CONE
      COMPLEX*16 ETEMP(MB2*MRC),GTEMP(MB2*MRC)
C
      DIMENSION PX(MB2),PY(MB2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     IF |MQN|.EQ.0 THEN NO INCREMENT IN (LQN,MQN) IS REQUIRED
      IF(IABS(MQN).EQ.0) RETURN
C
C     IMPORT GEOMETRIC VALUES FOR CENTRE OF INTEREST
      DO M=1,MAXM
        IF(IZ.EQ.1) THEN
          PX(M) = PAX(M)
          PY(M) = PAY(M)
        ELSEIF(IZ.EQ.2) THEN
          PX(M) = PBX(M)
          PY(M) = PBY(M)
        ENDIF
      ENDDO
C
C     PHASE TERM FOR SIGN OF MQN
      PHM = DFLOAT(ISIGN(1,MQN))
C
C     LOOP OVER ALL MAGNETIC NUMBERS UP TO THIS MQN (USE LQN AS COUNTER)
      DO LQN=0,IABS(MQN)-1
C
C       DEGENERACY COUNTER 2*LQN+1
        R2L1 = DFLOAT(2*LQN+1)
C
C       NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS LAM
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C       ISTART LABELS THE PREVIOUS LQN VALUE
C       JSTART LABELS THE CURRENT  LQN VALUE
C       KSTART LABELS THE NEXT     LQN VALUE
        JSTART = ISTART
        KSTART = JSTART + NTUV*MAXM
C
C**********************************************************************C
C                          INDEX MAPPINGS:                             C
C -------------------------------------------------------------------- C
C       J0-> E[LQN  ,LQN  ;IT  ,IU  ,IV]                               C
C       K0-> E[LQN+1,LQN+1;IT  ,IU  ,IV]                               C
C       K1-> E[LQN+1,LQN+1;IT+1,IU  ,IV]                               C
C       K2-> E[LQN+1,LQN+1;IT  ,IU+1,IV]                               C
C       K3-> E[LQN+1,LQN+1;IT-1,IU  ,IV]                               C
C       K4-> E[LQN+1,LQN+1;IT  ,IU-1,IV]                               C
C**********************************************************************C
C
C       LOOP OVER THE HGTF INDICES OF THE SEED LAYER
        DO IOUTER=0,LAM
          DO IT=0,IOUTER
            DO IU=0,IOUTER-IT
              IV = IOUTER-IT-IU
C
C             STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (LQN)
              J0 = JSTART + (IABC(IT  ,IU  ,IV)-1)*MAXM
              K0 = KSTART + (IABC(IT  ,IU  ,IV)-1)*MAXM
              K1 = KSTART + (IABC(IT+1,IU  ,IV)-1)*MAXM
              K2 = KSTART + (IABC(IT  ,IU+1,IV)-1)*MAXM
              IF(IT.NE.0) K3 = KSTART + (IABC(IT-1,IU  ,IV  )-1)*MAXM
              IF(IU.NE.0) K4 = KSTART + (IABC(IT  ,IU-1,IV  )-1)*MAXM
C
C             RECURRENCE RELATIONS ON LAYER (LQN+1) OF E'S AND G'S
C
C             TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
              DO M=1,MAXM
                T1 = R2L1/P2(M)
                TX = R2L1*PX(M)
                TY = R2L1*PY(M)
                ETEMP(K0+M) = ETEMP(K0+M) +          TX*ETEMP(J0+M)
     &                                    + PHM*CONE*TY*ETEMP(J0+M)
                GTEMP(K0+M) = GTEMP(K0+M) +          TX*GTEMP(J0+M)
     &                                    + PHM*CONE*TY*GTEMP(J0+M)
                ETEMP(K1+M) = ETEMP(K1+M) +          T1*ETEMP(J0+M)
                GTEMP(K1+M) = GTEMP(K1+M) +          T1*GTEMP(J0+M)
                ETEMP(K2+M) = ETEMP(K2+M) + PHM*CONE*T1*ETEMP(J0+M)
                GTEMP(K2+M) = GTEMP(K2+M) + PHM*CONE*T1*GTEMP(J0+M)
              ENDDO
C
C             SPECIAL CASE EXCLUDES IT=0
              IF(IT.NE.0) THEN
                FC = R2L1*DFLOAT(IT)
                DO M=1,MAXM
                  ETEMP(K3+M) = ETEMP(K3+M) +          FC*ETEMP(J0+M)
                  GTEMP(K3+M) = GTEMP(K3+M) +          FC*GTEMP(J0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IU=0
              IF(IU.NE.0) THEN
                FC = R2L1*DFLOAT(IU)
                DO M=1,MAXM
                  ETEMP(K4+M) = ETEMP(K4+M) + PHM*CONE*FC*ETEMP(J0+M)
                  GTEMP(K4+M) = GTEMP(K4+M) + PHM*CONE*FC*GTEMP(J0+M)
                ENDDO
              ENDIF
C
C             RECURRENCE RELATIONS ON G'S SOMETIMES HAVE CROSSING TERMS
              IF(NX.EQ.3) GOTO 100
              IF(IZ.EQ.1.AND.LR.EQ.'R') GOTO 100
              IF(IZ.EQ.2.AND.LR.EQ.'L') GOTO 100
C
              IF(NX.EQ.1) THEN
C             CROSSING TERM FROM X-DERIVATIVE
C
                DO M=1,MAXM
                  T0 = R2L1
                  GTEMP(K0+M) = GTEMP(K0+M) -          T0*ETEMP(J0+M)
                ENDDO
C
              ELSEIF(NX.EQ.2) THEN
C             CROSSING TERM FROM Y-DERIVATIVE
C
                DO M=1,MAXM
                  T0 = R2L1
                  GTEMP(K0+M) = GTEMP(K0+M) - PHM*CONE*T0*ETEMP(J0+M)
                ENDDO
C
              ENDIF
C
100           CONTINUE
C
C           END OF LOOPS OVER HGTF INDICES
            ENDDO
          ENDDO
        ENDDO
C
C       UPDATE 'PREVIOUS' START VALUE
        ISTART = ISTART + NTUV*MAXM
C
C       UPDATE LAMBDA VALUE
        LAM = LAM+1
C
C     END OF LOOP OVER MQN COUNTER
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ESTEPL(ETEMP,LAM,ISTART0,LQN,MQN,MAXM,IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         EEEEEEEE SSSSSS TTTTTTTT EEEEEEEE PPPPPPP  LL                C
C         EE      SS    SS   TT    EE       PP    PP LL                C
C         EE      SS         TT    EE       PP    PP LL                C
C         EEEEEE   SSSSSS    TT    EEEEEE   PP    PP LL                C
C         EE            SS   TT    EE       PPPPPPP  LL                C
C         EE      SS    SS   TT    EE       PP       LL                C
C         EEEEEEEE SSSSSS    TT    EEEEEEEE PP       LLLLLLLL          C
C                                                                      C
C -------------------------------------------------------------------- C
C  DECREMENT THE LQN, STARTING WITH E[MQN,±|MQN|], USING THE RECURSION C
C  ALGORITHM OF V.R.SAUNDERS IN `MOLECULAR INTEGRALS FOR GAUSSIAN TYPE C
C  FUNCTIONS' 1983 (EDITED BY G.H.F. DIERCKSEN AND S. WILSON).         C
C -------------------------------------------------------------------- C
C (1)  E[|MQN|,±MQN;IT,IU,IV] -> E[LQN+1,±MQN;IT,IU,IV]    Eq.(64)     C
C (2)  E[LQN  ,±MQN;IT,IU,IV] -> E[LQN+1,±MQN;IT,IU,IV]                C
C (3)  E[LQN-1,±MQN;IT,IU,IV] -> E[LQN+1,±MQN;IT,IU,IV]                C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ LAM     - LENGTH OF THE INPUT HGTF EXPANSION.                     C
C  ▶ ISTART0 - COUNTER TO TRACK THE START OF THIS PORTION OF THE LIST. C
C  ▶ LQN     - ORBITAL QUANTUM NUMBER.                                 C
C  ▶ MQN     - MAGNETIC QUANTUM NUMBER.                                C
C  ▶ MAXM    - NUMBER OF EXPONENT/DENSITY PAIRS.                       C
C  ▶ IZ      - CENTRE TO STEP UP.                                      C
C  OUTPUT:                                                             C
C  ▶ ETEMP   - UN-NORMALISED EXPANSION COEFFICIENTS FOR THIS BLOCK.    C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMPLEX*16 CONE
      COMPLEX*16 ETEMP(MB2*MRC)
C
      DIMENSION PP(MB2),PX(MB2),PY(MB2),PZ(MB2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     IF LQN.LE.|MQN| THEN NO INCREMENT IN LQN IS REQUIRED
      IF(LQN.LE.IABS(MQN)) RETURN
C
C     IMPORT GEOMETRIC VALUES FOR CENTRE OF INTEREST
      DO M=1,MAXM
        IF(IZ.EQ.1) THEN
          PP(M) = PA2(M)
          PX(M) = PAX(M)
          PY(M) = PAY(M)
          PZ(M) = PAZ(M)
        ELSEIF(IZ.EQ.2) THEN
          PP(M) = PB2(M)
          PX(M) = PBX(M)
          PY(M) = PBY(M)
          PZ(M) = PBZ(M)
        ENDIF
      ENDDO
C
C     NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS LAM
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C**********************************************************************C
C     STEP (1): E[|MQN|,MQN;IT,IU,IV] -> E[|MQN|+1,MQN;IT,IU,IV].      C
C               IT MAPS SOME INDEX SETS FROM DATA OBTAINED IN ESTEPLM. C
C -------------------------------------------------------------------- C
C     INDEX MAPPINGS: J0-> E[MQN  ,MQN;IT  ,IU  ,IV  ]                 C
C                     K0-> E[MQN+1,MQN;IT  ,IU  ,IV  ]                 C
C                     K6-> E[MQN+1,MQN;IT  ,IU  ,IV+1]                 C
C                     K9-> E[MQN+1,MQN;IT  ,IU  ,IV-1]                 C
C**********************************************************************C
C
C     ISTART0 LABELS THE GLOBAL STARTING VALUE
C     ISTART  LABELS THE PREVIOUS LQN VALUE
C     JSTART  LABELS THE CURRENT  LQN VALUE
C     KSTART  LABELS THE NEXT     LQN VALUE
      JSTART = ISTART0
      KSTART = JSTART + NTUV*MAXM
C
C     OVERALL LQN/MQN FACTOR SIMPLIFIES WHEN LQN = |MQN|
      RLM  = DFLOAT(2*IABS(MQN)+1)
C
C     LOOP OVER THE HGTF INDICES OF THE SEED LAYER
      DO IOUTER=0,LAM
        DO IT=0,IOUTER
          DO IU=0,IOUTER-IT
            IV = IOUTER-IT-IU
C
C           STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (LQN+1)
            J0 = JSTART + (IABC(IT  ,IU  ,IV  )-1)*MAXM
            K0 = KSTART + (IABC(IT  ,IU  ,IV  )-1)*MAXM
            K6 = KSTART + (IABC(IT  ,IU  ,IV+1)-1)*MAXM
            IF(IV.NE.0) K9 = KSTART + (IABC(IT  ,IU  ,IV-1)-1)*MAXM
C
C           INVOKE RECURRENCE RELATIONS ON LAYER (LQN+1)
C
C           TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
            DO M=1,MAXM
              TZ = RLM*PZ(M)
              TP = RLM/P2(M)
              ETEMP(K0+M) = ETEMP(K0+M) + TZ*ETEMP(J0+M)
              ETEMP(K6+M) = ETEMP(K6+M) + TP*ETEMP(J0+M)
            ENDDO
C
C           SPECIAL CASE EXCLUDES IV=0
            IF(IV.GE.1) THEN
              FC = RLM*DFLOAT(IV)
              DO M=1,MAXM
                ETEMP(K9+M) = ETEMP(K9+M) + FC*ETEMP(J0+M)
              ENDDO
            ENDIF
C
          ENDDO
        ENDDO
      ENDDO
C
C     UPDATE BLOCK LOCATORS
      ISTART  = ISTART0
      ISTART0 = ISTART0 + NTUV*MAXM
C
C     UPDATE LAMBDA VALUE
      LAM = LAM+1
C
C     IF LQN=|MQN|+1 THEN SET IS FINISHED
      IF(LQN.EQ.IABS(MQN)+1) RETURN
C
C**********************************************************************C
C   THIS AND SUBSEQUENT STEPS IN RECURRENCE INVOLVE THREE LAYERS:      C
C            E[LQN  ,MQN;IT,IU,IV] -> E[LQN+1,MQN;IT,IU,IV]            C
C            E[LQN-1,MQN;IT,IU,IV] -> E[LQN+1,MQN;IT,IU,IV]            C
C -------------------------------------------------------------------- C
C   INDEX MAPPINGS: I0 -> E[LQN-1,MQN;IT  ,IU  ,IV  ]                  C
C                   J0 -> E[LQN  ,MQN;IT  ,IU  ,IV  ]                  C
C                   K0 -> E[LQN+1,MQN;IT  ,IU  ,IV  ]                  C
C                   K1 -> E[LQN+1,MQN;IT+2,IU  ,IV  ]                  C
C                   K2 -> E[LQN+1,MQN;IT  ,IU+2,IV  ]                  C
C                   K3 -> E[LQN+1,MQN;IT  ,IU  ,IV+2]                  C
C                   K4 -> E[LQN+1,MQN;IT+1,IU  ,IV  ]                  C
C                   K5 -> E[LQN+1,MQN;IT  ,IU+1,IV  ]                  C
C                   K6 -> E[LQN+1,MQN;IT  ,IU  ,IV+1]                  C
C                   K7 -> E[LQN+1,MQN;IT-1,IU  ,IV  ]                  C
C                   K8 -> E[LQN+1,MQN;IT  ,IU-1,IV  ]                  C
C                   K9 -> E[LQN+1,MQN;IT  ,IU  ,IV-1]                  C
C                   K10-> E[LQN+1,MQN;IT-2,IU  ,IV  ]                  C
C                   K11-> E[LQN+1,MQN;IT  ,IU-2,IV  ]                  C
C                   K12-> E[LQN+1,MQN;IT  ,IU  ,IV-2]                  C
C**********************************************************************C
C
C     LOOP OVER LSTP FOR FIXED MQN DOWN TO LQN-1
      DO LSTP=IABS(MQN)+1,LQN-1

C       OVERALL LSTP/MQN FACTORS
        RLM1 = DFLOAT(2*LSTP+1)/DFLOAT(LSTP-IABS(MQN)+1)
        RLM2 =-DFLOAT(LSTP+IABS(MQN))/DBLE(LSTP-IABS(MQN)+1)
C
C       NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS LAM
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
        JSTART = ISTART0
        KSTART = ISTART0 + MAXM*NTUV
C
C**********************************************************************C
C     STEP (2): E[LSTP  ,MQN;IT,IU,IV] -> E[LSTP+1,MQN;IT,IU,IV].      C
C**********************************************************************C
C
C       LOOP OVER THE HGTF INDICES OF THE SEED LAYER
        DO IOUTER=0,LAM
          DO IT=0,IOUTER
            DO IU=0,IOUTER-IT
              IV = IOUTER-IT-IU
C
C             STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (LSTP)
              J0 = JSTART + (IABC(IT  ,IU  ,IV  )-1)*MAXM
              K0 = KSTART + (IABC(IT  ,IU  ,IV  )-1)*MAXM
              K6 = KSTART + (IABC(IT  ,IU  ,IV+1)-1)*MAXM
              IF(IV.NE.0) K9 = KSTART + (IABC(IT  ,IU  ,IV-1)-1)*MAXM
C
C             INVOKE RECURRENCE RELATIONS ON LAYER (LSTP)
C
C             TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
              DO M=1,MAXM
                TZ = PZ(M)
                TP = 1.0D0/P2(M)
                ETEMP(K0+M) = ETEMP(K0+M) + TZ*RLM1*ETEMP(J0+M)
                ETEMP(K6+M) = ETEMP(K6+M) + TP*RLM1*ETEMP(J0+M)
              ENDDO
C
C             SPECIAL CASE EXCLUDES IV=0
              IF(IV.GE.1) THEN
                FC = DFLOAT(IV)
                DO M=1,MAXM
                  ETEMP(K9+M) = ETEMP(K9+M) + FC*RLM1*ETEMP(J0+M)
                ENDDO
              ENDIF
C
            ENDDO
          ENDDO
        ENDDO
C
C**********************************************************************C
C     STEP (3): E[LSTP-1,MQN;IT,IU,IV] -> E[LSTP+1,MQN;IT,IU,IV].      C
C**********************************************************************C
C
C       LOOP OVER THE HGTF INDICES OF THE SEED LAYER
        DO IOUTER=0,LAM-1
          DO IT=0,IOUTER
            DO IU=0,IOUTER-IT
              IV = IOUTER-IT-IU
C
C             STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (LSTP)
              I0  = ISTART+(IABC(IT  ,IU  ,IV  )-1)*MAXM
              K0  = KSTART+(IABC(IT  ,IU  ,IV  )-1)*MAXM
              K1  = KSTART+(IABC(IT+2,IU  ,IV  )-1)*MAXM
              K2  = KSTART+(IABC(IT  ,IU+2,IV  )-1)*MAXM
              K3  = KSTART+(IABC(IT  ,IU  ,IV+2)-1)*MAXM
              K4  = KSTART+(IABC(IT+1,IU  ,IV  )-1)*MAXM
              K5  = KSTART+(IABC(IT  ,IU+1,IV  )-1)*MAXM
              K6  = KSTART+(IABC(IT  ,IU  ,IV+1)-1)*MAXM
              IF(IT.GT.0) K7  = KSTART + (IABC(IT-1,IU  ,IV  )-1)*MAXM
              IF(IU.GT.0) K8  = KSTART + (IABC(IT  ,IU-1,IV  )-1)*MAXM
              IF(IV.GT.0) K9  = KSTART + (IABC(IT  ,IU  ,IV-1)-1)*MAXM
              IF(IT.GT.1) K10 = KSTART + (IABC(IT-2,IU  ,IV  )-1)*MAXM
              IF(IU.GT.1) K11 = KSTART + (IABC(IT  ,IU-2,IV  )-1)*MAXM
              IF(IV.GT.1) K12 = KSTART + (IABC(IT  ,IU  ,IV-2)-1)*MAXM
C
C             INVOKE RECURRENCE RELATIONS ON LAYER (LSTP-1)
C
C             TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
              TI = DFLOAT(2*(IT+IU+IV)+3)
              DO M=1,MAXM
                T1 = 1.0D0/P22(M)
                TX = PX(M)/P(M)
                TY = PY(M)/P(M)
                TZ = PZ(M)/P(M)
                TT = PP(M) + TI/P2(M)
                ETEMP(K0+M) = ETEMP(K0+M) + RLM2*TT*ETEMP(I0+M)
                ETEMP(K1+M) = ETEMP(K1+M) + RLM2*T1*ETEMP(I0+M)
                ETEMP(K2+M) = ETEMP(K2+M) + RLM2*T1*ETEMP(I0+M)
                ETEMP(K3+M) = ETEMP(K3+M) + RLM2*T1*ETEMP(I0+M)
                ETEMP(K4+M) = ETEMP(K4+M) + RLM2*TX*ETEMP(I0+M)
                ETEMP(K5+M) = ETEMP(K5+M) + RLM2*TY*ETEMP(I0+M)
                ETEMP(K6+M) = ETEMP(K6+M) + RLM2*TZ*ETEMP(I0+M)
              ENDDO
C
C             SPECIAL CASE EXCLUDES IT=0
              IF(IT.GE.1) THEN
                DO M=1,MAXM
                  TX = DFLOAT(2*IT)*PX(M)
                  ETEMP(K7+M) = ETEMP(K7+M) + RLM2*TX*ETEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IU=0
              IF(IU.GE.1) THEN
                DO M=1,MAXM
                  TY = DFLOAT(2*IU)*PY(M)
                  ETEMP(K8+M) = ETEMP(K8+M) + RLM2*TY*ETEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IV=0
              IF(IV.GE.1) THEN
                DO M=1,MAXM
                  TZ = DFLOAT(2*IV)*PZ(M)
                  ETEMP(K9+M) = ETEMP(K9+M) + RLM2*TZ*ETEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IT=0,1
              IF(IT.GE.2) THEN
                TX = DFLOAT(IT*(IT-1))
                DO M=1,MAXM
                  ETEMP(K10+M) = ETEMP(K10+M) + RLM2*TX*ETEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IU=0,1
              IF(IU.GE.2) THEN
                TY = DFLOAT(IU*(IU-1))
                DO M=1,MAXM
                  ETEMP(K11+M) = ETEMP(K11+M) + RLM2*TY*ETEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IV=0,1
              IF(IV.GE.2) THEN
                TZ = DFLOAT(IV*(IV-1))
                DO M=1,MAXM
                  ETEMP(K12+M) = ETEMP(K12+M) + RLM2*TZ*ETEMP(I0+M)
                ENDDO
              ENDIF
C
            ENDDO
          ENDDO
        ENDDO
C
C       UPDATE BLOCK LOCATORS
        ISTART  = ISTART0
        ISTART0 = ISTART0 + MAXM*NTUV
C
C       UPDATE LAMBDA VALUE
        LAM = LAM+1
C
C     END OF LOOP OVER LSTP FOR FIXED MQN
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE GSTEPL(ETEMP,GTEMP,LAM,ISTART0,LQN,MQN,MAXM,IZ,NX,LR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          GGGGGG   SSSSSS TTTTTTTT EEEEEEEE PPPPPPP  LL               C
C         GG    GG SS    SS   TT    EE       PP    PP LL               C
C         GG       SS         TT    EE       PP    PP LL               C
C         GG        SSSSSS    TT    EEEEEE   PP    PP LL               C
C         GG   GGG       SS   TT    EE       PPPPPPP  LL               C
C         GG    GG SS    SS   TT    EE       PP       LL               C
C          GGGGGG   SSSSSS    TT    EEEEEEEE PP       LLLLLLLL         C
C                                                                      C
C -------------------------------------------------------------------- C
C  DECREMENT THE LQN, STARTING WITH G[MQN,±|MQN|], USING THE RECURSION C
C  ALGORITHM OF V.R.SAUNDERS IN `MOLECULAR INTEGRALS FOR GAUSSIAN TYPE C
C  FUNCTIONS' 1983 (EDITED BY G.H.F. DIERCKSEN AND S. WILSON).         C
C -------------------------------------------------------------------- C
C (1)  E/G[|MQN|,±MQN;IT,IU,IV] -> E/G[LQN+1,±MQN;IT,IU,IV]    Eq.(64) C
C (2)  E/G[LQN  ,±MQN;IT,IU,IV] -> E/G[LQN+1,±MQN;IT,IU,IV]            C
C (3)  E/G[LQN-1,±MQN;IT,IU,IV] -> E/G[LQN+1,±MQN;IT,IU,IV]            C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ LAM     - LENGTH OF THE INPUT HGTF EXPANSION.                     C
C  ▶ ISTART0 - COUNTER TO TRACK THE START OF THIS PORTION OF THE LIST. C
C  ▶ LQN     - ORBITAL QUANTUM NUMBER.                                 C
C  ▶ MQN     - MAGNETIC QUANTUM NUMBER.                                C
C  ▶ MAXM    - NUMBER OF EXPONENT/DENSITY PAIRS.                       C
C  ▶ IZ      - CENTRE TO STEP UP.                                      C
C  ▶ NX      - CARTESIAN DERIVATIVE DIRECTION.                         C
C  ▶ LR      - BASIS FUNCTION THAT PARTIAL DERIVATIVE APPLIES TO (L/R).C
C  OUTPUT:                                                             C
C  ▶ ETEMP   - UN-NORMALISED EXPANSION COEFFICIENTS FOR THIS BLOCK.    C
C  ▶ GTEMP - UN-NORMALISED DERIVATIVE COEFFICIENTS FOR THIS BLOCK.     C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 CONE
      COMPLEX*16 ETEMP(MB2*MRC),GTEMP(MB2*MRC)
C
      DIMENSION PP(MB2),PX(MB2),PY(MB2),PZ(MB2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     IF LQN.LE.|MQN| THEN NO INCREMENT IN LQN IS REQUIRED
      IF(LQN.LE.IABS(MQN)) RETURN
C
C     IDENTIFY ADDRESS FROM DERIVATIVE
      MX = KRONECK(1,NX)
      MY = KRONECK(2,NX)
      MZ = KRONECK(3,NX)
C
C     IMPORT GEOMETRIC VALUES FOR CENTRE OF INTEREST
      DO M=1,MAXM
        IF(IZ.EQ.1) THEN
          PP(M) = PA2(M)
          PX(M) = PAX(M)
          PY(M) = PAY(M)
          PZ(M) = PAZ(M)
        ELSEIF(IZ.EQ.2) THEN
          PP(M) = PB2(M)
          PX(M) = PBX(M)
          PY(M) = PBY(M)
          PZ(M) = PBZ(M)
        ENDIF
      ENDDO
C
C     NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS LAM
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C**********************************************************************C
C     STEP (1): E[|MQN|,MQN;IT,IU,IV] -> E[|MQN|+1,MQN;IT,IU,IV].      C
C               G[|MQN|,MQN;IT,IU,IV] -> G[|MQN|+1,MQN;IT,IU,IV].      C
C               IT MAPS SOME INDEX SETS FROM DATA OBTAINED IN GSTEPLM. C
C -------------------------------------------------------------------- C
C     INDEX MAPPINGS: J0-> E/G[MQN  ,MQN;IT  ,IU  ,IV  ]               C
C                     K0-> E/G[MQN+1,MQN;IT  ,IU  ,IV  ]               C
C                     K6-> E/G[MQN+1,MQN;IT  ,IU  ,IV+1]               C
C                     K9-> E/G[MQN+1,MQN;IT  ,IU  ,IV-1]               C
C**********************************************************************C
C
C     ISTART0 LABELS THE GLOBAL STARTING VALUE
C     ISTART  LABELS THE PREVIOUS LQN VALUE
C     JSTART  LABELS THE CURRENT  LQN VALUE
C     KSTART  LABELS THE NEXT     LQN VALUE
      JSTART = ISTART0
      KSTART = JSTART + NTUV*MAXM
C
C     OVERALL LQN/MQN FACTOR SIMPLIFIES WHEN LQN = |MQN|
      RLM  = DFLOAT(2*IABS(MQN)+1)
C
C     LOOP OVER THE HGTF INDICES OF THE SEED LAYER
      DO IOUTER=0,LAM
        DO IT=0,IOUTER
          DO IU=0,IOUTER-IT
            IV = IOUTER-IT-IU
C
C           STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (LQN+1)
            J0 = JSTART + (IABC(IT  ,IU  ,IV  )-1)*MAXM
            K0 = KSTART + (IABC(IT  ,IU  ,IV  )-1)*MAXM
            K6 = KSTART + (IABC(IT  ,IU  ,IV+1)-1)*MAXM
            IF(IV.NE.0) K9 = KSTART + (IABC(IT  ,IU  ,IV-1)-1)*MAXM
C
C           RECURRENCE RELATIONS ON LAYER (LQN+1) OF E'S AND G'S
C
C           TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
            DO M=1,MAXM
              TZ = RLM*PZ(M)
              TP = RLM/P2(M)
              ETEMP(K0+M) = ETEMP(K0+M) + TZ*ETEMP(J0+M)
              GTEMP(K0+M) = GTEMP(K0+M) + TZ*GTEMP(J0+M)
              ETEMP(K6+M) = ETEMP(K6+M) + TP*ETEMP(J0+M)
              GTEMP(K6+M) = GTEMP(K6+M) + TP*GTEMP(J0+M)
            ENDDO
C
C           SPECIAL CASE EXCLUDES IV=0
            IF(IV.GE.1) THEN
              FC = RLM*DFLOAT(IV)
              DO M=1,MAXM
                ETEMP(K9+M) = ETEMP(K9+M) + FC*ETEMP(J0+M)
                GTEMP(K9+M) = GTEMP(K9+M) + FC*GTEMP(J0+M)
              ENDDO
            ENDIF
C
C           RECURRENCE RELATIONS ON G'S SOMETIMES HAVE CROSSING TERMS
            IF(NX.NE.3) GOTO 100
            IF(IZ.EQ.1.AND.LR.EQ.'R') GOTO 100
            IF(IZ.EQ.2.AND.LR.EQ.'L') GOTO 100
C
C           CROSSING TERM FROM Z-DERIVATIVE
            T0 = RLM
            DO M=1,MAXM
              GTEMP(K0+M) = GTEMP(K0+M) - T0*ETEMP(J0+M)
            ENDDO
C
100         CONTINUE
C
          ENDDO
        ENDDO
      ENDDO
C
C     UPDATE BLOCK LOCATORS
      ISTART  = ISTART0
      ISTART0 = ISTART0 + NTUV*MAXM
C
C     UPDATE LAMBDA VALUE
      LAM = LAM+1
C
C     IF LQN=|MQN|+1 THEN SET IS FINISHED
      IF(LQN.EQ.IABS(MQN)+1) RETURN
C
C**********************************************************************C
C   THIS AND SUBSEQUENT STEPS IN RECURRENCE INVOLVE THREE LAYERS:      C
C            E[LQN  ,MQN;IT,IU,IV] -> E[LQN+1,MQN;IT,IU,IV]            C
C            E[LQN-1,MQN;IT,IU,IV] -> E[LQN+1,MQN;IT,IU,IV]            C
C -------------------------------------------------------------------- C
C   INDEX MAPPINGS: I0 -> E[LQN-1,MQN;IT  ,IU  ,IV  ]                  C
C                   J0 -> E[LQN  ,MQN;IT  ,IU  ,IV  ]                  C
C                   K0 -> E[LQN+1,MQN;IT  ,IU  ,IV  ]                  C
C                   K1 -> E[LQN+1,MQN;IT+2,IU  ,IV  ]                  C
C                   K2 -> E[LQN+1,MQN;IT  ,IU+2,IV  ]                  C
C                   K3 -> E[LQN+1,MQN;IT  ,IU  ,IV+2]                  C
C                   K4 -> E[LQN+1,MQN;IT+1,IU  ,IV  ]                  C
C                   K5 -> E[LQN+1,MQN;IT  ,IU+1,IV  ]                  C
C                   K6 -> E[LQN+1,MQN;IT  ,IU  ,IV+1]                  C
C                   K7 -> E[LQN+1,MQN;IT-1,IU  ,IV  ]                  C
C                   K8 -> E[LQN+1,MQN;IT  ,IU-1,IV  ]                  C
C                   K9 -> E[LQN+1,MQN;IT  ,IU  ,IV-1]                  C
C                   K10-> E[LQN+1,MQN;IT-2,IU  ,IV  ]                  C
C                   K11-> E[LQN+1,MQN;IT  ,IU-2,IV  ]                  C
C                   K12-> E[LQN+1,MQN;IT  ,IU  ,IV-2]                  C
C**********************************************************************C
C
C     LOOP OVER LQN FOR FIXED MQN
      DO LQN=IABS(MQN)+1,LQN-1

C       OVERALL LQN/MQN FACTORS
        RLM1 = DFLOAT(2*LQN+1)/DFLOAT(LQN-IABS(MQN)+1)
        RLM2 =-DFLOAT(LQN+IABS(MQN))/DBLE(LQN-IABS(MQN)+1)
C
C       NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS LAM
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
        JSTART = ISTART0
        KSTART = ISTART0 + MAXM*NTUV
C
C**********************************************************************C
C     STEP (2): E[LQN  ,MQN;IT,IU,IV] -> E[LQN+1,MQN;IT,IU,IV].        C
C**********************************************************************C
C
C       LOOP OVER THE HGTF INDICES OF THE SEED LAYER
        DO IOUTER=0,LAM
          DO IT=0,IOUTER
            DO IU=0,IOUTER-IT
              IV = IOUTER-IT-IU
C
C             STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (LQN)
              J0 = JSTART + (IABC(IT  ,IU  ,IV  )-1)*MAXM
              K0 = KSTART + (IABC(IT  ,IU  ,IV  )-1)*MAXM
              K6 = KSTART + (IABC(IT  ,IU  ,IV+1)-1)*MAXM
              IF(IV.NE.0) K9 = KSTART + (IABC(IT  ,IU  ,IV-1)-1)*MAXM
C
C             RECURRENCE RELATIONS ON LAYER (LQN) OF E'S AND G'S
C
C             TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
              DO M=1,MAXM
                TZ = PZ(M)
                TP = 1.0D0/P2(M)
                ETEMP(K0+M) = ETEMP(K0+M) + TZ*RLM1*ETEMP(J0+M)
                GTEMP(K0+M) = GTEMP(K0+M) + TZ*RLM1*GTEMP(J0+M)
                ETEMP(K6+M) = ETEMP(K6+M) + TP*RLM1*ETEMP(J0+M)
                GTEMP(K6+M) = GTEMP(K6+M) + TP*RLM1*GTEMP(J0+M)
              ENDDO
C
C             SPECIAL CASE EXCLUDES IV=0
              IF(IV.GE.1) THEN
                FC = DFLOAT(IV)
                DO M=1,MAXM
                  ETEMP(K9+M) = ETEMP(K9+M) + FC*RLM1*ETEMP(J0+M)
                  GTEMP(K9+M) = GTEMP(K9+M) + FC*RLM1*GTEMP(J0+M)
                ENDDO
              ENDIF
C
C             RECURRENCE RELATIONS ON G'S SOMETIMES HAVE CROSSING TERMS
              IF(NX.NE.3) GOTO 200
              IF(IZ.EQ.1.AND.LR.EQ.'R') GOTO 200
              IF(IZ.EQ.2.AND.LR.EQ.'L') GOTO 200
C
C             CROSSING TERM FROM Z-DERIVATIVE
              T0 = RLM1
              DO M=1,MAXM
                GTEMP(K0+M) = GTEMP(K0+M) - T0*ETEMP(J0+M)
              ENDDO
C
200           CONTINUE
C
            ENDDO
          ENDDO
        ENDDO
C
C**********************************************************************C
C     STEP (3): E[LQN-1,MQN;IT,IU,IV] -> E[LQN+1,MQN;IT,IU,IV].        C
C**********************************************************************C
C
C       LOOP OVER THE HGTF INDICES OF THE SEED LAYER
        DO IOUTER=0,LAM-1
          DO IT=0,IOUTER
            DO IU=0,IOUTER-IT
              IV = IOUTER-IT-IU
C
C             STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (LQN)
              I0  = ISTART+(IABC(IT  ,IU  ,IV  )-1)*MAXM
              K0  = KSTART+(IABC(IT  ,IU  ,IV  )-1)*MAXM
              K1  = KSTART+(IABC(IT+2,IU  ,IV  )-1)*MAXM
              K2  = KSTART+(IABC(IT  ,IU+2,IV  )-1)*MAXM
              K3  = KSTART+(IABC(IT  ,IU  ,IV+2)-1)*MAXM
              K4  = KSTART+(IABC(IT+1,IU  ,IV  )-1)*MAXM
              K5  = KSTART+(IABC(IT  ,IU+1,IV  )-1)*MAXM
              K6  = KSTART+(IABC(IT  ,IU  ,IV+1)-1)*MAXM
              KN  = KSTART+(IABC(IT+MX,IU+MY,IV+MZ)-1)*MAXM
              IF(IT.GT.0) K7  = KSTART+(IABC(IT-1,IU  ,IV  )-1)*MAXM
              IF(IU.GT.0) K8  = KSTART+(IABC(IT  ,IU-1,IV  )-1)*MAXM
              IF(IV.GT.0) K9  = KSTART+(IABC(IT  ,IU  ,IV-1)-1)*MAXM
              IF(IT.GT.1) K10 = KSTART+(IABC(IT-2,IU  ,IV  )-1)*MAXM
              IF(IU.GT.1) K11 = KSTART+(IABC(IT  ,IU-2,IV  )-1)*MAXM
              IF(IV.GT.1) K12 = KSTART+(IABC(IT  ,IU  ,IV-2)-1)*MAXM
              IF(IT-MX.GE.0.AND.IU-MY.GE.0.AND.IV-MZ.GE.0) THEN
                KM = KSTART+(IABC(IT-MX,IU-MY,IV-MZ)-1)*MAXM
              ENDIF
C
C             RECURRENCE RELATIONS ON LAYER (LQN-1) OF E'S AND G'S
C
C             TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
              TI = DFLOAT(2*(IT+IU+IV)+3)
              DO M=1,MAXM
                T1 = 1.0D0/P22(M)
                TX = PX(M)/P(M)
                TY = PY(M)/P(M)
                TZ = PZ(M)/P(M)
                TT = PP(M) + TI/P2(M)
                ETEMP(K0+M) = ETEMP(K0+M) + RLM2*TT*ETEMP(I0+M)
                GTEMP(K0+M) = GTEMP(K0+M) + RLM2*TT*GTEMP(I0+M)
                ETEMP(K1+M) = ETEMP(K1+M) + RLM2*T1*ETEMP(I0+M)
                GTEMP(K1+M) = GTEMP(K1+M) + RLM2*T1*GTEMP(I0+M)
                ETEMP(K2+M) = ETEMP(K2+M) + RLM2*T1*ETEMP(I0+M)
                GTEMP(K2+M) = GTEMP(K2+M) + RLM2*T1*GTEMP(I0+M)
                ETEMP(K3+M) = ETEMP(K3+M) + RLM2*T1*ETEMP(I0+M)
                GTEMP(K3+M) = GTEMP(K3+M) + RLM2*T1*GTEMP(I0+M)
                ETEMP(K4+M) = ETEMP(K4+M) + RLM2*TX*ETEMP(I0+M)
                GTEMP(K4+M) = GTEMP(K4+M) + RLM2*TX*GTEMP(I0+M)
                ETEMP(K5+M) = ETEMP(K5+M) + RLM2*TY*ETEMP(I0+M)
                GTEMP(K5+M) = GTEMP(K5+M) + RLM2*TY*GTEMP(I0+M)
                ETEMP(K6+M) = ETEMP(K6+M) + RLM2*TZ*ETEMP(I0+M)
                GTEMP(K6+M) = GTEMP(K6+M) + RLM2*TZ*GTEMP(I0+M)
              ENDDO
C
C             SPECIAL CASE EXCLUDES IT=0
              IF(IT.GE.1) THEN
                DO M=1,MAXM
                  TX = DFLOAT(2*IT)*PX(M)
                  ETEMP(K7+M) = ETEMP(K7+M) + RLM2*TX*ETEMP(I0+M)
                  GTEMP(K7+M) = GTEMP(K7+M) + RLM2*TX*GTEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IU=0
              IF(IU.GE.1) THEN
                DO M=1,MAXM
                  TY = DFLOAT(2*IU)*PY(M)
                  ETEMP(K8+M) = ETEMP(K8+M) + RLM2*TY*ETEMP(I0+M)
                  GTEMP(K8+M) = GTEMP(K8+M) + RLM2*TY*GTEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IV=0
              IF(IV.GE.1) THEN
                DO M=1,MAXM
                  TZ = DFLOAT(2*IV)*PZ(M)
                  ETEMP(K9+M) = ETEMP(K9+M) + RLM2*TZ*ETEMP(I0+M)
                  GTEMP(K9+M) = GTEMP(K9+M) + RLM2*TZ*GTEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IT=0,1
              IF(IT.GE.2) THEN
                TX = DFLOAT(IT*(IT-1))
                DO M=1,MAXM
                  ETEMP(K10+M) = ETEMP(K10+M) + RLM2*TX*ETEMP(I0+M)
                  GTEMP(K10+M) = GTEMP(K10+M) + RLM2*TX*GTEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IU=0,1
              IF(IU.GE.2) THEN
                TY = DFLOAT(IU*(IU-1))
                DO M=1,MAXM
                  ETEMP(K11+M) = ETEMP(K11+M) + RLM2*TY*ETEMP(I0+M)
                  GTEMP(K11+M) = GTEMP(K11+M) + RLM2*TY*GTEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IV=0,1
              IF(IV.GE.2) THEN
                TZ = DFLOAT(IV*(IV-1))
                DO M=1,MAXM
                  ETEMP(K12+M) = ETEMP(K12+M) + RLM2*TZ*ETEMP(I0+M)
                  GTEMP(K12+M) = GTEMP(K12+M) + RLM2*TZ*GTEMP(I0+M)
                ENDDO
              ENDIF
C
C             RECURRENCE RELATIONS ON G'S SOMETIMES HAVE CROSSING TERMS
              IF(IZ.EQ.1.AND.LR.EQ.'R') GOTO 300
              IF(IZ.EQ.2.AND.LR.EQ.'L') GOTO 300
C
C             TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
              DO M=1,MAXM
                T0 = 1.0D0/P(M)
                TN = 2.0D0*(MX*PX(M) + MY*PY(M) + MZ*PZ(M))
                GTEMP(KN+M) = GTEMP(KN+M) - RLM2*T0*ETEMP(I0+M)
                GTEMP(K0+M) = GTEMP(K0+M) - RLM2*TN*ETEMP(I0+M)
              ENDDO
C
C             SPECIAL CASE EXCLUDES IX=0
              IF(IT-MX.GE.0.AND.IU-MY.GE.0.AND.IV-MZ.GE.0) THEN
                TM = 2.0D0*DFLOAT(MX*IT + MY*IU + MZ*IV)
                DO M=1,MAXM
                  GTEMP(KM+M) = GTEMP(KM+M) - RLM2*TM*ETEMP(I0+M)
                ENDDO
              ENDIF
C
300           CONTINUE
C
            ENDDO
          ENDDO
        ENDDO
C
C       UPDATE BLOCK LOCATORS
        ISTART  = ISTART0
        ISTART0 = ISTART0 + MAXM*NTUV
C
C       UPDATE LAMBDA VALUE
        LAM = LAM+1
C
C     END OF LOOP OVER LQN FOR FIXED MQN
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ESTEPN(ESG,ENSG,LAM,MAXM,IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         EEEEEEEE SSSSSS TTTTTTTT EEEEEEEE PPPPPPP  NN    NN          C
C         EE      SS    SS   TT    EE       PP    PP NNN   NN          C
C         EE      SS         TT    EE       PP    PP NNNN  NN          C
C         EEEEEE   SSSSSS    TT    EEEEEE   PP    PP NN NN NN          C
C         EE            SS   TT    EE       PPPPPPP  NN  NNNN          C
C         EE      SS    SS   TT    EE       PP       NN   NNN          C
C         EEEEEEEE SSSSSS    TT    EEEEEEEE PP       NN    NN          C
C                                                                      C
C -------------------------------------------------------------------- C
C  INCREMENT THE QUANTUM NUMBER (NQN):                                 C
C                E[NQN  ,LQN,MQN] -> E[NQN+1,LQN,MQN].                 C
C -------------------------------------------------------------------- C
C  ▶ ONLY PERFORMS A SINGLE STEP IN NQN.                               C
C  ▶ LAM IS THE EFFECTIVE TOTAL ANGULAR MOMENTUM OF THE INPUT COEFFS.  C
C    EFFECTIVE TOTAL ANGULAR MOMENTUM OF THE OUTPUT COEFFS IS LAM+2.   C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ ESG  - EQ-COEFFICIENT BATCH.                                      C
C  ▶ LAM  - EFFECTIVE TOTAL ANGULAR MOMENTUM.                          C
C  ▶ MAXM - NUMBER OF EXPONENT/DENSITY PAIRS.                          C
C  ▶ IZ   - CENTRE TO STEP UP.                                         C
C  OUTPUT:                                                             C
C  ▶ ENSG - EQ-COEFFICIENT BATCH AFTER NQN HAS BEEN STEPPED UP.        C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION PP(MB2),PX(MB2),PY(MB2),PZ(MB2)
C
      COMPLEX*16 ESG(MB2,MEQ),ENSG(MB2,MEQ)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS LAM WITH N'=N+1
      NTUV = (LAM+3)*(LAM+4)*(LAM+5)/6
C
C     INITIALISE NEW ARRAY
      DO ITUV=1,NTUV
        DO M=1,MAXM
          ENSG(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     IMPORT GEOMETRIC VALUES FOR CENTRE OF INTEREST
      DO M=1,MAXM
        IF(IZ.EQ.1) THEN
          PP(M) = PA2(M)
          PX(M) = PAX(M)
          PY(M) = PAY(M)
          PZ(M) = PAZ(M)
        ELSEIF(IZ.EQ.2) THEN
          PP(M) = PB2(M)
          PX(M) = PBX(M)
          PY(M) = PBY(M)
          PZ(M) = PBZ(M)
        ENDIF
      ENDDO
C
C**********************************************************************C
C     INDEX MAPPINGS: I0 -> E[LQN-1,MQN;IT  ,IU  ,IV  ]                C
C                     K0 -> E[LQN+1,MQN;IT  ,IU  ,IV  ]                C
C                     K1 -> E[LQN+1,MQN;IT+2,IU  ,IV  ]                C
C                     K2 -> E[LQN+1,MQN;IT  ,IU+2,IV  ]                C
C                     K3 -> E[LQN+1,MQN;IT  ,IU  ,IV+2]                C
C                     K4 -> E[LQN+1,MQN;IT+1,IU  ,IV  ]                C
C                     K5 -> E[LQN+1,MQN;IT  ,IU+1,IV  ]                C
C                     K6 -> E[LQN+1,MQN;IT  ,IU  ,IV+1]                C
C                     K7 -> E[LQN+1,MQN;IT-1,IU  ,IV  ]                C
C                     K8 -> E[LQN+1,MQN;IT  ,IU-1,IV  ]                C
C                     K9 -> E[LQN+1,MQN;IT  ,IU  ,IV-1]                C
C                     K10-> E[LQN+1,MQN;IT-2,IU  ,IV  ]                C
C                     K11-> E[LQN+1,MQN;IT  ,IU-2,IV  ]                C
C                     K12-> E[LQN+1,MQN;IT  ,IU  ,IV-2]                C
C**********************************************************************C
C
      DO IOUTER=0,LAM
        DO IT=0,IOUTER
          DO IU=0,IOUTER-IT
            IV = IOUTER-IT-IU
C
C           STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (NQN)
            I0  = IABC(IT  ,IU  ,IV  )
            K0  = IABC(IT  ,IU  ,IV  )
            K1  = IABC(IT+2,IU  ,IV  )
            K2  = IABC(IT  ,IU+2,IV  )
            K3  = IABC(IT  ,IU  ,IV+2)
            K4  = IABC(IT+1,IU  ,IV  )
            K5  = IABC(IT  ,IU+1,IV  )
            K6  = IABC(IT  ,IU  ,IV+1)
            IF(IT.GT.0) K7  = IABC(IT-1,IU  ,IV  )
            IF(IU.GT.0) K8  = IABC(IT  ,IU-1,IV  )
            IF(IV.GT.0) K9  = IABC(IT  ,IU  ,IV-1)
            IF(IT.GT.1) K10 = IABC(IT-2,IU  ,IV  )
            IF(IU.GT.1) K11 = IABC(IT  ,IU-2,IV  )
            IF(IV.GT.1) K12 = IABC(IT  ,IU  ,IV-2)
C
C           INVOKE RECURRENCE RELATIONS ON LAYER (NQN+1)
C
C           TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
            TT = DFLOAT((2*(IT+IU+IV))+3)
            DO M=1,MAXM
              T1 = 1.0D0/P22(M)
              TX = PX(M)/P(M)
              TY = PY(M)/P(M)
              TZ = PZ(M)/P(M)
              TP = PP(M) + TT/P2(M)
              ENSG(M,K0) = ENSG(M,K0) + TP*ESG(M,I0)
              ENSG(M,K1) = ENSG(M,K1) + T1*ESG(M,I0)
              ENSG(M,K2) = ENSG(M,K2) + T1*ESG(M,I0)
              ENSG(M,K3) = ENSG(M,K3) + T1*ESG(M,I0)
              ENSG(M,K4) = ENSG(M,K4) + TX*ESG(M,I0)
              ENSG(M,K5) = ENSG(M,K5) + TY*ESG(M,I0)
              ENSG(M,K6) = ENSG(M,K6) + TZ*ESG(M,I0)
            ENDDO
C
C           SPECIAL CASE EXCLUDES IT=0
            IF(IT.GE.1) THEN
              RT2 = DFLOAT(2*IT)
              DO M=1,MAXM
                T0 = PX(M)*RT2
                ENSG(M,K7) = ENSG(M,K7) + T0*ESG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IU=0
            IF(IU.GE.1) THEN
              RU1 = DFLOAT(2*IU)
              DO M=1,MAXM
                T0 = PY(M)*RU1
                ENSG(M,K8) = ENSG(M,K8) + T0*ESG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IV=0
            IF(IV.GE.1) THEN
              RV1 = DFLOAT(2*IV)
              DO M=1,MAXM
                T0 = PZ(M)*RV1
                ENSG(M,K9) = ENSG(M,K9) + T0*ESG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IT=0,1
            IF(IT.GE.2) THEN
              RT2 = DFLOAT(IT*(IT-1))
              DO M=1,MAXM
                ENSG(M,K10) = ENSG(M,K10) + RT2*ESG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IU=0,1
            IF(IU.GE.2) THEN
              RU2 = DFLOAT(IU*(IU-1))
              DO M=1,MAXM
                ENSG(M,K11) = ENSG(M,K11) + RU2*ESG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IV=0,1
            IF(IV.GE.2) THEN
              RV2 = DFLOAT(IV*(IV-1))
              DO M=1,MAXM
                ENSG(M,K12) = ENSG(M,K12) + RV2*ESG(M,I0)
              ENDDO
            ENDIF

          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE GSTEPN(ESG,ENSG,GSG,GNSG,LAM,MAXM,IZ,NX,LR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          GGGGGG   SSSSSS TTTTTTTT EEEEEEEE PPPPPPP  NN    NN         C
C         GG    GG SS    SS   TT    EE       PP    PP NNN   NN         C
C         GG       SS         TT    EE       PP    PP NNNN  NN         C
C         GG        SSSSSS    TT    EEEEEE   PP    PP NN NN NN         C
C         GG   GGG       SS   TT    EE       PPPPPPP  NN  NNNN         C
C         GG    GG SS    SS   TT    EE       PP       NN   NNN         C
C          GGGGGG   SSSSSS    TT    EEEEEEEE PP       NN    NN         C
C                                                                      C
C -------------------------------------------------------------------- C
C  INCREMENT THE QUANTUM NUMBER (NQN):                                 C
C              E/G[NQN  ,LQN,MQN] -> E/G[NQN+1,LQN,MQN].               C
C -------------------------------------------------------------------- C
C  ▶ ONLY PERFORMS A SINGLE STEP IN NQN.                               C
C  ▶ LAM IS THE EFFECTIVE TOTAL ANGULAR MOMENTUM OF THE INPUT COEFFS.  C
C    EFFECTIVE TOTAL ANGULAR MOMENTUM OF THE OUTPUT COEFFS IS LAM+2.   C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ ESG  - EQ-COEFFICIENT BATCH.                                      C
C  ▶ GSG  - GQ-COEFFICIENT BATCH.                                      C
C  ▶ LAM  - EFFECTIVE TOTAL ANGULAR MOMENTUM.                          C
C  ▶ MAXM - NUMBER OF EXPONENT/DENSITY PAIRS.                          C
C  ▶ IZ   - CENTRE TO STEP UP.                                         C
C  ▶ NX   - CARTESIAN DERIVATIVE DIRECTION.                            C
C  ▶ LR   - BASIS FUNCTION THAT PARTIAL DERIVATIVE APPLIES TO (L/R).   C
C  OUTPUT:                                                             C
C  ▶ ENSG - EQ-COEFFICIENT BATCH AFTER NQN HAS BEEN STEPPED UP.        C
C  ▶ GNSG - GQ-COEFFICIENT BATCH AFTER NQN HAS BEEN STEPPED UP.        C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 LR
C
      COMPLEX*16 ESG(MB2,MEQ),ENSG(MB2,MEQ)
      COMPLEX*16 GSG(MB2,MEQ),GNSG(MB2,MEQ)
C
      DIMENSION PP(MB2),PX(MB2),PY(MB2),PZ(MB2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS LAM WITH N'=N+1
      NTUV = (LAM+3)*(LAM+4)*(LAM+5)/6
C
C     INITIALISE NEW ARRAY
      DO ITUV=1,NTUV
        DO M=1,MAXM
          ENSG(M,ITUV) = DCMPLX(0.0D0,0.0D0)
          GNSG(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     IDENTIFY ADDRESS FROM DERIVATIVE
      MX = KRONECK(1,NX)
      MY = KRONECK(2,NX)
      MZ = KRONECK(3,NX)
C
C     IMPORT GEOMETRIC VALUES FOR CENTRE OF INTEREST
      DO M=1,MAXM
        IF(IZ.EQ.1) THEN
          PP(M) = PA2(M)
          PX(M) = PAX(M)
          PY(M) = PAY(M)
          PZ(M) = PAZ(M)
        ELSEIF(IZ.EQ.2) THEN
          PP(M) = PB2(M)
          PX(M) = PBX(M)
          PY(M) = PBY(M)
          PZ(M) = PBZ(M)
        ENDIF
      ENDDO
C
C**********************************************************************C
C     INDEX MAPPINGS: I0 -> E[LQN-1,MQN;IT  ,IU  ,IV  ]                C
C                     K0 -> E[LQN+1,MQN;IT  ,IU  ,IV  ]                C
C                     K1 -> E[LQN+1,MQN;IT+2,IU  ,IV  ]                C
C                     K2 -> E[LQN+1,MQN;IT  ,IU+2,IV  ]                C
C                     K3 -> E[LQN+1,MQN;IT  ,IU  ,IV+2]                C
C                     K4 -> E[LQN+1,MQN;IT+1,IU  ,IV  ]                C
C                     K5 -> E[LQN+1,MQN;IT  ,IU+1,IV  ]                C
C                     K6 -> E[LQN+1,MQN;IT  ,IU  ,IV+1]                C
C                     K7 -> E[LQN+1,MQN;IT-1,IU  ,IV  ]                C
C                     K8 -> E[LQN+1,MQN;IT  ,IU-1,IV  ]                C
C                     K9 -> E[LQN+1,MQN;IT  ,IU  ,IV-1]                C
C                     K10-> E[LQN+1,MQN;IT-2,IU  ,IV  ]                C
C                     K11-> E[LQN+1,MQN;IT  ,IU-2,IV  ]                C
C                     K12-> E[LQN+1,MQN;IT  ,IU  ,IV-2]                C
C**********************************************************************C
C
      DO IOUTER=0,LAM
        DO IT=0,IOUTER
          DO IU=0,IOUTER-IT
            IV = IOUTER-IT-IU
C
C           STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (NQN)
            I0  = IABC(IT  ,IU  ,IV  )
            K0  = IABC(IT  ,IU  ,IV  )
            K1  = IABC(IT+2,IU  ,IV  )
            K2  = IABC(IT  ,IU+2,IV  )
            K3  = IABC(IT  ,IU  ,IV+2)
            K4  = IABC(IT+1,IU  ,IV  )
            K5  = IABC(IT  ,IU+1,IV  )
            K6  = IABC(IT  ,IU  ,IV+1)
            KN  = IABC(IT+MX,IU+MY,IV+MZ)
            IF(IT.GT.0) K7  = IABC(IT-1,IU  ,IV  )
            IF(IU.GT.0) K8  = IABC(IT  ,IU-1,IV  )
            IF(IV.GT.0) K9  = IABC(IT  ,IU  ,IV-1)
            IF(IT.GT.1) K10 = IABC(IT-2,IU  ,IV  )
            IF(IU.GT.1) K11 = IABC(IT  ,IU-2,IV  )
            IF(IV.GT.1) K12 = IABC(IT  ,IU  ,IV-2)
            IF(IT-MX.GE.0.AND.IU-MY.GE.0.AND.IV-MZ.GE.0) THEN
              KM = IABC(IT-MX,IU-MY,IV-MZ)
            ENDIF
C
C           RECURRENCE RELATIONS ON LAYER (LQN+1) OF E'S AND G'S
C
C           TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
            TT = DFLOAT((2*(IT+IU+IV))+3)
            DO M=1,MAXM
              T0 = 1.0D0/P(M)
              T1 = 1.0D0/P22(M)
              TX = PX(M)/P(M)
              TY = PY(M)/P(M)
              TZ = PZ(M)/P(M)
              TP = PP(M) + TT/P2(M)
              TN = 2.0D0*(MX*PX(M) + MY*PY(M) + MZ*PZ(M))
              ENSG(M,K0) = ENSG(M,K0) + TP*ESG(M,I0)
              GNSG(M,K0) = GNSG(M,K0) + TP*GSG(M,I0)
              ENSG(M,K1) = ENSG(M,K1) + T1*ESG(M,I0)
              GNSG(M,K1) = GNSG(M,K1) + T1*GSG(M,I0)
              ENSG(M,K2) = ENSG(M,K2) + T1*ESG(M,I0)
              GNSG(M,K2) = GNSG(M,K2) + T1*GSG(M,I0)
              ENSG(M,K3) = ENSG(M,K3) + T1*ESG(M,I0)
              GNSG(M,K3) = GNSG(M,K3) + T1*GSG(M,I0)
              ENSG(M,K4) = ENSG(M,K4) + TX*ESG(M,I0)
              GNSG(M,K4) = GNSG(M,K4) + TX*GSG(M,I0)
              ENSG(M,K5) = ENSG(M,K5) + TY*ESG(M,I0)
              GNSG(M,K5) = GNSG(M,K5) + TY*GSG(M,I0)
              ENSG(M,K6) = ENSG(M,K6) + TZ*ESG(M,I0)
              GNSG(M,K6) = GNSG(M,K6) + TZ*GSG(M,I0)
            ENDDO
C
C           SPECIAL CASE EXCLUDES IT=0
            IF(IT.GE.1) THEN
              RT2 = DFLOAT(2*IT)
              DO M=1,MAXM
                T0 = PX(M)*RT2
                ENSG(M,K7) = ENSG(M,K7) + T0*ESG(M,I0)
                GNSG(M,K7) = GNSG(M,K7) + T0*GSG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IU=0
            IF(IU.GE.1) THEN
              RU1 = DFLOAT(2*IU)
              DO M=1,MAXM
                T0 = PY(M)*RU1
                ENSG(M,K8) = ENSG(M,K8) + T0*ESG(M,I0)
                GNSG(M,K8) = GNSG(M,K8) + T0*GSG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IV=0
            IF(IV.GE.1) THEN
              RV1 = DFLOAT(2*IV)
              DO M=1,MAXM
                T0 = PZ(M)*RV1
                ENSG(M,K9) = ENSG(M,K9) + T0*ESG(M,I0)
                GNSG(M,K9) = GNSG(M,K9) + T0*GSG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IT=0,1
            IF(IT.GE.2) THEN
              RT2 = DFLOAT(IT*(IT-1))
              DO M=1,MAXM
                ENSG(M,K10) = ENSG(M,K10) + RT2*ESG(M,I0)
                GNSG(M,K10) = GNSG(M,K10) + RT2*GSG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IU=0,1
            IF(IU.GE.2) THEN
              RU2 = DFLOAT(IU*(IU-1))
              DO M=1,MAXM
                ENSG(M,K11) = ENSG(M,K11) + RU2*ESG(M,I0)
                GNSG(M,K11) = GNSG(M,K11) + RU2*GSG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IV=0,1
            IF(IV.GE.2) THEN
              RV2 = DFLOAT(IV*(IV-1))
              DO M=1,MAXM
                ENSG(M,K12) = ENSG(M,K12) + RV2*ESG(M,I0)
                GNSG(M,K12) = GNSG(M,K12) + RV2*GSG(M,I0)
              ENDDO
            ENDIF
C
C           RECURRENCE RELATIONS ON G'S SOMETIMES HAVE CROSSING TERMS
            IF(IZ.EQ.1.AND.LR.EQ.'R') GOTO 100
            IF(IZ.EQ.2.AND.LR.EQ.'L') GOTO 100
C
C           TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
            DO M=1,MAXM
              TN = 2.0D0*(MX*PX(M) + MY*PY(M) + MZ*PZ(M))
              T0 = 1.0D0/P(M)
              GNSG(M,KN) = GNSG(M,KN) - T0*ESG(M,I0)
              GNSG(M,K0) = GNSG(M,K0) - TN*ESG(M,I0)
            ENDDO
C
C           SPECIAL CASE EXCLUDES NX=0
            IF(IT-MX.GE.0.AND.IU-MY.GE.0.AND.IV-MZ.GE.0) THEN
              TM = 2.0D0*DFLOAT(MX*IT + MY*IU + MZ*IV)
              DO M=1,MAXM
                GNSG(M,KM) = GNSG(M,KM) - TM*ESG(M,I0)
              ENDDO
            ENDIF
C
100         CONTINUE
C
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE RNLL(RN,EXL,LQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 RRRRRRR  NN    NN LL       LL                        C
C                 RR    RR NNN   NN LL       LL                        C
C                 RR    RR NNNN  NN LL       LL                        C
C                 RR    RR NN NN NN LL       LL                        C
C                 RRRRRRR  NN  NNNN LL       LL                        C
C                 RR    RR NN   NNN LL       LL                        C
C                 RR    RR NN    NN LLLLLLLL LLLLLLLL                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  RNLL GENERATES THE LARGE-LARGE SGTF NORMALISATION CONSTANTS.        C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ LQN  - PAIR OF ORBITAL QUANTUM NUMBERS.                           C
C  ▶ NBAS - PAIR OF PARAMETER LIST LENGTHS.                            C
C  OUTPUT:                                                             C
C  ▶ RN   - BATCH OF NORMALISATION CONSTANTS.                          C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RN(MBS,2),EXL(MBS,2),LQN(2),NBAS(2)
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
      DO I=1,2
        T1  = PI12
        F1  = 0.5D0
        GML = DLOG(T1)
        DO N=2,LQN(I)+2
          GML = GML+DLOG(F1)
          F1  = F1 + 1.0D0
        ENDDO
        RLA = DFLOAT(LQN(I))
        GA1 = TWLG-GML
        RA1 = RLA+1.5D0
        DO IBAS=1,NBAS(I)
          ELOG       = DLOG(2.0D0*EXL(IBAS,I))
          RN(IBAS,I) = DEXP(0.5D0*(GA1+RA1*ELOG))
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE RNSS(RN,EXL,LQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 RRRRRRR  NN    NN  SSSSSS   SSSSSS                   C
C                 RR    RR NNN   NN SS    SS SS    SS                  C
C                 RR    RR NNNN  NN SS       SS                        C
C                 RR    RR NN NN NN  SSSSSS   SSSSSS                   C
C                 RRRRRRR  NN  NNNN       SS       SS                  C
C                 RR    RR NN   NNN SS    SS SS    SS                  C
C                 RR    RR NN    NN  SSSSSS   SSSSSS                   C
C                                                                      C
C -------------------------------------------------------------------- C
C  RNSS GENERATES THE SMALL-SMALL SGTF NORMALISATION CONSTANTS.        C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ LQN  - PAIR OF ORBITAL QUANTUM NUMBERS.                           C
C  ▶ NBAS - PAIR OF PARAMETER LIST LENGTHS.                            C
C  OUTPUT:                                                             C
C  ▶ RN - BATCH OF NORMALISATION CONSTANTS.                            C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RN(MBS,2),EXL(MBS,2),LQN(2),NBAS(2)
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
      DO I=1,2
        T1  = PI12
        F1  = 0.5D0
        GML = DLOG(T1)
        DO N=2,LQN(I)+3
          GML = GML+DLOG(F1)
          F1  = F1+1.0D0
        ENDDO
        RLA = DFLOAT(LQN(I))
        GA1 = TWLG-GML
        RA1 = RLA+0.5D0
        DO IBAS=1,NBAS(I)
          ELOG       = DLOG(2.0D0*EXL(IBAS,I))
          RN(IBAS,I) = DEXP(0.5D0*(GA1+RA1*ELOG))
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE RNLS(RN,EXL,LQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 RRRRRRR  NN    NN LL       SSSSSS                    C
C                 RR    RR NNN   NN LL      SS    SS                   C
C                 RR    RR NNNN  NN LL      SS                         C
C                 RR    RR NN NN NN LL       SSSSSS                    C
C                 RRRRRRR  NN  NNNN LL            SS                   C
C                 RR    RR NN   NNN LL      SS    SS                   C
C                 RR    RR NN    NN LLLLLLLL SSSSSS                    C
C                                                                      C
C -------------------------------------------------------------------- C
C  RNLS GENERATES THE LARGE-SMALL SGTF NORMALISATION CONSTANTS.        C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ LQN  - PAIR OF ORBITAL QUANTUM NUMBERS.                           C
C  ▶ NBAS - PAIR OF PARAMETER LIST LENGTHS.                            C
C  OUTPUT:                                                             C
C  ▶ RN   - BATCH OF NORMALISATION CONSTANTS.                          C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RN(MBS,2),EXL(MBS,2),LQN(2),NBAS(2)
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
      T1L  = PI12
      F1L  = 0.5D0
      GMLL = DLOG(T1L)
C
      DO N=2,LQN(1)+2
        GMLL = GMLL+DLOG(F1L)
        F1L  = F1L+1.0D0
      ENDDO
C
      RLAL = DFLOAT(LQN(1))
      GA1L = TWLG-GMLL
      RA1L = RLAL+1.5D0
C
      DO IBAS=1,NBAS(1)
        ELOGL      = DLOG(2.0D0*EXL(IBAS,1))
        RN(IBAS,1) = DEXP(0.5D0*(GA1L+RA1L*ELOGL))
      ENDDO
C
      T1S  = PI12
      F1S  = 0.5D0
      GMLS = DLOG(T1S)
C
      DO N=2,LQN(2)+3
        GMLS = GMLS+DLOG(F1S)
        F1S  = F1S+1.0D0
      ENDDO
C
      RLAS = DFLOAT(LQN(2))
      GA1S = TWLG-GMLS
      RA1S = RLAS+0.5D0
C
      DO JBAS=1,NBAS(2)
        ELOGS      = DLOG(2.0D0*EXL(JBAS,2))
        RN(JBAS,2) = DEXP(0.5D0*(GA1S+RA1S*ELOGS))
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE RNSL(RN,EXL,LQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 RRRRRRR  NN    NN  SSSSSS  LL                        C
C                 RR    RR NNN   NN SS    SS LL                        C
C                 RR    RR NNNN  NN SS       LL                        C
C                 RR    RR NN NN NN  SSSSSS  LL                        C
C                 RRRRRRR  NN  NNNN       SS LL                        C
C                 RR    RR NN   NNN SS    SS LL                        C
C                 RR    RR NN    NN  SSSSSS  LLLLLLLL                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  RNSL GENERATES THE SMALL-LARGE SGTF NORMALISATION CONSTANTS.        C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ LQN  - PAIR OF ORBITAL QUANTUM NUMBERS.                           C
C  ▶ NBAS - PAIR OF PARAMETER LIST LENGTHS.                            C
C  OUTPUT:                                                             C
C  ▶ RN   - BATCH OF NORMALISATION CONSTANTS.                          C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RN(MBS,2),EXL(MBS,2),LQN(2),NBAS(2)
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
      T1S  = PI12
      F1S  = 0.5D0
      GMSL = DLOG(T1S)
C
      DO N=2,LQN(1)+3
        GMSL = GMSL+DLOG(F1S)
        F1S  = F1S+1.0D0
      ENDDO
C
      RLAS = DFLOAT(LQN(1))
      GA1S = TWLG-GMSL
      RA1S = RLAS+0.5D0
C
      DO IBAS=1,NBAS(1)
        ELOGS      = DLOG(2.0D0*EXL(IBAS,1))
        RN(IBAS,1) = DEXP(0.5D0*(GA1S+RA1S*ELOGS))
      ENDDO
C
      T1L  = PI12
      F1L  = 0.5D0
      GMLL = DLOG(T1L)
C
      DO N=2,LQN(2)+2
        GMLL = GMLL+DLOG(F1L)
        F1L  = F1L+1.0D0
      ENDDO
C
      RLAL = DFLOAT(LQN(2))
      GA1L = TWLG-GMLL
      RA1L = RLAL+1.5D0
C
      DO JBAS=1,NBAS(2)
        ELOGL      = DLOG(2.0D0*EXL(JBAS,2))
        RN(JBAS,2) = DEXP(0.5D0*(GA1L+RA1L*ELOGL))
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE DNORM(NMAX,ECFF,ICMP,SCL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           DDDDDDD  NN    NN  OOOOOO  RRRRRRR  MM       MM            C
C           DD    DD NNN   NN OO    OO RR    RR MMM     MMM            C
C           DD    DD NNNN  NN OO    OO RR    RR MMMM   MMMM            C
C           DD    DD NN NN NN OO    OO RR    RR MM MM MM MM            C
C           DD    DD NN  NNNN OO    OO RRRRRRR  MM  MMM  MM            C
C           DD    DD NN   NNN OO    OO RR    RR MM   M   MM            C
C           DDDDDDD  NN    NN  OOOOOO  RR    RR MM       MM            C
C                                                                      C
C -------------------------------------------------------------------- C
C  DNORM CALCULATES A SCALE NORM FOR A REAL OR COMPLEX PART OF A LIST  C
C  ECFF OF LENGTH NMAX, AND STORES THE RESULT IN SCL.                  C
C**********************************************************************C
C
      DIMENSION ECMP(NMAX)
C
      COMPLEX*16 ECFF(NMAX)
C
C     IMPORT EITHER THE REAL OR COMPLEX COMPONENT FROM ECFF
      DO N=1,NMAX
        IF(ICMP.EQ.1) THEN
          ECMP(N) = DREAL(ECFF(N))
        ELSEIF(ICMP.EQ.2) THEN
          ECMP(N) = DIMAG(ECFF(N))
        ELSE
          WRITE(6, *) 'In DNORM: choose component 1 or 2.'
          WRITE(7, *) 'In DNORM: choose component 1 or 2.'
        ENDIF
      ENDDO
C
C     LOOP OVER ELEMENTS OF ECMP
      SSQ = 1.0D0
      SCL = 0.0D0
      DO N=1,NMAX
        IF(ECMP(N).NE.0.0D0) THEN
          ABN = DABS(ECMP(N))
          IF(SCL.LT.ABN) THEN
            SSQ = 1.0D0 + SSQ*(SCL/ABN)**2
          ELSE
            SSQ = SSQ   +     (ABN/SCL)**2
          ENDIF
        ENDIF
      ENDDO
      SCL = SCL*DSQRT(SSQ)
C
      RETURN
      END
