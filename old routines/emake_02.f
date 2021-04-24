C**********************************************************************C
C ==================================================================== C
C  [13] E-COEFFS: FINITE BASIS OVERLAP SPIN STRUCTURE FACTORS.         C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] EQFILE: MAIN ROUTINE FOR BUILDING A GLOBAL FILE OF EQ-COEFFS.  C
C   [B] ESETLL: CONSTRUCT AND SAVE ALL ELL0 COEFFICIENTS EXTERNALLY.   C
C   [C] ESETSS: CONSTRUCT AND SAVE ALL ESS0 COEFFICIENTS EXTERNALLY.   C
C   [D] ESETLS: CONSTRUCT AND SAVE ALL ELSQ COEFFICIENTS EXTERNALLY.   C
C   [E] EMAKELL: GENERATE A COMPLETE BLOCK OF EQLL COEFFICIENTS.       C
C   [F] EMAKESS: GENERATE A COMPLETE BLOCK OF EQSS COEFFICIENTS.       C
C   [G] EMAKELS: GENERATE A COMPLETE BLOCK OF EQLS COEFFICIENTS.       C
C   [H] EQLL: A RAW BLOCK OF EQLL COEFFICIENTS FOR EMAKELL.            C
C   [I] EQSS: A RAW BLOCK OF EQSS COEFFICIENTS FOR EMAKESS.            C
C   [J] EQLS: A RAW BLOCK OF EQLS COEFFICIENTS FOR EMAKELS.            C
C   [K] ESGTF: SET OF ES-COEFFS OVER SPHERICAL HARMONICS AND HGTFS.    C
C   [L] VRS: EXPANSION COEFFS IN HGTF OVERLAPS, CALLED IN ESGTF.       C
C   [M] STEPLM: SIMULTANEOUS INCREASE IN (L,M) FOR USE IN VRS.         C
C   [N] STEPL: INCREMENT IN L FOR USE IN VRS.                          C
C   [O] STEPN: INCREMENT IN N FOR USE IN VRS.                          C
C   [P] RNLL: A BLOCK OF LL NORMALISATION COEFFS.                      C
C   [Q] RNSS: A BLOCK OF SS NORMALISATION COEFFS.                      C
C   [R] RNLS: A BLOCK OF LS NORMALISATION COEFFS.                      C
C   [S] DNORM: NORM FOR A REAL OR COMPLEX PART OF EQ-COEFF LIST.       C
C**********************************************************************C
C
C
      SUBROUTINE EQFILE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          EEEEEEEE  QQQQQQ   FFFFFFFF IIII LL       EEEEEEEE          C
C          EE       QQ    QQ  FF        II  LL       EE                C
C          EE       QQ    QQ  FF        II  LL       EE                C
C          EEEEEE   QQ    QQ  FFFFFF    II  LL       EEEEEE            C
C          EE       QQ   QQQ  FF        II  LL       EE                C
C          EE       QQ    QQ  FF        II  LL       EE                C
C          EEEEEEEE  QQQQQQ Q FF       IIII LLLLLLLL EEEEEEEE          C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQFILE CONSTRUCTS A SET OF COMMON ARRAYS FOR ALL REQUIRED EQTT      C
C  COEFFICIENTS IN A CALCULATION THAT RESTS WITHIN QED.                C
C**********************************************************************C
      PARAMETER(MBS=26,MCT=6,MKP=9,MFL=10000000)
C
      CHARACTER*4 HMLTN
C
      DIMENSION NLL(0:MKP),NSS(0:MKP),NLS(0:MKP)
C
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
      COMMON/TMMD/TELL,TESS,TELS,TESL,TRLL,TRSS,TRLS,TRBR
C
C**********************************************************************C
C     LOOP OVER FOCK BLOCK AND COUNT ALL REQUIRED EQ-WORDS.            C
C**********************************************************************C
C
C     INITIALISE MAXIMUM LAMBDA
      LAMMX = 0
C
C     INITIALISE TOTAL COEFFICIENT COUNTERS
      NADE0LL =  0
      NADE0SS =  0
      NADEILS =  0
C
C     INITIALISE COEFFICIENT COUNTERS FOR LAMBDA CLASS
      DO ILAM=0,MKP
        NLL(ILAM) = 0
        NSS(ILAM) = 0
        NLS(ILAM) = 0
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
            IF(KVALS(KA,ICNTA).GT.0) THEN
              LQNA = KVALS(KA,ICNTA)
            ELSE
              LQNA =-KVALS(KA,ICNTA)-1
            ENDIF
            NBASA = NFUNCT(LQNA+1,ICNTA)
C
C           LOOP OVER KQNB VALUES
            DO KB=1,NKAP(ICNTB)
C
C             QUANTUM NUMBERS FOR BLOCK B
              IF(KVALS(KB,ICNTB).GT.0) THEN
                LQNB = KVALS(KB,ICNTB)
              ELSE
                LQNB =-KVALS(KB,ICNTB)-1
              ENDIF
              NBASB = NFUNCT(LQNB+1,ICNTB)
C
C             LOOP OVER |MQNA| VALUES
              DO MA=1,IABS(KVALS(KA,ICNTA))
                MJA = 2*MA-1
C
C               LOOP OVER |MQNB| VALUES
                DO MB=1,IABS(KVALS(KB,ICNTB))
                  MJB = 2*MB-1
C
C                 NUMBER OF BASIS FUNCTION OVERLAPS
                  MAXAB = NBASA*NBASB
C
C                 CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
                  LAMAB  = LQNA+LQNB
                  NTUVLL = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
                  NTUVSS = (LAMAB+3)*(LAMAB+4)*(LAMAB+5)/6
                  NTUVLS = (LAMAB+2)*(LAMAB+3)*(LAMAB+4)/6
C
C                 UPDATE LARGEST LAMBDA VALUE
                  IF(LAMAB.GT.LAMMX) THEN
                    LAMMX = LAMAB
                  ENDIF
C
C                 INCREASE NUMBER OF WORDS FOR THIS LAMBDA VALUE
                  NLL(LAMAB  ) = NLL(LAMAB  ) + NTUVLL*MAXAB
                  NSS(LAMAB+2) = NSS(LAMAB+2) + NTUVSS*MAXAB
                  NLS(LAMAB+1) = NLS(LAMAB+1) + NTUVLS*MAXAB*3
C
C                 INCREASE TOTAL NUMBER OF WORDS
                  NADE0LL =  NADE0LL + NTUVLL*MAXAB
                  NADE0SS =  NADE0SS + NTUVSS*MAXAB
                  NADEILS =  NADEILS + NTUVLS*MAXAB*3
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
C     DOUBLE LOOP OVER FOCK BLOCK COMPLETE
C
C**********************************************************************C
C     SUMMARY OF WORD ANALYSIS                                         C
C**********************************************************************C
C
C     SECTION TITLE
20    FORMAT(1X,A,4X,A,4X,A,6X,A,4X,A,8X,A,6X,A)
21    FORMAT(1X,A,8X,I2,4X,I5,3X,I9,7X,I2,3X,I10,5X,F10.3)
22    FORMAT(1X,A,20X,I10,12X,I10,5X,F10.3)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',22),'E-coefficient word analysis'
      WRITE(7, *) REPEAT(' ',22),'E-coefficient word analysis'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,20) 'Type','Lambda','Terms','Length','Mult.',
     &            'Words','Size (MB)'
      WRITE(7,20) 'Type','Lambda','Terms','Length','Mult.',
     &            'Words','Size (MB)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     INITIALISE OVERALL WORD COUNTER AND SIZE COUNTER
      NADETOT = 0
      NWRDTOT = 0
      SPCETOT = 0.0D0
C
C     E0LL ANALYSIS
      DO ILAM=0,LAMMX
        IF(NLL(ILAM).EQ.0) GOTO 200
        NTUVLL = (ILAM+1)*(ILAM+2)*(ILAM+3)/6
        SPCELL = NLL(ILAM)*8*8.0D-6
        WRITE(6,21) 'E0LL',ILAM,NTUVLL,NLL(ILAM),8,NLL(ILAM)*8,SPCELL
        WRITE(7,21) 'E0LL',ILAM,NTUVLL,NLL(ILAM),8,NLL(ILAM)*8,SPCELL
        NADETOT = NADETOT + NLL(ILAM)
        NWRDTOT = NWRDTOT + NLL(ILAM)*8
        SPCETOT = SPCETOT + SPCELL
200     CONTINUE
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
      IF(HMLTN.EQ.'NORL') GOTO 100
C
C     E0SS ANALYSIS
      DO ILAM=2,LAMMX+2
        IF(NSS(ILAM).EQ.0) GOTO 210
        NTUVSS = (ILAM+1)*(ILAM+2)*(ILAM+3)/6
        SPCESS = NSS(ILAM)*8*8.0D-6
        WRITE(6,21) 'E0SS',ILAM,NTUVSS,NSS(ILAM),8,NSS(ILAM)*8,SPCESS
        WRITE(7,21) 'E0SS',ILAM,NTUVSS,NSS(ILAM),8,NSS(ILAM)*8,SPCESS
        NADETOT = NADETOT + NSS(ILAM)
        NWRDTOT = NWRDTOT + NSS(ILAM)*8
        SPCETOT = SPCETOT + SPCESS
210     CONTINUE
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
      IF(HMLTN.EQ.'BARE'.OR.HMLTN.EQ.'DHFR') GOTO 100
C
C     EILS ANALYSIS
      DO ILAM=1,LAMMX+1
        IF(NLS(ILAM).EQ.0) GOTO 220
        NTUVLS = (ILAM+1)*(ILAM+2)*(ILAM+3)/6
        SPCELS = NSS(ILAM)*24*8.0D-6
        WRITE(6,21) 'EILS',ILAM,NTUVLS,NLS(ILAM),24,NLS(ILAM)*24,SPCELS
        WRITE(7,21) 'EILS',ILAM,NTUVLS,NLS(ILAM),24,NLS(ILAM)*24,SPCELS
        NADETOT = NADETOT + NLS(ILAM)
        NWRDTOT = NWRDTOT + NLS(ILAM)*24
        SPCETOT = SPCETOT + SPCELS
220     CONTINUE
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
100   CONTINUE
C
C     SUMMARY OF TOTALS
      WRITE(6,22) 'Total',NADETOT,NWRDTOT,SPCETOT
      WRITE(7,22) 'Total',NADETOT,NWRDTOT,SPCETOT
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     OPTION WHEN NUMBER OF WORDS EXCEEDS ALLOCATED SIZE LIMIT
      IF(NADE0LL.GT.MFL) THEN
        WRITE(6, *) 'In EQFILE: E0LL words exceed allocated limit.'
        WRITE(7, *) 'In EQFILE: E0LL words exceed allocated limit.'
        GOTO 150
      ENDIF
C
      IF(HMLTN.NE.'NORL') THEN
        IF(NADE0SS.GT.MFL) THEN
          WRITE(6, *) 'In EQFILE: E0SS words exceed allocated limit.'
          WRITE(7, *) 'In EQFILE: E0SS words exceed allocated limit.'
          GOTO 150
        ENDIF
      ENDIF
C
      IF(HMLTN.EQ.'DHFB'.OR.HMLTN.EQ.'DHFP'.OR.HMLTN.EQ.'DHFQ') THEN
        IF(NADEILS.GT.MFL) THEN
          WRITE(6, *) 'In EQFILE: EILS words exceed allocated limit.'
          WRITE(7, *) 'In EQFILE: EILS words exceed allocated limit.'
          GOTO 150
        ENDIF
      ENDIF

C     SIZE LIMITS ARE ALL OK -- SKIP TO BATCH GENERATION
      GOTO 250
C
C     ONE OF THE CLASSES EXCEEDS WORD LIMIT
150   CONTINUE
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     HAVE TO GENERATE E COEFFICIENTS BY BATCH
      WRITE(6, *) 'In EQFILE: E-coefficients to be generated by batch.'
      WRITE(7, *) 'In EQFILE: E-coefficients to be generated by batch.'
C
C     FLIP THE EQ-GENERATION TOGGLE AND EXIT
      IEQS = 0
      GOTO 300
C
250   CONTINUE
C
C**********************************************************************C
C     GENERATE COMPLETE BATCH OF EQ-COEFFS                             C
C**********************************************************************C
C
C     SECTION TITLE
      WRITE(6, *) REPEAT(' ',18),'Generating E-coefficient data files'
      WRITE(7, *) REPEAT(' ',18),'Generating E-coefficient data files'
C
C     E0LL COEFFICIENTS
      CALL CPU_TIME(TDM1)
      CALL ESETLL
      CALL CPU_TIME(TDM2)
      TELL = TELL+TDM2-TDM1
C
      IF(HMLTN.EQ.'NORL') GOTO 300
C
C     E0SS COEFFICIENTS
      CALL ESETSS
      CALL CPU_TIME(TDM3)
      TESS = TELL+TDM3-TDM2
C
      IF(HMLTN.EQ.'BARE'.OR.HMLTN.EQ.'DHFR') GOTO 300
C
C     EILS COEFFICIENTS
      CALL ESETLS
      CALL CPU_TIME(TDM4)
      TELS = TELS+TDM4-TDM3
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
      SUBROUTINE ESETLL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         EEEEEEEE SSSSSS  EEEEEEEE TTTTTTTT LL       LL               C
C         EE      SS    SS EE          TT    LL       LL               C
C         EE      SS       EE          TT    LL       LL               C
C         EEEEEE   SSSSSS  EEEEEE      TT    LL       LL               C
C         EE            SS EE          TT    LL       LL               C
C         EE      SS    SS EE          TT    LL       LL               C
C         EEEEEEEE SSSSSS  EEEEEEEE    TT    LLLLLLLL LLLLLLLL         C
C                                                                      C
C -------------------------------------------------------------------- C
C  ESETLL CONSTRUCTS ALL ELL0 COEFFICIENTS FOR A SYSTEM WITH CALLS TO  C
C  EMAKELL AND SAVES THEM TO AN EXTERNAL DATA FILE.                    C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
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
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFUNCT(LQN(1)+1,ICNTA)
        DO IBAS=1,NBAS(1)
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFUNCT(LQN(2)+1,ICNTB)
        DO JBAS=1,NBAS(2)
          EXPT(JBAS,2) = EXPSET(JBAS,LQN(2)+1,ICNTB)
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
      CALL EMAKELL(E11,E21,EXPT,XYZ,KQN,MQN,NBAS,+1,1,2,0)
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
C     GENERATE ELL0(CD) COEFFICIENTS
      CALL EMAKELL(E11,E21,EXPT,XYZ,KQN,MQN,NBAS,-1,1,2,0)
C
C     WRITE ELL0(CD) TO THE MAIN ARRAY
      DO ITUV=1,NTUVLL
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          E0LLFL(IAD+M,5) = DREAL(E11(M,ITUV))
          E0LLFL(IAD+M,6) = DIMAG(E11(M,ITUV))
          E0LLFL(IAD+M,7) = DREAL(E21(M,ITUV))
          E0LLFL(IAD+M,8) = DIMAG(E21(M,ITUV))
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
      SUBROUTINE ESETSS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         EEEEEEEE SSSSSS  EEEEEEEE TTTTTTTT SSSSSS   SSSSSS           C
C         EE      SS    SS EE          TT   SS    SS SS    SS          C
C         EE      SS       EE          TT   SS       SS                C
C         EEEEEE   SSSSSS  EEEEEE      TT    SSSSSS   SSSSSS           C
C         EE            SS EE          TT         SS       SS          C
C         EE      SS    SS EE          TT   SS    SS SS    SS          C
C         EEEEEEEE SSSSSS  EEEEEEEE    TT    SSSSSS   SSSSSS           C
C                                                                      C
C -------------------------------------------------------------------- C
C  ESETSS CONSTRUCTS ALL ESS0 COEFFICIENTS FOR A SYSTEM WITH CALLS TO  C
C  EMAKESS AND SAVES THEM TO AN EXTERNAL DATA FILE.                    C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
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
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFUNCT(LQN(1)+1,ICNTA)
        DO IBAS=1,NBAS(1)
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFUNCT(LQN(2)+1,ICNTB)
        DO JBAS=1,NBAS(2)
          EXPT(JBAS,2) = EXPSET(JBAS,LQN(2)+1,ICNTB)
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
      CALL EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NBAS,+1,1,2,0)
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
C     GENERATE ESS0(CD) COEFFICIENTS
      CALL EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NBAS,-1,1,2,0)
C
C     WRITE ESS0(CD) TO THE MAIN ARRAY
      DO ITUV=1,NTUVSS
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          E0SSFL(IAD+M,5) = DREAL(E11(M,ITUV))
          E0SSFL(IAD+M,6) = DIMAG(E11(M,ITUV))
          E0SSFL(IAD+M,7) = DREAL(E21(M,ITUV))
          E0SSFL(IAD+M,8) = DIMAG(E21(M,ITUV))
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
      SUBROUTINE ESETLS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         EEEEEEEE SSSSSS  EEEEEEEE TTTTTTTT LL       SSSSSS           C
C         EE      SS    SS EE          TT    LL      SS    SS          C
C         EE      SS       EE          TT    LL      SS                C
C         EEEEEE   SSSSSS  EEEEEE      TT    LL       SSSSSS           C
C         EE            SS EE          TT    LL            SS          C
C         EE      SS    SS EE          TT    LL      SS    SS          C
C         EEEEEEEE SSSSSS  EEEEEEEE    TT    LLLLLLLL SSSSSS           C
C                                                                      C
C -------------------------------------------------------------------- C
C  ESETLS CONSTRUCTS ALL ESS0 COEFFICIENTS FOR A SYSTEM WITH CALLS TO  C
C  EMAKELS AND SAVES THEM TO AN EXTERNAL DATA FILE.                    C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
C
      COMPLEX*16 E11(MB2,MEQ,3),E21(MB2,MEQ,3)
C
      COMMON/EILS/EILSFL(MFL,24),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
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
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFUNCT(LQN(1)+1,ICNTA)
        DO IBAS=1,NBAS(1)
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFUNCT(LQN(2)+1,ICNTB)
        DO JBAS=1,NBAS(2)
          EXPT(JBAS,2) = EXPSET(JBAS,LQN(2)+1,ICNTB)
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
C     GENERATE ELSI(AB) COEFFICIENTS
      CALL EMAKEB3(E11,E21,EXPT,XYZ,KQN,MQN,NBAS,+1,1,2)
C
C     WRITE ELSI(AB) TO THE MAIN ARRAY
      DO ITUV=1,NTUVLS
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          EILSFL(IAD+M, 1) = DREAL(E11(M,ITUV,1))
          EILSFL(IAD+M, 2) = DIMAG(E11(M,ITUV,1))
          EILSFL(IAD+M, 3) = DREAL(E21(M,ITUV,1))
          EILSFL(IAD+M, 4) = DIMAG(E21(M,ITUV,1))
          EILSFL(IAD+M, 5) = DREAL(E11(M,ITUV,2))
          EILSFL(IAD+M, 6) = DIMAG(E11(M,ITUV,2))
          EILSFL(IAD+M, 7) = DREAL(E21(M,ITUV,2))
          EILSFL(IAD+M, 8) = DIMAG(E21(M,ITUV,2))
          EILSFL(IAD+M, 9) = DREAL(E11(M,ITUV,3))
          EILSFL(IAD+M,10) = DIMAG(E11(M,ITUV,3))
          EILSFL(IAD+M,11) = DREAL(E21(M,ITUV,3))
          EILSFL(IAD+M,12) = DIMAG(E21(M,ITUV,3))
        ENDDO
      ENDDO
C
C     GENERATE ELSI(CD) COEFFICIENTS
      CALL EMAKEB3(E11,E21,EXPT,XYZ,KQN,MQN,NBAS,-1,1,2)
C
C     WRITE ELSI(CD) TO THE MAIN ARRAY
      DO ITUV=1,NTUVLS
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          EILSFL(IAD+M,13) = DREAL(E11(M,ITUV,1))
          EILSFL(IAD+M,14) = DIMAG(E11(M,ITUV,1))
          EILSFL(IAD+M,15) = DREAL(E21(M,ITUV,1))
          EILSFL(IAD+M,16) = DIMAG(E21(M,ITUV,1))
          EILSFL(IAD+M,17) = DREAL(E11(M,ITUV,2))
          EILSFL(IAD+M,18) = DIMAG(E11(M,ITUV,2))
          EILSFL(IAD+M,19) = DREAL(E21(M,ITUV,2))
          EILSFL(IAD+M,20) = DIMAG(E21(M,ITUV,2))
          EILSFL(IAD+M,21) = DREAL(E11(M,ITUV,3))
          EILSFL(IAD+M,22) = DIMAG(E11(M,ITUV,3))
          EILSFL(IAD+M,23) = DREAL(E21(M,ITUV,3))
          EILSFL(IAD+M,24) = DIMAG(E21(M,ITUV,3))
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
      SUBROUTINE EMAKELL(E11,E21,EXPT,XYZ,KQN,MQN,NBAS,IPHS,I1,I2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   EEEEEEEE MM       MM    AA    KK    KK EEEEEEEE LL      LL         C
C   EE       MMM     MMM   AAAA   KK   KK  EE       LL      LL         C
C   EE       MMMM   MMMM  AA  AA  KK  KK   EE       LL      LL         C
C   EEEEEE   MM MM MM MM AA    AA KKKKK    EEEEEE   LL      LL         C
C   EE       MM  MMM  MM AAAAAAAA KK  KK   EE       LL      LL         C
C   EE       MM   M   MM AA    AA KK   KK  EE       LL      LL         C
C   EEEEEEEE MM       MM AA    AA KK    KK EEEEEEEE LLLLLLL LLLLLLL    C
C                                                                      C
C -------------------------------------------------------------------- C
C  EMAKELL GENERATES BLOCKS OF SPHERICAL SPINOR E-COEFFICIENTS BY      C
C  CONTRACTING ON THE E-COEFFICIENTS FOR SCALAR SGTFS, USING A         C
C  DEVELOPMENT OF THE ALGORITHM OF V.R. SAUNDERS.                      C
C -------------------------------------------------------------------- C
C  H.M.QUINEY THE UNIVERSITY OF MELBOURNE (2008).                      C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4)
      DIMENSION EXL(MBS,2),COORD(3,2),KAP(2),MAG(2),NBS(2)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &             IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        COORD(IX,1) = XYZ(IX,I1)
        COORD(IX,2) = XYZ(IX,I2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KAP(1) = KQN(I1)
      KAP(2) = KQN(I2)
C
C     CALCULATE LQN VALUES
      IF(KAP(1).LT.0) THEN
        LA =-KAP(1)-1
      ELSE
        LA = KAP(1)
      ENDIF
      IF(KAP(2).LT.0) THEN
        LB =-KAP(2)-1
      ELSE
        LB = KAP(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NBS(1) = NBAS(I1)
      DO IBAS=1,NBS(1)
        EXL(IBAS,1) = EXPT(IBAS,I1)
      ENDDO
C
      NBS(2) = NBAS(I2)
      DO JBAS=1,NBS(2)
        EXL(JBAS,2) = EXPT(JBAS,I2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMLL  = LA+LB
      NTUVLL = (LAMLL+1)*(LAMLL+2)*(LAMLL+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR M PAIRS (-|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MAG(1) =-MQN(I1)
      MAG(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLL(E11,EXL,COORD,KAP,MAG,NBS,IQ)
C
C     MULTIPLY ELL0 BY PHASE TERM FOR R-INTEGRALS (SAVES COMP. TIME)
      DO ITUV=1,NTUVLL
        PHS = DFLOAT((IPHS)**(LAMVEC(ITUV)))
        DO M=1,NBS(1)*NBS(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR M PAIRS (+|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MAG(1) = MQN(I1)
      MAG(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLL(E21,EXL,COORD,KAP,MAG,NBS,IQ)
C
C     MULTIPLY ELL0 BY PHASE TERM FOR R-INTEGRALS (SAVES COMP. TIME)
      DO ITUV=1,NTUVLL
        PHS = DFLOAT((IPHS)**(LAMVEC(ITUV)))
        DO M=1,NBS(1)*NBS(2)
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
      SUBROUTINE EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NBAS,IPHS,I1,I2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   EEEEEEEE MM       MM    AA    KK    KK EEEEEEEE SSSSSS   SSSSSS    C
C   EE       MMM     MMM   AAAA   KK   KK  EE      SS    SS SS    SS   C
C   EE       MMMM   MMMM  AA  AA  KK  KK   EE      SS       SS         C
C   EEEEEE   MM MM MM MM AA    AA KKKKK    EEEEEE   SSSSSS   SSSSSS    C
C   EE       MM  MMM  MM AAAAAAAA KK  KK   EE            SS       SS   C
C   EE       MM   M   MM AA    AA KK   KK  EE      SS    SS SS    SS   C
C   EEEEEEEE MM       MM AA    AA KK    KK EEEEEEEE SSSSSS   SSSSSS    C
C                                                                      C
C -------------------------------------------------------------------- C
C  EMAKESS GENERATES BLOCKS OF SPHERICAL SPINOR E-COEFFICIENTS BY      C
C  CONTRACTING ON THE E-COEFFICIENTS FOR SCALAR SGTFS, USING A         C
C  DEVELOPMENT OF THE ALGORITHM OF V.R. SAUNDERS.                      C
C -------------------------------------------------------------------- C
C  H.M.QUINEY THE UNIVERSITY OF MELBOURNE (2008).                      C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4)
      DIMENSION EXL(MBS,2),COORD(3,2),KAP(2),MAG(2),NBS(2)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &             IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        COORD(IX,1) = XYZ(IX,I1)
        COORD(IX,2) = XYZ(IX,I2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KAP(1) = KQN(I1)
      KAP(2) = KQN(I2)
C
C     CALCULATE LQN VALUES
      IF(KAP(1).LT.0) THEN
        LA =-KAP(1)-1
      ELSE
        LA = KAP(1)
      ENDIF
      IF(KAP(2).LT.0) THEN
        LB =-KAP(2)-1
      ELSE
        LB = KAP(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NBS(1) = NBAS(I1)
      DO IBAS=1,NBS(1)
        EXL(IBAS,1) = EXPT(IBAS,I1)
      ENDDO
C
      NBS(2) = NBAS(I2)
      DO JBAS=1,NBS(2)
        EXL(JBAS,2) = EXPT(JBAS,I2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMSS  = LA+LB+2
      NTUVSS = (LAMSS+1)*(LAMSS+2)*(LAMSS+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR M PAIRS (-|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MAG(1) =-MQN(I1)
      MAG(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQSS(E11,EXL,COORD,KAP,MAG,NBS,IQ)
C
C     MULTIPLY ESS0 BY PHASE TERM FOR R-INTEGRALS (SAVES COMP. TIME)
      DO ITUV=1,NTUVSS
        PHS = DFLOAT((IPHS)**(LAMVEC(ITUV)))
        DO M=1,NBS(1)*NBS(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR M PAIRS (+|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MAG(1) = MQN(I1)
      MAG(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQSS(E21,EXL,COORD,KAP,MAG,NBS,IQ)
C
C     MULTIPLY ESS0 BY PHASE TERM FOR R-INTEGRALS (SAVES COMP. TIME)
      DO ITUV=1,NTUVSS
        PHS = DFLOAT((IPHS)**(LAMVEC(ITUV)))
        DO M=1,NBS(1)*NBS(2)
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
      SUBROUTINE EMAKELS(E11,E21,EXPT,XYZ,KQN,MQN,NBAS,IPHS,I1,I2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   EEEEEEEE MM       MM    AA    KK    KK EEEEEEEE LL      SSSSSS     C
C   EE       MMM     MMM   AAAA   KK   KK  EE       LL     SS    SS    C
C   EE       MMMM   MMMM  AA  AA  KK  KK   EE       LL     SS          C
C   EEEEEE   MM MM MM MM AA    AA KKKKK    EEEEEE   LL      SSSSSS     C
C   EE       MM  MMM  MM AAAAAAAA KK  KK   EE       LL           SS    C
C   EE       MM   M   MM AA    AA KK   KK  EE       LL     SS    SS    C
C   EEEEEEEE MM       MM AA    AA KK    KK EEEEEEEE LLLLLLL SSSSSS     C
C                                                                      C
C -------------------------------------------------------------------- C
C  EMAKELS GENERATES BLOCKS OF SPHERICAL SPINOR E-COEFFICIENTS BY      C
C  CONTRACTING ON THE E-COEFFICIENTS FOR VECTOR SGTFS, USING A         C
C  DEVELOPMENT OF THE ALGORITHM OF V.R. SAUNDERS.                      C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4)
      DIMENSION EXL(MBS,2),COORD(3,2),KAP(2),MAG(2),NBS(2)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &             IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        COORD(IX,1) = XYZ(IX,I1)
        COORD(IX,2) = XYZ(IX,I2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KAP(1) = KQN(I1)
      KAP(2) = KQN(I2)
C
C     CALCULATE LQN VALUES
      IF(KAP(1).LT.0) THEN
        LA =-KAP(1)-1
      ELSE
        LA = KAP(1)
      ENDIF
      IF(KAP(2).LT.0) THEN
        LB =-KAP(2)-1
      ELSE
        LB = KAP(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NBS(1) = NBAS(I1)
      DO IBAS=1,NBS(1)
        EXL(IBAS,1) = EXPT(IBAS,I1)
      ENDDO
C
      NBS(2) = NBAS(I2)
      DO JBAS=1,NBS(2)
        EXL(JBAS,2) = EXPT(JBAS,I2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMLS  = LA+LB+1
      NTUVLS = (LAMLS+1)*(LAMLS+2)*(LAMLS+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR M PAIRS (-|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MAG(1) =-MQN(I1)
      MAG(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLS(E11,EXL,COORD,KAP,MAG,NBS,IQ)
C
C     MULTIPLY ELSQ BY PHASE TERM FOR R-INTEGRALS (SAVES COMP. TIME)
      DO ITUV=1,NTUVLS
        PHS = DFLOAT((IPHS)**(LAMVEC(ITUV)))
        DO M=1,NBS(1)*NBS(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR M PAIRS (+|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MAG(1) = MQN(I1)
      MAG(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLS(E21,EXL,COORD,KAP,MAG,NBS,IQ)
C
C     MULTIPLY ELSQ BY PHASE TERM FOR R-INTEGRALS (SAVES COMP. TIME)
      DO ITUV=1,NTUVLS
        PHS = DFLOAT((IPHS)**(LAMVEC(ITUV)))
        DO M=1,NBS(1)*NBS(2)
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
      SUBROUTINE EMAKESL(E11,E21,EXPT,XYZ,KQN,MQN,NBAS,IPHS,I1,I2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   EEEEEEEE MM       MM    AA    KK    KK EEEEEEEE SSSSSS  LL         C
C   EE       MMM     MMM   AAAA   KK   KK  EE      SS    SS LL         C
C   EE       MMMM   MMMM  AA  AA  KK  KK   EE      SS       LL         C
C   EEEEEE   MM MM MM MM AA    AA KKKKK    EEEEEE   SSSSSS  LL         C
C   EE       MM  MMM  MM AAAAAAAA KK  KK   EE            SS LL         C
C   EE       MM   M   MM AA    AA KK   KK  EE      SS    SS LL         C
C   EEEEEEEE MM       MM AA    AA KK    KK EEEEEEEE SSSSSS  LLLLLLLL   C
C                                                                      C
C -------------------------------------------------------------------- C
C  EMAKESL GENERATES BLOCKS OF SPHERICAL SPINOR E-COEFFICIENTS BY      C
C  CONTRACTING ON THE E-COEFFICIENTS FOR VECTOR SGTFS, USING A         C
C  DEVELOPMENT OF THE ALGORITHM OF V.R. SAUNDERS.                      C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4)
      DIMENSION EXL(MBS,2),COORD(3,2),KAP(2),MAG(2),NBS(2)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &             IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        COORD(IX,1) = XYZ(IX,I1)
        COORD(IX,2) = XYZ(IX,I2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KAP(1) = KQN(I1)
      KAP(2) = KQN(I2)
C
C     CALCULATE LQN VALUES
      IF(KAP(1).LT.0) THEN
        LA =-KAP(1)-1
      ELSE
        LA = KAP(1)
      ENDIF
      IF(KAP(2).LT.0) THEN
        LB =-KAP(2)-1
      ELSE
        LB = KAP(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NBS(1) = NBAS(I1)
      DO IBAS=1,NBS(1)
        EXL(IBAS,1) = EXPT(IBAS,I1)
      ENDDO
C
      NBS(2) = NBAS(I2)
      DO JBAS=1,NBS(2)
        EXL(JBAS,2) = EXPT(JBAS,I2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMSL  = LA+LB+1
      NTUVSL = (LAMSL+1)*(LAMSL+2)*(LAMSL+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR M PAIRS (-|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MAG(1) =-MQN(I1)
      MAG(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
C      CALL EQSL(E11,EXL,COORD,KAP,MAG,NBS,IQ)
C
C     MULTIPLY ESLQ BY PHASE TERM FOR R-INTEGRALS (SAVES COMP. TIME)
      DO ITUV=1,NTUVSL
        PHS = DFLOAT((IPHS)**(LAMVEC(ITUV)))
        DO M=1,NBS(1)*NBS(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR M PAIRS (+|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MAG(1) = MQN(I1)
      MAG(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
C      CALL EQSL(E21,EXL,COORD,KAP,MAG,NBS,IQ)
C
C     MULTIPLY ESLQ BY PHASE TERM FOR R-INTEGRALS (SAVES COMP. TIME)
      DO ITUV=1,NTUVSL
        PHS = DFLOAT((IPHS)**(LAMVEC(ITUV)))
        DO M=1,NBS(1)*NBS(2)
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
      SUBROUTINE EMAKEB3(E11,E21,EXPT,XYZ,KQN,MQN,NBAS,IPHS,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  EEEEEEEE MM       MM    AA    KK    KK EEEEEEEE BBBBBBB   333333    C
C  EE       MMM     MMM   AAAA   KK   KK  EE       BB    BB 33    33   C
C  EE       MMMM   MMMM  AA  AA  KK  KK   EE       BB    BB       33   C
C  EEEEEE   MM MM MM MM AA    AA KKKKK    EEEEEE   BBBBBBB    33333    C
C  EE       MM  MMM  MM AAAAAAAA KK  KK   EE       BB    BB       33   C
C  EE       MM   M   MM AA    AA KK   KK  EE       BB    BB 33    33   C
C  EEEEEEEE MM       MM AA    AA KK    KK EEEEEEEE BBBBBBB   333333    C
C                                                                      C
C -------------------------------------------------------------------- C
C  EMAKEB3 GENERATES A VECTOR LIST OF ELSQ-COEFFICIENTS FOR A BATCH OF C
C  BREIT INTERACTION INTEGRALS USING A DEVELOPMENT OF THE ALGORITHM OF C
C  V.R. SAUNDERS.                                                      C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      COMPLEX*16 E11(MB2,MEQ,3),E21(MB2,MEQ,3)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4)
      DIMENSION EXL(MBS,2),COORD(3,2),KAP(2),MAG(2),NBS(2)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &             IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        COORD(IX,1) = XYZ(IX,I1)
        COORD(IX,2) = XYZ(IX,I2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KAP(1) = KQN(I1)
      KAP(2) = KQN(I2)
C
C     CALCULATE LQN VALUES
      IF(KAP(1).LT.0) THEN
        LA =-KAP(1)-1
      ELSE
        LA = KAP(1)
      ENDIF
      IF(KAP(2).LT.0) THEN
        LB =-KAP(2)-1
      ELSE
        LB = KAP(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NBS(1) = NBAS(I1)
      DO IBAS=1,NBS(1)
        EXL(IBAS,1) = EXPT(IBAS,I1)
      ENDDO
C
      NBS(2) = NBAS(I2)
      DO JBAS=1,NBS(2)
        EXL(JBAS,2) = EXPT(JBAS,I2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMLS  = LA+LB+1
      NTUVLS = (LAMLS+1)*(LAMLS+2)*(LAMLS+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR M PAIRS (-|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MAG(1) =-MQN(I1)
      MAG(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLS3(E11,EXL,COORD,KAP,MAG,NBS)
C
C     MULTIPLY ELSQ BY PHASE TERM FOR R-INTEGRALS (SAVES COMP. TIME)
      DO ITUV=1,NTUVLS
        PHS = DFLOAT((IPHS)**(LAMVEC(ITUV)))
        DO M=1,NBS(1)*NBS(2)
          DO IQ=1,3
            E11(M,ITUV,IQ) = PHS*E11(M,ITUV,IQ)
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR M PAIRS (+|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MAG(1) = MQN(I1)
      MAG(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLS3(E21,EXL,COORD,KAP,MAG,NBS)
C
C     MULTIPLY ELSQ BY PHASE TERM FOR R-INTEGRALS (SAVES COMP. TIME)
      DO ITUV=1,NTUVLS
        PHS = DFLOAT((IPHS)**(LAMVEC(ITUV)))
        DO M=1,NBS(1)*NBS(2)
          DO IQ=1,3
            E21(M,ITUV,IQ) = PHS*E21(M,ITUV,IQ)
          ENDDO
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
      SUBROUTINE EQLL(ELL,EXL,XYZ,KQN,MQN,NBAS,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                EEEEEEEE  QQQQQQ    LL       LL                       C
C                EE       QQ    QQ   LL       LL                       C
C                EE      QQ      QQ  LL       LL                       C
C                EEEEEE  QQ      QQ  LL       LL                       C
C                EE      QQ      QQ  LL       LL                       C
C                EE       QQ    QQ   LL       LL                       C
C                EEEEEEEE  QQQQQQ QQ LLLLLLLL LLLLLLLL                 C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQLL EVALUATES THE EQ-COEFFICIENTS FOR LARGE-LARGE BASIS FUNCTION   C
C  BLOCK OVERLAPS OF G-SPINOR FUNCTIONS FOR ANY CHOICE IQ = {0,1,2,3}. C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      COMPLEX*16 CONE
      COMPLEX*16 ELL(MB2,MEQ),ESG(MB2,MEQ)
C
      DIMENSION LLAB(2),MLAB(2)
      DIMENSION KQN(2),LQN(2),JQN(2),MQN(2),NBAS(2)
      DIMENSION EXL(MBS,2),RNORM(MBS,2),XYZ(3,2)
C
      COMMON/CTSN/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     THRESHOLD FOR CLEBSCH-GORDAN CEOFFICIENT MAGNITUDES
      SENS = 1.0D-10
C
C     SIGN MULTIPLIER FOR SIGMA COUPLING MATRICES
      IF(IQ.EQ.0.OR.IQ.EQ.1) THEN
        IPHS = 1
      ELSEIF(IQ.EQ.2.OR.IQ.EQ.3) THEN
        IPHS =-1
      ENDIF
C
C     INITIALISE COEFFICIENT STORAGE ARRAY TO ZERO
      DO ITUV=1,MEQ
        DO M=1,MB2
          ELL(M,ITUV) = DCMPLX(0.0D0,0.0D0)
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
      MAXM = NBAS(1)*NBAS(2)
C
C**********************************************************************C
C     INITIALISE COMMON GEOMETRICAL INFORMATION                        C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTRES A AND B
      DX  = (XYZ(1,1)-XYZ(1,2))**2
      DY  = (XYZ(2,1)-XYZ(2,2))**2
      DZ  = (XYZ(3,1)-XYZ(3,2))**2
      AB2 = DX+DY+DZ
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
C
          P(M)    = EXL(IBAS,1)+EXL(JBAS,2)
          P2(M)   = 2.0D0*P(M)
          P22(M)  = P2(M)*P2(M)
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
C     BASIS PAIR LQN LABELS
      LLAB(1) = LQN(1)
      LLAB(2) = LQN(2)
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ELL0 AND ELLZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         CONTRIBUTION TO ELL0 OR ELLZ: UPPER-UPPER
          DO ITUV=1,NTUV
            DO M=1,MAXM
              ELL(M,ITUV) = ELL(M,ITUV) +     CAB*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         CONTRIBUTION TO ELL0 OR ELLZ: LOWER-LOWER
          DO M=1,MAXM
            DO ITUV=1,NTUV
              ELL(M,ITUV) = ELL(M,ITUV) + IPHS*CAB*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ELLX AND ELLY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         CONTRIBUTION TO ELLX OR ELLY: UPPER-LOWER
          DO M=1,MAXM
            DO ITUV=1,NTUV
              ELL(M,ITUV) = ELL(M,ITUV) + IPHS*CAB*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         CONTRIBUTION TO ELLX OR ELLY: LOWER-UPPER
          DO M=1,MAXM
            DO ITUV=1,NTUV
              ELL(M,ITUV) = ELL(M,ITUV) +     CAB*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
C     GENERATE RNLL NORMALISATION CONSTANTS
      CALL RNLL(RNORM,EXL,LQN,NBAS)
C
C     NORMALISE THE ELLQ COEFFICIENT BLOCK
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M   = M+1
          RLL = RNORM(IBAS,1)*RNORM(JBAS,2)
          DO ITUV=1,NTUV
            ELL(M,ITUV) = RLL*ELL(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     SIGMA_Y SPECIAL CASE: MULTIPLY ELLY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV
            ELL(M,ITUV) = CONE*ELL(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE EQSS(ESS,EXL,XYZ,KQN,MQN,NBAS,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 EEEEEEEE  QQQQQQ    SSSSSS   SSSSSS                  C
C                 EE       QQ    QQ  SS    SS SS    SS                 C
C                 EE      QQ      QQ SS       SS                       C
C                 EEEEEE  QQ      QQ  SSSSSS   SSSSSS                  C
C                 EE      QQ      QQ       SS       SS                 C
C                 EE       QQ    QQ  SS    SS SS    SS                 C
C                 EEEEEEEE  QQQQQQ QQ SSSSSS   SSSSSS                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQSS EVALUATES THE EQ-COEFFICIENTS FOR SMALL-SMALL BASIS FUNCTION   C
C  BLOCK OVERLAPS OF G-SPINOR FUNCTIONS FOR ANY CHOICE IQ = {0,1,2,3}. C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      COMPLEX*16 CONE
      COMPLEX*16 ESS(MB2,MEQ),ESG(MB2,MEQ),ENSG(MB2,MEQ)
C
      DIMENSION LLAB(2),MLAB(2)
      DIMENSION KQN(2),LQN(2),JQN(2),MQN(2),NBAS(2)
      DIMENSION EXL(MBS,2),RNORM(MBS,2),XYZ(3,2)
      DIMENSION T22(MB2),T20(MB2),T02(MB2),T00(MB2)
C
      COMMON/CTSN/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     THRESHOLD FOR CLEBSCH-GORDAN CEOFFICIENT MAGNITUDES
      SENS = 1.0D-10
C
C     SIGN MULTIPLIER FOR SIGMA COUPLING MATRICES
      IF(IQ.EQ.0.OR.IQ.EQ.1) THEN
        IPHS = 1
      ELSEIF(IQ.EQ.2.OR.IQ.EQ.3) THEN
        IPHS =-1
      ENDIF
C
C     INITIALISE COEFFICIENT STORAGE ARRAY TO ZERO
      DO M=1,MB2
        DO ITUV=1,MEQ
          ESS(M,ITUV) = DCMPLX(0.0D0,0.0D0)
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
      MAXM = NBAS(1)*NBAS(2)
C
C     KINETIC PRE-FACTORS FOR THIS BLOCK
      RL1 = DFLOAT(2*LQN(1)+1)
      RL2 = DFLOAT(2*LQN(2)+1)
C
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
          T22(M) = 4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
          T20(M) =-2.0D0*RL2*EXL(IBAS,1)
          T02(M) =-2.0D0*RL1*EXL(JBAS,2)
          T00(M) = RL1*RL2
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INITIALISE COMMON GEOMETRICAL INFORMATION                        C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTRES A AND B
      DX  = (XYZ(1,1)-XYZ(1,2))**2
      DY  = (XYZ(2,1)-XYZ(2,2))**2
      DZ  = (XYZ(3,1)-XYZ(3,2))**2
      AB2 = DX+DY+DZ
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
C
          P(M)    = EXL(IBAS,1)+EXL(JBAS,2)
          P2(M)   = 2.0D0*P(M)
          P22(M)  = P2(M)*P2(M)
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
C     BASIS PAIR LQN LABELS
      LLAB(1) = LQN(1)+1
      LLAB(2) = LQN(2)+1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESS0 AND ESSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: UPPER-UPPER
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: LOWER-LOWER
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF

      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ESSX AND ESSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: UPPER-LOWER
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: LOWER-UPPER
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
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
C     BASIS PAIR LQN LABELS
      LLAB(1) = LQN(1)+1
      LLAB(2) = LQN(2)-1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESS0 AND ESSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T20(M) =-CAB*2.0D0*DFLOAT(2*LQN(2)+1)*EXL(IBAS,1)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: UPPER-UPPER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T20(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM,MAXM,2)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)+2
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: UPPER-UPPER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T20(M) =-CAB*2.0D0*DFLOAT(2*LQN(2)+1)*EXL(IBAS,1)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: LOWER-LOWER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T20(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM,MAXM,2)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)+2
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: LOWER-LOWER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ESSX AND ESSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          PREFAC =-CAB*2.0D0*DFLOAT(2*LQN(2)+1)
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T20(M) = PREFAC*EXL(IBAS,1)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: UPPER-LOWER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T20(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM,MAXM,2)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)+1
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: UPPER-LOWER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          PREFAC =-CAB*2.0D0*DFLOAT(2*LQN(2)+1)
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T20(M) = PREFAC*EXL(IBAS,1)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: LOWER-UPPER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T20(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM,MAXM,2)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)+1
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: LOWER-UPPER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
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
C     BASIS PAIR LQN LABELS
      LLAB(1) = LQN(1)-1
      LLAB(2) = LQN(2)+1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESS0 AND ESSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T02(M) =-CAB*2.0D0*DFLOAT(2*LQN(1)+1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: UPPER-UPPER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T02(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM,MAXM,1)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)+2
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: UPPER-UPPER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          PREFAC =-CAB*2.0D0*DFLOAT(2*LQN(1)+1)
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T02(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: LOWER-LOWER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T02(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM,MAXM,1)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)+2
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: LOWER-LOWER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ESSX AND ESSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          PREFAC =-CAB*2.0D0*DFLOAT(2*LQN(1)+1)
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T02(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: UPPER-LOWER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T02(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM,MAXM,1)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)+1
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: UPPER-LOWER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          PREFAC =-CAB*2.0D0*DFLOAT(2*LQN(1)+1)
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T02(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: LOWER-UPPER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T02(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM,MAXM,1)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)+1
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: LOWER-UPPER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
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
C     BASIS PAIR LQN LABELS
      LLAB(1) = LQN(1)-1
      LLAB(2) = LQN(2)-1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESS0 AND ESSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T00(M) = CAB*DFLOAT((2*LQN(1)+1)*(2*LQN(2)+1))
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: UPPER-UPPER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T00(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE FOR THE A CENTRE
          CALL STEPN(ESG,ENSG,LAM,MAXM,1)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION (DFNOTE: CHECK)
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+3)*(LAM+4)*(LAM+5)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          PREFAC =-CAB*2.0D0*DFLOAT(2*LQN(2)+1)
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T20(M) = PREFAC*EXL(IBAS,1)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: UPPER-UPPER FOR [N=1,N'=0]
          DO M=1,MAXM
            TK = T20(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE FOR THE B CENTRE
          CALL STEPN(ESG,ENSG,LAM,MAXM,2)
C
C         KINETIC AND ANGULAR PRE-FACTORS
          PREFAC =-CAB*2.0D0*DFLOAT(2*LQN(1)+1)
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T02(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          NTUV = (LAM+3)*(LAM+4)*(LAM+5)/6
C
C         CONTRIBUTION TO ESS0 OR ESSZ: UPPER-UPPER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T02(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE INDEX N IN ENSG ARRAY (CURRENTLY WITH [N=0,N'=1])
C         BY ONE ON THE A CENTRE, SO THAT IT OVERWRITES ESG VALUES AND
C         PROVIDES ES-COEFFICIENTS WITH [N=1,N'=1]. (LAM -> LAM + 2)
          CALL STEPN(ENSG,ESG,LAM+2,MAXM,1)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          NTUV  = (LAM+5)*(LAM+6)*(LAM+7)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: UPPER-UPPER FOR [N=1,N'=1]
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T00(M) = CAB*DFLOAT((2*LQN(1)+1)*(2*LQN(2)+1))
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: LOWER-LOWER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T00(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE FOR THE A CENTRE
          CALL STEPN(ESG,ENSG,LAM,MAXM,1)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          NTUV  = (LAM+3)*(LAM+4)*(LAM+5)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          PREFAC =-CAB*2.0D0*DFLOAT(2*LQN(2)+1)
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T20(M) = PREFAC*EXL(IBAS,1)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: LOWER-LOWER FOR [N=1,N'=0]
          DO M=1,MAXM
            TK = T20(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE FOR THE B CENTRE
          CALL STEPN(ESG,ENSG,LAM,MAXM,2)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          NTUV  = (LAM+3)*(LAM+4)*(LAM+5)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          PREFAC =-CAB*2.0D0*DFLOAT(2*LQN(1)+1)
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T02(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: LOWER-LOWER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T02(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE INDEX N IN ENSG ARRAY (CURRENTLY WITH [N=0,N'=1])
C         BY ONE ON THE A CENTRE, SO THAT IT OVERWRITES ESG VALUES AND
C         PROVIDES ES-COEFFICIENTS WITH [N=1,N'=1]. (LAM -> LAM + 2)
          CALL STEPN(ENSG,ESG,LAM+2,MAXM,1)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          NTUV  = (LAM+5)*(LAM+6)*(LAM+7)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESS0 OR ESSZ: LOWER-LOWER FOR [N=1,N'=1]
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ESSX AND ESSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T00(M) = CAB*DFLOAT((2*LQN(1)+1)*(2*LQN(2)+1))
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: UPPER-LOWER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T00(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE FOR THE A CENTRE FOR [N=1,N'=0]
          CALL STEPN(ESG,ENSG,LAM,MAXM,1)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)+2
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          PREFAC =-CAB*2.0D0*DFLOAT(2*LQN(2)+1)
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T20(M) = PREFAC*EXL(IBAS,1)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: UPPER-LOWER FOR [N=1,N'=0]
          DO M=1,MAXM
            TK = T20(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         DFNOTE: CHECK LAMBDA VALUE
C         INCREASE THE INDEX N BY ONE FOR THE B CENTRE FOR [N=0,N'=1]
          CALL STEPN(ESG,ENSG,LAM-2,MAXM,2)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION (DFNOTE: CHECK)
          LAM  = LLAB(1)+LLAB(2)+2
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          PREFAC =-CAB*2.0D0*DFLOAT(2*LQN(1)+1)
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T02(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: UPPER-LOWER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T02(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE INDEX N IN ENSG ARRAY (CURRENTLY WITH [N=0,N'=1])
C         BY ONE ON THE A CENTRE, SO THAT IT OVERWRITES ESG VALUES AND
C         PROVIDES ES-COEFFICIENTS WITH [N=1,N'=1]. (LAM -> LAM + 2)
          CALL STEPN(ENSG,ESG,LAM-2,MAXM,1)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)+4
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: UPPER-LOWER FOR [N=1,N'=1]
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) + IPHS*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T00(M) = CAB*DFLOAT((2*LQN(1)+1)*(2*LQN(2)+1))
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: LOWER-UPPER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T00(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE FOR THE A CENTRE
          CALL STEPN(ESG,ENSG,LAM,MAXM,1)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)+2
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          PREFAC =-CAB*2.0D0*DFLOAT(2*LQN(2)+1)
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T20(M) = PREFAC*EXL(IBAS,1)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: LOWER-UPPER FOR [N=1,N'=0]
          DO M=1,MAXM
            TK = T20(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE FOR THE B CENTRE
          CALL STEPN(ESG,ENSG,LAM-2,MAXM,2)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)+2
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          PREFAC =-CAB*2.0D0*DFLOAT(2*LQN(1)+1)
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M      = M+1
              T02(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: LOWER-UPPER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T02(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE INDEX N IN ENSG ARRAY (CURRENTLY WITH [N=0,N'=1])
C         BY ONE ON THE A CENTRE, SO THAT IT OVERWRITES ESG VALUES AND
C         PROVIDES ES-COEFFICIENTS WITH [N=1,N'=1]. (LAM -> LAM + 2)
          CALL STEPN(ENSG,ESG,LAM-2,MAXM,1)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)+4
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T22(M) = CAB*4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ESSX OR ESSY: LOWER-UPPER FOR [N=1,N'=1]
          DO M=1,MAXM
            TK = T22(M)
            DO ITUV=1,NTUV
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
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
      CALL RNSS(RNORM,EXL,LQN,NBAS)
C
C     MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
      LAM  = LQN(1)+LQN(2)+4
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     NORMALISE THE ESSQ COEFFICIENT BLOCK
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M   = M+1
          RSS = RNORM(IBAS,1)*RNORM(JBAS,2)
          DO ITUV=1,NTUV
            ESS(M,ITUV) = RSS*ESS(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     SIGMA_Y SPECIAL CASE: MULTIPLY ESSY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV
            ESS(M,ITUV) = CONE*ESS(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE EQLS(ELS,EXL,XYZ,KQN,MQN,NBAS,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                EEEEEEEE  QQQQQQ    LL       SSSSSS                   C
C                EE       QQ    QQ   LL      SS    SS                  C
C                EE      QQ      QQ  LL      SS                        C
C                EEEEEE  QQ      QQ  LL       SSSSSS                   C
C                EE      QQ      QQ  LL            SS                  C
C                EE       QQ    QQ   LL      SS    SS                  C
C                EEEEEEEE  QQQQQQ QQ LLLLLLLL SSSSSS                   C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQLS EVALUATES THE EQ-COEFFICIENTS FOR LARGE-SMALL CHARGE OVERLAP   C
C  OF G-SPINOR FUNCTIONS FOR ALL PAULI MATRICES IQ = {0,1,2,3}.        C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      COMPLEX*16 CONE
      COMPLEX*16 ELS(MB2,MEQ),ESG(MB2,MEQ),ENSG(MB2,MEQ)
C
      DIMENSION NBAS(2),KQN(2),JQN(2),LQN(2),MQN(2),LLAB(2),MLAB(2)
      DIMENSION EXL(MBS,2),RNORM(MBS,2),XYZ(3,2)
      DIMENSION T2(MB2),T0(MB2)
C
      COMMON/CTSN/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     THRESHOLD FOR CLEBSCH-GORDAN CEOFFICIENT MAGNITUDES
      SENS = 1.0D-10
C
C     SIGN MULTIPLIER FOR SIGMA COUPLING MATRICES
      IF(IQ.EQ.0.OR.IQ.EQ.1) THEN
        IPHS = 1
      ELSEIF(IQ.EQ.2.OR.IQ.EQ.3) THEN
        IPHS =-1
      ENDIF
C
C     INITIALISE COEFFICIENT STORAGE ARRAY TO ZERO
      DO M=1,MB2
        DO ITUV=1,MEQ
          ELS(M,ITUV) = DCMPLX(0.0D0,0.0D0)
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
        CBL = DSQRT(DFLOAT(JQN(2)+MQN(2)  )/DFLOAT(2*JQN(2)+4))
      ELSE
        CBU = DSQRT(DFLOAT(JQN(2)+MQN(2)  )/DFLOAT(2*JQN(2)  ))
        CBL = DSQRT(DFLOAT(JQN(2)-MQN(2)  )/DFLOAT(2*JQN(2)  ))
      ENDIF
C
C     DETERMINE THE NUMBER OF FUNCTIONS ON EACH CENTRE
      MAXM = NBAS(1)*NBAS(2)
C
C     KINETIC PRE-FACTORS FOR THIS BLOCK
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
          T2(M) =-2.0D0*EXL(JBAS,2)
          T0(M) = DFLOAT(2*LQN(2)+1)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INITIALISE COMMON GEOMETRICAL INFORMATION                        C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTRES A AND B
      DX  = (XYZ(1,1)-XYZ(1,2))**2
      DY  = (XYZ(2,1)-XYZ(2,2))**2
      DZ  = (XYZ(3,1)-XYZ(3,2))**2
      AB2 = DX+DY+DZ
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
C
          P(M)    = EXL(IBAS,1)+EXL(JBAS,2)
          P2(M)   = 2.0D0*P(M)
          P22(M)  = P2(M)*P2(M)
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
C     CASE 1: KQN(2).LT.0                                              C
C**********************************************************************C
C
      IF(KQN(2).GT.0) GOTO 100
C
C     BASIS PAIR LQN LABELS
      LLAB(1) = LQN(1)
      LLAB(2) = LQN(2)+1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ELS0 AND ELSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T2(M) =-CAB*2.0D0*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ELS0 OR ELSZ: UPPER-UPPER
          DO M=1,MAXM
            TK = T2(M)
            DO ITUV=1,NTUV
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T2(M) =-CAB*2.0D0*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ELS0 OR ELSZ: LOWER-LOWER
          DO M=1,MAXM
            TK = T2(M)
            DO ITUV=1,NTUV
              ELS(M,ITUV) = ELS(M,ITUV) + IPHS*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ELSX AND ELSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T2(M) =-CAB*2.0D0*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ELSX OR ELSY: UPPER-LOWER
          DO M=1,MAXM
            TK = T2(M)
            DO ITUV=1,NTUV
              ELS(M,ITUV) = ELS(M,ITUV) + IPHS*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T2(M) =-CAB*2.0D0*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ELSX OR ELSY: LOWER-UPPER
          DO M=1,MAXM
            TK = T2(M)
            DO ITUV=1,NTUV
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
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
C     BASIS PAIR LQN LABELS
      LLAB(1) = LQN(1)
      LLAB(2) = LQN(2)-1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ELS0 AND ELSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM   = LLAB(1)+LLAB(2)+1
          NTUV0 = (LAM+0)*(LAM+1)*(LAM+2)/6
          NTUV1 = (LAM+2)*(LAM+3)*(LAM+4)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T0(M) = CAB*DFLOAT(2*LQN(2)+1)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ELS0 OR ELSZ: UPPER-UPPER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T0(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T2(M) =-CAB*2.0D0*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ELS0 OR ELSZ: UPPER-UPPER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T2(M)
            DO ITUV=1,NTUV1
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM   = LLAB(1)+LLAB(2)+1
          NTUV0 = (LAM+0)*(LAM+1)*(LAM+2)/6
          NTUV1 = (LAM+2)*(LAM+3)*(LAM+4)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T0(M) = CAB*DFLOAT(2*LQN(2)+1)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ELS0 OR ELSZ: LOWER-LOWER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T0(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) + IPHS*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T2(M) =-CAB*2.0D0*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ELS0 OR ELSZ: LOWER-LOWER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T2(M)
            DO ITUV=1,NTUV1
              ELS(M,ITUV) = ELS(M,ITUV) + IPHS*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ELSX AND ELSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM   = LLAB(1)+LLAB(2)+1
          NTUV0 = (LAM+0)*(LAM+1)*(LAM+2)/6
          NTUV1 = (LAM+2)*(LAM+3)*(LAM+4)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T0(M) = CAB*DFLOAT(2*LQN(2)+1)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ELSX OR ELSY: UPPER-LOWER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T0(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) + IPHS*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T2(M) =-CAB*2.0D0*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ELSX OR ELSY: UPPER-LOWER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T2(M)
            DO ITUV=1,NTUV1
              ELS(M,ITUV) = ELS(M,ITUV) + IPHS*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM   = LLAB(1)+LLAB(2)+1
          NTUV0 = (LAM+0)*(LAM+1)*(LAM+2)/6
          NTUV1 = (LAM+2)*(LAM+3)*(LAM+4)/6
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T0(M) = CAB*DFLOAT(2*LQN(2)+1)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ELSX OR ELSY: LOWER-UPPER FOR [N=0,N'=0]
          DO M=1,MAXM
            TK = T0(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C         KINETIC AND ANGULAR PRE-FACTORS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M = M+1
              T2(M) =-CAB*2.0D0*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         CONTRIBUTION TO ELSX OR ELSY: LOWER-UPPER FOR [N=0,N'=1]
          DO M=1,MAXM
            TK = T2(M)
            DO ITUV=1,NTUV1
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF

200   CONTINUE
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
C     GENERATE RNLS NORMALISATION CONSTANTS
      CALL RNLS(RNORM,EXL,LQN,NBAS)
C
C     MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
      LAM  = LQN(1)+LQN(2)+1
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     NORMALISE THE ELSQ COEFFICIENT BLOCK
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M   = M+1
          RLS = RNORM(IBAS,1)*RNORM(JBAS,2)
          DO ITUV=1,NTUV
            ELS(M,ITUV) = RLS*ELS(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     SIGMA_Y SPECIAL CASE: MULTIPLY ELSY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV
            ELS(M,ITUV) = CONE*ELS(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE EQLS3(ELS,EXL,XYZ,KQN,MQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            EEEEEEEE  QQQQQQ    LL       SSSSSS   333333              C
C            EE       QQ    QQ   LL      SS    SS 33    33             C
C            EE      QQ      QQ  LL      SS             33             C
C            EEEEEE  QQ      QQ  LL       SSSSSS    33333              C
C            EE      QQ      QQ  LL            SS       33             C
C            EE       QQ    QQ   LL      SS    SS 33    33             C
C            EEEEEEEE  QQQQQQ QQ LLLLLLLL SSSSSS   333333              C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQLS EVALUATES THE EQ-COEFFICIENTS FOR LARGE-SMALL CHARGE OVERLAP   C
C  OF G-SPINOR FUNCTIONS FOR ALL PAULI MATRICES IQ = {0,1,2,3}.        C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      COMPLEX*16 CONE
      COMPLEX*16 ELS(MB2,MEQ,3),ESG(MB2,MEQ),ENSG(MB2,MEQ)
C
      DIMENSION NBAS(2),KQN(2),JQN(2),LQN(2),MQN(2),LLAB(2),MLAB(2)
      DIMENSION EXL(MBS,2),RNORM(MBS,2),XYZ(3,2)
      DIMENSION T2(MB2),T0(MB2)
C
      COMMON/CTSN/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     THRESHOLD FOR CLEBSCH-GORDAN CEOFFICIENT MAGNITUDES
      SENS = 1.0D-10
C
C     INITIALISE THE COEFFICIENTS TO ZERO
      DO M=1,MB2
        DO ITUV=1,MEQ
          DO IQ=1,3
            ELS(M,ITUV,IQ) = DCMPLX(0.0D0,0.0D0)
          ENDDO
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
C     CALCULATE THE APPROPRIATE CAB GORDAN FACTORS
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
C     DETERMINE THE NUMBER OF PAIRS IN THIS BLOCK
      MAXM  = NBAS(1)*NBAS(2)
C
C     KINETIC PRE-FACTORS FOR THIS BLOCK
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
          T2(M) =-2.0D0*EXL(JBAS,2)
          T0(M) = DFLOAT(2*LQN(2)+1)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INITIALISE COMMON GEOMETRICAL INFORMATION                        C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTRES A AND B
      DX  = (XYZ(1,1)-XYZ(1,2))**2
      DY  = (XYZ(2,1)-XYZ(2,2))**2
      DZ  = (XYZ(3,1)-XYZ(3,2))**2
      AB2 = DX+DY+DZ
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
C
          P(M)    = EXL(IBAS,1)+EXL(JBAS,2)
          P2(M)   = 2.0D0*P(M)
          P22(M)  = P2(M)*P2(M)
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
C     CASE 1: KQN(2).LT.0                                              C
C**********************************************************************C
C
      IF(KQN(2).GT.0) GOTO 100
C
C     BASIS PAIR LQN LABELS
      LLAB(1) = LQN(1)
      LLAB(2) = LQN(2)+1
C
C >>  TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MLAB(1) = (MQN(1)-1)/2
      MLAB(2) = (MQN(2)-1)/2
C
C     CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
      CAB = CAU*CBU
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM  = LLAB(1)+LLAB(2)
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C       KINETIC AND ANGULAR PRE-FACTORS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M = M+1
            T2(M) =-CAB*2.0D0*EXL(JBAS,2)
          ENDDO
        ENDDO
C
C       CONTRIBUTION TO ELSZ: UPPER-UPPER
        DO M=1,MAXM
          TK = T2(M)
          DO ITUV=1,NTUV
            ELS(M,ITUV,3) = ELS(M,ITUV,3) + TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
C >>  TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MLAB(1) = (MQN(1)+1)/2
      MLAB(2) = (MQN(2)+1)/2
C
C     CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
      CAB = CAL*CBL
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM  = LLAB(1)+LLAB(2)
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C       KINETIC AND ANGULAR PRE-FACTORS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M = M+1
            T2(M) = -CAB*2.0D0*EXL(JBAS,2)
          ENDDO
        ENDDO
C
C       CONTRIBUTION TO ELSZ: LOWER-LOWER
        DO M=1,MAXM
          TK = T2(M)
          DO ITUV=1,NTUV
            ELS(M,ITUV,3) = ELS(M,ITUV,3) - TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ELSX AND ELSY
C
C >>  TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MLAB(1) = (MQN(1)-1)/2
      MLAB(2) = (MQN(2)+1)/2
C
C     CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
      CAB = CAU*CBL
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM  = LLAB(1)+LLAB(2)
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C       KINETIC AND ANGULAR PRE-FACTORS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M = M+1
            T2(M) =-CAB*2.0D0*EXL(JBAS,2)
          ENDDO
        ENDDO
C
C       CONTRIBUTION TO ELSX AND ELSY: UPPER-LOWER
        DO M=1,MAXM
          TK = T2(M)
          DO ITUV=1,NTUV
            ELS(M,ITUV,1) = ELS(M,ITUV,1) + TK*ESG(M,ITUV)
            ELS(M,ITUV,2) = ELS(M,ITUV,2) - TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
C >>  TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MLAB(1) = (MQN(1)+1)/2
      MLAB(2) = (MQN(2)-1)/2
C
C     CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
      CAB = CAL*CBU
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM  = LLAB(1)+LLAB(2)
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C       KINETIC AND ANGULAR PRE-FACTORS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M = M+1
            T2(M) =-CAB*2.0D0*EXL(JBAS,2)
          ENDDO
        ENDDO
C
C       CONTRIBUTION TO ELSX AND ELSY: LOWER-UPPER
        DO M=1,MAXM
          TK = T2(M)
          DO ITUV=1,NTUV
            ELS(M,ITUV,1) = ELS(M,ITUV,1) + TK*ESG(M,ITUV)
            ELS(M,ITUV,2) = ELS(M,ITUV,2) + TK*ESG(M,ITUV)
          ENDDO
        ENDDO
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
C     BASIS PAIR LQN LABELS
      LLAB(1) = LQN(1)
      LLAB(2) = LQN(2)-1
C
C     INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
      LAM   = LLAB(1)+LLAB(2)+1
      NTUV0 = (LAM+0)*(LAM+1)*(LAM+2)/6
      NTUV1 = (LAM+2)*(LAM+3)*(LAM+4)/6
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ELS0 AND ELSZ
C
C >>  TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MLAB(1) = (MQN(1)-1)/2
      MLAB(2) = (MQN(2)-1)/2
C
C     CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
      CAB = CAU*CBU
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM  = LLAB(1)+LLAB(2)
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C       KINETIC AND ANGULAR PRE-FACTORS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M = M+1
            T0(M) = CAB*DFLOAT(2*LQN(2)+1)
          ENDDO
        ENDDO
C
C       CONTRIBUTION TO ELSZ: UPPER-UPPER FOR [N=0,N'=0]
        DO M=1,MAXM
          TK = T0(M)
          DO ITUV=1,NTUV0
            ELS(M,ITUV,3) = ELS(M,ITUV,3) + TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
C       INCREASE THE INDEX N BY ONE
        CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C       KINETIC AND ANGULAR PRE-FACTORS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M = M+1
            T2(M) =-CAB*2.0D0*EXL(JBAS,2)
          ENDDO
        ENDDO
C
C       CONTRIBUTION TO ELSZ: UPPER-UPPER FOR [N=0,N'=1]
        DO M=1,MAXM
          TK = T2(M)
          DO ITUV=1,NTUV1
            ELS(M,ITUV,3) = ELS(M,ITUV,3) + TK*ENSG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
C >>  TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MLAB(1) = (MQN(1)+1)/2
      MLAB(2) = (MQN(2)+1)/2
C
C     CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
      CAB = CAL*CBL
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM  = LLAB(1)+LLAB(2)
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C       KINETIC AND ANGULAR PRE-FACTORS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M = M+1
            T0(M) = CAB*DFLOAT(2*LQN(2)+1)
          ENDDO
        ENDDO
C
C       CONTRIBUTION TO ELSZ: LOWER-LOWER FOR [N=0,N'=0]
        DO M=1,MAXM
          TK = T0(M)
          DO ITUV=1,NTUV0
            ELS(M,ITUV,3) = ELS(M,ITUV,3) - TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
C       INCREASE THE INDEX N BY ONE
        CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C       KINETIC AND ANGULAR PRE-FACTORS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M = M+1
            T2(M) =-CAB*2.0D0*EXL(JBAS,2)
          ENDDO
        ENDDO
C
C       CONTRIBUTION TO ELSZ: LOWER-LOWER FOR [N=0,N'=1]
        DO M=1,MAXM
          TK = T2(M)
          DO ITUV=1,NTUV1
            ELS(M,ITUV,3) = ELS(M,ITUV,3) - TK*ENSG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ELSX AND ELSY
C
C >>  TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MLAB(1) = (MQN(1)-1)/2
      MLAB(2) = (MQN(2)+1)/2
C
C     CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
      CAB = CAU*CBL
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM  = LLAB(1)+LLAB(2)
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C       KINETIC AND ANGULAR PRE-FACTORS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M = M+1
            T0(M) = CAB*DFLOAT(2*LQN(2)+1)
          ENDDO
        ENDDO
C
C       CONTRIBUTION TO ELSX AND ELSY: UPPER-LOWER FOR [N=0,N'=0]
        DO M=1,MAXM
          TK = T0(M)
          DO ITUV=1,NTUV0
            ELS(M,ITUV,1) = ELS(M,ITUV,1) + TK*ESG(M,ITUV)
            ELS(M,ITUV,2) = ELS(M,ITUV,2) - TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
C       INCREASE THE INDEX N BY ONE
        CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C       KINETIC AND ANGULAR PRE-FACTORS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M = M+1
            T2(M) =-CAB*2.0D0*EXL(JBAS,2)
          ENDDO
        ENDDO
C
C       CONTRIBUTION TO ELSX AND ELSY: UPPER-LOWER FOR [N=0,N'=1]
        DO M=1,MAXM
          TK = T2(M)
          DO ITUV=1,NTUV1
            ELS(M,ITUV,1) = ELS(M,ITUV,1) + TK*ENSG(M,ITUV)
            ELS(M,ITUV,2) = ELS(M,ITUV,2) - TK*ENSG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
C >>  TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C     BASIS PAIR MQN LABELS
      MLAB(1) = (MQN(1)+1)/2
      MLAB(2) = (MQN(2)-1)/2
C
C     CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
      CAB = CAL*CBU
C
C     SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
      IF(DABS(CAB).GE.SENS) THEN
C
C       GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
        CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C       MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
        LAM  = LLAB(1)+LLAB(2)
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C       KINETIC AND ANGULAR PRE-FACTORS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M = M+1
            T0(M) = CAB*DFLOAT(2*LQN(2)+1)
          ENDDO
        ENDDO
C
C       CONTRIBUTION TO ELSX AND ELSY: LOWER-UPPER FOR [N=0,N'=0]
        DO M=1,MAXM
          TK = T0(M)
          DO ITUV=1,NTUV0
            ELS(M,ITUV,1) = ELS(M,ITUV,1) + TK*ESG(M,ITUV)
            ELS(M,ITUV,2) = ELS(M,ITUV,2) + TK*ESG(M,ITUV)
          ENDDO
        ENDDO
C
C       INCREASE THE INDEX N BY ONE
        CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C       KINETIC AND ANGULAR PRE-FACTORS
        M = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M = M+1
            T2(M) =-CAB*2.0D0*EXL(JBAS,2)
          ENDDO
        ENDDO
C
C       CONTRIBUTION TO ELSX AND ELSY: LOWER-UPPER FOR [N=0,N'=1]
        DO M=1,MAXM
          TK = T2(M)
          DO ITUV=1,NTUV1
            ELS(M,ITUV,1) = ELS(M,ITUV,1) + TK*ENSG(M,ITUV)
            ELS(M,ITUV,2) = ELS(M,ITUV,2) + TK*ENSG(M,ITUV)
          ENDDO
        ENDDO
C
      ENDIF
C
200   CONTINUE
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
C     GENERATE RNLS NORMALISATION CONSTANTS
      CALL RNLS(RNORM,EXL,LQN,NBAS)
C
C     MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
      LAM  = LQN(1)+LQN(2)+1
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     NORMALISE THE ELLQ COEFFICIENT BLOCK
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M   = M+1
          RLS = RNORM(IBAS,1)*RNORM(JBAS,2)
          DO ITUV=1,NTUV
            DO IQ=1,3
              ELS(M,ITUV,IQ) = RLS*ELS(M,ITUV,IQ)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     SIGMA_Y SPECIAL CASE: MULTIPLY ELSY RESULTS BY i.
      DO M=1,MAXM
        DO ITUV=1,NTUV
          ELS(M,ITUV,2) = CONE*ELS(M,ITUV,2)
        ENDDO
      ENDDO
CC
CC     MULTIPLY ALL RESULTS BY A FACTOR i
C      DO M=1,MAXM
C        DO ITUV=1,NTUV
C          DO IQ=1,3
C            ELS(M,ITUV,IQ) = CONE*ELS(M,ITUV,IQ)
C          ENDDO
C        ENDDO
C      ENDDO
CC
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
C  THE REQUIRED COEFFICIENTS ARE GENERATED BY A CALL TO VRS, WHICH IS  C
C  CONSTRUCTED ACCORDING TO THE RECURRENCE RELATIONS DEFINED BY        C
C  V.R. SAUNDERS. THE OUTPUT OF VRS IS THEN ADJUSTED TO INCLUDE THE    C
C  ANGULAR NORMALISATION CONSTANTS, AS WELL AS A PHASE FACTOR TO       C
C  CONVERT FROM THE SCHIFF TO THE CONDON-SHORTLEY PHASE CONVENTION.    C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C     LQN(I) - TARGET LQN VALUES ON CENTRES A AND B.                   C
C     MQN(I) - TARGET MQN VALUES ON CENTRES A AND B.                   C
C -------------------------------------------------------------------- C
C  H.M.QUINEY THE UNIVERSITY OF MELBOURNE (2008).                      C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      COMPLEX*16 ESG(MB2,MEQ)
C
      DIMENSION FACT(MKP),LQN(2),MQN(2),MQNLAB(2)
C
      COMMON/CTSN/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
      PI4 = 8.0D0*DASIN(1.0D0)
C
C     CALCULATE THE FACTORIAL FUNCTIONS
      LMAX = MAX(LQN(1),LQN(2))
      FACT(1) = 1.0D0
      DO M=1,2*LMAX
        FACT(M+1) = FACT(M)*DFLOAT(M)
      ENDDO
C
C     VRS IS CALLED WITH THE SIGN OF MQN(1) REVERSED
C     TO AFFECT COMPLEX CONJUGATION, ALONG WITH THE REQUISITE
C     PHASE, WHICH IS CALCULATED LATER.
      MQNLAB(1) =-MQN(1)
      MQNLAB(2) = MQN(2)
C
C     TRAP CASES FOR WHICH |MQN| EXCEEDS LQN (COULD BE CALLED BUT
C     WITH A ZERO MULTIPLICATIVE CONSTANT)
      IF(IABS(MQN(1)).GT.LQN(1).OR.IABS(MQN(2)).GT.LQN(2)) THEN
        LAMAB = LQN(1)+LQN(2)
        NTUV  = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
        DO ITUV=1,NTUV
          DO M=1,MAXM
            ESG(M,ITUV) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
        RETURN
      ELSE
        CALL VRS(ESG,LQN,MQNLAB,MAXM)
      ENDIF
C
C     IMPORT L AND BASIS PAIR MQNS
      L1 = LQN(1)
      L2 = LQN(2)
      M1 = IABS(MQN(1))
      M2 = IABS(MQN(2))
C
C     SPECIFY THE UPPER TERMINAL ON SUMMATION, AND CG COEFFICIENTS
      LAM    = LQN(1) + LQN(2)
      PREFAC = DFLOAT((2*L1+1)*(2*L2+1))
      PREFAC = PREFAC*FACT(L1-M1+1)/FACT(L1+M1+1)
      PREFAC = PREFAC*FACT(L2-M2+1)/FACT(L2+M2+1)
      PREFAC = DSQRT(PREFAC)
C     PHASE  = (-1.0D0)**(M1+M2+MQN(2))
      PHASE  = (-1.0D0)**((MQN(1)+MQN(2)+M1+M2)/2)
      PREFAC = PREFAC*PHASE/PI4
C
C     THERE ARE NTUV TOTAL TERMS IN THE SUM OVER A,B,C
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
      DO ITUV=1,NTUV
        DO M=1,MAXM
          ESG(M,ITUV) = PREFAC*ESG(M,ITUV)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE VRS(ESG,LQN,MQN,MAXM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                      VV    VV RRRRRRR   SSSSSS                       C
C                      VV    VV RR    RR SS    SS                      C
C                      VV    VV RR    RR SS                            C
C                      VV    VV RR    RR  SSSSSS                       C
C                       VV  VV  RRRRRRR        SS                      C
C                        VVVV   RR    RR SS    SS                      C
C                         VV    RR    RR  SSSSSS                       C
C                                                                      C
C -------------------------------------------------------------------- C
C  VRS EVALUATES THE EXPANSION COEFFICIENTS OF THE OVERLAP CHARGE      C
C  DENSITY OF SGTFS IN AN AUXILIARY HGTF. COEFFICIENTS ARE EVALUATED   C
C  USING THE RECURRENCE RELATIONS DEFINED BY VIC SAUNDERS IN:          C
C                                                                      C
C  V.R. SAUNDERS, "MOLECULAR INTEGRALS FOR GAUSSIAN-TYPE FUNCTIONS",   C
C  METHODS OF COMPUTATIONAL MOLECULAR PHYSICS, ED G.H.F.DIERCKSEN AND  C
C  S.WILSON, pp 1-26, REIDEL PUBLISHING, DORDRECHT (1983).             C
C                                                                      C
C  THE E-COEFFS IN THIS PROCEDURE ARE FOR AN UN-NORMALISED SGTF.       C
C  THE COEFFICIENTS ARE DETERMINED ACCORDING TO THE SAME RULES AS      C
C  DEFINED IN THE ABOVE ARTICLE. CONSEQUENTLY, IT SHOULD BE NOTED THAT C
C  THE COEFFICIENTS ARE THOSE OF SPHERICAL HARMONIC FUNCTIONS THAT ARE C
C    (*) UN-NORMALISED                                                 C
C    (*) SATISFY THE SCHIFF PHASE CONVENTION                           C
C                                                                      C
C  THE OUTLINE FOR THE GENERATION OF E-COEFFICIENTS IS TAKEN FROM p16  C
C  OF THE ABOVE ARTICLE. EQUATION NUMBERS ARE GIVEN IN COMMENTS.       C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C     LQN(I) - TARGET LQN VALUES ON CENTRES A AND B.                   C
C     MQN(I) - TARGET MQN VALUES ON CENTRES A AND B.                   C
C -------------------------------------------------------------------- C
C  H.M.QUINEY THE UNIVERSITY OF MELBOURNE (2008).                      C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=ML2*2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      COMPLEX*16 ESG(MB2,MEQ),ETEMP(MB2*MRC)
C
      DIMENSION NBAS(2),LQN(2),MQN(2)
C
      COMMON/CTSN/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     IMPORT LQN AND BASIS PAIR MQNS FOR LOCAL USE
      LQNA = LQN(1)
      LQNB = LQN(2)
      MQNA = MQN(1)
      MQNB = MQN(2)
      LMAX = LQNA + LQNB
C
C     CHECK THAT LMAX IS WITHIN THE BOUNDS OF MKP
      IF(LMAX.GT.MKP-1) THEN
        WRITE(6,20) LMAX,MKP-1
        WRITE(7,20) LMAX,MKP-1
        STOP
      ENDIF
20    FORMAT(2X,'Required value of LAMBDA = ',I3/
     &       2X,'Maximum allowed value of LAMBDA = ',I3//
     &       2X,'Reset MKP and recompile: terminating...'/)
C
C      SET INITIAL VALUES TO E[0,0;0,0;0,0,0,0] = RKAB
       DO M=1,MB2*MRC
         ETEMP(M) = DCMPLX(0.0D0,0.0D0)
       ENDDO
C
       DO M=1,MAXM
         ETEMP(M) = DCMPLX(RKAB(M),0.0D0)
       ENDDO
C
C      STEP 1:
C      GENERATE E[|MQNA|,MQNA;0,0] FROM E[0,0;0,0] USING
C      SIMULTANEOUS STEP OF LQN AND MQN ON CENTRE A
       ISTART = 0
       LAM    = 0
       IF(IABS(MQNA).NE.0) THEN
         CALL STEPLM(ETEMP,LAM,ISTART,MQNA,MAXM,1)
       ENDIF
C
C      STEP 2:
C      GENERATE E[LQNA,MQNA;0,0] FROM E[|MQNA|,MQNA;0,0]
C      USING THE STEP OF LQN ONLY ON CENTRE A
       IF(LQNA.GT.IABS(MQNA)) THEN
         CALL STEPL(ETEMP,LAM,ISTART,LQNA,MQNA,MAXM,1)
       ENDIF
C
C      STEP 3:
C      GENERATE E[LQNA,MQNA;|MQNB|,MQNB] FROM E[LQNA,MQNA;0,0]
C      USING SIMULTANEOUS STEP OF LQN AND MQN ON CENTRE B
       IF(IABS(MQNB).GT.0) THEN
         CALL STEPLM(ETEMP,LAM,ISTART,MQNB,MAXM,2)
       ENDIF
C
C      STEP 4:
C      GENERATE E[LQNA,MQNA;LQNB,MQNB] FROM E[LQNA,MQNA;|MQNB|,MQNB]
C      USING THE STEP OF LQN ONLY ON CENTRE B
       IF(LQNB.GT.IABS(MQNB)) THEN
         CALL STEPL(ETEMP,LAM,ISTART,LQNB,MQNB,MAXM,2)
       ENDIF
C
C      STEP 5:
C      COPY FINAL BLOCK OF ENTRIES AS THE REQUIRED OUTPUT
       ISTART0 = ISTART
       NTUV    = (LAM+1)*(LAM+2)*(LAM+3)/6
C
       K = 0
       DO ITUV=1,NTUV
         DO M=1,MAXM
           K = K+1
           ESG(M,ITUV) = ETEMP(ISTART0+K)
         ENDDO
       ENDDO
C
       RETURN
       END
C
C
      SUBROUTINE STEPLM(ETEMP,LAM,ISTART,MQN,MAXM,ICNT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        SSSSSS TTTTTTTT EEEEEEEE PPPPPPP  LL       MM       MM        C
C       SS    SS   TT    EE       PP    PP LL       MMM     MMM        C
C       SS         TT    EE       PP    PP LL       MMMM   MMMM        C
C        SSSSSS    TT    EEEEEE   PP    PP LL       MM MM MM MM        C
C             SS   TT    EE       PPPPPPP  LL       MM  MMM  MM        C
C       SS    SS   TT    EE       PP       LL       MM   M   MM        C
C        SSSSSS    TT    EEEEEEEE PP       LLLLLLLL MM       MM        C
C                                                                      C
C -------------------------------------------------------------------- C
C   SIMULTANEOUSLY INCREMENT/DECREMENT QUANTUM NUMBERS (LQN,MQN):      C
C               E[L, L;IT,IU,IV] -> E[L+1,L+1;IT,IU,IV]                C
C               E[L,-L;IT,IU,IV] -> E[L+1,L-1;IT,IU,IV]                C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    P2(M) - CONTAINS VALUES OF P*2, P=SUM OF EXPONENTS.               C
C    PX(M) - GEOMETRICAL VALUES OF X(M)-COORD(ICNT,X).                 C
C    PY(M) - GEOMETRICAL VALUES OF Y(M)-COORD(ICNT,Y).                 C
C    MAXM  - NUMBER OF EXPONENT/DENSITY PAIRS.                         C
C    LAM   - LENGTH OF THE INPUT HGTF EXPANSION.                       C
C    LQN   - L-QUANTUM NUMBER OF THE CENTRE TO BE INCREMENTED.         C
C    ICNT  - CENTRE TO STEP UP.                                        C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                      ML4=2*(MKP+1),MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      COMPLEX*16 CONE
      COMPLEX*16 ETEMP(*)
C
      DIMENSION PX(MB2),PY(MB2)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/CTSN/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     IMPORT GEOMETRICAL VALUES FOR CENTRE OF INTEREST
      DO M=1,MAXM
        IF(ICNT.EQ.1) THEN
          PX(M) = PAX(M)
          PY(M) = PAY(M)
        ELSEIF(ICNT.EQ.2) THEN
          PX(M) = PBX(M)
          PY(M) = PBY(M)
        ENDIF
      ENDDO
C
C     MAIN LOOP: FOR EACH M-QUANTUM NUMBER ON THIS CENTRE
      DO 100 MVAL=0,IABS(MQN)-1
C
C**********************************************************************C
C     COMPUTE THE BLOCK INDICES. THE RECURRENCE WILL RUN OVER          C
C     (LAM+1)*(LAM+2)*(LAM+3)/6 VALUES AND WILL GENERATE               C
C     (LAM+2)*(LAM+3)*(LAM+4)/6 VALUES IN THE NEXT LAYER               C
C**********************************************************************C
C
      RL1     = DFLOAT(2*IABS(MVAL)+1)
      NTUV    = (LAM+1)*(LAM+2)*(LAM+3)/6
      ISTART1 = ISTART
      ISTART2 = ISTART1+NTUV*MAXM
C
C**********************************************************************C
C                          INDEX MAPPINGS:                             C
C -------------------------------------------------------------------- C
C          I0-> E[LQN  ,LQN  ;IT  ,IU  ,IV]                            C
C          I1-> E[LQN+1,LQN+1;IT  ,IU  ,IV]                            C
C          I2-> E[LQN+1,LQN+1;IT+1,IU  ,IV]                            C
C          I3-> E[LQN+1,LQN+1;IT  ,IU+1,IV]                            C
C          I4-> E[LQN+1,LQN+1;IT-1,IU  ,IV]                            C
C          I5-> E[LQN+1,LQN+1;IT  ,IU-1,IV]                            C
C**********************************************************************C
C
C     INCREMENT THE M-QUANTUM NUMBER IF MQN > 0
      IF(MQN.GT.0) THEN
C       LOOP OVER THE HGTF INDICES OF THE SEED LAYER
        DO IOUTER=0,LAM
          DO IT=0,IOUTER
            DO IU=0,IOUTER-IT
              IV = IOUTER-IT-IU
              I0 = ISTART1 + (INABCD(IT  ,IU  ,IV)-1)*MAXM
              I1 = ISTART2 + (INABCD(IT  ,IU  ,IV)-1)*MAXM
              I2 = ISTART2 + (INABCD(IT+1,IU  ,IV)-1)*MAXM
              I3 = ISTART2 + (INABCD(IT  ,IU+1,IV)-1)*MAXM
C
              DO M=1,MAXM
                 T1 = RL1/P2(M)
                 TX = RL1*PX(M)
                 TY = RL1*PY(M)
                 ETEMP(I1+M) = ETEMP(I1+M) + TX*ETEMP(I0+M)
     &                                     + TY*ETEMP(I0+M)*CONE
                 ETEMP(I2+M) = ETEMP(I2+M) + T1*ETEMP(I0+M)
                 ETEMP(I3+M) = ETEMP(I3+M) + T1*ETEMP(I0+M)*CONE
              ENDDO
C
C             SPECIAL CASE EXCLUDES IT=0
              IF(IT.NE.0) THEN
                I4 = ISTART2 + (INABCD(IT-1,IU,IV)-1)*MAXM
                RT = DFLOAT(IT)
                FACTOR = RT*RL1
                DO M=1,MAXM
                  ETEMP(I4+M) = ETEMP(I4+M) + FACTOR*ETEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IU=0
              IF(IU.NE.0) THEN
                I5 = ISTART2 + (INABCD(IT,IU-1,IV)-1)*MAXM
                RU = DFLOAT(IU)
                FACTOR = RL1*RU
                DO M=1,MAXM
                  ETEMP(I5+M) = ETEMP(I5+M) + FACTOR*ETEMP(I0+M)*CONE
                ENDDO
              ENDIF
C
C           END OF LOOPS OVER HGTF INDICES
            ENDDO
          ENDDO
        ENDDO
C
C     DECREMENT THE M-QUANTUM NUMBER IF MQN < 0
      ELSE
C
C**********************************************************************C
C                          INDEX MAPPINGS:                             C
C -------------------------------------------------------------------- C
C          ISTART0 LABELS THE PREVIOUS LQN VALUE                       C
C          ISTART1 LABELS THE CURRENT  LQN VALUE                       C
C          ISTART2 LABELS THE NEXT     LQN VALUE                       C
C -------------------------------------------------------------------- C
C          I0-> E[LQN  ,LQN  ;IT  ,IU  ,IV]                            C
C          I1-> E[LQN+1,LQN+1;IT  ,IU  ,IV]                            C
C          I2-> E[LQN+1,LQN+1;IT+1,IU  ,IV]                            C
C          I3-> E[LQN+1,LQN+1;IT  ,IU+1,IV]                            C
C          I4-> E[LQN+1,LQN+1;IT-1,IU  ,IV]                            C
C          I5-> E[LQN+1,LQN+1;IT  ,IU-1,IV]                            C
C**********************************************************************C
C
C       LOOP OVER THE HGTF INDICES OF THE SEED LAYER
        DO IOUTER=0,LAM
          DO IT=0,IOUTER
            DO IU=0,IOUTER-IT
              IV = IOUTER-IT-IU
              I0 = ISTART1 + (INABCD(IT  ,IU  ,IV)-1)*MAXM
              I1 = ISTART2 + (INABCD(IT  ,IU  ,IV)-1)*MAXM
              I2 = ISTART2 + (INABCD(IT+1,IU  ,IV)-1)*MAXM
              I3 = ISTART2 + (INABCD(IT  ,IU+1,IV)-1)*MAXM
C
              DO M=1,MAXM
                T1 = RL1/P2(M)
                TX = RL1*PX(M)
                TY = RL1*PY(M)
                ETEMP(I1+M) = ETEMP(I1+M) + TX*ETEMP(I0+M)
     &                                    - TY*ETEMP(I0+M)*CONE
                ETEMP(I2+M) = ETEMP(I2+M) + T1*ETEMP(I0+M)
                ETEMP(I3+M) = ETEMP(I3+M) - T1*ETEMP(I0+M)*CONE
              ENDDO
C
C             SPECIAL CASE EXCLUDES IT=0
              IF(IT.NE.0) THEN
                I4 = ISTART2 + (INABCD(IT-1,IU  ,IV  )-1)*MAXM
                RT = DFLOAT(IT)
                FACTOR = RT*RL1
                DO M=1,MAXM
                  ETEMP(I4+M) = ETEMP(I4+M) + FACTOR*ETEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IU=0
              IF(IU.NE.0) THEN
                I5 = ISTART2 + (INABCD(IT  ,IU-1,IV  )-1)*MAXM
                RU = DFLOAT(IU)
                FACTOR = RL1*RU
                DO M=1,MAXM
                  ETEMP(I5+M) = ETEMP(I5+M) - FACTOR*ETEMP(I0+M)*CONE
                ENDDO
              ENDIF
C
C             END OF LOOPS OVER HGTF INDICES
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
C**********************************************************************C
C     END OF LOOP OVER MQN COUNTER. UPDATE THE VALUE OF LAM, AND       C
C     THE COUNTER THAT KEEPS TRACK OF THE BLOCKS OF E-COEFFICIENTS     C
C**********************************************************************C
C
      ISTART = ISTART + NTUV*MAXM
      LAM    = LAM + 1
C
100   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE STEPL(ETEMP,LAM,ISTART,LQN,MQN,MAXM,ICNT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C              SSSSSS TTTTTTTT EEEEEEEE PPPPPPP  LL                    C
C             SS    SS   TT    EE       PP    PP LL                    C
C             SS         TT    EE       PP    PP LL                    C
C              SSSSSS    TT    EEEEEE   PP    PP LL                    C
C                   SS   TT    EE       PPPPPPP  LL                    C
C             SS    SS   TT    EE       PP       LL                    C
C              SSSSSS    TT    EEEEEEEE PP       LLLLLLLL              C
C                                                                      C
C -------------------------------------------------------------------- C
C  INCREMENT THE LQN, STARTING AT E[L, L] OR E[L,-L].                  C
C -------------------------------------------------------------------- C
C  NOTE THAT THE FIRST APPLICATION OF EQ(64a) CAN ONLY MAP             C
C  E[L,L] -> E[L+1,M] OR E[L,-L] -> E[L+1,-L] AND IS TREATED SEP'TLY.  C
C  SUBSEQUENT STEPS MAP {E[L,L], E[L-1,L]} -> E[L+1,L].                C
C  FINAL APPLICATION OF THIS RULE GENERATES E[LMAX,L;0,0;T,U,V].       C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                      ML4=2*(MKP+1),MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      COMPLEX*16 CONE
      COMPLEX*16 ETEMP(*)
C
      DIMENSION PP(MB2),PX(MB2),PY(MB2),PZ(MB2)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/CTSN/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     IMPORT GEOMETRICAL VALUES FOR CENTRE OF INTEREST
      DO M=1,MAXM
        IF(ICNT.EQ.1) THEN
          PP(M) = PA2(M)
          PX(M) = PAX(M)
          PY(M) = PAY(M)
          PZ(M) = PAZ(M)
        ELSEIF(ICNT.EQ.2) THEN
          PP(M) = PB2(M)
          PX(M) = PBX(M)
          PY(M) = PBY(M)
          PZ(M) = PBZ(M)
        ENDIF
      ENDDO
C
C**********************************************************************C
C      THE FIRST STEP IS ALWAYS PERFORMED. IT MAPS THE INDEX SETS      C
C      E[MQN+1,MQN] <- E[MQN1,MQN1] FROM THE DATA OBTAINED IN STEPLM.  C
C**********************************************************************C
C
C     IF STEPL IS ENTERED WITH LQN.LE.|MQN| CONTROL IS RETURNED
C     AND NO COUNTERS ARE UPDATED
      IF(LQN.LE.IABS(MQN)) RETURN
C
      NTUV    = (LAM+1)*(LAM+2)*(LAM+3)/6
      ISTART1 = ISTART
      ISTART2 = ISTART1 + NTUV*MAXM
      RFACT1  = DFLOAT(2*IABS(MQN)+1)
C
C**********************************************************************C
C                          INDEX MAPPINGS:                             C
C -------------------------------------------------------------------- C
C       I0-> E[MQN  ,MQN;IT  ,IU  ,IV  ]                               C
C       I1-> E[MQN+1,MQN;IT  ,IU  ,IV  ]                               C
C       I2-> E[MQN+1,MQN;IT  ,IU  ,IV+1]                               C
C       I3-> E[MQN+1,MQN;IT  ,IU  ,IV-1]                               C
C**********************************************************************C
C
C     LOOP OVER THE HGTF INDICES OF THE SEED LAYER
      DO IOUTER=0,LAM
        DO IT=0,IOUTER
          DO IU=0,IOUTER-IT
            IV = IOUTER-IT-IU
            I0 = ISTART1 + (INABCD(IT  ,IU  ,IV  )-1)*MAXM
            I1 = ISTART2 + (INABCD(IT  ,IU  ,IV  )-1)*MAXM
            I2 = ISTART2 + (INABCD(IT  ,IU  ,IV+1)-1)*MAXM
            DO M=1,MAXM
              TZ = RFACT1*PZ(M)
              TP = RFACT1/P2(M)
              ETEMP(I1+M) = ETEMP(I1+M) + TZ*ETEMP(I0+M)
              ETEMP(I2+M) = ETEMP(I2+M) + TP*ETEMP(I0+M)
            ENDDO
            IF(IV.GE.1) THEN
              I3 = ISTART2 + (INABCD(IT  ,IU  ,IV-1)-1)*MAXM
              FACTOR = RFACT1*DFLOAT(IV)
              DO M=1,MAXM
                ETEMP(I3+M) = ETEMP(I3+M) + ETEMP(I0+M)*FACTOR
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
C     UPDATE LAM INDEX AND BLOCK LOCATOR
      ISTART0 = ISTART
      ISTART  = ISTART + NTUV*MAXM
      LAM     = LAM + 1
C
      IF(LQN.EQ.IABS(MQN)+1) RETURN
C
C**********************************************************************C
C     SECOND AND SUBSEQUENT STEPS IN THIS RECURRENCE INVOLVE THREE     C
C     LAYERS OF COEFFICIENTS:                                          C
C     E[LQN1+1,MQN1] <- {E[LQN1,MQN1],E[LQN-1,MQN1]}                   C
C**********************************************************************C
C
      DO LQN1=IABS(MQN)+1,LQN-1
        RL1M1   = DFLOAT(LQN1-IABS(MQN)+1)
        RFACT1  = DFLOAT(2*LQN1+1)/RL1M1
        NTUV    = (LAM+1)*(LAM+2)*(LAM+3)/6
        ISTART1 = ISTART
        ISTART2 = ISTART + MAXM*NTUV
C
C**********************************************************************C
C                          INDEX MAPPINGS:                             C
C -------------------------------------------------------------------- C
C     I0-> E[LQN  ,MQN;IT  ,IU  ,IV  ]                                 C
C     I1-> E[LQN+1,MQN;IT  ,IU  ,IV  ]                                 C
C     I2-> E[LQN+1,MQN;IT  ,IU  ,IV+1]                                 C
C     I3-> E[LQN+1,MQN;IT  ,IU  ,IV-1]                                 C
C**********************************************************************C
C
C       THE FIRST LOOP OVER ITUV INCLUDES ALL HGTF INDICES ON THE
C       LAYER CORRESPONDING TO THE CURRENT VALUE OF LQN1
        DO IOUTER=0,LAM
          DO IT=0,IOUTER
            DO IU=0,IOUTER-IT
              IV = IOUTER-IT-IU
              I0 = ISTART1 + (INABCD(IT  ,IU  ,IV  )-1)*MAXM
              I1 = ISTART2 + (INABCD(IT  ,IU  ,IV  )-1)*MAXM
              I2 = ISTART2 + (INABCD(IT  ,IU  ,IV+1)-1)*MAXM
              DO M=1,MAXM
                TZ = RFACT1*PZ(M)
                TP = RFACT1/P2(M)
                ETEMP(I1+M) = ETEMP(I1+M) + TZ*ETEMP(I0+M)
                ETEMP(I2+M) = ETEMP(I2+M) + TP*ETEMP(I0+M)
              ENDDO
              IF(IV.GE.1) THEN
                I3 = ISTART2 + (INABCD(IT  ,IU  ,IV-1)-1)*MAXM
                FACTOR = RFACT1*DFLOAT(IV)
                DO M=1,MAXM
                  ETEMP(I3+M) = ETEMP(I3+M) + ETEMP(I0+M)*FACTOR
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
C**********************************************************************C
C                          INDEX MAPPINGS:                             C
C -------------------------------------------------------------------- C
C     I0 -> E[LQN-1,MQN;IT  ,IU  ,IV  ]                                C
C     I1 -> E[LQN+1,MQN;IT+2,IU  ,IV  ]                                C
C     I2 -> E[LQN+1,MQN;IT  ,IU+2,IV  ]                                C
C     I3 -> E[LQN+1,MQN;IT  ,IU  ,IV+2]                                C
C     I4 -> E[LQN+1,MQN;IT+1,IU  ,IV  ]                                C
C     I5 -> E[LQN+1,MQN;IT  ,IU+1,IV  ]                                C
C     I6 -> E[LQN+1,MQN;IT  ,IU  ,IV+1]                                C
C     I7 -> E[LQN+1,MQN;IT  ,IU  ,IV  ]                                C
C     I8 -> E[LQN+1,MQN;IT-1,IU  ,IV  ]                                C
C     I9 -> E[LQN+1,MQN;IT  ,IU-1,IV  ]                                C
C     I10-> E[LQN+1,MQN;IT  ,IU  ,IV-1]                                C
C     I11-> E[LQN+1,MQN;IT-2,IU  ,IV  ]                                C
C     I12-> E[LQN+1,MQN;IT  ,IU-2,IV  ]                                C
C     I13-> E[LQN+1,MQN;IT  ,IU  ,IV-2]                                C
C**********************************************************************C
C
C       THE SECOND LOOP OVER ITUV INCLUDES ALL HGTF INDICES ON THE
C       LAYER CORRESPONDING TO (LQN1-1)
        RFACT1 =-DFLOAT(LQN1+IABS(MQN))/DBLE(LQN1-IABS(MQN)+1)
        DO IOUTER=0,LAM-1
          DO IT=0,IOUTER
            DO IU=0,IOUTER-IT
              IV = IOUTER-IT-IU
              I0 = ISTART0+(INABCD(IT  ,IU  ,IV  )-1)*MAXM
              I1 = ISTART2+(INABCD(IT+2,IU  ,IV  )-1)*MAXM
              I2 = ISTART2+(INABCD(IT  ,IU+2,IV  )-1)*MAXM
              I3 = ISTART2+(INABCD(IT  ,IU  ,IV+2)-1)*MAXM
              I4 = ISTART2+(INABCD(IT+1,IU  ,IV  )-1)*MAXM
              I5 = ISTART2+(INABCD(IT  ,IU+1,IV  )-1)*MAXM
              I6 = ISTART2+(INABCD(IT  ,IU  ,IV+1)-1)*MAXM
              I7 = ISTART2+(INABCD(IT  ,IU  ,IV  )-1)*MAXM
              TI = DFLOAT(2*(IT+IU+IV)+3)
              DO M=1,MAXM
                T1 = RFACT1/P22(M)
                T0 = RFACT1/P(M)
                TX = T0*PX(M)
                TY = T0*PY(M)
                TZ = T0*PZ(M)
                TT = RFACT1*(PP(M)+TI/P2(M))
                ETEMP(I1+M) = ETEMP(I1+M) + T1*ETEMP(I0+M)
                ETEMP(I2+M) = ETEMP(I2+M) + T1*ETEMP(I0+M)
                ETEMP(I3+M) = ETEMP(I3+M) + T1*ETEMP(I0+M)
                ETEMP(I4+M) = ETEMP(I4+M) + TX*ETEMP(I0+M)
                ETEMP(I5+M) = ETEMP(I5+M) + TY*ETEMP(I0+M)
                ETEMP(I6+M) = ETEMP(I6+M) + TZ*ETEMP(I0+M)
                ETEMP(I7+M) = ETEMP(I7+M) + TT*ETEMP(I0+M)
              ENDDO
              IF(IT.GE.1) THEN
                I8 = ISTART2 + (INABCD(IT-1,IU  ,IV  )-1)*MAXM
                T1 = RFACT1*DFLOAT(2*IT)
                DO M=1,MAXM
                  TX = T1*PX(M)
                  ETEMP(I8+M) = ETEMP(I8+M) + TX*ETEMP(I0+M)
                ENDDO
              ENDIF
              IF(IU.GE.1) THEN
                I9 = ISTART2 + (INABCD(IT  ,IU-1,IV  )-1)*MAXM
                T1 = RFACT1*DFLOAT(2*IU)
                DO M=1,MAXM
                  TY = T1*PY(M)
                  ETEMP(I9+M) = ETEMP(I9+M) + TY*ETEMP(I0+M)
                ENDDO
              ENDIF
              IF(IV.GE.1) THEN
                I10 = ISTART2 + (INABCD(IT  ,IU  ,IV-1)-1)*MAXM
                T1  = RFACT1*DFLOAT(2*IV)
                DO M=1,MAXM
                  TZ = T1*PZ(M)
                  ETEMP(I10+M) = ETEMP(I10+M) + TZ*ETEMP(I0+M)
                ENDDO
              ENDIF
              IF(IT.GE.2) THEN
                I11 = ISTART2 + (INABCD(IT-2,IU  ,IV  )-1)*MAXM
                T1  = RFACT1*DFLOAT(IT*(IT-1))
                DO M=1,MAXM
                  ETEMP(I11+M) = ETEMP(I11+M) + T1*ETEMP(I0+M)
                ENDDO
              ENDIF
              IF(IU.GE.2) THEN
                I12 = ISTART2 + (INABCD(IT  ,IU-2,IV  )-1)*MAXM
                T1  = RFACT1*DFLOAT(IU*(IU-1))
                DO M=1,MAXM
                  ETEMP(I12+M) = ETEMP(I12+M) + T1*ETEMP(I0+M)
                ENDDO
              ENDIF
              IF(IV.GE.2) THEN
                I13 = ISTART2 + (INABCD(IT  ,IU  ,IV-2)-1)*MAXM
                T1  = RFACT1*DFLOAT(IV*(IV-1))
                DO M=1,MAXM
                  ETEMP(I13+M) = ETEMP(I13+M) + T1*ETEMP(I0+M)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
C       END OF LOOP OVER LQN1 FOR FIXED MQN
C
        LAM     = LAM+1
        ISTART0 = ISTART
        ISTART  = ISTART + MAXM*NTUV
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE STEPN(ESG,ENSG,LAM,MAXM,ICNT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C              SSSSSS TTTTTTTT EEEEEEEE PPPPPPP  NN    NN              C
C             SS    SS   TT    EE       PP    PP NNN   NN              C
C             SS         TT    EE       PP    PP NNNN  NN              C
C              SSSSSS    TT    EEEEEE   PP    PP NN NN NN              C
C                   SS   TT    EE       PPPPPPP  NN  NNNN              C
C             SS    SS   TT    EE       PP       NN   NNN              C
C              SSSSSS    TT    EEEEEEEE PP       NN    NN              C
C                                                                      C
C -------------------------------------------------------------------- C
C  INCREMENT THE NQN, E[NQN,LQN,MQN] -> E[NQN+1,LQN,MQN].              C
C                                                                      C
C  NOTE THAT STEPN WILL ONLY PERFORM A SINGLE STEP IN NQN. IT USES AS  C
C  INPUT A SET OF PROCESSED E-COEFFICIENTS FROM VRS (ESGTFR,ESGTFI)    C
C  AND OUTPUTS THE INCREMENTED SET (ENSGTFR,ENSGTFI).                  C
C -------------------------------------------------------------------- C
C  LAM IS THE EFFECTIVE TOTAL ANGULAR MOMENTUM OF THE INPUT COEFFS.    C
C  THE EFFECTIVE TOTAL ANGULAR MOMENTUM OF THE OUTPUT COEFFS IS LAM+2. C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      COMPLEX*16 ESG(MB2,MEQ),ENSG(MB2,MEQ)
C
      DIMENSION PP(MB2),PX(MB2),PY(MB2),PZ(MB2)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/CTSN/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     IMPORT GEOMETRICAL VALUES FOR CENTRE OF INTEREST
      DO M=1,MAXM
        IF(ICNT.EQ.1) THEN
          PP(M) = PA2(M)
          PX(M) = PAX(M)
          PY(M) = PAY(M)
          PZ(M) = PAZ(M)
        ELSEIF(ICNT.EQ.2) THEN
          PP(M) = PB2(M)
          PX(M) = PBX(M)
          PY(M) = PBY(M)
          PZ(M) = PBZ(M)
        ENDIF
      ENDDO
C
C     SET THE TARGET COEFFICIENTS TO ZERO, TAKING INTO ACCOUNT THE
C     INCREMENT OF LAM BY TWO UNITS IN THE TARGET
      NTUV = (LAM+3)*(LAM+4)*(LAM+5)/6
      DO ITUV=1,NTUV
        DO M=1,MAXM
          ENSG(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C                          INDEX MAPPINGS                              C
C -------------------------------------------------------------------- C
C     I0 -> E[LQN-1,MQN;IT  ,IU  ,IV  ]                                C
C     I1 -> E[LQN+1,MQN;IT+2,IU  ,IV  ]                                C
C     I2 -> E[LQN+1,MQN;IT  ,IU+2,IV  ]                                C
C     I3 -> E[LQN+1,MQN;IT  ,IU  ,IV+2]                                C
C     I4 -> E[LQN+1,MQN;IT+1,IU  ,IV  ]                                C
C     I5 -> E[LQN+1,MQN;IT  ,IU+1,IV  ]                                C
C     I6 -> E[LQN+1,MQN;IT  ,IU  ,IV+1]                                C
C     I7 -> E[LQN+1,MQN;IT  ,IU  ,IV  ]                                C
C     I8 -> E[LQN+1,MQN;IT-1,IU  ,IV  ]                                C
C     I9 -> E[LQN+1,MQN;IT  ,IU-1,IV  ]                                C
C     I10-> E[LQN+1,MQN;IT  ,IU  ,IV-1]                                C
C     I11-> E[LQN+1,MQN;IT-2,IU  ,IV  ]                                C
C     I12-> E[LQN+1,MQN;IT  ,IU-2,IV  ]                                C
C     I13-> E[LQN+1,MQN;IT  ,IU  ,IV-2]                                C
C**********************************************************************C
C
      DO IOUTER=0,LAM
        DO IT=0,IOUTER
          DO IU=0,IOUTER-IT
            IV = IOUTER-IT-IU
C
            I0 = INABCD(IT  ,IU  ,IV  )
            I1 = INABCD(IT+2,IU  ,IV  )
            I2 = INABCD(IT  ,IU+2,IV  )
            I3 = INABCD(IT  ,IU  ,IV+2)
            I4 = INABCD(IT+1,IU  ,IV  )
            I5 = INABCD(IT  ,IU+1,IV  )
            I6 = INABCD(IT  ,IU  ,IV+1)
            I7 = INABCD(IT  ,IU  ,IV  )
C
            TT = DFLOAT((2*(IT+IU+IV))+3)
            DO M=1,MAXM
              T0 = 1.0D0/P22(M)
              TX = PX(M)/P(M)
              TY = PY(M)/P(M)
              TZ = PZ(M)/P(M)
              TP = PP(M) + TT/P2(M)
              ENSG(M,I1) = ENSG(M,I1) + T0*ESG(M,I0)
              ENSG(M,I2) = ENSG(M,I2) + T0*ESG(M,I0)
              ENSG(M,I3) = ENSG(M,I3) + T0*ESG(M,I0)
              ENSG(M,I4) = ENSG(M,I4) + TX*ESG(M,I0)
              ENSG(M,I5) = ENSG(M,I5) + TY*ESG(M,I0)
              ENSG(M,I6) = ENSG(M,I6) + TZ*ESG(M,I0)
              ENSG(M,I7) = ENSG(M,I7) + TP*ESG(M,I0)
            ENDDO
C
            IF(IT.GE.1) THEN
              RT2 = DFLOAT(2*IT)
              I8  = INABCD(IT-1,IU  ,IV  )
              DO M=1,MAXM
                T0 = PX(M)*RT2
                ENSG(M,I8) = ENSG(M,I8) + T0*ESG(M,I0)
              ENDDO
            ENDIF
C
            IF(IU.GE.1) THEN
              RU2 = DFLOAT(2*IU)
              I9   = INABCD(IT  ,IU-1,IV  )
              DO M=1,MAXM
                T0 = PY(M)*RU2
                ENSG(M,I9) = ENSG(M,I9) + T0*ESG(M,I0)
              ENDDO
            ENDIF
C
            IF(IV.GE.1) THEN
              RV2 = DFLOAT(2*IV)
              I10 = INABCD(IT  ,IU  ,IV-1)
              DO M=1,MAXM
                T0 = PZ(M)*RV2
                ENSG(M,I10) = ENSG(M,I10) + T0*ESG(M,I0)
              ENDDO
            ENDIF
C
            IF(IT.GE.2) THEN
              RT2 = DFLOAT(IT*(IT-1))
              I11 = INABCD(IT-2,IU  ,IV  )
              DO M=1,MAXM
                ENSG(M,I11) = ENSG(M,I11) + RT2*ESG(M,I0)
              ENDDO
            ENDIF
C
            IF(IU.GE.2) THEN
              RU2 = DFLOAT(IU*(IU-1))
              I12  = INABCD(IT  ,IU-2,IV  )
              DO M=1,MAXM
                ENSG(M,I12) = ENSG(M,I12) + RU2*ESG(M,I0)
              ENDDO
            ENDIF
C
            IF(IV.GE.2) THEN
              RV2 = DFLOAT(IV*(IV-1))
              I13 = INABCD(IT  ,IU  ,IV-2)
              DO M=1,MAXM
                ENSG(M,I13) = ENSG(M,I13) + RV2*ESG(M,I0)
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
      SUBROUTINE RNLL(RNORM,EXPT,LQN,NBAS)
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
C**********************************************************************C
      PARAMETER(MBS=26)
C
      DIMENSION RNORM(MBS,2),EXPT(MBS,2),LQN(2),NBAS(2)
      DATA PI,TWOLOG/3.1415926535897932D0,6.93147180559945309D-1/
C
      DO I=1,2
        T1  = DSQRT(PI)
        F1  = 0.5D0
        GML = DLOG(T1)
        DO N=2,LQN(I)+2
          GML = GML+DLOG(F1)
          F1  = F1 + 1.0D0
        ENDDO
        RLA = DFLOAT(LQN(I))
        GA1 = TWOLOG-GML
        RA1 = RLA+1.5D0
        DO IBAS=1,NBAS(I)
          ELOG          = DLOG(2.0D0*EXPT(IBAS,I))
          RNORM(IBAS,I) = DEXP(0.5D0*(GA1+RA1*ELOG))
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE RNSS(RNORM,EXPT,LQN,NBAS)
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
C**********************************************************************C
      PARAMETER(MBS=26)
C
      DIMENSION RNORM(MBS,2),EXPT(MBS,2),LQN(2),NBAS(2)
      DATA PI,TWOLOG/3.1415926535897932D0,6.93147180559945309D-1/
C
      DO I=1,2
        T1  = DSQRT(PI)
        F1  = 0.5D0
        GML = DLOG(T1)
        DO N=2,LQN(I)+3
          GML = GML+DLOG(F1)
          F1  = F1+1.0D0
        ENDDO
        RLA = DFLOAT(LQN(I))
        GA1 = TWOLOG-GML
        RA1 = RLA+0.5D0
        DO IBAS=1,NBAS(I)
          ELOG          = DLOG(2.0D0*EXPT(IBAS,I))
          RNORM(IBAS,I) = DEXP(0.5D0*(GA1+RA1*ELOG))
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE RNLS(RNORMLS,EXPT,LQN,NBAS)
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
C**********************************************************************C
      PARAMETER(MBS=26)
C
      DIMENSION RNORMLS(MBS,2),EXPT(MBS,2),LQN(2),NBAS(2)
      DATA PI,TWOLOG/3.1415926535897932D0,6.93147180559945309D-1/
C
      T1L  = DSQRT(PI)
      F1L  = 0.5D0
      GMLL = DLOG(T1L)
C
      DO N=2,LQN(1)+2
        GMLL = GMLL+DLOG(F1L)
        F1L  = F1L+1.0D0
      ENDDO
C
      RLAL = DFLOAT(LQN(1))
      GA1L = TWOLOG-GMLL
      RA1L = RLAL+1.5D0
C
      DO IBAS=1,NBAS(1)
        ELOGL           = DLOG(2.0D0*EXPT(IBAS,1))
        RNORMLS(IBAS,1) = DEXP(0.5D0*(GA1L+RA1L*ELOGL))
      ENDDO
C
      T1S  = DSQRT(PI)
      F1S  = 0.5D0
      GMLS = DLOG(T1S)
C
      DO N=2,LQN(2)+3
        GMLS = GMLS+DLOG(F1S)
        F1S  = F1S+1.0D0
      ENDDO
C
      RLAS = DFLOAT(LQN(2))
      GA1S = TWOLOG-GMLS
      RA1S = RLAS+0.5D0
C
      DO IBAS=1,NBAS(2)
        ELOGS           = DLOG(2.0D0*EXPT(IBAS,2))
        RNORMLS(IBAS,2) = DEXP(0.5D0*(GA1S+RA1S*ELOGS))
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
C     INITIATE LOOP OVER ELEMENTS OF ECMP
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

