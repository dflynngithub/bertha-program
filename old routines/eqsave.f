c     this version also produces the phase-adapted ~
c                                                  E[T,T';A,B,C]
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
      INCLUDE 'param.h'
C
      CHARACTER*4 HMLT
C
      DIMENSION NLL(0:MKP),NSS(0:MKP),NLS(0:MKP)
C
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/PRMS/HMLT,ICLC,INEW,IATM,ILIN
      COMMON/TGGL/ILEV,ICB1,ISML,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/TMMD/TELL,TESS,TELS,TESL,TRLL,TRSS,TRLS,TRSL,TRBR
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
      DO LAM=0,MKP
        NLL(LAM) = 0
        NSS(LAM) = 0
        NLS(LAM) = 0
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
      DO LAM=0,LAMMX
        IF(NLL(LAM).EQ.0) GOTO 200
        NTUVLL = (LAM+1)*(LAM+2)*(LAM+3)/6
        SPCELL = NLL(LAM)*8*8.0D-6
        WRITE(6,21) 'E0LL',LAM,NTUVLL,NLL(LAM),8,NLL(LAM)*8,SPCELL
        WRITE(7,21) 'E0LL',LAM,NTUVLL,NLL(LAM),8,NLL(LAM)*8,SPCELL
        NADETOT = NADETOT + NLL(LAM)
        NWRDTOT = NWRDTOT + NLL(LAM)*8
        SPCETOT = SPCETOT + SPCELL
200     CONTINUE
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
      IF(HMLT.EQ.'NORL') GOTO 100
C
C     E0SS ANALYSIS
      DO LAM=2,LAMMX+2
        IF(NSS(LAM).EQ.0) GOTO 210
        NTUVSS = (LAM+1)*(LAM+2)*(LAM+3)/6
        SPCESS = NSS(LAM)*8*8.0D-6
        WRITE(6,21) 'E0SS',LAM,NTUVSS,NSS(LAM),8,NSS(LAM)*8,SPCESS
        WRITE(7,21) 'E0SS',LAM,NTUVSS,NSS(LAM),8,NSS(LAM)*8,SPCESS
        NADETOT = NADETOT + NSS(LAM)
        NWRDTOT = NWRDTOT + NSS(LAM)*8
        SPCETOT = SPCETOT + SPCESS
210     CONTINUE
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
      IF(HMLT.EQ.'BARE'.OR.HMLT.EQ.'DHFR') GOTO 100
C
C     EILS ANALYSIS
      DO LAM=1,LAMMX+1
        IF(NLS(LAM).EQ.0) GOTO 220
        NTUVLS = (LAM+1)*(LAM+2)*(LAM+3)/6
        SPCELS = NSS(LAM)*24*8.0D-6
        WRITE(6,21) 'EILS',LAM,NTUVLS,NLS(LAM),24,NLS(LAM)*24,SPCELS
        WRITE(7,21) 'EILS',LAM,NTUVLS,NLS(LAM),24,NLS(LAM)*24,SPCELS
        NADETOT = NADETOT + NLS(LAM)
        NWRDTOT = NWRDTOT + NLS(LAM)*24
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
        WRITE(6, *) 'In EQSAVE: E0LL words exceed allocated limit.'
        WRITE(7, *) 'In EQSAVE: E0LL words exceed allocated limit.'
        GOTO 150
      ENDIF
C
      IF(HMLT.NE.'NORL') THEN
        IF(NADE0SS.GT.MFL) THEN
          WRITE(6, *) 'In EQSAVE: E0SS words exceed allocated limit.'
          WRITE(7, *) 'In EQSAVE: E0SS words exceed allocated limit.'
          GOTO 150
        ENDIF
      ENDIF
C
      IF(HMLT.EQ.'DHFB'.OR.HMLT.EQ.'DHFP'.OR.HMLT.EQ.'DHFQ') THEN
        IF(NADEILS.GT.MFL) THEN
          WRITE(6, *) 'In EQSAVE: EILS words exceed allocated limit.'
          WRITE(7, *) 'In EQSAVE: EILS words exceed allocated limit.'
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
      WRITE(6, *) 'In EQSAVE: E-coefficients to be generated by batch.'
      WRITE(7, *) 'In EQSAVE: E-coefficients to be generated by batch.'
C
C     FLIP THE EQ-GENERATION TOGGLE AND EXIT
      IEQS = 0
      GOTO 300
C
250   CONTINUE
C
C**********************************************************************C
C     GENERATE COMPLETE BATCHES OF EQ-COEFFS                           C
C**********************************************************************C
C
C     SECTION TITLE
      WRITE(6, *) REPEAT(' ',18),'Generating E-coefficient data files'
      WRITE(7, *) REPEAT(' ',18),'Generating E-coefficient data files'
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
