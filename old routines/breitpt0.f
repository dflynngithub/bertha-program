      SUBROUTINE BREITPT0(IZ,IA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  BBBBBBB  RRRRRRR  EEEEEEEE IIII TTTTTTTT PPPPPPP TTTTTTTT 000000    C
C  BB    BB RR    RR EE        II     TT    PP    PP   TT   00   000   C
C  BB    BB RR    RR EE        II     TT    PP    PP   TT   00  0000   C
C  BBBBBBB  RR    RR EEEEEE    II     TT    PP    PP   TT   00 00 00   C
C  BB    BB RRRRRRR  EE        II     TT    PPPPPPP    TT   0000  00   C
C  BB    BB RR    RR EE        II     TT    PP         TT   000   00   C
C  BBBBBBB  RR    RR EEEEEEEE IIII    TT    PP         TT    000000    C
C -------------------------------------------------------------------- C
C  BREITPT0 PERFORMS A FIRST-ORDER ANALYSIS OF THE BREIT INTERACTION   C
C  APPLIED TO AN AVERAGE-OVER-CONFIGURATION ATOMIC SOLUTION.           C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ IZ - ATOMIC CENTRE INDEX                                          C
C  ▶ IA - ANALYSIS MODE (1 FOR SYMMETRY TYPES, 2 FOR ATOMIC ORBITALS)  C
C**********************************************************************C
      INCLUDE 'param.h'
C
      CHARACTER*1 LLAB,QSGN
      CHARACTER*2 ELMT(120),ELNM,KLAB
      CHARACTER*4 HMLT
C
      DIMENSION QA(MKP),QE(MKP)
      DIMENSION CL(MBS,MKP,MKP+1),CS(MBS,MKP,MKP+1)
      DIMENSION DALL(MB2,MKP,MKP+1),DALS(MB2,MKP,MKP+1),
     &          DASL(MB2,MKP,MKP+1),DASS(MB2,MKP,MKP+1)
      DIMENSION DELL(MB2,MKP,MKP+1),DELS(MB2,MKP,MKP+1),
     &          DESL(MB2,MKP,MKP+1),DESS(MB2,MKP,MKP+1)
C
      DIMENSION DNALL(MB2,MKP),DNALS(MB2,MKP),
     &          DNASL(MB2,MKP),DNASS(MB2,MKP)
      DIMENSION DNELL(MB2,MKP),DNELS(MB2,MKP),
     &          DNESL(MB2,MKP),DNESS(MB2,MKP)
      DIMENSION EBSTT(MKP,MKP,4),EBS(4)
      DIMENSION EBOTT(MKP,MKP,MKP+1,MKP+1,4),EBO(4)
C
      COMPLEX*16 COEF(MDM,MDM)
C
      COMMON/ATMB/B11(MBD,MBD),B21(MBD,MBD),B12(MBD,MBD),B22(MBD,MBD)
      COMMON/ATMD/DLL1(MB2),DSL1(MB2),DSS1(MB2),DLS1(MB2),
     &            DLL2(MB2),DSL2(MB2),DSS2(MB2),DLS2(MB2)
      COMMON/AUFB/NORB(0:MEL,MKP+1),NUMOCC(0:MEL),LMXCONF
      COMMON/B0QN/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),CFMI(MCT),
     &            RNUC(MCT),FNUC(MCT,MFT),XNUC(MCT,MFT),NNUC(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EWDR,EWXC,EUEH
      COMMON/FILL/NCNF(MCT,0:MEL,MKP+1),NLVL(MCT,0:MEL),IAUF(MCT)
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
      COMMON/PRMS/HMLT,IOPT,IMOL,INEW,ILEV,ISML,ISWZ,IEQS,IERC,IPAR,ICOR
C
C**********************************************************************C
C     PREPARE DENSITY MATRIX FOR EACH ORBITAL                          C
C**********************************************************************C
C
C     LOOP OVER ALL OCCUPIED SYMMETRY TYPES
      DO IKP=1,2*LMXCONF+1
C
C       QUANTUM NUMBERS
        KQNA = KAPA(IKP,IZ)
        IF(KQNA.LT.0) THEN
          LQNA =-KQNA-1
        ELSE
          LQNA = KQNA
        ENDIF
C
C       AVERAGE FRACTIONAL SUBSHELL OCCUPANCY
        DO NQN=1,NUMOCC(LQNA)
          QA(NQN) = DFLOAT(NORB(LQNA,NQN))/DFLOAT(4*LQNA+2)
        ENDDO
C
C       EFFECTIVE FRACTIONAL SUBSHELL OCCUPANCY
        DO NQN=1,NUMOCC(LQNA)
          IF(NORB(LQNA,NQN).EQ.(4*LQNA+2)) THEN
            QE(NQN) = 1.0D0
          ELSE
            QE(NQN) = DFLOAT(NORB(LQNA,NQN)-1)/DFLOAT(4*LQNA+1)
          ENDIF
        ENDDO
C
C       LOOP OVER ALL OCCUPIED LEVELS
        DO NQN=1,NUMOCC(LQNA)
C
C         SKIP ANY UNOCCUPIED LEVELS
          IF(NORB(LQNA,NQN).NE.0) THEN
C
            QAF = DSQRT(QA(NQN))
C
C           COEFFICIENT MATRIX STARTING ADDRESSES
            IL1 = LRGE(IZ,IKP,1)
            IS1 = LRGE(IZ,IKP,1)+NSKP
          
            IROW = IOCTP(IZ,IKP,1,NQN)
C
C           LOOP OVER BASIS FUNCTIONS
            DO IBAS=1,NFNC(LQNA,IZ)
C
C             IMPORT COEFFICIENTS
              CL(IBAS,IKP,NQN) = DREAL(COEF(IL1+IBAS,IROW+NSKP))/QAF
              CS(IBAS,IKP,NQN) = DREAL(COEF(IS1+IBAS,IROW+NSKP))/QAF
C          
            ENDDO
C        
          ENDIF
C      
        ENDDO
C
C       BUILD ATOMIC CHARGE DENSITY MATRIX FOR THIS KQNA BLOCK
        M = 0
        DO IBAS=1,NFNC(LQNA,IZ)
          DO JBAS=1,NFNC(LQNA,IZ)
            M = M+1
C
C           SUM COUNTERS
            DNALL(M,IKP) = 0.0D0
            DNALS(M,IKP) = 0.0D0
            DNASL(M,IKP) = 0.0D0
            DNASS(M,IKP) = 0.0D0
C
C           SUM COUNTERS
            DNELL(M,IKP) = 0.0D0
            DNELS(M,IKP) = 0.0D0
            DNESL(M,IKP) = 0.0D0
            DNESS(M,IKP) = 0.0D0
C
C           LOOP OVER ALL LISTED SUBSHELLS OF THIS KQN TYPE
            DO NQN=1,NUMOCC(LQNA)
C
              QAF = QA(NQN)
              QEF = QE(NQN)
C
              DALL(M,IKP,NQN) = QAF*CL(IBAS,IKP,NQN)*CL(JBAS,IKP,NQN)
              DALS(M,IKP,NQN) = QAF*CL(IBAS,IKP,NQN)*CS(JBAS,IKP,NQN)
              DASL(M,IKP,NQN) = QAF*CS(IBAS,IKP,NQN)*CL(JBAS,IKP,NQN)
              DASS(M,IKP,NQN) = QAF*CS(IBAS,IKP,NQN)*CS(JBAS,IKP,NQN)
C
C             UPDATE NET DENSITY
              DNALL(M,IKP) = DNALL(M,IKP) + DALL(M,IKP,NQN)
              DNALS(M,IKP) = DNALS(M,IKP) + DALS(M,IKP,NQN)
              DNASL(M,IKP) = DNASL(M,IKP) + DASL(M,IKP,NQN)
              DNASS(M,IKP) = DNASS(M,IKP) + DASS(M,IKP,NQN)
C
              DELL(M,IKP,NQN) = QEF*CL(IBAS,IKP,NQN)*CL(JBAS,IKP,NQN)
              DELS(M,IKP,NQN) = QEF*CL(IBAS,IKP,NQN)*CS(JBAS,IKP,NQN)
              DESL(M,IKP,NQN) = QEF*CS(IBAS,IKP,NQN)*CL(JBAS,IKP,NQN)
              DESS(M,IKP,NQN) = QEF*CS(IBAS,IKP,NQN)*CS(JBAS,IKP,NQN)
C
C             UPDATE NET DENSITY
              DNELL(M,IKP) = DNELL(M,IKP) + DELL(M,IKP,NQN)
              DNELS(M,IKP) = DNELS(M,IKP) + DELS(M,IKP,NQN)
              DNESL(M,IKP) = DNESL(M,IKP) + DESL(M,IKP,NQN)
              DNESS(M,IKP) = DNESS(M,IKP) + DESS(M,IKP,NQN)
C
            ENDDO

          ENDDO
        ENDDO
C
C     END LOOP OVER KQNA VALUES
      ENDDO
C
C**********************************************************************C
C     BREIT ENERGY BY SYMMETRY TYPE PAIRS                              C
C**********************************************************************C
C
      IF(IA.EQ.2) GOTO 699
C
      WRITE(6,*) REPEAT('=',72)
      WRITE(7,*) REPEAT('=',72)
      WRITE(6,*) REPEAT(' ',17),'Atomic symmetry type Breit interaction'
      WRITE(7,*) REPEAT(' ',17),'Atomic symmetry type Breit interaction'
      WRITE(6,*) REPEAT('=',72)
      WRITE(7,*) REPEAT('=',72)
601   FORMAT(1X,'(κA|κB)',9X,
     &              'E(LL|SS)',8X,'E(SS|LL)',8X,'E(SL|LS)',10X,'E(TOT)')
602   FORMAT(1X,'(',A,'|',A,')',3X,
     &                             F14.10,2X,F14.10,2X,F14.10,2X,F14.10)
603   FORMAT(1X,'Atomic',4X,F14.10,2X,F14.10,2X,F14.10,2X,F14.10)
C
C     PRINT HEADER TO TERMINAL
      WRITE(6,601)
      WRITE(7,601)
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
C
C     INITIALISE BREIT COUNTERS
      DO ITT=1,4
        EBS(ITT) = 0.0D0
        DO KPA=1,MKP
          DO KPB=1,MKP
            EBSTT(KPA,KPB,ITT) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     LOOP OVER ALL OCCUPIED LQN VALUES
      DO LQNA=0,LMXCONF
C
C       RECORD LQNA VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
        NBASA = NFNC(LQNA,IZ)
        DO IBAS=1,NBASA
          EXLA(IBAS) = BEXL(IBAS,LQNA,IZ)
        ENDDO
C
C       LOOP OVER ALL OCCUPIED LQN VALUES
        DO LQNB=0,LMXCONF
C
C         RECORD LQNB VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
          NBASB = NFNC(LQNB,IZ)
          DO JBAS=1,NBASB
            EXLB(JBAS) = BEXL(JBAS,LQNB,IZ)
          ENDDO
C
C         NUMBER OF BASIS FUNCTION OVERLAPS IN THIS BLOCK
          MAXM = NBASB*NBASB
C
C         RELEVANT DENSITY MATRIX DEPENDS ON LQNA AND LQNB
          IF(LQNA.NE.LQNB) THEN
            DO M=1,MAXM
              DLL2(M) = DNALL(M,2*LQNB+1)
              DSL2(M) = DNASL(M,2*LQNB+1)
              DSS2(M) = DNASS(M,2*LQNB+1)
              DLS2(M) = DNALS(M,2*LQNB+1)
              IF(LQNB.GT.0) THEN
                DLL1(M) = DNALL(M,2*LQNB  )
                DSL1(M) = DNASL(M,2*LQNB  )
                DSS1(M) = DNASS(M,2*LQNB  )
                DLS1(M) = DNALS(M,2*LQNB  )
              ENDIF
            ENDDO
          ELSE
            DO M=1,MAXM
              DLL2(M) = DNELL(M,2*LQNB+1)
              DSL2(M) = DNESL(M,2*LQNB+1)
              DSS2(M) = DNESS(M,2*LQNB+1)
              DLS2(M) = DNELS(M,2*LQNB+1)
              IF(LQNB.GT.0) THEN
                DLL1(M) = DNELL(M,2*LQNB  )
                DSL1(M) = DNESL(M,2*LQNB  )
                DSS1(M) = DNESS(M,2*LQNB  )
                DLS1(M) = DNELS(M,2*LQNB  )
              ENDIF
            ENDDO
          ENDIF
C
C         GENERATE THE MEAN-FIELD ATOMIC BREIT MATRIX OVER DENSITIES
          CALL BREIT0
C
C         LOOP OVER SYMMETRY TYPES FOR THIS LQNA
          DO ISYMA=0,1
C
            KQNA = ((-1)**ISYMA)*LQNA - ISYMA
            RK2A = DFLOAT(2*IABS(KQNA))
            KPA  = 2*LQNA+ISYMA
C
            IF(KQNA.EQ.0) GOTO 651
C
C           LOOP OVER SYMMETRY TYPES FOR THIS LQNB
            DO ISYMB=0,1
C
              KQNB = ((-1)**ISYMB)*LQNB - ISYMB
              RK2B = DFLOAT(2*IABS(KQNB))
              KPB  = 2*LQNB+ISYMB
C
              IF(KQNB.EQ.0) GOTO 652
C
C             SUM OVER THIS BREIT MATRIX AND MULTIPLY BY DENSITY
              M = 0
              DO IBAS=1,NBASA
                DO JBAS=1,NBASA
                  M = M+1
C
                  IF(ISYMA.EQ.1.AND.ISYMB.EQ.1) THEN
                    BLL = B22(IBAS      ,JBAS      )
                    BSS = B22(IBAS+NBASA,JBAS+NBASA)
                    BSL = B22(IBAS+NBASA,JBAS      )
                  ELSEIF(ISYMA.EQ.1.AND.ISYMB.EQ.0) THEN
                    BLL = B21(IBAS      ,JBAS      )
                    BSS = B21(IBAS+NBASA,JBAS+NBASA)
                    BSL = B21(IBAS+NBASA,JBAS      )
                  ELSEIF(ISYMA.EQ.0.AND.ISYMB.EQ.1) THEN
                    BLL = B12(IBAS      ,JBAS      )
                    BSS = B12(IBAS+NBASA,JBAS+NBASA)
                    BSL = B12(IBAS+NBASA,JBAS      )
                  ELSEIF(ISYMA.EQ.0.AND.ISYMB.EQ.0) THEN
                    BLL = B11(IBAS      ,JBAS      )
                    BSS = B11(IBAS+NBASA,JBAS+NBASA)
                    BSL = B11(IBAS+NBASA,JBAS      )
                  ENDIF
C
                  EBSTT(KPA,KPB,1) = EBSTT(KPA,KPB,1) + BLL*DNALL(M,KPA)
                  EBSTT(KPA,KPB,2) = EBSTT(KPA,KPB,2) + BSS*DNASS(M,KPA)
                  EBSTT(KPA,KPB,3) = EBSTT(KPA,KPB,3) + BSL*DNASL(M,KPA)
C
                ENDDO
              ENDDO
C
C             MULTIPLY BY SUBSHELL DEGENERACY FACTOR
              DO ITT=1,3
                EBSTT(KPA,KPB,ITT) = 0.5D0*RK2A*RK2B*EBSTT(KPA,KPB,ITT)
              ENDDO
              EBSTT(KPA,KPB,4) = EBSTT(KPA,KPB,1) + EBSTT(KPA,KPB,2)
     &                                       +2.0D0*EBSTT(KPA,KPB,3)
C
              DO ITT=1,4
                EBS(ITT) = EBS(ITT) + EBSTT(KPA,KPB,ITT)
              ENDDO
C
652           CONTINUE
              ENDDO
C
651       CONTINUE
          ENDDO
C
C       END LOOP OVER LQNS FOR ORBITAL B
        ENDDO
C
C     END LOOP OVER LQNS FOR ORBITAL A
      ENDDO
C
C     QUOTE THE FINAL RESULTS IN A TABLE
      DO IKP=1,2*LMXCONF+1
        DO JKP=1,2*LMXCONF+1
          WRITE(6,602) KLAB(KAPA(IKP,IZ)),KLAB(KAPA(JKP,IZ)),
     &                                    (EBSTT(IKP,JKP,ITT),ITT=1,4)
          WRITE(7,602) KLAB(KAPA(IKP,IZ)),KLAB(KAPA(JKP,IZ)),
     &                                    (EBSTT(IKP,JKP,ITT),ITT=1,4)
        ENDDO
      ENDDO
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
      WRITE(6,603) (EBS(ITT),ITT=1,4)
      WRITE(7,603) (EBS(ITT),ITT=1,4)
      WRITE(6,*) REPEAT('=',72)
      WRITE(7,*) REPEAT('=',72)
      WRITE(6,*) ' '
      WRITE(7,*) ' '
C
C     RECORD TOTAL BREIT ENERGY
      EB = EBS(4)
C
699   CONTINUE
C
C**********************************************************************C
C     BREIT ENERGY BY ATOMIC ORBITAL PAIRS                             C
C**********************************************************************C
C
      IF(IA.EQ.1) GOTO 799
C
      WRITE(6,*) REPEAT('=',76)
      WRITE(7,*) REPEAT('=',76)
      WRITE(6,*) REPEAT(' ',22),'Atomic orbital Breit interaction'
      WRITE(7,*) REPEAT(' ',22),'Atomic orbital Breit interaction'
      WRITE(6,*) REPEAT('=',76)
      WRITE(7,*) REPEAT('=',76)
701   FORMAT(1X,'(n κA|n κB)',9X,
     &              'E(LL|SS)',8X,'E(SS|LL)',8X,'E(SL|LS)',10X,'E(TOT)')
702   FORMAT(1X,'(',I2,A,'|',I2,A,')',3X,
     &                             F14.10,2X,F14.10,2X,F14.10,2X,F14.10)
703   FORMAT(1X,'Atomic',8X,F14.10,2X,F14.10,2X,F14.10,2X,F14.10)
C
C     PRINT HEADER TO TERMINAL
      WRITE(6,701)
      WRITE(7,701)
      WRITE(6,*) REPEAT('-',76)
      WRITE(7,*) REPEAT('-',76)
C
C     INITIALISE BREIT COUNTERS
      DO ITT=1,4
        EBO(ITT) = 0.0D0
        DO KPA=1,MKP
          DO KPB=1,MKP
            DO IOCC=1,MKP+1
              DO JOCC=1,MKP+1
                EBOTT(KPA,KPB,IOCC,JOCC,ITT) = 0.0D0
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     LOOP OVER ALL OCCUPIED LQN VALUES
      DO LQNA=0,LMXCONF
C
C       RECORD LQNA VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
        NBASA = NFNC(LQNA,IZ)
        DO IBAS=1,NBASA
          EXLA(IBAS) = BEXL(IBAS,LQNA,IZ)
        ENDDO
C
C       LOOP OVER ALL OCCUPIED LQN VALUES
        DO LQNB=0,LMXCONF
C
C         RECORD LQNB VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
          NBASB = NFNC(LQNB,IZ)
          DO JBAS=1,NBASB
            EXLB(JBAS) = BEXL(JBAS,LQNB,IZ)
          ENDDO
          
          DO JOCC=1,NUMOCC(LQNB)
C
C         NUMBER OF BASIS FUNCTION OVERLAPS IN THIS BLOCK
          MAXM = NBASB*NBASB
C
C         RELEVANT DENSITY MATRIX DEPENDS ON LQNA AND LQNB
          IF(LQNA.NE.LQNB) THEN
            DO M=1,MAXM
              DLL2(M) = DALL(M,2*LQNB+1,JOCC)
              DSL2(M) = DASL(M,2*LQNB+1,JOCC)
              DSS2(M) = DASS(M,2*LQNB+1,JOCC)
              DLS2(M) = DALS(M,2*LQNB+1,JOCC)
              IF(LQNB.GT.0) THEN
                DLL1(M) = DALL(M,2*LQNB  ,JOCC)
                DSL1(M) = DASL(M,2*LQNB  ,JOCC)
                DSS1(M) = DASS(M,2*LQNB  ,JOCC)
                DLS1(M) = DALS(M,2*LQNB  ,JOCC)
              ENDIF
            ENDDO
          ELSE
            DO M=1,MAXM
              DLL2(M) = DELL(M,2*LQNB+1,JOCC)
              DSL2(M) = DESL(M,2*LQNB+1,JOCC)
              DSS2(M) = DESS(M,2*LQNB+1,JOCC)
              DLS2(M) = DELS(M,2*LQNB+1,JOCC)
              IF(LQNB.GT.0) THEN
                DLL1(M) = DELL(M,2*LQNB  ,JOCC)
                DSL1(M) = DESL(M,2*LQNB  ,JOCC)
                DSS1(M) = DESS(M,2*LQNB  ,JOCC)
                DLS1(M) = DELS(M,2*LQNB  ,JOCC)
              ENDIF
            ENDDO
          ENDIF
C
C         GENERATE THE MEAN-FIELD ATOMIC BREIT MATRIX OVER DENSITIES
          CALL BREIT0
C
          DO IOCC=1,NUMOCC(LQNA)
C
C         LOOP OVER SYMMETRY TYPES FOR THIS LQNA
          DO ISYMA=0,1
C
            KQNA = ((-1)**ISYMA)*LQNA - ISYMA
            RK2A = DFLOAT(2*IABS(KQNA))
            KPA  = 2*LQNA+ISYMA
C
            IF(KQNA.EQ.0) GOTO 751
C
C           LOOP OVER SYMMETRY TYPES FOR THIS LQNB
            DO ISYMB=0,1
C
              KQNB = ((-1)**ISYMB)*LQNB - ISYMB
              RK2B = DFLOAT(2*IABS(KQNB))
              KPB  = 2*LQNB+ISYMB
C
              IF(KQNB.EQ.0) GOTO 752
C
C             SUM OVER THIS BREIT MATRIX AND MULTIPLY BY DENSITY
              M = 0
              DO IBAS=1,NBASA
                DO JBAS=1,NBASA
                  M = M+1
C
                  IF(ISYMA.EQ.1.AND.ISYMB.EQ.1) THEN
                    BLL = B22(IBAS      ,JBAS      )
                    BSS = B22(IBAS+NBASA,JBAS+NBASA)
                    BSL = B22(IBAS+NBASA,JBAS      )
                  ELSEIF(ISYMA.EQ.1.AND.ISYMB.EQ.0) THEN
                    BLL = B21(IBAS      ,JBAS      )
                    BSS = B21(IBAS+NBASA,JBAS+NBASA)
                    BSL = B21(IBAS+NBASA,JBAS      )
                  ELSEIF(ISYMA.EQ.0.AND.ISYMB.EQ.1) THEN
                    BLL = B12(IBAS      ,JBAS      )
                    BSS = B12(IBAS+NBASA,JBAS+NBASA)
                    BSL = B12(IBAS+NBASA,JBAS      )
                  ELSEIF(ISYMA.EQ.0.AND.ISYMB.EQ.0) THEN
                    BLL = B11(IBAS      ,JBAS      )
                    BSS = B11(IBAS+NBASA,JBAS+NBASA)
                    BSL = B11(IBAS+NBASA,JBAS      )
                  ENDIF
C
                  EBOTT(KPA,KPB,IOCC,JOCC,1)
     &              = EBOTT(KPA,KPB,IOCC,JOCC,1) + BLL*DALL(M,KPA,IOCC)
                  EBOTT(KPA,KPB,IOCC,JOCC,2)
     &              = EBOTT(KPA,KPB,IOCC,JOCC,2) + BSS*DASS(M,KPA,IOCC)
                  EBOTT(KPA,KPB,IOCC,JOCC,3)
     &              = EBOTT(KPA,KPB,IOCC,JOCC,3) + BSL*DASL(M,KPA,IOCC)
C
                ENDDO
              ENDDO
C
C             MULTIPLY BY SUBSHELL DEGENERACY FACTOR
              DO ITT=1,3
                EBOTT(KPA,KPB,IOCC,JOCC,ITT)
     &                  = 0.5D0*RK2A*RK2B*EBOTT(KPA,KPB,IOCC,JOCC,ITT)
              ENDDO
              EBOTT(KPA,KPB,IOCC,JOCC,4) = EBOTT(KPA,KPB,IOCC,JOCC,1)
     &                                   + EBOTT(KPA,KPB,IOCC,JOCC,2)
     &                              +2.0D0*EBOTT(KPA,KPB,IOCC,JOCC,3)
C
              DO ITT=1,4
                EBO(ITT) = EBO(ITT) + EBOTT(KPA,KPB,IOCC,JOCC,ITT)
              ENDDO
C
752           CONTINUE
              ENDDO
C
751       CONTINUE
          ENDDO
C
C       NQNA
        ENDDO
C
C         NQNB
          ENDDO
C
C       END LOOP OVER LQNS FOR ORBITAL B
        ENDDO
C
C     END LOOP OVER LQNS FOR ORBITAL A
      ENDDO
C
C     LIST OF LQN VALUES FROM IKP
C
C     QUOTE THE FINAL RESULTS IN A TABLE
      DO IKP=1,2*LMXCONF+1
        KQNA = KAPA(IKP,IZ)
        IF(KQNA.LT.0) THEN
          LQNA =-KQNA-1
        ELSE
          LQNA = KQNA
        ENDIF
        DO IOCC=1,NUMOCC(LQNA)
          DO JKP=1,2*LMXCONF+1
            KQNB = KAPA(JKP,IZ)
            IF(KQNB.LT.0) THEN
              LQNB =-KQNB-1
            ELSE
              LQNB = KQNB
            ENDIF
            DO JOCC=1,NUMOCC(LQNB)
              WRITE(6,702) IOCC+LQNA,KLAB(KAPA(IKP,IZ)),
     &                     JOCC+LQNB,KLAB(KAPA(JKP,IZ)),
     &                     (EBOTT(IKP,JKP,IOCC,JOCC,ITT),ITT=1,4)
              WRITE(7,702) IOCC+LQNA,KLAB(KAPA(IKP,IZ)),
     &                     JOCC+LQNB,KLAB(KAPA(JKP,IZ)),
     &                     (EBOTT(IKP,JKP,IOCC,JOCC,ITT),ITT=1,4)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      WRITE(6,*) REPEAT('-',76)
      WRITE(7,*) REPEAT('-',76)
      WRITE(6,703) (EBO(ITT),ITT=1,4)
      WRITE(7,703) (EBO(ITT),ITT=1,4)
      WRITE(6,*) REPEAT('=',76)
      WRITE(7,*) REPEAT('=',76)
      WRITE(6,*) ' '
      WRITE(7,*) ' '
C
C     RECORD TOTAL BREIT ENERGY
      EB = EBO(4)
C
799   CONTINUE
C
C**********************************************************************C
C     ALL DONE.                                                        C
C**********************************************************************C
C
      RETURN
      END
