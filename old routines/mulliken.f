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
C**********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=MKP,MBS=26,MEL=100,MIT=30)
C
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM),C(MDM,MDM)
      COMPLEX*16 OLAP(MDM,MDM)
C
      DIMENSION FRAOCC1(MEL),FRAOCC2(MEL),BORDER(MEL)
C
      COMMON/COEF/C
      COMMON/ILAB/IADR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR,IOCCM0
C
C     GENERATE DIRECT OVERLAP MATRIX
      CALL MNPOLE(OLAP,0)
C
C     FIRST CENTER
      R1 = 0.0D0
      DO IOCC=1,NOCC
        MOCC = NSHIFT + IOCC
        DO ICOMP=1,2
          IADD=(ICOMP-1)*NSHIFT
          DO IR=1,IADR
            MR = IR + IADD
            DO IS=1,NSHIFT
              MS = IS + IADD
              R1 = R1 + C(MR,MOCC)*C(MS,MOCC)*OLAP(MR,MS)
            ENDDO
          ENDDO
        ENDDO
        FRAOCC1(IOCC) = R1
      ENDDO
C
C     SECOND CENTER
      R2 = 0.0D0
      DO IOCC=1,NOCC
        MOCC = NSHIFT + IOCC
        DO ICOMP=1,2
          IADD=(ICOMP-1)*NSHIFT
          DO IR=IADR+1,NSHIFT
            MR = IR + IADD
            DO IS=1,NSHIFT
              MS = IS + IADD
              R2 = R2 + C(MR,MOCC)*C(MS,MOCC)*OLAP(MR,MS)
            ENDDO
          ENDDO
        ENDDO
        FRAOCC2(IOCC) = R2
      ENDDO
C
C     FOR ALL OCCUPIED ELECTRONS, COMPUTE OVERLAP
      DO IOCC=1,NOCC
        MOCC = NSHIFT + IOCC
        TMP = 0.0D0
        DO ICOMP=1,2
          IADD=(ICOMP-1)*NSHIFT
          DO IR=1,IADR
            MR = IR + IADD
            DO IS=IADR+1,NSHIFT
              MS = IS + IADD
              TMP = TMP + C(MR,MOCC)*C(MS,MOCC)*OLAP(MR,MS)
            ENDDO
          ENDDO
        ENDDO
        BORDER(IOCC) = TMP
      ENDDO
C
C     PRINT RESULTS OF ANALYSIS
      WRITE(6, *) 'Mulliken population analysis:'
      WRITE(7, *) 'Mulliken population analysis:'
      WRITE(6,*) REPEAT('-',62)
      WRITE(7,*) REPEAT('-',62)
      WRITE(6, *) 'Total charge on center 1 =',ZNUC(1)-R1
      WRITE(7, *) 'Total charge on center 1 =',ZNUC(1)-R1
      WRITE(6, *) 'Total charge on center 2 =',ZNUC(2)-R2
      WRITE(7, *) 'Total charge on center 2 =',ZNUC(2)-R2
      WRITE(6, *) ''
      WRITE(7, *) ''
      WRITE(6, *) 'Distbtn of electrons in orbitals between centers'
      WRITE(7, *) 'Distbtn of electrons in orbitals between centers'
      DO IOCC=1,NOCC
        WRITE(6, *) IOCC,FRAOCC1(IOCC),FRAOCC2(IOCC)
        WRITE(7, *) IOCC,FRAOCC1(IOCC),FRAOCC2(IOCC)
      ENDDO
      WRITE(6, *) ''
      WRITE(7, *) ''
      WRITE(6, *) 'Overlap distribution of the orbitals:'
      WRITE(7, *) 'Overlap distribution of the orbitals:'
      DO IOCC=1,NOCC
        WRITE(6, *) IOCC,BORDER(IOCC)
        WRITE(7, *) IOCC,BORDER(IOCC)
      ENDDO
C
      RETURN
      END
