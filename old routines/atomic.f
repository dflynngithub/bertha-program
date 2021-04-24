      SUBROUTINE ATOMIC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             AA   TTTTTTTT OOOOOO  MM       MM IIII CCCCCC            C
C            AAAA     TT   OO    OO MMM     MMM  II CC    CC           C
C           AA  AA    TT   OO    OO MMMM   MMMM  II CC                 C
C          AA    AA   TT   OO    OO MM MM MM MM  II CC                 C
C          AAAAAAAA   TT   OO    OO MM  MMM  MM  II CC                 C
C          AA    AA   TT   OO    OO MM   M   MM  II CC    CC           C
C          AA    AA   TT    OOOOOO  MM       MM IIII CCCCCC            C
C                                                                      C
C                       CONTROLLING ROUTINE FOR                        C
C                          ** B E R T H A **                           C
C -------------------------------------------------------------------- C
C  ATOMIC PERFORMS A SINGLE-CENTRE SCF PROCEDURE FOR EACH ATOM IN THE  C
C  MOLECULE AND ASSEMBLES AN INITIAL COEFFICIENT MATRIX.               C
C**********************************************************************C
      INCLUDE 'param.h'
C
      CHARACTER*4  HMLT
      CHARACTER*16 HMS
      CHARACTER*20 STAMP
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 FOCK(MDM,MDM),OVAP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VUEH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),WDIR(MDM,MDM),
     &           WXCH(MDM,MDM),CPLE(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EWDR,EWXC,EUEH
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/MTRX/FOCK,OVAP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VUEH,QDIR,
     &            QXCH,WDIR,WXCH,CPLE
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
      COMMON/PRMS/HMLT,IOPT,IMOL,INEW,ILEV,ISML,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SHLL/ACFF,BCFF,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
      COMMON/TATM/TTOT
C
C     PRINT A BIG SECTION HEADER
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) REPEAT(' ',26),'ATOMIC HARTREE-FOCK SCF'
      WRITE(7, *) REPEAT(' ',26),'ATOMIC HARTREE-FOCK SCF'
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C     RECORD TIME AT START OF ATOMIC CALCULATION
      CALL CPU_TIME(TDUM)
C
C     INITIALISE MOLECULAR MATRICES AND RADIAL EXPECTATION VALUES
      DO I=1,NDIM
        DO J=1,NDIM
          COEF(I,J) = DCMPLX(0.0D0,0.0D0)
          OVAP(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     INITIALISE OCCUPATION COUNTER
      IOCCM0 = 0
C
C     SPECIAL EXIT FOR ZERO- AND ONE-ELECTRON PROBLEMS
      IF(NOCC.EQ.0) THEN
        WRITE(6, *) 'There are no electrons! Skip ATOMIC.'
        WRITE(7, *) 'There are no electrons! Skip ATOMIC.'
        RETURN
      ELSEIF(NOCC.EQ.1.AND.NCNT.EQ.1) THEN
        CALL ATOM1E(NCNT)
        CALL DENSTY0
        CALL SPECTRM0
        GOTO 300
      ENDIF
C
C     HARTREE-FOCK SCF PROCEDURE FOR EACH ISOLATED ATOM
      DO IZ=1,NCNT
        CALL HFSCF0(IZ)
        IF(HMLT.EQ.'DHFP') THEN
C         MODE 1 IS BY SYMMETRY TYPE PAIRS, MODE 2 IS BY ATOMIC ORBITALS
          CALL G2PAIR0(IZ,3,'BREIT')
        ENDIF
      ENDDO
C
C     GENERATE DENSITY MATRIX
      CALL DENSTY0
C
C     PRODUCE SPECTRUM SUMMARY
      CALL SPECTRM0
C
C     MOLECULAR ENERGIES
20    FORMAT(1X,A,31X,F21.12)
      WRITE(6, *) REPEAT(' ',19),'Molecular energies (Hartree units)'
      WRITE(7, *) REPEAT(' ',19),'Molecular energies (Hartree units)'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) 'Source',REPEAT(' ',60),'Energy'
      WRITE(7, *) 'Source',REPEAT(' ',60),'Energy'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) 'Nucleus-nucleus  (N)',ENUC
      WRITE(7,20) 'Nucleus-nucleus  (N)',ENUC
      WRITE(6,20) 'One-electron     (H)',EONE
      WRITE(7,20) 'One-electron     (H)',EONE
      IF(HMLT.EQ.'BARE') GOTO 202
      WRITE(6,20) 'Coulomb (closed) (G)',ECLG
      WRITE(7,20) 'Coulomb (closed) (G)',ECLG
      IF(HMLT.EQ.'NORL'.OR.HMLT.EQ.'DHFR') GOTO 202
      WRITE(6,20) 'Breit (closed)   (B)',EBRG
      WRITE(7,20) 'Breit (closed)   (B)',EBRG
      IF(HMLT.NE.'DHFQ') GOTO 202
      WRITE(6,20) 'Uehling          (U)',EUEH
      WRITE(7,20) 'Uehling          (U)',EUEH
202   CONTINUE
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) 'Molecule total   (F)',ETOT
      WRITE(7,20) 'Molecule total   (F)',ETOT
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
300   CONTINUE
C
C     SAVE EIGENVECTORS TO OUTPUT FILE
      OPEN(UNIT=8,FILE=TRIM(WFNFL),STATUS='UNKNOWN')
      REWIND(UNIT=10)
      DO I=1,NDIM
        WRITE(8, *) EIGN(I),(COEF(J,I),J=1,NDIM)
      ENDDO
      CLOSE(UNIT=8)
C
C     TIME TAKEN FOR ATOMIC CALCULATION
      CALL CPU_TIME(TTOT)
      TTOT = TTOT - TDUM
C
C     DATE AND TIME AT END OF ITERATION
      CALL TIMENOW(STAMP)
C
C     CALCULATION TIME
30    FORMAT(1X,A,26X,A)
      WRITE(6,30) 'Total atomic SCF time         ',HMS(TTOT)
      WRITE(7,30) 'Total atomic SCF time         ',HMS(TTOT)
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,30) 'Time at end of calculation',STAMP
      WRITE(7,30) 'Time at end of calculation',STAMP
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
      RETURN
      END
