      PROGRAM BERTHA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         BBBBBBB  EEEEEEEE RRRRRRR TTTTTTTT HH    HH    AA            C
C         BB    BB EE       RR    RR   TT    HH    HH   AAAA           C
C         BB    BB EE       RR    RR   TT    HH    HH  AA  AA          C
C         BBBBBBB  EEEEEE   RR    RR   TT    HHHHHHHH AA    AA         C
C         BB    BB EE       RRRRRRR    TT    HH    HH AAAAAAAA         C
C         BB    BB EE       RR    RR   TT    HH    HH AA    AA         C
C         BBBBBBB  EEEEEEEE RR    RR   TT    HH    HH AA    AA         C
C                                                                      C
C                 (THE PROGRAM FORMERLY KNOWN AS...)                   C
C                                                                      C
C     SSSSSS  WW         WW IIII RRRRRRR  LL      EEEEEEEE  SSSSSS     C
C    SS    SS WW         WW  II  RR    RR LL      EE       SS    SS    C
C    SS       WW         WW  II  RR    RR LL      EE       SS          C
C     SSSSSS  WW    W    WW  II  RR    RR LL      EEEEEE    SSSSSS     C
C          SS WW   WWW   WW  II  RRRRRRR  LL      EE             SS    C
C    SS    SS  WW WW WW WW   II  RR    RR LL      EE       SS    SS    C
C     SSSSSS    WW     WW   IIII RR    RR LLLLLLL EEEEEEEE  SSSSSS     C
C                                                                      C
C -------------------------------------------------------------------- C
C        A RELATIVISTIC MOLECULAR ELECTRONIC STRUCTURE PROGRAM         C
C            BASED ON THE ANALYTIC FINITE BASIS SET METHOD.            C
C                                                                      C
C        (c)   H.M.QUINEY, H.SKAANE, I.P.GRANT (OXFORD, 1996)          C
C              D. FLYNN (UNIMELB, 2017)                                C
C ******************************************************************** C
C                         HAMILTONIANS (HMLTN)                         C
C                         --------------------                         C
C   'NORL' NON-RELATIVISTIC HAMILTONIAN (PAULI EQUATION).              C
C   'BARE' BARE NUCLEUS DIRAC HAMILTONIAN (NO ELECTRON INTERACTION).   C
C   'DHFR' DIRAC-COULOMB HAMILTONIAN.                                  C
C   'DHFP' DIRAC-COULOMB HAMILTONIAN (+1ST ORDER BREIT).               C
C   'DHFB' DIRAC-COULOMB-BREIT HAMILTONIAN.                            C
C   'DHFQ' DIRAC-COULOMB-BREIT HAMILTONIAN WITH LEADING-ORDER QED.     C
C ******************************************************************** C
C                      CALCULATION TREES (ITREE)                       C
C                      -------------------------                       C
C    I: HARTREE-FOCK SCF CALCULATION.                                  C
C   II: MANY-BODY DIAGRAMMATIC PERTURBATION THEORY.                    C
C  III: MULTI-CONFIGURATIONAL SCF CALCULATION.                         C
C   IV: DENSITY MATRIX RENORMALISATION GROUP CALCULATION.              C
C    V: CALCULATION OF MOLECULAR EXPECTATION VALUES.                   C
C   VI: VISUALS (ELECTROMAGNETIC FIELDS, POTENTIALS, AMPLITUDES ETC).  C
C ******************************************************************** C
C                          TABLE OF CONTENTS                           C
C                          -----------------                           C
C   [1] INPUT/OUTPUT: READ FROM INPUT FILE AND SUMMARISE DATA.         C
C   [2] LABELS: MOLECULAR GEOMETRY AND FOCK MATRIX BLOCK ADDRESSES.    C
C   [3] DENSITIES: MOLECULAR DENSITIES, ENERGIES AND LEVEL SHIFTING.   C
C   [4] ATOMIC HARTREE-FOCK: AVERAGE OF CONFIG. ATOMIC SCF.            C
C   [5] MOLECULAR HARTREE-FOCK: MANY-CENTER HARTREE-FOCK SCF.          C
C   [6] MULTI-CONFIG: MANY-CETRE MULTICONFIG. SCF CALCULATIONS.        C
C   [7] DMRG: DENSITY MATRIX RENORMALISATION GROUP CALCULATIONS.       C
C   [8] MBPT: CORRELATION ENERGY CALCULATION ROUTINES.                 C
C   [9] EXPECTATION VALUES: OBSERVABLES FROM A CONVERGED SOLUTION.     C
C  [10] PLOTS: AMPLITUDES AND FIELDS/POTENTIALS IN DATA FILES.         C
C  [11] R-INTS: BOYS INTEGRALS R-INTEGRALS AND RELATED QUANTITIES.     C
C  [12] E-COEFFS: FINITE BASIS OVERLAP SPIN STRUCTURE FACTORS.         C
C**********************************************************************C
C
      CHARACTER*4  HMLTN
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
C
C     READ DATA FROM USER-SPECIFIED INPUT FILE
      CALL CARDIN
C
C     OPEN FILE FOR TERMINAL RECORD
      OPEN(UNIT=7,FILE=TRIM(OUTFL)//'.out',STATUS='UNKNOWN')
C
C     PRINT SUMMARY OF INPUT DATA
      CALL INPUT
C
C     PRINT MEMORY ALLOCATION SUMMARY
      CALL MEMORY
C
C     INTER-ATOMIC ANGLES AND NUCLEAR REPULSION ENERGY
      CALL NUCGEOM
C
C     FOCK MATRIX SYMMETRY TYPE INDICES
      CALL FOCKIND
C
C     CARTESIAN EXPANSION INDICES FOR BASIS FUNCTION OVERLAP PAIRS
      CALL CARTIND
C
C     ATOMIC HARTREE-FOCK SCF ROUTINE
      IF(INEW.EQ.0) THEN
        CALL ATOMIC
      ENDIF
C
C     MOLECULAR HARTREE-FOCK SCF ROUTINE
      IF(INEW.EQ.0.OR.ITREE.EQ.1) THEN
        IF(IMOL.NE.0) THEN
          CALL HFSCF
        ENDIF
      ENDIF
C
C     MANY-BODY DIAGRAMMATIC PERTURBATION THEORY
      IF(ITREE.EQ.2) THEN
        CALL MBPT
      ENDIF
C
C     MULTI-CONFIGURATIONAL SELF CONSISTENT FIELD CALCULATION
      IF(ITREE.EQ.3) THEN
        CALL MCSCF
      ENDIF
C
C     DENSITY MATRIX RENORMALISATION GROUP CALCULATION
      IF(ITREE.EQ.4) THEN
        CALL DMRG
      ENDIF
C
C     ONE-BODY HAMILTONIAN INTERACTIONS FROM A CONVERGED SOLUTION
      IF(ITREE.EQ.5) THEN
        CALL PT1BODY
      ENDIF
C
C     ELECTROMAGNETIC FIELDS AND POTENTIALS
      IF(ITREE.EQ.6) THEN
        CALL FIELDS
      ENDIF
C
C     PRINT SUMMARY OF OUTPUT DATA
      CALL OUTPUT
C
C     CLOSE FILE FOR TERMINAL RECORD
      CLOSE(UNIT=7)
C
C     SUCCESSFUL EXIT
      END PROGRAM
C
C
C**********************************************************************C
C ==================================================================== C
C   [1] INPUT/OUTPUT: READ FROM INPUT FILE AND SUMMARISE DATA.         C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] CARDIN: READ MOLECULAR DATA FROM A USER-SPECIFIED FILE.        C
C   [B] INPUT: WRITE A SUMMARY OF DATA INPUT OPTIONS TO TERMINAL.      C
C   [C] OUTPUT: WRITE A SUMMARY OF TOTAL CALCULATION STATS/DATA.       C
C**********************************************************************C
C
C
      SUBROUTINE CARDIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           CCCCCC     AA    RRRRRRR  DDDDDDD  IIII NN    NN           C
C          CC    CC   AAAA   RR    RR DD    DD  II  NNN   NN           C
C          CC        AA  AA  RR    RR DD    DD  II  NNNN  NN           C
C          CC       AA    AA RR    RR DD    DD  II  NN NN NN           C
C          CC       AAAAAAAA RRRRRRR  DD    DD  II  NN  NNNN           C
C          CC    CC AA    AA RR    RR DD    DD  II  NN   NNN           C
C           CCCCCC  AA    AA RR    RR DDDDDDD  IIII NN    NN           C
C                                                                      C
C -------------------------------------------------------------------- C
C  CARDIN READS AND PREPARES DATA FROM A USER-SPECIFIED INPUT FILE.    C
C  THIS IS ALSO WHERE ATOMIC ELEMENT NAMES AND CV ARE SPECIFIED.       C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*1  DUMLIN
      CHARACTER*2  ELMNT(120)
      CHARACTER*4  HMLTN
      CHARACTER*7  HMINT(10),PTYPE(10)
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 C(MDM,MDM)
C
      COMMON/ATOM/ELMNT
      COMMON/COEF/C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/FILL/NCNF(MCT,MKP,MKP+1),NLVL(MCT,MKP),IFILL(MCT)
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/LSHF/SHLEV(3),SHLV
      COMMON/PLOT/NPTYPE,PTYPE
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/PT1B/NHMINT,HMINT
      COMMON/SHLL/ACFF,BCFF,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN,NOELEC
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C
      DATA CV/1.370359898D2/
      DATA ELMNT/'H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne',
     &           'Na','Mg','Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca',
     &           'Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn',
     &           'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y' ,'Zr',
     &           'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     &           'Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd',
     &           'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     &           'Lu','Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg',
     &           'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     &           'Pa','U' ,'Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     &           'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds',
     &           'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og','Ue','Un'/
C
C**********************************************************************C
C     MOLECULE NAME AND CALCULATION TYPE                               C
C**********************************************************************C
C
C     MOLECULE LABEL
      READ(5,*) DUMLIN
      READ(5,*) MOLCL
C
C     CHOICE OF HAMILTONIAN HMLTN: NORL, BARE, DHFR, DHFP, DHFB OR DHFQ
10    FORMAT(A4)
      READ(5, *) DUMLIN
      READ(5,10) HMLTN
C
C     ALLOW VALID HAMILTONIANS TO PASS
      IF(HMLTN.EQ.'NORL'.OR.HMLTN.EQ.'BARE'.OR.HMLTN.EQ.'DHFR'.OR.
     &   HMLTN.EQ.'DHFP'.OR.HMLTN.EQ.'DHFB'.OR.HMLTN.EQ.'DHFQ') GOTO 20
C
C     UNKNOWN HAMILTONIAN - ABNORMAL EXIT
      WRITE(6, *) 'In CARDIN: unknown HMLTN value. HMLTN = ',HMLTN
      WRITE(7, *) 'In CARDIN: unknown HMLTN value. HMLTN = ',HMLTN
      STOP
C
20    CONTINUE
C
C     WAVE FUNCTION FILE NAME
      WFNFL = 'output/'//TRIM(MOLCL)//'_'//HMLTN//'.wfn'
C
C     HARTREE-FOCK SCF(1), MBPT(2), MCSCF(3), DMRG(4), HMINT(5), PLOTS(6)
      READ(5,*) DUMLIN
      READ(5,*) ITREE
C      
C     ENSURE USER HAS SELECTED VALID CHOICE OF ITREE
      IF(ITREE.LT.1.OR.ITREE.GT.6) THEN
        WRITE(6, *) 'In CARDIN: invalid calculation tree. ',ITREE
        WRITE(7, *) 'In CARDIN: invalid calculation tree. ',ITREE
        STOP
      ENDIF
C
C     OUTPUT FILE NAME
      IF(ITREE.EQ.1) THEN
        OUTFL = 'output/'//TRIM(MOLCL)//'_'//HMLTN//'_HFSCF'
      ELSEIF(ITREE.EQ.2) THEN
        OUTFL = 'output/'//TRIM(MOLCL)//'_'//HMLTN//'_MBPT'
      ELSEIF(ITREE.EQ.3) THEN
        OUTFL = 'output/'//TRIM(MOLCL)//'_'//HMLTN//'_MCSCF'
      ELSEIF(ITREE.EQ.4) THEN
        OUTFL = 'output/'//TRIM(MOLCL)//'_'//HMLTN//'_DMRG'
      ELSEIF(ITREE.EQ.5) THEN
        OUTFL = 'output/'//TRIM(MOLCL)//'_'//HMLTN//'_EXPVL'
      ELSEIF(ITREE.EQ.6) THEN
        OUTFL = 'output/'//TRIM(MOLCL)//'_'//HMLTN//'_PLOTS'
      ENDIF
C
C     NEW START(0), RESUME(1)
      READ(5,*) DUMLIN
      READ(5,*) INEW
C
C     READ IN A WAVE FUNCTION FILE IF PROMPTED
      IF(INEW.EQ.1) THEN
        OPEN(UNIT=8,FILE=WFNFL,STATUS='UNKNOWN')
        REWIND(UNIT=10)
        DO I=1,NDIM
          READ(8, *) EIGEN(I),(C(J,I),J=1,NDIM)
        ENDDO
        CLOSE(UNIT=8)
      ENDIF
C
C     E-COEFFICIENTS BY BATCH (0), TO LARGE EXTERNAL FILE (1)
      READ(5,*) DUMLIN
      READ(5,*) IEQS
C
C     SERIAL CALCULATION (0), OPENMP PARALLEL ENABLED (1)
      READ(5,*) DUMLIN
      READ(5,*) IPAR
C
C     MAX. NUMBER OF THREADS (UPDATE THIS TO TOTAL PROCESSORS MINUS ONE
      ICOR = 4
C
C**********************************************************************C
C     ATOMIC CENTERS AND BASIS FUNCTIONS                               C
C**********************************************************************C
C
C     BASIS SET TYPE: GEOMETRIC (1) OR OPTIMISED (2)
      READ(5,*) DUMLIN
      READ(5,*) INTYPE
C
C     NUMBER OF ATOMIC CENTERS
      READ(5,*) DUMLIN
      READ(5,*) NCNT
C
C     SPECIFY ATOM OR MOLECULE
      IF(NCNT.EQ.1) THEN
        IMOL = 0
      ELSE
        IMOL = 1
      ENDIF
C
C     INITIALISE MAXIMUM LQN AND DIMENSION COUNTERS
      LBIG = 0
      NDIM = 0
C
C     LOOP OVER NUCLEAR CENTERS
      DO ICNT=1,NCNT
C
C       CARTESIAN COORDINATES OF THIS CENTER       
        READ(5,*) DUMLIN
        READ(5,*) (COORD(J,ICNT),J=1,3)
C
C       ZNUC, ATOMIC MASS, MAXIMUM LQN AND ATOMIC CHARGE
        READ(5,*) DUMLIN
        READ(5,*) IZNUC(ICNT),AMASS(ICNT),LMAX(ICNT),IQNUC(ICNT)
C
C       AUFBAU FILLING FOR THIS CENTER: AUTOMATIC (0) OR MANUAL (1)
        READ(5,*) DUMLIN
        READ(5,*) IFILL(ICNT)
C
C       IF FILLING IS MANUAL, IMPORT ATOMIC ELECTRON CONFIGURATION
        IF(IFILL(ICNT).NE.0) THEN
          READ(5,*) DUMLIN
          DO L=1,LMAX(ICNT)+1
            READ(5,*) NLVL(ICNT,L),(NCNF(ICNT,L,N),N=1,NLVL(ICNT,L))
          ENDDO
        ENDIF
C
C       NUMBER OF KAPPA VALUES FOR THIS ATOM
        NKAP(ICNT) = 2*LMAX(ICNT)+1
C
C       NUCLEAR CHARGE AS A REAL VALUE
        ZNUC(ICNT) = DFLOAT(IZNUC(ICNT))
C
C       GAUSSIAN WIDTH PARAMETER FOR NUCLEAR CHARGE
        IF(IZNUC(ICNT).EQ.1) THEN
          CNUC(ICNT) = 2.1248239171D9
        ELSEIF(IZNUC(ICNT).EQ.8) THEN
          CNUC(ICNT) = 5.8631436655D8
        ELSE
          CDEN = AMASS(ICNT)**(1.0D0/3.0D0)
          CDEN = 0.836D0*CDEN + 0.57D0
          CDEN = 0.529177249D0/CDEN
          CNUC(ICNT) = 1.50D10*CDEN*CDEN
        ENDIF
C
C       UPDATE OVERALL MAXIMUM OCCURRING LQN
        IF(LMAX(ICNT).GT.LBIG) LBIG = LMAX(ICNT)
C
C ***   INITIATE IF STATEMENT FOR TYPE OF BASIS FUNCTION
        READ(5,*) DUMLIN
C
C >>>   GEOMETRIC BASIS FUNCTIONS
        IF(INTYPE.EQ.1) THEN
C
C         GENERATE THE EVEN TEMPERED ORBITAL EXPONENTS FOR EACH LQN
          DO LQN=0,LMAX(ICNT)
C
C           READ GENERATING PARAMETERS A, B AND NFUNCT
            READ(5,*) APARAM,BPARAM,NFUNCT(LQN+1,ICNT)
C
C           GENERATE NFUNCT BASIS EXPONENTS USING VARIABLE ZETA
            ZETA = APARAM
            DO IBAS=1,NFUNCT(LQN+1,ICNT)
              EXPSET(IBAS,LQN+1,ICNT) = ZETA
              ZETA = ZETA*BPARAM
            ENDDO
C
          ENDDO
C
C >>>   OPTIMISED EXPONENTS FROM A RECORDED LIST
        ELSEIF(INTYPE.EQ.2) THEN
C
C         READ IN THE OPTIMISED ORBITAL EXPONENTS FOR EACH LQN
          DO LQN=0,LMAX(ICNT)
C
C           READ NFUNCT
            READ(5,*) NFUNCT(LQN+1,ICNT)
C
C           READ BASIS EXPONENTS FROM A LIST
            DO IBAS=1,NFUNCT(LQN+1,ICNT)
              READ(5,*) EXPSET(IBAS,LQN+1,ICNT)
            ENDDO
C
          ENDDO
C
C ***   END IF STATEMENT FOR TYPE OF BASIS FUNCTION
        ENDIF
C
C       LOOP OVER ALL LQNS IN THIS CENTER AND ADD TO FOCK DIMENSION
        DO LQN=0,LMAX(ICNT)
C
C         EXTEND DIMENSION OF FOCK MATRIX
          NDIM = NDIM + 4*(2*LQN+1)*NFUNCT(LQN+1,ICNT)
C
C         ASSIGN KQN VALUES
          IF(LQN.NE.0) THEN
            KVALS(2*LQN  ,ICNT) = LQN
          ENDIF
          KVALS(2*LQN+1,ICNT) =-LQN-1
C
        ENDDO
C
C     END LOOP OVER NUCLEAR CENTERS        
      ENDDO
C
C     TOTAL DIMENSION DEPENDING ON CHOICE OF HAMILTONIAN
      IF(HMLTN.EQ.'NORL') THEN
        NDIM   = NDIM/2
        NSHIFT = 0
      ELSE
        NSHIFT = NDIM/2
      ENDIF
C
C**********************************************************************C
C     CLOSED/OPEN SHELL DETAILS                                        C
C**********************************************************************C
C
C     NUMBER OF CLOSED- AND OPEN-SHELL ELECTRONS AND TOTAL ELECTRONS
      READ(5,*) DUMLIN
      READ(5,*) NCLS,NOPN,NOELEC
C
C     TOTAL NUMBER OF ELECTRONS IN SYSTEM
C     DFNOTE: THIS IS A HACK FOR NOW
      NOCC = NCLS + NOPN
      NVIR = NDIM - NSHIFT - NOCC
C
C     NOCC = 0
C     DO IZ=1,NCNT
C       NOCC = NOCC + IZNUC(IZ) - IQNUC(IZ)
C     ENDDO
C
C *** INITIATE IF STATEMENT DEPENDING ON CLOSED/OPEN SHELLS

C >>> OPEN-SHELL MOLECULE
      goto 555
C     DFNOTE: ENABLE THIS OPTION
      IF(NOPN.NE.0) THEN
C
C       FRACTIONAL OCCUPANCY OF THE OPEN SHELL
        FOPEN = DFLOAT(NOELEC)/DFLOAT(NOPN)
C
C       LABELS FOR THE OPEN SHELL
        READ(5,*) DUMLIN
        READ(5,*) ACFF,BCFF,(IOPN(M),M=1,NOPN)
C
C       PRINT THE LABELS FOR THE CLOSED-SHELL SPINORS USING KNOWN
C       IDENTITY OF THE OPEN-SHELL SPINORS
        JCL = 1
        JOP = 1
        DO JCOUNT=1,NOCC
C
C         APPLY LABEL TO EACH ORBITAL
          IF(JCOUNT.NE.IOPN(JOP)) THEN
            ICLS(JCL) = JCOUNT
            JCL = JCL + 1
          ELSE
            JOP = JOP + 1
          ENDIF
C
        ENDDO
C
C >>> CLOSED-SHELL MOLECULE
      ELSE
C
C       LABEL THE CLOSED-SHELL ELECTRONS
        DO JCL=1,NCLS
          ICLS(JCL) = JCL
        ENDDO
C
C *** END IF STATEMENT FOR CLOSED/OPEN SHELLS
      ENDIF
C
555   continue
C
C     ALL SPINORS ARE CLOSED
      DO JCL=1,NCLS
        ICLS(JCL) = JCL
      ENDDO
C
C**********************************************************************C
C     LEVEL SHIFTING AND INTEGRAL INCLUSION                            C
C**********************************************************************C
C
C     LEVEL SHIFT PARAMETER FOR EACH INTEGRAL STAGE (SKAANE 4.4.3)
      READ(5,*) DUMLIN
      READ(5,*) (SHLEV(N),N=1,3)
C
C     STARTING STAGE OF INTEGRAL INCLUSION LEVEL (1-3)
      READ(5,*) DUMLIN
      READ(5,*) ILEV
C
C     REASONS TO CHANGE THE INTEGRAL INCLUSION LEVEL
      IF(HMLTN.EQ.'NORL') THEN
        ILEV = 1
      ENDIF
C
      IF(NCNT.EQ.1.OR.INEW.EQ.1) THEN
        ILEV = 3
      ENDIF
C
C     IMPLEMENT THE STARTING SHIFT FACTOR
      IF(ILEV.GE.1.AND.ILEV.LE.3) THEN
        SHLV = SHLEV(ILEV)
      ELSE
        WRITE(6, *) 'In CARDIN: invalid starting stage. ILEV = ',ILEV
        WRITE(7, *) 'In CARDIN: invalid starting stage. ILEV = ',ILEV
        STOP
      ENDIF
C
C**********************************************************************C
C     NON-HF CALCULATION DETAILS                                       C
C**********************************************************************C
C
C     EXPECTATION VALUE CALCULATIONS: ELCMNPL,MAGDIPL ETC
      IF(ITREE.EQ.5) THEN
C
C       READ NUMBER OF EXPECTATION VALUES
        READ(5,*) DUMLIN
        READ(5,*) NHMINT
C
C       PLACE LIMIT ON POSSIBLE NUMBER OF THEM
        IF(NHMINT.LT.1.OR.NHMINT.GT.10) THEN
          WRITE(6, *) 'In CARDIN: invalid number of expectation values.'
          WRITE(7, *) 'In CARDIN: invalid number of expectation values.'
          STOP
        ENDIF
C
C       READ IN EACH INTERACTION HAMILTONIAN
        DO N=1,NHMINT
          READ(5,*) HMINT(N)
        ENDDO
C
      ENDIF
C
C     DATA PLOTTING: AMPLITUDES, EM FIELDS AND POTENTIALS
      IF(ITREE.EQ.6) THEN
C
C       READ NUMBER OF PLOT TYPES
        READ(5,*) DUMLIN
        READ(5,*) NPTYPE
C
C       PLACE LIMIT ON POSSIBLE NUMBER OF THEM
        IF(NPTYPE.LT.1.OR.NPTYPE.GT.10) THEN
          WRITE(6, *) 'In CARDIN: invalid number of plot types.'
          WRITE(7, *) 'In CARDIN: invalid number of plot types.'
          STOP
        ENDIF
C
C       READ IN EACH PLOT TYPE
        DO N=1,NPTYPE
          READ(5,*) PTYPE(N)
        ENDDO
C
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE INPUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C               IIII NN    NN PPPPPPP  UU    UU TTTTTTTT               C
C                II  NNN   NN PP    PP UU    UU    TT                  C
C                II  NNNN  NN PP    PP UU    UU    TT                  C
C                II  NN NN NN PP    PP UU    UU    TT                  C
C                II  NN  NNNN PPPPPPP  UU    UU    TT                  C
C                II  NN   NNN PP       UU    UU    TT                  C
C               IIII NN    NN PP        UUUUUU     TT                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  INPUT PRINTS MOLECULAR DATA INPUT OPTIONS TO THE TERMINAL.          C
C**********************************************************************C
      PARAMETER(MDM=1200,MCT=6,MBS=26,MKP=9)
C
      CHARACTER*4  HMLTN
      CHARACTER*20 STAMP
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
C
C     PRINT A BIG SECTION HEADER
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) REPEAT(' ',30),'INPUT SUMMARY'
      WRITE(7, *) REPEAT(' ',30),'INPUT SUMMARY'
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C     TITLE FOR INPUT OPTIONS
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',21),'Calculation tree and data files'
      WRITE(7, *) REPEAT(' ',21),'Calculation tree and data files'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     CALCULATION TREE
      IF(ITREE.EQ.1) THEN
        WRITE(6, *) 'Calculation tree:',REPEAT(' ',43),'Hartree-Fock'
        WRITE(7, *) 'Calculation tree:',REPEAT(' ',43),'Hartree-Fock'
      ELSEIF(ITREE.EQ.2) THEN
        WRITE(6, *) 'Calculation tree:',REPEAT(' ',51),'MBPT'
        WRITE(7, *) 'Calculation tree:',REPEAT(' ',51),'MBPT'
      ELSEIF(ITREE.EQ.3) THEN
        WRITE(6, *) 'Calculation tree:',REPEAT(' ',50),'MCSCF'
        WRITE(7, *) 'Calculation tree:',REPEAT(' ',50),'MCSCF'
      ELSEIF(ITREE.EQ.4) THEN
        WRITE(6, *) 'Calculation tree:',REPEAT(' ',51),'DMRG'
        WRITE(7, *) 'Calculation tree:',REPEAT(' ',51),'DMRG'
      ELSEIF(ITREE.EQ.5) THEN
        WRITE(6, *) 'Calculation tree:',REPEAT(' ',44),'Expct. vals'
        WRITE(7, *) 'Calculation tree:',REPEAT(' ',44),'Expct. vals'
      ELSEIF(ITREE.EQ.6) THEN
        WRITE(6, *) 'Calculation tree:',REPEAT(' ',47),'Plotting'
        WRITE(7, *) 'Calculation tree:',REPEAT(' ',47),'Plotting'
      ENDIF
C
C     PRINT THE HAMILTONIAN OPTION
      WRITE(6, *) 'Hamiltonian type:',REPEAT(' ',51),HMLTN
      WRITE(7, *) 'Hamiltonian type:',REPEAT(' ',51),HMLTN
C
C     CONFIRM SOLUTION SPACE DIMENSION OR EXIT
      IF(NDIM.LE.MDM) THEN
        WRITE(6, *) 'Fock matrix dimension: ',REPEAT(' ',37),NDIM
        WRITE(7, *) 'Fock matrix dimension: ',REPEAT(' ',37),NDIM
      ELSEIF(NDIM.GT.MDM) THEN
        WRITE(6, *) 'In INPUT: matrix dimension too big. NDIM = ',NDIM
        WRITE(7, *) 'In INPUT: matrix dimension too big. NDIM = ',NDIM
        STOP
      ENDIF
C
C     NEW START OR READ IN HFSCF EXPANSION COEFFICIENTS
      IF(INEW.EQ.0) THEN
        WRITE(6, *) 'SCF density matrix:',REPEAT(' ',44),'New start'
        WRITE(7, *) 'SCF density matrix:',REPEAT(' ',44),'New start'
      ELSE
        WRITE(6, *) 'SCF density matrix:',REPEAT(' ',46),'Read in'
        WRITE(7, *) 'SCF density matrix:',REPEAT(' ',46),'Read in'
      ENDIF
C
C     E-COEFFICIENT CALCULATION
      IF(IEQS.EQ.0) THEN
        WRITE(6, *) 'E-coefficients:',REPEAT(' ',49),'By batch'
        WRITE(7, *) 'E-coefficients:',REPEAT(' ',49),'By batch'
      ELSE
        WRITE(6, *) 'E-coefficients:',REPEAT(' ',45),'Save to file'
        WRITE(7, *) 'E-coefficients:',REPEAT(' ',45),'Save to file'
      ENDIF
C
C     OPENMP PARALLEL OPTION
      IF(IPAR.EQ.0) THEN
        WRITE(6, *) 'OpenMP parallel option:',REPEAT(' ',41),'Disabled'
        WRITE(7, *) 'OpenMP parallel option:',REPEAT(' ',41),'Disabled'
      ELSE
        WRITE(6, *) 'OpenMP parallel option:',REPEAT(' ',42),'Enabled'
        WRITE(7, *) 'OpenMP parallel option:',REPEAT(' ',42),'Enabled'
        WRITE(6, *) 'OpenMP thread limit:'   ,REPEAT(' ',40),ICOR
        WRITE(7, *) 'OpenMP thread limit:'   ,REPEAT(' ',40),ICOR
      ENDIF
C
C     SECTION FOR FILE NAMES
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     PRINT INPUT FILE NAME
      LF = LEN(TRIM(MOLCL))
      WRITE(6, *) 'Molecule name:',REPEAT(' ',58-LF),TRIM(MOLCL)
      WRITE(7, *) 'Molecule name:',REPEAT(' ',58-LF),TRIM(MOLCL)
C
C     PRINT FILE OUTPUT NAMES
      LN = LEN(TRIM(OUTFL))
      WRITE(6, *) 'Output file name: ',REPEAT(' ',54-LN),TRIM(OUTFL)
      WRITE(7, *) 'Output file name: ',REPEAT(' ',54-LN),TRIM(OUTFL)
C
C     RECORD TIME AT BEGINNING OF CALCULATION
      CALL CPU_TIME(TBEG)
      CALL TIMENOW(STAMP)
      WRITE(6, *) 'Time at BERTHA initiation:',REPEAT(' ',26),STAMP
      WRITE(7, *) 'Time at BERTHA initiation:',REPEAT(' ',26),STAMP
C
C     END OF INPUT SUMMARY
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
      RETURN
      END
C
C
      SUBROUTINE MEMORY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     MM       MM EEEEEEEE MM       MM  OOOOOO  RRRRRRR  YY    YY      C
C     MMM     MMM EE       MMM     MMM OO    OO RR    RR YY    YY      C
C     MMMM   MMMM EE       MMMM   MMMM OO    OO RR    RR  YY  YY       C
C     MM MM MM MM EEEEEE   MM MM MM MM OO    OO RR    RR   YYYY        C
C     MM  MMM  MM EE       MM  MMM  MM OO    OO RRRRRRR     YY         C
C     MM   M   MM EE       MM   M   MM OO    OO RR    RR    YY         C
C     MM       MM EEEEEEEE MM       MM  OOOOOO  RR    RR    YY         C
C                                                                      C
C -------------------------------------------------------------------- C
C  MEMORY SUMMARISES THE SIZE AND MEMORY REQUIREMENTS OF BIG ARRAYS.   C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MFL=10000000,MKP=9,
     &                           MNU=MKP+1,MAB=2*MNU+5,NINT=1801,MLM=30)
C
      INTEGER*8 NAMEM,NCMEM,NDMEM,NEMEM,NMMEM,NRMEM,NTMEM
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 C(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)

      COMMON/BSKL/EK(MB2),EKL0(MB2),EIK(MB2,MAB),EJL(MB2,MAB),
     &            EL(MB2),EKL1(MB2),EKL(MB2,MAB),RNKL(MB2,4),
     &            B1(MB2,MAB,MAB),B2(MB2,MAB,MAB),B3(MB2,MAB,MAB),
     &            B4(MB2,MAB,MAB),B5(MB2,MAB,MAB),B6(MB2,MAB,MAB)
      COMMON/BTRE/BKLLSS(MB2,4),BKSLSL(MB2,4),BKSSLL(MB2,4),
     &            BMSLSL(MB2,4)
      COMMON/CLRE/CKLLLL(MB2,4),CKSSSS(MB2,4),CKSLSL(MB2,4),
     &            CJLLLL(MB2,4),CJSSSS(MB2,4),CJLLSS(MB2,4),
     &            CJSSLL(MB2,4)
      COMMON/COEF/C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EILS/EILSFL(MFL,24),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
C
C
C     TITLE FOR MEMORY SUMMARY
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',20),'Approx. system memory allocation'
      WRITE(7, *) REPEAT(' ',20),'Approx. system memory allocation'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)

C
C     APPROXIMATE THE MEMORY STORAGE ALLOCATION
      NAMEM = SIZE(B1) + SIZE(B2) + SIZE(B3) 
     &      + SIZE(B4) + SIZE(B5) + SIZE(B6)
     &      + SIZE(EK) + SIZE(EKL0) + SIZE(EIK) + SIZE(EJL) 
     &      + SIZE(EL) + SIZE(EKL1) + SIZE(EKL) + SIZE(RNKL)
      NCMEM = 2*SIZE(C)
      NDMEM = 2*SIZE(DENC) + 2*SIZE(DENO) + 2*SIZE(DENT)
      NEMEM = SIZE(E0LLFL) + SIZE(E0SSFL) + SIZE(EILSFL)
      NMMEM = 2*SIZE(OVAP) + 2*SIZE(HNUC) + 2*SIZE(HKIN)
     &      + 2*SIZE(VUEH) + 2*SIZE(GDIR) + 2*SIZE(GXCH)
     &      + 2*SIZE(QDIR) + 2*SIZE(QXCH) + 2*SIZE(BDIR)
     &      + 2*SIZE(BXCH) + 2*SIZE(FOCK)
      NRMEM = SIZE(CKLLLL) + SIZE(CKSSSS) + SIZE(CKSLSL) + SIZE(CJLLLL) 
     &      + SIZE(CJSSSS) + SIZE(CJLLSS) + SIZE(CJSSLL) 
     &      + SIZE(BKLLSS) + SIZE(BKSLSL) + SIZE(BKSSLL) + SIZE(BMSLSL)
      NTMEM = NAMEM + NCMEM + NDMEM + NEMEM + NMMEM + NRMEM
C
C     SIZES (IN GIGABYTES)
      SAMEM = NAMEM*8.0D-9
      SCMEM = NCMEM*8.0D-9
      SDMEM = NDMEM*8.0D-9
      SEMEM = NEMEM*8.0D-9
      SMMEM = NMMEM*8.0D-9
      SRMEM = NRMEM*8.0D-9
      STMEM = NTMEM*8.0D-9
C
20    FORMAT(1X,A,4X,A,19X,A,9X,A)
21    FORMAT(1X,A,6X,A,5X,I12,9X,F9.5)
      WRITE(6,20) 'COMMON','Description','Words (double)','Size (GB)'
      WRITE(7,20) 'COMMON','Description','Words (double)','Size (GB)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,21) 'COEF','Expansion coeffs.          ',NCMEM,SCMEM
      WRITE(7,21) 'COEF','Expansion coeffs.          ',NCMEM,SCMEM
      WRITE(6,21) 'DENS','Density matrices           ',NDMEM,SDMEM
      WRITE(7,21) 'DENS','Density matrices           ',NDMEM,SDMEM
      WRITE(6,21) 'CLRE','Atomic Coulomb/Breit ints. ',NRMEM,SRMEM
      WRITE(7,21) 'CLRE','Atomic Coulomb/Breit ints. ',NRMEM,SRMEM
      WRITE(6,21) 'BSKL','Atomic mean-field submatrx.',NAMEM,SAMEM
      WRITE(7,21) 'BSKL','Atomic mean-field submatrx.',NAMEM,SAMEM
      WRITE(6,21) 'MTRX','Molecular Fock matrices    ',NMMEM,SMMEM
      WRITE(7,21) 'MTRX','Molecular Fock matrices    ',NMMEM,SMMEM
      WRITE(6,21) 'EQTT','Eq-coefficient file        ',NEMEM,SEMEM
      WRITE(7,21) 'EQTT','Eq-coefficient file        ',NEMEM,SEMEM
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,21) '    ','Total                      ',NTMEM,STMEM
      WRITE(7,21) '    ','Total                      ',NTMEM,STMEM
C
C     END OF MEMORY SUMMARY
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
      RETURN
      END
C
C
      SUBROUTINE OUTPUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          OOOOOO  UU    UU TTTTTTTT PPPPPPP  UU    UU TTTTTTTT        C
C         OO    OO UU    UU    TT    PP    PP UU    UU    TT           C
C         OO    OO UU    UU    TT    PP    PP UU    UU    TT           C
C         OO    OO UU    UU    TT    PP    PP UU    UU    TT           C
C         OO    OO UU    UU    TT    PPPPPPP  UU    UU    TT           C
C         OO    OO UU    UU    TT    PP       UU    UU    TT           C
C          OOOOOO   UUUUUU     TT    PP        UUUUUU     TT           C
C                                                                      C
C -------------------------------------------------------------------- C
C  OUTPUT PRINTS A SUMMARY OF TOTAL CALCULATION DATA TO TERMINAL.      C
C**********************************************************************C
      PARAMETER(MCT=6,MBS=26,MKP=9)
C
      CHARACTER*4  HMLTN
      CHARACTER*16 HMS
      CHARACTER*20 STAMP
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
C
C     PRINT A BIG SECTION HEADER
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) REPEAT(' ',30),'OUTPUT SUMMARY'
      WRITE(7, *) REPEAT(' ',30),'OUTPUT SUMMARY'
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C     RECORD TIME AT BERTHA EXIT
      CALL TIMENOW(STAMP)
      CALL CPU_TIME(TEND)

      TTOT = TEND-TBEG
C
C     PRINT TABLE OF DATA
20    FORMAT(1X,A,26X,A)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',28),'CPU time summary'
      WRITE(7, *) REPEAT(' ',28),'CPU time summary'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      IF(INEW.EQ.0) THEN
        WRITE(6,20) 'Atomic Hartree-Fock SCF:      ',HMS(TATM)
        WRITE(7,20) 'Atomic Hartree-Fock SCF:      ',HMS(TATM)
      ENDIF
      IF(ITREE.EQ.1.OR.INEW.EQ.0) THEN
        WRITE(6,20) 'Molecular Hartree-Fock SCF:   ',HMS(TSCF)
        WRITE(7,20) 'Molecular Hartree-Fock SCF:   ',HMS(TSCF)
      ELSEIF(ITREE.EQ.2) THEN
        WRITE(6,20) 'Many-body perturbation theory:',HMS(TMPT)
        WRITE(7,20) 'Many-body perturbation theory:',HMS(TMPT)
      ELSEIF(ITREE.EQ.3) THEN
        WRITE(6,20) 'Multi-configurational SCF:    ',HMS(TMCF)
        WRITE(7,20) 'Multi-configurational SCF:    ',HMS(TMCF)
      ELSEIF(ITREE.EQ.4) THEN
        WRITE(6,20) 'Density matrix renorm. group: ',HMS(TDMG)
        WRITE(7,20) 'Density matrix renorm. group: ',HMS(TDMG)
      ELSEIF(ITREE.EQ.5) THEN
        WRITE(6,20) 'Property calculation:         ',HMS(TPRP)
        WRITE(7,20) 'Property calculation:         ',HMS(TPRP)
      ELSEIF(ITREE.EQ.6) THEN
        WRITE(6,20) 'Plotting:                     ',HMS(TPLT)
        WRITE(7,20) 'Plotting:                     ',HMS(TPLT)
      ENDIF
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) 'Total CPU time:               ',HMS(TTOT)
      WRITE(7,20) 'Total CPU time:               ',HMS(TTOT)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     SUCCESSFUL EXIT MESSAGE
      WRITE(6, *)
      WRITE(7, *)
      WRITE(6, *) 'Successful BERTHA exit at:',REPEAT(' ',26),STAMP
      WRITE(7, *) 'Successful BERTHA exit at:',REPEAT(' ',26),STAMP
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C   [2] LABELS: MOLECULAR GEOMETRY AND FOCK MATRIX BLOCK ADDRESSES.    C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] NUCGEOM: BOND DISTANCES AND NUCLEAR REPULSION ENERGY.          C
C   [B] FOCKIND: CALCULATE ADDRESSES OF FOCK MATRIX FOR BASIS QN'S.    C
C   [C] CARTIND: GENERATES INDICES FOR EQ-COEFFS AND R-INTEGRALS.      C
C   [D] AUFBAU: DETERMINES GROUND STATE ATOMIC ELECTRON CONFIG.        C
C   [E] SPECTRM: DISPLAYS EIGENVALES AND ATOMIC TERM SYMBOLS.          C
C   [F] ELLTERM: GIVES THE CHARACTER CORRESPONDING TO LQN VALUE L.     C
C   [G] ROTATE: PERFORM TWO EULER ROTATIONS ON ALL ATOMIC CENTERS.     C
C   [H] MMPROD: PRODUCT OF TWO SQUARE ARRAYS OF DOUBLES.               C
C   [I] MVPROD: PRODUCT OF A SQUARE MATRIX AND VECTOR OF DOUBLES.      C
C   [J] HMS: RETURNS A QUOTED TIME IN SECONDS AS 'HR-MIN-SEC'.         C
C   [K] TIMENOW: RETURNS A DATE STRING FOR THE CPU TIME.               C
C**********************************************************************C
C
C
      SUBROUTINE NUCGEOM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  NN    NN UU    UU  CCCCCC   GGGGGG  EEEEEEEE  OOOOOO  MM       MM   C
C  NNN   NN UU    UU CC    CC GG    GG EE       OO    OO MMM     MMM   C
C  NNNN  NN UU    UU CC       GG       EE       OO    OO MMMM   MMMM   C
C  NN NN NN UU    UU CC       GG       EEEEEE   OO    OO MM MM MM MM   C
C  NN  NNNN UU    UU CC       GG   GGG EE       OO    OO MM  MMM  MM   C
C  NN   NNN UU    UU CC    CC GG    GG EE       OO    OO MM   M   MM   C
C  NN    NN  UUUUUU   CCCCCC   GGGGGG  EEEEEEEE  OOOOOO  MM       MM   C
C                                                                      C
C -------------------------------------------------------------------- C
C  NUCGEOM TRANSLATES AND ROTATES A MOLECULE IN A WAY THAT IS BEST     C
C  SUITED TO EFFICIENT COMPUTATION, IDENTIFIES MOLECULAR SHAPE AND     C
C  BOND DISTANCES, AND CALCULATES NUCLEAR REPULSION ENERGY.            C
C -------------------------------------------------------------------- C
C  DFNOTE: THE TSYM PACKAGE (WERNER 1996) CAN BE DIRECTLY IMPLEMENTED  C
C          AT THIS POINT AND PERFORMS POINT GROUP SYMMETRY ANALYSIS.   C
C          THIS RESULTS IN A HAMILTONIAN MATRIX OF BLOCK STRUCTURE     C
C          ACCORDING TO IRREPS, AND SPEEDS UP SCF CALCULATIONS.        C
C**********************************************************************C
      PARAMETER(MBS=26,MCT=6,MKP=9)
C
      CHARACTER*2 ELMNT(120),ELA,ELB,ELC
      CHARACTER*8 SHAPE,SPCES
C
      DIMENSION XYZ(3,MCT),DIST(MCT),CENT(3),IZAD(MCT,3)
C
      COMMON/ATOM/ELMNT
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EMDR,EMXC,EUEH
      COMMON/GEOM/SHAPE
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR

      DATA PI/3.1415926535897932D0/
C
C
C     PRINT A BIG SECTION HEADER
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) REPEAT(' ',27),'MOLECULAR COORDINATES'
      WRITE(7, *) REPEAT(' ',27),'MOLECULAR COORDINATES'
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C**********************************************************************C
C     IDENTIFY LARGEST ATOMIC CENTERS                                  C
C**********************************************************************C
C
C     IDENTIFY THE THREE LARGEST NUCLEAR CHARGES
      LZ1 = 0
      LZ2 = 0
      LZ3 = 0
      DO N=1,NCNT
        IZ = IZNUC(N)
        IF(IZ.GT.LZ1) THEN
          LZ3 = LZ2
          LZ2 = LZ1
          LZ1 = IZ
        ELSEIF(IZ.GT.LZ2.AND.IZ.LT.LZ1) THEN
          LZ3 = LZ2
          LZ2 = IZ
        ELSEIF(IZ.GT.LZ3.AND.IZ.LT.LZ2) THEN
          LZ3 = IZ
        ENDIF
      ENDDO
C
C     COUNT CENTERS WITH TOP THREE NUCLEAR CHARGES
      NZ1 = 0
      NZ2 = 0
      NZ3 = 0
      DO N=1,NCNT
        IZ = IZNUC(N)
        IF(IZ.EQ.LZ1) THEN
          NZ1 = NZ1 + 1
          IZAD(NZ1,1) = N
        ELSEIF(IZ.EQ.LZ2) THEN
          NZ2 = NZ2 + 1
          IZAD(NZ2,2) = N
        ELSEIF(IZ.EQ.LZ3) THEN        
          NZ3 = NZ3 + 1
          IZAD(NZ3,3) = N
        ENDIF
      ENDDO
C
C     DESIGNATE SPECIES OF MOLECULE
      IF(NZ2.EQ.0) THEN
        SPCES = 'HOMONUC.'
      ELSE
        SPCES = 'HTRONUC.'
      ENDIF
C
C**********************************************************************C
C     TRANSLATION: ORIGIN COINCIDES WITH CENTER OF ALL HEAVY ATOMS     C
C**********************************************************************C
C
C     CALCULATE CENTER OF ALL CHARGES LZ1
      DCNT = 0.0D0
      DO I=1,3
        CENT(I) = 0.0D0
        DO N=1,NZ1
          CENT(I) = CENT(I) + COORD(I,IZAD(N,1))
        ENDDO
        CENT(I) = CENT(I)/NZ1
        DCNT    = DCNT + CENT(I)**2
      ENDDO
      DCNT = DSQRT(DCNT)
C
C     TRANSLATE MOLECULE BY CENT(I) AND CALCULATE DISTANCE TO ORIGIN
      DO N=1,NCNT
        DO I=1,3
          XYZ(I,N) = COORD(I,N) - CENT(I)
        ENDDO
        DIST(N) = DSQRT(XYZ(1,N)**2 + XYZ(2,N)**2 + XYZ(3,N)**2)
      ENDDO
C
C     ALERT THE USER IF TRANSLATION WAS SUBSTANTIAL
      IF(DCNT.GT.1.0D-4) THEN
        WRITE(6, *) 'In NUCGEOM: molecule has been translated.'
        WRITE(7, *) 'In NUCGEOM: molecule has been translated.'
      ENDIF
C
C     ONE-CENTER (MUST BE ATOMIC)
      IF(NCNT.EQ.1) THEN
C
C       DESIGNATE SHAPE
        SHAPE = 'ATOMIC'
C
C       NO ROTATIONS NECESSARY
        GOTO 300
C
      ENDIF
C
C**********************************************************************C
C     FIRST ROTATION: FIX ON A ROTATION CENTER AND ROTATE TO Z-AXIS.   C
C**********************************************************************C
C
C     IF AN LZ1 CENTER IS ON Z-AXIS BUT NOT THE ORIGIN, SKIP ROTATION
      DO NZ=1,NZ1
        N = IZAD(NZ,1)
        X = DIST(N) - DABS(XYZ(3,N))
        IF(DABS(X).LT.1.0D-6.AND.DIST(N).GT.1.0D-4) THEN
          GOTO 200
        ENDIF
      ENDDO
C
C     FIX ON AN LZ1 CENTER THAT DOES NOT LIE ON ORIGIN
      DO NZ=1,NZ1
        N = IZAD(NZ,1)
        IF(DIST(N).GT.1.0D-4) THEN
          DO I=1,3
            CENT(I) = XYZ(I,N)
          ENDDO
          DCNT = DIST(N)
          GOTO 35
        ENDIF
      ENDDO
C
C     IF WE HAVE MADE IT THIS FAR, THERE IS NO LZ1 CENTER ALREADY ON
C     THE Z-AXIS AND NO LZ1 CENTER OFF THE ORIGIN --> MOVE TO LZ2.
C
C     VECTOR CENTER OF ALL LZ2
      DCNT = 0.0D0
      DO I=1,3
        CENT(I) = 0.0D0
        DO N=1,NZ2
          CENT(I) = CENT(I) + XYZ(I,IZAD(N,2))
        ENDDO
        CENT(I) = CENT(I)/NZ2
        DCNT    = DCNT + CENT(I)**2
      ENDDO
      DCNT = DSQRT(DCNT)
C
C     IF CENTER OF LZ2 IS ON Z-AXIS BUT NOT ORIGIN, SKIP ROTATION
      X = DCNT - DABS(CENT(3))
      IF(DABS(X).LT.1.0D-4.AND.DCNT.GT.1.0D-4) THEN
        GOTO 200
      ENDIF
C
C     IF CENTER OF LZ2 IS ON ORIGIN, FIX ON AN LZ2 NOT ON THE ORIGIN
      IF(DCNT.LT.1.0D-4) THEN
        DO NZ=1,NZ2
          N = IZAD(NZ,2)
          IF(DIST(N).GT.1.0D-4) THEN
            DO I=1,3
              CENT(I) = XYZ(I,N)
            ENDDO
            DCNT = DIST(N)
            GOTO 35
          ENDIF
        ENDDO
      ENDIF
C
C     LZ1 AND LZ2 ALWAYS SUFFICIENT TO DEFINE A Z ORIENTATION
C
C     SKIP POINT (ROTATION CENTER HAS BEEN IDENTIFIED)
35    CONTINUE
C
C     CALCULATE FIRST ROTATION ANGLE
C     ALPHA IS ROTATION ANGLE AROUND XY TO BRING CENT TO Y=0
C     BETA  IS ROTATION ANGLE FROM (X,Z) TO (0,Z')
      RLT = DSQRT(CENT(1)**2 + CENT(2)**2)
      ALPHA = DACOS(CENT(1)/RLT)
      BETA  = DACOS(CENT(3)/DCNT)
      CALL ROTATE(XYZ,NCNT,ALPHA,BETA)
      WRITE(6, *) 'In NUCGEOM: molecule has been rotated.'
      WRITE(7, *) 'In NUCGEOM: molecule has been rotated.'
C
C     NO Z-AXIS ROTATION NECESSARY
200   CONTINUE
C
C     TWO-CENTER (MUST BE DIATOMIC)
      IF(NCNT.EQ.2) THEN
C
C       DESIGNATE SHAPE
        SHAPE = 'DIATOMIC'
C
C       NO FURTHER ROTATIONS NECESSARY
        GOTO 300
C        
      ENDIF
C
C**********************************************************************C
C     SECOND ROTATION: FIX ON A ROTATION CENTER AND ROTATE TO Y-AXIS.  C
C**********************************************************************C
C
C     CHECK WHETHER MOLECULE IS LINEAR
      DO N=1,NCNT
        TEST = XYZ(1,N)**2 + XYZ(2,N)**2
        IF(TEST.GT.1.0D-6) GOTO 40
      ENDDO
      SHAPE = 'LINEAR'
      GOTO 300
40    CONTINUE
C
C     IF AN LZ1 CENTER IS ON YZ-AXIS BUT NOT THE ORIGIN, SKIP ROTATION
      DO NZ=1,NZ1
        N = IZAD(NZ,1)
        X = DIST(N) - DSQRT(XYZ(2,N)**2 + XYZ(3,N)**2)
        IF(DABS(X).LT.1.0D-6.AND.DIST(N).GT.1.0D-4) THEN
          GOTO 300
        ENDIF
      ENDDO
C
C     FIX ON AN LZ1 CENTER THAT DOES NOT LIE ON ORIGIN OR Z-AXIS
      DO NZ=1,NZ1
        N = IZAD(NZ,1)
        X = DIST(N) - DABS(XYZ(2,N))
        IF(DABS(X).GT.1.0D-4.AND.DIST(N).GT.1.0D-4) THEN
          DO I=1,3
            CENT(I) = XYZ(I,N)
          ENDDO
          DCNT = DIST(N)
          GOTO 50
        ENDIF
      ENDDO
C
C     IF WE HAVE MADE IT THIS FAR, THERE IS NO LZ1 CENTER ALREADY ON
C     THE YZ-AXIS AND NO LZ1 CENTER OFF THE ORIGIN --> MOVE TO LZ2.
C
C     IF THERE ARE NO LZ2 CENTERS, SKIP ROTATION
      IF(NZ2.EQ.0) THEN
        GOTO 300
      ENDIF
C
C     VECTOR CENTER OF ALL LZ2
      DCNT = 0.0D0
      DO I=1,3
        CENT(I) = 0.0D0
        DO N=1,NZ2
          CENT(I) = CENT(I) + XYZ(I,IZAD(N,2))
        ENDDO
        CENT(I) = CENT(I)/NZ2
        DCNT    = DCNT + CENT(I)**2
      ENDDO
      DCNT = DSQRT(DCNT)
C
C     IF CENTER OF LZ2 IS ON YZ-AXIS BUT NOT ORIGIN, SKIP ROTATION
      X = DCNT - DSQRT(CENT(2)**2 + CENT(3)**2)
      IF(DABS(X).LT.1.0D-3.AND.DCNT.GT.1.0D-3) THEN
        GOTO 300
      ENDIF
C
C     IF CENTER OF LZ2 IS NOT ON Y-AXIS, MAKE IT THE CENTER
      X = DCNT - DABS(CENT(2))
      IF(DABS(X).GT.1.0D-3.AND.DCNT.GT.1.0D-3) THEN
        GOTO 50
      ENDIF
C
C     IF CENTER OF LZ2 IS ON Y-AXIS, FIX ON ONE OF THEM
      X = DCNT - DABS(CENT(2))
      IF(DABS(X).LT.1.0D-3.AND.DCNT.GT.1.0D-3) THEN
        DO NZ=1,NZ2
          N = IZAD(NZ,2)
          X = DIST(N) - DABS(XYZ(2,N))
          IF(DABS(X).GT.1.0D-3.AND.DIST(N).GT.1.0D-3) THEN
            DO I=1,3
              CENT(I) = XYZ(I,N)
            ENDDO
            DCNT = DIST(N)
            GOTO 50
          ENDIF
        ENDDO
      ENDIF
C
C     IF CENTER OF LZ2 IS ON ORIGIN, FIX ON ONE THAT IS OFF THE Z AXIS
      IF(DCNT.LT.1.0D-4) THEN
        DO NZ=1,NZ2
          N = IZAD(NZ,2)
          X = DIST(N) - DABS(CENT(3))
          IF(DABS(X).GT.1.0D-3.AND.DIST(N).GT.1.0D-3) THEN
            DO I=1,3
              CENT(I) = XYZ(I,N)
            ENDDO
            DCNT = DIST(N)
            GOTO 50
          ENDIF
        ENDDO
      ENDIF
C
C     IF WE HAVE MADE IT THIS FAR, THERE IS NO LZ2 CENTER ALREADY ON
C     THE YZ-AXIS AND NO LZ1 CENTER OFF THE ORIGIN --> MOVE TO LZ3.
C
C     IF THERE ARE NO LZ3 CENTERS, SKIP ROTATION
      IF(NZ3.EQ.0) THEN
        GOTO 300
      ENDIF
C
C     VECTOR CENTER OF ALL LZ3
      DCNT = 0.0D0
      DO I=1,3
        CENT(I) = 0.0D0
        DO N=1,NZ3
          CENT(I) = CENT(I) + XYZ(I,IZAD(N,2))
        ENDDO
        CENT(I) = CENT(I)/NZ3
        DCNT    = DCNT + CENT(I)**2
      ENDDO
      DCNT = DSQRT(DCNT)
C
C     IF CENTER OF LZ3 IS ON YZ-AXIS BUT NOT ORIGIN, SKIP ROTATION
      X = DCNT - DSQRT(CENT(2)**2 + CENT(3)**2)
      IF(DABS(X).LT.1.0D-3.AND.DCNT.GT.1.0D-3) THEN
        GOTO 300
      ENDIF
C
C     IF CENTER OF LZ3 IS NOT ON Y-AXIS, MAKE IT THE CENTER
      X = DCNT - DABS(CENT(2))
      IF(DABS(X).GT.1.0D-3.AND.DCNT.GT.1.0D-3) THEN
        GOTO 50
      ENDIF
C
C     IF CENTER OF LZ3 IS ON Y-AXIS, FIX ON ONE OF THEM
      X = DCNT - DABS(CENT(2))
      IF(DABS(X).LT.1.0D-3.AND.DCNT.GT.1.0D-3) THEN
        DO NZ=1,NZ3
          N = IZAD(NZ,3)
          X = DIST(N) - DABS(XYZ(2,N))
          IF(DABS(X).GT.1.0D-3.AND.DIST(N).GT.1.0D-3) THEN
            DO I=1,3
              CENT(I) = XYZ(I,N)
            ENDDO
            DCNT = DIST(N)
            GOTO 50
          ENDIF
        ENDDO
      ENDIF
C
C     IF CENTER OF LZ3 IS ON ORIGIN, FIX ON ONE THAT IS OFF THE Z AXIS
      IF(DCNT.LT.1.0D-4) THEN
        DO NZ=1,NZ3
          N = IZAD(NZ,3)
          X = DIST(N) - DABS(CENT(3))
          IF(DABS(X).GT.1.0D-3.AND.DIST(N).GT.1.0D-3) THEN
            DO I=1,3
              CENT(I) = XYZ(I,N)
            ENDDO
            DCNT = DIST(N)
            GOTO 50
          ENDIF
        ENDDO
      ENDIF
C
C     SKIP POINT (ROTATION CENTER HAS BEEN IDENTIFIED)
50    CONTINUE

C     CALCULATE ROTATION ANGLE
C     ALPHA IS ROTATION ANGLE AROUND XY TO BRING CENT TO Y=0
      RLT = DSQRT(CENT(1)**2 + CENT(2)**2)
      ALPHA = DACOS(CENT(1)/RLT)
      CALL ROTATE(XYZ,NCNT,ALPHA,0.0D0)
      WRITE(6, *) 'In NUCGEOM: molecule has been rotated.'
      WRITE(7, *) 'In NUCGEOM: molecule has been rotated.'
C
C     NO Y-AXIS ROTATION NECESSARY
300   CONTINUE
C
C     CHECK WHETHER MOLECULE IS PLANAR
      DO N=1,NCNT
        TEST = XYZ(1,N)**2
        IF(TEST.GT.1.0D-6) GOTO 60
      ENDDO
      SHAPE = 'PLANAR'
60    CONTINUE
C
C     NO MOLECULAR SYMMETRY
      SHAPE = 'NOSYMM'
C
C     TRANSFER ALL TEMPORARY COORDINATES TO THE COMMON ARRAY
      DO N=1,NCNT
        DO I=1,3
          COORD(I,N) = XYZ(I,N)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     MOLECULAR GEOMETRY, BOND DISTANCES AND NUCLEAR REPULSION ENERGY  C
C**********************************************************************C
C
C     ATOMIC COORDINATES
20    FORMAT(13X,A)
21    FORMAT(1X,A,12X,A,14X,A,14X,A)
22    FORMAT(1X,I2,' (',A,') ',6X,F16.10,5X,F16.10,5X,F16.10)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,20) '  Molecular geometry A: Cartesian coordinates'
      WRITE(7,20) '  Molecular geometry A: Cartesian coordinates'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,21) 'Center','x-coord','y-coord','z-coord'
      WRITE(7,21) 'Center','x-coord','y-coord','z-coord'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO ICNT=1,NCNT
        ELA = ELMNT(IZNUC(ICNT))
        WRITE(6,22) ICNT,ELA,COORD(1,ICNT),COORD(2,ICNT),COORD(3,ICNT)
        WRITE(7,22) ICNT,ELA,COORD(1,ICNT),COORD(2,ICNT),COORD(3,ICNT)
      ENDDO
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     BOND ANGLES AND DISTANCES
30    FORMAT(21X,A)
31    FORMAT(1X,A,8X,A,12X,A,12X,A)
32    FORMAT(1X,A,2X,A,7X,F14.6)
33    FORMAT(40X,A,2X,A,2X,A,9X,F14.6)
      IF(NCNT.EQ.1) THEN
C       NUCLEAR REPULSION ENERGY
        ENUC = 0.0D0     
      ELSEIF(NCNT.GT.1) THEN
        WRITE(6, *) REPEAT(' ',72)
        WRITE(7, *) REPEAT(' ',72)
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6,30) 'Molecular geometry B: R-matrix'
        WRITE(7,30) 'Molecular geometry B: R-matrix'
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6,31) 'C1  C2','Bond distance','C1  C2  C3','Angle (deg)'
        WRITE(7,31) 'C1  C2','Bond distance','C1  C2  C3','Angle'
        WRITE(6, *) REPEAT('-',72)
        WRITE(7, *) REPEAT('-',72)
C
        ICNT = 1
        DO JCNT=2,NCNT
          ELA = ELMNT(IZNUC(ICNT))
          ELB = ELMNT(IZNUC(JCNT))
          R1X = COORD(1,JCNT) - COORD(1,ICNT)
          R1Y = COORD(2,JCNT) - COORD(2,ICNT)
          R1Z = COORD(3,JCNT) - COORD(3,ICNT)
          D1  = DSQRT(R1X*R1X + R1Y*R1Y + R1Z*R1Z)
          WRITE(6,32) ELA,ELB,D1
          WRITE(7,32) ELA,ELB,D1
C
          DO KCNT=2,JCNT-1
            ELA = ELMNT(IZNUC(ICNT))
            ELB = ELMNT(IZNUC(JCNT))
            ELC = ELMNT(IZNUC(KCNT))
            R2X = COORD(1,KCNT) - COORD(1,ICNT)
            R2Y = COORD(2,KCNT) - COORD(2,ICNT)
            R2Z = COORD(3,KCNT) - COORD(3,ICNT)
            D2  = DSQRT(R2X*R2X + R2Y*R2Y + R2Z*R2Z)
            SP  = (R1X*R2X + R1Y*R2Y + R1Z*R2Z)
            D12 = D1*D2
            DEG = 180.0D0/PI
            ANG = DEG*DACOS(SP/D12)
            WRITE(6,33) ELB,ELA,ELC,ANG  
            WRITE(7,33) ELB,ELA,ELC,ANG  
          ENDDO
          IF(JCNT.NE.NCNT) THEN
            WRITE(6, *) REPEAT('-',72)
            WRITE(7, *) REPEAT('-',72)
          ENDIF
        ENDDO
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6, *) ' '
        WRITE(7, *) ' '
C
C       NUCLEAR REPULSION ENERGY
        ENUC = 0.0D0
        DO ICNT=1,NCNT
          DO JCNT=1,ICNT-1
            SEP  = DSQRT((COORD(1,ICNT) - COORD(1,JCNT))**2
     #                  +(COORD(2,ICNT) - COORD(2,JCNT))**2
     #                  +(COORD(3,ICNT) - COORD(3,JCNT))**2)
            EPNT = ZNUC(ICNT)*ZNUC(JCNT)/SEP
C
CC          THIS CODE INCLUDES GAUSSIAN CHARGE STRUCTURE EFFECTS, 
CC          BUT CORRECTIONS EXCEED DOUBLE FLOAT ACCURACY LIMITS...
C           EPRD = CNUC(ICNT)*CNUC(JCNT)
C           ESUM = CNUC(ICNT)+CNUC(JCNT)
C           EGAU = DSQRT(EPRD/ESUM)*SEP
C           EERF = DERF(EGAU)
C           ESEP = EPNT*EERF
C
            ESEP = EPNT
            ENUC = ENUC + ESEP
            
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE FOCKIND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       FFFFFFFF OOOOOO   CCCCCC  KK    KK IIII NN    NN DDDDDDD       C
C       FF      OO    OO CC    CC KK   KK   II  NNN   NN DD    DD      C
C       FF      OO    OO CC       KK  KK    II  NNNN  NN DD    DD      C
C       FFFFFF  OO    OO CC       KKKKK     II  NN NN NN DD    DD      C
C       FF      OO    OO CC       KK  KK    II  NN  NNNN DD    DD      C
C       FF      OO    OO CC    CC KK   KK   II  NN   NNN DD    DD      C
C       FF       OOOOOO   CCCCCC  KK    KK IIII NN    NN DDDDDDD       C
C                                                                      C
C -------------------------------------------------------------------- C
C  FOCKIND ASSIGNS INDICES FOR FOCK MATRIX BLOCKS DEPENDING ON         C
C  ICNT, KQN AND MQN QUANTUM NUMBERS OF EACH BASIS FUNCTION.           C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      COMMON/QNMS/ICNBAS(MDM),KQNBAS(MDM),MQNBAS(MDM)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     QUANTUM NUMBER LABELS
      ILST = 0
      DO ICNT=1,NCNT
        DO IMV=1,(MKP+1)/2
          MQN = 2*IMV-1
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
C
      RETURN
      END
C
C
      SUBROUTINE CARTIND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        CCCCCC     AA    RRRRRRR TTTTTTTT IIII NN    NN DDDDDDD       C
C       CC    CC   AAAA   RR    RR   TT     II  NNN   NN DD    DD      C
C       CC        AA  AA  RR    RR   TT     II  NNNN  NN DD    DD      C
C       CC       AA    AA RR    RR   TT     II  NN NN NN DD    DD      C
C       CC       AAAAAAAA RRRRRRR    TT     II  NN  NNNN DD    DD      C
C       CC    CC AA    AA RR    RR   TT     II  NN   NNN DD    DD      C
C        CCCCCC  AA    AA RR    RR   TT    IIII NN    NN DDDDDDD       C
C                                                                      C
C -------------------------------------------------------------------- C
C  CARTIND GENERATES INDICES FOR EQ-COEFFICIENTS AND R-INTEGRALS.      C
C**********************************************************************C
      PARAMETER(MKP=9,ML4=2*(MKP+1),MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &             IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
C
C     DETERMINE THE MAXIMUM POSSIBLE LAM VALUE, GIVEN MKP
      LAMMAX = MKP + 6

C     SCAN THROUGH ALL POSSIBLE COMBINATIONS ALPHA, BETA, GAMMA
C     LEADING TO A GIVEN LAM VALUE AND APPLY AN ADDRESS TO EACH
      IADD = 0
      DO LAM=0,LAMMAX
        DO IA=0,LAM
          DO IB=0,LAM
            DO IC=0,LAM
              IF(IA+IB+IC.NE.LAM) GOTO 10
              IADD             = IADD + 1
              IVEC(IADD)       = IA
              JVEC(IADD)       = IB
              KVEC(IADD)       = IC
              LAMVEC(IADD)     = LAM
              INABCD(IA,IB,IC) = IADD
10            CONTINUE
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE AUFBAU(IZNUC,IQNUC,NORB,NOCC,LMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            AA    UU    UU FFFFFFF BBBBBBB     AA    UU    UU         C
C           AAAA   UU    UU FF      BB    BB   AAAA   UU    UU         C
C          AA  AA  UU    UU FF      BB    BB  AA  AA  UU    UU         C
C         AA    AA UU    UU FFFFF   BBBBBBB  AA    AA UU    UU         C
C         AAAAAAAA UU    UU FF      BB    BB AAAAAAAA UU    UU         C
C         AA    AA UU    UU FF      BB    BB AA    AA UU    UU         C
C         AA    AA  UUUUUU  FF      BBBBBBB  AA    AA  UUUUUU          C
C                                                                      C
C -------------------------------------------------------------------- C
C  AUFBAU DETERMINES THE GROUND STATE ELECTRONIC CONFIGURATION OF A    C
C  NEUTRAL ATOM OF CHARGE IZNUC, UP TO COMPLETE OCCUPATION WITH THE    C
C  LIMIT LMAX = 4 (g-ORBITALS). 220 AVAILABLE ORBITALS.                C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C    LMAX IS THE HIGHEST LQN REQUIRED TO DESCRIBE THE GROUND STATE     C
C    NOCC SAVES THE NUMBER OF OCCUPIED NSHELLS FOR THIS LQN CLASS      C
C    NORB SAVES THE NUMBER OF ELECTRONS IN OF TYPE LQN IN SHELL N      C
C  PARAMETERS:                                                         C
C    IAUF STORES THE LQN OF ORBITALS IN ORDER OF HYDROGENIC ENERGY     C
C**********************************************************************C
      PARAMETER(MKP=9)
C
      DIMENSION NORB(MKP,12),NOCC(MKP),IAUF((MKP+1)*(MKP+3)/4)
C
C     INITIALISE THE COUNTER FOR NUMBER OF ELECTRONS IN EACH ORBTIAL
      DO LQN=0,(MKP-1)/2
        NOCC(LQN+1) = 0
      ENDDO
C
C     LQN OF ATOMIC ENERGY LEVELS IN AUFBAU ORDER
      ICT = 0
      DO L=0,(MKP-1)/2
        DO M=1,2
          DO J=L,0,-1
          ICT = ICT + 1
          IAUF(ICT) = J
          ENDDO
        ENDDO
      ENDDO
C
C     INITIALISE THE MAX LQN COUNTER
      LMAX = 0
C
C     INITIALISE THE NUMBER OF ELECTRONS LEFT TO FILL ORBITALS WITH
      ILEFT = IZNUC-IQNUC
C
C     INITIALISE LOOP OVER ORBITALS
      DO M=1,30
C
C       EXIT IF THERE ARE NO ELECTRONS LEFT TO COUNT
        IF(ILEFT.EQ.0) GOTO 20
C
C       THE LQN OF THIS ORBITAL IS STORED IN IAUF
        LQN = IAUF(M)
C
C       UPDATE THE MAX LQN COUNTER IF NECESSARY
        IF(LQN.GT.LMAX) THEN
          LMAX = LQN
        ENDIF
C
C       ADD TO THE NUMBER OF OCCUPIED NSHELLS FOR THIS LQN CLASS
        NOCC(LQN+1) = NOCC(LQN+1)+1
C
C       DETERMINE NO. OF ELECTRONS REQ'D TO FULLY OCCUPY THIS SUBSHELL
        IFULL = 4*LQN + 2
C
C ***   BEGIN IF STATEMENT TO DETERMINE THE SUBSHELL OCCUPATION
        IF(ILEFT.GT.IFULL) THEN
C >>>     IF THERE ARE MORE ELECTRONS LEFT THAN IFULL, FILL THE SUBSHELL
          NORB(LQN+1,NOCC(LQN+1)) = IFULL
          ILEFT = ILEFT-IFULL
        ELSE
C >>>     OTHERWISE, LEAVE ALL REMAINING ELECTRONS IN THIS NSHELL
          NORB(LQN+1,NOCC(LQN+1)) = ILEFT
          GOTO 20
C ****  END THE NSHELL IF STATEMENT
        ENDIF
C
C     END LOOP OVER ATOMIC ORBITALS
      ENDDO
20    CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE SPECTRM(IBND,IVIR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   SSSSSS  PPPPPPP  EEEEEEEE  CCCCCC  TTTTTTTT RRRRRRR  MM       MM   C
C  SS    SS PP    PP EE       CC    CC    TT    RR    RR MMM     MMM   C
C  SS       PP    PP EE       CC          TT    RR    RR MMMM   MMMM   C
C   SSSSSS  PP    PP EEEEEE   CC          TT    RR    RR MM MM MM MM   C
C        SS PPPPPPP  EE       CC          TT    RRRRRRR  MM  MMM  MM   C
C  SS    SS PP       EE       CC    CC    TT    RR    RR MM   M   MM   C
C   SSSSSS  PP       EEEEEEEE  CCCCCC     TT    RR    RR MM       MM   C
C                                                                      C
C -------------------------------------------------------------------- C
C  SPECTRM SUMMARISES THE ENERGY LEVELS OF AN SCF ITERATION, USING     C
C  ATOMIC TERM SYMBOLS, ATOMIC CENTERS AND FRACTIONAL POPULATIONS.     C
C  IT ALSO APPLIES SYMMETRY LABELS TO EACH STATE AND SORTS MQN M'FOLDS.C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    IBND - NUMBER OF BOUND ORBITALS TO INCLUDE.                       C
C    IVIR - NUMBER OF VIRTUAL ORBITALS TO INCLUDE.                     C
C -------------------------------------------------------------------- C
C  DFNOTE: FIGURE OUT HOW TO SORT THROUGH MANIFOLD IN ATOMIC CASE.     C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*1 CHL,ELLTERM
      CHARACTER*2 ELMNT(120),ELA
      CHARACTER*4 HMLTN
C
      COMPLEX*16 C(MDM,MDM),CTEMP,ROT1,ROT2
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
C
      DIMENSION FRAC(MDM,MCT,MKP,MKP+1),NPR(MCT,MKP,MKP+1)
      DIMENSION ICNLST(MDM),KQNLST(MDM),MQNLST(MDM),NQNLST(MDM),
     &          POPLST(MDM)
C
      COMMON/ATOM/ELMNT
      COMMON/COEF/C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/QNMS/ICNBAS(MDM),KQNBAS(MDM),MQNBAS(MDM)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
      DATA PI/3.1415926535897932D0/
C
C     SENSITIVITY TOLERANCE PARAMETER
      EPS = 1.0D-14
C
C     INITIALISE COUNTER ARRAYS
      DO ICNT=1,NCNT
        DO IKAP=1,NKAP(ICNT)
          IKQN = KVALS(IKAP,ICNT)
          IF(IKQN.LT.0) THEN
            ILQN =-IKQN-1
          ELSE
            ILQN = IKQN
          ENDIF
          NMV = 2*IABS(IKQN)
          DO IMV=1,NMV
            DO IOCC=1,IBND+IVIR
              FRAC(IOCC,ICNT,IKAP,IMV) = 0.0D0
            ENDDO
            NPR(ICNT,IKAP,IMV) = ILQN
          ENDDO
        ENDDO
      ENDDO     
C
C     INTERNAL ROTATION BETWEEN PAIRS OF STATES
      DO IPAIR=1,NDIM-NSHIFT,2
C
C       SKIP NEGATIVE ENERGY SPECTRUM
        MPAIR = IPAIR+NSHIFT
C
C       TEMPORARY LARGE VALUE
        RLRG = 10.0D10
C
C       INITIAL INCRIMENTAL RADIAN (SWEEP OVER ALL POSSIBLE ANGLES)
        RINC = PI/180.0D0
C
C       SEARCH FOR STARTING POINT BY SWEEPING ANGLES 0 <= PHI < PI
        DO NPHI=0,179
C
C         CALCULATE ROTATION ANGLE AND COS/SIN TRANSFORMATIONS
          PHI  = RINC*DFLOAT(NPHI)
          CPHI = DCOS(PHI)
          SPHI = DSIN(PHI)
          OLAP = 0.0D0
C
C         ROTATE ALL THE EXPANSION COEFFICIENT MQN PAIRS BY ANGLE PHI
          DO I=1,NDIM
            ROT1 = CPHI*C(I,MPAIR  ) + SPHI*C(I,MPAIR+1)
            ROT2 =-SPHI*C(I,MPAIR  ) + CPHI*C(I,MPAIR+1)
            OLAP = OLAP + CDABS(ROT1)*CDABS(ROT2)
          ENDDO
C
C         FIND PHI WHICH RESULTS IN SMALLEST SUM OF PRODUCTS
          IF(OLAP.LT.RLRG) THEN
            PHI0 = PHI
            RLRG = OLAP
          ENDIF
        ENDDO
C
C       NEW STARTING ROTATION ANGLE BASED ON THE ABOVE SEARCH
        PHI  = PHI0 - RINC
        SOLD = 1.0D11
C
C       SWEEP THROUGH INCREMENTAL ANGLES AND SEARCH AGAIN
        DO NPHI=1,4000
          PHI  = PHI + RINC
          CPHI = DCOS(PHI)
          SPHI = DSIN(PHI)
          OLAP = 0.0D0
C
C         ROTATE ALL THE EXPANSION COEFFICIENT MQN PAIRS BY ANGLE PHI
          DO I=1,NDIM
            ROT1 = CPHI*C(I,MPAIR  ) + SPHI*C(I,MPAIR+1)
            ROT2 =-SPHI*C(I,MPAIR  ) + CPHI*C(I,MPAIR+1)
            OLAP = OLAP + CDABS(ROT1)*CDABS(ROT2)
          ENDDO
C
C         IF THE NEW VALUE IS BELOW A TOLERANCE, FINISH.
          IF(DABS(OLAP-SOLD).LT.EPS) GOTO 1
C
C         IF SUM OF COEFFICIENT PRODUCTS IS BIGGER THAN COUNTER, REFINE.
          IF(OLAP.GT.SOLD) THEN
            RINC = -RINC/10.0D0
          ENDIF
C
C         DECREASE COUNTER VALUE
          SOLD = OLAP
C
        ENDDO
1       CONTINUE
C
C       PERFORM THE ACTUAL ROTATION USING THE BEST SOLUTION PHI
        CPHI = DCOS(PHI)
        SPHI = DSIN(PHI)
        DO I=1,NDIM
          ROT1 = CPHI*C(I,MPAIR  ) + SPHI*C(I,MPAIR+1)
          ROT2 =-SPHI*C(I,MPAIR  ) + CPHI*C(I,MPAIR+1)
C
          C(I,MPAIR  ) = ROT1
          C(I,MPAIR+1) = ROT2
        ENDDO
C      
      ENDDO
C
C     CALCULATE FRACTIONAL OCCUPATION FOR EACH SOLUTION VECTOR
C
C     LOOP OVER OCCUPIED ORBITALS
      DO IOCC=1,IBND+IVIR
C
C       IGNORE NEGATIVE SPECTRUM
        MOCC = IOCC + NSHIFT
C
C       LOOP OVER FOCK MATRIX ADDRESS BLOCKS
        DO I=1,NDIM-NSHIFT
          IS = I+NSHIFT
          ICNT = ICNBAS(I)
          IKQN = KQNBAS(I)
          IMQN = MQNBAS(I)
          IF(IKQN.LT.0) THEN
            IKAP =-2*IKQN-1
          ELSE
            IKAP = 2*IKQN
          ENDIF
          IF(IMQN.LT.0) THEN
            IMV = IABS(IMQN)
          ELSE
            IMV = IMQN+1
          ENDIF
C
          DO J=1,NDIM-NSHIFT
            JS = J+NSHIFT
            JCNT = ICNBAS(J)
            JKQN = KQNBAS(J)
            JMQN = MQNBAS(J)
C
C           LARGE AND SMALL CONTRIBUTIONS
            TMP = DREAL(DCONJG(C(I ,MOCC))*C(J ,MOCC)*OVAP(I ,J ))
     &          + DREAL(DCONJG(C(IS,MOCC))*C(JS,MOCC)*OVAP(IS,JS))
C
C           DECIDE WHERE TO PUT CONTRIBUTION
            IF(ICNT.EQ.JCNT.AND.IKQN.EQ.JKQN.AND.IMQN.EQ.JMQN) THEN
              FRAC(IOCC,ICNT,IKAP,IMV) = FRAC(IOCC,ICNT,IKAP,IMV) + TMP
            ENDIF
C
          ENDDO
        ENDDO
      ENDDO
C
C     IDENTIFY CONVENTIONAL DIRAC QUANTUM NUMBERS BASED ON FRAC.
C
C     LOOP OVER POSITIVE ENERGY SPECTRUM IN PAIRS
      DO IOCC=1,NDIM-NSHIFT
C
C       SEARCH FOR HIGHEST POPULATED ATOMIC STATE
        PLTN = 0.0D0
        DO ICNT=1,NCNT
          DO IKAP=1,NKAP(ICNT)
            IKQN = KVALS(IKAP,ICNT)
            NMV  = 2*IABS(IKQN)
            DO IMV=1,NMV
              IF(FRAC(IOCC,ICNT,IKAP,IMV).GT.PLTN) THEN
                PLTN = FRAC(IOCC,ICNT,IKAP,IMV)
                KCNT = ICNT
                KKQN = IKQN
                KKAP = IKAP
                KMV  = IMV
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        NPR(KCNT,KKAP,KMV) = NPR(KCNT,KKAP,KMV)+1
C
C       APPLY DIRAC LABELS ACCORDING TO HIGHEST POPULATED STATE
        NQNLST(IOCC) = NPR(KCNT,KKAP,KMV)
        ICNLST(IOCC) = KCNT
        KQNLST(IOCC) = KKQN
        POPLST(IOCC) = FRAC(IOCC,KCNT,KKAP,KMV)
        IF(MOD(KMV,2).EQ.1) THEN
          MQNLST(IOCC) =-KMV
        ELSE
          MQNLST(IOCC) = KMV-1
        ENDIF
C
C     END LOOP OVER ORBITALS
      ENDDO
C
C     SORTING: ORGANISE POSITIVE-ENERGY SOLUTIONS IN MQN PAIRS
C
      IF(NCNT.EQ.1) THEN
C
C       LOOP OVER ALL STATES
        DO IOCC=1,NDIM-NSHIFT
C
C         SKIP NEGATIVE ENERGY SPECTRUM
          MOCC = IOCC+NSHIFT
C
C         IDENTIFY KQN OF STATE
          IKQN = KQNLST(IOCC)
C
C         FOR THIS KQN VALUE, SEARCH NEXT 2*|KQN| ENTRIES AND ORDER
C
        ENDDO
C
      ELSEIF(NCNT.GT.1) THEN
C
C       LOOP OVER PAIRS OF STATES
        DO IPAIR=1,NDIM-NSHIFT,2
C
C         SKIP NEGATIVE ENERGY SPECTRUM
          MPAIR = IPAIR+NSHIFT
C
C         NO NEED TO SWAP IF MQN OF THE FIRST STATE IS NEGATIVE
          IF(MQNLST(IPAIR).LT.0) GOTO 100
C
C         SWAP EIGENVALUES, EXPANSION COEFFICIENTS AND MQN VALUES
          ETEMP          = EIGEN(MPAIR+1)
          EIGEN(MPAIR+1) = EIGEN(MPAIR  )
          EIGEN(MPAIR  ) = ETEMP
C
          DO I=1,NDIM
            CTEMP        = C(I,MPAIR+1)
            C(I,MPAIR+1) = C(I,MPAIR  )
            C(I,MPAIR  ) = CTEMP
          ENDDO
C
          MTEMP           = MQNLST(IPAIR+1)
          MQNLST(IPAIR+1) = MQNLST(IPAIR  )
          MQNLST(IPAIR  ) = MTEMP
C
          PTEMP           = POPLST(IPAIR+1)
          POPLST(IPAIR+1) = POPLST(IPAIR  )
          POPLST(IPAIR  ) = PTEMP
C
C         SKIP POINT FOR ALREADY-SORTED PAIRS
100       CONTINUE
C      
        ENDDO
C
      ENDIF
C
C     PRINT TITLE TO TERMINAL/FILE
20    FORMAT(1X,'Orb.',2X,'Center',4X,'Term sym.',3X,'m_j',14X,
     &                                    'Energy (au)',6X,'Population')
21    FORMAT(1X,I3,2X,I2,' (',A,')',4X,I2,A,'_',I1,'/2',4X,I2,'/2',7X,
     &                                                  F18.12,6X,F10.8)
C
      WRITE(6,20)
      WRITE(7,20)
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     SUMMARISE RESULTS
      DO IOCC=1,IBND+IVIR
C
        MOCC = IOCC+NSHIFT
C
        ICNT = ICNLST(IOCC)
        INQN = NQNLST(IOCC)
        IKQN = KQNLST(IOCC)
        IMQN = MQNLST(IOCC)
        ELA  = ELMNT(IZNUC(ICNT))
        IF(IKQN.LT.0) THEN
          ILQN =-IKQN-1
        ELSE
          ILQN = IKQN
        ENDIF
        CHL  = ELLTERM(ILQN)
        IJQN = 2*IABS(IKQN)-1
        
        IF(IOCC.LE.IBND) THEN
          PLTN = POPLST(IOCC)
        ELSE
          PLTN = 0.0D0
        ENDIF
C
        IF(HMLTN.EQ.'NORL') THEN
          PLTN = 0.5D0*PLTN
        ENDIF
C
C       OUTPUT TO TERMINAL
        WRITE(6,21) IOCC,ICNT,ELA,INQN,CHL,IJQN,IMQN,EIGEN(MOCC),PLTN
        WRITE(7,21) IOCC,ICNT,ELA,INQN,CHL,IJQN,IMQN,EIGEN(MOCC),PLTN
        IF(IOCC.EQ.IBND) THEN
          WRITE(6, *) REPEAT('-',72)
          WRITE(7, *) REPEAT('-',72)
        ENDIF
C
      ENDDO
C
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
      RETURN
      END
C
C
      FUNCTION ELLTERM(L)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   EEEEEEEE LL       LL     TTTTTTTT EEEEEEEE RRRRRRR  MM       MM    C
C   EE       LL       LL        TT    EE       RR    RR MMM     MMM    C
C   EE       LL       LL        TT    EE       RR    RR MMMM   MMMM    C
C   EEEEEE   LL       LL        TT    EEEEEE   RR    RR MM MM MM MM    C
C   EE       LL       LL        TT    EE       RRRRRRR  MM  MMM  MM    C
C   EE       LL       LL        TT    EE       RR    RR MM   M   MM    C
C   EEEEEEEE LLLLLLLL LLLLLLLL  TT    EEEEEEEE RR    RR MM       MM    C
C                                                                      C
C -------------------------------------------------------------------- C
C  ELLTERM(L) EVALUATES A CHARACTER CORRESPONDING TO TERM SYMBOL L.    C
C**********************************************************************C
C
      CHARACTER*1 ELLTERM
C
      IF(L.EQ.0) THEN
        ELLTERM = 's'
      ELSEIF(L.EQ.1) THEN
        ELLTERM = 'p'
      ELSEIF(L.EQ.2) THEN
        ELLTERM = 'd'
      ELSEIF(L.EQ.3) THEN
        ELLTERM = 'f'
      ELSEIF(L.EQ.4) THEN
        ELLTERM = 'g'
      ELSEIF(L.EQ.5) THEN
        ELLTERM = 'h'
      ELSEIF(L.EQ.6) THEN
        ELLTERM = 'i'
      ELSEIF(L.EQ.7) THEN
        ELLTERM = 'j'
      ELSEIF(L.EQ.8) THEN
        ELLTERM = 'k'
      ELSEIF(L.EQ.9) THEN
        ELLTERM = 'l'
      ELSEIF(L.EQ.10) THEN
        ELLTERM = 'm'
      ELSEIF(L.EQ.11) THEN
        ELLTERM = 'n'
      ELSEIF(L.EQ.12) THEN
        ELLTERM = 'o'
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE ROTATE(XYZ,NCNT,ALPHA,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          RRRRRRR   OOOOOO TTTTTTTT   AA   TTTTTTTT EEEEEEEE          C
C          RR    RR OO    OO   TT     AAAA     TT    EE                C
C          RR    RR OO    OO   TT    AA  AA    TT    EE                C
C          RR    RR OO    OO   TT   AA    AA   TT    EEEEEE            C
C          RRRRRRR  OO    OO   TT   AAAAAAAA   TT    EE                C
C          RR    RR OO    OO   TT   AA    AA   TT    EE                C
C          RR    RR  OOOOOO    TT   AA    AA   TT    EEEEEEEE          C
C                                                                      C
C -------------------------------------------------------------------- C
C  ROTATE PERFORMS TWO ROTATIONS ON ALL ATOMIC CENTERS, FIRST SO THAT  C
C  Y DOES NOT CHANGE, AND THEN SO THAT Z DOES NOT CHANGE.              C
C**********************************************************************C
C
      DIMENSION XYZ(3,NCNT),V(3),W(3)
      DIMENSION AR(3,3),BR(3,3),RR(3,3)
C
C     INITIALISE MATRICES
      DO I=1,3
        DO J=1,3
          AR(I,J) = 0.0D0
          BR(I,J) = 0.0D0
        ENDDO
      ENDDO
C
C     FIX ALL VALUES OF Z AND ROTATE BY ALPHA
      AR(3,3) = 1.0D0
      AR(1,1) = DCOS(ALPHA)
      AR(2,2) = DCOS(ALPHA)
      AR(1,2) = DSIN(ALPHA)
      AR(2,1) =-DSIN(ALPHA)
C
C     FIX ALL VALUES OF Y AND ROTATE BY BETA
      BR(2,2) = 1.0D0
      BR(1,1) = DCOS(BETA)
      BR(3,3) = DCOS(BETA)
      BR(3,1) = DSIN(BETA)
      BR(1,3) =-DSIN(BETA)
C
C     EVALUATE OVERALL ROTATION MATRIX
      CALL MMPROD(BR,AR,RR,3)     
C
C     ROTATE ALL VECTORS IN XYZ
      DO N=1,NCNT
        DO I=1,3
          V(I) = XYZ(I,N)
        ENDDO
        CALL MVPROD(RR,V,W,3)
        DO I=1,3
          XYZ(I,N) = W(I)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE MMPROD(A,B,X,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     MM       MM MM       MM PPPPPPP  RRRRRRR   OOOOOO  DDDDDDD       C
C     MMM     MMM MMM     MMM PP    PP RR    RR OO    OO DD    DD      C
C     MMMM   MMMM MMMM   MMMM PP    PP RR    RR OO    OO DD    DD      C
C     MM MM MM MM MM MM MM MM PP    PP RR    RR OO    OO DD    DD      C
C     MM  MMM  MM MM  MMM  MM PPPPPPP  RRRRRRR  OO    OO DD    DD      C
C     MM   M   MM MM   M   MM PP       RR    RR OO    OO DD    DD      C
C     MM       MM MM       MM PP       RR    RR  OOOOOO  DDDDDDD       C
C                                                                      C
C -------------------------------------------------------------------- C
C  MMPROD CALCULATES THE PRODUCT OF TWO SQUARE, DOUBLE-PRECISION       C
C  MATRIX ARRAYS OF DIMENSION N, AND OUTPUTS THE RESULT INTO X.        C
C**********************************************************************C
C
      DIMENSION A(N,N),B(N,N),X(N,N)
C
C     INITIALISE X MATRIX
      DO I=1,N
        DO J=1,N
          X(I,J) = 0.0D0
        ENDDO
      ENDDO
C
C     MATRIX PRODUCT
      DO I=1,N
        DO J=1,N
          DO K=1,N
            X(I,J) = X(I,J) + A(I,K)*B(K,J)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE MVPROD(A,V,W,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       MM       MM VV    VV PPPPPPP  RRRRRRR   OOOOOO  DDDDDDD        C
C       MMM     MMM VV    VV PP    PP RR    RR OO    OO DD    DD       C
C       MMMM   MMMM VV    VV PP    PP RR    RR OO    OO DD    DD       C
C       MM MM MM MM VV    VV PP    PP RR    RR OO    OO DD    DD       C
C       MM  MMM  MM  VV  VV  PPPPPPP  RRRRRRR  OO    OO DD    DD       C
C       MM   M   MM   VVVV   PP       RR    RR OO    OO DD    DD       C
C       MM       MM    VV    PP       RR    RR  OOOOOO  DDDDDDD        C
C                                                                      C
C -------------------------------------------------------------------- C
C  MVPROD CALCULATES THE PRODUCT OF A SQUARE MATRIX OF DIMENSION N     C
C  AND A VECTOR OF DIMENSION N, AND OUTPUTS THE RESULT INTO W.         C
C**********************************************************************C
C
      DIMENSION A(N,N),V(N),W(N)
C
C     INITIALISE W MATRIX
      DO I=1,N
        W(I) = 0.0D0
      ENDDO
C
C     MATRIX PRODUCT
      DO I=1,N
        DO J=1,N
          W(I) = W(I) + A(I,J)*V(J)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION HMS(TSEC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                     HH    HH MM       MM  SSSSSS                     C
C                     HH    HH MMM     MMM SS    SS                    C
C                     HH    HH MMMM   MMMM SS                          C
C                     HHHHHHHH MM MM MM MM  SSSSSS                     C
C                     HH    HH MM  MMM  MM       SS                    C
C                     HH    HH MM   M   MM SS    SS                    C
C                     HH    HH MM       MM  SSSSSS                     C
C                                                                      C
C -------------------------------------------------------------------- C
C   HMS RETURNS A QUOTED TIME IN SECONDS TO 'HR-MIN-SEC' FORMAT.       C
C**********************************************************************C
C
      CHARACTER*4 HOURS,MINUTES
      CHARACTER*8 SECONDS
      CHARACTER*16 HMS
C     
C     INITIALISE COUNTERS
      NMIN = 0
      NHRS = 0
C
C     BACKUP TSEC
      TEMP = TSEC
C
C     PERFORM MODULAR ARITHMETIC UNTIL 0.0D0 <= TSEC < 60.0D0
      DO WHILE (TSEC.GE.6.0D1)
        TSEC = TSEC - 60.0D0
        NMIN = NMIN + 1
      ENDDO
C
C     PERFORM MODULAR ARITHMETIC UNTIL 0 <= NMIN < 60
      DO WHILE (NMIN.GE.60)
        NMIN = NMIN - 60
        NHRS = NHRS + 1
      ENDDO
C
C     PRINT THE HOUR
      IF(NHRS.EQ.0) THEN
        WRITE(HOURS,'(A)') '    '
      ELSE
        WRITE(HOURS,'(I3,A)') NHRS,'h'
      ENDIF
C
C     PRINT THE MINUTE
      IF(NHRS.EQ.0) THEN
        IF(NMIN.EQ.0) THEN
          WRITE(MINUTES,'(1X,A)') '   '
        ELSE
          WRITE(MINUTES,'(1X,I2,A)') NMIN,'m'
        ENDIF
      ELSEIF(NHRS.NE.0) THEN
        IF(NMIN.LT.10) THEN
          WRITE(MINUTES,'(1X,I1,I1,A)') 0,NMIN,'m'
        ELSE
          WRITE(MINUTES,'(1X,I2,A)') NMIN,'m'
        ENDIF
      ENDIF
C
C     PRINT THE SECOND
      IF(NMIN.EQ.0) THEN
        IF(TSEC.LT.1.0D0) THEN
          WRITE(SECONDS,'(2X,F5.3,A)') TSEC,'s'
        ELSE
          WRITE(SECONDS,'(1X,F6.3,A)') TSEC,'s'
        ENDIF
      ELSEIF(NMIN.GT.0) THEN
        IF(TSEC.LT.10.0D0) THEN
          WRITE(SECONDS,'(1X,I1,F5.3,A)') 0,TSEC,'s'
        ELSE
          WRITE(SECONDS,'(1X,F6.3,A)') TSEC,'s'
        ENDIF
      ENDIF
C
      WRITE(HMS,'(A,A,A)') HOURS//MINUTES//SECONDS
C
C     RESTORE TSEC
      TSEC = TEMP
C
      RETURN
      END
C
C
      SUBROUTINE TIMENOW(STAMP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  TTTTTTTT IIII MM       MM EEEEEEEE NN    NN  OOOOOO  WW         WW  C
C     TT     II  MMM     MMM EE       NNN   NN OO    OO WW         WW  C
C     TT     II  MMMM   MMMM EE       NNNN  NN OO    OO WW         WW  C
C     TT     II  MM MM MM MM EEEEEE   NN NN NN OO    OO WW    W    WW  C
C     TT     II  MM  MMM  MM EE       NN  NNNN OO    OO WW   WWW   WW  C
C     TT     II  MM   M   MM EE       NN   NNN OO    OO  WW WW WW WW   C
C     TT    IIII MM       MM EEEEEEEE NN    NN  OOOOOO    WW     WW    C
C                                                                      C
C -------------------------------------------------------------------- C
C  TIMENOW CREATES A DATE STRING SPECIFYING CPU TIME WHEN CALLED.      C
C**********************************************************************C
C
      CHARACTER*5  ZONE
      CHARACTER*8  DATE
      CHARACTER*10 TIME
      CHARACTER*20 STAMP
C      
      DIMENSION IVL(8)
C
C     CALL TIME AND DATE ROUTINE
      CALL DATE_AND_TIME(DATE,TIME,ZONE,IVL)
C
C     PRINT THE DAY
      IF(IVL(3).LT.10) THEN
        WRITE(STAMP,'(1X,I1,I1,A)') 0,IVL(3),'/'
      ELSE
        WRITE(STAMP,'(1X,I2,A)') IVL(3),'/'
      ENDIF
C
C     PRINT THE MONTH
      IF(IVL(2).LT.10) THEN
        WRITE(STAMP,'(A,I1,I1,A)') TRIM(STAMP),0,IVL(2),'/'
      ELSE
        WRITE(STAMP,'(A,I2,A)') TRIM(STAMP),IVL(2),'/'
      ENDIF
C
C     PRINT THE YEAR
      WRITE(STAMP,'(A,I4)') TRIM(STAMP),IVL(1)
C
C     PRINT THE HOUR
      IF(IVL(5).LT.10) THEN
        WRITE(STAMP,'(A,A,I1,I1,A)') TRIM(STAMP),' ',0,IVL(5),':'
      ELSE
        WRITE(STAMP,'(A,A,I2,A)') TRIM(STAMP),' ',IVL(5),':'
      ENDIF
C
C     PRINT THE MINUTE
      IF(IVL(6).LT.10) THEN
        WRITE(STAMP,'(A,I1,I1,A)') TRIM(STAMP),0,IVL(6),':'
      ELSE
        WRITE(STAMP,'(A,I2,A)') TRIM(STAMP),IVL(6),':'
      ENDIF
C
C     PRINT THE SECOND
      IF(IVL(7).LT.10) THEN
        WRITE(STAMP,'(A,I1,I1)') TRIM(STAMP),0,IVL(7)
      ELSE
        WRITE(STAMP,'(A,I2)') TRIM(STAMP),IVL(7)
      ENDIF
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C   [3] DENSITIES: MOLECULAR DENSITIES, ENERGIES AND LEVEL SHIFTING.   C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] DENSTY0: GENERATES CLOSED- AND OPEN-SHELL MOLECULAR DENSITY.   C
C   [B] DENSTY: GENERATES CLOSED- AND OPEN-SHELL MOLECULAR DENSITY.    C
C   [C] LEVSHFT: APPLIES A LEVEL SHIFT TO OCCUPIED ORBITALS IN FOCK.   C
C   [D] ENERGIES: USE DENSITY MATRIX TO CALCULATE ENERGY TERMS.        C
C**********************************************************************C
C
C
      SUBROUTINE DENSTY0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    DDDDDDD  EEEEEEEE NN    NN  SSSSSS TTTTTTTT YY    YY  000000      C
C    DD    DD EE       NNN   NN SS    SS   TT    YY    YY 00    00     C
C    DD    DD EE       NNNN  NN SS         TT     YY  YY  00    00     C
C    DD    DD EEEEEE   NN NN NN  SSSSSS    TT      YYYY   00    00     C
C    DD    DD EE       NN  NNNN       SS   TT       YY    00    00     C
C    DD    DD EE       NN   NNN SS    SS   TT       YY    00    00     C
C    DDDDDDD  EEEEEEEE NN    NN  SSSSSS    TT       YY     000000      C
C                                                                      C
C -------------------------------------------------------------------- C
C  DENSTY0 IS A STARTING DENSITY ROUTINE FOR USE ONLY WHEN INEW = 0.   C
C  THIS IS BECAUSE 'ATOMIC' CALCULATES AVERAGE OVER SHELL OCCUPANCIES. C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM),C(MDM,MDM)
      COMPLEX*16 SUM
C
      COMMON/COEF/C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/OCPD/IOCPN(MDM),IOCCM0
      COMMON/SHLL/ALPHA,BETA,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN,NOELEC
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     CONSTRUCT THE CLOSED-SHELL AND TOTAL DENSITY ARRAYS
      DO I=1,NDIM
        DO J=1,NDIM
          SUM = DCMPLX(0.0D0,0.0D0)
          DO IOCC=1,IOCCM0
            SUM = SUM + DCONJG(C(I,IOCC+NSHIFT))*C(J,IOCC+NSHIFT)
          ENDDO
          DENC(I,J) = SUM
          DENO(I,J) = 0.0D0
          DENT(I,J) = DENC(I,J) + DENO(I,J)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE DENSTY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         DDDDDDD  EEEEEEEE NN    NN  SSSSSS TTTTTTTT YY    YY         C
C         DD    DD EE       NNN   NN SS    SS   TT    YY    YY         C
C         DD    DD EE       NNNN  NN SS         TT     YY  YY          C
C         DD    DD EEEEEE   NN NN NN  SSSSSS    TT      YYYY           C
C         DD    DD EE       NN  NNNN       SS   TT       YY            C
C         DD    DD EE       NN   NNN SS    SS   TT       YY            C
C         DDDDDDD  EEEEEEEE NN    NN  SSSSSS    TT       YY            C
C                                                                      C
C -------------------------------------------------------------------- C
C  DENSTY GENERATES DENSITY MATRICES FROM THE EXPANSION COEFFS C.      C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM),C(MDM,MDM)
      COMPLEX*16 SUM
C
      COMMON/COEF/C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/SHLL/ALPHA,BETA,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN,NOELEC
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     MAKE CLOSED-SHELL DENSITY AND EMPTY OPEN-SHELL DENSITY (RSCF 81)
      DO I=1,NDIM
        DO J=1,NDIM
          SUM = DCMPLX(0.0D0,0.0D0)
          DO IOCC=1,NCLS
            ICL = ICLS(IOCC)
            SUM = SUM + DCONJG(C(I,ICL+NSHIFT))*C(J,ICL+NSHIFT)
          ENDDO
          DENC(I,J) = SUM
          DENO(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     MAKE THE OPEN-SHELL DENSITY (RSCF 82)
      IF(NOPN.NE.0) THEN
        DO I=1,NDIM
          DO J=1,NDIM
            SUM = DCMPLX(0.0D0,0.0D0)
            DO IOCC=1,NOPN
              IOP = IOPN(IOCC)
              SUM = SUM + FOPEN*DCONJG(C(I,IOP+NSHIFT))*C(J,IOP+NSHIFT)
            ENDDO
            DENO(I,J) = SUM
          ENDDO
        ENDDO
      ENDIF
C
C     MAKE THE TOTAL DENSITY MATRIX (RSCF 83)
      DO I=1,NDIM
        DO J=1,NDIM
          DENT(I,J) = DENC(I,J) + DENO(I,J)
        ENDDO
      ENDDO
C     
      RETURN
      END
C
C
      SUBROUTINE LEVSHFT(SHLV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   LL       EEEEEEEE VV      VV  SSSSSS  HH    HH FFFFFFFF TTTTTTTT   C
C   LL       EE       VV      VV SS    SS HH    HH FF          TT      C
C   LL       EE       VV      VV SS       HH    HH FF          TT      C
C   LL       EEEEEE    VV    VV   SSSSSS  HHHHHHHH FFFFFF      TT      C
C   LL       EE         VV  VV         SS HH    HH FF          TT      C
C   LL       EE          VVVV    SS    SS HH    HH FF          TT      C
C   LLLLLLLL EEEEEEEE     VV      SSSSSS  HH    HH FF          TT      C
C                                                                      C
C -------------------------------------------------------------------- C
C  LEVSHFT APPLIES A LEVEL SHIFT OF SHLV TO UNOCCUPIED ORBITALS IN THE C
C  FOCK MATRIX.                                                        C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 C(MDM,MDM),A(MDM),SUM
C
      COMMON/COEF/C
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     LOOP OVER ALL NON-OCCUPIED POSITIVE-ENERGY (I.E., VIRTUAL) OCCS
      DO IVIR=NSHIFT+NOCC+1,NDIM
        DO I=1,NDIM
          SUM = DCMPLX(0.0D0,0.0D0)
          DO J=1,NDIM
            SUM = SUM + OVAP(I,J)*C(J,IVIR)
          ENDDO
          A(I) = SUM
        ENDDO
C
        DO I=1,NDIM
          DO J=1,NDIM
            FOCK(I,J) = FOCK(I,J) + SHLV*A(I)*DCONJG(A(J))
          ENDDO
        ENDDO
      ENDDO
C      
      RETURN
      END
C
C
      SUBROUTINE ENERGIES
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  EEEEEEEE NN    NN EEEEEEEE RRRRRRR   GGGGGG IIII EEEEEEEE SSSSSS    C
C  EE       NNN   NN EE       RR    RR GG    GG II  EE      SS    SS   C
C  EE       NNNN  NN EE       RR    RR GG       II  EE      SS         C
C  EEEEEE   NN NN NN EEEEEE   RR    RR GG       II  EEEEEE   SSSSSS    C
C  EE       NN  NNNN EE       RRRRRRR  GG   GGG II  EE            SS   C
C  EE       NN   NNN EE       RR    RR GG    GG II  EE      SS    SS   C
C  EEEEEEEE NN    NN EEEEEEEE RR    RR  GGGGGG IIII EEEEEEEE SSSSSS    C
C                                                                      C
C -------------------------------------------------------------------- C
C  ENERGIES CALCULATES INTERACTION ENERGIES OF THE CURRENT DENSITY     C
C  MATRIX OVER OCCUPIED SPINORS WITH THE MATRIX REP OF OPERATORS.      C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 ETMP(12)
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EMDR,EMXC,EUEH
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/SHLL/ALPHA,BETA,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN,NOELEC
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     INITIALISE MOLECULAR ENERGY COUNTERS
      DO N=1,12
        ETMP(N) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     CALCULATE MOLECULAR ENERGY COUNTERS OVER THE MOST RECENT DENSITY
      DO I=1,NDIM
        DO J=1,NDIM
          ETMP( 1) = ETMP( 1) +       DENT(I,J)*HNUC(I,J)
          ETMP( 2) = ETMP( 2) +       DENT(I,J)*HKIN(I,J)
          ETMP( 3) = ETMP( 3) +       DENT(I,J)*VUEH(I,J)
          ETMP( 4) = ETMP( 4) + 0.5D0*DENT(I,J)*GDIR(I,J)
          ETMP( 5) = ETMP( 5) + 0.5D0*DENT(I,J)*GXCH(I,J)
          ETMP( 6) = ETMP( 6) + 0.5D0*DENT(I,J)*BDIR(I,J)
          ETMP( 7) = ETMP( 7) + 0.5D0*DENT(I,J)*BXCH(I,J)
          ETMP( 8) = ETMP( 8) - 0.5D0*QDIR(I,J)*DENO(I,J)
     &                        + 0.5D0*QDIR(I,J)*DENT(I,J)*(FOPEN-1.0D0)
          ETMP( 9) = ETMP( 9) - 0.5D0*QXCH(I,J)*DENO(I,J)
     &                        + 0.5D0*QXCH(I,J)*DENT(I,J)*(FOPEN-1.0D0)
          ETMP(12) = ETMP(12) +       DENT(I,J)*FOCK(I,J)
        ENDDO
      ENDDO
C
C     REAL COMPONENTS ARE ACTUAL ENERGIES
      EHNC = DREAL(ETMP(1))
      EHKN = DREAL(ETMP(2))
      EUEH = DREAL(ETMP(3))
      EGDR = DREAL(ETMP(4))
      EGXC = DREAL(ETMP(5))
      EBDR = DREAL(ETMP(6))
      EBXC = DREAL(ETMP(7))
      EQDR = DREAL(ETMP(8))
      EQXC = DREAL(ETMP(9))
      EMDR = DREAL(ETMP(10))
      EMXC = DREAL(ETMP(11))
C
C     ADD ALL CONTRIBUTIONS TO THE TOTAL ENERGY
      EONE = EHNC + EHKN
      ECLG = EGDR - EGXC
      ECLQ = EQDR - EQXC
      EBRG = EBDR - EBXC
      EBRQ = EMDR - EMXC
      ETOT = ENUC + EONE +EUEH + ECLG + ECLQ + EBRG
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C   [4] ATOMIC HARTREE-FOCK: SINGLE-CENTER SCF CALCULATIONS.           C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] ATOMIC: MAIN ROUTNE FOR ATOMIC SCF CALCULATIONS.               C
C   [B] HFSCF0: ATOMIC SCF ROUTINE (AVERAGE OF CONFIGURATION MODEL.)   C
C   [C] OVRLP0: GENERATES ATOMIC OVERLAP MATRIX.                       C
C   [D] ONEEL0: GENERATES ATOMIC ONE-ELECTRON MATRIX (ALL HAMILS).     C
C   [E] BNUCINT0: CALCULATES AN ATOMIC BARE NUCLEUS INTEGRAL.          C
C   [F] COULOMB0: ATOMIC MEAN-FIELD COULOMB MATRIX (BETA INTEGRALS).   C
C   [G] BREIT0: CONSTRUCTION OF ATOMIC BREIT INTERACTION MATRIX.       C
C   [H] KLSET: ADDRESSES AND 'KL' EXPONENT POWERS FOR ATOMIC ERI/BII.  C
C   [I] BINT: BATCH OF BETA INTEGRALS FOR ATOMIC ERI/BII.              C
C   [J] ERI0: BATCH OF ATOMIC COULOMB INTERACTION INTEGRALS.           C
C   [K] BII0: BATCH OF ATOMIC BREIT INTERACTION INTEGRALS.             C
C   [L] ANGCOUL: ATOMIC ANGULAR COULOMB COEFFICIENTS.                  C
C   [M] ANGBREIT: ATOMIC ANGULAR BREIT COEFFICIENTS.                   C
C   [N] EXCHNG: EXCHANGE MAGNETIC COEFFICIENTS FOR ANGBREIT.           C
C   [O] ABC000: 3-J SYMBOLS FOR USE IN NON-REL ANGCOUL (LS COUPLING).  C
C   [P] SYM3JSQ: 3-J SYMBOLS FOR USE IN ANGCOUL AND ANGBREIT (JJ).     C
C   [Q] UEHLING0: GENERATES ATOMIC UEHLING INTERACTION MATRIX.         C
C   [R] UEHINT0: CALCULATES AN ATOMIC UEHLING INTEGRAL.                C
C   [S] RNORM0: GENERATE BATCHES OF ALL TT' NORMALISATION FACTORS.     C
C   [T] GAMGEN: LIST OF GAMMA FUNCTIONS FOR INT AND HALF-INT ARGS.     C
C   [U] FACTRLS: LIST OF FACTORIALS AND DOUBLE FACTORIALS.             C
C**********************************************************************C
C
C
      SUBROUTINE ATOMIC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             AA   TTTTTTTT  OOOOOO  MM       MM IIII  CCCCCC          C
C            AAAA     TT    OO    OO MMM     MMM  II  CC    CC         C
C           AA  AA    TT    OO    OO MMMM   MMMM  II  CC               C
C          AA    AA   TT    OO    OO MM MM MM MM  II  CC               C
C          AAAAAAAA   TT    OO    OO MM  MMM  MM  II  CC               C
C          AA    AA   TT    OO    OO MM   M   MM  II  CC    CC         C
C          AA    AA   TT     OOOOOO  MM       MM IIII  CCCCCC          C
C                                                                      C
C -------------------------------------------------------------------- C
C  ATOMIC PERFORMS A SINGLE-CENTER SCF PROCEDURE FOR EACH ATOM IN THE  C
C  MOLECULE AND ASSEMBLES AN INITIAL COEFFICIENT MATRIX.               C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4  HMLTN
      CHARACTER*16 HMS
      CHARACTER*20 STAMP
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 C(MDM,MDM)
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
C
      COMMON/COEF/C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EMDR,EMXC,EUEH
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/OCPD/IOCPN(MDM),IOCCM0
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SHLL/ACFF,BCFF,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN,NOELEC
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
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
C     SPECIAL EXIT FOR ONE-ELECTRON PROBLEMS
      IF(NOELEC.LE.1.AND.NCNT.NE.1) THEN
        WRITE(6, *) 'This is a one-electron problem. Skip ATOMIC.'
        WRITE(7, *) 'This is a one-electron problem. Skip ATOMIC.'
        RETURN
      ENDIF
C
C     RECORD TIME AT START OF ATOMIC CALCULATION
      CALL CPU_TIME(TDUM)
C
C     INITIALISE MOLECULAR MATRICES
      DO I=1,NDIM
        DO J=1,NDIM
          C(I,J)    = DCMPLX(0.0D0,0.0D0)
          OVAP(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     CALCULATE GAMMA FUNCTION VALUES FOR LATER USE
      CALL GAMGEN
C
C     INITIALISE OCCUPATION COUNTER
      IOCCM0 = 0
c
C     HARTREE-FOCK SCF PROCEDURE FOR EACH ISOLATED ATOM
      DO ICNT=1,NCNT
        CALL HFSCF0(ICNT)
      ENDDO
C
C     SAVE EIGENVECTORS TO OUTPUT FILE
      OPEN(UNIT=8,FILE=TRIM(WFNFL),STATUS='UNKNOWN')
      REWIND(UNIT=10)
      DO I=1,NDIM
        WRITE(8, *) EIGEN(I),(C(J,I),J=1,NDIM)
      ENDDO
      CLOSE(UNIT=8)
c
C     GENERATE DENSITY MATRIX
      CALL DENSTY0
C
C     TIME TAKEN FOR ATOMIC CALCULATION
      CALL CPU_TIME(TATM)
      TATM = TATM - TDUM
C
C     DATE AND TIME AT END OF ITERATION
      CALL TIMENOW(STAMP)
C
C     ATOMIC SCF SUMMARY HEADER
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
      WRITE(6, *) REPEAT(' ',27),'Atomic SCF summary'
      WRITE(7, *) REPEAT(' ',27),'Atomic SCF summary'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     PRODUCE SPECTRUM SUMMARY AND INCLUDE FIRST 6 VIRTUAL STATES
      CALL SPECTRM(IOCCM0,6)
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
      IF(HMLTN.EQ.'BARE') GOTO 202
      WRITE(6,20) 'Coulomb (closed) (G)',ECLG
      WRITE(7,20) 'Coulomb (closed) (G)',ECLG
      IF(HMLTN.EQ.'NORL'.OR.HMLTN.EQ.'DHFR') GOTO 202
      WRITE(6,20) 'Breit (closed)   (B)',EBRG
      WRITE(7,20) 'Breit (closed)   (B)',EBRG
202   CONTINUE
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) 'Molecule total   (F)',ETOT
      WRITE(7,20) 'Molecule total   (F)',ETOT
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     CALCULATION TIME
30    FORMAT(1X,A,26X,A)
      WRITE(6,30) 'Total atomic SCF time         ',HMS(TATM)
      WRITE(7,30) 'Total atomic SCF time         ',HMS(TATM)
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
C
C
      SUBROUTINE HFSCF0(ICNT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          HH    HH FFFFFFFF SSSSSS   CCCCCC  FFFFFFFF 000000          C
C          HH    HH FF      SS    SS CC    CC FF      00    00         C
C          HH    HH FF      SS       CC       FF      00    00         C
C          HHHHHHHH FFFFFF   SSSSSS  CC       FFFFFF  00    00         C
C          HH    HH FF            SS CC       FF      00    00         C
C          HH    HH FF      SS    SS CC    CC FF      00    00         C
C          HH    HH FF       SSSSSS   CCCCCC  FF       000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  HFSCF PERFORMS AN ATOMIC SINGLE-DETERMINANT ITERATIVE SELF-         C
C  CONSISTENT FIELD PROCEDURE OVER THE USER-SPECIFIED HAMILTONIAN.     C
C  USES THE CLOSED-SHELL AVERAGE OF CONFIGURATION MODEL, WITH SUBSHELL C
C  OCCUPATIONS DETERMINED EITHER MANUALLY OR BY THE AUFBAU ROUTINE.    C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    ICNT - ATOMIC CENTER INDEX                                        C
C -------------------------------------------------------------------- C
C  DFNOTE: UNFINISHED -- BREIT INTERACTION STILL A MESS.               C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MNU=MKP+1,
     &                                              LWK=128*MBS,MIT=100)
C
      CHARACTER*1 ELLTERM,QSGN
      CHARACTER*2 ELMNT(120),ELNM
      CHARACTER*4 HMLTN
      CHARACTER*8 ZWRT,QWRT,EWRT
      CHARACTER*80 TITLE
C
      COMPLEX*16 C(MDM,MDM)
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
C    
      DIMENSION BMAT(MDM,MDM)
C
      DIMENSION QE(MKP),QA(MKP),NORB(MKP,12),NUMOCC(MKP)
      DIMENSION W1(2*MBS),W2(2*MBS),T(LWK)
      DIMENSION O1(2*MBS,2*MBS),H1(2*MBS,2*MBS),U1(2*MBS,2*MBS),
     &          O2(2*MBS,2*MBS),H2(2*MBS,2*MBS),U2(2*MBS,2*MBS)
      DIMENSION G11(2*MBS,2*MBS),G21(2*MBS,2*MBS),
     &          G12(2*MBS,2*MBS),G22(2*MBS,2*MBS)
      DIMENSION B11(2*MBS,2*MBS),B21(2*MBS,2*MBS),
     &          B12(2*MBS,2*MBS),B22(2*MBS,2*MBS)
      DIMENSION DENLL(MB2,2*MKP+1),DFNLL(MB2,2*MKP+1),
     &          DENSL(MB2,2*MKP+1),DFNSL(MB2,2*MKP+1),
     &          DENSS(MB2,2*MKP+1),DFNSS(MB2,2*MKP+1)
      DIMENSION DEN1(MB2,3),DEN2(MB2,3)
C
      COMMON/ANGL/BK(MNU,4),DK(MNU,4),HK(MNU,4),FK(MNU,4),GM(MNU,4),
     &            NUS(MNU),NUNUM,NUMIN,NUMAX
      COMMON/ATOM/ELMNT
      COMMON/BSQN/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/BSXP/EXLA(MBS),EXLB(MBS)
      COMMON/COEF/C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EMDR,EMXC,EUEH
      COMMON/FILL/NCNF(MCT,MKP,MKP+1),NLVL(MCT,MKP),IFILL(MCT)
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/OCPD/IOCPN(MDM),IOCCM0
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     EMPTY BMAT ARRAY
      DO I=1,NDIM
        DO J=1,NDIM
          BMAT(I,J) = 0.0D0
        ENDDO
      ENDDO
C
C     CONVERGENCE TOLERANCE VALUE
      IF(HMLTN.EQ.'NORL') THEN
        EEPS = 5.0D-12
      ELSE
        EEPS = 5.0D-13
      ENDIF
C
C     IMPORT ATOMIC CHARGE DETAILS
      IZN  = IZNUC(ICNT)
      ZCRG = ZNUC(ICNT)
      ELNM = ELMNT(IZN)
      ICRG = IQNUC(ICNT)
      MLQN = LMAX(ICNT)
C
C     CONVERT IZN AND ICRG TO STRINGS AND WRITE TITLE
      IF(IZN.LT.10) THEN
        WRITE(ZWRT,'(A,I1)') 'Z = ',IZN
      ELSEIF(IZN.LT.100) THEN
        WRITE(ZWRT,'(A,I2)') 'Z = ',IZN
      ELSE
        WRITE(ZWRT,'(A,I3)') 'Z = ',IZN
      ENDIF
C
      IF(ICRG.LT.10) THEN
        WRITE(QWRT,'(A,I2)') 'Q = ',ICRG
      ELSEIF(IZN.LT.100) THEN
        WRITE(QWRT,'(A,I3)') 'Q = ',ICRG
      ELSE
        WRITE(QWRT,'(A,I4)') 'Q = ',ICRG
      ENDIF
C
      IF(ICRG.GT.0) THEN
        QSGN = '+'
      ELSEIF(ICRG.LT.0) THEN
        QSGN = '-'
      ENDIF
C
      ICMD = IABS(ICRG)
      IF(ICRG.EQ.0) THEN
        WRITE(EWRT,'(A,A,A)') '(',TRIM(ELNM),')'
      ELSEIF(ICMD.EQ.1) THEN
        WRITE(EWRT,'(A,A,A,A,A)') '(',TRIM(ELNM),'^',QSGN,')'
      ELSEIF(ICMD.LT.10) THEN
        WRITE(EWRT,'(A,A,A,I1,A,A)') '(',TRIM(ELNM),'^',ICMD,QSGN,')'
      ELSEIF(ICMD.LT.100) THEN
        WRITE(EWRT,'(A,A,A,I2,A,A)') '(',TRIM(ELNM),'^',ICMD,QSGN,')'
      ELSE
        WRITE(EWRT,'(A,A,A,I2,A,A)') '(',TRIM(ELNM),'^',ICMD,QSGN,')'
      ENDIF
C
C     PRINT TITLE SUMMARY FOR THIS ATOM
20    FORMAT(17X,'Center',I3,':',3X,A,3X,A,3X,A)
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,20) ICNT,ZWRT,QWRT,EWRT
      WRITE(7,20) ICNT,ZWRT,QWRT,EWRT
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     READ OR GENERATE AUFBAU OCCUPANCY FOR THIS ATOM
      IF(IFILL(ICNT).EQ.0) THEN
        CALL AUFBAU(IZN,ICRG,NORB,NUMOCC,LMXCONF)
      ELSE
        LMXCONF = LMAX(ICNT)
        DO LQN=1,LMXCONF+1
          NUMOCC(LQN) = NLVL(ICNT,LQN)
          DO N=1,NLVL(ICNT,LQN)
            NORB(LQN,N) = NCNF(ICNT,LQN,N)
          ENDDO
        ENDDO
      ENDIF
C
C     IDENTIFY THE HIGHEST OCCUPIED SHELL
      NMAX = 1
      DO LQN=1,LMXCONF+1
        IF(NUMOCC(LQN).GT.NMAX) THEN
          NMAX = NUMOCC(LQN)
        ENDIF
      ENDDO
C
C     CHECK WHETHER THERE ARE SUFFICIENT BASIS FUNCTION TYPES
      IF(MLQN.LT.LMXCONF) THEN
        WRITE(6, *) 'In HFSCF0: insufficient angular types in basis.'
        WRITE(7, *) 'In HFSCF0: insufficient angular types in basis.'
        WRITE(6, *) 'MLQN = ',MLQN,' and LMXCONF = ',LMXCONF
        WRITE(7, *) 'MLQN = ',MLQN,' and LMXCONF = ',LMXCONF
        STOP
      ENDIF
C
C     PRINT ATOMIC CONFIGURATION
30    FORMAT(1X,A,2X,A,1X,'|',2X'NSHELL ',12(2X,I2))
31    FORMAT(1X,'LQN = ',I1,3X,I2,2X,'|'1X,' OCC(',A,'):',A,12(2X,I2))
C
      IF(IFILL(ICNT).EQ.0) THEN
        WRITE(6,30) 'Aufbau:','#fns',(N,N=1,12)
        WRITE(7,30) 'Aufbau:','#fns',(N,N=1,12)
      ELSE
        WRITE(6,30) 'Manual:','#fns',(N,N=1,12)
        WRITE(7,30) 'Manual:','#fns',(N,N=1,12)
      ENDIF
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO LQN=0,LMAX(ICNT)
        WRITE(6,31) LQN,NFUNCT(LQN+1,ICNT),ELLTERM(LQN),
     &               REPEAT(' ',LQN*4),(NORB(LQN+1,J),J=1,NUMOCC(LQN+1))
        WRITE(7,31) LQN,NFUNCT(LQN+1,ICNT),ELLTERM(LQN),
     &               REPEAT(' ',LQN*4),(NORB(LQN+1,J),J=1,NUMOCC(LQN+1))
      ENDDO
C
C     RESULTS FOR EACH ITERATION
40    FORMAT(1X,A,10X,A,10X,A,13X,A,11X,A)
41    FORMAT(I3,2X,F15.6,2X,F15.6,2X,F18.6,4X,ES12.5)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,40) 'It','E1 (au)','E2 (au)','ET (au)','Ratio'
      WRITE(7,40) 'It','E1 (au)','E2 (au)','ET (au)','Ratio'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     IMPORT NUCLEAR RADIUS FOR THIS CENTER
      PNUC = CNUC(ICNT)
C
C     INITIALISE A STORAGE BIN FOR PREVIOUS ATOMIC ENERGY
      EPRV = 0.0D0
C
C**********************************************************************C
C     ZERO-BODY PROBLEM: IZN-ICRG = 0. (NO ELECTRONS AROUND CENTER.)   C
C**********************************************************************C
C
      IF(IZN.EQ.ICRG) THEN
C
C       NO OCCUPYING ELECTRON -> NO EIGENVALUE NEEDED
C
        EH = 0.0D0
        EG = 0.0D0
        EB = 0.0D0
        ENEW = EH-EG-EB

        WRITE(6,41) 1,EH,EG,ENEW,1.0D0
        WRITE(7,41) 1,EH,EG,ENEW,1.0D0
C
        GOTO 1001
C
C**********************************************************************C
C     ONE-BODY PROBLEM: IZN-ICRG = 1. (NO COULOMB ENERGY.)             C
C**********************************************************************C
C
      ELSEIF(IZN-ICRG.EQ.1) THEN
C
C       IMPORT ORDERED ELECTRON OCCUPATION NUMBER
        IOCCML = IOCCM0
C
C       GROUND STATE OF SINGLY-OCCUPIED ATOM: LQNA = 0
        LQNA = 0
C
C       EFFECTIVE AND AVERAGE OCCUPATION NUMBERS FOR THIS LQNA ORBITAL
C       A CLOSED SUBSHELL (NSHELL,LQNA) CONTAINS NCLS ELECTRONS
        NCLS = 4*LQNA+2
C
C       CORRESPONDING RELATIVISTIC QUANTUM NUMBER
        KAPA2 =-LQNA-1
        IF(HMLTN.EQ.'NORL') THEN
          RK2A2 = DFLOAT(NCLS)
        ELSE
          RK2A2 = DFLOAT(2*IABS(KAPA2))
        ENDIF
C
C       IMPORT BASIS FUNCTION EXPONENTS
        NFUNA = NFUNCT(LQNA+1,ICNT)
        DO IBAS=1,NFUNA
          EXLA(IBAS) = EXPSET(IBAS,LQNA+1,ICNT)
        ENDDO
C
C       MATRIX DIMENSIONS FOR THIS LQNA BLOCK
        IF(HMLTN.EQ.'NORL') THEN
          NBLC = 0
        ELSE
          NBLC = NFUNA
        ENDIF
        NMAT = NFUNA+NBLC
C
C       GENERATE OVERLAP AND BARE DIRAC MATRICES
        CALL OVRLP0(O2,EXLA,     KAPA2,NFUNA)
        CALL ONEEL0(H2,EXLA,ZCRG,KAPA2,NFUNA)
C
C       ATOMIC UEHLING INTERACTION
        IF(HMLTN.EQ.'DHFQ') THEN
C
C         GENERATE UEHLING MATRIX ELEMENTS
          CALL UEHLING0(U2,EXLA,ZCRG,KAPA2,NFUNA)
C
C         ADD UEHLING MATRIX ELEMENTS CONTRIBUTIONS TO FOCK MATRIX
          DO IBAS=1,NMAT
            DO JBAS=1,NMAT
              H2(IBAS,JBAS) = H2(IBAS,JBAS) + U2(IBAS,JBAS)
            ENDDO
          ENDDO
C
        ENDIF
C
C       DIAGONALISE MATRIX (THIS NEEDS LAPACK LIBRARY)
        CALL DSYGV(1,'V','U',NMAT,H2,2*MBS,O2,2*MBS,W2,T,LWK,INFO)
        IF(INFO.NE.0) THEN
          WRITE(6, *) 'In HFSCF0: eigenvalue solver DSYGV failed.',INFO
          WRITE(7, *) 'In HFSCF0: eigenvalue solver DSYGV failed.',INFO
          STOP
        ENDIF
C
C       COEFFICIENT MATRIX ADDRESSES
        IL1 = LARGE(ICNT,2*LQNA+1,1)
        IL2 = LARGE(ICNT,2*LQNA+1,2)
        IS1 = IL1 + NSHIFT
        IS2 = IL2 + NSHIFT
C
C       EFFECTIVE OCCUPATION NUMBER
        QF = 1.0D0/DSQRT(2.0D0)
C
C       COPY INTO MASTER COEFFICIENT LIST
        DO IBAS=1,NFUNA
C
C         LARGE COMPONENT OF KRAMERS PAIR
          CL = QF*H2(IBAS      ,NBLC+1)
          C(IL1+IBAS,IOCCML+NSHIFT+1) = DCMPLX(CL,0.0D0)
          C(IL2+IBAS,IOCCML+NSHIFT+2) = DCMPLX(CL,0.0D0)
C
C         SMALL COMPONENT OF KRAMERS PAIR
          IF(HMLTN.NE.'NORL') THEN
            CS = QF*H2(IBAS+NBLC,NBLC+1)
            C(IS1+IBAS,IOCCML+NSHIFT+1) = DCMPLX(CS,0.0D0)
            C(IS2+IBAS,IOCCML+NSHIFT+2) = DCMPLX(CS,0.0D0)
          ENDIF
C
        ENDDO
C
C       STORE LOWEST ENERGY EIGENVALUES TO MASTER LIST
        EIGEN(IOCCML+NSHIFT+1) = W2(NBLC+1)
        EIGEN(IOCCML+NSHIFT+2) = W2(NBLC+1)
C
C       INCREASE FOCK ADDRESS OF OCCUPIED ORBITALS (PAIR AT A TIME)
        IOCCML = IOCCML+2
C
C       ONE- AND TWO-BODY ENERGIES
        EH = W2(NBLC+1)
C
C       UEHLING INTERACTION ENERGIES FOR OCCUPIED ELECTRONS
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNA
            M = M+1
C
C           SMALL-COMPONENT ADDRESSES
            KBAS = IBAS+NBLC
            LBAS = JBAS+NBLC
C
C           ADD CONTRIBUTIONS TO UEHLING ENERGY
            EU = EU + U2(IBAS,JBAS)*H2(IBAS,NBLC+1)*H2(JBAS,NBLC+1)
            IF(HMLTN.NE.'NORL') THEN
              EU = EU + U2(KBAS,LBAS)*H2(KBAS,NBLC+1)*H2(LBAS,NBLC+1)
            ENDIF
C
          ENDDO
        ENDDO
C        
        EG = 0.0D0
        EB = 0.0D0
        ENEW = EH
C
C       WRITE RESULT
        WRITE(6,41) 1,EH,0.0D0,ENEW,1.0D0
        WRITE(7,41) 1,EH,0.0D0,ENEW,1.0D0
C
C       UPDATE FOCK LABEL FOR OCCUPATION COUNTER
        IOCCM0 = IOCCML
C
C       EXIT TO CONVERGENCE
        GOTO 1001
C
      ENDIF
C
C**********************************************************************C
C     TWO-BODY PROBLEM: INTERACTING ELECTRONS. (TREAT WITH SCF.)       C
C -------------------------------------------------------------------- C
C     ENTER ITERATIVE SELF-CONSISTENT FIELD PROCEDURE (USE INDEX 1000) C
C**********************************************************************C
C
      DO 1000 ITER=1,MIT

C       INITIALISE ONE-BODY AND TWO-BODY ENERGY COUNTERS
        EH = 0.0D0
        EU = 0.0D0
        EG = 0.0D0
        EB = 0.0D0
C
C       INITIALISE ELECTRON OCCUPATION COUNTER
        IOCCML = IOCCM0
C
C**********************************************************************C
C     FIRST LOOP: OVER BASIS FUNCTIONS I,J (USE INDEX 100)             C
C**********************************************************************C
C
C     LOOP OVER ALL OCCUPIED LQN VALUES
      DO 100 LQNA=0,LMXCONF
C
C     RECORD LQNA VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
      NFUNA = NFUNCT(LQNA+1,ICNT)
      DO IBAS=1,NFUNA
        EXLA(IBAS) = EXPSET(IBAS,LQNA+1,ICNT)
      ENDDO
C
C     MATRIX DIMENSIONS FOR THIS LQNA BLOCK
      IF(HMLTN.EQ.'NORL') THEN
        NBLC = 0
      ELSE
        NBLC = NFUNA
      ENDIF
      NMAT = NFUNA+NBLC
C
C     EFFECTIVE AND AVERAGE OCCUPATION NUMBERS FOR THIS LQNA ORBITAL
C     A CLOSED SUBSHELL (NSHELL,LQNA) CONTAINS NCLS ELECTRONS
      NCLS = 4*LQNA+2
C
C     FOR EACH OCCUPIED NSHELL OF THIS LQNA CLASS
      DO IOCC=1,NUMOCC(LQNA+1)
C
C       NUMBER OF CHARGES IN THIS SUBSHELL (NSHELL,LQNA)
        NQ = NORB(LQNA+1,IOCC)
C
C       IF SUBSHELL IS CLOSED THERE IS NO FRACTIONAL OCCUPANCY
        IF(NQ.EQ.NCLS) THEN
          QE(IOCC) = 1.0D0
C       IF SUBSHELL IS OPEN, CONSTRUCT FRACTION (GRANT 6.6.24)
        ELSE
          QE(IOCC) = DFLOAT(NQ-1)/DFLOAT(NCLS-1)
        ENDIF
C
C       ACTUAL FRACTIONAL SUBSHELL OCCUPANCY
        IF(NQ.GT.0) THEN
          QA(IOCC) = DFLOAT(NQ)/DFLOAT(NCLS)
        ELSE
          QA(IOCC) = 0.0D0
        ENDIF
C
      ENDDO
C
C >>> POSITIVE KAPPA(A) CHOICE (APPLIES ONLY FOR LQNA > 0)
      IF(LQNA.EQ.0.OR.HMLTN.EQ.'NORL') GOTO 130

      KAPA1 = LQNA
      RK2A1 = DFLOAT(2*IABS(KAPA1))
C
C     GENERATE OVERLAP AND BARE DIRAC MATRICES
      CALL OVRLP0(O1,EXLA,     KAPA1,NFUNA)
      CALL ONEEL0(H1,EXLA,ZCRG,KAPA1,NFUNA)
C
C     ATOMIC UEHLING INTERACTION
      IF(HMLTN.EQ.'DHFQ') THEN
C
C       GENERATE UEHLING MATRIX ELEMENTS
        CALL UEHLING0(U1,EXLA,ZCRG,KAPA1,NFUNA)
C
C       ADD UEHLING MATRIX ELEMENTS CONTRIBUTIONS TO FOCK MATRIX
        DO IBAS=1,NMAT
          DO JBAS=1,NMAT
            H1(IBAS,JBAS) = H1(IBAS,JBAS) + U1(IBAS,JBAS)
          ENDDO
        ENDDO
C
      ENDIF
C
130   CONTINUE
C
C >>> NEGATIVE KAPPA(A) CHOICE (APPLIES TO ALL LQNA VALUES)
      KAPA2 =-LQNA-1
      IF(HMLTN.EQ.'NORL') THEN
        RK2A2 = DFLOAT(NCLS)
      ELSE
        RK2A2 = DFLOAT(2*IABS(KAPA2))
      ENDIF
C
C     GENERATE OVERLAP AND BARE DIRAC MATRICES
      CALL OVRLP0(O2,EXLA,     KAPA2,NFUNA)
      CALL ONEEL0(H2,EXLA,ZCRG,KAPA2,NFUNA)
C
C     ATOMIC UEHLING INTERACTION
      IF(HMLTN.EQ.'DHFQ') THEN
C
C       GENERATE UEHLING MATRIX ELEMENTS
        CALL UEHLING0(U2,EXLA,ZCRG,KAPA2,NFUNA)
C
C       ADD UEHLING MATRIX ELEMENTS CONTRIBUTIONS TO FOCK MATRIX
        DO IBAS=1,NMAT
          DO JBAS=1,NMAT
            H2(IBAS,JBAS) = H2(IBAS,JBAS) + U2(IBAS,JBAS)
          ENDDO
        ENDDO
C
      ENDIF
C
C     SKIP SCF INTERACTIONS IN FIRST ITERATION
      IF(ITER.EQ.1) GOTO 150
C
C     INITIALISE RELEVANT COUNTERS AND ARRAYS
      RK2B1 = 0.0D0
      RK2B2 = 0.0D0
      DO IBAS=1,NMAT
        DO JBAS=1,NMAT
          G11(IBAS,JBAS) = 0.0D0
          G12(IBAS,JBAS) = 0.0D0
          G21(IBAS,JBAS) = 0.0D0
          G22(IBAS,JBAS) = 0.0D0
          B11(IBAS,JBAS) = 0.0D0
          B12(IBAS,JBAS) = 0.0D0
          B21(IBAS,JBAS) = 0.0D0
          B22(IBAS,JBAS) = 0.0D0
        ENDDO
      ENDDO
C
C**********************************************************************C
C     SECOND LOOP: OVER BASIS FUNCTIONS K,L (USE INDEX 200)            C
C**********************************************************************C
C
C     LOOP OVER ALL OCCUPIED LQN VALUES
      DO 200 LQNB=0,LMXCONF
C
C     RECORD LQNB VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
      NFUNB = NFUNCT(LQNB+1,ICNT)
      DO JBAS=1,NFUNB
        EXLB(JBAS) = EXPSET(JBAS,LQNB+1,ICNT)
      ENDDO
C
C     NUMBER OF BASIS FUNCTION OVERLAPS IN THIS BLOCK
      MAXM = NFUNB*NFUNB
C
C**********************************************************************C
C     GENERATE ATOMIC FOCK MATRIX                                      C
C**********************************************************************C
C
C     EVALUATE CLOSED-SHELL ELECTRON REPULSION ANGULAR INTEGRALS
      CALL ANGCOUL
C
C >>> POSITIVE KAPPA(B) CHOICE (APPLIES ONLY FOR LQNB > 0)
      IF(LQNB.EQ.0.OR.HMLTN.EQ.'NORL') GOTO 230
C
C     KAPPA(B) VALUE AND DEGENERACY
      KAPB1 = LQNB
      RK2B1 = DFLOAT(2*IABS(KAPB1))
C
C     RELEVANT DENSITY MATRIX DEPENDS ON LQNA AND LQNB
      IF(LQNA.EQ.LQNB) THEN
        DO M=1,MAXM
          DEN1(M,1) = DENLL(M,2*LQNB  )
          DEN1(M,2) = DENSL(M,2*LQNB  )
          DEN1(M,3) = DENSS(M,2*LQNB  )
        ENDDO        
      ELSEIF(LQNA.NE.LQNB) THEN
        DO M=1,MAXM
          DEN1(M,1) = DFNLL(M,2*LQNB  )
          DEN1(M,2) = DFNSL(M,2*LQNB  )
          DEN1(M,3) = DFNSS(M,2*LQNB  )
        ENDDO        
      ENDIF
C
230   CONTINUE
C
C >>> NEGATIVE KAPPA(B) CHOICE (APPLIES TO ALL LQNB VALUES)
C
C     KAPPA(B) VALUE AND DEGENERACY
      KAPB2 =-LQNB-1
      IF(HMLTN.EQ.'NORL') THEN
        RK2B2 = DFLOAT(4*LQNB+2)
      ELSE
        RK2B2 = DFLOAT(2*IABS(KAPB2))
      ENDIF
C
C     RELEVANT DENSITY MATRIX DEPENDS ON LQNA AND LQNB
      IF(LQNA.EQ.LQNB) THEN
        DO M=1,MAXM
          DEN2(M,1) = DENLL(M,2*LQNB+1)
          IF(HMLTN.NE.'NORL') THEN
            DEN2(M,2) = DENSL(M,2*LQNB+1)
            DEN2(M,3) = DENSS(M,2*LQNB+1)
          ENDIF
        ENDDO        
      ELSEIF(LQNA.NE.LQNB) THEN
        DO M=1,MAXM
          DEN2(M,1) = DFNLL(M,2*LQNB+1)
          IF(HMLTN.NE.'NORL') THEN
            DEN2(M,2) = DFNSL(M,2*LQNB+1)
            DEN2(M,3) = DFNSS(M,2*LQNB+1)
          ENDIF
        ENDDO
      ENDIF
C
C     GENERATE THE MEAN-FIELD ATOMIC COULOMB MATRIX OVER DENSITIES
      CALL COULOMB0(G11,G21,G12,G22,DEN1,DEN2)
C
C     ADD COULOMB MATRIX ELEMENTS CONTRIBUTIONS TO FOCK MATRIX
      DO IBAS=1,NMAT
        DO JBAS=1,NMAT
C
          IF(HMLTN.EQ.'NORL') THEN
C         NON-RELATIVISTIC COULOMB MATRIX
C
            H2(IBAS,JBAS) = H2(IBAS,JBAS) + RK2B2*G22(IBAS,JBAS)
C
          ELSE
C         RELATIVISTIC COULOMB MATRIX
C
            H1(IBAS,JBAS) = H1(IBAS,JBAS) + RK2B1*G11(IBAS,JBAS)
     &                                    + RK2B2*G12(IBAS,JBAS)
            H2(IBAS,JBAS) = H2(IBAS,JBAS) + RK2B1*G21(IBAS,JBAS)
     &                                    + RK2B2*G22(IBAS,JBAS)
C
          ENDIF
C
        ENDDO
      ENDDO
      
C      WRITE(*,*) '****',H2(1,1),RK2B1,G21(1,1),RK2B2,G22(1,1)
C
C     TWO-BODY EIGENVALUE ENERGIES FOR OCCUPIED ELECTRONS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNA
          M = M+1
C
          IF(HMLTN.EQ.'NORL') THEN
C         NON-RELATIVISTIC COULOMB ENERGY
C
            EG = EG + RK2A2*RK2B2*G22(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
C
          ELSE
C         RELATIVISTIC  COULOMB ENERGY
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NFUNA
            LBAS = JBAS+NFUNA
C
C           LL BLOCK
            EG = EG
     &         +       RK2A1*RK2B1*G11(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &         +       RK2A1*RK2B2*G12(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &         +       RK2A2*RK2B1*G21(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &         +       RK2A2*RK2B2*G22(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
C
C           SL BLOCK
            EG = EG
     &         + 2.0D0*RK2A1*RK2B1*G11(KBAS,JBAS)*DFNSL(M,2*LQNA  )
     &         + 2.0D0*RK2A1*RK2B2*G12(KBAS,JBAS)*DFNSL(M,2*LQNA  )
     &         + 2.0D0*RK2A2*RK2B1*G21(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
     &         + 2.0D0*RK2A2*RK2B2*G22(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
C
C           SS BLOCK
            EG = EG
     &         +       RK2A1*RK2B1*G11(KBAS,LBAS)*DFNSS(M,2*LQNA  )
     &         +       RK2A1*RK2B2*G12(KBAS,LBAS)*DFNSS(M,2*LQNA  )
     &         +       RK2A2*RK2B1*G21(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
     &         +       RK2A2*RK2B2*G22(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
          ENDIF
C
        ENDDO
      ENDDO
C
C     GENERATE THE MEAN-FIELD ATOMIC BREIT MATRIX
      IF(HMLTN.NE.'DHFB'.AND.HMLTN.NE.'DHFQ') GOTO 250
      IF(ITER.EQ.15) THEN
C
C     EVALUATE CLOSED-SHELL BREIT INTERACTION ANGULAR INTEGRALS
      CALL ANGBREIT
C
C     GENERATE THE MEAN-FIELD ATOMIC BREIT MATRIX OVER DENSITIES
      CALL BREIT0(B11,B21,B12,B22,DEN1,DEN2)
C
      GOTO 251
C     ADD TWO-PARTICLE CONTRIBUTIONS TO FOCK MATRIX
      DO IBAS=1,2*NFUNA
        DO JBAS=1,2*NFUNA
          H1(IBAS,JBAS) = H1(IBAS,JBAS) + RK2B1*B11(IBAS,JBAS)
     &                                  + RK2B2*B12(IBAS,JBAS)
          H2(IBAS,JBAS) = H2(IBAS,JBAS) + RK2B1*B21(IBAS,JBAS)
     &                                  + RK2B2*B22(IBAS,JBAS)
        ENDDO
      ENDDO
251   CONTINUE
C
C     TWO-BODY EIGENVALUE ENERGIES FOR OCCUPIED ELECTRONS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNA
          M = M+1
C
C         SMALL COMPONENT BLOCK ADDRESSES
          KBAS = IBAS+NFUNA
          LBAS = JBAS+NFUNA
C
C         LL BLOCK
          EB = EB
     &       +       RK2A1*RK2B1*B11(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &       +       RK2A1*RK2B2*B12(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &       +       RK2A2*RK2B1*B21(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &       +       RK2A2*RK2B2*B22(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
C
C         SL BLOCK
          EB = EB
     &       + 2.0D0*RK2A1*RK2B1*B11(KBAS,JBAS)*DFNSL(M,2*LQNA  )
     &       + 2.0D0*RK2A1*RK2B2*B12(KBAS,JBAS)*DFNSL(M,2*LQNA  )
     &       + 2.0D0*RK2A2*RK2B1*B21(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
     &       + 2.0D0*RK2A2*RK2B2*B22(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
C
C         SS BLOCK
          EB = EB
     &       +       RK2A1*RK2B1*B11(KBAS,LBAS)*DFNSS(M,2*LQNA  )
     &       +       RK2A1*RK2B2*B12(KBAS,LBAS)*DFNSS(M,2*LQNA  )
     &       +       RK2A2*RK2B1*B21(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
     &       +       RK2A2*RK2B2*B22(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
        ENDDO
      ENDDO
C
      WRITE(*,*) 'BENERGY = ',EB
      WRITE(*,*) 'CORRECT = ',6.3774D-5
      WRITE(*,*) '|ERROR| = ',DABS(EB-6.3774D-5)
      EB = 0.0D0
C
C     PUT BREIT MATRIX COMPONENTS INTO A BIGGER MATRIX
C
C     IGNORE ANYTHING THAT ISN'T AN s-TYPE OVERLAP TO START WITH...
      GOTO 160
C
      IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 161
C     CASE KQNA>0, KQNB>0
      IL1 = LARGE(ICNT,2*LQNA  ,1)
      IL2 = LARGE(ICNT,2*LQNA  ,2)
      IS1 = IL1 + NSHIFT
      IS2 = IL2 + NSHIFT
C
      JL1 = LARGE(ICNT,2*LQNB  ,1)
      JL2 = LARGE(ICNT,2*LQNB  ,2)
      JS1 = JL1 + NSHIFT
      JS2 = JL2 + NSHIFT
C
      RKK = RK2A1*RK2B1
C
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNA
C
C         SMALL COMPONENT BLOCK ADDRESSES
          KBAS = IBAS+NFUNA
          LBAS = JBAS+NFUNA
C
          BMAT(IL1+IBAS,JL1+JBAS) = B11(IBAS,JBAS)
          BMAT(IL2+IBAS,JL2+JBAS) = B11(IBAS,JBAS)

          BMAT(IL1+IBAS,JS1+JBAS) = B11(IBAS,LBAS)
          BMAT(IL2+IBAS,JS2+JBAS) = B11(IBAS,LBAS)

          BMAT(IS1+IBAS,JL1+JBAS) = B11(KBAS,JBAS)
          BMAT(IS2+IBAS,JL2+JBAS) = B11(KBAS,JBAS)

          BMAT(IS1+IBAS,JS1+JBAS) = B11(KBAS,LBAS)
          BMAT(IS2+IBAS,JS2+JBAS) = B11(KBAS,LBAS)
C
        ENDDO
      ENDDO
161   CONTINUE
C
      IF(LQNA.EQ.0) GOTO 162
C     CASE KQNA>0, KQNB<0
      IL1 = LARGE(ICNT,2*LQNA  ,1)
      IL2 = LARGE(ICNT,2*LQNA  ,2)
      IS1 = IL1 + NSHIFT
      IS2 = IL2 + NSHIFT
C
      JL1 = LARGE(ICNT,2*LQNB+1,1)
      JL2 = LARGE(ICNT,2*LQNB+1,2)
      JS1 = JL1 + NSHIFT
      JS2 = JL2 + NSHIFT
C
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNA
C
C         SMALL COMPONENT BLOCK ADDRESSES
          KBAS = IBAS+NFUNA
          LBAS = JBAS+NFUNA
C
          BMAT(IL1+IBAS,JL1+JBAS) = B12(IBAS,JBAS)
          BMAT(IL2+IBAS,JL2+JBAS) = B12(IBAS,JBAS)

          BMAT(IL1+IBAS,JS1+JBAS) = B12(IBAS,LBAS)
          BMAT(IL2+IBAS,JS2+JBAS) = B12(IBAS,LBAS)

          BMAT(IS1+IBAS,JL1+JBAS) = B12(KBAS,JBAS)
          BMAT(IS2+IBAS,JL2+JBAS) = B12(KBAS,JBAS)

          BMAT(IS1+IBAS,JS1+JBAS) = B12(KBAS,LBAS)
          BMAT(IS2+IBAS,JS2+JBAS) = B12(KBAS,LBAS)
C
        ENDDO
      ENDDO
162   CONTINUE
C
      IF(LQNB.EQ.0) GOTO 163
C     CASE KQNA<0, KQNB>0
      IL1 = LARGE(ICNT,2*LQNA+1,1)
      IL2 = LARGE(ICNT,2*LQNA+1,2)
      IS1 = IL1 + NSHIFT
      IS2 = IL2 + NSHIFT
C
      JL1 = LARGE(ICNT,2*LQNB  ,1)
      JL2 = LARGE(ICNT,2*LQNB  ,2)
      JS1 = JL1 + NSHIFT
      JS2 = JL2 + NSHIFT
C
      RKK = RK2A2*RK2B1
C
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNA
C
C         SMALL COMPONENT BLOCK ADDRESSES
          KBAS = IBAS+NFUNA
          LBAS = JBAS+NFUNA
C
          BMAT(IL1+IBAS,JL1+JBAS) = B21(IBAS,JBAS)
          BMAT(IL2+IBAS,JL2+JBAS) = B21(IBAS,JBAS)

          BMAT(IL1+IBAS,JS1+JBAS) = B21(IBAS,LBAS)
          BMAT(IL2+IBAS,JS2+JBAS) = B21(IBAS,LBAS)

          BMAT(IS1+IBAS,JL1+JBAS) = B21(KBAS,JBAS)
          BMAT(IS2+IBAS,JL2+JBAS) = B21(KBAS,JBAS)

          BMAT(IS1+IBAS,JS1+JBAS) = B21(KBAS,LBAS)
          BMAT(IS2+IBAS,JS2+JBAS) = B21(KBAS,LBAS)
C
        ENDDO
      ENDDO
163   CONTINUE

C     s-TYPE OVERLAPS...
160   CONTINUE
C
C     CASE KQNA<0, KQNB<0
      IF(LQNA.NE.0.OR.LQNB.NE.0) THEN
        WRITE(*,*) 'OH DEAR'
        GOTO 169
      ENDIF
      IL1 = LARGE(ICNT,2*LQNA+1,1)
      IL2 = LARGE(ICNT,2*LQNA+1,2)
      IS1 = IL1 + NSHIFT
      IS2 = IL2 + NSHIFT
C
      JL1 = LARGE(ICNT,2*LQNB+1,1)
      JL2 = LARGE(ICNT,2*LQNB+1,2)
      JS1 = JL1 + NSHIFT
      JS2 = JL2 + NSHIFT
C
      RKK = RK2A2*RK2B2
C
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNA
C
C         SMALL COMPONENT BLOCK ADDRESSES
          KBAS = IBAS+NFUNA
          LBAS = JBAS+NFUNA
C
          BMAT(IL1+IBAS,JL1+JBAS) = B22(IBAS,JBAS)
          BMAT(IL2+IBAS,JL2+JBAS) = B22(IBAS,JBAS)
C
          BMAT(IL1+IBAS,JS1+JBAS) = B22(IBAS,LBAS)
          BMAT(IL2+IBAS,JS2+JBAS) = B22(IBAS,LBAS)
C
          BMAT(IS1+IBAS,JL1+JBAS) = B22(KBAS,JBAS)
          BMAT(IS2+IBAS,JL2+JBAS) = B22(KBAS,JBAS)
C
          BMAT(IS1+IBAS,JS1+JBAS) = B22(KBAS,LBAS)
          BMAT(IS2+IBAS,JS2+JBAS) = B22(KBAS,LBAS)
          
          IF(IBAS.EQ.JBAS) THEN
            WRITE(*,*) IBAS,B22(KBAS,LBAS)
          ENDIF
C
        ENDDO
      ENDDO
169   CONTINUE
C
C     ITERATION 15 END IF
      ENDIF
C
250   CONTINUE
C
C     END LOOP OVER LQNS FOR ORBITAL B
200   CONTINUE
C
C       UEHLING INTERACTION ENERGIES FOR OCCUPIED ELECTRONS
        IF(HMLTN.EQ.'DHFQ') THEN
C
          M = 0
          DO IBAS=1,NFUNA
            DO JBAS=1,NFUNA
              M = M+1
C
              IF(HMLTN.EQ.'NORL') THEN
C             NON-RELATIVISTIC UEHLING ENERGY
C
                EU = EU + RK2A2*U2(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
C
              ELSE
C             RELATIVISTIC  UEHLING ENERGY
C
C               SMALL COMPONENT BLOCK ADDRESSES
                KBAS = IBAS+NFUNA
                LBAS = JBAS+NFUNA
C
C               LL BLOCK
                EU = EU + RK2A1*U1(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &                  + RK2A2*U2(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
C
C               SS BLOCK
                EU = EU + RK2A1*U1(KBAS,LBAS)*DFNSS(M,2*LQNA  )
     &                  + RK2A2*U2(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
              ENDIF
C
            ENDDO
          ENDDO
C
        ENDIF
C
C     FINISH GENERATING SCF CONTRIBUTIONS
150   CONTINUE
C
C     FINISHED CALCULATING OVERLAP COMBINATIONS BETWEEN THIS LQNA
C     VALUE AND ALL POSSIBLE LQNB VALUES
C
C**********************************************************************C
C     MATRIX DIAGONALISATION AND COEFFICIENT MATRIX UPDATES            C
C**********************************************************************C
C
C >>> POSITIVE KAPPA(A) CHOICE (APPLIES ONLY FOR LQNA > 0)
      IF(LQNA.EQ.0.OR.HMLTN.EQ.'NORL') GOTO 140
C
C     DIAGONALISE FOCK MATRIX (THIS NEEDS LAPACK LIBRARY)
      CALL DSYGV(1,'V','U',NMAT,H1,2*MBS,O1,2*MBS,W1,T,LWK,INFO)
      IF(INFO.NE.0) THEN
        WRITE(6, *) 'In HFSCF0: eigenvalue solver DSYGV failed.',INFO
        WRITE(7, *) 'In HFSCF0: eigenvalue solver DSYGV failed.',INFO
        STOP
      ENDIF
C
C     ATOMIC SELECTION RULE: ORTHOGONALITY IN BLOCKS OF KQN -> KA = KB
C
C     BEGIN LOOP OVER MQNA VALUES
      DO IMVAL=1,IABS(KAPA1)
C
C       COEFFICIENT MATRIX ADDRESSES
        IL1 = LARGE(ICNT,2*LQNA  ,IMVAL*2-1)
        IL2 = LARGE(ICNT,2*LQNA  ,IMVAL*2  )
        IS1 = IL1 + NSHIFT
        IS2 = IL2 + NSHIFT
C
C       COPY INTO MASTER COEFFICIENT LIST IF QA IS POSITIVE
        DO IOCC=1,NUMOCC(LQNA+1)
C
C         EFFECTIVE OCCUPATION NUMBER
          QF = DSQRT(QA(IOCC))
C
C         COPY INTO MASTER COEFFICIENT LIST
          DO IBAS=1,NFUNA
C
C           LARGE COMPONENT OF KRAMERS PAIR
            CL = QF*H1(IBAS      ,NFUNA+IOCC)
            C(IL1+IBAS,IOCCML+NSHIFT+1) = DCMPLX(CL,0.0D0)
            C(IL2+IBAS,IOCCML+NSHIFT+2) = DCMPLX(CL,0.0D0)
C
C           SMALL COMPONENT OF KRAMERS PAIR
            CS = QF*H1(IBAS+NFUNA,NFUNA+IOCC)
            C(IS1+IBAS,IOCCML+NSHIFT+1) = DCMPLX(CS,0.0D0)
            C(IS2+IBAS,IOCCML+NSHIFT+2) = DCMPLX(CS,0.0D0)
C
          ENDDO
C
C         STORE LOWEST ENERGY EIGENVALUES TO MASTER LIST
          EIGEN(IOCCML+NSHIFT+1) = W1(NFUNA+IOCC)
          EIGEN(IOCCML+NSHIFT+2) = W1(NFUNA+IOCC)
C
C         INCREASE FOCK ADDRESS OF OCCUPIED ORBITALS (PAIR AT A TIME)
          IOCCML = IOCCML+2
C
        ENDDO
      ENDDO
C
C
C     BUILD ATOMIC CHARGE DENSITY MATRIX FOR THIS KQNA BLOCK
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNA
          M = M+1
C
C         INITIALISE ATOMIC DENSITY LISTS FOR THIS BLOCK
          DENLL(M,2*LQNA  ) = 0.0D0
          DENSL(M,2*LQNA  ) = 0.0D0
          DENSS(M,2*LQNA  ) = 0.0D0
C
          DFNLL(M,2*LQNA  ) = 0.0D0
          DFNSL(M,2*LQNA  ) = 0.0D0
          DFNSS(M,2*LQNA  ) = 0.0D0
C
C         LOOP OVER ALL OCCUPIED SHELLS OF THIS KQN TYPE
          DO IOCC=1,NUMOCC(LQNA+1)
C
C           DENSITY OVERLAPS FROM EIGENVECTOR PRODUCTS
            DLL = H1(IBAS     ,NBLC+IOCC)*H1(JBAS     ,NBLC+IOCC)
            DSL = H1(IBAS+NBLC,NBLC+IOCC)*H1(JBAS     ,NBLC+IOCC)
            DSS = H1(IBAS+NBLC,NBLC+IOCC)*H1(JBAS+NBLC,NBLC+IOCC)
C
C           ADD DENSITY CONTRIBUTIONS TO ATOMIC LIST
            DENLL(M,2*LQNA  ) = DENLL(M,2*LQNA  ) + QE(IOCC)*DLL
            DENSL(M,2*LQNA  ) = DENSL(M,2*LQNA  ) + QE(IOCC)*DSL
            DENSS(M,2*LQNA  ) = DENSS(M,2*LQNA  ) + QE(IOCC)*DSS
C
            DFNLL(M,2*LQNA  ) = DFNLL(M,2*LQNA  ) + QA(IOCC)*DLL
            DFNSL(M,2*LQNA  ) = DFNSL(M,2*LQNA  ) + QA(IOCC)*DSL
            DFNSS(M,2*LQNA  ) = DFNSS(M,2*LQNA  ) + QA(IOCC)*DSS
C
          ENDDO
C
        ENDDO
      ENDDO
C
C     ONE-BODY EIGENVALUE ENERGIES FOR OCCUPIED ELECTRONS
      DO IOCC=1,NUMOCC(LQNA+1)
        EH = EH + QA(IOCC)*RK2A1*W1(NBLC+IOCC)
      ENDDO
C
140   CONTINUE
C
C
C >>> NEGATIVE KAPPA(A) CHOICE (APPLIES TO ALL LQNA VALUES)
C
C     DIAGONALISE FOCK MATRIX (THIS NEEDS LAPACK LIBRARY)
      CALL DSYGV(1,'V','U',NMAT,H2,2*MBS,O2,2*MBS,W2,T,LWK,INFO)
      IF(INFO.NE.0) THEN
        WRITE(6, *) 'In HFSCF0: eigenvalue solver DSYGV failed.',INFO
        WRITE(7, *) 'In HFSCF0: eigenvalue solver DSYGV failed.',INFO
        STOP
      ENDIF
C
C     ATOMIC SELECTION RULE: ORTHOGONALITY IN BLOCKS OF KQN -> KA = KB
C
C     BEGIN LOOP OVER MQNA VALUES
      DO IMVAL=1,IABS(KAPA2)
C
C       COEFFICIENT MATRIX ADDRESSES
        IL1 = LARGE(ICNT,2*LQNA+1,IMVAL*2-1)
        IL2 = LARGE(ICNT,2*LQNA+1,IMVAL*2  )
        IS1 = IL1 + NSHIFT
        IS2 = IL2 + NSHIFT
C
C       LOOP OVER ALL OCCUPIED PRINCIPAL SHELLS N FOR THIS KQNA
        DO IOCC=1,NUMOCC(LQNA+1)
C
C         EFFECTIVE OCCUPATION NUMBER
          QF = DSQRT(QA(IOCC))
C
C         COPY INTO MASTER COEFFICIENT LIST
          DO IBAS=1,NFUNA
C
C           LARGE COMPONENT OF KRAMERS PAIR
            CL = QF*H2(IBAS     ,NBLC+IOCC)
            C(IL1+IBAS,IOCCML+NSHIFT+1) = DCMPLX(CL,0.0D0)
            C(IL2+IBAS,IOCCML+NSHIFT+2) = DCMPLX(CL,0.0D0)
C
C           SMALL COMPONENT OF KRAMERS PAIR
            IF(HMLTN.NE.'NORL') THEN
              CS = QF*H2(IBAS+NBLC,NBLC+IOCC)
              C(IS1+IBAS,IOCCML+NSHIFT+1) = DCMPLX(CS,0.0D0)
              C(IS2+IBAS,IOCCML+NSHIFT+2) = DCMPLX(CS,0.0D0)
            ENDIF
C
          ENDDO
C
C         STORE LOWEST ENERGY EIGENVALUES TO MASTER LIST
          EIGEN(IOCCML+NSHIFT+1) = W2(NBLC+IOCC)
          EIGEN(IOCCML+NSHIFT+2) = W2(NBLC+IOCC)
C
C         INCREASE FOCK ADDRESS OF OCCUPIED ORBITALS (PAIR AT A TIME)
          IOCCML = IOCCML+2
C
        ENDDO
      ENDDO
C
C     NON-RELATIVISTIC SPECIAL CASE: ALSO FILL IN THE +KQNA BLOCK
      IF(HMLTN.EQ.'NORL'.AND.LQNA.GE.1) THEN
C
C       BEGIN LOOP OVER MQNA VALUES IN +KQNA BLOCK
        DO IMVAL=1,IABS(KVALS(2*LQNA  ,ICNT))
C
C         COEFFICIENT MATRIX ADDRESSES
          IL1 = LARGE(ICNT,2*LQNA  ,IMVAL*2-1)
          IL2 = LARGE(ICNT,2*LQNA  ,IMVAL*2  )
C
C         LOOP OVER ALL OCCUPIED PRINCIPAL SHELLS N FOR THIS KQNA
          DO IOCC=1,NUMOCC(LQNA+1)
C
C           EFFECTIVE OCCUPATION NUMBER
            QF = DSQRT(QA(IOCC))
C
C           COPY INTO MASTER COEFFICIENT LIST
            DO IBAS=1,NFUNA
C
C             LARGE COMPONENT OF KRAMERS PAIR
              CL = QF*H2(IBAS     ,NBLC+IOCC)
              C(IL1+IBAS,IOCCML+NSHIFT+1) = DCMPLX(CL,0.0D0)
              C(IL2+IBAS,IOCCML+NSHIFT+2) = DCMPLX(CL,0.0D0)
C
            ENDDO
C
C         STORE LOWEST ENERGY EIGENVALUES TO MASTER LIST
          EIGEN(IOCCML+NSHIFT+1) = W2(NBLC+IOCC)
          EIGEN(IOCCML+NSHIFT+2) = W2(NBLC+IOCC)
C
C         INCREASE FOCK ADDRESS OF OCCUPIED ORBITALS (PAIR AT A TIME)
          IOCCML = IOCCML+2
C
          ENDDO
C
        ENDDO    
      ENDIF
C
C     BUILD ATOMIC CHARGE DENSITY MATRIX FOR THIS KQNA BLOCK
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNA
          M = M+1
C
C         INITIALISE ATOMIC DENSITY LISTS FOR THIS BLOCK
          DENLL(M,2*LQNA+1) = 0.0D0
          DFNLL(M,2*LQNA+1) = 0.0D0
C
          IF(HMLTN.NE.'NORL') THEN
            DENSL(M,2*LQNA+1) = 0.0D0
            DFNSL(M,2*LQNA+1) = 0.0D0
C
            DENSS(M,2*LQNA+1) = 0.0D0
            DFNSS(M,2*LQNA+1) = 0.0D0
          ENDIF
C
C         LOOP OVER ALL OCCUPIED SHELLS OF THIS KQN TYPE
          DO IOCC=1,NUMOCC(LQNA+1)
C
C           LL DENSITY CONTRIBUTIONS
            DLL = H2(IBAS     ,NBLC+IOCC)*H2(JBAS     ,NBLC+IOCC)
            DENLL(M,2*LQNA+1) = DENLL(M,2*LQNA+1) + QE(IOCC)*DLL
            DFNLL(M,2*LQNA+1) = DFNLL(M,2*LQNA+1) + QA(IOCC)*DLL
C
            IF(HMLTN.NE.'NORL') THEN
C
C             SL DENSITY CONTRIBUTIONS
              DSL = H2(IBAS+NBLC,NBLC+IOCC)*H2(JBAS     ,NBLC+IOCC)
              DENSL(M,2*LQNA+1) = DENSL(M,2*LQNA+1) + QE(IOCC)*DSL
              DFNSL(M,2*LQNA+1) = DFNSL(M,2*LQNA+1) + QA(IOCC)*DSL
C
C             SS DENSITY CONTRIBUTIONS
              DSS = H2(IBAS+NBLC,NBLC+IOCC)*H2(JBAS+NBLC,NBLC+IOCC)
              DENSS(M,2*LQNA+1) = DENSS(M,2*LQNA+1) + QE(IOCC)*DSS
              DFNSS(M,2*LQNA+1) = DFNSS(M,2*LQNA+1) + QA(IOCC)*DSS
C
            ENDIF
C
          ENDDO
        ENDDO
      ENDDO
C
C     ONE-BODY EIGENVALUE ENERGIES FOR OCCUPIED ELECTRONS
      DO IOCC=1,NUMOCC(LQNA+1)
        EH = EH + QA(IOCC)*RK2A2*W2(NBLC+IOCC)
      ENDDO
C
C     END LOOP OVER LQNA VALUES
100   CONTINUE
C
C     COULOMB AND BREIT ENERGIES HAVE BEEN DOUBLE-COUNTED
      EG = EG/2.0D0
      EB = EB/2.0D0
C
      E2 = EG+EB
C
C     TOTAL ATOMIC ENERGY IN THIS ITERATION
      ENEW = EH-EG-EB
C
C     RELATIVE CHANGE IN ENERGY
      ETEST = DABS((EPRV-ENEW)/ENEW)
C
C     WRITE THE ITERATION NUMBER AND THE TOTAL ENERGY
      WRITE(6,41) ITER,EH,E2,ENEW,ETEST
      WRITE(7,41) ITER,EH,E2,ENEW,ETEST
C
C     SUCCESSFUL CONVERGENCE
      IF(ETEST.LE.EEPS) THEN
        GOTO 1001
      ELSE
        EPRV = ENEW
      ENDIF
C
C     BARE NUCLEUS APPROXIMATION
      IF(HMLTN.EQ.'BARE') GOTO 1001
C
C     END LOOP OVER ITERATIONS
1000  CONTINUE
C
C     WARN USER THAT ATOMIC SCF DID NOT CONVERGE
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6, *) 'In HFSCF0: no convergence in ',ITER,' iterations.'
      WRITE(7, *) 'In HFSCF0: no convergence in ',ITER,' iterations.'
C
C     COVERGENCE SUCCESSFUL
1001  CONTINUE
C
C**********************************************************************C
C     CONSTRUCT MOLECULAR OVERLAP MATRIX                               C
C**********************************************************************C
C
C     LOOP OVER ALL OCCUPIED LQN VALUES
      DO 101 LQNA=0,LMXCONF
C
C     IMPORT BASIS FUNCTION EXPONENTS
      NFUNA = NFUNCT(LQNA+1,ICNT)
      DO IBAS=1,NFUNA
        EXLA(IBAS) = EXPSET(IBAS,LQNA+1,ICNT)
      ENDDO
C
C     ONLY NEED 'KQN =+LQN' IN SOME CASES
      IF(LQNA.EQ.0) GOTO 131
C
C     FIRST RELATIVISTIC QUANTUM NUMBER
      KAPA1 = LQNA
C
C     GENERATE OVERLAP MATRIX
      CALL OVRLP0(O1,EXLA,KAPA1,NFUNA)
C
131   CONTINUE
C
C     NEED 'KQN =-LQN-1' ALWAYS
C
C     SECOND RELATIVISTIC QUANTUM NUMBER
      KAPA2 =-LQNA-1
C
C     GENERATE OVERLAP MATRIX
      CALL OVRLP0(O2,EXLA,KAPA2,NFUNA)
C
C     NOW FILL IN THESE OVERLAP MATRICES...
C
C     ONLY NEED 'KQN =+LQN' IN SOME CASES
      IF(LQNA.EQ.0) GOTO 141
C
C     BEGIN LOOP OVER MQNA VALUES
      DO IMVAL=1,IABS(KAPA1)
C
C       COEFFICIENT MATRIX ADDRESSES FOR 'KQN =+LQN' CHOICE
        IL1 = LARGE(ICNT,2*LQNA  ,IMVAL*2-1)
        IL2 = LARGE(ICNT,2*LQNA  ,IMVAL*2  )
        IS1 = IL1 + NSHIFT
        IS2 = IL2 + NSHIFT
C
C       TRANSFER OVERLAP MATRIX TO COMMON ARRAY
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNA
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NFUNA
            LBAS = JBAS+NFUNA
C
C           MATRIX ELEMENTS
            OVAP(IL1+IBAS,IL1+JBAS) = DCMPLX(O1(IBAS,JBAS),0.0D0)
            OVAP(IL2+IBAS,IL2+JBAS) = DCMPLX(O1(IBAS,JBAS),0.0D0)
            
            IF(HMLTN.NE.'NORL') THEN
              OVAP(IS1+IBAS,IS1+JBAS) = DCMPLX(O1(KBAS,LBAS),0.0D0)
              OVAP(IS2+IBAS,IS2+JBAS) = DCMPLX(O1(KBAS,LBAS),0.0D0)
            ENDIF
C
          ENDDO
        ENDDO
      ENDDO
C
141   CONTINUE
C
C     NEED 'KQN =-LQN-1' ALWAYS
C
C     BEGIN LOOP OVER MQNA VALUES
      DO IMVAL=1,IABS(KAPA2)
C
C       COEFFICIENT MATRIX ADDRESSES FOR 'KQN =+LQN' CHOICE
        IL1 = LARGE(ICNT,2*LQNA+1,IMVAL*2-1)
        IL2 = LARGE(ICNT,2*LQNA+1,IMVAL*2  )
        IS1 = IL1 + NSHIFT
        IS2 = IL2 + NSHIFT
C
C       TRANSFER OVERLAP MATRIX TO COMMON ARRAY
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNA
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NFUNA
            LBAS = JBAS+NFUNA
C
C           MATRIX ELEMENTS
            OVAP(IL1+IBAS,IL1+JBAS) = DCMPLX(O2(IBAS,JBAS),0.0D0)
            OVAP(IL2+IBAS,IL2+JBAS) = DCMPLX(O2(IBAS,JBAS),0.0D0)
C
            IF(HMLTN.NE.'NORL') THEN
              OVAP(IS1+IBAS,IS1+JBAS) = DCMPLX(O2(KBAS,LBAS),0.0D0)
              OVAP(IS2+IBAS,IS2+JBAS) = DCMPLX(O2(KBAS,LBAS),0.0D0)
            ENDIF
C
          ENDDO
        ENDDO
      ENDDO
C
C     END LOOP OVER LQN VALUES
101   CONTINUE
C
C**********************************************************************C
C     WRITTEN SUMMARY                                                  C
C**********************************************************************C
C
C     SUMMARY OF ENERGY CONTRIBUTIONS
50    FORMAT(1X,A,24X,F19.12)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',20),'Atomic energies (Hartree units)'
      WRITE(7, *) REPEAT(' ',20),'Atomic energies (Hartree units)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,50) 'One-electron energy          ',EH
      WRITE(7,50) 'One-electron energy          ',EH
      WRITE(6,50) 'Two-electron energy (Coulomb)',EG
      WRITE(7,50) 'Two-electron energy (Coulomb)',EG
      IF(HMLTN.NE.'DHFB'.AND.HMLTN.NE.'DHFQ') GOTO 500
      WRITE(6,50) 'Two-electron energy (Breit)  ',EB
      WRITE(7,50) 'Two-electron energy (Breit)  ',EB
500   CONTINUE
      WRITE(6,50) 'Total energy                 ',ENEW
      WRITE(7,50) 'Total energy                 ',ENEW
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
C
C     UPDATE COUNTER FOR HIGHEST OCCUPIED ATOMIC ORBITAL
      IOCCM0 = IOCCML
C
C     ADD RESULTS FROM THIS ATOM TO MOLECULAR ENERGY
      ETOT = ETOT + ENEW
      EONE = EONE + EH
      ECLG = ECLG + EG
      EBRG = EBRG + EB
C
C     DFNOTE: BREIT0 FIX
      IF(HMLTN.EQ.'DHFB') THEN
        TITLE = 'Atomic_BREIT0'
        CALL ARRYPLT(BMAT,TITLE,NDIM)
        STOP
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE OVRLP0(OVAP,EXL,KQN,NFUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          OOOOOO  VV    VV RRRRRRR  LL      PPPPPPP   000000          C
C         OO    OO VV    VV RR    RR LL      PP    PP 00    00         C
C         OO    OO VV    VV RR    RR LL      PP    PP 00    00         C
C         OO    OO VV    VV RR    RR LL      PP    PP 00    00         C
C         OO    OO  VV  VV  RRRRRRR  LL      PPPPPPP  00    00         C
C         OO    OO   VVVV   RR    RR LL      PP       00    00         C
C          OOOOOO     VV    RR    RR LLLLLLL PP        000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  OVRLP0 CALCULATES THE ATOMIC OVERLAP MATRIX FOR SYMMETRY TYPE KQN.  C
C**********************************************************************C
      PARAMETER(MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4 HMLTN
C
      DIMENSION RN(MBS*MBS,4),OVAP(2*MBS,2*MBS),EXL(MBS)
C
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     DETERMINE THE LQN
      IF(KQN.LT.0) THEN
        LQN =-KQN-1
      ELSE
        LQN = KQN
      ENDIF
      RL = DFLOAT(LQN)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NFUN,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NFUN
        EI = EXL(IBAS)
        DO JBAS=1,NFUN
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
          T32 = RL+1.5D0
          T52 = RL+2.5D0
          E32 = EIJ**T32
          E52 = EIJ**T52
          SLL = 0.5D0*RN(M,1)*GAMHLF(2*LQN+3)/E32
          SSS = 2.0D0*RN(M,3)*GAMHLF(2*LQN+5)*EPR/E52
C
C         OVERLAP MATRIX ELEMENTS
          OVAP(IBAS     ,JBAS     ) = SLL
          IF(HMLTN.EQ.'NORL') GOTO 50
          OVAP(IBAS+NFUN,JBAS     ) = 0.0D0
          OVAP(JBAS     ,IBAS+NFUN) = 0.0D0
          OVAP(IBAS+NFUN,JBAS+NFUN) = SSS
50        CONTINUE
C         
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ONEEL0(HMAT,EXL,ZCRG,KQN,NFUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          OOOOOO  NN    NN EEEEEEEE EEEEEEEE LL      000000           C
C         OO    OO NNN   NN EE       EE       LL     00    00          C
C         OO    OO NNNN  NN EE       EE       LL     00    00          C
C         OO    OO NN NN NN EEEEEE   EEEEEE   LL     00    00          C
C         OO    OO NN  NNNN EE       EE       LL     00    00          C
C         OO    OO NN   NNN EE       EE       LL     00    00          C
C          OOOOOO  NN    NN EEEEEEEE EEEEEEEE LLLLLLL 000000           C
C                                                                      C
C -------------------------------------------------------------------- C
C  ONEEL0 CALCULATES THE ATOMIC DIRAC AND OVERLAP MATRICES FOR         C
C  SYMMETRY TYPE KQN, USING EVEN-TEMPERED SGTFS.                       C
C**********************************************************************C
      PARAMETER(MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4 HMLTN
C
      DIMENSION RN(MBS*MBS,4),HMAT(2*MBS,2*MBS),EXL(MBS)
C
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     DETERMINE THE LQN
      IF(KQN.LT.0) THEN
        LQN =-KQN-1
      ELSE
        LQN = KQN
      ENDIF
      RL = DFLOAT(LQN)
      G  = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NFUN,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NFUN
        EI = EXL(IBAS)
        DO JBAS=1,NFUN
          M    = M+1
          EJ   = EXL(JBAS)
          EIJ  = EI+EJ
          EPR  = EI*EJ
          T52  = RL+2.5D0
          E52  = EIJ**T52
C
C         LL OVERLAP
          ULL =-ZCRG*RN(M,1)*BNUCINT0(2*LQN+1,EIJ)
          IF(HMLTN.EQ.'NORL') THEN
            PLL = RN(M,1)*GAMHLF(2*LQN+5)*EPR/E52
          ELSE
            PLL = 0.0D0
          ENDIF
C
C         TRANSFER INTO ARRAY
          HMAT(IBAS,JBAS) = ULL + PLL
C
C         LS,SL AND SS OVERLAPS
          IF(HMLTN.EQ.'NORL') GOTO 50
C
C         SMALL COMPONENT BLOCK ADDRESSES
          KBAS = IBAS+NFUN
          LBAS = JBAS+NFUN
C
C         OVERLAPS, KINETIC ELEMENTS AND SS POTENTIAL INTEGRALS
          SSS = 2.0D0*RN(M,3)*GAMHLF(2*LQN+5)*EPR/E52
          PSL = 2.0D0*RN(M,2)*GAMHLF(2*LQN+5)*EPR/E52
C
          VSA  = 4.0D0*EPR*BNUCINT0(2*LQN+3,EIJ)
          IF(KQN.GT.0) THEN
            VSN =-2.0D0*EIJ*G*BNUCINT0(2*LQN+1,EIJ)
     &                   +G*G*BNUCINT0(2*LQN-1,EIJ)
          ELSE
            VSN = 0.0D0
          ENDIF
          USS =-ZCRG*RN(M,3)*(VSA+VSN)
C
C         TRANSFER INTO ARRAY
          HMAT(KBAS,JBAS) = CV*PSL
          HMAT(JBAS,KBAS) = HMAT(KBAS,JBAS)
          HMAT(KBAS,LBAS) = USS-2.0D0*CV*CV*SSS
C
50        CONTINUE
C         
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION BNUCINT0(N,ZETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  BBBBBBB  NN    NN UU    UU  CCCCCC  IIII NN    NN TTTTTTTT 000000   C
C  BB    BB NNN   NN UU    UU CC    CC  II  NNN   NN    TT   00    00  C
C  BB    BB NNNN  NN UU    UU CC        II  NNNN  NN    TT   00    00  C
C  BBBBBBB  NN NN NN UU    UU CC        II  NN NN NN    TT   00    00  C
C  BB    BB NN  NNNN UU    UU CC        II  NN  NNNN    TT   00    00  C
C  BB    BB NN   NNN UU    UU CC    CC  II  NN   NNN    TT   00    00  C
C  BBBBBBB  NN    NN  UUUUUU   CCCCCC  IIII NN    NN    TT    000000   C
C                                                                      C
C -------------------------------------------------------------------- C
C  BNUCINT0 CALCULATES A ONE-CENTER BARE NUCLEAR ATTRACTION INTEGRAL   C
C  GIVEN EXPONENT SUM ZETA, ORDER N (RADIAL POWER) AND LOCAL NUCLEAR   C
C  WIDTH PNUC. MAXIMUM ORDER IS N=21, SO TREATS LQNMAX=8 (k-TYPE).     C
C -------------------------------------------------------------------- C
C  BNUCINT0(N,ZETA) = INT{R^N*EXP(-ZETA*R^2)*ERF(SQRT(PNUC)*R)}.       C
C**********************************************************************C
      PARAMETER(MBS=26,MCT=6,MKP=9)
C
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     ROUTINE ONLY ALLOWS ODD PARAMETERS N
      IF(MOD(N,2).NE.1) THEN
        WRITE(6, *) 'In BNUCINT0: order N must be odd. N = ',N
        WRITE(7, *) 'In BNUCINT0: order N must be odd. N = ',N
      ENDIF
C
C     FACTORS NEEDED FOR ALL PARAMETERS N
      X   = ZETA/PNUC
      X5  = X*X*X*X*X
      T0  = PNUC+ZETA
      RAT = PNUC/T0
      TRM = 0.5D0*DSQRT(PNUC)/ZETA/DSQRT(T0)
      DO I=1,(N-1)/2
        TRM = 0.5D0*TRM*RAT/ZETA
      ENDDO
C
      IF(N.EQ.1) THEN
        TRM = TRM
      ELSEIF(N.EQ.3) THEN
        VA  = 2.0D0 + 3.0D0*X
        TRM = TRM*VA
      ELSEIF(N.EQ.5) THEN
        VA  = 8.0D0 + 20.0D0*X + 15.0D0*X*X
        TRM = TRM*VA
      ELSEIF(N.EQ.7) THEN
        VA  = 16.0D0 + 56.0D0*X + 70.0D0*X*X + 35.0D0*X*X*X
        TRM = 3.0D0*TRM*VA
      ELSEIF(N.EQ.9) THEN
        VA  = 128.0D0 + 576.0D0*X + 1008.0D0*X*X + 840.0D0*X*X*X 
        VB  = 315.0D0*X*X*X*X
        TRM = 3.0D0*TRM*(VA+VB)
      ELSEIF(N.EQ.11) THEN
        VA  = 256.0D0 + 1408.0D0*X + 3168.0D0*X*X + 3696.0D0*X*X*X 
        VB  = 2310.0D0*X*X*X*X + 693.0D0*X*X*X*X*X
        TRM = 15.0D0*TRM*(VA+VB)
      ELSEIF(N.EQ.13) THEN
        VA  = 1024.0D0 + 6656.0D0*X + 18304.0D0*X*X 
        VB  = 27456.0D0*X*X*X + 24024.0D0*X*X*X*X 
        VC  = 12012.0D0*X5 + 3003.0D0*X5*X
        TRM = 45.0D0*TRM*(VA+VB+VC)
      ELSEIF(N.EQ.15) THEN
        VA  = 2048.0D0 + 15360.0D0*X + 49920.0D0*X*X 
        VB  = 91520.0D0*X*X*X+ 102960.0D0*X*X*X*X + 72072.0D0*X5 
        VC  = 30030.0D0*X5*X + 6435.0D0*X5*X*X
        TRM = 315.0D0*TRM*(VA+VB+VC)
      ELSEIF(N.EQ.17) THEN
        VA  = 32768.0D0 + 278528.0D0*X + 1044480.0D0*X*X 
        VB  = 2263040.0D0*X*X*X + 3111680.0D0*X*X*X*X 
        VC  = 2800512.0D0*X5 + 1633632.0D0*X5*X + 583440.0D0*X5*X*X
        VD  = 109395.0D0*X5*X*X*X
        TRM = 315.0D0*TRM*(VA+VB+VC+VD)
      ELSEIF(N.EQ.19) THEN
        VA  = 65536.0D0 + 6222592.0D0*X + 2646016.0D0*X*X 
        VB  = 6615040.0D0*X*X*X + 10749440.0D0*X*X*X*X 
        VC  = 11824384.0D0*X5 + 8868288.0D0*X5*X + 4434144.0D0*X5*X*X 
        VD  = 1385670.0D0*X5*X*X*X + 230945.0D0*X5*X*X*X*X
        TRM = 2835.0D0*TRM*(VA+VB+VC+VD)
      ELSEIF(N.EQ.21) THEN
        VA  = 262144.0D0 + 2752512.0D0*X + 13074432.0D0*X*X 
        VB  = 37044224.0D0*X*X*X + 69457920.0D0*X*X*X*X 
        VC  = 90295296.0D0*X5 + 82770688.0D0*X5*X 
        VD  = 53209728.0D0*X5*X*X + 23279256.0D0*X5*X*X*X 
        VE  = 6466460.0D0*X5*X*X*X*X + 969969.0D0*X5*X5
        TRM = 14175.0D0*TRM*(VA+VB+VC+VD+VE)
      ELSE 
        WRITE(6, *) 'In BNUCINT0: order N too large. N = ',N
        WRITE(7, *) 'In BNUCINT0: order N too large. N = ',N
      ENDIF
C
C     TRANSFER DATA TO BNUCINT0
      BNUCINT0 = TRM
C
      RETURN
      END
C
C
      SUBROUTINE COULOMB0(G11,G21,G12,G22,D1,D2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    CCCCCC   OOOOOO  UU    UU LL      MM       MM BBBBBBB   000000    C
C   CC    CC OO    OO UU    UU LL      MMM     MMM BB    BB 00    00   C
C   CC       OO    OO UU    UU LL      MMMM   MMMM BB    BB 00    00   C
C   CC       OO    OO UU    UU LL      MM MM MM MM BBBBBBB  00    00   C
C   CC       OO    OO UU    UU LL      MM  MMM  MM BB    BB 00    00   C
C   CC    CC OO    OO UU    UU LL      MM   M   MM BB    BB 00    00   C
C    CCCCCC   OOOOOO   UUUUUU  LLLLLLL MM       MM BBBBBBB   000000    C
C                                                                      C
C -------------------------------------------------------------------- C
C  COULOMB0 CONSTRUCTS THE ATOMIC COULOMB MATRIX FROM RADIAL DIRECT    C
C  AND EXCHANGE INTEGRALS AND A MEAN-FIELD CHARGE DENSITY.             C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MNU=MKP+1,MAB=2*MNU+5)
C
      CHARACTER*4 HMLTN
C
      DIMENSION G11(2*MBS,2*MBS),G21(2*MBS,2*MBS),
     &          G12(2*MBS,2*MBS),G22(2*MBS,2*MBS)
      DIMENSION D1(MB2,3),D2(MB2,3),RN(MB2,4)
C
      COMMON/BSIJ/EIJ(MAB),RNIJ(3),EI,EJ
      COMMON/BSKL/EK(MB2),EKL0(MB2),EIK(MB2,MAB),EJL(MB2,MAB),
     &            EL(MB2),EKL1(MB2),EKL(MB2,MAB),RNKL(MB2,4),
     &            B1(MB2,MAB,MAB),B2(MB2,MAB,MAB),B3(MB2,MAB,MAB),
     &            B4(MB2,MAB,MAB),B5(MB2,MAB,MAB),B6(MB2,MAB,MAB)
      COMMON/BSQN/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/BSXP/EXLA(MBS),EXLB(MBS)
      COMMON/CLRE/RKLLLL(MB2,4),RKSSSS(MB2,4),RKSLSL(MB2,4),
     &            RJLLLL(MB2,4),RJSSSS(MB2,4),RJLLSS(MB2,4),
     &            RJSSLL(MB2,4)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
C
C     INITIALISE COULOMB MATRIX
      DO IBAS=1,2*MBS
        DO JBAS=1,2*MBS
          G11(IBAS,JBAS) = 0.0D0
          G21(IBAS,JBAS) = 0.0D0
          G12(IBAS,JBAS) = 0.0D0
          G22(IBAS,JBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     GENERATE 'KL' EXPONENT POWERS FOR LATER B-INTEGRAL GENERATION
      CALL KLSET
C
C     GENERATE A BATCH OF NORMALISATION CONSTANTS
      CALL RNORM0(RN,EXLA,NFUNA,LQNA)
C
C     ITERATE OVER ALL MATRIX ELEMENTS
      IJ = 0
      DO IBAS=1,NFUNA
        EI = EXLA(IBAS)
        DO JBAS=1,NFUNA
C
C         BASIS EXPONENT COMBINATIONS FOR LATER B-INTEGRAL GENERATION
          IJ = IJ+1
C
          EJ = EXLA(JBAS)
          EIJ0 = EI+EJ
          EIJR = DSQRT(EIJ0)
          EIJA = EIJ0**(-LQNA)
          DO N=1,6
            EIJ(N) = EIJA
            EIJA   = EIJA/EIJR
          ENDDO
          RNIJ(1) = RN(IJ,1)
          RNIJ(2) = RN(IJ,2)
          RNIJ(3) = RN(IJ,3)
C
C         GENERATE A BATCH OF BETA INTEGRALS
          CALL BINT
C
C         GENERATE BATCH OF RADIAL INTEGRALS (J AND K MATRICES)
          CALL ERI0
C
          IF(HMLTN.EQ.'NORL') THEN
C         NON-RELATIVISTIC HAMILTONIAN
C
C           INITIALISE COUNTER
            GLL = 0.0D0
C
C           BUILD THE FOCK MATRIX
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,4)*D2(M,1) - RKLLLL(M,4)*D2(M,1)
            ENDDO
C
C           TRANSFER COUNTER VALUE TO COULOMB MATRIX
            G22(IBAS,JBAS) = GLL
C
          ELSE
C         RELATIVISTIC HAMILTONIAN
C
C           SMALL-COMPONENT MATRIX ADDRESSES
            KBAS = IBAS + NFUNA
            LBAS = JBAS + NFUNA
C
C    (11)   KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 10
C
C           INITIALISE COUNTERS          
            GLL = 0.0D0
            GSL = 0.0D0
            GSS = 0.0D0
C
C           SUM OVER MEAN FIELD CONTRIBUTIONS FOR THIS BASIS PAIR
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,1)*D1(M,1) - RKLLLL(M,1)*D1(M,1)
     &                  + RJLLSS(M,1)*D1(M,3)
              GSL = GSL                       - RKSLSL(M,1)*D1(M,2)
              GSS = GSS + RJSSSS(M,1)*D1(M,3) - RKSSSS(M,1)*D1(M,3)
     &                  + RJSSLL(M,1)*D1(M,1)
            ENDDO
C
C           TRANSFER COUNTER VALUES TO COULOMB MATRIX
            G11(IBAS,JBAS) = GLL
            G11(KBAS,JBAS) = GSL
            G11(JBAS,KBAS) = GSL
            G11(KBAS,LBAS) = GSS
C
10          CONTINUE
C
C    (21)   KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNB.EQ.0) GOTO 11
C
C           INITIALISE COUNTERS          
            GLL = 0.0D0
            GSL = 0.0D0
            GSS = 0.0D0
C
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,2)*D1(M,1) - RKLLLL(M,2)*D1(M,1)
     &                  + RJLLSS(M,2)*D1(M,3)
              GSL = GSL                       - RKSLSL(M,2)*D1(M,2)
              GSS = GSS + RJSSSS(M,2)*D1(M,3) - RKSSSS(M,2)*D1(M,3)
     &                  + RJSSLL(M,2)*D1(M,1)
            ENDDO
C
C           TRANSFER COUNTER VALUES TO COULOMB MATRIX
            G21(IBAS,JBAS) = GLL
            G21(KBAS,JBAS) = GSL
            G21(JBAS,KBAS) = GSL
            G21(KBAS,LBAS) = GSS
C
11          CONTINUE

C
C    (12)   KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNA.EQ.0) GOTO 12
C
C           INITIALISE COUNTERS          
            GLL = 0.0D0
            GSL = 0.0D0
            GSS = 0.0D0
C
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,3)*D2(M,1) - RKLLLL(M,3)*D2(M,1)
     &                  + RJLLSS(M,3)*D2(M,3)
              GSL = GSL                       - RKSLSL(M,3)*D2(M,2)
              GSS = GSS + RJSSSS(M,3)*D2(M,3) - RKSSSS(M,3)*D2(M,3)
     &                  + RJSSLL(M,3)*D2(M,1)
            ENDDO
C
C           TRANSFER COUNTER VALUES TO COULOMB MATRIX
            G12(IBAS,JBAS) = GLL
            G12(KBAS,JBAS) = GSL
            G12(JBAS,KBAS) = GSL
            G12(KBAS,LBAS) = GSS
C
12          CONTINUE
C
C    (22)   KQNA < 0 AND KQNB < 0  CONTRIBUTIONS (CANNOT SKIP)
C
C           INITIALISE COUNTERS          
            GLL = 0.0D0
            GSL = 0.0D0
            GSS = 0.0D0
C
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,4)*D2(M,1) - RKLLLL(M,4)*D2(M,1)
     &                  + RJLLSS(M,4)*D2(M,3)
              GSL = GSL                       - RKSLSL(M,4)*D2(M,2)
              GSS = GSS + RJSSSS(M,4)*D2(M,3) - RKSSSS(M,4)*D2(M,3)
     &                  + RJSSLL(M,4)*D2(M,1)
            ENDDO
C
C           TRANSFER COUNTER VALUES TO COULOMB MATRIX
            G22(IBAS,JBAS) = GLL
            G22(KBAS,JBAS) = GSL
            G22(JBAS,KBAS) = GSL
            G22(KBAS,LBAS) = GSS
C
          ENDIF
C
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE BREIT0(B11,B21,B12,B22,D1,D2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           BBBBBBB  RRRRRRR  EEEEEEEE IIII TTTTTTTT 000000            C
C           BB    BB RR    RR EE        II     TT   00    00           C
C           BB    BB RR    RR EE        II     TT   00    00           C
C           BBBBBBB  RR    RR EEEEEE    II     TT   00    00           C
C           BB    BB RRRRRRR  EE        II     TT   00    00           C
C           BB    BB RR    RR EE        II     TT   00    00           C
C           BBBBBBB  RR    RR EEEEEEEE IIII    TT    000000            C
C                                                                      C
C -------------------------------------------------------------------- C
C  BREIT0 CONSTRUCTS THE ATOMIC BREIT MATRIX FROM RADIAL DIRECT AND    C
C  EXCHANGE INTEGRALS AND A MEAN-FIELD CHARGE DENSITY.                 C
C -------------------------------------------------------------------- C
C  DFNOTE: NOT FINISHED YET. EXCHANGE TERMS INCORRECT FOR S-TYPE GTFS. C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MNU=MKP+1,MAB=2*MNU+5)
C
      DIMENSION B11(2*MBS,2*MBS),B21(2*MBS,2*MBS),
     &          B12(2*MBS,2*MBS),B22(2*MBS,2*MBS)
      DIMENSION D1(MB2,3),D2(MB2,3),RN(MB2,4)
C
      COMMON/BSIJ/EIJ(MAB),RNIJ(3),EI,EJ
      COMMON/BSKL/EK(MB2),EKL0(MB2),EIK(MB2,MAB),EJL(MB2,MAB),
     &            EL(MB2),EKL1(MB2),EKL(MB2,MAB),RNKL(MB2,4),
     &            B1(MB2,MAB,MAB),B2(MB2,MAB,MAB),B3(MB2,MAB,MAB),
     &            B4(MB2,MAB,MAB),B5(MB2,MAB,MAB),B6(MB2,MAB,MAB)
      COMMON/BSQN/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/BSXP/EXLA(MBS),EXLB(MBS)
      COMMON/BTRE/RKLLSS(MB2,4),RKSLSL(MB2,4),RKSSLL(MB2,4),
     &            RMSLSL(MB2,4)
C
C     GENERATE 'KL' EXPONENT POWERS FOR LATER B-INTEGRAL GENERATION
      CALL KLSET
C
C     GENERATE A BATCH OF NORMALISATION CONSTANTS
      CALL RNORM0(RN,EXLA,NFUNA,LQNA)
C
C     INITIALISE BREIT MATRIX
      DO IBAS=1,2*NFUNA
        DO JBAS=1,2*NFUNA
          B11(IBAS,JBAS) = 0.0D0
          B21(IBAS,JBAS) = 0.0D0
          B12(IBAS,JBAS) = 0.0D0
          B22(IBAS,JBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     ITERATE OVER ALL MATRIX ELEMENTS
      IJ = 0
      DO IBAS=1,NFUNA
        EI = EXLA(IBAS)
        DO JBAS=1,NFUNA
C         BASIS EXPONENT COMBINATIONS FOR LATER B-INTEGRAL GENERATION
          IJ = IJ+1
C
          KBAS = IBAS + NFUNA
          LBAS = JBAS + NFUNA
C
          EJ = EXLA(JBAS)
          EIJ0 = EI + EJ
          EIJR = DSQRT(EIJ0)
          EIJA = EIJ0**(-LQNA)
          DO N=1,6
            EIJ(N) = EIJA
            EIJA   = EIJA/EIJR
          ENDDO
          RNIJ(1) = RN(IJ,1)
          RNIJ(2) = RN(IJ,2)
          RNIJ(3) = RN(IJ,3)
C
C         GENERATE A BATCH OF BETA INTEGRALS
          CALL BINT
C
C         GENERATE BATCH OF RADIAL INTEGRALS (J AND K MATRICES)
          CALL BII0
C
C    (11) KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 10
C
C         INITIALISE COUNTERS
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0
C
C         SUM OVER MEAN FIELD CONTRIBUTIONS FOR THIS BASIS PAIR
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,1)*D1(M,3)
            BSL = BSL + RKSLSL(M,1)*D1(M,2) + RMSLSL(M,1)*D1(M,2)
            BSS = BSS + RKSSLL(M,1)*D1(M,1)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B11(IBAS,JBAS) = BLL
          B11(KBAS,JBAS) = BSL
          B11(JBAS,KBAS) = BSL
          B11(KBAS,LBAS) = BSS
C
10        CONTINUE
C
C    (21) KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNB.EQ.0) GOTO 11
C
C         INITIALISE COUNTERS
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0
C
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,2)*D1(M,3)
            BSL = BSL + RKSLSL(M,2)*D1(M,2) + RMSLSL(M,2)*D1(M,2)
            BSS = BSS + RKSSLL(M,2)*D1(M,1)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B21(IBAS,JBAS) = BLL
          B21(KBAS,JBAS) = BSL
          B21(JBAS,KBAS) = BSL
          B21(KBAS,LBAS) = BSS
C
11        CONTINUE

C
C    (12) KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0) GOTO 12
C
C         INITIALISE COUNTERS          
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0
C
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,3)*D2(M,3)
            BSL = BSL + RKSLSL(M,3)*D2(M,2) + RMSLSL(M,3)*D2(M,2)
            BSS = BSS + RKSSLL(M,3)*D2(M,1)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B12(IBAS,JBAS) = BLL
          B12(KBAS,JBAS) = BSL
          B12(JBAS,KBAS) = BSL
          B12(KBAS,LBAS) = BSS
C
12        CONTINUE
C
C    (22) KQNA < 0 AND KQNB < 0  CONTRIBUTIONS (CANNOT SKIP)
C
C         INITIALISE COUNTERS          
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0
C
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,4)*D2(M,3)
            BSL = BSL + RKSLSL(M,4)*D2(M,2) + RMSLSL(M,4)*D2(M,2)
            BSS = BSS + RKSSLL(M,4)*D2(M,1)
            IF(IBAS.EQ.1.AND.JBAS.EQ.1) THEN
              WRITE(*,*) M,RKSSLL(M,4),D2(M,1)
            ENDIF
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B22(IBAS,JBAS) = BLL
          B22(KBAS,JBAS) = BSL
          B22(JBAS,KBAS) = BSL
          B22(KBAS,LBAS) = BSS
C
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE KLSET
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             KK    KK LL       SSSSSS  EEEEEEEE TTTTTTTT              C
C             KK   KK  LL      SS    SS EE          TT                 C
C             KK  KK   LL      SS       EE          TT                 C
C             KKKKK    LL       SSSSSS  EEEEEE      TT                 C
C             KK  KK   LL            SS EE          TT                 C
C             KK   KK  LL      SS    SS EE          TT                 C
C             KK    KK LLLLLLLL SSSSSS  EEEEEEEE    TT                 C
C                                                                      C
C -------------------------------------------------------------------- C
C  KLSET GENERATES LISTS OF 'KL' EXPONENT POWERS FOR LATER USE         C
C  IN BETA INTEGRAL CONSTRUCTION AND ELECTRON REPULSION INTEGRALS.     C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MNU=MKP+1,MAB=2*MNU+5)
C
      COMMON/ANGL/BK(MNU,4),DK(MNU,4),HK(MNU,4),FK(MNU,4),GM(MNU,4),
     &            NUS(MNU),NUNUM,NUMIN,NUMAX
      COMMON/BSIJ/EIJ(MAB),RNIJ(3),EI,EJ
      COMMON/BSKL/EK(MB2),EKL0(MB2),EIK(MB2,MAB),EJL(MB2,MAB),
     &            EL(MB2),EKL1(MB2),EKL(MB2,MAB),RNKL(MB2,4),
     &            B1(MB2,MAB,MAB),B2(MB2,MAB,MAB),B3(MB2,MAB,MAB),
     &            B4(MB2,MAB,MAB),B5(MB2,MAB,MAB),B6(MB2,MAB,MAB)
      COMMON/BSQN/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/BSXP/EXLA(MBS),EXLB(MBS)
      COMMON/INDX/IKIND(MB2),JLIND(MB2)
C
C     GENERATE INDICES AND EXPONENT COMBINATIONS
      M = 0
      DO KBAS=1,NFUNB
        EK0 = EXLB(KBAS)
        DO LBAS=1, NFUNB
          M   = M+1
          EL0 = EXLB(LBAS)
          IKIND(M) = KBAS
          JLIND(M) = LBAS
          EK(M)    = EK0
          EL(M)    = EL0
          EKL0(M)  = EK0 + EL0
          EKL1(M)  = EK0*EL0
          EKLR     = DSQRT(EKL0(M))
          EKPW     = EKL0(M)**LQNB
          EKLA     = 1.0D0/EKPW
          DO N=1,6
            EKL(M,N) = EKLA
            EKLA     = EKLA/EKLR
          ENDDO
        ENDDO
      ENDDO
C
C     NORMALISATION CONSTANTS   
      CALL RNORM0(RNKL,EXLB,NFUNB,LQNB)
C
      RETURN
      END
C
C
      SUBROUTINE BINT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                   BBBBBBB  IIII NN    NN TTTTTTTT                    C
C                   BB    BB  II  NNN   NN    TT                       C
C                   BB    BB  II  NNNN  NN    TT                       C
C                   BBBBBBB   II  NN NN NN    TT                       C
C                   BB    BB  II  NN  NNNN    TT                       C
C                   BB    BB  II  NN   NNN    TT                       C
C                   BBBBBBB  IIII NN    NN    TT                       C
C                                                                      C
C -------------------------------------------------------------------- C
C  BINT GENERATES A BATCH OF BETA INTEGRALS FOR LATER USE IN THE       C
C  CONSTRUCTION OF TWO-ELECTRON INTERACTION INTEGRALS (G/B MATRIX).    C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MNU=MKP+1,MAB=2*MNU+5)
C
      CHARACTER*4 HMLTN
C
      DIMENSION XJ(MB2,2),XK(MB2,2),
     &          RTIK(MB2),RTJL(MB2),RTIK0(MBS),RTJL0(MBS),
     &          PTIK0(MBS),PTJL0(MBS),TTIK0(MBS),TTJL0(MBS)
      DIMENSION BETA(MB2),BTA1(MB2),XROOT(MB2),
     &          BIN(MB2),TRM(MB2),IAA(2),IBB(2)
C
      COMMON/ANGL/BK(MNU,4),DK(MNU,4),HK(MNU,4),FK(MNU,4),GM(MNU,4),
     &            NUS(MNU),NUNUM,NUMIN,NUMAX
      COMMON/BSIJ/EIJ(MAB),RNIJ(3),EI,EJ
      COMMON/BSKL/EK(MB2),EKL0(MB2),EIK(MB2,MAB),EJL(MB2,MAB),
     &            EL(MB2),EKL1(MB2),EKL(MB2,MAB),RNKL(MB2,4),
     &            B1(MB2,MAB,MAB),B2(MB2,MAB,MAB),B3(MB2,MAB,MAB),
     &            B4(MB2,MAB,MAB),B5(MB2,MAB,MAB),B6(MB2,MAB,MAB)
      COMMON/BSQN/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/BSXP/EXLA(MBS),EXLB(MBS)
      COMMON/INDX/IKIND(MB2),JLIND(MB2)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
C
C     TENSOR ORDER AND RELEVANT POWER
      NUMAX  = NUS(NUNUM)
      IPOWER = LQNA+LQNB-NUMAX
C
C     TEMPORARY ARRAYS USED ONLY IN NEXT LOOP
      DO KBAS=1,NFUNB
        TTIK0(KBAS) = EI+EXLB(KBAS)
        TTJL0(KBAS) = EJ+EXLB(KBAS)
        RTIK0(KBAS) = DSQRT(TTIK0(KBAS))
        RTJL0(KBAS) = DSQRT(TTJL0(KBAS))
        PTIK0(KBAS) = RTIK0(KBAS)**(-IPOWER)
        PTJL0(KBAS) = RTJL0(KBAS)**(-IPOWER)
      ENDDO
C
C     GENERATE ARRAYS XJ(MB2,2) AND XK(MB2,2) FOR LOCAL USE
C     AND SEED THE COMMON ARRAYS EIK(MAXM,MAB) AND EJL(MAX,MAB)
      TIJ0 = EI + EJ
      DO M=1,MAXM
        TIK0     = TTIK0(IKIND(M))
        TJL0     = TTJL0(JLIND(M))
        TIJKL    = TIJ0+EKL0(M)
        XJ(M,1)  = TIJ0/TIJKL
        XJ(M,2)  = EKL0(M)/TIJKL
        XK(M,1)  = TIK0/TIJKL
        XK(M,2)  = TJL0/TIJKL
        RTIK(M)  = RTIK0(IKIND(M))
        RTJL(M)  = RTJL0(JLIND(M))
        EIK(M,1) = PTIK0(IKIND(M))
        EJL(M,1) = PTJL0(JLIND(M))
      ENDDO
C
C     GENERATE COMMON ARRAYS EIK(MAXM,MAB) AND EJL(MAX,MAB)
      DO IV=2,2*NUMAX+6
        DO M=1,MAXM
          EIK(M,IV) = EIK(M,IV-1)/RTIK(M)
          EJL(M,IV) = EJL(M,IV-1)/RTJL(M)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     GENERATE ALL OF THE INCOMPLETE BETA FUNCTIONS FOR J-TYPE         C
C**********************************************************************C
C
C     PARAMETER NVALS USED TO DEFINE J- OR K- TYPE
      NVALS = 3
C
C     LOOP OVER ORDER FOR FIRST INDEX
      DO I1=1,NVALS
        NT1 = 2*(I1-1)
        IAA(1) = 2*LQNA+NT1+1
        IAA(2) = 2*LQNB+NT1+1
C
C       LOOP OVER ORDER FOR SECOND INDEX
        DO I2=1,NVALS
          NT2 = 2*(I2-1)
          IBB(1) = 2*LQNB+NT2
          IBB(2) = 2*LQNA+NT2
C
C         LOOP OVER BETA INTEGRAL TYPE
          DO IBETA=1,2
            IA =(IAA(IBETA)-1)/2
            IB = IBB(IBETA)   /2
C
C ***       BEGIN CONDITIONAL STATEMENT OVER IB VALUES
C >>>       CASE 1: IB > 1
            IF(IB.GT.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BTA1(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
              RA  = X
              RB  = DFLOAT(1-IB)
              RC  = X+1.0D0
              RD  = 1.0D0
              RCD = RC*RD
              FCT = RA*RB/RCD
              DO M=1,MAXM
                TRM(M) = FCT*XJ(M,IBETA)
                BIN(M) = 1.0D0+TRM(M)
              ENDDO
              RA = RA + 1.0D0
              RB = RB + 1.0D0
              RC = RC + 1.0D0
              RD = RD + 1.0D0
              DO IT=2,IB-1
                RCD = RC*RD
                FCT = RA*RB/RCD
                DO M=1,MAXM
                  TRM(M) = FCT*TRM(M)*XJ(M,IBETA)
                  BIN(M) = BIN(M)+TRM(M)
                ENDDO
                RA = RA + 1.0D0
                RB = RB + 1.0D0
                RC = RC + 1.0D0
                RD = RD + 1.0D0
              ENDDO
              DO M=1,MAXM
                BETA(M)     = BTA1(M)*BIN(M)
                B2(M,I1,I2) = BETA(M)
              ENDDO
C
C >>>       CASE 2: IB = 1
            ELSEIF(IB.EQ.1) THEN
              X  = DFLOAT(IA) + 0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BETA(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
C
C >>>       CASE 3: IB = 0
            ELSEIF(IB.EQ.0) THEN
              DO M=1,MAXM
                XROOT(M) = DSQRT(XJ(M,IBETA))
                DEN      = (1.0D0-XROOT(M))
                RAT      = (1.0D0+XROOT(M))/DEN
                BTA1(M)  = DLOG(RAT)
                BIN(M)   = 1.0D0
                TRM(M)   = XJ(M,IBETA)
              ENDDO
              IF(IA.GT.1) THEN
                DO K=2,IA
                  KK = 2*K-1
                  RK = DFLOAT(KK)
                  X  = 1.0D0/RK
                  DO M=1,MAXM
                    BIN(M) = BIN(M)+X*TRM(M)
                    TRM(M) = TRM(M)*XJ(M,IBETA)
                  ENDDO
                ENDDO
                DO M=1,MAXM
                  BETA(M) = BTA1(M) - 2.0D0*XROOT(M)*BIN(M)
                ENDDO
              ELSEIF(IA.EQ.1) THEN
                DO M=1,MAXM
                  BETA(M) = BTA1(M) - 2.0D0*XROOT(M)
                ENDDO
              ELSE
                DO M=1,MAXM
                  BETA(M) = BTA1(M)
                ENDDO
              ENDIF
C ***       END CONDITIONAL STATEMENT OVER IB VALUES
            ENDIF
C
            IF(IBETA.EQ.1) THEN
              DO M=1,MAXM
                B1(M,I1,I2) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXM
                B2(M,I1,I2) = BETA(M)
              ENDDO
            ENDIF
C
C         END LOOPS OVER IBETA AND INDICES I1, I2  
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     GENERATE ALL OF THE INCOMPLETE BETA FUNCTIONS FOR K-TYPE         C
C**********************************************************************C
C
C     PARAMETER NVALS USED TO DEFINE J- OR K- TYPE
      NVALS = ((NUS(NUNUM)-NUS(1))/2)+3
C
C     LOOP OVER ORDER FOR FIRST INDEX
      DO I1=1,NVALS
        NA = NUS(1)+2*(I1-1)
        IAA(1) = LQNA + LQNB + NA + 1
        IAA(2) = LQNA + LQNB + NA + 1
C
C       LOOP OVER ORDER FOR SECOND INDEX
        DO I2=1,NVALS
          NB = 2*(I2-1)-NUS(NUNUM)
          IBB(1) = LQNA + LQNB + NB
          IBB(2) = LQNA + LQNB + NB
C
C         LOOP OVER BETA INTEGRAL TYPE
          DO IBETA=1,2
            IA = (IAA(IBETA)-1)/2
            IB =  IBB(IBETA)   /2
C
C ***       BEGIN CONDITIONAL STATEMENT OVER IB VALUES
C >>>       CASE 1: IB > 1
            IF(IB.GT.1) THEN
              X  = DFLOAT(IA) + 0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BTA1(M) = (DSQRT(XK(M,IBETA))**IX)/X
              ENDDO
              RA  = X
              RB  = DFLOAT(1-IB)
              RC  = X + 1.0D0
              RD  = 1.0D0
              RCD  = RC*RD
              FCT = RA*RB/RCD
              DO M=1,MAXM
                TRM(M) = FCT*XK(M,IBETA)
                BIN(M) = 1.0D0+TRM(M)
              ENDDO
              RA  = RA + 1.0D0
              RB  = RB + 1.0D0
              RC  = RC + 1.0D0
              RD  = RD + 1.0D0
              DO IT=2,IB-1
                RCD = RC*RD
                FCT = RA*RB/RCD
                DO M=1,MAXM
                  TRM(M) = FCT*TRM(M)*XK(M,IBETA)
                  BIN(M) = BIN(M)+TRM(M)
                ENDDO
                RA = RA + 1.0D0
                RB = RB + 1.0D0
                RC = RC + 1.0D0
                RD = RD + 1.0D0
              ENDDO
              DO M=1,MAXM
                BETA(M) = BTA1(M)*BIN(M)
              ENDDO
C
C >>>       CASE 2: IB = 1
            ELSEIF(IB.EQ.1) THEN
              X  = DFLOAT(IA) + 0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BETA(M)  = (DSQRT(XK(M,IBETA))**IX)/X
              ENDDO
C
C >>>       CASE 3: IB = 0
            ELSEIF(IB.EQ.0) THEN
              DO M=1,MAXM
                XROOT(M) = DSQRT(XK(M,IBETA))
                DEN      = (1.0D0-XROOT(M))
                RAT      = (1.0D0+XROOT(M))/DEN
                BTA1(M)  = DLOG(RAT)
                BIN(M)   = 1.0D0
                TRM(M)   = XK(M,IBETA)
              ENDDO
C
              IF(IA.GT.1) THEN
                DO K=2,IA
                  KK = 2*K-1
                  RK = DFLOAT(KK)
                  X  = 1.0D0/RK
                  DO M=1,MAXM
                    BIN(M) = BIN(M)+X*TRM(M)
                    TRM(M) = TRM(M)*XK(M,IBETA)
                  ENDDO
                ENDDO
                DO M=1,MAXM
                  BETA(M) = BTA1(M) - 2.0D0*XROOT(M)*BIN(M)
                ENDDO
              ELSEIF(IA.EQ.1) THEN
                DO M=1,MAXM
                  BETA(M) = BTA1(M) - 2.0D0*XROOT(M)
                ENDDO
              ELSE
                DO M=1,MAXM
                  BETA(M) = BTA1(M)
                ENDDO
              ENDIF
C ***       END CONDITIONAL STATEMENT OVER IB VALUES
            ENDIF
C
            IF(IBETA.EQ.1) THEN
              DO M=1,MAXM
                B3(M,I1,I2) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXM
                B4(M,I1,I2) = BETA(M)
              ENDDO
            ENDIF
C
C         END LOOPS OVER IBETA AND INDICES I1, I2
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     GENERATE INCOMPLETE BETA FUNCTIONS FOR BREIT MATRIX ELEMENTS     C
C -------------------------------------------------------------------- C
C     DFNOTE: THIS SHOULD LOOK THE SAME AS B3 AND B4, EXCEPT THE       C
C             BETA INTEGRAL INDICES (IA,IB;Z) AND (IA',IB',Z') WILL    C
C             REQUIRE DIFFERENT POWERS.                                C
C**********************************************************************C
C
      IF(HMLTN.NE.'DHFB'.AND.HMLTN.NE.'DHFQ') GOTO 70
C
C     TOTAL NUMBER OF P,Q VALUES NEEDED (ACCOUNT FOR NU AND T=S CHOICE)
      NVALS = ((NUS(NUNUM)-NUS(1))/2) + 2
C
C     LOOP OVER ORDER FOR FIRST INDEX
      DO I1=1,NVALS
C
C       I1 PARAMETER VALUE FOR B5 AND B6
        IAA(1) = LQNA + LQNB + NUS(1) + 2*I1
        IAA(2) = LQNA + LQNB + NUS(1) + 2*I1
C
C       LOOP OVER ORDER FOR SECOND INDEX
        DO I2=1,NVALS
C
C         I2 PARAMETER VALUE FOR B5 AND B6
          IBB(1) = LQNA + LQNB - NUS(NUNUM) + 2*I2 - 1
          IBB(2) = LQNA + LQNB - NUS(NUNUM) + 2*I2 - 1
C
C         LOOP OVER BETA INTEGRAL TYPE (B5 OR B6)
          DO IBETA=1,2
            IA = (IAA(IBETA)-1)/2
            IB =  IBB(IBETA)   /2
C
C ***       BEGIN CONDITIONAL STATEMENT OVER IB VALUES
C >>>       CASE 1: IB > 1
            IF(IB.GT.1) THEN
              X  = DFLOAT(IA) + 0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BTA1(M) = (DSQRT(XK(M,IBETA))**IX)/X
              ENDDO
              RA  = X
              RB  = DFLOAT(1-IB)
              RC  = X + 1.0D0
              RD  = 1.0D0
              RCD = RC*RD
              FCT = RA*RB/RCD
              DO M=1,MAXM
                TRM(M) = FCT*XK(M,IBETA)
                BIN(M) = 1.0D0+TRM(M)
              ENDDO
              RA = RA + 1.0D0
              RB = RB + 1.0D0
              RC = RC + 1.0D0
              RD = RD + 1.0D0
              DO IT=2,IB-1
                RCD = RC*RD
                FCT = RA*RB/RCD
                DO M=1,MAXM
                  TRM(M) = FCT*TRM(M)*XK(M,IBETA)
                  BIN(M) = BIN(M)+TRM(M)
                ENDDO
                RA = RA + 1.0D0
                RB = RB + 1.0D0
                RC = RC + 1.0D0
                RD = RD + 1.0D0
              ENDDO
              DO M=1,MAXM
                BETA(M) = BTA1(M)*BIN(M)
              ENDDO
C
C >>>       CASE 2: IB = 1
            ELSEIF(IB.EQ.1) THEN
              X  = DFLOAT(IA) + 0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BETA(M)  = (DSQRT(XK(M,IBETA))**IX)/X
              ENDDO
C
C >>>       CASE 3: IB = 0
            ELSEIF(IB.EQ.0) THEN
              DO M=1,MAXM
                XROOT(M) = DSQRT(XK(M,IBETA))
                DEN      = (1.0D0-XROOT(M))
                RAT      = (1.0D0+XROOT(M))/DEN
                BTA1(M)  = DLOG(RAT)
                BIN(M)   = 1.0D0
                TRM(M)   = XK(M,IBETA)
              ENDDO
C
              IF(IA.GT.1) THEN
                DO K=2,IA
                  KK = 2*K-1
                  RK = DFLOAT(KK)
                  X  = 1.0D0/RK
                  DO M=1,MAXM
                    BIN(M) = BIN(M)+X*TRM(M)
                    TRM(M) = TRM(M)*XK(M,IBETA)
                  ENDDO
                ENDDO
                DO M=1,MAXM
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)*BIN(M)
                ENDDO
              ELSEIF(IA.EQ.1) THEN
                DO M=1,MAXM
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)
                ENDDO
              ELSE
                DO M=1,MAXM
                  BETA(M) = BTA1(M)
                ENDDO
              ENDIF
C ***       END CONDITIONAL STATEMENT OVER IB VALUES
            ENDIF
C
            IF(IBETA.EQ.1) THEN
              DO M=1,MAXM
                B5(M,I1,I2) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXM
                B6(M,I1,I2) = BETA(M)
              ENDDO
            ENDIF
C
C         END LOOPS OVER IBETA AND INDICES I1, I2
          ENDDO
        ENDDO
      ENDDO
C
70    CONTINUE
C
C     ALL BETA INTEGRAL LISTS COMPLETE
C
      RETURN
      END
C
C
      SUBROUTINE ERI0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                    EEEEEEEE RRRRRRR  IIII  000000                    C
C                    EE       RR    RR  II  00    00                   C
C                    EE       RR    RR  II  00    00                   C
C                    EEEEEE   RR    RR  II  00    00                   C
C                    EE       RRRRRRR   II  00    00                   C
C                    EE       RR    RR  II  00    00                   C
C                    EEEEEEEE RR    RR IIII  000000                    C
C                                                                      C
C -------------------------------------------------------------------- C
C  ERI0 EVALUATES A DIRECT AND EXCHANGE BATCH OF ELECTRON REPULSION    C
C  INTEGRALS OF ALL COMPONENT LABEL COMBINATIONS L AND S IN THE ATOMIC C
C  SCF PROCEDURE FOR A USER-SPECIFIED HAMILTONIAN.                     C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C    RJLLLL(M,N) - DIRECT ERI OVERLAP {LL,LL}                          C
C    RJSSSS(M,N) - DIRECT ERI OVERLAP {SS,SS}                          C
C    RJLLSS(M,N) - DIRECT ERI OVERLAP {LL,SS}                          C
C    RJSSLL(M,N) - DIRECT ERI OVERLAP {SS,LL}                          C
C    RKLLLL(M,N) - EXCHANGE ERI OVERLAP {LL,LL}                        C
C    RKSSSS(M,N) - EXCHANGE ERI OVERLAP {SS,SS}                        C
C    RKSLSL(M,N) - EXCHANGE ERI OVERLAP {SL,SL}                        C
C -------------------------------------------------------------------- C
C    N = 1 - KQN(A) > 0, KQN(B) > 0 (TYPICAL LABEL 11)                 C
C    N = 2 - KQN(A) < 0, KQN(B) > 0 (TYPICAL LABEL 12)                 C
C    N = 3 - KQN(A) > 0, KQN(B) < 0 (TYPICAL LABEL 21)                 C
C    N = 4 - KQN(A) < 0, KQN(B) < 0 (TYPICAL LABEL 22)                 C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MNU=MKP+1,MAB=2*MNU+5)
C
      CHARACTER*4 HMLTN
C
      COMMON/ANGL/BK(MNU,4),DK(MNU,4),HK(MNU,4),FK(MNU,4),GM(MNU,4),
     &            NUS(MNU),NUNUM,NUMIN,NUMAX
      COMMON/BSIJ/EIJ(MAB),RNIJ(3),EI,EJ
      COMMON/BSKL/EK(MB2),EKL0(MB2),EIK(MB2,MAB),EJL(MB2,MAB),
     &            EL(MB2),EKL1(MB2),EKL(MB2,MAB),RNKL(MB2,4),
     &            B1(MB2,MAB,MAB),B2(MB2,MAB,MAB),B3(MB2,MAB,MAB),
     &            B4(MB2,MAB,MAB),B5(MB2,MAB,MAB),B6(MB2,MAB,MAB)
      COMMON/BSQN/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/CLRE/RKLLLL(MB2,4),RKSSSS(MB2,4),RKSLSL(MB2,4),
     &            RJLLLL(MB2,4),RJSSSS(MB2,4),RJLLSS(MB2,4),
     &            RJSSLL(MB2,4)
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
C
C     EMPTY COUNTER ARRAYS FOR DIRECT AND EXCHANGE INTEGRALS
      DO M=1,MAXM
        DO N=1,4
          RJLLLL(M,N) = 0.0D0
          RJLLSS(M,N) = 0.0D0
          RJSSLL(M,N) = 0.0D0
          RJSSSS(M,N) = 0.0D0
          RKLLLL(M,N) = 0.0D0
          RKSLSL(M,N) = 0.0D0
          RKSSSS(M,N) = 0.0D0
        ENDDO
      ENDDO
C
      IF(HMLTN.EQ.'NORL') THEN
C     NON-RELATIVISTIC HAMILTONIAN
C
        C5 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+5)
C
      ELSE
C     RELATIVISTIC HAMILTONIAN
C
C       PREPARE VALUES FOR UPCOMING CALCULATIONS
        C1 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+1)
        C3 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+3)
        C5 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+5)
        C7 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+7)
        C9 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+9)
C
        V1 = 1.0D0
        V2 = 2.0D0
        V4 = 4.0D0
        V8 = 8.0D0
        VS = 1.6D1
C
        F  = DFLOAT(2*LQNA+1)
        G  = DFLOAT(2*LQNB+1)
C
        F0G0 = 1.0D0
        F1G0 = F
        F0G1 = G
        F1G1 = F*G
        F2G0 = F*F
        F0G2 = G*G
        F2G1 = F*F*G
        F1G2 = F*G*G
        F2G2 = F*F*G*G
C
      ENDIF
C
C**********************************************************************C
C     AN (LQNA,LQNB) COMBINATION HAS 1, 2 OR 4 (KQNA,KQNB) SUB-BLOCKS  C
C     SMALL-COMPONENT CONTRIBUTIONS DEPEND ON KQN SYMMETRY TYPE.       C
C**********************************************************************C
C
C     INITIATE LOOP OVER K,L BASIS FUNCTIONS
      DO M=1,MAXM
C
C       MORE VALUE PREPARATION
        E0000 = 1.0D0
        E1000 = EI
        E0100 = EJ
        E0010 = EK(M)
        E0001 = EL(M)
        E1100 = EI*EJ
        E1010 = EI*EK(M)
        E1001 = EI*EL(M)
        E0110 = EJ*EK(M)
        E0101 = EJ*EL(M)
        E0011 = EK(M)*EL(M)
        E1110 = EI*EJ*EK(M)
        E1101 = EI*EJ*EL(M)
        E1011 = EI*EK(M)*EL(M)
        E0111 = EJ*EK(M)*EL(M)
        E1111 = EI*EJ*EK(M)*EL(M)
C
C**********************************************************************C
C       DIRECT INTEGRAL MATRICES: RJSSSS, RJLLSS, RJSSLL, RJLLLL       C
C**********************************************************************C
C
        IF(HMLTN.EQ.'NORL') THEN
C       NON-RELATIVISTIC HAMILTONIAN
C
          B22 = EIJ(4)*EKL(M,3)*B1(M,2,2) + EIJ(3)*EKL(M,4)*B2(M,2,2)
          RJLLLL(M,4) = C5*B22
C
        ELSE
C       RELATIVISTIC HAMILTONIAN
C
C         TEMPORARY STORAGE OF VALUES
C         BXY MARKS 'B' BETA COMBINATION AND 'X' ONTO EFFECTIVE LQN STORE
          B00 = EIJ(2)*EKL(M,1)*B1(M,1,1) + EIJ(1)*EKL(M,2)*B2(M,1,1)
          B02 = EIJ(2)*EKL(M,3)*B1(M,1,2) + EIJ(1)*EKL(M,4)*B2(M,2,1)
          B04 = EIJ(2)*EKL(M,5)*B1(M,1,3) + EIJ(1)*EKL(M,6)*B2(M,3,1)
          B20 = EIJ(4)*EKL(M,1)*B1(M,2,1) + EIJ(3)*EKL(M,2)*B2(M,1,2)
          B22 = EIJ(4)*EKL(M,3)*B1(M,2,2) + EIJ(3)*EKL(M,4)*B2(M,2,2)
          B24 = EIJ(4)*EKL(M,5)*B1(M,2,3) + EIJ(3)*EKL(M,6)*B2(M,3,2)
          B40 = EIJ(6)*EKL(M,1)*B1(M,3,1) + EIJ(5)*EKL(M,2)*B2(M,1,3)
          B42 = EIJ(6)*EKL(M,3)*B1(M,3,2) + EIJ(5)*EKL(M,4)*B2(M,2,3)
          B44 = EIJ(6)*EKL(M,5)*B1(M,3,3) + EIJ(5)*EKL(M,6)*B2(M,3,3)
C
C         FILL RJ ARRAYS FOR THIS LQNA,LQNB BLOCK
C
C >>>     LQNA  =  0 AND LQNB  =  0 (REQUIRED FOR ALL BLOCKS)
          RJLLLL(M,4) = V1*F0G0*E0000*C5*B22
          IF(HMLTN.EQ.'NORL') GOTO 101
          RJLLSS(M,4) = V4*F0G0*E0011*C7*B24
          RJSSLL(M,4) = V4*F0G0*E1100*C7*B42
          RJSSSS(M,4) = VS*F0G0*E1111*C9*B44
C
C >>>     LQNA =/= 0                (NEED KQNA > 0 BLOCK)
          IF(LQNA.EQ.0) GOTO 103
          RJLLLL(M,3) = RJLLLL(M,4)
          RJLLSS(M,3) = RJLLSS(M,4)
          RJSSLL(M,3) = V4*F0G0*E1100*C7*B42
     &                - V2*F1G0*E1000*C5*B22 - V2*F1G0*E0100*C5*B22
     &                + V1*F2G0*E0000*C3*B02
          RJSSSS(M,3) = VS*F0G0*E1111*C9*B44
     &                - V8*F1G0*E1011*C7*B24 - V8*F1G0*E0111*C7*B24
     &                + V4*F2G0*E0011*C5*B04
103       CONTINUE
C
C >>>                    LQNB =/= 0 (NEED KQNB > 0 BLOCK)
          IF(LQNB.EQ.0) GOTO 102
          RJLLLL(M,2) = RJLLLL(M,4)
          RJLLSS(M,2) = V4*F0G0*E0011*C7*B24 
     &                - V2*F0G1*E0010*C5*B22 - V2*F0G1*E0001*C5*B22
     &                + V1*F0G2*E0000*C3*B20
          RJSSLL(M,2) = RJSSLL(M,4)
          RJSSSS(M,2) = VS*F0G0*E1111*C9*B44
     &                - V8*F0G1*E1110*C7*B42 - V8*F0G1*E1101*C7*B42
     &                + V4*F0G2*E1100*C5*B40
102       CONTINUE
C
C >>>     LQNA =/= 0 AND LQNB =/= 0 (NEED KQNA,KQNB > 0 BLOCK)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 101
          RJLLLL(M,1) = RJLLLL(M,4)
          RJLLSS(M,1) = RJLLSS(M,2)
          RJSSLL(M,1) = RJSSLL(M,3)
          RJSSSS(M,1) = VS*F0G0*E1111*C9*B44
     &                - V8*F0G1*E1110*C7*B42 - V8*F0G1*E1101*C7*B42
     &                - V8*F1G0*E1011*C7*B24 - V8*F1G0*E0111*C7*B24
     &                + V4*F2G0*E0011*C5*B04 + V4*F0G2*E1100*C5*B40
     &                + V4*F1G1*E1001*C5*B22 + V4*F1G1*E0110*C5*B22
     &                + V4*F1G1*E0101*C5*B22 + V4*F1G1*E1010*C5*B22
     &                - V2*F1G2*E1000*C3*B20 - V2*F1G2*E0100*C3*B20
     &                - V2*F2G1*E0010*C3*B02 - V2*F2G1*E0001*C3*B02
     &                + V1*F2G2*E0000*C1*B00
101       CONTINUE
C
        ENDIF
C
C**********************************************************************C
C       EXCHANGE INTEGRAL MATRICES: RKSSSS, RKSLSL, RKLLLL             C
C**********************************************************************C
C
C       LOOP OVER ORDERS FOR BETA INTEGRALS
        DO IV=1,NUNUM
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(IV)
C
C         INDEX OFFSETS FOR THE EXPONENT LISTS
          IA = NUMAX+NU+1
          IB = NUMAX-NU+1
C
C         INDEX OFFSETS FOR BETA INTEGRAL ARRAYS
          NX = (-NUMIN+NU)/2
          NY = ( NUMAX-NU)/2
C
          IF(HMLTN.EQ.'NORL') THEN
C         NON-RELATIVISTIC HAMILTONIAN
C
            B22 = EIK(M,IA+3)*EJL(M,IB+2)*B3(M,NX+2,NY+2)
     &          + EIK(M,IB+2)*EJL(M,IA+3)*B4(M,NX+2,NY+2)
            RKLLLL(M,4) = RKLLLL(M,4) + BK(IV,4)*C5*B22
C
          ELSE
C         RELATIVISTIC HAMILTONIAN
C
C           TEMPORARY STORAGE OF VALUES
C           BXY MARKS 'B' BETA COMBINATION AND 'X' TO EFFECTIVE LQN STORE
C           ONTO WHICH NU IS ADDED OR SUBTRACTED.         
            B00 = EIK(M,IA+1)*EJL(M,IB  )*B3(M,NX+1,NY+1)
     &          + EIK(M,IB  )*EJL(M,IA+1)*B4(M,NX+1,NY+1)
            B02 = EIK(M,IA+1)*EJL(M,IB+2)*B3(M,NX+1,NY+2)
     &          + EIK(M,IB  )*EJL(M,IA+3)*B4(M,NX+2,NY+1)
            B04 = EIK(M,IA+1)*EJL(M,IB+4)*B3(M,NX+1,NY+3)
     &          + EIK(M,IB  )*EJL(M,IA+5)*B4(M,NX+3,NY+1)
            B20 = EIK(M,IA+3)*EJL(M,IB  )*B3(M,NX+2,NY+1)
     &          + EIK(M,IB+2)*EJL(M,IA+1)*B4(M,NX+1,NY+2)
            B22 = EIK(M,IA+3)*EJL(M,IB+2)*B3(M,NX+2,NY+2)
     &          + EIK(M,IB+2)*EJL(M,IA+3)*B4(M,NX+2,NY+2)
            B24 = EIK(M,IA+3)*EJL(M,IB+4)*B3(M,NX+2,NY+3)
     &          + EIK(M,IB+2)*EJL(M,IA+5)*B4(M,NX+3,NY+2)
            B40 = EIK(M,IA+5)*EJL(M,IB  )*B3(M,NX+3,NY+1)
     &          + EIK(M,IB+4)*EJL(M,IA+1)*B4(M,NX+1,NY+3)
            B42 = EIK(M,IA+5)*EJL(M,IB+2)*B3(M,NX+3,NY+2)
     &          + EIK(M,IB+4)*EJL(M,IA+3)*B4(M,NX+2,NY+3)
            B44 = EIK(M,IA+5)*EJL(M,IB+4)*B3(M,NX+3,NY+3)
     &          + EIK(M,IB+4)*EJL(M,IA+5)*B4(M,NX+3,NY+3)
C
C >>>       LQNA  =  0 AND LQNB  =  0 (REQUIRED FOR ALL BLOCKS)
            RKLLLL(M,4) = RKLLLL(M,4) + BK(IV,4)*V1*F0G0*E0000*C5*B22
            RKSLSL(M,4) = RKSLSL(M,4) + BK(IV,4)*V4*F0G0*E1010*C7*B42
            RKSSSS(M,4) = RKSSSS(M,4) + BK(IV,4)*VS*F0G0*E1111*C9*B44
C
C >>>       LQNA =/= 0                (NEED KQNA > 0 BLOCK)
            IF(LQNA.EQ.0) GOTO 203
            RKLL = V1*F0G0*E0000*C5*B22
            RKSL = V4*F0G0*E1010*C7*B42 - V2*F1G0*E0010*C5*B22
            RKSS = VS*F0G0*E1111*C9*B44 - V8*F1G0*E1011*C7*B42 
     &           - V8*F1G0*E0111*C7*B24 + V4*F2G0*E0011*C5*B22
            RKLLLL(M,3) = RKLLLL(M,3) + BK(IV,3)*RKLL
            RKSLSL(M,3) = RKSLSL(M,3) + BK(IV,3)*RKSL
            RKSSSS(M,3) = RKSSSS(M,3) + BK(IV,3)*RKSS
203         CONTINUE
C
C >>>                    LQNB =/= 0 (NEED KQNB > 0 BLOCK)
            IF(LQNB.EQ.0) GOTO 202
            RKLL = V1*F0G0*E0000*C5*B22
            RKSL = V4*F0G0*E1010*C7*B42 - V2*F0G1*E1000*C5*B22
            RKSS = VS*F0G0*E1111*C9*B44 - V8*F0G1*E1110*C7*B42 
     &           - V8*F0G1*E1101*C7*B24 + V4*F0G2*E1100*C5*B22
            RKLLLL(M,2) = RKLLLL(M,2) + BK(IV,2)*RKLL
            RKSLSL(M,2) = RKSLSL(M,2) + BK(IV,2)*RKSL
            RKSSSS(M,2) = RKSSSS(M,2) + BK(IV,2)*RKSS
202         CONTINUE
C
C >>>       LQNA =/= 0 AND LQNB =/= 0 (NEED KQNA,KQNB > 0 BLOCK)
            IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 201
            RKLL = V1*F0G0*E0000*C5*B22
            RKSL = V4*F0G0*E1010*C7*B42 - V2*F1G0*E0010*C5*B22 
     &           - V2*F0G1*E1000*C5*B22 + V1*F1G1*E0000*C3*B02
            RKSS = VS*F0G0*E1111*C9*B44
     &           - V8*F0G1*E1110*C7*B42 - V8*F0G1*E1101*C7*B24 
     &           - V8*F1G0*E1011*C7*B42 - V8*F1G0*E0111*C7*B24
     &           + V4*F2G0*E0011*C5*B22 + V4*F0G2*E1100*C5*B22 
     &           + V4*F1G1*E0110*C5*B22 + V4*F1G1*E1001*C5*B22
     &           + V4*F1G1*E1010*C5*B40 + V4*F1G1*E0101*C5*B04
     &           - V2*F2G1*E0010*C3*B20 - V2*F1G2*E1000*C3*B20
     &           - V2*F2G1*E0001*C3*B02 - V2*F1G2*E0100*C3*B02
     &           + V1*F2G2*E0000*C1*B00
            RKLLLL(M,1) = RKLLLL(M,1) + BK(IV,1)*RKLL
            RKSLSL(M,1) = RKSLSL(M,1) + BK(IV,1)*RKSL
            RKSSSS(M,1) = RKSSSS(M,1) + BK(IV,1)*RKSS 
201         CONTINUE
C
          ENDIF
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     NORMALISE ACCUMULATED INTEGRALS (DIRECT AND EXCHANGE)            C
C**********************************************************************C
C
      IF(HMLTN.EQ.'NORL') THEN
C     NON-RELATIVISTIC HAMILTONIAN
C
        DO M=1,MAXM
          T0LLLL = RNIJ(1)*RNKL(M,1)
          RJLLLL(M,4) = RJLLLL(M,4)*T0LLLL
          RKLLLL(M,4) = RKLLLL(M,4)*T0LLLL
        ENDDO
C
      ELSE
C     RELATIVISTIC HAMILTONIAN
C
        DO M=1,MAXM
          T0LLLL = RNIJ(1)*RNKL(M,1)
          T0LLSS = RNIJ(1)*RNKL(M,3)
          T0SLSL = RNIJ(2)*RNKL(M,2)
          T0SSLL = RNIJ(3)*RNKL(M,1)
          T0SSSS = RNIJ(3)*RNKL(M,3)
          DO N=1,4
            RJLLLL(M,N) = RJLLLL(M,N)*T0LLLL
            RJLLSS(M,N) = RJLLSS(M,N)*T0LLSS
            RJSSLL(M,N) = RJSSLL(M,N)*T0SSLL
            RJSSSS(M,N) = RJSSSS(M,N)*T0SSSS
            RKLLLL(M,N) = RKLLLL(M,N)*T0LLLL
            RKSLSL(M,N) = RKSLSL(M,N)*T0SLSL
            RKSSSS(M,N) = RKSSSS(M,N)*T0SSSS
          ENDDO
        ENDDO
C
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE BII0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                     BBBBBBB  IIII IIII  000000                       C
C                     BB    BB  II   II  00    00                      C
C                     BB    BB  II   II  00    00                      C
C                     BBBBBBB   II   II  00    00                      C
C                     BB    BB  II   II  00    00                      C
C                     BB    BB  II   II  00    00                      C
C                     BBBBBBB  IIII IIII  000000                       C
C                                                                      C
C -------------------------------------------------------------------- C
C  BII0 EVALUATES A DIRECT AND EXCHANGE BATCH OF BREIT INTERACTION     C
C  INTEGRALS OF ALL COMPONENT LABEL COMBINATIONS L AND S IN THE ATOMIC C
C  (RELATIVISTIC) SCF PROCEDURE.                                       C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C    RKLLSS(M,N) - EXCHANGE BII OVERLAP {LL,SS}                        C
C    RKSLSL(M,N) - EXCHANGE BII OVERLAP {SL,SL}                        C
C    RKSSLL(M,N) - EXCHANGE BII OVERLAP {SS,LL}                        C
C    RMSLSL(M,N) - SEMI-RANGE BII OVERLAP {SL,SL}                      C
C -------------------------------------------------------------------- C
C    N = 1 - KQN(A) > 0, KQN(B) > 0 (TYPICAL LABEL 11)                 C
C    N = 2 - KQN(A) < 0, KQN(B) > 0 (TYPICAL LABEL 12)                 C
C    N = 3 - KQN(A) > 0, KQN(B) < 0 (TYPICAL LABEL 21)                 C
C    N = 4 - KQN(A) < 0, KQN(B) < 0 (TYPICAL LABEL 22)                 C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MNU=MKP+1,MAB=2*MNU+5)
C
      COMMON/ANGL/BK(MNU,4),DK(MNU,4),HK(MNU,4),FK(MNU,4),GM(MNU,4),
     &            NUS(MNU),NUNUM,NUMIN,NUMAX
      COMMON/BSIJ/EIJ(MAB),RNIJ(3),EI,EJ
      COMMON/BSKL/EK(MB2),EKL0(MB2),EIK(MB2,MAB),EJL(MB2,MAB),
     &            EL(MB2),EKL1(MB2),EKL(MB2,MAB),RNKL(MB2,4),
     &            B1(MB2,MAB,MAB),B2(MB2,MAB,MAB),B3(MB2,MAB,MAB),
     &            B4(MB2,MAB,MAB),B5(MB2,MAB,MAB),B6(MB2,MAB,MAB)
      COMMON/BSQN/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/BTRE/RKLLSS(MB2,4),RKSLSL(MB2,4),RKSSLL(MB2,4),
     &            RMSLSL(MB2,4)
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
C
C     EMPTY COUNTER ARRAYS FOR DIRECT AND EXCHANGE INTEGRALS
      DO M=1,MAXM
        DO N=1,4
          RKLLSS(M,N) = 0.0D0
          RKSLSL(M,N) = 0.0D0
          RKSSLL(M,N) = 0.0D0
          RMSLSL(M,N) = 0.0D0
        ENDDO
      ENDDO
C
C     PREPARE VALUES FOR UPCOMING CALCULATIONS
      C1 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+1)
      C3 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+3)
      C5 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+5)
      C7 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+7)
      C9 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+9)
C
      V1 = 1.0D0
      V2 = 2.0D0
      V4 = 4.0D0
      V8 = 8.0D0
      VS = 1.6D1
C
      F  = DFLOAT(2*LQNA+1)
      G  = DFLOAT(2*LQNB+1)
C
      F0G0 = 1.0D0
      F1G0 = F
      F0G1 = G
      F1G1 = F*G
      F2G0 = F*F
      F0G2 = G*G
      F2G1 = F*F*G
      F1G2 = F*G*G
      F2G2 = F*F*G*G
C
C**********************************************************************C
C     AN (LQNA,LQNB) COMBINATION HAS 1, 2 OR 4 (KQNA,KQNB) SUB-BLOCKS  C
C     SMALL-COMPONENT CONTRIBUTIONS DEPEND ON KQN SYMMETRY TYPE.       C
C**********************************************************************C
C
C     INITIATE LOOP OVER K,L BASIS FUNCTIONS
      DO M=1,MAXM
C
C       MORE VALUE PREPARATION
        E0000 = 1.0D0
        E1000 = EI
        E0100 = EJ
        E0010 = EK(M)
        E0001 = EL(M)
        E1100 = EI*EJ
        E1010 = EI*EK(M)
        E1001 = EI*EL(M)
        E0110 = EJ*EK(M)
        E0101 = EJ*EL(M)
        E0011 = EK(M)*EL(M)
        E1110 = EI*EJ*EK(M)
        E1101 = EI*EJ*EL(M)
        E1011 = EI*EK(M)*EL(M)
        E0111 = EJ*EK(M)*EL(M)
        E1111 = EI*EJ*EK(M)*EL(M)
C
C**********************************************************************C
C       EXCHANGE INTEGRAL MATRICES: RKLLSS, RKSLSL, RKSSLL             C
C**********************************************************************C
C
C       LOOP OVER ORDERS FOR BETA INTEGRALS
        DO IV=1,NUNUM
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(IV)
C
C         INDEX OFFSETS FOR THE EXPONENT LISTS
          IA = NUMAX+NU+1
          IB = NUMAX-NU+1
C
C         INDEX OFFSETS FOR BETA INTEGRAL ARRAYS
          NX = (-NUMIN+NU)/2
          NY = ( NUMAX-NU)/2
C
C         TEMPORARY STORAGE OF VALUES
C         BXY MARKS 'B' BETA COMBINATION AND 'X' TO EFFECTIVE LQN STORE
C         ONTO WHICH NU IS ADDED OR SUBTRACTED.
C
C         CONVENTIONAL EXCHANGE TYPE INTEGRALS
          B02 = EIK(M,IA+1)*EJL(M,IB+2)*B3(M,NX+1,NY+2)
     &        + EIK(M,IB  )*EJL(M,IA+1)*B4(M,NX+2,NY+1)
          B22 = EIK(M,IA+3)*EJL(M,IB+2)*B3(M,NX+2,NY+2)
     &        + EIK(M,IB+2)*EJL(M,IA+3)*B4(M,NX+2,NY+2)
          B42 = EIK(M,IA+5)*EJL(M,IB+2)*B3(M,NX+3,NY+2)
     &        + EIK(M,IB+4)*EJL(M,IA+3)*B4(M,NX+2,NY+3)
C
C         BREIT EXCHANGE TYPE INTEGRALS (FOR LL AND SS BLOCKS)
          B11 = EIK(M,IA+2)*EJL(M,IB+1)*B5(M,NX+1,NY+1)
     &        + EIK(M,IB+1)*EJL(M,IA+2)*B6(M,NX+1,NY+1)
          B13 = EIK(M,IA+2)*EJL(M,IB+3)*B5(M,NX+1,NY+2)
     &        + EIK(M,IB+1)*EJL(M,IA+4)*B6(M,NX+2,NY+1)
          B31 = EIK(M,IA+4)*EJL(M,IB+1)*B5(M,NX+2,NY+1)
     &        + EIK(M,IB+3)*EJL(M,IA+2)*B6(M,NX+1,NY+2)
C         BREIT0 FIX: RKSSLL(M,4) IN HE_GEO IS EMPTY -> B33 HAS NOTHING
          B33 = EIK(M,IA+4)*EJL(M,IB+3)*B5(M,NX+2,NY+3)
     &        + EIK(M,IB+5)*EJL(M,IA+2)*B6(M,NX+3,NY+2)
C
C
C >>>     LQNA  =  0 AND LQNB  =  0 (REQUIRED FOR ALL BLOCKS)
          RKLLSS(M,4) = RKLLSS(M,4) + DK(IV,4)*V4*F0G0*E0011*C7*B33
          RKSLSL(M,4) = RKSLSL(M,4) + HK(IV,4)*V4*F0G0*E1010*C7*B42
          RKSSLL(M,4) = RKSSLL(M,4) + FK(IV,4)*V4*F0G0*E1100*C7*B33
C
C >>>     LQNA =/= 0                (NEED KQNA > 0 BLOCK)
          IF(LQNA.EQ.0) GOTO 203
          RKLL = V4*F0G0*E0011*C7*B33
          RKSL = V4*F0G0*E1010*C7*B42 - V2*F1G0*E0010*C5*B22
          RKSS = V4*F0G0*E1100*C7*B33 - V2*F1G0*E1000*C5*B31
     &         - V2*F1G0*E0100*C5*B13 + V1*F2G0*E0000*C3*B11
          RKLLSS(M,3) = RKLLSS(M,3) + DK(IV,3)*RKLL
          RKSLSL(M,3) = RKSLSL(M,3) + HK(IV,3)*RKSL
          RKSSLL(M,3) = RKSSLL(M,3) + FK(IV,3)*RKSS
203       CONTINUE
C
C >>>     LQNB =/= 0                (NEED KQNB > 0 BLOCK)
          IF(LQNB.EQ.0) GOTO 202
          RKLL = V4*F0G0*E0011*C7*B33 - V2*F0G1*E0010*C5*B31
     &         - V2*F0G1*E0001*C5*B13 + V1*F1G1*E0000*C3*B11
          RKSL = V4*F0G0*E1010*C7*B42 - V2*F0G1*E1000*C5*B22
          RKSS = V4*F0G0*E1100*C7*B33
          RKLLSS(M,2) = RKLLSS(M,2) + DK(IV,2)*RKLL
          RKSLSL(M,2) = RKSLSL(M,2) + HK(IV,2)*RKSL
          RKSSLL(M,2) = RKSSLL(M,2) + FK(IV,2)*RKSS
202       CONTINUE
C
C >>>     LQNA =/= 0 AND LQNB =/= 0 (NEED KQNA,KQNB > 0 BLOCK)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 201
          RKLL = V4*F0G0*E0011*C7*B33 - V2*F0G1*E0010*C5*B31
     &         - V2*F0G1*E0001*C5*B13 + V1*F1G1*E0000*C3*B11
          RKSL = V4*F0G0*E1010*C7*B42 - V2*F1G0*E0010*C5*B22
     &         - V2*F0G1*E1000*C5*B22 + V1*F1G1*E0000*C3*B02
          RKSS = V4*F0G0*E1100*C7*B33 - V2*F1G0*E1000*C5*B31
     &         - V2*F1G0*E0100*C5*B13 + V1*F2G0*E0000*C3*B11
          RKLLSS(M,1) = RKLLSS(M,1) + DK(IV,1)*RKLL
          RKSLSL(M,1) = RKSLSL(M,1) + HK(IV,1)*RKSL
          RKSSLL(M,1) = RKSSLL(M,1) + FK(IV,1)*RKSS
201       CONTINUE
C
        ENDDO
C
C**********************************************************************C
C       BREIT INTEGRAL MATRICES: RMSLSL                                C
C**********************************************************************C
C
C       LOOP OVER ORDERS FOR BETA INTEGRALS
        DO IV=1,NUNUM
          NU = NUS(IV)
          IR = NUS(NUNUM)+NU+2
          IS = NUS(NUNUM)-NU+2
          NX = (-NUS(1  )+NU+2)/2
          NY = ( NUS(NUNUM)-NU+2)/2
C
C         TEMPORARY STORAGE OF VALUES
C         BXY MARKS 'B' BETA COMBINATION AND 'X' TO EFFECTIVE LQN STORE
C         ONTO WHICH NU IS ADDED OR SUBTRACTED.
cc          B11 = EIK(M,IR+1)*EJL(M,IS  )*B5(M,NX+1,NY+1)
cc     &        - EIK(M,IS  )*EJL(M,IR+1)*B6(M,NX+1,NY+1)
cc          B13 = EIK(M,IR+1)*EJL(M,IS+2)*B5(M,NX+1,NY+2)
cc     &        - EIK(M,IS  )*EJL(M,IR+3)*B6(M,NX+2,NY+1)
cc          B31 = EIK(M,IR+3)*EJL(M,IS  )*B5(M,NX+2,NY+1)
cc     &        - EIK(M,IS+2)*EJL(M,IR+1)*B6(M,NX+1,NY+2)
cc          B33 = EIK(M,IR+3)*EJL(M,IS+2)*B5(M,NX+2,NY+2)
cc     &        - EIK(M,IS+2)*EJL(M,IR+3)*B6(M,NX+2,NY+2)
C
          B11 = EIK(M,IR+1)*EJL(M,IS  )*B5(M,NX  ,NY  )
          B13 = EIK(M,IR+1)*EJL(M,IS+2)*B5(M,NX  ,NY+1)
          B31 = EIK(M,IR+3)*EJL(M,IS  )*B5(M,NX+1,NY  )
          B33 = EIK(M,IR+3)*EJL(M,IS+2)*B5(M,NX+1,NY+1)
C
C >>>     LQNA  =  0 AND LQNB  =  0 (REQUIRED FOR ALL BLOCKS)
          RMSLSL(M,4) = RMSLSL(M,4) + GM(IV,4)*V4*F0G0*E1100*C7*B33
C
C >>>     LQNA =/= 0                (NEED KQNA > 0 BLOCK)
          IF(LQNA.EQ.0) GOTO 303
          RMSL = V4*F0G0*E1100*C7*B33 - V2*F1G0*E1000*C5*B31
     &         - V2*F1G0*E0100*C5*B13 + V1*F2G0*E0000*C3*B11
          RMSLSL(M,3) = RMSLSL(M,3) + GM(IV,3)*RMSL
303       CONTINUE
C
C >>>     LQNB =/= 0                (NEED KQNB > 0 BLOCK)
          IF(LQNB.EQ.0) GOTO 302
          RMSL = V4*F0G0*E1100*C7*B33
          RMSLSL(M,2) = RMSLSL(M,2) + GM(IV,2)*RMSL
302       CONTINUE
C
C >>>     LQNA =/= 0 AND LQNB =/= 0 (NEED KQNA,KQNB > 0 BLOCK)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 301
          RMSL = V4*F0G0*E1100*C7*B33 - V2*F1G0*E1000*C5*B31
     &         - V2*F1G0*E0100*C5*B13 + V1*F2G0*E0000*C3*B11
          RMSLSL(M,1) = RMSLSL(M,1) + GM(IV,1)*RMSL
301       CONTINUE
C
        ENDDO
      ENDDO
C      
C**********************************************************************C
C     NORMALISE ACCUMULATED INTEGRALS (DIRECT AND EXCHANGE)            C
C**********************************************************************C
C
      DO M=1,MAXM
        T0LLSS = RNIJ(1)*RNKL(M,3)
        T0SLSL = RNIJ(2)*RNKL(M,2)
        T0SSLL = RNIJ(3)*RNKL(M,1)
        DO N=1,4
          RKLLSS(M,N) = RKLLSS(M,N)*T0LLSS
          RKSLSL(M,N) = RKSLSL(M,N)*T0SLSL
          RKSSLL(M,N) = RKSSLL(M,N)*T0SSLL
          RMSLSL(M,N) = RMSLSL(M,N)*T0SLSL
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ANGCOUL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       AA    NN    NN  GGGGGG   CCCCCC   OOOOOO  UU    UU LL          C
C      AAAA   NNN   NN GG    GG CC    CC OO    OO UU    UU LL          C
C     AA  AA  NNNN  NN GG       CC       OO    OO UU    UU LL          C
C    AA    AA NN NN NN GG       CC       OO    OO UU    UU LL          C
C    AAAAAAAA NN  NNNN GG   GGG CC       OO    OO UU    UU LL          C
C    AA    AA NN   NNN GG    GG CC    CC OO    OO UU    UU LL          C
C    AA    AA NN    NN  GGGGGG   CCCCCC   OOOOOO   UUUUUU  LLLLLLLL    C
C                                                                      C
C -------------------------------------------------------------------- C
C  ANGCOUL EVALUATES THE ANGULAR COEFFICIENTS OF THE COULOMB           C
C  INTERACTIONS FOR CLOSED SHELLS IN THE (L1,L2) MANIFOLD.             C
C**********************************************************************C
      PARAMETER(MKP=9,MNU=MKP+1)
C
      CHARACTER*4 HMLTN
C
      COMMON/ANGL/BK(MNU,4),DK(MNU,4),HK(MNU,4),FK(MNU,4),GM(MNU,4),
     &            NUS(MNU),NUNUM,NUMIN,NUMAX
      COMMON/BSQN/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
C
C     CALCULATE KQNA AND 2*JQNA VALUES FROM LQNA
      KLA =-LQNA-1
      KRA = LQNA
      JLA = 2*IABS(KLA)-1
      JRA = 2*IABS(KRA)-1
C
C     CALCULATE KQNB AND 2*JQNB VALUES FROM LQNB
      KLB =-LQNB-1
      KRB = LQNB
      JLB = 2*IABS(KLB)-1
      JRB = 2*IABS(KRB)-1
C
C     GENERATE LIST OF FACTORIALS
      CALL FACTRLS
C
C     START AND END PARAMETERS FROM TRIANGLE RULE
      NUMIN = IABS(LQNA-LQNB)
      NUMAX = LQNA+LQNB+1
      NUNUM = 0
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      DO NU=NUMIN,NUMAX
C
C       TEST WHETHER 'LQNA+LQNB+NU' ODD OR EVEN
        LTEST = LQNA+LQNB+NU
        LEVEN = 2*(LTEST/2)
C
C       ANGULAR COEFFICIENTS OF ODD PARITY ARE ZERO
        IF(LTEST.NE.LEVEN) GOTO 14
C
C       ANGULAR COEFFICIENTS OF EVEN PARITY ARE NON-ZERO
        NUNUM       = NUNUM+1
        NUS(NUNUM)  = NU
C
        IF(HMLTN.EQ.'NORL') THEN
          BK(NUNUM,4) = 0.5D0*ABC000(LQNA,LQNB,NU)
        ELSE
          BK(NUNUM,1) = SYM3JSQ(JRA,JRB,NU)
          BK(NUNUM,2) = SYM3JSQ(JLA,JRB,NU)
          BK(NUNUM,3) = SYM3JSQ(JRA,JLB,NU)
          BK(NUNUM,4) = SYM3JSQ(JLA,JLB,NU)
        ENDIF
C
14      CONTINUE
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ANGBREIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    AA    NN    NN  GGGGGG  BBBBBBB  RRRRRRR  EEEEEEEE IIII TTTTTTTT  C
C   AAAA   NNN   NN GG    GG BB    BB RR    RR EE        II     TT     C
C  AA  AA  NNNN  NN GG       BB    BB RR    RR EE        II     TT     C
C AA    AA NN NN NN GG       BBBBBBB  RR    RR EEEEEE    II     TT     C
C AAAAAAAA NN  NNNN GG   GGG BB    BB RRRRRRR  EE        II     TT     C
C AA    AA NN   NNN GG    GG BB    BB RR    RR EE        II     TT     C
C AA    AA NN    NN  GGGGGG  BBBBBBB  RR    RR EEEEEEEE IIII    TT     C
C                                                                      C
C -------------------------------------------------------------------- C
C  ANGBREIT EVALUATES ANGULAR COEFFICIENTS OF THE ATOMIC CLOSED SHELL  C
C  BREIT INTERACTION FOR ALL (K1,K2) VALUES IN THE MANIFOLD (L1,L2).   C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C    NUNUM - NUMBER OF NU VALUES THAT SATISFY PARITY RESTRICTION RULE. C
C    NUMIN - MINIMUM NU VALUE IN THIS MANIFOLD.                        C
C    NUMAX - MAXIMUM NU VALUE IN THIS MANIFOLD.                        C
C    DK(MNU,4) - (K1,K2) COEFFICIENTS FOR BREIT K^{NU,LL,SS} INTEGRALS C
C    HK(MNU,4) - (K1,K2) COEFFICIENTS FOR BREIT K^{NU,SL,SL} INTEGRALS C
C    FK(MNU,4) - (K1,K2) COEFFICIENTS FOR BREIT K^{NU,SS,LL} INTEGRALS C
C    GM(MNU,4) - (K1,K2) COEFFICIENTS FOR BREIT M^{NU,SL,SL} INTEGRALS C
C**********************************************************************C
      PARAMETER(MKP=9,MNU=MKP+1)
C
      DIMENSION S(4,2)
C
      COMMON/ANGL/BK(MNU,4),DK(MNU,4),HK(MNU,4),FK(MNU,4),GM(MNU,4),
     &            NUS(MNU),NUNUM,NUMIN,NUMAX
      COMMON/BSQN/NFUNA,NFUNB,LQNA,LQNB,MAXM
C
C     INITIALISE COEFFICIENT ARRAYS
      DO NU=1,MNU
        DO N=1,4
          DK(NU,N) = 0.0D0
          HK(NU,N) = 0.0D0
          FK(NU,N) = 0.0D0
          GM(NU,N) = 0.0D0
        ENDDO
        NUS(NU) = 0
      ENDDO
      NUNUM = 0
C
C     CALCULATE KQNA AND 2*JQNA VALUES FROM LQNA
      KLA =-LQNA-1
      KRA = LQNA
      JLA = 2*IABS(KLA)-1
      JRA = 2*IABS(KRA)-1
C
C     CALCULATE KQNB AND 2*JQNB VALUES FROM LQNB
      KLB =-LQNB-1
      KRB = LQNB
      JLB = 2*IABS(KLB)-1
      JRB = 2*IABS(KRB)-1
C
C     GENERATE LIST OF FACTORIALS
      CALL FACTRLS
C
C**********************************************************************C
C     (4) JLA AND JLB               (REQUIRED FOR ALL BLOCKS)          C
C**********************************************************************C
C
C     START AND END PARAMETERS FROM TRIANGLE RULE
      NUMIN = IABS(JLA-JLB)/2
      NUMAX =     (JLA+JLB)/2
      NNU   = 1
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      DO NU=NUMIN,NUMAX
C
C       GENERATE SQUARES OF 3J-SYMBOLS FOR THIS NU VALUE
        VAL = SYM3JSQ(JLA,JLB,NU)
C
C       DETERMINE PARITY OF COMBINATION LQNA+LQNB+NU
        LTEST = LQNA+LQNB+NU
        LEVEN = 2*(LTEST/2)
C
C       STEP PARAMETER FOR ADDITIONS TO COEFFICIENTS
        NSTEP = 1
        IF(NU.LE.0) THEN
          NSTEP = 0
        ENDIF
C
        IF(LTEST.NE.LEVEN.AND.NU.NE.0) THEN
C       CASE 1: ODD-PARITY COEFFICIENTS
C
          RNU = DFLOAT(NU*(NU+1))
          SYM = DFLOAT((KLA+KLB)**2)/RNU
          NUS(NNU)  = NU
          DK(NNU,4) = DK(NNU,4) + SYM*VAL
          HK(NNU,4) = HK(NNU,4) + SYM*VAL
          FK(NNU,4) = FK(NNU,4) + SYM*VAL
C
        ELSEIF(LTEST.EQ.LEVEN) THEN
C       CASE 2: EVEN-PARITY COEFFICIENTS
C
          CALL BRCOEFF(S,KLA,KLB,NU)
C
          NUS(NNU)  = NU-1
          DK(NNU,4) = DK(NNU,4) - S(1,1)*VAL
          HK(NNU,4) = HK(NNU,4) - S(2,1)*VAL
          FK(NNU,4) = FK(NNU,4) - S(3,1)*VAL
          GM(NNU,4) = GM(NNU,4) - S(4,1)*VAL
C
          NNU = NNU + NSTEP
          NUS(NNU)  = NU+1
          DK(NNU,4) = DK(NNU,4) - S(1,2)*VAL
          HK(NNU,4) = HK(NNU,4) - S(2,2)*VAL
          FK(NNU,4) = FK(NNU,4) - S(3,2)*VAL
          GM(NNU,4) = GM(NNU,4) - S(4,2)*VAL
C
        ENDIF
      ENDDO
C
      NUNUM = MAX(NNU,NUNUM)
C
C**********************************************************************C
C     (3) JRA AND JLB                (NEED KQNA > 0 BLOCK)             C
C**********************************************************************C
C
      IF(LQNA.EQ.0) GOTO 203
C
C     START AND END PARAMETERS FROM TRIANGLE RULE
      NUMIN = IABS(JRA-JLB)/2
      NUMAX =     (JRA+JLB)/2
      NNU = 1
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      DO NU=NUMIN,NUMAX
C
C       GENERATE SQUARES OF 3J-SYMBOLS FOR THIS NU VALUE
        VAL = SYM3JSQ(JRA,JLB,NU)
C
C       DETERMINE PARITY OF COMBINATION LQNA+LQNB+NU
        LTEST = LQNA+LQNB+NU
        LEVEN = 2*(LTEST/2)
C
C       STEP PARAMETER FOR ADDITIONS TO COEFFICIENTS
        NSTEP = 1
        IF(NU.LE.0) THEN
          NSTEP = 0
        ENDIF
C
        IF(LTEST.NE.LEVEN.AND.NU.NE.0) THEN
C       CASE 1: ODD-PARITY COEFFICIENTS
          RNU = DFLOAT(NU*(NU+1))
          SYM = DFLOAT((KRA+KLB)**2)/RNU
          NUS(NNU)  = NU
          DK(NNU,3) = DK(NNU,3) + SYM*VAL
          HK(NNU,3) = HK(NNU,3) + SYM*VAL
          FK(NNU,3) = FK(NNU,3) + SYM*VAL
        ELSEIF(LTEST.EQ.LEVEN) THEN
C       CASE 2: EVEN-PARITY COEFFICIENTS
C
          CALL BRCOEFF(S,KRA,KLB,NU)
C
          NUS(NNU) = NU-1
          DK(NNU,3) = DK(NNU,3) - S(1,1)*VAL
          HK(NNU,3) = HK(NNU,3) - S(2,1)*VAL
          FK(NNU,3) = FK(NNU,3) - S(3,1)*VAL
          GM(NNU,3) = GM(NNU,3) - S(4,1)*VAL
C
          NNU = NNU + NSTEP
          NUS(NNU)  = NU+1
          DK(NNU,3) = DK(NNU,3) - S(1,2)*VAL
          HK(NNU,3) = HK(NNU,3) - S(2,2)*VAL
          FK(NNU,3) = FK(NNU,3) - S(3,2)*VAL
          GM(NNU,3) = GM(NNU,3) - S(4,2)*VAL
        ENDIF
      ENDDO
C
      NUNUM = MAX(NNU,NUNUM)
C
203   CONTINUE
C
C**********************************************************************C
C     (2) JLA AND JRB               (NEED KQNB > 0 BLOCK)              C
C**********************************************************************C
C
      IF(LQNB.EQ.0) GOTO 202
C
C     START AND END PARAMETERS FROM TRIANGLE RULE
      NUMIN = IABS(JLA-JRB)/2
      NUMAX =     (JLA+JRB)/2
      NNU = 1
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      DO NU=NUMIN,NUMAX
C
C       GENERATE SQUARES OF 3J-SYMBOLS FOR THIS NU VALUE
        VAL = SYM3JSQ(JLA,JRB,NU)
C
C       DETERMINE PARITY OF COMBINATION LQNA+LQNB+NU
        LTEST = LQNA+LQNB+NU
        LEVEN = 2*(LTEST/2)
C
C       STEP PARAMETER FOR ADDITIONS TO COEFFICIENTS
        NSTEP = 1
        IF(NU.LE.0) THEN
          NSTEP = 0
        ENDIF
C
        IF(LTEST.NE.LEVEN.AND.NU.NE.0) THEN
C       CASE 1: ODD-PARITY COEFFICIENTS
          RNU = DFLOAT(NU*(NU+1))
          SYM = DFLOAT((KLA+KRB)**2)/RNU
          NUS(NNU)  = NU
          DK(NNU,2) = DK(NNU,2) + SYM*VAL
          HK(NNU,2) = HK(NNU,2) + SYM*VAL
          FK(NNU,2) = FK(NNU,2) + SYM*VAL
        ELSEIF(LTEST.EQ.LEVEN) THEN
C       CASE 2: EVEN-PARITY COEFFICIENTS
C
          CALL BRCOEFF(S,KLA,KRB,NU)
C
          NUS(NNU)  = NU-1
          DK(NNU,2) = DK(NNU,2) - S(1,1)*VAL
          HK(NNU,2) = HK(NNU,2) - S(2,1)*VAL
          FK(NNU,2) = FK(NNU,2) - S(3,1)*VAL
          GM(NNU,2) = GM(NNU,2) - S(4,1)*VAL
C
          NNU  = NNU + NSTEP
          NUS(NNU)  = NU+1
          DK(NNU,2) = DK(NNU,2) - S(1,2)*VAL
          HK(NNU,2) = HK(NNU,2) - S(2,2)*VAL
          FK(NNU,2) = FK(NNU,2) - S(3,2)*VAL
          GM(NNU,2) = GM(NNU,2) - S(4,2)*VAL
        ENDIF
      ENDDO
C
      NUNUM = MAX(NNU,NUNUM)
C
202   CONTINUE
C
C**********************************************************************C
C     (1) JRA AND JRB               (NEED KQNA,KQNB > 0 BLOCK)         C
C**********************************************************************C
C
      IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 201
C
C     START AND END PARAMETERS FROM TRIANGLE RULE
      NUMIN = IABS(JRA-JRB)/2
      NUMAX =     (JRA+JRB)/2
      NNU = 1
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      DO NU=NUMIN,NUMAX
C
C       GENERATE SQUARES OF 3J-SYMBOLS FOR THIS NU VALUE
        VAL = SYM3JSQ(JRA,JRB,NU)
C
C       DETERMINE PARITY OF COMBINATION LQNA+LQNB+NU
        LTEST = LQNA+LQNB+NU
        LEVEN = 2*(LTEST/2)
C
C       STEP PARAMETER FOR ADDITIONS TO COEFFICIENTS
        NSTEP = 1
        IF(NU.LE.0) THEN
          NSTEP = 0
        ENDIF
C
        IF(LTEST.NE.LEVEN.AND.NU.NE.0) THEN
C       CASE 1: ODD-PARITY COEFFICIENTS
          RNU = DFLOAT(NU*(NU+1))
          SYM = DFLOAT((KRA+KRB)*(KRA+KRB))/RNU
          NUS(NNU)  = NU
          DK(NNU,1) = DK(NNU,1) + SYM*VAL
          HK(NNU,1) = HK(NNU,1) + SYM*VAL
          FK(NNU,1) = FK(NNU,1) + SYM*VAL
        ELSEIF(LTEST.EQ.LEVEN) THEN
C       CASE 2: EVEN-PARITY COEFFICIENTS
C
          CALL BRCOEFF(S,KRA,KRB,NU)
C
          NUS(NNU)  = NU-1
          DK(NNU,1) = DK(NNU,1) - S(1,1)*VAL
          HK(NNU,1) = HK(NNU,1) - S(2,1)*VAL
          FK(NNU,1) = FK(NNU,1) - S(3,1)*VAL
          GM(NNU,1) = GM(NNU,1) - S(4,1)*VAL
C
          NNU       = NNU + NSTEP
          NUS(NNU)  = NU+1
          DK(NNU,1) = DK(NNU,1) - S(1,2)*VAL
          HK(NNU,1) = HK(NNU,1) - S(2,2)*VAL
          FK(NNU,1) = FK(NNU,1) - S(3,2)*VAL
          GM(NNU,1) = GM(NNU,1) - S(4,2)*VAL

        ENDIF
      ENDDO
C
      NUNUM = MAX(NNU,NUNUM)
C
201   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE BRCOEFF(S,KQNA,KQNB,NU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    BBBBBBB  RRRRRRR   CCCCCC   OOOOOO  EEEEEEEE FFFFFFFF FFFFFFFF    C
C    BB    BB RR    RR CC    CC OO    OO EE       FF       FF          C
C    BB    BB RR    RR CC       OO    OO EE       FF       FF          C
C    BBBBBBB  RR    RR CC       OO    OO EEEEEE   FFFFFF   FFFFFF      C
C    BB    BB RRRRRRR  CC       OO    OO EE       FF       FF          C
C    BB    BB RR    RR CC    CC OO    OO EE       FF       FF          C
C    BBBBBBB  RR    RR  CCCCCC   OOOOOO  EEEEEEEE FF       FF          C
C                                                                      C
C -------------------------------------------------------------------- C
C                          SWIRLES MODULE 24:                          C
C                                                                      C
C  BRCOEFF EVALUATES THE INTERMEDIATE COEFFICIENTS OF THE BREIT        C
C  INTERACTION (TABLE 3 OF GRANT AND PYPER 1976).                      C
C**********************************************************************C
      DIMENSION S(4,2)

      NU1 = NU-1
      NU2 = NU+1
      RU  = DFLOAT(NU)
      KK  = KQNB-KQNA
      RK  = DFLOAT(KK)
C
      IF(NU1.GE.0) THEN
        BD = DFLOAT(         2 *(NU1+NU1+1))
        CD = DFLOAT((NU1+NU1+1)*(NU1+NU1+2))
        B1 = DFLOAT(NU1+2)/BD
        C1 =-DFLOAT(NU1-1)/CD
        S(1,1) =-(RU+RK)*(B1+(C1*RK))
        S(2,1) = (B1*RU)-(C1*RK*RK)
        S(3,1) =-(RU-RK)*(B1-(C1*RK))
        S(4,1) = RK*(B1-(C1*RU))
      ELSE
        S(1,1) = 0.0D0
        S(2,1) = 0.0D0
        S(3,1) = 0.0D0
        S(4,1) = 0.0D0
      ENDIF
C
      IF(NU2.GE.1) THEN
        BD = DFLOAT(2*    (NU2+NU2+1))
        CD = DFLOAT(2*NU2*(NU2+NU2+1))
        B2 = DFLOAT(NU2-1)/BD
        C2 = DFLOAT(NU2+2)/CD
        S(1,2) =-( B2+(C2*RK))*(RK-RU-1.0D0)
        S(2,2) =-((B2*(RU+1.0D0))+(C2*RK*RK))
        S(3,2) = ( B2-(C2*RK))*(RK+RU+1.0D0)
        S(4,2) =-RK*(B2+(C2*(RU+1.0D0)))
      ELSE
        S(1,2) = 0.0D0
        S(2,2) = 0.0D0
        S(3,2) = 0.0D0
        S(4,2) = 0.0D0
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION ABC000(L1,L2,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           AA    BBBBBBB   CCCCCC   000000   000000   000000          C
C          AAAA   BB    BB CC    CC 00    00 00    00 00    00         C
C         AA  AA  BB    BB CC       00    00 00    00 00    00         C
C        AA    AA BBBBBBB  CC       00    00 00    00 00    00         C
C        AAAAAAAA BB    BB CC       00    00 00    00 00    00         C
C        AA    AA BB    BB CC    CC 00    00 00    00 00    00         C
C        AA    AA BBBBBBB   CCCCCC   000000   000000   000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  ABC000 EVALUATES THE NON-RELATIVISTIC 3-J SYMBOL FOR ATOMIC COULOMB C
C  ANGULAR COEFFICIENT ROUTINES, TAKEN FROM BRINK AND SATCHLER.        C
C                                                                      C
C  L1,L2, AND K MUST BE EQUAL TO THE ACTUAL (INTEGER) ANGULAR MOMENTA  C
C  OF THE ELECTRON AND PHOTON.                                         C
C**********************************************************************C
      COMMON/FCTS/RFACT(21),SFACT(21)
C
C     TRIANGLE INEQUALITY RESTRICTIONS
      IF(K.LT.IABS(L1-L2).OR.K.GT.(L1+L2)) THEN
        ABC000 = 0.0D0
        RETURN
      ENDIF
      LLK = L1+L2+K
C
C     PARITY SELECTION RULE
      IF((LLK/2)*2.NE.LLK) THEN
        ABC000 = 0.0D0
        RETURN
      ENDIF
C
      RF1 = RFACT(  L1+L2-K + 1)
      RF2 = RFACT(- L1+L2+K + 1)
      RF3 = RFACT(  L1-L2+K + 1)
      RF4 = RFACT(  L1+L2+K + 2)
      RF5 = RFACT(( L1+L2+K)/2 + 1)
      RF6 = RFACT(( L1+L2-K)/2 + 1)
      RF7 = RFACT(( L1-L2+K)/2 + 1)
      RF8 = RFACT((-L1+L2+K)/2 + 1)
C
      T1  = RF1*RF2*RF3
      T2  = T1/RF4
      T3  = RF6*RF7*RF8
      T4  = RF5/T3
C
      ABC000 = T2*T4*T4
C
      RETURN
      END
C
C
      FUNCTION SYM3JSQ(J1,J2,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  SSSSSS  YY    YY MM       MM  333333       JJJJ  SSSSSS   QQQQQQ    C
C SS    SS YY    YY MMM     MMM 33    33       JJ  SS    SS QQ    QQ   C
C SS       YY    YY MMMM   MMMM       33       JJ  SS       QQ    QQ   C
C  SSSSSS   YY  YY  MM MM MM MM    3333        JJ   SSSSSS  QQ    QQ   C
C       SS   YYYY   MM  MMM  MM       33       JJ        SS QQ    QQ   C
C SS    SS    YY    MM   M   MM  33   33 JJ    JJ  SS    SS QQ    QQ   C
C  SSSSSS     YY    MM       MM   33333   JJJJJJ    SSSSSS   QQQQQQ Q  C
C                                                                      C
C -------------------------------------------------------------------- C
C  SYM3JSQ EVALUATES THE SQUARE OF A 3-J SYMBOL,   /  j   K   j' \^2   C
C  WHERE j = J1/2 AND j' = J2/2, FOR THE           \-1/2  0  1/2 /     C
C  COULOMB/BREIT ANGULAR COEFFICIENT ROUTINES.                         C
C**********************************************************************C
      COMMON/FCTS/RFACT(21),SFACT(21)
C
C     TRIANGLE RULE RESTRICTIONS
      IF(K.LT.IABS((J1-J2)/2).OR.K.GT.(J1+J2)/2) THEN
        SYM3JSQ = 0.0D0
        RETURN
      ELSEIF(J1.LE.0.OR.J2.LE.0) THEN
        SYM3JSQ = 0.0D0
        RETURN
      ENDIF
C
C     VARIABLE WHICH DEPENDS ON PARITY OF ARGUMENTS
      JJK = (J1+J2)/2 + K
      IF((JJK/2)*2.EQ.JJK) THEN
        M = K
      ELSE
        M = K+1
      ENDIF
C
      RN1 = RFACT(( J1+J2)/2 - K + 1)
      RN2 = RFACT((-J1+J2)/2 + K + 1)
      RN3 = RFACT(( J1-J2)/2 + K + 1)
      RN4 = SFACT(( J1+J2)/2 + M + 1)
      RD1 = DFLOAT(J1+1)
      RD2 = DFLOAT(J2+1)
      RD3 = RFACT(( J1+J2)/2 + K + 2)
      RD4 = SFACT(( J1+J2)/2 - M + 1)
      RD5 = SFACT(( J1-J2)/2 + M    )
      RD6 = SFACT((-J1+J2)/2 + M    )
      PHI = (-1.0D0)**((J2-(3*J1))/2+M)
C       
      RNUM  = RN1*RN2*RN3*(RN4)**2
      RDEN  = RD1*RD2*RD3*(RD4*RD5*RD6)**2
C
      SYM3JSQ = PHI*RNUM/RDEN
C
      RETURN
      END
C
C
      SUBROUTINE UEHLING0(UMAT,EXL,ZCRG,KQN,NFUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  UU    UU EEEEEEEE HH    HH LL      IIII NN    NN  GGGGGG   000000   C
C  UU    UU EE       HH    HH LL       II  NNN   NN GG    GG 00    00  C
C  UU    UU EE       HH    HH LL       II  NNNN  NN GG       00    00  C
C  UU    UU EEEEEE   HHHHHHHH LL       II  NN NN NN GG       00    00  C
C  UU    UU EE       HH    HH LL       II  NN  NNNN GG   GGG 00    00  C
C  UU    UU EE       HH    HH LL       II  NN   NNN GG    GG 00    00  C
C   UUUUUU  EEEEEEEE HH    HH LLLLLLL IIII NN    NN  GGGGGG   000000   C
C                                                                      C
C -------------------------------------------------------------------- C
C  UEHLING0 GENERATES ATOMIC UEHLING INTERACTION MATRIX FOR KQNA.      C
C**********************************************************************C
      PARAMETER(MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4 HMLTN
C
      DIMENSION RN(MBS*MBS,4),UMAT(2*MBS,2*MBS),EXL(MBS)
C
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
      DATA PI,ROOTPI/3.1415926535897932D0,1.7724538509055160D0/
C
C     DETERMINE THE LQN
      IF(KQN.LT.0) THEN
        LQN =-KQN-1
      ELSE
        LQN = KQN
      ENDIF
      RL = DFLOAT(LQN)
      G  = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NFUN,LQN)
C
C     COMPTON WAVELENGTH FOR COUPLING TO THE ELECTRON FIELD
      CMP = 1.0D0/CV
C
C     PRE-FACTOR COMMON TO ALL MATRIX ELEMENTS
      VCF = 6.0D0*CV*PI
      VCF =-ZCRG/VCF
C
C     CHARGE MOMENT FACTOR
      FCT = 1.0D0/(PNUC*CMP*CMP)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NFUN
        EI = EXL(IBAS)
        DO JBAS=1,NFUN
          M    = M+1
          EJ   = EXL(JBAS)
          EIJ  = EI+EJ
          EPR  = EI*EJ
          T52  = RL+2.5D0
          E52  = EIJ**T52
C
C         LL OVERLAP
          UMAT(IBAS,JBAS) = VCF*RN(M,1)*UEHINT0(LQN+1,EIJ)/EIJ
C
C         SS OVERLAP
          IF(HMLTN.EQ.'NORL') GOTO 50
C
C         SMALL COMPONENT BLOCK ADDRESSES
          KBAS = IBAS+NFUN
          LBAS = JBAS+NFUN
C
          VSA = 4.0D0*EPR*UEHINT0(LQN+2,EIJ)
          IF(KQN.GT.0) THEN
            VSN =-2.0D0*EIJ*G*UEHINT0(LQN+1,EIJ) 
     &                  + G*G*UEHINT0(LQN  ,EIJ)
          ELSE
            VSN = 0.0D0
          ENDIF
          UMAT(KBAS,LBAS) = VCF*RN(M,3)*(VSA+VSN)
C
50        CONTINUE
C         
        ENDDO
      ENDDO
C      
      RETURN
      END
C
C
      FUNCTION UEHINT0(N,ZETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       UU    UU EEEEEEEE HH    HH IIII NN    NN TTTTTTTT 000000       C
C       UU    UU EE       HH    HH  II  NNN   NN    TT   00    00      C
C       UU    UU EE       HH    HH  II  NNNN  NN    TT   00    00      C
C       UU    UU EEEEEE   HHHHHHHH  II  NN NN NN    TT   00    00      C
C       UU    UU EE       HH    HH  II  NN  NNNN    TT   00    00      C
C       UU    UU EE       HH    HH  II  NN   NNN    TT   00    00      C
C        UUUUUU  EEEEEEEE HH    HH IIII NN    NN    TT    000000       C
C                                                                      C
C -------------------------------------------------------------------- C
C  UEHINT0 CALCULATES AN ATOMIC UEHLING POTENTIAL OVERLAP INTEGRAL     C
C  FOR A LOCAL GAUSSIAN NUCLEAR CHARGE DISTRIBUTION.                   C
C -------------------------------------------------------------------- C
C  DFNOTE: AT THE MOMENT THIS IS JUST A COPY OF BNUCINT0.              C
C**********************************************************************C
      PARAMETER(MBS=26,MCT=6,MKP=9)
C
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     ROUTINE ONLY ALLOWS ODD PARAMETERS N
      IF(MOD(N,2).NE.1) THEN
        WRITE(6, *) 'In UEHINT0: order N must be odd. N = ',N
        WRITE(7, *) 'In UEHINT0: order N must be odd. N = ',N
      ENDIF
C
C     FACTORS NEEDED FOR ALL PARAMETERS N
      X   = ZETA/PNUC
      X5  = X*X*X*X*X
      T0  = PNUC+ZETA
      RAT = PNUC/T0
      TRM = 0.5D0*DSQRT(PNUC)/ZETA/DSQRT(T0)
      DO I=1,(N-1)/2
        TRM = 0.5D0*TRM*RAT/ZETA
      ENDDO
C
      IF(N.EQ.1) THEN
        TRM = TRM
      ELSEIF(N.EQ.3) THEN
        VA  = 2.0D0 + 3.0D0*X
        TRM = TRM*VA
      ELSEIF(N.EQ.5) THEN
        VA  = 8.0D0 + 20.0D0*X + 15.0D0*X*X
        TRM = TRM*VA
      ELSEIF(N.EQ.7) THEN
        VA  = 16.0D0 + 56.0D0*X + 70.0D0*X*X + 35.0D0*X*X*X
        TRM = 3.0D0*TRM*VA
      ELSEIF(N.EQ.9) THEN
        VA  = 128.0D0 + 576.0D0*X + 1008.0D0*X*X + 840.0D0*X*X*X 
        VB  = 315.0D0*X*X*X*X
        TRM = 3.0D0*TRM*(VA+VB)
      ELSEIF(N.EQ.11) THEN
        VA  = 256.0D0 + 1408.0D0*X + 3168.0D0*X*X + 3696.0D0*X*X*X 
        VB  = 2310.0D0*X*X*X*X + 693.0D0*X*X*X*X*X
        TRM = 15.0D0*TRM*(VA+VB)
      ELSEIF(N.EQ.13) THEN
        VA  = 1024.0D0 + 6656.0D0*X + 18304.0D0*X*X 
        VB  = 27456.0D0*X*X*X + 24024.0D0*X*X*X*X 
        VC  = 12012.0D0*X5 + 3003.0D0*X5*X
        TRM = 45.0D0*TRM*(VA+VB+VC)
      ELSEIF(N.EQ.15) THEN
        VA  = 2048.0D0 + 15360.0D0*X + 49920.0D0*X*X 
        VB  = 91520.0D0*X*X*X+ 102960.0D0*X*X*X*X + 72072.0D0*X5 
        VC  = 30030.0D0*X5*X + 6435.0D0*X5*X*X
        TRM = 315.0D0*TRM*(VA+VB+VC)
      ELSEIF(N.EQ.17) THEN
        VA  = 32768.0D0 + 278528.0D0*X + 1044480.0D0*X*X 
        VB  = 2263040.0D0*X*X*X + 3111680.0D0*X*X*X*X 
        VC  = 2800512.0D0*X5 + 1633632.0D0*X5*X + 583440.0D0*X5*X*X
        VD  = 109395.0D0*X5*X*X*X
        TRM = 315.0D0*TRM*(VA+VB+VC+VD)
      ELSEIF(N.EQ.19) THEN
        VA  = 65536.0D0 + 6222592.0D0*X + 2646016.0D0*X*X 
        VB  = 6615040.0D0*X*X*X + 10749440.0D0*X*X*X*X 
        VC  = 11824384.0D0*X5 + 8868288.0D0*X5*X + 4434144.0D0*X5*X*X 
        VD  = 1385670.0D0*X5*X*X*X + 230945.0D0*X5*X*X*X*X
        TRM = 2835.0D0*TRM*(VA+VB+VC+VD)
      ELSEIF(N.EQ.21) THEN
        VA  = 262144.0D0 + 2752512.0D0*X + 13074432.0D0*X*X 
        VB  = 37044224.0D0*X*X*X + 69457920.0D0*X*X*X*X 
        VC  = 90295296.0D0*X5 + 82770688.0D0*X5*X 
        VD  = 53209728.0D0*X5*X*X + 23279256.0D0*X5*X*X*X 
        VE  = 6466460.0D0*X5*X*X*X*X + 969969.0D0*X5*X5
        TRM = 14175.0D0*TRM*(VA+VB+VC+VD+VE)
      ELSE 
        WRITE(6, *) 'In UEHINT0: order N too large. N = ',N
        WRITE(7, *) 'In UEHINT0: order N too large. N = ',N
      ENDIF
C
C     TRANSFER DATA TO UEHINT0
      UEHINT0 = TRM
C
      RETURN
      END
C
C
      SUBROUTINE RNORM0(RN,EXL,NFUN,LQN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        RRRRRRR  NN    NN  OOOOOO  RRRRRRR  MM     MM  000000         C
C        RR    RR NNN   NN OO    OO RR    RR MMM   MMM 00    00        C
C        RR    RR NNNN  NN OO    OO RR    RR MMMM MMMM 00    00        C
C        RR    RR NN NN NN OO    OO RR    RR MM MMM MM 00    00        C
C        RRRRRRR  NN  NNNN OO    OO RRRRRRR  MM  M  MM 00    00        C
C        RR    RR NN   NNN OO    OO RR    RR MM     MM 00    00        C
C        RR    RR NN    NN  OOOOOO  RR    RR MM     MM  000000         C
C                                                                      C
C -------------------------------------------------------------------- C
C  RNORM0 EVALUATES NORMALISATION CONSTANTS OF ALL VARIETIES.          C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS)
C
      DIMENSION RN(MB2,4),EXL(MBS),RNL(MBS),RNS(MBS)
C
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
C
      DATA TWOLOG/6.93147180559945309D-1/
C
      RLQN = DFLOAT(LQN)
      G1   = TWOLOG - GAMLOG(2*LQN+3)
      G2   = TWOLOG - GAMLOG(2*LQN+5)
      R1   = RLQN + 1.5D0
      R2   = RLQN + 0.5D0
      DO IBAS=1,NFUN
        ELOG      = DLOG(2.0D0*EXL(IBAS))
        RNL(IBAS) = DEXP(0.5D0*(G1+R1*ELOG))
        RNS(IBAS) = DEXP(0.5D0*(G2+R2*ELOG))
      ENDDO
C
C     RN(M,1) ARE THE LL NORMALISATION CONSTANTS
C     RN(M,2) ARE THE SL NORMALISATION CONSTANTS
C     RN(M,3) ARE THE SS NORMALISATION CONSTANTS
C     RN(M,4) ARE THE LS NORMALISATION CONSTANTS
C
      M = 0
      DO IBAS=1,NFUN
        DO JBAS=1,NFUN
          M = M+1
          RN(M,1) = RNL(IBAS)*RNL(JBAS)
          RN(M,2) = RNS(IBAS)*RNL(JBAS)
          RN(M,3) = RNS(IBAS)*RNS(JBAS)
          RN(M,4) = RNL(IBAS)*RNS(JBAS)
        ENDDO
      ENDDO
C      
      RETURN
      END
C
C
      SUBROUTINE GAMGEN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        GGGGGG     AA    MM       MM  GGGGGG  EEEEEEEE NN    NN       C
C       GG    GG   AAAA   MMM     MMM GG    GG EE       NNN   NN       C
C       GG        AA  AA  MMMM   MMMM GG       EE       NNNN  NN       C
C       GG       AA    AA MM MM MM MM GG       EEEEEE   NN NN NN       C
C       GG   GGG AAAAAAAA MM  MMM  MM GG   GGG EE       NN  NNNN       C
C       GG    GG AA    AA MM   M   MM GG    GG EE       NN   NNN       C
C        GGGGGG  AA    AA MM       MM  GGGGGG  EEEEEEEE NN    NN       C
C                                                                      C
C -------------------------------------------------------------------- C
C  GAMGEN EVALUATES INTEGER/HALF-INTEGER GAMMA VALUES AND THEIR LOGS.  C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C    GAMLOG(N) = DLOG(GAMMA(N/2))                                      C
C    GAMHLF(N) = GAMMA(N/2)                                            C
C**********************************************************************C 
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
C
      DATA ROOTPI,RTPILG/1.7724538509055160D0,5.7236494292470009D-1/
C
C     STARTING VALUES
      GAMLOG(1) = RTPILG
      GAMLOG(2) = 0.0D0
      GAMHLF(1) = ROOTPI
      GAMHLF(2) = 1.0D0
C
C     SEED VALUES FOR INCREMENT
      F1 = 0.5D0
      F2 = 1.0D0
C
C     FILL TABLE VALUES
      DO N=4,50,2
        GAMLOG(N-1) = GAMLOG(N-3) + DLOG(F1)
        GAMLOG(N  ) = GAMLOG(N-2) + DLOG(F2)
        GAMHLF(N-1) = GAMHLF(N-3)*F1
        GAMHLF(N  ) = GAMHLF(N-2)*F2
        F1 = F1 + 1.0D0
        F2 = F2 + 1.0D0
      ENDDO
C
      RETURN
      END
C
C
       SUBROUTINE FACTRLS
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     FFFFFFFF   AA     CCCCCC TTTTTTTT RRRRRRR  LL       SSSSSS       C
C     FF        AAAA   CC    CC   TT    RR    RR LL      SS    SS      C
C     FF       AA  AA  CC         TT    RR    RR LL      SS            C
C     FFFFFF  AA    AA CC         TT    RR    RR LL       SSSSSS       C
C     FF      AAAAAAAA CC         TT    RRRRRRR  LL            SS      C
C     FF      AA    AA CC    CC   TT    RR    RR LL      SS    SS      C
C     FF      AA    AA  CCCCCC    TT    RR    RR LLLLLLLL SSSSSS       C
C                                                                      C
C -------------------------------------------------------------------- C
C  FACTRLS GENERATES A SET OF N! AND N!! AS REAL NUMBERS.              C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C    RFACT - REGULAR FACTORIALS, N!                                    C
C    SFACT - SEMI-FACTORIALS, N!!                                      C
C**********************************************************************C
C
      COMMON/FCTS/RFACT(21),SFACT(21)
C
      RFACT(1) = 1.0D0
      RFACT(2) = 1.0D0
      SFACT(1) = 1.0D0
      SFACT(2) = 1.0D0
      DO I=3,21
        RNUMBER  = DFLOAT(I-1)
        RFACT(I) = RNUMBER*RFACT(I-1)
        SFACT(I) = RNUMBER*SFACT(I-2)
      ENDDO
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C   [5] MOLECULAR HARTREE-FOCK: MANY-CENTER SCF CALCULATIONS.          C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] HFSCF: MAIN ROUTINE FOR MOLECULAR HARTREE-FOCK SCF PROCEDURE.  C
C   [B] OVRLP: CONSTRUCTS THE ONE-ELECTRON OVERLAP MATRIX.             C
C   [C] ONEEL: ONE-ELECTRON MULTI-CENTER MATRIX OF INTEGRALS.          C
C   [D] UEHLING: CONSTRUCTS MULTI-CENTER UEHLING MATRIX ELEMENTS.      C
C   [E] COULOMB: MATRIX REP OF MEAN-FIELD COULOMB INTERACTION.         C
C   [F] ERI: GENERATES A BLOCK OF ELECTRON REPULSION INTEGRALS.        C
C   [G] BREIT: MATRIX REP OF MEAN-FIELD BREIT INTERACTION.             C
C   [H] BII: GENERATES A BLOCK OF BREIT INTERACTION INTEGRALS.         C
C   [I] RDASUM: SUM OF MAGNITUDES OF A REAL VECTOR.                    C
C   [J] CDASUM: SUM OF MAGNITUDES OF A COMPLEX VECTOR.                 C
C   [K] NCART: RETURNS THE CARTESIAN INDEX FROM A LOOP INDEX.          C
C**********************************************************************C
C
C
      SUBROUTINE HFSCF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C              HH    HH FFFFFFFF SSSSSS   CCCCCC  FFFFFFFF             C
C              HH    HH FF      SS    SS CC    CC FF                   C
C              HH    HH FF      SS       CC       FF                   C
C              HHHHHHHH FFFFFF   SSSSSS  CC       FFFFFF               C
C              HH    HH FF            SS CC       FF                   C
C              HH    HH FF      SS    SS CC    CC FF                   C
C              HH    HH FF       SSSSSS   CCCCCC  FF                   C
C                                                                      C
C -------------------------------------------------------------------- C
C  HFSCF PERFORMS A SINGLE-DETERMINANT ITERATIVE SELF-CONSISTENT FIELD C
C  PROCEDURE OVER THE USER-SPECIFIED HAMILTONIAN.                      C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9,LWK=64*MDM,MIT=200)
C
      CHARACTER*2  ELMNT(120)
      CHARACTER*4  HMLTN
      CHARACTER*16 HMS
      CHARACTER*20 STAMP
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 TESTER
C
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 DTMP(MDM,MDM),OTMP(MDM,MDM),C(MDM,MDM)
      COMPLEX*16 WORK(LWK)
C
      DIMENSION RWORK(3*MDM),ESAV(0:MIT),DNRM(MIT),WEDN(MIT)
      DIMENSION NMLEV(4),TMLEV(4)
      
      DIMENSION ARRAY(MDM,MDM)
C
      COMMON/ATOM/ELMNT
      COMMON/COEF/C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/EIGN/EIGEN(MDM)
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EMDR,EMXC,EUEH
      COMMON/FILL/NCNF(MCT,MKP,MKP+1),NLVL(MCT,MKP),IFILL(MCT)
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/LSHF/SHLEV(3),SHLV
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/QNMS/ICNBAS(MDM),KQNBAS(MDM),MQNBAS(MDM)
      COMMON/SHLL/ALPHA,BETA,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN,NOELEC
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
C
C     PRINT A BIG SECTION HEADER
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) REPEAT(' ',24),'MOLECULAR HARTREE-FOCK SCF'
      WRITE(7, *) REPEAT(' ',24),'MOLECULAR HARTREE-FOCK SCF'
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C     RECORD TIME AT START OF MOLECULAR SCF CALCULATION
      CALL CPU_TIME(TDM1)
C
C     IF READING IN PREVIOUS SOLUTION, CALCULATE MOLECULAR DENSITY
      IF(INEW.EQ.1) THEN
        CALL DENSTY
      ENDIF
C
C     EQ-COEFFICIENTS IN COMMON ARRAYS
      IF(IEQS.EQ.1) THEN
        CALL EQFILE
      ENDIF
C
C     TOLERANCE VALUE FOR VANISHING MATRIX ELEMENTS
      TOL = 1.0D-10
C
C     PARAMETERS FOR COMPLETING STAGES
      DEPS  = 1.0D-10
      EEPS  = 1.0D-12
      TRSH3 = 5.0D-11
      TRSH2 = 1.0D-08
C
C     INITIALISE ENERGY NORM STORAGE
      IF(INEW.EQ.0) THEN
        ESAV(0) = ETOT
      ELSE
        ESAV(0) = 1.0D0
      ENDIF
C
C     INITIALISE INTEGRAL INCLUSION LEVEL VALUES
      DO N=1,4
        NMLEV(N) = 0
        TMLEV(N) = 0.0D0
      ENDDO
C
C     INITIALISE ARRAYS AND TEMPORARY DENSITY MATRIX
      DO I=1,NDIM
        DO J=1,NDIM
          OVAP(I,J) = DCMPLX(0.0D0,0.0D0)
          HNUC(I,J) = DCMPLX(0.0D0,0.0D0)
          HKIN(I,J) = DCMPLX(0.0D0,0.0D0)
          VUEH(I,J) = DCMPLX(0.0D0,0.0D0)
          GDIR(I,J) = DCMPLX(0.0D0,0.0D0)
          GXCH(I,J) = DCMPLX(0.0D0,0.0D0)
          QDIR(I,J) = DCMPLX(0.0D0,0.0D0)
          QXCH(I,J) = DCMPLX(0.0D0,0.0D0)
          BDIR(I,J) = DCMPLX(0.0D0,0.0D0)
          BXCH(I,J) = DCMPLX(0.0D0,0.0D0)
          FOCK(I,J) = DCMPLX(0.0D0,0.0D0)
          DTMP(I,J) = DENT(I,J)
        ENDDO
      ENDDO
C
C     INITIALISE SCF TIME COUNTERS
      THMX = 0.0D0
      TGMX = 0.0D0
      TBMX = 0.0D0
      TEIG = 0.0D0
C
C     INITIALISE EQ-COEFF AND R-INT TIME COUNTERS
      TELL = 0.0D0
      TESS = 0.0D0
      TELS = 0.0D0
      TRLL = 0.0D0
      TRSS = 0.0D0
      TRLS = 0.0D0
      TRBR = 0.0D0
C
C     PRINT THE FIRST INTEGRAL INCLUSION LEVEL
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('*',72)
      WRITE(7, *) REPEAT('*',72)
      IF(ILEV.EQ.1) THEN
        WRITE(6,30)
        WRITE(7,30)
      ELSEIF(ILEV.EQ.2) THEN
        IF(HMLTN.NE.'DHFB'.AND.HMLTN.NE.'DHFQ') THEN
          WRITE(6,31)
          WRITE(7,31)
        ELSE
          WRITE(6,32)
          WRITE(7,32)
        ENDIF
      ELSEIF(ILEV.EQ.3) THEN
        WRITE(6,33)
        WRITE(7,33)     
      ENDIF
      WRITE(6, *) REPEAT('*',72)
      WRITE(7, *) REPEAT('*',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C       BREIT0 FIX
CC       GENERATE MATRIX REP OF BREIT INTERACTION
C        CALL BREIT
C
C        OPEN(UNIT=8,FILE='Molecular_BREIT.dat',STATUS='UNKNOWN')
C        REWIND(UNIT=8)
C        DO I=1,NDIM
C          WRITE(8, *) (dreal(bdir(i,j)-bxch(i,j)),J=1,NDIM)
C        ENDDO
C        CLOSE(UNIT=8)
C        STOP
C
C     SELF-CONSISTENT ITERATION PROCEDURE
      DO ITER=1,MIT
C
C       TIME AT START OF ITERATION
        CALL CPU_TIME(TMIT)
C
C       GENERATE ONE-BODY MATRIX AND OVERLAP MATRIX
        IF(ITER.EQ.1) THEN
          CALL OVRLP
          CALL ONEEL         
        ENDIF
        CALL CPU_TIME(T1EL)
C
C       GENERATE UEHLING INTERACTION MATRIX
        IF(HMLTN.EQ.'DHFQ') THEN
          IF(ITER.EQ.1) THEN
            CALL UEHLING
          ENDIF
        ENDIF
C
C       GENERATE MEAN-FIELD CLOSED- AND OPEN-SHELL COULOMB MATRIX
        IF(HMLTN.NE.'BARE'.OR.NOELEC.GT.1) THEN
          IF(NOELEC.GT.1) THEN
            CALL COULOMB
          ENDIF
        ENDIF
        CALL CPU_TIME(T2CL)
C
C       GENERATE MEAN-FIELD BREIT MATRIX
        IF(HMLTN.EQ.'DHFB'.OR.HMLTN.EQ.'DHFQ') THEN
          IF(NOELEC.GT.1) THEN
            CALL BREIT
          ENDIF
        ENDIF
        CALL CPU_TIME(T2BT)
C
C       CONSTRUCT FOCK MATRIX FROM ONE- AND TWO-BODY INTERACTIONS
        DO I=1,NDIM
          DO J=1,NDIM
            FOCK(I,J) = HNUC(I,J) + HKIN(I,J) + GDIR(I,J) - GXCH(I,J)
     &                - QDIR(I,J) + QXCH(I,J) + BDIR(I,J) - BXCH(I,J)
     &                + VUEH(I,J)
          ENDDO
        ENDDO
C
C       UPDATE MOLECULAR ENERGIES (BASED ON *PREVIOUS* DENSITY D^{N-1})
        CALL ENERGIES
C
C       SEARCH FOR MATRIX SPARSITY
        DO I=1,NDIM
          DO J=1,NDIM
C
            X = DREAL(FOCK(I,J))
            Y = DIMAG(FOCK(I,J))
C
C           ELIMINATE ANY VANISHINGLY SMALL FOCK MATRIX ELEMENTS
            IF(DABS(X).LT.TOL) THEN
              X = 0.0D0
            ENDIF
            IF(DABS(Y).LT.TOL) THEN
              Y = 0.0D0
            ENDIF
C
C           ALSO ELMINATE ALL DIAGONAL IMAGINARY FOCK MATRIX ELEMENTS
            IF(I.EQ.J) THEN
              Y = 0.0D0
            ENDIF
C
C           TRANSFER ELEMENT BACK TO FOCK MATRIX
            FOCK(I,J) = DCMPLX(X,Y)
C
          ENDDO
        ENDDO
C
C       LEVEL-SHIFT THE VIRTUAL SPACE TO MAKE ORBITALS LESS ACCESSIBLE
        IF(INEW.NE.0.OR.ITER.NE.1) THEN
          CALL LEVSHFT(SHLV)
        ENDIF
        NMLEV(ILEV) = NMLEV(ILEV)+1
C
C       SAVE OVERLAP MATRIX IN TEMPORARY MATRIX (ZHEGV OVERWRITES IT)
        DO I=1,NDIM
          DO J=1,NDIM
            OTMP(I,J) = OVAP(I,J)
          ENDDO
        ENDDO
C
C       DIAGONALISE FOCK MATRIX (REQUIRES LAPACK LIBRARY)
        CALL ZHEGV(1,'V','L',NDIM,FOCK,MDM,OVAP,MDM,
     &                                       EIGEN,WORK,LWK,RWORK,INFO)
        IF(INFO.NE.0) THEN
          WRITE(6, *) 'In HFSCF: eigenvalue solver ZHEGV failed.',INFO
          WRITE(7, *) 'In HFSCF: eigenvalue solver ZHEGV failed.',INFO
        ENDIF
        CALL CPU_TIME(TDGN)
C
C       TRANSFER EIGENVECTORS TO THE C ARRAY AND RESTORE OVAP ARRAY
        DO J=1,NDIM
          DO I=1,NDIM
            C(I,J)    = FOCK(I,J)
            OVAP(I,J) = OTMP(I,J)
          ENDDO
        ENDDO
C
C       DEDUCT LEVEL SHIFT VALUE FROM VIRTUAL ORBITAL EIGENVALUES
        IF(ITER.NE.1) THEN
          DO IVIR=NSHIFT+NOCC+1,NDIM
            EIGEN(IVIR) = EIGEN(IVIR) - SHLV
          ENDDO
        ENDIF
C
C       WRITE EIGENVECTORS TO OUTPUT FILE
        OPEN(UNIT=8,FILE=TRIM(WFNFL),STATUS='UNKNOWN')
        REWIND(UNIT=8)
        DO I=1,NDIM
          WRITE(8, *) EIGEN(I),(C(J,I),J=1,NDIM)
        ENDDO
        CLOSE(UNIT=8)
C
C       UPDATE THE DENSITY MATRIX
        CALL DENSTY
C
C       DENSITY DIFFERENCE NORM CALCULATION
        DNRM(ITER) = 0.0D0
        DO J=1,NDIM
          DO I=1,NDIM
            TMP        = ABS(DENT(I,J) - DTMP(I,J))
            DNRM(ITER) = DNRM(ITER) + TMP*TMP
            DTMP(I,J)  = DENT(I,J)
          ENDDO
        ENDDO
        RDM2 = DFLOAT(NDIM*NDIM)
        DNRM(ITER) = DSQRT(DNRM(ITER))/RDM2
C
C       IF DNRM IS SMALL ENOUGH, REDUCE REQUIREMENTS TO ENTER STAGE 3
        IF(DNRM(ITER).LE.1.0D-09) THEN
          TRSH3 = 1.0D2*DNRM(ITER)
        ENDIF
C      
C       WEIGHTED ENERGY DIFFERENCE NORM, WEDN
        DOFF = DABS(ETOT)+1.0D0
        WEDN(ITER) = DABS(ESAV(ITER-1)-ETOT)/DOFF
        ESAV(ITER) = ETOT
C
C       UPDATE TIME COUNTERS
        THMX = THMX + T1EL - TMIT
        TGMX = TGMX + T2CL - T1EL
        TBMX = TBMX + T2BT - T2CL
        TEIG = TEIG + TDGN - T2BT
        TMLEV(ILEV) = TMLEV(ILEV)+TDGN-TMIT
C
C       DATE AND TIME AT END OF ITERATION
        CALL TIMENOW(STAMP)
C
C       HEADER FOR ITERATION SUMMARY
20      FORMAT(27X,A,1X,I3)
        WRITE(6, *) ' '
        WRITE(7, *) ' '
        WRITE(6,20) 'Iteration number',ITER
        WRITE(7,20) 'Iteration number',ITER
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
C
C       PRODUCE SPECTRUM SUMMARY AND INCLUDE FIRST 6 VIRTUAL STATES
        LF = LEN(TRIM(MOLCL))+LEN(TRIM(HMLTN))
        LF = 23-LF/2
21      FORMAT(1X,A,'Molecular spectrum for ',A,' (',A,')')
        WRITE(6,21) REPEAT(' ',LF),TRIM(MOLCL),TRIM(HMLTN)
        WRITE(7,21) REPEAT(' ',LF),TRIM(MOLCL),TRIM(HMLTN)
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        CALL SPECTRM(NOCC,6)
C
C       MOLECULAR ENERGIES
22      FORMAT(1X,A,22X,F21.12)
        WRITE(6, *) REPEAT(' ',72)
        WRITE(7, *) REPEAT(' ',72)
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6, *) REPEAT(' ',19),'Molecular energies (Hartree units)'
        WRITE(7, *) REPEAT(' ',19),'Molecular energies (Hartree units)'
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6, *) 'Source',REPEAT(' ',60),'Energy'
        WRITE(7, *) 'Source',REPEAT(' ',60),'Energy'
        WRITE(6, *) REPEAT('-',72)
        WRITE(7, *) REPEAT('-',72)
        WRITE(6,22) 'Nucleus-nucleus           (N)',ENUC
        WRITE(7,22) 'Nucleus-nucleus           (N)',ENUC
        WRITE(6,22) 'Electron-nucleus          (V)',EHNC
        WRITE(7,22) 'Electron-nucleus          (V)',EHNC
        WRITE(6,22) 'Electron kinetic          (T)',EHKN
        WRITE(7,22) 'Electron kinetic          (T)',EHKN
        IF(HMLTN.NE.'DHFQ') GOTO 80
        WRITE(6,22) 'Uehling interaction       (U)',EUEH
        WRITE(7,22) 'Uehling interaction       (U)',EUEH
80      CONTINUE
        IF(HMLTN.EQ.'BARE') GOTO 100
        WRITE(6,22) 'Coulomb direct (closed)   (J)',EGDR
        WRITE(7,22) 'Coulomb direct (closed)   (J)',EGDR
        WRITE(6,22) 'Coulomb exchange (closed) (K)',EGXC
        WRITE(7,22) 'Coulomb exchange (closed) (K)',EGXC
        IF(NOPN.EQ.0) GOTO 110
        WRITE(6,22) 'Coulomb direct (open)     (Q)',EQDR
        WRITE(7,22) 'Coulomb direct (open)     (Q)',EQDR
        WRITE(6,22) 'Coulomb exchange (open)   (S)',EQXC
        WRITE(7,22) 'Coulomb exchange (open)   (S)',EQXC
110     CONTINUE
        IF(HMLTN.EQ.'NORL'.OR.HMLTN.EQ.'DHFR') GOTO 100
        WRITE(6,22) 'Breit direct (closed)     (B)',EBDR
        WRITE(7,22) 'Breit direct (closed)     (B)',EBDR
        WRITE(6,22) 'Breit exchange (closed)   (W)',EBXC
        WRITE(7,22) 'Breit exchange (closed)   (W)',EBXC
        IF(NOPN.EQ.0) GOTO 100
        WRITE(6,22) 'Breit direct (open)       (Y)',EMDR
        WRITE(7,22) 'Breit direct (open)       (Y)',EMDR
        WRITE(6,22) 'Breit exchange (open)     (Z)',EMXC
        WRITE(7,22) 'Breit exchange (open)     (Z)',EMXC
100     CONTINUE
        WRITE(6, *) REPEAT('-',72)
        WRITE(7, *) REPEAT('-',72)
        WRITE(6,22) 'Molecule total               ',ETOT
        WRITE(7,22) 'Molecule total               ',ETOT
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
C
C       MATRIX CONSTRUCTION STOPWATCH
23      FORMAT(1X,A,34X,A)
24      FORMAT(1X,A,1X,I3,24X,A)
        WRITE(6, *) REPEAT(' ',72)
        WRITE(7, *) REPEAT(' ',72)
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6, *) REPEAT(' ',22),'Matrix construction stopwatch'
        WRITE(7, *) REPEAT(' ',22),'Matrix construction stopwatch'
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6,23) 'One-electron          ',HMS(T1EL-TMIT)
        WRITE(7,23) 'One-electron          ',HMS(T1EL-TMIT)
        IF(HMLTN.EQ.'BARE') GOTO 207
        WRITE(6,23) 'Coulomb SCF           ',HMS(T2CL-T1EL)
        WRITE(7,23) 'Coulomb SCF           ',HMS(T2CL-T1EL)
        IF(HMLTN.EQ.'NORL'.OR.HMLTN.EQ.'DHFR') GOTO 207
        WRITE(6,23) 'Breit SCF             ',HMS(T2BT-T2CL)
        WRITE(7,23) 'Breit SCF             ',HMS(T2BT-T2CL)
207     CONTINUE
        WRITE(6,23) 'Matrix diagonalisation',HMS(TDGN-T2BT)
        WRITE(7,23) 'Matrix diagonalisation',HMS(TDGN-T2BT)
        WRITE(6, *) REPEAT('-',72)
        WRITE(7, *) REPEAT('-',72)
        WRITE(6,23) 'Total iteration time  ',HMS(TDGN-TMIT)
        WRITE(7,23) 'Total iteration time  ',HMS(TDGN-TMIT)
        WRITE(6,24) 'Time at end of iteration',ITER,STAMP
        WRITE(7,24) 'Time at end of iteration',ITER,STAMP
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
C
C       CONVERGENCE STATUS
25      FORMAT(1X,A,42X,I3)
26      FORMAT(1X,A,37X,F8.5)
27      FORMAT(1X,A,29X,1P,D16.9)
        WRITE(6, *) REPEAT(' ',72)
        WRITE(7, *) REPEAT(' ',72)
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6, *) REPEAT(' ',26),'Convergence analysis'
        WRITE(7, *) REPEAT(' ',26),'Convergence analysis'
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6,25) 'Iteration number           ',ITER
        WRITE(7,25) 'Iteration number           ',ITER
        WRITE(6,25) 'Integral inclusion level   ',ILEV
        WRITE(7,25) 'Integral inclusion level   ',ILEV
        WRITE(6,26) 'Level shift parameter      ',SHLV
        WRITE(7,26) 'Level shift parameter      ',SHLV
        WRITE(6,27) 'Density difference norm    ',DNRM(ITER)
        WRITE(7,27) 'Density difference norm    ',DNRM(ITER)
        WRITE(6,27) 'Weighted energy difference ',WEDN(ITER)
        WRITE(7,27) 'Weighted energy difference ',WEDN(ITER)
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6, *) ' '
        WRITE(7, *) ' '
C
C       TEST FOR CONVERGENCE OR INTEGRAL CLASS UPDATE
30      FORMAT(29X,'Stage 1: (LL|LL)')
31      FORMAT(23X,'Stage 2: (LL|SS) and (SS|LL)')
32      FORMAT(14X,'Stage 2: (LL|SS), (SS|LL), (LS|LS) and (SL|SL)')
33      FORMAT(29X,'Stage 3: (SS|SS)')
34      FORMAT(23X,'Stage 4: (LS|LS) and (SL|SL)')
C
C       IF WEIGHTED ENERGY DIFFERENCE IS NOT CONVERGING, EXIT
        IF(ITER.EQ.10.AND.WEDN(ITER).GT.1.0D-3) THEN
          WRITE(6, *) ' '
          WRITE(7, *) ' '
          WRITE(6, *) REPEAT('-',72)
          WRITE(7, *) REPEAT('-',72)
          WRITE(6, *) 'In HFSCF: energy not converging. Exit BERTHA.'
          WRITE(7, *) 'In HFSCF: energy not converging. Exit BERTHA.'
          WRITE(6, *) REPEAT('-',72)
          WRITE(7, *) REPEAT('-',72)
          STOP
        ENDIF
C
C       BARE NUCLEUS APPROXIMATION: NO COULOMB INTEGRALS AND NO SCF.
        IF(HMLTN.EQ.'BARE'.OR.NOELEC.LE.1) GOTO 300
C
C       CURRENTLY AT STAGE 1: (LL|LL)
        IF(ILEV.EQ.1) THEN
C
C         NON-RELATIVISTIC HAMILTONIAN ONLY INVOLVES (LL|LL)
          IF(HMLTN.EQ.'NORL') THEN
C
C           SATISFIES ALL CRITERIA - SUCCESSFUL CONVERGENCE
            IF(WEDN(ITER).LT.EEPS.AND.DNRM(ITER).LT.DEPS) THEN
              GOTO 300
            ENDIF
C
C         ALL OTHER HAMILTONIANS REQUIRE FURTHER INTEGRAL INCLUSIONS
          ELSE
C
C           IF STAGE 1 HAS NOT CONVERGED, ITERATE AGAIN
            IF(WEDN(ITER).GE.TRSH2) THEN
              SHLV = SHLEV(ILEV)
C           IF STAGE 1 HAS CONVERGED, PROCEED TO STAGE 2
            ELSEIF(WEDN(ITER).LT.TRSH2) THEN
              ILEV = 2
              SHLV = SHLEV(ILEV)
              WRITE(6, *) ' '
              WRITE(7, *) ' '
              WRITE(6, *) REPEAT('*',72)
              WRITE(7, *) REPEAT('*',72)
              IF(HMLTN.NE.'DHFB'.AND.HMLTN.NE.'DHFQ') THEN
                WRITE(6,31) 
                WRITE(7,31)
              ELSE
                WRITE(6,32) 
                WRITE(7,32)        
              ENDIF              
              WRITE(6, *) REPEAT('*',72)
              WRITE(7, *) REPEAT('*',72)
              WRITE(6, *) ' '
              WRITE(7, *) ' '
            ENDIF
C
          ENDIF
C
C       CURRENTLY AT STAGE 2: (LL|SS), (SS|LL), (LS|LS), (SL|SL)
        ELSEIF(ILEV.EQ.2) THEN

C         IF STAGE 2 HAS CONVERGED, PROCEED TO STAGE 3
          IF(WEDN(ITER).LT.TRSH3) THEN
            ILEV = 3
            SHLV = SHLEV(ILEV)
            WRITE(6, *) ' '
            WRITE(7, *) ' '
            WRITE(6, *) REPEAT('*',72)
            WRITE(7, *) REPEAT('*',72)
            WRITE(6,33)
            WRITE(7,33)
            WRITE(6, *) REPEAT('*',72)
            WRITE(7, *) REPEAT('*',72)
            WRITE(6, *) ' '
            WRITE(7, *) ' '
          ENDIF
C
C       CURRENTLY AT STAGE 3: (SS|SS)
        ELSEIF(ILEV.EQ.3) THEN
C
C         SATISFIES ALL CRITERIA - SUCCESSFUL CONVERGENCE
          IF(WEDN(ITER).LT.EEPS.AND.DNRM(ITER).LT.DEPS) THEN
            GOTO 300
          ENDIF
C
        ENDIF
C
C     END LOOP OVER ITERATIONS
      ENDDO
C
C     FORCED EXIT: UNSUCCESSFUL CONVERGENCE
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6, *) 'In HFSCF: convergence not attained. ITER = ',ITER
      WRITE(7, *) 'In HFSCF: convergence not attained. ITER = ',ITER
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      STOP
C
C     EARLY EXIT: SUCCESSFUL CONVERGENCE
300   CONTINUE
C
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('*',72)
      WRITE(7, *) REPEAT('*',72)
      WRITE(6, *) REPEAT(':',72)
      WRITE(7, *) REPEAT(':',72)
      WRITE(6, *) REPEAT(' ',25),'Successful convergence!'
      WRITE(7, *) REPEAT(' ',25),'Successful convergence!'
      WRITE(6, *) REPEAT(':',72)
      WRITE(7, *) REPEAT(':',72)
      WRITE(6, *) REPEAT('*',72)
      WRITE(7, *) REPEAT('*',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C
C**********************************************************************C
C     END OF SELF-CONSISTENT FIELD CALCULATIONS                        C
C**********************************************************************C
C
C     CALCULATE THE PERTUBATIVE VALUE OF THE BREIT ENERGY
      IF(HMLTN.EQ.'DHFP') THEN
        CALL CPU_TIME(TPRTI)
C
C       PRINT CALL TO NEW STAGE OF R-INTEGRALS
        WRITE(6, *) ' '
        WRITE(7, *) ' '
        WRITE(6, *) REPEAT('*',72)
        WRITE(7, *) REPEAT('*',72)
        WRITE(6,34) 
        WRITE(7,34)
        WRITE(6, *) REPEAT('*',72)
        WRITE(7, *) REPEAT('*',72)
        WRITE(6, *) ' '
        WRITE(7, *) ' '
C
C       TITLE FOR CALL TO BREIT ROUTINE
        WRITE(6, *) ' '
        WRITE(7, *) ' '
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6, *) REPEAT(' ',25),'Call to BREIT routine'
        WRITE(7, *) REPEAT(' ',25),'Call to BREIT routine'
        WRITE(6, *) REPEAT('-',72)
        WRITE(7, *) REPEAT('-',72)
C
C       GENERATE MATRIX REP OF BREIT INTERACTION
        CALL BREIT
C
C       ADD BREIT MATRIX TO MOST RECENT FOCK MATRIX
        DO I=1,NDIM
          DO J=1,NDIM
            FOCK(I,J) = FOCK(I,J) + BDIR(I,J) - BXCH(I,J)
          ENDDO
        ENDDO
C
C       DIAGONALISE FOCK MATRIX (REQUIRES LAPACK LIBRARY)
        CALL ZHEGV(1,'V','L',NDIM,FOCK,MDM,OVAP,MDM,
     &                                     EIGEN,WORK,LWK,RWORK,INFO)
        IF(INFO.NE.0) THEN
          WRITE(6, *) 'In HFSCF: eigenvalue solver ZHEGV failed.',INFO
          WRITE(7, *) 'In HFSCF: eigenvalue solver ZHEGV failed.',INFO
        ENDIF
C
C       TRANSFER EIGENVECTORS TO THE C ARRAY
        DO J=1,NDIM
          DO I=1,NDIM
            C(I,J) = FOCK(I,J)
          ENDDO
        ENDDO
C
C       DEDUCT LEVEL SHIFT VALUE FROM VIRTUAL ORBITAL EIGENVALUES
        IF(ITER.EQ.1) THEN
          DO IVIR=NSHIFT+NOCC+1,NDIM
            EIGEN(IVIR) = EIGEN(IVIR) - SHLV
          ENDDO
        ENDIF
C
C       UPDATE EIGENVALUES AND COEFFICIENTS
        OPEN(UNIT=8,FILE=TRIM(WFNFL)//'(+B).wfn',STATUS='UNKNOWN')
        REWIND(UNIT=8)
        DO I=1,NDIM
          WRITE(8, *) EIGEN(I),(C(J,I),J=1,NDIM)
        ENDDO
        CLOSE(UNIT=8)
C
C       RECALCULATE TOTAL ENERGY
        CALL ENERGIES
C
C       UPDATE TOTAL BREIT CALCULATION TIME
        CALL CPU_TIME(TPRTF)
        TBMX = TBMX + TPRTF - TPRTI
C
C       SUMMARISE DIRECT AND EXCHANGE ENERGIES
        WRITE(6,22) 'Breit direct (closed)     (B)',EBDR
        WRITE(7,22) 'Breit direct (closed)     (B)',EBDR
        WRITE(6,22) 'Breit exchange (closed)   (W)',EBXC
        WRITE(7,22) 'Breit exchange (closed)   (W)',EBXC
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6,23) 'Total BREIT time      ',HMS(TBMX)
        WRITE(7,23) 'Total BREIT time      ',HMS(TBMX)
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
C
C       LEVEL ITERATIONS AND TIME
        NMLEV(4) = 1
        TMLEV(4) = TBMX
        ITER     = ITER+1
C
      ENDIF
C
C     TIME AT END OF MOLECULAR CALCULATION
      CALL CPU_TIME(TDM2)
      TSCF = TDM2 - TDM1
C
C     DATE AND TIME AT END OF CALCULATION
      CALL TIMENOW(STAMP)
C
C     PRINT OUT FINAL SCF RESULTS
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',25),'Molecular SCF summary'
      WRITE(7, *) REPEAT(' ',25),'Molecular SCF summary'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     MOLECULAR ENERGIES
      EDIR = EGDR+EQDR+EBDR+EMDR
      EXCH = EGXC+EQXC+EBXC+EMXC
40    FORMAT(1X,A,22X,A,11X,A,14X,A)
41    FORMAT(1X,A,39X,F17.9)
42    FORMAT(1X,A,1X,F17.9,2X,F17.9,2X,F17.9)
      WRITE(6, *) REPEAT(' ',19),'Molecular energies (Hartree units)'
      WRITE(7, *) REPEAT(' ',19),'Molecular energies (Hartree units)'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,40) 'Source','Direct','Exchange','Total'
      WRITE(7,40) 'Source','Direct','Exchange','Total'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,41) 'Nucleus-nucleus ',ENUC
      WRITE(7,41) 'Nucleus-nucleus ',ENUC
      WRITE(6,41) 'Electron-nucleus',EHNC
      WRITE(7,41) 'Electron-nucleus',EHNC
      WRITE(6,41) 'Electron kinetic',EHKN
      WRITE(7,41) 'Electron kinetic',EHKN
      IF(HMLTN.NE.'DHFQ') GOTO 405
      WRITE(6,41) 'Uehling effects ',EUEH
      WRITE(6,41) 'Uehling effects ',EUEH
405   CONTINUE
      IF(HMLTN.EQ.'BARE') GOTO 400
      WRITE(6,42) 'Coulomb (closed)',EGDR,EGXC,ECLG
      WRITE(7,42) 'Coulomb (closed)',EGDR,EGXC,ECLG
      IF(NOPN.EQ.0) GOTO 410
      WRITE(6,42) 'Coulomb (open)  ',EQDR,EQXC,ECLQ
      WRITE(7,42) 'Coulomb (open)  ',EQDR,EQXC,ECLQ
410   CONTINUE
      IF(HMLTN.EQ.'NORL'.OR.HMLTN.EQ.'DHFR') GOTO 400
      WRITE(6,42) 'Breit (closed)  ',EBDR,EBXC,EBRG
      WRITE(7,42) 'Breit (closed)  ',EBDR,EBXC,EBRG
      IF(NOPN.EQ.0) GOTO 400
      WRITE(6,42) 'Breit (open)    ',EMDR,EMXC,EBRQ
      WRITE(7,42) 'Breit (open)    ',EMDR,EMXC,EBRQ
400   CONTINUE
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,42) 'Molecule        ',EDIR,EXCH,ETOT
      WRITE(7,42) 'Molecule        ',EDIR,EXCH,ETOT
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     E-COEFFICIENT AND R-INTEGRAL ANALYSIS
43    FORMAT(1X,A,5X,A,14X,A,20X,A)
44    FORMAT(1X,A,3X,A,14X,A,A)
45    FORMAT(38X,A,A)
46    FORMAT(1X,A,2X,A,33X,A)
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',21),'E-coefficients and R-integrals'
      WRITE(7, *) REPEAT(' ',21),'E-coefficients and R-integrals'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,43) 'E-coefficients','Time','R-integrals','Time'
      WRITE(7,43) 'E-coefficients','Time','R-integrals','Time'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,44) 'E0LL',HMS(TELL),'(LL|LL)            ',HMS(TRLL)
      WRITE(7,44) 'E0LL',HMS(TELL),'(LL|LL)            ',HMS(TRLL)
      IF(HMLTN.EQ.'NORL') GOTO 420
      WRITE(6,44) 'E0SS',HMS(TESS),'(LL|SS) and (SS|LL)',HMS(TRLS)
      WRITE(7,44) 'E0SS',HMS(TESS),'(LL|SS) and (SS|LL)',HMS(TRLS)
      WRITE(6,45)                  '(SS|SS)            ',HMS(TRSS)
      WRITE(7,45)                  '(SS|SS)            ',HMS(TRSS)
      IF(HMLTN.EQ.'BARE'.OR.HMLTN.EQ.'DHFR') GOTO 420
      WRITE(6,44) 'EILS',HMS(TELS),'(LS|LS) and (SL|SL)',HMS(TRBR)
      WRITE(7,44) 'EILS',HMS(TELS),'(LS|LS) and (SL|SL)',HMS(TRBR)
420   CONTINUE
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,46) 'Total',HMS(TELL+TESS+TELS),HMS(TRLL+TRLS+TRSS+TRBR)
      WRITE(7,46) 'Total',HMS(TELL+TESS+TELS),HMS(TRLL+TRLS+TRSS+TRBR)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     MATRIX CONSTRUCTION STOPWATCH
47    FORMAT(1X,A,34X,A)
48    FORMAT(1X,A,1X,25X,A)
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',22),'Matrix construction stopwatch'
      WRITE(7, *) REPEAT(' ',22),'Matrix construction stopwatch'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,47) 'One-electron          ',HMS(THMX)
      WRITE(7,47) 'One-electron          ',HMS(THMX)
      IF(HMLTN.EQ.'BARE') GOTO 430
      WRITE(6,47) 'Coulomb SCF           ',HMS(TGMX)
      WRITE(7,47) 'Coulomb SCF           ',HMS(TGMX)
      IF(HMLTN.EQ.'NORL'.OR.HMLTN.EQ.'DHFR') GOTO 430
      WRITE(6,47) 'Breit SCF             ',HMS(TBMX)
      WRITE(7,47) 'Breit SCF             ',HMS(TBMX)
430   CONTINUE
      WRITE(6,47) 'Matrix diagonalisation',HMS(TEIG)
      WRITE(7,47) 'Matrix diagonalisation',HMS(TEIG)
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,47) 'Total                 ',HMS(TSCF)
      WRITE(7,47) 'Total                 ',HMS(TSCF)
      WRITE(6,48) 'Time at end of calculation',STAMP
      WRITE(7,48) 'Time at end of calculation',STAMP
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     CONVERGENCE ANALYSIS
      T1 = TMLEV(1)
      T2 = TMLEV(2)
      T3 = TMLEV(3)
      T4 = TMLEV(4)
49    FORMAT(1X,A,3X,A,16X,A,6X,A,14X,A)
50    FORMAT(1X,I1,7X,A,2X,F8.5,13X,I3,2X,A)
51    FORMAT(1X,I1,7X,A)
52    FORMAT(9X,A,2X,F8.5,13X,I3,2X,A)
53    FORMAT(1X,A,46X,I3,2X,A)
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',25),'Convergence analysis'
      WRITE(7, *) REPEAT(' ',25),'Convergence analysis'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,49) 'Stage','Inclusion','Shift','Iterations','Time'
      WRITE(7,49) 'Stage','Inclusion','Shift','Iterations','Time'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,50) 1,'(LL|LL)             ',SHLEV(1),NMLEV(1),HMS(T1)
      WRITE(7,50) 1,'(LL|LL)             ',SHLEV(1),NMLEV(1),HMS(T1)
      IF(HMLTN.EQ.'NORL') GOTO 440
      IF(HMLTN.EQ.'DHFB'.OR.HMLTN.EQ.'DHFQ') THEN
        WRITE(6,51) 2,'(LL|SS) and (SS|LL),'
        WRITE(7,51) 2,'(LL|SS) and (SS|LL),'
        WRITE(6,52)   '(LS|LS) and (SL|SL) ',SHLEV(2),NMLEV(2),HMS(T2)
        WRITE(7,52)   '(LS|LS) and (SL|SL) ',SHLEV(2),NMLEV(2),HMS(T2)
      ELSE
        WRITE(6,50) 2,'(LL|SS) and (SS|LL) ',SHLEV(2),NMLEV(2),HMS(T2)
        WRITE(7,50) 2,'(LL|SS) and (SS|LL) ',SHLEV(2),NMLEV(2),HMS(T2)
      ENDIF
      WRITE(6,50) 3,'(SS|SS)             ',SHLEV(3),NMLEV(3),HMS(T3)
      WRITE(7,50) 3,'(SS|SS)             ',SHLEV(3),NMLEV(3),HMS(T3)
      IF(HMLTN.EQ.'DHFP') THEN
        WRITE(6,50) 4,'(LS|LS) and (SL|SL) ',SHLEV(3),NMLEV(4),HMS(T4)
        WRITE(7,50) 4,'(LS|LS) and (SL|SL) ',SHLEV(3),NMLEV(4),HMS(T4)
      ENDIF
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,53) 'Total',ITER,HMS(T1+T2+T3+T4)
      WRITE(7,53) 'Total',ITER,HMS(T1+T2+T3+T4)
440   CONTINUE
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
      RETURN
      END
C
C
      SUBROUTINE OVRLP   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C              OOOOOO  VV    VV RRRRRRR  LL      PPPPPPP               C
C             OO    OO VV    VV RR    RR LL      PP    PP              C
C             OO    OO VV    VV RR    RR LL      PP    PP              C
C             OO    OO VV    VV RR    RR LL      PP    PP              C
C             OO    OO  VV  VV  RRRRRRR  LL      PPPPPPP               C
C             OO    OO   VVVV   RR    RR LL      PP                    C
C              OOOOOO     VV    RR    RR LLLLLLL PP                    C
C                                                                      C
C -------------------------------------------------------------------- C
C  OVRLP CONSTRUCTS A MOLECULAR BASIS FUNCTION OVERLAP MATRIX.         C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 CONE
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 SLL(MBS,MBS,4),SSS(MBS,MBS,4)
C
      DIMENSION RC(MB2,MRC),EXPT(MBS,4),XYZ(3,4),APH(MB2),CP(MB2,3),
     &          PNC(MB2),KQN(4),LQN(4),MQN(4),NFUNS(4)
      DIMENSION INDEX(MCT,-MKP:MKP,2*(MKP+1)*MKP)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
      DATA PI/3.1415926535897932D0/
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     INITIALISE THE OVERLAP ARRAY
      DO I=1,NDIM
        DO J=1,NDIM
          OVAP(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTERS A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     LOOP OVER CENTER A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 0
        ENDIF
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNA    = NFUNCT(LQN(1)+1,ICNTA)
        NFUNS(1) = NFUNA
C
        DO IBAS=1,NFUNA
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNB    = NFUNCT(LQN(2)+1,ICNTB)
        NFUNS(2) = NFUNB
C
        DO IBAS=1,NFUNB
          EXPT(IBAS,2) = EXPSET(IBAS,LQN(2)+1,ICNTB)
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
C     IMPLEMENT ANY AVAILABLE SELECTION RULES HERE                     C
C**********************************************************************C
C
C     SELECTION RULES TO BE MADE BASED ON GEOMETRIC SYMMETRY,
C     ATOMIC COORDINATES AND QUANTUM NUMBER PAIRS. THE IDEA IS TO
C     SKIP CALCULATIONS THAT SHOULD PRODUCE NO OR NEGLIGIBLE EFFECT.
      IF(NCNT.LE.2) THEN
        IF(MQN(1).NE.MQN(2)) GOTO 2000
      ENDIF
C
C**********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS           C
C     OF (MA,MB). FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE          C
C     ORDERED                                                          C
C     11: = (-|MQN(A)|,-|MQN(B)|)                                      C
C     12: = (-|MQN(A)|,+|MQN(B)|)                                      C
C     21: = (+|MQN(A)|,-|MQN(B)|)                                      C
C     22: = (+|MQN(A)|,+|MQN(B)|)                                      C
C**********************************************************************C
C
C     CALCULATE COMPONENT OFFSETS
      IL1 = LARGE(ICNTA,KA,MJA  )
      IL2 = LARGE(ICNTA,KA,MJA+1)
      JL1 = LARGE(ICNTB,KB,MJB  )
      JL2 = LARGE(ICNTB,KB,MJB+1)
C
      IS1 = IL1 + NSHIFT
      IS2 = IL2 + NSHIFT
      JS1 = JL1 + NSHIFT
      JS2 = JL2 + NSHIFT
C
C     CONSTRUCTION OF ONE-ELECTRON OVERLAP MATRIX BY TT' BLOCKS...
C
C     THE PHASE GENERATES E022 AND E012 COEFFS FROM E011 AND E021
      FASE = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
     &                    *DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NFUNA*NFUNB
C
C**********************************************************************C
C     PART 1: THE LL MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)
      NTUVLL = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ELL0 COEFFICIENTS
      IALT = 1
      CALL CPU_TIME(TDM1)
      IF(IEQS.EQ.0) THEN
        CALL EMAKELL(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,1,2,0)
      ELSEIF(IEQS.EQ.1) THEN
        DO ITUV=1,NTUVLL
          IAD = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXAB
          DO M=1,MAXAB
            E11(M,ITUV) = DCMPLX(E0LLFL(IAD+M,1),E0LLFL(IAD+M,2))
            E21(M,ITUV) = DCMPLX(E0LLFL(IAD+M,3),E0LLFL(IAD+M,4))
          ENDDO
        ENDDO
      ENDIF
      CALL CPU_TIME(TDM2)
      TELL = TELL + TDM2 - TDM1
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          EIJ   = EXPT(IBAS,1)+EXPT(JBAS,2)
          EROOT = DSQRT(PI/EIJ)**3
          SLL(IBAS,JBAS,1) = EROOT*E11(M,1)
          SLL(IBAS,JBAS,3) = EROOT*E21(M,1)
          SLL(IBAS,JBAS,2) =-FASE*DCONJG(SLL(IBAS,JBAS,3))
          SLL(IBAS,JBAS,4) = FASE*DCONJG(SLL(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     NON-RELATIVISTIC OVERLAP CALCULATIONS COMPLETE
      IF(HMLTN.EQ.'NORL') GOTO 500
C
C**********************************************************************C
C     PART 2: THE SS MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)+2
      NTUVSS = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ESS0 COEFFICIENTS
      IALT = 1
      CALL CPU_TIME(TDM1)
      IF(IEQS.EQ.0) THEN
        CALL EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,1,2,0)
      ELSEIF(IEQS.EQ.1) THEN
        DO ITUV=1,NTUVSS
          IAD = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXAB
          DO M=1,MAXAB
            E11(M,ITUV) = DCMPLX(E0SSFL(IAD+M,1),E0SSFL(IAD+M,2))
            E21(M,ITUV) = DCMPLX(E0SSFL(IAD+M,3),E0SSFL(IAD+M,4))
          ENDDO
        ENDDO
      ENDIF
      CALL CPU_TIME(TDM2)
      TESS = TESS + TDM2 - TDM1
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          EIJ   = EXPT(IBAS,1)+EXPT(JBAS,2)
          EROOT = DSQRT(PI/EIJ)**3
          SSS(IBAS,JBAS,1) = EROOT*E11(M,1)
          SSS(IBAS,JBAS,3) = EROOT*E21(M,1)
          SSS(IBAS,JBAS,2) =-FASE*DCONJG(SSS(IBAS,JBAS,3))
          SSS(IBAS,JBAS,4) = FASE*DCONJG(SSS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
500   CONTINUE
C
C**********************************************************************C
C     WE NOW HAVE ALL PIECES OF THE OVERLAP MATRIX FOR THIS BLOCK OF   C
C     BASIS FUNCTIONS -- NOW OVERLAY THE RESULTS INTO OVAP.            C
C**********************************************************************C
C
C     CALCULATE COMPONENT OFFSETS
      IL1 = LARGE(ICNTA,KA,MJA  )
      IL2 = LARGE(ICNTA,KA,MJA+1)
      JL1 = LARGE(ICNTB,KB,MJB  )
      JL2 = LARGE(ICNTB,KB,MJB+1)
C
      IS1 = IL1 + NSHIFT
      IS2 = IL2 + NSHIFT
      JS1 = JL1 + NSHIFT
      JS2 = JL2 + NSHIFT
C
C     LL BLOCKS
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            OVAP(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            OVAP(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            OVAP(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            OVAP(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
            OVAP(JL1+JBAS,IL1+IBAS) = DCONJG(OVAP(IL1+IBAS,JL1+JBAS))
            OVAP(JL2+JBAS,IL1+IBAS) = DCONJG(OVAP(IL1+IBAS,JL2+JBAS))
            OVAP(JL1+JBAS,IL2+IBAS) = DCONJG(OVAP(IL2+IBAS,JL1+JBAS))
            OVAP(JL2+JBAS,IL2+IBAS) = DCONJG(OVAP(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            OVAP(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            OVAP(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            OVAP(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            OVAP(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
            OVAP(JL1+JBAS,IL1+IBAS) = DCONJG(OVAP(IL1+IBAS,JL1+JBAS))
            OVAP(JL2+JBAS,IL1+IBAS) = DCONJG(OVAP(IL1+IBAS,JL2+JBAS))
            OVAP(JL1+JBAS,IL2+IBAS) = DCONJG(OVAP(IL2+IBAS,JL1+JBAS))
            OVAP(JL2+JBAS,IL2+IBAS) = DCONJG(OVAP(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     NON-RELATIVISTIC OVERLAP MATRIX COMPLETE
      IF(HMLTN.EQ.'NORL') GOTO 600
C
C     SS BLOCKS
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            OVAP(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            OVAP(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            OVAP(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            OVAP(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
            OVAP(JS1+JBAS,IS1+IBAS) = DCONJG(OVAP(IS1+IBAS,JS1+JBAS))
            OVAP(JS2+JBAS,IS1+IBAS) = DCONJG(OVAP(IS1+IBAS,JS2+JBAS))
            OVAP(JS1+JBAS,IS2+IBAS) = DCONJG(OVAP(IS2+IBAS,JS1+JBAS))
            OVAP(JS2+JBAS,IS2+IBAS) = DCONJG(OVAP(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            OVAP(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            OVAP(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            OVAP(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            OVAP(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
            OVAP(JS1+JBAS,IS1+IBAS) = DCONJG(OVAP(IS1+IBAS,JS1+JBAS))
            OVAP(JS2+JBAS,IS1+IBAS) = DCONJG(OVAP(IS1+IBAS,JS2+JBAS))
            OVAP(JS1+JBAS,IS2+IBAS) = DCONJG(OVAP(IS2+IBAS,JS1+JBAS))
            OVAP(JS2+JBAS,IS2+IBAS) = DCONJG(OVAP(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
600   CONTINUE
C
2000  CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE ONEEL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C               OOOOOO  NN    NN EEEEEEE EEEEEEE LL                    C
C              OO    OO NNN   NN EE      EE      LL                    C
C              OO    OO NNNN  NN EE      EE      LL                    C
C              OO    OO NN NN NN EEEEEE  EEEEEE  LL                    C
C              OO    OO NN  NNNN EE      EE      LL                    C
C              OO    OO NN   NNN EE      EE      LL                    C
C               OOOOOO  NN    NN EEEEEEE EEEEEEE LLLLLLL               C
C                                                                      C
C -------------------------------------------------------------------- C
C  ONEEL CONSTRUCTS A FULL SET OF MULTI-CENTER OVERLAP, KINETIC AND    C
C  NUCLEAR ATTRACTION BASIS FUNCTION MATRIX ELEMENTS.                  C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
C
      COMPLEX*16 SLL(MBS,MBS,4),SSS(MBS,MBS,4),
     &           VLL(MBS,MBS,4),VSS(MBS,MBS,4),
     &           TLL(MBS,MBS,4),TLS(MBS,MBS,4),TSL(MBS,MBS,4)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      COMPLEX*16 E11A,E11B,E11C,TRM11,E21A,E21B,E21C,TRM21
      COMPLEX*16 CTMP1,CTMP2,CTMP3,CTMP4
C
      DIMENSION RC(MB2,MRC),CP(MB2,3),XYZ(3,4),APH(MB2),PNC(MB2),
     &          EXPT(MBS,4),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/DENS/DENC,DENO,DENT
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
      DATA PI/3.1415926535897932D0/
C
C     INITIALISE STORAGE MATRICES
      DO I=1,NDIM
        DO J=1,NDIM
          HNUC(I,J) = DCMPLX(0.0D0,0.0D0)
          HKIN(I,J) = DCMPLX(0.0D0,0.0D0)
          OVAP(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTERS A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     LOOP OVER CENTER A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 2000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 0
        ENDIF
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
         LQN(1) = KQN(1)
        ELSE
         LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNA    = NFUNCT(LQN(1)+1,ICNTA)
        NFUNS(1) = NFUNA
C
        DO IBAS=1,NFUNA
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNB    = NFUNCT(LQN(2)+1,ICNTB)
        NFUNS(2) = NFUNB
C
        DO IBAS=1,NFUNB
          EXPT(IBAS,2) = EXPSET(IBAS,LQN(2)+1,ICNTB)
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
C     IMPLEMENT ANY AVAILABLE SELECTION RULES HERE                     C
C**********************************************************************C
C
C     SELECTION RULES TO BE MADE BASED ON GEOMETRIC SYMMETRY,
C     ATOMIC COORDINATES AND QUANTUM NUMBER PAIRS. THE IDEA IS TO
C     SKIP CALCULATIONS THAT SHOULD PRODUCE NO OR NEGLIGIBLE EFFECT.
      IF(NCNT.LE.2) THEN
        IF(MQN(1).NE.MQN(2)) GOTO 2000
      ENDIF
C
C**********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS           C
C     OF (MA,MB). FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE          C
C     ORDERED                                                          C
C     11: = (-|MQN(A)|,-|MQN(B)|) -> 1                                 C
C     12: = (-|MQN(A)|,+|MQN(B)|) -> 2                                 C
C     21: = (+|MQN(A)|,-|MQN(B)|) -> 3                                 C
C     22: = (+|MQN(A)|,+|MQN(B)|) -> 4                                 C
C**********************************************************************C
C
C     CONSTRUCTION OF ONE-ELECTRON MATRICES BY TT' BLOCKS...
C
C     THE PHASE GENERATES E022 AND E012 COEFFS FROM E011 AND E021
      FASE = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
     &                    *DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXM = NFUNA*NFUNB
C
C     INITIALISE STORAGE ARRAYS
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          DO IB=1,4
            SLL(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
            SSS(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
            TLL(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
            TLS(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
            TSL(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
            VLL(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
            VSS(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     PART 1: THE LL MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)
      NTUVLL = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ELL0 COEFFICIENTS
      IALT = 1
      CALL CPU_TIME(TDM1)
      IF(IEQS.EQ.0) THEN
        CALL EMAKELL(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,1,2,0)
      ELSEIF(IEQS.EQ.1) THEN
        DO ITUV=1,NTUVLL
          IAD = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXM
          DO M=1,MAXM
            E11(M,ITUV) = DCMPLX(E0LLFL(IAD+M,1),E0LLFL(IAD+M,2))
            E21(M,ITUV) = DCMPLX(E0LLFL(IAD+M,3),E0LLFL(IAD+M,4))
          ENDDO
        ENDDO
      ENDIF
      CALL CPU_TIME(TDM2)
      TELL = TELL + TDM2 - TDM1
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          EIJ    = EXPT(IBAS,1)+EXPT(JBAS,2)
          EROOT  = DSQRT(PI/EIJ)**3
          SLL(IBAS,JBAS,1) = EROOT*E11(M,1)
          SLL(IBAS,JBAS,3) = EROOT*E21(M,1)
          SLL(IBAS,JBAS,2) =-FASE*DCONJG(SLL(IBAS,JBAS,3))
          SLL(IBAS,JBAS,4) = FASE*DCONJG(SLL(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     NUCLEAR ATTRACTION MATRIX ELEMENTS
      DO IZ=1,NCNT
C
C       NUCLEAR COORDINATES
        CX = COORD(1,IZ)
        CY = COORD(2,IZ)
        CZ = COORD(3,IZ)
C
C       GAUSSIAN PRODUCT THEOREM OVER BASIS FUNCTIONS
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M   = M+1
            EIJ = EXPT(IBAS,1)+EXPT(JBAS,2)
            ESM = CNUC(IZ)+EIJ
            PX  = (XYZ(1,1)*EXPT(IBAS,1) + XYZ(1,2)*EXPT(JBAS,2))/EIJ
            PY  = (XYZ(2,1)*EXPT(IBAS,1) + XYZ(2,2)*EXPT(JBAS,2))/EIJ
            PZ  = (XYZ(3,1)*EXPT(IBAS,1) + XYZ(3,2)*EXPT(JBAS,2))/EIJ
            APH(M) = (EIJ*CNUC(IZ))/ESM
            PNC(M) = 2.0D0*PI*DSQRT(CNUC(IZ)/ESM)*ZNUC(IZ)/EIJ
            CP(M,1) = CX - PX
            CP(M,2) = CY - PY
            CP(M,3) = CZ - PZ
          ENDDO
        ENDDO
C
C       GENERATE A BATCH OF R-INTEGRALS
        CALL CPU_TIME(TDM1)
        CALL RMAKE(RC,CP,APH,MAXM,LAM)
        CALL CPU_TIME(TDM2)
        TRLL = TRLL + TDM2 - TDM1
C
C       NUCLEAR ATTRACTION INTEGRALS AS A FINITE SUME OF ELL0 AND RC
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M = M+1
            DO ITUV=1,NTUVLL
              VLL(IBAS,JBAS,1) = VLL(IBAS,JBAS,1)
     &                         - PNC(M)*E11(M,ITUV)*RC(M,ITUV)
              VLL(IBAS,JBAS,3) = VLL(IBAS,JBAS,3) 
     &                         - PNC(M)*E21(M,ITUV)*RC(M,ITUV)
            ENDDO
            VLL(IBAS,JBAS,2) =-FASE*DCONJG(VLL(IBAS,JBAS,3))
            VLL(IBAS,JBAS,4) = FASE*DCONJG(VLL(IBAS,JBAS,1))
          ENDDO
        ENDDO
      ENDDO
C
C     CONSTRUCT NON-RELATIVISTIC KINETIC ENERGY INTEGRALS (IOS 91)
      IF(HMLTN.EQ.'NORL') THEN
        RL2 = DFLOAT(2*LQN(2)+3)
        M   = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M     = M+1
            EJ    = EXPT(JBAS,2)
            EIJ   = EXPT(IBAS,1) + EXPT(JBAS,2)
            EROOT = DSQRT(PI/EIJ)**3
            PX = (XYZ(1,1)*EXPT(IBAS,1) + XYZ(1,2)*EXPT(JBAS,2))/EIJ
            PY = (XYZ(2,1)*EXPT(IBAS,1) + XYZ(2,2)*EXPT(JBAS,2))/EIJ
            PZ = (XYZ(3,1)*EXPT(IBAS,1) + XYZ(3,2)*EXPT(JBAS,2))/EIJ
C
            PBX = PX - XYZ(1,2)
            PBY = PY - XYZ(2,2)
            PBZ = PZ - XYZ(3,2)
            PB2 = PBX*PBX + PBY*PBY + PBZ*PBZ
C       
            E0FC = EJ*RL2 - 2.0D0*EJ*EJ*PB2 - 3.00D0*EJ*EJ/EIJ
            E1FC = 4.0D0*EJ*EJ
C
C           TRUNCATE EXPRESSION DEPENDING ON LAM VALUE
C           ALL COMBINATIONS ALLOW FOR THE LAM = 0 MANIFOLD
            TRM11 = E0FC*E11(M,1)
            TRM21 = E0FC*E21(M,1)
C           IF LAM > 0 PROVIDE SECOND BUNCH OF TERMS
            IF(LAM.GE.1) THEN
              E11A = E11(M,INABCD(1,0,0))
              E21A = E21(M,INABCD(1,0,0))
              E11B = E11(M,INABCD(0,1,0))
              E21B = E21(M,INABCD(0,1,0))
              E11C = E11(M,INABCD(0,0,1))
              E21C = E21(M,INABCD(0,0,1))
              TRM11 = TRM11 - E1FC*(PBX*E11A + PBY*E11B + PBZ*E11C)
              TRM21 = TRM21 - E1FC*(PBX*E21A + PBY*E21B + PBZ*E21C)
            ENDIF
C           IF LAM > 1 PROVIDE FINAL BUNCH OF TERMS
            IF(LAM.GE.2) THEN
              E11A = E11(M,INABCD(2,0,0))
              E21A = E21(M,INABCD(2,0,0))
              E11B = E11(M,INABCD(0,2,0))
              E21B = E21(M,INABCD(0,2,0))
              E11C = E11(M,INABCD(0,0,2))
              E21C = E21(M,INABCD(0,0,2))
              TRM11 = TRM11 - E1FC*(E11A + E11B + E11C)
              TRM21 = TRM21 - E1FC*(E21A + E21B + E21C)
            ENDIF
            TLL(IBAS,JBAS,1) = EROOT*TRM11
            TLL(IBAS,JBAS,3) = EROOT*TRM21
            TLL(IBAS,JBAS,2) =-FASE*DCONJG(TLL(IBAS,JBAS,3))
            TLL(IBAS,JBAS,4) = FASE*DCONJG(TLL(IBAS,JBAS,1))
          ENDDO
        ENDDO
C       NON-RELATIVISTIC HAMILTONIAN MATRICES COMPLETE
        GOTO 500
      ENDIF
C     
C**********************************************************************C
C     PART 2: THE SS MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)+2
      NTUVSS = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ESS0 COEFFICIENTS
      IALT = 1
      CALL CPU_TIME(TDM1)
      IF(IEQS.EQ.0) THEN
        CALL EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,1,2,0)
      ELSEIF(IEQS.EQ.1) THEN
        DO ITUV=1,NTUVSS
          IAD = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXM
          DO M=1,MAXM
            E11(M,ITUV) = DCMPLX(E0SSFL(IAD+M,1),E0SSFL(IAD+M,2))
            E21(M,ITUV) = DCMPLX(E0SSFL(IAD+M,3),E0SSFL(IAD+M,4))
          ENDDO
        ENDDO
      ENDIF
      CALL CPU_TIME(TDM2)
      TESS = TESS + TDM2 - TDM1
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          EIJ     = EXPT(IBAS,1)+EXPT(JBAS,2)
          EROOT   = DSQRT(PI/EIJ)**3
          SSS(IBAS,JBAS,1) = EROOT*E11(M,1)
          SSS(IBAS,JBAS,3) = EROOT*E21(M,1)
          SSS(IBAS,JBAS,2) =-FASE*DCONJG(SSS(IBAS,JBAS,3))
          SSS(IBAS,JBAS,4) = FASE*DCONJG(SSS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     NUCLEAR ATTRACTION MATRIX ELEMENTS
      DO IZ=1,NCNT
C
C       NUCLEAR COORDINATES
        CX = COORD(1,IZ)
        CY = COORD(2,IZ)
        CZ = COORD(3,IZ)
C
C       GAUSSIAN PRODUCT THEOREM OVER BASIS FUNCTIONS
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M   = M+1
            EIJ = EXPT(IBAS,1)+EXPT(JBAS,2)
            ESM = CNUC(IZ) + EIJ
            PX  = (XYZ(1,1)*EXPT(IBAS,1) + XYZ(1,2)*EXPT(JBAS,2))/EIJ
            PY  = (XYZ(2,1)*EXPT(IBAS,1) + XYZ(2,2)*EXPT(JBAS,2))/EIJ
            PZ  = (XYZ(3,1)*EXPT(IBAS,1) + XYZ(3,2)*EXPT(JBAS,2))/EIJ
            APH(M) = (EIJ*CNUC(IZ))/ESM
            PNC(M) = 2.0D0*PI*DSQRT(CNUC(IZ)/ESM)*ZNUC(IZ)/EIJ
            CP(M,1) = CX - PX
            CP(M,2) = CY - PY
            CP(M,3) = CZ - PZ
          ENDDO
        ENDDO      
C
C       GENERATE A BATCH OF R-INTEGRALS
        CALL CPU_TIME(TDM1)
        CALL RMAKE(RC,CP,APH,MAXM,LAM)
        CALL CPU_TIME(TDM2)
        TRSS = TRSS + TDM2 - TDM1
C
C       NUCLEAR ATTRACTION INTEGRALS AS A FINITE SUME OF ESS0 AND RC
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M = M+1
            DO ITUV=1,NTUVSS
              VSS(IBAS,JBAS,1) = VSS(IBAS,JBAS,1) 
     &                           + PNC(M)*E11(M,ITUV)*RC(M,ITUV)
              VSS(IBAS,JBAS,3) = VSS(IBAS,JBAS,3) 
     &                           + PNC(M)*E21(M,ITUV)*RC(M,ITUV)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     SUBTRACT THE SS OVERLAP MATRIX AND FINISH CONSTRUCTION
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          VSS(IBAS,JBAS,1) =-VSS(IBAS,JBAS,1) 
     &                       - 2.0D0*CV*CV*SSS(IBAS,JBAS,1)
          VSS(IBAS,JBAS,3) =-VSS(IBAS,JBAS,3) 
     &                       - 2.0D0*CV*CV*SSS(IBAS,JBAS,3)
          VSS(IBAS,JBAS,2) =-FASE*DCONJG(VSS(IBAS,JBAS,3))
          VSS(IBAS,JBAS,4) = FASE*DCONJG(VSS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C**********************************************************************C
C     PART 3: THE SL MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)+2
      NTUVSS = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     KINETIC MATRIX ELEMENTS
      FACT = CV*DSQRT(DFLOAT(2*LQN(2)+3))
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          EIJ    = EXPT(IBAS,1) + EXPT(JBAS,2)
          EJRT   = FACT*DSQRT(EXPT(JBAS,2))
          EROOT  = DSQRT(PI/EIJ)**3
          TSL(IBAS,JBAS,1) = EJRT*EROOT*E11(M,1)
          TSL(IBAS,JBAS,3) = EJRT*EROOT*E21(M,1)
          TSL(IBAS,JBAS,2) =-FASE*DCONJG(TSL(IBAS,JBAS,3))
          TSL(IBAS,JBAS,4) = FASE*DCONJG(TSL(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C
C**********************************************************************C
C     PART 4: THE LS MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)+2
      NTUVSS = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ESS0 COEFFICIENTS
      IALT = 1
C
      CALL CPU_TIME(TDM1)
C
C     DFNOTE: THE INDEX SWAP '2 1' MAKES FILE IMPORT ANNOYING. I'M LAZY.
C             JUST GENERATE THESE EQ'S AS A BATCH (DOESN'T TAKE LONG).
C
C     IF(IEQS.EQ.0) THEN
        CALL EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,2,1,0)
C     ELSEIF(IEQS.EQ.1) THEN
C       DO ITUV=1,NTUVSS
C         IAD = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXM
C         M = 0
C         DO IBAS=1,NFUNA
C           DO JBAS=1,NFUNB
C             M = M+1
C             N = (JBAS-1)*NFUNB + IBAS
C             PI^{SL}_{IJ} = PI^{LS}_{JI}* (THE NEXT LINES ARE WRONG)
C             E11(M,ITUV) = DCMPLX(E0SSFL(IAD+N,1),-E0SSFL(IAD+N,2))
C             E21(M,ITUV) = DCMPLX(E0SSFL(IAD+N,3),-E0SSFL(IAD+N,4))
C           ENDDO
C         ENDDO
C       ENDDO
C     ENDIF
      CALL CPU_TIME(TDM2)
      TESS = TESS + TDM2 - TDM1
C
C     KINETIC MATRIX ELEMENTS
      FACT = CV*DSQRT(DFLOAT(2*LQN(1)+3))
      M = 0
      DO JBAS=1,NFUNB
        DO IBAS=1,NFUNA
          M = M+1
          EIJ    = EXPT(JBAS,2) + EXPT(IBAS,1)
          EIRT   = FACT*DSQRT(EXPT(IBAS,1))
          EROOT  = DSQRT(PI/EIJ)**3
          TLS(IBAS,JBAS,1) = EIRT*EROOT*E11(M,1)
          TLS(IBAS,JBAS,3) = EIRT*EROOT*E21(M,1)
          TLS(IBAS,JBAS,2) =-FASE*DCONJG(TLS(IBAS,JBAS,3))
          TLS(IBAS,JBAS,4) = FASE*DCONJG(TLS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     GENERATE LS MATRICES FROM THE ABOVE SL MATRICES
      M = 0
      DO JBAS=1,NFUNB
        DO IBAS=1,NFUNA
          M = M + 1
          CTMP1 = DCONJG(TLS(IBAS,JBAS,1))
          CTMP2 = DCONJG(TLS(IBAS,JBAS,2))
          CTMP3 = DCONJG(TLS(IBAS,JBAS,3))
          CTMP4 = DCONJG(TLS(IBAS,JBAS,4))
C
          TLS(IBAS,JBAS,1) = CTMP1
          TLS(IBAS,JBAS,2) = CTMP3
          TLS(IBAS,JBAS,3) = CTMP2
          TLS(IBAS,JBAS,4) = CTMP4
        ENDDO
      ENDDO
C
500   CONTINUE
C
C**********************************************************************C
C     WE NOW HAVE ALL PIECES OF HNUC AND HKIN FOR THIS BLOCK.          C
C**********************************************************************C
C
C     CALCULATE COMPONENT OFFSETS
      IL1 = LARGE(ICNTA,KA,MJA  )
      IL2 = LARGE(ICNTA,KA,MJA+1)
      JL1 = LARGE(ICNTB,KB,MJB  )
      JL2 = LARGE(ICNTB,KB,MJB+1)
C
      IS1 = IL1 + NSHIFT
      IS2 = IL2 + NSHIFT
      JS1 = JL1 + NSHIFT
      JS2 = JL2 + NSHIFT
C
C     LL OVERLAP BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            OVAP(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            OVAP(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            OVAP(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            OVAP(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
C
            OVAP(JL1+JBAS,IL1+IBAS) = DCONJG(OVAP(IL1+IBAS,JL1+JBAS))
            OVAP(JL2+JBAS,IL1+IBAS) = DCONJG(OVAP(IL1+IBAS,JL2+JBAS))
            OVAP(JL1+JBAS,IL2+IBAS) = DCONJG(OVAP(IL2+IBAS,JL1+JBAS))
            OVAP(JL2+JBAS,IL2+IBAS) = DCONJG(OVAP(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            OVAP(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            OVAP(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            OVAP(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            OVAP(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
C
            OVAP(JL1+JBAS,IL1+IBAS) = DCONJG(OVAP(IL1+IBAS,JL1+JBAS))
            OVAP(JL2+JBAS,IL1+IBAS) = DCONJG(OVAP(IL1+IBAS,JL2+JBAS))
            OVAP(JL1+JBAS,IL2+IBAS) = DCONJG(OVAP(IL2+IBAS,JL1+JBAS))
            OVAP(JL2+JBAS,IL2+IBAS) = DCONJG(OVAP(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     LL NUCLEAR POTENTIAL BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            HNUC(IL1+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,1)
            HNUC(IL1+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,2)
            HNUC(IL2+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,3)
            HNUC(IL2+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,4)
C
            HNUC(JL1+JBAS,IL1+IBAS) = DCONJG(HNUC(IL1+IBAS,JL1+JBAS))
            HNUC(JL2+JBAS,IL1+IBAS) = DCONJG(HNUC(IL1+IBAS,JL2+JBAS))
            HNUC(JL1+JBAS,IL2+IBAS) = DCONJG(HNUC(IL2+IBAS,JL1+JBAS))
            HNUC(JL2+JBAS,IL2+IBAS) = DCONJG(HNUC(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            HNUC(IL1+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,1)
            HNUC(IL1+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,2)
            HNUC(IL2+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,3)
            HNUC(IL2+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,4)
C
            HNUC(JL1+JBAS,IL1+IBAS) = DCONJG(HNUC(IL1+IBAS,JL1+JBAS))
            HNUC(JL2+JBAS,IL1+IBAS) = DCONJG(HNUC(IL1+IBAS,JL2+JBAS))
            HNUC(JL1+JBAS,IL2+IBAS) = DCONJG(HNUC(IL2+IBAS,JL1+JBAS))
            HNUC(JL2+JBAS,IL2+IBAS) = DCONJG(HNUC(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     NON-RELATIVISTIC HAMILTONIAN HAS A KINETIC MATRIX IN THE LL BLOCK
      IF(HMLTN.EQ.'NORL') THEN
C
C       LL KINETIC BLOCKS
        IF(IL1.GT.JL1) THEN
          DO JBAS=1,NFUNB
            DO IBAS=1,NFUNA
              HKIN(IL1+IBAS,JL1+JBAS) = TLL(IBAS,JBAS,1)
              HKIN(IL1+IBAS,JL2+JBAS) = TLL(IBAS,JBAS,2)
              HKIN(IL2+IBAS,JL1+JBAS) = TLL(IBAS,JBAS,3)
              HKIN(IL2+IBAS,JL2+JBAS) = TLL(IBAS,JBAS,4)
C
              HKIN(JL1+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JL1+JBAS))
              HKIN(JL2+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JL2+JBAS))
              HKIN(JL1+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JL1+JBAS))
              HKIN(JL2+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JL2+JBAS))
            ENDDO
          ENDDO
        ENDIF
C
        IF(IL1.EQ.JL1) THEN
          DO JBAS=1,NFUNB
            DO IBAS=JBAS,NFUNA
              HKIN(IL1+IBAS,JL1+JBAS) = TLL(IBAS,JBAS,1)
              HKIN(IL1+IBAS,JL2+JBAS) = TLL(IBAS,JBAS,2)
              HKIN(IL2+IBAS,JL1+JBAS) = TLL(IBAS,JBAS,3)
              HKIN(IL2+IBAS,JL2+JBAS) = TLL(IBAS,JBAS,4)
C
              HKIN(JL1+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JL1+JBAS))
              HKIN(JL2+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JL2+JBAS))
              HKIN(JL1+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JL1+JBAS))
              HKIN(JL2+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JL2+JBAS))
            ENDDO
          ENDDO
        ENDIF
C      
C       NON-RELATIVISTIC MATRIX CONSTRUCTION COMPLETE
        GOTO 600
C
      ENDIF
C
C     SS OVERLAP BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            OVAP(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            OVAP(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            OVAP(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            OVAP(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
C
            OVAP(JS1+JBAS,IS1+IBAS) = DCONJG(OVAP(IS1+IBAS,JS1+JBAS))
            OVAP(JS2+JBAS,IS1+IBAS) = DCONJG(OVAP(IS1+IBAS,JS2+JBAS))
            OVAP(JS1+JBAS,IS2+IBAS) = DCONJG(OVAP(IS2+IBAS,JS1+JBAS))
            OVAP(JS2+JBAS,IS2+IBAS) = DCONJG(OVAP(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            OVAP(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            OVAP(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            OVAP(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            OVAP(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
C
            OVAP(JS1+JBAS,IS1+IBAS) = DCONJG(OVAP(IS1+IBAS,JS1+JBAS))
            OVAP(JS2+JBAS,IS1+IBAS) = DCONJG(OVAP(IS1+IBAS,JS2+JBAS))
            OVAP(JS1+JBAS,IS2+IBAS) = DCONJG(OVAP(IS2+IBAS,JS1+JBAS))
            OVAP(JS2+JBAS,IS2+IBAS) = DCONJG(OVAP(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     SS NUCLEAR POTENTIAL BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            HNUC(IS1+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,1)
            HNUC(IS1+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,2)
            HNUC(IS2+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,3)
            HNUC(IS2+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,4)
C
            HNUC(JS1+JBAS,IS1+IBAS) = DCONJG(HNUC(IS1+IBAS,JS1+JBAS))
            HNUC(JS2+JBAS,IS1+IBAS) = DCONJG(HNUC(IS1+IBAS,JS2+JBAS))
            HNUC(JS1+JBAS,IS2+IBAS) = DCONJG(HNUC(IS2+IBAS,JS1+JBAS))
            HNUC(JS2+JBAS,IS2+IBAS) = DCONJG(HNUC(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            HNUC(IS1+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,1)
            HNUC(IS1+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,2)
            HNUC(IS2+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,3)
            HNUC(IS2+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,4)
C
            HNUC(JS1+JBAS,IS1+IBAS) = DCONJG(HNUC(IS1+IBAS,JS1+JBAS))
            HNUC(JS2+JBAS,IS1+IBAS) = DCONJG(HNUC(IS1+IBAS,JS2+JBAS))
            HNUC(JS1+JBAS,IS2+IBAS) = DCONJG(HNUC(IS2+IBAS,JS1+JBAS))
            HNUC(JS2+JBAS,IS2+IBAS) = DCONJG(HNUC(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     LS BLOCKS
      IF(IL1.GE.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            HKIN(IL1+IBAS,JS1+JBAS) = TLS(IBAS,JBAS,1)
            HKIN(IL1+IBAS,JS2+JBAS) = TLS(IBAS,JBAS,2)
            HKIN(IL2+IBAS,JS1+JBAS) = TLS(IBAS,JBAS,3)
            HKIN(IL2+IBAS,JS2+JBAS) = TLS(IBAS,JBAS,4)
C
            HKIN(JS1+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JS1+JBAS))
            HKIN(JS2+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JS2+JBAS))
            HKIN(JS1+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JS1+JBAS))
            HKIN(JS2+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     SL BLOCKS
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            HKIN(IS1+IBAS,JL1+JBAS) = TSL(IBAS,JBAS,1)
            HKIN(IS1+IBAS,JL2+JBAS) = TSL(IBAS,JBAS,2)
            HKIN(IS2+IBAS,JL1+JBAS) = TSL(IBAS,JBAS,3)
            HKIN(IS2+IBAS,JL2+JBAS) = TSL(IBAS,JBAS,4)
C
            HKIN(JL1+JBAS,IS1+IBAS) = DCONJG(HKIN(IS1+IBAS,JL1+JBAS))
            HKIN(JL2+JBAS,IS1+IBAS) = DCONJG(HKIN(IS1+IBAS,JL2+JBAS))
            HKIN(JL1+JBAS,IS2+IBAS) = DCONJG(HKIN(IS2+IBAS,JL1+JBAS))
            HKIN(JL2+JBAS,IS2+IBAS) = DCONJG(HKIN(IS2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
600   CONTINUE      
C
C     END LOOPS OVER BASIS PAIRS A,B
2000  CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE UEHLING
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       UU    UU EEEEEEEE HH    HH LL      IIII NN    NN  GGGGGG       C
C       UU    UU EE       HH    HH LL       II  NNN   NN GG    GG      C
C       UU    UU EE       HH    HH LL       II  NNNN  NN GG            C
C       UU    UU EEEEEE   HHHHHHHH LL       II  NN NN NN GG            C
C       UU    UU EE       HH    HH LL       II  NN  NNNN GG   GGG      C
C       UU    UU EE       HH    HH LL       II  NN   NNN GG    GG      C
C        UUUUUU  EEEEEEEE HH    HH LLLLLLL IIII NN    NN  GGGGGG       C
C                                                                      C
C -------------------------------------------------------------------- C
C  UEHLING CONSTRUCTS A FULL SET OF MULTI-CENTER UEHLING INTERACTION   C
C  MATRIX ELEMENTS.                                                    C
C -------------------------------------------------------------------- C
C  DFNOTE: AT THE MOMENT THIS IS JUST A COPY OF THE ONEEL NUCLEAR      C
C          ATTRACTION MATRIX MAKER.                                    C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
C
      COMPLEX*16 VLL(MBS,MBS,4),VSS(MBS,MBS,4)
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 E11A,E11B,E11C,TRM11,E21A,E21B,E21C,TRM21
      COMPLEX*16 CTMP1,CTMP2,CTMP3,CTMP4

      DIMENSION RC(MB2,MRC),CP(MB2,3),XYZ(3,4),APH(MB2),PNC(MB2),
     &          EXPT(MBS,4),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/DENS/DENC,DENO,DENT
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
      DATA PI/3.1415926535897932D0/
C
C     INITIALISE STORAGE MATRIX
      DO I=1,NDIM
        DO J=1,NDIM
          VUEH(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTERS A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     LOOP OVER CENTER A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 2000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 0
        ENDIF
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
         LQN(1) = KQN(1)
        ELSE
         LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNA    = NFUNCT(LQN(1)+1,ICNTA)
        NFUNS(1) = NFUNA
C
        DO IBAS=1,NFUNA
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNB    = NFUNCT(LQN(2)+1,ICNTB)
        NFUNS(2) = NFUNB
C
        DO IBAS=1,NFUNB
          EXPT(IBAS,2) = EXPSET(IBAS,LQN(2)+1,ICNTB)
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
C     IMPLEMENT ANY AVAILABLE SELECTION RULES HERE                     C
C**********************************************************************C
C
C     SELECTION RULES TO BE MADE BASED ON GEOMETRIC SYMMETRY,
C     ATOMIC COORDINATES AND QUANTUM NUMBER PAIRS. THE IDEA IS TO
C     SKIP CALCULATIONS THAT SHOULD PRODUCE NO OR NEGLIGIBLE EFFECT.
      IF(NCNT.LE.2) THEN
        IF(MQN(1).NE.MQN(2)) GOTO 2000
      ENDIF
C
C**********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS           C
C     OF (MA,MB). FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE          C
C     ORDERED                                                          C
C     11: = (-|MQN(A)|,-|MQN(B)|) -> 1                                 C
C     12: = (-|MQN(A)|,+|MQN(B)|) -> 2                                 C
C     21: = (+|MQN(A)|,-|MQN(B)|) -> 3                                 C
C     22: = (+|MQN(A)|,+|MQN(B)|) -> 4                                 C
C**********************************************************************C
C
C     CONSTRUCTION OF ONE-ELECTRON MATRICES BY TT' BLOCKS...
C
C     THE PHASE GENERATES E022 AND E012 COEFFS FROM E011 AND E021
      FASE = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
     &                    *DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NFUNA*NFUNB
C
C     INITIALISE STORAGE ARRAYS
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          DO IB=1,4
            VLL(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
            VSS(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     PART 1: THE LL MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)
      NTUVLL = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ELL0 COEFFICIENTS
      IALT = 1
      CALL CPU_TIME(TDM1)
      IF(IEQS.EQ.0) THEN
        CALL EMAKELL(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,1,2,0)
      ELSEIF(IEQS.EQ.1) THEN
        DO ITUV=1,NTUVLL
          IAD = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXAB
          DO M=1,MAXAB
            E11(M,ITUV) = DCMPLX(E0LLFL(IAD+M,1),E0LLFL(IAD+M,2))
            E21(M,ITUV) = DCMPLX(E0LLFL(IAD+M,3),E0LLFL(IAD+M,4))
          ENDDO
        ENDDO
      ENDIF
      CALL CPU_TIME(TDM2)
      TELL = TELL + TDM2 - TDM1
C
C     NUCLEAR ATTRACTION MATRIX ELEMENTS
      DO IZ=1,NCNT
C
C       NUCLEAR COORDINATES
        CX = COORD(1,IZ)
        CY = COORD(2,IZ)
        CZ = COORD(3,IZ)
C
C       GAUSSIAN PRODUCT THEOREM OVER BASIS FUNCTIONS
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M   = M+1
            EIJ = EXPT(IBAS,1)+EXPT(JBAS,2)
            ESM = CNUC(IZ)+EIJ
            PX  = (XYZ(1,1)*EXPT(IBAS,1) + XYZ(1,2)*EXPT(JBAS,2))/EIJ
            PY  = (XYZ(2,1)*EXPT(IBAS,1) + XYZ(2,2)*EXPT(JBAS,2))/EIJ
            PZ  = (XYZ(3,1)*EXPT(IBAS,1) + XYZ(3,2)*EXPT(JBAS,2))/EIJ
            APH(M) = (EIJ*CNUC(IZ))/ESM
            PNC(M) = 2.0D0*PI*DSQRT(CNUC(IZ)/ESM)*ZNUC(IZ)/EIJ
            CP(M,1) = CX - PX
            CP(M,2) = CY - PY
            CP(M,3) = CZ - PZ
          ENDDO
        ENDDO
C
C       GENERATE A BATCH OF R-INTEGRALS
        CALL CPU_TIME(TDM1)
        CALL RMAKE(RC,CP,APH,MAXAB,LAM)
        CALL CPU_TIME(TDM2)
        TRLL = TRLL + TDM2 - TDM1
C
C       NUCLEAR ATTRACTION INTEGRALS AS A FINITE SUME OF ELL0 AND RC
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M = M+1
            DO ITUV=1,NTUVLL
              VLL(IBAS,JBAS,1) = VLL(IBAS,JBAS,1)
     &                         - PNC(M)*E11(M,ITUV)*RC(M,ITUV)
              VLL(IBAS,JBAS,3) = VLL(IBAS,JBAS,3) 
     &                         - PNC(M)*E21(M,ITUV)*RC(M,ITUV)
            ENDDO
            VLL(IBAS,JBAS,2) =-FASE*DCONJG(VLL(IBAS,JBAS,3))
            VLL(IBAS,JBAS,4) = FASE*DCONJG(VLL(IBAS,JBAS,1))
          ENDDO
        ENDDO
      ENDDO
C
C     NON-RELATIVISTIC CALCULATIONS ARE COMPLETE
      IF(HMLTN.EQ.'NORL') GOTO 500
C     
C**********************************************************************C
C     PART 2: THE SS MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)+2
      NTUVSS = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ESS0 COEFFICIENTS
      IALT = 1
      CALL CPU_TIME(TDM1)
      IF(IEQS.EQ.0) THEN
        CALL EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,1,2,0)
      ELSEIF(IEQS.EQ.1) THEN
        DO ITUV=1,NTUVSS
          IAD = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXAB
          DO M=1,MAXAB
            E11(M,ITUV) = DCMPLX(E0SSFL(IAD+M,1),E0SSFL(IAD+M,2))
            E21(M,ITUV) = DCMPLX(E0SSFL(IAD+M,3),E0SSFL(IAD+M,4))
          ENDDO
        ENDDO
      ENDIF
      CALL CPU_TIME(TDM2)
      TESS = TESS + TDM2 - TDM1
C
C     NUCLEAR ATTRACTION MATRIX ELEMENTS
      DO IZ=1,NCNT
C
C       NUCLEAR COORDINATES
        CX = COORD(1,IZ)
        CY = COORD(2,IZ)
        CZ = COORD(3,IZ)
C
C       GAUSSIAN PRODUCT THEOREM OVER BASIS FUNCTIONS
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M   = M+1
            EIJ = EXPT(IBAS,1)+EXPT(JBAS,2)
            ESM = CNUC(IZ) + EIJ
            PX  = (XYZ(1,1)*EXPT(IBAS,1) + XYZ(1,2)*EXPT(JBAS,2))/EIJ
            PY  = (XYZ(2,1)*EXPT(IBAS,1) + XYZ(2,2)*EXPT(JBAS,2))/EIJ
            PZ  = (XYZ(3,1)*EXPT(IBAS,1) + XYZ(3,2)*EXPT(JBAS,2))/EIJ
            APH(M) = (EIJ*CNUC(IZ))/ESM
            PNC(M) = 2.0D0*PI*DSQRT(CNUC(IZ)/ESM)*ZNUC(IZ)/EIJ
            CP(M,1) = CX - PX
            CP(M,2) = CY - PY
            CP(M,3) = CZ - PZ
          ENDDO
        ENDDO      
C
C       GENERATE A BATCH OF R-INTEGRALS
        CALL CPU_TIME(TDM1)
        CALL RMAKE(RC,CP,APH,MAXAB,LAM)
        CALL CPU_TIME(TDM2)
        TRSS = TRSS + TDM2 - TDM1
C
C       NUCLEAR ATTRACTION INTEGRALS AS A FINITE SUME OF ELL0 AND RC
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M = M+1
            DO ITUV=1,NTUVSS
              VSS(IBAS,JBAS,1) = VSS(IBAS,JBAS,1) 
     &                           + PNC(M)*E11(M,ITUV)*RC(M,ITUV)
              VSS(IBAS,JBAS,3) = VSS(IBAS,JBAS,3) 
     &                           + PNC(M)*E21(M,ITUV)*RC(M,ITUV)
            ENDDO
            VSS(IBAS,JBAS,2) =-FASE*DCONJG(VSS(IBAS,JBAS,3))
            VSS(IBAS,JBAS,4) = FASE*DCONJG(VSS(IBAS,JBAS,1))
          ENDDO
        ENDDO
      ENDDO
C
500   CONTINUE
C
C**********************************************************************C
C     WE NOW HAVE ALL PIECES OF VUEH FOR THIS BLOCK.                   C
C**********************************************************************C
C
C     CALCULATE COMPONENT OFFSETS
      IL1 = LARGE(ICNTA,KA,MJA  )
      IL2 = LARGE(ICNTA,KA,MJA+1)
      JL1 = LARGE(ICNTB,KB,MJB  )
      JL2 = LARGE(ICNTB,KB,MJB+1)
C
      IS1 = IL1 + NSHIFT
      IS2 = IL2 + NSHIFT
      JS1 = JL1 + NSHIFT
      JS2 = JL2 + NSHIFT
C
C     LL UEHLING INTERACTION BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            VUEH(IL1+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,1)
            VUEH(IL1+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,2)
            VUEH(IL2+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,3)
            VUEH(IL2+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,4)
C
            VUEH(JL1+JBAS,IL1+IBAS) = DCONJG(VUEH(IL1+IBAS,JL1+JBAS))
            VUEH(JL2+JBAS,IL1+IBAS) = DCONJG(VUEH(IL1+IBAS,JL2+JBAS))
            VUEH(JL1+JBAS,IL2+IBAS) = DCONJG(VUEH(IL2+IBAS,JL1+JBAS))
            VUEH(JL2+JBAS,IL2+IBAS) = DCONJG(VUEH(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            VUEH(IL1+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,1)
            VUEH(IL1+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,2)
            VUEH(IL2+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,3)
            VUEH(IL2+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,4)
C
            VUEH(JL1+JBAS,IL1+IBAS) = DCONJG(VUEH(IL1+IBAS,JL1+JBAS))
            VUEH(JL2+JBAS,IL1+IBAS) = DCONJG(VUEH(IL1+IBAS,JL2+JBAS))
            VUEH(JL1+JBAS,IL2+IBAS) = DCONJG(VUEH(IL2+IBAS,JL1+JBAS))
            VUEH(JL2+JBAS,IL2+IBAS) = DCONJG(VUEH(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C      
C     NON-RELATIVISTIC MATRIX CONSTRUCTION COMPLETE
      GOTO 600
C
C     SS UEHLING INTERACTION BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            VUEH(IS1+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,1)
            VUEH(IS1+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,2)
            VUEH(IS2+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,3)
            VUEH(IS2+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,4)
C
            VUEH(JS1+JBAS,IS1+IBAS) = DCONJG(VUEH(IS1+IBAS,JS1+JBAS))
            VUEH(JS2+JBAS,IS1+IBAS) = DCONJG(VUEH(IS1+IBAS,JS2+JBAS))
            VUEH(JS1+JBAS,IS2+IBAS) = DCONJG(VUEH(IS2+IBAS,JS1+JBAS))
            VUEH(JS2+JBAS,IS2+IBAS) = DCONJG(VUEH(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            VUEH(IS1+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,1)
            VUEH(IS1+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,2)
            VUEH(IS2+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,3)
            VUEH(IS2+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,4)
C
            VUEH(JS1+JBAS,IS1+IBAS) = DCONJG(VUEH(IS1+IBAS,JS1+JBAS))
            VUEH(JS2+JBAS,IS1+IBAS) = DCONJG(VUEH(IS1+IBAS,JS2+JBAS))
            VUEH(JS1+JBAS,IS2+IBAS) = DCONJG(VUEH(IS2+IBAS,JS1+JBAS))
            VUEH(JS2+JBAS,IS2+IBAS) = DCONJG(VUEH(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
600   CONTINUE      
C
C     END LOOPS OVER BASIS PAIRS A,B
2000  CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE COULOMB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    CCCCCC   OOOOOO  UU    UU LL      OOOOOO  MM       MM BBBBBBB     C
C   CC    CC OO    OO UU    UU LL     OO    OO MMM     MMM BB    BB    C
C   CC       OO    OO UU    UU LL     OO    OO MMMM   MMMM BB    BB    C
C   CC       OO    OO UU    UU LL     OO    OO MM MM MM MM BBBBBBB     C
C   CC       OO    OO UU    UU LL     OO    OO MM  MMM  MM BB    BB    C
C   CC    CC OO    OO UU    UU LL     OO    OO MM   M   MM BB    BB    C
C    CCCCCC   OOOOOO   UUUUUU  LLLLLLL OOOOOO  MM       MM BBBBBBB     C
C                                                                      C
C -------------------------------------------------------------------- C
C  COULOMB GENERATES ELECTRON REPULSION INTEGRALS IN BATCHES AND       C
C  CALCULATES THE SCF COULOMB MATRIX (G), APPLYING IT DIRECTLY TO THE  C
C  FOCK MATRIX. IN THE CASE OF OPEN SUBSHELLS, THE Q MATRIX IS ALSO    C
C  INCLUDED. INTEGRAL SYMMETRY PARTIALLY EXPLOITED, BUT WITH ROOM FOR  C
C  IMPROVEMENT (GEOMETRIC SYMM, R-INT SYMM, E-COEFF SYMM).             C
C -------------------------------------------------------------------- C
C  DFNOTE: THIS ROUTINE COULD BENEFIT FROM PARALLELISATION -- OPENMP.  C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 T1(MDM),T2(MDM),T3(MDM)
      COMPLEX*16 RR(MB2,16)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NFUNS(4),LQN(4)
      DIMENSION ITQN(2),IFLG(11),ISCF(11,6)
      DIMENSION INDEX(MCT,-MKP:MKP,2*(MKP+1)*MKP)
C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/QNMS/ICNBAS(MDM),KQNBAS(MDM),MQNBAS(MDM)
      COMMON/SHLL/ACFF,BCFF,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN,NOELEC
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     ISCF TELLS WHICH INTEGRALS TO INCLUDE BASED ON OVERLAP COMBINATION
      DATA ISCF/1,0,1,0,1,0,1,0,1,1,0,
     &          1,0,0,0,1,0,0,0,1,1,0,
     &          0,1,1,0,0,0,0,0,1,1,0,
     &          1,0,0,1,1,0,0,0,1,0,0,
     &          0,1,0,1,0,0,0,0,1,0,0,
     &          0,1,0,0,0,0,0,0,1,0,0/
C
C     INITIALISE STORAGE MATRICES
      DO I=1,NDIM
        DO J=1,NDIM
          GDIR(I,J) = DCMPLX(0.0D0,0.0D0)
          GXCH(I,J) = DCMPLX(0.0D0,0.0D0)
          QDIR(I,J) = DCMPLX(0.0D0,0.0D0)
          QXCH(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     BARE NUCLEUS APPROXIMATION: NO MEAN-FIELD COULOMB MATRIX
      IF(HMLTN.EQ.'BARE') RETURN
C
C**********************************************************************C
C     INDEXING ROUTINE: SO WE CAN SET BASIS FUNCTION LABELS BY BLOCK   C
C**********************************************************************C
      ICOUNT = 0
C
C     LOOP OVER NUCLEAR CENTERS
      DO ICNT=1,NCNT
C
C       LOOP OVER KAPPA VALUES FOR THIS NUCLEAR CENTER
        DO KN=1,NKAP(ICNT)
C
C         IMPORT KAPPA, MAXIMUM MQN
          KAPPA = KVALS(KN,ICNT)
          MJMAX  = 2*IABS(KAPPA)-1
C
C         LOOP OVER MQN VALUES AND RECORD INDEX
          DO MJ=1,MJMAX,2
            ICOUNT = ICOUNT+1
            INDEX(ICNT,KAPPA,MJ) = ICOUNT
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     SET THE LIMITS FOR RUNNING OVER DENSITY COMBINATIONS             C
C**********************************************************************C
C
      IF(HMLTN.EQ.'NORL') THEN
        ITSTRT = 1
        ITSTOP = 1
        ITSKIP = 1
      ELSE
        ITSTRT = 1
        ITSTOP = 4
        ITSKIP = 3
      ENDIF
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTERS A AND B (USE INDEX 1000)      C
C**********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR A AND B: T = {L} OR {L,S}
      DO 1000 IT1=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITQN(1) = IT1
C
C       CALCULATE STARTING ADDRESS
        IF(IT1.EQ.1) THEN
          NADDAB = 0
        ELSE
          NADDAB = NSHIFT
        ENDIF
C
C     LOOP OVER CENTER A
      DO 1000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 1000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 0
        ENDIF
C
C     LOOP OVER KQN(A) VALUES
      DO 1000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNS(1) = NFUNCT(LQN(1)+1,ICNTA)
        DO IBAS=1,NFUNS(1)
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 1000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
         LQN(2) = KQN(2)
        ELSE
         LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNS(2) = NFUNCT(LQN(2)+1,ICNTB)
        DO JBAS=1, NFUNS(2)
          EXPT(JBAS,2) = EXPSET(JBAS,LQN(2)+1,ICNTB)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 1000 MA=1,IABS(KQN(1))
        MJA    = 2*MA-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 1000 MB=1,IABS(KQN(2))
        MJB    = 2*MB-1
        MQN(2) = MJB
C
C     CALCULATE NEW BLOCK OF E(AB) COEFFS AT NEXT OPPORTUNITY
      IEAB = 1
      IABLL = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB)
      IABSS = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB)
C
C**********************************************************************C
C      SECOND LAYER OF LOOPS, OVER CENTERS C AND D (INDEX 2000/2500)   C
C**********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR C AND D: T' = {L} OR {L,S}
      DO 2000 IT2=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITQN(2) = IT2
C
C       CALCULATE STARTING ADDRESS
        IF(IT2.EQ.1) THEN
          NADDCD = 0
        ELSE
          NADDCD = NSHIFT
        ENDIF

C     LOOP OVER CENTER C
      DO 2000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER C
        XYZ(1,3) = COORD(1,ICNTC)
        XYZ(2,3) = COORD(2,ICNTC)
        XYZ(3,3) = COORD(3,ICNTC)
C
C     LOOP OVER CENTER D
      DO 2000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER D
        XYZ(1,4) = COORD(1,ICNTD)
        XYZ(2,4) = COORD(2,ICNTD)
        XYZ(3,4) = COORD(3,ICNTD)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER C AND D
        IF(ICNTC.EQ.ICNTD) THEN
          INUCCD = 1
        ELSE
          INUCCD = 0
        ENDIF
C
C       PARAMETER FOR ATOMIC OR MULTICNTRE INTEGRAL
        IF(INUCAB*INUCCD.EQ.1.AND.ICNTA.EQ.ICNTC) THEN
          IATOM = 1
        ELSE
          IATOM = 0
        ENDIF
C
C ***   STAGES: DECISION TREE FOR SKIPPING MULTI-CENTER CONTRIBUTIONS
        IF(IATOM.EQ.0) THEN
C
C >>      STAGE 1: INCLUDE ONLY (LL|LL) REPULSION INTEGRALS
          IF(ILEV.EQ.1.AND.IT1+IT2.GT.2) THEN
            GOTO 2200
          ENDIF
C
C >>      STAGE 2: INCLUDE ONLY (LL|SS) AND (SS|LL) REPULSION INTEGRALS
          IF(ILEV.EQ.2.AND.IT1+IT2.GT.5) THEN
            GOTO 2200
          ENDIF
C
C ***   END OF ILEV DECISION TREE
        ENDIF
C
C     LOOP OVER KQN(C) VALUES
      DO 2500 KC=1,NKAP(ICNTC)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR C
        KQN(3) = KVALS(KC,ICNTC)
        IF(KQN(3).GT.0) THEN
          LQN(3) = KQN(3)
        ELSE
          LQN(3) =-KQN(3)-1
        ENDIF
C         
        NFUNS(3) = NFUNCT(LQN(3)+1,ICNTC)
        DO KBAS=1,NFUNS(3)
          EXPT(KBAS,3) = EXPSET(KBAS,LQN(3)+1,ICNTC)
        ENDDO
C
C     LOOP OVER KQN(D) VALUES
      DO 2500 KD=1,NKAP(ICNTD)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR D
        KQN(4) = KVALS(KD,ICNTD)
        IF(KQN(4).GT.0) THEN
          LQN(4) = KQN(4)
        ELSE
          LQN(4) =-KQN(4)-1
        ENDIF
C
        NFUNS(4) = NFUNCT(LQN(4)+1,ICNTD)
        DO LBAS=1,NFUNS(4)
          EXPT(LBAS,4) = EXPSET(LBAS,LQN(4)+1,ICNTD)
        ENDDO
C
C     LOOP OVER |MQN(C)| VALUES
      DO 2500 MC=1,IABS(KQN(3))
        MJC    = 2*MC-1
        MQN(3) = MJC
C
C     LOOP OVER |MQN(D)| VALUES
      DO 2500 MD=1,IABS(KQN(4))
        MJD    = 2*MD-1
        MQN(4) = MJD
C
C     CALCULATE NEW BLOCK OF E(CD) COEFFS AT NEXT OPPORTUNITY
      IECD = 1
      ICDLL = IAD0LL(ICNTC,ICNTD,KC,KD,MC,MD)
      ICDSS = IAD0SS(ICNTC,ICNTD,KC,KD,MC,MD)
C
C**********************************************************************C
C     FOR THIS CHOICE OF A,B,C AND D, COMPUTE ADDRESSES AND PHASES     C
C**********************************************************************C
C
C     CALCULATE BLOCK INDICES FOR {ABCD} COMBINATIONS
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
C     COMBINED BLOCK INDEX IN A TWO-FUNCTION LIST
      IQ12 = (IQ1*(IQ1-1))/2 + IQ2
      IQ34 = (IQ3*(IQ3-1))/2 + IQ4
C
C     FURTHER DEFINE STARTING ADDRESSES FOR {ABCD} BASIS LABELS
      IA1 = LARGE(ICNTA,KA,MJA  ) + NADDAB
      IB1 = LARGE(ICNTB,KB,MJB  ) + NADDAB
      IC1 = LARGE(ICNTC,KC,MJC  ) + NADDCD
      ID1 = LARGE(ICNTD,KD,MJD  ) + NADDCD
C
      IA2 = LARGE(ICNTA,KA,MJA+1) + NADDAB
      IB2 = LARGE(ICNTB,KB,MJB+1) + NADDAB
      IC2 = LARGE(ICNTC,KC,MJC+1) + NADDCD
      ID2 = LARGE(ICNTD,KD,MJD+1) + NADDCD
C
      JA1 = LARGE(ICNTA,KA,MJA  ) + NADDAB
      JB1 = LARGE(ICNTB,KB,MJB  ) + NADDAB
      JC1 = LARGE(ICNTC,KC,MJC  ) + NADDCD
      JD1 = LARGE(ICNTD,KD,MJD  ) + NADDCD
C
      JA2 = LARGE(ICNTA,KA,MJA+1) + NADDAB
      JB2 = LARGE(ICNTB,KB,MJB+1) + NADDAB
      JC2 = LARGE(ICNTC,KC,MJC+1) + NADDCD
      JD2 = LARGE(ICNTD,KD,MJD+1) + NADDCD
C
C     CALCULATE KQN PHASE FACTORS FOR PERMUTING INTEGRALS
      IF((KQN(1)*KQN(2)).GT.0) THEN 
        PKAB = 1.0D0
      ELSE
        PKAB =-1.0D0
      ENDIF
C        
      IF((KQN(3)*KQN(4)).GT.0) THEN 
        PKCD = 1.0D0
      ELSE
        PKCD =-1.0D0
      ENDIF
C
C     CALCULATE MQN PHASE FACTORS FOR PERMUTING INTEGRALS
      PMAB1 = DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PMAB2 = DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      PMCD1 = DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PMCD2 = DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
C     COMBINATIONS OF PHASE FACTORS TO BE USED FOR COULOMB
      PAB1 = PKAB*PMAB1
      PAB2 = PKAB*PMAB2
      PCD1 = PKCD*PMCD1
      PCD2 = PKCD*PMCD2
C
C**********************************************************************C
C     SKIP BATCHES WHICH CONFORM TO INTEGRAL SYMMETRIES                C
C**********************************************************************C
C
C     DIATOMIC MOLECULES CARRY STRICT SELECTION RULES ON MQNS
      IF(NCNT.LE.2) THEN
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) GOTO 2998
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) GOTO 2998
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) GOTO 2998
        GOTO 2999
      ENDIF
2998  CONTINUE

C     DECISION TREE FOR SKIPPING CONTRIBUTIONS DUE TO INTEGRAL SYMMETRY
      IF(IQ1.LT.IQ2) GOTO 2999
      IF(IQ3.LT.IQ4) GOTO 2999
      IF(IQ12.LT.IQ34) GOTO 2999
C
C     INDICATE BLOCKS TO BE INCLUDED AHEAD GIVEN A,B,C,D BASIS QNMS...
C     A =/= B AND C =/= D WITH AB LIST VALUE =/= CD LIST VALUE
      IF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 1
C     A=/=B AND C=/=D WITH IND(AB) = IND(CD)
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.EQ.IQ34) THEN
        ITSCF = 2
C     A=/=B AND C = D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 3
C     A = B AND C=/=D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 4
C     A = B AND C = D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 5
C     A = B AND C = D WITH IND(AB) = IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.EQ.IQ34) THEN
        ITSCF = 6
C     COMBINATION OF A,B,C,D NOT TO BE INCLUDED -- USE MATRIX CONJ LATER
      ELSE
        GO TO 2999
      ENDIF
C
C     READ IN FLAG VALUES FROM ISCF DATA BLOCK
      DO N=1,11
        IFLG(N) = ISCF(N,ITSCF)
      ENDDO
C
C     INCLUDE SPECIAL CASES...
C     A=/=B AND C=/=D WITH IND(AB)=/=IND(CD)
      IF(ITSCF.EQ.1) THEN
C       A = C
        IF(IQ1.EQ.IQ3) IFLG( 6) = 1
C       B = D
        IF(IQ2.EQ.IQ4) IFLG(11) = 1
C       B = C
        IF(IQ2.EQ.IQ3) IFLG( 8) = 1
C     A=/=B AND C = D WITH IND(AB)=/=IND(CD)
      ELSEIF(ITSCF.EQ.3) THEN
C       B = C
        IF(IQ2.EQ.IQ3) IFLG( 8) = 1
C     A = B AND C=/=D WITH IND(AB)=/=IND(CD)
      ELSEIF(ITSCF.EQ.4) THEN
C       B = C
        IF(IQ2.EQ.IQ3) IFLG( 8) = 1
      ENDIF
C
C**********************************************************************C
C     THIRD LAYER OF LOOPS, OVER BASIS FUNCTIONS A AND B (3000)        C
C -------------------------------------------------------------------- C
C     THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH         C
C     GENERATE CLOSED AND OPEN COULOMB MATRICES FROM SPINOR INTEGRALS. C
C     THESE INCLUDE PHASE FACTORS FOR THE PERMUTATION OF               C
C     KQN(1) <-> KQN(2) AND MQN(1) <-> MQN(2)                          C
C -------------------------------------------------------------------- C
C     (RSCF 86, 87)                                                    C
C**********************************************************************C
C
C     LOOP OVER ELEMENTS OF FOCK MATRIX BLOCK
      DO 3000 IBAS=1,NFUNS(1)
      DO 3000 JBAS=1,NFUNS(2)
C
C       GENERATE BATCH OF ELECTRON REPULSION INTEGRALS
        CALL ERI(RR,XYZ,KQN,MQN,EXPT,NFUNS,ITQN,IBAS,JBAS,IEAB,IECD)
C
C       THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH
C       GENERATE THE GDIR AND GXCH MATRIX FROM THE SPINOR INTEGRALS.
C
C**********************************************************************C
C     CONSTRUCT CLOSED-SHELL COULOMB INTERACTION MATRICES, QDIR/QXCH.  C
C**********************************************************************C
C
C       FIRST IFLG BATCH (DIRECT)
        IF(IFLG(1).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GDIR(IA1+IBAS,JB1+JBAS) = GDIR(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 1)*DENT(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M, 2)*DENT(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M, 3)*DENT(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M, 4)*DENT(IC2+KBAS,JD2+LBAS)
     &                     +      PCD1*RR(M, 4)*DENT(JD1+LBAS,IC1+KBAS)
     &                     +      PCD2*RR(M, 2)*DENT(JD1+LBAS,IC2+KBAS)
     &                     +      PCD2*RR(M, 3)*DENT(JD2+LBAS,IC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENT(JD2+LBAS,IC2+KBAS)
C
              GDIR(IA1+IBAS,JB2+JBAS) = GDIR(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M, 5)*DENT(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M, 6)*DENT(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M, 7)*DENT(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M, 8)*DENT(IC2+KBAS,JD2+LBAS)
     &                     +      PCD1*RR(M, 8)*DENT(JD1+LBAS,IC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENT(JD1+LBAS,IC2+KBAS)
     &                     +      PCD2*RR(M, 7)*DENT(JD2+LBAS,IC1+KBAS)
     &                     +      PCD1*RR(M, 5)*DENT(JD2+LBAS,IC2+KBAS)

C
              GDIR(IA2+IBAS,JB1+JBAS) = GDIR(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M, 9)*DENT(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M,10)*DENT(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M,11)*DENT(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M,12)*DENT(IC2+KBAS,JD2+LBAS)
     &                     +      PCD1*RR(M,12)*DENT(JD1+LBAS,IC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENT(JD1+LBAS,IC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENT(JD2+LBAS,IC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENT(JD2+LBAS,IC2+KBAS)
C
              GDIR(IA2+IBAS,JB2+JBAS) = GDIR(IA2+IBAS,JB2+JBAS)
     &                     +           RR(M,13)*DENT(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M,14)*DENT(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M,15)*DENT(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M,16)*DENT(IC2+KBAS,JD2+LBAS)
     &                     +      PCD1*RR(M,16)*DENT(JD1+LBAS,IC1+KBAS)
     &                     +      PCD2*RR(M,14)*DENT(JD1+LBAS,IC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENT(JD2+LBAS,IC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENT(JD2+LBAS,IC2+KBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SECOND IFLG BATCH (DIRECT)
        IF(IFLG(2).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GDIR(IA1+IBAS,JB1+JBAS) = GDIR(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 1)*DENT(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M, 2)*DENT(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M, 3)*DENT(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M, 4)*DENT(IC2+KBAS,JD2+LBAS)
C
              GDIR(IA1+IBAS,JB2+JBAS) = GDIR(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M, 5)*DENT(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M, 6)*DENT(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M, 7)*DENT(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M, 8)*DENT(IC2+KBAS,JD2+LBAS)
C
              GDIR(IA2+IBAS,JB1+JBAS) = GDIR(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M, 9)*DENT(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M,10)*DENT(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M,11)*DENT(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M,12)*DENT(IC2+KBAS,JD2+LBAS)
C
              GDIR(IA2+IBAS,JB2+JBAS) = GDIR(IA2+IBAS,JB2+JBAS)
     &                     +           RR(M,13)*DENT(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M,14)*DENT(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M,15)*DENT(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M,16)*DENT(IC2+KBAS,JD2+LBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       THIRD IFLG BATCH (DIRECT)
        IF(IFLG(3).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GDIR(IC1+KBAS,JD1+LBAS) = GDIR(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M, 1)*DENT(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 5)*DENT(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M, 9)*DENT(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,13)*DENT(IA2+IBAS,JB2+JBAS)
     &                     + PAB1*     RR(M,13)*DENT(JB1+JBAS,IA1+IBAS)
     &                     + PAB2*     RR(M, 5)*DENT(JB1+JBAS,IA2+IBAS)
     &                     + PAB2*     RR(M, 9)*DENT(JB2+JBAS,IA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENT(JB2+JBAS,IA2+IBAS)
C
              GDIR(IC1+KBAS,JD2+LBAS) = GDIR(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M, 2)*DENT(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 6)*DENT(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M,10)*DENT(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,14)*DENT(IA2+IBAS,JB2+JBAS)
     &                     + PAB1*     RR(M,14)*DENT(JB1+JBAS,IA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENT(JB1+JBAS,IA2+IBAS)
     &                     + PAB2*     RR(M,10)*DENT(JB2+JBAS,IA1+IBAS)
     &                     + PAB1*     RR(M, 2)*DENT(JB2+JBAS,IA2+IBAS)
C
              GDIR(IC2+KBAS,JD1+LBAS) = GDIR(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M, 3)*DENT(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 7)*DENT(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M,11)*DENT(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,15)*DENT(IA2+IBAS,JB2+JBAS)
     &                     + PAB1*     RR(M,15)*DENT(JB1+JBAS,IA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENT(JB1+JBAS,IA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENT(JB2+JBAS,IA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENT(JB2+JBAS,IA2+IBAS)
C
              GDIR(IC2+KBAS,JD2+LBAS) = GDIR(IC2+KBAS,JD2+LBAS)
     &                     +           RR(M, 4)*DENT(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 8)*DENT(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M,12)*DENT(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,16)*DENT(IA2+IBAS,JB2+JBAS)
     &                     + PAB1*     RR(M,16)*DENT(JB1+JBAS,IA1+IBAS)
     &                     + PAB2*     RR(M, 8)*DENT(JB1+JBAS,IA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENT(JB2+JBAS,IA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENT(JB2+JBAS,IA2+IBAS)
C 
            ENDDO
          ENDDO
        ENDIF
C
C       FOURTH IFLG BATCH (DIRECT)
        IF(IFLG(4).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GDIR(IC1+KBAS,JD1+LBAS) = GDIR(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M, 1)*DENT(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 5)*DENT(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M, 9)*DENT(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,13)*DENT(IA2+IBAS,JB2+JBAS)
C
              GDIR(IC1+KBAS,JD2+LBAS) = GDIR(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M, 2)*DENT(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 6)*DENT(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M,10)*DENT(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,14)*DENT(IA2+IBAS,JB2+JBAS)
C
              GDIR(IC2+KBAS,JD1+LBAS) = GDIR(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M, 3)*DENT(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 7)*DENT(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M,11)*DENT(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,15)*DENT(IA2+IBAS,JB2+JBAS)
C
              GDIR(IC2+KBAS,JD2+LBAS) = GDIR(IC2+KBAS,JD2+LBAS)
     &                     +           RR(M, 4)*DENT(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 8)*DENT(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M,12)*DENT(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,16)*DENT(IA2+IBAS,JB2+JBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       FIFTH IFLG BATCH (EXCHANGE)
        IF(IFLG(5).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GXCH(IA1+IBAS,JC1+KBAS) = GXCH(IA1+IBAS,JC1+KBAS)
     &                     +      PCD1*RR(M, 4)*DENT(ID1+LBAS,JB1+JBAS)
     &                     +      PCD1*RR(M, 8)*DENT(ID1+LBAS,JB2+JBAS)
     &                     +      PCD2*RR(M, 3)*DENT(ID2+LBAS,JB1+JBAS)
     &                     +      PCD2*RR(M, 7)*DENT(ID2+LBAS,JB2+JBAS)
C
              GXCH(IA1+IBAS,JC2+KBAS) = GXCH(IA1+IBAS,JC2+KBAS)
     &                     +      PCD2*RR(M, 2)*DENT(ID1+LBAS,JB1+JBAS)
     &                     +      PCD2*RR(M, 6)*DENT(ID1+LBAS,JB2+JBAS)
     &                     +      PCD1*RR(M, 1)*DENT(ID2+LBAS,JB1+JBAS)
     &                     +      PCD1*RR(M, 5)*DENT(ID2+LBAS,JB2+JBAS)
C
              GXCH(IA2+IBAS,JC1+KBAS) = GXCH(IA2+IBAS,JC1+KBAS)
     &                     +      PCD1*RR(M,12)*DENT(ID1+LBAS,JB1+JBAS)
     &                     +      PCD1*RR(M,16)*DENT(ID1+LBAS,JB2+JBAS)
     &                     +      PCD2*RR(M,11)*DENT(ID2+LBAS,JB1+JBAS)
     &                     +      PCD2*RR(M,15)*DENT(ID2+LBAS,JB2+JBAS)
C
              GXCH(IA2+IBAS,JC2+KBAS) = GXCH(IA2+IBAS,JC2+KBAS)
     &                     +      PCD2*RR(M,10)*DENT(ID1+LBAS,JB1+JBAS)
     &                     +      PCD2*RR(M,14)*DENT(ID1+LBAS,JB2+JBAS)
     &                     +      PCD1*RR(M, 9)*DENT(ID2+LBAS,JB1+JBAS)
     &                     +      PCD1*RR(M,13)*DENT(ID2+LBAS,JB2+JBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SIXTH IFLG BATCH (EXCHANGE)
        IF(IFLG(6).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GXCH(IC1+KBAS,JA1+IBAS) = GXCH(IC1+KBAS,JA1+IBAS)
     &                     + PAB1*     RR(M,13)*DENT(JB1+JBAS,ID1+LBAS)
     &                     + PAB1*     RR(M,14)*DENT(JB1+JBAS,ID2+LBAS)
     &                     + PAB2*     RR(M, 9)*DENT(JB2+JBAS,ID1+LBAS)
     &                     + PAB2*     RR(M,10)*DENT(JB2+JBAS,ID2+LBAS)
C
              GXCH(IC1+KBAS,JA2+IBAS) = GXCH(IC1+KBAS,JA2+IBAS)
     &                     + PAB2*     RR(M, 5)*DENT(JB1+JBAS,ID1+LBAS)
     &                     + PAB2*     RR(M, 6)*DENT(JB1+JBAS,ID2+LBAS)
     &                     + PAB1*     RR(M, 1)*DENT(JB2+JBAS,ID1+LBAS)
     &                     + PAB1*     RR(M, 2)*DENT(JB2+JBAS,ID2+LBAS)
C
              GXCH(IC2+KBAS,JA1+IBAS) = GXCH(IC2+KBAS,JA1+IBAS)
     &                     + PAB1*     RR(M,15)*DENT(JB1+JBAS,ID1+LBAS)
     &                     + PAB1*     RR(M,16)*DENT(JB1+JBAS,ID2+LBAS)
     &                     + PAB2*     RR(M,11)*DENT(JB2+JBAS,ID1+LBAS)
     &                     + PAB2*     RR(M,12)*DENT(JB2+JBAS,ID2+LBAS)
C
              GXCH(IC2+KBAS,JA2+IBAS) = GXCH(IC2+KBAS,JA2+IBAS)
     &                     + PAB2*     RR(M, 7)*DENT(JB1+JBAS,ID1+LBAS)
     &                     + PAB2*     RR(M, 8)*DENT(JB1+JBAS,ID2+LBAS)
     &                     + PAB1*     RR(M, 3)*DENT(JB2+JBAS,ID1+LBAS)
     &                     + PAB1*     RR(M, 4)*DENT(JB2+JBAS,ID2+LBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SEVENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(7).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GXCH(IB1+JBAS,JC1+KBAS) = GXCH(IB1+JBAS,JC1+KBAS)
     &                     + PAB1*PCD1*RR(M,16)*DENT(ID1+LBAS,JA1+IBAS)
     &                     + PAB2*PCD1*RR(M, 8)*DENT(ID1+LBAS,JA2+IBAS)
     &                     + PAB1*PCD2*RR(M,15)*DENT(ID2+LBAS,JA1+IBAS)
     &                     + PAB2*PCD2*RR(M, 7)*DENT(ID2+LBAS,JA2+IBAS)
C
              GXCH(IB1+JBAS,JC2+KBAS) = GXCH(IB1+JBAS,JC2+KBAS)
     &                     + PAB1*PCD2*RR(M,14)*DENT(ID1+LBAS,JA1+IBAS)
     &                     + PAB2*PCD2*RR(M, 6)*DENT(ID1+LBAS,JA2+IBAS)
     &                     + PAB1*PCD1*RR(M,13)*DENT(ID2+LBAS,JA1+IBAS)
     &                     + PAB2*PCD1*RR(M, 5)*DENT(ID2+LBAS,JA2+IBAS)
C
              GXCH(IB2+JBAS,JC1+KBAS) = GXCH(IB2+JBAS,JC1+KBAS)
     &                     + PAB2*PCD1*RR(M,12)*DENT(ID1+LBAS,JA1+IBAS)
     &                     + PAB1*PCD1*RR(M, 4)*DENT(ID1+LBAS,JA2+IBAS)
     &                     + PAB2*PCD2*RR(M,11)*DENT(ID2+LBAS,JA1+IBAS)
     &                     + PAB1*PCD2*RR(M, 3)*DENT(ID2+LBAS,JA2+IBAS)
C
              GXCH(IB2+JBAS,JC2+KBAS) = GXCH(IB2+JBAS,JC2+KBAS)
     &                     + PAB2*PCD2*RR(M,10)*DENT(ID1+LBAS,JA1+IBAS)
     &                     + PAB1*PCD2*RR(M, 2)*DENT(ID1+LBAS,JA2+IBAS)
     &                     + PAB2*PCD1*RR(M, 9)*DENT(ID2+LBAS,JA1+IBAS)
     &                     + PAB1*PCD1*RR(M, 1)*DENT(ID2+LBAS,JA2+IBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       EIGHTH IFLG BATCH (EXCHANGE)
        IF(IFLG(8).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GXCH(IC1+KBAS,JB1+JBAS) = GXCH(IC1+KBAS,JB1+JBAS)
     &                     +           RR(M, 1)*DENT(JA1+IBAS,ID1+LBAS)
     &                     +           RR(M, 2)*DENT(JA1+IBAS,ID2+LBAS)
     &                     +           RR(M, 9)*DENT(JA2+IBAS,ID1+LBAS)
     &                     +           RR(M,10)*DENT(JA2+IBAS,ID2+LBAS)
C
              GXCH(IC1+KBAS,JB2+JBAS) = GXCH(IC1+KBAS,JB2+JBAS)
     &                     +           RR(M, 5)*DENT(JA1+IBAS,ID1+LBAS)
     &                     +           RR(M, 6)*DENT(JA1+IBAS,ID2+LBAS)
     &                     +           RR(M,13)*DENT(JA2+IBAS,ID1+LBAS)
     &                     +           RR(M,14)*DENT(JA2+IBAS,ID2+LBAS)
C
              GXCH(IC2+KBAS,JB1+JBAS) = GXCH(IC2+KBAS,JB1+JBAS)
     &                     +           RR(M, 3)*DENT(JA1+IBAS,ID1+LBAS)
     &                     +           RR(M, 4)*DENT(JA1+IBAS,ID2+LBAS)
     &                     +           RR(M,11)*DENT(JA2+IBAS,ID1+LBAS)
     &                     +           RR(M,12)*DENT(JA2+IBAS,ID2+LBAS)
C
              GXCH(IC2+KBAS,JB2+JBAS) = GXCH(IC2+KBAS,JB2+JBAS)
     &                     +           RR(M, 7)*DENT(JA1+IBAS,ID1+LBAS)
     &                     +           RR(M, 8)*DENT(JA1+IBAS,ID2+LBAS)
     &                     +           RR(M,15)*DENT(JA2+IBAS,ID1+LBAS)
     &                     +           RR(M,16)*DENT(JA2+IBAS,ID2+LBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       NINTH IFLG BATCH (EXCHANGE)
        IF(IFLG(9).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C  
              GXCH(IA1+IBAS,JD1+LBAS) = GXCH(IA1+IBAS,JD1+LBAS)
     &                     +           RR(M, 1)*DENT(IC1+KBAS,JB1+JBAS)
     &                     +           RR(M, 5)*DENT(IC1+KBAS,JB2+JBAS)
     &                     +           RR(M, 3)*DENT(IC2+KBAS,JB1+JBAS)
     &                     +           RR(M, 7)*DENT(IC2+KBAS,JB2+JBAS)
C
              GXCH(IA1+IBAS,JD2+LBAS) = GXCH(IA1+IBAS,JD2+LBAS)
     &                     +           RR(M, 2)*DENT(IC1+KBAS,JB1+JBAS)
     &                     +           RR(M, 6)*DENT(IC1+KBAS,JB2+JBAS)
     &                     +           RR(M, 4)*DENT(IC2+KBAS,JB1+JBAS)
     &                     +           RR(M, 8)*DENT(IC2+KBAS,JB2+JBAS)
C
              GXCH(IA2+IBAS,JD1+LBAS) = GXCH(IA2+IBAS,JD1+LBAS)
     &                     +           RR(M, 9)*DENT(IC1+KBAS,JB1+JBAS)
     &                     +           RR(M,13)*DENT(IC1+KBAS,JB2+JBAS)
     &                     +           RR(M,11)*DENT(IC2+KBAS,JB1+JBAS)
     &                     +           RR(M,15)*DENT(IC2+KBAS,JB2+JBAS)
C
              GXCH(IA2+IBAS,JD2+LBAS) = GXCH(IA2+IBAS,JD2+LBAS)
     &                     +           RR(M,10)*DENT(IC1+KBAS,JB1+JBAS)
     &                     +           RR(M,14)*DENT(IC1+KBAS,JB2+JBAS)
     &                     +           RR(M,12)*DENT(IC2+KBAS,JB1+JBAS)
     &                     +           RR(M,16)*DENT(IC2+KBAS,JB2+JBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       TENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(10).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C
              GXCH(IB1+JBAS,JD1+LBAS) = GXCH(IB1+JBAS,JD1+LBAS)
     &                     + PAB1*     RR(M,13)*DENT(IC1+KBAS,JA1+IBAS)
     &                     + PAB2*     RR(M, 5)*DENT(IC1+KBAS,JA2+IBAS)
     &                     + PAB1*     RR(M,15)*DENT(IC2+KBAS,JA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENT(IC2+KBAS,JA2+IBAS)
C
              GXCH(IB1+JBAS,JD2+LBAS) = GXCH(IB1+JBAS,JD2+LBAS)
     &                     + PAB1*     RR(M,14)*DENT(IC1+KBAS,JA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENT(IC1+KBAS,JA2+IBAS)
     &                     + PAB1*     RR(M,16)*DENT(IC2+KBAS,JA1+IBAS)
     &                     + PAB2*     RR(M, 8)*DENT(IC2+KBAS,JA2+IBAS)
C
              GXCH(IB2+JBAS,JD1+LBAS) = GXCH(IB2+JBAS,JD1+LBAS)
     &                     + PAB2*     RR(M, 9)*DENT(IC1+KBAS,JA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENT(IC1+KBAS,JA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENT(IC2+KBAS,JA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENT(IC2+KBAS,JA2+IBAS)
C
              GXCH(IB2+JBAS,JD2+LBAS) = GXCH(IB2+JBAS,JD2+LBAS)
     &                     + PAB2*     RR(M,10)*DENT(IC1+KBAS,JA1+IBAS)
     &                     + PAB1*     RR(M, 2)*DENT(IC1+KBAS,JA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENT(IC2+KBAS,JA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENT(IC2+KBAS,JA2+IBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       ELEVENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(11).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C
              GXCH(ID1+LBAS,JB1+JBAS) = GXCH(ID1+LBAS,JB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENT(JA1+IBAS,IC1+KBAS)
     &                     +      PCD2*RR(M, 2)*DENT(JA1+IBAS,IC2+KBAS)
     &                     +      PCD1*RR(M,12)*DENT(JA2+IBAS,IC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENT(JA2+IBAS,IC2+KBAS)
C
              GXCH(ID1+LBAS,JB2+JBAS) = GXCH(ID1+LBAS,JB2+JBAS)
     &                     +      PCD1*RR(M, 8)*DENT(JA1+IBAS,IC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENT(JA1+IBAS,IC2+KBAS)
     &                     +      PCD1*RR(M,16)*DENT(JA2+IBAS,IC1+KBAS)
     &                     +      PCD2*RR(M,14)*DENT(JA2+IBAS,IC2+KBAS)
C
              GXCH(ID2+LBAS,JB1+JBAS) = GXCH(ID2+LBAS,JB1+JBAS)
     &                     +      PCD2*RR(M, 3)*DENT(JA1+IBAS,IC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENT(JA1+IBAS,IC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENT(JA2+IBAS,IC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENT(JA2+IBAS,IC2+KBAS)
C
              GXCH(ID2+LBAS,JB2+JBAS) = GXCH(ID2+LBAS,JB2+JBAS)
     &                     +      PCD2*RR(M, 7)*DENT(JA1+IBAS,IC1+KBAS)
     &                     +      PCD1*RR(M, 5)*DENT(JA1+IBAS,IC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENT(JA2+IBAS,IC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENT(JA2+IBAS,IC2+KBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C**********************************************************************C
C     CONSTRUCT OPEN-SHELL COULOMB INTERACTION MATRICES, QDIR/QXCH.    C
C**********************************************************************C
C
        IF(NOPN.EQ.0) GOTO 399
C
C       FIRST IFLG BATCH (DIRECT)
        IF(IFLG(1).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QDIR(IA1+IBAS,JB1+JBAS) = QDIR(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 1)*DENO(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M, 2)*DENO(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M, 3)*DENO(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M, 4)*DENO(IC2+KBAS,JD2+LBAS)
     &                     +      PCD1*RR(M, 4)*DENO(JD1+LBAS,IC1+KBAS)
     &                     +      PCD2*RR(M, 2)*DENO(JD1+LBAS,IC2+KBAS)
     &                     +      PCD2*RR(M, 3)*DENO(JD2+LBAS,IC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENO(JD2+LBAS,IC2+KBAS)
C
              QDIR(IA1+IBAS,JB2+JBAS) = QDIR(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M, 5)*DENO(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M, 6)*DENO(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M, 7)*DENO(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M, 8)*DENO(IC2+KBAS,JD2+LBAS)
     &                     +      PCD1*RR(M, 8)*DENO(JD1+LBAS,IC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENO(JD1+LBAS,IC2+KBAS)
     &                     +      PCD2*RR(M, 7)*DENO(JD2+LBAS,IC1+KBAS)
     &                     +      PCD1*RR(M, 5)*DENO(JD2+LBAS,IC2+KBAS)
C
              QDIR(IA2+IBAS,JB1+JBAS) = QDIR(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M, 9)*DENO(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M,10)*DENO(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M,11)*DENO(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M,12)*DENO(IC2+KBAS,JD2+LBAS)
     &                     +      PCD1*RR(M,12)*DENO(JD1+LBAS,IC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENO(JD1+LBAS,IC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENO(JD2+LBAS,IC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENO(JD2+LBAS,IC2+KBAS)
C
              QDIR(IA2+IBAS,JB2+JBAS) = QDIR(IA2+IBAS,JB2+JBAS)
     &                     +           RR(M,13)*DENO(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M,14)*DENO(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M,15)*DENO(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M,16)*DENO(IC2+KBAS,JD2+LBAS)
     &                     +      PCD1*RR(M,16)*DENO(JD1+LBAS,IC1+KBAS)
     &                     +      PCD2*RR(M,14)*DENO(JD1+LBAS,IC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENO(JD2+LBAS,IC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENO(JD2+LBAS,IC2+KBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SECOND IFLG BATCH (DIRECT)
        IF(IFLG(2).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QDIR(IA1+IBAS,JB1+JBAS) = QDIR(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 1)*DENO(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M, 2)*DENO(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M, 3)*DENO(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M, 4)*DENO(IC2+KBAS,JD2+LBAS)
C
              QDIR(IA1+IBAS,JB2+JBAS) = QDIR(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M, 5)*DENO(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M, 6)*DENO(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M, 7)*DENO(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M, 8)*DENO(IC2+KBAS,JD2+LBAS)
C
              QDIR(IA2+IBAS,JB1+JBAS) = QDIR(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M, 9)*DENO(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M,10)*DENO(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M,11)*DENO(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M,12)*DENO(IC2+KBAS,JD2+LBAS)
C
              QDIR(IA2+IBAS,JB2+JBAS) = QDIR(IA2+IBAS,JB2+JBAS)
     &                     +           RR(M,13)*DENO(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M,14)*DENO(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M,15)*DENO(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M,16)*DENO(IC2+KBAS,JD2+LBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       THIRD IFLG BATCH (DIRECT)
        IF(IFLG(3).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QDIR(IC1+KBAS,JD1+LBAS) = QDIR(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M, 1)*DENO(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 5)*DENO(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M, 9)*DENO(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,13)*DENO(IA2+IBAS,JB2+JBAS)
     &                     + PAB1*     RR(M,13)*DENO(JB1+JBAS,IA1+IBAS)
     &                     + PAB2*     RR(M, 5)*DENO(JB1+JBAS,IA2+IBAS)
     &                     + PAB2*     RR(M, 9)*DENO(JB2+JBAS,IA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENO(JB2+JBAS,IA2+IBAS)
C
              QDIR(IC1+KBAS,JD2+LBAS) = QDIR(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M, 2)*DENO(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 6)*DENO(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M,10)*DENO(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,14)*DENO(IA2+IBAS,JB2+JBAS)
     &                     + PAB1*     RR(M,14)*DENO(JB1+JBAS,IA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENO(JB1+JBAS,IA2+IBAS)
     &                     + PAB2*     RR(M,10)*DENO(JB2+JBAS,IA1+IBAS)
     &                     + PAB1*     RR(M, 2)*DENO(JB2+JBAS,IA2+IBAS)
C
              QDIR(IC2+KBAS,JD1+LBAS) = QDIR(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M, 3)*DENO(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 7)*DENO(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M,11)*DENO(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,15)*DENO(IA2+IBAS,JB2+JBAS)
     &                     + PAB1*     RR(M,15)*DENO(JB1+JBAS,IA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENO(JB1+JBAS,IA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENO(JB2+JBAS,IA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENO(JB2+JBAS,IA2+IBAS)
C
              QDIR(IC2+KBAS,JD2+LBAS) = QDIR(IC2+KBAS,JD2+LBAS)
     &                     +           RR(M, 4)*DENO(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 8)*DENO(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M,12)*DENO(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,16)*DENO(IA2+IBAS,JB2+JBAS)
     &                     + PAB1*     RR(M,16)*DENO(JB1+JBAS,IA1+IBAS)
     &                     + PAB2*     RR(M, 8)*DENO(JB1+JBAS,IA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENO(JB2+JBAS,IA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENO(JB2+JBAS,IA2+IBAS)
C 
            ENDDO
          ENDDO
        ENDIF
C
C       FOURTH IFLG BATCH (DIRECT)
        IF(IFLG(4).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QDIR(IC1+KBAS,JD1+LBAS) = QDIR(IC1+KBAS,JD1+LBAS)
     &                     +           RR(M, 1)*DENO(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 5)*DENO(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M, 9)*DENO(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,13)*DENO(IA2+IBAS,JB2+JBAS)
C
              QDIR(IC1+KBAS,JD2+LBAS) = QDIR(IC1+KBAS,JD2+LBAS)
     &                     +           RR(M, 2)*DENO(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 6)*DENO(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M,10)*DENO(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,14)*DENO(IA2+IBAS,JB2+JBAS)
C
              QDIR(IC2+KBAS,JD1+LBAS) = QDIR(IC2+KBAS,JD1+LBAS)
     &                     +           RR(M, 3)*DENO(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 7)*DENO(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M,11)*DENO(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,15)*DENO(IA2+IBAS,JB2+JBAS)
C
              QDIR(IC2+KBAS,JD2+LBAS) = QDIR(IC2+KBAS,JD2+LBAS)
     &                     +           RR(M, 4)*DENO(IA1+IBAS,JB1+JBAS)
     &                     +           RR(M, 8)*DENO(IA1+IBAS,JB2+JBAS)
     &                     +           RR(M,12)*DENO(IA2+IBAS,JB1+JBAS)
     &                     +           RR(M,16)*DENO(IA2+IBAS,JB2+JBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       FIFTH IFLG BATCH (EXCHANGE)
        IF(IFLG(5).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QXCH(IA1+IBAS,JC1+KBAS) = QXCH(IA1+IBAS,JC1+KBAS)
     &                     +      PCD1*RR(M, 4)*DENO(ID1+LBAS,JB1+JBAS)
     &                     +      PCD1*RR(M, 8)*DENO(ID1+LBAS,JB2+JBAS)
     &                     +      PCD2*RR(M, 3)*DENO(ID2+LBAS,JB1+JBAS)
     &                     +      PCD2*RR(M, 7)*DENO(ID2+LBAS,JB2+JBAS)
C
              QXCH(IA1+IBAS,JC2+KBAS) = QXCH(IA1+IBAS,JC2+KBAS)
     &                     +      PCD2*RR(M, 2)*DENO(ID1+LBAS,JB1+JBAS)
     &                     +      PCD2*RR(M, 6)*DENO(ID1+LBAS,JB2+JBAS)
     &                     +      PCD1*RR(M, 1)*DENO(ID2+LBAS,JB1+JBAS)
     &                     +      PCD1*RR(M, 5)*DENO(ID2+LBAS,JB2+JBAS)
C
              QXCH(IA2+IBAS,JC1+KBAS) = QXCH(IA2+IBAS,JC1+KBAS)
     &                     +      PCD1*RR(M,12)*DENO(ID1+LBAS,JB1+JBAS)
     &                     +      PCD1*RR(M,16)*DENO(ID1+LBAS,JB2+JBAS)
     &                     +      PCD2*RR(M,11)*DENO(ID2+LBAS,JB1+JBAS)
     &                     +      PCD2*RR(M,15)*DENO(ID2+LBAS,JB2+JBAS)
C
              QXCH(IA2+IBAS,JC2+KBAS) = QXCH(IA2+IBAS,JC2+KBAS)
     &                     +      PCD2*RR(M,10)*DENO(ID1+LBAS,JB1+JBAS)
     &                     +      PCD2*RR(M,14)*DENO(ID1+LBAS,JB2+JBAS)
     &                     +      PCD1*RR(M, 9)*DENO(ID2+LBAS,JB1+JBAS)
     &                     +      PCD1*RR(M,13)*DENO(ID2+LBAS,JB2+JBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SIXTH IFLG BATCH (EXCHANGE)
        IF(IFLG(6).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QXCH(IC1+KBAS,JA1+IBAS) = QXCH(IC1+KBAS,JA1+IBAS)
     &                     + PAB1*     RR(M,13)*DENO(JB1+JBAS,ID1+LBAS)
     &                     + PAB1*     RR(M,14)*DENO(JB1+JBAS,ID2+LBAS)
     &                     + PAB2*     RR(M, 9)*DENO(JB2+JBAS,ID1+LBAS)
     &                     + PAB2*     RR(M,10)*DENO(JB2+JBAS,ID2+LBAS)
C
              QXCH(IC1+KBAS,JA2+IBAS) = QXCH(IC1+KBAS,JA2+IBAS)
     &                     + PAB2*     RR(M, 5)*DENO(JB1+JBAS,ID1+LBAS)
     &                     + PAB2*     RR(M, 6)*DENO(JB1+JBAS,ID2+LBAS)
     &                     + PAB1*     RR(M, 1)*DENO(JB2+JBAS,ID1+LBAS)
     &                     + PAB1*     RR(M, 2)*DENO(JB2+JBAS,ID2+LBAS)
C
              QXCH(IC2+KBAS,JA1+IBAS) = QXCH(IC2+KBAS,JA1+IBAS)
     &                     + PAB1*     RR(M,15)*DENO(JB1+JBAS,ID1+LBAS)
     &                     + PAB1*     RR(M,16)*DENO(JB1+JBAS,ID2+LBAS)
     &                     + PAB2*     RR(M,11)*DENO(JB2+JBAS,ID1+LBAS)
     &                     + PAB2*     RR(M,12)*DENO(JB2+JBAS,ID2+LBAS)
C
              QXCH(IC2+KBAS,JA2+IBAS) = QXCH(IC2+KBAS,JA2+IBAS)
     &                     + PAB2*     RR(M, 7)*DENO(JB1+JBAS,ID1+LBAS)
     &                     + PAB2*     RR(M, 8)*DENO(JB1+JBAS,ID2+LBAS)
     &                     + PAB1*     RR(M, 3)*DENO(JB2+JBAS,ID1+LBAS)
     &                     + PAB1*     RR(M, 4)*DENO(JB2+JBAS,ID2+LBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SEVENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(7).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QXCH(IB1+JBAS,JC1+KBAS) = QXCH(IB1+JBAS,JC1+KBAS)
     &                     + PAB1*PCD1*RR(M,16)*DENO(ID1+LBAS,JA1+IBAS)
     &                     + PAB2*PCD1*RR(M, 8)*DENO(ID1+LBAS,JA2+IBAS)
     &                     + PAB1*PCD2*RR(M,15)*DENO(ID2+LBAS,JA1+IBAS)
     &                     + PAB2*PCD2*RR(M, 7)*DENO(ID2+LBAS,JA2+IBAS)
C
              QXCH(IB1+JBAS,JC2+KBAS) = QXCH(IB1+JBAS,JC2+KBAS)
     &                     + PAB1*PCD2*RR(M,14)*DENO(ID1+LBAS,JA1+IBAS)
     &                     + PAB2*PCD2*RR(M, 6)*DENO(ID1+LBAS,JA2+IBAS)
     &                     + PAB1*PCD1*RR(M,13)*DENO(ID2+LBAS,JA1+IBAS)
     &                     + PAB2*PCD1*RR(M, 5)*DENO(ID2+LBAS,JA2+IBAS)
C
              QXCH(IB2+JBAS,JC1+KBAS) = QXCH(IB2+JBAS,JC1+KBAS)
     &                     + PAB2*PCD1*RR(M,12)*DENO(ID1+LBAS,JA1+IBAS)
     &                     + PAB1*PCD1*RR(M, 4)*DENO(ID1+LBAS,JA2+IBAS)
     &                     + PAB2*PCD2*RR(M,11)*DENO(ID2+LBAS,JA1+IBAS)
     &                     + PAB1*PCD2*RR(M, 3)*DENO(ID2+LBAS,JA2+IBAS)

C
              QXCH(IB2+JBAS,JC2+KBAS) = QXCH(IB2+JBAS,JC2+KBAS)
     &                     + PAB2*PCD2*RR(M,10)*DENO(ID1+LBAS,JA1+IBAS)
     &                     + PAB1*PCD2*RR(M, 2)*DENO(ID1+LBAS,JA2+IBAS)
     &                     + PAB2*PCD1*RR(M, 9)*DENO(ID2+LBAS,JA1+IBAS)
     &                     + PAB1*PCD1*RR(M, 1)*DENO(ID2+LBAS,JA2+IBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       EIGHTH IFLG BATCH (EXCHANGE)
        IF(IFLG(8).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              QXCH(IC1+KBAS,JB1+JBAS) = QXCH(IC1+KBAS,JB1+JBAS)
     &                     +           RR(M, 1)*DENO(JA1+IBAS,ID1+LBAS)
     &                     +           RR(M, 2)*DENO(JA1+IBAS,ID2+LBAS)
     &                     +           RR(M, 9)*DENO(JA2+IBAS,ID1+LBAS)
     &                     +           RR(M,10)*DENO(JA2+IBAS,ID2+LBAS)
C
              QXCH(IC1+KBAS,JB2+JBAS) = QXCH(IC1+KBAS,JB2+JBAS)
     &                     +           RR(M, 5)*DENO(JA1+IBAS,ID1+LBAS)
     &                     +           RR(M, 6)*DENO(JA1+IBAS,ID2+LBAS)
     &                     +           RR(M,13)*DENO(JA2+IBAS,ID1+LBAS)
     &                     +           RR(M,14)*DENO(JA2+IBAS,ID2+LBAS)
C
              QXCH(IC2+KBAS,JB1+JBAS) = QXCH(IC2+KBAS,JB1+JBAS)
     &                     +           RR(M, 3)*DENO(JA1+IBAS,ID1+LBAS)
     &                     +           RR(M, 4)*DENO(JA1+IBAS,ID2+LBAS)
     &                     +           RR(M,11)*DENO(JA2+IBAS,ID1+LBAS)
     &                     +           RR(M,12)*DENO(JA2+IBAS,ID2+LBAS)
C
              QXCH(IC2+KBAS,JB2+JBAS) = QXCH(IC2+KBAS,JB2+JBAS)
     &                     +           RR(M, 7)*DENO(JA1+IBAS,ID1+LBAS)
     &                     +           RR(M, 8)*DENO(JA1+IBAS,ID2+LBAS)
     &                     +           RR(M,15)*DENO(JA2+IBAS,ID1+LBAS)
     &                     +           RR(M,16)*DENO(JA2+IBAS,ID2+LBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       NINTH IFLG BATCH (EXCHANGE)
        IF(IFLG(9).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C  
              QXCH(IA1+IBAS,JD1+LBAS) = QXCH(IA1+IBAS,JD1+LBAS)
     &                     +           RR(M, 1)*DENO(IC1+KBAS,JB1+JBAS)
     &                     +           RR(M, 5)*DENO(IC1+KBAS,JB2+JBAS)
     &                     +           RR(M, 3)*DENO(IC2+KBAS,JB1+JBAS)
     &                     +           RR(M, 7)*DENO(IC2+KBAS,JB2+JBAS)
C
              QXCH(IA1+IBAS,JD2+LBAS) = QXCH(IA1+IBAS,JD2+LBAS)
     &                     +           RR(M, 2)*DENO(IC1+KBAS,JB1+JBAS)
     &                     +           RR(M, 6)*DENO(IC1+KBAS,JB2+JBAS)
     &                     +           RR(M, 4)*DENO(IC2+KBAS,JB1+JBAS)
     &                     +           RR(M, 8)*DENO(IC2+KBAS,JB2+JBAS)
C
              QXCH(IA2+IBAS,JD1+LBAS) = QXCH(IA2+IBAS,JD1+LBAS)
     &                     +           RR(M, 9)*DENO(IC1+KBAS,JB1+JBAS)
     &                     +           RR(M,13)*DENO(IC1+KBAS,JB2+JBAS)
     &                     +           RR(M,11)*DENO(IC2+KBAS,JB1+JBAS)
     &                     +           RR(M,15)*DENO(IC2+KBAS,JB2+JBAS)
C
              QXCH(IA2+IBAS,JD2+LBAS) = QXCH(IA2+IBAS,JD2+LBAS)
     &                     +           RR(M,10)*DENO(IC1+KBAS,JB1+JBAS)
     &                     +           RR(M,14)*DENO(IC1+KBAS,JB2+JBAS)
     &                     +           RR(M,12)*DENO(IC2+KBAS,JB1+JBAS)
     &                     +           RR(M,16)*DENO(IC2+KBAS,JB2+JBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       TENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(10).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C
              QXCH(IB1+JBAS,JD1+LBAS) = QXCH(IB1+JBAS,JD1+LBAS)
     &                     + PAB1*     RR(M,13)*DENO(IC1+KBAS,JA1+IBAS)
     &                     + PAB2*     RR(M, 5)*DENO(IC1+KBAS,JA2+IBAS)
     &                     + PAB1*     RR(M,15)*DENO(IC2+KBAS,JA1+IBAS)
     &                     + PAB2*     RR(M, 7)*DENO(IC2+KBAS,JA2+IBAS)
C
              QXCH(IB1+JBAS,JD2+LBAS) = QXCH(IB1+JBAS,JD2+LBAS)
     &                     + PAB1*     RR(M,14)*DENO(IC1+KBAS,JA1+IBAS)
     &                     + PAB2*     RR(M, 6)*DENO(IC1+KBAS,JA2+IBAS)
     &                     + PAB1*     RR(M,16)*DENO(IC2+KBAS,JA1+IBAS)
     &                     + PAB2*     RR(M, 8)*DENO(IC2+KBAS,JA2+IBAS)
C
              QXCH(IB2+JBAS,JD1+LBAS) = QXCH(IB2+JBAS,JD1+LBAS)
     &                     + PAB2*     RR(M, 9)*DENO(IC1+KBAS,JA1+IBAS)
     &                     + PAB1*     RR(M, 1)*DENO(IC1+KBAS,JA2+IBAS)
     &                     + PAB2*     RR(M,11)*DENO(IC2+KBAS,JA1+IBAS)
     &                     + PAB1*     RR(M, 3)*DENO(IC2+KBAS,JA2+IBAS)
C
              QXCH(IB2+JBAS,JD2+LBAS) = QXCH(IB2+JBAS,JD2+LBAS)
     &                     + PAB2*     RR(M,10)*DENO(IC1+KBAS,JA1+IBAS)
     &                     + PAB1*     RR(M, 2)*DENO(IC1+KBAS,JA2+IBAS)
     &                     + PAB2*     RR(M,12)*DENO(IC2+KBAS,JA1+IBAS)
     &                     + PAB1*     RR(M, 4)*DENO(IC2+KBAS,JA2+IBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C       ELEVENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(11).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C
              QXCH(ID1+LBAS,JB1+JBAS) = QXCH(ID1+LBAS,JB1+JBAS)
     &                     +      PCD1*RR(M, 4)*DENO(JA1+IBAS,IC1+KBAS)
     &                     +      PCD2*RR(M, 2)*DENO(JA1+IBAS,IC2+KBAS)
     &                     +      PCD1*RR(M,12)*DENO(JA2+IBAS,IC1+KBAS)
     &                     +      PCD2*RR(M,10)*DENO(JA2+IBAS,IC2+KBAS)
C
              QXCH(ID1+LBAS,JB2+JBAS) = QXCH(ID1+LBAS,JB2+JBAS)
     &                     +      PCD1*RR(M, 8)*DENO(JA1+IBAS,IC1+KBAS)
     &                     +      PCD2*RR(M, 6)*DENO(JA1+IBAS,IC2+KBAS)
     &                     +      PCD1*RR(M,16)*DENO(JA2+IBAS,IC1+KBAS)
     &                     +      PCD2*RR(M,14)*DENO(JA2+IBAS,IC2+KBAS)
C
              QXCH(ID2+LBAS,JB1+JBAS) = QXCH(ID2+LBAS,JB1+JBAS)
     &                     +      PCD2*RR(M, 3)*DENO(JA1+IBAS,IC1+KBAS)
     &                     +      PCD1*RR(M, 1)*DENO(JA1+IBAS,IC2+KBAS)
     &                     +      PCD2*RR(M,11)*DENO(JA2+IBAS,IC1+KBAS)
     &                     +      PCD1*RR(M, 9)*DENO(JA2+IBAS,IC2+KBAS)
C
              QXCH(ID2+LBAS,JB2+JBAS) = QXCH(ID2+LBAS,JB2+JBAS)
     &                     +      PCD2*RR(M, 7)*DENO(JA1+IBAS,IC1+KBAS)
     &                     +      PCD1*RR(M, 5)*DENO(JA1+IBAS,IC2+KBAS)
     &                     +      PCD2*RR(M,15)*DENO(JA2+IBAS,IC1+KBAS)
     &                     +      PCD1*RR(M,13)*DENO(JA2+IBAS,IC2+KBAS)
C
            ENDDO
          ENDDO
        ENDIF
C
C     SKIPPING POINT FOR CLOSED SYSTEMS
399   CONTINUE
C
C     CLOSE ALL LOOPS OVER BASIS FUNCTIONS A,B,C,D
3000  CONTINUE
2999  CONTINUE
2500  CONTINUE
2200  CONTINUE
2000  CONTINUE
1000  CONTINUE
C
C**********************************************************************C
C     COMPLETE CONSTRUCTION OF ALL MATRICES BY CONJUGATION.            C
C**********************************************************************C
C
C     LOOP OVER LOWER TRIANGLE OF EACH TT' BLOCK
      DO J=1,NDIM-NSHIFT
        DO I=1,J
C
C         SMALL-COMPONENT ADDRESSES
          K = I + NSHIFT
          L = J + NSHIFT
C
C         SKIP DIAGONAL PARTS OF EACH SUB-BLOCK
          IF(ICNBAS(I).NE.ICNBAS(J)) GOTO 400
          IF(KQNBAS(I).NE.KQNBAS(J)) GOTO 400
          IF(IABS(MQNBAS(I)).NE.IABS(MQNBAS(J))) GOTO 400
          GOTO 401
400       CONTINUE
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF LL BLOCK
          GDIR(I,J) = GDIR(I,J) + DCONJG(GDIR(J,I))
          GDIR(J,I) =             DCONJG(GDIR(I,J))
          GXCH(I,J) = GXCH(I,J) + DCONJG(GXCH(J,I))
          GXCH(J,I) =             DCONJG(GXCH(I,J))
          QDIR(I,J) = QDIR(I,J) + DCONJG(QDIR(J,I))
          QDIR(J,I) =             DCONJG(QDIR(I,J))
          QXCH(I,J) = QXCH(I,J) + DCONJG(QXCH(J,I))
          QXCH(J,I) =             DCONJG(QXCH(I,J))
C
C         IF HMLTN = 'NORL' SKIP THE NEXT FEW CALCULATIONS
          IF(HMLTN.EQ.'NORL') GOTO 401
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF SS BLOCK
          GDIR(K,L) = GDIR(K,L) + DCONJG(GDIR(L,K))
          GDIR(L,K) =             DCONJG(GDIR(K,L))
          GXCH(K,L) = GXCH(K,L) + DCONJG(GXCH(L,K))
          GXCH(L,K) =             DCONJG(GXCH(K,L))
          QDIR(K,L) = QDIR(K,L) + DCONJG(QDIR(L,K))
          QDIR(L,K) =             DCONJG(QDIR(K,L))
          QXCH(K,L) = QXCH(K,L) + DCONJG(QXCH(L,K))
          QXCH(L,K) =             DCONJG(QXCH(K,L))
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF LS BLOCK
          GXCH(I,L) = GXCH(I,L) + DCONJG(GXCH(L,I))
          GXCH(L,I) =             DCONJG(GXCH(I,L))
          QXCH(I,L) = QXCH(I,L) + DCONJG(QXCH(L,I))
          QXCH(L,I) =             DCONJG(QXCH(I,L))
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF SL BLOCK
          GXCH(K,J) = GXCH(K,J) + DCONJG(GXCH(J,K))
          GXCH(J,K) =             DCONJG(GXCH(K,J))
          QXCH(K,J) = QXCH(K,J) + DCONJG(QXCH(J,K))
          QXCH(J,K) =             DCONJG(QXCH(K,J))
C
401       CONTINUE
        ENDDO
      ENDDO
C
C     MULTIPLY OPEN MATRIX BY ANGULAR COEFFICIENTS
      DO J=1,NDIM
        DO I=1,NDIM
          QDIR(I,J) = ACFF*QDIR(I,J)
          QXCH(I,J) = BCFF*QXCH(I,J)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     SPECIFY OPEN-SHELL COUPLING MATRIX R AND ADD TO COULOMB MATRIX   C
C -------------------------------------------------------------------- C
C               R = [S.D(O).Q + Q.D(O).S]  (RSCF 89)                   C
C**********************************************************************C
C         
      DO J=1,NDIM
        DO K=1,NDIM
          T1(K) = DCMPLX(0.0D0,0.0D0)
          T2(K) = DCMPLX(0.0D0,0.0D0)
          DO L=1,NDIM
            T1(K) = T1(K) + DENT(K,L)*QDIR(L,J) - DENT(K,L)*QDIR(L,J)
            T2(K) = T2(K) + DENT(K,L)*OVAP(L,J)
          ENDDO
        ENDDO
C
        DO I=1,NDIM
          T3(I) = DCMPLX(0.0D0,0.0D0)
          DO K=1,NDIM
            T3(I) = T3(I) + T1(K)*OVAP(I,K) 
     &                    + T2(K)*QDIR(I,K) - T2(K)*QXCH(I,K)
          ENDDO
        ENDDO
C
C       ADD THE PROJECTOR AND Q-MATRIX TO THE COULOMB MATRIX  (RSCF 91)
C       DO I=1,NDIM 
C         FOCK(I,J) = FOCK(I,J) - QMAT(I,J) + T3(I)
C       ENDDO
C       DFNOTE: OKAY... SO I HAVE TO FIX ENERGY TERM
C
        DO I=1,NDIM
          QDIR(I,J) = QDIR(I,J) - T3(I)
          QXCH(I,J) = QXCH(I,J) + T3(I)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ERI(RR,XYZ,KQN,MQN,EXPT,NFUNS,ITQN,I,J,IEAB,IECD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                        EEEEEEEE RRRRRRR  IIII                        C
C                        EE       RR    RR  II                         C
C                        EE       RR    RR  II                         C
C                        EEEEEE   RR    RR  II                         C
C                        EE       RRRRRRR   II                         C
C                        EE       RR    RR  II                         C
C                        EEEEEEEE RR    RR IIII                        C
C                                                                      C
C -------------------------------------------------------------------- C
C  ERI GENERATES BLOCKS OF ELECTRON REPULSION INTEGRALS OVER           C
C  KINETICALLY BALANCED G-SPINOR BASIS FUNCTIONS. THE DENSITIES ARE    C
C  EXPANDED IN A BASIS OF HERMITE GAUSSIANS AND THE INTEGRALS ARE      C
C  GENERATED USING THE MCMURCHIE-DAVIDSON ALGORITHM.                   C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    XYZ(3,4)  - COORDINATES OF THE 4 NUCLEAR CENTERS.                 C
C    KQN(4)    - KQN QUANTUM NUMBERS OF THE CENTERS.                   C
C    MQN(4)    - |MQN| QUANTUM NUMBERS OF THE CENTERS.                 C
C    EXPT(I,4) - EXPONENTS IN BASIS FUNCTION BLOCK FOR A,B,C,D         C
C    NFUNS(J)  - NUMBER OF FUNCTIONS ON CENTER J.                      C
C    ITQN(2)   - SPINOR LABEL PAIR FOR AB (I=1) AND CD (I=2).          C
C                ITQN(I) = {LL,LS,SL,SS}.                              C
C    I,J       - COMPONENT LABEL INDEX FOR BASIS FUNCTION PAIR ON AB . C
C    IEAB,IECD - 0 DON'T RECALCULATE E(AB/CD)-COEFFICIENTS             C
C                1 DO    RECALCULATE E(AB/CD)-COEFFICIENTS             C
C  OUTPUT:                                                             C
C    RR(MB2,16) - ERI'S FOR BLOCK AB, ALL 16 MQN PROJECTION COMBOS.    C
C -------------------------------------------------------------------- C
C  DFNOTE: THERE ARE NO SHORTCUTS FOR INUCAB OR INUCCD CONDITIONS.     C
C          IF IATOM = 1 (4 BASIS FUNCTIONS ON SAME CENTER), SHOULD     C
C          EVALUATE WITH RACAH ALGEBRA INSTEAD.                        C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 CONE
      COMPLEX*16 RR(MB2,16),Q1(MB2),Q2(MB2)
      COMPLEX*16 EAB11(MB2,MEQ),ECD11(MB2,MEQ),GAB11(MB2,MEQ),
     &           EAB21(MB2,MEQ),ECD21(MB2,MEQ),GAB21(MB2,MEQ)
C
      DIMENSION KQN(4),LQN(4),MQN(4),ITQN(2),NFUNS(4),EXPT(MBS,4)
      DIMENSION PQ(MB2,3),APH(MB2),PRE(MB2),RC(MB2,MRC),XYZ(3,4)
      DIMENSION IAB11(MEQ),IAB21(MEQ),ICD11(MEQ),ICD21(MEQ),IRC(MRC)
C
      SAVE EAB11,EAB21,ECD11,ECD21
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &             IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
      DATA ROOTPI,SENS/1.7724538509055160D0,1.0D-14/
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     EVALUATE LQNS FOR BASIS FUNCTIONS A, B, C, D
      DO N=1,4
        IF(KQN(N).LT.0) THEN
          LQN(N) =-KQN(N)-1
        ELSE
          LQN(N) = KQN(N)
        ENDIF
      ENDDO
C
C     COUNTER FOR BRA OR KET OVERLAPS
      IALTAB = 1
      IALTCD =-1
C
C     CALCULATE FINITE SUM TRUNCATION VALUES DEPENDING ON OVERLAP TYPE
      IF(ITQN(1).EQ.1.AND.ITQN(2).EQ.1) THEN
        LAMAB = LQN(1) + LQN(2)
        LAMCD = LQN(3) + LQN(4)
      ELSEIF(ITQN(1).EQ.1.AND.ITQN(2).EQ.4) THEN
        LAMAB = LQN(1) + LQN(2)
        LAMCD = LQN(3) + LQN(4) + 2
      ELSEIF(ITQN(1).EQ.4.AND.ITQN(2).EQ.1) THEN
        LAMAB = LQN(1) + LQN(2) + 2
        LAMCD = LQN(3) + LQN(4)
      ELSEIF(ITQN(1).EQ.4.AND.ITQN(2).EQ.4) THEN
        LAMAB = LQN(1) + LQN(2) + 2
        LAMCD = LQN(3) + LQN(4) + 2
      ELSE
        WRITE(6, *) 'In ERI: incorrect call.'
        WRITE(7, *) 'In ERI: incorrect call.'
        STOP
      ENDIF
C
C     NUMBER OF FINITE EXPANSION INDICES FOR EACH CENTER
      NTUVAB = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
      NTUVCD = (LAMCD+1)*(LAMCD+2)*(LAMCD+3)/6
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB  = NFUNS(1)*NFUNS(2)
      MAXCD  = NFUNS(3)*NFUNS(4)
C
C**********************************************************************C
C     IF ASKED TO RECALCULATE E(AB) COEFFICIENTS, DO THIS FIRST        C
C**********************************************************************C
C
      IF(IEAB.EQ.1) THEN
C
        IF(ITQN(1).EQ.1) THEN
          CALL CPU_TIME(TDM1)
          IF(IEQS.EQ.0) THEN
           CALL EMAKELL(EAB11,EAB21,EXPT,XYZ,KQN,MQN,NFUNS,IALTAB,1,2,0)
          ELSEIF(IEQS.EQ.1) THEN
            DO ITUV=1,NTUVAB
              IAD = IABLL + (ITUV-1)*MAXAB
              DO M=1,MAXAB
                EAB11(M,ITUV) = DCMPLX(E0LLFL(IAD+M,1),E0LLFL(IAD+M,2))
                EAB21(M,ITUV) = DCMPLX(E0LLFL(IAD+M,3),E0LLFL(IAD+M,4))
              ENDDO
            ENDDO
          ENDIF
          CALL CPU_TIME(TDM2)
          TELL = TELL + TDM2 - TDM1
        ELSEIF(ITQN(1).EQ.4) THEN
          CALL CPU_TIME(TDM1)
          IF(IEQS.EQ.0) THEN
           CALL EMAKESS(EAB11,EAB21,EXPT,XYZ,KQN,MQN,NFUNS,IALTAB,1,2,0)
          ELSEIF(IEQS.EQ.1) THEN
            DO ITUV=1,NTUVAB
              IAD = IABSS + (ITUV-1)*MAXAB
              DO M=1,MAXAB
                EAB11(M,ITUV) = DCMPLX(E0SSFL(IAD+M,1),E0SSFL(IAD+M,2))
                EAB21(M,ITUV) = DCMPLX(E0SSFL(IAD+M,3),E0SSFL(IAD+M,4))
              ENDDO
            ENDDO
          ENDIF
          CALL CPU_TIME(TDM2)
          TESS = TESS + TDM2 - TDM1
        ENDIF
C
C       SCREENING PROCEDURE: TEST MAGNITUDES OF E-COEFFICIENT LISTS
        DO IAB=1,NTUVAB
C
C         11 OVERLAP (AB PAIRS)
          TEST = CDASUM(MAXAB,EAB11(1,IAB))
          IF(TEST.LE.SENS) THEN
            IAB11(IAB) = 0
          ELSE
            IAB11(IAB) = 1
          ENDIF
C
C         21 OVERLAP (AB PAIRS)
          TEST = CDASUM(MAXAB,EAB21(1,IAB))
          IF(TEST.LE.SENS) THEN
            IAB21(IAB) = 0
          ELSE
            IAB21(IAB) = 1
          ENDIF
        
        ENDDO
C
C       DO NOT CALCULATE AGAIN UNTIL ASKED EXTERNALLY
        IEAB = 0
C
      ENDIF
C
C**********************************************************************C
C     IF ASKED TO RECALCULATE E(CD) COEFFICIENTS, DO THIS NEXT         C
C**********************************************************************C
C
      IF(IECD.EQ.1) THEN
C
        IF(ITQN(2).EQ.1) THEN
          CALL CPU_TIME(TDM1)
          IF(IEQS.EQ.0) THEN
           CALL EMAKELL(ECD11,ECD21,EXPT,XYZ,KQN,MQN,NFUNS,IALTCD,3,4,0)
          ELSEIF(IEQS.EQ.1) THEN
            DO ITUV=1,NTUVCD
              IAD = ICDLL + (ITUV-1)*MAXCD
              DO M=1,MAXCD
                ECD11(M,ITUV) = DCMPLX(E0LLFL(IAD+M,5),E0LLFL(IAD+M,6))
                ECD21(M,ITUV) = DCMPLX(E0LLFL(IAD+M,7),E0LLFL(IAD+M,8))
              ENDDO
            ENDDO
          ENDIF
          CALL CPU_TIME(TDM2)
          TELL = TELL + TDM2 - TDM1
        ELSEIF(ITQN(2).EQ.4) THEN
          CALL CPU_TIME(TDM1)
          IF(IEQS.EQ.0) THEN
           CALL EMAKESS(ECD11,ECD21,EXPT,XYZ,KQN,MQN,NFUNS,IALTCD,3,4,0)
          ELSEIF(IEQS.EQ.1) THEN
            DO ITUV=1,NTUVCD
              IAD = ICDSS + (ITUV-1)*MAXCD
              DO M=1,MAXCD
                ECD11(M,ITUV) = DCMPLX(E0SSFL(IAD+M,5),E0SSFL(IAD+M,6))
                ECD21(M,ITUV) = DCMPLX(E0SSFL(IAD+M,7),E0SSFL(IAD+M,8))
              ENDDO
            ENDDO
          ENDIF
          CALL CPU_TIME(TDM2)
          TESS = TESS + TDM2 - TDM1
        ENDIF
C
C       SCREENING PROCEDURE: TEST MAGNITUDES OF E-COEFFICIENT LISTS
        DO ICD=1,NTUVCD
C
C         11 OVERLAP (CD PAIRS)
          TEST = CDASUM(MAXCD,ECD11(1,ICD))
          IF(TEST.LE.SENS) THEN
            ICD11(ICD) = 0
          ELSE
            ICD11(ICD) = 1
          ENDIF
C
C         21 OVERLAP (CD PAIRS)
          TEST = CDASUM(MAXCD,ECD21(1,ICD))
          IF(TEST.LE.SENS) THEN
            ICD21(ICD) = 0
          ELSE
            ICD21(ICD) = 1
          ENDIF
C
        ENDDO
C
C       DO NOT CALCULATE AGAIN UNTIL ASKED EXTERNALLY
        IECD = 0
C
      ENDIF
C
C**********************************************************************C
C     R-INTEGRAL EVALUATION                                            C
C**********************************************************************C
C
C     GAUSSIAN OVERLAP VALUES  
      EIJ = EXPT(I,1) + EXPT(J,2)
      PX  = (XYZ(1,1)*EXPT(I,1) + XYZ(1,2)*EXPT(J,2))/EIJ
      PY  = (XYZ(2,1)*EXPT(I,1) + XYZ(2,2)*EXPT(J,2))/EIJ
      PZ  = (XYZ(3,1)*EXPT(I,1) + XYZ(3,2)*EXPT(J,2))/EIJ
C
      M = 0
      DO KBAS=1,NFUNS(3)
        DO LBAS=1,NFUNS(4)
          M = M+1         
          EKL = EXPT(KBAS,3) + EXPT(LBAS,4)
          QX  = (XYZ(1,3)*EXPT(KBAS,3) + XYZ(1,4)*EXPT(LBAS,4))/EKL
          QY  = (XYZ(2,3)*EXPT(KBAS,3) + XYZ(2,4)*EXPT(LBAS,4))/EKL
          QZ  = (XYZ(3,3)*EXPT(KBAS,3) + XYZ(3,4)*EXPT(LBAS,4))/EKL
          EIJKL   = EIJ+EKL
          APH(M)  = EIJ*EKL/EIJKL
          PQ(M,1) = QX-PX
          PQ(M,2) = QY-PY
          PQ(M,3) = QZ-PZ
          EMX     = DSQRT(EIJ+EKL)*EIJ*EKL
          PRE(M)  = 2.0D0*(ROOTPI**5)/EMX
        ENDDO
      ENDDO
C
      MAXCD = NFUNS(3)*NFUNS(4)
C
      CALL CPU_TIME (TDM1)
      CALL RMAKE(RC,PQ,APH,MAXCD,LAMAB+LAMCD)
      CALL CPU_TIME(TDM2)
      IF(ITQN(1).EQ.1.AND.ITQN(2).EQ.1) THEN
        TRLL = TRLL + TDM2 - TDM1
      ELSEIF(ITQN(1).EQ.4.AND.ITQN(2).EQ.4) THEN
        TRSS = TRSS + TDM2 - TDM1
      ELSE
        TRLS = TRLS + TDM2 - TDM1
      ENDIF
C
C     INITIALIZE ARRAY TO IMPLEMENT SPARSITY IN R-VECTOR
      LAMABCD  = LAMAB + LAMCD
      NTUVABCD = (LAMABCD+1)*(LAMABCD+2)*(LAMABCD+3)/6
C      
      DO NRC=1,NTUVABCD
        TEST = RDASUM(MAXCD,RC(1,NRC),1)
        IF(TEST.LE.SENS) THEN
          IRC(NRC) = 0
        ELSE
          IRC(NRC) = 1
        ENDIF
      ENDDO
C
C**********************************************************************C
C     CONSTRUCT INTERMEDIATE MATRICES FOR MCMURCHIE-DAVIDSON           C
C**********************************************************************C
C
      DO IAB=1,NTUVAB
C
        IAB1 = 0
        IAB2 = 0
C
C       SCREENING MARKERS
        IF(IAB11(IAB).EQ.1) THEN
          IAB1 = 1
        ENDIF
C
        IF(IAB21(IAB).EQ.1) THEN
          IAB2 = 1
        ENDIF
C
C       INITIALISE IF ANY E-COEFF. FOR THIS {IA,IB,IC} PASSES TEST
        IF(IAB1+IAB2.GT.0) THEN
          DO M=1,MAXCD
            GAB11(M,IAB) = DCMPLX(0.0D0,0.0D0)
            GAB21(M,IAB) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDIF
C
C       CALCULATE OVERALL ADDRESS AND ADD TO THE G-ARRAY BINS
        DO ICD=1,NTUVCD
          IRABCD = INABCD(IVEC(IAB)+IVEC(ICD),JVEC(IAB)+JVEC(ICD),
     &                                        KVEC(IAB)+KVEC(ICD))
C
C         SKIP THIS STEP IF THE R-INTEGRAL IS NOT BIG ENOUGH
          IF(IRC(IRABCD).EQ.0) GOTO 798
C
          IF(ICD11(ICD).EQ.1) THEN
            DO M=1,MAXCD
              GAB11(M,IAB) = GAB11(M,IAB) + ECD11(M,ICD)*RC(M,IRABCD)
            ENDDO
          ENDIF
C
          IF(ICD21(ICD).EQ.1) THEN
            DO M=1,MAXCD
              GAB21(M,IAB) = GAB21(M,IAB) + ECD21(M,ICD)*RC(M,IRABCD)
            ENDDO
          ENDIF
C 
798       CONTINUE
        ENDDO
      ENDDO
C
C**********************************************************************C
C     GENERATE ALL POSSIBLE TWO-ELECTRON INTEGRALS FROM THE            C
C     EAB COEFFICIENTS AND THE G-ARRAYS                                C
C**********************************************************************C
C
C     CALCULATE PHASES FOR BASIS FUNCTION OVERLAP COMBINATIONS
      P1 = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
      P2 = DFLOAT((-1)**((MQN(3)-MQN(4))/2))
      P1 = P1*DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
      P2 = P2*DFLOAT((KQN(3)*KQN(4))/IABS(KQN(3)*KQN(4)))
      P3 = P1*P2
C
C**********************************************************************C
C     INTEGRAL BATCH 1: ( - - || - - )                                 C
C**********************************************************************C
C
      NTUVAB = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
      IJ     = (I-1)*NFUNS(2) + J
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY
      DO M=1,MAXCD
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {IA,IB,IC}
      DO IAB=1,NTUVAB
        IF(IAB11(IAB).EQ.1) THEN
          DO M=1,MAXCD
            Q1(M) = Q1(M) +      EAB11(IJ,IAB)*DREAL(GAB11(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB11(IJ,IAB)*DIMAG(GAB11(M,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO M=1,MAXCD
        RR(M,1 ) =    (Q1(M)+Q2(M))*PRE(M)
        RR(M,4 ) = P2*(Q1(M)-Q2(M))*PRE(M)
        RR(M,13) = P3*DCONJG(RR(M,4))
        RR(M,16) = P3*DCONJG(RR(M,1))
      ENDDO
C
C**********************************************************************C
C     INTEGRAL BATCH 2: ( - - || + - )                                 C
C**********************************************************************C
C
      NTUVAB = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
      IJ     = (I-1)*NFUNS(2) + J
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY
      DO M=1,MAXCD
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {IA,IB,IC}
      DO IAB=1,NTUVAB
        IF(IAB11(IAB).EQ.1) THEN
          DO M=1,MAXCD
            Q1(M) = Q1(M) +      EAB11(IJ,IAB)*DREAL(GAB21(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB11(IJ,IAB)*DIMAG(GAB21(M,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO M=1,MAXCD
        RR(M,3 ) =    (Q1(M)+Q2(M))*PRE(M)
        RR(M,2 ) =-P2*(Q1(M)-Q2(M))*PRE(M)
        RR(M,15) =-P3*DCONJG(RR(M,2))
        RR(M,14) =-P3*DCONJG(RR(M,3))
      ENDDO
C
C**********************************************************************C
C     INTEGRAL BATCH 3: ( + - || - - )                                 C
C**********************************************************************C
C
      NTUVAB = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
      IJ     = (I-1)*NFUNS(2) + J
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY
      DO M=1,MAXCD
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {IA,IB,IC}
      DO IAB=1,NTUVAB
       IF(IAB21(IAB).EQ.1) THEN
          DO M=1,MAXCD
            Q1(M) = Q1(M) +      EAB21(IJ,IAB)*DREAL(GAB11(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB21(IJ,IAB)*DIMAG(GAB11(M,IAB))
          ENDDO
       ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO M=1,MAXCD
        RR(M,9 ) =    (Q1(M)+Q2(M))*PRE(M)
        RR(M,12) = P2*(Q1(M)-Q2(M))*PRE(M)
        RR(M,5 ) =-P3*DCONJG(RR(M,12))
        RR(M,8 ) =-P3*DCONJG(RR(M,9 ))
      ENDDO
C
C**********************************************************************C
C     INTEGRAL BATCH 4: ( + - || + - )                                 C
C**********************************************************************C
C
      NTUVAB = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
      IJ  = (I-1)*NFUNS(2) + J
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY
      DO M=1,MAXCD
        Q1(M) = DCMPLX(0.0D0,0.0D0)
        Q2(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {IA,IB,IC}
      DO IAB=1,NTUVAB
        IF(IAB21(IAB).EQ.1) THEN
          DO M=1,MAXCD
            Q1(M) = Q1(M) +      EAB21(IJ,IAB)*DREAL(GAB21(M,IAB))
            Q2(M) = Q2(M) + CONE*EAB21(IJ,IAB)*DIMAG(GAB21(M,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE RR ARRAY
      DO M=1,MAXCD
        RR(M,11) =    (Q1(M)+Q2(M))*PRE(M)
        RR(M,10) =-P2*(Q1(M)-Q2(M))*PRE(M)
        RR(M,7 ) = P3*DCONJG(RR(M,10))
        RR(M,6 ) = P3*DCONJG(RR(M,11))
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE BREIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C               BBBBBBB  RRRRRRR  EEEEEEEE IIII TTTTTTTT               C
C               BB    BB RR    RR EE        II     TT                  C
C               BB    BB RR    RR EE        II     TT                  C
C               BBBBBBB  RR    RR EEEEEE    II     TT                  C
C               BB    BB RRRRRRR  EE        II     TT                  C
C               BB    BB RR    RR EE        II     TT                  C
C               BBBBBBB  RR    RR EEEEEEEE IIII    TT                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  BREIT GENERATES ELECTRON INTERACTION INTEGRALS IN BATCHES AND       C
C  CALCULATES THE SCF BREIT MATRIX (B). INTEGRAL SYMMETRY IS PARTIALLY C
C  PARTIALLY EXPLOITED, BUT WITH ROOM FOR IMPROVEMENT.                 C
C  (GEOMETRIC SYMM, R-INT SYMM, E-COEFF SYMM).                         C
C -------------------------------------------------------------------- C
C  DFNOTE: THIS ROUTINE COULD BENEFIT FROM PARALLELISATION -- OPENMP.  C
C          LOOK INTO OPEN SHELL EXTENSIONS AND CASES WITH ZERO DIRECT. C
C**********************************************************************C    
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 BR(MB2,16)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NFUNS(4)
      DIMENSION INDEX(MCT,-MKP:MKP,2*(MKP+1)*MKP),RMX(MB2)
      DIMENSION IBINDX(MCT,MKP)
C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/EILS/EILSFL(MFL,24),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/IBSC/IBSCR(MB2),IBMAP(MB2)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/QNMS/ICNBAS(MDM),KQNBAS(MDM),MQNBAS(MDM)
      COMMON/SHLL/ALPHA,BETA,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN,NOELEC
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     INITIALISE STORAGE MATRICES
      DO I=1,NDIM
        DO J=1,NDIM
          BDIR(I,J) = DCMPLX(0.0D0,0.0D0)
          BXCH(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INDEXING ROUTINE: SO WE CAN SET BASIS FUNCTION LABELS BY BLOCK   C
C**********************************************************************C
      ICOUNT = 0
      IKOUNT = 0
C
C     LOOP OVER NUCLEAR CENTERS
      DO ICNT=1,NCNT
C
C       LOOP OVER KAPPA VALUES FOR THIS NUCLEAR CENTER
        DO KN=1,NKAP(ICNT)
C
C         IMPORT KAPPA, LQN, MAXIMUM MQN AND NUMBER OF BASIS FUNCTIONS
          KAPPA = KVALS(KN,ICNT)
          MJMAX = 2*IABS(KAPPA)-1
          IF(KAPPA.GT.0) THEN
            LQNN = KAPPA
          ELSE
            LQNN =-KAPPA-1
          ENDIF
          NFUN = NFUNCT(LQNN+1,ICNT)
C
C         RECORD CUMULATIVE INDEX COUNTER FOR NUMBER OF BASIS FUNCTIONS
          IBINDX(ICNT,KN) = IKOUNT
          IKOUNT          = IKOUNT + NFUN
C
C         LOOP OVER MQN VALUES AND RECORD INDEX (MORE THAN IKOUNT)
          DO MJ=1,MJMAX,2
            ICOUNT               = ICOUNT+1
            INDEX(ICNT,KAPPA,MJ) = ICOUNT
            LENIQ                = ICOUNT
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTERS A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     LOOP OVER CENTER A
      DO 1000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 1000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 0
        ENDIF
C
C     LOOP OVER KQN(A) VALUES
      DO 1000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNS(1) = NFUNCT(LQN(1)+1,ICNTA)
        DO IBAS=1,NFUNS(1)
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 1000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNS(2) = NFUNCT(LQN(2)+1,ICNTB)
        DO JBAS=1,NFUNS(2)
          EXPT(JBAS,2) = EXPSET(JBAS,LQN(2)+1,ICNTB)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 1000 MA=1,IABS(KQN(1))
        MJA    = 2*MA-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 1000 MB=1,IABS(KQN(2))
        MJB    = 2*MB-1
        MQN(2) = MJB
C
C     CALCULATE NEW BLOCK OF E(AB) COEFFS AT NEXT OPPORTUNITY
      IEAB  = 1
      IABLS = IADILS(ICNTA,ICNTB,KA,KB,MA,MB)
C
C**********************************************************************C
C     SECOND LAYER OF LOOPS, OVER CENTERS C AND D (USE INDEX 3000)     C
C**********************************************************************C
C
C     LOOP OVER CENTER C
      DO 2000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER C
        XYZ(1,3) = COORD(1,ICNTC)
        XYZ(2,3) = COORD(2,ICNTC)
        XYZ(3,3) = COORD(3,ICNTC)
C
C     LOOP OVER CENTER D
      DO 2000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER D
        XYZ(1,4) = COORD(1,ICNTD)
        XYZ(2,4) = COORD(2,ICNTD)
        XYZ(3,4) = COORD(3,ICNTD)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER C AND D
        IF(ICNTC.EQ.ICNTD) THEN
          INUCCD = 1
        ELSE
          INUCCD = 0
        ENDIF
C
C       PARAMETER FOR ATOMIC OR MULTICNTRE INTEGRAL
        IF(INUCAB*INUCCD.EQ.1.AND.ICNTA.EQ.ICNTC) THEN
          IATOM = 1
        ELSE
          IATOM = 0
        ENDIF
C
C       SKIP MULTI-CENTER CONTRIBUTIONS IN STAGE 1
        IF(IATOM.EQ.0.AND.ILEV.EQ.1) GOTO 2200
C
C     LOOP OVER KQN(C) VALUES
      DO 2500 KC=1,NKAP(ICNTC)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR C
        KQN(3) = KVALS(KC,ICNTC)
        IF(KQN(3).GT.0) THEN
          LQN(3) = KQN(3)
        ELSE
          LQN(3) =-KQN(3)-1
        ENDIF
C
        NFUNS(3) = NFUNCT(LQN(3)+1,ICNTC)
        DO KBAS=1,NFUNS(3)
          EXPT(KBAS,3) = EXPSET(KBAS,LQN(3)+1,ICNTC)
        ENDDO
C
C     LOOP OVER KQN(D) VALUES
      DO 2500 KD=1,NKAP(ICNTD)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR D
        KQN(4) = KVALS(KD,ICNTD)
        IF(KQN(4).GT.0) THEN
         LQN(4) = KQN(4)
        ELSE
         LQN(4) =-KQN(4)-1
        ENDIF
C
        NFUNS(4) = NFUNCT(LQN(4)+1,ICNTD)
        DO LBAS=1,NFUNS(4)
          EXPT(LBAS,4) = EXPSET(LBAS,LQN(4)+1,ICNTD)
        ENDDO
C
C     LOOP OVER |MQN(C)| VALUES
      DO 2500 MC=1,IABS(KQN(3))
        MJC    = 2*MC-1
        MQN(3) = MJC
C
C     LOOP OVER |MQN(D)| VALUES
      DO 2500 MD=1,IABS(KQN(4))
        MJD    = 2*MD-1
        MQN(4) = MJD
C
C     CALCULATE NEW BLOCK OF E(CD) COEFFS AT NEXT OPPORTUNITY
      IECD  = 1
      ICDLS = IADILS(ICNTC,ICNTD,KC,KD,MC,MD)
C
C**********************************************************************C
C     FOR THIS CHOICE OF A,B,C AND D, COMPUTE INDICES                  C
C**********************************************************************C
C
C     STARTING INDEX VALUES
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
C     LENIQ IS THE NUMBER OF BLOCKS TO BE CALCULATED
      IQ12 = (IQ1-1)*LENIQ + IQ2
      IQ34 = (IQ3-1)*LENIQ + IQ4
C
C     FURTHER DEFINE STARTING ADDRESSES FOR {ABCD} BASIS OVERLAPS
      IA1 = LARGE(ICNTA,KA,MJA  )
      IB1 = LARGE(ICNTB,KB,MJB  ) + NSHIFT
      IC1 = LARGE(ICNTC,KC,MJC  )
      ID1 = LARGE(ICNTD,KD,MJD  ) + NSHIFT
C
      IA2 = LARGE(ICNTA,KA,MJA+1)
      IB2 = LARGE(ICNTB,KB,MJB+1) + NSHIFT
      IC2 = LARGE(ICNTC,KC,MJC+1)
      ID2 = LARGE(ICNTD,KD,MJD+1) + NSHIFT
C
      JA1 = LARGE(ICNTA,KA,MJA  )
      JB1 = LARGE(ICNTB,KB,MJB  ) + NSHIFT
      JC1 = LARGE(ICNTC,KC,MJC  )
      JD1 = LARGE(ICNTD,KD,MJD  ) + NSHIFT
C
      JA2 = LARGE(ICNTA,KA,MJA+1)
      JB2 = LARGE(ICNTB,KB,MJB+1) + NSHIFT
      JC2 = LARGE(ICNTC,KC,MJC+1)
      JD2 = LARGE(ICNTD,KD,MJD+1) + NSHIFT
C
C     NOT SURE WHAT THESE ARE FOR, OR WHETHER WE REALLY NEED COMMON/BSCR
      IAA = IBINDX(ICNTA,KA)
      JBB = IBINDX(ICNTB,KB) + NSHIFT
C
      ICC = IBINDX(ICNTC,KC)
      JDD = IBINDX(ICNTD,KD) + NSHIFT
C
C**********************************************************************C
C     SKIP BATCHES WHICH CONFORM TO INTEGRAL SYMMETRIES                C
C**********************************************************************C
C
C     DIATOMIC MOLECULE SELECTION RULES
      IF(NCNT.LE.2) THEN
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) GOTO 2998
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) GOTO 2998
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) GOTO 2998
        GOTO 2999
      ENDIF
2998  CONTINUE
C
      IF(IQ12.LT.IQ34) GOTO 2999
C
C**********************************************************************C
C     FOR THIS CHOICE OF A,B,C AND D, COMPUTE PHASES                   C
C**********************************************************************C
C
C     CALCULATE KQN PHASE FACTORS FOR PERMUTING INTEGRALS.
C     NB! OPPOSITE PHASE AS IN THE LL/SS CASE SEEN IN SCF
      IF((KQN(1)*KQN(2)).GT.0) THEN 
        PKAB =-1.0D0
      ELSE
        PKAB = 1.0D0
      ENDIF
C
      IF((KQN(3)*KQN(4)).GT.0) THEN 
        PKCD =-1.0D0
      ELSE
        PKCD = 1.0D0
      ENDIF
C
C     CALCULATE MQN PHASE FACTORS FOR PERMUTING INTEGRALS
      PMAB1 = DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PMAB2 = DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      PMCD1 = DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PMCD2 = DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
      F1 = PKAB*PMAB1
      F2 = PKAB*PMAB2
      G1 = PKCD*PMCD1
      G2 = PKCD*PMCD2
C
C**********************************************************************C
C     THIRD LAYER OF LOOPS, OVER BASIS FUNCTIONS A AND B (USE 4000)    C
C**********************************************************************C
C
      DO 3000 IBAS=1,NFUNS(1)
      DO 3000 JBAS=1,NFUNS(2)
C
C       OVERRIDE SCREENING TESTS FOR THE SAKE OF CERTAINTY
        DO M=1,NFUNS(3)*NFUNS(4)
          IBMAP(M) = M
          IBSCR(M) = 1
        ENDDO
C
C       GENERATE A BATCH OF BREIT INTERACTION INTEGRALS
        CALL BII(BR,XYZ,KQN,MQN,EXPT,NFUNS,IBAS,JBAS,IEAB,IECD)
C
C       FIRST BATCH
        M = 0
        DO KBAS=1,NFUNS(3)
          DO LBAS=1,NFUNS(4)
            M = M+1
C
            BDIR(IA1+IBAS,JB1+JBAS) = BDIR(IA1+IBAS,JB1+JBAS)
     &                            +    BR(M, 1)*DENT(IC1+KBAS,JD1+LBAS)
     &                            +    BR(M, 2)*DENT(IC1+KBAS,JD2+LBAS)
     &                            +    BR(M, 3)*DENT(IC2+KBAS,JD1+LBAS)
     &                            +    BR(M, 4)*DENT(IC2+KBAS,JD2+LBAS)
     &                            + G1*BR(M, 4)*DENT(JD1+LBAS,IC1+KBAS)
     &                            + G2*BR(M, 2)*DENT(JD1+LBAS,IC2+KBAS)
     &                            + G2*BR(M, 3)*DENT(JD2+LBAS,IC1+KBAS)
     &                            + G1*BR(M, 1)*DENT(JD2+LBAS,IC2+KBAS)
C
            BDIR(IA1+IBAS,JB2+JBAS) = BDIR(IA1+IBAS,JB2+JBAS)
     &                            +    BR(M, 5)*DENT(IC1+KBAS,JD1+LBAS)
     &                            +    BR(M, 6)*DENT(IC1+KBAS,JD2+LBAS)
     &                            +    BR(M, 7)*DENT(IC2+KBAS,JD1+LBAS)
     &                            +    BR(M, 8)*DENT(IC2+KBAS,JD2+LBAS)
     &                            + G1*BR(M, 8)*DENT(JD1+LBAS,IC1+KBAS)
     &                            + G2*BR(M, 6)*DENT(JD1+LBAS,IC2+KBAS)
     &                            + G2*BR(M, 7)*DENT(JD2+LBAS,IC1+KBAS)
     &                            + G1*BR(M, 5)*DENT(JD2+LBAS,IC2+KBAS)
C
            BDIR(IA2+IBAS,JB1+JBAS) = BDIR(IA2+IBAS,JB1+JBAS)
     &                            +    BR(M, 9)*DENT(IC1+KBAS,JD1+LBAS)
     &                            +    BR(M,10)*DENT(IC1+KBAS,JD2+LBAS)
     &                            +    BR(M,11)*DENT(IC2+KBAS,JD1+LBAS)
     &                            +    BR(M,12)*DENT(IC2+KBAS,JD2+LBAS)
     &                            + G1*BR(M,12)*DENT(JD1+LBAS,IC1+KBAS)
     &                            + G2*BR(M,10)*DENT(JD1+LBAS,IC2+KBAS)
     &                            + G2*BR(M,11)*DENT(JD2+LBAS,IC1+KBAS)
     &                            + G1*BR(M, 9)*DENT(JD2+LBAS,IC2+KBAS)
C
            BDIR(IA2+IBAS,JB2+JBAS) = BDIR(IA2+IBAS,JB2+JBAS)
     &                            +    BR(M,13)*DENT(IC1+KBAS,JD1+LBAS)
     &                            +    BR(M,14)*DENT(IC1+KBAS,JD2+LBAS)
     &                            +    BR(M,15)*DENT(IC2+KBAS,JD1+LBAS)
     &                            +    BR(M,16)*DENT(IC2+KBAS,JD2+LBAS)
     &                            + G1*BR(M,16)*DENT(JD1+LBAS,IC1+KBAS)
     &                            + G2*BR(M,14)*DENT(JD1+LBAS,IC2+KBAS)
     &                            + G2*BR(M,15)*DENT(JD2+LBAS,IC1+KBAS)
     &                            + G1*BR(M,13)*DENT(JD2+LBAS,IC2+KBAS)
C
          ENDDO
        ENDDO
C
C       SECOND BATCH
        DO LBAS=1,NFUNS(4)
          DO KBAS=1,NFUNS(3)
            M = (KBAS-1)*NFUNS(4) + LBAS
C
            BXCH(IA1+IBAS,JD1+LBAS) = BXCH(IA1+IBAS,JD1+LBAS)
     &                            +    BR(M, 1)*DENT(IC1+KBAS,JB1+JBAS)
     &                            +    BR(M, 5)*DENT(IC1+KBAS,JB2+JBAS)
     &                            +    BR(M, 3)*DENT(IC2+KBAS,JB1+JBAS)
     &                            +    BR(M, 7)*DENT(IC2+KBAS,JB2+JBAS)
C
            BXCH(IA1+IBAS,JD2+LBAS) = BXCH(IA1+IBAS,JD2+LBAS)
     &                            +    BR(M, 2)*DENT(IC1+KBAS,JB1+JBAS)
     &                            +    BR(M, 6)*DENT(IC1+KBAS,JB2+JBAS)
     &                            +    BR(M, 4)*DENT(IC2+KBAS,JB1+JBAS)
     &                            +    BR(M, 8)*DENT(IC2+KBAS,JB2+JBAS)
C
            BXCH(IA2+IBAS,JD1+LBAS) = BXCH(IA2+IBAS,JD1+LBAS)
     &                            +    BR(M, 9)*DENT(IC1+KBAS,JB1+JBAS)
     &                            +    BR(M,13)*DENT(IC1+KBAS,JB2+JBAS)
     &                            +    BR(M,11)*DENT(IC2+KBAS,JB1+JBAS)
     &                            +    BR(M,15)*DENT(IC2+KBAS,JB2+JBAS)
C
            BXCH(IA2+IBAS,JD2+LBAS) = BXCH(IA2+IBAS,JD2+LBAS)
     &                            +    BR(M,10)*DENT(IC1+KBAS,JB1+JBAS)
     &                            +    BR(M,14)*DENT(IC1+KBAS,JB2+JBAS)
     &                            +    BR(M,12)*DENT(IC2+KBAS,JB1+JBAS)
     &                            +    BR(M,16)*DENT(IC2+KBAS,JB2+JBAS)
C
          ENDDO
        ENDDO
C
C       THIRD BATCH
        M = 0
        DO KBAS=1,NFUNS(3)
          DO LBAS=1,NFUNS(4)
            M = M+1
C
            BXCH(IA1+IBAS,JC1+KBAS) = BXCH(IA1+IBAS,JC1+KBAS)
     &                            + G1*BR(M, 4)*DENT(ID1+LBAS,JB1+JBAS)
     &                            + G1*BR(M, 8)*DENT(ID1+LBAS,JB2+JBAS)
     &                            + G2*BR(M, 3)*DENT(ID2+LBAS,JB1+JBAS)
     &                            + G2*BR(M, 7)*DENT(ID2+LBAS,JB2+JBAS)
C
            BXCH(IA1+IBAS,JC2+KBAS) = BXCH(IA1+IBAS,JC2+KBAS)
     &                            + G2*BR(M, 2)*DENT(ID1+LBAS,JB1+JBAS)
     &                            + G2*BR(M, 6)*DENT(ID1+LBAS,JB2+JBAS)
     &                            + G1*BR(M, 1)*DENT(ID2+LBAS,JB1+JBAS)
     &                            + G1*BR(M, 5)*DENT(ID2+LBAS,JB2+JBAS)
C
            BXCH(IA2+IBAS,JC1+KBAS) = BXCH(IA2+IBAS,JC1+KBAS)
     &                            + G1*BR(M,12)*DENT(ID1+LBAS,JB1+JBAS)
     &                            + G1*BR(M,16)*DENT(ID1+LBAS,JB2+JBAS)
     &                            + G2*BR(M,11)*DENT(ID2+LBAS,JB1+JBAS)
     &                            + G2*BR(M,15)*DENT(ID2+LBAS,JB2+JBAS)
C
            BXCH(IA2+IBAS,JC2+KBAS) = BXCH(IA2+IBAS,JC2+KBAS)
     &                            + G2*BR(M,10)*DENT(ID1+LBAS,JB1+JBAS)
     &                            + G2*BR(M,14)*DENT(ID1+LBAS,JB2+JBAS)
     &                            + G1*BR(M, 9)*DENT(ID2+LBAS,JB1+JBAS)
     &                            + G1*BR(M,13)*DENT(ID2+LBAS,JB2+JBAS)
C
          ENDDO
        ENDDO
C
C       FOURTH BATCH
        DO LBAS=1,NFUNS(4)
          DO KBAS=1,NFUNS(3)
            M = (KBAS-1)*NFUNS(4) + LBAS
C
            BXCH(IB1+JBAS,JD1+LBAS) = BXCH(IB1+JBAS,JD1+LBAS)
     &                            + F1*BR(M,13)*DENT(IC1+KBAS,JA1+IBAS)
     &                            + F2*BR(M, 5)*DENT(IC1+KBAS,JA2+IBAS)
     &                            + F1*BR(M,15)*DENT(IC2+KBAS,JA1+IBAS)
     &                            + F2*BR(M, 7)*DENT(IC2+KBAS,JA2+IBAS)
C
            BXCH(IB1+JBAS,JD2+LBAS) = BXCH(IB1+JBAS,JD2+LBAS)
     &                            + F1*BR(M,14)*DENT(IC1+KBAS,JA1+IBAS)
     &                            + F2*BR(M, 6)*DENT(IC1+KBAS,JA2+IBAS)
     &                            + F1*BR(M,16)*DENT(IC2+KBAS,JA1+IBAS)
     &                            + F2*BR(M, 8)*DENT(IC2+KBAS,JA2+IBAS)
C
            BXCH(IB2+JBAS,JD1+LBAS) = BXCH(IB2+JBAS,JD1+LBAS)
     &                            + F2*BR(M, 9)*DENT(IC1+KBAS,JA1+IBAS)
     &                            + F1*BR(M, 1)*DENT(IC1+KBAS,JA2+IBAS)
     &                            + F2*BR(M,11)*DENT(IC2+KBAS,JA1+IBAS)
     &                            + F1*BR(M, 3)*DENT(IC2+KBAS,JA2+IBAS)
C
            BXCH(IB2+JBAS,JD2+LBAS) = BXCH(IB2+JBAS,JD2+LBAS)
     &                            + F2*BR(M,10)*DENT(IC1+KBAS,JA1+IBAS)
     &                            + F1*BR(M, 2)*DENT(IC1+KBAS,JA2+IBAS)
     &                            + F2*BR(M,12)*DENT(IC2+KBAS,JA1+IBAS)
     &                            + F1*BR(M, 4)*DENT(IC2+KBAS,JA2+IBAS)
C
          ENDDO
        ENDDO
C
C ***   THERE IS AN EXTRA TERM FOR ELEMENTS SATISFYING IQ12 = IQ34
        IF(IQ12.EQ.IQ34) GOTO 100
C
C       FIFTH BATCH
        M = 0
        DO KBAS=1,NFUNS(3)
          DO LBAS=1,NFUNS(4)
            M = M+1
C
            BDIR(IC1+KBAS,JD1+LBAS) = BDIR(IC1+KBAS,JD1+LBAS)
     &                            +    BR(M, 1)*DENT(IA1+IBAS,JB1+JBAS)
     &                            +    BR(M, 5)*DENT(IA1+IBAS,JB2+JBAS)
     &                            +    BR(M, 9)*DENT(IA2+IBAS,JB1+JBAS)
     &                            +    BR(M,13)*DENT(IA2+IBAS,JB2+JBAS)
     &                            + F1*BR(M,13)*DENT(JB1+JBAS,IA1+IBAS)
     &                            + F2*BR(M, 5)*DENT(JB1+JBAS,IA2+IBAS)
     &                            + F2*BR(M, 9)*DENT(JB2+JBAS,IA1+IBAS)
     &                            + F1*BR(M, 1)*DENT(JB2+JBAS,IA2+IBAS)
C
            BDIR(IC1+KBAS,JD2+LBAS) = BDIR(IC1+KBAS,JD2+LBAS)
     &                            +    BR(M, 2)*DENT(IA1+IBAS,JB1+JBAS)
     &                            +    BR(M, 6)*DENT(IA1+IBAS,JB2+JBAS)
     &                            +    BR(M,10)*DENT(IA2+IBAS,JB1+JBAS)
     &                            +    BR(M,14)*DENT(IA2+IBAS,JB2+JBAS)
     &                            + F1*BR(M,14)*DENT(JB1+JBAS,IA1+IBAS)
     &                            + F2*BR(M, 6)*DENT(JB1+JBAS,IA2+IBAS)
     &                            + F2*BR(M,10)*DENT(JB2+JBAS,IA1+IBAS)
     &                            + F1*BR(M, 2)*DENT(JB2+JBAS,IA2+IBAS)
C
            BDIR(IC2+KBAS,JD1+LBAS) = BDIR(IC2+KBAS,JD1+LBAS)
     &                            +    BR(M, 3)*DENT(IA1+IBAS,JB1+JBAS)
     &                            +    BR(M, 7)*DENT(IA1+IBAS,JB2+JBAS)
     &                            +    BR(M,11)*DENT(IA2+IBAS,JB1+JBAS)
     &                            +    BR(M,15)*DENT(IA2+IBAS,JB2+JBAS)
     &                            + F1*BR(M,15)*DENT(JB1+JBAS,IA1+IBAS)
     &                            + F2*BR(M, 7)*DENT(JB1+JBAS,IA2+IBAS)
     &                            + F2*BR(M,11)*DENT(JB2+JBAS,IA1+IBAS)
     &                            + F1*BR(M, 3)*DENT(JB2+JBAS,IA2+IBAS)
C
            BDIR(IC2+KBAS,JD2+LBAS) = BDIR(IC2+KBAS,JD2+LBAS)
     &                            +    BR(M, 4)*DENT(IA1+IBAS,JB1+JBAS)
     &                            +    BR(M, 8)*DENT(IA1+IBAS,JB2+JBAS)
     &                            +    BR(M,12)*DENT(IA2+IBAS,JB1+JBAS)
     &                            +    BR(M,16)*DENT(IA2+IBAS,JB2+JBAS)
     &                            + F1*BR(M,16)*DENT(JB1+JBAS,IA1+IBAS)
     &                            + F2*BR(M, 8)*DENT(JB1+JBAS,IA2+IBAS)
     &                            + F2*BR(M,12)*DENT(JB2+JBAS,IA1+IBAS)
     &                            + F1*BR(M, 4)*DENT(JB2+JBAS,IA2+IBAS)
C
          ENDDO
        ENDDO
C
C       SIXTH BATCH
        M = 0
        DO KBAS=1,NFUNS(3)
          DO LBAS=1,NFUNS(4)
            M = M+1
C
            BXCH(IC1+KBAS,JA1+IBAS) = BXCH(IC1+KBAS,JA1+IBAS)
     &                            + F1*BR(M,13)*DENT(JB1+JBAS,ID1+LBAS)
     &                            + F1*BR(M,14)*DENT(JB1+JBAS,ID2+LBAS)
     &                            + F2*BR(M, 9)*DENT(JB2+JBAS,ID1+LBAS)
     &                            + F2*BR(M,10)*DENT(JB2+JBAS,ID2+LBAS)
C
            BXCH(IC1+KBAS,JA2+IBAS) = BXCH(IC1+KBAS,JA2+IBAS)
     &                            + F2*BR(M, 5)*DENT(JB1+JBAS,ID1+LBAS)
     &                            + F2*BR(M, 6)*DENT(JB1+JBAS,ID2+LBAS)
     &                            + F1*BR(M, 1)*DENT(JB2+JBAS,ID1+LBAS)
     &                            + F1*BR(M, 2)*DENT(JB2+JBAS,ID2+LBAS)
C
            BXCH(IC2+KBAS,JA1+IBAS) = BXCH(IC2+KBAS,JA1+IBAS)
     &                            + F1*BR(M,15)*DENT(JB1+JBAS,ID1+LBAS)
     &                            + F1*BR(M,16)*DENT(JB1+JBAS,ID2+LBAS)
     &                            + F2*BR(M,11)*DENT(JB2+JBAS,ID1+LBAS)
     &                            + F2*BR(M,12)*DENT(JB2+JBAS,ID2+LBAS)
C
            BXCH(IC2+KBAS,JA2+IBAS) = BXCH(IC2+KBAS,JA2+IBAS)
     &                            + F2*BR(M, 7)*DENT(JB1+JBAS,ID1+LBAS)
     &                            + F2*BR(M, 8)*DENT(JB1+JBAS,ID2+LBAS)
     &                            + F1*BR(M, 3)*DENT(JB2+JBAS,ID1+LBAS)
     &                            + F1*BR(M, 4)*DENT(JB2+JBAS,ID2+LBAS)
C
          ENDDO
        ENDDO
C
C       SEVENTH BATCH
        DO LBAS=1,NFUNS(4)
          DO KBAS=1,NFUNS(3)
            M = (KBAS-1)*NFUNS(4) + LBAS
C
            BXCH(ID1+LBAS,JB1+JBAS) = BXCH(ID1+LBAS,JB1+JBAS)
     &                            + G1*BR(M, 4)*DENT(JA1+IBAS,IC1+KBAS)
     &                            + G2*BR(M, 2)*DENT(JA1+IBAS,IC2+KBAS)
     &                            + G1*BR(M,12)*DENT(JA2+IBAS,IC1+KBAS)
     &                            + G2*BR(M,10)*DENT(JA2+IBAS,IC2+KBAS)
C
            BXCH(ID1+LBAS,JB2+JBAS) = BXCH(ID1+LBAS,JB2+JBAS)
     &                            + G1*BR(M, 8)*DENT(JA1+IBAS,IC1+KBAS)
     &                            + G2*BR(M, 6)*DENT(JA1+IBAS,IC2+KBAS)
     &                            + G1*BR(M,16)*DENT(JA2+IBAS,IC1+KBAS)
     &                            + G2*BR(M,14)*DENT(JA2+IBAS,IC2+KBAS)
C
            BXCH(ID2+LBAS,JB1+JBAS) = BXCH(ID2+LBAS,JB1+JBAS)
     &                            + G2*BR(M, 3)*DENT(JA1+IBAS,IC1+KBAS)
     &                            + G1*BR(M, 1)*DENT(JA1+IBAS,IC2+KBAS)
     &                            + G2*BR(M,11)*DENT(JA2+IBAS,IC1+KBAS)
     &                            + G1*BR(M, 9)*DENT(JA2+IBAS,IC2+KBAS)
C
            BXCH(ID2+LBAS,JB2+JBAS) = BXCH(ID2+LBAS,JB2+JBAS)
     &                            + G2*BR(M, 7)*DENT(JA1+IBAS,IC1+KBAS)
     &                            + G1*BR(M, 5)*DENT(JA1+IBAS,IC2+KBAS)
     &                            + G2*BR(M,15)*DENT(JA2+IBAS,IC1+KBAS)
     &                            + G1*BR(M,13)*DENT(JA2+IBAS,IC2+KBAS)
C
          ENDDO
        ENDDO
C
C       EIGHTH BATCH
        M = 0
        DO KBAS=1,NFUNS(3)
          DO LBAS=1,NFUNS(4)
            M = M+1
C
            BXCH(IC1+KBAS,JB1+JBAS) = BXCH(IC1+KBAS,JB1+JBAS)
     &                            +    BR(M, 1)*DENT(JA1+IBAS,ID1+LBAS)
     &                            +    BR(M, 2)*DENT(JA1+IBAS,ID2+LBAS)
     &                            +    BR(M, 9)*DENT(JA2+IBAS,ID1+LBAS)
     &                            +    BR(M,10)*DENT(JA2+IBAS,ID2+LBAS)
C
            BXCH(IC1+KBAS,JB2+JBAS) = BXCH(IC1+KBAS,JB2+JBAS)
     &                            +    BR(M, 5)*DENT(JA1+IBAS,ID1+LBAS)
     &                            +    BR(M, 6)*DENT(JA1+IBAS,ID2+LBAS)
     &                            +    BR(M,13)*DENT(JA2+IBAS,ID1+LBAS)
     &                            +    BR(M,14)*DENT(JA2+IBAS,ID2+LBAS)
C
            BXCH(IC2+KBAS,JB1+JBAS) = BXCH(IC2+KBAS,JB1+JBAS)
     &                            +    BR(M, 3)*DENT(JA1+IBAS,ID1+LBAS)
     &                            +    BR(M, 4)*DENT(JA1+IBAS,ID2+LBAS)
     &                            +    BR(M,11)*DENT(JA2+IBAS,ID1+LBAS)
     &                            +    BR(M,12)*DENT(JA2+IBAS,ID2+LBAS)
C
            BXCH(IC2+KBAS,JB2+JBAS) = BXCH(IC2+KBAS,JB2+JBAS)
     &                            +    BR(M, 7)*DENT(JA1+IBAS,ID1+LBAS)
     &                            +    BR(M, 8)*DENT(JA1+IBAS,ID2+LBAS)
     &                            +    BR(M,15)*DENT(JA2+IBAS,ID1+LBAS)
     &                            +    BR(M,16)*DENT(JA2+IBAS,ID2+LBAS)
C
          ENDDO
        ENDDO
C
100     CONTINUE
C
C**********************************************************************C
C     END OF BREIT MATRIX CONSTRUCTION                                 C
C**********************************************************************C
C
3000  CONTINUE
2999  CONTINUE
2500  CONTINUE
2200  CONTINUE
2000  CONTINUE
1000  CONTINUE
C
C**********************************************************************C
C     COMPLETE CONSTRUCTION OF BREIT MATRIX BY MATRIX CONJUGATION.     C
C**********************************************************************C
C
      DO I=1,NSHIFT
        DO J=NSHIFT+1,NDIM
          BDIR(J,I) = DCONJG(BDIR(I,J))
          BXCH(J,I) = DCONJG(BXCH(I,J))
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE BII(BR,XYZ,KQN,MQN,EXPT,NFUNS,IBAS,JBAS,IEAB,IECD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                          BBBBBBB  IIII IIII                          C
C                          BB    BB  II   II                           C
C                          BB    BB  II   II                           C
C                          BBBBBBB   II   II                           C
C                          BB    BB  II   II                           C
C                          BB    BB  II   II                           C
C                          BBBBBBB  IIII IIII                          C
C                                                                      C
C -------------------------------------------------------------------- C
C  BII GENERATES BLOCKS OF BREIT INTERACTION INTEGRALS OVER            C
C  KINETICALLY BALANCED G-SPINOR BASIS FUNCTIONS. THE DENSITIES ARE    C
C  EXPANDED IN A BASIS OF HERMITE GAUSSIANS AND THE INTEGRALS ARE      C
C  GENERATED USING THE MCMURCHIE-DAVIDSON ALGORITHM.                   C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    XYZ(3,4)  - COORDINATES OF THE 4 NUCLEAR CENTERS.                 C
C    KQN(4)    - KQN QUANTUM NUMBERS OF THE CENTERS.                   C
C    MQN(4)    - |MQN| QUANTUM NUMBERS OF THE CENTERS.                 C
C    EXPT(I,4) - EXPONENTS IN BASIS FUNCTION BLOCK FOR A,B,C,D         C
C    NFUNS(J)  - NUMBER OF FUNCTIONS ON CENTER J.                      C
C    I,J       - COMPONENT LABEL INDEX FOR BASIS FUNCTION PAIR ON AB . C
C    IEAB,IECD - 0 DON'T RECALCULATE E(AB/CD)-COEFFICIENTS             C
C                1 DO    RECALCULATE E(AB/CD)-COEFFICIENTS             C
C  OUTPUT:                                                             C
C    BR(MB2,16) - BII'S FOR BLOCK AB, ALL 16 MQN PROJECTION COMBOS.    C
C -------------------------------------------------------------------- C
C  DFNOTE: THERE ARE NO SHORTCUTS FOR INUCAB OR INUCCD CONDITIONS.     C
C          IF IATOM = 1 (4 BASIS FUNCTIONS ON SAME CENTER), SHOULD     C
C          EVALUATE WITH RACAH ALGEBRA INSTEAD!                        C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 CONE
      COMPLEX*16 BR(MB2,16),Q1(MB2),Q2(MB2)
      COMPLEX*16 EAB11(MB2,MEQ,3),EAB21(MB2,MEQ,3),BAB11(MB2,MEQ),
     &           ECD11(MB2,MEQ,3),ECD21(MB2,MEQ,3),BAB21(MB2,MEQ)
      COMPLEX*16 ELSAB11(MB2,4*MEQ*MEQ),ELSAB21(MB2,4*MEQ*MEQ),
     &           ELSCD11(MB2,4*MEQ*MEQ),ELSCD21(MB2,4*MEQ*MEQ)
C
      DIMENSION XYZ(3,4),KQN(4),LQN(4),MQN(4),NFUNS(4),IRC(MRC)
      DIMENSION RC(MB2,MRC),PQ(MB2,3),EXPT(MBS,4),APH(MB2),PRE(MB2)
      DIMENSION IAB11(MEQ,3),IAB21(MEQ,3),ICD11(MEQ,3),ICD21(MEQ,3)
C
      SAVE EAB11,EAB21,ECD11,ECD21
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/EILS/EILSFL(MFL,24),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/IBSC/IBSCR(MB2),IBMAP(MB2)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
      DATA ROOTPI5,SENS/1.7493418327624863D1,1.0D-14/
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     EVALUATE LQNS FOR BASIS FUNCTIONS A, B, C, D
      DO N=1,4
        IF(KQN(N).LT.0) THEN
         LQN(N) =-KQN(N)-1
        ELSE
         LQN(N) = KQN(N)
        ENDIF
      ENDDO
C      
C     MAXIMUM REQUIRED GTF CARTESIAN SUM
      LAMAB = LQN(1)+LQN(2)+1
      LAMCD = LQN(3)+LQN(4)+1
C
C     LENGTH OF EQ-COEFFICIENT LIST
      NTUVAB = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
      NTUVCD = (LAMCD+1)*(LAMCD+2)*(LAMCD+3)/6
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NFUNS(1)*NFUNS(2)
      MAXCD = NFUNS(3)*NFUNS(4)
C
C**********************************************************************C
C     IF ASKED TO RECALCULATE E(AB) COEFFICIENTS, DO THIS FIRST        C
C**********************************************************************C
C
      IF(IEAB.EQ.1) THEN
C
C       PHASE FACTOR FOR AB PAIRS
        IALTAB = 1
C
C       GENERATE ELS(AB) COEFFICIENTS
        CALL CPU_TIME(TDM1)
        IF(IEQS.EQ.0) THEN
          CALL EMAKEB3(EAB11,EAB21,EXPT,XYZ,KQN,MQN,NFUNS,IALTAB,1,2)
        ELSEIF(IEQS.EQ.1) THEN
          DO ITUV=1,NTUVAB
            IAD = IABLS + (ITUV-1)*MAXAB
            DO M=1,MAXAB
              EAB11(M,ITUV,1)=DCMPLX(EILSFL(IAD+M, 1),EILSFL(IAD+M, 2))
              EAB21(M,ITUV,1)=DCMPLX(EILSFL(IAD+M, 3),EILSFL(IAD+M, 4))
              EAB11(M,ITUV,2)=DCMPLX(EILSFL(IAD+M, 5),EILSFL(IAD+M, 6))
              EAB21(M,ITUV,2)=DCMPLX(EILSFL(IAD+M, 7),EILSFL(IAD+M, 8))
              EAB11(M,ITUV,3)=DCMPLX(EILSFL(IAD+M, 9),EILSFL(IAD+M,10))
              EAB21(M,ITUV,3)=DCMPLX(EILSFL(IAD+M,11),EILSFL(IAD+M,12))
            ENDDO
          ENDDO
        ENDIF
        CALL CPU_TIME(TDM2)
        TELS = TELS + TDM2 - TDM1
C
C       SCREENING PROCEDURE: TEST MAGNITUDES OF E-COEFFICIENT LISTS
        DO IQ=1,3
          DO IAB=1,NTUVAB
C
C           11 OVERLAP (AB PAIRS)
            TEST = CDASUM(MAXAB,EAB11(1,IAB,IQ))
            IF(TEST.LE.SENS) THEN
              IAB11(IAB,IQ) = 0
            ELSE
              IAB11(IAB,IQ) = 1
            ENDIF
C
C           21 OVERLAP (AB PAIRS)
            TEST = CDASUM(MAXAB,EAB21(1,IAB,IQ))
            IF(TEST.LE.SENS) THEN
              IAB21(IAB,IQ) = 0
            ELSE
              IAB21(IAB,IQ) = 1
            ENDIF
C
          ENDDO
        ENDDO
C
C       DO NOT CALCULATE AGAIN UNTIL ASKED EXTERNALLY
        IEAB = 0
C
      ENDIF
C
C**********************************************************************C
C     IF ASKED TO RECALCULATE E(CD) COEFFICIENTS, DO THIS NEXT         C
C**********************************************************************C
C
      IF(IECD.EQ.1) THEN
C
C       PHASE FACTOR FOR CD PAIRS
        IALTCD =-1
C
C       GENERATE ELS(CD) COEFFICIENTS
        CALL CPU_TIME(TDM1)
        IF(IEQS.EQ.0) THEN
          CALL EMAKEB3(ECD11,ECD21,EXPT,XYZ,KQN,MQN,NFUNS,IALTCD,3,4)
        ELSEIF(IEQS.EQ.1) THEN
          DO ITUV=1,NTUVCD
            IAD = ICDLS + (ITUV-1)*MAXCD
            DO M=1,MAXCD
              ECD11(M,ITUV,1)=DCMPLX(EILSFL(IAD+M,13),EILSFL(IAD+M,14))
              ECD21(M,ITUV,1)=DCMPLX(EILSFL(IAD+M,15),EILSFL(IAD+M,16))
              ECD11(M,ITUV,2)=DCMPLX(EILSFL(IAD+M,17),EILSFL(IAD+M,18))
              ECD21(M,ITUV,2)=DCMPLX(EILSFL(IAD+M,19),EILSFL(IAD+M,20))
              ECD11(M,ITUV,3)=DCMPLX(EILSFL(IAD+M,21),EILSFL(IAD+M,22))
              ECD21(M,ITUV,3)=DCMPLX(EILSFL(IAD+M,23),EILSFL(IAD+M,24))
            ENDDO
          ENDDO
        ENDIF
        CALL CPU_TIME(TDM2)
        TELS = TELS + TDM2 - TDM1
C 
C       SCREENING PROCEDURE: TEST MAGNITUDES OF E-COEFFICIENT LISTS
        DO IQ=1,3
          DO ICD=1,NTUVCD
C
C           11 OVERLAP (CD PAIRS)
            TEST = CDASUM(MAXCD,ECD11(1,ICD,IQ))
            IF(TEST.LE.SENS) THEN
              ICD11(ICD,IQ) = 0
            ELSE
              ICD11(ICD,IQ) = 1
            ENDIF
C
C           21 OVERLAP (CD PAIRS)
            TEST = CDASUM(MAXCD,ECD21(1,ICD,IQ))
            IF(TEST.LE.SENS) THEN
              ICD21(ICD,IQ) = 0
            ELSE
              ICD21(ICD,IQ) = 1
            ENDIF
C
          ENDDO
        ENDDO

C
C       DO NOT CALCULATE AGAIN UNTIL ASKED EXTERNALLY
        IECD = 0
C
      ENDIF
C
C**********************************************************************C
C     R-INTEGRAL EVALUATION                                            C
C**********************************************************************C
C
C     GAUSSIAN OVERLAP VALUES
      EIJ = EXPT(IBAS,1) + EXPT(JBAS,2)
      PX  = (XYZ(1,1)*EXPT(IBAS,1) + XYZ(1,2)*EXPT(JBAS,2))/EIJ
      PY  = (XYZ(2,1)*EXPT(IBAS,1) + XYZ(2,2)*EXPT(JBAS,2))/EIJ
      PZ  = (XYZ(3,1)*EXPT(IBAS,1) + XYZ(3,2)*EXPT(JBAS,2))/EIJ
C
C     M IS THE LONGER COUNTER, N IS SHORTENED DEPENDING ON IBSCR
      M = 0
      N = 0
      DO KBAS=1,NFUNS(3)
        DO LBAS=1,NFUNS(4)
          M = M+1
          IF(IBSCR(M).EQ.1) THEN
            N   = N+1
            EKL = EXPT(KBAS,3) + EXPT(LBAS,4)
            QX  = (XYZ(1,3)*EXPT(KBAS,3) + XYZ(1,4)*EXPT(LBAS,4))/EKL
            QY  = (XYZ(2,3)*EXPT(KBAS,3) + XYZ(2,4)*EXPT(LBAS,4))/EKL
            QZ  = (XYZ(3,3)*EXPT(KBAS,3) + XYZ(3,4)*EXPT(LBAS,4))/EKL
C
            EIJKL   = EIJ+EKL
            APH(N)  = EIJ*EKL/EIJKL
            PQ(N,1) = QX-PX
            PQ(N,2) = QY-PY
            PQ(N,3) = QZ-PZ
            EMX     = DSQRT(EIJ+EKL)*EIJ*EKL
            PRE(N)  = 2.0D0*ROOTPI5/EMX
          ENDIF
        ENDDO
      ENDDO
C
C     MAXN COUNTS # OVERLAPS THAT SURVIVE SCREENING
      MAXM = M
      MAXN = N
C
      CALL CPU_TIME(TDM1)
      CALL RMAKE(RC,PQ,APH,MAXN,LAMAB+LAMCD+2)
      CALL CPU_TIME(TDM2)
      TRBR = TRBR + TDM2 - TDM1
C
C     INITIALIZE ARRAY TO IMPLEMENT SPARSITY IN R-VECTOR
      LAMABCD  = LAMAB+LAMCD+2
      NTUVABCD = (LAMABCD+1)*(LAMABCD+2)*(LAMABCD+3)/6
C
      DO NRC=1,NTUVABCD
        TEST = RDASUM(MAXN,RC(1,NRC),1)
        IF(TEST.LE.SENS) THEN
          IRC(NRC) = 0
        ELSE
          IRC(NRC) = 1
        ENDIF
      ENDDO     
C
C     INITIALISE BR ARRAY
      DO ITG=1,16
        DO M=1,MAXCD
          BR(M,ITG) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CONSTRUCT INTERMEDIATE MATRICES FOR MCMURCHIE-DAVIDSON           C
C**********************************************************************C
C
C     BEGIN FIRST LOOP: CARTESIAN INDEX ICMP FOR CENTER AB (USE 6000)
      DO 6000 ICMP=1,3
C
C     CARTESIAN INDICES FOR {IA,IB,IC}
      IF(ICMP.EQ.1) THEN
        IDX = 1
        IDY = 0
        IDZ = 0
      ELSEIF(ICMP.EQ.2) THEN
        IDX = 0
        IDY = 1
        IDZ = 0
      ELSEIF(ICMP.EQ.3) THEN
        IDX = 0
        IDY = 0
        IDZ = 1
      ENDIF
C
C     INITIATE LOOP OVER ALL INDICES {IA,IB,IC} FOR AB
      DO IAB=1,NTUVAB
C
        IAB1 = 0
        IAB2 = 0
C
C       SCREENING MARKERS
        IF((IAB11(IAB,1)+IAB11(IAB,2)+IAB11(IAB,3)).NE.0) THEN
          IAB1 = 1
        ENDIF
C
        IF((IAB21(IAB,1)+IAB21(IAB,2)+IAB21(IAB,3)).NE.0) THEN
          IAB2 = 1
        ENDIF
C
C       INITIALISE IF ANY E-COEFF. FOR THIS {IA,IB,IC} PASSES TEST
        IF(IAB1+IAB2.GT.0) THEN
          DO N=1,MAXN
            BAB11(N,IAB) = DCMPLX(0.0D0,0.0D0)
            BAB21(N,IAB) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDIF
C
C       INITIATE LOOP OVER ALL INDICES {IA',IB',IC'} FOR CD
        DO JCMP=1,3
C
C         CARTESIAN INDICES FOR {IA',IB',IC'}
          IF(JCMP.EQ.1) THEN
            JDX = 1
            JDY = 0
            JDZ = 0
          ELSEIF(JCMP.EQ.2) THEN
            JDX = 0
            JDY = 1
            JDZ = 0
          ELSEIF(JCMP.EQ.3) THEN
            JDX = 0
            JDY = 0
            JDZ = 1
          ENDIF
C
C ***     SPECIAL CASE: CARTESIAN INDICES ARE EQUAL {BXX, BYY, BZZ}
          IF(ICMP.EQ.JCMP) THEN
C
C           BEGIN LOOP OVER ALL (CD) ADDRESSES FOR A GIVEN (AB) ADDRESS
            DO ICD=1,NTUVCD
C
C             OVERALL ADDRESS OF THIS {IA,IB,IC},{IA',IB',IC'}
              IRABCD = INABCD(IVEC(IAB)+IVEC(ICD),JVEC(IAB)+JVEC(ICD),
     &                                            KVEC(IAB)+KVEC(ICD))
C
C             SCREEN CALCULATION IF THE BIGGEST R-INTEGRAL IS TOO SMALL
              IF(IRC(IRABCD).EQ.0) GOTO 797
C
C             IF E11(CD) IS BIG ENOUGH, CONTINUE TO ADD ALL OVERLAPS
              IF(ICD11(ICD,JCMP).EQ.1) THEN
                DO N=1,MAXN
                  BAB11(N,IAB) = BAB11(N,IAB)
     &                         - ECD11(IBMAP(N),ICD,JCMP)*RC(N,IRABCD)
                ENDDO
              ENDIF
C
C             IF E21(CD) IS BIG ENOUGH, CONTINUE TO ADD ALL OVERLAPS
              IF(ICD21(ICD,JCMP).EQ.1) THEN
                DO N=1,MAXN
                  BAB21(N,IAB) = BAB21(N,IAB)
     &                         - ECD21(IBMAP(N),ICD,JCMP)*RC(N,IRABCD)
                ENDDO
              ENDIF
C
797         CONTINUE
            ENDDO
C
C ***     END CONDITIONAL OVER LIKE INDICES
          ENDIF
C
C         LOOP OVER FINITE SUM INDICES {IA',IB',IC'} FOR CD
          DO ICD=1,NTUVCD
C
            IF(JCMP.EQ.1) THEN
              RTP = DFLOAT(IVEC(IAB) + IVEC(ICD))
            ELSEIF(JCMP.EQ.2) THEN
              RTP = DFLOAT(JVEC(IAB) + JVEC(ICD))
            ELSE
              RTP = DFLOAT(KVEC(IAB) + KVEC(ICD))
            ENDIF
C
C           STARTING ADDRESSES FOR CARTESIAN COMPONENTS
            I1 = IVEC(IAB) + IVEC(ICD) + IDX + JDX
            J1 = JVEC(IAB) + JVEC(ICD) + IDY + JDY
            K1 = KVEC(IAB) + KVEC(ICD) + IDZ + JDZ
            IADR1 = INABCD(I1,J1,K1)
C
            I2 = IVEC(IAB) + IVEC(ICD) + IDX
            J2 = JVEC(IAB) + JVEC(ICD) + IDY
            K2 = KVEC(IAB) + KVEC(ICD) + IDZ
            IADR2 = INABCD(I2,J2,K2)
C
            I3 = IVEC(IAB) + IVEC(ICD) + IDX - JDX
            J3 = JVEC(IAB) + JVEC(ICD) + IDY - JDY
            K3 = KVEC(IAB) + KVEC(ICD) + IDZ - JDZ
            IF((I3.GE.0).AND.(J3.GE.0).AND.(K3.GE.0)) THEN
              IADR3 = INABCD(I3,J3,K3)
            ELSE
              IADR3 = 1
              RTP   = 0.0D0
            ENDIF
C
C           IF E11(CD) IS BIG ENOUGH, CONTINUE TO ADD ALL OVERLAPS
            IF(ICD11(ICD,JCMP).EQ.1) THEN
              DO N=1,MAXN
                T1 = 0.5D0*RC(N,IADR1)/APH(N)
                T2 = RC(N,IADR2)*PQ(N,JCMP)
                T3 = RC(N,IADR3)*RTP
                BAB11(N,IAB) = BAB11(N,IAB)
     &                         + ECD11(IBMAP(N),ICD,JCMP)*(T1-T2+T3)
              ENDDO
            ENDIF
C
C           IF E21(CD) IS BIG ENOUGH, CONTINUE TO ADD ALL OVERLAPS
            IF(ICD21(ICD,JCMP).EQ.1) THEN
              DO N=1,MAXN
                T1 = RC(N,IADR1)*0.5D0/APH(N)
                T2 = RC(N,IADR2)*PQ(N,JCMP)
                T3 = RC(N,IADR3)*RTP
                BAB21(N,IAB) = BAB21(N,IAB)
     &                       + ECD21(IBMAP(N),ICD,JCMP)*(T1-T2+T3)
              ENDDO
            ENDIF
C
C         END LOOP OVER INDICES {IA',IB',IC'} FOR CD
          ENDDO
C
C       END LOOP OVER CD INDEX {JX,JY,JZ}
        ENDDO
C
C     END LOOP OVER INDICES {IA,IB,IC} FOR AB
      ENDDO
C
C**********************************************************************C
C     GENERATE ALL POSSIBLE TWO-ELECTRON INTEGRALS FROM THE            C
C     EAB COEFFICIENTS AND THE G-ARRAYS                                C
C**********************************************************************C
C
C     CALCULATE PHASE FACTORS FOR MQN AND KQN COMBINATIONS (P1,P2,P3)
      P1 = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
      P2 = DFLOAT((-1)**((MQN(3)-MQN(4))/2))
C
      P1 =-(P1*DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2))))
      P2 =-(P2*DFLOAT((KQN(3)*KQN(4))/IABS(KQN(3)*KQN(4))))
C
      P3 = P1*P2
C
      IJ = (IBAS-1)*NFUNS(2) + JBAS
C
C**********************************************************************C
C     INTEGRAL BATCH 1: ( - - || - - )                                 C
C**********************************************************************C
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY 
      DO N=1,MAXN
        Q1(N) = DCMPLX(0.0D0,0.0D0)
        Q2(N) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {IA,IB,IC}
      DO IAB=1,NTUVAB
        IF(IAB11(IAB,ICMP).EQ.1) THEN
          DO N=1,MAXN
            Q1(N) = Q1(N) +      EAB11(IJ,IAB,ICMP)*DREAL(BAB11(N,IAB))
            Q2(N) = Q2(N) + CONE*EAB11(IJ,IAB,ICMP)*DIMAG(BAB11(N,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE BR ARRAY
      DO N=1,MAXN
        BR(N,1 ) = BR(N,1 ) +    (Q1(N)+Q2(N))*PRE(N)
        BR(N,4 ) = BR(N,4 ) + P2*(Q1(N)-Q2(N))*PRE(N)
        BR(N,13) = P3*DCONJG(BR(N,4))
        BR(N,16) = P3*DCONJG(BR(N,1))
      ENDDO
C
C**********************************************************************C
C     INTEGRAL BATCH 2: ( - - || + - )                                 C
C**********************************************************************C
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY 
      DO N=1,MAXN
        Q1(N) = DCMPLX(0.0D0,0.0D0)
        Q2(N) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {IA,IB,IC}
      DO IAB=1,NTUVAB
        IF(IAB11(IAB,ICMP).EQ.1) THEN
          DO N=1,MAXN
            Q1(N) = Q1(N) +      EAB11(IJ,IAB,ICMP)*DREAL(BAB21(N,IAB))
            Q2(N) = Q2(N) + CONE*EAB11(IJ,IAB,ICMP)*DIMAG(BAB21(N,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE BR ARRAY
      DO N=1,MAXN
        BR(N,3 ) = BR(N,3 ) +    (Q1(N)+Q2(N))*PRE(N)
        BR(N,2 ) = BR(N,2 ) - P2*(Q1(N)-Q2(N))*PRE(N)
        BR(N,15) =-P3*DCONJG(BR(N,2))
        BR(N,14) =-P3*DCONJG(BR(N,3))
      ENDDO
C
C**********************************************************************C
C     INTEGRAL BATCH 3: ( + - || - - )                                 C
C**********************************************************************C
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY 
      DO N=1,MAXN
        Q1(N) = DCMPLX(0.0D0,0.0D0)
        Q2(N) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {IA,IB,IC}
      DO IAB=1,NTUVAB
        IF(IAB21(IAB,ICMP).EQ.1) THEN
          DO N=1,MAXN
            Q1(N) = Q1(N) +      EAB21(IJ,IAB,ICMP)*DREAL(BAB11(N,IAB))
            Q2(N) = Q2(N) + CONE*EAB21(IJ,IAB,ICMP)*DIMAG(BAB11(N,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE BR ARRAY
      DO N=1,MAXN
        BR(N,9 ) = BR(N,9 ) +    (Q1(N)+Q2(N))*PRE(N)
        BR(N,12) = BR(N,12) + P2*(Q1(N)-Q2(N))*PRE(N)
        BR(N,5 ) =-P3*DCONJG(BR(N,12))
        BR(N,8 ) =-P3*DCONJG(BR(N,9 ))
      ENDDO
C
C**********************************************************************C
C     INTEGRAL BATCH 4: ( + - || + - )                                 C
C**********************************************************************C
C
C     EMPTY Q-ARRAYS TO IMPLEMENT SPARSITY 
      DO N=1,MAXN
        Q1(N) = DCMPLX(0.0D0,0.0D0)
        Q2(N) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     LOOP OVER ALL COMBINATIONS {IA,IB,IC}
      DO IAB=1,NTUVAB
        IF(IAB21(IAB,ICMP).EQ.1) THEN
          DO N=1,MAXN
            Q1(N) = Q1(N) +      EAB21(IJ,IAB,ICMP)*DREAL(BAB21(N,IAB))
            Q2(N) = Q2(N) + CONE*EAB21(IJ,IAB,ICMP)*DIMAG(BAB21(N,IAB))
          ENDDO
        ENDIF
      ENDDO
C
C     FILL THIS BATCH OF THE BR ARRAY
      DO N=1,MAXN
        BR(N,11) = BR(N,11) +    (Q1(N)+Q2(N))*PRE(N)
        BR(N,10) = BR(N,10) - P2*(Q1(N)-Q2(N))*PRE(N)
        BR(N,7 ) = P3*DCONJG(BR(N,10))
        BR(N,6 ) = P3*DCONJG(BR(N,11))
      ENDDO
C
C     END LOOP OVER INDICES {IX,IY,IZ}
6000  CONTINUE
C
C**********************************************************************C
C     BREIT OVERLAP ARRAY NOW FULLY CONSTRUCTED                        C
C**********************************************************************C
C
C     INCLUDE THE OUTSIDE FACTOR OF (1/2)
      DO ITG=1,16
        DO N=1,MAXN
          BR(IBMAP(N),ITG) = 0.5D0*BR(IBMAP(N),ITG)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION RDASUM(N,DX,INCX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       RRRRRRR  DDDDDDD     AA     SSSSSS  UU    UU MM       MM       C
C       RR    RR DD    DD   AAAA   SS    SS UU    UU MMM     MMM       C
C       RR    RR DD    DD  AA  AA  SS       UU    UU MMMM   MMMM       C
C       RR    RR DD    DD AA    AA  SSSSSS  UU    UU MM MM MM MM       C
C       RRRRRRR  DD    DD AAAAAAAA       SS UU    UU MM  MMM  MM       C
C       RR    RR DD    DD AA    AA SS    SS UU    UU MM   M   MM       C
C       RR    RR DDDDDDD  AA    AA  SSSSSS   UUUUUU  MM       MM       C
C                                                                      C
C -------------------------------------------------------------------- C
C  RDASUM RETURNS THE SUM OF MAGNITUDES OF A VECTOR OF DOUBLES DX(N).  C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    N     - NUMBER OF ELEMENTS IN INPUT VECTOR(S).                    C
C    DX    - DOUBLE PRECISION VECTOR WITH N ELEMENTS.                  C
C    INCX  - STORAGE SPACING BETWEEN ELEMENTS OF DX.                   C
C  OUTPUT:                                                             C
C    RDASUM - DOUBLE PRECISION RESULT (ZERO IF N < 0).                 C
C**********************************************************************C
C
      DIMENSION DX(N)
C
C     INITIALISE COUNTER
      RDASUM = 0.D0
C
C     IF VECTOR LENGTH IS ZERO, RESULT IS ZERO
      IF(N.LE.0) THEN
        RETURN
      ENDIF
C
C     DECISION TREE FOR INCREMENT STEP SIZE
      IF(INCX.EQ.1) THEN
        GOTO 20
      ENDIF

      NS = N*INCX
      DO I=1,NS,INCX
        RDASUM = RDASUM + DABS(DX(I))
      ENDDO
      RETURN
C
20    CONTINUE
C
C     CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
      M = MOD(N,6)
C
      IF(M.EQ.0) GOTO 40
C
      DO I=1,M
        RDASUM = RDASUM + DABS(DX(I))
      ENDDO
C      
      IF(N.LT.6) RETURN
C
   40 CONTINUE
C   
      MP1 = M+1
      DO I = MP1,N,6
        RDASUM = RDASUM + DABS(DX(I  )) + DABS(DX(I+1)) + DABS(DX(I+2))
     &                  + DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION CDASUM(N,DX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        CCCCCC  DDDDDDD     AA     SSSSSS  UU    UU MM       MM       C
C       CC    CC DD    DD   AAAA   SS    SS UU    UU MMM     MMM       C
C       CC       DD    DD  AA  AA  SS       UU    UU MMMM   MMMM       C
C       CC       DD    DD AA    AA  SSSSSS  UU    UU MM MM MM MM       C
C       CC       DD    DD AAAAAAAA       SS UU    UU MM  MMM  MM       C
C       CC    CC DD    DD AA    AA SS    SS UU    UU MM   M   MM       C
C        CCCCCC  DDDDDDD  AA    AA  SSSSSS   UUUUUU  MM       MM       C
C                                                                      C
C -------------------------------------------------------------------- C
C  CDASUM RETURNS THE SUM OF MAGNITUDES OF A COMPLEX VECTOR DX(N).     C
C**********************************************************************C
C
      COMPLEX*16 DX(N)
C
C     INITIALISE COUNTER
      CDASUM = 0.D0
C
C     IF VECTOR LENGTH IS ZERO, RESULT IS ZERO
      IF(N.LE.0) RETURN
C
C     CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
      M = MOD(N,6)
C
      IF(M.EQ.0) GOTO 40
C
      DO I=1,M
        CDASUM = CDASUM + ABS(DX(I))
      ENDDO
C
      IF(N.LT.6) RETURN
C
40    CONTINUE
C
      MP1 = M+1
      DO I=MP1,N,6
        CDASUM = CDASUM + ABS(DX(I  )) + ABS(DX(I+1)) + ABS(DX(I+2))
     &                  + ABS(DX(I+3)) + ABS(DX(I+4)) + ABS(DX(I+5))
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE NCART(IVECT,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             NN    NN  CCCCCC     AA    RRRRRRR TTTTTTTT              C
C             NNN   NN CC    CC   AAAA   RR    RR   TT                 C
C             NNNN  NN CC        AA  AA  RR    RR   TT                 C
C             NN NN NN CC       AA    AA RR    RR   TT                 C
C             NN  NNNN CC       AAAAAAAA RRRRRRR    TT                 C
C             NN   NNN CC    CC AA    AA RR    RR   TT                 C
C             NN    NN  CCCCCC  AA    AA RR    RR   TT                 C
C                                                                      C
C -------------------------------------------------------------------- C
C  NCART RETURNS THE CARTESIAN INDEX FROM THE INDEX VALUE IND.         C
C**********************************************************************C
C
      DIMENSION IVECT(3)
C
      IF(IND.EQ.1) THEN
        IVECT(1) = 1
        IVECT(2) = 0
        IVECT(3) = 0
      ELSEIF(IND.EQ.2) THEN
        IVECT(1) = 0
        IVECT(2) = 1
        IVECT(3) = 0
      ELSEIF(IND.EQ.3) THEN
        IVECT(1) = 0
        IVECT(2) = 0
        IVECT(3) = 1
      ELSE
        WRITE(6, *) 'In NCART: supplied index not valid',IND
        WRITE(7, *) 'In NCART: supplied index not valid',IND
        RETURN
      ENDIF
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C   [6] MULTI-CONFIG: MANY-CETRE MULTICONFIG. SCF CALCULATIONS.        C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] MCSCF: MAIN ROUTINE FOR MULTI-CONFIGURATIONAL SCF CALCULATION. C
C**********************************************************************C
C
C
      SUBROUTINE MCSCF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           MM       MM  CCCCCC   SSSSSS   CCCCCC  FFFFFFFF            C
C           MMM     MMM CC    CC SS    SS CC    CC FF                  C
C           MMMM   MMMM CC       SS       CC       FF                  C
C           MM MM MM MM CC        SSSSSS  CC       FFFFFF              C
C           MM  MMM  MM CC             SS CC       FF                  C
C           MM   M   MM CC    CC SS    SS CC    CC FF                  C
C           MM       MM  CCCCCC   SSSSSS   CCCCCC  FF                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  MCSCF PERFORMS A MULTI-CONFIGURATIONAL SCF CALCULATION USING THE    C
C  APPROACH OF KNOWLES AND WERNER (1985,1988).                         C
C**********************************************************************C
C
      CHARACTER*4  HMLTN
C
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
C
C     PRINT A BIG SECTION HEADER
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) REPEAT(' ',20),'MULTI-CONFIGURATIONAL MOLECULAR SCF'
      WRITE(7, *) REPEAT(' ',20),'MULTI-CONFIGURATIONAL MOLECULAR SCF'
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C     EQ-COEFFICIENTS IN COMMON ARRAYS
      IF(IEQS.EQ.1.AND.INEW.NE.0) THEN
        CALL EQFILE
      ENDIF
C
C     WARN USER THAT ROUTINE HASN'T BEEN WRITTEN YET
      WRITE(6, *) 'In BERTHA: MCSCF option not yet available.'
      WRITE(7, *) 'In BERTHA: MCSCF option not yet available.'
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C   [7] DMRG: DENSITY MATRIX RENORMALISATION GROUP CALCULATIONS.       C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] DMRG: DENSITY MATRIX RENORMALISATION GROUP CALCULATION.        C
C**********************************************************************C
C
C
      SUBROUTINE DMRG
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                DDDDDDD  MM       MM RRRRRRR   GGGGGG                 C
C                DD    DD MMM     MMM RR    RR GG    GG                C
C                DD    DD MMMM   MMMM RR    RR GG                      C
C                DD    DD MM MM MM MM RR    RR GG                      C
C                DD    DD MM  MMM  MM RRRRRRR  GG   GGG                C
C                DD    DD MM   M   MM RR    RR GG    GG                C
C                DDDDDDD  MM       MM RR    RR  GGGGGG                 C
C                                                                      C
C -------------------------------------------------------------------- C
C  DMRG PERFORMS A CALCULATION BASED ON THE FORMALISM OF THE DENSITY   C
C  MATRIX RENORMALISATION GROUP.                                       C
C**********************************************************************C
C
      CHARACTER*4  HMLTN
C
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
C     PRINT A BIG SECTION HEADER
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) REPEAT(' ',20),'DENSITY MATRIX RENORMALISATION GROUP'
      WRITE(7, *) REPEAT(' ',20),'DENSITY MATRIX RENORMALISATION GROUP'
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C     EQ-COEFFICIENTS IN COMMON ARRAYS
      IF(IEQS.EQ.1.AND.INEW.NE.0) THEN
        CALL EQFILE
      ENDIF
C
C     WARN USER THAT ROUTINE HASN'T BEEN WRITTEN YET
      WRITE(6, *) 'In BERTHA: DMRG option not yet available.'
      WRITE(7, *) 'In BERTHA: DMRG option not yet available.'
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C   [8] MBPT: CORRELATION ENERGY CALCULATION ROUTINES.                 C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] MBPT: MAIN ROUTINE FOR MANY-BODY DIAGRAMMATIC P.T.             C
C   [B] MP1: ZERO AND FIRST ORDER ENERGY ANALYSIS.                     C
C   [C] MP2: CORRELATION ENERGY BASED ON PRIOR CALCULATION.            C
C**********************************************************************C
C
C
      SUBROUTINE MBPT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                MM       MM BBBBBBB  PPPPPPP  TTTTTTTT                C
C                MMM     MMM BB    BB PP    PP    TT                   C
C                MMMM   MMMM BB    BB PP    PP    TT                   C
C                MM MM MM MM BBBBBBB  PP    PP    TT                   C
C                MM  MMM  MM BB    BB PPPPPPP     TT                   C
C                MM   M   MM BB    BB PP          TT                   C
C                MM       MM BBBBBBB  PP          TT                   C
C                                                                      C
C -------------------------------------------------------------------- C
C  MBPT PERFORMS MANY-BODY DIAGRAMMATIC EVALUATION ON A CONVERGED      C
C  MOLECULAR HARTREE-FOCK SOLUTION SPACE.                              C
C**********************************************************************C
C
      CHARACTER*4  HMLTN
C
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
C     PRINT A BIG SECTION HEADER
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) REPEAT(' ',23),'MANY-BODY PERTURBATION THEORY'
      WRITE(7, *) REPEAT(' ',23),'MANY-BODY PERTURBATION THEORY'
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C     RECORD TIME
      CALL CPU_TIME(TDUM)
C
C     EQ-COEFFICIENTS IN COMMON ARRAYS
      IF(IEQS.EQ.1.AND.INEW.NE.0) THEN
        CALL EQFILE
      ENDIF
C
C     CALL FIRST-ORDER MBPT ROUTINE
      CALL MP1
C
C     CALL SECOND-ORDER MBPT ROUTINE
      CALL MP2
C
C     TOTAL TIME TAKEN
      CALL CPU_TIME(TMPT)
      TMPT = TMPT-TDUM
C
      RETURN
      END
C
C
      SUBROUTINE MP1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                       MM       MM PPPPPPP   11                       C
C                       MMM     MMM PP    PP 111                       C
C                       MMMM   MMMM PP    PP  11                       C
C                       MM MM MM MM PP    PP  11                       C
C                       MM  MMM  MM PPPPPPP   11                       C
C                       MM   M   MM PP        11                       C
C                       MM       MM PP       1111                      C
C                                                                      C
C -------------------------------------------------------------------- C
C  MP1 EVALUATES ZERO AND FIRST ORDER ENERGIES FOR ALL OCCUPIED        C
C  SOLUTIONS TO A CONVERGED HARTREE-FOCK PROBLEM.                      C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000)
C
      CHARACTER*4  HMLTN
      CHARACTER*16 HMS
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 C(MDM,MDM),RR(MB2,16),ETMP1,ETMP2,ETMP3,ETMP4
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NFUNS(4),LQN(4)
      DIMENSION ITQN(2),IFLG(11),ISCF(11,6)
      DIMENSION INDEX(MCT,-MKP:MKP,2*(MKP+1)*MKP)
      DIMENSION T1(MDM),T2(MDM),T3(MDM)
C
      COMMON/COEF/C
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EMDR,EMXC,EUEH
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/QNMS/ICNBAS(MDM),KQNBAS(MDM),MQNBAS(MDM)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
      DIMENSION EAB1(NOCC,NOCC,6),EA1(NOCC,6)
C
C     ISCF TELLS WHICH INTEGRALS TO INCLUDE BASED ON OVERLAP COMBINATION
      DATA ISCF/1,0,1,0,1,0,1,0,1,1,0,
     &          1,0,0,0,1,0,0,0,1,1,0,
     &          0,1,1,0,0,0,0,0,1,1,0,
     &          1,0,0,1,1,0,0,0,1,0,0,
     &          0,1,0,1,0,0,0,0,1,0,0,
     &          0,1,0,0,0,0,0,0,1,0,0/
C
C
      CALL CPU_TIME(TDM1)
C
C     INITIALISE ENERGY COUNTERS
      DO N=1,6
        DO IOCCB=1,NOCC
          EA1(IOCCB,N) = 0.0D0
          DO IOCCA=1,NOCC
            EAB1(IOCCA,IOCCB,N) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     INDEXING ROUTINE: SO WE CAN SET BASIS FUNCTION LABELS BY BLOCK
      ICOUNT = 0
C
C     LOOP OVER NUCLEAR CENTERS
      DO ICNT=1,NCNT
C
C       LOOP OVER KAPPA VALUES FOR THIS NUCLEAR CENTER
        DO KN=1,NKAP(ICNT)
C
C         IMPORT KAPPA, MAXIMUM MQN
          KAPPA = KVALS(KN,ICNT)
          MJMAX  = 2*IABS(KAPPA)-1
C
C         LOOP OVER MQN VALUES AND RECORD INDEX
          DO MJ=1,MJMAX,2
            ICOUNT = ICOUNT+1
            INDEX(ICNT,KAPPA,MJ) = ICOUNT
          ENDDO
        ENDDO
      ENDDO
C
C     SET THE LIMITS FOR RUNNING OVER DENSITY COMBINATIONS
      IF(HMLTN.EQ.'NORL') THEN
        ITSTRT = 1
        ITSTOP = 1
        ITSKIP = 1
      ELSE
        ITSTRT = 1
        ITSTOP = 4
        ITSKIP = 3
      ENDIF
C
C     PRINT ENERGY CONTRIBUTIONS TO EACH ELECTRON PAIR
20    FORMAT(1X,A,12X,A,11X,A,11X,A,12X,A)
21    FORMAT(' (',I2,',',I2,')',4X,F13.7,5X,F11.7,5X,F11.7,5X,F13.7)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',24),'MP1 pair-wise summary'
      WRITE(7, *) REPEAT(' ',24),'MP1 pair-wise summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) '( a, b)','E1(H)','E1(J)','E1(K)','E1(ab)'
      WRITE(7,20) '( a, b)','E1(H)','E1(J)','E1(K)','E1(ab)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     CALCULATE ONE-BODY MATRIX REPS
      CALL ONEEL
C
C**********************************************************************C
C     ONE-BODY ENERGIES                                                C
C**********************************************************************C
C
      DO 1000 IOCCA=1,NOCC
      DO 1000 IOCCB=1,IOCCA
C
        IA = IOCCA + NSHIFT
        IB = IOCCB + NSHIFT
C
C       ONE-BODY ENERGY
        ETMP1 = DCMPLX(0.0D0,0.0D0)
        ETMP2 = DCMPLX(0.0D0,0.0D0)
        IF(IOCCA.NE.IOCCB) GOTO 90
        DO J=1,NDIM
          DO I=1,NDIM
            ETMP1 = ETMP1 + DCONJG(C(I,IA))*C(J,IA)*HNUC(I,J)
            ETMP2 = ETMP2 + DCONJG(C(I,IA))*C(J,IA)*HKIN(I,J)
          ENDDO
        ENDDO
90      CONTINUE
        EAB1(IOCCA,IOCCB,1) = DREAL(ETMP1)
        EAB1(IOCCA,IOCCB,2) = DREAL(ETMP2)
        EAB1(IOCCA,IOCCB,3) = DREAL(ETMP1 + ETMP2)
C
C       INITIALISE COULOMB DIRECT AND EXCHANGE ARRAYS
        DO J=1,NDIM
          DO I=1,NDIM
            GDIR(I,J) = DCMPLX(0.0D0,0.0D0)
            GXCH(I,J) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
cC
c100   CONTINUE
C
C**********************************************************************C
C     LOOP OVER UNIQUE OCCUPIED ORBITALS (A,B)       (INDEX 1000)      C
C**********************************************************************C
C
c      DO 1000 IOCCA=1,NOCC
c      DO 1000 IOCCB=1,IOCCA
cC
c        IA = IOCCA + NSHIFT
c        IB = IOCCB + NSHIFT
C
C**********************************************************************C
C     LOOP OVER SPINORS A, B BY BLOCK                (INDEX 3000)      C
C**********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR A AND B: T = {L} OR {L,S}
      DO 3000 IT1=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITQN(1) = IT1
C
C       CALCULATE STARTING ADDRESS
        IF(IT1.EQ.1) THEN
          NADDAB = 0
        ELSE
          NADDAB = NSHIFT
        ENDIF
C
C     LOOP OVER CENTER A
      DO 3000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 3000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
      DO 3000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNS(1) = NFUNCT(LQN(1)+1,ICNTA)
        DO IBAS=1,NFUNS(1)
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 3000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
         LQN(2) = KQN(2)
        ELSE
         LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNS(2) = NFUNCT(LQN(2)+1,ICNTB)
        DO JBAS=1,NFUNS(2)
          EXPT(JBAS,2) = EXPSET(JBAS,LQN(2)+1,ICNTB)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 3000 MA=1,IABS(KQN(1))
        MJA    = 2*MA-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 3000 MB=1,IABS(KQN(2))
        MJB    = 2*MB-1
        MQN(2) = MJB
C
C     CALCULATE NEW BLOCK OF E(AB) COEFFS AT NEXT OPPORTUNITY
      IEAB = 1
      IABLL = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB)
      IABSS = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB)
C
C**********************************************************************C
C     FOR THIS CHOICE OF A AND B, COMPUTE ADDRESSES AND PHASES         C
C**********************************************************************C
C
C     CALCULATE BLOCK INDICES FOR {AB} COMBINATIONS
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
C
      IQ12 = (IQ1*(IQ1-1))/2 + IQ2
C
C     FURTHER DEFINE STARTING ADDRESSES FOR {ABCD} BASIS LABELS
      IA1 = LARGE(ICNTA,KA,MJA  ) + NADDAB
      IB1 = LARGE(ICNTB,KB,MJB  ) + NADDAB
C
      IA2 = LARGE(ICNTA,KA,MJA+1) + NADDAB
      IB2 = LARGE(ICNTB,KB,MJB+1) + NADDAB
C
      JA1 = LARGE(ICNTA,KA,MJA  ) + NADDAB
      JB1 = LARGE(ICNTB,KB,MJB  ) + NADDAB
C
      JA2 = LARGE(ICNTA,KA,MJA+1) + NADDAB
      JB2 = LARGE(ICNTB,KB,MJB+1) + NADDAB
C
C     CALCULATE KQN PHASE FACTORS FOR PERMUTING INTEGRALS
      IF((KQN(1)*KQN(2)).GT.0) THEN 
        PKAB = 1.0D0
      ELSE
        PKAB =-1.0D0
      ENDIF
C
C     CALCULATE MQN PHASE FACTORS FOR PERMUTING INTEGRALS
      PMAB1 = DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PMAB2 = DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
C
      F1 = PKAB*PMAB1
      F2 = PKAB*PMAB2
C
C**********************************************************************C
C     LOOP OVER SPINORS C, D BY BLOCK                (INDEX 4000)      C
C**********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR C AND D: T = {L} OR {L,S}
      DO 4000 IT2=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITQN(2) = IT2
C
C       CALCULATE STARTING ADDRESS
        IF(IT2.EQ.1) THEN
          NADDCD = 0
        ELSE
          NADDCD = NSHIFT
        ENDIF
C
C     LOOP OVER CENTER C
      DO 4000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER C
        XYZ(1,3) = COORD(1,ICNTC)
        XYZ(2,3) = COORD(2,ICNTC)
        XYZ(3,3) = COORD(3,ICNTC)
C
C     LOOP OVER CENTER D
      DO 4000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER D
        XYZ(1,4) = COORD(1,ICNTD)
        XYZ(2,4) = COORD(2,ICNTD)
        XYZ(3,4) = COORD(3,ICNTD)
C
C     LOOP OVER KQN(C) VALUES
      DO 4000 KC=1,NKAP(ICNTC)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR C
        KQN(3) = KVALS(KC,ICNTC)
        IF(KQN(3).GT.0) THEN
          LQN(3) = KQN(3)
        ELSE
          LQN(3) =-KQN(3)-1
        ENDIF
C         
        NFUNS(3) = NFUNCT(LQN(3)+1,ICNTC)
        DO KBAS=1,NFUNS(3)
          EXPT(KBAS,3) = EXPSET(KBAS,LQN(3)+1,ICNTC)
        ENDDO
C
C     LOOP OVER KQN(D) VALUES
      DO 4000 KD=1,NKAP(ICNTD)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR D
        KQN(4) = KVALS(KD,ICNTD)
        IF(KQN(4).GT.0) THEN
          LQN(4) = KQN(4)
        ELSE
          LQN(4) =-KQN(4)-1
        ENDIF
C
        NFUNS(4) = NFUNCT(LQN(4)+1,ICNTD)
        DO LBAS=1,NFUNS(4)
          EXPT(LBAS,4) = EXPSET(LBAS,LQN(4)+1,ICNTD)
        ENDDO
C
C     LOOP OVER |MQN(C)| VALUES
      DO 4000 MC=1,IABS(KQN(3))
        MJC    = 2*MC-1
        MQN(3) = MJC
C
C     LOOP OVER |MQN(D)| VALUES
      DO 4000 MD=1,IABS(KQN(4))
        MJD    = 2*MD-1
        MQN(4) = MJD
C
C     CALCULATE NEW BLOCK OF E(CD) COEFFS AT NEXT OPPORTUNITY
      IECD = 1
      ICDLL = IAD0LL(ICNTC,ICNTD,KC,KD,MC,MD)
      ICDSS = IAD0SS(ICNTC,ICNTD,KC,KD,MC,MD)
C
C**********************************************************************C
C     FOR THIS CHOICE OF C AND D, COMPUTE ADDRESSES AND PHASES         C
C**********************************************************************C
C
C     CALCULATE BLOCK INDICES FOR {CD} COMBINATIONS
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
      IQ34 = (IQ3*(IQ3-1))/2 + IQ4
C
C     FURTHER DEFINE STARTING ADDRESSES FOR {CD} BASIS LABELS
      IC1 = LARGE(ICNTC,KC,MJC  ) + NADDCD
      ID1 = LARGE(ICNTD,KD,MJD  ) + NADDCD
C
      IC2 = LARGE(ICNTC,KC,MJC+1) + NADDCD
      ID2 = LARGE(ICNTD,KD,MJD+1) + NADDCD
C
      JC1 = LARGE(ICNTC,KC,MJC  ) + NADDCD
      JD1 = LARGE(ICNTD,KD,MJD  ) + NADDCD
C
      JC2 = LARGE(ICNTC,KC,MJC+1) + NADDCD
      JD2 = LARGE(ICNTD,KD,MJD+1) + NADDCD
C        
      IF((KQN(3)*KQN(4)).GT.0) THEN 
        PKCD = 1.0D0
      ELSE
        PKCD =-1.0D0
      ENDIF
C
C     CALCULATE MQN PHASE FACTORS FOR PERMUTING INTEGRALS
      PMCD1 = DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PMCD2 = DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
      G1 = PKCD*PMCD1
      G2 = PKCD*PMCD2
C
C**********************************************************************C
C     SKIP BATCHES WHICH CONFORM TO INTEGRAL SYMMETRIES                C
C**********************************************************************C
C
C     DIATOMIC MOLECULES CARRY STRICT SELECTION RULES ON MQNS
      IF(NCNT.LE.2) THEN
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) GOTO 3999
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) GOTO 3999
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) GOTO 3999
        GOTO 4000
      ENDIF
3999  CONTINUE
C
C     DECISION TREE FOR SKIPPING CONTRIBUTIONS DUE TO INTEGRAL SYMMETRY
      IF(IQ1.LT.IQ2)   GOTO 4000
      IF(IQ3.LT.IQ4)   GOTO 4000
      IF(IQ12.LT.IQ34) GOTO 4000
C
C     INDICATE BLOCKS TO BE INCLUDED AHEAD GIVEN A,B,C,D BASIS QNMS...
C     A =/= B AND C =/= D WITH AB LIST VALUE =/= CD LIST VALUE
      IF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 1
C     A=/=B AND C=/=D WITH IND(AB) = IND(CD)
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.EQ.IQ34) THEN
        ITSCF = 2
C     A=/=B AND C = D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 3
C     A = B AND C=/=D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 4
C     A = B AND C = D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 5
C     A = B AND C = D WITH IND(AB) = IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.EQ.IQ34) THEN
        ITSCF = 6
C     COMBINATION OF A,B,C,D NOT TO BE INCLUDED -- USE MATRIX CONJ LATER
      ELSE
        GO TO 4000
      ENDIF
C
C     READ IN FLAG VALUES FROM ISCF DATA BLOCK
      DO N=1,11
        IFLG(N) = ISCF(N,ITSCF)
      ENDDO
C
C     INCLUDE SPECIAL CASES...
C     A=/=B AND C=/=D WITH IND(AB)=/=IND(CD)
      IF(ITSCF.EQ.1) THEN
C       A = C
        IF(IQ1.EQ.IQ3) IFLG( 6) = 1
C       B = D
        IF(IQ2.EQ.IQ4) IFLG(11) = 1
C       B = C
        IF(IQ2.EQ.IQ3) IFLG( 8) = 1
C     A=/=B AND C = D WITH IND(AB)=/=IND(CD)
      ELSEIF(ITSCF.EQ.3) THEN
C       B = C
        IF(IQ2.EQ.IQ3) IFLG( 8) = 1
C     A = B AND C=/=D WITH IND(AB)=/=IND(CD)
      ELSEIF(ITSCF.EQ.4) THEN
C       B = C
        IF(IQ2.EQ.IQ3) IFLG( 8) = 1
      ENDIF
C
C
C**********************************************************************C
C     FOURTH LAYER OF LOOPS, OVER BASIS FUNCTIONS A AND B (5000)       C
C -------------------------------------------------------------------- C
C     THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH         C
C     GENERATE THE GDIR AND GXCH MATRICES FROM THE SPINOR INTEGRALS.   C
C     THESE INCLUDE PHASE FACTORS FOR THE PERMUTATION OF               C
C     KQN(1) <-> KQN(2) AND MQN(1) <-> MQN(2)                          C
C -------------------------------------------------------------------- C
C     (RSCF 86, 87)                                                    C
C**********************************************************************C
C
C     LOOP OVER ELEMENTS OF FOCK MATRIX BLOCK
      DO 5000 IBAS=1,NFUNS(1)
      DO 5000 JBAS=1,NFUNS(2)
C
C       GENERATE BATCH OF ELECTRON REPULSION INTEGRALS
        CALL ERI(RR,XYZ,KQN,MQN,EXPT,NFUNS,ITQN,IBAS,JBAS,IEAB,IECD)
C
C       THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH
C       GENERATE THE GDIR/GXCH MATRIX FROM THE SPINOR INTEGRALS.
C
C       FIRST IFLG BATCH (DIRECT)
        IF(IFLG(1).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GDIR(IA1+IBAS,JB1+JBAS) = GDIR(IA1+IBAS,JB1+JBAS)
     &           +       RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
     &           +    G1*RR(M, 4)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M, 2)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IB)
     &           +    G2*RR(M, 3)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M, 1)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IB)
C
              GDIR(IA1+IBAS,JB2+JBAS) = GDIR(IA1+IBAS,JB2+JBAS)
     &           +       RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
     &           +    G1*RR(M, 8)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M, 6)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IB)
     &           +    G2*RR(M, 7)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M, 5)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IB)
C
              GDIR(IA2+IBAS,JB1+JBAS) = GDIR(IA2+IBAS,JB1+JBAS)
     &           +       RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
     &           +    G1*RR(M,12)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M,10)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IB)
     &           +    G2*RR(M,11)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M, 9)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IB)
C
              GDIR(IA2+IBAS,JB2+JBAS) = GDIR(IA2+IBAS,JB2+JBAS)
     &           +       RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
     &           +    G1*RR(M,16)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M,14)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IB)
     &           +    G2*RR(M,15)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M,13)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IB)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SECOND IFLG BATCH (DIRECT)
        IF(IFLG(2).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GDIR(IA1+IBAS,JB1+JBAS) = GDIR(IA1+IBAS,JB1+JBAS)
     &           +       RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
C
              GDIR(IA1+IBAS,JB2+JBAS) = GDIR(IA1+IBAS,JB2+JBAS)
     &           +       RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
C
              GDIR(IA2+IBAS,JB1+JBAS) = GDIR(IA2+IBAS,JB1+JBAS)
     &           +       RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
C
              GDIR(IA2+IBAS,JB2+JBAS) = GDIR(IA2+IBAS,JB2+JBAS)
     &           +       RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IB)
     &           +       RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IB)
     &           +       RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IB)
C
            ENDDO
          ENDDO
        ENDIF

C       THIRD IFLG BATCH (DIRECT)
        IF(IFLG(3).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GDIR(IC1+KBAS,JD1+LBAS) = GDIR(IC1+KBAS,JD1+LBAS)
     &           +       RR(M, 1)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 5)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M, 9)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,13)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
     &           +    F1*RR(M,13)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F2*RR(M, 5)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IB)
     &           +    F2*RR(M, 9)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F1*RR(M, 1)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IB)
C
              GDIR(IC1+KBAS,JD2+LBAS) = GDIR(IC1+KBAS,JD2+LBAS)
     &           +       RR(M, 2)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 6)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,10)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,14)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
     &           +    F1*RR(M,14)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F2*RR(M, 6)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IB)
     &           +    F2*RR(M,10)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F1*RR(M, 2)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IB)
C
              GDIR(IC2+KBAS,JD1+LBAS) = GDIR(IC2+KBAS,JD1+LBAS)
     &           +       RR(M, 3)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 7)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,11)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,15)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
     &           +    F1*RR(M,15)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F2*RR(M, 7)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IB)
     &           +    F2*RR(M,11)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F1*RR(M, 3)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IB)
C
              GDIR(IC2+KBAS,JD2+LBAS) = GDIR(IC2+KBAS,JD2+LBAS)
     &           +       RR(M, 4)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 8)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,12)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,16)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
     &           +    F1*RR(M,16)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F2*RR(M, 8)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IB)
     &           +    F2*RR(M,12)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IB)
     &           +    F1*RR(M, 4)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IB)
C 
            ENDDO
          ENDDO
        ENDIF
C
C       FOURTH IFLG BATCH (DIRECT)
        IF(IFLG(4).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GDIR(IC1+KBAS,JD1+LBAS) = GDIR(IC1+KBAS,JD1+LBAS)
     &           +       RR(M, 1)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 5)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M, 9)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,13)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
C
              GDIR(IC1+KBAS,JD2+LBAS) = GDIR(IC1+KBAS,JD2+LBAS)
     &           +       RR(M, 2)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 6)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,10)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,14)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
C
              GDIR(IC2+KBAS,JD1+LBAS) = GDIR(IC2+KBAS,JD1+LBAS)
     &           +       RR(M, 3)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 7)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,11)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,15)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
C
              GDIR(IC2+KBAS,JD2+LBAS) = GDIR(IC2+KBAS,JD2+LBAS)
     &           +       RR(M, 4)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 8)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,12)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,16)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IB)
C
            ENDDO
          ENDDO
        ENDIF
C
C       FIFTH IFLG BATCH (EXCHANGE)
        IF(IFLG(5).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GXCH(IA1+IBAS,JC1+KBAS) = GXCH(IA1+IBAS,JC1+KBAS)
     &           +    G1*RR(M, 4)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G1*RR(M, 8)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IB)
     &           +    G2*RR(M, 3)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G2*RR(M, 7)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IB)
C
              GXCH(IA1+IBAS,JC2+KBAS) = GXCH(IA1+IBAS,JC2+KBAS)
     &           +    G2*RR(M, 2)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G2*RR(M, 6)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IB)
     &           +    G1*RR(M, 1)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G1*RR(M, 5)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IB)
C
              GXCH(IA2+IBAS,JC1+KBAS) = GXCH(IA2+IBAS,JC1+KBAS)
     &           +    G1*RR(M,12)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G1*RR(M,16)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IB)
     &           +    G2*RR(M,11)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G2*RR(M,15)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IB)
C
              GXCH(IA2+IBAS,JC2+KBAS) = GXCH(IA2+IBAS,JC2+KBAS)
     &           +    G2*RR(M,10)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G2*RR(M,14)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IB)
     &           +    G1*RR(M, 9)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IB)
     &           +    G1*RR(M,13)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IB)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SIXTH IFLG BATCH (EXCHANGE)
        IF(IFLG(6).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GXCH(IC1+KBAS,JA1+IBAS) = GXCH(IC1+KBAS,JA1+IBAS)
     &           +    F1*RR(M,13)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F1*RR(M,14)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IB)
     &           +    F2*RR(M, 9)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F2*RR(M,10)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IB)
C
              GXCH(IC1+KBAS,JA2+IBAS) = GXCH(IC1+KBAS,JA2+IBAS)
     &           +    F2*RR(M, 5)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F2*RR(M, 6)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IB)
     &           +    F1*RR(M, 1)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F1*RR(M, 2)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IB)
C
              GXCH(IC2+KBAS,JA1+IBAS) = GXCH(IC2+KBAS,JA1+IBAS)
     &           +    F1*RR(M,15)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F1*RR(M,16)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IB)
     &           +    F2*RR(M,11)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F2*RR(M,12)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IB)
C
              GXCH(IC2+KBAS,JA2+IBAS) = GXCH(IC2+KBAS,JA2+IBAS)
     &           +    F2*RR(M, 7)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F2*RR(M, 8)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IB)
     &           +    F1*RR(M, 3)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IB)
     &           +    F1*RR(M, 4)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IB)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SEVENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(7).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GXCH(IB1+JBAS,JC1+KBAS) = GXCH(IB1+JBAS,JC1+KBAS)
     &           + F1*G1*RR(M,16)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F2*G1*RR(M, 8)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IB)
     &           + F1*G2*RR(M,15)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F2*G2*RR(M, 7)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IB)
C
              GXCH(IB1+JBAS,JC2+KBAS) = GXCH(IB1+JBAS,JC2+KBAS)
     &           + F1*G2*RR(M,14)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F2*G2*RR(M, 6)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IB)
     &           + F1*G1*RR(M,13)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F2*G1*RR(M, 5)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IB)
C
              GXCH(IB2+JBAS,JC1+KBAS) = GXCH(IB2+JBAS,JC1+KBAS)
     &           + F2*G1*RR(M,12)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F1*G1*RR(M, 4)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IB)
     &           + F2*G2*RR(M,11)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F1*G2*RR(M, 3)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IB)
C
              GXCH(IB2+JBAS,JC2+KBAS) = GXCH(IB2+JBAS,JC2+KBAS)
     &           + F2*G2*RR(M,10)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F1*G2*RR(M, 2)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IB)
     &           + F2*G1*RR(M, 9)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IB)
     &           + F1*G1*RR(M, 1)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IB)
C
            ENDDO
          ENDDO
        ENDIF
C
C       EIGHTH IFLG BATCH (EXCHANGE)
        IF(IFLG(8).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNS(3)
            DO LBAS=1,NFUNS(4)
              M = M+1
C
              GXCH(IC1+KBAS,JB1+JBAS) = GXCH(IC1+KBAS,JB1+JBAS)
     &           +       RR(M, 1)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M, 2)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IB)
     &           +       RR(M, 9)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M,10)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IB)
C
              GXCH(IC1+KBAS,JB2+JBAS) = GXCH(IC1+KBAS,JB2+JBAS)
     &           +       RR(M, 5)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M, 6)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IB)
     &           +       RR(M,13)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M,14)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IB)
C
              GXCH(IC2+KBAS,JB1+JBAS) = GXCH(IC2+KBAS,JB1+JBAS)
     &           +       RR(M, 3)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M, 4)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IB)
     &           +       RR(M,11)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M,12)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IB)
C
              GXCH(IC2+KBAS,JB2+JBAS) = GXCH(IC2+KBAS,JB2+JBAS)
     &           +       RR(M, 7)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M, 8)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IB)
     &           +       RR(M,15)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IB)
     &           +       RR(M,16)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IB)
C
            ENDDO
          ENDDO
        ENDIF
C
C       NINTH IFLG BATCH (EXCHANGE)
        IF(IFLG(9).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C  
              GXCH(IA1+IBAS,JD1+LBAS) = GXCH(IA1+IBAS,JD1+LBAS)
     &           +       RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IB)
C
              GXCH(IA1+IBAS,JD2+LBAS) = GXCH(IA1+IBAS,JD2+LBAS)
     &           +       RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IB)
C
              GXCH(IA2+IBAS,JD1+LBAS) = GXCH(IA2+IBAS,JD1+LBAS)
     &           +       RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IB)
C
              GXCH(IA2+IBAS,JD2+LBAS) = GXCH(IA2+IBAS,JD2+LBAS)
     &           +       RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IB)
     &           +       RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IB)
     &           +       RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IB)
C
            ENDDO
          ENDDO
        ENDIF
C
C       TENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(10).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C
              GXCH(IB1+JBAS,JD1+LBAS) = GXCH(IB1+JBAS,JD1+LBAS)
     &           +    F1*RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F2*RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IB)
     &           +    F1*RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F2*RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IB)
C
              GXCH(IB1+JBAS,JD2+LBAS) = GXCH(IB1+JBAS,JD2+LBAS)
     &           +    F1*RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F2*RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IB)
     &           +    F1*RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F2*RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IB)
C
              GXCH(IB2+JBAS,JD1+LBAS) = GXCH(IB2+JBAS,JD1+LBAS)
     &           +    F2*RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F1*RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IB)
     &           +    F2*RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F1*RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IB)
C
              GXCH(IB2+JBAS,JD2+LBAS) = GXCH(IB2+JBAS,JD2+LBAS)
     &           +    F2*RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F1*RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IB)
     &           +    F2*RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IB)
     &           +    F1*RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IB)
C
            ENDDO
          ENDDO
        ENDIF
C
C       ELEVENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(11).EQ.1) THEN
          DO LBAS=1,NFUNS(4)
            DO KBAS=1,NFUNS(3)
              M = (KBAS-1)*NFUNS(4)+LBAS
C
              GXCH(ID1+LBAS,JB1+JBAS) = GXCH(ID1+LBAS,JB1+JBAS)
     &           +    G1*RR(M, 4)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M, 2)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IB)
     &           +    G1*RR(M,12)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M,10)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IB)
C
              GXCH(ID1+LBAS,JB2+JBAS) = GXCH(ID1+LBAS,JB2+JBAS)
     &           +    G1*RR(M, 8)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M, 6)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IB)
     &           +    G1*RR(M,16)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G2*RR(M,14)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IB)
C
              GXCH(ID2+LBAS,JB1+JBAS) = GXCH(ID2+LBAS,JB1+JBAS)
     &           +    G2*RR(M, 3)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M, 1)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IB)
     &           +    G2*RR(M,11)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M, 9)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IB)
C
              GXCH(ID2+LBAS,JB2+JBAS) = GXCH(ID2+LBAS,JB2+JBAS)
     &           +    G2*RR(M, 7)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M, 5)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IB)
     &           +    G2*RR(M,15)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IB)
     &           +    G1*RR(M,13)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IB)
C
            ENDDO
          ENDDO
        ENDIF
C
C     END LOOP OVER IAOCC,IBOCC BLOCK ADDRESSES
5000  CONTINUE
C     END LOOPS OVER A,B,C,D OVERLAP BLOCKS
4000  CONTINUE
3000  CONTINUE
C
C     COMPLETE CONSTRUCTION OF GDIR AND GXCH BY MATRIX CONJUGATION
      DO J=1,NDIM-NSHIFT
        DO I=1,J
C
C         SMALL-COMPONENT ADDRESSES
          K = I + NSHIFT
          L = J + NSHIFT
C
C         SKIP DIAGONAL PARTS OF EACH SUB-BLOCK
          IF(ICNBAS(I).NE.ICNBAS(J)) GOTO 400
          IF(KQNBAS(I).NE.KQNBAS(J)) GOTO 400
          IF(IABS(MQNBAS(I)).NE.IABS(MQNBAS(J))) GOTO 400
          GOTO 401
400       CONTINUE
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF LL BLOCK
          GDIR(I,J) = GDIR(I,J) + DCONJG(GDIR(J,I))
          GDIR(J,I) =             DCONJG(GDIR(I,J))
          GXCH(I,J) = GXCH(I,J) + DCONJG(GXCH(J,I))
          GXCH(J,I) =             DCONJG(GXCH(I,J))
C
C         IF HMLTN = 'NORL' SKIP THE NEXT FEW CALCULATIONS
          IF(HMLTN.EQ.'NORL') GOTO 401
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF LS BLOCK
          GDIR(I,L) = GDIR(I,L) + DCONJG(GDIR(L,I))
          GDIR(L,I) =             DCONJG(GDIR(I,L))
          GXCH(I,L) = GXCH(I,L) + DCONJG(GXCH(L,I))
          GXCH(L,I) =             DCONJG(GXCH(I,L))
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF SL BLOCK
          GDIR(K,J) = GDIR(K,J) + DCONJG(GDIR(J,K))
          GDIR(J,K) =             DCONJG(GDIR(K,J))
          GXCH(K,J) = GXCH(K,J) + DCONJG(GXCH(J,K))
          GXCH(J,K) =             DCONJG(GXCH(K,J))
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF SS BLOCK
          GDIR(K,L) = GDIR(K,L) + DCONJG(GDIR(L,K))
          GDIR(L,K) =             DCONJG(GDIR(K,L))
          GXCH(K,L) = GXCH(K,L) + DCONJG(GXCH(L,K))
          GXCH(L,K) =             DCONJG(GXCH(K,L))
C
401       CONTINUE
        ENDDO
      ENDDO
C
C     END LOOP OVER UNIQUE VIRTUAL ORBITAL COMBINATIONS (R,S)
C2000  CONTINUE
C
C**********************************************************************C
C     CALCULATE INTERACTION ENERGY FROM ORBITAL PAIR (A,B)             C
C**********************************************************************C
C
      ETMP3 = DCMPLX(0.0D0,0.0D0)
      ETMP4 = DCMPLX(0.0D0,0.0D0)
      DO J=1,NDIM
        DO I=1,NDIM
          ETMP3 = ETMP3 + DCONJG(C(I,IA))*C(J,IA)*GDIR(I,J)
          ETMP4 = ETMP4 + DCONJG(C(I,IA))*C(J,IA)*GXCH(I,J)
        ENDDO
      ENDDO
      EAB1(IOCCA,IOCCB,4) = DREAL(ETMP3)
      EAB1(IOCCA,IOCCB,5) = DREAL(ETMP4)
      EAB1(IOCCA,IOCCB,6) = DREAL(ETMP1 + ETMP2 + ETMP3 - ETMP4)
C
C     OUTPUT ENERGIES TO TERMINAL
      WRITE(6,21) IOCCA,IOCCB,(EAB1(IOCCA,IOCCB,N),N=3,6)
      WRITE(7,21) IOCCA,IOCCB,(EAB1(IOCCA,IOCCB,N),N=3,6)
C
      IF(IOCCA.NE.IOCCB) THEN
        DO N=1,6
          EAB1(IOCCB,IOCCA,N) = EAB1(IOCCA,IOCCB,N)
        ENDDO
      ENDIF
C
C     END LOOP OVER UNIQUE OCCUPIED ORBITAL COMBINATIONS (A,B)
1000  CONTINUE
C
C**********************************************************************C
C     END OF CALCULATION OVER OCCUPIED ORBITALS (A,B)                  C
C**********************************************************************C
C
C     RECORD TIME AT THE END OF MAIN CALCULATION
      CALL CPU_TIME(TDM2)
C
C     WRITE RESULTS OF EAB ENERGIES TO AN EXTERNAL FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_MP1.dat',STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO IOCCA=1,NOCC
        DO IOCCB=1,NOCC
          WRITE(8, *) (EAB1(IOCCA,IOCCB,N),N=1,6)
        ENDDO
      ENDDO
      CLOSE(UNIT=8)
C
C     CALCULATE MP1 SINGLE-PARTICLE ENERGIES AND MOLECULAR TOTALS
      EENC1 = 0.0D0
      EKIN1 = 0.0D0
      EBAR1 = 0.0D0
      EDIR1 = 0.0D0
      EXCH1 = 0.0D0
      DO IOCCA=1,NOCC
        DO N=1,5
          EA1(IOCCA,N) = 0.0D0
          DO IOCCB=1,NOCC
            EA1(IOCCA,N) = EA1(IOCCA,N) + EAB1(IOCCA,IOCCB,N)
          ENDDO
        ENDDO
        EA1(IOCCA,6) = EA1(IOCCA,3) + EA1(IOCCA,4) - EA1(IOCCA,5)
        EENC1 = EENC1 + EA1(IOCCA,1)
        EKIN1 = EKIN1 + EA1(IOCCA,2)
        EBAR1 = EBAR1 + EA1(IOCCA,3)
        EDIR1 = EDIR1 + EA1(IOCCA,4)*0.5D0
        EXCH1 = EXCH1 + EA1(IOCCA,5)*0.5D0
      ENDDO
      ETOT1 = EBAR1 + EDIR1 - EXCH1 + ENUC
C
C     MP1 SINGLE-PARTICLE SUMMARY
30    FORMAT(1X,A,12X,A,11X,A,11X,A,12X,A)
31    FORMAT('  ',I2,'    ',4X,F13.7,5X,F11.7,5X,F11.7,5X,F13.7)

      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',21),'MP1 single particle summary'
      WRITE(7, *) REPEAT(' ',21),'MP1 single particle summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,30) '  a    ','E1(H)','E1(J)','E1(K)',' E1(a)'
      WRITE(7,30) '  a    ','E1(H)','E1(J)','E1(K)',' E1(a)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)  
      DO IOCCA=1,NOCC
        WRITE(6,31) IOCCA,(EA1(IOCCA,N),N=3,6)
        WRITE(7,31) IOCCA,(EA1(IOCCA,N),N=3,6)
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     SUMMARY OF MOLECULAR ENERGY SOURCES
32    FORMAT(1X,A,4X,A,2X,'=',10X,F18.8,' au')
33    FORMAT(1X,A,5X,'=',23X,A)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',25),'MP1 molecular summary'
      WRITE(7, *) REPEAT(' ',25),'MP1 molecular summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,32) 'Hartree-Fock electron-nucleus','E1(V)',EENC1
      WRITE(7,32) 'Hartree-Fock electron-nucleus','E1(V)',EENC1
      WRITE(6,32) 'Hartree-Fock electron kinetic','E1(T)',EKIN1
      WRITE(7,32) 'Hartree-Fock electron kinetic','E1(T)',EKIN1
      WRITE(6,32) 'Hartree-Fock Coulomb direct  ','E1(J)',EDIR1
      WRITE(7,32) 'Hartree-Fock Coulomb direct  ','E1(J)',EDIR1
      WRITE(6,32) 'Hartree-Fock Coulomb exchange','E1(K)',EXCH1
      WRITE(7,32) 'Hartree-Fock Coulomb exchange','E1(K)',EXCH1
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,32) 'Nuclear repulsion            ','E0(N)',ENUC
      WRITE(7,32) 'Nuclear repulsion            ','E0(N)',ENUC
      WRITE(6,32) 'Hartree-Fock one-electron    ','E1(H)',EBAR1
      WRITE(7,32) 'Hartree-Fock one-electron    ','E1(H)',EBAR1
      WRITE(6,32) 'Hartree-Fock Coulomb total   ','E1(G)',EDIR1-EXCH1
      WRITE(7,32) 'Hartree-Fock Coulomb total   ','E1(G)',EDIR1-EXCH1
      WRITE(6,32) 'Hartree-Fock molecular energy','E1   ',ETOT1
      WRITE(7,32) 'Hartree-Fock molecular energy','E1   ',ETOT1
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,33) 'MP1 time                   ',HMS(TDM2-TDM1)
      WRITE(7,33) 'MP1 time                   ',HMS(TDM2-TDM1)
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
      RETURN
      END
C
C
      SUBROUTINE MP2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                    MM       MM PPPPPPP   222222                      C
C                    MMM     MMM PP    PP 22    22                     C
C                    MMMM   MMMM PP    PP       22                     C
C                    MM MM MM MM PP    PP     22                       C
C                    MM  MMM  MM PPPPPPP    22                         C
C                    MM   M   MM PP       22                           C
C                    MM       MM PP       22222222                     C
C                                                                      C
C -------------------------------------------------------------------- C
C  MP2 EVALUATES SECOND ORDER COULOMB DIRECT/EXCHANGE ENERGIES         C
C  FOR ALL OCCUPIED SOLUTIONS TO A CONVERGED HARTREE-FOCK PROBLEM.     C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000)
C
      CHARACTER*4  HMLTN
      CHARACTER*16 HMS
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 C(MDM,MDM),RR(MB2,16)
      COMPLEX*16 ETMP1,ETMP2,ETMP3,ETMP4,RNUMV,RNUMT,RNUMJ,RNUMK
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NFUNS(4),LQN(4)
      DIMENSION ITQN(2),IFLG(11),ISCF(11,6)
      DIMENSION INDEX(MCT,-MKP:MKP,2*(MKP+1)*MKP)
      DIMENSION T1(MDM),T2(MDM),T3(MDM)
C
      COMMON/COEF/C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EMDR,EMXC,EUEH
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/MTRX/OVAP,HNUC,HKIN,VUEH,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/QNMS/ICNBAS(MDM),KQNBAS(MDM),MQNBAS(MDM)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
      COMPLEX*16 ABRS(NVIR,NVIR),BARS(NVIR,NVIR),AVR(NVIR),ATR(NVIR)
C
      DIMENSION EAB2(NOCC,NOCC,4),EA2(NOCC,4)
C
C     ISCF TELLS WHICH INTEGRALS TO INCLUDE BASED ON OVERLAP COMBINATION
      DATA ISCF/1,0,1,0,1,0,1,0,1,1,0,
     &          1,0,0,0,1,0,0,0,1,1,0,
     &          0,1,1,0,0,0,0,0,1,1,0,
     &          1,0,0,1,1,0,0,0,1,0,0,
     &          0,1,0,1,0,0,0,0,1,0,0,
     &          0,1,0,0,0,0,0,0,1,0,0/
C
C
C     IMPORT MBPT1 PAIR RESULTS
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_MP1.dat',STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO IOCCA=1,NOCC
        DO IOCCB=1,NOCC
          READ(8, *) Q1,Q2,Q3,Q4,Q5,EAB2(IOCCA,IOCCB,1)
        ENDDO
      ENDDO
      CLOSE(UNIT=8)
C
C     INITIALISE MP2 ENERGY COUNTERS
      DO N=2,4
        DO IOCCA=1,NOCC
          EA2(IOCCA,N) = 0.0D0
          DO IOCCB=1,NOCC
            EAB2(IOCCA,IOCCB,N) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     INDEXING ROUTINE: SO WE CAN SET BASIS FUNCTION LABELS BY BLOCK
      ICOUNT = 0
C
C     LOOP OVER NUCLEAR CENTERS
      DO ICNT=1,NCNT
C
C       LOOP OVER KAPPA VALUES FOR THIS NUCLEAR CENTER
        DO KN=1,NKAP(ICNT)
C
C         IMPORT KAPPA, MAXIMUM MQN
          KAPPA = KVALS(KN,ICNT)
          MJMAX  = 2*IABS(KAPPA)-1
C
C         LOOP OVER MQN VALUES AND RECORD INDEX
          DO MJ=1,MJMAX,2
            ICOUNT = ICOUNT+1
            INDEX(ICNT,KAPPA,MJ) = ICOUNT
          ENDDO
        ENDDO
      ENDDO
C
C     SET THE LIMITS FOR RUNNING OVER DENSITY COMBINATIONS
      IF(HMLTN.EQ.'NORL') THEN
        ITSTRT = 1
        ITSTOP = 1
        ITSKIP = 1
      ELSE
        ITSTRT = 1
        ITSTOP = 4
        ITSKIP = 3
      ENDIF
C
20    FORMAT(1X,A,12X,A,11X,A,11X,A,12X,A)
21    FORMAT(' (',I2,',',I2,')',4X,F13.7,5X,F11.7,5X,F11.7,5X,F13.7)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',24),'MP2 pair-wise summary'
      WRITE(7, *) REPEAT(' ',24),'MP2 pair-wise summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) '( a, b)','E2(H)','E2(J)','E2(K)','E2(ab)'
      WRITE(7,20) '( a, b)','E2(H)','E2(J)','E2(K)','E2(ab)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C**********************************************************************C
C     LOOP OVER UNIQUE OCCUPIED ORBITALS (A,B)       (INDEX 1000)      C
C**********************************************************************C
C
      CALL CPU_TIME(TDM1)
      DO 1000 IOCCA=1,NOCC
      DO 1000 IOCCB=1,IOCCA
C
      IA = IOCCA + NSHIFT
      IB = IOCCB + NSHIFT
C
C     INITIALISE STORAGE ARRAYS FOR (AR|BS) AND (AR|SB)
      DO IR=1,NVIR
        DO IS=1,NVIR
          ABRS(IR,IS) = DCMPLX(0.0D0,0.0D0)
          BARS(IR,IS) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     LOOP OVER UNIQUE VIRTUAL ORBITALS (R,S)        (INDEX 2000)      C
C**********************************************************************C
C
C     GOTO 2001
      DO 2000 IR=1,NVIR
      DO 2000 IS=1,NVIR
C
        IOCCR = NSHIFT + NOCC + IR
        IOCCS = NSHIFT + NOCC + IS
C
C       INITIALISE COULOMB DIRECT AND EXCHANGE ARRAYS
        DO J=1,NDIM
          DO I=1,NDIM
            GDIR(I,J) = DCMPLX(0.0D0,0.0D0)
            GXCH(I,J) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
C
C**********************************************************************C
C     LOOP OVER SPINORS A, B BY BLOCK                (INDEX 3000)      C
C**********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR A AND B: T = {L} OR {L,S}
      DO 3000 IT1=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITQN(1) = IT1
C
C       CALCULATE STARTING ADDRESS
        IF(IT1.EQ.1) THEN
          NADDAB = 0
        ELSE
          NADDAB = NSHIFT
        ENDIF
C
C     LOOP OVER CENTER A
      DO 3000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 3000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
      DO 3000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNA    = NFUNCT(LQN(1)+1,ICNTA)
        NFUNS(1) = NFUNA
C
        DO IBAS=1,NFUNA
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 3000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
         LQN(2) = KQN(2)
        ELSE
         LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNB    = NFUNCT(LQN(2)+1,ICNTB)
        NFUNS(2) = NFUNB
C
        DO JBAS=1,NFUNB
          EXPT(JBAS,2) = EXPSET(JBAS,LQN(2)+1,ICNTB)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 3000 MA=1,IABS(KQN(1))
        MJA    = 2*MA-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 3000 MB=1,IABS(KQN(2))
        MJB    = 2*MB-1
        MQN(2) = MJB
C
C     CALCULATE NEW BLOCK OF E(AB) COEFFS AT NEXT OPPORTUNITY
      IEAB = 1
C
C**********************************************************************C
C     LOOP OVER SPINORS C, D BY BLOCK                (INDEX 4000)      C
C**********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR C AND D: T = {L} OR {L,S}
      DO 4000 IT2=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITQN(2) = IT2
C
C       CALCULATE STARTING ADDRESS
        IF(IT2.EQ.1) THEN
          NADDCD = 0
        ELSE
          NADDCD = NSHIFT
        ENDIF
C
C     LOOP OVER CENTER C
      DO 4000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER C
        XYZ(1,3) = COORD(1,ICNTC)
        XYZ(2,3) = COORD(2,ICNTC)
        XYZ(3,3) = COORD(3,ICNTC)
C
C     LOOP OVER CENTER D
      DO 4000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER D
        XYZ(1,4) = COORD(1,ICNTD)
        XYZ(2,4) = COORD(2,ICNTD)
        XYZ(3,4) = COORD(3,ICNTD)
C
C     LOOP OVER KQN(C) VALUES
      DO 4000 KC=1,NKAP(ICNTC)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR C
        KQN(3) = KVALS(KC,ICNTC)
        IF(KQN(3).GT.0) THEN
          LQN(3) = KQN(3)
        ELSE
          LQN(3) =-KQN(3)-1
        ENDIF
C         
        NFUNC    = NFUNCT(LQN(3)+1,ICNTC)
        NFUNS(3) = NFUNC
C
        DO KBAS=1,NFUNC
          EXPT(KBAS,3) = EXPSET(KBAS,LQN(3)+1,ICNTC)
        ENDDO
C
C     LOOP OVER KQN(D) VALUES
      DO 4000 KD=1,NKAP(ICNTD)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR D
        KQN(4) = KVALS(KD,ICNTD)
        IF(KQN(4).GT.0) THEN
          LQN(4) = KQN(4)
        ELSE
          LQN(4) =-KQN(4)-1
        ENDIF
C
        NFUND    = NFUNCT(LQN(4)+1,ICNTD)
        NFUNS(4) = NFUND
C
        DO LBAS=1,NFUND
          EXPT(LBAS,4) = EXPSET(LBAS,LQN(4)+1,ICNTD)
        ENDDO
C
C     LOOP OVER |MQN(C)| VALUES
      DO 4000 MC=1,IABS(KQN(3))
        MJC    = 2*MC-1
        MQN(3) = MJC
C
C     LOOP OVER |MQN(D)| VALUES
      DO 4000 MD=1,IABS(KQN(4))
        MJD    = 2*MD-1
        MQN(4) = MJD
C
C     CALCULATE NEW BLOCK OF E(CD) COEFFS AT NEXT OPPORTUNITY
      IECD = 1
C
C**********************************************************************C
C     FOR THIS CHOICE OF A,B,C AND D, COMPUTE ADDRESSES AND PHASES     C
C**********************************************************************C
C
C     CALCULATE BLOCK INDICES FOR {ABCD} COMBINATIONS
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
      IQ12 = (IQ1*(IQ1-1))/2 + IQ2
      IQ34 = (IQ3*(IQ3-1))/2 + IQ4
C
C     FURTHER DEFINE STARTING ADDRESSES FOR {ABCD} BASIS LABELS
      IA1 = LARGE(ICNTA,KA,MJA  ) + NADDAB
      IB1 = LARGE(ICNTB,KB,MJB  ) + NADDAB
      IC1 = LARGE(ICNTC,KC,MJC  ) + NADDCD
      ID1 = LARGE(ICNTD,KD,MJD  ) + NADDCD
C
      IA2 = LARGE(ICNTA,KA,MJA+1) + NADDAB
      IB2 = LARGE(ICNTB,KB,MJB+1) + NADDAB
      IC2 = LARGE(ICNTC,KC,MJC+1) + NADDCD
      ID2 = LARGE(ICNTD,KD,MJD+1) + NADDCD
C
      JA1 = LARGE(ICNTA,KA,MJA  ) + NADDAB
      JB1 = LARGE(ICNTB,KB,MJB  ) + NADDAB
      JC1 = LARGE(ICNTC,KC,MJC  ) + NADDCD
      JD1 = LARGE(ICNTD,KD,MJD  ) + NADDCD
C
      JA2 = LARGE(ICNTA,KA,MJA+1) + NADDAB
      JB2 = LARGE(ICNTB,KB,MJB+1) + NADDAB
      JC2 = LARGE(ICNTC,KC,MJC+1) + NADDCD
      JD2 = LARGE(ICNTD,KD,MJD+1) + NADDCD
C
C     CALCULATE KQN PHASE FACTORS FOR PERMUTING INTEGRALS
      IF((KQN(1)*KQN(2)).GT.0) THEN 
        PKAB = 1.0D0
      ELSE
        PKAB =-1.0D0
      ENDIF
C        
      IF((KQN(3)*KQN(4)).GT.0) THEN 
        PKCD = 1.0D0
      ELSE
        PKCD =-1.0D0
      ENDIF
C
C     CALCULATE MQN PHASE FACTORS FOR PERMUTING INTEGRALS
      PMAB1 = DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
      PMAB2 = DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      PMCD1 = DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
      PMCD2 = DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
C
      F1 = PKAB*PMAB1
      F2 = PKAB*PMAB2
      G1 = PKCD*PMCD1
      G2 = PKCD*PMCD2
C
C     INDICATE BLOCKS TO BE INCLUDED AHEAD GIVEN A,B,C,D BASIS QNMS...
C     A =/= B AND C =/= D WITH AB LIST VALUE =/= CD LIST VALUE
      IF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 1
C     A=/=B AND C=/=D WITH IND(AB) = IND(CD)
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.EQ.IQ34) THEN
        ITSCF = 2
C     A=/=B AND C = D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 3
C     A = B AND C=/=D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.GT.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 4
C     A = B AND C = D WITH IND(AB)=/=IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.GT.IQ34) THEN
        ITSCF = 5
C     A = B AND C = D WITH IND(AB) = IND(CD)
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQ12.EQ.IQ34) THEN
        ITSCF = 6
C     COMBINATION OF A,B,C,D NOT TO BE INCLUDED -- USE MATRIX CONJ LATER
      ELSE
        GO TO 4000
      ENDIF
C
C     READ IN FLAG VALUES FROM ISCF DATA BLOCK
      DO M=1,11
        IFLG(M) = ISCF(M,ITSCF)
      ENDDO
C
C     INCLUDE SPECIAL CASES
      IF(ITSCF.EQ.1.AND.IQ1.EQ.IQ3) IFLG(6)  = 1
      IF(ITSCF.EQ.1.AND.IQ2.EQ.IQ3) IFLG(8)  = 1
      IF(ITSCF.EQ.3.AND.IQ2.EQ.IQ3) IFLG(8)  = 1
      IF(ITSCF.EQ.4.AND.IQ2.EQ.IQ3) IFLG(8)  = 1
      IF(ITSCF.EQ.1.AND.IQ2.EQ.IQ4) IFLG(11) = 1
C
C**********************************************************************C
C     FOURTH LAYER OF LOOPS, OVER BASIS FUNCTIONS A AND B (5000)       C
C -------------------------------------------------------------------- C
C     THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH         C
C     GENERATE THE GDIR AND GXCH MATRICES FROM THE SPINOR INTEGRALS.   C
C     THESE INCLUDE PHASE FACTORS FOR THE PERMUTATION OF               C
C     KQN(1) <-> KQN(2) AND MQN(1) <-> MQN(2)                          C
C -------------------------------------------------------------------- C
C     (RSCF 86, 87)                                                    C
C**********************************************************************C
C
C     LOOP OVER ELEMENTS OF FOCK MATRIX BLOCK
      DO 5000 IBAS=1,NFUNA
      DO 5000 JBAS=1,NFUNB
C
C       GENERATE BATCH OF ELECTRON REPULSION INTEGRALS
        CALL ERI(RR,XYZ,KQN,MQN,EXPT,NFUNS,ITQN,IBAS,JBAS,IEAB,IECD)
C
C       THERE ARE ELEVEN DISTINCT PERMUTATIONAL ALGORITHMS WHICH
C       GENERATE THE GDIR/GXCH MATRIX FROM THE SPINOR INTEGRALS.
C
C       FIRST IFLG BATCH (DIRECT)
        IF(IFLG(1).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNC
            DO LBAS=1,NFUND
              M = M+1
C
              GDIR(IA1+IBAS,JB1+JBAS) = GDIR(IA1+IBAS,JB1+JBAS)
     &          +       RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IR)
     &          +       RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IR)
     &          +    G1*RR(M, 1)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M, 2)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M, 3)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IR)
     &          +    G1*RR(M, 4)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IR)
C
              GDIR(IA1+IBAS,JB2+JBAS) = GDIR(IA1+IBAS,JB2+JBAS)
     &          +       RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IR)
     &          +       RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IR)
     &          +    G1*RR(M, 5)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M, 6)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M, 7)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IR)
     &          +    G1*RR(M, 8)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IR)
C
              GDIR(IA2+IBAS,JB1+JBAS) = GDIR(IA2+IBAS,JB1+JBAS)
     &          +       RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IR)
     &          +       RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IR)
     &          +    G1*RR(M, 9)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M,10)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M,11)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IR)
     &          +    G1*RR(M,12)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IR)
C
              GDIR(IA2+IBAS,JB2+JBAS) = GDIR(IA2+IBAS,JB2+JBAS)
     &          +       RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IR)
     &          +       RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IR)
     &          +    G1*RR(M,13)*DCONJG(C(JD2+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M,14)*DCONJG(C(JD1+LBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M,15)*DCONJG(C(JD2+LBAS,IB))*C(IC1+KBAS,IR)
     &          +    G1*RR(M,16)*DCONJG(C(JD1+LBAS,IB))*C(IC1+KBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SECOND IFLG BATCH (DIRECT)
        IF(IFLG(2).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNC
            DO LBAS=1,NFUND
              M = M+1
C
              GDIR(IA1+IBAS,JB1+JBAS) = GDIR(IA1+IBAS,JB1+JBAS)
     &          +       RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IR)
     &          +       RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IR)
C
              GDIR(IA1+IBAS,JB2+JBAS) = GDIR(IA1+IBAS,JB2+JBAS)
     &          +       RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IR)
     &          +       RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IR)
C
              GDIR(IA2+IBAS,JB1+JBAS) = GDIR(IA2+IBAS,JB1+JBAS)
     &          +       RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IR)
     &          +       RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IR)
C
              GDIR(IA2+IBAS,JB2+JBAS) = GDIR(IA2+IBAS,JB2+JBAS)
     &          +       RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JD2+LBAS,IR)
     &          +       RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JD1+LBAS,IR)
     &          +       RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JD2+LBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       THIRD IFLG BATCH (DIRECT)
        IF(IFLG(3).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNC
            DO LBAS=1,NFUND
              M = M+1
C
              GDIR(IC1+KBAS,JD1+LBAS) = GDIR(IC1+KBAS,JD1+LBAS)
     &          +       RR(M, 1)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 5)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M, 9)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,13)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IR)
     &          +    F1*RR(M, 1)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M, 5)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M, 9)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IR)
     &          +    F1*RR(M,13)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IR)
C
              GDIR(IC1+KBAS,JD2+LBAS) = GDIR(IC1+KBAS,JD2+LBAS)
     &          +       RR(M, 2)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 6)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M,10)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,14)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IR)
     &          +    F1*RR(M, 2)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M, 6)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M,10)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IR)
     &          +    F1*RR(M,14)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IR)
C
              GDIR(IC2+KBAS,JD1+LBAS) = GDIR(IC2+KBAS,JD1+LBAS)
     &          +       RR(M, 3)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 7)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M,11)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,15)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IR)
     &          +    F1*RR(M, 3)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M, 7)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M,11)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IR)
     &          +    F1*RR(M,15)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IR)
C
              GDIR(IC2+KBAS,JD2+LBAS) = GDIR(IC2+KBAS,JD2+LBAS)
     &          +       RR(M, 4)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 8)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M,12)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,16)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IR)
     &          +    F1*RR(M, 4)*DCONJG(C(JB2+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M, 8)*DCONJG(C(JB1+JBAS,IB))*C(IA2+IBAS,IR)
     &          +    F2*RR(M,12)*DCONJG(C(JB2+JBAS,IB))*C(IA1+IBAS,IR)
     &          +    F1*RR(M,16)*DCONJG(C(JB1+JBAS,IB))*C(IA1+IBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       FOURTH IFLG BATCH (DIRECT)
        IF(IFLG(4).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNC
            DO LBAS=1,NFUND
              M = M+1
C
              GDIR(IC1+KBAS,JD1+LBAS) = GDIR(IC1+KBAS,JD1+LBAS)
     &          +       RR(M, 1)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 5)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M, 9)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,13)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IR)
C
              GDIR(IC1+KBAS,JD2+LBAS) = GDIR(IC1+KBAS,JD2+LBAS)
     &          +       RR(M, 2)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 6)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M,10)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,14)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IR)
C
              GDIR(IC2+KBAS,JD1+LBAS) = GDIR(IC2+KBAS,JD1+LBAS)
     &          +       RR(M, 3)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 7)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M,11)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,15)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IR)
C
              GDIR(IC2+KBAS,JD2+LBAS) = GDIR(IC2+KBAS,JD2+LBAS)
     &          +       RR(M, 4)*DCONJG(C(IA1+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 8)*DCONJG(C(IA1+IBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M,12)*DCONJG(C(IA2+IBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,16)*DCONJG(C(IA2+IBAS,IB))*C(JB2+JBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       FIFTH IFLG BATCH (EXCHANGE)
        IF(IFLG(5).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNC
            DO LBAS=1,NFUND
              M = M+1
C
              GXCH(IA1+IBAS,JC1+KBAS) = GXCH(IA1+IBAS,JC1+KBAS)
     &          +    G2*RR(M, 3)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G1*RR(M, 4)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G2*RR(M, 7)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IR)
     &          +    G1*RR(M, 8)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IR)
C
              GXCH(IA1+IBAS,JC2+KBAS) = GXCH(IA1+IBAS,JC2+KBAS)
     &          +    G1*RR(M, 1)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G2*RR(M, 2)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G1*RR(M, 5)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IR)
     &          +    G2*RR(M, 6)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IR)
C
              GXCH(IA2+IBAS,JC1+KBAS) = GXCH(IA2+IBAS,JC1+KBAS)
     &          +    G2*RR(M,11)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G1*RR(M,12)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G2*RR(M,15)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IR)
     &          +    G1*RR(M,16)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IR)
C
              GXCH(IA2+IBAS,JC2+KBAS) = GXCH(IA2+IBAS,JC2+KBAS)
     &          +    G1*RR(M, 9)*DCONJG(C(ID2+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G2*RR(M,10)*DCONJG(C(ID1+LBAS,IB))*C(JB1+JBAS,IR)
     &          +    G1*RR(M,13)*DCONJG(C(ID2+LBAS,IB))*C(JB2+JBAS,IR)
     &          +    G2*RR(M,14)*DCONJG(C(ID1+LBAS,IB))*C(JB2+JBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SIXTH IFLG BATCH (EXCHANGE)
        IF(IFLG(6).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNC
            DO LBAS=1,NFUND
              M = M+1
C
              GXCH(IC1+KBAS,JA1+IBAS) = GXCH(IC1+KBAS,JA1+IBAS)
     &          +    F2*RR(M, 9)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F2*RR(M,10)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IR)
     &          +    F1*RR(M,13)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F1*RR(M,14)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IR)
C
              GXCH(IC1+KBAS,JA2+IBAS) = GXCH(IC1+KBAS,JA2+IBAS)
     &          +    F1*RR(M, 1)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F1*RR(M, 2)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IR)
     &          +    F2*RR(M, 5)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F2*RR(M, 6)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IR)
C
              GXCH(IC2+KBAS,JA1+IBAS) = GXCH(IC2+KBAS,JA1+IBAS)
     &          +    F2*RR(M,11)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F2*RR(M,12)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IR)
     &          +    F1*RR(M,15)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F1*RR(M,16)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IR)
C
              GXCH(IC2+KBAS,JA2+IBAS) = GXCH(IC2+KBAS,JA2+IBAS)
     &          +    F1*RR(M, 3)*DCONJG(C(JB2+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F1*RR(M, 4)*DCONJG(C(JB2+JBAS,IB))*C(ID2+LBAS,IR)
     &          +    F2*RR(M, 7)*DCONJG(C(JB1+JBAS,IB))*C(ID1+LBAS,IR)
     &          +    F2*RR(M, 8)*DCONJG(C(JB1+JBAS,IB))*C(ID2+LBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       SEVENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(7).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNC
            DO LBAS=1,NFUND
              M = M+1
C
              GXCH(IB1+JBAS,JC1+KBAS) = GXCH(IB1+JBAS,JC1+KBAS)
     &          + F2*G2*RR(M, 7)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F2*G1*RR(M, 8)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F1*G2*RR(M,15)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IR)
     &          + F1*G1*RR(M,16)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IR)
C
              GXCH(IB1+JBAS,JC2+KBAS) = GXCH(IB1+JBAS,JC2+KBAS)
     &          + F2*G1*RR(M, 5)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F2*G2*RR(M, 6)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F1*G1*RR(M,13)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IR)
     &          + F1*G2*RR(M,14)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IR)
C
              GXCH(IB2+JBAS,JC1+KBAS) = GXCH(IB2+JBAS,JC1+KBAS)
     &          + F1*G2*RR(M, 3)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F1*G1*RR(M, 4)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F2*G2*RR(M,11)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IR)
     &          + F2*G1*RR(M,12)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IR)
C
              GXCH(IB2+JBAS,JC2+KBAS) = GXCH(IB2+JBAS,JC2+KBAS)
     &          + F1*G1*RR(M, 1)*DCONJG(C(ID2+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F1*G2*RR(M, 2)*DCONJG(C(ID1+LBAS,IB))*C(JA2+IBAS,IR)
     &          + F2*G1*RR(M, 9)*DCONJG(C(ID2+LBAS,IB))*C(JA1+IBAS,IR)
     &          + F2*G2*RR(M,10)*DCONJG(C(ID1+LBAS,IB))*C(JA1+IBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       EIGHTH IFLG BATCH (EXCHANGE)
        IF(IFLG(8).EQ.1) THEN
          M = 0
          DO KBAS=1,NFUNC
            DO LBAS=1,NFUND
              M = M+1
C
              GXCH(IC1+KBAS,JB1+JBAS) = GXCH(IC1+KBAS,JB1+JBAS)
     &          +       RR(M, 1)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M, 2)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IR)
     &          +       RR(M, 9)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M,10)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IR)
C
              GXCH(IC1+KBAS,JB2+JBAS) = GXCH(IC1+KBAS,JB2+JBAS)
     &          +       RR(M, 5)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M, 6)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IR)
     &          +       RR(M,13)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M,14)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IR)
C
              GXCH(IC2+KBAS,JB1+JBAS) = GXCH(IC2+KBAS,JB1+JBAS)
     &          +       RR(M, 3)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M, 4)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IR)
     &          +       RR(M,11)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M,12)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IR)
C
              GXCH(IC2+KBAS,JB2+JBAS) = GXCH(IC2+KBAS,JB2+JBAS)
     &          +       RR(M, 7)*DCONJG(C(JA1+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M, 8)*DCONJG(C(JA1+IBAS,IB))*C(ID2+LBAS,IR)
     &          +       RR(M,15)*DCONJG(C(JA2+IBAS,IB))*C(ID1+LBAS,IR)
     &          +       RR(M,16)*DCONJG(C(JA2+IBAS,IB))*C(ID2+LBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       NINTH IFLG BATCH (EXCHANGE)
        IF(IFLG(9).EQ.1) THEN
          DO LBAS=1,NFUND
            DO KBAS=1,NFUNC
              M = (KBAS-1)*NFUND+LBAS
C
              GXCH(IA1+IBAS,JD1+LBAS) = GXCH(IA1+IBAS,JD1+LBAS)
     &          +       RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IR)
C
              GXCH(IA1+IBAS,JD2+LBAS) = GXCH(IA1+IBAS,JD2+LBAS)
     &          +       RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IR)
C
              GXCH(IA2+IBAS,JD1+LBAS) = GXCH(IA2+IBAS,JD1+LBAS)
     &          +       RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IR)
C
              GXCH(IA2+IBAS,JD2+LBAS) = GXCH(IA2+IBAS,JD2+LBAS)
     &          +       RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JB1+JBAS,IR)
     &          +       RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JB2+JBAS,IR)
     &          +       RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JB2+JBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       TENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(10).EQ.1) THEN
          DO LBAS=1,NFUND
            DO KBAS=1,NFUNC
              M = (KBAS-1)*NFUND+LBAS
C
              GXCH(IB1+JBAS,JD1+LBAS) = GXCH(IB1+JBAS,JD1+LBAS)
     &          +    F2*RR(M, 5)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F2*RR(M, 7)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F1*RR(M,13)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IR)
     &          +    F1*RR(M,15)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IR)
C
              GXCH(IB1+JBAS,JD2+LBAS) = GXCH(IB1+JBAS,JD2+LBAS)
     &          +    F2*RR(M, 6)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F2*RR(M, 8)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F1*RR(M,14)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IR)
     &          +    F1*RR(M,16)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IR)
C
              GXCH(IB2+JBAS,JD1+LBAS) = GXCH(IB2+JBAS,JD1+LBAS)
     &          +    F1*RR(M, 1)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F1*RR(M, 3)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F2*RR(M, 9)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IR)
     &          +    F2*RR(M,11)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IR)
C
              GXCH(IB2+JBAS,JD2+LBAS) = GXCH(IB2+JBAS,JD2+LBAS)
     &          +    F1*RR(M, 2)*DCONJG(C(IC1+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F1*RR(M, 4)*DCONJG(C(IC2+KBAS,IB))*C(JA2+IBAS,IR)
     &          +    F2*RR(M,10)*DCONJG(C(IC1+KBAS,IB))*C(JA1+IBAS,IR)
     &          +    F2*RR(M,12)*DCONJG(C(IC2+KBAS,IB))*C(JA1+IBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C       ELEVENTH IFLG BATCH (EXCHANGE)
        IF(IFLG(11).EQ.1) THEN
          DO LBAS=1,NFUND
            DO KBAS=1,NFUNC
              M = (KBAS-1)*NFUND+LBAS
C
              GXCH(ID1+LBAS,JB1+JBAS) = GXCH(ID1+LBAS,JB1+JBAS)
     &          +    G2*RR(M, 2)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G1*RR(M, 4)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IR)
     &          +    G2*RR(M,10)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G1*RR(M,12)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IR)
C
              GXCH(ID1+LBAS,JB2+JBAS) = GXCH(ID1+LBAS,JB2+JBAS)
     &          +    G2*RR(M, 6)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G1*RR(M, 8)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IR)
     &          +    G2*RR(M,14)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G1*RR(M,16)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IR)
C
              GXCH(ID2+LBAS,JB1+JBAS) = GXCH(ID2+LBAS,JB1+JBAS)
     &          +    G1*RR(M, 1)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M, 3)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IR)
     &          +    G1*RR(M, 9)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M,11)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IR)
C
              GXCH(ID2+LBAS,JB2+JBAS) = GXCH(ID2+LBAS,JB2+JBAS)
     &          +    G1*RR(M, 5)*DCONJG(C(JA1+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M, 7)*DCONJG(C(JA1+IBAS,IB))*C(IC1+KBAS,IR)
     &          +    G1*RR(M,13)*DCONJG(C(JA2+IBAS,IB))*C(IC2+KBAS,IR)
     &          +    G2*RR(M,15)*DCONJG(C(JA2+IBAS,IB))*C(IC1+KBAS,IR)
C
            ENDDO
          ENDDO
        ENDIF
C
C     END LOOP OVER IAOCC,IBOCC BLOCK ADDRESSES
5000  CONTINUE
C     END LOOPS OVER A,B,C,D OVERLAP BLOCKS
4000  CONTINUE
3000  CONTINUE
C
C     COMPLETE CONSTRUCTION OF GDIR AND GXCH BY MATRIX CONJUGATION
      DO J=1,NDIM-NSHIFT
        DO I=1,J
C
C         SMALL-COMPONENT ADDRESSES
          K = I + NSHIFT
          L = J + NSHIFT
C
C         SKIP DIAGONAL PARTS OF EACH SUB-BLOCK
          IF(ICNBAS(I).NE.ICNBAS(J)) GOTO 400
          IF(KQNBAS(I).NE.KQNBAS(J)) GOTO 400
          IF(IABS(MQNBAS(I)).NE.IABS(MQNBAS(J))) GOTO 400
          GOTO 401
400       CONTINUE
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF LL BLOCK
          GDIR(I,J) = GDIR(I,J) + DCONJG(GDIR(J,I))
          GDIR(J,I) =             DCONJG(GDIR(I,J))
          GXCH(I,J) = GXCH(I,J) + DCONJG(GXCH(J,I))
          GXCH(J,I) =             DCONJG(GXCH(I,J))
C
C         IF HMLTN = 'NORL' SKIP THE NEXT FEW CALCULATIONS
          IF(HMLTN.EQ.'NORL') GOTO 401
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF LS BLOCK
          GDIR(I,L) = GDIR(I,L) + DCONJG(GDIR(L,I))
          GDIR(L,I) =             DCONJG(GDIR(I,L))
          GXCH(I,L) = GXCH(I,L) + DCONJG(GXCH(L,I))
          GXCH(L,I) =             DCONJG(GXCH(I,L))
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF SL BLOCK
          GDIR(K,J) = GDIR(K,J) + DCONJG(GDIR(J,K))
          GDIR(J,K) =             DCONJG(GDIR(K,J))
          GXCH(K,J) = GXCH(K,J) + DCONJG(GXCH(J,K))
          GXCH(J,K) =             DCONJG(GXCH(K,J))
C
C         COMPLETE LOWER AND THEN UPPER TRIANGLE OF SS BLOCK
          GDIR(K,L) = GDIR(K,L) + DCONJG(GDIR(L,K))
          GDIR(L,K) =             DCONJG(GDIR(K,L))
          GXCH(K,L) = GXCH(K,L) + DCONJG(GXCH(L,K))
          GXCH(L,K) =             DCONJG(GXCH(K,L))
C
401       CONTINUE
        ENDDO
      ENDDO
C
C     END LOOP OVER UNIQUE VIRTUAL ORBITAL COMBINATIONS (R,S)
2000  CONTINUE     
C2001  CONTINUE
C
C**********************************************************************C
C     CALCULATE INTERACTION ENERGY FROM ORBITAL PAIR (A,B)             C
C**********************************************************************C
C
C     ALL MATRIX ELEMENTS (AR|BS) AND (AR|SB) ARE STORED FOR SUMMATION
      ETMP3 = DCMPLX(0.0D0,0.0D0)
      ETMP4 = DCMPLX(0.0D0,0.0D0)
      DO IOCCR=1,NVIR
        DO IOCCS=1,NVIR
C
          IR = NSHIFT + NOCC + IOCCR
          IS = NSHIFT + NOCC + IOCCS
C
          RNUMJ = DCONJG(ABRS(IOCCR,IOCCS))*ABRS(IOCCR,IOCCS)
          RNUMK = DCONJG(BARS(IOCCR,IOCCS))*ABRS(IOCCR,IOCCS)
          RDEN  = EIGEN(IA)+EIGEN(IB)-EIGEN(IR)-EIGEN(IS)
          
          ETMP3 = ETMP3 + RNUMJ/RDEN
          ETMP4 = ETMP4 + RNUMK/RDEN
        ENDDO
      ENDDO
      EAB2(IOCCA,IOCCB,2) = DREAL(ETMP3)
      EAB2(IOCCA,IOCCB,3) = DREAL(ETMP4)
      EAB2(IOCCA,IOCCB,4) = DREAL(ETMP3 - ETMP4)
C
C     OUTPUT ENERGIES TO TERMINAL
      WRITE(6,21) IOCCA,IOCCB,(EAB2(IOCCA,IOCCB,N),N=1,4)
      WRITE(7,21) IOCCA,IOCCB,(EAB2(IOCCA,IOCCB,N),N=1,4)
C
      IF(IOCCA.NE.IOCCB) THEN
        DO N=1,4
          EAB2(IOCCB,IOCCA,N) = EAB2(IOCCA,IOCCB,N)
        ENDDO
      ENDIF
C
C     END LOOP OVER UNIQUE OCCUPIED ORBITAL COMBINATIONS (A,B)
1000  CONTINUE
C
C**********************************************************************C
C     END OF CALCULATION OVER OCCUPIED ORBITALS (A,B)                  C
C**********************************************************************C
C
C     RECORD TIME AT THE END OF MAIN CALCULATION
      CALL CPU_TIME(TDM2)
C
C     WRITE RESULTS OF MP2 ENERGIES TO AN EXTERNAL FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_MP2.dat',STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO IOCCA=1,NOCC
        DO IOCCB=1,NOCC
          WRITE(8, *) (EAB2(IOCCA,IOCCB,N),N=1,4)
        ENDDO
      ENDDO
      CLOSE(UNIT=8)
C
C     CALCULATE MBPT1 SINGLE-PARTICLE ENERGIES AND MOLECULAR TOTALS
      EDIR2 = 0.0D0
      EXCH2 = 0.0D0
      DO IOCCA=1,NOCC
        DO N=1,3
          EA2(IOCCA,N) = 0.0D0
          DO IOCCB=1,NOCC
            EA2(IOCCA,N) = EA2(IOCCA,N) + EAB2(IOCCA,IOCCB,N)
          ENDDO      
        ENDDO
        EA2(IOCCA,4) = EA2(IOCCA,2) - EA2(IOCCA,3)
        EDIR2 = EDIR2 + EA2(IOCCA,2)*0.5D0
        EXCH2 = EXCH2 + EA2(IOCCA,3)*0.5D0
      ENDDO
      ETOT2 = EDIR2 - EXCH2
C
C     ORBITAL SUMMARIES
30    FORMAT(1X,A,12X,A,11X,A,11X,A,12X,A)
31    FORMAT('  ',I2,'    ',4X,F13.7,5X,F11.7,5X,F11.7,5X,F13.7)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',21),'MP2 single particle summary'
      WRITE(7, *) REPEAT(' ',21),'MP2 single particle summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,30) '  a    ','E1(a)','E2(J)','E2(K)',' E2(a)'
      WRITE(7,30) '  a    ','E1(a)','E2(J)','E2(K)',' E2(a)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)  
      DO IOCCA=1,NOCC
        WRITE(6,31) IOCCA,(EA2(IOCCA,N),N=1,4)
        WRITE(7,31) IOCCA,(EA2(IOCCA,N),N=1,4)
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     TOTAL ENERGIES
32    FORMAT(1X,A,4X,A,2X,'=',10X,F18.8,' au')
33    FORMAT(1X,A,5X,'=',23X,A)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',20),'MP2 molecular summary'
      WRITE(7, *) REPEAT(' ',20),'MP2 molecular summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,32) 'Correlation Coulomb direct   ','E2(J)',EDIR2
      WRITE(7,32) 'Correlation Coulomb direct   ','E2(J)',EDIR2
      WRITE(6,32) 'Correlation Coulomb exchange ','E2(K)',EXCH2
      WRITE(7,32) 'Correlation Coulomb exchange ','E2(K)',EXCH2
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,32) 'Correlation Coulomb total    ','E2(G)',EDIR2-EXCH2
      WRITE(7,32) 'Correlation Coulomb total    ','E2(G)',EDIR2-EXCH2
      WRITE(6,32) 'Hartree-Fock molecular energy','E2   ',ETOT+ETOT2
      WRITE(7,32) 'Hartree-Fock molecular energy','E2   ',ETOT+ETOT2
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,33) 'MP2 time                   ',HMS(TDM2-TDM1)
      WRITE(7,33) 'MP2 time                   ',HMS(TDM2-TDM1)
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C   [9] EXPECTATION VALUES: OBSERVABLES FROM A CONVERGED SOLUTION.     C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] PT1BODY: MAIN ROUTINE FOR MOLECULAR EXPECTATION VALUES.        C
C   [B] PROPRTY: MOLECULAR EXPECTATION VALUE FROM DENSITY MATRIX.      C
C   [B] RS1: SET OF 1ST ORDER MOLECULAR MATRIX ELEMENTS AND E(2).      C
C   [C] RS2: SET OF 2ND ORDER MOLECULAR MATRIX ELEMENTS AND E(3).      C
C -------------------------------------------------------------------- C
C   [D] MULLIKN: MULLIKEN POPULATION ANALYSIS ON CONVERGED SOLUTION.   C
C   [E] ELCMNPL: MOLECULAR ELECTRIC MONOPOLE MOMENT ANALYSIS.          C
C   [F] ELCDIPL: MOLECULAR ELECTRIC DIPOLE MOMENT ANALYSIS.            C
C   [G] ELCQDPL: MOLECULAR ELECTRIC QUADRUPOLE MOMENT ANALYSIS.        C
C   [H] MAGDIPL: MOLECULAR MAGNETIC DIPOLE MOMENT ANALYSIS.            C
C   [I] MAGQDPL: MOLECULAR MAGNETIC QUADRUPOLE MOMENT ANALYSIS.        C
C   [J] STRKEFF: STARK EFFECT ANALYSIS, GIVEN ELECTRIC FIELD E.        C
C   [K] ZMANEFF: ZEEMAN EFFECT ANALYSIS, GIVEN MAGNETIC FIELD B.       C
C   [L] GTENSOR: MAGNETIC G-TENSOR CALCULATION.                        C
C   [M] HYPFINE: HYPERFINE INTERACTION ANALYSIS, GIVEN NUCLEAR MOMENT. C
C   [K] EEDMSML: ATOM-CENTERED PT-ODD EDM OPERATOR (WITH E-FIELD).     C
C   [J] EEDMEFF: ONE-BODY EFFECTIVE PT-ODD EDM OPERATOR.               C
C   [N] SCLPTEN: SCALAR PT-ODD ELECTRON-NUCLEAR INTERACTION ANALYSIS.  C
C   [O] VECPTEN: VECTOR PT-ODD ELECTRON-NUCLEAR INTERACTION ANALYSIS.  C
C   [P] PVIOLTN: P-ODD EFFECTIVE OPERATOR ANALYSIS.                    C
C   [Q] BETADCY: CORRECTIONS DUE TO THE NUCLEAR DECAY OF A CENTER.     C
C -------------------------------------------------------------------- C
C   [R] VMNPOLE: BASIS MONOPOLE MOMENT MATRIX OVER SIGMA_Q.            C
C   [S] VDIPOLE: BASIS DIPOLE MOMENT MATRIX OVER SIGMA_Q AND IX.       C
C   [T] VQDPOLE: BASIS QUADRUPOLE MOMENT MATRIX OVER SIGMA_Q, IX, JX.  C
C   [U] VEFIELD: BASIS ELECTRIC FIELD MATRIX OVER SIGMA_Q AND IX.      C
C   [V] VKNETIC: BASIS RELATIVISTIC KINETIC OVERLAP MATRIX.            C
C   [W] VLPLACE: BASIS NON-RELATIVISTIC KINETIC OVERLAP MATRIX.        C
C   [X] VNCATRC: BASIS NUCLEAR ATTRACTION OVERLAP MATRIX.              C
C**********************************************************************C
C
C
      SUBROUTINE PT1BODY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      PPPPPPP TTTTTTTT  11  BBBBBBB   OOOOOO  DDDDDDD  YY    YY       C
C      PP    PP   TT    111  BB    BB OO    OO DD    DD YY    YY       C
C      PP    PP   TT     11  BB    BB OO    OO DD    DD  YY  YY        C
C      PP    PP   TT     11  BBBBBBB  OO    OO DD    DD   YYYY         C
C      PPPPPPP    TT     11  BB    BB OO    OO DD    DD    YY          C
C      PP         TT     11  BB    BB OO    OO DD    DD    YY          C
C      PP         TT    1111 BBBBBBB   OOOOOO  DDDDDDD     YY          C
C                                                                      C
C -------------------------------------------------------------------- C
C  PT1BODY CALCULATES MATRIX ELEMENTS AND ENERGY CORRECTIONS GIVEN A   C
C  HARTREE-FOCK CALCULATION AND SET OF COEFFICIENTS AND ENERGIES.      C
C -------------------------------------------------------------------- C
C  EQFILE ONLY GENERATES ELL0, ESS0 AND MAYBE ELSI, AND IN GENERAL     C
C  THESE CALCULATIONS REQUIRE ETT'Q, SO GENERATE AS NEEDED INSTEAD.    C
C -------------------------------------------------------------------- C
C  BASIS OVERLAP MATRIX ELEMENTS AVAILABLE:                            C
C  (u,T|s_q|v,T')     - VMNPOLE(VIJ,IQ)      - DIRECT OVERLAP          C
C  (u,T|s.x|v,T')     - VDIPOLE(VIJ,IQ,IX)   - ELECTRIC DIPOLE MOMENT  C
C  (u,T|s.x.x'|v,T')  - VQDPOLE(VIJ,IQ,IX,JX)- ELECTRIC QUADRUPOLE     C
C  (u,T|s.x/r^3|v,T') - VEFIELD(VIJ,IQ,IX)   - E FIELD/B DIPOLE        C
C  (u,T|s.p|v,T')     - VKNETIC(VIJ)         - KINETIC                 C
C  (u,T|lap|v,T')     - VLPLACE(VIJ)         - NON-REL KINETIC         C
C  (u,T|nuc|v,T')     - VNCATRC(VIJ)         - NUCLEAR ATTRACTION      C
C**********************************************************************C
C
      CHARACTER*4  HMLTN
      CHARACTER*7  HMINT(10)
      CHARACTER*16 HMS
C
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/PT1B/NHMINT,HMINT
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
C     RECORD TIME
      CALL CPU_TIME(TDUM)
C
C     IF READING IN PREVIOUS SOLUTION, CALCULATE MOLECULAR DENSITY
      IF(INEW.EQ.1) THEN
        CALL DENSTY
      ENDIF
C
C     EQ-COEFF AND R-INT TIME INITIALISATION
      TELL = 0.0D0
      TESS = 0.0D0
      TELS = 0.0D0
C
C     LOOP OVER ALL REQUESTED INTERACTION HAMILTONIANS
      DO N=1,NHMINT

C       PRINT A TITLE
        WRITE(6, *) ' '
        WRITE(7, *) ' '
        WRITE(6, *) 'One-body H_{int} = ',HMINT(N)
        WRITE(7, *) 'One-body H_{int} = ',HMINT(N)
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
C
        IF(HMINT(N).EQ.'MULLIKN') THEN
          CALL MULLIKN(6)
        ELSEIF(HMINT(N).EQ.'ELCMNPL') THEN
          CALL ELCMNPL
        ELSEIF(HMINT(N).EQ.'ELCDIPL') THEN
          CALL ELCDIPL
        ELSEIF(HMINT(N).EQ.'ELCQDPL') THEN
          CALL ELCQDPL
        ELSEIF(HMINT(N).EQ.'MAGDIPL') THEN
          CALL MAGDIPL
        ELSEIF(HMINT(N).EQ.'MAGQDPL') THEN
          CALL MAGQDPL
        ELSEIF(HMINT(N).EQ.'STRKEFF') THEN
          CALL STRKEFF
        ELSEIF(HMINT(N).EQ.'ZMANEFF') THEN
          CALL ZMANEFF
        ELSEIF(HMINT(N).EQ.'GTENSOR') THEN
          CALL GTENSOR
        ELSEIF(HMINT(N).EQ.'HYPFINE') THEN
          CALL HYPFINE
        ELSEIF(HMINT(N).EQ.'EEDMSML') THEN
          CALL EEDMSML
        ELSEIF(HMINT(N).EQ.'EEDMEFF') THEN
          CALL EEDMEFF
        ELSEIF(HMINT(N).EQ.'SCLPTEN') THEN
          CALL SCLPTEN
        ELSEIF(HMINT(N).EQ.'VECPTEN') THEN
          CALL VECPTEN
        ELSEIF(HMINT(N).EQ.'PVIOLTN') THEN
          CALL PVIOLTN
        ELSEIF(HMINT(N).EQ.'BETADCY') THEN
          ICNT = 1
          ZDCY = 2.0D0
          CALL BETADCY(ICNT,ZDCY)
        ELSE
          WRITE(6, *) 'In PT1BODY: this operator is not available.'
          WRITE(7, *) 'In PT1BODY: this operator is not available.'
        ENDIF
C
C     END LOOP OVER INTERACTION HAMILTONIANS
      ENDDO
C
C     TOTAL TIME TAKEN
      CALL CPU_TIME(TPRP)
      TPRP = TPRP-TDUM
C
20    FORMAT(1X,A,37X,A)
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) 'Time in EMAKE (LL):',HMS(TELL)
      WRITE(7,20) 'Time in EMAKE (LL):',HMS(TELL)
      WRITE(6,20) 'Time in EMAKE (SS):',HMS(TESS)
      WRITE(7,20) 'Time in EMAKE (SS):',HMS(TESS)
      WRITE(6,20) 'Time in EMAKE (LS):',HMS(TELS)
      WRITE(7,20) 'Time in EMAKE (LS):',HMS(TELS)
      WRITE(6,*) REPEAT('=',72)
      WRITE(7,*) REPEAT('=',72)
C
      RETURN
      END
C
C
      SUBROUTINE PROPRTY(TOT,BLL,BSS,BLS,BSL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     PPPPPPP  RRRRRRR   OOOOOO  PPPPPPP  RRRRRRR TTTTTTTT YY    YY    C
C     PP    PP RR    RR OO    OO PP    PP RR    RR   TT    YY    YY    C
C     PP    PP RR    RR OO    OO PP    PP RR    RR   TT     YY  YY     C
C     PP    PP RR    RR OO    OO PP    PP RR    RR   TT      YYYY      C
C     PPPPPPP  RRRRRRR  OO    OO PPPPPPP  RRRRRRR    TT       YY       C
C     PP       RR    RR OO    OO PP       RR    RR   TT       YY       C
C     PP       RR    RR  OOOOOO  PP       RR    RR   TT       YY       C
C                                                                      C
C -------------------------------------------------------------------- C
C  PROPRTY CALCULATES A MOLECULAR EXPECTATION VALUE OVER THE DENSITY   C
C  MATRIX AND BASIS FUNCTION OVERLAPS IN THE B ARRAYS.                 C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)

      COMPLEX*16 SUMLL,SUMSS,SUMLS,SUMSL,TOT
      COMPLEX*16 BLL(MDM,MDM),BSS(MDM,MDM),BLS(MDM,MDM),BSL(MDM,MDM)
C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C
C     INITIALISE COUNTER ARRAYS FOR BASIS OVERLAP CONTRIBUTIONS
      SUMLL = DCMPLX(0.0D0,0.0D0)
      SUMSS = DCMPLX(0.0D0,0.0D0)
      SUMLS = DCMPLX(0.0D0,0.0D0)
      SUMSL = DCMPLX(0.0D0,0.0D0)
          
C     LOOP OVER ALL BASIS FUNCTIONS
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          SUMLL = SUMLL + BLL(I,J)*DENT(I       ,J       )
          IF(HMLTN.EQ.'NORL') GOTO 100
          SUMSS = SUMSS + BSS(I,J)*DENT(I+NSHIFT,J+NSHIFT)
          SUMLS = SUMLS + BLS(I,J)*DENT(I       ,J+NSHIFT)
          SUMSL = SUMSL + BSL(I,J)*DENT(I+NSHIFT,J       )
100       CONTINUE
        ENDDO
      ENDDO
C
C     MOLECULAR EXPECTATION VALUE
      TOT = SUMLL + SUMSS + SUMLS + SUMSL
C
      RETURN
      END
C
C
      SUBROUTINE RS1(V1,E2,BLL,BSS,BLS,BSL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                        RRRRRRR   SSSSSS   11                         C
C                        RR    RR SS    SS 111                         C
C                        RR    RR SS        11                         C
C                        RR    RR  SSSSSS   11                         C
C                        RRRRRRR        SS  11                         C
C                        RR    RR SS    SS  11                         C
C                        RR    RR  SSSSSS  1111                        C
C                                                                      C
C -------------------------------------------------------------------- C
C  RS1 ASSEMBLES AN ARRAY OF MOLECULAR MATRIX ELEMENTS FROM THE BASIS  C
C  FUNCTION OVERLAPS IN BTT BY FIRST-ORDER RAYLEIGH-SCHRODINGER THEORY.C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 SUMLL,SUMSS,SUMLS,SUMSL,TMP
      COMPLEX*16 C(MDM,MDM)
      COMPLEX*16 BLL(MDM,MDM),BSS(MDM,MDM),BLS(MDM,MDM),BSL(MDM,MDM)
      COMPLEX*16 V1(MDM,MDM),E2(MDM)
C
      COMMON/COEF/C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     TOLERANCE VALUE FOR VANISHING MATRIX ELEMENTS
      TOL = 1.0D-10
C
C     LOOP OVER ALL ORBITAL COMBINATIONS (NEGATIVE AND POSITIVE ENERGY)
      DO IOCC=1,NDIM
        DO JOCC=1,NDIM
C
C         INITIALISE COUNTER ARRAYS FOR BASIS OVERLAP CONTRIBUTIONS
          SUMLL = DCMPLX(0.0D0,0.0D0)
          SUMSS = DCMPLX(0.0D0,0.0D0)
          SUMLS = DCMPLX(0.0D0,0.0D0)
          SUMSL = DCMPLX(0.0D0,0.0D0)
C
C         LOOP OVER FOCK MATRIX ADDRESSES
          DO I=1,NDIM-NSHIFT
            DO J=1,NDIM-NSHIFT
C
              K = I+NSHIFT
              L = J+NSHIFT
C
C             LARGE AND SMALL CONTRIBUTIONS
              SUMLL = SUMLL + DCONJG(C(I,IOCC))*C(J,JOCC)*BLL(I,J)
              IF(HMLTN.EQ.'NORL') GOTO 100
              SUMSS = SUMSS + DCONJG(C(K,IOCC))*C(L,JOCC)*BSS(I,J)
              SUMLS = SUMLS + DCONJG(C(I,IOCC))*C(L,JOCC)*BLS(I,J)
              SUMSL = SUMSL + DCONJG(C(K,IOCC))*C(J,JOCC)*BSL(I,J)
100           CONTINUE
C
            ENDDO
          ENDDO
C
C         SAVE SUMS TO MATRIX
          V1(IOCC,JOCC) = SUMLL + SUMSS + SUMLS + SUMSL
C
C         VANISHING MATRIX ELEMENTS
          TMP1 = DREAL(V1(IOCC,JOCC))
          TMP2 = DIMAG(V1(IOCC,JOCC))
          IF(DABS(TMP1).LT.TOL) THEN
            TMP1 = 0.0D0
          ENDIF
          IF(DABS(TMP2).LT.TOL) THEN
            TMP2 = 0.0D0
          ENDIF
          V1(IOCC,JOCC) = DCMPLX(TMP1,TMP2)
C
        ENDDO
      ENDDO
C
C     CALCULATE ENERGY(2) CORRECTION FOR EACH ORBITAL
      DO IOCC=1,NDIM
        TMP = DCMPLX(0.0D0,0.0D0)
        DO JOCC=1,NDIM
          IF(JOCC.EQ.IOCC) GOTO 50
          IF(DABS(EIGEN(IOCC)-EIGEN(JOCC)).LT.1.0D-12) GOTO 50
          RN = DREAL(DCONJG(V1(IOCC,JOCC))*V1(IOCC,JOCC))
          RD = EIGEN(IOCC) - EIGEN(JOCC)
          TMP = TMP + RN/RD
50        CONTINUE
        ENDDO
        E2(IOCC) = TMP
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE RS2(V1,V2,E3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                      RRRRRRR   SSSSSS   222222                       C
C                      RR    RR SS    SS 22    22                      C
C                      RR    RR SS             22                      C
C                      RR    RR  SSSSSS      22                        C
C                      RRRRRRR        SS   22                          C
C                      RR    RR SS    SS 22                            C
C                      RR    RR  SSSSSS  22222222                      C
C                                                                      C
C -------------------------------------------------------------------- C
C  RS2 APPLIES 2ND ORDER RAYLEIGH-SCHRODINGER PERTURBATION THEORY TO A C
C  SET OF 1ST ORDER MATRIX ELEMENTS IN V1, AND OUTPUTS RESULTS TO V2.  C
C -------------------------------------------------------------------- C
C  DUE TO THE 2-FOLD MQN SYMMETRY, EITHER APPLY DEGENERATE RSPT OR     C
C  SIMPLY IGNORE THE RELEVANT PAIR ORBITAL.                            C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 V1(MDM,MDM),V2(MDM,MDM),E3(MDM)
C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     TOLERANCE VALUE FOR VANISHING MATRIX ELEMENTS
      TOL = 1.0D-10
C
C     LOOP OVER ALL ORBITAL COMBINATIONS (NEGATIVE AND POSITIVE ENERGY)
      DO IOCC=1,NDIM
        DO JOCC=1,NDIM
C
C         APPLY 2ND ORDER RSPT
          V1(IOCC,JOCC) = DCMPLX(0.0D0,0.0D0)
C
C         VANISHING MATRIX ELEMENTS
          TMP1 = DREAL(V2(IOCC,JOCC))
          TMP2 = DIMAG(V2(IOCC,JOCC))
          IF(DABS(TMP1).LT.TOL) THEN
            TMP1 = 0.0D0
          ENDIF
          IF(DABS(TMP2).LT.TOL) THEN
            TMP2 = 0.0D0
          ENDIF
          V2(IOCC,JOCC) = DCMPLX(TMP1,TMP2)
C
        ENDDO
      ENDDO
C      
C     THIRD ORDER ENERGY CORRECTION
      DO IOCC=1,NDIM
        E3(IOCC) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE MULLIKN(IVIR)
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
C  CONTRIBUTIONS FOR EACH ORBITAL ARE GROUPED BY ATOMIC CENTER, KQN    C
C  SYMMETRY AND MQN. ANY CHARGE DENSITY OTHER THAN THIS IS IN 'BRD'.   C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    IVIR - NUMBER OF VIRTUAL ORBITALS TO INCLUDE IN LIST.             C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*2 ELMNT(120),ELA
      CHARACTER*4 HMLTN
C
      COMPLEX*16 C(MDM,MDM),OLAPLL(MDM,MDM),OLAPSS(MDM,MDM)
C
      DIMENSION FRC(MDM,MCT,MKP,(MKP+1)/2,3)
      DIMENSION BRD(MDM,3),TOT(MDM,3),RHO(MCT,3)
C
      COMMON/ATOM/ELMNT
      COMMON/COEF/C
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/QNMS/ICNBAS(MDM),KQNBAS(MDM),MQNBAS(MDM)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     INITIALISE COUNTER MATRICES
      DO IOCC=1,NOCC+IVIR
        DO ICNT=1,NCNT
          DO IKAP=1,NKAP(ICNT)
            IKQN = KVALS(IKAP,ICNT)
            NMV  = IABS(IKQN)
            DO IMV=1,NMV
              DO IT=1,3
                FRC(IOCC,ICNT,IKAP,IMV,IT) = 0.0D0
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO IT=1,3
          BRD(IOCC,IT) = 0.0D0
          TOT(IOCC,IT) = 0.0D0
        ENDDO
      ENDDO
      DO ICNT=1,NCNT
        DO IT=1,3
          RHO(ICNT,IT) = 0.0D0
        ENDDO
      ENDDO
C
C     GENERATE DIRECT OVERLAP MATRICES
      CALL VMNPOLE(OLAPLL,1,0,1,2)
      CALL VMNPOLE(OLAPSS,4,0,1,2)
C
C     LOOP OVER OCCUPIED ORBITALS
      DO IOCC=1,NOCC+IVIR
C
C       IGNORE NEGATIVE SPECTRUM
        MOCC = IOCC + NSHIFT
C
C       LOOP OVER FOCK MATRIX ADDRESS BLOCKS
        DO I=1,NDIM-NSHIFT
          IS = I+NSHIFT
          ICNT = ICNBAS(I)
          IKQN = KQNBAS(I)
          IMQN = IABS(MQNBAS(I))
          IF(IKQN.LT.0) THEN
            IKAP =-2*IKQN-1
          ELSE
            IKAP = 2*IKQN
          ENDIF
          IMV = (IMQN+1)/2
          DO J=1,NDIM-NSHIFT
            JS = J+NSHIFT
            JCNT = ICNBAS(J)
            JKQN = KQNBAS(J)
            JMQN = IABS(MQNBAS(J))
C
C           LARGE AND SMALL CONTRIBUTIONS
            EL = DREAL(DCONJG(C(I ,MOCC))*C(J ,MOCC)*OLAPLL(I,J))
            ES = DREAL(DCONJG(C(IS,MOCC))*C(JS,MOCC)*OLAPSS(I,J))
C
C           UPDATE CHARGE ON CENTER ICNT (OCCUPIED ORBITALS ONLY)
            IF(IOCC.LE.NOCC) THEN
              RHO(ICNT,1) = RHO(ICNT,1) + EL
              RHO(ICNT,2) = RHO(ICNT,2) + ES
              RHO(ICNT,3) = RHO(ICNT,3) + EL + ES
            ENDIF

C           DECIDE WHERE TO PUT CONTRIBUTION
            IF(ICNT.EQ.JCNT.AND.IKQN.EQ.JKQN.AND.IMQN.EQ.JMQN) THEN
              FRC(IOCC,ICNT,IKAP,IMV,1) = FRC(IOCC,ICNT,IKAP,IMV,1) + EL
              FRC(IOCC,ICNT,IKAP,IMV,2) = FRC(IOCC,ICNT,IKAP,IMV,2) + ES
              FRC(IOCC,ICNT,IKAP,IMV,3) = FRC(IOCC,ICNT,IKAP,IMV,1) 
     &                                  + FRC(IOCC,ICNT,IKAP,IMV,2)
            ELSE
              BRD(IOCC,1) = BRD(IOCC,1) + EL
              BRD(IOCC,2) = BRD(IOCC,2) + ES
              BRD(IOCC,3) = BRD(IOCC,3) + EL + ES
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
                TOT(IOCC,1) = TOT(IOCC,1) + FRC(IOCC,ICNT,IKAP,IMV,1)
                TOT(IOCC,2) = TOT(IOCC,2) + FRC(IOCC,ICNT,IKAP,IMV,2)
                TOT(IOCC,3) = TOT(IOCC,3) + FRC(IOCC,ICNT,IKAP,IMV,3)
              ENDDO
            ENDDO
        ENDDO
        TOT(IOCC,1) = TOT(IOCC,1) + BRD(IOCC,1)
        TOT(IOCC,2) = TOT(IOCC,2) + BRD(IOCC,2)
        TOT(IOCC,3) = TOT(IOCC,3) + BRD(IOCC,3)
C
C     END LOOP OVER OCCUPIED ORBITALS
      ENDDO
C
C     RESULTS: CHARGES ON EACH CENTER
20    FORMAT(1X,'Total charge on center ',I2,' = ',F15.10)
21    FORMAT(1X,'Total charge on molecule  = ',F15.10)
      WRITE(6, *) 'Mulliken population analysis:'
      WRITE(7, *) 'Mulliken population analysis:'
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
      SUM = 0.0D0
      DO ICNT=1,NCNT
        WRITE(6,20) ICNT,ZNUC(ICNT)-RHO(ICNT,3)
        WRITE(7,20) ICNT,ZNUC(ICNT)-RHO(ICNT,3)
        SUM = SUM + ZNUC(ICNT)-RHO(ICNT,3)
      ENDDO
      WRITE(6,21) SUM
      WRITE(7,21) SUM
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     RESULTS: DIRAC BASIS DECOMPOSITION FOR EACH ORBITAL
22    FORMAT(' Orb.',2X,'Cent.',2X,'KQN',3X,'|MQN|',10X,
     &                                   'Q(L)',13X,'Q(S)',11X,'Q(TOT)')
23    FORMAT(1X,I3,2X,I2,'(',A,')',3X,I2,3X,I2,'/2',4X,F11.8,6X,
     &                                                   F11.8,6X,F11.8)
24    FORMAT(1X,I3,3X,I2,'(',A,')',2X,I2,3X,I2,'/2',4X,F11.8,6X,
     &                                                   F11.8,6X,F11.8)
25    FORMAT(1X,I3,3X,A ,15X,F11.8,6X,F11.8,6X,F11.8)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6,22)
      WRITE(7,22)
      DO IOCC=1,NOCC+IVIR
        WRITE(6, *) REPEAT('-',72)
        WRITE(7, *) REPEAT('-',72)
        IF(IOCC.EQ.NOCC+1) THEN
          WRITE(6, *) 'Virtual orbitals (not actually occupied):'
          WRITE(7, *) 'Virtual orbitals (not actually occupied):'
          WRITE(6, *) REPEAT('-',72)
          WRITE(7, *) REPEAT('-',72)        
        ENDIF
        DO ICNT=1,NCNT
          ELA = ELMNT(IZNUC(ICNT))
          DO IKAP=1,NKAP(ICNT)
            IKQN = KVALS(IKAP,ICNT)
            NMV  = IABS(IKQN)
            DO IMV=1,NMV
              IF(FRC(IOCC,ICNT,IKAP,IMV,3).GT.1.0D-9) THEN
                IF(NCNT.LT.10) THEN
                  WRITE(6,23) IOCC,ICNT,ELA,KVALS(IKAP,ICNT),2*IMV-1,
     &                               (FRC(IOCC,ICNT,IKAP,IMV,IT),IT=1,3)
                  WRITE(7,23) IOCC,ICNT,ELA,KVALS(IKAP,ICNT),2*IMV-1,
     &                               (FRC(IOCC,ICNT,IKAP,IMV,IT),IT=1,3)
                ELSE
                  WRITE(6,24) IOCC,ICNT,ELA,KVALS(IKAP,ICNT),2*IMV-1,
     &                               (FRC(IOCC,ICNT,IKAP,IMV,IT),IT=1,3)
                  WRITE(7,24) IOCC,ICNT,ELA,KVALS(IKAP,ICNT),2*IMV-1,
     &                               (FRC(IOCC,ICNT,IKAP,IMV,IT),IT=1,3)         
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        WRITE(6,25) IOCC,'Polar.',(BRD(IOCC,IT),IT=1,3)
        WRITE(7,25) IOCC,'Polar.',(BRD(IOCC,IT),IT=1,3)
        WRITE(6, *) REPEAT('-',72)
        WRITE(7, *) REPEAT('-',72)
        WRITE(6,25) IOCC,'Total ',(TOT(IOCC,IT),IT=1,3)
        WRITE(7,25) IOCC,'Total ',(TOT(IOCC,IT),IT=1,3)        
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ELCMNPL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   EEEEEEEE LL       CCCCCC  MM       MM NN    NN PPPPPPP  LL         C
C   EE       LL      CC    CC MMM     MMM NNN   NN PP    PP LL         C
C   EE       LL      CC       MMMM   MMMM NNNN  NN PP    PP LL         C
C   EEEEEE   LL      CC       MM MM MM MM NN NN NN PP    PP LL         C
C   EE       LL      CC       MM  MMM  MM NN  NNNN PPPPPPP  LL         C
C   EE       LL      CC    CC MM   M   MM NN   NNN PP       LL         C
C   EEEEEEEE LLLLLLLL CCCCCC  MM       MM NN    NN PP       LLLLLLLL   C
C                                                                      C
C -------------------------------------------------------------------- C
C  ELCMNPL PERFORMS A DIRECT OVERLAP ANALYSIS.                         C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 E1
      COMPLEX*16 V1(MDM,MDM),V2(MDM,MDM),E2(MDM),E3(MDM)
      COMPLEX*16 OLAPLL(MDM,MDM),OLAPSS(MDM,MDM),EMPTY(MDM,MDM)
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     INITIALISE DUMMY ARRAY
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          EMPTY(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     GENERATE DIRECT OVERLAP MATRIX
      CALL VMNPOLE(OLAPLL,1,0,1,2)
      CALL VMNPOLE(OLAPSS,4,0,1,2)
C
C     CALCULATE MOLECULAR EXPECTATION VALUE OVER THE DENSITY MATRIX
      CALL PROPRTY(E1,OLAPLL,OLAPSS,EMPTY,EMPTY)
C
      WRITE(6, *) 'Direct overlap analysis:'
      WRITE(7, *) 'Direct overlap analysis:'
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
C
C     WRITE OCCUPATION NUMBER EXPECTATION VALUE
20    FORMAT(1X,A,F15.10)
      WRITE(6,20) 'Occupation number = ',DREAL(E1)
      WRITE(7,20) 'Occupation number = ',DREAL(E1)
C
C     GENERATE ARRAY OF FIRST ORDER MATRIX ELEMENTS
      CALL RS1(V1,E2,OLAPLL,OLAPSS,EMPTY,EMPTY)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_ELCMNPL_1.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8, *) E2(IOCC),(V1(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
C     SECOND ORDER INTERACTION ELEMENTS
      CALL RS2(V1,V2,E3)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_ELCMNPL_2.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E3(IOCC),(V2(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
      RETURN
      END
C
C
      SUBROUTINE ELCDIPL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      EEEEEEEE LL       CCCCCC  DDDDDDD  IIII PPPPPPP  LL             C
C      EE       LL      CC    CC DD    DD  II  PP    PP LL             C
C      EE       LL      CC       DD    DD  II  PP    PP LL             C
C      EEEEEE   LL      CC       DD    DD  II  PP    PP LL             C
C      EE       LL      CC       DD    DD  II  PPPPPPP  LL             C
C      EE       LL      CC    CC DD    DD  II  PP       LL             C
C      EEEEEEEE LLLLLLLL CCCCCC  DDDDDDD  IIII PP       LLLLLLLL       C
C                                                                      C
C -------------------------------------------------------------------- C
C  ELCDIPL PERFORMS A MOLECULAR ELECTRIC DIPOLE ANALYSIS.              C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 EX1,EY1,EZ1
      COMPLEX*16 VX1(MDM,MDM),VY1(MDM,MDM),VZ1(MDM,MDM),
     &           VX2(MDM,MDM),VY2(MDM,MDM),VZ2(MDM,MDM),
     &           EX2(MDM),EY2(MDM),EZ2(MDM),EX3(MDM),EY3(MDM),EZ3(MDM)
      COMPLEX*16 DXLL(MDM,MDM),DYLL(MDM,MDM),DZLL(MDM,MDM),
     &           DXSS(MDM,MDM),DYSS(MDM,MDM),DZSS(MDM,MDM),
     &           EMPTY(MDM,MDM)
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     INITIALISE DUMMY ARRAY
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          EMPTY(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     GENERATE ELECTRIC DIPOLE BASIS FUNCTION OVERLAPS
      CALL VDIPOLE(DXLL,1,1,1,1,2)
      CALL VDIPOLE(DYLL,1,2,2,1,2)
      CALL VDIPOLE(DZLL,1,3,3,1,2)
C
      CALL VDIPOLE(DXSS,4,1,1,1,2)
      CALL VDIPOLE(DYSS,4,2,2,1,2)
      CALL VDIPOLE(DZSS,4,3,3,1,2)
C
C     CALCULATE MOLECULAR EXPECTATION VALUE OVER THE DENSITY MATRIX
      CALL PROPRTY(EX1,DXLL,DXSS,EMPTY,EMPTY)
      CALL PROPRTY(EY1,DYLL,DYSS,EMPTY,EMPTY)
      CALL PROPRTY(EZ1,DZLL,DZSS,EMPTY,EMPTY)
C
      RMUX = 0.5D0*(ZNUC(2)-ZNUC(1))*(COORD(1,2)-COORD(1,1))-DREAL(EX1)
      RMUY = 0.5D0*(ZNUC(2)-ZNUC(1))*(COORD(2,2)-COORD(2,1))-DREAL(EY1)
      RMUZ = 0.5D0*(ZNUC(2)-ZNUC(1))*(COORD(3,2)-COORD(3,1))-DREAL(EZ1)
      RMUX = RMUX*0.529D0*4.8D0
      RMUY = RMUY*0.529D0*4.8D0
      RMUZ = RMUZ*0.529D0*4.8D0
C
      WRITE(6, *) 'Electric dipole moment overlap analysis:'
      WRITE(7, *) 'Electric dipole moment overlap analysis:'
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
C
20    FORMAT(1X,A,3X,'(',F13.8,',',F13.8,',',F13.8,')')
      WRITE(6,20) 'Dipole moment = ',DREAL(EX1),DREAL(EY1),DREAL(EZ1)
      WRITE(7,20) 'Dipole moment = ',DREAL(EX1),DREAL(EY1),DREAL(EZ1)
      WRITE(6,20) 'Debye units   = ',RMUX,RMUY,RMUZ
      WRITE(7,20) 'Debye units   = ',RMUX,RMUY,RMUZ
C
C     GENERATE ARRAY OF FIRST ORDER MATRIX ELEMENTS
      CALL RS1(VX1,EX2,DXLL,DXSS,EMPTY,EMPTY)
      CALL RS1(VY1,EY2,DYLL,DYSS,EMPTY,EMPTY)
      CALL RS1(VZ1,EZ2,DZLL,DZSS,EMPTY,EMPTY)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=11,FILE=TRIM(OUTFL)//'_ELCDIPL_X1.dat',STATUS='UNKNOWN')
      OPEN(UNIT=12,FILE=TRIM(OUTFL)//'_ELCDIPL_Y1.dat',STATUS='UNKNOWN')
      OPEN(UNIT=13,FILE=TRIM(OUTFL)//'_ELCDIPL_Z1.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(11,*) EX2(IOCC),(VX1(IOCC,JOCC),JOCC=1,NDIM)
          WRITE(12,*) EY2(IOCC),(VY1(IOCC,JOCC),JOCC=1,NDIM)
          WRITE(13,*) EZ2(IOCC),(VZ1(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=13)
      CLOSE(UNIT=12)
      CLOSE(UNIT=11)
C
C     SECOND ORDER INTERACTION ELEMENTS
      CALL RS2(VX1,VX2,EX3)
      CALL RS2(VY1,VY2,EY3)
      CALL RS2(VZ1,VZ2,EZ3)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=11,FILE=TRIM(OUTFL)//'_ELCDIPL_X2.dat',STATUS='UNKNOWN')
      OPEN(UNIT=12,FILE=TRIM(OUTFL)//'_ELCDIPL_Y2.dat',STATUS='UNKNOWN')
      OPEN(UNIT=13,FILE=TRIM(OUTFL)//'_ELCDIPL_Z2.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(11,*) EX3(IOCC),(VX2(IOCC,JOCC),JOCC=1,NDIM)
          WRITE(12,*) EY3(IOCC),(VY2(IOCC,JOCC),JOCC=1,NDIM)
          WRITE(13,*) EZ3(IOCC),(VZ2(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=13)
      CLOSE(UNIT=12)
      CLOSE(UNIT=11)
C
      RETURN
      END
C
C
      SUBROUTINE ELCQDPL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    EEEEEEEE LL       CCCCCC   QQQQQQ   DDDDDDD  PPPPPPP  LL          C
C    EE       LL      CC    CC QQ    QQ  DD    DD PP    PP LL          C
C    EE       LL      CC       QQ    QQ  DD    DD PP    PP LL          C
C    EEEEEE   LL      CC       QQ    QQ  DD    DD PP    PP LL          C
C    EE       LL      CC       QQ   QQQ  DD    DD PPPPPPP  LL          C
C    EE       LL      CC    CC QQ    QQ  DD    DD PP       LL          C
C    EEEEEEEE LLLLLLLL CCCCCC   QQQQQQ Q DDDDDDD  PP       LLLLLLLL    C
C                                                                      C
C -------------------------------------------------------------------- C
C  ELCQDPL PERFORMS A MOLECULAR ELECTRIC QUADRUPOLE ANALYSIS.          C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM),C(MDM,MDM)
      COMPLEX*16 QXXLL(MDM,MDM),QYYLL(MDM,MDM),QZZLL(MDM,MDM),
     &           QXYLL(MDM,MDM),QYZLL(MDM,MDM),QZXLL(MDM,MDM),
     &           QXXSS(MDM,MDM),QYYSS(MDM,MDM),QZZSS(MDM,MDM),
     &           QXYSS(MDM,MDM),QYZSS(MDM,MDM),QZXSS(MDM,MDM)
C
      COMMON/COEF/C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     GENERATE ELECTRIC QUADRUPOLE OVERLAPS
      CALL VQDPOLE(QXXLL,1,0,1,1,1,2)
      CALL VQDPOLE(QYYLL,1,0,2,2,1,2)
      CALL VQDPOLE(QZZLL,1,0,3,3,1,2)
      CALL VQDPOLE(QXYLL,1,0,1,2,1,2)
      CALL VQDPOLE(QYZLL,1,0,2,3,1,2)
      CALL VQDPOLE(QZXLL,1,0,3,1,1,2)
C
      CALL VQDPOLE(QXXSS,4,0,1,1,1,2)
      CALL VQDPOLE(QYYSS,4,0,2,2,1,2)
      CALL VQDPOLE(QZZSS,4,0,3,3,1,2)
      CALL VQDPOLE(QXYSS,4,0,1,2,1,2)
      CALL VQDPOLE(QYZSS,4,0,2,3,1,2)
      CALL VQDPOLE(QZXSS,4,0,3,1,1,2)
C
      RETURN
      END
C
C
      SUBROUTINE MAGDIPL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     MM       MM    AA     GGGGGG  DDDDDDD  IIII PPPPPPP  LL          C
C     MMM     MMM   AAAA   GG    GG DD    DD  II  PP    PP LL          C
C     MMMM   MMMM  AA  AA  GG       DD    DD  II  PP    PP LL          C
C     MM MM MM MM AA    AA GG       DD    DD  II  PP    PP LL          C
C     MM  MMM  MM AAAAAAAA GG   GGG DD    DD  II  PPPPPPP  LL          C
C     MM   M   MM AA    AA GG    GG DD    DD  II  PP       LL          C
C     MM       MM AA    AA  GGGGGG  DDDDDDD  IIII PP       LLLLLLLL    C
C                                                                      C
C -------------------------------------------------------------------- C
C  MAGDIPL PERFORMS A MOLECULAR MAGNETIC DIPOLE ANALYSIS.              C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM),C(MDM,MDM)
      COMPLEX*16 DPXLL(MDM,MDM),DPYLL(MDM,MDM),DPZLL(MDM,MDM),
     &           DPXSS(MDM,MDM),DPYSS(MDM,MDM),DPZSS(MDM,MDM)
C
      COMMON/COEF/C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     GENERATE MAGNETIC DIPOLE OVERLAPS
      CALL VDIPOLE(DPXLL,1,1,2,1,2)
      CALL VDIPOLE(DPYLL,1,2,3,1,2)
      CALL VDIPOLE(DPZLL,1,3,1,1,2)
C
      CALL VDIPOLE(DPXLL,4,1,2,1,2)
      CALL VDIPOLE(DPYLL,4,2,3,1,2)
      CALL VDIPOLE(DPZLL,4,3,1,1,2)
C
      RETURN
      END
C
C
      SUBROUTINE MAGQDPL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  MM       MM    AA     GGGGGG   QQQQQQ   DDDDDDD  PPPPPPP  LL        C
C  MMM     MMM   AAAA   GG    GG QQ    QQ  DD    DD PP    PP LL        C
C  MMMM   MMMM  AA  AA  GG       QQ    QQ  DD    DD PP    PP LL        C
C  MM MM MM MM AA    AA GG       QQ    QQ  DD    DD PP    PP LL        C
C  MM  MMM  MM AAAAAAAA GG   GGG QQ   QQQ  DD    DD PPPPPPP  LL        C
C  MM   M   MM AA    AA GG    GG QQ    QQ  DD    DD PP       LL        C
C  MM       MM AA    AA  GGGGGG   QQQQQQ Q DDDDDDD  PP       LLLLLLLL  C
C                                                                      C
C -------------------------------------------------------------------- C
C  MAGQDPL PERFORMS A MOLECULAR MAGNETIC QUADRUPOLE ANALYSIS.          C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM),C(MDM,MDM)
      COMPLEX*16 QXXLL(MDM,MDM),QYYLL(MDM,MDM),QZZLL(MDM,MDM),
     &           QXYLL(MDM,MDM),QYZLL(MDM,MDM),QZXLL(MDM,MDM),
     &           QXXSS(MDM,MDM),QYYSS(MDM,MDM),QZZSS(MDM,MDM),
     &           QXYSS(MDM,MDM),QYZSS(MDM,MDM),QZXSS(MDM,MDM)
C
      COMMON/COEF/C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     GENERATE MAGNETIC QUADRUPOLE OVERLAPS
      CALL VQDPOLE(QXXLL,1,0,1,1,1,2)
      CALL VQDPOLE(QYYLL,1,0,2,2,1,2)
      CALL VQDPOLE(QZZLL,1,0,3,3,1,2)
      CALL VQDPOLE(QXYLL,1,0,1,2,1,2)
      CALL VQDPOLE(QYZLL,1,0,2,3,1,2)
      CALL VQDPOLE(QZXLL,1,0,3,1,1,2)
C
      CALL VQDPOLE(QXXSS,4,0,1,1,1,2)
      CALL VQDPOLE(QYYSS,4,0,2,2,1,2)
      CALL VQDPOLE(QZZSS,4,0,3,3,1,2)
      CALL VQDPOLE(QXYSS,4,0,1,2,1,2)
      CALL VQDPOLE(QYZSS,4,0,2,3,1,2)
      CALL VQDPOLE(QZXSS,4,0,3,1,1,2)
C
      RETURN
      END
C
C
      SUBROUTINE STRKEFF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      SSSSSS TTTTTTTT RRRRRRR  KK    KK EEEEEEEE FFFFFFFF FFFFFFFF    C
C     SS    SS   TT    RR    RR KK   KK  EE       FF       FF          C
C     SS         TT    RR    RR KK  KK   EE       FF       FF          C
C      SSSSSS    TT    RR    RR KKKKK    EEEEEE   FFFFFF   FFFFFF      C
C           SS   TT    RRRRRRR  KK  KK   EE       FF       FF          C
C     SS    SS   TT    RR    RR KK   KK  EE       FF       FF          C
C      SSSSSS    TT    RR    RR KK    KK EEEEEEEE FF       FF          C
C                                                                      C
C -------------------------------------------------------------------- C
C  STRKEFF PERFORMS A STARK EFFECT ANALYSIS, GIVEN APPLIED ELECTRIC    C
C  FIELD (EX,EY,EZ).                                                   C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 E1
      COMPLEX*16 V1(MDM,MDM),V2(MDM,MDM),E2(MDM),E3(MDM)
      COMPLEX*16 OLAPLL(MDM,MDM),OLAPSS(MDM,MDM),EMPTY(MDM,MDM)
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     INITIALISE DUMMY ARRAY
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          EMPTY(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     GENERATE DIRECT OVERLAP MATRIX
      CALL VMNPOLE(OLAPLL,1,0,1,2)
      CALL VMNPOLE(OLAPSS,4,0,1,2)
C
C     CALCULATE MOLECULAR EXPECTATION VALUE OVER THE DENSITY MATRIX
      CALL PROPRTY(E1,OLAPLL,OLAPSS,EMPTY,EMPTY)
C
      WRITE(6, *) 'Stark effect analysis:'
      WRITE(7, *) 'Stark effect analysis:'
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
C
C     WRITE OCCUPATION NUMBER EXPECTATION VALUE
20    FORMAT(1X,A,F15.10)
      WRITE(6,20) 'Molecular property = ',DREAL(E1)
      WRITE(7,20) 'Molecular property = ',DREAL(E1)
C
C     GENERATE ARRAY OF FIRST ORDER MATRIX ELEMENTS
      CALL RS1(V1,E2,OLAPLL,OLAPSS,EMPTY,EMPTY)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_STRKEFF_1.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E2(IOCC),(V1(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
C     SECOND ORDER INTERACTION ELEMENTS
      CALL RS2(V1,V2,E3)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_STRKEFF_2.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E3(IOCC),(V2(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
      RETURN
      END
C
C
      SUBROUTINE ZMANEFF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   ZZZZZZZZ MM       MM    AA    NN    NN EEEEEEEE FFFFFFFF FFFFFFFF  C
C        ZZ  MMM     MMM   AAAA   NNN   NN EE       FF       FF        C
C       ZZ   MMMM   MMMM  AA  AA  NNNN  NN EE       FF       FF        C
C      ZZ    MM MM MM MM AA    AA NN NN NN EEEEEE   FFFFFF   FFFFFF    C
C     ZZ     MM  MMM  MM AAAAAAAA NN  NNNN EE       FF       FF        C
C    ZZ      MM   M   MM AA    AA NN   NNN EE       FF       FF        C
C   ZZZZZZZZ MM       MM AA    AA NN    NN EEEEEEEE FF       FF        C
C                                                                      C
C -------------------------------------------------------------------- C
C  ZMANEFF PERFORMS A ZEEMAN EFFECT ANALYSIS, GIVEN APPLIED MAGNETIC   C
C  FIELD (BX,BY,BZ).                                                   C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 E1
      COMPLEX*16 V1(MDM,MDM),V2(MDM,MDM),E2(MDM),E3(MDM)
      COMPLEX*16 OLAPLL(MDM,MDM),OLAPSS(MDM,MDM),EMPTY(MDM,MDM)
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     INITIALISE DUMMY ARRAY
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          EMPTY(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     GENERATE DIRECT OVERLAP MATRIX
      CALL VMNPOLE(OLAPLL,1,0,1,2)
      CALL VMNPOLE(OLAPSS,4,0,1,2)
C
C     CALCULATE MOLECULAR EXPECTATION VALUE OVER THE DENSITY MATRIX
      CALL PROPRTY(E1,OLAPLL,OLAPSS,EMPTY,EMPTY)
C
      WRITE(6, *) 'Zeeman effect analysis:'
      WRITE(7, *) 'Zeeman effect analysis:'
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
C
C     WRITE OCCUPATION NUMBER EXPECTATION VALUE
20    FORMAT(1X,A,F15.10)
      WRITE(6,20) 'Molecular property = ',DREAL(E1)
      WRITE(7,20) 'Molecular property = ',DREAL(E1)
C
C     GENERATE ARRAY OF FIRST ORDER MATRIX ELEMENTS
      CALL RS1(V1,E2,OLAPLL,OLAPSS,EMPTY,EMPTY)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_ZMANEFF_1.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E2(IOCC),(V1(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
C     SECOND ORDER INTERACTION ELEMENTS
      CALL RS2(V1,V2,E3)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_ZMANEFF_2.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E3(IOCC),(V2(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
      RETURN
      END
C
C
      SUBROUTINE GTENSOR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     GGGGGG TTTTTTTT EEEEEEEE NN    NN  SSSSSS   OOOOOO  RRRRRRR      C
C    GG    GG   TT    EE       NNN   NN SS    SS OO    OO RR    RR     C
C    GG         TT    EE       NNNN  NN SS       OO    OO RR    RR     C
C    GG         TT    EEEEEE   NN NN NN  SSSSSS  OO    OO RR    RR     C
C    GG   GGG   TT    EE       NN  NNNN       SS OO    OO RRRRRRR      C
C    GG    GG   TT    EE       NN   NNN SS    SS OO    OO RR    RR     C
C     GGGGGG    TT    EEEEEEEE NN    NN  SSSSSS   OOOOOO  RR    RR     C
C                                                                      C
C -------------------------------------------------------------------- C
C  GTENSOR CALCULATES THE MAGNETIC G-TENSOR OF A MOLECULE.             C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 E1
      COMPLEX*16 V1(MDM,MDM),V2(MDM,MDM),E2(MDM),E3(MDM)
      COMPLEX*16 OLAPLL(MDM,MDM),OLAPSS(MDM,MDM),EMPTY(MDM,MDM)
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     INITIALISE DUMMY ARRAY
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          EMPTY(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     GENERATE DIRECT OVERLAP MATRIX
      CALL VMNPOLE(OLAPLL,1,0,1,2)
      CALL VMNPOLE(OLAPSS,4,0,1,2)
C
C     CALCULATE MOLECULAR EXPECTATION VALUE OVER THE DENSITY MATRIX
      CALL PROPRTY(E1,OLAPLL,OLAPSS,EMPTY,EMPTY)
C
      WRITE(6, *) 'G-tensor analysis:'
      WRITE(7, *) 'G-tensor analysis:'
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
C
C     WRITE OCCUPATION NUMBER EXPECTATION VALUE
20    FORMAT(1X,A,F15.10)
      WRITE(6,20) 'Molecular property = ',DREAL(E1)
      WRITE(7,20) 'Molecular property = ',DREAL(E1)
C
C     GENERATE ARRAY OF FIRST ORDER MATRIX ELEMENTS
      CALL RS1(V1,E2,OLAPLL,OLAPSS,EMPTY,EMPTY)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_GTENSOR_1.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E2(IOCC),(V1(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
C     SECOND ORDER INTERACTION ELEMENTS
      CALL RS2(V1,V2,E3)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_GTENSOR_2.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E3(IOCC),(V2(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
      RETURN
      END
C
C
      SUBROUTINE HYPFINE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      HH    HH YY    YY PPPPPPP  FFFFFFFF IIII NN    NN EEEEEEEE      C
C      HH    HH YY    YY PP    PP FF        II  NNN   NN EE            C
C      HH    HH YY    YY PP    PP FF        II  NNNN  NN EE            C
C      HHHHHHHH  YY  YY  PP    PP FFFFFF    II  NN NN NN EEEEEE        C
C      HH    HH   YYYY   PPPPPPP  FF        II  NN  NNNN EE            C
C      HH    HH    YY    PP       FF        II  NN   NNN EE            C
C      HH    HH    YY    PP       FF       IIII NN    NN EEEEEEEE      C
C                                                                      C
C -------------------------------------------------------------------- C
C  HYPFINE PERFORMS A HYPERFINE EFFECT ANALYSIS, GIVEN A MAGNETIC      C
C  NUCLEAR MOMENT.                                                     C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 E1
      COMPLEX*16 V1(MDM,MDM),V2(MDM,MDM),E2(MDM),E3(MDM)
      COMPLEX*16 OLAPLL(MDM,MDM),OLAPSS(MDM,MDM),EMPTY(MDM,MDM)
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     INITIALISE DUMMY ARRAY
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          EMPTY(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     GENERATE DIRECT OVERLAP MATRIX
      CALL VMNPOLE(OLAPLL,1,0,1,2)
      CALL VMNPOLE(OLAPSS,4,0,1,2)
C
C     CALCULATE MOLECULAR EXPECTATION VALUE OVER THE DENSITY MATRIX
      CALL PROPRTY(E1,OLAPLL,OLAPSS,EMPTY,EMPTY)
C
      WRITE(6, *) 'Hyperfine interaction analysis:'
      WRITE(7, *) 'Hyperfine interaction analysis:'
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
C
C     WRITE OCCUPATION NUMBER EXPECTATION VALUE
20    FORMAT(1X,A,F15.10)
      WRITE(6,20) 'Molecular property = ',DREAL(E1)
      WRITE(7,20) 'Molecular property = ',DREAL(E1)
C
C     GENERATE ARRAY OF FIRST ORDER MATRIX ELEMENTS
      CALL RS1(V1,E2,OLAPLL,OLAPSS,EMPTY,EMPTY)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_HYPFINE_1.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E2(IOCC),(V1(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
C     SECOND ORDER INTERACTION ELEMENTS
      CALL RS2(V1,V2,E3)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_HYPFINE_2.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E3(IOCC),(V2(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
      RETURN
      END
C
C
      SUBROUTINE EEDMSML
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C EEEEEEEE EEEEEEEE DDDDDDD  MM       MM  SSSSSS  MM       MM LL       C
C EE       EE       DD    DD MMM     MMM SS    SS MMM     MMM LL       C
C EE       EE       DD    DD MMMM   MMMM SS       MMMM   MMMM LL       C
C EEEEEE   EEEEEE   DD    DD MM MM MM MM  SSSSSS  MM MM MM MM LL       C
C EE       EE       DD    DD MM  MMM  MM       SS MM  MMM  MM LL       C
C EE       EE       DD    DD MM   M   MM SS    SS MM   M   MM LL       C
C EEEEEEEE EEEEEEEE DDDDDDD  MM       MM  SSSSSS  MM       MM LLLLLLLL C
C                                                                      C
C -------------------------------------------------------------------- C
C  EEDMSML PERFORMS AN ATOM-CENTERED PT-ODD ELECTRON EDM ANALYSIS,     C
C  USING THE SMALL-SMALL OVERLAP AND ELECTRIC FIELD OPERATOR.          C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 E1
      COMPLEX*16 V1(MDM,MDM),V2(MDM,MDM),E2(MDM),E3(MDM)
      COMPLEX*16 OLAPLL(MDM,MDM),OLAPSS(MDM,MDM),EMPTY(MDM,MDM)
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     INITIALISE DUMMY ARRAY
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          EMPTY(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     GENERATE DIRECT OVERLAP MATRIX
      CALL VMNPOLE(OLAPLL,1,0,1,2)
      CALL VMNPOLE(OLAPSS,4,0,1,2)
C
C     CALCULATE MOLECULAR EXPECTATION VALUE OVER THE DENSITY MATRIX
      CALL PROPRTY(E1,OLAPLL,OLAPSS,EMPTY,EMPTY)
C
      WRITE(6, *) 'Atom-centered PT-odd EDM operator analysis:'
      WRITE(7, *) 'Atom-centered PT-odd EDM operator analysis:'
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
C
C     WRITE OCCUPATION NUMBER EXPECTATION VALUE
20    FORMAT(1X,A,F15.10)
      WRITE(6,20) 'Molecular property = ',DREAL(E1)
      WRITE(7,20) 'Molecular property = ',DREAL(E1)
C
C     GENERATE ARRAY OF FIRST ORDER MATRIX ELEMENTS
      CALL RS1(V1,E2,OLAPLL,OLAPSS,EMPTY,EMPTY)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_EEDMSML_1.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E2(IOCC),(V1(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
C     SECOND ORDER INTERACTION ELEMENTS
      CALL RS2(V1,V2,E3)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_EEDMSML_2.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E3(IOCC),(V2(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
      RETURN
      END
C
C
      SUBROUTINE EEDMEFF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   EEEEEEEE EEEEEEEE DDDDDDD  MM       MM EEEEEEEE FFFFFFFF FFFFFFFF  C
C   EE       EE       DD    DD MMM     MMM EE       FF       FF        C
C   EE       EE       DD    DD MMMM   MMMM EE       FF       FF        C
C   EEEEEE   EEEEEE   DD    DD MM MM MM MM EEEEEE   FFFFFF   FFFFFF    C
C   EE       EE       DD    DD MM  MMM  MM EE       FF       FF        C
C   EE       EE       DD    DD MM   M   MM EE       FF       FF        C
C   EEEEEEEE EEEEEEEE DDDDDDD  MM       MM EEEEEEEE FF       FF        C
C                                                                      C
C -------------------------------------------------------------------- C
C  EEDMEFF PERFORMS ONE-BODY EFFECTIVE ELECTRON EDM ANALYSIS.          C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 E1
      COMPLEX*16 V1(MDM,MDM),V2(MDM,MDM),E2(MDM),E3(MDM)
      COMPLEX*16 OLAPLL(MDM,MDM),OLAPSS(MDM,MDM),EMPTY(MDM,MDM)
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     INITIALISE DUMMY ARRAY
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          EMPTY(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     GENERATE DIRECT OVERLAP MATRIX
      CALL VMNPOLE(OLAPLL,1,0,1,2)
      CALL VMNPOLE(OLAPSS,4,0,1,2)
C
C     CALCULATE MOLECULAR EXPECTATION VALUE OVER THE DENSITY MATRIX
      CALL PROPRTY(E1,OLAPLL,OLAPSS,EMPTY,EMPTY)
C
      WRITE(6, *) 'One-body effective PT-odd EDM operator analysis:'
      WRITE(7, *) 'One-body effective PT-odd EDM operator analysis:'
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
C
C     WRITE OCCUPATION NUMBER EXPECTATION VALUE
20    FORMAT(1X,A,F15.10)
      WRITE(6,20) 'Molecular property = ',DREAL(E1)
      WRITE(7,20) 'Molecular property = ',DREAL(E1)
C
C     GENERATE ARRAY OF FIRST ORDER MATRIX ELEMENTS
      CALL RS1(V1,E2,OLAPLL,OLAPSS,EMPTY,EMPTY)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_EEDMEFF_1.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E2(IOCC),(V1(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
C     SECOND ORDER INTERACTION ELEMENTS
      CALL RS2(V1,V2,E3)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_EEDMEFF_2.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E3(IOCC),(V2(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
      RETURN
      END
C
C
      SUBROUTINE SCLPTEN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      SSSSSS   CCCCCC  LL       PPPPPPP TTTTTTTT EEEEEEEE NN    NN    C
C     SS    SS CC    CC LL       PP    PP   TT    EE       NNN   NN    C
C     SS       CC       LL       PP    PP   TT    EE       NNNN  NN    C
C      SSSSSS  CC       LL       PP    PP   TT    EEEEEE   NN NN NN    C
C           SS CC       LL       PPPPPPP    TT    EE       NN  NNNN    C
C     SS    SS CC    CC LL       PP         TT    EE       NN   NNN    C
C      SSSSSS   CCCCCC  LLLLLLLL PP         TT    EEEEEEEE NN    NN    C
C                                                                      C
C -------------------------------------------------------------------- C
C  SCLPTEN SCALAR PT-ODD ELECTRON-NUCLEUS INTERACTION ANALYSIS.        C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 E1
      COMPLEX*16 V1(MDM,MDM),V2(MDM,MDM),E2(MDM),E3(MDM)
      COMPLEX*16 OLAPLL(MDM,MDM),OLAPSS(MDM,MDM),EMPTY(MDM,MDM)
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     INITIALISE DUMMY ARRAY
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          EMPTY(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     GENERATE DIRECT OVERLAP MATRIX
      CALL VMNPOLE(OLAPLL,1,0,1,2)
      CALL VMNPOLE(OLAPSS,4,0,1,2)
C
C     CALCULATE MOLECULAR EXPECTATION VALUE OVER THE DENSITY MATRIX
      CALL PROPRTY(E1,OLAPLL,OLAPSS,EMPTY,EMPTY)
C
      WRITE(6, *) 'Scalar PT-odd electron-nuclear operator analysis:'
      WRITE(7, *) 'Scalar PT-odd electron-nuclear operator analysis:'
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
C
C     WRITE OCCUPATION NUMBER EXPECTATION VALUE
20    FORMAT(1X,A,F15.10)
      WRITE(6,20) 'Molecular property = ',DREAL(E1)
      WRITE(7,20) 'Molecular property = ',DREAL(E1)
C
C     GENERATE ARRAY OF FIRST ORDER MATRIX ELEMENTS
      CALL RS1(V1,E2,OLAPLL,OLAPSS,EMPTY,EMPTY)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_SCLPTEN_1.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E2(IOCC),(V1(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
C     SECOND ORDER INTERACTION ELEMENTS
      CALL RS2(V1,V2,E3)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_SCLPTEN_2.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E3(IOCC),(V2(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
      RETURN
      END
C
C
      SUBROUTINE VECPTEN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    VV    VV EEEEEEEE  CCCCCC  PPPPPPP TTTTTTTT EEEEEEEE NN    NN     C
C    VV    VV EE       CC    CC PP    PP   TT    EE       NNN   NN     C
C    VV    VV EE       CC       PP    PP   TT    EE       NNNN  NN     C
C    VV    VV EEEEEE   CC       PP    PP   TT    EEEEEE   NN NN NN     C
C     VV  VV  EE       CC       PPPPPPP    TT    EE       NN  NNNN     C
C      VVVV   EE       CC    CC PP         TT    EE       NN   NNN     C
C       VV    EEEEEEEE  CCCCCC  PP         TT    EEEEEEEE NN    NN     C
C                                                                      C
C -------------------------------------------------------------------- C
C  VECPTEN VECTOR PT-ODD ELECTRON-NUCLEUS INTERACTION ANALYSIS.        C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 E1
      COMPLEX*16 V1(MDM,MDM),V2(MDM,MDM),E2(MDM),E3(MDM)
      COMPLEX*16 OLAPLL(MDM,MDM),OLAPSS(MDM,MDM),EMPTY(MDM,MDM)
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     INITIALISE DUMMY ARRAY
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          EMPTY(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     GENERATE DIRECT OVERLAP MATRIX
      CALL VMNPOLE(OLAPLL,1,0,1,2)
      CALL VMNPOLE(OLAPSS,4,0,1,2)
C
C     CALCULATE MOLECULAR EXPECTATION VALUE OVER THE DENSITY MATRIX
      CALL PROPRTY(E1,OLAPLL,OLAPSS,EMPTY,EMPTY)
C
      WRITE(6, *) 'Vector PT-odd electron-nuclear operator analysis:'
      WRITE(7, *) 'Vector PT-odd electron-nuclear operator analysis:'
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
C
C     WRITE OCCUPATION NUMBER EXPECTATION VALUE
20    FORMAT(1X,A,F15.10)
      WRITE(6,20) 'Molecular property = ',DREAL(E1)
      WRITE(7,20) 'Molecular property = ',DREAL(E1)
C
C     GENERATE ARRAY OF FIRST ORDER MATRIX ELEMENTS
      CALL RS1(V1,E2,OLAPLL,OLAPSS,EMPTY,EMPTY)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_VECPTEN_1.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E2(IOCC),(V1(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
C     SECOND ORDER INTERACTION ELEMENTS
      CALL RS2(V1,V2,E3)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_VECPTEN_2.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E3(IOCC),(V2(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
      RETURN
      END
C
C
      SUBROUTINE PVIOLTN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      PPPPPPP  VV    VV IIII  OOOOOO  LL      TTTTTTTT NN    NN       C
C      PP    PP VV    VV  II  OO    OO LL         TT    NNN   NN       C
C      PP    PP VV    VV  II  OO    OO LL         TT    NNNN  NN       C
C      PP    PP VV    VV  II  OO    OO LL         TT    NN NN NN       C
C      PPPPPPP   VV  VV   II  OO    OO LL         TT    NN  NNNN       C
C      PP         VVVV    II  OO    OO LL         TT    NN   NNN       C
C      PP          VV    IIII  OOOOOO  LLLLLLLL   TT    NN    NN       C
C                                                                      C
C -------------------------------------------------------------------- C
C  PVIOLTN PERFORMS A P-PODD EFFECTIVE OPERATOR ANALYSIS.              C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 E1
      COMPLEX*16 V1(MDM,MDM),V2(MDM,MDM),E2(MDM),E3(MDM)
      COMPLEX*16 OLAPLL(MDM,MDM),OLAPSS(MDM,MDM),EMPTY(MDM,MDM)
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     INITIALISE DUMMY ARRAY
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          EMPTY(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     GENERATE DIRECT OVERLAP MATRIX
      CALL VMNPOLE(OLAPLL,1,0,1,2)
      CALL VMNPOLE(OLAPSS,4,0,1,2)
C
C     CALCULATE MOLECULAR EXPECTATION VALUE OVER THE DENSITY MATRIX
      CALL PROPRTY(E1,OLAPLL,OLAPSS,EMPTY,EMPTY)
C
      WRITE(6, *) 'P-odd effective operator analysis:'
      WRITE(7, *) 'P-odd effective operator analysis:'
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
C
C     WRITE OCCUPATION NUMBER EXPECTATION VALUE
20    FORMAT(1X,A,F15.10)
      WRITE(6,20) 'Molecular property = ',DREAL(E1)
      WRITE(7,20) 'Molecular property = ',DREAL(E1)
C
C     GENERATE ARRAY OF FIRST ORDER MATRIX ELEMENTS
      CALL RS1(V1,E2,OLAPLL,OLAPSS,EMPTY,EMPTY)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_PVIOLTN_1.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E2(IOCC),(V1(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
C     SECOND ORDER INTERACTION ELEMENTS
      CALL RS2(V1,V2,E3)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_PVIOLTN_2.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E3(IOCC),(V2(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
      RETURN
      END
C
C
      SUBROUTINE BETADCY(ICNT,ZDCY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     BBBBBBB  EEEEEEEE TTTTTTTT   AA    DDDDDDD   CCCCCC  YY    YY    C
C     BB    BB EE          TT     AAAA   DD    DD CC    CC YY    YY    C
C     BB    BB EE          TT    AA  AA  DD    DD CC        YY  YY     C
C     BBBBBBB  EEEEEE      TT   AA    AA DD    DD CC         YYYY      C
C     BB    BB EE          TT   AAAAAAAA DD    DD CC          YY       C
C     BB    BB EE          TT   AA    AA DD    DD CC    CC    YY       C
C     BBBBBBB  EEEEEEEE    TT   AA    AA DDDDDDD   CCCCCC     YY       C
C                                                                      C
C -------------------------------------------------------------------- C
C  BETADCY TAKES THE RADIOACTIVE DECAY PROCESS IN WHICH ATOMIC CENTER  C
C  ICNT LOSES AN AMOUNT OF CHARGE ZDCY, AND PERFORMS AN ANALYSIS.      C
C -------------------------------------------------------------------- C
C  AT A FIRST PASS, GENERATE ONLY NUCLEAR ATTRACTION MATRIX ELEMENTS   C
C  DUE TO THE NEW CHARGE ON CENTER ICNT AND DON'T CHANGE ANY RADII.    C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      DIMENSION XYZ(3)
C
      COMPLEX*16 E1
      COMPLEX*16 V1(MDM,MDM),V2(MDM,MDM),E2(MDM),E3(MDM)
      COMPLEX*16 VNUCLL(MDM,MDM),VNUCSS(MDM,MDM),EMPTY(MDM,MDM)
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     ADJUST NUCLEAR CHARGE OF CENTER ICNT
      DO IX=1,3
        XYZ(IX) = COORD(IX,ICNT)
      ENDDO
C
C     NEW NUCLEAR CHARGE OF DECAYED CENTRE
      ZNEW = ZNUC(ICNT)-ZDCY
C
C     INITIALISE DUMMY ARRAY
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          EMPTY(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     GENERATE DIRECT OVERLAP MATRIX
C     DFNOTE: THIS ACTUALLY GENERATES MATRIX ELEMENTS FOR ALL ATOMIC
C             CENTERS IN THE MOLECULE. MIGHT BE WORTH CALCULATING 
C             MATRIX ELEMENTS OVER ICNT FROM SCRATCH HERE INSTEAD.
      CALL VNCATRC(VNUCLL,1,0,1,2)
      CALL VNCATRC(VNUCSS,4,0,1,2)
C
C     CALCULATE MOLECULAR EXPECTATION VALUE OVER THE DENSITY MATRIX
      CALL PROPRTY(E1,VNUCLL,VNUCSS,EMPTY,EMPTY)
C
      WRITE(6, *) 'Direct overlap analysis:'
      WRITE(7, *) 'Direct overlap analysis:'
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
C
C     WRITE ENERGY DIFFERENCE EXPECTATION VALUE
20    FORMAT(1X,A,F15.10)
      WRITE(6,20) 'Energy difference = ',DREAL(E1)
      WRITE(7,20) 'Energy difference = ',DREAL(E1)
C
C     GENERATE ARRAY OF FIRST ORDER MATRIX ELEMENTS
      CALL RS1(V1,E2,VNUCLL,VNUCSS,EMPTY,EMPTY)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_BETADCY_1.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E2(IOCC),(V1(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
C     SECOND ORDER INTERACTION ELEMENTS
      CALL RS2(V1,V2,E3)
C   
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_BETADCY_2.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8,*) E3(IOCC),(V2(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
      RETURN
      END
C
C
      SUBROUTINE VMNPOLE(VIJ,ITT,IQ,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  VV    VV MM       MM NN    NN PPPPPPP   OOOOOO  LL       EEEEEEEE   C
C  VV    VV MMM     MMM NNN   NN PP    PP OO    OO LL       EE         C
C  VV    VV MMMM   MMMM NNNN  NN PP    PP OO    OO LL       EE         C
C  VV    VV MM MM MM MM NN NN NN PP    PP OO    OO LL       EEEEEE     C
C   VV  VV  MM  MMM  MM NN  NNNN PPPPPPP  OO    OO LL       EE         C
C    VVVV   MM   M   MM NN   NNN PP       OO    OO LL       EE         C
C     VV    MM       MM NN    NN PP        OOOOOO  LLLLLLLL EEEEEEEE   C
C                                                                      C
C -------------------------------------------------------------------- C
C  VMNPOLE CONSTRUCTS A MATRIX OF (u,T|SIG_Q|v,T') OVERLAP INTEGRALS   C
C  OVER ALL BASIS FUNCTION PAIRS, AND SAVES THE RESULT TO VIJ.         C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    ITT   = {1,2,3,4} -> {LL,LS,SL,SS} COMPONENT COMBINATION.         C
C    IQ    = {0,1,2,3} -> {0,X,Y,Z} THE COUPLING SIGMA MATRIX.         C
C    I1,I2 = {1,2,3,4} -> BASIS FUNCTION OVERLAPS FOR EQ-COEFFS.       C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 VIJ(MDM,MDM),STR(MDM,MDM)
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 TMP(MBS,MBS,4)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
      DATA PI/3.1415926535897932D0/
C
C     INITIALISE THE OVERLAP ARRAY
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          VIJ(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     NON-RELATIVISTIC CONDITION
      IF(HMLTN.EQ.'NORL'.AND.ITT.GT.1) GOTO 100
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTERS A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     LOOP OVER CENTER A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 0
        ENDIF
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNA    = NFUNCT(LQN(1)+1,ICNTA)
        NFUNS(1) = NFUNA
C
        DO IBAS=1,NFUNA
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNB    = NFUNCT(LQN(2)+1,ICNTB)
        NFUNS(2) = NFUNB
C
        DO IBAS=1,NFUNB
          EXPT(IBAS,2) = EXPSET(IBAS,LQN(2)+1,ICNTB)
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
C     IMPLEMENT ANY AVAILABLE SELECTION RULES HERE                     C
C**********************************************************************C
C
C     SELECTION RULES TO BE MADE BASED ON GEOMETRIC SYMMETRY,
C     ATOMIC COORDINATES AND QUANTUM NUMBER PAIRS. THE IDEA IS TO
C     SKIP CALCULATIONS THAT SHOULD PRODUCE NO OR NEGLIGIBLE EFFECT.
      IF(NCNT.LE.2) THEN
        IF(MQN(1).NE.MQN(2)) GOTO 2000
      ENDIF
C
C**********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS           C
C     OF (MA,MB). FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE          C
C     ORDERED                                                          C
C     11: = (-|MQN(A)|,-|MQN(B)|)                                      C
C     12: = (-|MQN(A)|,+|MQN(B)|)                                      C
C     21: = (+|MQN(A)|,-|MQN(B)|)                                      C
C     22: = (+|MQN(A)|,+|MQN(B)|)                                      C
C**********************************************************************C
C
C     CALCULATE COMPONENT OFFSETS
      IL1 = LARGE(ICNTA,KA,MJA  )
      IL2 = LARGE(ICNTA,KA,MJA+1)
      JL1 = LARGE(ICNTB,KB,MJB  )
      JL2 = LARGE(ICNTB,KB,MJB+1)
C
      IS1 = IL1 + NSHIFT
      IS2 = IL2 + NSHIFT
      JS1 = JL1 + NSHIFT
      JS2 = JL2 + NSHIFT
C
C     CONSTRUCTION OF ONE-ELECTRON OVERLAP MATRIX BY TT' BLOCKS...
C
C     THE PHASE GENERATES E022 AND E012 COEFFS FROM E011 AND E021
      FASE = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
     &                    *DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXM = NFUNA*NFUNB
C
C**********************************************************************C
C     GENERATE OVERLAP MATRICES                                        C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)
      NTUVLL = (LAM+1)*(LAM+2)*(LAM+3)/6
C
      IALT = 1
C     GENERATE ELLQ COEFFICIENTS
      IF(ITT.EQ.1) THEN
        CALL CPU_TIME(TDM1)
        CALL EMAKELL(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,I1,I2,IQ)
        CALL CPU_TIME(TDM2)
        TELL = TELL + TDM2 - TDM1
C     GENERATE ESSQ COEFFICIENTS
      ELSEIF(ITT.EQ.4) THEN
        CALL CPU_TIME(TDM1)
        CALL EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,I1,I2,IQ)
        CALL CPU_TIME(TDM2)
        TESS = TESS + TDM2 - TDM1
C     GENERATE ELS1 COEFFICIENTS
      ELSEIF(ITT.EQ.2.OR.ITT.EQ.3) THEN
        CALL CPU_TIME(TDM1)
        CALL EMAKELS(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,I1,I2,IQ)
        CALL CPU_TIME(TDM2)
        TELS = TELS + TDM2 - TDM1
      ELSE
        WRITE(6,*) 'In OVERLAP: allowed values are ITT = {1,2,3,4}'
        WRITE(7,*) 'In OVERLAP: allowed values are ITT = {1,2,3,4}'
        WRITE(6,*) 'ITT = ',ITT
        WRITE(7,*) 'ITT = ',ITT
        STOP       
      ENDIF
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          EIJ    = EXPT(IBAS,1)+EXPT(JBAS,2)
          EROOT  = DSQRT(PI/EIJ)**3
          TMP(IBAS,JBAS,1) = EROOT*E11(M,1)
          TMP(IBAS,JBAS,3) = EROOT*E21(M,1)
          TMP(IBAS,JBAS,2) =-FASE*DCONJG(TMP(IBAS,JBAS,3))
          TMP(IBAS,JBAS,4) = FASE*DCONJG(TMP(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C**********************************************************************C
C     WE NOW HAVE ALL PIECES OF THE OVERLAP MATRIX FOR THIS BLOCK OF   C
C     BASIS FUNCTIONS -- NOW OVERLAY THE RESULTS INTO VIJ.             C
C**********************************************************************C
C
C     CALCULATE COMPONENT OFFSETS
      IL1 = LARGE(ICNTA,KA,MJA  )
      IL2 = LARGE(ICNTA,KA,MJA+1)
      JL1 = LARGE(ICNTB,KB,MJB  )
      JL2 = LARGE(ICNTB,KB,MJB+1)
C
C     ASSEMBLE TEMPORARY BLOCK INTO OVERLAP ARRAY
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            VIJ(IL1+IBAS,JL1+JBAS) = TMP(IBAS,JBAS,1)
            VIJ(IL1+IBAS,JL2+JBAS) = TMP(IBAS,JBAS,2)
            VIJ(IL2+IBAS,JL1+JBAS) = TMP(IBAS,JBAS,3)
            VIJ(IL2+IBAS,JL2+JBAS) = TMP(IBAS,JBAS,4)
            VIJ(JL1+JBAS,IL1+IBAS) = DCONJG(VIJ(IL1+IBAS,JL1+JBAS))
            VIJ(JL2+JBAS,IL1+IBAS) = DCONJG(VIJ(IL1+IBAS,JL2+JBAS))
            VIJ(JL1+JBAS,IL2+IBAS) = DCONJG(VIJ(IL2+IBAS,JL1+JBAS))
            VIJ(JL2+JBAS,IL2+IBAS) = DCONJG(VIJ(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            VIJ(IL1+IBAS,JL1+JBAS) = TMP(IBAS,JBAS,1)
            VIJ(IL1+IBAS,JL2+JBAS) = TMP(IBAS,JBAS,2)
            VIJ(IL2+IBAS,JL1+JBAS) = TMP(IBAS,JBAS,3)
            VIJ(IL2+IBAS,JL2+JBAS) = TMP(IBAS,JBAS,4)
            VIJ(JL1+JBAS,IL1+IBAS) = DCONJG(VIJ(IL1+IBAS,JL1+JBAS))
            VIJ(JL2+JBAS,IL1+IBAS) = DCONJG(VIJ(IL1+IBAS,JL2+JBAS))
            VIJ(JL1+JBAS,IL2+IBAS) = DCONJG(VIJ(IL2+IBAS,JL1+JBAS))
            VIJ(JL2+JBAS,IL2+IBAS) = DCONJG(VIJ(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
2000  CONTINUE
C
C     ESLQ[MU,NU;T,U,V] = ELSQ[NU,MU;T,U,V]*.
      IF(ITT.EQ.3) THEN
C       TRANSFER ITT=2 VIJ CASE TO A STORAGE ARRAY
        DO I=1,NDIM-NSHIFT
          DO J=1,NDIM-NSHIFT
            STR(I,J) = VIJ(I,J)
          ENDDO
        ENDDO
C       TRANSPOSE AND COMPLEX CONJUGATE
        DO I=1,NDIM-NSHIFT
          DO J=1,NDIM-NSHIFT
            VIJ(I,J) = DCONJG(STR(J,I))
          ENDDO
        ENDDO
      ENDIF
C
100   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE VDIPOLE(VIJ,ITT,IQ,IX,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      VV    VV DDDDDDD  IIII PPPPPPP   OOOOOO  LL      EEEEEEEE       C
C      VV    VV DD    DD  II  PP    PP OO    OO LL      EE             C
C      VV    VV DD    DD  II  PP    PP OO    OO LL      EE             C
C      VV    VV DD    DD  II  PP    PP OO    OO LL      EEEEEE         C
C       VV  VV  DD    DD  II  PPPPPPP  OO    OO LL      EE             C
C        VVVV   DD    DD  II  PP       OO    OO LL      EE             C
C         VV    DDDDDDD  IIII PP        OOOOOO  LLLLLLL EEEEEEEE       C
C                                                                      C
C -------------------------------------------------------------------- C
C  VDIPOLE CONSTRUCTS A MATRIX OF (u,T|SIG_Q.X|v,T') OVERLAP INTEGRALS C
C  OVER ALL BASIS FUNCTION PAIRS, AND SAVES THE RESULT TO VIJ.         C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    ITT   = {1,2,3,4} -> {LL,LS,SL,SS} COMPONENT COMBINATION.         C
C    IQ    = {0,1,2,3} -> {0,X,Y,Z} THE COUPLING SIGMA MATRIX.         C
C    IX    = {1,2,3}   -> {X,Y,Z} CARTESIAN WEIGHTING FACTOR.          C
C    I1,I2 = {1,2,3,4} -> BASIS FUNCTION OVERLAPS FOR EQ-COEFFS.       C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 CONE
      COMPLEX*16 VIJ(MDM,MDM)
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 SLL(MBS,MBS,4),SSS(MBS,MBS,4)
C
      DIMENSION RC(MB2,MRC),EXPT(MBS,4),XYZ(3,4),APH(MB2),CP(MB2,3),
     &          PNC(MB2),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
      DATA PI/3.1415926535897932D0/
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     INITIALISE THE OVERLAP ARRAY
      DO I=1,NDIM
        DO J=1,NDIM
          VIJ(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     NON-RELATIVISTIC CONDITION
      IF(HMLTN.EQ.'NORL'.AND.ITT.GT.1) GOTO 100
C
C     CENTER OF MASS COORDINATES
      CX = 0.0D0
      CY = 0.0D0
      CZ = 0.0D0
      CM = 0.0D0
      DO N=1,NCNT
        CX = CX + AMASS(N)*COORD(1,N)
        CY = CY + AMASS(N)*COORD(2,N)
        CZ = CZ + AMASS(N)*COORD(3,N)
        CM = CM + AMASS(N)
      ENDDO
      CX = CX/CM
      CY = CY/CM
      CZ = CZ/CM
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTERS A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     LOOP OVER CENTER A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 0
        ENDIF
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNA    = NFUNCT(LQN(1)+1,ICNTA)
        NFUNS(1) = NFUNA
C
        DO IBAS=1,NFUNA
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2)=KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNB    = NFUNCT(LQN(2)+1,ICNTB)
        NFUNS(2) = NFUNB
C
        DO IBAS=1,NFUNB
          EXPT(IBAS,2) = EXPSET(IBAS,LQN(2)+1,ICNTB)
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
C     IMPLEMENT ANY AVAILABLE SELECTION RULES HERE                     C
C**********************************************************************C
C
C     SELECTION RULES TO BE MADE BASED ON GEOMETRIC SYMMETRY,
C     ATOMIC COORDINATES AND QUANTUM NUMBER PAIRS. THE IDEA IS TO
C     SKIP CALCULATIONS THAT SHOULD PRODUCE NO OR NEGLIGIBLE EFFECT.
C     IF(MJA.NE.MJB) GOTO 2000
C
C**********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS           C
C     OF (MA,MB). FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE          C
C     ORDERED                                                          C
C     11: = (-|MQN(A)|,-|MQN(B)|)                                      C
C     12: = (-|MQN(A)|,+|MQN(B)|)                                      C
C     21: = (+|MQN(A)|,-|MQN(B)|)                                      C
C     22: = (+|MQN(A)|,+|MQN(B)|)                                      C
C**********************************************************************C
C
C     CONSTRUCTION OF ONE-ELECTRON OVERLAP MATRIX BY TT' BLOCKS...
C
C     THE PHASE GENERATES E022 AND E012 COEFFS FROM E011 AND E021
      FASE = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
     &                    *DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXM = NFUNA*NFUNB
C
C**********************************************************************C
C     PART 1: THE LL MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)
      NTUVLL = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ELL0 COEFFICIENTS
      IALT = 1
      CALL CPU_TIME(TDM1)
      CALL EMAKELL(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,1,2,0)
      CALL CPU_TIME(TDM2)
      TELL = TELL + TDM2 - TDM1
C
C     GAUSSIAN OVERLAP CENTERS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          EIJ = EXPT(IBAS,1)+EXPT(JBAS,2)
          PX  = (XYZ(1,1)*EXPT(IBAS,1) + XYZ(1,2)*EXPT(JBAS,2))/EIJ
          PY  = (XYZ(2,1)*EXPT(IBAS,1) + XYZ(2,2)*EXPT(JBAS,2))/EIJ
          PZ  = (XYZ(3,1)*EXPT(IBAS,1) + XYZ(3,2)*EXPT(JBAS,2))/EIJ
          CP(M,1) = PX - CX
          CP(M,2) = PY - CY
          CP(M,3) = PZ - CZ
        ENDDO
      ENDDO
C
C     CALCULATE OVERLAP MATRIX ELEMENTS
C
C *** BEGIN AN IF STATEMENT FOR LQN PAIRS
C >>> IF ONE OR MORE OF THE LQNS IS NONZERO
      IF(NTUVLL.GT.1) THEN
        IADR1 = NTUVLL-1
        IADR2 = NTUVLL
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M = M+1
            EIJ    = EXPT(IBAS,1)+EXPT(JBAS,2)
            EROOT  = DSQRT(PI/EIJ)**3
            SLL(IBAS,JBAS,1) = EROOT*(E11(M,IADR1)+CP(M,3)*E11(M,IADR2))
            SLL(IBAS,JBAS,3) = EROOT*(E21(M,IADR1)+CP(M,3)*E21(M,IADR2))
            SLL(IBAS,JBAS,2) =-FASE*DCONJG(SLL(IBAS,JBAS,3))
            SLL(IBAS,JBAS,4) = FASE*DCONJG(SLL(IBAS,JBAS,1))
          ENDDO
        ENDDO
C >>> IF BOTH LQNS ARE ZERO
      ELSE
        IADR2 = NTUVLL
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M = M+1
            EIJ    = EXPT(IBAS,1)+EXPT(JBAS,2)
            EROOT  = DSQRT(PI/EIJ)**3
            SLL(IBAS,JBAS,1) = EROOT*CP(M,3)*E11(M,IADR2)
            SLL(IBAS,JBAS,3) = EROOT*CP(M,3)*E21(M,IADR2)
            SLL(IBAS,JBAS,2) =-FASE*DCONJG(SLL(IBAS,JBAS,3))
            SLL(IBAS,JBAS,4) = FASE*DCONJG(SLL(IBAS,JBAS,1))
          ENDDO
        ENDDO     
C *** END IF STATEMENT OVER LQN PAIRS
      ENDIF
C
C**********************************************************************C
C     PART 2: THE SS MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)+2
      NTUVSS = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ESS0 COEFFICIENTS
      IALT = 1
      CALL CPU_TIME(TDM1)
      CALL EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,1,2,0)
      CALL CPU_TIME(TDM2)
      TESS = TESS + TDM2 - TDM1
C
C     GAUSSIAN OVERLAP CENTERS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          EIJ = EXPT(IBAS,1)+EXPT(JBAS,2)
          PX  = (XYZ(1,1)*EXPT(IBAS,1) + XYZ(1,2)*EXPT(JBAS,2))/EIJ
          PY  = (XYZ(2,1)*EXPT(IBAS,1) + XYZ(2,2)*EXPT(JBAS,2))/EIJ
          PZ  = (XYZ(3,1)*EXPT(IBAS,1) + XYZ(3,2)*EXPT(JBAS,2))/EIJ
          CP(M,1) = PX - CX
          CP(M,2) = PY - CY
          CP(M,3) = PZ - CZ
        ENDDO
      ENDDO
C
C     CALCULATE OVERLAP MATRIX ELEMENTS
C
C *** BEGIN AN IF STATEMENT FOR LQN PAIRS
C >>> IF ONE OR MORE OF THE LQNS IS NONZERO
      IF(NTUVSS.GT.1) THEN
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            IADR1 = NTUVSS-1
            IADR2 = NTUVSS
            M = M+1
            EIJ     = EXPT(IBAS,1)+EXPT(JBAS,2)
            EROOT   = DSQRT(PI/EIJ)**3
            SSS(IBAS,JBAS,1) = EROOT*E11(M,1)
            SSS(IBAS,JBAS,3) = EROOT*E21(M,1)
            SSS(IBAS,JBAS,2) =-FASE*DCONJG(SSS(IBAS,JBAS,3))
            SSS(IBAS,JBAS,4) = FASE*DCONJG(SSS(IBAS,JBAS,1))
          ENDDO
        ENDDO
C >>> IF BOTH LQNS ARE ZERO
      ELSE
        IADR2 = NTUVSS
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M = M+1
            EIJ    = EXPT(IBAS,1)+EXPT(JBAS,2)
            EROOT  = DSQRT(PI/EIJ)**3
            SLL(IBAS,JBAS,1) = EROOT*(E11(M,IADR1)+CP(M,3)*E11(M,IADR2))
            SLL(IBAS,JBAS,3) = EROOT*(E21(M,IADR1)+CP(M,3)*E21(M,IADR2))
            SLL(IBAS,JBAS,2) =-FASE*DCONJG(SLL(IBAS,JBAS,3))
            SLL(IBAS,JBAS,4) = FASE*DCONJG(SLL(IBAS,JBAS,1))
          ENDDO
        ENDDO     
C *** END IF STATEMENT OVER LQN PAIRS
      ENDIF
C
C**********************************************************************C
C     WE NOW HAVE ALL PIECES OF THE OVERLAP MATRIX FOR THIS BLOCK OF   C
C     BASIS FUNCTIONS -- NOW OVERLAY THE RESULTS INTO VIJ.             C
C**********************************************************************C
C
C     CALCULATE COMPONENT OFFSETS
      IL1 = LARGE(ICNTA,KA,MJA  )
      IL2 = LARGE(ICNTA,KA,MJA+1)
      JL1 = LARGE(ICNTB,KB,MJB  )
      JL2 = LARGE(ICNTB,KB,MJB+1)
C
      IS1 = IL1 + NSHIFT
      IS2 = IL2 + NSHIFT
      JS1 = JL1 + NSHIFT
      JS2 = JL2 + NSHIFT
C
C     LL BLOCKS
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            VIJ(IL1+IBAS,JL1+JBAS) = SLL(IBAS,J,1)
            VIJ(IL1+IBAS,JL2+JBAS) = SLL(IBAS,J,2)
            VIJ(IL2+IBAS,JL1+JBAS) = SLL(IBAS,J,3)
            VIJ(IL2+IBAS,JL2+JBAS) = SLL(IBAS,J,4)
            VIJ(JL1+JBAS,IL1+IBAS) = DCONJG(VIJ(IL1+IBAS,JL1+JBAS))
            VIJ(JL2+JBAS,IL1+IBAS) = DCONJG(VIJ(IL1+IBAS,JL2+JBAS))
            VIJ(JL1+JBAS,IL2+IBAS) = DCONJG(VIJ(IL2+IBAS,JL1+JBAS))
            VIJ(JL2+JBAS,IL2+IBAS) = DCONJG(VIJ(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            VIJ(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            VIJ(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            VIJ(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            VIJ(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
            VIJ(JL1+JBAS,IL1+IBAS) = DCONJG(VIJ(IL1+IBAS,JL1+JBAS))
            VIJ(JL2+JBAS,IL1+IBAS) = DCONJG(VIJ(IL1+IBAS,JL2+JBAS))
            VIJ(JL1+JBAS,IL2+IBAS) = DCONJG(VIJ(IL2+IBAS,JL1+JBAS))
            VIJ(JL2+JBAS,IL2+IBAS) = DCONJG(VIJ(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     SS BLOCKS
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            VIJ(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            VIJ(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            VIJ(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            VIJ(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
            VIJ(JS1+JBAS,IS1+IBAS) = DCONJG(VIJ(IS1+IBAS,JS1+JBAS))
            VIJ(JS2+JBAS,IS1+IBAS) = DCONJG(VIJ(IS1+IBAS,JS2+JBAS))
            VIJ(JS1+JBAS,IS2+IBAS) = DCONJG(VIJ(IS2+IBAS,JS1+JBAS))
            VIJ(JS2+JBAS,IS2+IBAS) = DCONJG(VIJ(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=J,NFUNA
            VIJ(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            VIJ(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            VIJ(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            VIJ(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
            VIJ(JS1+JBAS,IS1+IBAS) = DCONJG(VIJ(IS1+IBAS,JS1+JBAS))
            VIJ(JS2+JBAS,IS1+IBAS) = DCONJG(VIJ(IS1+IBAS,JS2+JBAS))
            VIJ(JS1+JBAS,IS2+IBAS) = DCONJG(VIJ(IS2+IBAS,JS1+JBAS))
            VIJ(JS2+JBAS,IS2+IBAS) = DCONJG(VIJ(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
2000  CONTINUE
C
100   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE VQDPOLE(VIJ,ITT,IQ,IX,JX,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    VV    VV  QQQQQQ   DDDDDDD  PPPPPPP   OOOOOO  LL      EEEEEEEE    C
C    VV    VV QQ    QQ  DD    DD PP    PP OO    OO LL      EE          C
C    VV    VV QQ    QQ  DD    DD PP    PP OO    OO LL      EE          C
C    VV    VV QQ    QQ  DD    DD PP    PP OO    OO LL      EEEEEE      C
C     VV  VV  QQ   QQQ  DD    DD PPPPPPP  OO    OO LL      EE          C
C      VVVV   QQ    QQ  DD    DD PP       OO    OO LL      EE          C
C       VV     QQQQQQ Q DDDDDDD  PP        OOOOOO  LLLLLLL EEEEEEEE    C
C                                                                      C
C -------------------------------------------------------------------- C
C  VQDPOLE CONSTRUCTS A MATRIX OF (u,T|SIG_Q.X.X'|v,T') OVERLAP        C
C  INTEGRALS OVER ALL BASIS FUNCTIONS, AND SAVES THE RESULT TO VIJ.    C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    ITT   = {1,2,3,4} -> {LL,LS,SL,SS} COMPONENT COMBINATION.         C
C    IQ    = {0,1,2,3} -> {0,X,Y,Z} THE COUPLING SIGMA MATRIX.         C
C    IX    = {1,2,3}   -> {X,Y,Z} FIRST CARTESIAN WEIGHTING FACTOR.    C
C    JX    = {1,2,3}   -> {X,Y,Z} SECOND CARTESIAN WEIGHTING FACTOR.   C
C    I1,I2 = {1,2,3,4} -> BASIS FUNCTION OVERLAPS FOR EQ-COEFFS.       C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 CONE
      COMPLEX*16 VIJ(MDM,MDM)
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 SLL(MBS,MBS,4),SSS(MBS,MBS,4)
C
      DIMENSION RC(MB2,MRC),EXPT(MBS,4),XYZ(3,4),APH(MB2),CP(MB2,3),
     &          PNC(MB2),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
      DATA PI/3.1415926535897932D0/
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     INITIALISE THE OVERLAP ARRAY
      DO I=1,NDIM
        DO J=1,NDIM
          VIJ(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     NON-RELATIVISTIC CONDITION
      IF(HMLTN.EQ.'NORL'.AND.ITT.GT.1) GOTO 100
C
100   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE VEFIELD(VIJ,ITT,IQ,IX,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      VV    VV EEEEEEEE FFFFFFFF IIII EEEEEEEE LL       DDDDDDD       C
C      VV    VV EE       FF        II  EE       LL       DD    DD      C
C      VV    VV EE       FF        II  EE       LL       DD    DD      C
C      VV    VV EEEEEE   FFFFFF    II  EEEEEE   LL       DD    DD      C
C       VV  VV  EE       FF        II  EE       LL       DD    DD      C
C        VVVV   EE       FF        II  EE       LL       DD    DD      C
C         VV    EEEEEEEE FF       IIII EEEEEEEE LLLLLLLL DDDDDDD       C
C                                                                      C
C -------------------------------------------------------------------- C
C  VEFIELD CONSTRUCTS A MATRIX OF (u,T|SIG.X/R^3|v,T') OVERLAP         C
C  INTEGRALS OVER ALL BASIS FUNCTIONS, AND SAVES THE RESULT TO VIJ.    C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    ITT   = {1,2,3,4} -> {LL,LS,SL,SS} COMPONENT COMBINATION.         C
C    IQ    = {0,1,2,3} -> {0,X,Y,Z} THE COUPLING SIGMA MATRIX.         C
C    IX    = {1,2,3}   -> {X,Y,Z} CARTESIAN WEIGHTING FACTOR.          C
C    I1,I2 = {1,2,3,4} -> BASIS FUNCTION OVERLAPS FOR EQ-COEFFS.       C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 CONE
      COMPLEX*16 VIJ(MDM,MDM)
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 SLL(MBS,MBS,4),SSS(MBS,MBS,4)
C
      DIMENSION RC(MB2,MRC),EXPT(MBS,4),XYZ(3,4),APH(MB2),CP(MB2,3),
     &          PNC(MB2),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
      DATA PI/3.1415926535897932D0/
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     INITIALISE THE OVERLAP ARRAY
      DO I=1,NDIM
        DO J=1,NDIM
          VIJ(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     NON-RELATIVISTIC CONDITION
      IF(HMLTN.EQ.'NORL'.AND.ITT.GT.1) GOTO 100
C
100   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE VKNETIC(VIJ,ITT,IQ,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      VV    VV KK    KK NN    NN EEEEEEEE TTTTTTTT IIII  CCCCCC       C
C      VV    VV KK   KK  NNN   NN EE          TT     II  CC    CC      C
C      VV    VV KK  KK   NNNN  NN EE          TT     II  CC            C
C      VV    VV KKKKK    NN NN NN EEEEEE      TT     II  CC            C
C       VV  VV  KK  KK   NN  NNNN EE          TT     II  CC            C
C        VVVV   KK   KK  NN   NNN EE          TT     II  CC    CC      C
C         VV    KK    KK NN    NN EEEEEEEE    TT    IIII  CCCCCC       C
C                                                                      C
C -------------------------------------------------------------------- C
C  VKNETIC CONSTRUCTS A MATRIX OF (u,T|SIG.X|v,T') OVERLAP             C
C  INTEGRALS OVER ALL BASIS FUNCTIONS, AND SAVES THE RESULT TO VIJ.    C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    ITT   = {1,2,3,4} -> {LL,LS,SL,SS} COMPONENT COMBINATION.         C
C    IQ    = {0,1,2,3} -> {0,X,Y,Z} THE COUPLING SIGMA MATRIX.         C
C    I1,I2 = {1,2,3,4} -> BASIS FUNCTION OVERLAPS FOR EQ-COEFFS.       C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 CONE
      COMPLEX*16 VIJ(MDM,MDM)
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 SLL(MBS,MBS,4),SSS(MBS,MBS,4)
C
      DIMENSION RC(MB2,MRC),EXPT(MBS,4),XYZ(3,4),APH(MB2),CP(MB2,3),
     &          PNC(MB2),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
      DATA PI/3.1415926535897932D0/
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     INITIALISE THE OVERLAP ARRAY
      DO I=1,NDIM
        DO J=1,NDIM
          VIJ(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     NON-RELATIVISTIC CONDITION
      IF(HMLTN.EQ.'NORL'.AND.ITT.GT.1) GOTO 100
C
100   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE VLPLACE(VIJ,ITT,IQ,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    VV    VV LL       PPPPPPP  LL          AA     CCCCCC  EEEEEEEE    C
C    VV    VV LL       PP    PP LL         AAAA   CC    CC EE          C
C    VV    VV LL       PP    PP LL        AA  AA  CC       EE          C
C    VV    VV LL       PP    PP LL       AA    AA CC       EEEEEE      C
C     VV  VV  LL       PPPPPPP  LL       AAAAAAAA CC       EE          C
C      VVVV   LL       PP       LL       AA    AA CC    CC EE          C
C       VV    LLLLLLLL PP       LLLLLLLL AA    AA  CCCCCC  EEEEEEEE    C
C                                                                      C
C -------------------------------------------------------------------- C
C  VLPLACE CONSTRUCTS A MATRIX OF (u,T|GRAD^2|v,T') OVERLAP            C
C  INTEGRALS OVER ALL BASIS FUNCTIONS, AND SAVES THE RESULT TO VIJ.    C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    ITT   = {1,2,3,4} -> {LL,LS,SL,SS} COMPONENT COMBINATION.         C
C    IQ    = {0,1,2,3} -> {0,X,Y,Z} THE COUPLING SIGMA MATRIX.         C
C    I1,I2 = {1,2,3,4} -> BASIS FUNCTION OVERLAPS FOR EQ-COEFFS.       C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 CONE
      COMPLEX*16 VIJ(MDM,MDM)
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 SLL(MBS,MBS,4),SSS(MBS,MBS,4)
C
      DIMENSION RC(MB2,MRC),EXPT(MBS,4),XYZ(3,4),APH(MB2),CP(MB2,3),
     &          PNC(MB2),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
      DATA PI/3.1415926535897932D0/
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     INITIALISE THE OVERLAP ARRAY
      DO I=1,NDIM
        DO J=1,NDIM
          VIJ(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     NON-RELATIVISTIC CONDITION
      IF(HMLTN.EQ.'NORL'.AND.ITT.GT.1) GOTO 100
C
100   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE VNCATRC(VIJ,ITT,IQ,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    VV    VV NN    NN  CCCCCC     AA   TTTTTTTT RRRRRRR   CCCCCC      C
C    VV    VV NNN   NN CC    CC   AAAA     TT    RR    RR CC    CC     C
C    VV    VV NNNN  NN CC        AA  AA    TT    RR    RR CC           C
C    VV    VV NN NN NN CC       AA    AA   TT    RR    RR CC           C
C     VV  VV  NN  NNNN CC       AAAAAAAA   TT    RRRRRRR  CC           C
C      VVVV   NN   NNN CC    CC AA    AA   TT    RR    RR CC    CC     C
C       VV    NN    NN  CCCCCC  AA    AA   TT    RR    RR  CCCCCC      C
C                                                                      C
C -------------------------------------------------------------------- C
C  VNCATRC CONSTRUCTS A MATRIX OF (u,T|nuc|v,T') OVERLAP               C
C  INTEGRALS OVER ALL BASIS FUNCTIONS, AND SAVES THE RESULT TO VIJ.    C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    ITT   = {1,2,3,4} -> {LL,LS,SL,SS} COMPONENT COMBINATION.         C
C    IQ    = {0,1,2,3} -> {0,X,Y,Z} THE COUPLING SIGMA MATRIX.         C
C    I1,I2 = {1,2,3,4} -> BASIS FUNCTION OVERLAPS FOR EQ-COEFFS.       C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 CONE
      COMPLEX*16 VIJ(MDM,MDM)
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 SLL(MBS,MBS,4),SSS(MBS,MBS,4)
C
C      DIMENSION RC(MB2,MRC),EXPT(MBS,4),XYZ(3,4),APH(MB2),CP(MB2,3),
C     &          PNC(MB2),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
      DATA PI/3.1415926535897932D0/
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     INITIALISE THE OVERLAP ARRAY
      DO I=1,NDIM
        DO J=1,NDIM
          VIJ(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     NON-RELATIVISTIC CONDITION
      IF(HMLTN.EQ.'NORL'.AND.ITT.GT.1) GOTO 100
C
100   CONTINUE
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C  [10] PLOTS: AMPLITUDES AND FIELDS/POTENTIALS IN DATA FILES.         C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] FIELDS: CALCULATE AMPLITUDES, FIELDS AND POTENTIALS.           C
C   [B] ARRYPLT: EXPORT ARRAY TO EXTERNAL DATA FILE AND PLOT.          C
C   [C] ORBCOEF: PRINT LIST OF EXPANSION COEFFICIENTS FOR AN ORBITAL.  C
C   [D] AMPLTDE: PLOT A DIRAC SPINOR AMPLITUDE ALONG ONE DIRECTION.    C
C   [E] J4EQSUM: PLOT 4-CURRENT ALONG ONE DIRECTION FOR ORBITAL PAIR.  C
C   [F] J4DIRCT: PLOT 4-CURRENT ALONG ONE DIRECTION FOR ORBITAL PAIR.  C
C   [G] POTENTL: PLOT 4-POTENTIAL ALONG ONE DIRECTION FOR ORBITAL PAIR.C
C   [H] ELCTRCF: PLOT E FIELD ALONG ONE DIRECTION FOR ORBITAL PAIR.    C
C   [I] MAGNTCF: PLOT B FIELD ALONG ONE DIRECTION FOR ORBITAL PAIR.    C
C   [J] GNUMAKE: GENERATE A GNUPLOT MAKE FILE FOR A DATA SET.          C
C   [K] CLEBSCH: CLEBSCH-GORDON COEFFICIENT FOR KQN,MQN.               C
C   [L] SPHHRM: VALUE OF Y_L^M AT TWO GIVEN ANGLES.                    C
C   [M] PLGNDR: RETURNS AN ASSOCIATED LEGENDRE POLYNOMIAL P_L^M(X).    C
C   [N] NFACT: INTEGER RESULT OF FACTORIAL N!                          C
C**********************************************************************C
C
C
      SUBROUTINE FIELDS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          FFFFFFFF IIII EEEEEEEE LL       DDDDDDD   SSSSSS            C
C          FF        II  EE       LL       DD    DD SS    SS           C
C          FF        II  EE       LL       DD    DD SS                 C
C          FFFFFF    II  EEEEEE   LL       DD    DD  SSSSSS            C
C          FF        II  EE       LL       DD    DD       SS           C
C          FF        II  EE       LL       DD    DD SS    SS           C
C          FF       IIII EEEEEEEE LLLLLLLL DDDDDDD   SSSSSS            C
C                                                                      C
C -------------------------------------------------------------------- C
C  FIELDS CALCULATES AMPLITUDES, EM FIELDS AND POTENTIALS, STORES IN   C
C  EXTERNAL DATA FILES AND PLOTS WHEN CALLED.                          C
C -------------------------------------------------------------------- C
C  EQFILE ONLY GENERATES ELL0, ESS0 AND MAYBE ELSI, AND IN GENERAL     C
C  THESE CALCULATIONS REQUIRE ETT'Q, SO GENERATE AS NEEDED INSTEAD.    C
C -------------------------------------------------------------------- C
C  FOR NOW PLOTS ARE GENERATED ALONG ONE DIRECTION ONLY, AND FOR IORB  C
C  AND JORB ORBITAL PAIRS. CAN EXTEND TO SURFACE AND DENSITY PLOTS, AS C
C  WELL AS OVERALL MOLECULAR AMPLITUDES AND FIELDS.                    C
C**********************************************************************C
      PARAMETER(MBS=26,MCT=6,MKP=9)
C
      CHARACTER*7  PTYPE(10)
      CHARACTER*16 HMS
C
      COMMON/PLOT/NPTYPE,PTYPE
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
C
C     RECORD TIME
      CALL CPU_TIME(TDUM)
C
C     EQ-COEFF AND R-INT TIME INITIALISATION
      TELL = 0.0D0
      TESS = 0.0D0
      TELS = 0.0D0
C
C     LOOP OVER NUMBER OF REQUESTED FIELD PLOTS
      DO N=1,NPTYPE
C
C      PRINT A TITLE
        WRITE(6, *) ' '
        WRITE(7, *) ' '
        WRITE(6, *) 'Plot for PTYPE = ',PTYPE(N)
        WRITE(7, *) 'Plot for PTYPE = ',PTYPE(N)
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
C
        IF(PTYPE(N).EQ.'ORBCOEF') THEN
          CALL ORBCOEF(NSHIFT+5)
        ELSEIF(PTYPE(N).EQ.'AMPLTDE') THEN
          CALL AMPLTDE(NSHIFT+1)
        ELSEIF(PTYPE(N).EQ.'J4EQSUM') THEN
          CALL J4EQSUM(NSHIFT+1,NSHIFT+1)
        ELSEIF(PTYPE(N).EQ.'J4DIRCT') THEN
          CALL J4DIRCT(NSHIFT+1,NSHIFT+1)
        ELSEIF(PTYPE(N).EQ.'POTENTL') THEN
          CALL POTENTL(NSHIFT+1,NSHIFT+1)
        ELSEIF(PTYPE(N).EQ.'ELCTRCF') THEN
          CALL ELCTRCF(NSHIFT+1,NSHIFT+1)
        ELSEIF(PTYPE(N).EQ.'MAGNTCF') THEN
          CALL MAGNTCF(NSHIFT+1,NSHIFT+1)
        ELSE
          WRITE(6, *) 'In FIELDS: plotting option is not available.'
          WRITE(7, *) 'In FIELDS: plotting option is not available.'
        ENDIF
C
C     END LOOP OVER FIELD PLOTS
      ENDDO
C
C     TOTAL TIME TAKEN
      CALL CPU_TIME(TPLT)
      TPLT = TPLT-TDUM
C
20    FORMAT(1X,A,37X,A)
      WRITE(6,*) REPEAT('-',72)
      WRITE(7,*) REPEAT('-',72)
      WRITE(6,20) 'Time in EMAKE (LL):',HMS(TELL)
      WRITE(7,20) 'Time in EMAKE (LL):',HMS(TELL)
      WRITE(6,20) 'Time in EMAKE (SS):',HMS(TESS)
      WRITE(7,20) 'Time in EMAKE (SS):',HMS(TESS)
      WRITE(6,20) 'Time in EMAKE (LS):',HMS(TELS)
      WRITE(7,20) 'Time in EMAKE (LS):',HMS(TELS)
      WRITE(6,*) REPEAT('=',72)
      WRITE(7,*) REPEAT('=',72)
C
      RETURN
      END
C
C
      SUBROUTINE ARRYPLT(ARRAY,TITLE,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       AA    RRRRRRR  RRRRRRR  YY    YY PPPPPPP  LL      TTTTTTTT     C
C      AAAA   RR    RR RR    RR YY    YY PP    PP LL         TT        C
C     AA  AA  RR    RR RR    RR YY    YY PP    PP LL         TT        C
C    AA    AA RR    RR RR    RR  YY  YY  PP    PP LL         TT        C
C    AAAAAAAA RRRRRRR  RRRRRRR    YYYY   PPPPPPP  LL         TT        C
C    AA    AA RR    RR RR    RR    YY    PP       LL         TT        C
C    AA    AA RR    RR RR    RR    YY    PP       LLLLLLLL   TT        C
C                                                                      C
C -------------------------------------------------------------------- C
C  ARRYPLT EXPORTS AN ARRAY TO AN EXTERNAL DATA FILE AND PLOTS IT.     C
C**********************************************************************C
      PARAMETER(MDM=1200)
C
      CHARACTER*80 TITLE
C
      DIMENSION ARRAY(MDM,MDM)
C
C     PRINT TO EXTERNAL DATA FILE
      OPEN(UNIT=8,FILE="plots/"//TRIM(TITLE)//".dat",STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO I=1,NDIM
        WRITE(8, *) (ARRAY(I,J),J=1,NDIM)
      ENDDO
      CLOSE(UNIT=8)
C
      XEND = DFLOAT(NDIM)-0.5D0
      YEND = DFLOAT(NDIM)-0.5D0
C
C     WRITE GNUPLOT MAKE FILE
      OPEN(UNIT=9,FILE='plots/'//TRIM(TITLE)//'.gnuplot',
     &                                                 STATUS='REPLACE')
      WRITE(9,'(A)') '#'//TRIM(TITLE)//'.gnuplot'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '#  Usage:'
      WRITE(9,'(A)') '#  gnuplot < '//TRIM(TITLE)//'.gnuplot'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Terminal output specs'
      WRITE(9,'(A)') 'set terminal pdf size 4,4'
      WRITE(9,'(A)') 'set output "plots/'//TRIM(TITLE)//'.pdf"'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') 'load "plots/pals/rdylgn.pal"'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Axes and title'
      WRITE(9,'(A)') 'set title sprintf("'//TRIM(TITLE)//'")'
      WRITE(9,'(A,F6.1,A)') 'set xrange [-0.5:',XEND,']'
      WRITE(9,'(A,F6.1,A)') 'set yrange [-0.5:',YEND,'] reverse'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Plot data to file'
      WRITE(9,'(A,I2,A,I2,A)') 'plot "plots/'//TRIM(TITLE)//'.dat"'
     &                                    //' matrix with image notitle'
      CLOSE(UNIT=9)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(TITLE)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(TITLE)//'.pdf')
C
      RETURN
      END
C
C
      SUBROUTINE ORBCOEF(IORB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      OOOOOO  RRRRRRR  BBBBBBB   CCCCCC   OOOOOO  EEEEEEE FFFFFFF     C
C     OO    OO RR    RR BB    BB CC    CC OO    OO EE      FF          C
C     OO    OO RR    RR BB    BB CC       OO    OO EE      FF          C
C     OO    OO RR    RR BBBBBBB  CC       OO    OO EEEEE   FFFFF       C
C     OO    OO RRRRRRR  BB    BB CC       OO    OO EE      FF          C
C     OO    OO RR    RR BB    BB CC    CC OO    OO EE      FF          C
C      OOOOOO  RR    RR BBBBBBB   CCCCCC   OOOOOO  EEEEEEE FF          C
C                                                                      C
C -------------------------------------------------------------------- C
C  ORBCOEF PRINTS OUT THE EXPANSION COEFFICIENTS FOR IORB, GROUPED     C
C  BY QUANTUM NUMBERS OF BASIS FUNCTIONS.                              C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*2 ELMNT(120),ELA
C
      COMPLEX*16 C(MDM,MDM)
C
      COMMON/ATOM/ELMNT
      COMMON/COEF/C
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     HEADER FOR COEFFICIENT LIST
20    FORMAT(1X,'****************',A,I3,'  ****************')
21    FORMAT(1X,'Center ',I2,4X,'(',A,')',5X,'#fns = ',I2,6X,'KQN =',
     &                           I2,6X,'LQN = ',I1,6X,'MQN =',A,I1,'/2')
22    FORMAT(1X,'IADD',9X,'Exponent',10X,A,7X,A)
23    FORMAT(1X,I4,4X,ES13.6,10X,ES17.10,7X,ES17.10)
      WRITE(6,20) '  Expansion coefficients for IORB =',IORB
      WRITE(7,20) '  Expansion coefficients for IORB =',IORB
C
      I = IORB
C
C     LOOP OVER NUCLEAR CENTERS
      DO ICNT=1,NCNT
C
C       SAVE ELEMENT LABEL
        ELA = ELMNT(IZNUC(ICNT))
C     
C       LOOP OVER ALL MQNS FOR THIS CENTER
        MMAX = LMAX(ICNT)+1
        DO IM=1,MMAX
C
C         MAGNITUDE OF THIS MQN
          MQN = 2*IM-1
C
C         MQN < 0: LOOP OVER ALL KQN WITH THIS |MQN|
          DO IK=MQN,NKAP(ICNT)
C
C           LARGE-COMPONENT ADDRESS OFFSET
            IAD = LARGE(ICNT,IK,MQN  )
C
C           LQN AND KQN VALUES
            IF(MOD(IK,2).EQ.0) THEN
              KQN = IK/2
              LQN = KQN
            ELSE
              KQN =-(IK+1)/2
              LQN =-KQN-1
            ENDIF
C
C           NUMBER OF BASIS FUNCTIONS FOR THIS LQN
            NBAS = NFUNCT(LQN+1,ICNT)
C
C           ATOMIC ADDRESS DETAILS
            WRITE(6, *) REPEAT('-',72)
            WRITE(7, *) REPEAT('-',72)
            WRITE(6,21) ICNT,ELA,NBAS,KQN,LQN,'-',MQN
            WRITE(7,21) ICNT,ELA,NBAS,KQN,LQN,'-',MQN
            WRITE(6, *) REPEAT('-',72)
            WRITE(7, *) REPEAT('-',72)
            WRITE(6,22) 'Re(CL(IADD,IORB))','Re(CS(IADD,IORB))'
            WRITE(7,22) 'Re(CL(IADD,IORB))','Re(CS(IADD,IORB))'
            WRITE(6, *) REPEAT('-',72)
            WRITE(7, *) REPEAT('-',72)
C
C           LIST THE EXPONENTS AND CORRESPONDING EXP. COEFFS
            DO IBAS=1,NBAS
              CL = DREAL(C(IAD+IBAS       ,IORB))
              CS = DREAL(C(IAD+IBAS+NSHIFT,IORB))
              WRITE(6,23) IAD+IBAS,EXPSET(IBAS,LQN+1,ICNT),CL,CS
              WRITE(7,23) IAD+IBAS,EXPSET(IBAS,LQN+1,ICNT),CL,CS
            ENDDO
          ENDDO
C
C         MQN > 0: LOOP OVER ALL KQN WITH THIS |MQN|
          DO IK=MQN,NKAP(ICNT)
C
C           LARGE-COMPONENT ADDRESS OFFSET
            IAD = LARGE(ICNT,IK,MQN+1)
C
C           LQN AND KQN VALUES
            IF(MOD(IK,2).EQ.0) THEN
              KQN = IK/2
              LQN = KQN
            ELSE
              KQN =-(IK+1)/2
              LQN =-KQN-1
            ENDIF
C
C           NUMBER OF BASIS FUNCTIONS FOR THIS LQN
            NBAS = NFUNCT(LQN+1,ICNT)
C
C           ATOMIC ADDRESS DETAILS
            WRITE(6, *) REPEAT('-',72)
            WRITE(7, *) REPEAT('-',72)
            WRITE(6,21) ICNT,ELA,NBAS,KQN,LQN,'+',MQN
            WRITE(7,21) ICNT,ELA,NBAS,KQN,LQN,'+',MQN
            WRITE(6, *) REPEAT('-',72)
            WRITE(7, *) REPEAT('-',72)
            WRITE(6,22) 'Re(CL(IADD,IORB))','Re(CS(IADD,IORB))'
            WRITE(7,22) 'Re(CL(IADD,IORB))','Re(CS(IADD,IORB))'
            WRITE(6, *) REPEAT('-',72)
            WRITE(7, *) REPEAT('-',72)
C
C           LIST THE EXPONENTS AND CORRESPONDING EXP. COEFFS
            DO IBAS=1,NBAS
              CL = DREAL(C(IAD+IBAS       ,IORB))
              CS = DREAL(C(IAD+IBAS+NSHIFT,IORB))
              WRITE(6,23) IAD+IBAS,EXPSET(IBAS,LQN+1,ICNT),CL,CS
              WRITE(7,23) IAD+IBAS,EXPSET(IBAS,LQN+1,ICNT),CL,CS
            ENDDO
          ENDDO       
        ENDDO    
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE AMPLTDE(IORB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      AA    MM       MM PPPPPPP  LL      TTTTTTTT DDDDDDD  EEEEEEEE   C
C     AAAA   MMM     MMM PP    PP LL         TT    DD    DD EE         C
C    AA  AA  MMMM   MMMM PP    PP LL         TT    DD    DD EE         C
C   AA    AA MM MM MM MM PP    PP LL         TT    DD    DD EEEEEE     C
C   AAAAAAAA MM  MMM  MM PPPPPPP  LL         TT    DD    DD EE         C
C   AA    AA MM   M   MM PP       LL         TT    DD    DD EE         C
C   AA    AA MM       MM PP       LLLLLLLL   TT    DDDDDDD  EEEEEEEE   C
C                                                                      C
C -------------------------------------------------------------------- C
C  AMPLTDE GENERATES A DATA SET AND GNUPLOT MAKE FILE FOR THE DIRAC    C
C  SPINOR AMPLITUDE OF ORBITAL IORB, ALONG THE Z-AXIS.                 C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4  HMLTN
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(4)
C
      COMPLEX*16 C(MDM,MDM)
C
      DIMENSION PSI(4)
C
      COMMON/COEF/C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     FILE NAME AND TITLES
      XOUT   = TRIM(MOLCL)//'_'//HMLTN//'_'//'AMPLTDE'
      TITLE  = 'AMPLTDE'//' for '//TRIM(MOLCL)//' of type '//HMLTN
      XAXIS  = 'z (a.u.)'
      YAXIS  = '{/Symbol y}(z)'
      KEY(1) = '{/Symbol y}^{L}_{u}(z)'
      KEY(2) = '{/Symbol y}^{L}_{d}(z)'
      KEY(3) = '{/Symbol y}^{S}_{u}(z)'
      KEY(4) = '{/Symbol y}^{S}_{d}(z)'
C
      I = IORB
C
C     DATA PLOTTING DETAILS
      NP = 1000
      Z0 = 0.0D0
      ZF = 5.0D0
      HS = (ZF-Z0)/NP
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C     BEGIN LOOP OVER ALL DATA POINTS
      DO IP=0,NP
C
C       SPECIFY CARTESIAN COORDINATES
        X = 0.0D0
        Y = 0.0D0
        Z = Z0 + HS*IP
C
C       DFNOTE: THIS IS DUMMY CODE JUST TO GET A DATA FILE
        PSI(1) = 1.0D0
        PSI(2) = Z*Z
        PSI(3) = DCOS(Z)
        PSI(4) = DEXP(-Z)
        
        WRITE(8,*) Z,(PSI(I),I=1,4)
C
C     CLOSE LOOP OVER DATA POINTS
      ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNUMAKE(XOUT,TITLE,XAXIS,YAXIS,4,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
      RETURN
      END
C
C
      SUBROUTINE J4EQSUM(IORB,JORB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      JJJJJ    44   EEEEEEEE  QQQQQQ    SSSSSS  UU    UU MM       MM  C
C        JJ    44    EE       QQ    QQ  SS    SS UU    UU MMM     MMM  C
C        JJ   44     EE       QQ    QQ  SS       UU    UU MMMM   MMMM  C
C        JJ  44 44   EEEEEE   QQ    QQ   SSSSSS  UU    UU MM MM MM MM  C
C        JJ 44444444 EE       QQ   QQQ        SS UU    UU MM  MMM  MM  C
C  JJ    JJ     44   EE       QQ    QQ  SS    SS UU    UU MM   M   MM  C
C   JJJJJJ      44   EEEEEEEE  QQQQQQ Q  SSSSSS   UUUUUU  MM       MM  C
C                                                                      C
C -------------------------------------------------------------------- C
C  J4EQSUM CREATES A PLOT OF THE 4-CURRENT OVERLAP BETWEEN IORB AND    C
C  JORB BY USE OF THE EQ-COEFFICIENTS AND GAUSSIAN PRODUCT THEOREM.    C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4  HMLTN
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(4)
C
      COMPLEX*16 C(MDM,MDM)
C
      COMMON/COEF/C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     FILE NAME AND TITLES
      XOUT   = TRIM(MOLCL)//'_'//HMLTN//'_'//'J4EQSUM'
      TITLE  = 'J4EQSUM'//' for '//TRIM(MOLCL)//' of type '//HMLTN
      XAXIS  = 'z (a.u.)'
      YAXIS  = '{j}_{/Symbol u}(z)'
      KEY(1) = '{/Symbol r}(z)'
      KEY(2) = 'j_{x}(z)'
      KEY(3) = 'j_{y}(z)'
      KEY(4) = 'j_{z}(z)'
C
      I = IORB
      J = JORB
C
C     DATA PLOTTING DETAILS
      NP = 1000
      Z0 = 0.0D0
      ZF = 5.0D0
      HS = (ZF-Z0)/NP
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C     BEGIN LOOP OVER ALL DATA POINTS
      DO IP=0,NP
C
C       SPECIFY CARTESIAN COORDINATES
        X = 0.0D0
        Y = 0.0D0
        Z = Z0 + HS*IP
C
C     CLOSE LOOP OVER DATA POINTS
      ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNUMAKE(XOUT,TITLE,XAXIS,YAXIS,4,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
      RETURN
      END
C
C
      SUBROUTINE J4DIRCT(IORB,JORB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      JJJJJ    44   DDDDDDD  IIII RRRRRRR  EEEEEEEE CCCCCC TTTTTTTT   C
C        JJ    44    DD    DD  II  RR    RR EE      CC    CC   TT      C
C        JJ   44     DD    DD  II  RR    RR EE      CC         TT      C
C        JJ  44 44   DD    DD  II  RR    RR EEEEEE  CC         TT      C
C        JJ 44444444 DD    DD  II  RRRRRRR  EE      CC         TT      C
C  JJ    JJ     44   DD    DD  II  RR    RR EE      CC    CC   TT      C
C   JJJJJJ      44   DDDDDDD  IIII RR    RR EEEEEEEE CCCCCC    TT      C
C                                                                      C
C -------------------------------------------------------------------- C
C  J4DIRCT CREATES A PLOT OF THE 4-CURRENT OVERLAP BETWEEN IORB AND    C
C  JORB BY DIRECT MULTIPLICATION OF GAUSSIANS AND SPHERICAL SPINORS.   C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4  HMLTN
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(4)
C
      COMPLEX*16 C(MDM,MDM)
C
      COMMON/COEF/C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     FILE NAME AND TITLES
      XOUT   = TRIM(MOLCL)//'_'//HMLTN//'_'//'J4DIRCT'
      TITLE  = 'J4DIRCT'//' for '//TRIM(MOLCL)//' of type '//HMLTN
      XAXIS  = 'z (a.u.)'
      YAXIS  = '{j}_{/Symbol u}(z)'
      KEY(1) = '{/Symbol r}(z)'
      KEY(2) = 'j_{x}(z)'
      KEY(3) = 'j_{y}(z)'
      KEY(4) = 'j_{z}(z)'
C
      I = IORB
      J = JORB
C
C     DATA PLOTTING DETAILS
      NP = 1000
      Z0 = 0.0D0
      ZF = 5.0D0
      HS = (ZF-Z0)/NP
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C     BEGIN LOOP OVER ALL DATA POINTS
      DO IP=0,NP
C
C       SPECIFY CARTESIAN COORDINATES
        X = 0.0D0
        Y = 0.0D0
        Z = Z0 + HS*IP
C
C     CLOSE LOOP OVER DATA POINTS
      ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNUMAKE(XOUT,TITLE,XAXIS,YAXIS,4,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
      RETURN
      END
C
C
      SUBROUTINE POTENTL(IORB,JORB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     PPPPPPP   OOOOOO TTTTTTTT EEEEEEEE NN    NN TTTTTTTT LL          C
C     PP    PP OO    OO   TT    EE       NNN   NN    TT    LL          C
C     PP    PP OO    OO   TT    EE       NNNN  NN    TT    LL          C
C     PP    PP OO    OO   TT    EEEEEE   NN NN NN    TT    LL          C
C     PPPPPPP  OO    OO   TT    EE       NN  NNNN    TT    LL          C
C     PP       OO    OO   TT    EE       NN   NNN    TT    LL          C
C     PP        OOOOOO    TT    EEEEEEEE NN    NN    TT    LLLLLLLL    C
C                                                                      C
C -------------------------------------------------------------------- C
C  POTENTL CREATES A PLOT OF THE 4-POTENTIAL OVERLAP BETWEEN IORB AND  C
C  JORB BY USE OF THE EQ-COEFFICIENTS AND BOYS INTEGRALS.              C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4  HMLTN
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(4)
C
      COMPLEX*16 C(MDM,MDM)
C
      COMMON/COEF/C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     FILE NAME AND TITLES
      XOUT   = TRIM(MOLCL)//'_'//HMLTN//'_'//'POTENTL'
      TITLE  = 'POTENTL'//' for '//TRIM(MOLCL)//' of type '//HMLTN
      XAXIS  = 'z (a.u.)'
      YAXIS  = '{A}_{/Symbol u}(z)'
      KEY(1) = '{/Symbol f}(z)'
      KEY(2) = 'A_{x}(z)'
      KEY(3) = 'A_{y}(z)'
      KEY(4) = 'A_{z}(z)'
C
      I = IORB
      J = JORB
C
C     DATA PLOTTING DETAILS
      NP = 1000
      Z0 = 0.0D0
      ZF = 5.0D0
      HS = (ZF-Z0)/NP
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C     BEGIN LOOP OVER ALL DATA POINTS
      DO IP=0,NP
C
C       SPECIFY CARTESIAN COORDINATES
        X = 0.0D0
        Y = 0.0D0
        Z = Z0 + HS*IP
C
C     CLOSE LOOP OVER DATA POINTS
      ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNUMAKE(XOUT,TITLE,XAXIS,YAXIS,4,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
      RETURN
      END
C
C
      SUBROUTINE ELCTRCF(IORB,JORB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     EEEEEEEE LL       CCCCCC TTTTTTTT RRRRRRR   CCCCCC  FFFFFFFF     C
C     EE       LL      CC    CC   TT    RR    RR CC    CC FF           C
C     EE       LL      CC         TT    RR    RR CC       FF           C
C     EEEEEE   LL      CC         TT    RR    RR CC       FFFFFF       C
C     EE       LL      CC         TT    RRRRRRR  CC       FF           C
C     EE       LL      CC    CC   TT    RR    RR CC    CC FF           C
C     EEEEEEEE LLLLLLLL CCCCCC    TT    RR    RR  CCCCCC  FF           C
C                                                                      C
C -------------------------------------------------------------------- C
C  ELCTRCF CREATES A PLOT OF THE ELECTRIC FIELD OVERLAP BETWEEN IORB   C
C  AND JORB BY USE OF THE EQ-COEFFICIENTS AND BOYS INTEGRALS.          C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4  HMLTN
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(4)
C
      COMPLEX*16 C(MDM,MDM)
C
      COMMON/COEF/C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     FILE NAME AND TITLES
      XOUT   = TRIM(MOLCL)//'_'//HMLTN//'_'//'ELCTRCF'
      TITLE  = 'ELCTRCF'//' for '//TRIM(MOLCL)//' of type '//HMLTN
      XAXIS  = 'z (a.u.)'
      YAXIS  = '{E}_{i}(z)'
      KEY(1) = 'E_{x}(z)'
      KEY(2) = 'E_{y}(z)'
      KEY(3) = 'E_{z}(z)'
      KEY(4) = '|E|(z)'
C
      I = IORB
      J = JORB
C
C     DATA PLOTTING DETAILS
      NP = 1000
      Z0 = 0.0D0
      ZF = 5.0D0
      HS = (ZF-Z0)/NP
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C     BEGIN LOOP OVER ALL DATA POINTS
      DO IP=0,NP
C
C       SPECIFY CARTESIAN COORDINATES
        X = 0.0D0
        Y = 0.0D0
        Z = Z0 + HS*IP
C
C     CLOSE LOOP OVER DATA POINTS
      ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNUMAKE(XOUT,TITLE,XAXIS,YAXIS,4,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
      RETURN
      END
C
C
      SUBROUTINE MAGNTCF(IORB,JORB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   MM       MM    AA     GGGGGG  NN    NN TTTTTTTT CCCCCC  FFFFFFFF   C
C   MMM     MMM   AAAA   GG    GG NNN   NN    TT   CC    CC FF         C
C   MMMM   MMMM  AA  AA  GG       NNNN  NN    TT   CC       FF         C
C   MM MM MM MM AA    AA GG       NN NN NN    TT   CC       FFFFFF     C
C   MM  MMM  MM AAAAAAAA GG   GGG NN  NNNN    TT   CC       FF         C
C   MM   M   MM AA    AA GG    GG NN   NNN    TT   CC    CC FF         C
C   MM       MM AA    AA  GGGGGG  NN    NN    TT    CCCCCC  FF         C
C                                                                      C
C -------------------------------------------------------------------- C
C  MAGNTCF CREATES A PLOT OF THE MAGNETIC FIELD OVERLAP BETWEEN IORB   C
C  AND JORB BY USE OF THE EQ-COEFFICIENTS AND BOYS INTEGRALS.          C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=6,MKP=9)
C
      CHARACTER*4  HMLTN
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(4)
C
      COMPLEX*16 C(MDM,MDM)
C
      COMMON/COEF/C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,ILEV,ISWZ,IEQS,IERC,IPAR,ICOR
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C     FILE NAME AND TITLES
      XOUT   = TRIM(MOLCL)//'_'//HMLTN//'_'//'MAGNTCF'
      TITLE  = 'MAGNTCF'//' for '//TRIM(MOLCL)//' of type '//HMLTN
      XAXIS  = 'z (a.u.)'
      YAXIS  = '{B}_{i}(z)'
      KEY(1) = 'B_{x}(z)'
      KEY(2) = 'B_{y}(z)'
      KEY(3) = 'B_{z}(z)'
      KEY(4) = '|B|(z)'
C
      I = IORB
      J = JORB
C
C     DATA PLOTTING DETAILS
      NP = 1000
      Z0 = 0.0D0
      ZF = 5.0D0
      HS = (ZF-Z0)/NP
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C     BEGIN LOOP OVER ALL DATA POINTS
      DO IP=0,NP
C
C       SPECIFY CARTESIAN COORDINATES
        X = 0.0D0
        Y = 0.0D0
        Z = Z0 + HS*IP
C
C     CLOSE LOOP OVER DATA POINTS
      ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNUMAKE(XOUT,TITLE,XAXIS,YAXIS,4,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
      RETURN
      END
C
C
      SUBROUTINE GNUMAKE(XOUT,TITLE,XAXIS,YAXIS,NDAT,KEY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   GGGGGG  NN    NN UU    UU MM       MM    AA    KK    KK EEEEEEEE   C
C  GG    GG NNN   NN UU    UU MMM     MMM   AAAA   KK   KK  EE         C
C  GG       NNNN  NN UU    UU MMMM   MMMM  AA  AA  KK  KK   EE         C
C  GG       NN NN NN UU    UU MM MM MM MM AA    AA KKKKK    EEEEEE     C
C  GG   GGG NN  NNNN UU    UU MM  MMM  MM AAAAAAAA KK  KK   EE         C
C  GG    GG NN   NNN UU    UU MM   M   MM AA    AA KK   KK  EE         C
C   GGGGGG  NN    NN  UUUUUU  MM       MM AA    AA KK    KK EEEEEEEE   C
C                                                                      C
C -------------------------------------------------------------------- C
C  GNUMAKE IS A CONTROLLING ROUTINE THAT GENERATES A GNUPLOT MAKE FILE C
C  FOR A SET OF DATA POINTS.                                           C
C**********************************************************************C
C
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(NDAT)
C
      OPEN(UNIT=9,FILE='plots/'//TRIM(XOUT)//'.gnuplot',
     &                                                 STATUS='REPLACE')
      WRITE(9,'(A)') '#'//TRIM(XOUT)//'.gnuplot'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '#  Usage:'
      WRITE(9,'(A)') '#  gnuplot < '//TRIM(XOUT)//'.gnuplot'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Plot raw data'
      DO N=1,NDAT
        IF(NDAT.EQ.1) THEN
          WRITE(9,'(A,I2,A)') 'plot "plots/'//TRIM(XOUT)//'.dat"'
     &                        //' using 1:',N+1,' with lines'
        ELSEIF(NDAT.GT.1.AND.N.EQ.1) THEN
          WRITE(9,'(A,I2,A)') 'plot "plots/'//TRIM(XOUT)//'.dat"'
     &                        //' using 1:',N+1,' with lines,\'
        ELSEIF(NDAT.GT.1.AND.N.GT.1.AND.N.LT.NDAT) THEN
          WRITE(9,'(A,I2,A)') '     "plots/'//TRIM(XOUT)//'.dat"'
     &                        //' using 1:',N+1,' with lines,\'
        ELSEIF(NDAT.GT.1.AND.N.EQ.NDAT) THEN
          WRITE(9,'(A,I2,A)') '     "plots/'//TRIM(XOUT)//'.dat"'
     &                        //' using 1:',N+1,' with lines'
        ENDIF
      ENDDO
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Terminal output specs'
      WRITE(9,'(A)') 'set terminal pdf enhance font "palatino,10"'
      WRITE(9,'(A)') 'set output "plots/'//TRIM(XOUT)//'.pdf"'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Axes and title'
      WRITE(9,'(A)') 'set title sprintf("'//TRIM(TITLE)//'")'
      WRITE(9,'(A)') 'set xlabel "'//TRIM(XAXIS)//'"'
      WRITE(9,'(A)') 'set ylabel "'//TRIM(YAXIS)//'"'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Plotting range'
      WRITE(9,'(A)') 'xmin = GPVAL_X_MIN'
      WRITE(9,'(A)') 'xmax = GPVAL_X_MAX'
      WRITE(9,'(A)') 'ymin = GPVAL_Y_MIN'
      WRITE(9,'(A)') 'ymax = GPVAL_Y_MAX'
      WRITE(9,'(A)') 'set xrange [xmin:xmax] noreverse nowriteback'
      WRITE(9,'(A)') 'set yrange [ymin:ymax] noreverse nowriteback'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Grid style'
      WRITE(9,'(A)') 'set style line 102 lc rgb"#808080" lt 0 lw 1'
      WRITE(9,'(A)') 'set grid back ls 102'  
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Plot data to file'
      DO N=1,NDAT
        IF(NDAT.EQ.1) THEN
          WRITE(9,'(A,I2,A,I2,A)') 'plot "plots/'//TRIM(XOUT)//'.dat"'
     &                        //' using 1:',N+1,' with lines ls ',N+1,
     &                                 ' title "'//TRIM(KEY(N))//'"'
        ELSEIF(NDAT.GT.1.AND.N.EQ.1) THEN
          WRITE(9,'(A,I2,A,I2,A)') 'plot "plots/'//TRIM(XOUT)//'.dat"'
     &                        //' using 1:',N+1,' with lines ls ',N+1,
     &                                 ' title "'//TRIM(KEY(N))//'",\'
        ELSEIF(NDAT.GT.1.AND.N.GT.1.AND.N.LT.NDAT) THEN
          WRITE(9,'(A,I2,A,I2,A)') '     "plots/'//TRIM(XOUT)//'.dat"'
     &                        //' using 1:',N+1,' with lines ls ',N+1,
     &                                 ' title "'//TRIM(KEY(N))//'",\'
        ELSEIF(NDAT.GT.1.AND.N.EQ.NDAT) THEN
          WRITE(9,'(A,I2,A,I2,A)') '     "plots/'//TRIM(XOUT)//'.dat"'
     &                        //' using 1:',N+1,' with lines ls ',N+1,
     &                                 ' title "'//TRIM(KEY(N))//'"'
        ENDIF
      ENDDO
      CLOSE(UNIT=9)
C
      WRITE(6, *) 'Created command file "'//TRIM(XOUT)//'.gnuplot."'
      WRITE(7, *) 'Created command file "'//TRIM(XOUT)//'.gnuplot."'
C
      RETURN
      END
C
C
      FUNCTION CLEBSCH(KQN,MQN,NSGN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     CCCCCC  LL       EEEEEEEE BBBBBBB   SSSSSS   CCCCCC  HH    HH    C
C    CC    CC LL       EE       BB    BB SS    SS CC    CC HH    HH    C
C    CC       LL       EE       BB    BB SS       CC       HH    HH    C
C    CC       LL       EEEEEE   BBBBBBB   SSSSSS  CC       HHHHHHHH    C
C    CC       LL       EE       BB    BB       SS CC       HH    HH    C
C    CC    CC LL       EE       BB    BB SS    SS CC    CC HH    HH    C
C     CCCCCC  LLLLLLLL EEEEEEEE BBBBBBB   SSSSSS   CCCCCC  HH    HH    C
C                                                                      C
C -------------------------------------------------------------------- C
C     CLEBSCH DETERMINES A CLEBSCH-GORDON COEFFICIENT FOR KQN,MQN.     C
C -------------------------------------------------------------------- C
C     KQN IS FINE BUT MQN IS DOUBLE THE ACTUAL VALUE, AND NSGN IS      C
C     THE PARITY NUMBER.                                               C
C**********************************************************************C
C
      IF(KQN.LT.0) THEN
        LQN =-KQN-1
      ELSE
        LQN = KQN
      ENDIF
C
      T1 = 2.0D0*LQN + 1.0D0 + MQN*NSGN
      T2 = 4.0D0*LQN + 2.0D0
      
      CLEBSCH = DSQRT(T1/T2)
C
      RETURN
      END
C
C
      FUNCTION SPHHRM(THETA,PHI,L,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        SSSSSS  PPPPPPP  HH    HH HH    HH RRRRRRR  MM       MM       C
C       SS    SS PP    PP HH    HH HH    HH RR    RR MMM     MMM       C
C       SS       PP    PP HH    HH HH    HH RR    RR MMMM   MMMM       C
C        SSSSSS  PP    PP HHHHHHHH HHHHHHHH RR    RR MM MM MM MM       C
C             SS PPPPPPP  HH    HH HH    HH RRRRRRR  MM  MMM  MM       C
C       SS    SS PP       HH    HH HH    HH RR    RR MM   M   MM       C
C        SSSSSS  PP       HH    HH HH    HH RR    RR MM       MM       C
C                                                                      C
C -------------------------------------------------------------------- C
C  SPHHRM RETURNS A COMPLEX*16 SPHERICAL HARMONIC Y_L^M (THETA,PHI).   C
C -------------------------------------------------------------------- C
C  PARAMETERS:                                                         C
C   INPUT  THETA - ZENITH ANGLE (FROM POSITIVE Z DOWN).                C
C          PHI   - AZIMUTH ANGLE (FROM POSITIVE X TOWARD POSITIVE Y).  C
C          L     - ORBITAL QUANTUM NUMBER (0,1,...).                   C
C          M     - MAGNETIC QUANTUM NUMBER (-L,...,0,...,+L).          C
C  OUTPUT  SPHHRM - DOUBLE COMPLEX NUMBER.                             C
C**********************************************************************C
C
      COMPLEX*16 SPHHRM
C
      DATA PI/3.1415926535897932D0/
C      
C     FOR CONDON-SHORTLEY PHASE CONVENTION
      MTEMP = M
      IF(M.LT.0) THEN
        M  = -M
        IP = (-1)**M
      ENDIF
C
      DEN = 4.0D0*PI*NFACT(L+M)
      FCT = NFACT(L-M)*(2.0D0*L+1.0D0)/DEN
      FCT = DSQRT(FCT)
      AZR = DCOS(M*PHI)
      AZI = DSIN(M*PHI)
      ARG = DCOS(THETA)
      PLM = PLGNDR(L,M,ARG)
C
      IF(M.GE.0) THEN
        SPHHRM = DCMPLX(   FCT*PLM*AZR,    FCT*PLM*AZI)
      ELSEIF(M.LT.0) THEN
        SPHHRM = DCMPLX(IP*FCT*PLM*AZR,-IP*FCT*PLM*AZI)
      ENDIF      
      
      M = MTEMP
C
      RETURN
      END
C
C
      FUNCTION PLGNDR(L,M,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        PPPPPPP  LL       GGGGGG  NN    NN DDDDDDD  RRRRRRR           C
C        PP    PP LL      GG    GG NNN   NN DD    DD RR    RR          C
C        PP    PP LL      GG    GG NNNN  NN DD    DD RR    RR          C
C        PP    PP LL      GG       NN NN NN DD    DD RR    RR          C
C        PPPPPPP  LL      GG   GGG NN  NNNN DD    DD RRRRRRR           C
C        PP       LL      GG    GG NN   NNN DD    DD RR    RR          C
C        PP       LLLLLLLL GGGGGG  NN    NN DDDDDDD  RR    RR          C
C                                                                      C
C -------------------------------------------------------------------- C
C  PLGNDR RETURNS AN ASSOCIATED LEGENDRE POLYNOMIAL P_L^M(X).          C
C -------------------------------------------------------------------- C
C  PARAMETERS:                                                         C
C   INPUT  X - ARGUMENT WITH -1 <= X <= 1.                             C
C          L - ORBITAL QUANTUM NUMBER (0,1,...).                       C
C          M - MAGNETIC QUANTUM NUMBER (0,...,+L).                     C
C  OUTPUT  PLGNDR - ASSOCIATED LEGENDRE POLYNOMIAL P_L^M(X).           C
C**********************************************************************C
C
C     CHECK VALIDITY OF ARGUMENTS
      IF(M.LT.0) THEN
        WRITE(6, *) 'In PLGNDR: bad argument - M < 0.'
        WRITE(7, *) 'In PLGNDR: bad argument - M < 0.'
        STOP
      ELSEIF(M.GT.L) THEN
        WRITE(6, *) 'In PLGNDR: bad argument - L < M: ',L,'<',M
        WRITE(7, *) 'In PLGNDR: bad argument - L < M: ',L,'<',M
        STOP
      ELSEIF(DABS(X).GT.1) THEN
        WRITE(6, *) 'In PLGNDR: bad argument - X NOTIN [-1,+1].'
        WRITE(7, *) 'In PLGNDR: bad argument - X NOTIN [-1,+1].'
        STOP
      ENDIF
C
C     COMPUTE P^M_M (X)
      PMM = 1.0D0
      IF(M.GT.0) THEN
        SOMX2 = DSQRT((1.0D0-X)*(1.0D0+X))
        FACT  = 1.0D0
        DO I=1,M
          PMM  = -PMM*FACT*SOMX2
          FACT =  FACT + 2.0D0
        ENDDO
      ENDIF
C
C     APPLY RECURRENCE RELATION UNTIL WE HAVE P^M_L
      IF(L.EQ.M) THEN
        PLGNDR = PMM
      ELSE
        PMMP1 = X*PMM*(2*M+1)
        IF(L.EQ.M+1) THEN
C         COMPUTE P^M_M+1
          PLGNDR = PMMP1
        ELSE
C         COMPUTE P^M_L, L>M+1
          DO LL=M+2,L
            LLM   = LL-M
            PLL   = (X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/LLM
            PMM   = PMMP1
            PMMP1 = PLL
          ENDDO
          PLGNDR = PLL
        ENDIF
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION NFACT(N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             NN    NN FFFFFFFF   AA     CCCCCC  TTTTTTTT              C
C             NNN   NN FF        AAAA   CC    CC    TT                 C
C             NNNN  NN FF       AA  AA  CC          TT                 C
C             NN NN NN FFFFFF  AA    AA CC          TT                 C
C             NN  NNNN FF      AAAAAAAA CC          TT                 C
C             NN   NNN FF      AA    AA CC    CC    TT                 C
C             NN    NN FF      AA    AA  CCCCCC     TT                 C
C                                                                      C
C -------------------------------------------------------------------- C
C  NFACT(N) RETURNS THE INTEGER RESTULT OF THE FACTORIAL N!.           C
C**********************************************************************C     
      PARAMETER(NMAX=30)
C
C     CHECK THAT ARGUMENT IS VALID
      IF(N.LT.0) THEN
        WRITE(6, *) 'In NFACT: negative input number.',N
        WRITE(7, *) 'In NFACT: negative input number.',N
        STOP
      ELSEIF(N.GT.NMAX) THEN
        WRITE(6,*) 'In NFACT: N too large. N = ',N
        WRITE(7,*) 'In NFACT: N too large. N = ',N
        STOP
      ENDIF
C
C     FACTORIAL AS A PRODUCT OF INTEGERS
      NFACT = 1
      DO M=1,N
        NFACT = M*NFACT
      ENDDO
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C  [11] R-INTS: BOYS INTEGRALS R-INTEGRALS AND RELATED QUANTITIES.     C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] RMAKE: BATCH OF R-INTEGRALS FOR BASIS FUNCTION OVERLAPS.       C
C   [B] WMAKE: BATCH OF W-INTEGRALS FOR BASIS FUNCTION OVERLAPS.       C
C   [C] FUNFM: LIST OF BOYS INTEGRALS FOR USE IN RMAKE.                C
C   [D] BOYSPLOT: OUTPUT DATA FILE WITH FAMILY OF BOYS FUNCTIONS.      C
C   [E] HCOEFF: HERMITE POLYNOMIAL COEFFICIENT FOR USE IN WMAKE.       C
C   [F] HGTF: ARRAY OF HGTFS EVAULATED AT (X,Y,Z).                     C
C   [G] HERMITE: EVALUATION OF H_I (P,X) BY RECURRENCE.                C
C**********************************************************************C

      SUBROUTINE RMAKE(RC,QP,APH,MAXM,LAM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           RRRRRRR  MM       MM    AA    KK    KK EEEEEEEE            C
C           RR    RR MMM     MMM   AAAA   KK   KK  EE                  C
C           RR    RR MMMM   MMMM  AA  AA  KK  KK   EE                  C
C           RR    RR MM MM MM MM AA    AA KKKKK    EEEEEE              C
C           RRRRRRR  MM  MMM  MM AAAAAAAA KK  KK   EE                  C
C           RR    RR MM   M   MM AA    AA KK   KK  EE                  C
C           RR    RR MM       MM AA    AA KK    KK EEEEEEEE            C
C                                                                      C
C -------------------------------------------------------------------- C
C  RMAKE GENERATES A COMPLETE SET OF R-INTEGRALS REQUIRED IN THE       C
C  FINITE SUM REPRESENTATION OF A MULTI-CENTER GAUSSIAN OVERLAP.       C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MLM=26,MKP=9,
     &                      ML4=2*(MKP+1),MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      DIMENSION FS(MB2,MLM),APH(MB2),QP(MB2,3),RC(MB2,MRC),RC2(MB2,MRC)
      DIMENSION F1(MB2,MLM),F2(MB2,MLM),F3(MB2,MLM),F4(MB2,MLM)
      DIMENSION X1(MB2),X2(MB2),X3(MB2),X4(MB2),
     &          I1(MB2),I2(MB2),I3(MB2),I4(MB2)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
C
C
C**********************************************************************C
C     THE FIRST STEP OF THIS ROUTINE IS TO EVALUATE THE REQUIRED       C
C     BOYS INTEGRALS. THIS IS DIVIDED INTO CASES, DEPENDING ON THE     C
C     MAGNITUDE OF THE ARGUMENT.                                       C
C -------------------------------------------------------------------- C
C        FS_M (X) = INT_{0}^{1} T^{2M} EXP(-X*T^{2}) DT                C
C -------------------------------------------------------------------- C
C   EVALUATED FOR ALL VALUES OF M IN THE RANGE 0 < M < LAM.            C
C             FOR ALL VALUES OF X IN THE RANGE X > 0.                  C
C**********************************************************************C
C
      N1 = 0
      N2 = 0
      N3 = 0
      N4 = 0
C
C     FOR EACH PAIR OF BASIS FUNCTIONS (EXPONENTS EI AND EJ IN 'M'),
C     DETERMINE THE BEST WAY TO EVALUATE THE BOYS FUNCTION
      DO M=1,MAXM
C
        X = APH(M)*(QP(M,1)*QP(M,1)+QP(M,2)*QP(M,2)+QP(M,3)*QP(M,3))
C
C       CASE 1: IF X IS ALMOST ZERO (SO WHEN Q=P OR EIJ<<1)
        IF(X.LE.1.0D-11) THEN
          N1     = N1+1
          X1(N1) = X
          I1(N1) = M
C
C       CASE 2: IF X IS SMALLER THAN 17.0D0
        ELSEIF(X.GT.1.0D-11.AND.X.LE.17.0D0) THEN
          N2     = N2+1
          X2(N2) = X
          I2(N2) = M
C
C       CASE 3: IF X IS SMALLER THAN 30.0D0
        ELSEIF(X.GT.17.0D0.AND.X.LE.30.0D0) THEN
          N3     = N3+1
          X3(N3) = X
          I3(N3) = M
C
C       CASE 4: IF X IS LARGER THAN 30.0D0
        ELSE
          N4     = N4+1
          X4(N4) = X
          I4(N4) = M
        ENDIF
C
      ENDDO
C
C     EVALUATE THE BOYS INTEGRALS -- A BATCH FOR EACH ITYPE
C
C     CASE 1: ARGUMENT OF THE BOYS FUNCTION IS X=0.0D0.
C     THE VALUE OF THIS FUNCTION IS 2N+1 (DONE IN FUNFM).
      IF(N1.GT.0) THEN
        CALL FUNFM(F1,X1,N1,LAM,1)
        DO K=1,LAM+1
          DO M=1,N1
            FS(I1(M),K) = F1(M,K)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 2: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=17.0D0.
C     EVALUATE WITH LOCAL POLYNOMIAL EXPANSION OF ORDER 5,
C     AND RECURRENCE IN DIRECTION OF DECREASING M.
      IF(N2.GT.0) THEN
        CALL FUNFM(F2,X2,N2,LAM,2)
        DO K=1,LAM+1
          DO M=1,N2
            FS(I2(M),K) = F2(M,K)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 3: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=30.0D0.
C     EVALUATE USING ASYMPTOTIC FORMULA WITH EXPONENTIAL.
      IF(N3.GT.0) THEN
        CALL FUNFM(F3,X3,N3,LAM,3)
        DO K=1,LAM+1
          DO M=1,N3
            FS(I3(M),K) = F3(M,K)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 4: ARGUMENT OF THE BOYS FUNCTION IS LARGER THAN X=30.0D0.
C     EVALUATE USING ASYMPTOTIC FORMULA WITHOUT EXPONENTIAL.
      IF(N4.GT.0) THEN
        CALL FUNFM(F4,X4,N4,LAM,4)
        DO K=1,LAM+1
          DO M=1,N4
            FS(I4(M),K) = F4(M,K)
          ENDDO
        ENDDO
      ENDIF
C
C**********************************************************************C
C     THE SECOND STEP OF THIS ROUTINE IS TO EVALUATE THE ACTUAL        C
C     R-COEFFICIENTS, BASED ON THE CORRESPONDING BOYS FUNCTIONS.       C
C**********************************************************************C
C
C     CONSTRUCT TOP LEVEL (FOR MAXIMUM LAM VALUE)
      DO M=1,MAXM
        RC(M,1)=((-2.0D0*APH(M))**(LAM))*FS(M,LAM+1)
      ENDDO
C
C     MINIMUM LEVEL ILEV BASED ON LAM VALUE
      IF(MOD(LAM,2).EQ.0) THEN
       ITUVMIN = 1
      ELSE
       ITUVMIN = 2
      ENDIF
C
C     INITIALISE ITUV COUNTER (RELATES TO # CARTESIAN INDICES FOR lam)
      ITUV=-1
C
C     MAIN LOOP: LEVEL 'ILEV' STARTING AT LAM-1 AND WORKING BACKWARDS
      DO ILEV=LAM-1,ITUVMIN,-2
C
C       UPDATE ITUV COUNTER
        ITUV = ITUV+1
C
C       LOOP OVER ALL UNIQUE (IT,IU,IV)
        DO IT=0,ITUV
          RIT = DFLOAT(IT)
C
          DO IU=0,ITUV-IT
            RIU=DFLOAT(IU)
C
            DO IV=0,ITUV-IT-IU
              RIV=DFLOAT(IV)
C
C             R-INTEGRAL CARTESIAN DESTINATION ADDRESSES
              N1 = INABCD(IT+1,IU  ,IV  )
              N2 = INABCD(IT  ,IU+1,IV  )
              N3 = INABCD(IT  ,IU  ,IV+1)
C
C             RECURRENCE RELATIONS DIFFERENT IF ANY (IT,IU,IV) ARE ZERO

C             CASE (IT,IU,IV)
              IF(IT.NE.0.AND.IU.NE.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE (IT,IU, 0)
              ELSEIF(IT.NE.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
               ENDDO
C
C             CASE (IT, 0,IV)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE (IT, 0, 0)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE ( 0,IU,IV)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE ( 0,IU, 0)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE ( 0, 0,IV)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE ( 0, 0, 0)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
                ENDDO
C
C             ALL POSSIBILITIES ACCOUNTED FOR -- END LOOP
              ENDIF
C
C           END LOOPS OVER (IT,IU,IV) ADDRESSES FOR THIS ILEV
            ENDDO
          ENDDO
        ENDDO
C
C       ADD IN (IT=0,IU=0,IV=0) CASE
        DO M=1,MAXM
          RC2(M,1) = ((-2.0D0*APH(M))**(ILEV))*FS(M,ILEV+1)
        ENDDO
C
C
C       UPDATE ITUV COUNTER
        ITUV = ITUV+1
C
C       LOOP OVER ALL UNIQUE (IT,IU,IV)
        DO IT=0,ITUV
          RIT=DFLOAT(IT)
          DO IU=0,ITUV-IT
            RIU=DFLOAT(IU)
            DO IV=0,ITUV-IT-IU
              RIV=DFLOAT(IV)
C
C             R-INTEGRAL CARTESIAN DESTINATION ADDRESSES
              N1 = INABCD(IT+1,IU  ,IV  )
              N2 = INABCD(IT  ,IU+1,IV  )
              N3 = INABCD(IT  ,IU  ,IV+1)
C
C             RECURRENCE RELATIONS DIFFERENT IF ANY (IT,IU,IV) ARE ZERO

C             CASE (IT,IU,IV)
              IF(IT.NE.0.AND.IU.NE.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
C
C             CASE (IT,IU, 0)
              ELSEIF(IT.NE.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1)
                ENDDO
C
C             CASE (IT, 0,IV)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
C
C             CASE (IT, 0, 0)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                DO M=1,MAXM
                 RC(M,N1) = -QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                 RC(M,N2) = -QP(M,2)*RC2(M,K1)
                 RC(M,N3) = -QP(M,3)*RC2(M,K1)
               ENDDO
C
C             CASE ( 0,IU,IV)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
C
C             CASE ( 0,IU, 0)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1)
                ENDDO
C
C             CASE ( 0, 0,IV)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
C
C             CASE ( 0, 0, 0)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                DO M=1,MAXM
                  RC(M,N1) = -QP(M,1)*RC2(M,K1)
                  RC(M,N2) = -QP(M,2)*RC2(M,K1)
                  RC(M,N3) = -QP(M,3)*RC2(M,K1)
                ENDDO
C
C             ALL POSSIBILITIES ACCOUNTED FOR -- END LOOP
              ENDIF
C
C           END LOOPS OVER (IT,IU,IV) ADDRESSES FOR THIS ILEV
            ENDDO
          ENDDO
        ENDDO
C
C       ADD IN (IT=0,IU=0,IV=0) CASE
C
        DO M=1,MAXM
          RC(M,1) = ((-2.0D0*APH(M))**(ILEV-1))*FS(M,ILEV)
        ENDDO
C
      ENDDO
C
C
C     AN ADDITIONAL LOOP OVER ADDRESSES (WHEN LAM IS ODD)
      IF(MOD(LAM,2).EQ.1) THEN
C
C       UPDATE ITUV COUNTER
        ITUV = ITUV+1
C
C       LOOP OVER ALL UNIQUE (IT,IU,IV)
        DO IT=0,ITUV
          RIT=DFLOAT(IT)
          DO IU=0,ITUV-IT
            RIU=DFLOAT(IU)
            DO IV=0,ITUV-IT-IU
              RIV=DFLOAT(IV)
C
C             R-INTEGRAL CARTESIAN DESTINATION ADDRESSES
              N1 = INABCD(IT+1,IU  ,IV  )
              N2 = INABCD(IT  ,IU+1,IV  )
              N3 = INABCD(IT  ,IU  ,IV+1)
C
C             RECURRENCE RELATIONS DIFFERENT IF ANY (IT,IU,IV) ARE ZERO

C             CASE (IT,IU,IV)
              IF(IT.NE.0.AND.IU.NE.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE (IT,IU, 0)
              ELSEIF(IT.NE.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE (IT, 0,IV)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE (IT, 0, 0)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M1 = INABCD(IT-1,IU  ,IV  )
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE ( 0,IU,IV)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE ( 0,IU, 0)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M2 = INABCD(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE ( 0, 0,IV)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                M3 = INABCD(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE ( 0, 0, 0)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
                K1 = INABCD(IT  ,IU  ,IV  )
                DO M=1,MAXM
                  RC2(M,N1) = -QP(M,1)*RC(M,K1)
                  RC2(M,N2) = -QP(M,2)*RC(M,K1)
                  RC2(M,N3) = -QP(M,3)*RC(M,K1)
               ENDDO
C
C             ALL POSSIBILITIES ACCOUNTED FOR -- END LOOP
              ENDIF
C
C           END LOOPS OVER (IT,IU,IV) ADDRESSES FOR THIS ILEV
            ENDDO
          ENDDO
        ENDDO
C
C       ADD IN (IT=0,IU=0,IV=0) CASE
C
        DO M=1,MAXM
          RC2(M,1) = FS(M,1)
        ENDDO
C
C       MOVE THE RC2 ARRAY INTO RC
C
        ITMAX = (LAM+1)*(LAM+2)*(LAM+3)/6
        DO IT=1,ITMAX
          DO M=1,MAXM
            RC(M,IT) = RC2(M,IT)
          ENDDO
        ENDDO
C
C     END IF STATEMENT FOR THE ODD LAM CASE
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE WMAKE(WC,RP,APH,MAXM,LAM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         WW         WW MM       MM    AA    KK    KK EEEEEEEE         C
C         WW         WW MMM     MMM   AAAA   KK   KK  EE               C
C         WW         WW MMMM   MMMM  AA  AA  KK  KK   EE               C
C         WW    W    WW MM MM MM MM AA    AA KKKKK    EEEEEE           C
C          WW  WWW  WW  MM  MMM  MM AAAAAAAA KK  KK   EE               C
C           WWWW WWWW   MM   M   MM AA    AA KK   KK  EE               C
C            WW   WW    MM       MM AA    AA KK    KK EEEEEEEE         C
C                                                                      C
C -------------------------------------------------------------------- C
C  FOR A PARTICULAR BLOCK OF BASIS FUNCTION OVERLAPS, MAKE THE W       C
C  FUNCTIONS W(APH,RP;A,B,C) AS A FUNCTION OF (X,Y,Z) -- THIS          C
C  COMES TO THE POTENTIAL FROM A HGTF OVERLAP CHARGE SOURCE.           C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MLM=30,
     &                      ML4=2*(MKP+1),MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      DIMENSION FS(MB2,MLM),APH(MB2),RP(MB2,3),WC(MB2,MRC)
      DIMENSION X1(MB2),X2(MB2),X3(MB2),X4(MB2),
     &          I1(MB2),I2(MB2),I3(MB2),I4(MB2)
      DIMENSION F1(MB2,MLM),F2(MB2,MLM),F3(MB2,MLM),F4(MB2,MLM)     
      DIMENSION H(0:MKP-1,0:(MKP-1)/2,3)
C
      COMMON/ACCSS/INABCD(0:ML4,0:ML4,0:ML4),
     &             IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
C
      DATA PI/3.141592653589793D0/      
C
C**********************************************************************C
C     THE FIRST STEP OF THIS ROUTINE IS TO EVALUATE THE REQUIRED       C
C     BOYS INTEGRALS. THIS IS DIVIDED INTO CASES, DEPENDING ON THE     C
C     MAGNITUDE OF THE ARGUMENT.                                       C
C -------------------------------------------------------------------- C
C        FS_M (X) = INT_{0}^{1} T^{2M} EXP(-X*T^{2}) DT                C
C -------------------------------------------------------------------- C
C   EVALUATED FOR ALL VALUES OF M IN THE RANGE 0 < M < LAM.            C
C             FOR ALL VALUES OF X IN THE RANGE X > 0.                  C
C**********************************************************************C
C
C     CALCULATE GAMMA FUNCTION VALUES FOR LATER USE
      CALL GAMGEN
C
      N1 = 0
      N2 = 0
      N3 = 0
      N4 = 0
C
C     FOR EACH PAIR OF BASIS FUNCTIONS (EXPONENTS EI AND EJ IN 'M'),
C     DETERMINE THE BEST WAY TO EVALUATE THE BOYS FUNCTION
      DO M=1,MAXM
C       X = EIJ*(P-C)^2
        X = APH(M)*(RP(M,1)*RP(M,1)+RP(M,2)*RP(M,2)+RP(M,3)*RP(M,3))
C       CASE 1: IF X IS ALMOST ZERO (SO WHEN Q=P OR EIJ<<1)
        IF(X.LE.1.00D-11) THEN
          N1     = N1 + 1
          X1(N1) = X
          I1(N1) = M
C       CASE 2: IF X IS SMALLER THAN X=17.0D0
        ELSEIF(X.GT.1.00D-11.AND.X.LE.1.70D+01) THEN
          N2     = N2 + 1
          X2(N2) = X
          I2(N2) = M
C       CASE 3: IF X IS SMALLER THAN ABOUT 30.0D0
        ELSEIF(X.GT.1.70D+01.AND.X.LE.3.00D+01) THEN
          N3     = N3 + 1
          X3(N3) = X
          I3(N3) = M
C       CASE 4: IF X IS LARGER THAN ABOUT 30.0D0
        ELSE
          N4     = N4 + 1
          X4(N4) = X
          I4(N4) = M
        ENDIF
      ENDDO
C
C     CASE 1: ARGUMENT OF THE BOYS FUNCTION IS X=0.0D0.
C     THE VALUE OF THIS FUNCTION IS 2N+1 (DONE IN FUNFM).
      IF(N1.NE.0) THEN
        CALL FUNFM(F1,X1,N1,LAM,1)     
        DO JJ=1,LAM+1
          DO M=1,N1
            FS(I1(M),JJ) = F1(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 2: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=17.0D0.
C     EVALUATE WITH LOCAL POLYNOMIAL EXPANSION OF ORDER 5,
C     AND RECURRENCE IN DIRECTION OF DECREASING M.
      IF(N2.NE.0) THEN
        CALL FUNFM(F2,X2,N2,LAM,2)
        DO JJ=1,LAM+1
          DO M=1,N2
            FS(I2(M),JJ) = F2(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 3: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=30.0D0.
C     EVALUATE USING ASYMPTOTIC FORMULA WITH EXPONENTIAL.
      IF(N3.NE.0) THEN
        CALL FUNFM(F3,X3,N3,LAM,3)
        DO JJ=1,LAM+1
          DO M=1,N3
            FS(I3(M),JJ) = F3(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C     CASE 4: ARGUMENT OF THE BOYS FUNCTION IS LARGER THAN X=30.0D0.
C     EVALUATE USING ASYMPTOTIC FORMULA WITHOUT EXPONENTIAL.
      IF(N4.NE.0) THEN
        CALL FUNFM(F4,X4,N4,LAM,4)
        DO JJ=1,LAM+1
          DO M=1,N4
            FS(I4(M),JJ) = F4(M,JJ)
          ENDDO
        ENDDO
      ENDIF
C
C**********************************************************************C
C     THE SECOND STEP OF THIS SUBROUTINE IS TO EVALUATE THE ACTUAL     C
C     W-FUNCTIONS, BASED ON THE CORRESPONDING BOYS INTEGRALS.          C
C**********************************************************************C
C
C     INITIALISE THE WC ARRAY
      DO M=1,MAXM
        DO ITUV=1,LAM
          WC(M,ITUV) = 0.0D0
        ENDDO
      ENDDO
C
C     NUMBER OF ADDRESSES IADD FOR LAM VALUE
      NTUV = ((LAM+1)*(LAM+2)*(LAM+3))/6
C
C     LOOP OVER BASIS FUNCTION PAIRS
      DO M=1,MAXM
C       
C       GENERATE ALL NECESSARY HERMITE POLYNOMIALS FOR LATER REFERENCE
        DO ILAM=0,LAM
C         TOTAL NUMBER OF TERMS IN HERMITE EXPANSION
          NTERMS = (ILAM-MOD(ILAM,2))/2
          DO N=0,NTERMS
            H(ILAM,N,1) = HCOEFF(APH(M),RP(M,1),ILAM,N)
            H(ILAM,N,2) = HCOEFF(APH(M),RP(M,2),ILAM,N)
            H(ILAM,N,3) = HCOEFF(APH(M),RP(M,3),ILAM,N)
          ENDDO
        ENDDO
C
C       START A LOOP THAT RUNS OVER ALL COMBINATIONS ITUV={T,U,V}
        DO ITUV=1,NTUV
C
C         ESTABLISH PARAMETERS T,U,V AND LAM FOR THIS VALUE ITUV
          IT   = IVEC(ITUV)
          IU   = JVEC(ITUV)
          IV   = KVEC(ITUV)
          ILAM = LAMVEC(ITUV)
C
C         FOR GIVEN ITUV, START THREE NESTED LOOPS FOR X,Y,Z DECOMP.
          STORE = 0.0D0
          DO IA=0,(IT-MOD(IT,2))/2
            IPOW = (IT + MOD(IT,2))/2 + IA
            DO IB=0,(IU-MOD(IU,2))/2
              JPOW = (IU + MOD(IU,2))/2 + IB
              DO IC=0,(IV-MOD(IV,2))/2
                KPOW = (IV + MOD(IV,2))/2 + IC
C               
C               PRODUCT OF HERMITE EXPANSION COEFFICIENTS
                HTMS = H(IT,IA,1)*H(IU,IB,2)*H(IV,IC,3)
C
C               BOYS FUNCTION OF RELEVANT INDEX NBOYS
                NBOYS = IPOW + JPOW + KPOW
                FBOYS = FS(M,NBOYS+1)
C               ADD TO W FUNCTION BIN THE RELEVANT PRODUCT
                STORE = STORE + HTMS*FBOYS
C
              ENDDO
            ENDDO
          ENDDO
C
C         MULTIPLY RESULT BY 2*PI/LAM
          WC(M,ITUV) = 2.0D0*PI*STORE/APH(M)
C
        ENDDO
C     
      ENDDO
C      
      RETURN
      END
C
C
      SUBROUTINE FUNFM(FM,T,N,LAM,ITYPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           FFFFFFFF UU    UU NN    NN FFFFFFFF MM       MM            C
C           FF       UU    UU NNN   NN FF       MMM     MMM            C
C           FF       UU    UU NNNN  NN FF       MMMM   MMMM            C
C           FFFFFF   UU    UU NN NN NN FFFFFF   MM MM MM MM            C
C           FF       UU    UU NN  NNNN FF       MM  MMM  MM            C
C           FF       UU    UU NN   NNN FF       MM   M   MM            C
C           FF        UUUUUU  NN    NN FF       MM       MM            C
C                                                                      C
C -------------------------------------------------------------------- C
C  FUNFM EVALUATES INTEGRAL [INT_{0}^{1} U^{2M} EXP(-T*U^{2}) dU]      C
C  FOR VARIABLE T > 0 FOR ALL ORDERS 0 < M < LAM.                      C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C    ITYPE = 1 - SPECIAL CASE X = 0.0D0.                               C
C    ITYPE = 2 - POWER SERIES AND REVERSE RECURRENCE.                  C
C                (ONLY MSER TERMS WILL BE USED, SO USE MUST SUPPLY A   C
C                 VALUE APPROPRIATE TO THE MAX VALUE OF X IN BATCH).   C
C    ITYPE = 3 - ASYMPTOTIC EXPANSION AND FORWARD RECURRENCE.          C
C    ITYPE = 4 - ASYMPTOTIC EXPANSION AND FORWARD RECURRENCE.          C
C                ALL TERMS DEPENDING ON EXP(-X) ARE OMITTED TO AVOID   C
C                NUMERICAL UNDERFLOW PROBLEMS. MSER NOT REQUIRED.      C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MLM=26,MSER=60)
C
      DIMENSION FM(MB2,MLM),T(MB2),TLAM(MB2),TT2(MB2),
     &          TXP(MB2),TRT(MB2)
C
      DATA PIROOT,A0,B0/8.862269254527580D-1,4.994501191201870D-1,
     &                                       4.551838436668326D-1/
C
C**********************************************************************C
C     ITYPE = 1: SPECIAL CASE FOR T = 0.0D0                            C
C**********************************************************************C
C
      IF(ITYPE.EQ.1) THEN
        DO K=1,LAM+1
          MVAL  = K-1
          VALUE = 1.0D0/DFLOAT(MVAL+MVAL+1)
          DO M=1,N
            FM(M,K) = VALUE
          ENDDO
        ENDDO
        RETURN
C
C**********************************************************************C
C     ITYPE = 2: POWER SERIES EVALUATION                               C
C**********************************************************************C
C
      ELSEIF(ITYPE.EQ.2) THEN
C
C       INITIALIZE THE POWER SERIES FOR M = LAM
        DO M=1,N
          TXP(M)      = DEXP(-T(M))
          TT2(M)      = 2.0D0*T(M)
          TLAM(M)     = 1.0D0
          FM(M,LAM+1) = 1.0D0
        ENDDO
C
C       LOOP OVER TERMS IN THE POWER SERIES
        DO K=1,MSER
          DLAM = DFLOAT(2*(LAM+K)+1)
          DO M=1,N
            TLAM(M)     = TLAM(M)*(TT2(M)/DLAM)
            FM(M,LAM+1) = FM(M,LAM+1) + TLAM(M)
          ENDDO
        ENDDO
C
C       RESCALE BY THE PREFACTOR
        DEN = DFLOAT((2*LAM)+1)
        DO M=1,N
          FM(M,LAM+1) = FM(M,LAM+1)*TXP(M)/DEN
        ENDDO
C
C       NOW COMPLETE TABLE BY BACKWARDS RECURRENCE
        DO I=1,LAM
          MIND  = LAM-I+1
          MVAL  = MIND-1
          COEFF = DFLOAT(MVAL+MVAL+1)
          DO M=1,N
            FM(M,MIND) = (TT2(M)*FM(M,MIND+1) + TXP(M))/COEFF
          ENDDO
        ENDDO
C
C**********************************************************************C
C     ITYPE = 3: ASYMPTOTIC EXPANSION WITH VERY LARGE ARGUMENT.        C
C**********************************************************************C
C
      ELSEIF(ITYPE.EQ.3) THEN
C
C       INITIALIZE THE ASYMPTOTIC EXPANSION
        DO M=1,N
          TXP(M) = DEXP(-T(M))
          TT2(M) = 2.0D0*T(M)
          TRT(M) = DSQRT(T(M)) 
        ENDDO
C
C       SEED VALUES
        DO M=1,N
          FM(M,1) = A0/(B0+T(M))
        ENDDO
C
C       RESCALE BY THE PREFACTOR
        DO M=1,N
          FM(M,1) = (PIROOT/TRT(M)) - (TXP(M)*FM(M,1))
        ENDDO
C
C       NOW COMPLETE TABLE BY FORWARD RECURRENCE
        DO MIND=1,LAM
          MVAL  = MIND-1
          COEFF = DFLOAT(MVAL+MVAL+1)
          DO M=1,N
            FM(M,MIND+1) = (COEFF*FM(M,MIND) - TXP(M))/TT2(M)
          ENDDO
        ENDDO
C
C**********************************************************************C
C     ITYPE = 4: ASYMPTOTIC EXPANSION WITH VERY LARGE ARGUMENT         C
C**********************************************************************C
C
      ELSEIF(ITYPE.EQ.4) THEN
C
C     INITIALIZE THE ASYMPTOTIC EXPANSION
      DO M=1,N
        TT2(M)  = 2.0D0*T(M)
        FM(M,1) = PIROOT/DSQRT(T(M))
      ENDDO
C
C     NOW COMPLETE TABLE BY FORWARD RECURRENCE
      DO MIND=1,LAM
        MVAL  = MIND-1
        COEFF = DFLOAT(MVAL+MVAL+1)
        DO M=1,N
          FM(M,MIND+1) = (COEFF*FM(M,MIND))/TT2(M)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     ITYPE OUT OF RANGE: INVALID INPUT TO FUNFM                       C
C**********************************************************************C
C
      ELSE
91      FORMAT(2X,'In FUNFM: invalid type (must be 1-4)',I4)
        WRITE(6,91) ITYPE
        WRITE(7,91) ITYPE
        STOP
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE BOYSGEN(LAMBDA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    BBBBBBB   OOOOOO  YY    YY  SSSSSS   GGGGGG  EEEEEEEE NN    NN    C
C    BB    BB OO    OO YY    YY SS    SS GG    GG EE       NNN   NN    C
C    BB    BB OO    OO  YY  YY  SS       GG       EE       NNNN  NN    C
C    BBBBBBB  OO    OO   YYYY    SSSSSS  GG       EEEEEE   NN NN NN    C
C    BB    BB OO    OO    YY          SS GG   GGG EE       NN  NNNN    C
C    BB    BB OO    OO    YY    SS    SS GG    GG EE       NN   NNN    C
C    BBBBBBB   OOOOOO     YY     SSSSSS   GGGGGG  EEEEEEEE NN    NN    C
C                                                                      C
C -------------------------------------------------------------------- C
C  BOYSGEN PRODUCES A DATA FILE WHICH CONTAINS A FAMILY OF BOYS        C
C  FUNCTIONS OVER A SPECIFIED REGION, GIVEN A MAXIMUM FAMILY           C
C  PARAMETER DETERMINED BY LAMBDA.                                     C
C**********************************************************************C
      PARAMETER (MBS=26,MB2=MBS*MBS,MLM=30)
C
      DIMENSION X1(MB2),X2(MB2),X3(MB2),X4(MB2)
      DIMENSION F1(MB2,MLM),F2(MB2,MLM),F3(MB2,MLM),
     &          F4(MB2,MLM),FS(MB2,MLM)
C
C     EVALUATION PARAMETERS
      NTOT = 400
      XMIN = 0.00D0
      XMAX = 4.00D1
      HSTP = (XMAX-XMIN)/NTOT
C
C**********************************************************************C
C     THE FIRST STEP OF THIS SUBROUTINE IS TO EVALUATE THE REQUIRED    C
C     BOYS INTEGRALS. THIS IS DIVIDED INTO CASES, DEPENDING ON THE     C
C     MAGNITUDE OF THE ARGUMENT.                                       C
C -------------------------------------------------------------------- C
C        FS_M (X) = INT_{0}^{1} T^{2M} EXP(-X*T^{2}) DT                C
C -------------------------------------------------------------------- C
C   EVALUATED FOR ALL VALUES OF M IN THE RANGE 0 < M < LAMBDA.         C
C             FOR ALL VALUES OF X IN THE RANGE X > 0.                  C
C**********************************************************************C
C
      OPEN(UNIT=8,FILE='plots/boysfunction.dat',STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO NX=0,NTOT
C
      X = XMIN + HSTP*NX
C      
C     DETERMINE THE BEST WAY TO EVALUATE THE BOYS FUNCTION
      IF(X.LE.1.00D-11) THEN
C     CASE 0: ARGUMENT OF THE BOYS FUNCTION IS X = 0.
C             THE VALUE OF THIS FUNCTION IS (2N+1).
        N1     = 1
        X1(N1) = X
        CALL FUNFM(F1,X1,N1,LAMBDA,1)
        DO JJ=1,LAMBDA+1
          FS(N1,JJ) = F1(N1,JJ)
        ENDDO
      ELSEIF(X.GT.1.00D-11.AND.X.LE.1.70D+01) THEN
C     CASE 1: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=17.
C             EVALUATE WITH LOCAL POLYNOMIAL EXPANSION OF ORDER 5,
C             AND RECURRENCE IN DIRECTION OF DECREASING M.
        N2     = 1
        X2(N2) = X
        CALL FUNFM(F2,X2,N2,LAMBDA,2)
        DO JJ=1,LAMBDA+1
          FS(N2,JJ) = F2(N2,JJ)
        ENDDO
      ELSEIF(X.GT.1.70D+01.AND.X.LE.3.00D+01) THEN
C     CASE 2: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=30.
C             EVALUATE USING ASYMPTOTIC FORMULA WITH EXPONENTIAL.
        N3     = 1
        X3(N3) = X
        CALL FUNFM(F3,X3,N3,LAMBDA,3)
        DO JJ=1,LAMBDA+1
          FS(N3,JJ) = F3(N3,JJ)
        ENDDO
      ELSE
C     CASE 3: ARGUMENT OF THE BOYS FUNCTION IS LARGER THAN X=30.
C             EVALUATE USING ASYMPTOTIC FORMULA WITHOUT EXPONENTIAL.
        N4     = 1
        X4(N4) = X
        CALL FUNFM(F4,X4,N4,LAMBDA,4)
        DO JJ=1,LAMBDA+1
          FS(N4,JJ) = F4(N4,JJ)
        ENDDO
      ENDIF
C      
      WRITE(8,*) X,(FS(1,L),L=1,LAMBDA+1)
C      
      ENDDO
      CLOSE(UNIT=8)
C      
      RETURN
      END
C
C
      FUNCTION HCOEFF(P,X,LAM,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        HH    HH  CCCCCC   OOOOOO  EEEEEEEE FFFFFFFF FFFFFFFF         C
C        HH    HH CC    CC OO    OO EE       FF       FF               C
C        HH    HH CC       OO    OO EE       FF       FF               C
C        HHHHHHHH CC       OO    OO EEEEEE   FFFFFF   FFFFFF           C
C        HH    HH CC       OO    OO EE       FF       FF               C
C        HH    HH CC    CC OO    OO EE       FF       FF               C
C        HH    HH  CCCCCC   OOOOOO  EEEEEEEE FF       FF               C
C                                                                      C
C -------------------------------------------------------------------- C
C  HCOEFF DETERMINES THE NTH COEFFICIENT OF A HERMITE POLYNOMIAL       C
C  EXPANSION FOR H_LAM (P,X) FOR USE BY THE W FUNCTIONS, WHICH IN      C
C  TURN GENERATES THE POTENTIAL OF A WAVE FUNCTION.                    C
C**********************************************************************C
      PARAMETER (LAMMAX=40)
C
      IF(LAM.LT.0) THEN
        WRITE(6, *) 'In HCOEFF: index less than zero. LAM = ',LAM
        WRITE(7, *) 'In HCOEFF: index less than zero. LAM = ',LAM
        STOP
      ENDIF
C
      IF(LAM.GT.LAMMAX) THEN
        WRITE(6, *) 'In HCOEFF: index greater than max. LAM = ',LAM
        WRITE(6, *) 'In HCOEFF: index greater than max. LAM = ',LAM
        STOP
      ENDIF      
C
C     EVALUATE SOME QUANTITIES COMMON TO BOTH CLASSES
      NPRE  = NFACT(LAM)
      NPOW  = (LAM - MOD(LAM,2))/2 - N
      TOP   = (-P)**(NPOW)
      NDEN2 = NFACT(NPOW)
      NOTHR = 2*N + MOD(LAM,2)
      NDEN1 = NFACT(NOTHR)
      TRM   = 2.0D0*P*X
      POLY  = TRM**NOTHR
C
      ND12   = NDEN1*NDEN2
      HCOEFF = NPRE*TOP*POLY/ND12
C
      RETURN
      END
C
C
      SUBROUTINE HGTF(HABC,XYZEVAL,XYZ,EXL,KAPPA,NFUNS,I1,I2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 HH    HH  GGGGGG  TTTTTTTT FFFFFFFF                  C
C                 HH    HH GG    GG    TT    FF                        C
C                 HH    HH GG          TT    FF                        C
C                 HHHHHHHH GG          TT    FFFFFF                    C
C                 HH    HH GG   GGG    TT    FF                        C
C                 HH    HH GG    GG    TT    FF                        C
C                 HH    HH  GGGGGG     TT    FF                        C
C                                                                      C
C -------------------------------------------------------------------- C
C  HGTF GENERATES AN ARRAY OF HGTF FUNCTIONS EVALUATED AT A SET OF     C
C  COORDINATES XYZEVAL(3), USING THE GAUSSIAN PRODUCT THEOREM TO       C
C  TRANSFORM AN ARBITRARY PRODUCT OF TWO BASIS FUNCTIONS INTO A        C
C  FINITE SUM OF EQ-COEFFICIENTS AND HGTFS ON A SINGLE CENTER, RP.     C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=2*ML2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &             IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
C      
      DIMENSION HABC(MB2,MEQ)
      DIMENSION RLOC(3,2),EXL(MBS,4),KAPPA(4),NFUNS(4)
      DIMENSION XYZEVAL(3),XYZ(3,4)
C
C     TRANSFER KAPPA TO LOCAL LQNS
      IF(KAPPA(I1).LT.0) THEN
        LQN1 =-KAPPA(I1)-1
      ELSE
        LQN1 = KAPPA(I1) 
      ENDIF
C      
      IF(KAPPA(I2).LT.0) THEN
        LQN2 =-KAPPA(I2)-1
      ELSE
        LQN2 = KAPPA(I2)
      ENDIF
C
C     DETERMINE THE NUMBER OF FUNCTIONS ON EACH CENTER
      NFUNA = NFUNS(I1)
      NFUNB = NFUNS(I2)
      MAXM  = NFUNS(I1)*NFUNS(I2)
C
C     ASSUME MAXIMAL CASE (NECESSARY FOR SS OVERLAP)
      LAMBDA = LQN1 + LQN2 + 2
      NTUV   = ((LAMBDA+1)*(LAMBDA+2)*(LAMBDA+3))/6
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NFUNA
        EXPA = EXL(IBAS,I1)
        DO JBAS=1,NFUNB
          M    = M + 1
          EXPB = EXL(JBAS,I2)
C         SUM OF EXPONENTS
          PAB  = EXPA + EXPB
C         HGTF CENTER COORDINATES
          PX = (EXPA*XYZ(1,I1) + EXPB*XYZ(1,I2))/PAB
          PY = (EXPA*XYZ(2,I1) + EXPB*XYZ(2,I2))/PAB
          PZ = (EXPA*XYZ(3,I1) + EXPB*XYZ(3,I2))/PAB
C         COORDINATES WRT LOCAL ORIGIN
          RPX = XYZEVAL(1) - PX
          RPY = XYZEVAL(2) - PY
          RPZ = XYZEVAL(3) - PZ
C         GAUSSIAN COMPONENT OF HGTF
          GSS  = DEXP(-PAB*(RPX*RPX + RPY*RPY + RPZ*RPZ))
          DO ITUV=1,NTUV
C           HERMITE POLYNOMIALS FOR THIS BASIS PAIR OVER ALL ITUV
            HALPH = HERMITE(PAB,RPX,IVEC(ITUV))
            HBETA = HERMITE(PAB,RPY,JVEC(ITUV))
            HGAMA = HERMITE(PAB,RPZ,KVEC(ITUV))
C           ADD TO BASIS FUNCTION PRODUCT
            HABC(M,ITUV) = HALPH*HBETA*HGAMA*GSS
          ENDDO          
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION HERMITE(P,X,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    HH    HH EEEEEEEE RRRRRRR  MM       MM IIII TTTTTTTT EEEEEEEE     C
C    HH    HH EE       RR    RR MMM     MMM  II     TT    EE           C
C    HH    HH EE       RR    RR MMMM   MMMM  II     TT    EE           C
C    HHHHHHHH EEEEEE   RR    RR MM MM MM MM  II     TT    EEEEEE       C
C    HH    HH EE       RRRRRRR  MM  MMM  MM  II     TT    EE           C
C    HH    HH EE       RR    RR MM   M   MM  II     TT    EE           C
C    HH    HH EEEEEEEE RR    RR MM       MM IIII    TT    EEEEEEEE     C
C                                                                      C
C -------------------------------------------------------------------- C
C  HERMITE EVALUATES HERMITE POLYNOMIAL H_I (P,X) BY RECURRENCE.       C
C**********************************************************************C
      PARAMETER (IMAX=40)
C
      IF(I.LT.0) THEN
        WRITE(6, *) 'In HERMITE: index less than zero. I = ',I
        WRITE(7, *) 'In HERMITE: index less than zero. I = ',I
        STOP
      ENDIF
C
      IF(I.GT.IMAX) THEN
        WRITE(6, *) 'In HERMITE: index greater than max. I = ',I
        WRITE(7, *) 'In HERMITE: index greater than max. I = ',I
        STOP
      ENDIF      
C
C     NEED FIRST TWO VALUES TO ESTABLISH RECURRENCE RELATION      
      TEMP1 = 1.0D0
      TEMP2 = 2.0D0*P*X
C
      IF(I.EQ.0) THEN
        HERMITE = TEMP1
      ELSEIF(I.EQ.1) THEN
        HERMITE = TEMP2
      ELSEIF(I.GT.1) THEN
        DO N=2,I
          TEMP3 = 2.0D0*P*(X*TEMP2 - (N-1)*TEMP1)
          TEMP1 = TEMP2
          TEMP2 = TEMP3
        ENDDO
        HERMITE = TEMP2
      ENDIF
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C  [12] E-COEFFS: FINITE BASIS OVERLAP SPIN STRUCTURE FACTORS.         C
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
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
      COMMON/TIME/TATM,TSCF,TMPT,TMCF,TDMG,TPRP,TPLT,TTOT,T1EL,T2CL,
     &            T2BR,TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR,TBEG,TEND
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
C     LOOP OVER CENTERS A AND B
      DO ICNTA=1,NCNT
        DO ICNTB=1,NCNT
C
C         LOOP OVER KQNA VALUES
          DO KA=1,NKAP(ICNTA)
C
C           QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
            IF(KVALS(KA,ICNTA).GT.0) THEN
              LQNA = KVALS(KA,ICNTA)
            ELSE
              LQNA =-KVALS(KA,ICNTA)-1
            ENDIF
            NFUNA = NFUNCT(LQNA+1,ICNTA)
C
C           LOOP OVER KQNB VALUES
            DO KB=1,NKAP(ICNTB)
C
C             QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
              IF(KVALS(KB,ICNTB).GT.0) THEN
                LQNB = KVALS(KB,ICNTB)
              ELSE
                LQNB =-KVALS(KB,ICNTB)-1
              ENDIF
              NFUNB = NFUNCT(LQNB+1,ICNTB)
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
                  MAXAB = NFUNA*NFUNB
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
        WRITE(6,*) 'In EQFILE: E0LL words exceed allocated limit.'
        WRITE(7,*) 'In EQFILE: E0LL words exceed allocated limit.'
        GOTO 150
      ENDIF
C
      IF(HMLTN.NE.'NORL') THEN
        IF(NADE0SS.GT.MFL) THEN
          WRITE(6,*) 'In EQFILE: E0SS words exceed allocated limit.'
          WRITE(7,*) 'In EQFILE: E0SS words exceed allocated limit.'
          GOTO 150
        ENDIF
      ENDIF
C
      IF(HMLTN.EQ.'DHFB'.OR.HMLTN.EQ.'DHFP'.OR.HMLTN.EQ.'DHFQ') THEN
        IF(NADEILS.GT.MFL) THEN
          WRITE(6,*) 'In EQFILE: EILS words exceed allocated limit.'
          WRITE(7,*) 'In EQFILE: EILS words exceed allocated limit.'
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
      WRITE(6,*) 'In EQFILE: E-coefficients will be generated by batch.'
      WRITE(7,*) 'In EQFILE: E-coefficients will be generated by batch.'
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
C        EEEEEEEE  SSSSSS  EEEEEEEE TTTTTTTT LL       LL               C
C        EE       SS    SS EE          TT    LL       LL               C
C        EE       SS       EE          TT    LL       LL               C
C        EEEEEE    SSSSSS  EEEEEE      TT    LL       LL               C
C        EE             SS EE          TT    LL       LL               C
C        EE       SS    SS EE          TT    LL       LL               C
C        EEEEEEEE  SSSSSS  EEEEEEEE    TT    LLLLLLLL LLLLLLLL         C
C                                                                      C
C -------------------------------------------------------------------- C
C  ESETLL CONSTRUCTS ALL ELL0 COEFFICIENTS FOR A SYSTEM WITH CALLS TO  C
C  EMAKELL AND SAVES THEM TO AN EXTERNAL DATA FILE.                    C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTERS A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     INITIALISE COUNT INDEX
      ICOUNT = 0
C
C     LOOP OVER CENTER A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 0
        ENDIF
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNS(1) = NFUNCT(LQN(1)+1,ICNTA)
        DO IBAS=1,NFUNS(1)
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNS(2) = NFUNCT(LQN(2)+1,ICNTB)
        DO JBAS=1,NFUNS(2)
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
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS           C
C     OF (MA,MB). FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE          C
C     ORDERED                                                          C
C     11: = (-|MQN(A)|,-|MQN(B)|)                                      C
C     12: = (-|MQN(A)|,+|MQN(B)|)                                      C
C     21: = (+|MQN(A)|,-|MQN(B)|)                                      C
C     22: = (+|MQN(A)|,+|MQN(B)|)                                      C
C**********************************************************************C
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NFUNS(1)*NFUNS(2)
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMAB  = LQN(1)+LQN(2)
      NTUVLL = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
C
C     INDEX TO TRACK START OF THIS BLOCK OF COEFFICIENTS
      IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB) = ICOUNT
C
C     GENERATE ELL0(AB) COEFFICIENTS
      CALL EMAKELL(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,+1,1,2,0)
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
      CALL EMAKELL(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,-1,1,2,0)
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
C        EEEEEEEE  SSSSSS  EEEEEEEE TTTTTTTT  SSSSSS   SSSSSS          C
C        EE       SS    SS EE          TT    SS    SS SS    SS         C
C        EE       SS       EE          TT    SS       SS               C
C        EEEEEE    SSSSSS  EEEEEE      TT     SSSSSS   SSSSSS          C
C        EE             SS EE          TT          SS       SS         C
C        EE       SS    SS EE          TT    SS    SS SS    SS         C
C        EEEEEEEE  SSSSSS  EEEEEEEE    TT     SSSSSS   SSSSSS          C
C                                                                      C
C -------------------------------------------------------------------- C
C  ESETSS CONSTRUCTS ALL ESS0 COEFFICIENTS FOR A SYSTEM WITH CALLS TO  C
C  EMAKESS AND SAVES THEM TO AN EXTERNAL DATA FILE.                    C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTERS A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     INITIALISE COUNT INDEX
      ICOUNT = 0
C
C     LOOP OVER CENTER A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 0
        ENDIF
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNS(1) = NFUNCT(LQN(1)+1,ICNTA)
        DO IBAS=1,NFUNS(1)
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNS(2) = NFUNCT(LQN(2)+1,ICNTB)
        DO JBAS=1,NFUNS(2)
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
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS           C
C     OF (MA,MB). FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE          C
C     ORDERED                                                          C
C     11: = (-|MQN(A)|,-|MQN(B)|)                                      C
C     12: = (-|MQN(A)|,+|MQN(B)|)                                      C
C     21: = (+|MQN(A)|,-|MQN(B)|)                                      C
C     22: = (+|MQN(A)|,+|MQN(B)|)                                      C
C**********************************************************************C
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NFUNS(1)*NFUNS(2)
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMAB  = LQN(1)+LQN(2)+2
      NTUVSS = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
C
C     INDEX TO TRACK START OF THIS BLOCK OF COEFFICIENTS
      IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB) = ICOUNT
C
C     GENERATE ESS0(AB) COEFFICIENTS
      CALL EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,+1,1,2,0)
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
      CALL EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,-1,1,2,0)
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
C        EEEEEEEE  SSSSSS  EEEEEEEE TTTTTTTT LL        SSSSSS          C
C        EE       SS    SS EE          TT    LL       SS    SS         C
C        EE       SS       EE          TT    LL       SS               C
C        EEEEEE    SSSSSS  EEEEEE      TT    LL        SSSSSS          C
C        EE             SS EE          TT    LL             SS         C
C        EE       SS    SS EE          TT    LL       SS    SS         C
C        EEEEEEEE  SSSSSS  EEEEEEEE    TT    LLLLLLLL  SSSSSS          C
C                                                                      C
C -------------------------------------------------------------------- C
C  ESETLS CONSTRUCTS ALL ESS0 COEFFICIENTS FOR A SYSTEM WITH CALLS TO  C
C  EMAKELS AND SAVES THEM TO AN EXTERNAL DATA FILE.                    C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MCT=6,MKP=9,MFL=10000000,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      COMPLEX*16 E11(MB2,MEQ,3),E21(MB2,MEQ,3)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      COMMON/EILS/EILSFL(MFL,24),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTERS A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     INITIALISE COUNT INDEX
      ICOUNT = 0
C
C     LOOP OVER CENTER A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTER B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTER B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTER/MULTI-CENTER OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 0
        ENDIF
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
          LQN(1) = KQN(1)
        ELSE
          LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNS(1) = NFUNCT(LQN(1)+1,ICNTA)
        DO IBAS=1,NFUNS(1)
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNS(2) = NFUNCT(LQN(2)+1,ICNTB)
        DO JBAS=1,NFUNS(2)
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
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS           C
C     OF (MA,MB). FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE          C
C     ORDERED                                                          C
C     11: = (-|MQN(A)|,-|MQN(B)|)                                      C
C     12: = (-|MQN(A)|,+|MQN(B)|)                                      C
C     21: = (+|MQN(A)|,-|MQN(B)|)                                      C
C     22: = (+|MQN(A)|,+|MQN(B)|)                                      C
C**********************************************************************C
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NFUNS(1)*NFUNS(2)
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMAB  = LQN(1)+LQN(2)+1
      NTUVLS = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
C
C     INDEX TO TRACK START OF THIS BLOCK OF COEFFICIENTS
      IADILS(ICNTA,ICNTB,KA,KB,MA,MB) = ICOUNT
C
C     GENERATE ELSI(AB) COEFFICIENTS
      CALL EMAKEB3(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,+1,1,2)
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
      CALL EMAKEB3(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,-1,1,2)
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
      SUBROUTINE EMAKELL(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,I1,I2,IQ)
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
C  DEVELOPMENT OF THE ALGORITHM OF V.R.SAUNDERS.                       C
C -------------------------------------------------------------------- C
C  H.M.QUINEY THE UNIVERSITY OF MELBOURNE (2008).                      C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NFUNS(4)
      DIMENSION EXL(MBS,2),KQNLAB(2),MQNLAB(2),NFLAB(2),COORD(3,2)
C
C     DETERMINE LQNS FROM KQNS
      IF(KQN(I1).GT.0) THEN
        LA = KQN(I1)
      ELSE
        LA =-KQN(I1)-1
      ENDIF
      IF(KQN(I2).GT.0) THEN
        LB = KQN(I2)
      ELSE
        LB =-KQN(I2)-1
      ENDIF
C
C     SUMMATION TERMINATES AT LAM = LA + LB FOR LL PAIRS
      LAMAB = LA + LB
C
C     KAPPA LABELS AND NUMBER OF FUNCTIONS ON A GIVEN CENTER
      KQNLAB(1) = KQN(I1)
      KQNLAB(2) = KQN(I2)
      NFLAB(1)  = NFUNS(I1)
      NFLAB(2)  = NFUNS(I2)
C
C     CARTESIAN COORDINATES OF EACH CENTER
      DO I=1,3
        COORD(I,1) = XYZ(I,I1)
        COORD(I,2) = XYZ(I,I2)
      ENDDO
C
C     BASIS SET EXPONENTS FOR CENTER A
      DO IBAS=1,NFUNS(I1)
       EXL(IBAS,1) = EXPT(IBAS,I1)       
      ENDDO
C
C     BASIS SET EXPONENTS FOR CENTER B
      DO JBAS=1,NFUNS(I2)
       EXL(JBAS,2) = EXPT(JBAS,I2)
      ENDDO                  
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR M PAIRS (-|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQNLAB(1) =-MQN(I1)
      MQNLAB(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLL(E11,EXL,COORD,KQNLAB,MQNLAB,NFLAB,IQ)
C
C     MULTIPLY SOME OF THE ELL0 COEFFICIENTS BY A PHASE TERM
      ITUV=0
      DO MTUV=0,LAMAB
        PHASE = DFLOAT((IALT)**(MTUV))
        DO MDUM=1,((MTUV+1)*(MTUV+2))/2
          ITUV = ITUV+1
          DO M=1,NFUNS(I1)*NFUNS(I2)
            E11(M,ITUV) = E11(M,ITUV)*PHASE
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR M PAIRS (+|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQNLAB(1) = MQN(I1)
      MQNLAB(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLL(E21,EXL,COORD,KQNLAB,MQNLAB,NFLAB,IQ)
C
C     MULTIPLY SOME OF THE ELL0 COEFFICIENTS BY A PHASE TERM
      ITUV=0
      DO MTUV=0,LAMAB
        PHASE = DFLOAT((IALT)**(MTUV))
        DO MDUM=1,((MTUV+1)*(MTUV+2))/2
          ITUV = ITUV+1
          DO M=1,NFUNS(I1)*NFUNS(I2)
            E21(M,ITUV) = E21(M,ITUV)*PHASE
          ENDDO
        ENDDO
      ENDDO
C
C     NOTE THAT E22 AND E12 ARE RELATED TO THESE BY PHASE FACTORS.
C
      RETURN
      END
C
C
      SUBROUTINE EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,I1,I2,IQ)
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
C  DEVELOPMENT OF THE ALGORITHM OF V.R.SAUNDERS.                       C
C -------------------------------------------------------------------- C
C  H.M.QUINEY THE UNIVERSITY OF MELBOURNE (2008).                      C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NFUNS(4)
      DIMENSION EXL(MBS,2),KQNLAB(2),MQNLAB(2),NFLAB(2),COORD(3,2)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
C     DETERMINE LQNS FROM KQNS
      IF(KQN(I1).GT.0) THEN
        LA = KQN(I1)
      ELSE
        LA =-KQN(I1)-1
      ENDIF
      IF(KQN(I2).GT.0) THEN
        LB = KQN(I2)
      ELSE
        LB =-KQN(I2)-1
      ENDIF
C
C     SUMMATION TERMINATES AT LAMAB = LA + LB + 2 FOR SS PAIRS
      LAMAB = LA + LB + 2
C
C     KAPPA LABELS AND NUMBER OF FUNCTIONS ON A GIVEN CENTER
      KQNLAB(1) = KQN(I1)
      KQNLAB(2) = KQN(I2)
      NFLAB(1)  = NFUNS(I1)
      NFLAB(2)  = NFUNS(I2)
C
C     CARTESIAN COORDINATES OF EACH CENTER
      DO I=1,3
        COORD(I,1) = XYZ(I,I1)
        COORD(I,2) = XYZ(I,I2)
      ENDDO
C
C     BASIS SET EXPONENTS FOR CENTER A
      DO IBAS=1,NFUNS(I1)
        EXL(IBAS,1) = EXPT(IBAS,I1)
      ENDDO
C
C     BASIS SET EXPONENTS FOR CENTER B
      DO JBAS=1,NFUNS(I2)
        EXL(JBAS,2) = EXPT(JBAS,I2)
      ENDDO      
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR M PAIRS (-|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQNLAB(1) =-MQN(I1)
      MQNLAB(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQSS(E11,EXL,COORD,KQNLAB,MQNLAB,NFLAB,IQ)
C
C     MULTIPLY SOME OF THE ESS0 COEFFICIENTS BY A PHASE TERM
      ITUV=0
      DO MTUV=0,LAMAB
        PHASE = DFLOAT((IALT)**(MTUV))
        DO MDUM=1,((MTUV+1)*(MTUV+2))/2
          ITUV = ITUV+1
          DO M=1,NFUNS(I1)*NFUNS(I2)
            E11(M,ITUV) = E11(M,ITUV)*PHASE
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR M PAIRS (+|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQNLAB(1) = MQN(I1)
      MQNLAB(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQSS(E21,EXL,COORD,KQNLAB,MQNLAB,NFLAB,IQ)
C
C     MULTIPLY SOME OF THE ESS0 COEFFICIENTS BY A PHASE TERM
      ITUV=0
      DO MTUV=0,LAMAB
        PHASE = DFLOAT((IALT)**(MTUV))
        DO MDUM=1,((MTUV+1)*(MTUV+2))/2
         ITUV = ITUV+1
          DO M=1,NFUNS(I1)*NFUNS(I2)
            E21(M,ITUV) = E21(M,ITUV)*PHASE
          ENDDO
        ENDDO
      ENDDO
C
C     NOTE THAT E22 AND E12 ARE RELATED TO THESE BY SIMPLE FACTORS.
C
      RETURN
      END
C
C
      SUBROUTINE EMAKELS(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,I1,I2,IQ)
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
C  DEVELOPMENT OF THE ALGORITHM OF V.R.SAUNDERS.                       C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NFUNS(4)
      DIMENSION EXL(MBS,2),KQNLAB(2),MQNLAB(2),NFLAB(2),COORD(3,2)
C
C     DETERMINE LQNS FROM KQNS
      IF(KQN(I1).GT.0) THEN
        LA = KQN(I1)
      ELSE
        LA =-KQN(I1)-1
      ENDIF
      IF(KQN(I2).GT.0) THEN
        LB = KQN(I2)
      ELSE
        LB =-KQN(I2)-1
      ENDIF
C
C     THE SUMMATION TERMINATES AT LAMAB = LA + LB + 1 FOR LS PAIRS
      LAMAB = LA + LB + 1
C
C     KAPPA LABELS AND NUMBER OF FUNCTIONS ON A GIVEN CENTER
      KQNLAB(1) = KQN(I1)
      KQNLAB(2) = KQN(I2)
      NFLAB(1)  = NFUNS(I1)
      NFLAB(2)  = NFUNS(I2)
C
C     CARTESIAN COORDINATES OF EACH CENTER
      DO I=1,3
        COORD(I,1) = XYZ(I,I1)
        COORD(I,2) = XYZ(I,I2)
      ENDDO
C
C     BASIS SET EXPONENTS FOR CENTER A
      DO IBAS=1,NFUNS(I1)
       EXL(IBAS,1) = EXPT(IBAS,I1)
      ENDDO
C
C     BASIS SET EXPONENTS FOR CENTER B
      DO JBAS=1,NFUNS(I2)
       EXL(JBAS,2) = EXPT(JBAS,I2)
      ENDDO    
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR M PAIRS (-|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQNLAB(1) =-MQN(I1)
      MQNLAB(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLS(E11,EXL,COORD,KQNLAB,MQNLAB,NFLAB,IQ)
C
C     MULTIPLY SOME OF THE ELSQ COEFFICIENTS BY A PHASE TERM
      ITUV = 0
      DO MTUV=0,LAMAB
        PHASE = DFLOAT((IALT)**(MTUV))
        DO MDUM=1,(MTUV+1)*(MTUV+2)/2
         ITUV = ITUV+1
           DO M=1,NFUNS(I1)*NFUNS(I2)
             E11(M,ITUV) = E11(M,ITUV)*PHASE
           ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR M PAIRS (+|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQNLAB(1) = MQN(I1)
      MQNLAB(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLS(E21,EXL,COORD,KQNLAB,MQNLAB,NFLAB,IQ)
C
C     MULTIPLY SOME OF THE ELS1 COEFFICIENTS BY A PHASE TERM
      ITUV = 0
      DO MTUV=0,LAMAB
        PHASE = DFLOAT((IALT)**(MTUV))
        DO MDUM=1,(MTUV+1)*(MTUV+2)/2
          ITUV = ITUV + 1
          DO M=1,NFUNS(I1)*NFUNS(I2)
            E21(M,ITUV) = E21(M,ITUV)*PHASE
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
      SUBROUTINE EMAKEB3(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,I1,I2)
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
C  EMAKEC3 GENERATES A VECTOR LIST OF E-COEFFICIENTS FOR A BATCH OF    C
C  BREIT INTERACTION INTEGRALS USING A DEVELOPMENT OF THE ALGORITHM OF C
C  V.R.SAUNDERS.                                                       C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      COMPLEX*16 E11(MB2,MEQ,3),E21(MB2,MEQ,3)
C
      DIMENSION EXPT(MBS,4),XYZ(3,4),KQN(4),MQN(4),NFUNS(4)
      DIMENSION EXL(MBS,2),KQNLAB(2),MQNLAB(2),NFLAB(2),COORD(3,2)
C
C     DETERMINE LQNS FROM KQNS
      IF(KQN(I1).GT.0) THEN
        LA = KQN(I1)
      ELSE
        LA =-KQN(I1)-1
      ENDIF
      IF(KQN(I2).GT.0) THEN
        LB = KQN(I2)
      ELSE
        LB =-KQN(I2)-1
      ENDIF
C
C     THE SUMMATION TERMINATES AT LAMAB = LA + LB + 1 FOR LS PAIRS
      LAMAB = LA + LB + 1
C
C     KAPPA LABELS AND NUMBER OF FUNCTIONS ON A GIVEN CENTER
      KQNLAB(1) = KQN(I1)
      KQNLAB(2) = KQN(I2)
      NFLAB(1)  = NFUNS(I1)
      NFLAB(2)  = NFUNS(I2)
C
C     CARTESIAN COORDINATES OF EACH CENTER
      DO I=1,3
        COORD(I,1) = XYZ(I,I1)
        COORD(I,2) = XYZ(I,I2)
      ENDDO
C
C     BASIS SET EXPONENTS FOR CENTER A
      DO IBAS=1,NFUNS(I1)
       EXL(IBAS,1) = EXPT(IBAS,I1)
      ENDDO
C
C     BASIS SET EXPONENTS FOR CENTER B
      DO JBAS=1,NFUNS(I2)
       EXL(JBAS,2) = EXPT(JBAS,I2)
      ENDDO
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR M PAIRS (-|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQNLAB(1) =-MQN(I1)
      MQNLAB(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLS3(E11,EXL,COORD,KQNLAB,MQNLAB,NFLAB)
C
C     MULTIPLY SOME OF THE ELSQ COEFFICIENTS BY A PHASE TERM
      ITUV = 0
      DO MTUV=0,LAMAB
        PHASE = DFLOAT((IALT)**(MTUV))
        DO MDUM=1,(MTUV+1)*(MTUV+2)/2
         ITUV = ITUV+1
           DO M=1,NFUNS(I1)*NFUNS(I2)
             DO IQ=1,3
               E11(M,ITUV,IQ) = PHASE*E11(M,ITUV,IQ)
             ENDDO
           ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR M PAIRS (+|MA|,-|MB|)             C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQNLAB(1) = MQN(I1)
      MQNLAB(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLS3(E21,EXL,COORD,KQNLAB,MQNLAB,NFLAB)
C
C     MULTIPLY SOME OF THE ELS1 COEFFICIENTS BY A PHASE TERM
      ITUV = 0
      DO MTUV=0,LAMAB
        PHASE = DFLOAT((IALT)**(MTUV))
        DO MDUM=1,(MTUV+1)*(MTUV+2)/2
          ITUV = ITUV + 1
          DO M=1,NFUNS(I1)*NFUNS(I2)
            DO IQ=1,3
               E21(M,ITUV,IQ) = PHASE*E21(M,ITUV,IQ)
            ENDDO
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
      SUBROUTINE EQLL(ELL,EXL,COORD,KQN,MQN,NFUN,IQ)
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
C  EQLL EVALUATES THE EQ-COEFFICIENTS FOR LARGE-LARGE CHARGE OVERLAP   C
C  OF G-SPINOR FUNCTIONS FOR ALL PAULI MATRICES IQ = {0,1,2,3}.        C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      COMPLEX*16 CONE
      COMPLEX*16 ELL(MB2,MEQ),ESG(MB2,MEQ)
C
      DIMENSION EXL(MBS,2),RNORM(MBS,2),TEMP(MB2),COORD(3,2),
     &          NFUN(2),KQN(2),JQN(2),LQN(2),MQN(2),LLAB(2),MLAB(2)
C
      COMMON/GSPR/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     THRESHOLD FOR CG-SCREENING PROCESS
      SENS = 1.0D-10
C
C     TOTAL ANGULAR MOMENTUM QUANTUM NUMBER
      JQN(1) = 2*IABS(KQN(1))-1
      JQN(2) = 2*IABS(KQN(2))-1
      RJ10   = DFLOAT(2*JQN(1))
      RJ20   = DFLOAT(2*JQN(2))
      RJ12   = DFLOAT(2*(JQN(1)+2))
      RJ22   = DFLOAT(2*(JQN(2)+2))
C
C     MAP KAPPA QUANTUM NUMBERS ONTO LQN QUANTUM NUMBERS, AND
C     CALCULATE THE APPROPRIATE CLEBSCH GORDAN FACTORS
      IF(KQN(1).LT.0) THEN 
        LQN(1) = -KQN(1)-1
        CAU    =  DSQRT(DFLOAT(JQN(1)+MQN(1))/RJ10)
        CAL    =  DSQRT(DFLOAT(JQN(1)-MQN(1))/RJ10)
      ELSE
        LQN(1) =  KQN(1)
        CAU    = -DSQRT(DFLOAT(JQN(1)+2-MQN(1))/RJ12)
        CAL    =  DSQRT(DFLOAT(JQN(1)+2+MQN(1))/RJ12)
      ENDIF
C      
      IF(KQN(2).LT.0) THEN 
        LQN(2) = -KQN(2)-1
        CBU    =  DSQRT(DFLOAT(JQN(2)+MQN(2))/RJ20)
        CBL    =  DSQRT(DFLOAT(JQN(2)-MQN(2))/RJ20)
      ELSE
        LQN(2) =  KQN(2)
        CBU    = -DSQRT(DFLOAT(JQN(2)+2-MQN(2))/RJ22)
        CBL    =  DSQRT(DFLOAT(JQN(2)+2+MQN(2))/RJ22)
      ENDIF
C
C     DETERMINE THE NUMBER OF FUNCTIONS ON EACH CENTER
      NFUNA = NFUN(1)
      NFUNB = NFUN(2)
      MAXM  = NFUN(1)*NFUN(2)
C
C**********************************************************************C
C     INITIALIZE COMMON GEOMETRICAL INFORMATION                        C
C -------------------------------------------------------------------- C
C     RKAB.EQ.1 -> INCORPORATE RKAB(M) INTO COEFFICIENTS               C
C     RKAB.NE.1 -> SET ALL RKAB(M) = (1.0D0,0.0D0)                     C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTERS A AND B
      AB2 = (COORD(1,1)-COORD(1,2))**2 + (COORD(2,1)-COORD(2,2))**2
     &                                 + (COORD(3,1)-COORD(3,2))**2
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NFUNA
        EXPA = EXL(IBAS,1)
        DO JBAS=1,NFUNB
          EXPB    = EXL(JBAS,2)
          M       = M + 1
          PAB     = EXPA + EXPB
          PX      = (EXPA*COORD(1,1) + EXPB*COORD(1,2))/PAB
          PY      = (EXPA*COORD(2,1) + EXPB*COORD(2,2))/PAB
          PZ      = (EXPA*COORD(3,1) + EXPB*COORD(3,2))/PAB
          P(M)    = PAB
          P2(M)   = PAB*2.0D0
          P22(M)  = P2(M)*P2(M)
          PAX(M)  = PX - COORD(1,1)
          PAY(M)  = PY - COORD(2,1)
          PAZ(M)  = PZ - COORD(3,1)
          PBX(M)  = PX - COORD(1,2)
          PBY(M)  = PY - COORD(2,2)
          PBZ(M)  = PZ - COORD(3,2)
          PA2(M)  = PAX(M)*PAX(M) + PAY(M)*PAY(M) + PAZ(M)*PAZ(M)
          PB2(M)  = PBX(M)*PBX(M) + PBY(M)*PBY(M) + PBZ(M)*PBZ(M)
          RKAB(M) = DEXP(-(EXPA*EXPB*AB2)/PAB)
        ENDDO
      ENDDO
C
C     INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
      LAM  = LQN(1)+LQN(2)
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     INITIALIZE THE COEFFICIENTS TO ZERO
      DO ITUV=1,NTUV
        DO M=1,MAXM
          ELL(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CALCULATION OF CONTRIBUTIONS FROM MQN PAIRS                      C
C**********************************************************************C
C 
C     BASIS PAIR LQNS
      LLAB(1) = LQN(1)
      LLAB(2) = LQN(2)
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR SIGMA_0,Z
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (M-1/2, M'-1/2) CONTRIBUTIONS
C         
C       BASIS PAIR MQNS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)-1)/2
C     
C       DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER PRODUCT)
        CLEBSCH = CAU*CBU
C
C       IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELLQ
        IF(DABS(CLEBSCH).GE.SENS) THEN
C
C         GENERATE THE ES-COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         FIRST CONTRIBUTION TO SIGMA_0
          IF(IQ.EQ.0) THEN
            DO ITUV=1,NTUV
              DO M=1,MAXM
                ELL(M,ITUV) = ELL(M,ITUV) + CLEBSCH*ESG(M,ITUV)
              ENDDO
            ENDDO
          ENDIF
C
C         FIRST CONTRIBUTION TO SIGMA_Z
          IF(IQ.EQ.3) THEN
            DO ITUV=1,NTUV
              DO M=1,MAXM
                ELL(M,ITUV) = ELL(M,ITUV) + CLEBSCH*ESG(M,ITUV)
              ENDDO
            ENDDO
          ENDIF
C
        ENDIF
C
C >>    TERM 22: (M+1/2, M'+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQNS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)+1)/2
C
C       DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER PRODUCT)
        CLEBSCH = CAL*CBL
C
C       IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELLQ
        IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C         GENERATE THE ES-COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         SECOND CONTRIBUTION TO SIGMA_0
          IF(IQ.EQ.0) THEN
            DO ITUV=1,NTUV
              DO M=1,MAXM
                ELL(M,ITUV) = ELL(M,ITUV) + CLEBSCH*ESG(M,ITUV)
              ENDDO
            ENDDO
          ENDIF
C
C         SECOND CONTRIBUTION TO SIGMA_0
          IF(IQ.EQ.3) THEN
            DO ITUV=1,NTUV
              DO M=1,MAXM
                ELL(M,ITUV) = ELL(M,ITUV) - CLEBSCH*ESG(M,ITUV)
              ENDDO
            ENDDO
          ENDIF
C          
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR SIGMA_X,Y
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (M-1/2, M'+1/2) CONTRIBUTIONS
C         
C       BASIS PAIR MQNS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)+1)/2
C     
C       DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER PRODUCT)
        CLEBSCH = CAU*CBL
C
C       IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELLQ
        IF(DABS(CLEBSCH).GE.SENS) THEN
C
C         GENERATE THE ES-COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         FIRST CONTRIBUTION TO SIGMA_X
          IF(IQ.EQ.1) THEN
            DO ITUV=1,NTUV
              DO M=1,MAXM
                ELL(M,ITUV) = ELL(M,ITUV) + CLEBSCH*ESG(M,ITUV)
              ENDDO
            ENDDO
          ENDIF
C
C         FIRST CONTRIBUTION TO SIGMA_Y
          IF(IQ.EQ.2) THEN
            DO ITUV=1,NTUV
              DO M=1,MAXM
                ELL(M,ITUV) = ELL(M,ITUV) - CLEBSCH*ESG(M,ITUV)
              ENDDO
            ENDDO
          ENDIF
C
        ENDIF
C
C >>    TERM 21: (M+1/2, M'-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQNS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER PRODUCT)
        CLEBSCH = CAL*CBU
C
C       IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELLQ
        IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C         GENERATE THE ES-COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         SECOND CONTRIBUTION TO SIGMA_X
          IF(IQ.EQ.1) THEN
            DO ITUV=1,NTUV
              DO M=1,MAXM
                ELL(M,ITUV) = ELL(M,ITUV) + CLEBSCH*ESG(M,ITUV)
              ENDDO
            ENDDO
          ENDIF
C
C         SECOND CONTRIBUTION TO SIGMA_Y
          IF(IQ.EQ.2) THEN
            DO ITUV=1,NTUV
              DO M=1,MAXM
                ELL(M,ITUV) = ELL(M,ITUV) + CLEBSCH*ESG(M,ITUV)
              ENDDO
            ENDDO
          ENDIF
C          
        ENDIF
      ENDIF
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
C     GENERATE LARGE-LARGE NORMALISATION CONSTANTS
      CALL RNLL(RNORM,EXL,LQN,NFUN)
C
C     PRODUCT OF NORMALISATION CONSTANTS
      M=0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M       = M+1
          TEMP(M) = RNORM(IBAS,1)*RNORM(JBAS,2)
        ENDDO
      ENDDO
C
C     MULTIPLY CONSTANTS BY THE ELL0 COEFFICIENTS
C     REMEMBER THAT SIGMA_Y NEEDS AN EXTRA FACTOR OF i
      DO ITUV=1,NTUV
        DO M=1,MAXM
          ELL(M,ITUV) = TEMP(M)*ELL(M,ITUV)
          IF(IQ.EQ.2) THEN
            ELL(M,ITUV) = ELL(M,ITUV)*CONE
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE EQSS(ESS,EXL,COORD,KQN,MQN,NFUN,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                EEEEEEEE  QQQQQQ     SSSSSS   SSSSSS                  C
C                EE       QQ    QQ   SS    SS SS    SS                 C
C                EE      QQ      QQ  SS       SS                       C
C                EEEEEE  QQ      QQ   SSSSSS   SSSSSS                  C
C                EE      QQ      QQ        SS       SS                 C
C                EE       QQ    QQ   SS    SS SS    SS                 C
C                EEEEEEEE  QQQQQQ QQ  SSSSSS   SSSSSS                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQSS EVALUATES THE EQ-COEFFICIENTS FOR SMALL-SMALL CHARGE OVERLAP   C
C  OF G-SPINOR FUNCTIONS FOR ALL PAULI MATRICES IQ = {0,1,2,3}.        C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      COMPLEX*16 CONE
      COMPLEX*16 ESS(MB2,MEQ),ESG(MB2,MEQ),ENSG(MB2,MEQ)
C
      DIMENSION EXL(MBS,2),RNORM(MBS,2),TEMP(MB2),COORD(3,2),
     &          NFUN(2),KQN(2),JQN(2),LQN(2),MQN(2),LLAB(2),MLAB(2)
C
      COMMON/GSPR/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     THRESHOLD FOR CG-SCREENING PROCESS
      SENS = 1.0D-10
C
C     TOTAL ANGULAR MOMENTUM QUANTUM NUMBER
      JQN(1) = 2*IABS(KQN(1))-1
      JQN(2) = 2*IABS(KQN(2))-1
      RJ10   = DFLOAT(2*JQN(1))
      RJ20   = DFLOAT(2*JQN(2))
      RJ12   = DFLOAT(2*(JQN(1)+2))
      RJ22   = DFLOAT(2*(JQN(2)+2))
C
C     MAP KAPPA QUANTUM NUMBERS ONTO LQN QUANTUM NUMBERS
C     CALCULATE THE APPROPRIATE CLEBSCH GORDAN FACTORS
      IF(KQN(1).LT.0) THEN 
        LQN(1) = -KQN(1)-1
        CAU    = -DSQRT(DFLOAT(JQN(1)+2-MQN(1))/RJ12)
        CAL    =  DSQRT(DFLOAT(JQN(1)+2+MQN(1))/RJ12)
      ELSE
        LQN(1) =  KQN(1)
        CAU    =  DSQRT(DFLOAT(JQN(1)+MQN(1))/RJ10)
        CAL    =  DSQRT(DFLOAT(JQN(1)-MQN(1))/RJ10)
      ENDIF
C
      IF(KQN(2).LT.0) THEN 
        LQN(2) = -KQN(2)-1
        CBU    = -DSQRT(DFLOAT(JQN(2)+2-MQN(2))/RJ22)
        CBL    =  DSQRT(DFLOAT(JQN(2)+2+MQN(2))/RJ22)
      ELSE
        LQN(2) =  KQN(2)
        CBU    =  DSQRT(DFLOAT(JQN(2)+MQN(2))/RJ20)
        CBL    =  DSQRT(DFLOAT(JQN(2)-MQN(2))/RJ20)
      ENDIF
C
C     DETERMINE THE NUMBER OF FUNCTIONS ON EACH CENTER
      NFUNA = NFUN(1)
      NFUNB = NFUN(2)
      MAXM  = NFUN(1)*NFUN(2)
C
C**********************************************************************C
C     INITIALIZE COMMON GEOMETRICAL INFORMATION                        C
C -------------------------------------------------------------------- C
C     RKAB.EQ.1 -> INCORPORATE RKAB(M) INTO COEFFICIENTS               C
C     RKAB.NE.1 -> SET ALL RKAB(M) = (1.0D0,0.0D0)                     C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTERS A AND B
      AB2 = (COORD(1,1)-COORD(1,2))**2 + (COORD(2,1)-COORD(2,2))**2
     &                                 + (COORD(3,1)-COORD(3,2))**2
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NFUNA
        EXPA = EXL(IBAS,1)
        DO JBAS=1,NFUNB
          M       = M + 1
          EXPB    = EXL(JBAS,2)
          PAB     = EXPA + EXPB
          PX      = (EXPA*COORD(1,1) + EXPB*COORD(1,2))/PAB
          PY      = (EXPA*COORD(2,1) + EXPB*COORD(2,2))/PAB
          PZ      = (EXPA*COORD(3,1) + EXPB*COORD(3,2))/PAB
          P(M)    = PAB
          P2(M)   = PAB*2.0D0
          P22(M)  = P2(M)*P2(M)
          PAX(M)  = PX - COORD(1,1)
          PAY(M)  = PY - COORD(2,1)
          PAZ(M)  = PZ - COORD(3,1)
          PBX(M)  = PX - COORD(1,2)
          PBY(M)  = PY - COORD(2,2)
          PBZ(M)  = PZ - COORD(3,2)
          PA2(M)  = PAX(M)*PAX(M) + PAY(M)*PAY(M) + PAZ(M)*PAZ(M)
          PB2(M)  = PBX(M)*PBX(M) + PBY(M)*PBY(M) + PBZ(M)*PBZ(M)
          RKAB(M) = DEXP(-(EXPA*EXPB*AB2)/PAB)
        ENDDO
      ENDDO
C
C     MAXIMAL INDEX SUMMATION TERMINAL USED IN THIS ROUTINE
      LAM  = LQN(1)+LQN(2)+2
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     INITIALIZE THE COEFFICIENTS TO ZERO
      DO ITUV=1,NTUV
        DO M=1,MAXM
          ESS(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO      
C
C**********************************************************************C
C     CASE 1: KQN(1).LT.0 AND KQN(2).LT.0                              C
C**********************************************************************C
C
      IF(KQN(1).LT.0.AND.KQN(2).LT.0) THEN
C 
C       BASIS PAIR LQNS
        LLAB(1) = LQN(1)+1
        LLAB(2) = LQN(2)+1
C
C       TERMS 11 AND 22 ARE ONLY NECESSARY FOR SIGMA_0,Z
        IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>      TERM 11: (M-1/2, M'-1/2) CONTRIBUTIONS
C                
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)-1)/2
          MLAB(2) = (MQN(2)-1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER PRODUCT)
          CLEBSCH = CAU*CBU
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6 
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M+1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_0
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Z
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
C
C >>      TERM 22: (M+1/2, M'+1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)+1)/2
          MLAB(2) = (MQN(2)+1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER PRODUCT)
          CLEBSCH = CAL*CBL
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN  
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M+1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_0
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Z
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF

        ENDIF
C
C       TERMS 12 AND 21 ARE ONLY NECESSARY FOR SIGMA_X,Y
        IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>      TERM 12: (M-1/2, M'+1/2) CONTRIBUTIONS
C                
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)-1)/2
          MLAB(2) = (MQN(2)+1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER PRODUCT)
          CLEBSCH = CAU*CBL
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M+1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_X
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Y
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
C
C >>      TERM 21: (M+1/2, M'-1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)+1)/2
          MLAB(2) = (MQN(2)-1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER PRODUCT)
          CLEBSCH = CAL*CBU
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN  
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M+1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_X
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Y
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
        ENDIF
C        
      ENDIF
C
C**********************************************************************C
C     CASE 2: KQN(1).LT.0 AND KQN(2).GT.0                              C
C**********************************************************************C
C
      IF(KQN(1).LT.0.AND.KQN(2).GT.0) THEN 
C
C       BASIS PAIR LQNS
        LLAB(1) = LQN(1)+1
        LLAB(2) = LQN(2)-1
C
C       TERMS 11 AND 22 ARE ONLY NECESSARY FOR SIGMA_0,Z
        IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>      TERM 11: (M-1/2, M'-1/2) CONTRIBUTIONS
C                
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)-1)/2
          MLAB(2) = (MQN(2)-1)/2
C
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER PRODUCT)
          CLEBSCH = CAU*CBU
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(2)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(IBAS,1)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=0]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=0]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INCREASE THE INDEX N BY ONE
            CALL STEPN(ESG,ENSG,LAM,MAXM,2)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M = M+1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)+2
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           FIRST CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=1]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=1]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
C
C >>      TERM 22: (M+1/2, M'+1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)+1)/2
          MLAB(2) = (MQN(2)+1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER PRODUCT)
          CLEBSCH = CAL*CBL
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN  
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(2)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(IBAS,1)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=0]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=0]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INCREASE THE INDEX N BY ONE
            CALL STEPN(ESG,ENSG,LAM,MAXM,2)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)+2
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=1]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=1]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
        ENDIF
C
C       TERMS 12 AND 21 ARE ONLY NECESSARY FOR SIGMA_X,Y
        IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>      TERM 12: (M-1/2, M'+1/2) CONTRIBUTIONS
C                
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)-1)/2
          MLAB(2) = (MQN(2)+1)/2
C
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER PRODUCT)
          CLEBSCH = CAU*CBL
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(2)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(IBAS,1)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_X FOR [N=0,N'=0]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=0]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)+1
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCREASE THE INDEX N BY ONE
            CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M = M+1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_X FOR [N=0,N'=1]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=1]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
C
C >>      TERM 21: (M+1/2, M'-1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)+1)/2
          MLAB(2) = (MQN(2)-1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER PRODUCT)
          CLEBSCH = CAL*CBU
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN  
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(2)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(IBAS,1)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_X FOR [N=0,N'=0]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=0]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)+1
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCREASE THE INDEX N BY ONE
            CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_X FOR [N=0,N'=1]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=1]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
        ENDIF
C
      ENDIF
C
C**********************************************************************C
C     CASE 3: KQN(1).GT.0 AND KQN(2).LT.0                              C
C**********************************************************************C
C
      IF(KQN(1).GT.0.AND.KQN(2).LT.0) THEN
C
C       BASIS PAIR LQNS
        LLAB(1) = LQN(1)-1
        LLAB(2) = LQN(2)+1
C
C       TERMS 11 AND 22 ARE ONLY NECESSARY FOR SIGMA_0,Z
        IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>      TERM 11: (M-1/2, M'-1/2) CONTRIBUTIONS
C                
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)-1)/2
          MLAB(2) = (MQN(2)-1)/2
C
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER PRODUCT)
          CLEBSCH = CAU*CBU
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(1)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=0]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=0]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INCREASE THE INDEX N BY ONE
            CALL STEPN(ESG,ENSG,LAM,MAXM,1)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)+2
            NTUV  = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M = M+1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=1]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=1]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
C
C >>      TERM 22: (M+1/2, M'+1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)+1)/2
          MLAB(2) = (MQN(2)+1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER PRODUCT)
          CLEBSCH = CAL*CBL
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN  
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(1)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=0]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=0]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INCREASE THE INDEX N BY ONE
            CALL STEPN(ESG,ENSG,LAM,MAXM,1)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)+2
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=1]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=1]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
        ENDIF
C
C       TERMS 12 AND 21 ARE ONLY NECESSARY FOR SIGMA_X,Y
        IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>      TERM 12: (M-1/2, M'+1/2) CONTRIBUTIONS
C                
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)-1)/2
          MLAB(2) = (MQN(2)+1)/2
C
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER PRODUCT)
          CLEBSCH = CAU*CBL
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(1)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_X FOR [N=0,N'=0]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=0]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)+1
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCREASE THE INDEX N BY ONE
            CALL STEPN(ESG,ENSG,LAM-1,MAXM,1)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M = M+1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_X FOR [N=0,N'=1]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=1]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
C
C >>      TERM 21: (M+1/2, M'-1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)+1)/2
          MLAB(2) = (MQN(2)-1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER PRODUCT)
          CLEBSCH = CAL*CBU
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN  
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(1)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_X FOR [N=0,N'=0]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=0]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)+1
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCREASE THE INDEX N BY ONE
            CALL STEPN(ESG,ENSG,LAM-1,MAXM,1)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_X FOR [N=0,N'=1]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=1]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
        ENDIF
      ENDIF
C
C**********************************************************************C
C     CASE 4: KQN(1).GT.0 AND KQN(2).GT.0                              C
C**********************************************************************C
C
      IF(KQN(1).GT.0.AND.KQN(2).GT.0) THEN 
C
C       BASIS PAIR LQNS
        LLAB(1) = LQN(1)-1
        LLAB(2) = LQN(2)-1
C
C       TERMS 11 AND 22 ARE ONLY NECESSARY FOR SIGMA_0,Z
        IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>      TERM 11: (M-1/2, M'-1/2) CONTRIBUTIONS
C                
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)-1)/2
          MLAB(2) = (MQN(2)-1)/2
C
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER PRODUCT)
          CLEBSCH = CAU*CBU
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)-2
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = DFLOAT((2*LQN(1)+1)*(2*LQN(2)+1))*CLEBSCH
C
C           FIRST CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=0]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + PREFAC*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=0]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + PREFAC*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INCREASE THE INDEX N BY ONE FOR THE A CENTER
            CALL STEPN(ESG,ENSG,LAM,MAXM,1)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            NTUV  = (LAM+3)*(LAM+4)*(LAM+5)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(2)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(IBAS,1)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_0 FOR [N=1,N'=0]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Z FOR [N=1,N'=0]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INCREASE THE INDEX N BY ONE FOR THE B CENTER
            CALL STEPN(ESG,ENSG,LAM,MAXM,2)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(1)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            NTUV  = (LAM+3)*(LAM+4)*(LAM+5)/6
C
C           FIRST CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=1]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=1]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INCREASE INDEX N IN ENSG ARRAY (CURRENTLY WITH [N=0,N'=1])
C           BY ONE ON THE A CENTER, SO THAT IT OVERWRITES ESG VALUES AND
C           PROVIDES ES-COEFFICIENTS WITH [N=1,N'=1]. (LAM -> LAM + 2)
            CALL STEPN(ENSG,ESG,LAM+2,MAXM,1)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            NTUV  = (LAM+5)*(LAM+6)*(LAM+7)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M = M+1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_0 FOR [N=1,N'=1]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Z FOR [N=1,N'=1]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
C
C >>      TERM 22: (M+1/2, M'+1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)+1)/2
          MLAB(2) = (MQN(2)+1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER PRODUCT)
          CLEBSCH = CAL*CBL
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)-2
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = DFLOAT((2*LQN(1)+1)*(2*LQN(2)+1))*CLEBSCH
C
C           SECOND CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=0]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + PREFAC*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=0]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - PREFAC*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INCREASE THE INDEX N BY ONE FOR THE A CENTER
            CALL STEPN(ESG,ENSG,LAM,MAXM,1)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            NTUV  = (LAM+3)*(LAM+4)*(LAM+5)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(2)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(IBAS,1)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_0 FOR [N=1,N'=0]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Z FOR [N=1,N'=0]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INCREASE THE INDEX N BY ONE FOR THE B CENTER
            CALL STEPN(ESG,ENSG,LAM,MAXM,2)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            NTUV  = (LAM+3)*(LAM+4)*(LAM+5)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(1)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=1]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=1]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INCREASE INDEX N IN ENSG ARRAY (CURRENTLY WITH [N=0,N'=1])
C           BY ONE ON THE A CENTER, SO THAT IT OVERWRITES ESG VALUES AND
C           PROVIDES ES-COEFFICIENTS WITH [N=1,N'=1]. (LAM -> LAM + 2)
            CALL STEPN(ENSG,ESG,LAM+2,MAXM,1)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            NTUV  = (LAM+5)*(LAM+6)*(LAM+7)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M = M+1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_0 FOR [N=1,N'=1]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Z FOR [N=1,N'=1]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
        ENDIF
C
C       TERMS 12 AND 21 ARE ONLY NECESSARY FOR SIGMA_X,Y
        IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>      TERM 12: (M-1/2, M'+1/2) CONTRIBUTIONS
C                
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)-1)/2
          MLAB(2) = (MQN(2)+1)/2
C
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER PRODUCT)
          CLEBSCH = CAU*CBL
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = DFLOAT((2*LQN(1)+1)*(2*LQN(2)+1))*CLEBSCH
C
C           FIRST CONTRIBUTION TO SIGMA_X FOR [N=0,N'=0]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + PREFAC*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=0]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - PREFAC*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)+2
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCREASE THE INDEX N BY ONE FOR THE A CENTER FOR [N=1,N'=0]
            CALL STEPN(ESG,ENSG,LAM-2,MAXM,1)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(2)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(IBAS,1)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_X FOR [N=1,N'=0]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Y FOR [N=1,N'=0]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)+2
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCREASE THE INDEX N BY ONE FOR THE B CENTER FOR [N=0,N'=1]
            CALL STEPN(ESG,ENSG,LAM-2,MAXM,2)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(1)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_X FOR [N=0,N'=1]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=1]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)+4
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCREASE INDEX N IN ENSG ARRAY (CURRENTLY WITH [N=0,N'=1])
C           BY ONE ON THE A CENTER, SO THAT IT OVERWRITES ESG VALUES AND
C           PROVIDES ES-COEFFICIENTS WITH [N=1,N'=1]. (LAM -> LAM + 2)
            CALL STEPN(ENSG,ESG,LAM-4,MAXM,1)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M = M+1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_X FOR [N=1,N'=1]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Y FOR [N=1,N'=1]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) - TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
C
C >>      TERM 21: (M+1/2, M'-1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)+1)/2
          MLAB(2) = (MQN(2)-1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER PRODUCT)
          CLEBSCH = CAL*CBU
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ESS0
          IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = DFLOAT((2*LQN(1)+1)*(2*LQN(2)+1))*CLEBSCH
C
C           SECOND CONTRIBUTION TO SIGMA_X FOR [N=0,N'=0]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + PREFAC*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=0]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + PREFAC*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)+2
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCREASE THE INDEX N BY ONE FOR THE A CENTER
            CALL STEPN(ESG,ENSG,LAM-2,MAXM,1)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(2)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(IBAS,1)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_X FOR [N=1,N'=0]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Y FOR [N=1,N'=0]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)+2
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCREASE THE INDEX N BY ONE FOR THE B CENTER
            CALL STEPN(ESG,ENSG,LAM-2,MAXM,2)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*DFLOAT(2*LQN(1)+1)*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_X FOR [N=0,N'=1]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=1]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)+4
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCREASE INDEX N IN ENSG ARRAY (CURRENTLY WITH [N=0,N'=1])
C           BY ONE ON THE A CENTER, SO THAT IT OVERWRITES ESG VALUES AND
C           PROVIDES ES-COEFFICIENTS WITH [N=1,N'=1]. (LAM -> LAM + 2)
            CALL STEPN(ENSG,ESG,LAM-4,MAXM,1)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = 4.0D0*CLEBSCH
            M = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M = M+1
                TEMP(M) = PREFAC*EXL(IBAS,1)*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_X FOR [N=1,N'=1]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Y FOR [N=1,N'=1]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ESS(M,ITUV) = ESS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
        ENDIF
C
      ENDIF
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
C     GENERATE SMALL-SMALL NORMALISATION CONSTANTS
      CALL RNSS(RNORM,EXL,LQN,NFUN)
C
C     MAXIMAL INDEX SUMMATION TERMINAL USED IN THIS ROUTINE
      LAM  = LQN(1)+LQN(2)+4
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     PRODUCT OF NORMALISATION CONSTANTS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M       = M + 1
          TEMP(M) = RNORM(IBAS,1)*RNORM(JBAS,2)
        ENDDO
      ENDDO
C
C     BRING THESE COEFFICIENTS INTO THE ESSQ VALUES AND ALSO FACTOR i
C     REMEMBER THAT SIGMA_Y NEEDS AN EXTRA FACTOR OF i
      DO ITUV=1,NTUV
        DO M=1,MAXM
          ESS(M,ITUV) = TEMP(M)*ESS(M,ITUV)
          IF(IQ.EQ.2) THEN
            ESS(M,ITUV) = ESS(M,ITUV)*CONE
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE EQLS(ELS,EXL,COORD,KQN,MQN,NFUN,IQ)
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
      DIMENSION EXL(MBS,2),RNORM(MBS,2),TEMP(MB2),COORD(3,2),
     &          NFUN(2),KQN(2),JQN(2),LQN(2),MQN(2),LLAB(2),MLAB(2)
C
      COMMON/GSPR/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     SET CLEBSCH-GORDAN SENSITIVITY
      SENS = 1.0D-10
C
C     TOTAL ANGULAR MOMENTUM QUANTUM NUMBER
      JQN(1) = 2*IABS(KQN(1))-1
      JQN(2) = 2*IABS(KQN(2))-1
      RJ10   = DFLOAT(2*JQN(1))
      RJ20   = DFLOAT(2*JQN(2))
      RJ12   = DFLOAT(2*(JQN(1)+2))
      RJ22   = DFLOAT(2*(JQN(2)+2))
C
C     MAP KAPPA QUANTUM NUMBERS ONTO LQN QUANTUM NUMBERS
C     CALCULATE THE APPROPRIATE CLEBSCH GORDAN FACTORS
      IF(KQN(1).LT.0) THEN 
        LQN(1) = -KQN(1)-1
        CAU    =  DSQRT(DFLOAT(JQN(1)+MQN(1))/RJ10)
        CAL    =  DSQRT(DFLOAT(JQN(1)-MQN(1))/RJ10)
      ELSE
        LQN(1) =  KQN(1)
        CAU    = -DSQRT(DFLOAT(JQN(1)+2-MQN(1))/RJ12)
        CAL    =  DSQRT(DFLOAT(JQN(1)+2+MQN(1))/RJ12)
      ENDIF
C
      IF(KQN(2).LT.0) THEN 
        LQN(2) = -KQN(2)-1
        CBU    = -DSQRT(DFLOAT(JQN(2)+2-MQN(2))/RJ22)
        CBL    =  DSQRT(DFLOAT(JQN(2)+2+MQN(2))/RJ22)
      ELSE
        LQN(2) = KQN(2)
        CBU    = DSQRT(DFLOAT(JQN(2)+MQN(2))/RJ20)
        CBL    = DSQRT(DFLOAT(JQN(2)-MQN(2))/RJ20)
      ENDIF
C
C     DETERMINE THE NUMBER OF FUNCTIONS ON EACH CENTER
      NFUNA = NFUN(1)
      NFUNB = NFUN(2)
      MAXM  = NFUN(1)*NFUN(2)
C
C**********************************************************************C
C     INITIALIZE COMMON GEOMETRICAL INFORMATION                        C
C -------------------------------------------------------------------- C
C     RKAB.EQ.1 -> INCORPORATE RKAB(M) INTO COEFFICIENTS               C
C     RKAB.NE.1 -> SET ALL RKAB(M) = (1.0D0,0.0D0)                     C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTERS A AND B
      AB2 = (COORD(1,1)-COORD(1,2))**2 + (COORD(2,1)-COORD(2,2))**2
     &                                 + (COORD(3,1)-COORD(3,2))**2
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NFUNA
        EXPA = EXL(IBAS,1)
        DO JBAS=1,NFUNB
          M       = M + 1
          EXPB    = EXL(JBAS,2)
          PAB     = EXPA + EXPB
          PX      = (EXPA*COORD(1,1) + EXPB*COORD(1,2))/PAB
          PY      = (EXPA*COORD(2,1) + EXPB*COORD(2,2))/PAB
          PZ      = (EXPA*COORD(3,1) + EXPB*COORD(3,2))/PAB
          P(M)    = PAB
          P2(M)   = PAB*2.0D0
          P22(M)  = P2(M)*P2(M)
          PAX(M)  = PX - COORD(1,1)
          PAY(M)  = PY - COORD(2,1)
          PAZ(M)  = PZ - COORD(3,1)
          PBX(M)  = PX - COORD(1,2)
          PBY(M)  = PY - COORD(2,2)
          PBZ(M)  = PZ - COORD(3,2)
          PA2(M)  = PAX(M)*PAX(M) + PAY(M)*PAY(M) + PAZ(M)*PAZ(M)
          PB2(M)  = PBX(M)*PBX(M) + PBY(M)*PBY(M) + PBZ(M)*PBZ(M)
          RKAB(M) = DEXP(-(EXPA*EXPB*AB2)/PAB)
        ENDDO
      ENDDO
C
C     MAXIMAL INDEX SUMMATION TERMINAL USED IN THIS ROUTINE
C     WHY IS THIS +3 AND NOT +1?
      LAM  = LQN(1)+LQN(2)+3
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     INITIALIZE THE COEFFICIENTS TO ZERO
      DO ITUV=1,NTUV
        DO M=1,MAXM
          ELS(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CASE 1: KQN(2).LT.0                                              C
C**********************************************************************C
C
      IF(KQN(2).LT.0) THEN
C 
C       BASIS PAIR LQNS
        LLAB(1) = LQN(1)
        LLAB(2) = LQN(2)+1
C
C       TERMS 11 AND 22 ARE ONLY NECESSARY FOR SIGMA_0,Z
        IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>      TERM 11: (M-1/2, M'-1/2) CONTRIBUTIONS
C         
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)-1)/2
          MLAB(2) = (MQN(2)-1)/2
C
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER PRODUCT)
          CLEBSCH = CAU*CBU
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
          IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)+1
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*CLEBSCH
            M      = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_0
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Z
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
C
C >>      TERM 22: (M+1/2, M'+1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)+1)/2
          MLAB(2) = (MQN(2)+1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER PRODUCT)
          CLEBSCH = CAL*CBL
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
          IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LLAB(1)+LLAB(2)
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*CLEBSCH
            M      = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_0
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Z
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) - TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
        ENDIF
C
C       TERMS 12 AND 21 ARE ONLY NECESSARY FOR SIGMA_X,Y
        IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>      TERM 12: (M-1/2, M'+1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)-1)/2
          MLAB(2) = (MQN(2)+1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER-LOWER PRODUCT)
          CLEBSCH = CAU*CBL
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
          IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)+1
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*CLEBSCH
            M      = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_X
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Y
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) - TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
C
C >>      TERM 21: (M+1/2, M'-1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)+1)/2
          MLAB(2) = (MQN(2)-1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER-UPPER PRODUCT)
          CLEBSCH = CAL*CBU
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
          IF(DABS(CLEBSCH).GE.SENS) THEN
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)+1
            NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*CLEBSCH
            M      = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_X
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Y
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + TEMP(M)*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
        ENDIF
C
C**********************************************************************C
C     CASE 2: KQN(2).GT.0                                              C
C**********************************************************************C
C
      ELSE
C 
C       BASIS PAIR LQNS
        LLAB(1) = LQN(1)
        LLAB(2) = LQN(2)-1
C
C       TERMS 11 AND 22 ARE ONLY NECESSARY FOR SIGMA_0,Z
        IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>      TERM 11: (M-1/2, M'-1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)-1)/2
          MLAB(2) = (MQN(2)-1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER-UPPER PRODUCT)
          CLEBSCH = CAU*CBU
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
          IF(DABS(CLEBSCH).GE.SENS) THEN
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)
            NTUV = LAM*(LAM+1)*(LAM+2)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            FACTOR = DFLOAT(2*LQN(2)+1)*CLEBSCH
C
C           FIRST CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=0]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + FACTOR*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=0]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + FACTOR*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INCREASE THE INDEX N BY ONE
            CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)
            NTUV = (LAM+2)*(LAM+3)*(LAM+4)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*CLEBSCH
            M      = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=1]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=1]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
C
C >>      TERM 22: (M+1/2, M'+1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)+1)/2
          MLAB(2) = (MQN(2)+1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER-LOWER PRODUCT)
          CLEBSCH = CAL*CBL
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
          IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)
            NTUV = LAM*(LAM+1)*(LAM+2)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            FACTOR = DFLOAT(2*LQN(2)+1)*CLEBSCH
C
C           SECOND CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=0]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + FACTOR*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=0]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) - FACTOR*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INCREASE THE INDEX N BY ONE
            CALL STEPN(ESG,ENSG,LAM,MAXM,2)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)
            NTUV = (LAM+2)*(LAM+3)*(LAM+4)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*CLEBSCH
            M      = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_0 FOR [N=0,N'=1]
            IF(IQ.EQ.0) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=1]
            IF(IQ.EQ.3) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) - TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
        ENDIF
C
C       TERMS 12 AND 21 ARE ONLY NECESSARY FOR SIGMA_X,Y
        IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>      TERM 12: (M-1/2, M'+1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)-1)/2
          MLAB(2) = (MQN(2)+1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER-LOWER PRODUCT)
          CLEBSCH = CAU*CBL
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
          IF(DABS(CLEBSCH).GE.SENS) THEN
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)
            NTUV = LAM*(LAM+1)*(LAM+2)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            FACTOR = DFLOAT(2*LQN(2)+1)*CLEBSCH
C
C           FIRST CONTRIBUTION TO SIGMA_X FOR [N=0,N'=0]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + FACTOR*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=0]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) - FACTOR*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INCREASE THE INDEX N BY ONE
            CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)
            NTUV = (LAM+2)*(LAM+3)*(LAM+4)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*CLEBSCH
            M      = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M + 1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           FIRST CONTRIBUTION TO SIGMA_X FOR [N=0,N'=1]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           FIRST CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=1]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) - TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
C
C >>      TERM 21: (M+1/2, M'-1/2) CONTRIBUTIONS
C
C         BASIS PAIR MQNS
          MLAB(1) = (MQN(1)+1)/2
          MLAB(2) = (MQN(2)-1)/2
C      
C         DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER-UPPER PRODUCT)
          CLEBSCH = CAL*CBU
C
C         IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
          IF(DABS(CLEBSCH).GE.SENS) THEN
C
C           GENERATE THE ES-COEFFICIENTS
            CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)
            NTUV = LAM*(LAM+1)*(LAM+2)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            FACTOR = DFLOAT(2*LQN(2)+1)*CLEBSCH
C
C           SECOND CONTRIBUTION TO SIGMA_X FOR [N=0,N'=0]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + FACTOR*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=0]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + FACTOR*ESG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           INCREASE THE INDEX N BY ONE
            CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C           INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
            LAM  = LQN(1)+LQN(2)
            NTUV = (LAM+2)*(LAM+3)*(LAM+4)/6
C
C           INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
            PREFAC = -2.0D0*CLEBSCH
            M      = 0
            DO IBAS=1,NFUNA
              DO JBAS=1,NFUNB
                M       = M+1
                TEMP(M) = PREFAC*EXL(JBAS,2)
              ENDDO
            ENDDO
C
C           SECOND CONTRIBUTION TO SIGMA_X FOR [N=0,N'=1]
            IF(IQ.EQ.1) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
C           SECOND CONTRIBUTION TO SIGMA_Y FOR [N=0,N'=1]
            IF(IQ.EQ.2) THEN
              DO ITUV=1,NTUV
                DO M=1,MAXM
                  ELS(M,ITUV) = ELS(M,ITUV) + TEMP(M)*ENSG(M,ITUV)
                ENDDO
              ENDDO
            ENDIF
C
          ENDIF
        ENDIF
      ENDIF
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
      CALL RNLS(RNORM,EXL,LQN,NFUN)
C
C     MAXIMAL INDEX SUMMATION TERMINAL USED IN THIS ROUTINE
      LAM  = LQN(1)+LQN(2)+3
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE THE LS NORMALISATION FACTORS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          TEMP(M) = RNORM(IBAS,1)*RNORM(JBAS,2)
        ENDDO
      ENDDO
C
C     BRING THESE COEFFICIENTS INTO THE ELSQ VALUES AND ALSO FACTOR i
C     REMEMBER THAT SIGMA_Y NEEDS AN EXTRA FACTOR OF i
      DO ITUV=1,NTUV
        DO M=1,MAXM
          ELS(M,ITUV) = TEMP(M)*ELS(M,ITUV)*CONE
          IF(IQ.EQ.2) THEN
            ELS(M,ITUV) = ELS(M,ITUV)*CONE
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE EQLS3(ELS,EXL,COORD,KQN,MQN,NFUN)
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
      DIMENSION EXL(MBS,2),RNORM(MBS,2),TEMP(MB2),COORD(3,2),
     &          NFUN(2),KQN(2),JQN(2),LQN(2),MQN(2),LLAB(2),MLAB(2)
C
      COMMON/GSPR/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
      DATA SENS/1.0D-10/
C
C     IMAGINARY UNIT
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     TOTAL ANGULAR MOMENTUM QUANTUM NUMBER
      JQN(1) = 2*IABS(KQN(1))-1
      JQN(2) = 2*IABS(KQN(2))-1
      RJ10   = DFLOAT(2*JQN(1))
      RJ20   = DFLOAT(2*JQN(2))
      RJ12   = DFLOAT(2*(JQN(1)+2))
      RJ22   = DFLOAT(2*(JQN(2)+2))
C
C     MAP KAPPA QUANTUM NUMBERS ONTO LQN QUANTUM NUMBERS
C     CALCULATE THE APPROPRIATE CLEBSCH GORDAN FACTORS
      IF(KQN(1).LT.0) THEN 
        LQN(1) =-KQN(1)-1
        CAU    = DSQRT(DFLOAT(JQN(1)+MQN(1))/RJ10)
        CAL    = DSQRT(DFLOAT(JQN(1)-MQN(1))/RJ10)
      ELSE
        LQN(1) = KQN(1)
        CAU    =-DSQRT(DFLOAT(JQN(1)+2-MQN(1))/RJ12)
        CAL    = DSQRT(DFLOAT(JQN(1)+2+MQN(1))/RJ12)
      ENDIF
C
      IF(KQN(2).LT.0) THEN 
        LQN(2) =-KQN(2)-1
        CBU    =-DSQRT(DFLOAT(JQN(2)+2-MQN(2))/RJ22)
        CBL    = DSQRT(DFLOAT(JQN(2)+2+MQN(2))/RJ22)
      ELSE
        LQN(2) = KQN(2)
        CBU    = DSQRT(DFLOAT(JQN(2)+MQN(2))/RJ20)
        CBL    = DSQRT(DFLOAT(JQN(2)-MQN(2))/RJ20)
      ENDIF
C
C     DETERMINE THE NUMBER OF FUNCTIONS ON EACH CENTER
      NFUNA = NFUN(1)
      NFUNB = NFUN(2)
      MAXM  = NFUN(1)*NFUN(2)
C
C**********************************************************************C
C     INITIALIZE COMMON GEOMETRICAL INFORMATION                        C
C -------------------------------------------------------------------- C
C     RKAB.EQ.1 -> INCORPORATE RKAB(M) INTO COEFFICIENTS               C
C     RKAB.NE.1 -> SET ALL RKAB(M) = (1.0D0,0.0D0)                     C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTERS A AND B
      AB2 = (COORD(1,1)-COORD(1,2))**2 + (COORD(2,1)-COORD(2,2))**2
     &                                 + (COORD(3,1)-COORD(3,2))**2
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NFUNA
        EXPA = EXL(IBAS,1)
        DO JBAS=1,NFUNB
          M       = M + 1
          EXPB    = EXL(JBAS,2)
          PAB     = EXPA + EXPB
          PX      = (EXPA*COORD(1,1) + EXPB*COORD(1,2))/PAB
          PY      = (EXPA*COORD(2,1) + EXPB*COORD(2,2))/PAB
          PZ      = (EXPA*COORD(3,1) + EXPB*COORD(3,2))/PAB
          P(M)    = PAB
          P2(M)   = PAB*2.0D0
          P22(M)  = P2(M)*P2(M)
          PAX(M)  = PX - COORD(1,1)
          PAY(M)  = PY - COORD(2,1)
          PAZ(M)  = PZ - COORD(3,1)
          PBX(M)  = PX - COORD(1,2)
          PBY(M)  = PY - COORD(2,2)
          PBZ(M)  = PZ - COORD(3,2)
          PA2(M)  = PAX(M)*PAX(M) + PAY(M)*PAY(M) + PAZ(M)*PAZ(M)
          PB2(M)  = PBX(M)*PBX(M) + PBY(M)*PBY(M) + PBZ(M)*PBZ(M)
          RKAB(M) = DEXP(-(EXPA*EXPB*AB2)/PAB)
        ENDDO
      ENDDO
C
C     MAXIMAL INDEX SUMMATION TERMINAL USED IN THIS ROUTINE
      LAM  = LQN(1)+LQN(2)+1
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     INITIALIZE THE COEFFICIENTS TO ZERO
      DO IQ=1,3
        DO ITUV=1,NTUV
          DO M=1,MAXM
            ELS(M,ITUV,IQ) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CASE 1: KQN(2).LT.0                                              C
C**********************************************************************C
C
      IF(KQN(2).LT.0) THEN
C 
C       BASIS PAIR LQNS
        LLAB(1) = LQN(1)
        LLAB(2) = LQN(2)+1
C
C >>    TERM 11: (M-1/2, M'-1/2) CONTRIBUTIONS
C         
C       BASIS PAIR MQNS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)-1)/2
C
C       DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER PRODUCT)
        CLEBSCH = CAU*CBU
C
C       IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
        IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C         GENERATE THE ES-COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM  = LQN(1)+LQN(2)+1
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
          M = 0
          PREFAC =-2.0D0*CLEBSCH
          DO IBAS=1,NFUNA
            DO JBAS=1,NFUNB
              M = M+1
              TEMP(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         FIRST CONTRIBUTION TO SIGMA_Z
          DO ITUV=1,NTUV
            DO M=1,MAXM
              ELS(M,ITUV,3) = ELS(M,ITUV,3) + TEMP(M)*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (M+1/2, M'+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQNS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)+1)/2
C      
C       DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER PRODUCT)
        CLEBSCH = CAL*CBL
C
C       IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
        IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM  = LLAB(1)+LLAB(2)
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         GENERATE THE ES-COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS   
          M = 0
          PREFAC =-2.0D0*CLEBSCH
          DO IBAS=1,NFUNA
            DO JBAS=1,NFUNB
              M = M+1
              TEMP(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         SECOND CONTRIBUTION TO SIGMA_Z
          DO ITUV=1,NTUV
            DO M=1,MAXM
              ELS(M,ITUV,3) = ELS(M,ITUV,3) - TEMP(M)*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C       TERMS 12 AND 21 ARE ONLY NECESSARY FOR SIGMA_X AND SIGMA_Y
C
C >>    TERM 12: (M-1/2, M'+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQNS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)+1)/2
C      
C       DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER-LOWER PRODUCT)
        CLEBSCH = CAU*CBL
C
C       IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
        IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C         GENERATE THE ES-COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM  = LQN(1)+LQN(2)+1
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
          M = 0
          PREFAC =-2.0D0*CLEBSCH
          DO IBAS=1,NFUNA
            DO JBAS=1,NFUNB
              M = M+1
              TEMP(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         FIRST CONTRIBUTION TO SIGMA_X AND SIGMA_Y
          DO ITUV=1,NTUV
            DO M=1,MAXM
              ELS(M,ITUV,1) = ELS(M,ITUV,1) + TEMP(M)*ESG(M,ITUV)
              ELS(M,ITUV,2) = ELS(M,ITUV,2) - TEMP(M)*ESG(M,ITUV)
            ENDDO
          ENDDO          
C
        ENDIF
C
C >>    TERM 21: (M+1/2, M'-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQNS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)-1)/2
C      
C       DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER-UPPER PRODUCT)
        CLEBSCH = CAL*CBU
C
C       IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
        IF(DABS(CLEBSCH).GE.SENS) THEN
C
C         GENERATE THE ES-COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM  = LQN(1)+LQN(2)+1
          NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C         INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
          M = 0
          PREFAC =-2.0D0*CLEBSCH
          DO IBAS=1,NFUNA
            DO JBAS=1,NFUNB
              M = M+1
              TEMP(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         SECOND CONTRIBUTION TO SIGMA_X AND SIGMA_Y
          DO ITUV=1,NTUV
            DO M=1,MAXM
              ELS(M,ITUV,1) = ELS(M,ITUV,1) + TEMP(M)*ESG(M,ITUV)
              ELS(M,ITUV,2) = ELS(M,ITUV,2) + TEMP(M)*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C**********************************************************************C
C     CASE 2: KQN(2).GT.0                                              C
C**********************************************************************C
C
      ELSE
C 
C       BASIS PAIR LQNS
        LLAB(1) = LQN(1)
        LLAB(2) = LQN(2)-1
C
C       TERMS 11 AND 22 ARE ONLY NECESSARY FOR SIGMA_0,Z
C
C >>    TERM 11: (M-1/2, M'-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQNS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)-1)/2
C      
C       DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER-UPPER PRODUCT)
        CLEBSCH = CAU*CBU
C
C       IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
        IF(DABS(CLEBSCH).GE.SENS) THEN
C
C         GENERATE THE ES-COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM  = LQN(1)+LQN(2)
          NTUV = LAM*(LAM+1)*(LAM+2)/6
C
C         INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
          FACTOR = DFLOAT(2*LQN(2)+1)*CLEBSCH
C
C         FIRST CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=0]
          DO ITUV=1,NTUV
            DO M=1,MAXM
              ELS(M,ITUV,3) = ELS(M,ITUV,3) + FACTOR*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM  = LQN(1)+LQN(2)
          NTUV = (LAM+2)*(LAM+3)*(LAM+4)/6
C
C         INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
          M = 0
          PREFAC =-2.0D0*CLEBSCH
          DO IBAS=1,NFUNA
            DO JBAS=1,NFUNB
              M = M+1
              TEMP(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         FIRST CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=1]
          DO ITUV=1,NTUV
            DO M=1,MAXM
              ELS(M,ITUV,3) = ELS(M,ITUV,3) + TEMP(M)*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (M+1/2, M'+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQNS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)+1)/2
C      
C       DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER-LOWER PRODUCT)
        CLEBSCH = CAL*CBL
C
C       IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
        IF(DABS(CLEBSCH).GE.SENS) THEN 
C
C         GENERATE THE ES-COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM  = LQN(1)+LQN(2)
          NTUV = LAM*(LAM+1)*(LAM+2)/6
C
C         INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
          FACTOR = DFLOAT(2*LQN(2)+1)*CLEBSCH
C
C         SECOND CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=0]
          DO ITUV=1,NTUV
            DO M=1,MAXM
              ELS(M,ITUV,3) = ELS(M,ITUV,3) - FACTOR*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM,MAXM,2)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM  = LQN(1)+LQN(2)
          NTUV = (LAM+2)*(LAM+3)*(LAM+4)/6
C
C         INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
          M = 0
          PREFAC =-2.0D0*CLEBSCH
          DO IBAS=1,NFUNA
            DO JBAS=1,NFUNB
              M = M+1
              TEMP(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         SECOND CONTRIBUTION TO SIGMA_Z FOR [N=0,N'=1]
          DO ITUV=1,NTUV
            DO M=1,MAXM
              ELS(M,ITUV,3) = ELS(M,ITUV,3) - TEMP(M)*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C       TERMS 12 AND 21 ARE ONLY NECESSARY FOR SIGMA_X,Y
C
C >>    TERM 12: (M-1/2, M'+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQNS
        MLAB(1) = (MQN(1)-1)/2
        MLAB(2) = (MQN(2)+1)/2
C    
C       DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (UPPER-LOWER PRODUCT)
        CLEBSCH = CAU*CBL
C
C       IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
        IF(DABS(CLEBSCH).GE.SENS) THEN
C
C         GENERATE THE ES-COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM  = LQN(1)+LQN(2)
          NTUV = LAM*(LAM+1)*(LAM+2)/6
C
C         INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
          FACTOR = DFLOAT(2*LQN(2)+1)*CLEBSCH
C
C         FIRST CONTRIBUTION TO SIGMA_X AND SIGMA_Y FOR [N=0,N'=0]
          DO ITUV=1,NTUV
            DO M=1,MAXM
              ELS(M,ITUV,1) = ELS(M,ITUV,1) + FACTOR*ESG(M,ITUV)
              ELS(M,ITUV,2) = ELS(M,ITUV,2) - FACTOR*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM  = LQN(1)+LQN(2)
          NTUV = (LAM+2)*(LAM+3)*(LAM+4)/6
C
C         INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
          M = 0
          PREFAC =-2.0D0*CLEBSCH
          DO IBAS=1,NFUNA
            DO JBAS=1,NFUNB
              M = M+1
              TEMP(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         FIRST CONTRIBUTION TO SIGMA_X AND SIGMA_Y FOR [N=0,N'=1]
          DO ITUV=1,NTUV
            DO M=1,MAXM
              ELS(M,ITUV,1) = ELS(M,ITUV,1) + TEMP(M)*ENSG(M,ITUV)
              ELS(M,ITUV,2) = ELS(M,ITUV,2) - TEMP(M)*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (M+1/2, M'-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQNS
        MLAB(1) = (MQN(1)+1)/2
        MLAB(2) = (MQN(2)-1)/2
C      
C       DEFINE PRODUCT OF CLEBSCH-GORDON FACTORS (LOWER-UPPER PRODUCT)
        CLEBSCH = CAL*CBU
C
C       IF THE CG COEFFICIENT IS LARGE ENOUGH, CALCULATE ELSQ
        IF(DABS(CLEBSCH).GE.SENS) THEN
C
C         GENERATE THE ES-COEFFICIENTS
          CALL ESGTF(ESG,LLAB,MLAB,MAXM)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM  = LQN(1)+LQN(2)
          NTUV = LAM*(LAM+1)*(LAM+2)/6
C
C         INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
          FACTOR = DFLOAT(2*LQN(2)+1)*CLEBSCH
C
C         SECOND CONTRIBUTION TO SIGMA_X AND SIGMA_Y FOR [N=0,N'=0]
          DO ITUV=1,NTUV
            DO M=1,MAXM
              ELS(M,ITUV,1) = ELS(M,ITUV,1) + FACTOR*ESG(M,ITUV)
              ELS(M,ITUV,2) = ELS(M,ITUV,2) + FACTOR*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE
          CALL STEPN(ESG,ENSG,LAM-1,MAXM,2)
C
C         INDEX SUMMATION TERMINAL BASED ON MAX DEGREE OF HGTF
          LAM  = LQN(1)+LQN(2)
          NTUV = (LAM+2)*(LAM+3)*(LAM+4)/6
C
C         INCORPORATE KINETIC BALANCE CONTRACTION COEFFICIENTS
          M = 0
          PREFAC =-2.0D0*CLEBSCH
          DO IBAS=1,NFUNA
            DO JBAS=1,NFUNB
              M = M+1
              TEMP(M) = PREFAC*EXL(JBAS,2)
            ENDDO
          ENDDO
C
C         SECOND CONTRIBUTION TO SIGMA_X AND SIGMA_Y FOR [N=0,N'=1]
          DO ITUV=1,NTUV
            DO M=1,MAXM
              ELS(M,ITUV,1) = ELS(M,ITUV,1) + TEMP(M)*ENSG(M,ITUV)
              ELS(M,ITUV,2) = ELS(M,ITUV,2) + TEMP(M)*ENSG(M,ITUV)
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
      CALL RNLS(RNORM,EXL,LQN,NFUN)
C
C     MAXIMAL INDEX SUMMATION TERMINAL USED IN THIS ROUTINE
      LAM  = LQN(1)+LQN(2)+3
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE THE LS NORMALISATION FACTORS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          TEMP(M) = RNORM(IBAS,1)*RNORM(JBAS,2)
        ENDDO
      ENDDO
C
C     BRING THESE COEFFICIENTS INTO THE ELSQ VALUES AND ALSO FACTOR i
C     REMEMBER THAT SIGMA_Y NEEDS AN EXTRA FACTOR OF i
      DO IQ=1,3
        DO ITUV=1,NTUV
          DO M=1,MAXM
            ELS(M,ITUV,IQ) = CONE*TEMP(M)*ELS(M,ITUV,IQ)
            IF(IQ.EQ.2) THEN
              ELS(M,ITUV,IQ) = CONE*ELS(M,ITUV,IQ)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
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
C  THE REQUIRED COEFFICIENTS ARE GENERATED BY A CALL TO VRS, WHICH IS  C
C  CONSTRUCTED ACCORDING TO THE RECURRENCE RELATIONS DEFINED BY        C
C  V.R.SAUNDERS. THE OUTPUT OF VRS IS THEN ADJUSTED TO INCLUDE THE     C
C  ANGULAR NORMALISATION CONSTANTS, AS WELL AS A PHASE FACTOR TO       C
C  CONVERT FROM THE SCHIFF TO THE CONDON-SHORTLEY PHASE CONVENTION.    C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C     LQN(I) - TARGET LQN VALUES ON CENTERS A AND B.                   C
C     MQN(I) - TARGET MQN VALUES ON CENTERS A AND B.                   C
C -------------------------------------------------------------------- C
C  H.M.QUINEY THE UNIVERSITY OF MELBOURNE (2008).                      C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6)
C
      COMPLEX*16 ESG(MB2,MEQ)
      DIMENSION FACT(MKP),LQN(2),MQN(2),MQNLAB(2)
C
      COMMON/GSPR/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
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
      MQNLAB(1) = -MQN(1)
      MQNLAB(2) =  MQN(2)
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
C  V.R.SAUNDERS,"MOLECULAR INTEGRALS FOR GAUSSIAN-TYPE FUNCTIONS",     C
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
C     LQN(I) - TARGET LQN VALUES ON CENTERS A AND B.                   C
C     MQN(I) - TARGET MQN VALUES ON CENTERS A AND B.                   C
C -------------------------------------------------------------------- C
C  H.M.QUINEY THE UNIVERSITY OF MELBOURNE (2008).                      C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,ITUVRS=716,
     &                          ML2=MKP+1,MEQ=(ML2+1)*(ML2+2)*(ML2+3)/6,
     &                          ML4=ML2*2,MRC=(ML4+1)*(ML4+2)*(ML4+3)/6)
C
      COMPLEX*16 ESG(MB2,MEQ),ETEMP(MB2*MRC)
C
      DIMENSION NFUN(2),LQN(2),MQN(2)
C
      COMMON/GSPR/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     IMPORT L AND BASIS PAIR MQNS FOR LOCAL USE
      LQNA = LQN(1)
      LQNB = LQN(2)
      MQNA = MQN(1)
      MQNB = MQN(2)
      LMAX = LQNA + LQNB
C
C     AT THE MOMENT, THESE ROUTINES HANDLE AT MOST G-TYPE FUNCTIONS.
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
C      SIMULTANEOUS STEP OF LQN AND MQN ON CENTER A
       ISTART = 0
       LAM    = 0
       IF(IABS(MQNA).NE.0) THEN 
         CALL STEPLM(ETEMP,LAM,ISTART,MQNA,MAXM,1)
       ENDIF
C
C      STEP 2:
C      GENERATE E[LQNA,MQNA;0,0] FROM E[|MQNA|,MQNA;0,0]
C      USING THE STEP OF LQN ONLY ON CENTER A
       IF(LQNA.GT.IABS(MQNA)) THEN
         CALL STEPL(ETEMP,LAM,ISTART,LQNA,MQNA,MAXM,1)
       ENDIF
C
C      STEP 3:
C      GENERATE E[LQNA,MQNA;|MQNB|,MQNB] FROM E[LQNA,MQNA;0,0] 
C      USING SIMULTANEOUS STEP OF LQN AND MQN ON CENTER B
       IF(IABS(MQNB).GT.0) THEN 
         CALL STEPLM(ETEMP,LAM,ISTART,MQNB,MAXM,2)
       ENDIF
C
C      STEP 4:
C      GENERATE E[LQNA,MQNA;LQNB,MQNB] FROM E[LQNA,MQNA;|MQNB|,MQNB]
C      USING THE STEP OF LQN ONLY ON CENTER B 
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
C    LQN   - L-QUANTUM NUMBER OF THE CENTER TO BE INCREMENTED.         C
C    ICNT  - CENTER TO STEP UP.                                        C
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
      COMMON/GSPR/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     IMPORT GEOMETRICAL VALUES FOR CENTER OF INTEREST
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
C     MAIN LOOP: FOR EACH M-QUANTUM NUMBER ON THIS CENTER
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
      COMMON/GSPR/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     IMPORT GEOMETRICAL VALUES FOR CENTER OF INTEREST
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
        RFACT1 = -DFLOAT(LQN1+IABS(MQN))/DBLE(LQN1-IABS(MQN)+1)
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
      DIMENSION PP(MB2),PX(MB2),PY(MB2),PZ(MB2)
C
      COMMON/ACSS/INABCD(0:ML4,0:ML4,0:ML4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/GSPR/P(MB2),P2(MB2),P22(MB2),RKAB(MB2),PA2(MB2),PB2(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2)
C
C     IMPORT GEOMETRICAL VALUES FOR CENTER OF INTEREST
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
      SUBROUTINE RNLL(RNORM,EXPT,LQN,NFUN)
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
      DIMENSION RNORM(MBS,2),EXPT(MBS,2),LQN(2),NFUN(2)
      DATA PI,TWOLOG/3.1415926535897932D0,6.93147180559945309D-1/
C
      DO ICNT=1,2
        T1  = DSQRT(PI)
        F1  = 0.5D0
        GML = DLOG(T1)
        DO M=2,LQN(ICNT)+2
          GML = GML+DLOG(F1)
          F1  = F1 + 1.0D0
        ENDDO
        RLA = DFLOAT(LQN(ICNT))
        GA1 = TWOLOG - GML
        RA1 = RLA + 1.50D0
        DO M=1,NFUN(ICNT)
          ELOG          = DLOG(2.0D0*EXPT(M,ICNT))
          RNORM(M,ICNT) = DEXP(0.5D0*(GA1 + RA1*ELOG))
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE RNSS(RNORM,EXPT,LQN,NFUN)
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
      DIMENSION RNORM(MBS,2),EXPT(MBS,2),LQN(2),NFUN(2)
      DATA PI,TWOLOG/3.1415926535897932D0,6.93147180559945309D-1/
C
      DO ICNT=1,2
        T1  = DSQRT(PI)
        F1  = 0.5D0
        GML = DLOG(T1)
        DO M=2,LQN(ICNT)+3
          GML = GML+DLOG(F1)
          F1  = F1+1.0D0
        ENDDO
        RLA = DFLOAT(LQN(ICNT))
        GA1 = TWOLOG - GML
        RA1 = RLA + 0.50D0
        DO M=1,NFUN(ICNT)
          ELOG           = DLOG(2.0D0*EXPT(M,ICNT))
          RNORM(M,ICNT) = DEXP(0.5D0*(GA1 + RA1*ELOG))
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE RNLS(RNORMLS,EXPT,LQN,NFUN)
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
      DIMENSION RNORMLS(MBS,2),EXPT(MBS,2),LQN(2),NFUN(2)
      DATA PI,TWOLOG/3.1415926535897932D0,6.93147180559945309D-1/
C
      T1L  = DSQRT(PI)
      F1L  = 0.5D0
      GMLL = DLOG(T1L)
C
      DO M=2,LQN(1)+2
        GMLL = GMLL+DLOG(F1L)
        F1L  = F1L + 1.0D0
      ENDDO
C
      RLAL = DFLOAT(LQN(1))
      GA1L = TWOLOG - GMLL
      RA1L = RLAL + 1.5D0
C
      DO M=1,NFUN(1)
        ELOGL        = DLOG(2.0D0*EXPT(M,1))
        RNORMLS(M,1) = DEXP(0.5D0*(GA1L + RA1L*ELOGL))
      ENDDO
C
      T1S  = DSQRT(PI)
      F1S  = 0.5D0
      GMLS = DLOG(T1S)
C
      DO M=2,LQN(2)+3
        GMLS = GMLS + DLOG(F1S)
        F1S  = F1S + 1.0D0
      ENDDO
C
      RLAS = DFLOAT(LQN(2))
      GA1S = TWOLOG - GMLS
      RA1S = RLAS + 0.5D0
C
      DO M=1,NFUN(2)
        ELOGS        = DLOG(2.0D0*EXPT(M,2))
        RNORMLS(M,2) = DEXP(0.5D0*(GA1S+RA1S*ELOGS))
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
C            DDDDDDD  NN    NN  OOOOOO  RRRRRRR  MM     MM             C
C            DD    DD NNN   NN OO    OO RR    RR MMM   MMM             C
C            DD    DD NNNN  NN OO    OO RR    RR MMMM MMMM             C
C            DD    DD NN NN NN OO    OO RR    RR MM MMM MM             C
C            DD    DD NN  NNNN OO    OO RRRRRRR  MM  M  MM             C
C            DD    DD NN   NNN OO    OO RR    RR MM     MM             C
C            DDDDDDD  NN    NN  OOOOOO  RR    RR MM     MM             C
C                                                                      C
C -------------------------------------------------------------------- C
C  DNORM CALCULATES A SCALE NORM FOR A REAL OR COMPLEX PART OF A LIST  C
C  ECFF OF LENGTH NMAX, AND STORES THE RESULT IN SCL.                  C
C**********************************************************************C
C
      COMPLEX*16 ECFF(NMAX)
C
      DIMENSION ECMP(NMAX)
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
C
