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
C      SSSSSS  WW        WW IIII RRRRRRR  LL      EEEEEEEE SSSSSS      C
C     SS    SS WW        WW  II  RR    RR LL      EE      SS    SS     C
C     SS       WW   WW   WW  II  RR    RR LL      EE      SS           C
C      SSSSSS  WW  WWWW  WW  II  RR    RR LL      EEEEEE   SSSSSS      C
C           SS WW WW  WW WW  II  RRRRRRR  LL      EE            SS     C
C     SS    SS WWWW    WWWW  II  RR    RR LL      EE      SS    SS     C
C      SSSSSS  WW        WW IIII RR    RR LLLLLLL EEEEEEEE SSSSSS      C
C                                                                      C
C ******************************************************************** C
C                                                                      C
C  *   *   *   *  MM       MM BBBBBBB  PPPPPPP TTTTTTTT *   *   *   *  C
C    *   *   *    MMM     MMM BB    BB PP    PP   TT      *   *   *    C
C  *   *   *   *  MMMM   MMMM BB    BB PP    PP   TT    *   *   *   *  C
C    *   *   *    MM MM MM MM BBBBBBB  PP    PP   TT      *   *   *    C
C  *   *   *   *  MM  MMM  MM BB    BB PPPPPPP    TT    *   *   *   *  C
C    *   *   *    MM   M   MM BB    BB PP         TT      *   *   *    C
C  *   *   *   *  MM       MM BBBBBBB  PP         TT    *   *   *   *  C
C                                                                      C
C ******************************************************************** C
C                                                                      C
C        A RELATIVISTIC MOLECULAR ELECTRONIC STRUCTURE PROGRAM         C
C            BASED ON THE ANALYTIC FINITE BASIS SET METHOD.            C
C                                                                      C
C                 H.M.QUINEY, H.SKAANE (OXFORD, 1997)                  C
C                       D. FLYNN (UNIMELB, 2020)                       C
C                                                                      C
C -------------------------------------------------------------------- C
C                          HAMILTONIANS (HMLT)                         C
C                          -------------------                         C
C ▶ 'NORL' NON-RELATIVISTIC HAMILTONIAN (PAULI EQUATION).              C
C ▶ 'BARE' BARE NUCLEUS DIRAC HAMILTONIAN (NO ELECTRON INTERACTION).   C
C ▶ 'DHFR' DIRAC-COULOMB HAMILTONIAN.                                  C
C ▶ 'DHFP' DIRAC-COULOMB HAMILTONIAN (+1ST ORDER BREIT AND QED).       C
C ▶ 'DHFB' DIRAC-COULOMB-BREIT HAMILTONIAN.                            C
C ▶ 'DHFQ' DIRAC-COULOMB-BREIT HAMILTONIAN WITH LEADING-ORDER QED.     C
C -------------------------------------------------------------------- C
C                   CALCULATION TREE OPTIONS (TREE)                    C
C                   -------------------------------                    C
C  HFSCF: HARTREE-FOCK SCF CALCULATION.                                C
C  MBPTN: NTH-ORDER MANY-BODY DIAGRAMMATIC PERTURBATION THEORY.        C
C  MCSCF: MULTI-CONFIGURATIONAL SCF CALCULATION.                       C
C  EXPVL: CALCULATION OF MOLECULAR EXPECTATION VALUES.                 C
C  PLOTS: VISUALS (ELECTROMAGNETIC FIELDS, POTENTIALS, FORM FACTORS).  C
C -------------------------------------------------------------------- C
C                          TABLE OF CONTENTS                           C
C                          -----------------                           C
C   [1] INPUT/OUTPUT: READ FROM INPUT FILE AND SUMMARISE DATA.         C
C   [2] LABELS: MOLECULAR GEOMETRY AND FOCK MATRIX BLOCK ADDRESSES.    C
C   [3] DENSITIES: MOLECULAR DENSITIES, ENERGIES AND LEVEL SHIFTING.   C
C   [4] ATOMIC HARTREE-FOCK: AVERAGE OF CONFIG. ATOMIC SCF.            C
C   [5] MOLECULAR HARTREE-FOCK: MANY-CENTRE HARTREE-FOCK SCF.          C
C   [8] MBPTN: CORRELATION ENERGY CALCULATION ROUTINES.                C
C  [11] R-INTS: BOYS INTEGRALS R-INTEGRALS AND RELATED QUANTITIES.     C
C  [12] EQ-COEFFS: BASIS FUNCTION OVERLAP SPIN-STRUCTURE FACTORS.      C
C  [13] GQ-COEFFS: ANALYTIC DERIVS OF OVERLAP SPIN-STRUCTURE FACTORS.  C
C  [14] SCREENING: ROUTINES TO ESTIMATE FOCK MATRIX CONTRIBUTIONS.     C
C                                                                      C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/TCPU/TBEG,TNUC,TATM,TSCF,TMPT,TPRP,TPLT,TTOT
C
C     READ DATA FROM USER-SPECIFIED INPUT FILE
      CALL CARDIN
C
C     THIS IS EXPLICITLY AN MBPT PROGRAM
C     IF(TREE.NE.'MBPTN') RETURN
C
C     OPEN FILE FOR TERMINAL RECORD
      OPEN(UNIT=7,FILE=TRIM(OUTFL)//'_mbpt.out',STATUS='UNKNOWN')
C
C     PRINT SUMMARY OUTPUTS IF READING IN
      IF(READIN) THEN
C
C       PRINT SUMMARY OF INPUT DATA
        CALL INPUT
C
C       PRINT MEMORY ALLOCATION SUMMARY
        CALL MEMORY
C
      ENDIF
C
C     INTER-ATOMIC ANGLES AND NUCLEAR REPULSION ENERGY
      CALL NUCGEOM
C
C     NUCLEAR POTENTIALS
      CALL CPU_TIME(T0)
      CALL NUCCOUL
      CALL NUCVCPL
      CALL CPU_TIME(TNUC)
      TNUC = TNUC-T0
C
C     FOCK MATRIX SYMMETRY TYPE INDICES
      CALL FOCKIND
C
C     CARTESIAN EXPANSION INDICES FOR BASIS FUNCTION OVERLAP PAIRS
      CALL CARTIND
C
C     MANY-BODY DIAGRAMMATIC PERTURBATION THEORY
      CALL CPU_TIME(T0)
      CALL MBPTN
      CALL CPU_TIME(TMPT)
      TMPT = TMPT-T0
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
C   [C] MEMORY: WRITE A SUMMARY OF MEMORY REQUIREMENTS OF BIG ARRAYS.  C
C   [D] OUTPUT: WRITE A SUMMARY OF TOTAL CALCULATION STATS/DATA.       C
C   [E] PHYSPRM: PHYSICAL CONSTANTS, TRANSCENDENTALS AND ELEMENT NAMES.C
C   [F] GAMGEN: LIST OF GAMMA FUNCTIONS FOR INT AND HALF-INT ARGS.     C
C   [G] FACTRL: LIST OF FACTORIALS AND DOUBLE FACTORIALS.              C
C**********************************************************************C
C
C
      SUBROUTINE CARDIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            CCCCCC     AA    RRRRRRR  DDDDDDD IIII NN    NN           C
C           CC    CC   AAAA   RR    RR DD    DD II  NNN   NN           C
C           CC        AA  AA  RR    RR DD    DD II  NNNN  NN           C
C           CC       AA    AA RR    RR DD    DD II  NN NN NN           C
C           CC       AAAAAAAA RRRRRRR  DD    DD II  NN  NNNN           C
C           CC    CC AA    AA RR    RR DD    DD II  NN   NNN           C
C            CCCCCC  AA    AA RR    RR DDDDDDD IIII NN    NN           C
C                                                                      C
C                          INPUT ROUTINE FOR                           C
C                          ** B E R T H A **                           C
C -------------------------------------------------------------------- C
C  CARDIN READS AND PREPARES DATA FROM A USER-SPECIFIED INPUT FILE.    C
C  THIS IS ALSO WHERE ATOMIC ELEMENT NAMES AND CV ARE SPECIFIED.       C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C     INCLUDE 'omp_lib.h'
C
      LOGICAL FILETHERE
C
      CHARACTER*1  DUMLIN
      CHARACTER*5  NMDL
      CHARACTER*6  CNFG
      CHARACTER*7  HMINT(50),PTYPE(10)
      CHARACTER*9  BTYP
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 COEF(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/FILL/NCNF(MCT,0:MEL,MKP+1),NLVL(MCT,0:MEL),CNFG(MCT)
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/LSHF/SHLEV(4),SHLV,ILEV
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/PLOT/NPTYPE,PTYPE
      COMMON/PT1B/NHMINT,HMINT
      COMMON/SHLL/ACFF,BCFF,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
C
C     GENERATE SOME MATHEMATICAL LISTS
      CALL FACTRL
      CALL GAMGEN
C
C     MANUAL CHOICES
      EQFILE = .TRUE.
      OPENMP = .FALSE.
C
C**********************************************************************C
C     MOLECULE NAME AND CALCULATION TYPE                               C
C**********************************************************************C
C
C     MOLECULE OUTPUT STRING
      READ(5, *) DUMLIN
      READ(5, *) MOLCL
C
C     CALCULATION TREE: HFSCF, MBPTN, MCSCF, EXPVL, PLOTS
      READ(5, *) DUMLIN
      READ(5, *) TREE
      TREE = 'MBPTN'
C
C     CALCULATION TREE CHECK -- IF UNKNOWN, EXIT.
      IF(TREE.NE.'HFSCF'.AND.TREE.NE.'MBPTN'.AND.TREE.NE.'MCSCF'.AND.
     &   TREE.NE.'EXPVL'.AND.TREE.NE.'PLOTS') THEN
        WRITE(6, *) 'In CARDIN: invalid calculation tree. ',TREE
        WRITE(7, *) 'In CARDIN: invalid calculation tree. ',TREE
        STOP
      ENDIF
C
C     HAMILTONIAN: NORL, BARE, DHFR, DHFP, DHFB OR DHFQ
10    FORMAT(A4)
      READ(5, *) DUMLIN
      READ(5,10) HMLT
C
C     HAMILTONIAN CHECK -- IF UNKNOWN, EXIT.
      IF(HMLT.NE.'NORL'.AND.HMLT.NE.'BARE'.AND.HMLT.NE.'DHFR'.AND.
     &   HMLT.NE.'DHFP'.AND.HMLT.NE.'DHFB'.AND.HMLT.NE.'DHFQ') THEN
        WRITE(6, *) 'In CARDIN: unknown HMLT value. HMLT = ',HMLT
        WRITE(7, *) 'In CARDIN: unknown HMLT value. HMLT = ',HMLT
        STOP
      ENDIF
C
C     WAVE FUNCTION FILE NAME
      WFNFL = 'output/'//TRIM(MOLCL)//'_'//HMLT//'.wfn'
C
C     OUTPUT FILE NAME
      OUTFL = 'output/'//TRIM(MOLCL)//'_'//HMLT//'_'//TREE
C
C     CONTINUING CALCULATION: READ-IN (TRUE), NEW START (FALSE)
      READ(5, *) DUMLIN
      READ(5, *) READIN
C
C**********************************************************************C
C     ATOMIC CENTRES AND BASIS FUNCTIONS                               C
C**********************************************************************C
C
C     NUMBER OF ATOMIC CENTRES
      READ(5, *) DUMLIN
      READ(5, *) NCNT
C
C     CHECK THAT NCNT CAN BE SUPPORTED BY SYSTEM PARAMETERS
      IF(NCNT.GT.MCT) THEN
        WRITE(6, *) 'In CARDIN: number of centres exceeds MCT.',MCT
        WRITE(7, *) 'In CARDIN: number of centres exceeds MCT.',MCT
        STOP
      ENDIF
C
C     INITIALISE MAXIMUM LQN AND DIMENSION COUNTERS
      LBIG = 0
      NDIM = 0
C
C     LOOP OVER ATOMIC CENTRES
      DO IZ=1,NCNT
C
C       CARTESIAN COORDINATES OF THIS CENTRE
        READ(5, *) DUMLIN
        READ(5, *) DUMLIN
        READ(5, *) (BXYZ(J,IZ),J=1,3)
C
C       NUCLEUS: MODEL, CHARGE, MASS
        READ(5, *) DUMLIN
        READ(5, *) NMDL(IZ),ZNUC(IZ),ANUC(IZ)
C
C       CHECK THAT NUCLEAR MODEL IS ALLOWED
        IF(NMDL(IZ).NE.'POINT'.AND.NMDL(IZ).NE.'GAUSS'.AND.
     &     NMDL(IZ).NE.'FERMI'.AND.NMDL(IZ).NE.'UNIFM') THEN
          WRITE(6, *) 'In CARDIN: illegal nuclear model.',NMDL(IZ)
          WRITE(7, *) 'In CARDIN: illegal nuclear model.',NMDL(IZ)
          STOP
        ENDIF
C
C       ATOMIC BASIS: TYPE AND LMAX
        READ(5, *) DUMLIN
        READ(5, *) BTYP,LMAX
C
C       CHECK THAT BASIS TYPE IS ALLOWED
        IF(BTYP.NE.'OPTIMISED'.AND.BTYP.NE.'GEOMETRIC'
     &                        .AND.BTYP.NE.'EVENTEMPR') THEN
          WRITE(6, *) 'In CARDIN: illegal basis type.',BTYP
          WRITE(7, *) 'In CARDIN: illegal basis type.',BTYP
          STOP
        ENDIF
C
C       NUMBER OF KAPPA VALUES FOR THIS ATOM
        NKAP(IZ) = 2*LMAX+1
C
C       CHECK THAT LMAX CAN BE SUPPORTED BY SYSTEM PARAMETERS
        IF(2*LMAX+1.GT.MKP) THEN
          WRITE(6, *) 'In CARDIN: LMAX runs outside MKP storage.'
          WRITE(7, *) 'In CARDIN: LMAX runs outside MKP storage.'
          STOP
        ENDIF
C
C       ELECTRON CONFIGURATION CHOICE AND NUCLEAR ELECTRONIC CHARGE
        READ(5, *) DUMLIN
        READ(5, *) CNFG(IZ),IQNC(IZ)
C
C       CHECK THAT CNFG IS ALLOWED
        IF(CNFG(IZ).NE.'AUFBAU'.AND.CNFG(IZ).NE.'MANUAL') THEN
          WRITE(6, *) 'In CARDIN: illegal CNFG option.',CNFG(IZ)
          WRITE(7, *) 'In CARDIN: illegal CNFG option.',CNFG(IZ)
          STOP
        ENDIF
C
C       IF FILLING IS MANUAL, IMPORT ATOMIC ELECTRON CONFIGURATION
        IF(CNFG(IZ).EQ.'MANUAL') THEN
          READ(5, *) DUMLIN
          DO L=0,LMAX
            READ(5, *) NLVL(IZ,L),(NCNF(IZ,L,N),N=1,NLVL(IZ,L))
          ENDDO
        ENDIF
C
C       UPDATE OVERALL MAXIMUM OCCURRING LQN
        IF(LMAX.GT.LBIG) LBIG = LMAX
C
C       OPTIMISED EXPONENTS FROM A RECORDED LIST
        IF(BTYP.EQ.'OPTIMISED') THEN
C
C         TITLE FILLER
          READ(5, *) DUMLIN
C
C         READ IN THE OPTIMISED ORBITAL EXPONENTS FOR EACH LQN
          DO LQN=0,LMAX
C
C           READ NUMBER OF BASIS FUNCTIONS
            READ(5, *) NFNC(LQN,IZ)
C
C           CHECK THAT THIS CAN BE SUPPORTED BY SYSTEM PARAMETERS
            IF(NFNC(LQN,IZ).GT.MBS) THEN
              NOOPS = NFNC(LQN,IZ)
              WRITE(6, *) 'In CARDIN: too many basis functions.',NOOPS
              WRITE(7, *) 'In CARDIN: too many basis functions.',NOOPS
              STOP
            ENDIF
C
C           READ BASIS EXPONENTS FROM A LIST
            DO IBAS=1,NFNC(LQN,IZ)
              READ(5, *) BEXL(IBAS,LQN,IZ)
            ENDDO
C
          ENDDO
C
C       GEOMETRIC BASIS FUNCTIONS
        ELSEIF(BTYP.EQ.'GEOMETRIC') THEN
C
C         TITLE FILLER
          READ(5, *) DUMLIN
C
C         GENERATE THE EVEN TEMPERED ORBITAL EXPONENTS FOR EACH LQN
          DO LQN=0,LMAX
C
C           READ GENERATING PARAMETERS A, B AND NFNC
            READ(5, *) ALPH,BETA,NFNC(LQN,IZ)
C
C           GENERATE NFNC BASIS EXPONENTS USING VARIABLE ZETA
            ZETA = ALPH
            DO IBAS=1,NFNC(LQN,IZ)
              BEXL(IBAS,LQN,IZ) = ZETA
              ZETA = ZETA*BETA
            ENDDO
C
C           CHECK THAT NFNC CAN BE SUPPORTED BY SYSTEM PARAMETERS
            IF(NFNC(LQN,IZ).GT.MBS) THEN
              NOOPS = NFNC(LQN,IZ)
              WRITE(6, *) 'In CARDIN: too many basis functions.',NOOPS
              WRITE(7, *) 'In CARDIN: too many basis functions.',NOOPS
              STOP
            ENDIF
C
          ENDDO
C
C       EVEN-TEMPERED GEOMETRIC BASIS FUNCTIONS
        ELSEIF(BTYP.EQ.'EVENTEMPR') THEN
C
C         TITLE FILLER
          READ(5, *) DUMLIN
C
C         GENERATE THE EVEN TEMPERED ORBITAL EXPONENTS FOR EACH LQN
          DO LQN=0,LMAX
C
C           READ GENERATING PARAMETERS A, B, NFNC AND STARTING POINT
            READ(5, *) ALPH,BETA,NFNC(LQN,IZ),N1
C
C           GENERATE NFNC BASIS EXPONENTS USING VARIABLE ZETA
            ZETA = ALPH
            DO IBAS=1,N1-1
              ZETA = ZETA*BETA
            ENDDO
            DO IBAS=1,NFNC(LQN,IZ)
              BEXL(IBAS,LQN,IZ) = ZETA
              ZETA = ZETA*BETA
            ENDDO
C
C           CHECK THAT NFNC CAN BE SUPPORTED BY SYSTEM PARAMETERS
            IF(NFNC(LQN,IZ).GT.MBS) THEN
              NOOPS = NFNC(LQN,IZ)
              WRITE(6, *) 'In CARDIN: too many basis functions.',NOOPS
              WRITE(7, *) 'In CARDIN: too many basis functions.',NOOPS
              STOP
            ENDIF
C
          ENDDO
C
C       END IF STATEMENT FOR TYPE OF BASIS FUNCTION
        ENDIF
C
C       LOOP OVER ALL LQNS IN THIS CENTRE AND ADD TO FOCK DIMENSION
        DO LQN=0,LMAX
C
C         EXTEND DIMENSION OF FOCK MATRIX
          NDIM = NDIM + 4*(2*LQN+1)*NFNC(LQN,IZ)
C
C         ASSIGN KQN VALUES
          IF(LQN.NE.0) THEN
            KAPA(2*LQN  ,IZ) = LQN
          ENDIF
          KAPA(2*LQN+1,IZ) =-LQN-1
C
        ENDDO
C
C     END LOOP OVER ATOMIC CENTRES
      ENDDO
C
C     TOTAL DIMENSION DEPENDING ON CHOICE OF HAMILTONIAN
      IF(HMLT.EQ.'NORL') THEN
        NDIM = NDIM/2
        NSKP = 0
      ELSE
        NSKP = NDIM/2
      ENDIF
C
C     CHECK THAT SYSTEM PARAMETERS CAN SUPPORT NDIM
      IF(NDIM.GT.MDM) THEN
        WRITE(6, *) 'In CARDIN: matrix dimension exceeds maximum.',NDIM
        WRITE(7, *) 'In CARDIN: matrix dimension exceeds maximum.',NDIM
        STOP
      ENDIF
C
C**********************************************************************C
C     CLOSED/OPEN SHELL DETAILS                                        C
C**********************************************************************C
C
C     NUMBER OF CLOSED- AND OPEN-SHELL ORBITALS AND OPEN SHELL ELECTRONS
      READ(5, *) DUMLIN
      READ(5, *) NCLS,NOPN,NOELEC
C
C     TOTAL NUMBER OF 'OCCUPIED' AND 'VIRTUAL' ORBITALS IN SYSTEM
      NOCC = NCLS+NOPN
      NVRT = NDIM-NSKP-NOCC
C
C     INVALID CHOICE OF NOPN
      IF(NOPN.LT.0) THEN
        WRITE(6, *) 'In CARDIN: invalid value NOPN.',NOPN
        WRITE(7, *) 'In CARDIN: invalid value NOPN.',NOPN
        STOP
      ENDIF
C
C     CLOSED-SHELL MOLECULAR CONFIGURATION
      IF(NOPN.EQ.0) THEN
C
C       LABEL THE CLOSED-SHELL ELECTRON ORBITALS
        DO JCL=1,NCLS
          ICLS(JCL) = JCL
        ENDDO
C
C     OPEN-SHELL MOLECULAR CONFIGURATION
      ELSEIF(NOPN.GT.0) THEN
C
C       FRACTIONAL OCCUPANCY OF THE OPEN-SHELL ORBITALS
        FOPN = DFLOAT(NOELEC)/DFLOAT(NOPN)
C
C       SYMMETRY-ADAPTED CONFIGURATION OF OPEN-SHELL ORBITALS
        READ(5, *) DUMLIN
        READ(5, *) ACFF,BCFF,(IOPN(M),M=1,NOPN)
C
C       LABEL THE CLOSED-SHELL ELECTRON ORBITALS
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
      ENDIF
C
C**********************************************************************C
C     LEVEL SHIFTING AND INTEGRAL INCLUSION STARTING POINT             C
C**********************************************************************C
C
C     LEVEL SHIFT PARAMETER FOR EACH INTEGRAL STAGE (SKAANE 4.4.3)
      READ(5, *) DUMLIN
      READ(5, *) (SHLEV(N),N=1,3)
C
C     STARTING STAGE OF INTEGRAL INCLUSION LEVEL (1-3)
      READ(5, *) DUMLIN
      READ(5, *) ILEV
C
C     CHECK FOR INVALID CHOICE OF ILEV
      IF(ILEV.LT.1.AND.ILEV.GT.4) THEN
        WRITE(6, *) 'In CARDIN: invalid starting stage. ILEV = ',ILEV
        WRITE(7, *) 'In CARDIN: invalid starting stage. ILEV = ',ILEV
        STOP
      ENDIF
C
C     REASONS TO CHANGE THE INTEGRAL INCLUSION LEVEL (STAGE)
      IF(HMLT.EQ.'NORL') THEN
        ILEV = 1
      ENDIF
C
      IF(NCNT.EQ.1.OR.READIN) THEN
        ILEV = 3
      ENDIF
C
C     IMPLEMENT THE STARTING SHIFT FACTOR
      IF(ILEV.GE.1.AND.ILEV.LE.3) THEN
        SHLV = SHLEV(ILEV)
      ELSEIF(ILEV.EQ.4) THEN
        SHLV = 0.0D0
      ENDIF
C
C**********************************************************************C
C     NON-HF CALCULATION DETAILS                                       C
C**********************************************************************C
C
C     EXPECTATION VALUE CALCULATIONS: ORTHGNL,MAGDIPL ETC
      IF(TREE.EQ.'EXPVL') THEN
C
C       READ NUMBER OF EXPECTATION VALUES
        READ(5, *) DUMLIN
        READ(5, *) NHMINT
C
C       PLACE LIMIT ON POSSIBLE NUMBER OF THEM
        IF(NHMINT.LT.1.OR.NHMINT.GT.50) THEN
          WRITE(6, *) 'In CARDIN: invalid number of expectation values.'
          WRITE(7, *) 'In CARDIN: invalid number of expectation values.'
          STOP
        ENDIF
C
C       READ IN EACH INTERACTION HAMILTONIAN
        DO N=1,NHMINT
          READ(5, *) HMINT(N)
        ENDDO
C
      ENDIF
C
C     DATA PLOTTING: AMPLITUDES, EM FIELDS AND POTENTIALS
      IF(TREE.EQ.'PLOTS') THEN
C
C       READ NUMBER OF PLOT TYPES
        READ(5, *) DUMLIN
        READ(5, *) NPTYPE
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
          READ(5, *) PTYPE(N)
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C     READ IN ANY EXTERNAL DATA FILES                                  C
C**********************************************************************C
C
C     READ IN A WAVE FUNCTION FILE IF PROMPTED
C      IF(READIN) THEN
        INQUIRE(FILE=WFNFL,EXIST=FILETHERE)
        IF(FILETHERE) THEN
          OPEN(UNIT=8,FILE=WFNFL,STATUS='UNKNOWN')
          REWIND(UNIT=8)
          DO I=1,NDIM
            READ(8, *) EIGN(I),(COEF(J,I),J=1,NDIM)
          ENDDO
          CLOSE(UNIT=8)
        ELSE
          READIN = .FALSE.
        ENDIF
C      ENDIF
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
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*20 STAMP
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/TCPU/TCPUS(8)
      COMMON/TMMD/TMMDS(9)
      COMMON/TSCF/TSCFS(37)
C
C     INITIALISE TIME COUNTERS
      DO NT=1,8
        TCPUS(NT) = 0.0D0
      ENDDO
C
      DO NT=1,37
        TSCFS(NT) = 0.0D0
      ENDDO
C
      DO NT=1,9
        TMMDS(NT) = 0.0D0
      ENDDO
C
C     PRINT A BIG SECTION HEADER
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) REPEAT(' ',30),'BERTHA-MBPT'
      WRITE(7, *) REPEAT(' ',30),'BERTHA-MBPT'
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',29),'Input summary'
      WRITE(7, *) REPEAT(' ',29),'Input summary'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     PRINT THE HAMILTONIAN OPTION
      WRITE(6, *) 'Hamiltonian type:',REPEAT(' ',51),HMLT
      WRITE(7, *) 'Hamiltonian type:',REPEAT(' ',51),HMLT
C
C     PRINT INPUT FILE NAME
      LF = LEN(TRIM(MOLCL))
      WRITE(6, *) 'Molecule name:',REPEAT(' ',58-LF),TRIM(MOLCL)
      WRITE(7, *) 'Molecule name:',REPEAT(' ',58-LF),TRIM(MOLCL)
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
      IF(READIN) THEN
        WRITE(6, *) 'SCF density matrix:',REPEAT(' ',46),'Read in'
        WRITE(7, *) 'SCF density matrix:',REPEAT(' ',46),'Read in'
      ELSE
        WRITE(6, *) 'SCF density matrix:',REPEAT(' ',44),'New start'
        WRITE(7, *) 'SCF density matrix:',REPEAT(' ',44),'New start'
      ENDIF
C
C     EQ-COEFFICIENT CALCULATION
      IF(EQFILE) THEN
        WRITE(6, *) 'Eq-coefficients:',REPEAT(' ',44),'Save to file'
        WRITE(7, *) 'Eq-coefficients:',REPEAT(' ',44),'Save to file'
      ELSE
        WRITE(6, *) 'Eq-coefficients:',REPEAT(' ',48),'By batch'
        WRITE(7, *) 'Eq-coefficients:',REPEAT(' ',48),'By batch'
      ENDIF
C
C     OPENMP PARALLEL OPTION
      IF(OPENMP) THEN
        WRITE(6, *) 'OpenMP parallel option:',REPEAT(' ',42),'Enabled'
        WRITE(7, *) 'OpenMP parallel option:',REPEAT(' ',42),'Enabled'
        WRITE(6, *) 'OpenMP thread limit:'   ,REPEAT(' ',40),NPRCSR
        WRITE(7, *) 'OpenMP thread limit:'   ,REPEAT(' ',40),NPRCSR
      ELSE
        WRITE(6, *) 'OpenMP parallel option:',REPEAT(' ',41),'Disabled'
        WRITE(7, *) 'OpenMP parallel option:',REPEAT(' ',41),'Disabled'
      ENDIF
C
C     SECTION FOR FILE NAMES
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     PRINT FILE OUTPUT NAMES
      LN = LEN(TRIM(OUTFL))
      WRITE(6, *) 'Output file name: ',REPEAT(' ',54-LN),TRIM(OUTFL)
      WRITE(7, *) 'Output file name: ',REPEAT(' ',54-LN),TRIM(OUTFL)
C
C     RECORD TIME AT BEGINNING OF CALCULATION
      CALL CPU_TIME(TCPUS(1))
      CALL TIMENOW(STAMP)
      WRITE(6, *) 'Time at BERTHA initiation:',REPEAT(' ',26),STAMP
      WRITE(7, *) 'Time at BERTHA initiation:',REPEAT(' ',26),STAMP
C
C     SECTION FOR PHYSICAL PARAMETERS
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     PRINT SPEED OF LIGHT
20    FORMAT(1X,A,43X,F14.10)
      WRITE(6,20) 'Speed of light:',CV
      WRITE(7,20) 'Speed of light:',CV
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
      INCLUDE 'parameters.h'
C
      INTEGER*16 NCMEM,NDMEM,NEMEM,NMMEM,NRMEM,NTMEM
      INTEGER*16 NHMEM
C
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VANM(MDM,MDM),
     &           VSLF(MDM,MDM),VUEH(MDM,MDM),VWKR(MDM,MDM),
     &           VKSB(MDM,MDM),QDIR(MDM,MDM),QXCH(MDM,MDM),
     &           WDIR(MDM,MDM),WXCH(MDM,MDM),CPLE(MDM,MDM)
C
      COMMON/DENS/DENC,DENO,DENT
      COMMON/EIGC/COEF
      COMMON/E0LL/E0LLFL(MFL,4),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,4),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EILS/EILSFL(MFL,12),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EISL/EISLFL(MFL,12),IADISL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VANM,VSLF,
     &            VUEH,VWKR,VKSB,QDIR,QXCH,WDIR,WXCH,CPLE
C
C     TITLE FOR MEMORY SUMMARY
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',24),'Stack memory allocation'
      WRITE(7, *) REPEAT(' ',24),'Stack memory allocation'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     APPROXIMATE THE MEMORY STORAGE ALLOCATION
      NCMEM = 2*SIZE(COEF)
      NDMEM = 2*SIZE(DENC) + 2*SIZE(DENO) + 2*SIZE(DENT)
      NEMEM = SIZE(E0LLFL) + SIZE(E0SSFL) + SIZE(EILSFL) + SIZE(EISLFL)
      NHMEM = 20*MFL
      NMMEM = 2*SIZE(FOCK) + 2*SIZE(OVLP) + 2*SIZE(HNUC) + 2*SIZE(HKIN)
     &      + 2*SIZE(GDIR) + 2*SIZE(GXCH) + 2*SIZE(BDIR) + 2*SIZE(BXCH)
     &      + 2*SIZE(VANM) + 2*SIZE(VSLF) + 2*SIZE(VUEH) + 2*SIZE(VWKR)
     &      + 2*SIZE(VKSB) + 2*SIZE(QDIR) + 2*SIZE(QXCH) + 2*SIZE(QDIR)
     &      + 2*SIZE(QXCH) + 2*SIZE(WDIR) + 2*SIZE(WXCH) + 2*SIZE(CPLE)
      NTMEM = NCMEM + NDMEM + NEMEM + NHMEM + NMMEM
C
C     SIZES (IN GIGABYTES)
      SCMEM = DFLOAT(NCMEM)*8.0D-9
      SDMEM = DFLOAT(NDMEM)*8.0D-9
      SEMEM = DFLOAT(NEMEM)*8.0D-9
      SHMEM = DFLOAT(NHMEM)*8.0D-9
      SMMEM = DFLOAT(NMMEM)*8.0D-9
      STMEM = DFLOAT(NTMEM)*8.0D-9
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
      WRITE(6,21) 'MTRX','Molecular Fock matrices    ',NMMEM,SMMEM
      WRITE(7,21) 'MTRX','Molecular Fock matrices    ',NMMEM,SMMEM
      WRITE(6,21) 'EQTT','Eq-coefficient file        ',NEMEM,SEMEM
      WRITE(7,21) 'EQTT','Eq-coefficient file        ',NEMEM,SEMEM
      WRITE(6,21) 'RCTT','R-integral class file      ',NHMEM,SHMEM
      WRITE(7,21) 'RCTT','R-integral class file      ',NHMEM,SHMEM
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
C                      EXIT SUMMARY ROUTINE FOR                        C
C                          ** B E R T H A **                           C
C -------------------------------------------------------------------- C
C  OUTPUT PRINTS A SUMMARY OF TOTAL CALCULATION DATA TO TERMINAL.      C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*16 HMS
      CHARACTER*20 STAMP
C
      COMMON/TCPU/TBEG,TNUC,TATM,TSCF,TMPT,TPRP,TPLT,TTOT
C
C     OVERALL CALCULATION TIME
      CALL CPU_TIME(TTOT)
      TTOT = TTOT-TBEG
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
      WRITE(6,20) 'Nuclear fitting:              ',HMS(TNUC)
      WRITE(7,20) 'Nuclear fitting:              ',HMS(TNUC)
      WRITE(6,20) 'Many-body perturbation theory:',HMS(TMPT)
      WRITE(7,20) 'Many-body perturbation theory:',HMS(TMPT)
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) 'Total CPU time:               ',HMS(TTOT)
      WRITE(7,20) 'Total CPU time:               ',HMS(TTOT)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     RECORD TIME AT BERTHA EXIT
      CALL TIMENOW(STAMP)
C
C     SUCCESSFUL EXIT MESSAGE
30    FORMAT(1X,A,21(' '),A)
      WRITE(6, *)
      WRITE(7, *)
      WRITE(6,30) 'Successful BERTHA-MBPT exit at:',STAMP
      WRITE(7,30) 'Successful BERTHA-MBPT exit at:',STAMP
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
      RETURN
      END
C
C
      BLOCK DATA PHYSPRM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  PPPPPPP  HH    HH YY    YY  SSSSSS  PPPPPPP  RRRRRRR  MM       MM   C
C  PP    PP HH    HH YY    YY SS    SS PP    PP RR    RR MMM     MMM   C
C  PP    PP HH    HH YY    YY SS       PP    PP RR    RR MMMM   MMMM   C
C  PP    PP HHHHHHHH  YY  YY   SSSSSS  PP    PP RR    RR MM MM MM MM   C
C  PPPPPPP  HH    HH   YYYY         SS PPPPPPP  RRRRRRR  MM  MMM  MM   C
C  PP       HH    HH    YY    SS    SS PP       RR    RR MM   M   MM   C
C  PP       HH    HH    YY     SSSSSS  PP       RR    RR MM       MM   C
C                                                                      C
C -------------------------------------------------------------------- C
C  BLOCK DATA FOR PHYSICAL CONSTANTS, CALCULATED TRANSCENDENTALS AND   C
C  PERIODIC TABLE ELEMENT NAMES.                                       C
C**********************************************************************C
C
      CHARACTER*2 ELMT(120)
C
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/MDLV/ELMT
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     PERIODIC TABLE ELEMENT NAMES
      DATA ELMT/'H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne',
     &          'Na','Mg','Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca',
     &          'Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn',
     &          'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y' ,'Zr',
     &          'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     &          'Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd',
     &          'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     &          'Lu','Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg',
     &          'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     &          'Pa','U' ,'Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     &          'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds',
     &          'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og','Ue','Un'/
C
C     PI AND ROOTS
      DATA PI,PI12,PI32,PI52/3.1415926535897932D0,1.7724538509055160D0,
     &                       5.5683279968317078D0,1.7493418327624863D1/
C
C     NATURAL LOGARITHM DATA
      DATA PILG,TWLG,THLG/1.1447298858494002D0,0.6931471805599453D0,
     &                                         1.0986122886681097D0/
C
C     NATURAL CONSTANTS
      DATA EULR/0.5772156649015329D0/
C
C     OTHER ROOTS
      DATA TW12/1.4142135623730950D0/
C
C     SPEED OF LIGHT (CODATA 2018 MEASUREMENT - UNCERTAINTY 0.000000021)
      DATA CV/137.035999084D0/
C
C     SPEED OF LIGHT (CODATA 2010 MEASUREMENT)
c     DATA CV/137.0359990D0/
C
C     SPEED OF LIGHT (CODATA 1999 MEASUREMENT)
C     DATA CV/137.03599976D0/
C
C     SPEED OF LIGHT (USED IN HAAKON'S THESIS)
c     DATA CV/137.0359898D0/
C
C     RELATIVE ELECTRON, MUON AND TAUON MASSES (ATOMIC UNITS)
      DATA EMSS,UMSS,TMSS/1.0000D0,206.7683D0,3477.1429D0/
C
C     RATIO OF PROTON MASS TO ELECTRON MASS AND PROTON RADIUS (IN fm)
      DATA PMSS,PRAD/1836.153D0,0.842D0/
C
C     COMPTON WAVELENGTH (1.0D0/CV)
      DATA CMPW/7.297353064D-3/
C
C     FREE ELECTRON G-FACTOR
      DATA GFREE/2.0023193043622D0/
C
C     FERMI COUPLING CONSTANT (IN eV^-2)
C     DATA GFRMI/1.1663787D-23/
C
C     FERMI COUPLING CONSTANT (IN au -- USED IN HAAKON'S THESIS)
C     DATA GFRMI/2.22D-14/
C
C     FERMI COUPLING CONSTANT (IN au)
      DATA GFRMI/2.22255D-14/
C
C     WEINBERG MIXING ANGLE SIN(θ_W)^2
      DATA WEIN/0.23120D0/
C
C     CONVERSION FACTORS: HARTREE TO HZ, HARTREE TO eV, HARTREE TO cm^-1
      DATA CHZ,CEV,CCM/6.579684D+15,2.7211332663D+01,2.1947422215D+05/
C
C     CONVERSION FACTORS: BOHR RADIUS TO FM OR ÅNGSTROM
      DATA CFM,CNG/5.29177211D+4,5.29177211D-1/
C
C     CONVERSION FACTOR: ELECTRON CHARGE ÅNGSTROM TO DEBYE (eÅ -> D)
      DATA CDB/2.54174623D+0/
C
C     http://www.fileformat.info/info/unicode/font/gnu_unifont/grid.htm
C     https://coolsymbol.com/
C     WRITE(6, *) 'αβγδεζηθικλμνξοπρϱςστυφχψω'
C     WRITE(6, *) 'ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡ  ΣΤΥΦΧΨΩ'
C     WRITE(6, *) 'ⅈℵℎℏℓℒℋℛℐℜℑℕℙℚℝℤℂ℘' 
C     WRITE(6, *) '×·∕−±∓≠∞—–…≪≫≤≥'
C     WRITE(6, *) '⅟½⅓⅕⅙⅛⅔⅖⅚⅜¾⅗⅝⅞⅘'
C     WRITE(6, *) 'π∞Σ√∛'∜∫∬∭∮∯∰∀∂∃∄∅∆∇∈∉∋∌∎∏∐∑'
C     WRITE(6, *) '∗∝∠∡∢∧∨∩∪∴∵∶∷∼∽∿⊕⊖⊗⊘⊙⊚⊛⊜⊝␢Å'
C     WRITE(6, *) '‗‘’‚‛“”„‟•‣❛❜❝❞¿'
C     WRITE(6, *) '✓✔✗✘☓√☑☐☒☼❄❆♤♠♧♣♡♥♢♦'
C     WRITE(6, *) '¢$€£¥®™☎⌨✁✂✎✏✐§☛'
C     WRITE(6, *) '↕↖↗↘↙↚↛↥↦↧↰↱↲↳↴↺↻↼↽↾↿⇀⇁⇂⇃⇄⇅⇆⇇⇠⇡⇢⏎▶➔➘➙➚➛➜➞↵⇑⇓'
C     WRITE(6, *) '★☆✡✦✧✩✪✰✢✣✤✥✱✲✳✴✵✶✷✸✹✺✻✼✽✾✿❀❁❂❃❇❈❉❊❋❄❆❅≛'
C     WRITE(6, *) 'äÄéÉ'
C
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
C  ▶ GAMLOG(N) = DLOG(GAMMA(N/2))                                      C
C  ▶ GAMHLF(N) = GAMMA(N/2)                                            C
C**********************************************************************C
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     STARTING VALUES
      GAMLOG(1) = PILG*0.5D0
      GAMLOG(2) = 0.0D0
      GAMHLF(1) = PI12
      GAMHLF(2) = 1.0D0
C
C     SEED VALUES FOR INCREMENT
      F1 = 0.5D0
      F2 = 1.0D0
C
C     FILL TABLE VALUES
      DO N=4,300,2
        GAMLOG(N-1) = GAMLOG(N-3)+DLOG(F1)
        GAMLOG(N  ) = GAMLOG(N-2)+DLOG(F2)
        GAMHLF(N-1) = GAMHLF(N-3)*F1
        GAMHLF(N  ) = GAMHLF(N-2)*F2
        F1 = F1+1.0D0
        F2 = F2+1.0D0
      ENDDO
C
      RETURN
      END
C
C
       SUBROUTINE FACTRL
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         FFFFFFFF   AA     CCCCCC TTTTTTTT RRRRRRR  LL                C
C         FF        AAAA   CC    CC   TT    RR    RR LL                C
C         FF       AA  AA  CC         TT    RR    RR LL                C
C         FFFFFF  AA    AA CC         TT    RR    RR LL                C
C         FF      AAAAAAAA CC         TT    RRRRRRR  LL                C
C         FF      AA    AA CC    CC   TT    RR    RR LL                C
C         FF      AA    AA  CCCCCC    TT    RR    RR LLLLLLLL          C
C                                                                      C
C -------------------------------------------------------------------- C
C  FACTRL GENERATES A SET OF N! AND N!! AS REAL NUMBERS.               C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C  ▶ RFACT - REGULAR FACTORIALS, RFACT(N) = N!                         C
C  ▶ SFACT - SEMI-FACTORIALS,    SFACT(N) = N!!                        C
C**********************************************************************C
C
      COMMON/FCTS/RFACT(0:80),SFACT(0:80)
C
      RFACT(0) = 1.0D0
      RFACT(1) = 1.0D0
      SFACT(0) = 1.0D0
      SFACT(1) = 1.0D0
      DO I=2,80
        RNUMBER  = DFLOAT(I)
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
C   [2] LABELS: MOLECULAR GEOMETRY AND FOCK MATRIX BLOCK ADDRESSES.    C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] NUCCOUL: NUCLEAR CHARGE DISTRIBUTION DETAILS.                  C
C   [B] VNFERMI: NORMALISED FERMI NUCLEAR POTENTIAL AT RADIUS R.       C
C   [C] POLYLOG: THE POLYLOGARITHM OF NEGATIVE EXPONENTIAL ARGUMENT.   C
C   [D] NUCGEOM: BOND DISTANCES AND NUCLEAR REPULSION ENERGY.          C
C   [E] FOCKIND: CALCULATE ADDRESSES OF FOCK MATRIX FOR BASIS QN'S.    C
C   [F] CARTIND: GENERATES INDICES FOR EQ-COEFFS AND R-INTEGRALS.      C
C   [G] AUFBAU: DETERMINES GROUND STATE ATOMIC ELECTRON CONFIG.        C
C   [H] SPECTRM0: ATOMIC SPECTRUM W/ EIGENVALUES AND RADIAL MOMENTS.   C
C   [I] SPECTRM: MOLECULAR SPECTRUM W/ EIGENVALUES AND TERM SYMBOLS.   C
C   [J] LLAB: GIVES THE CHARACTER CORRESPONDING TO LQN VALUE LQN.      C
C   [K] ROTATE: PERFORM TWO EULER ROTATIONS ON ALL ATOMIC CENTRES.     C
C   [L] MMPROD: PRODUCT OF TWO SQUARE ARRAYS OF DOUBLES.               C
C   [M] MVPROD: PRODUCT OF A SQUARE MATRIX AND VECTOR OF DOUBLES.      C
C   [N] HMS: RETURNS A QUOTED TIME IN SECONDS AS 'MIN-SEC'.            C
C   [O] MS: RETURNS A QUOTED TIME IN SECONDS AS 'HR-MIN-SEC'.          C
C   [P] TIMENOW: RETURNS A DATE STRING FOR THE CPU TIME.               C
C**********************************************************************C
C
C
      SUBROUTINE NUCCOUL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    NN    NN UU    UU  CCCCCC   CCCCCC   OOOOOO  UU    UU LL          C
C    NNN   NN UU    UU CC    CC CC    CC OO    OO UU    UU LL          C
C    NNNN  NN UU    UU CC       CC       OO    OO UU    UU LL          C
C    NN NN NN UU    UU CC       CC       OO    OO UU    UU LL          C
C    NN  NNNN UU    UU CC       CC       OO    OO UU    UU LL          C
C    NN   NNN UU    UU CC    CC CC    CC OO    OO UU    UU LL          C
C    NN    NN  UUUUUU   CCCCCC   CCCCCC   OOOOOO   UUUUUU  LLLLLLLL    C
C                                                                      C
C -------------------------------------------------------------------- C
C  NUCCOUL CONTROLS THE CHARGE POTENTIALS THAT APPLY TO EACH NUCLEAR   C
C  CENTRE, GIVEN A RANGE OF NUCLEAR CHARGE MODELS AND QED OPTIONS.     C
C -------------------------------------------------------------------- C
C  THIS VERSION OF THE ROUTINE READS IN RESULTS FROM BERTHA-NUCLEAR.   C
C -------------------------------------------------------------------- C
C  ▶ NCOUL = 0 (POINT), NCOUL = 1 (GAUSS), NCOUL > 1 (FERMI).          C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL FILETHERE
      CHARACTER*2 ELMT(120)
      CHARACTER*5 NMDL
      CHARACTER*40 MOLCL,WFNFL,OUTFL,NUCFL
C
      DIMENSION RSQ(MCT),NCOUL(MCT)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/MDLV/ELMT
C
C     NUCLEUS FILE NAME
      NUCFL = 'output/'//TRIM(MOLCL)//'_nuclear.dat'
C
C     TRY TO READ THE FILE IN
      INQUIRE(FILE=TRIM(NUCFL),EXIST=FILETHERE)
C
C     READ IN AN EXISTING FILE
      IF(FILETHERE) THEN
        OPEN(UNIT=8,FILE=TRIM(NUCFL),STATUS='UNKNOWN')
C
C         LOOP OVER NUCLEAR CENTERS
          DO IZ=1,NCNT
C
C           NUMBER OF FITTING FUNCTIONS AND RMS RADIUS
            READ(8, *) NNUC(IZ),RNUC(IZ)
C
C           AMPLITUDES AND PARAMETERS FOR EACH NUCLEAR FITTING FUNCTION
            DO IFT=0,NNUC(IZ)
              READ(8, *) FNUC(IZ,IFT),XNUC(IZ,IFT)
            ENDDO
C
C           FERMI MODEL
            TFMI(IZ) = 2.30D0/CFM
            AFMI(IZ) = 0.25D0*TFMI(IZ)/THLG
C
          ENDDO
C
C       CLOSE DATA FILE
        CLOSE(UNIT=8)
C
C     CANNOT OPEN A FILE SO ASSUME GAUSSIAN
      ELSE
C
C       ALERT THE USER THAT THE NUCLEAR FILE COULD NOT BE FOUND
        WRITE(6, *) 'In NUCCOUL: could not find nuclear data.'
        WRITE(7, *) 'In NUCCOUL: could not find nuclear data.'
C
C       GAUSSIAN FITTING FOR ALL NUCLEI
        DO IZ=1,NCNT
C
C         NUCLEAR RADIUS
          IF(INT(ZNUC(IZ)).EQ.1) THEN
C           PROTON HAS A VERY SMALL RADIUS (IN FM)
            RFM = PRAD
          ELSEIF(ZNUC(IZ).GT.1.0D0.AND.ZNUC(IZ).LT.90.0D0) THEN
C           EMPIRICAL FORMULA FOR NUCLEAR RMS RADIUS IF Z < 90 (IN FM)
            RFM = 0.836D0*(ANUC(IZ)**(1.0D0/3.0D0)) + 0.57D0
C           RFM = 1.20D0*(ANUC(IZ)**(1.0D0/3.0D0))
          ELSE
C           EMPIRICAL FORMULA FOR NUCLEAR RMS RADIUS IF Z > 90 (IN FM)
            RFM = 0.770D0*(ANUC(IZ)**(1.0D0/3.0D0)) + 0.98D0
          ENDIF
C
C         CONVERT RESULT TO ATOMIC UNITS
          RNUC(IZ) = RFM/CFM
C
C         GAUSSIAN MODEL
          NNUC(IZ)   = 0
          FNUC(IZ,0) = 1.0D0
          XNUC(IZ,0) = 1.5D0/(RNUC(IZ)*RNUC(IZ))
C
C         FERMI MODEL
          TFMI(IZ) = 2.30D0/CFM
          AFMI(IZ) = 0.25D0*TFMI(IZ)/THLG
C
        ENDDO
C
C       ALERT THE USER THAT NUCLEI ARE GAUSSIAN NOW
21      FORMAT(26X,A)
22      FORMAT(1X,A,2X,'|',3X,A,4X,A,5X,A,1X,'|',1X,A,3X,A,3X,A)
23      FORMAT(1X,I2,' (',A,')',1X,'|',5X,I3,4X,F8.4,5X,F6.4,1X,'|',
     &                                             4X,I2,7X,A,2X,F10.8)
        WRITE(6, *) ' '
        WRITE(7, *) ' '
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6,21) 'Nuclear charge summary'
        WRITE(7,21) 'Nuclear charge summary'
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6,22) 'Centre','Z (e)','A (m_p+)','R (fm)',
     &              'N_fit','Potential','R-squared'
        WRITE(7,22) 'Centre','Z (e)','A (m_p+)','R (fm)',
     &              'N_fit','Potential','R-squared'
        WRITE(6, *) REPEAT('-',72)
        WRITE(7, *) REPEAT('-',72)
        DO IZ=1,NCNT
          IZNC = INT(ZNUC(IZ))
          WRITE(6,23) IZ,ELMT(IZNC),IZNC,ANUC(IZ),RNUC(IZ)*CFM,NNUC(IZ),
     &                                                 NMDL(IZ),RSQ(IZ)
          WRITE(7,23) IZ,ELMT(IZNC),IZNC,ANUC(IZ),RNUC(IZ)*CFM,NNUC(IZ),
     &                                                 NMDL(IZ),RSQ(IZ)
        ENDDO
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6, *) ' '
        WRITE(7, *) ' '
C
      ENDIF
C
      RETURN
      END
C
C
C
C
      SUBROUTINE NUCVCPL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    NN    NN UU    UU  CCCCCC  VV    VV  CCCCCC  PPPPPPP  LL          C
C    NNN   NN UU    UU CC    CC VV    VV CC    CC PP    PP LL          C
C    NNNN  NN UU    UU CC       VV    VV CC       PP    PP LL          C
C    NN NN NN UU    UU CC       VV    VV CC       PP    PP LL          C
C    NN  NNNN UU    UU CC        VV  VV  CC       PPPPPPP  LL          C
C    NN   NNN UU    UU CC    CC   VVVV   CC    CC PP       LL          C
C    NN    NN  UUUUUU   CCCCCC     VV     CCCCCC  PP       LLLLLLLL    C
C                                                                      C
C -------------------------------------------------------------------- C
C  NUCVCPL CONTROLS THE VACUUM POLARISATION POTENTIALS FOR ALL ATOMIC  C
C  CENTRES, GIVEN A RANGE OF NUCLEAR CHARGE MODELS AND QED OPTIONS.    C
C -------------------------------------------------------------------- C
C  THIS VERSION OF THE ROUTINE READS IN RESULTS FROM BERTHA-NUCLEAR.   C
C -------------------------------------------------------------------- C
C  ▶ NFT > 26 (FOR WHICHEVER NUCLEAR CHARGE MODEL.)                    C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL FILETHERE
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*60 VACFL
C
      COMMON/BUEH/RUEH(MCT,3),FUEH(MCT,MFT),XUEH(MCT,MFT),NUEH(MCT)
      COMMON/BKSB/RKSB(MCT,3),FKSB(MCT,MFT),XKSB(MCT,MFT),NKSB(MCT)
      COMMON/BWKR/RWKR(MCT,3),FWKR(MCT,MFT),XWKR(MCT,MFT),NWKR(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/QUEH/RU(0:NRAD),VU(MCT,0:NRAD),R0U,RMU,RIU,NUU,NEU
      COMMON/QWKR/RW(0:NRAD),VW(MCT,0:NRAD),R0W,RMW,RIW,NUW,NEW
      COMMON/QKSB/RK(0:NRAD),VK(MCT,0:NRAD),R0K,RMK,RIK,NUK,NEK
C
C     UEHLING POTENTIAL
      IF(HMLT.NE.'DHFQ'.AND.HMLT.NE.'DHFP') RETURN
C
C     UEHLING FITTING SET
      VACFL = 'output/'//TRIM(MOLCL)//'_vpueh.dat'
      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
      IF(FILETHERE) THEN
C       READ FITTING DATA IN FROM AN EXISTING FILE
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          DO IZ=1,NCNT
            READ(8, *) NUEH(IZ)
            DO IFT=1,NUEH(IZ)
              READ(8, *) FUEH(IZ,IFT),XUEH(IZ,IFT)
            ENDDO
          ENDDO
        CLOSE(UNIT=8)
      ELSE
C       ALERT THE USER THAT THE UEHLING FILE COULD NOT BE FOUND
        WRITE(6, *) 'In NUCVCPL: could not find Uehling data.'
        WRITE(7, *) 'In NUCVCPL: could not find Uehling data.'
        DO IZ=1,NCNT
          NUEH(IZ) = 0
        ENDDO
      ENDIF
CC
CC     UEHLING RADIAL GRID DATA
C      VACFL = 'output/'//TRIM(MOLCL)//'_vpueh_grid.dat'
C      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
C      IF(FILETHERE) THEN
CC       READ FITTING DATA IN FROM AN EXISTING FILE
C        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
C          READ(8, *) R0U,RMU,RIU,NUU,NEU
C          DO N=0,NRAD
C            READ(8, *) RU(N),(VU(IZ,N),IZ=1,NCNT)
C          ENDDO
C        CLOSE(UNIT=8)
C      ELSE
CC       ALERT THE USER THAT THE UEH GRID FILE COULD NOT BE FOUND
C        WRITE(6, *) 'In NUCVCPL: could not find UEH grid data.'
C        WRITE(7, *) 'In NUCVCPL: could not find UEH grid data.'
C        DO N=0,NRAD
C          RU(N) = 0.0D0
C          DO IZ=1,NCNT
C            VU(IZ,N) = 0.0D0
C          ENDDO
C        ENDDO
C      ENDIF
C
C     WICHMANN-KROLL FITTING SET
      VACFL = 'output/'//TRIM(MOLCL)//'_vpwkr.dat'
      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
      IF(FILETHERE) THEN
C       READ FITTING DATA IN FROM AN EXISTING FILE
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          DO IZ=1,NCNT
            READ(8, *) NWKR(IZ)
            DO IFT=1,NWKR(IZ)
              READ(8, *) FWKR(IZ,IFT),XWKR(IZ,IFT)
            ENDDO
          ENDDO
        CLOSE(UNIT=8)
      ELSE
C       ALERT THE USER THAT THE WICHMANN-KROLL FILE COULD NOT BE FOUND
        WRITE(6, *) 'In NUCVCPL: could not find Wichmann-Kroll data.'
        WRITE(7, *) 'In NUCVCPL: could not find Wichmann-Kroll data.'
        DO IZ=1,NCNT
          NWKR(IZ) = 0
        ENDDO
      ENDIF
CC
CC     WICHMANN-KROLL RADIAL GRID DATA
C      VACFL = 'output/'//TRIM(MOLCL)//'_vpwkr_grid.dat'
C      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
C      IF(FILETHERE) THEN
CC       READ FITTING DATA IN FROM AN EXISTING FILE
C        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
C          READ(8, *) R0W,RMW,RIW,NUW,NEW
C          DO N=0,NRAD
C            READ(8, *) RW(N),(VW(IZ,N),IZ=1,NCNT)
C          ENDDO
C        CLOSE(UNIT=8)
C      ELSE
CC       ALERT THE USER THAT THE WKR GRID FILE COULD NOT BE FOUND
C        WRITE(6, *) 'In NUCVCPL: could not find WKR grid data.'
C        WRITE(7, *) 'In NUCVCPL: could not find WKR grid data.'
C        DO N=0,NRAD
C          RW(N) = 0.0D0
C          DO IZ=1,NCNT
C            VW(IZ,N) = 0.0D0
C          ENDDO
C        ENDDO
C      ENDIF
C
C     KÄLLÉN-SABRY FITTING SET
      VACFL = 'output/'//TRIM(MOLCL)//'_vpksb.dat'
      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
      IF(FILETHERE) THEN
C       READ FITTING DATA IN FROM AN EXISTING FILE
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          DO IZ=1,NCNT
            READ(8, *) NKSB(IZ)
            DO IFT=1,NKSB(IZ)
              READ(8, *) FKSB(IZ,IFT),XKSB(IZ,IFT)
            ENDDO
          ENDDO
        CLOSE(UNIT=8)
      ELSE
C       ALERT THE USER THAT THE KÄLLÉN-SABRY FILE COULD NOT BE FOUND
        WRITE(6, *) 'In NUCVCPL: could not find Källén-Sabry data.'
        WRITE(7, *) 'In NUCVCPL: could not find Källén-Sabry data.'
        DO IZ=1,NCNT
          NKSB(IZ) = 0
        ENDDO
      ENDIF
CC
CC     KÄLLÉN-SABRY RADIAL GRID DATA
C      VACFL = 'output/'//TRIM(MOLCL)//'_vpwkr_grid.dat'
C      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
C      IF(FILETHERE) THEN
CC       READ FITTING DATA IN FROM AN EXISTING FILE
C        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
C          READ(8, *) R0K,RMK,RIK,NUK,NEK
C          DO N=0,NRAD
C            READ(8, *) RK(N),(VK(IZ,N),IZ=1,NCNT)
C          ENDDO
C        CLOSE(UNIT=8)
C      ELSE
CC       ALERT THE USER THAT THE KSB GRID FILE COULD NOT BE FOUND
C        WRITE(6, *) 'In NUCVCPL: could not find KSB grid data.'
C        WRITE(7, *) 'In NUCVCPL: could not find KSB grid data.'
C        DO N=0,NRAD
C          RW(N) = 0.0D0
C          DO IZ=1,NCNT
C            VW(IZ,N) = 0.0D0
C          ENDDO
C        ENDDO
C      ENDIF
C
      RETURN
      END
C
C
      FUNCTION RHONUC(MODEL,IZ,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        RRRRRRR  HH    HH  OOOOOO  NN    NN UU    UU  CCCCCC          C
C        RR    RR HH    HH OO    OO NNN   NN UU    UU CC    CC         C
C        RR    RR HH    HH OO    OO NNNN  NN UU    UU CC               C
C        RR    RR HHHHHHHH OO    OO NN NN NN UU    UU CC               C
C        RRRRRRR  HH    HH OO    OO NN  NNNN UU    UU CC               C
C        RR    RR HH    HH OO    OO NN   NNN UU    UU CC    CC         C
C        RR    RR HH    HH  OOOOOO  NN    NN  UUUUUU   CCCCCC          C
C                                                                      C
C -------------------------------------------------------------------- C
C  RHONUC EVALUATES THE VALUE OF THE NORMALISED NUCLEAR CHARGE DENSITY C
C  FOR NUCLEUS IZ AT RADIUS R, BASED ON A NUCLEAR CHARGE MODEL.        C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5 NMDL,MODEL
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     RADIUS MUST BE NON-NEGATIVE
      IF(R.LT.0.0D0) THEN
        WRITE(6, *) 'In RHONUC: radius R must be non-negative. R = ',R
        WRITE(7, *) 'In RHONUC: radius R must be non-negative. R = ',R
        STOP
      ENDIF
C
C     REASONS TO USE GAUSSIAN NUCLEAR CHARGE MODEL
      EPS = DSQRT(1.4D0)*PI*0.25D0*TFMI(IZ)/THLG
      IF(EPS.GT.RNUC(IZ)) THEN
C
C       NORMALISATION FACTOR
        XI = 3.0D0/(2.0D0*RNUC(IZ)*RNUC(IZ))
        YN = XI*DSQRT(XI)/PI32
        RHONUC = YN*DEXP(-R*R*XI)
        RETURN
C
      ENDIF
C
C     GAUSSIAN NUCLEAR CHARGE MODEL
      IF(MODEL.EQ.'GAUSS') THEN
C
C       NORMALISATION FACTOR
        ZT = 3.0D0/(2.0D0*RNUC(IZ)*RNUC(IZ))
        YN = ZT*DSQRT(ZT)/(PI*PI12)
        RHONUC = YN*DEXP(-R*R*ZT)
C      
      ENDIF
C
C     FERMI NUCLEAR CHARGE MODEL
      IF(MODEL.EQ.'FERMI') THEN
C
        A  = AFMI(IZ)
        C = 5.0D0*RNUC(IZ)*RNUC(IZ)/3.0D0
     &    - 7.0D0*PI*PI*AFMI(IZ)*AFMI(IZ)/3.0D0
        C = DSQRT(C)
        U  = A/C
        C3 = C*C*C
C
C       NORMALISATION CONSTANT
        Y1 = 1.0D0
        Y2 = PI*PI*U*U
        Y3 =-6.0D0*U*U*U*POLYLOG(3,-1.0D0/U)
        YN = 1.0D0/(Y1+Y2+Y3)
        YN = 3.0D0*YN/(4.0D0*PI*C3)
C
C       EXACT FERMI DENSITY FORMULA
        PR = 1.0D0 + DEXP((R-C)/A)
        RHONUC = YN/PR
C
      ENDIF
C
C     HOMOGENEOUS BALL MODEL
      IF(MODEL.EQ.'UNIFM') THEN
C
        B = DSQRT(5.0D0/3.0D0)*RNUC(IZ)
        IF(R.LE.B) THEN
          RHONUC = 0.75D0/(PI*B*B*B)
        ELSE
          RHONUC = 0.00D0
        ENDIF
C
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION VNFERMI(R,C,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    VV    VV NN    NN FFFFFFFF EEEEEEEE RRRRRRR  MM       MM IIII     C
C    VV    VV NNN   NN FF       EE       RR    RR MMM     MMM  II      C
C    VV    VV NNNN  NN FF       EE       RR    RR MMMM   MMMM  II      C
C    VV    VV NN NN NN FFFFFF   EEEEEE   RR    RR MM MM MM MM  II      C
C     VV  VV  NN  NNNN FF       EE       RRRRRRR  MM  MMM  MM  II      C
C      VVVV   NN   NNN FF       EE       RR    RR MM   M   MM  II      C
C       VV    NN    NN FF       EEEEEEEE RR    RR MM       MM IIII     C
C                                                                      C
C -------------------------------------------------------------------- C
C  VNFERMI EVALUATES THE NORMALISED SCALAR POTENTIAL ARISING FROM THE  C
C  TWO-PARAMETER FERMI NUCLEAR CHARGE DISTRIBUTION AT RADIUS R, GIVEN  C
C  HALF-DENSITY RADIUS C AND THE NUCLEAR SKIN THICKNESS T.             C
C  THE POTENTIAL IS EVALUATED IN CASES DEPENDING ON R AND C, AND       C
C  INVOLVES THE USE OF AN AUXILLIARY FUNCTION POLYLOG AS WELL.         C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     RADIUS MUST BE NON-NEGATIVE
      IF(R.LT.0.0D0) THEN
        WRITE(6, *) 'In VNFERMI: radius R must be non-negative. R = ',R
        WRITE(7, *) 'In VNFERMI: radius R must be non-negative. R = ',R
        STOP
      ENDIF
C
C     WARN THE USER IF C<A BUT DON'T STOP
      IF(C.LT.0.25D0*T/THLG) THEN
        WRITE(6, *) 'In VNFERMI: A/C > 1 so expressions diverge.'
        WRITE(7, *) 'In VNFERMI: A/C > 1 so expressions diverge.'
      ENDIF
C
C     THESE PARAMETERS MAKE FOR NEATER EXPRESSIONS
      A   = 0.25D0*T/THLG
      PIA = PI*A
      QAC = A/C
      PAC = PI*A/C
      C3  = C*C*C
C
C     NORMALISATION CONSTANT
      Y1 = 1.0D0
      Y2 = PAC*PAC
      Y3 = 6.0D0*POLYLOG(3,-C/A)*(QAC**3)
      YN = 1.0D0/(Y1+Y2-Y3)
C
C     POTENTIAL DEPENDS ON RADIUS
      IF(R.LE.C) THEN
        V1 = 1.5D0/C
        V2 = 0.5D0*R*R/C3
        V3 = 0.5D0*PIA*PIA/C3
        V4 = 3.0D0*A*A*POLYLOG(2,(R-C)/A)/C3
        IF(R.LT.1.0D-10) THEN
          V5 = (6.0D0*A*A*POLYLOG(2,-C/A) + 3.0D0*A*R*POLYLOG(1,-C/A)
     &                                        + R*R*POLYLOG(0,-C/A))/C3
          V6 = 0.0D0
        ELSE
          V5 = 6.0D0*(QAC**3)*POLYLOG(3,(R-C)/A)/R
          V6 = 6.0D0*(QAC**3)*POLYLOG(3,  -C /A)/R
        ENDIF
        VR = V1-V2+V3-V4+V5-V6
      ELSEIF(R.GT.C) THEN
        V1 = 1.0D0
        V2 = PAC*PAC
        V3 = 0.0D0
        V4 = 3.0D0*R*A*A*POLYLOG(2,-(R-C)/A)/C3
        V5 = 6.0D0*(QAC**3)*POLYLOG(3,-(R-C)/A)
        V6 = 6.0D0*(QAC**3)*POLYLOG(3,   -C /A)
        VR = (V1+V2+V3+V4+V5-V6)/R
      ENDIF
C
C     FINAL VALUE OF POTENTIAL
      VNFERMI =-YN*VR
C
      RETURN
      END
C
C
      FUNCTION POLYLOG(K,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     PPPPPPP   OOOOOO  LL      YY    YY LL       OOOOOO   GGGGGG      C
C     PP    PP OO    OO LL      YY    YY LL      OO    OO GG    GG     C
C     PP    PP OO    OO LL      YY    YY LL      OO    OO GG           C
C     PP    PP OO    OO LL       YY  YY  LL      OO    OO GG           C
C     PPPPPPP  OO    OO LL        YYYY   LL      OO    OO GG   GGG     C
C     PP       OO    OO LL         YY    LL      OO    OO GG    GG     C
C     PP        OOOOOO  LLLLLLLL   YY    LLLLLLLL OOOOOO   GGGGGG      C
C                                                                      C
C -------------------------------------------------------------------- C
C  POLYLOG EVALUATES THE POLYLOGARITHM FUNCTION OF ORDER K AND         C
C  ARGUMENT -DEXP(X) USING A TRUNCATED SERIES EXPANSION OF LENGTH NTR. C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
C     TRUNCATION OPTION SHOULD BE GREATER THAN 1
      IF(NTR.LE.1) THEN
        WRITE(6, *) 'In POLYLOG: series trunction too low. NTR = ',NTR
        WRITE(7, *) 'In POLYLOG: series trunction too low. NTR = ',NTR
        STOP
      ENDIF
C
C     ARGUMENT X MUST BE NEGATIVE
      IF(X.GT.0.0D0) THEN
        WRITE(6, *) 'In POLYLOG: argument X must be negative. X = ',X
        WRITE(7, *) 'In POLYLOG: argument X must be negative. X = ',X
        STOP
      ENDIF
C
C     INITIAL VALUES
      PHS = 1.0D0
      VAL = 0.0D0
      DO N=1,NTR
        PHS =-PHS
        RNK = DFLOAT(N)**K
        VAL = VAL + PHS*DEXP(X*N)/RNK
      ENDDO
C
C     VALUE OF POLYLOG
      POLYLOG = VAL
C
      RETURN
      END
C
C
      FUNCTION CLEBSCH(LQN,MQN,NSGN)
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
C  CLEBSCH DETERMINES A CLEBSCH-GORDAN COEFFICIENT FOR LQN, MQN AND    C
C  NSGN -- THESE ARE THE COEFFICIENTS IN FRONT OF THE SPIN-ANGULAR     C
C  FUNCTIONS FOR A GIVEN COMPONENT TYPE.                               C
C**********************************************************************C
C
      T1 = 2.0D0*LQN + 1.0D0 + MQN*NSGN
      T2 = 4.0D0*LQN + 2.0D0
C
      CLEBSCH = DSQRT(T1/T2)
C
      RETURN
      END
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
C  ▶ THE TSYM PACKAGE (WERNER 1996) CAN BE DIRECTLY IMPLEMENTED HERE - C
C    IT PERFORMS POINT GROUP SYMMETRY ANALYSIS. THIS ALLOWS A BLOCK    C
C    STRUCTURE ACCORDING TO IRREPS FOR THE HAMILTONIAN MATRIX (FASTER).C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL TRANSFORM
C
      CHARACTER*2 ELMT(120),ELA,ELB,ELC
      CHARACTER*4 UCHR
      CHARACTER*5 NMDL
      CHARACTER*8 SPCES
C
      DIMENSION XYZ(3,MCT),DIST(MCT),CENT(3),IZAD(MCT,3)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EWDR,EWXC,EANM,ESLF,EUEH,
     &            EWKR,EKSB
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/MDLV/ELMT
C
C     CHOICE OF UNITS -- 1 FOR AU, 2 FOR ÅNGSTROM
      IUNIT = 1
C
      IF(IUNIT.EQ.1) THEN
        UCHR = '(a0)'
        CFCT = 1.0D0
      ELSEIF(IUNIT.EQ.2) THEN
        UCHR = '(Å)'
        CFCT = CNG
      ENDIF
C
C     PRINT A BIG SECTION HEADER
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) REPEAT(' ',19),'NUCLEAR COORDINATES AND POTENTIALS'
      WRITE(7, *) REPEAT(' ',19),'NUCLEAR COORDINATES AND POTENTIALS'
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C     LOGICAL OPTION TO DISABLE ROTATION/TRANSLATION
      TRANSFORM = .TRUE.
C
      IF(.NOT.TRANSFORM) GOTO 500
C
C**********************************************************************C
C     IDENTIFY LARGEST ATOMIC CENTRES                                  C
C**********************************************************************C
C
C     IDENTIFY THE THREE LARGEST NUCLEAR CHARGES
      LZ1 = 0
      LZ2 = 0
      LZ3 = 0
      DO N=1,NCNT
        IZ = INT(ZNUC(N))
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
C     COUNT CENTRES WITH TOP THREE NUCLEAR CHARGES
      NZ1 = 0
      NZ2 = 0
      NZ3 = 0
      DO N=1,NCNT
        IZ = INT(ZNUC(N))
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
C     TRANSLATION: ORIGIN COINCIDES WITH CENTRE OF ALL HEAVY ATOMS     C
C**********************************************************************C
C
C     CALCULATE CENTRE OF ALL CHARGES LZ1
      DCNT = 0.0D0
      DO I=1,3
        CENT(I) = 0.0D0
        DO N=1,NZ1
          CENT(I) = CENT(I) + BXYZ(I,IZAD(N,1))
        ENDDO
        CENT(I) = CENT(I)/NZ1
        DCNT    = DCNT + CENT(I)**2
      ENDDO
      DCNT = DSQRT(DCNT)
C
C     TRANSLATE MOLECULE BY CENT(I) AND CALCULATE DISTANCE TO ORIGIN
      DO N=1,NCNT
        DO I=1,3
          XYZ(I,N) = BXYZ(I,N) - CENT(I)
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
C     ONE-CENTRE (MUST BE ATOMIC)
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
C     FIRST ROTATION: FIX ON A ROTATION CENTRE AND ROTATE TO Z-AXIS.   C
C**********************************************************************C
C
C     IF AN LZ1 CENTRE IS ON Z-AXIS BUT NOT THE ORIGIN, SKIP ROTATION
      DO NZ=1,NZ1
        N = IZAD(NZ,1)
        X = DIST(N) - DABS(XYZ(3,N))
        IF(DABS(X).LT.1.0D-6.AND.DIST(N).GT.1.0D-4) THEN
          GOTO 200
        ENDIF
      ENDDO
C
C     FIX ON AN LZ1 CENTRE THAT DOES NOT LIE ON ORIGIN
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
C     IF WE HAVE MADE IT THIS FAR, THERE IS NO LZ1 CENTRE ALREADY ON
C     THE Z-AXIS AND NO LZ1 CENTRE OFF THE ORIGIN --> MOVE TO LZ2.
C
C     VECTOR CENTRE OF ALL LZ2
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
C     IF CENTRE OF LZ2 IS ON Z-AXIS BUT NOT ORIGIN, SKIP ROTATION
      X = DCNT - DABS(CENT(3))
      IF(DABS(X).LT.1.0D-4.AND.DCNT.GT.1.0D-4) THEN
        GOTO 200
      ENDIF
C
C     IF CENTRE OF LZ2 IS ON ORIGIN, FIX ON AN LZ2 NOT ON THE ORIGIN
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
C     SKIP POINT (ROTATION CENTRE HAS BEEN IDENTIFIED)
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
C     TWO-CENTRE (MUST BE DIATOMIC)
      IF(NCNT.EQ.2) THEN
C
C       DESIGNATE SHAPE
        SHAPE = 'DIATOM'
C
C       NO FURTHER ROTATIONS NECESSARY
        GOTO 300
C
      ENDIF
C
C**********************************************************************C
C     SECOND ROTATION: FIX ON A ROTATION CENTRE AND ROTATE TO Y-AXIS.  C
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
C     IF AN LZ1 CENTRE IS ON YZ-AXIS BUT NOT THE ORIGIN, SKIP ROTATION
      DO NZ=1,NZ1
        N = IZAD(NZ,1)
        X = DIST(N) - DSQRT(XYZ(2,N)**2 + XYZ(3,N)**2)
        IF(DABS(X).LT.1.0D-6.AND.DIST(N).GT.1.0D-4) THEN
          GOTO 300
        ENDIF
      ENDDO
C
C     FIX ON AN LZ1 CENTRE THAT DOES NOT LIE ON ORIGIN OR Z-AXIS
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
C     IF WE HAVE MADE IT THIS FAR, THERE IS NO LZ1 CENTRE ALREADY ON
C     THE YZ-AXIS AND NO LZ1 CENTRE OFF THE ORIGIN --> MOVE TO LZ2.
C
C     IF THERE ARE NO LZ2 CENTRES, SKIP ROTATION
      IF(NZ2.EQ.0) THEN
        GOTO 300
      ENDIF
C
C     VECTOR CENTRE OF ALL LZ2
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
C     IF CENTRE OF LZ2 IS ON YZ-AXIS BUT NOT ORIGIN, SKIP ROTATION
      X = DCNT - DSQRT(CENT(2)**2 + CENT(3)**2)
      IF(DABS(X).LT.1.0D-3.AND.DCNT.GT.1.0D-3) THEN
        GOTO 300
      ENDIF
C
C     IF CENTRE OF LZ2 IS NOT ON Y-AXIS, MAKE IT THE CENTRE
      X = DCNT - DABS(CENT(2))
      IF(DABS(X).GT.1.0D-3.AND.DCNT.GT.1.0D-3) THEN
        GOTO 50
      ENDIF
C
C     IF CENTRE OF LZ2 IS ON Y-AXIS, FIX ON ONE OF THEM
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
C     IF CENTRE OF LZ2 IS ON ORIGIN, FIX ON ONE THAT IS OFF THE Z AXIS
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
C     IF WE HAVE MADE IT THIS FAR, THERE IS NO LZ2 CENTRE ALREADY ON
C     THE YZ-AXIS AND NO LZ1 CENTRE OFF THE ORIGIN --> MOVE TO LZ3.
C
C     IF THERE ARE NO LZ3 CENTRES, SKIP ROTATION
      IF(NZ3.EQ.0) THEN
        GOTO 300
      ENDIF
C
C     VECTOR CENTRE OF ALL LZ3
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
C     IF CENTRE OF LZ3 IS ON YZ-AXIS BUT NOT ORIGIN, SKIP ROTATION
      X = DCNT - DSQRT(CENT(2)**2 + CENT(3)**2)
      IF(DABS(X).LT.1.0D-3.AND.DCNT.GT.1.0D-3) THEN
        GOTO 300
      ENDIF
C
C     IF CENTRE OF LZ3 IS NOT ON Y-AXIS, MAKE IT THE CENTRE
      X = DCNT - DABS(CENT(2))
      IF(DABS(X).GT.1.0D-3.AND.DCNT.GT.1.0D-3) THEN
        GOTO 50
      ENDIF
C
C     IF CENTRE OF LZ3 IS ON Y-AXIS, FIX ON ONE OF THEM
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
C     IF CENTRE OF LZ3 IS ON ORIGIN, FIX ON ONE THAT IS OFF THE Z AXIS
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
C     SKIP POINT (ROTATION CENTRE HAS BEEN IDENTIFIED)
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
C     IF MOLECULE ALREADY HAS LABEL, SKIP
      IF(SHAPE.EQ.'ATOMIC') GOTO 400
      IF(SHAPE.EQ.'DIATOM'.OR.SHAPE.EQ.'LINEAR') GOTO 400
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
400   CONTINUE
C
C     TRANSFER ALL TEMPORARY COORDINATES TO THE COMMON ARRAY
      DO N=1,NCNT
        DO I=1,3
          BXYZ(I,N) = XYZ(I,N)
        ENDDO
      ENDDO
C
C     ROTATION/TRANSLATION SKIP
500   CONTINUE
C
C**********************************************************************C
C     MOLECULAR GEOMETRY, BOND DISTANCES AND NUCLEAR REPULSION ENERGY  C
C**********************************************************************C
C
C     ATOMIC COORDINATES
20    FORMAT(11X,A,1X,A)
21    FORMAT(1X,A,12X,A,14X,A,14X,A)
22    FORMAT(1X,I2,' (',A,') ',6X,F16.10,5X,F16.10,5X,F16.10)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,20) '  Molecular geometry A: Cartesian coordinates',UCHR
      WRITE(7,20) '  Molecular geometry A: Cartesian coordinates',UCHR
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,21) 'Centre','x-coord','y-coord','z-coord'
      WRITE(7,21) 'Centre','x-coord','y-coord','z-coord'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO IZ=1,NCNT
        ELA = ELMT(INT(ZNUC(IZ)))
        WRITE(6,22) IZ,ELA,(CFCT*BXYZ(IX,IZ),IX=1,3)
        WRITE(7,22) IZ,ELA,(CFCT*BXYZ(IX,IZ),IX=1,3)
      ENDDO
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C     BOND ANGLES AND DISTANCES
30    FORMAT(21X,A)
31    FORMAT(1X,A,6X,A,1X,A,11X,A,12X,A)
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
        WRITE(6,31) 'C1  C2','Bond length',UCHR,'C1  C2  C3',
     &                                                    'Angle (deg)'
        WRITE(7,31) 'C1  C2','Bond length',UCHR,'C1  C2  C3',
     &                                                    'Angle (deg)'
        WRITE(6, *) REPEAT('-',72)
        WRITE(7, *) REPEAT('-',72)
C
        ICNT = 1
        DO JCNT=2,NCNT
          ELA = ELMT(INT(ZNUC(ICNT)))
          ELB = ELMT(INT(ZNUC(JCNT)))
          R1X = BXYZ(1,JCNT) - BXYZ(1,ICNT)
          R1Y = BXYZ(2,JCNT) - BXYZ(2,ICNT)
          R1Z = BXYZ(3,JCNT) - BXYZ(3,ICNT)
          D1  = DSQRT(R1X*R1X + R1Y*R1Y + R1Z*R1Z)
          WRITE(6,32) ELA,ELB,CFCT*D1
          WRITE(7,32) ELA,ELB,CFCT*D1
C
          DO KCNT=2,JCNT-1
            ELA = ELMT(INT(ZNUC(ICNT)))
            ELB = ELMT(INT(ZNUC(JCNT)))
            ELC = ELMT(INT(ZNUC(KCNT)))
            R2X = BXYZ(1,KCNT) - BXYZ(1,ICNT)
            R2Y = BXYZ(2,KCNT) - BXYZ(2,ICNT)
            R2Z = BXYZ(3,KCNT) - BXYZ(3,ICNT)
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
            SEP  = DSQRT((BXYZ(1,ICNT) - BXYZ(1,JCNT))**2
     #                  +(BXYZ(2,ICNT) - BXYZ(2,JCNT))**2
     #                  +(BXYZ(3,ICNT) - BXYZ(3,JCNT))**2)
            EPNT = ZNUC(ICNT)*ZNUC(JCNT)/SEP
C
CC          TODO: MAKE SURE THIS IS RIGHT -- UPDATED WITH FERMI CHARGE
CC          THIS CODE INCLUDES GAUSSIAN CHARGE STRUCTURE EFFECTS,
CC          BUT CORRECTIONS EXCEED DOUBLE FLOAT ACCURACY LIMITS...
C           EERF = 0.0D0
C           DO IFT=1,NNUC(ICNT)
C             DO JFT=1,NNUC(JCNT)
C               EPRD = XNUC(ICNT,IFT)*XNUC(JCNT,JFT)
C               ESUM = XNUC(ICNT,IFT)+XNUC(JCNT,JFT)
C               EGAU = DSQRT(EPRD/ESUM)*SEP
C               EERF = EERF + FNUC(ICNT,IFT)*FNUC(JCNT,JFT)*DERF(EGAU)
C             ENDDO
C           ENDDO
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
C  FOCKIND ASSIGNS INDICES FOR FOCK MATRIX BLOCKS DEPENDING ON BASIS   C
C  NUMBERS IZ, KQN AND MQN QUANTUM NUMBERS OF EACH BASIS FUNCTION.     C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/QNMS/LABICN(MDM),LABKQN(MDM),LABMQN(MDM)
C
C     QUANTUM NUMBER LABELS
      ILST = 0
      DO IZ=1,NCNT
        DO IMV=0,MEL
          MQN = 2*IMV+1
C
C         LABEL MQN<0 BLOCKS
          DO KA=1,NKAP(IZ)
            KQN = KAPA(KA,IZ)
            LQN = LVAL(KQN)
            NBAS  = NFNC(LQN,IZ)
            MQMAX = 2*IABS(KQN)-1
            IF(MQMAX.GE.MQN) THEN
              LRGE(IZ,KA,MQN) = ILST
              DO IBAS=1,NBAS
                LABICN(ILST+IBAS) = IZ
                LABKQN(ILST+IBAS) = KQN
                LABMQN(ILST+IBAS) =-MQN
              ENDDO
              ILST = ILST+NBAS
            ENDIF
          ENDDO
C
C         LABEL MQN>0 BLOCKS
          DO KA=1,NKAP(IZ)
            KQN = KAPA(KA,IZ)
            LQN = LVAL(KQN)
            NBAS  = NFNC(LQN,IZ)
            MQMAX = 2*IABS(KQN)-1
            IF(MQMAX.GE.MQN) THEN
              LRGE(IZ,KA,MQN+1) = ILST
              DO IBAS=1,NBAS
                LABICN(ILST+IBAS) = IZ
                LABKQN(ILST+IBAS) = KQN
                LABMQN(ILST+IBAS) = MQN
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
      INCLUDE 'parameters.h'
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
C     LOOP OVER ALL (A,B,C) SO THAT A+B+C=LAM AND APPLY ORDERED ADDRESS
C     (USE T, U AND V AS INTERMEDIATE COUNTERS RESPECTIVELY)
      N = 0
      DO LAM=0,ML4
        DO IT=0,LAM
          DO IU=0,LAM
            DO IV=0,LAM
C
C             TEST WHETER A+B+C=LAM
              IF(IT+IU+IV.NE.LAM) GOTO 10
C
C             UPDATE ADDRESS
              N = N+1
C
C             CARTESIAN INDICES (A,B,C) AND LAM VALUE FOR THIS ADDRESS
              IA(N) = IT
              IB(N) = IU
              IC(N) = IV
              ILAM(N) = LAM
C
C             GLOBAL ADDRESS FOR THIS (A,B,C) INDEX
              IABC(IT,IU,IV) = N
C
C             SKIP POINT FOR (A,B,C) THAT DO NOT ADD TO LAM
10            CONTINUE
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE AUFBAU(IZNC,IQNC,NORB,NOCC,LMAX)
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
C  AUFBAU DETERMINES THE GROUND STATE ELECTRONIC CONFIGURATION OF AN   C
C  ATOM WITH IQNC ELECTRONS. LIMITS ON # ELECTRONS IN THE ATOM:        C
C -------------------------------------------------------------------- C
C  LMAX |  0    1    2    3    4    5 |   OR IN GENERAL FOR LMAX,      C
C  NMAX |  4   20   56  120  220  364 |   NMAX = (2L+2)(2L+3)(2L+4)/6. C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C  ▶ LMAX IS THE HIGHEST LQN REQUIRED TO DESCRIBE THE GROUND STATE     C
C  ▶ NOCC SAVES THE NUMBER OF OCCUPIED NSHELLS FOR THIS LQN CLASS      C
C  ▶ NORB SAVES THE NUMBER OF ELECTRONS IN OF TYPE LQN IN SHELL N      C
C  PARAMETERS:                                                         C
C  ▶ IORD STORES THE LQN OF ORBITALS IN ORDER OF HYDROGENIC ENERGY     C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION NORB(0:MEL,MKP+1),NOCC(0:MEL),IORD((MEL+1)*(MEL+2))
C
C     INITIALISE THE OCCUPIED NQN SHELL COUNTER FOR EACH LQN
      DO LQN=0,MEL
        NOCC(LQN) = 0
      ENDDO
C
C     RECORD LQNS AS THEY APPEAR IN A FULL AUFBAU COUNT UP TO LMAX
      ICT = 0
C
C     EACH LQN AS THE HIGHEST
      DO LHIGH=0,MEL
C
C       GET TWO DIAGONAL STRIKES FOR THIS LHIGH
        DO MDIAG=1,2
C
C         GO BACKWARDS FROM THIS LQN DOWN TO ZERO AND RECORD LQN
          DO LQN=LHIGH,0,-1
            ICT = ICT + 1
            IORD(ICT) = LQN
          ENDDO
C
        ENDDO
C
      ENDDO
C
C     INITIALISE THE MAX LQN COUNTER
      LMAX = 0
C
C     INITIALISE THE NUMBER OF ELECTRONS LEFT TO FILL ORBITALS WITH
      ILEFT = IQNC
C
C     INITIALISE LOOP OVER ORBITALS
      DO M=1,(MEL+1)*(MEL+2)
C
C       EXIT IF THERE ARE NO ELECTRONS LEFT TO COUNT
        IF(ILEFT.EQ.0) GOTO 20
C
C       THE LQN OF THIS ORBITAL IS STORED IN IORD
        LQN = IORD(M)
C
C       UPDATE THE MAX LQN COUNTER IF NECESSARY
        IF(LQN.GT.LMAX) THEN
          LMAX = LQN
        ENDIF
C
C       ADD TO THE NUMBER OF OCCUPIED NSHELLS FOR THIS LQN CLASS
        NOCC(LQN) = NOCC(LQN)+1
C
C       DETERMINE NO. OF ELECTRONS REQ'D TO FULLY OCCUPY THIS SUBSHELL
        IFULL = 4*LQN + 2
C
C       BEGIN IF STATEMENT TO DETERMINE THE SUBSHELL OCCUPATION
        IF(ILEFT.GT.IFULL) THEN
C
C         IF THERE ARE MORE ELECTRONS LEFT THAN IFULL, FILL THE SUBSHELL
          NORB(LQN,NOCC(LQN)) = IFULL
          ILEFT = ILEFT-IFULL
C
        ELSE
C
C         OTHERWISE, LEAVE ALL REMAINING ELECTRONS IN THIS NSHELL
          NORB(LQN,NOCC(LQN)) = ILEFT
          GOTO 20
C
C       END THE NSHELL IF STATEMENT
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
C  ROTATE PERFORMS TWO ROTATIONS ON ALL ATOMIC CENTRES, FIRST SO THAT  C
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
      FUNCTION MS(TSEC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                          MM       MM  SSSSSS                         C
C                          MMM     MMM SS    SS                        C
C                          MMMM   MMMM SS                              C
C                          MM MM MM MM  SSSSSS                         C
C                          MM  MMM  MM       SS                        C
C                          MM   M   MM SS    SS                        C
C                          MM       MM  SSSSSS                         C
C                                                                      C
C -------------------------------------------------------------------- C
C   MS RETURNS A QUOTED TIME IN SECONDS USING 'MIN-SEC' FORMAT.        C
C**********************************************************************C
      CHARACTER*4 MINUTES
      CHARACTER*7 SECONDS
      CHARACTER*11 MS
C
C     INITIALISE COUNTERS
      NMIN = 0
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
C     PRINT THE MINUTE
      IF(NMIN.EQ.0) THEN
        WRITE(MINUTES,'(A)') '    '
      ELSE
        WRITE(MINUTES,'(I3,A)') NMIN,'m'
      ENDIF
C
C     PRINT THE SECOND
      IF(NMIN.EQ.0) THEN
        IF(TSEC.LT.1.0D0) THEN
          WRITE(SECONDS,'(2X,F4.2,A)') TSEC,'s'
        ELSE
          WRITE(SECONDS,'(1X,F5.2,A)') TSEC,'s'
        ENDIF
      ELSEIF(NMIN.GT.0) THEN
        IF(TSEC.LT.10.0D0) THEN
          WRITE(SECONDS,'(1X,I1,F4.2,A)') 0,TSEC,'s'
        ELSE
          WRITE(SECONDS,'(1X,F5.2,A)') TSEC,'s'
        ENDIF
      ENDIF
C
      WRITE(MS,'(A,A)') MINUTES//SECONDS
C
C     RESTORE TSEC
      TSEC = TEMP
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
C   HMS RETURNS A QUOTED TIME IN SECONDS USING 'HR-MIN-SEC' FORMAT.    C
C**********************************************************************C
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
C  TTTTTTTT IIII MM       MM EEEEEEEE NN    NN  OOOOOO  WW        WW   C
C     TT     II  MMM     MMM EE       NNN   NN OO    OO WW        WW   C
C     TT     II  MMMM   MMMM EE       NNNN  NN OO    OO WW   WW   WW   C
C     TT     II  MM MM MM MM EEEEEE   NN NN NN OO    OO WW  WWWW  WW   C
C     TT     II  MM  MMM  MM EE       NN  NNNN OO    OO WW WW  WW WW   C
C     TT     II  MM   M   MM EE       NN   NNN OO    OO WWWW    WWWW   C
C     TT    IIII MM       MM EEEEEEEE NN    NN  OOOOOO  WW        WW   C
C                                                                      C
C -------------------------------------------------------------------- C
C  TIMENOW CREATES A DATE STRING SPECIFYING CPU TIME WHEN CALLED.      C
C**********************************************************************C
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
      SUBROUTINE PROGRSS(I,N,TBEG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    PPPPPPP  RRRRRRR   OOOOOO   GGGGGG  RRRRRRR   SSSSSS   SSSSSS     C
C    PP    PP RR    RR OO    OO GG    GG RR    RR SS    SS SS    SS    C
C    PP    PP RR    RR OO    OO GG       RR    RR SS       SS          C
C    PP    PP RR    RR OO    OO GG       RR    RR  SSSSSS   SSSSSS     C
C    PPPPPPP  RRRRRRR  OO    OO GG   GGG RRRRRRR        SS       SS    C
C    PP       RR    RR OO    OO GG    GG RR    RR SS    SS SS    SS    C
C    PP       RR    RR  OOOOOO   GGGGGG  RR    RR  SSSSSS   SSSSSS     C
C                                                                      C
C -------------------------------------------------------------------- C
C  PROGRSS TRACKS HOW MUCH WORK HAS BEEN DONE AND EXPECTED COMPLETION. C
C  REQUIRES NUMBER OF PROCESSES HAVE BEEN FINISHED, HOW MANY ARE TO    C
C  BE DONE IN TOTAL, AND START TIME OF JOB AS A WHOLE.                 C
C**********************************************************************C
C
      CHARACTER*16 HMS
C
      INTEGER*16 I,N
C
C     COMPLETION RATIO
      RAT = DFLOAT(I)/DFLOAT(N)
C
C     TIME NOW
      CALL CPU_TIME(TNOW)
C
C     DIFFERENCE BETWEEN TIME NOW AND START TIME
      TDIF = TNOW-TBEG
C
C     PROJECTED TOTAL COMPLETION TIME
      TTOT = (TNOW-TBEG)/RAT
C
C     PROJECTED END TIME
      TEND = TBEG + TTOT
C
C     EXPECTED TIME LEFT
      TLFT = TTOT-TNOW+TBEG
C
C     FOR COMPLETION MARKERS
10    FORMAT('Completion: ',F8.4,' %',8X,'Time left: ',1X,A)
      WRITE(6,10) 100.0D0*RAT,HMS(TLFT)
      WRITE(7,10) 100.0D0*RAT,HMS(TLFT)
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
C**********************************************************************C
C
C
      SUBROUTINE EIGSORT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       EEEEEEEE IIII GGGGGG   SSSSSS   OOOOOO  RRRRRRR TTTTTTTT       C
C       EE        II GG    GG SS    SS OO    OO RR    RR   TT          C
C       EE        II GG       SS       OO    OO RR    RR   TT          C
C       EEEEEE    II GG        SSSSSS  OO    OO RR    RR   TT          C
C       EE        II GG   GGG       SS OO    OO RRRRRRR    TT          C
C       EE        II GG    GG SS    SS OO    OO RR    RR   TT          C
C       EEEEEEEE IIII GGGGGG   SSSSSS   OOOOOO  RR    RR   TT          C
C                                                                      C
C -------------------------------------------------------------------- C
C  EIGSORT RE-ORGANISES THE ENERGY EIGENVALUES AND COEFFCIENT ARRAY    C
C  IN PREPARATION FOR MOLECULAR CALCULATIONS (IOCC: 1->NOCC).          C
C  CONTEXT: ONLY ACTUALLY NEED THIS WHEN DEALING WITH SELF-ENERGY.     C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMPLEX*16 CTMP(MDM)
      COMPLEX*16 COEF(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
C
C     LOOP OVER OCCUPIED LEVELS AS INDICATED IN IOCPN
      DO IOCC=1,IOCCM0
        IOCAD = IOCPN(IOCC)
C
C       SKIP IF IOCC AND IOCAD ARE THE SAME
        IF(IOCC.EQ.IOCAD) GOTO 101
C
C       EXECUTE A SWAP
        ETMP = EIGN(IOCAD+NSKP)
        DO I=1,NDIM
          CTMP(I) = COEF(I,IOCAD+NSKP)
        ENDDO
C
        EIGN(IOCAD+NSKP) = EIGN(IOCC+NSKP)
        DO I=1,NDIM
          COEF(I,IOCAD+NSKP) = COEF(I,IOCC+NSKP)
        ENDDO
        
        EIGN(IOCC+NSKP) = ETMP
        DO I=1,NDIM
          COEF(I,IOCC+NSKP) = CTMP(I)
        ENDDO
C
101     CONTINUE
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE DENSTY0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    DDDDDDD  EEEEEEEE NN    NN  SSSSSS TTTTTTTT YY    YY  000000      C
C    DD    DD EE       NNN   NN SS    SS   TT    YY    YY 00   000     C
C    DD    DD EE       NNNN  NN SS         TT     YY  YY  00  0000     C
C    DD    DD EEEEEE   NN NN NN  SSSSSS    TT      YYYY   00 00 00     C
C    DD    DD EE       NN  NNNN       SS   TT       YY    0000  00     C
C    DD    DD EE       NN   NNN SS    SS   TT       YY    000   00     C
C    DDDDDDD  EEEEEEEE NN    NN  SSSSSS    TT       YY     000000      C
C                                                                      C
C -------------------------------------------------------------------- C
C  DENSTY0 IS A STARTING DENSITY ROUTINE FOR USE ONLY WHEN READIN =.F. C
C  THIS IS BECAUSE 'ATOMIC' CALCULATES AVERAGE OVER SHELL OCCUPANCIES. C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMPLEX*16 SUM
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/DENS/DENC,DENO,DENT
      COMMON/EIGC/COEF
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
      COMMON/SHLL/ACFF,BCFF,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
C
C     CONSTRUCT THE CLOSED-SHELL AND TOTAL DENSITY ARRAYS
      DO I=1,NDIM
        DO J=1,NDIM
          SUM = DCMPLX(0.0D0,0.0D0)
          DO IOCC=1,IOCCM0
            IOCAD = IOCPN(IOCC)
            SUM = SUM + DCONJG(COEF(I,IOCAD+NSKP))*COEF(J,IOCAD+NSKP)
          ENDDO
          DENC(I,J) = SUM
          DENT(I,J) = DENC(I,J)
          DENO(I,J) = DCMPLX(0.0D0,0.0D0)
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
C  DENSTY GENERATES DENSITY MATRICES FROM THE EXPANSION COEFFS COEF.   C
C**********************************************************************C
      INCLUDE 'parameters.h'
c
      COMPLEX*16 SUM
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/DENS/DENC,DENO,DENT
      COMMON/EIGC/COEF
      COMMON/SHLL/ACFF,BCFF,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
C
C     MAKE CLOSED-SHELL DENSITY AND EMPTY OPEN-SHELL DENSITY (RSCF 81)
      DO I=1,NDIM
        DO J=1,NDIM
          SUM = DCMPLX(0.0D0,0.0D0)
          DO IOCC=1,NCLS
            ICL = ICLS(IOCC)
            SUM = SUM +      DCONJG(COEF(I,ICL+NSKP))*COEF(J,ICL+NSKP)
          ENDDO
          DENC(I,J) = SUM
          DENT(I,J) = DENC(I,J)
          DENO(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
      IF(NOPN.EQ.0) GOTO 100
C
C     MAKE THE OPEN-SHELL DENSITY (RSCF 82)
      DO I=1,NDIM
        DO J=1,NDIM
          SUM = DCMPLX(0.0D0,0.0D0)
          DO IOCC=1,NOPN
            IOP = IOPN(IOCC)
            SUM = SUM + FOPN*DCONJG(COEF(I,IOP+NSKP))*COEF(J,IOP+NSKP)
          ENDDO
          DENO(I,J) = SUM
          DENT(I,J) = DENT(I,J)+DENO(I,J)
        ENDDO
      ENDDO
C
100   CONTINUE
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
      INCLUDE 'parameters.h'
C
      COMPLEX*16 A(MDM),SUM
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VANM(MDM,MDM),
     &           VSLF(MDM,MDM),VUEH(MDM,MDM),VWKR(MDM,MDM),
     &           VKSB(MDM,MDM),QDIR(MDM,MDM),QXCH(MDM,MDM),
     &           WDIR(MDM,MDM),WXCH(MDM,MDM),CPLE(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/EIGC/COEF
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VANM,VSLF,
     &            VUEH,VWKR,VKSB,QDIR,QXCH,WDIR,WXCH,CPLE
C
C     LOOP OVER ALL NON-OCCUPIED POSITIVE-ENERGY (I.E., VIRTUAL) OCCS
      DO IVIR=NSKP+NOCC+1,NDIM
C
C       LIST OF FOCK ADDRESSES WITH 'I' WITH AN AMOUNT FROM EACH IVIR
        DO I=1,NDIM
          SUM = DCMPLX(0.0D0,0.0D0)
          DO J=1,NDIM
            SUM = SUM + OVLP(I,J)*COEF(J,IVIR)
          ENDDO
          A(I) = SUM
        ENDDO
C
C       ADD A LEVEL SHIFT'S WORTH OF INTERACTION ENERGY TO ALL FOCK(I,J)
        DO I=1,NDIM
          DO J=1,NDIM
            FOCK(I,J) = FOCK(I,J) + SHLV*A(I)*DCONJG(A(J))
          ENDDO
        ENDDO
C
      ENDDO
C
      RETURN
      END
C
C**********************************************************************C
C ==================================================================== C
C   [4] ATOMIC HARTREE-FOCK: SINGLE-CENTRE SCF CALCULATIONS.           C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [N] GAMLWR: LOWER INCOMPLETE GAMMA FUNCTION gamma(A,X).            C
C   [O] GAMUPR: UPPER INCOMPLETE GAMMA FUNCTION GAMMA(A,X).            C
C**********************************************************************C
C
C
      FUNCTION ERFINT0(L,EIJ,ZTA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      EEEEEEEE RRRRRRR  FFFFFFFF IIII NN    NN TTTTTTTT 000000        C
C      EE       RR    RR FF        II  NNN   NN    TT   00   000       C
C      EE       RR    RR FF        II  NNNN  NN    TT   00  0000       C
C      EEEEEE   RR    RR FFFFFF    II  NN NN NN    TT   00 00 00       C
C      EE       RRRRRRR  FF        II  NN  NNNN    TT   0000  00       C
C      EE       RR    RR FF        II  NN   NNN    TT   000   00       C
C      EEEEEEEE RR    RR FF       IIII NN    NN    TT    000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  ERFINT0 CALCULATES THE VALUE OF AN INTEGRAL OVER A GAUSSIAN WITH    C
C  EXPONENT EIJ, AN ERROR FUNCTION WITH NUCLEAR WIDTH PARAMETER ZTA    C
C  AND AN ODD POLYNOMIAL ORDER, 2L+1. SUPPORTS UP TO RELATIVISTIC L=9. C
C                           ∞                                          C
C          ERFINT(L,λ,ζ) = ∫ r^2L+1 exp(-λ r^2) erf(√ζ r) dr.          C
C                           0                                          C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In ERFINT0: order L must be positive. L = ',L
        WRITE(7, *) 'In ERFINT0: order L must be positive. L = ',L
        STOP
      ENDIF
C
C     FACTORS NEEDED FOR ALL PARAMETERS L
      X   = EIJ/ZTA
      X5  = X*X*X*X*X
      T0  = ZTA+EIJ
      RAT = ZTA/T0
      TRM = 0.5D0*DSQRT(ZTA)/EIJ/DSQRT(T0)
      DO I=1,L
        TRM = 0.5D0*TRM*RAT/EIJ
      ENDDO
C
      IF(L.EQ.0) THEN
        TRM = TRM
      ELSEIF(L.EQ.1) THEN
        VA  = 2.0D0 + 3.0D0*X
        TRM = TRM*VA
      ELSEIF(L.EQ.2) THEN
        VA  = 8.0D0 + 20.0D0*X + 15.0D0*X*X
        TRM = TRM*VA
      ELSEIF(L.EQ.3) THEN
        VA  = 16.0D0 + 56.0D0*X + 70.0D0*X*X + 35.0D0*X*X*X
        TRM = 3.0D0*TRM*VA
      ELSEIF(L.EQ.4) THEN
        VA  = 128.0D0 + 576.0D0*X + 1008.0D0*X*X + 840.0D0*X*X*X
        VB  = 315.0D0*X*X*X*X
        TRM = 3.0D0*TRM*(VA+VB)
      ELSEIF(L.EQ.5) THEN
        VA  = 256.0D0 + 1408.0D0*X + 3168.0D0*X*X + 3696.0D0*X*X*X
        VB  = 2310.0D0*X*X*X*X + 693.0D0*X*X*X*X*X
        TRM = 15.0D0*TRM*(VA+VB)
      ELSEIF(L.EQ.6) THEN
        VA  = 1024.0D0 + 6656.0D0*X + 18304.0D0*X*X
        VB  = 27456.0D0*X*X*X + 24024.0D0*X*X*X*X
        VC  = 12012.0D0*X5 + 3003.0D0*X5*X
        TRM = 45.0D0*TRM*(VA+VB+VC)
      ELSEIF(L.EQ.7) THEN
        VA  = 2048.0D0 + 15360.0D0*X + 49920.0D0*X*X
        VB  = 91520.0D0*X*X*X+ 102960.0D0*X*X*X*X + 72072.0D0*X5
        VC  = 30030.0D0*X5*X + 6435.0D0*X5*X*X
        TRM = 315.0D0*TRM*(VA+VB+VC)
      ELSEIF(L.EQ.8) THEN
        VA  = 32768.0D0 + 278528.0D0*X + 1044480.0D0*X*X
        VB  = 2263040.0D0*X*X*X + 3111680.0D0*X*X*X*X
        VC  = 2800512.0D0*X5 + 1633632.0D0*X5*X + 583440.0D0*X5*X*X
        VD  = 109395.0D0*X5*X*X*X
        TRM = 315.0D0*TRM*(VA+VB+VC+VD)
      ELSEIF(L.EQ.9) THEN
        VA  = 65536.0D0 + 6222592.0D0*X + 2646016.0D0*X*X
        VB  = 6615040.0D0*X*X*X + 10749440.0D0*X*X*X*X
        VC  = 11824384.0D0*X5 + 8868288.0D0*X5*X + 4434144.0D0*X5*X*X
        VD  = 1385670.0D0*X5*X*X*X + 230945.0D0*X5*X*X*X*X
        TRM = 2835.0D0*TRM*(VA+VB+VC+VD)
      ELSEIF(L.EQ.10) THEN
        VA  = 262144.0D0 + 2752512.0D0*X + 13074432.0D0*X*X
        VB  = 37044224.0D0*X*X*X + 69457920.0D0*X*X*X*X
        VC  = 90295296.0D0*X5 + 82770688.0D0*X5*X
        VD  = 53209728.0D0*X5*X*X + 23279256.0D0*X5*X*X*X
        VE  = 6466460.0D0*X5*X*X*X*X + 969969.0D0*X5*X5
        TRM = 14175.0D0*TRM*(VA+VB+VC+VD+VE)
      ELSE
        WRITE(6, *) 'In ERFINT0: order L too large. L = ',L
        WRITE(7, *) 'In ERFINT0: order L too large. L = ',L
      ENDIF
C
C     TRANSFER DATA TO ERFINT0
      ERFINT0 = TRM
C
      RETURN
      END
C
C
      FUNCTION GAMUPR(L,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        GGGGGG     AA    MM       MM UU    UU PPPPPPP  RRRRRRR        C
C       GG    GG   AAAA   MMM     MMM UU    UU PP    PP RR    RR       C
C       GG        AA  AA  MMMM   MMMM UU    UU PP    PP RR    RR       C
C       GG       AA    AA MM MM MM MM UU    UU PP    PP RR    RR       C
C       GG   GGG AAAAAAAA MM  MMM  MM UU    UU PPPPPPP  RRRRRRR        C
C       GG    GG AA    AA MM   M   MM UU    UU PP       RR    RR       C
C        GGGGGG  AA    AA MM       MM  UUUUUU  PP       RR    RR       C
C                                                                      C
C -------------------------------------------------------------------- C
C  GAMUPR RETURNS THE UPPER INCOMPLETE GAMMA FUNCTION FOR INTEGER OR   C
C  HALF-INTEGER ARGUMENTS. INPUT L IS TWICE THE ACTUAL ARGUMENT.       C
C  SOLUTIONS ARE FINITE AND ALGEBRAIC -- NO APPROXIMATIONS NEEDED.     C
C  NECESSARY -- AN ALGEBRAIC SOLUTION HAS BEEN DEDUCED.                C
C                                                                      C
C                GAMUPR(ℓ,x) = Γ(ℓ/2,x) with ℓ>0 and x>0.              C
C                                                                      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION PIM(L)
C
      COMMON/FCTS/RFACT(0:80),SFACT(0:80)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     ROUTINE ONLY ALLOWS PARAMETERS L>0
      IF(L.LT.1) THEN
        WRITE(6, *) 'In GAMUPR: order L must be at least 1. L = ',L
        WRITE(7, *) 'In GAMUPR: order L must be at least 1. L = ',L
        STOP
      ENDIF
C
C     INTEGER ARGUMENTS L
      IF(MOD(L,2).EQ.0) THEN
C
C       EXPONENTIAL FACTOR
        FC1 = DEXP(-X)
C
C       INITIALISE POLYNOMIAL COUNTERS
        POLY = 0.0D0
        XPOW = 1.0D0
C
C       LOOP OVER REQUIRED POLYNOMIAL DEGREES
        DO I=0,(L-2)/2
          POLY = POLY + XPOW*RFACT((L-2)/2)/RFACT(I)
          XPOW = XPOW*X
        ENDDO
C
C       VALUE OF INTEGER-ORDER UPPER INCOMPLETE GAMMA FUNCTION
        GAMUPR = FC1*POLY
C
      ELSE
C
C       FACTORS REQUIRED FOR ALL ORDERS
        X12 = DSQRT(X)
        FC1 = PI12*(1.0D0-DERF(X12))
        FC2 = DEXP(-X)*X12
C
C       POCHHAMMER SYMBOL
        P32 = 1.0D0
        DO K=0,(L-3)/2
          P32 = P32*(0.5D0+DFLOAT(K))
        ENDDO
C
C       FACTOR FOR EACH POWER TERM
        LST = (L-1)/2
        IF(LST.GT.0) THEN
          PSD = 1.0D0
          DO I=1,LST
            PIM(I) = PSD
            PSD    = 0.5D0*PSD*DFLOAT(L-2*I)
          ENDDO
        ENDIF
C
C       INITIALISE POLYNOMIAL COUNTERS
        POLY = 0.0D0
        XPOW = 1.0D0
C
C       LOOP OVER REQUIRED POLYNOMIAL DEGREES
        DO I=0,(L-3)/2
C
          POLY = POLY + XPOW*PIM(LST-I)
          XPOW = XPOW*X
C
        ENDDO
C
C       VALUE OF INTEGER-ORDER UPPER INCOMPLETE GAMMA FUNCTION
        GAMUPR = P32*FC1 + POLY*FC2
C
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION GAMLWR(L,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      GGGGGG     AA    MM       MM LL       WW        WW RRRRRRR      C
C     GG    GG   AAAA   MMM     MMM LL       WW        WW RR    RR     C
C     GG        AA  AA  MMMM   MMMM LL       WW   WW   WW RR    RR     C
C     GG       AA    AA MM MM MM MM LL       WW  WWWW  WW RR    RR     C
C     GG   GGG AAAAAAAA MM  MMM  MM LL       WW WW  WW WW RRRRRRR      C
C     GG    GG AA    AA MM   M   MM LL       WWWW    WWWW RR    RR     C
C      GGGGGG  AA    AA MM       MM LLLLLLLL WW        WW RR    RR     C
C                                                                      C
C -------------------------------------------------------------------- C
C  GAMLWR RETURNS THE LOWER INCOMPLETE GAMMA FUNCTION FOR INTEGER OR   C
C  HALF-INTEGER ARGUMENTS. INPUT L IS TWICE THE ACTUAL ARGUMENT.       C
C  SOLUTIONS ARE FINITE AND ALGEBRAIC -- NO APPROXIMATIONS NEEDED.     C
C                                                                      C
C                GAMLWR(ℓ,x) = γ(ℓ/2,x) with ℓ>0 and x>0.              C
C                                                                      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION PIM(L)
C
      COMMON/FCTS/RFACT(0:80),SFACT(0:80)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     ROUTINE ONLY ALLOWS PARAMETERS L>0
      IF(L.LT.1) THEN
        WRITE(6, *) 'In GAMLWR: order L must be at least 1. L = ',L
        WRITE(7, *) 'In GAMLWR: order L must be at least 1. L = ',L
        STOP
      ENDIF
C
C     INTEGER ARGUMENTS L
      IF(MOD(L,2).EQ.0) THEN
C
C       EXPONENTIAL FACTOR
        FC1 = DEXP(-X)
C
C       INITIALISE POLYNOMIAL COUNTERS
        POLY = 0.0D0
        XPOW = 1.0D0
C
C       LOOP OVER REQUIRED POLYNOMIAL DEGREES
        DO I=0,(L-2)/2
          POLY = POLY + XPOW*RFACT((L-2)/2)/RFACT(I)
          XPOW = XPOW*X
        ENDDO
C
C       VALUE OF INTEGER-ORDER UPPER INCOMPLETE GAMMA FUNCTION
        GAMLWR = RFACT((L-2)/2) - FC1*POLY
C
      ELSE
C
C       FACTORS REQUIRED FOR ALL ORDERS
        X12 = DSQRT(X)
        FC1 = PI12*DERF(X12)
        FC2 = DEXP(-X)*X12
C
C       POCHHAMMER SYMBOL
        P32 = 1.0D0
        DO K=0,(L-3)/2
          P32 = P32*(0.5D0+DFLOAT(K))
        ENDDO
C
C       FACTOR FOR EACH POWER TERM
        LST = (L-1)/2
        IF(LST.GT.0) THEN
          PSD = 1.0D0
          DO I=1,LST
            PIM(I) = PSD
            PSD    = 0.5D0*PSD*DFLOAT(L-2*I)
          ENDDO
        ENDIF
C
C       INITIALISE POLYNOMIAL COUNTERS
        POLY = 0.0D0
        XPOW = 1.0D0
C
C       LOOP OVER REQUIRED POLYNOMIAL DEGREES
        DO I=0,(L-3)/2
C
          POLY = POLY + XPOW*PIM(LST-I)
          XPOW = XPOW*X
C
        ENDDO
C
C       VALUE OF INTEGER-ORDER UPPER INCOMPLETE GAMMA FUNCTION
        GAMLWR = P32*FC1 - POLY*FC2
C
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION FMIINT0(L,EIJ,A,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       FFFFFFFF MM       MM IIII IIII NN    NN TTTTTTTT 000000        C
C       FF       MMM     MMM  II   II  NNN   NN    TT   00   000       C
C       FF       MMMM   MMMM  II   II  NNNN  NN    TT   00  0000       C
C       FFFFFF   MM MM MM MM  II   II  NN NN NN    TT   00 00 00       C
C       FF       MM  MMM  MM  II   II  NN  NNNN    TT   0000  00       C
C       FF       MM   M   MM  II   II  NN   NNN    TT   000   00       C
C       FF       MM       MM IIII IIII NN    NN    TT    000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  FMIINT0 CALCULATES AN AUXILLIARY FERMI POTENTIAL INTEGRAL OVER      C
C  A PAIR OF ATOMIC BASIS FUNCTIONS, WHERE L IS A POSITIVE INTEGER.    C
C                        ∞                                             C
C      FMIINT0(L,λ,ζ) = ∫ r^2L+1 exp(-λ r^2) r*V_fermi(r) dr.          C
C                        0                                             C
C  THIS ROUTINE USES MY OWN RECIPE FOR THE ANNOYING INTEGRAL CLASSES,  C
C  FOR THE FIRST FEW L CASES AND ON THE CONDITION THAT A√λ < 0.05D0.   C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/FCTS/RFACT(0:80),SFACT(0:80)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
      DIMENSION U(0:60),P(24)
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In FMIINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In FMIINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
C     MULTIPLES OF THE FACTOR U
      U(0) = 1.0D0
      DO K=1,60
        U(K) = U(K-1)*A/C
      ENDDO
C
C     SPECIAL VALUES OF THE RIEMANN-ZETA FUNCTION
      P( 4) =-(PI**4)*7.0D0/720.0D0
      P( 6) =-(PI**6)*31.0D0/30240.0D0
      P( 8) =-(PI**8)*127.0D0/12096.0D2
      P(10) =-(PI**10)*73.0D0/684288.0D1
      P(12) =-(PI**12)*1414477.0D0/1307674368.0D3
      P(14) =-(PI**14)*8191.0D0/747242496.0D2
      P(16) =-(PI**16)*16931177.0D0/152437569184.0D4
      P(18) =-(PI**18)*5749691557.0D0/5109094217170944.0D3
      P(20) =-(PI**20)*91546277357.0D0/8028576626982912.0D5
      P(22) =-(PI**22)*3324754717.0D0/28777755182432256.0D4
      P(24) =-(PI**24)*1982765468311237.0D0/16938241367317436694528.0D5
C
C     CORRECT DIMENSIONS ARISE FROM THESE TERMS
      EL1 = 1.0D0/(EIJ**(L+1))
      C2L = C**(2*L+2)
C
C     DIMENSIONLESS GAUSSIAN PARAMETER
      SIGMA = EIJ*C*C
C
C     FERMI NORMALISATION CONSTANT
      FNRM = 1.0D0 + PI*PI*U(2) - 6.0D0*U(3)*POLYLOG(3,-1.0D0/U(1))
      FNRM = 1.0D0/FNRM
C
C     INTEGRAL TYPE: INCOMPLETE GAMMA FUNCTIONS
      FA  = 1.5D0 + PI*PI*U(2)*0.5D0
      FB  =-0.5D0
      FC  = 1.0D0 + PI*PI*U(2)
      Y1A = 0.5D0*GAMLWR(2*L+3,SIGMA)/(DSQRT(SIGMA))
      Y1B = 0.5D0*GAMLWR(2*L+5,SIGMA)/(SIGMA*DSQRT(SIGMA))
      Y1C = 0.5D0*GAMUPR(2*L+2,SIGMA)
      Y1  = FA*Y1A + FB*Y1B + FC*Y1C
C
C     INTEGRAL TYPE: SEMI-INFINITE GAMMA FUNCTIONS
      FD  =-6.0D0*U(3)*POLYLOG(3,-1.0D0/U(1))
      Y2D = 0.5D0*RFACT(L)
      
      Y2  = FD*Y2D
C
C     INTEGRAL TYPE: POLYLOG FUNCTIONS WITH GAUSSIANS AND POLYNOMIALS
C                    (EVALUATION METHOD DEPENDS ON SIGMA PARAMETER)
      Y34 = 0.0D0
C
C     EXPAND THE GAUSSIAN FACTOR IN TAYLOR SERIES AND ASSUME LARGE N/U
C      IF(SIGMA.LT.1.0D0) THEN
       GOTO 123
C
C       MAXIMUM EXPANSION ORDER FOR GAUSSIAN
        KMAX = 20
C
C       INITIALISE BINS FOR GAUSSIAN EXPANSION SERIES     
        PHK = 1.0D0
        SGK = 1.0D0
        FTK = 1.0D0
C
C       TAYLOR EXPANSION OVER THE GAUSSIAN FOR SMALL SIGMA
        DO K=0,KMAX
C
          BIN = U(1)*U(2*L+2*K)*POLYLOG(2*L+2*K+5,-1.0D0/U(1))
          DO J=0,MIN(L+K,10)
C
            PR = 1.0D0/RFACT(2*L+2*K-2*J+1)
            BIN = BIN + 2.0D0*U(2*J)*P(2*J+4)*PR
C
          ENDDO
C             
          Y34 = Y34 + PHK*SGK*DFLOAT(L+K+2)*RFACT(2*L+2*K+1)*BIN/FTK
C
C         UPDATE PHASES, FACTORIALS AND POLYNOMIALS IN K
          PHK =-PHK
          SGK = SGK*SIGMA
          FTK = FTK*DFLOAT(K+1)
C          
        ENDDO
C     
C       INCLUDE MULTIPLICATIVE FACTOR
        Y34 = 6.0D0*U(4)*Y34
C
C     NUMERICAL QUADRATURE (WORKS FOR ANY CASE BUT COMPUTATIONALLY SLOW)
C     ELSE
123     CONTINUE
C
C       INTEGRATION GRID DETAILS
        NLIN = 100
        NEXP = 1000
        RMAX = 15.0D0
        ELIN = 1.0D0/DFLOAT(NLIN)
        EEXP = DLOG(RMAX/C)/DFLOAT(NEXP)
C
C       NUMERICALLY INTEGRATE REMAINING TERMS FROM 0 TO C (LINEAR)
        X0C = 0.0D0
        DO N=0,NLIN
          ZN  = ELIN*DFLOAT(N)
          Z1  = ZN**(2*L+1)
          Z2  = DEXP(-SIGMA*ZN*ZN)
          Z3  =-ZN*POLYLOG(2,(ZN-1.0D0)/U(1))
          Z4  = 2.0D0*U(1)*POLYLOG(3,(ZN-1.0D0)/U(1))
          X0C = X0C + EXTINT11(Z1*Z2*(Z3+Z4),N,NLIN)
        ENDDO
        X0C = 15.0D0*U(2)*ELIN*X0C/299376.0D0
C
C       NUMERICALLY INTEGRATE REMAINING TERMS FROM C TO INFINITY (EXP.)
        XCI = 0.0D0
        DO N=0,NEXP
          TN  = EEXP*DFLOAT(N)
          ZN  = DEXP(TN)
          Z1  = ZN**(2*L+2)
          Z2  = DEXP(-SIGMA*ZN*ZN)
          Z3  = ZN*POLYLOG(2,(1.0D0-ZN)/U(1))
          Z4  = 2.0D0*U(1)*POLYLOG(3,(1.0D0-ZN)/U(1))
          XCI = XCI + EXTINT11(Z1*Z2*(Z3+Z4),N,NEXP)
        ENDDO
        XCI = 15.0D0*U(2)*EEXP*XCI/299376.0D0
C        
        Y34 = X0C+XCI
C
C      ENDIF
C
C     COMBINE ALL TERMS AND APPLY NORMALISATION FACTOR
      FMIINT0 = FNRM*(EL1*Y1 + EL1*Y2 + C2L*Y34)
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C   [5] MOLECULAR HARTREE-FOCK: MANY-CENTRE SCF CALCULATIONS.          C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [B] OVRLP: CONSTRUCTS THE ONE-ELECTRON OVERLAP MATRIX.             C
C   [C] ONEEL: ONE-ELECTRON MULTI-CENTRE MATRIX OF INTEGRALS.          C
C   [F] COULOMB: CONSTRUCTS ALL MULTI-CENTRE COULOMB SCF ELEMENTS.     C
C   [G] ERI: GENERATES A BLOCK OF ELECTRON REPULSION INTEGRALS.        C
C   [I] BREIT: MATRIX REP OF MEAN-FIELD BREIT INTERACTION.             C
C   [J] BII: GENERATES A BLOCK OF BREIT INTERACTION INTEGRALS.         C
C   [N] NCART: RETURNS THE CARTESIAN INDEX FROM A LOOP INDEX.          C
C   [O] NCNTRS: RETURNS NUMBER OF UNIQUE NUCLEAR CENTRES FROM INPUT.   C
C   [P] ZPROJ: RETURNS A METRIC FOR Z-AXIS DEVIATION OF SOME COORDS.   C
C   [Q] KLAB: RETURNS THE CONVENTIONAL DIRAC SYMMETRY TYPE STRING.     C
C   [R] MLAB: RETURNS THE CONVENTIONAL DIRAC MAGNETIC STATE STRING.    C
C**********************************************************************C
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
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      DIMENSION RC(MB2,MRC),EXL(MBS,4),XYZ(3,4)
      DIMENSION KQN(4),LQN(4),MQN(4),NBAS(4)
C
      COMPLEX*16 CONE
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 SLL(MBS,MBS,4),SSS(MBS,MBS,4)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VANM(MDM,MDM),
     &           VSLF(MDM,MDM),VUEH(MDM,MDM),VWKR(MDM,MDM),
     &           VKSB(MDM,MDM),QDIR(MDM,MDM),QXCH(MDM,MDM),
     &           WDIR(MDM,MDM),WXCH(MDM,MDM),CPLE(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/E0LL/E0LLFL(MFL,4),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,4),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VANM,VSLF,
     &            VUEH,VWKR,VKSB,QDIR,QXCH,WDIR,WXCH,CPLE
      COMMON/TMMD/TELL,TESS,TELS,TESL,TRLL,TRSS,TRLS,TRSL,TRBR
      COMMON/TSCF/TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRW,TCC1,TCC2,TCMC,
     &            TB1B,TB1R,TB1F,TB1M,TBEC,TBRM,TBRW,TBC1,TBC2,TBMC,
     &            TSMX,TUMX,THMX,TAMX,TC1T,TC2T,TCVT,TB1T,TB2T,TACC,
     &            TEIG,TSCR,TTOT,TC1S,TC2S,TB1S,TB2S
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     INITIALISE THE OVERLAP ARRAY
      DO I=1,NDIM
        DO J=1,NDIM
          OVLP(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)      C
C**********************************************************************C
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
        LQN(1) = LVAL(KQN(1))
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
        LQN(2) = LVAL(KQN(2))
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
C     CALCULATE COMPONENT OFFSETS
      IL1 = LRGE(ICNTA,KA,MJA  )
      IL2 = LRGE(ICNTA,KA,MJA+1)
      JL1 = LRGE(ICNTB,KB,MJB  )
      JL2 = LRGE(ICNTB,KB,MJB+1)
C
      IS1 = IL1+NSKP
      IS2 = IL2+NSKP
      JS1 = JL1+NSKP
      JS2 = JL2+NSKP
C
C     CONSTRUCTION OF ONE-ELECTRON OVERLAP MATRIX BY TT' BLOCKS...
C
C     THIS PHASE RELATES EQ22 AND EQ12 COEFFS TO EQ11 AND EQ21
      PHS = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
     &                    *DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NBAS(1)*NBAS(2)
C
C**********************************************************************C
C     PART 1: THE LL MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)
      NTUVLL = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ELL0 COEFFICIENTS (IPHS=+1)
      CALL CPU_TIME(TDM1)
      IF(EQFILE) THEN
        DO ITUV=1,NTUVLL
          IAD = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXAB
          DO M=1,MAXAB
            E11(M,ITUV) = DCMPLX(E0LLFL(IAD+M,1),E0LLFL(IAD+M,2))
            E21(M,ITUV) = DCMPLX(E0LLFL(IAD+M,3),E0LLFL(IAD+M,4))
          ENDDO
        ENDDO
      ELSE
        CALL EQLLMK(E11,E21,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,0)
      ENDIF
      CALL CPU_TIME(TDM2)
      TELL = TELL+TDM2-TDM1
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
          EIJ   = EXL(IBAS,1)+EXL(JBAS,2)
          EROOT = DSQRT(PI/EIJ)**3
          SLL(IBAS,JBAS,1) = EROOT*E11(M,1)
          SLL(IBAS,JBAS,3) = EROOT*E21(M,1)
          SLL(IBAS,JBAS,2) =-PHS*DCONJG(SLL(IBAS,JBAS,3))
          SLL(IBAS,JBAS,4) = PHS*DCONJG(SLL(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     NON-RELATIVISTIC OVERLAP CALCULATIONS COMPLETE
      IF(HMLT.EQ.'NORL') GOTO 500
C
C**********************************************************************C
C     PART 2: THE SS MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)+2
      NTUVSS = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ESS0 COEFFICIENTS (IPHS=+1)
      CALL CPU_TIME(TDM1)
      IF(EQFILE) THEN
        DO ITUV=1,NTUVSS
          IAD = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXAB
          DO M=1,MAXAB
            E11(M,ITUV) = DCMPLX(E0SSFL(IAD+M,1),E0SSFL(IAD+M,2))
            E21(M,ITUV) = DCMPLX(E0SSFL(IAD+M,3),E0SSFL(IAD+M,4))
          ENDDO
        ENDDO
      ELSE
        CALL EQSSMK(E11,E21,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,0)
      ENDIF
      CALL CPU_TIME(TDM2)
      TESS = TESS+TDM2-TDM1
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
          EIJ   = EXL(IBAS,1)+EXL(JBAS,2)
          EROOT = DSQRT(PI/EIJ)**3
          SSS(IBAS,JBAS,1) = EROOT*E11(M,1)
          SSS(IBAS,JBAS,3) = EROOT*E21(M,1)
          SSS(IBAS,JBAS,2) =-PHS*DCONJG(SSS(IBAS,JBAS,3))
          SSS(IBAS,JBAS,4) = PHS*DCONJG(SSS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
500   CONTINUE
C
C**********************************************************************C
C     WE NOW HAVE ALL PIECES OF THE OVERLAP MATRIX FOR THIS BLOCK OF   C
C     BASIS FUNCTIONS -- NOW OVERLAY THE RESULTS INTO OVLP.            C
C**********************************************************************C
C
C     CALCULATE COMPONENT OFFSETS
      IL1 = LRGE(ICNTA,KA,MJA  )
      IL2 = LRGE(ICNTA,KA,MJA+1)
      JL1 = LRGE(ICNTB,KB,MJB  )
      JL2 = LRGE(ICNTB,KB,MJB+1)
C
      IS1 = IL1+NSKP
      IS2 = IL2+NSKP
      JS1 = JL1+NSKP
      JS2 = JL2+NSKP
C
C     LL BLOCKS
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=1,NBAS(1)
            OVLP(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            OVLP(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            OVLP(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            OVLP(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
            OVLP(JL1+JBAS,IL1+IBAS) = DCONJG(OVLP(IL1+IBAS,JL1+JBAS))
            OVLP(JL2+JBAS,IL1+IBAS) = DCONJG(OVLP(IL1+IBAS,JL2+JBAS))
            OVLP(JL1+JBAS,IL2+IBAS) = DCONJG(OVLP(IL2+IBAS,JL1+JBAS))
            OVLP(JL2+JBAS,IL2+IBAS) = DCONJG(OVLP(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=JBAS,NBAS(1)
            OVLP(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            OVLP(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            OVLP(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            OVLP(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
            OVLP(JL1+JBAS,IL1+IBAS) = DCONJG(OVLP(IL1+IBAS,JL1+JBAS))
            OVLP(JL2+JBAS,IL1+IBAS) = DCONJG(OVLP(IL1+IBAS,JL2+JBAS))
            OVLP(JL1+JBAS,IL2+IBAS) = DCONJG(OVLP(IL2+IBAS,JL1+JBAS))
            OVLP(JL2+JBAS,IL2+IBAS) = DCONJG(OVLP(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     NON-RELATIVISTIC OVERLAP MATRIX COMPLETE
      IF(HMLT.EQ.'NORL') GOTO 600
C
C     SS BLOCKS
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=1,NBAS(1)
            OVLP(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            OVLP(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            OVLP(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            OVLP(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
            OVLP(JS1+JBAS,IS1+IBAS) = DCONJG(OVLP(IS1+IBAS,JS1+JBAS))
            OVLP(JS2+JBAS,IS1+IBAS) = DCONJG(OVLP(IS1+IBAS,JS2+JBAS))
            OVLP(JS1+JBAS,IS2+IBAS) = DCONJG(OVLP(IS2+IBAS,JS1+JBAS))
            OVLP(JS2+JBAS,IS2+IBAS) = DCONJG(OVLP(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=JBAS,NBAS(1)
            OVLP(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            OVLP(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            OVLP(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            OVLP(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
            OVLP(JS1+JBAS,IS1+IBAS) = DCONJG(OVLP(IS1+IBAS,JS1+JBAS))
            OVLP(JS2+JBAS,IS1+IBAS) = DCONJG(OVLP(IS1+IBAS,JS2+JBAS))
            OVLP(JS1+JBAS,IS2+IBAS) = DCONJG(OVLP(IS2+IBAS,JS1+JBAS))
            OVLP(JS2+JBAS,IS2+IBAS) = DCONJG(OVLP(IS2+IBAS,JS2+JBAS))
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
C  ONEEL CONSTRUCTS A FULL SET OF MULTI-CENTRE OVERLAP, KINETIC AND    C
C  NUCLEAR ATTRACTION BASIS FUNCTION MATRIX ELEMENTS.                  C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
      DIMENSION RC(MB2,MRC),CP(MB2,3),APH(MB2),PNC(MB2)
      DIMENSION RN(MBS)
C
      COMPLEX*16 TA1,TA2
      COMPLEX*16 CTMP1,CTMP2,CTMP3,CTMP4
      COMPLEX*16 E11A,E11B,E11C,TRM11,E21A,E21B,E21C,TRM21
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
      COMPLEX*16 SLL(MBS,MBS,4),SSS(MBS,MBS,4),
     &           VLL(MBS,MBS,4),VSS(MBS,MBS,4),
     &           TLL(MBS,MBS,4),TLS(MBS,MBS,4),TSL(MBS,MBS,4)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VANM(MDM,MDM),
     &           VSLF(MDM,MDM),VUEH(MDM,MDM),VWKR(MDM,MDM),
     &           VKSB(MDM,MDM),QDIR(MDM,MDM),QXCH(MDM,MDM),
     &           WDIR(MDM,MDM),WXCH(MDM,MDM),CPLE(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/E0LL/E0LLFL(MFL,4),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,4),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VANM,VSLF,
     &            VUEH,VWKR,VKSB,QDIR,QXCH,WDIR,WXCH,CPLE
      COMMON/TMMD/TELL,TESS,TELS,TESL,TRLL,TRSS,TRLS,TRSL,TRBR
C
C     INITIALISE STORAGE MATRICES
      DO I=1,NDIM
        DO J=1,NDIM
          HNUC(I,J) = DCMPLX(0.0D0,0.0D0)
          HKIN(I,J) = DCMPLX(0.0D0,0.0D0)
          OVLP(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)      C
C**********************************************************************C
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
      DO 2000 ICNTB=1,ICNTA
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
        LQN(1) = LVAL(KQN(1))
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
        LQN(2) = LVAL(KQN(2))
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
C     CONSTRUCTION OF ONE-ELECTRON MATRICES BY TT' BLOCKS...
C
C     THIS PHASE RELATES EQ22 AND EQ12 COEFFS TO EQ11 AND EQ21
      PHS = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
     &                    *DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXM = NBAS(1)*NBAS(2)
C
C     INITIALISE STORAGE ARRAYS
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          DO IM=1,4
            SLL(IBAS,JBAS,IM) = DCMPLX(0.0D0,0.0D0)
            SSS(IBAS,JBAS,IM) = DCMPLX(0.0D0,0.0D0)
            TLL(IBAS,JBAS,IM) = DCMPLX(0.0D0,0.0D0)
            TLS(IBAS,JBAS,IM) = DCMPLX(0.0D0,0.0D0)
            TSL(IBAS,JBAS,IM) = DCMPLX(0.0D0,0.0D0)
            VLL(IBAS,JBAS,IM) = DCMPLX(0.0D0,0.0D0)
            VSS(IBAS,JBAS,IM) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     PART 1: THE LL MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMLL  = LQN(1)+LQN(2)
      NTUVLL = (LAMLL+1)*(LAMLL+2)*(LAMLL+3)/6
C
C     GENERATE ELL0 COEFFICIENTS (IPHS = +1)
      CALL CPU_TIME(TDM1)
      IF(EQFILE) THEN
        DO ITUV=1,NTUVLL
          IAD = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXM
          DO M=1,MAXM
            E11(M,ITUV) = DCMPLX(E0LLFL(IAD+M,1),E0LLFL(IAD+M,2))
            E21(M,ITUV) = DCMPLX(E0LLFL(IAD+M,3),E0LLFL(IAD+M,4))
          ENDDO
        ENDDO
      ELSE
        CALL EQLLMK(E11,E21,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,0)
      ENDIF
      CALL CPU_TIME(TDM2)
      TELL = TELL+TDM2-TDM1
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
          EIJ   = EXL(IBAS,1)+EXL(JBAS,2)
          EROOT = DSQRT(PI/EIJ)**3
          SLL(IBAS,JBAS,1) = EROOT*E11(M,1)
          SLL(IBAS,JBAS,3) = EROOT*E21(M,1)
          SLL(IBAS,JBAS,2) =-PHS*DCONJG(SLL(IBAS,JBAS,3))
          SLL(IBAS,JBAS,4) = PHS*DCONJG(SLL(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     NUCLEAR ATTRACTION MATRIX ELEMENTS
      DO IZ=1,NCNT
C
C       ATOM-CENTRED SPECIAL CODE
        IF(ICNTA.EQ.IZ.AND.ICNTB.EQ.IZ) THEN
C
C         ANGULAR SYMMETRY SELECTION RULES
          IF(KA.NE.KB) GOTO 2500
          IF(MA.NE.MB) GOTO 2500
C
C         GAUSSIAN WIDTH PARAMETER
          GPRM = 1.5D0/(RNUC(IZ)*RNUC(IZ))
C
C         FERMI DISTRIBUTION PARAMETER CFMI
          CFMI = 5.0D0*RNUC(IZ)*RNUC(IZ)/3.0D0
     &         - 7.0D0*PI*PI*AFMI(IZ)*AFMI(IZ)/3.0D0
          CFMI = DSQRT(CFMI)
C
C         NORMALISATION CONSTANTS
          G1 = TWLG-GAMLOG(2*LQN(1)+3)
          R1 = DFLOAT(LQN(1))+1.5D0
          DO IBAS=1,NBAS(1)
            ELOG     = DLOG(2.0D0*EXL(IBAS,1))
            RN(IBAS) = DEXP(0.5D0*(G1+R1*ELOG))
          ENDDO
C
C         LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
          DO IBAS=1,NBAS(1)
            EI = EXL(IBAS,1)
            DO JBAS=1,NBAS(2)
              EJ  = EXL(JBAS,2)
              EIJ = EI+EJ
              EL1 = EIJ**(LQN(1)+1)
C
C             RAW INTEGRAL DEPENDS ON NUCLEAR CHARGE MODEL
              IF(NMDL(IZ).EQ.'POINT') THEN
                ULL = 0.5D0*GAMHLF(2*LQN(1)+2)/EL1
              ELSEIF(NMDL(IZ).EQ.'GAUSS') THEN
                ULL = ERFINT0(LQN(1),EIJ,GPRM)
              ELSEIF(NMDL(IZ).EQ.'FERMI') THEN
                ULL = FMIINT0(LQN(1),EIJ,AFMI(IZ),CFMI)
              ELSEIF(NMDL(IZ).EQ.'UNIFM') THEN
                X   = 5.0D0*EIJ*RNUC(IZ)*RNUC(IZ)/3.0D0
                FAC = 2.0D0*DSQRT(X)*X
                GL3 = GAMLWR(2*LQN(1)+3,X)
                GL5 = GAMLWR(2*LQN(1)+5,X)
                GU2 = GAMUPR(2*LQN(1)+2,X)
                ULL = 0.5D0*(GU2 + (3.0D0*X*GL3-GL5)/FAC)/EL1
              ENDIF
C
C             TRANSFER INTO NUCLEAR ATTRACTION MATRIX
              VLL(IBAS,JBAS,1) = VLL(IBAS,JBAS,1)
     &                         - ZNUC(IZ)*RN(IBAS)*RN(JBAS)*ULL
              VLL(IBAS,JBAS,4) = VLL(IBAS,JBAS,4)
     &                         - ZNUC(IZ)*RN(IBAS)*RN(JBAS)*ULL
C
            ENDDO
          ENDDO
C
C       NON-ATOM-CENTRED CODE
        ELSE
C
C         NUCLEAR COORDINATES
          CX = BXYZ(1,IZ)
          CY = BXYZ(2,IZ)
          CZ = BXYZ(3,IZ)
C
C         GAUSSIAN PRODUCT THEOREM OVER BASIS FUNCTIONS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M   = M+1
              EIJ = EXL(IBAS,1)+EXL(JBAS,2)
              PX  = (XYZ(1,1)*EXL(IBAS,1) + XYZ(1,2)*EXL(JBAS,2))/EIJ
              PY  = (XYZ(2,1)*EXL(IBAS,1) + XYZ(2,2)*EXL(JBAS,2))/EIJ
              PZ  = (XYZ(3,1)*EXL(IBAS,1) + XYZ(3,2)*EXL(JBAS,2))/EIJ
              CP(M,1) = CX-PX
              CP(M,2) = CY-PY
              CP(M,3) = CZ-PZ
            ENDDO
          ENDDO
C
          IF(NMDL(IZ).EQ.'POINT') THEN
C         POINT-NUCLEUS APPROXIMATION
C
C           PREPARE ELEMENTS FOR RMAKE
            M = 0
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                M = M+1
C
C               POINT-NUCLEUS EFFECTIVE PARAMETERS
                EIJ    = EXL(IBAS,1)+EXL(JBAS,2)
                APH(M) = EIJ
                PNC(M) = 2.0D0*ZNUC(IZ)*PI/EIJ
C
              ENDDO
            ENDDO
C
C           GENERATE A BATCH OF R-INTEGRALS
            CALL CPU_TIME(TDM1)
            CALL RMAKE(RC,CP,APH,MAXM,LAMLL)
            CALL CPU_TIME(TDM2)
            TRLL = TRLL + TDM2 - TDM1
C
C           NUCLEAR ATTRACTION INTEGRALS (FINITE SUM OF ELL0 AND RC)
            M = 0
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                M = M+1
                DO ITUV=1,NTUVLL
                  VLL(IBAS,JBAS,1) = VLL(IBAS,JBAS,1)
     &                             - PNC(M)*E11(M,ITUV)*RC(M,ITUV)
                  VLL(IBAS,JBAS,3) = VLL(IBAS,JBAS,3)
     &                             - PNC(M)*E21(M,ITUV)*RC(M,ITUV)
                ENDDO
                VLL(IBAS,JBAS,2) =-PHS*DCONJG(VLL(IBAS,JBAS,3))
                VLL(IBAS,JBAS,4) = PHS*DCONJG(VLL(IBAS,JBAS,1))
              ENDDO
            ENDDO
C
          ELSE
C         BEST-FIT EXPANSION
C
C           GAUSSIAN NUCLEAR ATTRACTION
            XI = XNUC(IZ,0)
            FC = FNUC(IZ,0)
C
C           PREPARE ELEMENTS FOR RMAKE
            M = 0
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                M = M+1
C
C               FINITE-NUCLEUS BOYS EXPONENT AND MULTIPLIER
                EIJ    = EXL(IBAS,1)+EXL(JBAS,2)
                ESM    = EIJ+XI
                APH(M) = EIJ*XI/ESM
                PNC(M) = 2.0D0*PI*ZNUC(IZ)*FC*DSQRT(XI/ESM)/EIJ
C
              ENDDO
            ENDDO
C
C           GENERATE A BATCH OF R-INTEGRALS
            CALL CPU_TIME(TDM1)
            CALL RMAKE(RC,CP,APH,MAXM,LAMLL)
            CALL CPU_TIME(TDM2)
            TRLL = TRLL + TDM2 - TDM1
C
C           NUCLEAR ATTRACTION INTEGRALS (FINITE SUM OF ELL0 AND RC)
            M = 0
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                M = M+1
                DO ITUV=1,NTUVLL
                  VLL(IBAS,JBAS,1) = VLL(IBAS,JBAS,1)
     &                             - PNC(M)*E11(M,ITUV)*RC(M,ITUV)
                  VLL(IBAS,JBAS,3) = VLL(IBAS,JBAS,3)
     &                             - PNC(M)*E21(M,ITUV)*RC(M,ITUV)
                ENDDO
                VLL(IBAS,JBAS,2) =-PHS*DCONJG(VLL(IBAS,JBAS,3))
                VLL(IBAS,JBAS,4) = PHS*DCONJG(VLL(IBAS,JBAS,1))
              ENDDO
            ENDDO
C
C           CORRECTION FUNCTIONS FOR DETAILED NUCLEI
            IF(NNUC(IZ).LE.0) GOTO 100
            M = 0
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                M = M+1
C
C               EXPONENT COMBINATIONS
                EIJ = EXL(IBAS,1)+EXL(JBAS,2)
C
C               GAUSSIAN PRODUCT THEOREM CENTRE OF BASIS FUNCTIONS
                PX = (XYZ(1,1)*EXL(IBAS,1) + XYZ(1,2)*EXL(JBAS,2))/EIJ
                PY = (XYZ(2,1)*EXL(IBAS,1) + XYZ(2,2)*EXL(JBAS,2))/EIJ
                PZ = (XYZ(3,1)*EXL(IBAS,1) + XYZ(3,2)*EXL(JBAS,2))/EIJ
C
C               COORDINATES OF EVALUATION
                RPX = PX-BXYZ(1,IZ)
                RPY = PY-BXYZ(2,IZ)
                RPZ = PZ-BXYZ(3,IZ)
                RSQ = RPX*RPX + RPY*RPY + RPZ*RPZ
C
C               INITIALISE TEMPORARY STORAGE ARRAYS
                TA1 = DCMPLX(0.0D0,0.0D0)
                TA2 = DCMPLX(0.0D0,0.0D0)
C
C               CORRECTION FUNCTIONS FOR DETAILED NUCLEI
                DO IFT=1,NNUC(IZ)
C
C                 GAUSSIAN EXPONENT AND FRACTIONAL CONTRIBUTION
                  XI = XNUC(IZ,IFT)
                  FC = FNUC(IZ,IFT)
C
C                 GAUSSIAN COMPONENT OF HGTF
                  ESM = EIJ+XI
                  ECP = EIJ*XI/ESM
                  GSS = DEXP(-ECP*RSQ)
                  AMP =-0.5D0*XI/PI
C
C                 LOOP OVER ALL {A,B,C}
                  DO ITUV=1,NTUVLL
C
C                   CARTESIAN HERMITE POLYNOMIAL VALUES
                    HALPH = HERMITE(ECP,RPX,IA(ITUV))
                    HBETA = HERMITE(ECP,RPY,IB(ITUV))
                    HGAMA = HERMITE(ECP,RPZ,IC(ITUV))
C
C                   HGTF IS A PRODUCT OF THESE
                    HABC = HALPH*HBETA*HGAMA*GSS
C
C                   MULTIPLY BY EQ-COEFF AND ADD TO TOTAL
                    TA1 = TA1 + ZNUC(IZ)*FC*AMP*E11(M,ITUV)*HABC
                    TA2 = TA2 + ZNUC(IZ)*FC*AMP*E21(M,ITUV)*HABC
C
                  ENDDO
C
C               END LOOP OVER NUCLEAR BASIS FOR IZ
                ENDDO
C
C               MATRIX ELEMENTS
                VLL(IBAS,JBAS,1) = VLL(IBAS,JBAS,1) + TA1
                VLL(IBAS,JBAS,3) = VLL(IBAS,JBAS,3) + TA2
                VLL(IBAS,JBAS,2) =-PHS*DCONJG(VLL(IBAS,JBAS,3))
                VLL(IBAS,JBAS,4) = PHS*DCONJG(VLL(IBAS,JBAS,1))
C
C             END LOOP OVER BASIS SET PAIR
              ENDDO
            ENDDO
C
100         CONTINUE
C
C         END CONDITIONAL ON POINT/FINITE NUCLEAR TYPE
          ENDIF
C
C       END CONDITIONAL ON ATOM-CENTRED OR OTHER
        ENDIF
C
C       ATOM-CENTRED SELECTION RULE SKIP POINT
2500    CONTINUE
C
C     END LOOP OVER CENTRES IZ
      ENDDO
C
C     CONSTRUCT NON-RELATIVISTIC KINETIC ENERGY INTEGRALS (IOS 91)
      IF(HMLT.EQ.'NORL') THEN
        RL2 = DFLOAT(2*LQN(2)+3)
        M   = 0
        DO IBAS=1,NBAS(1)
          DO JBAS=1,NBAS(2)
            M   = M+1
C
C           GAUSSIAN PRODUCT THEOREM DETAILS
            EJ  = EXL(JBAS,2)
            EIJ = EXL(IBAS,1) + EXL(JBAS,2)
            ERT = DSQRT(PI/EIJ)**3
            PX = (XYZ(1,1)*EXL(IBAS,1) + XYZ(1,2)*EXL(JBAS,2))/EIJ
            PY = (XYZ(2,1)*EXL(IBAS,1) + XYZ(2,2)*EXL(JBAS,2))/EIJ
            PZ = (XYZ(3,1)*EXL(IBAS,1) + XYZ(3,2)*EXL(JBAS,2))/EIJ
C
C           COORDINATES RELATIVE TO NUCLEAR CENTRE
            PBX = PX-XYZ(1,2)
            PBY = PY-XYZ(2,2)
            PBZ = PZ-XYZ(3,2)
            PB2 = PBX*PBX + PBY*PBY + PBZ*PBZ
C
            E0FC = EJ*RL2 - 2.0D0*EJ*EJ*PB2 - 3.0D0*EJ*EJ/EIJ
            E1FC = 4.0D0*EJ*EJ
C
C           TRUNCATE EXPRESSION DEPENDING ON LAMLL VALUE
C           ALL COMBINATIONS ALLOW FOR THE LAMLL = 0 MANIFOLD
            TRM11 = E0FC*E11(M,IABC(0,0,0))
            TRM21 = E0FC*E21(M,IABC(0,0,0))
C           IF LAMLL > 0 PROVIDE SECOND BUNCH OF TERMS
            IF(LAMLL.GE.1) THEN
              E11A = E11(M,IABC(1,0,0))
              E21A = E21(M,IABC(1,0,0))
              E11B = E11(M,IABC(0,1,0))
              E21B = E21(M,IABC(0,1,0))
              E11C = E11(M,IABC(0,0,1))
              E21C = E21(M,IABC(0,0,1))
              TRM11 = TRM11 - E1FC*(PBX*E11A + PBY*E11B + PBZ*E11C)
              TRM21 = TRM21 - E1FC*(PBX*E21A + PBY*E21B + PBZ*E21C)
            ENDIF
C           IF LAMLL > 1 PROVIDE FINAL BUNCH OF TERMS
            IF(LAMLL.GE.2) THEN
              E11A = E11(M,IABC(2,0,0))
              E21A = E21(M,IABC(2,0,0))
              E11B = E11(M,IABC(0,2,0))
              E21B = E21(M,IABC(0,2,0))
              E11C = E11(M,IABC(0,0,2))
              E21C = E21(M,IABC(0,0,2))
              TRM11 = TRM11 - E1FC*(E11A + E11B + E11C)
              TRM21 = TRM21 - E1FC*(E21A + E21B + E21C)
            ENDIF
            TLL(IBAS,JBAS,1) = ERT*TRM11
            TLL(IBAS,JBAS,3) = ERT*TRM21
            TLL(IBAS,JBAS,2) =-PHS*DCONJG(TLL(IBAS,JBAS,3))
            TLL(IBAS,JBAS,4) = PHS*DCONJG(TLL(IBAS,JBAS,1))
C            
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
      LAMSS  = LQN(1)+LQN(2)+2
      NTUVSS = (LAMSS+1)*(LAMSS+2)*(LAMSS+3)/6
C
C     GENERATE ESS0 COEFFICIENTS (IPHS = +1)
      CALL CPU_TIME(TDM1)
      IF(EQFILE) THEN
        DO ITUV=1,NTUVSS
          IAD = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXM
          DO M=1,MAXM
            E11(M,ITUV) = DCMPLX(E0SSFL(IAD+M,1),E0SSFL(IAD+M,2))
            E21(M,ITUV) = DCMPLX(E0SSFL(IAD+M,3),E0SSFL(IAD+M,4))
          ENDDO
        ENDDO
      ELSE
        CALL EQSSMK(E11,E21,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,0)
      ENDIF
      CALL CPU_TIME(TDM2)
      TESS = TESS+TDM2-TDM1
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
          EIJ   = EXL(IBAS,1)+EXL(JBAS,2)
          EROOT = DSQRT(PI/EIJ)**3
          SSS(IBAS,JBAS,1) = EROOT*E11(M,1)
          SSS(IBAS,JBAS,3) = EROOT*E21(M,1)
          SSS(IBAS,JBAS,2) =-PHS*DCONJG(SSS(IBAS,JBAS,3))
          SSS(IBAS,JBAS,4) = PHS*DCONJG(SSS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     NUCLEAR ATTRACTION MATRIX ELEMENTS
      DO IZ=1,NCNT
C
C       ATOM-CENTRED SPECIAL CODE
        IF(ICNTA.EQ.IZ.AND.ICNTB.EQ.IZ) THEN
C
C         ANGULAR SYMMETRY SELECTION RULES
          IF(KA.NE.KB) GOTO 2600
          IF(MA.NE.MB) GOTO 2600
C
C         ORBITAL TERM
          TL = DFLOAT(2*LQN(1)+1)
C
C         GAUSSIAN WIDTH PARAMETER
          GPRM = 1.5D0/(RNUC(IZ)*RNUC(IZ))
C
C         FERMI DISTRIBUTION PARAMETER CFMI
          CFMI = 5.0D0*RNUC(IZ)*RNUC(IZ)/3.0D0
     &         - 7.0D0*PI*PI*AFMI(IZ)*AFMI(IZ)/3.0D0
          CFMI = DSQRT(CFMI)
C
C         NORMALISATION CONSTANTS
          G2 = TWLG-GAMLOG(2*LQN(1)+5)
          R2 = LQN(1)+0.5D0
          DO IBAS=1,NBAS(1)
            ELOG     = DLOG(2.0D0*EXL(IBAS,1))
            RN(IBAS) = DEXP(0.5D0*(G2+R2*ELOG))
          ENDDO
C
C         LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
          DO IBAS=1,NBAS(1)
            EI = EXL(IBAS,1)
            DO JBAS=1,NBAS(2)
              EJ  = EXL(JBAS,2)
              EIJ = EI+EJ
              EL0 = EIJ**(LQN(1))
              EL2 = EIJ**(LQN(1)+2)
              EPR = EI*EJ
C
C             FIRST RAW INTEGRAL DEPENDS ON NUCLEAR CHARGE MODEL
              IF(NMDL(IZ).EQ.'POINT') THEN
                VSA = 2.0D0*EPR*GAMHLF(2*LQN(1)+4)/EL2
              ELSEIF(NMDL(IZ).EQ.'GAUSS') THEN
                VSA = 4.0D0*EPR*ERFINT0(LQN(1)+1,EIJ,GPRM)
              ELSEIF(NMDL(IZ).EQ.'FERMI') THEN
                VSA = 4.0D0*EPR*FMIINT0(LQN(1)+1,EIJ,AFMI(IZ),CFMI)
              ELSEIF(NMDL(IZ).EQ.'UNIFM') THEN
                X   = 5.0D0*EIJ*RNUC(IZ)*RNUC(IZ)/3.0D0
                FAC = 2.0D0*DSQRT(X)*X
                GL5 = GAMLWR(2*LQN(1)+5,X)
                GL7 = GAMLWR(2*LQN(1)+7,X)
                GU4 = GAMUPR(2*LQN(1)+4,X)
                VSA = 2.0D0*EPR*(GU4 + (3.0D0*X*GL5 - GL7)/FAC)/EL2
              ENDIF
C
C             SECOND RAW INTEGRAL DEPENDS ON NUCLEAR CHARGE AND PARITY
              IF(KQN(1).LT.0) THEN
                VSN = 0.0D0
              ELSE
                IF(NMDL(IZ).EQ.'POINT') THEN
                  VSN = 0.5D0*GAMHLF(2*LQN(1))*TL/EL0
                ELSEIF(NMDL(IZ).EQ.'GAUSS') THEN
                  VSN =-2.0D0*EIJ*TL*ERFINT0(LQN(1)  ,EIJ,GPRM)
     &                        +TL*TL*ERFINT0(LQN(1)-1,EIJ,GPRM)
                ELSEIF(NMDL(IZ).EQ.'FERMI') THEN
                  VSN =-2.0D0*EIJ*TL*FMIINT0(LQN(1)  ,EIJ,AFMI(IZ),CFMI)
     &                        +TL*TL*FMIINT0(LQN(1)-1,EIJ,AFMI(IZ),CFMI)
                ELSEIF(NMDL(IZ).EQ.'UNIFM') THEN
                  X   = 5.0D0*EIJ*RNUC(IZ)*RNUC(IZ)/3.0D0
                  FAC = 2.0D0*DSQRT(X)*X
                  GL1 = GAMLWR(2*LQN(1)+1,X)
                  GL3 = GAMLWR(2*LQN(1)+3,X)
                  GL5 = GAMLWR(2*LQN(1)+5,X)
                  GU0 = GAMUPR(2*LQN(1)  ,X)
                  GU2 = GAMUPR(2*LQN(1)+2,X)
                  VSN =-GU2 - (3.0D0*X*GL3-GL5)/FAC
     &                + 0.5D0*TL*(GU0 + (3.0D0*X*GL1 - GL3)/FAC)
                  VSN = TL*VSN/EL0
                ENDIF
              ENDIF
C
C             SUM OF RAW INTEGRALS
              USS = VSA+VSN
C
C             TRANSFER INTO NUCLEAR ATTRACTION MATRIX
              VSS(IBAS,JBAS,1) = VSS(IBAS,JBAS,1)
     &                         - ZNUC(IZ)*RN(IBAS)*RN(JBAS)*USS
              VSS(IBAS,JBAS,4) = VSS(IBAS,JBAS,4)
     &                         - ZNUC(IZ)*RN(IBAS)*RN(JBAS)*USS
C
            ENDDO
          ENDDO
C
C       NON-ATOM-CENTRED CODE
        ELSE
C
C         NUCLEAR COORDINATES
          CX = BXYZ(1,IZ)
          CY = BXYZ(2,IZ)
          CZ = BXYZ(3,IZ)
C
C         GAUSSIAN PRODUCT THEOREM OVER BASIS FUNCTIONS
          M = 0
          DO IBAS=1,NBAS(1)
            DO JBAS=1,NBAS(2)
              M   = M+1
              EIJ = EXL(IBAS,1)+EXL(JBAS,2)
              PX  = (XYZ(1,1)*EXL(IBAS,1) + XYZ(1,2)*EXL(JBAS,2))/EIJ
              PY  = (XYZ(2,1)*EXL(IBAS,1) + XYZ(2,2)*EXL(JBAS,2))/EIJ
              PZ  = (XYZ(3,1)*EXL(IBAS,1) + XYZ(3,2)*EXL(JBAS,2))/EIJ
              CP(M,1) = CX-PX
              CP(M,2) = CY-PY
              CP(M,3) = CZ-PZ
            ENDDO
          ENDDO
C
          IF(NMDL(IZ).EQ.'POINT') THEN
C         POINT-NUCLEUS APPROXIMATION
C
C           PREPARE ELEMENTS FOR RMAKE
            M = 0
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                M = M+1
C
C               POINT-NUCLEUS EFFECTIVE PARAMETERS
                EIJ = EXL(IBAS,1)+EXL(JBAS,2)
                APH(M) = EIJ
                PNC(M) = 2.0D0*ZNUC(IZ)*PI/EIJ
C
              ENDDO
            ENDDO
C
C           GENERATE A BATCH OF R-INTEGRALS
            CALL CPU_TIME(TDM1)
            CALL RMAKE(RC,CP,APH,MAXM,LAMSS)
            CALL CPU_TIME(TDM2)
            TRSS = TRSS + TDM2 - TDM1
C
C           NUCLEAR ATTRACTION INTEGRALS (FINITE SUM OF ESS0 AND RC)
            M = 0
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                M = M+1
                DO ITUV=1,NTUVSS
                  VSS(IBAS,JBAS,1) = VSS(IBAS,JBAS,1)
     &                             - PNC(M)*E11(M,ITUV)*RC(M,ITUV)
                  VSS(IBAS,JBAS,3) = VSS(IBAS,JBAS,3)
     &                             - PNC(M)*E21(M,ITUV)*RC(M,ITUV)
                ENDDO
                VSS(IBAS,JBAS,2) =-PHS*DCONJG(VSS(IBAS,JBAS,3))
                VSS(IBAS,JBAS,4) = PHS*DCONJG(VSS(IBAS,JBAS,1))
              ENDDO
            ENDDO
C
          ELSE
C         BEST-FIT EXPANSION
C
C           GAUSSIAN NUCLEAR ATTRACTION
            XI = XNUC(IZ,0)
            FC = FNUC(IZ,0)
C
C           PREPARE ELEMENTS FOR RMAKE
            M = 0
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                M = M+1
C
C               FINITE-NUCLEUS BOYS EXPONENT AND MULTIPLIER
                EIJ    = EXL(IBAS,1)+EXL(JBAS,2)
                ESM    = EIJ+XI
                APH(M) = EIJ*XI/ESM
                PNC(M) = 2.0D0*PI*ZNUC(IZ)*FC*DSQRT(XI/ESM)/EIJ
C
              ENDDO
            ENDDO
C
C           GENERATE A BATCH OF R-INTEGRALS
            CALL CPU_TIME(TDM1)
            CALL RMAKE(RC,CP,APH,MAXM,LAMSS)
            CALL CPU_TIME(TDM2)
            TRSS = TRSS + TDM2 - TDM1
C
C           NUCLEAR ATTRACTION INTEGRALS (FINITE SUM OF ESS0 AND RC)
            M = 0
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                M = M+1
                DO ITUV=1,NTUVSS
                  VSS(IBAS,JBAS,1) = VSS(IBAS,JBAS,1)
     &                             - PNC(M)*E11(M,ITUV)*RC(M,ITUV)
                  VSS(IBAS,JBAS,3) = VSS(IBAS,JBAS,3)
     &                             - PNC(M)*E21(M,ITUV)*RC(M,ITUV)
                ENDDO
                VSS(IBAS,JBAS,2) =-PHS*DCONJG(VSS(IBAS,JBAS,3))
                VSS(IBAS,JBAS,4) = PHS*DCONJG(VSS(IBAS,JBAS,1))
              ENDDO
            ENDDO
C
C           CORRECTION FUNCTIONS FOR DETAILED NUCLEI
            IF(NNUC(IZ).LE.0) GOTO 200
            M = 0
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
                M = M+1
C
C               EXPONENT COMBINATIONS
                EIJ = EXL(IBAS,1)+EXL(JBAS,2)
C
C               GAUSSIAN PRODUCT THEOREM CENTRE OF BASIS FUNCTIONS
                PX = (XYZ(1,1)*EXL(IBAS,1) + XYZ(1,2)*EXL(JBAS,2))/EIJ
                PY = (XYZ(2,1)*EXL(IBAS,1) + XYZ(2,2)*EXL(JBAS,2))/EIJ
                PZ = (XYZ(3,1)*EXL(IBAS,1) + XYZ(3,2)*EXL(JBAS,2))/EIJ
C
C               COORDINATES OF EVALUATION
                RPX = PX-BXYZ(1,IZ)
                RPY = PY-BXYZ(2,IZ)
                RPZ = PZ-BXYZ(3,IZ)
                RSQ = RPX*RPX + RPY*RPY + RPZ*RPZ
C
C               INITIALISE TEMPORARY STORAGE ARRAYS
                TA1 = DCMPLX(0.0D0,0.0D0)
                TA2 = DCMPLX(0.0D0,0.0D0)
C
C               CORRECTION FUNCTIONS FOR DETAILED NUCLEI
                DO IFT=1,NNUC(IZ)
C
C                 GAUSSIAN EXPONENT AND FRACTIONAL CONTRIBUTION
                  XI = XNUC(IZ,IFT)
                  FC = FNUC(IZ,IFT)
C
C                 GAUSSIAN COMPONENT OF HGTF
                  ESM = EIJ+XI
                  ECP = EIJ*XI/ESM
                  GSS = DEXP(-ECP*RSQ)
                  AMP =-0.5D0*XI/PI
C
C                 LOOP OVER ALL {A,B,C}
                  DO ITUV=1,NTUVSS
C
C                   CARTESIAN HERMITE POLYNOMIAL VALUES
                    HALPH = HERMITE(ECP,RPX,IA(ITUV))
                    HBETA = HERMITE(ECP,RPY,IB(ITUV))
                    HGAMA = HERMITE(ECP,RPZ,IC(ITUV))
C
C                   HGTF IS A PRODUCT OF THESE
                    HABC = HALPH*HBETA*HGAMA*GSS
C
C                   MULTIPLY BY EQ-COEFF AND ADD TO TOTAL
                    TA1 = TA1 + ZNUC(IZ)*FC*AMP*E11(M,ITUV)*HABC
                    TA2 = TA2 + ZNUC(IZ)*FC*AMP*E21(M,ITUV)*HABC
C
                  ENDDO
C
C               END LOOP OVER NUCLEAR BASIS FOR IZ
                ENDDO
C
C               MATRIX ELEMENTS
                VSS(IBAS,JBAS,1) = VSS(IBAS,JBAS,1) + TA1
                VSS(IBAS,JBAS,3) = VSS(IBAS,JBAS,3) + TA2
                VSS(IBAS,JBAS,2) =-PHS*DCONJG(VSS(IBAS,JBAS,3))
                VSS(IBAS,JBAS,4) = PHS*DCONJG(VSS(IBAS,JBAS,1))
C
C             END LOOP OVER BASIS SET PAIR
              ENDDO
            ENDDO
C
200         CONTINUE
C
C         END CONDITIONAL ON POINT/FINITE NUCLEAR TYPE
          ENDIF
C
C       END CONDITIONAL ON ATOM-CENTRED OR OTHER
        ENDIF
C
C       ATOM-CENTRED SELECTION RULE SKIP POINT
2600    CONTINUE
C
C     END LOOP OVER CENTRES IZ
      ENDDO
C
C     SUBTRACT THE SS OVERLAP MATRIX AND FINISH CONSTRUCTION
      CV2 = 2.0D0*EMSS*CV*CV
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          VSS(IBAS,JBAS,1) = VSS(IBAS,JBAS,1) - CV2*SSS(IBAS,JBAS,1)
          VSS(IBAS,JBAS,3) = VSS(IBAS,JBAS,3) - CV2*SSS(IBAS,JBAS,3)
          VSS(IBAS,JBAS,2) =-PHS*DCONJG(VSS(IBAS,JBAS,3))
          VSS(IBAS,JBAS,4) = PHS*DCONJG(VSS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C**********************************************************************C
C     PART 3: THE SL MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMSS  = LQN(1)+LQN(2)+2
      NTUVSS = (LAMSS+1)*(LAMSS+2)*(LAMSS+3)/6
C
C     KINETIC MATRIX ELEMENTS
      FACT = CV*DSQRT(DFLOAT(2*LQN(2)+3))
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
          EIJ    = EXL(IBAS,1) + EXL(JBAS,2)
          EJRT   = FACT*DSQRT(EXL(JBAS,2))
          EROOT  = DSQRT(PI/EIJ)**3
          TSL(IBAS,JBAS,1) = EJRT*EROOT*E11(M,1)
          TSL(IBAS,JBAS,3) = EJRT*EROOT*E21(M,1)
          TSL(IBAS,JBAS,2) =-PHS*DCONJG(TSL(IBAS,JBAS,3))
          TSL(IBAS,JBAS,4) = PHS*DCONJG(TSL(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C
C**********************************************************************C
C     PART 4: THE LS MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMSS  = LQN(1)+LQN(2)+2
      NTUVSS = (LAMSS+1)*(LAMSS+2)*(LAMSS+3)/6
C
C     GENERATE ESS0 COEFFICIENTS (IPHS = +1)
C
      CALL CPU_TIME(TDM1)
C
C     TODO: THE INDEX SWAP '2 1' MAKES FILE IMPORT ANNOYING. I'M LAZY.
C             JUST GENERATE THESE EQ'S AS A BATCH (DOESN'T TAKE LONG).
C
C     IF(EQFILE) THEN
C       DO ITUV=1,NTUVSS
C         IAD = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB) + (ITUV-1)*MAXM
C         M = 0
C         DO IBAS=1,NBAS(1)
C           DO JBAS=1,NBAS(2)
C             M = M+1
C             N = (JBAS-1)*NBAS(2) + IBAS
C             PI^{SL}_{IJ} = PI^{LS}_{JI}* (THE NEXT LINES ARE WRONG)
C             E11(M,ITUV) = DCMPLX(E0SSFL(IAD+N,1),-E0SSFL(IAD+N,2))
C             E21(M,ITUV) = DCMPLX(E0SSFL(IAD+N,3),-E0SSFL(IAD+N,4))
C           ENDDO
C         ENDDO
C       ENDDO
C     ELSE
        CALL EQSSMK(E11,E21,EXL,XYZ,KQN,MQN,NBAS,+1,2,1,0)
C     ENDIF
      CALL CPU_TIME(TDM2)
      TESS = TESS+TDM2-TDM1
C
C     KINETIC MATRIX ELEMENTS
      FACT = CV*DSQRT(DFLOAT(2*LQN(1)+3))
      M = 0
      DO JBAS=1,NBAS(2)
        DO IBAS=1,NBAS(1)
          M = M+1
          EIJ    = EXL(JBAS,2) + EXL(IBAS,1)
          EIRT   = FACT*DSQRT(EXL(IBAS,1))
          EROOT  = DSQRT(PI/EIJ)**3
          TLS(IBAS,JBAS,1) = EIRT*EROOT*E11(M,1)
          TLS(IBAS,JBAS,3) = EIRT*EROOT*E21(M,1)
          TLS(IBAS,JBAS,2) =-PHS*DCONJG(TLS(IBAS,JBAS,3))
          TLS(IBAS,JBAS,4) = PHS*DCONJG(TLS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     GENERATE LS MATRICES FROM THE ABOVE SL MATRICES
      M = 0
      DO JBAS=1,NBAS(2)
        DO IBAS=1,NBAS(1)
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
      IL1 = LRGE(ICNTA,KA,MJA  )
      IL2 = LRGE(ICNTA,KA,MJA+1)
      JL1 = LRGE(ICNTB,KB,MJB  )
      JL2 = LRGE(ICNTB,KB,MJB+1)
C
      IS1 = IL1+NSKP
      IS2 = IL2+NSKP
      JS1 = JL1+NSKP
      JS2 = JL2+NSKP
C
C     LL OVERLAP BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=1,NBAS(1)
            OVLP(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            OVLP(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            OVLP(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            OVLP(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
C
            OVLP(JL1+JBAS,IL1+IBAS) = DCONJG(OVLP(IL1+IBAS,JL1+JBAS))
            OVLP(JL2+JBAS,IL1+IBAS) = DCONJG(OVLP(IL1+IBAS,JL2+JBAS))
            OVLP(JL1+JBAS,IL2+IBAS) = DCONJG(OVLP(IL2+IBAS,JL1+JBAS))
            OVLP(JL2+JBAS,IL2+IBAS) = DCONJG(OVLP(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=JBAS,NBAS(1)
            OVLP(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            OVLP(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            OVLP(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            OVLP(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
C
            OVLP(JL1+JBAS,IL1+IBAS) = DCONJG(OVLP(IL1+IBAS,JL1+JBAS))
            OVLP(JL2+JBAS,IL1+IBAS) = DCONJG(OVLP(IL1+IBAS,JL2+JBAS))
            OVLP(JL1+JBAS,IL2+IBAS) = DCONJG(OVLP(IL2+IBAS,JL1+JBAS))
            OVLP(JL2+JBAS,IL2+IBAS) = DCONJG(OVLP(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     LL NUCLEAR POTENTIAL BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=1,NBAS(1)
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
        DO JBAS=1,NBAS(2)
          DO IBAS=JBAS,NBAS(1)
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
      IF(HMLT.EQ.'NORL') THEN
C
C       LL KINETIC BLOCKS
        IF(IL1.GT.JL1) THEN
          DO JBAS=1,NBAS(2)
            DO IBAS=1,NBAS(1)
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
          DO JBAS=1,NBAS(2)
            DO IBAS=JBAS,NBAS(1)
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
        DO JBAS=1,NBAS(2)
          DO IBAS=1,NBAS(1)
            OVLP(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            OVLP(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            OVLP(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            OVLP(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
C
            OVLP(JS1+JBAS,IS1+IBAS) = DCONJG(OVLP(IS1+IBAS,JS1+JBAS))
            OVLP(JS2+JBAS,IS1+IBAS) = DCONJG(OVLP(IS1+IBAS,JS2+JBAS))
            OVLP(JS1+JBAS,IS2+IBAS) = DCONJG(OVLP(IS2+IBAS,JS1+JBAS))
            OVLP(JS2+JBAS,IS2+IBAS) = DCONJG(OVLP(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=JBAS,NBAS(1)
            OVLP(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            OVLP(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            OVLP(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            OVLP(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
C
            OVLP(JS1+JBAS,IS1+IBAS) = DCONJG(OVLP(IS1+IBAS,JS1+JBAS))
            OVLP(JS2+JBAS,IS1+IBAS) = DCONJG(OVLP(IS1+IBAS,JS2+JBAS))
            OVLP(JS1+JBAS,IS2+IBAS) = DCONJG(OVLP(IS2+IBAS,JS1+JBAS))
            OVLP(JS2+JBAS,IS2+IBAS) = DCONJG(OVLP(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     SS NUCLEAR POTENTIAL BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NBAS(2)
          DO IBAS=1,NBAS(1)
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
        DO JBAS=1,NBAS(2)
          DO IBAS=JBAS,NBAS(1)
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
        DO JBAS=1,NBAS(2)
          DO IBAS=1,NBAS(1)
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
        DO JBAS=1,NBAS(2)
          DO IBAS=1,NBAS(1)
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
      SUBROUTINE ERI(RR,XYZ,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN)
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
C  ERI GENERATES A BATCH OF MOLECULAR ELECTRON REPULSION INTEGRALS BY  C
C  MEANS OF THE MCMURCHIE-DAVIDSION ALGORITHM (DOUBLE FINITE SUM OVER  C
C  EQ-COEFFICIENTS AND INTEGRALS OVER A PAIR OF HGTFS.)                C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ XYZ  - FULL SET OF CARTESIAN BASIS CENTRES.                       C
C  ▶ KQN  - FULL SET OF RELATIVISTIC LABELS.                           C
C  ▶ MQN  - FULL SET OF MAGNETIC QUANTUM NUMBERS (MAGNITUDE).          C
C  ▶ NBAS - FULL SET OF EXPONENT LIST LENGTHS.                         C
C  ▶ EXL  - FULL LISTS OF EXPONENTS IN THE BLOCK.                      C
C  ▶ IBAS - 1ST BASIS FUNCTION (HELD CONSTANT DURING ROUTINE).         C
C  ▶ JBAS - 2ND BASIS FUNCTION (HELD CONSTANT DURING ROUTINE).         C
C  ▶ ITN  - COMPONENT OVERLAP (T,T') FOR CD. ITN(I) = {LL,LS,SL,SS}.   C
C  OUTPUT:                                                             C
C  ▶ RR   - ERI'S FOR BLOCK AB, ALL 16 MQN SIGN COMBINATIONS.          C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL L0CASE
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4),LQN(4),ITN(2)
      DIMENSION PQ(MB2,3),APH(MB2),PRE(MB2),RC(MB2,MRC)
      DIMENSION IABR11(MEQ),IABR21(MEQ),IABI11(MEQ),IABI21(MEQ)
      DIMENSION ICDR11(MEQ),ICDR21(MEQ),ICDI11(MEQ),ICDI21(MEQ)
      DIMENSION IRC(MRC)
      DIMENSION RCTTFL(20*MFL),IRCTTFL(MFL)
      DIMENSION GABR11(MB2,MEQ),GABI11(MB2,MEQ),
     &          GABR21(MB2,MEQ),GABI21(MB2,MEQ)
      DIMENSION QR1(MB2),QI1(MB2),QR2(MB2),QI2(MB2)
C
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 EAB11(MB2,MEQ),ECD11(MB2,MEQ),
     &           EAB21(MB2,MEQ),ECD21(MB2,MEQ)
C
      COMMON/E0LL/E0LLFL(MFL,4),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,4),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS,IABSL,ICDSL
      COMMON/IRCM/IEAB,IECD,IRIJ(MBS,MBS)
      COMMON/ISCR/IMTX(MB2,11),ISCR(MB2),IMAP(MB2),IBCH,ITOG,MAXN
      COMMON/LSHF/SHLEV(4),SHLV,ILEV
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/TMMD/TELL,TESS,TELS,TESL,TRLL,TRSS,TRLS,TRSL,TRBR
      COMMON/TSCF/TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRW,TCC1,TCC2,TCMC,
     &            TB1B,TB1R,TB1F,TB1M,TBEC,TBRM,TBRW,TBC1,TBC2,TBMC,
     &            TSMX,TUMX,THMX,TAMX,TC1T,TC2T,TCVT,TB1T,TB2T,TACC,
     &            TEIG,TSCR,TTOT,TC1S,TC2S,TB1S,TB2S
C
C     EQ-COEFFICIENT SENSITIVITY PARAMETER
      DATA SENS/1.0D-10/
C
C     ILLEGAL COMPONENT OVERLAP CHECKER
      DO IT=1,2
        IF(ITN(IT).NE.1.AND.ITN(IT).NE.4) THEN
          WRITE(6, *) 'In ERI: illegal component overlaps in ITN.'
          WRITE(7, *) 'In ERI: illegal component overlaps in ITN.'
          STOP
        ENDIF
      ENDDO
C
C     EVALUATE LQNS FOR BASIS FUNCTIONS (A,B,C,D)
      DO N=1,4
        LQN(N) = LVAL(KQN(N))
      ENDDO
C
C     SPECIAL CASE FOR S-TYPE OVERLAPS (ONLY EVER NEEDED ONCE)
      IF(LQN(1)+LQN(2)+LQN(3)+LQN(4).EQ.0) THEN
        L0CASE = .TRUE.
      ELSE
        L0CASE = .FALSE.
      ENDIF
C
C     INTEGRAL SKIPPING ON MOLECULAR GROUP SYMMETRY CLASS BASIS
      IF(SHAPE.EQ.'ATOMIC') THEN
        ISYM = 2
      ELSEIF(SHAPE.EQ.'DIATOM'.OR.SHAPE.EQ.'LINEAR') THEN
        ISYM = 1
      ELSE
        ISYM = 0
      ENDIF
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NBAS(1)*NBAS(2)
      MAXCD = NBAS(3)*NBAS(4)
C
C     PHASE FACTORS FOR AB AND CD PAIR OVERLAPS
      IPHSAB = 1
      IPHSCD =-1
C
C     MCMURCHIE-DAVIDSON MAXINUM ORDER FOR EQ(AB)-COEFFICIENTS
      IF(ITN(1).EQ.1) THEN
        LAMAB = LQN(1)+LQN(2)
      ELSEIF(ITN(1).EQ.4) THEN
        LAMAB = LQN(1)+LQN(2)+2
      ENDIF
C
C     MCMURCHIE-DAVIDSON MAXIMUM ORDER FOR EQ(CD)-COEFFICIENTS
      IF(ITN(2).EQ.1) THEN
        LAMCD = LQN(3)+LQN(4)
      ELSEIF(ITN(2).EQ.4) THEN
        LAMCD = LQN(3)+LQN(4)+2
      ENDIF
C
C     MCMURCHIE-DAVIDSON MAXIMUM ORDER FOR CONTRACTED R-INTEGRAL BATCH
      LAMABCD = LAMAB+LAMCD
C
C     MAXIMUM LAMBDA ORDER FOR RCTTFL SAVED LIST (ONLY IF MCNT.GT.2)
      IF(.NOT.RCFILE.OR.HMLT.EQ.'NORL') THEN
        LAMABCDFL = LAMABCD
      ELSEIF(RCFILE.AND.HMLT.NE.'NORL') THEN
        IF(ILEV.EQ.1) THEN
          LAMABCDFL = LQN(1)+LQN(2)+LQN(3)+LQN(4)
        ELSEIF(ILEV.EQ.2) THEN
          LAMABCDFL = LQN(1)+LQN(2)+LQN(3)+LQN(4)+2
        ELSEIF(ILEV.EQ.3) THEN
          LAMABCDFL = LQN(1)+LQN(2)+LQN(3)+LQN(4)+4
        ENDIF
      ENDIF
C
C     UTILISE EXPANSION SUBSET (WHEN RCFILE.EQ.FALSE)
      IF((ITN(1)+ITN(2))/3.EQ.ILEV) THEN
        ISKIP = 0
      ELSE
        ISKIP = 1
      ENDIF
C
C     PAIR AND GROUP EQ-COEFFICIENT LIST AND R-INTEGRAL BATCH LENGTHS
      NTUVAB     = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
      NTUVCD     = (LAMCD+1)*(LAMCD+2)*(LAMCD+3)/6
      NTUVABCD   = (LAMABCD+1)*(LAMABCD+2)*(LAMABCD+3)/6
      NTUVABCDFL = (LAMABCDFL+1)*(LAMABCDFL+2)*(LAMABCDFL+3)/6
C
C     LIST ADDRESS FOR (AB|  ) AND GAUSSIAN EXPONENT FOR AB OVERLAP
      IJ  = (IBAS-1)*NBAS(2)+JBAS
      EIJ = EXL(IBAS,1)+EXL(JBAS,2)
C
C     INITIALISE RR ARRAY
      DO M=1,MAXCD
        DO ITG=1,16
          RR(M,ITG) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     GENERATE NEW BATCH OF E(AB|  ) COEFFICIENTS IF PROMPTED          C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
      IF(IEAB.EQ.0) GOTO 100
C
C     START TIME
      CALL CPU_TIME(TDM1)
C
      IF(EQFILE) THEN
C
C       OPTION 1: READ FROM LOCAL EQ-COEFFICIENT FILE
        IF(ITN(1).EQ.1) THEN
          DO ITUV=1,NTUVAB
            IAD = IABLL + (ITUV-1)*MAXAB
            DO M=1,MAXAB
              EAB11(M,ITUV) = DCMPLX(E0LLFL(IAD+M,1),E0LLFL(IAD+M,2))
              EAB21(M,ITUV) = DCMPLX(E0LLFL(IAD+M,3),E0LLFL(IAD+M,4))
            ENDDO
          ENDDO
        ELSEIF(ITN(1).EQ.4) THEN
          DO ITUV=1,NTUVAB
            IAD = IABSS + (ITUV-1)*MAXAB
            DO M=1,MAXAB
              EAB11(M,ITUV) = DCMPLX(E0SSFL(IAD+M,1),E0SSFL(IAD+M,2))
              EAB21(M,ITUV) = DCMPLX(E0SSFL(IAD+M,3),E0SSFL(IAD+M,4))
            ENDDO
          ENDDO
        ENDIF
C
      ELSE
C
C       OPTION 2: CALCULATE FROM SCRATCH
        IF(ITN(1).EQ.1) THEN
          CALL EQLLMK(EAB11,EAB21,EXL,XYZ,KQN,MQN,NBAS,IPHSAB,1,2,0)
        ELSEIF(ITN(1).EQ.4) THEN
          CALL EQSSMK(EAB11,EAB21,EXL,XYZ,KQN,MQN,NBAS,IPHSAB,1,2,0)
        ENDIF
C
      ENDIF
C
C     RECORD THE TIME TAKEN TO GENERATE/READ THE E(AB|  ) COEFFICIENTS
      CALL CPU_TIME(TDM2)
      IF(ITN(1).EQ.1) THEN
        TELL = TELL+TDM2-TDM1
      ELSEIF(ITN(1).EQ.4) THEN
        TESS = TESS+TDM2-TDM1
      ENDIF
C
C     SCREENING: TEST E(AB| -) COLUMNS OF CARTESIAN INDEX (T ,U ,V )
      DO IAB=1,NTUVAB
C
C       Re{E(AB|--)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXAB
          ERAB11 = DREAL(EAB11(M,IAB))
          SUM = SUM + DABS(ERAB11)
          IF(SUM.GT.SENS) THEN
            IABR11(IAB) = 1
            GOTO 101
          ENDIF
        ENDDO
        IABR11(IAB) = 0
101     CONTINUE
C
C       Im{E(AB|--)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXAB
          EIAB11 = DIMAG(EAB11(M,IAB))
          SUM = SUM + DABS(EIAB11)
          IF(SUM.GT.SENS) THEN
            IABI11(IAB) = 1
            GOTO 102
          ENDIF
        ENDDO
        IABI11(IAB) = 0
102     CONTINUE
C
C       Re{E(AB|+-)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXAB
          ERAB21 = DREAL(EAB21(M,IAB))
          SUM = SUM + DABS(ERAB21)
          IF(SUM.GT.SENS) THEN
            IABR21(IAB) = 1
            GOTO 103
          ENDIF
        ENDDO
        IABR21(IAB) = 0
103     CONTINUE
C
C       Im{E(AB|+-)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXAB
          EIAB21 = DIMAG(EAB21(M,IAB))
          SUM = SUM + DABS(EIAB21)
          IF(SUM.GT.SENS) THEN
            IABI21(IAB) = 1
            GOTO 104
          ENDIF
        ENDDO
        IABI21(IAB) = 0
104     CONTINUE
C
      ENDDO
C
C     DO NOT CALCULATE AGAIN UNTIL PROMPTED EXTERNALLY
      IEAB = 0
C
100   CONTINUE
C
C**********************************************************************C
C     GENERATE NEW BATCH OF E(CD| -) COEFFICIENTS IF PROMPTED          C
C**********************************************************************C
C
      IF(IECD.EQ.0) GOTO 200
C
C     START TIME
      CALL CPU_TIME(TDM1)
C
C     OPTION 1: READ FROM LOCAL EQ-COEFFICIENT FILE
      IF(EQFILE) THEN
        IF(ITN(2).EQ.1) THEN
          DO ITUV=1,NTUVCD
            IAD = ICDLL + (ITUV-1)*MAXCD
            Z   = DFLOAT((-1)**(ILAM(ITUV)))
            DO M=1,MAXCD
              ECD11(M,ITUV) = Z*DCMPLX(E0LLFL(IAD+M,1),E0LLFL(IAD+M,2))
              ECD21(M,ITUV) = Z*DCMPLX(E0LLFL(IAD+M,3),E0LLFL(IAD+M,4))
            ENDDO
          ENDDO
        ELSEIF(ITN(2).EQ.4) THEN
          DO ITUV=1,NTUVCD
            IAD = ICDSS + (ITUV-1)*MAXCD
            Z   = DFLOAT((-1)**(ILAM(ITUV)))
            DO M=1,MAXCD
              ECD11(M,ITUV) = Z*DCMPLX(E0SSFL(IAD+M,1),E0SSFL(IAD+M,2))
              ECD21(M,ITUV) = Z*DCMPLX(E0SSFL(IAD+M,3),E0SSFL(IAD+M,4))
            ENDDO
          ENDDO
        ENDIF
C
C     OPTION 2: CALCULATE FROM SCRATCH
      ELSE
        IF(ITN(2).EQ.1) THEN
          CALL EQLLMK(ECD11,ECD21,EXL,XY Z,KQN,MQN,NBAS,IPHSCD,3,4,0)
        ELSEIF(ITN(2).EQ.4) THEN
          CALL EQSSMK(ECD11,ECD21,EXL,XYZ,KQN,MQN,NBAS,IPHSCD,3,4,0)
        ENDIF
      ENDIF
C
C     RECORD THE TIME TAKEN TO GENERATE/READ THE E(CD|  ) COEFFICIENTS
      CALL CPU_TIME(TDM2)
      IF(ITN(2).EQ.1) THEN
        TELL = TELL+TDM2-TDM1
      ELSEIF(ITN(2).EQ.4) THEN
        TESS = TESS+TDM2-TDM1
      ENDIF
C
C     SCREENING: TEST E(CD| -) COLUMNS OF CARTESIAN INDEX (T',U',V')
      DO ICD=1,NTUVCD
C
C       Re{E(CD|--)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXCD
          ERCD11 = DREAL(ECD11(M,ICD))
          SUM = SUM + DABS(ERCD11)
          IF(SUM.GT.SENS) THEN
            ICDR11(ICD) = 1
            GOTO 201
          ENDIF
        ENDDO
        ICDR11(ICD) = 0
201     CONTINUE
C
C       Im{E(CD|--)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXCD
          EICD11 = DIMAG(ECD11(M,ICD))
          SUM = SUM + DABS(EICD11)
          IF(SUM.GT.SENS) THEN
            ICDI11(ICD) = 1
            GOTO 202
          ENDIF
        ENDDO
        ICDI11(ICD) = 0
202     CONTINUE
C
C       Re{E(CD|+-)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXCD
          ERCD21 = DREAL(ECD21(M,ICD))
          SUM = SUM + DABS(ERCD21)
          IF(SUM.GT.SENS) THEN
            ICDR21(ICD) = 1
            GOTO 203
          ENDIF
        ENDDO
        ICDR21(ICD) = 0
203     CONTINUE
C
C       Im{E(CD|+-)} COEFFICIENTS
        SUM = 0.0D0
        DO M=1,MAXCD
          EICD21 = DIMAG(ECD21(M,ICD))
          SUM = SUM + DABS(EICD21)
          IF(SUM.GT.SENS) THEN
            ICDI21(ICD) = 1
            GOTO 204
          ENDIF
        ENDDO
        ICDI21(ICD) = 0
204     CONTINUE
C
      ENDDO
C
C     DO NOT CALCULATE AGAIN UNTIL PROMPTED EXTERNALLY
      IECD = 0
C
200   CONTINUE
      CALL CPU_TIME(T2)
      TCEC = TCEC+T2-T1
C
C**********************************************************************C
C     GENERATE NEW BATCH OF RC(AB|CD) INTEGRALS IF PROMPTED            C
C**********************************************************************C
C
C     START TIME
      CALL CPU_TIME(TDM1)
C
C     SKIP IF A HIGHER LEVEL OF THIS BATCH EXISTS IN RC
      IF(.NOT.RCFILE.AND.ISKIP.EQ.0) GOTO 300
C
C     SKIP IF A HIGHER LEVEL OF THIS BATCH EXISTS IN RC
      IF(L0CASE.AND.ISKIP.EQ.0) GOTO 300
C
C     SKIP IF INTEGRAL BATCH EXISTS IN FILE
      IF(RCFILE.AND.IRIJ(IBAS,JBAS).EQ.0) GOTO 300
C
      CALL CPU_TIME(T1)
C
C     GAUSSIAN OVERLAP CENTRE
      PX = (XYZ(1,1)*EXL(IBAS,1)+XYZ(1,2)*EXL(JBAS,2))/EIJ
      PY = (XYZ(2,1)*EXL(IBAS,1)+XYZ(2,2)*EXL(JBAS,2))/EIJ
      PZ = (XYZ(3,1)*EXL(IBAS,1)+XYZ(3,2)*EXL(JBAS,2))/EIJ
C
C     AUXILLIARY DATA FOR RMAKE ROUTINE
      M = 0
      N = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
          IF(RCFILE.AND.L0CASE.AND.ISCR(M).EQ.0) GOTO 301
          IF(.NOT.RCFILE.AND.ISCR(M).EQ.0) GOTO 301
          N   = N+1
          EKL = EXL(KBAS,3)+EXL(LBAS,4)
          QX  = (XYZ(1,3)*EXL(KBAS,3)+XYZ(1,4)*EXL(LBAS,4))/EKL
          QY  = (XYZ(2,3)*EXL(KBAS,3)+XYZ(2,4)*EXL(LBAS,4))/EKL
          QZ  = (XYZ(3,3)*EXL(KBAS,3)+XYZ(3,4)*EXL(LBAS,4))/EKL
          PQ(N,1) = QX-PX
          PQ(N,2) = QY-PY
          PQ(N,3) = QZ-PZ
          APH(N)  = EIJ*EKL/(EIJ+EKL)
301       CONTINUE
        ENDDO
      ENDDO
C
C     BATCH SIZE AND EXPANSION LENGTH DEPENDS ON MODE
      IF(RCFILE.AND..NOT.L0CASE) THEN
        MBCH = MAXCD
        MLAM = LAMABCDFL
        MNTV = NTUVABCDFL
      ELSE
        MBCH = MAXN
        MLAM = LAMABCD
        MNTV = NTUVABCD
      ENDIF
C
C     GENERATE R-INTEGRALS
      CALL RMAKE(RC,PQ,APH,MBCH,MLAM)
C
C     SCREENING: TEST RC(AB|CD) COLUMNS WITH INDEX (T+T',U+U',V+V')
      DO INTV=1,MNTV
C
C       SUM OF RC(AB|CD) MAGNITUDES
        SUM = 0.0D0
        DO N=1,MBCH
          SUM = SUM + DABS(RC(N,INTV))
          IF(SUM.GT.SENS) THEN
            IRC(INTV) = 1
            GOTO 302
          ENDIF
        ENDDO
        IRC(INTV) = 0
302     CONTINUE
C
      ENDDO
C
      CALL CPU_TIME(T2)
      TCRM = TCRM+T2-T1
C
C     CONTINUE ONLY IF INTEGRALS ARE TO BE SAVED TO LARGE FILE
      IF(L0CASE) GOTO 300
      IF(.NOT.RCFILE) GOTO 300
C
      CALL CPU_TIME(T3)
C
C     TEST WHETHER FINAL ADDRESS IS STILL INSIDE ARRAY BOUNDS
      IF(20*MFL.LT.IJ*MAXCD*NTUVABCDFL) THEN
C       OUT OF BOUNDS: PRINT WARNING BUT KEEP GOING
        WRITE(6, *) 'In ERI: RCTT words exceed allocated limit.'
        WRITE(7, *) 'In ERI: RCTT words exceed allocated limit.'
        GOTO 300
      ELSE
C       DO NOT CALCULATE AGAIN UNTIL PROMPTED EXTERNALLY
        IRIJ(IBAS,JBAS) = 0
      ENDIF
C
C     STARTING ADDRESS FOR THIS BATCH OF SAVED R(AB|CD) INTEGRALS
      IADRTT = (IJ-1)*MAXCD*NTUVABCDFL
C
C     COPY THIS BATCH OF INTEGRALS TO A SAVED LIST
      DO IABCDFL=1,NTUVABCDFL
        IAD = IADRTT + MAXCD*(IABCDFL-1)
        DO M=1,MAXCD
          RCTTFL(IAD+M) = RC(M,IABCDFL)
        ENDDO
      ENDDO
C
C     STARTING ADDRESS FOR THIS BATCH OF SCREENING FLAGS
      IADSCR = (IJ-1)*NTUVABCDFL
C
C     COPY SCREENING MARKERS TO A SAVED LIST
      DO IABCDFL=1,NTUVABCDFL
        IRCTTFL(IADSCR+IABCDFL) = IRC(IABCDFL)
      ENDDO
C
      CALL CPU_TIME(T4)
      TCRW = TCRW+T4-T3
C
300   CONTINUE
C
C     RECORD THE TIME TAKEN TO GENERATE THE RC(AB|CD) BATCH
      CALL CPU_TIME(TDM2)
      IF(ITN(1).EQ.1.AND.ITN(2).EQ.1) THEN
        TRLL = TRLL+TDM2-TDM1
      ELSEIF(ITN(1).EQ.4.AND.ITN(2).EQ.4) THEN
        TRSS = TRSS+TDM2-TDM1
      ELSE
        TRLS = TRLS+TDM2-TDM1
      ENDIF
C
C**********************************************************************C
C     PERFORM FIRST CONTRACTION: G(AB| -) = E(CD| -)*RC(AB|CD).        C
C     THIS YIELDS ALL MQN SIGN POSSIBILITIES FOR C AND D.              C
C**********************************************************************C
C
C     TIME AT START OF FIRST CONTRACTION
      CALL CPU_TIME(T1)
C
C     LOOP OVER ALL ADDRESSES FOR E(AB| -) FINITE EXPANSION
      DO IAB=1,NTUVAB
C
C       RESET CONTRACTION STORAGE ARRAYS G(AB| -)
        DO N=1,MAXN
          GABR11(N,IAB) = 0.0D0
          GABI11(N,IAB) = 0.0D0
          GABR21(N,IAB) = 0.0D0
          GABI21(N,IAB) = 0.0D0
        ENDDO
C
C       SKIP ENTIRE PROCESS IF E(AB| -) PASSES SCREENING CONDITION
        IABALL = IABR11(IAB)+IABI11(IAB)+IABR21(IAB)+IABI21(IAB)
        IF(IABALL.EQ.0) GOTO 401
C
C       CALCULATE DIRECTLY FROM SMALL RC(AB|CD) LOCAL ARRAY
        IF(.NOT.RCFILE.OR.L0CASE) THEN
C
C         LOOP OVER ALL FINITE EXPANSION ADDRESSES FOR E(CD| -)
          DO ICD=1,NTUVCD
C
C           CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
            IRABCD = IABC(IA(IAB)+IA(ICD),IB(IAB)+IB(ICD),
     &                                    IC(IAB)+IC(ICD))
C
C           SKIP THIS STEP IF THE RC(AB|CD) PASSES SCREENING CONDITION
            IF(IRC(IRABCD).EQ.0) GOTO 402
C
C           SKIP THIS STEP IF THE E(CD) PASSES SCREENING CONDITION
            ICDALL = ICDR11(ICD)+ICDI11(ICD)+ICDR21(ICD)+ICDI21(ICD)
            IF(ICDALL.EQ.0) GOTO 402
C
C           CONTRIBUTIONS TO Re{G(AB|--)} FROM EACH Re{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABR11(IAB).EQ.0) GOTO 411
            IF(ISYM.EQ.2.AND.IABR11(IAB).EQ.0) GOTO 411
            IF(ICDR11(ICD).EQ.1) THEN
              DO N=1,MAXN
                GABR11(N,IAB) = GABR11(N,IAB)
     &                         + DREAL(ECD11(IMAP(N),ICD))*RC(N,IRABCD)
              ENDDO
            ENDIF
411         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|--)} FROM EACH Im{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABI11(IAB).EQ.0) GOTO 412
            IF(ISYM.EQ.2.AND.IABI11(IAB).EQ.0) GOTO 412
            IF(ICDI11(ICD).EQ.1) THEN
              DO N=1,MAXN
                GABI11(N,IAB) = GABI11(N,IAB)
     &                         + DIMAG(ECD11(IMAP(N),ICD))*RC(N,IRABCD)
              ENDDO
            ENDIF
412         CONTINUE
C
C           CONTRIBUTIONS TO Re{G(AB|+-)} FROM EACH Re{E(CD|+-)} ADDRESS
            IF(ISYM.EQ.1.AND.IABR21(IAB).EQ.0) GOTO 413
            IF(ISYM.EQ.2.AND.IABR21(IAB).EQ.0) GOTO 413
            IF(ICDR21(ICD).EQ.1) THEN
              DO N=1,MAXN
                GABR21(N,IAB) = GABR21(N,IAB)
     &                         + DREAL(ECD21(IMAP(N),ICD))*RC(N,IRABCD)
              ENDDO
            ENDIF
413         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|+-)} FROM EACH Im{E(CD|+-)} ADDRESS
            IF(ISYM.EQ.1.AND.IABI21(IAB).EQ.0) GOTO 414
            IF(ISYM.EQ.2.AND.IABI21(IAB).EQ.0) GOTO 414
            IF(ICDI21(ICD).EQ.1) THEN
              DO N=1,MAXN
                GABI21(N,IAB) = GABI21(N,IAB)
     &                         + DIMAG(ECD21(IMAP(N),ICD))*RC(N,IRABCD)
              ENDDO
            ENDIF
414         CONTINUE
C
C           SKIP POINT FOR RC(AB|CD) SCREENING
402         CONTINUE
C
C         END LOOP OVER E(CD|  ) FINITE EXPANSION ADDRESSES
          ENDDO
C
C       CALCULATE DIRECTLY FROM LARGE RC(AB|CD) SCRATCH ARRAY
        ELSE
C
C         LOOP OVER ALL FINITE EXPANSION ADDRESSES FOR E(CD| -)
          DO ICD=1,NTUVCD
C
C           CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
            IRABCD = IABC(IA(IAB)+IA(ICD),IB(IAB)+IB(ICD),
     &                                    IC(IAB)+IC(ICD))
C
C           STARTING ADDRESS FOR THIS BATCH OF SAVED R(AB|CD) INTEGRALS
            IAD = (IJ-1)*MAXCD*NTUVABCDFL + MAXCD*(IRABCD-1)
C
C           SKIP THIS STEP IF THE RC(AB|CD) PASSES SCREENING CONDITION
            IF(IRCTTFL((IJ-1)*NTUVABCDFL+IRABCD).EQ.0) GOTO 432
C
C           SKIP THIS STEP IF THE E(CD) PASSES SCREENING CONDITION
            ICDALL = ICDR11(ICD)+ICDI11(ICD)+ICDR21(ICD)+ICDI21(ICD)
            IF(ICDALL.EQ.0) GOTO 432
C
C           CONTRIBUTIONS TO Re{G(AB|--)} FROM EACH Re{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABR11(IAB).EQ.0) GOTO 421
            IF(ISYM.EQ.2.AND.IABR11(IAB).EQ.0) GOTO 421
            IF(ICDR11(ICD).EQ.1) THEN
              DO N=1,MAXN
                GABR11(N,IAB) = GABR11(N,IAB)
     &                  + DREAL(ECD11(IMAP(N),ICD))*RCTTFL(IAD+IMAP(N))
              ENDDO
            ENDIF
421         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|--)} FROM EACH Im{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABI11(IAB).EQ.0) GOTO 422
            IF(ISYM.EQ.2.AND.IABI11(IAB).EQ.0) GOTO 422
            IF(ICDI11(ICD).EQ.1) THEN
              DO N=1,MAXN
                GABI11(N,IAB) = GABI11(N,IAB)
     &                  + DIMAG(ECD11(IMAP(N),ICD))*RCTTFL(IAD+IMAP(N))
              ENDDO
            ENDIF
422         CONTINUE
C
C           CONTRIBUTIONS TO Re{G(AB|+-)} FROM EACH Re{E(CD|+-)} ADDRESS
            IF(ISYM.EQ.1.AND.IABR21(IAB).EQ.0) GOTO 423
            IF(ISYM.EQ.2.AND.IABR21(IAB).EQ.0) GOTO 423
            IF(ICDR21(ICD).EQ.1) THEN
              DO N=1,MAXN
                GABR21(N,IAB) = GABR21(N,IAB)
     &                  + DREAL(ECD21(IMAP(N),ICD))*RCTTFL(IAD+IMAP(N))
              ENDDO
            ENDIF
423         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|+-)} FROM EACH Im{E(CD|+-)} ADDRESS
            IF(ISYM.EQ.1.AND.IABI21(IAB).EQ.0) GOTO 424
            IF(ISYM.EQ.2.AND.IABI21(IAB).EQ.0) GOTO 424
            IF(ICDI21(ICD).EQ.1) THEN
              DO N=1,MAXN
                GABI21(N,IAB) = GABI21(N,IAB)
     &                  + DIMAG(ECD21(IMAP(N),ICD))*RCTTFL(IAD+IMAP(N))
              ENDDO
            ENDIF
424         CONTINUE
C
C           SKIP POINT FOR RC(AB|CD) SCREENING
432         CONTINUE
C
C         END LOOP OVER E(CD|  ) FINITE EXPANSION ADDRESSES
          ENDDO
C
        ENDIF
C
C       SKIP POINT FOR E(AB|  ) SCREENING
401     CONTINUE
C
C     END LOOP OVER E(AB|  ) FINITE EXPANSION ADDRESSES
      ENDDO
C
C     TIME AT END OF FIRST CONTRACTION
      CALL CPU_TIME(T2)
      TCC1 = TCC1+T2-T1
C
C**********************************************************************C
C     PERFORM SECOND CONTRACTION: ( -| -) = E(AB| -)*G(AB| -).         C
C     THIS YIELDS A FULL BATCH OF TWO-ELECTRON INTEGRALS (16 PERM'NS). C
C**********************************************************************C
C
C     CALCULATE PHASES FOR BASIS FUNCTION OVERLAP COMBINATIONS
      PAB = ISIGN(1,KQN(1)*KQN(2))*(-1)**((MQN(1)-MQN(2))/2)
      PCD = ISIGN(1,KQN(3)*KQN(4))*(-1)**((MQN(3)-MQN(4))/2)
C
      PABCD = PAB*PCD
C
C     1ST SET: ( 1) = (--|--)   ( 4) = (--|++)
C              (16) = (++|++)   (13) = (++|--)
C
C     RESET CONTRACTION STORAGE LISTS
      DO N=1,MAXN
        QR1(N) = 0.0D0
        QI1(N) = 0.0D0
        QR2(N) = 0.0D0
        QI2(N) = 0.0D0
      ENDDO
C
C     RAW CONTRACTION (--|--) = E(AB|--)*(Re{G(AB|--)} + i*Im{G(AB|--)})
      DO IAB=1,NTUVAB
        IF(IABR11(IAB).EQ.1) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB11(IJ,IAB))*GABR11(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB11(IJ,IAB))*GABI11(N,IAB)
          ENDDO
        ENDIF
        IF(IABI11(IAB).EQ.1) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB11(IJ,IAB))*GABR11(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB11(IJ,IAB))*GABI11(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     APPLY PHASE RELATIONS AND NORMALISATION FACTORS TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N, 1) =     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N, 4) = PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
        RR(N,16) = PABCD*DCONJG(RR(N, 1))
        RR(N,13) = PABCD*DCONJG(RR(N, 4))
      ENDDO
C
      IF(ISYM.EQ.1) GOTO 501
      IF(ISYM.EQ.2) GOTO 501
C
C     2ND SET: ( 3) = (--|+-)   ( 2) = (--|-+)
C              (14) = (++|-+)   (15) = (++|+-)
C
C     RESET CONTRACTION STORAGE LISTS
      DO N=1,MAXN
        QR1(N) = 0.0D0
        QI1(N) = 0.0D0
        QR2(N) = 0.0D0
        QI2(N) = 0.0D0
      ENDDO
C
C     RAW CONTRACTION (--|+-) = E(AB|--)*(Re{G(AB|+-)} + i*Im{G(AB|+-)})
      DO IAB=1,NTUVAB
        IF(IABR11(IAB).EQ.1) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB11(IJ,IAB))*GABR21(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB11(IJ,IAB))*GABI21(N,IAB)
          ENDDO
        ENDIF
        IF(IABI11(IAB).EQ.1) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB11(IJ,IAB))*GABR21(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB11(IJ,IAB))*GABI21(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     APPLY PHASE RELATIONS AND NORMALISATION FACTORS TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N, 3) =     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N, 2) =-PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
        RR(N,14) =-PABCD*DCONJG(RR(N, 3))
        RR(N,15) =-PABCD*DCONJG(RR(N, 2))
      ENDDO
C
C     3RD SET: ( 9) = (+-|--)   (12) = (+-|++)
C              ( 8) = (-+|++)   ( 5) = (-+|--)
C
C     RESET CONTRACTION STORAGE LISTS
      DO N=1,MAXN
        QR1(N) = 0.0D0
        QI1(N) = 0.0D0
        QR2(N) = 0.0D0
        QI2(N) = 0.0D0
      ENDDO
C
C     RAW CONTRACTION (+-|--) = E(AB|+-)*(Re{G(AB|--)} + i*Im{G(AB|--)})
      DO IAB=1,NTUVAB
        IF(IABR21(IAB).EQ.1) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB21(IJ,IAB))*GABR11(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB21(IJ,IAB))*GABI11(N,IAB)
          ENDDO
        ENDIF
        IF(IABI21(IAB).EQ.1) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB21(IJ,IAB))*GABR11(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB21(IJ,IAB))*GABI11(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     APPLY PHASE RELATIONS AND NORMALISATION FACTORS TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N, 9) =     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N,12) = PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
        RR(N, 8) =-PABCD*DCONJG(RR(N, 9))
        RR(N, 5) =-PABCD*DCONJG(RR(N,12))
      ENDDO
C
501   CONTINUE
C
C     4TH SET: (11) = (+-|+-)   (10) = (+-|-+)
C              ( 6) = (-+|-+)   ( 7) = (-+|+-)
C
C     RESET CONTRACTION STORAGE LISTS
      DO N=1,MAXN
        QR1(N) = 0.0D0
        QI1(N) = 0.0D0
        QR2(N) = 0.0D0
        QI2(N) = 0.0D0
      ENDDO
C
C     RAW CONTRACTION (+-|+-) = E(AB|+-)*(Re{G(AB|+-)} + i*Im{G(AB|+-)})
      DO IAB=1,NTUVAB
        IF(IABR21(IAB).EQ.1) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB21(IJ,IAB))*GABR21(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB21(IJ,IAB))*GABI21(N,IAB)
          ENDDO
        ENDIF
        IF(IABI21(IAB).EQ.1) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB21(IJ,IAB))*GABR21(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB21(IJ,IAB))*GABI21(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     APPLY PHASE RELATIONS AND NORMALISATION FACTORS TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N,11) =     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N,10) =-PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
        RR(N, 6) = PABCD*DCONJG(RR(N,11))
        RR(N, 7) = PABCD*DCONJG(RR(N,10))
      ENDDO
C
C     TIME AT END OF SECOND CONTRACTION
      CALL CPU_TIME(T3)
      TCC2 = TCC2+T3-T2
C
C**********************************************************************C
C     COULOMB INTEGRAL BATCH NOW FULLY CONSTRUCTED                     C
C**********************************************************************C
C
C     CALCULATE THE R-INTEGRAL NORMALISATION FACTOR
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
          IF(ISCR(M).EQ.1) THEN
            EKL = EXL(KBAS,3)+EXL(LBAS,4)
            EMX = DSQRT(EIJ+EKL)*EIJ*EKL
            PRE(M) = 2.0D0*PI52/EMX
          ENDIF
        ENDDO
      ENDDO
C
C     INCLUDE THE R-INTEGRAL NORMALISATION FACTOR
      DO N=1,MAXN
        DO ITG=1,16
          RR(N,ITG) = PRE(IMAP(N))*RR(N,ITG)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE BII(RR,XYZ,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN,BGAUNT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                          BBBBBBB IIII IIII                           C
C                          BB    BB II   II                            C
C                          BB    BB II   II                            C
C                          BBBBBBB  II   II                            C
C                          BB    BB II   II                            C
C                          BB    BB II   II                            C
C                          BBBBBBB IIII IIII                           C
C                                                                      C
C -------------------------------------------------------------------- C
C  ERI GENERATES A BATCH OF MOLECULAR ELECTRON REPULSION INTEGRALS BY  C
C  MEANS OF THE MCMURCHIE-DAVIDSION ALGORITHM (DOUBLE FINITE SUM OVER  C
C  EQ-COEFFICIENTS AND INTEGRALS OVER A PAIR OF HGTFS.)                C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ XYZ    - FULL SET OF CARTESIAN BASIS CENTRES.                     C
C  ▶ KQN    - FULL SET OF RELATIVISTIC LABELS.                         C
C  ▶ MQN    - FULL SET OF MAGNETIC QUANTUM NUMBERS (MAGNITUDE).        C
C  ▶ NBAS   - FULL SET OF EXPONENT LIST LENGTHS.                       C
C  ▶ EXL    - FULL LISTS OF EXPONENTS IN THE BLOCK.                    C
C  ▶ IBAS   - 1ST BASIS FUNCTION (HELD CONSTANT DURING ROUTINE).       C
C  ▶ JBAS   - 2ND BASIS FUNCTION (HELD CONSTANT DURING ROUTINE).       C
C  ▶ ITN    - COMPONENT OVERLAP COMBINATION.                           C
C  ▶ BGAUNT - GAUNT INTERACTION OVERRIDE OPTION.                       C
C  OUTPUT:                                                             C
C  ▶ RR     - BII'S FOR BLOCK AB, ALL 16 MQN SIGN COMBINATIONS.        C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL L0CASE,BGAUNT
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4),ITN(2)
      DIMENSION PQ(MB2,3),APH(MB2),PRE(MB2),T(MB2),RC(MB2,MRC)
      DIMENSION IABR11(MEQ,3),IABI11(MEQ,3),IABR21(MEQ,3),IABI21(MEQ,3)
      DIMENSION ICDR11(MEQ,3),ICDI11(MEQ,3),ICDR21(MEQ,3),ICDI21(MEQ,3)
      DIMENSION IRC(MRC)
      DIMENSION RCTTFL(20*MFL),IRCTTFL(MFL)
      DIMENSION GABR11(MB2,MEQ),GABI11(MB2,MEQ),
     &          GABR21(MB2,MEQ),GABI21(MB2,MEQ)
      DIMENSION QR1(MB2),QI1(MB2),QR2(MB2),QI2(MB2)
      DIMENSION IDX(3),JDX(3)
C
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 EAB11(MB2,MEQ,3),EAB21(MB2,MEQ,3),
     &           ECD11(MB2,MEQ,3),ECD21(MB2,MEQ,3)
C
      COMMON/EILS/EILSFL(MFL,12),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EISL/EISLFL(MFL,12),IADISL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS,IABSL,ICDSL
      COMMON/IRCM/IEAB,IECD,IRIJ(MBS,MBS)
      COMMON/ISCR/IMTX(MB2,11),ISCR(MB2),IMAP(MB2),IBCH,ITOG,MAXN
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/TMMD/TELL,TESS,TELS,TESL,TRLL,TRSS,TRLS,TRSL,TRBR
      COMMON/TSCF/TC1B,TC1R,TC1F,TC1M,TCEC,TCRM,TCRW,TCC1,TCC2,TCMC,
     &            TB1B,TB1R,TB1F,TB1M,TBEC,TBRM,TBRW,TBC1,TBC2,TBMC,
     &            TSMX,TUMX,THMX,TAMX,TC1T,TC2T,TCVT,TB1T,TB2T,TACC,
     &            TEIG,TSCR,TTOT,TC1S,TC2S,TB1S,TB2S
C
C     EQ-COEFFICIENT SENSITIVITY PARAMETER
      DATA SENS/1.0D-10/
C
C     ILLEGAL COMPONENT OVERLAP CHECKER
      DO IT=1,2
        IF(ITN(IT).NE.2.AND.ITN(IT).NE.3) THEN
          WRITE(6, *) 'In BII: illegal component overlaps in ITN.'
          WRITE(7, *) 'In BII: illegal component overlaps in ITN.'
          STOP
        ENDIF
      ENDDO
C
C     EVALUATE LQNS FOR BASIS FUNCTIONS (A,B,C,D)
      DO N=1,4
        LQN(N) = LVAL(KQN(N))
      ENDDO
C
C     SPECIAL CASE FOR S-TYPE OVERLAPS (ONLY EVER NEEDED ONCE)
      IF(LQN(1)+LQN(2)+LQN(3)+LQN(4).EQ.0) THEN
        L0CASE = .TRUE.
      ELSE
        L0CASE = .FALSE.
      ENDIF
C
C     INTEGRAL SKIPPING ON MOLECULAR GROUP SYMMETRY CLASS BASIS
      IF(SHAPE.EQ.'ATOMIC') THEN
        ISYM = 2
      ELSEIF(SHAPE.EQ.'DIATOM'.OR.SHAPE.EQ.'LINEAR') THEN
        ISYM = 1
      ELSE
        ISYM = 0
      ENDIF
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NBAS(1)*NBAS(2)
      MAXCD = NBAS(3)*NBAS(4)
C
C     PHASE FACTOR FOR AB AND CD PAIR OVERLAPS
      IPHSAB = 1
      IPHSCD =-1
C
C     MCMURCHIE-DAVIDSON MAXIMUM ORDER FOR EQ-COEFFICIENTS
      LAMAB = LQN(1)+LQN(2)+1
      LAMCD = LQN(3)+LQN(4)+1
C
C     MCMURCHIE-DAVIDSON MAXIMUM ORDER FOR CONTRACTED R-INTEGRAL BATCH
      IF(BGAUNT) THEN
C       GAUNT INTEGRALS ONLY
        LAMABCD = LAMAB+LAMCD
      ELSE
C       FULL BREIT INTERACTION
        LAMABCD = LAMAB+LAMCD+2
      ENDIF
C
C     MCMURCHIE-DAVIDSON EQ-COEFFICIENT AND R-INTEGRAL LIST LENGTHS
      NTUVAB   = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
      NTUVCD   = (LAMCD+1)*(LAMCD+2)*(LAMCD+3)/6
      NTUVABCD = (LAMABCD+1)*(LAMABCD+2)*(LAMABCD+3)/6
C
C     LIST ADDRESS FOR (AB|  ) AND GAUSSIAN EXPONENT FOR AB OVERLAP
      IJ  = (IBAS-1)*NBAS(2)+JBAS
      EIJ = EXL(IBAS,1)+EXL(JBAS,2)
C
C     INITIALISE RR ARRAY
      DO M=1,MAXCD
        DO ITG=1,16
          RR(M,ITG) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     GENERATE NEW BATCH OF E(AB|  ) COEFFICIENTS IF PROMPTED          C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
      IF(IEAB.EQ.0) GOTO 100
C
C     START TIME
      CALL CPU_TIME(TDM1)
C
      IF(EQFILE) THEN
C
C       OPTION 1: READ FROM LOCAL EQ-COEFFICIENT FILE
        IF(ITN(1).EQ.2) THEN
          DO ITUV=1,NTUVAB
            IAD = IABLS + (ITUV-1)*MAXAB
            DO M=1,MAXAB
            EAB11(M,ITUV,1) = DCMPLX(EILSFL(IAD+M, 1),EILSFL(IAD+M, 2))
            EAB21(M,ITUV,1) = DCMPLX(EILSFL(IAD+M, 3),EILSFL(IAD+M, 4))
            EAB11(M,ITUV,2) = DCMPLX(EILSFL(IAD+M, 5),EILSFL(IAD+M, 6))
            EAB21(M,ITUV,2) = DCMPLX(EILSFL(IAD+M, 7),EILSFL(IAD+M, 8))
            EAB11(M,ITUV,3) = DCMPLX(EILSFL(IAD+M, 9),EILSFL(IAD+M,10))
            EAB21(M,ITUV,3) = DCMPLX(EILSFL(IAD+M,11),EILSFL(IAD+M,12))
            ENDDO
          ENDDO
        ELSEIF(ITN(1).EQ.3) THEN
          DO ITUV=1,NTUVAB
            IAD = IABSL + (ITUV-1)*MAXAB
            DO M=1,MAXAB
            EAB11(M,ITUV,1) = DCMPLX(EISLFL(IAD+M, 1),EISLFL(IAD+M, 2))
            EAB21(M,ITUV,1) = DCMPLX(EISLFL(IAD+M, 3),EISLFL(IAD+M, 4))
            EAB11(M,ITUV,2) = DCMPLX(EISLFL(IAD+M, 5),EISLFL(IAD+M, 6))
            EAB21(M,ITUV,2) = DCMPLX(EISLFL(IAD+M, 7),EISLFL(IAD+M, 8))
            EAB11(M,ITUV,3) = DCMPLX(EISLFL(IAD+M, 9),EISLFL(IAD+M,10))
            EAB21(M,ITUV,3) = DCMPLX(EISLFL(IAD+M,11),EISLFL(IAD+M,12))
            ENDDO
          ENDDO
        ENDIF
C
      ELSE
C
C       OPTION 2: CALCULATE FROM SCRATCH
        IF(ITN(1).EQ.2) THEN
          CALL EILSB3(EAB11,EAB21,EXL,XYZ,KQN,MQN,NBAS,IPHSAB,1,2)
        ELSEIF(ITN(1).EQ.3) THEN
          CALL EISLB3(EAB11,EAB21,EXL,XYZ,KQN,MQN,NBAS,IPHSAB,1,2)
        ENDIF
C
      ENDIF
C
C     RECORD THE TIME TAKEN TO GENERATE/READ THE E(AB|  ) COEFFICIENTS
      CALL CPU_TIME(TDM2)
      IF(ITN(1).EQ.2) THEN
        TELS = TELS+TDM2-TDM1
      ELSEIF(ITN(1).EQ.3) THEN
        TESL = TESL+TDM2-TDM1
      ENDIF
C
C     SCREENING PROCEDURE: NORM SUM OF EQ-COEFFICIENT LIST FOR EACH IAB
      DO ICMP=1,3
        DO IAB=1,NTUVAB
C
C         Re{E(AB|--)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXAB
            ER = DREAL(EAB11(M,IAB,ICMP))
            SUM = SUM + DABS(ER)
            IF(SUM.GT.SENS) THEN
              IABR11(IAB,ICMP) = 1
              GOTO 101
            ENDIF
          ENDDO
          IABR11(IAB,ICMP) = 0
101       CONTINUE
C
C         Im{E(AB|--)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXAB
            EI = DIMAG(EAB11(M,IAB,ICMP))
            SUM = SUM + DABS(EI)
            IF(SUM.GT.SENS) THEN
              IABI11(IAB,ICMP) = 1
              GOTO 102
            ENDIF
          ENDDO
          IABI11(IAB,ICMP) = 0
102       CONTINUE
C
C         Re{E(AB|+-)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXAB
            ER = DREAL(EAB21(M,IAB,ICMP))
            SUM = SUM + DABS(ER)
            IF(SUM.GT.SENS) THEN
              IABR21(IAB,ICMP) = 1
              GOTO 103
            ENDIF
          ENDDO
          IABR21(IAB,ICMP) = 0
103       CONTINUE
C
C         Im{E(AB|+-)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXAB
            EI = DIMAG(EAB21(M,IAB,ICMP))
            SUM = SUM + DABS(EI)
            IF(SUM.GT.SENS) THEN
              IABI21(IAB,ICMP) = 1
              GOTO 104
            ENDIF
          ENDDO
          IABI21(IAB,ICMP) = 0
104       CONTINUE
C
        ENDDO
      ENDDO
C
C     DO NOT CALCULATE AGAIN UNTIL PROMPTED EXTERNALLY
      IEAB = 0
C
100   CONTINUE
C
C**********************************************************************C
C     GENERATE NEW BATCH OF E(CD| -) COEFFICIENTS IF PROMPTED          C
C**********************************************************************C
C
      IF(IECD.EQ.0) GOTO 200
C
C     START TIME
      CALL CPU_TIME(TDM1)
C
C     OPTION 1: READ FROM LOCAL EQ-COEFFICIENT FILE
      IF(EQFILE) THEN
        IF(ITN(2).EQ.2) THEN
          DO ITUV=1,NTUVCD
            IAD = ICDLS + (ITUV-1)*MAXCD
            Z   = DFLOAT((-1)**(ILAM(ITUV)))
            DO M=1,MAXCD
            ECD11(M,ITUV,1)=Z*DCMPLX(EILSFL(IAD+M, 1),EILSFL(IAD+M, 2))
            ECD21(M,ITUV,1)=Z*DCMPLX(EILSFL(IAD+M, 3),EILSFL(IAD+M, 4))
            ECD11(M,ITUV,2)=Z*DCMPLX(EILSFL(IAD+M, 5),EILSFL(IAD+M, 6))
            ECD21(M,ITUV,2)=Z*DCMPLX(EILSFL(IAD+M, 7),EILSFL(IAD+M, 8))
            ECD11(M,ITUV,3)=Z*DCMPLX(EILSFL(IAD+M, 9),EILSFL(IAD+M,10))
            ECD21(M,ITUV,3)=Z*DCMPLX(EILSFL(IAD+M,11),EILSFL(IAD+M,12))
            ENDDO
          ENDDO
        ELSEIF(ITN(2).EQ.3) THEN
          DO ITUV=1,NTUVCD
            IAD = ICDSL + (ITUV-1)*MAXCD
            Z   = DFLOAT((-1)**(ILAM(ITUV)))
            DO M=1,MAXCD
            ECD11(M,ITUV,1)=Z*DCMPLX(EISLFL(IAD+M, 1),EISLFL(IAD+M, 2))
            ECD21(M,ITUV,1)=Z*DCMPLX(EISLFL(IAD+M, 3),EISLFL(IAD+M, 4))
            ECD11(M,ITUV,2)=Z*DCMPLX(EISLFL(IAD+M, 5),EISLFL(IAD+M, 6))
            ECD21(M,ITUV,2)=Z*DCMPLX(EISLFL(IAD+M, 7),EISLFL(IAD+M, 8))
            ECD11(M,ITUV,3)=Z*DCMPLX(EISLFL(IAD+M, 9),EISLFL(IAD+M,10))
            ECD21(M,ITUV,3)=Z*DCMPLX(EISLFL(IAD+M,11),EISLFL(IAD+M,12))
            ENDDO
          ENDDO
        ENDIF
C
C     OPTION 2: CALCULATE FROM SCRATCH
      ELSE
        IF(ITN(2).EQ.2) THEN
          CALL EILSB3(ECD11,ECD21,EXL,XYZ,KQN,MQN,NBAS,IPHSCD,3,4)
        ELSEIF(ITN(2).EQ.3) THEN
          CALL EISLB3(ECD11,ECD21,EXL,XYZ,KQN,MQN,NBAS,IPHSCD,3,4)
        ENDIF
      ENDIF
C
C     RECORD THE TIME TAKEN TO GENERATE/READ THE E(CD|  ) COEFFICIENTS
      CALL CPU_TIME(TDM2)
      IF(ITN(2).EQ.2) THEN
        TELS = TELS+TDM2-TDM1
      ELSEIF(ITN(2).EQ.3) THEN
        TESL = TESL+TDM2-TDM1
      ENDIF
C
C     SCREENING PROCEDURE: NORM SUM OF EQ-COEFFICIENT LIST FOR EACH ICD
      DO JCMP=1,3
        DO ICD=1,NTUVCD
C
C         Re{E(CD|--)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXCD
            ER = DREAL(ECD11(M,ICD,JCMP))
            SUM = SUM + DABS(ER)
            IF(SUM.GT.SENS) THEN
              ICDR11(ICD,JCMP) = 1
              GOTO 201
            ENDIF
          ENDDO
          ICDR11(ICD,JCMP) = 0
201       CONTINUE
C
C         Im{E(CD|--)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXCD
            EI = DIMAG(ECD11(M,ICD,JCMP))
            SUM = SUM + DABS(EI)
            IF(SUM.GT.SENS) THEN
              ICDI11(ICD,JCMP) = 1
              GOTO 202
            ENDIF
          ENDDO
          ICDI11(ICD,JCMP) = 0
202       CONTINUE
C
C         Re{E(CD|+-)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXCD
            ER = DREAL(ECD21(M,ICD,JCMP))
            SUM = SUM + DABS(ER)
            IF(SUM.GT.SENS) THEN
              ICDR21(ICD,JCMP) = 1
              GOTO 203
            ENDIF
          ENDDO
          ICDR21(ICD,JCMP) = 0
203       CONTINUE
C
C         Im{E(CD|+-)} COEFFICIENTS
          SUM = 0.0D0
          DO M=1,MAXCD
            EI = DIMAG(ECD21(M,ICD,JCMP))
            SUM = SUM + DABS(EI)
            IF(SUM.GT.SENS) THEN
              ICDI21(ICD,JCMP) = 1
              GOTO 204
            ENDIF
          ENDDO
          ICDI21(ICD,JCMP) = 0
204       CONTINUE
C
        ENDDO
      ENDDO
C
C     DO NOT CALCULATE AGAIN UNTIL ASKED EXTERNALLY
      IECD = 0
C
200   CONTINUE
      CALL CPU_TIME(T2)
      TBEC = TBEC+T2-T1
C
C**********************************************************************C
C     GENERATE NEW BATCH OF RC(AB|CD) INTEGRALS IF PROMPTED            C
C**********************************************************************C
C
C     START TIME
      CALL CPU_TIME(TDM1)
C
C     GAUSSIAN OVERLAP CENTRE
      PX = (XYZ(1,1)*EXL(IBAS,1)+XYZ(1,2)*EXL(JBAS,2))/EIJ
      PY = (XYZ(2,1)*EXL(IBAS,1)+XYZ(2,2)*EXL(JBAS,2))/EIJ
      PZ = (XYZ(3,1)*EXL(IBAS,1)+XYZ(3,2)*EXL(JBAS,2))/EIJ
C
C     AUXILLIARY DATA FOR RMAKE ROUTINE
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M   = M+1
          EKL = EXL(KBAS,3)+EXL(LBAS,4)
          QX  = (XYZ(1,3)*EXL(KBAS,3)+XYZ(1,4)*EXL(LBAS,4))/EKL
          QY  = (XYZ(2,3)*EXL(KBAS,3)+XYZ(2,4)*EXL(LBAS,4))/EKL
          QZ  = (XYZ(3,3)*EXL(KBAS,3)+XYZ(3,4)*EXL(LBAS,4))/EKL
          PQ(M,1) = QX-PX
          PQ(M,2) = QY-PY
          PQ(M,3) = QZ-PZ
          APH(M)  = EIJ*EKL/(EIJ+EKL)
        ENDDO
      ENDDO
C
C     SKIP IF INTEGRAL BATCH EXISTS IN FILE
      IF(RCFILE.AND.IRIJ(IBAS,JBAS).EQ.0) GOTO 300
C
      CALL CPU_TIME(T1)
C
C     BATCH SIZE AND EXPANSION LENGTH DEPENDS ON MODE
      IF(RCFILE.AND..NOT.L0CASE) THEN
        MBCH = MAXCD
      ELSE
        MBCH = MAXN
      ENDIF
C
C     SHORTEN RMAKE DATA IF POSSIBLE
      IF(ITOG.NE.0) THEN
        IF(.NOT.RCFILE.OR.L0CASE) THEN
          DO N=1,MAXN
            PQ(N,1) = PQ(IMAP(N),1)
            PQ(N,2) = PQ(IMAP(N),2)
            PQ(N,3) = PQ(IMAP(N),3)
            APH(N)  = APH(IMAP(N))
          ENDDO
        ENDIF
      ENDIF
C
C     GENERATE R-INTEGRALS
      CALL RMAKE(RC,PQ,APH,MBCH,LAMABCD)
C
C     SCREENING: TEST RC(AB|CD) COLUMNS WITH INDEX (T+T',U+U',V+V')
      DO IABCD=1,NTUVABCD
C
C       SUM OF RC(AB|CD) MAGNITUDES
        SUM = 0.0D0
        DO N=1,MBCH
          SUM = SUM + DABS(RC(N,IABCD))
          IF(SUM.GT.SENS) THEN
            IRC(IABCD) = 1
            GOTO 301
          ENDIF
        ENDDO
        IRC(IABCD) = 0
301     CONTINUE
C
      ENDDO
C
      CALL CPU_TIME(T2)
      TBRM = TBRM+T2-T1
C
C     CONTINUE ONLY IF INTEGRALS ARE TO BE SAVED TO LARGE FILE
      IF(.NOT.RCFILE.OR.L0CASE) GOTO 300
C
      CALL CPU_TIME(T3)
C
C     TEST WHETHER FINAL ADDRESS IS STILL INSIDE ARRAY BOUNDS
      IF(20*MFL.LT.IJ*MAXCD*NTUVABCD) THEN
C       OUT OF BOUNDS: PRINT WARNING BUT KEEP GOING
        WRITE(6, *) 'In BII: RCTT words exceed allocated limit.'
        WRITE(7, *) 'In BII: RCTT words exceed allocated limit.'
        GOTO 300
      ELSE
C       DO NOT CALCULATE AGAIN UNTIL PROMPTED EXTERNALLY
        IRIJ(IBAS,JBAS) = 0
      ENDIF
C
C     STARTING ADDRESS FOR THIS BATCH OF SAVED R(AB|CD) INTEGRALS
      IADRTT = (IJ-1)*MAXCD*NTUVABCD
C
C     COPY THIS BATCH OF INTEGRALS TO A SAVED LIST
      DO IABCD=1,NTUVABCD
        IAD = IADRTT + MAXCD*(IABCD-1)
        DO M=1,MAXCD
          RCTTFL(IAD+M) = RC(M,IABCD)
        ENDDO
      ENDDO
C
C     STARTING ADDRESS FOR SCREENING FLAGS
      IADSCR = (IJ-1)*NTUVABCD
C
C     COPY SCREENING MARKERS TO A SAVED LIST
      DO IABCD=1,NTUVABCD
        IRCTTFL(IADSCR+IABCD) = IRC(IABCD)
      ENDDO
C
      CALL CPU_TIME(T4)
      TBRW = TBRW+T4-T3
C
300   CONTINUE
C
C     RECORD THE TIME TAKEN TO GENERATE THE RC(AB|CD) BATCH
      CALL CPU_TIME(TDM2)
      TRBR = TRBR+TDM2-TDM1
C
C**********************************************************************C
C     PERFORM FIRST CONTRACTION: G(AB| -) = E(CD| -)*RC(AB|CD).        C
C     THIS YIELDS ALL MQN SIGN POSSIBILITIES FOR C AND D.              C
C**********************************************************************C
C
C     LOOP OVER CARTESIAN INDEX ICMP FOR CENTRE AB (USE INDEX 6000)
      DO 6000 ICMP=1,3
C
C     TIME AT START OF FIRST CONTRACTION FOR THIS ICMP INDEX
      CALL CPU_TIME(T1I)
C
C     CARTESIAN INDEX ICMP AS A VECTOR, IDX
      CALL NCART(IDX,ICMP)
C
C     LOOP OVER ALL ADDRESSES FOR E(AB| -) FINITE EXPANSION
      DO IAB=1,NTUVAB
C
C       RESET CONTRACTION STORAGE ARRAYS G(AB| -)
        DO N=1,MAXN
          GABR11(N,IAB) = 0.0D0
          GABI11(N,IAB) = 0.0D0
          GABR21(N,IAB) = 0.0D0
          GABI21(N,IAB) = 0.0D0
        ENDDO
C
C       SKIP ENTIRE PROCESS IF E(AB| -) PASSES SCREENING CONDITION
        IABALL = IABR11(IAB,ICMP)+IABI11(IAB,ICMP)+IABR21(IAB,ICMP)
     &                                            +IABI21(IAB,ICMP)
        IF(IABALL.EQ.0) GOTO 401
C
C >>>>> GAUNT INTERACTION
C
C       LOOP OVER ALL FINITE EXPANSION ADDRESSES FOR E(CD| -)
        DO ICD=1,NTUVCD
C
C         CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
          IRABCD = IABC(IA(IAB)+IA(ICD),IB(IAB)+IB(ICD),IC(IAB)+IC(ICD))
C
C         CALCULATE DIRECTLY FROM RC(AB|CD) LOCAL ARRAY
          IF(.NOT.RCFILE.OR.L0CASE) THEN
C
C           SKIP THIS STEP IF THE RC(AB|CD) PASSES SCREENING CONDITION
            IF(IRC(IRABCD).EQ.0) GOTO 402
C
C           SKIP THIS STEP IF THE E(CD) PASSES SCREENING CONDITION
            ICDALL = ICDR11(ICD,ICMP)+ICDI11(ICD,ICMP)+ICDR21(ICD,ICMP)
     &                                                +ICDI21(ICD,ICMP)
            IF(ICDALL.EQ.0) GOTO 402
C
C           CONTRIBUTIONS TO Re{G(AB|--)} FROM EACH Re{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABR11(IAB,ICMP).EQ.0) GOTO 411
            IF(ISYM.EQ.2.AND.IABR11(IAB,ICMP).EQ.0) GOTO 411
            IF(ICDR11(ICD,ICMP).EQ.1) THEN
              DO N=1,MAXN
                GABR11(N,IAB) = GABR11(N,IAB)
     &                    + DREAL(ECD11(IMAP(N),ICD,ICMP))*RC(N,IRABCD)
              ENDDO
            ENDIF
411         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|--)} FROM EACH Im{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABI11(IAB,ICMP).EQ.0) GOTO 412
            IF(ISYM.EQ.2.AND.IABI11(IAB,ICMP).EQ.0) GOTO 412
            IF(ICDI11(ICD,ICMP).EQ.1) THEN
              DO N=1,MAXN
                GABI11(N,IAB) = GABI11(N,IAB)
     &                    + DIMAG(ECD11(IMAP(N),ICD,ICMP))*RC(N,IRABCD)
              ENDDO
            ENDIF
412         CONTINUE
C
C           CONTRIBUTIONS TO Re{G(AB|+-)} FROM EACH Re{E(CD|+-)} ADDRESS
            IF(ISYM.EQ.1.AND.IABR21(IAB,ICMP).EQ.0) GOTO 413
            IF(ISYM.EQ.2.AND.IABR21(IAB,ICMP).EQ.0) GOTO 413
            IF(ICDR21(ICD,ICMP).EQ.1) THEN
              DO N=1,MAXN
                GABR21(N,IAB) = GABR21(N,IAB)
     &                    + DREAL(ECD21(IMAP(N),ICD,ICMP))*RC(N,IRABCD)
              ENDDO
            ENDIF
413         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|+-)} FROM EACH Im{E(CD|+-)} ADDRESS
            IF(ISYM.EQ.1.AND.IABI21(IAB,ICMP).EQ.0) GOTO 414
            IF(ISYM.EQ.2.AND.IABI21(IAB,ICMP).EQ.0) GOTO 414
            IF(ICDI21(ICD,ICMP).EQ.1) THEN
              DO N=1,MAXN
                GABI21(N,IAB) = GABI21(N,IAB)
     &                    + DIMAG(ECD21(IMAP(N),ICD,ICMP))*RC(N,IRABCD)
              ENDDO
            ENDIF
414         CONTINUE
C
C           SKIP POINT FOR RC(AB|CD) AND E(CD) SCREENING
402         CONTINUE
C
          ELSE
C
C           STARTING ADDRESS FOR THIS BATCH OF SAVED R(AB|CD) INTEGRALS
            IAD = (IJ-1)*MAXCD*NTUVABCD + MAXCD*(IRABCD-1)
C
C           SKIP THIS STEP IF THE RC(AB|CD) PASSES SCREENING CONDITION
            IF(IRCTTFL((IJ-1)*NTUVABCD+IRABCD).EQ.0) GOTO 462
C
C           SKIP THIS STEP IF THE E(CD) PASSES SCREENING CONDITION
            ICDALL = ICDR11(ICD,ICMP)+ICDI11(ICD,ICMP)+ICDR21(ICD,ICMP)
     &                                                +ICDI21(ICD,ICMP)
            IF(ICDALL.EQ.0) GOTO 462
C
C           CONTRIBUTIONS TO Re{G(AB|--)} FROM EACH Re{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABR11(IAB,ICMP).EQ.0) GOTO 451
            IF(ISYM.EQ.2.AND.IABR11(IAB,ICMP).EQ.0) GOTO 451
            IF(ICDR11(ICD,ICMP).EQ.1) THEN
              DO N=1,MAXN
                GABR11(N,IAB) = GABR11(N,IAB)
     &             + DREAL(ECD11(IMAP(N),ICD,ICMP))*RCTTFL(IAD+IMAP(N))
              ENDDO
            ENDIF
451         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|--)} FROM EACH Im{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABI11(IAB,ICMP).EQ.0) GOTO 452
            IF(ISYM.EQ.2.AND.IABI11(IAB,ICMP).EQ.0) GOTO 452
            IF(ICDI11(ICD,ICMP).EQ.1) THEN
              DO N=1,MAXN
                GABI11(N,IAB) = GABI11(N,IAB)
     &             + DIMAG(ECD11(IMAP(N),ICD,ICMP))*RCTTFL(IAD+IMAP(N))
              ENDDO
            ENDIF
452         CONTINUE
C
C           CONTRIBUTIONS TO Re{G(AB|+-)} FROM EACH Re{E(CD|+-)} ADDRESS
            IF(ISYM.EQ.1.AND.IABR21(IAB,ICMP).EQ.0) GOTO 453
            IF(ISYM.EQ.2.AND.IABR21(IAB,ICMP).EQ.0) GOTO 453
            IF(ICDR21(ICD,ICMP).EQ.1) THEN
              DO N=1,MAXN
                GABR21(N,IAB) = GABR21(N,IAB)
     &             + DREAL(ECD21(IMAP(N),ICD,ICMP))*RCTTFL(IAD+IMAP(N))
              ENDDO
            ENDIF
453         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|+-)} FROM EACH Im{E(CD|+-)} ADDRESS
            IF(ISYM.EQ.1.AND.IABI21(IAB,ICMP).EQ.0) GOTO 454
            IF(ISYM.EQ.2.AND.IABI21(IAB,ICMP).EQ.0) GOTO 454
            IF(ICDI21(ICD,ICMP).EQ.1) THEN
              DO N=1,MAXN
                GABI21(N,IAB) = GABI21(N,IAB)
     &             + DIMAG(ECD21(IMAP(N),ICD,ICMP))*RCTTFL(IAD+IMAP(N))
              ENDDO
            ENDIF
454         CONTINUE
C
C           SKIP POINT FOR RC(AB|CD) AND E(CD) SCREENING
462         CONTINUE
C
          ENDIF
C
C >>>>>   GAUGE TERM (REQUIRES ADDITIONAL CARTESIAN SUM Q')
          IF(BGAUNT) GOTO 480
C
C         LOOP OVER CARTESIAN INDEX JCMP FOR CENTRE CD
          DO JCMP=1,3
C
C           SKIP THIS STEP IF THE E(CD) PASSES SCREENING CONDITION
            IF(ICDR11(ICD,JCMP)+ICDI11(ICD,JCMP)
     &                    +ICDR21(ICD,JCMP)+ICDI21(ICD,JCMP).EQ.0) THEN
              GOTO 403
            ENDIF
C
C           CARTESIAN INDEX JCMP AS A VECTOR, JDX
            CALL NCART(JDX,JCMP)
C
C           NEW ADDRESS DEPENDING ON JCMP CARTESIAN INDEX
            IF(JCMP.EQ.1) THEN
              RTP = DFLOAT(IA(IAB)+IA(ICD))
            ELSEIF(JCMP.EQ.2) THEN
              RTP = DFLOAT(IB(IAB)+IB(ICD))
            ELSEIF(JCMP.EQ.3) THEN
              RTP = DFLOAT(IC(IAB)+IC(ICD))
            ENDIF
C
C           FIRST CONTRIBUTION ADDRESS
            I1 = IA(IAB)+IA(ICD)+IDX(1)+JDX(1)
            J1 = IB(IAB)+IB(ICD)+IDX(2)+JDX(2)
            K1 = IC(IAB)+IC(ICD)+IDX(3)+JDX(3)
C
C           CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
            IADR1 = IABC(I1,J1,K1)
C
C           SECOND CONTRIBUTION ADDRESS
            I2 = IA(IAB)+IA(ICD)+IDX(1)
            J2 = IB(IAB)+IB(ICD)+IDX(2)
            K2 = IC(IAB)+IC(ICD)+IDX(3)
C
C           CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
            IADR2 = IABC(I2,J2,K2)
C
C           THIRD CONTRIBUTION ADDRESS
            I3 = IA(IAB)+IA(ICD)+IDX(1)-JDX(1)
            J3 = IB(IAB)+IB(ICD)+IDX(2)-JDX(2)
            K3 = IC(IAB)+IC(ICD)+IDX(3)-JDX(3)
C
C           CALCULATE RC ADDRESS FOR THIS PARTICULAR AB/CD OVERLAP
            IF(I3.GE.0.AND.J3.GE.0.AND.K3.GE.0) THEN
              IADR3 = IABC(I3,J3,K3)
            ELSE
              IADR3 = 0
            ENDIF
C
C           CALCULATE DIRECTLY FROM RC(AB|CD) LOCAL ARRAY
            IF(.NOT.RCFILE.OR.L0CASE) THEN
C
C             SKIP THIS STEP IF THE RC(AB|CD) PASSES SCREENING CONDITION
              IF(IADR3.NE.0) THEN
                IF(IRC(IADR1)+IRC(IADR2)+IRC(IADR3).EQ.0) GOTO 403
              ELSE
                IF(IRC(IADR1)+IRC(IADR2).EQ.0) GOTO 403
              ENDIF
C
C             PRE-FACTORS FOR THE UPCOMING CONTRACTION
              DO N=1,MAXN
                T1 = RC(N,IADR1)*0.5D0/APH(N)
                T2 = RC(N,IADR2)*PQ(N,JCMP)
                IF(I3.GE.0.AND.J3.GE.0.AND.K3.GE.0) THEN
                  T3 = RC(N,IADR3)*RTP
                ELSE
                  T3 = 0.0D0
                ENDIF
                T(N) = T1-T2+T3
              ENDDO
C
            ELSE
C
C             SKIP THIS STEP IF THE RC(AB|CD) PASSES SCREENING CONDITION
              IA1 = (IJ-1)*MAXCD*NTUVABCD + MAXCD*(IADR1-1)
              IA2 = (IJ-1)*MAXCD*NTUVABCD + MAXCD*(IADR2-1)
              IA3 = (IJ-1)*MAXCD*NTUVABCD + MAXCD*(IADR3-1)
              IF(IADR3.NE.0) THEN
                IF(IA1+IA2+IA3.EQ.0) GOTO 403
              ELSE
                IF(IA1+IA2.EQ.0) GOTO 403
              ENDIF
C
C             PRE-FACTORS FOR THE UPCOMING CONTRACTION
              DO N=1,MAXN
                T1 = RCTTFL(IA1+IMAP(N))*0.5D0/APH(IMAP(N))
                T2 = RCTTFL(IA2+IMAP(N))*PQ(IMAP(N),JCMP)
                IF(I3.GE.0.AND.J3.GE.0.AND.K3.GE.0) THEN
                  T3 = RCTTFL(IA3+IMAP(N))*RTP
                ELSE
                  T3 = 0.0D0
                ENDIF
                T(N) = T1-T2+T3
              ENDDO
C
            ENDIF
C
C           CONTRIBUTIONS TO Re{G(AB|--)} FROM EACH Re{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABR11(IAB,ICMP).EQ.0) GOTO 415
            IF(ISYM.EQ.2.AND.IABR11(IAB,ICMP).EQ.0) GOTO 415
            IF(ICDR11(ICD,JCMP).EQ.1) THEN
              DO N=1,MAXN
                GABR11(N,IAB) = GABR11(N,IAB)
     &                            - DREAL(ECD11(IMAP(N),ICD,JCMP))*T(N)
              ENDDO
            ENDIF
415         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|--)} FROM EACH Im{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABI11(IAB,ICMP).EQ.0) GOTO 416
            IF(ISYM.EQ.2.AND.IABI11(IAB,ICMP).EQ.0) GOTO 416
            IF(ICDI11(ICD,JCMP).EQ.1) THEN
              DO N=1,MAXN
                GABI11(N,IAB) = GABI11(N,IAB)
     &                            - DIMAG(ECD11(IMAP(N),ICD,JCMP))*T(N)
              ENDDO
            ENDIF
416         CONTINUE
C
C           CONTRIBUTIONS TO Re{G(AB|--)} FROM EACH Re{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABR21(IAB,ICMP).EQ.0) GOTO 417
            IF(ISYM.EQ.2.AND.IABR21(IAB,ICMP).EQ.0) GOTO 417
            IF(ICDR21(ICD,JCMP).EQ.1) THEN
              DO N=1,MAXN
                GABR21(N,IAB) = GABR21(N,IAB)
     &                            - DREAL(ECD21(IMAP(N),ICD,JCMP))*T(N)
              ENDDO
            ENDIF
417         CONTINUE
C
C           CONTRIBUTIONS TO Im{G(AB|--)} FROM EACH Im{E(CD|--)} ADDRESS
            IF(ISYM.EQ.1.AND.IABI21(IAB,ICMP).EQ.0) GOTO 418
            IF(ISYM.EQ.2.AND.IABI21(IAB,ICMP).EQ.0) GOTO 418
            IF(ICDI21(ICD,JCMP).EQ.1) THEN
              DO N=1,MAXN
                GABI21(N,IAB) = GABI21(N,IAB)
     &                            - DIMAG(ECD21(IMAP(N),ICD,JCMP))*T(N)
              ENDDO
            ENDIF
418         CONTINUE
C
C           SKIP POINT FOR E(CD) SCREENING
403         CONTINUE
C
C           END LOOP OVER CARTESIAN INDEX JCMP FOR CENTRE CD
            ENDDO
C
C           SKIP POINT FOR GAUNT INTERACTION ONLY
480         CONTINUE
C
C         END LOOP OVER E(CD|  ) FINITE EXPANSION ADDRESSES
          ENDDO
C
C       SKIP POINT FOR E(AB|  ) SCREENING
401     CONTINUE
C
C     END LOOP OVER E(AB|  ) FINITE EXPANSION ADDRESSES
      ENDDO
C
C     TIME AT END OF FIRST CONTRACTION FOR THIS ICMP INDEX
      CALL CPU_TIME(T1F)
      TBC1 = TBC1+T1F-T1I
C
C**********************************************************************C
C     PERFORM SECOND CONTRACTION: ( -| -) = E(AB| -)*G(AB| -).         C
C     THIS YIELDS A FULL BATCH OF TWO-ELECTRON INTEGRALS (16 PERM'NS). C
C**********************************************************************C
C
C     CALCULATE PHASES FOR BASIS FUNCTION OVERLAP COMBINATIONS
      PAB =-ISIGN(1,KQN(1)*KQN(2))*(-1)**((MQN(1)-MQN(2))/2)
      PCD =-ISIGN(1,KQN(3)*KQN(4))*(-1)**((MQN(3)-MQN(4))/2)
C
      PABCD = PAB*PCD
C
C     1ST SET: ( 1) = (--|--)   ( 4) = (--|++)
C              (16) = (++|++)   (13) = (++|--)
C
C     RESET CONTRACTION STORAGE LISTS
      DO N=1,MAXN
        QR1(N) = 0.0D0
        QI1(N) = 0.0D0
        QR2(N) = 0.0D0
        QI2(N) = 0.0D0
      ENDDO
C
C     RAW CONTRACTION (--|--) = E(AB|--)*(Re{G(AB|--)} + i*Im{G(AB|--)})
      DO IAB=1,NTUVAB
        IF(IABR11(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB11(IJ,IAB,ICMP))*GABR11(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB11(IJ,IAB,ICMP))*GABI11(N,IAB)
          ENDDO
        ENDIF
        IF(IABI11(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB11(IJ,IAB,ICMP))*GABR11(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB11(IJ,IAB,ICMP))*GABI11(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     ADD THIS ICMP TERM TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N,1 ) = RR(N,1 ) +     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N,4 ) = RR(N,4 ) + PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
      ENDDO
C
      IF(ISYM.EQ.1) GOTO 501
      IF(ISYM.EQ.2) GOTO 501
C
C     2ND SET: ( 3) = (--|+-)   ( 2) = (--|-+)
C              (14) = (++|-+)   (15) = (++|+-)
C
C     RESET CONTRACTION STORAGE LISTS
      DO N=1,MAXN
        QR1(N) = 0.0D0
        QI1(N) = 0.0D0
        QR2(N) = 0.0D0
        QI2(N) = 0.0D0
      ENDDO
C
C     RAW CONTRACTION (--|+-) = E(AB|--)*(Re{G(AB|+-)} + i*Im{G(AB|+-)})
      DO IAB=1,NTUVAB
        IF(IABR11(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB11(IJ,IAB,ICMP))*GABR21(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB11(IJ,IAB,ICMP))*GABI21(N,IAB)
          ENDDO
        ENDIF
        IF(IABI11(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB11(IJ,IAB,ICMP))*GABR21(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB11(IJ,IAB,ICMP))*GABI21(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     ADD THIS ICMP TERM TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N,3 ) = RR(N,3 ) +     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N,2 ) = RR(N,2 ) - PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
      ENDDO
C
C     3RD SET: ( 9) = (+-|--)   (12) = (+-|++)
C              ( 8) = (-+|++)   ( 5) = (-+|--)
C
C     RESET CONTRACTION STORAGE LISTS
      DO N=1,MAXN
        QR1(N) = 0.0D0
        QI1(N) = 0.0D0
        QR2(N) = 0.0D0
        QI2(N) = 0.0D0
      ENDDO
C
C     RAW CONTRACTION (+-|--) = E(AB|+-)*(Re{G(AB|--)} + i*Im{G(AB|--)})
      DO IAB=1,NTUVAB
        IF(IABR21(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB21(IJ,IAB,ICMP))*GABR11(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB21(IJ,IAB,ICMP))*GABI11(N,IAB)
          ENDDO
        ENDIF
        IF(IABI21(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB21(IJ,IAB,ICMP))*GABR11(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB21(IJ,IAB,ICMP))*GABI11(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     ADD THIS ICMP TERM TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N,9 ) = RR(N,9 ) +     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N,12) = RR(N,12) + PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
      ENDDO
C
501   CONTINUE
C
C     4TH SET: (11) = (+-|+-)   (10) = (+-|-+)
C              ( 6) = (-+|-+)   ( 7) = (-+|+-)
C
C     RESET CONTRACTION STORAGE LISTS
      DO N=1,MAXN
        QR1(N) = 0.0D0
        QI1(N) = 0.0D0
        QR2(N) = 0.0D0
        QI2(N) = 0.0D0
      ENDDO
C
C     RAW CONTRACTION (+-|+-) = E(AB|+-)*(Re{G(AB|+-)} + i*Im{G(AB|+-)})
      DO IAB=1,NTUVAB
        IF(IABR21(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QR1(N) = QR1(N) + DREAL(EAB21(IJ,IAB,ICMP))*GABR21(N,IAB)
            QI2(N) = QI2(N) + DREAL(EAB21(IJ,IAB,ICMP))*GABI21(N,IAB)
          ENDDO
        ENDIF
        IF(IABI21(IAB,ICMP).NE.0) THEN
          DO N=1,MAXN
            QI1(N) = QI1(N) + DIMAG(EAB21(IJ,IAB,ICMP))*GABR21(N,IAB)
            QR2(N) = QR2(N) - DIMAG(EAB21(IJ,IAB,ICMP))*GABI21(N,IAB)
          ENDDO
        ENDIF
      ENDDO
C
C     ADD THIS ICMP TERM TO RAW CONTRACTION
      DO N=1,MAXN
        RR(N,11) = RR(N,11) +     DCMPLX(QR1(N)+QR2(N),QI1(N)+QI2(N))
        RR(N,10) = RR(N,10) - PCD*DCMPLX(QR1(N)-QR2(N),QI1(N)-QI2(N))
      ENDDO
C
C     TIME AT END OF SECOND CONTRACTION FOR THIS ICMP INDEX
      CALL CPU_TIME(T2F)
      TBC2 = TBC2+T2F-T1F
C
C     END LOOP OVER CARTESIAN INDICES {IX,IY,IZ}
6000  CONTINUE
C
C     HALF OF THE RR ARRAY CAN BE GENERATED WITH PHASE RELATIONS
      DO N=1,MAXN
        RR(N,16) = PABCD*DCONJG(RR(N,1 ))
        RR(N,13) = PABCD*DCONJG(RR(N,4 ))
        RR(N,14) =-PABCD*DCONJG(RR(N,3 ))
        RR(N,15) =-PABCD*DCONJG(RR(N,2 ))
        RR(N,8 ) =-PABCD*DCONJG(RR(N,9 ))
        RR(N,5 ) =-PABCD*DCONJG(RR(N,12))
        RR(N,6 ) = PABCD*DCONJG(RR(N,11))
        RR(N,7 ) = PABCD*DCONJG(RR(N,10))
      ENDDO     
C
C**********************************************************************C
C     BREIT INTEGRAL BATCH NOW FULLY CONSTRUCTED                       C
C**********************************************************************C
C
C     CALCULATE THE R-INTEGRAL NORMALISATION FACTOR
      M = 0
      DO KBAS=1,NBAS(3)
        DO LBAS=1,NBAS(4)
          M = M+1
          IF(ISCR(M).EQ.1) THEN
            EKL = EXL(KBAS,3)+EXL(LBAS,4)
            EMX = DSQRT(EIJ+EKL)*EIJ*EKL
            PRE(M) = 2.0D0*PI52/EMX
          ENDIF
        ENDDO
      ENDDO
C
C     INCLUDE THE OUTSIDE FACTOR OF (-1/2) AND MOVE TO FULL ARRAY
      DO N=1,MAXN
        DO ITG=1,16
          RR(N,ITG) =-0.5D0*PRE(IMAP(N))*RR(N,ITG)
        ENDDO
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
      FUNCTION NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        NN    NN  CCCCCC  NN    NN TTTTTTTT RRRRRRR   SSSSSS          C
C        NNN   NN CC    CC NNN   NN    TT    RR    RR SS    SS         C
C        NNNN  NN CC       NNNN  NN    TT    RR    RR SS               C
C        NN NN NN CC       NN NN NN    TT    RR    RR  SSSSSS          C
C        NN  NNNN CC       NN  NNNN    TT    RRRRRRR        SS         C
C        NN   NNN CC    CC NN   NNN    TT    RR    RR SS    SS         C
C        NN    NN  CCCCCC  NN    NN    TT    RR    RR  SSSSSS          C
C                                                                      C
C -------------------------------------------------------------------- C
C  NCNTRS RETURNS NUMBER OF UNIQUE NUCLEAR CENTRES FROM INPUT VALUES.  C
C**********************************************************************C
C
C     STORE ICNTA IN FIRST PLACE
      NCNTRS = 1
C
C     CHECK ICNTB AGAINST STORED VALUES
      IF(ICNTB.NE.ICNTA) THEN
        NCNTRS = NCNTRS + 1
      ENDIF
C
C     CHECK ICNTC AGAINST STORED VALUES
      IF(ICNTC.NE.ICNTA.AND.ICNTC.NE.ICNTB) THEN
        NCNTRS = NCNTRS + 1
      ENDIF
C
C     CHECK ICNTD AGAINST STORED VALUES
      IF(ICNTD.NE.ICNTA.AND.ICNTD.NE.ICNTB.AND.ICNTD.NE.ICNTC) THEN
        NCNTRS = NCNTRS + 1
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION ZPROJ(XYZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            ZZZZZZZZ PPPPPPP  RRRRRRR   OOOOOO     JJJJJJ             C
C                 ZZ  PP    PP RR    RR OO    OO       JJ              C
C                ZZ   PP    PP RR    RR OO    OO       JJ              C
C               ZZ    PP    PP RR    RR OO    OO       JJ              C
C              ZZ     PPPPPPP  RRRRRRR  OO    OO       JJ              C
C             ZZ      PP       RR    RR OO    OO JJ    JJ              C
C            ZZZZZZZZ PP       RR    RR  OOOOOO   JJJJJJ               C
C                                                                      C
C -------------------------------------------------------------------- C
C  ZPROJ RETURNS THE SUM OF ABSOLUTE DIFFERENCES BETWEEN X AND Y       C
C  LOCATIONS OF FOUR SUPPLIED COORDINATES.                             C
C**********************************************************************C
C
       DIMENSION XYZ(3,4)
C
       ZPROJ = 0.0D0
C
       DO IX=1,2
         DO M=1,3
           DO N=M+1,4
             ZPROJ = ZPROJ + DABS(XYZ(IX,M)-XYZ(IX,N))
           ENDDO
         ENDDO
       ENDDO
C
       RETURN
       END
C
C
      FUNCTION LVAL(KQN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                  LL      VV    VV    AA    LL                        C
C                  LL      VV    VV   AAAA   LL                        C
C                  LL      VV    VV  AA  AA  LL                        C
C                  LL      VV    VV AA    AA LL                        C
C                  LL       VV  VV  AAAAAAAA LL                        C
C                  LL        VVVV   AA    AA LL                        C
C                  LLLLLLLL   VV    AA    AA LLLLLLLL                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  LVAL RETURNS THE LQN FROM A KQN VALUE.                              C
C**********************************************************************C
C
      IF(KQN.LT.0) THEN
        LVAL =-KQN-1
      ELSE
        LVAL = KQN
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION LLAB(LQN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 LL       LL          AA    BBBBBBB                   C
C                 LL       LL         AAAA   BB    BB                  C
C                 LL       LL        AA  AA  BB    BB                  C
C                 LL       LL       AA    AA BBBBBBB                   C
C                 LL       LL       AAAAAAAA BB    BB                  C
C                 LL       LL       AA    AA BB    BB                  C
C                 LLLLLLLL LLLLLLLL AA    AA BBBBBBB                   C
C                                                                      C
C -------------------------------------------------------------------- C
C  LLAB RETURNS THE CONVENTIONAL ATOMIC ORBITAL TYPE LABEL FOR THE     C
C  NON-RELATIVISTIC QUANTUM NUMBER LQN OF LENGTH 1.                    C
C**********************************************************************C
C
      CHARACTER*1 LLAB
C
      IF(LQN.EQ.0) THEN
        LLAB = 's'
      ELSEIF(LQN.EQ.1) THEN
        LLAB = 'p'
      ELSEIF(LQN.EQ.2) THEN
        LLAB = 'd'
      ELSEIF(LQN.EQ.3) THEN
        LLAB = 'f'
      ELSEIF(LQN.EQ.4) THEN
        LLAB = 'g'
      ELSEIF(LQN.EQ.5) THEN
        LLAB = 'h'
      ELSEIF(LQN.EQ.6) THEN
        LLAB = 'i'
      ELSEIF(LQN.EQ.7) THEN
        LLAB = 'j'
      ELSEIF(LQN.EQ.8) THEN
        LLAB = 'k'
      ELSEIF(LQN.EQ.9) THEN
        LLAB = 'l'
      ELSEIF(LQN.EQ.10) THEN
        LLAB = 'm'
      ELSEIF(LQN.EQ.11) THEN
        LLAB = 'n'
      ELSEIF(LQN.EQ.12) THEN
        LLAB = 'o'
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION KLAB(KQN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 KK    KK LL          AA    BBBBBBB                   C
C                 KK   KK  LL         AAAA   BB    BB                  C
C                 KK  KK   LL        AA  AA  BB    BB                  C
C                 KKKKK    LL       AA    AA BBBBBBB                   C
C                 KK  KK   LL       AAAAAAAA BB    BB                  C
C                 KK   KK  LL       AA    AA BB    BB                  C
C                 KK    KK LLLLLLLL AA    AA BBBBBBB                   C
C                                                                      C
C -------------------------------------------------------------------- C
C  KLAB RETURNS THE CONVENTIONAL SYMMETRY TYPE LABEL FOR DIRAC         C
C  ORBITAL TYPE KQN AS A STRING OF LENGTH 2.                           C
C**********************************************************************C
C
      CHARACTER*1 LLAB,CHS
      CHARACTER*2 KLAB
C
      LQN = LVAL(KQN)
C
      IF(KQN.EQ.-1) THEN
        CHS = ' '
      ELSEIF(KQN.LT.0) THEN
        CHS = '-'
      ELSE
        CHS = '+'
      ENDIF
C
C     STITCH TOGETHER THE LQN TITLE AND THE PARITY LABEL
      WRITE(KLAB,'(A,A)') LLAB(LQN),CHS
C
      RETURN
      END
C
C
      FUNCTION MLAB(MQN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                MM       MM LL          AA    BBBBBBB                 C
C                MMM     MMM LL         AAAA   BB    BB                C
C                MMMM   MMMM LL        AA  AA  BB    BB                C
C                MM MM MM MM LL       AA    AA BBBBBBB                 C
C                MM  MMM  MM LL       AAAAAAAA BB    BB                C
C                MM   M   MM LL       AA    AA BB    BB                C
C                MM       MM LLLLLLLL AA    AA BBBBBBB                 C
C                                                                      C
C -------------------------------------------------------------------- C
C  MLAB RETURNS THE CONVENTIONAL SYMMETRY TYPE LABEL FOR DIRAC         C
C  MAGNETIC NUMBER MQN AS A STRING OF LENGTH 5.                        C
C**********************************************************************C
C
       CHARACTER*1 CHS
       CHARACTER*5 MLAB
C
       IF(MOD(MQN,2).EQ.0) THEN
         CHS = '+'
       ELSE
         CHS = '-'
       ENDIF
C
C      ODD INTEGER FOR NUMERATOR
       MAG = MQN-MOD(MQN+1,2)
C
C      STITCH TOGETHER THE MQN TITLE
       IF(MAG.LT.10) THEN
         WRITE(MLAB,'(A,I1,A)') CHS,MAG,'/2 '
       ELSE
         WRITE(MLAB,'(A,I2,A)') CHS,MAG,'/2'
       ENDIF
C
       RETURN
       END
C
C
C**********************************************************************C
C ==================================================================== C
C   [8] MBPTN: CORRELATION ENERGY CALCULATION ROUTINES.                C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] MBPTN: MAIN ROUTINE FOR MANY-BODY DIAGRAMMATIC P.T.            C
C   [B] MBPT1: ZERO- AND FIRST-ORDER ENERGY ANALYSIS.                  C
C   [C] MBPT2: SECOND-ORDER CORRELATION ENERGY ANALYSIS.               C
C   [D] MBPTCB: COULOMB + GAUNT/BREIT SECOND ORDER REDUCTION.          C
C**********************************************************************C
C
C
      SUBROUTINE MBPTN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            MM       MM BBBBBBB  PPPPPPP TTTTTTTT NN    NN            C
C            MMM     MMM BB    BB PP    PP   TT    NNN   NN            C
C            MMMM   MMMM BB    BB PP    PP   TT    NNNN  NN            C
C            MM MM MM MM BBBBBBB  PP    PP   TT    NN NN NN            C
C            MM  MMM  MM BB    BB PPPPPPP    TT    NN  NNNN            C
C            MM   M   MM BB    BB PP         TT    NN   NNN            C
C            MM       MM BBBBBBB  PP         TT    NN    NN            C
C                                                                      C
C                       CONTROLLING ROUTINE FOR                        C
C                          ** B E R T H A **                           C
C -------------------------------------------------------------------- C
C  MBPTN PERFORMS NTH-ORDER MANY-BODY DIAGRAMMATIC EVALUATION ON A     C
C  CONVERGED MOLECULAR HARTREE-FOCK SOLUTION SPACE.                    C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 G2INT
      CHARACTER*6 IOSAV
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/TMPT/EMTY
C
C     PRINT A BIG SECTION HEADER
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) REPEAT(' ',21),'MANY-BODY PERTURBATION THEORY'
      WRITE(7, *) REPEAT(' ',21),'MANY-BODY PERTURBATION THEORY'
      WRITE(6, *) REPEAT('#',72)
      WRITE(7, *) REPEAT('#',72)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C     RECORD TIME
      CALL CPU_TIME(TDUM)
C
C     EQ-COEFFICIENTS IN COMMON ARRAYS
      IF(EQFILE) THEN
        CALL EQSAVE
      ENDIF
C
C     GENERATE DENSITY MATRIX -- FOR ONE-BODY TERMS
      CALL DENSTY
C
C     RESTRICT SET OF OCCUPIED AND VIRTUAL STATES
      MINO = 1
      NUMO = NOCC-MINO+1
      MINV = NOCC+1
      NUMV = 100
C     NUMV = NVRT-200
C
C     LABEL THE TWO-BODY INTERACTION ('COULM', 'GAUNT' OR 'BREIT')
      G2INT = 'COULM'
c     G2INT = 'BREIT'
C
C     LABEL THE IO SAVE/WRITE OPTION
C     IOSAV = 'WRTOUT'
C     IOSAV = 'NOTHNG'
      IOSAV = 'READIN'
C     IOSAV = 'SILENT'
C
C     SUMMARY OF INPUT CHOICES
      WRITE(6, *) 'First occupied orbital:     ',MINO
      WRITE(7, *) 'First occupied orbital:     ',MINO
      WRITE(6, *) 'Number of occupied orbitals:',NUMO
      WRITE(7, *) 'Number of occupied orbitals:',NUMO
      WRITE(6, *) 'First virtual orbital:      ',MINV
      WRITE(7, *) 'First virtual orbital:      ',MINV
      WRITE(6, *) 'Number of virtual orbitals: ',NUMV
      WRITE(7, *) 'Number of virtual orbitals: ',NUMV
      IF(G2INT.EQ.'COULM') THEN
        WRITE(6, *) '2-body interaction: Coulomb.'
        WRITE(7, *) '2-body interaction: Coulomb.'
      ELSEIF(G2INT.EQ.'GAUNT') THEN
        WRITE(6, *) '2-body interaction: Gaunt.'
        WRITE(7, *) '2-body interaction: Gaunt.'
      ELSEIF(G2INT.EQ.'BREIT') THEN
        WRITE(6, *) '2-body interaction: Breit.'
        WRITE(7, *) '2-body interaction: Breit.'
      ENDIF
C
C     CALL FIRST-ORDER MBPT ROUTINE (SHOULD CONFIRM SCF RESULTS)
c     CALL MBPT1(MINO,NUMO,G2INT,IOSAV)
C
C     CALL SECOND-ORDER MBPT ROUTINE
      CALL MBPT2(MINO,NUMO,MINV,NUMV,G2INT,IOSAV)
C
C     COULOMB + GAUNT/BREIT E1 AND E2 ANALYSIS
c     CALL MBPTNPR(MINO,NUMO,MINV,NUMV,'COULM','BREIT')
C
C     TOTAL TIME TAKEN
      CALL CPU_TIME(TMPT)
      TTOT = TTOT-TDUM
C
      RETURN
      END
C
C
      SUBROUTINE MBPT1(MINO,NUMO,G2INT,IOSAV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C              MM       MM BBBBBBB  PPPPPPP TTTTTTTT 11                C
C              MMM     MMM BB    BB PP    PP   TT   111                C
C              MMMM   MMMM BB    BB PP    PP   TT    11                C
C              MM MM MM MM BBBBBBB  PP    PP   TT    11                C
C              MM  MMM  MM BB    BB PPPPPPP    TT    11                C
C              MM   M   MM BB    BB PP         TT    11                C
C              MM       MM BBBBBBB  PP         TT   1111               C
C                                                                      C
C -------------------------------------------------------------------- C
C  MBPT1 EVALUATES ZERO- AND FIRST-ORDER ENERGIES FOR ALL OCCUPIED     C
C  SOLUTIONS TO A CONVERGED HARTREE-FOCK PROBLEM.                      C
C -------------------------------------------------------------------- C
C INPUT:                                                               C
C  ▶ MINO  - LOWEST OCCUPIED STATE TO ACCOUNT FOR. (FULL: 1)           C
C  ▶ NUMO  - NUMBER OF OCCUPIED STATES TO ACCOUNT FOR. (FULL: NOCC)    C
C  ▶ G2INT - NAME OF TWO-BODY OPERATOR ('COULM', 'GAUNT' OR 'BREIT').  C
C  ▶ IOSAV - I/O OPTION (READ RESULTS, WRITE RESULTS, DO NOTHING).     C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL FILETHERE,BGAUNT
C
      CHARACTER*5  G2INT
      CHARACTER*6  IOSAV
      CHARACTER*16 HMS
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*60 CADBOUT,EAB1OUT
      CHARACTER*80 TITLE
C
      INTEGER*16 INEST,NNEST,NCADB
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4),LQN(4),JQN(4)
      DIMENSION INDEX(MCT,-(MEL+1):MEL,MKP),ITN(2)
      DIMENSION EAB1(NUMO,NUMO,6),EA1(NUMO,6)
C
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 B1(MBS*NUMO,8),B2(MBS*NUMO,8)
      COMPLEX*16 DB(MB2,NUMO*NUMO,4)
      COMPLEX*16 ADB1(MBS,NUMO*NUMO*NUMO,2),ADB2(MBS,NUMO*NUMO*NUMO,2)
      COMPLEX*16 CADB((NUMO+1)*NUMO*NUMO*NUMO/2)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM),HNUC(MDM,MDM),
     &           HKIN(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           BDIR(MDM,MDM),BXCH(MDM,MDM),VANM(MDM,MDM),
     &           VSLF(MDM,MDM),VUEH(MDM,MDM),VWKR(MDM,MDM),
     &           VKSB(MDM,MDM),QDIR(MDM,MDM),QXCH(MDM,MDM),
     &           WDIR(MDM,MDM),WXCH(MDM,MDM),CPLE(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/E0LL/E0LLFL(MFL,4),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,4),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EILS/EILSFL(MFL,12),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EISL/EISLFL(MFL,12),IADISL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS,IABSL,ICDSL
      COMMON/IRCM/IEAB,IECD,IRIJ(MBS,MBS)
      COMMON/ISCR/IMTX(MB2,11),ISCR(MB2),IMAP(MB2),IBCH,ITOG,MAXN
      COMMON/MTRX/FOCK,OVLP,HNUC,HKIN,GDIR,GXCH,BDIR,BXCH,VANM,VSLF,
     &            VUEH,VWKR,VKSB,QDIR,QXCH,WDIR,WXCH,CPLE
C
C     TURN OFF RC(AB|CD) LOCAL FILE PROCESS
      RCFILE = .FALSE.
C
C     WARNINGS BASED ON INVALID HMLT VS. G2INT COMBINATIONS
1     FORMAT(1X,'In MBPT1: HMLT = ',A,' but G2INT = ',A,'.')
      IF(HMLT.EQ.'BARE') THEN
        WRITE(6,1) HMLT,G2INT
        WRITE(7,1) HMLT,G2INT
      ELSEIF(HMLT.EQ.'NORL'.OR.HMLT.EQ.'DHFR'.OR.HMLT.EQ.'DHFP') THEN
        IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') THEN
          WRITE(6,1) HMLT,G2INT
          WRITE(7,1) HMLT,G2INT
        ENDIF
      ENDIF
C
C     INTEGRAL SKIPPING ON MOLECULAR GROUP SYMMETRY CLASS BASIS
      IF(SHAPE.EQ.'ATOMIC') THEN
        ISYM = 2
      ELSEIF(SHAPE.EQ.'DIATOM'.OR.SHAPE.EQ.'LINEAR') THEN
        ISYM = 1
      ELSE
        ISYM = 0
      ENDIF
C
C     RELATIVISTIC COMPONENT LABELS FOR LOOPING
      IF(HMLT.EQ.'NORL') THEN
        IF(G2INT.EQ.'COULM') THEN
          ITSTRT = 1
          ITSTOP = 1
          ITSKIP = 1
        ENDIF
      ELSE
        IF(G2INT.EQ.'COULM') THEN
          ITSTRT = 4
          ITSTOP = 1
          ITSKIP =-3
        ELSEIF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') THEN
          ITSTRT = 2
          ITSTOP = 3
          ITSKIP = 1        
        ENDIF
      ENDIF
C
C     GAUNT REPLACEMENT TOGGLE
      IF(G2INT.EQ.'GAUNT') THEN
        BGAUNT = .TRUE.
      ELSE
        BGAUNT = .FALSE.
      ENDIF
C
C     INITIALISE TIME COUNTERS
      T1EL = 0.0D0
      TERI = 0.0D0
      TCN1 = 0.0D0
      TCN2 = 0.0D0
      TCN3 = 0.0D0
      TCN4 = 0.0D0
      TSUM = 0.0D0
C
C     START TIME COUNTER      
      CALL CPU_TIME(TBEG)
C
C     INITIALISE ENERGY COUNTERS
      DO N=1,6
        DO IOCCB=1,NUMO
          EA1(IOCCB,N) = 0.0D0
          DO IOCCA=1,NUMO
            EAB1(IOCCA,IOCCB,N) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     INTIALISE THE ARRAY FOR (AC|BD) VALUES
      M = 0
      DO IOCCC=1,NUMO
        DO IOCCA=1,NUMO
          DO IOCCB=1,IOCCA
            DO IOCCD=1,NUMO
              M = M+1
              CADB(M) = DCMPLX(0.0D0,0.0D0)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     TRY TO READ CABD ARRAY FROM EXISTING FILE
      IF(IOSAV.EQ.'READIN') THEN
        CADBOUT = TRIM(OUTFL)//'_CADB_'//G2INT//'.dat'
        INQUIRE(FILE=TRIM(CADBOUT),EXIST=FILETHERE)
        IF(FILETHERE) THEN
          OPEN(UNIT=8,FILE=TRIM(CADBOUT),STATUS='OLD')
          REWIND(UNIT=8)
C         READ MINO AND NUMO
          READ(8, *) MINOTMP,NUMOTMP
C         IF DIFFERENT TO INPUT, WRITE TO LOG AND UPDATE
          IF(MINOTMP.NE.MINO) THEN
            WRITE(6, *) 'In MBPT1: input MINO = ',MINO
            WRITE(6, *) '           file MINO = ',MINOTMP
          ENDIF
          IF(NUMOTMP.NE.NUMO) THEN
            WRITE(6, *) 'In MBPT1: input NUMO = ',NUMO
            WRITE(6, *) '           file NUMO = ',NUMOTMP
          ENDIF
          MINO = MINOTMP
          NUMO = NUMOTMP
C         LOOP OVER OCCUPIED STATES IOCCA AND VIRTUAL STATES IVRTR
          DO IOCCC=1,NUMO
            DO IOCCA=1,NUMO
C             LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
              DO IOCCD=1,NUMO
                DO IOCCB=1,IOCCA
C                 MAIN LIST ADDRESS FOR THIS IOCCA,IVRTR,IOCCB,IVRTS
                  MCADB = (IOCCA-1)*NUMO*IOCCA*NUMO/2
     &                  + (IOCCC-1)*NUMO*IOCCA + (IOCCB-1)*NUMO + IOCCD
                  READ(8, *) CADB(MCADB)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CLOSE(UNIT=8)
C         SKIP CALCULATION PHASE
          GOTO 80
        ELSE
C       IF THE FILE DOESN'T YET EXIST, CALCULATE AND SAVE IT
          WRITE(6, *) 'In MBPT1: tried to read CADB but failed.'
          WRITE(7, *) 'In MBPT1: tried to read CADB but failed.'
          IOSAV = 'WRTOUT'
        ENDIF
      ENDIF
C
C**********************************************************************C
C     ORDERED INDEX OF (ICNT,KQN,MQN) COMBINATIONS                     C
C**********************************************************************C
C
      ICOUNT = 0
C
C     LOOP OVER NUCLEAR CENTRES
      DO ICT=1,NCNT
C
C       LOOP OVER KAPPA VALUES FOR THIS NUCLEAR CENTRE
        DO KN=1,NKAP(ICT)
C
C         IMPORT KAPPA, MAXIMUM MQN
          KAPPA = KAPA(KN,ICT)
          MJMAX = 2*IABS(KAPPA)-1
C
C         LOOP OVER MQN VALUES AND RECORD INDEX
          DO MJ=1,MJMAX,2
            ICOUNT              = ICOUNT+1
            INDEX(ICT,KAPPA,MJ) = ICOUNT
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     ORDERED COUNTER FOR NESTED LOOP PROGRESS (ASSUMES ERI/BII FAST)  C
C**********************************************************************C
C
      NNEST = 0
C
C     BLOCKS A AND B
      DO 100 ICNTA=1,NCNT
      DO 100 ICNTB=1,ICNTA
      DO 200 KA=1,NKAP(ICNTA)
      KQN(1) = KAPA(KA,ICNTA)
      LQN(1) = LVAL(KQN(1))
      JQN(1) = 2*IABS(KQN(1))-1
      NBAS(1) = NFNC(LQN(1),ICNTA)
      DO 200 KB=1,NKAP(ICNTB)
      KQN(2) = KAPA(KB,ICNTB)
      LQN(2) = LVAL(KQN(2))
      JQN(2) = 2*IABS(KQN(2))-1
      NBAS(2) = NFNC(LQN(2),ICNTB)
      DO 300 MA=1,IABS(KAPA(KA,ICNTA))
      MQN(1) = 2*MA-1
      IQ1 = INDEX(ICNTA,KQN(1),2*MA-1)
      DO 300 MB=1,IABS(KAPA(KB,ICNTB))
      MQN(2) = 2*MB-1
      IQ2 = INDEX(ICNTB,KQN(2),2*MB-1)
C       REASONS TO SKIP
        IF(G2INT.EQ.'COULM') THEN
          IF(IQ1.LT.IQ2) GOTO 301
        ENDIF
C       COMPONENT IT1
        DO 400 IT1=ITSTRT,ITSTOP,ITSKIP
C       BLOCKS C AND D
        DO 500 ICNTC=1,NCNT
        DO 500 ICNTD=1,NCNT
        MCNT = NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
        DO 600 KC=1,NKAP(ICNTC)
        KQN(3) = KAPA(KC,ICNTC)
        LQN(3) = LVAL(KQN(3))
        JQN(3) = 2*IABS(KQN(3))-1
        NBAS(3) = NFNC(LQN(3),ICNTC)
        DO 600 KD=1,NKAP(ICNTD)
        KQN(4) = KAPA(KD,ICNTD)
        LQN(4) = LVAL(KQN(4))
        JQN(4) = 2*IABS(KQN(4))-1
        NBAS(4) = NFNC(LQN(4),ICNTD)
C       SELECTION RULES BASED ON KQN
        IF(MCNT.EQ.1) THEN
          IPARAB = MOD(LQN(1)+LQN(2),2)
          IPARCD = MOD(LQN(3)+LQN(4),2)
          IF(IPARAB.NE.IPARCD) THEN
            GOTO 601
          ENDIF
          NUI = MAX0(IABS(JQN(1)-JQN(2))/2,IABS(JQN(3)-JQN(4))/2)
          NUF = MIN0(    (JQN(1)+JQN(2))/2,    (JQN(3)+JQN(4))/2)
          IF(NUI.GT.NUF) THEN
            GOTO 601
          ENDIF
          ISELK = 0
          DO NU=NUI,NUF
            IPARAB = MOD(LQN(1)+LQN(2)+NU,2)
            IPARCD = MOD(LQN(3)+LQN(4)+NU,2)
            IF(IPARAB.EQ.1.AND.IPARCD.EQ.1.AND.NU.NE.0) ISELK = 1
            IF(IPARAB.EQ.0.AND.IPARCD.EQ.0) ISELK = 1
          ENDDO
          IF(ISELK.EQ.0) GOTO 601
        ENDIF
        DO 700 MC=1,IABS(KAPA(KC,ICNTC))
        MQN(3) = 2*MC-1
        IQ3 = INDEX(ICNTC,KQN(3),2*MC-1)
        DO 700 MD=1,IABS(KAPA(KD,ICNTD))
        MQN(4) = 2*MD-1
        IQ4 = INDEX(ICNTD,KQN(4),2*MD-1)
C         REASONS TO SKIP
          IF(G2INT.EQ.'COULM') THEN
            IF(IQ3.LT.IQ4) GOTO 701
          ENDIF
C         SELECTION RULES BASED ON MQN
          IF(MCNT.EQ.1) THEN
            ISELM = 0
            DO ISGN1=1,2
              DO ISGN2=1,2
                DO ISGN3=1,2
                  DO ISGN4=1,2
                    MMJA = MQN(1)*((-1)**ISGN1)
                    MMJB = MQN(2)*((-1)**ISGN2)
                    MMJC = MQN(3)*((-1)**ISGN3)
                    MMJD = MQN(4)*((-1)**ISGN4)
                    IF(MMJA-MMJB.EQ.MMJD-MMJC) ISELM = 1
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            IF(ISELM.EQ.0) GOTO 701
          ENDIF
          IF(ISYM.EQ.1.OR.ISYM.EQ.2) THEN
            ISELM = 0
            IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) ISELM = 1
            IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) ISELM = 1
            IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) ISELM = 1
            IF(ISELM.EQ.0) GOTO 701
          ENDIF
C         COMPONENT IT2
          DO 800 IT2=ITSTRT,ITSTOP,ITSKIP
C           FIRST CONTRACTION
            NNEST = NNEST + NBAS(1)*NBAS(2)*NBAS(3)*NBAS(4)*NUMO
            IF(G2INT.EQ.'COULM'.AND.IQ3.EQ.IQ4) GOTO 901
            IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') GOTO 901
            NNEST = NNEST + NBAS(1)*NBAS(2)*NBAS(3)*NBAS(4)*NUMO
901         CONTINUE
C           SECOND CONTRACTION
            NNEST = NNEST + NBAS(1)*NBAS(2)*NBAS(3)*NUMO*NUMO
            IF(G2INT.EQ.'COULM'.AND.IQ3.EQ.IQ4) GOTO 902
            IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') GOTO 902
            NNEST = NNEST + NBAS(1)*NBAS(2)*NBAS(4)*NUMO*NUMO
902         CONTINUE
800       CONTINUE
701       CONTINUE
700       CONTINUE
601       CONTINUE
600       CONTINUE
500       CONTINUE
C         THIRD CONTRACTION
          NNEST = NNEST + NBAS(1)*NBAS(2)*NUMO*(NUMO+1)*NUMO/2
          IF(G2INT.EQ.'COULM'.AND.IQ1.EQ.IQ2) GOTO 401
          IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') GOTO 401
          NNEST = NNEST + NBAS(1)*NBAS(2)*NUMO*(NUMO+1)*NUMO/2
401       CONTINUE
C         FOURTH CONTRACTION
          NNEST = NNEST + NBAS(1)*NUMO*NUMO*(NUMO+1)*NUMO/2 
          IF(G2INT.EQ.'COULM'.AND.IQ1.EQ.IQ2) GOTO 402
          IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') GOTO 402
          NNEST = NNEST + NBAS(2)*NUMO*NUMO*(NUMO+1)*NUMO/2 
402       CONTINUE
400     CONTINUE
301     CONTINUE
300     CONTINUE
200     CONTINUE
100     CONTINUE
C
C     START CONTRACTION COUNTER
      INEST = 0
C
C**********************************************************************C
C     ONE-BODY ENERGIES (INSIGNIFICANT COMPUTATIONAL COST)             C
C**********************************************************************C
C
C     RECORD TIME AT THE START OF THIS PROCESS
      CALL CPU_TIME(T1)
C
C     CALCULATE ONE-BODY MATRIX REPS
      CALL ONEEL
C
C     LOOP OVER ALL OCCUPIED ORBITAL PAIRS AND CONTRACT
      E1H = 0.0D0
      DO IOCCA=1,NUMO
        DO IOCCB=1,IOCCA
C
C         FOCK MATRIX ADDRESS FOR IOCCA AND IOCCB
          IA = MINO-1+IOCCA+NSKP
          IB = MINO-1+IOCCB+NSKP
C
C         ONE-BODY ENERGY
          EN = 0.0D0
          ET = 0.0D0
          IF(IOCCA.EQ.IOCCB) THEN
            DO J=1,NDIM
              DO I=1,NDIM
                EN = EN + DREAL(HNUC(I,J)*DCONJG(COEF(I,IA))*COEF(J,IA))
                ET = ET + DREAL(HKIN(I,J)*DCONJG(COEF(I,IA))*COEF(J,IA))
              ENDDO
            ENDDO
          ENDIF
          EAB1(IOCCA,IOCCB,1) = EN
          EAB1(IOCCA,IOCCB,2) = ET
          EAB1(IOCCA,IOCCB,3) = EN+ET
          E1H = E1H + EAB1(IOCCA,IOCCB,3)
C
        ENDDO
      ENDDO
C
C     RECORD TIME AT THE END OF THIS PROCESS
      CALL CPU_TIME(T2)
      T1EL = T1EL+T2-T1
C
C**********************************************************************C
C     LOOP OVER ATOMIC CENTRES A AND B (USE INDEX 1000)                C
C**********************************************************************C
C
C     LOOP OVER CENTRE A
      DO 1000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 1000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C**********************************************************************C
C     LOOP OVER KQN SYMMETRY TYPES A AND B (USE INDEX 2000)            C
C**********************************************************************C
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,ICNTA)
        LQN(1) = LVAL(KQN(1))
        JQN(1) = 2*IABS(KQN(1))-1
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
        LQN(2) = LVAL(KQN(2))
        JQN(2) = 2*IABS(KQN(2))-1
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1, NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
        ENDDO
C
C**********************************************************************C
C     LOOP OVER |MQN| PROJECTIONS A AND B (INDEX 3000)                 C
C**********************************************************************C
C
C     LOOP OVER |MQN(A)| VALUES
      DO 3000 MA=1,IABS(KQN(1))
        MQN(1) = 2*MA-1
C
C     LOOP OVER |MQN(B)| VALUES
      DO 3000 MB=1,IABS(KQN(2))
        MQN(2) = 2*MB-1
C
C     INDEX ASSIGNMENT
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
C
C     SKIP CONTRIBUTIONS THAT ARISE BY PERMUTATION OF INTEGRALS
      IF(G2INT.EQ.'COULM') THEN
        IF(IQ1.LT.IQ2) GOTO 3001
      ENDIF
C     IF(IQ1.LT.IQ2) GOTO 3001
C
C     EQ-COEFFICIENT PHASE FACTORS FOR PERMUTATION OF R-INTEGRALS
      IF(G2INT.EQ.'COULM') THEN
        PAB1 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
        PAB2 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      ELSEIF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') THEN
        PAB1 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
        PAB2 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      ENDIF
C
C**********************************************************************C
C     LOOP OVER COMPONENT OVERLAP LABELS A AND B (INDEX 4000)          C
C**********************************************************************C
C
C     COMPONENT LABEL FOR A AND B: TT = LL(1) or SS(4) <- COULOMB
C                                  TT = LS(2) or SL(3) <- BREIT
      DO 4000 IT1=ITSTRT,ITSTOP,ITSKIP
C
C       EQ-COEFFICIENT STARTING ADDRESSES FOR (AB) PAIR
        IF(IT1.EQ.1) THEN
          IABLL = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB)
        ELSEIF(IT1.EQ.2) THEN
          IABLS = IADILS(ICNTA,ICNTB,KA,KB,MA,MB)
        ELSEIF(IT1.EQ.3) THEN
          IABSL = IADISL(ICNTA,ICNTB,KA,KB,MA,MB)
        ELSEIF(IT1.EQ.4) THEN
          IABSS = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB)
        ENDIF
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITN(1) = IT1
C
C       CALCULATE STARTING ADDRESS
        IF(ITN(1).EQ.1) THEN
          NADDA = 0
          NADDB = 0
        ELSEIF(ITN(1).EQ.2) THEN
          NADDA = 0
          NADDB = NSKP
        ELSEIF(ITN(1).EQ.3) THEN
          NADDA = NSKP
          NADDB = 0
        ELSEIF(ITN(1).EQ.4) THEN
          NADDA = NSKP
          NADDB = NSKP
        ENDIF
C
C       FLAG READ-IN OF E0(AB) COEFFICIENTS FOR THIS COMPONENT LABEL
        IEAB = 1
C
C**********************************************************************C
C     FOCK MATRIX STARTING ADDRESSES                                   C
C**********************************************************************C
C
C     FOCK ADDRESS FOR EACH BASIS FUNCTION (WITH SPIN PROJECTION)
      NA1 = LRGE(ICNTA,KA,2*MA-1) + NADDA
      NA2 = LRGE(ICNTA,KA,2*MA  ) + NADDA
C
      NB1 = LRGE(ICNTB,KB,2*MB-1) + NADDB
      NB2 = LRGE(ICNTB,KB,2*MB  ) + NADDB
C
C     CLEAR ARRAY FOR THE COMPLETED CONTRACTION OVER BLOCKS C AND D
      DO MDB=1,NUMO*NUMO
        DO MIJ=1,NBAS(1)*NBAS(2)
          DO IJSPIN=1,4
            DB(MIJ,MDB,IJSPIN) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     LOOP OVER ATOMIC CENTRES C AND D (USE INDEX 5000)                C
C**********************************************************************C
C
C     LOOP OVER CENTRE C
      DO 5000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE C
        XYZ(1,3) = BXYZ(1,ICNTC)
        XYZ(2,3) = BXYZ(2,ICNTC)
        XYZ(3,3) = BXYZ(3,ICNTC)
C
C     LOOP OVER CENTRE D
      DO 5000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE D
        XYZ(1,4) = BXYZ(1,ICNTD)
        XYZ(2,4) = BXYZ(2,ICNTD)
        XYZ(3,4) = BXYZ(3,ICNTD)
C
C     NUMBER OF NUCLEAR CENTRES INVOLVED IN THIS OVERLAP
      MCNT = NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
C
C**********************************************************************C
C     LOOP OVER KQN SYMMETRY TYPES C AND D (USE INDEX 6000)            C
C**********************************************************************C
C
C     LOOP OVER KQN(C) VALUES
      DO 6000 KC=1,NKAP(ICNTC)
C
C       QUANTUM NUMBERS FOR BLOCK C
        KQN(3) = KAPA(KC,ICNTC)
        LQN(3) = LVAL(KQN(3))
        JQN(3) = 2*IABS(KQN(3))-1
C
C       BASIS EXPONENTS FOR BLOCK C
        NBAS(3) = NFNC(LQN(3),ICNTC)
        DO KBAS=1,NBAS(3)
          EXL(KBAS,3) = BEXL(KBAS,LQN(3),ICNTC)
        ENDDO
C
C     LOOP OVER KQN(D) VALUES
      DO 6000 KD=1,NKAP(ICNTD)
C
C       QUANTUM NUMBERS FOR BLOCK D
        KQN(4) = KAPA(KD,ICNTD)
        LQN(4) = LVAL(KQN(4))
        JQN(4) = 2*IABS(KQN(4))-1
C
C       BASIS EXPONENTS FOR BLOCK D
        NBAS(4) = NFNC(LQN(4),ICNTD)
        DO LBAS=1,NBAS(4)
          EXL(LBAS,4) = BEXL(LBAS,LQN(4),ICNTD)
        ENDDO
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON KQN                           C
C**********************************************************************C
C
C     ATOM-CENTRED SELECTION RULES
      IF(MCNT.EQ.1) THEN
C
C       LQN PAIR PARITY (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQN(1)+LQN(2),2)
        IPARCD = MOD(LQN(3)+LQN(4),2)
C
C       LQN PAIR PARITY SELECTION RULE
        IF(IPARAB.NE.IPARCD) THEN
          GOTO 6001
        ENDIF
C
C       JQN TRIANGLE RULE CHECK FOR MULTIPOLE EXPANSION (ATOM-CENTRED)
        NUI = MAX0(IABS(JQN(1)-JQN(2))/2,IABS(JQN(3)-JQN(4))/2)
        NUF = MIN0(    (JQN(1)+JQN(2))/2,    (JQN(3)+JQN(4))/2)
        IF(NUI.GT.NUF) THEN
          GOTO 6001
        ENDIF
C
C       ADDITIONAL LQN SELECTION RULE PARITY ANALYSIS
        ISELK = 0
        DO NU=NUI,NUF
C
C         A AND B: LQN(1)+LQN(2)+NU EVEN OR ODD (0 IF EVEN, 1 IF ODD)
          IPARAB = MOD(LQN(1)+LQN(2)+NU,2)
C
C         C AND D: LQN(3)+LQN(4)+NU EVEN OR ODD (0 IF EVEN, 1 IF ODD)
          IPARCD = MOD(LQN(3)+LQN(4)+NU,2)
C
C         CASE 1: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH ODD (UNLESS NU=0)
          IF(IPARAB.EQ.1.AND.IPARCD.EQ.1.AND.NU.NE.0) ISELK = 1
C
C         CASE 2: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH EVEN
          IF(IPARAB.EQ.0.AND.IPARCD.EQ.0) ISELK = 1
C
        ENDDO
        IF(ISELK.EQ.0) GOTO 6001
C
      ENDIF
C
C**********************************************************************C
C     LOOP OVER |MQN| PROJECTIONS C AND D (INDEX 7000)                 C
C**********************************************************************C
C
C     LOOP OVER |MQN(C)| VALUES
      DO 7000 MC=1,IABS(KQN(3))
        MQN(3) = 2*MC-1
C
C     LOOP OVER |MQN(D)| VALUES
      DO 7000 MD=1,IABS(KQN(4))
        MQN(4) = 2*MD-1
C
C     INDEX ASSIGNMENT
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
C     SKIP CONTRIBUTIONS THAT ARISE BY PERMUTATION OF INTEGRALS
      IF(G2INT.EQ.'COULM') THEN
        IF(IQ3.LT.IQ4) GOTO 7001
      ENDIF
C      IF(IQ3.LT.IQ4) GOTO 7001
C
C     EQ-COEFFICIENT PHASE FACTORS FOR PERMUTATION OF R-INTEGRALS
      IF(G2INT.EQ.'COULM') THEN
        PCD1 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
        PCD2 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
      ELSEIF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') THEN
        PCD1 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
        PCD2 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
      ENDIF
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON MQN                           C
C**********************************************************************C
C
C     SPIN PROJECTION CONSERVED ALONG Z-AXIS FOR LINEAR MOLECULES
      IF(ISYM.EQ.1.OR.ISYM.EQ.2) THEN
        ISELM = 0
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) ISELM = 1
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) ISELM = 1
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) ISELM = 1
        IF(ISELM.EQ.0) GOTO 7001
      ENDIF
C
C     ATOM-CENTRED SELECTION RULES (ONLY APPLIES IF RACAH1 SWITCHED OFF)
      IF(MCNT.EQ.1) THEN
        ISELM = 0
        DO ISGN1=1,2
          DO ISGN2=1,2
            DO ISGN3=1,2
              DO ISGN4=1,2
                MMJA = MQN(1)*((-1)**ISGN1)
                MMJB = MQN(2)*((-1)**ISGN2)
                MMJC = MQN(3)*((-1)**ISGN3)
                MMJD = MQN(4)*((-1)**ISGN4)
                IF(MMJA-MMJB.EQ.MMJD-MMJC) ISELM = 1
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IF(ISELM.EQ.0) GOTO 7001
      ENDIF
C
C**********************************************************************C
C     LOOP OVER COMPONENT OVERLAP LABELS C AND D (INDEX 8000)          C
C**********************************************************************C
C
C     COMPONENT LABEL FOR C AND D: TT = LL(1) or SS(4) <- COULOMB
C                                  TT = LS(2) or SL(3) <- BREIT
      DO 8000 IT2=ITSTRT,ITSTOP,ITSKIP
C
C       EQ-COEFFICIENT STARTING ADDRESSES FOR (CD) PAIR
        IF(IT2.EQ.1) THEN
          ICDLL = IAD0LL(ICNTC,ICNTD,KC,KD,MC,MD)
        ELSEIF(IT2.EQ.2) THEN
          ICDLS = IADILS(ICNTC,ICNTD,KC,KD,MC,MD)
        ELSEIF(IT2.EQ.3) THEN
          ICDSL = IADISL(ICNTC,ICNTD,KC,KD,MC,MD)
        ELSEIF(IT2.EQ.4) THEN
          ICDSS = IAD0SS(ICNTC,ICNTD,KC,KD,MC,MD)
        ENDIF
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITN(2) = IT2
C
C       CALCULATE STARTING ADDRESS
        IF(ITN(2).EQ.1) THEN
          NADDC = 0
          NADDD = 0
        ELSEIF(ITN(2).EQ.2) THEN
          NADDC = 0
          NADDD = NSKP
        ELSEIF(ITN(2).EQ.3) THEN
          NADDC = NSKP
          NADDD = 0
        ELSEIF(ITN(2).EQ.4) THEN
          NADDC = NSKP
          NADDD = NSKP
        ENDIF
C
C       FLAG READ-IN OF E0(CD) COEFFICIENTS FOR THIS COMPONENT LABEL
        IECD = 1
C
C**********************************************************************C
C     FOCK MATRIX STARTING ADDRESSES                                   C
C**********************************************************************C
C
C     FOCK ADDRESS FOR EACH BASIS FUNCTION (WITH SPIN PROJECTION)
      NC1 = LRGE(ICNTC,KC,2*MC-1) + NADDC
      NC2 = LRGE(ICNTC,KC,2*MC  ) + NADDC
C
      ND1 = LRGE(ICNTD,KD,2*MD-1) + NADDD
      ND2 = LRGE(ICNTD,KD,2*MD  ) + NADDD
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS IN BLOCKS A AND B (INDEX 9000)         C
C**********************************************************************C
C
      DO 9000 IBAS=1,NBAS(1)
      DO 9000 JBAS=1,NBAS(2)
C
C     LIST ADDRESS FOR THIS IBAS,JBAS COMBINATION
      MIJ = (IBAS-1)*NBAS(2)+JBAS
C
C     RESET SCREENING COUNTERS
      DO M=1,NBAS(3)*NBAS(4)
        IMAP(M) = M
        ISCR(M) = 1
      ENDDO
      MAXN = NBAS(3)*NBAS(4)
C
C     BATCH OF TWO-BODY INTEGRALS (IJ|KL) FOR FIXED (IJ)
      CALL CPU_TIME(T1)
      IF(G2INT.EQ.'COULM') THEN
        CALL ERI(RR,XYZ,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN)
      ELSEIF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') THEN
        CALL BII(RR,XYZ,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN,BGAUNT)
      ENDIF
      CALL CPU_TIME(T2)
      TERI = TERI+T2-T1
C
C     CLEAR ARRAY FOR FIRST CONTRACTION (DIRECT)
      DO MKB=1,NBAS(3)*NUMO
        DO IJKSPIN=1,8
          B1(MKB,IJKSPIN) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     CLEAR ARRAY FOR FIRST CONTRACTION (SWAP)
      DO MLB=1,NBAS(4)*NUMO
        DO IJLSPIN=1,8
          B2(MLB,IJLSPIN) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     FIRST CONTRACTION:                                               C
C     (IJ;T|KL;T') -> (IJ;T|KB;T')  AND  (IJ;T|LK;T') -> (IJ;T|LB;T')  C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     FIRST CONTRACTION (NORMAL): (IJ;T|KL;T') -> (IJ;T|KB;T')
C
C     LOOP OVER BASIS FUNCTIONS IN BLOCK C AND OCCUPIED STATES IOCCB
      DO KBAS=1,NBAS(3)
        DO IOCCB=1,NUMO
C
C         FOCK MATRIX ADDRESS FOR IOCCB
          IB = MINO-1+IOCCB+NSKP
C
C         LIST ADDRESS FOR THIS KBAS,IOCCB COMBINATION
          MKB = (KBAS-1)*NUMO + IOCCB
C
C         LOOP OVER BASIS FUNCTIONS IN BLOCK D AND CONTRACT OVER ERI
          DO LBAS=1,NBAS(4)
C
C           LIST ADDRESS FOR THIS KBAS,LBAS COMBINATION
            M = (KBAS-1)*NBAS(4) + LBAS
C
C           (--|-B) = (--|--) + (--|-+)
            B1(MKB,1) = B1(MKB,1) +      RR(M, 1)*COEF(ND1+LBAS,IB)
     &                            +      RR(M, 2)*COEF(ND2+LBAS,IB)
C           (+-|-B) = (+-|--) + (+-|-+)
            B1(MKB,2) = B1(MKB,2) +      RR(M, 9)*COEF(ND1+LBAS,IB)
     &                            +      RR(M,10)*COEF(ND2+LBAS,IB)
C           (-+|-B) = (-+|--) + (-+|-+)
            B1(MKB,3) = B1(MKB,3) +      RR(M, 5)*COEF(ND1+LBAS,IB)
     &                            +      RR(M, 6)*COEF(ND2+LBAS,IB)
C           (++|-B) = (++|--) + (++|-+)
            B1(MKB,4) = B1(MKB,4) +      RR(M,13)*COEF(ND1+LBAS,IB)
     &                            +      RR(M,14)*COEF(ND2+LBAS,IB)
C           (--|+B) = (--|+-) + (--|++)
            B1(MKB,5) = B1(MKB,5) +      RR(M, 3)*COEF(ND1+LBAS,IB)
     &                            +      RR(M, 4)*COEF(ND2+LBAS,IB)
C           (+-|+B) = (+-|+-) + (+-|++)
            B1(MKB,6) = B1(MKB,6) +      RR(M,11)*COEF(ND1+LBAS,IB)
     &                            +      RR(M,12)*COEF(ND2+LBAS,IB)
C           (-+|+B) = (-+|+-) + (-+|++)
            B1(MKB,7) = B1(MKB,7) +      RR(M, 7)*COEF(ND1+LBAS,IB)
     &                            +      RR(M, 8)*COEF(ND2+LBAS,IB)
C           (++|+B) = (++|+-) + (++|++)
            B1(MKB,8) = B1(MKB,8) +      RR(M,15)*COEF(ND1+LBAS,IB)
     &                            +      RR(M,16)*COEF(ND2+LBAS,IB)
C
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NBAS(4)*NUMO*NBAS(3)
C
C     FIRST CONTRACTION (SWAP): (IJ;T|LK;T') -> (IJ;T|LB;T')
      IF(G2INT.EQ.'COULM'.AND.IQ3.EQ.IQ4) GOTO 9001
      IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') GOTO 9001
C
C     LOOP OVER BASIS FUNCTIONS IN BLOCK D AND OCCUPIED STATES IOCCB
      DO LBAS=1,NBAS(4)
        DO IOCCB=1,NUMO
C
C         FOCK MATRIX ADDRESS FOR IOCCB
          IB = MINO-1+IOCCB+NSKP
C
C         LIST ADDRESS FOR THIS LBAS,IOCCB COMBINATION
          MLB = (LBAS-1)*NUMO + IOCCB
C
C         LOOP OVER BASIS FUNCTIONS IN BLOCK C AND CONTRACT OVER ERI
          DO KBAS=1,NBAS(3)
C
C           LIST ADDRESS FOR THIS KBAS,LBAS COMBINATION
            M = (KBAS-1)*NBAS(4) + LBAS
C
C           (--|+B) = PAB*{(--|B+)} = PAB*{(--|++) + ((--|-+))}
            B2(MLB,1) = B2(MLB,1) + PCD1*RR(M, 4)*COEF(NC1+KBAS,IB)
     &                            + PCD2*RR(M, 2)*COEF(NC2+KBAS,IB)
C
C           (+-|+B) = PAB*{(+-|B+)} = PAB*{(+-|++) + ((+-|-+))}
            B2(MLB,2) = B2(MLB,2) + PCD1*RR(M,12)*COEF(NC1+KBAS,IB)
     &                            + PCD2*RR(M,10)*COEF(NC2+KBAS,IB)
C
C           (-+|+B) = PAB*{(-+|B+)} = PAB*{(-+|++) + ((-+|-+))}
            B2(MLB,3) = B2(MLB,3) + PCD1*RR(M, 8)*COEF(NC1+KBAS,IB)
     &                            + PCD2*RR(M, 6)*COEF(NC2+KBAS,IB)
C
C           (++|+B) = PAB*{(++|B+)} = PAB*{(++|++) + ((++|-+))}
            B2(MLB,4) = B2(MLB,4) + PCD1*RR(M,16)*COEF(NC1+KBAS,IB)
     &                            + PCD2*RR(M,14)*COEF(NC2+KBAS,IB)
C
C           (--|-B) = PAB*{(--|B-)} = PAB*{(--|+-) + ((--|--))}
            B2(MLB,5) = B2(MLB,5) + PCD2*RR(M, 3)*COEF(NC1+KBAS,IB)
     &                            + PCD1*RR(M, 1)*COEF(NC2+KBAS,IB)
C
C           (+-|-B) = PAB*{(+-|B-)} = PAB*{(+-|+-) + ((+-|--))}
            B2(MLB,6) = B2(MLB,6) + PCD2*RR(M,11)*COEF(NC1+KBAS,IB)
     &                            + PCD1*RR(M, 9)*COEF(NC2+KBAS,IB)
C
C           (-+|-B) = PAB*{(-+|B-)} = PAB*{(-+|+-) + ((-+|--))}
            B2(MLB,7) = B2(MLB,7) + PCD2*RR(M, 7)*COEF(NC1+KBAS,IB)
     &                            + PCD1*RR(M, 5)*COEF(NC2+KBAS,IB)
C
C           (++|-B) = PAB*{(++|B-)} = PAB*{(++|+-) + ((++|--))}
            B2(MLB,8) = B2(MLB,8) + PCD2*RR(M,15)*COEF(NC1+KBAS,IB)
     &                            + PCD1*RR(M,13)*COEF(NC2+KBAS,IB)
C
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NBAS(4)*NUMO*NBAS(3)
C
C     SKIP POINT FOR IQ3=IQ4
9001  CONTINUE
C
      CALL CPU_TIME(T2)
      TCN1 = TCN1+T2-T1
C
C**********************************************************************C
C     SECOND CONTRACTION:                                   ~          C
C     (IJ;T|KB;T') -> (IJ;T|DB)  AND  (IJ;T|LB;T') -> (IJ;T|DB)        C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     SECOND CONTRACTION (NORMAL): (IJ;T|KB) -> (IJ;T|DB)
C
C     LOOP OVER OCCUPIED STATES IOCCB AND IOCCD
      DO IOCCB=1,NUMO
        DO IOCCD=1,NUMO
C
C         FOCK MATRIX ADDRESS FOR IOCCD
          ID = MINO-1+IOCCD+NSKP
C
C         LIST ADDRESS FOR THIS IOCCD,IOCCB COMBINATION IN B1
          MDB = (IOCCB-1)*NUMO + IOCCD
C
C         LOOP OVER BASIS FUNCTIONS IN BLOCK C AND CONTRACT OVER B1
          DO KBAS=1,NBAS(3)
C
C           LIST ADDRESS FOR THIS KBAS,IOCCB COMBINATION
            MKB = (KBAS-1)*NUMO + IOCCB
C
C           (--|DB) = (--|-B) + (--|+B)
            DB(MIJ,MDB,1) = DB(MIJ,MDB,1)
     &                    + B1(MKB,1)*DCONJG(COEF(NC1+KBAS,ID))
     &                    + B1(MKB,5)*DCONJG(COEF(NC2+KBAS,ID))
C
C           (-+|DB) = (-+|-B) + (-+|+B)
            DB(MIJ,MDB,2) = DB(MIJ,MDB,2)
     &                    + B1(MKB,3)*DCONJG(COEF(NC1+KBAS,ID))
     &                    + B1(MKB,7)*DCONJG(COEF(NC2+KBAS,ID))
C
C           (+-|DB) = (+-|-B) + (+-|+B)
            DB(MIJ,MDB,3) = DB(MIJ,MDB,3)
     &                    + B1(MKB,2)*DCONJG(COEF(NC1+KBAS,ID))
     &                    + B1(MKB,6)*DCONJG(COEF(NC2+KBAS,ID))
C
C           (++|DB) = (++|-B) + (++|+B)
            DB(MIJ,MDB,4) = DB(MIJ,MDB,4)
     &                    + B1(MKB,4)*DCONJG(COEF(NC1+KBAS,ID))
     &                    + B1(MKB,8)*DCONJG(COEF(NC2+KBAS,ID))
C
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NUMO*NUMO*NBAS(3)
C
C     SECOND CONTRACTION (SWAP): (IJ;T|LB) -> (IJ;T|DB)
      IF(G2INT.EQ.'COULM'.AND.IQ3.EQ.IQ4) GOTO 9002
      IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') GOTO 9002
C
C     LOOP OVER OCCUPIED STATES IOCCB AND IOCCD
      DO IOCCB=1,NUMO
        DO IOCCD=1,NUMO
C
C         FOCK MATRIX ADDRESS FOR IOCCD
          ID = MINO-1+IOCCD+NSKP
C
C         LIST ADDRESS FOR THIS IOCCD,IOCCB COMBINATION IN B1
          MDB = (IOCCB-1)*NUMO + IOCCD
C
C         LOOP OVER BASIS FUNCTIONS IN BLOCK D AND CONTRACT OVER B2
          DO LBAS=1,NBAS(4)
C
C           LIST ADDRESS FOR THIS LBAS,IOCCB COMBINATION
            MLB = (LBAS-1)*NUMO + IOCCB
C
C           (--|DB) = (--|+B) + (--|-B)
            DB(MIJ,MDB,1) = DB(MIJ,MDB,1)
     &                    + B2(MLB,1)*DCONJG(COEF(ND1+LBAS,ID))
     &                    + B2(MLB,5)*DCONJG(COEF(ND2+LBAS,ID))
C
C           (-+|DB) = (-+|+B) + (-+|-B)
            DB(MIJ,MDB,2) = DB(MIJ,MDB,2)
     &                    + B2(MLB,3)*DCONJG(COEF(ND1+LBAS,ID))
     &                    + B2(MLB,7)*DCONJG(COEF(ND2+LBAS,ID))
C
C           (+-|DB) = (+-|+B) + (+-|-B)
            DB(MIJ,MDB,3) = DB(MIJ,MDB,3)
     &                    + B2(MLB,2)*DCONJG(COEF(ND1+LBAS,ID))
     &                    + B2(MLB,6)*DCONJG(COEF(ND2+LBAS,ID))
C
C           (++|DB) = (++|+B) + (++|-B)
            DB(MIJ,MDB,4) = DB(MIJ,MDB,4)
     &                    + B2(MLB,4)*DCONJG(COEF(ND1+LBAS,ID))
     &                    + B2(MLB,8)*DCONJG(COEF(ND2+LBAS,ID))
C
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NUMO*NUMO*NBAS(4)
C
C     SKIP POINT FOR IQ3=IQ4
9002  CONTINUE
C
      CALL CPU_TIME(T2)
      TCN2 = TCN2+T2-T1
C
C     END LOOP OVER BASIS PAIR (IBAS,JBAS)
9000  CONTINUE
C
C     END LOOP OVER COMPONENT OVERLAP T'T'
8000  CONTINUE
C
C     SKIP POINT FOR IQ3.LT.IQ4
7001  CONTINUE
C     ALL CONTRIBUTIONS FROM BLOCK (C,D) NOW ACCOUNTED FOR
7000  CONTINUE
C     SKIP POINT FOR KQN SELECTION RULES
6001  CONTINUE
6000  CONTINUE
5000  CONTINUE
C
C**********************************************************************C
C     THIRD CONTRACTION:                                               C
C     (IJ;T|DB) -> (IA;T|DB)  AND  (JI;T|DB) -> (JA;T|DB)              C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     THIRD CONTRACTION (NORMAL): (IJ;T|DB) -> (IA;T|DB)
C
C     CLEAR ARRAY FOR THIRD CONTRACTION (NORMAL)
      DO MADB=1,NUMO*NUMO*NUMO
        DO IBAS=1,NBAS(1)
          DO ISPIN=1,2
            ADB1(IBAS,MADB,ISPIN) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C     LOOP OVER OCCUPIED STATES IOCCA
      DO IOCCA=1,NUMO
C
C       FOCK MATRIX ADDRESS FOR IOCCA
        IA = MINO-1+IOCCA+NSKP
C
C       LOOP OVER OCCUPIED STATES IOCCB AND IOCCD
        DO IOCCB=1,IOCCA
          DO IOCCD=1,NUMO
C
C           LIST ADDRESS FOR THIS IOCCD,IOCCB COMBINATION
            MDB = (IOCCB-1)*NUMO+IOCCD
C
C           LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MDB
            MADB = (IOCCA-1)*NUMO*NUMO + MDB
C
C           LOOP OVER BASIS FUNCTIONS IN A AND B, CONTRACT OVER DB
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
C
C               LIST ADDRESS FOR THIS IBAS,JBAS COMBINATION
                MIJ = (IBAS-1)*NBAS(2)+JBAS
C
C               (-A|DB) = (--|DB) + (-+|DB)
                ADB1(IBAS,MADB,1) = ADB1(IBAS,MADB,1)
     &                      +      DB(MIJ,MDB,1)*COEF(NB1+JBAS,IA)
     &                      +      DB(MIJ,MDB,2)*COEF(NB2+JBAS,IA)
C               (+A|DB) = (+-|DB) + (++|DB)
                ADB1(IBAS,MADB,2) = ADB1(IBAS,MADB,2)
     &                      +      DB(MIJ,MDB,3)*COEF(NB1+JBAS,IA)
     &                      +      DB(MIJ,MDB,4)*COEF(NB2+JBAS,IA)
C
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NUMO*(NUMO+1)*NUMO*NBAS(1)*NBAS(2)/2
C
C     THIRD CONTRACTION (SWAP): (JI;T|DB) -> (JA;T|DB)
      IF(G2INT.EQ.'COULM'.AND.IQ1.EQ.IQ2) GOTO 4001
      IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') GOTO 4001
C
C     CLEAR ARRAY FOR THIRD CONTRACTION (SWAP)
      DO MADB=1,NUMO*NUMO*NUMO
        DO JBAS=1,NBAS(2)
          DO JSPIN=1,2
            ADB2(JBAS,MADB,JSPIN) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C     LOOP OVER OCCUPIED STATES IOCCA
      DO IOCCA=1,NUMO
C
C       FOCK MATRIX ADDRESS FOR IOCCA
        IA = MINO-1+IOCCA+NSKP
C
C       LOOP OVER OCCUPIED STATES IOCCB AND IOCCD
        DO IOCCB=1,IOCCA
          DO IOCCD=1,NUMO
C
C           LIST ADDRESS FOR THIS IOCCD,IOCCB COMBINATION
            MDB = (IOCCB-1)*NUMO+IOCCD
C
C           LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MDB
            MADB = (IOCCA-1)*NUMO*NUMO + MDB
C
C           LOOP OVER BASIS FUNCTIONS IN A AND B, CONTRACT OVER DB
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
C
C               LIST ADDRESS FOR THIS IBAS,JBAS COMBINATION
                MIJ = (IBAS-1)*NBAS(2)+JBAS
C
C               (+A|DB) = PCD*{(+A|DB)} = PCD*{(++|DB) + (-+|DB)}
                ADB2(JBAS,MADB,1) = ADB2(JBAS,MADB,1)
     &                      + PAB1*DB(MIJ,MDB,4)*COEF(NA1+IBAS,IA)
     &                      + PAB2*DB(MIJ,MDB,2)*COEF(NA2+IBAS,IA)
C               (-A|DB) = PCD*{(-A|DB)} = PCD*{(+-|DB) + (--|DB)}
                ADB2(JBAS,MADB,2) = ADB2(JBAS,MADB,2)
     &                      + PAB2*DB(MIJ,MDB,3)*COEF(NA1+IBAS,IA)
     &                      + PAB1*DB(MIJ,MDB,1)*COEF(NA2+IBAS,IA)
C
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NUMO*(NUMO+1)*NUMO*NBAS(1)*NBAS(2)/2
C
C     SKIP POINT FOR IQ1=IQ2
4001  CONTINUE
C
      CALL CPU_TIME(T2)
      TCN3 = TCN3+T2-T1
C
C**********************************************************************C
C     FOURTH CONTRACTION:                      ~                       C
C     (IA;T|DB) -> (CA|DB)  AND  (JA;T|DB) -> (CA|DB)                  C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     FOURTH CONTRACTION (NORMAL): (IA;T|DB) -> (CA|DB)
C
C     LOOP OVER OCCUPIED STATES IOCCA AND IOCCC
      DO IOCCA=1,NUMO
        DO IOCCC=1,NUMO
C
C         FOCK MATRIX ADDRESS FOR IOCCC
          IC = MINO-1+IOCCC+NSKP
C
C         LOOP OVER OCCUPIED STATES IOCCB AND IOCCD
          DO IOCCB=1,IOCCA
            DO IOCCD=1,NUMO
C
C             LIST ADDRESS FOR THIS IOCCD,IOCCB COMBINATION
              MDB = (IOCCB-1)*NUMO + IOCCD
C
C             LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MDB
              MADB = (IOCCA-1)*NUMO*NUMO + MDB
C
C             LIST ADDRESS FOR THIS IOCCC,IOCCA AND THE ABOVE MDB
              MCADB = (IOCCA-1)*NUMO*IOCCA*NUMO/2
     &              + (IOCCC-1)*NUMO*IOCCA + MDB
C
C             LOOP OVER BASIS FUNCTIONS IN BLOCK A, CONTRACT OVER ADB1
              DO IBAS=1,NBAS(1)
C
C               (CA|DB) = (-A|DB) + (+A|DB)
                CADB(MCADB) = CADB(MCADB)
     &                     + ADB1(IBAS,MADB,1)*DCONJG(COEF(NA1+IBAS,IC))
     &                     + ADB1(IBAS,MADB,2)*DCONJG(COEF(NA2+IBAS,IC))
C
              ENDDO
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NUMO*NUMO*(NUMO+1)*NUMO*NBAS(1)/2 
C                                              ~
C     FOURTH CONTRACTION (SWAP): (JA;T|DB) -> (CA|DB)
      IF(G2INT.EQ.'COULM'.AND.IQ1.EQ.IQ2) GOTO 4002
      IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') GOTO 4002
C
C     LOOP OVER OCCUPIED STATES IOCCA AND IOCCC
      DO IOCCA=1,NUMO
        DO IOCCC=1,NUMO
C
C         FOCK MATRIX ADDRESS FOR IOCCC
          IC = MINO-1+IOCCC+NSKP
C
C         LOOP OVER OCCUPIED STATES IOCCB AND IOCCD
          DO IOCCB=1,IOCCA
            DO IOCCD=1,NUMO
C
C             LIST ADDRESS FOR THIS IOCCD,IOCCB COMBINATION
              MDB = (IOCCB-1)*NUMO + IOCCD
C
C             LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MDB
              MADB = (IOCCA-1)*NUMO*NUMO + MDB
C
C             LIST ADDRESS FOR THIS IOCCC,IOCCA AND THE ABOVE MDB
              MCADB = (IOCCA-1)*NUMO*IOCCA*NUMO/2
     &              + (IOCCC-1)*NUMO*IOCCA + MDB
C
C             LOOP OVER BASIS FUNCTIONS IN BLOCK B, CONTRACT OVER ACB2
              DO JBAS=1,NBAS(2)
C
C               (CA|DB) = (+A|DB) + (-A|DB)
                CADB(MCADB) = CADB(MCADB)
     &                     + ADB2(JBAS,MADB,1)*DCONJG(COEF(NB1+JBAS,IC))
     &                     + ADB2(JBAS,MADB,2)*DCONJG(COEF(NB2+JBAS,IC))
C
              ENDDO
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NUMO*NUMO*(NUMO+1)*NUMO*NBAS(2)/2 
C
C     SKIP POINT FOR IQ1=IQ2
4002  CONTINUE
C
      CALL CPU_TIME(T2)
      TCN4 = TCN4+T2-T1
C
C     END LOOP OVER COMPONENT OVERLAP TT
4000  CONTINUE
C
C     WRITE PERCENTAGE COMPLETION
      IF(IOSAV.NE.'SILENT') THEN
        CALL PROGRSS(INEST,NNEST,TBEG)
      ENDIF
C
C     SKIP POINT FOR IQ1.LT.IQ2
3001  CONTINUE
C
C     ALL CONTRIBUTIONS FROM BLOCK (A,B) NOW ACCOUNTED FOR
3000  CONTINUE
2000  CONTINUE
1000  CONTINUE
C
C**********************************************************************C
C     CALCULATE SECOND-ORDER PAIR CORRELATION ENERGY FROM (CA|DB)      C
C**********************************************************************C
C
C     WRITE CADB ARRAY TO EXTERNAL FILE IF PROMPTED
10    FORMAT(/1X,72('=')/24(' '),'Writing data file in MBPT1'
     &       /1X,72('-')/' Operator: ',A/1X,'File name: ',A
     &       /1X,'File size: ',F7.3,A,26X,'# words: ',I16/1X,72('='))
      IF(IOSAV.EQ.'WRTOUT'.OR.IOSAV.EQ.'SILENT') THEN
        CADBOUT = TRIM(OUTFL)//'_CADB_'//G2INT//'.dat'
        NCADB   = (NUMO+1)*NUMO*NUMO*NUMO/2
        IF(IOSAV.NE.'SILENT') THEN
          SCADB = DFLOAT(NCADB)*16.0D-6
          IF(SCADB.LT.1.0D0) THEN
            WRITE(6,10) G2INT,TRIM(CADBOUT),1.0D+3*SCADB,' KB',NCADB
            WRITE(7,10) G2INT,TRIM(CADBOUT),1.0D+3*SCADB,' KB',NCADB
          ELSEIF(SCADB.LT.1.0D3) THEN
            WRITE(6,10) G2INT,TRIM(CADBOUT),1.0D+0*SCADB,' MB',NCADB
            WRITE(7,10) G2INT,TRIM(CADBOUT),1.0D+0*SCADB,' MB',NCADB
          ELSE
            WRITE(6,10) G2INT,TRIM(CADBOUT),1.0D-3*SCADB,' GB',NCADB
            WRITE(7,10) G2INT,TRIM(CADBOUT),1.0D-3*SCADB,' GB',NCADB
          ENDIF
        ENDIF
        OPEN(UNIT=8,FILE=TRIM(CADBOUT),STATUS='UNKNOWN')
        REWIND(UNIT=8)
C       WRITE MINO AND NUMO
        WRITE(8, *) MINO,NUMO
C       LOOP OVER OCCUPIED STATES IOCCA AND VIRTUAL STATES IVRTR
        DO IOCCC=1,NUMO
          DO IOCCA=1,NUMO
C           LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
            DO IOCCD=1,NUMO
              DO IOCCB=1,IOCCA
C               MAIN LIST ADDRESS FOR THIS IOCCA,IVRTR,IOCCB,IVRTS
                MCADB = (IOCCA-1)*NUMO*IOCCA*NUMO/2
     &                + (IOCCC-1)*NUMO*IOCCA + (IOCCB-1)*NUMO + IOCCD
                WRITE(8, *) CADB(MCADB)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CLOSE(UNIT=8)
      ENDIF
C
C     SKIP POINT FOR READIN OPTION
80    CONTINUE
C
      CALL CPU_TIME(T1)
C
C     FOR EACH IOCCA,IOCCB PAIR, SUM OVER IVRTR AND IVRTS CONTRIBUTIONS
C
C     LOOP OVER OCCUPIED STATES IOCCA AND VIRTUAL STATES IVRTR
      DO IOCCC=1,NUMO
        DO IOCCA=1,NUMO
C
C         LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
          DO IOCCD=1,NUMO
            DO IOCCB=1,IOCCA
C
C             MAIN LIST ADDRESS FOR THIS IOCCA,IVRTR,IOCCB,IVRTS
              MCADB = (IOCCA-1)*NUMO*IOCCA*NUMO/2
     &              + (IOCCC-1)*NUMO*IOCCA + (IOCCB-1)*NUMO + IOCCD
C
              IF(IOCCA.EQ.IOCCC.AND.IOCCB.EQ.IOCCD) THEN
                EAB1(IOCCA,IOCCB,4) = DREAL(CADB(MCADB))
              ENDIF
C
              IF(IOCCB.EQ.IOCCC.AND.IOCCA.EQ.IOCCD) THEN
                EAB1(IOCCA,IOCCB,5) =-DREAL(CADB(MCADB))
              ENDIF
C
C             ADD TO DIRECT AND EXCHANGE BINS
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     FILL IN THE OTHER HALF OF THE ARRAY AND CALCULATE TOTALS
      E1D = 0.0D0
      E1X = 0.0D0
      DO IOCCA=1,NUMO
        DO IOCCB=1,IOCCA
C
C         INTERMEDIATE VALUES
          EAB1ONE = EAB1(IOCCA,IOCCB,3)
          EAB1DIR = EAB1(IOCCA,IOCCB,4)
          EAB1XCH = EAB1(IOCCA,IOCCB,5)
          EAB1TWO = EAB1DIR + EAB1XCH
C
C         PUT THESE INTO EAB1 AND ADD CONTRIBUTION TO E1
          EAB1(IOCCA,IOCCB,6) = EAB1TWO
          IF(IOCCA.NE.IOCCB) THEN
            EAB1(IOCCB,IOCCA,4) = EAB1DIR
            EAB1(IOCCB,IOCCA,5) = EAB1XCH
            EAB1(IOCCB,IOCCA,6) = EAB1TWO
            E1D = E1D +       EAB1DIR
            E1X = E1X +       EAB1XCH
          ELSE
            E1D = E1D + 0.5D0*EAB1DIR
            E1X = E1X + 0.5D0*EAB1XCH
          ENDIF
        ENDDO
      ENDDO
      E1S = E1H+E1D+E1X
C
C     WRITE RESULTS OF EAB ENERGIES TO AN EXTERNAL FILE
      EAB1OUT = TRIM(OUTFL)//'_EAB1_'//G2INT//'.dat'
      OPEN(UNIT=8,FILE=TRIM(EAB1OUT),STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO IOCCA=1,NUMO
        DO IOCCB=1,NUMO
          WRITE(8, *) (EAB1(IOCCA,IOCCB,N),N=1,6)
        ENDDO
      ENDDO
      CLOSE(UNIT=8)
C
C**********************************************************************C
C     CALCULATE SECOND-ORDER SINGLE ORBITAL ENERGY                     C
C**********************************************************************C
C
C     FOR EACH IOCCA, SUM OVER THE IOCCB CONTRIBUTIONS
      DO IOCCA=1,NUMO
        DO N=1,6
          EA1(IOCCA,N) = 0.0D0
          DO IOCCB=1,NUMO
            EA1(IOCCA,N) = EA1(IOCCA,N) + EAB1(IOCCA,IOCCB,N)
          ENDDO
        ENDDO
C       SIXTH COUNTER CAN STORE THE TOTAL SCF ENERGY
        EA1(IOCCA,6) = EA1(IOCCA,6) + EA1(IOCCA,3)
      ENDDO
C
      CALL CPU_TIME(T2)
      TSUM = TSUM+T2-T1
C
C**********************************************************************C
C     TERMINAL OUTPUT SUMMARY                                          C
C**********************************************************************C
C
      IF(IOSAV.EQ.'SILENT') RETURN
C
C     MBPT1 PAIRWISE SUMMARY
20    FORMAT(1X,A,11X,A,11X,A,11X,A,12X,A)
21    FORMAT(' (',I2,',',I2,')',F16.7,F16.7,F16.7,F17.7)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',25),'MBPT1 pairwise summary'
      WRITE(7, *) REPEAT(' ',25),'MBPT1 pairwise summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) '( a, b)','E1(H)','E1(J)','E1(K)','E1(G)'
      WRITE(7,20) '( a, b)','E1(H)','E1(J)','E1(K)','E1(G)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO IOCCA=1,NUMO
        IANUM = IOCCA+MINO-1
        DO IOCCB=1,IOCCA
          IBNUM = IOCCB+MINO-1
          WRITE(6,21) IANUM,IBNUM,(EAB1(IOCCA,IOCCB,N),N=3,6)
          WRITE(7,21) IANUM,IBNUM,(EAB1(IOCCA,IOCCB,N),N=3,6)
        ENDDO
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     MBPT1 SINGLE-PARTICLE SUMMARY
30    FORMAT(1X,A,11X,A,11X,A,11X,A,10X,A)
31    FORMAT('  ',I2,'    ',F16.7,F16.7,F16.7,F17.7)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',21),'MBPT1 single particle summary'
      WRITE(7, *) REPEAT(' ',21),'MBPT1 single particle summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,30) '  a    ','E1(H)','E1(J)','E1(K)','E1(tot)'
      WRITE(7,30) '  a    ','E1(H)','E1(J)','E1(K)','E1(tot)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO IOCCA=1,NUMO
        IANUM = IOCCA+MINO-1
        WRITE(6,31) IANUM,(EA1(IOCCA,N),N=3,6)
        WRITE(7,31) IANUM,(EA1(IOCCA,N),N=3,6)
      ENDDO
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     MBPT1 TOTAL FIRST-ORDER INTERACTION
32    FORMAT(' total  ',F16.7,F16.7,F16.7,F17.7)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',19),'MBPT1 first order molecular energy'
      WRITE(7, *) REPEAT(' ',19),'MBPT1 first order molecular energy'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,32) E1H,E1D,E1X,E1S
      WRITE(7,32) E1H,E1D,E1X,E1S
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6, *) E1D,E1X
      WRITE(7, *) E1D,E1X
C
C     MBPT1 LABOUR ANALYSIS
      CALL CPU_TIME(TFIN)
      TTOT = TFIN - TBEG
      TOTH = TTOT - (T1EL + TERI + TCN1 + TCN2 + TSUM)

40    FORMAT(1X,A,15X,A)
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',26),'MBPT1 labour analysis'
      WRITE(7, *) REPEAT(' ',26),'MBPT1 labour analysis'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
      WRITE(6,40) 'One-body terms - EH(A,B)                 ',HMS(T1EL)
      WRITE(7,40) 'One-body terms - EH(A,B)                 ',HMS(T1EL)
      WRITE(6,40) 'ERI construction - (IJ|KL)               ',HMS(TERI)
      WRITE(7,40) 'ERI construction - (IJ|KL)               ',HMS(TERI)
      WRITE(6,40) '1st contraction  - (IJ|KB)               ',HMS(TCN1)
      WRITE(7,40) '1st contraction  - (IJ|KB)               ',HMS(TCN1)
      WRITE(6,40) '2nd contraction  - (IJ|BB) and (IJ|AB)   ',HMS(TCN2)
      WRITE(7,40) '2nd contraction  - (IJ|BB) and (IJ|AB)   ',HMS(TCN2)
      WRITE(6,40) '3rd contraction  - (IA|BB) and (IA|AB)   ',HMS(TCN3)
      WRITE(7,40) '3rd contraction  - (IA|BB) and (IA|AB)   ',HMS(TCN3)
      WRITE(6,40) '4th contraction  - (AA|BB) and (BA|AB)   ',HMS(TCN4)
      WRITE(7,40) '4th contraction  - (AA|BB) and (BA|AB)   ',HMS(TCN4)
      WRITE(6,40) 'Other                                    ',HMS(TOTH)
      WRITE(7,40) 'Other                                    ',HMS(TOTH)
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,40) 'Total MBPT1 time                         ',HMS(TTOT)
      WRITE(7,40) 'Total MBPT1 time                         ',HMS(TTOT)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
      RETURN
      END
C
C
      SUBROUTINE MBPT2(MINO,NUMO,MINV,NUMV,G2INT,IOSAV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            MM       MM BBBBBBB  PPPPPPP TTTTTTTT 222222              C
C            MMM     MMM BB    BB PP    PP   TT   22    22             C
C            MMMM   MMMM BB    BB PP    PP   TT         22             C
C            MM MM MM MM BBBBBBB  PP    PP   TT       22               C
C            MM  MMM  MM BB    BB PPPPPPP    TT     22                 C
C            MM   M   MM BB    BB PP         TT   22                   C
C            MM       MM BBBBBBB  PP         TT   22222222             C
C                                                                      C
C -------------------------------------------------------------------- C
C  MBPT2 CALCULATES SECOND-ORDER PAIR CORRELATION CORRECTIONS OVER A   C
C  USER-SPECIFIED TWO-BODY OPERATOR. IT USES A RELATIVISTIC ADAPTION   C
C  OF THE DIRECT MP2 ALGORITHM OF HEAD-GORDON, POPLE AND FRISCH (1988).C
C  THIS IS NOT MOLLER-PLESSET PERTURBATION THEORY, SO CALL IT MBPT2.   C
C  ELECTRON REPULATION INTEGRALS ARE GENERATED ONCE, AT THE EXPENSE OF C
C  ADDITIONAL MEMORY STORAGE -- THEY ARE CONTRACTED IN FOUR STEPS.     C
C  THIS ALGORITHM EXPLOITS SYMMETRIES BETWEEN (A,B) AND (C,D) IN THE   C
C  INTEGRALS (AB|CD), BUT NOT THE SWAP (AB|CD)<->(CD|AB). STRUCTURE    C
C  MIMICS THAT OF 'COULOMB', BUT DOES NOT IMPLEMENT SELECTION RULES.   C
C  CONTRACTION ADDRESSES ARE CALCULATED EXPLICITLY RATHER THAN JUST    C
C  WITH UPDATING COUNTERS -- THIS IS TO GUIDE THE USER IN TRACKING.    C
C -------------------------------------------------------------------- C
C INPUT:                                                               C
C ▶ MINO  - LOWEST OCCUPIED STATE TO ACCOUNT FOR. (FULL: 1)            C
C ▶ NUMO  - NUMBER OF OCCUPIED STATES TO ACCOUNT FOR. (FULL: NOCC)     C
C ▶ MINV  - LOWEST VIRTUAL STATE TO ACCOUNT FOR. (FULL: NOCC+1)        C
C ▶ NUMV  - NUMBER OF VIRTUAL STATES TO ACCOUNT FOR. (FULL: NVRT)      C
C ▶ G2INT - NAME OF TWO-BODY OPERATOR ('COULM', 'GAUNT' OR 'BREIT').   C
C ▶ IOSAV - I/O OPTION (READ RESULTS, WRITE RESULTS, DO NOTHING).      C
C -------------------------------------------------------------------- C
C ▶ LOOP STRUCTURE NOT LIKE 'COULOMB' -- MUST CONTRACT C AND D ASAP.   C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL FILETHERE,BGAUNT
C
      CHARACTER*5  G2INT
      CHARACTER*6  IOSAV
      CHARACTER*16 HMS
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*60 RASBOUT,EAB1OUT,EAB2OUT
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4),LQN(4),JQN(4)
      DIMENSION INDEX(MCT,-(MEL+1):MEL,MKP),ITN(2)
      DIMENSION EAB2(NUMO,NUMO,6),EA2(NUMO,6)
C
      INTEGER*16 INEST,NNEST,NRASB
C
      COMPLEX*16 RR(MB2,16)
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 B(MBS*NUMO,8),SB(MB2,NUMO*NUMV,4)
      COMPLEX*16 ASB1(MBS,NUMO*NUMO*NUMV,2),ASB2(MBS,NUMO*NUMO*NUMV,2)
      COMPLEX*16 RASB((NUMO+1)*NUMO*NUMV*NUMV/2)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/E0LL/E0LLFL(MFL,4),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/E0SS/E0SSFL(MFL,4),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EILS/EILSFL(MFL,12),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/EISL/EISLFL(MFL,12),IADISL(MCT,MCT,MKP,MKP,MKP,MKP)
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/IQTT/IABLL,ICDLL,IABSS,ICDSS,IABLS,ICDLS,IABSL,ICDSL
      COMMON/IRCM/IEAB,IECD,IRIJ(MBS,MBS)
      COMMON/ISCR/IMTX(MB2,11),ISCR(MB2),IMAP(MB2),IBCH,ITOG,MAXN
C
C     TURN OFF RC(AB|CD) LOCAL FILE PROCESS
      RCFILE = .FALSE.
C
C     WARNINGS BASED ON INVALID HMLT VS. G2INT COMBINATIONS
1     FORMAT(1X,'In MBPT2: HMLT = ',A,' but G2INT = ',A,'.')
      IF(HMLT.EQ.'BARE') THEN
        WRITE(6,1) HMLT,G2INT
        WRITE(7,1) HMLT,G2INT
      ELSEIF(HMLT.EQ.'NORL'.OR.HMLT.EQ.'DHFR'.OR.HMLT.EQ.'DHFP') THEN
        IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') THEN
          WRITE(6,1) HMLT,G2INT
          WRITE(7,1) HMLT,G2INT
        ENDIF
      ENDIF
C
C     INTEGRAL SKIPPING ON MOLECULAR GROUP SYMMETRY CLASS BASIS
      IF(SHAPE.EQ.'ATOMIC') THEN
        ISYM = 2
      ELSEIF(SHAPE.EQ.'DIATOM'.OR.SHAPE.EQ.'LINEAR') THEN
        ISYM = 1
      ELSE
        ISYM = 0
      ENDIF
C
C     COMPONENT OVERLAP LABELS TO LOOP OVER
      IF(HMLT.EQ.'NORL') THEN
        IF(G2INT.EQ.'COULM') THEN
          ITSTRT = 1
          ITSTOP = 1
          ITSKIP = 1
        ENDIF
      ELSE
        IF(G2INT.EQ.'COULM') THEN
          ITSTRT = 4
          ITSTOP = 1
          ITSKIP =-3
        ELSEIF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') THEN
          ITSTRT = 2
          ITSTOP = 3
          ITSKIP = 1        
        ENDIF
      ENDIF
C
C     GAUNT REPLACEMENT TOGGLE
      IF(G2INT.EQ.'GAUNT') THEN
        BGAUNT = .TRUE.
      ELSE
        BGAUNT = .FALSE.
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
        RASBOUT = TRIM(OUTFL)//'_RASB_'//G2INT//'.dat'
        INQUIRE(FILE=TRIM(RASBOUT),EXIST=FILETHERE)
        IF(FILETHERE) THEN
          OPEN(UNIT=8,FILE=TRIM(RASBOUT),STATUS='OLD')
          REWIND(UNIT=8)
C         READ MINO AND NUMO
          READ(8, *) MINOTMP,NUMOTMP,MINVTMP,NUMVTMP
C         IF DIFFERENT TO INPUT, WRITE TO LOG AND UPDATE
          IF(MINOTMP.NE.MINO) THEN
            WRITE(6, *) 'In MBPT2: MINO clash.',MINO,MINOTMP,G2INT
            WRITE(7, *) 'In MBPT2: MINO clash.',MINO,MINOTMP,G2INT
          ENDIF
          IF(NUMOTMP.NE.NUMO) THEN
            WRITE(6, *) 'In MBPT2: NUMO clash.',NUMO,NUMOTMP,G2INT
            WRITE(7, *) 'In MBPT2: NUMO clash.',NUMO,NUMOTMP,G2INT
          ENDIF
          IF(MINVTMP.NE.MINV) THEN
            WRITE(6, *) 'In MBPT2: MINV clash.',MINV,MINVTMP,G2INT
            WRITE(7, *) 'In MBPT2: MINV clash.',MINV,MINVTMP,G2INT
          ENDIF
          IF(NUMVTMP.NE.NUMV) THEN
            WRITE(6, *) 'In MBPT2: NUMV clash.',NUMV,NUMVTMP,G2INT
            WRITE(7, *) 'In MBPT2: NUMV clash.',NUMV,NUMVTMP,G2INT
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
          GOTO 80
        ELSE
C       IF THE FILE DOESN'T YET EXIST, CALCULATE AND SAVE IT
          WRITE(6, *) 'In MBPT2: tried to read RASB but failed.'
          WRITE(7, *) 'In MBPT2: tried to read RASB but failed.'
          IOSAV = 'WRTOUT'
        ENDIF
      ENDIF
C
C**********************************************************************C
C     ORDERED INDEX OF (ICNT,KQN,MQN) COMBINATIONS                     C
C**********************************************************************C
C
      ICOUNT = 0
C
C     LOOP OVER NUCLEAR CENTRES
      DO ICT=1,NCNT
C
C       LOOP OVER KAPPA VALUES FOR THIS NUCLEAR CENTRE
        DO KN=1,NKAP(ICT)
C
C         IMPORT KAPPA, MAXIMUM MQN
          KAPPA = KAPA(KN,ICT)
          MJMAX = 2*IABS(KAPPA)-1
C
C         LOOP OVER MQN VALUES AND RECORD INDEX
          DO MJ=1,MJMAX,2
            ICOUNT              = ICOUNT+1
            INDEX(ICT,KAPPA,MJ) = ICOUNT
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     ORDERED COUNTER FOR NESTED LOOP PROGRESS                         C
C**********************************************************************C
C
      NNEST = 0
C
C     COMPONENT IT1, BLOCKS A AND B
      DO 100 IT1=ITSTRT,ITSTOP,ITSKIP
      DO 200 ICNTA=1,NCNT
      DO 200 KA=1,NKAP(ICNTA)
      KQN(1) = KAPA(KA,ICNTA)
      LQN(1) = LVAL(KQN(1))
      JQN(1) = 2*IABS(KQN(1))-1
      NBAS(1) = NFNC(LQN(1),ICNTA)
      DO 200 MA=1,IABS(KAPA(KA,ICNTA))
      MQN(1) = 2*MA-1
      IQ1 = INDEX(ICNTA,KQN(1),2*MA-1)
      DO 300 ICNTB=1,ICNTA
      DO 300 KB=1,NKAP(ICNTB)
      KQN(2) = KAPA(KB,ICNTB)
      LQN(2) = LVAL(KQN(2))
      JQN(2) = 2*IABS(KQN(2))-1
      NBAS(2) = NFNC(LQN(2),ICNTB)
      DO 300 MB=1,IABS(KAPA(KB,ICNTB))
      MQN(2) = 2*MB-1
      IQ2 = INDEX(ICNTB,KQN(2),2*MB-1)
C       REASONS TO SKIP
        IF(G2INT.EQ.'COULM') THEN
          IF(IQ1.LT.IQ2) GOTO 301
        ENDIF
C       COMPONENT IT2, BLOCK C
        DO 500 IT2=ITSTRT,ITSTOP,ITSKIP
        DO 600 ICNTC=1,NCNT
        DO 600 KC=1,NKAP(ICNTC)
        KQN(3) = KAPA(KC,ICNTC)
        LQN(3) = LVAL(KQN(3))
        JQN(3) = 2*IABS(KQN(3))-1
        NBAS(3) = NFNC(LQN(3),ICNTC)
        DO 600 MC=1,IABS(KAPA(KC,ICNTC))
        MQN(3) = 2*MC-1
        IQ3 = INDEX(ICNTC,KQN(3),2*MC-1)
C         BLOCK D
          DO 700 ICNTD=1,NCNT
          DO 700 KD=1,NKAP(ICNTD)
          KQN(4) = KAPA(KD,ICNTD)
          LQN(4) = LVAL(KQN(4))
          JQN(4) = 2*IABS(KQN(4))-1
          NBAS(4) = NFNC(LQN(4),ICNTD)
          DO 700 MD=1,IABS(KAPA(KD,ICNTD))
          MQN(4) = 2*MD-1
          IQ4 = INDEX(ICNTD,KQN(4),2*MD-1)
          MCNT = NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
          IF(MCNT.EQ.1) THEN
            IPARAB = MOD(LQN(1)+LQN(2),2)
            IPARCD = MOD(LQN(3)+LQN(4),2)
            IF(IPARAB.NE.IPARCD) THEN
              GOTO 701
            ENDIF
            NUI = MAX0(IABS(JQN(1)-JQN(2))/2,IABS(JQN(3)-JQN(4))/2)
            NUF = MIN0(    (JQN(1)+JQN(2))/2,    (JQN(3)+JQN(4))/2)
            IF(NUI.GT.NUF) THEN
              GOTO 701
            ENDIF
            ISELK = 0
            DO NU=NUI,NUF
              IPARAB = MOD(LQN(1)+LQN(2)+NU,2)
              IPARCD = MOD(LQN(3)+LQN(4)+NU,2)
              IF(IPARAB.EQ.1.AND.IPARCD.EQ.1.AND.NU.NE.0) ISELK = 1
              IF(IPARAB.EQ.0.AND.IPARCD.EQ.0) ISELK = 1
            ENDDO
            IF(ISELK.EQ.0) GOTO 701
            ISELM = 0
            DO ISGN1=1,2
              DO ISGN2=1,2
                DO ISGN3=1,2
                  DO ISGN4=1,2
                    MMJA = MQN(1)*((-1)**ISGN1)
                    MMJB = MQN(2)*((-1)**ISGN2)
                    MMJC = MQN(3)*((-1)**ISGN3)
                    MMJD = MQN(4)*((-1)**ISGN4)
                    IF(MMJA-MMJB.EQ.MMJD-MMJC) ISELM = 1
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            IF(ISELM.EQ.0) GOTO 701
          ENDIF
          IF(ISYM.EQ.1.OR.ISYM.EQ.2) THEN
            ISELM = 0
            IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) ISELM = 1
            IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) ISELM = 1
            IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) ISELM = 1
            IF(ISELM.EQ.0) GOTO 701
          ENDIF
C           FIRST CONTRACTION
            NNEST = NNEST + NBAS(1)*NBAS(2)*NBAS(3)*NBAS(4)*NUMO
701       CONTINUE
700       CONTINUE
C         SECOND CONTRACTION
          NNEST = NNEST + NBAS(1)*NBAS(2)*NBAS(3)*NUMO*NUMV
600     CONTINUE
500     CONTINUE
C       THIRD CONTRACTION
        NNEST = NNEST + NUMO*(NUMO+1)*NUMV*NBAS(1)*NBAS(2)/2
        IF(G2INT.EQ.'COULM'.AND.IQ1.EQ.IQ2) GOTO 302
        IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') GOTO 302
        NNEST = NNEST + NUMO*(NUMO+1)*NUMV*NBAS(1)*NBAS(2)/2
302     CONTINUE
C       FOURTH CONTRACTION
        NNEST = NNEST + NUMO*NUMV*(NUMO+1)*NUMV*NBAS(1)/2 
        IF(G2INT.EQ.'COULM'.AND.IQ1.EQ.IQ2) GOTO 303
        IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') GOTO 303
        NNEST = NNEST + NUMO*NUMV*(NUMO+1)*NUMV*NBAS(2)/2 
303   CONTINUE
301   CONTINUE
300   CONTINUE
200   CONTINUE
100   CONTINUE
C
C     START CONTRACTION COUNTER
      INEST = 0
C
C**********************************************************************C
C     LOOP OVER COMPONENT OVERLAP LABELS TT(AB) (USE INDEX 1000)       C
C**********************************************************************C
C
C     COMPONENT LABEL FOR A AND B: TT = LL(1) or SS(4) <- COULOMB
C                                  TT = LS(2) or SL(3) <- BREIT
      DO 1000 IT1=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITN(1) = IT1
C
C       CALCULATE STARTING ADDRESS
        IF(ITN(1).EQ.1) THEN
          NADDA = 0
          NADDB = 0
        ELSEIF(ITN(1).EQ.2) THEN
          NADDA = 0
          NADDB = NSKP
        ELSEIF(ITN(1).EQ.3) THEN
          NADDA = NSKP
          NADDB = 0
        ELSEIF(ITN(1).EQ.4) THEN
          NADDA = NSKP
          NADDB = NSKP
        ENDIF
C
C**********************************************************************C
C     LOOP OVER BLOCKS OF A (USE INDEX 2000)                           C
C**********************************************************************C
C
C     LOOP OVER CENTRE A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,ICNTA)
        LQN(1) = LVAL(KQN(1))
        JQN(1) = 2*IABS(KQN(1))-1
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),ICNTA)
        DO IBAS=1,NBAS(1)
          EXL(IBAS,1) = BEXL(IBAS,LQN(1),ICNTA)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 2000 MA=1,IABS(KQN(1))
        MQN(1) = 2*MA-1
C
C     INDEX ASSIGNMENT
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
C
C     STARTING FOCK ADDRESS FOR BASIS FUNCTIONS ON A
      NA1 = LRGE(ICNTA,KA,2*MA-1) + NADDA
      NA2 = LRGE(ICNTA,KA,2*MA  ) + NADDA
C
C**********************************************************************C
C     LOOP OVER BLOCKS OF B (USE INDEX 3000)                           C
C**********************************************************************C
C
C     LOOP OVER CENTRE B
      DO 3000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C     LOOP OVER KQN(B) VALUES
      DO 3000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KAPA(KB,ICNTB)
        LQN(2) = LVAL(KQN(2))
        JQN(2) = 2*IABS(KQN(2))-1
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
        ENDDO
C
C     LOOP OVER |MQN(B)| VALUES
      DO 3000 MB=1,IABS(KQN(2))
        MQN(2) = 2*MB-1
C
C     INDEX ASSIGNMENT
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
C
C     STARTING FOCK ADDRESS FOR BASIS FUNCTIONS ON B
      NB1 = LRGE(ICNTB,KB,2*MB-1) + NADDB
      NB2 = LRGE(ICNTB,KB,2*MB  ) + NADDB
C
C     EQ-COEFFICIENT STARTING ADDRESSES FOR (AB) PAIR
      IF(IT1.EQ.1) THEN
        IABLL = IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB)
      ELSEIF(IT1.EQ.2) THEN
        IABLS = IADILS(ICNTA,ICNTB,KA,KB,MA,MB)
      ELSEIF(IT1.EQ.3) THEN
        IABSL = IADISL(ICNTA,ICNTB,KA,KB,MA,MB)
      ELSEIF(IT1.EQ.4) THEN
        IABSS = IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB)
      ENDIF
C
C     FLAG READ-IN OF E0(AB) COEFFICIENTS FOR THIS COMPONENT LABEL
      IEAB = 1
C
C     EQ-COEFFICIENT PHASE FACTORS FOR PERMUTATION OF R-INTEGRALS
      IF(G2INT.EQ.'COULM') THEN
        PAB1 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
        PAB2 = ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      ELSEIF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') THEN
        PAB1 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**((-MQN(1)+MQN(2))/2))
        PAB2 =-ISIGN(1,KQN(1)*KQN(2))*DFLOAT((-1)**(( MQN(1)+MQN(2))/2))
      ENDIF
C
C     SKIP CONTRIBUTIONS THAT ARISE BY PERMUTATION OF INTEGRALS
      IF(G2INT.EQ.'COULM') THEN
        IF(IQ1.LT.IQ2) GOTO 3001
      ENDIF
C
C     CLEAR ARRAY FOR SECOND CONTRACTION (KL->SB)
      DO MIJ=1,NBAS(1)*NBAS(2)
        DO MSB=1,NUMO*NUMV
          DO IJSPIN=1,4
            SB(MIJ,MSB,IJSPIN) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS (IBAS,JBAS) (INDEX 4000)               C
C**********************************************************************C
C
      DO 4000 IBAS=1,NBAS(1)
      DO 4000 JBAS=1,NBAS(2)
C
C     LIST ADDRESS FOR THIS IBAS,JBAS COMBINATION
      MIJ = (IBAS-1)*NBAS(2)+JBAS
C
C**********************************************************************C
C     LOOP OVER COMPONENT OVERLAP LABELS TT(CD) (USE INDEX 5000)       C
C**********************************************************************C
C
C     COMPONENT LABEL FOR C AND D: TT = LL(1) or SS(4) <- COULOMB
C                                  TT = LS(2) or SL(3) <- BREIT
      DO 5000 IT2=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITN(2) = IT2
C
C       CALCULATE STARTING ADDRESS
        IF(ITN(2).EQ.1) THEN
          NADDC = 0
          NADDD = 0
        ELSEIF(ITN(2).EQ.2) THEN
          NADDC = 0
          NADDD = NSKP
        ELSEIF(ITN(2).EQ.3) THEN
          NADDC = NSKP
          NADDD = 0
        ELSEIF(ITN(2).EQ.4) THEN
          NADDC = NSKP
          NADDD = NSKP
        ENDIF
C
C**********************************************************************C
C     LOOP OVER BLOCKS OF C (USE INDEX 6000)                           C
C**********************************************************************C
C
C     LOOP OVER CENTRE C
      DO 6000 ICNTC=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE C
        XYZ(1,3) = BXYZ(1,ICNTC)
        XYZ(2,3) = BXYZ(2,ICNTC)
        XYZ(3,3) = BXYZ(3,ICNTC)
C
C     LOOP OVER KQN(C) VALUES
      DO 6000 KC=1,NKAP(ICNTC)
C
C       QUANTUM NUMBERS FOR BLOCK C
        KQN(3) = KAPA(KC,ICNTC)
        LQN(3) = LVAL(KQN(3))
        JQN(3) = 2*IABS(KQN(3))-1
C
C       BASIS EXPONENTS FOR BLOCK C
        NBAS(3) = NFNC(LQN(3),ICNTC)
        DO KBAS=1,NBAS(3)
          EXL(KBAS,3) = BEXL(KBAS,LQN(3),ICNTC)
        ENDDO
C
C     LOOP OVER |MQN(C)| VALUES
      DO 6000 MC=1,IABS(KQN(3))
        MQN(3) = 2*MC-1
C
C     INDEX ASSIGNMENT
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
C
C     STARTING FOCK ADDRESS FOR BASIS FUNCTIONS ON C
      NC1 = LRGE(ICNTC,KC,2*MC-1) + NADDC
      NC2 = LRGE(ICNTC,KC,2*MC  ) + NADDC
C
C     CLEAR ARRAY FOR FIRST CONTRACTION (L->B)
      DO MKB=1,NBAS(3)*NUMO
        DO IJKSPIN=1,8
          B(MKB,IJKSPIN) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     LOOP OVER BLOCKS OF D (USE INDEX 7000)                           C
C**********************************************************************C
C
C     LOOP OVER CENTRE D
      DO 7000 ICNTD=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE D
        XYZ(1,4) = BXYZ(1,ICNTD)
        XYZ(2,4) = BXYZ(2,ICNTD)
        XYZ(3,4) = BXYZ(3,ICNTD)
C
C     NUMBER OF NUCLEAR CENTRES INVOLVED IN THIS OVERLAP
      MCNT = NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
C
C     LOOP OVER KQN(D) VALUES
      DO 7000 KD=1,NKAP(ICNTD)
C
C       QUANTUM NUMBERS FOR BLOCK D
        KQN(4) = KAPA(KD,ICNTD)
        LQN(4) = LVAL(KQN(4))
        JQN(4) = 2*IABS(KQN(4))-1
C
C       BASIS EXPONENTS FOR BLOCK D
        NBAS(4) = NFNC(LQN(4),ICNTD)
        DO LBAS=1,NBAS(4)
          EXL(LBAS,4) = BEXL(LBAS,LQN(4),ICNTD)
        ENDDO
C
C     LOOP OVER |MQN(D)| VALUES
      DO 7000 MD=1,IABS(KQN(4))
        MQN(4) = 2*MD-1
C
C     INDEX ASSIGNMENT
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
C     STARTING FOCK ADDRESS FOR BASIS FUNCTIONS ON D
      ND1 = LRGE(ICNTD,KD,2*MD-1) + NADDD
      ND2 = LRGE(ICNTD,KD,2*MD  ) + NADDD
C
C     EQ-COEFFICIENT STARTING ADDRESSES FOR (CD) PAIR
      IF(IT2.EQ.1) THEN
        ICDLL = IAD0LL(ICNTC,ICNTD,KC,KD,MC,MD)
      ELSEIF(IT2.EQ.2) THEN
        ICDLS = IADILS(ICNTC,ICNTD,KC,KD,MC,MD)
      ELSEIF(IT2.EQ.3) THEN
        ICDSL = IADISL(ICNTC,ICNTD,KC,KD,MC,MD)
      ELSEIF(IT2.EQ.4) THEN
        ICDSS = IAD0SS(ICNTC,ICNTD,KC,KD,MC,MD)
      ENDIF
C
C     FLAG READ-IN OF E0(CD) COEFFICIENTS FOR THIS COMPONENT LABEL
      IECD = 1
C
C     EQ-COEFFICIENT PHASE FACTORS FOR PERMUTATION OF R-INTEGRALS
      IF(G2INT.EQ.'COULM') THEN
        PCD1 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
        PCD2 = ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
      ELSEIF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') THEN
        PCD1 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**((-MQN(3)+MQN(4))/2))
        PCD2 =-ISIGN(1,KQN(3)*KQN(4))*DFLOAT((-1)**(( MQN(3)+MQN(4))/2))
      ENDIF
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON KQN                           C
C**********************************************************************C
C
C     ATOM-CENTRED SELECTION RULES
      IF(MCNT.EQ.1) THEN
C
C       LQN PAIR PARITY (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQN(1)+LQN(2),2)
        IPARCD = MOD(LQN(3)+LQN(4),2)
C
C       LQN PAIR PARITY SELECTION RULE
        IF(IPARAB.NE.IPARCD) THEN
          GOTO 7001
        ENDIF
C
C       JQN TRIANGLE RULE CHECK FOR MULTIPOLE EXPANSION (ATOM-CENTRED)
        NUI = MAX0(IABS(JQN(1)-JQN(2))/2,IABS(JQN(3)-JQN(4))/2)
        NUF = MIN0(    (JQN(1)+JQN(2))/2,    (JQN(3)+JQN(4))/2)
        IF(NUI.GT.NUF) THEN
          GOTO 7001
        ENDIF
C
C       ADDITIONAL LQN SELECTION RULE PARITY ANALYSIS
        ISELK = 0
        DO NU=NUI,NUF
C
C         A AND B: LQN(1)+LQN(2)+NU EVEN OR ODD (0 IF EVEN, 1 IF ODD)
          IPARAB = MOD(LQN(1)+LQN(2)+NU,2)
C
C         C AND D: LQN(3)+LQN(4)+NU EVEN OR ODD (0 IF EVEN, 1 IF ODD)
          IPARCD = MOD(LQN(3)+LQN(4)+NU,2)
C
C         CASE 1: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH ODD (UNLESS NU=0)
          IF(IPARAB.EQ.1.AND.IPARCD.EQ.1.AND.NU.NE.0) ISELK = 1
C
C         CASE 2: LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH EVEN
          IF(IPARAB.EQ.0.AND.IPARCD.EQ.0) ISELK = 1
C
        ENDDO
        IF(ISELK.EQ.0) GOTO 7001
C
      ENDIF
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON MQN                           C
C**********************************************************************C
C
C     SPIN PROJECTION CONSERVED ALONG Z-AXIS FOR LINEAR MOLECULES
      IF(ISYM.EQ.1.OR.ISYM.EQ.2) THEN
        ISELM = 0
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) ISELM = 1
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) ISELM = 1
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) ISELM = 1
        IF(ISELM.EQ.0) GOTO 7001
      ENDIF
C
C     ATOM-CENTRED SELECTION RULES (ONLY APPLIES IF RACAH1 SWITCHED OFF)
      IF(MCNT.EQ.1) THEN
        ISELM = 0
        DO ISGN1=1,2
          DO ISGN2=1,2
            DO ISGN3=1,2
              DO ISGN4=1,2
                MMJA = MQN(1)*((-1)**ISGN1)
                MMJB = MQN(2)*((-1)**ISGN2)
                MMJC = MQN(3)*((-1)**ISGN3)
                MMJD = MQN(4)*((-1)**ISGN4)
                IF(MMJA-MMJB.EQ.MMJD-MMJC) ISELM = 1
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IF(ISELM.EQ.0) GOTO 7001
      ENDIF
C
C**********************************************************************C
C     ALL LOOPS NOW COMPLETE -- GENERATE BATCH OF ERIs AND CONTRACT    C
C**********************************************************************C
C
C     RESET SCREENING COUNTERS
      DO M=1,NBAS(3)*NBAS(4)
        IMAP(M) = M
        ISCR(M) = 1
      ENDDO
      MAXN = NBAS(3)*NBAS(4)
C
C     BATCH OF ELECTRON INTERACTION INTEGRALS (IJ|KL) FOR FIXED (IJ)
      CALL CPU_TIME(T1)
      IF(G2INT.EQ.'COULM') THEN
        CALL ERI(RR,XYZ,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN)
      ELSEIF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') THEN
        CALL BII(RR,XYZ,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN,BGAUNT)
      ENDIF
      CALL CPU_TIME(T2)
      TERI = TERI+T2-T1
C
C**********************************************************************C
C     FIRST CONTRACTION:                                               C
C     (IJ;T|KL;T') -> (IJ;T|KB;T')                                     C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     FIRST CONTRACTION (NORMAL): (IJ;T|KL;T') -> (IJ;T|KB;T')
C
C     LOOP OVER BASIS FUNCTIONS IN BLOCK C AND OCCUPIED STATES IOCCB
      DO KBAS=1,NBAS(3)
        DO IOCCB=1,NUMO
C
C         FOCK MATRIX ADDRESS FOR IOCCB
          IB = MINO-1+IOCCB+NSKP
C
C         LIST ADDRESS FOR THIS KBAS,IOCCB COMBINATION
          MKB = (KBAS-1)*NUMO + IOCCB
C
C         CONTRACT OVER ALL LBAS IN BLOCK D
          DO LBAS=1,NBAS(4)
C
C           LIST ADDRESS FOR THIS KBAS,LBAS COMBINATION
            M = (KBAS-1)*NBAS(4) + LBAS
C
C           (--|-B) = (--|--) + (--|-+)
            B(MKB,1) = B(MKB,1) + RR(M, 1)*COEF(ND1+LBAS,IB)
     &                          + RR(M, 2)*COEF(ND2+LBAS,IB)
C
C           (+-|-B) = (+-|--) + (+-|-+)
            B(MKB,2) = B(MKB,2) + RR(M, 9)*COEF(ND1+LBAS,IB)
     &                          + RR(M,10)*COEF(ND2+LBAS,IB)
C
C           (-+|-B) = (-+|--) + (-+|-+)
            B(MKB,3) = B(MKB,3) + RR(M, 5)*COEF(ND1+LBAS,IB)
     &                          + RR(M, 6)*COEF(ND2+LBAS,IB)
C
C           (++|-B) = (++|--) + (++|-+)
            B(MKB,4) = B(MKB,4) + RR(M,13)*COEF(ND1+LBAS,IB)
     &                          + RR(M,14)*COEF(ND2+LBAS,IB)
C
C           (--|+B) = (--|+-) + (--|++)
            B(MKB,5) = B(MKB,5) + RR(M, 3)*COEF(ND1+LBAS,IB)
     &                          + RR(M, 4)*COEF(ND2+LBAS,IB)
C
C           (+-|+B) = (+-|+-) + (+-|++)
            B(MKB,6) = B(MKB,6) + RR(M,11)*COEF(ND1+LBAS,IB)
     &                          + RR(M,12)*COEF(ND2+LBAS,IB)
C
C           (-+|+B) = (-+|+-) + (-+|++)
            B(MKB,7) = B(MKB,7) + RR(M, 7)*COEF(ND1+LBAS,IB)
     &                          + RR(M, 8)*COEF(ND2+LBAS,IB)
C
C           (++|+B) = (++|+-) + (++|++)
            B(MKB,8) = B(MKB,8) + RR(M,15)*COEF(ND1+LBAS,IB)
     &                          + RR(M,16)*COEF(ND2+LBAS,IB)
C
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NBAS(4)*NUMO*NBAS(3)
C
      CALL CPU_TIME(T2)
      TCN1 = TCN1+T2-T1
C
C     SKIP POINT FOR KQN AND MQN SELECTION RULES
7001  CONTINUE
C
C     FIRST CONTRACTION COMPLETE FOR THIS D BLOCK
7000  CONTINUE
C
C**********************************************************************C
C     SECOND CONTRACTION:                                              C
C     (IJ;T|KB;T') -> (IJ;T|SB)                                        C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     SECOND CONTRACTION (NORMAL): (IJ;T|KB;T') -> (IJ;T|SB;T')
C
C     LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
      DO IOCCB=1,NUMO
        DO IVRTS=1,NUMV
C
C         FOCK MATRIX ADDRESS FOR IVRTS
          IS = MINV-1+IVRTS+NSKP
C
C         LIST ADDRESS FOR THIS IVRTS,IOCCB COMBINATION IN B1
          MSB = (IOCCB-1)*NUMV + IVRTS
C
C         CONTRACT OVER ALL KBAS IN BLOCK C
          DO KBAS=1,NBAS(3)
C
C           LIST ADDRESS FOR THIS KBAS,IOCCB COMBINATION
            MKB = (KBAS-1)*NUMO + IOCCB
C
C           (--|SB) = (--|-B) + (--|+B)
            SB(MIJ,MSB,1) = SB(MIJ,MSB,1)
     &                             + B(MKB,1)*DCONJG(COEF(NC1+KBAS,IS))
     &                             + B(MKB,5)*DCONJG(COEF(NC2+KBAS,IS))
C
C           (-+|SB) = (-+|-B) + (-+|+B)
            SB(MIJ,MSB,2) = SB(MIJ,MSB,2)
     &                             + B(MKB,3)*DCONJG(COEF(NC1+KBAS,IS))
     &                             + B(MKB,7)*DCONJG(COEF(NC2+KBAS,IS))
C
C           (+-|SB) = (+-|-B) + (+-|+B)
            SB(MIJ,MSB,3) = SB(MIJ,MSB,3)
     &                             + B(MKB,2)*DCONJG(COEF(NC1+KBAS,IS))
     &                             + B(MKB,6)*DCONJG(COEF(NC2+KBAS,IS))
C
C           (++|SB) = (++|-B) + (++|+B)
            SB(MIJ,MSB,4) = SB(MIJ,MSB,4)
     &                             + B(MKB,4)*DCONJG(COEF(NC1+KBAS,IS))
     &                             + B(MKB,8)*DCONJG(COEF(NC2+KBAS,IS))
C
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NUMO*NUMV*NBAS(3)
C
      CALL CPU_TIME(T2)
      TCN2 = TCN2+T2-T1
C
C     SECOND CONTRACTION COMPLETE FOR THIS C BLOCK AND T'T' OVERLAP
6000  CONTINUE
5000  CONTINUE
C
C     SECOND CONTRACTION COMPLETE FOR ALL (IBAS,JBAS)
4000  CONTINUE
C
C     WRITE PERCENTAGE COMPLETION
      IF(IOSAV.NE.'SILENT') THEN
        CALL PROGRSS(INEST,NNEST,TBEG)
      ENDIF
C
C**********************************************************************C
C     THIRD CONTRACTION:                                               C
C     (IJ;T|SB) -> (IA;T|SB)  AND  (JI;T|SB) -> (JA;T|SB)              C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     THIRD CONTRACTION (DIRECT): (IJ;T|SB) -> (IA;T|SB)
C
C     CLEAR ARRAY FOR THIRD CONTRACTION (DIRECT)
      DO MASB=1,NUMV*NUMO*NUMO
        DO IBAS=1,NBAS(1)
          DO ISPIN=1,2
            ASB1(IBAS,MASB,ISPIN) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C     LOOP OVER OCCUPIED STATES IOCCA
      DO IOCCA=1,NUMO
C
C       FOCK MATRIX ADDRESS FOR IOCCA
        IA = MINO-1+IOCCA+NSKP
C
C       LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
        DO IOCCB=1,IOCCA
          DO IVRTS=1,NUMV
C
C           LIST ADDRESS FOR THIS IVRTS,IOCCB COMBINATION
            MSB = (IOCCB-1)*NUMV+IVRTS
C
C           LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MSB
            MASB = (IOCCA-1)*NUMV*NUMO + MSB
C
C           CONTRACT OVER ALL (IBAS,JBAS) IN BLOCKS A AND B
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
C
C               LIST ADDRESS FOR THIS IBAS,JBAS COMBINATION
                MIJ = (IBAS-1)*NBAS(2)+JBAS
C
C               (-A|SB) = (--|SB) + (-+|SB)
                ASB1(IBAS,MASB,1) = ASB1(IBAS,MASB,1)
     &                            +     SB(MIJ,MSB,1)*COEF(NB1+JBAS,IA)
     &                            +     SB(MIJ,MSB,2)*COEF(NB2+JBAS,IA)
C
C               (+A|SB) = (+-|SB) + (++|SB)
                ASB1(IBAS,MASB,2) = ASB1(IBAS,MASB,2)
     &                            +     SB(MIJ,MSB,3)*COEF(NB1+JBAS,IA)
     &                            +     SB(MIJ,MSB,4)*COEF(NB2+JBAS,IA)
C
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NUMO*(NUMO+1)*NUMV*NBAS(1)*NBAS(2)/2 
C
C     THIRD CONTRACTION (SWAP): (JI;T|SB) -> (JA;T|SB)
      IF(G2INT.EQ.'COULM'.AND.IQ1.EQ.IQ2) GOTO 3002
      IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') GOTO 3002
C
C     CLEAR ARRAY FOR THIRD CONTRACTION (SWAP)
      DO MASB=1,NUMV*NUMO*NUMO
        DO JBAS=1,NBAS(2)
          DO JSPIN=1,2
            ASB2(JBAS,MASB,JSPIN) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C     LOOP OVER OCCUPIED STATES IOCCA
      DO IOCCA=1,NUMO
C
C       FOCK MATRIX ADDRESS FOR IOCCA
        IA = MINO-1+IOCCA+NSKP
C
C       LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
        DO IOCCB=1,IOCCA
          DO IVRTS=1,NUMV
C
C           LIST ADDRESS FOR THIS IVRTS,IOCCB COMBINATION
            MSB = (IOCCB-1)*NUMV+IVRTS
C
C           LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MSB
            MASB = (IOCCA-1)*NUMV*NUMO + MSB
C
C           CONTRACT OVER ALL (IBAS,JBAS) IN BLOCKS A AND B
            DO IBAS=1,NBAS(1)
              DO JBAS=1,NBAS(2)
C
C               LIST ADDRESS FOR THIS IBAS,JBAS COMBINATION
                MIJ = (IBAS-1)*NBAS(2)+JBAS
C
C               (+A|SB) = PCD*{(+A|SB)} = PCD*{(++|SB) + (-+|SB)}
                ASB2(JBAS,MASB,1) = ASB2(JBAS,MASB,1)
     &                           + PAB1*SB(MIJ,MSB,4)*COEF(NA1+IBAS,IA)
     &                           + PAB2*SB(MIJ,MSB,2)*COEF(NA2+IBAS,IA)
C
C               (-A|SB) = PCD*{(-A|SB)} = PCD*{(+-|SB) + (--|SB)}
                ASB2(JBAS,MASB,2) = ASB2(JBAS,MASB,2)
     &                           + PAB2*SB(MIJ,MSB,3)*COEF(NA1+IBAS,IA)
     &                           + PAB1*SB(MIJ,MSB,1)*COEF(NA2+IBAS,IA)
C
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NUMO*(NUMO+1)*NUMV*NBAS(1)*NBAS(2)/2 
C
C     SKIP POINT FOR IQ1=IQ2
3002  CONTINUE
C
      CALL CPU_TIME(T2)
      TCN3 = TCN3+T2-T1
C
C**********************************************************************C
C     FOURTH CONTRACTION:                      ~                       C
C     (IA;T|SB) -> (RA|SB)  AND  (JA;T|SB) -> (RA|SB)                  C
C**********************************************************************C
C
      CALL CPU_TIME(T1)
C
C     FOURTH CONTRACTION (DIRECT): (IA;T|SB) -> (RA|SB)
C
C     LOOP OVER OCCUPIED STATES IOCCA AND VIRTUAL STATES IVRTR
      DO IOCCA=1,NUMO
        DO IVRTR=1,NUMV
C
C         FOCK MATRIX ADDRESS FOR IVRTR
          IR = MINV-1+IVRTR+NSKP
C
C         LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
          DO IOCCB=1,IOCCA
            DO IVRTS=1,NUMV
C
C             LIST ADDRESS FOR THIS IVRTS,IOCCB COMBINATION
              MSB = (IOCCB-1)*NUMV + IVRTS
C
C             LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MSB
              MASB = (IOCCA-1)*NUMV*NUMO + MSB
C
C             LIST ADDRESS FOR THIS IVRTR,IOCCA AND THE ABOVE MSB
              MRASB = (IOCCA-1)*NUMV*IOCCA*NUMV/2
     &              + (IVRTR-1)*NUMV*IOCCA + MSB
C
C             CONTRACT OVER ALL IBAS IN BLOCK A
              DO IBAS=1,NBAS(1)
C
C               (RA|SB) = (-A|SB) + (+A|SB)
                RASB(MRASB) = RASB(MRASB)
     &                    + ASB1(IBAS,MASB,1)*DCONJG(COEF(NA1+IBAS,IR))
     &                    + ASB1(IBAS,MASB,2)*DCONJG(COEF(NA2+IBAS,IR))
C
              ENDDO
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NUMO*NUMV*(NUMO+1)*NUMV*NBAS(1)/2 
C                                              ~
C     FOURTH CONTRACTION (SWAP): (JA;T|SB) -> (RA|SB)
      IF(G2INT.EQ.'COULM'.AND.IQ1.EQ.IQ2) GOTO 3003
      IF(G2INT.EQ.'GAUNT'.OR.G2INT.EQ.'BREIT') GOTO 3003
C
C     LOOP OVER OCCUPIED STATES IOCCA AND VIRTUAL STATES IVRTR
      DO IOCCA=1,NUMO
        DO IVRTR=1,NUMV
C
C         FOCK MATRIX ADDRESS FOR IVRTR
          IR = MINV-1+IVRTR+NSKP
C
C         LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
          DO IOCCB=1,IOCCA
            DO IVRTS=1,NUMV
C
C             LIST ADDRESS FOR THIS IVRTS,IOCCB COMBINATION
              MSB = (IOCCB-1)*NUMV+IVRTS
C
C             LIST ADDRESS FOR THIS IOCCA AND THE ABOVE MSB
              MASB = (IOCCA-1)*NUMV*NUMO + MSB
C
C             LIST ADDRESS FOR THIS IVRTR,IOCCA AND THE ABOVE MSB
              MRASB = (IOCCA-1)*NUMV*IOCCA*NUMV/2
     &              + (IVRTR-1)*NUMV*IOCCA + MSB
C
C             CONTRACT OVER ALL JBAS IN BLOCK B
              DO JBAS=1,NBAS(2)
C
C               (PA|SB) = (+A|SB) + (-A|SB)
                RASB(MRASB) = RASB(MRASB)
     &                    + ASB2(JBAS,MASB,1)*DCONJG(COEF(NB1+JBAS,IR))
     &                    + ASB2(JBAS,MASB,2)*DCONJG(COEF(NB2+JBAS,IR))
C
              ENDDO
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      INEST = INEST + NUMO*NUMV*(NUMO+1)*NUMV*NBAS(2)/2 
C
C     SKIP POINT FOR IQ1=IQ2
3003  CONTINUE
C
      CALL CPU_TIME(T2)
      TCN4 = TCN4+T2-T1
C
C     SKIP POINT FOR IQ1.LT.IQ2
3001  CONTINUE
C
C     ALL CONTRIBUTIONS FROM THIS CLASS (A,B) BLOCK NOW ACCOUNTED FOR
3000  CONTINUE
2000  CONTINUE
1000  CONTINUE
C
C
C**********************************************************************C
C     CALCULATE SECOND-ORDER PAIR CORRELATION ENERGY                   C
C**********************************************************************C
C
C     WRITE RASB ARRAY TO EXTERNAL FILE IF PROMPTED
10    FORMAT(/1X,72('=')/24(' '),'Writing data file in MBPT2'
     &       /1X,72('-')/1X,' Operator: ',A/1X,'File name: ',A
     &       /1X,'File size: ',F7.3,A,26X,'# words: ',I16/1X,72('='))
      IF(IOSAV.EQ.'WRTOUT'.OR.IOSAV.EQ.'SILENT') THEN
c       RASBOUT = TRIM(OUTFL)//'_RASB_'//G2INT//'500.dat'
        RASBOUT = TRIM(OUTFL)//'_RASB_'//G2INT//'.dat'
        NRASB   = (NUMO+1)*NUMO*NUMV*NUMV/2
        IF(IOSAV.NE.'SILENT') THEN
          SRASB = DFLOAT(NRASB)*16.0D-6
          IF(SRASB.LT.1.0D0) THEN
            WRITE(6,10) G2INT,TRIM(RASBOUT),1.0D+3*SRASB,' KB',NRASB
            WRITE(7,10) G2INT,TRIM(RASBOUT),1.0D+3*SRASB,' KB',NRASB
          ELSEIF(SRASB.LT.1.0D3) THEN
            WRITE(6,10) G2INT,TRIM(RASBOUT),1.0D+0*SRASB,' MB',NRASB
            WRITE(7,10) G2INT,TRIM(RASBOUT),1.0D+0*SRASB,' MB',NRASB
          ELSE
            WRITE(6,10) G2INT,TRIM(RASBOUT),1.0D-3*SRASB,' GB',NRASB
            WRITE(7,10) G2INT,TRIM(RASBOUT),1.0D-3*SRASB,' GB',NRASB
          ENDIF
        ENDIF
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
80    CONTINUE
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
      IF(IOSAV.EQ.'SILENT') RETURN
C
C     MBPT2 PAIRWISE SUMMARY
20    FORMAT(1X,A,11X,A,11X,A,11X,A,12X,A)
21    FORMAT(' (',I2,',',I2,')',3X,F13.7,5X,F11.7,5X,F11.7,4X,F13.7)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',25),'MBPT2 pairwise summary'
      WRITE(7, *) REPEAT(' ',25),'MBPT2 pairwise summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) '( a, b)','E1(G)','E2(J)','E2(K)','E2(G)'
      WRITE(7,20) '( a, b)','E1(G)','E2(J)','E2(K)','E2(G)'
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
C     MBPT2 SINGLE-PARTICLE SUMMARY
30    FORMAT(1X,A,11X,A,11X,A,11X,A,11X,A)
31    FORMAT('  ',I2,'    ',3X,F13.7,5X,F11.7,5X,F11.7,4X,F13.7)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',21),'MBPT2 single particle summary'
      WRITE(7, *) REPEAT(' ',21),'MBPT2 single particle summary'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,30) '  a    ','E1(G)','E2(J)','E2(K)',' E2(G)'
      WRITE(7,30) '  a    ','E1(G)','E2(J)','E2(K)',' E2(G)'
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
C     MBPT2 TOTAL SECOND ORDER CORRELATION ENERGY
32    FORMAT(' total  ',3X,F13.7,5X,F11.7,5X,F11.7,4X,F13.7)

      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',17),'MBPT2 second order correlation energy'
      WRITE(7, *) REPEAT(' ',17),'MBPT2 second order correlation energy'
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

40    FORMAT(1X,A,24X,A)
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',26),'MBPT2 labour analysis'
      WRITE(7, *) REPEAT(' ',26),'MBPT2 labour analysis'
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
      WRITE(6,40) 'ERI construction - (IJ|KL)      ',HMS(TERI)
      WRITE(7,40) 'ERI construction - (IJ|KL)      ',HMS(TERI)
      WRITE(6,40) '1st contraction  - (IJ|KB)      ',HMS(TCN1)
      WRITE(7,40) '1st contraction  - (IJ|KB)      ',HMS(TCN1)
      WRITE(6,40) '2nd contraction  - (IJ|SB)      ',HMS(TCN2)
      WRITE(7,40) '2nd contraction  - (IJ|SB)      ',HMS(TCN2)
      WRITE(6,40) '3rd contraction  - (IA|SB)      ',HMS(TCN3)
      WRITE(7,40) '3rd contraction  - (IA|SB)      ',HMS(TCN3)
      WRITE(6,40) '4th contraction  - (RA|SB)      ',HMS(TCN4)
      WRITE(7,40) '4th contraction  - (RA|SB)      ',HMS(TCN4)
      WRITE(6,40) 'Virtual orbital sum - E2(A,B)   ',HMS(TSUM)
      WRITE(7,40) 'Virtual orbital sum - E2(A,B)   ',HMS(TSUM)
      WRITE(6,40) 'Other                           ',HMS(TOTH)
      WRITE(7,40) 'Other                           ',HMS(TOTH)
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,40) 'Total MBPT2 time                ',HMS(TTOT)
      WRITE(7,40) 'Total MBPT2 time                ',HMS(TTOT)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
      RETURN
      END
C
C
      SUBROUTINE MBPTNPR(MINO,NUMO,MINV,NUMV,G2INT1,G2INT2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   MM       MM BBBBBBB  PPPPPPP TTTTTTTT NN    NN PPPPPPP  RRRRRRR    C
C   MMM     MMM BB    BB PP    PP   TT    NNN   NN PP    PP RR    RR   C
C   MMMM   MMMM BB    BB PP    PP   TT    NNNN  NN PP    PP RR    RR   C
C   MM MM MM MM BBBBBBB  PP    PP   TT    NN NN NN PP    PP RR    RR   C
C   MM  MMM  MM BB    BB PPPPPPP    TT    NN  NNNN PPPPPPP  RRRRRRR    C
C   MM   M   MM BB    BB PP         TT    NN   NNN PP       RR    RR   C
C   MM       MM BBBBBBB  PP         TT    NN    NN PP       RR    RR   C
C                                                                      C
C -------------------------------------------------------------------- C
C  MBPTNPR READS CALCULATED MBPT1 AND MBPT2 RESULTS UNDER THE COULOMB  C
C  AND BREIT/GAUNT INTERACTIONS (OR ANY OTHER TWO-BODY OPERATORS), AND C
C  THEN CALCULATES DIRECT/EXCHANGE CONTRIBUTIONS FOR THEIR SUM.        C
C  MBPT1 ENERGY COUNTERS ARE READ IN FROM FILE (CAN OPTIONALLY READ IN C
C  (AC|BD) MATRIX ELEMENTS IN BOTH CASES IF YOU HAVE A USE FOR THEM)   C
C  BUT MBPT2 ENERGY COUNTERS ARE RE-CALCULATED, SINCE WE WILL NEED TO  C
C  EVALUATE THE CROSS-TERMS ANYWAY.                                    C
C -------------------------------------------------------------------- C
C INPUT:                                                               C
C ▶ MINO   - LOWEST OCCUPIED STATE TO ACCOUNT FOR. (FULL: 1)           C
C ▶ NUMO   - NUMBER OF OCCUPIED STATES TO ACCOUNT FOR. (FULL: NOCC)    C
C ▶ MINV   - LOWEST VIRTUAL STATE TO ACCOUNT FOR. (FULL: NOCC+1)       C
C ▶ NUMV   - NUMBER OF VIRTUAL STATES TO ACCOUNT FOR. (FULL: NVRT)     C
C ▶ G2INT1 - NAME OF 1ST 2-BODY OPERATOR ('COULM', 'GAUNT' OR 'BREIT').C
C ▶ G2INT2 - NAME OF 2ND 2-BODY OPERATOR ('COULM', 'GAUNT' OR 'BREIT').C
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
      DIMENSION EAB1(NUMO,NUMO,8),EA1(NUMO,8)
      DIMENSION EAB2(NUMO,NUMO,5),EA2(NUMO,5)
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
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EWDR,EWXC,EANM,ESLF,EUEH,
     &            EWKR,EKSB
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
C
C     WARNINGS BASED ON INVALID HMLT VS. G2INT COMBINATIONS
10    FORMAT(1X,'In MBPTNPR: HMLT = ',A,' but G2INT1 = ',A,
     &                                          ' and G2INT2 = ',A,'.')
      IF(HMLT.EQ.'NORL'.OR.HMLT.EQ.'DHFR'.OR.HMLT.EQ.'DHFP') THEN
        IF(G2INT1.EQ.'GAUNT'.OR.G2INT1.EQ.'BREIT'.OR.
     &     G2INT2.EQ.'GAUNT'.OR.G2INT2.EQ.'BREIT') THEN
          WRITE(6,10) HMLT,G2INT1,G2INT2
          WRITE(7,10) HMLT,G2INT1,G2INT2
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
C     INITIALISE ENERGY COUNTER ARRAYS
      DO IOCCA=1,NUMO
        DO IOCCB=1,IOCCA
          DO N=1,8
            EAB1(IOCCA,IOCCB,N) = 0.0D0
          ENDDO
          DO N=1,5
            EAB2(IOCCA,IOCCB,N) = 0.0D0
          ENDDO
        ENDDO
        DO N=1,8
          EA1(IOCCA,N) = 0.0D0
        ENDDO
        DO N=1,5
          EA2(IOCCA,N) = 0.0D0
        ENDDO
      ENDDO
C
C     CHECK FOR EXISTENCE OF MBPT1 DATA GIVEN G2INT1
      EAB1OUT = TRIM(OUTFL)//'_EAB1_'//G2INT1//'.dat'
      INQUIRE(FILE=TRIM(EAB1OUT),EXIST=FILETHERE)
C
C     IF NO FILE EXISTS, NEED MBPT1 (CAN USE ANY IOSAV HERE)
      IF(.NOT.FILETHERE) THEN
        WRITE(6, *) 'In MBPTNPR: calling MBPT1 for G2INT = ',G2INT1
        WRITE(7, *) 'In MBPTNPR: calling MBPT1 for G2INT = ',G2INT1
C       CALL MBPT1(MINO,NUMO,G2INT1,'SILENT')
      ENDIF
C
C     IMPORT MBPT1 ENERGY RESULTS FOR G2INT1
      OPEN(UNIT=8,FILE=TRIM(EAB1OUT),STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO IOCCA=1,NUMO
        DO IOCCB=1,NUMO
          READ(8, *) T1,T2,(EAB1(IOCCA,IOCCB,N),N=1,4)
        ENDDO
      ENDDO
      CLOSE(UNIT=8)
C
C     CHECK FOR EXISTENCE OF MBPT1 DATA GIVEN G2INT2
      EAB1OUT = TRIM(OUTFL)//'_EAB1_'//G2INT2//'.dat'
      INQUIRE(FILE=TRIM(EAB1OUT),EXIST=FILETHERE)
C
C     IF NO FILE EXISTS, NEED MBPT1 (CAN USE ANY IOSAV HERE)
      IF(.NOT.FILETHERE) THEN
        WRITE(6, *) 'In MBPTNPR: calling MBPT1 for G2INT = ',G2INT2
        WRITE(7, *) 'In MBPTNPR: calling MBPT1 for G2INT = ',G2INT2
C       CALL MBPT1(MINO,NUMO,G2INT2,'SILENT')
      ENDIF
C
C     IMPORT MBPT1 ENERGY RESULTS FOR G2INT2
      OPEN(UNIT=8,FILE=TRIM(EAB1OUT),STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO IOCCA=1,NUMO
        DO IOCCB=1,NUMO
          READ(8, *) T1,T2,T3,(EAB1(IOCCA,IOCCB,N),N=5,7)
        ENDDO
      ENDDO
      CLOSE(UNIT=8)
C
C     TOTALS FOR EAB1
      DO IOCCA=1,NUMO
        DO IOCCB=1,NUMO
          EAB1(IOCCA,IOCCB,8) = EAB1(IOCCA,IOCCB,1)
     &                        + EAB1(IOCCA,IOCCB,4)+EAB1(IOCCA,IOCCB,7)
        ENDDO
      ENDDO
C
C     INITIALISE THE ARRAY OF (AR|BS) VALUES
      M = 0
      DO IVRTR=1,NUMV
        DO IOCCA=1,NUMO
          DO IOCCB=1,IOCCA
            DO IVRTS=1,NUMV
              M = M+1
              RASB1(M) = DCMPLX(0.0D0,0.0D0)
              RASB2(M) = DCMPLX(0.0D0,0.0D0)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     CHECK FOR EXISTENCE OF MBPT1 DATA GIVEN G2INT1
c     RASBOUT = TRIM(OUTFL)//'_RASB_'//G2INT1//'500.dat'
      RASBOUT = TRIM(OUTFL)//'_RASB_'//G2INT1//'.dat'
      INQUIRE(FILE=TRIM(RASBOUT),EXIST=FILETHERE)
C
C     IF NO FILE EXISTS, NEED MBPT2
      IF(.NOT.FILETHERE) THEN
        WRITE(6, *) 'In MBPTNPR: calling MBPT2 for G2INT = ',G2INT1
        WRITE(7, *) 'In MBPTNPR: calling MBPT2 for G2INT = ',G2INT1
        CALL MBPT2(MINO,NUMO,MINV,NUMV,G2INT1,'SILENT')
      ENDIF
C
C     IMPORT (AR|BS) MATRIX ELEMENTS FOR G2INT1
      OPEN(UNIT=8,FILE=TRIM(RASBOUT),STATUS='OLD')
      REWIND(UNIT=8)
C     READ MINO AND NUMO
      READ(8, *) MINOTMP,NUMOTMP,MINVTMP,NUMVTMP
C     IF DIFFERENT TO INPUT, WRITE TO LOG AND UPDATE
      IF(MINOTMP.NE.MINO) THEN
        WRITE(6, *) 'In MBPTNPR: MINO clash.',MINO,MINOTMP,G2INT1
        WRITE(7, *) 'In MBPTNPR: MINO clash.',MINO,MINOTMP,G2INT1
      ENDIF
      IF(NUMOTMP.NE.NUMO) THEN
        WRITE(6, *) 'In MBPTNPR: NUMO clash.',NUMO,NUMOTMP,G2INT1
        WRITE(7, *) 'In MBPTNPR: NUMO clash.',NUMO,NUMOTMP,G2INT1
      ENDIF
      IF(MINVTMP.NE.MINV) THEN
        WRITE(6, *) 'In MBPTNPR: MINV clash.',MINV,MINVTMP,G2INT1
        WRITE(7, *) 'In MBPTNPR: MINV clash.',MINV,MINVTMP,G2INT1
      ENDIF
      IF(NUMVTMP.NE.NUMV) THEN
        WRITE(6, *) 'In MBPTNPR: NUMV clash.',NUMV,NUMVTMP,G2INT1
        WRITE(7, *) 'In MBPTNPR: NUMV clash.',NUMV,NUMVTMP,G2INT1
      ENDIF
      MINO = MINOTMP
      NUMO = NUMOTMP
      MINV = MINVTMP
      NUMV = NUMVTMP
C     LOOP OVER OCCUPIED STATES IOCCA AND VIRTUAL STATES IVRTR
      DO IOCCA=1,NUMO
        DO IVRTR=1,NUMV
C         LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
          DO IOCCB=1,IOCCA
            DO IVRTS=1,NUMV
C             MAIN LIST ADDRESS FOR THIS IOCCA,IVRTR,IOCCB,IVRTS
              MRASB = (IOCCA-1)*NUMV*IOCCA*NUMV/2
     &              + (IVRTR-1)*NUMV*IOCCA + (IOCCB-1)*NUMV + IVRTS
              READ(8, *) RASB1(MRASB)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CLOSE(UNIT=8)
C
C     CHECK FOR EXISTENCE OF MBPT1 DATA GIVEN G2INT2
c     RASBOUT = TRIM(OUTFL)//'_RASB_'//G2INT2//'500.dat'
      RASBOUT = TRIM(OUTFL)//'_RASB_'//G2INT2//'.dat'
      INQUIRE(FILE=TRIM(RASBOUT),EXIST=FILETHERE)
C
C     IF NO FILE EXISTS, NEED MBPT2
      IF(.NOT.FILETHERE) THEN
        WRITE(6, *) 'In MBPTNPR: calling MBPT2 for G2INT = ',G2INT2
        WRITE(7, *) 'In MBPTNPR: calling MBPT2 for G2INT = ',G2INT2
        CALL MBPT2(MINO,NUMO,MINV,NUMV,G2INT2,'SILENT')
      ENDIF
C
C     IMPORT (AR|BS) MATRIX ELEMENTS FOR G2INT2
      OPEN(UNIT=8,FILE=TRIM(RASBOUT),STATUS='OLD')
      REWIND(UNIT=8)
C     READ MINO AND NUMO
      READ(8, *) MINOTMP,NUMOTMP,MINVTMP,NUMVTMP
C     IF DIFFERENT TO INPUT, WRITE TO LOG AND UPDATE
      IF(MINOTMP.NE.MINO) THEN
        WRITE(6, *) 'In MBPTNPR: MINO clash.',MINO,MINOTMP,G2INT2
        WRITE(7, *) 'In MBPTNPR: MINO clash.',MINO,MINOTMP,G2INT2
      ENDIF
      IF(NUMOTMP.NE.NUMO) THEN
        WRITE(6, *) 'In MBPTNPR: NUMO clash.',NUMO,NUMOTMP,G2INT2
        WRITE(7, *) 'In MBPTNPR: NUMO clash.',NUMO,NUMOTMP,G2INT2
      ENDIF
      IF(MINVTMP.NE.MINV) THEN
        WRITE(6, *) 'In MBPTNPR: MINV clash.',MINV,MINVTMP,G2INT2
        WRITE(7, *) 'In MBPTNPR: MINV clash.',MINV,MINVTMP,G2INT2
      ENDIF
      IF(NUMVTMP.NE.NUMV) THEN
        WRITE(6, *) 'In MBPTNPR: NUMV clash.',NUMV,NUMVTMP,G2INT2
        WRITE(7, *) 'In MBPTNPR: NUMV clash.',NUMV,NUMVTMP,G2INT2
      ENDIF
      MINO = MINOTMP
      NUMO = NUMOTMP
      MINV = MINVTMP
      NUMV = NUMVTMP
C     LOOP OVER OCCUPIED STATES IOCCA AND VIRTUAL STATES IVRTR
      DO IOCCA=1,NUMO
        DO IVRTR=1,NUMV
C         LOOP OVER OCCUPIED STATES IOCCB AND VIRTUAL STATES IVRTS
          DO IOCCB=1,IOCCA
            DO IVRTS=1,NUMV
C             MAIN LIST ADDRESS FOR THIS IOCCA,IVRTR,IOCCB,IVRTS
              MRASB = (IOCCA-1)*NUMV*IOCCA*NUMV/2
     &              + (IVRTR-1)*NUMV*IOCCA + (IOCCB-1)*NUMV + IVRTS
              READ(8, *) RASB2(MRASB)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CLOSE(UNIT=8)
C
C**********************************************************************C
C     CALCULATE SECOND-ORDER PAIR CORRELATION ENERGY                   C
C**********************************************************************C
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
C             DENOMINATOR FOR BOTH CONTRIBUTIONS
              EABRS = EIGN(IA) + EIGN(IB) - EIGN(IR) - EIGN(IS)
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
              RNUM11D = DREAL(RASB1(MRASB)*DCONJG(RASB1(MRASB)))
              RNUM11X =-DREAL(RASB1(MRASB)*DCONJG(RASB1(MSARB)))
              V11     = RNUM11D+RNUM11X
C
              RNUM12D = DREAL(RASB1(MRASB)*DCONJG(RASB2(MRASB)))
              RNUM12X =-DREAL(RASB1(MRASB)*DCONJG(RASB2(MSARB)))
              V12     = RNUM12D+RNUM12X
C
              RNUM21D = DREAL(RASB2(MRASB)*DCONJG(RASB1(MRASB)))
              RNUM21X =-DREAL(RASB2(MRASB)*DCONJG(RASB1(MSARB)))
              V21     = RNUM21D+RNUM21X
C
              RNUM22D = DREAL(RASB2(MRASB)*DCONJG(RASB2(MRASB)))
              RNUM22X =-DREAL(RASB2(MRASB)*DCONJG(RASB2(MSARB)))
              V22     = RNUM22D+RNUM22X
C
C             ADD TO OCCUPIED COUNTER BINS
              EAB2(IOCCA,IOCCB,1) = EAB2(IOCCA,IOCCB,1) + V11/EABRS
              EAB2(IOCCA,IOCCB,2) = EAB2(IOCCA,IOCCB,2) + V12/EABRS
              EAB2(IOCCA,IOCCB,3) = EAB2(IOCCA,IOCCB,3) + V21/EABRS
              EAB2(IOCCA,IOCCB,4) = EAB2(IOCCA,IOCCB,4) + V22/EABRS
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     THE OTHER HALF OF THE ARRAY
      DO IOCCA=1,NUMO
        DO IOCCB=1,IOCCA-1
          DO N=1,4
            EAB2(IOCCB,IOCCA,N) = EAB2(IOCCA,IOCCB,N)
          ENDDO
        ENDDO
      ENDDO
C
C     COMBINED ENERGY
      DO IOCCA=1,NUMO
        DO IOCCB=1,NUMO
          ETMP = 0.0D0
          DO N=1,4
            ETMP = ETMP + EAB2(IOCCA,IOCCB,N)
          ENDDO
          EAB2(IOCCA,IOCCB,5) = ETMP
        ENDDO
      ENDDO
C
C     CALCULATE TOTALS
      E11 = 0.0D0
      E12 = 0.0D0
      E21 = 0.0D0
      E22 = 0.0D0
      DO IOCCA=1,NUMO
        DO IOCCB=1,IOCCA
          E11 = E11 + 0.5D0*EAB2(IOCCA,IOCCB,1)
          E12 = E12 + 0.5D0*EAB2(IOCCA,IOCCB,2)
          E21 = E21 + 0.5D0*EAB2(IOCCA,IOCCB,3)
          E22 = E22 + 0.5D0*EAB2(IOCCA,IOCCB,4)
        ENDDO
      ENDDO
      ETT = E11+E12+E21+E22
C
C**********************************************************************C
C     CALCULATE FIRST- AND SECOND-ORDER SINGLE ORBITAL ENERGY          C
C**********************************************************************C
C
C     FOR EACH IOCCA, SUM OVER THE IOCCB CONTRIBUTIONS
      DO IOCCA=1,NUMO
        DO N=1,8
          EA1(IOCCA,N) = 0.0D0
          DO IOCCB=1,NUMO
            EA1(IOCCA,N) = EA1(IOCCA,N) + EAB1(IOCCA,IOCCB,N)
          ENDDO
        ENDDO
      ENDDO
C
C     CALCULATE TOTALS
      E1VK = 0.0D0
      E1GD = 0.0D0
      E2GX = 0.0D0
      E2BD = 0.0D0
      E2BX = 0.0D0
      DO IOCCA=1,NUMO
        DO IOCCB=1,IOCCA
          E1VK = E1VK + 0.5D0*EAB1(IOCCA,IOCCB,1)
          E1GD = E1GD + 0.5D0*EAB1(IOCCA,IOCCB,2)
          E1GX = E1GX + 0.5D0*EAB1(IOCCA,IOCCB,3)
          E1BD = E1BD + 0.5D0*EAB1(IOCCA,IOCCB,5)
          E1BX = E1BX + 0.5D0*EAB1(IOCCA,IOCCB,6)
        ENDDO
      ENDDO
      E1GT = E1GD+E1GX
      E1BT = E1BD+E1BX
      E1T  = 2.0D0*(E1VK+E1GT+E1BT)
C
C     FOR EACH IOCCA, SUM OVER THE IOCCB CONTRIBUTIONS
      DO IOCCA=1,NUMO
        DO N=1,5
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
C     MBPT1(PR) PAIRWISE SUMMARY
20    FORMAT(1X,A,7X,'E1(V+K)',9X,'E1(C)',9X,'E1(B)',7X,'E1(tot)',7X,A)
21    FORMAT(' (',I2,',',I2,')',5(1X,F13.7))
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',25),'MBPT1(PR) pairwise summary'
      WRITE(7, *) REPEAT(' ',25),'MBPT1(PR) pairwise summary'
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
      WRITE(6,20) '( a, b)','E1(a,b)'
      WRITE(7,20) '( a, b)','E1(a,b)'
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
      DO IA=1,NUMO
        IANUM = MINO-1+IA+NSKP
        DO IB=1,IA
          IBNUM = MINO-1+IB+NSKP
          IF(IA.EQ.IB) THEN
            E1 = EIGN(IANUM)
          ELSE
            E1 = 0.0D0
          ENDIF
          WRITE(6,21) IA,IB,EAB1(IA,IB,1),EAB1(IA,IB,4),
     &                      EAB1(IA,IB,7),EAB1(IA,IB,8),E1
          WRITE(7,21) IA,IB,EAB1(IA,IB,1),EAB1(IA,IB,4),
     &                      EAB1(IA,IB,7),EAB1(IA,IB,8),E1
        ENDDO
      ENDDO
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
C
C     MBPT1(PR) SINGLE-PARTICLE SUMMARY
22    FORMAT('  ',I2,'    ',5(1X,F13.7))
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',22),'MBPT1(PR) single particle summary'
      WRITE(7, *) REPEAT(' ',22),'MBPT1(PR) single particle summary'
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
      WRITE(6,20) '  a    ','  E1(a)'
      WRITE(7,20) '  a    ','  E1(a)'
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
      DO IA=1,NUMO
        IANUM = MINO-1+IA+NSKP
        E1 = EIGN(IANUM)
        WRITE(6,22) IA,EA1(IA,1),EA1(IA,4),EA1(IA,7),EA1(IA,8),E1
        WRITE(7,22) IA,EA1(IA,1),EA1(IA,4),EA1(IA,7),EA1(IA,8),E1
      ENDDO
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
C
C     MBPT1(PR) TOTAL SECOND ORDER CORRELATION ENERGY
23    FORMAT(' total  ',5(1X,F13.7))
C
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',24),'MBPT1(PR) correlation energy'
      WRITE(7, *) REPEAT(' ',24),'MBPT1(PR) correlation energy'
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
      WRITE(6,23) E1VK,E1GT,E1BT,E1T,ETOT
      WRITE(7,23) E1VK,E1GT,E1BT,E1T,ETOT
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
C
C     MBPT2(PR) PAIRWISE SUMMARY
30    FORMAT(1X,A,8X,'E2(CC)',8X,'E2(CB)',8X,'E2(BC)',8X,'E2(BB)',7X,A)
31    FORMAT(' (',I2,',',I2,')',5(1X,F13.7))
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',25),'MBPT2(PR) pairwise summary'
      WRITE(7, *) REPEAT(' ',25),'MBPT2(PR) pairwise summary'
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
      WRITE(6,30) '( a, b)','E2(a,b)'
      WRITE(7,30) '( a, b)','E2(a,b)'
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
      DO IA=1,NUMO
        IANUM = IA+MINO-1
        DO IB=1,IA
          IBNUM = IB+MINO-1
          WRITE(6,31) IANUM,IBNUM,(EAB2(IA,IB,N),N=1,5)
          WRITE(7,31) IANUM,IBNUM,(EAB2(IA,IB,N),N=1,5)
        ENDDO
      ENDDO
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
C
C     MBPT2(PR) SINGLE-PARTICLE SUMMARY
32    FORMAT('  ',I2,'    ',5(1X,F13.7))
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',22),'MBPT2(PR) single particle summary'
      WRITE(7, *) REPEAT(' ',22),'MBPT2(PR) single particle summary'
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
      WRITE(6,30) '  a    ','  E2(a)'
      WRITE(7,30) '  a    ','  E2(a)'
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
      DO IA=1,NUMO
        IANUM = IA+MINO-1
        WRITE(6,32) IANUM,(EA2(IA,N),N=1,5)
        WRITE(7,32) IANUM,(EA2(IA,N),N=1,5)
      ENDDO
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
C
C     MBPT2(PR) TOTAL SECOND ORDER CORRELATION ENERGY
33    FORMAT(' total  ',5(1X,F13.7))
C
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',24),'MBPT2(PR) correlation energy'
      WRITE(7, *) REPEAT(' ',24),'MBPT2(PR) correlation energy'
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
      WRITE(6,33) E11,E12,E21,E22,ETT
      WRITE(7,33) E11,E12,E21,E22,ETT
      WRITE(6, *) REPEAT('-',77)
      WRITE(7, *) REPEAT('-',77)
C
C     TIME TAKEN IN THIS ROUTINE
      CALL CPU_TIME(TFIN)
      TTOT = TFIN - TBEG
      TOTH = TTOT - (TERI + TCN1 + TCN2 + TCN3 + TCN4 + TSUM)
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
C   [B] FUNFM: LIST OF BOYS INTEGRALS FOR USE IN RMAKE.                C
C   [C] BOYSGEN: OUTPUT DATA FILE WITH FAMILY OF BOYS FUNCTIONS.       C
C   [D] HGTFS: BATCH OF HGTF AMPLITUDES EVALUATED AT (X,Y,Z).          C
C   [E] HERMITE: EVALUATION OF H_I (P,X) BY RECURRENCE.                C
C   [F] HTFRMS: BATCH OF HGTF FOURIER TRANSFORMS AT (QX,QY,QZ).        C
C**********************************************************************C
C
C
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
C  FINITE SUM REPRESENTATION OF A MULTI-CENTRE GAUSSIAN OVERLAP.       C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION FS(MB2,ML4),APH(MB2),QP(MB2,3),RC(MB2,MRC),RC2(MB2,MRC)
      DIMENSION F1(MB2,ML4),F2(MB2,ML4),F3(MB2,ML4),F4(MB2,ML4)
      DIMENSION X1(MB2),X2(MB2),X3(MB2),X4(MB2),
     &          I1(MB2),I2(MB2),I3(MB2),I4(MB2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
C     VALUES WHICH DETERMINE EVALUATION METHOD
      EPSZR = 1.0D-11
      EPSPL = 1.7D+01
      EPSAS = 3.0D+01
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
C     INITIALISE COUNTERS FOR EVALUATION METHODS IN BATCH
      N1 = 0
      N2 = 0
      N3 = 0
      N4 = 0
C
C     FOR EACH PAIR OF BASIS FUNCTIONS (EXPONENTS EI AND EJ IN 'M'),
C     DETERMINE THE BEST WAY TO EVALUATE THE BOYS FUNCTION
      DO M=1,MAXM
C
        X = QP(M,1)*QP(M,1) + QP(M,2)*QP(M,2) + QP(M,3)*QP(M,3)
        X = X*APH(M)
C
C       CASE 1: IF X=0.0D0 USE EXACT FORMULA
        IF(X.LE.EPSZR) THEN
C
          N1     = N1+1
          X1(N1) = X
          I1(N1) = M
C
C       CASE 2: IF X<EPSPL USE A POLYNOMIAL EXPANSION
        ELSEIF(X.GT.EPSZR.AND.X.LE.EPSPL) THEN
C
          N2     = N2+1
          X2(N2) = X
          I2(N2) = M
C
C       CASE 3: IF X<EPSAS USE ASYMPTOIC FORMULA WITH EXPONENTIAL
        ELSEIF(X.GT.EPSPL.AND.X.LE.EPSAS) THEN
C
          N3     = N3+1
          X3(N3) = X
          I3(N3) = M
C
C       CASE 4: IF X>EPSAS USE ASYMPTOIC FORMULA WITHOUT EXPONENTIAL
        ELSE
C
          N4     = N4+1
          X4(N4) = X
          I4(N4) = M
C
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
C     CASE 2: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=EPSPL.
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
C     CASE 3: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=EPSAS.
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
C     CASE 4: ARGUMENT OF THE BOYS FUNCTION IS LARGER THAN X=EPSAS.
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
C     WITH THE FULL SET OF BOYS' INTEGRALS WE NOW APPLY RECURRENCE     C
C     RELATIONS TO F_N (X) AND EVALUATE THE R-INTEGRALS.               C
C**********************************************************************C
C
C     CONSTRUCT TOP LEVEL (FOR MAXIMUM LAM VALUE)
      DO M=1,MAXM
        RC(M,1)= (-2.0D0*APH(M))**LAM
        RC(M,1) = RC(M,1)*FS(M,LAM+1)
      ENDDO
C
C     MINIMUM LEVEL ILV BASED ON LAM VALUE
      IF(MOD(LAM,2).EQ.0) THEN
       ITUVMIN = 1
      ELSE
       ITUVMIN = 2
      ENDIF
C
C     INITIALISE ITUV COUNTER (RELATES TO # CARTESIAN INDICES FOR lam)
      ITUV=-1
C
C     MAIN LOOP: LEVEL 'ILV' STARTING AT LAM-1 AND WORKING BACKWARDS
      DO ILV=LAM-1,ITUVMIN,-2
C
C       UPDATE ITUV COUNTER
        ITUV = ITUV+1
C
C       LOOP OVER ALL UNIQUE (IT,IU,IV)
        DO IT=0,ITUV
          RIT = DFLOAT(IT)
C
          DO IU=0,ITUV-IT
            RIU = DFLOAT(IU)
C
            DO IV=0,ITUV-IT-IU
              RIV = DFLOAT(IV)
C
C             R-INTEGRAL CARTESIAN DESTINATION ADDRESSES
              N1 = IABC(IT+1,IU  ,IV  )
              N2 = IABC(IT  ,IU+1,IV  )
              N3 = IABC(IT  ,IU  ,IV+1)
C
C             RECURRENCE RELATIONS DIFFERENT IF ANY (IT,IU,IV) ARE ZERO

C             CASE (IT,IU,IV)
              IF(IT.NE.0.AND.IU.NE.0.AND.IV.NE.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M1 = IABC(IT-1,IU  ,IV  )
                M2 = IABC(IT  ,IU-1,IV  )
                M3 = IABC(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE (IT,IU, 0)
              ELSEIF(IT.NE.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M1 = IABC(IT-1,IU  ,IV  )
                M2 = IABC(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
               ENDDO
C
C             CASE (IT, 0,IV)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M1 = IABC(IT-1,IU  ,IV  )
                M3 = IABC(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE (IT, 0, 0)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M1 = IABC(IT-1,IU  ,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE ( 0,IU,IV)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.NE.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M2 = IABC(IT  ,IU-1,IV  )
                M3 = IABC(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE ( 0,IU, 0)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M2 = IABC(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE ( 0, 0,IV)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M3 = IABC(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE ( 0, 0, 0)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
                ENDDO
C
C             ALL POSSIBILITIES ACCOUNTED FOR -- END LOOP
              ENDIF
C
C           END LOOPS OVER (IT,IU,IV) ADDRESSES FOR THIS ILV
            ENDDO
C
          ENDDO
C
        ENDDO
C
C       ADD IN (IT=0,IU=0,IV=0) CASE
        DO M=1,MAXM
          RC2(M,1) = (-2.0D0*APH(M))**ILV
          RC2(M,1) = RC2(M,1)*FS(M,ILV+1)
        ENDDO
C
C       UPDATE ITUV COUNTER
        ITUV = ITUV+1
C
C       LOOP OVER ALL UNIQUE (IT,IU,IV)
        DO IT=0,ITUV
          RIT = DFLOAT(IT)
C
          DO IU=0,ITUV-IT
            RIU = DFLOAT(IU)
C
            DO IV=0,ITUV-IT-IU
              RIV = DFLOAT(IV)
C
C             R-INTEGRAL CARTESIAN DESTINATION ADDRESSES
              N1 = IABC(IT+1,IU  ,IV  )
              N2 = IABC(IT  ,IU+1,IV  )
              N3 = IABC(IT  ,IU  ,IV+1)
C
C             RECURRENCE RELATIONS DIFFERENT IF ANY (IT,IU,IV) ARE ZERO

C             CASE (IT,IU,IV)
              IF(IT.NE.0.AND.IU.NE.0.AND.IV.NE.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M1 = IABC(IT-1,IU  ,IV  )
                M2 = IABC(IT  ,IU-1,IV  )
                M3 = IABC(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC(M,N1) =-QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                  RC(M,N2) =-QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) =-QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
C
C             CASE (IT,IU, 0)
              ELSEIF(IT.NE.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M1 = IABC(IT-1,IU  ,IV  )
                M2 = IABC(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC(M,N1) =-QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                  RC(M,N2) =-QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) =-QP(M,3)*RC2(M,K1)
                ENDDO
C
C             CASE (IT, 0,IV)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M1 = IABC(IT-1,IU  ,IV  )
                M3 = IABC(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC(M,N1) =-QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                  RC(M,N2) =-QP(M,2)*RC2(M,K1)
                  RC(M,N3) =-QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
C
C             CASE (IT, 0, 0)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M1 = IABC(IT-1,IU  ,IV  )
                DO M=1,MAXM
                 RC(M,N1) =-QP(M,1)*RC2(M,K1) + RIT*RC2(M,M1)
                 RC(M,N2) =-QP(M,2)*RC2(M,K1)
                 RC(M,N3) =-QP(M,3)*RC2(M,K1)
               ENDDO
C
C             CASE ( 0,IU,IV)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.NE.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M2 = IABC(IT  ,IU-1,IV  )
                M3 = IABC(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC(M,N1) =-QP(M,1)*RC2(M,K1)
                  RC(M,N2) =-QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) =-QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
C
C             CASE ( 0,IU, 0)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M2 = IABC(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC(M,N1) =-QP(M,1)*RC2(M,K1)
                  RC(M,N2) =-QP(M,2)*RC2(M,K1) + RIU*RC2(M,M2)
                  RC(M,N3) =-QP(M,3)*RC2(M,K1)
                ENDDO
C
C             CASE ( 0, 0,IV)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M3 = IABC(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC(M,N1) =-QP(M,1)*RC2(M,K1)
                  RC(M,N2) =-QP(M,2)*RC2(M,K1)
                  RC(M,N3) =-QP(M,3)*RC2(M,K1) + RIV*RC2(M,M3)
                ENDDO
C
C             CASE ( 0, 0, 0)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                DO M=1,MAXM
                  RC(M,N1) =-QP(M,1)*RC2(M,K1)
                  RC(M,N2) =-QP(M,2)*RC2(M,K1)
                  RC(M,N3) =-QP(M,3)*RC2(M,K1)
                ENDDO
C
C             ALL POSSIBILITIES ACCOUNTED FOR -- END LOOP
              ENDIF
C
C           END LOOPS OVER (IT,IU,IV) ADDRESSES FOR THIS ILV
            ENDDO
C
          ENDDO
C
        ENDDO
C
C       ADD IN (IT=0,IU=0,IV=0) CASE
C
        DO M=1,MAXM
          RC(M,1) = (-2.0D0*APH(M))**(ILV-1)
          RC(M,1) = RC(M,1)*FS(M,ILV)
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
C
          DO IU=0,ITUV-IT
            RIU=DFLOAT(IU)
C
            DO IV=0,ITUV-IT-IU
              RIV=DFLOAT(IV)
C
C             R-INTEGRAL CARTESIAN DESTINATION ADDRESSES
              N1 = IABC(IT+1,IU  ,IV  )
              N2 = IABC(IT  ,IU+1,IV  )
              N3 = IABC(IT  ,IU  ,IV+1)
C
C             RECURRENCE RELATIONS DIFFERENT IF ANY (IT,IU,IV) ARE ZERO

C             CASE (IT,IU,IV)
              IF(IT.NE.0.AND.IU.NE.0.AND.IV.NE.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M1 = IABC(IT-1,IU  ,IV  )
                M2 = IABC(IT  ,IU-1,IV  )
                M3 = IABC(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE (IT,IU, 0)
              ELSEIF(IT.NE.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M1 = IABC(IT-1,IU  ,IV  )
                M2 = IABC(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE (IT, 0,IV)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M1 = IABC(IT-1,IU  ,IV  )
                M3 = IABC(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE (IT, 0, 0)
              ELSEIF(IT.NE.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M1 = IABC(IT-1,IU  ,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1) + RIT*RC(M,M1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE ( 0,IU,IV)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.NE.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M2 = IABC(IT  ,IU-1,IV  )
                M3 = IABC(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE ( 0,IU, 0)
              ELSEIF(IT.EQ.0.AND.IU.NE.0.AND.IV.EQ.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M2 = IABC(IT  ,IU-1,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1) + RIU*RC(M,M2)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
                ENDDO
C
C             CASE ( 0, 0,IV)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.NE.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                M3 = IABC(IT  ,IU  ,IV-1)
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1) + RIV*RC(M,M3)
                ENDDO
C
C             CASE ( 0, 0, 0)
              ELSEIF(IT.EQ.0.AND.IU.EQ.0.AND.IV.EQ.0) THEN
C
                K1 = IABC(IT  ,IU  ,IV  )
                DO M=1,MAXM
                  RC2(M,N1) =-QP(M,1)*RC(M,K1)
                  RC2(M,N2) =-QP(M,2)*RC(M,K1)
                  RC2(M,N3) =-QP(M,3)*RC(M,K1)
               ENDDO
C
C             ALL POSSIBILITIES ACCOUNTED FOR -- END LOOP
              ENDIF
C
C           END LOOPS OVER (IT,IU,IV) ADDRESSES FOR THIS ILV
            ENDDO
C
          ENDDO
C
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
C  ▶ ITYPE = 1 - SPECIAL CASE X = 0.0D0.                               C
C  ▶ ITYPE = 2 - POWER SERIES AND REVERSE RECURRENCE.                  C
C                (60 TERMS WILL BE USED, SO USE MUST SUPPLY A VALUE    C
C                 APPROPRIATE TO THE MAX VALUE OF X IN BATCH).         C
C  ▶ ITYPE = 3 - ASYMPTOTIC EXPANSION AND FORWARD RECURRENCE.          C
C  ▶ ITYPE = 4 - ASYMPTOTIC EXPANSION AND FORWARD RECURRENCE.          C
C                ALL TERMS DEPENDING ON EXP(-X) ARE OMITTED TO AVOID   C
C                NUMERICAL UNDERFLOW PROBLEMS.                         C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION FM(MB2,ML4),T(MB2),TLAM(MB2),TT2(MB2),
     &          TXP(MB2),TRT(MB2)
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
      DATA A0,B0/4.994501191201870D-1,4.551838436668326D-1/
C
C**********************************************************************C
C     ITYPE = 1: SPECIAL CASE FOR T = 0.0D0                            C
C**********************************************************************C
C
      IF(ITYPE.EQ.1) THEN
C
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
C       INITIALISE THE POWER SERIES FOR M = LAM
        DO M=1,N
          TXP(M)      = DEXP(-T(M))
          TT2(M)      = 2.0D0*T(M)
          TLAM(M)     = 1.0D0
          FM(M,LAM+1) = 1.0D0
        ENDDO
C
C       LOOP OVER TERMS IN THE POWER SERIES
        DO K=1,60
          DLAM = DFLOAT(2*LAM+2*K+1)
          DO M=1,N
            TLAM(M)     = TLAM(M)*TT2(M)/DLAM
            FM(M,LAM+1) = FM(M,LAM+1) + TLAM(M)
          ENDDO
        ENDDO
C
C       RESCALE BY THE PREFACTOR
        DEN = DFLOAT(2*LAM+1)
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
C     ITYPE = 3: ASYMPTOTIC EXPANSION WITH EXPONENTIAL ARGUMENT.       C
C**********************************************************************C
C
      ELSEIF(ITYPE.EQ.3) THEN
C
C       INITIALISE THE ASYMPTOTIC EXPANSION
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
          FM(M,1) = 0.5D0*PI12/TRT(M) - TXP(M)*FM(M,1)
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
C       INITIALISE THE ASYMPTOTIC EXPANSION
        DO M=1,N
          TT2(M)  = 2.0D0*T(M)
          FM(M,1) = 0.5D0*PI12/DSQRT(T(M))
        ENDDO
C
C       NOW COMPLETE TABLE BY FORWARD RECURRENCE
        DO MIND=1,LAM
          MVAL  = MIND-1
          COEFF = DFLOAT(MVAL+MVAL+1)
          DO M=1,N
            FM(M,MIND+1) = COEFF*FM(M,MIND)/TT2(M)
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
      SUBROUTINE BOYSGEN(LAM)
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
C  PARAMETER DETERMINED BY LAM.                                        C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION X1(MB2),X2(MB2),X3(MB2),X4(MB2)
      DIMENSION F1(MB2,ML4),F2(MB2,ML4),F3(MB2,ML4),
     &          F4(MB2,ML4),FS(MB2,ML4)
C
C     EVALUATION PARAMETERS
      NTOT = 400
      XMIN = 0.00D0
      XMAX = 4.00D1
      HSTP = (XMAX-XMIN)/NTOT
C
C     VALUES WHICH DETERMINE EVALUATION METHOD
      EPSZR = 1.0D-11
      EPSPL = 1.7D+01
      EPSAS = 3.0D+01
C
C**********************************************************************C
C     THE FIRST STEP OF THIS SUBROUTINE IS TO EVALUATE THE REQUIRED    C
C     BOYS INTEGRALS. THIS IS DIVIDED INTO CASES, DEPENDING ON THE     C
C     MAGNITUDE OF THE ARGUMENT.                                       C
C -------------------------------------------------------------------- C
C        FS_M (X) = INT_{0}^{1} T^{2M} EXP(-X*T^{2}) DT                C
C -------------------------------------------------------------------- C
C   EVALUATED FOR ALL VALUES OF M IN THE RANGE 0 < M < LAM.            C
C             FOR ALL VALUES OF X IN THE RANGE X > 0.                  C
C**********************************************************************C
C
      OPEN(UNIT=8,FILE='plots/boysfunction.dat',STATUS='UNKNOWN')
      REWIND(UNIT=8)
      DO NX=0,NTOT
C
      X = XMIN + HSTP*NX
C
C     CASE 0: ARGUMENT OF THE BOYS FUNCTION IS X = 0.
C             THE VALUE OF THIS FUNCTION IS (2N+1).
      IF(X.LE.EPSZR) THEN
C
        N1     = 1
        X1(N1) = X
        CALL FUNFM(F1,X1,N1,LAM,1)
        DO JJ=1,LAM+1
          FS(N1,JJ) = F1(N1,JJ)
        ENDDO
C
C     CASE 1: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=17.
C             EVALUATE WITH LOCAL POLYNOMIAL EXPANSION OF ORDER 5,
C             AND RECURRENCE IN DIRECTION OF DECREASING M.
      ELSEIF(X.GT.EPSZR.AND.X.LE.EPSPL) THEN
C
        N2     = 1
        X2(N2) = X
        CALL FUNFM(F2,X2,N2,LAM,2)
        DO JJ=1,LAM+1
          FS(N2,JJ) = F2(N2,JJ)
        ENDDO
C
C     CASE 2: ARGUMENT OF THE BOYS FUNCTION IS SMALLER THAN X=30.
C             EVALUATE USING ASYMPTOTIC FORMULA WITH EXPONENTIAL.
      ELSEIF(X.GT.EPSPL.AND.X.LE.EPSAS) THEN
C
        N3     = 1
        X3(N3) = X
        CALL FUNFM(F3,X3,N3,LAM,3)
        DO JJ=1,LAM+1
          FS(N3,JJ) = F3(N3,JJ)
        ENDDO
C
C     CASE 3: ARGUMENT OF THE BOYS FUNCTION IS LARGER THAN X=30.
C             EVALUATE USING ASYMPTOTIC FORMULA WITHOUT EXPONENTIAL.
      ELSE
C
        N4     = 1
        X4(N4) = X
        CALL FUNFM(F4,X4,N4,LAM,4)
        DO JJ=1,LAM+1
          FS(N4,JJ) = F4(N4,JJ)
        ENDDO
C
      ENDIF
C
      WRITE(8, *) X,(FS(1,L),L=1,LAM+1)
C
      ENDDO
      CLOSE(UNIT=8)
C
      RETURN
      END
C
C
      SUBROUTINE HGTFS(HABC,XYZEVAL,EXL,XYZ,NBAS,LAM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C              HH    HH  GGGGGG TTTTTTTT FFFFFFFF SSSSSS               C
C              HH    HH GG    GG   TT    FF      SS    SS              C
C              HH    HH GG         TT    FF      SS                    C
C              HHHHHHHH GG         TT    FFFFFF   SSSSSS               C
C              HH    HH GG   GGG   TT    FF            SS              C
C              HH    HH GG    GG   TT    FF      SS    SS              C
C              HH    HH  GGGGGG    TT    FF       SSSSSS               C
C                                                                      C
C -------------------------------------------------------------------- C
C  HGTFS GENERATES AN ARRAY OF HGTF FUNCTIONS EVALUATED AT A SET OF    C
C  COORDINATES XYZEVAL(3), USING THE GAUSSIAN PRODUCT THEOREM TO       C
C  TRANSFORM AN ARBITRARY PRODUCT OF TWO BASIS FUNCTIONS INTO A        C
C  FINITE SUM OF EQ-COEFFICIENTS AND HGTFS ON A SINGLE CENTRE, RP.     C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION HABC(MB2,MEQ)
      DIMENSION EXL(MBS,2),NBAS(2)
      DIMENSION XYZEVAL(3),XYZ(3,2)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
C     CALCULATE NUMBER OF UNIQUE CONTRIBUTIONS TO FINITE EXPANSION
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
C
C         SUM OF EXPONENTS
          PAB  = EXL(IBAS,1)+EXL(JBAS,2)
C
C         HGTF CENTRE COORDINATES
          PX = (EXL(IBAS,1)*XYZ(1,1)+EXL(JBAS,2)*XYZ(1,2))/PAB
          PY = (EXL(IBAS,1)*XYZ(2,1)+EXL(JBAS,2)*XYZ(2,2))/PAB
          PZ = (EXL(IBAS,1)*XYZ(3,1)+EXL(JBAS,2)*XYZ(3,2))/PAB
C
C         COORDINATES WRT LOCAL ORIGIN
          RPX = XYZEVAL(1)-PX
          RPY = XYZEVAL(2)-PY
          RPZ = XYZEVAL(3)-PZ
C
C         GAUSSIAN COMPONENT OF HGTF
          GSS = DEXP(-PAB*(RPX*RPX+RPY*RPY+RPZ*RPZ))
C
C         CALCULATE HERMITE POLYNOMIAL PRODUCTS FOR ALL {A,B,C}
          DO ITUV=1,NTUV
C
C           HERMITE POLYNOMIALS FOR THIS BASIS PAIR OVER ALL ITUV
            HALPH = HERMITE(PAB,RPX,IA(ITUV))
            HBETA = HERMITE(PAB,RPY,IB(ITUV))
            HGAMA = HERMITE(PAB,RPZ,IC(ITUV))
C
C           ADD TO BASIS FUNCTION PRODUCT
            HABC(M,ITUV) = HALPH*HBETA*HGAMA*GSS
C
          ENDDO
C
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
C
      IF(I.LT.0) THEN
        WRITE(6, *) 'In HERMITE: index less than zero. I = ',I
        WRITE(7, *) 'In HERMITE: index less than zero. I = ',I
        STOP
      ENDIF
C
      IF(I.GT.40) THEN
        WRITE(6, *) 'In HERMITE: index greater than 40. I = ',I
        WRITE(7, *) 'In HERMITE: index greater than 40. I = ',I
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
      SUBROUTINE HTFRMS(HABC,XYZEVAL,EXL,XYZ,NBAS,LAM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       HH    HH TTTTTTTT FFFFFFFF RRRRRRR  MM       MM  SSSSSS        C
C       HH    HH    TT    FF       RR    RR MMM     MMM SS    SS       C
C       HH    HH    TT    FF       RR    RR MMMM   MMMM SS             C
C       HHHHHHHH    TT    FFFFFF   RR    RR MM MM MM MM  SSSSSS        C
C       HH    HH    TT    FF       RRRRRRR  MM  MMM  MM       SS       C
C       HH    HH    TT    FF       RR    RR MM   M   MM SS    SS       C
C       HH    HH    TT    FF       RR    RR MM       MM  SSSSSS        C
C                                                                      C
C -------------------------------------------------------------------- C
C  HGTFS GENERATES AN ARRAY OF HGTF FOURIER TRANSFORMS EVALUATED AT A  C
C  SET OF COORDINATES XYZEVAL(3), USING THE GAUSSIAN PRODUCT THEOREM   C
C  TO TRANSFORM AN ARBITRARY PRODUCT OF TWO BASIS FUNCTIONS INTO A     C
C  FINITE SUM OF EQ-COEFFICIENTS AND HGTFS ON A SINGLE CENTRE, RP.     C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION EXL(MBS,2),NBAS(2)
      DIMENSION XYZEVAL(3),XYZ(3,2)
C
      COMPLEX*16 CONE,PHS,PHSQP
      COMPLEX*16 HABC(MB2,MEQ)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     UNIT COMPLEX NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     CALCULATE NUMBER OF UNIQUE CONTRIBUTIONS TO FINITE EXPANSION
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          M = M+1
C
C         SUM OF EXPONENTS
          PAB  = EXL(IBAS,1)+EXL(JBAS,2)
C
C         HGTF CENTRE COORDINATES
          PX = (EXL(IBAS,1)*XYZ(1,1)+EXL(JBAS,2)*XYZ(1,2))/PAB
          PY = (EXL(IBAS,1)*XYZ(2,1)+EXL(JBAS,2)*XYZ(2,2))/PAB
          PZ = (EXL(IBAS,1)*XYZ(3,1)+EXL(JBAS,2)*XYZ(3,2))/PAB
C
C         COORDINATES WRT LOCAL ORIGIN
          QPX = XYZEVAL(1)*PX
          QPY = XYZEVAL(2)*PY
          QPZ = XYZEVAL(3)*PZ
C
C         DOT PRODUCT FACTOR (COMPLEX WITH UNIT MAGNITUDE)
          QPVEC = QPX+QPY+QPZ
C         PHSQP = DEXP(-CONE*QPVEC)
          PHSQP = DCOS(QPVEC) + CONE*DSIN(QPVEC)
C
C         GAUSSIAN COMPONENT OF HTFRM
          QX2 = XYZEVAL(1)*XYZEVAL(1)
          QY2 = XYZEVAL(2)*XYZEVAL(2)
          QZ2 = XYZEVAL(3)*XYZEVAL(3)
          GSS = DEXP(-0.75D0*(QX2*QX2+QY2*QY2+QZ2*QZ2)/PAB)
C
C         INTEGRAL TRANSFORM TERM
          PRT = PAB*DSQRT(PAB)
          SCL = PI*PI12/PRT
C
C         CALCULATE HERMITE POLYNOMIAL PRODUCTS FOR ALL {A,B,C}
          DO ITUV=1,NTUV
C
C           TEST TO DETERMINE MULTIPLICATIVE FACTOR
            IF(MOD(ILAM(ITUV),4).EQ.0) THEN
              PHS = DCMPLX(1.0D0,0.0D0)
            ELSEIF(MOD(ILAM(ITUV),4).EQ.1) THEN
              PHS =-DCMPLX(0.0D0,1.0D0)
            ELSEIF(MOD(ILAM(ITUV),4).EQ.2) THEN
              PHS =-DCMPLX(1.0D0,0.0D0)
            ELSEIF(MOD(ILAM(ITUV),4).EQ.3) THEN
              PHS = DCMPLX(0.0D0,1.0D0)
            ENDIF
C
C           HERMITE POLYNOMIALS FOR THIS BASIS PAIR OVER ALL ITUV
            QXA = XYZEVAL(1)**(IA(ITUV))
            QXB = XYZEVAL(2)**(IB(ITUV))
            QXC = XYZEVAL(3)**(IC(ITUV))
C
C           ADD TO BASIS FUNCTION PRODUCT
            HABC(M,ITUV) = SCL*PHS*QXA*QXB*QXC*PHSQP*GSS
C
          ENDDO
C
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C  [12] EQ-COEFFS: BASIS FUNCTION OVERLAP SPIN-STRUCTURE FACTORS.      C
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
C   [E] EILSB3: GENERATE A VECTOR BATCH OF EILS COEFFICIENTS FOR BREIT.C
C   [F] EISLB3: GENERATE A VECTOR BATCH OF EISL COEFFICIENTS FOR BREIT.C
C   [G] EQLL: A RAW BLOCK OF EQLL COEFFICIENTS FOR EQLLMK.             C
C   [H] EQLS: A RAW BLOCK OF EQLS COEFFICIENTS FOR EQLSMK.             C
C   [I] EQSL: A RAW BLOCK OF EQSL COEFFICIENTS FOR EQSLMK.             C
C   [J] EQSS: A RAW BLOCK OF EQSS COEFFICIENTS FOR EQSSMK.             C
C   [K] ESGTF: SET OF ES-COEFFS OVER SPHERICAL HARMONICS AND HGTFS.    C
C   [L] EVRS: EXPANSION COEFFS IN HGTF OVERLAPS, CALLED IN ESGTF.      C
C   [M] ESTEPLM: SIMULTANEOUS INCREASE IN (L,M) FOR USE IN EVRS.       C
C   [N] ESTEPL: INCREMENT IN L FOR USE IN EVRS.                        C
C   [O] ESTEPN: INCREMENT IN N FOR USE IN EVRS.                        C
C -------------------------------------------------------------------- C
C   [A] RNORM1: A BLOCK OF TT' NORMALISATION COEFFS.                   C
C   [B] DNORM: NORM FOR A REAL OR COMPLEX PART OF EQ-COEFF LIST.       C
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
      COMPLEX*16 G11(MB2,MEQ),G21(MB2,MEQ)
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
        LQN(1) = LVAL(KQN(1))
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
        LQN(2) = LVAL(KQN(2))
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
      COMPLEX*16 G11(MB2,MEQ),G21(MB2,MEQ)
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
        LQN(1) = LVAL(KQN(1))
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
        LQN(2) = LVAL(KQN(2))
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
        LQN(1) = LVAL(KQN(1))
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
        LQN(2) = LVAL(KQN(2))
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
        LQN(1) = LVAL(KQN(1))
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
        LQN(2) = LVAL(KQN(2))
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
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
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
C     CALCULATE LQN VALUES
      LA = LVAL(KQ2(1))
      LB = LVAL(KQ2(2))
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
      CALL EQLL(E11,EX2,XY2,KQ2,MQ2,NB2,IQ)
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
      CALL EQLL(E21,EX2,XY2,KQ2,MQ2,NB2,IQ)
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
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
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
C     CALCULATE LQN VALUES
      LA = LVAL(KQ2(1))
      LB = LVAL(KQ2(2))
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
      CALL EQLS(E11,EX2,XY2,KQ2,MQ2,NB2,IQ)
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
      CALL EQLS(E21,EX2,XY2,KQ2,MQ2,NB2,IQ)
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
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
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
C     CALCULATE LQN VALUES
      LA = LVAL(KQ2(1))
      LB = LVAL(KQ2(2))
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
      CALL EQSL(E11,EX2,XY2,KQ2,MQ2,NB2,IQ)
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
      CALL EQSL(E21,EX2,XY2,KQ2,MQ2,NB2,IQ)
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
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
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
C     CALCULATE LQN VALUES
      LA = LVAL(KQ2(1))
      LB = LVAL(KQ2(2))
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
      CALL EQSS(E11,EX2,XY2,KQ2,MQ2,NB2,IQ)
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
      CALL EQSS(E21,EX2,XY2,KQ2,MQ2,NB2,IQ)
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
C  THIS ACTUALLY MAKES A VECTOR LIST OF THE EQLS NEEDED FOR BREIT.     C
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
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMPLEX*16 E11(MB2,MEQ,3),E21(MB2,MEQ,3)
      COMPLEX*16 E11X(MB2,MEQ),E11Y(MB2,MEQ),E11Z(MB2,MEQ),
     &           E21X(MB2,MEQ),E21Y(MB2,MEQ),E21Z(MB2,MEQ)
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
      CALL EQLS(E11X,EX2,XY2,KQ2,MQ2,NB2,1)
      CALL EQLS(E11Y,EX2,XY2,KQ2,MQ2,NB2,2)
      CALL EQLS(E11Z,EX2,XY2,KQ2,MQ2,NB2,3)
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
      CALL EQLS(E21X,EX2,XY2,KQ2,MQ2,NB2,1)
      CALL EQLS(E21Y,EX2,XY2,KQ2,MQ2,NB2,2)
      CALL EQLS(E21Z,EX2,XY2,KQ2,MQ2,NB2,3)
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
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMPLEX*16 E11(MB2,MEQ,3),E21(MB2,MEQ,3)
      COMPLEX*16 E11X(MB2,MEQ),E11Y(MB2,MEQ),E11Z(MB2,MEQ),
     &           E21X(MB2,MEQ),E21Y(MB2,MEQ),E21Z(MB2,MEQ)
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
      CALL EQSL(E11X,EX2,XY2,KQ2,MQ2,NB2,1)
      CALL EQSL(E11Y,EX2,XY2,KQ2,MQ2,NB2,2)
      CALL EQSL(E11Z,EX2,XY2,KQ2,MQ2,NB2,3)
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
      CALL EQSL(E21X,EX2,XY2,KQ2,MQ2,NB2,1)
      CALL EQSL(E21Y,EX2,XY2,KQ2,MQ2,NB2,2)
      CALL EQSL(E21Z,EX2,XY2,KQ2,MQ2,NB2,3)
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
      SUBROUTINE EQLL(ELL,EXL,XYZ,KQN,MQN,NBS,IQ)
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
C  EQLL GENERATES A BLOCK OF RAW EQLL/GQLL-COEFFICIENTS FOR A GIVEN    C
C  COMBINATION OF MQNS (WITH SIGN), TO BE USED BY EQLLMK.              C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ XYZ  - PAIR OF BASIS FUNCTION CARTESIAN CENTRES.                  C
C  ▶ KQN  - PAIR OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.              C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ NBS  - PAIR OF BASIS FUNCTION BLOCK LENGTHS.                      C
C  ▶ IQ   - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.            C
C  OUTPUT:                                                             C
C  ▶ ELL  - RAW EQLL-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION LQLAB(2),MQLAB(2)
      DIMENSION EXL(MBS,2),XYZ(3,2),KQN(2),LQN(2),JQN(2),MQN(2),NBS(2)
      DIMENSION RNLL(MBS,MBS)
C
      COMPLEX*16 CONE
      COMPLEX*16 ELL(MB2,MEQ),ESG(MB2,MEQ)
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
        ENDDO
      ENDDO
C
C     CALCULATE LQN VALUES
      LQN(1) = LVAL(KQN(1))
      LQN(2) = LVAL(KQN(2))
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
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ELL0 AND ELLZ
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLL
          DO M=1,MAXM
            DO ITUV=1,NTUV0
              ELL(M,ITUV) = ELL(M,ITUV) +     CAB*ESG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLL
          DO M=1,MAXM
            DO ITUV=1,NTUV0
              ELL(M,ITUV) = ELL(M,ITUV) + SIG*CAB*ESG(M,ITUV)
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLL
          DO M=1,MAXM
            DO ITUV=1,NTUV0
              ELL(M,ITUV) = ELL(M,ITUV) + SIG*CAB*ESG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLL
          DO M=1,MAXM
            DO ITUV=1,NTUV0
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
C     INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
      LAM0  = LQLAB(1)+LQLAB(2)
      NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C     GENERATE RNLL NORMALISATION CONSTANTS
      CALL RNORM1(RNLL,EXL,LQN,NBS,1)
C
C     NORMALISE THE ELLQ COEFFICIENT BLOCK
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
          DO ITUV=1,NTUV0
            ELL(M,ITUV) = RNLL(IBAS,JBAS)*ELL(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     σ_Y SPECIAL CASE: MULTIPLY ELLY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV0
            ELL(M,ITUV) = CONE*ELL(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE EQLS(ELS,EXL,XYZ,KQN,MQN,NBS,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 EEEEEEEE  QQQQQQ    LL       SSSSSS                  C
C                 EE       QQ    QQ   LL      SS    SS                 C
C                 EE      QQ      QQ  LL      SS                       C
C                 EEEEEE  QQ      QQ  LL       SSSSSS                  C
C                 EE      QQ      QQ  LL            SS                 C
C                 EE       QQ    QQ   LL      SS    SS                 C
C                 EEEEEEEE  QQQQQQ QQ LLLLLLLL SSSSSS                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQLS GENERATES A BLOCK OF RAW EQLS-COEFFICIENTS FOR A GIVEN         C
C  COMBINATION OF MQNS (WITH SIGN), TO BE USED BY EQLSMK.              C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ XYZ  - PAIR OF BASIS FUNCTION CARTESIAN CENTRES.                  C
C  ▶ KQN  - PAIR OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.              C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ NBS  - PAIR OF BASIS FUNCTION BLOCK LENGTHS.                      C
C  ▶ IQ   - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.            C
C  OUTPUT:                                                             C
C  ▶ ELS  - RAW EQLS-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION LQLAB(2),MQLAB(2)
      DIMENSION EXL(MBS,2),XYZ(3,2),KQN(2),LQN(2),JQN(2),MQN(2),NBS(2)
      DIMENSION RNLS(MBS,MBS)
      DIMENSION T2(MB2),T0(MB2)
C
      COMPLEX*16 CONE
      COMPLEX*16 ELS(MB2,MEQ),ESG(MB2,MEQ),ENSG(MB2,MEQ)
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
        ENDDO
      ENDDO
C
C     CALCULATE LQN VALUES
      LQN(1) = LVAL(KQN(1))
      LQN(2) = LVAL(KQN(2))
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
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ELS0 AND ELSZ
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
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
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)
      LQLAB(2) = LQN(2)-1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ELS0 AND ELSZ
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
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
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
              ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
              ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
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
C     MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
      LAM2  = LQN(1)+LQN(2)+2
      NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C     GENERATE RNLS NORMALISATION CONSTANTS
      CALL RNORM1(RNLS,EXL,LQN,NBS,2)
C
C     NORMALISE THE ELSQ COEFFICIENT BLOCK AND MULTIPLY BY FACTOR i
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
          DO ITUV=1,NTUV2
            ELS(M,ITUV) = CONE*RNLS(IBAS,JBAS)*ELS(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     σ_Y SPECIAL CASE: MULTIPLY ELSY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV2
            ELS(M,ITUV) = CONE*ELS(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE EQSL(ESL,EXL,XYZ,KQN,MQN,NBS,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 EEEEEEEE  QQQQQQ    SSSSSS  LL                       C
C                 EE       QQ    QQ  SS    SS LL                       C
C                 EE      QQ      QQ SS       LL                       C
C                 EEEEEE  QQ      QQ  SSSSSS  LL                       C
C                 EE      QQ      QQ       SS LL                       C
C                 EE       QQ    QQ  SS    SS LL                       C
C                 EEEEEEEE  QQQQQQ QQ SSSSSS  LLLLLLLL                 C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQSL GENERATES A BLOCK OF RAW EQSL-COEFFICIENTS FOR A GIVEN         C
C  COMBINATION OF MQNS (WITH SIGN), TO BE USED BY EQSLMK.              C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ XYZ  - PAIR OF BASIS FUNCTION CARTESIAN CENTRES.                  C
C  ▶ KQN  - PAIR OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.              C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ NBS  - PAIR OF BASIS FUNCTION BLOCK LENGTHS.                      C
C  ▶ IQ   - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.            C
C  OUTPUT:                                                             C
C  ▶ ESL  - RAW EQSL-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION LQLAB(2),MQLAB(2)
      DIMENSION EXL(MBS,2),XYZ(3,2),KQN(2),LQN(2),JQN(2),MQN(2),NBS(2)
      DIMENSION RNSL(MBS,MBS)
      DIMENSION T2(MB2),T0(MB2)
C
      COMPLEX*16 CONE
      COMPLEX*16 ESL(MB2,MEQ),ESG(MB2,MEQ),ENSG(MB2,MEQ)
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
        ENDDO
      ENDDO
C
C     CALCULATE LQN VALUES
      LQN(1) = LVAL(KQN(1))
      LQN(2) = LVAL(KQN(2))
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
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESL0 AND ESLZ
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ESLX AND ESLY
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
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
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESL0 AND ESLZ
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
              ESL(M,ITUV) = ESL(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
              ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ESLX AND ESLY
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
              ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
              ESL(M,ITUV) = ESL(M,ITUV) +     TK*ENSG(M,ITUV)
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
C     MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
      LAM2  = LQN(1)+LQN(2)+2
      NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C     GENERATE RNSL NORMALISATION CONSTANTS
      CALL RNORM1(RNSL,EXL,LQN,NBS,3)
C
C     NORMALISE THE ESLQ COEFFICIENT BLOCK AND MULTIPLY BY FACTOR -i
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
          DO ITUV=1,NTUV2
            ESL(M,ITUV) =-CONE*RNSL(IBAS,JBAS)*ESL(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     σ_Y SPECIAL CASE: MULTIPLY ESLY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV2
            ESL(M,ITUV) = CONE*ESL(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE EQSS(ESS,EXL,XYZ,KQN,MQN,NBS,IQ)
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
C  EQSS GENERATES A BLOCK OF RAW EQSS-COEFFICIENTS FOR A GIVEN         C
C  COMBINATION OF MQNS (WITH SIGN), TO BE USED BY EQSSMK.              C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ XYZ  - PAIR OF BASIS FUNCTION CARTESIAN CENTRES.                  C
C  ▶ KQN  - PAIR OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.              C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ NBS  - PAIR OF BASIS FUNCTION BLOCK LENGTHS.                      C
C  ▶ IQ   - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.            C
C  OUTPUT:                                                             C
C  ▶ ESS  - RAW EQSS-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION LQLAB(2),MQLAB(2)
      DIMENSION EXL(MBS,2),XYZ(3,2),KQN(2),LQN(2),JQN(2),MQN(2),NBS(2)
      DIMENSION RNSS(MBS,MBS)
      DIMENSION T22(MB2),T20(MB2),T02(MB2),T00(MB2)
C
      COMPLEX*16 CONE
      COMPLEX*16 ESS(MB2,MEQ),ESG(MB2,MEQ),ENSG(MB2,MEQ)
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
        ENDDO
      ENDDO
C
C     CALCULATE LQN VALUES
      LQN(1) = LVAL(KQN(1))
      LQN(2) = LVAL(KQN(2))
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
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESS0 AND ESSZ
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV0
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
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)+1
      LQLAB(2) = LQN(2)-1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESS0 AND ESSZ
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) +    TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
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
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)-1
      LQLAB(2) = LQN(2)+1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESS0 AND ESSZ
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
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
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)-1
      LQLAB(2) = LQN(2)-1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESS0 AND ESSZ
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          LAM4  = LQLAB(1)+LQLAB(2)+4
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
          NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T00(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ESG)
          CALL ESTEPN(ENSG,ESG,LAM2,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV4
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          LAM4  = LQLAB(1)+LQLAB(2)+4
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
          NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T00(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ESG)
          CALL ESTEPN(ENSG,ESG,LAM2,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV4
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          LAM4  = LQLAB(1)+LQLAB(2)+4
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
          NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T00(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ESG)
          CALL ESTEPN(ENSG,ESG,LAM2,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV4
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
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
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          LAM4  = LQLAB(1)+LQLAB(2)+4
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
          NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T00(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ESG)
          CALL ESTEPN(ENSG,ESG,LAM2,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV4
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
      CALL RNORM1(RNSS,EXL,LQN,NBS,4)
C
C     MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
      LAM4  = LQN(1)+LQN(2)+4
      NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C     NORMALISE THE ESSQ COEFFICIENT BLOCK
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
          DO ITUV=1,NTUV4
            ESS(M,ITUV) = RNSS(IBAS,JBAS)*ESS(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     σ_Y SPECIAL CASE: MULTIPLY ESSY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV4
            ESS(M,ITUV) = CONE*ESS(M,ITUV)
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
C     INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
      LAM  = LQN(1)+LQN(2)
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
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
      DIMENSION NBAS(2),LQN(2),MQN(2)
C
      COMPLEX*16 ESG(MB2,MEQ),ETEMP(MB2*MRC)
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
      DIMENSION PP(MB2),PX(MB2),PY(MB2),PZ(MB2)
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
      SUBROUTINE RNORM1(RNTT,EXL,LQN,NBAS,ITT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         RRRRRRR  NN    NN  OOOOOO  RRRRRRR  MM       MM  11          C
C         RR    RR NNN   NN OO    OO RR    RR MMM     MMM 111          C
C         RR    RR NNNN  NN OO    OO RR    RR MMMM   MMMM  11          C
C         RR    RR NN NN NN OO    OO RR    RR MM MM MM MM  11          C
C         RRRRRRR  NN  NNNN OO    OO RRRRRRR  MM  MMM  MM  11          C
C         RR    RR NN   NNN OO    OO RR    RR MM   M   MM  11          C
C         RR    RR NN    NN  OOOOOO  RR    RR MM       MM 1111         C
C                                                                      C
C -------------------------------------------------------------------- C
C  RNORM1 GENERATES TT' RADIAL BASIS FUNCTION NORMALISATION CONSTANTS, C
C  BUT KEEPS THEM ISOLATED IN A PAIR OF COLUMNS.                       C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ LQN  - PAIR OF ORBITAL QUANTUM NUMBERS.                           C
C  ▶ NBAS - PAIR OF PARAMETER LIST LENGTHS.                            C
C  ▶ ITT  - COMPONENT-TYPE OVERLAPS.                                   C
C  OUTPUT:                                                             C
C  ▶ RNTT - BATCH OF NORMALISATION CONSTANTS.                          C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RNTT(MBS,MBS)
      DIMENSION RN1(MBS),RN2(MBS)
      DIMENSION EXL(MBS,2),LQN(2),NBAS(2)
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     LEFT-HAND NORMALISATION LIST
      GL = TWLG-GAMLOG(2*LQN(1)+3)
      GS = TWLG-GAMLOG(2*LQN(1)+5)
      UL = DFLOAT(LQN(1))+1.5D0
      US = DFLOAT(LQN(1))+0.5D0
      DO IBAS=1,NBAS(1)
        ELOG = DLOG(2.0D0*EXL(IBAS,1))
        IF(ITT.EQ.1.OR.ITT.EQ.2) THEN
          RN1(IBAS) = DEXP(0.5D0*(GL+UL*ELOG))
        ELSEIF(ITT.EQ.3.OR.ITT.EQ.4) THEN
          RN1(IBAS) = DEXP(0.5D0*(GS+US*ELOG))
        ENDIF
      ENDDO
C
C     RIGHT-HAND NORMALISATION LIST
      GL = TWLG-GAMLOG(2*LQN(2)+3)
      GS = TWLG-GAMLOG(2*LQN(2)+5)
      UL = DFLOAT(LQN(2))+1.5D0
      US = DFLOAT(LQN(2))+0.5D0
      DO JBAS=1,NBAS(2)
        ELOG = DLOG(2.0D0*EXL(JBAS,2))
        IF(ITT.EQ.1.OR.ITT.EQ.3) THEN
          RN2(JBAS) = DEXP(0.5D0*(GL+UL*ELOG))
        ELSEIF(ITT.EQ.2.OR.ITT.EQ.4) THEN
          RN2(JBAS) = DEXP(0.5D0*(GS+US*ELOG))
        ENDIF
      ENDDO
C
C     FULL MATRIX SET OF NORMALISATION CONSTANTS
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          RNTT(IBAS,JBAS) = RN1(IBAS)*RN2(JBAS)
        ENDDO
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
C     SENSITIVITY TOLERANCE PARAMETER
      EPS = 1.0D-10
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
        IF(DABS(ECMP(N)).GT.EPS) THEN
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
C
      FUNCTION EXTINT11(Y,I,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     EEEEEEEE XX     XX TTTTTTTT IIII NN    NN TTTTTTTT 11   11       C
C     EE        XX   XX     TT     II  NNN   NN    TT   111  111       C
C     EE         XX XX      TT     II  NNNN  NN    TT    11   11       C
C     EEEEEE      XXX       TT     II  NN NN NN    TT    11   11       C
C     EE         XX XX      TT     II  NN  NNNN    TT    11   11       C
C     EE        XX   XX     TT     II  NN   NNN    TT    11   11       C
C     EEEEEEEE XX     XX    TT    IIII NN    NN    TT   1111 1111      C
C                                                                      C
C -------------------------------------------------------------------- C
C  EXTINT11 DETERMINES THE VALUE OF THE CONTRIBUTION TO AN INTEGRAL    C
C  BASED ON THE REPEATED 11-POINT NEWTON-COTES FORMULA, WHERE Y IS     C
C  THE VALUE F(X_I), I IS THE INDEX OF INTEREST, AND THE START AND     C
C  END POINTS OF THE INTEGRAL ARE F(X_0) AND F(X_N) RESPECTIVELY.      C
C  IT IS NUMERICALLY ADVANTAGEOUS TO MULTIPLY BY 5.0D+0/2.99376D+5     C
C  AFTER THE SUM IS CALCULATED.                                        C
C**********************************************************************C
C
C     FUNCTION ONLY ALLOWS I BETWEEN 0 AND N
      IF(I.LT.0.OR.I.GT.N) THEN
        WRITE(6, *) 'In EXTINT11: invalid index I. I = ',I
        WRITE(7, *) 'In EXTINT11: invalid index I. I = ',I
        STOP
      ENDIF
C
C     FOR THIS FORMULA TO APPLY, N MUST BE A MULTIPLE OF 10
      IF(MOD(N,10).NE.0.OR.N.LT.10) THEN
        M = MOD(N,10)
        WRITE(6, *) 'In EXTINT11: invalid discretisation N. N = ',M
        WRITE(7, *) 'In EXTINT11: invalid discretisation N. N = ',M
        STOP
      ENDIF
C
      IF(I.EQ.0.OR.I.EQ.N) THEN
        EXTINT11 = 16067.0D0*Y
      ELSEIF(MOD(I,10).EQ.0) THEN
        EXTINT11 = 32134.0D0*Y
      ELSEIF(MOD(I,10).EQ.1.OR.MOD(I,10).EQ.9) THEN
        EXTINT11 = 106300.0D0*Y
      ELSEIF(MOD(I,10).EQ.2.OR.MOD(I,10).EQ.8) THEN
        EXTINT11 =-48525.0D0*Y
      ELSEIF(MOD(I,10).EQ.3.OR.MOD(I,10).EQ.7) THEN
        EXTINT11 = 272400.0D0*Y
      ELSEIF(MOD(I,10).EQ.4.OR.MOD(I,10).EQ.6) THEN
        EXTINT11 =-260550.0D0*Y
      ELSEIF(MOD(I,10).EQ.5) THEN
        EXTINT11 = 427368.0D0*Y
      ENDIF
C
      RETURN
      END

