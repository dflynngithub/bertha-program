C     THIS IS JUST A BERTHA SHELL THAT BUDGETS INTEGRAL NUMBERS AND COSTS
C
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
C     SSSSSS  WW         WW IIII RRRRRRR  LL      EEEEEEEE SSSSSS      C
C    SS    SS WW         WW  II  RR    RR LL      EE      SS    SS     C
C    SS       WW         WW  II  RR    RR LL      EE      SS           C
C     SSSSSS  WW    W    WW  II  RR    RR LL      EEEEEE   SSSSSS      C
C          SS WW   WWW   WW  II  RRRRRRR  LL      EE            SS     C
C    SS    SS  WW WW WW WW   II  RR    RR LL      EE      SS    SS     C
C     SSSSSS    WW     WW   IIII RR    RR LLLLLLL EEEEEEEE SSSSSS      C
C                                                                      C
C -------------------------------------------------------------------- C
C        A RELATIVISTIC MOLECULAR ELECTRONIC STRUCTURE PROGRAM         C
C            BASED ON THE ANALYTIC FINITE BASIS SET METHOD.            C
C                                                                      C
C                 H.M.QUINEY, H.SKAANE (OXFORD, 1997)                  C
C                       D. FLYNN (UNIMELB, 2019)                       C
C -------------------------------------------------------------------- C
C                          HAMILTONIANS (HMLT)                         C
C                          -------------------                         C
C ▶ 'NORL' NON-RELATIVISTIC HAMILTONIAN (PAULI EQUATION).              C
C ▶ 'BARE' BARE NUCLEUS DIRAC HAMILTONIAN (NO ELECTRON INTERACTION).   C
C ▶ 'DHFR' DIRAC-COULOMB HAMILTONIAN.                                  C
C ▶ 'DHFP' DIRAC-COULOMB HAMILTONIAN (+1ST ORDER BREIT).               C
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
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
C     READ DATA FROM USER-SPECIFIED INPUT FILE
      CALL CARDIN
C
C     PRINT SUMMARY OF INPUT DATA
      CALL INPUT
C
C     PRINT MEMORY ALLOCATION SUMMARY
      CALL MEMORY
C
C     REDUCE ALL THE INFO TO A BUDGET
      CALL REDUCE
C
C     BLANK EQ-COEFFICIENT GENERATOR
      IF(EQFILE) THEN
        CALL EQSAVE
      ENDIF
C
C     SUCCESSFUL EXIT
      END PROGRAM
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
      CHARACTER*7  HMINT(10),PTYPE(10)
      CHARACTER*9  BTYP
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 COEF(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,MFT),XNUC(MCT,MFT),NNUC(MCT),NMDL(MCT)
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
C
C     CALCULATION TREE CHECK -- IF UNKNOWN, EXIT.
      IF(TREE.NE.'HFSCF'.AND.TREE.NE.'MBPTN'.AND.TREE.NE.'MCSCF'.AND.
     &   TREE.NE.'EXPVL'.AND.TREE.NE.'PLOTS') THEN
        WRITE(6, *) 'In CARDIN: invalid calculation tree. ',TREE
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
        ENDIF
C
C       NUMBER OF KAPPA VALUES FOR THIS ATOM
        NKAP(IZ) = 2*LMAX+1
C
C       CHECK THAT LMAX CAN BE SUPPORTED BY SYSTEM PARAMETERS
        IF(2*LMAX+1.GT.MKP) THEN
          WRITE(6, *) 'In CARDIN: LMAX runs outside MKP storage.'
        ENDIF
C
C       ELECTRON CONFIGURATION CHOICE AND NUCLEAR ELECTRONIC CHARGE
        READ(5, *) DUMLIN
        READ(5, *) CNFG(IZ),IQNC(IZ)
C
C       CHECK THAT CNFG IS ALLOWED
        IF(CNFG(IZ).NE.'AUFBAU'.AND.CNFG(IZ).NE.'MANUAL') THEN
          WRITE(6, *) 'In CARDIN: illegal CNFG option.',CNFG(IZ)
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
C
        WRITE(6, *) 'In CARDIN: invalid value NOPN.',NOPN
C
C     CLOSED-SHELL MOLECULAR CONFIGURATION
      ELSEIF(NOPN.EQ.0) THEN
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
      IF(READIN) THEN
        INQUIRE(FILE=WFNFL,EXIST=FILETHERE)
        IF(.NOT.FILETHERE) THEN
          WRITE(6, *) 'Expected existing wfn file but none exist!',WFNFL
        ENDIF
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
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/TMMD/TMMDS(9)
      COMMON/TSCF/TSCFS(37)
C
C     INITIALISE TIME COUNTERS
      DO NT=1,35
        TSCFS(NT) = 0.0D0
      ENDDO
C
      DO NT=1,9
        TMMDS(NT) = 0.0D0
      ENDDO
C
C     PRINT A BIG SECTION HEADER
      WRITE(6, *) ' '
      WRITE(6, *) REPEAT('#',72)
      WRITE(6, *) REPEAT(' ',30),'INPUT SUMMARY'
      WRITE(6, *) REPEAT('#',72)
      WRITE(6, *) ' '
C
C     TITLE FOR INPUT OPTIONS
      WRITE(6, *) ' '
      WRITE(6, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',21),'Calculation tree and data files'
      WRITE(6, *) REPEAT('=',72)
C
C     CALCULATION TREE
      IF(TREE.EQ.'HFSCF') THEN
        WRITE(6, *) 'Calculation tree:',REPEAT(' ',43),'Hartree-Fock'
      ELSEIF(TREE.EQ.'MBPTN') THEN
        WRITE(6, *) 'Calculation tree:',REPEAT(' ',50),'MBPTN'
      ELSEIF(TREE.EQ.'MCSCF') THEN
        WRITE(6, *) 'Calculation tree:',REPEAT(' ',50),'MCSCF'
      ELSEIF(TREE.EQ.'EXPVL') THEN
        WRITE(6, *) 'Calculation tree:',REPEAT(' ',44),'Expct. vals'
      ELSEIF(TREE.EQ.'PLOTS') THEN
        WRITE(6, *) 'Calculation tree:',REPEAT(' ',47),'Plotting'
      ENDIF
C
C     PRINT THE HAMILTONIAN OPTION
      WRITE(6, *) 'Hamiltonian type:',REPEAT(' ',51),HMLT
C
C     CONFIRM SOLUTION SPACE DIMENSION OR EXIT
      IF(NDIM.LE.MDM) THEN
        WRITE(6, *) 'Fock matrix dimension: ',REPEAT(' ',37),NDIM
      ELSEIF(NDIM.GT.MDM) THEN
        WRITE(6, *) 'In INPUT: matrix dimension too big. NDIM = ',NDIM
        STOP
      ENDIF
C
C     NEW START OR READ IN HFSCF EXPANSION COEFFICIENTS
      IF(READIN) THEN
        WRITE(6, *) 'SCF density matrix:',REPEAT(' ',46),'Read in'
      ELSE
        WRITE(6, *) 'SCF density matrix:',REPEAT(' ',44),'New start'
      ENDIF
C
C     EQ-COEFFICIENT CALCULATION
      IF(EQFILE) THEN
        WRITE(6, *) 'Eq-coefficients:',REPEAT(' ',44),'Save to file'
      ELSE
        WRITE(6, *) 'Eq-coefficients:',REPEAT(' ',48),'By batch'
      ENDIF
C
C     OPENMP PARALLEL OPTION
      IF(OPENMP) THEN
        WRITE(6, *) 'OpenMP parallel option:',REPEAT(' ',42),'Enabled'
        WRITE(6, *) 'OpenMP thread limit:'   ,REPEAT(' ',40),NPRCSR
      ELSE
        WRITE(6, *) 'OpenMP parallel option:',REPEAT(' ',41),'Disabled'
      ENDIF
C
C     SECTION FOR FILE NAMES
      WRITE(6, *) REPEAT('-',72)
C
C     PRINT INPUT FILE NAME
      LF = LEN(TRIM(MOLCL))
      WRITE(6, *) 'Molecule name:',REPEAT(' ',58-LF),TRIM(MOLCL)
C
C     PRINT FILE OUTPUT NAMES
      LN = LEN(TRIM(OUTFL))
      WRITE(6, *) 'Output file name: ',REPEAT(' ',54-LN),TRIM(OUTFL)
C
C     SECTION FOR PHYSICAL PARAMETERS
      WRITE(6, *) REPEAT('-',72)
C
C     PRINT SPEED OF LIGHT
20    FORMAT(1X,A,43X,F14.10)
      WRITE(6,20) 'Speed of light:',CV
C
C     END OF INPUT SUMMARY
      WRITE(6, *) REPEAT('=',72)
      WRITE(6, *) ' '
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
C     TITLE FOR MEMORY SUMMARY
      WRITE(6, *) ' '
      WRITE(6, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',24),'Stack memory allocation'
      WRITE(6, *) REPEAT('=',72)
C
C     APPROXIMATE MEMORY STORAGE ALLOCATION
      NCMEM = 2*MDM*MDM
      NDMEM = 2*3*MDM*MDM
      NEMEM = 4*4*MFL
      NHMEM = 20*MFL
      NMMEM = 2*18*MDM*MDM
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
      WRITE(6, *) REPEAT('-',72)
      WRITE(6,21) 'COEF','Expansion coeffs.          ',NCMEM,SCMEM
      WRITE(6,21) 'DENS','Density matrices           ',NDMEM,SDMEM
      WRITE(6,21) 'MTRX','Molecular Fock matrices    ',NMMEM,SMMEM
      WRITE(6,21) 'EQTT','Eq-coefficient file        ',NEMEM,SEMEM
      WRITE(6,21) 'RCTT','R-integral class file      ',NHMEM,SHMEM
      WRITE(6, *) REPEAT('-',72)
      WRITE(6,21) '    ','Total                      ',NTMEM,STMEM
C
C     END OF MEMORY SUMMARY
      WRITE(6, *) REPEAT('=',72)
      WRITE(6, *) ' '
C
      RETURN
      END
C
C
      SUBROUTINE REDUCE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        RRRRRRR  EEEEEEEE DDDDDDD  UU    UU  CCCCCC  EEEEEEEE         C
C        RR    RR EE       DD    DD UU    UU CC    CC EE               C
C        RR    RR EE       DD    DD UU    UU CC       EE               C
C        RR    RR EEEEEE   DD    DD UU    UU CC       EEEEEE           C
C        RRRRRRR  EE       DD    DD UU    UU CC       EE               C
C        RR    RR EE       DD    DD UU    UU CC    CC EE               C
C        RR    RR EEEEEEEE DDDDDDD   UUUUUU   CCCCCC  EEEEEEEE         C
C                                                                      C
C -------------------------------------------------------------------- C
C     REDUCE MAKES A SERIES OF INTEGRAL COST BUDGETS.                  C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      LOGICAL SELECTIONRULES
C
      INTEGER*16 N2EB(5,7),N2EI(5,7),N2EM(5,7),N2EC(5,7),N2ES(5,7)
      DIMENSION  F2EB(5,7),F2EI(5,7),F2EM(5,7),F2EC(5,7),F2ES(5,7)
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBAS(4),LQN(4),JQN(4)
      DIMENSION INDEX(MCT,-(MEL+1):MEL,MKP),ICNT(4),ITN(2)
      DIMENSION MAPTTTT(4,4)
      DIMENSION NMLEV(5),TMLEV(5)
C
      COMPLEX*16 RR(MB2,16)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/IRCM/IEAB,IECD,IRIJ(MBS,MBS)
      COMMON/LSHF/SHLEV(4),SHLV,ILEV
      COMMON/QNMS/LABICN(MDM),LABKQN(MDM),LABMQN(MDM)
C
C     TWO-ELECTRON COMPONENT OVERLAP ADDRESSES
      DATA MAPTTTT/1,0,0,2,0,5,6,0,0,7,8,0,3,0,0,4/
C
      SELECTIONRULES = .TRUE.
C
C**********************************************************************C
C     SCF CALCULATION CHOICES AND ARRAY INITIALISATION                 C
C**********************************************************************C
C
C     INITIALISE INTEGRAL INCLUSION LEVEL VALUES
      DO N=1,5
        NMLEV(N) = 0
        TMLEV(N) = 0.0D0
      ENDDO
C
C     RESET TWO-ELECTRON INTEGRAL SCREENING COUNTERS
      DO MCNT=1,5
        DO ITT=1,7
          N2EB(MCNT,ITT) = 0
          N2EI(MCNT,ITT) = 0
          N2EM(MCNT,ITT) = 0
          N2EC(MCNT,ITT) = 0
          N2ES(MCNT,ITT) = 0
        ENDDO
      ENDDO
C
C**********************************************************************C
C     PRELIMINARY STORAGE OF ARRAYS AND TITLES                         C
C**********************************************************************C
C
C     COMPONENT OVERLAP LABELS TO LOOP OVER
      ITSTRT = 4
      ITSTOP = 1
      ITSKIP =-3
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
C     LOOP OVER ALL ATOMIC CENTRES (USE INDEX 1000)                    C
C**********************************************************************C
C
C     LOOP OVER CENTRE A
      DO 1000 ICNTA=1,NCNT
        ICNT(1) = ICNTA
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 1000 ICNTB=1,ICNTA
        ICNT(2) = ICNTB
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C     LOOP OVER CENTRE C
      DO 1000 ICNTC=1,NCNT
        ICNT(3) = ICNTC
C
C       CARTESIAN COORDINATES OF CENTRE C
        XYZ(1,3) = BXYZ(1,ICNTC)
        XYZ(2,3) = BXYZ(2,ICNTC)
        XYZ(3,3) = BXYZ(3,ICNTC)
C
C     LOOP OVER CENTRE D
      DO 1000 ICNTD=1,NCNT
        ICNT(4) = ICNTD
C
C       CARTESIAN COORDINATES OF CENTRE D
        XYZ(1,4) = BXYZ(1,ICNTD)
        XYZ(2,4) = BXYZ(2,ICNTD)
        XYZ(3,4) = BXYZ(3,ICNTD)
C
C     NUMBER OF NUCLEAR CENTRES INVOLVED IN THIS OVERLAP
      MCNT = NCNTRS(ICNTA,ICNTB,ICNTC,ICNTD)
C
C     SKIP ONE-CENTRE CONTRIBUTIONS (DEFER TO RACAH ALGEBRA ROUTINE)
      IF(SELECTIONRULES) THEN
        IF(MCNT.EQ.1.AND.RACAH1) THEN
          GOTO 1001
        ENDIF
      ENDIF
C
C**********************************************************************C
C     LOOP OVER ALL KQN SYMMETRY TYPES (USE INDEX 2000)                C
C**********************************************************************C
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
        JQN(1) = 2*IABS(KQN(1))-1
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),ICNTA)
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
        JQN(2) = 2*IABS(KQN(2))-1
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
C
C     LOOP OVER KQN(C) VALUES
      DO 2000 KC=1,NKAP(ICNTC)
C
C       QUANTUM NUMBERS FOR BLOCK C
        KQN(3) = KAPA(KC,ICNTC)
        IF(KQN(3).LT.0) THEN
          LQN(3) =-KQN(3)-1
        ELSE
          LQN(3) = KQN(3)
        ENDIF
        JQN(3) = 2*IABS(KQN(3))-1
C
C       BASIS EXPONENTS FOR BLOCK C
        NBAS(3) = NFNC(LQN(3),ICNTC)
C
C     LOOP OVER KQN(D) VALUES
      DO 2000 KD=1,NKAP(ICNTD)
C
C       QUANTUM NUMBERS FOR BLOCK D
        KQN(4) = KAPA(KD,ICNTD)
        IF(KQN(4).LT.0) THEN
          LQN(4) =-KQN(4)-1
        ELSE
          LQN(4) = KQN(4)
        ENDIF
        JQN(4) = 2*IABS(KQN(4))-1
C
C       BASIS EXPONENTS FOR BLOCK D
        NBAS(4) = NFNC(LQN(4),ICNTD)
C
C     THIS UNIQUELY DEFINES A FULL SET OF RC(AB|CD) INTEGRALS -- RESET
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
          IRIJ(IBAS,JBAS) = 1
        ENDDO
      ENDDO
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON KQN                           C
C**********************************************************************C
C
      IF(SELECTIONRULES) THEN
C     ATOM-CENTRED SELECTION RULES (ONLY APPLIES IF RACAH1 SWITCHED OFF)
      IF(MCNT.EQ.1) THEN
C
C       LQN PAIR PARITY (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQN(1)+LQN(2),2)
        IPARCD = MOD(LQN(3)+LQN(4),2)
C
C       LQN PAIR PARITY SELECTION RULE
        IF(IPARAB.NE.IPARCD) THEN
          GOTO 2001
        ENDIF
C
C       JQN TRIANGLE RULE CHECK FOR MULTIPOLE EXPANSION (ATOM-CENTRED)
        NUI = MAX0(IABS(JQN(1)-JQN(2))/2,IABS(JQN(3)-JQN(4))/2)
        NUF = MIN0(    (JQN(1)+JQN(2))/2,    (JQN(3)+JQN(4))/2)
        IF(NUI.GT.NUF) THEN
          GOTO 2001
        ENDIF
C
C       ADDITIONAL LQN SELECTION RULE PARITY ANALYSIS
        ISELK = 0
        DO NU=NUI,NUF
C
C         A AND B: LQN(1)+LQN(2)+NU EVEN OR ODD (0 IF EVEN, 1 IF ODD)
          IPARAB = MOD(LQN(1)+LQN(2)+NU,2)
          IPARCD = MOD(LQN(3)+LQN(4)+NU,2)
C
C         LQNA+LQNB+NU AND LQNC+LQND+NU ARE BOTH EVEN
          IF(IPARAB.EQ.0.AND.IPARCD.EQ.0) ISELK = 1
C
        ENDDO
        IF(ISELK.EQ.0) GOTO 2001
C
      ENDIF
      ENDIF
C
C**********************************************************************C
C     LOOP OVER ALL |MQN| PROJECTIONS (INDEX 3000)                     C
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
C     LOOP OVER |MQN(C)| VALUES
      DO 3000 MC=1,IABS(KQN(3))
        MQN(3) = 2*MC-1
C
C     LOOP OVER |MQN(D)| VALUES
      DO 3000 MD=1,IABS(KQN(4))
        MQN(4) = 2*MD-1
C
C**********************************************************************C
C     MOLECULAR SELECTION RULES BASED ON MQN                           C
C**********************************************************************C
C
C     SPIN PROJECTION CONSERVED ALONG Z-AXIS FOR LINEAR MOLECULES
      IF(SELECTIONRULES) THEN
      IF(ISYM.EQ.1.OR.ISYM.EQ.2) THEN
        IF(MQN(1).EQ.MQN(2).AND.MQN(3).EQ.MQN(4)) GOTO 3003
        IF(MQN(1).EQ.MQN(3).AND.MQN(2).EQ.MQN(4)) GOTO 3003
        IF(MQN(1).EQ.MQN(4).AND.MQN(2).EQ.MQN(3)) GOTO 3003
        GOTO 3001
      ENDIF
3003  CONTINUE
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
        IF(ISELM.EQ.0) GOTO 3001
      ENDIF
      ENDIF
C
C**********************************************************************C
C     IDENTIFICATION OF ERI SYMMETRIES AVAILABLE TO THIS BLOCK         C
C**********************************************************************C
C
C     CALCULATE BLOCK INDICES FOR {ABCD} COMBINATIONS
      IQ1 = INDEX(ICNTA,KQN(1),MQN(1))
      IQ2 = INDEX(ICNTB,KQN(2),MQN(2))
      IQ3 = INDEX(ICNTC,KQN(3),MQN(3))
      IQ4 = INDEX(ICNTD,KQN(4),MQN(4))
C
C     COMBINED BLOCK INDEX IN A TWO-FUNCTION LIST
      IQL = (IQ1*(IQ1-1))/2 + IQ2
      IQR = (IQ3*(IQ3-1))/2 + IQ4
C
C     SKIP CONTRIBUTIONS THAT ARISE BY PERMUTATION OF INTEGRALS
      IF(IQ1.LT.IQ2) GOTO 3001
      IF(IQ3.LT.IQ4) GOTO 3001
      IF(IQL.LT.IQR) GOTO 3001
C
      IF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQL.GT.IQR) THEN
C       IQ1 > IQ2, IQ3 > IQ4, IQL > IQR
        ITSCF = 1
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.GT.IQ4.AND.IQL.EQ.IQR) THEN
C       IQ1 > IQ2, IQ3 > IQ4, IQL = IQR
        ITSCF = 2
      ELSEIF(IQ1.GT.IQ2.AND.IQ3.EQ.IQ4.AND.IQL.GT.IQR) THEN
C       IQ1 > IQ2, IQ3 = IQ4, IQL > IQR
        ITSCF = 3
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.GT.IQ4.AND.IQL.GT.IQR) THEN
C       IQ1 = IQ2, IQ3 > IQ4, IQL > IQR
        ITSCF = 4
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQL.GT.IQR) THEN
C       IQ1 = IQ2, IQ3 = IQ4, IQL > IQR
        ITSCF = 5
      ELSEIF(IQ1.EQ.IQ2.AND.IQ3.EQ.IQ4.AND.IQL.EQ.IQR) THEN
C       IQ1 = IQ2, IQ3 = IQ4, IQL = IQR
        ITSCF = 6
      ELSE
C       ALL OTHER CASES GENERATED BY THE ABOVE, SO SKIP
        GOTO 3001
      ENDIF
C
C**********************************************************************C
C     LOOP OVER COMPONENT OVERLAP OPTIONS (INDEX 4000)                 C
C**********************************************************************C
C
C     LOOP OVER COMPONENT LABEL FOR A AND B: TT = LL (1) or SS (4)
      DO 4000 IT1=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITN(1) = IT1
C
C       FLAG READ-IN OF E0(AB) COEFFICIENTS FOR THIS COMPONENT LABEL
        IEAB = 1
C
C     LOOP OVER COMPONENT LABEL FOR C AND D: T'T' = LL (1) or SS (4)
      DO 4000 IT2=ITSTRT,ITSTOP,ITSKIP
C
C       PUT BLOCK VALUES INTO AN ARRAY FOR ERI ROUTINE LATER
        ITN(2) = IT2
C
C       FLAG READ-IN OF E0(CD) COEFFICIENTS FOR THIS COMPONENT LABEL
        IECD = 1
C
C     COMPONENT OVERLAP INDEX {(LL|LL)=1,(LL|SS)=2,(SS|LL)=3,(SS|SS)=4}
      ITT = MAPTTTT(IT1,IT2)
C
C     STAGE 1: INCLUDE ONLY (LL|LL) REPULSION INTEGRALS
C      IF(ILEV.EQ.1.AND.ITT.GT.1) THEN
C        GOTO 4001
C      ENDIF
C
C     STAGE 2: INCLUDE ONLY (LL|SS) AND (SS|LL) REPULSION INTEGRALS
C      IF(ILEV.EQ.2.AND.ITT.GT.3) THEN
C        GOTO 4001
C      ENDIF
C
C     STAGE 3: INCLUDE ONLY TWO-CENTRE (SS|SS) REPULSION INTEGRALS
C      IF(ILEV.EQ.3.AND.ITT.EQ.4) THEN
CC
C        IF(SSSSI3.AND.MCNT.EQ.3) THEN
C          GOTO 4001
C        ENDIF
C        IF(SSSSI4.AND.MCNT.EQ.4) THEN
C          GOTO 4001
C        ENDIF
C
C        OPTIONS FOR VISSCHER SMALL COMPONENT APPROXIMATION
C
c        IF(MCNT.EQ.2) THEN
cC         DIRECT CHOICE
c          IF(ICNTA.EQ.ICNTB.AND.ICNTC.EQ.ICNTD) THEN
c            GOTO 4002
cC         EXCHANGE CHOICE
c          ELSEIF(ICNTA.EQ.ICNTD.AND.ICNTB.EQ.ICNTC) THEN
c            GOTO 4002
c          ENDIF
c          GOTO 4001
cC         TODO: THIS IS WHERE I WOULD EVALUATE THE POINT-COULOMB RESULTS
c4002      CONTINUE
C        ENDIF
C
C      ENDIF
      
C      IF(ITT.NE.4) GOTO 4001
C
C     UPDATE COUNTER FOR NUMBER OF CLASSES
      N2EB(MCNT,ITT) = N2EB(MCNT,ITT)+1
C
C     HANDLE BREIT CASE WITH ITT=4
      IF(ITT.EQ.4) THEN
        N2EB(MCNT,5) = N2EB(MCNT,5)+1
        N2EB(MCNT,6) = N2EB(MCNT,6)+1
      ENDIF      
C
C**********************************************************************C
C     LOOP OVER BASIS FUNCTIONS (IBAS,JBAS) TO CONSTRUCT GMAT/QMAT     C
C**********************************************************************C
C
C     START OF PARALLEL REGION
C!$OMP PARALLEL DO COLLAPSE(2)
C!$OMP&  PRIVATE(RR,IBCH,IFLG)
C!$OMP&  SHARED(XYZ,KQN,MQN,EXL,NBAS,ITN)
      DO IBAS=1,NBAS(1)
        DO JBAS=1,NBAS(2)
C
C         SCHWARZ SCREENING (ECONOMIC ONLY WHEN SCREENING FRACTION BIG)
          IF(SCHWRZ.AND.ITT.GT.1) THEN
            ITOG = 1
          ELSE
            ITOG = 0
          ENDIF
C
c         CALL SCHWARZ(GDSC,SENS,TC2S)
          if(IT1.eq.1) then
            lamabc = lqn(1)+lqn(2)
          else
            lamabc = lqn(1)+lqn(2)+2
            lamabb = lqn(1)+lqn(2)+1
          endif
          
          if(IT2.eq.1) then
            lamcdc = lqn(3)+lqn(4)
          else
            lamcdc = lqn(3)+lqn(4)+2
            lamcdb = lqn(3)+lqn(4)+1
          endif
          ntuvabc =   (lamabc+1)*(lamabc+2)*(lamabc+3)/6
          ntuvcdc =   (lamcdc+1)*(lamcdc+2)*(lamcdc+3)/6

          ntuvabg = 3*(lamabb+1)*(lamabb+2)*(lamabb+3)/6
          ntuvcdg = 3*(lamcdb+1)*(lamcdb+2)*(lamcdb+3)/6

          ntuvabb = 9*(lamabb+1)*(lamabb+2)*(lamabb+3)/6
          ntuvcdb = 9*(lamcdb+1)*(lamcdb+2)*(lamcdb+3)/6

          ncntc = ntuvabc*(ntuvcdc+1)
          ncntg = ntuvabg*(ntuvcdg+1)
          ncntb = ntuvabb*(ntuvcdb+1)
c
          if(itt.eq.1) then
            lamabcdc = lqn(1)+lqn(2)+lqn(3)+lqn(4)
          elseif(itt.eq.2.or.itt.eq.3) then
            lamabcdc = lqn(1)+lqn(2)+lqn(3)+lqn(4)+2
          elseif(itt.eq.4) then
            lamabcdc = lqn(1)+lqn(2)+lqn(3)+lqn(4)+4
            lamabcdg = lqn(1)+lqn(2)+lqn(3)+lqn(4)+2
            lamabcdb = lqn(1)+lqn(2)+lqn(3)+lqn(4)+4
          endif
          ntuvc = (lamabcdc+1)*(lamabcdc+2)*(lamabcdc+3)/6
          ntuvg = (lamabcdg+1)*(lamabcdg+2)*(lamabcdg+3)/6
          ntuvb = (lamabcdb+1)*(lamabcdb+2)*(lamabcdb+3)/6
C
C         UPDATE COUNTER FOR NUMBER OF INTEGRALS AND SCREENED INTEGRALS
          N2EI(MCNT,ITT) = N2EI(MCNT,ITT)+NBAS(3)*NBAS(4)
          N2EM(MCNT,ITT) = N2EM(MCNT,ITT)+NBAS(3)*NBAS(4)*ntuvc
          N2EC(MCNT,ITT) = N2EC(MCNT,ITT)+NBAS(3)*NBAS(4)*ncntc
          N2ES(MCNT,ITT) = N2ES(MCNT,ITT)+NBAS(3)*NBAS(4)-MAXN
C
C         HANDLE BREIT CASE WITH ITT=4
          IF(ITT.EQ.4) THEN
            N2EI(MCNT,5) = N2EI(MCNT,5)+NBAS(3)*NBAS(4)
            N2EM(MCNT,5) = N2EM(MCNT,5)+NBAS(3)*NBAS(4)*ntuvg
            N2EC(MCNT,5) = N2EC(MCNT,5)+NBAS(3)*NBAS(4)*ncntg
            N2ES(MCNT,5) = N2ES(MCNT,5)+NBAS(3)*NBAS(4)-MAXN

            N2EI(MCNT,6) = N2EI(MCNT,6)+NBAS(3)*NBAS(4)
            N2EM(MCNT,6) = N2EM(MCNT,6)+NBAS(3)*NBAS(4)*ntuvb
            N2EC(MCNT,6) = N2EC(MCNT,6)+NBAS(3)*NBAS(4)*ncntb
            N2ES(MCNT,6) = N2ES(MCNT,6)+NBAS(3)*NBAS(4)-MAXN
          ENDIF      
C
C         CONDITIONAL TO SKIP THIS BATCH
          IF(IBCH.EQ.1) THEN
C
C           GENERATE BATCH OF ELECTRON REPULSION INTEGRALS
C           CALL ERI(RR,XYZ,KQN,MQN,NBAS,EXL,IBAS,JBAS,ITN)
C
          ENDIF
C
        ENDDO
      ENDDO
C     END OF PARALLEL REGION
C!$OMP END PARALLEL DO
C
C     RECORD TIME AT END OF BATCH
C
4001  CONTINUE
4000  CONTINUE
3001  CONTINUE
3000  CONTINUE
2001  CONTINUE
2000  CONTINUE
1001  CONTINUE
1000  CONTINUE
C
C**********************************************************************C
C     START OF SELF-CONSISTENT FIELD CALCULATIONS                      C
C**********************************************************************C
C
C     ADD TWO- AND MANY-CENTRE BINS TO (TT|TT) TOTALS
      DO MCNT=1,4
        DO ITT=1,6
          N2EB(5,ITT) = N2EB(5,ITT) + N2EB(MCNT,ITT)
          N2EI(5,ITT) = N2EI(5,ITT) + N2EI(MCNT,ITT)
          N2EM(5,ITT) = N2EM(5,ITT) + N2EM(MCNT,ITT)
          N2EC(5,ITT) = N2EC(5,ITT) + N2EC(MCNT,ITT)
          N2ES(5,ITT) = N2ES(5,ITT) + N2ES(MCNT,ITT)
        ENDDO
      ENDDO
C
C     ADD ALL (TT|TT) RESULTS TO TOTALS
      DO MCNT=1,5
        DO ITT=1,6
          N2EB(MCNT,7) = N2EB(MCNT,7) + N2EB(MCNT,ITT)
          N2EI(MCNT,7) = N2EI(MCNT,7) + N2EI(MCNT,ITT)
          N2EM(MCNT,7) = N2EM(MCNT,7) + N2EM(MCNT,ITT)
          N2EC(MCNT,7) = N2EC(MCNT,7) + N2EC(MCNT,ITT)
          N2ES(MCNT,7) = N2ES(MCNT,7) + N2ES(MCNT,ITT)
        ENDDO
      ENDDO
C
C     REDUCE TO PERCENTAGES
      DO MCNT=1,5
        DO ITT=1,7
          F2EB(MCNT,ITT) = 100*DFLOAT(N2EB(MCNT,ITT))/DFLOAT(N2EB(5,7))
          F2EI(MCNT,ITT) = 100*DFLOAT(N2EI(MCNT,ITT))/DFLOAT(N2EI(5,7))
          F2EM(MCNT,ITT) = 100*DFLOAT(N2EM(MCNT,ITT))/DFLOAT(N2EM(5,7))
          F2EC(MCNT,ITT) = 100*DFLOAT(N2EC(MCNT,ITT))/DFLOAT(N2EC(5,7))
        ENDDO
      ENDDO
C
C     FRACTION OF BLOCKS SCREENED
      DO MCNT=1,5
        DO ITT=1,6
          IF(N2EI(MCNT,ITT).NE.0) THEN
            RATIO = DFLOAT(N2ES(MCNT,ITT))/DFLOAT(N2EI(MCNT,ITT))
            F2ES(MCNT,ITT) = 100.0D0*RATIO
          ELSE
            F2ES(MCNT,ITT) = 0.0D0
          ENDIF
        ENDDO
      ENDDO
C
C     TWO-ELECTRON INTEGRAL SUMMARY
23    FORMAT(1X,A,3X,'1-centre',6X,'2-centre',6X,'3-centre',
     &                                         6X,'4-centre',9X,'Total')
24    FORMAT(1X,A,3X,I11,3X,I11,3X,I11,3X,I11,3X,I11)
25    FORMAT(1X,A,3X,F11.2,3X,F11.2,3X,F11.2,3X,F11.2,3X,F11.2)
c
      NLNES = 82
      WRITE(6, *) REPEAT('=',NLNES)
      WRITE(6, *) REPEAT(' ',24),'Two-electron calculation breakdown'
      WRITE(6, *) REPEAT('=',NLNES)
      WRITE(6,23) 'Blocks  (#)    '
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,24) '(LL|LL)     ',(N2EB(MCNT,1),MCNT=1,5)
      WRITE(6,24) '(LL|SS)     ',(N2EB(MCNT,2),MCNT=1,5)
      WRITE(6,24) '(SS|LL)     ',(N2EB(MCNT,3),MCNT=1,5)
      WRITE(6,24) '(SS|SS)     ',(N2EB(MCNT,4),MCNT=1,5)
      WRITE(6,24) '(SL|LS)gaunt',(N2EB(MCNT,5),MCNT=1,5)
      WRITE(6,24) '(SL|LS)gauge',(N2EB(MCNT,6),MCNT=1,5)
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,24) 'Total       ',(N2EB(MCNT,7),MCNT=1,5)

      WRITE(6, *) REPEAT('=',NLNES)
      WRITE(6,23) 'Integrals (M)  '
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,25) '(LL|LL)     ',(1.0D-6*N2EI(MCNT,1),MCNT=1,5)
      WRITE(6,25) '(LL|SS)     ',(1.0D-6*N2EI(MCNT,2),MCNT=1,5)
      WRITE(6,25) '(SS|LL)     ',(1.0D-6*N2EI(MCNT,3),MCNT=1,5)
      WRITE(6,25) '(SS|SS)     ',(1.0D-6*N2EI(MCNT,4),MCNT=1,5)
      WRITE(6,25) '(SL|LS)gaunt',(1.0D-6*N2EI(MCNT,5),MCNT=1,5)
      WRITE(6,25) '(SL|LS)gauge',(1.0D-6*N2EI(MCNT,6),MCNT=1,5)
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,25) 'Total       ',(1.0D-6*N2EI(MCNT,7),MCNT=1,5)
      WRITE(6, *) REPEAT('=',NLNES)

      WRITE(6,23) 'HGTFs olp (M)  '
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,25) '(LL|LL)     ',(1.0D-6*N2EM(MCNT,1),MCNT=1,5)
      WRITE(6,25) '(LL|SS)     ',(1.0D-6*N2EM(MCNT,2),MCNT=1,5)
      WRITE(6,25) '(SS|LL)     ',(1.0D-6*N2EM(MCNT,3),MCNT=1,5)
      WRITE(6,25) '(SS|SS)     ',(1.0D-6*N2EM(MCNT,4),MCNT=1,5)
      WRITE(6,25) '(SL|LS)gaunt',(1.0D-6*N2EM(MCNT,5),MCNT=1,5)
      WRITE(6,25) '(SL|LS)gauge',(1.0D-6*N2EM(MCNT,6),MCNT=1,5)
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,25) 'Total       ',(1.0D-6*N2EM(MCNT,7),MCNT=1,5)
      WRITE(6, *) REPEAT('=',NLNES)

      WRITE(6,23) 'contractions(B)'
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,25) '(LL|LL)     ',(1.0D-9*N2EC(MCNT,1),MCNT=1,5)
      WRITE(6,25) '(LL|SS)     ',(1.0D-9*N2EC(MCNT,2),MCNT=1,5)
      WRITE(6,25) '(SS|LL)     ',(1.0D-9*N2EC(MCNT,3),MCNT=1,5)
      WRITE(6,25) '(SS|SS)     ',(1.0D-9*N2EC(MCNT,4),MCNT=1,5)
      WRITE(6,25) '(SL|LS)gaunt',(1.0D-9*N2EC(MCNT,5),MCNT=1,5)
      WRITE(6,25) '(SL|LS)gauge',(1.0D-9*N2EC(MCNT,6),MCNT=1,5)
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,25) 'Total       ',(1.0D-9*N2EC(MCNT,7),MCNT=1,5)
      WRITE(6, *) REPEAT('=',NLNES)

      WRITE(6,23) 'Blocks  (%)    '
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,25) '(LL|LL)     ',(F2EB(MCNT,1),MCNT=1,5)
      WRITE(6,25) '(LL|SS)     ',(F2EB(MCNT,2),MCNT=1,5)
      WRITE(6,25) '(SS|LL)     ',(F2EB(MCNT,3),MCNT=1,5)
      WRITE(6,25) '(SS|SS)     ',(F2EB(MCNT,4),MCNT=1,5)
      WRITE(6,25) '(SL|LS)gaunt',(F2EB(MCNT,5),MCNT=1,5)
      WRITE(6,25) '(SL|LS)gauge',(F2EB(MCNT,6),MCNT=1,5)
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,25) 'Total       ',(F2EB(MCNT,7),MCNT=1,5)

      WRITE(6, *) REPEAT('=',NLNES)
      WRITE(6,23) 'Integrals(#)   '
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,25) '(LL|LL)     ',(F2EI(MCNT,1),MCNT=1,5)
      WRITE(6,25) '(LL|SS)     ',(F2EI(MCNT,2),MCNT=1,5)
      WRITE(6,25) '(SS|LL)     ',(F2EI(MCNT,3),MCNT=1,5)
      WRITE(6,25) '(SS|SS)     ',(F2EI(MCNT,4),MCNT=1,5)
      WRITE(6,25) '(SL|LS)gaunt',(F2EI(MCNT,5),MCNT=1,5)
      WRITE(6,25) '(SL|LS)gauge',(F2EI(MCNT,6),MCNT=1,5)
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,25) 'Total       ',(F2EI(MCNT,7),MCNT=1,5)
      WRITE(6, *) REPEAT('=',NLNES)

      WRITE(6,23) 'HGTFs olp(%)   '
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,25) '(LL|LL)     ',(F2EM(MCNT,1),MCNT=1,5)
      WRITE(6,25) '(LL|SS)     ',(F2EM(MCNT,2),MCNT=1,5)
      WRITE(6,25) '(SS|LL)     ',(F2EM(MCNT,3),MCNT=1,5)
      WRITE(6,25) '(SS|SS)     ',(F2EM(MCNT,4),MCNT=1,5)
      WRITE(6,25) '(SL|LS)gaunt',(F2EM(MCNT,5),MCNT=1,5)
      WRITE(6,25) '(SL|LS)gauge',(F2EM(MCNT,6),MCNT=1,5)
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,25) 'Total       ',(F2EM(MCNT,7),MCNT=1,5)
      WRITE(6, *) REPEAT('=',NLNES)

      WRITE(6,23) 'contractions   '
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,25) '(LL|LL)     ',(F2EC(MCNT,1),MCNT=1,5)
      WRITE(6,25) '(LL|SS)     ',(F2EC(MCNT,2),MCNT=1,5)
      WRITE(6,25) '(SS|LL)     ',(F2EC(MCNT,3),MCNT=1,5)
      WRITE(6,25) '(SS|SS)     ',(F2EC(MCNT,4),MCNT=1,5)
      WRITE(6,25) '(SL|LS)gaunt',(F2EC(MCNT,5),MCNT=1,5)
      WRITE(6,25) '(SL|LS)gauge',(F2EC(MCNT,6),MCNT=1,5)
      WRITE(6, *) REPEAT('-',NLNES)
      WRITE(6,25) 'Total       ',(F2EC(MCNT,7),MCNT=1,5)
      WRITE(6, *) REPEAT('=',NLNES)
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
      WRITE(6, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',22),'Eq-coefficient word analysis'
      WRITE(6, *) REPEAT('=',72)
      WRITE(6,20) 'Type','Λ','Terms','Length','#','Words','Size (MB)'
C
C     E0LL ANALYSIS
      WRITE(6, *) REPEAT('-',72)
      DO LAM=0,LAMMX
        IF(LAM.EQ.0) THEN
          WRITE(6,21) 'E0LL',LAM,NTUV(LAM,1),NTRM(LAM,1), 4,NWRD(LAM,1),
     &                                                      SPCE(LAM,1)
        ELSE
          WRITE(6,21) '    ',LAM,NTUV(LAM,1),NTRM(LAM,1), 4,NWRD(LAM,1),
     &                                                      SPCE(LAM,1)
        ENDIF
      ENDDO
C
      IF(HMLT.EQ.'NORL') GOTO 100
C
C     E0SS ANALYSIS
      WRITE(6, *) REPEAT('-',72)
      DO LAM=2,LAMMX+2
        IF(LAM.EQ.2) THEN
          WRITE(6,21) 'E0SS',LAM,NTUV(LAM,2),NTRM(LAM,2), 4,NWRD(LAM,2),
     &                                                      SPCE(LAM,2)
        ELSE
          WRITE(6,21) '    ',LAM,NTUV(LAM,2),NTRM(LAM,2), 4,NWRD(LAM,2),
     &                                                      SPCE(LAM,2)
        ENDIF
      ENDDO
C
      IF(HMLT.EQ.'BARE'.OR.HMLT.EQ.'DHFR') GOTO 100
C
C     E0SS ANALYSIS
      WRITE(6, *) REPEAT('-',72)
      DO LAM=1,LAMMX+1
        IF(LAM.EQ.1) THEN
          WRITE(6,21) 'EILS',LAM,NTUV(LAM,3),NTRM(LAM,3),12,NWRD(LAM,3),
     &                                                      SPCE(LAM,3)
        ELSE
          WRITE(6,21) '    ',LAM,NTUV(LAM,3),NTRM(LAM,3),12,NWRD(LAM,3),
     &                                                      SPCE(LAM,3)
        ENDIF
      ENDDO
C
      IF(TREE.NE.'MBPTN') GOTO 100
C
C     EISL ANALYSIS
      WRITE(6, *) REPEAT('-',72)
      DO LAM=1,LAMMX+1
        IF(LAM.EQ.1) THEN
          WRITE(6,21) 'EISL',LAM,NTUV(LAM,4),NTRM(LAM,4),12,NWRD(LAM,4),
     &                                                      SPCE(LAM,4)
        ELSE
          WRITE(6,21) '    ',LAM,NTUV(LAM,4),NTRM(LAM,4),12,NWRD(LAM,4),
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
      WRITE(6,24) 'Total','Λ_max','Terms','Length','#',
     &                                              'Words','Size (MB)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(6,23) 'E0LL',LAMMX  ,NTUVT(1),NTRMT(1), 4,NWRDT(1),SPCET(1)
      IF(HMLT.EQ.'NORL') GOTO 200
      WRITE(6,23) 'E0SS',LAMMX+2,NTUVT(2),NTRMT(2), 4,NWRDT(2),SPCET(2)
      IF(HMLT.EQ.'DHFR') GOTO 200
      WRITE(6,23) 'EILS',LAMMX+1,NTUVT(3),NTRMT(3),12,NWRDT(3),SPCET(3)
      IF(TREE.NE.'MBPTN') GOTO 200
      IF(HMLT.EQ.'NORL'.OR.HMLT.EQ.'DHFR') GOTO 200
      WRITE(6,23) 'EISL',LAMMX+1,NTUVT(4),NTRMT(4),12,NWRDT(4),SPCET(4)
200   CONTINUE
C      
      WRITE(6, *) REPEAT('-',72)
      WRITE(6,25) '       ',NWRDNET,SPCENET
      WRITE(6, *) REPEAT('-',72)
      WRITE(6,22) 'parameters.h',MFL,MWRD,MWRD*MFL,SPCEMFL
      WRITE(6, *) REPEAT('=',72)
C
C     OPTION WHEN NUMBER OF WORDS EXCEEDS ALLOCATED SIZE LIMIT
      IF(NTRMT(1).GT.MFL) THEN
        WRITE(6, *) 'In EQSAVE: E0LL words exceed allocated limit.'
        GOTO 150
      ENDIF
C
      IF(HMLT.NE.'NORL') THEN
        IF(NTRMT(2).GT.MFL) THEN
          WRITE(6, *) 'In EQSAVE: E0SS words exceed allocated limit.'
          GOTO 150
        ENDIF
      ENDIF
C     NO NEED TO ASK ABOUT EILS OR EISL! E0SS IS LARGER THAN EITHER.
C
      GOTO 250
C
      WRITE(6, *) REPEAT('-',72)
C
C     HAVE TO GENERATE E COEFFICIENTS BY BATCH
      WRITE(6, *) 'In EQSAVE: Eq-coefficients to be generated by batch.'
C
C     ONE OF THE CLASSES EXCEEDS WORD LIMIT
150   CONTINUE
C
C     END OF SECTION
      WRITE(6, *) REPEAT('=',72)
      WRITE(6, *) ' '
C
250   CONTINUE
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
C
      END
