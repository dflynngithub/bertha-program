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
C             aa   tttttttt oooooo  mm       mm iiii cccccc            C
C            aaaa     tt   oo    oo mmm     mmm  ii cc    cc           C
C           aa  aa    tt   oo    oo mmmm   mmmm  ii cc                 C
C          aa    aa   tt   oo    oo mm mm mm mm  ii cc                 C
C          aaaaaaaa   tt   oo    oo mm  mmm  mm  ii cc                 C
C          aa    aa   tt   oo    oo mm   m   mm  ii cc    cc           C
C          aa    aa   tt    oooooo  mm       mm iiii cccccc            C
C                                                                      C
C ******************************************************************** C
C                                                                      C
C    THIS MODIFIED VERSION OF BERTHA JUST INCLUDES THE ATOMIC CODE.    C
C    (FACILITATES BIGGER BASIS SETS -- GOOD FOR SUPERHEAVY TESTING).   C
C                                                                      C
C -------------------------------------------------------------------- C
C        A RELATIVISTIC MOLECULAR ELECTRONIC STRUCTURE PROGRAM         C
C            BASED ON THE ANALYTIC FINITE BASIS SET METHOD.            C
C                                                                      C
C                 H.M.QUINEY, H.SKAANE (OXFORD, 1997)                  C
C                       D. FLYNN (UNIMELB, 2020)                       C
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
C                          TABLE OF CONTENTS                           C
C                          -----------------                           C
C   [1] INPUT/OUTPUT: READ FROM INPUT FILE AND SUMMARISE DATA.         C
C   [2] LABELS: MOLECULAR GEOMETRY AND FOCK MATRIX BLOCK ADDRESSES.    C
C   [3] DENSITIES: MOLECULAR DENSITIES, ENERGIES AND LEVEL SHIFTING.   C
C   [4] ATOMIC HARTREE-FOCK: AVERAGE OF CONFIG. ATOMIC SCF.            C
C   [5] QED: UEHLING ROUTINE DEVELOPMENT.                              C
C                                                                      C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/TCPU/TBEG,TNUC,TATM,TSCF,TMPT,TPRP,TPLT,TTOT
C
C     READ DATA FROM USER-SPECIFIED INPUT FILE
      CALL SYSTEM_CLOCK(ICL1,RATE)
      TBEG = DFLOAT(ICL1)/RATE
      CALL CARDIN
C
C     OPEN FILE FOR TERMINAL RECORD
      OPEN(UNIT=7,FILE=TRIM(OUTFL)//'_atomic.out',STATUS='UNKNOWN')
C
C     PRINT SUMMARY OF INPUT DATA
      CALL INPUT
C
C     NUCLEAR POTENTIALS
      CALL SYSTEM_CLOCK(ICL3,RATE)
      CALL NUCCOUL
      CALL NUCVCPL
      CALL SYSTEM_CLOCK(ICL4)
      TNUC = DFLOAT(ICL4-ICL3)/RATE
C
C     FOCK MATRIX SYMMETRY TYPE INDICES
      CALL FOCKIND
C
C     ATOMIC HARTREE-FOCK SCF ROUTINE
C     IF(.NOT.READIN) THEN
        CALL SYSTEM_CLOCK(ICL3,RATE)
        CALL ATOMIC
        CALL SYSTEM_CLOCK(ICL4)
        TATM = DFLOAT(ICL4-ICL3)/RATE
C     ENDIF
C
C     PLOT THE RADIAL WAVEFUNCTION NEAR THE NUCLEUS
C      CALL RADPLOT(1,0.0D0,5.0D0/CFM,1000)
C      CALL RATPLOT(1,0.0D0,5.0D0/CFM,1000)
C      CALL RADPLOT(1,0.01D0,0.10D0,1000)
C      CALL RATPLOT(1,0.01D0,0.10D0,1000)
C
C     PRINT SUMMARY OF OUTPUT DATA
C
C     OVERALL CALCULATION TIME
      CALL SYSTEM_CLOCK(ICL2,RATE)
      TTOT = DFLOAT(ICL2-ICL1)/RATE
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
      CHARACTER*7  HMINT(10),PTYPE(10)
      CHARACTER*9  BTYP
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 COEF(MDM,MDM),OVLP(MDM,MDM)
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
C
        WRITE(6, *) 'In CARDIN: invalid value NOPN.',NOPN
        WRITE(7, *) 'In CARDIN: invalid value NOPN.',NOPN
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
      IF(ILEV.LT.1.AND.ILEV.GT.5) THEN
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
      ELSEIF(ILEV.EQ.4.OR.ILEV.EQ.5) THEN
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
        IF(NHMINT.LT.1.OR.NHMINT.GT.10) THEN
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
      IF(READIN) THEN
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
      ENDIF
C
C     INITIALISE COEFFICIENT AND OVERLAP MATRICES IF NO FILE SPECIFIED
      IF(.NOT.READIN) THEN
        DO I=1,NDIM
          DO J=1,NDIM
            COEF(I,J) = DCMPLX(0.0D0,0.0D0)
            OVLP(I,J) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
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
      WRITE(6, *) REPEAT(' ',29),'BERTHA-ATOMIC'
      WRITE(7, *) REPEAT(' ',29),'BERTHA-ATOMIC'
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
      CALL SYSTEM_CLOCK(ICL3,RATE)
      TCPUS(1) = DFLOAT(ICL3)/RATE
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
      WRITE(6,20) 'Atomic Hartree-Fock SCF:      ',HMS(TATM)
      WRITE(7,20) 'Atomic Hartree-Fock SCF:      ',HMS(TATM)
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
30    FORMAT(1X,A,19(' '),A)
      WRITE(6, *)
      WRITE(7, *)
      WRITE(6,30) 'Successful BERTHA-ATOMIC exit at:',STAMP
      WRITE(7,30) 'Successful BERTHA-ATOMIC exit at:',STAMP
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
c     DATA CV/137.03599976D0/
C
C     SPEED OF LIGHT (USED IN HAAKON'S THESIS)
C     DATA CV/137.0359898D0/
C
C     SPEED OF LIGHT (GRASP 1992)
C     DATA CV/137.03599976D0/
C
C     SPEED OF LIGHT (PARPIA AND MOHANTY 1992)
C     DATA CV/137.0359895D0/
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
C     WRITE(6,*) 'αβγδεζηθικλμνξοπρϱςστυφχψω'
C     WRITE(6,*) 'ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡ  ΣΤΥΦΧΨΩ'
C     WRITE(6,*) 'ⅈℵℎℏℓℒℋℛℐℜℑℕℙℚℝℤℂ℘'
C     WRITE(6,*) '×·∕−±∓≠∞—–…≪≫≤≥'
C     WRITE(6,*) '⅟½⅓⅕⅙⅛⅔⅖⅚⅜¾⅗⅝⅞⅘'
C     WRITE(6,*) 'π∞Σ√∛'∜∫∬∭∮∯∰∀∂∃∄∅∆∇∈∉∋∌∎∏∐∑'
C     WRITE(6,*) '∗∝∠∡∢∧∨∩∪∴∵∶∷∼∽∿⊕⊖⊗⊘⊙⊚⊛⊜⊝␢Å'
C     WRITE(6,*) '‗‘’‚‛“”„‟•‣❛❜❝❞¿'
C     WRITE(6,*) '✓✔✗✘☓√☑☐☒☼❄❆♤♠♧♣♡♥♢♦'
C     WRITE(6,*) '¢$€£¥®™☎⌨✁✂✎✏✐§☛'
C     WRITE(6,*) '↕↖↗↘↙↚↛↥↦↧↰↱↲↳↴↺↻↼↽↾↿⇀⇁⇂⇃⇄⇅⇆⇇⇠⇡⇢⏎▶➔➘➙➚➛➜➞↵⇑⇓'
C     WRITE(6,*) '★☆✡✦✧✩✪✰✢✣✤✥✱✲✳✴✵✶✷✸✹✺✻✼✽✾✿❀❁❂❃❇❈❉❊❋❄❆❅≛'
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
C   [F] CARTIND: GENERATES INDICES FOR EQ-COEFFS AND R-INTEGRALS.      C
C   [G] AUFBAU: DETERMINES GROUND STATE ATOMIC ELECTRON CONFIG.        C
C   [H] SPECTRM0: ATOMIC SPECTRUM W/ EIGENVALUES AND RADIAL MOMENTS.   C
C   [I] SPECTRM: MOLECULAR SPECTRUM W/ EIGENVALUES AND TERM SYMBOLS.   C
C   [J] LLAB: GIVES THE CHARACTER CORRESPONDING TO LQN VALUE LQN.      C
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
C     GAUSSIAN FITTING FOR ALL NUCLEI
      DO IZ=1,NCNT
C
C       NUCLEAR RADIUS
        IF(INT(ZNUC(IZ)).EQ.1) THEN
C         PROTON HAS A VERY SMALL RADIUS (IN FM)
          RFM = PRAD
        ELSEIF(ZNUC(IZ).GT.1.0D0.AND.ZNUC(IZ).LT.90.0D0) THEN
C         EMPIRICAL FORMULA FOR NUCLEAR RMS RADIUS IF Z < 90 (IN FM)
          RFM = 0.836D0*(ANUC(IZ)**(1.0D0/3.0D0)) + 0.57D0
C         RFM = 1.20D0*(ANUC(IZ)**(1.0D0/3.0D0))
        ELSE
C         EMPIRICAL FORMULA FOR NUCLEAR RMS RADIUS IF Z > 90 (IN FM)
          RFM = 0.770D0*(ANUC(IZ)**(1.0D0/3.0D0)) + 0.98D0
        ENDIF
C
C       CONVERT RESULT TO ATOMIC UNITS
        RNUC(IZ) = RFM/CFM
C
C       GAUSSIAN MODEL
        NNUC(IZ)   = 0
        FNUC(IZ,0) = 1.0D0
        XNUC(IZ,0) = 1.5D0/(RNUC(IZ)*RNUC(IZ))
C
C       FERMI MODEL
        TFMI(IZ) = 2.30D0/CFM
        AFMI(IZ) = 0.25D0*TFMI(IZ)/THLG
C
      ENDDO
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
C          LOCK OUT ACCESS TO FUNCTIONS FOR GAUSSIAN NUCLEUS
           IF(NMDL(IZ).EQ.'GAUSS') THEN
             NNUC(IZ) = 0
           ENDIF
C
C         FERMI MODEL PARAMETERS
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
      COMMON/BUEE/RUEE(MCT,3),FUEE(MCT,MFT),XUEE(MCT,MFT),NUEE(MCT)
      COMMON/BUEU/RUEU(MCT,3),FUEU(MCT,MFT),XUEU(MCT,MFT),NUEU(MCT)
      COMMON/BUEP/RUEP(MCT,3),FUEP(MCT,MFT),XUEP(MCT,MFT),NUEP(MCT)
      COMMON/BWKR/RWKR(MCT,3),FWKR(MCT,MFT),XWKR(MCT,MFT),NWKR(MCT)
      COMMON/BKSB/RKSB(MCT,3),FKSB(MCT,MFT),XKSB(MCT,MFT),NKSB(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/QUEE/RE(0:NRAD),VE(MCT,0:NRAD),R0E,RME,RIE,NLE,NEE
      COMMON/QUEU/RU(0:NRAD),VU(MCT,0:NRAD),R0U,RMU,RIU,NLU,NEU
      COMMON/QUEP/RP(0:NRAD),VP(MCT,0:NRAD),R0P,RMP,RIP,NLP,NEP
      COMMON/QWKR/RW(0:NRAD),VW(MCT,0:NRAD),R0W,RMW,RIW,NLZ,NEW
      COMMON/QKSB/RK(0:NRAD),VK(MCT,0:NRAD),R0K,RMK,RIK,NLK,NEK
C
C     EARLY EXIT IF VACPOL MEAN FIELDS AREN'T NEEDED
      IF(HMLT.NE.'DHFQ'.AND.HMLT.NE.'DHFP') RETURN
C
C**********************************************************************C
C     ORDER (α/π) (Zα) : UEHLING POTENTIAL (e+e-)                      C
C**********************************************************************C
C
C     UEHLING FITTING SET
      VACFL = 'output/'//TRIM(MOLCL)//'_vpuee.dat'
      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
C
      IF(FILETHERE) THEN
C
C       READ FITTING DATA IN FROM AN EXISTING FILE
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          DO IZ=1,NCNT
            READ(8, *) NUEE(IZ)
            DO IFT=1,NUEE(IZ)
              READ(8, *) FUEE(IZ,IFT),XUEE(IZ,IFT)
            ENDDO
          ENDDO
        CLOSE(UNIT=8)
C
      ELSE
C
C       ALERT THE USER THAT THE UEHLING FILE COULD NOT BE FOUND
        WRITE(6, *) 'In NUCVCPL: could not find Uehling (e+e-) data.'
        WRITE(7, *) 'In NUCVCPL: could not find Uehling (e+e-) data.'
        DO IZ=1,NCNT
          NUEE(IZ) = 0
        ENDDO
C
      ENDIF
C
C     UEHLING RADIAL GRID DATA
      VACFL = 'output/'//TRIM(MOLCL)//'_vpuee_grid.dat'
      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
C
      IF(FILETHERE) THEN
C
C       READ FITTING DATA IN FROM AN EXISTING FILE
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          READ(8, *) R0E,REU,RIE,NLE,NEE
          DO N=0,NRAD
            READ(8, *) RE(N),(VE(IZ,N),IZ=1,NCNT)
          ENDDO
        CLOSE(UNIT=8)
C
      ELSE
C
C       ALERT THE USER THAT THE UEH GRID FILE COULD NOT BE FOUND
        WRITE(6, *) 'In NUCVCPL: could not find UEH (e+e-) grid data.'
        WRITE(7, *) 'In NUCVCPL: could not find UEH (e+e-) grid data.'
        DO N=0,NRAD
          RE(N) = 0.0D0
          DO IZ=1,NCNT
            VE(IZ,N) = 0.0D0
          ENDDO
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C     ORDER (α/π) (Zα) : UEHLING POTENTIAL (μ+μ-)                      C
C**********************************************************************C
C
C     UEHLING FITTING SET
      VACFL = 'output/'//TRIM(MOLCL)//'_vpueu.dat'
      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
C
      IF(FILETHERE) THEN
C
C       READ FITTING DATA IN FROM AN EXISTING FILE
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          DO IZ=1,NCNT
            READ(8, *) NUEU(IZ)
            DO IFT=1,NUEU(IZ)
              READ(8, *) FUEU(IZ,IFT),XUEU(IZ,IFT)
            ENDDO
          ENDDO
        CLOSE(UNIT=8)
C
      ELSE
C
C       ALERT THE USER THAT THE UEHLING FILE COULD NOT BE FOUND
        WRITE(6, *) 'In NUCVCPL: could not find Uehling (μ+μ-) data.'
        WRITE(7, *) 'In NUCVCPL: could not find Uehling (μ+μ-) data.'
        DO IZ=1,NCNT
          NUEU(IZ) = 0
        ENDDO
C
      ENDIF
C
C     UEHLING RADIAL GRID DATA
      VACFL = 'output/'//TRIM(MOLCL)//'_vpueu_grid.dat'
      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
C
      IF(FILETHERE) THEN
C
C       READ FITTING DATA IN FROM AN EXISTING FILE
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          READ(8, *) R0U,RMU,RIU,NLU,NEU
          DO N=0,NRAD
            READ(8, *) RU(N),(VU(IZ,N),IZ=1,NCNT)
          ENDDO
        CLOSE(UNIT=8)
C
      ELSE
C
C       ALERT THE USER THAT THE UEH GRID FILE COULD NOT BE FOUND
        WRITE(6, *) 'In NUCVCPL: could not find UEH (μ+μ-) grid data.'
        WRITE(7, *) 'In NUCVCPL: could not find UEH (μ+μ-) grid data.'
        DO N=0,NRAD
          RU(N) = 0.0D0
          DO IZ=1,NCNT
            VP(IZ,N) = 0.0D0
          ENDDO
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C     ORDER (α/π) (Zα) : UEHLING POTENTIAL (π+π-)                      C
C**********************************************************************C
C
C     UEHLING FITTING SET
      VACFL = 'output/'//TRIM(MOLCL)//'_vpuep.dat'
      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
C
      IF(FILETHERE) THEN
C
C       READ FITTING DATA IN FROM AN EXISTING FILE
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          DO IZ=1,NCNT
            READ(8, *) NUEP(IZ)
            DO IFT=1,NUEP(IZ)
              READ(8, *) FUEP(IZ,IFT),XUEP(IZ,IFT)
            ENDDO
          ENDDO
        CLOSE(UNIT=8)
C
      ELSE
C
C       ALERT THE USER THAT THE UEHLING FILE COULD NOT BE FOUND
        WRITE(6, *) 'In NUCVCPL: could not find Uehling (π+π-) data.'
        WRITE(7, *) 'In NUCVCPL: could not find Uehling (π+π-) data.'
        DO IZ=1,NCNT
          NUEP(IZ) = 0
        ENDDO
C
      ENDIF
C
C     UEHLING RADIAL GRID DATA
      VACFL = 'output/'//TRIM(MOLCL)//'_vpuep_grid.dat'
      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
C
      IF(FILETHERE) THEN
C
C       READ FITTING DATA IN FROM AN EXISTING FILE
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          READ(8, *) R0P,RMP,RIP,NLP,NEP
          DO N=0,NRAD
            READ(8, *) RP(N),(VU(IZ,N),IZ=1,NCNT)
          ENDDO
        CLOSE(UNIT=8)
C
      ELSE
C
C       ALERT THE USER THAT THE UEH GRID FILE COULD NOT BE FOUND
        WRITE(6, *) 'In NUCVCPL: could not find UEH (π+π-) grid data.'
        WRITE(7, *) 'In NUCVCPL: could not find UEH (π+π-) grid data.'
        DO N=0,NRAD
          RP(N) = 0.0D0
          DO IZ=1,NCNT
            VU(IZ,N) = 0.0D0
          ENDDO
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C     ORDER (α/π) (Zα)³: WICHMANN-KROLL POTENTIAL                      C
C**********************************************************************C
C
C     WICHMANN-KROLL FITTING SET
      VACFL = 'output/'//TRIM(MOLCL)//'_vpwkr.dat'
      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
C
      IF(FILETHERE) THEN
C
C       READ FITTING DATA IN FROM AN EXISTING FILE
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          DO IZ=1,NCNT
            READ(8, *) NWKR(IZ)
            DO IFT=1,NWKR(IZ)
              READ(8, *) FWKR(IZ,IFT),XWKR(IZ,IFT)
            ENDDO
          ENDDO
        CLOSE(UNIT=8)
C
      ELSE
C
C       ALERT THE USER THAT THE WICHMANN-KROLL FILE COULD NOT BE FOUND
        WRITE(6, *) 'In NUCVCPL: could not find Wichmann-Kroll data.'
        WRITE(7, *) 'In NUCVCPL: could not find Wichmann-Kroll data.'
        DO IZ=1,NCNT
          NWKR(IZ) = 0
        ENDDO
C
      ENDIF
C
C     WICHMANN-KROLL RADIAL GRID DATA
      VACFL = 'output/'//TRIM(MOLCL)//'_vpwkr_grid.dat'
C
      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
      IF(FILETHERE) THEN
C
C       READ FITTING DATA IN FROM AN EXISTING FILE
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          READ(8, *) R0W,RMW,RIW,NLZ,NEW
          DO N=0,NRAD
            READ(8, *) RW(N),(VW(IZ,N),IZ=1,NCNT)
          ENDDO
        CLOSE(UNIT=8)
C
      ELSE
C
C       ALERT THE USER THAT THE WKR GRID FILE COULD NOT BE FOUND
        WRITE(6, *) 'In NUCVCPL: could not find WKR grid data.'
        WRITE(7, *) 'In NUCVCPL: could not find WKR grid data.'
        DO N=0,NRAD
          RW(N) = 0.0D0
          DO IZ=1,NCNT
            VW(IZ,N) = 0.0D0
          ENDDO
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C     ORDER (α/π)²(Zα) : KÄLLÉN-SABRY POTENTIAL                        C
C**********************************************************************C
C
C     KÄLLÉN-SABRY FITTING SET
      VACFL = 'output/'//TRIM(MOLCL)//'_vpksb.dat'
      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
C
      IF(FILETHERE) THEN
C
C       READ FITTING DATA IN FROM AN EXISTING FILE
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          DO IZ=1,NCNT
            READ(8, *) NKSB(IZ)
            DO IFT=1,NKSB(IZ)
              READ(8, *) FKSB(IZ,IFT),XKSB(IZ,IFT)
            ENDDO
          ENDDO
        CLOSE(UNIT=8)
C
      ELSE
C
C       ALERT THE USER THAT THE KÄLLÉN-SABRY FILE COULD NOT BE FOUND
        WRITE(6, *) 'In NUCVCPL: could not find Källén-Sabry data.'
        WRITE(7, *) 'In NUCVCPL: could not find Källén-Sabry data.'
        DO IZ=1,NCNT
          NKSB(IZ) = 0
        ENDDO
C
      ENDIF
C
C     KÄLLÉN-SABRY RADIAL GRID DATA
      VACFL = 'output/'//TRIM(MOLCL)//'_vpksb_grid.dat'
      INQUIRE(FILE=TRIM(VACFL),EXIST=FILETHERE)
C
      IF(FILETHERE) THEN
C
C       READ FITTING DATA IN FROM AN EXISTING FILE
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          READ(8, *) R0K,RMK,RIK,NLK,NEK
          DO N=0,NRAD
            READ(8, *) RK(N),(VK(IZ,N),IZ=1,NCNT)
          ENDDO
        CLOSE(UNIT=8)
C
      ELSE
C
C       ALERT THE USER THAT THE KSB GRID FILE COULD NOT BE FOUND
        WRITE(6, *) 'In NUCVCPL: could not find KSB grid data.'
        WRITE(7, *) 'In NUCVCPL: could not find KSB grid data.'
        DO N=0,NRAD
          RK(N) = 0.0D0
          DO IZ=1,NCNT
            VK(IZ,N) = 0.0D0
          ENDDO
        ENDDO
C
      ENDIF
C
C     NUMBER OF DATA POINTS IN UNIFORMLY-SPACED AND EXPONENTIAL REGION
      NLIN = NRAD/10 - MOD(NRAD/10,10)
      NEXP = NRAD-NLIN
C
C     GENERATE RADIAL GRID (LINEAR FROM FEMTOMETERS, EXPONENTIAL IN AU)
      RORI =  0.0D0/CFM
      RMID = 10.0D0/CFM
      RMAX =  0.5D0
C
      HL = (RMID-RORI)/DFLOAT(NLIN)
      HE = DLOG(RMAX/RMID)/DFLOAT(NEXP)
C
      RETURN
      END
C
C
      FUNCTION VNGAUSS(R,XI,FRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    VV    VV NN    NN  GGGGGG     AA    UU    UU  SSSSSS   SSSSSS     C
C    VV    VV NNN   NN GG    GG   AAAA   UU    UU SS    SS SS    SS    C
C    VV    VV NNNN  NN GG        AA  AA  UU    UU SS       SS          C
C    VV    VV NN NN NN GG       AA    AA UU    UU  SSSSSS   SSSSSS     C
C     VV  VV  NN  NNNN GG   GGG AAAAAAAA UU    UU       SS       SS    C
C      VVVV   NN   NNN GG    GG AA    AA UU    UU SS    SS SS    SS    C
C       VV    NN    NN  GGGGGG  AA    AA  UUUUUU   SSSSSS   SSSSSS     C
C                                                                      C
C -------------------------------------------------------------------- C
C  VNGAUSS EVALUATES THE POTENTIAL ARISING FROM ONE TERM IN THE SET OF C
C  GAUSSIAN DISTRIBUTIONS AT A PARTICULAR RADIUS R.                    C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
      IF(R.LT.1.0D-10) THEN
        XR2 = XI*R*R
        V0  = 2.0D0-2.0D0*XR2/3.0D0+XR2*XR2/5.0D0-XR2*XR2*XR2/21.0D0
        VNGAUSS = -V0*DSQRT(XI)/PI12
      ELSE
        VNGAUSS = -DERF(DSQRT(XI)*R)/R
      ENDIF
C
C     FINAL VALUE FOR THE POTENTIAL
      VNGAUSS = FRAC*VNGAUSS
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
C     FERMI NUCLEAR CHARGE MODEL
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
      FUNCTION ERNUC(MODEL,IZ,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             EEEEEEEE RRRRRRR  NN    NN UU    UU  CCCCCC              C
C             EE       RR    RR NNN   NN UU    UU CC    CC             C
C             EE       RR    RR NNNN  NN UU    UU CC                   C
C             EEEEEE   RR    RR NN NN NN UU    UU CC                   C
C             EE       RRRRRRR  NN  NNNN UU    UU CC                   C
C             EE       RR    RR NN   NNN UU    UU CC    CC             C
C             EEEEEEEE RR    RR NN    NN  UUUUUU   CCCCCC              C
C                                                                      C
C -------------------------------------------------------------------- C
C  ERNUC EVALUATES THE VALUE OF THE NORMALISED NUCLEAR ELECTRIC FIELD  C
C  (RADIAL COMPONENT) AT RADIUS R, BASED ON A NUCLEAR CHARGE MODEL.    C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5 NMDL
      CHARACTER*5 MODEL
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     RADIUS MUST BE NON-NEGATIVE
      IF(R.LT.0.0D0) THEN
        WRITE(6, *) 'In ERNUC: radius R must be non-negative. R = ',R
        WRITE(7, *) 'In ERNUC: radius R must be non-negative. R = ',R
        STOP
      ENDIF
C
C     REASONS TO USE GAUSSIAN NUCLEAR CHARGE MODEL
      EPS = DSQRT(1.4D0)*PI*0.25D0*TFMI(IZ)/THLG
      IF(EPS.GT.RNUC(IZ)) THEN
C
        XI = 3.0D0/(2.0D0*RNUC(IZ)*RNUC(IZ))
        IF(R.LT.0.2D0*RNUC(IZ)) THEN
          X1 = DSQRT(XI)*R
          X3 = XI*R*R*X1
          X5 = XI*R*R*X3
          X7 = XI*R*R*X5
          X9 = XI*R*R*X7
          ERNUC = 4.0D0*X1/3.0D0 - 4.0D0*X3/5.0D0 + 2.0D0*X5/7.0D0
     &          - 2.0D0*X7/27.0D0 + X9/66.0D0
          ERNUC = XI*ERNUC/PI12
        ELSE
          X1 = DSQRT(XI)*R
          ERNUC =-2.0D0*X1*DEXP(-X1*X1)/PI12 + DERF(X1)
          ERNUC = ERNUC/(R*R)
        ENDIF
        RETURN
C
      ENDIF
C
C     GAUSSIAN NUCLEAR CHARGE MODEL
      IF(MODEL.EQ.'GAUSS') THEN
C
        XI = 3.0D0/(2.0D0*RNUC(IZ)*RNUC(IZ))
        IF(R.LT.0.2D0*RNUC(IZ)) THEN
          X1 = DSQRT(XI)*R
          X3 = XI*R*R*X1
          X5 = XI*R*R*X3
          X7 = XI*R*R*X5
          X9 = XI*R*R*X7
          ERNUC = 4.0D0*X1/3.0D0 - 4.0D0*X3/5.0D0 + 2.0D0*X5/7.0D0
     &          - 2.0D0*X7/27.0D0 + X9/66.0D0
          ERNUC = XI*ERNUC/PI12
        ELSE
          X1 = DSQRT(XI)*R
          ERNUC =-2.0D0*X1*DEXP(-X1*X1)/PI12 + DERF(X1)
          ERNUC = ERNUC/(R*R)
        ENDIF
C      
      ENDIF
C
C     FERMI NUCLEAR CHARGE MODEL
      IF(MODEL.EQ.'FERMI') THEN
C
        A = AFMI(IZ)
        C = 5.0D0*RNUC(IZ)*RNUC(IZ)/3.0D0
     &    - 7.0D0*PI*PI*AFMI(IZ)*AFMI(IZ)/3.0D0
        C = DSQRT(C)
        U = A/C
C
C       NORMALISATION CONSTANT
        Y1 = 1.0D0
        Y2 = PI*PI*U*U
        Y3 =-6.0D0*U*U*U*POLYLOG(3,-1.0D0/U)
        YN = 1.0D0/(Y1+Y2+Y3)
C
C       ELECTRIC FIELD BY EXPANSION
        IF(R.LT.0.2D0*C) THEN
          ERNUC =-YN*R/(C*C*C)
        ELSEIF(R.LT.C) THEN
          E1 =-R*R*R/(C*C*C)
          E3 = 6.0D0*U*U*U*POLYLOG(3,-C/A)
          E4 =-6.0D0*U*U*U*POLYLOG(3,(R-C)/A)
          E5 =-3.0D0*U*R*R*POLYLOG(1,(R-C)/A)/(C*C)
          E6 = 6.0D0*U*U*R*POLYLOG(2,(R-C)/A)/C
          ERNUC = YN*(E1+   E3+E4+E5+E6)/(R*R)
        ELSE
          E1 =-1.0D0
          E2 =-(PI*U)**2
          E3 = 6.0D0*U*U*U*POLYLOG(3,-C/A)
          E4 =-6.0D0*U*U*U*POLYLOG(3,(C-R)/A)
          E5 =-3.0D0*U*R*R*POLYLOG(1,(C-R)/A)/(C*C)
          E6 =-6.0D0*U*U*R*POLYLOG(2,(C-R)/A)/C
          ERNUC = YN*(E1+E2+E3+E4+E5+E6)/(R*R)
        ENDIF
C
      ENDIF
C
C     HOMOGENEOUS BALL MODEL
      IF(MODEL.EQ.'UNIFM') THEN
C
        B = DSQRT(5.0D0/3.0D0)*RNUC(IZ)
        IF(R.LE.B) THEN
          RN3   = RNUC(IZ)*RNUC(IZ)*RNUC(IZ)
          ERNUC = R/(B*B*B)
        ELSE
          ERNUC = 1.0D0/(R*R)
        ENDIF
C
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION ERNUCDNS(MODEL,IZ,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    EEEEEEEE RRRRRRR  NN    NN  CCCCCC  DDDDDDD  NN    NN  SSSSSS     C
C    EE       RR    RR NNN   NN CC    CC DD    DD NNN   NN SS    SS    C
C    EE       RR    RR NNNN  NN CC       DD    DD NNNN  NN SS          C
C    EEEEEE   RR    RR NN NN NN CC       DD    DD NN NN NN  SSSSSS     C
C    EE       RRRRRRR  NN  NNNN CC       DD    DD NN  NNNN       SS    C
C    EE       RR    RR NN   NNN CC    CC DD    DD NN   NNN SS    SS    C
C    EEEEEEEE RR    RR NN    NN  CCCCCC  DDDDDDD  NN    NN  SSSSSS     C
C                                                                      C
C -------------------------------------------------------------------- C
C  ERNUCDNS EVALUATES THE VALUE OF THE NORMALISED NUCLEAR E_R FIELD    C
C  AT RADIUS R, MULTIPLIED BY R**2, BASED ON A NUCLEAR CHARGE MODEL.   C
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
        WRITE(6, *) 'In ERNUCDNS: radius R must be non-negative. R = ',R
        WRITE(7, *) 'In ERNUCDNS: radius R must be non-negative. R = ',R
        STOP
      ENDIF
C
C     POINT-LINE NUCLEAR CHARGE MODEL
      IF(MODEL.EQ.'POINT') THEN
        ERNUCDNS = 1.0D0
      ENDIF
C
C     GAUSSIAN NUCLEAR CHARGE MODEL
      IF(MODEL.EQ.'GAUSS') THEN
        XI = 3.0D0/(2.0D0*RNUC(IZ)*RNUC(IZ))
        X1 = DSQRT(XI)*R
        ERNUCDNS =-2.0D0*X1*DEXP(-X1*X1)/PI12 + DERF(X1)
      ENDIF
C
C     FERMI NUCLEAR CHARGE MODEL
      IF(MODEL.EQ.'FERMI') THEN
C
        A = AFMI(IZ)
        C = 5.0D0*RNUC(IZ)*RNUC(IZ)/3.0D0
     &    - 7.0D0*PI*PI*AFMI(IZ)*AFMI(IZ)/3.0D0
        C = DSQRT(C)
        U = A/C
C
C       NORMALISATION CONSTANT
        Y1 = 1.0D0
        Y2 = PI*PI*U*U
        Y3 =-6.0D0*U*U*U*POLYLOG(3,-1.0D0/U)
        YN = 1.0D0/(Y1+Y2+Y3)
C
C       ELECTRIC FIELD BY EXPANSION
        IF(R.LT.C) THEN
          E1 =-R*R*R/(C*C*C)
          E3 = 6.0D0*U*U*U*POLYLOG(3,-C/A)
          E4 =-6.0D0*U*U*U*POLYLOG(3,(R-C)/A)
          E5 =-3.0D0*U*R*R*POLYLOG(1,(R-C)/A)/(C*C)
          E6 = 6.0D0*U*U*R*POLYLOG(2,(R-C)/A)/C
          ERNUCDNS = YN*(E1+   E3+E4+E5+E6)
        ELSE
          E1 =-1.0D0
          E2 =-(PI*U)**2
          E3 = 6.0D0*U*U*U*POLYLOG(3,-C/A)
          E4 =-6.0D0*U*U*U*POLYLOG(3,(C-R)/A)
          E5 =-3.0D0*U*R*R*POLYLOG(1,(C-R)/A)/(C*C)
          E6 =-6.0D0*U*U*R*POLYLOG(2,(C-R)/A)/C
          ERNUCDNS = YN*(E1+E2+E3+E4+E5+E6)
        ENDIF
C
      ENDIF
C
C     HOMOGENEOUS BALL MODEL
      IF(MODEL.EQ.'UNIFM') THEN
        B = DSQRT(5.0D0/3.0D0)*RNUC(IZ)
        IF(R.LE.B) THEN
          RN3   = RNUC(IZ)*RNUC(IZ)*RNUC(IZ)
          ERNUCDNS = R*R*R/(B*B*B)
        ELSE
          ERNUCDNS = 1.0D0
        ENDIF
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
            IF(KQN.GT.0) THEN
              LQN = KQN
            ELSE
              LQN =-KQN-1
            ENDIF
            NBAS  = NFNC(LQN,IZ)
            MQMAX = 2*IABS(KQN)-1
            IF(MQMAX.GE.MQN) THEN
              LRGE(IZ,KA,MQN) = ILST
              ILST = ILST+NBAS
            ENDIF
          ENDDO
C
C         LABEL MQN>0 BLOCKS
          DO KA=1,NKAP(IZ)
            KQN = KAPA(KA,IZ)
            IF(KQN.GT.0) THEN
              LQN = KQN
            ELSE
              LQN =-KQN-1
            ENDIF
            NBAS  = NFNC(LQN,IZ)
            MQMAX = 2*IABS(KQN)-1
            IF(MQMAX.GE.MQN) THEN
              LRGE(IZ,KA,MQN+1) = ILST
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
        IFULL = 4*LQN+2
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
      SUBROUTINE ORIBHV0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       OOOOOO  RRRRRRR  IIII BBBBBBB  HH    HH VV    VV  000000       C
C      OO    OO RR    RR  II  BB    BB HH    HH VV    VV 00   000      C
C      OO    OO RR    RR  II  BB    BB HH    HH VV    VV 00  0000      C
C      OO    OO RR    RR  II  BBBBBBB  HHHHHHHH VV    VV 00 00 00      C
C      OO    OO RRRRRRR   II  BB    BB HH    HH  VV  VV  0000  00      C
C      OO    OO RR    RR  II  BB    BB HH    HH   VVVV   000   00      C
C       OOOOOO  RR    RR IIII BBBBBBB  HH    HH    VV     000000       C
C                                                                      C
C -------------------------------------------------------------------- C
C  ORIBHV ASSESSES THE ORIGIN BEHAVIOUR OF ELECTRONIC ORBITALS.        C
C  THIS SHOULD NOT BE USED FOR THE POINT-NUCLEAR MODEL.                C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*1 LLAB,CHS
      CHARACTER*2 ELMT(120),ELNM
      CHARACTER*5 NMDL
      CHARACTER*6 CNFG
C
      COMPLEX*16 COEF(MDM,MDM)
C
      DIMENSION QA(MKP),NUMOCC(0:MEL),NORB(0:MEL,MKP+1)
      DIMENSION RGLB(MCT,MKP+1,MKP,10),VE(MCT)
      DIMENSION PLTN(MCT,MKP+1,MKP),IGLB(MCT,MKP+1,MKP)
      DIMENSION R0(MBD),R2(MBD),R4(MBD),EXL(MBS)
      DIMENSION RRB(MBD,MBD)
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
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/MDLV/ELMT
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C**********************************************************************C
C     SORTING SECTION: FIND ADDRESSES FOR USUAL ATOMIC ORDER           C
C**********************************************************************C
C
C     COUNTER FOR OCCUPIED ORBITALS
      IOCCML = 0
C
C     LOOP OVER NUCLEAR CENTRES IN THE MOLECULE
      DO IZ=1,NCNT
C
C       IMPORT ATOMIC CHARGE DETAILS
        IZNV = INT(ZNUC(IZ))
        ELNM = ELMT(IZNV)
        ICRG = IQNC(IZ)
        MLQN = (NKAP(IZ)-1)/2
C
C       INITIALISE ELECTRONIC POTENTIAL AT THE ORIGIN
        VE(IZ) = 0.0D0
C
C       INITIALISE THE ORBITAL FILLING ARRAY
        DO LQN=0,MEL
          DO NQN=1,MKP+1
            NORB(LQN,NQN) = 0
          ENDDO
        ENDDO
C
C       READ OR GENERATE AUFBAU OCCUPANCY FOR THIS ATOM
        IF(CNFG(IZ).EQ.'MANUAL') THEN
          LMXCONF = (NKAP(IZ)-1)/2
          DO LQN=0,LMXCONF
            NUMOCC(LQN) = NLVL(IZ,LQN)
            DO NQN=1,NLVL(IZ,LQN)
              NORB(LQN,NQN) = NCNF(IZ,LQN,NQN)
            ENDDO
          ENDDO
        ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
          CALL AUFBAU(IZNV,ICRG,NORB,NUMOCC,LMXCONF)
        ELSE
          WRITE(6, *) 'In ORIBHV0: invalid configuration choice.'
          WRITE(7, *) 'In ORIBHV0: invalid configuration choice.'
        ENDIF
C
C       LOOP OVER LQNS
        DO LQN=0,LMXCONF
C
C         RECORD LQNA VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
          NBAS = NFNC(LQN,IZ)
          DO IBAS=1,NBAS
            EXL(IBAS) = BEXL(IBAS,LQN,IZ)
          ENDDO
C
C         FRACTIONAL SUBSHELL OCCUPANCY
          DO NQN=1,NUMOCC(LQN)
            IF(NOCC.EQ.1) THEN
              QA(NQN) = 1.0D0
            ELSE
              QA(NQN) = DFLOAT(NORB(LQN,NQN))/DFLOAT(4*LQN+2)
            ENDIF
          ENDDO
C
C         POSITIVE KQN CHOICE (APPLIES ONLY FOR LQN>0)
C
          IF(LQN.EQ.0.OR.HMLT.EQ.'NORL') GOTO 110
C
C         CALCULATE BATCH OF ORIGIN AMPLITUDES
          CALL AMPLT0(R0,R2,R4,EXL, LQN  ,NBAS)
C
C         CALCULATE BATCH OF <1/R> INTEGRALS
          CALL MOMNT0(RRB,EXL, LQN  ,NBAS,-1)
C
C         LOOP OVER |MQN| CHOICES FOR THIS KQN
          DO MQN=1,LQN
C
C           LOOP OVER LISTED NQNS
            DO NQN=1,NUMOCC(LQN)
C
C             SKIP THIS STEP IF ORBITAL IS UNOCCUPIED
              IF(NORB(LQN,NQN).EQ.0) GOTO 101
              IF(IOCTP(IZ,2*LQN  ,1,NQN).EQ.0) GOTO 101
C
C             SPHERICAL SYMMETRY: FOCUS ON FIRST MQN ONLY
              IF(MQN.EQ.1) THEN
C
C               ADDRESS
                IGLB(IZ,NQN+LQN,2*LQN  ) = IOCPN(IOCCML+1)
                PLTN(IZ,NQN+LQN,2*LQN  ) = QA(NQN)
C
C               FOCK MATRIX ADDRESSES
                IL = LRGE(IZ,2*LQN  ,MQN)
                IS = LRGE(IZ,2*LQN  ,MQN)+NSKP
C
C               COEFFICIENT MATRIX OFFSET
                MOCC = IOCPN(IOCCML+1)+NSKP
C
C               POLYNOMIAL COEFFICIENTS
                P0 = 0.0D0
                P2 = 0.0D0
                P4 = 0.0D0
                Q0 = 0.0D0
                Q2 = 0.0D0
                Q4 = 0.0D0
                DO IBAS=1,NBAS
C
                  P0 = P0 + DREAL(COEF(IL+IBAS,MOCC))*R0(IBAS     )
                  P2 = P2 + DREAL(COEF(IL+IBAS,MOCC))*R2(IBAS     )
                  P4 = P4 + DREAL(COEF(IL+IBAS,MOCC))*R4(IBAS     )
C
                  Q0 = Q0 + DREAL(COEF(IS+IBAS,MOCC))*R0(IBAS+NBAS)
                  Q2 = Q2 + DREAL(COEF(IS+IBAS,MOCC))*R2(IBAS+NBAS)
                  Q4 = Q4 + DREAL(COEF(IS+IBAS,MOCC))*R4(IBAS+NBAS)
C
                ENDDO
                RGLB(IZ,NQN+LQN,2*LQN  ,1) = P0/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN  ,2) = P2/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN  ,3) = P4/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN  ,4) = Q0/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN  ,5) = Q2/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN  ,6) = Q4/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN  ,7) = P0/Q0
                RGLB(IZ,NQN+LQN,2*LQN  ,8) = P2/Q2
                RGLB(IZ,NQN+LQN,2*LQN  ,9) = P4/Q4
C
C               RECIPROCAL MOMENT
                RB = 0.0D0
                DO IBAS=1,NBAS
                  DO JBAS=1,NBAS
C
                    DAA = DREAL(COEF(IL+IBAS,MOCC)*COEF(IL+JBAS,MOCC))
                    DBB = DREAL(COEF(IS+IBAS,MOCC)*COEF(IS+JBAS,MOCC))
C
                    RB = RB + DAA*RRB(IBAS     ,JBAS     ) 
     &                      + DBB*RRB(IBAS+NBAS,JBAS+NBAS)
C
                  ENDDO
                ENDDO
                RGLB(IZ,NQN+LQN,2*LQN  ,10) = RB/QA(NQN)
C
C               DEGENERACY
                IDGN = 2*LQN
C
C               UPDATE ELECTRONIC POTENTIAL AT THE ORIGIN
                VE(IZ) = VE(IZ) + DFLOAT(IDGN)*RB/QA(NQN)
C
              ENDIF
C
C             UPDATE NUMBER OF OCCUPIED PAIRS
              IOCCML = IOCCML+2
C
C             SKIP STEP FOR UNOCCUPIED ORBITALS
101           CONTINUE
C
C           END LOOP OVER LISTED NQNS
            ENDDO
C
          ENDDO
C
110       CONTINUE
C
C         NEGATIVE KQN CHOICE
C
C         CALCULATE BATCH OF ORIGIN AMPLITUDES
          CALL AMPLT0(R0,R2,R4,EXL,-LQN-1,NBAS)
C
C         CALCULATE BATCH OF <1/R> INTEGRALS
          CALL MOMNT0(RRB,EXL,-LQN-1,NBAS,-1)
C
C         LOOP OVER |MQN| CHOICES
          DO MQN=1,LQN+1
C
C           LOOP OVER LISTED NQNS
            DO NQN=1,NUMOCC(LQN)
C
C             SKIP THIS STEP IF ORBITAL IS UNOCCUPIED
              IF(NORB(LQN,NQN).EQ.0) GOTO 102
              IF(IOCTP(IZ,2*LQN+1,1,NQN).EQ.0) GOTO 102
C
C             SPHERICAL SYMMETRY: FOCUS ON FIRST MQN ONLY
              IF(MQN.EQ.1) THEN
C
C               ADDRESS
                IGLB(IZ,NQN+LQN,2*LQN+1) = IOCPN(IOCCML+1)
                PLTN(IZ,NQN+LQN,2*LQN+1) = QA(NQN)
C
C               FOCK MATRIX ADDRESSES
                IL = LRGE(IZ,2*LQN+1,MQN)
                IS = LRGE(IZ,2*LQN+1,MQN)+NSKP
C
C               COEFFICIENT MATRIX OFFSET
                MOCC = IOCPN(IOCCML+1)+NSKP
C
C               RADIAL MOMENTS
                P0 = 0.0D0
                P2 = 0.0D0
                P4 = 0.0D0
                Q0 = 0.0D0
                Q2 = 0.0D0
                Q4 = 0.0D0
                DO IBAS=1,NBAS
C
                  P0 = P0 + DREAL(COEF(IL+IBAS,MOCC))*R0(IBAS     )
                  P2 = P2 + DREAL(COEF(IL+IBAS,MOCC))*R2(IBAS     )
                  P4 = P4 + DREAL(COEF(IL+IBAS,MOCC))*R4(IBAS     )
C
                  Q0 = Q0 + DREAL(COEF(IS+IBAS,MOCC))*R0(IBAS+NBAS)
                  Q2 = Q2 + DREAL(COEF(IS+IBAS,MOCC))*R2(IBAS+NBAS)
                  Q4 = Q4 + DREAL(COEF(IS+IBAS,MOCC))*R4(IBAS+NBAS)
C
                ENDDO
                RGLB(IZ,NQN+LQN,2*LQN+1,1) = P0/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN+1,2) = P2/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN+1,3) = P4/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN+1,4) = Q0/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN+1,5) = Q2/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN+1,6) = Q4/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN+1,7) = Q0/P0
                RGLB(IZ,NQN+LQN,2*LQN+1,8) = Q2/P2
                RGLB(IZ,NQN+LQN,2*LQN+1,9) = Q4/P4
C
C               RECIPROCAL MOMENT
                RB = 0.0D0
                DO IBAS=1,NBAS
                  DO JBAS=1,NBAS
C
                    DAA = DREAL(COEF(IL+IBAS,MOCC)*COEF(IL+JBAS,MOCC))
                    RB = RB + DAA*RRB(IBAS     ,JBAS     ) 
C
                    IF(HMLT.EQ.'NORL') GOTO 115
C
                    DBB = DREAL(COEF(IS+IBAS,MOCC)*COEF(IS+JBAS,MOCC))
                    RB = RB + DBB*RRB(IBAS+NBAS,JBAS+NBAS)
C
115                 CONTINUE
C
                  ENDDO
                ENDDO
                RGLB(IZ,NQN+LQN,2*LQN+1,10) = RB/QA(NQN)
C
C               DEGENERACY
                IF(HMLT.EQ.'NORL') THEN
                  IDGN = 4*LQN+2
                ELSE
                  IDGN = 2*LQN+2
                ENDIF
C
C               UPDATE ELECTRONIC POTENTIAL AT THE ORIGIN
                VE(IZ) = VE(IZ) + DFLOAT(IDGN)*RB/QA(NQN)
C
              ENDIF
C
C             UPDATE NUMBER OF OCCUPIED PAIRS
              IOCCML = IOCCML+2
C
C             SKIP STEP FOR UNOCCUPIED ORBITALS
102           CONTINUE
C
C           END LOOP OVER LISTED NQNS
            ENDDO
C
          ENDDO
C
C         NON-RELATIVISTIC SPECIAL CASE: ALSO FILL IN THE +KQN BLOCK
          IF(HMLT.EQ.'NORL'.AND.LQN.GE.1) THEN
C
C           LOOP OVER |MQN| CHOICES
            DO MQN=1,LQN
C
C             LOOP OVER LISTED NQNS
              DO NQN=1,NUMOCC(LQN)
C
C               SKIP THIS STEP IF ORBITAL IS UNOCCUPIED
                IF(NORB(LQN,NQN).GT.0) THEN
                  IOCCML = IOCCML+2
                ENDIF
C
              ENDDO
C
            ENDDO
C
          ENDIF
C
C       END LOOP OVER LQNS
        ENDDO
C
C     END LOOP OVER NUCLEAR CENTRES
      ENDDO     
C
C**********************************************************************C
C     PRINT SPECTRUM                                                   C
C**********************************************************************C
C
C     PRINT TITLE TO TERMINAL/FILE
20    FORMAT(1X,A,': nuclear model = ',A,3X,A,I3,9X,A,F15.6,A)
21    FORMAT(29X,A,F8.4,4X,A,F15.6,A)
22    FORMAT(1X,' nl',11X,'E (au)',3X,'f(E,V_0,l)',
     &                                 13X,'p0',13X,'q0',6X,'f(p0,q0)')
23    FORMAT(1X,I2,A,A,' | ',F13.6,1X,F12.6,' |',
     &                                       ES13.6,2X,ES13.6,1X,F13.6)
24    FORMAT(1X,A,15X,F13.6,A)
C
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT(' ',19),'Quality of wavefunctions at origin'
      WRITE(7, *) REPEAT(' ',19),'Quality of wavefunctions at origin'
      WRITE(6, *) REPEAT(' ',22),'V(r)/Z = V_0 + ½V_0" r^2 + ... '
      WRITE(7, *) REPEAT(' ',22),'V(r)/Z = V_0 + ½V_0" r^2 + ... '
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C     LOOP OVER NUCLEAR CENTRES IN THE MOLECULE
      DO IZ=1,NCNT
C
C       POTENTIAL POWER SERIES
        X12 = DSQRT(XNUC(IZ,0))
        V0  =-2.0D0*FNUC(IZ,0)*X12/PI12
        V2  = 4.0D0*FNUC(IZ,0)*XNUC(IZ,0)*X12/(3.0D0*PI12)
        IF(NNUC(IZ).GT.0) THEN
          DO IFT=1,NNUC(IZ)
            V0  = V0 + FNUC(IZ,IFT)
            V2  = V2 - 6.0D0*XNUC(IZ,IFT)
          ENDDO
        ENDIF
C       CORRECT ANSWERS SHOULD BE
C       V2 = 4.0D0*PI*RHONUC(NMDL(IZ),IZ,0.0D0)
C
C       IN CASE OF HOMOGENEOUS DISTRIBUTION, OVERRIDE V0 WITH EXACT VALUE
        IF(NMDL(IZ).EQ.'UNIFM') THEN
          V0 =-1.5D0/RNUC(IZ) 
          V2 = 4.0D0*PI*RHONUC('UNIFM',IZ,0.0D0)
        ENDIF
C
C       IN CASE OF FERMI DISTRIBUTION, OVERRIDE V0 WITH EXACT VALUE
        IF(NMDL(IZ).EQ.'FERMI') THEN
          CFMI = 5.0D0*RNUC(IZ)*RNUC(IZ)/3.0D0
     &             - 7.0D0*PI*PI*AFMI(IZ)*AFMI(IZ)/3.0D0
          CFMI = DSQRT(CFMI)
          V0 = VNFERMI(0.0D0,CFMI,TFMI(IZ))
          V2 = 4.0D0*PI*RHONUC('FERMI',IZ,0.0D0)
        ENDIF
C
C       GIVEN DHFQ HAMILTONIAN, INCLUDE SELF-INTERACTION TERM
c       dfnote: wrong in existing form -- figure out
        IF(HMLT.EQ.'DHFQ') THEN
C
C         LOW- AND HIGH-ENERGY PARTITION CUTOFF
          CUTK = CV
C
C         AMPLITUDE FOR FREE-WAVE SELF-INTERACTION TERM (B&S 19.3)
          AG2 = (DLOG(EMSS*CV/CUTK)-TWLG+11.0D0/24.0D0)/(3.0D0*CV*PI)
          AG2 =-AG2/(EMSS*EMSS*CV*CV)
C
C         CONTRIBUTION TO V0
c         V0 = V0 - AG2*4.0D0*PI*RHONUC(NMDL(IZ),IZ,0.0D0)
C
        ENDIF
C
C       SCALE TO APROPRIATE UNITS (FM)
        V2 = V2/(CFM*CFM)
C
C       PRINT HEADER
        WRITE(6,20) ELMT(INT(ZNUC(IZ))),NMDL(IZ),
     &              'Z = ',INT(ZNUC(IZ)),'V_0 =',ZNUC(IZ)*V0,' au'
        WRITE(7,20) ELMT(INT(ZNUC(IZ))),NMDL(IZ),
     &              'Z = ',INT(ZNUC(IZ)),'V_0 =',ZNUC(IZ)*V0,' au'
        WRITE(6,21) 'A = ',    ANUC(IZ) ,'V_0"=',ZNUC(IZ)*V2,' au/fm^2'
        WRITE(7,21) 'A = ',    ANUC(IZ) ,'V_0"=',ZNUC(IZ)*V2,' au/fm^2'
        WRITE(6, *)
        WRITE(7, *)
        WRITE(6,24) 'Electronic potential at the origin =',VE(IZ),' au'
        WRITE(7,24) 'Electronic potential at the origin =',VE(IZ),' au'
        WRITE(6, *) REPEAT('=',77)
        WRITE(7, *) REPEAT('=',77)
        WRITE(6,22)
        WRITE(7,22)
        WRITE(6, *) REPEAT('-',77)
        WRITE(7, *) REPEAT('-',77)
C
C       IMPORT ATOMIC CHARGE DETAILS
        IZNV = INT(ZNUC(IZ))
        ICRG = IQNC(IZ)
C
C       INITIALISE THE ORBITAL FILLING ARRAY
        DO LQN=0,MEL
          DO NQN=1,MKP+1
            NORB(LQN,NQN) = 0
          ENDDO
        ENDDO
C
C       READ OR GENERATE AUFBAU OCCUPANCY FOR THIS ATOM
        IF(CNFG(IZ).EQ.'MANUAL') THEN
          LMXCONF = (NKAP(IZ)-1)/2
          DO L=0,LMXCONF
            NUMOCC(L) = NLVL(IZ,L)
            DO N=1,NLVL(IZ,L)
              NORB(L,N) = NCNF(IZ,L,N)
            ENDDO
          ENDDO
        ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
          CALL AUFBAU(IZNV,ICRG,NORB,NUMOCC,LMXCONF)
        ELSE
          WRITE(6, *) 'In ORIBHV0: invalid configuration choice.'
          WRITE(7, *) 'In ORIBHV0: invalid configuration choice.'
        ENDIF
C
C       LOOP OVER ORBITALS IN THE USUAL AUFBAU ORDER
        DO NQN=1,MKP
          DO LQN=0,NQN-1
C
            IF(LQN.GT.MEL) GOTO 150
            IF(NORB(LQN,NQN-LQN).EQ.0) GOTO 150
C
C           POSITIVE KQN CHOICE (APPLIES ONLY FOR LQN>0)
            IF(LQN.EQ.0.OR.HMLT.EQ.'NORL') GOTO 120
            IF(IOCTP(IZ,2*LQN  ,1,NQN-LQN).EQ.0) GOTO 120
C
C           SYMMETRY LABEL
            CHS = '+'
C
C           IDENTIFY THE RIGHT INDEX AND PROPERTIES
            MOCC = IGLB(IZ,NQN,2*LQN  )+NSKP
            FRAC = PLTN(IZ,NQN,2*LQN  )
            P0   = RGLB(IZ,NQN,2*LQN  , 1)
            P2   = RGLB(IZ,NQN,2*LQN  , 2)
            P4   = RGLB(IZ,NQN,2*LQN  , 3)
            Q0   = RGLB(IZ,NQN,2*LQN  , 4)
            Q2   = RGLB(IZ,NQN,2*LQN  , 5)
            Q4   = RGLB(IZ,NQN,2*LQN  , 6)
            R00  = RGLB(IZ,NQN,2*LQN  , 7)
            R22  = RGLB(IZ,NQN,2*LQN  , 8)
            R44  = RGLB(IZ,NQN,2*LQN  , 9)
            RAVB = RGLB(IZ,NQN,2*LQN  ,10)
C
C           DEGENERACY
            IDGN = 2*LQN
C
C           EFFECTIVE TOTAL POTENTIAL AT THE ORIGIN FOR THIS ELECTRON
C           VT = ZNUC(IZ)*V0 + VE(IZ) - RAVB
            VT = ZNUC(IZ)*V0
C
C           EXPECTED RATIO IN THIS ASYMPTOTIC LIMIT
            VRT =-(VT-2.0D0*CV*CV-EIGN(MOCC))/((2*LQN+1)*CV)
C
C           OUTPUT TO TERMINAL
            WRITE(6,23) NQN,LLAB(LQN),CHS,EIGN(MOCC),VRT,P0,Q0,R00
            WRITE(7,23) NQN,LLAB(LQN),CHS,EIGN(MOCC),VRT,P0,Q0,R00
C
120         CONTINUE
C
C           NEGATIVE KQN CHOICE
            IF(IOCTP(IZ,2*LQN+1,1,NQN-LQN).EQ.0) GOTO 130
C
C           SYMMETRY LABEL IF NEEDED
            IF(LQN.EQ.0.OR.HMLT.EQ.'NORL') THEN
              CHS = ' '
            ELSE
              CHS = '-'
            ENDIF
C
C           IDENTIFY THE RIGHT INDEX AND PROPERTIES
            MOCC = IGLB(IZ,NQN,2*LQN+1)+NSKP
            FRAC = PLTN(IZ,NQN,2*LQN+1)
            P0   = RGLB(IZ,NQN,2*LQN+1, 1)
            P2   = RGLB(IZ,NQN,2*LQN+1, 2)
            P4   = RGLB(IZ,NQN,2*LQN+1, 3)
            Q0   = RGLB(IZ,NQN,2*LQN+1, 4)
            Q2   = RGLB(IZ,NQN,2*LQN+1, 5)
            Q4   = RGLB(IZ,NQN,2*LQN+1, 6)
            R00  = RGLB(IZ,NQN,2*LQN+1, 7)
            R22  = RGLB(IZ,NQN,2*LQN+1, 8)
            R44  = RGLB(IZ,NQN,2*LQN+1, 9)
            RAVB = RGLB(IZ,NQN,2*LQN+1,10)
C
C           DEGENERACY
            IF(HMLT.EQ.'NORL') THEN
              IDGN = 4*LQN+2
            ELSE
              IDGN = 2*LQN+2
            ENDIF
C
C           EFFECTIVE TOTAL POTENTIAL AT THE ORIGIN FOR THIS ELECTRON
C           VT = ZNUC(IZ)*V0 + VE(IZ)
C           VT = ZNUC(IZ)*V0 + VE(IZ) - RAVB
C           VT = ZNUC(IZ)*V0 + VE(IZ)
            VT = ZNUC(IZ)*V0
C           VT = ZNUC(IZ)*V0 + VE(IZ) - RAVB
C
C           EXPECTED RATIO IN THIS ASYMPTOTIC LIMIT
            IF(HMLT.EQ.'NORL') THEN
              VRT = LQN+3.0D0
            ELSE
              VRT = (VT-EIGN(MOCC))/((2*LQN+3)*CV)
            ENDIF
C
C           OUTPUT TO TERMINAL
            WRITE(6,23) NQN,LLAB(LQN),CHS,EIGN(MOCC),VRT,P0,Q0,R00
            WRITE(7,23) NQN,LLAB(LQN),CHS,EIGN(MOCC),VRT,P0,Q0,R00
C
130         CONTINUE
C
150         CONTINUE
C
          ENDDO
C
        ENDDO
C
C       DELIMETER BETWEEN NUCLEAR CENTRES
        IF(IZ.NE.NCNT) THEN
          WRITE(6, *) REPEAT('=',77)
          WRITE(7, *) REPEAT('=',77)
          WRITE(6, *) ' '
          WRITE(7, *) ' '
        ENDIF

C     END LOOP OVER NUCLEAR CENTRES
      ENDDO
C
      WRITE(6, *) REPEAT('=',77)
      WRITE(7, *) REPEAT('=',77)
C
      RETURN
      END
C
C
      SUBROUTINE SPECTRM0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    SSSSSS  PPPPPPP   CCCCCC TTTTTTTT RRRRRRR  MM       MM  000000    C
C   SS    SS PP    PP CC    CC   TT    RR    RR MMM     MMM 00   000   C
C   SS       PP    PP CC         TT    RR    RR MMMM   MMMM 00  0000   C
C    SSSSSS  PP    PP CC         TT    RR    RR MM MM MM MM 00 00 00   C
C         SS PPPPPPP  CC         TT    RRRRRRR  MM  MMM  MM 0000  00   C
C   SS    SS PP       CC    CC   TT    RR    RR MM   M   MM 000   00   C
C    SSSSSS  PP        CCCCCC    TT    RR    RR MM       MM  000000    C
C                                                                      C
C -------------------------------------------------------------------- C
C  SPECTRM0 SUMMARISES RESULTS OF THE AVERAGE-OVER-CONFIGURATION       C
C  ATOMIC HARTREE-FOCK SOLUTIONS IN THE USUAL DIRAC BASIS, WITH SOME   C
C  RADIAL EXPECTATION VALUES FOR EACH DEGENERATE ENERGY LEVEL.         C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*1 LLAB,CHS
      CHARACTER*2 ELMT(120),ELNM
      CHARACTER*5 NMDL
      CHARACTER*6 CNFG
C
      COMPLEX*16 COEF(MDM,MDM)
C
      DIMENSION QA(MKP),NUMOCC(0:MEL),NORB(0:MEL,MKP+1)
      DIMENSION PLTN(MCT,MKP+1,MKP),IGLB(MCT,MKP+1,MKP)
      DIMENSION RRA(MBD,MBD),RRB(MBD,MBD),RRC(MBD,MBD),
     &          RRD(MBD,MBD),RRE(MBD,MBD)
      
      DIMENSION EXL(MBS)
      DIMENSION RGLB(MCT,MKP+1,MKP,5)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/FILL/NCNF(MCT,0:MEL,MKP+1),NLVL(MCT,0:MEL),CNFG(MCT)
      COMMON/MDLV/ELMT
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
      COMMON/SHLL/ACFF,BCFF,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
C
C**********************************************************************C
C     SORTING SECTION: FIND ADDRESSES FOR USUAL ATOMIC ORDER           C
C**********************************************************************C
C
C     COUNTER FOR OCCUPIED ORBITALS
      IOCCML = 0
C
C     LOOP OVER NUCLEAR CENTRES IN THE MOLECULE
      DO IZ=1,NCNT
C
C       IMPORT ATOMIC CHARGE DETAILS
        IZNV = INT(ZNUC(IZ))
        ELNM = ELMT(IZNV)
        ICRG = IQNC(IZ)
        MLQN = (NKAP(IZ)-1)/2
C
C       INITIALISE THE ORBITAL FILLING ARRAY
        DO LQN=0,MEL
          DO NQN=1,MKP+1
            NORB(LQN,NQN) = 0
          ENDDO
        ENDDO
C
C       READ OR GENERATE AUFBAU OCCUPANCY FOR THIS ATOM
        IF(CNFG(IZ).EQ.'MANUAL') THEN
          LMXCONF = (NKAP(IZ)-1)/2
          DO LQN=0,LMXCONF
            NUMOCC(LQN) = NLVL(IZ,LQN)
            DO NQN=1,NLVL(IZ,LQN)
              NORB(LQN,NQN) = NCNF(IZ,LQN,NQN)
            ENDDO
          ENDDO
        ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
          CALL AUFBAU(IZNV,ICRG,NORB,NUMOCC,LMXCONF)
        ELSE
          WRITE(6, *) 'In SPECTRM0: invalid configuration choice.'
          WRITE(7, *) 'In SPECTRM0: invalid configuration choice.'
        ENDIF
C
C       LOOP OVER LQNS
        DO LQN=0,LMXCONF
C
C         RECORD LQNA VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
          NBAS = NFNC(LQN,IZ)
          DO IBAS=1,NBAS
            EXL(IBAS) = BEXL(IBAS,LQN,IZ)
          ENDDO
C
C         FRACTIONAL SUBSHELL OCCUPANCY
          DO NQN=1,NUMOCC(LQN)
            IF(NOCC.EQ.1) THEN
              QA(NQN) = 1.0D0
            ELSE
              QA(NQN) = DFLOAT(NORB(LQN,NQN))/DFLOAT(4*LQN+2)
            ENDIF
          ENDDO
C
C         POSITIVE KQN CHOICE (APPLIES ONLY FOR LQN>0)
C
          IF(LQN.EQ.0.OR.HMLT.EQ.'NORL') GOTO 110
C
C         CALCULATE BATCH OF <1/R>, <R> AND OTHER INTEGRALS
          CALL MOMNT0(RRA,EXL, LQN  ,NBAS,-2)
          CALL MOMNT0(RRB,EXL, LQN  ,NBAS,-1)
          CALL MOMNT0(RRC,EXL, LQN  ,NBAS, 1)
          CALL MOMNT0(RRD,EXL, LQN  ,NBAS, 2)
          CALL MOMNT0(RRE,EXL, LQN  ,NBAS, 3)
C
C         LOOP OVER |MQN| CHOICES FOR THIS KQN
          DO MQN=1,LQN
C
C           LOOP OVER LISTED NQNS
            DO NQN=1,NUMOCC(LQN)
C
C             SKIP THIS STEP IF ORBITAL IS UNOCCUPIED
              IF(NORB(LQN,NQN).EQ.0) GOTO 101
              IF(IOCTP(IZ,2*LQN  ,1,NQN).EQ.0) GOTO 101
C
C             SPHERICAL SYMMETRY: FOCUS ON FIRST MQN ONLY
              IF(MQN.EQ.1) THEN
C
C               ADDRESS AND AVERAGE OCCUPATION
                IGLB(IZ,NQN+LQN,2*LQN  ) = IOCPN(IOCCML+1)
                PLTN(IZ,NQN+LQN,2*LQN  ) = QA(NQN)
C
C               FOCK MATRIX ADDRESSES
                IL = LRGE(IZ,2*LQN  ,MQN)
                IS = LRGE(IZ,2*LQN  ,MQN)+NSKP
C
C               COEFFICIENT MATRIX OFFSET
                MOCC = IOCPN(IOCCML+1)+NSKP
C
C               RADIAL MOMENTS
                RA = 0.0D0
                RB = 0.0D0
                RC = 0.0D0
                RD = 0.0D0
                RE = 0.0D0
                DO IBAS=1,NBAS
                  DO JBAS=1,NBAS
C
C                   LARGE-COMPONENT OFFSETS
                    IA = IL+IBAS
                    JA = IL+JBAS
C
C                   SMALL-COMPONENT OFFSETS
                    IB = IS+IBAS
                    JB = IS+JBAS
C
                    KBAS = IBAS+NBAS
                    LBAS = JBAS+NBAS
                    
                    DAA = DREAL(DCONJG(COEF(IA,MOCC))*COEF(JA,MOCC))
                    DBB = DREAL(DCONJG(COEF(IB,MOCC))*COEF(JB,MOCC))
C
                    RA = RA + DAA*RRA(IBAS,JBAS) + DBB*RRA(KBAS,LBAS)
                    RB = RB + DAA*RRB(IBAS,JBAS) + DBB*RRB(KBAS,LBAS)
                    RC = RC + DAA*RRC(IBAS,JBAS) + DBB*RRC(KBAS,LBAS)
                    RD = RD + DAA*RRD(IBAS,JBAS) + DBB*RRD(KBAS,LBAS)
                    RE = RE + DAA*RRE(IBAS,JBAS) + DBB*RRE(KBAS,LBAS)
C
                  ENDDO
                ENDDO
                RGLB(IZ,NQN+LQN,2*LQN  ,1) = RA/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN  ,2) = RB/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN  ,3) = RC/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN  ,4) = RD/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN  ,5) = RE/QA(NQN)
C
              ENDIF
C
C             UPDATE NUMBER OF OCCUPIED PAIRS
              IOCCML = IOCCML+2
C
C             SKIP STEP FOR UNOCCUPIED ORBITALS
101           CONTINUE
C
C           END LOOP OVER LISTED NQNS
            ENDDO
C
          ENDDO
C
110       CONTINUE
C
C         NEGATIVE KQN CHOICE
C
C         CALCULATE BATCH OF <1/R>, <R> AND OTHER INTEGRALS
          CALL MOMNT0(RRA,EXL,-LQN-1,NBAS,-2)
          CALL MOMNT0(RRB,EXL,-LQN-1,NBAS,-1)
          CALL MOMNT0(RRC,EXL,-LQN-1,NBAS, 1)
          CALL MOMNT0(RRD,EXL,-LQN-1,NBAS, 2)
          CALL MOMNT0(RRE,EXL,-LQN-1,NBAS, 3)
C
C         LOOP OVER |MQN| CHOICES
          DO MQN=1,LQN+1
C
C           LOOP OVER LISTED NQNS
            DO NQN=1,NUMOCC(LQN)
C
C             SKIP THIS STEP IF ORBITAL IS UNOCCUPIED
              IF(NORB(LQN,NQN).EQ.0) GOTO 102
              IF(IOCTP(IZ,2*LQN+1,1,NQN).EQ.0) GOTO 102
C
C             SPHERICAL SYMMETRY: FOCUS ON FIRST MQN ONLY
              IF(MQN.EQ.1) THEN
C
C               ADDRESS AND AVERAGE OCCUPATION
                IGLB(IZ,NQN+LQN,2*LQN+1) = IOCPN(IOCCML+1)
                PLTN(IZ,NQN+LQN,2*LQN+1) = QA(NQN)
C
C               FOCK MATRIX ADDRESSES
                IL = LRGE(IZ,2*LQN+1,MQN)
                IS = LRGE(IZ,2*LQN+1,MQN)+NSKP
C
C               COEFFICIENT MATRIX OFFSET
                MOCC = IOCPN(IOCCML+1)+NSKP
C
C               RADIAL MOMENTS
                RA = 0.0D0
                RB = 0.0D0
                RC = 0.0D0
                RD = 0.0D0
                RE = 0.0D0
                DO IBAS=1,NBAS
                  DO JBAS=1,NBAS
C
C                   LARGE-COMPONENT OFFSETS
                    IA = IL+IBAS
                    JA = IL+JBAS
C
C                   SMALL-COMPONENT OFFSETS
                    IB = IS+IBAS
                    JB = IS+JBAS
C
                    KBAS = IBAS+NBAS
                    LBAS = JBAS+NBAS
C
                    DAA = DREAL(DCONJG(COEF(IA,MOCC))*COEF(JA,MOCC))
                    RA  = RA + DAA*RRA(IBAS,JBAS)
                    RB  = RB + DAA*RRB(IBAS,JBAS)
                    RC  = RC + DAA*RRC(IBAS,JBAS)
                    RD  = RD + DAA*RRD(IBAS,JBAS)
                    RE  = RE + DAA*RRE(IBAS,JBAS)
C
                    IF(HMLT.EQ.'NORL') GOTO 115
C
                    DBB = DREAL(DCONJG(COEF(IB,MOCC))*COEF(JB,MOCC))
                    RA = RA + DBB*RRA(KBAS,LBAS)
                    RB = RB + DBB*RRB(KBAS,LBAS)
                    RC = RC + DBB*RRC(KBAS,LBAS)
                    RD = RD + DBB*RRD(KBAS,LBAS)
                    RE = RE + DBB*RRE(KBAS,LBAS)
C
115                 CONTINUE
C
                  ENDDO
                ENDDO
                RGLB(IZ,NQN+LQN,2*LQN+1,1) = RA/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN+1,2) = RB/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN+1,3) = RC/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN+1,4) = RD/QA(NQN)
                RGLB(IZ,NQN+LQN,2*LQN+1,5) = RE/QA(NQN)
C
              ENDIF
C
C             UPDATE NUMBER OF OCCUPIED PAIRS
              IOCCML = IOCCML+2
C
C             SKIP STEP FOR UNOCCUPIED ORBITALS
102           CONTINUE
C
C           END LOOP OVER LISTED NQNS
            ENDDO
C
          ENDDO
C
C         NON-RELATIVISTIC SPECIAL CASE: ALSO FILL IN THE +KQN BLOCK
          IF(HMLT.EQ.'NORL'.AND.LQN.GE.1) THEN
C
C           LOOP OVER |MQN| CHOICES
            DO MQN=1,LQN
C
C             LOOP OVER LISTED NQNS
              DO NQN=1,NUMOCC(LQN)
C
C               SKIP THIS STEP IF ORBITAL IS UNOCCUPIED
                IF(NORB(LQN,NQN).GT.0) THEN
                  IOCCML = IOCCML+2
                ENDIF
C
              ENDDO
C
            ENDDO
C
          ENDIF
C
C       END LOOP OVER LQNS
        ENDDO
C
C     END LOOP OVER NUCLEAR CENTRES
      ENDDO     
C
C**********************************************************************C
C     PRINT SPECTRUM                                                   C
C**********************************************************************C
C
C     PRINT TITLE TO TERMINAL/FILE
20    FORMAT(1X,'Z  | nl',9X,'Energy (au)',2X,'Ω',2X,
     & 'Frac |',6X,'<1/r^2>',8X,'<1/r>',10X,'<r>',8X,'<r^2>',8X,'<r^3>')
21    FORMAT(1X,A,1X,'|',I2,A,A,F19.12,1X,I2,1X,F5.3,' |',5(1X,ES12.6))
C
      WRITE(6, *) REPEAT(' ',44),'Atomic summary'
      WRITE(7, *) REPEAT(' ',44),'Atomic summary'
      WRITE(6, *) REPEAT('=',103)
      WRITE(7, *) REPEAT('=',103)
      WRITE(6,20)
      WRITE(7,20)
      WRITE(6, *) REPEAT('-',103)
      WRITE(7, *) REPEAT('-',103)
C
C     LOOP OVER NUCLEAR CENTRES IN THE MOLECULE
      DO IZ=1,NCNT
C
C       IMPORT ATOMIC CHARGE DETAILS
        IZNV = INT(ZNUC(IZ))
        ELNM = ELMT(IZNV)
        ICRG = IQNC(IZ)
C
C       INITIALISE THE ORBITAL FILLING ARRAY
        DO LQN=0,MEL
          DO NQN=1,MKP+1
            NORB(LQN,NQN) = 0
          ENDDO
        ENDDO
C
C       READ OR GENERATE AUFBAU OCCUPANCY FOR THIS ATOM
        IF(CNFG(IZ).EQ.'MANUAL') THEN
          LMXCONF = (NKAP(IZ)-1)/2
          DO L=0,LMXCONF
            NUMOCC(L) = NLVL(IZ,L)
            DO N=1,NLVL(IZ,L)
              NORB(L,N) = NCNF(IZ,L,N)
            ENDDO
          ENDDO
        ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
          CALL AUFBAU(IZNV,ICRG,NORB,NUMOCC,LMXCONF)
        ELSE
          WRITE(6, *) 'In SPECTRM0: invalid configuration choice.'
          WRITE(7, *) 'In SPECTRM0: invalid configuration choice.'
        ENDIF
C
C       LOOP OVER ORBITALS IN THE USUAL AUFBAU ORDER
        DO NQN=1,MKP
          DO LQN=0,NQN-1
C
            IF(LQN.GT.MEL) GOTO 150
            IF(NORB(LQN,NQN-LQN).EQ.0) GOTO 150
C
C           POSITIVE KQN CHOICE (APPLIES ONLY FOR LQN>0)
            IF(LQN.EQ.0.OR.HMLT.EQ.'NORL') GOTO 120
            IF(IOCTP(IZ,2*LQN  ,1,NQN-LQN).EQ.0) GOTO 120
C
C           SYMMETRY LABEL
            CHS = '+'
C
C           DEGENERACY
            IDGN = 2*LQN
C
C           IDENTIFY THE RIGHT INDEX AND PROPERTIES
            MOCC = IGLB(IZ,NQN,2*LQN  )+NSKP
            FRAC = PLTN(IZ,NQN,2*LQN  )
            RAVA = RGLB(IZ,NQN,2*LQN  ,1)
            RAVB = RGLB(IZ,NQN,2*LQN  ,2)
            RAVC = RGLB(IZ,NQN,2*LQN  ,3)
            RAVD = RGLB(IZ,NQN,2*LQN  ,4)
            RAVE = RGLB(IZ,NQN,2*LQN  ,5)
C
C           OUTPUT TO TERMINAL
            WRITE(6,21) ELNM,NQN,LLAB(LQN),CHS,EIGN(MOCC),IDGN,FRAC,
     &                                         RAVA,RAVB,RAVC,RAVD,RAVE
            WRITE(7,21) ELNM,NQN,LLAB(LQN),CHS,EIGN(MOCC),IDGN,FRAC,
     &                                         RAVA,RAVB,RAVC,RAVD,RAVE
C
120         CONTINUE
C
C           NEGATIVE KQN CHOICE
            IF(IOCTP(IZ,2*LQN+1,1,NQN-LQN).EQ.0) GOTO 130
C
C           SYMMETRY LABEL IF NEEDED
            IF(LQN.EQ.0.OR.HMLT.EQ.'NORL') THEN
              CHS = ' '
            ELSE
              CHS = '-'
            ENDIF
C
C           DEGENERACY
            IF(HMLT.EQ.'NORL') THEN
              IDGN = 4*LQN+2
            ELSE
              IDGN = 2*LQN+2
            ENDIF
C
C           IDENTIFY THE RIGHT INDEX AND PROPERTIES
            MOCC = IGLB(IZ,NQN,2*LQN+1)+NSKP
            FRAC = PLTN(IZ,NQN,2*LQN+1)
            RAVA = RGLB(IZ,NQN,2*LQN+1,1)
            RAVB = RGLB(IZ,NQN,2*LQN+1,2)
            RAVC = RGLB(IZ,NQN,2*LQN+1,3)
            RAVD = RGLB(IZ,NQN,2*LQN+1,4)
            RAVE = RGLB(IZ,NQN,2*LQN+1,5)
C
CC          IDENTIFY ADDRESS LOCATION OF OCCUPIED STATE
C           IOCPN(IOCCML+1) = LRGE(IZ,KP,MADD)+NQN
C           IOCTP(IZ,KP,MADD,NQN) = IOCPN(IOCCML+1)
C
C           OUTPUT TO TERMINAL
            WRITE(6,21) ELNM,NQN,LLAB(LQN),CHS,EIGN(MOCC),IDGN,FRAC,
     &                                         RAVA,RAVB,RAVC,RAVD,RAVE
            WRITE(7,21) ELNM,NQN,LLAB(LQN),CHS,EIGN(MOCC),IDGN,FRAC,
     &                                         RAVA,RAVB,RAVC,RAVD,RAVE
C
130         CONTINUE
C
150         CONTINUE
C
          ENDDO
C
        ENDDO
C
C       DELIMETER BETWEEN NUCLEAR CENTRES
        IF(IZ.NE.NCNT) THEN
          WRITE(6, *) REPEAT('-',103)
          WRITE(7, *) REPEAT('-',103)
        ENDIF

C     END LOOP OVER NUCLEAR CENTRES
      ENDDO
C
      WRITE(6, *) REPEAT('=',103)
      WRITE(7, *) REPEAT('=',103)
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
      CALL SYSTEM_CLOCK(ICL1,RATE)
C
C     DIFFERENCE BETWEEN TIME NOW AND START TIME
      TDIF = DFLOAT(ICL1)/RATE-TBEG
C
C     PROJECTED TOTAL COMPLETION TIME
      TTOT = (DFLOAT(ICL1)/RATE-TBEG)/RAT
C
C     PROJECTED END TIME
      TEND = TBEG + TTOT
C
C     EXPECTED TIME LEFT
      TLFT = TTOT-DFLOAT(ICL1)/RATE+TBEG
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
C   [A] EIGSORT: ORGANISE EIGENVALUES BY OCCUPATION AND ORDER.         C
C   [B] DENSTY0: GENERATES CLOSED- AND OPEN-SHELL MOLECULAR DENSITY.   C
C   [C] LEVSHFT: APPLIES A LEVEL SHIFT TO OCCUPIED ORBITALS IN FOCK.   C
C   [D] ENERGIES: USE DENSITY MATRIX TO CALCULATE ENERGY TERMS.        C
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
      COMPLEX*16 CTMP(MDM,MDM),ETMP(MDM)
      COMPLEX*16 COEF(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
C
C     OCCUPIED SPINORS AT THE START OF THE POSITIVE ENERGY SPECTRUM
      DO IOCC=1,IOCCM0
        ETMP(IOCC) = EIGN(IOCPN(IOCC)+NSKP)
        DO I=1,NDIM
          CTMP(I,IOCC) = COEF(I,IOCPN(IOCC)+NSKP)
        ENDDO
      ENDDO
C
C     VIRTUAL SPINORS IN THE REMAINING SPECTRUM
      ITRC = IOCCM0
      DO IVRT=1,NDIM-NSKP
        DO IOCC=1,IOCCM0
          IF(IOCPN(IOCC).EQ.IVRT) GOTO 101
        ENDDO
        ITRC = ITRC+1
        ETMP(ITRC) = EIGN(IVRT+NSKP)
        DO I=1,NDIM
          CTMP(I,ITRC) = COEF(I,IVRT+NSKP)
        ENDDO
101     CONTINUE
      ENDDO
C
C     TRANSFER ALL VALUES BACK INTO THE GLOBAL ARRAYS
      DO ISPN=1,NDIM-NSKP
        EIGN(ISPN+NSKP) = ETMP(ISPN)
        DO I=1,NDIM
          COEF(I,ISPN+NSKP) = CTMP(I,ISPN)
        ENDDO
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
      COMPLEX*16 DENT(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/DENS/DENT
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
          DENT(I,J) = SUM
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C   [4] ATOMIC HARTREE-FOCK: SINGLE-CENTRE SCF CALCULATIONS.           C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] ATOMIC: MAIN ROUTNE FOR ATOMIC SCF CALCULATIONS.               C
C   [B] ATOM1E: HYDROGENIC ATOMIC SOLVER (ONE-BODY ONLY, SO NO SCF.)   C
C   [C] HFSCF0: ATOMIC SCF ROUTINE (AVERAGE OF CONFIGURATION MODEL.)   C
C -------------------------------------------------------------------- C
C   [A] OVRLP0: GENERATES ATOMIC OVERLAP MATRIX FOR A GIVEN KQN.       C
C   [B] ONEEL0: GENERATES ATOMIC ONE-ELECTRON MATRIX (ALL HAMILS).     C
C   [C] KINTC0: GENERATES ATOMIC KINETIC MATRIX FOR A GIVEN KQN.       C
C   [D] PNTNC0: GENERATES ATOMIC POINT-NUCLEAR MATRIX FOR A GIVEN KQN. C
C   [E] UNINC0: GENERATES ATOMIC UNIFORM-NUCLEAR MATRIX   "         ". C
C   [F] GSSNC0: GENERATES ATOMIC GAUSSIAN-NUCLEAR MATRIX  "         ". C
C   [G] NCOLP0: GENERATES ATOMIC NUCLEAR OVERLAP MATRIX (GAUSSIAN).    C
C   [H] MOMNT0: GENERATES ATOMIC RADIAL INTEGRALS FOR A GIVEN KQN.     C
C   [I] UEHNC0: GENERATES ATOMIC NUCLEAR UEHLING MATRIX.               C
C   [J] WKRNC0: GENERATES ATOMIC NUCLEAR WICHMANN-KROLL MATRIX.        C
C   [K] KSBNC0: GENERATES ATOMIC NUCLEAR KÄLLÉN-SABRY MATRIX.          C
C   [L] COULM0: ATOMIC MEAN-FIELD COULOMB MATRIX (BETA INTEGRALS).     C
C   [M] RKCLM0: BATCH OF ATOMIC COULOMB INTERACTION INTEGRALS.         C
C   [N] BREIT0: CONSTRUCTION OF ATOMIC BREIT INTERACTION MATRIX.       C
C   [O] RKBRT0: BATCH OF ATOMIC BREIT INTERACTION INTEGRALS.           C
C -------------------------------------------------------------------- C
C   [A] ERFINT0: INTEGRAL OVER A GAUSSIAN AND ERROR FUNCTION.          C
C   [B] UEHINT0: CALCULATES AN ATOMIC UEHLING INTEGRAL.                C
C   [C] WKRINT0: CALCULATES AN ATOMIC WICHMANN-KROLL INTEGRAL.         C
C   [D] KSBINT0: CALCULATES AN ATOMIC KÄLLÉN-SABRY INTEGRAL.           C
C   [E] ANGCLM0: ATOMIC ANGULAR COULOMB COEFFICIENTS.                  C
C   [F] ANGBRT0: ATOMIC ANGULAR BREIT COEFFICIENTS.                    C
C   [G] BRCOEF0: EXCHANGE MAGNETIC COEFFICIENTS FOR CLOSED-SHELL BREIT.C
C   [H] ANGSQLS: SQUARE OF A WIGNER-3J SYMBOL (LS COUPLING).           C
C   [I] ANGSQJJ: SQUARE OF A WIGNER-3J SYMBOL (JJ COUPLING).           C
C   [J] IJSET0: BASIS SET INTERMEDIATES FOR IJ-PAIRS IN RKCLM0/RKBRT0. C
C   [K] KLSET0: BASIS SET INTERMEDIATES FOR KL-PAIRS IN RKCLM0/RKBRT0. C
C   [L] RNORM0: GENERATE BATCHES OF ALL TT' NORMALISATION FACTORS.     C
C   [M] GAMLWR: LOWER INCOMPLETE GAMMA FUNCTION gamma(A,X).            C
C   [N] GAMUPR: UPPER INCOMPLETE GAMMA FUNCTION GAMMA(A,X).            C
C   [O] VMOMNT: EVEN MOMENTS OF THE FERMI NUCLEAR CHARGE.              C
C   [P] VNFERMI: NORMALISED FERMI POTENTIAL EVALUATED AT RADIUS R.     C
C   [Q] POLYLOG: EVALUATE POLYLOG FUNCTION OF A PARTICULAR ARGUMENT.   C
C**********************************************************************C
C
C
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
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*16 HMS
      CHARACTER*20 STAMP
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 FOCK(MDM,MDM),OVLP(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EWDR,EWXC,EANM,ESLF,EUEE,
     &            EUEU,EUEP,EWKR,EKSB,ESCF
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/MTRX/FOCK,OVLP
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
      COMMON/SHLL/ACFF,BCFF,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
      COMMON/TATM/TTOT
C
      ENUC = 0.0D0
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
      CALL SYSTEM_CLOCK(ICL3,RATE)
      CALL SYSTEM_CLOCK(ICL4)
      TNUC = DFLOAT(ICL4-ICL3)/RATE
      CALL CPU_TIME(TDUM)
C
C     INITIALISE OCCUPATION COUNTER
      IOCCM0 = 0
C
C     IF READING IN A WAVEFUNCTION, BUILD CONFIGURATION AND SKIP
      IF(READIN) THEN
        DO IZ=1,NCNT
          CALL CONFIG0(IZ)
C         CALL H1INT0(IZ,'ANOMLS',EANM)
C         CALL H1INT0(IZ,'SLFLW0',ESLF)
C         CALL H1INT0(IZ,'SLFLWN',ESLF)
C         CALL H1INT0(IZ,'SLFHI0',ESLF)
c         CALL H1INT0(IZ,'UEENUC',EUEE)
c         CALL H1INT0(IZ,'UEUNUC',EUEU)
c         CALL H1INT0(IZ,'UEPNUC',EUEP)
C          CALL H1INT0(IZ,'WKRNUC',EWKR)
C          CALL H1INT0(IZ,'KSBNUC',EKSB)
c          CALL H1INT0(IZ,'UEHSCF',ESCF)
CC         CALL H2INT0(IZ,'UEENUC')
Cc         CALL H1INT0(IZ,'KSBNUC',EKSB)
C          CALL G1PAIR0(IZ,1,'BREIT')
        ENDDO
        GOTO 15
      ENDIF
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
C       BOUNDARY VALUES
        CALL ORIBHV0
        WRITE(6, *) REPEAT(' ',72)
        WRITE(7, *) REPEAT(' ',72)
        IF(HMLT.EQ.'DHFP') THEN
          CALL H1INT0(1,'UEENUC',EUEE)
          CALL H1INT0(1,'UEUNUC',EUEU)
          CALL H1INT0(1,'UEPNUC',EUEP)
          CALL H1INT0(1,'WKRNUC',EWKR)
          CALL H1INT0(1,'KSBNUC',EKSB)
          CALL H1INT0(1,'UEHSCF',ESCF)
C         CALL H1INT0(1,'ANOMLS',EANM)
C         CALL H1INT0(1,'SLFLW0',ESLF)
C         CALL H1INT0(1,'SLFLWN',ESLF)
C         CALL H1INT0(1,'SLFHI0',ESLF)
        ENDIF
        GOTO 300
      ENDIF
C
C     HARTREE-FOCK SCF PROCEDURE FOR EACH ISOLATED ATOM
      DO IZ=1,NCNT
C
C       AVERAGE-OVER-CONFIGURATION ATOMIC SCF PROCEDURE
        CALL HFSCF0(IZ)
C
C       BREIT INTERACTION ORBITAL-BY-ORBITAL
        IF(HMLT.EQ.'DHFB') THEN
          IF(NOCC.NE.1) THEN
            CALL G1PAIR0(IZ,1,'BREIT')
          ENDIF
        ENDIF
C
C       PERTURBATIVE BREIT AND QED INTERACTIONS
        IF(HMLT.EQ.'DHFP') THEN
          EANM = 0.0D0
          ESLF = 0.0D0
          EUEE = 0.0D0
          EUEU = 0.0D0
          EUEP = 0.0D0
          EWKR = 0.0D0
          EKSB = 0.0D0
          CALL H1INT0(IZ,'UEENUC',EUEE)
          CALL H1INT0(IZ,'UEUNUC',EUEU)
          CALL H1INT0(IZ,'UEPNUC',EUEP)
          CALL H1INT0(IZ,'WKRNUC',EWKR)
          CALL H1INT0(IZ,'KSBNUC',EKSB)
          CALL H1INT0(IZ,'UEHSCF',ESCF)
c         CALL H1INT0(IZ,'ANOMLS',EANM)
c         CALL H1INT0(IZ,'SLFLW0',ESLF)
c         CALL H1INT0(IZ,'SLFLWN',ESLF)
c         CALL H1INT0(IZ,'SLFHI0',ESLF)
C         CALL H2INT0(IZ,'UEENUC')
          IF(NOCC.NE.1) THEN
            CALL G1PAIR0(IZ,1,'BREIT')
          ENDIF
        ENDIF
C
C       PERTURBATIVE COULOMB INTERACTION
C       CALL G1PAIR0(IZ,1,'COULM')
C
      ENDDO
C
15    CONTINUE
C
C     GENERATE DENSITY MATRIX
      CALL DENSTY0
C
C     PRODUCE SPECTRUM SUMMARY
      CALL SPECTRM0
C
      IF(READIN) GOTO 16
C
C     MOLECULAR ENERGIES
20    FORMAT(1X,A,29X,F23.14)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
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
      IF(HMLT.NE.'DHFQ'.AND.HMLT.NE.'DHFP') GOTO 202
      WRITE(6,20) 'Anomalous moment (A)',EANM
      WRITE(7,20) 'Anomalous moment (A)',EANM
      WRITE(6,20) 'Self-interaction (S)',ESLF
      WRITE(7,20) 'Self-interaction (S)',ESLF
      WRITE(6,20) 'Uehling e+e-     (E)',EUEE
      WRITE(7,20) 'Uehling e+e-     (E)',EUEE
      WRITE(6,20) 'Wichmann-Kroll   (W)',EWKR
      WRITE(7,20) 'Wichmann-Kroll   (W)',EWKR
      WRITE(6,20) 'Källén-Sabry     (R)',EKSB
      WRITE(7,20) 'Källén-Sabry     (R)',EKSB
      WRITE(6,20) 'Uehling μ+μ-     (U)',EUEU
      WRITE(7,20) 'Uehling μ+μ-     (U)',EUEU
      WRITE(6,20) 'Uehling π+π-     (P)',EUEP
      WRITE(7,20) 'Uehling π+π-     (P)',EUEP
      WRITE(6,20) 'Uehling SCF      (D)',ESCF
      WRITE(7,20) 'Uehling SCF      (D)',ESCF
202   CONTINUE
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,20) 'Molecule total   (F)',ETOT
      WRITE(7,20) 'Molecule total   (F)',ETOT
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
16    CONTINUE
C
C     BOUNDARY VALUES
C      CALL ORIBHV0
C
300   CONTINUE
C
C     SORT THE OCCUPIED STATES (IOCC: 1 -> NOCC)
      CALL EIGSORT
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
      CALL SYSTEM_CLOCK(ICL3,RATE)
      CALL SYSTEM_CLOCK(ICL4)
      TNUC = DFLOAT(ICL4-ICL3)/RATE
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
C
C
      SUBROUTINE ATOM1E(IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             AA   TTTTTTTT OOOOOO  MM       MM  11  EEEEEEEE          C
C            AAAA     TT   OO    OO MMM     MMM 111  EE                C
C           AA  AA    TT   OO    OO MMMM   MMMM  11  EE                C
C          AA    AA   TT   OO    OO MM MM MM MM  11  EEEEEE            C
C          AAAAAAAA   TT   OO    OO MM  MMM  MM  11  EE                C
C          AA    AA   TT   OO    OO MM   M   MM  11  EE                C
C          AA    AA   TT    OOOOOO  MM       MM 1111 EEEEEEEE          C
C                                                                      C
C -------------------------------------------------------------------- C
C  ATOM1E CALCUATES THE SOLUTION TO A ONE-ELECTRON ATOMIC PROBLEM.     C
C  THERE IS NO NEED FOR AVERAGE-OVER CONFIGURATION TREATMENT HERE, SO  C
C  RESULTS ARE STORED IN THE MASTER LIST OF A SINGLE ORBITAL.          C
C  ALL SOLUTIONS OF THE EIGENVALUE PROBLEM ARE EXPORTED TO COEF.       C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ IZ - ATOMIC CENTRE INDEX                                          C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*1 LLAB,QSGN
      CHARACTER*2 ELMT(120),ELNM
      CHARACTER*5 NMDL
      CHARACTER*6 CNFG
      CHARACTER*8 ZWRT,QWRT,EWRT
C
      DIMENSION O(MBD,MBD),H(MBD,MBD),C(MBD,MBD),T(MBD,MBD),V(MBD,MBD),
     &          S(MBD,MBD),Y(MBD,MBD),E(MBD,MBD),U(MBD,MBD),P(MBD,MBD),
     &          W(MBD,MBD),R(MBD,MBD),A(MBD,MBD),D(MBD,MBD)
      DIMENSION EIG(MBD),TE(LWK0),EXL(MBS)
      DIMENSION NORB(0:MEL,MKP+1),NUMOCC(0:MEL)
C
      COMPLEX*16 COEF(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/FILL/NCNF(MCT,0:MEL,MKP+1),NLVL(MCT,0:MEL),CNFG(MCT)
      COMMON/MDLV/ELMT
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     TEMPORARY ELECTRON MASS MODIFIER (DECOUPLES NUCLEAR MOTION)
      XMSS = EMSS
c     EMSS = (XMSS*PMSS*ANUC(IZ))/(XMSS+PMSS*ANUC(IZ))
C
C     IMPORT ATOMIC CHARGE DETAILS
      IZNV = INT(ZNUC(IZ))
      ELNM = ELMT(IZNV)
      ICRG = IQNC(IZ)
      MLQN =(NKAP(IZ)-1)/2
C
C     CONVERT IZNV AND ICRG TO STRINGS AND WRITE TITLE
      IF(IZNV.LT.10) THEN
        WRITE(ZWRT,'(A,I1)') 'Z = ',IZNV
      ELSEIF(IZNV.LT.100) THEN
        WRITE(ZWRT,'(A,I2)') 'Z = ',IZNV
      ELSE
        WRITE(ZWRT,'(A,I3)') 'Z = ',IZNV
      ENDIF
C
      IF(IZNV-ICRG.LT.10) THEN
        WRITE(QWRT,'(A,I2)') 'Q = ',IZNV-ICRG
      ELSEIF(IZNV.LT.100) THEN
        WRITE(QWRT,'(A,I3)') 'Q = ',IZNV-ICRG
      ELSE
        WRITE(QWRT,'(A,I4)') 'Q = ',IZNV-ICRG
      ENDIF
C
      IF(IZNV-ICRG.GT.0) THEN
        QSGN = '+'
      ELSEIF(IZNV-ICRG.LT.0) THEN
        QSGN = '-'
      ENDIF
C
      ICMD = IABS(IZNV-ICRG)
      IF(IZNV-ICRG.EQ.0) THEN
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
20    FORMAT(17X,'Centre',I3,':',3X,A,3X,A,3X,A)
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,20) IZ,ZWRT,QWRT,EWRT
      WRITE(7,20) IZ,ZWRT,QWRT,EWRT
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     READ OR GENERATE AUFBAU OCCUPANCY FOR THIS ATOM
      IF(CNFG(IZ).EQ.'MANUAL') THEN
        LMXCONF = (NKAP(IZ)-1)/2
        DO L=0,LMXCONF
          NUMOCC(L) = NLVL(IZ,L)
          DO N=1,NLVL(IZ,L)
            NORB(L,N) = NCNF(IZ,L,N)
          ENDDO
        ENDDO
      ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
        CALL AUFBAU(IZNV,ICRG,NORB,NUMOCC,LMXCONF)
      ELSE
        WRITE(6, *) 'In ATOM1E: invalid configuration choice.'
        WRITE(7, *) 'In ATOM1E: invalid configuration choice.'
      ENDIF
C
C     CHECK WHETHER THERE ARE SUFFICIENT BASIS FUNCTION TYPES
      IF(MLQN.LT.LMXCONF) THEN
        WRITE(6, *) 'In ATOM1E: insufficient angular types in basis.'
        WRITE(7, *) 'In ATOM1E: insufficient angular types in basis.'
        WRITE(6, *) 'MLQN = ',MLQN,' and LMXCONF = ',LMXCONF
        WRITE(7, *) 'MLQN = ',MLQN,' and LMXCONF = ',LMXCONF
        STOP
      ENDIF
C
C     PRINT ATOMIC CONFIGURATION
31    FORMAT(1X,A,2X,A,1X,'|',2X'NSHELL ',12(2X,I2))
32    FORMAT(1X,'LQN = ',I1,2X,I3,2X,'|'1X,' OCC(',A,'):',A,12(2X,I2))
C
      WRITE(6, *) 'Reduced electron (orbiter) mass = ',EMSS
      WRITE(7, *) 'Reduced electron (orbiter) mass = ',EMSS
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      IF(CNFG(IZ).EQ.'MANUAL') THEN
        WRITE(6,31) 'Manual:','#fns',(N,N=1,12)
        WRITE(7,31) 'Manual:','#fns',(N,N=1,12)
      ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
        WRITE(6,31) 'Aufbau:','#fns',(N,N=1,12)
        WRITE(7,31) 'Aufbau:','#fns',(N,N=1,12)
      ENDIF
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO L=0,(NKAP(IZ)-1)/2
        WRITE(6,32) L,NFNC(L,IZ),LLAB(L),
     &                        REPEAT(' ',L*4),(NORB(L,J),J=1,NUMOCC(L))
        WRITE(7,32) L,NFNC(L,IZ),LLAB(L),
     &                        REPEAT(' ',L*4),(NORB(L,J),J=1,NUMOCC(L))
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
C**********************************************************************C
C     ONE-BODY PROBLEM: ICRG = 1. (NO COULOMB ENERGY.)                 C
C**********************************************************************C
C
C     IMPORT ORDERED ELECTRON OCCUPATION NUMBER
      IOCCML = IOCCM0
C
C     IDENTIFY LQNA AND NSHELL FOR OCCUPIED ELECTRON
      DO L=0,LMXCONF
        DO N=1,NUMOCC(L)
          IF(NORB(L,N).EQ.1) THEN
            LQNOCC = L
            NSHELL = N
          ENDIF
        ENDDO
      ENDDO
C
C     RELATIVISTIC QUANTUM NUMBER
      KQNOCC = LQNOCC
C     KQNOCC =-(LQNOCC+1)
C
C     OVERRIDE FOR LQN=0 WITH POSITIVE SYMMETRY
      IF(LQNOCC.EQ.0) THEN
        KQNOCC =-1
      ENDIF
C
C     MAGNETIC NUMBER (DOUBLE THE ACTUAL VALUE)
      MQNOCC =-1
C
C     CORRESPONDING FOCK ADDRESS FOR THIS MQN
      IF(MQNOCC.LT.0) THEN
        MADD =-MQNOCC
      ELSE
        MADD = MQNOCC+1
      ENDIF
C
C     LOOP OVER ALL AVAILABLE KQN VALUES
      DO KA=1,NKAP(IZ)
C
C       ACTUAL KQN AND LQN VALUES
        KQN = KAPA(KA,IZ)
        LQN = LVAL(KQN)
C
C       IMPORT BASIS FUNCTION EXPONENTS
        NBAS = NFNC(LQN,IZ)
        DO IBAS=1,NBAS
          EXL(IBAS) = BEXL(IBAS,LQN,IZ)
        ENDDO
C
C       MATRIX DIMENSIONS FOR THIS LQN BLOCK
        IF(HMLT.EQ.'NORL') THEN
          NBLC = 0
        ELSE
          NBLC = NBAS
        ENDIF
        NMAT = NBAS+NBLC
C
C       CLEAR HAMILTONIAN AND OVERLAP ARRAYS
        DO IBAS=1,NMAT
          DO JBAS=1,NMAT
            O(IBAS,JBAS) = 0.0D0
            H(IBAS,JBAS) = 0.0D0
            C(IBAS,JBAS) = 0.0D0
            T(IBAS,JBAS) = 0.0D0
            V(IBAS,JBAS) = 0.0D0
            S(IBAS,JBAS) = 0.0D0
            Y(IBAS,JBAS) = 0.0D0
            E(IBAS,JBAS) = 0.0D0
            U(IBAS,JBAS) = 0.0D0
            P(IBAS,JBAS) = 0.0D0
            W(IBAS,JBAS) = 0.0D0
            R(IBAS,JBAS) = 0.0D0
            A(IBAS,JBAS) = 0.0D0
            D(IBAS,JBAS) = 0.0D0
          ENDDO
        ENDDO
C
C       BASIS FUNCTION OVERLAPS
        CALL OVRLP0(O,EXL,KQN,NBAS)
C
C       KINETIC ENERGY
        CALL KINTC0(T,EXL,KQN,NBAS)
C
C       NUCLEAR ATTRACTION
        IF(NMDL(IZ).EQ.'POINT') THEN
           CALL PNTNC0(V,EXL,IZ,KQN,NBAS)
         ELSEIF(NMDL(IZ).EQ.'GAUSS') THEN
           CALL GSSNC0(V,EXL,IZ,KQN,NBAS)
         ELSEIF(NMDL(IZ).EQ.'FERMI') THEN
           CALL FMINC0(V,EXL,IZ,KQN,NBAS)
         ELSEIF(NMDL(IZ).EQ.'UNIFM') THEN
           CALL UNINC0(V,EXL,IZ,KQN,NBAS)
         ENDIF
C
C       FOR NUCLEAR BEST-FIT TESTING PURPOSES
C       CALL GSMNC0(V,EXL,IZ,KQN,NBAS)
C       CALL ONEEL0(T,V,EXL,IZ,KQN,NBAS)
C
C       SELF-INTERACTION (HIGH-ENERGY) AND NUCLEAR VACUUM POLARISATION
        IF(HMLT.EQ.'DHFQ') THEN
c         CALL ANMLS0(A,EXL,IZ,KQN,NBAS)
c         CALL SLFLW0(Y,    IZ,KQN     )
c         CALL SLFHI0(S,EXL,IZ,KQN,NBAS)
          CALL UEHNC0(E,EXL,IZ,KQN,NBAS,'ELEC')
          CALL UEHNC0(U,EXL,IZ,KQN,NBAS,'MUON')
          CALL UEHNC0(P,EXL,IZ,KQN,NBAS,'PION')
          CALL WKRNC0(W,EXL,IZ,KQN,NBAS)
          CALL KSBNC0(R,EXL,IZ,KQN,NBAS)
          CALL UEHEL0(D,EXL,IZ,KQN,NBAS)
        ENDIF
C
C       ONE-ELECTRON HAMILTONIAN MATRIX AND WORKING EIGENVECTOR MATRIX
        DO IBAS=1,NMAT
          DO JBAS=1,NMAT
            H(IBAS,JBAS) = T(IBAS,JBAS)+V(IBAS,JBAS)+A(IBAS,JBAS)
     &                   + S(IBAS,JBAS)+Y(IBAS,JBAS)+E(IBAS,JBAS)
     &                   + U(IBAS,JBAS)+P(IBAS,JBAS)+W(IBAS,JBAS)
     &                   + R(IBAS,JBAS)+D(IBAS,JBAS)
            C(IBAS,JBAS) = H(IBAS,JBAS)
          ENDDO
        ENDDO
C
C       DIAGONALISE MATRIX (THIS NEEDS LAPACK LIBRARY)
        CALL DSYGV(1,'V','U',NMAT,C,MBD,O,MBD,EIG,TE,LWK0,INFO)
        IF(INFO.NE.0) THEN
          WRITE(6, *) 'In ATOM1E: eigenvalue solver DSYGV failed.',INFO
          WRITE(7, *) 'In ATOM1E: eigenvalue solver DSYGV failed.',INFO
          STOP
        ENDIF
C
C       EXPORT ALL SOLUTIONS TO THE GLOBAL COEF ARRAY
C       BEGIN LOOP OVER ALL MQN VALUES
        DO IMVAL=1,IABS(KQN)
C
C         FOCK MATRIX STARTING ADDRESS
          IL1 = LRGE(IZ,KA,2*IMVAL-1)
          IL2 = LRGE(IZ,KA,2*IMVAL  )
          IF(HMLT.NE.'NORL') THEN
            IS1 = LRGE(IZ,KA,2*IMVAL-1)+NSKP
            IS2 = LRGE(IZ,KA,2*IMVAL  )+NSKP
          ENDIF
C
C         LOOP OVER OCCUPIED SOLUTIONS OF THE EIGENVALUE PROBLEM
C         DO IOCC=1,NUMOCC(LQN)
C         LOOP OVER ALL SOLUTIONS OF THE EIGENVALUE PROBLEM
          DO IOCC=1,NBAS
C
C           COPY ENERGY EIGENVALUES TO MASTER LIST
            EIGN(NSKP+IL1+IOCC) = EIG(NBLC+IOCC)
            EIGN(NSKP+IL2+IOCC) = EIG(NBLC+IOCC)
            IF(HMLT.NE.'NORL') THEN
              EIGN(     IL1+IOCC) = EIG(     IOCC)
              EIGN(     IL2+IOCC) = EIG(     IOCC)
            ENDIF
C
C           LOOP OVER BASIS FUNCTIONS
            DO IBAS=1,NBAS
C
C             SPINOR COEFFICIENT FOR THIS BASIS FUNCTION
              CLP = C(     IBAS,NBLC+IOCC)
              IF(HMLT.NE.'NORL') THEN
                CSP = C(NBLC+IBAS,NBLC+IOCC)
                CLN = C(     IBAS,     IOCC)
                CSN = C(NBLC+IBAS,     IOCC)
              ENDIF
C
C             COPY INTO MASTER COEFFICIENT LIST
              COEF(IL1+IBAS,NSKP+IL1+IOCC) = DCMPLX(CLP,0.0D0)
              COEF(IL2+IBAS,NSKP+IL2+IOCC) = DCMPLX(CLP,0.0D0)
              IF(HMLT.NE.'NORL') THEN
                COEF(IS1+IBAS,NSKP+IL1+IOCC) = DCMPLX(CSP,0.0D0)
                COEF(IS2+IBAS,NSKP+IL2+IOCC) = DCMPLX(CSP,0.0D0)
                COEF(IL1+IBAS,     IL1+IOCC) = DCMPLX(CLN,0.0D0)
                COEF(IL2+IBAS,     IL2+IOCC) = DCMPLX(CLN,0.0D0)
                COEF(IS1+IBAS,     IL1+IOCC) = DCMPLX(CSN,0.0D0)
                COEF(IS2+IBAS,     IL2+IOCC) = DCMPLX(CSN,0.0D0)
              ENDIF
C
            ENDDO
C
          ENDDO
        ENDDO
C
C       IDENTIFY OCCUPIED STATE AND FILL IN DENSITY DETAILS
        IF(KQN.NE.KQNOCC) GOTO 150
C
C       IDENTIFY ADDRESS LOCATION OF OCCUPIED STATE
        IOCPN(IOCCML+1) = LRGE(IZ,KA,MADD)+NSHELL
        IOCTP(IZ,KA,MADD,NSHELL) = IOCPN(IOCCML+1)
C
C       INCREASE STORED OCCUPATION NUMBER
        IOCCML = IOCCML+1
C
C       ONE- AND TWO-BODY ENERGIES
        EH = EIG(NBLC+NSHELL)
        EK = 0.0D0
        EZ = 0.0D0
        ES = 0.0D0
        EE = 0.0D0
        EU = 0.0D0
        EP = 0.0D0
        EW = 0.0D0
        ER = 0.0D0
        ED = 0.0D0
C
C       ONE-BODY ENERGIES FOR OCCUPIED ELECTRONS
        M = 0
        DO IBAS=1,NBAS
          DO JBAS=1,NBAS
            M = M+1
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NBLC
            LBAS = JBAS+NBLC
C
            DLL = C(IBAS,NBLC+NSHELL)*C(JBAS,NBLC+NSHELL)
            DLS = C(IBAS,NBLC+NSHELL)*C(LBAS,NBLC+NSHELL)
            DSL = C(KBAS,NBLC+NSHELL)*C(JBAS,NBLC+NSHELL)
            DSS = C(KBAS,NBLC+NSHELL)*C(LBAS,NBLC+NSHELL)
C
C           KINETIC ENERGY: LL BLOCK
            IF(HMLT.EQ.'NORL') THEN
C
              EK = EK + T(IBAS,JBAS)*DLL
C
C           KINETIC ENERGY: SL BLOCK
            ELSE
C
              EK = EK + 2.0D0*T(KBAS,JBAS)*DSL
C
            ENDIF
C
C           NUCLEAR ATTRACTION
            EZ = EZ + V(IBAS,JBAS)*DLL + V(KBAS,LBAS)*DSS
C
C           ANOMALOUS ELECTRON MAGNETIC DIPOLE MOMENT
            EA = EA + A(IBAS,LBAS)*DLS + A(KBAS,JBAS)*DSL
C
C           SELF-INTERACTION (HIGH-ENERGY)
            ES = ES + S(IBAS,JBAS)*DLL + S(KBAS,LBAS)*DSS
C
C           SELF-INTERACTION (LOW-ENERGY)
            EY = EY + Y(IBAS,JBAS)*DLL + Y(IBAS,LBAS)*DLS
     &              + Y(KBAS,JBAS)*DSL + Y(KBAS,LBAS)*DSS
C
C           UEHLING NUCLEAR VACUUM POLARISATION (e+e-)
            EE = EE + E(IBAS,JBAS)*DLL + E(KBAS,LBAS)*DSS
C
C           UEHLING NUCLEAR VACUUM POLARISATION (μ+μ-)
            EU = EU + U(IBAS,JBAS)*DLL + U(KBAS,LBAS)*DSS
C
C           UEHLING NUCLEAR VACUUM POLARISATION (π+π-)
            EP = EP + P(IBAS,JBAS)*DLL + P(KBAS,LBAS)*DSS
C
C           WICHMANN-KROLL NUCLEAR VACUUM POLARISATION
            EW = EW + W(IBAS,JBAS)*DLL + W(KBAS,LBAS)*DSS
C
C           KÄLLÉN-SABRY NUCLEAR VACUUM POLARISATION
            ER = ER + R(IBAS,JBAS)*DLL + R(KBAS,LBAS)*DSS
C
C           UEHLING ELECTRONIC VACUUM POLARISATION
            ED = ED + D(IBAS,JBAS)*DLL + D(KBAS,LBAS)*DSS
C
          ENDDO
        ENDDO
C
150     CONTINUE
C
C     END LOOP OVER ALL AVAILABLE KQN VALUES
      ENDDO
C
      EG = 0.0D0
      EB = 0.0D0
      E1 = EK+EN+EZ+EA+ES+EY+EE+EU+EP+EW+ER+ED
      ENEW = EH
C
C     WRITE RESULT
      WRITE(6,41) 1,EH,0.0D0,ENEW,1.0D0
      WRITE(7,41) 1,EH,0.0D0,ENEW,1.0D0
C
C     STUPID RATIO THAT GRASP USES
      EVRL = (ENEW-EK)/EK
C
C**********************************************************************C
C     WRITTEN SUMMARY                                                  C
C**********************************************************************C
C
C     SUMMARY OF ENERGY CONTRIBUTIONS
50    FORMAT(1X,A,21X,F21.12)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',20),'Atomic energies (Hartree units)'
      WRITE(7, *) REPEAT(' ',20),'Atomic energies (Hartree units)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,50) 'Electron kinetic energy       ',EK
      WRITE(7,50) 'Electron kinetic energy       ',EK
      WRITE(6,50) 'Nuclear attraction energy     ',EZ
      WRITE(7,50) 'Nuclear attraction energy     ',EZ
      WRITE(6,50) 'Two-electron energy (Coulomb) ',EG
      WRITE(7,50) 'Two-electron energy (Coulomb) ',EG
      IF(HMLT.NE.'DHFB'.AND.HMLT.NE.'DHFQ') GOTO 500
      WRITE(6,50) 'Two-electron energy (Breit)   ',EB
      WRITE(7,50) 'Two-electron energy (Breit)   ',EB
      IF(HMLT.NE.'DHFQ') GOTO 500
      WRITE(6,50) 'Anomalous moment energy       ',EA
      WRITE(7,50) 'Anomalous moment energy       ',EA
      WRITE(6,50) 'Self-interaction energy (high)',ES
      WRITE(7,50) 'Self-interaction energy (high)',ES
      WRITE(6,50) 'Self-interaction energy (low) ',EY
      WRITE(7,50) 'Self-interaction energy (low) ',EY
      WRITE(6,50) 'Uehling energy (e+e-)         ',EE
      WRITE(7,50) 'Uehling energy (e+e-)         ',EE
      WRITE(6,50) 'Wichmann-Kroll energy         ',EW
      WRITE(7,50) 'Wichmann-Kroll energy         ',EW
      WRITE(6,50) 'Källén-Sabry energy           ',ER
      WRITE(7,50) 'Källén-Sabry energy           ',ER
      WRITE(6,50) 'Uehling energy (μ+μ-)         ',EU
      WRITE(7,50) 'Uehling energy (μ+μ-)         ',EU
      WRITE(6,50) 'Uehling energy (π+π-)         ',EP
      WRITE(7,50) 'Uehling energy (π+π-)         ',EP
      WRITE(6,50) 'Uehling SCF energy            ',ED
      WRITE(7,50) 'Uehling SCF energy            ',ED
500   CONTINUE
      WRITE(6,50) 'Total energy                  ',ENEW
      WRITE(7,50) 'Total energy                  ',ENEW
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,50) 'Virial nuclear/kinetic ratio  ',EVRL
      WRITE(7,50) 'Virial nuclear/kinetic ratio  ',EVRL
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
C
C     CONVERT BACK TO ORIGINAL ELECTRON MASS
      EMSS = XMSS
C
      RETURN
      END
C
C
      SUBROUTINE HFSCF0(IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          HH    HH FFFFFFFF SSSSSS   CCCCCC  FFFFFFFF 000000          C
C          HH    HH FF      SS    SS CC    CC FF      00   000         C
C          HH    HH FF      SS       CC       FF      00  0000         C
C          HHHHHHHH FFFFFF   SSSSSS  CC       FFFFFF  00 00 00         C
C          HH    HH FF            SS CC       FF      0000  00         C
C          HH    HH FF      SS    SS CC    CC FF      000   00         C
C          HH    HH FF       SSSSSS   CCCCCC  FF       000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  HFSCF PERFORMS AN ATOMIC SINGLE-DETERMINANT ITERATIVE SELF-         C
C  CONSISTENT FIELD PROCEDURE OVER THE USER-SPECIFIED HAMILTONIAN.     C
C  USES THE CLOSED-SHELL AVERAGE OF CONFIGURATION MODEL, WITH SUBSHELL C
C  OCCUPATIONS DETERMINED EITHER MANUALLY OR BY THE AUFBAU ROUTINE.    C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ IZ - ATOMIC CENTRE INDEX                                          C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*1 LLAB,QSGN
      CHARACTER*2 ELMT(120),ELNM
      CHARACTER*5 NMDL
      CHARACTER*6 CNFG
      CHARACTER*8 ZWRT,QWRT,EWRT
C
      DIMENSION QE(MKP),QA(MKP)
      DIMENSION EIG1(MBD),EIG2(MBD),TE(LWK0)
      DIMENSION O1(MBD,MBD),H1(MBD,MBD),C1(MBD,MBD),T1(MBD,MBD),
     &          V1(MBD,MBD),E1(MBD,MBD),U1(MBD,MBD),P1(MBD,MBD),
     &          W1(MBD,MBD),R1(MBD,MBD),D1(MBD,MBD),S1(MBD,MBD),
     &          Y1(MBD,MBD),A1(MBD,MBD)
      DIMENSION O2(MBD,MBD),H2(MBD,MBD),C2(MBD,MBD),T2(MBD,MBD),
     &          V2(MBD,MBD),E2(MBD,MBD),U2(MBD,MBD),P2(MBD,MBD),
     &          W2(MBD,MBD),R2(MBD,MBD),D2(MBD,MBD),S2(MBD,MBD),
     &          Y2(MBD,MBD),A2(MBD,MBD)
      DIMENSION DENLL(MB2,MKP),DFNLL(MB2,MKP),
     &          DENSL(MB2,MKP),DFNSL(MB2,MKP),
     &          DENSS(MB2,MKP),DFNSS(MB2,MKP),
     &          DENLS(MB2,MKP),DFNLS(MB2,MKP)
C
      COMPLEX*16 COEF(MDM,MDM)
C
      COMMON/ATMB/B11(MBD,MBD),B21(MBD,MBD),B12(MBD,MBD),B22(MBD,MBD)
      COMMON/ATMC/G11(MBD,MBD),G21(MBD,MBD),G12(MBD,MBD),G22(MBD,MBD)
      COMMON/ATMD/DLL1(MB2),DSL1(MB2),DSS1(MB2),DLS1(MB2),
     &            DLL2(MB2),DSL2(MB2),DSS2(MB2),DLS2(MB2)
      COMMON/AUFB/NORB(0:MEL,MKP+1),NUMOCC(0:MEL),LMXCONF
      COMMON/B0QN/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EWDR,EWXC,EANM,ESLF,EUEE,
     &            EUEU,EUEP,EWKR,EKSB,ESCF
      COMMON/FILL/NCNF(MCT,0:MEL,MKP+1),NLVL(MCT,0:MEL),CNFG(MCT)
      COMMON/MDLV/ELMT
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
C
C     CONVERGENCE TOLERANCE VALUE
      IF(HMLT.EQ.'NORL') THEN
        ENRGTOL = 5.0D-12
      ELSE
        ENRGTOL = 1.0D-12
C       ENRGTOL = 5.0D-13
      ENDIF
C
C     IMPORT ATOMIC CHARGE DETAILS
      IZNV = INT(ZNUC(IZ))
      ELNM = ELMT(IZNV)
      ICRG = IQNC(IZ)
      MLQN =(NKAP(IZ)-1)/2
C
C     CONVERT IZNV AND ICRG TO STRINGS AND WRITE TITLE
      IF(IZNV.LT.10) THEN
        WRITE(ZWRT,'(A,I1)') 'Z = ',IZNV
      ELSEIF(IZNV.LT.100) THEN
        WRITE(ZWRT,'(A,I2)') 'Z = ',IZNV
      ELSE
        WRITE(ZWRT,'(A,I3)') 'Z = ',IZNV
      ENDIF
C
      IF(IZNV-ICRG.LT.10) THEN
        WRITE(QWRT,'(A,I2)') 'Q = ',IZNV-ICRG
      ELSEIF(IZNV.LT.100) THEN
        WRITE(QWRT,'(A,I3)') 'Q = ',IZNV-ICRG
      ELSE
        WRITE(QWRT,'(A,I4)') 'Q = ',IZNV-ICRG
      ENDIF
C
      IF(IZNV-ICRG.GT.0) THEN
        QSGN = '+'
      ELSEIF(IZNV-ICRG.LT.0) THEN
        QSGN = '-'
      ENDIF
C
      ICMD = IABS(IZNV-ICRG)
      IF(IZNV-ICRG.EQ.0) THEN
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
20    FORMAT(17X,'Centre',I3,':',3X,A,3X,A,3X,A)
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,20) IZ,ZWRT,QWRT,EWRT
      WRITE(7,20) IZ,ZWRT,QWRT,EWRT
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     READ OR GENERATE AUFBAU OCCUPANCY FOR THIS ATOM
      IF(CNFG(IZ).EQ.'MANUAL') THEN
        LMXCONF = MLQN
        DO LQN=0,LMXCONF
          NUMOCC(LQN) = NLVL(IZ,LQN)
          DO NQN=1,NLVL(IZ,LQN)
            NORB(LQN,NQN) = NCNF(IZ,LQN,NQN)
          ENDDO
        ENDDO
      ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
        CALL AUFBAU(IZNV,ICRG,NORB,NUMOCC,LMXCONF)
      ELSE
        WRITE(6, *) 'In HFSCF0: invalid configuration choice.'
        WRITE(7, *) 'In HFSCF0: invalid configuration choice.'
      ENDIF
C
C     IDENTIFY THE HIGHEST OCCUPIED SHELL
      NMAX = 1
      DO LQN=0,LMXCONF
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
31    FORMAT(1X,'LQN = ',I1,2X,I3,2X,'|'1X,' OCC(',A,'):',A,12(2X,I2))
C
      IF(CNFG(IZ).EQ.'MANUAL') THEN
        WRITE(6,30) 'Manual:','#fns',(N,N=1,12)
        WRITE(7,30) 'Manual:','#fns',(N,N=1,12)
      ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
        WRITE(6,30) 'Aufbau:','#fns',(N,N=1,12)
        WRITE(7,30) 'Aufbau:','#fns',(N,N=1,12)
      ENDIF
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO LQN=0,MLQN
        WRITE(6,31) LQN,NFNC(LQN,IZ),LLAB(LQN),
     &                   REPEAT(' ',LQN*4),(NORB(LQN,J),J=1,NUMOCC(LQN))
        WRITE(7,31) LQN,NFNC(LQN,IZ),LLAB(LQN),
     &                   REPEAT(' ',LQN*4),(NORB(LQN,J),J=1,NUMOCC(LQN))
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
C     EARLY EXIT OPTION (NO ELECTRONS IN THIS ATOM)
      IF(ICRG.EQ.0) THEN
        WRITE(6,41) 1,0.0D0,0.0D0,0.0D0,1.0D0
        WRITE(7,41) 1,0.0D0,0.0D0,0.0D0,1.0D0
        WRITE(6, *) REPEAT('=',72)
        WRITE(7, *) REPEAT('=',72)
        WRITE(6, *) REPEAT(' ',72)
        WRITE(7, *) REPEAT(' ',72)
        RETURN
      ENDIF
C
C     INITIALISE A STORAGE BIN FOR PREVIOUS ATOMIC ENERGY
      EPRV = 0.0D0
C
C     TWO-BODY PROBLEM: INTERACTING ELECTRONS. (TREAT WITH SCF.)
      DO 1000 ITER=1,MIT0
C
C       INITIALISE ONE-BODY AND TWO-BODY ENERGY COUNTERS
        EH = 0.0D0
        EK = 0.0D0
        EZ = 0.0D0
        EA = 0.0D0
        ES = 0.0D0
        EY = 0.0D0
        EE = 0.0D0
        EU = 0.0D0
        EP = 0.0D0
        EW = 0.0D0
        ER = 0.0D0
        ED = 0.0D0
        EG = 0.0D0
        EB = 0.0D0
C
C       INITIALISE ELECTRON OCCUPATION COUNTER
        IOCCML = IOCCM0
C
C**********************************************************************C
C     ONE-ELECTRON PART: LOOP OVER BASIS FUNCTIONS I,J (USE INDEX 100) C
C**********************************************************************C
C
C     LOOP OVER ALL OCCUPIED LQN VALUES
C     DO 100 LQNA=0,LMXCONF
C     LOOP OVER ALL AVAILABLE LQN VALUES
      DO 100 LQNA=0,MLQN
C
C     RECORD LQNA VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
      NBASA = NFNC(LQNA,IZ)
      DO IBAS=1,NBASA
        EXLA(IBAS) = BEXL(IBAS,LQNA,IZ)
      ENDDO
C
C     MATRIX DIMENSIONS FOR THIS LQNA BLOCK
      IF(HMLT.EQ.'NORL') THEN
        NBLC = 0
      ELSE
        NBLC = NBASA
      ENDIF
      NMAT = NBASA+NBLC
C
C     EFFECTIVE AND AVERAGE OCCUPATION NUMBERS FOR THIS LQNA ORBITAL
C     A CLOSED SUBSHELL (NSHELL,LQNA) CONTAINS NCLS ELECTRONS
      NCLS = 4*LQNA+2
C
C     FOR EACH LISTED SUBSHELLS OF LQN TYPE
      DO IOCC=1,NUMOCC(LQNA)
C
C       NUMBER OF CHARGES IN THIS SUBSHELL (NSHELL=IOCC+LQNA)
        NQ = NORB(LQNA,IOCC)
C
        IF(NQ.EQ.NCLS) THEN
C         IF SUBSHELL IS CLOSED THERE IS NO FRACTIONAL OCCUPANCY
          QE(IOCC) = 1.0D0
        ELSE
C         IF SUBSHELL IS OPEN, CONSTRUCT FRACTION (GRANT 6.6.24)
          QE(IOCC) = DFLOAT(NQ-1)/DFLOAT(NCLS-1)
        ENDIF
C
C       ACTUAL FRACTIONAL SUBSHELL OCCUPANCY
        IF(NQ.GT.0) THEN
          QA(IOCC) = DFLOAT(NQ)/DFLOAT(NCLS)
        ELSE
          QA(IOCC) = 1.0D0
        ENDIF
C
      ENDDO
C
C     POSITIVE KAPPA(A) CHOICE (APPLIES ONLY FOR LQNA > 0)
      IF(LQNA.EQ.0.OR.HMLT.EQ.'NORL') GOTO 130

      KAPA1 = LQNA
      RK2A1 = DFLOAT(2*IABS(KAPA1))
C
C     INITIALISE HAMILTONIAN AND OVERLAP ARRAYS
      DO IBAS=1,NMAT
        DO JBAS=1,NMAT
          O1(IBAS,JBAS) = 0.0D0
          H1(IBAS,JBAS) = 0.0D0
          C1(IBAS,JBAS) = 0.0D0
          T1(IBAS,JBAS) = 0.0D0
          V1(IBAS,JBAS) = 0.0D0
          S1(IBAS,JBAS) = 0.0D0
          Y1(IBAS,JBAS) = 0.0D0
          E1(IBAS,JBAS) = 0.0D0
          U1(IBAS,JBAS) = 0.0D0
          P1(IBAS,JBAS) = 0.0D0
          W1(IBAS,JBAS) = 0.0D0
          R1(IBAS,JBAS) = 0.0D0
          D1(IBAS,JBAS) = 0.0D0
          A1(IBAS,JBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     BASIS FUNCTION OVERLAPS
      CALL OVRLP0(O1,EXLA,KAPA1,NBASA)
C
C     KINETIC ENERGY
      CALL KINTC0(T1,EXLA,KAPA1,NBASA)
C
C     NUCLEAR ATTRACTION
      IF(NMDL(IZ).EQ.'POINT') THEN
        CALL PNTNC0(V1,EXLA,IZ,KAPA1,NBASA)
      ELSEIF(NMDL(IZ).EQ.'GAUSS') THEN
        CALL GSSNC0(V1,EXLA,IZ,KAPA1,NBASA)
      ELSEIF(NMDL(IZ).EQ.'FERMI') THEN
        CALL FMINC0(V1,EXLA,IZ,KAPA1,NBASA)
        CALL FMINC0temp(V1,EXLA,IZ,KAPA1,NBASA)
      ELSEIF(NMDL(IZ).EQ.'UNIFM') THEN
        CALL UNINC0(V1,EXLA,IZ,KAPA1,NBASA)
      ENDIF
C
C     FOR NUCLEAR BEST-FIT TESTING PURPOSES
C     CALL GSMNC0(V1,EXLA,IZ,KAPA1,NBASA)
C     CALL ONEEL0(T1,V1,EXLA,IZ,KAPA1,NBASA)
C
C     QED-LEVEL INTERACTIONS (ONE-BODY ONLY)
      IF(HMLT.EQ.'DHFQ') THEN
c       CALL ANMLS0(A1,EXLA,IZ,KAPA1,NBASA)
c       CALL SLFLW0(Y1,     IZ,KAPA1      )
c       CALL SLFHI0(S1,EXLA,IZ,KAPA1,NBASA)
        CALL UEHNC0(E1,EXLA,IZ,KAPA1,NBASA,'ELEC')
C       CALL UEHNC0(U1,EXLA,IZ,KAPA1,NBASA,'MUON')
C       CALL UEHNC0(P1,EXLA,IZ,KAPA1,NBASA,'PION')
        CALL WKRNC0(W1,EXLA,IZ,KAPA1,NBASA)
        CALL KSBNC0(R1,EXLA,IZ,KAPA1,NBASA)
C       CALL UEHEL0(D1,EXLA,IZ,KAPA1,NBASA)
      ENDIF
C
C     ONE-ELECTRON HAMILTONIAN MATRIX
      DO IBAS=1,NMAT
        DO JBAS=1,NMAT
          H1(IBAS,JBAS) = T1(IBAS,JBAS)+V1(IBAS,JBAS)+A1(IBAS,JBAS)
     &                  + Y1(IBAS,JBAS)+S1(IBAS,JBAS)+E1(IBAS,JBAS)
     &                  + U1(IBAS,JBAS)+P1(IBAS,JBAS)+W1(IBAS,JBAS)
     &                  + R1(IBAS,JBAS)+D1(IBAS,JBAS)
        ENDDO
      ENDDO
C
130   CONTINUE
C
C     NEGATIVE KAPPA(A) CHOICE (APPLIES TO ALL LQNA VALUES)
      KAPA2 =-LQNA-1
      IF(HMLT.EQ.'NORL') THEN
        RK2A2 = DFLOAT(NCLS)
      ELSE
        RK2A2 = DFLOAT(2*IABS(KAPA2))
      ENDIF
C
C     INITIALISE HAMILTONIAN AND OVERLAP ARRAYS
      DO IBAS=1,NMAT
        DO JBAS=1,NMAT
          O2(IBAS,JBAS) = 0.0D0
          H2(IBAS,JBAS) = 0.0D0
          C2(IBAS,JBAS) = 0.0D0
          T2(IBAS,JBAS) = 0.0D0
          V2(IBAS,JBAS) = 0.0D0
          S2(IBAS,JBAS) = 0.0D0
          Y2(IBAS,JBAS) = 0.0D0
          E2(IBAS,JBAS) = 0.0D0
          U2(IBAS,JBAS) = 0.0D0
          P2(IBAS,JBAS) = 0.0D0
          W2(IBAS,JBAS) = 0.0D0
          R2(IBAS,JBAS) = 0.0D0
          D2(IBAS,JBAS) = 0.0D0
          A2(IBAS,JBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     BASIS FUNCTION OVERLAPS
      CALL OVRLP0(O2,EXLA,KAPA2,NBASA)
C
C     KINETIC ENERGY
      CALL KINTC0(T2,EXLA,KAPA2,NBASA)
C
C     NUCLEAR ATTRACTION
      IF(NMDL(IZ).EQ.'POINT') THEN
        CALL PNTNC0(V2,EXLA,IZ,KAPA2,NBASA)
      ELSEIF(NMDL(IZ).EQ.'GAUSS') THEN
        CALL GSSNC0(V2,EXLA,IZ,KAPA2,NBASA)
      ELSEIF(NMDL(IZ).EQ.'FERMI') THEN
        CALL FMINC0(V2,EXLA,IZ,KAPA2,NBASA)
        CALL FMINC0temp(V2,EXLA,IZ,KAPA2,NBASA)
      ELSEIF(NMDL(IZ).EQ.'UNIFM') THEN
        CALL UNINC0(V2,EXLA,IZ,KAPA2,NBASA)
      ENDIF
C
C     FOR NUCLEAR BEST-FIT TESTING PURPOSES
C     CALL GSMNC0(V2,EXLA,IZ,KAPA2,NBASA)
C     CALL ONEEL0(T2,V2,EXLA,IZ,KAPA2,NBASA)
C
C     QED-LEVEL INTERACTIONS (ONE-BODY ONLY)
      IF(HMLT.EQ.'DHFQ') THEN
C       CALL ANMLS0(A2,EXLA,IZ,KAPA2,NBASA)
C       CALL SLFLW0(Y2,     IZ,KAPA2      )
C       CALL SLFHI0(S2,EXLA,IZ,KAPA2,NBASA)
        CALL UEHNC0(E2,EXLA,IZ,KAPA2,NBASA,'ELEC')
C       CALL UEHNC0(U2,EXLA,IZ,KAPA2,NBASA,'MUON')
C       CALL UEHNC0(P2,EXLA,IZ,KAPA2,NBASA,'PION')
        CALL WKRNC0(W2,EXLA,IZ,KAPA2,NBASA)
        CALL KSBNC0(R2,EXLA,IZ,KAPA2,NBASA)
C       CALL UEHEL0(D2,EXLA,IZ,KAPA2,NBASA)
      ENDIF
C
C     ONE-ELECTRON HAMILTONIAN MATRIX
      DO IBAS=1,NMAT
        DO JBAS=1,NMAT
          H2(IBAS,JBAS) = T2(IBAS,JBAS)+V2(IBAS,JBAS)+A2(IBAS,JBAS)
     &                  + S2(IBAS,JBAS)+Y2(IBAS,JBAS)+E2(IBAS,JBAS)
     &                  + U2(IBAS,JBAS)+P2(IBAS,JBAS)+W2(IBAS,JBAS)
     &                  + R2(IBAS,JBAS)+D2(IBAS,JBAS)
        ENDDO
      ENDDO
C
C     REASONS TO SKIP SCF INTERACTIONS
      IF(ITER.EQ.1.OR.ICRG.EQ.1) GOTO 150
C
C     IF SELF-ENERGY IS NEEDED, ALWAYS CALCULATE SCF ELEMENTS
      IF(HMLT.EQ.'DHFP'.OR.HMLT.EQ.'DHFQ') GOTO 180

C     IF MBPT2 IS NEEDED, ALWAYS CALCULATE SCF ELEMENTS
      IF(TREE.EQ.'MBPTN') GOTO 180
C
C     SKIP SCF INTERACTIONS IF LQN LEVEL IS UNOCCUPIED
      IF(LQNA.GT.LMXCONF) GOTO 150
C
C     SKIP POINT WHEN SCF ELEMENTS ARE REQUIRED
180   CONTINUE
C
C     INITIALISE RELEVANT COUNTERS AND ARRAYS
      RK2B1 = 0.0D0
      RK2B2 = 0.0D0
C
C**********************************************************************C
C     TWO-ELECTRON PART: LOOP OVER BASIS FUNCTIONS K,L (USE INDEX 200) C
C**********************************************************************C
C
C     LOOP OVER ALL OCCUPIED LQN VALUES
      DO 200 LQNB=0,LMXCONF
C
C     RECORD LQNB VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
      NBASB = NFNC(LQNB,IZ)
      DO JBAS=1,NBASB
        EXLB(JBAS) = BEXL(JBAS,LQNB,IZ)
      ENDDO
C
C     NUMBER OF BASIS FUNCTION OVERLAPS IN THIS BLOCK
      MAXM = NBASB*NBASB
C
C     POSITIVE KAPPA(B) CHOICE (APPLIES ONLY FOR LQNB > 0)
      IF(LQNB.EQ.0.OR.HMLT.EQ.'NORL') GOTO 230
C
C     ANGULAR QUANTUM NUMBER AND DEGENERACY
      KAPB1 = LQNB
      RK2B1 = DFLOAT(2*IABS(KAPB1))
C
C     RELEVANT DENSITY MATRIX DEPENDS ON LQNA AND LQNB
      IF(LQNA.EQ.LQNB) THEN
        DO M=1,MAXM
          DLL1(M) = DENLL(M,2*LQNB  )
          DSL1(M) = DENSL(M,2*LQNB  )
          DSS1(M) = DENSS(M,2*LQNB  )
          DLS1(M) = DENLS(M,2*LQNB  )
        ENDDO
      ELSEIF(LQNA.NE.LQNB) THEN
        DO M=1,MAXM
          DLL1(M) = DFNLL(M,2*LQNB  )
          DSL1(M) = DFNSL(M,2*LQNB  )
          DSS1(M) = DFNSS(M,2*LQNB  )
          DLS1(M) = DFNLS(M,2*LQNB  )
        ENDDO
      ENDIF
C
230   CONTINUE
C
C     NEGATIVE KAPPA(B) CHOICE (APPLIES TO ALL LQNB VALUES)
C
C     ANGULAR QUANTUM NUMBER AND DEGENERACY
      KAPB2 =-LQNB-1
      IF(HMLT.EQ.'NORL') THEN
        RK2B2 = DFLOAT(4*LQNB+2)
      ELSE
        RK2B2 = DFLOAT(2*IABS(KAPB2))
      ENDIF
C
C     RELEVANT DENSITY MATRIX DEPENDS ON LQNA AND LQNB
      IF(LQNA.EQ.LQNB) THEN
        DO M=1,MAXM
          DLL2(M) = DENLL(M,2*LQNB+1)
          IF(HMLT.NE.'NORL') THEN
            DSL2(M) = DENSL(M,2*LQNB+1)
            DSS2(M) = DENSS(M,2*LQNB+1)
            DLS2(M) = DENLS(M,2*LQNB+1)
          ENDIF
        ENDDO
      ELSEIF(LQNA.NE.LQNB) THEN
        DO M=1,MAXM
          DLL2(M) = DFNLL(M,2*LQNB+1)
          IF(HMLT.NE.'NORL') THEN
            DSL2(M) = DFNSL(M,2*LQNB+1)
            DSS2(M) = DFNSS(M,2*LQNB+1)
            DLS2(M) = DFNLS(M,2*LQNB+1)
          ENDIF
        ENDDO
      ENDIF
C
C**********************************************************************C
C     GENERATE ATOMIC MEAN-FIELD COULOMB MATRIX                        C
C**********************************************************************C
C
C     GENERATE THE MEAN-FIELD ATOMIC COULOMB MATRIX OVER DENSITIES
      CALL COULM0
C
C     ADD COULOMB MATRIX ELEMENTS CONTRIBUTIONS TO FOCK MATRIX
      DO IBAS=1,NMAT
        DO JBAS=1,NMAT
C
          IF(HMLT.EQ.'NORL') THEN
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
C
C     TWO-BODY EIGENVALUE ENERGIES FOR OCCUPIED ELECTRONS
      M = 0
      DO IBAS=1,NBASA
        DO JBAS=1,NBASA
          M = M+1
C
          IF(HMLT.EQ.'NORL') THEN
C         NON-RELATIVISTIC COULOMB ENERGY
C
            EG = EG + RK2A2*RK2B2*G22(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
C
          ELSE
C         RELATIVISTIC COULOMB ENERGY
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NBASA
            LBAS = JBAS+NBASA
C
C           CONTRIBUTIONS THAT ALWAYS MATTER
            EG = EG
     &         +       RK2A2*RK2B1*G21(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &         +       RK2A2*RK2B2*G22(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &         + 2.0D0*RK2A2*RK2B1*G21(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
     &         + 2.0D0*RK2A2*RK2B2*G22(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
     &         +       RK2A2*RK2B1*G21(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
     &         +       RK2A2*RK2B2*G22(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
C           CONTRIBUTIONS THAT ONLY MATTER WHEN LQNA>0
            IF(LQNA.EQ.0) GOTO 235
C
            EG = EG
     &         +       RK2A1*RK2B1*G11(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &         +       RK2A1*RK2B2*G12(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &         + 2.0D0*RK2A1*RK2B1*G11(KBAS,JBAS)*DFNSL(M,2*LQNA  )
     &         + 2.0D0*RK2A1*RK2B2*G12(KBAS,JBAS)*DFNSL(M,2*LQNA  )
     &         +       RK2A1*RK2B1*G11(KBAS,LBAS)*DFNSS(M,2*LQNA  )
     &         +       RK2A1*RK2B2*G12(KBAS,LBAS)*DFNSS(M,2*LQNA  )
C
235         CONTINUE
C
          ENDIF
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     GENERATE ATOMIC MEAN-FIELD BREIT MATRIX                          C
C**********************************************************************C
C
C     GENERATE THE MEAN-FIELD ATOMIC BREIT MATRIX
      IF(HMLT.NE.'DHFB'.AND.HMLT.NE.'DHFQ') GOTO 240
C
C     GENERATE THE MEAN-FIELD ATOMIC BREIT MATRIX OVER DENSITIES
      CALL BREIT0
C
C     ADD TWO-PARTICLE CONTRIBUTIONS TO FOCK MATRIX
      DO IBAS=1,NMAT
        DO JBAS=1,NMAT
C
          H1(IBAS,JBAS) = H1(IBAS,JBAS) + RK2B1*B11(IBAS,JBAS)
     &                                  + RK2B2*B12(IBAS,JBAS)
          H2(IBAS,JBAS) = H2(IBAS,JBAS) + RK2B1*B21(IBAS,JBAS)
     &                                  + RK2B2*B22(IBAS,JBAS)
C
        ENDDO
      ENDDO
C
C     TWO-BODY EIGENVALUE ENERGIES FOR OCCUPIED ELECTRONS
      M = 0
      DO IBAS=1,NBASA
        DO JBAS=1,NBASA
          M = M+1
C
C         SMALL COMPONENT BLOCK ADDRESSES
          KBAS = IBAS+NBASA
          LBAS = JBAS+NBASA
C
C         CONTRIBUTIONS THAT ALWAYS MATTER
          EB = EB
     &       +       RK2A2*RK2B1*B21(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &       +       RK2A2*RK2B2*B22(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &       + 2.0D0*RK2A2*RK2B1*B21(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
     &       + 2.0D0*RK2A2*RK2B2*B22(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
     &       +       RK2A2*RK2B1*B21(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
     &       +       RK2A2*RK2B2*B22(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
C         CONTRIBUTIONS THAT ONLY MATTER WHEN LQNA>0
          IF(LQNA.EQ.0) GOTO 245
C
C         LL BLOCK
          EB = EB
     &       +       RK2A1*RK2B1*B11(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &       +       RK2A1*RK2B2*B12(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &       + 2.0D0*RK2A1*RK2B1*B11(KBAS,JBAS)*DFNSL(M,2*LQNA  )
     &       + 2.0D0*RK2A1*RK2B2*B12(KBAS,JBAS)*DFNSL(M,2*LQNA  )
     &       +       RK2A1*RK2B1*B11(KBAS,LBAS)*DFNSS(M,2*LQNA  )
     &       +       RK2A1*RK2B2*B12(KBAS,LBAS)*DFNSS(M,2*LQNA  )
C
245       CONTINUE
C
        ENDDO
      ENDDO
C
240   CONTINUE
C
C     END LOOP OVER LQNS FOR ORBITAL B
200   CONTINUE
C
C     FINISHED CALCULATING OVERLAP COMBINATIONS BETWEEN THIS LQNA
C     VALUE AND ALL POSSIBLE LQNB VALUES
C
C**********************************************************************C
C     FINISHED GENERATING TWO-ELECTRON INTEGRALS                       C
C**********************************************************************C
C
C     SKIP POINT OUT OF SCF CODE
150   CONTINUE
C
C     ONE-BODY ENERGIES FOR OCCUPIED ELECTRONS
      M = 0
      DO IBAS=1,NBASA
        DO JBAS=1,NBASA
          M = M+1
C
C         SMALL COMPONENT BLOCK ADDRESSES
          KBAS = IBAS+NBASA
          LBAS = JBAS+NBASA
C
C         KINETIC ENERGY: LL BLOCK
          IF(HMLT.EQ.'NORL') THEN
C
            EK = EK + RK2A2*T2(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
C
C         KINETIC ENERGY: SL BLOCK
          ELSE
C
            EK = EK + 2.0D0*RK2A2*T2(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
C
            IF(LQNA.NE.0) THEN
              EK = EK + 2.0D0*RK2A1*T1(KBAS,JBAS)*DFNSL(M,2*LQNA  )
            ENDIF
C
          ENDIF
C
C         NUCLEAR ATTRACTION
          EZ = EZ + RK2A2*V2(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &            + RK2A2*V2(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
C         ANOMALOUS ELECTRON MAGNETIC DIPOLE MOMENT
          EA = EA + RK2A2*A2(IBAS,LBAS)*DFNLS(M,2*LQNA+1)
     &            + RK2A2*A2(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
C
C         SELF-INTERACTION (HIGH-ENERGY)
          ES = ES + RK2A2*S2(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &            + RK2A2*S2(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
C         SELF-INTERACTION (LOW-ENERGY)
          EY = EY + RK2A2*Y2(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &            + RK2A2*Y2(IBAS,LBAS)*DFNLS(M,2*LQNA+1)
     &            + RK2A2*Y2(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
     &            + RK2A2*Y2(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
C         UEHLING NUCLEAR VACUUM POLARISATION (e+e-)
          EE = EE + RK2A2*E2(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &            + RK2A2*E2(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
C         UEHLING NUCLEAR VACUUM POLARISATION (μ+μ-)
          EU = EU + RK2A2*U2(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &            + RK2A2*U2(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
C         UEHLING NUCLEAR VACUUM POLARISATION (π+π-)
          EP = EP + RK2A2*P2(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &            + RK2A2*P2(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
C         WICHMANN-KROLL NUCLEAR VACUUM POLARISATION
          EW = EW + RK2A2*W2(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &            + RK2A2*W2(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
C         KÄLLÉN-SABRY NUCLEAR VACUUM POLARISATION
          ER = ER + RK2A2*R2(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &            + RK2A2*R2(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
C         UEHLING ELECTRONIC VACUUM POLARISATION
          ED = ED + RK2A2*D2(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &            + RK2A2*D2(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
C         CONTRIBUTIONS THAT ONLY MATTER WHEN LQNA>0
          IF(LQNA.EQ.0) GOTO 205
C
C         NUCLEAR ATTRACTION
          EZ = EZ + RK2A1*V1(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &            + RK2A1*V1(KBAS,LBAS)*DFNSS(M,2*LQNA  )
C
C         ANOMALOUS ELECTRON MAGNETIC DIPOLE MOMENT
          EA = EA + RK2A1*A1(IBAS,LBAS)*DFNLS(M,2*LQNA  )
     &            + RK2A1*A1(KBAS,JBAS)*DFNSL(M,2*LQNA  )
C
C         SELF-INTERACTION (HIGH-ENERGY)
          ES = ES + RK2A1*S1(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &            + RK2A1*S1(KBAS,LBAS)*DFNSS(M,2*LQNA  )
C
C         SELF-INTERACTION (LOW-ENERGY)
          EY = EY + RK2A1*Y1(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &            + RK2A1*Y1(IBAS,LBAS)*DFNLS(M,2*LQNA  )
     &            + RK2A1*Y1(KBAS,JBAS)*DFNSL(M,2*LQNA  )
     &            + RK2A1*Y1(KBAS,LBAS)*DFNSS(M,2*LQNA  )
C
C         UEHLING NUCLEAR VACUUM POLARISATION (e+e-)
          EE = EE + RK2A1*E1(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &            + RK2A1*E1(KBAS,LBAS)*DFNSS(M,2*LQNA  )
C
C         UEHLING NUCLEAR VACUUM POLARISATION (μ+μ-)
          EU = EU + RK2A1*U1(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &            + RK2A1*U1(KBAS,LBAS)*DFNSS(M,2*LQNA  )
C
C         UEHLING NUCLEAR VACUUM POLARISATION (π+π-)
          EP = EP + RK2A1*P1(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &            + RK2A1*P1(KBAS,LBAS)*DFNSS(M,2*LQNA  )
C
C         WICHMANN-KROLL NUCLEAR VACUUM POLARISATION
          EW = EW + RK2A1*W1(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &            + RK2A1*W1(KBAS,LBAS)*DFNSS(M,2*LQNA  )
C
C         KÄLLÉN-SABRY NUCLEAR VACUUM POLARISATION
          ER = ER + RK2A1*R1(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &            + RK2A1*R1(KBAS,LBAS)*DFNSS(M,2*LQNA  )
C
C         UEHLING ELECTRONIC VACUUM POLARISATION
          ED = ED + RK2A1*D1(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &            + RK2A1*D1(KBAS,LBAS)*DFNSS(M,2*LQNA  )
C
205       CONTINUE
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     MATRIX DIAGONALISATION AND COEFFICIENT MATRIX UPDATES            C
C**********************************************************************C
C
C     POSITIVE KAPPA(A) CHOICE (APPLIES ONLY FOR LQNA > 0)
      IF(LQNA.EQ.0.OR.HMLT.EQ.'NORL') GOTO 140
C
C     WORKING EIGENVECTOR MATRIX
      DO IBAS=1,NMAT
        DO JBAS=1,NMAT
          C1(IBAS,JBAS) = H1(IBAS,JBAS)
        ENDDO
      ENDDO
C
C     DIAGONALISE FOCK MATRIX (THIS NEEDS LAPACK LIBRARY)
      CALL DSYGV(1,'V','U',NMAT,C1,MBD,O1,MBD,EIG1,TE,LWK0,INFO)
      IF(INFO.NE.0) THEN
        WRITE(6, *) 'In HFSCF0: eigenvalue solver DSYGV(1) failed.',INFO
        WRITE(7, *) 'In HFSCF0: eigenvalue solver DSYGV(1) failed.',INFO
        STOP
      ENDIF
c
c      diagonalise the overlap matrix to inspect linear dependence      
c      CALL DSYEV('V','U',NMAT,O1,MBD,EIG1,TE,LWK0,INFO)
c      IF(INFO.NE.0) THEN
c        WRITE(6, *) 'In HFSCF0: eigenvalue solver DSYEV(1) failed.',INFO
c        WRITE(7, *) 'In HFSCF0: eigenvalue solver DSYEV(1) failed.',INFO
c        STOP
c      ENDIF
c      write(*,*) '----'
c333   format(i3,1x,es24.14,1x,es24.14)
c      write(*,*) KAPA1
c      do n=1,NBLC
c        write(*,333) n,EIG1(n),EIG1(n+nblc)
c      enddo
C
C     ATOMIC SELECTION RULE: ORTHOGONALITY IN BLOCKS OF KQN -> KA = KB
C
C     BEGIN LOOP OVER ALL MQNA VALUES
      DO IMVAL=1,IABS(KAPA1)
C
C       COEFFICIENT MATRIX ADDRESSES
        IL1 = LRGE(IZ,2*LQNA  ,2*IMVAL-1)
        IL2 = LRGE(IZ,2*LQNA  ,2*IMVAL  )
        IS1 = LRGE(IZ,2*LQNA  ,2*IMVAL-1)+NSKP
        IS2 = LRGE(IZ,2*LQNA  ,2*IMVAL  )+NSKP
C
C       LOOP OVER OCCUPIED SOLUTIONS OF THE EIGENVALUE PROBLEM
c       DO IOCC=1,NUMOCC(LQNA)
C       LOOP OVER ALL SOLUTIONS OF THE EIGENVALUE PROBLEM
        DO IOCC=1,NBASA
C
C         COPY ENERGY EIGENVALUES TO MASTER LIST
          EIGN(NSKP+IL1+IOCC) = EIG1(NBLC+IOCC)
          EIGN(NSKP+IL2+IOCC) = EIG1(NBLC+IOCC)
          EIGN(     IL1+IOCC) = EIG1(     IOCC)
          EIGN(     IL2+IOCC) = EIG1(     IOCC)
C
C         EFFECTIVE FRACTIONAL OCCUPANCY FOR THIS ORBITAL
          IF(IOCC.LE.NUMOCC(LQNA)) THEN
            QF = DSQRT(QA(IOCC))
          ELSE
            QF = 1.0D0
          ENDIF
C
C         LOOP OVER BASIS FUNCTIONS
          DO IBAS=1,NBASA
C
C           KRAMERS PAIR LARGE- AND SMALL-COMPONENT COEFFICIENTS
            CLP = QF*C1(     IBAS,NBLC+IOCC)
            CSP = QF*C1(NBLC+IBAS,NBLC+IOCC)
            CLN = QF*C1(     IBAS,     IOCC)
            CSN = QF*C1(NBLC+IBAS,     IOCC)
C
C           COPY INTO MASTER COEFFICIENT LIST
            COEF(IL1+IBAS,NSKP+IL1+IOCC) = DCMPLX(CLP,0.0D0)
            COEF(IL2+IBAS,NSKP+IL2+IOCC) = DCMPLX(CLP,0.0D0)
            COEF(IS1+IBAS,NSKP+IL1+IOCC) = DCMPLX(CSP,0.0D0)
            COEF(IS2+IBAS,NSKP+IL2+IOCC) = DCMPLX(CSP,0.0D0)
            COEF(IL1+IBAS,     IL1+IOCC) = DCMPLX(CLN,0.0D0)
            COEF(IL2+IBAS,     IL2+IOCC) = DCMPLX(CLN,0.0D0)
            COEF(IS1+IBAS,     IL1+IOCC) = DCMPLX(CSN,0.0D0)
            COEF(IS2+IBAS,     IL2+IOCC) = DCMPLX(CSN,0.0D0)
C          
          ENDDO
        ENDDO
C
C       LOOP OVER ALL OCCUPIED SUBSHELLS OF THIS KQN TYPE GIVEN MQN
        DO IOCC=1,NUMOCC(LQNA)
C
C         SKIP UNOCCUPIED ORBITALS
          IF(NORB(LQNA,IOCC).NE.0) THEN
C
C           FOCK ADDRESS FOR THIS PAIR
            IOCPN(IOCCML+1) = LRGE(IZ,2*LQNA  ,2*IMVAL-1)+IOCC
            IOCPN(IOCCML+2) = LRGE(IZ,2*LQNA  ,2*IMVAL  )+IOCC
C
C           SAVE THESE STARTING ADDRESSES
            IOCTP(IZ,2*LQNA  ,2*IMVAL-1,IOCC) = IOCPN(IOCCML+1)
            IOCTP(IZ,2*LQNA  ,2*IMVAL  ,IOCC) = IOCPN(IOCCML+2)
C
C           INCREASE FOCK ADDRESS OF OCCUPIED ORBITALS (PAIR AT A TIME)
            IOCCML = IOCCML+2
C
          ENDIF
C
        ENDDO
C
      ENDDO
C
C     BUILD ATOMIC CHARGE DENSITY MATRIX FOR THIS KQNA BLOCK
      M = 0
      DO IBAS=1,NBASA
        DO JBAS=1,NBASA
          M = M+1
C
C         INITIALISE ATOMIC DENSITY LISTS FOR THIS BLOCK
          DENLL(M,2*LQNA  ) = 0.0D0
          DENSL(M,2*LQNA  ) = 0.0D0
          DENSS(M,2*LQNA  ) = 0.0D0
          DENLS(M,2*LQNA  ) = 0.0D0
C
          DFNLL(M,2*LQNA  ) = 0.0D0
          DFNSL(M,2*LQNA  ) = 0.0D0
          DFNSS(M,2*LQNA  ) = 0.0D0
          DFNLS(M,2*LQNA  ) = 0.0D0
C
C         LOOP OVER ALL LISTED SUBSHELLS OF THIS KQN TYPE
          DO IOCC=1,NUMOCC(LQNA)
C
C           SKIP THIS STEP SUBSHELL IS UNOCCUPIED
            IF(NORB(LQNA,IOCC).GT.0) THEN
            
C             DENSITY OVERLAPS FROM EIGENVECTOR PRODUCTS
              DLL = C1(IBAS     ,NBLC+IOCC)*C1(JBAS     ,NBLC+IOCC)
              DSL = C1(IBAS+NBLC,NBLC+IOCC)*C1(JBAS     ,NBLC+IOCC)
              DSS = C1(IBAS+NBLC,NBLC+IOCC)*C1(JBAS+NBLC,NBLC+IOCC)
              DLS = C1(IBAS     ,NBLC+IOCC)*C1(JBAS+NBLC,NBLC+IOCC)
C
C             ADD DENSITY CONTRIBUTIONS TO ATOMIC LIST
              DENLL(M,2*LQNA  ) = DENLL(M,2*LQNA  ) + QE(IOCC)*DLL
              DENSL(M,2*LQNA  ) = DENSL(M,2*LQNA  ) + QE(IOCC)*DSL
              DENSS(M,2*LQNA  ) = DENSS(M,2*LQNA  ) + QE(IOCC)*DSS
              DENLS(M,2*LQNA  ) = DENLS(M,2*LQNA  ) + QE(IOCC)*DLS
C
              DFNLL(M,2*LQNA  ) = DFNLL(M,2*LQNA  ) + QA(IOCC)*DLL
              DFNSL(M,2*LQNA  ) = DFNSL(M,2*LQNA  ) + QA(IOCC)*DSL
              DFNSS(M,2*LQNA  ) = DFNSS(M,2*LQNA  ) + QA(IOCC)*DSS
              DFNLS(M,2*LQNA  ) = DFNLS(M,2*LQNA  ) + QA(IOCC)*DLS
C
             ENDIF
C
          ENDDO
C
        ENDDO
      ENDDO
C
C     ONE-BODY EIGENVALUE ENERGIES FOR OCCUPIED ELECTRONS
      DO IOCC=1,NUMOCC(LQNA)
        IF(NORB(LQNA,IOCC).GT.0) THEN
          EH = EH + QA(IOCC)*RK2A1*EIG1(NBLC+IOCC)
        ENDIF
      ENDDO
C
140   CONTINUE
C
C     NEGATIVE KAPPA(A) CHOICE (APPLIES TO ALL LQNA VALUES)
C
C     WORKING EIGENVECTOR MATRIX
      DO IBAS=1,NMAT
        DO JBAS=1,NMAT
          C2(IBAS,JBAS) = H2(IBAS,JBAS)
        ENDDO
      ENDDO
C
C     DIAGONALISE FOCK MATRIX (THIS NEEDS LAPACK LIBRARY)
      CALL DSYGV(1,'V','U',NMAT,C2,MBD,O2,MBD,EIG2,TE,LWK0,INFO)
      IF(INFO.NE.0) THEN
        WRITE(6, *) 'In HFSCF0: eigenvalue solver DSYGV(2) failed.',INFO
        WRITE(7, *) 'In HFSCF0: eigenvalue solver DSYGV(2) failed.',INFO
        STOP
      ENDIF
c
c      diagonalise the overlap matrix to inspect linear dependence      
c      CALL DSYEV('V','U',NMAT,O2,MBD,EIG2,TE,LWK0,INFO)
c      IF(INFO.NE.0) THEN
c        WRITE(6, *) 'In HFSCF0: eigenvalue solver DSYEV(2) failed.',INFO
c        WRITE(7, *) 'In HFSCF0: eigenvalue solver DSYEV(2) failed.',INFO
c        STOP
c      ENDIF
c      write(*,*) '----'
c      write(*,*) KAPA2
c      do n=1,NBLC
c        write(*,333) n,EIG2(n),EIG2(n+nblc)
c      enddo
C
C     ATOMIC SELECTION RULE: ORTHOGONALITY IN BLOCKS OF KQN -> KA = KB
C
C     BEGIN LOOP OVER ALL MQNA VALUES
      DO IMVAL=1,IABS(KAPA2)
C
C       COEFFICIENT MATRIX ADDRESSES
        IL1 = LRGE(IZ,2*LQNA+1,2*IMVAL-1)
        IL2 = LRGE(IZ,2*LQNA+1,2*IMVAL  )
        IS1 = LRGE(IZ,2*LQNA+1,2*IMVAL-1)+NSKP
        IS2 = LRGE(IZ,2*LQNA+1,2*IMVAL  )+NSKP
C
C       LOOP OVER OCCUPIED SOLUTIONS OF THE EIGENVALUE PROBLEM
c       DO IOCC=1,NUMOCC(LQNA)
C       LOOP OVER ALL SOLUTIONS OF THE EIGENVALUE PROBLEM
        DO IOCC=1,NBASA
C
C         COPY ENERGY EIGENVALUES TO MASTER LIST
          EIGN(NSKP+IL1+IOCC) = EIG2(NBLC+IOCC)
          EIGN(NSKP+IL2+IOCC) = EIG2(NBLC+IOCC)
          EIGN(     IL1+IOCC) = EIG2(     IOCC)
          EIGN(     IL2+IOCC) = EIG2(     IOCC)
C
C         EFFECTIVE FRACTIONAL OCCUPANCY FOR THIS ORBITAL
          IF(IOCC.LE.NUMOCC(LQNA)) THEN
            QF = DSQRT(QA(IOCC))
          ELSE
            QF = 1.0D0
          ENDIF
C
C         LOOP OVER BASIS FUNCTIONS
          DO IBAS=1,NBASA
C
C           KRAMERS PAIR LARGE- AND SMALL-COMPONENT COEFFICIENTS
            CLP = QF*C2(     IBAS,NBLC+IOCC)
            IF(HMLT.NE.'NORL') THEN
              CSP = QF*C2(NBLC+IBAS,NBLC+IOCC)
              CLN = QF*C2(     IBAS,     IOCC)
              CSN = QF*C2(NBLC+IBAS,     IOCC)
            ENDIF
C
C           COPY INTO MASTER COEFFICIENT LIST
            COEF(IL1+IBAS,NSKP+IL1+IOCC) = DCMPLX(CLP,0.0D0)
            COEF(IL2+IBAS,NSKP+IL2+IOCC) = DCMPLX(CLP,0.0D0)
            IF(HMLT.NE.'NORL') THEN
              COEF(IS1+IBAS,NSKP+IL1+IOCC) = DCMPLX(CSP,0.0D0)
              COEF(IS2+IBAS,NSKP+IL2+IOCC) = DCMPLX(CSP,0.0D0)
              COEF(IL1+IBAS,     IL1+IOCC) = DCMPLX(CLN,0.0D0)
              COEF(IL2+IBAS,     IL2+IOCC) = DCMPLX(CLN,0.0D0)
              COEF(IS1+IBAS,     IL1+IOCC) = DCMPLX(CSN,0.0D0)
              COEF(IS2+IBAS,     IL2+IOCC) = DCMPLX(CSN,0.0D0)
            ENDIF
C          
          ENDDO
        ENDDO
C
C       LOOP OVER ALL OCCUPIED SUBSHELLS OF THIS KQN TYPE GIVEN MQN
        DO IOCC=1,NUMOCC(LQNA)
C
C         SKIP UNOCCUPIED ORBITALS
          IF(NORB(LQNA,IOCC).NE.0) THEN
C
C           FOCK ADDRESS FOR THIS PAIR
            IOCPN(IOCCML+1) = LRGE(IZ,2*LQNA+1,2*IMVAL-1)+IOCC
            IOCPN(IOCCML+2) = LRGE(IZ,2*LQNA+1,2*IMVAL  )+IOCC
C
C           SAVE THESE STARTING ADDRESSES
            IOCTP(IZ,2*LQNA+1,2*IMVAL-1,IOCC) = IOCPN(IOCCML+1)
            IOCTP(IZ,2*LQNA+1,2*IMVAL  ,IOCC) = IOCPN(IOCCML+2)
C
C           INCREASE FOCK ADDRESS OF OCCUPIED ORBITALS (PAIR AT A TIME)
            IOCCML = IOCCML+2
C
          ENDIF
C
        ENDDO
C
      ENDDO
C
C     NON-RELATIVISTIC SPECIAL CASE: ALSO FILL IN THE +KQNA BLOCK
      IF(HMLT.EQ.'NORL'.AND.LQNA.GT.0) THEN
C
        DO IMVAL=1,IABS(KAPA(2*LQNA  ,IZ))
C
C         COEFFICIENT MATRIX ADDRESSES
          IL1 = LRGE(IZ,2*LQNA  ,2*IMVAL-1)
          IL2 = LRGE(IZ,2*LQNA  ,2*IMVAL  )
C
C         LOOP OVER OCCUPIED SOLUTIONS OF THE EIGENVALUE PROBLEM
C         DO IOCC=1,NUMOCC(LQNA)
C         LOOP OVER ALL SOLUTIONS OF THE EIGENVALUE PROBLEM
          DO IOCC=1,NUMOCC(LQNA)
C
C           COPY ENERGY EIGENVALUES TO MASTER LIST
            EIGN(NSKP+IL1+IOCC) = EIG2(NBLC+IOCC)
            EIGN(NSKP+IL2+IOCC) = EIG2(NBLC+IOCC)
C
C           EFFECTIVE FRACTIONAL OCCUPANCY FOR THIS ORBITAL
            IF(IOCC.LE.NUMOCC(LQNA)) THEN
              QF = DSQRT(QA(IOCC))
            ELSE
              QF = 1.0D0
            ENDIF
C
C           LOOP OVER BASIS FUNCTIONS
            DO IBAS=1,NBASA
C
C             KRAMERS PAIR LARGE- AND SMALL-COMPONENT COEFFICIENTS
              CLP = QF*C2(     IBAS,NBLC+IOCC)
C
C             COPY INTO MASTER COEFFICIENT LIST
              COEF(IL1+IBAS,NSKP+IL1+IOCC) = DCMPLX(CLP,0.0D0)
              COEF(IL2+IBAS,NSKP+IL2+IOCC) = DCMPLX(CLP,0.0D0)
C          
            ENDDO
          ENDDO
C
C         LOOP OVER ALL OCCUPIED SUBSHELLS OF THIS KQN TYPE GIVEN MQN
          DO IOCC=1,NUMOCC(LQNA)
C
C           SKIP UNOCCUPIED ORBITALS
            IF(NORB(LQNA,IOCC).NE.0) THEN
C
C             FOCK ADDRESS FOR THIS PAIR
              IOCPN(IOCCML+1) = LRGE(IZ,2*LQNA  ,2*IMVAL-1)+IOCC
              IOCPN(IOCCML+2) = LRGE(IZ,2*LQNA  ,2*IMVAL  )+IOCC
C
C             SAVE THESE STARTING ADDRESSES
              IOCTP(IZ,2*LQNA  ,2*IMVAL-1,IOCC) = IOCPN(IOCCML+1)
              IOCTP(IZ,2*LQNA  ,2*IMVAL  ,IOCC) = IOCPN(IOCCML+2)
C
C             INCREASE FOCK ADDRESS OF OCCUPIED ORBITALS (PAIR AT A TIME)
              IOCCML = IOCCML+2
C
            ENDIF
C
          ENDDO
C
        ENDDO
C
      ENDIF
C
C     BUILD ATOMIC CHARGE DENSITY MATRIX FOR THIS KQNA BLOCK
      M = 0
      DO IBAS=1,NBASA
        DO JBAS=1,NBASA
          M = M+1
C
C         INITIALISE ATOMIC DENSITY LISTS FOR THIS BLOCK
          DENLL(M,2*LQNA+1) = 0.0D0
          DFNLL(M,2*LQNA+1) = 0.0D0
C
          IF(HMLT.NE.'NORL') THEN
            DENSL(M,2*LQNA+1) = 0.0D0
            DFNSL(M,2*LQNA+1) = 0.0D0
C
            DENSS(M,2*LQNA+1) = 0.0D0
            DFNSS(M,2*LQNA+1) = 0.0D0
C
            DENLS(M,2*LQNA+1) = 0.0D0
            DFNLS(M,2*LQNA+1) = 0.0D0
          ENDIF
C
C         LOOP OVER ALL LISTED SUBSHELLS OF THIS KQN TYPE
          DO IOCC=1,NUMOCC(LQNA)
C
C           SKIP THIS STEP IF SUBSHELL IS UNOCCUPIED
            IF(NORB(LQNA,IOCC).GT.0) THEN
C
C             LL DENSITY CONTRIBUTIONS
              DLL = C2(IBAS     ,NBLC+IOCC)*C2(JBAS     ,NBLC+IOCC)
              DENLL(M,2*LQNA+1) = DENLL(M,2*LQNA+1) + QE(IOCC)*DLL
              DFNLL(M,2*LQNA+1) = DFNLL(M,2*LQNA+1) + QA(IOCC)*DLL
C
              IF(HMLT.NE.'NORL') THEN
C
C               SL DENSITY CONTRIBUTIONS
                DSL = C2(IBAS+NBLC,NBLC+IOCC)*C2(JBAS     ,NBLC+IOCC)
                DENSL(M,2*LQNA+1) = DENSL(M,2*LQNA+1) + QE(IOCC)*DSL
                DFNSL(M,2*LQNA+1) = DFNSL(M,2*LQNA+1) + QA(IOCC)*DSL
C
C               SS DENSITY CONTRIBUTIONS
                DSS = C2(IBAS+NBLC,NBLC+IOCC)*C2(JBAS+NBLC,NBLC+IOCC)
                DENSS(M,2*LQNA+1) = DENSS(M,2*LQNA+1) + QE(IOCC)*DSS
                DFNSS(M,2*LQNA+1) = DFNSS(M,2*LQNA+1) + QA(IOCC)*DSS
C
C               LS DENSITY CONTRIBUTIONS
                DLS = C2(IBAS     ,NBLC+IOCC)*C2(JBAS+NBLC,NBLC+IOCC)
                DENLS(M,2*LQNA+1) = DENLS(M,2*LQNA+1) + QE(IOCC)*DLS
                DFNLS(M,2*LQNA+1) = DFNLS(M,2*LQNA+1) + QA(IOCC)*DLS
C
              ENDIF
C
            ENDIF
C
          ENDDO
        ENDDO
      ENDDO
C
C     ONE-BODY EIGENVALUE ENERGIES FOR OCCUPIED ELECTRONS
      DO IOCC=1,NUMOCC(LQNA)
        IF(NORB(LQNA,IOCC).GT.0) THEN
          EH = EH + QA(IOCC)*RK2A2*EIG2(NBLC+IOCC)
        ENDIF
      ENDDO
C
C     END LOOP OVER LQNA VALUES
100   CONTINUE
C
C     COULOMB AND BREIT ENERGIES HAVE BEEN DOUBLE-COUNTED
      EG = EG/2.0D0
      EB = EB/2.0D0
C
      EG2 = EG+EB
C
C     TOTAL ATOMIC ENERGY IN THIS ITERATION
C     ENEW = EK+EZ+EE+EU+EP+EW-EG-EB
      ENEW = EH-EG-EB
C
C     RELATIVE CHANGE IN ENERGY
      ETEST = DABS((EPRV-ENEW)/ENEW)
C
C     WRITE THE ITERATION NUMBER AND THE TOTAL ENERGY
      WRITE(6,41) ITER,EH,EG2,ENEW,ETEST
      WRITE(7,41) ITER,EH,EG2,ENEW,ETEST
C
C     CONVERGENCE CRITERIA NOT RELEVANT (ONE-BODY PROBLEM)
      IF(ICRG.EQ.1.AND.HMLT.NE.'DHFQ') THEN
        GOTO 1001
      ENDIF
C
C     SUCCESSFUL CONVERGENCE
      IF(ETEST.LE.ENRGTOL) THEN
        GOTO 1001
      ELSE
        EPRV = ENEW
      ENDIF
C
C     BARE NUCLEUS APPROXIMATION
      IF(HMLT.EQ.'BARE') GOTO 1001
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
C     WRITTEN SUMMARY                                                  C
C**********************************************************************C
C
C     SKIP SUMMARY FOR TRIVIAL SOLUTIONS
      IF(ITER.EQ.1) GOTO 550
C
C     SUMMARY OF ENERGY CONTRIBUTIONS
50    FORMAT(1X,A,17X,F25.12)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',20),'Atomic energies (Hartree units)'
      WRITE(7, *) REPEAT(' ',20),'Atomic energies (Hartree units)'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,50) 'Electron kinetic energy       ',EK
      WRITE(7,50) 'Electron kinetic energy       ',EK
      WRITE(6,50) 'Nuclear attraction energy     ',EZ
      WRITE(7,50) 'Nuclear attraction energy     ',EZ
      WRITE(6,50) 'Two-electron energy (Coulomb) ',EG
      WRITE(7,50) 'Two-electron energy (Coulomb) ',EG
      IF(HMLT.NE.'DHFB'.AND.HMLT.NE.'DHFQ') GOTO 500
      WRITE(6,50) 'Two-electron energy (Breit)   ',EB
      WRITE(7,50) 'Two-electron energy (Breit)   ',EB
      IF(HMLT.NE.'DHFQ') GOTO 500
      WRITE(6,50) 'Anomalous moment              ',EA
      WRITE(7,50) 'Anomalous moment              ',EA
      WRITE(6,50) 'Self-interaction energy (low) ',EY
      WRITE(7,50) 'Self-interaction energy (low) ',EY
      WRITE(6,50) 'Self-interaction energy (high)',ES
      WRITE(7,50) 'Self-interaction energy (high)',ES
      WRITE(6,50) 'Nuclear Uehling energy (e+e-) ',EE
      WRITE(7,50) 'Nuclear Uehling energy (e+e-) ',EE
      WRITE(6,50) 'Nuclear Wichmann-Kroll energy ',EW
      WRITE(7,50) 'Nuclear Wichmann-Kroll energy ',EW
      WRITE(6,50) 'Nuclear Källén-Sabry energy   ',ER
      WRITE(7,50) 'Nuclear Källén-Sabry energy   ',ER
      WRITE(6,50) 'Nuclear Uehling energy (μ+μ-) ',EU
      WRITE(7,50) 'Nuclear Uehling energy (μ+μ-) ',EU
      WRITE(6,50) 'Nuclear Uehling energy (π+π-) ',EP
      WRITE(7,50) 'Nuclear Uehling energy (π+π-) ',EP
      WRITE(6,50) 'Electronic Uehling energy     ',ED
      WRITE(7,50) 'Electronic Uehling energy     ',ED
500   CONTINUE
      WRITE(6,50) 'Total energy                  ',ENEW
      WRITE(7,50) 'Total energy                  ',ENEW
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     CALCULATE THE SEEMINGLY POINTLESS RATIO THAT GRASP USES
      WRITE(6,50) 'Virial nuclear/kinetic ratio  ',(ENEW-EK)/EK
      WRITE(7,50) 'Virial nuclear/kinetic ratio  ',(ENEW-EK)/EK
550   CONTINUE
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
C
C**********************************************************************C
C     ORBITAL AND ENERGY COUNTERS FOR WHOLE MOLECULE                   C
C**********************************************************************C
C
C     UPDATE COUNTER FOR HIGHEST OCCUPIED ATOMIC ORBITAL
      IOCCM0 = IOCCML
C
C     ADD RESULTS FROM THIS ATOM TO MOLECULAR ENERGY
      ETOT = ETOT + ENEW
      EONE = EONE + EH
      ECLG = ECLG + EG
      EBRG = EBRG + EB
      ESLF = ESLF + ES+EY
      EANM = EANM + EA
      EUEE = EUEE + EE
      EUEU = EUEU + EU
      EUEP = EUEP + EP
      EWKR = EWKR + EW
      EKSB = EKSB + ER
      ESCF = ESCF + ED
C
      RETURN
      END
C
C
      SUBROUTINE CONFIG0(IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       CCCCCC   OOOOOO  NN    NN FFFFFFFF IIII GGGGGG   000000        C
C      CC    CC OO    OO NNN   NN FF        II GG    GG 00   000       C
C      CC       OO    OO NNNN  NN FF        II GG       00  0000       C
C      CC       OO    OO NN NN NN FFFFFF    II GG       00 00 00       C
C      CC       OO    OO NN  NNNN FF        II GG   GGG 0000  00       C
C      CC    CC OO    OO NN   NNN FF        II GG    GG 000   00       C
C       CCCCCC   OOOOOO  NN    NN FF       IIII GGGGGG   000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  CONFIG0 DETERMINES THE ELECTRONIC CONFIGURATION OF AN ATOM.         C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ IZ - ATOMIC CENTRE INDEX                                          C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*1 LLAB,QSGN
      CHARACTER*2 ELMT(120),ELNM
      CHARACTER*5 NMDL
      CHARACTER*6 CNFG
      CHARACTER*8 ZWRT,QWRT,EWRT
C
      COMMON/AUFB/NORB(0:MEL,MKP+1),NUMOCC(0:MEL),LMXCONF
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/FILL/NCNF(MCT,0:MEL,MKP+1),NLVL(MCT,0:MEL),CNFG(MCT)
      COMMON/MDLV/ELMT
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
C
C     IMPORT ATOMIC CHARGE DETAILS
      IZNV = INT(ZNUC(IZ))
      ELNM = ELMT(IZNV)
      ICRG = IQNC(IZ)
      MLQN =(NKAP(IZ)-1)/2
C
C     CONVERT IZNV AND ICRG TO STRINGS AND WRITE TITLE
      IF(IZNV.LT.10) THEN
        WRITE(ZWRT,'(A,I1)') 'Z = ',IZNV
      ELSEIF(IZNV.LT.100) THEN
        WRITE(ZWRT,'(A,I2)') 'Z = ',IZNV
      ELSE
        WRITE(ZWRT,'(A,I3)') 'Z = ',IZNV
      ENDIF
C
      IF(IZNV-ICRG.LT.10) THEN
        WRITE(QWRT,'(A,I2)') 'Q = ',IZNV-ICRG
      ELSEIF(IZNV.LT.100) THEN
        WRITE(QWRT,'(A,I3)') 'Q = ',IZNV-ICRG
      ELSE
        WRITE(QWRT,'(A,I4)') 'Q = ',IZNV-ICRG
      ENDIF
C
      IF(IZNV-ICRG.GT.0) THEN
        QSGN = '+'
      ELSEIF(IZNV-ICRG.LT.0) THEN
        QSGN = '-'
      ENDIF
C
      ICMD = IABS(IZNV-ICRG)
      IF(IZNV-ICRG.EQ.0) THEN
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
20    FORMAT(17X,'Centre',I3,':',3X,A,3X,A,3X,A)
      WRITE(6, *) REPEAT(' ',72)
      WRITE(7, *) REPEAT(' ',72)
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6,20) IZ,ZWRT,QWRT,EWRT
      WRITE(7,20) IZ,ZWRT,QWRT,EWRT
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C     READ OR GENERATE AUFBAU OCCUPANCY FOR THIS ATOM
      IF(CNFG(IZ).EQ.'MANUAL') THEN
        LMXCONF = MLQN
        DO LQN=0,LMXCONF
          NUMOCC(LQN) = NLVL(IZ,LQN)
          DO NQN=1,NLVL(IZ,LQN)
            NORB(LQN,NQN) = NCNF(IZ,LQN,NQN)
          ENDDO
        ENDDO
      ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
        CALL AUFBAU(IZNV,ICRG,NORB,NUMOCC,LMXCONF)
      ELSE
        WRITE(6, *) 'In CONFIG0: invalid configuration choice.'
        WRITE(7, *) 'In CONFIG0: invalid configuration choice.'
      ENDIF
C
C     IDENTIFY THE HIGHEST OCCUPIED SHELL
      NMAX = 1
      DO LQN=0,LMXCONF
        IF(NUMOCC(LQN).GT.NMAX) THEN
          NMAX = NUMOCC(LQN)
        ENDIF
      ENDDO
C
C     CHECK WHETHER THERE ARE SUFFICIENT BASIS FUNCTION TYPES
      IF(MLQN.LT.LMXCONF) THEN
        WRITE(6, *) 'In CONFIG0: insufficient angular types in basis.'
        WRITE(7, *) 'In CONFIG0: insufficient angular types in basis.'
        WRITE(6, *) 'MLQN = ',MLQN,' and LMXCONF = ',LMXCONF
        WRITE(7, *) 'MLQN = ',MLQN,' and LMXCONF = ',LMXCONF
        STOP
      ENDIF
C
C     PRINT ATOMIC CONFIGURATION
30    FORMAT(1X,A,2X,A,1X,'|',2X'NSHELL ',12(2X,I2))
31    FORMAT(1X,'LQN = ',I1,2X,I3,2X,'|'1X,' OCC(',A,'):',A,12(2X,I2))
C
      IF(CNFG(IZ).EQ.'MANUAL') THEN
        WRITE(6,30) 'Manual:','#fns',(N,N=1,12)
        WRITE(7,30) 'Manual:','#fns',(N,N=1,12)
      ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
        WRITE(6,30) 'Aufbau:','#fns',(N,N=1,12)
        WRITE(7,30) 'Aufbau:','#fns',(N,N=1,12)
      ENDIF
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      DO LQN=0,MLQN
        WRITE(6,31) LQN,NFNC(LQN,IZ),LLAB(LQN),
     &                   REPEAT(' ',LQN*4),(NORB(LQN,J),J=1,NUMOCC(LQN))
        WRITE(7,31) LQN,NFNC(LQN,IZ),LLAB(LQN),
     &                   REPEAT(' ',LQN*4),(NORB(LQN,J),J=1,NUMOCC(LQN))
      ENDDO
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
      WRITE(6, *) ''
      WRITE(7, *) ''
C
C     INITIALISE ELECTRON OCCUPATION COUNTER
      IOCCML = IOCCM0
C
C     LOOP OVER ALL AVAILABLE LQN VALUES
      DO 100 LQNA=0,MLQN
C
C     EFFECTIVE AND AVERAGE OCCUPATION NUMBERS FOR THIS LQNA ORBITAL
C     A CLOSED SUBSHELL (NSHELL,LQNA) CONTAINS NCLS ELECTRONS
      NCLS = 4*LQNA+2
C
C     POSITIVE KAPPA(A) CHOICE
      KAPA1 = LQNA
C
C     NEGATIVE KAPPA(A) CHOICE (APPLIES TO ALL LQNA VALUES)
      KAPA2 =-LQNA-1
C
C**********************************************************************C
C     MATRIX DIAGONALISATION AND COEFFICIENT MATRIX UPDATES            C
C**********************************************************************C
C
C     POSITIVE KAPPA(A) CHOICE (APPLIES ONLY FOR LQNA > 0)
      IF(LQNA.EQ.0.OR.HMLT.EQ.'NORL') GOTO 140
C
C     ATOMIC SELECTION RULE: ORTHOGONALITY IN BLOCKS OF KQN -> KA = KB
C
C     BEGIN LOOP OVER ALL MQNA VALUES
      DO IMVAL=1,IABS(KAPA1)
C
C       LOOP OVER ALL OCCUPIED SUBSHELLS OF THIS KQN TYPE GIVEN MQN
        DO IOCC=1,NUMOCC(LQNA)
C
C         SKIP UNOCCUPIED ORBITALS
          IF(NORB(LQNA,IOCC).NE.0) THEN
C
C           FOCK ADDRESS FOR THIS PAIR
            IOCPN(IOCCML+1) = LRGE(IZ,2*LQNA  ,2*IMVAL-1)+IOCC
            IOCPN(IOCCML+2) = LRGE(IZ,2*LQNA  ,2*IMVAL  )+IOCC
C
C           SAVE THESE STARTING ADDRESSES
            IOCTP(IZ,2*LQNA  ,2*IMVAL-1,IOCC) = IOCPN(IOCCML+1)
            IOCTP(IZ,2*LQNA  ,2*IMVAL  ,IOCC) = IOCPN(IOCCML+2)
C
C           INCREASE FOCK ADDRESS OF OCCUPIED ORBITALS (PAIR AT A TIME)
            IOCCML = IOCCML+2
C
          ENDIF
C
        ENDDO
C
      ENDDO
C
140   CONTINUE
C
C     NEGATIVE KAPPA(A) CHOICE (APPLIES TO ALL LQNA VALUES)
C
C     BEGIN LOOP OVER ALL MQNA VALUES
      DO IMVAL=1,IABS(KAPA2)
C
C       LOOP OVER ALL OCCUPIED SUBSHELLS OF THIS KQN TYPE GIVEN MQN
        DO IOCC=1,NUMOCC(LQNA)
C
C         SKIP UNOCCUPIED ORBITALS
          IF(NORB(LQNA,IOCC).NE.0) THEN
C
C           FOCK ADDRESS FOR THIS PAIR
            IOCPN(IOCCML+1) = LRGE(IZ,2*LQNA+1,2*IMVAL-1)+IOCC
            IOCPN(IOCCML+2) = LRGE(IZ,2*LQNA+1,2*IMVAL  )+IOCC
C
C           SAVE THESE STARTING ADDRESSES
            IOCTP(IZ,2*LQNA+1,2*IMVAL-1,IOCC) = IOCPN(IOCCML+1)
            IOCTP(IZ,2*LQNA+1,2*IMVAL  ,IOCC) = IOCPN(IOCCML+2)
C
C           INCREASE FOCK ADDRESS OF OCCUPIED ORBITALS (PAIR AT A TIME)
            IOCCML = IOCCML+2
C
          ENDIF
C
        ENDDO
C
      ENDDO
C
C     NON-RELATIVISTIC SPECIAL CASE: ALSO FILL IN THE +KQNA BLOCK
      IF(HMLT.EQ.'NORL'.AND.LQNA.GT.0) THEN
C
        DO IMVAL=1,IABS(KAPA(2*LQNA  ,IZ))
C
C         LOOP OVER ALL OCCUPIED SUBSHELLS OF THIS KQN TYPE GIVEN MQN
          DO IOCC=1,NUMOCC(LQNA)
C
C           SKIP UNOCCUPIED ORBITALS
            IF(NORB(LQNA,IOCC).NE.0) THEN
C
C             FOCK ADDRESS FOR THIS PAIR
              IOCPN(IOCCML+1) = LRGE(IZ,2*LQNA  ,2*IMVAL-1)+IOCC
              IOCPN(IOCCML+2) = LRGE(IZ,2*LQNA  ,2*IMVAL  )+IOCC
C
C             SAVE THESE STARTING ADDRESSES
              IOCTP(IZ,2*LQNA  ,2*IMVAL-1,IOCC) = IOCPN(IOCCML+1)
              IOCTP(IZ,2*LQNA  ,2*IMVAL  ,IOCC) = IOCPN(IOCCML+2)
C
C             INCREASE FOCK ADDRESS OF OCCUPIED ORBITALS (PAIR AT A TIME)
              IOCCML = IOCCML+2
C
            ENDIF
C
          ENDDO
C
        ENDDO
C
      ENDIF
C
C     END LOOP OVER LQNA VALUES
100   CONTINUE
C
C     UPDATE COUNTER FOR HIGHEST OCCUPIED ATOMIC ORBITAL
      IOCCM0 = IOCCML
C
      RETURN
      END
C
C
      SUBROUTINE H1INT0(IZ,H1INT,E1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             HH    HH  11  IIII NN    NN TTTTTTTT 000000              C
C             HH    HH 111   II  NNN   NN    TT   00   000             C
C             HH    HH  11   II  NNNN  NN    TT   00  0000             C
C             HHHHHHHH  11   II  NN NN NN    TT   00 00 00             C
C             HH    HH  11   II  NN  NNNN    TT   0000  00             C
C             HH    HH  11   II  NN   NNN    TT   000   00             C
C             HH    HH 1111 IIII NN    NN    TT    000000              C
C                                                                      C
C -------------------------------------------------------------------- C
C  H1INT0 CALCULATES THE ATOMIC PERTURBATIVE EFFECTS ARISING FROM AN   C
C  INTERACTION, H1INT. THE RESULTS ARE GIVEN BY ELECTRONIC ORBITAL.    C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ IZ    - ATOMIC CENTRE INDEX.                                      C
C  ▶ H1INT - ONE-BODY INTERACTION HAMILTONIAN.                         C
C  OUTPUT:                                                             C
C  ▶ E1    - TOTAL ATOMIC ENERGY UNDER THIS INTERACTION HAMILTONIAN.   C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*1 LLAB,CHS
      CHARACTER*5 NMDL
      CHARACTER*6 CNFG,H1INT
C
      COMPLEX*16 COEF(MDM,MDM)
C
      DIMENSION QA(MKP)
      DIMENSION PLTN(MCT,MKP+1,MKP),IGLB(MCT,MKP+1,MKP)
      DIMENSION V1(MBD,MBD)
      
      DIMENSION EXL(MBS)
      DIMENSION RGLB(MCT,MKP+1,MKP,5)
C
      COMMON/AUFB/NORB(0:MEL,MKP+1),NUMOCC(0:MEL),LMXCONF
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/FILL/NCNF(MCT,0:MEL,MKP+1),NLVL(MCT,0:MEL),CNFG(MCT)
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
      COMMON/SHLL/ACFF,BCFF,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
C
C     INITIALISE TOTAL ENERGIES
      ELLT = 0.0D0
      ELST = 0.0D0
      ESLT = 0.0D0
      ESST = 0.0D0
C
C**********************************************************************C
C     SORTING SECTION: FIND ADDRESSES FOR USUAL ATOMIC ORDER           C
C**********************************************************************C
C
C     COUNTER FOR OCCUPIED ORBITALS
      IOCCML = 0
C
C     IMPORT ATOMIC CHARGE DETAILS
      IZNV = INT(ZNUC(IZ))
      ICRG = IQNC(IZ)
      MLQN = (NKAP(IZ)-1)/2
C
C     INITIALISE THE ORBITAL FILLING ARRAY
      DO LQN=0,MEL
        DO NQN=1,MKP+1
          NORB(LQN,NQN) = 0
        ENDDO
      ENDDO
C
C     READ OR GENERATE AUFBAU OCCUPANCY FOR THIS ATOM
      IF(CNFG(IZ).EQ.'MANUAL') THEN
        LMXCONF = (NKAP(IZ)-1)/2
        DO LQN=0,LMXCONF
          NUMOCC(LQN) = NLVL(IZ,LQN)
          DO NQN=1,NLVL(IZ,LQN)
            NORB(LQN,NQN) = NCNF(IZ,LQN,NQN)
          ENDDO
        ENDDO
      ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
        CALL AUFBAU(IZNV,ICRG,NORB,NUMOCC,LMXCONF)
      ELSE
        WRITE(6, *) 'In H1INT0: invalid configuration choice.'
        WRITE(7, *) 'In H1INT0: invalid configuration choice.'
      ENDIF
C
C     LOOP OVER LQNS
      DO LQN=0,LMXCONF
C
C       RECORD LQNA VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
        NBAS = NFNC(LQN,IZ)
        DO IBAS=1,NBAS
          EXL(IBAS) = BEXL(IBAS,LQN,IZ)
        ENDDO
C
C       FRACTIONAL SUBSHELL OCCUPANCY
        DO NQN=1,NUMOCC(LQN)
          QA(NQN) = DFLOAT(NORB(LQN,NQN))/DFLOAT(4*LQN+2)
          IF(NORB(LQN,NQN).EQ.1) THEN
            QA(NQN) = 1.0D0
          ENDIF
        ENDDO
C
C       POSITIVE KQN CHOICE (APPLIES ONLY FOR LQN>0)
        IF(LQN.EQ.0) GOTO 110
        KQN = LQN
C
C       DEGENERACY
        IDGN = 2*LQN
C
C       SPECIAL CASE: ONE-ELECTRON SYSTEM
        IF(NOCC.EQ.1) THEN
          DO NQN=1,NUMOCC(LQN)
            IDGN    = 1
            QA(NQN) = 1.0D0
          ENDDO
        ENDIF
C
C       CALCULATE BATCH OF MATRIX ELEMENT INTEGRALS
        IF(H1INT.EQ.'ANOMLS') THEN
          CALL ANMLS0(V1,EXL,IZ,KQN,NBAS)
        ELSEIF(H1INT.EQ.'SLFLW0') THEN
          CALL SLFLW0(V1,    IZ,KQN     )
        ELSEIF(H1INT.EQ.'SLFHI0') THEN
          CALL SLFHI0(V1,EXL,IZ,KQN,NBAS)
        ELSEIF(H1INT.EQ.'UEENUC') THEN
          CALL UEHNC0(V1,EXL,IZ,KQN,NBAS,'ELEC')
        ELSEIF(H1INT.EQ.'UEUNUC') THEN
          CALL UEHNC0(V1,EXL,IZ,KQN,NBAS,'MUON')
        ELSEIF(H1INT.EQ.'UEPNUC') THEN
          CALL UEHNC0(V1,EXL,IZ,KQN,NBAS,'PION')
        ELSEIF(H1INT.EQ.'WKRNUC') THEN
          CALL WKRNC0(V1,EXL,IZ,KQN,NBAS)
        ELSEIF(H1INT.EQ.'KSBNUC') THEN
          CALL KSBNC0(V1,EXL,IZ,KQN,NBAS)
        ELSEIF(H1INT.EQ.'UEHSCF') THEN
          CALL UEHEL0(V1,EXL,IZ,KQN,NBAS)
        ENDIF
C
C       LOOP OVER |MQN| CHOICES FOR THIS KQN
        DO MQN=1,LQN
C
C         LOOP OVER LISTED NQNS
          DO NQN=1,NUMOCC(LQN)
C
C           SKIP THIS STEP IF ORBITAL IS UNOCCUPIED
            IF(NORB(LQN,NQN).EQ.0) GOTO 101
C
C           SPHERICAL SYMMETRY: FOCUS ON FIRST MQN ONLY
            IF(MQN.EQ.1) THEN
C
C             ADDRESS AND AVERAGE OCCUPATION
              IGLB(IZ,NQN+LQN,2*LQN  ) = IOCPN(IOCCML+1)
              PLTN(IZ,NQN+LQN,2*LQN  ) = QA(NQN)
C
C             FOCK MATRIX ADDRESSES
              IL = LRGE(IZ,2*LQN  ,MQN)
              IS = LRGE(IZ,2*LQN  ,MQN)+NSKP
C
C             COEFFICIENT MATRIX OFFSET
              MOCC = IOCTP(IZ,2*LQN  ,MQN,NQN)+NSKP
C
C             SPECIAL OPTION FOR ENERGY-DEPENDENT OPERATOR
              IF(H1INT.EQ.'SLFLWN') THEN
                CALL SLFLWN(V1,IZ,KQN,MOCC)
              ENDIF
C
C             RADIAL MOMENTS
              ELL = 0.0D0
              ELS = 0.0D0
              ESL = 0.0D0
              ESS = 0.0D0
              DO IBAS=1,NBAS
                DO JBAS=1,NBAS
C
C                 LARGE-COMPONENT OFFSETS
                  IA = IL+IBAS
                  JA = IL+JBAS
C
C                 SMALL-COMPONENT OFFSETS
                  IB = IS+IBAS
                  JB = IS+JBAS
C
                  KBAS = IBAS+NBAS
                  LBAS = JBAS+NBAS
                    
                  DAA = DREAL(DCONJG(COEF(IA,MOCC))*COEF(JA,MOCC))
                  DAB = DREAL(DCONJG(COEF(IA,MOCC))*COEF(JB,MOCC))
                  DBA = DREAL(DCONJG(COEF(IB,MOCC))*COEF(JA,MOCC))
                  DBB = DREAL(DCONJG(COEF(IB,MOCC))*COEF(JB,MOCC))
C
                  ELL = ELL + DAA*V1(IBAS,JBAS)
                  ELS = ELS + DAB*V1(IBAS,LBAS)
                  ESL = ESL + DBA*V1(KBAS,JBAS)
                  ESS = ESS + DBB*V1(KBAS,LBAS)
C
                ENDDO
              ENDDO
              ETOT = ELL+ELS+ESL+ESS
              RGLB(IZ,NQN+LQN,2*LQN  ,1) = ELL/QA(NQN)
              RGLB(IZ,NQN+LQN,2*LQN  ,2) = ELS/QA(NQN)
              RGLB(IZ,NQN+LQN,2*LQN  ,3) = ESL/QA(NQN)
              RGLB(IZ,NQN+LQN,2*LQN  ,4) = ESS/QA(NQN)
              RGLB(IZ,NQN+LQN,2*LQN  ,5) = ETOT/QA(NQN)
              
              ELLT = ELLT + IDGN*ELL
              ELST = ELST + IDGN*ELS
              ESLT = ESLT + IDGN*ESL
              ESST = ESST + IDGN*ESS
C
            ENDIF
C
C           UPDATE NUMBER OF OCCUPIED PAIRS
            IOCCML = IOCCML+2
C
C           SKIP STEP FOR UNOCCUPIED ORBITALS
101         CONTINUE
C
C         END LOOP OVER LISTED NQNS
          ENDDO
C
        ENDDO
C
110     CONTINUE
C
C       NEGATIVE KQN CHOICE
        KQN =-LQN-1
C
C       DEGENERACY
        IDGN = 2*LQN+2
C
C       SPECIAL CASE: ONE-ELECTRON SYSTEM
        IF(NOCC.EQ.1) THEN
          DO NQN=1,NUMOCC(LQN)
            IDGN    = 1
            QA(NQN) = 1.0D0
          ENDDO
        ENDIF
C
C       CALCULATE BATCH OF MATRIX ELEMENT INTEGRALS
        IF(H1INT.EQ.'ANOMLS') THEN
          CALL ANMLS0(V1,EXL,IZ,KQN,NBAS)
        ELSEIF(H1INT.EQ.'SLFLW0') THEN
          CALL SLFLW0(V1,    IZ,KQN     )
        ELSEIF(H1INT.EQ.'SLFHI0') THEN
          CALL SLFHI0(V1,EXL,IZ,KQN,NBAS)
        ELSEIF(H1INT.EQ.'UEENUC') THEN
          CALL UEHNC0(V1,EXL,IZ,KQN,NBAS,'ELEC')
        ELSEIF(H1INT.EQ.'UEUNUC') THEN
          CALL UEHNC0(V1,EXL,IZ,KQN,NBAS,'MUON')
        ELSEIF(H1INT.EQ.'UEPNUC') THEN
          CALL UEHNC0(V1,EXL,IZ,KQN,NBAS,'PION')
        ELSEIF(H1INT.EQ.'WKRNUC') THEN
          CALL WKRNC0(V1,EXL,IZ,KQN,NBAS)
        ELSEIF(H1INT.EQ.'KSBNUC') THEN
          CALL KSBNC0(V1,EXL,IZ,KQN,NBAS)
        ELSEIF(H1INT.EQ.'UEHSCF') THEN
          CALL UEHEL0(V1,EXL,IZ,KQN,NBAS)
        ENDIF
C
C       LOOP OVER |MQN| CHOICES
        DO MQN=1,LQN+1
C
C         LOOP OVER LISTED NQNS
          DO NQN=1,NUMOCC(LQN)
C
C           SKIP THIS STEP IF ORBITAL IS UNOCCUPIED
            IF(NORB(LQN,NQN).EQ.0) GOTO 102
C
C           SPHERICAL SYMMETRY: FOCUS ON FIRST MQN ONLY
            IF(MQN.EQ.1) THEN
C
C             ADDRESS AND AVERAGE OCCUPATION
              IGLB(IZ,NQN+LQN,2*LQN+1) = IOCPN(IOCCML+1)
              PLTN(IZ,NQN+LQN,2*LQN+1) = QA(NQN)
C
C             FOCK MATRIX ADDRESSES
              IL = LRGE(IZ,2*LQN+1,MQN)
              IS = LRGE(IZ,2*LQN+1,MQN)+NSKP
C
C             COEFFICIENT MATRIX OFFSET
              MOCC = IOCTP(IZ,2*LQN+1,MQN,NQN)+NSKP
C
C             SPECIAL OPTION FOR ENERGY-DEPENDENT OPERATOR
              IF(H1INT.EQ.'SLFLWN') THEN
                CALL SLFLWN(V1,IZ,KQN,MOCC)
              ENDIF
C
C             RADIAL MOMENTS
              ELL = 0.0D0
              ELS = 0.0D0
              ESL = 0.0D0
              ESS = 0.0D0
              DO IBAS=1,NBAS
                DO JBAS=1,NBAS
C
C                 LARGE-COMPONENT OFFSETS
                  IA = IL+IBAS
                  JA = IL+JBAS
C
C                 SMALL-COMPONENT OFFSETS
                  IB = IS+IBAS
                  JB = IS+JBAS
C
                  KBAS = IBAS+NBAS
                  LBAS = JBAS+NBAS
C
                  DAA = DREAL(DCONJG(COEF(IA,MOCC))*COEF(JA,MOCC))
                  DAB = DREAL(DCONJG(COEF(IA,MOCC))*COEF(JB,MOCC))
                  DBA = DREAL(DCONJG(COEF(IB,MOCC))*COEF(JA,MOCC))
                  DBB = DREAL(DCONJG(COEF(IB,MOCC))*COEF(JB,MOCC))
C
                  ELL = ELL + DAA*V1(IBAS,JBAS)
                  ELS = ELS + DAB*V1(IBAS,LBAS)
                  ESL = ESL + DBA*V1(KBAS,JBAS)
                  ESS = ESS + DBB*V1(KBAS,LBAS)
C
                ENDDO
              ENDDO

              ETOT = ELL+ELS+ESL+ESS
              RGLB(IZ,NQN+LQN,2*LQN+1,1) = ELL/QA(NQN)
              RGLB(IZ,NQN+LQN,2*LQN+1,2) = ELS/QA(NQN)
              RGLB(IZ,NQN+LQN,2*LQN+1,3) = ESL/QA(NQN)
              RGLB(IZ,NQN+LQN,2*LQN+1,4) = ESS/QA(NQN)
              RGLB(IZ,NQN+LQN,2*LQN+1,5) = ETOT/QA(NQN)
C
              ELLT = ELLT + IDGN*ELL
              ELST = ELST + IDGN*ELS
              ESLT = ESLT + IDGN*ESL
              ESST = ESST + IDGN*ESS
C
            ENDIF
C
C           UPDATE NUMBER OF OCCUPIED PAIRS
            IOCCML = IOCCML+2
C
C           SKIP STEP FOR UNOCCUPIED ORBITALS
102         CONTINUE
C
C         END LOOP OVER LISTED NQNS
          ENDDO
C
        ENDDO
C
C     END LOOP OVER LQNS
      ENDDO
C
C**********************************************************************C
C     PRINT SPECTRUM                                                   C
C**********************************************************************C
C
C     PRINT TITLE TO TERMINAL/FILE
20    FORMAT(1X,' nl  |',11X,'E(LL)',11X,
     &                 'E(LS)',11X,'E(SL)',11X,'E(SS) |',9X,'E(TOT)')
21    FORMAT(1X,I2,A,A,' |',4(1X,F15.10),' |',F15.10)
22    FORMAT(25X,'Interaction ',A,' by atomic orbitals')
23    FORMAT(31X,'First order in RSPT: ',A)
24    FORMAT(1X,A,' |',4(1X,F15.10),' |',F15.10)
C
      WRITE(6,22) H1INT
      WRITE(7,22) H1INT
      WRITE(6,23) 'E(1)'
      WRITE(7,23) 'E(1)'
      WRITE(6, *) REPEAT('=',87)
      WRITE(7, *) REPEAT('=',87)
      WRITE(6,20)
      WRITE(7,20)
      WRITE(6, *) REPEAT('-',87)
      WRITE(7, *) REPEAT('-',87)
C
C     IMPORT ATOMIC CHARGE DETAILS
      IZNV = INT(ZNUC(IZ))
      ICRG = IQNC(IZ)
C
C     INITIALISE THE ORBITAL FILLING ARRAY
      DO LQN=0,MEL
        DO NQN=1,MKP+1
          NORB(LQN,NQN) = 0
        ENDDO
      ENDDO
C
C     READ OR GENERATE AUFBAU OCCUPANCY FOR THIS ATOM
      IF(CNFG(IZ).EQ.'MANUAL') THEN
        LMXCONF = (NKAP(IZ)-1)/2
        DO L=0,LMXCONF
          NUMOCC(L) = NLVL(IZ,L)
          DO N=1,NLVL(IZ,L)
            NORB(L,N) = NCNF(IZ,L,N)
          ENDDO
        ENDDO
      ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
        CALL AUFBAU(IZNV,ICRG,NORB,NUMOCC,LMXCONF)
      ELSE
        WRITE(6, *) 'In H1INT0: invalid configuration choice.'
        WRITE(7, *) 'In H1INT0: invalid configuration choice.'
      ENDIF
C
C     LOOP OVER ORBITALS IN THE USUAL AUFBAU ORDER
      DO NQN=1,MKP
        DO LQN=0,NQN-1
C
          IF(LQN.GT.MEL) GOTO 150
          IF(NORB(LQN,NQN-LQN).EQ.0) GOTO 150
C
C         POSITIVE KQN CHOICE (APPLIES ONLY FOR LQN>0)
          IF(LQN.EQ.0) GOTO 120
C
C         SYMMETRY LABEL
          CHS = '-'
C
C         DEGENERACY
          IDGN = 2*LQN
C
C         IDENTIFY THE RIGHT INDEX AND PROPERTIES
          MOCC = IGLB(IZ,NQN,2*LQN  )+NSKP
          FRAC = PLTN(IZ,NQN,2*LQN  )
          RAVA = RGLB(IZ,NQN,2*LQN  ,1)
          RAVB = RGLB(IZ,NQN,2*LQN  ,2)
          RAVC = RGLB(IZ,NQN,2*LQN  ,3)
          RAVD = RGLB(IZ,NQN,2*LQN  ,4)
          RAVE = RGLB(IZ,NQN,2*LQN  ,5)
C
C         OUTPUT TO TERMINAL
          WRITE(6,21) NQN,LLAB(LQN),CHS,RAVA,RAVB,RAVC,RAVD,RAVE
          WRITE(7,21) NQN,LLAB(LQN),CHS,RAVA,RAVB,RAVC,RAVD,RAVE
C
120       CONTINUE
C
C         NEGATIVE KQN CHOICE (APPLIES ONLY FOR LQN>0)
C
C         SYMMETRY LABEL IF NEEDED
          IF(LQN.EQ.0) THEN
            CHS = ' '
          ELSE
            CHS = '+'
          ENDIF
C
C         DEGENERACY
          IDGN = 2*LQN+2
C
C         IDENTIFY THE RIGHT INDEX AND PROPERTIES
          MOCC = IGLB(IZ,NQN,2*LQN+1)+NSKP
          FRAC = PLTN(IZ,NQN,2*LQN+1)
          RAVA = RGLB(IZ,NQN,2*LQN+1,1)
          RAVB = RGLB(IZ,NQN,2*LQN+1,2)
          RAVC = RGLB(IZ,NQN,2*LQN+1,3)
          RAVD = RGLB(IZ,NQN,2*LQN+1,4)
          RAVE = RGLB(IZ,NQN,2*LQN+1,5)
C
C         OUTPUT TO TERMINAL
          WRITE(6,21) NQN,LLAB(LQN),CHS,RAVA,RAVB,RAVC,RAVD,RAVE
          WRITE(7,21) NQN,LLAB(LQN),CHS,RAVA,RAVB,RAVC,RAVD,RAVE
C
150       CONTINUE
C
        ENDDO
C
      ENDDO
C
      ETTT = ELLT+ELST+ESLT+ESST
      WRITE(6, *) REPEAT('-',87)
      WRITE(7, *) REPEAT('-',87)
      WRITE(6,24) 'Atom',ELLT,ELST,ESLT,ESST,ETTT
      WRITE(7,24) 'Atom',ELLT,ELST,ESLT,ESST,ETTT
      WRITE(6, *) REPEAT('=',87)
      WRITE(7, *) REPEAT('=',87)
      WRITE(6, *) ''
      WRITE(7, *) ''
C
C     ADD TO THE TALLIED TOTAL ATOMIC ENERGY
      E1 = E1+ETTT
C
      RETURN
      END
C
C
      SUBROUTINE H2INT0(IZ,H2INT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           HH    HH  222222  IIII NN    NN TTTTTTTT 000000            c
C           HH    HH 22    22  II  NNN   NN    TT   00   000           C
C           HH    HH       22  II  NNNN  NN    TT   00  0000           C
C           HHHHHHHH      22   II  NN NN NN    TT   00 00 00           C
C           HH    HH    22     II  NN  NNNN    TT   0000  00           C
C           HH    HH  22       II  NN   NNN    TT   000   00           C
C           HH    HH 22222222 IIII NN    NN    TT    000000            C
C                                                                      C
C -------------------------------------------------------------------- C
C  H2INT0 CALCULATES THE ATOMIC PERTURBATIVE EFFECTS ARISING FROM AN   C
C  INTERACTION, H2INT, AT SECOND ORDER IN RAYLEIGH-SCHRODINGER         C
C  PERTURBATION THEORY. RESULTS ARE GIVEN BY ELECTRONIC ORBITAL, AND   C
C  DECOMPOSED INTO NEGATIVE- AND POSITIVE-ENERGY CONTRIBUTIONS.        C
C -------------------------------------------------------------------- C
C  SO FAR THIS ROUTINE ONLY WORKS FOR OPERATORS THAT HAVE A SIMPLE     C
C  SPIN-ANGULAR STRUCTURE, TO TAKE ADVANTAGE OF SELECTION RULES.       C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ IZ    - ATOMIC CENTRE INDEX.                                      C
C  ▶ H2INT - ONE-BODY INTERACTION HAMILTONIAN.                         C
C  OUTPUT:                                                             C
C  ▶ E1    - TOTAL ATOMIC ENERGY UNDER THIS INTERACTION HAMILTONIAN.   C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*1 LLAB,CHS
      CHARACTER*5 NMDL
      CHARACTER*6 CNFG,H2INT
C
      COMPLEX*16 COEF(MDM,MDM)
C
      DIMENSION QA(MKP)
      DIMENSION PLTN(MKP+1,MKP),IGLB(MKP+1,MKP)
      DIMENSION V1(MBD,MBD)
C
      DIMENSION EXL(MBS)
      DIMENSION RGLB(MKP+1,MKP,5)
      DIMENSION SGLB(MKP+1,MKP,MBD,5)
C
      COMMON/AUFB/NORB(0:MEL,MKP+1),NUMOCC(0:MEL),LMXCONF
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/FILL/NCNF(MCT,0:MEL,MKP+1),NLVL(MCT,0:MEL),CNFG(MCT)
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
      COMMON/SHLL/ACFF,BCFF,FOPN,ICLS(MDM),IOPN(MDM),NCLS,NOPN,NOELEC
C
C**********************************************************************C
C     SORTING SECTION: FIND ADDRESSES FOR USUAL ATOMIC ORDER           C
C**********************************************************************C
C
C     COUNTER FOR OCCUPIED ORBITALS
      IOCCML = 0
C
C     IMPORT ATOMIC CHARGE DETAILS
      IZNV = INT(ZNUC(IZ))
      ICRG = IQNC(IZ)
      MLQN = (NKAP(IZ)-1)/2
C
C     INITIALISE THE ORBITAL FILLING ARRAY
      DO LQN=0,MEL
        DO NQN=1,MKP+1
          NORB(LQN,NQN) = 0
        ENDDO
      ENDDO
C
C     READ OR GENERATE AUFBAU OCCUPANCY FOR THIS ATOM
      IF(CNFG(IZ).EQ.'MANUAL') THEN
        LMXCONF = (NKAP(IZ)-1)/2
        DO LQN=0,LMXCONF
          NUMOCC(LQN) = NLVL(IZ,LQN)
          DO NQN=1,NLVL(IZ,LQN)
            NORB(LQN,NQN) = NCNF(IZ,LQN,NQN)
          ENDDO
        ENDDO
      ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
        CALL AUFBAU(IZNV,ICRG,NORB,NUMOCC,LMXCONF)
      ELSE
        WRITE(6, *) 'In H2INT0: invalid configuration choice.'
        WRITE(7, *) 'In H2INT0: invalid configuration choice.'
      ENDIF
C
C     LOOP OVER LQNS
      DO LQN=0,LMXCONF
C
C       RECORD LQNA VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
        NBAS = NFNC(LQN,IZ)
        DO IBAS=1,NBAS
          EXL(IBAS) = BEXL(IBAS,LQN,IZ)
        ENDDO
C
C       FRACTIONAL SUBSHELL OCCUPANCY
        DO NQN=1,NUMOCC(LQN)
          QA(NQN) = DFLOAT(NORB(LQN,NQN))/DFLOAT(4*LQN+2)
          IF(NORB(LQN,NQN).EQ.1) THEN
            QA(NQN) = 1.0D0
          ENDIF
        ENDDO
C
C       POSITIVE KQN CHOICE (APPLIES ONLY FOR LQN>0)
        IF(LQN.EQ.0) GOTO 110
        KQN = LQN
C
C       DEGENERACY
        IDGN = 2*LQN
C
C       SPECIAL CASE: ONE-ELECTRON SYSTEM
        IF(NOCC.EQ.1) THEN
          DO NQN=1,NUMOCC(LQN)
            IDGN    = 1
            QA(NQN) = 1.0D0
          ENDDO
        ENDIF
C
C       CALCULATE BATCH OF MATRIX ELEMENT INTEGRALS
        IF(H2INT.EQ.'SLFHI0') THEN
          CALL SLFHI0(V1,EXL,IZ,KQN,NBAS)
        ELSEIF(H2INT.EQ.'UEENUC') THEN
          CALL UEHNC0(V1,EXL,IZ,KQN,NBAS,'ELEC')
        ELSEIF(H2INT.EQ.'UEUNUC') THEN
          CALL UEHNC0(V1,EXL,IZ,KQN,NBAS,'MUON')
        ELSEIF(H2INT.EQ.'UEPNUC') THEN
          CALL UEHNC0(V1,EXL,IZ,KQN,NBAS,'PION')
        ELSEIF(H2INT.EQ.'WKRNUC') THEN
          CALL WKRNC0(V1,EXL,IZ,KQN,NBAS)
        ELSEIF(H2INT.EQ.'KSBNUC') THEN
          CALL KSBNC0(V1,EXL,IZ,KQN,NBAS)
        ELSEIF(H2INT.EQ.'UEHSCF') THEN
          CALL UEHEL0(V1,EXL,IZ,KQN,NBAS)
        ENDIF
C
C       LOOP OVER |MQN| CHOICES FOR THIS KQN
        DO MQN=1,LQN
C
C         LOOP OVER LISTED NQNS
          DO NQN=1,NUMOCC(LQN)
C
C           SKIP THIS STEP IF ORBITAL IS UNOCCUPIED
            IF(NORB(LQN,NQN).EQ.0) GOTO 101
C
C           SPHERICAL SYMMETRY: FOCUS ON FIRST MQN ONLY
            IF(MQN.EQ.1) THEN
C
C             ADDRESS AND AVERAGE OCCUPATION
              IGLB(NQN+LQN,2*LQN  ) = IOCPN(IOCCML+1)
              PLTN(NQN+LQN,2*LQN  ) = QA(NQN)
C
C             FOCK MATRIX ADDRESSES
              IL = LRGE(IZ,2*LQN  ,MQN)
              IS = LRGE(IZ,2*LQN  ,MQN)+NSKP
C
C             COEFFICIENT MATRIX OFFSET
              MOCC = IOCTP(IZ,2*LQN  ,MQN,NQN)+NSKP
C
C             LOOP OVER ALL STATES OF THIS TYPE
              DO N=1,NBAS
C
C               RADIAL MOMENTS
                ELL = 0.0D0
                ELS = 0.0D0
                ESL = 0.0D0
                ESS = 0.0D0
                WLL = 0.0D0
                WLS = 0.0D0
                WSL = 0.0D0
                WSS = 0.0D0
                DO IBAS=1,NBAS
                  DO JBAS=1,NBAS
C
C                   LARGE-COMPONENT OFFSETS
                    IA = IL+IBAS
                    JA = IL+JBAS
C
C                   SMALL-COMPONENT OFFSETS
                    IB = IS+IBAS
                    JB = IS+JBAS
C
C                   COEFFICIENT MATRIX OFFSET
                    NOCC = IOCTP(IZ,2*LQN  ,MQN,1)-1+N
C
                    DAA = DREAL(DCONJG(COEF(IA,MOCC))*COEF(JA,NOCC))
                    DAB = DREAL(DCONJG(COEF(IA,MOCC))*COEF(JB,NOCC))
                    DBA = DREAL(DCONJG(COEF(IB,MOCC))*COEF(JA,NOCC))
                    DBB = DREAL(DCONJG(COEF(IB,MOCC))*COEF(JB,NOCC))
C
                    KBAS = IBAS+NBAS
                    LBAS = JBAS+NBAS
C
                    ELL = ELL + DAA*V1(IBAS,JBAS)
                    ELS = ELS + DAB*V1(IBAS,LBAS)
                    ESL = ESL + DBA*V1(KBAS,JBAS)
                    ESS = ESS + DBB*V1(KBAS,LBAS)
C
C                   COEFFICIENT MATRIX OFFSET
                    NOCC = IOCTP(IZ,2*LQN  ,MQN,1)-1+NSKP+N
C
                    DAA = DREAL(DCONJG(COEF(IA,MOCC))*COEF(JA,NOCC))
                    DAB = DREAL(DCONJG(COEF(IA,MOCC))*COEF(JB,NOCC))
                    DBA = DREAL(DCONJG(COEF(IB,MOCC))*COEF(JA,NOCC))
                    DBB = DREAL(DCONJG(COEF(IB,MOCC))*COEF(JB,NOCC))
C
                    KBAS = IBAS+NBAS
                    LBAS = JBAS+NBAS
C
                    WLL = WLL + DAA*V1(IBAS,JBAS)
                    WLS = WLS + DAB*V1(IBAS,LBAS)
                    WSL = WSL + DBA*V1(KBAS,JBAS)
                    WSS = WSS + DBB*V1(KBAS,LBAS)
C
                  ENDDO
                ENDDO
                ETOT = ELL+ELS+ESL+ESS
                WTOT = WLL+WLS+WSL+WSS
                RGLB(NQN+LQN,2*LQN  ,1) = ELL/QA(NQN)
                RGLB(NQN+LQN,2*LQN  ,2) = ELS/QA(NQN)
                RGLB(NQN+LQN,2*LQN  ,3) = ESL/QA(NQN)
                RGLB(NQN+LQN,2*LQN  ,4) = ESS/QA(NQN)
                RGLB(NQN+LQN,2*LQN  ,5) = ETOT/QA(NQN)
                SGLB(NQN+LQN,2*LQN  ,N     ,1) = ELL/QA(NQN)
                SGLB(NQN+LQN,2*LQN  ,N     ,2) = ELS/QA(NQN)
                SGLB(NQN+LQN,2*LQN  ,N     ,3) = ESL/QA(NQN)
                SGLB(NQN+LQN,2*LQN  ,N     ,4) = ESS/QA(NQN)
                SGLB(NQN+LQN,2*LQN  ,N     ,5) = ETOT/QA(NQN)
                SGLB(NQN+LQN,2*LQN  ,N+NBAS,1) = WLL/QA(NQN)
                SGLB(NQN+LQN,2*LQN  ,N+NBAS,2) = WLS/QA(NQN)
                SGLB(NQN+LQN,2*LQN  ,N+NBAS,3) = WSL/QA(NQN)
                SGLB(NQN+LQN,2*LQN  ,N+NBAS,4) = WSS/QA(NQN)
                SGLB(NQN+LQN,2*LQN  ,N+NBAS,5) = WTOT/QA(NQN)
C
              ENDDO
C
            ENDIF
C
C           UPDATE NUMBER OF OCCUPIED PAIRS
            IOCCML = IOCCML+2
C
C           SKIP STEP FOR UNOCCUPIED ORBITALS
101         CONTINUE
C
C         END LOOP OVER LISTED NQNS
          ENDDO
C
        ENDDO
C
110     CONTINUE
C
C       NEGATIVE KQN CHOICE
        KQN =-LQN-1
C
C       DEGENERACY
        IDGN = 2*LQN+2
C
C       SPECIAL CASE: ONE-ELECTRON SYSTEM
        IF(NOCC.EQ.1) THEN
          DO NQN=1,NUMOCC(LQN)
            IDGN    = 1
            QA(NQN) = 1.0D0
          ENDDO
        ENDIF
C
C       CALCULATE BATCH OF MATRIX ELEMENT INTEGRALS
        IF(H2INT.EQ.'SLFHI0') THEN
          CALL SLFHI0(V1,EXL,IZ,KQN,NBAS)
        ELSEIF(H2INT.EQ.'UEENUC') THEN
          CALL UEHNC0(V1,EXL,IZ,KQN,NBAS,'ELEC')
        ELSEIF(H2INT.EQ.'UEUNUC') THEN
          CALL UEHNC0(V1,EXL,IZ,KQN,NBAS,'MUON')
        ELSEIF(H2INT.EQ.'UEPNUC') THEN
          CALL UEHNC0(V1,EXL,IZ,KQN,NBAS,'PION')
        ELSEIF(H2INT.EQ.'WKRNUC') THEN
          CALL WKRNC0(V1,EXL,IZ,KQN,NBAS)
        ELSEIF(H2INT.EQ.'KSBNUC') THEN
          CALL KSBNC0(V1,EXL,IZ,KQN,NBAS)
        ELSEIF(H2INT.EQ.'UEHSCF') THEN
          CALL UEHEL0(V1,EXL,IZ,KQN,NBAS)
        ENDIF
C
C       LOOP OVER |MQN| CHOICES
        DO MQN=1,LQN+1
C
C         COEFFICIENT MATRIX OFFSET
          NOCC = IOCTP(IZ,2*LQN+1,MQN,1)
C
C         LOOP OVER LISTED NQNS
          DO NQN=1,NUMOCC(LQN)
C
C           SKIP THIS STEP IF ORBITAL IS UNOCCUPIED
            IF(NORB(LQN,NQN).EQ.0) GOTO 102
C
C           SPHERICAL SYMMETRY: FOCUS ON FIRST MQN ONLY
            IF(MQN.EQ.1) THEN
C
C             ADDRESS AND AVERAGE OCCUPATION
              IGLB(NQN+LQN,2*LQN+1) = IOCPN(IOCCML+1)
              PLTN(NQN+LQN,2*LQN+1) = QA(NQN)
C
C             FOCK MATRIX ADDRESSES
              IL = LRGE(IZ,2*LQN+1,MQN)
              IS = LRGE(IZ,2*LQN+1,MQN)+NSKP
C
C             COEFFICIENT MATRIX OFFSET
              MOCC = IOCTP(IZ,2*LQN+1,MQN,NQN)+NSKP
C
C             LOOP OVER ALL STATES OF THIS TYPE
              DO N=1,NBAS
C
C               RADIAL MOMENTS
                ELL = 0.0D0
                ELS = 0.0D0
                ESL = 0.0D0
                ESS = 0.0D0
                WLL = 0.0D0
                WLS = 0.0D0
                WSL = 0.0D0
                WSS = 0.0D0
                DO IBAS=1,NBAS
                  DO JBAS=1,NBAS
C
C                   LARGE-COMPONENT OFFSETS
                    IA = IL+IBAS
                    JA = IL+JBAS
C
C                   SMALL-COMPONENT OFFSETS
                    IB = IS+IBAS
                    JB = IS+JBAS
C
C                   COEFFICIENT MATRIX OFFSET
                    NOCC = IOCTP(IZ,2*LQN+1,MQN,1)-1+N
C
                    DAA = DREAL(DCONJG(COEF(IA,MOCC))*COEF(JA,NOCC))
                    DAB = DREAL(DCONJG(COEF(IA,MOCC))*COEF(JB,NOCC))
                    DBA = DREAL(DCONJG(COEF(IB,MOCC))*COEF(JA,NOCC))
                    DBB = DREAL(DCONJG(COEF(IB,MOCC))*COEF(JB,NOCC))
C
                    KBAS = IBAS+NBAS
                    LBAS = JBAS+NBAS
C
                    ELL = ELL + DAA*V1(IBAS,JBAS)
                    ELS = ELS + DAB*V1(IBAS,LBAS)
                    ESL = ESL + DBA*V1(KBAS,JBAS)
                    ESS = ESS + DBB*V1(KBAS,LBAS)
C
C                   COEFFICIENT MATRIX OFFSET
                    NOCC = IOCTP(IZ,2*LQN+1,MQN,1)-1+NSKP+N
C
                    DAA = DREAL(DCONJG(COEF(IA,MOCC))*COEF(JA,NOCC))
                    DAB = DREAL(DCONJG(COEF(IA,MOCC))*COEF(JB,NOCC))
                    DBA = DREAL(DCONJG(COEF(IB,MOCC))*COEF(JA,NOCC))
                    DBB = DREAL(DCONJG(COEF(IB,MOCC))*COEF(JB,NOCC))
C
                    KBAS = IBAS+NBAS
                    LBAS = JBAS+NBAS
C
                    WLL = WLL + DAA*V1(IBAS,JBAS)
                    WLS = WLS + DAB*V1(IBAS,LBAS)
                    WSL = WSL + DBA*V1(KBAS,JBAS)
                    WSS = WSS + DBB*V1(KBAS,LBAS)
C
                  ENDDO
                ENDDO

                ETOT = ELL+ELS+ESL+ESS
                WTOT = WLL+WLS+WSL+WSS
                RGLB(NQN+LQN,2*LQN+1,1) = ELL/QA(NQN)
                RGLB(NQN+LQN,2*LQN+1,2) = ELS/QA(NQN)
                RGLB(NQN+LQN,2*LQN+1,3) = ESL/QA(NQN)
                RGLB(NQN+LQN,2*LQN+1,4) = ESS/QA(NQN)
                RGLB(NQN+LQN,2*LQN+1,5) = ETOT/QA(NQN)
                SGLB(NQN+LQN,2*LQN+1,N     ,1) = ELL/QA(NQN)
                SGLB(NQN+LQN,2*LQN+1,N     ,2) = ELS/QA(NQN)
                SGLB(NQN+LQN,2*LQN+1,N     ,3) = ESL/QA(NQN)
                SGLB(NQN+LQN,2*LQN+1,N     ,4) = ESS/QA(NQN)
                SGLB(NQN+LQN,2*LQN+1,N     ,5) = ETOT/QA(NQN)
                SGLB(NQN+LQN,2*LQN+1,N+NBAS,1) = WLL/QA(NQN)
                SGLB(NQN+LQN,2*LQN+1,N+NBAS,2) = WLS/QA(NQN)
                SGLB(NQN+LQN,2*LQN+1,N+NBAS,3) = WSL/QA(NQN)
                SGLB(NQN+LQN,2*LQN+1,N+NBAS,4) = WSS/QA(NQN)
                SGLB(NQN+LQN,2*LQN+1,N+NBAS,5) = WTOT/QA(NQN)
C
              ENDDO
C
            ENDIF
C
C           UPDATE NUMBER OF OCCUPIED PAIRS
            IOCCML = IOCCML+2
C
C           SKIP STEP FOR UNOCCUPIED ORBITALS
102         CONTINUE
C
C         END LOOP OVER LISTED NQNS
          ENDDO
C
        ENDDO
C
C     END LOOP OVER LQNS
      ENDDO
C
C**********************************************************************C
C     NOW ALL REQUIRED V(1) HAVE BEEN CALCULATED.                      C
C**********************************************************************C
C
C     CALCULATE W+(I) = SUM_J |<I|V|J>|^2/(E(I)-E(J))
C      DO N=1,NBAS
C        VIJ = SGLB(NQN+LQN,2*LQN+1,N     ,1)
C        IS  = LRGE(IZ,KA,1)+NSKP
C        EIJ = EIGN(IS+IOCC)-EIGN(IS+N)
C        RGLB(NQN+LQN,2*LQN+1,1) = VIJ*VIJ/EIJ
C      ENDDO
       
C
C**********************************************************************C
C     PRINT SPECTRUM (POSITIVE ENERGY)                                 C
C**********************************************************************C
C
C     PRINT TITLE TO TERMINAL/FILE
20    FORMAT(1X,' nl  |',11X,'E(LL)',11X,
     &                 'E(LS)',11X,'E(SL)',11X,'E(SS) |',9X,'E(TOT)')
21    FORMAT(1X,I2,A,A,' |',4(1X,F15.10),' |',F15.10)
22    FORMAT(25X,'Interaction ',A,' by atomic orbitals')
23    FORMAT(30X,'Second order in RSPT: ',A)
24    FORMAT(1X,A,' |',4(1X,F15.10),' |',F15.10)
C
      WRITE(6,22) H2INT
      WRITE(7,22) H2INT
      WRITE(6,23) 'E(2+)'
      WRITE(7,23) 'E(2+)'
      WRITE(6, *) REPEAT('=',87)
      WRITE(7, *) REPEAT('=',87)
      WRITE(6,20)
      WRITE(7,20)
      WRITE(6, *) REPEAT('-',87)
      WRITE(7, *) REPEAT('-',87)
C
C     IMPORT ATOMIC CHARGE DETAILS
      IZNV = INT(ZNUC(IZ))
      ICRG = IQNC(IZ)
C
C     INITIALISE THE ORBITAL FILLING ARRAY
      DO LQN=0,MEL
        DO NQN=1,MKP+1
          NORB(LQN,NQN) = 0
        ENDDO
      ENDDO
C
C     READ OR GENERATE AUFBAU OCCUPANCY FOR THIS ATOM
      IF(CNFG(IZ).EQ.'MANUAL') THEN
        LMXCONF = (NKAP(IZ)-1)/2
        DO L=0,LMXCONF
          NUMOCC(L) = NLVL(IZ,L)
          DO N=1,NLVL(IZ,L)
            NORB(L,N) = NCNF(IZ,L,N)
          ENDDO
        ENDDO
      ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
        CALL AUFBAU(IZNV,ICRG,NORB,NUMOCC,LMXCONF)
      ELSE
        WRITE(6, *) 'In H2INT0: invalid configuration choice.'
        WRITE(7, *) 'In H2INT0: invalid configuration choice.'
      ENDIF
C
C     LOOP OVER ORBITALS IN THE USUAL AUFBAU ORDER
      DO NQN=1,MKP
        DO LQN=0,NQN-1
C
          IF(LQN.GT.MEL) GOTO 150
          IF(NORB(LQN,NQN-LQN).EQ.0) GOTO 150
C
C         POSITIVE KQN CHOICE (APPLIES ONLY FOR LQN>0)
          IF(LQN.EQ.0) GOTO 120
C
C         SYMMETRY LABEL
          CHS = '-'
C
C         DEGENERACY
          IDGN = 2*LQN
C
C         IDENTIFY THE RIGHT INDEX AND PROPERTIES
          MOCC = IGLB(NQN,2*LQN  )+NSKP
          FRAC = PLTN(NQN,2*LQN  )
          RAVA = RGLB(NQN,2*LQN  ,1)
          RAVB = RGLB(NQN,2*LQN  ,2)
          RAVC = RGLB(NQN,2*LQN  ,3)
          RAVD = RGLB(NQN,2*LQN  ,4)
          RAVE = RGLB(NQN,2*LQN  ,5)
C
C         OUTPUT TO TERMINAL
          WRITE(6,21) NQN,LLAB(LQN),CHS,RAVA,RAVB,RAVC,RAVD,RAVE
          WRITE(7,21) NQN,LLAB(LQN),CHS,RAVA,RAVB,RAVC,RAVD,RAVE
C
120       CONTINUE
C
C         NEGATIVE KQN CHOICE (APPLIES ONLY FOR LQN>0)
C
C         SYMMETRY LABEL IF NEEDED
          IF(LQN.EQ.0) THEN
            CHS = ' '
          ELSE
            CHS = '+'
          ENDIF
C
C         DEGENERACY
          IDGN = 2*LQN+2
C
C         IDENTIFY THE RIGHT INDEX AND PROPERTIES
          MOCC = IGLB(NQN,2*LQN+1)+NSKP
          FRAC = PLTN(NQN,2*LQN+1)
          RAVA = RGLB(NQN,2*LQN+1,1)
          RAVB = RGLB(NQN,2*LQN+1,2)
          RAVC = RGLB(NQN,2*LQN+1,3)
          RAVD = RGLB(NQN,2*LQN+1,4)
          RAVE = RGLB(NQN,2*LQN+1,5)
C
C         OUTPUT TO TERMINAL
          WRITE(6,21) NQN,LLAB(LQN),CHS,RAVA,RAVB,RAVC,RAVD,RAVE
          WRITE(7,21) NQN,LLAB(LQN),CHS,RAVA,RAVB,RAVC,RAVD,RAVE
C
150       CONTINUE
C
        ENDDO
C
      ENDDO
C
C     INITIALISE TOTAL ENERGIES
      ELLT = 0.0D0
      ELST = 0.0D0
      ESLT = 0.0D0
      ESST = 0.0D0
C
      ETTT = ELLT+ELST+ESLT+ESST
      WRITE(6, *) REPEAT('-',87)
      WRITE(7, *) REPEAT('-',87)
      WRITE(6,24) 'Atom',ELLT,ELST,ESLT,ESST,ETTT
      WRITE(7,24) 'Atom',ELLT,ELST,ESLT,ESST,ETTT
      WRITE(6, *) REPEAT('=',87)
      WRITE(7, *) REPEAT('=',87)
      WRITE(6, *) ''
      WRITE(7, *) ''
C
C**********************************************************************C
C     PRINT SPECTRUM (NEGATIVE ENERGY)                                 C
C**********************************************************************C
C
      WRITE(6,22) H2INT
      WRITE(7,22) H2INT
      WRITE(6,23) 'E(2-)'
      WRITE(7,23) 'E(2-)'
      WRITE(6, *) REPEAT('=',87)
      WRITE(7, *) REPEAT('=',87)
      WRITE(6,20)
      WRITE(7,20)
      WRITE(6, *) REPEAT('-',87)
      WRITE(7, *) REPEAT('-',87)
C
C     IMPORT ATOMIC CHARGE DETAILS
      IZNV = INT(ZNUC(IZ))
      ICRG = IQNC(IZ)
C
C     INITIALISE THE ORBITAL FILLING ARRAY
      DO LQN=0,MEL
        DO NQN=1,MKP+1
          NORB(LQN,NQN) = 0
        ENDDO
      ENDDO
C
C     READ OR GENERATE AUFBAU OCCUPANCY FOR THIS ATOM
      IF(CNFG(IZ).EQ.'MANUAL') THEN
        LMXCONF = (NKAP(IZ)-1)/2
        DO L=0,LMXCONF
          NUMOCC(L) = NLVL(IZ,L)
          DO N=1,NLVL(IZ,L)
            NORB(L,N) = NCNF(IZ,L,N)
          ENDDO
        ENDDO
      ELSEIF(CNFG(IZ).EQ.'AUFBAU') THEN
        CALL AUFBAU(IZNV,ICRG,NORB,NUMOCC,LMXCONF)
      ELSE
        WRITE(6, *) 'In H2INT0: invalid configuration choice.'
        WRITE(7, *) 'In H2INT0: invalid configuration choice.'
      ENDIF
C
C     LOOP OVER ORBITALS IN THE USUAL AUFBAU ORDER
      DO NQN=1,MKP
        DO LQN=0,NQN-1
C
          IF(LQN.GT.MEL) GOTO 155
          IF(NORB(LQN,NQN-LQN).EQ.0) GOTO 155
C
C         POSITIVE KQN CHOICE (APPLIES ONLY FOR LQN>0)
          IF(LQN.EQ.0) GOTO 125
C
C         SYMMETRY LABEL
          CHS = '-'
C
C         DEGENERACY
          IDGN = 2*LQN
C
C         IDENTIFY THE RIGHT INDEX AND PROPERTIES
          MOCC = IGLB(NQN,2*LQN  )+NSKP
          FRAC = PLTN(NQN,2*LQN  )
          RAVA = RGLB(NQN,2*LQN  ,1)
          RAVB = RGLB(NQN,2*LQN  ,2)
          RAVC = RGLB(NQN,2*LQN  ,3)
          RAVD = RGLB(NQN,2*LQN  ,4)
          RAVE = RGLB(NQN,2*LQN  ,5)
C
C         OUTPUT TO TERMINAL
          WRITE(6,21) NQN,LLAB(LQN),CHS,RAVA,RAVB,RAVC,RAVD,RAVE
          WRITE(7,21) NQN,LLAB(LQN),CHS,RAVA,RAVB,RAVC,RAVD,RAVE
C
125       CONTINUE
C
C         NEGATIVE KQN CHOICE (APPLIES ONLY FOR LQN>0)
C
C         SYMMETRY LABEL IF NEEDED
          IF(LQN.EQ.0) THEN
            CHS = ' '
          ELSE
            CHS = '+'
          ENDIF
C
C         DEGENERACY
          IDGN = 2*LQN+2
C
C         IDENTIFY THE RIGHT INDEX AND PROPERTIES
          MOCC = IGLB(NQN,2*LQN+1)+NSKP
          FRAC = PLTN(NQN,2*LQN+1)
          RAVA = RGLB(NQN,2*LQN+1,1)
          RAVB = RGLB(NQN,2*LQN+1,2)
          RAVC = RGLB(NQN,2*LQN+1,3)
          RAVD = RGLB(NQN,2*LQN+1,4)
          RAVE = RGLB(NQN,2*LQN+1,5)
C
C         OUTPUT TO TERMINAL
          WRITE(6,21) NQN,LLAB(LQN),CHS,RAVA,RAVB,RAVC,RAVD,RAVE
          WRITE(7,21) NQN,LLAB(LQN),CHS,RAVA,RAVB,RAVC,RAVD,RAVE
C
155       CONTINUE
C
        ENDDO
C
      ENDDO
C
C     INITIALISE TOTAL ENERGIES
      ELLT = 0.0D0
      ELST = 0.0D0
      ESLT = 0.0D0
      ESST = 0.0D0
C
      ETTT = ELLT+ELST+ESLT+ESST
      WRITE(6, *) REPEAT('-',87)
      WRITE(7, *) REPEAT('-',87)
      WRITE(6,24) 'Atom',ELLT,ELST,ESLT,ESST,ETTT
      WRITE(7,24) 'Atom',ELLT,ELST,ESLT,ESST,ETTT
      WRITE(6, *) REPEAT('=',87)
      WRITE(7, *) REPEAT('=',87)
      WRITE(6, *) ''
      WRITE(7, *) ''
C
      RETURN
      END
C
C
      SUBROUTINE G1PAIR0(IZ,IA,G2INT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        GGGGGG   11  PPPPPPPP     AA    IIII RRRRRRR   000000         C
C       GG    GG 111  PP     PP   AAAA    II  RR    RR 00   000        C
C       GG        11  PP     PP  AA  AA   II  RR    RR 00  0000        C
C       GG        11  PP     PP AA    AA  II  RR    RR 00 00 00        C
C       GG   GGG  11  PPPPPPPP  AAAAAAAA  II  RRRRRRR  0000  00        C
C       GG    GG  11  PP        AA    AA  II  RR    RR 000   00        C
C        GGGGGG  1111 PP        AA    AA IIII RR    RR  000000         C
C                                                                      C
C -------------------------------------------------------------------- C
C  G1PAIR0 PERFORMS A FIRST-ORDER ANALYSIS OF A 2-BODY INTERACTION     C
C  APPLIED TO AN AVERAGE-OVER-CONFIGURATION ATOMIC SOLUTION.           C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ IZ    - ATOMIC CENTRE INDEX.                                      C
C  ▶ IA    - ANALYSIS MODE (1 - SYMMETRY TYPES, 2 - ATOMIC ORBITALS).  C
C  ▶ G2INT - NAME OF TWO-BODY INTERACTION (COULM OR BREIT).            C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*1 LLAB,QSGN
      CHARACTER*2 ELMT(120),ELNM,KLAB
      CHARACTER*5 NMDL,G2INT
      CHARACTER*6 CNFG
C
      DIMENSION QA(MKP),QE(MKP)
      DIMENSION CL(MBS,MKP,MKP+1),CS(MBS,MKP,MKP+1)
      DIMENSION DALL(MB2,MKP,MKP+1),DALS(MB2,MKP,MKP+1),
     &          DASL(MB2,MKP,MKP+1),DASS(MB2,MKP,MKP+1)
      DIMENSION DELL(MB2,MKP,MKP+1),DELS(MB2,MKP,MKP+1),
     &          DESL(MB2,MKP,MKP+1),DESS(MB2,MKP,MKP+1)
      DIMENSION DNALL(MB2,MKP),DNALS(MB2,MKP),
     &          DNASL(MB2,MKP),DNASS(MB2,MKP)
      DIMENSION DNELL(MB2,MKP),DNELS(MB2,MKP),
     &          DNESL(MB2,MKP),DNESS(MB2,MKP)
      DIMENSION GSET(MBD,MBD,0:1,0:1)
      DIMENSION EGSTT(MKP,MKP,5),EGS(5)
      DIMENSION EGOTT(MKP,MKP,MKP+1,MKP+1,5),EGO(5)
C
      COMPLEX*16 COEF(MDM,MDM)
C
      COMMON/ATMB/B11(MBD,MBD),B21(MBD,MBD),B12(MBD,MBD),B22(MBD,MBD)
      COMMON/ATMC/G11(MBD,MBD),G21(MBD,MBD),G12(MBD,MBD),G22(MBD,MBD)
      COMMON/ATMD/DLL1(MB2),DSL1(MB2),DSS1(MB2),DLS1(MB2),
     &            DLL2(MB2),DSL2(MB2),DSS2(MB2),DLS2(MB2)
      COMMON/AUFB/NORB(0:MEL,MKP+1),NUMOCC(0:MEL),LMXCONF
      COMMON/B0QN/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EWDR,EWXC,EANM,ESLF,EUEE,
     &            EUEU,EUEP,EWKR,EKSB,ESCF
      COMMON/FILL/NCNF(MCT,0:MEL,MKP+1),NLVL(MCT,0:MEL),CNFG(MCT)
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
C
C     CONDITIONS UNDER WHICH THIS MAKES NO SENSE
      IF(HMLT.EQ.'NORL'.AND.G2INT.EQ.'BREIT') THEN
        WRITE(6, *) 'In G1PAIR0: cannot apply Breit to NORL.'
        WRITE(7, *) 'In G1PAIR0: cannot apply Breit to NORL.'
      ELSEIF(G2INT.NE.'COULM'.AND.G2INT.NE.'BREIT') THEN
        WRITE(6, *) 'In G1PAIR0: select a valid operator.'
        WRITE(7, *) 'In G1PAIR0: select a valid operator.'
      ENDIF
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
600   FORMAT(1X,88('='),/,27X,A,1X,A,1X,A,/,1X,88('='))
601   FORMAT(1X,'(κA|κB)',9X,'E(LL|LL)',8X,'E(LL|SS)',8X,
     &                       'E(SS|LL)',8X,'E(SS|SS)',10X,'E(TOT)')
602   FORMAT(1X,'(κA|κB)',9X,'E(LL|SS)',8X,'E(LS|LS)',8X,
     &                       'E(SL|LS)',8X,'E(SS|LL)',10X,'E(TOT)')
603   FORMAT(1X,'(',A,'|',A,')',3X,5(F14.10,2X))
604   FORMAT(1X,'Atomic',2X,5(2X,F14.10))
C
C     PRINT HEADER TO TERMINAL
      WRITE(6,600) 'Interaction',G2INT,'by symmetry pairs'
      WRITE(7,600) 'Interaction',G2INT,'by symmetry pairs'
      IF(G2INT.EQ.'COULM') THEN
        WRITE(6,601)
        WRITE(7,601)
      ELSEIF(G2INT.EQ.'BREIT') THEN
        WRITE(6,602)
        WRITE(7,602)
      ENDIF
      WRITE(6,*) REPEAT('-',88)
      WRITE(7,*) REPEAT('-',88)
C
C     INITIALISE BREIT COUNTERS
      DO ITT=1,5
        EGS(ITT) = 0.0D0
        DO KPA=1,MKP
          DO KPB=1,MKP
            EGSTT(KPA,KPB,ITT) = 0.0D0
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
C       FULL MATRIX DIMENSION DEPENDS ON HMLT
        IF(HMLT.EQ.'NORL') THEN
          NMATA =   NBASA
        ELSE
          NMATA = 2*NBASA
        ENDIF
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
C         BATCH OF TWO-BODY INTEGRALS
          IF(G2INT.EQ.'COULM') THEN
C
C           GENERATE THE MEAN-FIELD ATOMIC COULOMB MATRIX OVER DENSITIES
            CALL COULM0
C
C           TRANSFER RESULTS INTO LOCAL ARRAY
            DO IBAS=1,NMATA
              DO JBAS=1,NMATA
                GSET(IBAS,JBAS,1,1) = G22(IBAS,JBAS)
                GSET(IBAS,JBAS,1,0) = G21(IBAS,JBAS)
                GSET(IBAS,JBAS,0,1) = G12(IBAS,JBAS)
                GSET(IBAS,JBAS,0,0) = G11(IBAS,JBAS)
              ENDDO
            ENDDO
C
          ELSEIF(G2INT.EQ.'BREIT') THEN
C
C           GENERATE THE MEAN-FIELD ATOMIC BREIT MATRIX OVER DENSITIES
            CALL BREIT0
C
C           TRANSFER RESULTS INTO LOCAL ARRAY
            DO IBAS=1,NMATA
              DO JBAS=1,NMATA
                GSET(IBAS,JBAS,1,1) = B22(IBAS,JBAS)
                GSET(IBAS,JBAS,1,0) = B21(IBAS,JBAS)
                GSET(IBAS,JBAS,0,1) = B12(IBAS,JBAS)
                GSET(IBAS,JBAS,0,0) = B11(IBAS,JBAS)
              ENDDO
            ENDDO
C          
          ENDIF
C
C         LOOP OVER SYMMETRY TYPES FOR THIS LQNA
          DO 650 ISA=0,1
C
          KQNA = ((-1)**ISA)*LQNA - ISA
          RK2A = DFLOAT(2*IABS(KQNA))
          KPA  = 2*LQNA+ISA
C
          IF(KQNA.EQ.0.OR.HMLT.EQ.'NORL') GOTO 651
C
C         LOOP OVER SYMMETRY TYPES FOR THIS LQNB
          DO 660 ISB=0,1
C
          KQNB = ((-1)**ISB)*LQNB - ISB
          RK2B = DFLOAT(2*IABS(KQNB))
          KPB  = 2*LQNB+ISB
C
          IF(KQNB.EQ.0.OR.HMLT.EQ.'NORL') GOTO 661
C
C         INITIALISE SOME COUNTERS
          GLL = 0.0D0
          GLS = 0.0D0
          GSL = 0.0D0
          GSS = 0.0D0
C
C         SUM OVER THIS BREIT MATRIX AND MULTIPLY BY DENSITY
C         TODO: I SWAPPED THESE OUT TO MATCH UP WITH BREIT1...
          M = 0
          DO IBAS=1,NBASA
            DO JBAS=1,NBASA
              M = M+1
C
C             SMALL COMPONENT BLOCK ADDRESSES
              KBAS = IBAS+NBASA
              LBAS = JBAS+NBASA
C
              GSS = GSS + GSET(IBAS,JBAS,ISA,ISB)*DNALL(M,KPA)
              GSL = GSL + GSET(IBAS,LBAS,ISA,ISB)*DNALS(M,KPA)
              GLS = GLS + GSET(KBAS,JBAS,ISA,ISB)*DNASL(M,KPA)
              GLL = GLL + GSET(KBAS,LBAS,ISA,ISB)*DNASS(M,KPA)
C
            ENDDO
          ENDDO
C
C         MULTIPLY BY SUBSHELL DEGENERACY FACTOR
          EGSTT(KPA,KPB,1) = 0.5D0*RK2A*RK2B*GLL
          EGSTT(KPA,KPB,2) = 0.5D0*RK2A*RK2B*GLS
          EGSTT(KPA,KPB,3) = 0.5D0*RK2A*RK2B*GSL
          EGSTT(KPA,KPB,4) = 0.5D0*RK2A*RK2B*GSS
C
C         TOTAL VALUE              
          EGSTT(KPA,KPB,5) = EGSTT(KPA,KPB,1) + EGSTT(KPA,KPB,2)
     &                     + EGSTT(KPA,KPB,3) + EGSTT(KPA,KPB,4)
C
C         ADD TO NET COMPONENT-TYPE COUNTER
          DO ITT=1,5
            EGS(ITT) = EGS(ITT) + EGSTT(KPA,KPB,ITT)
          ENDDO
C
661       CONTINUE
660       CONTINUE
C
651       CONTINUE
650       CONTINUE
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
          WRITE(6,603) KLAB(KAPA(IKP,IZ)),KLAB(KAPA(JKP,IZ)),
     &                                    (EGSTT(IKP,JKP,ITT),ITT=1,5)
          WRITE(7,603) KLAB(KAPA(IKP,IZ)),KLAB(KAPA(JKP,IZ)),
     &                                    (EGSTT(IKP,JKP,ITT),ITT=1,5)
        ENDDO
      ENDDO
      WRITE(6,*) REPEAT('-',88)
      WRITE(7,*) REPEAT('-',88)
      WRITE(6,604) (EGS(ITT),ITT=1,5)
      WRITE(7,604) (EGS(ITT),ITT=1,5)
      WRITE(6,*) REPEAT('=',88)
      WRITE(7,*) REPEAT('=',88)
      WRITE(6,*) ' '
      WRITE(7,*) ' '
C
C     ATOMIC TWO-BODY INTERACTION ENERGY
      E2 = EGS(5)
C
699   CONTINUE
C
C**********************************************************************C
C     BREIT ENERGY BY ATOMIC ORBITAL PAIRS                             C
C**********************************************************************C
C
      IF(IA.EQ.1) GOTO 799
C
700   FORMAT(1X,92('='),/,26X,A,1X,A,1X,A,/,1X,92('='))
701   FORMAT(1X,'(n κA|n κB)',9X,'E(LL|LL)',8X,'E(LL|SS)',8X,
     &                       'E(SS|LL)',8X,'E(SS|SS)',10X,'E(TOT)')
702   FORMAT(1X,'(n κA|n κB)',9X,'E(LL|SS)',8X,'E(LS|LS)',8X,
     &                       'E(SL|LS)',8X,'E(SS|LL)',10X,'E(TOT)')
703   FORMAT(1X,'(',I2,A,'|',I2,A,')',3X,5(F14.10,2X))
704   FORMAT(1X,'Atomic',6X,5(2X,F14.10))
C
C     PRINT HEADER TO TERMINAL
      WRITE(6,700) 'Interaction',G2INT,'by atomic orbital pairs'
      WRITE(7,700) 'Interaction',G2INT,'by atomic orbital pairs'
      IF(G2INT.EQ.'COULM') THEN
        WRITE(6,701)
        WRITE(7,701)
      ELSEIF(G2INT.EQ.'BREIT') THEN
        WRITE(6,702)
        WRITE(7,702)
      ENDIF
      WRITE(6,*) REPEAT('-',92)
      WRITE(7,*) REPEAT('-',92)
C
C     INITIALISE BREIT COUNTERS
      DO ITT=1,5
        EGO(ITT) = 0.0D0
        DO KPA=1,MKP
          DO KPB=1,MKP
            DO IOCC=1,MKP+1
              DO JOCC=1,MKP+1
                EGOTT(KPA,KPB,IOCC,JOCC,ITT) = 0.0D0
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     LOOP OVER ALL OCCUPIED LQNA VALUES
      DO LQNA=0,LMXCONF
C
C       RECORD LQNA VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
        NBASA = NFNC(LQNA,IZ)
        DO IBAS=1,NBASA
          EXLA(IBAS) = BEXL(IBAS,LQNA,IZ)
        ENDDO
C
C       FULL MATRIX DIMENSION DEPENDS ON HMLT
        IF(HMLT.EQ.'NORL') THEN
          NMATA =   NBASA
        ELSE
          NMATA = 2*NBASA
        ENDIF
C
C       LOOP OVER ALL OCCUPIED LQNB VALUES
        DO LQNB=0,LMXCONF
C
C         RECORD LQNB VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
          NBASB = NFNC(LQNB,IZ)
          DO JBAS=1,NBASB
            EXLB(JBAS) = BEXL(JBAS,LQNB,IZ)
          ENDDO
C
C         LOOP OVER ALL OCCUPIED NQNS FOR LQNB TYPE
          DO JOCC=1,NUMOCC(LQNB)
C
C           NUMBER OF BASIS FUNCTION OVERLAPS IN THIS BLOCK
            MAXM = NBASB*NBASB
C
C           RELEVANT DENSITY MATRIX DEPENDS ON LQNA AND LQNB
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
C           BATCH OF TWO-BODY INTEGRALS
            IF(G2INT.EQ.'COULM') THEN
C
C             GENERATE THE MEAN-FIELD ATOMIC COULOMB MATRIX
              CALL COULM0
C
C             TRANSFER RESULTS INTO LOCAL ARRAY
              DO IBAS=1,NMATA
                DO JBAS=1,NMATA
                  GSET(IBAS,JBAS,1,1) = G22(IBAS,JBAS)
                  GSET(IBAS,JBAS,1,0) = G21(IBAS,JBAS)
                  GSET(IBAS,JBAS,0,1) = G12(IBAS,JBAS)
                  GSET(IBAS,JBAS,0,0) = G11(IBAS,JBAS)
                ENDDO
              ENDDO
C
            ELSEIF(G2INT.EQ.'BREIT') THEN
C
C             GENERATE THE MEAN-FIELD ATOMIC BREIT MATRIX
              CALL BREIT0
C
C             TRANSFER RESULTS INTO LOCAL ARRAY
              DO IBAS=1,NMATA
                DO JBAS=1,NMATA
                  GSET(IBAS,JBAS,1,1) = B22(IBAS,JBAS)
                  GSET(IBAS,JBAS,1,0) = B21(IBAS,JBAS)
                  GSET(IBAS,JBAS,0,1) = B12(IBAS,JBAS)
                  GSET(IBAS,JBAS,0,0) = B11(IBAS,JBAS)
                ENDDO
              ENDDO
C          
            ENDIF
C
C           LOOP OVER ALL OCCUPIED NQNS FOR LQNA TYPE
            DO IOCC=1,NUMOCC(LQNA)
C
C             LOOP OVER SYMMETRY TYPES FOR THIS LQNA
              DO 750 ISA=0,1
C
              KQNA = ((-1)**ISA)*LQNA - ISA
              RK2A = DFLOAT(2*IABS(KQNA))
              KPA  = 2*LQNA+ISA
C
              IF(KQNA.EQ.0.OR.HMLT.EQ.'NORL') GOTO 751
C
C             LOOP OVER SYMMETRY TYPES FOR THIS LQNB
              DO 760 ISB=0,1
C
              KQNB = ((-1)**ISB)*LQNB - ISB
              RK2B = DFLOAT(2*IABS(KQNB))
              KPB  = 2*LQNB+ISB
C
              IF(KQNB.EQ.0.OR.HMLT.EQ.'NORL') GOTO 761
C
C             INITIALISE SOME COUNTERS
              GLL = 0.0D0
              GLS = 0.0D0
              GSL = 0.0D0
              GSS = 0.0D0
C
C             SUM OVER THIS BREIT MATRIX AND MULTIPLY BY DENSITY
              M = 0
              DO IBAS=1,NBASA
                DO JBAS=1,NBASA
                  M = M+1
C
C                 SMALL COMPONENT BLOCK ADDRESSES
                  KBAS = IBAS+NBASA
                  LBAS = JBAS+NBASA
C
                  GLL = GLL + GSET(IBAS,JBAS,ISA,ISB)*DALL(M,KPA,IOCC)
                  GLS = GLS + GSET(IBAS,LBAS,ISA,ISB)*DALS(M,KPA,IOCC)
                  GSL = GSL + GSET(KBAS,JBAS,ISA,ISB)*DASL(M,KPA,IOCC)
                  GSS = GSS + GSET(KBAS,LBAS,ISA,ISB)*DASS(M,KPA,IOCC)
C
                ENDDO
              ENDDO
C
C             MULTIPLY BY SUBSHELL DEGENERACY FACTOR
              EGOTT(KPA,KPB,IOCC,JOCC,1) = 0.5D0*RK2A*RK2B*GLL
              EGOTT(KPA,KPB,IOCC,JOCC,2) = 0.5D0*RK2A*RK2B*GLS
              EGOTT(KPA,KPB,IOCC,JOCC,3) = 0.5D0*RK2A*RK2B*GSL
              EGOTT(KPA,KPB,IOCC,JOCC,4) = 0.5D0*RK2A*RK2B*GSS
C
C             TOTAL VALUE              
              EGOTT(KPA,KPB,IOCC,JOCC,5)
     &        = EGOTT(KPA,KPB,IOCC,JOCC,1) + EGOTT(KPA,KPB,IOCC,JOCC,2)
     &        + EGOTT(KPA,KPB,IOCC,JOCC,3) + EGOTT(KPA,KPB,IOCC,JOCC,4)
C
C             ADD TO NET COMPONENT-TYPE COUNTER
              DO ITT=1,5
                EGO(ITT) = EGO(ITT) + EGOTT(KPA,KPB,IOCC,JOCC,ITT)
              ENDDO
C
761           CONTINUE
760           CONTINUE
C
751           CONTINUE
750           CONTINUE
C
            ENDDO
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
              WRITE(6,703) IOCC+LQNA,KLAB(KAPA(IKP,IZ)),
     &                     JOCC+LQNB,KLAB(KAPA(JKP,IZ)),
     &                     (EGOTT(IKP,JKP,IOCC,JOCC,ITT),ITT=1,5)
              WRITE(7,703) IOCC+LQNA,KLAB(KAPA(IKP,IZ)),
     &                     JOCC+LQNB,KLAB(KAPA(JKP,IZ)),
     &                     (EGOTT(IKP,JKP,IOCC,JOCC,ITT),ITT=1,5)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      WRITE(6,*) REPEAT('-',92)
      WRITE(7,*) REPEAT('-',92)
      WRITE(6,704) (EGO(ITT),ITT=1,5)
      WRITE(7,704) (EGO(ITT),ITT=1,5)
      WRITE(6,*) REPEAT('=',92)
      WRITE(7,*) REPEAT('=',92)
      WRITE(6,*) ' '
      WRITE(7,*) ' '
C
C     ATOMIC TWO-BODY INTERACTION ENERGY
      E2 = EGO(5)
C
799   CONTINUE
C
C**********************************************************************C
C     UPDATE CORRECT ENERGY COUNTER AND EXIT                           C
C**********************************************************************C
C
C     RECORD TOTAL BREIT ENERGY
      IF(HMLT.EQ.'BARE') THEN
        IF(G2INT.EQ.'COULM') THEN
          ECLG = ECLG + E2
        ENDIF
      ENDIF
      
      IF(HMLT.NE.'DHFB'.AND.HMLT.NE.'DHFQ') THEN
        IF(G2INT.EQ.'BREIT') THEN
          EBRG = EBRG + E2
        ENDIF
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE OVRLP0(OMAT,EXL,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          OOOOOO  VV    VV RRRRRRR  LL      PPPPPPP   000000          C
C         OO    OO VV    VV RR    RR LL      PP    PP 00   000         C
C         OO    OO VV    VV RR    RR LL      PP    PP 00  0000         C
C         OO    OO VV    VV RR    RR LL      PP    PP 00 00 00         C
C         OO    OO  VV  VV  RRRRRRR  LL      PPPPPPP  0000  00         C
C         OO    OO   VVVV   RR    RR LL      PP       000   00         C
C          OOOOOO     VV    RR    RR LLLLLLL PP        000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  OVRLP0 CALCULATES THE ATOMIC OVERLAP MATRIX FOR SYMMETRY TYPE KQN.  C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      DIMENSION RN(MBS*MBS,4),OMAT(MBD,MBD),EXL(MBS)
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
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
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
          T32 = RL+1.5D0
          T52 = RL+2.5D0
          E32 = EIJ**T32
          E52 = EIJ**T52
          SLL = 0.5D0*RN(M,1)*GAMHLF(2*LQN+3)/E32
          SSS = 2.0D0*RN(M,4)*GAMHLF(2*LQN+5)*EPR/E52
C
C         OVERLAP MATRIX ELEMENTS
          OMAT(IBAS     ,JBAS     ) = SLL
          IF(HMLT.EQ.'NORL') GOTO 50
          OMAT(IBAS+NBAS,JBAS     ) = 0.0D0
          OMAT(JBAS     ,IBAS+NBAS) = 0.0D0
          OMAT(IBAS+NBAS,JBAS+NBAS) = SSS
50        CONTINUE
C
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ONEEL0(TMAT,VMAT,EXL,IZ,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          OOOOOO  NN    NN EEEEEEEE EEEEEEEE LL      000000           C
C         OO    OO NNN   NN EE       EE       LL     00   000          C
C         OO    OO NNNN  NN EE       EE       LL     00  0000          C
C         OO    OO NN NN NN EEEEEE   EEEEEE   LL     00 00 00          C
C         OO    OO NN  NNNN EE       EE       LL     0000  00          C
C         OO    OO NN   NNN EE       EE       LL     000   00          C
C          OOOOOO  NN    NN EEEEEEEE EEEEEEEE LLLLLLL 000000           C
C                                                                      C
C -------------------------------------------------------------------- C
C  ONEEL0 CALCULATES THE ATOMIC DIRAC AND OVERLAP MATRICES FOR         C
C  SYMMETRY TYPE KQN, USING EVEN-TEMPERED SGTFS.                       C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RN(MBS*MBS,4),TMAT(MBD,MBD),VMAT(MBD,MBD),EXL(MBS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/FCTS/RFACT(0:80),SFACT(0:80)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      RL  = DFLOAT(LQN)
      TL  = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M    = M+1
          EJ   = EXL(JBAS)
          EIJ  = EI+EJ
          EL1  = EIJ**(LQN+1)
          EPR  = EI*EJ
          T52  = RL+2.5D0
          E52  = EIJ**T52
C
C         LL OVERLAP
          ULL = 0.0D0
          IF(NMDL(IZ).EQ.'POINT') THEN
C           POINT-NUCLEUS APPROXIMATION
            ULL = 0.5D0*GAMHLF(2*LQN+2)/EL1
          ELSE
C           NUCLEAR BEST-FIT EXPANSION
            ULL = FNUC(IZ,0)*ERFINT0(LQN  ,EIJ,XNUC(IZ,0))
            IF(NNUC(IZ).GT.0) THEN
              ULL = ULL + ZCVINT0(IZ,LQN  ,EIJ)
            ENDIF
          ENDIF

          IF(HMLT.EQ.'NORL') THEN
            PLL = GAMHLF(2*LQN+5)*EPR/E52
          ELSE
            PLL = 0.0D0
          ENDIF
C
C         TRANSFER INTO KINETIC AND NUCLEAR ATTRACTION MATRICES
          TMAT(IBAS,JBAS) =           RN(M,1)*PLL
          VMAT(IBAS,JBAS) = -ZNUC(IZ)*RN(M,1)*ULL
C
C         LS,SL AND SS OVERLAPS
          IF(HMLT.EQ.'NORL') GOTO 50
C
C         OVERLAPS, KINETIC ELEMENTS AND SS POTENTIAL INTEGRALS
          SSS = 2.0D0*RN(M,4)*GAMHLF(2*LQN+5)*EPR/E52
          PSL = 2.0D0*RN(M,3)*GAMHLF(2*LQN+5)*EPR/E52
C
          USS = 0.0D0
          IF(NMDL(IZ).EQ.'POINT') THEN
C           POINT-NUCLEUS APPROXIMATION
            VSA = 2.0D0*GAMHLF(2*LQN+4)*EPR/(EL1*EIJ)
            IF(KQN.LT.0) THEN
              VSN = 0.0D0
            ELSE
              VSN = 0.5D0*GAMHLF(2*LQN  )*TL*EIJ/EL1
            ENDIF
            USS = VSA+VSN
          ELSE
C           NUCLEAR BEST-FIT EXPANSION
            VSA = 4.0D0*EPR*ERFINT0(LQN+1,EIJ,XNUC(IZ,0))
            IF(KQN.LT.0) THEN
              VSN = 0.0D0
            ELSE
              VSN =-2.0D0*EIJ*TL*ERFINT0(LQN  ,EIJ,XNUC(IZ,0))
     &                    +TL*TL*ERFINT0(LQN-1,EIJ,XNUC(IZ,0))
            ENDIF
            USS = USS + FNUC(IZ,0)*(VSA+VSN)
            IF(NNUC(IZ).GT.0) THEN
              VSA = 4.0D0*EPR*ZCVINT0(IZ,LQN+1,EIJ)
              IF(KQN.LT.0) THEN
                VSN = 0.0D0
              ELSE
                VSN =-2.0D0*EIJ*TL*ZCVINT0(IZ,LQN  ,EIJ)
     &                      +TL*TL*ZCVINT0(IZ,LQN-1,EIJ)
              ENDIF
              USS = USS + VSA+VSN
            ENDIF
          ENDIF
C
C         SMALL COMPONENT BLOCK ADDRESSES
          KBAS = IBAS+NBAS
          LBAS = JBAS+NBAS
C
C         TRANSFER INTO KINETIC AND NUCLEAR ATTRACTION MATRICES
          TMAT(KBAS,JBAS) = CV*PSL
          TMAT(JBAS,KBAS) = CV*PSL
          VMAT(KBAS,LBAS) =-ZNUC(IZ)*RN(M,4)*USS-2.0D0*EMSS*CV*CV*SSS
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
      SUBROUTINE KINTC0(TMAT,EXL,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           KK    KK IIII NN    NN TTTTTTTT CCCCCC   000000            C
C           KK   KK   II  NNN   NN    TT   CC    CC 00   000           C
C           KK  KK    II  NNNN  NN    TT   CC       00  0000           C
C           KKKKK     II  NN NN NN    TT   CC       00 00 00           C
C           KK  KK    II  NN  NNNN    TT   CC       0000  00           C
C           KK   KK   II  NN   NNN    TT   CC    CC 000   00           C
C           KK    KK IIII NN    NN    TT    CCCCCC   000000            C
C                                                                      C
C -------------------------------------------------------------------- C
C  KINTC0 CALCULATES THE ATOM-CENTRED KINETIC MATRIX ELEMENTS.         C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      DIMENSION RN(MBS*MBS,4),TMAT(MBD,MBD),EXL(MBS)
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      RL  = DFLOAT(LQN)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M    = M+1
          EJ   = EXL(JBAS)
          EIJ  = EI+EJ
          EPR  = EI*EJ
          T52  = RL+2.5D0
          E52  = EIJ**T52
C
C         MATRIX ELEMENT DEPENDS ON HAMILTONIAN
          IF(HMLT.EQ.'NORL') THEN
C
C           NON-RELATIVISTIC CASE
            PLL = RN(M,1)*GAMHLF(2*LQN+5)*EPR/E52
            TMAT(IBAS,JBAS) = PLL
C
          ELSE
C
C           OVERLAPS, KINETIC ELEMENTS AND SS POTENTIAL INTEGRALS
            PSL = 2.0D0*RN(M,3)*GAMHLF(2*LQN+5)*EPR/E52
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NBAS
            LBAS = JBAS+NBAS
C
C           TRANSFER INTO KINETIC MATRIX
            TMAT(KBAS,JBAS) = CV*PSL
            TMAT(JBAS,KBAS) = CV*PSL
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
      SUBROUTINE PNTNC0(VMAT,EXL,IZ,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        PPPPPPP  NN    NN TTTTTTTT NN    NN  CCCCCC   000000          C
C        PP    PP NNN   NN    TT    NNN   NN CC    CC 00   000         C
C        PP    PP NNNN  NN    TT    NNNN  NN CC       00  0000         C
C        PP    PP NN NN NN    TT    NN NN NN CC       00 00 00         C
C        PPPPPPP  NN  NNNN    TT    NN  NNNN CC       0000  00         C
C        PP       NN   NNN    TT    NN   NNN CC    CC 000   00         C
C        PP       NN    NN    TT    NN    NN  CCCCCC   000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  PNTNC0 CALCULATES ATOM-CENTRED NUCLEAR ATTRACTION MATRIX ELEMENTS   C
C  FOR A POINT-NUCLEAR CHARGE MODEL.                                   C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RN(MBS*MBS,4),VMAT(MBD,MBD),EXL(MBS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      RL  = DFLOAT(LQN)
      TL  = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
          EL1 = EIJ**(LQN+1)
          T52 = RL+2.5D0
          E52 = EIJ**T52
C
C         MATRIX ELEMENT DEPENDS ON HAMILTONIAN
          IF(HMLT.EQ.'NORL') THEN
C
C           NON-RELATIVISTIC CASE
            ULL = 0.5D0*GAMHLF(2*LQN+2)/EL1
C
C           TRANSFER INTO NUCLEAR ATTRACTION MATRIX
            VMAT(IBAS,JBAS) = -ZNUC(IZ)*RN(M,1)*ULL
C
          ELSE
C
C           RELATIVISTIC CASE
            SSS = 2.0D0*GAMHLF(2*LQN+5)*EPR/E52
            ULL = 0.5D0*GAMHLF(2*LQN+2)/EL1
C
            VSA = 2.0D0*GAMHLF(2*LQN+4)*EPR/(EL1*EIJ)
            IF(KQN.LT.0) THEN
              VSN = 0.0D0
            ELSE
              VSN = 0.5D0*GAMHLF(2*LQN  )*TL*EIJ/EL1
            ENDIF
            USS = VSA+VSN
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NBAS
            LBAS = JBAS+NBAS
C
C           TRANSFER INTO NUCLEAR ATTRACTION MATRIX
            VMAT(IBAS,JBAS) = -ZNUC(IZ)*RN(M,1)*ULL
            VMAT(KBAS,LBAS) = -ZNUC(IZ)*RN(M,4)*USS
     &                        -2.0D0*EMSS*CV*CV*RN(M,4)*SSS
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
      SUBROUTINE UNINC0(VMAT,EXL,IZ,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          UU    UU NN    NN IIII NN    NN  CCCCCC   000000            C
C          UU    UU NNN   NN  II  NNN   NN CC    CC 00   000           C
C          UU    UU NNNN  NN  II  NNNN  NN CC       00  0000           C
C          UU    UU NN NN NN  II  NN NN NN CC       00 00 00           C
C          UU    UU NN  NNNN  II  NN  NNNN CC       0000  00           C
C          UU    UU NN   NNN  II  NN   NNN CC    CC 000   00           C
C           UUUUUU  NN    NN IIII NN    NN  CCCCCC   000000            C
C                                                                      C
C -------------------------------------------------------------------- C
C  UNINC0 CALCULATES ATOM-CENTRED NUCLEAR ATTRACTION MATRIX ELEMENTS   C
C  FOR A UNIFORMLY-CHARGED SPHERICAL NUCLEAR MODEL.                    C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RN(MBS*MBS,4),VMAT(MBD,MBD),EXL(MBS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      TL  = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
          E02 = EIJ**(LQN)
          E22 = E02*EIJ
          E42 = E22*EIJ
          E52 = E42*DSQRT(EIJ)
C
          X   = 5.0D0*EIJ*RNUC(IZ)*RNUC(IZ)/3.0D0
          FAC = 2.0D0*DSQRT(X)*X
C
C         MATRIX ELEMENT DEPENDS ON HAMILTONIAN
          IF(HMLT.EQ.'NORL') THEN
C
C           NON-RELATIVISTIC CASE
            T1  = GAMUPR(2*LQN+2,X)
            T2  = 3.0D0*X*GAMLWR(2*LQN+3,X)-GAMLWR(2*LQN+5,X)
            ULL = T1 + T2/FAC
            ULL = 0.5D0*ULL/E22
C
C           TRANSFER INTO NUCLEAR ATTRACTION MATRIX
            VMAT(IBAS,JBAS) = -ZNUC(IZ)*RN(M,1)*ULL
C
          ELSE
C
C           RELATIVISTIC CASE
            T1  = GAMUPR(2*LQN+2,X)
            T2  = 3.0D0*X*GAMLWR(2*LQN+3,X)-GAMLWR(2*LQN+5,X)
            ULL = T1 + T2/FAC
            ULL = 0.5D0*ULL/E22
C
            SSS = 2.0D0*GAMHLF(2*LQN+5)*EPR/E52
C
            T1  = GAMUPR(2*LQN+4,X)
            T2  = 3.0D0*X*GAMLWR(2*LQN+5,X)-GAMLWR(2*LQN+7,X)
            VS2 = T1 + T2/FAC
            VS2 = 2.0D0*EPR*VS2/E42
C
            IF(KQN.LT.0) THEN
C
              VS0 = 0.0D0
              VS1 = 0.0D0
C
            ELSE
C
              T1  = GAMUPR(2*LQN+2,X)
              T2  = 3.0D0*X*GAMLWR(2*LQN+3,X)-GAMLWR(2*LQN+5,X)
              VS1 = T1 + T2/FAC
              VS1 =-TL*VS1/E02
C
              T1  = GAMUPR(2*LQN  ,X)
              T2  = 3.0D0*X*GAMLWR(2*LQN+1,X)-GAMLWR(2*LQN+3,X)
              VS0 = T1 + T2/FAC
              VS0 = 0.5D0*TL*TL*VS0/E02
C
            ENDIF
            USS = VS0+VS1+VS2
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NBAS
            LBAS = JBAS+NBAS
C
C           TRANSFER INTO NUCLEAR ATTRACTION MATRIX
            VMAT(IBAS,JBAS) = -ZNUC(IZ)*RN(M,1)*ULL
            VMAT(KBAS,LBAS) = -ZNUC(IZ)*RN(M,4)*USS
     &                        -2.0D0*EMSS*CV*CV*RN(M,4)*SSS
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
      SUBROUTINE GSSNC0(VMAT,EXL,IZ,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         GGGGGG   SSSSSS   SSSSSS  NN    NN  CCCCCC   000000          C
C        GG    GG SS    SS SS    SS NNN   NN CC    CC 00   000         C
C        GG       SS       SS       NNNN  NN CC       00  0000         C
C        GG        SSSSSS   SSSSSS  NN NN NN CC       00 00 00         C
C        GG   GGG       SS       SS NN  NNNN CC       0000  00         C
C        GG    GG SS    SS SS    SS NN   NNN CC    CC 000   00         C
C         GGGGGG   SSSSSS   SSSSSS  NN    NN  CCCCCC   000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  GSSNC0 CALCULATES ATOM-CENTRED NUCLEAR ATTRACTION MATRIX ELEMENTS   C
C  FOR A GAUSSIAN NUCLEAR CHARGE MODEL.                                C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RN(MBS*MBS,4),VMAT(MBD,MBD),EXL(MBS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      RL  = DFLOAT(LQN)
      TL  = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     GAUSSIAN WIDTH PARAMETER
      GPRM = 1.5D0/(RNUC(IZ)*RNUC(IZ))
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
          T52 = RL+2.5D0
          E52 = EIJ**T52
C
C         MATRIX ELEMENT DEPENDS ON HAMILTONIAN
          IF(HMLT.EQ.'NORL') THEN
C
C           NON-RELATIVISTIC CASE
            ULL = ERFINT0(LQN  ,EIJ,GPRM)
C
C           TRANSFER INTO NUCLEAR ATTRACTION MATRIX
            VMAT(IBAS,JBAS) = -ZNUC(IZ)*RN(M,1)*ULL
C
          ELSE
C
C           RELATIVISTIC CASE
            SSS = 2.0D0*GAMHLF(2*LQN+5)*EPR/E52
            ULL = ERFINT0(LQN  ,EIJ,GPRM)
            VSA = 4.0D0*EPR*ERFINT0(LQN+1,EIJ,GPRM)
            IF(KQN.LT.0) THEN
              VSN = 0.0D0
            ELSE
              VSN =-2.0D0*EIJ*TL*ERFINT0(LQN  ,EIJ,GPRM)
     &                    +TL*TL*ERFINT0(LQN-1,EIJ,GPRM)
            ENDIF
            USS = VSA+VSN
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NBAS
            LBAS = JBAS+NBAS
C
C           TRANSFER INTO ARRAY
            VMAT(IBAS,JBAS) = -ZNUC(IZ)*RN(M,1)*ULL
            VMAT(KBAS,LBAS) = -ZNUC(IZ)*RN(M,4)*USS
     &                        -2.0D0*EMSS*CV*CV*RN(M,4)*SSS
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
      SUBROUTINE GSMNC0(VMAT,EXL,IZ,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        GGGGGG   SSSSSS  MM       MM NN    NN  CCCCCC   000000        C
C       GG    GG SS    SS MMM     MMM NNN   NN CC    CC 00   000       C
C       GG       SS       MMMM   MMMM NNNN  NN CC       00  0000       C
C       GG        SSSSSS  MM MM MM MM NN NN NN CC       00 00 00       C
C       GG   GGG       SS MM  MMM  MM NN  NNNN CC       0000  00       C
C       GG    GG SS    SS MM   M   MM NN   NNN CC    CC 000   00       C
C        GGGGGG   SSSSSS  MM       MM NN    NN  CCCCCC   000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  GSMNC0 CALCULATES ATOM-CENTRED NUCLEAR ATTRACTION MATRIX ELEMENTS   C
C  FOR A BASIS OF GAUSSIAN NUCLEAR CHARGES, FRACTIONALLY DISTRIBUTED.  C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RN(MBS*MBS,4),VMAT(MBD,MBD),EXL(MBS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      RL  = DFLOAT(LQN)
      TL  = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
          T52 = RL+2.5D0
          E52 = EIJ**T52
C
C         MATRIX ELEMENT DEPENDS ON HAMILTONIAN
          IF(HMLT.EQ.'NORL') THEN
C
C           NON-RELATIVISTIC CASE
            ULL =-FNUC(IZ,0)*ERFINT0(LQN  ,EIJ,XNUC(IZ,0))
            IF(NNUC(IZ).GT.0) THEN
              ULL = ULL + ZCVINT0(IZ,LQN  ,EIJ)
            ENDIF
C
C           TRANSFER INTO NUCLEAR ATTRACTION MATRIX
            VMAT(IBAS,JBAS) = -ZNUC(IZ)*RN(M,1)*ULL
C
          ELSE
C
C           RELATIVISTIC CASE
            SSS = 2.0D0*GAMHLF(2*LQN+5)*EPR/E52
C
            ULL =-FNUC(IZ,0)*ERFINT0(LQN  ,EIJ,XNUC(IZ,0))
            IF(NNUC(IZ).GT.0) THEN
              ULL = ULL + ZCVINT0(IZ,LQN  ,EIJ)
            ENDIF
C
C           NUCLEAR BEST-FIT EXPANSION
            VSA = 4.0D0*EPR*ERFINT0(LQN+1,EIJ,XNUC(IZ,0))
            IF(KQN.LT.0) THEN
              VSN = 0.0D0
            ELSE
              VSN =-2.0D0*EIJ*TL*ERFINT0(LQN  ,EIJ,XNUC(IZ,0))
     &                    +TL*TL*ERFINT0(LQN-1,EIJ,XNUC(IZ,0))
            ENDIF
            USS =-FNUC(IZ,0)*(VSA+VSN)
C
            IF(NNUC(IZ).GT.0) THEN
              VSA = 4.0D0*EPR*ZCVINT0(IZ,LQN+1,EIJ)
              IF(KQN.LT.0) THEN
                VSN = 0.0D0
              ELSE
                VSN =-2.0D0*EIJ*TL*ZCVINT0(IZ,LQN  ,EIJ)
     &                      +TL*TL*ZCVINT0(IZ,LQN-1,EIJ)
              ENDIF
              USS = USS + VSA+VSN
            ENDIF
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NBAS
            LBAS = JBAS+NBAS
C
C           TRANSFER INTO NUCLEAR ATTRACTION MATRIX
            VMAT(IBAS,JBAS) = ZNUC(IZ)*RN(M,1)*ULL
            VMAT(KBAS,LBAS) = ZNUC(IZ)*RN(M,4)*USS
     &                        -2.0D0*EMSS*CV*CV*RN(M,4)*SSS
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
      SUBROUTINE FMINC0(VMAT,EXL,IZ,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         FFFFFFFF MM       MM IIII NN    NN  CCCCCC   000000          C
C         FF       MMM     MMM  II  NNN   NN CC    CC 00   000         C
C         FF       MMMM   MMMM  II  NNNN  NN CC       00  0000         C
C         FFFFFF   MM MM MM MM  II  NN NN NN CC       00 00 00         C
C         FF       MM  MMM  MM  II  NN  NNNN CC       0000  00         C
C         FF       MM   M   MM  II  NN   NNN CC    CC 000   00         C
C         FF       MM       MM IIII NN    NN  CCCCCC   000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  FMINC0 CALCULATES ATOM-CENTRED NUCLEAR ATTRACTION MATRIX ELEMENTS   C
C  FOR A FERMI NUCLEAR CHARGE MODEL.                                   C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RN(MBS*MBS,4),VMAT(MBD,MBD),EXL(MBS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      RL  = DFLOAT(LQN)
      TL  = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     TWO-PARAMETER FERMI DETAILS
      APRM = AFMI(IZ)
      CPRM = 5.0D0*RNUC(IZ)*RNUC(IZ)/3.0D0
     &     - 7.0D0*PI*PI*AFMI(IZ)*AFMI(IZ)/3.0D0
      CPRM = DSQRT(CPRM)
C
C     WARN USER IF FERMI MODEL IS NOT PHYSICAL
      IF(CPRM.LT.APRM) THEN
        WRITE(6, *) 'In FMINC0: nucleus is too small for this model.'
        WRITE(7, *) 'In FMINC0: nucleus is too small for this model.'
      ENDIF
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
          T52 = RL+2.5D0
          E52 = EIJ**T52
C
C         MATRIX ELEMENT DEPENDS ON HAMILTONIAN
          IF(HMLT.EQ.'NORL') THEN
C
C           NON-RELATIVISTIC CASE
            ULL = FMIINT0(LQN  ,EIJ,APRM,CPRM)
C
C           TRANSFER INTO NUCLEAR ATTRACTION MATRIX
            VMAT(IBAS,JBAS) = -ZNUC(IZ)*RN(M,1)*ULL
C
          ELSE
C
C           RELATIVISTIC CASE
            ULL = FMIINT0(LQN  ,EIJ,APRM,CPRM)
C
            SSS = 2.0D0*GAMHLF(2*LQN+5)*EPR/E52
C
            VSA = 4.0D0*EPR*FMIINT0(LQN+1,EIJ,APRM,CPRM)
            IF(KQN.LT.0) THEN
              VSN = 0.0D0
            ELSE
              VSN =-2.0D0*EIJ*TL*FMIINT0(LQN  ,EIJ,APRM,CPRM)
     &                    +TL*TL*FMIINT0(LQN-1,EIJ,APRM,CPRM)
            ENDIF
            USS = VSA+VSN
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NBAS
            LBAS = JBAS+NBAS
C
C           TRANSFER INTO NUCLEAR ATTRACTION MATRIX
            VMAT(IBAS,JBAS) = -ZNUC(IZ)*RN(M,1)*ULL
            VMAT(KBAS,LBAS) = -ZNUC(IZ)*RN(M,4)*USS
     &                        -2.0D0*EMSS*CV*CV*RN(M,4)*SSS
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
      SUBROUTINE FMINC0temp(VMAT,EXL,IZ,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         FFFFFFFF MM       MM IIII NN    NN  CCCCCC   000000          C
C         FF       MMM     MMM  II  NNN   NN CC    CC 00   000         C
C         FF       MMMM   MMMM  II  NNNN  NN CC       00  0000         C
C         FFFFFF   MM MM MM MM  II  NN NN NN CC       00 00 00         C
C         FF       MM  MMM  MM  II  NN  NNNN CC       0000  00         C
C         FF       MM   M   MM  II  NN   NNN CC    CC 000   00         C
C         FF       MM       MM IIII NN    NN  CCCCCC   000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  FMINC0Temp CALCULATES ATOM-CENTRED NUCLEAR ATTRACTION MATRIX        C
C  ELEMENTS FOR A FERMI NUCLEAR CHARGE MODEL BY STARTING WITH THE      C
C  ANALYTIC GAUSSIAN FORMULA AND THEN TAKING A QUADRATURE DIFFERENCE.  C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RN(MBS*MBS,4),VMAT(MBD,MBD),EXL(MBS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      RL  = DFLOAT(LQN)
      TL  = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     GAUSSIAN WIDTH PARAMETER
      GPRM = 1.5D0/(RNUC(IZ)*RNUC(IZ))
C
C     TWO-PARAMETER FERMI DETAILS
      TPRM = TFMI(IZ)
      CPRM = 5.0D0*RNUC(IZ)*RNUC(IZ)/3.0D0
     &     - 7.0D0*PI*PI*AFMI(IZ)*AFMI(IZ)/3.0D0
      CPRM = DSQRT(CPRM)
      RPRM = RNUC(IZ)
C
C     WARN USER IF FERMI MODEL IS NOT PHYSICAL
      IF(CPRM.LT.APRM) THEN
        WRITE(6, *) 'In FRMINC0: nucleus is too small for this model.'
        WRITE(7, *) 'In FRMINC0: nucleus is too small for this model.'
      ENDIF
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
          T52 = RL+2.5D0
          E52 = EIJ**T52
C
C         MATRIX ELEMENT DEPENDS ON HAMILTONIAN
          IF(HMLT.EQ.'NORL') THEN
C
C           NON-RELATIVISTIC CASE
            ULLG = ERFINT0(LQN  ,EIJ,GPRM)
            ULLD = FMIINT0temp(LQN  ,EIJ,TPRM,CPRM,RPRM)
C
C           TRANSFER INTO NUCLEAR ATTRACTION MATRIX
            VMAT(IBAS,JBAS) = -ZNUC(IZ)*RN(M,1)*(ULLG+ULLD)
C
          ELSE
C
C           RELATIVISTIC CASE
            SSS  = 2.0D0*GAMHLF(2*LQN+5)*EPR/E52
            ULLG = ERFINT0(LQN  ,EIJ,GPRM)
            ULLD = FMIINT0temp(LQN  ,EIJ,TPRM,CPRM,RPRM)
            VSAG = 4.0D0*EPR*ERFINT0(LQN+1,EIJ,GPRM)
            VSAD = 4.0D0*EPR*FMIINT0temp(LQN+1,EIJ,TPRM,CPRM,RPRM)
            IF(KQN.LT.0) THEN
              VSNG = 0.0D0
              VSND = 0.0D0
            ELSE
              VSNG =-2.0D0*EIJ*TL*ERFINT0(LQN  ,EIJ,GPRM)
     &                     +TL*TL*ERFINT0(LQN-1,EIJ,GPRM)
              VSND =-2.0D0*EIJ*TL*FMIINT0temp(LQN  ,EIJ,TPRM,CPRM,RPRM)
     &                     +TL*TL*FMIINT0temp(LQN-1,EIJ,TPRM,CPRM,RPRM)
            ENDIF
            USSG = VSAG+VSNG
            USSD = VSAD+VSND
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NBAS
            LBAS = JBAS+NBAS
C
C           TRANSFER INTO ARRAY
            VMAT(IBAS,JBAS) = -ZNUC(IZ)*RN(M,1)*(ULLG+ULLD)
            VMAT(KBAS,LBAS) = -ZNUC(IZ)*RN(M,4)*(USSG+USSD)
     &                        -2.0D0*EMSS*CV*CV*RN(M,4)*SSS
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
      SUBROUTINE NCOLP0(OMAT,EXL,IZ,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        NN    NN  CCCCCC   OOOOOO  LL       PPPPPPP   000000          C
C        NNN   NN CC    CC OO    OO LL       PP    PP 00   000         C
C        NNNN  NN CC       OO    OO LL       PP    PP 00  0000         C
C        NN NN NN CC       OO    OO LL       PP    PP 00 00 00         C
C        NN  NNNN CC       OO    OO LL       PPPPPPP  0000  00         C
C        NN   NNN CC    CC OO    OO LL       PP       000   00         C
C        NN    NN  CCCCCC   OOOOOO  LLLLLLLL PP        000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  NCOLP0 CALCULATES ATOM-CENTRED NUCLEAR OVERLAP MATRIX ELEMENTS FOR  C
C  A GAUSSIAN NUCLEAR CHARGE MODEL (NORMALISED CHARGE DENSITY).        C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RN(MBS*MBS,4),OMAT(MBD,MBD),EXL(MBS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     INITIALISE ARRAY
      DO IBAS=1,NBAS
        DO JBAS=1,NBAS
          OMAT(IBAS     ,JBAS     ) = 0.0D0
          OMAT(IBAS+NBAS,JBAS+NBAS) = 0.0D0
          OMAT(IBAS+NBAS,JBAS     ) = 0.0D0
          OMAT(IBAS     ,JBAS+NBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      RL  = DFLOAT(LQN)
      TL  = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
          T52 = RL+2.5D0
          E52 = EIJ**T52
C
C         MATRIX ELEMENT DEPENDS ON HAMILTONIAN
          IF(HMLT.EQ.'NORL') THEN
C
C           NON-RELATIVISTIC CASE
            OLL = FNUC(IZ,0)*GSSINT0(LQN  ,EIJ,XNUC(IZ,0))
            IF(NNUC(IZ).GT.0) THEN
              OLL = OLL + ZCDINT0(IZ,LQN  ,EIJ)
            ENDIF
C
C           TRANSFER INTO DIRAC MATRIX
            OMAT(IBAS,JBAS) = RN(M,1)*OLL
C
          ELSE
C
C           RELATIVISTIC CASE
C
            OLL = FNUC(IZ,0)*GSSINT0(LQN  ,EIJ,XNUC(IZ,0))
            IF(NNUC(IZ).GT.0) THEN
              OLL = OLL + ZCDINT0(IZ,LQN  ,EIJ)
            ENDIF

            VSA = 4.0D0*EPR*GSSINT0(LQN+1,EIJ,XNUC(IZ,0))
            IF(KQN.LT.0) THEN
              VSN = 0.0D0
            ELSE
              VSN = -2.0D0*EIJ*TL*GSSINT0(LQN  ,EIJ,XNUC(IZ,0))
     &                     +TL*TL*GSSINT0(LQN-1,EIJ,XNUC(IZ,0))
            ENDIF
            OSS = OSS + FNUC(IZ,0)*(VSA+VSN)
            IF(NNUC(IZ).GT.0) THEN
              VSA = 4.0D0*EPR*ZCDINT0(IZ,LQN+1,EIJ)
              IF(KQN.LT.0) THEN
                VSN = 0.0D0
              ELSE
                VSN = -2.0D0*EIJ*TL*ZCDINT0(IZ,LQN  ,EIJ)
     &                       +TL*TL*ZCDINT0(IZ,LQN-1,EIJ)
              ENDIF
              OSS = OSS + (VSA+VSN)
            ENDIF
C
C           TRANSFER INTO NUCLEAR OVERLAP MATRIX
            OMAT(IBAS     ,JBAS     ) = RN(M,1)*OLL
            OMAT(IBAS+NBAS,JBAS+NBAS) = RN(M,4)*OSS
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
      SUBROUTINE MOMNT0(RR,EXL,KQN,NBAS,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      MM       MM  OOOOOO  MM       MM NN    NN TTTTTTTT 000000       C
C      MMM     MMM OO    OO MMM     MMM NNN   NN    TT   00   000      C
C      MMMM   MMMM OO    OO MMMM   MMMM NNNN  NN    TT   00  0000      C
C      MM MM MM MM OO    OO MM MM MM MM NN NN NN    TT   00 00 00      C
C      MM  MMM  MM OO    OO MM  MMM  MM NN  NNNN    TT   0000  00      C
C      MM   M   MM OO    OO MM   M   MM NN   NNN    TT   000   00      C
C      MM       MM  OOOOOO  MM       MM NN    NN    TT    000000       C
C                                                                      C
C -------------------------------------------------------------------- C
C  MOMNT0 CALCULATES RADIAL MATRIX ELEMENTS <r^K>^{TT} FOR INTEGER K.  C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      DIMENSION RN(MBS*MBS,4),RR(MBD,MBD),EXL(MBS)
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      RL  = DFLOAT(LQN)
      RK  = DFLOAT(K)
      TL  = DFLOAT(2*LQN+1)
C
C     INITIALISE ARRAYS
      DO IBAS=1,NBAS
        DO JBAS=1,NBAS
          RR(IBAS     ,JBAS     ) = 0.0D0
          RR(IBAS     ,JBAS+NBAS) = 0.0D0
          RR(IBAS+NBAS,JBAS     ) = 0.0D0
          RR(IBAS+NBAS,JBAS+NBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     CONDITIONAL FOR K THAT PRODUCES DIVERGING INTEGRAL
      IF(2*LQN+K.LE.-3) THEN
        RETURN
      ENDIF
C
C     NORMALISATION FACTORS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        DO JBAS=1,NBAS
          M = M+1
C
C         EXPONENT COMBINATIONS
          EI  = EXL(IBAS)
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
C
C         POWERS AND FACTORS
          TK3 = RL+0.5D0*RK+1.5D0
          TK5 = RL+0.5D0*RK+2.5D0
          EK3 = EIJ**TK3
          EK5 = EIJ**TK5
C
          IF(KQN.LT.0) THEN
            TM1 = 0.0D0
          ELSE
            TM1 =-RK*TL*EIJ*EIJ + 4.0D0*EPR*(RL+0.5D0+0.5D0*RK)*TK3         
          ENDIF
C
C         PRODUCT OF EACH PART
          RLL = 0.5D0*RN(M,1)*GAMHLF(2*LQN+K+3)/EK3
          IF(HMLT.EQ.'NORL') THEN
            RSS = 0.0D0
          ELSEIF(KQN.LT.0) THEN
            RSS = 2.0D0*RN(M,4)*GAMHLF(2*LQN+K+5)*EPR/EK5
          ELSEIF(KQN.GT.0) THEN
            RSS = 0.5D0*RN(M,4)*GAMHLF(2*LQN+K+1)*TM1/EK5
          ENDIF
C
C         SMALL COMPONENT BLOCK ADDRESSES
          KBAS = IBAS+NBAS
          LBAS = JBAS+NBAS
C
C         MATRIX ELEMENTS
          RR(IBAS,JBAS) = RLL
          IF(HMLT.NE.'NORL') THEN
            RR(IBAS,LBAS) = 0.0D0
            RR(KBAS,JBAS) = 0.0D0
            RR(KBAS,LBAS) = RSS
          ENDIF
C
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE AMPLT0(R0,R2,R4,EXL,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           AA    MM       MM PPPPPPP  LL      TTTTTTTT 000000         C
C          AAAA   MMM     MMM PP    PP LL         TT   00   000        C
C         AA  AA  MMMM   MMMM PP    PP LL         TT   00  0000        C
C        AA    AA MM MM MM MM PP    PP LL         TT   00 00 00        C
C        AAAAAAAA MM  MMM  MM PPPPPPP  LL         TT   0000  00        C
C        AA    AA MM   M   MM PP       LL         TT   000   00        C
C        AA    AA MM       MM PP       LLLLLLLL   TT    000000         C
C                                                                      C
C -------------------------------------------------------------------- C
C  AMPLT0 EVALUATES RADIAL FUNCTION POWER SERIES AMPLITUDES.           C
C -------------------------------------------------------------------- C
C  P(r) = r^(LQN+1)[p_0 + p_2 r^2 + ... ]   COLLECT p_0, p_2           C
C  Q(r) = r^(LQN  )[q_0 + q_2 r^2 + ... ]                              C
C      OR r^(LQN+2)[q_0 + q_2 r^2 + ... ]   COLLECT q_0, q_2           C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      DIMENSION R0(MBD),R2(MBD),R4(MBD)
      DIMENSION RNL(MBS),RNS(MBS),EXL(MBS)
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     INITIALISE ARRAY
      DO IBAS=1,2*NBAS
        R0(IBAS) = 0.0D0
        R2(IBAS) = 0.0D0
      ENDDO
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      RL  = DFLOAT(LQN)
C
C     NORMALISATION CONSTANTS
      GA = TWLG-GAMLOG(2*LQN+3)
      GB = TWLG-GAMLOG(2*LQN+5)
      RA = RL+1.5D0
      RB = RL+0.5D0
      DO IBAS=1,NBAS
        ELOG      = DLOG(2.0D0*EXL(IBAS))
        RNL(IBAS) = DEXP(0.5D0*(GA+RA*ELOG))
        RNS(IBAS) = DEXP(0.5D0*(GB+RB*ELOG))
      ENDDO
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      DO IBAS=1,NBAS
C
C       EXPONENTS
        EI = EXL(IBAS)
C
C       LARGE-COMPONENT POWER SERIES AMPLITUDES
        RL0 = RNL(IBAS)
        RL2 =-RNL(IBAS)*EI
        RL4 = RNL(IBAS)*0.5D0*EI*EI
C
C       SMALL-COMPONENT POWER SERIES AMPLITUDES
        IF(HMLT.EQ.'NORL') THEN
          RS0 = RNL(IBAS)*(RL+3.0D0)
          RS2 =-RNL(IBAS)*(RL+5.0D0)*EI
          RS4 = RNL(IBAS)*(RL+7.0D0)*0.5D0*EI*EI
        ELSEIF(KQN.LT.0) THEN
          RS0 =-RNS(IBAS)*2.0D0*EI
          RS2 = RNS(IBAS)*2.0D0*EI*EI
          RS4 =-RNS(IBAS)*EI*EI*EI
        ELSEIF(KQN.GT.0) THEN
          RS0 = RNS(IBAS)*(2.0D0*RL+1.0D0)
          RS2 =-RNS(IBAS)*(2.0D0*RL+3.0D0)*EI
          RS4 =-RNS(IBAS)*(2.0D0*RL+5.0D0)*0.5D0*EI*EI
        ENDIF
C
C       MATRIX ELEMENTS
        R0(IBAS) = RL0
        R2(IBAS) = RL2
        R4(IBAS) = RL4

        R0(IBAS+NBAS) = RS0
        R2(IBAS+NBAS) = RS2
        R4(IBAS+NBAS) = RS4
C
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ANMLS0(AMAT,EXL,IZ,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          AA    NN    NN MM       MM LL       SSSSSS   000000         C
C         AAAA   NNN   NN MMM     MMM LL      SS    SS 00   000        C
C        AA  AA  NNNN  NN MMMM   MMMM LL      SS       00  0000        C
C       AA    AA NN NN NN MM MM MM MM LL       SSSSSS  00 00 00        C
C       AAAAAAAA NN  NNNN MM  MMM  MM LL            SS 0000  00        C
C       AA    AA NN   NNN MM   M   MM LL      SS    SS 000   00        C
C       AA    AA NN    NN MM       MM LLLLLLLL SSSSSS   000000         C
C                                                                      C
C -------------------------------------------------------------------- C
C  ANMLS0 GENERATES ATOMIC ELECTRON MAGNETIC MOMENT INTERACTION MATRIX C
C  ELEMENTS FOR SYMMETRY TYPE KQNA.                                    C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RTT(MBS*MBS,4),AMAT(MBD,MBD),EXL(MBS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/FCTS/RFACT(0:80),SFACT(0:80)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     INITIALISE ARRAY
      DO IBAS=1,2*NBAS
        DO JBAS=1,2*NBAS
          AMAT(IBAS,JBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     AMPLITUDE FOR ELECTRON ANOMALOUS MAGNETIC MOMENT TERM (B&S 19.3)
      AG1 =-1.0D0/(4.0D0*PI*EMSS*CV*CV)
C
C     INTEGRATION GRID DETAILS (WHEN QUADRATURE IS INVOKED)
      NLIN = 100
      NEXP = 5000
      RMAX = 25.0D0
      ELIN = RNUC(IZ)/DFLOAT(NLIN)
      EEXP = DLOG(RMAX/RNUC(IZ))/DFLOAT(NEXP)
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      TL  = DFLOAT(LQN+KQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RTT,EXL,NBAS,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
C
C         INITIALISE RAW MATRIX ELEMENTS
          ULS = 0.0D0
          USL = 0.0D0
C
C         RAW MATRIX ELEMENTS FOR A POINT-NUCLEUS
          IF(NMDL(IZ).EQ.'POINT') THEN
C
            IF(LQN.EQ.0) THEN
              ULS =-EJ/EIJ
              USL =-EI/EIJ
            ELSE
              AFC = 0.5D0*RFACT(LQN-1)/(EIJ**(LQN+1))
              ULS = AFC*((KQN+LQN+1)*EI+(KQN-LQN+1)*EJ)
              USL = AFC*((KQN-LQN+1)*EI+(KQN+LQN+1)*EJ)
            ENDIF
C
C         RAW MATRIX ELEMENTS FOR A GAUSSIAN (OR BEST-FIT) NUCLEUS
          ELSEIF(NMDL(IZ).EQ.'GAUSS'.OR.NMDL(IZ).EQ.'FERMI') THEN
C
            XI  = XNUC(IZ,0)
            FC  = FNUC(IZ,0)
C
            XRT =-2.0D0*DSQRT(XI)/PI12
            PES = XI+EIJ
            P12 = DSQRT(PES)
            P32 = PES*P12
            XEL = PES**LQN
C
            ULS = ULS - FC*XRT*EJ*GAMHLF(2*LQN+3)/(XEL*P32)
            USL = USL - FC*XRT*EI*GAMHLF(2*LQN+3)/(XEL*P32)
            IF(KQN.GT.0) THEN
              ULS = ULS + 0.5D0*FC*XRT*TL*GAMHLF(2*LQN+1)/(XEL*P12)
              USL = USL + 0.5D0*FC*XRT*TL*GAMHLF(2*LQN+1)/(XEL*P12)
            ENDIF
            ULS = ULS - 2.0D0*FC*EJ*ERFINT0(LQN  ,EIJ,XI)
            USL = USL - 2.0D0*FC*EI*ERFINT0(LQN  ,EIJ,XI)
            IF(KQN.GT.0) THEN
              ULS = ULS + FC*TL*ERFINT0(LQN-1,EIJ,XI)
              USL = USL + FC*TL*ERFINT0(LQN-1,EIJ,XI)
            ENDIF
            
            IF(NNUC(IZ).GT.0) THEN
              ULS = ULS + 2.0D0*TL*ZCEINT0(IZ,LQN  ,EIJ)
     &                  - 4.0D0*EJ*ZCEINT0(IZ,LQN+1,EIJ)
              USL = USL + 2.0D0*TL*ZCEINT0(IZ,LQN  ,EIJ)
     &                  - 4.0D0*EI*ZCEINT0(IZ,LQN+1,EIJ)
            ENDIF
C
C         RAW MATRIX ELEMENTS FOR A GENERIC NUCLEUS (BY QUADRATURE)
          ELSE
C
C           USE A LINEAR GRID TO INTEGRATE FROM 0 TO RNUC
            DO N=0,NLIN
              RN  = ELIN*DFLOAT(N)
              IF(KQN.LT.0) THEN
                Z1  = RN**(2*LQN+1)
                ZLS =-2.0D0*EJ
                ZSL =-2.0D0*EI
              ELSEIF(KQN.GT.0) THEN
                Z1  = RN**(2*LQN-1)
                ZLS = TL-2.0D0*EJ*RN*RN
                ZSL = TL-2.0D0*EI*RN*RN
              ENDIF
              Z2  = DEXP(-EIJ*RN*RN)
              Z3  = ERNUCDNS(NMDL(IZ),IZ,RN)
              ULS = ULS - ELIN*EXTINT11(Z1*Z2*Z3*ZLS,N,NLIN)
              USL = USL - ELIN*EXTINT11(Z1*Z2*Z3*ZSL,N,NLIN)
            ENDDO
C
C           USE AN EXP GRID TO INTEGRATE FROM RNUC TO RMAX
            DO N=0,NEXP
              TI  = EEXP*DFLOAT(N)
              RN  = RNUC(IZ)*DEXP(TI)
              Z1  = RN**(2*LQN)
              Z2  = DEXP(-EIJ*RN*RN)
              Z3  = ERNUCDNS(NMDL(IZ),IZ,RN)
              ZLS = TL-2.0D0*EJ*RN*RN
              ZSL = TL-2.0D0*EI*RN*RN
              ULS = ULS - EEXP*EXTINT11(Z1*Z2*Z3*ZLS,N,NEXP)
              USL = USL - EEXP*EXTINT11(Z1*Z2*Z3*ZSL,N,NEXP)
            ENDDO
C
C           INCLUDE NEWTON-COTES MULTIPLICATIVE FACTOR
            ULS = 5.0D0*ULS/299376.0D0
            USL = 5.0D0*USL/299376.0D0
C
          ENDIF
C
C         SKIP POINT FOR ULTRAVIOLET CUTOFF PARAMETER
101       CONTINUE
C
C         TRANSFER INTO ARRAY
          AMAT(IBAS     ,JBAS+NBAS) = ZNUC(IZ)*RTT(M,2)*AG1*ULS
          AMAT(IBAS+NBAS,JBAS     ) = ZNUC(IZ)*RTT(M,3)*AG1*USL
C
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE SLFLW0(YMAT,IZ,KQNA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       SSSSSS  LL       FFFFFFFF LL       WW        WW  000000        C
C      SS    SS LL       FF       LL       WW        WW 00   000       C
C      SS       LL       FF       LL       WW   WW   WW 00  0000       C
C       SSSSSS  LL       FFFFFF   LL       WW  WWWW  WW 00 00 00       C
C            SS LL       FF       LL       WW WW  WW WW 0000  00       C
C      SS    SS LL       FF       LL       WWWW    WWWW 000   00       C
C       SSSSSS  LLLLLLLL FF       LLLLLLLL WW        WW  000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  SLFLW0 GENERATES ATOMIC ELECTRON SELF-INTERACTION (LOW-ENERGY)      C
C  MATRIX ELEMENTS FOR SYMMETRY TYPE KQNA USING THE BETHE FORMULATION. C
C -------------------------------------------------------------------- C
C  ▶ THE `LL' COMPONENT OF THIS CAN BE USED IN THE 'NORL' OPTION BUT   C
C    THE INPUT DECK ASSUMES THAT 'DHFQ' MEANS RELATIVISTIC.            C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      COMPLEX*16 SPINDOT
      COMPLEX*16 COEF(MDM,MDM)
C
      DIMENSION YMAT(MBD,MBD)
      DIMENSION RNLL(MBS,MBS),RNLS(MBS,MBS),RNSL(MBS,MBS),RNSS(MBS,MBS)
      DIMENSION EX2(MBS,2),LQ2(2),NBS2(2)
      DIMENSION PIJLS(MBS,MBS),PIJSL(MBS,MBS)
      DIMENSION PINLS(MBS,MBS),PINSL(MBS,MBS)
      DIMENSION PMNLS(MBS,MBS),PMNSL(MBS,MBS)
      DIMENSION OPRJ(MBD,MBS),ELWM(MBS)
      DIMENSION OLL(MBS,MBS),OSS(MBS,MBS)
C
      COMMON/AUFB/NORB(0:MEL,MKP+1),NUMOCC(0:MEL),LMXCONF
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
      DATA DEPS/1.0D-2/
C
C     LOW- AND HIGH-ENERGY PARTITION CUTOFF
      CUTK = CV
C     CUTK = CV*ZNUC(IZ)*ZNUC(IZ)
C
C     AMPLITUDE FOR LOW-ENERGY CONTRIBUTION
      ALW =-2.0D0/(3.0D0*PI*CV)
C
C     DETERMINE THE ORBITAL VALUE LQNA
      LQNA = LVAL(KQNA)
      TLA  = DFLOAT(LQNA+KQNA+1)
C
C     NUMBER OF BASIS FUNCTIONS
      NBASA = NFNC(LQNA,IZ)
      IF(KQNA.LT.0) THEN
        KA =-2*KQNA-1
      ELSE
        KA = 2*KQNA
      ENDIF
C
C     INITIALISE LIST OF ENERGY VALUES
      DO IOCC=1,NUMOCC(LQNA)
        IF(NORB(LQNA,IOCC).NE.0) THEN
          ELWM(IOCC) = 0.0D0
        ENDIF
      ENDDO
C
C     PUT CURRENT BLOCK DETAILS IN SOME ARRAYS
      LQ2(1)  = LQNA
      NBS2(1) = NBASA
      DO IBAS=1,NBASA
        EX2(IBAS,1) = BEXL(IBAS,LQNA,IZ)
      ENDDO
C
C     LOOP OVER SYMMETRY TYPES FOR THE VIRTUAL STATES
      DO KB=1,NKAP(IZ)
C
C       DETERMINE THE ORBITAL VALUE LQNB
        KQNB = KAPA(KB,IZ)
        LQNB = LVAL(KQNB)
        TLB  = DFLOAT(LQNB+KQNB+1)
C
C       NUMBER OF BASIS FUNCTIONS
        NBASB = NFNC(LQNB,IZ)
C
C       PUT CURRENT BLOCK DETAILS INTO SOME ARRAYS
        LQ2(2)  = LQNB
        NBS2(2) = NBASB
        DO JBAS=1,NBASB
          EX2(JBAS,2) = BEXL(JBAS,LQNB,IZ)
        ENDDO
C
C       SELECTION RULE: |LA-LB|=1
        IF(IABS(LQNA-LQNB).NE.1) GOTO 101
C       
C       ANGULAR FACTORS OVER THIS KQNA/KQNB BLOCK (WILL BE REAL)
        ALL = 0.0D0
        ALS = 0.0D0
        ASL = 0.0D0
        ASS = 0.0D0
        DO MQNA=-2*IABS(KQNA)+1,2*IABS(KQNA)-1,2
          DO MQNB=-2*IABS(KQNB)+1,2*IABS(KQNB)-1,2
            ALL = ALL + DREAL(SPINDOT(2,3,KQNA,KQNB,MQNA,MQNB))
            ALS = ALS - DREAL(SPINDOT(2,2,KQNA,KQNB,MQNA,MQNB))
            ASL = ASL - DREAL(SPINDOT(3,3,KQNA,KQNB,MQNA,MQNB))
            ASS = ASS + DREAL(SPINDOT(3,2,KQNA,KQNB,MQNA,MQNB))
          ENDDO
        ENDDO
C
C       SKIP OPTION (ALL ANGULAR FACTORS ARE ZERO)
        IF(DABS(ALL)+DABS(ALS)+DABS(ASL)+DABS(ASS).LT.1.0D-8) GOTO 101
C
C       RNTT NORMALISATION CONSTANTS FOR THIS BLOCK
        CALL RNORM1(RNLS,EX2,LQ2,NBS2,2)
        CALL RNORM1(RNSL,EX2,LQ2,NBS2,3)
C
C       A BLOCK OF RADIAL INTEGRALS
        LAB = LQNA+LQNB
        GAM = GAMHLF(LAB+2)
        DO IBAS=1,NBASA
          EI = EX2(IBAS,1)
          DO JBAS=1,NBASB
            EJ = EX2(JBAS,2)
C
C           GAUSSIAN PARAMETER COMBINATIONS
            EIJ = EI+EJ
            EP4 = EIJ**(0.5D0*LAB+2.0D0)
C
C           NORMALISED MATRIX ELEMENTS
            PLS = 0.5D0*GAM*(TLB*EIJ-EJ*(LAB+2.0D0))/EP4
            PSL = 0.5D0*GAM*(TLA*EIJ-EI*(LAB+2.0D0))/EP4
            PIJLS(IBAS,JBAS) = RNLS(IBAS,JBAS)*PLS
            PIJSL(IBAS,JBAS) = RNSL(IBAS,JBAS)*PSL
          ENDDO
        ENDDO
C
C       INITIALISE SOME NEW ARRAYS
        DO IBAS=1,NBASA
          DO NPSV=1,NBASB
            PINLS(IBAS,NPSV) = 0.0D0
            PINSL(IBAS,NPSV) = 0.0D0
          ENDDO
        ENDDO
C
C       CONTRACT THE RADIAL INTEGRALS OVER POSITIVE-ENERGY K' STATES
        DO NPSV=1,NBASB
C
C         STARTING ADDRESS FOR THIS BLOCK
          NL = LRGE(IZ,KB,1)
          NS = LRGE(IZ,KB,1)+NSKP
C
C         LOOP OVER BASIS FUNCTIONS IBAS OF BLOCK K
          DO IBAS=1,NBASA
C
C           CONTRACT OVER INDEX JBAS
            DO JBAS=1,NBASB
C
C             COUPLE PLS ELEMENTS WITH SMALL-COMPONENT COEFFICIENTS
              PINLS(IBAS,NPSV) = PINLS(IBAS,NPSV)
     &                  + DREAL(COEF(NS+JBAS,NS+NPSV))*PIJLS(IBAS,JBAS)
C
C             COUPLE PSL ELEMENTS WITH LARGE-COMPONENT COEFFICIENTS
              PINSL(IBAS,NPSV) = PINSL(IBAS,NPSV)
     &                  + DREAL(COEF(NL+JBAS,NS+NPSV))*PIJSL(IBAS,JBAS)
C            
            ENDDO
C
C         END LOOP OVER BASIS FUNCTIONS ON BLOCK K          
          ENDDO
C
C       END LOOP OVER POSITIVE-ENERGY STATES
        ENDDO
C
C       INITIALISE SOME NEW ARRAYS
        DO IOCC=1,NUMOCC(LQNA)
          IF(NORB(LQNA,IOCC).NE.0) THEN
            DO NPSV=1,NBASB
              PMNLS(IOCC,NPSV) = 0.0D0
              PMNSL(IOCC,NPSV) = 0.0D0
            ENDDO
          ENDIF
        ENDDO
C
C       CONTRACT THE RADIAL INTEGRALS OVER POSITIVE-ENERGY K STATES
        DO IOCC=1,NUMOCC(LQNA)
          IF(NORB(LQNA,IOCC).NE.0) THEN
C
C           SUBSHELL OCCUPANCY FACTOR
            SBFC = 1.0D0/DSQRT(DFLOAT(2*IABS(KQNA)))
C
C           STARTING ADDRESS FOR THIS BLOCK
            IL = LRGE(IZ,KA,1)
            IS = LRGE(IZ,KA,1)+NSKP
C
C           LOOP OVER POSITIVE-ENERGY STATES OF BLOCK K'
            DO NPSV=1,NBASB
C
C             CONTRACT OVER INDEX JBAS
              DO IBAS=1,NBASA
C
C               COUPLE PLS ELEMENTS WITH SMALL-COMPONENT COEFFICIENTS
                PMNLS(IOCC,NPSV) = PMNLS(IOCC,NPSV)
     &             + SBFC*DREAL(COEF(IL+IBAS,IS+IOCC))*PINLS(IBAS,NPSV)
C
C               COUPLE PSL ELEMENTS WITH LARGE-COMPONENT COEFFICIENTS
                PMNSL(IOCC,NPSV) = PMNSL(IOCC,NPSV)
     &             + SBFC*DREAL(COEF(IS+IBAS,IS+IOCC))*PINSL(IBAS,NPSV)
C            
              ENDDO
C
C           END LOOP OVER BASIS FUNCTIONS ON BLOCK K          
            ENDDO
C
C         END LOOP OVER OCCUPIED STATES
          ENDIF
        ENDDO
C
C
C       MATRIX ELEMENTS  <IOCC|Γ|IOCC>
        DO IOCC=1,NUMOCC(LQNA)
          IF(NORB(LQNA,IOCC).NE.0) THEN
C
C           STARTING ADDRESS FOR THIS BLOCK
            IS = LRGE(IZ,KA,1)+NSKP
C
C           REFERENCE ENERGY
            EI = EIGN(IS+IOCC)
C
            DO NPSV=1,NBASB
C
C             STARTING ADDRESS FOR THIS BLOCK
              NS = LRGE(IZ,KB,1)+NSKP
C
C             STARTING ADDRESS IN FOCK MATRIX AND ENERGY LEVEL
              EN = EIGN(NS+NPSV)
C
C             ENERGY DIFFERENCE
              DEIN = EI-EN
C
C             MANUALLY SKIP CASES OF DEGENERATE ENERGY LEVELS
              IF(DABS(DEIN).GT.DEPS) THEN
C
C               LOGARITHMIC ENERGY TERM
                DTRM = (CV*CUTK-DEIN)/DABS(DEIN)
                DTRM = DEIN*DLOG(DTRM)
C
C               VECTOR DOT PRODUCT OF TRANSITION CURRENTS
                DOT = ALL*PMNLS(IOCC,NPSV)*PMNLS(IOCC,NPSV)
     &              + ALS*PMNLS(IOCC,NPSV)*PMNSL(IOCC,NPSV)
     &              + ASL*PMNSL(IOCC,NPSV)*PMNLS(IOCC,NPSV)
     &              + ASS*PMNSL(IOCC,NPSV)*PMNSL(IOCC,NPSV)
C
C               MATRIX ELEMENT
                ELWM(IOCC) = ELWM(IOCC) + ALW*DOT*DTRM
C
              ENDIF
C
            ENDDO
          ENDIF
        ENDDO
C
C       SKIP POINT FOR ANGULAR SELECTION RULE
101     CONTINUE
C
C     END LOOP OVER KQNB
      ENDDO
C
C     WE NOW HAVE OUR LIST OF <IOCC|Γ|IOCC>
C
C**********************************************************************C
C     APPLY PROJECTION OPERATOR TO GENERATE VSLF (MEAN-FIELD)          C
C**********************************************************************C
C
C     PUT CURRENT BLOCK DETAILS IN SOME ARRAYS
      LQ2(2)  = LQNA
      NBS2(2) = NBASA
      DO JBAS=1,NBASA
        EX2(JBAS,2) = BEXL(JBAS,LQNA,IZ)
      ENDDO
C
C     RNTT NORMALISATION CONSTANTS FOR THIS BLOCK
      CALL RNORM1(RNLL,EX2,LQ2,NBS2,1)
      CALL RNORM1(RNSS,EX2,LQ2,NBS2,4)
C
C     GENERATE OVERLAP MATRICES
      GLL = GAMHLF(2*LQNA+3)
      GSS = GAMHLF(2*LQNA+5)
      DO IBAS=1,NBASA
        EI = EX2(IBAS,1)
        DO JBAS=1,NBASA
          EJ  = EX2(JBAS,2)
          EIJ = EI+EJ
          EPR = EI*EJ
          T32 = DFLOAT(LQNA)+1.5D0
          T52 = DFLOAT(LQNA)+2.5D0
          E32 = EIJ**T32
          E52 = EIJ**T52
          OLL(IBAS,JBAS) = 0.5D0*RNLL(IBAS,JBAS)*GLL/E32
          OSS(IBAS,JBAS) = 2.0D0*RNSS(IBAS,JBAS)*GSS*EPR/E52
        ENDDO
      ENDDO
C
C     STARTING ADDRESS FOR THIS BLOCK
      IL = LRGE(IZ,KA,1)
      IS = LRGE(IZ,KA,1)+NSKP
C
C     LOOP OVER BASIS FUNCTIONS IN THIS LQNA BLOCK
      DO IBAS=1,NBASA
C       LOOP OVER ALL OCCUPIED SUBSHELLS OF THIS KQN TYPE GIVEN MQN
        DO IOCC=1,NUMOCC(LQNA)
          IF(NORB(LQNA,IOCC).NE.0) THEN
            OPRJ(IBAS      ,IOCC) = 0.0D0
            OPRJ(IBAS+NBASA,IOCC) = 0.0D0
            DO JBAS=1,NBASA
              OPRJ(IBAS      ,IOCC) = OPRJ(IBAS      ,IOCC)
     &                        + COEF(IL+JBAS,IS+IOCC)*OLL(IBAS,JBAS)
              OPRJ(IBAS+NBASA,IOCC) = OPRJ(IBAS+NBASA,IOCC)
     &                        + COEF(IS+JBAS,IS+IOCC)*OSS(IBAS,JBAS)
            ENDDO
          ENDIF
        ENDDO
C
      ENDDO
C
C     GENERATE THE MEAN-FIELD-COMPATIBLE MATRIX ELEMENTS
      DO IBAS=1,2*NBASA
        DO JBAS=1,2*NBASA
          YMAT(IBAS,JBAS) = 0.0D0
          DO IOCC=1,NUMOCC(LQNA)
            IF(NORB(LQNA,IOCC).NE.0) THEN
              YMAT(IBAS,JBAS) = YMAT(IBAS,JBAS)
     &                     + OPRJ(IBAS,IOCC)*ELWM(IOCC)*OPRJ(JBAS,IOCC)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE SLFLWN(YMAT,IZ,KQNA,MOCC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       SSSSSS  LL       FFFFFFFF LL       WW        WW NN    NN       C
C      SS    SS LL       FF       LL       WW        WW NNN   NN       C
C      SS       LL       FF       LL       WW   WW   WW NNNN  NN       C
C       SSSSSS  LL       FFFFFF   LL       WW  WWWW  WW NN NN NN       C
C            SS LL       FF       LL       WW WW  WW WW NN  NNNN       C
C      SS    SS LL       FF       LL       WWWW    WWWW NN   NNN       C
C       SSSSSS  LLLLLLLL FF       LLLLLLLL WW        WW NN    NN       C
C                                                                      C
C -------------------------------------------------------------------- C
C  SLFLWN GENERATES ATOMIC ELECTRON SELF-INTERACTION (LOW-ENERGY)      C
C  MATRIX ELEMENTS FOR SYMMETRY TYPE KQNA USING THE BETHE FORMULATION. C
C  IT DOES *NOT* PROJECT -- THIS IS FOR OCCUPIED ORBITAL MOCC ONLY.    C
C -------------------------------------------------------------------- C
C  ▶ THE `LL' COMPONENT OF THIS CAN BE USED IN THE 'NORL' OPTION BUT   C
C    THE INPUT DECK ASSUMES THAT 'DHFQ' MEANS RELATIVISTIC.            C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      COMPLEX*16 SPINDOT
      COMPLEX*16 COEF(MDM,MDM)
C
      DIMENSION RKNT(MBS)
      DIMENSION YMAT(MBD,MBD)
      DIMENSION RNLS(MBS,MBS),RNSL(MBS,MBS)
      DIMENSION EX2(MBS,2),LQ2(2),NBS2(2)
      DIMENSION PIJLS(MBS,MBS),PIJSL(MBS,MBS)
      DIMENSION PINLS(MBS,MBS),PINSL(MBS,MBS)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
      DATA DEPS/1.0D-2/
C
C     LOW- AND HIGH-ENERGY PARTITION CUTOFF
      CUTK = CV
C     CUTK = CV*ZNUC(IZ)*ZNUC(IZ)
C
C     AMPLITUDE FOR LOW-ENERGY CONTRIBUTION
      ALW =-2.0D0/(3.0D0*PI*CV)
C
C     DETERMINE THE ORBITAL VALUE LQNA
      LQNA = LVAL(KQNA)
      TLA  = DFLOAT(LQNA+KQNA+1)
C
C     NUMBER OF BASIS FUNCTIONS
      NBASA = NFNC(LQNA,IZ)
      IF(KQNA.LT.0) THEN
        KA =-2*KQNA-1
      ELSE
        KA = 2*KQNA
      ENDIF
C
C     INITIALISE MATRIX
      DO IBAS=1,2*NBASA
        DO JBAS=1,2*NBASA
          YMAT(IBAS,JBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     PUT CURRENT BLOCK DETAILS IN SOME ARRAYS
      LQ2(1)  = LQNA
      NBS2(1) = NBASA
      DO IBAS=1,NBASA
        EX2(IBAS,1) = BEXL(IBAS,LQNA,IZ)
      ENDDO
C
C     LOOP OVER SYMMETRY TYPES FOR THE VIRTUAL STATES
      DO KB=1,NKAP(IZ)
C
C       DETERMINE THE ORBITAL VALUE LQNB
        KQNB = KAPA(KB,IZ)
        LQNB = LVAL(KQNB)
        TLB  = DFLOAT(LQNB+KQNB+1)
C
C       NUMBER OF BASIS FUNCTIONS
        NBASB = NFNC(LQNB,IZ)
C
C       PUT CURRENT BLOCK DETAILS INTO SOME ARRAYS
        LQ2(2)  = LQNB
        NBS2(2) = NBASB
        DO JBAS=1,NBASB
          EX2(JBAS,2) = BEXL(JBAS,LQNB,IZ)
        ENDDO
C
C       SELECTION RULE: |LA-LB|=1
        IF(IABS(LQNA-LQNB).NE.1) GOTO 101
C       
C       ANGULAR FACTORS OVER THIS KQNA/KQNB BLOCK (WILL BE REAL)
        ALL = 0.0D0
        ALS = 0.0D0
        ASL = 0.0D0
        ASS = 0.0D0
        DO MQNA=-2*IABS(KQNA)+1,2*IABS(KQNA)-1,2
          DO MQNB=-2*IABS(KQNB)+1,2*IABS(KQNB)-1,2
            ALL = ALL + DREAL(SPINDOT(2,3,KQNA,KQNB,MQNA,MQNB))
            ALS = ALS - DREAL(SPINDOT(2,2,KQNA,KQNB,MQNA,MQNB))
            ASL = ASL - DREAL(SPINDOT(3,3,KQNA,KQNB,MQNA,MQNB))
            ASS = ASS + DREAL(SPINDOT(3,2,KQNA,KQNB,MQNA,MQNB))
          ENDDO
        ENDDO
C
C       SERIES OF ENERGY DIFFERENCES
        EM = EIGN(MOCC)
        DO NPSV=1,NBASB
          NS   = LRGE(IZ,KB,1)+NSKP
          EDIF = EM-EIGN(NS+NPSV)
          ARGL = (CUTK*CV-EDIF)/DABS(EDIF)
          IF(DABS(EDIF).GT.DEPS) THEN
            RKNT(NPSV) = EDIF*DLOG(ARGL)
          ELSE
            RKNT(NPSV) = 0.0D0
          ENDIF
        ENDDO
C
C       SKIP OPTION (ALL ANGULAR FACTORS ARE ZERO)
        IF(DABS(ALL)+DABS(ALS)+DABS(ASL)+DABS(ASS).LT.1.0D-8) GOTO 101
C
C       RNTT NORMALISATION CONSTANTS FOR THIS BLOCK
        CALL RNORM1(RNLS,EX2,LQ2,NBS2,2)
        CALL RNORM1(RNSL,EX2,LQ2,NBS2,3)
C
C       A BLOCK OF RADIAL INTEGRALS
        LAB = LQNA+LQNB
        GAM = GAMHLF(LAB+2)
        DO IBAS=1,NBASA
          EI = EX2(IBAS,1)
          DO JBAS=1,NBASB
            EJ = EX2(JBAS,2)
C
C           GAUSSIAN PARAMETER COMBINATIONS
            EIJ = EI+EJ
            EP4 = EIJ**(0.5D0*LAB+2.0D0)
C
C           NORMALISED MATRIX ELEMENTS
            PLS = 0.5D0*GAM*(TLB*EIJ-EJ*(LAB+2.0D0))/EP4
            PSL = 0.5D0*GAM*(TLA*EIJ-EI*(LAB+2.0D0))/EP4
            PIJLS(IBAS,JBAS) = RNLS(IBAS,JBAS)*PLS
            PIJSL(IBAS,JBAS) = RNSL(IBAS,JBAS)*PSL
          ENDDO
        ENDDO
C
C       INITIALISE SOME NEW ARRAYS
        DO IBAS=1,NBASA
          DO NPSV=1,NBASB
            PINLS(IBAS,NPSV) = 0.0D0
            PINSL(IBAS,NPSV) = 0.0D0
          ENDDO
        ENDDO
C
C       CONTRACT THE RADIAL INTEGRALS OVER POSITIVE-ENERGY K' STATES
        DO NPSV=1,NBASB
C
C         STARTING ADDRESS FOR THIS BLOCK
          NL = LRGE(IZ,KB,1)
          NS = LRGE(IZ,KB,1)+NSKP
C
C         LOOP OVER BASIS FUNCTIONS IBAS OF BLOCK K
          DO IBAS=1,NBASA
C
C           CONTRACT OVER INDEX JBAS
            DO JBAS=1,NBASB
C
C             COUPLE PLS ELEMENTS WITH SMALL-COMPONENT COEFFICIENTS
              PINLS(IBAS,NPSV) = PINLS(IBAS,NPSV)
     &                  + DREAL(COEF(NS+JBAS,NS+NPSV))*PIJLS(IBAS,JBAS)
C
C             COUPLE PSL ELEMENTS WITH LARGE-COMPONENT COEFFICIENTS
              PINSL(IBAS,NPSV) = PINSL(IBAS,NPSV)
     &                  + DREAL(COEF(NL+JBAS,NS+NPSV))*PIJSL(IBAS,JBAS)
C            
            ENDDO
C
C         END LOOP OVER BASIS FUNCTIONS ON BLOCK K          
          ENDDO
C
C       END LOOP OVER POSITIVE-ENERGY STATES
        ENDDO
C
C       SUBSHELL OCCUPANCY FACTOR
        SBFC = 1.0D0/DFLOAT(2*IABS(KQNA))
C
C       FORM THE DOT PRODUCT OF MATRICES AND MULTIPLY BY ENERGY TERM
        DO IBAS=1,NBASA
          DO JBAS=1,NBASA
            DO NPSV=1,NBASB
C
              YMAT(IBAS      ,JBAS      ) = YMAT(IBAS      ,JBAS      )
     &      + SBFC*ALL*ALW*RKNT(NPSV)*PINLS(IBAS,NPSV)*PINLS(JBAS,NPSV)
C
              YMAT(IBAS      ,JBAS+NBASA) = YMAT(IBAS      ,JBAS+NBASA)
     &      + SBFC*ALS*ALW*RKNT(NPSV)*PINLS(IBAS,NPSV)*PINSL(JBAS,NPSV)
C
              YMAT(IBAS+NBASA,JBAS      ) = YMAT(IBAS+NBASA,JBAS      )
     &      + SBFC*ASL*ALW*RKNT(NPSV)*PINSL(IBAS,NPSV)*PINLS(JBAS,NPSV)
C
              YMAT(IBAS+NBASA,JBAS+NBASA) = YMAT(IBAS+NBASA,JBAS+NBASA)
     &      + SBFC*ASS*ALW*RKNT(NPSV)*PINSL(IBAS,NPSV)*PINSL(JBAS,NPSV)
C
            ENDDO
          ENDDO
        ENDDO
C
C       SKIP POINT FOR ANGULAR SELECTION RULE
101     CONTINUE
C
C     END LOOP OVER KQNB
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE SLFHI0(SMAT,EXL,IZ,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           SSSSSS  LL       FFFFFFFF HH    HH IIII  000000            C
C          SS    SS LL       FF       HH    HH  II  00   000           C
C          SS       LL       FF       HH    HH  II  00  0000           C
C           SSSSSS  LL       FFFFFF   HHHHHHHH  II  00 00 00           C
C                SS LL       FF       HH    HH  II  0000  00           C
C          SS    SS LL       FF       HH    HH  II  000   00           C
C           SSSSSS  LLLLLLLL FF       HH    HH IIII  000000            C
C                                                                      C
C -------------------------------------------------------------------- C
C  SLFHI0 GENERATES ATOMIC ELECTRON SELF-INTERACTION (HIGH-ENERGY)     C
C  MATRIX ELEMENTS FOR SYMMETRY TYPE KQNA.                             C
C -------------------------------------------------------------------- C
C  ▶ THE `LL' COMPONENT OF THIS CAN BE USED IN THE 'NORL' OPTION BUT   C
C    THE INPUT DECK ASSUMES THAT 'DHFQ' MEANS RELATIVISTIC.            C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RTT(MBS*MBS,4),SMAT(MBD,MBD),EXL(MBS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/FCTS/RFACT(0:80),SFACT(0:80)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     INITIALISE ARRAY
      DO IBAS=1,2*NBAS
        DO JBAS=1,2*NBAS
          SMAT(IBAS,JBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     LOW- AND HIGH-ENERGY PARTITION CUTOFF
      CUTK = CV
C     CUTK = ZNUC(IZ)*ZNUC(IZ)/CV
C
C     AMPLITUDE FOR FREE-WAVE SELF-INTERACTION TERM (B&S 19.3)
      AG2 = (DLOG(EMSS*CV/CUTK)-TWLG+11.0D0/24.0D0)/(3.0D0*CV*PI)
      AG2 =-AG2/(EMSS*EMSS*CV*CV)
C
C     INTEGRATION GRID DETAILS (WHEN QUADRATURE IS INVOKED)
      NRHO = 6000
      ERHO = 4.0D0*RNUC(IZ)/DFLOAT(NRHO)
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      TL  = DFLOAT(LQN+KQN+1)
C
C     UPPER CUTOFF PARAMETER (TRIAL PHASE)
C     CUTINF = 1.0D-4/(RNUC(IZ)*RNUC(IZ))
      CUTINF = 2476.0D0
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RTT,EXL,NBAS,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
C
C         INITIALISE RAW MATRIX ELEMENTS
          ULL = 0.0D0
          USS = 0.0D0
C
C         ULTRAVIOLET CUTOFF (RESEARCH PHASE)
C         IF(EI.GT.CUTINF.OR.EJ.GT.CUTINF) GOTO 101
C         IF(EIJ.GT.CUTINF) GOTO 101
C
C         RAW MATRIX ELEMENTS FOR A POINT-NUCLEUS
          IF(NMDL(IZ).EQ.'POINT') THEN
C
            IF(LQN.EQ.0) THEN
              ULL = 1.0D0/(4.0D0*PI)
            ELSEIF(LQN.EQ.1) THEN
              USS = TL*TL/(4.0D0*PI)
            ENDIF
C
C         RAW MATRIX ELEMENTS FOR A GAUSSIAN (OR BEST-FIT) NUCLEUS
          ELSEIF(NMDL(IZ).EQ.'GAUSS'.OR.NMDL(IZ).EQ.'FERMI') THEN
C
            ULL = FNUC(IZ,0)*GSSINT0(LQN  ,EIJ,XNUC(IZ,0))
            IF(NNUC(IZ).GT.0) THEN
              ULL = ULL + ZCDINT0(IZ,LQN  ,EIJ)
            ENDIF

            VSA = 4.0D0*EPR*GSSINT0(LQN+1,EIJ,XNUC(IZ,0))
            IF(KQN.LT.0) THEN
              VSN = 0.0D0
            ELSE
              VSN = -2.0D0*EIJ*TL*GSSINT0(LQN  ,EIJ,XNUC(IZ,0))
     &                     +TL*TL*GSSINT0(LQN-1,EIJ,XNUC(IZ,0))
            ENDIF
            USS = USS + FNUC(IZ,0)*(VSA+VSN)
            IF(NNUC(IZ).GT.0) THEN
              VSA = 4.0D0*EPR*ZCDINT0(IZ,LQN+1,EIJ)
              IF(KQN.LT.0) THEN
                VSN = 0.0D0
              ELSE
                VSN = -2.0D0*EIJ*TL*ZCDINT0(IZ,LQN  ,EIJ)
     &                       +TL*TL*ZCDINT0(IZ,LQN-1,EIJ)
              ENDIF
              USS = USS + (VSA+VSN)
            ENDIF
C
C         RAW MATRIX ELEMENTS FOR A GENERIC NUCLEUS (BY QUADRATURE)
          ELSE
C
C           USE A LINEAR GRID TO INTEGRATE FROM 0 TO 4*RNUC
            DO N=0,NRHO
              RN  = ERHO*DFLOAT(N)
              ZLL = RN**(2*LQN+2)
              ZSS = RN**(2*LQN  )
              ZSS = ZSS*(TL-2.0D0*EI*RN*RN)*(TL-2.0D0*EJ*RN*RN)
              RSM = DEXP(-EIJ*RN*RN)
              PZR = RHONUC(NMDL(IZ),IZ,RN)
              ULL = ULL + ERHO*EXTINT11(ZLL*RSM*PZR,N,NRHO)
              USS = USS + ERHO*EXTINT11(ZSS*RSM*PZR,N,NRHO)
            ENDDO
C
C           INCLUDE NEWTON-COTES MULTIPLICATIVE FACTOR
            ULL = 5.0D0*ULL/299376.0D0
            USS = 5.0D0*USS/299376.0D0
C
          ENDIF
C
C         SKIP POINT FOR ULTRAVIOLET CUTOFF PARAMETER
101       CONTINUE
C
C         TRANSFER INTO ARRAY
          SMAT(IBAS     ,JBAS     ) =-AG2*4.0D0*PI*ZNUC(IZ)*RTT(M,1)*ULL
          SMAT(IBAS+NBAS,JBAS+NBAS) =-AG2*4.0D0*PI*ZNUC(IZ)*RTT(M,4)*USS
C
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE UEHNC0(UMAT,EXL,IZ,KQN,NBAS,COUPLING)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        UU    UU EEEEEEEE HH    HH NN    NN  CCCCCC   000000          C
C        UU    UU EE       HH    HH NNN   NN CC    CC 00   000         C
C        UU    UU EE       HH    HH NNNN  NN CC       00  0000         C
C        UU    UU EEEEEE   HHHHHHHH NN NN NN CC       00 00 00         C
C        UU    UU EE       HH    HH NN  NNNN CC       0000  00         C
C        UU    UU EE       HH    HH NN   NNN CC    CC 000   00         C
C         UUUUUU  EEEEEEEE HH    HH NN    NN  CCCCCC   000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  UEHNC0 GENERATES ATOMIC NUCLEAR VACUUM POLARISATION INTERACTION     C
C  MATRIX ELEMENTS FOR KQNA (UEHLING INTERACTION TO FIRST ORDER).      C
C -------------------------------------------------------------------- C
C  ▶ THE `LL' COMPONENT OF THIS CAN BE USED IN THE 'NORL' OPTION BUT   C
C    THE INPUT DECK ASSUMES THAT 'DHFQ' MEANS RELATIVISTIC.            C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*4 COUPLING
      CHARACTER*5 NMDL
C
      DIMENSION RN(MBS*MBS,4),UMAT(MBD,MBD),EXL(MBS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      RL  = DFLOAT(LQN)
      TL  = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
C
C         MATRIX ELEMENT DEPENDS ON HAMILTONIAN
          IF(HMLT.EQ.'NORL') THEN
C
C           NON-RELATIVISTIC CASE
C           (HAVE TO SELECT MANUALLY -- CONFLICT IN CARDIN INPUT DECK.)
            ULL = UEHINT0(IZ,LQN  ,EIJ,COUPLING)
            UMAT(IBAS,JBAS) = ZNUC(IZ)*RN(M,1)*ULL
C
          ELSE
C
C           RELATIVISTIC CASE
            ULL = UEHINT0(IZ,LQN  ,EIJ,COUPLING)
C
            VSA = 4.0D0*EPR*UEHINT0(IZ,LQN+1,EIJ,COUPLING)
            IF(KQN.LT.0) THEN
              VSN = 0.0D0
            ELSE
              VSN =-2.0D0*EIJ*TL*UEHINT0(IZ,LQN  ,EIJ,COUPLING)
     &                    +TL*TL*UEHINT0(IZ,LQN-1,EIJ,COUPLING)
            ENDIF
            USS = VSA+VSN
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NBAS
            LBAS = JBAS+NBAS
C
C           TRANSFER INTO ARRAY
            UMAT(IBAS,JBAS) = ZNUC(IZ)*RN(M,1)*ULL
            UMAT(IBAS,LBAS) = 0.0D0
            UMAT(KBAS,JBAS) = 0.0D0
            UMAT(KBAS,LBAS) = ZNUC(IZ)*RN(M,4)*USS
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
      SUBROUTINE WKRNC0(UMAT,EXL,IZ,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      WW        WW KK    KK RRRRRRR  NN    NN  CCCCCC   000000        C
C      WW        WW KK   KK  RR    RR NNN   NN CC    CC 00   000       C
C      WW   WW   WW KK  KK   RR    RR NNNN  NN CC       00  0000       C
C      WW  WWWW  WW KKKKK    RR    RR NN NN NN CC       00 00 00       C
C      WW WW  WW WW KK  KK   RRRRRRR  NN  NNNN CC       0000  00       C
C      WWWW    WWWW KK   KK  RR    RR NN   NNN CC    CC 000   00       C
C      WW        WW KK    KK RR    RR NN    NN  CCCCCC   000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  WKRNC0 GENERATES ATOMIC NUCLEAR VACUUM POLARISATION INTERACTION     C
C  MATRIX ELEMENTS FOR KQNA (WICHMANN-KROLL INTERACTION).              C
C -------------------------------------------------------------------- C
C  ▶ THE `LL' COMPONENT OF THIS CAN BE USED IN THE 'NORL' OPTION BUT   C
C    THE INPUT DECK ASSUMES THAT 'DHFQ' MEANS RELATIVISTIC.            C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RN(MBS*MBS,4),UMAT(MBD,MBD),EXL(MBS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      RL  = DFLOAT(LQN)
      TL  = DFLOAT(2*LQN+1)
C
C     NUCLEAR CHARGE FACTOR
      Z3 = ZNUC(IZ)*ZNUC(IZ)*ZNUC(IZ)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
C
C         MATRIX ELEMENT DEPENDS ON HAMILTONIAN
          IF(HMLT.EQ.'NORL') THEN
C
C           NON-RELATIVISTIC CASE
C           (HAVE TO SELECT MANUALLY -- CONFLICT IN CARDIN INPUT DECK.)
            ULL = WKRINT0(IZ,LQN  ,EIJ)
            UMAT(IBAS,JBAS) = ZNUC(IZ)*RN(M,1)*ULL
C
          ELSE
C
C           RELATIVISTIC CASE
            ULL = WKRINT0(IZ,LQN  ,EIJ)
C
            VSA = 4.0D0*EPR*WKRINT0(IZ,LQN+1,EIJ)
            IF(KQN.LT.0) THEN
              VSN = 0.0D0
            ELSE
              VSN =-2.0D0*EIJ*TL*WKRINT0(IZ,LQN  ,EIJ)
     &                    +TL*TL*WKRINT0(IZ,LQN-1,EIJ)
            ENDIF
            USS = VSA+VSN
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NBAS
            LBAS = JBAS+NBAS
C
C           TRANSFER INTO ARRAY
            UMAT(IBAS,JBAS) = Z3*RN(M,1)*ULL
            UMAT(IBAS,LBAS) = 0.0D0
            UMAT(KBAS,JBAS) = 0.0D0
            UMAT(KBAS,LBAS) = Z3*RN(M,4)*USS
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
      SUBROUTINE KSBNC0(UMAT,EXL,IZ,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         KK    KK SSSSSS  BBBBBBB  NN    NN  CCCCCC   000000          C
C         KK   KK SS    SS BB    BB NNN   NN CC    CC 00   000         C
C         KK  KK  SS       BB    BB NNNN  NN CC       00  0000         C
C         KKKKK    SSSSSS  BBBBBBB  NN NN NN CC       00 00 00         C
C         KK  KK        SS BB    BB NN  NNNN CC       0000  00         C
C         KK   KK SS    SS BB    BB NN   NNN CC    CC 000   00         C
C         KK    KK SSSSSS  BBBBBBB  NN    NN  CCCCCC   000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  KSBNC0 GENERATES ATOMIC NUCLEAR KÄLLÉN-SABRY INTERACTION MATRIX     C
C  ELEMENTS FOR KQNA SYMMETRY TYPE.                                    C
C -------------------------------------------------------------------- C
C  ▶ THE `LL' COMPONENT OF THIS CAN BE USED IN THE 'NORL' OPTION BUT   C
C    THE INPUT DECK ASSUMES THAT 'DHFQ' MEANS RELATIVISTIC.            C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RN(MBS*MBS,4),UMAT(MBD,MBD),EXL(MBS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      RL  = DFLOAT(LQN)
      TL  = DFLOAT(2*LQN+1)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
C
C         MATRIX ELEMENT DEPENDS ON HAMILTONIAN
          IF(HMLT.EQ.'NORL') THEN
C
C           NON-RELATIVISTIC CASE
C           (HAVE TO SELECT MANUALLY -- CONFLICT IN CARDIN INPUT DECK.)
            ULL = RKSBINT0(IZ,LQN  ,EIJ)
            UMAT(IBAS,JBAS) = ZNUC(IZ)*RN(M,1)*ULL
C
          ELSE
C
C           RELATIVISTIC CASE
            ULL = RKSBINT0(IZ,LQN  ,EIJ)
C
            VSA = 4.0D0*EPR*RKSBINT0(IZ,LQN+1,EIJ)
            IF(KQN.LT.0) THEN
              VSN = 0.0D0
            ELSE
              VSN =-2.0D0*EIJ*TL*RKSBINT0(IZ,LQN  ,EIJ)
     &                    +TL*TL*RKSBINT0(IZ,LQN-1,EIJ)
            ENDIF
            USS = VSA+VSN
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NBAS
            LBAS = JBAS+NBAS
C
C           TRANSFER INTO ARRAY
            UMAT(IBAS,JBAS) = ZNUC(IZ)*RN(M,1)*ULL
            UMAT(IBAS,LBAS) = 0.0D0
            UMAT(KBAS,JBAS) = 0.0D0
            UMAT(KBAS,LBAS) = ZNUC(IZ)*RN(M,4)*USS
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
      SUBROUTINE UEHEL0(UMAT,EXL,IZ,KQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         UU    UU EEEEEEEE HH    HH EEEEEEEE LL       000000          C
C         UU    UU EE       HH    HH EE       LL      00   000         C
C         UU    UU EE       HH    HH EE       LL      00  0000         C
C         UU    UU EEEEEE   HHHHHHHH EEEEEE   LL      00 00 00         C
C         UU    UU EE       HH    HH EE       LL      0000  00         C
C         UU    UU EE       HH    HH EE       LL      000   00         C
C          UUUUUU  EEEEEEEE HH    HH EEEEEEEE LLLLLLLL 000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  UEHEL0 GENERATES ATOMIC ELECTRONIC VACUUM POLARISATION INTERACTION  C
C  MATRIX ELEMENTS FOR KQNA (UEHLING INTERACTION TO FIRST ORDER).      C
C -------------------------------------------------------------------- C
C  I GET PRETTY SIMPLE WHEN EVALUATING THE DENSITY HERE. ASSUME IT'S   C
C  JUST THE LARGE-COMPONENT S-LIKE ORBITALS, AND THEY'RE ALL FULL.     C
C -------------------------------------------------------------------- C
C  ▶ THE `LL' COMPONENT OF THIS CAN BE USED IN THE 'NORL' OPTION BUT   C
C    THE INPUT DECK ASSUMES THAT 'DHFQ' MEANS RELATIVISTIC.            C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*5 NMDL
C
      DIMENSION RN(MBS*MBS,4),UMAT(MBD,MBD),EXL(MBS)
      DIMENSION RS(MBS*MBS,4),DS(MBS*MBS),EXS(MBS)
      DIMENSION CL(MBS,MKP+1)
C
      COMPLEX*16 COEF(MDM,MDM)
C
      COMMON/AUFB/NORB(0:MEL,MKP+1),NUMOCC(0:MEL),LMXCONF
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/OCPD/IOCPN(MDM),IOCCM0,IOCTP(MCT,MKP,MKP+1,MKP+1)
C
C     CALL IN THE S-TYPE FUNCTIONS ON THIS CENTRE
      NS = NFNC(0,IZ)
      DO K=1,NS
        EXS(K) = BEXL(K,0,IZ)
      ENDDO
C
C     NORMALISATION
      CALL RNORM0(RS,EXS,NS,0)
C
C     S-TYPE ORBITAL COEFFICIENTS
      DO NQN=1,NUMOCC(0)
C
C       SKIP ANY UNOCCUPIED LEVELS
        IF(NORB(0,NQN).NE.0) THEN
C
C         COEFFICIENT MATRIX STARTING ADDRESSES
          IL1  = LRGE(IZ,1,1)
          IROW = IOCTP(IZ,1,1,NQN)
C
C         LOOP OVER BASIS FUNCTIONS
          DO IBAS=1,NS
            CL(IBAS,NQN) = DREAL(COEF(IL1+IBAS,IROW+NSKP))
          ENDDO
C        
        ENDIF
C      
      ENDDO
C
C     BUILD ATOMIC CHARGE DENSITY MATRIX FOR S-STATES
      N = 0
      TOT = 0.0D0
      DO KBAS=1,NS
        EK  = EXS(KBAS)
        DO LBAS=1,NS
          N = N+1
          EL  = EXS(LBAS)
          EKL = EK+EL
          DS(N) = 0.0D0
          DO NQN=1,NUMOCC(0)
            DS(N) = DS(N) + 2.0D0*CL(KBAS,NQN)*CL(LBAS,NQN)
          ENDDO
          TOT = TOT + DS(N)*RS(N,1)*PI12*0.25D0/(EKL*DSQRT(EKL))
        ENDDO
      ENDDO
C
C     DETERMINE THE LQN
      LQN = LVAL(KQN)
      RL  = DFLOAT(LQN)
      TL  = DFLOAT(2*LQN+1)
      AMP =-CMPW*CMPW/(240.0D0*PI*PI*PI*CV)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NBAS,LQN)
C
C     LOOP OVER PAIRS OF BASIS FUNCTIONS WITHIN THIS KQN BLOCK
      M = 0
      DO IBAS=1,NBAS
        EI = EXL(IBAS)
        DO JBAS=1,NBAS
          M   = M+1
          EJ  = EXL(JBAS)
          EIJ = EI+EJ
          EPR = EI*EJ
C
C         MATRIX ELEMENT DEPENDS ON HAMILTONIAN
          IF(HMLT.EQ.'NORL') THEN
C
C           NON-RELATIVISTIC CASE
            ULL = DNSINT0(LQN  ,EIJ,DS,RS,EXS,NS)
C
C           TRANSFER INTO DIRAC MATRIX
            UMAT(IBAS,JBAS) = AMP*RN(M,1)*ULL
C
          ELSE
C
C           RELATIVISTIC CASE
            ULL = DNSINT0(LQN  ,EIJ,DS,RS,EXS,NS)
C
            VSA = 4.0D0*EPR*DNSINT0(LQN+1,EIJ,DS,RS,EXS,NS)
            IF(KQN.LT.0) THEN
              VSN = 0.0D0
            ELSE
              VSN = -2.0D0*EIJ*TL*DNSINT0(LQN  ,EIJ,DS,RS,EXS,NS)
     &                     +TL*TL*DNSINT0(LQN-1,EIJ,DS,RS,EXS,NS)
            ENDIF
            USS = VSA+VSN
C
C           TRANSFER INTO NUCLEAR OVERLAP MATRIX
            UMAT(IBAS     ,JBAS     ) = AMP*RN(M,1)*ULL
            UMAT(IBAS     ,JBAS+NBAS) = 0.0D0
            UMAT(IBAS+NBAS,JBAS     ) = 0.0D0
            UMAT(IBAS+NBAS,JBAS+NBAS) = AMP*RN(M,4)*USS
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
      SUBROUTINE COULM0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        CCCCCC   OOOOOO  UU    UU LL      MM       MM  000000         C
C       CC    CC OO    OO UU    UU LL      MMM     MMM 00   000        C
C       CC       OO    OO UU    UU LL      MMMM   MMMM 00  0000        C
C       CC       OO    OO UU    UU LL      MM MM MM MM 00 00 00        C
C       CC       OO    OO UU    UU LL      MM  MMM  MM 0000  00        C
C       CC    CC OO    OO UU    UU LL      MM   M   MM 000   00        C
C        CCCCCC   OOOOOO   UUUUUU  LLLLLLL MM       MM  000000         C
C                                                                      C
C -------------------------------------------------------------------- C
C  COULM0 CONSTRUCTS THE ATOMIC COULOMB MATRIX FROM RADIAL DIRECT      C
C  AND EXCHANGE INTEGRALS AND A MEAN-FIELD CHARGE DENSITY.             C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      DIMENSION RJLLLL(MB2,4),RJSSSS(MB2,4),RJLLSS(MB2,4),RJSSLL(MB2,4),
     &          RKLLLL(MB2,4),RKSSSS(MB2,4),RKSLSL(MB2,4)
C
      COMMON/ATMC/G11(MBD,MBD),G21(MBD,MBD),G12(MBD,MBD),G22(MBD,MBD)
      COMMON/ATMD/DLL1(MB2),DSL1(MB2),DSS1(MB2),DLS1(MB2),
     &            DLL2(MB2),DSL2(MB2),DSS2(MB2),DLS2(MB2)
      COMMON/B0IJ/EIJ(-MTN:MTN),RNIJ(4),EI,EJ,IBAS,JBAS
      COMMON/B0QN/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/RCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
C
C     INITIALISE COULOMB MATRIX
      DO IBAS=1,MBD
        DO JBAS=1,MBD
          G11(IBAS,JBAS) = 0.0D0
          G21(IBAS,JBAS) = 0.0D0
          G12(IBAS,JBAS) = 0.0D0
          G22(IBAS,JBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     COEFFICIENTS FOR INTEGRAL ASSEMBLY
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
      TI = DFLOAT(2*LQNA+1)
      TJ = DFLOAT(2*LQNA+1)
      TK = DFLOAT(2*LQNB+1)
      TL = DFLOAT(2*LQNB+1)
C
      T0000 = 1.0D0
      T1000 = TI
      T0100 = TJ
      T0010 = TK
      T0001 = TL
      T1100 = TI*TJ
      T1010 = TI*TK
      T1001 = TI*TL
      T0110 = TJ*TK
      T0101 = TJ*TL
      T0011 = TK*TL
      T1110 = TI*TJ*TK
      T1101 = TI*TJ*TL
      T1011 = TI*TK*TL
      T0111 = TJ*TK*TL
      T1111 = TI*TJ*TK*TL
C
C     EVALUATE CLOSED-SHELL ELECTRON REPULSION ANGULAR INTEGRALS
      CALL ANGCLM0
C
C     BASIS SET INTERMEDIATES FOR THE KL-PAIRS
      CALL KLSET0
C
C     ITERATE OVER ALL MATRIX ELEMENTS
      DO IBAS=1,NBASA
        DO JBAS=1,NBASA
C
C         GAUSSIAN EXPONENTS FOR THIS PAIR
          EI = EXLA(IBAS)
          EJ = EXLA(JBAS)
C
C         BASIS SET INTERMEDIATES FOR THE IJ-PAIRS
          CALL IJSET0
C
C         GENERATE BATCH OF RADIAL INTEGRALS (J AND K MATRICES)
          CALL RKCLM0(RJLLLL,RJSSSS,RJLLSS,RJSSLL,RKLLLL,RKSSSS,RKSLSL)
C
          IF(HMLT.EQ.'NORL') THEN
C         NON-RELATIVISTIC HAMILTONIAN
C
C           INITIALISE COUNTER
            GLL = 0.0D0
C
C           BUILD THE FOCK MATRIX
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,1)*DLL2(M) - RKLLLL(M,1)*DLL2(M)
            ENDDO
C
C           TRANSFER COUNTER VALUE TO COULOMB MATRIX
            G22(IBAS,JBAS) = GLL
C
          ELSE
C         RELATIVISTIC HAMILTONIAN
C
C           SMALL-COMPONENT MATRIX ADDRESSES
            KBAS = IBAS + NBASA
            LBAS = JBAS + NBASA
C
C    (22)   KQNA < 0 AND KQNB < 0  CONTRIBUTIONS (CANNOT SKIP)
C
C           INITIALISE COUNTERS
            GLL = 0.0D0
            GSL = 0.0D0
            GSS = 0.0D0
C
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,1)*DLL2(M) - RKLLLL(M,1)*DLL2(M)
     &                  + RJLLSS(M,1)*DSS2(M)
              GSL = GSL                       - RKSLSL(M,1)*DSL2(M)
              GSS = GSS + RJSSSS(M,1)*DSS2(M) - RKSSSS(M,1)*DSS2(M)
     &                  + RJSSLL(M,1)*DLL2(M)
            ENDDO
C
C           TRANSFER COUNTER VALUES TO COULOMB MATRIX
            G22(IBAS,JBAS) = GLL
            G22(KBAS,JBAS) = GSL
            G22(JBAS,KBAS) = GSL
            G22(KBAS,LBAS) = GSS
C
C    (21)   KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNB.EQ.0) GOTO 200
C
C           INITIALISE COUNTERS
            GLL = 0.0D0
            GSL = 0.0D0
            GSS = 0.0D0
C
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,2)*DLL1(M) - RKLLLL(M,2)*DLL1(M)
     &                  + RJLLSS(M,2)*DSS1(M)
              GSL = GSL                       - RKSLSL(M,2)*DSL1(M)
              GSS = GSS + RJSSSS(M,2)*DSS1(M) - RKSSSS(M,2)*DSS1(M)
     &                  + RJSSLL(M,2)*DLL1(M)
            ENDDO
C
C           TRANSFER COUNTER VALUES TO COULOMB MATRIX
            G21(IBAS,JBAS) = GLL
            G21(KBAS,JBAS) = GSL
            G21(JBAS,KBAS) = GSL
            G21(KBAS,LBAS) = GSS
C
200         CONTINUE
C
C    (12)   KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNA.EQ.0) GOTO 300
C
C           INITIALISE COUNTERS
            GLL = 0.0D0
            GSL = 0.0D0
            GSS = 0.0D0
C
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,3)*DLL2(M) - RKLLLL(M,3)*DLL2(M)
     &                  + RJLLSS(M,3)*DSS2(M)
              GSL = GSL                       - RKSLSL(M,3)*DSL2(M)
              GSS = GSS + RJSSSS(M,3)*DSS2(M) - RKSSSS(M,3)*DSS2(M)
     &                  + RJSSLL(M,3)*DLL2(M)
            ENDDO
C
C           TRANSFER COUNTER VALUES TO COULOMB MATRIX
            G12(IBAS,JBAS) = GLL
            G12(KBAS,JBAS) = GSL
            G12(JBAS,KBAS) = GSL
            G12(KBAS,LBAS) = GSS
C
300         CONTINUE
C
C    (11)   KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 400
C
C           INITIALISE COUNTERS
            GLL = 0.0D0
            GSL = 0.0D0
            GSS = 0.0D0
C
C           SUM OVER MEAN FIELD CONTRIBUTIONS FOR THIS BASIS PAIR
            DO M=1,MAXM
              GLL = GLL + RJLLLL(M,4)*DLL1(M) - RKLLLL(M,4)*DLL1(M)
     &                  + RJLLSS(M,4)*DSS1(M)
              GSL = GSL                       - RKSLSL(M,4)*DSL1(M)
              GSS = GSS + RJSSSS(M,4)*DSS1(M) - RKSSSS(M,4)*DSS1(M)
     &                  + RJSSLL(M,4)*DLL1(M)
            ENDDO
C
C           TRANSFER COUNTER VALUES TO COULOMB MATRIX
            G11(IBAS,JBAS) = GLL
            G11(KBAS,JBAS) = GSL
            G11(JBAS,KBAS) = GSL
            G11(KBAS,LBAS) = GSS
C
400         CONTINUE
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
      SUBROUTINE RKCLM0(RJLLLL,RJSSSS,RJLLSS,RJSSLL,
     &                                       RKLLLL,RKSSSS,RKSLSL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       RRRRRRR  KK    KK CCCCCC  LL       MM       MM  000000         C
C       RR    RR KK   KK CC    CC LL       MMM     MMM 00   000        C
C       RR    RR KK  KK  CC       LL       MMMM   MMMM 00  0000        C
C       RR    RR KKKKK   CC       LL       MM MM MM MM 00 00 00        C
C       RRRRRRR  KK  KK  CC       LL       MM  MMM  MM 0000  00        C
C       RR    RR KK   KK CC    CC LL       MM   M   MM 000   00        C
C       RR    RR KK    KK CCCCCC  LLLLLLLL MM       MM  000000         C
C                                                                      C
C -------------------------------------------------------------------- C
C  RKCLM0 EVALUATES A DIRECT AND EXCHANGE BATCH OF ELECTRON REPULSION  C
C  INTEGRALS OF ALL COMPONENT LABEL COMBINATIONS L AND S IN THE ATOMIC C
C  SCF PROCEDURE FOR A USER-SPECIFIED HAMILTONIAN.                     C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C  ▶ RJLLLL(M,N) - DIRECT ERI OVERLAP {LL,LL}                          C
C  ▶ RJSSSS(M,N) - DIRECT ERI OVERLAP {SS,SS}                          C
C  ▶ RJLLSS(M,N) - DIRECT ERI OVERLAP {LL,SS}                          C
C  ▶ RJSSLL(M,N) - DIRECT ERI OVERLAP {SS,LL}                          C
C  ▶ RKLLLL(M,N) - EXCHANGE ERI OVERLAP {LL,LL}                        C
C  ▶ RKSSSS(M,N) - EXCHANGE ERI OVERLAP {SS,SS}                        C
C  ▶ RKSLSL(M,N) - EXCHANGE ERI OVERLAP {SL,SL}                        C
C -------------------------------------------------------------------- C
C  COULOMB MATRIX LIST INDEX:                                          C
C  ▶ N=1 - KQN(A)<0, KQN(B)<0 (TYPICAL LABEL 22)                       C
C  ▶ N=2 - KQN(A)<0, KQN(B)>0 (TYPICAL LABEL 12)                       C
C  ▶ N=3 - KQN(A)>0, KQN(B)<0 (TYPICAL LABEL 21)                       C
C  ▶ N=4 - KQN(A)>0, KQN(B)>0 (TYPICAL LABEL 11)                       C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      DIMENSION XJ(MB2,2),XK(MB2,2),IAA(2),IBB(2)
      DIMENSION BETA(MB2),BTA1(MB2),XROOT(MB2),BIN(MB2),TRM(MB2)
      DIMENSION BDU(MB2,-MTN:MTN,-MTN:MTN),BDL(MB2,-MTN:MTN,-MTN:MTN)
      DIMENSION BXU(MB2,-MTN:MTN,-MTN:MTN),BXL(MB2,-MTN:MTN,-MTN:MTN)
C
      DIMENSION RTIK0(MBS),RTJL0(MBS),PTIK0(MBS),PTJL0(MBS)
      DIMENSION RJLLLL(MB2,4),RJSSSS(MB2,4),RJLLSS(MB2,4),RJSSLL(MB2,4),
     &          RKLLLL(MB2,4),RKSSSS(MB2,4),RKSLSL(MB2,4)
C
      COMMON/B0IJ/EIJ(-MTN:MTN),RNIJ(4),EI,EJ,IBAS,JBAS
      COMMON/B0IK/EIK(MB2,-MTN:MTN),IKIND(MB2)
      COMMON/B0JL/EJL(MB2,-MTN:MTN),JLIND(MB2)
      COMMON/B0KL/EKL(MB2,-MTN:MTN),RNKL(MB2,4),EK(MB2),EL(MB2)
      COMMON/B0QN/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/RCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
      COMMON/XTNS/BK(MNU,4),ELL(MNU,4),ESS(MNU,4),ESL(MNU,4),GSL(MNU,4)
C
C**********************************************************************C
C     PREPARATION OF EXPONENT POWERS AND BETA FUNCTION ARGUMENTS.      C
C**********************************************************************C
C
C     TENSOR ORDER LIMITS
      NUI = NUS(1)
      NUF = NUS(NUNUM)
C
C     GAUSSIAN EXPONENT FOR (IBAS,JBAS) OVERLAP
      TIJ0 = EI+EJ
C
C     INITIALISE THE ARRAY XJ(M,2) FOR INCOMPLETE BETA FUNCTION ARGS.
      DO M=1,MAXM
        TKL0    = EK(M)+EL(M)
        TIJKL   = TIJ0+TKL0
        XJ(M,1) = TIJ0/TIJKL
        XJ(M,2) = TKL0/TIJKL
      ENDDO
C
C     LOWEST EXPONENT POWER
      IPOWER = LQNA+LQNB-NUF
C
C     A BLOCK OF BASIS EXPONENT PRODUCTS
      DO KBAS=1,NBASB
        RTIK0(KBAS) = DSQRT(EI+EXLB(KBAS))
        RTJL0(KBAS) = DSQRT(EJ+EXLB(KBAS))
        PTIK0(KBAS) = RTIK0(KBAS)**(-IPOWER)
        PTJL0(KBAS) = RTJL0(KBAS)**(-IPOWER)
      ENDDO
C
C     CALCULATE A FULL SET OF EXPONENT OVERLAPS FOR EXCHANGE
      DO M=1,MAXM
        RTIK  = RTIK0(IKIND(M))
        RTJL  = RTJL0(JLIND(M))
        EIK(M,-NUF) = PTIK0(IKIND(M))
        EJL(M,-NUF) = PTJL0(JLIND(M))
        DO IPOW=-NUF+1,NUF+5
          EIK(M,IPOW) = EIK(M,IPOW-1)/RTIK
          EJL(M,IPOW) = EJL(M,IPOW-1)/RTJL
        ENDDO
      ENDDO
C
C     CALCULATE LIST OF BETA FUNCTION ARGUMENTS, XK(MB2,2) (Z AND Z')
      DO M=1,MAXM
        TKL0    = EK(M)+EL(M)
        TIJKL   = EI+EJ+TKL0
        XK(M,1) = (EI+EK(M))/TIJKL
        XK(M,2) = (EJ+EL(M))/TIJKL
      ENDDO
C
C**********************************************************************C
C     INLINE BETA FUNCTION CODE FOR DIRECT TERMS                       C
C**********************************************************************C
C
C     NUMBER OF LEVELS NEEDED TO ACCOUNT FOR ALL DIRECT INTEGRALS
      NVALS = 3
C
C     LOOP OVER ORDER FOR FIRST INDEX
      DO NX=1,NVALS
C
C       TWICE THE ACTUAL FAMILY VALUE
        IAA(1) = 2*LQNA+2*NX-1
        IAA(2) = 2*LQNB+2*NX-1
C
C       LOOP OVER ORDER FOR SECOND INDEX
        DO NY=1,NVALS
          IBB(1) = 2*LQNB+2*NY-2
          IBB(2) = 2*LQNA+2*NY-2
C
C         LOOP OVER BETA INTEGRAL TYPE
          DO IBETA=1,2
            IA =(IAA(IBETA)-1)/2
            IB = IBB(IBETA)   /2
C
C           BEGIN CONDITIONAL STATEMENT OVER IB VALUES
C           CASE 1: IB > 1
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
              RA = RA+1.0D0
              RB = RB+1.0D0
              RC = RC+1.0D0
              RD = RD+1.0D0
              DO J=2,IB-1
                RCD = RC*RD
                FCT = RA*RB/RCD
                DO M=1,MAXM
                  TRM(M) = FCT*TRM(M)*XJ(M,IBETA)
                  BIN(M) = BIN(M)+TRM(M)
                ENDDO
                RA = RA+1.0D0
                RB = RB+1.0D0
                RC = RC+1.0D0
                RD = RD+1.0D0
              ENDDO
              DO M=1,MAXM
                BETA(M) = BTA1(M)*BIN(M)
              ENDDO
C
C           CASE 2: IB = 1
            ELSEIF(IB.EQ.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BETA(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
C
C           CASE 3: IB = 0
            ELSEIF(IB.EQ.0) THEN
              DO M=1,MAXM
                XROOT(M) = DSQRT(XJ(M,IBETA))
                DEN      =  1.0D0-XROOT(M)
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
C
            ENDIF
C
C           SORT THE BETA INTEGRAL INTO THE CORRECT ARRAY
            IF(IBETA.EQ.1) THEN
              DO M=1,MAXM
                BDU(M,2*NX-1,2*NY-2) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXM
                BDL(M,2*NX-1,2*NY-2) = BETA(M)
              ENDDO
            ENDIF
C
C         END LOOPS OVER IBETA AND INDICES NX, NY
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INLINE BETA FUNCTION CODE FOR EXCHANGE TERMS                     C
C**********************************************************************C
C
C     NUMBER OF LEVELS NEEDED TO ACCOUNT FOR ALL EXCHANGE INTEGRALS
      NVALS = (NUF-NUI)/2+3
C
C     LOOP OVER ORDER FOR FIRST INDEX
      DO NX=1,NVALS
        IAA(1) = LQNA+LQNB+NUI+2*NX-1
        IAA(2) = LQNA+LQNB+NUI+2*NX-1
C
C       LOOP OVER ORDER FOR SECOND INDEX
        DO NY=1,NVALS
          IBB(1) = LQNA+LQNB-NUF+2*NY-2
          IBB(2) = LQNA+LQNB-NUF+2*NY-2
C
C         LOOP OVER BETA INTEGRAL TYPE
          DO IBETA=1,2
            IA = (IAA(IBETA)-1)/2
            IB =  IBB(IBETA)   /2
C
C           BEGIN CONDITIONAL STATEMENT OVER IB VALUES
C           CASE 1: IB > 1
            IF(IB.GT.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BTA1(M) = (DSQRT(XK(M,IBETA))**IX)/X
              ENDDO
              RA = X
              RB = DFLOAT(1-IB)
              RC = 1.0D0+X
              RD = 1.0D0
              RCD = RC*RD
              FCT = RA*RB/RCD
              DO M=1,MAXM
                TRM(M) = FCT*XK(M,IBETA)
                BIN(M) = 1.0D0+TRM(M)
              ENDDO
              RA = RA+1.0D0
              RB = RB+1.0D0
              RC = RC+1.0D0
              RD = RD+1.0D0
              DO J=2,IB-1
                RCD = RC*RD
                FCT = RA*RB/RCD
                DO M=1,MAXM
                  TRM(M) = FCT*TRM(M)*XK(M,IBETA)
                  BIN(M) = BIN(M)+TRM(M)
                ENDDO
                RA = RA+1.0D0
                RB = RB+1.0D0
                RC = RC+1.0D0
                RD = RD+1.0D0
              ENDDO
              DO M=1,MAXM
                BETA(M) = BTA1(M)*BIN(M)
              ENDDO
C
C           CASE 2: IB = 1
            ELSEIF(IB.EQ.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BETA(M) = (DSQRT(XK(M,IBETA))**IX)/X
              ENDDO
C
C           CASE 3: IB = 0
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
C
            ENDIF
C
C           SORT THE BETA INTEGRAL INTO THE UPPER/LOWER ARRAYS
            IF(IBETA.EQ.1) THEN
              DO M=1,MAXM
                BXU(M, NUI+2*NX-1,-NUF+2*NY-2) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXM
                BXL(M, NUI+2*NX-1,-NUF+2*NY-2) = BETA(M)
              ENDDO
            ENDIF
C
C         END LOOPS OVER IBETA AND INDICES NX, NY
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     RADIAL INTEGRALS OVER SPINORS (SEPARATED BY TENSOR ORDER)        C
C**********************************************************************C
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
C     VALUE PREPARATION
      E0000 = 1.0D0
      E1000 = EI
      E0100 = EJ
      E1100 = EI*EJ
C
C     LOOP OVER K,L BASIS FUNCTIONS
      DO M=1,MAXM
C
C       MORE VALUE PREPARATION
        E0010 = EK(M)
        E0001 = EL(M)
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
        IF(HMLT.EQ.'NORL') THEN
C       NON-RELATIVISTIC HAMILTONIAN
C
          B32 = EIJ(2)*EKL(M,3)*BDL(M,3,2) + EIJ(3)*EKL(M,2)*BDU(M,3,2)
C
          RJLLLL(M,1) = V1*T0000*E0000*C5*B32
C
        ELSE
C       RELATIVISTIC HAMILTONIAN
C
C         TEMPORARY STORAGE OF RAW RJ(1,M) (LTEN=1 BECAUSE NU=0 ONLY)
          B10 = EIJ(0)*EKL(M,1)*BDL(M,1,0) + EIJ(1)*EKL(M,0)*BDU(M,1,0)
          B12 = EIJ(2)*EKL(M,1)*BDL(M,1,2) + EIJ(3)*EKL(M,0)*BDU(M,3,0)
          B14 = EIJ(4)*EKL(M,1)*BDL(M,1,4) + EIJ(5)*EKL(M,0)*BDU(M,5,0)
          B30 = EIJ(0)*EKL(M,3)*BDL(M,3,0) + EIJ(1)*EKL(M,2)*BDU(M,1,2)
          B32 = EIJ(2)*EKL(M,3)*BDL(M,3,2) + EIJ(3)*EKL(M,2)*BDU(M,3,2)
          B34 = EIJ(4)*EKL(M,3)*BDL(M,3,4) + EIJ(5)*EKL(M,2)*BDU(M,5,2)
          B50 = EIJ(0)*EKL(M,5)*BDL(M,5,0) + EIJ(1)*EKL(M,4)*BDU(M,1,4)
          B52 = EIJ(2)*EKL(M,5)*BDL(M,5,2) + EIJ(3)*EKL(M,4)*BDU(M,3,4)
          B54 = EIJ(4)*EKL(M,5)*BDL(M,5,4) + EIJ(5)*EKL(M,4)*BDU(M,5,4)
C
C         FILL RJ ARRAYS FOR THIS LQNA,LQNB BLOCK
C
C         KQNA < 0 AND KQNB < 0 CONTRIBUTIONS (CANNOT SKIP)
          RJLLLL(M,1) = V1*T0000*E0000*C5*B32
          RJLLSS(M,1) = V4*T0000*E0011*C7*B52
          RJSSLL(M,1) = V4*T0000*E1100*C7*B34
          RJSSSS(M,1) = VS*T0000*E1111*C9*B54
C
C         KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNB.EQ.0) GOTO 102
          RJLLLL(M,2) = V1*T0000*E0000*C5*B32
          RJLLSS(M,2) = V4*T0000*E0011*C7*B52
     &                - V2*T0001*E0010*C5*B32 - V2*T0010*E0001*C5*B32
     &                + V1*T0011*E0000*C3*B12
          RJSSLL(M,2) = V4*T0000*E1100*C7*B34
          RJSSSS(M,2) = VS*T0000*E1111*C9*B54
     &                - V8*T0001*E1110*C7*B34 - V8*T0010*E1101*C7*B34
     &                + V4*T0011*E1100*C5*B14
102       CONTINUE
C
C         KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0) GOTO 103
          RJLLLL(M,3) = V1*T0000*E0000*C5*B32
          RJLLSS(M,3) = V4*T0000*E0011*C7*B52
          RJSSLL(M,3) = V4*T0000*E1100*C7*B34
     &                - V2*T0100*E1000*C5*B32 - V2*T1000*E0100*C5*B32
     &                + V1*T1100*E0000*C3*B30
          RJSSSS(M,3) = VS*T0000*E1111*C9*B54
     &                - V8*T0100*E1011*C7*B52 - V8*T1000*E0111*C7*B52
     &                + V4*T1100*E0011*C5*B50
103       CONTINUE
C
C         KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 104
          RJLLLL(M,4) = V1*T0000*E0000*C5*B32
          RJLLSS(M,4) = V4*T0000*E0011*C7*B52
     &                - V2*T0001*E0010*C5*B32 - V2*T0010*E0001*C5*B32
     &                + V1*T0011*E0000*C3*B12
          RJSSLL(M,4) = V4*T0000*E1100*C7*B34
     &                - V2*T0100*E1000*C5*B32 - V2*T1000*E0100*C5*B32
     &                + V1*T1100*E0000*C3*B30
          RJSSSS(M,4) = VS*T0000*E1111*C9*B54
     &                - V8*T0001*E1110*C7*B34 - V8*T0010*E1101*C7*B34
     &                - V8*T0100*E1011*C7*B52 - V8*T1000*E0111*C7*B52
     &                + V4*T1100*E0011*C5*B50 + V4*T0011*E1100*C5*B14
     &                + V4*T0110*E1001*C5*B32 + V4*T1001*E0110*C5*B32
     &                + V4*T1010*E0101*C5*B32 + V4*T0101*E1010*C5*B32
     &                - V2*T0111*E1000*C3*B12 - V2*T1011*E0100*C3*B12
     &                - V2*T1101*E0010*C3*B30 - V2*T1110*E0001*C3*B30
     &                + V1*T1111*E0000*C1*B10
104       CONTINUE
C
        ENDIF
C
C**********************************************************************C
C       EXCHANGE INTEGRAL MATRICES: RKSSSS, RKSLSL, RKLLLL             C
C**********************************************************************C
C
C       LOOP OVER THE TENSOR ORDERS OF THE COULOMB INTERACTION
        DO LTEN=1,NUNUM
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(LTEN)
C
          IF(HMLT.EQ.'NORL') THEN
C         NON-RELATIVISTIC HAMILTONIAN
C
            B32 = EIK(M,-NU+2)*EJL(M, NU+3)*BXL(M, NU+3,-NU+2)
     &          + EIK(M, NU+3)*EJL(M,-NU+2)*BXU(M, NU+3,-NU+2)
C
            RKLL = V1*T0000*E0000*C5*B32
C
            RKLLLL(M,1) = RKLLLL(M,1) + BK(LTEN,1)*RKLL
C
          ELSE
C         RELATIVISTIC HAMILTONIAN
C
C           TEMPORARY STORAGE OF RAW RK(LTEN,M)
            B10 = EIK(M,-NU  )*EJL(M, NU+1)*BXL(M, NU+1,-NU  )
     &          + EIK(M, NU+1)*EJL(M,-NU  )*BXU(M, NU+1,-NU  )
            B12 = EIK(M,-NU+2)*EJL(M, NU+1)*BXL(M, NU+1,-NU+2)
     &          + EIK(M, NU+3)*EJL(M,-NU  )*BXU(M, NU+3,-NU  )
            B14 = EIK(M,-NU+4)*EJL(M, NU+1)*BXL(M, NU+1,-NU+4)
     &          + EIK(M, NU+5)*EJL(M,-NU  )*BXU(M, NU+5,-NU  )
            B30 = EIK(M,-NU  )*EJL(M, NU+3)*BXL(M, NU+3,-NU  )
     &          + EIK(M, NU+1)*EJL(M,-NU+2)*BXU(M, NU+1,-NU+2)
            B32 = EIK(M,-NU+2)*EJL(M, NU+3)*BXL(M, NU+3,-NU+2)
     &          + EIK(M, NU+3)*EJL(M,-NU+2)*BXU(M, NU+3,-NU+2)
            B34 = EIK(M,-NU+4)*EJL(M, NU+3)*BXL(M, NU+3,-NU+4)
     &          + EIK(M, NU+5)*EJL(M,-NU+2)*BXU(M, NU+5,-NU+2)
            B50 = EIK(M,-NU  )*EJL(M, NU+5)*BXL(M, NU+5,-NU  )
     &          + EIK(M, NU+1)*EJL(M,-NU+4)*BXU(M, NU+1,-NU+4)
            B52 = EIK(M,-NU+2)*EJL(M, NU+5)*BXL(M, NU+5,-NU+2)
     &          + EIK(M, NU+3)*EJL(M,-NU+4)*BXU(M, NU+3,-NU+4)
            B54 = EIK(M,-NU+4)*EJL(M, NU+5)*BXL(M, NU+5,-NU+4)
     &          + EIK(M, NU+5)*EJL(M,-NU+4)*BXU(M, NU+5,-NU+4)
C
C           KQNA < 0 AND KQNB < 0 CONTRIBUTIONS (CANNOT SKIP)
            RKLL = V1*T0000*E0000*C5*B32
            RKSL = V4*T0000*E1010*C7*B34
            RKSS = VS*T0000*E1111*C9*B54
C
            RKLLLL(M,1) = RKLLLL(M,1) + BK(LTEN,1)*RKLL
            RKSLSL(M,1) = RKSLSL(M,1) + BK(LTEN,1)*RKSL
            RKSSSS(M,1) = RKSSSS(M,1) + BK(LTEN,1)*RKSS
C
C           KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNB.EQ.0) GOTO 202
            RKLL = V1*T0000*E0000*C5*B32
            RKSL = V4*T0000*E1010*C7*B34 - V2*T0010*E1000*C5*B32
            RKSS = VS*T0000*E1111*C9*B54 - V8*T0001*E1110*C7*B34
     &           - V8*T0010*E1101*C7*B52 + V4*T0011*E1100*C5*B32
C
            RKLLLL(M,2) = RKLLLL(M,2) + BK(LTEN,2)*RKLL
            RKSLSL(M,2) = RKSLSL(M,2) + BK(LTEN,2)*RKSL
            RKSSSS(M,2) = RKSSSS(M,2) + BK(LTEN,2)*RKSS
202         CONTINUE
C
C           KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNA.EQ.0) GOTO 203
            RKLL = V1*T0000*E0000*C5*B32
            RKSL = V4*T0000*E1010*C7*B34 - V2*T1000*E0010*C5*B32
            RKSS = VS*T0000*E1111*C9*B54 - V8*T0100*E1011*C7*B34
     &           - V8*T1000*E0111*C7*B52 + V4*T1100*E0011*C5*B32
C
            RKLLLL(M,3) = RKLLLL(M,3) + BK(LTEN,3)*RKLL
            RKSLSL(M,3) = RKSLSL(M,3) + BK(LTEN,3)*RKSL
            RKSSSS(M,3) = RKSSSS(M,3) + BK(LTEN,3)*RKSS
203         CONTINUE
C
C           KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
            IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 204
            RKLL = V1*T0000*E0000*C5*B32
            RKSL = V4*T0000*E1010*C7*B34 - V2*T1000*E0010*C5*B32
     &           - V2*T0010*E1000*C5*B32 + V1*T1010*E0000*C3*B30
            RKSS = VS*T0000*E1111*C9*B54
     &           - V8*T0001*E1110*C7*B34 - V8*T0010*E1101*C7*B52
     &           - V8*T0100*E1011*C7*B34 - V8*T1000*E0111*C7*B52
     &           + V4*T1100*E0011*C5*B32 + V4*T0011*E1100*C5*B32
     &           + V4*T1001*E0110*C5*B32 + V4*T0110*E1001*C5*B32
     &           + V4*T0101*E1010*C5*B14 + V4*T1010*E0101*C5*B50
     &           - V2*T1101*E0010*C3*B12 - V2*T0111*E1000*C3*B12
     &           - V2*T1110*E0001*C3*B30 - V2*T1011*E0100*C3*B30
     &           + V1*T1111*E0000*C1*B10
C
            RKLLLL(M,4) = RKLLLL(M,4) + BK(LTEN,4)*RKLL
            RKSLSL(M,4) = RKSLSL(M,4) + BK(LTEN,4)*RKSL
            RKSSSS(M,4) = RKSSSS(M,4) + BK(LTEN,4)*RKSS
204         CONTINUE
C
          ENDIF
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     APPLY NORMALISATION FACTORS                                      C
C**********************************************************************C
C
      IF(HMLT.EQ.'NORL') THEN
C     NON-RELATIVISTIC HAMILTONIAN
C
        DO M=1,MAXM
          RNLLLL = RNIJ(1)*RNKL(M,1)
          RJLLLL(M,1) = RJLLLL(M,1)*RNLLLL
          RKLLLL(M,1) = RKLLLL(M,1)*RNLLLL
        ENDDO
C
      ELSE
C     RELATIVISTIC HAMILTONIAN
C
        DO M=1,MAXM
          RNLLLL = RNIJ(1)*RNKL(M,1)
          RNLLSS = RNIJ(1)*RNKL(M,4)
          RNSLSL = RNIJ(3)*RNKL(M,3)
          RNSSLL = RNIJ(4)*RNKL(M,1)
          RNSSSS = RNIJ(4)*RNKL(M,4)
          DO N=1,4
            RJLLLL(M,N) = RJLLLL(M,N)*RNLLLL
            RJLLSS(M,N) = RJLLSS(M,N)*RNLLSS
            RJSSLL(M,N) = RJSSLL(M,N)*RNSSLL
            RJSSSS(M,N) = RJSSSS(M,N)*RNSSSS
            RKLLLL(M,N) = RKLLLL(M,N)*RNLLLL
            RKSLSL(M,N) = RKSLSL(M,N)*RNSLSL
            RKSSSS(M,N) = RKSSSS(M,N)*RNSSSS
          ENDDO
        ENDDO
C
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE BREIT0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           BBBBBBB  RRRRRRR  EEEEEEEE IIII TTTTTTTT 000000            C
C           BB    BB RR    RR EE        II     TT   00   000           C
C           BB    BB RR    RR EE        II     TT   00  0000           C
C           BBBBBBB  RR    RR EEEEEE    II     TT   00 00 00           C
C           BB    BB RRRRRRR  EE        II     TT   0000  00           C
C           BB    BB RR    RR EE        II     TT   000   00           C
C           BBBBBBB  RR    RR EEEEEEEE IIII    TT    000000            C
C                                                                      C
C -------------------------------------------------------------------- C
C  BREIT0 CONSTRUCTS THE ATOMIC BREIT MATRIX FROM RADIAL DIRECT AND    C
C  EXCHANGE INTEGRALS AND A MEAN-FIELD CHARGE DENSITY.                 C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RKLLSS(MB2,4),RKSLLS(MB2,4),RKSSLL(MB2,4),RMSLLS(MB2,4)
C
      COMMON/ATMB/B11(MBD,MBD),B21(MBD,MBD),B12(MBD,MBD),B22(MBD,MBD)
      COMMON/ATMD/DLL1(MB2),DSL1(MB2),DSS1(MB2),DLS1(MB2),
     &            DLL2(MB2),DSL2(MB2),DSS2(MB2),DLS2(MB2)
      COMMON/B0IJ/EIJ(-MTN:MTN),RNIJ(4),EI,EJ,IBAS,JBAS
      COMMON/B0QN/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/RCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
C
C     INITIALISE BREIT MATRIX
      DO IBAS=1,MBD
        DO JBAS=1,MBD
          B11(IBAS,JBAS) = 0.0D0
          B21(IBAS,JBAS) = 0.0D0
          B12(IBAS,JBAS) = 0.0D0
          B22(IBAS,JBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     COEFFICIENTS FOR INTEGRAL ASSEMBLY
      C3 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+3)
      C5 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+5)
      C7 = 0.25D0*GAMHLF(2*LQNA+2*LQNB+7)
C
      V1 = 1.0D0
      V2 = 2.0D0
      V4 = 4.0D0
C
      TI = DFLOAT(2*LQNA+1)
      TJ = DFLOAT(2*LQNA+1)
      TK = DFLOAT(2*LQNB+1)
      TL = DFLOAT(2*LQNB+1)
C
      T0000 = 1.0D0
      T1000 = TI
      T0100 = TJ
      T0010 = TK
      T0001 = TL
      T1100 = TI*TJ
      T1001 = TI*TL
      T0011 = TK*TL
C
C     EVALUATE CLOSED-SHELL BREIT INTERACTION ANGULAR INTEGRALS
      CALL ANGBRT0
C
C     BASIS SET INTERMEDIATES FOR THE KL-PAIRS
      CALL KLSET0
C
C     ITERATE OVER ALL MATRIX ELEMENTS
      DO IBAS=1,NBASA
        DO JBAS=1,NBASA
C
C         GAUSSIAN EXPONENTS FOR THIS PAIR
          EI = EXLA(IBAS)
          EJ = EXLA(JBAS)
C
C         BASIS SET INTERMEDIATES FOR THE IJ-PAIRS
          CALL IJSET0
C
C         GENERATE BATCH OF RADIAL INTEGRALS (J AND K MATRICES)
          CALL RKBRT0(RKLLSS,RKSLLS,RKSSLL,RMSLLS)
C
C         SMALL-COMPONENT MATRIX ADDRESSES
          KBAS = IBAS + NBASA
          LBAS = JBAS + NBASA
C
C    (22) KQNA < 0 AND KQNB < 0 CONTRIBUTIONS (CANNOT SKIP)
C
C         INITIALISE COUNTERS
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0
          
          KAPA2 =-LQNA-1
          KAPB2 =-LQNB-1
          RK2A2 = DFLOAT(2*IABS(KAPA2))
          RK2B2 = DFLOAT(2*IABS(KAPB2))
C
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,1)*DSS2(M)
            BSL = BSL + RKSLLS(M,1)*DLS2(M) + RMSLLS(M,1)*DLS2(M)
            BSS = BSS + RKSSLL(M,1)*DLL2(M)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B22(IBAS,JBAS) = BLL
          B22(KBAS,JBAS) = BSL
          B22(JBAS,KBAS) = BSL
          B22(KBAS,LBAS) = BSS
C
C    (21) KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNB.EQ.0) GOTO 200
C
C         INITIALISE COUNTERS
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0

          KAPA2 =-LQNA-1
          KAPB1 = LQNB
          RK2A2 = DFLOAT(2*IABS(KAPA2))
          RK2B1 = DFLOAT(2*IABS(KAPB1))
C
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,2)*DSS1(M)
            BSL = BSL + RKSLLS(M,2)*DLS1(M) + RMSLLS(M,2)*DLS1(M)
            BSS = BSS + RKSSLL(M,2)*DLL1(M)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B21(IBAS,JBAS) = BLL
          B21(KBAS,JBAS) = BSL
          B21(JBAS,KBAS) = BSL
          B21(KBAS,LBAS) = BSS
C
200       CONTINUE
C
C    (12) KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0) GOTO 300
C
C         INITIALISE COUNTERS
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0

          KAPA1 = LQNA
          KAPB2 =-LQNB-1
          RK2A1 = DFLOAT(2*IABS(KAPA1))
          RK2B2 = DFLOAT(2*IABS(KAPB2))
C
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,3)*DSS2(M)
            BSL = BSL + RKSLLS(M,3)*DLS2(M) + RMSLLS(M,3)*DLS2(M)
            BSS = BSS + RKSSLL(M,3)*DLL2(M)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B12(IBAS,JBAS) = BLL
          B12(KBAS,JBAS) = BSL
          B12(JBAS,KBAS) = BSL
          B12(KBAS,LBAS) = BSS
C
300       CONTINUE
C
C    (11) KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 400
C
C         INITIALISE COUNTERS
          BLL = 0.0D0
          BSL = 0.0D0
          BSS = 0.0D0

          KAPA1 = LQNA
          KAPB1 = LQNB
          RK2A1 = DFLOAT(2*IABS(KAPA1))
          RK2B1 = DFLOAT(2*IABS(KAPB1))
C
C         SUM OVER MEAN FIELD CONTRIBUTIONS FOR THIS BASIS PAIR
          DO M=1,MAXM
            BLL = BLL + RKLLSS(M,4)*DSS1(M)
            BSL = BSL + RKSLLS(M,4)*DLS1(M) + RMSLLS(M,4)*DLS1(M)
            BSS = BSS + RKSSLL(M,4)*DLL1(M)
          ENDDO
C
C         TRANSFER COUNTER VALUES TO BREIT MATRIX
          B11(IBAS,JBAS) = BLL
          B11(KBAS,JBAS) = BSL
          B11(JBAS,KBAS) = BSL
          B11(KBAS,LBAS) = BSS
C
400       CONTINUE
C
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE RKBRT0(RKLLSS,RKSLLS,RKSSLL,RMSLLS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         RRRRRRR  KK    KK BBBBBBB  RRRRRRR TTTTTTTT 000000           C
C         RR    RR KK   KK  BB    BB RR    RR   TT   00   000          C
C         RR    RR KK  KK   BB    BB RR    RR   TT   00  0000          C
C         RR    RR KKKKK    BBBBBBB  RR    RR   TT   00 00 00          C
C         RRRRRRR  KK  KK   BB    BB RRRRRRR    TT   0000  00          C
C         RR    RR KK   KK  BB    BB RR    RR   TT   000   00          C
C         RR    RR KK    KK BBBBBBB  RR    RR   TT    000000           C
C                                                                      C
C -------------------------------------------------------------------- C
C  RKBRT0 EVALUATES A DIRECT AND EXCHANGE BATCH OF BREIT INTERACTION   C
C  INTEGRALS OF ALL COMPONENT LABEL COMBINATIONS L AND S IN THE ATOMIC C
C  (RELATIVISTIC) SCF PROCEDURE.                                       C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C  ▶ RKLLSS(M,N) - EXCHANGE BII OVERLAP {LL,SS}                        C
C  ▶ RKSLLS(M,N) - EXCHANGE BII OVERLAP {SL,SL}                        C
C  ▶ RKSSLL(M,N) - EXCHANGE BII OVERLAP {SS,LL}                        C
C  ▶ RMSLLS(M,N) - SEMI-RANGE BII OVERLAP {SL,SL}                      C
C -------------------------------------------------------------------- C
C  BREIT MATRIX LIST INDEX:                                            C
C  ▶ N=1 - KQN(A)<0, KQN(B)<0 (TYPICAL LABEL 22)                       C
C  ▶ N=2 - KQN(A)<0, KQN(B)>0 (TYPICAL LABEL 12)                       C
C  ▶ N=3 - KQN(A)>0, KQN(B)<0 (TYPICAL LABEL 21)                       C
C  ▶ N=4 - KQN(A)>0, KQN(B)>0 (TYPICAL LABEL 11)                       C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION XK(MB2,2),IAA(2),IBB(2)
      DIMENSION RTIK0(MBS),RTJL0(MBS),PTIK0(MBS),PTJL0(MBS)
      DIMENSION BETA(MB2),BTA1(MB2),XROOT(MB2),BIN(MB2),TRM(MB2)
      DIMENSION BXU(MB2,-MTN:MTN,-MTN:MTN),BXL(MB2,-MTN:MTN,-MTN:MTN)
      DIMENSION RKLLSS(MB2,4),RKSLLS(MB2,4),RKSSLL(MB2,4),RMSLLS(MB2,4)
C
      COMMON/B0IJ/EIJ(-MTN:MTN),RNIJ(4),EI,EJ,IBAS,JBAS
      COMMON/B0IK/EIK(MB2,-MTN:MTN),IKIND(MB2)
      COMMON/B0JL/EJL(MB2,-MTN:MTN),JLIND(MB2)
      COMMON/B0KL/EKL(MB2,-MTN:MTN),RNKL(MB2,4),EK(MB2),EL(MB2)
      COMMON/B0QN/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/RCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
      COMMON/XTNS/BK(MNU,4),ELL(MNU,4),ESS(MNU,4),ESL(MNU,4),GSL(MNU,4)
C
C**********************************************************************C
C     PREPARATION OF EXPONENT POWERS AND BETA FUNCTION ARGUMENTS.      C
C**********************************************************************C
C
C     TENSOR ORDER LIMITS
      NUI = NUS(1)
      NUF = NUS(NUNUM)
C
C     GAUSSIAN EXPONENT FOR (IBAS,JBAS) OVERLAP
      TIJ0 = EI+EJ
C
C     INITIALISE THE ARRAY XK(M,2) FOR INCOMPLETE BETA FUNCTION ARGS.
      DO M=1,MAXM
        TKL0  = EK(M)+EL(M)
        TIJKL = TIJ0+TKL0
        XK(M,1) = (EI+EK(M))/TIJKL
        XK(M,2) = (EJ+EL(M))/TIJKL
      ENDDO
C
C     LOWEST EXPONENT POWER
      IPOWER = LQNA+LQNB-NUF
C
C     A BLOCK OF BASIS EXPONENT PRODUCTS
      DO KBAS=1,NBASB
        RTIK0(KBAS) = DSQRT(EI+EXLB(KBAS))
        RTJL0(KBAS) = DSQRT(EJ+EXLB(KBAS))
        PTIK0(KBAS) = RTIK0(KBAS)**(-IPOWER)
        PTJL0(KBAS) = RTJL0(KBAS)**(-IPOWER)
      ENDDO
C
C     CALCULATE A FULL SET OF EXPONENT OVERLAPS FOR EXCHANGE
      DO M=1,MAXM
        RTIK = RTIK0(IKIND(M))
        RTJL = RTJL0(JLIND(M))
        EIK(M,-NUF) = PTIK0(IKIND(M))
        EJL(M,-NUF) = PTJL0(JLIND(M))
        DO IPOW=-NUF+1,NUF+5
          EIK(M,IPOW) = EIK(M,IPOW-1)/RTIK
          EJL(M,IPOW) = EJL(M,IPOW-1)/RTJL
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INLINE BETA FUNCTION CODE FOR EXCHANGE TERMS                     C
C**********************************************************************C
C
C     NUMBER OF LEVELS NEEDED TO ACCOUNT FOR ALL EXCHANGE INTEGRALS
      NVALS = (NUF-NUI)/2+2
C
C     LOOP OVER ORDER FOR FIRST INDEX
      DO NX=1,NVALS
        IAA(1) = LQNA+LQNB+NUI+2*NX
        IAA(2) = LQNA+LQNB+NUI+2*NX
C
C       LOOP OVER ORDER FOR SECOND INDEX
        DO NY=1,NVALS
          IBB(1) = LQNA+LQNB-NUF+2*NY-1
          IBB(2) = LQNA+LQNB-NUF+2*NY-1
C
C         LOOP OVER BETA INTEGRAL TYPE
          DO IBETA=1,2
            IA = (IAA(IBETA)-1)/2
            IB =  IBB(IBETA)   /2
C
C           CASE 1: IB > 1
            IF(IB.GT.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BTA1(M) = (DSQRT(XK(M,IBETA))**IX)/X
              ENDDO
              RA = X
              RB = DFLOAT(1-IB)
              RC = 1.0D0+X
              RD = 1.0D0
              RCD = RC*RD
              FCT = RA*RB/RCD
              DO M=1,MAXM
                TRM(M) = FCT*XK(M,IBETA)
                BIN(M) = 1.0D0+TRM(M)
              ENDDO
              RA = RA+1.0D0
              RB = RB+1.0D0
              RC = RC+1.0D0
              RD = RD+1.0D0
              DO J=2,IB-1
                RCD = RC*RD
                FCT = RA*RB/RCD
                DO M=1,MAXM
                  TRM(M) = FCT*TRM(M)*XK(M,IBETA)
                  BIN(M) = BIN(M)+TRM(M)
                ENDDO
                RA = RA+1.0D0
                RB = RB+1.0D0
                RC = RC+1.0D0
                RD = RD+1.0D0
              ENDDO
              DO M=1,MAXM
                BETA(M) = BTA1(M)*BIN(M)
              ENDDO
C
C           CASE 2: IB = 1
            ELSEIF(IB.EQ.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BETA(M) = (DSQRT(XK(M,IBETA))**IX)/X
              ENDDO
C
C           CASE 3: IB = 0
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
C             CASE A: IA > 1
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
C             CASE B: IA = 1
              ELSEIF(IA.EQ.1) THEN
                DO M=1,MAXM
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)
                ENDDO
C             CASE C: IA = 0
              ELSE
                DO M=1,MAXM
                  BETA(M) = BTA1(M)
                ENDDO
              ENDIF
C
            ENDIF
C
C           SORT THE BETA INTEGRAL INTO THE UPPER/LOWER ARRAYS
            IF(IBETA.EQ.1) THEN
              DO M=1,MAXM
                BXU(M, NUI+2*NX  ,-NUF+2*NY-1) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXM
                BXL(M, NUI+2*NX  ,-NUF+2*NY-1) = BETA(M)
              ENDDO
            ENDIF
C
C         END LOOPS OVER IBETA AND INDICES NX, NY
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     RADIAL INTEGRALS OVER SPINORS (SEPARATED BY TENSOR ORDER)        C
C**********************************************************************C
C
C     EMPTY COUNTER ARRAYS FOR DIRECT AND EXCHANGE INTEGRALS
      DO M=1,MAXM
        DO N=1,4
          RKLLSS(M,N) = 0.0D0
          RKSLLS(M,N) = 0.0D0
          RKSSLL(M,N) = 0.0D0
          RMSLLS(M,N) = 0.0D0
        ENDDO
      ENDDO
C
      E0000 = 1.0D0
      E1000 = EI
      E0100 = EJ
      E1100 = EI*EJ
C
C     LOOP OVER K,L BASIS FUNCTIONS
      DO M=1,MAXM
C
C       MORE VALUE PREPARATION
        E0010 = EK(M)
        E0001 = EL(M)
        E1001 = EI*EL(M)
        E0011 = EK(M)*EL(M)
C
C**********************************************************************C
C       EXCHANGE INTEGRAL MATRICES: RKLLSS, RKSLLS, RKSSLL             C
C**********************************************************************C
C
C       LOOP OVER THE TENSOR ORDERS OF THE BREIT INTERACTION
        DO LTEN=1,NUNUM
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(LTEN)
C
C         TEMPORARY STORAGE OF RAW RK(LTEN,M)
          B21 = EIK(M,-NU+1)*EJL(M, NU+2)*BXL(M, NU+2,-NU+1)
     &        + EIK(M, NU+2)*EJL(M,-NU+1)*BXU(M, NU+2,-NU+1)
          B23 = EIK(M,-NU+3)*EJL(M, NU+2)*BXL(M, NU+2,-NU+3)
     &        + EIK(M, NU+4)*EJL(M,-NU+1)*BXU(M, NU+4,-NU+1)
          B41 = EIK(M,-NU+1)*EJL(M, NU+4)*BXL(M, NU+4,-NU+1)
     &        + EIK(M, NU+2)*EJL(M,-NU+3)*BXU(M, NU+2,-NU+3)
          B43 = EIK(M,-NU+3)*EJL(M, NU+4)*BXL(M, NU+4,-NU+3)
     &        + EIK(M, NU+4)*EJL(M,-NU+3)*BXU(M, NU+4,-NU+3)
C
C         FILL RK ARRAYS FOR THIS LQNA,LQNB BLOCK
C
C         KQNA < 0 AND KQNB < 0 CONTRIBUTIONS (CANNOT SKIP)
          RKLL = V4*T0000*E0011*C7*B43
          RKSL = V4*T0000*E1001*C7*B43
          RKSS = V4*T0000*E1100*C7*B43
C
          RKLLSS(M,1) = RKLLSS(M,1) + ELL(LTEN,1)*RKLL
          RKSLLS(M,1) = RKSLLS(M,1) + ESL(LTEN,1)*RKSL
          RKSSLL(M,1) = RKSSLL(M,1) + ESS(LTEN,1)*RKSS
C
C         KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNB.EQ.0) GOTO 202
          RKLL = V4*T0000*E0011*C7*B43 - V2*T0001*E0010*C5*B23
     &         - V2*T0010*E0001*C5*B41 + V1*T0011*E0000*C3*B21
          RKSL = V4*T0000*E1001*C7*B43 - V2*T0001*E1000*C5*B23
          RKSS = V4*T0000*E1100*C7*B43
C
          RKLLSS(M,2) = RKLLSS(M,2) + ELL(LTEN,2)*RKLL
          RKSLLS(M,2) = RKSLLS(M,2) + ESL(LTEN,2)*RKSL
          RKSSLL(M,2) = RKSSLL(M,2) + ESS(LTEN,2)*RKSS
202       CONTINUE
C
C         KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0) GOTO 203
          RKLL = V4*T0000*E0011*C7*B43
          RKSL = V4*T0000*E1001*C7*B43 - V2*T1000*E0001*C5*B41
          RKSS = V4*T0000*E1100*C7*B43 - V2*T0100*E1000*C5*B23
     &         - V2*T1000*E0100*C5*B41 + V1*T1100*E0000*C3*B21
C
          RKLLSS(M,3) = RKLLSS(M,3) + ELL(LTEN,3)*RKLL
          RKSLLS(M,3) = RKSLLS(M,3) + ESL(LTEN,3)*RKSL
          RKSSLL(M,3) = RKSSLL(M,3) + ESS(LTEN,3)*RKSS
203       CONTINUE
C
C         KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 204
          RKLL = V4*T0000*E0011*C7*B43 - V2*T0001*E0010*C5*B23
     &         - V2*T0010*E0001*C5*B41 + V1*T0011*E0000*C3*B21
          RKSL = V4*T0000*E1001*C7*B43 - V2*T1000*E0001*C5*B41
     &         - V2*T0001*E1000*C5*B23 + V1*T1001*E0000*C3*B21
          RKSS = V4*T0000*E1100*C7*B43 - V2*T0100*E1000*C5*B23
     &         - V2*T1000*E0100*C5*B41 + V1*T1100*E0000*C3*B21
C
          RKLLSS(M,4) = RKLLSS(M,4) + ELL(LTEN,4)*RKLL
          RKSLLS(M,4) = RKSLLS(M,4) + ESL(LTEN,4)*RKSL
          RKSSLL(M,4) = RKSSLL(M,4) + ESS(LTEN,4)*RKSS
204       CONTINUE
C
        ENDDO
C
C**********************************************************************C
C       HALF-RANGE EXCHANGE INTEGRAL MATRICES: RMSLLS                  C
C**********************************************************************C
C
C       LOOP OVER THE TENSOR ORDERS OF THE COULOMB INTERACTION
        DO LTEN=1,NUNUM
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(LTEN)
C
C         TEMPORARY STORAGE OF RAW RM(LTEN,M) FOR RMSLLS
          B21 = EIK(M,-NU+1)*EJL(M, NU+2)*BXL(M, NU+2,-NU+1)
     &        - EIK(M, NU+2)*EJL(M,-NU+1)*BXU(M, NU+2,-NU+1)
          B23 = EIK(M,-NU+3)*EJL(M, NU+2)*BXL(M, NU+2,-NU+3)
     &        - EIK(M, NU+4)*EJL(M,-NU+1)*BXU(M, NU+4,-NU+1)
          B41 = EIK(M,-NU+1)*EJL(M, NU+4)*BXL(M, NU+4,-NU+1)
     &        - EIK(M, NU+2)*EJL(M,-NU+3)*BXU(M, NU+2,-NU+3)
          B43 = EIK(M,-NU+3)*EJL(M, NU+4)*BXL(M, NU+4,-NU+3)
     &        - EIK(M, NU+4)*EJL(M,-NU+3)*BXU(M, NU+4,-NU+3)
C
C         FILL RM ARRAYS FOR THIS LQNA,LQNB BLOCK
C
C         KQNA < 0 AND KQNB < 0 CONTRIBUTIONS (CANNOT SKIP)
          RMSL = V4*T0000*E1001*C7*B43

          RMSLLS(M,1) = RMSLLS(M,1) + GSL(LTEN,1)*RMSL
C
C         KQNA < 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNB.EQ.0) GOTO 302
          RMSL = V4*T0000*E1001*C7*B43 - V2*T0001*E1000*C5*B23
C
          RMSLLS(M,2) = RMSLLS(M,2) + GSL(LTEN,2)*RMSL
302       CONTINUE
C
C         KQNA > 0 AND KQNB < 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0) GOTO 303
          RMSL = V4*T0000*E1001*C7*B43 - V2*T1000*E0001*C5*B41
C
          RMSLLS(M,3) = RMSLLS(M,3) + GSL(LTEN,3)*RMSL
303       CONTINUE
C
C         KQNA > 0 AND KQNB > 0 CONTRIBUTIONS (SKIP IF POSSIBLE)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 304
          RMSL = V4*T0000*E1001*C7*B43 - V2*T1000*E0001*C5*B41
     &         - V2*T0001*E1000*C5*B23 + V1*T1001*E0000*C3*B21
C
          RMSLLS(M,4) = RMSLLS(M,4) + GSL(LTEN,4)*RMSL
304       CONTINUE
C
        ENDDO
C
      ENDDO
C
C**********************************************************************C
C     APPLY NORMALISATION FACTORS                                      C
C**********************************************************************C
C
      DO M=1,MAXM
        RNLLSS = RNIJ(1)*RNKL(M,4)
        RNSSLL = RNIJ(4)*RNKL(M,1)
        RNSLLS = RNIJ(3)*RNKL(M,2)
        DO N=1,4
          RKLLSS(M,N) = RNLLSS*RKLLSS(M,N)
          RKSSLL(M,N) = RNSSLL*RKSSLL(M,N)
          RKSLLS(M,N) = RNSLLS*RKSLLS(M,N)
          RMSLLS(M,N) = RNSLLS*RMSLLS(M,N)
        ENDDO
      ENDDO
C
      RETURN
      END
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
      FUNCTION ZCVINT0(IZ,L,EIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       ZZZZZZZZ  CCCCCC  VV    VV IIII NN    NN TTTTTTTT 000000       C
C            ZZ  CC    CC VV    VV  II  NNN   NN    TT   00   000      C
C           ZZ   CC       VV    VV  II  NNNN  NN    TT   00  0000      C
C          ZZ    CC       VV    VV  II  NN NN NN    TT   00 00 00      C
C         ZZ     CC        VV  VV   II  NN  NNNN    TT   0000  00      C
C        ZZ      CC    CC   VVVV    II  NN   NNN    TT   000   00      C
C       ZZZZZZZZ  CCCCCC     VV    IIII NN    NN    TT    000000       C
C                                                                      C
C -------------------------------------------------------------------- C
C  ZCVINT0 CALCULATES AN NUCLEAR POTENTIAL CORRECTION MATRIX ELEMENT   C
C  OVER AN ATOMIC BASIS FUNCTION PRODUCT, WHERE L>-1. THE POTENTIAL IS C
C  A COMBINATION OF ZERO-CHARGE POTENTIAL BASIS FUNCTIONS.             C
C                        ∞                                             C
C        ZCVINT0(IZ,L,λ) = ∫ r^2L+2 exp(-λ r^2) V_nuc(r;IZ) dr         C
C                        0                                             C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5 NMDL
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In ZCVINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In ZCVINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
      ZCVINT0 = 0.0D0
      IF(NNUC(IZ).EQ.0) RETURN
C
      DO IFT=1,NNUC(IZ)
C
        XI  = XNUC(IZ,IFT)
        FC  = FNUC(IZ,IFT)
C
        PES = XI+EIJ
        RGN = 1.0D0/(PES*DSQRT(PES))
        XEL = PES**L
C
        ZCVINT0 = ZCVINT0 + 0.5D0*FC*GAMHLF(2*L+3)*RGN/XEL
C
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION GSSINT0(L,EIJ,XI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       GGGGGG   SSSSSS   SSSSSS  IIII NN    NN TTTTTTTT 000000        C
C      GG    GG SS    SS SS    SS  II  NNN   NN    TT   00   000       C
C      GG       SS       SS        II  NNNN  NN    TT   00  0000       C
C      GG        SSSSSS   SSSSSS   II  NN NN NN    TT   00 00 00       C
C      GG   GGG       SS       SS  II  NN  NNNN    TT   0000  00       C
C      GG    GG SS    SS SS    SS  II  NN   NNN    TT   000   00       C
C       GGGGGG   SSSSSS   SSSSSS  IIII NN    NN    TT    000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  GSSINT0 CALCULATES A RADIAL INTEGRAL OVER A POLYNOMIAL IN R AND     C
C  TWO GAUSSIANS -- EXPONENTS EIJ AND XI.                              C
C                              ∞                                       C
C             GSSINT(L,λ,ζ) = ∫ r^2L+2 exp(-(λ+ζ) r^2) dr.             C
C                              0                                       C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In GSSINT0: order L must be positive. L = ',L
        WRITE(7, *) 'In GSSINT0: order L must be positive. L = ',L
        STOP
      ENDIF
C
      XPR = XI*DSQRT(XI)/PI32
      PES = XI+EIJ
      RGN = XPR/(PES*DSQRT(PES))
      XEL = PES**L
C
C     TRANSFER DATA TO GSSINT0
      GSSINT0 = 0.5D0*GAMHLF(2*L+3)*RGN/XEL
C
      RETURN
      END
C
C
      FUNCTION GSSINT02(L,EIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       GGGGGG   SSSSSS   SSSSSS  IIII NN    NN TTTTTTTT 000000        C
C      GG    GG SS    SS SS    SS  II  NNN   NN    TT   00   000       C
C      GG       SS       SS        II  NNNN  NN    TT   00  0000       C
C      GG        SSSSSS   SSSSSS   II  NN NN NN    TT   00 00 00       C
C      GG   GGG       SS       SS  II  NN  NNNN    TT   0000  00       C
C      GG    GG SS    SS SS    SS  II  NN   NNN    TT   000   00       C
C       GGGGGG   SSSSSS   SSSSSS  IIII NN    NN    TT    000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  GSSINT0 CALCULATES A RADIAL INTEGRAL OVER A POLYNOMIAL IN R AND     C
C  TWO GAUSSIANS -- EXPONENTS EIJ AND XI.                              C
C                              ∞                                       C
C             GSSINT(L,λ,ζ) = ∫ r^2L+2 exp(-(λ+ζ) r^2) dr.             C
C                              0                                       C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In GSSINT02: order L must be positive. L = ',L
        WRITE(7, *) 'In GSSINT02: order L must be positive. L = ',L
        STOP
      ENDIF
C
      RGN = 1.0D0/(EIJ*DSQRT(EIJ))
      XEL = EIJ**L
C
C     TRANSFER DATA TO GSSINT0
      GSSINT02 = 0.5D0*GAMHLF(2*L+3)*RGN/XEL
C
      RETURN
      END
C
C
      FUNCTION DNSINT0(L,EIJ,DS,RS,EXS,NS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      DDDDDDD  NN    NN  SSSSSS  IIII NN    NN TTTTTTTT 000000        C
C      DD    DD NNN   NN SS    SS  II  NNN   NN    TT   00   000       C
C      DD    DD NNNN  NN SS        II  NNNN  NN    TT   00  0000       C
C      DD    DD NN NN NN  SSSSSS   II  NN NN NN    TT   00 00 00       C
C      DD    DD NN  NNNN       SS  II  NN  NNNN    TT   0000  00       C
C      DD    DD NN   NNN SS    SS  II  NN   NNN    TT   000   00       C
C      DDDDDDD  NN    NN  SSSSSS  IIII NN    NN    TT    000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  DNSINT0 CALCULATES A RADIAL INTEGRAL OVER A POLYNOMIAL IN R AND     C
C  TWO GAUSSIANS, WHERE THE SECOND IS ACTUALLY A SUM OF S-FUNCTIONS.   C
C                                 ∞                                    C
C          GSSINT(L,λ,ζ) = Σ D N ∫ r^2L+2 exp(-(λ+ζ) r^2) dr.          C
C                                 0                                    C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RS(MBS*MBS,4),DS(MBS*MBS),EXS(MBS)
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In DNSINT0: order L must be positive. L = ',L
        WRITE(7, *) 'In DNSINT0: order L must be positive. L = ',L
        STOP
      ENDIF
C
C     EMPTY THE BIN
      DNSINT0 = 0.0D0
C
C     LOOP OVER S-TYPE DENSITY
      N = 0
      DO KBAS=1,NS
        EK = EXS(KBAS)
        DO LBAS=1,NS
          N   = N+1
          EL  = EXS(LBAS)
          EKL = EK+EL
          DNSINT0 = DNSINT0 + DS(N)*RS(N,1)*GSSINT02(L,EIJ+EKL)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION ZCDINT0(IZ,L,EIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       ZZZZZZZZ  CCCCCC  DDDDDDD  IIII NN    NN TTTTTTTT 000000       C
C            ZZ  CC    CC DD    DD  II  NNN   NN    TT   00   000      C
C           ZZ   CC       DD    DD  II  NNNN  NN    TT   00  0000      C
C          ZZ    CC       DD    DD  II  NN NN NN    TT   00 00 00      C
C         ZZ     CC       DD    DD  II  NN  NNNN    TT   0000  00      C
C        ZZ      CC    CC DD    DD  II  NN   NNN    TT   000   00      C
C       ZZZZZZZZ  CCCCCC  DDDDDDD  IIII NN    NN    TT    000000       C
C                                                                      C
C -------------------------------------------------------------------- C
C  ZCDINT0 CALCULATES AN NUCLEAR DENSITY CORRECTION MATRIX ELEMENT     C
C  OVER AN ATOMIC BASIS FUNCTION PRODUCT, WHERE L>-1. THE DENSITY IS   C
C  A COMBINATION OF ZERO-CHARGE DENSITY BASIS FUNCTIONS.               C
C                        ∞                                             C
C       ZCDINT0(IZ,L,λ) = ∫ r^2L exp(-λ r^2) r^2*V_nuc(r;IZ) dr        C
C                        0                                             C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5 NMDL
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In ZCDINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In ZCDINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
      ZCDINT0 = 0.0D0
      IF(NNUC(IZ).EQ.0) RETURN
C
      DO IFT=1,NNUC(IZ)
C
        XI  = XNUC(IZ,IFT)
        FC  = FNUC(IZ,IFT)
C
        XPR = XI/PI
        PES = XI+EIJ
        RGN = XPR/(PES*PES*DSQRT(PES))
        XEL = PES**L
        BRK =-3.0D0*EIJ + 2.0D0*L*XI
C
        ZCDINT0 = ZCDINT0 + 0.25D0*FC*GAMHLF(2*L+3)*RGN/XEL
C
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION ZCEINT0(IZ,L,EIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       ZZZZZZZZ  CCCCCC  EEEEEEEE IIII NN    NN TTTTTTTT 000000       C
C            ZZ  CC    CC EE        II  NNN   NN    TT   00   000      C
C           ZZ   CC       EE        II  NNNN  NN    TT   00  0000      C
C          ZZ    CC       EEEEEE    II  NN NN NN    TT   00 00 00      C
C         ZZ     CC       EE        II  NN  NNNN    TT   0000  00      C
C        ZZ      CC    CC EE        II  NN   NNN    TT   000   00      C
C       ZZZZZZZZ  CCCCCC  EEEEEEEE IIII NN    NN    TT    000000       C
C                                                                      C
C -------------------------------------------------------------------- C
C  ZCEINT0 CALCULATES AN NUCLEAR ELECTRIC FIELD CORRECTION MATRIX      C
C  ELEMENT OVER AN ATOMIC BASIS FUNCTION PRODUCT, WHERE L>-1. THE      C
C  Er-FIELD IS A COMBINATION OF ZERO-CHARGE DENSITY BASIS FUNCTIONS.   C
C                        ∞                                             C
C      ZCEINT0(IZ,L,λ) = ∫ r^2L exp(-λ r^2) r^2*Er_nuc(r;IZ) dr        C
C                        0                                             C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5 NMDL
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In ZCEINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In ZCEINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
      ZCEINT0 = 0.0D0
      IF(NNUC(IZ).EQ.0) RETURN
C
      DO IFT=1,NNUC(IZ)
C
        XI  = XNUC(IZ,IFT)
        FC  = FNUC(IZ,IFT)
C
        XPR = XI/PI
        PES = XI+EIJ
        RGN = XPR/(PES*PES*DSQRT(PES))
        XEL = PES**L
        BRK =-3.0D0*EIJ + 2.0D0*L*XI
C
        ZCEINT0 = ZCEINT0 + 0.25D0*FC*XI*GAMHLF(2*L+3)*RGN/XEL
C
      ENDDO
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
      IF(SIGMA.LT.0.25D0) THEN
c       GOTO 123
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
      ELSE
c123     CONTINUE
        NNUMR = NNUMR+1
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
      ENDIF
C
C     COMBINE ALL TERMS AND APPLY NORMALISATION FACTOR
      FMIINT0 = FNRM*(EL1*Y1 + EL1*Y2 + C2L*Y34)
C
      RETURN
      END
C
C
      FUNCTION FMIINT0PM(L,EIJ,A,C)
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
C  THIS ROUTINE USES THE PARPIA AND MOHANTY RECIPE.                    C
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
      FUNCTION FMIINT0temp(L,EIJ,TPRM,CPRM,RPRM)
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
C  FMIINT0temp CALCULATES AN AUXILLIARY FERMI POTENTIAL INTEGRAL OVER  C
C  A PAIR OF ATOMIC BASIS FUNCTIONS, WHERE L IS A POSITIVE INTEGER.    C
C  IT DOES THIS BY QUADRATURE, SUBTRACTING OFF THE GAUSSIAN POTENTIAL. C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In FMIINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In FMIINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
C     NUMBER OF DATA POINTS IN UNIFORMLY-SPACED AND EXPONENTIAL REGION
      NLIN = 2000
C
C     GENERATE RADIAL GRID (LINEAR FROM FEMTOMETERS, EXPONENTIAL IN AU)
      RORI = 0.0D0
      RMID = 6.0D0*RPRM
C
      HL = (RMID-RORI)/DFLOAT(NLIN)
C
C     UNIFORMLY-SPACED GRID VALUES
      X0M = 0.0D0
      DO N=0,NLIN
        R   = RORI + HL*DFLOAT(N)
        Z1  = R**(2*L+2)
        Z2  = DEXP(-EIJ*R*R)
        VF  = VNFERMI(R,CPRM,TPRM)
        VG  = VNGAUSS(R,1.5D0/(RPRM*RPRM),1.0D0)
        X0M = X0M + EXTINT11(Z1*Z2*(VF-VG),N,NLIN)
      ENDDO
      X0M = 5.0D0*HL*X0M/299376.0D0
C        
      FMIINT0temp =-X0M
C
      RETURN
      END
C
C
      FUNCTION EXTINT2(Y,I,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      EEEEEEEE XX     XX TTTTTTTT IIII NN    NN TTTTTTTT 222222       C
C      EE        XX   XX     TT     II  NNN   NN    TT   22    22      C
C      EE         XX XX      TT     II  NNNN  NN    TT         22      C
C      EEEEEE      XXX       TT     II  NN NN NN    TT       22        C
C      EE         XX XX      TT     II  NN  NNNN    TT     22          C
C      EE        XX   XX     TT     II  NN   NNN    TT   22            C
C      EEEEEEEE XX     XX    TT    IIII NN    NN    TT   22222222      C
C                                                                      C
C -------------------------------------------------------------------- C
C  EXTINT2 DETERMINES THE VALUE OF THE CONTRIBUTION TO AN INTEGRAL     C
C  BASED ON THE TRAPEZOIDAL RULE.                                      C
C**********************************************************************C
C
C     FUNCTION ONLY ALLOWS I BETWEEN 0 AND N
      IF(I.LT.0.OR.I.GT.N) THEN
        WRITE(6, *) 'In EXTINT2: invalid index I. I = ',I
        WRITE(7, *) 'In EXTINT2: invalid index I. I = ',I
        STOP
      ENDIF
C
C     FOR THIS FORMULA TO APPLY, N MUST BE A MULTIPLE OF 2
      IF(MOD(N,2).NE.0.OR.N.LT.2) THEN
        WRITE(6, *) 'In EXTINT2: invalid discretisation N. N = ',N
        WRITE(7, *) 'In EXTINT2: invalid discretisation N. N = ',N
        STOP
      ENDIF
C
      IF(I.EQ.0.OR.I.EQ.N) THEN
        EXTINT2 = 0.5D0*Y
      ELSE
        EXTINT2 = Y
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION EXTINT5(Y,I,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     EEEEEEEE XX     XX TTTTTTTT IIII NN    NN TTTTTTTT 55555555      C
C     EE        XX   XX     TT     II  NNN   NN    TT    55            C
C     EE         XX XX      TT     II  NNNN  NN    TT    55            C
C     EEEEEE      XXX       TT     II  NN NN NN    TT    5555555       C
C     EE         XX XX      TT     II  NN  NNNN    TT          55      C
C     EE        XX   XX     TT     II  NN   NNN    TT    55    55      C
C     EEEEEEEE XX     XX    TT    IIII NN    NN    TT     555555       C
C                                                                      C
C -------------------------------------------------------------------- C
C  EXTINT5 DETERMINES THE VALUE OF THE CONTRIBUTION TO AN INTEGRAL     C
C  BASED ON THE REPEATED 5-POINT NEWTON-COTES FORMULA (ALSO KNOWN AS   C
C  BODE'S RULE), WHERE Y IS THE VALUE F(X_I), I IS THE INDEX OF        C
C  INTEREST, AND THE START AND END POINTS OF THE INTEGRAL ARE F(X_0)   C
C  AND F(X_N) RESPECTIVELY. IT IS NUMERICALLY ADVANTAGEOUS TO MULTIPLY C
C  BY 2.0D0/45.0D0 AFTER THE SUM IS CALCULATED.                        C
C**********************************************************************C
C
C     FUNCTION ONLY ALLOWS I BETWEEN 0 AND N
      IF(I.LT.0.OR.I.GT.N) THEN
        WRITE(6, *) 'In EXTINT5: invalid index I. I = ',I
        WRITE(7, *) 'In EXTINT5: invalid index I. I = ',I
        STOP
      ENDIF
C
C     FOR THIS FORMULA TO APPLY, N MUST BE A MULTIPLE OF 4
      IF(MOD(N,4).NE.0.OR.N.LT.4) THEN
        M = MOD(N,4)
        WRITE(6, *) 'In EXTINT5: invalid discretisation N. N = ',M
        WRITE(7, *) 'In EXTINT5: invalid discretisation N. N = ',M
        STOP
      ENDIF
C
      IF(I.EQ.0.OR.I.EQ.N) THEN
        EXTINT5 = 7.0D0*Y
      ELSEIF(MOD(I,4).EQ.0) THEN
        EXTINT5 = 14.0D0*Y
      ELSEIF(MOD(I,4).EQ.1.OR.MOD(I,4).EQ.3) THEN
        EXTINT5 = 32.0D0*Y
      ELSEIF(MOD(I,4).EQ.2) THEN
        EXTINT5 = 12.0D0*Y
      ENDIF
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
C
C
      FUNCTION GCF(A,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                       GGGGGG   CCCCCC  FFFFFFFF                      C
C                      GG    GG CC    CC FF                            C
C                      GG       CC       FF                            C
C                      GG       CC       FFFFFF                        C
C                      GG   GGG CC       FF                            C
C                      GG    GG CC    CC FF                            C
C                       GGGGGG   CCCCCC  FF                            C
C                                                                      C
C -------------------------------------------------------------------- C
C  GCF EVALUATES THE CONFLUENT HYPERGEOMETRIC FUNCTION U(1,1+A;X) BY   C
C  CONTINUED FRACTION EXPANSION -- (6.2.7) IN NUMERICAL RECIPES, USING C
C  THE MODIFIED LENZ'S METHOD WITH B0 = 0.0D0.                         C
C  REQUIRE X > 0.0D0, AND CONVERGENCE IS RAPID FOR X > A+1.0D0.        C
C                                                                      C
C                   1      1.(1-A)   2.(2-A)                           C
C  U(1,1+A;X) =  ------- . ------- . ------- ...                       C
C                X+1-A-    X+3-A-    X+5-A-                            C
C                                                                      C
C**********************************************************************C
      DATA DLT,EPS/1.0D-30,1.0D-14/
      DATA ITM/100/
C
C     FUNCTION ONLY ALLOWS PARAMETERS X>0
      IF(X.LT.0.0D0) THEN
        WRITE(6, *) 'In GCF: argument X must be positive. X = ',X
        WRITE(7, *) 'In GCF: argument X must be positive. X = ',X
      ENDIF
C
C     FUNCTION CONVERGES ONLY IF X>A+1
      IF(X.LT.A+1.0D0) THEN
        WRITE(6, *) 'In GCF: use a series expansion for this order.',A,X
        WRITE(7, *) 'In GCF: use a series expansion for this order.',A,X
      ENDIF
C
      DEN = X+1.0D0-A
      SAF = 1.0D0/DLT
      RAT = 1.0D0/DEN
      GCF = RAT
C
C     ITERATE TO CONVERGENCE
      DO I=1,ITM
        
        RNM  =-I*(I-A)
        DEN = DEN + 2.0D0
C
C       UPDATE RAT AND FIX IF RESULT UNDERFLOWS
        RAT = RNM*RAT + DEN
        IF(DABS(RAT).LT.DLT) THEN
          RAT = DLT
        ENDIF
C
C       UPDATE SAF AND FIX IF RESULT UNDERFLOWS
        SAF = DEN + RNM/SAF
        IF(DABS(SAF).LT.DLT) THEN
          SAF = DLT
        ENDIF
C
        RAT = 1.0D0/RAT
        DEL = RAT*SAF
C
C       UPDATE GCF AND CHECK AGREEMENT
        GCF = GCF*DEL
        IF(DABS(DEL-1.0D0).LT.EPS) THEN
          GOTO 10
        ENDIF
        
      ENDDO
C
C     METHOD FAILURE
      WRITE(6, *) 'In GCF: A too large or ITM too small. A = ',A
      WRITE(7, *) 'In GCF: A too large or ITM too small. A = ',A
      STOP
C
C     METHOD SUCCESS
10    CONTINUE
C
      RETURN
      END
C
C
      FUNCTION PHIX(L,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                   PPPPPPP  HH    HH IIII XX    XX                    C
C                   PP    PP HH    HH  II   XX  XX                     C
C                   PP    PP HH    HH  II    XXXX                      C
C                   PP    PP HHHHHHHH  II     XX                       C
C                   PPPPPPP  HH    HH  II    XXXX                      C
C                   PP       HH    HH  II   XX  XX                     C
C                   PP       HH    HH IIII XX    XX                    C
C                                                                      C
C -------------------------------------------------------------------- C
C  PHIX IS AN AUXILLIARY FUNCTION WHICH EVALUATES AS A POWER SERIES    C
C  OF TRUNCTATED LENGTH N:0->20 THE COMBINATION OF LOWER INCOMPLETE    C
C  GAMMA FUNCTIONS NEEDED FOR ATOMIC UNIFORM NUCLEUS MATRIX ELEMENTS:  C
C                                                                      C
C  Φ(ℓ,x) = [3*x*γ(ℓ+1/2,x)-γ(ℓ+3/2,x)]/8*π*x^(3/2) with ℓ>0 and x>0.  C
C                                                                      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/FCTS/RFACT(0:80),SFACT(0:80)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     ROUTINE ONLY ALLOWS PARAMETERS L>0
      IF(L.LT.1) THEN
        WRITE(6, *) 'In PHIX: order L must be at least 1. L = ',L
        WRITE(7, *) 'In PHIX: order L must be at least 1. L = ',L
      ENDIF
C
C     PRE-MULTIPLYING FACTOR
      FC1 = X**(L)/PI
C
C     INITIALISE POLYNOMIAL COUNTERS
      PLY = 0.0D0
      PHS = 1.0D0
      XPW = 1.0D0
C
C     LOOP OVER REQUIRED POLYNOMIAL DEGREES (CHOOSE 20)
      DO N=0,20
        RAT = DFLOAT(L+N+2)/DFLOAT((2*L+2*N+1)*(2*L+2*N+3))
        PLY = PLY + PHS*XPW*RAT/RFACT(N)
        XPW = XPW*X
        PHS =-PHS
      ENDDO
C
C     VALUE OF INTEGER-ORDER UPPER INCOMPLETE GAMMA FUNCTION
      PHIX = FC1*PLY
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
      FUNCTION UEHINT0(IZ,L,EIJ,COUPLING)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       UU    UU EEEEEEEE HH    HH IIII NN    NN TTTTTTTT 000000       C
C       UU    UU EE       HH    HH  II  NNN   NN    TT   00   000      C
C       UU    UU EE       HH    HH  II  NNNN  NN    TT   00  0000      C
C       UU    UU EEEEEE   HHHHHHHH  II  NN NN NN    TT   00 00 00      C
C       UU    UU EE       HH    HH  II  NN  NNNN    TT   0000  00      C
C       UU    UU EE       HH    HH  II  NN   NNN    TT   000   00      C
C        UUUUUU  EEEEEEEE HH    HH IIII NN    NN    TT    000000       C
C                                                                      C
C -------------------------------------------------------------------- C
C  UEHINT0 CALCULATES AN ELEMENTARY UEHLING MATRIX ELEMENT OVER A      C
C  PAIR OF ATOMIC BASIS FUNCTIONS, WHERE L IS A NON-NEGATIVE INTEGER.  C
C  THE UEHLING POTENTIAL IS REPRESENTED BY A FITTING SET.              C
C                        ∞                                             C
C       UEHINT0(IZ,L,λ) = ∫ r^2L exp(-λ r^2) r^2*V_ueh(r;IZ) dr        C
C                        0                                             C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*4 COUPLING
C
      COMMON/BUEE/RUEE(MCT,3),FUEE(MCT,MFT),XUEE(MCT,MFT),NUEE(MCT)
      COMMON/BUEU/RUEU(MCT,3),FUEU(MCT,MFT),XUEU(MCT,MFT),NUEU(MCT)
      COMMON/BUEP/RUEP(MCT,3),FUEP(MCT,MFT),XUEP(MCT,MFT),NUEP(MCT)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In UEHINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In UEHINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
      UEHINT0 = 0.0D0
      DO IFT=1,NUEE(IZ)
C
        IF(COUPLING.EQ.'ELEC') THEN
          XI  = XUEE(IZ,IFT)
          FC  = FUEE(IZ,IFT)
        ELSEIF(COUPLING.EQ.'MUON') THEN
          XI  = XUEU(IZ,IFT)
          FC  = FUEU(IZ,IFT)
        ELSEIF(COUPLING.EQ.'PION') THEN
          XI  = XUEP(IZ,IFT)
          FC  = FUEP(IZ,IFT)
        ENDIF
C
        PES = XI+EIJ
        RGN = 1.0D0/(PES*DSQRT(PES))
        XEL = PES**L
C
        UEHINT0 = UEHINT0 + 0.5D0*FC*GAMHLF(2*L+3)*RGN/XEL
C
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION UEHINT0QE(IZ,L,EIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       UU    UU EEEEEEEE HH    HH IIII NN    NN TTTTTTTT 000000       C
C       UU    UU EE       HH    HH  II  NNN   NN    TT   00   000      C
C       UU    UU EE       HH    HH  II  NNNN  NN    TT   00  0000      C
C       UU    UU EEEEEE   HHHHHHHH  II  NN NN NN    TT   00 00 00      C
C       UU    UU EE       HH    HH  II  NN  NNNN    TT   0000  00      C
C       UU    UU EE       HH    HH  II  NN   NNN    TT   000   00      C
C        UUUUUU  EEEEEEEE HH    HH IIII NN    NN    TT    000000       C
C                                                                      C
C -------------------------------------------------------------------- C
C  UEHINT0Q CALCULATES AN ELEMENTARY UEHLING MATRIX ELEMENT OVER A     C
C  PAIR OF ATOMIC BASIS FUNCTIONS, WHERE L IS A NON-NEGATIVE INTEGER.  C
C  THE UEHLING POTENTIAL HAS ALREADY BEEN CALCULATED ON A GRID.        C
C                        ∞                                             C
C        UEHINT0Q(Z,L,λ) = ∫ r^2L exp(-λ r^2) r^2*V_ueh(r;Z) dr        C
C                        0                                             C
C -------------------------------------------------------------------- C
C ▶ IF YOU DON'T HAVE ACCESS TO VVAC, GO BACK AND ENABLE THAT CHUNK OF C
C   CODE IN THE VUEHGEN SUBROUTINE IN BERTHA-NUCLEAR.                  C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/QUEE/RAD(0:NRAD),VVAC(MCT,0:NRAD),RORI,RMID,RMAX,NLIN,NEXP
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In UEHINT0QE: invalid parameter L. L = ',L
        WRITE(7, *) 'In UEHINT0QE: invalid parameter L. L = ',L
        STOP
      ENDIF
C
C     STEP SIZES
      HL = (RMID-RORI)/DFLOAT(NLIN)
      HE = DLOG(RMAX/RMID)/DFLOAT(NEXP)
C
C     UNIFORMLY-SPACED GRID VALUES
      X0M = 0.0D0
      DO N=0,NLIN
        RN  = RAD(N)
        Z1  = RN**(2*L)
        Z2  = DEXP(-EIJ*RN*RN)
        Z3  = VVAC(IZ,N)
        X0M = X0M + EXTINT11(Z1*Z2*Z3,N,NLIN)
      ENDDO
      X0M = 5.0D0*HL*X0M/299376.0D0
C
C     EXPONENTIALLY-SPACED GRID VALUES
      XMI = 0.0D0
      DO N=NLIN,NRAD
        RN  = RAD(N)
        Z1  = RN**(2*L+1)
        Z2  = DEXP(-EIJ*RN*RN)
        Z3  = VVAC(IZ,N)
        XMI = XMI + EXTINT11(Z1*Z2*Z3,N-NLIN,NEXP)
      ENDDO
      XMI = 5.0D0*HE*XMI/299376.0D0
C        
      UEHINT0QE = X0M+XMI
C
      RETURN
      END
C
C
      FUNCTION UEHINT0QU(IZ,L,EIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       UU    UU EEEEEEEE HH    HH IIII NN    NN TTTTTTTT 000000       C
C       UU    UU EE       HH    HH  II  NNN   NN    TT   00   000      C
C       UU    UU EE       HH    HH  II  NNNN  NN    TT   00  0000      C
C       UU    UU EEEEEE   HHHHHHHH  II  NN NN NN    TT   00 00 00      C
C       UU    UU EE       HH    HH  II  NN  NNNN    TT   0000  00      C
C       UU    UU EE       HH    HH  II  NN   NNN    TT   000   00      C
C        UUUUUU  EEEEEEEE HH    HH IIII NN    NN    TT    000000       C
C                                                                      C
C -------------------------------------------------------------------- C
C  UEHINT0Q CALCULATES AN ELEMENTARY UEHLING MATRIX ELEMENT OVER A     C
C  PAIR OF ATOMIC BASIS FUNCTIONS, WHERE L IS A NON-NEGATIVE INTEGER.  C
C  THE UEHLING POTENTIAL HAS ALREADY BEEN CALCULATED ON A GRID.        C
C                        ∞                                             C
C        UEHINT0Q(Z,L,λ) = ∫ r^2L exp(-λ r^2) r^2*V_ueh(r;Z) dr        C
C                        0                                             C
C -------------------------------------------------------------------- C
C ▶ IF YOU DON'T HAVE ACCESS TO VVAC, GO BACK AND ENABLE THAT CHUNK OF C
C   CODE IN THE VUEHGEN SUBROUTINE IN BERTHA-NUCLEAR.                  C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/QUEU/RAD(0:NRAD),VVAC(MCT,0:NRAD),RORI,RMID,RMAX,NLIN,NEXP
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In UEHINT0QU: invalid parameter L. L = ',L
        WRITE(7, *) 'In UEHINT0QU: invalid parameter L. L = ',L
        STOP
      ENDIF
C
C     STEP SIZES
      HL = (RMID-RORI)/DFLOAT(NLIN)
      HE = DLOG(RMAX/RMID)/DFLOAT(NEXP)
C
C     UNIFORMLY-SPACED GRID VALUES
      X0M = 0.0D0
      DO N=0,NLIN
        RN  = RAD(N)
        Z1  = RN**(2*L)
        Z2  = DEXP(-EIJ*RN*RN)
        Z3  = VVAC(IZ,N)
        X0M = X0M + EXTINT11(Z1*Z2*Z3,N,NLIN)
      ENDDO
      X0M = 5.0D0*HL*X0M/299376.0D0
C
C     EXPONENTIALLY-SPACED GRID VALUES
      XMI = 0.0D0
      DO N=NLIN,NRAD
        RN  = RAD(N)
        Z1  = RN**(2*L+1)
        Z2  = DEXP(-EIJ*RN*RN)
        Z3  = VVAC(IZ,N)
        XMI = XMI + EXTINT11(Z1*Z2*Z3,N-NLIN,NEXP)
      ENDDO
      XMI = 5.0D0*HE*XMI/299376.0D0
C        
      UEHINT0QU = X0M+XMI
C
      RETURN
      END
C
C
      FUNCTION UEHINT0QP(IZ,L,EIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       UU    UU EEEEEEEE HH    HH IIII NN    NN TTTTTTTT 000000       C
C       UU    UU EE       HH    HH  II  NNN   NN    TT   00   000      C
C       UU    UU EE       HH    HH  II  NNNN  NN    TT   00  0000      C
C       UU    UU EEEEEE   HHHHHHHH  II  NN NN NN    TT   00 00 00      C
C       UU    UU EE       HH    HH  II  NN  NNNN    TT   0000  00      C
C       UU    UU EE       HH    HH  II  NN   NNN    TT   000   00      C
C        UUUUUU  EEEEEEEE HH    HH IIII NN    NN    TT    000000       C
C                                                                      C
C -------------------------------------------------------------------- C
C  UEHINT0Q CALCULATES AN ELEMENTARY UEHLING MATRIX ELEMENT OVER A     C
C  PAIR OF ATOMIC BASIS FUNCTIONS, WHERE L IS A NON-NEGATIVE INTEGER.  C
C  THE UEHLING POTENTIAL HAS ALREADY BEEN CALCULATED ON A GRID.        C
C                        ∞                                             C
C        UEHINT0Q(Z,L,λ) = ∫ r^2L exp(-λ r^2) r^2*V_ueh(r;Z) dr        C
C                        0                                             C
C -------------------------------------------------------------------- C
C ▶ IF YOU DON'T HAVE ACCESS TO VVAC, GO BACK AND ENABLE THAT CHUNK OF C
C   CODE IN THE VUEHGEN SUBROUTINE IN BERTHA-NUCLEAR.                  C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/QUEP/RAD(0:NRAD),VVAC(MCT,0:NRAD),RORI,RMID,RMAX,NLIN,NEXP
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In UEHINT0QP: invalid parameter L. L = ',L
        WRITE(7, *) 'In UEHINT0QP: invalid parameter L. L = ',L
        STOP
      ENDIF
C
C     STEP SIZES
      HL = (RMID-RORI)/DFLOAT(NLIN)
      HE = DLOG(RMAX/RMID)/DFLOAT(NEXP)
C
C     UNIFORMLY-SPACED GRID VALUES
      X0M = 0.0D0
      DO N=0,NLIN
        RN  = RAD(N)
        Z1  = RN**(2*L)
        Z2  = DEXP(-EIJ*RN*RN)
        Z3  = VVAC(IZ,N)
        X0M = X0M + EXTINT11(Z1*Z2*Z3,N,NLIN)
      ENDDO
      X0M = 5.0D0*HL*X0M/299376.0D0
C
C     EXPONENTIALLY-SPACED GRID VALUES
      XMI = 0.0D0
      DO N=NLIN,NRAD
        RN  = RAD(N)
        Z1  = RN**(2*L+1)
        Z2  = DEXP(-EIJ*RN*RN)
        Z3  = VVAC(IZ,N)
        XMI = XMI + EXTINT11(Z1*Z2*Z3,N-NLIN,NEXP)
      ENDDO
      XMI = 5.0D0*HE*XMI/299376.0D0
C        
      UEHINT0QP = X0M+XMI
C
      RETURN
      END
C
C
      FUNCTION WKRINT0(IZ,L,EIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    WW        WW KK    KK RRRRRRR  IIII NN    NN TTTTTTTT 000000      C
C    WW        WW KK   KK  RR    RR  II  NNN   NN    TT   00   000     C
C    WW   WW   WW KK  KK   RR    RR  II  NNNN  NN    TT   00  0000     C
C    WW  WWWW  WW KKKKK    RR    RR  II  NN NN NN    TT   00 00 00     C
C    WW WW  WW WW KK  KK   RRRRRRR   II  NN  NNNN    TT   0000  00     C
C    WWWW    WWWW KK   KK  RR    RR  II  NN   NNN    TT   000   00     C
C    WW        WW KK    KK RR    RR IIII NN    NN    TT    000000      C
C                                                                      C
C -------------------------------------------------------------------- C
C  WKRINT0 CALCULATES AN ELEMENTARY WICHMANN-KROLL MATRIX ELEMENT OVER C
C  A PAIR OF ATOMIC BASIS FUNCTIONS, WHERE L IS A NON-NEGATIVE INTEGER.C
C  THE WICHMANN-KROLL POTENTIAL IS REPRESENTED BY A FITTING SET.       C
C                        ∞                                             C
C       WKRINT0(IZ,L,λ) = ∫ r^2L exp(-λ r^2) r^2*V_wkr(r;IZ) dr        C
C                        0                                             C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/BWKR/RWKR(MCT,3),FWKR(MCT,MFT),XWKR(MCT,MFT),NWKR(MCT)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In WKRINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In WKRINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
      WKRINT0 = 0.0D0
      DO IFT=1,NWKR(IZ)
C
        XI  = XWKR(IZ,IFT)
        FC  = FWKR(IZ,IFT)
C
        PES = XI+EIJ
        RGN = 1.0D0/(PES*DSQRT(PES))
        XEL = PES**L
C
        WKRINT0 = WKRINT0 + 0.5D0*FC*GAMHLF(2*L+3)*RGN/XEL
C
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION WKRINT0Q(IZ,L,EIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    WW        WW KK    KK RRRRRRR  IIII NN    NN TTTTTTTT 000000      C
C    WW        WW KK   KK  RR    RR  II  NNN   NN    TT   00   000     C
C    WW   WW   WW KK  KK   RR    RR  II  NNNN  NN    TT   00  0000     C
C    WW  WWWW  WW KKKKK    RR    RR  II  NN NN NN    TT   00 00 00     C
C    WW WW  WW WW KK  KK   RRRRRRR   II  NN  NNNN    TT   0000  00     C
C    WWWW    WWWW KK   KK  RR    RR  II  NN   NNN    TT   000   00     C
C    WW        WW KK    KK RR    RR IIII NN    NN    TT    000000      C
C                                                                      C
C -------------------------------------------------------------------- C
C  WKRINT0 CALCULATES AN ELEMENTARY WICHMANN-KROLL MATRIX ELEMENT OVER C
C  A PAIR OF ATOMIC BASIS FUNCTIONS, WHERE L IS A NON-NEGATIVE INTEGER.C
C  THE WICHMANN-KROLL POTENTIAL HAS ALREADY BEEN CALCULATED ON A GRID. C
C                        ∞                                             C
C        WKRINT0(Z,L,λ) = ∫ r^2L exp(-λ r^2) r^2*V_wkr(r;Z) dr         C
C                        0                                             C
C -------------------------------------------------------------------- C
C ▶ IF YOU DON'T HAVE ACCESS TO VVAC, GO BACK AND ENABLE THAT CHUNK OF C
C   CODE IN THE VWKRGEN SUBROUTINE IN BERTHA-NUCLEAR.                  C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/QWKR/RAD(0:NRAD),VVAC(MCT,0:NRAD),RORI,RMID,RMAX,NLIN,NEXP
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In WKRINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In WKRINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
C     STEP SIZES
      HL = (RMID-RORI)/DFLOAT(NLIN)
      HE = DLOG(RMAX/RMID)/DFLOAT(NEXP)
C
C     UNIFORMLY-SPACED GRID VALUES
      X0M = 0.0D0
      DO N=0,NLIN
        RN  = RAD(N)
        Z1  = RN**(2*L)
        Z2  = DEXP(-EIJ*RN*RN)
        Z3  = VVAC(IZ,N)
        X0M = X0M + EXTINT11(Z1*Z2*Z3,N,NLIN)
      ENDDO
      X0M = 5.0D0*HL*X0M/299376.0D0
C
C     EXPONENTIALLY-SPACED GRID VALUES
      XMI = 0.0D0
      DO N=NLIN,NRAD
        RN  = RAD(N)
        Z1  = RN**(2*L+1)
        Z2  = DEXP(-EIJ*RN*RN)
        Z3  = VVAC(IZ,N)
        XMI = XMI + EXTINT11(Z1*Z2*Z3,N-NLIN,NEXP)
      ENDDO
      XMI = 5.0D0*HE*XMI/299376.0D0
C        
      WKRINT0Q = X0M+XMI
C
      RETURN
      END
C
C
      FUNCTION RKSBINT0(IZ,L,EIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       KK    KK SSSSSS  BBBBBBB  IIII NN    NN TTTTTTTT 000000        C
C       KK   KK SS    SS BB    BB  II  NNN   NN    TT   00   000       C
C       KK  KK  SS       BB    BB  II  NNNN  NN    TT   00  0000       C
C       KKKKK    SSSSSS  BBBBBBB   II  NN NN NN    TT   00 00 00       C
C       KK  KK        SS BB    BB  II  NN  NNNN    TT   0000  00       C
C       KK   KK SS    SS BB    BB  II  NN   NNN    TT   000   00       C
C       KK    KK SSSSSS  BBBBBBB  IIII NN    NN    TT    000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  KSBINT0 CALCULATES AN ELEMENTARY KÄLLÉN-SABRY MATRIX ELEMENT OVER A C
C  PAIR OF ATOMIC BASIS FUNCTIONS, WHERE L IS A NON-NEGATIVE INTEGER.  C
C  THE KÄLLÉN-SABRY POTENTIAL IS REPRESENTED BY A FITTING SET.         C
C                        ∞                                             C
C       KSBINT0(IZ,L,λ) = ∫ r^2L exp(-λ r^2) r^2*V_ksb(r;IZ) dr        C
C                        0                                             C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/BKSB/RKSB(MCT,3),FKSB(MCT,MFT),XKSB(MCT,MFT),NKSB(MCT)
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In KSBINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In KSBINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
      RKSBINT0 = 0.0D0
      DO IFT=1,NKSB(IZ)
C
        XI  = XKSB(IZ,IFT)
        FC  = FKSB(IZ,IFT)
C
        PES = XI+EIJ
        RGN = 1.0D0/(PES*DSQRT(PES))
        XEL = PES**L
C
        RKSBINT0 = RKSBINT0 + 0.5D0*FC*GAMHLF(2*L+3)*RGN/XEL
C
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION RKSBINT0Q(IZ,L,EIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       KK    KK SSSSSS  BBBBBBB  IIII NN    NN TTTTTTTT 000000        C
C       KK   KK SS    SS BB    BB  II  NNN   NN    TT   00   000       C
C       KK  KK  SS       BB    BB  II  NNNN  NN    TT   00  0000       C
C       KKKKK    SSSSSS  BBBBBBB   II  NN NN NN    TT   00 00 00       C
C       KK  KK        SS BB    BB  II  NN  NNNN    TT   0000  00       C
C       KK   KK SS    SS BB    BB  II  NN   NNN    TT   000   00       C
C       KK    KK SSSSSS  BBBBBBB  IIII NN    NN    TT    000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  KSBINT0 CALCULATES AN ELEMENTARY KÄLLÉN-SABRY MATRIX ELEMENT OVER A C
C  PAIR OF ATOMIC BASIS FUNCTIONS, WHERE L IS A NON-NEGATIVE INTEGER.  C
C  THE KÄLLÉN-SABRY POTENTIAL HAS ALREADY BEEN CALCULATED ON A GRID.   C
C                        ∞                                             C
C        KSBINT0(Z,L,λ) = ∫ r^2L exp(-λ r^2) r^2*V_ksb(r;Z) dr         C
C                        0                                             C
C -------------------------------------------------------------------- C
C ▶ IF YOU DON'T HAVE ACCESS TO VVAC, GO BACK AND ENABLE THAT CHUNK OF C
C   CODE IN THE VKSBGEN SUBROUTINE IN BERTHA-NUCLEAR.                  C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/QKSB/RAD(0:NRAD),VVAC(MCT,0:NRAD),RORI,RMID,RMAX,NLIN,NEXP
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In KSBINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In KSBINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
C     STEP SIZES
      HL = (RMID-RORI)/DFLOAT(NLIN)
      HE = DLOG(RMAX/RMID)/DFLOAT(NEXP)
C
C     UNIFORMLY-SPACED GRID VALUES
      X0M = 0.0D0
      DO N=0,NLIN
        RN  = RAD(N)
        Z1  = RN**(2*L)
        Z2  = DEXP(-EIJ*RN*RN)
        Z3  = VVAC(IZ,N)
        X0M = X0M + EXTINT11(Z1*Z2*Z3,N,NLIN)
      ENDDO
      X0M = 5.0D0*HL*X0M/299376.0D0
C
C     EXPONENTIALLY-SPACED GRID VALUES
      XMI = 0.0D0
      DO N=NLIN,NRAD
        RN  = RAD(N)
        Z1  = RN**(2*L+1)
        Z2  = DEXP(-EIJ*RN*RN)
        Z3  = VVAC(IZ,N)
        XMI = XMI + EXTINT11(Z1*Z2*Z3,N-NLIN,NEXP)
      ENDDO
      XMI = 5.0D0*HE*XMI/299376.0D0
C        
      RKSBINT0Q = X0M+XMI
C
      RETURN
      END
C
C
      SUBROUTINE ANGCLM0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     AA    NN    NN  GGGGGG   CCCCCC  LL       MM       MM  000000    C
C    AAAA   NNN   NN GG    GG CC    CC LL       MMM     MMM 00   000   C
C   AA  AA  NNNN  NN GG       CC       LL       MMMM   MMMM 00  0000   C
C  AA    AA NN NN NN GG       CC       LL       MM MM MM MM 00 00 00   C
C  AAAAAAAA NN  NNNN GG   GGG CC       LL       MM  MMM  MM 0000  00   C
C  AA    AA NN   NNN GG    GG CC    CC LL       MM   M   MM 000   00   C
C  AA    AA NN    NN  GGGGGG   CCCCCC  LLLLLLLL MM       MM  000000    C
C                                                                      C
C -------------------------------------------------------------------- C
C  ANGCLM0 EVALUATES THE ANGULAR COEFFICIENTS OF THE COULOMB           C
C  INTERACTIONS FOR CLOSED SHELLS IN THE (L1,L2) MANIFOLD.             C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      COMMON/B0QN/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
      COMMON/XTNS/BK(MNU,4),ELL(MNU,4),ESS(MNU,4),ESL(MNU,4),GSL(MNU,4)
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
C     START AND END PARAMETERS FROM TRIANGLE RULE
      NUI = IABS(LQNA-LQNB)
      NUF = LQNA+LQNB+1
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      LTEN = 0
      DO NU=NUI,NUF
C
C       PARITY OF 'LQNA+LQNB+NU' (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQNA+LQNB+NU,2)
C
        IF(IPARAB.EQ.0) THEN
C       ONLY ANGULAR COEFFICIENTS OF EVEN PARITY ARE NON-ZERO
C
C         SAVE THIS TENSOR ORDER
          LTEN      = LTEN+1
          NUS(LTEN) = NU
C
          IF(HMLT.EQ.'NORL') THEN
            BK(LTEN,1) = 0.5D0*ANGSQLS(LQNA,LQNB,NU)
          ELSE
            BK(LTEN,1) = ANGSQJJ(JLA,JLB,NU)
            BK(LTEN,2) = ANGSQJJ(JLA,JRB,NU)
            BK(LTEN,3) = ANGSQJJ(JRA,JLB,NU)
            BK(LTEN,4) = ANGSQJJ(JRA,JRB,NU)
          ENDIF
        ENDIF
C
      ENDDO
C
C     NUMBER OF SURVIVING TENSOR ORDERS
      NUNUM = LTEN
C
      RETURN
      END
C
C
      SUBROUTINE ANGBRT0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        AA    NN    NN  GGGGGG  BBBBBBB  RRRRRRR TTTTTTTT 000000      C
C       AAAA   NNN   NN GG    GG BB    BB RR    RR   TT   00   000     C
C      AA  AA  NNNN  NN GG       BB    BB RR    RR   TT   00  0000     C
C     AA    AA NN NN NN GG       BBBBBBB  RR    RR   TT   00 00 00     C
C     AAAAAAAA NN  NNNN GG   GGG BB    BB RRRRRRR    TT   0000  00     C
C     AA    AA NN   NNN GG    GG BB    BB RR    RR   TT   000   00     C
C     AA    AA NN    NN  GGGGGG  BBBBBBB  RR    RR   TT    000000      C
C                                                                      C
C -------------------------------------------------------------------- C
C  ANGBRT0 EVALUATES ANGULAR COEFFICIENTS OF THE ATOMIC CLOSED SHELL   C
C  BREIT INTERACTION FOR ALL (K1,K2) VALUES IN THE MANIFOLD (L1,L2).   C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C  ▶ NUNUM - NUMBER OF NU VALUES THAT SATISFY PARITY RESTRICTION RULE. C
C  ▶ NUI - MINIMUM NU VALUE IN THIS MANIFOLD.                          C
C  ▶ NUF - MAXIMUM NU VALUE IN THIS MANIFOLD.                          C
C  ▶ ELL(MNU,4) - ANGULAR TERMS FOR CLOSED BREIT EK(JA,JB;LL) TERMS    C
C  ▶ ESS(MNU,4) - ANGULAR TERMS FOR CLOSED BREIT EK(JA,JB;SS) TERMS    C
C  ▶ ESL(MNU,4) - ANGULAR TERMS FOR CLOSED BREIT EK(JA,JB;SL) TERMS    C
C  ▶ GSL(MNU,4) - ANGULAR TERMS FOR CLOSED BREIT EK(JA,JB;SL) TERMS    C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION SCOEFF(4,2)
C
      COMMON/B0QN/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
      COMMON/XTNS/BK(MNU,4),ELL(MNU,4),ESS(MNU,4),ESL(MNU,4),GSL(MNU,4)
C
C     INITIALISE COEFFICIENT ARRAYS
      DO LTEN=1,MNU
        DO N=1,4
          ELL(LTEN,N) = 0.0D0
          ESS(LTEN,N) = 0.0D0
          ESL(LTEN,N) = 0.0D0
          GSL(LTEN,N) = 0.0D0
        ENDDO
        NUS(LTEN) = 0
      ENDDO
      NUNUM = 0
C
C     SPECIFY ALLOWED NU VALUES (PARITY LQNA+LQNB+NU MUST BE ODD)
      NUNUM = 0
      DO NU=0,IABS(LQNA+LQNB+1)
        IF(NU.GE.0.AND.MOD(LQNA+LQNB+NU,2).EQ.1) THEN
          NUNUM = NUNUM+1
          NUS(NUNUM) = NU
        ENDIF
      ENDDO
C
C     BREIT TENSOR ORDER LIMITS BASED ON PARITY CHECK
      NUI = NUS(1)
      NUF = NUS(NUNUM)
C
C**********************************************************************C
C     (1) KQNA < 0 AND KQNB < 0   (CANNOT SKIP)                        C
C**********************************************************************C
C
C     KQNA AND KQNB
      KQNA =-LQNA-1
      KQNB =-LQNB-1
C
C     JQNA AND JQNB
      JQNA = 2*IABS(KQNA)-1
      JQNB = 2*IABS(KQNB)-1
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      LTEN = 1
      DO NU=NUI,NUF
C
C       SKIP IF THIS VALUE OF NU EXCEEDS TRIANGLE RULE FOR KQN BLOCK
        IF(NU.GT.(JQNA+JQNB)/2) GOTO 101
C
C       ANGULAR FACTOR FOR THIS (KQNA,KQNB,NU)
        AJJN0 = ANGSQJJ(JQNA,JQNB,NU)
C
C       PARITY OF 'LQNA+LQNB+NU' (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQNA+LQNB+NU,2)
C
C       CASE 1: ODD-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.1.AND.NU.NE.0) THEN
C
C         INTERMEDIATE (KQN,NU) COEFFICIENTS
          RCOEFF = DFLOAT((KQNA+KQNB)*(KQNA+KQNB))/DFLOAT(NU*(NU+1))
C
C         THIS CONTRIBUTION BELONGS IN "NU" BIN
          NUS(LTEN) = NU
C
C         CONTRIBUTIONS TO TENSOR ORDER "NU"
          ELL(LTEN,1) = ELL(LTEN,1) + AJJN0*RCOEFF
          ESS(LTEN,1) = ESS(LTEN,1) + AJJN0*RCOEFF
          ESL(LTEN,1) = ESL(LTEN,1) + AJJN0*RCOEFF
C
        ENDIF
C
C       CASE 2: EVEN-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.0) THEN
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          CALL BRCOEF0(SCOEFF,KQNA,KQNB,NU)
C
C         FIRST CONTRIBUTION BELONGS IN "NU-1" BIN
          NUS(LTEN) = NU-1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,1) = ELL(LTEN,1) + AJJN0*SCOEFF(1,1)
          ESL(LTEN,1) = ESL(LTEN,1) + AJJN0*SCOEFF(2,1)
          ESS(LTEN,1) = ESS(LTEN,1) + AJJN0*SCOEFF(3,1)
          GSL(LTEN,1) = GSL(LTEN,1) + AJJN0*SCOEFF(4,1)
C
C         INCREASE TENSOR LIST LENGTH IF ENTRIES ARE NON-ZERO
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SECOND CONTRIBUTION BELONGS IN "NU+1" BIN
          NUS(LTEN) = NU+1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,1) = ELL(LTEN,1) + AJJN0*SCOEFF(1,2)
          ESL(LTEN,1) = ESL(LTEN,1) + AJJN0*SCOEFF(2,2)
          ESS(LTEN,1) = ESS(LTEN,1) + AJJN0*SCOEFF(3,2)
          GSL(LTEN,1) = GSL(LTEN,1) + AJJN0*SCOEFF(4,2)
C
        ENDIF
C
101     CONTINUE
C
      ENDDO
C
C**********************************************************************C
C     (2) KQNA < 0 AND KQNB > 0   (SKIP IF POSSIBLE)                   C
C**********************************************************************C
C
      IF(LQNB.EQ.0) GOTO 200
C
C     KQNA AND KQNB
      KQNA =-LQNA-1
      KQNB = LQNB
C
C     JQNA AND JQNB
      JQNA = 2*IABS(KQNA)-1
      JQNB = 2*IABS(KQNB)-1
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      LTEN = 1
      DO NU=NUI,NUF
C
C       SKIP IF THIS VALUE OF NU EXCEEDS TRIANGLE RULE FOR KQN BLOCK
        IF(NU.GT.(JQNA+JQNB)/2) GOTO 201
C
C       ANGULAR FACTOR FOR THIS (KQNA,KQNB,NU)
        AJJN0 = ANGSQJJ(JQNA,JQNB,NU)
C
C       PARITY OF 'LQNA+LQNB+NU' (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQNA+LQNB+NU,2)
C
C       CASE 1: ODD-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.1.AND.NU.NE.0) THEN
C
C         INTERMEDIATE (KQN,NU) COEFFICIENTS
          RCOEFF = DFLOAT((KQNA+KQNB)*(KQNA+KQNB))/DFLOAT(NU*(NU+1))
C
C         THIS CONTRIBUTION BELONGS IN "NU" BIN
          NUS(LTEN) = NU
C
C         CONTRIBUTIONS TO TENSOR ORDER "NU"
          ELL(LTEN,2) = ELL(LTEN,2) + AJJN0*RCOEFF
          ESS(LTEN,2) = ESS(LTEN,2) + AJJN0*RCOEFF
          ESL(LTEN,2) = ESL(LTEN,2) + AJJN0*RCOEFF
C
        ENDIF
C
C       CASE 2: EVEN-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.0) THEN
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          CALL BRCOEF0(SCOEFF,KQNA,KQNB,NU)
C
C         FIRST CONTRIBUTION BELONGS IN "NU-1" BIN
          NUS(LTEN) = NU-1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,2) = ELL(LTEN,2) + AJJN0*SCOEFF(1,1)
          ESL(LTEN,2) = ESL(LTEN,2) + AJJN0*SCOEFF(2,1)
          ESS(LTEN,2) = ESS(LTEN,2) + AJJN0*SCOEFF(3,1)
          GSL(LTEN,2) = GSL(LTEN,2) + AJJN0*SCOEFF(4,1)
C
C         INCREASE TENSOR LIST LENGTH IF ENTRIES ARE NON-ZERO
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SECOND CONTRIBUTION BELONGS IN "NU+1" BIN
          NUS(LTEN) = NU+1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,2) = ELL(LTEN,2) + AJJN0*SCOEFF(1,2)
          ESL(LTEN,2) = ESL(LTEN,2) + AJJN0*SCOEFF(2,2)
          ESS(LTEN,2) = ESS(LTEN,2) + AJJN0*SCOEFF(3,2)
          GSL(LTEN,2) = GSL(LTEN,2) + AJJN0*SCOEFF(4,2)
C
        ENDIF
C
201     CONTINUE
C
      ENDDO
C
200   CONTINUE
C
C**********************************************************************C
C     (3) KQNA > 0 AND KQNB < 0   (SKIP IF POSSIBLE)                   C
C**********************************************************************C
C
      IF(LQNA.EQ.0) GOTO 300
C
C     KQNA AND KQNB
      KQNA = LQNA
      KQNB =-LQNB-1
C
C     JQNA AND JQNB
      JQNA = 2*IABS(KQNA)-1
      JQNB = 2*IABS(KQNB)-1
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      LTEN = 1
      DO NU=NUI,NUF
C
C       SKIP IF THIS VALUE OF NU EXCEEDS TRIANGLE RULE FOR KQN BLOCK
        IF(NU.GT.(JQNA+JQNB)/2) GOTO 301
C
C       ANGULAR FACTOR FOR THIS (KQNA,KQNB,NU)
        AJJN0 = ANGSQJJ(JQNA,JQNB,NU)
C
C       PARITY OF 'LQNA+LQNB+NU' (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQNA+LQNB+NU,2)
C
C       CASE 1: ODD-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.1.AND.NU.NE.0) THEN
C
C         INTERMEDIATE (KQN,NU) COEFFICIENTS
          RCOEFF = DFLOAT((KQNA+KQNB)*(KQNA+KQNB))/DFLOAT(NU*(NU+1))
C
C         THIS CONTRIBUTION BELONGS IN "NU" BIN
          NUS(LTEN) = NU
C
C         CONTRIBUTIONS TO TENSOR ORDER "NU"
          ELL(LTEN,3) = ELL(LTEN,3) + AJJN0*RCOEFF
          ESS(LTEN,3) = ESS(LTEN,3) + AJJN0*RCOEFF
          ESL(LTEN,3) = ESL(LTEN,3) + AJJN0*RCOEFF
C
        ENDIF
C
C       CASE 2: EVEN-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.0) THEN
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          CALL BRCOEF0(SCOEFF,KQNA,KQNB,NU)
C
C         FIRST CONTRIBUTION BELONGS IN "NU-1" BIN
          NUS(LTEN) = NU-1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,3) = ELL(LTEN,3) + AJJN0*SCOEFF(1,1)
          ESL(LTEN,3) = ESL(LTEN,3) + AJJN0*SCOEFF(2,1)
          ESS(LTEN,3) = ESS(LTEN,3) + AJJN0*SCOEFF(3,1)
          GSL(LTEN,3) = GSL(LTEN,3) + AJJN0*SCOEFF(4,1)
C
C         INCREASE TENSOR LIST LENGTH IF ENTRIES ARE NON-ZERO
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SECOND CONTRIBUTION BELONGS IN "NU+1" BIN
          NUS(LTEN) = NU+1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,3) = ELL(LTEN,3) + AJJN0*SCOEFF(1,2)
          ESL(LTEN,3) = ESL(LTEN,3) + AJJN0*SCOEFF(2,2)
          ESS(LTEN,3) = ESS(LTEN,3) + AJJN0*SCOEFF(3,2)
          GSL(LTEN,3) = GSL(LTEN,3) + AJJN0*SCOEFF(4,2)
C
        ENDIF
C
301     CONTINUE
C
      ENDDO
C
300   CONTINUE
C
C**********************************************************************C
C     (4) KQNA > 0 AND KQNB > 0   (SKIP IF POSSIBLE)                   C
C**********************************************************************C
C
      IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 400
C
C     KQNA AND KQNB
      KQNA = LQNA
      KQNB = LQNB
C
C     JQNA AND JQNB
      JQNA = 2*IABS(KQNA)-1
      JQNB = 2*IABS(KQNB)-1
C
C     LOOP OVER ALL NU VALUES WITHIN TRIANGLE RULE
      LTEN = 1
      DO NU=NUI,NUF
C
C       SKIP IF THIS VALUE OF NU EXCEEDS TRIANGLE RULE FOR KQN BLOCK
        IF(NU.GT.(JQNA+JQNB)/2) GOTO 401
C
C       ANGULAR FACTOR FOR THIS (KQNA,KQNB,NU)
        AJJN0 = ANGSQJJ(JQNA,JQNB,NU)
C
C       PARITY OF 'LQNA+LQNB+NU' (0 IF EVEN, 1 IF ODD)
        IPARAB = MOD(LQNA+LQNB+NU,2)
C
C       CASE 1: ODD-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.1.AND.NU.NE.0) THEN
C
C         INTERMEDIATE (KQN,NU) COEFFICIENTS
          RCOEFF = DFLOAT((KQNA+KQNB)*(KQNA+KQNB))/DFLOAT(NU*(NU+1))
C
C         THIS CONTRIBUTION BELONGS IN "NU" BIN
          NUS(LTEN) = NU
C
C         CONTRIBUTIONS TO TENSOR ORDER "NU"
          ELL(LTEN,4) = ELL(LTEN,4) + AJJN0*RCOEFF
          ESS(LTEN,4) = ESS(LTEN,4) + AJJN0*RCOEFF
          ESL(LTEN,4) = ESL(LTEN,4) + AJJN0*RCOEFF
C
        ENDIF
C
C       CASE 2: EVEN-PARITY CONTRIBUTIONS
        IF(IPARAB.EQ.0) THEN
C
C         CALCULATE INTERMEDIATE COEFFICIENTS
          CALL BRCOEF0(SCOEFF,KQNA,KQNB,NU)
C
C         FIRST CONTRIBUTION BELONGS IN "NU-1" BIN
          NUS(LTEN) = NU-1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,4) = ELL(LTEN,4) + AJJN0*SCOEFF(1,1)
          ESL(LTEN,4) = ESL(LTEN,4) + AJJN0*SCOEFF(2,1)
          ESS(LTEN,4) = ESS(LTEN,4) + AJJN0*SCOEFF(3,1)
          GSL(LTEN,4) = GSL(LTEN,4) + AJJN0*SCOEFF(4,1)
C
C         INCREASE TENSOR LIST LENGTH IF ENTRIES ARE NON-ZERO
          IF(NU.GT.0) THEN
            LTEN = LTEN+1
          ENDIF
C
C         SECOND CONTRIBUTION BELONGS IN "NU+1" BIN
          NUS(LTEN) = NU+1
C
C         CONTRIBUTIONS TO ANGULAR TERMS
          ELL(LTEN,4) = ELL(LTEN,4) + AJJN0*SCOEFF(1,2)
          ESL(LTEN,4) = ESL(LTEN,4) + AJJN0*SCOEFF(2,2)
          ESS(LTEN,4) = ESS(LTEN,4) + AJJN0*SCOEFF(3,2)
          GSL(LTEN,4) = GSL(LTEN,4) + AJJN0*SCOEFF(4,2)
C
        ENDIF
C
401     CONTINUE
C
      ENDDO
C
400   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE BRCOEF0(SCOEFF,KQNA,KQNB,NU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    BBBBBBB  RRRRRRR   CCCCCC   OOOOOO  EEEEEEEE FFFFFFFF 000000      C
C    BB    BB RR    RR CC    CC OO    OO EE       FF      00   000     C
C    BB    BB RR    RR CC       OO    OO EE       FF      00  0000     C
C    BBBBBBB  RR    RR CC       OO    OO EEEEEE   FFFFFF  00 00 00     C
C    BB    BB RRRRRRR  CC       OO    OO EE       FF      0000  00     C
C    BB    BB RR    RR CC    CC OO    OO EE       FF      000   00     C
C    BBBBBBB  RR    RR  CCCCCC   OOOOOO  EEEEEEEE FF       000000      C
C                                                                      C
C -------------------------------------------------------------------- C
C  BRCOEF0 EVALUATES THE INTERMEDIATE COEFFICIENTS OF THE BREIT        C
C  INTERACTION FOR CLOSED SHELLS (TABLE 3 OF GRANT AND PYPER 1976).    C
C**********************************************************************C
      DIMENSION SCOEFF(4,2)
C

      RU  = DFLOAT(NU)
      RK  = DFLOAT(KQNB-KQNA)
C
      IF(NU.GT.0) THEN
        RM = RU-1.0D0
        B1 = (RM+2.0D0)/(2.0D0*   (2.0D0*RM+1.0D0))
        C1 =-(RM-1.0D0)/(2.0D0*RU*(2.0D0*RM+1.0D0))
        SCOEFF(1,1) = (RK+RU)*(C1*RK+B1)
        SCOEFF(2,1) = -B1*RU + C1*RK*RK
        SCOEFF(3,1) = (RK-RU)*(C1*RK-B1)
        SCOEFF(4,1) =  RK    *(C1*RU-B1)
      ELSE
        SCOEFF(1,1) = 0.0D0
        SCOEFF(2,1) = 0.0D0
        SCOEFF(3,1) = 0.0D0
        SCOEFF(4,1) = 0.0D0
      ENDIF
      IF(NU+1.GT.1) THEN
        RP = RU+1.0D0
        B2 = (RP-1.0D0)/(2.0D0*   (2.0D0*RP+1.0D0))
        C2 = (RP+2.0D0)/(2.0D0*RP*(2.0D0*RP+1.0D0))
        SCOEFF(1,2) = (RK-RP)*(C2*RK+B2)
        SCOEFF(2,2) =  B2*RP + C2*RK*RK 
        SCOEFF(3,2) = (RK+RP)*(C2*RK-B2)
        SCOEFF(4,2) =  RK    *(C2*RP+B2)
      ELSE
        SCOEFF(1,2) = 0.0D0
        SCOEFF(2,2) = 0.0D0
        SCOEFF(3,2) = 0.0D0
        SCOEFF(4,2) = 0.0D0
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION ANGSQLS(L1,L2,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      AA    NN    NN  GGGGGG   SSSSSS    QQQQQQ    LL       SSSSSS    C
C     AAAA   NNN   NN GG    GG SS    SS  QQ    QQ   LL      SS    SS   C
C    AA  AA  NNNN  NN GG       SS       QQ      QQ  LL      SS         C
C   AA    AA NN NN NN GG        SSSSSS  QQ      QQ  LL       SSSSSS    C
C   AAAAAAAA NN  NNNN GG   GGG       SS QQ      QQ  LL            SS   C
C   AA    AA NN   NNN GG    GG SS    SS  QQ    QQ   LL      SS    SS   C
C   AA    AA NN    NN  GGGGGG   SSSSSS    QQQQQQ QQ LLLLLLLL SSSSSS    C
C                                                                      C
C -------------------------------------------------------------------- C
C  ANGSQLS EVALUATES THE NON-RELATIVISTIC 3-J SYMBOL FOR ATOMIC        C
C  COULOMB ANGULAR COEFFICIENT ROUTINES, TAKEN FROM BRINK AND SATCHLER.C
C  L1,L2, AND K MUST BE EQUAL TO THE ACTUAL (INTEGER) ANGULAR MOMENTA  C
C  OF THE ELECTRON AND PHOTON.                                         C
C**********************************************************************C
      COMMON/FCTS/RFACT(0:80),SFACT(0:80)
C
C     TRIANGLE INEQUALITY RESTRICTIONS
      IF(K.LT.IABS(L1-L2).OR.K.GT.(L1+L2)) THEN
        ANGSQLS = 0.0D0
        RETURN
      ENDIF
      LLK = L1+L2+K
C
C     PARITY SELECTION RULE
      IF((LLK/2)*2.NE.LLK) THEN
        ANGSQLS = 0.0D0
        RETURN
      ENDIF
C
      RF1 = RFACT(  L1+L2-K   )
      RF2 = RFACT(- L1+L2+K   )
      RF3 = RFACT(  L1-L2+K   )
      RF4 = RFACT(  L1+L2+K +1)
      RF5 = RFACT(( L1+L2+K)/2)
      RF6 = RFACT(( L1+L2-K)/2)
      RF7 = RFACT(( L1-L2+K)/2)
      RF8 = RFACT((-L1+L2+K)/2)
C
      T1 = RF1*RF2*RF3
      T2 = T1/RF4
      T3 = RF6*RF7*RF8
      T4 = RF5/T3
C
      ANGSQLS = T2*T4*T4
C
      RETURN
      END
C
C
      FUNCTION ANGSQJJ(J1,J2,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      AA    NN    NN  GGGGGG   SSSSSS    QQQQQQ       JJJJ     JJJJ   C
C     AAAA   NNN   NN GG    GG SS    SS  QQ    QQ        JJ       JJ   C
C    AA  AA  NNNN  NN GG       SS       QQ      QQ       JJ       JJ   C
C   AA    AA NN NN NN GG        SSSSSS  QQ      QQ       JJ       JJ   C
C   AAAAAAAA NN  NNNN GG   GGG       SS QQ      QQ       JJ       JJ   C
C   AA    AA NN   NNN GG    GG SS    SS  QQ    QQ  JJ    JJ JJ    JJ   C
C   AA    AA NN    NN  GGGGGG   SSSSSS    QQQQQQ QQ JJJJJJ   JJJJJJ    C
C                                                                      C
C -------------------------------------------------------------------- C
C  ANGSQJJ EVALUATES THE SQUARE OF A 3-J SYMBOL,   /  j   K   j' \^2   C
C  WHERE j = J1/2 AND j' = J2/2, FOR THE           \-1/2  0  1/2 /     C
C  COULOMB/BREIT ANGULAR COEFFICIENT ROUTINES.                         C
C**********************************************************************C
      COMMON/FCTS/RFACT(0:80),SFACT(0:80)
C
C     TRIANGLE RULE RESTRICTIONS
      IF(K.LT.IABS((J1-J2)/2).OR.K.GT.(J1+J2)/2) THEN
        ANGSQJJ = 0.0D0
        RETURN
      ELSEIF(J1.LE.0.OR.J2.LE.0) THEN
        ANGSQJJ = 0.0D0
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
      RN1 = RFACT(( J1+J2)/2 - K)
      RN2 = RFACT((-J1+J2)/2 + K)
      RN3 = RFACT(( J1-J2)/2 + K)
      RN4 = SFACT(( J1+J2)/2 + M)
      RD1 = DFLOAT(J1+1)
      RD2 = DFLOAT(J2+1)
      RD3 = RFACT(( J1+J2)/2 + K + 1)
      RD4 = SFACT(( J1+J2)/2 - M    )
      RD5 = SFACT(( J1-J2)/2 + M - 1)
      RD6 = SFACT((-J1+J2)/2 + M - 1)
      PHS = (-1.0D0)**((J2-(3*J1))/2+M)
C
      RNUM  = RN1*RN2*RN3*(RN4)**2
      RDEN  = RD1*RD2*RD3*(RD4*RD5*RD6)**2
C
      ANGSQJJ = PHS*RNUM/RDEN
C
      RETURN
      END
C
C
      SUBROUTINE IJSET0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           IIII     JJJJ SSSSSS  EEEEEEEE TTTTTTTT 000000             C
C            II       JJ SS    SS EE          TT   00   000            C
C            II       JJ SS       EE          TT   00  0000            C
C            II       JJ  SSSSSS  EEEEEE      TT   00 00 00            C
C            II       JJ       SS EE          TT   0000  00            C
C            II JJ    JJ SS    SS EE          TT   000   00            C
C           IIII JJJJJJ   SSSSSS  EEEEEEEE    TT    000000             C
C                                                                      C
C -------------------------------------------------------------------- C
C  IJSET0 GENERATES BASIS SET INTERMEDIATES FOR IJ-PAIRS TO BE USED    C
C  IN THE CONSTRUCTION OF ATOMIC TWO-ELECTRON INTEGRALS.               C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/B0IJ/EIJ(-MTN:MTN),RNIJ(4),EI,EJ,IBAS,JBAS
      COMMON/B0QN/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
C     NORMALISATION CONSTANTS FOR EXPONENTS EI AND EJ
      RL = DFLOAT(LQNA)
      G1 = TWLG-GAMLOG(2*LQNA+3)
      G2 = TWLG-GAMLOG(2*LQNA+5)
      R1 = RL+1.5D0
      R2 = RL+0.5D0
C
      ELOG = DLOG(2.0D0*EI)
      RNLI = DEXP(0.5D0*(G1+R1*ELOG))
      RNSI = DEXP(0.5D0*(G2+R2*ELOG))
C
      ELOG = DLOG(2.0D0*EJ)
      RNLJ = DEXP(0.5D0*(G1+R1*ELOG))
      RNSJ = DEXP(0.5D0*(G2+R2*ELOG))
C
C     NORMALISATION PAIRS
      RNIJ(1) = RNLI*RNLJ
      RNIJ(2) = RNLI*RNSJ
      RNIJ(3) = RNSI*RNLJ
      RNIJ(4) = RNSI*RNSJ
C
C     POWERS OF THE EXPONENT SUM
      EIJ0 = EI+EJ
      EIJR = DSQRT(EIJ0)
      EIJA = EIJ0**(-LQNA)
      DO N=0,MTN
        EIJ(N) = EIJA
        EIJA   = EIJA/EIJR
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE KLSET0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         KK    KK LL       SSSSSS  EEEEEEEE TTTTTTTT 000000           C
C         KK   KK  LL      SS    SS EE          TT   00   000          C
C         KK  KK   LL      SS       EE          TT   00  0000          C
C         KKKKK    LL       SSSSSS  EEEEEE      TT   00 00 00          C
C         KK  KK   LL            SS EE          TT   0000  00          C
C         KK   KK  LL      SS    SS EE          TT   000   00          C
C         KK    KK LLLLLLLL SSSSSS  EEEEEEEE    TT    000000           C
C                                                                      C
C -------------------------------------------------------------------- C
C  KLSET0 GENERATES BASIS SET INTERMEDIATES FOR KL-PAIRS TO BE USED    C
C  IN THE CONSTRUCTION OF ATOMIC TWO-ELECTRON INTEGRALS.               C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/B0QN/EXLA(MBS),EXLB(MBS),NBASA,NBASB,LQNA,LQNB,MAXM
      COMMON/B0IK/EIK(MB2,-MTN:MTN),IKIND(MB2)
      COMMON/B0JL/EJL(MB2,-MTN:MTN),JLIND(MB2)
      COMMON/B0KL/EKL(MB2,-MTN:MTN),RNKL(MB2,4),EK(MB2),EL(MB2)
C
C     GENERATE INDICES AND EXPONENT COMBINATIONS
      M = 0
      DO KBAS=1,NBASB
        EK0 = EXLB(KBAS)
        DO LBAS=1, NBASB
          M   = M+1
          EL0 = EXLB(LBAS)
          IKIND(M) = KBAS
          JLIND(M) = LBAS
          EK(M)    = EK0
          EL(M)    = EL0
          EKL0     = EK0+EL0
          EKLR     = DSQRT(EKL0)
          EKPW     = EKL0**LQNB
          EKLA     = 1.0D0/EKPW
          DO N=0,MTN
            EKL(M,N) = EKLA
            EKLA     = EKLA/EKLR
          ENDDO
        ENDDO
      ENDDO
C
C     NORMALISATION CONSTANTS
      CALL RNORM0(RNKL,EXLB,NBASB,LQNB)
C
      RETURN
      END
C
C
      SUBROUTINE RNORM0(RN,EXL,NBAS,LQN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       RRRRRRR  NN    NN  OOOOOO  RRRRRRR  MM       MM  000000        C
C       RR    RR NNN   NN OO    OO RR    RR MMM     MMM 00   000       C
C       RR    RR NNNN  NN OO    OO RR    RR MMMM   MMMM 00  0000       C
C       RR    RR NN NN NN OO    OO RR    RR MM MM MM MM 00 00 00       C
C       RRRRRRR  NN  NNNN OO    OO RRRRRRR  MM  MMM  MM 0000  00       C
C       RR    RR NN   NNN OO    OO RR    RR MM   M   MM 000   00       C
C       RR    RR NN    NN  OOOOOO  RR    RR MM       MM  000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  RNORM0 EVALUATES NORMALISATION CONSTANTS OF ALL VARIETIES.          C
C  USE THIS ROUTINE IF YOU'RE IN A BLOCK WITH THE SAME IZ AND KQN.     C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RN(MB2,4),EXL(MBS),RNL(MBS),RNS(MBS)
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
      RL = DFLOAT(LQN)
      G1 = TWLG-GAMLOG(2*LQN+3)
      G2 = TWLG-GAMLOG(2*LQN+5)
      R1 = RL+1.5D0
      R2 = RL+0.5D0
      DO IBAS=1,NBAS
        ELOG      = DLOG(2.0D0*EXL(IBAS))
        RNL(IBAS) = DEXP(0.5D0*(G1+R1*ELOG))
        RNS(IBAS) = DEXP(0.5D0*(G2+R2*ELOG))
      ENDDO
C
C     RN(M,1) ARE THE LL NORMALISATION CONSTANTS
C     RN(M,2) ARE THE LS NORMALISATION CONSTANTS
C     RN(M,3) ARE THE SL NORMALISATION CONSTANTS
C     RN(M,4) ARE THE SS NORMALISATION CONSTANTS
C
      M = 0
      DO IBAS=1,NBAS
        DO JBAS=1,NBAS
          M = M+1
          RN(M,1) = RNL(IBAS)*RNL(JBAS)
          RN(M,2) = RNL(IBAS)*RNS(JBAS)
          RN(M,3) = RNS(IBAS)*RNL(JBAS)
          RN(M,4) = RNS(IBAS)*RNS(JBAS)
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
C
C
      FUNCTION SPINDOT(ITT1,ITT2,K1,K2,M1,M2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        SSSSSS  PPPPPPP IIII NN    NN DDDDDDD   OOOOOO TTTTTTTT       C
C       SS    SS PP    PP II  NNN   NN DD    DD OO    OO   TT          C
C       SS       PP    PP II  NNNN  NN DD    DD OO    OO   TT          C
C        SSSSSS  PP    PP II  NN NN NN DD    DD OO    OO   TT          C
C             SS PPPPPPP  II  NN  NNNN DD    DD OO    OO   TT          C
C       SS    SS PP       II  NN   NNN DD    DD OO    OO   TT          C
C        SSSSSS  PP      IIII NN    NN DDDDDDD   OOOOOO    TT          C
C                                                                      C
C -------------------------------------------------------------------- C
C  SPINDOT DETERMINES A DOT PRODUCT OVER ANGULAR SPINOR INTEGRALS:     C
C                                                                      C
C       3                                                              C
C       Σ  A(T1,T1',IQ,κ1,κ2,m1,m2) • A(T2,T2',IQ,κ1,κ2,m1,m2).        C
C     IQ=1                                                             C
C                                                                      C
C  THIS FUNCTION WAS WRITTEN SPECIFICALLY FOR THE SLFLW0 BETHE FORMULA.C
C -------------------------------------------------------------------- C
C  ▶ ITT1  = {1,2,3,4} -> {LL,LS,SL,SS} COMPONENT COMBINATION.         C
C  ▶ ITT2  = {1,2,3,4} -> {LL,LS,SL,SS} COMPONENT COMBINATION.         C
C  ▶ K1,K2 = RELATIVISTIC SYMMETRY (KAPPA) NUMBERS.                    C
C  ▶ M1,M2 = MAGNETIC QUANTUM NUMBERS (INPUT IS DOUBLE THE TRUE VALUE).C
C**********************************************************************C
C
      COMPLEX*16 SPINDOT,SPINANG
C
      SPINDOT = SPINANG(ITT1,1,K1,K2,M1,M2)*SPINANG(ITT2,1,K2,K1,M2,M1)
     &        + SPINANG(ITT1,2,K1,K2,M1,M2)*SPINANG(ITT2,2,K2,K1,M2,M1)
     &        + SPINANG(ITT1,3,K1,K2,M1,M2)*SPINANG(ITT2,3,K2,K1,M2,M1)
C
      RETURN
      END
C
C
      FUNCTION SPINANG(ITT,IQ,K1,K2,M1,M2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       SSSSSS  PPPPPPP IIII NN    NN    AA    NN    NN  GGGGGG        C
C      SS    SS PP    PP II  NNN   NN   AAAA   NNN   NN GG    GG       C
C      SS       PP    PP II  NNNN  NN  AA  AA  NNNN  NN GG             C
C       SSSSSS  PP    PP II  NN NN NN AA    AA NN NN NN GG             C
C            SS PPPPPPP  II  NN  NNNN AAAAAAAA NN  NNNN GG   GGG       C
C      SS    SS PP       II  NN   NNN AA    AA NN   NNN GG    GG       C
C       SSSSSS  PP      IIII NN    NN AA    AA NN    NN  GGGGGG        C
C                                                                      C
C -------------------------------------------------------------------- C
C  SPINANG DETERMINES THE INTEGRAL OVER TWO ANGULAR SPINORS AND A      C
C  PAULI MATRIX, INVOLVING COMPONENT TYPES AND CARTESIAN MATRICES.     C
C                                                                      C
C                                     ↥T            T'                 C
C       A(T,T',IQ,κ1,κ2,m1,m2) = ∬   χ   (θ,φ) σ   χ   (θ,φ) dΩ        C
C                                 4π  κ1,m1     IQ  κ2,m2              C
C                                                                      C
C -------------------------------------------------------------------- C
C  ▶ ITT   = {1,2,3,4} -> {LL,LS,SL,SS} COMPONENT COMBINATION.         C
C  ▶ IQ    = {0,1,2,3} THE PAULI COUPLING MATRIX.                      C
C  ▶ K1,K2 = RELATIVISTIC SYMMETRY (KAPPA) NUMBERS.                    C
C  ▶ M1,M2 = MAGNETIC QUANTUM NUMBERS (INPUT IS DOUBLE THE TRUE VALUE).C
C**********************************************************************C
C
      COMPLEX*16 SPINANG
      COMPLEX*16 CONE
C
C     UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     INITIAL VALUE FOR THIS INTEGRAL
      SPINANG = DCMPLX(0.0D0,0.0D0)
C
C     INITIAL VALUES FOR COEFFICIENTS
C
C     ORBITAL AND PARITY NUMBERS
      IF(ITT.EQ.1.OR.ITT.EQ.2) THEN
        L1 = LVAL( K1)
        N1 = K1/IABS(K1)
      ELSE
        L1 = LVAL(-K1)
        N1 =-K1/IABS(K1)
      ENDIF
C
      IF(ITT.EQ.1.OR.ITT.EQ.3) THEN
        L2 = LVAL( K2)
        N2 = K2/IABS(K2)
      ELSE
        L2 = LVAL(-K2)
        N2 =-K2/IABS(K2)
      ENDIF
C
C     LQN SELECTION RULES
      IF(L1.NE.L2) RETURN
C
C     COEFFICIENTS BASED ON MQN SELECTION RULES
      C11 = DCMPLX(0.0D0,0.0D0)
      C22 = DCMPLX(0.0D0,0.0D0)
      IF(M1.EQ.M2) THEN
        C11 = N1*N2*CLEBSCH(L1,M1,-N1)*CLEBSCH(L2,M2,-N2)
        C22 =       CLEBSCH(L1,M1, N1)*CLEBSCH(L2,M2, N2)
      ENDIF
C
      C12 = DCMPLX(0.0D0,0.0D0)
      IF(M1.EQ.M2+2) THEN
        C12 =-N1*   CLEBSCH(L1,M1,-N1)*CLEBSCH(L2,M2, N2)
      ENDIF
C
      C21 = DCMPLX(0.0D0,0.0D0)
      IF(M1.EQ.M2-2) THEN
        C21 =-   N2*CLEBSCH(L1,M1, N1)*CLEBSCH(L2,M2,-N2)
      ENDIF
C
C     IMPLEMENT PAULI MATRICES
      IF(IQ.EQ.0) THEN
        SPINANG = C11+C22
      ELSEIF(IQ.EQ.1) THEN
        SPINANG = C12+C21
      ELSEIF(IQ.EQ.2) THEN
        SPINANG =-C12+C21
        SPINANG = CONE*SPINANG
      ELSEIF(IQ.EQ.3) THEN
        SPINANG = C11-C22
      ENDIF
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
      IF(KQN.LT.0) THEN
        LQN =-KQN-1
      ELSE
        LQN = KQN
      ENDIF
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
      SUBROUTINE GNULINE(XOUT,TITLE,XAXIS,YAXIS,NDAT,KEY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       GGGGGG  NN    NN UU    UU LL       IIII NN    NN EEEEEEEE      C
C      GG    GG NNN   NN UU    UU LL        II  NNN   NN EE            C
C      GG       NNNN  NN UU    UU LL        II  NNNN  NN EE            C
C      GG       NN NN NN UU    UU LL        II  NN NN NN EEEEEE        C
C      GG   GGG NN  NNNN UU    UU LL        II  NN  NNNN EE            C
C      GG    GG NN   NNN UU    UU LL        II  NN   NNN EE            C
C       GGGGGG  NN    NN  UUUUUU  LLLLLLLL IIII NN    NN EEEEEEEE      C
C                                                                      C
C -------------------------------------------------------------------- C
C  GNULINE IS A CONTROLLING ROUTINE THAT GENERATES A GNUPLOT MAKE FILE C
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
      WRITE(9,'(A)') '# Load line style definitions'
      WRITE(9,'(A)') 'load "plots/plotstyles.pal"'
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
      WRITE(6, *) 'Created command file "'//TRIM(XOUT)//'.gnuplot".'
      WRITE(7, *) 'Created command file "'//TRIM(XOUT)//'.gnuplot".'
C
      RETURN
      END
C
C
      SUBROUTINE DGNUMAP(ARRAY,TITLE,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  DDDDDDD   GGGGGG  NN    NN UU    UU MM       MM    AA    PPPPPPP    C
C  DD    DD GG    GG NNN   NN UU    UU MMM     MMM   AAAA   PP    PP   C
C  DD    DD GG       NNNN  NN UU    UU MMMM   MMMM  AA  AA  PP    PP   C
C  DD    DD GG       NN NN NN UU    UU MM MM MM MM AA    AA PP    PP   C
C  DD    DD GG   GGG NN  NNNN UU    UU MM  MMM  MM AAAAAAAA PPPPPPP    C
C  DD    DD GG    GG NN   NNN UU    UU MM   M   MM AA    AA PP         C
C  DDDDDDD   GGGGGG  NN    NN  UUUUUU  MM       MM AA    AA PP         C
C                                                                      C
C -------------------------------------------------------------------- C
C  DGNUMAP EXPORTS AN ARRAY OF DOUBLE-PRECISION NUMBERS TO AN EXTERNAL C
C  DATA FILE AND PLOTS IT AS ONE HEAT MAP.                             C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*80 TITLE
C
      DIMENSION ARRAY(NDIM,NDIM)
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
      WRITE(9,'(A)') 'load "plots/pals/jet2.pal"'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Axes and title'
      WRITE(9,'(A)') 'set title sprintf("'//TRIM(TITLE)//'")'
      WRITE(9,'(A,F6.1,A)') 'set xrange [-0.5:',XEND,']'
      WRITE(9,'(A,F6.1,A)') 'set yrange [',YEND,':-0.5] reverse'
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
      SUBROUTINE RADPLOT(IOCC,R0,RF,NTOT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     RRRRRRR     AA    DDDDDDD  PPPPPPP  LL       OOOOOO TTTTTTTT     C
C     RR    RR   AAAA   DD    DD PP    PP LL      OO    OO   TT        C
C     RR    RR  AA  AA  DD    DD PP    PP LL      OO    OO   TT        C
C     RR    RR AA    AA DD    DD PP    PP LL      OO    OO   TT        C
C     RRRRRRR  AAAAAAAA DD    DD PPPPPPP  LL      OO    OO   TT        C
C     RR    RR AA    AA DD    DD PP       LL      OO    OO   TT        C
C     RR    RR AA    AA DDDDDDD  PP       LLLLLLLL OOOOOO    TT        C
C                                                                      C
C -------------------------------------------------------------------- C
C  RADPLOT TELLS YOU ABOUT THE RADIAL WAVEFUNCTION OF AN ORBITAL.      C
C  INPUT IS IN ATOMIC UNITS. I WROTE THIS ASSUMING ATOMIC CODE.        C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(2)
C
      COMPLEX*16 COEF(MDM,MDM)
C
      DIMENSION EXL(MBS)
      DIMENSION RAD(0:NTOT),PR(0:NTOT),QR(0:NTOT)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/EIGC/COEF
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
C
C     GENERATE RADIAL GRID (LINEAR SCALE) AND INITIALISE AMPLITUDES
      HSPC = (RF-R0)/DFLOAT(NTOT)
      DO N=0,NTOT
        RAD(N) = R0 + HSPC*DFLOAT(N)
        PR(N)  = 0.0D0
        QR(N)  = 0.0D0
      ENDDO
C
C     LOOP OVER KQN VALUES
      DO 2000 KA=1,NKAP(1)
C
C       ORBITAL QUANTUM NUMBERS
        KQN = KAPA(KA,1)
        JQN = 2*IABS(KQN)-1
        IF(KQN.LT.0) THEN
          LQN =-KQN-1
        ELSE
          LQN = KQN
        ENDIF
C
C       BASIS EXPONENTS
        NBAS = NFNC(LQN,1)
        DO IBAS=1,NBAS
          EXL(IBAS) = BEXL(IBAS,LQN,1)
        ENDDO
C
C     LOOP OVER |MQN| VALUES
      DO 3000 MA=1,IABS(KQN)
        MJA = 2*MA-1
        MQN = MJA
C
C     EXPANSION COEFFICIENT MATRIX STARTING ADDRESSES
      NA1 = LRGE(1,KA,2*MA-1)
      NA2 = LRGE(1,KA,2*MA  )
C
C     BEGIN LOOP OVER ALL GRID COORDINATES
      DO N=0,NTOT
C
C       RADIAL COORDINATE
        R = RAD(N)
C
C       PRE-FACTORS FOR NORMALISATION CONSTANTS
        RL = DFLOAT(LQN)
        G1 = TWLG-GAMLOG(2*LQN+3)
        G2 = TWLG-GAMLOG(2*LQN+5)
        R1 = RL+1.5D0
        R2 = RL+0.5D0
C
C       LOOP OVER BASIS FUNCTIONS
        DO IBAS=1,NBAS
C
C         NORMALISATION CONSTANTS
          ELOG = DLOG(2.0D0*EXL(IBAS))
          RNL  = DEXP(0.5D0*(G1+R1*ELOG))
          RNS  = DEXP(0.5D0*(G2+R2*ELOG))
C
C         LARGE AND SMALL RADIAL AMPLITUDES
CC          FL = RNL*(R**(LQN+1))*DEXP(-EXL(IBAS)*R*R)
CC          FS = RNS*(KQN+LQN+1.0D0-2.0D0*EXL(IBAS)*R*R)
CC     &                          *(R**(LQN))*DEXP(-EXL(IBAS)*R*R)
          
          IF(LQN.EQ.0) THEN
            FL = RNL*DEXP(-EXL(IBAS)*R*R)
            FS = RNS*(-2.0D0*EXL(IBAS))*DEXP(-EXL(IBAS)*R*R)
          ELSE
            FL = RNL*(R**(LQN))*DEXP(-EXL(IBAS)*R*R)
            FS = RNS*(KQN+LQN+1.0D0-2.0D0*EXL(IBAS)*R*R)
     &                            *(R**(LQN-1))*DEXP(-EXL(IBAS)*R*R)
          ENDIF
C
C         ADDRESS OFFSETS FOR SMALL-COMPONENTS
          KBAS = IBAS+NSKP
C
C         MULTIPLY BY CORRESPONDING EXPANSION COEFFICIENT ELEMENTS
          PR(N) = PR(N) + DREAL(COEF(NA1+IBAS,NSKP+IOCC))*FL
     &                  + DREAL(COEF(NA2+IBAS,NSKP+IOCC))*FL
C
          IF(HMLT.EQ.'NORL') GOTO 100
C
          QR(N) = QR(N) + DREAL(COEF(NA1+KBAS,NSKP+IOCC))*FS
     &                  + DREAL(COEF(NA2+KBAS,NSKP+IOCC))*FS
C          
100       CONTINUE
C
        ENDDO
C
C     CLOSE LOOP OVER GRID COORDINATES
      ENDDO
C
C     END LOOP OVER BASIS QUANTUM NUMBERS
3000  CONTINUE
2000  CONTINUE
C
C     FILE NAME AND TITLES
      XOUT   = 'RadialAmplitude'
      TITLE  = 'Radial Dirac amplitude'
      XAXIS  = 'r [fm]'
      YAXIS  = 'R(r)'
      KEY(1) = 'Large'
      KEY(2) = 'Small'
C
C     PRINT DATA TO EXTERNAL FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
        DO N=0,NTOT
          WRITE(8, *) RAD(N)*CFM,PR(N)/CFM,QR(N)/CFM
        ENDDO
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
      RETURN
      END
C
C
      SUBROUTINE RATPLOT(IOCC,R0,RF,NTOT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     RRRRRRR     AA    TTTTTTTT PPPPPPP  LL       OOOOOO TTTTTTTT     C
C     RR    RR   AAAA      TT    PP    PP LL      OO    OO   TT        C
C     RR    RR  AA  AA     TT    PP    PP LL      OO    OO   TT        C
C     RR    RR AA    AA    TT    PP    PP LL      OO    OO   TT        C
C     RRRRRRR  AAAAAAAA    TT    PPPPPPP  LL      OO    OO   TT        C
C     RR    RR AA    AA    TT    PP       LL      OO    OO   TT        C
C     RR    RR AA    AA    TT    PP       LLLLLLLL OOOOOO    TT        C
C                                                                      C
C -------------------------------------------------------------------- C
C  RATPLOT TELLS YOU ABOUT THE RADIAL WAVEFUNCTION RATIOS OF IOCC.     C
C  INPUT IS IN ATOMIC UNITS. I WROTE THIS ASSUMING ATOMIC CODE.        C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(1)
C
      COMPLEX*16 COEF(MDM,MDM)
C
      DIMENSION EXL(MBS)
      DIMENSION RAD(0:NTOT),PR(0:NTOT),QR(0:NTOT)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/EIGC/COEF
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
C
C     GENERATE RADIAL GRID (LINEAR SCALE) AND INITIALISE AMPLITUDES
      HSPC = (RF-R0)/DFLOAT(NTOT)
      DO N=0,NTOT
        RAD(N) = R0 + HSPC*DFLOAT(N)
        PR(N)  = 0.0D0
        QR(N)  = 0.0D0
      ENDDO
C
C     LOOP OVER KQN VALUES
      DO 2000 KA=1,NKAP(1)
C
C       ORBITAL QUANTUM NUMBERS
        KQN = KAPA(KA,1)
        JQN = 2*IABS(KQN)-1
        IF(KQN.LT.0) THEN
          LQN =-KQN-1
        ELSE
          LQN = KQN
        ENDIF
C
C       BASIS EXPONENTS
        NBAS = NFNC(LQN,1)
        DO IBAS=1,NBAS
          EXL(IBAS) = BEXL(IBAS,LQN,1)
        ENDDO
C
C     LOOP OVER |MQN| VALUES
      DO 3000 MA=1,IABS(KQN)
        MJA = 2*MA-1
        MQN = MJA
C
C     EXPANSION COEFFICIENT MATRIX STARTING ADDRESSES
      NA1 = LRGE(1,KA,2*MA-1)
      NA2 = LRGE(1,KA,2*MA  )
C
C     BEGIN LOOP OVER ALL GRID COORDINATES
      DO N=0,NTOT
C
C       RADIAL COORDINATE
        R = RAD(N)
C
C       PRE-FACTORS FOR NORMALISATION CONSTANTS
        RL = DFLOAT(LQN)
        G1 = TWLG-GAMLOG(2*LQN+3)
        G2 = TWLG-GAMLOG(2*LQN+5)
        R1 = RL+1.5D0
        R2 = RL+0.5D0
C
C       LOOP OVER BASIS FUNCTIONS
        DO IBAS=1,NBAS
C
C         NORMALISATION CONSTANTS
          ELOG = DLOG(2.0D0*EXL(IBAS))
          RNL  = DEXP(0.5D0*(G1+R1*ELOG))
          RNS  = DEXP(0.5D0*(G2+R2*ELOG))
C
C         LARGE AND SMALL RADIAL AMPLITUDES
CC          FL = RNL*(R**(LQN+1))*DEXP(-EXL(IBAS)*R*R)
CC          FS = RNS*(KQN+LQN+1.0D0-2.0D0*EXL(IBAS)*R*R)
CC     &                          *(R**(LQN))*DEXP(-EXL(IBAS)*R*R)
C
C         LARGE- AND SMALL-AMPLITUDE P0 AND Q0
          GSS = DEXP(-EXL(IBAS)*R*R)
          IF(KQN.LT.0) THEN
            FL = RNL*DEXP(-EXL(IBAS)*R*R)
            FS = RNS*(-2.0D0*EXL(IBAS))*GSS
          ELSE
            FL = RNL*DEXP(-EXL(IBAS)*R*R)
            FS = RNS*(KQN+LQN+1.0D0-2.0D0*EXL(IBAS)*R*R)*GSS
          ENDIF
C
C         ADDRESS OFFSETS FOR SMALL-COMPONENTS
          KBAS = IBAS+NSKP
C
C         MULTIPLY BY CORRESPONDING EXPANSION COEFFICIENT ELEMENTS
          PR(N) = PR(N) + DREAL(COEF(NA1+IBAS,NSKP+IOCC))*FL
     &                  + DREAL(COEF(NA2+IBAS,NSKP+IOCC))*FL
          QR(N) = QR(N) + DREAL(COEF(NA1+KBAS,NSKP+IOCC))*FS
     &                  + DREAL(COEF(NA2+KBAS,NSKP+IOCC))*FS
C
        ENDDO
C
C     CLOSE LOOP OVER GRID COORDINATES
      ENDDO
C
C     END LOOP OVER BASIS QUANTUM NUMBERS
3000  CONTINUE
2000  CONTINUE
C
C     FILE NAME AND TITLES
      XOUT   = 'RadialAmplitudeRatio'
      TITLE  = 'Radial Dirac ratio'
      XAXIS  = 'r [fm]'
      YAXIS  = 'R(r)'
      KEY(1) = 'Ratio'
C
C     PRINT DATA TO EXTERNAL FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
        DO N=0,NTOT
          WRITE(8, *) RAD(N)*CFM,QR(N)/PR(N)
        ENDDO
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,1,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
      RETURN
      END
