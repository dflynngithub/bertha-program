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
C  * NN    NN UU    UU  CCCCCC  LL       EEEEEEEE    AA    RRRRRRR  *  C
C    NNN   NN UU    UU CC    CC LL       EE         AAAA   RR    RR    C
C  * NNNN  NN UU    UU CC       LL       EE        AA  AA  RR    RR *  C
C    NN NN NN UU    UU CC       LL       EEEEEE   AA    AA RR    RR    C
C  * NN  NNNN UU    UU CC       LL       EE       AAAAAAAA RRRRRRR  *  C
C    NN   NNN UU    UU CC    CC LL       EE       AA    AA RR    RR    C
C  * NN    NN  UUUUUU   CCCCCC  LLLLLLLL EEEEEEEE AA    AA RR    RR *  C
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
C     OPEN FILE FOR TERMINAL RECORD
      OPEN(UNIT=7,FILE=TRIM(OUTFL)//'_nuclear.out',STATUS='UNKNOWN')
C
C     PRINT SUMMARY OF INPUT DATA
      CALL INPUT
C
C     NUCLEAR POTENTIALS
      CALL CPU_TIME(T0)
      CALL NUCCOUL
      CALL NUCVCPL
      CALL CPU_TIME(TNUC)
      TNUC = TNUC-T0
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
C
C     CALCULATION TREE CHECK -- IF UNKNOWN, EXIT.
      IF(TREE.NE.'HFSCF'.AND.TREE.NE.'MBPTN'.AND.TREE.NE.'MCSCF'.AND.
     &   TREE.NE.'EXPVL'.AND.TREE.NE.'PLOTS') THEN
        WRITE(6, *) 'In CARDIN: invalid calculation tree. ',TREE
        WRITE(7, *) 'In CARDIN: invalid calculation tree. ',TREE
c       STOP
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
c         STOP
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
c         STOP
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
c             STOP
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
c             STOP
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
c             STOP
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
c       STOP
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
c       STOP
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
c       STOP
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
c         STOP
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
c         STOP
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
cC     READ IN A WAVE FUNCTION FILE IF PROMPTED
c      IF(READIN) THEN
c        INQUIRE(FILE=WFNFL,EXIST=FILETHERE)
c        IF(FILETHERE) THEN
c          OPEN(UNIT=8,FILE=WFNFL,STATUS='UNKNOWN')
c          REWIND(UNIT=8)
c          DO I=1,NDIM
c            READ(8, *) EIGN(I),(COEF(J,I),J=1,NDIM)
c          ENDDO
c          CLOSE(UNIT=8)
c        ELSE
c          READIN = .FALSE.
c        ENDIF
c      ENDIF
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
      WRITE(6, *) REPEAT(' ',29),'BERTHA-NUCLEAR'
      WRITE(7, *) REPEAT(' ',29),'BERTHA-NUCLEAR'
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
c       STOP
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
30    FORMAT(1X,A,18(' '),A)
      WRITE(6, *)
      WRITE(7, *)
      WRITE(6,30) 'Successful BERTHA-NUCLEAR exit at:',STAMP
      WRITE(7,30) 'Successful BERTHA-NUCLEAR exit at:',STAMP
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
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/MDLV/ELMT
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/BMSS/RMPC,RMP0,RMKC,RMK0
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
C     RIEMANN ZETA FUNCTION
      DATA ZTA3/1.2020569031595943D0/
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
C     RELATIVE MESON MASSES (ATOMIC UNITS)
      DATA RMPC,RMP0,RMKC,RMK0/263.35D0,264.14D0,973.9D0,973.8D0/
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
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
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
C   [B] VFERMI: NORMALISED FERMI NUCLEAR POTENTIAL AT RADIUS R.        C
C   [C] POLYLOG: THE POLYLOGARITHM OF NEGATIVE EXPONENTIAL ARGUMENT.   C
C   [D] NUCGEOM: BOND DISTANCES AND NUCLEAR REPULSION ENERGY.          C
C   [E] FOCKIND: CALCULATE ADDRESSES OF FOCK MATRIX FOR BASIS QN'S.    C
C   [F] CARTIND: GENERATES INDICES FOR EQ-COEFFS AND R-INTEGRALS.      C
C   [G] AUFBAU: DETERMINES GROUND STATE ATOMIC ELECTRON CONFIG.        C
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
C  ▶ NCOUL = 0 (POINT), NCOUL = 1 (GAUSS), NCOUL > 1 (FERMI).          C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
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
      COMMON/MDLV/ELMT
C
C     TOGGLE FOR MORE DETAILED SUMMARIES (0 FOR OFF, 1 FOR ON)
      IWRT = 0
C
C     NUMBER OF FITTING FUNCTIONS FOR COULOMB PROCEDURE
      DO IZ=1,NCNT
        IF(NMDL(IZ).EQ.'POINT'.OR.NMDL(IZ).EQ.'GAUSS') THEN
          NCOUL(IZ) = 0
        ELSEIF(NMDL(IZ).EQ.'UNIFM'.OR.NMDL(IZ).EQ.'FERMI') THEN
          NCOUL(IZ) = 12
        ELSE
          WRITE(6, *) 'In NUCCOUL: invalid charge model. ',NMDL(IZ)
          WRITE(7, *) 'In NUCCOUL: invalid charge model. ',NMDL(IZ)
        ENDIF
      ENDDO
C
C     NUCLEAR COULOMB POTENTIAL BEST-FIT PROCEDURE
      DO IZ=1,NCNT
C
C       GEOMETRIC PARAMETER SEARCH DETAILS
        IF(ZNUC(IZ).LE.10.0D0) THEN
          ALPH = 1.00D-1
          BETA = 1.40D+0
          GAMA =-0.30D+0
          RMAX = 2.00D+0
        ELSEIF(ZNUC(IZ).GT.10.0D0.AND.ZNUC(IZ).LE.90.0D0) THEN
          ALPH = 2.00D-1
          BETA = 1.40D+0
          GAMA =-0.40D+0
          RMAX = 2.50D+0
c          ALPH = 2.00D+0
c          BETA = 1.60D+0
c          GAMA = 0.21D+0
c          RMAX = 2.50D+0
        ELSE
          ALPH = 3.00D-1
          BETA = 1.30D+0
          GAMA =-0.21D+0
          RMAX = 2.50D+0
        ENDIF
C
C       CALL THE FITTING PROCEDURE
        CALL VCLMGEN(IZ,NCOUL(IZ),IWRT,RSQ(IZ),RMAX,ALPH,BETA,GAMA)
C
      ENDDO
CC
CC     EIGENVALUE SPECTRUM
C      DO I=0,100
C        ALPH = 2.0D-2 + DFLOAT(I)*1.0D-2/DFLOAT(100)
C        CALL EIGBASIS(10,ALPH,BETA,GAMA)
C      ENDDO
C
C     WRITE RESULTS OF COULOMB BEST-FIT PROCEDURE
21    FORMAT(26X,A)
22    FORMAT(1X,A,2X,'|',3X,A,4X,A,5X,A,1X,'|',1X,A,3X,A,3X,A)
23    FORMAT(1X,I2,' (',A,')',1X,'|',5X,I3,4X,F8.4,5X,F6.4,1X,'|',
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
     &            'N_fit','Potential','R-squared'
      WRITE(7,22) 'Centre','Z (e)','A (m_p+)','R (fm)',
     &            'N_fit','Potential','R-squared'
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
C     NUCLEUS FILE NAME
      NUCFL = 'output/'//TRIM(MOLCL)//'_nuclear.dat'
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE=TRIM(NUCFL),STATUS='UNKNOWN')
C
C       LOOP OVER NUCLEAR CENTERS
        DO IZ=1,NCNT
C
C         NUMBER OF FITTING FUNCTIONS AND RMS RADIUS
          WRITE(8, *) NNUC(IZ),RNUC(IZ)
C
C         AMPLITUDES AND PARAMETERS FOR EACH NUCLEAR FITTING FUNCTION
          DO IFT=0,NNUC(IZ)
            WRITE(8, *) FNUC(IZ,IFT),XNUC(IZ,IFT)
          ENDDO
C
        ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
      RETURN
      END
C
C
      SUBROUTINE EIGBASIS(NFT,ALPH,BETA,GAMA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  VV    VV  CCCCCC  LL       MM       MM  GGGGGG  EEEEEEEE NN    NN   C
C  VV    VV CC    CC LL       MMM     MMM GG    GG EE       NNN   NN   C
C  VV    VV CC       LL       MMMM   MMMM GG       EE       NNNN  NN   C
C  VV    VV CC       LL       MM MM MM MM GG       EEEEEE   NN NN NN   C
C   VV  VV  CC       LL       MM  MMM  MM GG   GGG EE       NN  NNNN   C
C    VVVV   CC    CC LL       MM   M   MM GG    GG EE       NN   NNN   C
C     VV     CCCCCC  LLLLLLLL MM       MM  GGGGGG  EEEEEEEE NN    NN   C
C                                                                      C
C -------------------------------------------------------------------- C
C  VCLMGEN EVALUATES THE DETAILS OF NFT GAUSSIAN NUCLEAR FUNCTIONS     C
C  FOR EACH CENTRE IN THE MOLECULE, SO AS TO BEST FIT THE FERMI MODEL. C
C  EXPONENTS XNUC ARE SELECTED MANUALLY, AND FRACTIONS FNUC CALCULATED.C
C -------------------------------------------------------------------- C
C  MATCHING CRITERIA AVAILABLE FOR GAUSSIAN AMPLITUDES:                C
C   ▶ TOTAL NORMALISED CHARGE OF UNITY <R^0>.                          C
C   ▶ SECOND MOMENT <R^2> SAME AS EMPIRICAL FORMULA OR AVAILABLE DATA. C
C   ▶ ANY OTHER INTEGER MOMENT CAN BE MATCHED WITH FERMI MODEL.        C
C   ▶ ANY OTHER REAL VALUED MOMENT (K > -3) CAN BE MATCHED WITH DATA.  C
C   ▶ FORCED AGREEMENT WITH FERMI POTENTIAL VFERMI(R) AT ANY RADIUS R. C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5  NMDL
C
      DIMENSION X(NFT,NFT),Z(NFT),Y(NFT),XI(NFT)
      DIMENSION IPIV(NFT)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
C     GEOMETRIC SEED VALUES
      AF = ALPH
      BF = BETA
      GF = GAMA
C
C     SET OF GAUSSIAN EXPONENTS
      XI(1) = AF/(RNUC(1)**2)
      DO I=1,NFT-1
        RT = DFLOAT(I+1)/DFLOAT(NFT+1)
        XI(I+1) = BF*XI(I)*(1.0D0 + GAMA*RT*RT)
      ENDDO
CC
CC     PRINT EXPONENTS      
C      WRITE(*,*) 'EXPONENTS'
C      WRITE(*,*) RNUC(1)
C      DO I=1,NFT
C        WRITE(*,*) I,XI(I)*(RNUC(1)**2)
C      ENDDO
C
C     POTENTIAL IN NUCLEAR REGION
      DO I=1,NFT
        DO J=1,NFT
          XPR = XI(I)*XI(J)
          XSM = XI(I)+XI(J)
          XPR = XPR**(1.75D0)
          XSM = XSM**(3.50D0)
          X(I,J) = 8.0D0*TW12*XPR/XSM
C          XPR = XPR*XPR
C          XSM = XSM**(3.50D0)
C          X(I,J) = 15.0D0*XPR/(4.0D0*PI12*XSM)
        ENDDO
      ENDDO
C
111   FORMAT(10(F4.1,1X))
C
C     SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
      CALL dgetrf(NFT,NFT,X,NFT,IPIV,INFO)
      DET = 1.0D0
      DO I=1,NFT
        DET = DET*X(I,I)
      ENDDO
      WRITE(*,*) 'DET = ',AF,DET
C     CALL dsyevd('V','U',1,X,NFT,Z,Y,NFT,IPIV,NFT,INFO)
C     CALL DSYGV(1,'V','U',NMAT,C,MBD,O,MBD,E,TE,LWK0,INFO)
C     CALL DGESV(NFT,1,X,NFT,IPIV,Z,NFT,INFO)
C
      RETURN
      END
C
C
      SUBROUTINE VCLMGEN(IZ,NFT,IWRT,RSQBIG,RMAX,ALPH,BETA,GAMA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  VV    VV  CCCCCC  LL       MM       MM  GGGGGG  EEEEEEEE NN    NN   C
C  VV    VV CC    CC LL       MMM     MMM GG    GG EE       NNN   NN   C
C  VV    VV CC       LL       MMMM   MMMM GG       EE       NNNN  NN   C
C  VV    VV CC       LL       MM MM MM MM GG       EEEEEE   NN NN NN   C
C   VV  VV  CC       LL       MM  MMM  MM GG   GGG EE       NN  NNNN   C
C    VVVV   CC    CC LL       MM   M   MM GG    GG EE       NN   NNN   C
C     VV     CCCCCC  LLLLLLLL MM       MM  GGGGGG  EEEEEEEE NN    NN   C
C                                                                      C
C -------------------------------------------------------------------- C
C  VCLMGEN EVALUATES THE DETAILS OF NFT GAUSSIAN NUCLEAR FUNCTIONS     C
C  FOR EACH CENTRE IN THE MOLECULE, SO AS TO BEST FIT THE FERMI MODEL. C
C  EXPONENTS XNUC ARE SELECTED MANUALLY, AND FRACTIONS FNUC CALCULATED.C
C -------------------------------------------------------------------- C
C  MATCHING CRITERIA AVAILABLE FOR GAUSSIAN AMPLITUDES:                C
C   ▶ TOTAL NORMALISED CHARGE OF UNITY <R^0>.                          C
C   ▶ SECOND MOMENT <R^2> SAME AS EMPIRICAL FORMULA OR AVAILABLE DATA. C
C   ▶ ANY OTHER INTEGER MOMENT CAN BE MATCHED WITH FERMI MODEL.        C
C   ▶ ANY OTHER REAL VALUED MOMENT (K > -3) CAN BE MATCHED WITH DATA.  C
C   ▶ FORCED AGREEMENT WITH FERMI POTENTIAL VFERMI(R) AT ANY RADIUS R. C
C**********************************************************************C
      INCLUDE 'parameters.h'
      PARAMETER(NPTS=6000)
C
      CHARACTER*5  NMDL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(4),KEYNFT(NFT)
C
      DIMENSION CFMI(MCT)
      DIMENSION RADS(0:NPTS)
      DIMENSION VNUC(4,0:NPTS),PNUC(4,0:NPTS),ENUC(4,0:NPTS)
      DIMENSION PAC(11),RAC(11),PRC(11),PLG(11)
      DIMENSION RF(-2:8),RG(-2:8),R0(-2:8),RU(-2:8)
      DIMENSION VU(0:NPTS),VF(0:NPTS),VG(0:NPTS),V0(0:NPTS)
      DIMENSION RM(NFT),PM(NFT)
C
      DIMENSION X(NFT,NFT),Z(NFT),Y(NFT)
      DIMENSION IPIV(NFT)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     MULTIPLES OF CONVERSION FACTOR
      C2 = CFM*CFM
      C4 = C2*CFM*CFM
      C6 = C4*CFM*CFM
C
C     NUMBER OF FITTING FUNCTIONS
      NNUC(IZ) = NFT
C
C     NUCLEAR RADIUS
      IF(INT(ZNUC(IZ)).EQ.1) THEN
C       PROTON HAS A VERY SMALL RADIUS (IN FM)
        RFM = PRAD
      ELSEIF(ZNUC(IZ).GT.1.0D0.AND.ZNUC(IZ).LT.190.0D0) THEN
C       EMPIRICAL FORMULA FOR NUCLEAR RMS RADIUS IF Z < 90 (IN FM)
        RFM = 0.836D0*(ANUC(IZ)**(1.0D0/3.0D0)) + 0.57D0
C       RFM = 1.20D0*(ANUC(IZ)**(1.0D0/3.0D0))
      ELSE
C       EMPIRICAL FORMULA FOR NUCLEAR RMS RADIUS IF Z > 90 (IN FM)
        RFM = 0.770D0*(ANUC(IZ)**(1.0D0/3.0D0)) + 0.98D0
      ENDIF
C
C     CONVERT RESULT TO ATOMIC UNITS
      RNUC(IZ) = RFM/CFM
C
C     POINT-NUCLEUS MODEL
      IF(NMDL(IZ).EQ.'POINT') THEN
        RNUC(IZ) = 0.0D0
        RSQBIG   = 1.0D0
        GOTO 20
      ENDIF
C
C**********************************************************************C
C     SET UP THE RADIAL TESTING GRID                                   C
C**********************************************************************C
C
C     COULOMB POTENTIAL MATCHING VALUES
      HD = RMAX/DFLOAT(NFT-4)
      DO K=4,NFT
        RM(K) = RNUC(IZ)*HD*(K-4)
      ENDDO
C
C**********************************************************************C
C     GAUSSIAN NUCLEAR MODEL DETAILS                                   C
C**********************************************************************C
C
100   CONTINUE
C
C     ELEMENTARY GAUSSIAN FITTING
      FNUC(IZ,0) = 1.0D0
      XNUC(IZ,0) = 1.5D0/(RNUC(IZ)*RNUC(IZ))
C
C     POTENTIAL ON THE GRID FOR GAUSSIAN MODEL
      DO IPTS=0,NPTS
        R = 5.0D0*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
        V0(IPTS) = VGAUSS(R,XNUC(IZ,0),1.0D0)
      ENDDO
C
C     RADIAL MOMENTS FOR A GAUSSIAN MODEL
      DO I=-2,8
        FC1 = 4.0D0/PI12
        FC2 = (2.0D0/3.0D0)**(0.5D0*I)
        FC3 = RNUC(IZ)**I
        FC4 = GAMHLF(5+I)
        FC5 = 3.0D0+I
        R0(I) = FC1*FC2*FC3*FC4/FC5
      ENDDO
C
C     GAUSSIAN MODEL
      IF(NMDL(IZ).EQ.'GAUSS') THEN
        RSQBIG   = 1.0D0
        NNUC(IZ) = 0
        GOTO 20
      ENDIF
C
C**********************************************************************C
C     FERMI NUCLEAR MODEL DETAILS                                      C
C**********************************************************************C
C
C     FERMI SKIN THICKNESS AND HALF-DENSITY PARAMETERS (IN A0)
      TFMI(IZ) = 2.30D0/CFM
      AFMI(IZ) = 0.25D0*TFMI(IZ)/THLG
      CFMI(IZ) = 5.0D0*RNUC(IZ)*RNUC(IZ)/3.0D0
     &         - 7.0D0*PI*PI*AFMI(IZ)*AFMI(IZ)/3.0D0
      CFMI(IZ) = DSQRT(CFMI(IZ))
C
C     CHECK WHETHER FERMI MODEL IS POSSIBLE
      EPS = DSQRT(1.4D0)*PI*0.25D0*TFMI(IZ)/THLG
      IF(EPS.GT.RNUC(IZ)) THEN
        NMDL(IZ) = 'GAUSS'
        GOTO 100
      ENDIF
C
C     POTENTIAL ON THE GRID FOR FERMI MODEL
      DO IPTS=0,NPTS
        R = 5.0D0*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
        VF(IPTS) = VFERMI(R,CFMI(IZ),TFMI(IZ))
      ENDDO
C
C     THESE PARAMETERS MAKE FOR NEATER EXPRESSIONS
      PST = 1.0D0
      RST = 1.0D0
      DO I=1,11
        PST = PST*PI*AFMI(IZ)/CFMI(IZ)
        RST = RST*AFMI(IZ)/CFMI(IZ)
        PAC(I) = PST
        RAC(I) = RST
        PLG(I) = POLYLOG(I,-1.0D0/RAC(1))
      ENDDO
C
C     FERMI MODEL NORMALISATION CONSTANT RHO0
      RHO0 = 1.0D0 + PAC(2) - 6.0D0*RAC(3)*PLG(3)
      RHO0 = 3.0D0/(RHO0*4.0D0*PI*(CFMI(IZ)**3))
C
C     FERMI NORMALISATION CHARGE DENSITY
      DO I=1,11
        PRC(I) = 4.0D0*PI*RHO0*(CFMI(IZ)**(I))/DFLOAT(I)
      ENDDO
C
C     FERMI MODEL RADIAL MOMENTS
      RF(-2) = PRC(1)*(1.0D0 - RAC(1)*PLG(1))
      RF(-1) = PRC(2)*(1.0D0 + PAC(2)/3.0D0 + 2.0D0*RAC(2)*PLG(2))
      RF( 0) = PRC(3)*(1.0D0 + PAC(2) - 6.0D0*RAC(3)*PLG(3))
      RF( 1) = PRC(4)*(1.0D0 + 2.0D0*PAC(2) + 7.0D0*PAC(4)/15.0D0
     &       + 24.0D0*RAC(4)*PLG(4))
      RF( 2) = PRC(5)*(1.0D0 + 10.0D0*PAC(2)/3.0D0
     &       + 7.0D0*PAC(4)/3.0D0 - 120.0D0*RAC(5)*PLG(5))
      RF( 3) = PRC(6)*(1.0D0 + 5.0D0*PAC(2) + 7.0D0*PAC(4)
     &       + 31.0D0*PAC(6)/21.0D0 + 720.0D0*RAC(6)*PLG(6))
      RF( 4) = PRC(7)*(1.0D0 + 7.0D0*PAC(2) + 49.0D0*PAC(4)/3.0D0
     &       + 31.0D0*PAC(6)/3.0D0 - 5040.0D0*RAC(7)*PLG(7))
      RF( 5) = PRC(8)*(1.0D0 + 28.0D0*PAC(2) + 98.0D0*PAC(4)/3.0D0
     &       + 124.0D0*PAC(6)/3.0D0 + 127.0D0*PAC(8)/15.0D0 
     &       + 40320.0D0*RAC(8)*PLG(8))
      RF( 6) = PRC(9)*(1.0D0 + 12.0D0*PAC(2) + 294.0D0*PAC(4)/5.0D0
     &       + 124.0D0*PAC(6)/3.0D0 + 381.0D0*PAC(8)
     &       - 362880.0D0*RAC(9)*PLG(9))
      RF( 7) = PRC(10)*(1.0D0 + 150.0D0*PAC(2) + 98.0D0*PAC(4)
     &       + 310.0D0*PAC(6) + 381.0D0*PAC(8)  + 2555.0D0*PAC(10)
     &       + 3628800.0D0*RAC(10)*PLG(10))
      RF( 8) = PRC(11)*(1.0D0 + 55.0D0*PAC(2) + 154.0D0*PAC(4)
     &       + 682.0D0*PAC(6) + 1397.0D0*PAC(8) + 2555.0D0*PAC(10)
     &       - 39916800.0D0*RAC(11)*PLG(11))
C
C     SET UP MATRIX EQUATIONS FOR FERMI DISTRIBUTION
      IF(NMDL(IZ).EQ.'FERMI') THEN
        Y(1) = 0.0D0
        Y(2) = C4*(RF(4)-R0(4))
        Y(3) = C6*(RF(6)-R0(6))
C
C       POTENTIAL IN NUCLEAR REGION
        DO K=4,NFT
          R    = RM(K)
          Y(K) = VFERMI(R,CFMI(IZ),TFMI(IZ))
     &         - VGAUSS(R,XNUC(IZ,0),1.0D0)
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C     UNIFORM NUCLEAR MODEL DETAILS                                    C
C**********************************************************************C
C
C     POTENTIAL ON THE GRID FOR UNIFORM MODEL
      DO IPTS=0,NPTS
        R = 5.0D0*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
        VU(IPTS) = VUNIFM(R,RNUC(IZ))
      ENDDO
C
C     MOMENTS FOR UNIFORMLY-CHARGED NUCLEUS
      DO I=-2,8
        FC1 = 3.0D0
        FC2 = (5.0D0/3.0D0)**(0.5D0*I)
        FC3 = RNUC(IZ)**I
        FC5 = 3.0D0+I
        RU(I) = FC1*FC2*FC3/FC5
      ENDDO
C
C     SET UP MATRIX EQUATIONS FOR FERMI DISTRIBUTION
      IF(NMDL(IZ).EQ.'UNIFM') THEN
        Y(1) = 0.0D0
        Y(2) = C4*(RU(4)-R0(4))
        Y(3) = C6*(RU(6)-R0(6))
C
C       POTENTIAL IN NUCLEAR REGION
        DO K=4,NFT
          R    = RM(K)
          Y(K) = VUNIFM(R,RNUC(IZ))
     &         - VGAUSS(R,XNUC(IZ,0),1.0D0)
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C     GIVEN GEOMETRIC BETA, FIND OPTIMAL ALPHA VALUE.                  C
C**********************************************************************C
C
C     NUMBER OF INCREMENTS IN SEARCH
      NAF = 1200
C
C     INITIALISE BEST-FIT R-SQUARED AND INDEX FOR ALPHA
      RSQBIG = 0.0D0
      IABG   = 0
C
C     ITERATE OVER ALPHA SEARCH
      DO IA=0,NAF
C
C       GEOMETRIC SEED VALUES
        AF = ALPH*(1.0D0 + 3.0D0*DFLOAT(IA)/DFLOAT(NAF))/(RNUC(IZ)**2)
        BF = BETA
C
C       STORE GEOMETRIC SET OF PARAMETERS IN XNUC
        XI = AF
        DO IFT=1,NFT
          XNUC(IZ,IFT) = XI
          RT = DFLOAT(IFT-1)/DFLOAT(NFT+1)
          XI = BF*XI*(1.0D0 + GAMA*RT*RT)
        ENDDO
C
C       RESET FERMI Z MATRIX
        DO IFT=1,NFT
          Z(IFT) = Y(IFT)
        ENDDO
C
C       EVEN RADIAL MOMENTS OF EACH BASIS FUNCTION
        DO JFT=1,NFT
          X12 = DSQRT(XNUC(IZ,JFT))
          X32 = X12*XNUC(IZ,JFT)
          X52 = X32*XNUC(IZ,JFT)
          X72 = X52*XNUC(IZ,JFT)
          X(1,JFT) = C2*(1.5D0*PI12/X32)
          X(2,JFT) = C4*(7.5D0*PI12/X52)
          X(3,JFT) = C6*(315.0D0*PI12/(8.0D0*X72))
        ENDDO
C
C       POTENTIAL IN NUCLEAR REGION
        DO K=4,NFT
          R = RM(K)
          DO JFT=1,NFT
            X(K,JFT) = DEXP(-XNUC(IZ,JFT)*R*R)
          ENDDO
        ENDDO
C
C       SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
        CALL DGESV(NFT,1,X,NFT,IPIV,Z,NFT,INFO)
C
C       TRANSFER THE X VALUES TO AN ARRAY
        DO IFT=1,NFT
          FNUC(IZ,IFT) = Z(IFT)
        ENDDO
C
C       BEST-FIT GAUSSIAN POTENTIAL
        DO IPTS=0,NPTS
          R = 5.0D0*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
          VG(IPTS) = VGAUSS(R,XNUC(IZ,0),1.0D0)
          DO IFT=1,NFT
            VG(IPTS) = VG(IPTS) + FNUC(IZ,IFT)*DEXP(-XNUC(IZ,IFT)*R*R)
          ENDDO
        ENDDO
C
C       AVERAGE EVALUATED POTENTIAL OVER ALL RADII
        YB = 0.0D0
        DO IPTS=0,NPTS
          YB = YB + VG(IPTS)
        ENDDO
        YB = YB/DFLOAT(NPTS+1)
C
C       BITS AND PIECES FOR THE R-SQUARED VALUE
        SRES = 0.0D0
        STOT = 0.0D0
        IF(NMDL(IZ).EQ.'FERMI') THEN
          DO IPTS=0,NPTS
            SRES = SRES + (VG(IPTS)-VF(IPTS))**2
            STOT = STOT + (VG(IPTS)-YB      )**2
          ENDDO
        ELSEIF(NMDL(IZ).EQ.'UNIFM') THEN
          DO IPTS=0,NPTS
            SRES = SRES + (VG(IPTS)-VU(IPTS))**2
            STOT = STOT + (VG(IPTS)-YB      )**2
          ENDDO
        ENDIF
C
C       DECIDE WHETHER THIS IS THE BEST FIT SO FAR
        RSQ = 1.0D0 - SRES/STOT
        IF(RSQ.GT.RSQBIG) THEN
          IABG = IA
          RSQBIG = RSQ
        ENDIF
C
      ENDDO
C
C**********************************************************************C
C     END OF BEST-FIT DETERMINATION -- NOW IMPLEMENT THE BEST ALPHA.   C
C**********************************************************************C
C
C     GEOMETRIC SEED VALUES
      AF = ALPH*(1.0D0 + 3.0D0*DFLOAT(IABG)/DFLOAT(NAF))
      BF = BETA
      GF = GAMA
C
C     SET OF GAUSSIAN EXPONENTS
      XI = AF/(RNUC(IZ)**2)
      DO IFT=1,NFT
        XNUC(IZ,IFT) = XI
        RT = DFLOAT(IFT-1)/DFLOAT(NFT+1)
        XI = BF*XI*(1.0D0 + GAMA*RT*RT)
      ENDDO
C
C     RESET FERMI Z MATRIX
      DO IFT=1,NFT
        Z(IFT) = Y(IFT)
      ENDDO
C
C     EVEN RADIAL MOMENTS OF EACH BASIS FUNCTION
      DO JFT=1,NFT
        X12 = DSQRT(XNUC(IZ,JFT))
        X32 = X12*XNUC(IZ,JFT)
        X52 = X32*XNUC(IZ,JFT)
        X72 = X52*XNUC(IZ,JFT)
        X(1,JFT) = C2*(1.5D0*PI12/X32)
        X(2,JFT) = C4*(7.5D0*PI12/X52)
        X(3,JFT) = C6*(315.0D0*PI12/(8.0D0*X72))
      ENDDO
C
C     POTENTIAL IN NUCLEAR REGION
      DO K=4,NFT
        R = RM(K)
        DO JFT=1,NFT
          X(K,JFT) = DEXP(-XNUC(IZ,JFT)*R*R)
        ENDDO
      ENDDO
C
C     SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
      CALL DGESV(NFT,1,X,NFT,IPIV,Z,NFT,INFO)
C
C     TRANSFER THE X VALUES TO FRACTIONAL ARRAY
      DO IFT=1,NFT
        FNUC(IZ,IFT) = Z(IFT)
      ENDDO
C
C     QUIT HERE IF NO WRITTEN RESULTS ARE REQUIRED
      IF(IWRT.EQ.0) GOTO 20
C
C     WRITE THE BEST-FIT SOLUTION
35    FORMAT(1X,A,I2)
36    FORMAT(1X,A,F10.8,1X,A)
37    FORMAT(1X,I2,2X,F19.8,2X,F11.8)
      WRITE(6, *) ''
      WRITE(7, *) ''
      WRITE(6,35) 'Nuclear gaussian basis set for IZ = ',IZ
      WRITE(7,35) 'Nuclear gaussian basis set for IZ = ',IZ
      WRITE(6,36) 'Here ξ0 = 3/2R_n^2 and R_n  = ',RNUC(IZ),'a0'
      WRITE(7,36) 'Here ξ0 = 3/2R_n^2 and R_n  = ',RNUC(IZ),'a0'
      WRITE(6,36) 'Least-squares best fit: R^2 = ',RSQBIG
      WRITE(7,36) 'Least-squares best fit: R^2 = ',RSQBIG
      WRITE(6, *) REPEAT('=',36)
      WRITE(7, *) REPEAT('=',36)
      WRITE(6,36) ' i                A_{i}    ξ_{i}*R^2'
      WRITE(7,36) ' i                A_{i}    ξ_{i}*R^2'
      WRITE(6, *) REPEAT('-',36)
      WRITE(7, *) REPEAT('-',36)
      WRITE(6,37) 0,1.0D0,1.5D0
      WRITE(7,37) 0,1.0D0,1.5D0
      WRITE(6, *) REPEAT('-',36)
      WRITE(7, *) REPEAT('-',36)
      DO IFT=1,NFT
        WRITE(6,37) IFT,FNUC(IZ,IFT),RNUC(IZ)*RNUC(IZ)*XNUC(IZ,IFT)
        WRITE(7,37) IFT,FNUC(IZ,IFT),RNUC(IZ)*RNUC(IZ)*XNUC(IZ,IFT)
      ENDDO
      WRITE(6, *) REPEAT('=',36)
      WRITE(7, *) REPEAT('=',36)
C
C**********************************************************************C
C     RADIAL MOMENTS FOR THIS AND OTHER MODELS.                        C
C**********************************************************************C
C
C     FINITE SET OF GAUSSIANS
      DO I=-2,8
        FC1 = 2.0D0/PI12
        FC2 = GAMHLF(3+I)
        FC3 = XNUC(IZ,0)**(0.5D0*I)
        RG(I) = FC1*FC2/FC3
        DO IFT=1,NFT
          XI  = XNUC(IZ,IFT)
          FC  = FNUC(IZ,IFT)
          FC1 = FC*DFLOAT(I)
          FC2 = GAMHLF(3+I)
          FC3 = XI**(I+1)
          RG(I) = RG(I) + FC1*FC2/DSQRT(FC3)
        ENDDO
      ENDDO
C
C     CREATE TABLE FOR THE MOMENTS
49    FORMAT(2X,A,7X,A,12X,A,12X,A,12X,A)
52    FORMAT(2X,'<r^',I2,'>^1/',I2,4(2X,ES16.8))
53    FORMAT(2X,'<r^',I2,'>',5X,4(2X,ES16.8))
      WRITE(6, *) ''
      WRITE(7, *) ''
      WRITE(6,49) 'Moment (fm or e)','Gauss0','Unifrm',' Fermi','Bestft'
      WRITE(7,49) 'Moment (fm or e)','Gauss0','Unifrm',' Fermi','Bestft'
      WRITE(6, *) REPEAT('-',84)
      WRITE(7, *) REPEAT('-',84)
      CF = 1.0D0/(CFM**2)
      DO I=-2,8
        IF(I.NE.0) THEN
          RU(I) = (RU(I)*CF)**(1.0D0/DFLOAT(I))
          RF(I) = (RF(I)*CF)**(1.0D0/DFLOAT(I))
          RG(I) = (RG(I)*CF)**(1.0D0/DFLOAT(I))
          R0(I) = (R0(I)*CF)**(1.0D0/DFLOAT(I))
          WRITE(6,52) I,I,R0(I),RU(I),RF(I),RG(I)
          WRITE(7,52) I,I,R0(I),RU(I),RF(I),RG(I)
        ELSE
          RU(I) = RU(I)*CF
          RF(I) = RF(I)*CF
          RG(I) = RG(I)*CF
          R0(I) = R0(I)*CF
          WRITE(6,53) I,  R0(I),RU(I),RF(I),RG(I)
          WRITE(7,53) I,  R0(I),RU(I),RF(I),RG(I)
        ENDIF
        CF = CF*CFM
      ENDDO
      WRITE(6, *) REPEAT('-',84)
      WRITE(7, *) REPEAT('-',84)
      WRITE(6, *) ''
      WRITE(7, *) ''
C
C**********************************************************************C
C     SCALAR POTENTIALS FOR THIS AND OTHER MODELS.                     C
C**********************************************************************C
C
C     WRITE ALL THIS INTO A BIG ARRAY
      DO IPTS=0,NPTS
        RADS(IPTS) = 5.0D0*DFLOAT(IPTS)/DFLOAT(NPTS)
        VNUC(1,IPTS) = VU(IPTS)
        VNUC(2,IPTS) = VF(IPTS)
        VNUC(3,IPTS) = V0(IPTS)
        VNUC(4,IPTS) = VG(IPTS)
      ENDDO
C
C     CREATE A TABLE FOR THE VNUC
C     GOTO 23
50    FORMAT(1X,F8.4,'R',3X,F20.10,2X,F20.10,2X,F20.10)
51    FORMAT(2X,A,13X,A,11X,A,11X,A)
      WRITE(6,51) 'Radius r','V_Gauss0(r)','V_Fermi(r)','V_GaussN(r)'
      WRITE(7,51) 'Radius r','V_Gauss0(r)','V_Fermi(r)','V_GaussN(r)'
      WRITE(6, *) REPEAT('-',76)
      WRITE(7, *) REPEAT('-',76)
      DO IPTS=0,NPTS,100
        R = 10.0D0*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
        WRITE(6,50) R/RNUC(IZ),V0(IPTS),VF(IPTS),VG(IPTS)
        WRITE(7,50) R/RNUC(IZ),V0(IPTS),VF(IPTS),VG(IPTS)
      ENDDO
      WRITE(6,50) 15.0D0,
     &            VGAUSS(15.0D0*RNUC(IZ),XNUC(IZ,0),1.0D0),
     &            VFERMI(15.0D0*RNUC(IZ),CFMI(IZ),TFMI(IZ)),
     &            VFITGS(15.0D0*RNUC(IZ),IZ)
      WRITE(6,50) 20.0D0,
     &            VGAUSS(20.0D0*RNUC(IZ),XNUC(IZ,0),1.0D0),
     &            VFERMI(20.0D0*RNUC(IZ),CFMI(IZ),TFMI(IZ)),
     &            VFITGS(20.0D0*RNUC(IZ),IZ)
      WRITE(6, *) REPEAT('-',76)
      WRITE(7, *) REPEAT('-',76)
23    CONTINUE
C
C     FILE NAME AND TITLES
      XOUT   = 'Potentials'
      TITLE  = 'Potentials'
      XAXIS  = 'r/RNUC'
      YAXIS  = '-V(r)/Z'
      KEY(1) = 'Uniform'
      KEY(2) = 'Fermi'
      KEY(3) = 'Gaussian'
      KEY(4) = 'Fitting'
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C       BEGIN LOOP OVER ALL POINTS
        DO IPTS=0,NPTS
C
C         PRINT DISTANCE ALONG PATH AND FOUR-CURRENT
          WRITE(8, *) RADS(IPTS),(-VNUC(ITP,IPTS),ITP=1,4)
C
        ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,4,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C**********************************************************************C
C     RADIAL CHARGE DENSITIES FOR THIS AND OTHER MODELS.               C
C**********************************************************************C
C
C     CHARGE DENSITIES FOR ALL EXACT NUCLEAR MODELS
      RAT = XNUC(IZ,0)/PI
      DO IPTS=0,NPTS
C
C       RADIAL VALUE
        R = 5.0D0*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C       CHARGE DENSITY AT THIS RADIUS
        RADS(IPTS) = R/RNUC(IZ)
        PNUC(1,IPTS) = 0.0D0
        PNUC(3,IPTS) = 0.0D0
        PNUC(4,IPTS) = 0.0D0
        PNUC(1,IPTS) = RHONUC('UNIFM',IZ,R)
        PNUC(2,IPTS) = RHONUC('FERMI',IZ,R)
        PNUC(3,IPTS) = RHONUC('GAUSS',IZ,R)
        PNUC(4,IPTS) = RAT*DSQRT(RAT)*DEXP(-XNUC(IZ,0)*R*R)
        DO IFT=1,NFT
          XI   = XNUC(IZ,IFT)
          RHOI = XI*(3.0D0-2.0D0*XI*R*R)*DEXP(-XI*R*R)
          PNUC(4,IPTS) = PNUC(4,IPTS) - FNUC(IZ,IFT)*RHOI/(2.0D0*PI)
        ENDDO
      ENDDO
C
C     FILE NAME AND TITLES
      XOUT   = 'Densities'
      TITLE  = 'Densities'
      XAXIS  = 'r/RNUC'
      YAXIS  = '-rho(r)/Z'
      KEY(1) = 'Uniform'
      KEY(2) = 'Fermi'
      KEY(3) = 'Gaussian'
      KEY(4) = 'Fitting'
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C       BEGIN LOOP OVER ALL POINTS
        DO IPTS=0,NPTS
C
C         PRINT DISTANCE ALONG PATH AND FOUR-CURRENT
          WRITE(8, *) RADS(IPTS),(PNUC(ITP,IPTS),ITP=1,4)
C
        ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,4,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C**********************************************************************C
C     WEIGHTED RADIAL CHARGE DENSITIES FOR THIS AND OTHER MODELS.      C
C**********************************************************************C
C
C     FILE NAME AND TITLES
      XOUT   = 'WeightedDensities'
      TITLE  = 'WeightedDensities'
      XAXIS  = 'r/RNUC'
      YAXIS  = '4*pi*r*r*rho(r)/Z'
      KEY(1) = 'Uniform'
      KEY(2) = 'Fermi'
      KEY(3) = 'Gaussian'
      KEY(4) = 'Fitting'
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C       BEGIN LOOP OVER ALL POINTS
        DO IPTS=0,NPTS
C
C         PRINT DISTANCE ALONG PATH AND FOUR-CURRENT
          R = RADS(IPTS)
          W = R*(RNUC(IZ)*RNUC(IZ))
          WRITE(8, *) R,(4.0D0*PI*W*W*PNUC(ITP,IPTS),ITP=1,4)
C
        ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,4,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C**********************************************************************C
C     RADIAL ELECTRIC FIELD FOR THIS AND OTHER MODELS.                 C
C**********************************************************************C
C
C     CHARGE DENSITIES FOR ALL EXACT NUCLEAR MODELS
      DO IPTS=0,NPTS
C
C       RADIAL VALUE
        R = 5.0D0*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C       CHARGE DENSITY AT THIS RADIUS
        RADS(IPTS) = R/RNUC(IZ)
        ENUC(1,IPTS) = ERNUC('UNIFM',IZ,R)
        ENUC(2,IPTS) = ERNUC('FERMI',IZ,R)
        ENUC(3,IPTS) = ERNUC('GAUSS',IZ,R)
        ENUC(4,IPTS) = 0.0D0
        DO IFT=1,NFT
          XI = XNUC(IZ,IFT)
          X1 = DSQRT(XI)*R
          IF(R.LT.0.2D0*RNUC(IZ)) THEN
            X3 = XI*R*R*X1
            X5 = XI*R*R*X3
            X7 = XI*R*R*X5
            X9 = XI*R*R*X7
            ERRI = 4.0D0*X1/3.0D0 - 4.0D0*X3/5.0D0 + 2.0D0*X5/7.0D0
     &            - 2.0D0*X7/27.0D0 + X9/66.0D0
            ERRI = XI*ERRI/PI12
          ELSE
            ERRI =-2.0D0*X1*DEXP(-X1*X1)/PI12 + DERF(X1)
            ERRI = ERRI/(R*R)
          ENDIF
          ENUC(4,IPTS) = ENUC(4,IPTS) + FNUC(IZ,IFT)*ERRI
        ENDDO
C
      ENDDO
C
C     FILE NAME AND TITLES
      XOUT   = 'RadialEfield'
      TITLE  = 'Radial electric field component E_{r}(r)'
      XAXIS  = 'r/RNUC'
      YAXIS  = 'E_{r}(r)/Z'
      KEY(1) = 'Uniform'
      KEY(2) = 'Fermi'
      KEY(3) = 'Gaussian'
      KEY(4) = 'Fitting'
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C       BEGIN LOOP OVER ALL POINTS
        DO IPTS=0,NPTS
C
C         PRINT DISTANCE ALONG PATH AND FOUR-CURRENT
          WRITE(8, *) RADS(IPTS),(ABS(ENUC(ITP,IPTS)),ITP=1,4)
C
        ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
cC
cC     GENERATE A GNUPLOT MAKE FILE
c      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,4,KEY)
cC
cC     EXECUTE GNUPLOT COMMAND IN TERMINAL
c      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
c      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C**********************************************************************C
C     BASIS FUNCTION POTENTIALS (PRE-NORMALISED AT THE ORIGIN)         C
C**********************************************************************C
C
C     FILE NAME AND TITLES
      XOUT   = 'Basis-potentials'
      TITLE  = 'Basis-potenitals'
      XAXIS  = 'r/RNUC'
      YAXIS  = '-V(r)'
      DO IFT=1,NFT
        KEYNFT(IFT) = 'X'
      ENDDO
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C       BEGIN LOOP OVER ALL POINTS
        DO IPTS=0,NPTS
C
C         PRINT DISTANCE ALONG PATH AND FOUR-CURRENT
          R = RADS(IPTS)/CFM
          WRITE(8, *) RADS(IPTS),(DEXP(-XNUC(IZ,IFT)*R*R),IFT=1,NFT)
C
        ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,NFT,KEYNFT)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C**********************************************************************C
C     BASIS FUNCTION DENSITIES                                         C
C**********************************************************************C
C
C     FILE NAME AND TITLES
      XOUT   = 'Basis-densities'
      TITLE  = 'Basis-densities'
      XAXIS  = 'r/RNUC'
      YAXIS  = 'rho(r)'
      DO IFT=1,NFT
        KEYNFT(IFT) = 'X'
      ENDDO
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C       BEGIN LOOP OVER ALL POINTS
        DO IPTS=0,NPTS
C
C         PRINT DISTANCE ALONG PATH AND FOUR-CURRENT
          R = RADS(IPTS)/CFM
          DO IFT=1,NFT
            PM(IFT) = (3.0D0-2.0D0*XNUC(IZ,IFT)*R*R)
            PM(IFT) = XNUC(IZ,IFT)*PM(IFT)*DEXP(-XNUC(IZ,IFT)*R*R)
          ENDDO
          WRITE(8, *) RADS(IPTS),(R*R*PM(IFT),IFT=1,NFT)
C
        ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,NFT,KEYNFT)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
20    CONTINUE
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
C  ▶ NFT > 26 (FOR WHICHEVER NUCLEAR CHARGE MODEL.)                    C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*2 ELMT(120)
      CHARACTER*5 NMDL
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*60 VACFL
C
      DIMENSION RMS(MCT,5),RZR(MCT,5),QPL(MCT,5),NFT(MCT,5)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BUEE/RUEE(MCT,3),FUEE(MCT,MFT),XUEE(MCT,MFT),NUEE(MCT)
      COMMON/BUEU/RUEU(MCT,3),FUEU(MCT,MFT),XUEU(MCT,MFT),NUEU(MCT)
      COMMON/BUEP/RUEP(MCT,3),FUEP(MCT,MFT),XUEP(MCT,MFT),NUEP(MCT)
      COMMON/BKSB/RKSB(MCT,3),FKSB(MCT,MFT),XKSB(MCT,MFT),NKSB(MCT)
      COMMON/BWKR/RWKR(MCT,3),FWKR(MCT,MFT),XWKR(MCT,MFT),NWKR(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/MDLV/ELMT
      COMMON/QUEE/RE(0:NRAD),VE(MCT,0:NRAD),R0E,RME,RIE,NUE,NEE
      COMMON/QUEU/RU(0:NRAD),VU(MCT,0:NRAD),R0U,RMU,RIU,NUU,NEU
      COMMON/QUEP/RP(0:NRAD),VP(MCT,0:NRAD),R0P,RMP,RIP,MUP,NEP
      COMMON/QWKR/RW(0:NRAD),VW(MCT,0:NRAD),R0W,RMW,RIW,NUW,NEW
      COMMON/QKSB/RK(0:NRAD),VK(MCT,0:NRAD),R0K,RMK,RIK,NUK,NEK
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/BMSS/RMPC,RMP0,RMKC,RMK0
C
C     TOGGLE FOR MORE DETAILED SUMMARIES (0 FOR OFF, 1 FOR ON)
      IWRT = 0
C
C     EARLY EXIT IF VACPOL MEAN FIELDS AREN'T NEEDED
      IF(HMLT.NE.'DHFQ'.AND.HMLT.NE.'DHFP') RETURN
C
C     VIRTUAL FIELD MASS MODIFIER (1.0D0 FOR ELECTRON-POSITRON FIELD)
      RV = EMSS
C
C     NUMBER OF FITTING FUNCTIONS FOR ALL PROCEDURES
      DO IZ=1,NCNT
        DO K=1,5
          NFT(IZ,K) = 26
        ENDDO
        NUEE(IZ) = 0
        NUEU(IZ) = 0
        NUEP(IZ) = 0
        NWKR(IZ) = 0
        NKSB(IZ) = 0
      ENDDO
C
C     EMPTY SOME QUANTITATIVE COUNTERS
      DO IZ=1,NCNT
        DO K=1,5
          RMS(IZ,K) = 0.0D0
          RZR(IZ,K) = 0.0D0
          QPL(IZ,K) = 0.0D0
        ENDDO
      ENDDO
C
C     INITIALISE INTEGRATION GRID
      CALL GENGRID
C
C**********************************************************************C
C     PERFORM FITTING PROCEDURE FOR ALL VP POTENTIALS.                 C
C**********************************************************************C
C
C     WRITE RESULTS OF VACPOL BEST-FIT PROCEDURE
31    FORMAT(33X,A)
32    FORMAT(1X,A,2X,'|',1X,A,7X,A,3X,'| ',
     &                         A,3X,A,3X,A,' |   ',A,2X,A,' | ',A,2X,A)
33    FORMAT(1X,I2,' (',A,')',1X,'|',1X,A,2X,A,' |',1X,F7.4,2X,F8.4,
     &        2X,F8.4,' | ',F8.5,2X,F7.4,' | ',I2,4X,F10.8)
34    FORMAT(9X,'|',1X,A,2X,A,' |',1X,F7.4,2X,F8.4,
     &        2X,F8.4,' | ',F8.5,2X,F7.4,' | ',I2,4X,F10.8)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('=',102)
      WRITE(7, *) REPEAT('=',102)
      WRITE(6,31) 'Nuclear vacuum polarisation summary'
      WRITE(7,31) 'Nuclear vacuum polarisation summary'
      WRITE(6, *) REPEAT('=',102)
      WRITE(7, *) REPEAT('=',102)
      WRITE(6,32) 'Centre','Potential','Order',
     &            'R0 (fm)','R2 (fm)','R4 (fm)',
     &            'Q± (e)','RZ (fm)','N_fit','R-squared'
      WRITE(7,32) 'Centre','Potential','Order',
     &            'R0 (fm)','R2 (fm)','R4 (fm)',
     &            'Q± (e)','RZ (fm)','N_fit','R-squared'
      WRITE(6, *) REPEAT('-',102)
      WRITE(7, *) REPEAT('-',102)
      DO IZ=1,NCNT
C
C       ELEMENT LABEL
        IZNC = INT(ZNUC(IZ))
C
c       GOTO 100
C       GENERATE UEHLING POTENTIAL (e+e-)
        CALL VUEHGEN(IZ,NFT(IZ,1),IWRT,RV,RMS(IZ,1),RZR(IZ,1),
     &   QPL(IZ,1),'ELEC',RUEE,FUEE,XUEE,NUEE,RE,VE,R0E,RME,RIE,NUE,NEE)
        QNET = RUEE(IZ,1)
        QRMS = DSQRT(RUEE(IZ,2))
        QRMQ = DSQRT(DSQRT(RUEE(IZ,3)))
        WRITE(6,33) IZ,ELMT(IZNC),
     &              'Uehling (e+e-)','α (Zα) ',
     &              QNET*CFM,QRMS*CFM,QRMQ*CFM,
     &              ZNUC(IZ)*QPL(IZ,1),RZR(IZ,1)*CFM,
     &              NFT(IZ,1),RMS(IZ,1)
        WRITE(7,33) IZ,ELMT(IZNC),
     &              'Uehling (e+e-)','α (Zα) ',
     &              QNET*CFM,QRMS*CFM,QRMQ*CFM,
     &              ZNUC(IZ)*QPL(IZ,1),RZR(IZ,1)*CFM,
     &              NFT(IZ,1),RMS(IZ,1)
100     CONTINUE
C
c       GOTO 200
C       GENERATE WICHMANN-KROLL POTENTIAL
        CALL VWKRGEN(IZ,NFT(IZ,2),IWRT,RV,RMS(IZ,2),RZR(IZ,2),QPL(IZ,2))
        QNET = RWKR(IZ,1)
        QRMS =-DSQRT(-RWKR(IZ,2))
        QRMQ =-DSQRT(DSQRT(-RWKR(IZ,3)))
        WRITE(6,34) 'Wichmann-Kroll','α (Zα)³',
     &              QNET*CFM,QRMS*CFM,QRMQ*CFM,
     &             -ZNUC(IZ)*ZNUC(IZ)*ZNUC(IZ)*QPL(IZ,2),RZR(IZ,2)*CFM,
     &              NFT(IZ,2),RMS(IZ,2)
        WRITE(7,34) 'Wichmann-Kroll','α (Zα)³',
     &              QNET*CFM,QRMS*CFM,QRMQ*CFM,
     &             -ZNUC(IZ)*ZNUC(IZ)*ZNUC(IZ)*QPL(IZ,2),RZR(IZ,2)*CFM,
     &              NFT(IZ,2),RMS(IZ,2)
200     CONTINUE
C
c       GOTO 300
C       GENERATE KÄLLÉN-SABRY POTENTIAL
        CALL KSAUXFN
        CALL VKSBGEN(IZ,NFT(IZ,3),IWRT,RV,RMS(IZ,3),RZR(IZ,3),QPL(IZ,3))
        QNET = RKSB(IZ,1)
        QRMS = DSQRT(RKSB(IZ,2))
        QRMQ = DSQRT(DSQRT(RKSB(IZ,3)))
        WRITE(6,34) 'Källén-Sabry  ','α²(Zα) ',
     &              QNET*CFM,QRMS*CFM,QRMQ*CFM,
     &              ZNUC(IZ)*QPL(IZ,3),RZR(IZ,3)*CFM,
     &              NFT(IZ,3),RMS(IZ,3)
        WRITE(7,34) 'Källén-Sabry  ','α²(Zα) ',
     &              QNET*CFM,QRMS*CFM,QRMQ*CFM,
     &              ZNUC(IZ)*QPL(IZ,3),RZR(IZ,3)*CFM,
     &              NFT(IZ,3),RMS(IZ,3)
300     CONTINUE
C
c       GOTO 400
C       GENERATE UEHLING POTENTIAL (μ+μ-)
        CALL VUEHGEN(IZ,NFT(IZ,4),IWRT,UMSS,RMS(IZ,4),RZR(IZ,4),
     &   QPL(IZ,4),'MUON',RUEU,FUEU,XUEU,NUEU,RU,VU,R0U,RMU,RIU,NUU,NEU)
        QNET = RUEU(IZ,1)
        QRMS = DSQRT(RUEU(IZ,2))
        QRMQ = DSQRT(DSQRT(RUEU(IZ,3)))
        WRITE(6,34) 'Uehling (μ+μ-)','α (Zα) ',
     &              QNET*CFM,QRMS*CFM,QRMQ*CFM,
     &              ZNUC(IZ)*QPL(IZ,4),RZR(IZ,4)*CFM,
     &              NFT(IZ,2),RMS(IZ,4)
        WRITE(7,34) 'Uehling (μ+μ-)','α (Zα) ',
     &              QNET*CFM,QRMS*CFM,QRMQ*CFM,
     &              ZNUC(IZ)*QPL(IZ,4),RZR(IZ,4)*CFM,
     &              NFT(IZ,4),RMS(IZ,4)
400     CONTINUE
C
C       GOTO 500
C       GENERATE UEHLING POTENTIAL (π+π-)
        CALL VUEHGEN(IZ,NFT(IZ,5),IWRT,RMPC,RMS(IZ,5),RZR(IZ,5),
     &   QPL(IZ,5),'PION',RUEP,FUEP,XUEP,NUEP,RP,VP,R0P,RMP,RIP,MUP,NEP)
        QNET = RUEP(IZ,1)
        QRMS = DSQRT(RUEP(IZ,2))
        QRMQ = DSQRT(DSQRT(RUEP(IZ,3)))
        WRITE(6,34) 'Uehling (π+π-)','α (Zα) ',
     &              QNET*CFM,QRMS*CFM,QRMQ*CFM,
     &              ZNUC(IZ)*QPL(IZ,5),RZR(IZ,5)*CFM,
     &              NFT(IZ,5),RMS(IZ,5)
        WRITE(7,34) 'Uehling (π+π-)','α (Zα) ',
     &              QNET*CFM,QRMS*CFM,QRMQ*CFM,
     &              ZNUC(IZ)*QPL(IZ,5),RZR(IZ,5)*CFM,
     &              NFT(IZ,5),RMS(IZ,5)
500     CONTINUE
C
        IF(IZ.NE.NCNT) THEN
          WRITE(6, *) REPEAT('-',102)
        ENDIF
C
      ENDDO
      WRITE(6, *) REPEAT('=',102)
      WRITE(7, *) REPEAT('=',102)
      WRITE(6, *) ' '
      WRITE(7, *) ' '
C
C**********************************************************************C
C     SAVE THE EFFECTIVE POTENTIALS TO EXTERNAL FILES AS NEEDED.       C
C**********************************************************************C
C
C     SAVE UEHLING FITTING SET (e+e-)
C     GOTO 101
      VACFL = 'output/'//TRIM(MOLCL)//'_vpuee.dat'
      OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
        DO IZ=1,NCNT
          WRITE(8, *) NUEE(IZ)
          DO IFT=1,NUEE(IZ)
            WRITE(8, *) FUEE(IZ,IFT),XUEE(IZ,IFT)
          ENDDO
        ENDDO
      CLOSE(UNIT=8)
C
C     SAVE UEHLING RADIAL GRID DATA (e+e-)
      IF(IWRT.EQ.1) THEN
        VACFL = 'output/'//TRIM(MOLCL)//'_vpuee_grid.dat'
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          WRITE(8, *) R0E,RME,RIE,NUE,NEE
          DO N=0,NRAD
            WRITE(8, *) RE(N),(VE(IZ,N),IZ=1,NCNT)
          ENDDO
        CLOSE(UNIT=8)
      ENDIF
101   CONTINUE
C
C     SAVE WICHMANN-KROLL FITTING SET
C     GOTO 201
      VACFL = 'output/'//TRIM(MOLCL)//'_vpwkr.dat'
      OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
        DO IZ=1,NCNT
          WRITE(8, *) NWKR(IZ)
          DO IFT=1,NWKR(IZ)
            WRITE(8, *) FWKR(IZ,IFT),XWKR(IZ,IFT)
          ENDDO
        ENDDO
      CLOSE(UNIT=8)
C
C     SAVE WICHMANN-KROLL RADIAL GRID DATA
      IF(IWRT.EQ.1) THEN
        VACFL = 'output/'//TRIM(MOLCL)//'_vpwkr_grid.dat'
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          WRITE(8, *) R0W,RMW,RIW,NUW,NEW
          DO N=0,NRAD
            WRITE(8, *) RW(N),(VW(IZ,N),IZ=1,NCNT)
          ENDDO
        CLOSE(UNIT=8)
      ENDIF
201   CONTINUE
C
C     SAVE KÄLLÉN-SABRY FITTING SET
C     GOTO 301
      VACFL = 'output/'//TRIM(MOLCL)//'_vpksb.dat'
      OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
        DO IZ=1,NCNT
          WRITE(8, *) NKSB(IZ)
          DO IFT=1,NKSB(IZ)
            WRITE(8, *) FKSB(IZ,IFT),XKSB(IZ,IFT)
          ENDDO
        ENDDO
      CLOSE(UNIT=8)
C
C     SAVE KÄLLÉN-SABRY RADIAL GRID DATA
      IF(IWRT.EQ.1) THEN
        VACFL = 'output/'//TRIM(MOLCL)//'_vpksb_grid.dat'
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          WRITE(8, *) R0K,RMK,RIK,NUK,NEK
          DO N=0,NRAD
            WRITE(8, *) RK(N),(VK(IZ,N),IZ=1,NCNT)
          ENDDO
        CLOSE(UNIT=8)
      ENDIF
301   CONTINUE
C
C     SAVE UEHLING FITTING SET (μ+μ-)
C     GOTO 401
      VACFL = 'output/'//TRIM(MOLCL)//'_vpueu.dat'
      OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
        DO IZ=1,NCNT
          WRITE(8, *) NUEU(IZ)
          DO IFT=1,NUEU(IZ)
            WRITE(8, *) FUEU(IZ,IFT),XUEU(IZ,IFT)
          ENDDO
        ENDDO
      CLOSE(UNIT=8)
C
C     SAVE UEHLING RADIAL GRID DATA (μ+μ-)
      IF(IWRT.EQ.1) THEN
        VACFL = 'output/'//TRIM(MOLCL)//'_vpueu_grid.dat'
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          WRITE(8, *) R0U,RMU,RIU,NUU,NEU
          DO N=0,NRAD
            WRITE(8, *) RU(N),(VU(IZ,N),IZ=1,NCNT)
          ENDDO
        CLOSE(UNIT=8)
      ENDIF
401   CONTINUE
C
C     SAVE UEHLING FITTING SET (π+π-)
C     GOTO 501
      VACFL = 'output/'//TRIM(MOLCL)//'_vpuep.dat'
      OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
        DO IZ=1,NCNT
          WRITE(8, *) NUEP(IZ)
          DO IFT=1,NUEP(IZ)
            WRITE(8, *) FUEP(IZ,IFT),XUEP(IZ,IFT)
          ENDDO
        ENDDO
      CLOSE(UNIT=8)
C
C     SAVE UEHLING RADIAL GRID DATA (π+π-)
      IF(IWRT.EQ.1) THEN
        VACFL = 'output/'//TRIM(MOLCL)//'_vpuep_grid.dat'
        OPEN(UNIT=8,FILE=TRIM(VACFL),STATUS='UNKNOWN')
          WRITE(8, *) R0P,RMP,RIP,MUP,NEP
          DO N=0,NRAD
            WRITE(8, *) RP(N),(VP(IZ,N),IZ=1,NCNT)
          ENDDO
        CLOSE(UNIT=8)
      ENDIF
501   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE VUEHGEN(IZ,NFT,IWRT,RV,RSQBIG,RZER,QPOL,COUPLING,
     &           RUEH,FUEH,XUEH,NUEH,RAD,VVAC,RORI,RMID,RMAX,NLIN,NEXP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    VV    VV UU    UU EEEEEEEE HH    HH  GGGGGG  EEEEEEEE NN    NN    C
C    VV    VV UU    UU EE       HH    HH GG    GG EE       NNN   NN    C
C    VV    VV UU    UU EE       HH    HH GG       EE       NNNN  NN    C
C    VV    VV UU    UU EEEEEE   HHHHHHHH GG       EEEEEE   NN NN NN    C
C     VV  VV  UU    UU EE       HH    HH GG   GGG EE       NN  NNNN    C
C      VVVV   UU    UU EE       HH    HH GG    GG EE       NN   NNN    C
C       VV     UUUUUU  EEEEEEEE HH    HH  GGGGGG  EEEEEEEE NN    NN    C
C                                                                      C
C -------------------------------------------------------------------- C
C  VUEHGEN GENERATES A BEST-FIT GAUSSIAN SET FOR THE UEHLING POTENTIAL C
C  ARISING FROM NUCLEUS IZ, USING NFT TOTAL GAUSSIANS.                 C
C -------------------------------------------------------------------- C
C  MATCHING CRITERIA AVAILABLE FOR GAUSSIAN AMPLITUDES:                C
C   ▶ RADIAL ARGUMENTS SCALED BY NUCLEAR RMS RADIUS RNUC(IZ).          C
C   ▶ UEHLING POTENTIAL FIRST SCALED TO R*R*V(R).                      C
C   ▶ ZEROTH MOMENT SATISFIED AUTOMATICALLY BUT SECOND MOMENT ENFORCED.C
C   ▶ POTENTIAL MATCHED TO NFT POINTS, FROM 3.0D0*RNUC TO 200.0D0*RNUC.C
C   ▶ MATCHING POINTS ARE EXPONENTIALLY SPACED.                        C
C   ▶ GAUSSIAN EXPONENTS FROM MODIFIED GEOMETRIC SERIES WITH THE SAME  C
C     BETA AND GAMMA, BUT ALPHA IS SCALED BY RNUC*RNUC.                C
C   ▶ WHEN THIS IS COMPLETE, POTENTIAL V(R) NEAR R=0.0D0 IS AUGMENTED  C
C     WITH SOME MORE GAUSSIANS WITH LARGER EXPONENTS IF V(R) IS NEEDED.C
C -------------------------------------------------------------------- C
C INPUT:                                                               C
C   ▶ IZ:   NUCLEAR CENTRE.                                            C
C   ▶ NFT:  DESIRED NUMBER OF FITTING FUNCTIONS IN THE SET.            C
C   ▶ IWRT: EXTENDED WRITTEN OUTPUT OPTION.                            C
C   ▶ RV:   MASS MODIFIER (1.0D0 FOR ELECTRONS, MASS RATIO OTHERWISE). C
C**********************************************************************C
      INCLUDE 'parameters.h'
      PARAMETER(NPTS=6000)
C
      EXTERNAL CHIUEHF,CHIUEHB
C
      LOGICAL FILETHERE
C
      CHARACTER*2  ELMT(120)
      CHARACTER*3  ATRM
      CHARACTER*4  COUPLING
      CHARACTER*5  NMDL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(2)
C
      DIMENSION RUEH(MCT,3),FUEH(MCT,MFT),XUEH(MCT,MFT),NUEH(MCT)
      DIMENSION RAD(0:NRAD),VVAC(MCT,0:NRAD)
      DIMENSION XS(0:NLW),Y2S(0:NLW),D2S(0:NLW),Y1S(0:NLW),D1S(0:NLW),
     &          XB(0:NUP),Y2B(0:NUP),D2B(0:NUP)
      DIMENSION VF(0:NPTS),VG(0:NPTS),RU(0:NPTS)
      DIMENSION X(NFT,NFT),Z(NFT),Y(NFT),RM(NFT),IPIV(NFT)
      DIMENSION RHO(0:NSRC),RHOTMP(0:NRVP)
      DIMENSION POL(0:NPTS),RK(-2:5)
      DIMENSION FMQ(0:NQVP),PUEH(0:NQVP),GUEH(0:NQVP)
      DIMENSION VTMP(0:NRAD)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/MDLV/ELMT
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/QGRD/QFM(0:NQVP),QEP(0:NQVP),QMIN,QMID,QMAX
      COMMON/RGRD/RFM(0:NRVP),REP(0:NRVP),RMIN,RMD,RMX
C
C     MULTIPLIER FOR COMPTON WAVELENGTH (USE FOR MUON OR TAUON FIELD)
      CMPF = CMPW/RV
C
C     CHECK THAT THERE ARE ENOUGH FITTING FUNCTIONS
      IF(NFT.LT.15) THEN
        WRITE(6, *) 'In VUEHGEN: need more fitting functions. NFT =',NFT
        WRITE(7, *) 'In VUEHGEN: need more fitting functions. NFT =',NFT
        STOP
      ELSEIF(NFT+1.GT.MFT) THEN
        WRITE(6, *) 'In VUEHGEN: too many fitting functions. NFT =',NFT
        WRITE(7, *) 'In VUEHGEN: too many fitting functions. NFT =',NFT
        STOP
      ENDIF
C
C     NUMBER OF FITTING FUNCTIONS
      NUEH(IZ) = NFT
C
C     ZEROTH, SECOND AND FOURTH MOMENTS OF POLARISED CHARGE DENSITY
      FAC = CMPF*CMPF/(PI*CV)
      FNT = RNUC(IZ)/CMPF
      RUEH(IZ,1) = 0.0D0
      IF(COUPLING.EQ.'ELEC'.OR.COUPLING.EQ.'MUON') THEN
        RUEH(IZ,2) = 2.0D0*FAC/5.0D0
        RUEH(IZ,3) = 6.0D0*FAC*CMPF*CMPF/7.0D0
        RUEH(IZ,3) = RUEH(IZ,3)*(1.0D0 + 14.0D0*FNT*FNT/9.0D0)
      ELSEIF(COUPLING.EQ.'PION') THEN
        RUEH(IZ,2) = FAC/10.0D0
        RUEH(IZ,3) = 7.0D0*FAC*CMPF*CMPF
        RUEH(IZ,3) = RUEH(IZ,3)*(1.0D0 + 7.0D0*FNT*FNT/3.0D0)
      ELSE
        WRITE(6, *) 'In VUEHGEN: unrecognised coupling field.',COUPLING
        WRITE(7, *) 'In VUEHGEN: unrecognised coupling field.',COUPLING
        RETURN
      ENDIF
C
C     NUMERICAL CALCULATION OF THE UEHLING POTENTIAL
C     1: FOURIER, 2: DIRECT INTEGRATION, 3: F+R CHEBYSHEV
      IPOT = 1
C
C**********************************************************************C
C     GENERATE POLARISATION FORM FACTOR                                C
C**********************************************************************C
C
      IF(IPOT.NE.1) GOTO 100
C
C     READ SPECTRAL FUNCTION FROM FILE OR RECALCULATE
      IF(COUPLING.EQ.'ELEC') THEN
C
        INQUIRE(FILE='spectral/pi_ueh_electron.dat',EXIST=FILETHERE)
        IF(FILETHERE) THEN
          OPEN(UNIT=8,FILE='spectral/pi_ueh_electron.dat',
     &                                               STATUS='UNKNOWN')
          REWIND(UNIT=8)
          DO I=0,NQVP
            READ(8, *) M,Q,PUEH(I)
            IF(M.NE.I.OR.Q.NE.QFM(I)) THEN
              FILETHERE = .FALSE.
            GOTO 111
            ENDIF
          ENDDO
111       CONTINUE
          CLOSE(UNIT=8)
        ENDIF
        IF(.NOT.FILETHERE) THEN
          OPEN(UNIT=8,FILE='spectral/pi_ueh_electron.dat',
     &                                               STATUS='UNKNOWN')
          REWIND(UNIT=8)
          DO I=0,NQVP
            PUEH(I) = PIUEHF(QFM(I),CMPF)
            WRITE(8, *) I,QFM(I),PUEH(I)
          ENDDO
          CLOSE(UNIT=8)
        ENDIF
C
      ELSEIF(COUPLING.EQ.'MUON') THEN
C
        INQUIRE(FILE='spectral/pi_ueh_muon.dat',EXIST=FILETHERE)
        IF(FILETHERE) THEN
          OPEN(UNIT=8,FILE='spectral/pi_ueh_muon.dat',STATUS='UNKNOWN')
          REWIND(UNIT=8)
          DO I=0,NQVP
            READ(8, *) M,Q,PUEH(I)
            IF(M.NE.I.OR.Q.NE.QFM(I)) THEN
              FILETHERE = .FALSE.
            GOTO 121
            ENDIF
          ENDDO
121       CONTINUE
          CLOSE(UNIT=8)
        ENDIF
        IF(.NOT.FILETHERE) THEN
          OPEN(UNIT=8,FILE='spectral/pi_ueh_muon.dat',STATUS='UNKNOWN')
          REWIND(UNIT=8)
          DO I=0,NQVP
            PUEH(I) = PIUEHF(QFM(I),CMPF)
            WRITE(8, *) I,QFM(I),PUEH(I)
          ENDDO
          CLOSE(UNIT=8)
        ENDIF
C
      ELSEIF(COUPLING.EQ.'PION') THEN
C
        INQUIRE(FILE='spectral/pi_ueh_pion+.dat',EXIST=FILETHERE)
        IF(FILETHERE) THEN
          OPEN(UNIT=8,FILE='spectral/pi_ueh_pion+.dat',STATUS='UNKNOWN')
          REWIND(UNIT=8)
          DO I=0,NQVP
            READ(8, *) M,Q,PUEH(I)
            IF(M.NE.I.OR.Q.NE.QFM(I)) THEN
              FILETHERE = .FALSE.
            GOTO 131
            ENDIF
          ENDDO
131       CONTINUE
          CLOSE(UNIT=8)
        ENDIF
        IF(.NOT.FILETHERE) THEN
          OPEN(UNIT=8,FILE='spectral/pi_ueh_pion+.dat',STATUS='UNKNOWN')
          REWIND(UNIT=8)
          DO I=0,NQVP
            PUEH(I) = PIUEHB(QFM(I),CMPF)
            WRITE(8, *) I,QFM(I),PUEH(I)
          ENDDO
          CLOSE(UNIT=8)
        ENDIF
C
      ENDIF
CC
CC     NUCLEAR FORM FACTOR
C      DO I=0,NQVP
C        FMQ(I) = 0.0D0
C        Q = QFM(I)
C        FAC = 4.0D0*XNUC(IZ,0)
C        FMQ(I) = FMQ(I) + FNUC(IZ,0)*DEXP(-Q*Q/FAC)
C        DO K=1,NNUC(IZ)
C          FAC = 4.0D0*XNUC(IZ,K)
C          X32 = XNUC(IZ,K)*DSQRT(XNUC(IZ,K))
C          PRE = 0.25D0*PI12*Q*Q/X32
C          FMQ(I) = FMQ(I) - PRE*FNUC(IZ,K)*DEXP(-Q*Q/FAC)
C        ENDDO
C      ENDDO
CC
C     NUCLEAR CHARGE DISTRIBUTION ALONG RADIAL GRID
      DO N=0,NRVP
        R = RFM(N)
        RHOTMP(N) = RHONUC(NMDL(IZ),IZ,R)
      ENDDO

C     CORRESPONDING NUCLEAR FORM FACTOR
      DO I=0,NQVP
        Q = QFM(I)
        FMQ(I) = 0.0d0
        DO N=0,NRVP
          R = RFM(N)
          FMQ(I) = FMQ(I)
     &           + EXTINT11(R*R*RHOTMP(N)*SPHBJ0(Q*R)*REP(N),N,NRVP)
        ENDDO
        FMQ(I) = 5.0D0*FMQ(I)/2.99376D+5
        FMQ(I) = 4.0D0*PI*FMQ(I)
      ENDDO
C
C     CALCULATE THE POLARISATION FORM FACTOR AS A PRODUCT
      DO I=0,NQVP
        GUEH(I) = FMQ(I)*PUEH(I)
      ENDDO
C
C     NOW THE POTENTIAL AT ANY R IS A WEIGHTED FIT
C
C     EXPLICITLY CALCULATE V0 (SPECIAL CASE)
      V0 = 0.0D0
      DO I=0,NQVP
        Q = QFM(I)
        V0 = V0 + EXTINT11(QEP(I)*GUEH(I),I,NQVP)
      ENDDO
      V0 =-2.0D0*5.0D0*V0/(PI*2.99376D+5)
C
100   CONTINUE
C
C**********************************************************************C
C     GENERATE AND IMPORT THE INTERPOLATED CHI FUNCTIONS               C
C**********************************************************************C
C
      IF(IPOT.NE.2) GOTO 200
C
C     ONLY SUPPORTS FERMIONIC FIELDS
      IF(COUPLING.NE.'ELEC'.AND.COUPLING.NE.'MUON') THEN
        WRITE(6, *) 'In VUEHGEN: unrecognised coupling field.',COUPLING
        WRITE(7, *) 'In VUEHGEN: unrecognised coupling field.',COUPLING
        RETURN
      ENDIF
C
C     IMPORT SPLINE DATA FOR X*CHI_1(X) GIVEN 0 <= X <= XSPL
      INQUIRE(FILE='spline/chi1x.dat',EXIST=FILETHERE)
      IF(.NOT.FILETHERE) CALL CHIGEN
      OPEN(UNIT=50,FILE='spline/chi1x.dat',STATUS='UNKNOWN')
      REWIND(UNIT=50)
      READ(50, *) XSPL
      DO N=0,NLW
        READ(50, *) XS(N),Y1S(N),D1S(N)
      ENDDO
      CLOSE(UNIT=50)
C
C     IMPORT SPLINE DATA FOR CHI_2(X) GIVEN 0 <= X <= XSPL
      INQUIRE(FILE='spline/chi2_small.dat',EXIST=FILETHERE)
      IF(.NOT.FILETHERE) CALL CHIGEN
      OPEN(UNIT=52,FILE='spline/chi2_small.dat',STATUS='UNKNOWN')
      REWIND(UNIT=52)
      READ(52, *) XSPL
      DO N=0,NLW
        READ(52, *) XS(N),Y2S(N),D2S(N)
      ENDDO
      CLOSE(UNIT=52)
C
C     IMPORT SPLINE DATA FOR CHI_2(X) GIVEN X >= XSPL
      INQUIRE(FILE='spline/chi2_big.dat',EXIST=FILETHERE)
      IF(.NOT.FILETHERE) CALL CHIGEN
      OPEN(UNIT=53,FILE='spline/chi2_big.dat',STATUS='UNKNOWN')
      REWIND(UNIT=53)
      READ(53, *) XSPL
      DO N=0,NUP
        READ(53, *) XB(N),Y2B(N),D2B(N)
      ENDDO
      CLOSE(UNIT=53)
C
C**********************************************************************C
C     INTEGRATION PARAMETERS AND UEHLING POTENTIAL AT THE ORIGIN       C
C**********************************************************************C
C
C     SEARCH FOR CHARGE RADIUS SMAX FOR WHICH RHO(SMAX) < 1.0D-16
      SMAX = 0.0D0
60    SMAX = SMAX + 0.1D0/CFM
      PMAX = SMAX*RHONUC(NMDL(IZ),IZ,SMAX)
      IF(DABS(PMAX).GT.1.0D-16) GOTO 60
C
C     SOURCE CHARGE STEP SIZE
      HS = SMAX/DFLOAT(NSRC)
C
C     EVALUATE CHARGE DENSITY ON UNIFORMLY-SPACED GRID
      DO M=0,NSRC
        S = HS*DFLOAT(M)
        RHO(M) = RHONUC(NMDL(IZ),IZ,S)
      ENDDO
C
C     UEHLING POTENTIAL ORIGIN VALUE
      V0 = 0.0D0
      DO M=1,NSRC
C
C       SOURCE CHARGE RADIUS
        S  = HS*DFLOAT(M)
C
C       X*CHI_1(X) VALUE (ARGUMENT WILL ALWAYS BE LESS THAN 0.6D0)
        X0 = 2.0D0*S/CMPF
        CALL SPLNINT(XS,Y1S,D1S,X0,C0,NLW)
C
C       CONTRIBUTION TO INTEGRAND
        V0 = V0 + EXTINT11(RHO(M)*C0,M,NSRC)
C
      ENDDO
C
C     MULTIPLICATIVE FACTORS FOR VUEH(0)
      V0 = 5.0D0*HS*V0/2.99376D+5
      V0 =-4.0D0*CMPF*V0/(3.0D0*CV)
C
200   CONTINUE
C
C**********************************************************************C
C     UEHLING POTENTIAL VALUES ON THE TESTING GRID.                    C
C**********************************************************************C
C
C     UEHLING POTENTIAL SAMPLING VALUES (FOR R-SQUARED CALCULATION)
      RLIM = 250.0D0
C
C     STORE R*R*V(R) ON THE UNIFORM GRID IN VF(0:NPTS)
      RU(0) = 0.0D0
      VF(0) = 0.0D0
C
C     POTENTIAL BY CONVOLUTION
      IF(IPOT.EQ.1) THEN
C
        DO IPTS=1,NPTS
C
C         RADIUS ON UNIFORM SCALE
          R = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C         STORE FOR LATER
          RU(IPTS) = R
C
C         NUMERICAL INVERSE HANKEL TRANSFORM (FOR A POTENTIAL)
          CALL VRIHNKT(GUEH,R,VF(IPTS))
C
        ENDDO
C
C     POTENTIAL BY DIRECT INTEGRATION
      ELSEIF(IPOT.EQ.2) THEN
C
        DO IPTS=1,NPTS
C
C         RADIUS ON UNIFORM SCALE
          R = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C         STORE FOR LATER
          RU(IPTS) = R
C
C         INITIALISE POTENTIAL VALUE
          VF(IPTS) = 0.0D0
C
C         INTEGRATE OVER CHARGE SOURCE
          DO M=0,NSRC
C
C           SOURCE RADIUS
            S  = HS*DFLOAT(M)
C
C           CHI FUNCTION ARGUMENTS
            XM = 2.0D0*DABS(R-S)/CMPF
            XP = 2.0D0*DABS(R+S)/CMPF
C
C           COMPONENTS OF INTEGRAND
            IF(0.5D0*(XM+XP).LT.XSPL) THEN
              CALL SPLNINT(XS,Y2S,D2S,XM,CM,NLW)
              CALL SPLNINT(XS,Y2S,D2S,XP,CP,NLW)    
            ELSE
              CALL SPLNINT(XB,Y2B,D2B,XM,CM,NUP)
              CALL SPLNINT(XB,Y2B,D2B,XP,CP,NUP)
            ENDIF
C
C           PERFORM THE MAPPING CHI(X) = SPLINE(X)*EXP(-X)
            CM = CM*DEXP(-XM)
            CP = CP*DEXP(-XP)
C
C           CONTRIBUTION TO INTEGRAND
            VF(IPTS) = VF(IPTS) + EXTINT11(S*RHO(M)*(CM-CP),M,NSRC)
C
          ENDDO
C
C         INTEGRATION WEIGHTING FACTORS
          VF(IPTS) = 5.0D0*HS*VF(IPTS)/2.99376D+5
C
C         OTHER FACTORS NEEDED FOR V(R)
          VF(IPTS) =-2.0D0*CMPF*VF(IPTS)/(3.0D0*CV*R)
C
        ENDDO
C
C     FULLERTON/RINKER METHOD
      ELSEIF(IPOT.EQ.3) THEN
C
        ISWITCH = 1
C
        DO IPTS=1,NPTS
C
C         RADIUS ON UNIFORM SCALE
          R = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C         STORE FOR LATER
          RU(IPTS) = R
C
C         EVALUATE V AS IF IN ASYMPTOTIC REGION
          IF(R.NE.0.0D0) THEN
            V = FUNK(2.0D0*R/CMPF,1)/(PI*CMPF)
          ENDIF
C
C         EVALUATE V AS IF IN NUCLEAR REGION
          IF(ISWITCH.EQ.0) GOTO 202
C
C         INTIALISE POTENTIAL COUNTER
          VN = 0.0D0
C
C         INTEGRATE OVER CHARGE SOURCE
          DO M=0,NSRC
C
C           SET SOURCE RADIUS
            S = HS*DFLOAT(M)
C
C           FUNCTION ARGUMENTS
            XM = 2.0D0*DABS(R-S)/CMPF
            XP = 2.0D0*DABS(R+S)/CMPF
C
C           CHEBYSHEV APPROXIMATION FIT (FULLERTON+RINKER)
            CM = FUNK(XM,0)
            CP = FUNK(XP,0)
C
C           CONTRIBUTION TO INTEGRAND
            VN = VN + EXTINT11(S*RHO(M)*(CM-CP),M,NSRC)
C
          ENDDO
C
C         INTEGRATION WEIGHTING FACTORS
          VN = 5.0D0*HS*VN/2.99376D+5
C
C         IF FULL AND ASYMPTOTIC RESULTS ARE SIMILAR, FLIP SWITCH
          IF(DABS((VN-V)/VN)*ZNUC(IZ).LT.1.0D-5) THEN
            ISWITCH = 0
          ENDIF
C
C         OVERRIDE POTENTIAL NOW
          V = VN
C
202       CONTINUE
C
C         OTHER FACTORS AND FINAL VALUE V(R)
          VF(IPTS) =-2.0D0*CMPF*V/(3.0D0*CV*R)
C
        ENDDO
C      
      ENDIF
C
C**********************************************************************C
C     UEHLING POTENTIAL VALUES ON THE MATCHING GRID                    C
C**********************************************************************C
C
C     FIRST CONDITION MATCHES RMS RADIUS
      Y(1) = RUEH(IZ,2)
C
C     SECOND CONDITION MATCHES FOURTH RADIAL MOMENT
      Y(2) = RUEH(IZ,3)
C
C     DIVIDE NFT POINTS INTO DIRECT AND WEIGHTED REGIONS
C     NDCT = (NFT)/4
      NDCT = 6
      NWGT = NFT-NDCT-2
C
C     UEHLING POTENTIAL MATCHING VALUES
      R0 = 0.0D0
      RC = 2.0D0
C     RC = 3.0D0
      RN = 240.0D0
C     RN = 200.0D0
      HD = (RC-R0)/DFLOAT(NDCT)
      HW = DLOG(RN/RC)/DFLOAT(NWGT-1)
C
C     RADIAL MATCHING POINTS
      DO IDCT=1,NDCT
        RM(IDCT     ) = RNUC(IZ)*(R0 + DFLOAT(IDCT-1)*HD)
      ENDDO
      DO IWGT=1,NWGT
        RM(IWGT+NDCT) = RNUC(IZ)*RC*DEXP(DFLOAT(IWGT-1)*HW)
      ENDDO
C
C     STORE R*R*V(R) MATCHING VALUES IN Y(NFT) -- MATRIX EQN SOLUTIONS
      Y(3) = V0
C
      DO IFT=4,NFT
C
C       RADIAL MATCHING POINT
        R = RM(IFT-2)
C
C       INITIALISE VALUE FOR THE POTENTIAL HERE
        Y(IFT) = 0.0D0
C
C       POTENTIAL BY CONVOLUTION
        IF(IPOT.EQ.1) THEN
C
C         NUMERICAL INVERSE HANKEL TRANSFORM (FOR A POTENTIAL)
          CALL VRIHNKT(GUEH,R,Y(IFT))
C
C       POTENTIAL BY DIRECT INTEGRATION
        ELSEIF(IPOT.EQ.2) THEN
C
C         INTEGRATE OVER CHARGE SOURCE
          DO M=0,NSRC
C
C           SET SOURCE RADIUS
            S  = HS*DFLOAT(M)
C
C           CHI FUNCTION ARGUMENTS
            XM = 2.0D0*DABS(R-S)/CMPF
            XP = 2.0D0*DABS(R+S)/CMPF
C
C           COMPONENTS OF INTEGRAND
            IF(0.5D0*(XM+XP).LT.XSPL) THEN
              CALL SPLNINT(XS,Y2S,D2S,XM,CM,NLW)
              CALL SPLNINT(XS,Y2S,D2S,XP,CP,NLW)
            ELSE
              CALL SPLNINT(XB,Y2B,D2B,XM,CM,NUP)
              CALL SPLNINT(XB,Y2B,D2B,XP,CP,NUP)
            ENDIF
C
C           PERFORM THE MAPPING CHI(X) = SPLINE(X)*EXP(-X)
            CM = CM*DEXP(-XM)
            CP = CP*DEXP(-XP)
C
C           CONTRIBUTION TO INTEGRAND
            Y(IFT) = Y(IFT) + EXTINT11(S*RHO(M)*(CM-CP),M,NSRC)
C
          ENDDO
C
C         INTEGRATION WEIGHTING FACTORS
          Y(IFT) = 5.0D0*HS*Y(IFT)/2.99376D+5
C
C         OTHER FACTORS AND UEHLING POTENTIAL V(R)
          Y(IFT) =-2.0D0*CMPF*Y(IFT)/(3.0D0*CV*R)
C
C       FULLERTON/RINKER METHOD
        ELSEIF(IPOT.EQ.3) THEN
C
          ISWITCH = 1
C
C         EVALUATE V AS IF IN ASYMPTOTIC REGION
          IF(R.NE.0.0D0) THEN
            V = FUNK(2.0D0*R/CMPF,1)/(PI*CMPF)
          ENDIF
C
C         EVALUATE V AS IF IN NUCLEAR REGION
          IF(ISWITCH.EQ.0) GOTO 212
C
C         INTIALISE POTENTIAL COUNTER
          VN = 0.0D0
C
C         INTEGRATE OVER CHARGE SOURCE
          DO M=0,NSRC
C
C           SET SOURCE RADIUS
            S = HS*DFLOAT(M)
C
C           FUNCTION ARGUMENTS
            XM = 2.0D0*DABS(R-S)/CMPF
            XP = 2.0D0*DABS(R+S)/CMPF
C
C           CHEBYSHEV APPROXIMATION FIT (FULLERTON+RINKER)
            CM = FUNK(XM,0)
            CP = FUNK(XP,0)
C
C           CONTRIBUTION TO INTEGRAND
            VN = VN + EXTINT11(S*RHO(M)*(CM-CP),M,NSRC)
C
          ENDDO
C
C         INTEGRATION WEIGHTING FACTORS
          VN = 5.0D0*HS*VN/2.99376D+5
C
C         IF FULL AND ASYMPTOTIC RESULTS ARE SIMILAR, FLIP SWITCH
          IF(DABS((VN-V)/VN)*ZNUC(IZ).LT.1.0D-5) THEN
            ISWITCH = 0
          ENDIF
C
C         OVERRIDE POTENTIAL NOW
          V = VN
C
212       CONTINUE
C
C         OTHER FACTORS AND FINAL VALUE V(R)
          Y(IFT) =-2.0D0*CMPF*V/(3.0D0*CV*R)
C
        ENDIF
C
C       DATA POINTS IN THE WEIGHTED REGION MUST BE MULTIPLIED BY R*R
        IF(IFT-2.GT.NDCT) THEN
          Y(IFT) = R*R*Y(IFT)
        ENDIF
C
      ENDDO
C
C**********************************************************************C
C     GIVEN MODIFIED GEOMETRIC SET OF NFT PARAMETERS, FIND BEST ALPHA. C
C**********************************************************************C
C
C     GEOMETRIC PARAMETER SEARCH DETAILS
      ALPH = 2.00D-5
      BETA = 1.60D+0
      GAMA = 0.21D+0
C
C     NUMBER OF INCREMENTS IN SEARCH
      NAF = 1200
C
C     INITIALISE BEST-FIT R-SQUARED AND INDEX FOR ALPHA
      RSQBIG = 0.0D0
      IABG   = 0
C
C     ITERATE OVER ALPHA SEARCH
      DO IA=0,NAF
C
C       TRIAL ALPHA PARAMETER AND FROZEN BETA PARAMETER
        AF = ALPH*(1.0D0 + 3.0D0*DFLOAT(IA)/DFLOAT(NAF))/(RNUC(IZ)**2)
        BF = BETA
C
C       STORE GEOMETRIC SET OF PARAMETERS IN XUEH
        XI = AF
        DO IFT=1,NFT
          XUEH(IZ,IFT) = XI
          RT = DFLOAT(IFT)/DFLOAT(NFT+1)
          XI = BF*XI*(1.0D0 + GAMA*RT*RT)
        ENDDO
C
C       TRANSFER MATCHING VALUES V(R) AND R*R*V(R) INTO Z VECTOR
        DO IFT=1,NFT
          Z(IFT) = Y(IFT)
        ENDDO
C
C       SET UP MATRIX EQUATIONS FOR SET OF FITTING GAUSSIANS
C
C       MATCH THE SECOND MOMENT
        DO JFT=1,NFT
          RAT = PI/XUEH(IZ,JFT)
          X(1,JFT) =-6.0D0*RAT*DSQRT(RAT)/(4.0D0*PI)
        ENDDO
C
C       MATCH THE FOURTH MOMENT
        DO JFT=1,NFT
          RAT = PI/XUEH(IZ,JFT)
          X(2,JFT) =-30.0D0*RAT*DSQRT(RAT)/(4.0D0*PI*XUEH(IZ,JFT))
        ENDDO
C
C       MATCH THE POTENTIAL ITSELF
        DO IFT=3,NFT
C
C         MATCHING RADIUS
          R = RM(IFT-2)
C
C         GAUSSIAN POTENTIALS AT THIS RADIUS
          DO JFT=1,NFT
            IF(IFT-2.LE.NDCT) THEN
              X(IFT,JFT) =     DEXP(-XUEH(IZ,JFT)*R*R)
            ELSE
              X(IFT,JFT) = R*R*DEXP(-XUEH(IZ,JFT)*R*R)
            ENDIF
          ENDDO
C
        ENDDO
C
C       SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
        CALL DGESV(NFT,1,X,NFT,IPIV,Z,NFT,INFO)
C
C       TRANSFER AMPLITUDES VALUES INTO FUEH ARRAY
        DO IFT=1,NFT
          FUEH(IZ,IFT) = Z(IFT)
        ENDDO
C
C       BEST-FIT GAUSSIAN POTENTIAL VALUES ACROSS UNIFORM GRID
        DO IPTS=0,NPTS
C
C         RADIUS ON UNIFORM SCALE
          R = RU(IPTS)
C
C         INITIALISE BEST-FIT POTENTIAL AT THIS RADIUS
          VG(IPTS) = 0.0D0
          DO IFT=1,NFT
            VG(IPTS) = VG(IPTS) + FUEH(IZ,IFT)*DEXP(-XUEH(IZ,IFT)*R*R)
          ENDDO
C
        ENDDO
C
C       AVERAGE POTENTIAL VALUE
        YB = 0.0D0
        DO IPTS=0,NPTS
          YB = YB + VG(IPTS)
        ENDDO
        YB = YB/DFLOAT(NPTS+1)
C
C       PREPARATION FOR R-SQUARED VALUE
        SRES = 0.0D0
        STOT = 0.0D0
        DO IPTS=0,NPTS
          SRES = SRES + RU(IPTS)*RU(IPTS)*(VG(IPTS)-VF(IPTS))**2
          STOT = STOT + RU(IPTS)*RU(IPTS)*(VG(IPTS)-YB      )**2
        ENDDO
C
C       R-SQUARED VALUE AND BEST-FIT CHECK
        RSQ = 1.0D0 - SRES/STOT
        IF(RSQ.GT.RSQBIG) THEN
          IABG   = IA
          RSQBIG = RSQ
        ENDIF
C
      ENDDO
C
C     IF NO SOLUTION CAN BE FOUND, ALERT THE USER
      IF(RSQBIG.LT.0.0D0) THEN
        WRITE(6, *) 'In VUEHGEN: best-fit search failed. IZ =',IZ
        WRITE(7, *) 'In VUEHGEN: best-fit search failed. IZ =',IZ
        STOP
      ENDIF
C
C**********************************************************************C
C     BEST-FIT ALPHA HAS NOW BEEN FOUND -- RESTORE THESE RESULTS.      C
C**********************************************************************C
C
C     BEST-FIT VALUES FOR ALPHA, BETA AND GAMMA
      AF = ALPH*(1.0D0 + 3.0D0*DFLOAT(IABG)/DFLOAT(NAF))
      BF = BETA
      GF = GAMA
C
C     STORE GEOMETRIC SET OF PARAMETERS IN XUEH
      XI = AF/(RNUC(IZ)**2)
      DO IFT=1,NFT
        XUEH(IZ,IFT) = XI
        RT = DFLOAT(IFT)/DFLOAT(NFT+1)
        XI = BF*XI*(1.0D0 + GF*RT*RT)
      ENDDO
C
C     TRANSFER MATCHING VALUES V(R) AND R*R*V(R) INTO Z VECTOR
      DO IFT=1,NFT
        Z(IFT) = Y(IFT)
      ENDDO
C
C     SET UP MATRIX EQUATIONS FOR SET OF FITTING GAUSSIANS
C
C     MATCH THE SECOND MOMENT
      DO JFT=1,NFT
        RAT = PI/XUEH(IZ,JFT)
        X(1,JFT) =-6.0D0*RAT*DSQRT(RAT)/(4.0D0*PI)
      ENDDO
C
C     MATCH THE FOURTH MOMENT
      DO JFT=1,NFT
        RAT = PI/XUEH(IZ,JFT)
        X(2,JFT) =-30.0D0*RAT*DSQRT(RAT)/(4.0D0*PI*XUEH(IZ,JFT))
      ENDDO
C
C     MATCH THE POTENTIAL ITSELF
      DO IFT=3,NFT
C
C       MATCHING RADIUS
        R = RM(IFT-2)
C
C       GAUSSIAN POTENTIALS AT THIS RADIUS
        DO JFT=1,NFT
          IF(IFT-2.LE.NDCT) THEN
            X(IFT,JFT) =     DEXP(-XUEH(IZ,JFT)*R*R)
          ELSE
            X(IFT,JFT) = R*R*DEXP(-XUEH(IZ,JFT)*R*R)
          ENDIF
        ENDDO
C
      ENDDO
C
C     SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
      CALL DGESV(NFT,1,X,NFT,IPIV,Z,NFT,INFO)
C
C     TRANSFER AMPLITUDES VALUES INTO FUEH ARRAY
      DO IFT=1,NFT
        FUEH(IZ,IFT) = Z(IFT)
      ENDDO
C
C**********************************************************************C
C     RESULTS: RADIAL MOMENTS                                          C
C**********************************************************************C
C
C     BEST-FIT NUCLEAR UEHLING RADIAL MOMENTS (IN HARTREE)
      DO K=-2,5
        BIN = 0.0D0
        DO IFT=1,NFT
          FR = FUEH(IZ,IFT)
          XI = XUEH(IZ,IFT)
          XR = XI**(K+1)
          BIN = BIN-FR*DFLOAT(K)*GAMHLF(K+3)/DSQRT(XR)
        ENDDO
        RK(K) = BIN
      ENDDO
C
C**********************************************************************C
C     RESULTS: POLARISATION CHARGE DENSITY                             C
C**********************************************************************C
C
C     VACUUM POLARISATION CHARGE DENSITY
      DO IPTS=0,NPTS
C
C       RADIUS R
        R = RU(IPTS)
C
C       INTIALISE CHARGE COUNTER
        Q = 0.0D0
        DO IFT=1,NFT
          FR = FUEH(IZ,IFT)
          XI = XUEH(IZ,IFT)
          Q  = Q + 2.0D0*FR*XI*R*R*(3.0D0-2.0D0*R*R*XI)*DEXP(-XI*R*R)
        ENDDO
        POL(IPTS) = Q
C
      ENDDO
C
C     LOCATION OF CHARGE ZERO (LINE OF BEST FIT BETWEEN TWO POINTS)
      RZER = 0.0D0
      DO IPTS=NPTS,0,-1
        IF(POL(IPTS).LT.0.0D0.AND.RU(IPTS).LT.10.0D0/CFM) THEN
          PD   = POL(IPTS)-POL(IPTS-1)
          RZER = (RU(IPTS-1)*POL(IPTS)-RU(IPTS)*POL(IPTS-1))/PD
          GOTO 51
        ENDIF
      ENDDO
51    CONTINUE
C
C     INTEGRATED POLARISED CHARGE UP TO THIS ZERO
      QPOL = 0.0D0
      DO IFT=1,NFT
        FR = FUEH(IZ,IFT)
        XI = XUEH(IZ,IFT)
        QPOL = QPOL - 2.0D0*(RZER**3)*XI*FR*DEXP(-XI*RZER*RZER)
      ENDDO
C
C**********************************************************************C
C     PRINT RESULTS IF PROMPTED                                        C
C**********************************************************************C
C
C     SKIP OUTPUT SECTION UNLESS PROMPTED
      IF(IWRT.EQ.0) GOTO 20
C
C     SOLUTION SET
35    FORMAT(1X,A,I2,A,I3,A,A,A)
36    FORMAT(1X,A,F17.15,1X,A)
37    FORMAT(1X,A,F14.12,A,F8.6,A,F8.6)
38    FORMAT(1X,A,I2,A,ES15.8,A,ES15.8,A)
39    FORMAT(1X,A,ES17.10,1X,A)
      WRITE(6, *) ''
      WRITE(7, *) ''
      WRITE(6, *) 'Coupling field: ',COUPLING
      WRITE(7, *) 'Coupling field: ',COUPLING
      WRITE(6, *) 'Mass of virtual particle: ',RV
      WRITE(7, *) 'Mass of virtual particle: ',RV
      WRITE(6,35) 'Uehling best-fit for centre IZ = ',IZ,
     &                  '  (',INT(ANUC(IZ)),'^',ELMT(INT(ZNUC(IZ))),')'
      WRITE(7,35) 'Uehling best-fit for centre IZ = ',IZ,
     &                  '  (',INT(ANUC(IZ)),'^',ELMT(INT(ZNUC(IZ))),')'
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
      WRITE(6,36) 'Here ξ0 = 1/R_n^2   and   R_n  = ',RNUC(IZ),'a0'
      WRITE(7,36) 'Here ξ0 = 1/R_n^2   and   R_n  = ',RNUC(IZ),'a0'
      WRITE(6,36) '                                 ',RNUC(IZ)*CFM,'fm'
      WRITE(7,36) '                                 ',RNUC(IZ)*CFM,'fm'
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
      WRITE(6,37) 'α = ',AF,' ξ0,   β = ',BF,',   γ = ',GF
      WRITE(7,37) 'α = ',AF,' ξ0,   β = ',BF,',   γ = ',GF
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
      WRITE(6,36) 'Least-squares best fit:       R^2 = ',RSQBIG
      WRITE(7,36) 'Least-squares best fit:       R^2 = ',RSQBIG
      WRITE(6,39) '1-R^2                             = ',1.0D0-RSQBIG
      WRITE(7,39) '1-R^2                             = ',1.0D0-RSQBIG
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
      DO IFT=1,NFT
        FR = FUEH(IZ,IFT)
        XI = RNUC(IZ)*RNUC(IZ)*XUEH(IZ,IFT)
        WRITE(6,38) 'f_',IFT,'(r): ',FR,'*exp(-',XI,' ξ0 r^2)'
        WRITE(7,38) 'f_',IFT,'(r): ',FR,'*exp(-',XI,' ξ0 r^2)'
      ENDDO
C
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
C
C     LOCATION OF CHARGE ZERO
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('-',50)
      WRITE(7, *) REPEAT('-',50)
      WRITE(6, *) 'Charge zero        :',RZER*CFM,' fm'
      WRITE(7, *) 'Charge zero        :',RZER*CFM,' fm'
      WRITE(6, *) 'Charge polarisation:',QPOL
      WRITE(7, *) 'Charge polarisation:',QPOL
      WRITE(6, *) REPEAT('-',50)
      WRITE(7, *) REPEAT('-',50)
C
C**********************************************************************C
C     BEST-FIT UEHLING POTENTIAL ON THE TESTING GRID.                  C
C**********************************************************************C
C
C     THIS TAKES TIME -- SKIP UNLESS YOU WANT QUADRATURE/PLOTTING
C     GOTO 20
C
C     BEST-FIT GAUSSIAN POTENTIAL VALUES ACROSS UNIFORM GRID
      DO IPTS=0,NPTS
C
C       RADIUS ON UNIFORM SCALE
        R = RU(IPTS)
C
C       INITIALISE BEST-FIT POTENTIAL AT THIS RADIUS
        VG(IPTS) = 0.0D0
        DO IFT=1,NFT
          VG(IPTS) = VG(IPTS) + FUEH(IZ,IFT)*DEXP(-XUEH(IZ,IFT)*R*R)
        ENDDO
C
      ENDDO
C
C**********************************************************************C
C     EXACT UEHLING POTENTIAL (WEIGHTED BY R*R) ON COMPOSITE GRID.     C
C**********************************************************************C
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
      DO N=0,NLIN
        RAD(N) = RORI + HL*DFLOAT(N)
      ENDDO
C
      DO N=0,NEXP
        RAD(N+NLIN) = RMID*DEXP(HE*DFLOAT(N))
      ENDDO
C
C     ORIGIN VALUE IS A SPECIAL CASE
      VVAC(IZ,0) = 0.0D0
C
C     POTENTIAL BY CONVOLUTION
      IF(IPOT.EQ.1) THEN
C
C       APPROXIMATE POTENTIAL EVERYWHERE
        AMP =-1.0D0/(PI*CV)
        RLONG = 200.0D0/(CFM*RV)
        IF(COUPLING.EQ.'ELEC'.OR.COUPLING.EQ.'MUON') THEN
          CALL VPOLPOT(CHIUEHF,CMPF,RAD,VTMP,AMP,RNUC(IZ),RLONG)
        ELSEIF(COUPLING.EQ.'PION') THEN
          CALL VPOLPOT(CHIUEHB,CMPF,RAD,VTMP,AMP,RNUC(IZ),RLONG)
        ENDIF
        DO N=1,NRAD
          R = RAD(N)
          VVAC(IZ,N) = R*R*VTMP(N)
        ENDDO
C
C       OVERRIDE WHEN R < 200 fm
        DO N=1,NRAD
C
C         RADIUS R AND INITIALISE COUNTER FOR POTENTIAL
          R = RAD(N)
C
C         NUMERICAL INVERSE HANKEL TRANSFORM (FOR A POTENTIAL)
          IF(R.LE.200.0D0/CFM) THEN
            CALL VRIHNKT(GUEH,R,V)
            VVAC(IZ,N) = R*R*V
          ENDIF
C
        ENDDO
C
C     POTENTIAL BY DIRECT INTEGRATION
      ELSEIF(IPOT.EQ.2) THEN
C
C       VALUE OF R*R*V(R) AT EACH OF THESE RADII
        DO N=1,NRAD
C
C         RADIUS R AND INITIALISE COUNTER FOR POTENTIAL
          R = RAD(N)
C
C         INTIALISE POTENTIAL COUNTER
          V = 0.0D0
C
C         INTEGRATE OVER CHARGE SOURCE
          DO M=0,NSRC
C
C           SET SOURCE RADIUS
            S  = HS*DFLOAT(M)
C
C           CHI FUNCTION ARGUMENTS
            XM = 2.0D0*DABS(R-S)/CMPF
            XP = 2.0D0*DABS(R+S)/CMPF
C
C           COMPONENTS OF INTEGRAND
            IF(0.5D0*(XM+XP).LT.XSPL) THEN
              CALL SPLNINT(XS,Y2S,D2S,XM,CM,NLW)
              CALL SPLNINT(XS,Y2S,D2S,XP,CP,NLW)    
            ELSE
              CALL SPLNINT(XB,Y2B,D2B,XM,CM,NUP)
              CALL SPLNINT(XB,Y2B,D2B,XP,CP,NUP)
            ENDIF
C
C           PERFORM THE MAPPING CHI(X) = SPLINE(X)*EXP(-X)
            CM = CM*DEXP(-XM)
            CP = CP*DEXP(-XP)
C
C           CONTRIBUTION TO INTEGRAND
            V = V + EXTINT11(S*RHO(M)*(CM-CP),M,NSRC)
C
          ENDDO
C
C         INTEGRATION WEIGHTING FACTORS
          V = 5.0D0*HS*V/2.99376D+5
C
C         OTHER FACTORS AND FINAL VALUE R*R*V(R)
          VVAC(IZ,N) =-2.0D0*CMPF*R*V/(3.0D0*CV)
C
        ENDDO
C
C     FULLERTON/RINKER METHOD
      ELSEIF(IPOT.EQ.3) THEN
C
        ISWITCH = 1
C
C       VALUE OF R*R*V(R) AT EACH OF THESE RADII
        DO N=1,NRAD
C
C         RADIUS R AND INITIALISE COUNTER FOR POTENTIAL
          R = RAD(N)
C
C         EVALUATE V AS IF IN ASYMPTOTIC REGION
          IF(R.NE.0.0D0) THEN
            V = FUNK(2.0D0*R/CMPF,1)/(PI*CMPF)
          ENDIF
C
C         EVALUATE V AS IF IN NUCLEAR REGION
          IF(ISWITCH.EQ.0) GOTO 222
C
C         INTIALISE POTENTIAL COUNTER
          VN = 0.0D0
C
C         INTEGRATE OVER CHARGE SOURCE
          DO M=0,NSRC
C
C           SET SOURCE RADIUS
            S = HS*DFLOAT(M)
C
C           FUNCTION ARGUMENTS
            XM = 2.0D0*DABS(R-S)/CMPF
            XP = 2.0D0*DABS(R+S)/CMPF
C
C           CHEBYSHEV APPROXIMATION FIT (FULLERTON+RINKER)
            CM = FUNK(XM,0)
            CP = FUNK(XP,0)
C
C           CONTRIBUTION TO INTEGRAND
            VN = VN + EXTINT11(S*RHO(M)*(CM-CP),M,NSRC)
C           VN = VN + EXTINT5(S*RHO(M)*(CM-CP),M,NSRC)
C
          ENDDO
C
C         INTEGRATION WEIGHTING FACTORS
          VN = 5.0D0*HS*VN/2.99376D+5
C         VN = 2.0D0*HS*VN/45.0D0
C
C         IF FULL AND ASYMPTOTIC RESULTS ARE SIMILAR, FLIP SWITCH
          IF(DABS((VN-V)/VN)*ZNUC(IZ).LT.1.0D-5) THEN
            ISWITCH = 0
          ENDIF
C
C         OVERRIDE POTENTIAL NOW
          V = VN
C
222       CONTINUE
C
C         OTHER FACTORS AND FINAL VALUE R*R*V(R)
          VVAC(IZ,N) =-2.0D0*CMPF*R*V/(3.0D0*CV)
C
        ENDDO
C      
      ENDIF
C
C**********************************************************************C
C     PLOTTING SECTION                                                 C
C**********************************************************************C
C
      IF(IWRT.NE.1) GOTO 20
C
C     NAME TAG FOR PDF DOCUMENTS
      IF(INT(ANUC(IZ)).LT.10) THEN
        WRITE(ATRM,'(A,I1)') '00',INT(ANUC(IZ))
      ELSEIF(INT(ANUC(IZ)).LT.100) THEN
        WRITE(ATRM,'(A,I2)') '0',INT(ANUC(IZ))
      ELSE
        WRITE(ATRM,'(I3)') INT(ANUC(IZ))
      ENDIF
C
C     DETAILS COMMON TO ALL PLOTS
      XAXIS  = 'r/RNUC'
      KEY(1) = 'Exact Uehling'
      KEY(2) = 'Best fit'
C
C     UNIFORM GRID VALUES R*R*V(R)
      XOUT   = 'Uehling-uniform'//TRIM(ATRM)
      TITLE  = 'Uehling r^{2}*V(r) on uniform grid'
      YAXIS  = 'r^{2}*V(r)/Z'
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
      DO IPTS=0,NPTS
        R = RU(IPTS)
        WRITE(8, *) R/RNUC(IZ),R*R*VF(IPTS),R*R*VG(IPTS)
      ENDDO
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C     NUMERICAL EVALUATION FOR V(R)
      XOUT   = 'Uehling-pointwise'//TRIM(ATRM)
      TITLE  = 'Uehling V(r) on piecewise grid'
      YAXIS  = 'V(r)/Z'
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
      DO N=0,NRAD
        R = RAD(N)
        IF(R/RNUC(IZ).GT.5.0D0) GOTO 123
C       TRUE UEHLING POTENTIAL
        IF(N.EQ.0) THEN
          VTR = V0
        ELSE
          VTR = VVAC(IZ,N)/(R*R)
        ENDIF
C       APPROXIMATE VALUES
        V = 0.0D0
        DO IFT=1,NFT
          V = V + FUEH(IZ,IFT)*DEXP(-XUEH(IZ,IFT)*R*R)
        ENDDO
        WRITE(8, *) R*cfm,VTR,V
      ENDDO
123   CONTINUE
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C     UNIFORM GRID VALUES R*R*V(R)
      XOUT   = 'Uehling-charge'//TRIM(ATRM)
      TITLE  = 'Uehling rho(r) on uniform grid'
      YAXIS  = 'rho(r)/Z'
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
      DO IPTS=0,NPTS/25
        R = RU(IPTS)
C       NUMERICAL INVERSE HANKEL TRANSFORM (FOR A CHARGE DENSITY)
        CALL PRIHNKT(GUEH,R,RHOUEH)
        WRITE(8, *) R/rnuc(iz),-4.0D0*R*R*PI*RHOUEH,POL(IPTS)
      ENDDO
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
20    CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE VWKRGEN(IZ,NFT,IWRT,RV,RSQBIG,RZER,QPOL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  VV    VV WW        WW KK    KK RRRRRRR   GGGGGG  EEEEEEEE NN    NN  C
C  VV    VV WW        WW KK   KK  RR    RR GG    GG EE       NNN   NN  C
C  VV    VV WW   WW   WW KK  KK   RR    RR GG       EE       NNNN  NN  C
C  VV    VV WW  WWWW  WW KKKKK    RR    RR GG       EEEEEE   NN NN NN  C
C   VV  VV  WW WW  WW WW KK  KK   RRRRRRR  GG   GGG EE       NN  NNNN  C
C    VVVV   WWWW    WWWW KK   KK  RR    RR GG    GG EE       NN   NNN  C
C     VV    WW        WW KK    KK RR    RR  GGGGGG  EEEEEEEE NN    NN  C
C                                                                      C
C -------------------------------------------------------------------- C
C  VWKRGEN GENERATES A BEST-FIT GAUSSIAN SET FOR THE WICHMANN-KROLL    C
C  POTENTIAL ARISING FROM NUCLEUS IZ, USING NFT TOTAL GAUSSIANS.       C
C -------------------------------------------------------------------- C
C  MATCHING CRITERIA AVAILABLE FOR GAUSSIAN AMPLITUDES:                C
C   ▶ RADIAL ARGUMENTS SCALED BY NUCLEAR RMS RADIUS RNUC(IZ).          C
C   ▶ WICHMANN-KROLL POTENTIAL FIRST SCALED TO R*R*V(R).               C
C   ▶ ZEROTH MOMENT SATISFIED AUTOMATICALLY BUT SECOND MOMENT ENFORCED.C
C   ▶ POTENTIAL MATCHED TO NFT POINTS, FROM 3.0D0*RNUC TO 200.0D0*RNUC.C
C   ▶ MATCHING POINTS ARE EXPONENTIALLY SPACED.                        C
C   ▶ GAUSSIAN EXPONENTS FROM MODIFIED GEOMETRIC SERIES WITH THE SAME  C
C     BETA AND GAMMA, BUT ALPHA IS SCALED BY RNUC*RNUC.                C
C   ▶ WHEN THIS IS COMPLETE, POTENTIAL V(R) NEAR R=0.0D0 IS AUGMENTED  C
C     WITH SOME MORE GAUSSIANS WITH LARGER EXPONENTS IF V(R) IS NEEDED.C
C -------------------------------------------------------------------- C
C INPUT:                                                               C
C   ▶ IZ:   NUCLEAR CENTRE.                                            C
C   ▶ NFT:  DESIRED NUMBER OF FITTING FUNCTIONS IN THE SET.            C
C   ▶ IWRT: EXTENDED WRITTEN OUTPUT OPTION.                            C
C   ▶ RV:   MASS MODIFIER (1.0D0 FOR ELECTRONS, MASS RATIO OTHERWISE). C
C**********************************************************************C
      INCLUDE 'parameters.h'
      PARAMETER(NPTS=6000)
C
      LOGICAL FILETHERE
C
      EXTERNAL CHIWKRF
C
      CHARACTER*2  ELMT(120)
      CHARACTER*3  ATRM
      CHARACTER*5  NMDL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(2)
C
      DIMENSION VF(0:NPTS),VG(0:NPTS),RU(0:NPTS)
      DIMENSION X(NFT,NFT),Z(NFT),Y(NFT),RM(NFT),IPIV(NFT)
      DIMENSION RHO(0:NSRC),RHOTMP(0:NRVP)
      DIMENSION POL(0:NPTS),RK(-2:5)
      DIMENSION FMQ(0:NQVP),PWKR(0:NQVP),GWKR(0:NQVP)
      DIMENSION VTMP(0:NRAD)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BWKR/RWKR(MCT,3),FWKR(MCT,MFT),XWKR(MCT,MFT),NWKR(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/MDLV/ELMT
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/QGRD/QFM(0:NQVP),QEP(0:NQVP),QMIN,QMID,QMAX
      COMMON/RGRD/RFM(0:NRVP),REP(0:NRVP),RMIN,RMD,RMX
      COMMON/QWKR/RAD(0:NRAD),VVAC(MCT,0:NRAD),RORI,RMID,RMAX,NLIN,NEXP
C
C     MULTIPLIER FOR COMPTON WAVELENGTH (USE FOR MUON OR TAUON FIELD)
      CMPF = CMPW/RV
C
C     CHECK THAT THERE ARE ENOUGH FITTING FUNCTIONS
      IF(NFT.LT.15) THEN
        WRITE(6, *) 'In VWKRGEN: need more fitting functions. NFT =',NFT
        WRITE(7, *) 'In VWKRGEN: need more fitting functions. NFT =',NFT
        STOP
      ELSEIF(NFT+1.GT.MFT) THEN
        WRITE(6, *) 'In VWKRGEN: too many fitting functions. NFT =',NFT
        WRITE(7, *) 'In VWKRGEN: too many fitting functions. NFT =',NFT
        STOP
      ENDIF
C
C     NUMBER OF FITTING FUNCTIONS
      NWKR(IZ) = NFT
C
C     ZEROTH, SECOND AND FOURTH MOMENTS OF POLARISED CHARGE DENSITY
      FAC = CMPF*CMPF/(PI*CV*CV*CV)
      FNT = RNUC(IZ)/CMPF
      RWKR(IZ,1) = 0.0D0
      RWKR(IZ,2) =-9.0D0*FAC/106.0D0
      RWKR(IZ,3) =-46.0D0*FAC*CMPF*CMPF/55.0D0
      RWKR(IZ,3) = RWKR(IZ,3)*(1.0D0 + 825.0D0*FNT*FNT/2438.0D0)
C
C     NUMERICAL CALCULATION OF THE WICHMANN-KROLL POTENTIAL
C     1: FOURIER, 2: ------------------, 3: -------------
      IPOT = 1
C
C**********************************************************************C
C     GENERATE POLARISATION FORM FACTOR                                C
C**********************************************************************C
C
      IF(IPOT.NE.1) GOTO 100
C
C     READ SPECTRAL FUNCTION FROM FILE OR RECALCULATE
      INQUIRE(FILE='spectral/pi_wkr_electron.dat',EXIST=FILETHERE)
      IF(FILETHERE) THEN
        OPEN(UNIT=8,FILE='spectral/pi_wkr_electron.dat',
     &                                               STATUS='UNKNOWN')
        REWIND(UNIT=8)
        DO I=0,NQVP
          READ(8, *) M,Q,PWKR(I)
          IF(M.NE.I) THEN
c         IF(M.NE.I.OR.Q.NE.QFM(I)) THEN
          	FILETHERE = .FALSE.
            GOTO 111
          ENDIF
        ENDDO
111     CONTINUE
        CLOSE(UNIT=8)
      ENDIF
      IF(.NOT.FILETHERE) THEN
        WRITE(6, *) 'In VWKRGEN: cannot locate spectral function.'
        WRITE(7, *) 'In VWKRGEN: cannot locate spectral function.'
        RETURN
      ENDIF
110   CONTINUE
C
C     NUCLEAR FORM FACTOR
      DO I=0,NQVP
        FMQ(I) = 0.0D0
        Q = QFM(I)
        FAC = 4.0D0*XNUC(IZ,0)
        FMQ(I) = FMQ(I) + FNUC(IZ,0)*DEXP(-Q*Q/FAC)
        DO K=1,NNUC(IZ)
          FAC = 4.0D0*XNUC(IZ,K)
          X32 = XNUC(IZ,K)*DSQRT(XNUC(IZ,K))
          PRE = 0.25D0*PI12*Q*Q/X32
          FMQ(I) = FMQ(I) - PRE*FNUC(IZ,K)*DEXP(-Q*Q/FAC)
        ENDDO
      ENDDO
CC
C     NUCLEAR CHARGE DISTRIBUTION ALONG RADIAL GRID
      DO N=0,NRVP
        R = RFM(N)
        RHOTMP(N) = RHONUC(NMDL(IZ),IZ,R)
      ENDDO
C
C     CORRESPONDING NUCLEAR FORM FACTOR
      DO I=0,NQVP
        Q = QFM(I)
        DO N=0,NRVP
          R = RFM(N)
          FMQ(I) = FMQ(I)
     &           + EXTINT11(R*R*RHOTMP(N)*SPHBJ0(Q*R)*REP(N),N,NRVP)
        ENDDO
        FMQ(I) = 5.0D0*FMQ(I)/2.99376D+5
        FMQ(I) = 4.0D0*PI*FMQ(I)
      ENDDO
C
C     CALCULATE THE POLARISATION FORM FACTOR AS A PRODUCT
      DO I=0,NQVP
        GWKR(I) = FMQ(I)*PWKR(I)
      ENDDO
C
C     EXPLICITLY CALCULATE V0 (SPECIAL CASE)
      V0 = 0.0D0
      DO I=0,NQVP
        Q = QFM(I)
        V0 = V0 + EXTINT11(QEP(I)*GWKR(I),I,NQVP)
      ENDDO
      V0 =-2.0D0*5.0D0*V0/(PI*2.99376D+5)
C
100   CONTINUE
C
C**********************************************************************C
C     NUCLEAR CHARGE INTEGRATION PARAMETERS.                           C
C**********************************************************************C
C
C     SEARCH FOR CHARGE RADIUS SMAX FOR WHICH RHO(SMAX) < 1.0D-16
      SMAX = 0.0D0
60    SMAX = SMAX + 0.1D0/CFM
      PMAX = SMAX*RHONUC(NMDL(IZ),IZ,SMAX)
      IF(DABS(PMAX).GT.1.0D-16) GOTO 60
C
C     SOURCE CHARGE STEP SIZE
      HS = SMAX/DFLOAT(NSRC)
C
C     EVALUATE CHARGE DENSITY ON UNIFORMLY-SPACED GRID
      DO M=0,NSRC
        S = HS*DFLOAT(M)
        RHO(M) = RHONUC(NMDL(IZ),IZ,S)
      ENDDO
C
C**********************************************************************C
C     WICHMANN-KROLL POTENTIAL VALUES ON THE TESTING GRID.             C
C**********************************************************************C
C
C     WICHMANN-KROLL POTENTIAL SAMPLING VALUES (FOR R-SQUARED EVAL)
      RLIM = 250.0D0
C
C     STORE R*R*V(R) ON THE UNIFORM GRID IN VF(0:NPTS)
      RU(0) = 0.0D0
      VF(0) = 0.0D0
C
C     POTENTIAL BY CONVOLUTION
      IF(IPOT.EQ.1) THEN
C
        DO IPTS=1,NPTS
C
C         RADIUS ON UNIFORM SCALE
          R = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C         STORE FOR LATER
          RU(IPTS) = R
C
C         NUMERICAL INVERSE HANKEL TRANSFORM (FOR A POTENTIAL)
          CALL VRIHNKT(GWKR,R,VF(IPTS))
C
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C     WICHMANN-KROLL POTENTIAL VALUES ON THE MATCHING GRID.            C
C**********************************************************************C
C
C     FIRST CONDITION MATCHES RMS RADIUS
      Y(1) = RWKR(IZ,2)
C
C     SECOND CONDITION MATCHES FOURTH RADIAL MOMENT
      Y(2) = RWKR(IZ,3)
C
C     DIVIDE NFT POINTS INTO DIRECT AND WEIGHTED REGIONS
C     NDCT = (NFT)/4
      NDCT = 6
      NWGT = NFT-NDCT-2
C
C     WICHMANN-KROLL POTENTIAL MATCHING VALUES
      R0 = 0.0D0
      RC = 2.0D0
C     RC = 3.0D0
      RN = 200.0D0
C     RN = 200.0D0
      HD = (RC-R0)/DFLOAT(NDCT)
      HW = DLOG(RN/RC)/DFLOAT(NWGT-1)
C
C     RADIAL MATCHING POINTS
      DO IDCT=1,NDCT
        RM(IDCT     ) = RNUC(IZ)*(R0 + DFLOAT(IDCT-1)*HD)
      ENDDO
      DO IWGT=1,NWGT
        RM(IWGT+NDCT) = RNUC(IZ)*RC*DEXP(DFLOAT(IWGT-1)*HW)
      ENDDO
C
C     STORE R*R*V(R) MATCHING VALUES IN Y(NFT) -- MATRIX EQN SOLUTIONS
      Y(3) = V0
C
      DO IFT=4,NFT
C
C       RADIAL MATCHING POINT
        R = RM(IFT-2)
C
C       INITIALISE VALUE FOR THE POTENTIAL HERE
        Y(IFT) = 0.0D0
C
C       POTENTIAL BY CONVOLUTION
        IF(IPOT.EQ.1) THEN
C
C         NUMERICAL INVERSE HANKEL TRANSFORM (FOR A POTENTIAL)
          CALL VRIHNKT(GWKR,R,Y(IFT))
C
        ENDIF
C
C       DATA POINTS IN THE WEIGHTED REGION MUST BE MULTIPLIED BY R*R
        IF(IFT-2.GT.NDCT) THEN
          Y(IFT) = R*R*Y(IFT)
        ENDIF
C
      ENDDO
C
C**********************************************************************C
C     GIVEN MODIFIED GEOMETRIC SET OF NFT PARAMETERS, FIND BEST ALPHA. C
C**********************************************************************C
C
C     GEOMETRIC PARAMETER SEARCH DETAILS
      ALPH = 2.00D-5
      BETA = 1.60D+0
      GAMA = 0.21D+0
C
C     NUMBER OF INCREMENTS IN SEARCH
      NAF = 1200
C
C     INITIALISE BEST-FIT R-SQUARED AND INDEX FOR ALPHA
      RSQBIG = 0.0D0
      IABG   = 0
C
C     ITERATE OVER ALPHA SEARCH
      DO IA=0,NAF
C
C       TRIAL ALPHA PARAMETER AND FROZEN BETA PARAMETER
        AF = ALPH*(1.0D0 + 3.0D0*DFLOAT(IA)/DFLOAT(NAF))/(RNUC(IZ)**2)
        BF = BETA
C
C       STORE GEOMETRIC SET OF PARAMETERS IN XWKR
        XI = AF
        DO IFT=1,NFT
          XWKR(IZ,IFT) = XI
          RT = DFLOAT(IFT)/DFLOAT(NFT+1)
          XI = BF*XI*(1.0D0 + GAMA*RT*RT)
        ENDDO
C
C       TRANSFER MATCHING VALUES V(R) AND R*R*V(R) INTO Z VECTOR
        DO IFT=1,NFT
          Z(IFT) = Y(IFT)
        ENDDO
C
C       SET UP MATRIX EQUATIONS FOR SET OF FITTING GAUSSIANS
C
C       MATCH THE SECOND MOMENT
        DO JFT=1,NFT
          RAT = PI/XWKR(IZ,JFT)
          X(1,JFT) =-6.0D0*RAT*DSQRT(RAT)/(4.0D0*PI)
        ENDDO
C
C       MATCH THE FOURTH MOMENT
        DO JFT=1,NFT
          RAT = PI/XWKR(IZ,JFT)
          X(2,JFT) =-30.0D0*RAT*DSQRT(RAT)/(4.0D0*PI*XWKR(IZ,JFT))
        ENDDO
C
C       MATCH THE POTENTIAL ITSELF
        DO IFT=3,NFT
C
C         MATCHING RADIUS
          R = RM(IFT-2)
C
C         GAUSSIAN POTENTIALS AT THIS RADIUS
          DO JFT=1,NFT
            IF(IFT-2.LE.NDCT) THEN
              X(IFT,JFT) =     DEXP(-XWKR(IZ,JFT)*R*R)
            ELSE
              X(IFT,JFT) = R*R*DEXP(-XWKR(IZ,JFT)*R*R)
            ENDIF
          ENDDO
C
        ENDDO
C
C       SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
        CALL DGESV(NFT,1,X,NFT,IPIV,Z,NFT,INFO)
C
C       TRANSFER AMPLITUDES VALUES INTO FWKR ARRAY
        DO IFT=1,NFT
          FWKR(IZ,IFT) = Z(IFT)
        ENDDO
C
C       BEST-FIT GAUSSIAN POTENTIAL VALUES ACROSS UNIFORM GRID
        DO IPTS=0,NPTS
C
C         RADIUS ON UNIFORM SCALE
          R = RU(IPTS)
C
C         INITIALISE BEST-FIT POTENTIAL AT THIS RADIUS
          VG(IPTS) = 0.0D0
          DO IFT=1,NFT
            VG(IPTS) = VG(IPTS) + FWKR(IZ,IFT)*DEXP(-XWKR(IZ,IFT)*R*R)
          ENDDO
C
        ENDDO
C
C       AVERAGE POTENTIAL VALUE
        YB = 0.0D0
        DO IPTS=0,NPTS
          YB = YB + VG(IPTS)
        ENDDO
        YB = YB/DFLOAT(NPTS+1)
C
C       PREPARATION FOR R-SQUARED VALUE
        SRES = 0.0D0
        STOT = 0.0D0
        DO IPTS=0,NPTS
          SRES = SRES + RU(IPTS)*RU(IPTS)*(VG(IPTS)-VF(IPTS))**2
          STOT = STOT + RU(IPTS)*RU(IPTS)*(VG(IPTS)-YB      )**2
        ENDDO
C
C       R-SQUARED VALUE AND BEST-FIT CHECK
        RSQ = 1.0D0 - SRES/STOT
        IF(RSQ.GT.RSQBIG) THEN
          IABG   = IA
          RSQBIG = RSQ
        ENDIF
C
      ENDDO
C
C     IF NO SOLUTION CAN BE FOUND, ALERT THE USER
      IF(RSQBIG.LT.0.0D0) THEN
        WRITE(6, *) 'In VWKRGEN: best-fit search failed. IZ =',IZ
        WRITE(7, *) 'In VWKRGEN: best-fit search failed. IZ =',IZ
        STOP
      ENDIF
C
C**********************************************************************C
C     BEST-FIT ALPHA HAS NOW BEEN FOUND -- RESTORE THESE RESULTS.      C
C**********************************************************************C
C
C     BEST-FIT VALUES FOR ALPHA, BETA AND GAMMA
      AF = ALPH*(1.0D0 + 3.0D0*DFLOAT(IABG)/DFLOAT(NAF))
      BF = BETA
      GF = GAMA
C
C     STORE GEOMETRIC SET OF PARAMETERS IN XWKR
      XI = AF/(RNUC(IZ)**2)
      DO IFT=1,NFT
        XWKR(IZ,IFT) = XI
        RT = DFLOAT(IFT)/DFLOAT(NFT+1)
        XI = BF*XI*(1.0D0 + GF*RT*RT)
      ENDDO
C
C     TRANSFER MATCHING VALUES V(R) AND R*R*V(R) INTO Z VECTOR
      DO IFT=1,NFT
        Z(IFT) = Y(IFT)
      ENDDO
C
C     SET UP MATRIX EQUATIONS FOR SET OF FITTING GAUSSIANS
C
C     MATCH THE SECOND MOMENT
      DO JFT=1,NFT
        RAT = PI/XWKR(IZ,JFT)
        X(1,JFT) =-6.0D0*RAT*DSQRT(RAT)/(4.0D0*PI)
      ENDDO
C
C     MATCH THE FOURTH MOMENT
      DO JFT=1,NFT
        RAT = PI/XWKR(IZ,JFT)
        X(2,JFT) =-30.0D0*RAT*DSQRT(RAT)/(4.0D0*PI*XWKR(IZ,JFT))
      ENDDO
C
C     MATCH THE POTENTIAL ITSELF
      DO IFT=3,NFT
C
C       MATCHING RADIUS
        R = RM(IFT-2)
C
C       GAUSSIAN POTENTIALS AT THIS RADIUS
        DO JFT=1,NFT
          IF(IFT-2.LE.NDCT) THEN
            X(IFT,JFT) =     DEXP(-XWKR(IZ,JFT)*R*R)
          ELSE
            X(IFT,JFT) = R*R*DEXP(-XWKR(IZ,JFT)*R*R)
          ENDIF
        ENDDO
C
      ENDDO
C
C     SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
      CALL DGESV(NFT,1,X,NFT,IPIV,Z,NFT,INFO)
C
C     TRANSFER AMPLITUDES VALUES INTO FWKR ARRAY
      DO IFT=1,NFT
        FWKR(IZ,IFT) = Z(IFT)
      ENDDO
C
C**********************************************************************C
C     RESULTS: RADIAL MOMENTS                                          C
C**********************************************************************C
C
C     BEST-FIT NUCLEAR WICHMANN-KROLL RADIAL MOMENTS (IN HARTREE)
      DO K=-2,5
        BIN = 0.0D0
        DO IFT=1,NFT
          FR = FWKR(IZ,IFT)
          XI = XWKR(IZ,IFT)
          XR = XI**(K+1)
          BIN = BIN-FR*DFLOAT(K)*GAMHLF(K+3)/DSQRT(XR)
        ENDDO
        RK(K) = BIN
      ENDDO
C
C**********************************************************************C
C     RESULTS: POLARISATION CHARGE DENSITY                             C
C**********************************************************************C
C
C     VACUUM POLARISATION CHARGE DENSITY
      DO IPTS=0,NPTS
C
C       RADIUS R
        R = RU(IPTS)
C
C       INTIALISE CHARGE COUNTER
        Q = 0.0D0
        DO IFT=1,NFT
          FR = FWKR(IZ,IFT)
          XI = XWKR(IZ,IFT)
          Q  = Q + 2.0D0*FR*XI*R*R*(3.0D0-2.0D0*R*R*XI)*DEXP(-XI*R*R)
        ENDDO
        POL(IPTS) = Q
C
      ENDDO
C
C     LOCATION OF CHARGE ZERO (LINE OF BEST FIT BETWEEN TWO POINTS)
      RZER = 0.0D0
      DO IPTS=NPTS,0,-1
        IF(POL(IPTS).GT.0.0D0.AND.RU(IPTS).LT.20.0D0/CFM) THEN
          PD   = POL(IPTS)-POL(IPTS-1)
          RZER = (RU(IPTS-1)*POL(IPTS)-RU(IPTS)*POL(IPTS-1))/PD
          GOTO 51
        ENDIF
      ENDDO
51    CONTINUE
C
C     INTEGRATED POLARISED CHARGE UP TO THIS ZERO
      QPOL = 0.0D0
      DO IFT=1,NFT
        FR = FWKR(IZ,IFT)
        XI = XWKR(IZ,IFT)
        QPOL = QPOL - 2.0D0*(RZER**3)*XI*FR*DEXP(-XI*RZER*RZER)
      ENDDO
C
C**********************************************************************C
C     PRINT RESULTS IF PROMPTED                                        C
C**********************************************************************C
C
C     SKIP OUTPUT SECTION UNLESS PROMPTED
      IF(IWRT.EQ.0) GOTO 20
C
C     SOLUTION SET
35    FORMAT(1X,A,I2,A,I3,A,A,A)
36    FORMAT(1X,A,F17.15,1X,A)
37    FORMAT(1X,A,F14.12,A,F8.6,A,F8.6)
38    FORMAT(1X,A,I2,A,ES15.8,A,ES15.8,A)
39    FORMAT(1X,A,ES17.10,1X,A)
      WRITE(6, *) ''
      WRITE(7, *) ''
      WRITE(6,35) 'Wichmann-Kroll best-fit for centre IZ = ',IZ,
     &                  '  (',INT(ANUC(IZ)),'^',ELMT(INT(ZNUC(IZ))),')'
      WRITE(7,35) 'Wichmann-Kroll best-fit for centre IZ = ',IZ,
     &                  '  (',INT(ANUC(IZ)),'^',ELMT(INT(ZNUC(IZ))),')'
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
      WRITE(6,36) 'Here ξ0 = 1/R_n^2   and   R_n  = ',RNUC(IZ),'a0'
      WRITE(7,36) 'Here ξ0 = 1/R_n^2   and   R_n  = ',RNUC(IZ),'a0'
      WRITE(6,36) '                                 ',RNUC(IZ)*CFM,'fm'
      WRITE(7,36) '                                 ',RNUC(IZ)*CFM,'fm'
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
      WRITE(6,37) 'α = ',AF,' ξ0,   β = ',BF,',   γ = ',GF
      WRITE(7,37) 'α = ',AF,' ξ0,   β = ',BF,',   γ = ',GF
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
      WRITE(6,36) 'Least-squares best fit:       R^2 = ',RSQBIG
      WRITE(7,36) 'Least-squares best fit:       R^2 = ',RSQBIG
      WRITE(6,39) '1-R^2                             = ',1.0D0-RSQBIG
      WRITE(7,39) '1-R^2                             = ',1.0D0-RSQBIG
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
      DO IFT=1,NFT
        FR = FWKR(IZ,IFT)
        XI = RNUC(IZ)*RNUC(IZ)*XWKR(IZ,IFT)
        WRITE(6,38) 'f_',IFT,'(r): ',FR,'*exp(-',XI,' ξ0 r^2)'
        WRITE(7,38) 'f_',IFT,'(r): ',FR,'*exp(-',XI,' ξ0 r^2)'
      ENDDO
C
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
C
C     LOCATION OF CHARGE ZERO
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('-',50)
      WRITE(7, *) REPEAT('-',50)
      WRITE(6, *) 'Charge zero        :',RZER*CFM,' fm'
      WRITE(7, *) 'Charge zero        :',RZER*CFM,' fm'
      WRITE(6, *) 'Charge polarisation:',QPOL
      WRITE(7, *) 'Charge polarisation:',QPOL
      WRITE(6, *) REPEAT('-',50)
      WRITE(7, *) REPEAT('-',50)
C
C**********************************************************************C
C     BEST-FIT WICHMANN-KROLL POTENTIAL ON THE TESTING GRID.           C
C**********************************************************************C
C
C     THIS TAKES TIME -- SKIP UNLESS YOU WANT QUADRATURE/PLOTTING
c     GOTO 20
C
C     BEST-FIT GAUSSIAN POTENTIAL VALUES ACROSS UNIFORM GRID
      DO IPTS=0,NPTS
C
C       RADIUS ON UNIFORM SCALE
        R = RU(IPTS)
C
C       INITIALISE BEST-FIT POTENTIAL AT THIS RADIUS
        VG(IPTS) = 0.0D0
        DO IFT=1,NFT
          VG(IPTS) = VG(IPTS) + FWKR(IZ,IFT)*DEXP(-XWKR(IZ,IFT)*R*R)
        ENDDO
C
      ENDDO
C
C**********************************************************************C
C     EXACT WICHMANN-KROLL POTENTIAL (WEIGHTED BY R*R) ON GRID.        C
C**********************************************************************C
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
      DO N=0,NLIN
        RAD(N) = RORI + HL*DFLOAT(N)
      ENDDO
C
      DO N=0,NEXP
        RAD(N+NLIN) = RMID*DEXP(HE*DFLOAT(N))
      ENDDO
C
C     ORIGIN VALUE IS A SPECIAL CASE
      VVAC(IZ,0) = 0.0D0
C
C     POTENTIAL BY CONVOLUTION
      IF(IPOT.EQ.1) THEN
C
CC       THIS COMMENTED CODE TURNED OUT TO BE UNRELIABLY VOLATILE
C
C        ALPI = 1.0D0/(CV*PI)
C        ZALF = 1.0D0/CV
CC
CC       APPROXIMATE POTENTIAL EVERYWHERE
CC       AMP =-1.0D0/(PI*CV*CV*CV)
Cc       CALL VPOLPOT(CHIWKRF,CMPF,RAD,VTMP,AMP,RNUC(IZ),200.0D0/CFM)
C        DO N=1,NRAD
C          R = RAD(N)
Cc         VVAC(IZ,N) = R*R*VTMP(N)/CNG
C          VVAC(IZ,N) =-CV*ALPI*ZALF*ZALF*ZALF*R*VWKRFPT(R)
C          WRITE(*,*) N,R,VVAC(IZ,N),VWKRFPT(R)
C        ENDDO
C
C       FAR FROM THE NUCLEUS, JUST EVALUATE THE STRUCTURELESS POTENTIAL
        CALL VWKRFPT(IZ)
C
C       OVERRIDE WHEN R < 200 fm
        DO N=0,NRAD
C
C         RADIUS R AND INITIALISE COUNTER FOR POTENTIAL
          R = RAD(N)
C
C         NUMERICAL INVERSE HANKEL TRANSFORM (FOR A POTENTIAL)
          IF(R.LE.200.0D0/CFM) THEN
            CALL VRIHNKT(GWKR,R,V)
            VVAC(IZ,N) = R*R*V
          ENDIF
C
        ENDDO
C
      ENDIF
C
C**********************************************************************C
C     PLOTTING SECTION                                                 C
C**********************************************************************C
C
      IF(IWRT.NE.1) GOTO 20
C
C     NAME TAG FOR PDF DOCUMENTS
      IF(INT(ANUC(IZ)).LT.10) THEN
        WRITE(ATRM,'(A,I1)') '00',INT(ANUC(IZ))
      ELSEIF(INT(ANUC(IZ)).LT.100) THEN
        WRITE(ATRM,'(A,I2)') '0',INT(ANUC(IZ))
      ELSE
        WRITE(ATRM,'(I3)') INT(ANUC(IZ))
      ENDIF
C
C     DETAILS COMMON TO ALL PLOTS
      XAXIS  = 'r/RNUC'
      KEY(1) = 'Exact Wichmann-Kroll'
      KEY(2) = 'Best fit'
C
C     UNIFORM GRID VALUES R*R*V(R)
      XOUT   = 'Wichmann-Kroll-uniform'//TRIM(ATRM)
      TITLE  = 'Wichmann-Kroll r^{2}*V(r) on uniform grid'
      YAXIS  = 'r^{2}*V(r)/Z'
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
      DO IPTS=0,NPTS
        R = RU(IPTS)
        WRITE(8, *) R/RNUC(IZ),R*R*VF(IPTS),R*R*VG(IPTS)
      ENDDO
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C     NUMERICAL EVALUATION FOR V(R)
      XOUT   = 'Wichmann-Kroll-pointwise'//TRIM(ATRM)
      TITLE  = 'Wichmann-Kroll V(r) on piecewise grid'
      YAXIS  = 'V(r)/Z'
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
      DO N=0,NRAD
        R = RAD(N)
        IF(R/RNUC(IZ).GT.5.0D0) GOTO 123
C       TRUE WICHMANN-KROLL POTENTIAL
        IF(N.EQ.0) THEN
          VTR = V0
        ELSE
          VTR = VVAC(IZ,N)/(R*R)
        ENDIF
C       APPROXIMATE VALUES
        V = 0.0D0
        DO IFT=1,NFT
          V = V + FWKR(IZ,IFT)*DEXP(-XWKR(IZ,IFT)*R*R)
        ENDDO
        WRITE(8, *) R*cfm,VTR,V
      ENDDO
123   CONTINUE
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C     UNIFORM GRID VALUES R*R*V(R)
      XOUT   = 'Wichmann-Kroll-charge'//TRIM(ATRM)
      TITLE  = 'Wichmann-Kroll rho(r) on uniform grid'
      YAXIS  = 'rho(r)/Z'
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
      DO IPTS=0,NPTS/25
        R = RU(IPTS)
C       IF(R*CFM.LE.16.0D0) THEN
C         NUMERICAL INVERSE HANKEL TRANSFORM (FOR A CHARGE DENSITY)
          CALL PRIHNKT(GWKR,R,RHOWKR)
          WRITE(8, *) R/rnuc(iz),-4.0D0*R*R*PI*RHOWKR,POL(IPTS)
C       ENDIF
      ENDDO
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
20    CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE VKSBGEN(IZ,NFT,IWRT,RV,RSQBIG,RZER,QPOL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    VV    VV KK    KK  SSSSSS  BBBBBBB   GGGGGG  EEEEEEEE NN    NN    C
C    VV    VV KK   KK  SS    SS BB    BB GG    GG EE       NNN   NN    C
C    VV    VV KK  KK   SS       BB    BB GG       EE       NNNN  NN    C
C    VV    VV KKKKK     SSSSSS  BBBBBBB  GG       EEEEEE   NN NN NN    C
C     VV  VV  KK  KK         SS BB    BB GG   GGG EE       NN  NNNN    C
C      VVVV   KK   KK  SS    SS BB    BB GG    GG EE       NN   NNN    C
C       VV    KK    KK  SSSSSS  BBBBBBB   GGGGGG  EEEEEEEE NN    NN    C
C                                                                      C
C -------------------------------------------------------------------- C
C  VKSBGEN GENERATES A BEST-FIT GAUSSIAN SET FOR THE KÄLLÉN-SABRY      C
C  POTENTIAL ARISING FROM NUCLEUS IZ, USING NFT TOTAL GAUSSIANS.       C
C -------------------------------------------------------------------- C
C  MATCHING CRITERIA AVAILABLE FOR GAUSSIAN AMPLITUDES:                C
C   ▶ RADIAL ARGUMENTS SCALED BY NUCLEAR RMS RADIUS RNUC(IZ).          C
C   ▶ KÄLLÉN-SABRY POTENTIAL FIRST SCALED TO R*R*V(R).                 C
C   ▶ ZEROTH MOMENT SATISFIED AUTOMATICALLY BUT SECOND MOMENT ENFORCED.C
C   ▶ POTENTIAL MATCHED TO NFT POINTS, FROM 3.0D0*RNUC TO 200.0D0*RNUC.C
C   ▶ MATCHING POINTS ARE EXPONENTIALLY SPACED.                        C
C   ▶ GAUSSIAN EXPONENTS FROM MODIFIED GEOMETRIC SERIES WITH THE SAME  C
C     BETA AND GAMMA, BUT ALPHA IS SCALED BY RNUC*RNUC.                C
C   ▶ WHEN THIS IS COMPLETE, POTENTIAL V(R) NEAR R=0.0D0 IS AUGMENTED  C
C     WITH SOME MORE GAUSSIANS WITH LARGER EXPONENTS IF V(R) IS NEEDED.C
C -------------------------------------------------------------------- C
C INPUT:                                                               C
C   ▶ IZ:   NUCLEAR CENTRE.                                            C
C   ▶ NFT:  DESIRED NUMBER OF FITTING FUNCTIONS IN THE SET.            C
C   ▶ IWRT: EXTENDED WRITTEN OUTPUT OPTION.                            C
C   ▶ RV:   MASS MODIFIER (1.0D0 FOR ELECTRONS, MASS RATIO OTHERWISE). C
C**********************************************************************C
      INCLUDE 'parameters.h'
      PARAMETER(NPTS=6000)
C
      EXTERNAL CHIKSBF
C
      LOGICAL FILETHERE
C
      CHARACTER*2  ELMT(120)
      CHARACTER*3  ATRM
      CHARACTER*5  NMDL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(2)
C
      DIMENSION VF(0:NPTS),VG(0:NPTS),RU(0:NPTS)
      DIMENSION X(NFT,NFT),Z(NFT),Y(NFT),RM(NFT),IPIV(NFT)
      DIMENSION RHO(0:NSRC),RHOTMP(0:NRVP)
      DIMENSION POL(0:NPTS),RK(-2:5)
      DIMENSION FMQ(0:NQVP),PKSB(0:NQVP),GKSB(0:NQVP)
      DIMENSION VTMP(0:NRAD)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BKSB/RKSB(MCT,3),FKSB(MCT,MFT),XKSB(MCT,MFT),NKSB(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/MDLV/ELMT
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/QGRD/QFM(0:NQVP),QEP(0:NQVP),QMIN,QMID,QMAX
      COMMON/RGRD/RFM(0:NRVP),REP(0:NRVP),RMIN,RMD,RMX
      COMMON/QKSB/RAD(0:NRAD),VVAC(MCT,0:NRAD),RORI,RMID,RMAX,NLIN,NEXP
C
C     MULTIPLIER FOR COMPTON WAVELENGTH (USE FOR MUON OR TAUON FIELD)
      CMPF = CMPW/(RV)
C
C     CHECK THAT THERE ARE ENOUGH FITTING FUNCTIONS
      IF(NFT.LT.15) THEN
        WRITE(6, *) 'In VKSBGEN: need more fitting functions. NFT =',NFT
        WRITE(7, *) 'In VKSBGEN: need more fitting functions. NFT =',NFT
        STOP
      ELSEIF(NFT+1.GT.MFT) THEN
        WRITE(6, *) 'In VKSBGEN: too many fitting functions. NFT =',NFT
        WRITE(7, *) 'In VKSBGEN: too many fitting functions. NFT =',NFT
        STOP
      ENDIF
C
C     NUMBER OF FITTING FUNCTIONS
      NKSB(IZ) = NFT
C
C     ZEROTH, SECOND AND FOURTH MOMENTS OF POLARISED CHARGE DENSITY
      FAC = CMPF*CMPF/(PI*PI*CV*CV)
      FNT = RNUC(IZ)/CMPF
      RKSB(IZ,1) = 0.0D0
      RKSB(IZ,2) = 41.0D0*FAC/27.0D0
      RKSB(IZ,3) = 401.0D0*FAC*CMPF*CMPF/90.0D0
      RKSB(IZ,3) = RKSB(IZ,3)*(1.0D0 + 4100.0D0*FNT*FNT/3609.0D0)
C
C     NUMERICAL CALCULATION OF THE KÄLLÉN-SABRY POTENTIAL
C     1: FOURIER, 2: ------------------, 3: F+R CHEBYSHEV
      IPOT = 1
C
C**********************************************************************C
C     GENERATE POLARISATION FORM FACTOR                                C
C**********************************************************************C
C
      IF(IPOT.NE.1) GOTO 100
C
C     READ SPECTRAL FUNCTION FROM FILE OR RECALCULATE
      INQUIRE(FILE='spectral/pi_ksb_electron.dat',EXIST=FILETHERE)
      IF(FILETHERE) THEN
        OPEN(UNIT=8,FILE='spectral/pi_ksb_electron.dat',
     &                                               STATUS='UNKNOWN')
        REWIND(UNIT=8)
        DO I=0,NQVP
          READ(8, *) M,Q,PKSB(I)
          IF(M.NE.I.OR.Q.NE.QFM(I)) THEN
          	FILETHERE = .FALSE.
            GOTO 111
          ENDIF
        ENDDO
111     CONTINUE
        CLOSE(UNIT=8)
      ENDIF
      IF(.NOT.FILETHERE) THEN
        OPEN(UNIT=8,FILE='spectral/pi_ksb_electron.dat',
     &                                               STATUS='UNKNOWN')
        REWIND(UNIT=8)
        DO I=0,NQVP
          PKSB(I) = PIKSBF(QFM(I),CMPF)
          WRITE(8, *) I,QFM(I),PKSB(I)
        ENDDO
        CLOSE(UNIT=8)
      ENDIF
110   CONTINUE
CC
CC     NUCLEAR FORM FACTOR
C      DO I=0,NQVP
C        FMQ(I) = 0.0D0
C        Q = QFM(I)
C        FAC = 4.0D0*XNUC(IZ,0)
C        FMQ(I) = FMQ(I) + FNUC(IZ,0)*DEXP(-Q*Q/FAC)
C        DO K=1,NNUC(IZ)
C          FAC = 4.0D0*XNUC(IZ,K)
C          X32 = XNUC(IZ,K)*DSQRT(XNUC(IZ,K))
C          PRE = 0.25D0*PI12*Q*Q/X32
C          FMQ(I) = FMQ(I) - PRE*FNUC(IZ,K)*DEXP(-Q*Q/FAC)
C        ENDDO
C      ENDDO
CC
C     NUCLEAR CHARGE DISTRIBUTION ALONG RADIAL GRID
      DO N=0,NRVP
        R = RFM(N)
        RHOTMP(N) = RHONUC(NMDL(IZ),IZ,R)
      ENDDO

C     CORRESPONDING NUCLEAR FORM FACTOR
      DO I=0,NQVP
        Q = QFM(I)
        DO N=0,NRVP
          R = RFM(N)
          FMQ(I) = FMQ(I)
     &           + EXTINT11(R*R*RHOTMP(N)*SPHBJ0(Q*R)*REP(N),N,NRVP)
        ENDDO
        FMQ(I) = 5.0D0*FMQ(I)/2.99376D+5
        FMQ(I) = 4.0D0*PI*FMQ(I)
      ENDDO
C
C     CALCULATE THE POLARISATION FORM FACTOR AS A PRODUCT
      DO I=0,NQVP
        GKSB(I) = FMQ(I)*PKSB(I)
      ENDDO
C
C     EXPLICITLY CALCULATE V0 (SPECIAL CASE)
      V0 = 0.0D0
      DO I=0,NQVP
        Q = QFM(I)
        V0 = V0 + EXTINT11(QEP(I)*GKSB(I),I,NQVP)
      ENDDO
      V0 =-2.0D0*5.0D0*V0/(PI*2.99376D+5)
C
100   CONTINUE
C
C**********************************************************************C
C     INTEGRATION PARAMETERS AND UEHLING POTENTIAL AT THE ORIGIN       C
C**********************************************************************C
C
C     SEARCH FOR CHARGE RADIUS SMAX FOR WHICH RHO(SMAX) < 1.0D-16
      SMAX = 0.0D0
60    SMAX = SMAX + 0.1D0/CFM
      PMAX = SMAX*RHONUC(NMDL(IZ),IZ,SMAX)
      IF(DABS(PMAX).GT.1.0D-16) GOTO 60
C
C     SOURCE CHARGE STEP SIZE
      HS = SMAX/DFLOAT(NSRC)
C
C     EVALUATE CHARGE DENSITY ON UNIFORMLY-SPACED GRID
      DO M=0,NSRC
        S = HS*DFLOAT(M)
        RHO(M) = RHONUC(NMDL(IZ),IZ,S)
      ENDDO
C
C     POTENTIAL AT THE ORIGIN USING CHEBYSHEV
      IF(IPOT.EQ.3) THEN
C
C       INTEGRATE OVER CHARGE SOURCE
        V0 = 0.0D0
        DO M=0,NSRC
C
C         SET SOURCE RADIUS
          S = HS*DFLOAT(M)
C
C         FUNCTION ARGUMENTS
          XS = 2.0D0*DABS(S)/CMPF
C
C         CHEBYSHEV APPROXIMATION FIT (FULLERTON+RINKER)
          CS = FUNL(XS,1)
C
C         CONTRIBUTION TO INTEGRAND
          V0 = V0 + EXTINT11(S*RHO(M)*CS,M,NSRC)
C
        ENDDO
        V0 = 5.0D0*HS*V0/2.99376D+5
        V0 =-4.0D0*V0/(PI*CV*CV)
C
      ENDIF
C
C**********************************************************************C
C     KÄLLÉN-SABRY POTENTIAL VALUES ON THE TESTING GRID.               C
C**********************************************************************C
C
C     KÄLLÉN-SABRY POTENTIAL SAMPLING VALUES (FOR R-SQUARED CALCULATION)
      RLIM = 250.0D0
C
C     STORE R*R*V(R) ON THE UNIFORM GRID IN VF(0:NPTS)
      RU(0) = 0.0D0
      VF(0) = 0.0D0
C
C     POTENTIAL BY CONVOLUTION
      IF(IPOT.EQ.1) THEN
C
        DO IPTS=1,NPTS
C
C         RADIUS ON UNIFORM SCALE
          R = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C         STORE FOR LATER
          RU(IPTS) = R
C
C         NUMERICAL INVERSE HANKEL TRANSFORM (FOR A POTENTIAL)
          CALL VRIHNKT(GKSB,R,VF(IPTS))
C
        ENDDO
C
C     FULLERTON/RINKER METHOD
      ELSEIF(IPOT.EQ.3) THEN
C
        ISWITCH = 1
C
        DO IPTS=1,NPTS
C
C         RADIUS ON UNIFORM SCALE
          R = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
C
C         STORE FOR LATER
          RU(IPTS) = R
C
C         EVALUATE V AS IF IN ASYMPTOTIC REGION
          IF(R.NE.0.0D0) THEN
            V = FUNL(2.0D0*R/CMPF,1)/(PI*CMPF)
          ENDIF
C
C         EVALUATE V AS IF IN NUCLEAR REGION
          IF(ISWITCH.EQ.0) GOTO 202
C
C         INTIALISE POTENTIAL COUNTER
          VN = 0.0D0
C
C         INTEGRATE OVER CHARGE SOURCE
          DO M=0,NSRC
C
C           SET SOURCE RADIUS
            S = HS*DFLOAT(M)
C
C           FUNCTION ARGUMENTS
            XM = 2.0D0*DABS(R-S)/CMPF
            XP = 2.0D0*DABS(R+S)/CMPF
C
C           CHEBYSHEV APPROXIMATION FIT (FULLERTON+RINKER)
            CM = FUNL(XM,0)
            CP = FUNL(XP,0)
C
C           CONTRIBUTION TO INTEGRAND
            VN = VN + EXTINT11(S*RHO(M)*(CM-CP),M,NSRC)
C
          ENDDO
C
C         INTEGRATION WEIGHTING FACTORS
          VN = 5.0D0*HS*VN/2.99376D+5
C
C         IF FULL AND ASYMPTOTIC RESULTS ARE SIMILAR, FLIP SWITCH
          IF(DABS((VN-V)/VN)*ZNUC(IZ).LT.1.0D-5) THEN
            ISWITCH = 0
          ENDIF
C
C         OVERRIDE POTENTIAL NOW
          V = VN
C
202       CONTINUE
C
C         OTHER FACTORS AND FINAL VALUE V(R)
          VF(IPTS) =-CMPF*V/(PI*CV*CV*R)
C
        ENDDO
C      
      ENDIF
C
C**********************************************************************C
C     KÄLLÉN-SABRY POTENTIAL VALUES ON THE MATCHING GRID.              C
C**********************************************************************C
C
C     FIRST CONDITION MATCHES RMS RADIUS
      Y(1) = RKSB(IZ,2)
C
C     SECOND CONDITION MATCHES FOURTH RADIAL MOMENT
      Y(2) = RKSB(IZ,3)
C
C     DIVIDE NFT POINTS INTO DIRECT AND WEIGHTED REGIONS
C     NDCT = (NFT)/4
      NDCT = 6
      NWGT = NFT-NDCT-2
C
C     KÄLLÉN-SABRY POTENTIAL MATCHING VALUES
      R0 = 0.0D0
      RC = 2.0D0
C     RC = 3.0D0
      RN = 240.0D0
C     RN = 200.0D0
      HD = (RC-R0)/DFLOAT(NDCT)
      HW = DLOG(RN/RC)/DFLOAT(NWGT-1)
C
C     RADIAL MATCHING POINTS
      DO IDCT=1,NDCT
        RM(IDCT     ) = RNUC(IZ)*(R0 + DFLOAT(IDCT-1)*HD)
      ENDDO
      DO IWGT=1,NWGT
        RM(IWGT+NDCT) = RNUC(IZ)*RC*DEXP(DFLOAT(IWGT-1)*HW)
      ENDDO
C
C     STORE R*R*V(R) MATCHING VALUES IN Y(NFT) -- MATRIX EQN SOLUTIONS
      Y(3) = V0
C
      DO IFT=4,NFT
C
C       RADIAL MATCHING POINT
        R = RM(IFT-2)
C
C       INITIALISE VALUE FOR THE POTENTIAL HERE
        Y(IFT) = 0.0D0
C
C       POTENTIAL BY CONVOLUTION
        IF(IPOT.EQ.1) THEN
C
C         NUMERICAL INVERSE HANKEL TRANSFORM (FOR A POTENTIAL)
          CALL VRIHNKT(GKSB,R,Y(IFT))
C
C       FULLERTON/RINKER METHOD
        ELSEIF(IPOT.EQ.3) THEN
C
          ISWITCH = 1
C
C         EVALUATE V AS IF IN ASYMPTOTIC REGION
          IF(R.NE.0.0D0) THEN
            V = FUNL(2.0D0*R/CMPF,1)/(PI*CMPF)
          ENDIF
C
C         EVALUATE V AS IF IN NUCLEAR REGION
          IF(ISWITCH.EQ.0) GOTO 212
C
C         INTIALISE POTENTIAL COUNTER
          VN = 0.0D0
C
C         INTEGRATE OVER CHARGE SOURCE
          DO M=0,NSRC
C
C           SET SOURCE RADIUS
            S = HS*DFLOAT(M)
C
C           FUNCTION ARGUMENTS
            XM = 2.0D0*DABS(R-S)/CMPF
            XP = 2.0D0*DABS(R+S)/CMPF
C
C           CHEBYSHEV APPROXIMATION FIT (FULLERTON+RINKER)
            CM = FUNL(XM,0)
            CP = FUNL(XP,0)
C
C           CONTRIBUTION TO INTEGRAND
            VN = VN + EXTINT11(S*RHO(M)*(CM-CP),M,NSRC)
C
          ENDDO
C
C         INTEGRATION WEIGHTING FACTORS
          VN = 5.0D0*HS*VN/2.99376D+5
C
C         IF FULL AND ASYMPTOTIC RESULTS ARE SIMILAR, FLIP SWITCH
          IF(DABS((VN-V)/VN)*ZNUC(IZ).LT.1.0D-5) THEN
            ISWITCH = 0
          ENDIF
C
C         OVERRIDE POTENTIAL NOW
          V = VN
C
212       CONTINUE
C
C         OTHER FACTORS AND FINAL VALUE V(R)
          Y(IFT) =-CMPF*V/(PI*CV*CV*R)
C
        ENDIF
C
C       DATA POINTS IN THE WEIGHTED REGION MUST BE MULTIPLIED BY R*R
        IF(IFT-2.GT.NDCT) THEN
          Y(IFT) = R*R*Y(IFT)
        ENDIF
C
      ENDDO
C
C**********************************************************************C
C     GIVEN MODIFIED GEOMETRIC SET OF NFT PARAMETERS, FIND BEST ALPHA. C
C**********************************************************************C
C
C     GEOMETRIC PARAMETER SEARCH DETAILS
      ALPH = 2.00D-5
      BETA = 1.60D+0
      GAMA = 0.21D+0
C
C     NUMBER OF INCREMENTS IN SEARCH
      NAF = 1200
C
C     INITIALISE BEST-FIT R-SQUARED AND INDEX FOR ALPHA
      RSQBIG = 0.0D0
      IABG   = 0
C
C     ITERATE OVER ALPHA SEARCH
      DO IA=0,NAF
C
C       TRIAL ALPHA PARAMETER AND FROZEN BETA PARAMETER
        AF = ALPH*(1.0D0 + 3.0D0*DFLOAT(IA)/DFLOAT(NAF))/(RNUC(IZ)**2)
        BF = BETA
C
C       STORE GEOMETRIC SET OF PARAMETERS IN XKSB
        XI = AF
        DO IFT=1,NFT
          XKSB(IZ,IFT) = XI
          RT = DFLOAT(IFT)/DFLOAT(NFT+1)
          XI = BF*XI*(1.0D0 + GAMA*RT*RT)
        ENDDO
C
C       TRANSFER MATCHING VALUES V(R) AND R*R*V(R) INTO Z VECTOR
        DO IFT=1,NFT
          Z(IFT) = Y(IFT)
        ENDDO
C
C       SET UP MATRIX EQUATIONS FOR SET OF FITTING GAUSSIANS
C
C       MATCH THE SECOND MOMENT
        DO JFT=1,NFT
          RAT = PI/XKSB(IZ,JFT)
          X(1,JFT) =-6.0D0*RAT*DSQRT(RAT)/(4.0D0*PI)
        ENDDO
C
C       MATCH THE FOURTH MOMENT
        DO JFT=1,NFT
          RAT = PI/XKSB(IZ,JFT)
          X(2,JFT) =-30.0D0*RAT*DSQRT(RAT)/(4.0D0*PI*XKSB(IZ,JFT))
        ENDDO
C
C       MATCH THE POTENTIAL ITSELF
        DO IFT=3,NFT
C
C         MATCHING RADIUS
          R = RM(IFT-2)
C
C         GAUSSIAN POTENTIALS AT THIS RADIUS
          DO JFT=1,NFT
            IF(IFT-2.LE.NDCT) THEN
              X(IFT,JFT) =     DEXP(-XKSB(IZ,JFT)*R*R)
            ELSE
              X(IFT,JFT) = R*R*DEXP(-XKSB(IZ,JFT)*R*R)
            ENDIF
          ENDDO
C
        ENDDO
C
C       SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
        CALL DGESV(NFT,1,X,NFT,IPIV,Z,NFT,INFO)
C
C       TRANSFER AMPLITUDES VALUES INTO FKSB ARRAY
        DO IFT=1,NFT
          FKSB(IZ,IFT) = Z(IFT)
        ENDDO
C
C       BEST-FIT GAUSSIAN POTENTIAL VALUES ACROSS UNIFORM GRID
        DO IPTS=0,NPTS
C
C         RADIUS ON UNIFORM SCALE
          R = RU(IPTS)
C
C         INITIALISE BEST-FIT POTENTIAL AT THIS RADIUS
          VG(IPTS) = 0.0D0
          DO IFT=1,NFT
            VG(IPTS) = VG(IPTS) + FKSB(IZ,IFT)*DEXP(-XKSB(IZ,IFT)*R*R)
          ENDDO
C
        ENDDO
C
C       AVERAGE POTENTIAL VALUE
        YB = 0.0D0
        DO IPTS=0,NPTS
          YB = YB + VG(IPTS)
        ENDDO
        YB = YB/DFLOAT(NPTS+1)
C
C       PREPARATION FOR R-SQUARED VALUE
        SRES = 0.0D0
        STOT = 0.0D0
        DO IPTS=0,NPTS
          SRES = SRES + RU(IPTS)*RU(IPTS)*(VG(IPTS)-VF(IPTS))**2
          STOT = STOT + RU(IPTS)*RU(IPTS)*(VG(IPTS)-YB      )**2
        ENDDO
C
C       R-SQUARED VALUE AND BEST-FIT CHECK
        RSQ = 1.0D0 - SRES/STOT
        IF(RSQ.GT.RSQBIG) THEN
          IABG   = IA
          RSQBIG = RSQ
        ENDIF
C
      ENDDO
C
C     IF NO SOLUTION CAN BE FOUND, ALERT THE USER
      IF(RSQBIG.LT.0.0D0) THEN
        WRITE(6, *) 'In VKSBGEN: best-fit search failed. IZ =',IZ
        WRITE(7, *) 'In VKSBGEN: best-fit search failed. IZ =',IZ
        STOP
      ENDIF
C
C**********************************************************************C
C     BEST-FIT ALPHA HAS NOW BEEN FOUND -- RESTORE THESE RESULTS.      C
C**********************************************************************C
C
C     BEST-FIT VALUES FOR ALPHA, BETA AND GAMMA
      AF = ALPH*(1.0D0 + 3.0D0*DFLOAT(IABG)/DFLOAT(NAF))
      BF = BETA
      GF = GAMA
C
C     STORE GEOMETRIC SET OF PARAMETERS IN XKSB
      XI = AF/(RNUC(IZ)**2)
      DO IFT=1,NFT
        XKSB(IZ,IFT) = XI
        RT = DFLOAT(IFT)/DFLOAT(NFT+1)
        XI = BF*XI*(1.0D0 + GF*RT*RT)
      ENDDO
C
C     TRANSFER MATCHING VALUES V(R) AND R*R*V(R) INTO Z VECTOR
      DO IFT=1,NFT
        Z(IFT) = Y(IFT)
      ENDDO
C
C     SET UP MATRIX EQUATIONS FOR SET OF FITTING GAUSSIANS
C
C     MATCH THE SECOND MOMENT
      DO JFT=1,NFT
        RAT = PI/XKSB(IZ,JFT)
        X(1,JFT) =-6.0D0*RAT*DSQRT(RAT)/(4.0D0*PI)
      ENDDO
C
C     MATCH THE FOURTH MOMENT
      DO JFT=1,NFT
        RAT = PI/XKSB(IZ,JFT)
        X(2,JFT) =-30.0D0*RAT*DSQRT(RAT)/(4.0D0*PI*XKSB(IZ,JFT))
      ENDDO
C
C     MATCH THE POTENTIAL ITSELF
      DO IFT=3,NFT
C
C       MATCHING RADIUS
        R = RM(IFT-2)
C
C       GAUSSIAN POTENTIALS AT THIS RADIUS
        DO JFT=1,NFT
          IF(IFT-2.LE.NDCT) THEN
            X(IFT,JFT) =     DEXP(-XKSB(IZ,JFT)*R*R)
          ELSE
            X(IFT,JFT) = R*R*DEXP(-XKSB(IZ,JFT)*R*R)
          ENDIF
        ENDDO
C
      ENDDO
C
C     SOLVE THE MATRIX EQUATION X.A = Z FOR AMLPITUDES A
      CALL DGESV(NFT,1,X,NFT,IPIV,Z,NFT,INFO)
C
C     TRANSFER AMPLITUDES VALUES INTO FKSB ARRAY
      DO IFT=1,NFT
        FKSB(IZ,IFT) = Z(IFT)
      ENDDO
C
C**********************************************************************C
C     RESULTS: RADIAL MOMENTS                                          C
C**********************************************************************C
C
C     BEST-FIT NUCLEAR KÄLLÉN-SABRY RADIAL MOMENTS (IN HARTREE)
      DO K=-2,5
        BIN = 0.0D0
        DO IFT=1,NFT
          FR = FKSB(IZ,IFT)
          XI = XKSB(IZ,IFT)
          XR = XI**(K+1)
          BIN = BIN-FR*DFLOAT(K)*GAMHLF(K+3)/DSQRT(XR)
        ENDDO
        RK(K) = BIN
      ENDDO
C
C**********************************************************************C
C     RESULTS: POLARISATION CHARGE DENSITY                             C
C**********************************************************************C
C
C     VACUUM POLARISATION CHARGE DENSITY
      DO IPTS=0,NPTS
C
C       RADIUS R
        R = RU(IPTS)
C
C       INTIALISE CHARGE COUNTER
        Q = 0.0D0
        DO IFT=1,NFT
          FR = FKSB(IZ,IFT)
          XI = XKSB(IZ,IFT)
          Q  = Q + 2.0D0*FR*XI*R*R*(3.0D0-2.0D0*R*R*XI)*DEXP(-XI*R*R)
        ENDDO
        POL(IPTS) = Q
C
      ENDDO
C
C     LOCATION OF CHARGE ZERO (LINE OF BEST FIT BETWEEN TWO POINTS)
      RZER = 0.0D0
      DO IPTS=NPTS,0,-1
        IF(POL(IPTS).LT.0.0D0.AND.RU(IPTS).LT.10.0D0/CFM) THEN
          PD   = POL(IPTS)-POL(IPTS-1)
          RZER = (RU(IPTS-1)*POL(IPTS)-RU(IPTS)*POL(IPTS-1))/PD
          GOTO 51
        ENDIF
      ENDDO
51    CONTINUE
C
C     INTEGRATED POLARISED CHARGE UP TO THIS ZERO
      QPOL = 0.0D0
      DO IFT=1,NFT
        FR = FKSB(IZ,IFT)
        XI = XKSB(IZ,IFT)
        QPOL = QPOL - 2.0D0*(RZER**3)*XI*FR*DEXP(-XI*RZER*RZER)
      ENDDO
C
C**********************************************************************C
C     PRINT RESULTS IF PROMPTED                                        C
C**********************************************************************C
C
C     SKIP OUTPUT SECTION UNLESS PROMPTED
      IF(IWRT.EQ.0) GOTO 20
C
C     SOLUTION SET
35    FORMAT(1X,A,I2,A,I3,A,A,A)
36    FORMAT(1X,A,F17.15,1X,A)
37    FORMAT(1X,A,F14.12,A,F8.6,A,F8.6)
38    FORMAT(1X,A,I2,A,ES15.8,A,ES15.8,A)
39    FORMAT(1X,A,ES17.10,1X,A)
      WRITE(6, *) ''
      WRITE(7, *) ''
      WRITE(6,35) 'Källén-Sabry best-fit for centre IZ = ',IZ,
     &                  '  (',INT(ANUC(IZ)),'^',ELMT(INT(ZNUC(IZ))),')'
      WRITE(7,35) 'Källén-Sabry best-fit for centre IZ = ',IZ,
     &                  '  (',INT(ANUC(IZ)),'^',ELMT(INT(ZNUC(IZ))),')'
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
      WRITE(6,36) 'Here ξ0 = 1/R_n^2   and   R_n  = ',RNUC(IZ),'a0'
      WRITE(7,36) 'Here ξ0 = 1/R_n^2   and   R_n  = ',RNUC(IZ),'a0'
      WRITE(6,36) '                                 ',RNUC(IZ)*CFM,'fm'
      WRITE(7,36) '                                 ',RNUC(IZ)*CFM,'fm'
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
      WRITE(6,37) 'α = ',AF,' ξ0,   β = ',BF,',   γ = ',GF
      WRITE(7,37) 'α = ',AF,' ξ0,   β = ',BF,',   γ = ',GF
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
      WRITE(6,36) 'Least-squares best fit:       R^2 = ',RSQBIG
      WRITE(7,36) 'Least-squares best fit:       R^2 = ',RSQBIG
      WRITE(6,39) '1-R^2                             = ',1.0D0-RSQBIG
      WRITE(7,39) '1-R^2                             = ',1.0D0-RSQBIG
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
      DO IFT=1,NFT
        FR = FKSB(IZ,IFT)
        XI = RNUC(IZ)*RNUC(IZ)*XKSB(IZ,IFT)
        WRITE(6,38) 'f_',IFT,'(r): ',FR,'*exp(-',XI,' ξ0 r^2)'
        WRITE(7,38) 'f_',IFT,'(r): ',FR,'*exp(-',XI,' ξ0 r^2)'
      ENDDO
C
      WRITE(6, *) REPEAT('-',53)
      WRITE(7, *) REPEAT('-',53)
C
C     LOCATION OF CHARGE ZERO
      WRITE(6, *) ' '
      WRITE(7, *) ' '
      WRITE(6, *) REPEAT('-',50)
      WRITE(7, *) REPEAT('-',50)
      WRITE(6, *) 'Charge zero        :',RZER*CFM,' fm'
      WRITE(7, *) 'Charge zero        :',RZER*CFM,' fm'
      WRITE(6, *) 'Charge polarisation:',QPOL
      WRITE(7, *) 'Charge polarisation:',QPOL
      WRITE(6, *) REPEAT('-',50)
      WRITE(7, *) REPEAT('-',50)
C
C**********************************************************************C
C     BEST-FIT KÄLLÉN-SABRY POTENTIAL ON THE TESTING GRID.             C
C**********************************************************************C
C
C     THIS TAKES TIME -- SKIP UNLESS YOU WANT QUADRATURE/PLOTTING
c     GOTO 20
C
C     BEST-FIT GAUSSIAN POTENTIAL VALUES ACROSS UNIFORM GRID
      DO IPTS=0,NPTS
C
C       RADIUS ON UNIFORM SCALE
        R = RU(IPTS)
C
C       INITIALISE BEST-FIT POTENTIAL AT THIS RADIUS
        VG(IPTS) = 0.0D0
        DO IFT=1,NFT
          VG(IPTS) = VG(IPTS) + FKSB(IZ,IFT)*DEXP(-XKSB(IZ,IFT)*R*R)
        ENDDO
C
      ENDDO
C
C**********************************************************************C
C     EXACT KÄLLÉN-SABRY POTENTIAL (WEIGHTED BY R*R) ON COMPOSITE GRID.C
C**********************************************************************C
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
      DO N=0,NLIN
        RAD(N) = RORI + HL*DFLOAT(N)
      ENDDO
C
      DO N=0,NEXP
        RAD(N+NLIN) = RMID*DEXP(HE*DFLOAT(N))
      ENDDO
C
C     ORIGIN VALUE IS A SPECIAL CASE
      VVAC(IZ,0) = 0.0D0
C
C     POTENTIAL BY CONVOLUTION
      IF(IPOT.EQ.1) THEN
C
C       APPROXIMATE POTENTIAL EVERYWHERE
        AMP =-1.0D0/(PI*PI*CV*CV)
        CALL VPOLPOT(CHIKSBF,CMPF,RAD,VTMP,AMP,RNUC(IZ),200.0D0/CFM)
        DO N=1,NRAD
          R = RAD(N)
          VVAC(IZ,N) = R*R*VTMP(N)
        ENDDO
C
C       OVERRIDE WHEN R < 200 fm
        DO N=1,NRAD
C
C         RADIUS R AND INITIALISE COUNTER FOR POTENTIAL
          R = RAD(N)
C
C         NUMERICAL INVERSE HANKEL TRANSFORM (FOR A POTENTIAL)
          IF(R.LE.200.0D0/CFM) THEN
            CALL VRIHNKT(GKSB,R,V)
            VVAC(IZ,N) = R*R*V
          ENDIF
C
        ENDDO
C
C     FULLERTON/RINKER METHOD
      ELSEIF(IPOT.EQ.3) THEN
C
        ISWITCH = 1
C
C       VALUE OF R*R*V(R) AT EACH OF THESE RADII
        DO N=1,NRAD
C
C         RADIUS R AND INITIALISE COUNTER FOR POTENTIAL
          R = RAD(N)
C
C         EVALUATE V AS IF IN ASYMPTOTIC REGION
          IF(R.NE.0.0D0) THEN
            V = FUNL(2.0D0*R/CMPF,1)/(PI*CMPF)
          ENDIF
C
C         EVALUATE V AS IF IN NUCLEAR REGION
          IF(ISWITCH.EQ.0) GOTO 222
C
C         INTIALISE POTENTIAL COUNTER
          VN = 0.0D0
C
C         INTEGRATE OVER CHARGE SOURCE
          DO M=0,NSRC
C
C           SET SOURCE RADIUS
            S = HS*DFLOAT(M)
C
C           FUNCTION ARGUMENTS
            XM = 2.0D0*DABS(R-S)/CMPF
            XP = 2.0D0*DABS(R+S)/CMPF
C
C           CHEBYSHEV APPROXIMATION FIT (FULLERTON+RINKER)
            CM = FUNL(XM,0)
            CP = FUNL(XP,0)
C
C           CONTRIBUTION TO INTEGRAND
            VN = VN + EXTINT11(S*RHO(M)*(CM-CP),M,NSRC)
C           VN = VN + EXTINT5(S*RHO(M)*(CM-CP),M,NSRC)
C
          ENDDO
C
C         INTEGRATION WEIGHTING FACTORS
          VN = 5.0D0*HS*VN/2.99376D+5
C         VN = 2.0D0*HS*VN/45.0D0
C
C         IF FULL AND ASYMPTOTIC RESULTS ARE SIMILAR, FLIP SWITCH
          IF(DABS((VN-V)/VN)*ZNUC(IZ).LT.1.0D-5) THEN
            ISWITCH = 0
          ENDIF
C
C         OVERRIDE POTENTIAL NOW
          V = VN
C
222       CONTINUE
C
C         OTHER FACTORS AND FINAL VALUE R*R*V(R)
          VVAC(IZ,N) =-CMPF*R*V/(PI*CV*CV)
C
        ENDDO
C      
      ENDIF
C
C**********************************************************************C
C     PLOTTING SECTION                                                 C
C**********************************************************************C
C
      IF(IWRT.NE.1) GOTO 20
C
C     NAME TAG FOR PDF DOCUMENTS
      IF(INT(ANUC(IZ)).LT.10) THEN
        WRITE(ATRM,'(A,I1)') '00',INT(ANUC(IZ))
      ELSEIF(INT(ANUC(IZ)).LT.100) THEN
        WRITE(ATRM,'(A,I2)') '0',INT(ANUC(IZ))
      ELSE
        WRITE(ATRM,'(I3)') INT(ANUC(IZ))
      ENDIF
C
C     DETAILS COMMON TO ALL PLOTS
      XAXIS  = 'r/RNUC'
      KEY(1) = 'Exact Källén-Sabry'
      KEY(2) = 'Best fit'
C
C     UNIFORM GRID VALUES R*R*V(R)
      XOUT   = 'Källén-Sabry-uniform'//TRIM(ATRM)
      TITLE  = 'Källén-Sabry r^{2}*V(r) on uniform grid'
      YAXIS  = 'r^{2}*V(r)/Z'
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
      DO IPTS=0,NPTS
        R = RU(IPTS)
        WRITE(8, *) R/RNUC(IZ),R*R*VF(IPTS),R*R*VG(IPTS)
      ENDDO
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C     NUMERICAL EVALUATION FOR V(R)
      XOUT   = 'Källén-Sabry-pointwise'//TRIM(ATRM)
      TITLE  = 'Källén-Sabry V(r) on piecewise grid'
      YAXIS  = 'V(r)/Z'
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
      DO N=0,NRAD
        R = RAD(N)
        IF(R/RNUC(IZ).GT.5.0D0) GOTO 123
C       TRUE KÄLLÉN-SABRY POTENTIAL
        IF(N.EQ.0) THEN
          VTR = V0
        ELSE
          VTR = VVAC(IZ,N)/(R*R)
        ENDIF
C       APPROXIMATE VALUES
        V = 0.0D0
        DO IFT=1,NFT
          V = V + FKSB(IZ,IFT)*DEXP(-XKSB(IZ,IFT)*R*R)
        ENDDO
        WRITE(8, *) R*cfm,VTR,V
      ENDDO
123   CONTINUE
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C     UNIFORM GRID VALUES R*R*V(R)
      XOUT   = 'Källén-Sabry-charge'//TRIM(ATRM)
      TITLE  = 'Källén-Sabry rho(r) on uniform grid'
      YAXIS  = 'rho(r)/Z'
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
      DO IPTS=0,NPTS/25
        R = RU(IPTS)
C       NUMERICAL INVERSE HANKEL TRANSFORM (FOR A CHARGE DENSITY)
        CALL PRIHNKT(GKSB,R,RHOKSB)
        WRITE(8, *) R/rnuc(iz),-4.0D0*R*R*PI*RHOKSB,POL(IPTS)
      ENDDO
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,2,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
20    CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE UEHPLOT(IZ,RZER,QPOL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     UU    UU EEEEEEEE HH    HH PPPPPPP  LL       OOOOOO TTTTTTTT     C
C     UU    UU EE       HH    HH PP    PP LL      OO    OO   TT        C
C     UU    UU EE       HH    HH PP    PP LL      OO    OO   TT        C
C     UU    UU EEEEEE   HHHHHHHH PP    PP LL      OO    OO   TT        C
C     UU    UU EE       HH    HH PPPPPPP  LL      OO    OO   TT        C
C     UU    UU EE       HH    HH PP       LL      OO    OO   TT        C
C      UUUUUU  EEEEEEEE HH    HH PP       LLLLLLLL OOOOOO    TT        C
C                                                                      C
C -------------------------------------------------------------------- C
C  UEHPLOT GENERATES AND PLOTS SOME POLARISED CHARGED DENSITIES,       C
C  COMPARING AGAINST THE ACTUAL (NORMALISED) DENSITY AND APPROX. FORMS.C
C**********************************************************************C
      INCLUDE 'parameters.h'
      PARAMETER(NPTS=500,MPTS=3000)
C
      CHARACTER*3  ATRM
      CHARACTER*5  NMDL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(4)
C
      DIMENSION RU(0:NPTS),PN(0:NPTS)
      DIMENSION DN(0:NPTS),DP(0:NPTS),DU(0:NPTS),DB(0:NPTS)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BUEE/RUEE(MCT,3),FUEE(MCT,MFT),XUEE(MCT,MFT),NUEE(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     MULTIPLIER FOR COMPTON WAVELENGTH (USE FOR MUON OR TAUON FIELD)
      CMPF = CMPW/1.0D0
C
C     NUMBER OF FITTING FUNCTIONS
      NFT = NUEE(IZ)
C
C     ZEROTH, SECOND AND FOURTH MOMENTS OF POLARISED CHARGE DENSITY
      RUEE(IZ,1) = 0.0D0
      RUEE(IZ,2) = 2.0D0*CMPF*CMPF/(5.0D0*PI*CV)
      RUEE(IZ,3) = 4.0D0*CMPF*CMPF*RNUC(IZ)*RNUC(IZ)/(3.0D0*PI*CV) 
     &           + 6.0D0*CMPF*CMPF*CMPF*CMPF/(7.0D0*PI*CV)
C
C     EXTENDED INTEGRATION FACTOR
      FI = 5.0D0/2.99376D5
C
C     GENERATE UNIFORM RADIAL GRID
      RLIM = 6.0D0
      DO IPTS=0,NPTS
        RU(IPTS) = RLIM*DFLOAT(IPTS)*RNUC(IZ)/DFLOAT(NPTS)
      ENDDO
      HR = RLIM*RNUC(IZ)/DFLOAT(NPTS)
C
C**********************************************************************C
C     NUCLEAR CHARGE DENSITY (NORMALISED AND RESCALED BY FACTOR CMPF)  C
C**********************************************************************C
C
C     NOW EVALUATE THE WEIGHTED DENSITY
      DO IPTS=0,NPTS
        R = RU(IPTS)
        DN(IPTS) = 4.0D0*PI*R*R*RHONUC(NMDL(IZ),IZ,R)
      ENDDO
C
C**********************************************************************C
C     NUMERICAL POLARISED CHARGE DENSITY                               C
C**********************************************************************C
C
C     GOTO 100
C     PRE-CALCULATE NUCLEAR CHARGE DENSITY ALONG SUITABLE GRID
      DO KPTS=0,NPTS
        S = RU(KPTS)
        PN(KPTS) = RHONUC(NMDL(IZ),IZ,S)
      ENDDO
C
C     VACUUM POLARISATION CHARGE DENSITY
      DO IPTS=0,NPTS
C
C       RADIUS R
        R = RU(IPTS)
C
C       INTIALISE CHARGE COUNTER HERE
        Q = 0.0D0
C
C       START INTEGRATION OVER ENERGY VARIABLE T
        T0 = 1.0D0
        TN = 1.0D2
        HT = (TN-T0)/DFLOAT(MPTS)
        DO JPTS=0,MPTS
C
          T = T0 + (TN-T0)*DFLOAT(JPTS)/DFLOAT(MPTS)
          U = 1.0D0/T
          E = U*DSQRT(1.0D0-U*U)*(1.0D0+0.5D0*U*U)
C
C         START INTEGRATION OVER SOURCE CHARGE VARIABLE S
          P = 0.0D0
          DO KPTS=0,NPTS
            S  = RU(KPTS)
            T1 = DEXP(-2.0D0*DABS(R-S)*T/CMPF)
            T2 = DEXP(-2.0D0*DABS(R+S)*T/CMPF)
            P  = P + EXTINT11(S*PN(KPTS)*(T1-T2),KPTS,NPTS)
          ENDDO
C
          G = FI*T*R*HR*P/CMPF - 0.25D0*DN(IPTS)/PI
          Q = Q + EXTINT11(E*G,JPTS,MPTS)
C          
        ENDDO
C
        DP(IPTS) = 8.0D0*FI*HT*Q/(3.0D0*CV)
C
      ENDDO
C
C     LOCATION OF CHARGE ZERO (LINE OF BEST FIT BETWEEN TWO POINTS)
      ZP = 0.0D0
      DO IPTS=NPTS,0,-1
        IF(DP(IPTS).LT.0.0D0) THEN
          PD = DP(IPTS)-DP(IPTS-1)
          ZP = (RU(IPTS-1)*DP(IPTS)-RU(IPTS)*DP(IPTS-1))/PD
          GOTO 51
        ENDIF
      ENDDO
51    CONTINUE
C
C     INTEGRATED POLARISED CHARGE UP TO THIS ZERO
      NPOL = 0
      DO IPTS=1,NPTS
        NPOL = NPOL+1
        IF(RU(IPTS).GT.ZP.AND.MOD(IPTS,10).EQ.0) GOTO 52
      ENDDO
52    CONTINUE

      QP = 0.0D0
      DO IPTS=1,NPOL
        QP = QP - EXTINT11(DP(IPTS),IPTS,NPOL)
      ENDDO
      QP = FI*HR*QP
100   CONTINUE
C
C**********************************************************************C
C     UEHLING'S APPROXIMATE CHARGE DENSITY                             C
C**********************************************************************C
C
C     CAN DO THIS EXACTLY IN CASE OF GAUSSIAN AND FERMI MODEL
      WRITE(6, *) cmpf,cv
      WRITE(7, *) cmpf,cv
      FR =-CMPF*CMPF*CMPF/(15.0D0*CV*PI*PI)
      IF(NMDL(IZ).EQ.'GAUSS') THEN
        DO IPTS=0,NPTS
          R = RU(IPTS)
          B = 9.0D0/(RNUC(IZ)*RNUC(IZ)*RNUC(IZ)*RNUC(IZ))
          DU(IPTS) = FR*B*(RNUC(IZ)-R)*(RNUC(IZ)+R)*DN(IPTS)
        ENDDO
      ELSEIF(NMDL(IZ).EQ.'FERMI') THEN
        DO IPTS=0,NPTS
          R = RU(IPTS)
C
C         FITTING PARAMETERS
          A = AFMI(IZ)
          C = 5.0D0*RNUC(IZ)*RNUC(IZ)/3.0D0
     &      - 7.0D0*PI*PI*AFMI(IZ)*AFMI(IZ)/3.0D0
          C = DSQRT(C)
          U  = A/C
          C3 = C*C*C
C
C         NORMALISATION CONSTANT
          Y1 = 1.0D0
          Y2 = PI*PI*U*U
          Y3 =-6.0D0*U*U*U*POLYLOG(3,-1.0D0/U)
          YN = 1.0D0/(Y1+Y2+Y3)
          YN = 3.0D0*YN/(4.0D0*PI*C3)
C
C         BITS AND PIECES
          X1 = PI*R*YN/(A*A)
          X2 = DCOSH(0.5D0*(R-C)/A)
          X3 = 2.0D0*A - R*DTANH(0.5D0*(R-C)/A)
C
          DU(IPTS) = FR*X1*X3/(X2*X2)
C
        ENDDO
      ELSE
        DO IPTS=0,NPTS
          DU(IPTS) = 0.0D0
        ENDDO
      ENDIF
      DO IPTS=0,NPTS
        DU(IPTS) = 0.0D0
      ENDDO
C
C     LOCATION OF CHARGE ZERO (LINE OF BEST FIT BETWEEN TWO POINTS)
      ZU = 0.0D0
      DO IPTS=NPTS,0,-1
        IF(DU(IPTS).LT.0.0D0) THEN
          PD = DU(IPTS)-DU(IPTS-1)
          ZU = (RU(IPTS-1)*DU(IPTS)-RU(IPTS)*DU(IPTS-1))/PD
          GOTO 53
        ENDIF
      ENDDO
53    CONTINUE
C
C     INTEGRATED POLARISED CHARGE UP TO THIS ZERO
      NPOL = 0
      DO IPTS=1,NPTS
        NPOL = NPOL+1
        IF(RU(IPTS).GT.ZU.AND.MOD(IPTS,10).EQ.0) GOTO 54
      ENDDO
54    CONTINUE

      QU = 0.0D0
      DO IPTS=1,NPOL
        QU = QU - EXTINT11(DU(IPTS),IPTS,NPOL)
      ENDDO
      QU = FI*HR*QU
C
C**********************************************************************C
C     BEST-FIT CHARGE DENSITY                                          C
C**********************************************************************C
C
C     VACUUM POLARISATION CHARGE DENSITY
      DO IPTS=0,NPTS
C
C       RADIUS R
        R = RU(IPTS)
C
C       INTIALISE CHARGE COUNTER
        Q = 0.0D0
        DO IFT=1,NFT
          FR = FUEE(IZ,IFT)
          XI = XUEE(IZ,IFT)
          Q  = Q + 2.0D0*FR*XI*R*R*(3.0D0-2.0D0*R*R*XI)*DEXP(-XI*R*R)
        ENDDO
        DB(IPTS) = Q
C
      ENDDO
C
C     LOCATION OF CHARGE ZERO (LINE OF BEST FIT BETWEEN TWO POINTS)
      ZB = RZER
C
C     INTEGRATED POLARISED CHARGE UP TO THIS ZERO
      QB = QPOL
C
C**********************************************************************C
C     PLOTTING SECTION                                                 C
C**********************************************************************C
C
      WRITE(6, *) IZ
      WRITE(7, *) IZ
      WRITE(6, *) ZP/RNUC(IZ),QP
      WRITE(7, *) ZP/RNUC(IZ),QP
      WRITE(6, *) ZU/RNUC(IZ),QU
      WRITE(7, *) ZU/RNUC(IZ),QU
      WRITE(6, *) ZB/RNUC(IZ),QB
      WRITE(7, *) ZB/RNUC(IZ),QB
C
C     NAME TAG FOR PDF DOCUMENTS
      IF(INT(ANUC(IZ)).LT.10) THEN
        WRITE(ATRM,'(A,I1)') '00',INT(ANUC(IZ))
      ELSEIF(INT(ANUC(IZ)).LT.100) THEN
        WRITE(ATRM,'(A,I2)') '0',INT(ANUC(IZ))
      ELSE
        WRITE(ATRM,'(I3)') INT(ANUC(IZ))
      ENDIF
C
C     DETAILS COMMON TO ALL PLOTS
      XAXIS  = 'r/RNUC'
      KEY(1) = 'Underlying nuclear density'
      KEY(2) = 'Numerical Uehling density'
      KEY(3) = 'Approx. Uehling density'
      KEY(4) = 'Best fit Uehling density'
C
C     UNIFORM GRID VALUES R*R*P(R)
      XOUT   = 'Uehling-polarised'//TRIM(ATRM)
      TITLE  = 'Uehling polarised radial charge density'
      YAXIS  = 'rho(r)/Z'
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
      DO IPTS=0,NPTS
        R = RU(IPTS)
        WRITE(8, *) R/RNUC(IZ),CMPF*DN(IPTS),DP(IPTS),DU(IPTS),DB(IPTS)
      ENDDO
      CLOSE(UNIT=8)
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,4,KEY)
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
      RETURN
      END
C
C
      FUNCTION FUNL(X,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 FFFFFFFF UU    UU NN    NN LL                        C
C                 FF       UU    UU NNN   NN LL                        C
C                 FF       UU    UU NNNN  NN LL                        C
C                 FFFFFF   UU    UU NN NN NN LL                        C
C                 FF       UU    UU NN  NNNN LL                        C
C                 FF       UU    UU NN   NNN LL                        C
C                 FF        UUUUUU  NN    NN LLLLLLLL                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  FUNL IS AN AUXILLIARY FUNCTION EVALUATOR FOR THE KÄLLÉN-SABRY       C
C  POTENTIAL BASED ON THE CHEBYSHEV POLYNOMIAL APPROXIMATION CARRIED   C
C  OUT IN FULLTERTON & RINKER (1975).                                  C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION F(0:5,0:1),G(0:2,0:1),H(0:1,0:1)
C
      DATA F/ 1.990159D+00,-2.397605D+00, 1.046471D+00,
     &       -0.367066D+00, 0.063740D+00,-0.037058D+00, 
     &        1.646407D+00,-2.092942D+00, 0.962310D+00,
     &       -0.254960D+00, 0.164404D+00, 0.000000D+00/
      DATA G/ 0.751198D+00, 0.138889D+00, 0.020886D+00, 
     &        0.137691D+00,-0.416667D+00,-0.097486D+00/
      DATA H/-0.444444D+00,-0.003472D+00,
     &        0.444444D+00, 0.017361D+00/
C
C     I FOUND THIS IN GRASP, NOT THE ORIGINAL PAPER
      DATA A,B/2.2D+00,-1.72D+00/
C     DATA A,B/6.720D+00,-6.647D+00/
C
C     INITIALISE THE FUNCTION VALUE
      FUNL = 0.0D0
C
C     THIS ROUTINE ONLY ACCEPTS K=0 AND K=1
      IF(K.NE.0.AND.K.NE.1) THEN
        WRITE(6, *) 'In FUNL: can only choose K=0 or K=1. K=',K
        WRITE(7, *) 'In FUNL: can only choose K=0 or K=1. K=',K
        RETURN
      ENDIF
C
C     THIS ROUTINE ONLY ACCEPTS X >= 0.0D0
      IF(X.LT.0.0D0) THEN
        WRITE(6, *) 'In FUNL: can only choose X.GE.0.0D0. X=',X
        WRITE(7, *) 'In FUNL: can only choose X.GE.0.0D0. X=',X
        RETURN
      ENDIF
C
C     SPECIAL CASE: X=0.0D0
      IF(X.EQ.0.0D0) THEN
        IF(K.EQ.0) FUNL = F(0,K)
        RETURN
      ENDIF
C
      IF(X.LE.2.0D0) THEN
C     CASE A: 0.0D+00 < X < 2.0D+00 -- RATIONAL APPROXIMATION
C
        XPOW = 1.0D0
        DO I=0,5
          FUNL = FUNL + F(I,K)*XPOW
          XPOW = XPOW*X
        ENDDO
C
        XPOW = X**(1-K)
        DO I=0,2
          FUNL = FUNL + G(I,K)*XPOW*DLOG(X)
          XPOW = XPOW*X*X
        ENDDO
C
        XPOW = X**(1-K)
        DO I=0,1
          FUNL = FUNL + H(I,K)*XPOW*DLOG(X)*DLOG(X)
          XPOW = XPOW*X*X*X*X
        ENDDO
C
      ELSE
C     CASE B: X > 2.0D+00 -- ASYMPTOTIC EXPANSION
C
        IF(K.EQ.0) THEN
          FUNL = A + B/DSQRT(X)
        ELSEIF(K.EQ.1) THEN
          FUNL = A*DSQRT(X)*(X+2.0D0)+0.5D0*B*(2.0D0*X+5.0D0)
        ENDIF
        XPOW = X**(4+3*K)
        FUNL = FUNL*DEXP(-X)/DSQRT(XPOW)
C
C        I FOUND THIS IN GRASP, NOT THE ORIGINAL PAPER
         FUNL = A + B/X
         IF(K.NE.0) THEN
           FUNL = FUNL + (FUNL + B/X)/X
         ENDIF
         FUNL = FUNL*DEXP(-X)/X
C
      ENDIF
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
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
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
      CHARACTER*5 NMDL,MODEL
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
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
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
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
      FUNCTION VGAUSS(R,XI,FRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        VV    VV  GGGGGG     AA    UU    UU  SSSSSS   SSSSSS          C
C        VV    VV GG    GG   AAAA   UU    UU SS    SS SS    SS         C
C        VV    VV GG        AA  AA  UU    UU SS       SS               C
C        VV    VV GG       AA    AA UU    UU  SSSSSS   SSSSSS          C
C         VV  VV  GG   GGG AAAAAAAA UU    UU       SS       SS         C
C          VVVV   GG    GG AA    AA UU    UU SS    SS SS    SS         C
C           VV     GGGGGG  AA    AA  UUUUUU   SSSSSS   SSSSSS          C
C                                                                      C
C -------------------------------------------------------------------- C
C  VGAUSS EVALUATES THE POTENTIAL ARISING FROM ONE TERM IN THE SET OF  C
C  GAUSSIAN DISTRIBUTIONS AT A PARTICULAR RADIUS R.                    C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
      IF(R.LT.1.0D-10) THEN
        XR2 = XI*R*R
        V0  = 2.0D0-2.0D0*XR2/3.0D0+XR2*XR2/5.0D0-XR2*XR2*XR2/21.0D0
        VGAUSS = -V0*DSQRT(XI)/PI12
      ELSE
        VGAUSS = -DERF(DSQRT(XI)*R)/R
      ENDIF
C
C     FINAL VALUE FOR THE POTENTIAL
      VGAUSS = FRAC*VGAUSS
C
      RETURN
      END
C
C
      FUNCTION VFERMI(R,C,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         VV    VV FFFFFFFF EEEEEEEE RRRRRRR  MM       MM IIII         C
C         VV    VV FF       EE       RR    RR MMM     MMM  II          C
C         VV    VV FF       EE       RR    RR MMMM   MMMM  II          C
C         VV    VV FFFFFF   EEEEEE   RR    RR MM MM MM MM  II          C
C          VV  VV  FF       EE       RRRRRRR  MM  MMM  MM  II          C
C           VVVV   FF       EE       RR    RR MM   M   MM  II          C
C            VV    FF       EEEEEEEE RR    RR MM       MM IIII         C
C                                                                      C
C -------------------------------------------------------------------- C
C  VFERMI EVALUATES THE NORMALISED SCALAR POTENTIAL ARISING FROM THE   C
C  TWO-PARAMETER FERMI NUCLEAR CHARGE DISTRIBUTION AT RADIUS R, GIVEN  C
C  HALF-DENSITY RADIUS C AND THE NUCLEAR SKIN THICKNESS T.             C
C  THE POTENTIAL IS EVALUATED IN CASES DEPENDING ON R AND C, AND       C
C  INVOLVES THE USE OF AN AUXILLIARY FUNCTION POLYLOG AS WELL.         C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
C     RADIUS MUST BE NON-NEGATIVE
      IF(R.LT.0.0D0) THEN
        WRITE(6, *) 'In VFERMI: radius R must be non-negative. R = ',R
        WRITE(7, *) 'In VFERMI: radius R must be non-negative. R = ',R
        STOP
      ENDIF
C
C     WARN THE USER IF C<A BUT DON'T STOP
      IF(C.LT.0.25D0*T/THLG) THEN
        WRITE(6, *) 'In VFERMI: A/C > 1 so expressions diverge.'
        WRITE(7, *) 'In VFERMI: A/C > 1 so expressions diverge.'
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
      VFERMI =-YN*VR
C
      RETURN
      END
C
C
      FUNCTION VUNIFM(R,RNUC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         VV    VV UU    UU NN    NN IIII FFFFFFFF MM       MM         C
C         VV    VV UU    UU NNN   NN  II  FF       MMM     MMM         C
C         VV    VV UU    UU NNNN  NN  II  FF       MMMM   MMMM         C
C         VV    VV UU    UU NN NN NN  II  FFFFFF   MM MM MM MM         C
C          VV  VV  UU    UU NN  NNNN  II  FF       MM  MMM  MM         C
C           VVVV   UU    UU NN   NNN  II  FF       MM   M   MM         C
C            VV     UUUUUU  NN    NN IIII FF       MM       MM         C
C                                                                      C
C -------------------------------------------------------------------- C
C  VUNIFM IS THE NORMALISED COULOMB POTENTIAL FROM A UNIFORM CHARGE    C
C  DENSITY WITH ROOT MEAN SQUARE VALUE RNUC.                           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
C     RADIUS MUST BE NON-NEGATIVE
      IF(R.LT.0.0D0) THEN
        WRITE(6, *) 'In VUNIFM: radius R must be non-negative. R = ',R
        WRITE(7, *) 'In VUNIFM: radius R must be non-negative. R = ',R
        STOP
      ENDIF
C
C     UNIFORM POTENTIAL
      CU = DSQRT(5.0D0/3.0D0)*RNUC
      IF(R.LT.CU) THEN
        VUNIFM =-0.5D0*(3.0D0 - (R/CU)**2)/CU
      ELSE
        VUNIFM =-1.0D0/R
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION VFITGS(R,IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           VV    VV FFFFFFFF IIII TTTTTTTT GGGGGG   SSSSSS            C
C           VV    VV FF        II     TT   GG    GG SS    SS           C
C           VV    VV FF        II     TT   GG       SS                 C
C           VV    VV FFFFFF    II     TT   GG        SSSSSS            C
C            VV  VV  FF        II     TT   GG   GGG       SS           C
C             VVVV   FF        II     TT   GG    GG SS    SS           C
C              VV    FF       IIII    TT    GGGGGG   SSSSSS            C
C                                                                      C
C -------------------------------------------------------------------- C
C  VFITGS EVALUATES THE COULOMB POTENTIAL ARISING FROM A FITTED SET    C
C  OF NUCLEAR FUNCTIONS.                                               C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5  NMDL
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
C     RADIUS MUST BE NON-NEGATIVE
      IF(R.LT.0.0D0) THEN
        WRITE(6, *) 'In VFITGS: radius R must be non-negative. R = ',R
        WRITE(7, *) 'In VFITGS: radius R must be non-negative. R = ',R
        STOP
      ENDIF
C
C     START WITH A GAUSSIAN POTENTIAL
      VFITGS = VGAUSS(R,XNUC(IZ,0),1.0D0)
C
C     ADD A SERIES OF FITTED FUNCTIONS
      DO IFT=1,NNUC(IZ)
        VFITGS = VFITGS + FNUC(IZ,IFT)*DEXP(-XNUC(IZ,IFT)*R*R)
      ENDDO
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
C   TTTTTTTT IIII MM       MM EEEEEEEE NN    NN  OOOOOO  WW      WW    C
C      TT     II  MMM     MMM EE       NNN   NN OO    OO WW      WW    C
C      TT     II  MMMM   MMMM EE       NNNN  NN OO    OO WW      WW    C
C      TT     II  MM MM MM MM EEEEEE   NN NN NN OO    OO WW  WW  WW    C
C      TT     II  MM  MMM  MM EE       NN  NNNN OO    OO WW WWWW WW    C
C      TT     II  MM   M   MM EE       NN   NNN OO    OO WWWW  WWWW    C
C      TT    IIII MM       MM EEEEEEEE NN    NN  OOOOOO  WW      WW    C
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
C**********************************************************************C
C ==================================================================== C
C   [4] ATOMIC HARTREE-FOCK: SINGLE-CENTRE SCF CALCULATIONS.           C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] ERFINT0: INTEGRAL OVER A GAUSSIAN AND ERROR FUNCTION.          C
C   [N] GAMLWR: LOWER INCOMPLETE GAMMA FUNCTION gamma(A,X).            C
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
C**********************************************************************C
C ==================================================================== C
C  [10] PLOTS: AMPLITUDES AND FIELDS/POTENTIALS IN DATA FILES.         C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] FIELDS: CALCULATE AMPLITUDES, FIELDS AND POTENTIALS.           C
C   [B] DGNUMAP: EXPORT AND PLOT A SQUARE ARRAY OF DOUBLES.            C
C   [C] ZGNUMAP: EXPORT AND PLOT A SQUARE ARRAY OF COMPLEX DOUBLES.    C
C   [D] ORBCOEF: PRINT LIST OF EXPANSION COEFFICIENTS FOR AN ORBITAL.  C
C   [E] AMPLTDE: PLOT A DIRAC SPINOR AMPLITUDE ALONG ONE DIRECTION.    C
C   [F] DENSMAP: HEAT MAP OF ELECTRONIC CHARGE DENSITY (ON A 2D GRID). C
C   [G] J4CRRNT: PLOT 4-CURRENT ALONG ONE DIRECTION.                   C
C   [H] POTENTL: PLOT 4-POTENTIAL ALONG ONE DIRECTION.                 C
C   [I] ELCTRCF: PLOT E FIELD ALONG ONE DIRECTION.                     C
C   [J] MAGNTCF: PLOT B FIELD ALONG ONE DIRECTION.                     C
C   [K] FRMFCTR: PLOT THE ELECTRON SCATTERING FORM FACTOR.             C
C   [L] GNULINE: GENERATE A GNUPLOT LINE PLOT MAKE FILE FOR DATA SET.  C
C   [M] GNUDENS: GENERATE A GNUPLOT MATRIX PLOT MAKE FILE FOR DENSITY. C
C   [N] CLEBSCH: CLEBSCH-GORDAN COEFFICIENT FOR KQN,MQN.               C
C**********************************************************************C
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
      IF(NDAT.LT.6) THEN
        WRITE(9,'(A)') 'load "plots/plotstyles.pal"'
      ELSE
        WRITE(9,'(A)') 'load "plots/plotstyles2.pal"'
      ENDIF
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
     &                        //' using 1:',N+1,' with lines ls ',N,
     &                                 ' title "'//TRIM(KEY(N))//'"'
        ELSEIF(NDAT.GT.1.AND.N.EQ.1) THEN
          WRITE(9,'(A,I2,A,I2,A)') 'plot "plots/'//TRIM(XOUT)//'.dat"'
     &                        //' using 1:',N+1,' with lines ls ',N,
     &                                 ' title "'//TRIM(KEY(N))//'",\'
        ELSEIF(NDAT.GT.1.AND.N.GT.1.AND.N.LT.NDAT) THEN
          WRITE(9,'(A,I2,A,I2,A)') '     "plots/'//TRIM(XOUT)//'.dat"'
     &                        //' using 1:',N+1,' with lines ls ',N,
     &                                 ' title "'//TRIM(KEY(N))//'",\'
        ELSEIF(NDAT.GT.1.AND.N.EQ.NDAT) THEN
          WRITE(9,'(A,I2,A,I2,A)') '     "plots/'//TRIM(XOUT)//'.dat"'
     &                        //' using 1:',N+1,' with lines ls ',N,
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
      SUBROUTINE GNUDENS(XOUT,TITLE,NHRZ,NVTC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     GGGGGG  NN    NN UU    UU DDDDDDD  EEEEEEEE NN    NN  SSSSSS     C
C    GG    GG NNN   NN UU    UU DD    DD EE       NNN   NN SS    SS    C
C    GG       NNNN  NN UU    UU DD    DD EE       NNNN  NN SS          C
C    GG       NN NN NN UU    UU DD    DD EEEEEE   NN NN NN  SSSSSS     C
C    GG   GGG NN  NNNN UU    UU DD    DD EE       NN  NNNN       SS    C
C    GG    GG NN   NNN UU    UU DD    DD EE       NN   NNN SS    SS    C
C     GGGGGG  NN    NN  UUUUUU  DDDDDDD  EEEEEEEE NN    NN  SSSSSS     C
C                                                                      C
C -------------------------------------------------------------------- C
C  GNUDENS IS A CONTROLLING ROUTINE THAT GENERATES A GNUPLOT MAKE FILE C
C  FOR A DENSITY MAP.                                                  C
C**********************************************************************C
C
      CHARACTER*80 XOUT,TITLE
C
      RATIO = DFLOAT(NVTC)/DFLOAT(NHRZ)
      SZHRZ = 10.0D0
      SZVTC = 10.0D0*RATIO
C
      OPEN(UNIT=9,FILE='plots/'//TRIM(XOUT)//'.gnuplot',
     &                                                 STATUS='REPLACE')
      WRITE(9,'(A)') '#'//TRIM(XOUT)//'.gnuplot'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '#  Usage:'
      WRITE(9,'(A)') '#  gnuplot < '//TRIM(XOUT)//'.gnuplot'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Terminal output specs'
      WRITE(9,'(A,F5.2)') 'set size ratio ',RATIO
      WRITE(9,'(A,F5.2,A,F5.2,A)') 'set terminal pdf size ',SZHRZ,'cm,',
     &                                     SZVTC,'cm font "palatino,10"'
      WRITE(9,'(A)') 'set output "plots/'//TRIM(XOUT)//'.pdf"'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Load line style definitions'
      WRITE(9,'(A)') 'load "plots/pals/jet.pal"'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Axes and title'
      WRITE(9,'(A)') 'set title sprintf("'//TRIM(TITLE)//'")'
      WRITE(9,'(A)') 'set xlabel "NHRZ"'
      WRITE(9,'(A)') 'set ylabel "NVTC"'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Grid style'
      WRITE(9,'(A)') 'set style line 102 lc rgb"#808080" lt 0 lw 1'
      WRITE(9,'(A)') 'set grid back ls 102'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Plotting options'
      WRITE(9,'(A)') 'set view map'
      WRITE(9,'(A)') 'set palette model RGB'
      WRITE(9,'(A)') 'set palette defined'
      WRITE(9,'(A)') '#'
      WRITE(9,'(A)') '# Plot data to file'
      WRITE(9,'(A,I2,A,I2)') 'splot "plots/'//TRIM(XOUT)//'.dat"'
     &                                    //' matrix with image notitle'
      CLOSE(UNIT=9)
C
      WRITE(6, *) 'Created command file "'//TRIM(XOUT)//'.gnuplot".'
      WRITE(7, *) 'Created command file "'//TRIM(XOUT)//'.gnuplot".'
C
      RETURN
      END
C
C
C**********************************************************************C
C ==================================================================== C
C  [15] QED: UEHLING ROUTINE DEVELOPMENT.                              C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] RHOASYM: ASYMPTOTIC FORM FOR UEHLING DENSITY AT LARGE R.       C
C**********************************************************************C
C
C
      SUBROUTINE RHOASYM(IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  RRRRRRR  HH    HH  OOOOOO     AA     SSSSSS  YY    YY MM       MM   C
C  RR    RR HH    HH OO    OO   AAAA   SS    SS YY    YY MMM     MMM   C
C  RR    RR HH    HH OO    OO  AA  AA  SS       YY    YY MMMM   MMMM   C
C  RR    RR HHHHHHHH OO    OO AA    AA  SSSSSS   YY  YY  MM MM MM MM   C
C  RRRRRRR  HH    HH OO    OO AAAAAAAA       SS   YYYY   MM  MMM  MM   C
C  RR    RR HH    HH OO    OO AA    AA SS    SS    YY    MM   M   MM   C
C  RR    RR HH    HH  OOOOOO  AA    AA  SSSSSS     YY    MM       MM   C
C                                                                      C
C -------------------------------------------------------------------- C
C  RHOASYM CALCULATES THE ASYMPTOTIC FORM FOR THE UEHLING CHARGE       C
C  DENSITY AT A LARGE RADIUS R, GIVEN A SINGLE GAUSSIAN CHARGE MODEL.  C
C  REQUIRES A NUMERICAL INTEGRATION OVER ENERGY DOMAIN T.              C
C  THIS CALCULATES R*RHO(R) BUT THEN PLOTS R*R*RHO(R).                 C
C  ------------------------------------------------------------------- C
C  TODO: FOR NOW, IGNORE FERMI MODEL TREATMENT AND GUESS THE ZETA      C
C        THAT MATCHES BEST WITH NUCLEAR RMS VALUE.                     C
C**********************************************************************C
      INCLUDE 'parameters.h'
      PARAMETER(NGRD=1000,NQDT=500)
C
      CHARACTER*5  NMDL
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(1)
C
      DIMENSION RAD(0:NGRD),RHO(0:NGRD)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     COMPTON WAVELENGTH MODIFIER
      CMPF = CMPW/1.0D0
C
C     WRITE NUCLEAR RADIUS (IN FM)
      RFM = CFM*RNUC(IZ)
      WRITE(6, *) 'Nuclear radius = ',RFM,'fm'
C
C     NUCLEAR WIDTH PARAMETER
      ZETA = 1.5D0/(RNUC(IZ)*RNUC(IZ))
C
C     LOOP OVER ALL RADIUS VALUES
      H = 200.0D0*RFM/(DFLOAT(NGRD)*CFM)
      DO N=0,NGRD
C
C       RADIUS AT THIS POINT
        RAD(N) = H*DFLOAT(N)
        RHO(N) = 0.0D0
C
C       LOOP OVER INTEGRAND POINTS T
        U = (5.0D0-1.0D0)/DFLOAT(NQDT)
        DO J=0,NQDT
C
          T  = 1.0D0 + U*DFLOAT(J)
C
          Y1 = 1.0D0/(T*T*T)
          Y2 = DSQRT(1.0D0 - 1.0D0/(T*T))
          Y3 = 1.0D0 + 0.5D0/(T*T)
          Y4 = DEXP(-2.0D0*T*RAD(N)*CV)
          Y5 = DEXP(EMSS*CV*CV/(ZETA*T*T))

          RHO(N) = RHO(N) + EXTINT11(Y1*Y2*Y3*Y4*Y5,J,NQDT)
          
        ENDDO
C
      ENDDO
C
C     MULTIPLICATIVE FACTOR
      FAC = 5.0D0*U/2.99376D5
      FAC = 8.0D0*FAC*CMPF/(3.0D0*PI*CV*CV)
      DO N=0,NGRD
        RHO(N) = FAC*RHO(N)
      ENDDO
C
C     FILE NAME AND TITLES
      XOUT   = 'Asymptotic'
      TITLE  = 'r*r*rho(r)'
      XAXIS  = 'r/RNUC'
      YAXIS  = 'r*r*rho(r)'
      KEY(1) = 'Gaussian'
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C       BEGIN LOOP OVER ALL POINTS
        DO N=0,NGRD
C
C         PRINT DISTANCE ALONG PATH AND FOUR-CURRENT
          WRITE(8, *) RAD(N)*CFM/RFM,RAD(N)*RHO(N)*CFM*1.0D3*4.0D0*PI
C
        ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,1,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
C
      RETURN
      END
C
C
      SUBROUTINE CHIGEN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            CCCCCC  HH    HH IIII GGGGGG  EEEEEEEE NN    NN           C
C           CC    CC HH    HH  II GG    GG EE       NNN   NN           C
C           CC       HH    HH  II GG       EE       NNNN  NN           C
C           CC       HHHHHHHH  II GG       EEEEEE   NN NN NN           C
C           CC       HH    HH  II GG   GGG EE       NN  NNNN           C
C           CC    CC HH    HH  II GG    GG EE       NN   NNN           C
C            CCCCCC  HH    HH IIII GGGGGG  EEEEEEEE NN    NN           C
C                                                                      C
C -------------------------------------------------------------------- C
C  CHIGEN EXPLICITLY GENERATES THE FUNCTION CHI_2(X) AND SAVES IT TO   C
C  FILE, FOR LATER SPLINE INTERPOLATION IN UEHLING CALCULATIONS.       C
C  IT USES THE EXPONENTIAL INTEGRAL EXPANSION METHOD.                  C
C  THE SPLINE FUNCTION IS ACTUALLY χ_2(X)*EXP(X), LEAVING A POLYNOMIAL.C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION XS(0:NLW),Y2S(0:NLW),D2S(0:NLW),Y1S(0:NLW),D1S(0:NLW),
     &          XB(0:NUP),Y2B(0:NUP),D2B(0:NUP)
C
C     JOINING POINT FOR BOTH SPLINE FUNCTIONS
      XSPL = 0.60D0
C
C**********************************************************************C
C     INTERPOLATING FUNCTIONS IN REGION 0.0D0 < X <~ 0.6D0             C
C**********************************************************************C
C
C     REGION GRID LENGTHS AND ENDPOINTS
      NCOR = 800
      XBEG = 0.00D0
      XCOR = 1.00D-4
      XEND = 1.25D0*XSPL
C
C     VERY SMALL X: UNIFORMLY-SPACED GRID
      DO N=0,NCOR
        XS(N)  = XBEG + (XCOR-XBEG)*DFLOAT(N)/DFLOAT(NCOR)
        Y1S(N) = CHI1SF(XS(N))
        Y2S(N) = CHIFNC(XS(N),2,0)*DEXP(XS(N))
      ENDDO
C
C     LARGER X: EXPONENTIALLY-SPACED GRID
      SPC = DLOG(XEND/XCOR)/DFLOAT(NLW-NCOR)
      DO N=0,NLW-NCOR
        XS(N+NCOR)  = XCOR*DEXP(N*SPC)
        Y1S(N+NCOR) = CHI1SF(XS(N+NCOR))
        Y2S(N+NCOR) = CHIFNC(XS(N+NCOR),2,0)*DEXP(XS(N+NCOR))
      ENDDO
C
C     DERIVATIVES AT SPLINE ENDS
      D10 = 1.0D30
      D1N = CHI1SD(XEND)
      D20 =-1.0D30
      D2N = (CHIDRV(XEND,2,0)+CHIFNC(XEND,2,0))*DEXP(XEND)
C
C     GENERATE SPLINES
      CALL SPLNGEN(XS,Y1S,D1S,D10,D1N,NLW)
      CALL SPLNGEN(XS,Y2S,D2S,D20,D2N,NLW)
C
C     WRITE SPLINE DATA TO FILE
      OPEN(UNIT=50,FILE='spline/chi1x.dat',STATUS='UNKNOWN')
      REWIND(UNIT=50)
        WRITE(50, *) XSPL
        DO N=0,NLW
          WRITE(50, *) XS(N),Y1S(N),D1S(N)
        ENDDO
      CLOSE(UNIT=50)
C
C     WRITE SPLINE DATA TO FILE
      OPEN(UNIT=52,FILE='spline/chi2_small.dat',STATUS='UNKNOWN')
      REWIND(UNIT=52)
        WRITE(52, *) XSPL
        DO N=0,NLW
          WRITE(52, *) XS(N),Y2S(N),D2S(N)
        ENDDO
      CLOSE(UNIT=52)
C
C**********************************************************************C
C     INTERPOLATING FUNCTIONS IN REGION X >~ 0.6D0                     C
C**********************************************************************C
C
C     LARGE X ENDPOINTS
      XBEG =   0.90D0*XSPL
      XEND = 250.00D0
C
C     USE EXPONENTIALLY-SPACED GRID
      SPC = DLOG(XEND/XBEG)/DFLOAT(NUP)
      DO N=0,NUP
        XB(N)  = XBEG*DEXP(N*SPC)
        Y2B(N) = CHIFNC(XB(N),2,1)*DEXP(XB(N))
      ENDDO
C
C     DERIVATIVES AT SPLINE ENDS
      D20 = (CHIDRV(XBEG,2,1)+CHIFNC(XBEG,2,1))*DEXP(XBEG)
      D2N = (CHIDRV(XEND,2,1)+CHIFNC(XEND,2,1))*DEXP(XEND)
C
C     GENERATE SPLINE
      CALL SPLNGEN(XB,Y2B,D2B,D20,D2N,NUP)
C
C     WRITE SPLINE DATA TO FILE
      OPEN(UNIT=53,FILE='spline/chi2_big.dat',STATUS='UNKNOWN')
      REWIND(UNIT=53)
        WRITE(53, *) XSPL
        DO N=0,NUP
          WRITE(53, *) XB(N),Y2B(N),D2B(N)
        ENDDO
      CLOSE(UNIT=53)
C
      RETURN
      END
C
C
      FUNCTION CHIFNCQUAD(X,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  CCCCCC  HH    HH IIII FFFFFFFF NN    NN  CCCCCC   QQQQQQ   DDDDDDD  C
C CC    CC HH    HH  II  FF       NNN   NN CC    CC QQ    QQ  DD    DD C          
C CC       HH    HH  II  FF       NNNN  NN CC      QQ      QQ DD    DD C
C CC       HHHHHHHH  II  FFFFFF   NN NN NN CC      QQ      QQ DD    DD C
C CC       HH    HH  II  FF       NN  NNNN CC      QQ      QQ DD    DD C
C CC    CC HH    HH  II  FF       NN   NNN CC    CC QQ    QQ  DD    DD C
C  CCCCCC  HH    HH IIII FF       NN    NN  CCCCCC   QQQQQQ Q DDDDDDD  C
C                                                                      C
C -------------------------------------------------------------------- C      
C  EVALUATE THE FUNCTION CHI_N (X) FOR PARTICULAR VALUE OF X,          C
C  USING A TRUNCATED SERIES WITH NTR TERMS.                            C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
C     SENSITIVITY TOLERANCE PARAMETER
      EPS = 1.0D-8
C
C     CHECK THAT N IS WITHIN ALLOWED VALUES
      IF(N.LE.1.AND.DABS(X).LT.EPS) THEN
        WRITE(6, *) 'In CHIFNCQUAD: divergent behaviour when N<1. N=',N
        WRITE(7, *) 'In CHIFNCQUAD: divergent behaviour when N<1. N=',N
        CHIFNCQUAD = 1.0D30
      ENDIF
C
C     INTEGRATION PARAMETERS (UNIFORM GRID)
      TBEG = 1.0D0
      TEND = 1.0D4
      HT   = (TEND/TBEG)/DFLOAT(KTR)
C
      CHIFNCQUAD = 0.0D0
      DO K=0,KTR
        T = TBEG + HT*DFLOAT(K)
        TR = 1.0D0/T
        Y1 = TR**N
        Y2 = DSQRT(1.0D0 - TR*TR)
        Y3 = 1.0D0 + 0.5D0*TR*TR
        Y4 = DEXP(-X*T)
        Z  = Y1*Y2*Y3*Y4
        CHIFNCQUAD = CHIFNCQUAD + 5.0D0*HT*EXTINT11(Z,K,KTR)/2.99376D+5
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION CHIFNC(X,N,MODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           CCCCCC  HH    HH IIII FFFFFFFF NN    NN  CCCCCC            C
C          CC    CC HH    HH  II  FF       NNN   NN CC    CC           C          
C          CC       HH    HH  II  FF       NNNN  NN CC                 C
C          CC       HHHHHHHH  II  FFFFFF   NN NN NN CC                 C
C          CC       HH    HH  II  FF       NN  NNNN CC                 C
C          CC    CC HH    HH  II  FF       NN   NNN CC    CC           C
C           CCCCCC  HH    HH IIII FF       NN    NN  CCCCCC            C
C                                                                      C
C -------------------------------------------------------------------- C      
C  EVALUATE THE FUNCTION CHI_N (X) FOR PARTICULAR VALUE OF X,          C
C  USING A TRUNCATED SERIES WITH NTR TERMS.                            C
C                   ∞      _______                                     C
C         χ_N(X) = ∫ t^-N √1-1/t^2 (1 + 1/2t^2) exp(-x t) dt.          C
C                   1                                                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  ▶  MODE = 0: USE EXACT AND APPROX. CHI_N (0.0D0) TO MIN. ERROR      C
C  ▶  MODE = 1: DO NOT USE EXACT VALUE (LARGE X).                      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
C     SENSITIVITY TOLERANCE PARAMETER
      EPS = 1.0D-8
C
C     CHECK THAT N IS WITHIN ALLOWED VALUES
      IF(N.LE.1.AND.DABS(X).LT.EPS) THEN
        WRITE(6, *) 'In CHIFNC: divergent behaviour when N < 1. N = ',N
        WRITE(7, *) 'In CHIFNC: divergent behaviour when N < 1. N = ',N
        STOP
      ENDIF
C
C     SPECIAL CASE: X = 0.0D0
      IF(DABS(X).LE.EPS) THEN
        CHIFNC = 0.75D0*PI12*GAMHLF(N+3)/(DFLOAT(N-1)*GAMHLF(N+4))
        RETURN
      ENDIF
C
C     χ_N(X) WITH SUM OF EXPONENTIAL FUNCTIONS
C
C     INITIALISE THE CHIFNC COUNTER
      CHIFNC = 0.0D0
C
C     LOOP FROM K=0 TO TRUNCATED INFINITY
      PHS =-1.0D0
      DO K=0,KTR
C
        RK = DFLOAT(K)
        IF(K.EQ.0) THEN
          BNM = 1.0D0
        ELSE
          BNM = (1.5D0-RK)*BNM/RK
        ENDIF
C
        PHS =-PHS
        ENK = EXPINTE(2*K+N,X) + 0.5D0*EXPINTE(2*K+N+2,X)
C
        CHIFNC = CHIFNC + PHS*BNM*ENK
C
      ENDDO
C
C     χ_N(X) NEAR ZERO, TAKING A CLEVER DIFFERENCE FROM χ_N(0)
      IF(MODE.EQ.0.AND.N.GT.1) THEN
C
C       EXACT VALUE AT X = 0.0D0
        TRU = 0.75D0*PI12*GAMHLF(N+3)/(DFLOAT(N-1)*GAMHLF(N+4))
C
C       INITIALISE THE COUNTER FOR RESIDUAL
        RSD = 0.0D0
C
        PHS =-1.0D0
        DO K=0,KTR
C
          RK = DFLOAT(K)
          IF(K.EQ.0) THEN
            BNM = 1.0D0
          ELSE
            BNM = (1.5D0-RK)*BNM/RK
          ENDIF
C
          PHS =-PHS
          RAT = DFLOAT(6*K+3*N+1)/DFLOAT(2*(2*K+N-1)*(2*K+N+1))
C
          RSD = RSD + PHS*BNM*RAT
C
        ENDDO
C
        CHIFNC = CHIFNC+TRU-RSD
C
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION CHIDRV(X,N,MODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           CCCCCC  HH    HH IIII DDDDDDD  RRRRRRR  VV    VV           C
C          CC    CC HH    HH  II  DD    DD RR    RR VV    VV           C
C          CC       HH    HH  II  DD    DD RR    RR VV    VV           C
C          CC       HHHHHHHH  II  DD    DD RR    RR VV    VV           C
C          CC       HH    HH  II  DD    DD RRRRRRR   VV  VV            C
C          CC    CC HH    HH  II  DD    DD RR    RR   VVVV             C
C           CCCCCC  HH    HH IIII DDDDDDD  RR    RR    VV              C
C                                                                      C
C -------------------------------------------------------------------- C
C  CHI2DERIV EVALUATES THE INSTANTANEOUS DERIVATIVE OF CHI_N(X) USING  C
C  A SERIES OF EXPONENTIAL INTEGRALS, OF LENGTH NTR.                   C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
C     SENSITIVITY TOLERANCE PARAMETER
      EPS = 1.0D-8
C
C     CHECK THAT N IS WITHIN ALLOWED VALUES
      IF(N.LE.2.AND.DABS(X).LT.EPS) THEN
        WRITE(6, *) 'In CHIDRV: divergent behaviour when N < 2. N = ',N
        WRITE(7, *) 'In CHIDRV: divergent behaviour when N < 2. N = ',N
        STOP
      ENDIF
C
C     SPECIAL CASE: X = 0.0D0
      IF(DABS(X).LT.EPS) THEN
        CHIDRV =-0.75D0*PI12*GAMHLF(N+2)/(DFLOAT(N-2)*GAMHLF(N+3))
        RETURN
      ENDIF
C
C     χ'_N(X) WITH SUM OF EXPONENTIAL FUNCTIONS
C
C     INITIALISE THE CHIDRV COUNTER
      CHIDRV = 0.0D0
C
C     LOOP FROM K=0 TO TRUNCATED INFINITY
      PHS =-1.0D0
      DO K=0,KTR
C
        RK = DFLOAT(K)
        IF(K.EQ.0) THEN
          BNM = 1.0D0
        ELSE
          BNM = (1.5D0-RK)*BNM/RK
        ENDIF
C
        PHS =-PHS
        ENK = EXPINTE(2*K+N-1,X) + 0.5D0*EXPINTE(2*K+N+1,X)
C
        CHIDRV = CHIDRV - PHS*BNM*ENK
C
      ENDDO
C
C     χ'_N(X) NEAR ZERO, TAKING A CLEVER DIFFERENCE FROM χ'_N(0)
      IF(MODE.EQ.0.AND.N.GT.2) THEN
C
C       EXACT VALUE AT X = 0.0D0
        TRU =-0.75D0*PI12*GAMHLF(N+2)/(DFLOAT(N-2)*GAMHLF(N+3))
C
C       INITIALISE THE COUNTER FOR RESIDUAL
        RSD = 0.0D0
C
        PHS =-1.0D0
        DO K=0,KTR
C
          RK = DFLOAT(K)
          IF(K.EQ.0) THEN
            BNM = 1.0D0
          ELSE
            BNM = (1.5D0-RK)*BNM/RK
          ENDIF
C
          PHS =-PHS
          RAT = DFLOAT(6*K+3*N-2)/DFLOAT(2*(2*K+N-2)*(2*K+N))
C
          RSD = RSD - PHS*BNM*RAT
C
        ENDDO
C
        CHIDRV = CHIDRV+TRU-RSD
C
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION CHI1SF(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             CCCCCC  HH    HH IIII  11   SSSSSS  FFFFFFFF             C
C            CC    CC HH    HH  II  111  SS    SS FF                   C
C            CC       HH    HH  II   11  SS       FF                   C
C            CC       HHHHHHHH  II   11   SSSSSS  FFFFFF               C
C            CC       HH    HH  II   11        SS FF                   C
C            CC    CC HH    HH  II   11  SS    SS FF                   C
C             CCCCCC  HH    HH IIII 1111  SSSSSS  FF                   C
C                                                                      C
C -------------------------------------------------------------------- C      
C  EVALUATE THE FUNCTION X*CHI_1 (X) FOR PARTICULAR VALUE OF X, USING  C
C  A TRUNCATED SERIES WITH NTR TERMS.                                  C
C                    ∞      _______                                    C
C      X*χ_N(X) = X*∫ t^-N √1-1/t^2 (1 + 1/2t^2) exp(-x t) dt.         C
C                    1                                                 C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
C     X*χ_1(X) WITH SUM OF EXPONENTIAL FUNCTIONS
C
C     INITIALISE THE CHI1SF COUNTER
      CHI1SF = 0.0D0
C
C     LOOP FROM K=0 TO TRUNCATED INFINITY
      PHS =-1.0D0
      DO K=0,KTR
C
        RK = DFLOAT(K)
        IF(K.EQ.0) THEN
          BNM = 1.0D0
        ELSE
          BNM = (1.5D0-RK)*BNM/RK
        ENDIF
C
        PHS =-PHS
        E10 = 1.5D0*DEXP(-X)
        E11 =-DFLOAT(2*K+1)*EXPINTE(2*K+2,X)
        E12 =-DFLOAT(2*K+3)*EXPINTE(2*K+4,X)*0.5D0
        ENK = E10 + E11 + E12
C
        CHI1SF = CHI1SF + PHS*BNM*ENK
C
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION CHI1SD(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             CCCCCC  HH    HH IIII  11   SSSSSS  DDDDDDD              C
C            CC    CC HH    HH  II  111  SS    SS DD    DD             C
C            CC       HH    HH  II   11  SS       DD    DD             C
C            CC       HHHHHHHH  II   11   SSSSSS  DD    DD             C
C            CC       HH    HH  II   11        SS DD    DD             C
C            CC    CC HH    HH  II   11  SS    SS DD    DD             C
C             CCCCCC  HH    HH IIII 1111  SSSSSS  DDDDDDD              C
C                                                                      C
C -------------------------------------------------------------------- C      
C  EVALUATE THE DERIVATIVE OF THE FUNCTION X*CHI_1 (X) FOR PARTICULAR  C
C  VALUE OF X, USING A TRUNCATED SERIES WITH NTR TERMS.                C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
C     SENSITIVITY TOLERANCE PARAMETER
      EPS = 1.0D-8
C
C     CHECK THAT X IS WITHIN ALLOWED VALUES
      IF(DABS(X).LT.EPS) THEN
        WRITE(6, *) 'In CHI1SD: divergent behaviour when X = 0.0D0.'
        WRITE(7, *) 'In CHI1SD: divergent behaviour when X = 0.0D0.'
        STOP
      ENDIF
C
C     (X*χ_1(X))' WITH SUM OF EXPONENTIAL FUNCTIONS
C
C     INITIALISE THE CHI1SD COUNTER
      CHI1SD = 0.0D0
C
C     LOOP FROM K=0 TO TRUNCATED INFINITY
      PHS =-1.0D0
      DO K=0,KTR
C
        RK = DFLOAT(K)
        IF(K.EQ.0) THEN
          BNM = 1.0D0
        ELSE
          BNM = (1.5D0-RK)*BNM/RK
        ENDIF
C
        PHS =-PHS
        E10 = 1.5D0*DEXP(-X)
        E11 = DFLOAT(2*K+1)*EXPINTE(2*K+1,X)
        E12 = DFLOAT(2*K+3)*EXPINTE(2*K+3,X)*0.5D0
        ENK = E10 + E11 + E12
C
        CHI1SD = CHI1SD + PHS*BNM*ENK
C
      ENDDO
C
      RETURN
      END
C
C
      FUNCTION EXPINTE(N,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      EEEEEEEE XX     XX PPPPPPP IIII NN    NN TTTTTTTT EEEEEEEE      C
C      EE        XX   XX  PP    PP II  NNN   NN    TT    EE            C
C      EE         XX XX   PP    PP II  NNNN  NN    TT    EE            C
C      EEEEEE      XXX    PP    PP II  NN NN NN    TT    EEEEEE        C
C      EE         XX XX   PPPPPPP  II  NN  NNNN    TT    EE            C
C      EE        XX   XX  PP       II  NN   NNN    TT    EE            C
C      EEEEEEEE XX     XX PP      IIII NN    NN    TT    EEEEEEEE      C
C                                                                      C
C -------------------------------------------------------------------- C
C  EXPINTE EVALUATES THE GENERALISED EXPONENTIAL INTEGRAL, E_N(X).     C
C**********************************************************************C
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
C     SENSITIVITY TOLERANCE PARAMETER
      EPS = 1.0D-8
C
C     CHECK FOR ILLEGAL INPUT
      IF(X.LT.0.0D0) THEN
        WRITE(6, *) 'In EXPINTE: X < 0.0D0. X = ',X 
        WRITE(7, *) 'In EXPINTE: X < 0.0D0. X = ',X 
        RETURN
      ELSEIF(N.LT.2.AND.DABS(X).LT.EPS) THEN
        WRITE(6, *) 'In EXPINTE: cannot have X = 0.0D0 and N < 2.'
        WRITE(7, *) 'In EXPINTE: cannot have X = 0.0D0 and N < 2.'
        RETURN
      ENDIF
C
      IF(N.LE.0.AND.DABS(X).GT.EPS) THEN
C     CASE N < 0, X >= 0.0D0
C
        EXPINTE = DEXP(-X)/X
        DO K=1,N
          EXPINTE = DEXP(-X)/X + DFLOAT(K)*EXPINTE/X
        ENDDO
C
      ELSEIF(N.GE.2.AND.DABS(X).GT.EPS) THEN
C     CASE N < 2, X  = 0.0D0
C
        EXPINTE = 1.0D0/DFLOAT(N-1)
C
      ELSEIF(X.GT.0.0D0.AND.X.LE.1.0D0) THEN
C     CASE 0.0D0 <= X <= 1.0D0 (SMALL ARGUMENT, SERIES EVALUATION)
C
C       STARTING VALUE
        IF(N.EQ.1) THEN
          EXPINTE =-DLOG(X)-EULR
        ELSE
          EXPINTE = 1.0D0/DFLOAT(N-1)
        ENDIF
C
C       LOOP UNTIL CORRECTIONS LEVEL OFF
        SER = 1.0D0
C
        DO K=1,100
C
          SER =-SER*X/DFLOAT(K)
C
C         CORRECTION FACTOR
          IF(K.NE.N-1) THEN
            DEL =-SER/DFLOAT(K-N+1)
          ELSE
C         COMPUTE PSI(N)
            PSI =-EULR
            DO KI=1,N-1
              PSI = PSI + 1.0D0/DFLOAT(KI)
            ENDDO
            DEL = SER*(-DLOG(X) + PSI)
          ENDIF
C
          EXPINTE = EXPINTE + DEL
          IF(DABS(DEL).LT.DABS(EXPINTE)*1.0D-9) THEN
            RETURN
          ENDIF
C
        ENDDO
C
C       FAILURE TO EXIT
        WRITE(6, *) 'In EXPINTE: series method failed.'
        WRITE(7, *) 'In EXPINTE: series method failed.'
        RETURN
C
C
      ELSEIF(X.GT.1.0D0) THEN
C     CASE X >= 1.0D0 (LARGE ARGUMENT, LENZ'S ALGORITHM)
C
        B = X+N
        C = 1.0D+30
        D = 1.0D0/B
        H = D
C
        DO K=1,100
          A   =-K*DFLOAT(N+K-1)
          B   = B+2.0D0
          D   = 1.0D0/(A*D+B)
          C   = B + A/C
          DEL = C*D
          H   = H*DEL
C
          EXPINTE = H*DEXP(-X)
          IF(DABS(DEL-1.0D0).LT.1.0D-9) THEN
            RETURN
          ENDIF
C
        ENDDO
C
C       FAILURE TO EXIT
        WRITE(6, *) 'In EXPINTE: continued fraction method failed.'
        WRITE(7, *) 'In EXPINTE: continued fraction method failed.'
        RETURN
C
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE SPLNGEN(X,Y,D2,D10,D1N,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     SSSSSS  PPPPPPP  LL       NN    NN  GGGGGG  EEEEEEEE NN    NN    C
C    SS    SS PP    PP LL       NNN   NN GG    GG EE       NNN   NN    C
C    SS       PP    PP LL       NNNN  NN GG       EE       NNNN  NN    C
C     SSSSSS  PP    PP LL       NN NN NN GG       EEEEEE   NN NN NN    C
C          SS PPPPPPP  LL       NN  NNNN GG   GGG EE       NN  NNNN    C
C    SS    SS PP       LL       NN   NNN GG    GG EE       NN   NNN    C
C     SSSSSS  PP       LLLLLLLL NN    NN  GGGGGG  EEEEEEEE NN    NN    C
C                                                                      C
C -------------------------------------------------------------------- C
C  SPLNGEN IS A CUBIC SPLINE GENERATOR, TAKING AN ORDERED LIST OF      C
C  ARGUMENTS X(0:N) AND THE CORRESPONDING FUNCTION VALUES Y(0:N) AS    C
C  WELL AS ITS FIRST DERIVATIVE AT THE SPLINE BOUNDARIES, AND GIVING   C
C  AN ARRAY D2(0:N) OF ITS SECOND DERIVATIVES AT ALL ARGUMENT VALUES.  C
C  THERE ARE SPECIAL CASES WHEN THE BOUNDARY DERIVATIVES AREN'T KNOWN. C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C   ▶   X - ARRAY OF INPUT VALUES, X_I.                                C
C   ▶   Y - ARRAY OF FUNCTION VALUES, Y_I.                             C
C   ▶   N - LENGTH OF ARRAYS FOR X, Y AND D2.                          C
C   ▶ D10 - FIRST DERIVATIVE OF INTERPOLATING FUNCTION AT X_0.         C
C   ▶ D1N - FIRST DERIVATIVE OF INTERPOLATING FUNCTION AT X_N.         C
C OUTPUT:                                                              C
C   ▶  D2 - SECOND DERIVATIVES OF THE INTERPOLATING FUNCTION AT X_I.   C
C**********************************************************************C
C
      DIMENSION X(0:N),Y(0:N),D2(0:N),UT(0:N)
C
C     FORWARD EVALUATION STARTING WITH LOWER BOUNDARY CONDITION D10
C
      IF(DABS(D10).GT.1.0D+29) THEN
C       NATURAL OPTION (D1 LIKELY DIVERGES HERE)
        D2(0) = 0.0D0
        UT(0) = 0.0D0
      ELSE
C       SPECIFIED DERIVATIVE
        D2(0) =-0.5D0
        UT(0) = (3.0D0/(X(1)-X(0)))*((Y(1)-Y(0))/(X(1)-X(0))-D10)
      ENDIF
C
C     DECOMPOSITION LOOP OF THE TRIDIAGONAL ALGORITHM
      DO I=1,N-1
C
        XDP = X(I+1)-X(I  )
        XDM = X(I  )-X(I-1)
        XDB = X(I+1)-X(I-1)
        YDP = Y(I+1)-Y(I  )
        YDM = Y(I  )-Y(I-1)
C
        P     = XDM*D2(I-1)/XDB + 2.0D0
        D2(I) = (XDM/XDB-1.0D0)/P
        UT(I) = 6.0D0*(YDP/XDP-YDM/XDM) - XDM*UT(I-1)
        UT(I) = UT(I)/(P*XDB)
C
      ENDDO
C
C     BACKWARD CORRECTION STARTING WITH UPPER BOUNDARY CONDITION D1N
C
      IF(DABS(D1N).GT.1.0D+29) THEN
C
C       NATURAL OPTION (D1 LIKELY DIVERGES HERE)
        QN = 0.0D0
        UN = 0.0D0
C
      ELSE
C
C       SPECIFIED DERIVATIVE
        QN = 0.5D0
        UN = (3.0D0/(X(N)-X(N-1)))*(D1N-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
C
      ENDIF
C
C     BACK-SUBSTITUTION LOOP USING THE TRIDIAGONAL ALGORITHM
      D2(N) = (UN - QN*UT(N-1))/(QN*D2(N-1)+1.0D0)
      DO I=0,N-1
        K     = N-I
        D2(K) = D2(K)*D2(K+1) + UT(K)
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE SPLNINT(XI,YI,D2I,X,Y,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       SSSSSS  PPPPPPP  LL       NN    NN IIII NN    NN TTTTTTTT      C
C      SS    SS PP    PP LL       NNN   NN  II  NNN   NN    TT         C
C      SS       PP    PP LL       NNNN  NN  II  NNNN  NN    TT         C
C       SSSSSS  PP    PP LL       NN NN NN  II  NN NN NN    TT         C
C            SS PPPPPPP  LL       NN  NNNN  II  NN  NNNN    TT         C
C      SS    SS PP       LL       NN   NNN  II  NN   NNN    TT         C
C       SSSSSS  PP       LLLLLLLL NN    NN IIII NN    NN    TT         C
C                                                                      C
C -------------------------------------------------------------------- C
C  SPLNINT IS AN INTERPOLATOR ROUTINE THAT TAKES AN EXISTING SPLINE    C
C  (WITH ARGUMENTS, FUNCTION VALUES AND SECOND DERIVATIVES IN ARRAYS   C
C  XI, YI AND D2I OF LENGTH N) AND AN ARGUMENT X, THEN RETURNS THE     C
C  CUBIC SPLINE INTERPOLATED FUNCTION VALUE Y.                         C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C   ▶  XI - ARRAY OF SPLINE INPUT VALUES X_I.                          C
C   ▶  YI - ARRAY OF SPLINE FUNCTION VALUES Y_I = Y(X_I).              C
C   ▶ D2I - ARRAY OF SPLINE SECOND DERIVATIVE VALUES D2I = Y''(X_I).   C
C   ▶   N - LENGTH OF ARRAYS FOR XI, YI AND D2I.                       C
C   ▶   X - INPUT VALUE OF INTEREST.                                   C
C OUTPUT:                                                              C
C   ▶   Y - CUBIC-SPLINE INTERPOLATED FUNCTION VALUE OF ARGUMENT X.    C
C**********************************************************************C
C
      DIMENSION XI(0:N),YI(0:N),D2I(0:N)
C
C     SENSITIVITY TOLERANCE PARAMETER
      EPS = 1.0D-10
C
C     INITIALISE ENDPOINTS OF CORRECT SPLINE ADDRESS
      M1 = 0
      M2 = N
C
C     USE BISECTION TO SEARCH FOR SPLINE WINDOW.
C     IF SPLNINT ARGUMENTS ARE IN ORDER, STORE M1 AND M2 BETWEEN CALLS.
1     IF(M2-M1.GT.1) THEN
        K = (M2+M1)/2
        IF(XI(K).GT.X) THEN
          M2 = K
        ELSE
          M1 = K
        ENDIF
        GOTO 1
      ENDIF
C
C     SPACING BETWEEN THE NEIGHBOURING XI VALUES
      H = XI(M2)-XI(M1)
C
C     CHECK THAT SPLINE VALUES XI ARE DISTINCT
      IF(DABS(H).LE.EPS) THEN
        WRITE(6, *) 'In SPLNINT: bad XI input values.'
        WRITE(7, *) 'In SPLNINT: bad XI input values.'
        STOP
      ENDIF
C
C     EVALUATE CUBIC-SPLINE POLYNOMIAL
      R = (XI(M2)-X     )/H
      S = (X     -XI(M1))/H
      T = R*(R*R-1.0D0)*D2I(M1)
      U = S*(S*S-1.0D0)*D2I(M1)
C
      Y = R*YI(M1) + S*YI(M2) + H*H*(T+U)/6.0D0
C
      RETURN
      END
C
C
      SUBROUTINE VUEHGRD(IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    VV    VV UU    UU EEEEEEEE HH    HH  GGGGGG  RRRRRRR  DDDDDDD     C
C    VV    VV UU    UU EE       HH    HH GG    GG RR    RR DD    DD    C
C    VV    VV UU    UU EE       HH    HH GG       RR    RR DD    DD    C
C    VV    VV UU    UU EEEEEE   HHHHHHHH GG       RR    RR DD    DD    C
C     VV  VV  UU    UU EE       HH    HH GG   GGG RRRRRRR  DD    DD    C
C      VVVV   UU    UU EE       HH    HH GG    GG RR    RR DD    DD    C
C       VV     UUUUUU  EEEEEEEE HH    HH  GGGGGG  RR    RR DDDDDDD     C
C                                                                      C
C -------------------------------------------------------------------- C
C  VUEHGRD EVALUATES THE UEHLING POTENTIAL DUE TO A NUCLEAR CHARGE     C
C  SOURCE (GAUSSIAN SHAPE) FOR A SERIES OF RADII.                      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5  NMDL
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(1)
C
      DIMENSION XS(0:NLW),Y2S(0:NLW),D2S(0:NLW),Y1S(0:NLW),D1S(0:NLW),
     &          XB(0:NUP),Y2B(0:NUP),D2B(0:NUP)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/VPOL/RAD(0:NRAD),VUEH(0:NRAD),EUEH(0:NRAD),PUEH(0:NRAD)
C
C     MULTIPLIER FOR COMPTON WAVELENGTH (USE FOR MUON OR TAUON FIELD)
      CMPF = 1.0D0*CMPW
C
C     IMPORT SPLINE DATA FOR X*CHI_1(X) GIVEN 0 <= X <= XSPL
      OPEN(UNIT=50,FILE='spline/chi1x.dat',STATUS='UNKNOWN')
      REWIND(UNIT=50)
        READ(50, *) XSPL
        DO N=0,NLW
          READ(50, *) XS(N),Y1S(N),D1S(N)
        ENDDO
      CLOSE(UNIT=50)
C
C     IMPORT SPLINE DATA FOR CHI_2(X) GIVEN 0 <= X <= XSPL
      OPEN(UNIT=52,FILE='spline/chi2_small.dat',STATUS='UNKNOWN')
      REWIND(UNIT=52)
        READ(52, *) XSPL
        DO N=0,NLW
          READ(52, *) XS(N),Y2S(N),D2S(N)
        ENDDO
      CLOSE(UNIT=52)
C
C     IMPORT SPLINE DATA FOR CHI_2(X) GIVEN X >= XSPL
      OPEN(UNIT=53,FILE='spline/chi2_big.dat',STATUS='UNKNOWN')
      REWIND(UNIT=53)
        READ(53, *) XSPL
        DO N=0,NUP
          READ(53, *) XB(N),Y2B(N),D2B(N)
        ENDDO
      CLOSE(UNIT=53)
C
C     RADIAL PLOTTING INFORMATION
      R0 = 0.0D0/CFM
      RN = 2500.0D0/CFM
      HR = (RN-R0)/DFLOAT(NRAD)
C
C     SEARCH FOR CHARGE RADIUS SMAX FOR WHICH RHO(SMAX) < 1.0D-16
      S0 = 0.0D0
      SM = 0.0D0
31    SM = SM + 0.001D0/CFM
      PM = SM*rhonuc(NMDL(IZ),IZ,SM)
      IF(DABS(PM).GT.1.0D-10) GOTO 31
c      SM = DSQRT(5.0D0/3.0D0)*RNUC(IZ)
      HS = (SM-S0)/DFLOAT(NSRC)
C      
      RHOTOT = 0.0D0
      DO M=0,NSRC
        S = S0 + HS*DFLOAT(M)
        Z = S*S*rhonuc(NMDL(IZ),IZ,S)
        RHOTOT = RHOTOT + 5.0D0*HS*EXTINT11(Z,M,NSRC)/2.99376D+5
      ENDDO
C
C     GENERATE RADIAL GRID AND POTENTIALS
      DO N=0,NRAD
C
C       SET RADIUS
        R = R0 + HR*DFLOAT(N)
        RAD(N) = R
C
C       SPECIAL CASE: R = 0.0D0 (UEHLING POTENTIAL AT THE ORIGIN)
        IF(N.EQ.0) THEN
C
C         INITIALISE COUNTER FOR V(R)
          V = 0.0D0
C
C         INTEGRATE OVER CHARGE SOURCE (DON'T INCLUDE S = 0.0D0 CASE)
          DO M=1,NSRC
C
C           SET SOURCE RADIUS
            S  = S0 + HS*DFLOAT(M)
C
C           CHI FUNCTION ARGUMENTS
            X0 = 2.0D0*S/CMPF
C
C           COMPONENTS OF INTEGRAND
            CALL SPLNINT(XS,Y1S,D1S,X0,C0,NLW)
C
C           CONTRIBUTION TO INTEGRAND
            Z = rhonuc(NMDL(IZ),IZ,S)*C0
            V = V + 5.0D0*HS*EXTINT11(Z,M,NSRC)/2.99376D+5
C
          ENDDO
C
C         MULTIPLICATIVE FACTORS FOR VUEH(N)
          VUEH(N) =-4.0D0*CMPF*V/(3.0D0*CV)
C
        ELSE
C       ALL NON-ZERO RADII
C
C         INITIALISE COUNTER FOR V(R)
          V = 0.0D0
C
C         INTEGRATE OVER CHARGE SOURCE
          DO M=0,NSRC
C
C           SET SOURCE RADIUS
            S  = S0 + HS*DFLOAT(M)
C
C           CHI FUNCTION ARGUMENTS
            XM = 2.0D0*DABS(R-S)/CMPF
            XP = 2.0D0*DABS(R+S)/CMPF
C
C           COMPONENTS OF INTEGRAND
            IF(0.5D0*(XM+XP).LT.XSPL) THEN
              CALL SPLNINT(XS,Y2S,D2S,XM,CM,NLW)
              CALL SPLNINT(XS,Y2S,D2S,XP,CP,NLW)    
            ELSE
              CALL SPLNINT(XB,Y2B,D2B,XM,CM,NUP)
              CALL SPLNINT(XB,Y2B,D2B,XP,CP,NUP)
            ENDIF
C
C           PERFORM THE MAPPING CHI(X) = SPLINE(X)*EXP(-X)
            CM = CM*DEXP(-XM)
            CP = CP*DEXP(-XP)
C
C           CONTRIBUTION TO INTEGRAND
            Z = S*rhonuc(NMDL(IZ),IZ,S)*(CM-CP)
            V = V + 5.0D0*HS*EXTINT11(Z,M,NSRC)/2.99376D+5
C          
          ENDDO
C
C         MULTIPLICATIVE FACTORS FOR VUEH(N)
          VUEH(N) =-2.0D0*CMPF*V/(3.0D0*CV*R)
C
        ENDIF
C
      ENDDO
      
      RETURN
C
C     FILE NAME AND TITLES
      XOUT   = 'Vueh'
      TITLE  = 'Uehling potential'
      XAXIS  = 'r [fm]'
      YAXIS  = 'V(r) [C/fm]'
      KEY(1) = 'Gaussian'
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C       BEGIN LOOP OVER ALL POINTS
        DO N=0,300
C
C         PRINT DISTANCE ALONG PATH AND FOUR-CURRENT
          WRITE(8, *) RAD(N)*CFM,VUEH(N)/CFM
C
        ENDDO
C
C     CLOSE DATA FILE
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
C
C
      SUBROUTINE EUEHGRD(IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    EEEEEEEE UU    UU EEEEEEEE HH    HH  GGGGGG  RRRRRRR  DDDDDDD     C
C    EE       UU    UU EE       HH    HH GG    GG RR    RR DD    DD    C
C    EE       UU    UU EE       HH    HH GG       RR    RR DD    DD    C
C    EEEEEE   UU    UU EEEEEE   HHHHHHHH GG       RR    RR DD    DD    C
C    EE       UU    UU EE       HH    HH GG   GGG RRRRRRR  DD    DD    C
C    EE       UU    UU EE       HH    HH GG    GG RR    RR DD    DD    C
C    EEEEEEEE  UUUUUU  EEEEEEEE HH    HH  GGGGGG  RR    RR DDDDDDD     C
C                                                                      C
C -------------------------------------------------------------------- C
C  EUEHGRD EVALUATES THE UEHLING ELECTRIC FIELD DUE TO A NUCLEAR       C
C  CHARGE SOURCE (GAUSSIAN SHAPE) FOR A SERIES OF RADII.               C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(1)
C
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/VPOL/RAD(0:NRAD),VUEH(0:NRAD),EUEH(0:NRAD),PUEH(0:NRAD)
      
C     TAKE NUMERICAL DERIVATIVE
      H = RAD(1)-RAD(0)
      DO N=1,NRAD-1
        EUEH(N) =-0.5D0*(VUEH(N+1)-VUEH(N-1))/H
      ENDDO
C
C     SPECIAL CASES (PHYSICALLY-MOTIVATED)
      EUEH(0   ) = 0.0D0
      EUEH(NRAD) = 0.0D0
C
C     THIS FIXES THE INCONSISTENCY BETWEEN V(0) AND V(R) METHODS
      EUEH(1) = 0.5D0*EUEH(2)
C
      RETURN
C
C     FILE NAME AND TITLES
      XOUT   = 'Eueh'
      TITLE  = 'Uehling electric field E_{r}(r)'
      XAXIS  = 'r [fm]'
      YAXIS  = 'E_{r}(r) [C/fm^2]'
      KEY(1) = 'Gaussian'
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C       BEGIN LOOP OVER ALL POINTS
        DO N=0,300
C
C         PRINT DISTANCE ALONG PATH AND FOUR-CURRENT
          WRITE(8, *) RAD(N)*CFM,EUEH(N)/(CFM*CFM)
C
        ENDDO
C
C     CLOSE DATA FILE
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
C
C
      SUBROUTINE PUEHGRD(IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    PPPPPPP  UU    UU EEEEEEEE HH    HH  GGGGGG  RRRRRRR  DDDDDDD     C
C    PP    PP UU    UU EE       HH    HH GG    GG RR    RR DD    DD    C
C    PP    PP UU    UU EE       HH    HH GG       RR    RR DD    DD    C
C    PPPPPPP  UU    UU EEEEEE   HHHHHHHH GG       RR    RR DD    DD    C
C    PP       UU    UU EE       HH    HH GG   GGG RRRRRRR  DD    DD    C
C    PP       UU    UU EE       HH    HH GG    GG RR    RR DD    DD    C
C    PP        UUUUUU  EEEEEEEE HH    HH  GGGGGG  RR    RR DDDDDDD     C
C                                                                      C
C -------------------------------------------------------------------- C
C  PUEHGRD EVALUATES THE UEHLING CHARGE DISTRIBUTION FIELD DUE TO A    C
C  NUCLEAR CHARGE SOURCE (GAUSSIAN SHAPE) FOR A SERIES OF RADII.       C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5  NMDL
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(1)
C
      DIMENSION ETMP(0:NRAD),RK(7)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/VPOL/RAD(0:NRAD),VUEH(0:NRAD),EUEH(0:NRAD),PUEH(0:NRAD)
C
C     PRE-MULTIPLICATION (WEIGH BY R^2)
      DO N=0,NRAD
        ETMP(N) = RAD(N)*RAD(N)*EUEH(N)
      ENDDO
C
C     FINITE DIFFERENCE FIRST DERIVATIVE
      HR = RAD(1)-RAD(0)
      DO N=1,NRAD-1
        PUEH(N) = 0.5D0*(ETMP(N+1)-ETMP(N-1))/HR
      ENDDO
C
C     SPECIAL CASES (PHYSICALLY-MOTIVATED)
      PUEH(0   ) = 0.0D0
      PUEH(NRAD) = 0.0D0
C
C     SEARCH FOR START OF POSITIVE CHARGE
      RZR = 0.0D0
      NPS = 0
      DO N=0,NRAD
        IF(PUEH(N).GT.0.0D0) THEN
          NPS = N
          RZR = (PUEH(N)*RAD(N-1)-PUEH(N-1)*RAD(N))/(PUEH(N)-PUEH(N-1))
          RZR = RZR*CFM
          GOTO 10
        ENDIF
      ENDDO
10    CONTINUE
C      
C     NUMBER OF INTEGRATION STEPS MUST BE DIVISIBLE BY 10
      NCLN = NRAD-NPS-MOD(NRAD-NPS,10)
C
C     INTEGRATE UP OVER POSITIVE PARTS
      QP = 0.0D0
      DO N=NPS,NPS+NCLN
        QP = QP + 5.0D0*HR*EXTINT11(PUEH(N),N-NPS,NCLN)/2.99376D+5
      ENDDO
C
C     CALCULATION OF MOMENTS
C     WRITE(6, *) 'Radial moments...'
      RK(1) = VUEH(0)/CFM
      DO K=-2,3
C
        HR = RAD(1)-RAD(0)
        RK(K+4) = 0.0D0
        DO N=0,NRAD
          IF(N.EQ.0) THEN
            Z = 0.0D0
          ELSE
            Z  = (RAD(N)**K)*PUEH(N)
          ENDIF
          RK(K+4) = RK(K+4) + 5.0D0*HR*EXTINT11(Z,N,NRAD)/2.99376D+5
        ENDDO
C
        RK(K+4) = RK(K+4)*(CFM**K)
C
      ENDDO
C
      AVAL = (RNUC(IZ)*CFM - 0.570)/0.836D0
      AVAL = AVAL**3
      
      WRITE(6, *) RNUC(IZ)*CFM,AVAL,(RK(I),I=1,7),QP,RZR
C
C     FILE NAME AND TITLES
      XOUT   = 'Pueh_frm'
      TITLE  = 'Weighted Uehling radial charge density rho(r)'
      XAXIS  = 'r [fm]'
      YAXIS  = 'r^2rho(r) [C/fm^3]'
      KEY(1) = 'Gaussian'
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C       BEGIN LOOP OVER ALL POINTS
        DO N=0,300
C
C         PRINT DISTANCE ALONG PATH AND FOUR-CURRENT
          WRITE(8, *) RAD(N)*CFM,PUEH(N)/CFM
C
        ENDDO
C
C     CLOSE DATA FILE
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
C
C
      SUBROUTINE PUEHUNI(IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      PPPPPPP  UU    UU EEEEEEEE HH    HH UU    UU NN    NN IIII      C
C      PP    PP UU    UU EE       HH    HH UU    UU NNN   NN  II       C
C      PP    PP UU    UU EE       HH    HH UU    UU NNNN  NN  II       C
C      PPPPPPP  UU    UU EEEEEE   HHHHHHHH UU    UU NN NN NN  II       C
C      PP       UU    UU EE       HH    HH UU    UU NN  NNNN  II       C
C      PP       UU    UU EE       HH    HH UU    UU NN   NNN  II       C
C      PP        UUUUUU  EEEEEEEE HH    HH  UUUUUU  NN    NN IIII      C
C                                                                      C
C -------------------------------------------------------------------- C
C  PUEHLNG EVALUATES THE UEHLING CHARGE DISTRIBUTION FIELD DUE TO A    C
C  NUCLEAR CHARGE SOURCE (GAUSSIAN SHAPE) FOR A SERIES OF RADII.       C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5  NMDL
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(1)
C
      DIMENSION ETMP(0:NRAD),RK(7)
      DIMENSION C0(7),CNS(7),CND(7)
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/VPOL/RAD(0:NRAD),VUEH(0:NRAD),EUEH(0:NRAD),PUEH(0:NRAD)
C
C     MULTIPLIER FOR COMPTON WAVELENGTH (USE FOR MUON OR TAUON FIELD)
      CMPF = 1.0D0*CMPW
C
C     RADIAL PLOTTING INFORMATION
      R0 = 0.0D0/CFM
      RN = 2000.0D0/CFM
      HR = (RN-R0)/DFLOAT(NRAD)
C
C     CHARGE SOURCE EXTENT AND UNITLESS PARAMETER
      SM = DSQRT(5.0D0/3.0D0)*RNUC(IZ)
      X  = 2.0D0*SM/CMPF
C
C     SKIP THIS STEP
      GOTO 50
C
C     GENERATE RADIAL GRID AND CHARGE DENSITY
      DO N=0,NRAD
C
C       SET RADIUS
        R = R0 + HR*DFLOAT(N)
        RAD(N) = R
C
C       INITIALISE COUNTER FOR RHO(R)
        RHO = 0.0D0
C
C       MULTIPLICATIVE FACTOR
        RMLT =-CMPF*R/(2.0D0*PI*CV*SM*SM*SM)
        IF(R.LE.SM) THEN
C
          XM = 2.0D0*(SM-R)/CMPF
          XP = 2.0D0*(SM+R)/CMPF
C
C         CHI_2(X) TERMS
          CHIM = CHIFNC(XM,2,0)
          CHIP = CHIFNC(XP,2,0)
C
          RHO = RHO + CHIM-CHIP
C
C         CHI_1(X) TERMS
          CHIM = CHIFNC(XM,1,0)
          CHIP = CHIFNC(XP,1,0)
C
          RHO = RHO + X*(CHIM-CHIP)
C
        ELSE
C
          XM = 2.0D0*(R-SM)/CMPF
          XP = 2.0D0*(R+SM)/CMPF
C
C         CHI_2(X) TERMS
          CHIM = CHIFNC(XM,2,0)
          CHIP = CHIFNC(XP,2,0)
C
          RHO = RHO + CHIM-CHIP
C
C         CHI_1(X) TERMS
          CHIM = CHIFNC(XM,1,0)
          CHIP = CHIFNC(XP,1,0)
C
          RHO = RHO - X*(CHIM+CHIP)
C
        ENDIF
        
        PUEH(N) = RMLT*RHO
C        
      ENDDO
C
C     FILE NAME AND TITLES
      XOUT   = 'Pueh_uni_exct'
      TITLE  = 'Weighted Uehling radial charge density rho(r)'
      XAXIS  = 'r [fm]'
      YAXIS  = 'r^2rho(r) [C/fm^3]'
      KEY(1) = 'Gaussian'
C
C     OPEN DATA FILE
      OPEN(UNIT=8,FILE='plots/'//TRIM(XOUT)//'.dat',STATUS='REPLACE')
C
C       BEGIN LOOP OVER ALL POINTS
        DO N=0,300
C
C         PRINT DISTANCE ALONG PATH AND FOUR-CURRENT
          WRITE(8, *) RAD(N)*CFM,PUEH(N)/CFM
C
        ENDDO
C
C     CLOSE DATA FILE
      CLOSE(UNIT=8)
C
C     GENERATE A GNUPLOT MAKE FILE
      CALL GNULINE(XOUT,TITLE,XAXIS,YAXIS,1,KEY)
C
C     EXECUTE GNUPLOT COMMAND IN TERMINAL
      CALL SYSTEM('gnuplot plots/'//TRIM(XOUT)//'.gnuplot')
      CALL SYSTEM('xdg-open plots/'//TRIM(XOUT)//'.pdf')
C
50    CONTINUE
C
C     CHARGE SOURCE EXTENT AND UNITLESS PARAMETER
      SM = DSQRT(5.0D0/3.0D0)*RNUC(IZ)
      X  = 2.0D0*SM/CMPF
C
C     CHI_N(0) VALUES
      C0(2) = 9.0D0*PI/32.0D0
      C0(3) = 0.4D0
      C0(4) = 5.0D0*PI/64.0D0
      C0(5) = 6.0D0/35.0D0
      C0(6) = 21.0D0*PI/512.0D0
      C0(7) = 32.0D0/315.0D0
C
C     CHI_N(X) AND CHI_N(2X) VALUES
      DO N=2,7
        CNS(N) = CHIFNC(      X,N,0)
        CND(N) = CHIFNC(2.0D0*X,N,0)
      ENDDO
C
C     SEARCH FOR START OF POSITIVE CHARGE
      RZR = DSQRT(5.0D0/3.0D0)*RNUC(IZ)*CFM
C
C     TOTAL CHARGE IN EACH REGION (MAGNITUDE)
c     ***
      QP =-C0(4) + CND(4) + 2.0D0*X*CND(3) + X*X*C0(2) + X*X*CND(2)
      QP = QP/(CV*PI*X*X*X)
C
C     CALCULATION OF MOMENTS
C
C     <R^-2> UNDEFINED FOR NOW
      RK(2) = 0.0D0
C
C     <R^-1>
      RK(3) =-C0(3) + CNS(3) + X*CNS(2)
      RK(3) = 4.0D0*RK(3)/(PI*X*X*X*CFM)
C
C     <R^0>
      RK(4) = 0.0D0
C
C     <R^1>
      RK(5) =-C0(5) + CNS(5) + X*CNS(4) + 0.5D0*X*X*C0(3)
      RK(5) = 2.0D0*RK(5)*CMPF*CFM/(CV*PI*X*X*X)
C
C     <R^2>
      RK(6) = 2.0D0*((CMPF*CFM)**2)/(5.0D0*CV*PI)
C
C     <R^3>
c     ***
      RK(7) =-C0(7) + CNS(7) + X*CNS(6) + 0.5D0*X*X*C0(5)
     &      + 0.125D0*X*X*X*X*C0(3)
C
C      WRITE(6, *) RNUC(IZ)*CFM,-C0(7),CNS(7),X*CNS(6),0.5D0*X*X*C0(5),
C     &           0.125D0*X*X*X*X*C0(3),rk(7)
C      WRITE(7, *) RNUC(IZ)*CFM,-C0(7),CNS(7),X*CNS(6),0.5D0*X*X*C0(5),
C     &           0.125D0*X*X*X*X*C0(3),rk(7)
C      IF(RNUC(IZ)*CFM.EQ.2.0D0) STOP

      RK(7) = 6.0D0*RK(7)/(CV*PI)
      RK(7) = RK(7)*((CMPF*CFM/X)**3)
C
C     VUEH(0) DEFINED AS IDENTICAL TO <R^-1>
      RK(1) = RK(3)
C
      AVAL = (RNUC(IZ)*CFM - 0.570D0)/0.836D0
      AVAL = AVAL**3
      
      WRITE(6, *) RNUC(IZ)*CFM,AVAL,(RK(I),I=1,7),QP,RZR
      WRITE(7, *) RNUC(IZ)*CFM,AVAL,(RK(I),I=1,7),QP,RZR
C
      RETURN
      END
C
C
      FUNCTION VASYMP(IZ,R,CMPF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       VV    VV    AA     SSSSSS  YY    YY MM       MM PPPPPPP        C
C       VV    VV   AAAA   SS    SS YY    YY MMM     MMM PP    PP       C
C       VV    VV  AA  AA  SS        YY  YY  MMMM   MMMM PP    PP       C
C       VV    VV AA    AA  SSSSSS    YYYY   MM MM MM MM PP    PP       C
C        VV  VV  AAAAAAAA       SS    YY    MM  MMM  MM PPPPPPP        C
C         VVVV   AA    AA SS    SS    YY    MM   M   MM PP             C
C          VV    AA    AA  SSSSSS     YY    MM       MM PP             C
C                                                                      C
C -------------------------------------------------------------------- C
C  VASYMP IS THE ASYMPTOTIC POTENTIAL FOR NUCLEAR CENTRE IZ.           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5 NMDL
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
C
      RP2  = RNUC(IZ)*RNUC(IZ)
      ZETA = 1.5D0/RP2
C
      GSS = 1.0D0/(ZETA*CMPF*CMPF)
      EFC =-2.0D0*R/CMPF
C
C     INTEGRATION PARAMETERS
      NTRM = 200000
      TBEG = 1.0D0
      TEND = 4.0D3
      HT   = (TEND-TBEG)/DFLOAT(NTRM)
C
      Y = 0.0D0
      DO N=0,NTRM
        T  = TBEG + HT*DFLOAT(N)
        X1 = 1.0D0/T
        X2 = DSQRT(1.0D0 - X1*X1)
        X3 = 1.0D0 + 0.5D0*X1*X1
        X4 = DEXP(GSS*X1*X1)
        X5 = DEXP(EFC*T)
        Y  = Y + EXTINT11(X1*X2*X3*X4*X5,N,NTRM)
      ENDDO
      VASYMP = 5.0D0*HT*Y/2.99376D+5
C
      RETURN
      END
C
C
      FUNCTION QUNIFRM(IZ,S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     QQQQQQ   UU    UU NN    NN IIII FFFFFFFF RRRRRRR  MM       MM    C
C    QQ    QQ  UU    UU NNN   NN  II  FF       RR    RR MMM     MMM    C
C   QQ      QQ UU    UU NNNN  NN  II  FF       RR    RR MMMM   MMMM    C
C   QQ      QQ UU    UU NN NN NN  II  FFFFFF   RR    RR MM MM MM MM    C
C   QQ      QQ UU    UU NN  NNNN  II  FF       RRRRRRR  MM  MMM  MM    C
C    QQ    QQ  UU    UU NN   NNN  II  FF       RR    RR MM   M   MM    C
C     QQQQQQ QQ UUUUUU  NN    NN IIII FF       RR    RR MM       MM    C
C                                                                      C
C -------------------------------------------------------------------- C
C  QUNIFRM IS THE RADIAL CHARGE DENSITY GENERATED BY A PROTON WITH     C
C  POSITIVE UNIT CHARGE, MEAN SQUARE RADIUS RPAU, AND WHOSE CHARGE     C
C  DISTRIBUTION IS A UNIFORM BALL.                                     C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5 NMDL
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
      B   = DSQRT(5.0D0/3.0D0)*RNUC(IZ)
      IF(S.LE.B) THEN
        B3      = B*B*B
        QUNIFRM = 0.75/(PI*B3)
      ELSE
        QUNIFRM = 0.0D0
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION QGAUSS0(IZ,S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     QQQQQQ    GGGGGG     AA    UU    UU  SSSSSS   SSSSSS   000000    C
C    QQ    QQ  GG    GG   AAAA   UU    UU SS    SS SS    SS 00   000   C
C   QQ      QQ GG        AA  AA  UU    UU SS       SS       00  0000   C
C   QQ      QQ GG       AA    AA UU    UU  SSSSSS   SSSSSS  00 00 00   C
C   QQ      QQ GG   GGG AAAAAAAA UU    UU       SS       SS 0000  00   C
C    QQ    QQ  GG    GG AA    AA UU    UU SS    SS SS    SS 000   00   C
C     QQQQQQ QQ GGGGGG  AA    AA  UUUUUU   SSSSSS   SSSSSS   000000    C
C                                                                      C
C -------------------------------------------------------------------- C
C  QGAUSS0 IS THE RADIAL CHARGE DENSITY GENERATED BY A PROTON WITH     C
C  POSITIVE UNIT CHARGE, MEAN SQUARE RADIUS RPAU, AND WHOSE CHARGE     C
C   DISTRIBUTION IS A SIMPLE GAUSSIAN.                                 C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*5 NMDL
C
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,0:MFT),XNUC(MCT,0:MFT),NNUC(MCT),NMDL(MCT)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
      RP2 = RNUC(IZ)*RNUC(IZ)
      AGS = 1.5D0/RP2
      TRM = AGS/PI     
C
      QGAUSS0 = TRM*DSQRT(TRM)*DEXP(-AGS*S*S)
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
C**********************************************************************C
C ==================================================================== C
C  [16] VACPOL: POTENTIALS BY CONVOLUTION                              C
C ==================================================================== C
C  ROUTINES AND FUNCTIONS:                                             C
C -------------------------------------------------------------------- C
C   [A] GENGRID: ASYMPTOTIC FORM FOR UEHLING DENSITY AT LARGE R.       C
C   [B] CLMDIAG: GENERATES ARRAY OF COULOMB SELF-OVERLAPS.             C
C   [C] BRTDIAG: GENERATES ARRAY OF BREIT SELF-OVERLAPS.               C
C**********************************************************************C
C
C
      SUBROUTINE GENGRID
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       GGGGGG  EEEEEEEE NN    NN  GGGGGG  RRRRRRR  IIII DDDDDDD       C
C      GG    GG EE       NNN   NN GG    GG RR    RR  II  DD    DD      C
C      GG       EE       NNNN  NN GG       RR    RR  II  DD    DD      C
C      GG       EEEEEE   NN NN NN GG       RR    RR  II  DD    DD      C
C      GG   GGG EE       NN  NNNN GG   GGG RRRRRRR   II  DD    DD      C
C      GG    GG EE       NN   NNN GG    GG RR    RR  II  DD    DD      C
C       GGGGGG  EEEEEEEE NN    NN  GGGGGG  RR    RR IIII DDDDDDD       C
C                                                                      C
C -------------------------------------------------------------------- C
C  GENGRID GENERATES RADIAL AND MOMENTUM GRIDS FOR VACPOL QUANTITIES.  C
C  ▶ RADIUS R IS GIVEN IN FM                                           C
C  ▶ MOMENTUM Q IS GIVEN IN FM^-1                                      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/RGRD/RFM(0:NRVP),REP(0:NRVP),RMIN,RMID,RMAX
      COMMON/QGRD/QFM(0:NQVP),QEP(0:NQVP),QMIN,QMID,QMAX
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
C
C     DIVIDE THE RADIAL GRID INTO A UNIFORM AND EXPONENTIAL SECTION
      RMIN = 0.0D0/CFM
      RMID = 4.0D1/CFM
      RMAX = 6.0D3/CFM
C
      NRLN = NRVP/3
      NREX = NRVP-NRLN
C
C     UNIFORM RADIAL GRID (IN AU)
      REPUN = (RMID-RMIN)/DFLOAT(NRLN)
      DO N=0,NRLN
        RFM(N) = REPUN*DFLOAT(N)
        REP(N) = REPUN
      ENDDO
C
C     EXPONENTIAL RADIAL GRID (IN AU)
      HEPEX = DLOG(RMAX/RMID)/DFLOAT(NREX)
      DO N=NRLN+1,NRVP
        I = N-NRLN
        RFM(N) = RMID*DEXP(HEPEX*DFLOAT(I))
        REP(N) = HEPEX*RFM(N)
      ENDDO
C
C     DIVIDE THE MOMENTUM GRID INTO A UNIFORM AND EXPONENTIAL SECTION
      QMIN = 0.0D0*CFM
      QMID = 1.5D1*CFM
      QMAX = 3.0D1*CFM
C
C     UNIFORM MOMENTUM GRID (IN AU^-1)
      QEPUN = (QMID-QMIN)/DFLOAT(NQVP)
      DO I=0,NQVP
        QFM(I) = QEPUN*DFLOAT(I)
        QEP(I) = QEPUN
      ENDDO
C
C      COULD USE AN EXPONENTIALLY TAILING GRID AS WELL
C      I DECIDED TO ABANDON THIS BECAUSE FERMI ORIGIN WAS THE PROBLEM
CC
C      NQLN = 2*NQ/3
C      NQEX = NQ-NQLN
CC
CC     UNIFORM MOMENTUM GRID (IN FM^-1)
C      QEPUN = (QMID-QMIN)/DFLOAT(NQLN)
C      DO N=0,NQLN
C        QFM(N) = QEPUN*DFLOAT(N)
C        QEP(N) = QEPUN
C      ENDDO
CC
CC     EXPONENTIAL MOMENTUM GRID (IN FM^-1)
C      HEPEX = DLOG(QMAX/QMID)/DFLOAT(NQEX)
C      DO N=NQLN+1,NQ
C        I = N-NQLN
C        QFM(N) = QMID*DEXP(HEPEX*DFLOAT(I))
C        QEP(N) = HEPEX*QFM(N)
C      ENDDO
C
      RETURN
      END
C
C
      FUNCTION PIUEHF(Q,CMPF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          PPPPPPP  IIII UU    UU EEEEEEEE HH    HH FFFFFFFF           C
C          PP    PP  II  UU    UU EE       HH    HH FF                 C
C          PP    PP  II  UU    UU EE       HH    HH FF                 C
C          PP    PP  II  UU    UU EEEEEE   HHHHHHHH FFFFFF             C
C          PPPPPPP   II  UU    UU EE       HH    HH FF                 C
C          PP        II  UU    UU EE       HH    HH FF                 C
C          PP       IIII  UUUUUU  EEEEEEEE HH    HH FF                 C
C                                                                      C
C -------------------------------------------------------------------- C
C  PIUEHF IS THE SPECTRAL FUNCTION FOR THE UEHLING INTERACTION THAT    C
C  INVOLVES A FERMION WITH ITS REDUCED COMPTON WAVELENGTH CMPF.        C
C**********************************************************************C
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DIMENSIONLESS VARIABLE U
      U  = Q*CMPF
      U2 = U*U
C
C     SWITCHING PARAMETERS
      UPOW = 1.0D+0
      UPAS = 3.0D+5
C
      IF(U.LT.UPOW) THEN
C     POWER SERIES FOR SMALL U
        U4 = U2*U2
        U6 = U4*U2
        U8 = U6*U2
        PIUEHF = U2/1.0D1 - 3.0D0*U4/2.8D2 + U6/6.3D2 - U8/369.6D1
      ELSEIF(U.GE.UPOW.AND.U.LT.UPAS) THEN
C     FULL EXPRESSION FOR MEDIUM U
        TRM = (U2-2.0D0)*DSQRT(U2+4.0D0)*DASINH(0.5D0*U)/(U2*U)
        PIUEHF =-5.0D0/6.0D0 + 2.0D0/U2 + TRM
      ELSEIF(U.GE.UPAS) THEN
C     ULTRA-RELATIVISTIC CASE FOR LARGE U
        A =-5.0D0/6.0D0 + DLOG(U)
        B = 3.0D0
        C =-3.0D0/2.0D0 + 6.0D0*DLOG(U)
        PIUEHF = A + B/U2 + C/(U2*U2)
        PIUEHF = 0.0D0
      ENDIF
C
C     AMPLITUDE
      AMP = 2.0D0/(3.0D0*PI*CV)
      PIUEHF = AMP*PIUEHF
C
      RETURN
      END
C
C
      FUNCTION PIUEHB(Q,CMPF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          PPPPPPP  IIII UU    UU EEEEEEEE HH    HH BBBBBBB            C
C          PP    PP  II  UU    UU EE       HH    HH BB    BB           C
C          PP    PP  II  UU    UU EE       HH    HH BB    BB           C
C          PP    PP  II  UU    UU EEEEEE   HHHHHHHH BBBBBBB            C
C          PPPPPPP   II  UU    UU EE       HH    HH BB    BB           C
C          PP        II  UU    UU EE       HH    HH BB    BB           C
C          PP       IIII  UUUUUU  EEEEEEEE HH    HH BBBBBBB            C
C                                                                      C
C -------------------------------------------------------------------- C
C  PIUEHV IS THE SPECTRAL FUNCTION FOR THE UEHLING INTERACTION THAT    C
C  INVOLVES A BOSON WITH ITS REDUCED COMPTON WAVELENGTH CMPF.          C
C -------------------------------------------------------------------- C
C  CHARGE ZNUC IS NORMALISED IN THIS FORMULATION (RESTORE IN BERTHA).  C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DIMENSIONLESS VARIABLE U
      U  = Q*CMPF
      U2 = U*U
C
C     SWITCHING PARAMETERS
      UPOW = 1.0D-1
      UPAS = 3.0D+5
C
      IF(U.LT.UPOW) THEN
C     POWER SERIES FOR SMALL U
        U4 = U2*U2
        U6 = U4*U2
        U8 = U6*U2
        PIUEHB = U2/4.0D1 - 3.0D0*U4/5.6D2 + U6/50.4D2 - U8/369.6D2
      ELSEIF(U.GE.UPOW.AND.U.LT.UPAS) THEN
C     FULL EXPRESSION FOR MEDIUM U
        TRM = (U2+4.0D0)*DSQRT(U2+4.0D0)*DASINH(0.5D0*U)/(2.0D0*U2*U)
        PIUEHB = -2.0D0/3.0D0 - 2.0D0/U2 + TRM
      ELSEIF(U.GE.UPAS) THEN
C     ULTRA-RELATIVISTIC CASE FOR LARGE U
        A =-2.0D0/3.0D0 + 0.5D0*DLOG(U)
        B =-3.0D0/2.0D0 + 3.0D0*DLOG(U)
        C = 9.0D0/4.0D0 + 3.0D0*DLOG(U)
        PIUEHB = A + B/U2 + C/(U2*U2)
        PIUEHB = 0.0D0
      ENDIF
C
C     AMPLITUDE
      AMP    = 2.0D0/(3.0D0*PI*CV)
      PIUEHB = AMP*PIUEHB
C
      RETURN
      END
C
C
      FUNCTION PIKSBF(Q,CMPF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          PPPPPPP  IIII KK    KK  SSSSSS  BBBBBBB  FFFFFFFF           C
C          PP    PP  II  KK   KK  SS    SS BB    BB FF                 C
C          PP    PP  II  KK  KK   SS       BB    BB FF                 C
C          PP    PP  II  KKKKK     SSSSSS  BBBBBBB  FFFFFF             C
C          PPPPPPP   II  KK  KK         SS BB    BB FF                 C
C          PP        II  KK   KK  SS    SS BB    BB FF                 C
C          PP       IIII KK    KK  SSSSSS  BBBBBBB  FF                 C
C                                                                      C
C -------------------------------------------------------------------- C
C  PIKSBF IS THE SPECTRAL FUNCTION FOR THE KÄLLÉN-SABRY INTERACTION    C
C  THAT INVOLVES A FERMION WITH ITS REDUCED COMPTON WAVELENGTH CMPF.   C
C**********************************************************************C
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     DIMENSIONLESS VARIABLE U
      U  = Q*CMPF
      U2 = U*U
C
C     SWITCHING PARAMETERS
      UPOW = 1.20D+0
      UPAS = 3.00D+3
C
      IF(U.LT.UPOW) THEN
C     POWER SERIES FOR SMALL U
        U4 = U2*U2
        U6 = U4*U2
        U8 = U6*U2
        PIKSBF = 41.0D0*U2/162.0D0
     &         - 401.0D0*U4/10800.0D0
     &         + 54919.0D0*U6/7938000.0D0
     &         - 25993.0D0*U8/16329600.0D0
      ELSEIF(U.GE.UPOW.AND.U.LT.UPAS) THEN
C     FULL EXPRESSION FOR MEDIUM U
        DL2 = 1.0D0 + 4.0D0/U2
        DLT = DSQRT(DL2)
        DLM = 1.0D0-DLT
        DLP = 1.0D0+DLT
        DRT = DABS(DLM)/DLP
        AD3 = DABS(1.0D0-DL2)*DABS(1.0D0-DL2)*DABS(1.0D0-DL2)
        T1  =-13.0D0/108.0D0 + 11.0D0*DL2/72.0D0 - DL2*DL2/3.0D0
        T2  = 19.0D0/24.0D0 - 55.0D0*DL2/72.0D0 + DL2*DL2/3.0D0
        T3  = DLT*DLOG(1.0D0/DRT)
        T4  = 33.0D0/32.0D0 + 23.0D0*DL2/16.0D0 
     &      - 23.0D0*DL2*DL2/32.0D0 + DL2*DL2*DL2/12.0D0
        T5  = DLOG(1.0D0/DRT)*DLOG(1.0D0/DRT) - PI*PI*THETAKS(DLM)
        T6  = PHIKS(DLM/DLP) + 2.0D0*PHIKS(-DLM/DLP)
     &      + 0.25D0*PI*PI
     &      - 0.75D0*PI*PI*THETAKS(DLM)
     &      - 0.75D0*DLOG(1.0D0/DRT)*DLOG(1.0D0/DRT)
     &      + 0.50D0*DLOG(1.0D0/DRT)*DLOG(64.0D0*DL2*DL2/AD3)
        T7  = DLT*(3.0D0-DL2)
        T8  = 3.0D0 + 2.0D0*DL2 - DL2*DL2
        T9  = FKS(DL2) + 1.5D0*GKS(DL2) - HKS(DL2)
        PIKSBF = T1 + T2*T3 - T4*T5 + T6*T7 + T8*T9
        PIKSBF =-PIKSBF/3.0D0
      ELSEIF(U.GE.UPAS) THEN
C     ULTRA-RELATIVISTIC CASE FOR LARGE U
        A =-5.0D0/6.0D0 + DLOG(U)
        B = 3.0D0
        C =-3.0D0/2.0D0 + 6.0D0*DLOG(U)
        PIKSBF = A + B/U2 + C/(U2*U2)
        PIKSBF = 0.0D0
      ENDIF
C
C     AMPLITUDE
      PIKSBF = PIKSBF/(PI*PI*CV*CV)
C
      RETURN
      END
C
C
      FUNCTION THETAKS(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                       θθθθθθ  KK    KK  SSSSSS                       C
C                      θθ    θθ KK   KK  SS    SS                      C
C                      θθ    θθ KK  KK   SS                            C
C                      θθθθθθθθ KKKKK     SSSSSS                       C
C                      θθ    θθ KK  KK         SS                      C
C                      θθ    θθ KK   KK  SS    SS                      C
C                       θθθθθθ  KK    KK  SSSSSS                       C
C                                                                      C
C -------------------------------------------------------------------- C
C  THETAKS IS A SPECIAL FUNCTION NEEDED FOR THE KÄLLÉN-SABRY FUNCTION. C
C**********************************************************************C
C
      THETAKS = 0.5D0*(1.0D0 + X/DABS(X))
C
      RETURN
      END
C
C
      FUNCTION PHIKS(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                         ΦΦ    KK    KK  SSSSSS                       C
C                       ΦΦΦΦΦΦ  KK   KK  SS    SS                      C
C                      ΦΦ ΦΦ ΦΦ KK  KK   SS                            C
C                      ΦΦ ΦΦ ΦΦ KKKKK     SSSSSS                       C
C                      ΦΦ ΦΦ ΦΦ KK  KK         SS                      C
C                       ΦΦΦΦΦΦ  KK   KK  SS    SS                      C
C                         ΦΦ    KK    KK  SSSSSS                       C
C                                                                      C
C -------------------------------------------------------------------- C
C  PHIKS IS A SPECIAL FUNCTION NEEDED FOR THE KÄLLÉN-SABRY FUNCTION.   C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
C     EARLY EXIT CONDITION
      IF(DABS(X-1.0D0).LT.1.0D-6) RETURN
C
      TMIN = 1.0D0
      TMAX = X
C
      HSPC = (TMAX-TMIN)/DFLOAT(NTVP)
C
C     LOOP OVER NTVP EQUALLY-SPACED INTEGRATION POINTS
      PHIKS = 0.0D0
      DO K=0,NTVP
        RK = DFLOAT(K)
C
C       VALUE OF T
        T = TMIN + HSPC*RK
C
C       INTEGRAND TERMS
        C1 = DLOG(DABS(1.0D0+T))
        C2 = 1.0D0/T
C
C       ADD TO THE COUNTER
        PHIKS = PHIKS + EXTINT11(C1*C2,K,NTVP)
C
      ENDDO
C
      PHIKS = 5.0D0*HSPC*PHIKS/2.99376D+5
C
      RETURN
      END
C
C
      FUNCTION FKS(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                      FFFFFFFF KK    KK  SSSSSS                       C
C                      FF       KK   KK  SS    SS                      C
C                      FF       KK  KK   SS                            C
C                      FFFFFF   KKKKK     SSSSSS                       C
C                      FF       KK  KK         SS                      C
C                      FF       KK   KK  SS    SS                      C
C                      FF       KK    KK  SSSSSS                       C
C                                                                      C
C -------------------------------------------------------------------- C
C  FKS IS A SPECIAL FUNCTION NEEDED FOR THE KÄLLÉN-SABRY FUNCTION.     C
C  ENCOUNTERS A SINGULARITY AT EXACTLY T=-1.0D0, SO FOR NOW, SKIP IT.  C
C**********************************************************************C
C
      TMIN =-1.0D0
      TMAX =+1.0D0
      NTOT = 10000
      
      HSPC = (TMAX-TMIN)/DFLOAT(NTOT)
C
C     LOOP OVER NTVP EQUALLY-SPACED INTEGRATION POINTS
      FKS = 0.0D0
      DO K=1,NTOT
        RK = DFLOAT(K)
C
C       VALUE OF T
        T = TMIN + HSPC*RK
C
C       INTEGRAND TERMS
        IF(DABS(T).LT.1.0D-3) THEN
          C1 = 0.0D0
          TP = 1.0D0
          DO I=0,6
            C1 = C1 + DFLOAT(I+1)*TP
            TP = TP*T
          ENDDO
        ELSE
          C1 = DLOG(1.0D0+T)/T
        ENDIF
        C2 = DLOG(DABS(1.0D0-T*T/X))
C
C       ADD TO THE COUNTER
        FKS = FKS + EXTINT11(C1*C2,K,NTOT)
C
      ENDDO
C
      FKS = 5.0D0*HSPC*FKS/2.99376D+5
C
      RETURN
      END
C
C
      FUNCTION GKS(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                       GGGGGG  KK    KK  SSSSSS                       C
C                      GG    GG KK   KK  SS    SS                      C
C                      GG       KK  KK   SS                            C
C                      GG       KKKKK     SSSSSS                       C
C                      GG   GGG KK  KK         SS                      C
C                      GG    GG KK   KK  SS    SS                      C
C                       GGGGGG  KK    KK  SSSSSS                       C
C                                                                      C
C -------------------------------------------------------------------- C
C  GKS IS A SPECIAL FUNCTION NEEDED FOR THE KÄLLÉN-SABRY FUNCTION.     C
C  ENCOUNTERS A SINGULARITY AT EXACTLY T=-1.0D0, SO FOR NOW, SKIP IT.  C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      TMIN =-1.0D0
      TMAX =+1.0D0-1.0D-8
      NTOT = 10000
      
      HSPC = (TMAX-TMIN)/DFLOAT(NTOT)
C
C     LOOP OVER NTOT EQUALLY-SPACED INTEGRATION POINTS
      GKS = 0.0D0
      DO K=1,NTOT
        RK = DFLOAT(K)
C
C       VALUE OF T
        T = TMIN + HSPC*RK
C
C       INTEGRAND TERMS
        C1 = DLOG(0.5D0*DABS(1.0D0-T))
        C2 = 1.0D0/(T+1.0D0)
        IF(DABS(T+1.0D0).LT.1.0D-4) THEN
          T1 = T+1.0D0
          C1 = 0.5D0 + T1/8.0D0 + T1*T1/24.0D0 + T1*T1*T1/64.0D0
          C2 =-1.0D0
        ENDIF
        IF(DABS(T-1.0D0).LT.1.0D-4) THEN
          T1 = T+1.0D0
          C1 = 0.5D0 + T1/8.0D0 + T1*T1/24.0D0 + T1*T1*T1/64.0D0
          C2 =-1.0D0
        ENDIF
        C3 = DLOG(DABS(1.0D0 -T*T/X))
C
C       ADD TO THE COUNTER
        GKS = GKS + EXTINT11(C1*C2*C3,K,NTOT)
C
      ENDDO
C
      GKS = 5.0D0*HSPC*GKS/2.99376D+5
C
      RETURN
      END
C
C
      FUNCTION HKS(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                      HH    HH KK    KK  SSSSSS                       C
C                      HH    HH KK   KK  SS    SS                      C
C                      HH    HH KK  KK   SS                            C
C                      HHHHHHHH KKKKK     SSSSSS                       C
C                      HH    HH KK  KK         SS                      C
C                      HH    HH KK   KK  SS    SS                      C
C                      HH    HH KK    KK  SSSSSS                       C
C                                                                      C
C -------------------------------------------------------------------- C
C  HKS IS A SPECIAL FUNCTION NEEDED FOR THE KÄLLÉN-SABRY FUNCTION.     C
C  ENCOUNTERS A SINGULARITY AT EXACTLY T=-1.0D0, SO FOR NOW, SKIP IT.  C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      TMIN =-1.0D0
      TMAX =+1.0D0
      NTOT = 10000
C      
      HSPC = (TMAX-TMIN)/DFLOAT(NTOT)
C
C     LOOP OVER NTOT EQUALLY-SPACED INTEGRATION POINTS
      HKS = 0.0D0
      DO K=1,NTOT
        RK = DFLOAT(K)
C
C       VALUE OF T
        T = TMIN + HSPC*RK
C
C       INTEGRAND TERMS
        C1 = DLOG(DABS(T))
        C2 = 1.0D0/(T+1.0D0)
        IF(DABS(T).LT.1.0D-4) THEN
          C1 = 0.0D0
          C2 = 1.0D0
        ENDIF
        C3 = DLOG(DABS(1.0D0 -T*T/X))
C
C       ADD TO THE COUNTER
        HKS = HKS + EXTINT11(C1*C2*C3,K,NTOT)
C
      ENDDO
C
      HKS = 5.0D0*HSPC*HKS/2.99376D+5
C
      RETURN
      END
C
C
      FUNCTION SPHBJ0(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          SSSSSS  PPPPPPP  HH    HH BBBBBBB  JJJJJJJJ 000000          C
C         SS    SS PP    PP HH    HH BB    BB      JJ 00   000         C
C         SS       PP    PP HH    HH BB    BB      JJ 00  0000         C
C          SSSSSS  PP    PP HHHHHHHH BBBBBBB       JJ 00 00 00         C
C               SS PPPPPPP  HH    HH BB    BB      JJ 0000  00         C
C         SS    SS PP       HH    HH BB    BB JJ   JJ 000   00         C
C          SSSSSS  PP       HH    HH BBBBBBB   JJJJJ   000000          C
C                                                                      C
C -------------------------------------------------------------------- C
C  SHPHJ0 EVALUATES THE ZEROTH-ORDER REGULAR BESSEL FUNCTION, J0(X).   X
C**********************************************************************C
C
      IF(DABS(X).LT.1.0D-2) THEN
C     SMALL X
        X2 = X*X
        SPHBJ0 = 1.0D0 - X2/6.0D0 + X2*X2/120.0D0 - X2*X2*X2/5040.0D0
     &         + X2*X2*X2*X2/362880.0D0
      ELSE
C     LARGE X
        SPHBJ0 = DSIN(X)/X
      ENDIF
C
      RETURN
      END
C
C     
      SUBROUTINE PRIHNKT(GQ,R,RHO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      PPPPPPP  RRRRRRR  IIII HH    HH NN    NN KK    KK TTTTTTTT      C
C      PP    PP RR    RR  II  HH    HH NNN   NN KK   KK     TT         C
C      PP    PP RR    RR  II  HH    HH NNNN  NN KK  KK      TT         C
C      PP    PP RR    RR  II  HHHHHHHH NN NN NN KKKKK       TT         C
C      PPPPPPP  RRRRRRR   II  HH    HH NN  NNNN KK  KK      TT         C
C      PP       RR    RR  II  HH    HH NN   NNN KK   KK     TT         C
C      PP       RR    RR IIII HH    HH NN    NN KK    KK    TT         C
C                                                                      C
C -------------------------------------------------------------------- C
C  PRIHNKT CALCULATES A CHARGE DENSITY BY INVERSE HANKEL TRANSFORM,    C
C  STARTING WITH A RADIAL FORM FACTOR STORED IN GQ.                    C
C -------------------------------------------------------------------- C
C  ▶ RADIUS R IS GIVEN IN AU                                           C
C  ▶ MOMENTUM Q IS GIVEN IN AU^-1                                      C
C  ▶ NEED TO USE THE COMMON-STORED RADIAL GRID.                        C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION GQ(0:NQVP)
C
      COMMON/QGRD/QFM(0:NQVP),QEP(0:NQVP),QMIN,QMID,QMAX
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
      RHO = 0.0D0
      DO I=0,NQVP
        Q   = QFM(I)
        A   = SPHBJ0(Q*R)*QEP(I)
        RHO = RHO + EXTINT11(A*Q*Q*GQ(I),I,NQVP)
      ENDDO
      RHO = 0.5D0*5.0D0*RHO/(PI*PI*2.99376D+5)
C
      RETURN
      END
C
C     
      SUBROUTINE VRIHNKT(GQ,R,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      VV    VV RRRRRRR  IIII HH    HH NN    NN KK    KK TTTTTTTT      C
C      VV    VV RR    RR  II  HH    HH NNN   NN KK   KK     TT         C
C      VV    VV RR    RR  II  HH    HH NNNN  NN KK  KK      TT         C
C      VV    VV RR    RR  II  HHHHHHHH NN NN NN KKKKK       TT         C
C       VV  VV  RRRRRRR   II  HH    HH NN  NNNN KK  KK      TT         C
C        VVVV   RR    RR  II  HH    HH NN   NNN KK   KK     TT         C
C         VV    RR    RR IIII HH    HH NN    NN KK    KK    TT         C
C                                                                      C
C -------------------------------------------------------------------- C
C  VRIHNKT CALCULATES A POTENTIAL BY INVERSE HANKEL TRANSFORM,         C
C  STARTING WITH A RADIAL FORM FACTOR STORED IN GQ.                    C
C -------------------------------------------------------------------- C
C  ▶ RADIUS R IS GIVEN IN AU                                           C
C  ▶ MOMENTUM Q IS GIVEN IN AU^-1                                      C
C  ▶ NEED TO USE THE COMMON-STORED RADIAL GRID.                        C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION GQ(0:NQVP)
C
      COMMON/QGRD/QFM(0:NQVP),QEP(0:NQVP),QMIN,QMID,QMAX
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
      V = 0.0D0
      DO I=0,NQVP
        Q = QFM(I)
        A = SPHBJ0(Q*R)*QEP(I)
        V = V + EXTINT11(A*GQ(I),I,NQVP)
      ENDDO
      V =-2.0D0*5.0D0*V/(PI*2.99376D+5)
C
      RETURN
      END
C
C
      SUBROUTINE KSAUXFN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    KK    KK SSSSSS     AA    UU    UU XX    XX FFFFFFFF NN    NN     C
C    KK   KK SS    SS   AAAA   UU    UU  XX  XX  FF       NNN   NN     C
C    KK  KK  SS        AA  AA  UU    UU   XXXX   FF       NNNN  NN     C
C    KKKKK    SSSSSS  AA    AA UU    UU    XX    FFFFFF   NN NN NN     C
C    KK  KK        SS AAAAAAAA UU    UU   XXXX   FF       NN  NNNN     C
C    KK   KK SS    SS AA    AA UU    UU  XX  XX  FF       NN   NNN     C
C    KK    KK SSSSSS  AA    AA  UUUUUU  XX    XX FF       NN    NN     C
C                                                                      C
C -------------------------------------------------------------------- C
C  KSAUXFN GENERATES THE AUXILLIARY FUNCTION NEEDED FOR FERMIONIC      C
C  KÄLLÉN-SABRY POTENTIALS AND CHARGE DENSITIES.                       C
C  STRATEGY: MARCH BACKWARDS THROUGH ALL T-INTEGRATION ARGUMENTS,      C
C            KNOWING THE INTEGRAL AT INFINITY IS ZERO, AND CUT UP THE  C
C            T-VALUES FOR A FULL 11-POINT INTEGRAL SWEEP.              C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/AUXK/YFNC(0:NTVP)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
      DIMENSION TGRD(0:NTVP)
C
C     INTEGRATION GRID DETAILS
      TMIN = 1.0D0
      TMAX = 1.0D5
      HSPC = DLOG(TMAX/TMIN)/DFLOAT(NTVP)
C
C     GENERATE THE RIGHT GRID OF T VALUES AND DIFFERENTIAL ELEMENTS
      DO K=0,NTVP
        TGRD(K) = TMIN*DEXP(HSPC*DFLOAT(K))
      ENDDO
C
C     INITIALISE THE INTEGRAL AT PRACTICAL INFINITY
      YFNC(NTVP) = 0.0D0
C
C     LOOP BACKWARDS OVER THE T-VARIABLE INTEGRATION POINTS
      DO K=NTVP,2,-1
C
C       MIN AND MAX Y VALUES
        YMAX = TGRD(K  )
        YMIN = TGRD(K-1)
        YDIF = (YMAX-YMIN)/DFLOAT(10)
C
        YBIN = 0.0D0
C
C       LOOP OVER 10 UNIFORMLY-SPACED INTEGRATION POINTS
        DO J=0,10
C
C         VALUE OF Y ALONG INTEGRATION GRID AND DIFF. ELEMENT
          Y = YMIN + YDIF*DFLOAT(J)
C
C         INTEGRAND TERMS
          Y2 = Y*Y
          C1 = (3.0D0*Y2-1.0D0)/(Y*(Y2-1.0D0))
          C2 = DLOG(Y + DSQRT(Y2-1.0D0))
          C3 = 1.0D0/DSQRT(Y2-1.0D0)
          C4 = DLOG(8.0D0*Y*(Y2-1.0D0))
C
C         ADD TO THE COUNTER
          YBIN = YBIN + 5.0D0*EXTINT11(C1*C2-C3*C4,J,10)/2.99376D+5
C
        ENDDO
C
C       ADDITIONAL CONTRIBUTION TO THE FUNCTION
        YFNC(K-1) = YFNC(K) + YDIF*YBIN
C
      ENDDO
C
C     SPECIAL VALUE AT T=1
      YFNC(0) = PI*PI/4.0D0
C
      RETURN
      END
C
C
      FUNCTION CHIUEHF(X,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       CCCCCC  HH    HH IIII UU    UU EEEEEEEE HH    HH FFFFFFFF      C
C      CC    CC HH    HH  II  UU    UU EE       HH    HH FF            C
C      CC       HH    HH  II  UU    UU EE       HH    HH FF            C
C      CC       HHHHHHHH  II  UU    UU EEEEEE   HHHHHHHH FFFFFF        C
C      CC       HH    HH  II  UU    UU EE       HH    HH FF            C
C      CC    CC HH    HH  II  UU    UU EE       HH    HH FF            C
C       CCCCCC  HH    HH IIII  UUUUUU  EEEEEEEE HH    HH FF            C
C                                                                      C
C -------------------------------------------------------------------- C
C  CHIUEHF GENERATES THE AUXILLIARY FUNCTION NEEDED FOR FERMIONIC      C
C  UEHLING POTENTIALS AND CHARGE DENSITIES.                            C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
      TMIN = 1.0D0
      TMAX = 1.0D4
C
      HSPC = DLOG(TMAX/TMIN)/DFLOAT(NTVP)
C
C     LOOP OVER NTVP EXPONENTIALLY-SPACED INTEGRATION POINTS
      CHIUEHF = 0.0D0
      DO K=0,NTVP
C
        T    = TMIN*DEXP(HSPC*DFLOAT(K))
        TSPC = HSPC*T
C
C       INTEGRAND TERMS
        T2 = T*T
        C1 = 1.0D0/(T**N)
        C2 = DSQRT(1.0D0 - 1.0D0/T2)
        C3 = 1.0D0 + 0.5D0/T2
        C4 = DEXP(-X*T)
C
C       ADD TO THE COUNTER
        CHIUEHF = CHIUEHF + TSPC*EXTINT11(C1*C2*C3*C4,K,NTVP)
C
      ENDDO
C
      CHIUEHF = 5.0D0*CHIUEHF/2.99376D+5
      CHIUEHF = 2.0D0*CHIUEHF/3.0D0
C
      RETURN
      END
C
C
      FUNCTION CHIUEHB(X,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       CCCCCC  HH    HH IIII UU    UU EEEEEEEE HH    HH BBBBBBB       C
C      CC    CC HH    HH  II  UU    UU EE       HH    HH BB    BB      C
C      CC       HH    HH  II  UU    UU EE       HH    HH BB    BB      C
C      CC       HHHHHHHH  II  UU    UU EEEEEE   HHHHHHHH BBBBBBB       C
C      CC       HH    HH  II  UU    UU EE       HH    HH BB    BB      C
C      CC    CC HH    HH  II  UU    UU EE       HH    HH BB    BB      C
C       CCCCCC  HH    HH IIII  UUUUUU  EEEEEEEE HH    HH BBBBBBB       C
C                                                                      C
C -------------------------------------------------------------------- C
C  CHIUEHB GENERATES THE AUXILLIARY FUNCTION NEEDED FOR FERMIONIC      C
C  UEHLING POTENTIALS AND CHARGE DENSITIES.                            C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
      TMIN = 1.0D0
      TMAX = 1.0D4
C
      HSPC = DLOG(TMAX/TMIN)/DFLOAT(NTVP)
C
C     LOOP OVER NT EXPONENTIALLY-SPACED INTEGRATION POINTS
      CHIUEHB = 0.0D0
      DO K=0,NTVP
C
        T    = TMIN*DEXP(HSPC*DFLOAT(K))
        TSPC = HSPC*T
C
C       INTEGRAND TERMS
        T2 = T*T
        C1 = 0.5D0/(T**N)
        C2 = 1.0D0 - 1.0D0/T2
        C3 = DSQRT(C2)
        C4 = DEXP(-X*T)
C
C       ADD TO THE COUNTER
        CHIUEHB = CHIUEHB + TSPC*EXTINT11(C1*C2*C3*C4,K,NTVP)
C
      ENDDO
C
      CHIUEHB = 5.0D0*CHIUEHB/2.99376D+5
      CHIUEHB = 2.0D0*CHIUEHB/3.0D0
C
      RETURN
      END
C
C
      SUBROUTINE VWKRFPT(IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   VV    VV WW      WW KK    KK RRRRRRR  FFFFFFFF PPPPPPP TTTTTTTT    C
C   VV    VV WW      WW KK   KK  RR    RR FF       PP    PP   TT       C
C   VV    VV WW  WW  WW KK  KK   RR    RR FF       PP    PP   TT       C
C   VV    VV WW WWWW WW KKKKK    RR    RR FFFFFF   PP    PP   TT       C
C    VV  VV  WWWW  WWWW KK  KK   RRRRRRR  FF       PPPPPPP    TT       C
C     VVVV   WWW    WWW KK   KK  RR    RR FF       PP         TT       C
C      VV    WW      WW KK    KK RR    RR FF       PP         TT       C
C                                                                      C
C -------------------------------------------------------------------- C
C  VWKRFPT GENERATES THE WICHMANN-KROLL POTENTIAL AS IF FROM A POINT.  C
C**********************************************************************C
      PARAMETER(LIMIT=1000000,LENW=5*LIMIT)
C
      INCLUDE 'parameters.h'
C
      EXTERNAL VKERNEL
C
      LOGICAL FILETHERE
C
      DIMENSION IWORK(LIMIT)
      DIMENSION WORK(LENW)
C
      COMMON/TVAL/T
      COMMON/RVAL/R
      COMMON/QWKR/RAD(0:NRAD),VVAC(MCT,0:NRAD),RORI,RMID,RMAX,NLIN,NEXP
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
C
C     ATTEMPT TO READ STRUCTURELESS POTENTIAL FROM EXTERNAL FILE
      INQUIRE(FILE='spectral/VR2_wkr_electron.dat',EXIST=FILETHERE)
      IF(FILETHERE) THEN
        OPEN(UNIT=8,FILE='spectral/VR2_wkr_electron.dat',
     &                                               STATUS='UNKNOWN')
        REWIND(UNIT=8)
        DO N=0,NRAD
          READ(8, *) M,R,VVAC(IZ,N)
C         IF(M.NE.N) THEN
          IF(M.NE.N.OR.R.NE.RAD(N)) THEN
          	FILETHERE = .FALSE.
            GOTO 111
          ENDIF
        ENDDO
111     CONTINUE
        CLOSE(UNIT=8)
C       POTENTIAL SUCCESSFULLY READ
        RETURN
      ENDIF
      IF(.NOT.FILETHERE) THEN
        WRITE(6, *) 'In VWKRFPT: problem reading structureless VWK(r).'
        WRITE(7, *) 'In VWKRFPT: problem reading structureless VWK(r).'
        DO N=0,NRAD
          VVAC(IZ,N) = 0.0D0
        ENDDO
      ENDIF
110   CONTINUE
C
C     MULTIPLIER FOR COMPTON WAVELENGTH IN FM (FOR MUON OR TAUON FIELD)
      CMPF = CMPW/EMSS
C
C     PERTURBATIVE ORDERS (ALPHA/PI) AND (Z*ALPHA) -- IGNORE Z
      ALPI = 1.0D0/(CV*PI)
      ZALF = 1.0D0/CV
C
C     LOOP OVER THE WHOLE RADIAL GRID
      DO N=0,NRAD
C
C       MODIFIED RADIUS, RAD/CMPF
        R = RAD(N)/CMPF
C
C       INITIALISE RAW POTENTIAL (IT'S ACTUALLY V(R)*R)
        VR = 0.0D0
C
C        POTENTIAL NEAR THE ORIGIN
         IF(RAD(N).LE.0.001D0) THEN
C        CASE 1: USE THE BLOMQVIST FORMULA FOR SMALL R
C
          ICASE = 1
          VR = VWKRPWR(R)
C
C       NUMERICAL INTEGRATION REQUIRED
        ELSEIF(RAD(N).GT.0.001D0.AND.RAD(N).LT.10.0D0) THEN
C
C         FIRST INTEGRAL ON [0:1]
          CALL DQAGS(VKERNEL,0.0D0,1.0D0,1.0D-12,1.0D-12,RINT1,ABSERR,
     &                             NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
C
C         SECOND INTEGRAL ON [1:∞)
          CALL DQAGS(VKERNEL,1.0D0,1.0D3,1.0D-12,1.0D-12,RINT2,ABSERR,
     &                             NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
C
          ICASE = 2
          VR = RINT1+RINT2
C
C       ASYMPTOTIC FORMULA
        ELSE
C
          ICASE = 3
          VR = VWKRASY(R)
C
        ENDIF
C
C       APPLY PERTURBATIVE ORDERS AND RADIAL WEIGHTING FACTOR
        VVAC(IZ,N) =-CV*ALPI*ZALF*ZALF*ZALF*RAD(N)*VR
C
      ENDDO
C
C     WRITE TO FILE
      OPEN(UNIT=8,FILE='spectral/VR2_wkr_electron.dat',STATUS='NEW')
      DO N=0,NRAD
        WRITE(8, *) N,RAD(N),VVAC(IZ,N)
      ENDDO
      CLOSE(UNIT=8)
C
      RETURN
      END
C
C
      FUNCTION VWKRPWR(R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C VV    VV WW       WW KK    KK RRRRRRR  PPPPPPP  WW       WW RRRRRRR  C
C VV    VV WW   W   WW KK   KK  RR    RR PP    PP WW   W   WW RR    RR C
C VV    VV WW  WWW  WW KK  KK   RR    RR PP    PP WW  WWW  WW RR    RR C
C VV    VV WW WW WW WW KKKKK    RR    RR PP    PP WW WW WW WW RR    RR C
C  VV  VV  WWWW   WWWW KK  KK   RRRRRRR  PPPPPPP  WWWW   WWWW RRRRRRR  C
C   VVVV   WWW     WWW KK   KK  RR    RR PP       WWW     WWW RR    RR C
C    VV    WW       WW KK    KK RR    RR PP       WW       WW RR    RR C
C                                                                      C
C -------------------------------------------------------------------- C
C  VWKRPWR EVALUATES THE TERMS {...} IN THE POWER SERIES EXPANSION FOR C
C  THE POINT-NUCLEAR WICHMANN-KROLL POTENTIAL, AS GIVEN BY BLOMQVIST   C
C  NUCLEAR PHYSICS B45 (95-103) IN 1972. A FACTOR 1/R HAS BEEN REMOVED.C
C**********************************************************************C
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
      PI2 = PI*PI
      PI3 = PI2*PI
      PI4 = PI2*PI2
C
C     ARGUMENT MUST BE NON-NEGATIVE
      IF(R.LT.0.0D0) THEN
        VWKRPWR = 0.0D0
        RETURN
      ENDIF
C
C     SPECIAL VALUE AT THE ORIGIN      
      T0 =-2.0D0*ZTA3/3.0D0 + PI2/6.0D0 - 7.0D0/9.0D0
      IF(R.LE.0.0D0) THEN
        VWKRPWR = T0
        RETURN
      ENDIF
C
C     VALUE FOR R>0
      R2 = R*R
      RLOGG = DLOG(R)+EULR
C
      T1 = 2.0D0*PI*ZTA3 - PI3/4.0D0
      T2 =-6.0D0*ZTA3 + PI4/16.D0 + PI2/6.0D0
      T3 = 2.0D0*PI/9.0D0
      T4 = 2.0D0*PI*ZTA3/3.0D0 + 4.0D0*PI*TWLG/9.0D0 - 31.0D0*PI/27.0D0
      T5 = 1.0D0/12.0D0
      T6 = 5.0D0*PI2/54.0D0 - 19.0D0/36.0D0
      T7 = 13.0D0*ZTA3/18.0D0 - 109.0D0*PI2/432.0D0 + 859.0D0/864.0D0
C
      VWKRPWR = T0 + T1*R + T2*R2 + T3*RLOGG*R*R2 + T4*R*R2 +
     &          T5*RLOGG*RLOGG*R2*R2 + T6*RLOGG*R2*R2 + T7*R2*R2
C
      RETURN
      END
C
C
      FUNCTION VWKRASY(R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  VV    VV WW       WW KK    KK RRRRRRR     AA     SSSSSS  YY    YY   C
C  VV    VV WW   W   WW KK   KK  RR    RR   AAAA   SS    SS YY    YY   C 
C  VV    VV WW  WWW  WW KK  KK   RR    RR  AA  AA  SS        YY  YY    C 
C  VV    VV WW WW WW WW KKKKK    RR    RR AA    AA  SSSSSS    YYYY     C 
C   VV  VV  WWWW   WWWW KK  KK   RRRRRRR  AAAAAAAA       SS    YY      C 
C    VVVV   WWW     WWW KK   KK  RR    RR AA    AA SS    SS    YY      C 
C     VV    WW       WW KK    KK RR    RR AA    AA  SSSSSS     YY      C
C                                                                      C
C -------------------------------------------------------------------- C
C  VWKRASY IS AN ASYMPTOTIC EXPANSION OF VWK, VALID FOR R>10.          C
C  THIS IS A QUINEY-ORIGINAL POWER SERIES EXPANSION IN 1/R^2.          C
C  OVERALL CONSTANTS AND A FACTOR 1/R HAVE BEEN REMOVED.               C
C**********************************************************************C
C
C     ARGUMENT MUST BE NON-NEGATIVE
      IF(R.LT.0.0D0) THEN
        VWKRASY = 0.0D0
        RETURN
      ENDIF
C
C     LARGE VALUES OF R
      R2  = 1.0D0/(R**2)
      R4  = R2*R2
      R6  = R2*R4
      R8  = R6*R2
      R10 = R8*R2
      R12 = R10*R2
      R14 = R12*R2
C
      VWKRASY = (2.0D0/225.0D0)*R4
     &        + (59.0D0/1323.0D0)*R6
     &        + (659.0D0/1575.0D0)*R8
     &        + (8.3912D0/1.2705D0)*R10
     &        + (217.824448D0/1.366365D0)*R12
     &        + (27479.392D0/5.005D0)*R14
C
      RETURN
      END
C
C
      FUNCTION VKERNEL(T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    VV    VV KK    KK EEEEEEEE RRRRRRR  NN    NN EEEEEEEE LL          C
C    VV    VV KK   KK  EE       RR    RR NNN   NN EE       LL          C
C    VV    VV KK  KK   EE       RR    RR NNNN  NN EE       LL          C
C    VV    VV KKKKK    EEEEEE   RR    RR NN NN NN EEEEEE   LL          C
C     VV  VV  KK  KK   EE       RRRRRRR  NN  NNNN EE       LL          C
C      VVVV   KK   KK  EE       RR    RR NN   NNN EE       LL          C
C       VV    KK    KK EEEEEEEE RR    RR NN    NN EEEEEEEE LLLLLLLL    C
C                                                                      C
C -------------------------------------------------------------------- C
C  VKERNEL APPLIES THE EXPONENTIAL FACTOR TO GG(T).                    C
C  GG(T) IS THE FUNCTION SPECIFICALLY FOR THE WICHMANN-KROLL DIAGRAM.  C
C**********************************************************************C
      COMMON/RVAL/R
C
      VKERNEL = DEXP(-2.0D0*R*T)*GG(T,1)
C
      RETURN
      END
C
C
      FUNCTION CHIWKRF(X,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      CCCCCC  HH    HH IIII WW      WW KK    KK RRRRRRR  FFFFFFFF     C
C     CC    CC HH    HH  II  WW      WW KK   KK  RR    RR FF           C
C     CC       HH    HH  II  WW  WW  WW KK  KK   RR    RR FF           C
C     CC       HHHHHHHH  II  WW WWWW WW KKKKK    RR    RR FFFFFF       C
C     CC       HH    HH  II  WWWW  WWWW KK  KK   RRRRRRR  FF           C
C     CC    CC HH    HH  II  WWW    WWW KK   KK  RR    RR FF           C
C      CCCCCC  HH    HH IIII WW      WW KK    KK RR    RR FF           C
C                                                                      C
C -------------------------------------------------------------------- C
C  CHIKSBF GENERATES THE AUXILLIARY FUNCTION NEEDED FOR FERMIONIC      C
C  WICHMANN-KROLL POTENTIALS AND CHARGE DENSITIES.                     C
C**********************************************************************C
      PARAMETER(LIMIT=1000000,LENW=5*LIMIT)
C
      INCLUDE 'parameters.h'
C
      EXTERNAL CKERNEL
C
      DIMENSION IWORK(LIMIT)
      DIMENSION WORK(LENW)
C
      COMMON/NVAL/NN
      COMMON/XVAL/XX
C
      NN = N
      XX = X
C
C     FIRST INTEGRAL ON [0:1]
      CALL DQAGS(CKERNEL,0.0D0,1.0D0,1.0D-12,1.0D-12,RINT1,ABSERR,NEVAL,
     &                                   IER,LIMIT,LENW,LAST,IWORK,WORK)
C
C     SECOND INTEGRAL ON [1:∞)
      CALL DQAGS(CKERNEL,1.0D0,1.0D3,1.0D-12,1.0D-12,RINT2,ABSERR,NEVAL,
     &                                   IER,LIMIT,LENW,LAST,IWORK,WORK)
C
      CHIWKRF = RINT1+RINT2
C
      RETURN
      END
C
C
      FUNCTION CKERNEL(T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     CCCCCC  KK    KK EEEEEEEE RRRRRRR  NN    NN EEEEEEEE LL          C
C    CC    CC KK   KK  EE       RR    RR NNN   NN EE       LL          C
C    CC       KK  KK   EE       RR    RR NNNN  NN EE       LL          C
C    CC       KKKKK    EEEEEE   RR    RR NN NN NN EEEEEE   LL          C
C    CC       KK  KK   EE       RRRRRRR  NN  NNNN EE       LL          C
C    CC    CC KK   KK  EE       RR    RR NN   NNN EE       LL          C
C     CCCCCC  KK    KK EEEEEEEE RR    RR NN    NN EEEEEEEE LLLLLLLL    C
C                                                                      C
C -------------------------------------------------------------------- C
C  RKERNEL SELECTS OUT A PARTICULAR GG FAMILY.                         C
C**********************************************************************C
      COMMON/NVAL/NN
      COMMON/XVAL/XX
C
      RKERNEL = GG(T,1-NN)*DEXP(-XX*T)
C
      RETURN
      END
C
C
      FUNCTION GG(TT,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                           GGGGGG   GGGGGG                            C
C                          GG    GG GG    GG                           C
C                          GG       GG                                 C
C                          GG       GG                                 C
C                          GG   GGG GG   GGG                           C
C                          GG    GG GG    GG                           C
C                           GGGGGG   GGGGGG                            C
C                                                                      C
C -------------------------------------------------------------------- C
C  GG(T) IS THE FUNCTION SPECIFICALLY FOR THE WICHMANN-KROLL DIAGRAM.  C
C  OVERALL CONSTANTS AND A FACTOR 1/R HAVE BEEN REMOVED.               C
C  USE K=1 FOR THE STRUCTURELESS POTENTIAL AND OTHER VALUES FOR CHI(K).C
C**********************************************************************C
      PARAMETER(LIMIT=1000000,LENW=5*LIMIT)
C
      EXTERNAL GKERNEL
C
      DIMENSION IWORK(LIMIT)
      DIMENSION WORK(LENW)
C
      COMMON/TVAL/T
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
      T  = TT
      T2 = T*T
C
      IF(T.LE.0.0D0) THEN 
C     CASE 0: RETURN G=0.0D0 FOR T=<0.0D0
C
        G = 0.0D0
C
      ELSEIF(T.GT.0.0D0.AND.T.LE.0.1D0) THEN      
C     CASE 1: EXPANSION METHOD FOR SMALL T (T<0.1D0)
C
C       LOWEST-ORDER POLYNOMIAL IN THE EXPANSION
        TK = T**(4-K)
C
        G  = (1.6D0/67.5D0)*TK
     &     + (4.72D0/198.45D0)*TK*T2
     &     + (1.0544D0/49.6125D0)*TK*T2*T2
     &     + (6.71296D0/360.18675D0)*TK*T2*T2*T2
     &     + (3.485191168D0/213.050462625D0)*TK*T2*T2*T2*T2
     &     + (4.39670272D0/304.35780375D0)*TK*T2*T2*T2*T2*T2
C
      ELSEIF(T.GT.0.1D0.AND.T.LE.(1.0D0-1.0D-12)) THEN 
C     CASE 2: DIRECT NUMERICAL INTEGRATION FOR 0.1D0 < T < 1.0D0
C
        A = 0.0D0
        B = DASIN(1.0D0)
        EPSABS = 1.0D-12
        EPSREL = 1.0D-12
        CALL DQAGS(GKERNEL,A,B,EPSABS,EPSREL,RINT1,ABSERR,NEVAL,IER,
     &                                      LIMIT,LENW,LAST,IWORK,WORK)
        G = RINT1/(T**(K-1))
C
      ELSEIF(DABS(T-1.0D0).LE.1.0D-12) THEN
C     CASE 3: VALUE VERY CLOSE TO T=1
C
        G = 0.45985669359783D0/(T**(K-1))
C
      ELSEIF(T.GT.(1.0D0+1.0D-12)) THEN
C     CASE 4: DIRECT NUMERICAL INTEGRATION FOR T>1 (CORRECTION NEAR X=1)
C
        DELTA1 = DMIN1(1.0D-3*(T-1.0D0),1.0D-3)
C
C       RINT1
        A = 0.0D0
        B = DASIN((1.0D0-DELTA1)/T)
        EPSABS = 1.0D-10
        EPSREL = 1.0D-10
        CALL DQAGS(GKERNEL,A,B,EPSABS,EPSREL,RINT1,ABSERR1,NEVAL,IER,
     &                                      LIMIT,LENW,LAST,IWORK,WORK)
C
C       RINT2
        A = DASIN((1.0D0+DELTA1)/T)
        B = DASIN(1.0D0)
        EPSABS = 1.0D-10
        EPSREL = 1.0D-10
        CALL DQAGS(GKERNEL,A,B,EPSABS,EPSREL,RINT2,ABSERR2,NEVAL,IER,
     &                                      LIMIT,LENW,LAST,IWORK,WORK)
C
C       RINT3
        RINT3 = BRIDGE(T,DELTA1)
C
C       RINT4: THERE IS A CONSTANT TO BE ADDED FOR T>1
        PI2   = PI*PI
        T4    = T**4
        RINT4 =-PI2*DSQRT(T*T-1.0D0)/(12.0D0*T4)
C
C       FINAL SUM
        G = (RINT1+RINT2+RINT3+RINT4)/(T**(K-1))
C
      ENDIF
C
      GG = G
C
      RETURN
      END
C
C
      FUNCTION GKERNEL(THETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     GGGGGG  KK    KK EEEEEEEE RRRRRRR  NN    NN EEEEEEEE LL          C
C    GG    GG KK   KK  EE       RR    RR NNN   NN EE       LL          C
C    GG       KK  KK   EE       RR    RR NNNN  NN EE       LL          C
C    GG       KKKKK    EEEEEE   RR    RR NN NN NN EEEEEE   LL          C
C    GG   GGG KK  KK   EE       RRRRRRR  NN  NNNN EE       LL          C
C    GG    GG KK   KK  EE       RR    RR NN   NNN EE       LL          C
C     GGGGGG  KK    KK EEEEEEEE RR    RR NN    NN EEEEEEEE LLLLLLLL    C
C                                                                      C
C -------------------------------------------------------------------- C
C  GKERNEL IS THE KERENEL OF THE INTEGRAL GG(T), AS IN BLOMQVIST (1972)C
C  THE VARIABLE X IS DEFINED BY X = T*SIN(THETA) SO DX = T*COS(THETA)  C
C  AND D(THETA) = SQRT(T^2-X^2). A FACTOR OF 1/T^4 IS ABSORBED AT      C
C  THIS LEVEL FOR STABILITY OF THE INTEGRATION FOR LARGE T.            C
C**********************************************************************C
      COMMON/TVAL/T
C
      T4 = T**4
      X  = T*DSIN(THETA)
      DX = T*DCOS(THETA)
      GKERNEL = DX*DX*FWK(X)/T4
C
      RETURN
      END
C
C
      FUNCTION FWK(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                    FFFFFFFF WW       WW KK    KK                     C
C                    FF       WW   W   WW KK   KK                      C
C                    FF       WW  WWW  WW KK  KK                       C
C                    FFFFFF   WW WW WW WW KKKKK                        C
C                    FF       WWWW   WWWW KK  KK                       C
C                    FF       WWW     WWW KK   KK                      C
C                    FF       WW       WW KK    KK                     C
C                                                                      C
C -------------------------------------------------------------------- C
C FWK EVALUATES THE FUNCTION F(X) DEFINED BY BLOMQVIST (1972).         C
C**********************************************************************C
C
      IF(X.LT.0.0D0) THEN
C     CASE 1: REJECT NEGATIVE VALUES OF X AND SET FUNCTION EQUAL TO ZERO
C
        WRITE(6,*) 'ILLEGAL VALUE OF X IN FWK(X)'
        FWK = 0.0D0
        RETURN
C
      ELSEIF(X.GE.0.0D0.AND.X.LE.0.1D0) THEN
C     CASE 2: EMPLOY SERIES EXPANSION OF FWK FOR SMALL X, 0<X<X0
C
        X2 = X*X
        X4 = X2*X2
        X5 = X4*X
        X6 = X4*X2
C
        FWK = X5*(14.0D0/45.0D0 + X2*(59.0D0/126.0D0)
     &      + X4*(7249.0D0/12600.0D0)+X6*(136357.0D0/207900.0D0))
C
        RETURN
C
      ELSEIF(X.GT.0.1D0.AND.X.LT.1.0D0) THEN
C     CASE 3: USE COMPLETE F(X) EXPRESSION FOR X0<X<1 
C
        X2 = X*X
        T1 = DLOG(1.0D0-X2)
        T2 = DLOG((1.0D0+X)/(1.0D0-X))
C
        FWK =-2.0D0*X*POLYLOG2(X2) - X*T1*T1 + (1.0D0-X2)*T1*T2/X2
     &      + (1.0D0-X2)*T2*T2/(4.0D0*X)
     &      + (2.0D0-X2)*T1/(X*(1.0D0-X2))
     &      + (3.0D0-2.0D0*X2)*T2/(1.0D0-X2) - 3.0D0*X
C
        RETURN
C
      ELSEIF(X.EQ.1.0D0) THEN
C     CASE 4: F(X) IS NOT DEFINED AT X=1
C             SET F(1)=0 AND HANDLE SEPARATELY
C
        FWK = 0.0D0
C
      ELSEIF(X.GT.1.0D0) THEN
C     CASE 5: USE COMPLETE F(X) EXPRESSION FOR X>1
C
        X2 = X*X
        A0 = 1.0D0/X
        A1 = 1.0D0/X2
        A2 = 1.0D0-A1
        T0 = DLOG(X)
        T1 = DLOG(A2)
        T2 = DLOG((X+1)/(X-1))
C
        FWK = A1*POLYLOG2(A1) 
     &      - ((3.0D0*X2+1.0D0)/(2.0D0*X))*(POLYLOG2(A0)-POLYLOG2(-A0))
     &      - ((2.0D0*X2-1.0D0)/(2.0D0*X2))*(T1*T1 + T2*T2)
     &      - (2.0D0*X-1.0D0)*T1*T2
     &      + ((3.0D0*X2+1.0D0)/(4.0D0*X))*T2*T2
     &      - 2.0D0*T0*T1
     &      - ((3.0D0*X2+1.0D0)/(2.0D0*X))*T0*T2
     &      + (5.0D0-(X*(3.0D0*X2-2.0D0)/(X2-1.0D0)))*T1
     &      + (((3.0D0*X2+2.0D0)/X)-((3.0D0*X2-2.0D0)/(X2-1.0D0)))*T2
     &      + 3.0D0*T0-3.0D0
C
        RETURN
C
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION BRIDGE(T,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          BBBBBBB  RRRRRRR  IIII DDDDDDD   GGGGGG  EEEEEEEE           C
C          BB    BB RR    RR  II  DD    DD GG    GG EE                 C
C          BB    BB RR    RR  II  DD    DD GG       EE                 C
C          BBBBBBB  RR    RR  II  DD    DD GG       EEEEEE             C
C          BB    BB RRRRRRR   II  DD    DD GG   GGG EE                 C
C          BB    BB RR    RR  II  DD    DD GG    GG EE                 C
C          BBBBBBB  RR    RR IIII DDDDDDD   GGGGGG  EEEEEEEE           C
C                                                                      C
C -------------------------------------------------------------------- C
C  BRIDGE ANALYTICALLY HANDLES THE DISCONTINUOUS PART OF THE WK KERNEL C
C  FOR ARGUMENTS (1+DELTA), USING TERM-BY-TERM CANCELLATIONS OF TERMS. C
C**********************************************************************C
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
C
      DLG = DLOG(D)
C
      TR = T*T-1.0D0
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
C
C     COEFFICIENT OF D^1
      CTD1 = (1.0D0/3.0D0)*(TR**3.0D0)*(2.0D0*PI*PI*TR
     &     + 3.0D0*(-5.0D0-2.0D0*TWLG*TWLG+T2*(5.0D0-12.0D0*TWLG
     &     + 2.0D0*TWLG*TWLG)+10.0D0*TWLG)
     &     + 3.0D0*TR*(1.0D0+4.0D0*TWLG)*DLG)
      CTD1 =-CTD1/T4
C
C     COEFFICIENT OF D^2
      CTD2 = (1.0D0/8.0D0)*(TR**3.0D0)*(44.0D0-36.0D0*T2+PI*PI*TR
     &     + 8.0D0*(3.0D0*T2-4.0D0)*DLG+8.0D0*DLG*DLG)
      CTD2 =-CTD2/T4
C
C     COEFFICIENT OF D^3
      CTD3 = (1.0D0/216.0D0)*TR*(3.0D0*PI*PI*(13.0D0
     &     - 25.0D0*T2+3.0D0*T4+9.0D0*T6)
     &     + 6.0D0*(-65.0D0-54.0D0*TWLG+6.0D0*TWLG*TWLG)
     &     + 4.0D0*T2*(284.0D0+297.0D0*TWLG+81.0D0*TWLG*TWLG)
     &     + 4.0D0*T6*(140.0D0+33.0D0*TWLG+117.0D0*TWLG*TWLG)
     &     - 2.0D0*T4*(653.0D0+462.0D0*TWLG+414.0D0*TWLG*TWLG)
     &     - 3.0D0*TR*(174.0D0+72.0D0*TWLG+2.0D0*T4*(37.0D0+12.0D0*TWLG)
     &     - 2.0D0*T2*(130.0D0+72.0D0*TWLG))*DLG
     &     + 36.0D0*(TR**2.0D0)*(T2-3.0D0)*DLG*DLG)
      CTD3 = CTD3/T4
C
C     COEFFICIENT OF D^4
      CTD4 = (1.0D0/2.88D2)*TR*(4.46D2-2.97D2*T2-9.24D2*T4
     &     + 8.2D2*T6+3.0D0*PI*PI*(-4.0D0+2.4D1*T2-3.3D1*T4+1.3D1*T6)
     &     - 1.2D1*(1.4D1+2.1D1*T2-9.6D1*T4+6.4D1*T6)*DLG
     &     + 7.2D1*(-2.0D0+1.3D1*T2-1.8D1*T4+8.0D0*T6)*DLG*DLG)
      CTD4 =-CTD4/T4
C
      PREFAC = TR**(-7.0D0/2.0D0)
      BRIDGE = PREFAC*(CTD1*D + CTD2*D*D + CTD3*D*D*D + CTD4*D*D*D*D)
C
      RETURN
      END
C
C
      FUNCTION POLYLOG2(X)
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
C POLYLOG2 EVALUATES LI_2(X) for |X| =< 1.                             C
C**********************************************************************C
      DATA BASEL,EPS/1.644934066848226D0,1.0D-15/
C
C     1. REJECT VALUES OUTSIDE |X|=< 1
C
      IF(DABS(X).GT.1.0D0) THEN
C
        WRITE(6,*) 'ILLEGAL VALUE OF X in POLYLOG2. X= ',X
        POLYLOG2=0.0D0
        RETURN
C
C     2. THE CASE X=1 IS THE BASEL PROBLEM, PI^2/6
C
      ELSEIF(X.EQ.1.0D0) THEN
C
        POLYLOG2=BASEL
        RETURN
C
C     3. THE CASE X=0 IS TRIVIAL
C
      ELSEIF(X.EQ.0.0D0) THEN
C
        POLYLOG2=0.0D0
        RETURN
C
C     4. THE CASE X=-1 IS -PI^2/12
C
      ELSEIF(X.EQ.-1.0D0) THEN
C
        POLYLOG2=-BASEL/2.0D0
C
C     5. EVALUATE USING THE SERIES REPRESENTATION FOR X< 0.5. THE
C        EXIT ACKNOWLEDGES THAT |X| COULD BE CLOSE TO ZERO
C
      ELSEIF(X.LE.5.0D-1.AND.X.GE.-5.0D-1) THEN
C
        SUM=0.0D0
        XI=X
        DO I=1,1000000000
          RI=DBLE(I)
          TERM=XI/(RI*RI)
          SUM=SUM+TERM
          XI=XI*X
          IF(DABS(TERM/SUM).LE.EPS.OR.DABS(TERM).LE.EPS) GOTO 10
        ENDDO
10      CONTINUE
        POLYLOG2=SUM
        RETURN
C
C     6. EVALUATE FOR -1<X<-0.5 USING THE REFLECTION FORMULA  
C
      ELSEIF(X.GT.-1.0D0.AND.X.LT.-5.0D-1) THEN
C
        XM=-X
        SUM1=0.0D0
        X0=1.0D0-XM
        XI=X0
        DO I=1,1000000000
          RI=DBLE(I)
          TERM=XI/(RI*RI)
          SUM1=SUM1+TERM
          XI=XI*X0
          IF(DABS(TERM/SUM1).LE.EPS.OR.DABS(TERM).LE.EPS) GOTO 11
        ENDDO
11      CONTINUE     
        SUM2=0.0D0
        XP=1.0D0-(XM*XM)
        XI=XP
        DO I=1,1000000000
          RI=DBLE(I)
          TERM=XI/(RI*RI)
          SUM2=SUM2+TERM
          XI=XI*XP
          IF(DABS(TERM/SUM2).LE.EPS.OR.DABS(TERM).LE.EPS) GOTO 12
        ENDDO
12      CONTINUE
c       write(6,*) X,-BASEL/2.0D0, -DLOG(XM)*DLOG(1.0D0+XM),
c     & SUM1,-SUM2/2.0D0
        POLYLOG2=-(BASEL/2.0D0)-DLOG(XM)*DLOG(1.0D0+XM)+SUM1
     &           -(SUM2/2.0D0)      
        RETURN   
C
C     7. REMAINING RANGE IS 0.5 < X < 1. USE THE EULER REFLECTION FORMULA
C
      ELSE
C
        X0=1.0D0-X
        SUM=0.0D0
        XI=X0
        DO I=1,100000000
          RI=DBLE(I)
          TERM=XI/(RI*RI)
          SUM=SUM+TERM
          XI=XI*X0
          IF(DABS(TERM/SUM).LE.EPS.OR.DABS(TERM).LE.EPS) GOTO 20
        ENDDO
20      CONTINUE
        POLYLOG2=BASEL-DLOG(X0)*DLOG(X)-SUM
        RETURN
C
      ENDIF
C
      RETURN
      END
C
C
      FUNCTION CHIKSBF(X,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       CCCCCC  HH    HH IIII KK    KK  SSSSSS  BBBBBBB  FFFFFFFF      C
C      CC    CC HH    HH  II  KK   KK  SS    SS BB    BB FF            C
C      CC       HH    HH  II  KK  KK   SS       BB    BB FF            C
C      CC       HHHHHHHH  II  KKKKK     SSSSSS  BBBBBBB  FFFFFF        C
C      CC       HH    HH  II  KK  KK         SS BB    BB FF            C
C      CC    CC HH    HH  II  KK   KK  SS    SS BB    BB FF            C
C       CCCCCC  HH    HH IIII KK    KK  SSSSSS  BBBBBBB  FF            C
C                                                                      C
C -------------------------------------------------------------------- C
C  CHIKSBF GENERATES THE AUXILLIARY FUNCTION NEEDED FOR FERMIONIC      C
C  KÄLLÉN-SABRY POTENTIALS AND CHARGE DENSITIES.                       C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/AUXK/YFNC(0:NTVP)
C
      TMIN = 1.0D0
      TMAX = 1.0D5
C
      HSPC = DLOG(TMAX/TMIN)/DFLOAT(NTVP)
C
C     LOOP OVER NTVP EXPONENTIALLY-SPACED INTEGRATION POINTS
      CHIKSBF = 0.0D0
      DO K=0,NTVP
C
        T    = TMIN*DEXP(HSPC*DFLOAT(K))
        TSPC = HSPC*T
C
C       INTEGRAND TERMS
        T2 = T*T
        T4 = T2*T2
        T6 = T4*T2
        T8 = T6*T2
        Z0 = T**(N-1)
        A1 = 13.0D0/(54.0D0*T2) + 7.0D0/(108.0D0*T4) + 2.0D0/(9.0D0*T6)
        A2 = DSQRT(T2-1.0D0)
        B1 =-44.0D0/(9.0D0*T) + 2.0D0/(3.0D0*T*T2)
     &       + 5.0D0/(4.0D0*T*T4) + 2.0D0/(9.0D0*T*T6)
        B2 = DLOG(T+DSQRT(T2-1.0D0))
        C1 = 4.0D0/(3.0D0*T2) + 2.0D0/(3.0D0*T4)
        C2 = DSQRT(T2-1.0D0)
        C3 = DLOG(8.0D0*T*(T2-1.0D0))
        D1 =-8.0D0/(3.0D0*T) + 2.0D0/(3.0D0*T*T4)
        D2 = YFNC(K)
        Z5 = DEXP(-X*T)
        
        Z1 = A1*A2
        Z2 = B1*B2
        Z3 = C1*C2*C3
        IF(K.EQ.0) Z3 = 0.0D0
        Z4 = D1*D2
C
C       ADD TO THE COUNTER
        CHIKSBF = CHIKSBF + TSPC*EXTINT11((Z1+Z2+Z3+Z4)*Z5/Z0,K,NTVP)
C
      ENDDO
C
      CHIKSBF =-5.0D0*CHIKSBF/2.99376D+5
C
      RETURN
      END
C
C
      SUBROUTINE VPOLPOT(FUNC,CMPF,RFM,FPTNT,AMP,RNUC,RLONG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    VV    VV PPPPPPP   OOOOOO  LL       PPPPPPP   OOOOOO  TTTTTTTT    C
C    VV    VV PP    PP OO    OO LL       PP    PP OO    OO    TT       C
C    VV    VV PP    PP OO    OO LL       PP    PP OO    OO    TT       C
C    VV    VV PP    PP OO    OO LL       PP    PP OO    OO    TT       C
C     VV  VV  PPPPPPP  OO    OO LL       PPPPPPP  OO    OO    TT       C
C      VVVV   PP       OO    OO LL       PP       OO    OO    TT       C
C       VV    PP        OOOOOO  LLLLLLLL PP        OOOOOO     TT       C
C                                                                      C
C -------------------------------------------------------------------- C
C  VPOLPOT CALCULATES THE POLARISED POTENTIAL BY THE FULLERTON         C
C  AND RINKER (1975) MULTIPOLE EXPANSION METHOD. THIS METHOD CAN BE    C
C  EXTENDED TO APPLY TO FERMIONIC/BOSONIC INTERACTIONS AS WELL AS      C
C  HIGHER-ORDER INTERACTIONS.                                          C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RFM(0:NRAD),FPTNT(0:NRAD)
C
      EXTERNAL FUNC
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR,ZTA3
      COMMON/CONV/CHZ,CEV,CCM,CFM,CNG,CDB
C
C     UNITLESS NUCLEAR RATIO
      RAT = RNUC/CMPF
C
C     NUCLEAR CHARGE MOMENT DETAILS (USE A GAUSSIAN)
      R0 = 1.0D0
      R2 = 2.0D0*RAT*RAT/3.0D0
      R4 = 2.0D0*RAT*RAT*RAT*RAT/9.0D0
C
C     IF R IS LARGE ENOUGH, EVALUATE THE POTENTIAL
      DO N=0,NRAD
        R = RFM(N)
        IF(R.LT.RLONG) GOTO 10
        X   = 2.0D0*R/CMPF
        SER = FUNC(X,1) + R2*FUNC(X,-1) + R4*FUNC(X,-3)
C       SER = FUNC(X,1)
        FPTNT(N) = AMP*SER/R
        IF(DABS(FPTNT(N)).LT.1.0D-16) GOTO 50
10      CONTINUE
      ENDDO
50    CONTINUE
C
      RETURN
      END
C
C
