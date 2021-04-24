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
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26)
C
      CHARACTER*1 DUMLIN
      CHARACTER*2 ELMNT(120)
      CHARACTER*4 HMLTN
      CHARACTER*7 HMINT(10),PTYPE(10)
      CHARACTER*40 FILNAM,STRING
C
      COMPLEX*16 C(MDM,MDM)
C
      COMMON/ATOM/ELMNT
      COMMON/COEF/C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/FILL/NCNF(MCT,MKP,MKP+1),NLVL(MCT,MKP),IFILL(MCT)
      COMMON/FLNM/STRING,FILNAM,LN,LF
      COMMON/LSHF/SHLV1,SHLV2,SHLV3,SHLV
      COMMON/PLOT/NPTYPE,PTYPE
      COMMON/PRMS/CV,HMLTN,ITER,IIL,ITREE,IEQS
      COMMON/PT1B/NHMINT,HMINT
      COMMON/SHLL/ALPHA,BETA,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR,IOCCM0
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
C     OUTPUT STRING
      READ(5,*) DUMLIN
      READ(5,*) FILNAM
C
C     CHOICE OF HAMILTONIAN HMLTN: NORL, BARE, DHFR, DHFP, OR DHFB
9     FORMAT(A4)
      READ(5, *) DUMLIN
      READ(5, 9) HMLTN
C
C     MAKE SURE THE HAMILTONIAN INPUT IS VALID
      IF(HMLTN.NE.'NORL'.AND.HMLTN.NE.'DHFR'.AND.HMLTN.NE.'DHFP'
     &                  .AND.HMLTN.NE.'DHFB'.AND.HMLTN.NE.'BARE') THEN
        WRITE(6, *) 'In BERTHA: unknown HMLTN value. Abnormal exit.'
        WRITE(7, *) 'In BERTHA: unknown HMLTN value. Abnormal exit.'
        WRITE(6, *) 'HMLTN = ',HLMTN
        WRITE(7, *) 'HMLTN = ',HLMTN
        STOP
      ENDIF
C
C     NEW SCF(1), OLD SCF(2), MBPT(3), MCSCF(4), DMRG(5), HMINT(6), PLOTS(7)
      READ(5,*) DUMLIN
      READ(5,*) ITREE
C      
C     ENSURE USER HAS SELECTED VALID CHOICE OF ITREE
      IF(ITREE.LT.1.OR.ITREE.GT.7) THEN
        WRITE(6, *) 'In BERTHA: invalid calculation tree. ',ITREE
        STOP
      ENDIF
C
      LF = LEN(TRIM(FILNAM))
C     DETERMINE LENGTH OF THIS FILNAM STRING, EXCLUDING TRAILING BLANKS
C      DO I=LEN(FILNAM),1,-1
C        IF(FILNAM(I:I).NE.' ') GOTO 40
C      ENDDO
C40    CONTINUE
C      LF = I
C
C     ADJUST OUTPUT FILE NAME AND SPECIFY DIRECTORY
      IF(ITREE.EQ.1.OR.ITREE.EQ.2) THEN
        STRING = 'output/'//FILNAM(:LF)//'_'//HMLTN//'_SCF'
      ELSEIF(ITREE.EQ.3) THEN
        STRING = 'output/'//FILNAM(:LF)//'_'//HMLTN//'_MBPT'
      ELSEIF(ITREE.EQ.4) THEN
        STRING = 'output/'//FILNAM(:LF)//'_'//HMLTN//'_MCSCF'
      ELSEIF(ITREE.EQ.5) THEN
        STRING = 'output/'//FILNAM(:LF)//'_'//HMLTN//'_DMRG'
      ELSEIF(ITREE.EQ.6) THEN
        STRING = 'output/'//FILNAM(:LF)//'_'//HMLTN//'_EXPVL'
      ELSEIF(ITREE.EQ.7) THEN
        STRING = 'output/'//FILNAM(:LF)//'_'//HMLTN//'_PLOTS'
      ENDIF
C
      LN = LEN(TRIM(STRING))
C     DETERMINE LENGTH OF THIS STRING STRING, EXCLUDING TRAILING BLANKS
      DO I=LEN(STRING),1,-1
        IF(STRING(I:I).NE.' ') GOTO 41
      ENDDO
41    CONTINUE
      LN = I
C
C     BASIS SET TYPE: GEOMETRIC (1) OR OPTIMISED (2)
      READ(5,*) DUMLIN
      READ(5,*) INTYPE
C
C     NUMBER OF ATOMIC CENTERS
      READ(5,*) DUMLIN
      READ(5,*) NCNT
C
C     INITIATE LOOP OVER NUCLEAR CENTERS
      MLQN = 0
      NDIM = 0
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
C       CALCULATE SOME INTERMEDIATE ATOMIC PARAMETERS (Z, NKAP, RN)
        NKAP(ICNT) = 2*LMAX(ICNT) + 1
        ZNUC(ICNT) = DFLOAT(IZNUC(ICNT))
        IF(IZNUC(ICNT).EQ.1) THEN
          CNUC(ICNT) = 0.21248239171D+10
        ELSEIF(IZNUC(ICNT).EQ.8) THEN
          CNUC(ICNT) = 0.58631436655D+09
        ELSE
          CNUC(ICNT) = 8.36D-1*(AMASS(ICNT)**(1.0D0/3.0D0))
          CNUC(ICNT) = 1.50D+10*(5.29177249D-1/(CNUC(ICNT)+5.7D-1))**2
        ENDIF
C
C       UPDATE OVERALL MAXIMUM OCCURRING LQN
        IF(LMAX(ICNT).GT.MLQN) MLQN = LMAX(ICNT)
C
        READ(5,*) DUMLIN
C ***   INITIATE IF STATEMENT FOR TYPE OF BASIS FUNCTION
C >>>   GEOMETRIC BASIS FUNCTIONS
        IF(INTYPE.EQ.1) THEN
C         GENERATE THE EVEN TEMPERED ORBITAL EXPONENTS FOR EACH LQN
          DO LQN=0,LMAX(ICNT)
C           READ GENERATING PARAMETERS A, B AND NFUNCT
            READ(5,*) APARAM,BPARAM,NFUNCT(LQN+1,ICNT)
C
C           GENERATE NFUNCT BASIS EXPONENTS USING VARIABLE ZETA
            ZETA = APARAM
            DO IBAS=1,NFUNCT(LQN+1,ICNT)
              EXPSET(IBAS,LQN+1,ICNT) = ZETA
              ZETA = ZETA*BPARAM
            ENDDO
          ENDDO
C >>>   OPTIMISED EXPONENTS FROM A RECORDED LIST
        ELSEIF(INTYPE.EQ.2) THEN
C         READ IN THE OPTIMISED ORBITAL EXPONENTS FOR EACH LQN
          DO LQN=0,LMAX(ICNT)
C           READ NFUNCT
            READ(5,*) NFUNCT(LQN+1,ICNT)
C           READ BASIS EXPONENTS FROM A LIST
            DO IBAS=1,NFUNCT(LQN+1,ICNT)
              READ(5,*) EXPSET(IBAS,LQN+1,ICNT)
            ENDDO
          ENDDO
C ***   END IF STATEMENT FOR TYPE OF BASIS FUNCTION
        ENDIF
C
        DO LQN=0,LMAX(ICNT)
C         EXTEND DIMENSION OF FOCK MATRIX
          NDIM = NDIM + 4*(2*LQN+1)*NFUNCT(LQN+1,ICNT)
C         ASSIGN KAPPA VALUES
          IF(LQN.EQ.0) THEN
            KVALS(1      ,ICNT) =-1
          ELSE
            KVALS(2*LQN  ,ICNT) = LQN
            KVALS(2*LQN+1,ICNT) =-LQN-1
          ENDIF
        ENDDO
C
C     END LOOP OVER NUCLEAR CENTERS        
      ENDDO
C
C     GENERATE MINIMAL LIST OF BOYS INTEGRALS FOR USE IN FUNFX/RMAKE
      CALL GFINIT(4*MLQN+10)
C
C     TOTAL DIMENSION DEPENDING ON CHOICE OF HAMILTONIAN
      IF(HMLTN.EQ.'NORL') THEN
        NDIM   = NDIM/2
        NSHIFT = 0
      ELSE
        NSHIFT = NDIM/2
      ENDIF
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
C      NOCC = 0
C      DO IZ=1,NCNT
C        NOCC = NOCC + IZNUC(IZ) - IQNUC(IZ)
C      ENDDO
C
C *** INITIATE IF STATEMENT DEPENDING ON CLOSED/OPEN SHELLS
C >>> OPEN-SHELL MOLECULE
C     DFNOTE: ENABLE THIS OPTION
      goto 555
      IF(NOPN.NE.0) THEN
C       FRACTIONAL OCCUPANCY OF THE OPEN SHELL
        FOPEN = DFLOAT(NOELEC)/DFLOAT(NOPN)
C       LABELS FOR THE OPEN SHELL
        READ(5,*) DUMLIN
        READ(5,*) ALPHA,BETA,(IOPN(M),M=1,NOPN)
C       PRINT THE LABELS FOR THE CLOSED-SHELL SPINORS USING KNOWN
C       IDENTITY OF THE OPEN-SHELL SPINORS
        JCL = 1
        JOP = 1
        DO JCOUNT=1,NOCC
          IF(JCOUNT.NE.IOPN(JOP)) THEN
            ICLS(JCL) = JCOUNT
            JCL = JCL + 1
          ELSE
            JOP = JOP + 1
          ENDIF
        ENDDO
C >>> CLOSED-SHELL MOLECULE
      ELSE
C       LABEL THE CLOSED-SHELL ELECTRONS
        DO JCL=1,NCLS
          ICLS(JCL) = JCL
        ENDDO
C *** END IF STATEMENT FOR CLOSED/OPEN SHELLS
      ENDIF
555   continue
C
C     ALL SPINORS ARE CLOSED
      DO JCL=1,NCLS
        ICLS(JCL) = JCL
      ENDDO
C
C     LEVEL SHIFT PARAMETER FOR EACH INTEGRAL INCLUSION LEVEL
      READ(5,*) DUMLIN
      READ(5,*) SHLV1,SHLV2,SHLV3
C
C     STARTING STAGE OF INTEGRAL INCLUSION LEVEL (1-3)
      READ(5,*) DUMLIN
      READ(5,*) IIL
C
C     REASONS TO SKIP TO FINAL INTEGRAL INCLUSION LEVEL
      IF(NCNT.EQ.1.OR.HMLTN.EQ.'NORL'.OR.ITREE.NE.1) THEN
        IIL = 3
      ENDIF
C
C     USE THE STARTING STAGE TO DETERMINE THE CURRENT SHIFT FACTOR
      IF(IIL.EQ.1) THEN
        SHLV = SHLV1
      ELSEIF(IIL.EQ.2) THEN
        SHLV = SHLV2
      ELSEIF(IIL.EQ.3) THEN
        SHLV = SHLV3
      ELSE
        WRITE(6, *) 'In BERTHA: invalid starting stage. IIL = ',IIL
        STOP
      ENDIF
C
C     DAMPING FACTOR AND RELATIVE TRESHOLD FOR INITIATION OF DAMPING
      READ(5,*) DUMLIN
      READ(5,*) DAMP,DTHRESH
C
C     E-COEFFICIENTS BY BATCH (0), TO LARGE EXTERNAL FILE (1)
      READ(5,*) DUMLIN
      READ(5,*) IEQS
C
C     READ IN MOST RECENT EIGEN-SOLUTIONS UNLESS A NEW SCF CALCULATION
      IF(ITREE.NE.1) THEN
        OPEN(UNIT=10,FILE='output/'//FILNAM(:LF)//'_'//HMLTN//'_SCF.wfn'
     &                                                ,STATUS='UNKNOWN')
        REWIND(UNIT=10)
        DO I=1,NDIM
          READ(10, *) EIGEN(I),(C(J,I),J=1,NDIM)
        ENDDO
        CLOSE(UNIT=10)
      ENDIF
C
C     EXPECTATION VALUE CALCULATIONS: ELCMNPL,MAGDIPL ETC
      IF(ITREE.EQ.6) THEN
        READ(5,*) DUMLIN
        READ(5,*) NHMINT
        IF(NHMINT.LT.1.OR.NHMINT.GT.10) THEN
          WRITE(6, *) 'In BERTHA: invalid number of expectation values.'
          STOP
        ENDIF
        DO N=1,NHMINT
          READ(5,*) HMINT(N)
        ENDDO
      ENDIF
C
C     DATA PLOTTING: AMPLITUDES, EM FIELDS AND POTENTIALS
      IF(ITREE.EQ.7) THEN
        READ(5,*) DUMLIN
        READ(5,*) NPTYPE
        IF(NPTYPE.LT.1.OR.NPTYPE.GT.10) THEN
          WRITE(6, *) 'In BERTHA: invalid number of plot types.'
          STOP
        ENDIF
        DO N=1,NPTYPE
          READ(5,*) PTYPE(N)
        ENDDO
      ENDIF
C
      RETURN
      END

