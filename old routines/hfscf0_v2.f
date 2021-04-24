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
      PARAMETER(MDM=2500,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26,M2B=2*MBS,
     &                        MB2=MBS*MBS,LWK=128*MBS,MNU=MKP+1,MIT=100)
C
      CHARACTER*1 ELLTERM,QSGN
      CHARACTER*2 ELMNT(120),ELNM
      CHARACTER*4 HMLTN
      CHARACTER*8 ZWRT,QWRT,EWRT
C
      COMPLEX*16 C(MDM,MDM)
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           VUEH(MDM,MDM),GDIR(MDM,MDM),GXCH(MDM,MDM),
     &           QDIR(MDM,MDM),QXCH(MDM,MDM),BDIR(MDM,MDM),
     &           BXCH(MDM,MDM),FOCK(MDM,MDM)
C    
      COMPLEX*16 BMAT(MDM,MDM)
C
      DIMENSION QE(MKP),QA(MKP),NORB(MKP,12),NUMOCC(MKP)
      DIMENSION W1(M2B),W2(M2B),T(LWK)
      DIMENSION O1(M2B,M2B),O2(M2B,M2B),H1(M2B,M2B),H2(M2B,M2B),
     &          U1(M2B,M2B),U2(M2B,M2B)
      DIMENSION G11(M2B,M2B),G21(M2B,M2B),G12(M2B,M2B),G22(M2B,M2B)
      DIMENSION B11(M2B,M2B),B21(M2B,M2B),B12(M2B,M2B),B22(M2B,M2B)
      DIMENSION DENLL(MB2,2*MKP+1),DFNLL(MB2,2*MKP+1),
     &          DENSL(MB2,2*MKP+1),DFNSL(MB2,2*MKP+1),
     &          DENSS(MB2,2*MKP+1),DFNSS(MB2,2*MKP+1)
      DIMENSION DEN1(MB2,3),DEN2(MB2,3)
C
      COMMON/ANGL/BK(MNU,4),DK(MNU,4),HK(MNU,4),FK(MNU,4),GM(MNU,4),
     &            NUS(MNU),NUNUM,NUMIN,NUMAX
      COMMON/ATOM/ELMNT
      COMMON/BLOC/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/BSIS/EXLA(MBS),EXLB(MBS)
      COMMON/COEF/C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/ENRG/ETOT,ENUC,EONE,ECLG,ECLQ,EBRG,EBRQ,EHNC,EHKN,EGDR,
     &            EGXC,EQDR,EQXC,EBDR,EBXC,EMDR,EMXC,EUEH
      COMMON/FILL/NCNF(MCT,MKP,MKP+1),NLVL(MCT,MKP),IFILL(MCT)
      COMMON/MTRX/OVAP,HNUC,HKIN,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/OCPD/IOCPN(MDM),IOCCM0
      COMMON/PRMS/CV,HMLTN,INEW,ITREE,IMOL,ILEV,IEQS
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR
C
      DATA EEPS/5.0D-13/
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
41    FORMAT(1X,I2,2X,F15.6,2X,F15.6,2X,F18.6,4X,ES12.5)
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
        RETURN
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
C       COEFFICIENT MATRIX ADDRESSES
        IL1 = LARGE(ICNT,2*LQNA+1,1)
        IL2 = LARGE(ICNT,2*LQNA+1,2)
        IS1 = IL1 + NSHIFT
        IS2 = IL2 + NSHIFT
C
C       TRANSFER OVERLAP MATRIX TO COMMON ARRAY
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNA
C
C           MATRIX ELEMENTS
            OVAP(IL1+IBAS,IL1+JBAS) = DCMPLX(O2(IBAS,JBAS),0.0D0)
            OVAP(IL2+IBAS,IL2+JBAS) = DCMPLX(O2(IBAS,JBAS),0.0D0)
C
            IF(HMLTN.NE.'NORL') THEN
C
C             SMALL COMPONENT BLOCK ADDRESSES
              KBAS = IBAS+NFUNA
              LBAS = JBAS+NFUNA
C
C             MATRIX ELEMENTS
              OVAP(IS1+IBAS,IS1+JBAS) = DCMPLX(O2(KBAS,LBAS),0.0D0)
              OVAP(IS2+IBAS,IS2+JBAS) = DCMPLX(O2(KBAS,LBAS),0.0D0)
C
            ENDIF
C
          ENDDO
        ENDDO
C
C       DIAGONALISE MATRIX (THIS NEEDS LAPACK LIBRARY)
        CALL DSYGV(1,'V','U',NMAT,H2,M2B,O2,M2B,W2,T,LWK,INFO)
        IF(INFO.NE.0) THEN
          WRITE(6, *) 'In HFSCF0: eigenvalue solver DSYGV failed.',INFO
          WRITE(7, *) 'In HFSCF0: eigenvalue solver DSYGV failed.',INFO
          STOP
        ENDIF
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
C**********************************************************************C
C     FIRST ITERATION ONLY INVOLVES BARE NUCLEUS PROBLEM               C
C**********************************************************************C
C
      IF(ITER.NE.1) THEN
C
C**********************************************************************C
C     SECOND LOOP: OVER BASIS FUNCTIONS K,L (USE INDEX 200)            C
C**********************************************************************C
C
C       LOOP OVER ALL OCCUPIED LQN VALUES
        DO 200 LQNB=0,LMXCONF
C
C       RECORD LQNB VALUE AND READ BASIS FUNCTIONS FOR THIS LQN
        NFUNB = NFUNCT(LQNB+1,ICNT)
        DO JBAS=1,NFUNB
          EXLB(JBAS) = EXPSET(JBAS,LQNB+1,ICNT)
        ENDDO
C
C       NUMBER OF BASIS FUNCTION OVERLAPS IN THIS BLOCK
        MAXM = NFUNB*NFUNB
C
C**********************************************************************C
C     GENERATE ATOMIC FOCK MATRIX                                      C
C**********************************************************************C
C
C       EVALUATE CLOSED-SHELL ELECTRON REPULSION ANGULAR INTEGRALS
        CALL ANGCOUL
C
C >>>   POSITIVE KAPPA(B) CHOICE (APPLIES ONLY FOR LQNB > 0)
        IF(LQNB.EQ.0.OR.HMLTN.EQ.'NORL') GOTO 230
C
C       KAPPA(B) VALUE AND DEGENERACY
        KAPB1 = LQNB
        RK2B1 = DFLOAT(2*IABS(KAPB1))
C
C       RELEVANT DENSITY MATRIX DEPENDS ON LQNA AND LQNB
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
230     CONTINUE
C
C >>>   NEGATIVE KAPPA(B) CHOICE (APPLIES TO ALL LQNB VALUES)
C
C       KAPPA(B) VALUE AND DEGENERACY
        KAPB2 =-LQNB-1
        IF(HMLTN.EQ.'NORL') THEN
          RK2B2 = DFLOAT(4*LQNB+2)
        ELSE
          RK2B2 = DFLOAT(2*IABS(KAPB2))
        ENDIF
C
C       RELEVANT DENSITY MATRIX DEPENDS ON LQNA AND LQNB
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
C       GENERATE THE MEAN-FIELD ATOMIC COULOMB MATRIX OVER DENSITIES
        CALL COULOMB0(G11,G21,G12,G22,DEN1,DEN2)
C
C       ADD COULOMB MATRIX ELEMENTS CONTRIBUTIONS TO FOCK MATRIX
        DO IBAS=1,NMAT
          DO JBAS=1,NMAT
C
            IF(HMLTN.EQ.'NORL') THEN
C           NON-RELATIVISTIC COULOMB MATRIX
C
              H2(IBAS,JBAS) = H2(IBAS,JBAS) + RK2B2*G22(IBAS,JBAS)
C
            ELSE
C           RELATIVISTIC COULOMB MATRIX
C
              H1(IBAS,JBAS) = H1(IBAS,JBAS) + RK2B1*G11(IBAS,JBAS)
     &                                      + RK2B2*G12(IBAS,JBAS)
              H2(IBAS,JBAS) = H2(IBAS,JBAS) + RK2B1*G21(IBAS,JBAS)
     &                                      + RK2B2*G22(IBAS,JBAS)
C
            ENDIF
C
          ENDDO
        ENDDO
C
C       TWO-BODY EIGENVALUE ENERGIES FOR OCCUPIED ELECTRONS
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNA
            M = M+1
C
            IF(HMLTN.EQ.'NORL') THEN
C           NON-RELATIVISTIC COULOMB ENERGY
C
              EG = EG + RK2A2*RK2B2*G22(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
C
            ELSE
C           RELATIVISTIC  COULOMB ENERGY
C
C             SMALL COMPONENT BLOCK ADDRESSES
              KBAS = IBAS+NFUNA
              LBAS = JBAS+NFUNA
C
C             LL BLOCK
              EG = EG
     &           +       RK2A1*RK2B1*G11(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &           +       RK2A1*RK2B2*G12(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &           +       RK2A2*RK2B1*G21(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &           +       RK2A2*RK2B2*G22(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
C
C             SL BLOCK
              EG = EG
     &           + 2.0D0*RK2A1*RK2B1*G11(KBAS,JBAS)*DFNSL(M,2*LQNA  )
     &           + 2.0D0*RK2A1*RK2B2*G12(KBAS,JBAS)*DFNSL(M,2*LQNA  )
     &           + 2.0D0*RK2A2*RK2B1*G21(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
     &           + 2.0D0*RK2A2*RK2B2*G22(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
C
C             SS BLOCK
              EG = EG
     &           +       RK2A1*RK2B1*G11(KBAS,LBAS)*DFNSS(M,2*LQNA  )
     &           +       RK2A1*RK2B2*G12(KBAS,LBAS)*DFNSS(M,2*LQNA  )
     &           +       RK2A2*RK2B1*G21(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
     &           +       RK2A2*RK2B2*G22(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
            ENDIF
C
          ENDDO
        ENDDO
C
C       GENERATE THE MEAN-FIELD ATOMIC BREIT MATRIX
        goto 250
        IF(HMLTN.NE.'DHFB'.AND.HMLTN.NE.'DHFQ') GOTO 250
C
C       EVALUATE CLOSED-SHELL BREIT INTERACTION ANGULAR INTEGRALS
        CALL ANGBREIT
C
C       GENERATE THE MEAN-FIELD ATOMIC BREIT MATRIX OVER DENSITIES
        CALL BREIT0(B11,B22,B12,B22,DEN1,DEN2)
C
C       ADD TWO-PARTICLE CONTRIBUTIONS TO FOCK MATRIX
        DO IBAS=1,2*NFUNA
          DO JBAS=1,2*NFUNA
            H1(IBAS,JBAS) = H1(IBAS,JBAS) + RK2B1*B11(IBAS,JBAS)
     &                                    + RK2B2*B12(IBAS,JBAS)
            H2(IBAS,JBAS) = H2(IBAS,JBAS) + RK2B1*B21(IBAS,JBAS)
     &                                    + RK2B2*B22(IBAS,JBAS)
          ENDDO
        ENDDO
C
C       TWO-BODY EIGENVALUE ENERGIES FOR OCCUPIED ELECTRONS
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNA
            M = M+1
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NFUNA
            LBAS = JBAS+NFUNA
C
C           LL BLOCK
            EB = EB
     &         +       RK2A1*RK2B1*B11(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &         +       RK2A1*RK2B2*B12(IBAS,JBAS)*DFNLL(M,2*LQNA  )
     &         +       RK2A2*RK2B1*B21(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
     &         +       RK2A2*RK2B2*B22(IBAS,JBAS)*DFNLL(M,2*LQNA+1)
C
C           SL BLOCK
            EB = EB
     &         + 2.0D0*RK2A1*RK2B1*B11(KBAS,JBAS)*DFNSL(M,2*LQNA  )
     &         + 2.0D0*RK2A1*RK2B2*B12(KBAS,JBAS)*DFNSL(M,2*LQNA  )
     &         + 2.0D0*RK2A2*RK2B1*B21(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
     &         + 2.0D0*RK2A2*RK2B2*B22(KBAS,JBAS)*DFNSL(M,2*LQNA+1)
C
C           SS BLOCK
            EB = EB
     &         +       RK2A1*RK2B1*B11(KBAS,LBAS)*DFNSS(M,2*LQNA  )
     &         +       RK2A1*RK2B2*B12(KBAS,LBAS)*DFNSS(M,2*LQNA  )
     &         +       RK2A2*RK2B1*B21(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
     &         +       RK2A2*RK2B2*B22(KBAS,LBAS)*DFNSS(M,2*LQNA+1)
C
          ENDDO
        ENDDO
C
C       PUT BREIT MATRIX COMPONENTS INTO A BIGGER MATRIX
        IL1 = LARGE(ICNT,2*LQNA+1,1)
        IL2 = LARGE(ICNT,2*LQNA+1,2)
        IS1 = IL1 + NSHIFT
        IS2 = IL2 + NSHIFT

        JL1 = LARGE(ICNT,2*LQNA+1,1)
        JL2 = LARGE(ICNT,2*LQNA+1,2)
        JS1 = JL1 + NSHIFT
        JS2 = JL2 + NSHIFT
C
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNA
C
C           SMALL COMPONENT BLOCK ADDRESSES
            KBAS = IBAS+NFUNA
            LBAS = JBAS+NFUNA
C
cc          BMAT(IL1+IBAS,JL1+JBAS) = 2.0D0*B22(IBAS      ,JBAS      )
cc          BMAT(IL2+IBAS,JL2+JBAS) = 2.0D0*B22(IBAS      ,JBAS      )
cc          BMAT(IL1+IBAS,JS1+JBAS) = 0.8d0*B22(IBAS      ,JBAS+NFUNA)
cc          BMAT(IL2+IBAS,JS2+JBAS) = 0.8d0*B22(IBAS      ,JBAS+NFUNA)
cc          BMAT(IS1+IBAS,JL1+JBAS) = 2.0D0*B22(IBAS+NFUNA,JBAS      )
cc          BMAT(IS2+IBAS,JL2+JBAS) = 2.0D0*B22(IBAS+NFUNA,JBAS      )
cc          BMAT(IS1+IBAS,JS1+JBAS) = 2.0D0*B22(IBAS+NFUNA,JBAS+NFUNA)
cc          BMAT(IS2+IBAS,JS2+JBAS) = 2.0D0*B22(IBAS+NFUNA,JBAS+NFUNA)
C
            BMAT(IL1+IBAS,JL1+JBAS) = DCMPLX(RK2B2*B22(IBAS,JBAS),0.0D0)
            BMAT(IL2+IBAS,JL2+JBAS) = DCMPLX(RK2B2*B22(IBAS,JBAS),0.0D0)
            BMAT(IL1+IBAS,JS1+JBAS) = DCMPLX(RK2B2*B22(IBAS,LBAS),0.0D0)
            BMAT(IL2+IBAS,JS2+JBAS) = DCMPLX(RK2B2*B22(IBAS,LBAS),0.0D0)
            BMAT(IS1+IBAS,JL1+JBAS) = DCMPLX(RK2B2*B22(KBAS,JBAS),0.0D0)
            BMAT(IS2+IBAS,JL2+JBAS) = DCMPLX(RK2B2*B22(KBAS,JBAS),0.0D0)
            BMAT(IS1+IBAS,JS1+JBAS) = DCMPLX(RK2B2*B22(KBAS,LBAS),0.0D0)
            BMAT(IS2+IBAS,JS2+JBAS) = DCMPLX(RK2B2*B22(KBAS,LBAS),0.0D0)
C
          ENDDO
        ENDDO
C
C
250     CONTINUE
C
C       END LOOP OVER LQNS FOR ORBITAL B
200     CONTINUE
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
C     FINISH GENERATING ATOMIC FOCK MATRIX, END CONDITIONAL OVER ITER
      ENDIF
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
C     BEGIN LOOP OVER MQNA VALUES
      DO IMVAL=1,IABS(KAPA1)
C
C       COEFFICIENT MATRIX ADDRESSES
        IL1 = LARGE(ICNT,2*LQNA  ,IMVAL*2-1)
        IL2 = LARGE(ICNT,2*LQNA  ,IMVAL*2  )
        IS1 = IL1 + NSHIFT
        IS2 = IL2 + NSHIFT
C
C
C       TRANSFER OVERLAP MATRIX TO COMMON ARRAY
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNA
C
C           MATRIX ELEMENTS
            OVAP(IL1+IBAS,IL1+JBAS) = DCMPLX(O1(IBAS,JBAS),0.0D0)
            OVAP(IL2+IBAS,IL2+JBAS) = DCMPLX(O1(IBAS,JBAS),0.0D0)
C
            IF(HMLTN.NE.'NORL') THEN
C
C             SMALL COMPONENT BLOCK ADDRESSES
              KBAS = IBAS+NFUNA
              LBAS = JBAS+NFUNA
C
C             MATRIX ELEMENTS
              OVAP(IS1+IBAS,IS1+JBAS) = DCMPLX(O1(KBAS,LBAS),0.0D0)
              OVAP(IS2+IBAS,IS2+JBAS) = DCMPLX(O1(KBAS,LBAS),0.0D0)
C
            ENDIF
C
          ENDDO
        ENDDO
      ENDDO
C
C     DIAGONALISE FOCK MATRIX (THIS NEEDS LAPACK LIBRARY)
      CALL DSYGV(1,'V','U',NMAT,H1,M2B,O1,M2B,W1,T,LWK,INFO)
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
C >>> NEGATIVE KAPPA(A) CHOICE (APPLIES TO ALL LQNA VALUES)
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
C       TRANSFER OVERLAP MATRIX TO COMMON ARRAY
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNA
C
C           MATRIX ELEMENTS
            OVAP(IL1+IBAS,IL1+JBAS) = DCMPLX(O2(IBAS,JBAS),0.0D0)
            OVAP(IL2+IBAS,IL2+JBAS) = DCMPLX(O2(IBAS,JBAS),0.0D0)
C
            IF(HMLTN.NE.'NORL') THEN
C
C             SMALL COMPONENT BLOCK ADDRESSES
              KBAS = IBAS+NFUNA
              LBAS = JBAS+NFUNA
C
C             MATRIX ELEMENTS
              OVAP(IS1+IBAS,IS1+JBAS) = DCMPLX(O2(KBAS,LBAS),0.0D0)
              OVAP(IS2+IBAS,IS2+JBAS) = DCMPLX(O2(KBAS,LBAS),0.0D0)
C
            ENDIF
C
          ENDDO
        ENDDO
      ENDDO
C
C     DIAGONALISE FOCK MATRIX (THIS NEEDS LAPACK LIBRARY)
      CALL DSYGV(1,'V','U',NMAT,H2,M2B,O2,M2B,W2,T,LWK,INFO)
      IF(INFO.NE.0) THEN
        WRITE(6, *) 'In HFSCF0: eigenvalue solver DSYGV failed.',INFO
        WRITE(7, *) 'In HFSCF0: eigenvalue solver DSYGV failed.',INFO
        STOP
      ENDIF
C
C     ATOMIC SELECTION RULE: ORTHOGONALITY IN BLOCKS OF KQN -> KA = KB
C
C     BEGIN LOOP OVER MQNA VALUES
      DO IMVAL=1,IABS(KVALS(2*LQNA+1,ICNT))
C
C       COEFFICIENT MATRIX ADDRESSES
        IL1 = LARGE(ICNT,2*LQNA+1,IMVAL*2-1)
        IL2 = LARGE(ICNT,2*LQNA+1,IMVAL*2  )
        IS1 = IL1 + NSHIFT
        IS2 = IL2 + NSHIFT
C
C       LOOP OVER ALL OCCUPIED PRINCIPLE SHELLS N FOR THIS KQNA
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
C         LOOP OVER ALL OCCUPIED PRINCIPLE SHELLS N FOR THIS KQNA
          DO IOCC=1,NUMOCC(LQNA+1)
C
C           EFFECTIVE OCCUPATION NUMBER
            QF = DSQRT(QA(IOCC))
C
C           LARGE COMPONENT OF KRAMERS PAIR
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
      ENDIF
C
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
        OPEN(UNIT=8,FILE='fock1.dat',STATUS='UNKNOWN')
        REWIND(UNIT=8)
        DO I=1,NDIM
          WRITE(8, *) (DREAL(BMAT(I,J)),J=1,NDIM)
        ENDDO
        CLOSE(UNIT=8)
        STOP
      ENDIF
C
      RETURN
      END