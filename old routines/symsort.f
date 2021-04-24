c       do this after scf iteration is done, but before printing to terminal
C       REWRITE COEFFICIENT MATRIX IN SYMMETRY-ADAPTED LIST
C       IF(NCNT.LE.2) THEN
C          NMVALS = 0
C          CALL SYMSORT(NSYMOC,MLABEL,NMVALS)
C       ENDIF

      SUBROUTINE SYMSORT(NSYMOC,MLABEL,NSYM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    SSSSSS  YY    YY MM       MM  SSSSSS   OOOOOO  RRRRRRR TTTTTTTT   C
C   SS    SS YY    YY MMM     MMM SS    SS OO    OO RR    RR   TT      C
C   SS        YY  YY  MMMM   MMMM SS       OO    OO RR    RR   TT      C
C    SSSSSS    YYYY   MM MMMMM MM  SSSSSS  OO    OO RR    RR   TT      C
C         SS    YY    MM  MMM  MM       SS OO    OO RRRRRRR    TT      C
C   SS    SS    YY    MM   M   MM SS    SS OO    OO RR    RR   TT      C
C    SSSSSS     YY    MM       MM  SSSSSS   OOOOOO  RR    RR   TT      C
C                                                                      C
C -------------------------------------------------------------------- C
C  HAAKON WROTE THIS ROUTINE. THE THEORY GOES LIKE THIS: IN AN         C
C  ATOMIC OR DIATOMIC MOLECULE, MJ IS A 'GOOD' QUANTUM NUMBER, AND     C
C  BERTHA GENERATES A DEGENERATE MANIFOLD UPON WHICH A SET OF          C
C  EIGENVALUE ENERGIES ARE THE SAME. THEN THE EXPANSION COEFFICIENTS   C
C  FORM A LINEAR COMBINATION WITHIN THAT SET (SAY, THE MJ = +/- 3/2    C
C  AND +/- 1/2 STATES OF A THE 2P_3/2 ORBITAL). THEREFORE TO OBTAIN    C
C  NICE, CLEAN ORBITALS OF PURE CHARACTER, WE NEED TO ROTATE THE SET   C
C  OF STATES BY AN ANGLE. ONCE THAT ANGLE IS DETERMINED, THE ROTATION  C
C  WILL BE THE SAME FOR ALL MATRIX ELEMENTS -- SO EXPECT THAT THIS     C
C  ROUTINE CAN BE WRITTEN EVEN MORE SIMPLY.                            C
C**********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26)
C
      COMPLEX*16 C(MDM,MDM)
C
      DIMENSION ISYM(MDM*2,MMV*2),ICOUNT(MMV*2),
     &          NSYMOC(MMV*2),MLABEL(MDM),JLABEL(MDM)
C
      COMMON/COEF/C
      COMMON/EIGN/EIGEN(MDM)
      COMMON/SHLL/ALPHA,BETA,FOPEN,ICLS(500),IOPN(6),NCLS,NOPN
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR,IOCCM0
C
      DATA PI/3.1415926535897932D0/
C
C     SENSITIVITY TOLERANCE PARAMETERS
      TOL1 = 1.0D-14
      TOL2 = 1.0D-10
C
C**********************************************************************C
C     TAKE MQN PAIRS OF C AND ROTATE BETWEEN THEM TO SEPARATE OUT +/-. C
C**********************************************************************C
C
C     LOOP OVER PAIRS OF STATES
      DO IPAIR=NSHIFT+1,NDIM,2
C
C       TEMPORARY LARGE VALUE
        RLRG = 10.0D10
C
C       INITIAL INCRIMENTAL RADIAN (SWEEP OVER ALL POSSIBLE ANGLES)
        RINC = 2*PI/DFLOAT(360)
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
            ROT1 = CPHI*C(I,IPAIR  ) + SPHI*C(I,IPAIR+1)
            ROT2 =-SPHI*C(I,IPAIR  ) + CPHI*C(I,IPAIR+1)
            OLAP = OLAP + DABS(ROT1)*DABS(ROT2)
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
            ROT1 = CPHI*C(I,IPAIR  ) + SPHI*C(I,IPAIR+1)
            ROT2 =-SPHI*C(I,IPAIR  ) + CPHI*C(I,IPAIR+1)
            OLAP = OLAP + DABS(ROT1)*DABS(ROT2)
          ENDDO
C
C         IF THE NEW VALUE IS BELOW A TOLERANCE, FINISH.
          IF(DABS(OLAP-SOLD).LT.TOL1) GOTO 1
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
          ROT1 = CPHI*C(I,IPAIR  ) + SPHI*C(I,IPAIR+1)
          ROT2 =-SPHI*C(I,IPAIR  ) + CPHI*C(I,IPAIR+1)
C
          C(I,IPAIR  ) = ROT1
          C(I,IPAIR+1) = ROT2
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CHARACTERISE ALL OCCUPIED ORBITALS BY MQN TYPE WITH SIGN         C
C**********************************************************************C
C
C     EMPTY THE ICOUNT AND ISYM ARRAYS
      DO MA=1,MMV*2
        ICOUNT(MA) = 0
        DO N=1,MDM*2
          ISYM(N,MA) = 0
        ENDDO
      ENDDO
C
C     LOOP OVER ALL OCCUPIED ORBITALS (NOTE: USED TO GO UP TO NDIM)
      DO IOCC=NSHIFT+1,NDIM
C
C       LOOP OVER ALL MQNS
        DO MA=1,MMV
C
C         THE ACTUAL MQN IS HALF OF THIS
          MQN = 2*MA-1
C
C         INITIALISE NEGATIVE SPIN AND POSITIVE SPIN COUNTERS
          BINN = 0.0D0
          BINP = 0.0D0
C
C         LOOP OVER ALL NUCLEAR CENTERS
          DO ICNT=1,NCNT
C
C           LOOP OVER THE KAPPA VALUES OF EACH CENTER
            DO KA=1,NKAP(ICNT)
C
C             QUANTUM NUMBERS AND BASIS EXPONENTS FOR THIS SPINOR
              KQN = KVALS(KA,ICNT)
              IF(KQN.GT.0) THEN
                LQN = KQN
              ELSE
                LQN =-KQN-1
              ENDIF
C
              NFUN = NFUNCT(LQN+1,ICNT)
C
C             IF THE CURRENT MQN EXCEEDS KQN POSSIBILITIES MOVE ON
              IF(MQN.GT.2*IABS(KQN)-1) GOTO 2
C
C             FIND THE STARTING POINT IN THE FOCK MATRIX FOR THESE QNMS
              INDX1 = LARGE(ICNT,KA,2*MA-1)
              INDX2 = LARGE(ICNT,KA,2*MA  )
C
C             COUNT UP TOTAL CONTRIBUTIONS TO -VE AND +VE SPIN MQN
              DO M=1,NFUN
                BINN = BINN + DCONJG(C(INDX1+M,IOCC))*C(INDX1+M,IOCC)
                BINP = BINP + DCONJG(C(INDX2+M,IOCC))*C(INDX2+M,IOCC)
              ENDDO
C
2             CONTINUE
            ENDDO
          ENDDO
C
C         IF -SPIN COUNTER IS SMALL IT DOESN'T QUALIFY FOR -MQN STATUS
          IF(BINN.LE.TOL2) THEN
            GOTO 3
C
C         IF -SPIN COUNTER IS ABOVE THE TOLERANCE IT DOES QUALIFY
          ELSE
            ICOUNT(MA)          = ICOUNT(MA)+1
            ISYM(ICOUNT(MA),MA) = IOCC
            JLABEL(IOCC)        =-MQN
            GOTO 5
          ENDIF
3         CONTINUE
C
C         IF +SPIN COUNTER IS SMALL IT DOESN'T QUALIFY FOR +MQN STATUS
          IF(BINP.LE.TOL2) THEN
            GOTO 4
C
C         IF +SPIN COUNTER IS ABOVE THE TOLERANCE IT DOES QUALIFY
          ELSE
            ICOUNT(MA)          = ICOUNT(MA)+1
            ISYM(ICOUNT(MA),MA) = IOCC
            JLABEL(IOCC)        = MQN
            GOTO 5
          ENDIF
4         CONTINUE
C
C       END LOOP OVER MQNS
        ENDDO
C
C       IF ORBITAL CAN'T BE CLASSIFIED THEN ROTATION MUST HAVE FAILED
        WRITE(6, *) 'In SYMSORT: couldnt classify orbital number',IOCC
        WRITE(7, *) 'In SYMSORT: couldnt classify orbital number',IOCC
        STOP
C
5       CONTINUE
      ENDDO
C
C**********************************************************************C
C     REARRANGE THE MATRIX APPROPRIATELY                               C
C**********************************************************************C
C
      IF(NSYM.GT.0) THEN
C       REARRANGE COEFFICIENT MATRIX SO THAT THE FIRST NOCC ORBITALS
C       CONFORM TO INPUT SYMMETRY
        IOCC = 0
        DO MSYM=1,NSYM
          MQN = 2*MSYM - 1
          DO ISYMOC=1,NSYMOC(MSYM)
            LABVEC       = ISYM(ISYMOC,MSYM)
            IOCC         = IOCC+1
            MLABEL(IOCC) = MQN
            EIGEN(IOCC)  = EIGEN(LABVEC)
            DO M=1,NDIM
              C(M,IOCC) = C(M,LABVEC)
            ENDDO
          ENDDO
        ENDDO
C
        DO MVEC=1,IOCC
          EIGEN(MVEC+NSHIFT) = EIGEN(MVEC)
          DO M=1,NDIM
            C(M,MVEC+NSHIFT) = C(M,MVEC)
          ENDDO
        ENDDO
        IF(IOCC.NE.NOCC) THEN
          WRITE(6,*) 'In SYMSORT: some eigenvectors have gone missing.'
          WRITE(7,*) 'In SYMSORT: some eigenvectors have gone missing.'
          STOP
        ENDIF
      ELSE
C       NSYM=0 REQUIRES ONLY THAT THE STATES ARE LABELLED BUT NOT
C       BLOCKED BY SYMMETRY IN THE COEFFICIENT MATRIX
        DO IOCC=NSHIFT+1,NDIM
          MLABEL(IOCC-NSHIFT) = JLABEL(IOCC)
        ENDDO
      ENDIF
C     
      RETURN
      END
