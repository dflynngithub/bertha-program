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
      COMMON/B0IJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ,IBAS,JBAS
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
