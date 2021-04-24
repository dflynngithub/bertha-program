      SUBROUTINE COULOMBNR0(G11,DEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C    CCCCCC  LL      MM       MM BBBBBBB  NN    NN RRRRRRR   000000    C
C   CC    CC LL      MMM     MMM BB    BB NNN   NN RR    RR 00    00   C
C   CC       LL      MMMM   MMMM BB    BB NNNN  NN RR    RR 00    00   C
C   CC       LL      MM MM MM MM BBBBBBB  NN NN NN RR    RR 00    00   C
C   CC       LL      MM  MMM  MM BB    BB NN  NNNN RRRRRRR  00    00   C
C   CC    CC LL      MM   M   MM BB    BB NN   NNN RR    RR 00    00   C
C    CCCCCC  LLLLLLL MM       MM BBBBBBB  NN    NN RR    RR  000000    C
C                                                                      C
C -------------------------------------------------------------------- C
C  COULOMBNR0 CONSTRUCTS THE PAULI ATOMIC COULOMB MATRIX FROM LISTS    C
C  OF INTEGRALS AND THE DENSITY MATRIX.                                C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MNU=MKP+1,MAB=2*MNU+5)
C
      DIMENSION G11(2*MBS,2*MBS),DEN(MB2),RN(MB2,4)
C
      COMMON/BLOC/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/BSIS/EXLA(MBS),EXLB(MBS)
      COMMON/CLNR/RJLLLL(MB2),RKLLLL(MB2)
      COMMON/IJ/EIJ(MAB),RNIJ(3),EI,EJ
      COMMON/KL/EK(MB2),EKL0(MB2),EIK(MB2,MAB),EJL(MB2,MAB),
     &          EL(MB2),EKL1(MB2),EKL(MB2,MAB),RNKL(MB2,4),
     &          B1(MB2,MAB,MAB),B2(MB2,MAB,MAB),B3(MB2,MAB,MAB),
     &          B4(MB2,MAB,MAB),B5(MB2,MAB,MAB),B6(MB2,MAB,MAB)
C
C     INITIALISE FOCK MATRIX
      DO IBAS=1,2*MBS
        DO JBAS=1,2*MBS
          G11(IBAS,JBAS) = 0.0D0
        ENDDO
      ENDDO
C
C     GENERATE 'KL' EXPONENT POWERS FOR USE IN B INTEGRAL GENERATION
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
          IJ = IJ+1
C
          EJ   = EXLA(JBAS)
          EIJ0 = EI+EJ
          EIJR = DSQRT(EIJ0)
          EIJA = EIJ0**(-LQNA)
          DO N=1,6
            EIJ(N) = EIJA
            EIJA   = EIJA/EIJR
          ENDDO

          RNIJ(1) = RN(IJ,1)
C
C         GENERATE A BATCH OF BETA INTEGRALS
          CALL BINT
C
C         GENERATE A BATCH OF NON-REL COULOMB INTEGRALS (J AND K)
          CALL ERINR0
C
C         BUILD THE FOCK MATRIX
          DO M=1,MAXM
            G11(IBAS,JBAS) = G11(IBAS,JBAS) 
     &                     + RJLLLL(M)*DEN(M) + RKLLLL(M)*DEN(M)
          ENDDO
C
        ENDDO
      ENDDO
C
      RETURN
      END

