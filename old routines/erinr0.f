      SUBROUTINE ERINR0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          EEEEEEEE RRRRRRR  IIII NN    NN RRRRRRR   000000            C
C          EE       RR    RR  II  NNN   NN RR    RR 00    00           C
C          EE       RR    RR  II  NNNN  NN RR    RR 00    00           C
C          EEEEEE   RR    RR  II  NN NN NN RR    RR 00    00           C
C          EE       RRRRRRR   II  NN  NNNN RRRRRRR  00    00           C
C          EE       RR    RR  II  NN   NNN RR    RR 00    00           C
C          EEEEEEEE RR    RR IIII NN    NN RR    RR  000000            C
C                                                                      C
C -------------------------------------------------------------------- C
C  ERINR0 EVALUATES A BATCH OF TWO ELECTRON INTEGRALS FOR THE          C
C  NON-RELATIVISTIC SCF PROCEDURE.                                     C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MNU=MKP+1,MAB=2*MNU+5)
C
      COMMON/ANGL/BK(MNU,4),DK(MNU,4),HK(MNU,4),FK(MNU,4),GM(MNU,4),
     &            NUS(MNU),NUNUM,NUMIN,NUMAX
      COMMON/BLOC/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/CLNR/RJLLLL(MB2),RKLLLL(MB2)
      COMMON/GMFN/GAMMAL(100),GAM(100)
      COMMON/IJ/EIJ(MAB),RNIJ(3),EI,EJ
      COMMON/KL/EK(MB2),EKL0(MB2),EIK(MB2,MAB),EJL(MB2,MAB),
     &          EL(MB2),EKL1(MB2),EKL(MB2,MAB),RNKL(MB2,4),
     &          B1(MB2,MAB,MAB),B2(MB2,MAB,MAB),B3(MB2,MAB,MAB),
     &          B4(MB2,MAB,MAB),B5(MB2,MAB,MAB),B6(MB2,MAB,MAB)
C
C     EMPTY COUNTER ARRAYS FOR DIRECT AND EXCHANGE INTEGRALS
      DO M=1,MAXM
        RJLLLL(M) = 0.0D0
        RKLLLL(M) = 0.0D0
      ENDDO
C
C     PREPARE VALUES FOR UPCOMING LISTS
      C5 = GAM(2*LQNA+2*LQNB+5)
C
C     INITIATE LOOP OVER K,L BASIS FUNCTIONS
      DO M=1,MAXM
C
C**********************************************************************C
C     DIRECT INTEGRAL MATRIX: RJLLLL                                   C
C**********************************************************************C
C
C       TEMPORARY STORAGE OF VALUES
        V1111 = EIJ(4)*EKL(M,3)*B1(M,2,2) + EIJ(3)*EKL(M,4)*B2(M,2,2)
C
C       FILL RJ ARRAY FOR THIS LQNA,LQNB BLOCK
        RJLLLL(M) = C5*V1111
C
C**********************************************************************C
C      EXCHANGE INTEGRAL MATRIX: RKLLLL                                C
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
          W1111 = EIK(M,IR+2)*EJL(M,IS+1)*B3(M,NX+1,NY+1) +
     &            EIK(M,IS+1)*EJL(M,IR+2)*B4(M,NX+1,NY+1)
C
C         FILL RK ARRAY FOR THIS LQNA,LQNB BLOCK
          RKLLLL(M) = RKLLLL(M) - BK(IV,4)*C5*W1111
C
        ENDDO
C
      ENDDO
C
C**********************************************************************C
C     NORMALISE ACCUMULATED INTEGRALS (DIRECT AND EXCHANGE)            C
C**********************************************************************C
C
      DO M=1,MAXM
        T0LLLL = RNIJ(1)*RNKL(M,1)
        RJLLLL(M) = RJLLLL(M)*T0LLLL
        RKLLLL(M) = RKLLLL(M)*T0LLLL
      ENDDO
C
      RETURN
      END

