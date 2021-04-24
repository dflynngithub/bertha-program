      SUBROUTINE ERIRE0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           EEEEEEEE RRRRRRR  IIII RRRRRRR  EEEEEEEE 000000            C
C           EE       RR    RR  II  RR    RR EE      00    00           C
C           EE       RR    RR  II  RR    RR EE      00    00           C
C           EEEEEE   RR    RR  II  RR    RR EEEEEE  00    00           C
C           EE       RRRRRRR   II  RRRRRRR  EE      00    00           C
C           EE       RR    RR  II  RR    RR EE      00    00           C
C           EEEEEEEE RR    RR IIII RR    RR EEEEEEEE 000000            C
C                                                                      C
C -------------------------------------------------------------------- C
C  ERIRE0 EVALUATES A DIRECT AND EXCHANGE BATCH OF ELECTRON REPULSION  C
C  INTEGRALS OF ALL COMPONENT LABEL COMBINATIONS L AND S IN THE ATOMIC C
C  (RELATIVISTIC) SCF PROCEDURE.                                       C
C -------------------------------------------------------------------- C
C  OUTPUT:                                                             C
C    RJLLLL(M,N) - DIRECT ERI OVERLAP {LL,LL}                          C
C    RJSSSS(M,N) - DIRECT ERI OVERLAP {SS,SS}                          C
C    RJLLSS(M,N) - DIRECT ERI OVERLAP {LL,SS}                          C
C    RJSSLL(M,N) - DIRECT ERI OVERLAP {SS,LL}                          C
C    RKLLLL(M,N) - EXCHANGE ERI OVERLAP {LL,LL}                        C
C    RKSSSS(M,N) - EXCHANGE ERI OVERLAP {SS,SS}                        C
C    RKSLSL(M,N) - EXCHANGE ERI OVERLAP {SL,SL}                        C
C -------------------------------------------------------------------- C
C    N = 1 - KQN(A) > 0, KQN(B) > 0 (TYPICAL LABEL 11)                 C
C    N = 2 - KQN(A) < 0, KQN(B) > 0 (TYPICAL LABEL 12)                 C
C    N = 3 - KQN(A) > 0, KQN(B) < 0 (TYPICAL LABEL 21)                 C
C    N = 4 - KQN(A) < 0, KQN(B) < 0 (TYPICAL LABEL 22)                 C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MKP=9,MNU=MKP+1,MAB=2*MNU+5)
C
      COMMON/ANGL/BK(MNU,4),DK(MNU,4),HK(MNU,4),FK(MNU,4),GM(MNU,4),
     &            NUS(MNU),NUNUM,NUMIN,NUMAX
      COMMON/BLOC/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/CLRE/RKLLLL(MB2,4),RKSSSS(MB2,4),RKSLSL(MB2,4),
     &            RJLLLL(MB2,4),RJSSSS(MB2,4),RJLLSS(MB2,4),
     &            RJSSLL(MB2,4)
      COMMON/GMFN/GAMMAL(100),GAM(100)
      COMMON/IJ/EIJ(MAB),RNIJ(3),EI,EJ
      COMMON/KL/EK(MB2),EKL0(MB2),EIK(MB2,MAB),EJL(MB2,MAB),
     &          EL(MB2),EKL1(MB2),EKL(MB2,MAB),RNKL(MB2,4),
     &          B1(MB2,MAB,MAB),B2(MB2,MAB,MAB),B3(MB2,MAB,MAB),
     &          B4(MB2,MAB,MAB),B5(MB2,MAB,MAB),B6(MB2,MAB,MAB)
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
C     PREPARE VALUES FOR UPCOMING CALCULATIONS
      C1 = GAM(2*LQNA+2*LQNB+1)
      C3 = GAM(2*LQNA+2*LQNB+3)
      C5 = GAM(2*LQNA+2*LQNB+5)
      C7 = GAM(2*LQNA+2*LQNB+7)
      C9 = GAM(2*LQNA+2*LQNB+9)
C
      V1 = 1.0D0
      V2 = 2.0D0
      V4 = 4.0D0
      V8 = 8.0D0
      VS = 1.6D1
C
      F  = DFLOAT(2*LQNA+1)
      G  = DFLOAT(2*LQNB+1)
C
      F0G0 = 1.0D0
      F1G0 = F
      F0G1 = G
      F1G1 = F*G
      F2G0 = F*F
      F0G2 = G*G
      F2G1 = F*F*G
      F1G2 = F*G*G
      F2G2 = F*F*G*G
C
C**********************************************************************C
C     AN (LQNA,LQNB) COMBINATION HAS 1, 2 OR 4 (KQNA,KQNB) SUB-BLOCKS  C
C     SMALL-COMPONENT CONTRIBUTIONS DEPEND ON KQN SYMMETRY TYPE.       C
C**********************************************************************C
C
C     INITIATE LOOP OVER K,L BASIS FUNCTIONS
      DO M=1,MAXM
C
C       MORE VALUE PREPARATION
        E0000 = 1.0D0
        E1000 = EI
        E0100 = EJ
        E0010 = EK(M)
        E0001 = EL(M)
        E1100 = EI*EJ
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
C       TEMPORARY STORAGE OF VALUES
C       BXY MARKS 'B' BETA COMBINATION AND 'X' ONTO EFFECTIVE LQN STORE
        B00 = EIJ(2)*EKL(M,1)*B1(M,1,1) + EIJ(1)*EKL(M,2)*B2(M,1,1)
        B02 = EIJ(2)*EKL(M,3)*B1(M,1,2) + EIJ(1)*EKL(M,4)*B2(M,2,1)
        B04 = EIJ(2)*EKL(M,5)*B1(M,1,3) + EIJ(1)*EKL(M,6)*B2(M,3,1)
        B20 = EIJ(4)*EKL(M,1)*B1(M,2,1) + EIJ(3)*EKL(M,2)*B2(M,1,2)
        B22 = EIJ(4)*EKL(M,3)*B1(M,2,2) + EIJ(3)*EKL(M,4)*B2(M,2,2)
        B24 = EIJ(4)*EKL(M,5)*B1(M,2,3) + EIJ(3)*EKL(M,6)*B2(M,3,2)
        B40 = EIJ(6)*EKL(M,1)*B1(M,3,1) + EIJ(5)*EKL(M,2)*B2(M,1,3)
        B42 = EIJ(6)*EKL(M,3)*B1(M,3,2) + EIJ(5)*EKL(M,4)*B2(M,2,3)
        B44 = EIJ(6)*EKL(M,5)*B1(M,3,3) + EIJ(5)*EKL(M,6)*B2(M,3,3)
C
C       FILL RJ ARRAYS FOR THIS LQNA,LQNB BLOCK
C
C >>>   LQNA  =  0 AND LQNB  =  0 (REQUIRED FOR ALL BLOCKS)
        RJLLLL(M,4) = V1*F0G0*E0000*C5*B22
        RJLLSS(M,4) = V4*F0G0*E0011*C7*B24
        RJSSLL(M,4) = V4*F0G0*E1100*C7*B42
        RJSSSS(M,4) = VS*F0G0*E1111*C9*B44
C
C >>>   LQNA =/= 0                (NEED KQNA > 0 BLOCK)
        IF(LQNA.EQ.0) GOTO 103
        RJLLLL(M,3) = RJLLLL(M,4)
        RJLLSS(M,3) = RJLLSS(M,4)
        RJSSLL(M,3) = V4*F0G0*E1100*C7*B42
     &              - V2*F1G0*E1000*C5*B22 - V2*F1G0*E0100*C5*B22
     &              + V1*F2G0*E0000*C3*B02
        RJSSSS(M,3) = VS*F0G0*E1111*C9*B44
     &              - V8*F1G0*E1011*C7*B24 - V8*F1G0*E0111*C7*B24
     &              + V4*F2G0*E0011*C5*B04
103     CONTINUE
C
C >>>                  LQNB =/= 0 (NEED KQNB > 0 BLOCK)
        IF(LQNB.EQ.0) GOTO 102
        RJLLLL(M,2) = RJLLLL(M,4)
        RJLLSS(M,2) = V4*F0G0*E0011*C7*B24 
     &              - V2*F0G1*E0010*C5*B22 - V2*F0G1*E0001*C5*B22
     &              + V1*F0G2*E0000*C3*B20
        RJSSLL(M,2) = RJSSLL(M,4)
        RJSSSS(M,2) = VS*F0G0*E1111*C9*B44
     &              - V8*F0G1*E1110*C7*B42 - V8*F0G1*E1101*C7*B42
     &              + V4*F0G2*E1100*C5*B40
102     CONTINUE
C
C >>>   LQNA =/= 0 AND LQNB =/= 0 (NEED KQNA,KQNB > 0 BLOCK)
        IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 101
        RJLLLL(M,1) = RJLLLL(M,4)
        RJLLSS(M,1) = RJLLSS(M,2)
        RJSSLL(M,1) = RJSSLL(M,3)
        RJSSSS(M,1) = VS*F0G0*E1111*C9*B44
     &              - V8*F0G1*E1110*C7*B42 - V8*F0G1*E1101*C7*B42
     &              - V8*F1G0*E1011*C7*B24 - V8*F1G0*E0111*C7*B24
     &              + V4*F2G0*E0011*C5*B04 + V4*F0G2*E1100*C5*B40
     &              + V4*F1G1*E1001*C5*B22 + V4*F1G1*E0110*C5*B22
     &              + V4*F1G1*E0101*C5*B22 + V4*F1G1*E1010*C5*B22
     &              - V2*F1G2*E1000*C3*B20 - V2*F1G2*E0100*C3*B20
     &              - V2*F2G1*E0010*C3*B02 - V2*F2G1*E0001*C3*B02
     &              + V1*F2G2*E0000*C1*B00
101     CONTINUE
C
C**********************************************************************C
C       EXCHANGE INTEGRAL MATRICES: RKSSSS, RKSLSL, RKLLLL             C
C**********************************************************************C
C
C       LOOP OVER ORDERS FOR BETA INTEGRALS
        DO IV=1,NUNUM
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(IV)
C
C         INDEX OFFSETS FOR THE EXPONENT LISTS
          IA = NUMAX+NU+1
          IB = NUMAX-NU+1
C
C         INDEX OFFSETS FOR BETA INTEGRAL ARRAYS
          NX = (-NUMIN+NU)/2
          NY = ( NUMAX-NU)/2
C
C         TEMPORARY STORAGE OF VALUES
C         BXY MARKS 'B' BETA COMBINATION AND 'X' TO EFFECTIVE LQN STORE
C         ONTO WHICH NU IS ADDED OR SUBTRACTED.
          
          B00 = EIK(M,IA+1)*EJL(M,IB  )*B3(M,NX+1,NY+1)
     &        + EIK(M,IB  )*EJL(M,IA+1)*B4(M,NX+1,NY+1)
          B02 = EIK(M,IA+1)*EJL(M,IB+2)*B3(M,NX+1,NY+2)
     &        + EIK(M,IB  )*EJL(M,IA+3)*B4(M,NX+2,NY+1)
          B04 = EIK(M,IA+1)*EJL(M,IB+4)*B3(M,NX+1,NY+3)
     &        + EIK(M,IB  )*EJL(M,IA+5)*B4(M,NX+3,NY+1)
          B20 = EIK(M,IA+3)*EJL(M,IB  )*B3(M,NX+2,NY+1)
     &        + EIK(M,IB+2)*EJL(M,IA+1)*B4(M,NX+1,NY+2)
          B22 = EIK(M,IA+3)*EJL(M,IB+2)*B3(M,NX+2,NY+2)
     &        + EIK(M,IB+2)*EJL(M,IA+3)*B4(M,NX+2,NY+2)
          B24 = EIK(M,IA+3)*EJL(M,IB+4)*B3(M,NX+2,NY+3)
     &        + EIK(M,IB+2)*EJL(M,IA+5)*B4(M,NX+3,NY+2)
          B40 = EIK(M,IA+5)*EJL(M,IB  )*B3(M,NX+3,NY+1)
     &        + EIK(M,IB+4)*EJL(M,IA+1)*B4(M,NX+1,NY+3)
          B42 = EIK(M,IA+5)*EJL(M,IB+2)*B3(M,NX+3,NY+2)
     &        + EIK(M,IB+4)*EJL(M,IA+3)*B4(M,NX+2,NY+3)
          B44 = EIK(M,IA+5)*EJL(M,IB+4)*B3(M,NX+3,NY+3)
     &        + EIK(M,IB+4)*EJL(M,IA+5)*B4(M,NX+3,NY+3)
C
C >>>     LQNA  =  0 AND LQNB  =  0 (REQUIRED FOR ALL BLOCKS)
          RKLLLL(M,4) = RKLLLL(M,4) + BK(IV,4)*V1*F0G0*E0000*C5*B22
          RKSLSL(M,4) = RKSLSL(M,4) + BK(IV,4)*V4*F0G0*E1010*C7*B42
          RKSSSS(M,4) = RKSSSS(M,4) + BK(IV,4)*VS*F0G0*E1111*C9*B44
C
C >>>     LQNA =/= 0                (NEED KQNA > 0 BLOCK)
          IF(LQNA.EQ.0) GOTO 203
          RKLL = V1*F0G0*E0000*C5*B22
          RKSL = V4*F0G0*E1010*C7*B42 - V2*F1G0*E0010*C5*B22
          RKSS = VS*F0G0*E1111*C9*B44 - V8*F1G0*E1011*C7*B42 
     &         - V8*F1G0*E0111*C7*B24 + V4*F2G0*E0011*C5*B22
          RKLLLL(M,3) = RKLLLL(M,3) + BK(IV,3)*RKLL
          RKSLSL(M,3) = RKSLSL(M,3) + BK(IV,3)*RKSL
          RKSSSS(M,3) = RKSSSS(M,3) + BK(IV,3)*RKSS
203       CONTINUE
C
C >>>                    LQNB =/= 0 (NEED KQNB > 0 BLOCK)
          IF(LQNB.EQ.0) GOTO 202
          RKLL = V1*F0G0*E0000*C5*B22
          RKSL = V4*F0G0*E1010*C7*B42 - V2*F0G1*E1000*C5*B22
          RKSS = VS*F0G0*E1111*C9*B44 - V8*F0G1*E1110*C7*B42 
     &         - V8*F0G1*E1101*C7*B24 + V4*F0G2*E1100*C5*B22
          RKLLLL(M,2) = RKLLLL(M,2) + BK(IV,2)*RKLL
          RKSLSL(M,2) = RKSLSL(M,2) + BK(IV,2)*RKSL
          RKSSSS(M,2) = RKSSSS(M,2) + BK(IV,2)*RKSS
202       CONTINUE
C
C >>>     LQNA =/= 0 AND LQNB =/= 0 (NEED KQNA,KQNB > 0 BLOCK)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 201
          RKLL = V1*F0G0*E0000*C5*B22
          RKSL = V4*F0G0*E1010*C7*B42 - V2*F1G0*E0010*C5*B22 
     &         - V2*F0G1*E1000*C5*B22 + V1*F1G1*E0000*C3*B02
          RKSS = VS*F0G0*E1111*C9*B44
     &         - V8*F0G1*E1110*C7*B42 - V8*F0G1*E1101*C7*B24 
     &         - V8*F1G0*E1011*C7*B42 - V8*F1G0*E0111*C7*B24
     &         + V4*F2G0*E0011*C5*B22 + V4*F0G2*E1100*C5*B22 
     &         + V4*F1G1*E0110*C5*B22 + V4*F1G1*E1001*C5*B22
     &         + V4*F1G1*E1010*C5*B40 + V4*F1G1*E0101*C5*B04
     &         - V2*F2G1*E0010*C3*B20 - V2*F1G2*E1000*C3*B20
     &         - V2*F2G1*E0001*C3*B02 - V2*F1G2*E0100*C3*B02
     &         + V1*F2G2*E0000*C1*B00
          RKLLLL(M,1) = RKLLLL(M,1) + BK(IV,1)*RKLL
          RKSLSL(M,1) = RKSLSL(M,1) + BK(IV,1)*RKSL
          RKSSSS(M,1) = RKSSSS(M,1) + BK(IV,1)*RKSS 
201       CONTINUE
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     NORMALISE ACCUMULATED INTEGRALS (DIRECT AND EXCHANGE)            C
C**********************************************************************C
C
      DO M=1,MAXM
        T0LLLL = RNIJ(1)*RNKL(M,1)
        T0LLSS = RNIJ(1)*RNKL(M,3)
        T0SLSL = RNIJ(2)*RNKL(M,2)
        T0SSLL = RNIJ(3)*RNKL(M,1)
        T0SSSS = RNIJ(3)*RNKL(M,3)
        DO N=1,4
          RJLLLL(M,N) = RJLLLL(M,N)*T0LLLL
          RJLLSS(M,N) = RJLLSS(M,N)*T0LLSS
          RJSSLL(M,N) = RJSSLL(M,N)*T0SSLL
          RJSSSS(M,N) = RJSSSS(M,N)*T0SSSS
          RKLLLL(M,N) = RKLLLL(M,N)*T0LLLL
          RKSLSL(M,N) = RKSLSL(M,N)*T0SLSL
          RKSSSS(M,N) = RKSSSS(M,N)*T0SSSS
        ENDDO
      ENDDO
C
      RETURN
      END

