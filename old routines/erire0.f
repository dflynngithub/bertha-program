      SUBROUTINE ERIRE0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C         EEEEEEEE RRRRRRR  IIIII RRRRRRR  EEEEEEEE  000000           C
C         EE       RR    RR  II   RR    RR EE       00    00          C
C         EE       RR    RR  II   RR    RR EE       00    00          C
C         EEEEEE   RR    RR  II   RR    RR EEEEEE   00    00          C
C         EE       RRRRRRR   II   RRRRRRR  EE       00    00          C
C         EE       RR    RR  II   RR    RR EE       00    00          C
C         EEEEEEEE RR    RR IIIII RR    RR EEEEEEEE  000000           C
C                                                                     C
C ------------------------------------------------------------------- C
C     ERIRE0 EVALUATES BATCHES OF TWO ELECTRON INTEGRALS FOR THE      C
C     RELATIVISTIC SCF PROCEDURE, WITH THE FOLLOWING ENTRIES:         C
C ------------------------------------------------------------------- C
C       11: KQN(A) > 0, KQN(B) > 0 -> R(M,1)                          C
C       12: KQN(A) < 0, KQN(B) > 0 -> R(M,2)                          C
C       21: KQN(A) > 0, KQN(B) < 0 -> R(M,3)                          C
C       22: KQN(A) < 0, KQN(B) < 0 -> R(M,4)                          C
C*********************************************************************C
      PARAMETER(MAXB=26,MAXB2=MAXB*MAXB,MAXAB=20,NUMAX=10)
C
      COMMON/ANGL/BK1(NUMAX),BK2(NUMAX),BK3(NUMAX),BK4(NUMAX),
     &            BK(NUMAX),NUS(NUMAX),NNUS
      COMMON/BLOC/NFUNA,NFUNB,LQNA,LQNB,MAXM
      COMMON/CLRE/RKLLLL(MAXB2,4),RKSSSS(MAXB2,4),RKSLSL(MAXB2,4),
     &            RJLLLL(MAXB2,4),RJSSSS(MAXB2,4),RJLLSS(MAXB2,4),
     &            RJSSLL(MAXB2,4)
      COMMON/GMT2/GAMMAL(100),GAM(100)
      COMMON/IJ/EIJ(MAXAB),RNIJ(3),EI,EJ
      COMMON/KL/EK(MAXB2),EKL0(MAXB2),EIK(MAXB2,MAXAB),EJL(MAXB2,MAXAB),
     &          EL(MAXB2),EKL1(MAXB2),EKL(MAXB2,MAXAB),RNKL(MAXB2,4),
     &          B1(MAXB2,MAXAB,MAXAB),B2(MAXB2,MAXAB,MAXAB),
     &          B3(MAXB2,MAXAB,MAXAB),B4(MAXB2,MAXAB,MAXAB)    
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
C     GAMMA FUNCTION LQN ARGUMENT
      LT = 2*(LQNA+LQNB)
      C1 = GAM(LT+1)
      C3 = GAM(LT+3)
      C5 = GAM(LT+5)
      C7 = GAM(LT+7)
      C9 = GAM(LT+9)
C
C     G REFERENCES TO THE (I,J) PAIR AND F = REFERENCES TO THE (K,L) PAIR
      G = DFLOAT(2*LQNA+1)
      F = DFLOAT(2*LQNB+1)     
C
C*********************************************************************C
C     SMALL COMPONENT IS SENSITIVE TO SYMMETRY TYPE FOR AN LQN VALUE  C
C     SO CALCULATE ALL CHOICES FOR LQN VALUES OTHER THAN LQN = 0      C
C*********************************************************************C
C
C *** BEGIN CONDITIONAL STATEMENT OVER LQNA AND LQNB TYPES
C
C >>> LQNA =/= 0 AND LQNB =/= 0
      IF(LQNA.NE.0.AND.LQNB.NE.0) THEN
C
C       PREPARE VALUES FOR UPCOMING LISTS
        T13  = G*G*F*2.0D0
        T45  = 4.0D0*G*F
        T58  = 8.0D0*G
        T68  = 2.0D0*EI
        T73  = 4.0D0*EI
        T90  = 8.0D0*EI*EJ
        T102 = 1.6D1*EI
        T150 = C5*G*G*4.0D0
        T155 = T58*EJ
        T160 = C3*F*F
        T165 = T58*EI
        T170 = C5*T73*EJ*F*F
        T175 = T90*F 
        T180 = C9*T102*EJ
        T185 = C1*G*G*F*F
        T190 = C3*T13
        T195 = C3*2.0D0*G*F*F*EJ
        T200 = T45*EJ 
        T205 = 2.0D0*T200  
        T210 = T68*G
        T215 = T73*G
        T220 = C5*F
        T225 = C5*G
        T230 = C9*EJ*T102
        T235 = C7*4.0D0
        T240 = C3*G*G
        T245 = C7*EJ
        T250 = 4.0D0*EI*EJ
        T255 = 2.0D0*EJ
C
C       DIRECT INTEGRAL MATRICES: RJSSSS, RJLLSS, RJSSLL, RJLLLL
        DO M=1,MAXM
          T19  = EIJ(2)*EKL(M,3)*B1(M,1,2) + EIJ(1)*EKL(M,4)*B2(M,2,1)
          T33  = EKL1(M)*T150*(EIJ(2)*EKL(M,5)*B1(M,1,3)
     &           + EIJ(1)*EKL(M,6)*B2(M,3,1))
          T41  = EIJ(4)*EKL(M,1)*B1(M,2,1) + EIJ(3)*EKL(M,2)*B2(M,1,2)
          T51  = EIJ(4)*EKL(M,3)*B1(M,2,2) + EIJ(3)*EKL(M,4)*B2(M,2,2)
          T52  = C5*T51
          T65  = EIJ(4)*EKL(M,5)*B1(M,2,3) + EIJ(3)*EKL(M,6)*B2(M,3,2)
          T66  = C7*EKL1(M)*T65
          T67  = T155*T66
          T71  = T160*T41
          T75  = F*EL(M)
          T79  = T165*T66
          T88  = T170*(EIJ(6)*EKL(M,1)*B1(M,3,1)
     &               + EIJ(5)*EKL(M,2)*B2(M,1,3))
          T95  = EIJ(6)*EKL(M,3)*B1(M,3,2) 
     &         + EIJ(5)*EKL(M,4)*B2(M,2,3)
          T96  = C7*T95
          T98  = T90*T75*T96
          T101 = T175*EK(M)*T96
          T109 = EIJ(6)*EKL(M,5)*B1(M,3,3) 
     &         + EIJ(5)*EKL(M,6)*B2(M,3,3)
          T111 = EKL1(M)*T180*T109
          T112 = T185*(EIJ(2)*EKL(M,1)*B1(M,1,1)
     &               + EIJ(1)*EKL(M,2)*B2(M,1,1))
     &                 - EL(M)*T190*T19-EK(M)*T190*T19 + T33 - T195*T41 
     &            + T200*EL(M)*T52+2*T45*EJ*EK(M)*T52 - T67 - T210*T71
     &            + T215*T75*T52 - T79 + T88 - T98 - T101 + T111
          T126 = T220*T51
          T142 = T225*T51
C
          RJLLLL(M,1) = T52
          RJLLSS(M,1) = T71-2.0D0*EK(M)*T126
     &                     -2.0D0*EL(M)*T126 + T235*EKL1(M)*T65
          RJSSLL(M,1) = T240*T19 - T255*T142 - T68*T142 + T73*T245*T95
          RJSSSS(M,1) = T112
C
          RJLLLL(M,2) = RJLLLL(M,1)
          RJLLSS(M,2) = RJLLSS(M,1)
          RJSSLL(M,2) = T250*T96
          RJSSSS(M,2) = T88-T98-T101+T111
C
          RJLLLL(M,3) = RJLLLL(M,1)
          RJLLSS(M,3) = 4.0D0*T66
          RJSSLL(M,3) = RJSSLL(M,1)
          RJSSSS(M,3) = T33-T67-T79+T111
C
          RJLLLL(M,4) = RJLLLL(M,1)
          RJLLSS(M,4) = RJLLSS(M,3)
          RJSSLL(M,4) = RJSSLL(M,2)
          RJSSSS(M,4) = EKL1(M)*T230*T109
        ENDDO
C
C       LOOP OVER ORDERS FOR BETA INTEGRALS
        DO INU=1,NNUS
          NU   = NUS(INU)
          IR1  = NU + NUS(NNUS) + 2
          IR2  = NUS(NNUS) - NU + 2
          NX   = (NU - NUS(1) + 2)/2
          NY   = (NUS(NNUS) - NU + 2)/2
C
C         PREPARE VALUES FOR UPCOMING LIST
          T1   = G**2
          T2   = F**2
          T12  = T1*F*2.0D0
          T44  = 4.0D0*G*F
          T58  = 8.0D0*G
          T68  = 2.0D0*EI
          T73  = 4.0D0*EI
          T100 = 8.0D0*EI*EJ
          T150 = C5*T1*4.0D0
          T155 = T58*EJ
          T160 = T58*EI
          T165 = C5*T73*EJ*G*G
          T170 = C7*T100
          T175 = C7*T100*F
          T180 = C9*1.6D1*EI*EJ
          T185 = C1*T1*G*G
          T190 = C3*T12
          T195 = C3*2.0D0*G*G*G*EJ
          T200 = C5*T44*EJ
          T205 = T44*EJ
          T210 = C3*T68*G*G*G
          T215 = C5*2.0D0*G
          T220 = C5*T68*F
          T225 = C3*G*F
          T230 = T73*G
          T235 = C5*T44*EI
          T240 = 4.0D0*EI
C
C         EXCHANGE INTEGRAL MATRICES: RKSSSS, RKSLSL, RKLLLL
          DO M=1,MAXM
            T18  = EIK(M,IR1  )*EJL(M,IR2+1)*B3(M,  NX,1+NY)
     &           + EIK(M,IR2-1)*EJL(M,IR1+2)*B4(M,1+NX,  NY)
            T26  = EIK(M,IR1+2)*EJL(M,IR2-1)*B3(M,1+NX,  NY)
     &           + EIK(M,IR2+1)*EJL(M,IR1  )*B4(M,  NX,1+NY)
            T35  = EIK(M,IR1+2)*EJL(M,IR2+1)*B3(M,1+NX,1+NY)
     &           + EIK(M,IR2+1)*EJL(M,IR1+2)*B4(M,1+NX,1+NY)
            T37  = EKL1(M)*T150*T35
            T55  = C5*T35
            T60  = C7*EKL1(M)
            T65  = EIK(M,IR1+2)*EJL(M,IR2+3)*B3(M,1+NX,2+NY)
     &           + EIK(M,IR2+1)*EJL(M,IR1+4)*B4(M,2+NX,1+NY)
            T67  = T155*T60*T65
            T75  = F*EL(M)
            T92  = EIK(M,IR1+4)*EJL(M,IR2+1)*B3(M,2+NX,1+NY)
     &           + EIK(M,IR2+3)*EJL(M,IR1+2)*B4(M,1+NX,2+NY)
            T94  = T160*T60*T92
            T98  = T165*T35
            T103 = T75*T170*T65
            T107 = EK(M)*T175*T92
            T115 = EIK(M,IR1+4)*EJL(M,IR2+3)*B3(M,2+NX,2+NY)
     &           + EIK(M,IR2+3)*EJL(M,IR1+4)*B4(M,2+NX,2+NY)
            T117 = EKL1(M)*T180*T115
            T118 = T185*(EIK(M,IR1)*EJL(M,IR2-1)*B3(M,NX,NY)
     &           + EIK(M,IR2-1)*EJL(M,IR1  )*B4(M,NX,NY)) 
     &           - EL(M)*T190*T18 - EK(M)*T190*T26 + T37 - T195*T18
     &           + EL(M)*T200*(EIK(M,IR1  )*EJL(M,IR2+3)*B3(M,NX,2+NY)
     &                       + EIK(M,IR2-1)*EJL(M,IR1+4)*B4(M,2+NX,NY))
     &           + EK(M)*T235*(EIK(M,IR1+4)*EJL(M,IR2-1)*B3(M,2+NX,NY)
     &                       + EIK(M,IR2+3)*EJL(M,IR1  )*B4(M,NX,2+NY))  
     &           + EK(M)*T205*T55 - T67 - T210*T26 + T230*T75*T55
     &                         - T94 + T98 - T103 - T107 + T117
            T136 = EK(M)*T215*T35
            T139 = T220*T35
            T141 = C7*EK(M)*T92
            T142 = T73*T141
C
            RKLLLL(M,1) = RKLLLL(M,1)+BK1(INU)*C5*T35
            RKSLSL(M,1) = RKSLSL(M,1)+BK1(INU)*(T225*T18-T136-T139+T142)
            RKSSSS(M,1) = RKSSSS(M,1)+BK1(INU)*T118
C 
            RKLLLL(M,2) = RKLLLL(M,2)+BK2(INU)*C5*T35
            RKSLSL(M,2) = RKSLSL(M,2)+BK2(INU)*(-T139+T142)
            RKSSSS(M,2) = RKSSSS(M,2)+BK2(INU)*(T98-T103-T107+T117)
C
            RKLLLL(M,3) = RKLLLL(M,3)+BK3(INU)*C5*T35
            RKSLSL(M,3) = RKSLSL(M,3)+BK3(INU)*(-T136+T142)
            RKSSSS(M,3) = RKSSSS(M,3)+BK3(INU)*(T37-T67-T94+T117)
C
            RKLLLL(M,4) = RKLLLL(M,4)+BK4(INU)*C5*T35
            RKSLSL(M,4) = RKSLSL(M,4)+BK4(INU)*T240*T141
            RKSSSS(M,4) = RKSSSS(M,4)+BK4(INU)*EKL1(M)*T180*T115
          ENDDO
        ENDDO
C
C >>> LQNA =/= 0 AND LQNB  =  0
      ELSEIF(LQNA.NE.0.AND.LQNB.EQ.0) THEN
C
C       PREPARE VALUES FOR UPCOMING LIST
        T12  = 8.0D0*G
        T24  = 1.6D1*EI
        T150 = C5*G
        T155 = C5*G*G*4.0D0
        T160 = T12*EJ
        T165 = T12*EI
        T170 = C9*T24*EJ
        T180 = 2.0D0*EJ
        T185 = 2.0D0*EI
        T190 = C7*4.0D0*EI*EJ
C
C       DIRECT INTEGRAL MATRICES: RJSSSS, RJLLSS, RJSSLL, RJLLLL
        DO M=1,MAXM
          T20 = C7*EKL1(M)
     &         *(EIJ(4)*EKL(M,5)*B1(M,2,3) + EIJ(3)*EKL(M,6)*B2(M,3,2))
          T31 =  EIJ(6)*EKL(M,5)*B1(M,3,3) + EIJ(5)*EKL(M,6)*B2(M,3,3)
          T58 =                              EIJ(4)*EKL(M,3)*B1(M,2,2)  
     &                                     + EIJ(3)*EKL(M,4)*B2(M,2,2)
          T59 = T150*T58
          T69 =  EIJ(6)*EKL(M,3)*B1(M,3,2) + EIJ(5)*EKL(M,4)*B2(M,2,3)
C
          RJLLLL(M,3) = C5*T58
          RJLLSS(M,3) = 4.0D0*T20 
          RJSSLL(M,3) = C3*G*G*(EIJ(2)*EKL(M,3)*B1(M,1,2)
     &                                 +EIJ(1)*EKL(M,4)*B2(M,2,1))
     &                                 -T180*T59 - T185*T59 + T190*T69
          RJSSSS(M,3) = EKL1(M)*T155*(EIJ(2)*EKL(M,5)*B1(M,1,3)
     &                              + EIJ(1)*EKL(M,6)*B2(M,3,1))
     &                 - T160*T20-T165*T20+T170*EKL1(M)*T31
C
          RJLLLL(M,4) = RJLLLL(M,3)
          RJLLSS(M,4) = RJLLSS(M,3)
          RJSSLL(M,4) = T190*T69
          RJSSSS(M,4) = EKL1(M)*T170*T31 
        ENDDO
C
C       LOOP OVER ORDERS FOR BETA INTEGRALS
        DO INU=1,NNUS
          NU   = NUS(INU)
          IR1  = NUS(NNUS) + NU + 2
          IR2  = NUS(NNUS) - NU + 2
          NX   = (NU-NUS(1)+2)/2
          NY   = (NUS(NNUS)-NU+2)/2
C
          T1   = G**2
          T12  = 8.0D0*G
          T150 = C5*T1*4.0D0
          T155 = T12*EJ
          T160 = T12*EI
          T165 = C9*1.6D1*EI*EJ
          T170 =-C5*2.0D0*G
          T175 = 4.0D0*EI 
C
C         EXCHANGE INTEGRAL MATRICES: RKSSSS, RKSLSL, RKLLLL
          DO M=1,MAXM
            T9  = EIK(M,IR1+2)*EJL(M,IR2+1)*B3(M,NX+1,NY+1)
     &          + EIK(M,IR2+1)*EJL(M,IR1+2)*B4(M,NX+1,NY+1)
            T14 = C7*EKL1(M)
            T27 = EIK(M,IR1+4)*EJL(M,IR2+1)*B3(M,NX+2,NY+1)
     &          + EIK(M,IR2+3)*EJL(M,IR1+2)*B4(M,NX+1,NY+2)
            T37 = EIK(M,IR1+4)*EJL(M,IR2+3)*B3(M,NX+2,NY+2)
     &          + EIK(M,IR2+3)*EJL(M,IR1+4)*B4(M,NX+2,NY+2)
            T54 = C7*EK(M)*T27
C
            RKLLLL(M,3) = RKLLLL(M,3) +BK3(INU)*C5*T9
            RKSLSL(M,3) = RKSLSL(M,3) +BK3(INU)*(EK(M)*T170*T9+T175*T54)
            RKSSSS(M,3) = RKSSSS(M,3) +BK3(INU)*(-T160*T14*T27+EKL1(M)*
     &                          T165*T37+EK(M)*EL(M)*T150*T9-T155*T14
     &                     *(EIK(M,IR1+2)*EJL(M,IR2+3)*B3(M,1+NX,2+NY)
     &                     +EIK(M,IR2+1)*EJL(M,IR1+4)*B4(M,2+NX,1+NY)))
C
            RKLLLL(M,4) = RKLLLL(M,4) + BK4(INU)*C5*T9
            RKSLSL(M,4) = RKSLSL(M,4) + BK4(INU)*T175*T54
            RKSSSS(M,4) = RKSSSS(M,4) + BK4(INU)*EKL1(M)*T165*T37
          ENDDO
        ENDDO
C
C >>> LQNA  =  0 AND LQNB =/= 0
      ELSEIF(LQNA.EQ.0.AND.LQNB.NE.0) THEN
C
C       PREPARE VALUES FOR UPCOMING LISTS
        F31 = C3*F*F
        F51 = C5
        F52 =-C5*2.0D0*F*EI
        F54 = C5*4.0D0*F*F*EI*EJ
        F71 = C7*4.0D0*EI*EJ
C
C       INITIATE LOOP OVER K,L BASIS FUNCTIONS
        DO M=1,MAXM
C
C         MORE VALUE PREPARATION
          F53 =-C5*2.0D0*F*EKL0(M)
          F72 = C7*4.0D0*EI*EK(M)
          F73 =-C7*8.0D0*F*EI*EJ*EL(M)
          F74 =-C7*8.0D0*F*EI*EJ*EK(M)
          F75 =-C7*8.0D0*F*EI*EJ*EKL0(M)
          F76 = C7*4.0D0*EKL1(M)
          F91 = C9*1.6D1*EI*EJ*EKL1(M)
C
C         DIRECT INTEGRAL MATRICES: RJSSSS, RJLLSS, RJSSLL, RJLLLL
C
C         TEMPORARY STORAGE OF VALUES
          R4121 = EIJ(4)*EKL(M,1)*B1(M,2,1) + EIJ(3)*EKL(M,2)*B2(M,1,2)
          R4322 = EIJ(4)*EKL(M,3)*B1(M,2,2) + EIJ(3)*EKL(M,4)*B2(M,2,2)
          R4523 = EIJ(4)*EKL(M,5)*B1(M,2,3) + EIJ(3)*EKL(M,6)*B2(M,3,2)
          R6131 = EIJ(6)*EKL(M,1)*B1(M,3,1) + EIJ(5)*EKL(M,2)*B2(M,1,3)
          R6332 = EIJ(6)*EKL(M,3)*B1(M,3,2) + EIJ(5)*EKL(M,4)*B2(M,2,3)
          R6533 = EIJ(6)*EKL(M,5)*B1(M,3,3) + EIJ(5)*EKL(M,6)*B2(M,3,3)         
C        
C         TRANSFER TO RJ ARRAYS
          RJLLLL(M,2) = F51*R4322
          RJLLSS(M,2) = F31*R4121 + F53*R4322 + F76*R4523
          RJSSLL(M,2) = F71*R6332
          RJSSSS(M,2) = F54*R6131 + F75*R6332 + F91*R6533
C
          RJLLLL(M,4) = F51*R4322
          RJLLSS(M,4) = F76*R4523
          RJSSLL(M,4) = F71*R6332
          RJSSSS(M,4) = F91*R6533
C
C         EXCHANGE INTEGRAL MATRICES: RKSSSS, RKSLSL, RKLLLL
C
C         LOOP OVER ORDERS FOR BETA INTEGRALS
          DO INU=1,NNUS
            NU   = NUS(INU)
            IR1  = NUS(NNUS) + NU + 2
            IR2  = NUS(NNUS) - NU + 2
            NX   = (-NUS(1)   +NU + 2)/2
            NY   = ( NUS(NNUS)-NU + 2)/2
C
C           TEMPORARY STORAGE OF VALUES
            R2111 = EIK(M,IR1+2)*EJL(M,IR2+1)*B3(M,NX+1,NY+1)
     &            + EIK(M,IR2+1)*EJL(M,IR1+2)*B4(M,NX+1,NY+1)
            R4121 = EIK(M,IR1+4)*EJL(M,IR2+1)*B3(M,NX+2,NY+1)
     &            + EIK(M,IR2+3)*EJL(M,IR1+2)*B4(M,NX+1,NY+2)
            R2312 = EIK(M,IR1+2)*EJL(M,IR2+3)*B3(M,NX+1,NY+2)
     &            + EIK(M,IR2+1)*EJL(M,IR1+4)*B4(M,NX+2,NY+1)
            R4322 = EIK(M,IR1+4)*EJL(M,IR2+3)*B3(M,NX+2,NY+2)
     &            + EIK(M,IR2+3)*EJL(M,IR1+4)*B4(M,NX+2,NY+2)
C
C           TRANSFER TO RK ARRAYS
            RKLLLL(M,2) = RKLLLL(M,2) + BK2(INU)*F51*R2111
            RKSLSL(M,2) = RKSLSL(M,2) + BK2(INU)*(F52*R2111 + F72*R4121)
            RKSSSS(M,2) = RKSSSS(M,2) + BK2(INU)*(F54*R2111 + F73*R2312
     &                                           +F74*R4121 + F91*R4322)
     
            RKLLLL(M,4) = RKLLLL(M,4) + BK4(INU)*F51*R2111
            RKSLSL(M,4) = RKSLSL(M,4) + BK4(INU)*F72*R4121
            RKSSSS(M,4) = RKSSSS(M,4) + BK4(INU)*F91*R4322
          ENDDO
C
        ENDDO      
C       
C >>> LQNA  =  0 AND LQNB  =  0
      ELSEIF(LQNA.EQ.0.AND.LQNB.EQ.0) THEN
C
C       PREPARE VALUES FOR UPCOMING LISTS
        F51 = C5
        F71 = C7*4.0D0*EI*EJ
        F91 = C9*1.6D1*EI*EJ
C
C       INITIATE LOOP OVER K,L BASIS FUNCTIONS
        DO M=1,MAXM
C
C         MORE VALUE PREPARATION
          F72 = C7*4.0D0*EKL1(M)
          F73 = C7*4.0D0*EI*EK(M)
          F92 = C9*1.6D1*EI*EJ*EKL1(M)
C
C         DIRECT INTEGRAL MATRICES: RJSSSS, RJLLSS, RJSSLL, RJLLLL
C
C         TEMPORARY STORAGE OF VALUES
          R4322 = EIJ(4)*EKL(M,3)*B1(M,2,2) + EIJ(3)*EKL(M,4)*B2(M,2,2)
          R4523 = EIJ(4)*EKL(M,5)*B1(M,2,3) + EIJ(3)*EKL(M,6)*B2(M,3,2)
          R6332 = EIJ(6)*EKL(M,3)*B1(M,3,2) + EIJ(5)*EKL(M,4)*B2(M,2,3)
          R6533 = EIJ(6)*EKL(M,5)*B1(M,3,3) + EIJ(5)*EKL(M,6)*B2(M,3,3)
C        
C         TRANSFER TO RJ ARRAYS
          RJLLLL(M,4) = F51*R4322
          RJLLSS(M,4) = F72*R4523
          RJSSLL(M,4) = F71*R6332
          RJSSSS(M,4) = F92*R6533
C
C         EXCHANGE INTEGRAL MATRICES: RKSSSS, RKSLSL, RKLLLL
C
C         LOOP OVER ORDERS FOR BETA INTEGRALS
          DO INU=1,NNUS
            NU   = NUS(INU)
            IR1  = NUS(NNUS) + NU + 2
            IR2  = NUS(NNUS) - NU + 2
            NX   = (-NUS(1)   +NU + 2)/2
            NY   = ( NUS(NNUS)-NU + 2)/2
C
C           TEMPORARY STORAGE OF VALUES
            R2111 = EIK(M,IR1+2)*EJL(M,IR2+1)*B3(M,NX+1,NY+1)
     &            + EIK(M,IR2+1)*EJL(M,IR1+2)*B4(M,NX+1,NY+1)
            R4121 = EIK(M,IR1+4)*EJL(M,IR2+1)*B3(M,NX+2,NY+1)
     &            + EIK(M,IR2+3)*EJL(M,IR1+2)*B4(M,NX+1,NY+2)
            R4322 = EIK(M,IR1+4)*EJL(M,IR2+3)*B3(M,NX+2,NY+2)
     &            + EIK(M,IR2+3)*EJL(M,IR1+4)*B4(M,NX+2,NY+2)
C
C           TRANSFER TO RK ARRAYS
            RKLLLL(M,4) = RKLLLL(M,4) + BK4(INU)*F51*R2111
            RKSLSL(M,4) = RKSLSL(M,4) + BK4(INU)*F73*R4121
            RKSSSS(M,4) = RKSSSS(M,4) + BK4(INU)*F92*R4322
          ENDDO
C
        ENDDO
C
C *** END CONDITIONAL STATEMENT OVER LQNA AND LQNB COMBINATIONS
      ENDIF
C
C*********************************************************************C
C     NORMALISE ACCUMULATED INTEGRALS (DIRECT AND EXCHANGE)           C
C*********************************************************************C
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
      
      GOTO 111
      IF(LQNA.NE.0.AND.LQNB.NE.0) THEN
        WRITE(*,*) RJLLLL(1,1),RJLLSS(1,1),RJSSLL(1,1),RJSSSS(1,1)
        WRITE(*,*) RJLLLL(1,2),RJLLSS(1,2),RJSSLL(1,2),RJSSSS(1,2)
        WRITE(*,*) RJLLLL(1,3),RJLLSS(1,3),RJSSLL(1,3),RJSSSS(1,3)
        WRITE(*,*) RJLLLL(1,4),RJLLSS(1,4),RJSSLL(1,4),RJSSSS(1,4)
        WRITE(*,*) RKLLLL(1,1),RKSLSL(1,1),RKSSSS(1,1)
        WRITE(*,*) RKLLLL(1,2),RKSLSL(1,2),RKSSSS(1,2)
        WRITE(*,*) RKLLLL(1,3),RKSLSL(1,3),RKSSSS(1,3)
        WRITE(*,*) RKLLLL(1,4),RKSLSL(1,4),RKSSSS(1,4)
        STOP
      ENDIF
111   CONTINUE      
C
C     END OF INTEGRAL BATCH GENERATION
C
      RETURN
      END



C
C >>>     LQNA  =  0 AND LQNB  =  0 (REQUIRED FOR ALL BLOCKS)
          RKLLSS(M,4) = RKLLSS(M,4) + DK(IV,4)*V4*F0G0*E0011*C7*B33
          RKSLSL(M,4) = RKSLSL(M,4) + HK(IV,4)*V4*F0G0*E1010*C7*B42
          RKSSLL(M,4) = RKSSLL(M,4) + FK(IV,4)*V4*F0G0*E1100*C7*B33         
C
C >>>     LQNA =/= 0                (NEED KQNA > 0 BLOCK)
          IF(LQNA.EQ.0) GOTO 203
          RKLL = V4*F0G0*E0011*C7*B33
          RKSL = V4*F0G0*E1010*C7*B42 - V2*F1G0*E0010*C5*B22
          RKSS = V4*F0G0*E1100*C7*B33 - V2*F1G0*E1000*C5*B31
     &         - V2*F1G0*E0100*C5*B13 + V1*F2G0*E0000*C3*B11
          RKLLSS(M,3) = RKLLSS(M,3) + DK(IV,3)*RKLL
          RKSLSL(M,3) = RKSLSL(M,3) + HK(IV,3)*RKSL
          RKSSLL(M,3) = RKSSLL(M,3) + FK(IV,3)*RKSS
203       CONTINUE
C
C >>>     LQNB =/= 0                (NEED KQNB > 0 BLOCK)
          IF(LQNB.EQ.0) GOTO 202
          RKLL = V4*F0G0*E0011*C7*B33 - V2*F0G1*E0010*C5*B31
     &         - V2*F0G1*E0001*C5*B13 + V1*F1G1*E0000*C3*B11
          RKSL = V4*F0G0*E1010*C7*B42 - V2*F0G1*E1000*C5*B22
          RKSS = V4*F0G0*E1100*C7*B33
          RKLLSS(M,2) = RKLLSS(M,2) + DK(IV,2)*RKLL
          RKSLSL(M,2) = RKSLSL(M,2) + HK(IV,2)*RKSL
          RKSSLL(M,2) = RKSSLL(M,2) + FK(IV,2)*RKSS
202       CONTINUE
C
C >>>     LQNA =/= 0 AND LQNB =/= 0 (NEED KQNA,KQNB > 0 BLOCK)
          IF(LQNA.EQ.0.OR.LQNB.EQ.0) GOTO 201
          RKLL = V4*F0G0*E0011*C7*B33 - V2*F0G1*E0010*C5*B31
     &         - V2*F0G1*E0001*C5*B13 + V1*F1G1*E0000*C3*B11
          RKSL = V4*F0G0*E1010*C7*B42 - V2*F1G0*E0010*C5*B22
     &         - V2*F0G1*E1000*C5*B22 + V1*F1G1*E0000*C3*B02
          RKSS = V4*F0G0*E1100*C7*B33 - V2*F1G0*E1000*C5*B31
     &         - V2*F1G0*E0100*C5*B13 + V1*F2G0*E0000*C3*B11
          RKLLSS(M,1) = RKLLSS(M,1) + DK(IV,1)*RKLL
          RKSLSL(M,1) = RKSLSL(M,1) + HK(IV,1)*RKSL
          RKSSLL(M,1) = RKSSLL(M,1) + FK(IV,1)*RKSS
201       CONTINUE


