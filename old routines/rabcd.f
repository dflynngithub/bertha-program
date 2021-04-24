      SUBROUTINE RABCD(KQN,LQN,MAXM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C             RRRRRRR     AA    BBBBBBB   CCCCCC  DDDDDDD              C
C             RR    RR   AAAA   BB    BB CC    CC DD    DD             C
C             RR    RR  AA  AA  BB    BB CC       DD    DD             C
C             RR    RR AA    AA BBBBBBB  CC       DD    DD             C
C             RRRRRRR  AAAAAAAA BB    BB CC       DD    DD             C
C             RR    RR AA    AA BB    BB CC    CC DD    DD             C
C             RR    RR AA    AA BBBBBBB   CCCCCC  DDDDDDD              C
C                                                                      C
C -------------------------------------------------------------------- C
C  RABCD EVALUATES THE EFFECTIVE INTERACTION STRENGTHS OF THE COULOMB  C
C  AND BREIT INTERACTIONS IN A G-SPINOR BASIS SET. THE ALGORITHM IS    C
C  SEGMENTED TO TAKE INTO ACCOUNT THE SIGNS OF KQN IN THE T=S CASE.    C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS,MAB=10,MNU=10)
C
      CHARACTER*4 HMLTN
C
      COMMON/GAMA/GAMLOG(50),GAMHLF(50)
      COMMON/PRMS/CV,HMLTN,ITREE,IMOL,INEW,IEQS,IPAR,ICOR,ILEV
C
      COMMON/XCAL1/XLLLL(MB2,MNU),XLLSS(MB2,MNU),
     &             XSSLL(MB2,MNU),XSSSS(MB2,MNU)
      COMMON/XCAL2/EK(MB2),EKL0(MB2),
     &             EL(MB2),EKL1(MB2,MAB),EKL2(MB2,MAB), 
     &             RNLLKL(MB2),RNSSKL(MB2),
     &             RNSLKL(MB2),RNLSKL(MB2),
     &             B1(MB2,MAB,MAB),B2(MB2,MAB,MAB)
      COMMON/XCAL3/AI,AJ,AK,AL,IATEST(16),
     &             NUS(MNU),NUCMAX,NUCMIN,NUNUM
      COMMON/XCAL4/EIJ1(MAB),EIJ2(MAB),RNLLIJ,RNSSIJ,RNSLIJ,RNLSIJ,EI,EJ
C
      DIMENSION XJ(MB2,2),XROOT(MB2),BETA(MB2),BTA1(MB2),
     &          BIN(MB2),TRM(MB2),KQN(4),LQN(4),IAA(2),IBB(2)  
C
      DIMENSION KSET(16)
C
C     GAUSSIAN EXPONENT FOR (IBAS,JBAS) OVERLAP
      TIJ0 = EI+EJ
C
C     INITIALIZE THE ARRAY XJ(M,2) FOR INCOMPLETE BETA FUNCTION ARGS.
      DO M=1,MAXM
        TIJKL   = TIJ0+EKL0(M)
        XJ(M,1) = TIJ0/TIJKL
        XJ(M,2) = EKL0(M)/TIJKL
      ENDDO
C
C**********************************************************************C
C     BEGIN INLINE BETA FUNCTION CODE                                  C
C**********************************************************************C
C
C     GENERATE A MASTER TABLE OF INCOMPLETE BETA FUNCTIONS FOR THIS
C     BATCH OF RADIAL INTEGRALS. (THE SAME BATCH CAN BE USED
C     FOR BOTH THE COULOMB AND BREIT INTEGRALS)
C
      NVALS = ((NUS(NUNUM)-NUS(1))/2)+3
C
C     LOOP OVER EXPANSION TERMINALS FOR FIRST PAIR
      DO I1=1,NVALS
        NT1 = 2*(I1-1) + NUS(1)
        IAA(1) = LQN(1)+LQN(2)+NT1+1
        IAA(2) = LQN(3)+LQN(4)+NT1+1
C
C       LOOP OVER EXPANSION TERMINALS FOR SECOND PAIR
        DO I2=1,NVALS
          NT2 = 2*(I2-1) - NUS(NUNUM)
          IBB(1) = LQN(3)+LQN(4)+NT2
          IBB(2) = LQN(1)+LQN(2)+NT2
C
C         LOOP OVER BETA INTEGRAL COMBINATIONS
          DO IBETA=1,2
            IA =(IAA(IBETA)-1)/2
            IB = IBB(IBETA)   /2
C
C ***       BEGIN CONDITIONAL STATEMENT OVER IB VALUES
C >>>       CASE 1: IB > 1
            IF(IB.GT.1) THEN
              X  = DFLOAT(IA) + 0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BTA1(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
              RA  = X
              RB  = DFLOAT(1-IB)
              RC  = X + 1.0D0
              RD  = 1.0D0
              RCD = RC*RD
              FCT = RA*RB/RCD
              DO M=1,MAXM
                TRM(M) = FCT*XJ(M,IBETA)
                BIN(M) = 1.0D0 + TRM(M)
              ENDDO
              RA = RA + 1.0D0
              RB = RB + 1.0D0
              RC = RC + 1.0D0
              RD = RD + 1.0D0
              DO IT=2,IB-1
                RCD = RC*RD
                FCT = RA*RB/RCD
                DO M=1,MAXM
                  TRM(M) = FCT*TRM(M)*XJ(M,IBETA)
                  BIN(M) = BIN(M) + TRM(M)
                ENDDO
                RA = RA + 1.0D0
                RB = RB + 1.0D0
                RC = RC + 1.0D0
                RD = RD + 1.0D0
              ENDDO
              DO M=1,MAXM
                BETA(M) = BTA1(M)*BIN(M)
              ENDDO
C
C >>>       CASE 2: IB = 1
C
            ELSEIF(IB.EQ.1) THEN
              X  = DFLOAT(IA) + 0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BETA(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
C
C >>>       CASE 3: IB = 0
            ELSEIF(IB.EQ.0) THEN
              DO M=1,MAXM
                XROOT(M) = DSQRT(XJ(M,IBETA))
                DEN      = (1.0D0-XROOT(M))
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
                    BIN(M) = BIN(M) + X*TRM(M)
                    TRM(M) = TRM(M)*XJ(M,IBETA)
                  ENDDO
                ENDDO
                DO M=1,MAXM
                  BETA(M) = BTA1(M) - 2.0D0*XROOT(M)*BIN(M)
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
C ***       END CONDITIONAL STATEMENT OVER IB VALUES
            ENDIF
C
            IF(IBETA.EQ.1) THEN
              DO M=1,MAXM
                B1(M,I1,I2) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXM
                B2(M,I1,I2) = BETA(M)
              ENDDO
            ENDIF

          ENDDO

        ENDDO
      ENDDO
C
C**********************************************************************C
C     END INLINE BETA FUNCTION CODE                                    C
C**********************************************************************C
C
C     PREPARE VALUES FOR UPCOMING CALCULATIONS
      IF(HMLTN.EQ.'NORL') THEN
C     NON-RELATIVISTIC HAMILTONIAN
C
        C5 = GAMHLF(LQN(1)+LQN(2)+LQN(3)+LQN(4)+5)
C
      ELSE
C     RELATIVISTIC HAMILTONIAN
C
        C1 = 0.25D0*GAMHLF(LQN(1)+LQN(2)+LQN(3)+LQN(4)+1)
        C3 = 0.25D0*GAMHLF(LQN(1)+LQN(2)+LQN(3)+LQN(4)+3)
        C5 = 0.25D0*GAMHLF(LQN(1)+LQN(2)+LQN(3)+LQN(4)+5)
        C7 = 0.25D0*GAMHLF(LQN(1)+LQN(2)+LQN(3)+LQN(4)+7)
        C9 = 0.25D0*GAMHLF(LQN(1)+LQN(2)+LQN(3)+LQN(4)+9)     
C
        V1 = 1.0D0
        V2 = 2.0D0
        V4 = 4.0D0
        V8 = 8.0D0
        VS = 1.6D1
C
        KTI = LQN(1)+KQN(1)+1
        KTJ = LQN(2)+KQN(2)+1
        KTK = LQN(3)+KQN(3)+1
        KTL = LQN(4)+KQN(4)+1
C
        TI = DFLOAT(KTI)
        TJ = DFLOAT(KTJ)
        TK = DFLOAT(KTK)
        TL = DFLOAT(KTL)
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
C       CASES FOR INTEGRAL INCLUSION
        DO N=1,16
          KSET(N) = 0
        ENDDO
C
        IF(KQN(4).GT.0)                 KSET( 2) = 1
        IF(KQN(3).GT.0)                 KSET( 3) = 1
        IF(KQN(3).GT.0.AND.KQN(4).GT.0) KSET( 4) = 1
        IF(KQN(2).GT.0)                 KSET( 5) = 1
        IF(KQN(2).GT.0.AND.KQN(4).GT.0) KSET( 6) = 1
        IF(KQN(2).GT.0.AND.KQN(4).GT.0) KSET( 7) = 1
        IF(KQN(2).GT.0.AND.KQN(4).GT.0) KSET( 8) = 1
        IF(KQN(2).GT.0.AND.KQN(4).GT.0) KSET( 9) = 1
        IF(KQN(2).GT.0.AND.KQN(4).GT.0) KSET(10) = 1
        IF(KQN(2).GT.0.AND.KQN(4).GT.0) KSET(11) = 1
        IF(KQN(2).GT.0.AND.KQN(4).GT.0) KSET(12) = 1
        IF(KQN(2).GT.0.AND.KQN(4).GT.0) KSET(13) = 1
        IF(KQN(2).GT.0.AND.KQN(4).GT.0) KSET(14) = 1
        IF(KQN(2).GT.0.AND.KQN(4).GT.0) KSET(15) = 1
        IF(KQN(2).GT.0.AND.KQN(4).GT.0) KSET(16) = 1
C
      ENDIF
C
C**********************************************************************C
C     A (KQN(1),KQN(2),KQN(3),KQN(4)) COMBINATION HAS BETWEEN 1 AND 16 C
C     CONTRIBUTING TERMS DUE TO VANISHING PARTS OF SMALL COMPONENTS.   C
C**********************************************************************C
C
C     INITIATE LOOP OVER K,L BASIS FUNCTIONS
      DO M=1,MAXM
C
C       MORE VALUE PREPARATION
        IF(HMLTN.NE.'NORL') THEN
C         NON-RELATIVISTIC HAMILTONIAN
C
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
        ENDIF
C
C       LOOP OVER THE TENSOR ORDERS OF THE COULOMB INTERACTION
        DO IV=1,NUNUM
C
C         IMPORT CURRENT NU VALUE FROM THE NUS ARRAY
          NUVAL = NUS(IV)
C
C         INDEX OFFSETS FOR TENSOR ORDER
          IQ1 = (-NUS(    1)+NUVAL)/2
          IQ2 = (+NUS(NUNUM)-NUVAL)/2
C
C**********************************************************************C
C         SPECIAL CASE: NON-RELATIVISTIC HAMILTONIAN                   C
C**********************************************************************C
C
          IF(HMLTN.EQ.'NORL') THEN
C         NON-RELATIVISTIC HAMILTONIAN
C
C           BETA INTEGRALS WITH POWERS OF EXPONENTIALS
            B11 = EIJ1(IQ1+2)*EKL1(M,IQ2+2)*B1(M,IQ1+2,IQ2+2)
     &          + EIJ2(IQ2+2)*EKL2(M,IQ1+2)*B2(M,IQ1+2,IQ2+2)
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XLLLL(M,IV) = C5*B11
C
C           SKIP PAST THE RELATIVISTIC STAGES
            GOTO 999
C
          ENDIF
C
C**********************************************************************C
C         LARGER CASE: RELATIVISTIC HAMILTONIAN                        C
C**********************************************************************C
C
C         BETA INTEGRALS WITH POWERS OF EXPONENTIALS
          B00 = EIJ1(IQ1+1)*EKL1(M,IQ2+1)*B1(M,IQ1+1,IQ2+1)
     &        + EIJ2(IQ2+1)*EKL2(M,IQ1+1)*B2(M,IQ1+1,IQ2+1)
          B01 = EIJ1(IQ1+1)*EKL1(M,IQ2+2)*B1(M,IQ1+1,IQ2+2) 
     &        + EIJ2(IQ2+1)*EKL2(M,IQ1+2)*B2(M,IQ1+2,IQ2+1)
          B02 = EIJ1(IQ1+1)*EKL1(M,IQ2+3)*B1(M,IQ1+1,IQ2+3) 
     &        + EIJ2(IQ2+1)*EKL2(M,IQ1+3)*B2(M,IQ1+3,IQ2+1)
          B10 = EIJ1(IQ1+2)*EKL1(M,IQ2+1)*B1(M,IQ1+2,IQ2+1) 
     &        + EIJ2(IQ2+2)*EKL2(M,IQ1+1)*B2(M,IQ1+1,IQ2+2)
          B11 = EIJ1(IQ1+2)*EKL1(M,IQ2+2)*B1(M,IQ1+2,IQ2+2) 
     &        + EIJ2(IQ2+2)*EKL2(M,IQ1+2)*B2(M,IQ1+2,IQ2+2)
          B12 = EIJ1(IQ1+2)*EKL1(M,IQ2+3)*B1(M,IQ1+2,IQ2+3) 
     &        + EIJ2(IQ2+2)*EKL2(M,IQ1+3)*B2(M,IQ1+3,IQ2+2)
          B20 = EIJ1(IQ1+3)*EKL1(M,IQ2+1)*B1(M,IQ1+3,IQ2+1) 
     &        + EIJ2(IQ2+3)*EKL2(M,IQ1+1)*B2(M,IQ1+1,IQ2+3)
          B21 = EIJ1(IQ1+3)*EKL1(M,IQ2+2)*B1(M,IQ1+3,IQ2+2) 
     &        + EIJ2(IQ2+3)*EKL2(M,IQ1+2)*B2(M,IQ1+2,IQ2+3)
          B22 = EIJ1(IQ1+3)*EKL1(M,IQ2+3)*B1(M,IQ1+3,IQ2+3) 
     &        + EIJ2(IQ2+3)*EKL2(M,IQ1+3)*B2(M,IQ1+3,IQ2+3)
C
C
C         STAGE 1 (APPLIES TO ALL KQN POSSIBILITIES):
C
C         EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
          XLLLL(M,IV) = V1*T0000*E0000*C5*B11
          XLLSS(M,IV) = V4*T0000*E0011*C7*B12
          XSSLL(M,IV) = V4*T0000*E1100*C7*B21
          XSSSS(M,IV) = VS*T0000*E1111*C9*B22
C
C         STAGE 2: KQN(4)>0
          IF(KQN(4).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XLLSS(M,IV) = XLLSS(M,IV) - V2*T0001*E0010*C5*B11
            XSSSS(M,IV) = XSSSS(M,IV) - V8*T0001*E1110*C7*B21
C
          ENDIF
C
C         STAGE 3: KQN(3)>0
          IF(KQN(3).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XLLSS(M,IV) = XLLSS(M,IV) - V2*T0010*E0001*C5*B11
            XSSSS(M,IV) = XSSSS(M,IV) - V8*T0010*E1101*C7*B21
C
          ENDIF
C
C         STAGE 4: KQN(3)>0 AND KQN(4)>0
          IF(KQN(3).GT.0.AND.KQN(4).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XLLSS(M,IV) = XLLSS(M,IV) + V1*T0011*E0000*C3*B10
            XSSSS(M,IV) = XSSSS(M,IV) + V4*T0011*E1100*C5*B20
C
          ENDIF
C
C         STAGE 5: KQN(2)>0
          IF(KQN(2).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XSSLL(M,IV) = XSSLL(M,IV) - V2*T0100*E1000*C5*B11
            XSSSS(M,IV) = XSSSS(M,IV) - V8*T0100*E1011*C7*B12
C
          ENDIF
C
C         STAGE 6: KQN(2)>0 AND KQN(4)>0
          IF(KQN(2).GT.0.AND.KQN(4).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XSSSS(M,IV) = XSSSS(M,IV)+ V4*T0101*E1010*C5*B11
C
          ENDIF
C
C         STAGE 7: KQN(2)>0 AND KQN(3)>0
          IF(KQN(2).GT.0.AND.KQN(3).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XSSSS(M,IV) = XSSSS(M,IV)+ V4*T0110*E1001*C5*B11
C
          ENDIF
C
C         STAGE 8: KQN(2)>0 AND KQN(3)>0 AND KQN(4)>0
          IF(KQN(2).GT.0.AND.KQN(3).GT.0.AND.KQN(4).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XSSSS(M,IV) = XSSSS(M,IV) - V2*T0111*E1000*C3*B10
C
          ENDIF
C
C         STAGE 9: KQN(1)>0
          IF(KQN(1).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XSSLL(M,IV) = XSSLL(M,IV) - V2*T1000*E0100*C5*B11
            XSSSS(M,IV) = XSSSS(M,IV) - V8*T1000*E0111*C7*B12
C
          ENDIF
C
C         STAGE 10: KQN(1)>0 AND KQN(4)>0
          IF(KQN(1).GT.0.AND.KQN(4).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XSSSS(M,IV) = XSSSS(M,IV) + V4*T1001*E0110*C5*B11
C
          ENDIF
C
C         STAGE 11: KQN(1)>0 AND KQN(3)>0
          IF(KQN(1).GT.0.AND.KQN(3).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XSSSS(M,IV) = XSSSS(M,IV) + V4*T1010*E0101*C5*B11
C
          ENDIF
C
C         STAGE 12: KQN(1)>0 AND KQN(3)>0 AND KQN(4)>0
          IF(KQN(1).GT.0.AND.KQN(3).GT.0.AND.KQN(4).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XSSSS(M,IV) = XSSSS(M,IV)- V2*T1011*E0100*C3*B10
C
          ENDIF
C
C         STAGE 13: KQN(1)>0 AND KQN(2)>0
          IF(KQN(1).GT.0.AND.KQN(2).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XSSLL(M,IV) = XSSLL(M,IV)+ V1*T1100*E0000*C3*B01
            XSSSS(M,IV) = XSSSS(M,IV)+ V4*T1100*E0011*C5*B02
C
          ENDIF
C
C         STAGE 14: KQN(1)>0 AND KQN(2)>0 AND KQN(4)>0
          IF(KQN(1).GT.0.AND.KQN(2).GT.0.AND.KQN(4).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XSSSS(M,IV) = XSSSS(M,IV)- V2*T1101*E0010*C3*B01
C
          ENDIF
C
C         STAGE 15: KQN(1)>0 AND KQN(2)>0 AND KQN(3)>0
          IF(KQN(1).GT.0.AND.KQN(2).GT.0.AND.KQN(3).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XSSSS(M,IV) = XSSSS(M,IV) - V2*T1110*E0001*C3*B01
C
          ENDIF
C
C         STAGE 16: KQN(1)>0 AND KQN(2)>0 AND KQN(3)>0 AND KQN(4)>0
          IF(KQN(1).GT.0.AND.KQN(2).GT.0.AND.KQN(3).GT.0.AND. 
     &                                                KQN(4).GT.0) THEN
C
C           EFFECTIVE INTERACTION STRENGTH FOR THIS TENSOR ORDER
            XSSSS(M,IV) = XSSSS(M,IV)+ V1*T1111*E0000*C1*B00
C
          ENDIF
C
C         SKIP POINT FOR NON-RELATIVISTIC HAMILTONIANS
999       CONTINUE
C
C       END LOOP OVER TENSOR ORDERS
        ENDDO
C
C     END LOOP OVER K,L BASIS FUNCTIONS
      ENDDO
C
C**********************************************************************C
C     NORMALISE ACCUMULATED INTEGRALS (DIRECT AND EXCHANGE)            C
C**********************************************************************C
C
      IF(HMLTN.EQ.'NORL') THEN
C     NON-RELATIVISTIC HAMILTONIAN
C
        DO M=1,MAXM
          DO IV=1,NUNUM      
            XLLLL(M,IV) = RNLLIJ*RNLLKL(M)*XLLLL(M,IV)
            XLLSS(M,IV) = 0.0D0
            XSSLL(M,IV) = 0.0D0
            XSSSS(M,IV) = 0.0D0
          ENDDO
        ENDDO
C
      ELSE
C     RELATIVISTIC HAMILTONIAN
C
        DO M=1,MAXM
          DO IV=1,NUNUM      
            XLLLL(M,IV) = RNLLIJ*RNLLKL(M)*XLLLL(M,IV)
            XLLSS(M,IV) = RNLLIJ*RNSSKL(M)*XLLSS(M,IV)
            XSSLL(M,IV) = RNSSIJ*RNLLKL(M)*XSSLL(M,IV)
            XSSSS(M,IV) = RNSSIJ*RNSSKL(M)*XSSSS(M,IV)
          ENDDO
        ENDDO
C
      ENDIF
C
      RETURN
      END
