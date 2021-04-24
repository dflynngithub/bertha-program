      SUBROUTINE RKBRT1(RJLSSL,RJSLLS,RJSLSL,RJLSLS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C            RRRRRRR  KK    KK BBBBBBB  RRRRRRR TTTTTTTT 11            C
C            RR    RR KK   KK  BB    BB RR    RR   TT   111            C
C            RR    RR KK  KK   BB    BB RR    RR   TT    11            C
C            RR    RR KKKKK    BBBBBBB  RR    RR   TT    11            C
C            RRRRRRR  KK  KK   BB    BB RRRRRRR    TT    11            C
C            RR    RR KK   KK  BB    BB RR    RR   TT    11            C
C            RR    RR KK    KK BBBBBBB  RR    RR   TT   1111           C
C                                                                      C
C -------------------------------------------------------------------- C
C  RKBRT1 EVALUATES A BATCH OF GENERAL (OPEN-SHELL) ONE-CENTRE RADIAL  C
C  INTEGRALS OVER THE BREIT INTERACTION, SEPARATING RESULTS BY THE     C
C  ALLOWED TENSOR ORDERS RK(ABCD).                                     C
C -------------------------------------------------------------------- C
C  DFNOTE: NOT WORKING YET...                                          C
C**********************************************************************C
      INCLUDE 'param.h'
C
      DIMENSION XJ(MB2,2),IAA(2),IBB(2)
      DIMENSION BETA(MB2),BTA1(MB2),XROOT(MB2),BIN(MB2),TRM(MB2)
      DIMENSION BDU(MB2,-MAB:MAB,-MAB:MAB),BDL(MB2,-MAB:MAB,-MAB:MAB)
C
      DIMENSION RJLSSL(MB2,MNU,2),RJSLLS(MB2,MNU,2),
     &          RJSLSL(MB2,MNU,2),RJLSLS(MB2,MNU,2)
C
      COMMON/B0IJ/EIJ(-MAB:MAB),RNIJ(4),EI,EJ,IBAS,JBAS
      COMMON/B0KL/EKL(MB2,-MAB:MAB),RNKL(MB2,4),EK(MB2),EL(MB2)
      COMMON/B1QN/EXL(MBS,4),MQN(4),KQN(4),LQN(4),NBAS(4),MAXM
      COMMON/RCFF/T0000,T1000,T0100,T0010,T0001,T1100,T1010,T1001,
     &            T0110,T0101,T0011,T1110,T1101,T1011,T0111,T1111,
     &            C1,C3,C5,C7,C9,V1,V2,V4,V8,VS
      COMMON/XNUS/NUS(MNU),NUI,NUF,NUNUM
C
C**********************************************************************C
C     PREPARATION OF EXPONENT POWERS AND BETA FUNCTION ARGUMENTS.      C
C**********************************************************************C
C
C     GAUSSIAN EXPONENT FOR (IBAS,JBAS) OVERLAP
      TIJ0 = EI+EJ
C
C     INITIALISE THE ARRAY XJ(M,2) FOR INCOMPLETE BETA FUNCTION ARGS.
      DO M=1,MAXM
        TKL0    = EK(M)+EL(M)
        TIJKL   = TIJ0+TKL0
        XJ(M,1) = TIJ0/TIJKL
        XJ(M,2) = TKL0/TIJKL
      ENDDO
C
C**********************************************************************C
C     INLINE BETA FUNCTION CODE FOR DIRECT TERMS                       C
C**********************************************************************C
C
C     NUMBER OF LEVELS NEEDED TO ACCOUNT FOR ALL BETA INTEGRALS
      NVALS = (NUF-NUI)/2+2
C
C     LOOP OVER EXPANSION TERMINALS FOR FIRST PAIR
      DO NX=1,NVALS
        IAA(1) = LQN(1)+LQN(2)+NUI+2*NX
        IAA(2) = LQN(3)+LQN(4)+NUI+2*NX
C
C       LOOP OVER EXPANSION TERMINALS FOR SECOND PAIR
        DO NY=1,NVALS
          IBB(1) = LQN(3)+LQN(4)-NUF+2*NY-1
          IBB(2) = LQN(1)+LQN(2)-NUF+2*NY-1
C
C         LOOP OVER BETA INTEGRAL COMBINATIONS
          DO IBETA=1,2
            IA =(IAA(IBETA)-1)/2
            IB = IBB(IBETA)   /2
C
C           BEGIN CONDITIONAL STATEMENT OVER IB VALUES
C           CASE 1: IB > 1
            IF(IB.GT.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BTA1(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
              RA  = X
              RB  = DFLOAT(1-IB)
              RC  = X+1.0D0
              RD  = 1.0D0
              RCD = RC*RD
              FCT = RA*RB/RCD
              DO M=1,MAXM
                TRM(M) = FCT*XJ(M,IBETA)
                BIN(M) = 1.0D0+TRM(M)
              ENDDO
              RA = RA+1.0D0
              RB = RB+1.0D0
              RC = RC+1.0D0
              RD = RD+1.0D0
              DO IT=2,IB-1
                RCD = RC*RD
                FCT = RA*RB/RCD
                DO M=1,MAXM
                  TRM(M) = FCT*TRM(M)*XJ(M,IBETA)
                  BIN(M) = BIN(M)+TRM(M)
                ENDDO
                RA = RA+1.0D0
                RB = RB+1.0D0
                RC = RC+1.0D0
                RD = RD+1.0D0
              ENDDO
              DO M=1,MAXM
                BETA(M) = BTA1(M)*BIN(M)
              ENDDO
C
C           CASE 2: IB = 1
            ELSEIF(IB.EQ.1) THEN
              X  = DFLOAT(IA)+0.5D0
              IX = 2*IA+1
              DO M=1,MAXM
                BETA(M) = (DSQRT(XJ(M,IBETA))**IX)/X
              ENDDO
C
C           CASE 3: IB = 0
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
                    BIN(M) = BIN(M)+X*TRM(M)
                    TRM(M) = TRM(M)*XJ(M,IBETA)
                  ENDDO
                ENDDO
                DO M=1,MAXM
                  BETA(M) = BTA1(M)-2.0D0*XROOT(M)*BIN(M)
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
C           END CONDITIONAL STATEMENT OVER IB VALUES
            ENDIF
C
            IF(IBETA.EQ.1) THEN
              DO M=1,MAXM
                BDU(M, NUI+2*NX  ,-NUF+2*NY-1) = BETA(M)
              ENDDO
            ELSE
              DO M=1,MAXM
                BDL(M, NUI+2*NX  ,-NUF+2*NY-1) = BETA(M)
              ENDDO
            ENDIF

          ENDDO

        ENDDO
      ENDDO
C
C**********************************************************************C
C     RADIAL INTEGRALS OVER SPINORS (SEPARATED BY TENSOR ORDER)        C
C**********************************************************************C
C
C     EMPTY COUNTER ARRAYS FOR DIRECT AND EXCHANGE INTEGRALS
      DO M=1,MAXM
        DO LTEN=1,NUNUM
          DO ILU=1,2
            RJLSLS(M,LTEN,ILU) = 0.0D0
            RJSLSL(M,LTEN,ILU) = 0.0D0
            RJSLLS(M,LTEN,ILU) = 0.0D0
            RJLSSL(M,LTEN,ILU) = 0.0D0
          ENDDO
        ENDDO
      ENDDO
C
C     VALUE PREPARATION
      E0000 = 1.0D0
      E1000 = EI
      E0100 = EJ
      E1100 = EI*EJ
C
C     LOOP OVER K,L BASIS FUNCTIONS
      DO M=1,MAXM
C
C       MORE VALUE PREPARATION
        E0010 = EK(M)
        E0001 = EL(M)
        E1010 = EI*EK(M)
        E1001 = EI*EL(M)
        E0110 = EJ*EK(M)
        E0101 = EJ*EL(M)
        E0011 = EK(M)*EL(M)
C
C       LOOP OVER THE TENSOR ORDERS OF THE BREIT INTERACTION
C       (U -> Z, L -> Z')
        DO LTEN=1,NUNUM
C
C         IMPORT NU VALUE FROM ARRAY
          NU = NUS(LTEN)
C
C         TEMPORARY STORAGE OF RAW RJ(LTEN,M)
          B21L = EIJ(-NU+1)*EKL(M, NU+2)*BDL(M, NU+2,-NU+1)
          B21U = EIJ( NU+2)*EKL(M,-NU+1)*BDU(M, NU+2,-NU+1)
          B23L = EIJ(-NU+3)*EKL(M, NU+2)*BDL(M, NU+2,-NU+3)
          B41U = EIJ( NU+4)*EKL(M,-NU+1)*BDU(M, NU+4,-NU+1)
          B41L = EIJ(-NU+1)*EKL(M, NU+4)*BDL(M, NU+4,-NU+1)
          B23U = EIJ( NU+2)*EKL(M,-NU+3)*BDU(M, NU+2,-NU+3)
          B43L = EIJ(-NU+3)*EKL(M, NU+4)*BDL(M, NU+4,-NU+3)
          B43U = EIJ( NU+4)*EKL(M,-NU+3)*BDU(M, NU+4,-NU+3)
C
C         EFFECTIVE INTERACTION STRENGTH RJLSSL(M,LTEN)
          RJLSSL(M,LTEN,1)
     &         = V4*T0000*E0110*C7*B43L - V2*T0010*E0100*C5*B23L
     &         - V2*T0100*E0010*C5*B41L + V1*T0110*E0000*C3*B21L
C
          RJLSSL(M,LTEN,2) 
     &         = V4*T0000*E0110*C7*B43U - V2*T0010*E0100*C5*B41U
     &         - V2*T0100*E0010*C5*B23U + V1*T0110*E0000*C3*B21U
C
C         EFFECTIVE INTERACTION STRENGTH RJSLLS(M,LTEN)
          RJSLLS(M,LTEN,1) 
     &         = V4*T0000*E1001*C7*B43L - V2*T0001*E1000*C5*B23L
     &         - V2*T1000*E0001*C5*B41L + V1*T1001*E0000*C3*B21L
C
          RJSLLS(M,LTEN,2) 
     &         = V4*T0000*E1001*C7*B43U - V2*T0001*E1000*C5*B41U
     &         - V2*T1000*E0001*C5*B23U + V1*T1001*E0000*C3*B21U
C
C         EFFECTIVE INTERACTION STRENGTH RJSLSL(M,LTEN)
          RJSLSL(M,LTEN,1)
     &         = V4*T0000*E1010*C7*B43L - V2*T0010*E1000*C5*B23L
     &         - V2*T1000*E0010*C5*B41L + V1*T1010*E0000*C3*B21L
C
          RJSLSL(M,LTEN,2)
     &         = V4*T0000*E1010*C7*B43U - V2*T0010*E1000*C5*B41U
     &         - V2*T1000*E0010*C5*B23U + V1*T1010*E0000*C3*B21U
C
C         EFFECTIVE INTERACTION STRENGTH RJLSLS(M,LTEN)
          RJLSLS(M,LTEN,1) 
     &         = V4*T0000*E0101*C7*B43L - V2*T0001*E0100*C5*B23L
     &         - V2*T0100*E0001*C5*B41L + V1*T0101*E0000*C3*B21L
C
          RJLSLS(M,LTEN,2) 
     &         = V4*T0000*E0101*C7*B43U - V2*T0001*E0100*C5*B41U
     &         - V2*T0100*E0001*C5*B23U + V1*T0101*E0000*C3*B21U
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
C     APPLY NORMALISATION FACTORS                                      C
C**********************************************************************C
C
      DO M=1,MAXM
        RNLSLS = RNIJ(4)*RNKL(M,4)
        RNSLSL = RNIJ(2)*RNKL(M,2)
        RNSLLS = RNIJ(2)*RNKL(M,4)
        RNLSSL = RNIJ(4)*RNKL(M,2)
        DO LTEN=1,NUNUM
          DO ILU=1,2
            RJLSLS(M,LTEN,ILU) = RNLSLS*RJLSLS(M,LTEN,ILU)
            RJSLSL(M,LTEN,ILU) = RNSLSL*RJSLSL(M,LTEN,ILU)
            RJSLLS(M,LTEN,ILU) = RNSLLS*RJSLLS(M,LTEN,ILU)
            RJLSSL(M,LTEN,ILU) = RNLSSL*RJLSSL(M,LTEN,ILU)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END

