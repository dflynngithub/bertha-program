      SUBROUTINE RNORMF(EXPT,LQN,NFUNS,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       RRRRRRR  NN    NN  OOOOOO  RRRRRRR  MM       MM FFFFFFFF       C      
C       RR    RR NNN   NN OO    OO RR    RR MMM     MMM FF             C 
C       RR    RR NNNN  NN OO    OO RR    RR MMMM   MMMM FF             C
C       RR    RR NN NN NN OO    OO RR    RR MM MM MM MM FFFFFF         C
C       RRRRRRR  NN  NNNN OO    OO RRRRRRR  MM  MMM  MM FF             C
C       RR    RR NN   NNN OO    OO RR    RR MM   M   MM FF             C
C       RR    RR NN    NN  OOOOOO  RR    RR MM       MM FF             C
C                                                                      C
C -------------------------------------------------------------------- C
C     RNORMF EVALUATES NORMALISATION CONSTANTS OF ALL VARIETIES.       C
C**********************************************************************C
      PARAMETER(MBS=26,MB2=MBS*MBS)
C
      DIMENSION EXPT(MBS,4),NFUNS(4),LQN(4)
      DIMENSION RNAL(MBS),RNAS(MBS),RNBL(MBS),RNBS(MBS)
      DIMENSION EXPA(MB2),EXPB(MB2),EXPAB(MB2)
C
      COMMON/GMFN/GAMMAL(100),GAMMAF(100)
      COMMON/RNRM/RNLL(MB2),RNSL(MB2),RNLS(MB2),RNSS(MB2)
C
      DATA TWOLOG/6.93147180559945309D-1/
C
      LA    = LQN(I1)
      LB    = LQN(I2)
      NFUNA = NFUNS(I1)
      NFUNB = NFUNS(I2)
      MAXM  = NFUNA*NFUNB
      RLA   = DFLOAT(LA)
      RLB   = DFLOAT(LB)
      GA1   = TWOLOG - GAMMAL(2*LA + 3)
      GA2   = TWOLOG - GAMMAL(2*LA + 5)
      GB1   = TWOLOG - GAMMAL(2*LB + 3)
      GB2   = TWOLOG - GAMMAL(2*LB + 5)
      RA1   = RLA + 1.5D0
      RA2   = RLA + 0.5D0
      RB1   = RLB + 1.5D0
      RB2   = RLB + 0.5D0
      DO IBAS=1,NFUNA
        ELOG       = DLOG(2.0D0*EXPT(IBAS,I1))
        RNAL(IBAS) = DEXP(0.5D0*(GA1+RA1*ELOG))
        RNAS(IBAS) = DEXP(0.5D0*(GA2+RA2*ELOG))
      ENDDO
      DO JBAS=1,NFUNB
        ELOG       = DLOG(2.0D0*EXPT(JBAS,I2))
        RNBL(JBAS) = DEXP(0.5D0*(GB1+RB1*ELOG))
        RNBS(JBAS) = DEXP(0.5D0*(GB2+RB2*ELOG))
      ENDDO
C
C     RNLL(M) ARE THE LL NORMALISATION CONSTANTS
C     RNSL(M) ARE THE SL NORMALISATION CONSTANTS
C     RNSS(M) ARE THE SS NORMALISATION CONSTANTS
C
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M       = M + 1
          RNLL(M) = RNAL(IBAS)*RNBL(JBAS)
          RNSL(M) = RNAS(IBAS)*RNBL(JBAS)
          RNLS(M) = RNAL(IBAS)*RNBS(JBAS)
          RNSS(M) = RNAS(IBAS)*RNBS(JBAS)
          EXPA(M) = EXPT(IBAS,I1)
          EXPB(M) = EXPT(JBAS,I2)
        ENDDO
      ENDDO
C
      MAXM = M
      DO M=1,MAXM
        EXPAB(M) = EXPA(M)*EXPB(M)
      ENDDO
C
      RETURN
      END
