C
C
      SUBROUTINE RNLL(RN,EXL,LQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 RRRRRRR  NN    NN LL       LL                        C
C                 RR    RR NNN   NN LL       LL                        C
C                 RR    RR NNNN  NN LL       LL                        C
C                 RR    RR NN NN NN LL       LL                        C
C                 RRRRRRR  NN  NNNN LL       LL                        C
C                 RR    RR NN   NNN LL       LL                        C
C                 RR    RR NN    NN LLLLLLLL LLLLLLLL                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  RNLL GENERATES THE LARGE-LARGE SGTF NORMALISATION CONSTANTS.        C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ LQN  - PAIR OF ORBITAL QUANTUM NUMBERS.                           C
C  ▶ NBAS - PAIR OF PARAMETER LIST LENGTHS.                            C
C  OUTPUT:                                                             C
C  ▶ RN   - BATCH OF NORMALISATION CONSTANTS.                          C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RN(MBS,2),EXL(MBS,2),LQN(2),NBAS(2)
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
      DO I=1,2
        T1  = PI12
        F1  = 0.5D0
        GML = DLOG(T1)
        DO N=2,LQN(I)+2
          GML = GML+DLOG(F1)
          F1  = F1 + 1.0D0
        ENDDO
        RLA = DFLOAT(LQN(I))
        GA1 = TWLG-GML
        RA1 = RLA+1.5D0
        DO IBAS=1,NBAS(I)
          ELOG       = DLOG(2.0D0*EXL(IBAS,I))
          RN(IBAS,I) = DEXP(0.5D0*(GA1+RA1*ELOG))
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE RNSS(RN,EXL,LQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 RRRRRRR  NN    NN  SSSSSS   SSSSSS                   C
C                 RR    RR NNN   NN SS    SS SS    SS                  C
C                 RR    RR NNNN  NN SS       SS                        C
C                 RR    RR NN NN NN  SSSSSS   SSSSSS                   C
C                 RRRRRRR  NN  NNNN       SS       SS                  C
C                 RR    RR NN   NNN SS    SS SS    SS                  C
C                 RR    RR NN    NN  SSSSSS   SSSSSS                   C
C                                                                      C
C -------------------------------------------------------------------- C
C  RNSS GENERATES THE SMALL-SMALL SGTF NORMALISATION CONSTANTS.        C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ LQN  - PAIR OF ORBITAL QUANTUM NUMBERS.                           C
C  ▶ NBAS - PAIR OF PARAMETER LIST LENGTHS.                            C
C  OUTPUT:                                                             C
C  ▶ RN - BATCH OF NORMALISATION CONSTANTS.                            C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RN(MBS,2),EXL(MBS,2),LQN(2),NBAS(2)
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
      DO I=1,2
        T1  = PI12
        F1  = 0.5D0
        GML = DLOG(T1)
        DO N=2,LQN(I)+3
          GML = GML+DLOG(F1)
          F1  = F1+1.0D0
        ENDDO
        RLA = DFLOAT(LQN(I))
        GA1 = TWLG-GML
        RA1 = RLA+0.5D0
        DO IBAS=1,NBAS(I)
          ELOG       = DLOG(2.0D0*EXL(IBAS,I))
          RN(IBAS,I) = DEXP(0.5D0*(GA1+RA1*ELOG))
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE RNLS(RN,EXL,LQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 RRRRRRR  NN    NN LL       SSSSSS                    C
C                 RR    RR NNN   NN LL      SS    SS                   C
C                 RR    RR NNNN  NN LL      SS                         C
C                 RR    RR NN NN NN LL       SSSSSS                    C
C                 RRRRRRR  NN  NNNN LL            SS                   C
C                 RR    RR NN   NNN LL      SS    SS                   C
C                 RR    RR NN    NN LLLLLLLL SSSSSS                    C
C                                                                      C
C -------------------------------------------------------------------- C
C  RNLS GENERATES THE LARGE-SMALL SGTF NORMALISATION CONSTANTS.        C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ LQN  - PAIR OF ORBITAL QUANTUM NUMBERS.                           C
C  ▶ NBAS - PAIR OF PARAMETER LIST LENGTHS.                            C
C  OUTPUT:                                                             C
C  ▶ RN   - BATCH OF NORMALISATION CONSTANTS.                          C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RN(MBS,2),EXL(MBS,2),LQN(2),NBAS(2)
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
      T1L  = PI12
      F1L  = 0.5D0
      GMLL = DLOG(T1L)
C
      DO N=2,LQN(1)+2
        GMLL = GMLL+DLOG(F1L)
        F1L  = F1L+1.0D0
      ENDDO
C
      RLAL = DFLOAT(LQN(1))
      GA1L = TWLG-GMLL
      RA1L = RLAL+1.5D0
C
      DO IBAS=1,NBAS(1)
        ELOGL      = DLOG(2.0D0*EXL(IBAS,1))
        RN(IBAS,1) = DEXP(0.5D0*(GA1L+RA1L*ELOGL))
      ENDDO
C
      T1S  = PI12
      F1S  = 0.5D0
      GMLS = DLOG(T1S)
C
      DO N=2,LQN(2)+3
        GMLS = GMLS+DLOG(F1S)
        F1S  = F1S+1.0D0
      ENDDO
C
      RLAS = DFLOAT(LQN(2))
      GA1S = TWLG-GMLS
      RA1S = RLAS+0.5D0
C
      DO JBAS=1,NBAS(2)
        ELOGS      = DLOG(2.0D0*EXL(JBAS,2))
        RN(JBAS,2) = DEXP(0.5D0*(GA1S+RA1S*ELOGS))
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE RNSL(RN,EXL,LQN,NBAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 RRRRRRR  NN    NN  SSSSSS  LL                        C
C                 RR    RR NNN   NN SS    SS LL                        C
C                 RR    RR NNNN  NN SS       LL                        C
C                 RR    RR NN NN NN  SSSSSS  LL                        C
C                 RRRRRRR  NN  NNNN       SS LL                        C
C                 RR    RR NN   NNN SS    SS LL                        C
C                 RR    RR NN    NN  SSSSSS  LLLLLLLL                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  RNSL GENERATES THE SMALL-LARGE SGTF NORMALISATION CONSTANTS.        C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ LQN  - PAIR OF ORBITAL QUANTUM NUMBERS.                           C
C  ▶ NBAS - PAIR OF PARAMETER LIST LENGTHS.                            C
C  OUTPUT:                                                             C
C  ▶ RN   - BATCH OF NORMALISATION CONSTANTS.                          C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION RN(MBS,2),EXL(MBS,2),LQN(2),NBAS(2)
C
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
C
      T1S  = PI12
      F1S  = 0.5D0
      GMSL = DLOG(T1S)
C
      DO N=2,LQN(1)+3
        GMSL = GMSL+DLOG(F1S)
        F1S  = F1S+1.0D0
      ENDDO
C
      RLAS = DFLOAT(LQN(1))
      GA1S = TWLG-GMSL
      RA1S = RLAS+0.5D0
C
      DO IBAS=1,NBAS(1)
        ELOGS      = DLOG(2.0D0*EXL(IBAS,1))
        RN(IBAS,1) = DEXP(0.5D0*(GA1S+RA1S*ELOGS))
      ENDDO
C
      T1L  = PI12
      F1L  = 0.5D0
      GMLL = DLOG(T1L)
C
      DO N=2,LQN(2)+2
        GMLL = GMLL+DLOG(F1L)
        F1L  = F1L+1.0D0
      ENDDO
C
      RLAL = DFLOAT(LQN(2))
      GA1L = TWLG-GMLL
      RA1L = RLAL+1.5D0
C
      DO JBAS=1,NBAS(2)
        ELOGL      = DLOG(2.0D0*EXL(JBAS,2))
        RN(JBAS,2) = DEXP(0.5D0*(GA1L+RA1L*ELOGL))
      ENDDO
C
      RETURN
      END

