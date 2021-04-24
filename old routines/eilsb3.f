C
C
      SUBROUTINE EILSB3(E11,E21,EXL,XYZ,KQN,MQN,NBS,IPHS,IA1,IA2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           EEEEEEEE IIII LL       SSSSSS  BBBBBBB   333333            C
C           EE        II  LL      SS    SS BB    BB 33    33           C
C           EE        II  LL      SS       BB    BB       33           C
C           EEEEEE    II  LL       SSSSSS  BBBBBBB    33333            C
C           EE        II  LL            SS BB    BB       33           C
C           EE        II  LL      SS    SS BB    BB 33    33           C
C           EEEEEEEE IIII LLLLLLLL SSSSSS  BBBBBBB   333333            C
C                                                                      C
C -------------------------------------------------------------------- C
C  EILSB3 GENERATES A BACTCH OF EQLL COEFFICIENTS FOR BASIS FUNCTION   C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB+1 -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.      C
C  THIS ACTUALLY MAKES A VECTOR LIST OF EQLS NEEDED FOR BREIT.         C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'param.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMPLEX*16 E11(MB2,MEQ,3),E21(MB2,MEQ,3)
      COMPLEX*16 E11X(MB2,MEQ),E11Y(MB2,MEQ),E11Z(MB2,MEQ),
     ^           E21X(MB2,MEQ),E21Y(MB2,MEQ),E21Z(MB2,MEQ)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
C
C     TRANSFER NUCLEAR CENTRE COORDINATES TO LOCAL ARRAY
      DO IX=1,3
        XY2(IX,1) = XYZ(IX,IA1)
        XY2(IX,2) = XYZ(IX,IA2)
      ENDDO
C
C     TRANSFER KQN ENTRIES TO LOCAL ARRAY
      KQ2(1) = KQN(IA1)
      KQ2(2) = KQN(IA2)
C
C     CALCULATE LQN VALUES
      IF(KQ2(1).LT.0) THEN
        LA =-KQ2(1)-1
      ELSE
        LA = KQ2(1)
      ENDIF
      IF(KQ2(2).LT.0) THEN
        LB =-KQ2(2)-1
      ELSE
        LB = KQ2(2)
      ENDIF
C
C     TRANSFER BASIS BLOCKS TO LOCAL ARRAY
      NB2(1) = NBS(IA1)
      DO IBAS=1,NB2(1)
        EX2(IBAS,1) = EXL(IBAS,IA1)
      ENDDO
C
      NB2(2) = NBS(IA2)
      DO JBAS=1,NB2(2)
        EX2(JBAS,2) = EXL(JBAS,IA2)
      ENDDO
C
C     MAXIMUM POLYNOMIAL DEGREE IN FINITE EXPANSION
      LAMLS  = LA+LB+1
      NTUVLS = (LAMLS+1)*(LAMLS+2)*(LAMLS+3)/6
C
C**********************************************************************C
C     1: GENERATE AND STORE E11, FOR MQN PAIRS (-|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) =-MQN(IA1)
      MQ2(2) =-MQN(IA2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLS(E11X,EX2,XY2,KQ2,MQ2,NB2,1)
      CALL EQLS(E11Y,EX2,XY2,KQ2,MQ2,NB2,2)
      CALL EQLS(E11Z,EX2,XY2,KQ2,MQ2,NB2,3)
C
C     MULTIPLY ELSQ BY PHASE TERM FOR R-INTEGRALS (SAVES COMP. TIME)
      DO ITUV=1,NTUVLS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV,1) = PHS*E11X(M,ITUV)
          E11(M,ITUV,2) = PHS*E11Y(M,ITUV)
          E11(M,ITUV,3) = PHS*E11Z(M,ITUV)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     2: GENERATE AND STORE E21, FOR MQN PAIRS (+|MA|,-|MB|)           C
C**********************************************************************C
C
C     BASIS PAIR MQNS
      MQ2(1) = MQN(IA1)
      MQ2(2) =-MQN(IA2)

C     GENERATE THE RAW COEFFICIENTS
      CALL EQLS(E21X,EX2,XY2,KQ2,MQ2,NB2,1)
      CALL EQLS(E21Y,EX2,XY2,KQ2,MQ2,NB2,2)
      CALL EQLS(E21Z,EX2,XY2,KQ2,MQ2,NB2,3)
C
C     MULTIPLY ELSQ BY PHASE TERM FOR R-INTEGRALS (SAVES COMP. TIME)
      DO ITUV=1,NTUVLS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV,1) = PHS*E21X(M,ITUV)
          E21(M,ITUV,2) = PHS*E21Y(M,ITUV)
          E21(M,ITUV,3) = PHS*E21Z(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT E22 AND E12 ARE RELATED TO THESE BY PHASE FACTORS.
C     THERE IS ALSO A RELATION WHICH ALLOWS US TO OBTAIN ESLQ FROM ELSQ.
C
      RETURN
      END

