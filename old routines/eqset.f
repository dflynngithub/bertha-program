C
C
      SUBROUTINE EQLLMK(E11,E21,EXL,XYZ,KQN,MQN,NBS,IPHS,IA1,IA2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      EEEEEEEE  QQQQQQ    LL       LL       MM       MM KK    KK      C
C      EE       QQ    QQ   LL       LL       MMM     MMM KK   KK       C
C      EE      QQ      QQ  LL       LL       MMMM   MMMM KK  KK        C
C      EEEEEE  QQ      QQ  LL       LL       MM MM MM MM KKKKK         C
C      EE      QQ      QQ  LL       LL       MM  MMM  MM KK  KK        C
C      EE       QQ    QQ   LL       LL       MM   M   MM KK   KK       C
C      EEEEEEEE  QQQQQQ QQ LLLLLLLL LLLLLLLL MM       MM KK    KK      C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQLLMK GENERATES A BATCH OF EQLL COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.        C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - BASIS INDICES TO CONSTRUCT EQ-BLOCK FROM.               C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,4),XY2(3,4),KQ2(4),MQ2(4),NB2(4)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
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
      LAMLL  = LA+LB
      NTUVLL = (LAMLL+1)*(LAMLL+2)*(LAMLL+3)/6
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
      CALL EQLL(E11,EX2,XY2,KQ2,MQ2,NB2,IQ)
C                                    ~
C     MULTIPLY EQLL BY PHASE TERM IF EQLL ARE NEEDED
      DO ITUV=1,NTUVLL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
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
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLL(E21,EX2,XY2,KQ2,MQ2,NB2,IQ)
C                                    ~
C     MULTIPLY EQLL BY PHASE TERM IF EQLL ARE NEEDED
      DO ITUV=1,NTUVLL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV) = PHS*E21(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT E22 AND E12 ARE RELATED TO THESE BY PHASE FACTORS.
C
      RETURN
      END
C
C
      SUBROUTINE EQLSMK(E11,E21,EXL,XYZ,KQN,MQN,NBS,IPHS,IA1,IA2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      EEEEEEEE  QQQQQQ    LL       SSSSSS  MM       MM KK    KK       C
C      EE       QQ    QQ   LL      SS    SS MMM     MMM KK   KK        C
C      EE      QQ      QQ  LL      SS       MMMM   MMMM KK  KK         C
C      EEEEEE  QQ      QQ  LL       SSSSSS  MM MM MM MM KKKKK          C
C      EE      QQ      QQ  LL            SS MM  MMM  MM KK  KK         C
C      EE       QQ    QQ   LL      SS    SS MM   M   MM KK   KK        C
C      EEEEEEEE  QQQQQQ QQ LLLLLLLL SSSSSS  MM       MM KK    KK       C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQLSMK GENERATES A BATCH OF EQLS COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB+1 -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.      C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - BASIS INDICES TO CONSTRUCT EQ-BLOCK FROM.               C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
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
      CALL EQLS(E11,EX2,XY2,KQ2,MQ2,NB2,IQ)
C                                    ~
C     MULTIPLY EQLS BY PHASE TERM IF EQLS ARE NEEDED
      DO ITUV=1,NTUVLS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
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
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLS(E21,EX2,XY2,KQ2,MQ2,NB2,IQ)
C                                    ~
C     MULTIPLY EQLS BY PHASE TERM IF EQLS ARE NEEDED
      DO ITUV=1,NTUVLS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV) = PHS*E21(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT E22 AND E12 ARE RELATED TO THESE BY PHASE FACTORS.
C     THERE IS ALSO A RELATION WHICH ALLOWS US TO OBTAIN ESLQ FROM ELSQ.
C
      RETURN
      END
C
C
      SUBROUTINE EQSLMK(E11,E21,EXL,XYZ,KQN,MQN,NBS,IPHS,IA1,IA2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      EEEEEEEE  QQQQQQ    SSSSSS  LL       MM       MM KK    KK       C
C      EE       QQ    QQ  SS    SS LL       MMM     MMM KK   KK        C
C      EE      QQ      QQ SS       LL       MMMM   MMMM KK  KK         C
C      EEEEEE  QQ      QQ  SSSSSS  LL       MM MM MM MM KKKKK          C
C      EE      QQ      QQ       SS LL       MM  MMM  MM KK  KK         C
C      EE       QQ    QQ  SS    SS LL       MM   M   MM KK   KK        C
C      EEEEEEEE  QQQQQQ QQ SSSSSS  LLLLLLLL MM       MM KK    KK       C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQSLMK GENERATES A BATCH OF EQSL COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB+1 -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.      C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - BASIS INDICES TO CONSTRUCT EQ-BLOCK FROM.               C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
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
      LAMSL  = LA+LB+1
      NTUVSL = (LAMSL+1)*(LAMSL+2)*(LAMSL+3)/6
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
      CALL EQSL(E11,EX2,XY2,KQ2,MQ2,NB2,IQ)
C                                    ~
C     MULTIPLY EQSL BY PHASE TERM IF EQSL ARE NEEDED
      DO ITUV=1,NTUVSL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
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
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQSL(E21,EX2,XY2,KQ2,MQ2,NB2,IQ)
C                                    ~
C     MULTIPLY EQSL BY PHASE TERM IF EQSL ARE NEEDED
      DO ITUV=1,NTUVSL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV) = PHS*E21(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT E22 AND E12 ARE RELATED TO THESE BY PHASE FACTORS.
C
      RETURN
      END
C
C
      SUBROUTINE EQSSMK(E11,E21,EXL,XYZ,KQN,MQN,NBS,IPHS,IA1,IA2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      EEEEEEEE  QQQQQQ    SSSSSS   SSSSSS  MM       MM KK    KK       C
C      EE       QQ    QQ  SS    SS SS    SS MMM     MMM KK   KK        C
C      EE      QQ      QQ SS       SS       MMMM   MMMM KK  KK         C
C      EEEEEE  QQ      QQ  SSSSSS   SSSSSS  MM MM MM MM KKKKK          C
C      EE      QQ      QQ       SS       SS MM  MMM  MM KK  KK         C
C      EE       QQ    QQ  SS    SS SS    SS MM   M   MM KK   KK        C
C      EEEEEEEE  QQQQQQ QQ SSSSSS   SSSSSS  MM       MM KK    KK       C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQSSMK GENERATES A BATCH OF EQSS COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB+2 -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.      C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL     - FULL LIST OF BASIS SET PARAMETERS.                      C
C  ▶ XYZ     - FULL LIST OF BASIS FUNCTION CARTESIAN CENTRES.          C
C  ▶ KQN     - FULL LIST OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.      C
C  ▶ MQN     - FULL LIST OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.   C
C  ▶ NBS     - FULL LIST OF BASIS FUNCTION BLOCK LENGTHS.              C
C  ▶ IPHS    - EQ-COEFFICIENT PHASE (NON-TRIVIAL FOR COULOMB/BREIT).   C
C  ▶ IA1/IA2 - BASIS INDICES TO CONSTRUCT EQ-BLOCK FROM.               C
C  ▶ IQ      - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.         C
C  OUTPUT:                                                             C
C  ▶ E11     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (-|MA|,-|MB|)      C
C  ▶ E21     - UNIQUE EQ-COEFFICIENTS FOR MQN PAIRS (+|MA|,-|MB|)      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
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
      LAMSS  = LA+LB+2
      NTUVSS = (LAMSS+1)*(LAMSS+2)*(LAMSS+3)/6
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
      CALL EQSS(E11,EX2,XY2,KQ2,MQ2,NB2,IQ)
C                                    ~
C     MULTIPLY EQSS BY PHASE TERM IF EQSS ARE NEEDED
      DO ITUV=1,NTUVSS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E11(M,ITUV) = PHS*E11(M,ITUV)
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
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQSS(E21,EX2,XY2,KQ2,MQ2,NB2,IQ)
C                                    ~
C     MULTIPLY EQSS BY PHASE TERM IF EQSS ARE NEEDED
      DO ITUV=1,NTUVSS
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV) = PHS*E21(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT E22 AND E12 ARE RELATED TO THESE BY SIMPLE FACTORS.
C
      RETURN
      END
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
C  EILSB3 GENERATES A BATCH OF EILS COEFFICIENTS FOR BASIS FUNCTION    C
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
      INCLUDE 'parameters.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMPLEX*16 E11(MB2,MEQ,3),E21(MB2,MEQ,3)
      COMPLEX*16 E11X(MB2,MEQ),E11Y(MB2,MEQ),E11Z(MB2,MEQ),
     &           E21X(MB2,MEQ),E21Y(MB2,MEQ),E21Z(MB2,MEQ)
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
C                                    ~
C     MULTIPLY EQLS BY PHASE TERM IF EQLS ARE NEEDED
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
C                                    ~
C     MULTIPLY EQLS BY PHASE TERM IF EQLS ARE NEEDED
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
C
C
      SUBROUTINE EISLB3(E11,E21,EXL,XYZ,KQN,MQN,NBS,IPHS,IA1,IA2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          EEEEEEEE IIII  SSSSSS  LL       BBBBBBB   333333            C
C          EE        II  SS    SS LL       BB    BB 33    33           C
C          EE        II  SS       LL       BB    BB       33           C
C          EEEEEE    II   SSSSSS  LL       BBBBBBB    33333            C
C          EE        II        SS LL       BB    BB       33           C
C          EE        II  SS    SS LL       BB    BB 33    33           C
C          EEEEEEEE IIII  SSSSSS  LLLLLLLL BBBBBBB   333333            C
C                                                                      C
C -------------------------------------------------------------------- C
C  EISLB3 GENERATES A BATCH OF EISL COEFFICIENTS FOR BASIS FUNCTION    C
C  OVERLAPS, WITH A FINITE EXPANSION OF LENGTH Λ = (λ+1)(λ+2)(λ+3)/6   C
C  WHERE Λ = LA+LB+1 -- THIS USES THE ALGORITHM OF V.R. SAUNDERS.      C
C  THIS ACTUALLY MAKES A VECTOR LIST OF EQSL NEEDED FOR BREIT (MBPT).  C
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
      INCLUDE 'parameters.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),MQN(4),NBS(4)
      DIMENSION EX2(MBS,2),XY2(3,2),KQ2(2),MQ2(2),NB2(2)
C
      COMPLEX*16 E11(MB2,MEQ,3),E21(MB2,MEQ,3)
      COMPLEX*16 E11X(MB2,MEQ),E11Y(MB2,MEQ),E11Z(MB2,MEQ),
     &           E21X(MB2,MEQ),E21Y(MB2,MEQ),E21Z(MB2,MEQ)
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
      LAMSL  = LA+LB+1
      NTUVSL = (LAMSL+1)*(LAMSL+2)*(LAMSL+3)/6
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
      CALL EQSL(E11X,EX2,XY2,KQ2,MQ2,NB2,1)
      CALL EQSL(E11Y,EX2,XY2,KQ2,MQ2,NB2,2)
      CALL EQSL(E11Z,EX2,XY2,KQ2,MQ2,NB2,3)
C                                    ~
C     MULTIPLY EQSL BY PHASE TERM IF EQSL ARE NEEDED
      DO ITUV=1,NTUVSL
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
      CALL EQSL(E21X,EX2,XY2,KQ2,MQ2,NB2,1)
      CALL EQSL(E21Y,EX2,XY2,KQ2,MQ2,NB2,2)
      CALL EQSL(E21Z,EX2,XY2,KQ2,MQ2,NB2,3)
C                                    ~
C     MULTIPLY EQSL BY PHASE TERM IF EQSL ARE NEEDED
      DO ITUV=1,NTUVSL
        PHS = DFLOAT((IPHS)**(ILAM(ITUV)))
        DO M=1,NB2(1)*NB2(2)
          E21(M,ITUV,1) = PHS*E21X(M,ITUV)
          E21(M,ITUV,2) = PHS*E21Y(M,ITUV)
          E21(M,ITUV,3) = PHS*E21Z(M,ITUV)
        ENDDO
      ENDDO
C
C     NOTE THAT E22 AND E12 ARE RELATED TO THESE BY PHASE FACTORS.
C     THERE IS ALSO A RELATION WHICH ALLOWS US TO OBTAIN ELSQ FROM ESLQ.
C
      RETURN
      END
C
C
      SUBROUTINE EQLL(ELL,EXL,XYZ,KQN,MQN,NBS,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                EEEEEEEE  QQQQQQ    LL       LL                       C
C                EE       QQ    QQ   LL       LL                       C
C                EE      QQ      QQ  LL       LL                       C
C                EEEEEE  QQ      QQ  LL       LL                       C
C                EE      QQ      QQ  LL       LL                       C
C                EE       QQ    QQ   LL       LL                       C
C                EEEEEEEE  QQQQQQ QQ LLLLLLLL LLLLLLLL                 C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQLL GENERATES A BLOCK OF RAW EQLL/GQLL-COEFFICIENTS FOR A GIVEN    C
C  COMBINATION OF MQNS (WITH SIGN), TO BE USED BY EQLLMK.              C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ XYZ  - PAIR OF BASIS FUNCTION CARTESIAN CENTRES.                  C
C  ▶ KQN  - PAIR OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.              C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ NBS  - PAIR OF BASIS FUNCTION BLOCK LENGTHS.                      C
C  ▶ IQ   - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.            C
C  OUTPUT:                                                             C
C  ▶ ELL  - RAW EQLL-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION LQLAB(2),MQLAB(2)
      DIMENSION EXL(MBS,2),XYZ(3,2),KQN(2),LQN(2),JQN(2),MQN(2),NBS(2)
      DIMENSION RN(MBS,2)
C
      COMPLEX*16 CONE
      COMPLEX*16 ELL(MB2,MEQ),ESG(MB2,MEQ)
C
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     CLEBSCH-GORDAN SENSITIVITY PARAMETER
      DATA SENS/1.0D-14/
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     SIGN MULITPLIER FOR σ COUPLING MATRICES
      IF(IQ.EQ.0.OR.IQ.EQ.1) THEN
        SIG = 1.0D0
      ELSEIF(IQ.EQ.2.OR.IQ.EQ.3) THEN
        SIG =-1.0D0
      ENDIF
C
C     INITIALISE COEFFICIENT STORAGE ARRAY TO ZERO
      DO ITUV=1,MEQ
        DO M=1,MB2
          ELL(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     CALCULATE LQN VALUES
      IF(KQN(1).LT.0) THEN
        LQN(1) =-KQN(1)-1
      ELSE
        LQN(1) = KQN(1)
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        LQN(2) =-KQN(2)-1
      ELSE
        LQN(2) = KQN(2)
      ENDIF
C
C     CALCULATE JQN VALUES
      JQN(1) = 2*IABS(KQN(1))-1
      JQN(2) = 2*IABS(KQN(2))-1
C
C     CALCULATE THE APPROPRIATE CLEBSCH-GORDAN FACTORS
      IF(KQN(1).LT.0) THEN
        CAU = DSQRT(DFLOAT(JQN(1)+MQN(1)  )/DFLOAT(2*JQN(1)  ))
        CAL = DSQRT(DFLOAT(JQN(1)-MQN(1)  )/DFLOAT(2*JQN(1)  ))
      ELSE
        CAU =-DSQRT(DFLOAT(JQN(1)-MQN(1)+2)/DFLOAT(2*JQN(1)+4))
        CAL = DSQRT(DFLOAT(JQN(1)+MQN(1)+2)/DFLOAT(2*JQN(1)+4))
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        CBU = DSQRT(DFLOAT(JQN(2)+MQN(2)  )/DFLOAT(2*JQN(2)  ))
        CBL = DSQRT(DFLOAT(JQN(2)-MQN(2)  )/DFLOAT(2*JQN(2)  ))
      ELSE
        CBU =-DSQRT(DFLOAT(JQN(2)-MQN(2)+2)/DFLOAT(2*JQN(2)+4))
        CBL = DSQRT(DFLOAT(JQN(2)+MQN(2)+2)/DFLOAT(2*JQN(2)+4))
      ENDIF
C
C     BASIS FUNCTION OVERLAP LIST LENGTH
      MAXM = NBS(1)*NBS(2)
C
C**********************************************************************C
C     INITIALISE COMMON GEOMETRIC INFORMATION                          C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTRES A AND B
      ABDST(1) = XYZ(1,1)-XYZ(1,2)
      ABDST(2) = XYZ(2,1)-XYZ(2,2)
      ABDST(3) = XYZ(3,1)-XYZ(3,2)
      AB2 = ABDST(1)*ABDST(1) + ABDST(2)*ABDST(2) + ABDST(3)*ABDST(3)
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
C
          EIBS(M) = EXL(IBAS,1)
          EJBS(M) = EXL(JBAS,2)
          P(M)    = EXL(IBAS,1)+EXL(JBAS,2)
          P2(M)   = 2.0D0*P(M)
          P22(M)  = P2(M)*P2(M)
          UV(M)   = EXL(IBAS,1)*EXL(JBAS,2)
          PX      = (EXL(IBAS,1)*XYZ(1,1)+EXL(JBAS,2)*XYZ(1,2))/P(M)
          PY      = (EXL(IBAS,1)*XYZ(2,1)+EXL(JBAS,2)*XYZ(2,2))/P(M)
          PZ      = (EXL(IBAS,1)*XYZ(3,1)+EXL(JBAS,2)*XYZ(3,2))/P(M)
          PAX(M)  = PX-XYZ(1,1)
          PAY(M)  = PY-XYZ(2,1)
          PAZ(M)  = PZ-XYZ(3,1)
          PBX(M)  = PX-XYZ(1,2)
          PBY(M)  = PY-XYZ(2,2)
          PBZ(M)  = PZ-XYZ(3,2)
          PA2(M)  = PAX(M)*PAX(M)+PAY(M)*PAY(M)+PAZ(M)*PAZ(M)
          PB2(M)  = PBX(M)*PBX(M)+PBY(M)*PBY(M)+PBZ(M)*PBZ(M)
          RKAB(M) = DEXP(-(EXL(IBAS,1)*EXL(JBAS,2)*AB2)/P(M))
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CASE 1: ALL KQN(1) AND KQN(2) TYPES.                             C
C**********************************************************************C
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)
      LQLAB(2) = LQN(2)
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ELL0 AND ELLZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLL
          DO M=1,MAXM
            DO ITUV=1,NTUV0
              ELL(M,ITUV) = ELL(M,ITUV) +     CAB*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLL
          DO M=1,MAXM
            DO ITUV=1,NTUV0
              ELL(M,ITUV) = ELL(M,ITUV) + SIG*CAB*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ELLX AND ELLY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLL
          DO M=1,MAXM
            DO ITUV=1,NTUV0
              ELL(M,ITUV) = ELL(M,ITUV) + SIG*CAB*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLL
          DO M=1,MAXM
            DO ITUV=1,NTUV0
              ELL(M,ITUV) = ELL(M,ITUV) +     CAB*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
C     INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
      LAM0  = LQLAB(1)+LQLAB(2)
      NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C     GENERATE RNLL NORMALISATION CONSTANTS
      CALL RNLL(RN,EXL,LQN,NBS)
C
C     NORMALISE THE ELLQ COEFFICIENT BLOCK
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M   = M+1
          RLL = RN(IBAS,1)*RN(JBAS,2)
          DO ITUV=1,NTUV0
            ELL(M,ITUV) = RLL*ELL(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     σ_Y SPECIAL CASE: MULTIPLY ELLY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV0
            ELL(M,ITUV) = CONE*ELL(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE EQLS(ELS,EXL,XYZ,KQN,MQN,NBS,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 EEEEEEEE  QQQQQQ    LL       SSSSSS                  C
C                 EE       QQ    QQ   LL      SS    SS                 C
C                 EE      QQ      QQ  LL      SS                       C
C                 EEEEEE  QQ      QQ  LL       SSSSSS                  C
C                 EE      QQ      QQ  LL            SS                 C
C                 EE       QQ    QQ   LL      SS    SS                 C
C                 EEEEEEEE  QQQQQQ QQ LLLLLLLL SSSSSS                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQLS GENERATES A BLOCK OF RAW EQLS-COEFFICIENTS FOR A GIVEN         C
C  COMBINATION OF MQNS (WITH SIGN), TO BE USED BY EQLSMK.              C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ XYZ  - PAIR OF BASIS FUNCTION CARTESIAN CENTRES.                  C
C  ▶ KQN  - PAIR OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.              C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ NBS  - PAIR OF BASIS FUNCTION BLOCK LENGTHS.                      C
C  ▶ IQ   - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.            C
C  OUTPUT:                                                             C
C  ▶ ELS  - RAW EQLS-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION LQLAB(2),MQLAB(2)
      DIMENSION EXL(MBS,2),XYZ(3,2),KQN(2),LQN(2),JQN(2),MQN(2),NBS(2)
      DIMENSION RN(MBS,2)
      DIMENSION T2(MB2),T0(MB2)
C
      COMPLEX*16 CONE
      COMPLEX*16 ELS(MB2,MEQ),ESG(MB2,MEQ),ENSG(MB2,MEQ)
C
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     CLEBSCH-GORDAN SENSITIVITY PARAMETER
      DATA SENS/1.0D-14/
C
C     SIGN MULITPLIER FOR σ COUPLING MATRICES
      IF(IQ.EQ.0.OR.IQ.EQ.1) THEN
        SIG = 1.0D0
      ELSEIF(IQ.EQ.2.OR.IQ.EQ.3) THEN
        SIG =-1.0D0
      ENDIF
C
C     INITIALISE COEFFICIENT STORAGE ARRAY TO ZERO
      DO M=1,MB2
        DO ITUV=1,MEQ
          ELS(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     CALCULATE LQN VALUES
      IF(KQN(1).LT.0) THEN
        LQN(1) =-KQN(1)-1
      ELSE
        LQN(1) = KQN(1)
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        LQN(2) =-KQN(2)-1
      ELSE
        LQN(2) = KQN(2)
      ENDIF
C
C     CALCULATE JQN VALUES
      JQN(1) = 2*IABS(KQN(1))-1
      JQN(2) = 2*IABS(KQN(2))-1
C
C     CALCULATE THE APPROPRIATE CLEBSCH-GORDAN FACTORS
      IF(KQN(1).LT.0) THEN
        CAU = DSQRT(DFLOAT(JQN(1)+MQN(1)  )/DFLOAT(2*JQN(1)  ))
        CAL = DSQRT(DFLOAT(JQN(1)-MQN(1)  )/DFLOAT(2*JQN(1)  ))
      ELSE
        CAU =-DSQRT(DFLOAT(JQN(1)-MQN(1)+2)/DFLOAT(2*JQN(1)+4))
        CAL = DSQRT(DFLOAT(JQN(1)+MQN(1)+2)/DFLOAT(2*JQN(1)+4))
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        CBU =-DSQRT(DFLOAT(JQN(2)-MQN(2)+2)/DFLOAT(2*JQN(2)+4))
        CBL = DSQRT(DFLOAT(JQN(2)+MQN(2)+2)/DFLOAT(2*JQN(2)+4))
      ELSE
        CBU = DSQRT(DFLOAT(JQN(2)+MQN(2)  )/DFLOAT(2*JQN(2)  ))
        CBL = DSQRT(DFLOAT(JQN(2)-MQN(2)  )/DFLOAT(2*JQN(2)  ))
      ENDIF
C
C     DETERMINE THE NUMBER OF FUNCTIONS ON EACH CENTRE
      MAXM = NBS(1)*NBS(2)
C
C     KINETIC PRE-FACTORS FOR THIS BLOCK
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
          T2(M) =-2.0D0*EXL(JBAS,2)
          T0(M) = DFLOAT(2*LQN(2)+1)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INITIALISE COMMON GEOMETRIC INFORMATION                          C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTRES A AND B
      ABDST(1) = XYZ(1,1)-XYZ(1,2)
      ABDST(2) = XYZ(2,1)-XYZ(2,2)
      ABDST(3) = XYZ(3,1)-XYZ(3,2)
      AB2 = ABDST(1)*ABDST(1) + ABDST(2)*ABDST(2) + ABDST(3)*ABDST(3)
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
C
          EIBS(M) = EXL(IBAS,1)
          EJBS(M) = EXL(JBAS,2)
          P(M)    = EXL(IBAS,1)+EXL(JBAS,2)
          P2(M)   = 2.0D0*P(M)
          P22(M)  = P2(M)*P2(M)
          UV(M)   = EXL(IBAS,1)*EXL(JBAS,2)
          PX      = (EXL(IBAS,1)*XYZ(1,1)+EXL(JBAS,2)*XYZ(1,2))/P(M)
          PY      = (EXL(IBAS,1)*XYZ(2,1)+EXL(JBAS,2)*XYZ(2,2))/P(M)
          PZ      = (EXL(IBAS,1)*XYZ(3,1)+EXL(JBAS,2)*XYZ(3,2))/P(M)
          PAX(M)  = PX-XYZ(1,1)
          PAY(M)  = PY-XYZ(2,1)
          PAZ(M)  = PZ-XYZ(3,1)
          PBX(M)  = PX-XYZ(1,2)
          PBY(M)  = PY-XYZ(2,2)
          PBZ(M)  = PZ-XYZ(3,2)
          PA2(M)  = PAX(M)*PAX(M) + PAY(M)*PAY(M) + PAZ(M)*PAZ(M)
          PB2(M)  = PBX(M)*PBX(M) + PBY(M)*PBY(M) + PBZ(M)*PBZ(M)
          RKAB(M) = DEXP(-EXL(IBAS,1)*EXL(JBAS,2)*AB2/P(M))
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CASE 1: KQN(2).LT.0                                              C
C**********************************************************************C
C
      IF(KQN(2).GT.0) GOTO 100
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)
      LQLAB(2) = LQN(2)+1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ELS0 AND ELSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ELSX AND ELSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
      ENDIF
C
100   CONTINUE
C
C**********************************************************************C
C     CASE 2: KQN(2).GT.0                                              C
C**********************************************************************C
C
      IF(KQN(2).LT.0) GOTO 200
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)
      LQLAB(2) = LQN(2)-1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ELS0 AND ELSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
              ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ELSX AND ELSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
              ELS(M,ITUV) = ELS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQLS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
              ELS(M,ITUV) = ELS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF

200   CONTINUE
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
C     MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
      LAM2  = LQN(1)+LQN(2)+2
      NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C     GENERATE RNLS NORMALISATION CONSTANTS
      CALL RNLS(RN,EXL,LQN,NBS)
C
C     NORMALISE THE ELSQ COEFFICIENT BLOCK AND MULTIPLY BY FACTOR i
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M   = M+1
          RLS = RN(IBAS,1)*RN(JBAS,2)
          DO ITUV=1,NTUV2
            ELS(M,ITUV) = CONE*RLS*ELS(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     σ_Y SPECIAL CASE: MULTIPLY ELSY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV2
            ELS(M,ITUV) = CONE*ELS(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE EQSL(ESL,EXL,XYZ,KQN,MQN,NBS,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 EEEEEEEE  QQQQQQ    SSSSSS  LL                       C
C                 EE       QQ    QQ  SS    SS LL                       C
C                 EE      QQ      QQ SS       LL                       C
C                 EEEEEE  QQ      QQ  SSSSSS  LL                       C
C                 EE      QQ      QQ       SS LL                       C
C                 EE       QQ    QQ  SS    SS LL                       C
C                 EEEEEEEE  QQQQQQ QQ SSSSSS  LLLLLLLL                 C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQSL GENERATES A BLOCK OF RAW EQSL-COEFFICIENTS FOR A GIVEN         C
C  COMBINATION OF MQNS (WITH SIGN), TO BE USED BY EQSLMK.              C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ XYZ  - PAIR OF BASIS FUNCTION CARTESIAN CENTRES.                  C
C  ▶ KQN  - PAIR OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.              C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ NBS  - PAIR OF BASIS FUNCTION BLOCK LENGTHS.                      C
C  ▶ IQ   - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.            C
C  OUTPUT:                                                             C
C  ▶ ESL  - RAW EQSL-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION LQLAB(2),MQLAB(2)
      DIMENSION EXL(MBS,2),XYZ(3,2),KQN(2),LQN(2),JQN(2),MQN(2),NBS(2)
      DIMENSION RN(MBS,2)
      DIMENSION T2(MB2),T0(MB2)
C
      COMPLEX*16 CONE
      COMPLEX*16 ESL(MB2,MEQ),ESG(MB2,MEQ),ENSG(MB2,MEQ)
C
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     CLEBSCH-GORDAN SENSITIVITY PARAMETER
      DATA SENS/1.0D-14/
C
C     SIGN MULITPLIER FOR σ COUPLING MATRICES
      IF(IQ.EQ.0.OR.IQ.EQ.1) THEN
        SIG = 1.0D0
      ELSEIF(IQ.EQ.2.OR.IQ.EQ.3) THEN
        SIG =-1.0D0
      ENDIF
C
C     INITIALISE COEFFICIENT STORAGE ARRAY TO ZERO
      DO M=1,MB2
        DO ITUV=1,MEQ
          ESL(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     CALCULATE LQN VALUES
      IF(KQN(1).LT.0) THEN
        LQN(1) =-KQN(1)-1
      ELSE
        LQN(1) = KQN(1)
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        LQN(2) =-KQN(2)-1
      ELSE
        LQN(2) = KQN(2)
      ENDIF
C
C     CALCULATE JQN VALUES
      JQN(1) = 2*IABS(KQN(1))-1
      JQN(2) = 2*IABS(KQN(2))-1
C
C     CALCULATE THE APPROPRIATE CLEBSCH-GORDAN FACTORS
      IF(KQN(1).LT.0) THEN
        CAU =-DSQRT(DFLOAT(JQN(1)-MQN(1)+2)/DFLOAT(2*JQN(1)+4))
        CAL = DSQRT(DFLOAT(JQN(1)+MQN(1)+2)/DFLOAT(2*JQN(1)+4))
      ELSE
        CAU = DSQRT(DFLOAT(JQN(1)+MQN(1)  )/DFLOAT(2*JQN(1)  ))
        CAL = DSQRT(DFLOAT(JQN(1)-MQN(1)  )/DFLOAT(2*JQN(1)  ))
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        CBU = DSQRT(DFLOAT(JQN(2)+MQN(2)  )/DFLOAT(2*JQN(2)  ))
        CBL = DSQRT(DFLOAT(JQN(2)-MQN(2)  )/DFLOAT(2*JQN(2)  ))
      ELSE
        CBU =-DSQRT(DFLOAT(JQN(2)-MQN(2)+2)/DFLOAT(2*JQN(2)+4))
        CBL = DSQRT(DFLOAT(JQN(2)+MQN(2)+2)/DFLOAT(2*JQN(2)+4))
      ENDIF
C
C     DETERMINE THE NUMBER OF FUNCTIONS ON EACH CENTRE
      MAXM = NBS(1)*NBS(2)
C
C     KINETIC PRE-FACTORS FOR THIS BLOCK
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
          T2(M) =-2.0D0*EXL(IBAS,1)
          T0(M) = DFLOAT(2*LQN(1)+1)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INITIALISE COMMON GEOMETRIC INFORMATION                          C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTRES A AND B
      ABDST(1) = XYZ(1,1)-XYZ(1,2)
      ABDST(2) = XYZ(2,1)-XYZ(2,2)
      ABDST(3) = XYZ(3,1)-XYZ(3,2)
      AB2 = ABDST(1)*ABDST(1) + ABDST(2)*ABDST(2) + ABDST(3)*ABDST(3)
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
C
          EIBS(M) = EXL(IBAS,1)
          EJBS(M) = EXL(JBAS,2)
          P(M)    = EXL(IBAS,1)+EXL(JBAS,2)
          P2(M)   = 2.0D0*P(M)
          P22(M)  = P2(M)*P2(M)
          UV(M)   = EXL(IBAS,1)*EXL(JBAS,2)
          PX      = (EXL(IBAS,1)*XYZ(1,1)+EXL(JBAS,2)*XYZ(1,2))/P(M)
          PY      = (EXL(IBAS,1)*XYZ(2,1)+EXL(JBAS,2)*XYZ(2,2))/P(M)
          PZ      = (EXL(IBAS,1)*XYZ(3,1)+EXL(JBAS,2)*XYZ(3,2))/P(M)
          PAX(M)  = PX-XYZ(1,1)
          PAY(M)  = PY-XYZ(2,1)
          PAZ(M)  = PZ-XYZ(3,1)
          PBX(M)  = PX-XYZ(1,2)
          PBY(M)  = PY-XYZ(2,2)
          PBZ(M)  = PZ-XYZ(3,2)
          PA2(M)  = PAX(M)*PAX(M) + PAY(M)*PAY(M) + PAZ(M)*PAZ(M)
          PB2(M)  = PBX(M)*PBX(M) + PBY(M)*PBY(M) + PBZ(M)*PBZ(M)
          RKAB(M) = DEXP(-(EXL(IBAS,1)*EXL(JBAS,2)*AB2)/P(M))
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CASE 1: KQN(2).LT.0                                              C
C**********************************************************************C
C
      IF(KQN(1).GT.0) GOTO 100
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)+1
      LQLAB(2) = LQN(2)
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESL0 AND ESLZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ESLX AND ESLY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
100   CONTINUE
C
C**********************************************************************C
C     CASE 2: KQN(1).GT.0                                              C
C**********************************************************************C
C
      IF(KQN(1).LT.0) GOTO 200
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)-1
      LQLAB(2) = LQN(2)
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESL0 AND ESLZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
              ESL(M,ITUV) = ESL(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
              ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ESLX AND ESLY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
              ESL(M,ITUV) = ESL(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL
          DO M=1,MAXM
            TK = CAB*T0(M)
            DO ITUV=1,NTUV0
              ESL(M,ITUV) = ESL(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSL [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T2(M)
            DO ITUV=1,NTUV2
              ESL(M,ITUV) = ESL(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF

200   CONTINUE
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
C     MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
      LAM2  = LQN(1)+LQN(2)+2
      NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C     GENERATE RNSL NORMALISATION CONSTANTS
      CALL RNSL(RN,EXL,LQN,NBS)
C
C     NORMALISE THE ESLQ COEFFICIENT BLOCK AND MULTIPLY BY FACTOR -i
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M   = M+1
          RSL = RN(IBAS,1)*RN(JBAS,2)
          DO ITUV=1,NTUV2
            ESL(M,ITUV) =-CONE*RSL*ESL(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     σ_Y SPECIAL CASE: MULTIPLY ESLY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV2
            ESL(M,ITUV) = CONE*ESL(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE EQSS(ESS,EXL,XYZ,KQN,MQN,NBS,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 EEEEEEEE  QQQQQQ    SSSSSS   SSSSSS                  C
C                 EE       QQ    QQ  SS    SS SS    SS                 C
C                 EE      QQ      QQ SS       SS                       C
C                 EEEEEE  QQ      QQ  SSSSSS   SSSSSS                  C
C                 EE      QQ      QQ       SS       SS                 C
C                 EE       QQ    QQ  SS    SS SS    SS                 C
C                 EEEEEEEE  QQQQQQ QQ SSSSSS   SSSSSS                  C
C                                                                      C
C -------------------------------------------------------------------- C
C  EQSS GENERATES A BLOCK OF RAW EQSS-COEFFICIENTS FOR A GIVEN         C
C  COMBINATION OF MQNS (WITH SIGN), TO BE USED BY EQSSMK.              C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ EXL  - PAIR OF BASIS SET PARAMETER LISTS.                         C
C  ▶ XYZ  - PAIR OF BASIS FUNCTION CARTESIAN CENTRES.                  C
C  ▶ KQN  - PAIR OF BASIS FUNCTION KAPPA QUANTUM NUMBERS.              C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ NBS  - PAIR OF BASIS FUNCTION BLOCK LENGTHS.                      C
C  ▶ IQ   - PAULI SIGMA COUPLING MATRIX FOR SPINOR OVERLAP.            C
C  OUTPUT:                                                             C
C  ▶ ESS  - RAW EQSS-COEFFICIENTS FOR THIS MQN PAIR (MA,MB).           C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION LQLAB(2),MQLAB(2)
      DIMENSION EXL(MBS,2),XYZ(3,2),KQN(2),LQN(2),JQN(2),MQN(2),NBS(2)
      DIMENSION RN(MBS,2)
      DIMENSION T22(MB2),T20(MB2),T02(MB2),T00(MB2)
C
      COMPLEX*16 CONE
      COMPLEX*16 ESS(MB2,MEQ),ESG(MB2,MEQ),ENSG(MB2,MEQ)
C
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     CLEBSCH-GORDAN SENSITIVITY PARAMETER
      DATA SENS/1.0D-14/
C
C     SIGN MULITPLIER FOR σ COUPLING MATRICES
      IF(IQ.EQ.0.OR.IQ.EQ.1) THEN
        SIG = 1.0D0
      ELSEIF(IQ.EQ.2.OR.IQ.EQ.3) THEN
        SIG =-1.0D0
      ENDIF
C
C     INITIALISE COEFFICIENT STORAGE ARRAY TO ZERO
      DO M=1,MB2
        DO ITUV=1,MEQ
          ESS(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     CALCULATE LQN VALUES
      IF(KQN(1).LT.0) THEN
        LQN(1) =-KQN(1)-1
      ELSE
        LQN(1) = KQN(1)
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        LQN(2) =-KQN(2)-1
      ELSE
        LQN(2) = KQN(2)
      ENDIF
C
C     CALCULATE JQN VALUES
      JQN(1) = 2*IABS(KQN(1))-1
      JQN(2) = 2*IABS(KQN(2))-1
C
C     CALCULATE THE APPROPRIATE CLEBSCH-GORDAN FACTORS
      IF(KQN(1).LT.0) THEN
        CAU =-DSQRT(DFLOAT(JQN(1)-MQN(1)+2)/DFLOAT(2*JQN(1)+4))
        CAL = DSQRT(DFLOAT(JQN(1)+MQN(1)+2)/DFLOAT(2*JQN(1)+4))
      ELSE
        CAU = DSQRT(DFLOAT(JQN(1)+MQN(1)  )/DFLOAT(2*JQN(1)  ))
        CAL = DSQRT(DFLOAT(JQN(1)-MQN(1)  )/DFLOAT(2*JQN(1)  ))
      ENDIF
C
      IF(KQN(2).LT.0) THEN
        CBU =-DSQRT(DFLOAT(JQN(2)-MQN(2)+2)/DFLOAT(2*JQN(2)+4))
        CBL = DSQRT(DFLOAT(JQN(2)+MQN(2)+2)/DFLOAT(2*JQN(2)+4))
      ELSE
        CBU = DSQRT(DFLOAT(JQN(2)+MQN(2)  )/DFLOAT(2*JQN(2)  ))
        CBL = DSQRT(DFLOAT(JQN(2)-MQN(2)  )/DFLOAT(2*JQN(2)  ))
      ENDIF
C
C     DETERMINE THE NUMBER OF FUNCTIONS ON EACH CENTRE
      MAXM = NBS(1)*NBS(2)
C
C     KINETIC PRE-FACTORS FOR THIS BLOCK
      RL1 = DFLOAT(2*LQN(1)+1)
      RL2 = DFLOAT(2*LQN(2)+1)
C
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
          T22(M) = 4.0D0*EXL(IBAS,1)*EXL(JBAS,2)
          T20(M) =-2.0D0*RL2*EXL(IBAS,1)
          T02(M) =-2.0D0*RL1*EXL(JBAS,2)
          T00(M) = RL1*RL2
        ENDDO
      ENDDO
C
C**********************************************************************C
C     INITIALISE COMMON GEOMETRIC INFORMATION                          C
C**********************************************************************C
C
C     EUCLIDIAN DISTANCE BETWEEN CENTRES A AND B
      ABDST(1) = XYZ(1,1)-XYZ(1,2)
      ABDST(2) = XYZ(2,1)-XYZ(2,2)
      ABDST(3) = XYZ(3,1)-XYZ(3,2)
      AB2 = ABDST(1)*ABDST(1) + ABDST(2)*ABDST(2) + ABDST(3)*ABDST(3)
C
C     GAUSSIAN PRODUCT THEOREM IMPLEMENTATION
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M = M+1
C
          EIBS(M) = EXL(IBAS,1)
          EJBS(M) = EXL(JBAS,2)
          P(M)    = EXL(IBAS,1)+EXL(JBAS,2)
          P2(M)   = 2.0D0*P(M)
          P22(M)  = P2(M)*P2(M)
          UV(M)   = EXL(IBAS,1)*EXL(JBAS,2)
          PX      = (EXL(IBAS,1)*XYZ(1,1)+EXL(JBAS,2)*XYZ(1,2))/P(M)
          PY      = (EXL(IBAS,1)*XYZ(2,1)+EXL(JBAS,2)*XYZ(2,2))/P(M)
          PZ      = (EXL(IBAS,1)*XYZ(3,1)+EXL(JBAS,2)*XYZ(3,2))/P(M)
          PAX(M)  = PX-XYZ(1,1)
          PAY(M)  = PY-XYZ(2,1)
          PAZ(M)  = PZ-XYZ(3,1)
          PBX(M)  = PX-XYZ(1,2)
          PBY(M)  = PY-XYZ(2,2)
          PBZ(M)  = PZ-XYZ(3,2)
          PA2(M)  = PAX(M)*PAX(M)+PAY(M)*PAY(M)+PAZ(M)*PAZ(M)
          PB2(M)  = PBX(M)*PBX(M)+PBY(M)*PBY(M)+PBZ(M)*PBZ(M)
          RKAB(M) = DEXP(-(EXL(IBAS,1)*EXL(JBAS,2)*AB2)/P(M))
C
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CASE 1: KQN(1).LT.0 AND KQN(2).LT.0                              C
C**********************************************************************C
C
      IF(KQN(1).GT.0.OR.KQN(2).GT.0) GOTO 100
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)+1
      LQLAB(2) = LQN(2)+1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESS0 AND ESSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ESSX AND ESSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
100   CONTINUE
C
C**********************************************************************C
C     CASE 2: KQN(1).LT.0 AND KQN(2).GT.0                              C
C**********************************************************************C
C
      IF(KQN(1).GT.0.OR.KQN(2).LT.0) GOTO 200
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)+1
      LQLAB(2) = LQN(2)-1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESS0 AND ESSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) +    TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ESSX AND ESSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
200   CONTINUE
C
C**********************************************************************C
C     CASE 3: KQN(1).GT.0 AND KQN(2).LT.0                              C
C**********************************************************************C
C
      IF(KQN(1).LT.0.OR.KQN(2).GT.0) GOTO 300
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)-1
      LQLAB(2) = LQN(2)+1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESS0 AND ESSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ESSX AND ESSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF

300   CONTINUE      
C
C**********************************************************************C
C     CASE 4: KQN(1).GT.0 AND KQN(2).GT.0                              C
C**********************************************************************C
C
      IF(KQN(1).LT.0.OR.KQN(2).LT.0) GOTO 400
C
C     GENERATING LQN LABELS
      LQLAB(1) = LQN(1)-1
      LQLAB(2) = LQN(2)-1
C
C     TERMS 11 AND 22 ARE ONLY NECESSARY FOR ESS0 AND ESSZ
      IF(IQ.EQ.0.OR.IQ.EQ.3) THEN
C
C >>    TERM 11: (MA-1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-UPPER
        CAB = CAU*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          LAM4  = LQLAB(1)+LQLAB(2)+4
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
          NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T00(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ESG)
          CALL ESTEPN(ENSG,ESG,LAM2,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV4
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 22: (MA+1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-LOWER
        CAB = CAL*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          LAM4  = LQLAB(1)+LQLAB(2)+4
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
          NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T00(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ESG)
          CALL ESTEPN(ENSG,ESG,LAM2,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV4
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
C     TERMS 12 AND 21 ARE ONLY NECESSARY FOR ESSX AND ESSY
      IF(IQ.EQ.1.OR.IQ.EQ.2) THEN
C
C >>    TERM 12: (MA-1/2, MB+1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)-1)/2
        MQLAB(2) = (MQN(2)+1)/2
C
C       CLEBSCH-GORDAN PRODUCT: UPPER-LOWER
        CAB = CAU*CBL
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          LAM4  = LQLAB(1)+LQLAB(2)+4
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
          NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T00(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ESG)
          CALL ESTEPN(ENSG,ESG,LAM2,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV4
              ESS(M,ITUV) = ESS(M,ITUV) + SIG*TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
C
C >>    TERM 21: (MA+1/2, MB-1/2) CONTRIBUTIONS
C
C       BASIS PAIR MQN LABELS
        MQLAB(1) = (MQN(1)+1)/2
        MQLAB(2) = (MQN(2)-1)/2
C
C       CLEBSCH-GORDAN PRODUCT: LOWER-UPPER
        CAB = CAL*CBU
C
C       SCREEN CONTRIBUTION BY ANGULAR PRODUCT PAIR
        IF(DABS(CAB).GE.SENS) THEN
C
C         INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
          LAM0  = LQLAB(1)+LQLAB(2)
          LAM2  = LQLAB(1)+LQLAB(2)+2
          LAM4  = LQLAB(1)+LQLAB(2)+4
          NTUV0 = (LAM0+1)*(LAM0+2)*(LAM0+3)/6
          NTUV2 = (LAM2+1)*(LAM2+2)*(LAM2+3)/6
          NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C         GENERATE RAW SPHERICAL GAUSSIAN (ES) COEFFICIENTS
          CALL ESGTF(ESG,LQLAB,MQLAB,MAXM)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS
          DO M=1,MAXM
            TK = CAB*T00(M)
            DO ITUV=1,NTUV0
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=0]
          DO M=1,MAXM
            TK = CAB*T20(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE B (WRITE TO ENSG)
          CALL ESTEPN(ESG,ENSG,LAM0,MAXM,2)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=0,NB=1]
          DO M=1,MAXM
            TK = CAB*T02(M)
            DO ITUV=1,NTUV2
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ENSG(M,ITUV)
            ENDDO
          ENDDO
C
C         INCREASE THE INDEX N BY ONE ON CENTRE A (WRITE TO ESG)
          CALL ESTEPN(ENSG,ESG,LAM2,MAXM,1)
C
C         ADD THIS ANGULAR-(ES) CONTRIBUTION TO EQSS [NA=1,NB=1]
          DO M=1,MAXM
            TK = CAB*T22(M)
            DO ITUV=1,NTUV4
              ESS(M,ITUV) = ESS(M,ITUV) +     TK*ESG(M,ITUV)
            ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C
400   CONTINUE
C
C**********************************************************************C
C     GAUSSIAN NORMALISATION FACTORS                                   C
C**********************************************************************C
C
C     GENERATE RNSS NORMALISATION CONSTANTS
      CALL RNSS(RN,EXL,LQN,NBS)
C
C     MAX POLYNOMIAL DEGREE OF HGTFS IN FINITE EXPANSION
      LAM4  = LQN(1)+LQN(2)+4
      NTUV4 = (LAM4+1)*(LAM4+2)*(LAM4+3)/6
C
C     NORMALISE THE ESSQ COEFFICIENT BLOCK
      M = 0
      DO IBAS=1,NBS(1)
        DO JBAS=1,NBS(2)
          M   = M+1
          RSS = RN(IBAS,1)*RN(JBAS,2)
          DO ITUV=1,NTUV4
            ESS(M,ITUV) = RSS*ESS(M,ITUV)
          ENDDO
        ENDDO
      ENDDO
C
C     σ_Y SPECIAL CASE: MULTIPLY ESSY RESULTS BY i.
      IF(IQ.EQ.2) THEN
        DO M=1,MAXM
          DO ITUV=1,NTUV4
            ESS(M,ITUV) = CONE*ESS(M,ITUV)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE ESGTF(ESG,LQN,MQN,MAXM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C              EEEEEEEE SSSSSS   GGGGGG TTTTTTTT FFFFFFFF              C
C              EE      SS    SS GG    GG   TT    FF                    C
C              EE      SS       GG         TT    FF                    C
C              EEEEEE   SSSSSS  GG         TT    FFFFFF                C
C              EE            SS GG   GGG   TT    FF                    C
C              EE      SS    SS GG    GG   TT    FF                    C
C              EEEEEEEE SSSSSS   GGGGGG    TT    FF                    C
C                                                                      C
C -------------------------------------------------------------------- C
C  ESGTF CONSTRUCTS THE EXPANSION COEFFICIENTS OF THE OVERLAP DENSITY  C
C  OF TWO SPHERICAL HARMONIC FUNCTIONS IN AN AUXILIARY HGTF BASIS.     C
C                                                                      C
C  THE OVERLAP DENSITY IS DEFINED BY Y*[L,M]Y[L',M'], WHERE Y[L,M] ARE C
C  SPHERICAL HARMONICS FOLLOWING THE CONDON-SHORTLEY PHASE CONVENTION. C
C                                                                      C
C  THE REQUIRED COEFFICIENTS ARE GENERATED BY A CALL TO EVRS, WHICH IS C
C  CONSTRUCTED ACCORDING TO THE RECURRENCE RELATIONS DEFINED BY        C
C  V.R.SAUNDERS. THE OUTPUT OF EVRS IS THEN ADJUSTED TO INCLUDE THE    C
C  ANGULAR NORMALISATION CONSTANTS, AS WELL AS A PHASE FACTOR TO       C
C  CONVERT FROM THE SCHIFF TO THE CONDON-SHORTLEY PHASE CONVENTION.    C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ LQN  - PAIR OF BASIS SET ORBITAL QUANTUM NUMBERS.                 C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ MAXM - SIZE OF THIS BLOCK.                                        C
C  OUTPUT:                                                             C
C  ▶ ESG  - EXPANSION COEFFICIENTS FOR EACH OVERLAP IN THE BLOCK.      C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION LQN(2),MQN(2),MQNLAB(2)
C
      COMPLEX*16 ESG(MB2,MEQ)
C
      COMMON/FCTS/RFACT(0:80),SFACT(0:80)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     EVRS IS CALLED WITH THE SIGN OF MQN(1) REVERSED
C     TO AFFECT COMPLEX CONJUGATION, ALONG WITH THE REQUISITE
C     PHASE, WHICH IS CALCULATED LATER.
      MQNLAB(1) =-MQN(1)
      MQNLAB(2) = MQN(2)
C
C     GENERATE RAW COEFFICIENTS WITH EVRS
      CALL EVRS(ESG,LQN,MQNLAB,MAXM)
C
C     MQN PHASES
      PHS1 = (-1.0D0)**((MQN(1)+IABS(MQN(1)))/2)
      PHS2 = (-1.0D0)**((MQN(2)+IABS(MQN(2)))/2)
C
C     CG COEFFICIENTS
      PI4 = 0.25D0/PI
      DGL = DFLOAT((2*LQN(1)+1)*(2*LQN(2)+1))
      CG1 = RFACT(LQN(1)-IABS(MQN(1)))/RFACT(LQN(1)+IABS(MQN(1)))
      CG2 = RFACT(LQN(2)-IABS(MQN(2)))/RFACT(LQN(2)+IABS(MQN(2)))
C
C     ANGULAR NORMALISATION CONSTANT
      ANG = PHS1*PHS2*PI4*DSQRT(DGL*CG1*CG2)
C
C     INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
      LAM  = LQN(1)+LQN(2)
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     APPLY ANGULAR FACTOR TO RAW COEFFICIENTS
      DO M=1,MAXM
        DO ITUV=1,NTUV
          ESG(M,ITUV) = ANG*ESG(M,ITUV)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE EVRS(ESG,LQN,MQN,MAXM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C                 EEEEEEEE VV    VV RRRRRRR   SSSSSS                   C
C                 EE       VV    VV RR    RR SS    SS                  C
C                 EE       VV    VV RR    RR SS                        C
C                 EEEEEE   VV    VV RR    RR  SSSSSS                   C
C                 EE        VV  VV  RRRRRRR        SS                  C
C                 EE         VVVV   RR    RR SS    SS                  C
C                 EEEEEEEE    VV    RR    RR  SSSSSS                   C
C                                                                      C
C -------------------------------------------------------------------- C
C  EVRS EVALUATES THE EXPANSION COEFFICIENTS OF THE OVERLAP CHARGE     C
C  DENSITY OF SGTFS IN AN AUXILIARY HGTF. COEFFICIENTS ARE EVALUATED   C
C  USING THE RECURRENCE RELATIONS DEFINED BY VIC SAUNDERS IN:          C
C                                                                      C
C  V.R.SAUNDERS, "MOLECULAR INTEGRALS FOR GAUSSIAN-TYPE FUNCTIONS",    C
C  METHODS OF COMPUTATIONAL MOLECULAR PHYSICS, DIERCKSEN AND WILSON,   C
C  pp 1-26, REIDEL PUBLISHING, DORDRECHT (1983).                       C
C                                                                      C
C  THE EQ-COEFFS IN THIS PROCEDURE ARE FOR AN UN-NORMALISED SGTF.      C
C  THE COEFFICIENTS ARE DETERMINED ACCORDING TO THE SAME RULES AS      C
C  DEFINED IN THE ABOVE ARTICLE. CONSEQUENTLY, IT SHOULD BE NOTED THAT C
C  THE COEFFICIENTS ARE THOSE OF SPHERICAL HARMONIC FUNCTIONS THAT ARE C
C   ▶ UN-NORMALISED                                                    C
C   ▶ SATISFY THE SCHIFF PHASE CONVENTION.                             C
C                                                                      C
C  THE OUTLINE FOR THE GENERATION OF EQ-COEFFICIENTS IS TAKEN FROM p16 C
C  OF THE ABOVE ARTICLE. EQUATION NUMBERS ARE GIVEN IN COMMENTS.       C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ LQN  - PAIR OF BASIS SET ORBITAL QUANTUM NUMBERS.                 C
C  ▶ MQN  - PAIR OF BASIS FUNCTION MAGNETIC QUANTUM NUMBERS.           C
C  ▶ MAXM - SIZE OF THIS BLOCK.                                        C
C  OUTPUT:                                                             C
C  ▶ ESG  - UN-NORMALISED EXPANSION COEFFICIENTS FOR THIS BLOCK.       C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION NBAS(2),LQN(2),MQN(2)
C
      COMPLEX*16 ESG(MB2,MEQ),ETEMP(MB2*MRC)
C
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     IMPORT LQN AND BASIS PAIR MQNS FOR LOCAL USE
      LQNA = LQN(1)
      LQNB = LQN(2)
      MQNA = MQN(1)
      MQNB = MQN(2)
C
C     INDEX SUMMATION TERMINALS BASED ON MAX DEGREE OF HGTF
      LMAX = LQN(1)+LQN(2)
      NTUV = (LMAX+1)*(LMAX+2)*(LMAX+3)/6
C
C     INITIALISE ESG AND ARRAY
      DO M=1,MAXM
        DO ITUV=1,NTUV
          ESG(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     EXIT IF LQN,MQN ORDERS YIELD ZERO CG-COEFFICIENTS
      IF(IABS(MQN(1)).GT.LQN(1).OR.IABS(MQN(2)).GT.LQN(2)) RETURN
C
C     CHECK THAT LMAX IS WITHIN THE BOUNDS OF MKP
      IF(LMAX.GT.MKP+1) THEN
        WRITE(6, *) 'In EVRS: LMAX exceeds MKP+1 parameter.',LMAX,MKP+1
        WRITE(7, *) 'In EVRS: LMAX exceeds MKP+1 parameter.',LMAX,MKP+1
        STOP
      ENDIF
C
c     INITIALISE TEMPORARY ARRAY
      DO M=1,MB2*MRC
        ETEMP(M) = DCMPLX(0.0D0,0.0D0)
      ENDDO
C
C     SET INITIAL VALUES TO E[0,0;0,0;0,0,0,0] = RKAB
      DO M=1,MAXM
        ETEMP(M) = DCMPLX(RKAB(M),0.0D0)
      ENDDO
C
C     STEP 1:
C     GENERATE E[|MQNA|,MQNA;0,0] FROM E[0,0;0,0] USING
C     SIMULTANEOUS STEP OF LQN AND MQN ON CENTRE A
      ISTART = 0
      LAM    = 0
      CALL ESTEPLM(ETEMP,LAM,ISTART,MQNA,MAXM,1)
C
C     STEP 2:
C     GENERATE E[LQNA,MQNA;0,0] FROM E[|MQNA|,MQNA;0,0]
C     USING THE STEP OF LQN ONLY ON CENTRE A
      CALL ESTEPL(ETEMP,LAM,ISTART,LQNA,MQNA,MAXM,1)
C
C     STEP 3:
C     GENERATE E[LQNA,MQNA;|MQNB|,MQNB] FROM E[LQNA,MQNA;0,0]
C     USING SIMULTANEOUS STEP OF LQN AND MQN ON CENTRE B
      CALL ESTEPLM(ETEMP,LAM,ISTART,MQNB,MAXM,2)
C
C     STEP 4:
C     GENERATE E[LQNA,MQNA;LQNB,MQNB] FROM E[LQNA,MQNA;|MQNB|,MQNB]
C     USING THE STEP OF LQN ONLY ON CENTRE B
      CALL ESTEPL(ETEMP,LAM,ISTART,LQNB,MQNB,MAXM,2)
C
C     NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS LAM
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     STEP 5:
C     COPY FINAL BLOCK OF ENTRIES AS THE REQUIRED OUTPUT
      K = 0
      DO ITUV=1,NTUV
        DO M=1,MAXM
          K = K+1
          ESG(M,ITUV) = ETEMP(ISTART+K)
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ESTEPLM(ETEMP,LAM,ISTART,MQN,MAXM,IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C   EEEEEEEE SSSSSS TTTTTTTT EEEEEEEE PPPPPPP  LL       MM       MM    C
C   EE      SS    SS   TT    EE       PP    PP LL       MMM     MMM    C
C   EE      SS         TT    EE       PP    PP LL       MMMM   MMMM    C
C   EEEEEE   SSSSSS    TT    EEEEEE   PP    PP LL       MM MM MM MM    C
C   EE            SS   TT    EE       PPPPPPP  LL       MM  MMM  MM    C
C   EE      SS    SS   TT    EE       PP       LL       MM   M   MM    C
C   EEEEEEEE SSSSSS    TT    EEEEEEEE PP       LLLLLLLL MM       MM    C
C                                                                      C
C -------------------------------------------------------------------- C
C  SIMULTANEOUSLY INCREMENT THE QUANTUM NUMBERS LQN & MQN, STARTING    C
C  WITH E[0,0], USING THE RECURSION ALGORITHM OF V.R.SAUNDERS IN       C
C  `MOLECULAR INTEGRALS FOR GAUSSIAN TYPE FUNCTIONS' 1983.             C
C -------------------------------------------------------------------- C
C       E[0,0;IT,IU,IV] -> E[|MQN|+1,±(|MQN|+1);IT,IU,IV]  Eq.(64)     C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ LAM    - LENGTH OF THE INPUT HGTF EXPANSION.                      C
C  ▶ ISTART - COUNTER TO TRACK THE START OF THIS PORTION OF THE LIST.  C
C  ▶ MQN    - MAGNETIC QUANTUM NUMBER.                                 C
C  ▶ MAXM   - NUMBER OF EXPONENT/DENSITY PAIRS.                        C
C  ▶ IZ     - CENTRE TO STEP UP.                                       C
C  OUTPUT:                                                             C
C  ▶ ETEMP  - UN-NORMALISED EXPANSION COEFFICIENTS FOR THIS BLOCK.     C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION PX(MB2),PY(MB2)
C
      COMPLEX*16 CONE
      COMPLEX*16 ETEMP(MB2*MRC)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     IF |MQN|.EQ.0 THEN NO INCREMENT IN (LQN,MQN) IS REQUIRED
      IF(IABS(MQN).EQ.0) RETURN
C
C     IMPORT GEOMETRIC VALUES FOR CENTRE OF INTEREST
      DO M=1,MAXM
        IF(IZ.EQ.1) THEN
          PX(M) = PAX(M)
          PY(M) = PAY(M)
        ELSEIF(IZ.EQ.2) THEN
          PX(M) = PBX(M)
          PY(M) = PBY(M)
        ENDIF
      ENDDO
C
C     PHASE TERM FOR SIGN OF MQN
      PHM = DFLOAT(ISIGN(1,MQN))
C
C     LOOP OVER ALL MAGNETIC NUMBERS UP TO THIS MQN (USE LQN AS COUNTER)
      DO LQN=0,IABS(MQN)-1
C
C       DEGENERACY COUNTER 2*LQN+1
        R2L1 = DFLOAT(2*LQN+1)
C
C       NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS LAM
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C       ISTART LABELS THE PREVIOUS LQN VALUE
C       JSTART LABELS THE CURRENT  LQN VALUE
C       KSTART LABELS THE NEXT     LQN VALUE
        JSTART = ISTART
        KSTART = JSTART + NTUV*MAXM
C
C**********************************************************************C
C                          INDEX MAPPINGS:                             C
C -------------------------------------------------------------------- C
C       J0-> E[LQN  ,LQN  ;IT  ,IU  ,IV]                               C
C       K0-> E[LQN+1,LQN+1;IT  ,IU  ,IV]                               C
C       K1-> E[LQN+1,LQN+1;IT+1,IU  ,IV]                               C
C       K2-> E[LQN+1,LQN+1;IT  ,IU+1,IV]                               C
C       K3-> E[LQN+1,LQN+1;IT-1,IU  ,IV]                               C
C       K4-> E[LQN+1,LQN+1;IT  ,IU-1,IV]                               C
C**********************************************************************C
C
C       LOOP OVER THE HGTF INDICES OF THE SEED LAYER
        DO IOUTER=0,LAM
          DO IT=0,IOUTER
            DO IU=0,IOUTER-IT
              IV = IOUTER-IT-IU
C
C             STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (LQN)
              J0 = JSTART + (IABC(IT  ,IU  ,IV)-1)*MAXM
              K0 = KSTART + (IABC(IT  ,IU  ,IV)-1)*MAXM
              K1 = KSTART + (IABC(IT+1,IU  ,IV)-1)*MAXM
              K2 = KSTART + (IABC(IT  ,IU+1,IV)-1)*MAXM
              IF(IT.NE.0) K3 = KSTART + (IABC(IT-1,IU  ,IV  )-1)*MAXM
              IF(IU.NE.0) K4 = KSTART + (IABC(IT  ,IU-1,IV  )-1)*MAXM
C
C             INVOKE RECURRENCE RELATIONS ON LAYER (LQN+1)
C
C             TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
              DO M=1,MAXM
                T1 = R2L1/P2(M)
                TX = R2L1*PX(M)
                TY = R2L1*PY(M)
                ETEMP(K0+M) = ETEMP(K0+M) +          TX*ETEMP(J0+M)
     &                                    + PHM*CONE*TY*ETEMP(J0+M)
                ETEMP(K1+M) = ETEMP(K1+M) +          T1*ETEMP(J0+M)
                ETEMP(K2+M) = ETEMP(K2+M) + PHM*CONE*T1*ETEMP(J0+M)
              ENDDO
C
C             SPECIAL CASE EXCLUDES IT=0
              IF(IT.NE.0) THEN
                FC = R2L1*DFLOAT(IT)
                DO M=1,MAXM
                  ETEMP(K3+M) = ETEMP(K3+M) +          FC*ETEMP(J0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IU=0
              IF(IU.NE.0) THEN
                FC = R2L1*DFLOAT(IU)
                DO M=1,MAXM
                  ETEMP(K4+M) = ETEMP(K4+M) + PHM*CONE*FC*ETEMP(J0+M)
                ENDDO
              ENDIF
C
C           END OF LOOPS OVER HGTF INDICES
            ENDDO
          ENDDO
        ENDDO
C
C       UPDATE 'PREVIOUS' START VALUE
        ISTART = ISTART + NTUV*MAXM
C
C       UPDATE LAMBDA VALUE
        LAM = LAM+1
C
C     END OF LOOP OVER MQN COUNTER
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ESTEPL(ETEMP,LAM,ISTART0,LQN,MQN,MAXM,IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         EEEEEEEE SSSSSS TTTTTTTT EEEEEEEE PPPPPPP  LL                C
C         EE      SS    SS   TT    EE       PP    PP LL                C
C         EE      SS         TT    EE       PP    PP LL                C
C         EEEEEE   SSSSSS    TT    EEEEEE   PP    PP LL                C
C         EE            SS   TT    EE       PPPPPPP  LL                C
C         EE      SS    SS   TT    EE       PP       LL                C
C         EEEEEEEE SSSSSS    TT    EEEEEEEE PP       LLLLLLLL          C
C                                                                      C
C -------------------------------------------------------------------- C
C  DECREMENT THE LQN, STARTING WITH E[MQN,±|MQN|], USING THE RECURSION C
C  ALGORITHM OF V.R.SAUNDERS IN `MOLECULAR INTEGRALS FOR GAUSSIAN TYPE C
C  FUNCTIONS' 1983 (EDITED BY G.H.F. DIERCKSEN AND S. WILSON).         C
C -------------------------------------------------------------------- C
C (1)  E[|MQN|,±MQN;IT,IU,IV] -> E[LQN+1,±MQN;IT,IU,IV]    Eq.(64)     C
C (2)  E[LQN  ,±MQN;IT,IU,IV] -> E[LQN+1,±MQN;IT,IU,IV]                C
C (3)  E[LQN-1,±MQN;IT,IU,IV] -> E[LQN+1,±MQN;IT,IU,IV]                C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ LAM     - LENGTH OF THE INPUT HGTF EXPANSION.                     C
C  ▶ ISTART0 - COUNTER TO TRACK THE START OF THIS PORTION OF THE LIST. C
C  ▶ LQN     - ORBITAL QUANTUM NUMBER.                                 C
C  ▶ MQN     - MAGNETIC QUANTUM NUMBER.                                C
C  ▶ MAXM    - NUMBER OF EXPONENT/DENSITY PAIRS.                       C
C  ▶ IZ      - CENTRE TO STEP UP.                                      C
C  OUTPUT:                                                             C
C  ▶ ETEMP   - UN-NORMALISED EXPANSION COEFFICIENTS FOR THIS BLOCK.    C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION PP(MB2),PX(MB2),PY(MB2),PZ(MB2)
C
      COMPLEX*16 CONE
      COMPLEX*16 ETEMP(MB2*MRC)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     DEFINE THE UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     IF LQN.LE.|MQN| THEN NO INCREMENT IN LQN IS REQUIRED
      IF(LQN.LE.IABS(MQN)) RETURN
C
C     IMPORT GEOMETRIC VALUES FOR CENTRE OF INTEREST
      DO M=1,MAXM
        IF(IZ.EQ.1) THEN
          PP(M) = PA2(M)
          PX(M) = PAX(M)
          PY(M) = PAY(M)
          PZ(M) = PAZ(M)
        ELSEIF(IZ.EQ.2) THEN
          PP(M) = PB2(M)
          PX(M) = PBX(M)
          PY(M) = PBY(M)
          PZ(M) = PBZ(M)
        ENDIF
      ENDDO
C
C     NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS LAM
      NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C**********************************************************************C
C     STEP (1): E[|MQN|,MQN;IT,IU,IV] -> E[|MQN|+1,MQN;IT,IU,IV].      C
C               IT MAPS SOME INDEX SETS FROM DATA OBTAINED IN ESTEPLM. C
C -------------------------------------------------------------------- C
C     INDEX MAPPINGS: J0-> E[MQN  ,MQN;IT  ,IU  ,IV  ]                 C
C                     K0-> E[MQN+1,MQN;IT  ,IU  ,IV  ]                 C
C                     K6-> E[MQN+1,MQN;IT  ,IU  ,IV+1]                 C
C                     K9-> E[MQN+1,MQN;IT  ,IU  ,IV-1]                 C
C**********************************************************************C
C
C     ISTART0 LABELS THE GLOBAL STARTING VALUE
C     ISTART  LABELS THE PREVIOUS LQN VALUE
C     JSTART  LABELS THE CURRENT  LQN VALUE
C     KSTART  LABELS THE NEXT     LQN VALUE
      JSTART = ISTART0
      KSTART = JSTART + NTUV*MAXM
C
C     OVERALL LQN/MQN FACTOR SIMPLIFIES WHEN LQN = |MQN|
      RLM  = DFLOAT(2*IABS(MQN)+1)
C
C     LOOP OVER THE HGTF INDICES OF THE SEED LAYER
      DO IOUTER=0,LAM
        DO IT=0,IOUTER
          DO IU=0,IOUTER-IT
            IV = IOUTER-IT-IU
C
C           STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (LQN+1)
            J0 = JSTART + (IABC(IT  ,IU  ,IV  )-1)*MAXM
            K0 = KSTART + (IABC(IT  ,IU  ,IV  )-1)*MAXM
            K6 = KSTART + (IABC(IT  ,IU  ,IV+1)-1)*MAXM
            IF(IV.NE.0) K9 = KSTART + (IABC(IT  ,IU  ,IV-1)-1)*MAXM
C
C           INVOKE RECURRENCE RELATIONS ON LAYER (LQN+1)
C
C           TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
            DO M=1,MAXM
              TZ = RLM*PZ(M)
              TP = RLM/P2(M)
              ETEMP(K0+M) = ETEMP(K0+M) + TZ*ETEMP(J0+M)
              ETEMP(K6+M) = ETEMP(K6+M) + TP*ETEMP(J0+M)
            ENDDO
C
C           SPECIAL CASE EXCLUDES IV=0
            IF(IV.GE.1) THEN
              FC = RLM*DFLOAT(IV)
              DO M=1,MAXM
                ETEMP(K9+M) = ETEMP(K9+M) + FC*ETEMP(J0+M)
              ENDDO
            ENDIF
C
          ENDDO
        ENDDO
      ENDDO
C
C     UPDATE BLOCK LOCATORS
      ISTART  = ISTART0
      ISTART0 = ISTART0 + NTUV*MAXM
C
C     UPDATE LAMBDA VALUE
      LAM = LAM+1
C
C     IF LQN=|MQN|+1 THEN SET IS FINISHED
      IF(LQN.EQ.IABS(MQN)+1) RETURN
C
C**********************************************************************C
C   THIS AND SUBSEQUENT STEPS IN RECURRENCE INVOLVE THREE LAYERS:      C
C            E[LQN  ,MQN;IT,IU,IV] -> E[LQN+1,MQN;IT,IU,IV]            C
C            E[LQN-1,MQN;IT,IU,IV] -> E[LQN+1,MQN;IT,IU,IV]            C
C -------------------------------------------------------------------- C
C   INDEX MAPPINGS: I0 -> E[LQN-1,MQN;IT  ,IU  ,IV  ]                  C
C                   J0 -> E[LQN  ,MQN;IT  ,IU  ,IV  ]                  C
C                   K0 -> E[LQN+1,MQN;IT  ,IU  ,IV  ]                  C
C                   K1 -> E[LQN+1,MQN;IT+2,IU  ,IV  ]                  C
C                   K2 -> E[LQN+1,MQN;IT  ,IU+2,IV  ]                  C
C                   K3 -> E[LQN+1,MQN;IT  ,IU  ,IV+2]                  C
C                   K4 -> E[LQN+1,MQN;IT+1,IU  ,IV  ]                  C
C                   K5 -> E[LQN+1,MQN;IT  ,IU+1,IV  ]                  C
C                   K6 -> E[LQN+1,MQN;IT  ,IU  ,IV+1]                  C
C                   K7 -> E[LQN+1,MQN;IT-1,IU  ,IV  ]                  C
C                   K8 -> E[LQN+1,MQN;IT  ,IU-1,IV  ]                  C
C                   K9 -> E[LQN+1,MQN;IT  ,IU  ,IV-1]                  C
C                   K10-> E[LQN+1,MQN;IT-2,IU  ,IV  ]                  C
C                   K11-> E[LQN+1,MQN;IT  ,IU-2,IV  ]                  C
C                   K12-> E[LQN+1,MQN;IT  ,IU  ,IV-2]                  C
C**********************************************************************C
C
C     LOOP OVER LSTP FOR FIXED MQN DOWN TO LQN-1
      DO LSTP=IABS(MQN)+1,LQN-1

C       OVERALL LSTP/MQN FACTORS
        RLM1 = DFLOAT(2*LSTP+1)/DFLOAT(LSTP-IABS(MQN)+1)
        RLM2 =-DFLOAT(LSTP+IABS(MQN))/DBLE(LSTP-IABS(MQN)+1)
C
C       NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS LAM
        NTUV = (LAM+1)*(LAM+2)*(LAM+3)/6
C
        JSTART = ISTART0
        KSTART = ISTART0 + MAXM*NTUV
C
C**********************************************************************C
C     STEP (2): E[LSTP  ,MQN;IT,IU,IV] -> E[LSTP+1,MQN;IT,IU,IV].      C
C**********************************************************************C
C
C       LOOP OVER THE HGTF INDICES OF THE SEED LAYER
        DO IOUTER=0,LAM
          DO IT=0,IOUTER
            DO IU=0,IOUTER-IT
              IV = IOUTER-IT-IU
C
C             STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (LSTP)
              J0 = JSTART + (IABC(IT  ,IU  ,IV  )-1)*MAXM
              K0 = KSTART + (IABC(IT  ,IU  ,IV  )-1)*MAXM
              K6 = KSTART + (IABC(IT  ,IU  ,IV+1)-1)*MAXM
              IF(IV.NE.0) K9 = KSTART + (IABC(IT  ,IU  ,IV-1)-1)*MAXM
C
C             INVOKE RECURRENCE RELATIONS ON LAYER (LSTP)
C
C             TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
              DO M=1,MAXM
                TZ = PZ(M)
                TP = 1.0D0/P2(M)
                ETEMP(K0+M) = ETEMP(K0+M) + TZ*RLM1*ETEMP(J0+M)
                ETEMP(K6+M) = ETEMP(K6+M) + TP*RLM1*ETEMP(J0+M)
              ENDDO
C
C             SPECIAL CASE EXCLUDES IV=0
              IF(IV.GE.1) THEN
                FC = DFLOAT(IV)
                DO M=1,MAXM
                  ETEMP(K9+M) = ETEMP(K9+M) + FC*RLM1*ETEMP(J0+M)
                ENDDO
              ENDIF
C
            ENDDO
          ENDDO
        ENDDO
C
C**********************************************************************C
C     STEP (3): E[LSTP-1,MQN;IT,IU,IV] -> E[LSTP+1,MQN;IT,IU,IV].      C
C**********************************************************************C
C
C       LOOP OVER THE HGTF INDICES OF THE SEED LAYER
        DO IOUTER=0,LAM-1
          DO IT=0,IOUTER
            DO IU=0,IOUTER-IT
              IV = IOUTER-IT-IU
C
C             STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (LSTP)
              I0  = ISTART+(IABC(IT  ,IU  ,IV  )-1)*MAXM
              K0  = KSTART+(IABC(IT  ,IU  ,IV  )-1)*MAXM
              K1  = KSTART+(IABC(IT+2,IU  ,IV  )-1)*MAXM
              K2  = KSTART+(IABC(IT  ,IU+2,IV  )-1)*MAXM
              K3  = KSTART+(IABC(IT  ,IU  ,IV+2)-1)*MAXM
              K4  = KSTART+(IABC(IT+1,IU  ,IV  )-1)*MAXM
              K5  = KSTART+(IABC(IT  ,IU+1,IV  )-1)*MAXM
              K6  = KSTART+(IABC(IT  ,IU  ,IV+1)-1)*MAXM
              IF(IT.GT.0) K7  = KSTART + (IABC(IT-1,IU  ,IV  )-1)*MAXM
              IF(IU.GT.0) K8  = KSTART + (IABC(IT  ,IU-1,IV  )-1)*MAXM
              IF(IV.GT.0) K9  = KSTART + (IABC(IT  ,IU  ,IV-1)-1)*MAXM
              IF(IT.GT.1) K10 = KSTART + (IABC(IT-2,IU  ,IV  )-1)*MAXM
              IF(IU.GT.1) K11 = KSTART + (IABC(IT  ,IU-2,IV  )-1)*MAXM
              IF(IV.GT.1) K12 = KSTART + (IABC(IT  ,IU  ,IV-2)-1)*MAXM
C
C             INVOKE RECURRENCE RELATIONS ON LAYER (LSTP-1)
C
C             TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
              TI = DFLOAT(2*(IT+IU+IV)+3)
              DO M=1,MAXM
                T1 = 1.0D0/P22(M)
                TX = PX(M)/P(M)
                TY = PY(M)/P(M)
                TZ = PZ(M)/P(M)
                TT = PP(M) + TI/P2(M)
                ETEMP(K0+M) = ETEMP(K0+M) + RLM2*TT*ETEMP(I0+M)
                ETEMP(K1+M) = ETEMP(K1+M) + RLM2*T1*ETEMP(I0+M)
                ETEMP(K2+M) = ETEMP(K2+M) + RLM2*T1*ETEMP(I0+M)
                ETEMP(K3+M) = ETEMP(K3+M) + RLM2*T1*ETEMP(I0+M)
                ETEMP(K4+M) = ETEMP(K4+M) + RLM2*TX*ETEMP(I0+M)
                ETEMP(K5+M) = ETEMP(K5+M) + RLM2*TY*ETEMP(I0+M)
                ETEMP(K6+M) = ETEMP(K6+M) + RLM2*TZ*ETEMP(I0+M)
              ENDDO
C
C             SPECIAL CASE EXCLUDES IT=0
              IF(IT.GE.1) THEN
                DO M=1,MAXM
                  TX = DFLOAT(2*IT)*PX(M)
                  ETEMP(K7+M) = ETEMP(K7+M) + RLM2*TX*ETEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IU=0
              IF(IU.GE.1) THEN
                DO M=1,MAXM
                  TY = DFLOAT(2*IU)*PY(M)
                  ETEMP(K8+M) = ETEMP(K8+M) + RLM2*TY*ETEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IV=0
              IF(IV.GE.1) THEN
                DO M=1,MAXM
                  TZ = DFLOAT(2*IV)*PZ(M)
                  ETEMP(K9+M) = ETEMP(K9+M) + RLM2*TZ*ETEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IT=0,1
              IF(IT.GE.2) THEN
                TX = DFLOAT(IT*(IT-1))
                DO M=1,MAXM
                  ETEMP(K10+M) = ETEMP(K10+M) + RLM2*TX*ETEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IU=0,1
              IF(IU.GE.2) THEN
                TY = DFLOAT(IU*(IU-1))
                DO M=1,MAXM
                  ETEMP(K11+M) = ETEMP(K11+M) + RLM2*TY*ETEMP(I0+M)
                ENDDO
              ENDIF
C
C             SPECIAL CASE EXCLUDES IV=0,1
              IF(IV.GE.2) THEN
                TZ = DFLOAT(IV*(IV-1))
                DO M=1,MAXM
                  ETEMP(K12+M) = ETEMP(K12+M) + RLM2*TZ*ETEMP(I0+M)
                ENDDO
              ENDIF
C
            ENDDO
          ENDDO
        ENDDO
C
C       UPDATE BLOCK LOCATORS
        ISTART  = ISTART0
        ISTART0 = ISTART0 + MAXM*NTUV
C
C       UPDATE LAMBDA VALUE
        LAM = LAM+1
C
C     END OF LOOP OVER LSTP FOR FIXED MQN
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE ESTEPN(ESG,ENSG,LAM,MAXM,IZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         EEEEEEEE SSSSSS TTTTTTTT EEEEEEEE PPPPPPP  NN    NN          C
C         EE      SS    SS   TT    EE       PP    PP NNN   NN          C
C         EE      SS         TT    EE       PP    PP NNNN  NN          C
C         EEEEEE   SSSSSS    TT    EEEEEE   PP    PP NN NN NN          C
C         EE            SS   TT    EE       PPPPPPP  NN  NNNN          C
C         EE      SS    SS   TT    EE       PP       NN   NNN          C
C         EEEEEEEE SSSSSS    TT    EEEEEEEE PP       NN    NN          C
C                                                                      C
C -------------------------------------------------------------------- C
C  INCREMENT THE QUANTUM NUMBER (NQN):                                 C
C                E[NQN  ,LQN,MQN] -> E[NQN+1,LQN,MQN].                 C
C -------------------------------------------------------------------- C
C  ▶ ONLY PERFORMS A SINGLE STEP IN NQN.                               C
C  ▶ LAM IS THE EFFECTIVE TOTAL ANGULAR MOMENTUM OF THE INPUT COEFFS.  C
C    EFFECTIVE TOTAL ANGULAR MOMENTUM OF THE OUTPUT COEFFS IS LAM+2.   C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ ESG  - EQ-COEFFICIENT BATCH.                                      C
C  ▶ LAM  - EFFECTIVE TOTAL ANGULAR MOMENTUM.                          C
C  ▶ MAXM - NUMBER OF EXPONENT/DENSITY PAIRS.                          C
C  ▶ IZ   - CENTRE TO STEP UP.                                         C
C  OUTPUT:                                                             C
C  ▶ ENSG - EQ-COEFFICIENT BATCH AFTER NQN HAS BEEN STEPPED UP.        C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      DIMENSION PP(MB2),PX(MB2),PY(MB2),PZ(MB2)
C
      COMPLEX*16 ESG(MB2,MEQ),ENSG(MB2,MEQ)
C
      COMMON/ICRT/IABC(0:ML4,0:ML4,0:ML4),IA(MRC),IB(MRC),IC(MRC),
     &            ILAM(MRC)
      COMMON/MCMD/P(MB2),P2(MB2),P22(MB2),UV(MB2),RKAB(MB2),
     &            PAX(MB2),PAY(MB2),PAZ(MB2),PBX(MB2),PBY(MB2),PBZ(MB2),
     &            PA2(MB2),PB2(MB2),EIBS(MB2),EJBS(MB2),ABDST(3)
C
C     NUMBER OF TERMS IN CARTESIAN EXPANSION FOR THIS LAM WITH N'=N+1
      NTUV = (LAM+3)*(LAM+4)*(LAM+5)/6
C
C     INITIALISE NEW ARRAY
      DO ITUV=1,NTUV
        DO M=1,MAXM
          ENSG(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     IMPORT GEOMETRIC VALUES FOR CENTRE OF INTEREST
      DO M=1,MAXM
        IF(IZ.EQ.1) THEN
          PP(M) = PA2(M)
          PX(M) = PAX(M)
          PY(M) = PAY(M)
          PZ(M) = PAZ(M)
        ELSEIF(IZ.EQ.2) THEN
          PP(M) = PB2(M)
          PX(M) = PBX(M)
          PY(M) = PBY(M)
          PZ(M) = PBZ(M)
        ENDIF
      ENDDO
C
C**********************************************************************C
C     INDEX MAPPINGS: I0 -> E[LQN-1,MQN;IT  ,IU  ,IV  ]                C
C                     K0 -> E[LQN+1,MQN;IT  ,IU  ,IV  ]                C
C                     K1 -> E[LQN+1,MQN;IT+2,IU  ,IV  ]                C
C                     K2 -> E[LQN+1,MQN;IT  ,IU+2,IV  ]                C
C                     K3 -> E[LQN+1,MQN;IT  ,IU  ,IV+2]                C
C                     K4 -> E[LQN+1,MQN;IT+1,IU  ,IV  ]                C
C                     K5 -> E[LQN+1,MQN;IT  ,IU+1,IV  ]                C
C                     K6 -> E[LQN+1,MQN;IT  ,IU  ,IV+1]                C
C                     K7 -> E[LQN+1,MQN;IT-1,IU  ,IV  ]                C
C                     K8 -> E[LQN+1,MQN;IT  ,IU-1,IV  ]                C
C                     K9 -> E[LQN+1,MQN;IT  ,IU  ,IV-1]                C
C                     K10-> E[LQN+1,MQN;IT-2,IU  ,IV  ]                C
C                     K11-> E[LQN+1,MQN;IT  ,IU-2,IV  ]                C
C                     K12-> E[LQN+1,MQN;IT  ,IU  ,IV-2]                C
C**********************************************************************C
C
      DO IOUTER=0,LAM
        DO IT=0,IOUTER
          DO IU=0,IOUTER-IT
            IV = IOUTER-IT-IU
C
C           STARTING LOCATION FOR GIVEN (IA,IB,IC) AND EFF. (NQN)
            I0  = IABC(IT  ,IU  ,IV  )
            K0  = IABC(IT  ,IU  ,IV  )
            K1  = IABC(IT+2,IU  ,IV  )
            K2  = IABC(IT  ,IU+2,IV  )
            K3  = IABC(IT  ,IU  ,IV+2)
            K4  = IABC(IT+1,IU  ,IV  )
            K5  = IABC(IT  ,IU+1,IV  )
            K6  = IABC(IT  ,IU  ,IV+1)
            IF(IT.GT.0) K7  = IABC(IT-1,IU  ,IV  )
            IF(IU.GT.0) K8  = IABC(IT  ,IU-1,IV  )
            IF(IV.GT.0) K9  = IABC(IT  ,IU  ,IV-1)
            IF(IT.GT.1) K10 = IABC(IT-2,IU  ,IV  )
            IF(IU.GT.1) K11 = IABC(IT  ,IU-2,IV  )
            IF(IV.GT.1) K12 = IABC(IT  ,IU  ,IV-2)
C
C           INVOKE RECURRENCE RELATIONS ON LAYER (NQN+1)
C
C           TERMS THAT ALWAYS APPLY (NO SPECIAL CASES)
            TT = DFLOAT((2*(IT+IU+IV))+3)
            DO M=1,MAXM
              T1 = 1.0D0/P22(M)
              TX = PX(M)/P(M)
              TY = PY(M)/P(M)
              TZ = PZ(M)/P(M)
              TP = PP(M) + TT/P2(M)
              ENSG(M,K0) = ENSG(M,K0) + TP*ESG(M,I0)
              ENSG(M,K1) = ENSG(M,K1) + T1*ESG(M,I0)
              ENSG(M,K2) = ENSG(M,K2) + T1*ESG(M,I0)
              ENSG(M,K3) = ENSG(M,K3) + T1*ESG(M,I0)
              ENSG(M,K4) = ENSG(M,K4) + TX*ESG(M,I0)
              ENSG(M,K5) = ENSG(M,K5) + TY*ESG(M,I0)
              ENSG(M,K6) = ENSG(M,K6) + TZ*ESG(M,I0)
            ENDDO
C
C           SPECIAL CASE EXCLUDES IT=0
            IF(IT.GE.1) THEN
              RT2 = DFLOAT(2*IT)
              DO M=1,MAXM
                T0 = PX(M)*RT2
                ENSG(M,K7) = ENSG(M,K7) + T0*ESG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IU=0
            IF(IU.GE.1) THEN
              RU1 = DFLOAT(2*IU)
              DO M=1,MAXM
                T0 = PY(M)*RU1
                ENSG(M,K8) = ENSG(M,K8) + T0*ESG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IV=0
            IF(IV.GE.1) THEN
              RV1 = DFLOAT(2*IV)
              DO M=1,MAXM
                T0 = PZ(M)*RV1
                ENSG(M,K9) = ENSG(M,K9) + T0*ESG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IT=0,1
            IF(IT.GE.2) THEN
              RT2 = DFLOAT(IT*(IT-1))
              DO M=1,MAXM
                ENSG(M,K10) = ENSG(M,K10) + RT2*ESG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IU=0,1
            IF(IU.GE.2) THEN
              RU2 = DFLOAT(IU*(IU-1))
              DO M=1,MAXM
                ENSG(M,K11) = ENSG(M,K11) + RU2*ESG(M,I0)
              ENDDO
            ENDIF
C
C           SPECIAL CASE EXCLUDES IV=0,1
            IF(IV.GE.2) THEN
              RV2 = DFLOAT(IV*(IV-1))
              DO M=1,MAXM
                ENSG(M,K12) = ENSG(M,K12) + RV2*ESG(M,I0)
              ENDDO
            ENDIF

          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END

