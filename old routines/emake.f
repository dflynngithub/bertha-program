      SUBROUTINE EMAKELL(ELL11,ELL21,EXPT,KQN,MQN,NFUNS,IALT,I1,I2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C    EEEEEEEE MM     MM    AA    KK    KK EEEEEEEE LL      LL         C
C    EE       MMM   MMM   AAAA   KK   KK  EE       LL      LL         C
C    EE       MMMM MMMM  AA  AA  KK  KK   EE       LL      LL         C
C    EEEEEE   MM MMM MM AA    AA KKKKK    EEEEEE   LL      LL         C
C    EE       MM  M  MM AAAAAAAA KK  KK   EE       LL      LL         C
C    EE       MM     MM AA    AA KK   KK  EE       LL      LL         C
C    EEEEEEEE MM     MM AA    AA KK    KK EEEEEEEE LLLLLLL LLLLLLL    C
C                                                                     C
C ------------------------------------------------------------------- C
C     EMAKELL GENERATES BLOCKS OF SPHERICAL SPINOR E-COEFFICIENTS     C
C     BY CONTRACTING ON THE E-COEFFICIENTS FOR SCALAR SPHERICAL       C
C     GAUSSIAN FUNCTIONS, USING A DEVELOPMENT OF THE ALGORITHM OF     C
C     V.R.SAUNDERS                                                    C
C ------------------------------------------------------------------- C
C     INPUT --                                                        C
C      I1: INDEX 1 (OF FOUR OPTIONS)                                  C
C      I2: INDEX 2 (OF FOUR OPTIONS)                                  C
C ------------------------------------------------------------------- C
C     H.M.QUINEY THE UNIVERSITY OF MELBOURNE  (2008)                  C 
C*********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=MKP,MBS=26,MB2=MBS*MBS,
     &          MLL=(MKP+8)*(MKP+9)*(MKP+10)/6)
C
      COMPLEX*16 ELL11(MB2,MLL),ELL21(MB2,MLL)
C
      DIMENSION EXPT(MBS,4),KQN(4),MQN(4),NFUNS(4)
      DIMENSION EXL(MBS,2),KQNLAB(2),MQNLAB(2),NFLAB(2),COORD(3,2)
C
      COMMON/XYZ4/XYZ(3,4)
C
C     DETERMINE LQNS FROM KQNS
      IF(KQN(I1).GT.0) THEN
        LA = KQN(I1)
      ELSE
        LA =-KQN(I1)-1
      ENDIF
      IF(KQN(I2).GT.0) THEN
        LB = KQN(I2)
      ELSE
        LB =-KQN(I2)-1
      ENDIF
C
C     SUMMATION TERMINATES AT LAMBDA = LA + LB FOR LL PAIRS
      LAMAB = LA + LB
C
C     KAPPA LABELS AND NUMBER OF FUNCTIONS ON A GIVEN CENTRE
      KQNLAB(1) = KQN(I1)
      KQNLAB(2) = KQN(I2)
      NFLAB(1)  = NFUNS(I1)
      NFLAB(2)  = NFUNS(I2)
C
C     CARTESIAN COORDINATES OF EACH CENTRE
      DO I=1,3
        COORD(I,1) = XYZ(I,I1)
        COORD(I,2) = XYZ(I,I2)
      ENDDO
C
C     BASIS SET EXPONENTS FOR CENTRE A
      DO IBAS=1,NFUNS(I1)
       EXL(IBAS,1) = EXPT(IBAS,I1)       
      ENDDO
C
C     BASIS SET EXPONENTS FOR CENTRE B
      DO JBAS=1,NFUNS(I2)
       EXL(JBAS,2) = EXPT(JBAS,I2)
      ENDDO                  
C
C*********************************************************************C
C     1: GENERATE AND STORE ELL11, FOR M PAIRS (-|MA|,-|MB|)          C
C*********************************************************************C
C
C     M QUANTUM NUMBERS
      MQNLAB(1) =-MQN(I1)
      MQNLAB(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLL(ELL11,EXL,COORD,KQNLAB,MQNLAB,NFLAB,IQ)
C
C     MULTIPLY SOME OF THE ELL0 COEFFICIENTS BY A PHASE TERM
      ITUV=0
      DO MTUV=0,LAMAB
        PHASE = DFLOAT((IALT)**(MTUV))
        DO MDUM=1,((MTUV+1)*(MTUV+2))/2
          ITUV = ITUV+1
          DO M=1,NFUNS(I1)*NFUNS(I2)
            ELL11(M,ITUV) = ELL11(M,ITUV)*PHASE
          ENDDO
        ENDDO
      ENDDO
C
C*********************************************************************C
C     2: GENERATE AND STORE ELL21, FOR M PAIRS (+|MA|,-|MB|)          C
C*********************************************************************C
C
C     M QUANTUM NUMBERS
      MQNLAB(1) = MQN(I1)
      MQNLAB(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLL(ELL21,EXL,COORD,KQNLAB,MQNLAB,NFLAB,IQ)
C
C     MULTIPLY SOME OF THE ELL0 COEFFICIENTS BY A PHASE TERM
      ITUV=0
      DO MTUV=0,LAMAB
        PHASE = DFLOAT((IALT)**(MTUV))
        DO MDUM=1,((MTUV+1)*(MTUV+2))/2
          ITUV = ITUV+1
          DO M=1,NFUNS(I1)*NFUNS(I2)
            ELL21(M,ITUV) = ELL21(M,ITUV)*PHASE
          ENDDO
        ENDDO
      ENDDO
C
C     NOTE THAT ELL22 AND ELL12 ARE RELATED TO THESE BY PHASE FACTORS.
C
      RETURN
      END
C
C
      SUBROUTINE EMAKESS(ESS11,ESS21,EXPT,KQN,MQN,NFUNS,IALT,I1,I2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C   EEEEEEEE MM     MM    AA    KK    KK EEEEEEEE  SSSSSS   SSSSSS    C
C   EE       MMM   MMM   AAAA   KK   KK  EE       SS    SS SS    SS   C
C   EE       MMMM MMMM  AA  AA  KK  KK   EE       SS       SS         C
C   EEEEEE   MM MMM MM AA    AA KKKKK    EEEEEE    SSSSSS   SSSSSS    C
C   EE       MM  M  MM AAAAAAAA KK  KK   EE             SS       SS   C
C   EE       MM     MM AA    AA KK   KK  EE       SS    SS SS    SS   C
C   EEEEEEEE MM     MM AA    AA KK    KK EEEEEEEE  SSSSSS   SSSSSS    C
C                                                                     C
C ------------------------------------------------------------------- C
C     EMAKESS GENERATES BLOCKS OF SPHERICAL SPINOR E-COEFFICIENTS     C
C     BY CONTRACTING ON THE E-COEFFICIENTS FOR SCALAR SPHERICAL       C
C     GAUSSIAN FUNCTIONS, USING A DEVELOPMENT OF THE ALGORITHM OF     C
C     V.R.SAUNDERS                                                    C
C ------------------------------------------------------------------- C
C     INPUT --                                                        C
C      I1: INDEX 1 (OF FOUR OPTIONS)                                  C
C      I2: INDEX 2 (OF FOUR OPTIONS)                                  C
C ------------------------------------------------------------------- C
C     H.M.QUINEY THE UNIVERSITY OF MELBOURNE  (2008)                  C 
C*********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=MKP,MBS=26,MB2=MBS*MBS,
     &          MLL=(MKP+8)*(MKP+9)*(MKP+10)/6)
C
      COMMON/XYZ4/XYZ(3,4)
C
      DIMENSION EXPT(MBS,4),KQN(4),MQN(4),NFUNS(4)
      DIMENSION EXL(MBS,2),KQNLAB(2),MQNLAB(2),NFLAB(2),COORD(3,2)
C
      COMPLEX*16 ESS11(MB2,MLL),ESS21(MB2,MLL)
C
C     DETERMINE LQNS FROM KQNS
      IF(KQN(I1).GT.0) THEN
        LA = KQN(I1)
      ELSE
        LA =-KQN(I1)-1
      ENDIF
      IF(KQN(I2).GT.0) THEN
        LB = KQN(I2)
      ELSE
        LB =-KQN(I2)-1
      ENDIF
C
C     SUMMATION TERMINATES AT LAMBDA = LA + LB + 2 FOR SS PAIRS
      LAMAB = LA + LB + 2
C
C     KAPPA LABELS AND NUMBER OF FUNCTIONS ON A GIVEN CENTRE
      KQNLAB(1) = KQN(I1)
      KQNLAB(2) = KQN(I2)
      NFLAB(1)  = NFUNS(I1)
      NFLAB(2)  = NFUNS(I2)
C
C     CARTESIAN COORDINATES OF EACH CENTRE
      DO I=1,3
        COORD(I,1) = XYZ(I,I1)
        COORD(I,2) = XYZ(I,I2)
      ENDDO
C
C     BASIS SET EXPONENTS FOR CENTRE A
      DO IBAS=1,NFUNS(I1)
        EXL(IBAS,1) = EXPT(IBAS,I1)
      ENDDO
C
C     BASIS SET EXPONENTS FOR CENTRE B
      DO JBAS=1,NFUNS(I2)
        EXL(JBAS,2) = EXPT(JBAS,I2)
      ENDDO      
C
C*********************************************************************C
C     1: GENERATE AND STORE ESS11, FOR M PAIRS (-|MA|,-|MB|)          C
C*********************************************************************C
C
C     M QUANTUM NUMBERS
      MQNLAB(1) =-MQN(I1)
      MQNLAB(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQSS(ESS11,EXL,COORD,KQNLAB,MQNLAB,NFLAB,IQ)
C
C     MULTIPLY SOME OF THE ESS0 COEFFICIENTS BY A PHASE TERM
      ITUV=0
      DO MTUV=0,LAMAB
        PHASE = DFLOAT((IALT)**(MTUV))
        DO MDUM=1,((MTUV+1)*(MTUV+2))/2
          ITUV = ITUV+1
          DO M=1,NFUNS(I1)*NFUNS(I2)
            ESS11(M,ITUV) = ESS11(M,ITUV)*PHASE
          ENDDO
        ENDDO
      ENDDO
C
C*********************************************************************C
C     2: GENERATE AND STORE ESS21, FOR M PAIRS (+|MA|,-|MB|)          C
C*********************************************************************C
C
C     M QUANTUM NUMBERS
      MQNLAB(1) = MQN(I1)
      MQNLAB(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQSS(ESS21,EXL,COORD,KQNLAB,MQNLAB,NFLAB,IQ)
C
C     MULTIPLY SOME OF THE ES0 COEFFICIENTS BY A PHASE TERM
      ITUV=0
      DO MTUV=0,LAMAB
        PHASE = DFLOAT((IALT)**(MTUV))
        DO MDUM=1,((MTUV+1)*(MTUV+2))/2
         ITUV = ITUV+1
          DO M=1,NFUNS(I1)*NFUNS(I2)
            ESS21(M,ITUV) = ESS21(M,ITUV)*PHASE
          ENDDO
        ENDDO
      ENDDO
C
C     NOTE THAT ESS22 AND ESS12 ARE RELATED TO THESE BY SIMPLE FACTORS.
C
      RETURN
      END
C
C
      SUBROUTINE EMAKELS(ELS11,ELS21,EXPT,KQN,MQN,NFUNS,IALT,I1,I2,IQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C    EEEEEEEE MM     MM    AA    KK    KK EEEEEEEE LL      SSSSSS     C
C    EE       MMM   MMM   AAAA   KK   KK  EE       LL     SS    SS    C
C    EE       MMMM MMMM  AA  AA  KK  KK   EE       LL     SS          C
C    EEEEEE   MM MMM MM AA    AA KKKKK    EEEEEE   LL      SSSSSS     C
C    EE       MM  M  MM AAAAAAAA KK  KK   EE       LL           SS    C
C    EE       MM     MM AA    AA KK   KK  EE       LL     SS    SS    C
C    EEEEEEEE MM     MM AA    AA KK    KK EEEEEEEE LLLLLLL SSSSSS     C
C                                                                     C
C ------------------------------------------------------------------- C
C     EMAKELS GENERATES BLOCKS OF SPHERICAL SPINOR E-COEFFICIENTS     C
C     BY CONTRACTING ON THE E-COEFFICIENTS FOR VECTOR SPHERICAL       C
C     GAUSSIAN FUNCTIONS, USING A DEVELOPMENT OF THE ALGORITHM OF     C
C     V.R.SAUNDERS                                                    C
C ------------------------------------------------------------------- C
C     INPUT --                                                        C
C      I1: INDEX 1 (OF FOUR OPTIONS)                                  C
C      I2: INDEX 2 (OF FOUR OPTIONS)                                  C
C*********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=MKP,MBS=26,MB2=MBS*MBS,
     &          MLL=(MKP+8)*(MKP+9)*(MKP+10)/6)
C
      COMPLEX*16 ELS11(MB2,MLL),ELS21(MB2,MLL)
C
      DIMENSION EXPT(MBS,4),KQN(4),MQN(4),NFUNS(4)
      DIMENSION EXL(MBS,2),KQNLAB(2),MQNLAB(2),NFLAB(2),COORD(3,2)
C
      COMMON/XYZ4/XYZ(3,4)
C
C     DETERMINE LQNS FROM KQNS
      IF(KQN(I1).GT.0) THEN
        LA = KQN(I1)
      ELSE
        LA =-KQN(I1)-1
      ENDIF
      IF(KQN(I2).GT.0) THEN
        LB = KQN(I2)
      ELSE
        LB =-KQN(I2)-1
      ENDIF
C
C     THE SUMMATION TERMINATES AT LAMBDA = LA + LB + 1 FOR LS PAIRS
      LAMAB = LA + LB + 1
C
C     KAPPA LABELS AND NUMBER OF FUNCTIONS ON A GIVEN CENTRE
      KQNLAB(1) = KQN(I1)
      KQNLAB(2) = KQN(I2)
      NFLAB(1)  = NFUNS(I1)
      NFLAB(2)  = NFUNS(I2)
C
C     CARTESIAN COORDINATES OF EACH CENTRE
      DO I=1,3
        COORD(I,1) = XYZ(I,I1)
        COORD(I,2) = XYZ(I,I2)
      ENDDO
C
C     BASIS SET EXPONENTS FOR CENTRE A
      DO IBAS=1,NFUNS(I1)
       EXL(IBAS,1) = EXPT(IBAS,I1)
      ENDDO
C
C     BASIS SET EXPONENTS FOR CENTRE B
      DO JBAS=1,NFUNS(I2)
       EXL(JBAS,2) = EXPT(JBAS,I2)
      ENDDO    
C
C     EMPTY THE E-COEFFICIENT ARRAYS
      DO ITUV=1,MLL
        DO M=1,MB2
          ELS11(M,ITUV) = DCMPLX(0.0D0,0.0D0)
          ELS21(M,ITUV) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C*********************************************************************C
C     1: GENERATE AND STORE ELS11, FOR M PAIRS (-|MA|,-|MB|)          C
C*********************************************************************C
C
C     M QUANTUM NUMBERS
      MQNLAB(1) =-MQN(I1)
      MQNLAB(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLS(ELS11,EXL,COORD,KQNLAB,MQNLAB,NFLAB,IQ)
C
C     MULTIPLY SOME OF THE ELSQ COEFFICIENTS BY A PHASE TERM
      ITUV = 0
      DO MTUV=0,LAMAB
        PHASE = DFLOAT((IALT)**(MTUV))
        DO MDUM=1,((MTUV+1)*(MTUV+2))/2
         ITUV = ITUV+1
           DO M=1,NFUNS(I1)*NFUNS(I2)
             ELS11(M,ITUV) = ELS11(M,ITUV)*PHASE
           ENDDO
        ENDDO
      ENDDO
C
C*********************************************************************C
C     2: GENERATE AND STORE ELS21, FOR M PAIRS (+|MA|,-|MB|)          C
C*********************************************************************C
C
C     M QUANTUM NUMBERS
      MQNLAB(1) = MQN(I1)
      MQNLAB(2) =-MQN(I2)
C
C     GENERATE THE RAW COEFFICIENTS
      CALL EQLS(ELS21,EXL,COORD,KQNLAB,MQNLAB,NFLAB,IQ)
C
C     MULTIPLY SOME OF THE ELS1 COEFFICIENTS BY A PHASE TERM
      ITUV = 0
      DO MTUV=0,LAMAB
        PHASE = DFLOAT((IALT)**(MTUV))
        DO MDUM=1,((MTUV+1)*(MTUV+2))/2
          ITUV = ITUV + 1
          DO M=1,NFUNS(I1)*NFUNS(I2)
            ELS21(M,ITUV) = ELS21(M,ITUV)*PHASE
          ENDDO
        ENDDO
      ENDDO
C
C     NOTE THAT ELS22 AND ELS12 ARE RELATED TO THESE BY PHASE FACTORS.
C     THERE IS ALSO A RELATION WHICH ALLOWS US TO OBTAIN ESLQ FROM ELSQ.
C
      RETURN
      END

