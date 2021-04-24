      SUBROUTINE ONEEL0(HMAT,OVAP,EXL,ZCRG,KQN,NFUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          OOOOOO  NN    NN EEEEEEEE EEEEEEEE LL      000000           C
C         OO    OO NNN   NN EE       EE       LL     00    00          C
C         OO    OO NNNN  NN EE       EE       LL     00    00          C
C         OO    OO NN NN NN EEEEEE   EEEEEE   LL     00    00          C
C         OO    OO NN  NNNN EE       EE       LL     00    00          C
C         OO    OO NN   NNN EE       EE       LL     00    00          C
C          OOOOOO  NN    NN EEEEEEEE EEEEEEEE LLLLLLL 000000           C
C                                                                      C
C -------------------------------------------------------------------- C
C  ONEELRE0 CALCULATES THE DIRAC AND OVERLAP MATRICES FOR SYMMETRY     C
C  TYPE KQN, USING EVEN-TEMPERED SGTFS.                                C
C**********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=(MKP+1)/2,MBS=26)
C
      CHARACTER*4 HMLTN
C
      DIMENSION HMAT(2*MBS,2*MBS),OVAP(2*MBS,2*MBS),EXL(MBS)
      DIMENSION RN(MBS*MBS,4)
C
      COMMON/GMFN/GAMMAL(100),GAMMAF(100)
      COMMON/PRMS/CV,HMLTN,ITER,IIL,ITREE,IEQS
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR,IOCCM0
C
C     GENERATE GAMMA FUNCTIONS
      CALL GAMMAS
C
C     DETERMINE THE LQN
      IF(KQN.LT.0) THEN
        LQN =-(KQN+1)
      ELSE
        LQN =  KQN
      ENDIF
      RL  = DFLOAT(LQN)
      RL3 = DFLOAT(2*LQN+3)
C
C     GENERATE NORMALISATION FACTORS FOR THESE EXPONENTS
      CALL RNORM0(RN,EXL,NFUN,LQN)
C
C     CONSTRUCT THE DIRAC MATRIX
      IF(KQN.GT.0) THEN
C     DIRAC MATRIX FOR KQN > 0
        G = DFLOAT(2*LQN+1)
        M = 0
        DO IBAS=1,NFUN
          EI = EXL(IBAS)
          DO JBAS=1,NFUN
            M    = M+1
            EJ   = EXL(JBAS)
            EIJ  = EI + EJ
            EIJP = EI*EJ
            T3   = RINT(2*LQN+1,EIJ)
            T6   = RL + 0.5D0
            T7   = EIJ**T6
            T8   = 1.0D0/T7
            T10  = G**2
            T14  = RL + 1.5D0
            T16  = EIJ**2
            T20  = T10*0.5D0 - G*T6 + 2.0D0*EIJP*T6*T14/T16
            T21  = GAMMAF(2*LQN+1)*T8*T20*4.0D0
            T34  = CV*CV
            T40  = EIJ**T14
            T52  = EIJ**(RL + 2.5D0)
            F1   =-ZCRG*RN(M,1)*T3
            F2   = CV*RN(M,2)*T21
            F3   =-ZCRG*RN(M,3)*(T10*RINT(2*LQN-1,EIJ) 
     &              - 2.0D0*EIJ*G*T3 + 4.0D0*EIJP*RINT(2*LQN+3,EIJ)) 
     &              - 2.0D0*T34*RN(M,3)*T21
            F4   = RN(M,1)*2.0D0*GAMMAF(2*LQN+3)/T40
            F5   = RN(M,3)*GAMMAF(2*LQN+1)*4.0D0*T8*T20
C
C           OVERLAP MATRIX
            OVAP(IBAS     ,JBAS     ) = F4
            OVAP(IBAS+NFUN,JBAS     ) = 0.0D0
            OVAP(JBAS     ,IBAS+NFUN) = 0.0D0
            OVAP(IBAS+NFUN,JBAS+NFUN) = F5
C
C           DIRAC MATRIX
            HMAT(IBAS     ,JBAS     ) = F1
            HMAT(IBAS+NFUN,JBAS     ) = F2
            HMAT(JBAS     ,IBAS+NFUN) = HMAT(IBAS+NFUN,JBAS)
            HMAT(IBAS+NFUN,JBAS+NFUN) = F3
          ENDDO
        ENDDO
      ELSE
C     DIRAC MATRIX FOR KQN < 0
        M = 0
        DO IBAS=1,NFUN
          EI = EXL(IBAS)
          DO JBAS=1,NFUN
            M    = M + 1
            EJ   = EXL(JBAS)
            EIJ  = EI + EJ
            EIJP = EI*EJ
            T6   = RL + 0.5D0
            T7   = EIJ**T6
            T8   = 1.0D0/T7
            T12  = RL + 1.5D0
            T14  = EIJ**2
            T15  = 1.0D0/T14
            T16  = T6*T12*T15
            T24  = 2.0D0**2
            T25  = CV*CV
            T27  = RN(M,3)*GAMMAF(2*LQN+1)*4.0D0
            T34  = EIJ**T12
            F1   =-ZCRG*RN(M,1)*RINT(2*LQN+1,EIJ)
            F2   = CV*RN(M,2)*GAMMAF(2*LQN+1)*8.0D0*T8*EIJP*T16
            F3   =-4.0D0*ZCRG*RN(M,3)*EIJP*RINT(2*LQN+3,EIJ)
     &              - T24*T25*T27*T8*EIJP*T16
            F4   = RN(M,1)*2.0D0*GAMMAF(2*LQN+3)/T34
            F5   = T27*T8*2.0D0*EIJP*T6*T12*T15
            F6   = F4*RL3*EIJP/EIJ
C
C           OVERLAP MATRIX
            IF(HMLTN.EQ.'NORL') THEN
              OVAP(IBAS     ,JBAS     ) = F4
            ELSE
              OVAP(IBAS     ,JBAS     ) = F4
              OVAP(IBAS+NFUN,JBAS     ) = 0.0D0
              OVAP(JBAS     ,IBAS+NFUN) = 0.0D0
              OVAP(IBAS+NFUN,JBAS+NFUN) = F5          
            ENDIF
C
C           DIRAC MATRIX
            IF(HMLTN.EQ.'NORL') THEN
              HMAT(IBAS     ,JBAS     ) = F1 + F6
            ELSE
              HMAT(IBAS     ,JBAS     ) = F1
              HMAT(IBAS+NFUN,JBAS     ) = F2
              HMAT(JBAS     ,IBAS+NFUN) = HMAT(IBAS+NFUN,JBAS)
              HMAT(IBAS+NFUN,JBAS+NFUN) = F3     
            ENDIF
C         
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
