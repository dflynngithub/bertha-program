      SUBROUTINE ONEEL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C               OOOOOO  NN    NN EEEEEEE EEEEEEE LL                   C
C              OO    OO NNN   NN EE      EE      LL                   C
C              OO    OO NNNN  NN EE      EE      LL                   C
C              OO    OO NN NN NN EEEEEE  EEEEEE  LL                   C
C              OO    OO NN  NNNN EE      EE      LL                   C
C              OO    OO NN   NNN EE      EE      LL                   C
C               OOOOOO  NN    NN EEEEEEE EEEEEEE LLLLLLL              C
C                                                                     C
C ------------------------------------------------------------------- C
C     ONEEL CONSTRUCTS ONE-ELECTRON MULTI-CENTRE HMAT INTEGRALS       C
C     AND CALCULATES THE ONE-ELECTRON ENERGY OVER THE ATOMIC DENSITY. C
C*********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=MKP,MBS=26,MB2=MBS*MBS,
     &          IL4=2*(MKP-1),MLL=(MKP+8)*(MKP+9)*(MKP+10)/6,MRC=969)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 HMAT(MDM,MDM),GMAT(MDM,MDM),QMAT(MDM,MDM),
     &           BMAT(MDM,MDM),OVAP(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
      COMPLEX*16 HLL(MBS,MBS,4),HLS(MBS,MBS,4),SLL(MBS,MBS,4),
     &           HSL(MBS,MBS,4),HSS(MBS,MBS,4),SSS(MBS,MBS,4)
      COMPLEX*16 E11(MB2,MLL),E21(MB2,MLL)
      COMPLEX*16 E11A,E11B,E11C,TRM11,E21A,E21B,E21C,TRM21,ETMP,CTEMP(4)
C
      COMMON/ACSS/INABCD(0:IL4,0:IL4,0:IL4),
     &            IVEC(MLL),JVEC(MLL),KVEC(MLL),LAMVEC(MLL)
      COMMON/DENS/DENC,DENO,DENT
      COMMON/ENRG/ETOT,ENUC,EONE,ECLM,EBRT
      COMMON/MTRX/HMAT,GMAT,QMAT,BMAT,OVAP,FOCK
      COMMON/PRMS/CV,HMLTN,ITER,IALL,IRUN
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,IOCCM0
      COMMON/XYZ4/XYZ(3,4)
C
      DIMENSION RC(MB2,MRC),CP(MB2,3),APH(MB2),PNC(MB2),
     &          EXPT(MBS,4),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      DATA PI/3.1415926535897932D0/
C
C     INITIALISE ONE-ELECTRON ENERGY COUNTER
      EONE  = 0.0D0
C
C*********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)     C
C*********************************************************************C
C
C     LOOP OVER CENTRE A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 2000 ICNTB=1,ICNTA
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTRE/MULTI-CENTRE OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 2
        ENDIF
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
         LQN(1) = KQN(1)
        ELSE
         LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNA    = NFUNCT(LQN(1)+1,ICNTA)
        NFUNS(1) = NFUNA
C
        DO IBAS=1,NFUNA
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2)=KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNB    = NFUNCT(LQN(2)+1,ICNTB)
        NFUNS(2) = NFUNB
C
        DO IBAS=1,NFUNB
          EXPT(IBAS,2) = EXPSET(IBAS,LQN(2)+1,ICNTB)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 2000 MA=1,IABS(KQN(1))
        MJA    = 2*MA-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 2000 MB=1,IABS(KQN(2))
        MJB    = 2*MB-1
        MQN(2) = MJB
C
C*********************************************************************C
C     IMPLEMENT ANY AVAILABLE SELECTION RULES HERE                    C
C*********************************************************************C
C
C     SELECTION RULES TO BE MADE BASED ON GEOMETRIC SYMMETRY,
C     ATOMIC COORDINATES AND QUANTUM NUMBER PAIRS. THE IDEA IS TO
C     SKIP CALCULATIONS THAT SHOULD PRODUCE NO OR NEGLIGIBLE EFFECT.
      IF(NCNT.LE.2) THEN
        IF(MQN(1).NE.MQN(2)) GOTO 2000
      ENDIF
C
C*********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS          C
C     OF (MA,MB). FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE         C
C     ORDERED                                                         C
C     11: = (-|MQN(A)|,-|MQN(B)|) -> 1                                C
C     12: = (-|MQN(A)|,+|MQN(B)|) -> 2                                C
C     21: = (+|MQN(A)|,-|MQN(B)|) -> 3                                C
C     22: = (+|MQN(A)|,+|MQN(B)|) -> 4                                C
C*********************************************************************C
C
C     CONSTRUCTION OF ONE-ELECTRON HMAT MATRIX BY TT' BLOCKS...
C
C     THE PHASE GENERATES E022 AND E012 COEFFS FROM E011 AND E021
      FASE = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
     &                    *DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXM = NFUNA*NFUNB
C
C*********************************************************************C
C     PART 1: THE LL MATRICES                                         C
C*********************************************************************C
C
C     CALCULATE LAMBDA VALUES FOR THIS OVERLAP CHOICE
      LAMBDA = LQN(1)+LQN(2)
      NTUVLL = ((LAMBDA+1)*(LAMBDA+2)*(LAMBDA+3))/6
C
C     GENERATE ELL0 COEFFICIENTS
      IALT = 1
      CALL EMAKELL(E11,E21,EXPT,KQN,MQN,NFUNS,IALT,1,2,0)
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO I=1,NFUNA
        DO J=1,NFUNB
          M = M+1
          EIJ    = EXPT(I,1)+EXPT(J,2)
          EROOT  = DSQRT(PI/EIJ)**3
          SLL(I,J,1) = EROOT*E11(M,1)
          SLL(I,J,3) = EROOT*E21(M,1)
          SLL(I,J,2) =-FASE*DCONJG(SLL(I,J,3))
          SLL(I,J,4) = FASE*DCONJG(SLL(I,J,1))
        ENDDO
      ENDDO
C
C     INITIALISE KINETIC ARRAYS
      DO I=1,NFUNA
        DO J=1,NFUNB
          DO IB=1,4
            HLL(I,J,IB) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C     CONSTRUCT NON-REL KINETIC ENERGY INTEGRALS (IOS 91)
      IF(HMLTN.NE.'NORL') GOTO 400
      RL2 = DFLOAT(2*LQN(2)+3)
      M   = 0
      DO I=1,NFUNA
        DO J=1,NFUNB
          M     = M+1
          EJ    = EXPT(J,2)
          EIJ   = EXPT(I,1) + EXPT(J,2)
          EROOT = DSQRT(PI/EIJ)**3
          PX = (XYZ(1,1)*EXPT(I,1) + XYZ(1,2)*EXPT(J,2))/EIJ
          PY = (XYZ(2,1)*EXPT(I,1) + XYZ(2,2)*EXPT(J,2))/EIJ
          PZ = (XYZ(3,1)*EXPT(I,1) + XYZ(3,2)*EXPT(J,2))/EIJ
C
          PBX = PX - XYZ(1,2)
          PBY = PY - XYZ(2,2)
          PBZ = PZ - XYZ(3,2)
          PB2 = (PBX*PBX) + (PBY*PBY) + (PBZ*PBZ)
C       
          E0FC = EJ*RL2 - 2.0D0*EJ*EJ*(PB2 + 1.50D0/EIJ)
          E1FC = 4.0D0*EJ*EJ
C
C         TRUNCATE EXPRESSION DEPENDING ON LAMBDA VALUE
C         ALL COMBINATIONS ALLOW FOR THE LAMBDA = 0 MANIFOLD           
          TRM11 = E0FC*E11(M,1)
          TRM21 = E0FC*E21(M,1)
C         IF LAMBDA > 0 PROVIDE SECOND BUNCH OF TERMS
          IF(LAMBDA.GE.1) THEN
            E11A = E11(M,INABCD(1,0,0))
            E21A = E21(M,INABCD(1,0,0))
            E11B = E11(M,INABCD(0,1,0))
            E21B = E21(M,INABCD(0,1,0))
            E11C = E11(M,INABCD(0,0,1))
            E21C = E21(M,INABCD(0,0,1))
            TRM11 = TRM11 - E1FC*(PBX*E11A + PBY*E11B + PBZ*E11C)
            TRM21 = TRM21 - E1FC*(PBX*E21A + PBY*E21B + PBZ*E21C)
          ENDIF
C         IF LAMBDA > 1 PROVIDE FINAL BUNCH OF TERMS
          IF(LAMBDA.GE.2) THEN
            E11A = E11(M,INABCD(2,0,0))
            E21A = E21(M,INABCD(2,0,0))
            E11B = E11(M,INABCD(0,2,0))
            E21B = E21(M,INABCD(0,2,0))
            E11C = E11(M,INABCD(0,0,2))
            E21C = E21(M,INABCD(0,0,2))
            TRM11 = TRM11 - E1FC*(E11A + E11B + E11C)
            TRM21 = TRM21 - E1FC*(E21A + E21B + E21C)
          ENDIF
          HLL(I,J,1) = EROOT*TRM11
          HLL(I,J,3) = EROOT*TRM21
        ENDDO
      ENDDO
400   CONTINUE
C
C     LOOP OVER ALL NUCLEAR CENTRES FOR NUCLEAR ATTRACTION INTEGRALS
      DO IZ=1,NCNT
C
C       NUCLEAR COORDINATES
        CX = COORD(1,IZ)
        CY = COORD(2,IZ)
        CZ = COORD(3,IZ)
C
C       GAUSSIAN PRODUCT THEOREM OVER BASIS FUNCTIONS
        M = 0
        DO I=1,NFUNA
          DO J=1,NFUNB
            M   = M+1
            EIJ = EXPT(I,1)+EXPT(J,2)
            ESM = CNUC(IZ)+EIJ
            PX  = (XYZ(1,1)*EXPT(I,1) + XYZ(1,2)*EXPT(J,2))/EIJ
            PY  = (XYZ(2,1)*EXPT(I,1) + XYZ(2,2)*EXPT(J,2))/EIJ
            PZ  = (XYZ(3,1)*EXPT(I,1) + XYZ(3,2)*EXPT(J,2))/EIJ
            APH(M) = (EIJ*CNUC(IZ))/ESM
            PNC(M) = 2.0D0*PI*DSQRT(CNUC(IZ)/ESM)*ZNUC(IZ)/EIJ
            CP(M,1) = CX - PX
            CP(M,2) = CY - PY
            CP(M,3) = CZ - PZ
          ENDDO
        ENDDO
C
C       GENERATE A BATCH OF R-INTEGRALS
        CALL RMAKE(RC,CP,APH,MAXM,LAMBDA)
C
C       NUCLEAR ATTRACTION INTEGRALS AS A FINITE SUME OF ELL0 AND RC
        M = 0
        DO I=1,NFUNA
          DO J=1,NFUNB
            M = M+1
            DO ITUV=1,NTUVLL
              HLL(I,J,1) = HLL(I,J,1) - PNC(M)*E11(M,ITUV)*RC(M,ITUV)
              HLL(I,J,3) = HLL(I,J,3) - PNC(M)*E21(M,ITUV)*RC(M,ITUV)
            ENDDO
            HLL(I,J,2) =-FASE*DCONJG(HLL(I,J,3))
            HLL(I,J,4) = FASE*DCONJG(HLL(I,J,1))
          ENDDO
        ENDDO
      ENDDO
C
C     NON-REL HAMILTONIAN MATRICES COMPLETE
      IF(HMLTN.EQ.'NORL') GOTO 500
C     
C*********************************************************************C
C     PART 2: THE SS MATRICES                                         C
C*********************************************************************C
C
C     CALCULATE LAMBDA VALUES FOR THIS OVERLAP CHOICE
      LAMBDA = LQN(1)+LQN(2)+2
      NTUVSS = ((LAMBDA+1)*(LAMBDA+2)*(LAMBDA+3))/6
C
C     GENERATE ESS0 COEFFICIENTS
      IALT = 1
      CALL EMAKESS(E11,E21,EXPT,KQN,MQN,NFUNS,IALT,1,2,0)
C
C     OVERLAP MATRIX ELEMENTS
      M=0
      DO I=1,NFUNA
        DO J=1,NFUNB
          M = M+1
          EIJ     = EXPT(I,1)+EXPT(J,2)
          EROOT   = DSQRT(PI/EIJ)**3
          SSS(I,J,1) = EROOT*E11(M,1)
          SSS(I,J,3) = EROOT*E21(M,1)
          SSS(I,J,2) =-FASE*DCONJG(SSS(I,J,3))
          SSS(I,J,4) = FASE*DCONJG(SSS(I,J,1))
        ENDDO
      ENDDO
C
C     INITIALISE KINETIC ARRAYS
      DO I=1,NFUNA
        DO J=1,NFUNB
          DO IB=1,4
            HSS(I,J,IB) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C     LOOP OVER ALL NUCLEAR CENTRES FOR NUCLEAR ATTRACTION INTEGRALS
      DO IZ=1,NCNT
C
C       NUCLEAR COORDINATES
        CX = COORD(1,IZ)
        CY = COORD(2,IZ)
        CZ = COORD(3,IZ)
C
C       GAUSSIAN PRODUCT THEOREM OVER BASIS FUNCTIONS
        M = 0
        DO I=1,NFUNA
          DO J=1,NFUNB
            M   = M+1
            EIJ = EXPT(I,1)+EXPT(J,2)
            ESM = CNUC(IZ) + EIJ
            PX  = (XYZ(1,1)*EXPT(I,1) + XYZ(1,2)*EXPT(J,2))/EIJ
            PY  = (XYZ(2,1)*EXPT(I,1) + XYZ(2,2)*EXPT(J,2))/EIJ
            PZ  = (XYZ(3,1)*EXPT(I,1) + XYZ(3,2)*EXPT(J,2))/EIJ
            APH(M) = (EIJ*CNUC(IZ))/ESM
            PNC(M) = 2.0D0*PI*DSQRT(CNUC(IZ)/ESM)*ZNUC(IZ)/EIJ
            CP(M,1) = CX - PX
            CP(M,2) = CY - PY
            CP(M,3) = CZ - PZ
          ENDDO
        ENDDO      
C
C       GENERATE A BATCH OF R-INTEGRALS
        CALL RMAKE(RC,CP,APH,MAXM,LAMBDA)
C
C       NUCLEAR ATTRACTION INTEGRALS AS A FINITE SUME OF ESS0 AND RC
        M = 0
        DO I=1,NFUNA
          DO J=1,NFUNB
            M = M+1
            DO ITUV=1,NTUVSS
              HSS(I,J,1) = HSS(I,J,1) + PNC(M)*E11(M,ITUV)*RC(M,ITUV)
              HSS(I,J,3) = HSS(I,J,3) + PNC(M)*E21(M,ITUV)*RC(M,ITUV)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     SUBTRACT THE SS OVERLAP MATRIX AND FINISH CONSTRUCTION
      DO I=1,NFUNA
        DO J=1,NFUNB
          HSS(I,J,1) =-HSS(I,J,1) - 2.0D0*CV*CV*SSS(I,J,1)
          HSS(I,J,3) =-HSS(I,J,3) - 2.0D0*CV*CV*SSS(I,J,3)
          HSS(I,J,2) =-FASE*DCONJG(HSS(I,J,3))
          HSS(I,J,4) = FASE*DCONJG(HSS(I,J,1))
        ENDDO
      ENDDO
C
C*********************************************************************C
C     PART 3: THE SL MATRICES                                         C
C*********************************************************************C
C
C     CALCULATE LAMBDA VALUES FOR THIS OVERLAP CHOICE
      LAMBDA = LQN(1)+LQN(2)+2
      NTUVSS = ((LAMBDA+1)*(LAMBDA+2)*(LAMBDA+3))/6
C
C     INITIALISE KINETIC ARRAYS
      DO I=1,NFUNA
        DO J=1,NFUNB
          DO IB=1,4
            HSL(I,J,IB) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C     KINETIC MATRIX ELEMENTS
      FACT = CV*DSQRT(DFLOAT(2*LQN(2)+3))
      M = 0
      DO I=1,NFUNA
        DO J=1,NFUNB
          M = M+1
          EIJ    = EXPT(I,1) + EXPT(J,2)
          EJRT   = FACT*DSQRT(EXPT(J,2))
          EROOT  = DSQRT(PI/EIJ)**3
          HSL(I,J,1) = EJRT*EROOT*E11(M,1)
          HSL(I,J,3) = EJRT*EROOT*E21(M,1)
          HSL(I,J,2) =-FASE*DCONJG(HSL(I,J,3))
          HSL(I,J,4) = FASE*DCONJG(HSL(I,J,1))
        ENDDO
      ENDDO
C
C
C*********************************************************************C
C     PART 4: THE LS MATRICES                                         C
C*********************************************************************C
C
C     CALCULATE LAMBDA VALUES FOR THIS OVERLAP CHOICE
      LAMBDA = LQN(1)+LQN(2)+2
      NTUVSS = ((LAMBDA+1)*(LAMBDA+2)*(LAMBDA+3))/6
C
C     INITIALISE KINETIC ARRAYS
      DO I=1,NFUNA
        DO J=1,NFUNB
          DO IB=1,4
            HLS(I,J,IB) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C     GENERATE ESS0 COEFFICIENTS
      IALT = 1
      CALL RNORMF(EXPT,LQN,NFUNS,2,1)
      CALL EMAKESS(E11,E21,EXPT,KQN,MQN,NFUNS,IALT,2,1,0)
C
C     KINETIC MATRIX ELEMENTS
      FACT = CV*DSQRT(DFLOAT(2*LQN(1)+3))
      M = 0
      DO J=1,NFUNB
        DO I=1,NFUNA
          M = M+1
          EIJ    = EXPT(J,2) + EXPT(I,1)
          EIRT   = FACT*DSQRT(EXPT(I,1))
          EROOT  = DSQRT(PI/EIJ)**3
          HLS(I,J,1) = EIRT*EROOT*E11(M,1)
          HLS(I,J,3) = EIRT*EROOT*E21(M,1)
          HLS(I,J,2) =-FASE*DCONJG(HLS(I,J,3))
          HLS(I,J,4) = FASE*DCONJG(HLS(I,J,1))
        ENDDO
      ENDDO
C
C     GENERATE LS MATRICES FROM THE ABOVE SL MATRICES
      M = 0
      DO J=1,NFUNB
        DO I=1,NFUNA
          M = M + 1
          CTEMP(1) = DCONJG(HLS(I,J,1))
          CTEMP(2) = DCONJG(HLS(I,J,2))
          CTEMP(3) = DCONJG(HLS(I,J,3))
          CTEMP(4) = DCONJG(HLS(I,J,4))
C
          HLS(I,J,1) = CTEMP(1)
          HLS(I,J,2) = CTEMP(3)
          HLS(I,J,3) = CTEMP(2)
          HLS(I,J,4) = CTEMP(4)
        ENDDO
      ENDDO
C
500   CONTINUE
C
C*********************************************************************C
C     WE NOW HAVE ALL PIECES OF THE HMAT MATRIX FOR THIS BLOCK OF     C
C     BASIS FUNCTIONS -- NOW OVERLAY THE RESULTS INTO HMAT AND OVAP.  C
C*********************************************************************C
C
C     CALCULATE COMPONENT OFFSETS
      IL1 = LARGE(ICNTA,KA,MJA  )
      IL2 = LARGE(ICNTA,KA,MJA+1)
      JL1 = LARGE(ICNTB,KB,MJB  )
      JL2 = LARGE(ICNTB,KB,MJB+1)
C
      IS1 = IL1 + NSHIFT
      IS2 = IL2 + NSHIFT
      JS1 = JL1 + NSHIFT
      JS2 = JL2 + NSHIFT
C
C     LL BLOCKS
      IF(IL1.GT.JL1) THEN
        DO J=1,NFUNB
          DO I=1,NFUNA
            HMAT(IL1+I,JL1+J) = HMAT(IL1+I,JL1+J) + HLL(I,J,1)
            HMAT(IL1+I,JL2+J) = HMAT(IL1+I,JL2+J) + HLL(I,J,2)
            HMAT(IL2+I,JL1+J) = HMAT(IL2+I,JL1+J) + HLL(I,J,3)
            HMAT(IL2+I,JL2+J) = HMAT(IL2+I,JL2+J) + HLL(I,J,4)
C
            HMAT(JL1+J,IL1+I) = DCONJG(HMAT(IL1+I,JL1+J))
            HMAT(JL2+J,IL1+I) = DCONJG(HMAT(IL1+I,JL2+J))
            HMAT(JL1+J,IL2+I) = DCONJG(HMAT(IL2+I,JL1+J))
            HMAT(JL2+J,IL2+I) = DCONJG(HMAT(IL2+I,JL2+J))
C
            OVAP(IL1+I,JL1+J) = SLL(I,J,1)
            OVAP(IL1+I,JL2+J) = SLL(I,J,2)
            OVAP(IL2+I,JL1+J) = SLL(I,J,3)
            OVAP(IL2+I,JL2+J) = SLL(I,J,4)
C
            OVAP(JL1+J,IL1+I) = DCONJG(OVAP(IL1+I,JL1+J))
            OVAP(JL2+J,IL1+I) = DCONJG(OVAP(IL1+I,JL2+J))
            OVAP(JL1+J,IL2+I) = DCONJG(OVAP(IL2+I,JL1+J))
            OVAP(JL2+J,IL2+I) = DCONJG(OVAP(IL2+I,JL2+J))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO J=1,NFUNB
          DO I=J,NFUNA
            HMAT(IL1+I,JL1+J) = HMAT(IL1+I,JL1+J) + HLL(I,J,1)
            HMAT(IL1+I,JL2+J) = HMAT(IL1+I,JL2+J) + HLL(I,J,2)
            HMAT(IL2+I,JL1+J) = HMAT(IL2+I,JL1+J) + HLL(I,J,3)
            HMAT(IL2+I,JL2+J) = HMAT(IL2+I,JL2+J) + HLL(I,J,4)
C
            HMAT(JL1+J,IL1+I) = DCONJG(HMAT(IL1+I,JL1+J))
            HMAT(JL2+J,IL1+I) = DCONJG(HMAT(IL1+I,JL2+J))
            HMAT(JL1+J,IL2+I) = DCONJG(HMAT(IL2+I,JL1+J))
            HMAT(JL2+J,IL2+I) = DCONJG(HMAT(IL2+I,JL2+J))
C
            OVAP(IL1+I,JL1+J) = SLL(I,J,1)
            OVAP(IL1+I,JL2+J) = SLL(I,J,2)
            OVAP(IL2+I,JL1+J) = SLL(I,J,3)
            OVAP(IL2+I,JL2+J) = SLL(I,J,4)
C
            OVAP(JL1+J,IL1+I) = DCONJG(OVAP(IL1+I,JL1+J))
            OVAP(JL2+J,IL1+I) = DCONJG(OVAP(IL1+I,JL2+J))
            OVAP(JL1+J,IL2+I) = DCONJG(OVAP(IL2+I,JL1+J))
            OVAP(JL2+J,IL2+I) = DCONJG(OVAP(IL2+I,JL2+J))
          ENDDO
        ENDDO
      ENDIF
C
C     NON-REL HAMILTONIAN MATRIX COMPLETE
      IF(HMLTN.EQ.'NORL') GOTO 600
C
C     LS BLOCKS
      IF(IL1.GE.JL1) THEN
        DO J=1,NFUNB
          DO I=1,NFUNA
            HMAT(IL1+I,JS1+J) = HMAT(IL1+I,JS1+J) + HLS(I,J,1)
            HMAT(IL1+I,JS2+J) = HMAT(IL1+I,JS2+J) + HLS(I,J,2)
            HMAT(IL2+I,JS1+J) = HMAT(IL2+I,JS1+J) + HLS(I,J,3)
            HMAT(IL2+I,JS2+J) = HMAT(IL2+I,JS2+J) + HLS(I,J,4)
C
            HMAT(JS1+J,IL1+I) = DCONJG(HMAT(IL1+I,JS1+J))
            HMAT(JS2+J,IL1+I) = DCONJG(HMAT(IL1+I,JS2+J))
            HMAT(JS1+J,IL2+I) = DCONJG(HMAT(IL2+I,JS1+J))
            HMAT(JS2+J,IL2+I) = DCONJG(HMAT(IL2+I,JS2+J))
          ENDDO
        ENDDO
      ENDIF
C
C     SL BLOCKS
      IF(IL1.GT.JL1) THEN
        DO J=1,NFUNB
          DO I=1,NFUNA
            HMAT(IS1+I,JL1+J) = HMAT(IS1+I,JL1+J) + HSL(I,J,1)
            HMAT(IS1+I,JL2+J) = HMAT(IS1+I,JL2+J) + HSL(I,J,2)
            HMAT(IS2+I,JL1+J) = HMAT(IS2+I,JL1+J) + HSL(I,J,3)
            HMAT(IS2+I,JL2+J) = HMAT(IS2+I,JL2+J) + HSL(I,J,4)
C
            HMAT(JL1+J,IS1+I) = DCONJG(HMAT(IS1+I,JL1+J))
            HMAT(JL2+J,IS1+I) = DCONJG(HMAT(IS1+I,JL2+J))
            HMAT(JL1+J,IS2+I) = DCONJG(HMAT(IS2+I,JL1+J))
            HMAT(JL2+J,IS2+I) = DCONJG(HMAT(IS2+I,JL2+J))
          ENDDO
        ENDDO
      ENDIF
C
C     SS BLOCKS
      IF(IL1.GT.JL1) THEN
        DO J=1,NFUNB
          DO I=1,NFUNA
            HMAT(IS1+I,JS1+J) = HMAT(IS1+I,JS1+J) + HSS(I,J,1)
            HMAT(IS1+I,JS2+J) = HMAT(IS1+I,JS2+J) + HSS(I,J,2)
            HMAT(IS2+I,JS1+J) = HMAT(IS2+I,JS1+J) + HSS(I,J,3)
            HMAT(IS2+I,JS2+J) = HMAT(IS2+I,JS2+J) + HSS(I,J,4)
C
            HMAT(JS1+J,IS1+I) = DCONJG(HMAT(IS1+I,JS1+J))
            HMAT(JS2+J,IS1+I) = DCONJG(HMAT(IS1+I,JS2+J))
            HMAT(JS1+J,IS2+I) = DCONJG(HMAT(IS2+I,JS1+J))
            HMAT(JS2+J,IS2+I) = DCONJG(HMAT(IS2+I,JS2+J))
C
            OVAP(IS1+I,JS1+J) = SSS(I,J,1)
            OVAP(IS1+I,JS2+J) = SSS(I,J,2)
            OVAP(IS2+I,JS1+J) = SSS(I,J,3)
            OVAP(IS2+I,JS2+J) = SSS(I,J,4)
C
            OVAP(JS1+J,IS1+I) = DCONJG(OVAP(IS1+I,JS1+J))
            OVAP(JS2+J,IS1+I) = DCONJG(OVAP(IS1+I,JS2+J))
            OVAP(JS1+J,IS2+I) = DCONJG(OVAP(IS2+I,JS1+J))
            OVAP(JS2+J,IS2+I) = DCONJG(OVAP(IS2+I,JS2+J))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO J=1,NFUNB
          DO I=J,NFUNA
            HMAT(IS1+I,JS1+J) = HMAT(IS1+I,JS1+J) + HSS(I,J,1)
            HMAT(IS1+I,JS2+J) = HMAT(IS1+I,JS2+J) + HSS(I,J,2)
            HMAT(IS2+I,JS1+J) = HMAT(IS2+I,JS1+J) + HSS(I,J,3)
            HMAT(IS2+I,JS2+J) = HMAT(IS2+I,JS2+J) + HSS(I,J,4)
C
            HMAT(JS1+J,IS1+I) = DCONJG(HMAT(IS1+I,JS1+J))
            HMAT(JS2+J,IS1+I) = DCONJG(HMAT(IS1+I,JS2+J))
            HMAT(JS1+J,IS2+I) = DCONJG(HMAT(IS2+I,JS1+J))
            HMAT(JS2+J,IS2+I) = DCONJG(HMAT(IS2+I,JS2+J))
C
            OVAP(IS1+I,JS1+J) = SSS(I,J,1)
            OVAP(IS1+I,JS2+J) = SSS(I,J,2)
            OVAP(IS2+I,JS1+J) = SSS(I,J,3)
            OVAP(IS2+I,JS2+J) = SSS(I,J,4)
C
            OVAP(JS1+J,IS1+I) = DCONJG(OVAP(IS1+I,JS1+J))
            OVAP(JS2+J,IS1+I) = DCONJG(OVAP(IS1+I,JS2+J))
            OVAP(JS1+J,IS2+I) = DCONJG(OVAP(IS2+I,JS1+J))
            OVAP(JS2+J,IS2+I) = DCONJG(OVAP(IS2+I,JS2+J))
          ENDDO
        ENDDO
      ENDIF 
C
600   CONTINUE      
C
C*********************************************************************C
C     CALCULATE THE ONE-ELECTRON CONTRIBUTION TO THE TOTAL ENERGY     C
C*********************************************************************C
C
C     TEMPORARY ENERGY CONTER FOR THIS BLOCK
      ETMP = DCMPLX(0.0D0,0.0D0)
C
C     ALL CONTRIBUTIONS TO ONE-ELECTRON ENERGY IN THIS BLOCK
      M = 0
      DO I=1,NFUNA
        DO J=1,NFUNB
          M = M+1
          ETMP = ETMP
     &    + DENC(IL1+I,JL1+J)*HLL(I,J,1) + DENC(IL1+I,JL2+J)*HLL(I,J,2) 
     &    + DENC(IL2+I,JL1+J)*HLL(I,J,3) + DENC(IL2+I,JL2+J)*HLL(I,J,4)    
          IF(HMLTN.EQ.'NORL') GOTO 700
          ETMP = ETMP
     &    + DENC(IL1+I,JS1+J)*HLS(I,J,1) + DENC(IL1+I,JS2+J)*HLS(I,J,2)
     &    + DENC(IL2+I,JS1+J)*HLS(I,J,3) + DENC(IL2+I,JS2+J)*HLS(I,J,4)
     &    + DENC(IS1+I,JL1+J)*HSL(I,J,1) + DENC(IS1+I,JL2+J)*HSL(I,J,2) 
     &    + DENC(IS2+I,JL1+J)*HSL(I,J,3) + DENC(IS2+I,JL2+J)*HSL(I,J,4) 
     &    + DENC(IS1+I,JS1+J)*HSS(I,J,1) + DENC(IS1+I,JS2+J)*HSS(I,J,2)
     &    + DENC(IS2+I,JS1+J)*HSS(I,J,3) + DENC(IS2+I,JS2+J)*HSS(I,J,4)
700       CONTINUE
        ENDDO
      ENDDO
C
C     DECIDE WHETHER THIS BLOCK IS ATOMIC OR MULTI-CENTRE
      IF(ICNTA.EQ.ICNTB) THEN
        SYMFAC = 1.0D0
      ELSE
        SYMFAC = 2.0D0
      ENDIF
C
C     UPDATE ONE-ELECTRON ENERGY
      EONE = EONE + SYMFAC*DREAL(ETMP)
C
2000  CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE OVRLAP   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C          OOOOOO  VV    VV RRRRRRR  LL         AA    PPPPPPP         C
C         OO    OO VV    VV RR    RR LL        AAAA   PP    PP        C
C         OO    OO VV    VV RR    RR LL       AA  AA  PP    PP        C
C         OO    OO VV    VV RR    RR LL      AA    AA PP    PP        C
C         OO    OO  VV  VV  RRRRRRR  LL      AAAAAAAA PPPPPPP         C
C         OO    OO   VVVV   RR    RR LL      AA    AA PP              C
C          OOOOOO     VV    RR    RR LLLLLLL AA    AA PP              C
C                                                                     C
C ------------------------------------------------------------------- C
C     OVRLAP CONSTRUCTS THE ONE-ELECTRON OVERLAP MATRIX.              C
C*********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=MKP,MBS=26,MB2=MBS*MBS,
     &          IL4=2*(MKP-1),MLL=(MKP+8)*(MKP+9)*(MKP+10)/6,MRC=969)
C
      CHARACTER*4 HMLTN
C   
      COMPLEX*16 CONE
      COMPLEX*16 HMAT(MDM,MDM),GMAT(MDM,MDM),QMAT(MDM,MDM),
     &           BMAT(MDM,MDM),OVAP(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 E11(MB2,MLL),E21(MB2,MLL)
      COMPLEX*16 SLL(MBS,MBS,4),SSS(MBS,MBS,4)
C
      COMMON/ACSS/INABCD(0:IL4,0:IL4,0:IL4),
     &            IVEC(MLL),JVEC(MLL),KVEC(MLL),LAMVEC(MLL)
      COMMON/PRMS/CV,HMLTN,ITER,IALL,IRUN
      COMMON/MTRX/HMAT,GMAT,QMAT,BMAT,OVAP,FOCK
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,IOCCM0
      COMMON/XYZ4/XYZ(3,4)
C
      DIMENSION RC(MB2,MRC),EXPT(MBS,4),APH(MB2),CP(MB2,3),
     &          PNC(MB2),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      DATA PI/3.1415926535897932D0/
C
      CONE = DCMPLX(0.0D0,1.0D0)
C
C     INITIALISE THE OVERLAP ARRAY
      DO I=1,NDIM
        DO J=1,NDIM
          OVAP(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C*********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)     C
C*********************************************************************C
C
C     LOOP OVER CENTRE A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = COORD(1,ICNTA)
        XYZ(2,1) = COORD(2,ICNTA)
        XYZ(3,1) = COORD(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = COORD(1,ICNTB)
        XYZ(2,2) = COORD(2,ICNTB)
        XYZ(3,2) = COORD(3,ICNTB)
C
C       PARAMETER FOR SINGLE-CENTRE/MULTI-CENTRE OVERLAP OVER A AND B
        IF(ICNTA.EQ.ICNTB) THEN
          INUCAB = 1
        ELSE
          INUCAB = 2
        ENDIF
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR A
        KQN(1) = KVALS(KA,ICNTA)
        IF(KQN(1).GT.0) THEN
         LQN(1) = KQN(1)
        ELSE
         LQN(1) =-KQN(1)-1
        ENDIF
C
        NFUNA    = NFUNCT(LQN(1)+1,ICNTA)
        NFUNS(1) = NFUNA
C
        DO IBAS=1,NFUNA
          EXPT(IBAS,1) = EXPSET(IBAS,LQN(1)+1,ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS AND BASIS EXPONENTS FOR B
        KQN(2)=KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNB    = NFUNCT(LQN(2)+1,ICNTB)
        NFUNS(2) = NFUNB
C
        DO IBAS=1,NFUNB
          EXPT(IBAS,2) = EXPSET(IBAS,LQN(2)+1,ICNTB)
        ENDDO
C
C     LOOP OVER |MQN(A)| VALUES
      DO 2000 MA=1,IABS(KQN(1))
        MJA    = (2*MA)-1
        MQN(1) = MJA
C
C     LOOP OVER |MQN(B)| VALUES
      DO 2000 MB=1,IABS(KQN(2))
        MJB    = (2*MB)-1
        MQN(2) = MJB
C
C*********************************************************************C
C     IMPLEMENT ANY AVAILABLE SELECTION RULES HERE                    C
C*********************************************************************C
C
C     SELECTION RULES TO BE MADE BASED ON GEOMETRIC SYMMETRY,
C     ATOMIC COORDINATES AND QUANTUM NUMBER PAIRS. THE IDEA IS TO
C     SKIP CALCULATIONS THAT SHOULD PRODUCE NO OR NEGLIGIBLE EFFECT.
C     IF(MJA.NE.MJB) GOTO 2000
C
C*********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS          C
C     OF (MA,MB). FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE         C
C     ORDERED                                                         C
C     11: = (-|MQN(A)|,-|MQN(B)|)                                     C
C     12: = (-|MQN(A)|,+|MQN(B)|)                                     C
C     21: = (+|MQN(A)|,-|MQN(B)|)                                     C
C     22: = (+|MQN(A)|,+|MQN(B)|)                                     C
C*********************************************************************C
C
C     CONSTRUCTION OF ONE-ELECTRON OVERLAP MATRIX BY TT' BLOCKS...
C
C     THE PHASE GENERATES E022 AND E012 COEFFS FROM E011 AND E021
      FASE = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
     &                    *DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXM = NFUNA*NFUNB
C
C*********************************************************************C
C     PART 1: THE LL MATRICES                                         C
C*********************************************************************C
C
C     CALCULATE LAMBDA VALUES FOR THIS OVERLAP CHOICE
      LAMBDA = LQN(1)+LQN(2)
      NTUVLL = ((LAMBDA+1)*(LAMBDA+2)*(LAMBDA+3))/6
C
C     GENERATE ELL0 COEFFICIENTS
      IALT = 1
      CALL EMAKELL(E11,E21,EXPT,KQN,MQN,NFUNS,IALT,1,2,0)
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO I=1,NFUNA
        DO J=1,NFUNB
          M = M+1
          EIJ    = EXPT(I,1)+EXPT(J,2)
          EROOT  = DSQRT(PI/EIJ)**3
          SLL(I,J,1) = EROOT*E11(M,1)
          SLL(I,J,3) = EROOT*E21(M,1)
          SLL(I,J,2) =-FASE*DCONJG(SLL(I,J,3))
          SLL(I,J,4) = FASE*DCONJG(SLL(I,J,1))
        ENDDO
      ENDDO
C
C*********************************************************************C
C     PART 2: THE SS MATRICES                                         C
C*********************************************************************C
C
C     CALCULATE LAMBDA VALUES FOR THIS OVERLAP CHOICE
      LAMBDA = LQN(1)+LQN(2)+2
      NTUVSS = ((LAMBDA+1)*(LAMBDA+2)*(LAMBDA+3))/6
C
C     GENERATE ESS0 COEFFICIENTS
      IALT = 1
      CALL EMAKESS(E11,E21,EXPT,KQN,MQN,NFUNS,IALT,1,2,0)
C
C     OVERLAP MATRIX ELEMENTS
      M=0
      DO I=1,NFUNA
        DO J=1,NFUNB
          M = M+1
          EIJ     = EXPT(I,1)+EXPT(J,2)
          EROOT   = DSQRT(PI/EIJ)**3
          SSS(I,J,1) = EROOT*E11(M,1)
          SSS(I,J,3) = EROOT*E21(M,1)
          SSS(I,J,2) =-FASE*DCONJG(SSS(I,J,3))
          SSS(I,J,4) = FASE*DCONJG(SSS(I,J,1))
        ENDDO
      ENDDO
C
C*********************************************************************C
C     WE NOW HAVE ALL PIECES OF THE OVERLAP MATRIX FOR THIS BLOCK OF  C
C     BASIS FUNCTIONS -- NOW OVERLAY THE RESULTS INTO OVAP.           C
C*********************************************************************C
C
C     CALCULATE COMPONENT OFFSETS
      IL1 = LARGE(ICNTA,KA,MJA  )
      IL2 = LARGE(ICNTA,KA,MJA+1)
      JL1 = LARGE(ICNTB,KB,MJB  )
      JL2 = LARGE(ICNTB,KB,MJB+1)
C
      IS1 = IL1 + NSHIFT
      IS2 = IL2 + NSHIFT
      JS1 = JL1 + NSHIFT
      JS2 = JL2 + NSHIFT
C
C     LL BLOCKS
      IF(IL1.GT.JL1) THEN
        DO J=1,NFUNB
          DO I=1,NFUNA
            OVAP(IL1+I,JL1+J) = SLL(I,J,1)
            OVAP(IL1+I,JL2+J) = SLL(I,J,2)
            OVAP(IL2+I,JL1+J) = SLL(I,J,3)
            OVAP(IL2+I,JL2+J) = SLL(I,J,4)
            OVAP(JL1+J,IL1+I) = DCONJG(OVAP(IL1+I,JL1+J))
            OVAP(JL2+J,IL1+I) = DCONJG(OVAP(IL1+I,JL2+J))
            OVAP(JL1+J,IL2+I) = DCONJG(OVAP(IL2+I,JL1+J))
            OVAP(JL2+J,IL2+I) = DCONJG(OVAP(IL2+I,JL2+J))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO J=1,NFUNB
          DO I=J,NFUNA
            OVAP(IL1+I,JL1+J) = SLL(I,J,1)
            OVAP(IL1+I,JL2+J) = SLL(I,J,2)
            OVAP(IL2+I,JL1+J) = SLL(I,J,3)
            OVAP(IL2+I,JL2+J) = SLL(I,J,4)
            OVAP(JL1+J,IL1+I) = DCONJG(OVAP(IL1+I,JL1+J))
            OVAP(JL2+J,IL1+I) = DCONJG(OVAP(IL1+I,JL2+J))
            OVAP(JL1+J,IL2+I) = DCONJG(OVAP(IL2+I,JL1+J))
            OVAP(JL2+J,IL2+I) = DCONJG(OVAP(IL2+I,JL2+J))
          ENDDO
        ENDDO
      ENDIF
C
C     SS BLOCKS
      IF(IL1.GT.JL1) THEN
        DO J=1,NFUNB
          DO I=1,NFUNA
            OVAP(IS1+I,JS1+J) = SSS(I,J,1)
            OVAP(IS1+I,JS2+J) = SSS(I,J,2)
            OVAP(IS2+I,JS1+J) = SSS(I,J,3)
            OVAP(IS2+I,JS2+J) = SSS(I,J,4)
            OVAP(JS1+J,IS1+I) = DCONJG(OVAP(IS1+I,JS1+J))
            OVAP(JS2+J,IS1+I) = DCONJG(OVAP(IS1+I,JS2+J))
            OVAP(JS1+J,IS2+I) = DCONJG(OVAP(IS2+I,JS1+J))
            OVAP(JS2+J,IS2+I) = DCONJG(OVAP(IS2+I,JS2+J))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO J=1,NFUNB
          DO I=J,NFUNA
            OVAP(IS1+I,JS1+J) = SSS(I,J,1)
            OVAP(IS1+I,JS2+J) = SSS(I,J,2)
            OVAP(IS2+I,JS1+J) = SSS(I,J,3)
            OVAP(IS2+I,JS2+J) = SSS(I,J,4)
            OVAP(JS1+J,IS1+I) = DCONJG(OVAP(IS1+I,JS1+J))
            OVAP(JS2+J,IS1+I) = DCONJG(OVAP(IS1+I,JS2+J))
            OVAP(JS1+J,IS2+I) = DCONJG(OVAP(IS2+I,JS1+J))
            OVAP(JS2+J,IS2+I) = DCONJG(OVAP(IS2+I,JS2+J))
          ENDDO
        ENDDO
      ENDIF
C
2000  CONTINUE
C
      RETURN
      END
