      SUBROUTINE ONEEL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C               OOOOOO  NN    NN EEEEEEE EEEEEEE LL                    C
C              OO    OO NNN   NN EE      EE      LL                    C
C              OO    OO NNNN  NN EE      EE      LL                    C
C              OO    OO NN NN NN EEEEEE  EEEEEE  LL                    C
C              OO    OO NN  NNNN EE      EE      LL                    C
C              OO    OO NN   NNN EE      EE      LL                    C
C               OOOOOO  NN    NN EEEEEEE EEEEEEE LLLLLLL               C
C                                                                      C
C -------------------------------------------------------------------- C
C     ONEEL CONSTRUCTS THE OVERLAP AND ONE-ELECTRON MULTI-CENTRE       C
C     MATRICES AND CALCULATES THE ONE-ELECTRON ENERGY.                 C
C**********************************************************************C
      PARAMETER(MDM=1600,MCT=15,MKP=9,MMV=MKP,MBS=26,MB2=MBS*MBS,
     &          MLL=MKP*(MKP+1)*(MKP+2)/6,IL4=2*(MKP-1),
     &          MRC=(IL4+1)*(IL4+2)*(IL4+3)/6)
C
      CHARACTER*4 HMLTN
C
      COMPLEX*16 OVAP(MDM,MDM),HNUC(MDM,MDM),HKIN(MDM,MDM),
     &           GDIR(MDM,MDM),GXCH(MDM,MDM),QDIR(MDM,MDM),
     &           QXCH(MDM,MDM),BDIR(MDM,MDM),BXCH(MDM,MDM),FOCK(MDM,MDM)
      COMPLEX*16 DENC(MDM,MDM),DENO(MDM,MDM),DENT(MDM,MDM)
C
      COMPLEX*16 ELLAB11(MB2,4*MLL*MLL),ELLAB21(MB2,4*MLL*MLL),
     &           ESSAB11(MB2,4*MLL*MLL),ESSAB21(MB2,4*MLL*MLL)
C
      COMPLEX*16 SLL(MBS,MBS,4),SSS(MBS,MBS,4),
     &           VLL(MBS,MBS,4),VSS(MBS,MBS,4),
     &           TLL(MBS,MBS,4),TLS(MBS,MBS,4),TSL(MBS,MBS,4)
      COMPLEX*16 E11(MB2,MLL),E21(MB2,MLL)
      COMPLEX*16 E11A,E11B,E11C,TRM11,E21A,E21B,E21C,TRM21
      COMPLEX*16 CTMP1,CTMP2,CTMP3,CTMP4

      DIMENSION RC(MB2,MRC),CP(MB2,3),XYZ(3,4),APH(MB2),PNC(MB2),
     &          EXPT(MBS,4),KQN(4),LQN(4),MQN(4),NFUNS(4)
C
      COMMON/ABLL/ELLAB11,ELLAB21
      COMMON/ABSS/ESSAB11,ESSAB21
      COMMON/ACSS/INABCD(0:IL4,0:IL4,0:IL4),
     &            IVEC(MRC),JVEC(MRC),KVEC(MRC),LAMVEC(MRC)
      COMMON/DENS/DENC,DENO,DENT
      COMMON/ILLM/ILLAD(MCT,MCT,MKP,MKP,MKP,MKP),IABLL,ICDLL
      COMMON/ISSM/ISSAD(MCT,MCT,MKP,MKP,MKP,MKP),IABSS,ICDSS
      COMMON/MTRX/OVAP,HNUC,HKIN,GDIR,GXCH,QDIR,QXCH,BDIR,BXCH,FOCK
      COMMON/PRMS/CV,HMLTN,ITER,IALL,IRUN,IEQS
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,2*MMV),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVIR,IOCCM0
      COMMON/TIME/TELL,TESS,TELS,TRLL,TRSS,TRLS,TRBR
C
      DATA PI/3.1415926535897932D0/
C
C     INITIALISE STORAGE MATRICES
      DO I=1,NDIM
        DO J=1,NDIM
          HNUC(I,J) = DCMPLX(0.0D0,0.0D0)
          HKIN(I,J) = DCMPLX(0.0D0,0.0D0)
          OVAP(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)      C
C**********************************************************************C
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
          INUCAB = 0
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
C**********************************************************************C
C     IMPLEMENT ANY AVAILABLE SELECTION RULES HERE                     C
C**********************************************************************C
C
C     SELECTION RULES TO BE MADE BASED ON GEOMETRIC SYMMETRY,
C     ATOMIC COORDINATES AND QUANTUM NUMBER PAIRS. THE IDEA IS TO
C     SKIP CALCULATIONS THAT SHOULD PRODUCE NO OR NEGLIGIBLE EFFECT.
      IF(NCNT.LE.2) THEN
        IF(MQN(1).NE.MQN(2)) GOTO 2000
      ENDIF
C
C**********************************************************************C
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS           C
C     OF (MA,MB). FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE          C
C     ORDERED                                                          C
C     11: = (-|MQN(A)|,-|MQN(B)|) -> 1                                 C
C     12: = (-|MQN(A)|,+|MQN(B)|) -> 2                                 C
C     21: = (+|MQN(A)|,-|MQN(B)|) -> 3                                 C
C     22: = (+|MQN(A)|,+|MQN(B)|) -> 4                                 C
C**********************************************************************C
C
C     CONSTRUCTION OF ONE-ELECTRON MATRICES BY TT' BLOCKS...
C
C     THE PHASE GENERATES E022 AND E012 COEFFS FROM E011 AND E021
      FASE = DFLOAT((-1)**((MQN(1)-MQN(2))/2))
     &                    *DFLOAT((KQN(1)*KQN(2))/IABS(KQN(1)*KQN(2)))
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXM = NFUNA*NFUNB
C
C     INITIALISE STORAGE ARRAYS
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          DO IB=1,4
            SLL(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
            SSS(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
            TLL(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
            TLS(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
            TSL(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
            VLL(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
            VSS(IBAS,JBAS,IB) = DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     PART 1: THE LL MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)
      NTUVLL = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ELL0 COEFFICIENTS
      IALT = 1
      CALL CPU_TIME(TBEG)
      IF(IEQS.EQ.0) THEN
        CALL EMAKELL(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,1,2,0)
      ELSEIF(IEQS.EQ.1) THEN
        IABLL = ILLAD(ICNTA,ICNTB,KA,KB,MA,MB)
        DO ITUV=1,NTUVLL
          DO M=1,MAXM
            E11(M,ITUV) = ELLAB11(M,IABLL+ITUV)
            E21(M,ITUV) = ELLAB21(M,IABLL+ITUV)
          ENDDO
        ENDDO
      ENDIF
      CALL CPU_TIME(TFIN)
      TELL = TELL + TFIN - TBEG
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          EIJ    = EXPT(IBAS,1)+EXPT(JBAS,2)
          EROOT  = DSQRT(PI/EIJ)**3
          SLL(IBAS,JBAS,1) = EROOT*E11(M,1)
          SLL(IBAS,JBAS,3) = EROOT*E21(M,1)
          SLL(IBAS,JBAS,2) =-FASE*DCONJG(SLL(IBAS,JBAS,3))
          SLL(IBAS,JBAS,4) = FASE*DCONJG(SLL(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     NUCLEAR ATTRACTION MATRIX ELEMENTS
      DO IZ=1,NCNT
C
C       NUCLEAR COORDINATES
        CX = COORD(1,IZ)
        CY = COORD(2,IZ)
        CZ = COORD(3,IZ)
C
C       GAUSSIAN PRODUCT THEOREM OVER BASIS FUNCTIONS
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M   = M+1
            EIJ = EXPT(IBAS,1)+EXPT(JBAS,2)
            ESM = CNUC(IZ)+EIJ
            PX  = (XYZ(1,1)*EXPT(IBAS,1) + XYZ(1,2)*EXPT(JBAS,2))/EIJ
            PY  = (XYZ(2,1)*EXPT(IBAS,1) + XYZ(2,2)*EXPT(JBAS,2))/EIJ
            PZ  = (XYZ(3,1)*EXPT(IBAS,1) + XYZ(3,2)*EXPT(JBAS,2))/EIJ
            APH(M) = (EIJ*CNUC(IZ))/ESM
            PNC(M) = 2.0D0*PI*DSQRT(CNUC(IZ)/ESM)*ZNUC(IZ)/EIJ
            CP(M,1) = CX - PX
            CP(M,2) = CY - PY
            CP(M,3) = CZ - PZ
          ENDDO
        ENDDO
C
C       GENERATE A BATCH OF R-INTEGRALS
        CALL CPU_TIME(TBEG)
        CALL RMAKE(RC,CP,APH,MAXM,LAM)
        CALL CPU_TIME(TFIN)
        TRLL = TRLL + TFIN - TBEG
C
C       NUCLEAR ATTRACTION INTEGRALS AS A FINITE SUME OF ELL0 AND RC
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M = M+1
            DO ITUV=1,NTUVLL
              VLL(IBAS,JBAS,1) = VLL(IBAS,JBAS,1)
     &                         - PNC(M)*E11(M,ITUV)*RC(M,ITUV)
              VLL(IBAS,JBAS,3) = VLL(IBAS,JBAS,3) 
     &                         - PNC(M)*E21(M,ITUV)*RC(M,ITUV)
            ENDDO
            VLL(IBAS,JBAS,2) =-FASE*DCONJG(VLL(IBAS,JBAS,3))
            VLL(IBAS,JBAS,4) = FASE*DCONJG(VLL(IBAS,JBAS,1))
          ENDDO
        ENDDO
      ENDDO
C
C     CONSTRUCT NON-REL KINETIC ENERGY INTEGRALS (IOS 91)
      IF(HMLTN.EQ.'NORL') THEN
        RL2 = DFLOAT(2*LQN(2)+3)
        M   = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M     = M+1
            EJ    = EXPT(JBAS,2)
            EIJ   = EXPT(IBAS,1) + EXPT(JBAS,2)
            EROOT = DSQRT(PI/EIJ)**3
            PX = (XYZ(1,1)*EXPT(IBAS,1) + XYZ(1,2)*EXPT(JBAS,2))/EIJ
            PY = (XYZ(2,1)*EXPT(IBAS,1) + XYZ(2,2)*EXPT(JBAS,2))/EIJ
            PZ = (XYZ(3,1)*EXPT(IBAS,1) + XYZ(3,2)*EXPT(JBAS,2))/EIJ
C
            PBX = PX - XYZ(1,2)
            PBY = PY - XYZ(2,2)
            PBZ = PZ - XYZ(3,2)
            PB2 = PBX*PBX + PBY*PBY + PBZ*PBZ
C       
            E0FC = EJ*RL2 - 2.0D0*EJ*EJ*PB2 - 3.00D0*EJ*EJ/EIJ
            E1FC = 4.0D0*EJ*EJ
C
C           TRUNCATE EXPRESSION DEPENDING ON LAM VALUE
C           ALL COMBINATIONS ALLOW FOR THE LAM = 0 MANIFOLD           
            TRM11 = E0FC*E11(M,1)
            TRM21 = E0FC*E21(M,1)
C           IF LAM > 0 PROVIDE SECOND BUNCH OF TERMS
            IF(LAM.GE.1) THEN
              E11A = E11(M,INABCD(1,0,0))
              E21A = E21(M,INABCD(1,0,0))
              E11B = E11(M,INABCD(0,1,0))
              E21B = E21(M,INABCD(0,1,0))
              E11C = E11(M,INABCD(0,0,1))
              E21C = E21(M,INABCD(0,0,1))
              TRM11 = TRM11 - E1FC*(PBX*E11A + PBY*E11B + PBZ*E11C)
              TRM21 = TRM21 - E1FC*(PBX*E21A + PBY*E21B + PBZ*E21C)
            ENDIF
C           IF LAM > 1 PROVIDE FINAL BUNCH OF TERMS
            IF(LAM.GE.2) THEN
              E11A = E11(M,INABCD(2,0,0))
              E21A = E21(M,INABCD(2,0,0))
              E11B = E11(M,INABCD(0,2,0))
              E21B = E21(M,INABCD(0,2,0))
              E11C = E11(M,INABCD(0,0,2))
              E21C = E21(M,INABCD(0,0,2))
              TRM11 = TRM11 - E1FC*(E11A + E11B + E11C)
              TRM21 = TRM21 - E1FC*(E21A + E21B + E21C)
            ENDIF
            TLL(IBAS,JBAS,1) = EROOT*TRM11
            TLL(IBAS,JBAS,3) = EROOT*TRM21
            TLL(IBAS,JBAS,2) =-FASE*DCONJG(TLL(IBAS,JBAS,3))
            TLL(IBAS,JBAS,4) = FASE*DCONJG(TLL(IBAS,JBAS,1))
          ENDDO
        ENDDO
C       NON-REL HAMILTONIAN MATRICES COMPLETE
        GOTO 500
      ENDIF
C     
C**********************************************************************C
C     PART 2: THE SS MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)+2
      NTUVSS = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ESS0 COEFFICIENTS
      IALT = 1
      CALL CPU_TIME(TBEG)
      IF(IEQS.EQ.0) THEN
        CALL EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,1,2,0)
      ELSEIF(IEQS.EQ.1) THEN
        IABSS = ISSAD(ICNTA,ICNTB,KA,KB,MA,MB)
        DO ITUV=1,NTUVSS
          DO M=1,MAXM
            E11(M,ITUV) = ESSAB11(M,IABSS+ITUV)
            E21(M,ITUV) = ESSAB21(M,IABSS+ITUV)
          ENDDO
        ENDDO
      ENDIF
      CALL CPU_TIME(TFIN)
      TESS = TESS + TFIN - TBEG
C
C     OVERLAP MATRIX ELEMENTS
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          EIJ     = EXPT(IBAS,1)+EXPT(JBAS,2)
          EROOT   = DSQRT(PI/EIJ)**3
          SSS(IBAS,JBAS,1) = EROOT*E11(M,1)
          SSS(IBAS,JBAS,3) = EROOT*E21(M,1)
          SSS(IBAS,JBAS,2) =-FASE*DCONJG(SSS(IBAS,JBAS,3))
          SSS(IBAS,JBAS,4) = FASE*DCONJG(SSS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     NUCLEAR ATTRACTION MATRIX ELEMENTS
      DO IZ=1,NCNT
C
C       NUCLEAR COORDINATES
        CX = COORD(1,IZ)
        CY = COORD(2,IZ)
        CZ = COORD(3,IZ)
C
C       GAUSSIAN PRODUCT THEOREM OVER BASIS FUNCTIONS
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M   = M+1
            EIJ = EXPT(IBAS,1)+EXPT(JBAS,2)
            ESM = CNUC(IZ) + EIJ
            PX  = (XYZ(1,1)*EXPT(IBAS,1) + XYZ(1,2)*EXPT(JBAS,2))/EIJ
            PY  = (XYZ(2,1)*EXPT(IBAS,1) + XYZ(2,2)*EXPT(JBAS,2))/EIJ
            PZ  = (XYZ(3,1)*EXPT(IBAS,1) + XYZ(3,2)*EXPT(JBAS,2))/EIJ
            APH(M) = (EIJ*CNUC(IZ))/ESM
            PNC(M) = 2.0D0*PI*DSQRT(CNUC(IZ)/ESM)*ZNUC(IZ)/EIJ
            CP(M,1) = CX - PX
            CP(M,2) = CY - PY
            CP(M,3) = CZ - PZ
          ENDDO
        ENDDO      
C
C       GENERATE A BATCH OF R-INTEGRALS
        CALL CPU_TIME(TBEG)
        CALL RMAKE(RC,CP,APH,MAXM,LAM)
        CALL CPU_TIME(TFIN)
        TRSS = TRSS + TFIN - TBEG
C
C       NUCLEAR ATTRACTION INTEGRALS AS A FINITE SUME OF ESS0 AND RC
        M = 0
        DO IBAS=1,NFUNA
          DO JBAS=1,NFUNB
            M = M+1
            DO ITUV=1,NTUVSS
              VSS(IBAS,JBAS,1) = VSS(IBAS,JBAS,1) 
     &                           + PNC(M)*E11(M,ITUV)*RC(M,ITUV)
              VSS(IBAS,JBAS,3) = VSS(IBAS,JBAS,3) 
     &                           + PNC(M)*E21(M,ITUV)*RC(M,ITUV)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C     SUBTRACT THE SS OVERLAP MATRIX AND FINISH CONSTRUCTION
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          VSS(IBAS,JBAS,1) =-VSS(IBAS,JBAS,1) 
     &                       - 2.0D0*CV*CV*SSS(IBAS,JBAS,1)
          VSS(IBAS,JBAS,3) =-VSS(IBAS,JBAS,3) 
     &                       - 2.0D0*CV*CV*SSS(IBAS,JBAS,3)
          VSS(IBAS,JBAS,2) =-FASE*DCONJG(VSS(IBAS,JBAS,3))
          VSS(IBAS,JBAS,4) = FASE*DCONJG(VSS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C**********************************************************************C
C     PART 3: THE SL MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)+2
      NTUVSS = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     KINETIC MATRIX ELEMENTS
      FACT = CV*DSQRT(DFLOAT(2*LQN(2)+3))
      M = 0
      DO IBAS=1,NFUNA
        DO JBAS=1,NFUNB
          M = M+1
          EIJ    = EXPT(IBAS,1) + EXPT(JBAS,2)
          EJRT   = FACT*DSQRT(EXPT(JBAS,2))
          EROOT  = DSQRT(PI/EIJ)**3
          TSL(IBAS,JBAS,1) = EJRT*EROOT*E11(M,1)
          TSL(IBAS,JBAS,3) = EJRT*EROOT*E21(M,1)
          TSL(IBAS,JBAS,2) =-FASE*DCONJG(TSL(IBAS,JBAS,3))
          TSL(IBAS,JBAS,4) = FASE*DCONJG(TSL(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C
C**********************************************************************C
C     PART 4: THE LS MATRICES                                          C
C**********************************************************************C
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAM    = LQN(1)+LQN(2)+2
      NTUVSS = (LAM+1)*(LAM+2)*(LAM+3)/6
C
C     GENERATE ESS0 COEFFICIENTS
      IALT = 1
C
      CALL CPU_TIME(TBEG)
      IF(IEQS.EQ.0) THEN
        CALL EMAKESS(E11,E21,EXPT,XYZ,KQN,MQN,NFUNS,IALT,2,1,0)
      ELSEIF(IEQS.EQ.1) THEN
        IABSS = ISSAD(ICNTB,ICNTA,KB,KA,MB,MA)
        DO ITUV=1,NTUVSS
          DO M=1,MAXM
            E11(M,ITUV) = ESSAB11(M,IABSS+ITUV)
            E21(M,ITUV) = ESSAB21(M,IABSS+ITUV)
          ENDDO
        ENDDO
      ENDIF
      CALL CPU_TIME(TFIN)
      TESS = TESS + TFIN - TBEG
C
C     KINETIC MATRIX ELEMENTS
      FACT = CV*DSQRT(DFLOAT(2*LQN(1)+3))
      M = 0
      DO JBAS=1,NFUNB
        DO IBAS=1,NFUNA
          M = M+1
          EIJ    = EXPT(JBAS,2) + EXPT(IBAS,1)
          EIRT   = FACT*DSQRT(EXPT(IBAS,1))
          EROOT  = DSQRT(PI/EIJ)**3
          TLS(IBAS,JBAS,1) = EIRT*EROOT*E11(M,1)
          TLS(IBAS,JBAS,3) = EIRT*EROOT*E21(M,1)
          TLS(IBAS,JBAS,2) =-FASE*DCONJG(TLS(IBAS,JBAS,3))
          TLS(IBAS,JBAS,4) = FASE*DCONJG(TLS(IBAS,JBAS,1))
        ENDDO
      ENDDO
C
C     GENERATE LS MATRICES FROM THE ABOVE SL MATRICES
      M = 0
      DO JBAS=1,NFUNB
        DO IBAS=1,NFUNA
          M = M + 1
          CTMP1 = DCONJG(TLS(IBAS,JBAS,1))
          CTMP2 = DCONJG(TLS(IBAS,JBAS,2))
          CTMP3 = DCONJG(TLS(IBAS,JBAS,3))
          CTMP4 = DCONJG(TLS(IBAS,JBAS,4))
C
          TLS(IBAS,JBAS,1) = CTMP1
          TLS(IBAS,JBAS,2) = CTMP3
          TLS(IBAS,JBAS,3) = CTMP2
          TLS(IBAS,JBAS,4) = CTMP4
        ENDDO
      ENDDO
C
500   CONTINUE
C
C**********************************************************************C
C     WE NOW HAVE ALL PIECES OF HNUC AND HKIN FOR THIS BLOCK.          C
C**********************************************************************C
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
C     LL OVERLAP BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            OVAP(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            OVAP(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            OVAP(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            OVAP(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
C
            OVAP(JL1+JBAS,IL1+IBAS) = DCONJG(OVAP(IL1+IBAS,JL1+JBAS))
            OVAP(JL2+JBAS,IL1+IBAS) = DCONJG(OVAP(IL1+IBAS,JL2+JBAS))
            OVAP(JL1+JBAS,IL2+IBAS) = DCONJG(OVAP(IL2+IBAS,JL1+JBAS))
            OVAP(JL2+JBAS,IL2+IBAS) = DCONJG(OVAP(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            OVAP(IL1+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,1)
            OVAP(IL1+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,2)
            OVAP(IL2+IBAS,JL1+JBAS) = SLL(IBAS,JBAS,3)
            OVAP(IL2+IBAS,JL2+JBAS) = SLL(IBAS,JBAS,4)
C
            OVAP(JL1+JBAS,IL1+IBAS) = DCONJG(OVAP(IL1+IBAS,JL1+JBAS))
            OVAP(JL2+JBAS,IL1+IBAS) = DCONJG(OVAP(IL1+IBAS,JL2+JBAS))
            OVAP(JL1+JBAS,IL2+IBAS) = DCONJG(OVAP(IL2+IBAS,JL1+JBAS))
            OVAP(JL2+JBAS,IL2+IBAS) = DCONJG(OVAP(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     LL NUCLEAR POTENTIAL BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            HNUC(IL1+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,1)
            HNUC(IL1+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,2)
            HNUC(IL2+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,3)
            HNUC(IL2+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,4)
C
            HNUC(JL1+JBAS,IL1+IBAS) = DCONJG(HNUC(IL1+IBAS,JL1+JBAS))
            HNUC(JL2+JBAS,IL1+IBAS) = DCONJG(HNUC(IL1+IBAS,JL2+JBAS))
            HNUC(JL1+JBAS,IL2+IBAS) = DCONJG(HNUC(IL2+IBAS,JL1+JBAS))
            HNUC(JL2+JBAS,IL2+IBAS) = DCONJG(HNUC(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            HNUC(IL1+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,1)
            HNUC(IL1+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,2)
            HNUC(IL2+IBAS,JL1+JBAS) = VLL(IBAS,JBAS,3)
            HNUC(IL2+IBAS,JL2+JBAS) = VLL(IBAS,JBAS,4)
C
            HNUC(JL1+JBAS,IL1+IBAS) = DCONJG(HNUC(IL1+IBAS,JL1+JBAS))
            HNUC(JL2+JBAS,IL1+IBAS) = DCONJG(HNUC(IL1+IBAS,JL2+JBAS))
            HNUC(JL1+JBAS,IL2+IBAS) = DCONJG(HNUC(IL2+IBAS,JL1+JBAS))
            HNUC(JL2+JBAS,IL2+IBAS) = DCONJG(HNUC(IL2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     NON-REL HAMILTONIAN HAS A KINETIC MATRIX IN THE LL BLOCK
      IF(HMLTN.EQ.'NORL') THEN
C
C       LL KINETIC BLOCKS
        IF(IL1.GT.JL1) THEN
          DO JBAS=1,NFUNB
            DO IBAS=1,NFUNA
              HKIN(IL1+IBAS,JL1+JBAS) = TLL(IBAS,JBAS,1)
              HKIN(IL1+IBAS,JL2+JBAS) = TLL(IBAS,JBAS,2)
              HKIN(IL2+IBAS,JL1+JBAS) = TLL(IBAS,JBAS,3)
              HKIN(IL2+IBAS,JL2+JBAS) = TLL(IBAS,JBAS,4)
C
              HKIN(JL1+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JL1+JBAS))
              HKIN(JL2+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JL2+JBAS))
              HKIN(JL1+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JL1+JBAS))
              HKIN(JL2+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JL2+JBAS))
            ENDDO
          ENDDO
        ENDIF
C
        IF(IL1.EQ.JL1) THEN
          DO JBAS=1,NFUNB
            DO IBAS=JBAS,NFUNA
              HKIN(IL1+IBAS,JL1+JBAS) = TLL(IBAS,JBAS,1)
              HKIN(IL1+IBAS,JL2+JBAS) = TLL(IBAS,JBAS,2)
              HKIN(IL2+IBAS,JL1+JBAS) = TLL(IBAS,JBAS,3)
              HKIN(IL2+IBAS,JL2+JBAS) = TLL(IBAS,JBAS,4)
C
              HKIN(JL1+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JL1+JBAS))
              HKIN(JL2+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JL2+JBAS))
              HKIN(JL1+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JL1+JBAS))
              HKIN(JL2+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JL2+JBAS))
            ENDDO
          ENDDO
        ENDIF
C      
C       NON-REL MATRIX CONSTRUCTION COMPLETE
        GOTO 600
C
      ENDIF
C
C     SS OVERLAP BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            OVAP(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            OVAP(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            OVAP(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            OVAP(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
C
            OVAP(JS1+JBAS,IS1+IBAS) = DCONJG(OVAP(IS1+IBAS,JS1+JBAS))
            OVAP(JS2+JBAS,IS1+IBAS) = DCONJG(OVAP(IS1+IBAS,JS2+JBAS))
            OVAP(JS1+JBAS,IS2+IBAS) = DCONJG(OVAP(IS2+IBAS,JS1+JBAS))
            OVAP(JS2+JBAS,IS2+IBAS) = DCONJG(OVAP(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            OVAP(IS1+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,1)
            OVAP(IS1+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,2)
            OVAP(IS2+IBAS,JS1+JBAS) = SSS(IBAS,JBAS,3)
            OVAP(IS2+IBAS,JS2+JBAS) = SSS(IBAS,JBAS,4)
C
            OVAP(JS1+JBAS,IS1+IBAS) = DCONJG(OVAP(IS1+IBAS,JS1+JBAS))
            OVAP(JS2+JBAS,IS1+IBAS) = DCONJG(OVAP(IS1+IBAS,JS2+JBAS))
            OVAP(JS1+JBAS,IS2+IBAS) = DCONJG(OVAP(IS2+IBAS,JS1+JBAS))
            OVAP(JS2+JBAS,IS2+IBAS) = DCONJG(OVAP(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     SS NUCLEAR POTENTIAL BLOCK
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            HNUC(IS1+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,1)
            HNUC(IS1+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,2)
            HNUC(IS2+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,3)
            HNUC(IS2+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,4)
C
            HNUC(JS1+JBAS,IS1+IBAS) = DCONJG(HNUC(IS1+IBAS,JS1+JBAS))
            HNUC(JS2+JBAS,IS1+IBAS) = DCONJG(HNUC(IS1+IBAS,JS2+JBAS))
            HNUC(JS1+JBAS,IS2+IBAS) = DCONJG(HNUC(IS2+IBAS,JS1+JBAS))
            HNUC(JS2+JBAS,IS2+IBAS) = DCONJG(HNUC(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
      IF(IL1.EQ.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=JBAS,NFUNA
            HNUC(IS1+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,1)
            HNUC(IS1+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,2)
            HNUC(IS2+IBAS,JS1+JBAS) = VSS(IBAS,JBAS,3)
            HNUC(IS2+IBAS,JS2+JBAS) = VSS(IBAS,JBAS,4)
C
            HNUC(JS1+JBAS,IS1+IBAS) = DCONJG(HNUC(IS1+IBAS,JS1+JBAS))
            HNUC(JS2+JBAS,IS1+IBAS) = DCONJG(HNUC(IS1+IBAS,JS2+JBAS))
            HNUC(JS1+JBAS,IS2+IBAS) = DCONJG(HNUC(IS2+IBAS,JS1+JBAS))
            HNUC(JS2+JBAS,IS2+IBAS) = DCONJG(HNUC(IS2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     LS BLOCKS
      IF(IL1.GE.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            HKIN(IL1+IBAS,JS1+JBAS) = TLS(IBAS,JBAS,1)
            HKIN(IL1+IBAS,JS2+JBAS) = TLS(IBAS,JBAS,2)
            HKIN(IL2+IBAS,JS1+JBAS) = TLS(IBAS,JBAS,3)
            HKIN(IL2+IBAS,JS2+JBAS) = TLS(IBAS,JBAS,4)
C
            HKIN(JS1+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JS1+JBAS))
            HKIN(JS2+JBAS,IL1+IBAS) = DCONJG(HKIN(IL1+IBAS,JS2+JBAS))
            HKIN(JS1+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JS1+JBAS))
            HKIN(JS2+JBAS,IL2+IBAS) = DCONJG(HKIN(IL2+IBAS,JS2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
C     SL BLOCKS
      IF(IL1.GT.JL1) THEN
        DO JBAS=1,NFUNB
          DO IBAS=1,NFUNA
            HKIN(IS1+IBAS,JL1+JBAS) = TSL(IBAS,JBAS,1)
            HKIN(IS1+IBAS,JL2+JBAS) = TSL(IBAS,JBAS,2)
            HKIN(IS2+IBAS,JL1+JBAS) = TSL(IBAS,JBAS,3)
            HKIN(IS2+IBAS,JL2+JBAS) = TSL(IBAS,JBAS,4)
C
            HKIN(JL1+JBAS,IS1+IBAS) = DCONJG(HKIN(IS1+IBAS,JL1+JBAS))
            HKIN(JL2+JBAS,IS1+IBAS) = DCONJG(HKIN(IS1+IBAS,JL2+JBAS))
            HKIN(JL1+JBAS,IS2+IBAS) = DCONJG(HKIN(IS2+IBAS,JL1+JBAS))
            HKIN(JL2+JBAS,IS2+IBAS) = DCONJG(HKIN(IS2+IBAS,JL2+JBAS))
          ENDDO
        ENDDO
      ENDIF
C
600   CONTINUE      
C
C     END LOOPS OVER BASIS PAIRS A,B
2000  CONTINUE
C
      RETURN
      END

