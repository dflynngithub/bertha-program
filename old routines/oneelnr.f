      SUBROUTINE ONEELNR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C     OOOOOO  NN    NN EEEEEEEE EEEEEEEE LL      NN    NN RRRRRRR     C
C    OO    OO NNN   NN EE       EE       LL      NNN   NN RR    RR    C
C    OO    OO NNNN  NN EE       EE       LL      NNNN  NN RR    RR    C
C    OO    OO NN NN NN EEEEEE   EEEEEE   LL      NN NN NN RR    RR    C
C    OO    OO NN  NNNN EE       EE       LL      NN  NNNN RRRRRRR     C
C    OO    OO NN   NNN EE       EE       LL      NN   NNN RR    RR    C
C     OOOOOO  NN    NN EEEEEEEE EEEEEEEE LLLLLLL NN    NN RR    RR    C
C                                                                     C
C ------------------------------------------------------------------- C
C     ONEELNR CONSTRUCTS NON-REL ONE-ELECTRON FOCK MULTI-CENTRE       C
C     INTEGRALS AND GIVES THE ONE-ELECTRON ENERGY OVER ATOMIC DENSITY.C
C*********************************************************************C
      PARAMETER(MXDM=1600,NCNTM=15,MAXB=26,MAXB2=MAXB*MAXB,NKAPM=9,
     &          LMAX=(NKAPM-1)/2,LAMMX=2*LMAX,IL4=2*LAMMX,MAXMV=NKAPM,
     &          MAXMV2=2*MAXMV,MLL=((NKAPM+8)*(NKAPM+9)*(NKAPM+10))/6,
     &          MAXR=969)
C
      COMPLEX*16 HMAT(MXDM,MXDM),GMAT(MXDM,MXDM),QMAT(MXDM,MXDM),
     &           BMAT(MXDM,MXDM),OVAP(MXDM,MXDM),FOCK(MXDM,MXDM)
      COMPLEX*16 DENC(MXDM,MXDM),DENO(MXDM,MXDM),DENT(MXDM,MXDM)
      COMPLEX*16 HLL(MAXB,MAXB,4),SLL(MAXB,MAXB,4)
      COMPLEX*16 E11(MAXB2,MLL),E21(MAXB2,MLL)
      COMPLEX*16 E11A,E11B,E11C,TRM11,E21A,E21B,E21C,TRM21,ETMP
C
      COMMON/ACSS/INABCD(0:IL4,0:IL4,0:IL4),
     &            IVEC(MLL),JVEC(MLL),KVEC(MLL),LAMVEC(MLL)
      COMMON/DENS/DENC,DENO,DENT
      COMMON/ENRG/ETOT,ENUC,EONE,ECLM,EBRT
      COMMON/LBLS/LARGE(NCNTM,NKAPM,MAXMV2)
      COMMON/LGHT/CV
      COMMON/MTRX/HMAT,GMAT,QMAT,BMAT,OVAP,FOCK
      COMMON/SPEC/EXPSET(MAXB,NKAPM,NCNTM),COORD(3,NCNTM),ZNUC(NCNTM),
     &            AMASS(NCNTM),CNUC(NCNTM),NFUNCT(NKAPM,NCNTM),
     &            IZNUC(NCNTM),IQNUC(NCNTM),LMAXX(NCNTM),NKAP(NCNTM),
     &            KVALS(NKAPM,NCNTM),NCNT,NDIM,NSHIFT,NOCC
      COMMON/XYZ4/XYZ(3,4)
C
      DIMENSION RC(MAXB2,MAXR),EXPT(MAXB,4),APH(MAXB2),CP(MAXB2,3),
     &          PNC(MAXB2),KQN(4),LQN(4),MQN(4),NFUNS(4)
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
        KQN(2) = KVALS(KB,ICNTB)
        IF(KQN(2).GT.0) THEN
          LQN(2) = KQN(2)
        ELSE
          LQN(2) =-KQN(2)-1
        ENDIF
C
        NFUNB    = NFUNCT(LQN(2)+1,ICNTB)
        NFUNS(2) = NFUNB
C
        DO JBAS=1,NFUNB
          EXPT(JBAS,2) = EXPSET(JBAS,LQN(2)+1,ICNTB)
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
        IF(MQN(1).NE.MQN(2)) GOTO 1999
      ENDIF
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
          EIJ    = EXPT(I,1) + EXPT(J,2)
          EROOT  = DSQRT(PI/EIJ)**3
          SLL(I,J,1) = EROOT*E11(M,1)
          SLL(I,J,3) = EROOT*E21(M,1)
          SLL(I,J,2) =-FASE*DCONJG(SLL(I,J,3))
          SLL(I,J,4) = FASE*DCONJG(SLL(I,J,1))
        ENDDO
      ENDDO
C
C     CONSTRUCT OVERLAP AND KINETIC ENERGY INTEGRALS (IOS 91)
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
            EIJ = EXPT(I,1) + EXPT(J,2)
            ESM = CNUC(IZ) + EIJ
            PX  = (XYZ(1,1)*EXPT(I,1) + XYZ(1,2)*EXPT(J,2))/EIJ
            PY  = (XYZ(2,1)*EXPT(I,1) + XYZ(2,2)*EXPT(J,2))/EIJ
            PZ  = (XYZ(3,1)*EXPT(I,1) + XYZ(3,2)*EXPT(J,2))/EIJ
            APH(M) = EIJ*CNUC(IZ)/ESM
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
            HMAT(JL1+J,IL2+I) = DCONJG(HMAT(IL2+I,JL1+J))
            HMAT(JL2+J,IL1+I) = DCONJG(HMAT(IL1+I,JL2+J))
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
      IF(IL1.EQ.JL1) THEN
        DO J=1,NFUNB
          DO I=J,NFUNA
            HMAT(IL1+I,JL1+J) = HMAT(IL1+I,JL1+J) + HLL(I,J,1)
            HMAT(IL1+I,JL2+J) = HMAT(IL1+I,JL2+J) + HLL(I,J,2)
            HMAT(IL2+I,JL1+J) = HMAT(IL2+I,JL1+J) + HLL(I,J,3)
            HMAT(IL2+I,JL2+J) = HMAT(IL2+I,JL2+J) + HLL(I,J,4)
C
            HMAT(JL1+J,IL1+I) = DCONJG(HMAT(IL1+I,JL1+J))
            HMAT(JL1+J,IL2+I) = DCONJG(HMAT(IL2+I,JL1+J))
            HMAT(JL2+J,IL1+I) = DCONJG(HMAT(IL1+I,JL2+J))
            HMAT(JL2+J,IL2+I) = DCONJG(HMAT(IL2+I,JL2+J))
C
            OVAP(IL1+I,JL1+J) = SLL(I,J,1)
            OVAP(IL1+I,JL2+J) = SLL(I,J,2)
            OVAP(IL2+I,JL1+J) = SLL(I,J,3)
            OVAP(IL2+I,JL2+J) = SLL(I,J,4)
C
            OVAP(JL1+J,IL1+I) = DCONJG(OVAP(IL1+I,JL1+J))
            OVAP(JL1+J,IL2+I) = DCONJG(OVAP(IL2+I,JL1+J))
            OVAP(JL2+J,IL1+I) = DCONJG(OVAP(IL1+I,JL2+J))
            OVAP(JL2+J,IL2+I) = DCONJG(OVAP(IL2+I,JL2+J))
          ENDDO
        ENDDO
      ENDIF
C
C*********************************************************************C
C     CALCULATE THE ONE ELECTRON CONTRIBUTION TO THE TOTAL ENERGY     C
C*********************************************************************C
C
C     INITIALISE TEMPORARY ENERGY COUNTER FOR THIS BLOCK
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
1999  CONTINUE
C
2000  CONTINUE
C
      RETURN
      END

