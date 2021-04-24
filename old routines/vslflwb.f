C     THIS VERSION INCLUDES A TESTING ROUTINE WHICH ORDERS THE POSITIVE-
C     ENERGY SPECTRUM CONTRIBUTIONS BY BASIS FUNCTION, AND FOR A HE BASIS
C     SET WITH 48s 48p, COLLAPSES MQN CONTRIBUTIONS INTO THE SAME KQN GROUP.
C

      SUBROUTINE VSLFLWB(VLL,VLS,VSL,VSS,IOCC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C  VV    VV  SSSSSS  LL       FFFFFFFF LL      WW         WW BBBBBBB   C
C  VV    VV SS    SS LL       FF       LL      WW         WW BB    BB  C
C  VV    VV SS       LL       FF       LL      WW         WW BB    BB  C
C  VV    VV  SSSSSS  LL       FFFFFF   LL      WW    W    WW BBBBBBB   C
C   VV  VV        SS LL       FF       LL       WW  WWW  WW  BB    BB  C
C    VVVV   SS    SS LL       FF       LL        WWWW WWWW   BB    BB  C
C     VV     SSSSSS  LLLLLLLL FF       LLLLLLLL   WW   WW    BBBBBBB   C
C                                                                      C
C -------------------------------------------------------------------- C
C  VSLFLWB CONSTRUCTS A MATRIX OF (μ,T|σ_Q.Vslf|ν,T') ELECTRON SELF-   C
C  INTERACTION INTEGRALS OVER ALL BASIS FUNCTION PAIRS (BETHE METHOD). C
C -------------------------------------------------------------------- C
C  INPUT:                                                              C
C  ▶ IOCC    = OCCUPIED STATE (AUTOMATICALLY SHIFTED).                 C
C  ▶ ITT     = {1,2,3,4} -> {LL,LS,SL,SS} COMPONENT COMBINATION.       C
C  ▶ IQ      = {0,1,2,3} -> {0,X,Y,Z} THE COUPLING σ MATRIX.           C
C  ▶ IA1,IA2 = {1,2,3,4} -> BASIS FUNCTION OVERLAPS FOR EQ-COEFFS.     C
C**********************************************************************C
      INCLUDE 'parameters.h'
C
      CHARACTER*1 CHL,LLAB
      CHARACTER*5 NMDL
C
      DIMENSION RKNT(MDM)
      DIMENSION IORD(MDM)
      dimension tally(mbs)
C
      COMPLEX*16 TSM
      COMPLEX*16 CL,CS
      COMPLEX*16 EX(MDM,4),EY(MDM,4),EZ(MDM,4),ET(MDM,4)
      COMPLEX*16 VLL(MDM,MDM),VLS(MDM,MDM),VSL(MDM,MDM),VSS(MDM,MDM)
      COMPLEX*16 COEF(MDM,MDM)
      COMPLEX*16 OLSX(MDM,MDM),OLSY(MDM,MDM),OLSZ(MDM,MDM),
     &           OSLX(MDM,MDM),OSLY(MDM,MDM),OSLZ(MDM,MDM)
      COMPLEX*16 PINLSX(MDM,MDM),PINLSY(MDM,MDM),PINLSZ(MDM,MDM)
      COMPLEX*16 PINSLX(MDM,MDM),PINSLY(MDM,MDM),PINSLZ(MDM,MDM)
      COMPLEX*16 WLS(MDM,3),WSL(MDM,3)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,MFT),XNUC(MCT,MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EIGC/COEF
      COMMON/EIGE/EIGN(MDM)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/PHYS/CV,EMSS,UMSS,TMSS,PMSS,PRAD,CMPW,GFREE,GFRMI,WEIN
      COMMON/QNMS/LABICN(MDM),LABKQN(MDM),LABMQN(MDM)
      COMMON/SPEC/KQNLST(MDM),MQNLST(MDM),NQNLST(MDM)
      COMMON/OTTI/OLSX,OLSY,OLSZ,OSLX,OSLY,OSLZ
C
      DATA EPS/1.0D-2/
C
C     AMPLITUDE FOR LOW-ENERGY CONTRIBUTION
      ALW =-2.0D0/(3.0D0*PI*CV)
C
C     SET THE LOW-ENERGY CUTOFF TO M*CV FOR NOW
      CUTK = CV
      
      do i=1,mbs
        tally(i) = 0.0d0
      enddo
C
C     INITIALISE STORAGE ARRAY
      DO I=1,NDIM-NSKP
        DO J=1,NDIM-NSKP
          VLL(I,J) = DCMPLX(0.0D0,0.0D0)
          VLS(I,J) = DCMPLX(0.0D0,0.0D0)
          VSL(I,J) = DCMPLX(0.0D0,0.0D0)
          VSS(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
      
      if(iocc.ne.1) return
C
C     SERIES OF ENERGY DIFFERENCES
      EM = EIGN(IOCC+NSKP)
      DO N=1,NDIM-NSKP
        EDIF = EM-EIGN(N+NSKP)
        ARGL = (CUTK*CV-EDIF)/DABS(EDIF)
        IF(DABS(EDIF).GT.EPS) THEN
          RKNT(N) = EDIF*DLOG(ARGL)
        ELSE
          RKNT(N) = 0.0D0
        ENDIF
      ENDDO
C
      DO I=1,NDIM-NSKP
        DO J=1,NDIM-NSKP
C          IF(LABKQN(J).NE.-2) THEN
C            OLSX(I,J) = DCMPLX(0.0D0,0.0D0)
C            OLSY(I,J) = DCMPLX(0.0D0,0.0D0)
C            OLSZ(I,J) = DCMPLX(0.0D0,0.0D0)
C            OSLX(I,J) = DCMPLX(0.0D0,0.0D0)
C            OSLY(I,J) = DCMPLX(0.0D0,0.0D0)
C            OSLZ(I,J) = DCMPLX(0.0D0,0.0D0)
C          ENDIF
        ENDDO
      ENDDO
C
C     GENERATE MATRIX ELEMENTS <I|SIG|N> FOR ALL N
      DO I=1,NDIM-NSKP
C
C       LOOP OVER ALL POSITIVE-ENERGY STATES N
        DO N=1,NDIM-NSKP
C
C         INITIALISE THESE MATRIX ELEMENTS
          PINLSX(I,N) = DCMPLX(0.0D0,0.0D0)
          PINLSY(I,N) = DCMPLX(0.0D0,0.0D0)
          PINLSZ(I,N) = DCMPLX(0.0D0,0.0D0)
          PINSLX(I,N) = DCMPLX(0.0D0,0.0D0)
          PINSLY(I,N) = DCMPLX(0.0D0,0.0D0)
          PINSLZ(I,N) = DCMPLX(0.0D0,0.0D0)
C
C         CONTRACT OVER COEFFICIENTS AND RAW MATRIX ELEMENTS
          DO L=1,NDIM-NSKP
C
C           SMALL-COMPONENT COEFFICIENT
            CL = COEF(L     ,N+NSKP)
            CS = COEF(L+NSKP,N+NSKP)
C
C           <I,L|SIG_X|N,S>
            PINLSX(I,N) = PINLSX(I,N) + CS*OLSX(I,L)
            PINLSY(I,N) = PINLSY(I,N) + CS*OLSY(I,L)
            PINLSZ(I,N) = PINLSZ(I,N) + CS*OLSZ(I,L)
C
C           <I,S|SIG_X|N,L>
            PINSLX(I,N) = PINSLX(I,N) + CL*OSLX(I,L)
            PINSLY(I,N) = PINSLY(I,N) + CL*OSLY(I,L)
            PINSLZ(I,N) = PINSLZ(I,N) + CL*OSLZ(I,L)
C
          ENDDO
C
        ENDDO
      ENDDO
C
C     FORM THE DOT PRODUCT OF MATRICES AND MULTIPLY BY ENERGY TERM
      DO I=1,NDIM-NSKP
        DO J=1,NDIM-NSKP
C
          VLL(I,J) = DCMPLX(0.0D0,0.0D0)
          VLS(I,J) = DCMPLX(0.0D0,0.0D0)
          VSL(I,J) = DCMPLX(0.0D0,0.0D0)
          VSS(I,J) = DCMPLX(0.0D0,0.0D0)
C
C         LOOP OVER ALL POSITIVE-ENERGY STATES N
          DO N=1,NDIM-NSKP
C
C           ITT = 1: <I,L|SIG|N,S><N,S|SIG|I,L>
            TXX = PINLSX(I,N)*DCONJG(PINLSX(J,N))
            TYY = PINLSY(I,N)*DCONJG(PINLSY(J,N))
            TZZ = PINLSZ(I,N)*DCONJG(PINLSZ(J,N))
C
            VLL(I,J) = VLL(I,J) + ALW*RKNT(N)*(TXX+TYY+TZZ)
C
C           ITT = 2: <I,L|SIG|N,S><N,L|SIG|I,S>
            TXX = PINLSX(I,N)*DCONJG(PINSLX(J,N))
            TYY = PINLSY(I,N)*DCONJG(PINSLY(J,N))
            TZZ = PINLSZ(I,N)*DCONJG(PINSLZ(J,N))
C
            VLS(I,J) = VLS(I,J) + ALW*RKNT(N)*(TXX+TYY+TZZ)
C
C           ITT = 3: <I,S|SIG|N,L><N,S|SIG|I,L>
            TXX = PINSLX(I,N)*DCONJG(PINLSX(J,N))
            TYY = PINSLY(I,N)*DCONJG(PINLSY(J,N))
            TZZ = PINSLZ(I,N)*DCONJG(PINLSZ(J,N))
C
            VSL(I,J) = VSL(I,J) + ALW*RKNT(N)*(TXX+TYY+TZZ)
C
C           ITT = 4: <I,S|SIG|N,L><N,L|SIG|I,S>
            TXX = PINSLX(I,N)*DCONJG(PINSLX(J,N))
            TYY = PINSLY(I,N)*DCONJG(PINSLY(J,N))
            TZZ = PINSLZ(I,N)*DCONJG(PINSLZ(J,N))
C
            VSS(I,J) = VSS(I,J) + ALW*RKNT(N)*(TXX+TYY+TZZ)
C
          ENDDO
        ENDDO
      ENDDO
C
C**********************************************************************C
C     CODE FOR TESTING PURPOSES (II)                                   C
C**********************************************************************C
C
C     GENERATE MATRIX ELEMENTS <IOCC|SIG|N> FOR ALL N
      IF(IOCC.NE.1) RETURN
C
C     GENERATE MATRIX ELEMENTS <M|SIG|N> FOR ALL N
C
C     LOOP OVER STATES N
      DO N=1,NDIM-NSKP
C
C       INITIALISE MATRIX ELEMENT STORAGE
        WLS(N,1) = DCMPLX(0.0D0,0.0D0)
        WLS(N,2) = DCMPLX(0.0D0,0.0D0)
        WLS(N,3) = DCMPLX(0.0D0,0.0D0)
        WSL(N,1) = DCMPLX(0.0D0,0.0D0)
        WSL(N,2) = DCMPLX(0.0D0,0.0D0)
        WSL(N,3) = DCMPLX(0.0D0,0.0D0)
C
        DO I=1,NDIM-NSKP
C
C         SMALL-COMPONENT COEFFICIENT
          CL = COEF(I     ,IOCC+NSKP)
          CS = COEF(I+NSKP,IOCC+NSKP)
C
C         <IOCC,L|SIG_X|N,S>
          WLS(N,1) = WLS(N,1) + DCONJG(CL)*PINLSX(I,N)
          WLS(N,2) = WLS(N,2) + DCONJG(CL)*PINLSY(I,N)
          WLS(N,3) = WLS(N,3) + DCONJG(CL)*PINLSZ(I,N)
C
C         <IOCC,L|SIG_X|N,S>
          WSL(N,1) = WSL(N,1) + DCONJG(CS)*PINSLX(I,N)
          WSL(N,2) = WSL(N,2) + DCONJG(CS)*PINSLY(I,N)
          WSL(N,3) = WSL(N,3) + DCONJG(CS)*PINSLZ(I,N)
C
        ENDDO
      ENDDO
C
C     CALCULATE CONTRIBUTIONS TO THE SELF-INTERACTION ENERGY
      WRITE(6, *) 'IOCC = ',IOCC
      WRITE(7, *) 'IOCC = ',IOCC
      NLNES = 143
      DO ITT=4,4
        WRITE(6, *) 'ITT  = ',ITT
        WRITE(7, *) 'ITT  = ',ITT
        WRITE(6, *) REPEAT('=',NLNES)
        WRITE(7, *) REPEAT('=',NLNES)
        WRITE(6,81)
        WRITE(7,81)
        WRITE(6, *) REPEAT('=',NLNES)
        WRITE(7, *) REPEAT('=',NLNES)
        TSM = DCMPLX(0.0D0,0.0D0)
        DO N=1,NDIM-NSKP
          IF(ITT.EQ.1) THEN
            EX(N,ITT) = WLS(N,1)*DCONJG(WLS(N,1))
            EY(N,ITT) = WLS(N,2)*DCONJG(WLS(N,2))
            EZ(N,ITT) = WLS(N,3)*DCONJG(WLS(N,3))
          ELSEIF(ITT.EQ.2) THEN
            EX(N,ITT) = WLS(N,1)*DCONJG(WSL(N,1))
            EY(N,ITT) = WLS(N,2)*DCONJG(WSL(N,2))
            EZ(N,ITT) = WLS(N,3)*DCONJG(WSL(N,3))
          ELSEIF(ITT.EQ.3) THEN
            EX(N,ITT) = WSL(N,1)*DCONJG(WLS(N,1))
            EY(N,ITT) = WSL(N,2)*DCONJG(WLS(N,2))
            EZ(N,ITT) = WSL(N,3)*DCONJG(WLS(N,3))
          ELSEIF(ITT.EQ.4) THEN
            EX(N,ITT) = WSL(N,1)*DCONJG(WSL(N,1))
            EY(N,ITT) = WSL(N,2)*DCONJG(WSL(N,2))
            EZ(N,ITT) = WSL(N,3)*DCONJG(WSL(N,3))
          ENDIF          
          ET(N,ITT) = ALW*RKNT(N)*(EX(N,ITT)+EY(N,ITT)+EZ(N,ITT))
          TSM = TSM + ET(N,ITT)
        ENDDO
C
        TSM = DCMPLX(0.0D0,0.0D0)
        DO N=1,NDIM-NSKP
c         IF(ABS(ET(N,ITT)).LT.1.0D-14) GOTO 86
c          if(LABKQN(N).ne.-2) goto 86
          ital = mod(n,48)
          if(ital.eq.0) ital = 48
          tally(ital) = tally(ital) + DREAL(ET(N,ITT))
          IF(LABKQN(N).LT.0) THEN
            ILQN =-LABKQN(N)-1
          ELSE
            ILQN = LABKQN(N)
          ENDIF
          CHL = LLAB(ILQN)
          TSM = TSM + ET(N,ITT)
          WRITE(6,80) N,NQNLST(N),CHL,2*IABS(LABKQN(N))-1,LABMQN(N),
     &                EIGN(N+NSKP),RKNT(N),DREAL(EX(N,ITT)),
     &                DREAL(EY(N,ITT)),DREAL(EZ(N,ITT)),
     &                DREAL(ET(N,ITT)),DREAL(TSM)
86        continue
        ENDDO
        WRITE(6, *) REPEAT('=',NLNES)
        WRITE(7, *) REPEAT('=',NLNES)
        
C
      ENDDO
      
      bigsum = 0.0d0
      do ital=1,48
        write(*,*) ital,tally(ital)
        bigsum = bigsum + tally(ital)
      enddo
      WRITE(*,*) 'tot',bigsum
C
C**********************************************************************C
C     CODE FOR TESTING PURPOSES                                        C
C**********************************************************************C
C
      goto 888
C     GENERATE MATRIX ELEMENTS <IOCC|SIG|N> FOR ALL N
C
C     LOOP OVER STATES N
      DO N=1,NDIM-NSKP
C
C       INITIALISE MATRIX ELEMENT STORAGE
        WLS(N,1) = DCMPLX(0.0D0,0.0D0)
        WLS(N,2) = DCMPLX(0.0D0,0.0D0)
        WLS(N,3) = DCMPLX(0.0D0,0.0D0)
        WSL(N,1) = DCMPLX(0.0D0,0.0D0)
        WSL(N,2) = DCMPLX(0.0D0,0.0D0)
        WSL(N,3) = DCMPLX(0.0D0,0.0D0)
C
        DO I=1,NDIM-NSKP
C
C         SMALL-COMPONENT COEFFICIENT
          CL = COEF(I     ,IOCC+NSKP)
          CS = COEF(I+NSKP,IOCC+NSKP)
C
C         <IOCC,L|SIG_X|N,S>
          WLS(N,1) = WLS(N,1) + DCONJG(CL)*PINLSX(I,N)
          WLS(N,2) = WLS(N,2) + DCONJG(CL)*PINLSY(I,N)
          WLS(N,3) = WLS(N,3) + DCONJG(CL)*PINLSZ(I,N)
C
C         <IOCC,L|SIG_X|N,S>
          WSL(N,1) = WSL(N,1) + DCONJG(CS)*PINSLX(I,N)
          WSL(N,2) = WSL(N,2) + DCONJG(CS)*PINSLY(I,N)
          WSL(N,3) = WSL(N,3) + DCONJG(CS)*PINSLZ(I,N)
C
        ENDDO
      ENDDO
C
C     CALCULATE CONTRIBUTIONS TO THE SELF-INTERACTION ENERGY
      WRITE(6, *) 'IOCC = ',IOCC
      WRITE(7, *) 'IOCC = ',IOCC
80    FORMAT(1X,I3,2X,I4,A,'_',I1,'/2',2X,I2,'/2',2X,F19.10,2X,F17.10,
     &                         2X,'| ',3(F13.10,1X),'|',2(1X,ES17.10))
81    FORMAT(3X,'N',5X,'spinor',3X,'m_j',17X,'E(N)',12X,'RKNT(N)',2X,
     &     '|',4X,'<I|sX|N>^2',4X,'<I|sY|N>^2',4X,'<I|sZ|N>^2',1X,
     &     '|',13X,'SE(N)',13X,'TALLY')
      NLNES = 142
      DO ITT=1,4
        WRITE(6, *) 'ITT  = ',ITT
        WRITE(7, *) 'ITT  = ',ITT
        WRITE(6, *) REPEAT('=',NLNES)
        WRITE(7, *) REPEAT('=',NLNES)
        WRITE(6,81)
        WRITE(7,81)
        WRITE(6, *) REPEAT('=',NLNES)
        WRITE(7, *) REPEAT('=',NLNES)
        TSM = DCMPLX(0.0D0,0.0D0)
        DO N=1,NDIM-NSKP
          IF(ITT.EQ.1) THEN
            EX(N,ITT) = WLS(N,1)*DCONJG(WLS(N,1))
            EY(N,ITT) = WLS(N,2)*DCONJG(WLS(N,2))
            EZ(N,ITT) = WLS(N,3)*DCONJG(WLS(N,3))
          ELSEIF(ITT.EQ.2) THEN
            EX(N,ITT) = WLS(N,1)*DCONJG(WSL(N,1))
            EY(N,ITT) = WLS(N,2)*DCONJG(WSL(N,2))
            EZ(N,ITT) = WLS(N,3)*DCONJG(WSL(N,3))
          ELSEIF(ITT.EQ.3) THEN
            EX(N,ITT) = WSL(N,1)*DCONJG(WLS(N,1))
            EY(N,ITT) = WSL(N,2)*DCONJG(WLS(N,2))
            EZ(N,ITT) = WSL(N,3)*DCONJG(WLS(N,3))
          ELSEIF(ITT.EQ.4) THEN
            EX(N,ITT) = WSL(N,1)*DCONJG(WSL(N,1))
            EY(N,ITT) = WSL(N,2)*DCONJG(WSL(N,2))
            EZ(N,ITT) = WSL(N,3)*DCONJG(WSL(N,3))
          ENDIF          
          ET(N,ITT) = ALW*RKNT(N)*(EX(N,ITT)+EY(N,ITT)+EZ(N,ITT))
          TSM = TSM + ET(N,ITT)
        ENDDO
C      
C       ORGANISE RESULTS IN ORDER OF CONTRIBUTIONS
        DO N=1,NDIM-NSKP
          IORD(N) = 0
        ENDDO
C
        TSM = DCMPLX(0.0D0)
        DO M=1,NDIM-NSKP
          EBIG = 0.0D0
          DO N=1,NDIM-NSKP
            goto 74
            IF(IORD(N).EQ.1) GOTO 74
            IF(ABS(ET(N,ITT)).GT.EBIG) THEN
              EBIG = ABS(ET(N,ITT))
              K = N
            ENDIF
74          CONTINUE
          ENDDO
          k = m
          IF(ABS(ET(K,ITT)).LT.1.0D-9) GOTO 76
          if(LABKQN(K).ne.-2) goto 76
c          IF(ABS(ET(K,ITT)).LT.1.0D-9) GOTO 75
          IORD(K) = 1
          IF(LABKQN(K).LT.0) THEN
            ILQN =-LABKQN(K)-1
          ELSE
            ILQN = LABKQN(K)
          ENDIF
          CHL  = LLAB(ILQN)
          TSM = TSM + ET(K,ITT)
          WRITE(6,80) K,NQNLST(K),CHL,2*IABS(LABKQN(K))-1,LABMQN(K),
     &                EIGN(K+NSKP),RKNT(K),DREAL(EX(K,ITT)),
     &                DREAL(EY(K,ITT)),DREAL(EZ(K,ITT)),
     &                DREAL(ET(K,ITT)),DREAL(TSM)
          WRITE(7,80) K,NQNLST(K),CHL,2*IABS(LABKQN(K))-1,LABMQN(K),
     &                EIGN(K+NSKP),RKNT(K),DREAL(EX(K,ITT)),
     &                DREAL(EY(K,ITT)),DREAL(EZ(K,ITT)),
     &                DREAL(ET(K,ITT)),DREAL(TSM)
76        continue
        ENDDO
75      CONTINUE
        WRITE(6, *) REPEAT('=',NLNES)
        WRITE(7, *) REPEAT('=',NLNES)
C
      ENDDO
888   continue
C
      RETURN
      END

