      SUBROUTINE E119FILE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C     I USED THIS ROUTINE TO IMPORT A GRASP AMPLITUDE/UEHLING/DENSITY  C
C     FILE FOR E119+, THEN AT EACH RADIUS TO WRITE SOME EQUIVALENT     C
C     BERTHA OUTPUT.                                                   C
C**********************************************************************C
      INCLUDE 'parameters.h'
      INCLUDE 'scfoptions.h'
C
      CHARACTER*2  ELMT(120)
      CHARACTER*3  ATRM
      CHARACTER*5  NMDL
      CHARACTER*40 MOLCL,WFNFL,OUTFL
      CHARACTER*80 XOUT,TITLE,XAXIS,YAXIS,KEY(2)
C
      DIMENSION EXL(MBS),XYZ(3)
      DIMENSION RADII(19990),VPOT(19990),DTOT(19990)
      DIMENSION DGSP(19990)
      DIMENSION PGSP(19990,31),QGSP(19990,31)
      DIMENSION DLG(31),DSG(31)
      DIMENSION DLB(31),DSB(31)
      DIMENSION P(19990,31),Q(19990,31)
      DIMENSION IMAP(31),IDGN(31)
C
      COMPLEX*16 CONE
      COMPLEX*16 COEF(MDM,MDM)
C
      COMMON/BDIM/NDIM,NSKP,NOCC,NVRT
      COMMON/BNUC/ZNUC(MCT),ANUC(MCT),TFMI(MCT),AFMI(MCT),RNUC(MCT),
     &            FNUC(MCT,MFT),XNUC(MCT,MFT),NNUC(MCT),NMDL(MCT)
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/BUEH/RUEH(MCT,3),FUEH(MCT,MFT),XUEH(MCT,MFT),NUEH(MCT)
      COMMON/EIGC/COEF
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12,EULR
      COMMON/UEHQ/RAD(0:NRAD),VVAC(MCT,0:NRAD),RORI,RMID,RMAX,NLIN,NEXP
C
C     READ IN THE RADIAL INFORMATION
      OPEN(UNIT=8,FILE='gaussian_amplitudes_grasp.dat',STATUS='UNKNOWN')
        READ(8,*) 
        DO I=1,19990
          READ(8,*) J,RADII(I),V,DGSP(I),(PGSP(I,N),QGSP(I,N),N=1,31)
        ENDDO
      CLOSE(UNIT=8)
C
C     UNIT IMAGINARY NUMBER
      CONE = DCMPLX(0.0D0,1.0D0)
C
C**********************************************************************C
C     INITIALISE THE NECESSARY ARRAY ELEMENTS                          C
C**********************************************************************C
C
      DO I=1,19990
        DTOT(I) = 0.0D0
        DO N=1,31
          P(I,N) = 0.0D0
          Q(I,N) = 0.0D0
        ENDDO
      ENDDO
      DO I=1,19990
        DTOT(I) = 0.0D0
        DO N=1,31
          P(I,N) = 0.0D0
          Q(I,N) = 0.0D0
        ENDDO
      ENDDO
C
      IMAP( 1) =   1
      IMAP( 2) =   3
      IMAP( 3) =   5
      IMAP( 4) =   7
      IMAP( 5) =  11
      IMAP( 6) =  13
      IMAP( 7) =  15
      IMAP( 8) =  19
      IMAP( 9) =  23
      IMAP(10) =  29
      IMAP(11) =  31
      IMAP(12) =  33
      IMAP(13) =  37
      IMAP(14) =  41
      IMAP(15) =  47
      IMAP(16) =  53
      IMAP(17) =  61
      IMAP(18) =  63
      IMAP(19) =  65
      IMAP(20) =  69
      IMAP(21) =  73
      IMAP(22) =  83
      IMAP(23) =  89
      IMAP(24) =  80
      IMAP(25) =  81
      IMAP(26) =  97
      IMAP(27) = 101
      IMAP(28) = 105
      IMAP(29) = 111
      IMAP(30) = 113
      IMAP(31) = 115
C
      IDGN( 1) = 2
      IDGN( 2) = 2
      IDGN( 3) = 2
      IDGN( 4) = 4
      IDGN( 5) = 2
      IDGN( 6) = 2
      IDGN( 7) = 4
      IDGN( 8) = 4
      IDGN( 9) = 6
      IDGN(10) = 2
      IDGN(11) = 2
      IDGN(12) = 4
      IDGN(13) = 4
      IDGN(14) = 6
      IDGN(15) = 6
      IDGN(16) = 8
      IDGN(17) = 2
      IDGN(18) = 2
      IDGN(19) = 4
      IDGN(20) = 4
      IDGN(21) = 6
C      IDGN(22) = 6
      IDGN(22) = 2
C      IDGN(23) = 8
      IDGN(23) = 2
C      IDGN(24) = 2
      IDGN(24) = 4
C      IDGN(25) = 2
      IDGN(25) = 6
C      IDGN(26) = 4
      IDGN(26) = 8
      IDGN(27) = 4
      IDGN(28) = 6
      IDGN(29) = 2
      IDGN(30) = 2
      IDGN(31) = 4
      
      IEL = 0
      DO N=1,31
        IEL = IEL + IDGN(N)
      ENDDO
      WRITE(*,*) 'NUMBER OF ELECTRONS: ',IEL
C
C**********************************************************************C
C      CALCULATE CONTRIBUTIONS TO RADIAL FUNCTIONS                     C
C**********************************************************************C
C
C     LOOP OVER NUCLEAR CENTRES
      DO 1000 ICNTA=1,NCNT
C
C     LOOP OVER KQN VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       ORBITAL QUANTUM NUMBERS
        KQN = KAPA(KA,ICNTA)
        JQN = 2*IABS(KQN)-1
        LQN = LVAL(KQN)
C
C       BASIS EXPONENTS
        NBAS = NFNC(LQN,ICNTA)
        DO IBAS=1,NBAS
          EXL(IBAS) = BEXL(IBAS,LQN,ICNTA)
        ENDDO
C
C     LOOP OVER |MQN| VALUES
      DO 3000 MA=1,IABS(KQN)
C
C     EXPANSION COEFFICIENT MATRIX STARTING ADDRESSES
      NA1 = LRGE(ICNTA,KA,2*MA-1)
C
C     PRE-FACTORS FOR NORMALISATION CONSTANTS
      RL = DFLOAT(LQN)
      G1 = TWLG-GAMLOG(2*LQN+3)
      G2 = TWLG-GAMLOG(2*LQN+5)
      R1 = RL+1.5D0
      R2 = RL+0.5D0
C
C     LOOP OVER BASIS FUNCTIONS
      DO IBAS=1,NBAS
C
C       NORMALISATION CONSTANTS
        ELOG = DLOG(2.0D0*EXL(IBAS))
        RNL  = DEXP(0.5D0*(G1+R1*ELOG))
        RNS  = DEXP(0.5D0*(G2+R2*ELOG))
C
C       LOOP OVER ALL RADII
        DO I=1,19990
C
          R = RADII(I)
C
C         LARGE AND SMALL RADIAL AMPLITUDES
          DR = R**(LQN)
          DE = DEXP(-EXL(IBAS)*R*R)
          FL = RNL*DR*R*DE
          FS = RNS*DR*(KQN+LQN+1.0D0-2.0D0*EXL(IBAS)*R*R)*DE
C
C         ADDRESS OFFSETS FOR SMALL-COMPONENTS
          KBAS = IBAS+NSKP
C
C         LOOP OVER OCCUPIED ORBITALS
          DO N=1,31
C
C           MULTIPLY BY CORRESPONDING EXPANSION COEFFICIENT ELEMENTS
            P(I,N) = P(I,N) + COEF(NA1+IBAS,IMAP(N)+NSKP)*FL
            Q(I,N) = Q(I,N) + COEF(NA1+KBAS,IMAP(N)+NSKP)*FS
C
          ENDDO
C          
        ENDDO
C
      ENDDO
C
C     END LOOP OVER BASIS QUANTUM NUMBERS
3000  CONTINUE
2000  CONTINUE
1000  CONTINUE
C
C     OPEN A NEW FILE AND LOOP OVER RADII
      OPEN(UNIT=9,FILE='gaussian_amplitudes_bertha.dat',
     &                                               STATUS='UNKNOWN')
      PTOT = 0.0D0
      QTOT = 0.0D0
      GTOT = 0.0D0
      DO I=1,19990
        IF(I.EQ.1) THEN
          HS = 0.0D0
        ELSE
          HS = RADII(I)-RADII(I-1)
        ENDIF
        VU = 0.0D0
        DO N=1,26
          VU = VU + FUEH(1,N)*DEXP(-XUEH(1,N)*RADII(I)*RADII(I))
        ENDDO
        DO N=1,31
          DTOT(I) = DTOT(I) + IDGN(N)*(P(I,N)*P(I,N)+Q(I,N)*Q(I,N))
          PTOT = PTOT + IDGN(N)*P(I,N)*P(I,N)*HS
          QTOT = QTOT + IDGN(N)*Q(I,N)*Q(I,N)*HS
        ENDDO
        GTOT = GTOT + DGSP(I)*HS
C
C       WRITE OUT RESULTS
        WRITE(9,*) J,RADII(I),VU*119.0D0,DTOT(I),(P(I,N),Q(I,N),N=1,31)
        IF(MOD(I,100).EQ.0) THEN
          WRITE(*,*) RADII(I),P(I,1),PGSP(I,1)
        ENDIF
C
C     END LOOP OVER RADII AND CLOSE FILE
      ENDDO
      CLOSE(UNIT=9)
C
C     BERTHA ENCLOSED CHARGE IN DENSITIES
      DO N=1,31
        DLB(N) = 0.0D0
        DSB(N) = 0.0D0
        DO I=1,19990
          IF(I.EQ.1) THEN
            HS = 0.0D0
          ELSE
            HS = RADII(I)-RADII(I-1)
          ENDIF
          DLB(N) = DLB(N) + P(I,N)*P(I,N)*HS
          DSB(N) = DSB(N) + Q(I,N)*Q(I,N)*HS
        ENDDO
      ENDDO
C
C     BERTHA NET CHARGES
      BLT = 0.0D0
      BST = 0.0D0
      DO N=1,31
        BLT = BLT + IDGN(N)*DLB(N)
        BST = BST + IDGN(N)*DSB(N)
        WRITE(*,*) 'BERTHA:',N,DLB(N),DSB(N),DLB(N)+DSB(N)
        WRITE(*, *)
      ENDDO
C
C     GRASP ENCLOSED CHARGE IN DENSITIES
      DO N=1,31
        DLG(N) = 0.0D0
        DSG(N) = 0.0D0
        DO I=1,19990
          IF(I.EQ.1) THEN
            HS = 0.0D0
          ELSE
            HS = RADII(I)-RADII(I-1)
          ENDIF
          DLG(N) = DLG(N) + PGSP(I,N)*PGSP(I,N)*HS
          DSG(N) = DSG(N) + QGSP(I,N)*QGSP(I,N)*HS
        ENDDO
      ENDDO
C
C     GRASP NET CHARGES
      GLT = 0.0D0
      GST = 0.0D0
      DO N=1,31
        GLT = GLT + IDGN(N)*DLG(N)
        GST = GST + IDGN(N)*DSG(N)
        WRITE(*,*) 'GRASP: ',N,DLG(N),DSG(N),DLG(N)+DSG(N)
        WRITE(*, *)
      ENDDO
C
C     READ IN THE RADIAL INFORMATION
c      OPEN(UNIT=8,FILE='gaussian_amplitudes_grasp.dat',
c     &                                                STATUS='UNKNOWN')
c        WRITE(8,*) '#structure: i r(i) v(i) dens(i) p1(i) q1(i) p2(i)..'
c        DO I=1,19990
c          DGSP(I) = 0.0D0
c          DO N=1,31
c            D = PGSP(I,N)*PGSP(I,N) + QGSP(I,N)*QGSP(I,N)
c            DGSP(I) = DGSP(I) + IDGN(N)*D
c          ENDDO
c          WRITE(8,*) I,RADII(I),V,DGSP(I),(PGSP(I,N),QGSP(I,N),N=1,31)
c        ENDDO
c      CLOSE(UNIT=8)
C
C     SUMMARISE RESULTS
      WRITE(*,*) 'TOTALS'
      WRITE(*,*) 'BERTHA:',BLT,BST,BLT+BST
      WRITE(*,*) 'GRASP: ',GLT,GST,GLT+GST,GTOT
C
      RETURN
      END

