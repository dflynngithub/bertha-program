      SUBROUTINE ORTHGNL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C     OOOOOO  RRRRRRR TTTTTTTT HH    HH  GGGGGG  NN    NN LL           C
C    OO    OO RR    RR   TT    HH    HH GG    GG NNN   NN LL           C
C    OO    OO RR    RR   TT    HH    HH GG       NNNN  NN LL           C
C    OO    OO RR    RR   TT    HHHHHHHH GG       NN NN NN LL           C
C    OO    OO RRRRRRR    TT    HH    HH GG   GGG NN  NNNN LL           C
C    OO    OO RR    RR   TT    HH    HH GG    GG NN   NNN LL           C
C     OOOOOO  RR    RR   TT    HH    HH  GGGGGG  NN    NN LLLLLLLL     C
C                                                                      C
C -------------------------------------------------------------------- C
C  ORTHGNL PERFORMS A MOLECULAR ORTHOGONALITY ANALYSIS.                C
C**********************************************************************C
      PARAMETER(MDM=1200,MBS=26,MCT=12,MKP=9)
C
      CHARACTER*40 MOLCL,WFNFL,OUTFL
C
      COMPLEX*16 E1(5)
      COMPLEX*16 V1(MDM,MDM)
      COMPLEX*16 OLAPLL(MDM,MDM),OLAPSS(MDM,MDM),EMPTY(MDM,MDM)
C
      DIMENSION EPOS(0:21),ENEG(0:21),ETOT(0:21)
C
      COMMON/FLNM/MOLCL,WFNFL,OUTFL
      COMMON/SPEC/EXPSET(MBS,MKP,MCT),COORD(3,MCT),ZNUC(MCT),AMASS(MCT),
     &            CNUC(MCT),PNUC,LARGE(MCT,MKP,MKP+1),NFUNCT(MKP,MCT),
     &            KVALS(MKP,MCT),IZNUC(MCT),IQNUC(MCT),LMAX(MCT),
     &            NKAP(MCT),NCNT,NDIM,NSHIFT,NOCC,NVRT
C
C     INITIALISE DUMMY ARRAY
      DO I=1,NDIM-NSHIFT
        DO J=1,NDIM-NSHIFT
          EMPTY(I,J) = DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO
C
C     GENERATE DIRECT OVERLAP MATRIX (ZEROTH MOMENT)
      CALL VMOMNT0(OLAPLL,1,0,1,2)
      CALL VMOMNT0(OLAPSS,4,0,1,2)
C
C**********************************************************************C
C     MOLECULAR EXPECTATION VALUE                                      C
C**********************************************************************C
C
      WRITE(6, *) 'Direct overlap analysis:'
      WRITE(7, *) 'Direct overlap analysis:'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     CALCULATE MOLECULAR EXPECTATION VALUE OVER THE DENSITY MATRIX
      CALL PROPRTY(E1,OLAPLL,OLAPSS,EMPTY,EMPTY)
C
C     WRITE OCCUPATION NUMBER EXPECTATION VALUE
20    FORMAT(1X,A,F15.10)
      WRITE(6,20) 'Occupation number = ',DREAL(E1(5))
      WRITE(7,20) 'Occupation number = ',DREAL(E1(5))
C
C     END SECTION
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C**********************************************************************C
C     ORBITAL EXPECTATION VALUES                                       C
C**********************************************************************C
C
C     MAXIMUM ORDER IN PERTURBATIVE EXPANSION
      NORD = 1
C
30    FORMAT(1X,A,12X,A,I1,A,13X,A,I1,A,13X,A,I1,A)
31    FORMAT(1X,I3,2X,F21.14,2X,F21.14,2X,F21.14)
      WRITE(6, *) 'Expectation value for each orbital IOCC:'
      WRITE(7, *) 'Expectation value for each orbital IOCC:'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
      WRITE(6,30) 'Orb.','E^(',NORD,')[pos]','E^{',NORD,'}[neg]',
     &                                       'E^{',NORD,'}[tot]'
      WRITE(7,30) 'Orb.','E^(',NORD,')[pos]','E^{',NORD,'}[neg]',
     &                                       'E^{',NORD,'}[tot]'
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     GENERATE ARRAY OF FIRST ORDER MATRIX ELEMENTS
      CALL RSPT1(V1,OLAPLL,OLAPSS,EMPTY,EMPTY,NORD)
C
C     PRINT MATRIX ELEMENTS TO FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_ELCMNPL_V1.dat',STATUS='UNKNOWN')
        DO IOCC=1,NDIM
          WRITE(8, *) (V1(IOCC,JOCC),JOCC=1,NDIM)
        ENDDO
      CLOSE(UNIT=8)
C
C     MAPPING BETWEEN NORD AND TOTAL ENERGY AT THAT ORDER
      IF(NORD.EQ.0) THEN
        NAD = 0
      ELSEIF(NORD.EQ.1) THEN
        NAD = 1
      ELSEIF(NORD.EQ.2) THEN
        NAD = 2
      ELSEIF(NORD.EQ.3) THEN
        NAD = 5
      ELSEIF(NORD.EQ.4) THEN
        NAD = 10
      ELSEIF(NORD.EQ.5) THEN
        NAD = 20
      ENDIF
C
C     OPEN EXTERNAL FILE
      OPEN(UNIT=8,FILE=TRIM(OUTFL)//'_ELCMNPL_EN.dat',STATUS='UNKNOWN')
C
C       LOOP OVER ALL OCCUPIED ORBITALS
        DO IOCC=1,NOCC
C
C         DIAGRAMATIC PERTURBATION THEORY ON V1
          CALL DIAGRMTC(V1,EPOS,ENEG,ETOT,NORD,NSHIFT+IOCC,0)
C
C         WRITE RESULTS TO FILE
          WRITE(6,31) IOCC,EPOS(NAD),ENEG(NAD),ETOT(NAD)
          WRITE(7,31) IOCC,EPOS(NAD),ENEG(NAD),ETOT(NAD)
          WRITE(8,31) IOCC,EPOS(NAD),ENEG(NAD),ETOT(NAD)
C
C       END LOOP OVER OCCUPIED ORBITALS
        ENDDO
C
C     CLOSE EXTERNAL FILE
      CLOSE(UNIT=8)
C
C     END SECTION
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
C**********************************************************************C
C     SET OF GOLDSTONE DIAGRAMS FOR ONE ORBITAL                        C
C**********************************************************************C
C
C     SPECIFY THE ORBITAL
      IOCC = 1
C
C     MAXIMUM PERTURBATIVE ORDER
      NORD = 4
C
40    FORMAT(1X,A,I3)
      WRITE(6,40) 'Goldstone diagram values for IOCC = ',IOCC
      WRITE(7,40) 'Goldstone diagram values for IOCC = ',IOCC
      WRITE(6, *) REPEAT('-',72)
      WRITE(7, *) REPEAT('-',72)
C
C     DIAGRAMATIC PERTURBATION THEORY ON IOCC DUE TO V1
      CALL DIAGRMTC(V1,EPOS,ENEG,ETOT,NORD,NSHIFT+IOCC,1)
C
      WRITE(6, *) REPEAT('=',72)
      WRITE(7, *) REPEAT('=',72)
C
      RETURN
      END
