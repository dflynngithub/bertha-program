c     these versions also produce the phase-adapted ~ 
c                                                   E[T,T';A,B,C]

      SUBROUTINE E0LLGN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         EEEEEEEE 000000  LL       LL       GGGGGG  NN    NN          C
C         EE      00   000 LL       LL      GG    GG NNN   NN          C
C         EE      00  0000 LL       LL      GG       NNNN  NN          C
C         EEEEEE  00 00 00 LL       LL      GG       NN NN NN          C
C         EE      0000  00 LL       LL      GG   GGG NN  NNNN          C
C         EE      000   00 LL       LL      GG    GG NN   NNN          C
C         EEEEEEEE 000000  LLLLLLLL LLLLLLLL GGGGGG  NN    NN          C
C                                                                      C
C -------------------------------------------------------------------- C
C  E0LLGN GENERATES A FULL SET OF E0LL COEFFICIENTS AND SAVES TO A     C
C  COMMON ARRAY E0LL, INCLUDING AN ADDRESS INDEX.                      C
C**********************************************************************C
      INCLUDE 'param.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/E0LL/E0LLFL(MFL,8),IAD0LL(MCT,MCT,MKP,MKP,MKP,MKP)
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     INITIALISE COUNT INDEX
      ICOUNT = 0
C
C     LOOP OVER CENTRE A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,ICNTA)
        IF(KQN(1).LT.0) THEN
          LQN(1) =-KQN(1)-1
        ELSE
          LQN(1) = KQN(1)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),ICNTA)
        DO IBAS=1,NBAS(1)
          EXL(IBAS,1) = BEXL(IBAS,LQN(1),ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KAPA(KB,ICNTB)
        IF(KQN(2).LT.0) THEN
          LQN(2) =-KQN(2)-1
        ELSE
          LQN(2) = KQN(2)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
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
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS (MA,MB).  C
C     FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE ORDERED:             C
C -------------------------------------------------------------------- C
C  11: = (-|MQN(A)|,-|MQN(B)|) -> 1  12: = (-|MQN(A)|,+|MQN(B)|) -> 2  C
C  21: = (+|MQN(A)|,-|MQN(B)|) -> 3  22: = (+|MQN(A)|,+|MQN(B)|) -> 4  C
C**********************************************************************C
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NBAS(1)*NBAS(2)
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMAB  = LQN(1)+LQN(2)
      NTUVLL = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
C
C     INDEX TO TRACK START OF THIS BLOCK OF COEFFICIENTS
      IAD0LL(ICNTA,ICNTB,KA,KB,MA,MB) = ICOUNT
C
C     GENERATE ELL0(AB) COEFFICIENTS
      CALL EQLLMK(E11,E21,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,0)
C
C     WRITE ELL0(AB) TO THE MAIN ARRAY
      DO ITUV=1,NTUVLL
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          E0LLFL(IAD+M,1) = DREAL(E11(M,ITUV))
          E0LLFL(IAD+M,2) = DIMAG(E11(M,ITUV))
          E0LLFL(IAD+M,3) = DREAL(E21(M,ITUV))
          E0LLFL(IAD+M,4) = DIMAG(E21(M,ITUV))
        ENDDO
      ENDDO
C
C     GENERATE ELL0(CD) COEFFICIENTS
      CALL EQLLMK(E11,E21,EXL,XYZ,KQN,MQN,NBAS,-1,1,2,0)
C
C     WRITE ELL0(CD) TO THE MAIN ARRAY
      DO ITUV=1,NTUVLL
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          E0LLFL(IAD+M,5) = DREAL(E11(M,ITUV))
          E0LLFL(IAD+M,6) = DIMAG(E11(M,ITUV))
          E0LLFL(IAD+M,7) = DREAL(E21(M,ITUV))
          E0LLFL(IAD+M,8) = DIMAG(E21(M,ITUV))
        ENDDO
      ENDDO
C
C     INCREASE COUNT INDEX
      ICOUNT = ICOUNT + NTUVLL*MAXAB
C
C     END LOOPS OVER FOCK BLOCK
2000  CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE E0SSGN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         EEEEEEEE 000000   SSSSSS   SSSSSS   GGGGGG  NN    NN         C
C         EE      00   000 SS    SS SS    SS GG    GG NNN   NN         C
C         EE      00  0000 SS       SS       GG       NNNN  NN         C
C         EEEEEE  00 00 00  SSSSSS   SSSSSS  GG       NN NN NN         C
C         EE      0000  00       SS       SS GG   GGG NN  NNNN         C
C         EE      000   00 SS    SS SS    SS GG    GG NN   NNN         C
C         EEEEEEEE 000000   SSSSSS   SSSSSS   GGGGGG  NN    NN         C
C                                                                      C
C -------------------------------------------------------------------- C
C  E0SSGN GENERATES A FULL SET OF E0SS COEFFICIENTS AND SAVES TO A     C
C  COMMON ARRAY E0SS, INCLUDING AN ADDRESS INDEX.                      C
C**********************************************************************C
      INCLUDE 'param.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
C
      COMPLEX*16 E11(MB2,MEQ),E21(MB2,MEQ)
C
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/E0SS/E0SSFL(MFL,8),IAD0SS(MCT,MCT,MKP,MKP,MKP,MKP)
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     INITIALISE COUNT INDEX
      ICOUNT = 0
C
C     LOOP OVER CENTRE A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,ICNTA)
        IF(KQN(1).LT.0) THEN
          LQN(1) =-KQN(1)-1
        ELSE
          LQN(1) = KQN(1)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),ICNTA)
        DO IBAS=1,NBAS(1)
          EXL(IBAS,1) = BEXL(IBAS,LQN(1),ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KAPA(KB,ICNTB)
        IF(KQN(2).LT.0) THEN
          LQN(2) =-KQN(2)-1
        ELSE
          LQN(2) = KQN(2)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
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
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS (MA,MB).  C
C     FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE ORDERED:             C
C -------------------------------------------------------------------- C
C  11: = (-|MQN(A)|,-|MQN(B)|) -> 1  12: = (-|MQN(A)|,+|MQN(B)|) -> 2  C
C  21: = (+|MQN(A)|,-|MQN(B)|) -> 3  22: = (+|MQN(A)|,+|MQN(B)|) -> 4  C
C**********************************************************************C
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NBAS(1)*NBAS(2)
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMAB  = LQN(1)+LQN(2)+2
      NTUVSS = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
C
C     INDEX TO TRACK START OF THIS BLOCK OF COEFFICIENTS
      IAD0SS(ICNTA,ICNTB,KA,KB,MA,MB) = ICOUNT
C
C     GENERATE ESS0(AB) COEFFICIENTS
      CALL EQSSMK(E11,E21,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,0)
C
C     WRITE ESS0(AB) TO THE MAIN ARRAY
      DO ITUV=1,NTUVSS
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          E0SSFL(IAD+M,1) = DREAL(E11(M,ITUV))
          E0SSFL(IAD+M,2) = DIMAG(E11(M,ITUV))
          E0SSFL(IAD+M,3) = DREAL(E21(M,ITUV))
          E0SSFL(IAD+M,4) = DIMAG(E21(M,ITUV))
        ENDDO
      ENDDO
C
C     GENERATE ESS0(CD) COEFFICIENTS
      CALL EQSSMK(E11,E21,EXL,XYZ,KQN,MQN,NBAS,-1,1,2,0)
C
C     WRITE ESS0(CD) TO THE MAIN ARRAY
      DO ITUV=1,NTUVSS
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          E0SSFL(IAD+M,5) = DREAL(E11(M,ITUV))
          E0SSFL(IAD+M,6) = DIMAG(E11(M,ITUV))
          E0SSFL(IAD+M,7) = DREAL(E21(M,ITUV))
          E0SSFL(IAD+M,8) = DIMAG(E21(M,ITUV))
        ENDDO
      ENDDO
C
C     INCREASE COUNT INDEX
      ICOUNT = ICOUNT + NTUVSS*MAXAB
C
C     END LOOPS OVER FOCK BLOCK
2000  CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE EILSGN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           EEEEEEEE IIII LL       SSSSSS   GGGGGG  NN    NN           C
C           EE        II  LL      SS    SS GG    GG NNN   NN           C
C           EE        II  LL      SS       GG       NNNN  NN           C
C           EEEEEE    II  LL       SSSSSS  GG       NN NN NN           C
C           EE        II  LL            SS GG   GGG NN  NNNN           C
C           EE        II  LL      SS    SS GG    GG NN   NNN           C
C           EEEEEEEE IIII LLLLLLLL SSSSSS   GGGGGG  NN    NN           C
C                                                                      C
C -------------------------------------------------------------------- C
C  EILSGN GENERATES A FULL SET OF EILS COEFFICIENTS AND SAVES TO A     C
C  COMMON ARRAY EILS, INCLUDING AN ADDRESS INDEX.                      C
C**********************************************************************C
      INCLUDE 'param.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
C
      COMPLEX*16 E11X(MB2,MEQ),E11Y(MB2,MEQ),E11Z(MB2,MEQ),
     &           E21X(MB2,MEQ),E21Y(MB2,MEQ),E21Z(MB2,MEQ)
C
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EILS/EILSFL(MFL,24),IADILS(MCT,MCT,MKP,MKP,MKP,MKP)
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     INITIALISE COUNT INDEX
      ICOUNT = 0
C
C     LOOP OVER CENTRE A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,ICNTA)
        IF(KQN(1).LT.0) THEN
          LQN(1) =-KQN(1)-1
        ELSE
          LQN(1) = KQN(1)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),ICNTA)
        DO IBAS=1,NBAS(1)
          EXL(IBAS,1) = BEXL(IBAS,LQN(1),ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KAPA(KB,ICNTB)
        IF(KQN(2).LT.0) THEN
          LQN(2) =-KQN(2)-1
        ELSE
          LQN(2) = KQN(2)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
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
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS (MA,MB).  C
C     FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE ORDERED:             C
C -------------------------------------------------------------------- C
C  11: = (-|MQN(A)|,-|MQN(B)|) -> 1  12: = (-|MQN(A)|,+|MQN(B)|) -> 2  C
C  21: = (+|MQN(A)|,-|MQN(B)|) -> 3  22: = (+|MQN(A)|,+|MQN(B)|) -> 4  C
C**********************************************************************C
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NBAS(1)*NBAS(2)
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMAB  = LQN(1)+LQN(2)+1
      NTUVLS = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
C
C     INDEX TO TRACK START OF THIS BLOCK OF COEFFICIENTS
      IADILS(ICNTA,ICNTB,KA,KB,MA,MB) = ICOUNT
C
C     GENERATE EILS(AB) COEFFICIENTS
      CALL EQLSMK(E11X,E21X,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,1)
      CALL EQLSMK(E11Y,E21Y,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,2)
      CALL EQLSMK(E11Z,E21Z,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,3)
C
C     WRITE EILS(AB) TO THE MAIN ARRAY
      DO ITUV=1,NTUVLS
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          EILSFL(IAD+M, 1) = DREAL(E11X(M,ITUV))
          EILSFL(IAD+M, 2) = DIMAG(E11X(M,ITUV))
          EILSFL(IAD+M, 3) = DREAL(E21X(M,ITUV))
          EILSFL(IAD+M, 4) = DIMAG(E21X(M,ITUV))
          EILSFL(IAD+M, 5) = DREAL(E11Y(M,ITUV))
          EILSFL(IAD+M, 6) = DIMAG(E11Y(M,ITUV))
          EILSFL(IAD+M, 7) = DREAL(E21Y(M,ITUV))
          EILSFL(IAD+M, 8) = DIMAG(E21Y(M,ITUV))
          EILSFL(IAD+M, 9) = DREAL(E11Z(M,ITUV))
          EILSFL(IAD+M,10) = DIMAG(E11Z(M,ITUV))
          EILSFL(IAD+M,11) = DREAL(E21Z(M,ITUV))
          EILSFL(IAD+M,12) = DIMAG(E21Z(M,ITUV))
        ENDDO
      ENDDO
C
C     GENERATE EILS(CD) COEFFICIENTS
      CALL EQLSMK(E11X,E21X,EXL,XYZ,KQN,MQN,NBAS,-1,1,2,1)
      CALL EQLSMK(E11Y,E21Y,EXL,XYZ,KQN,MQN,NBAS,-1,1,2,2)
      CALL EQLSMK(E11Z,E21Z,EXL,XYZ,KQN,MQN,NBAS,-1,1,2,3)
C
C     WRITE EILS(CD) TO THE MAIN ARRAY
      DO ITUV=1,NTUVLS
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          EILSFL(IAD+M,13) = DREAL(E11X(M,ITUV))
          EILSFL(IAD+M,14) = DIMAG(E11X(M,ITUV))
          EILSFL(IAD+M,15) = DREAL(E21X(M,ITUV))
          EILSFL(IAD+M,16) = DIMAG(E21X(M,ITUV))
          EILSFL(IAD+M,17) = DREAL(E11Y(M,ITUV))
          EILSFL(IAD+M,18) = DIMAG(E11Y(M,ITUV))
          EILSFL(IAD+M,19) = DREAL(E21Y(M,ITUV))
          EILSFL(IAD+M,20) = DIMAG(E21Y(M,ITUV))
          EILSFL(IAD+M,21) = DREAL(E11Z(M,ITUV))
          EILSFL(IAD+M,22) = DIMAG(E11Z(M,ITUV))
          EILSFL(IAD+M,23) = DREAL(E21Z(M,ITUV))
          EILSFL(IAD+M,24) = DIMAG(E21Z(M,ITUV))
        ENDDO
      ENDDO
C
C     INCREASE COUNT INDEX
      ICOUNT = ICOUNT + NTUVLS*MAXAB
C
C     END LOOPS OVER FOCK BLOCK
2000  CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE EISLGN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C           EEEEEEEE IIII  SSSSSS  LL       GGGGGG  NN    NN           C
C           EE        II  SS    SS LL      GG    GG NNN   NN           C
C           EE        II  SS       LL      GG       NNNN  NN           C
C           EEEEEE    II   SSSSSS  LL      GG       NN NN NN           C
C           EE        II        SS LL      GG   GGG NN  NNNN           C
C           EE        II  SS    SS LL      GG    GG NN   NNN           C
C           EEEEEEEE IIII  SSSSSS  LLLLLLLL GGGGGG  NN    NN           C
C                                                                      C
C -------------------------------------------------------------------- C
C  EISLGN GENERATES A FULL SET OF EISL COEFFICIENTS AND SAVES TO A     C
C  COMMON ARRAY EISL, INCLUDING AN ADDRESS INDEX.                      C
C**********************************************************************C
      INCLUDE 'param.h'
C
      DIMENSION EXL(MBS,4),XYZ(3,4),KQN(4),LQN(4),MQN(4),NBAS(4)
C
      COMPLEX*16 E11X(MB2,MEQ),E11Y(MB2,MEQ),E11Z(MB2,MEQ),
     &           E21X(MB2,MEQ),E21Y(MB2,MEQ),E21Z(MB2,MEQ)
C
      COMMON/BSET/BEXL(MBS,0:MEL,MCT),BXYZ(3,MCT),LRGE(MCT,MKP,MKP+1),
     &            KAPA(MKP,MCT),NFNC(0:MEL,MCT),NKAP(MCT),IQNC(MCT),NCNT
      COMMON/EISL/EISLFL(MFL,24),IADISL(MCT,MCT,MKP,MKP,MKP,MKP)
C
C**********************************************************************C
C     FIRST LAYER OF LOOPS, OVER CENTRES A AND B (USE INDEX 2000)      C
C**********************************************************************C
C
C     INITIALISE COUNT INDEX
      ICOUNT = 0
C
C     LOOP OVER CENTRE A
      DO 2000 ICNTA=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE A
        XYZ(1,1) = BXYZ(1,ICNTA)
        XYZ(2,1) = BXYZ(2,ICNTA)
        XYZ(3,1) = BXYZ(3,ICNTA)
C
C     LOOP OVER CENTRE B
      DO 2000 ICNTB=1,NCNT
C
C       CARTESIAN COORDINATES OF CENTRE B
        XYZ(1,2) = BXYZ(1,ICNTB)
        XYZ(2,2) = BXYZ(2,ICNTB)
        XYZ(3,2) = BXYZ(3,ICNTB)
C
C     LOOP OVER KQN(A) VALUES
      DO 2000 KA=1,NKAP(ICNTA)
C
C       QUANTUM NUMBERS FOR BLOCK A
        KQN(1) = KAPA(KA,ICNTA)
        IF(KQN(1).LT.0) THEN
          LQN(1) =-KQN(1)-1
        ELSE
          LQN(1) = KQN(1)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK A
        NBAS(1) = NFNC(LQN(1),ICNTA)
        DO IBAS=1,NBAS(1)
          EXL(IBAS,1) = BEXL(IBAS,LQN(1),ICNTA)
        ENDDO
C
C     LOOP OVER KQN(B) VALUES
      DO 2000 KB=1,NKAP(ICNTB)
C
C       QUANTUM NUMBERS FOR BLOCK B
        KQN(2) = KAPA(KB,ICNTB)
        IF(KQN(2).LT.0) THEN
          LQN(2) =-KQN(2)-1
        ELSE
          LQN(2) = KQN(2)
        ENDIF
C
C       BASIS EXPONENTS FOR BLOCK B
        NBAS(2) = NFNC(LQN(2),ICNTB)
        DO JBAS=1,NBAS(2)
          EXL(JBAS,2) = BEXL(JBAS,LQN(2),ICNTB)
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
C     AT THIS POINT, WE ARE WITHIN A BLOCK OF 4 COMBINATIONS (MA,MB).  C
C     FOR GIVEN (|MA|,|MB|), THE COMBINATIONS ARE ORDERED:             C
C -------------------------------------------------------------------- C
C  11: = (-|MQN(A)|,-|MQN(B)|) -> 1  12: = (-|MQN(A)|,+|MQN(B)|) -> 2  C
C  21: = (+|MQN(A)|,-|MQN(B)|) -> 3  22: = (+|MQN(A)|,+|MQN(B)|) -> 4  C
C**********************************************************************C
C
C     NUMBER OF BASIS FUNCTION OVERLAPS
      MAXAB = NBAS(1)*NBAS(2)
C
C     CALCULATE LAM VALUES FOR THIS OVERLAP CHOICE
      LAMAB  = LQN(1)+LQN(2)+1
      NTUVSL = (LAMAB+1)*(LAMAB+2)*(LAMAB+3)/6
C
C     INDEX TO TRACK START OF THIS BLOCK OF COEFFICIENTS
      IADISL(ICNTA,ICNTB,KA,KB,MA,MB) = ICOUNT
C
C     GENERATE EISL(AB) COEFFICIENTS
      CALL EQSLMK(E11X,E21X,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,1)
      CALL EQSLMK(E11Y,E21Y,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,2)
      CALL EQSLMK(E11Z,E21Z,EXL,XYZ,KQN,MQN,NBAS,+1,1,2,3)
C
C     WRITE EISL(AB) TO THE MAIN ARRAY
      DO ITUV=1,NTUVSL
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          EISLFL(IAD+M, 1) = DREAL(E11X(M,ITUV))
          EISLFL(IAD+M, 2) = DIMAG(E11X(M,ITUV))
          EISLFL(IAD+M, 3) = DREAL(E21X(M,ITUV))
          EISLFL(IAD+M, 4) = DIMAG(E21X(M,ITUV))
          EISLFL(IAD+M, 5) = DREAL(E11Y(M,ITUV))
          EISLFL(IAD+M, 6) = DIMAG(E11Y(M,ITUV))
          EISLFL(IAD+M, 7) = DREAL(E21Y(M,ITUV))
          EISLFL(IAD+M, 8) = DIMAG(E21Y(M,ITUV))
          EISLFL(IAD+M, 9) = DREAL(E11Z(M,ITUV))
          EISLFL(IAD+M,10) = DIMAG(E11Z(M,ITUV))
          EISLFL(IAD+M,11) = DREAL(E21Z(M,ITUV))
          EISLFL(IAD+M,12) = DIMAG(E21Z(M,ITUV))
        ENDDO
      ENDDO
C
C     GENERATE EISL(CD) COEFFICIENTS
      CALL EQSLMK(E11X,E21X,EXL,XYZ,KQN,MQN,NBAS,-1,1,2,1)
      CALL EQSLMK(E11Y,E21Y,EXL,XYZ,KQN,MQN,NBAS,-1,1,2,2)
      CALL EQSLMK(E11Z,E21Z,EXL,XYZ,KQN,MQN,NBAS,-1,1,2,3)
C
C     WRITE EISL(CD) TO THE MAIN ARRAY
      DO ITUV=1,NTUVSL
        IAD = ICOUNT + (ITUV-1)*MAXAB
        DO M=1,MAXAB
          EISLFL(IAD+M,13) = DREAL(E11X(M,ITUV))
          EISLFL(IAD+M,14) = DIMAG(E11X(M,ITUV))
          EISLFL(IAD+M,15) = DREAL(E21X(M,ITUV))
          EISLFL(IAD+M,16) = DIMAG(E21X(M,ITUV))
          EISLFL(IAD+M,17) = DREAL(E11Y(M,ITUV))
          EISLFL(IAD+M,18) = DIMAG(E11Y(M,ITUV))
          EISLFL(IAD+M,19) = DREAL(E21Y(M,ITUV))
          EISLFL(IAD+M,20) = DIMAG(E21Y(M,ITUV))
          EISLFL(IAD+M,21) = DREAL(E11Z(M,ITUV))
          EISLFL(IAD+M,22) = DIMAG(E11Z(M,ITUV))
          EISLFL(IAD+M,23) = DREAL(E21Z(M,ITUV))
          EISLFL(IAD+M,24) = DIMAG(E21Z(M,ITUV))
        ENDDO
      ENDDO
C
C     INCREASE COUNT INDEX
      ICOUNT = ICOUNT + NTUVSL*MAXAB
C
C     END LOOPS OVER FOCK BLOCK
2000  CONTINUE
C
      RETURN
      END
