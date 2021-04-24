      SUBROUTINE FUNFM(FM,T,LENGTH,MMAX,ITYPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C           FFFFFFFF UU     UU NN    NN FFFFFFFF MM      MM           C
C           FF       UU     UU NNN   NN FF       MMM    MMM           C
C           FF       UU     UU NNNN  NN FF       MMMM  MMMM           C
C           FFFFFF   UU     UU NN NN NN FFFFFF   MM MMMM MM           C
C           FF       UU     UU NN  NNNN FF       MM  MM  MM           C
C           FF       UU     UU NN   NNN FF       MM      MM           C
C           FF        UUUUUUU  NN    NN FF       MM      MM           C
C                                                                     C
C ------------------------------------------------------------------- C
C     FUNFM EVALUATES INTEGRAL [INT_{0}^{1} U^{2M} EXP(-TU^{2}) DU]   C
C     FOR VARIABLE T > 0 FOR ALL ORDERS 0 < M < MMAX.                 C
C ------------------------------------------------------------------- C
C     ITYPE = 0  SPECIAL CASE T = 0.0D0                               C
C     ITYPE = 1  POWER SERIES AND REVERSE RECURRENCE                  C
C                (ONLY MSERIES TERMS WILL BE USED, SO USE MUST SUPPLY C
C                A VALUE APPROPRIATE TO THE MAX VALUE OF T IN BATCH). C
C     ITYPE = 2  ASYMPTOTIC EXPANSION AND FORWARD RECURRENCE.         C        
C     ITYPE = 3  ASYMPTOTIC EXPANSION AND FORWARD RECURRENCE.         C
C                ALL TERMS DEPENDING ON EXP(-T) ARE OMITTED TO AVOID  C
C                NUMERICAL UNDERFLOW PROBLEMS. MSERIES NOT REQUIRED.  C
C*********************************************************************C
      PARAMETER(MAXB=26,MAXB2=MAXB*MAXB,MAXL=26,MSERIES=60)

      DIMENSION FM(MAXB2,MAXL),T(MAXB2),TMMAX(MAXB2),TT2(MAXB2),
     &          TEXP(MAXB2),TROOT(MAXB2)
     
      DATA PIROOT,A0,B0/8.862269254527580D-1,4.994501191201870D-1,
     &                  4.551838436668326D-1/
C
C     ITYPE = 0: SPECIAL CASE FOR T = 0.0D0
      IF(ITYPE.EQ.0) THEN
        DO K=1,MMAX+1
          MVAL  = K-1
          VALUE = 1.0D0/DFLOAT(MVAL+MVAL+1)
          DO M=1,LENGTH
            FM(M,K) = VALUE
          ENDDO
        ENDDO
        RETURN
C
C     ITYPE = 1: POWER SERIES EVALUATION (INITIALIZE SERIES FOR M = MMAX)
      ELSEIF(ITYPE.EQ.1) THEN
        DO M=1,LENGTH
          TEXP(M)      = DEXP(-T(M))
          TT2(M)       = 2.0D0*T(M)
          TMMAX(M)     = 1.0D0
          FM(M,MMAX+1) = 1.0D0
        ENDDO
C       LOOP OVER TERMS IN THE POWER SERIES
C       (CONVERGENCE IS ACHIEVED WHEN EXP(-T)*TERM(K)< 1.0D-14 WHERE
C       TERM(K) IS THE KTH TERM IN THE POWER SERIES)
        DO K=1,MSERIES
          DMMAX = DFLOAT(2*(MMAX+K)+1)
          DO M=1,LENGTH
            TMMAX(M)     = TMMAX(M)*(TT2(M)/DMMAX)
            FM(M,MMAX+1) = FM(M,MMAX+1)+TMMAX(M)
          ENDDO
        ENDDO
C       RESCALE BY THE PREFACTOR
        DEN = DFLOAT((2*MMAX)+1)
        DO M=1,LENGTH
          FM(M,MMAX+1) = FM(M,MMAX+1)*TEXP(M)/DEN
        ENDDO
C       NOW COMPLETE TABLE BY BACKWARDS RECURRENCE
        DO I=1,MMAX
          MIND  = MMAX-I+1
          MVAL  = MIND-1
          COEFF = DFLOAT(MVAL+MVAL+1)
          DO M=1,LENGTH
            FM(M,MIND) = (TT2(M)*FM(M,MIND+1)+TEXP(M))/COEFF
          ENDDO
        ENDDO
C
C     ITYPE = 2: ASYMPTOTIC EXPANSION WITH VERY LARGE ARGUMENT
      ELSEIF(ITYPE.EQ.2) THEN
C       INITIALIZE THE ASYMPTOTIC EXPANSION
        DO M=1,LENGTH
          TEXP(M)  = DEXP(-T(M))
          TT2(M)   = 2.0D0*T(M)
          TROOT(M) = DSQRT(T(M)) 
        ENDDO
        DO M=1,LENGTH
          FM(M,1) = A0/(B0+T(M))
        ENDDO
C       RESCALE BY THE PREFACTOR
        DO M=1,LENGTH
          FM(M,1) = (PIROOT/TROOT(M))-(TEXP(M)*FM(M,1))
        ENDDO
C       NOW COMPLETE TABLE BY FORWARD RECURRENCE
        DO MIND=1,MMAX
          MVAL  = MIND-1
          COEFF = DFLOAT(MVAL+MVAL+1)
          DO M=1,LENGTH
            FM(M,MIND+1) = (COEFF*FM(M,MIND)-TEXP(M))/TT2(M)
          ENDDO
        ENDDO
C
C     ITYPE = 3: ASYMPTOTIC EXPANSION WITH VERY LARGE ARGUMENT
      ELSEIF(ITYPE.EQ.3) THEN
C       INITIALIZE THE ASYMPTOTIC EXPANSION
        DO M=1,LENGTH
          TT2(M)  = 2.0D0*T(M)
          FM(M,1) = PIROOT/DSQRT(T(M))
        ENDDO
C       NOW COMPLETE TABLE BY FORWARD RECURRENCE
        DO MIND=1,MMAX
          MVAL  = MIND-1
          COEFF = DFLOAT(MVAL+MVAL+1)
          DO M=1,LENGTH
            FM(M,MIND+1) = (COEFF*FM(M,MIND))/TT2(M)
          ENDDO
        ENDDO
C
C     ITYPE OUT OF RANGE: INVALID INPUT TO FUNFM
      ELSE
        WRITE(6,901) ITYPE
901     FORMAT(2X,'In FUNFM: invalid type (must be 0-3)',I4)
        STOP
      ENDIF
C
      RETURN
      END



      SUBROUTINE GFINIT(MMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C*********************************************************************C
C                                                                     C
C          GGGGG   FFFFFF  IIIIII  NN    NN  IIIIII  TTTTTTTT         C
C         GG   GG  FF        II    NNN   NN    II       TT            C
C         GG       FF        II    NNNN  NN    II       TT            C
C         GG  GGG  FFFFFF    II    NN NN NN    II       TT            C
C         GG   GG  FF        II    NN  NNNN    II       TT            C
C         GG   GG  FF        II    NN   NNN    II       TT            C
C          GGGGG   FF      IIIIII  NN    NN  IIIIII     TT            C
C                                                                     C
C     GFINIT EVALUATES THE BOYS INCOMPLETE GAMMA FUNCTION INTEGRAL    C
C             F_M(T)=INT_{0}^{1} U^{2M} EXP(-TU^{2}) DU               C
C                                                                     C
C     FOR 0 < M < MMAX AND 0 < T < 18.0 in STEPS OF 0.01              C
C     USING A POWER SERIES REPRESENTATION. VALUES ARE STORED FOR USE  C
C     IN SUBSEQUENT EVALUATIONS USING THE TAYLOR SERIES               C
C     F_N(X+H) = F_N(X) - F_{N+1}(X)H + (1/2!)F_{N+2}(X)H^2 - ...     C
C ------------------------------------------------------------------- C
C     THE VALUES X ARE CHOSEN FOR A STEP 0 < H < 0.01                 C
C*********************************************************************C
      PARAMETER(NINTER=1801,ISTEP=100,MAXL=30,MSERIES=60)
      COMMON/FUNFM/FM(NINTER,MAXL),T(NINTER)
      DIMENSION TMMAX(NINTER),TT2(NINTER),TEXP(NINTER)
C*********************************************************************C
C     SPECIAL CASE FOR T=0                                            C
C*********************************************************************C
      T(1) = 0.0D0
      DO K=1,MMAX+1
        MVAL    = K-1
        VALUE   = 1.0D0/DFLOAT(MVAL+MVAL+1)
        FM(1,K) = VALUE
      ENDDO
C*********************************************************************C
C     POWER SERIES EVALUATION                                         C
C     INITIALIZE THE POWER SERIES FOR M = MMAX                        C
C*********************************************************************C
      DO M=2,NINTER
        T(M)         = DFLOAT(M-1)/DFLOAT(ISTEP)
        TEXP(M)      = DEXP(-T(M))
        TT2(M)       = 2.0D0*T(M)
        TMMAX(M)     = 1.0D0
        FM(M,MMAX+1) = 1.0D0
      ENDDO
C*********************************************************************C
C     LOOP OVER TERMS IN THE POWER SERIES                             C
C                                                                     C
C     CONVERGENCE IS ACHIEVED WHEN EXP(-T)*TERM(K)< 1.0E-14           C
C     WHERE TERM(K) IS THE KTH TERM IN THE POWER SERIES               C
C     NOTE THAT THE TERMS ARE ALWAYS POSITIVE SO THAT THERE IS        C
C     NO NEED TO TEST FOR ABSOLUTE VALUE                              C
C                                                                     C
C*********************************************************************C
      DO K=1,MSERIES
        DMMAX = DFLOAT(2*(MMAX+K)+1)
        DO M=2,NINTER
          TMMAX(M)     = TMMAX(M)*(TT2(M)/DMMAX)
          FM(M,MMAX+1) = FM(M,MMAX+1)+TMMAX(M)
        ENDDO
      ENDDO
C*********************************************************************C
C     RESCALE BY THE PREFACTOR                                        C
C*********************************************************************C
      DEN = DFLOAT((2*MMAX)+1)
      DO M=2,NINTER
        FM(M,MMAX+1) = FM(M,MMAX+1)*TEXP(M)/DEN
      ENDDO
C*********************************************************************C
C     NOW COMPLETE TABLE BY BACKWARDS RECURRENCE                      C
C*********************************************************************C
      DO I=1,MMAX
        MIND  = MMAX-I+1
        MVAL  = MIND-1
        COEFF = DFLOAT(MVAL+MVAL+1)
        DO M=2,NINTER
          FM(M,MIND) = (TT2(M)*FM(M,MIND+1)+TEXP(M))/COEFF
        ENDDO
      ENDDO
      RETURN
      END

