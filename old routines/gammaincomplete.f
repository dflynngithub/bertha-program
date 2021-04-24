C
C
      FUNCTION GAMMLN(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      GGGGGG     AA    MM       MM MM       MM LL       NN    NN      C
C     GG    GG   AAAA   MMM     MMM MMM     MMM LL       NNN   NN      C
C     GG        AA  AA  MMMM   MMMM MMMM   MMMM LL       NNNN  NN      C
C     GG       AA    AA MM MM MM MM MM MM MM MM LL       NN NN NN      C
C     GG   GGG AAAAAAAA MM  MMM  MM MM  MMM  MM LL       NN  NNNN      C
C     GG    GG AA    AA MM   M   MM MM   M   MM LL       NN   NNN      C
C      GGGGGG  AA    AA MM       MM MM       MM LLLLLLLL NN    NN      C
C                                                                      C
C -------------------------------------------------------------------- C
C  GAMMLN RETURNS LN(GAMMA(X)) FOR X>0. FULL ACCURACY IS OBTAINED FOR  C
C  X>1. FOR 0<X<1, THE REFLECTION FORMULA CAN BE USED FIRST.           C
C**********************************************************************C
      DIMENSION COF(6)
C
      DATA COF/76.18009173D0,-86.50532033D0,24.01409822D0,
     &        -1.231739516D0,0.120858003D-2,-0.536382D-5/
      DATA STP/2.50662827465D0/
C
      Y   = X-1.0D0
      TMP = Y+5.5D0
      TMP = TMP*DLOG(TMP)-TMP
      SER = 1.0D0
C
      DO J=1,6
        Y   = Y+1.0D0
        SER = SER+COF(J)/Y
      ENDDO
C
      GAMMLN = TMP + DLOG(STP*SER)
C
      RETURN
      END
C
C
      FUNCTION GAMLWR(A,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      GGGGGG     AA    MM       MM LL      WW         WW RRRRRRR      C
C     GG    GG   AAAA   MMM     MMM LL      WW         WW RR    RR     C
C     GG        AA  AA  MMMM   MMMM LL      WW         WW RR    RR     C
C     GG       AA    AA MM MM MM MM LL      WW    W    WW RR    RR     C
C     GG   GGG AAAAAAAA MM  MMM  MM LL       WW  WWW  WW  RRRRRRR      C
C     GG    GG AA    AA MM   M   MM LL        WWWW WWWW   RR    RR     C
C      GGGGGG  AA    AA MM       MM LLLLLLLL   WW   WW    RR    RR     C
C                                                                      C
C -------------------------------------------------------------------- C
C  GAMLWR RETURNS THE LOWER INCOMPLETE GAMMA FUNCTION GAM(A)*P(A,X).   C
C**********************************************************************C
C
C     INPUT CHECK: MAKE SURE X>0.
      IF(X.LT.0.0D0) THEN
        WRITE(6, *) 'In GAMLWR: X<0.0D0'
        WRITE(7, *) 'In GAMLWR: X<0.0D0'
      ENDIF
C
C     INPUT CHECK: MAKE SURE A >0.
      IF(A.LE.0.0D0) THEN
        WRITE(6, *) 'In GAMLWR: A<=0.0D0'
        WRITE(7, *) 'In GAMLWR: A<=0.0D0'
      ENDIF
C
C     GAMMA FUNCTION VALUE
      GLN = GAMMLN(A)
C
C     REGULARISED INCOMPLETE GAMMA FUNCTION VALUE
      IF(X.LT.A+1.0D0) THEN
C
C       USE THE SERIES REPRESENTATION
        GAMLWR = CDFSER(A,X)
C
      ELSE
C
C       USE THE CONTINUED FRACTION REPRESENTATION
        GAMLWR = 1.0D0-CDFFRC(A,X)
C
      ENDIF
C
C     gamma(A,X) = GAMMA(A)*P(A,X)
      GAMLWR = DEXP(GLN)*GAMLWR
C
      RETURN
      END
C
C
      FUNCTION GAMUPR(A,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        GGGGGG     AA    MM       MM UU    UU PPPPPPP  RRRRRRR        C
C       GG    GG   AAAA   MMM     MMM UU    UU PP    PP RR    RR       C
C       GG        AA  AA  MMMM   MMMM UU    UU PP    PP RR    RR       C
C       GG       AA    AA MM MM MM MM UU    UU PP    PP RR    RR       C
C       GG   GGG AAAAAAAA MM  MMM  MM UU    UU PPPPPPP  RRRRRRR        C
C       GG    GG AA    AA MM   M   MM UU    UU PP       RR    RR       C
C        GGGGGG  AA    AA MM       MM  UUUUUU  PP       RR    RR       C
C                                                                      C
C -------------------------------------------------------------------- C
C  GAMMP RETURNS THE UPPER INCOMPLETE GAMMA FUNCTION GAM(A)*Q(A,X).    C
C**********************************************************************C
C
C     INPUT CHECK: MAKE SURE X>0.
      IF(X.LT.0.0D0) THEN
        WRITE(6, *) 'In GAMUPR: X<0.0D0'
        WRITE(7, *) 'In GAMUPR: X<0.0D0'
      ENDIF
C
C     INPUT CHECK: MAKE SURE A >0.
      IF(A.LE.0.0D0) THEN
        WRITE(6, *) 'In GAMUPR: A<=0.0D0'
        WRITE(7, *) 'In GAMUPR: A<=0.0D0'
      ENDIF
C
C     GAMMA FUNCTION VALUE
      GLN = GAMMLN(A)
C
C     REGULARISED INCOMPLETE GAMMA FUNCTION VALUE
      IF(X.LT.A+1.0D0) THEN
C
C       USE THE SERIES REPRESENTATION
        GAMUPR = 1.0D0-CDFSER(A,X)
C
      ELSE
C
C       USE THE CONTINUED FRACTION REPRESENTATION
        GAMUPR = CDFFRC(A,X)
C
      ENDIF
C
C     GAMMA(A,X) = GAMMA(A)*Q(A,X)
      GAMUPR = DEXP(GLN)*GAMUPR
C
      RETURN
      END
C
C
      FUNCTION CDFSER(A,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          CCCCCC  DDDDDDD  FFFFFFFF SSSSSS  EEEEEEEE RRRRRRR          C
C         CC    CC DD    DD FF      SS    SS EE       RR    RR         C
C         CC       DD    DD FF      SS       EE       RR    RR         C
C         CC       DD    DD FFFFFF   SSSSSS  EEEEEE   RR    RR         C
C         CC       DD    DD FF            SS EE       RRRRRRR          C
C         CC    CC DD    DD FF      SS    SS EE       RR    RR         C
C          CCCCCC  DDDDDDD  FF       SSSSSS  EEEEEEEE RR    RR         C
C                                                                      C
C -------------------------------------------------------------------- C
C  CDFSER RETURNS THE CUMULATIVE DISTRIBUTION FUNCTION (ALSO CALLED    C
C  THE REGULARISED INCOMPLETE GAMMA FUNCTION) BY SERIES REPRESENTATION.C
C**********************************************************************C
C
C     PARAMETER WARNINGS
      IF(X.LE.0.0D0) THEN
        WRITE(6, *) 'In CDFSER: X < 0.0'
        WRITE(7, *) 'In CDFSER: X < 0.0'
        GAMSER = 0.0D0
        RETURN
      ENDIF
C
C     LOG OF GAMMA FUNCTION
      GLN = GAMMLN(A)
C
C     PREPARE INITIAL VALUES
      AP  = A
      SUM = 1.0D0/A
      DEL = SUM
C
C     LOOP UNTIL RESULT IS CONSISTENT
      DO N=1,100
C
C       UPDATE VALUE OF FUNCTION
        AP  = AP+1.0D0
        DEL = DEL*X/AP
        SUM = SUM + DEL
C
C       CHECK FOR AGREEMENT
        IF(DABS(DEL).LT.DABS(SUM)*3.0D-7) THEN
          GOTO 10
        ENDIF
C
      ENDDO
C
C     FAILURE TO CONVERGE TO A GOOD ANSWER
      WRITE(6, *) 'In CDFSER: A too large. Need more iterations.'
      WRITE(7, *) 'In CDFSER: A too large. Need more iterations.'
      STOP
C
C     SUCCESSFUL EXIT
10    CONTINUE
      CDFSER = SUM*DEXP(-X+A*DLOG(X)-GLN)
C
      RETURN
      END
C
C
      FUNCTION CDFFRC(A,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         CCCCCC  DDDDDDD  FFFFFFFF FFFFFFFF RRRRRRR   CCCCCC          C
C        CC    CC DD    DD FF       FF       RR    RR CC    CC         C
C        CC       DD    DD FF       FF       RR    RR CC               C
C        CC       DD    DD FFFFFF   FFFFFF   RR    RR CC               C
C        CC       DD    DD FF       FF       RRRRRRR  CC               C
C        CC    CC DD    DD FF       FF       RR    RR CC    CC         C
C         CCCCCC  DDDDDDD  FF       FF       RR    RR  CCCCCC          C
C                                                                      C
C -------------------------------------------------------------------- C
C  CDFSER RETURNS THE CUMULATIVE DISTRIBUTION FUNCTION (ALSO CALLED    C
C  THE REGULARISED INCOMPLETE GAMMA FUNCTION) BY CONTINUED FRACTION.   C
C -------------------------------------------------------------------- C
C  DFNOTE: THIS IS EQUIVALENT TO THE ROUTINE IN `NUMERICAL RECIPES',   C
C          BUT BOTH PRODUCE THE WRONG RESULT.                          C
C**********************************************************************C
C
C     LOG OF GAMMA FUNCTION
      GLN = GAMMLN(A)
C
C     INITIALISE VALUES
      A0  = 1.0D0
      A1  = X
      B0  = 0.0D0
      B1  = 1.0D0
      FAC = 1.0D0
      TST = 0.0D0
C
C     LOOP UNTIL RESULT IS CONSISTENT
      DO N=1,100
C
C       INTERMEDIATE VALUES
        RN = DFLOAT(N)
        AN = RN-A
        A0 = FAC*(A1+A0*AN)
        B0 = FAC*(B1+B0*AN)
        FN = FAC*RN
        A1 = X*A0 + FN*A1
        B1 = X*B0 + FN*B1
C
C       TEST WHETHER TO RENORMALISE
        IF(A1.NE.0.0D0) THEN
C
          FAC = 1.0D0/A1
          SUM = B1*FAC
C
C         CHECK FOR AGREEMENT
          IF(DABS((SUM-TST)/SUM).LT.3.0D-7) THEN
            GOTO 10
          ENDIF
C
C         SAVE MOST RECENT GAMMA VALUE
          TST = SUM
C
        ENDIF
C
      ENDDO
C
C     FAILURE TO CONVERGE TO A GOOD ANSWER
      WRITE(6, *) 'In CDFFRC: A too large. Need more iterations.'
      WRITE(7, *) 'In CDFFRC: A too large. Need more iterations.'
C
C     SUCCESSFUL EXIT
10    CONTINUE
      CDFFRC = SUM*DEXP(-X+A*DLOG(X)-GLN)
C
      RETURN
      END

