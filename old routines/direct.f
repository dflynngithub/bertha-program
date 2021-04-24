C
C
      SUBROUTINE OPENFL(IRECL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C         OOOOOO  PPPPPPP  EEEEEEEE NN    NN FFFFFFFF LL               C
C        OO    OO PP    PP EE       NNN   NN FF       LL               C
C        OO    OO PP    PP EE       NNNN  NN FF       LL               C
C        OO    OO PP    PP EEEEEE   NN NN NN FFFFFF   LL               C
C        OO    OO PPPPPPP  EE       NN  NNNN FF       LL               C
C        OO    OO PP       EE       NN   NNN FF       LL               C
C         OOOOOO  PP       EEEEEEEE NN    NN FF       LLLLLLLL         C
C                                                                      C
C -------------------------------------------------------------------- C
C  OPENFL OPENS A NEW DIRECT-ACCESS FILE FOR RC-INTEGRALS.             C
C**********************************************************************C
C
      OPEN(UNIT=60,FILE='RCRAW',ACCESS='DIRECT',STATUS='NEW',
     &                                   RECL=IRECL,FORM='UNFORMATTED')
C
      RETURN
      END
C
C
      SUBROUTINE CLSEFL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C          CCCCCC  LL       SSSSSS  EEEEEEEE FFFFFFFF LL               C
C         CC    CC LL      SS    SS EE       FF       LL               C
C         CC       LL      SS       EE       FF       LL               C
C         CC       LL       SSSSSS  EEEEEE   FFFFFF   LL               C
C         CC       LL            SS EE       FF       LL               C
C         CC    CC LL      SS    SS EE       FF       LL               C
C          CCCCCC  LLLLLLLL SSSSSS  EEEEEEEE FF       LLLLLLLL         C
C                                                                      C
C -------------------------------------------------------------------- C
C  CLSEFL CLOSES AND DELETES A DIRECT-ACCESS FILE FOR RC-INTEGRALS.    C
C**********************************************************************C
C
      CLOSE(UNIT=60,STATUS='DELETE')
C
      RETURN
      END
C
C
      SUBROUTINE WRTEFL(RC,ISTRT,NREC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C      WW         WW RRRRRRR TTTTTTTT EEEEEEEE FFFFFFFF LL             C
C      WW         WW RR    RR   TT    EE       FF       LL             C
C      WW         WW RR    RR   TT    EE       FF       LL             C
C      WW    W    WW RR    RR   TT    EEEEEE   FFFFFF   LL             C
C      WW   WWW   WW RRRRRRR    TT    EE       FF       LL             C
C       WW WW WW WW  RR    RR   TT    EE       FF       LL             C
C        WWW   WWW   RR    RR   TT    EEEEEEEE FF       LLLLLLLL       C
C                                                                      C
C -------------------------------------------------------------------- C
C  WRTEFL WRITES TO A DIRECT-ACCESS FILE FOR RC-INTEGRALS.             C
C**********************************************************************C
C
      ISTOP = ISTRT+NREC
      WRITE(UNIT=60,REC=NREC) (RC(I),I=ISTRT,ISTOP)
C
      RETURN
      END
C
C
      SUBROUTINE READFL(RC,ISTRT,NREC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C        RRRRRRR  EEEEEEEE    AA    DDDDDDD  FFFFFFFF LL               C
C        RR    RR EE         AAAA   DD    DD FF       LL               C
C        RR    RR EE        AA  AA  DD    DD FF       LL               C
C        RR    RR EEEEEE   AA    AA DD    DD FFFFFF   LL               C
C        RRRRRRR  EE       AAAAAAAA DD    DD FF       LL               C
C        RR    RR EE       AA    AA DD    DD FF       LL               C
C        RR    RR EEEEEEEE AA    AA DDDDDDD  FF       LLLLLLLL         C
C                                                                      C
C -------------------------------------------------------------------- C
C  READFL READS FROM A DIRECT-ACCESS FILE FOR RC-INTEGRALS.            C
C**********************************************************************C
C
      ISTOP = ISTRT+NREC
      READ(UNIT=60,REC=NREC) (RC(I),I=ISTRT,ISTOP)
C
      RETURN
      END

