      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NMAX=26)
C
      DIMENSION ES(NMAX),EP(NMAX),ED(NMAX),EF(NMAX),EG(NMAX),EH(NMAX)
C
      OPEN(UNIT=11,FILE='store.dat',STATUS='OLD')
      REWIND(UNIT=11)

      READ(11,*) 
      READ(11,*) 
      READ(11,*) 
      READ(11,*) NS
      DO N=1,NS
        READ(11,*) ES(N)
      ENDDO

      READ(11,*) 
      READ(11,*) NP
      IF(NP.EQ.0) GOTO 100
      DO N=1,NP
        READ(11,*) EP(N)
      ENDDO

      READ(11,*) 
      READ(11,*) ND
      IF(ND.EQ.0) GOTO 100
      DO N=1,ND
        READ(11,*) ED(N)
      ENDDO

      READ(11,*) 
      READ(11,*) NF
      IF(NF.EQ.0) GOTO 100
      DO N=1,NF
        READ(11,*) EF(N)
      ENDDO

      READ(11,*) 
      READ(11,*) NG
      IF(NG.EQ.0) GOTO 100
      DO N=1,NG
        READ(11,*) EG(N)
      ENDDO

      READ(11,*) 
      READ(11,*) NH
      IF(NH.EQ.0) GOTO 100
      DO N=1,NH
        READ(11,*) EH(N)
      ENDDO

100   CONTINUE
      CLOSE(UNIT=11)

1     FORMAT(1P,D14.8)
2     FORMAT(I1)
3     FORMAT(I2)

      IF(NS.LT.10) THEN
        WRITE(*,2) NS
      ELSEIF(NS.GE.10) THEN
        WRITE(*,3) NS
      ENDIF
      DO N=1,NS
        WRITE(*,1) ES(N)
      ENDDO

      IF(NP.EQ.0) THEN
        STOP
      ELSEIF(NP.GT.0.AND.NP.LT.10) THEN
        WRITE(*,2) NP
      ELSEIF(NP.GE.10) THEN
        WRITE(*,3) NP
      ENDIF
      DO N=1,NP
        WRITE(*,1) EP(N)
      ENDDO

      IF(ND.EQ.0) THEN
        STOP
      ELSEIF(ND.GT.0.AND.ND.LT.10) THEN
        WRITE(*,2) ND
      ELSEIF(ND.GE.10) THEN
        WRITE(*,3) ND
      ENDIF
      DO N=1,ND
        WRITE(*,1) ED(N)
      ENDDO
      
      IF(NF.EQ.0) THEN 
        STOP
      ELSEIF(NF.GT.0.AND.NF.LT.10) THEN
        WRITE(*,2) NF
      ELSEIF(NF.GE.10) THEN
        WRITE(*,3) NF
      ENDIF
      DO N=1,NF
        WRITE(*,1) EF(N)
      ENDDO

      IF(NG.EQ.0) THEN
        STOP
      ELSEIF(NG.GT.0.AND.NG.LT.10) THEN
        WRITE(*,2) NG
      ELSEIF(NG.GE.10) THEN
        WRITE(*,3) NG
      ENDIF
      DO N=1,NG
        WRITE(*,1) EG(N)
      ENDDO

      IF(NH.EQ.0) THEN
        STOP
      ELSEIF(NH.GT.0.AND.NH.LT.10) THEN
        WRITE(*,2) NH
      ELSEIF(NH.GE.10) THEN
        WRITE(*,3) NH
      ENDIF
      DO N=1,NH
        WRITE(*,1) EH(N)
      ENDDO
C
      STOP
      END
