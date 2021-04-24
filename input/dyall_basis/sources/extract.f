      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NMAX=26)
C
      DIMENSION ES(NMAX),EP(NMAX),ED(NMAX),EF(NMAX)
C
      OPEN(UNIT=11,FILE='store.dat',STATUS='OLD')
      REWIND(UNIT=11)
      READ(11,*) NS,NP,ND,NF
      DO N=1,NS
        READ(11,*) NTMP,ES(N)
      ENDDO
      DO N=1,NP
        READ(11,*) NTMP,EP(N)
      ENDDO
      DO N=1,ND
        READ(11,*) NTMP,ED(N)
      ENDDO
      DO N=1,NF
        READ(11,*) NTMP,EF(N)
      ENDDO
      CLOSE(UNIT=11)

1     FORMAT(1P,D13.7)
2     FORMAT(I2)

      WRITE(*,2) NS
      DO N=1,NS
        WRITE(*,1) ES(N)
      ENDDO

      WRITE(*,2) NP
      DO N=1,NP
        WRITE(*,1) EP(N)
      ENDDO

      WRITE(*,2) ND
      DO N=1,ND
        WRITE(*,1) ED(N)
      ENDDO
      
      IF(NF.EQ.0) STOP
      
      WRITE(*,2) NF
      DO N=1,NF
        WRITE(*,1) EF(N)
      ENDDO
C
      STOP
      END
