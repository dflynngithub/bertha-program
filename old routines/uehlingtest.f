C
C     RE-CALCULATE CHI_N(X)
C     CALL CHIGEN
C
C     LIMITS OF NUCLEAR CHARGE RADIUS (FM TO ATOMIC UNITS)
      R0 = 1.0D0/CFM
      RN = 6.0D0/CFM
      NT = 1000
      HN = (RN-R0)/DFLOAT(NT)
      
      WRITE(*,*) 'unifrm_exct'
      
      k = 0
C
C     LOOP OVER NUCLEAR RADIUS VALUES
      DO I=0,NT
C
C       NUCLEAR RADIUS
        RNUC(1) = R0 + HN*DFLOAT(I)
C
C       ATOMIC MASS
        ANUC(1) = RNUC(1)*CFM - 0.57D0
        ANUC(1) = (ANUC(1)/0.836D0)**3
C
        IF(DSQRT(1.4D0)*PI*0.25D0*2.3D0/THLG.GT.RNUC(1)*CFM) THEN
C       KICK OUT NUCLEI THAT DON'T SATISFY MINIMUM CONDITIONS
          GOTO 20
        ELSEIF(ANUC(1).GT.300.0D0) THEN
C       KICK OUT NUCLEI THAT ARE A BIT TOO BIG
          GOTO 20
        ENDIF
C
C       SKIN THICKNESS AND RELATED PARAMETER (FM TO ATOMIC UNITS)
        TFMI(1) = 2.30D0/CFM
        AFMI(1) = 0.25D0*TFMI(1)/THLG
C
C       HALF-THICKNESS
        CFMI(1) = 5.0D0*RNUC(1)*RNUC(1)/3.0D0
     &          - 7.0D0*PI*PI*AFMI(1)*AFMI(1)/3.0D0
        CFMI(1) = DSQRT(CFMI(1))
C
        IF(AFMI(1)/CFMI(1).GT.1.0D0) GOTO 20
        
        k = k+1

c        if(k.ge.1675.and.k.le.1696) then
C        CALL VUEHLNG(1)
C        CALL EUEHLNG(1)
C        CALL PUEHLNG(1)
        CALL PUEHUNI(1)

C       CALL RHOASYM(1)
c        endif

C
20      CONTINUE
C
      ENDDO
      STOP
      
      
      
C
C
c      OPEN(UNIT=10,FILE='plots/abgparams.dat',STATUS='REPLACE')
c      OPEN(UNIT=11,FILE='plots/exponents.dat',STATUS='REPLACE')
c      OPEN(UNIT=12,FILE='plots/amplituds.dat',STATUS='REPLACE')
c
c        DO IA=1,300
c          
c        ANUC(1) = DFLOAT(IA)
c
cC       NUCLEAR CHARGE FITTING ROUTINE, GIVEN # FUNCTIONS AND IWRT (0/1)
        CALL NUCCHRG(ICNT,10,0)
c
c        EPS = DSQRT(1.4D0)*PI*0.25D0*2.3D0/THLG
c        IF(EPS.GT.RNUC(1)*CFM) THEN
c          GOTO 200
c        ENDIF
cC
c        IF(AFMI(1)/CFMI(1).GT.1.0D0) GOTO 200
c        IF(NNUC(1).EQ.1) GOTO 200
c          WRITE(*,*) '----> WORKING ON ',IA,ANUC(1),RNUC(1)*CFM
cC
cc         OPTIONAL UEHLING FITTING ROUTINE
C          IF(HMLT.EQ.'DHFQ') THEN
C            CALL NUCUEHL(ICNT,26,1,AR,BF,GF,RS)
C          ENDIF
cC
c90        FORMAT(I3,3X,F10.8,3X,ES15.8,3X,F8.4,3X,F12.8,3X,F16.14)
c          WRITE(10,90) IA,RNUC(1)*CFM,AR,BF,GF,RS
c          WRITE(* ,90) IA,RNUC(1)*CFM,AR,BF,GF,RS
c          WRITE(11, *) IA,(XUEH(1,IFT)*(RNUC(1)**2),IFT=1,26)
c          WRITE(* , *) IA,(XUEH(1,IFT)*(RNUC(1)**2),IFT=1,26)
c          WRITE(12, *) IA,(FUEH(1,IFT),IFT=1,26)
c          WRITE(* , *) IA,(FUEH(1,IFT),IFT=1,26)
c
c200       CONTINUE
c        
c        ENDDO
cC
c      CLOSE(UNIT=12)
c      CLOSE(UNIT=11)
c      CLOSE(UNIT=10)
cc
c      STOP

