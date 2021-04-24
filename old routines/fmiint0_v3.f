      FUNCTION FMIINT0(L,EIJ,A,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**********************************************************************C
C                                                                      C
C       FFFFFFFF MM       MM IIII IIII NN    NN TTTTTTTT 000000        C
C       FF       MMM     MMM  II   II  NNN   NN    TT   00   000       C
C       FF       MMMM   MMMM  II   II  NNNN  NN    TT   00  0000       C
C       FFFFFF   MM MM MM MM  II   II  NN NN NN    TT   00 00 00       C
C       FF       MM  MMM  MM  II   II  NN  NNNN    TT   0000  00       C
C       FF       MM   M   MM  II   II  NN   NNN    TT   000   00       C
C       FF       MM       MM IIII IIII NN    NN    TT    000000        C
C                                                                      C
C -------------------------------------------------------------------- C
C  FMIINT0 CALCULATES AN AUXILLIARY FERMI POTENTIAL INTEGRAL OVER      C
C  A PAIR OF ATOMIC BASIS FUNCTIONS, WHERE L IS A POSITIVE INTEGER.    C
C                                                                      C
C    FMIINT0(L,λ,ζ) = ∫^{∞}_{0} r^2L+1 exp(-λ r^2) r*V_fermi(r) dr.    C
C                                                                      C
C  DFNOTE: THIS VERSION WAS USED FOR TESTING PURPOSES, SINCE THE       C
C          CUMBERSOME TERMS NEED A SERIES EXPANSION OF SOME KIND.      C
C**********************************************************************C
      INCLUDE 'param.h'
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In FMIINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In FMIINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
C     FACTORS NEEDED FOR ALL PARAMETERS N
      U   = A/C
      U3  = U*U*U
      X   = EIJ*C*C
      X12 = DSQRT(X)
      C2L = C*C
      C2L = C2L**(L+1)
C
C     FERMI NORMALISATION CONSTANT
      FNRM = 1.0D0 + PI*PI*U*U - 6.0D0*U3*POLYLOG(3,-1.0D0/U)
      FNRM = 1.0D0/FNRM
C
C     INTEGRAL TYPE: INCOMPLETE GAMMA FUNCTION
      X1A = GAMLWR(2*L+3,X)*(3.0D0+PI*PI*U*U)/(2.0D0*X12)
      X1B =-GAMLWR(2*L+5,X)/(2.0D0*X12*X12*X12)
      X1C = GAMUPR(2*L+2,X)*(1.0D0+PI*PI*U*U)
C
C     INTEGRAL TYPE: LOGARITHMIC INTEGRAL
      X2A =-6.0D0*GAMHLF(2*L)*U3*POLYLOG(3,-1.0D0/U)
C
C     TRANSFER DATA TO FMIINT0
      X1 = X1A+X1B+X1C
      X2 = X2A
      FMIINT0 = 0.5D0*FNRM*(X1+X2)/(EIJ**(L+1))
C
C     ALGORITHMS AVAILABLE FOR THE REMAINING SET OF INTEGRALS
c      X3A = (DEXP(-X)/X - 1.0D0/(X*DSQRT(X)))*POLYLOG(2,0.0D0)
c      X3B = 2.0D0*U*U*POLYLOG(4,0.0D0)/DSQRT(X)
c      X3C = (1.0D0/DSQRT(X)-1.0D0)*POLYLOG(1,-1.0D0/U)/(4.0D0*U*X*X)
c      X3D = U3*POLYLOG(5,-1.0D0/U)/DSQRT(X)
c       
c      X3  = 3.0D0*A*A*(X3A - X3B + X3C - X3D)
cC
c      X4A = 2.0D0*U*POLYLOG(4,0.0D0)
c      X4B = U*U*POLYLOG(5,-1.0D0/U)
c       
c      X4  = 6.0D0*A*A*U*(X4A + X4B)
cC      
c      FMIINT0 = FMIINT0 + FNRM*(X3+X4)

cC
cC     ALGORITHMS AVAILABLE FOR THE REMAINING SET OF INTEGRALS
c
c      X3A = (DEXP(-X)/X - 1.0D0/(X*DSQRT(X)))*POLYLOG(2,0.0D0)
c      X3B = 2.0D0*U*U*POLYLOG(4,0.0D0)/DSQRT(X)
c      X3C = (1.0D0/DSQRT(X)-1.0D0)*POLYLOG(1,-1.0D0/U)/(4.0D0*U*X*X)
c      X3D = U3*POLYLOG(5,-1.0D0/U)/DSQRT(X)
c       
c      X3  = 3.0D0*A*A*(X3A - X3B + X3C - X3D)
cC
c      X4A = 2.0D0*U*POLYLOG(4,0.0D0)
c      X4B = U*U*POLYLOG(5,-1.0D0/U)
c       
c      X4  = 6.0D0*A*A*U*(X4A + X4B)
cC      
c      FMIINT0 = FMIINT0 + FNRM*(X3+X4)
c
C     ALGORITHMS AVAILABLE FOR THE REMAINING SET OF INTEGRALS

c      Y1 = 4.0D0*U3*POLYLOG(6,-1.0D0/U)
c      Y2 = 36.0D0*X*U3*U*U*POLYLOG(8,-1.0D0/U)
c      Y3 = (DEXP(-X) - 1.0D0 + X - 0.5D0*X*X)*POLYLOG(2,0.0D0)/X
c      Y4 = 4.0D0*U*U*(2.0D0-3.0D0*X)*POLYLOG(4,0.0D0)
c      Y5 = 72.0D0*U3*U*X*POLYLOG(6,0.0D0)

c      y1  =-72.0d0*u3*u*x*polylog(6,0.0d0)
c      y2  =  4.0d0*u*u*(2.0d0-3.0d0*x)*polylog(4,0.0d0)
c      y3  =  4.0d0*u3*polylog(5,-1.0d0/u)
c      y4  =-36.0d0*u3*u*u*x*polylog(7,-1.0d0/u)
c
c      y1  = 960.0d0*u3*u3*x*x*polylog(8,0.0d0)
c      y2  =   8.0d0*u3*u*x*(-9.0d0+20.0d0*x)*polylog(6,0.0d0)
c      y3  =   4.0d0*u*u*(2.0d0-3.0d0*x+2.0d0*x*x)*polylog(4,0.0d0)
c      y4  =   4.0d0*u3*polylog(5,-1.0d0/u)
c      y5  = -36.0d0*u3*u*u*x*polylog(7,-1.0d0/u)
c      y6  = 480.0d0*u3*u3*u*x*x*polylog(9,-1.0d0/u)
c
c      Y34  = 3.0D0*A*A*(Y1+Y2+Y3+Y4+Y5+Y6)

c      if(l.eq.0) then
c        y1  =     4.0d0*u*u*(2.0d0-3.0d0*x+2.0d0*x*x-1.0d0*x*x*x/12.0d0)
c     &                                                *polylog(4,0.0d0)
c        y2  =     8.0d0*u3*u*x*(-9.0d0+20.0d0*x-1.75d0*x*x)
c     &                                                *polylog(6,0.0d0)
c        y3  =    64.0d0*u3*u3*x*x*(15.0d0-35.0d0*x/8.0d0)
c     &                                                *polylog(8,0.0d0)
c        y4  =-16800.0d0*u3*u3*u*u*x*x*x*polylog(10,0.0d0)
c        y5  =     4.0d0*u3*polylog(5,-1.0d0/u)
c        y6  =   -36.0d0*u3*u*u*x*polylog(7,-1.0d0/u)
c        y7  =   480.0d0*u3*u3*u*x*x*polylog(9,-1.0d0/u)
c        y8  = -8400.0d0*u3*u3*u3*x*x*x*polylog(11,-1.0d0/u)
c
c        Y34  = 3.0D0*A*A*(Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8)
      u02 = u*u
      u04 = u*u*u02
      u06 = u*u*u04
      u08 = u*u*u06
      u10 = u*u*u08
      u12 = u*u*u10
      u14 = u*u*u12
      u16 = u*u*u14
      u18 = u*u*u16

      x0  = 1.0d0
      x1  = x*x0
      x2  = x*x1
      x3  = x*x2
      x4  = x*x3
      x5  = x*x4
      x6  = x*x5
      x7  = x*x6

      p04 = polylog( 4,0.0d0)
      p06 = polylog( 6,0.0d0)
      p08 = polylog( 8,0.0d0)
      p10 = polylog(10,0.0d0)
      p12 = polylog(12,0.0d0)
      p14 = polylog(14,0.0d0)
      p16 = polylog(16,0.0d0)
      p18 = polylog(18,0.0d0)
      
      s05 = polylog( 5,-1.0d0/u)
      s07 = polylog( 7,-1.0d0/u)
      s09 = polylog( 9,-1.0d0/u)
      s11 = polylog(11,-1.0d0/u)
      s13 = polylog(13,-1.0d0/u)
      s15 = polylog(15,-1.0d0/u)
      s17 = polylog(17,-1.0d0/u)
      s19 = polylog(19,-1.0d0/u)
      s21 = polylog(21,-1.0d0/u)
      s23 = polylog(23,-1.0d0/u)
      
      y04 = 0.0d0
      y06 = 0.0d0
      y08 = 0.0d0
      y10 = 0.0d0
      y12 = 0.0d0
      y14 = 0.0d0
      y16 = 0.0d0
      y18 = 0.0d0

      e05 = 0.0d0
      e07 = 0.0d0
      e09 = 0.0d0
      e11 = 0.0d0
      e13 = 0.0d0
      e15 = 0.0d0
      e17 = 0.0d0
      e19 = 0.0d0
      e21 = 0.0d0
      e23 = 0.0d0

      Y34 = 0.0D0
      if(l.eq.0.and.u*dsqrt(x).lt.0.1d0) then
c
        y04 = 8.0d0 - 12.0d0*x1 + 8.0d0*x2 - 10.0d0*x3/3.0d0
     &      + x4 - 7.0d0*x5/30.0d0 + 2.0d0*x6/45.0d0 - x7/140.0d0
        y06 =-72.0d0 + 160.0d0*x1 - 140.0d0*x2 + 72.0d0*x3 
     &      - 77.0d0*x4/3.0d0 + 104.0d0*x5/15.0d0 - 3.0d0*x6/2.0d0
        y08 = 960.0d0 - 2800.0d0*x1 + 3024.0d0*x2 - 1848.0d0*x3 
     &      + 2288.0d0*x4/3.0d0 - 234.0d0*x5
        y10 =-16800.0d0 + 60480.0d0*x1 - 77616.0d0*x2 + 54912.0d0*x3 
     &      - 25740.0d0*x4
        y12 = 362880.0d0 - 1552320.0d0*x1 + 2306304.0d0*x2 
     &      - 1853280.0d0*x3
        y14 =-9313920.0d0 + 46126080.0d0*x1 - 77837760.0d0*x2
        y16 = 276756480.0d0 - 1556755200.0d0*x1
        y18 =-9340531200.0d0
        
        yp   = u02*x0*y04*p04 + u04*x1*y06*p06 + u06*x2*y08*p08
     &       + u08*x3*y10*p10 + u10*x4*y12*p12 + u12*x5*y14*p14
     &       + u14*x6*y16*p16 + u16*x7*y18*p18

        e05 = 1.0d0
        e07 =-9.0d0*x1*u02
        e09 = 120.0d0*x2*u04
        e11 =-2100.0d0*x3*u06
        e13 = 45360.0d0*x4*u08
        e15 =-1164240.0d0*x5*u10
        e17 = 34594560.0d0*x6*u12
        e19 =-1167566400.0d0*x7*u14
        
        ye  = e05*s05 + e07*s07 + e09*s09 + e11*s11
     &      + e13*s13 + e15*s15 + e17*s17 + e19*s19

        Y34  = 3.0D0*A*A*(yp + 4.0d0*u02*u*ye)
C
      elseif(l.eq.1.and.u*dsqrt(x).lt.0.1d0) then
c
        y04 = 12.0d0 - 16.0d0*x1 + 10.0d0*x2 - 4.0d0*x3
     &      + 7.0d0*x4/6.0d0 - 4.0d0*x5/15.0d0 + x6/20.0d0 - x7/126.0d0
        y06 = 72.0d0 - 320.0d0*x1 + 420.0d0*x2 - 288.0d0*x3 
     &      + 385.0d0*x4/3.0d0 - 208.0d0*x5/5.0d0 + 21.0d0*x6/2.0d0
        y08 =-1920.0d0 + 8400.0d0*x1 - 12096.0d0*x2 + 9240.0d0*x3 
     &      - 4576.0d0*x4 + 1638.0d0*x5
        y10 = 50400.0d0 - 241920.0d0*x1 + 388080.0d0*x2 - 329472.0d0*x3 
     &      + 180180.0d0*x4
        y12 =-1451520.0d0 + 7761600.0d0*x1 - 13837824.0d0*x2 
     &      + 12972960.0d0*x3
        y14 = 46569600.0d0 - 276756480.0d0*x1 + 544864320.0d0*x2
        y16 =-1660538880.0d0 + 10897286400.0d0*x1
        y18 = 65383718400.0d0
        
        yp   = u02*x0*y04*p04 + u04*x0*y06*p06 + u06*x1*y08*p08
     &       + u08*x2*y10*p10 + u10*x3*y12*p12 + u12*x4*y14*p14
     &       + u14*x5*y16*p16 + u16*x6*y18*p18

        e07 = 3.0d0
        e09 =-80.0d0*x1*u02
        e11 = 2100.0d0*x2*u04
        e13 =-60480.0d0*x3*u06
        e15 = 1940400.0d0*x4*u08
        e17 =-69189120.0d0*x5*u10
        e19 = 2724321600.0d0*x6*u12
        e21 =-117621505000.0d0*x7*u14
        
        ye  = e07*s07 + e09*s09 + e11*s11 + e13*s13
     &      + e15*s15 + e17*s17 + e19*s19 + e21*s21

        Y34  = 3.0D0*U*U*C2L*(yp + 12.0d0*u04*u*ye)
C
      elseif(l.eq.2.and.u*dsqrt(x).lt.0.1d0) then
c
        y04 = 16.0d0 - 20.0d0*x1 + 12.0d0*x2 - 14.0d0*x3/3.0d0
     &      + 3.0d0*x4/10.0d0 - 3.0d0*x5/10.0d0 + x6/18.0d0
     &      - 11.0d0*x7/1260.0d0
        y06 = 320.0d0 - 840.0d0*x1 + 864.0d0*x2 - 1540.0d0*x3/3.0d0
     &      + 208.0d0*x4 - 63.0d0*x5 + 136.0d0*x6/9.0d0
        y08 = 1920.0d0 - 16800.0d0*x1 + 36288.0d0*x2 - 36960.0d0*x3 
     &      + 22880.0d0*x4 - 9828.0d0*x5
        y10 =-100800.0d0 + 825760.0d0*x1 - 1552320.0d0*x2
     &      + 1647360.0d0*x3 - 1081080.0d0*x4
        y12 = 4354560.0d0 - 31046400.0d0*x1 + 69189120.0d0*x2 
     &      - 77837760.0d0*x3
        y14 =-186278400.0d0 + 1383782400.0d0*x1 - 3269185920.0d0*x2
        y16 = 8302694400.0d0 - 65383708400.0d0*x1
        y18 =-392302310400.0d0
        
        yp   = u02*x0*y04*p04 + u04*x0*y06*p06 + u06*x0*y08*p08
     &       + u08*x1*y10*p10 + u10*x2*y12*p12 + u12*x3*y14*p14
     &       + u14*x4*y16*p16 + u16*x5*y18*p18

        e09 = 2.0d0
        e11 =-105.0d0*x1*u02
        e13 = 4536.0d0*x2*u04
        e15 =-194040.0d0*x3*u06
        e17 = 8648640.0d0*x4*u08
        e19 =-408648240.0d0*x5*u10
        e21 = 20583763200.0d0*x6*u12
        e23 =-1106230245120.0d0*x7*u14
        
        ye  = e09*s09 + e11*s11 + e13*s13 + e15*s15
     &      + e17*s17 + e19*s19 + e21*s21 + e23*s23

        Y34  = 3.0D0*U*U*C2L*C2L*(yp + 480.0d0*u06*u*ye)
c
      ELSE
        Y34 = 0.0D0
      endif

C      
      FMIINT0 = FMIINT0 + FNRM*Y34

C
C     NOTHING   -51891.322203
C     QUADRTR   -51891.186189
C     GAUSIAN   -51891.361960
C
      RETURN
      END
C
C FERMI   A+B      34864.988757
C FERMI W/C+D  0.5 34864.979133
C FERMI W/C+D  0.1 34864.979133
C FERMI 20000 15.0 34864.863362
C GAUSSMN          34864.857069

