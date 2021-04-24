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
C  THIS ROUTINE USES MY OWN RECIPE FOR THE ANNOYING INTEGRAL CLASSES,  C
C  FOR THE FIRST FEW L CASES AND ON THE CONDITION THAT A√λ < 0.05D0.   C
C  DFNOTE: I'VE RETAINED TERMS FROM THE POLYNOMIAL EXPRESSIONS HERE,   C
C          WHICH DON'T ACTUALLY CONTRIBUTE TO NUMERICAL ACCURACY IN A  C
C          MEASURABLE WAY. THIS ROUTINE ALSO RETAINS THE S_K(X) TERMS. C
C**********************************************************************C
      INCLUDE 'param.h'
C
      COMMON/GAMA/GAMLOG(300),GAMHLF(300)
      COMMON/MATH/PI,PI12,PI32,PI52,PILG,TWLG,THLG,TW12
C
      DIMENSION U(0:17),X(0:9),P(20),S(21)
C
C     ROUTINE ONLY ALLOWS NON-NEGATIVE INTEGER L
      IF(L.LT.0) THEN
        WRITE(6, *) 'In FMIINT0: invalid parameter L. L = ',L
        WRITE(7, *) 'In FMIINT0: invalid parameter L. L = ',L
        STOP
      ENDIF
C
C     FACTORS NEEDED FOR ALL PARAMETERS N
      U(0) = 1.0D0
      DO K=1,17
        U(K) = U(K-1)*A/C
      ENDDO
      
      X(0) = 1.0D0
      DO K=1,9
        X(K) = X(K-1)*EIJ*C*C
      ENDDO

      DO K=2,20,2
        P(K) = POLYLOG(K,0.0D0)
      ENDDO
      
      P( 2) =-(PI**2)/12.0D0
      P( 4) =-(PI**4)*7.0D0/720.0D0
      P( 6) =-(PI**6)*31.0D0/30240.0D0
      P( 8) =-(PI**8)*127.0D0/12096.0D2
      P(10) =-(PI**10)*73.0D0/6842880.0D0
      P(12) =-(PI**12)*1414477.0D0/13076743688.0D2
      P(14) =-(PI**14)*8191.0D0/747242496.0D2
      P(16) =-(PI**16)*16931177.0D0/152437469184.0D4
      P(18) =-(PI**18)*5749691557.0D0/5109094217170944.0D3
      P(20) =-(PI**20)*91546277357.0D0/8028576626982912.0D5
      
      DO K=3,25,2
        S(K) = POLYLOG(K,-1.0D0/U(1))
      ENDDO
C
C     FERMI NORMALISATION CONSTANT
      FNRM = 1.0D0 + PI*PI*U(2) - 6.0D0*U(3)*S(3)
      FNRM = 1.0D0/FNRM
C
C     INTEGRAL TYPE: INCOMPLETE GAMMA FUNCTION
      X1A = GAMLWR(2*L+3,X(1))*(3.0D0+PI*PI*U(2))/(2.0D0*DSQRT(X(1)))
      X1B =-GAMLWR(2*L+5,X(1))/(2.0D0*X(1)*DSQRT(X(1)))
      X1C = GAMUPR(2*L+2,X(1))*(1.0D0+PI*PI*U(2))
C
C     INTEGRAL TYPE: LOGARITHMIC INTEGRAL
      X2A =-6.0D0*GAMHLF(2*L)*U(3)*S(3)
C
C     TRANSFER DATA TO FMIINT0
      X1 = X1A+X1B+X1C
      X2 = X2A
      FMIINT0 = 0.5D0*FNRM*(X1+X2)/(EIJ**(L+1))
C
      TOL = 0.05D0
C
C     ALGORITHMS AVAILABLE FOR THE REMAINING SET OF INTEGRALS
      IF(L.EQ.0.AND.A*DSQRT(EIJ).LT.TOL) THEN
C
        Y02 = 1.0D0/36288.0D1 - X(1)/3628.0D2
        Y04 = 8.0D0 - 12.0D0*X(1) + 8.0D0*X(2) - 1.0D1*X(3)/3.0D0
     &      + X(4) - 7.0D0*X(5)/3.0D1 + 2.0D0*X(6)/45.0D0 
     &      - X(7)/14.0D1 + X(8)/2016.0D1
        Y06 =-72.0D0 + 16.0D1*X(1) - 14.0D1*X(2) + 72.0D0*X(3) 
     &      - 77.0D0*X(4)/3.0D0 + 104.0D0*X(5)/15.0D0 - 1.5D0*X(6)
     &      + 17.0D0*X(7)/126.0D1
        Y08 = 96.0D1 - 28.0D2*X(1) + 3024.0D0*X(2) - 1848.0D0*X(3) 
     &      + 2288.0D0*X(4)/3.0D0 - 234.0D0*X(5) + 17.0D0*X(6)/6.0D0
        Y10 =-168.0D2 + 6048.0D1*X(1) - 77616.0D0*X(2)+ 54912.0D0*X(3)
     &      - 2574.0D1*X(4) + 442.0D0*X(5)
        Y12 = 36288.0D1 - 155232.0D1*X(1) + 2306304.0D0*X(2) 
     &      - 185328.0D1*X(3) + 4862.0D1*X(4)
        Y14 =-931392.0D1 + 4612608.0D1*X(1) - 7783776.0D1*X(2)
     &      + 350064.0D1*X(3)
        Y16 = 27675648.0D1 - 15567552.0D2*X(1) + 14702688.0D1*X(2)
        Y18 =-93405312.0D2 + 29405376.0D2*X(1)
        Y20 = 176432256.0D2
C        
        YP   = Y02*P( 2)*X(8)/U( 1) + Y04*P( 4)*X(0)*U( 1)
     &       + Y06*P( 6)*X(1)*U( 3) + Y08*P( 8)*X(2)*U( 5)
     &       + Y10*P(10)*X(3)*U( 7) + Y12*P(12)*X(4)*U( 9)
     &       + Y14*P(14)*X(5)*U(11) + Y16*P(16)*X(6)*U(13)
     &       + Y18*P(18)*X(7)*U(15) + Y20*P(20)*X(8)*U(17)
C
        E05 = 1.0D0
        E07 =-9.0D0*X(1)*U(2)
        E09 = 12.0D1*X(2)*U(4)
        E11 =-21.0D2*X(3)*U(6)
        E13 = 4536.0D1*X(4)*U(8)
        E15 =-116424.0D1*X(5)*U(10)
        E17 = 3459456.0D1*X(6)*U(12)
        E19 =-11675664.0D2*X(7)*U(14)
        E21 = 22054032.0D2*X(8)*U(16)
C    
        YE  = E05*S( 5) + E07*S( 7) + E09*S( 9) + E11*S(11)
     &      + E13*S(13) + E15*S(15) + E17*S(17) + E19*S(19) + E21*S(21)
C
        Y34  = 3.0D0*U(3)*C*C*(YP + 4.0D0*U(2)*YE)
C
      ELSEIF(L.EQ.1.AND.A*DSQRT(EIJ).LT.TOL) THEN
C
        Y02 =-19.0D0/8064.0D1 + 17.0D0*X(1)/72576.0D1
     &      - 19.0D0*X(2)/72576.0D2 + X(3)/36288.0D2
        Y04 = 12.0D0 - 16.0D0*X(1) + 1.0D1*X(2) - 4.0D0*X(3)
     &      + 7.0D0*X(4)/6.0D0 - 4.0D0*X(5)/15.0D0 - 71.0D0*X(6)/504.0D1
     &      + X(7)/1344.0D1
        Y06 = 72.0D0 - 32.0D1*X(1) + 42.0D1*X(2) - 288.0D0*X(3) 
     &      + 385.0D0*X(4)/3.0D0 - 41.6D0*X(5) - 71.0D0*X(6)/24.0D0
     &      + 17.0D0*X(7)/84.0D1
        Y08 =-192.0D1 + 84.0D2*X(1) - 12096.0D0*X(2) + 924.0D1*X(3) 
     &      - 4576.0D0*X(4) - 461.5D0*X(5) + 4.25D0*X(6)
        Y10 = 504.0D2 - 24192.0D1*X(1) + 38808.0D1*X(2)
     &      - 329472.0D0*X(3) - 50765.0D0*X(4) + 663.0D0*X(5)
        Y12 =-145152.0D1 + 77616.0D2*X(1) - 13837824.0D0*X(2) 
     &      - 365508.0D1*X(3) + 7293.0D1*X(4)
        Y14 = 465696.0D2 - 27675648.0D1*X(1) - 15351336.0D1*X(2)
     &      + 525096.0D1*X(3)
        Y16 =-166053888.0D1 - 30702672.0D2*X(1) + 22054032.0D1*X(2)
        Y18 =-184216032.0D2 + 44108064.0D2*X(1)
        Y20 = 264648384.0D2
C
        YP   = Y02*P( 2)*X(6)/U( 1) + Y04*P( 4)*X(0)*U( 1)
     &       + Y06*P( 6)*X(0)*U( 3) + Y08*P( 8)*X(1)*U( 5)
     &       + Y10*P(10)*X(2)*U( 7) + Y12*P(12)*X(3)*U( 9)
     &       + Y14*P(14)*X(4)*U(11) + Y16*P(16)*X(5)*U(13)
     &       + Y18*P(18)*X(6)*U(15) + Y20*P(20)*X(7)*U(17)
C
        E07 = 1.0D0
        E09 =-8.0D1*X(1)*U(2)/3.0D0
        E11 = 7.0D2*X(2)*U(4)
        E13 =-2016.0D1*X(3)*U(6)
        E15 = 6468.0D2*X(4)*U(8)
        E17 =-2306304.0D1*X(5)*U(10)
        E19 =-2558556.0D2*X(6)*U(12)
        E21 = 3675672.0D2*X(7)*U(14)
C   
        YE  = E07*S( 7) + E09*S( 9) + E11*S(11) + E13*S(13)
     &      + E15*S(15) + E17*S(17) + E19*S(19) + E21*S(21)
C
        Y34 = 3.0D0*U(3)*(C**4)*(YP + 36.0D0*U(4)*YE)
C
      ELSEIF(L.EQ.2.AND.A*DSQRT(EIJ).LT.TOL) THEN
C
        Y02 =-323.0D0/2016.0D1 + 19.0D0*X(1)/2016.0D1
     &      - 47.0D0*X(2)/290304.0D0 + 43.0D0*X(3)/20736.0D2
     &      - 17.0D0*X(4)/72576.0D2 + X(5)/36288.0D2
        Y04 = 16.0D0 - 2.0D1*X(1) + 12.0D0*X(2) - 14.0D0*X(3)/3.0D0
     &      - 2.03125D0*X(4) - 11.0D0*X(5)/252.0D0 + X(6)/5376.0D0
        Y06 = 32.0D1 - 84.0D1*X(1) + 864.0D0*X(2) - 154.0D1*X(3)/3.0D0
     &      - 316.875D0*X(4) - 55.0D0*X(5)/6.0D0 + 17.0D0*X(6)/336.0D0
        Y08 = 192.0D1 - 168.0D2*X(1) + 36288.0D0*X(2) - 3696.0D1*X(3)
     &      - 34856.25D0*X(4) - 143.0D1*X(5) + 10.625D0*X(6)
        Y10 =-1008.0D2 + 72576.0D1*X(1) - 155232.0D1*X(2)
     &      - 250965.0D1*X(3) - 1573.0D2*X(4) + 1657.5D0*X(5)
        Y12 = 435456.0D1 - 310464.0D2*X(1) - 1054053.0D2*X(2) 
     &      - 113256.0D2*X(3) + 182325.0D0*X(4)
        Y14 =-1862784.0D2 - 2108106.0D3*X(1) - 4756752.0D2*X(2)
     &      + 131274.0D2*X(3)
        Y16 =-12648636.0D3 - 9513504.0D3*X(1) + 5513508.0D2*X(2)
        Y18 =-57081024.0D3 + 11027016.0D3*X(1)
        Y20 = 66162096.0D3
C
        YP   = Y02*P( 2)*X(4)/U( 1) + Y04*P( 4)*X(0)*U( 1)
     &       + Y06*P( 6)*X(0)*U( 3) + Y08*P( 8)*X(0)*U( 5)
     &       + Y10*P(10)*X(1)*U( 7) + Y12*P(12)*X(2)*U( 9)
     &       + Y14*P(14)*X(3)*U(11) + Y16*P(16)*X(4)*U(13)
     &       + Y18*P(18)*X(5)*U(15) + Y20*P(20)*X(6)*U(17)
C
        E09 = 1.0D0
        E11 =-52.5D0*X(1)*U(2)
        E13 = 2268.0D0*X(2)*U(4)
        E15 =-9702.0D1*X(3)*U(6)
        E17 =-6587831.25D0*X(4)*U(8)
        E19 =-297297.0D2*X(5)*U(10)
        E21 = 34459425.0D0*X(6)*U(12)
C
        YE  = E09*S( 9) + E11*S(11) + E13*S(13)
     &      + E15*S(15) + E17*S(17) + E19*S(19) + E21*S(21)
C
        Y34 = 3.0D0*U(3)*(C**6)*(YP + 96.0D1*U(6)*YE)
C
      ELSEIF(L.EQ.3.AND.A*DSQRT(EIJ).LT.TOL) THEN
C
        Y02 =-323.0D0/384.0D0 - 323.0D0*X(1)/8064.0D0
     &      - 19.0D0*X(2)/3072.0D0 + X(3)/1024.0D0
     &      - 4229.0D0*X(4)/290304.0D2 + 13.0D0*X(5)/6912.0D2
     &      - X(6)/48384.0D1 + X(7)/36288.0D2
        Y04 = 2.0D1 - 24.0D0*X(1) - 117.21875D0*X(2)
     &      - 2639.0D0*X(3)/192.0D0 - 35.0D0*X(4)/192.0D0
     &      + X(5)/1536.0D0
        Y06 = 84.0D1 - 1728.0D0*X(1) - 12894.0625D0*X(2)
     &      - 2144.1875D0*X(3) - 38.28125D0*X(4) + 17.0D0*X(5)/96.0D0
        Y08 = 168.0D2 - 72576.0D0*X(1) - 928372.5D0*X(2)
     &      - 235860.625D0*X(3) - 5971.875D0*X(4) + 37.1875D0*X(5)
        Y10 = 1008.0D2 - 145152.0D1*X(1) - 38991645.0D0*X(2)
     &      - 16981965.0D0*X(3) - 656906.25D0*X(4) + 4801.25D0*X(5)
        Y12 =-870912.0D1 - 7798329.0D2*X(1) - 71324253.0D1*X(2) 
     &      - 4729725.0D1*X(3) + 638137.5D0*X(4)
        Y14 =-46789974.0D2 - 142648506.0D2*X(1) - 19864845.0D2*X(2)
     &      + 459459.0D2*X(3)
        Y16 =-855891036.0D2 - 3972969.0D4*X(1) + 19297278.0D2*X(2)
        Y18 =-23837814.0D4 + 38594556.0D3*X(1)
        Y20 = 231567336.0D3
C
        YP   = Y02*P( 2)*X(2)/U( 1) + Y04*P( 4)*X(0)*U( 1)
     &       + Y06*P( 6)*X(0)*U( 3) + Y08*P( 8)*X(0)*U( 5)
     &       + Y10*P(10)*X(0)*U( 7) + Y12*P(12)*X(1)*U( 9)
     &       + Y14*P(14)*X(2)*U(11) + Y16*P(16)*X(3)*U(13)
     &       + Y18*P(18)*X(4)*U(15) + Y20*P(20)*X(5)*U(17)
C
        E11 = 1.0D0
        E13 =-86.4D0*X(1)*U( 2)
        E15 =-46418.625D0*X(2)*U( 4)
        E17 =-849098.25D0*X(3)*U( 6)
        E19 =-2364862.5D0*X(4)*U( 8)
        E21 = 2297295.0D0*X(5)*U(10)
C
        YE  = E11*S(11) + E13*S(13) + E15*S(15)
     &      + E17*S(17) + E19*S(19) + E21*S(21)
C
        Y34 = 3.0D0*U(3)*(C**8)*(YP + 504.0D2*U(8)*YE)
C
      ELSEIF(L.EQ.4.AND.A*DSQRT(EIJ).LT.TOL) THEN
C
        Y02 =-32.8046875D0 - 3553.0D0*X(1)/384.0D0
     &      - 7429.0D0*X(2)/1344.0D1 + 0.018554687D0*X(3)
     &      - 269.0D0*X(4)/55296.0D0 + 5627.0D0*X(5)/64512.0D2
     &      - 1277.0D0*X(6)/96768.0D2 + 253.0D0*X(7)/145152.0D2
     &      - 13.0D0*X(8)/72576.0D2 + X(9)/36288.0D2
        Y04 =-3584.515625D0 - 1471.40625D0*X(1) - 100.078125D0*X(2)
     &      - 0.953125D0*X(3) + 0.002929688D0*X(4)
        Y06 =-258085.125D0 - 161854.6875D0*X(1) - 15612.1875D0*X(2)
     &      - 200.15625D0*X(3) + 0.796875D0*X(4)
        Y08 =-10839575.25D0 - 11653537.5D0*X(1) - 1717340.625D0*X(2)
     &      - 31224.375D0*X(3) + 167.34375D0*X(4)
        Y10 =-216791505.0D0 - 489448575.0D0*X(1) - 123648525.0D0*X(2)
     &      - 3434681.25D0*X(3) + 26105.625D0*X(4)
        Y12 =-130074903.0D1 - 97889715.0D2*X(1) - 519323805.0D1*X(2) 
     &      - 24729705.0D1*X(3) + 2871618.75D0*X(4)
        Y14 =-58733829.0D3 - 103864761.0D3*X(1) - 103864761.0D2*X(2)
     &      + 20675655.0D1*X(3)
        Y16 =-623188566.0D3 - 207729522.0D3*X(1) + 86837751.0D2*X(2)
        Y18 =-1246377132.0D3 + 173675502.0D3*X(1)
        Y20 = 1042053012.0D3
C
        YP   = Y02*P( 2)*X(0)/U( 1) + Y04*P( 4)*X(0)*U( 1)
     &       + Y06*P( 6)*X(0)*U( 3) + Y08*P( 8)*X(0)*U( 5)
     &       + Y10*P(10)*X(0)*U( 7) + Y12*P(12)*X(0)*U( 9)
     &       + Y14*P(14)*X(1)*U(11) + Y16*P(16)*X(2)*U(13)
     &       + Y18*P(18)*X(3)*U(15) + Y20*P(20)*X(4)*U(17)
C
        E13 =-1.0D0
        E15 =-103587.0D2*X(1)*U(2)/229409.0D0
        E17 =-1099098.0D2*X(2)*U(4)/229409.0D0
        E19 =-2198196.0D2*X(3)*U(6)/229409.0D0
        E21 = 1837836.0D2*X(4)*U(8)/229409.0D0
C
        YE  = E13*S(13) + E15*S(15)
     &      + E17*S(17) + E19*S(19) + E21*S(21)
C
        Y34 = 3.0D0*U(3)*(C**10)*(YP + 650374515.0D0*U(10)*YE)
C
      ELSE
C
C       INTEGRATION GRID DETAILS
        NLIN = 100
        NEXP = 1000
        RMAX = 10.0D0
        ELIN = 1.0D0/DFLOAT(NLIN)
        EEXP = DLOG(RMAX/C)/DFLOAT(NEXP)
C
C       NUMERICALLY INTEGRATE REMAINING TERMS FROM 0 TO C
        X0C = 0.0D0
        DO N=0,NLIN
          ZN  = ELIN*DFLOAT(N)
          Y1  = ZN**(2*L+1)
          Y2  = DEXP(-X(1)*ZN*ZN)
          Y3  =-ZN*POLYLOG(2,(ZN-1.0D0)/U(1))
          Y4  = 2.0D0*U(1)*POLYLOG(3,(ZN-1.0D0)/U(1))
          X0C = X0C + EXTINT11(Y1*Y2*(Y3+Y4),N,NLIN)
        ENDDO
        X0C = 15.0D0*U(2)*(C**(2*L+2))*ELIN*X0C/299376.0D0
C
C       NUMERICALLY INTEGRATE REMAINING TERMS FROM C TO INFINITY (EXP.)
        XCI = 0.0D0
        DO N=0,NEXP
          TN  = EEXP*DFLOAT(N)
          ZN  = DEXP(TN)
          Y1  = ZN**(2*L+2)
          Y2  = DEXP(-X(1)*ZN*ZN)
          Y3  = ZN*POLYLOG(2,(1.0D0-ZN)/U(1))
          Y4  = 2.0D0*U(1)*POLYLOG(3,(1.0D0-ZN)/U(1))
          XCI = XCI + EXTINT11(Y1*Y2*(Y3+Y4),N,NEXP)
        ENDDO
        XCI = 15.0D0*U(2)*(C**(2*L+2))*EEXP*XCI/299376.0D0
        
        Y34 = X0C+XCI
C
      ENDIF
C      
      FMIINT0 = FMIINT0 + FNRM*Y34
C
      RETURN
      END

