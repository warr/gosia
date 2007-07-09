 
C----------------------------------------------------------------------
 
      SUBROUTINE SETIN
      IMPLICIT NONE
      REAL*8 ADB , CH , EPS , EROot , FIEx , SH
      INTEGER*4 IAXs , IEXp , k , LP1 , LP10 , LP11 , LP12 , LP13 , 
     &          LP14 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /HIPER / SH(365) , CH(365)
      COMMON /ADX   / ADB(365)
      COMMON /KIN   / EPS(50) , EROot(50) , FIEx(50,2) , IEXp , IAXs(50)
      DO k = 1 , LP12
         ADB(k) = EPS(IEXp)*SH(k) + .03*(k-1)
      ENDDO
      END
