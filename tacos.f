 
C----------------------------------------------------------------------
C FUNCTION TACOS
C
C Called by: ARCCOS, CEGRY, COORD, GOSIA2
C Calls:     TASIN
C
C Purpose: evaluate arccosine(x)
C
C We use: arccos(x) = pi/2 - arcsin(x)
 
      REAL*8 FUNCTION TACOS(X)
      IMPLICIT NONE
      REAL*8 TASIN , X
      
      TACOS = 1.570796327 - TASIN(X) ! 1.570796327 = pi / 2
      END
