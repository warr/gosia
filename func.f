 
C----------------------------------------------------------------------
C FUNCTION FUNC
C
C Called by: LAGRAN
C
C Purpose: evaluates f(y) = y, f(y) = log_e(y) or f(y) = sqrt(y), depending
C on the flag I
      
      REAL*8 FUNCTION FUNC(Y,I)
      IMPLICIT NONE
      INTEGER*4 I
      REAL*8 Y
      
      IF ( I.EQ.2 ) THEN
         IF ( Y.LT.1.E-12 ) Y = 1.E-12
         FUNC = LOG(Y)
         RETURN
      ELSEIF ( I.EQ.3 ) THEN
         FUNC = SQRT(Y)
         GOTO 99999
      ENDIF
      FUNC = Y
99999 END
