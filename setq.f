 
C----------------------------------------------------------------------
C SUBROUTINE SETQ
C
C Called by: SNAKE
C
C Purpose: to store the collision functions in the QFUNC common block, so
C that they can be retrieved by the function GETQ.
C
C Uses global variables:
C      QQQ    - the collision functions Q
C
C Formal parameters:
C      Ind    - Index of Q function to store
C
C Return value:
C      Q-function value
C
      SUBROUTINE SETQ(Ind, Value)

      IMPLICIT NONE
      INTEGER*4 Ind
      REAL*8 QQQ, Value
      COMMON /QFUNC / QQQ(4900)
      
C     Check the index is valid

      IF ( Ind.LE.0 .OR. Ind.GE. 4900) THEN
         WRITE(*,*) 'Invalid Q function number ', Ind
         STOP 'Bug!'
      ENDIF

      QQQ(Ind) = Value

      RETURN
      END
