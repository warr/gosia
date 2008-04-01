 
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
C      Iomega - Index of omega value for which we want the Q value
C      Lambda - Lambda for which we want  the Q value
C      Mu     - Mu for which we want the Q value
C      Value  - the value of the Q function we want to store
C
C The index is Iomega + LOCQ(lambda, mu), where icnt is an index to the
C appropriate value of omega. The steps of omega are in 0.03.
      
      SUBROUTINE SETQ(Iomega, Lambda, Mu, Value)

      IMPLICIT NONE
      INTEGER*4 ind, Iomega, Lambda, Mu, LOCQ
      REAL*8 QQQ, Value
      COMMON /ALLC  / LOCQ(8,7)
      COMMON /QFUNC / QQQ(4900)

C     Check lambda is valid

      IF ( Lambda.LT.1 .OR. Lambda.GT.8 ) THEN
         WRITE(*,*) 'Invalid value for lambda ', Lambda
         STOP 'Bug!'
      ENDIF

C     Check mu is valid

      IF ( Mu.LT.1 .OR. Mu.GT.Lambda+1 ) THEN 
         WRITE(*,*) 'Invalid value for lamda = ', Lambda, ' mu =', Mu
         STOP 'Bug!'
      ENDIF
      
C     Calculate index

      ind = LOCQ(Lambda, Mu) + Iomega

C     Check the index is valid

      IF ( ind.LE.0 .OR. ind.GE. 4900 ) THEN
         WRITE(*,*) 'Invalid Q function number ', ind
         STOP 'Bug!'
      ENDIF

C     Set the value
      
      QQQ(ind) = Value

      RETURN
      END
