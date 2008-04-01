 
C----------------------------------------------------------------------
C FUNCTION GETQ
C
C Called by: LAISUM
C
C Purpose: to retrieve the collision functions from the QFUNC common block,
C that were set using the subroutine SETQ.
C
C Uses global variables:
C      LOCQ   - the location of the Q function for each lambda, mu combination
C      QQQ    - the collision functions Q
C
C Formal parameters:
C      Iomega - Index of omega value for which we want the Q value
C      Lambda - Lambda for which we want  the Q value
C      Mu     - Mu for which we want the Q value
C
C Return value:
C      Q-function value
C
C The index is Iomega + LOCQ(lambda, mu), where icnt is an index to the
C appropriate value of omega. The steps of omega are in 0.03.

      REAL*8 FUNCTION GETQ(Iomega, Lambda, Mu)

      IMPLICIT NONE
      INTEGER*4 ind, Iomega, Lambda, Mu, LOCQ
      REAL*8 QQQ
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

C     Get the value
      
      GETQ = QQQ(ind)

      RETURN
      END
