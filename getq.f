 
C----------------------------------------------------------------------
C FUNCTION GETQ
C
C Called by: LAISUM
C
C Purpose: to retrieve the collision functions from the QFUNC common block,
C that were set using the subroutine SETQ.
C
C Uses global variables:
C      QQQ    - the collision functions Q
C
C Formal parameters:
C      Ind    - Index of Q function to store
C
      REAL*8 FUNCTION GETQ(Ind)

      IMPLICIT NONE
      INTEGER*4 Ind
      REAL*8 QQQ
      COMMON /QFUNC / QQQ(4900)

C     Check the index is valid

      IF ( Ind.LE.0 .OR. Ind.GE. 4900) THEN
         WRITE(*,*) 'Invalid Q function number ', Ind
         STOP 'Bug!'
      ENDIF

      GETQ = QQQ(Ind)

      RETURN
      END
