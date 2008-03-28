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
