 
C----------------------------------------------------------------------
 
      SUBROUTINE PODZIEL(I,J)
      IMPLICIT NONE
      INTEGER*4 I , IAPr , IDIve , ISEx , J , k , l , l1 , l2 , LERf , 
     &          LP1 , LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , LP3 , 
     &          LP4 , LP6
      INTEGER*4 LP7 , LP8 , LP9
      REAL*8 QAPr
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /APRCAT/ QAPr(500,2,7) , IAPr(500,2) , ISEx(75)
      COMMON /APRX  / LERf , IDIve(50,2)
      IF ( I.NE.3 ) THEN
         IF ( I.EQ.1 ) THEN
            l1 = IDIve(J,1)
            IDIve(J,1) = l1 + 1
            GOTO 100
         ELSE
            l1 = IDIve(J,2)
            IDIve(J,2) = l1 + 1
         ENDIF
      ENDIF
      l2 = IDIve(J,2)
      IF ( I.EQ.3 ) l1 = 1
      DO k = 1 , LP2
         DO l = 1 , 7
            QAPr(k,2,l) = QAPr(k,2,l)*l1/l2
         ENDDO
      ENDDO
      IF ( I.EQ.2 ) WRITE (22,99001) J , IDIve(J,1) , l2
      IF ( I.NE.3 ) RETURN
 100  l2 = IDIve(J,1)
      IF ( I.EQ.3 ) l1 = 1
      DO k = 1 , LP2
         DO l = 1 , 7
            QAPr(k,1,l) = QAPr(k,1,l)*l1/l2
         ENDDO
      ENDDO
      IF ( I.EQ.1 ) WRITE (22,99001) J , l2 , IDIve(J,2)
      RETURN
99001 FORMAT (5X,'*****',1X,'EXP(A) EXPANSION FAILURE!',1X,'*****'/5X,
     &        'EXPERIMENT',1X,1I2,3X,'NEW SUBDIVISION',1X,'(',1I1,',',
     &        1I1,')')
      END
