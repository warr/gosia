 
C----------------------------------------------------------------------
 
      SUBROUTINE HALF(Iso)
      IMPLICIT NONE
      REAL*8 CAT
      INTEGER*4 ir , ISMax , Iso , j , NDIm , NMAx , NMAx1 , NSTart , 
     &          NSTop
      COMPLEX*16 ARM , fpom
      COMMON /CEXC0 / NSTart(76) , NSTop(75)
      COMMON /AZ    / ARM(600,7)
      COMMON /CLCOM8/ CAT(600,3) , ISMax
      COMMON /COEX2 / NMAx , NDIm , NMAx1
      IF ( Iso.EQ.0 ) THEN
         DO j = 1 , NMAx
            ir = NSTart(j) - 1
 20         ir = ir + 1
            fpom = ARM(ir,3)
            ARM(ir,1) = -.0625*(ARM(ir,1)+ARM(ir,4))
     &                  + .5625*(ARM(ir,2)+ARM(ir,3))
            ARM(ir,3) = ARM(ir,3)*.75 + .375*ARM(ir,4) - ARM(ir,2)/8.
            ARM(ir,2) = fpom
            IF ( CAT(ir,3).LT.-.1 ) GOTO 20
         ENDDO
         GOTO 99999
      ENDIF
      DO j = 1 , ISMax
         fpom = ARM(j,3)
         ARM(j,1) = -.0625*(ARM(j,4)+ARM(j,1))
     &              + .5625*(ARM(j,2)+ARM(j,3))
         ARM(j,3) = ARM(j,3)*.75 + .375*ARM(j,4) - ARM(j,2)/8.
         ARM(j,2) = fpom
      ENDDO
99999 END
