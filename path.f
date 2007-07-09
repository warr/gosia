 
C----------------------------------------------------------------------
 
      SUBROUTINE PATH(Irld)
      IMPLICIT NONE
      REAL*8 CAT , spm , vl
      INTEGER*4 i , IPAth , Irld , ISMax , isp , ist , j , MAGa , NDIm , 
     &          NMAx , NMAx1 , NSTart , NSTop
      COMMON /CEXC0 / NSTart(76) , NSTop(75)
      COMMON /PTH   / IPAth(75) , MAGa(75)
      COMMON /COEX2 / NMAx , NDIm , NMAx1
      COMMON /CLCOM8/ CAT(600,3) , ISMax
      spm = CAT(Irld,3)
      DO i = 2 , NMAx
         IPAth(i) = 0
         ist = NSTart(i)
         IF ( ist.NE.0 ) THEN
            isp = NSTop(i)
            DO j = ist , isp
               vl = CAT(j,3)
               IF ( ABS(vl-spm).LT.1.E-6 ) GOTO 50
            ENDDO
         ENDIF
         GOTO 100
 50      IPAth(i) = j
 100  ENDDO
      IPAth(1) = Irld
      END
