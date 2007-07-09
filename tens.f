 
C----------------------------------------------------------------------
 
      SUBROUTINE TENS(Bten)
      IMPLICIT NONE
      REAL*8 arg , Bten , DJMM , DSIgs , EPS , EROot , FIEx , TETacm , 
     &       TREp , ZETa
      INTEGER*4 i , IAXs , IEXp , ind , inz , iph , ix , k , k1 , kp , 
     &          l , lp , lpp , lx , lxx , LZEta , NDIm , NMAx , NMAx1
      DIMENSION Bten(1200)
      COMMON /KIN   / EPS(50) , EROot(50) , FIEx(50,2) , IEXp , IAXs(50)
      COMMON /TCM   / TETacm(50) , TREp(50) , DSIgs(50)
      COMMON /CCOUP / ZETa(50000) , LZEta(8)
      COMMON /COEX2 / NMAx , NDIm , NMAx1
      ix = NMAx*28
      arg = 1.570796327 + TETacm(IEXp)/2.
      DO i = 1 , ix
         ZETa(i) = 0.
      ENDDO
      DO i = 2 , NMAx
         DO kp = 1 , 7 , 2
            k = kp - 1
            k1 = INT(DBLE(k)/2.+.01)
            IF ( k.EQ.0 ) THEN
               ind = (i-2)*16 + 1
               inz = (i-1)*28 + 1
               ZETa(inz) = Bten(ind)
            ELSE
               DO lp = 1 , kp
                  IF ( IAXs(IEXp).NE.0 .OR. lp.EQ.1 ) THEN
                     inz = (i-1)*28 + k1*7 + lp
                     l = lp - 1
                     DO lpp = 1 , kp
                        ind = k*k/4 + lpp + (i-2)*16
                        lx = lpp - 1
                        lxx = lx
 2                      iph = (-1)**(l+INT(DBLE(lxx)/2.))
                        ZETa(inz) = ZETa(inz) + Bten(ind)
     &                              *iph*DJMM(arg,k,lx,l)
                        IF ( lpp.NE.1 ) THEN
                           IF ( lx.GE.0 ) THEN
                              lx = -lx
                              lxx = lx - 1
                              GOTO 2
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      END
