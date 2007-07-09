 
C----------------------------------------------------------------------
 
      SUBROUTINE NEWCAT(Iexp,Jidim)
      IMPLICIT NONE
      REAL*8 a , b , FXIS1 , FXIS2 , PARx , PARxm , q1 , q2 , QAPr , 
     &       wg , wl , XI , XIR , xp , xx , zt
      INTEGER*4 IAPr , Iexp , IPAth , ISEx , ist , istop , Jidim , k , 
     &          kk , LAMda , LAMmax , LDNum , LEAd , MAGa , MULti , n , 
     &          NDIm , ng , nl , NMAx
      INTEGER*4 NMAx1
      COMMON /MAP   / PARx(50,12,5) , PARxm(50,4,10,6) , XIR(6,50)
      COMMON /CXI   / XI(500)
      COMMON /CLCOM / LAMda(8) , LEAd(2,500) , LDNum(8,75) , LAMmax , 
     &                MULti(8)
      COMMON /COEX2 / NMAx , NDIm , NMAx1
      COMMON /PTH   / IPAth(75) , MAGa(75)
      COMMON /APRCAT/ QAPr(500,2,7) , IAPr(500,2) , ISEx(75)
      Jidim = NMAx + 1
      IF ( MAGa(Iexp).NE.0 ) Jidim = 3*NMAx + 1
      ist = 1
      DO kk = 1 , 6
         IF ( MULti(kk).NE.0 ) THEN
            istop = MULti(kk) - 1 + ist
            DO k = ist , istop
               xx = ABS(XI(k))
               xx = xx/XIR(kk,Iexp)
               DO n = 1 , 7 , 3
                  IF ( MAGa(Iexp).NE.0 .OR. n.EQ.4 ) THEN
                     zt = QAPr(k,1,n)
                     zt = ABS(zt)
                     xp = 9.*xx
                     nl = INT(xp) + 1
                     wg = xp - DBLE(nl-1)
                     ng = nl + 1
                     wl = DBLE(nl) - xp
                     a = wg*PARxm(Iexp,1,ng,kk) + wl*PARxm(Iexp,1,nl,kk)
                     b = wg*PARxm(Iexp,2,ng,kk) + wl*PARxm(Iexp,2,nl,kk)
                     q1 = a*zt + b
                     a = wg*PARxm(Iexp,3,ng,kk) + wl*PARxm(Iexp,3,nl,kk)
                     b = wg*PARxm(Iexp,4,ng,kk) + wl*PARxm(Iexp,4,nl,kk)
                     q2 = a*zt + b
                     QAPr(k,2,n) = QAPr(k,1,n)*q2*FXIS2(k,n)
                     QAPr(k,1,n) = QAPr(k,1,n)*q1*FXIS1(k,n)
                     IF ( IAPr(k,1).EQ.IAPr(k,2) ) THEN
                        QAPr(k,1,n) = 0.
                        QAPr(k,2,n) = QAPr(k,2,n)/2.
                     ENDIF
                  ENDIF
               ENDDO
               IF ( MAGa(Iexp).NE.0 ) THEN
                  DO n = 2 , 6
                     IF ( n.NE.4 ) THEN
                        zt = QAPr(k,1,n)
                        zt = ABS(zt)
                        xp = 4.*xx
                        nl = INT(xp) + 1
                        wg = xp - DBLE(nl-1)
                        ng = nl + 1
                        wl = DBLE(nl) - xp
                        q1 = wg*PARx(Iexp,2*kk-1,ng)
     &                       + wl*PARx(Iexp,2*kk-1,nl)
                        q2 = wg*PARx(Iexp,2*kk,ng)
     &                       + wl*PARx(Iexp,2*kk,nl)
                        QAPr(k,2,n) = QAPr(k,1,n)*q2*FXIS2(k,n)
                        QAPr(k,1,n) = QAPr(k,1,n)*q1*FXIS1(k,n)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            ist = istop + 1
         ENDIF
      ENDDO
      END
