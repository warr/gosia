 
C----------------------------------------------------------------------
 
      SUBROUTINE LAIAMP(Ir,W0)
      IMPLICIT NONE
      REAL*8 CAT , D2W , ELM , ELMl , ELMu , EPS , epsi , EROot , errt , 
     &       FIEx , pm , ppp , rmir , rmis , rmu , SA , TCABS , W0 , 
     &       XI , xiv
      REAL*8 z , ZETa
      INTEGER*4 i1 , i2 , i3 , IAXs , IEXp , indx , Ir , is , is1 , 
     &          is2 , ISG , ISG1 , ISMax , ismin , isplus , KDIv , la , 
     &          lam , LAMda , LAMmax
      INTEGER*4 LAMr , ld , LDNum , LEAd , LEADF , LZEta , m , MEM , 
     &          mrange , mua , MULti , NDIv , NPT , NSTart , NSTop , 
     &          NSW , nz
      COMPLEX*16 ARM , STAMP , dis , uhuj
      COMMON /CLCOM / LAMda(8) , LEAd(2,500) , LDNum(8,75) , LAMmax , 
     &                MULti(8)
      COMMON /AZ    / ARM(600,7)
      COMMON /CAUX  / NPT , NDIv , KDIv , LAMr(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /CCOUP / ZETa(50000) , LZEta(8)
      COMMON /CLCOM8/ CAT(600,3) , ISMax
      COMMON /COMME / ELM(500) , ELMu(500) , ELMl(500) , SA(500)
      COMMON /CEXC0 / NSTart(76) , NSTop(75)
      COMMON /KIN   / EPS(50) , EROot(50) , FIEx(50,2) , IEXp , IAXs(50)
      COMMON /CXI   / XI(500)
      ppp = 0.
      epsi = EPS(IEXp)
      errt = EROot(IEXp)
      rmir = CAT(Ir,3)
      DO i1 = 1 , LAMmax
         lam = LAMda(i1)
         nz = LZEta(lam)
         IF ( LAMr(lam).NE.0 ) THEN
            la = lam
            IF ( lam.GT.6 ) lam = lam - 6
            ld = LDNum(la,1)
            IF ( ld.NE.0 ) THEN
               DO i2 = 1 , ld
                  m = LEADF(1,i2,la)
                  indx = MEM(1,m,la)
                  xiv = XI(indx)
                  ismin = 0
                  is1 = NSTart(m)
                  IF ( NSTart(m).NE.0 ) THEN
                     isplus = INT(rmir-CAT(is1,3)) - lam
                     IF ( isplus.LT.0 ) THEN
                        ismin = isplus
                        isplus = 0
                     ENDIF
                     is2 = is1 + isplus - 1
                     mrange = 2*lam + 1 + ismin
                     IF ( is2+mrange.GT.NSTop(m) ) mrange = NSTop(m)
     &                    - is2
                     IF ( mrange.GT.0 ) THEN
                        DO i3 = 1 , mrange
                           is = is2 + i3
                           nz = nz + 1
                           z = ZETa(nz)
                           rmis = CAT(is,3)
                           rmu = rmis - rmir
                           mua = ABS(rmu) + 1.1
                           IF ( lam.LE.6 .OR. mua.NE.1 ) THEN
                              CALL FAZA1(la,mua,rmir,rmis,dis,rmu)
                              pm = ELM(indx)*z
                              uhuj = STAMP(epsi,errt,xiv,.03D0,W0,lam,
     &                               mua)
                              ARM(is,5) = dis*pm*uhuj
                              ppp = ppp + TCABS(ARM(is,5))
     &                              *TCABS(ARM(is,5))
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      ARM(Ir,5) = CMPLX(SQRT(1.-ppp),0.)
      END
