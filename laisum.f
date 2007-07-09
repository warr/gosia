 
C----------------------------------------------------------------------
 
      SUBROUTINE LAISUM(Ir,N,Rsg,Lam,Ld,Nz,I57)
      IMPLICIT NONE
      REAL*8 ACCa , ACCur , CAT , D2W , DIPol , ELM , ELMl , ELMu , EN , 
     &       q , rmir , rmis , rmu , Rsg , SA , SPIn , z , ZETa , ZPOl
      INTEGER*4 i2 , i3 , I57 , iii , indq , indx , Ir , irs , is , 
     &          is1 , is2 , ISG , ISG1 , ISHa , ISMax , ismin , ISO , 
     &          isplus , ISStar , ISSto
      INTEGER*4 KDIv , la , Lam , LAMr , Ld , LOCq , LP1 , LP10 , LP11 , 
     &          LP12 , LP13 , LP14 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , 
     &          LP9 , LZEta
      INTEGER*4 m , mrange , MSTore , mua , N , NDIv , NPT , NSTart , 
     &          NSTop , NSW , Nz
      COMPLEX*16 ARM , FAZA , pamp , EXPo , pamp1
      COMMON /PSPIN / ISHa
      COMMON /AZ    / ARM(600,7)
      COMMON /COEX  / EN(75) , SPIn(75) , ACCur , DIPol , ZPOl , ACCa , 
     &                ISO
      COMMON /CAUX  / NPT , NDIv , KDIv , LAMr(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /PINT  / ISStar(76) , ISSto(75) , MSTore(2,75)
      COMMON /ADBXI / EXPo(500)
      COMMON /CCOUP / ZETa(50000) , LZEta(8)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /CLCOM8/ CAT(600,3) , ISMax
      COMMON /COMME / ELM(500) , ELMu(500) , ELMl(500) , SA(500)
      COMMON /ALLC  / LOCq(8,7)
      COMMON /CEXC0 / NSTart(76) , NSTop(75)
      rmir = CAT(Ir,3)
      iii = 0
      IF ( Lam.GT.6 ) iii = 1
      la = Lam
      IF ( Lam.GT.6 ) Lam = Lam - 6
      DO i2 = 1 , Ld
         pamp = (0.,0.)
         m = MSTore(1,i2)
         indx = MSTore(2,i2)
         ismin = 0
         is1 = NSTart(m)
         IF ( is1.NE.0 ) THEN
            isplus = INT(rmir-CAT(is1,3)) - Lam
            IF ( isplus.LT.0 ) THEN
               ismin = isplus
               isplus = 0
            ENDIF
            is2 = is1 + isplus - 1
            mrange = 2*Lam + 1 + ismin
            IF ( is2+mrange.GT.NSTop(m) ) mrange = NSTop(m) - is2
            IF ( mrange.GT.0 ) THEN
               DO i3 = 1 , mrange
                  is = is2 + i3
                  rmis = CAT(is,3)
                  IF ( ISO.NE.0 .OR. rmir.LE..1 .OR. rmis.LE..1 ) THEN
                     rmu = rmis - rmir
                     mua = ABS(rmu) + 1.1
                     IF ( la.LE.6 .OR. mua.NE.1 ) THEN
                        indq = LOCq(Lam,mua) + NPT
                        Nz = Nz + 1
                        z = ZETa(Nz)
                        q = ZETa(indq+LP7)
                        IF ( NDIv.NE.0 ) q = ZETa(indq+LP7) + DBLE(KDIv)
     &                       *(ZETa(indq+LP7+ISG1)-ZETa(indq+LP7))
     &                       /DBLE(NDIv)
                        pamp1 = FAZA(la,mua,rmu,Rsg)*q*z
                        IF ( ISO.NE.0 .OR. rmir.LE..1 ) THEN
                           pamp = pamp1*ARM(is,I57) + pamp
                           IF ( ISO.EQ.0 .AND. rmis.GT..1 ) GOTO 10
                        ENDIF
                        IF ( N.NE.m ) THEN
                           irs = (-1)**(INT(rmir+rmis)-ISHa+iii)
                           ARM(is,6) = ARM(is,6) + irs*pamp1*ARM(Ir,I57)
                           ISStar(i2) = MIN(is,ISStar(i2))
                           ISSto(i2) = MAX(is,ISSto(i2))
                        ENDIF
                     ENDIF
                  ENDIF
 10            ENDDO
               IF ( N.EQ.m ) THEN
                  ARM(Ir,4) = ARM(Ir,4) + pamp*ELM(indx)
               ELSE
                  ARM(Ir,4) = ARM(Ir,4) + pamp*ELM(indx)*EXPo(indx)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      Lam = la
      END
