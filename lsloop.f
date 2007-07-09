 
C----------------------------------------------------------------------
 
      SUBROUTINE LSLOOP(Ir,N,Nz,Ld,Lam,La,Ssqrt,Icg,Iexp)
      IMPLICIT NONE
      REAL*8 ACCa , ACCur , CAT , DIPol , ELM , ELMl , ELMu , EN , phz , 
     &       PSI , QAPr , rmir , rmis , SA , SPIn , Ssqrt , WTHREJ , 
     &       ZETa , ZPOl
      INTEGER*4 i2 , i3 , IAPr , Icg , Iexp , IFAc , iiex , indx , 
     &          inqa , inr , ins , IPAth , Ir , is , is1 , is2 , ISEx , 
     &          ISMax , ismin , ISO
      INTEGER*4 isplus , jg1 , jg2 , jrmir , La , Lam , lam2 , Ld , 
     &          LEADF , LP1 , LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , 
     &          LP3 , LP4 , LP6 , LP7
      INTEGER*4 LP8 , LP9 , LZEta , m , MAGa , MEM , mrange , mt , N , 
     &          NSTart , NSTop , Nz
      COMMON /COEX  / EN(75) , SPIn(75) , ACCur , DIPol , ZPOl , ACCa , 
     &                ISO
      COMMON /PCOM  / PSI(500)
      COMMON /CCOUP / ZETa(50000) , LZEta(8)
      COMMON /CLCOM8/ CAT(600,3) , ISMax
      COMMON /CEXC0 / NSTart(76) , NSTop(75)
      COMMON /APRCAT/ QAPr(500,2,7) , IAPr(500,2) , ISEx(75)
      COMMON /PTH   / IPAth(75) , MAGa(75)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /COMME / ELM(500) , ELMu(500) , ELMl(500) , SA(500)
      COMMON /CLCOM0/ IFAc(75)
      lam2 = 2*Lam
      inr = CAT(Ir,2)*2.
      rmir = CAT(Ir,3)
      jrmir = 2.*rmir
      DO i2 = 1 , Ld
         m = LEADF(N,i2,La)
         indx = MEM(N,m,La)
         IAPr(indx,1) = N
         IAPr(indx,2) = m
         ismin = 0
         ins = SPIn(m)*2.
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
                  IF ( ISO.NE.0 .OR. rmis.LE..1 .OR. rmir.LE..1 ) THEN
                     jg1 = -rmis*2.
                     jg2 = (rmis-rmir)*2.
                     IF ( Icg.NE.2 .OR. ABS(jg2).LE.2*MAGa(Iexp) ) THEN
                        IF ( La.LE.6 .OR. jg2.NE.0 ) THEN
                           Nz = Nz + 1
                           IF ( Nz.LE.LP7 ) THEN
                              iiex = (ins+jg1)/2
                              phz = (-1.0)**iiex
                              ZETa(Nz) = phz*PSI(indx)
     &                           *Ssqrt*WTHREJ(ins,lam2,inr,jg1,jg2,
     &                           jrmir)
                              IF ( Icg.NE.1 ) THEN
                                 mt = CAT(is,1)
                                 CALL CODE7(Ir,is,N,mt,inqa,indx)
                                 IF ( ABS(ELM(indx)).LT.1.E-6 )
     &                                ELM(indx) = 1.E-6
                                 IF ( inqa.NE.-1 ) THEN
                                    QAPr(indx,1,inqa) = ZETa(Nz)
     &                                 *ELM(indx)
                                    IF ( ISO.EQ.0 .AND. inqa.EQ.1 )
     &                                 QAPr(indx,1,7) = QAPr(indx,1,1)
     &                                 *IFAc(m)
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      END
