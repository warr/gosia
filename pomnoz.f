 
C----------------------------------------------------------------------
 
      SUBROUTINE POMNOZ(Acca,L,Iw,Ktoto,Img,Jidim)
      IMPLICIT NONE
      REAL*8 Acca , QAPr , sig , TCABS , test , u
      INTEGER*4 IAPr , IDIve , Img , INHb , IPAth , ISEx , IVAr , Iw , 
     &          Jidim , k , kk , Ktoto , L , LERf , LMAxe , m , MAGa , 
     &          MAGexc , mc , mc1
      INTEGER*4 MEMax , MEMx6 , mw , mw1
      COMPLEX*16 ARM , ci
      COMMON /INHI  / INHb
      COMMON /APRCAT/ QAPr(500,2,7) , IAPr(500,2) , ISEx(75)
      COMMON /PTH   / IPAth(75) , MAGa(75)
      COMMON /CEXC  / MAGexc , MEMax , LMAxe , MEMx6 , IVAr(500)
      COMMON /AZ    / ARM(600,7)
      COMMON /APRX  / LERf , IDIve(50,2)
      DATA ci/(0.,-1.)/
      sig = 1.
      IF ( L.NE.2 ) sig = -1.
      DO kk = 1 , Jidim
         ARM(kk,1) = ARM(kk,Iw)
      ENDDO
      DO k = 1 , 100
         Ktoto = Ktoto + 1
         DO m = 1 , MEMx6
            mw1 = IAPr(m,1)
            mc1 = IAPr(m,2)
            IF ( IPAth(mw1).NE.0 .AND. IPAth(mc1).NE.0 ) THEN
               mw = IPAth(mw1) + 1
               mc = IPAth(mc1) + 1
               IF ( Ktoto.GE.ISEx(mc1) ) THEN
                  IF ( Img.EQ.1 ) THEN
                     ARM(mw,2) = ARM(mw,2) + QAPr(m,L,4)*ARM(mc,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPr(m,L,4)*ARM(mw,1)
                  ELSE
                     ARM(mw,2) = ARM(mw,2) + QAPr(m,L,4)*ARM(mc,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPr(m,L,4)*ARM(mw,1)
                     ARM(mw-1,2) = ARM(mw-1,2) + QAPr(m,L,2)*ARM(mc,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPr(m,L,2)*ARM(mw-1,1)
                     ARM(mw-1,2) = ARM(mw-1,2) + QAPr(m,L,1)*ARM(mc-1,1)
                     ARM(mc-1,2) = ARM(mc-1,2) + sig*QAPr(m,L,1)
     &                             *ARM(mw-1,1)
                     ARM(mw,2) = ARM(mw,2) + QAPr(m,L,3)*ARM(mc-1,1)
                     ARM(mc-1,2) = ARM(mc-1,2) + sig*QAPr(m,L,3)
     &                             *ARM(mw,1)
                     ARM(mw,2) = ARM(mw,2) + QAPr(m,L,5)*ARM(mc+1,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPr(m,L,6)*ARM(mw+1,1)
                     ARM(mc+1,2) = ARM(mc+1,2) + sig*QAPr(m,L,5)
     &                             *ARM(mw,1)
                     ARM(mw+1,2) = ARM(mw+1,2) + QAPr(m,L,6)*ARM(mc,1)
                     ARM(mw+1,2) = ARM(mw+1,2) + QAPr(m,L,7)*ARM(mc+1,1)
                     ARM(mc+1,2) = ARM(mc+1,2) + sig*QAPr(m,L,7)
     &                             *ARM(mw+1,1)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         test = 0.
         DO m = 1 , Jidim
            ARM(m,1) = ARM(m,2)/k
            ARM(m,2) = (0.,0.)
            IF ( L.NE.1 ) ARM(m,1) = ARM(m,1)*ci
            ARM(m,Iw) = ARM(m,Iw) + ARM(m,1)
            IF ( k.GT.5 ) THEN
               u = TCABS(ARM(m,Iw))
               test = test + u*u
            ENDIF
         ENDDO
         IF ( ABS(test-1.).LT.Acca ) GOTO 99999
      ENDDO
      IF ( INHb.NE.1 ) LERf = 1
99999 END
