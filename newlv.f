 
C----------------------------------------------------------------------
 
      SUBROUTINE NEWLV(N,Ld,La)
      IMPLICIT NONE
      REAL*8 D2W
      INTEGER*4 i2 , IFLg , indx , ISG , ISG1 , ISStar , ISSto , KDIv , 
     &          La , LAMda , LAMmax , LAMr , Ld , LDNum , LEAd , LEADF , 
     &          m , MEM , MSTore , MULti
      INTEGER*4 N , NDIv , NPT , NSTart , NSTop , NSW
      COMPLEX*16 EXPo , EXPON
      COMMON /CLCOM / LAMda(8) , LEAd(2,500) , LDNum(8,75) , LAMmax , 
     &                MULti(8)
      COMMON /CAUX  / NPT , NDIv , KDIv , LAMr(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /PINT  / ISStar(76) , ISSto(75) , MSTore(2,75)
      COMMON /ADBXI / EXPo(500)
      COMMON /FLA   / IFLg
      COMMON /CEXC0 / NSTart(76) , NSTop(75)
      Ld = LDNum(La,N)
      IF ( Ld.EQ.0 ) RETURN
      DO i2 = 1 , Ld
         m = LEADF(N,i2,La)
         ISStar(i2) = NSTop(m)
         ISSto(i2) = NSTart(m)
         MSTore(1,i2) = m
         indx = MEM(N,m,La)
         MSTore(2,i2) = indx
         IF ( IFLg.NE.0 ) THEN
            IF ( m.NE.N ) EXPo(indx) = EXPON(indx,NPT,ISG,ISG1,NDIv,KDIv
     &                                 )
         ENDIF
      ENDDO
      END
