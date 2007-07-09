 
C----------------------------------------------------------------------
 
      SUBROUTINE PRELM(Iop)
      IMPLICIT NONE
      REAL*8 ACCa , ACCur , b , DIPol , ELM , ELMl , ELMu , EN , HLM , 
     &       pv , SA , SPIn , ste , ZPOl
      INTEGER*4 inx , Iop , ISO , isp , IVAr , j , k , kk , l , LAMda , 
     &          LAMmax , LDNum , LEAd , LMAxe , m , MAGexc , MEMax , 
     &          MEMx6 , MULti , NDIm
      INTEGER*4 NMAx , NMAx1
      CHARACTER*3 wrn
      COMMON /HHH   / HLM(500)
      COMMON /COMME / ELM(500) , ELMu(500) , ELMl(500) , SA(500)
      COMMON /CEXC  / MAGexc , MEMax , LMAxe , MEMx6 , IVAr(500)
      COMMON /COEX  / EN(75) , SPIn(75) , ACCur , DIPol , ZPOl , ACCa , 
     &                ISO
      COMMON /CLCOM / LAMda(8) , LEAd(2,500) , LDNum(8,75) , LAMmax , 
     &                MULti(8)
      COMMON /COEX2 / NMAx , NDIm , NMAx1
      inx = 0
      WRITE (22,99001)
99001 FORMAT (2X/40X,'MATRIX ELEMENTS',//)
      DO j = 1 , 8
         m = MULti(j)
         IF ( m.NE.0 ) THEN
            WRITE (22,99002) j
99002       FORMAT (5X,'MULTIPOLARITY=',1I1)
            IF ( Iop.EQ.1 ) WRITE (22,99003)
99003       FORMAT (4X,'INDEX',3X,'NF',5X,'NS',10X,'ME')
            IF ( Iop.EQ.2 ) WRITE (22,99004)
99004       FORMAT (4X,'INDEX',3X,'NF',5X,'NS',10X,'ME',15X,'LIMITS')
            IF ( Iop.EQ.3 ) WRITE (22,99005)
99005       FORMAT (4X,'INDEX',3X,'NF',5X,'NS',10X,'ME',10X,'PC CHANGE',
     &              5X,'RED. TRANS. PROB.')
            DO k = 1 , NMAx
               l = LDNum(j,k)
               IF ( l.NE.0 ) THEN
                  DO kk = 1 , l
                     inx = inx + 1
                     IF ( Iop.EQ.2 ) THEN
                        IF ( IVAr(inx).EQ.0 ) THEN
                           WRITE (22,99006) inx , LEAd(1,inx) , 
     &                            LEAd(2,inx) , ELM(inx)
99006                      FORMAT (5X,1I3,5X,1I2,5X,1I2,5X,1F10.5,5X,
     &                             'FIXED')
                        ELSEIF ( IVAr(inx).GT.1000 ) THEN
                           WRITE (22,99007) inx , LEAd(1,inx) , 
     &                            LEAd(2,inx) , ELM(inx) , 
     &                            (IVAr(inx)-1000)
99007                      FORMAT (5X,1I3,5X,1I2,5X,1I2,5X,1F10.5,5X,
     &                             'COUPLED TO',1X,1I3)
                        ELSE
                           WRITE (22,99009) inx , LEAd(1,inx) , 
     &                            LEAd(2,inx) , ELM(inx) , ELMl(inx) , 
     &                            ELMu(inx)
                        ENDIF
                     ELSEIF ( Iop.EQ.3 ) THEN
                        isp = LEAd(2,inx)
                        pv = (ELMu(inx)-ELMl(inx))/100.
                        wrn = '   '
                        IF ( (ELM(inx)-ELMl(inx)).LT.pv ) wrn = '*?*'
                        IF ( (ELMu(inx)-ELM(inx)).LT.pv ) wrn = '*?*'
                        ste = HLM(inx)
                        b = ELM(inx)*ELM(inx)/(2.*SPIn(isp)+1.)
                        IF ( LEAd(1,inx).EQ.LEAd(2,inx) ) b = 9999999.
                        WRITE (22,99009) inx , LEAd(1,inx) , LEAd(2,inx)
     &                         , ELM(inx) , 100.*(ELM(inx)-ste)/ste , 
     &                         b , wrn
                     ELSE
                        WRITE (22,99008) inx , LEAd(1,inx) , LEAd(2,inx)
     &                         , ELM(inx)
99008                   FORMAT (5X,1I3,5X,1I2,5X,1I2,5X,1F10.5)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
99009 FORMAT (5X,1I3,5X,1I2,5X,1I2,3(5X,1F10.5),1A3)
      END
