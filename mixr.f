 
C----------------------------------------------------------------------
 
      SUBROUTINE MIXR(Nw,Ipsw,Chi,Chilo)
      IMPLICIT NONE
      REAL*8 Chi , Chilo , dl , DMIx , DMIxe , ELM , ELMl , ELMu , SA , 
     &       TAU
      INTEGER*4 i , IMIx , INTr , inx , inx1 , IPS1 , Ipsw , it , KSEq , 
     &          LNY , NDL , Nw
      COMMON /LEV   / TAU(75) , KSEq(500,4)
      COMMON /COMME / ELM(500) , ELMu(500) , ELMl(500) , SA(500)
      COMMON /MIXD  / DMIxe(20,2) , DMIx(20) , IMIx(20) , NDL
      COMMON /LOGY  / LNY , INTr , IPS1
      IF ( NDL.EQ.0 ) RETURN
      Nw = Nw + NDL
      DO i = 1 , NDL
         it = IMIx(i)
         inx = KSEq(it,1)
         inx1 = KSEq(it,2)
         IF ( ABS(ELM(inx1)).LT.1.E-5 ) ELM(inx1) = 1.E-5
         dl = DMIx(i)*ELM(inx)/ELM(inx1)
         IF ( Ipsw.EQ.1 ) DMIx(i) = dl
         Chi = Chi + (dl-DMIxe(i,1))**2/DMIxe(i,2)/DMIxe(i,2)
         IF ( LNY.EQ.1 ) Chilo = Chilo + 
     &                           (DMIxe(i,1)*LOG(ABS(dl/DMIxe(i,1)))
     &                           /DMIxe(i,2))**2
      ENDDO
      IF ( Ipsw.EQ.0 ) RETURN
      WRITE (22,99001)
99001 FORMAT (1X//10X,'E2/M1 MIXING RATIOS'/10X,'TRANSITION',10X,
     &        'EXP.DELTA',10X,'CALC.DELTA',10X,'SIGMA'/)
      DO i = 1 , NDL
         dl = (DMIx(i)-DMIxe(i,1))/DMIxe(i,2)
         it = IMIx(i)
         WRITE (22,99002) KSEq(it,3) , KSEq(it,4) , DMIxe(i,1) , DMIx(i)
     &                    , dl
99002    FORMAT (10X,1I2,'---',1I2,14X,1F7.2,12X,1F7.2,13X,1F5.2)
      ENDDO
      END
