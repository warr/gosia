 
C----------------------------------------------------------------------
 
      SUBROUTINE PTICC(Idr)
      IMPLICIT NONE
      REAL*8 ACCa , ACCur , cone1 , cone2 , conm1 , CONV , DIPol , EN , 
     &       enet , SPIn , TAU , ZPOl
      INTEGER*4 Idr , iinx , ISO , KSEq , l , LAMda , LAMmax , LDNum , 
     &          LEAd , MULti , nf , ni
      COMMON /CLCOM / LAMda(8) , LEAd(2,500) , LDNum(8,75) , LAMmax , 
     &                MULti(8)
      COMMON /COEX  / EN(75) , SPIn(75) , ACCur , DIPol , ZPOl , ACCa , 
     &                ISO
      COMMON /LEV   / TAU(75) , KSEq(500,4)
      WRITE (22,99001)
99001 FORMAT (1X//20X,'CALCULATED INTERNAL CONVERSION ',
     &        'COEFFICIENTS FOR E1,E2 AND M1'//5X,'NI',5X,'NF',7X,'II',
     &        8X,'IF',9X,'ENERGY(MEV)',6X,'ICC(E1)',8X,'ICC(E2)',8X,
     &        'ICC(M1)')
      DO l = 1 , Idr
         iinx = KSEq(l,1)
         ni = KSEq(l,3)
         nf = KSEq(l,4)
         enet = EN(ni) - EN(nf)
         cone2 = CONV(enet,2)
         IF ( ABS(SPIn(ni)-SPIn(nf)).GT.2. ) cone2 = 0.
         conm1 = 0.
         cone1 = 0.
         IF ( iinx.LE.MULti(1) ) cone1 = CONV(enet,1)
         IF ( ABS(SPIn(ni)-SPIn(nf)).LT.2. ) conm1 = CONV(enet,4)
         WRITE (22,99002) ni , nf , SPIn(ni) , SPIn(nf) , enet , cone1 , 
     &                    cone2 , conm1
99002    FORMAT (5X,I2,5X,I2,7X,F4.1,6X,F4.1,9X,F6.4,8X,E9.4,6X,E9.4,6X,
     &           E9.4)
      ENDDO
      END
