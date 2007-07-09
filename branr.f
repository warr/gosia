 
C----------------------------------------------------------------------
 
      SUBROUTINE BRANR(Chisq,Nwyr,Chilo)
      IMPLICIT NONE
      REAL*8 ACCa , ACCur , BRAt , ch1 , ch2 , Chilo , Chisq , CONV , 
     &       DELta , DIPol , ELM , ELMl , ELMu , EN , ENDec , eng1 , 
     &       eng2 , ENZ , SA , SPIn
      REAL*8 TAU , u , ZPOl
      INTEGER*4 i1 , i2 , IBRc , iflg , iout , IPRm , ISO , ITMa , itt , 
     &          j1 , j2 , k , KSEq , lab1 , lab2 , LAMda , LAMmax , 
     &          LDNum , LEAd , mul2
      INTEGER*4 MULti , n1 , n2 , NBRa , Nwyr
      COMMON /CLCOM / LAMda(8) , LEAd(2,500) , LDNum(8,75) , LAMmax , 
     &                MULti(8)
      COMMON /COEX  / EN(75) , SPIn(75) , ACCur , DIPol , ZPOl , ACCa , 
     &                ISO
      COMMON /BRNCH / BRAt(50,2) , IBRc(2,50) , NBRa
      COMMON /TRA   / DELta(500,3) , ENDec(500) , ITMa(50,200) , 
     &                ENZ(200)
      COMMON /PRT   / IPRm(20)
      COMMON /COMME / ELM(500) , ELMu(500) , ELMl(500) , SA(500)
      COMMON /LEV   / TAU(75) , KSEq(500,4)
      IF ( NBRa.EQ.0 ) RETURN
      IF ( IPRm(3).EQ.-1 ) WRITE (22,99001)
99001 FORMAT (1X,///10X,'EXP. AND CALCULATED BRANCHING RATIOS',//5X,
     &        'NS1',5X,'NF1',5X,'NS2',5X,'NF2',5X,'RATIO(1:2)',9X,
     &        'ERROR',7X,'CALC.RATIO',5X,'(EXP-CAL)/ERROR',//)
      Nwyr = Nwyr + NBRa
      mul2 = MULti(1) + MULti(2)
      DO k = 1 , NBRa
         ch1 = 0.
         ch2 = 0.
         iflg = 1
         itt = 1
         iout = 0
         n1 = IBRc(1,k)
         n2 = IBRc(2,k)
         i1 = KSEq(n1,1)
         i2 = KSEq(n2,1)
         eng1 = EN(KSEq(n1,3)) - EN(KSEq(n1,4))
         eng2 = EN(KSEq(n2,3)) - EN(KSEq(n2,4))
         IF ( i1.NE.0 ) THEN
            IF ( i1.LE.MULti(1) ) lab1 = 1
            IF ( i1.GT.MULti(1) .AND. i1.LE.mul2 ) lab1 = 2
            IF ( i1.GT.mul2 ) lab1 = 3
         ENDIF
         IF ( i2.NE.0 ) THEN
            IF ( i2.LE.MULti(1) ) lab2 = 1
            IF ( i2.GT.MULti(1) .AND. i2.LE.mul2 ) lab2 = 2
            IF ( i2.GT.mul2 ) lab2 = 3
         ENDIF
         IF ( i1.NE.0 ) ch1 = ELM(i1)*ELM(i1)*DELta(n1,1)
     &                        /(1.+CONV(eng1,lab1))
         IF ( i2.NE.0 ) ch2 = ELM(i2)*ELM(i2)*DELta(n2,1)
     &                        /(1.+CONV(eng2,lab2))
         j1 = KSEq(n1,2)
         IF ( j1.NE.0 ) THEN
            iflg = iflg + 1
            lab1 = lab1 + 2
            ch1 = ch1 + ELM(j1)*ELM(j1)*DELta(n1,2)/(1.+CONV(eng1,lab1))
         ENDIF
         j2 = KSEq(n2,2)
         IF ( j2.NE.0 ) THEN
            iflg = iflg + 1
            lab2 = lab2 + 2
            ch2 = ch2 + ELM(j2)*ELM(j2)*DELta(n2,2)/(1.+CONV(eng2,lab2))
         ENDIF
         u = (ch1/ch2-BRAt(k,1))/BRAt(k,2)
         Chisq = Chisq + u*u
         Chilo = Chilo + (BRAt(k,1)*LOG(ch1/ch2/BRAt(k,1))/BRAt(k,2))**2
         IF ( IPRm(3).EQ.-1 ) WRITE (22,99002) KSEq(n1,3) , KSEq(n1,4) , 
     &                               KSEq(n2,3) , KSEq(n2,4) , BRAt(k,1)
     &                               , BRAt(k,2) , ch1/ch2 , -u
99002    FORMAT (5X,3(1I2,6X),1I2,5X,3(1F10.5,5X),5X,1F4.1)
      ENDDO
      IF ( IPRm(3).EQ.-1 ) IPRm(3) = 0
      END
