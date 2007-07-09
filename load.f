 
C----------------------------------------------------------------------
 
      SUBROUTINE LOAD(Iexp,Ient,Icg,Polm,Joj)
      IMPLICIT NONE
      REAL*8 a1 , a2 , aaz2 , aaz3 , aazz , ACCa , ACCur , ah , CAT , 
     &       cpsi , dep , DIPol , EMMa , EN , EP , eta , etan , Polm , 
     &       pp1 , pp2
      REAL*8 ppp , PSI , QAPr , rlam , SPIn , ssqrt , szet , TLBdg , 
     &       VINf , wrt , wrtm , XA , XA1 , XI , z1 , z2 , zet , ZETa , 
     &       ZPOl , zsqa
      INTEGER*4 i , i1 , i2 , i3 , IAPr , Icg , Ient , Iexp , IPAth , 
     &          ir , is , ISEx , ISHa , ISMax , ISO , ispi , ispo , 
     &          IVAr , IZ , IZ1
      INTEGER*4 jj , jjj , Joj , la , lam , lam1 , LAMda , LAMmax , ld , 
     &          LDNum , LEAd , LMAx , LMAxe , LP1 , LP10 , LP11 , LP12 , 
     &          LP13 , LP14 , LP2
      INTEGER*4 LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , LZEta , m , m1 , 
     &          m2 , MAGa , MAGexc , MEMax , MEMx6 , mstop , MULti , n , 
     &          n2 , n3 , NCM
      INTEGER*4 NDIm , NEXpt , NMAx , NMAx1 , nn , NSTart , NSTop , nz
      LOGICAL ERR
      COMMON /CLCOM / LAMda(8) , LEAd(2,500) , LDNum(8,75) , LAMmax , 
     &                MULti(8)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /COEX  / EN(75) , SPIn(75) , ACCur , DIPol , ZPOl , ACCa , 
     &                ISO
      COMMON /PSPIN / ISHa
      COMMON /CEXC  / MAGexc , MEMax , LMAxe , MEMx6 , IVAr(500)
      COMMON /PCOM  / PSI(500)
      COMMON /CCOUP / ZETa(50000) , LZEta(8)
      COMMON /CLM   / LMAx
      COMMON /CLCOM8/ CAT(600,3) , ISMax
      COMMON /CLCOM9/ ERR
      COMMON /COEX2 / NMAx , NDIm , NMAx1
      COMMON /CX    / NEXpt , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBdg(50) , VINf(50)
      COMMON /CEXC0 / NSTart(76) , NSTop(75)
      COMMON /CXI   / XI(500)
      COMMON /CAUX0 / EMMa(75) , NCM
      COMMON /APRCAT/ QAPr(500,2,7) , IAPr(500,2) , ISEx(75)
      COMMON /PTH   / IPAth(75) , MAGa(75)
      DIMENSION etan(75) , cpsi(8)
      LMAx = INT(SPIn(1)+1.1)
      IF ( Ient.EQ.1 ) THEN
         ISHa = 0
         ispi = INT(SPIn(1)+.51)
         ispo = INT(SPIn(1)+.49)
         IF ( ispi.NE.ispo ) ISHa = 1
         z1 = DBLE(ABS(IZ1(Iexp)))
         z2 = DBLE(IZ)
         a1 = XA1(Iexp)
         a2 = XA
         ZPOl = DIPol*EP(Iexp)*a2/(z2*z2*(1.+a1/a2))
         IF ( IZ1(Iexp).LT.0 ) ZPOl = DIPol*EP(Iexp)
     &                                *a1/(z1*z1*(1.+a2/a1))
         IF ( IZ1(Iexp).LE.0 ) THEN
            ah = a1
            a1 = a2
            a2 = ah
         ENDIF
         eta = z1*z2*SQRT(a1/EP(Iexp))/6.349770
         DO m = 1 , NMAx
            dep = (1.0+a1/a2)*EN(m)
            zet = dep/EP(Iexp)
            szet = SQRT(1.0-zet)
            etan(m) = eta/szet
         ENDDO
         DO n = 1 , MEMax
            i1 = LEAd(1,n)
            i2 = LEAd(2,n)
            XI(n) = etan(i1) - etan(i2)
         ENDDO
         aazz = 1./(1.+a1/a2)/z1/z2
         cpsi(1) = 5.169286*aazz
         IF ( LMAxe.NE.1 ) THEN
            aaz2 = aazz*aazz
            cpsi(2) = 14.359366*aaz2
            IF ( LMAxe.NE.2 ) THEN
               aaz3 = aazz*aaz2
               cpsi(3) = 56.982577*aaz3
               IF ( LMAxe.NE.3 ) THEN
                  aazz = aaz2*aaz2
                  cpsi(4) = 263.812653*aazz
                  IF ( LMAxe.NE.4 ) THEN
                     aaz2 = aaz3*aaz2
                     cpsi(5) = 1332.409500*aaz2
                     IF ( LMAxe.NE.5 ) THEN
                        aazz = aaz3*aaz3
                        cpsi(6) = 7117.691577*aazz
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF ( MAGexc.NE.0 ) THEN
            aazz = VINf(Iexp)/95.0981942
            cpsi(7) = aazz*cpsi(1)
            IF ( LAMmax.NE.8 ) cpsi(8) = aazz*cpsi(2)
         ENDIF
         zsqa = z1*SQRT(a1)
         i3 = 1
         ppp = 1. + a1/a2
         DO i1 = 1 , LAMmax
            lam = LAMda(i1)
            lam1 = lam
            IF ( lam.GT.6 ) lam1 = lam - 6
            DO n2 = 1 , NMAx
               nn = LDNum(lam,n2)
               IF ( nn.NE.0 ) THEN
                  n3 = LEAd(1,i3)
                  pp1 = EP(Iexp) - ppp*EN(n3)
                  DO m1 = 1 , nn
                     m2 = LEAd(2,i3)
                     i2 = i3
                     i3 = i3 + 1
                     pp2 = EP(Iexp) - ppp*EN(m2)
                     PSI(i2) = cpsi(lam)*zsqa*(pp1*pp2)
     &                         **((2.*DBLE(lam1)-1.)/4.)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         IF ( Ient.EQ.1 ) RETURN
      ENDIF
      DO n = 1 , NMAx
         NSTart(n) = 0
         NSTop(n) = 0
      ENDDO
      is = 1
      NSTart(1) = 1
      DO n = 1 , NMAx
         wrt = Polm - EMMa(Iexp)
         wrtm = Polm + EMMa(Iexp)
         IF ( Icg.EQ.2 ) wrt = Polm - DBLE(MAGa(Iexp))
         IF ( Icg.EQ.2 ) wrtm = Polm + DBLE(MAGa(Iexp))
         IF ( wrtm.LT.-SPIn(n) ) THEN
            NSTart(n) = 0
         ELSE
            IF ( ABS(wrt).GT.SPIn(n) ) wrt = -SPIn(n)
            IF ( wrtm.GT.SPIn(n) ) wrtm = SPIn(n)
            mstop = INT(wrtm-wrt+1.01)
            DO i = 1 , mstop
               CAT(is,1) = n
               CAT(is,2) = SPIn(n)
               CAT(is,3) = wrt + DBLE(i-1)
               IF ( n.EQ.1 .AND. ABS(CAT(is,3)-Polm).LT.1.E-6 ) Joj = is
               is = is + 1
            ENDDO
         ENDIF
         NSTart(n+1) = is
         NSTop(n) = is - 1
      ENDDO
      ISMax = is - 1
      IF ( ISMax.LE.LP10 ) THEN
         IF ( Ient.EQ.3 ) RETURN
         nz = 0
         DO jj = 1 , 7
            DO jjj = 1 , MEMax
               QAPr(jjj,1,jj) = 0.
               QAPr(jjj,2,jj) = 0.
            ENDDO
         ENDDO
         DO i = 1 , 8
            LZEta(i) = 0
         ENDDO
         DO i1 = 1 , LAMmax
            lam = LAMda(i1)
            IF ( Icg.NE.2 .OR. lam.LE.6 ) THEN
               la = lam
               IF ( lam.GT.6 ) lam = lam - 6
               rlam = DBLE(lam)
               ssqrt = SQRT(2.*rlam+1.)
               LZEta(la) = nz
               ir = 0
 10            ir = ir + 1
               IF ( ir.LE.ISMax ) THEN
                  n = CAT(ir,1)
                  IF ( Icg.NE.1 ) THEN
                     IF ( MAGa(Iexp).EQ.0 .AND. ir.NE.IPAth(n) ) GOTO 10
                     IF ( ABS(ir-IPAth(n)).GT.1 ) GOTO 10
                  ENDIF
                  ld = LDNum(la,n)
                  IF ( ld.EQ.0 ) THEN
                     ir = ir + NSTop(n) - NSTart(n)
                  ELSE
                     CALL LSLOOP(ir,n,nz,ld,lam,la,ssqrt,Icg,Iexp)
                  ENDIF
                  GOTO 10
               ENDIF
            ENDIF
         ENDDO
         IF ( nz.GT.LP7 ) THEN
            WRITE (22,99001) LP7
99001       FORMAT (1x,
     &              'ERROR - NUMBER OF ELEMENTS IN ZETA ARRAY EXCEEDS',
     &              'ZEMAX',5X,'(ZEMAX =',I6,')')
         ELSE
            RETURN
         ENDIF
      ELSE
         WRITE (22,99002) LP10
99002    FORMAT (' ERROR-ISMAX EXCEEDS MAGMAX',5X,'(MAGMAX =',I4,')')
      ENDIF
      ERR = .TRUE.
      END
