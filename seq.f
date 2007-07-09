 
C----------------------------------------------------------------------
 
      SUBROUTINE SEQ(Idr)
      IMPLICIT NONE
      REAL*8 ACCa , ACCur , CONV , DELta , DIPol , ega , egs , emax , 
     &       EN , ENDec , ENZ , F , FP , GF , GKP , SPIn , spinf , 
     &       spini , TAU , twoi
      REAL*8 ZPOl
      INTEGER*4 idecay , Idr , indx , inx , inx1 , ir , is , ISO , 
     &          istr1 , istr2 , ITMa , j , js , jsave , k , KLEc , kpa , 
     &          KSEq , l , la
      INTEGER*4 la1 , LAMda , LAMmax , ld , LDNum , LEAd , LEADF , LP1 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , LP3 , LP4 , 
     &          LP6 , LP7 , LP8 , LP9
      INTEGER*4 m , m1 , m6 , MEM , mk , mule , mulm , MULti , n , n1 , 
     &          NDIm , NMAx , NMAx1 , nob
      COMMON /COEX2 / NMAx , NDIm , NMAx1
      COMMON /TRA   / DELta(500,3) , ENDec(500) , ITMa(50,200) , 
     &                ENZ(200)
      COMMON /CLCOM / LAMda(8) , LEAd(2,500) , LDNum(8,75) , LAMmax , 
     &                MULti(8)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /COEX  / EN(75) , SPIn(75) , ACCur , DIPol , ZPOl , ACCa , 
     &                ISO
      COMMON /LEV   / TAU(75) , KSEq(500,4)
      COMMON /CATLF / FP(4,500,3) , GKP(4,500,2) , KLEc(75)
      m6 = 0
      DO l = 1 , 6
         m6 = m6 + MULti(l)
      ENDDO
      idecay = 0
      Idr = 0
      DO l = 1 , LP3
         KLEc(l) = 0
      ENDDO
      DO k = 1 , LP2
         DO j = 1 , 3
            DO l = 1 , 4
               FP(l,k,j) = 0.
               IF ( j.NE.3 ) GKP(l,k,j) = 0.
            ENDDO
            DELta(k,j) = 0.
         ENDDO
      ENDDO
      DO n = 1 , NMAx
         TAU(n) = EN(n)
      ENDDO
      DO n = 1 , NMAx
         emax = 0.
         DO j = 1 , NMAx
            IF ( TAU(j).GE.emax ) THEN
               emax = TAU(j)
               jsave = j
            ENDIF
         ENDDO
         DO is = 1 , NMAx
            DO la = 1 , 8
               IF ( la.LE.3 .OR. la.EQ.7 .OR. la.EQ.8 ) THEN
                  ld = LDNum(la,is)
                  IF ( ld.NE.0 ) THEN
                     DO ir = 1 , ld
                        m = LEADF(is,ir,la)
                        IF ( m.EQ.jsave .OR. is.EQ.jsave ) THEN
                           IF ( is.NE.jsave .OR. EN(m).LT.EN(is) ) THEN
                              IF ( m.NE.jsave .OR. EN(is).LT.EN(m) )
     &                             THEN
                                 indx = MEM(is,m,la)
                                 idecay = idecay + 1
                                 KSEq(idecay,1) = m
                                 KSEq(idecay,2) = is
                                 KSEq(idecay,3) = indx
                                 KSEq(idecay,4) = la + 10
                                 IF ( EN(m).LE.EN(is) ) THEN
                                    KSEq(idecay,1) = is
                                    KSEq(idecay,2) = m
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         TAU(jsave) = -1.
      ENDDO
      DO l = 1 , idecay
         istr1 = 0
         IF ( KSEq(l,4).LT.10 ) GOTO 200
         istr2 = 0
         n = KSEq(l,1)
         m = KSEq(l,2)
         inx = KSEq(l,3)
         la = KSEq(l,4) - 10
         ega = EN(n) - EN(m)
         twoi = 1./SQRT(2.*SPIn(n)+1.)
         spini = SPIn(n) + .001
         spinf = SPIn(m) + .001
         egs = SQRT(ega)*twoi
         js = l + 1
         la1 = 0
         inx1 = 0
         DO j = js , idecay
            IF ( KSEq(j,4).GE.10 ) THEN
               n1 = KSEq(j,1)
               m1 = KSEq(j,2)
               IF ( n1.EQ.n .AND. m1.EQ.m ) THEN
                  inx1 = KSEq(j,3)
                  la1 = KSEq(j,4) - 10
                  KSEq(j,4) = KSEq(j,4) - 10
               ENDIF
            ENDIF
         ENDDO
         KSEq(l,4) = KSEq(l,4) - 10
         Idr = Idr + 1
         mule = 0
         mulm = 0
         nob = 1
 50      IF ( la.LE.3 ) THEN
            IF ( la.EQ.1 ) THEN
               DELta(Idr,1) = 399.05*ega*egs
               mule = 1
               istr1 = 1
            ELSEIF ( la.EQ.2 ) THEN
               DELta(Idr,1) = 3.4928*egs*ega*ega
               mule = 2
               istr1 = 2
            ELSEIF ( la.EQ.3 ) THEN
               DELta(Idr,1) = .02391*ega*ega*ega*egs
               mule = 3
               istr1 = 3
            ELSE
               GOTO 100
            ENDIF
            GOTO 150
         ENDIF
 100     la = la - 6
         IF ( la.EQ.2 ) THEN
            DELta(Idr,2) = .0368*ega*ega*egs
            mulm = 2
            istr2 = 5
         ELSE
            DELta(Idr,2) = 4.1952*ega*egs
            mulm = 1
            istr2 = 4
         ENDIF
 150     IF ( nob.NE.2 ) THEN
            IF ( mule.NE.1 ) THEN
               nob = nob + 1
               IF ( la.GT.3 ) inx1 = inx
               IF ( la1.NE.0 ) THEN
                  la = la1
                  GOTO 50
               ENDIF
            ENDIF
            inx1 = 0
         ENDIF
         DELta(Idr,3) = DELta(Idr,1)*DELta(Idr,2)
         DELta(Idr,1) = DELta(Idr,1)*DELta(Idr,1)
         DELta(Idr,2) = DELta(Idr,2)*DELta(Idr,2)
         KSEq(Idr,1) = inx
         KSEq(Idr,2) = inx1
         KSEq(Idr,3) = n
         KSEq(Idr,4) = m
         IF ( inx.GT.m6 ) THEN
            KSEq(Idr,2) = inx
            KSEq(Idr,1) = 0
         ENDIF
         ENDec(Idr) = EN(n) - EN(m)
         DO mk = 1 , 7 , 2
            kpa = mk/2 + 1
            k = mk - 1
            IF ( mule.GE.3 .OR. k.NE.6 ) THEN
               GKP(kpa,Idr,1) = GF(k,spini,spinf,mule)*DELta(Idr,1)
     &                          *(1.+CONV(ega,istr1))
               GKP(kpa,Idr,2) = GF(k,spini,spinf,mulm)*DELta(Idr,2)
     &                          *(1.+CONV(ega,istr2))
               FP(kpa,Idr,1) = F(k,spini,spinf,mule,mule)*DELta(Idr,1)
               FP(kpa,Idr,3) = F(k,spini,spinf,mulm,mule)*DELta(Idr,3)
               FP(kpa,Idr,2) = F(k,spini,spinf,mulm,mulm)*DELta(Idr,2)
            ENDIF
         ENDDO
         DELta(Idr,1) = DELta(Idr,1)*(1.+CONV(ega,istr1))
         DELta(Idr,2) = DELta(Idr,2)*(CONV(ega,istr2)+1.)
         KLEc(n) = KLEc(n) + 1
 200  ENDDO
      NMAx1 = 0
      DO n = 1 , NMAx
         IF ( KLEc(n).NE.0 ) NMAx1 = NMAx1 + 1
      ENDDO
      END
