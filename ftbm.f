 
C----------------------------------------------------------------------
 
      SUBROUTINE FTBM(Icll,Chisq,Idr,Ncall,Chilo,Bten)
      IMPLICIT NONE
      REAL*8 ACCa , ACCur , AGEli , aval , Bten , CAT , CC , Chilo , 
     &       chis1 , CHIs11 , chish , Chisq , chisx , chx , CORf , 
     &       DIPol , DYEx , EG , ELM , ELMl
      REAL*8 ELMu , EMH , EN , EP , EPS , EROot , fc , FIEx , fx , 
     &       polm , pr , prop , Q , SA , SPIn , TAU , TLBdg , UPL , 
     &       val , VINf
      REAL*8 wz , XA , XA1 , YEXp , YNRm , ZETa , ZPOl
      INTEGER*4 i1 , i11 , iapx , IAXs , Icll , idec , Idr , IDRn , 
     &          IEXp , iflg , IGRd , ii , ILE , ile1 , ile2 , ile3 , 
     &          ilin , indx , inko , INM
      INTEGER*4 inp , inpo , inpx , INTr , inzz , inzzz , IPAth , IPRm , 
     &          IPS1 , ISMax , ISO , issp , ITAk2 , itemp , IVAr , ixx , 
     &          IY , IZ , IZ1 , izzz
      INTEGER*4 j , jj , jjgg , jjj , jk , jkl , jm , jmf , jmt , jmte , 
     &          jpp , jpz , JSKip , jy , k , karm , kk , kk6 , kkx , kmt
      INTEGER*4 knm , KSEq , kx , larm , lcc , lcou , LFL , LFL1 , 
     &          LFL2 , licz , lix , llx , lm , LMAx , LMAxe , lmh , 
     &          LNY , loc , loch , loct
      INTEGER*4 lp , LP1 , LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , 
     &          LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , lpit , lput , lpx , 
     &          lpxd , ls , lst
      INTEGER*4 luu , lx , LZEta , MAGa , MAGexc , MEMax , MEMx6 , 
     &          NANg , Ncall , NDIm , NEXpt , NICc , NLIft , nlin , 
     &          NMAx , NMAx1 , nowr , npoz , nrest , NSTart
      INTEGER*4 NSTop , NWR , nwyr , NYLde
      COMPLEX*16 ARM
      DIMENSION jmte(6) , prop(6) , Bten(1200)
      COMMON /CX    / NEXpt , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBdg(50) , VINf(50)
      COMMON /CEXC0 / NSTart(76) , NSTop(75)
      COMMON /CCC   / EG(50) , CC(50,5) , AGEli(50,200,2) , Q(3,200,8) , 
     &                NICc , NANg(200)
      COMMON /ILEWY / NWR
      COMMON /CH1T  / CHIs11
      COMMON /IGRAD / IGRd
      COMMON /LCZP  / EMH , INM , LFL1 , LFL2 , LFL
      COMMON /UWAGA / ITAk2
      COMMON /LEV   / TAU(75) , KSEq(500,4)
      COMMON /CCOUP / ZETa(50000) , LZEta(8)
      COMMON /KIN   / EPS(50) , EROot(50) , FIEx(50,2) , IEXp , IAXs(50)
      COMMON /YEXPT / YEXp(32,1500) , IY(1500,32) , CORf(1500,32) , 
     &                DYEx(32,1500) , NYLde(50,32) , UPL(32,50) , 
     &                YNRm(32,50) , IDRn , ILE(32)
      COMMON /COMME / ELM(500) , ELMu(500) , ELMl(500) , SA(500)
      COMMON /CLM   / LMAx
      COMMON /COEX  / EN(75) , SPIn(75) , ACCur , DIPol , ZPOl , ACCa , 
     &                ISO
      COMMON /CLCOM8/ CAT(600,3) , ISMax
      COMMON /AZ    / ARM(600,7)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /COEX2 / NMAx , NDIm , NMAx1
      COMMON /PTH   / IPAth(75) , MAGa(75)
      COMMON /PRT   / IPRm(20)
      COMMON /CEXC  / MAGexc , MEMax , LMAxe , MEMx6 , IVAr(500)
      COMMON /SKP   / JSKip(50)
      COMMON /LIFE  / NLIft
      COMMON /LOGY  / LNY , INTr , IPS1
      issp = 0
      Chilo = 0.
      fx = 2.*SPIn(1) + 1.
      Chisq = 0.
      LFL = 0
      chis1 = 0.
      ixx = NDIm*MEMax + LP11
      DO i1 = 1 , ixx
         ZETa(i1) = 0.
      ENDDO
      DO ii = 1 , LP6
         ILE(ii) = 1
      ENDDO
      itemp = 0
      NWR = 0
      iapx = 1
      DO jkl = 1 , NEXpt
         IEXp = jkl
         IGRd = 0
         LFL2 = 1
         IF ( ITAk2.EQ.-1 ) THEN
            DO larm = 1 , 4
               DO karm = 1 , LP10
                  ARM(karm,larm) = (0.,0.)
               ENDDO
            ENDDO
         ENDIF
         iflg = 0
         IF ( IEXp.NE.1 ) THEN
            kk = NANg(IEXp)
            DO jjj = 1 , LP6
               ILE(jjj) = ILE(jjj) + NYLde(IEXp-1,jjj)
            ENDDO
         ENDIF
         lp = 3
         IF ( JSKip(jkl).EQ.0 ) GOTO 200
         IF ( MAGa(IEXp).EQ.0 ) lp = 1
         IF ( Ncall.EQ.0 ) GOTO 150
         IF ( Icll.EQ.4 ) GOTO 100
 50      loch = LP3*(MEMax-1) + NMAx + LP11
         DO k = 1 , loch
            ZETa(k) = 0.
         ENDDO
         CALL LOAD(IEXp,1,2,0.D0,jj)
         DO k = 1 , LMAx
            fc = 2.
            IF ( k.EQ.LMAx ) fc = 1.
            IF ( DBLE(INT(SPIn(1))).LT.SPIn(1) ) fc = 2.
            loc = 0
            polm = DBLE(k-1) - SPIn(1)
            CALL LOAD(IEXp,3,2,polm,jj)
            CALL PATH(jj)
            CALL LOAD(IEXp,2,2,polm,jj)
            CALL APRAM(IEXp,1,1,jj,ACCa)
            IF ( Ncall.NE.0 ) THEN
               IF ( Icll.NE.3 ) THEN
                  DO indx = 1 , MEMx6
                     CALL APRAM(IEXp,0,indx,jj,ACCa)
                     kx = 0
                     DO i11 = 1 , NMAx
                        IF ( NSTart(i11).NE.0 ) THEN
                           loc = LP3*(indx-1) + i11 + LP11
                           jpp = INT(2.*SPIn(i11)+1.)
                           lpx = MIN(lp,jpp)
                           IF ( ISO.NE.0 ) lpx = NSTop(i11)
     &                          - NSTart(i11) + 1
                           DO lpxd = 1 , lpx
                              kx = kx + 1
                              ZETa(loc) = ZETa(loc) + fc*DBLE(ARM(kx,5))
     &                           *DBLE(ARM(kx,6))
     &                           /fx + fc*IMAG(ARM(kx,5))
     &                           *IMAG(ARM(kx,6))/fx
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
            CALL TENB(k,Bten,LMAx)
         ENDDO
         IF ( loc.NE.0 ) THEN
            REWIND 14
            WRITE (14,*) (ZETa(i11),i11=LP8,loch)
         ENDIF
         CALL TENS(Bten)
         IF ( Ncall.EQ.0 ) GOTO 200
         IF ( Icll.GE.2 ) GOTO 200
         llx = 28*NMAx
         DO lx = 1 , llx
            ZETa(LP9+lx) = ZETa(lx)
         ENDDO
         IF ( Icll.NE.1 ) GOTO 200
 100     iapx = 0
         issp = 1
         CALL LOAD(IEXp,1,1,0.D0,jj)
         CALL ALLOC(ACCur)
         CALL SNAKE(IEXp,ZPOl)
         CALL SETIN
         DO k = 1 , LMAx
            polm = DBLE(k-1) - SPIn(1)
            CALL LOAD(IEXp,2,1,polm,kk)
            IF ( IPRm(7).EQ.-1 ) WRITE (22,99001) polm , IEXp
99001       FORMAT (1X//40X,'EXCITATION AMPLITUDES'//10X,'M=',1F4.1,5X,
     &              'EXPERIMENT',1X,1I2//5X,'LEVEL',2X,'SPIN',2X,'M',5X,
     &              'REAL AMPLITUDE',2X,'IMAGINARY AMPLITUDE'//)
            CALL STING(kk)
            CALL PATH(kk)
            CALL INTG(IEXp)
            CALL TENB(k,Bten,LMAx)
            IF ( IPRm(7).EQ.-1 ) THEN
               DO j = 1 , ISMax
                  WRITE (22,99002) INT(CAT(j,1)) , CAT(j,2) , CAT(j,3) , 
     &                             DBLE(ARM(j,5)) , IMAG(ARM(j,5))
99002             FORMAT (7X,1I2,3X,1F4.1,2X,1F4.1,2X,1E14.6,2X,1E14.6)
               ENDDO
            ENDIF
         ENDDO
         CALL TENS(Bten)
         IF ( IPRm(7).EQ.-1 ) THEN
            DO jjgg = 2 , NMAx
               loct = (jjgg-1)*28 + 1
               WRITE (22,99003) jjgg , ZETa(loct)
99003          FORMAT (2X,'LEVEL',1X,1I2,10X,'POPULATION',1X,1E14.6)
            ENDDO
         ENDIF
         GOTO 200
 150     IF ( iflg.EQ.1 ) THEN
            itemp = 1
            iflg = 2
            GOTO 50
         ELSE
            IF ( iflg.EQ.2 ) GOTO 300
            itemp = 2
            iflg = 1
            GOTO 100
         ENDIF
 200     CALL CEGRY(Chisq,itemp,Chilo,Idr,nwyr,Icll,issp,0)
         issp = 0
         IF ( Ncall.EQ.0 .AND. JSKip(jkl).NE.0 ) THEN
            IF ( Ncall.NE.0 ) GOTO 200
            GOTO 150
         ELSE
            NWR = NWR + nwyr
            IF ( Icll.LE.2 .AND. JSKip(jkl).NE.0 ) THEN
               IF ( IEXp.EQ.1 ) chish = CHIs11
               IF ( Icll.EQ.1 ) chis1 = CHIs11
               IF ( Icll.EQ.0 ) chis1 = Chisq
               LFL2 = 0
               IGRd = 1
               IF ( ITAk2.EQ.-1 ) LFL = 1
               REWIND 14
               READ (14,*) (ZETa(i11),i11=LP8,loch)
               DO larm = 1 , 4
                  DO karm = 1 , LP10
                     ARM(karm,larm) = (0.,0.)
                  ENDDO
               ENDDO
               chisx = 0.
               llx = 28*NMAx
               DO lix = 1 , llx
                  ZETa(LP9+lix) = ZETa(lix)
               ENDDO
               CALL CEGRY(chisx,itemp,Chilo,Idr,nwyr,0,0,1)
               DO knm = 1 , MEMax
                  INM = knm
                  chisx = 0.
                  EMH = ELM(INM)
                  ELM(INM) = 1.05*EMH
                  lcc = LP3*(INM-1) + LP11
                  DO lst = 2 , NMAx
                     wz = ZETa(lst+lcc)
                     inpx = (lst-1)*28
                     DO jy = 1 , 4
                        inp = inpx + (jy-1)*7
                        IF ( jy.EQ.1 ) pr = ZETa(LP13+inp) + 1.E-12
                        jmf = 2*jy - 1
                        IF ( IAXs(IEXp).EQ.0 ) jmf = 1
                        DO jm = 1 , jmf
                           inp = inp + 1
                           ZETa(inp) = ZETa(inp+LP9)*(1.+.1*EMH*wz/pr)
                        ENDDO
                     ENDDO
                  ENDDO
                  CALL CEGRY(chisx,itemp,Chilo,Idr,nwyr,0,0,0)
                  ELM(INM) = EMH
               ENDDO
               IF ( ITAk2.EQ.-1 .AND. LFL1.NE.0 ) THEN
                  IF ( IPRm(17).NE.0 ) THEN
                     kmt = ABS(IPRm(17))
                     WRITE (22,99004) IEXp
99004                FORMAT (1X///20X,'EXPERIMENT',11X,1I2,5X,
     &                       'D(LOG(P))/D(LOG(ME)) MAP'/20X,52('-')///)
                     nlin = (NMAx-2)/6 + 1
                     nrest = NMAx - 1 - 6*(nlin-1)
                     DO ilin = 1 , nlin
                        npoz = 6
                        IF ( ilin.EQ.nlin ) npoz = nrest
                        inpo = (ilin-1)*6 + 2
                        inko = inpo + npoz - 1
                        lpit = 0
                        DO lm = inpo , inko
                           lpit = lpit + 1
                           jmte(lpit) = lm
                        ENDDO
                        WRITE (22,99005) (jmte(lm),lm=1,lpit)
99005                   FORMAT (5X,'LEVEL',6(8X,1I2,9X))
                        WRITE (22,99006)
     &                         (ZETa(LP13+(jpz-1)*28),jpz=inpo,inko)
99006                   FORMAT (1X,'EXC.PROB.',6(5X,1E10.4,4X))
                        DO jmt = 1 , kmt
                           lput = 0
                           DO ls = inpo , inko
                              lput = lput + 1
                              prop(lput) = 0.
                              DO lm = 1 , MEMx6
                                 inzz = ls + LP3*(lm-1) + LP11
                                 inzzz = LP13 + (ls-1)*28
                                 IF ( ABS(ZETa(inzzz)).LT.1.E-20 )
     &                                ZETa(inzzz) = 1.E-20
                                 val = 2.*ELM(lm)*ZETa(inzz)/ZETa(inzzz)
                                 aval = ABS(val)
                                 IF ( aval.GT.ABS(prop(lput)) ) THEN
                                    prop(lput) = val
                                    lmh = lm
                                    jmte(lput) = lm
                                 ENDIF
                              ENDDO
                              izzz = (lmh-1)*LP3 + LP11 + ls
                              ZETa(izzz) = 0.
                           ENDDO
                           WRITE (22,99007)
     &                            (jmte(lcou),prop(lcou),lcou=1,npoz)
99007                      FORMAT (10X,6(2X,'(',1X,1I3,1X,1E8.2,')',2X))
                        ENDDO
                     ENDDO
                     REWIND 14
                     READ (14,*) (ZETa(i11),i11=LP8,loch)
                     IF ( IPRm(17).LT.0 ) GOTO 300
                  ENDIF
                  LFL = 0
                  WRITE (22,99008) IEXp
99008             FORMAT (10X,'EXPERIMENT',1X,1I2/10X,
     &                    'D(LOG(Y)/D(LOG(ME))',//)
                  ile1 = ILE(1) + NYLde(IEXp,1) - 1
                  ile3 = ILE(1)
                  licz = 0
                  DO ile2 = ile3 , ile1
                     licz = licz + 1
                     idec = IY(ile2,1)
                     IF ( idec.GT.1000 ) idec = idec/1000
                     luu = 6*licz - 5
                     jk = (luu-1)/LP10 + 1
                     kk = luu - LP10*(jk-1)
                     kk6 = kk + 5
                     WRITE (22,99009) KSEq(idec,3) , KSEq(idec,4) , 
     &                                (INT(DBLE(ARM(kkx,jk))),
     &                                IMAG(ARM(kkx,jk)),kkx=kk,kk6)
99009                FORMAT (2X,1I2,'--',1I2,5X,
     &                       6('(',1I3,2X,1E8.2,')',3X))
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
 300  ENDDO
      IF ( ITAk2.EQ.-1 .AND. Icll.LT.2 ) ITAk2 = 0
      IF ( Ncall.NE.0 ) THEN
         IF ( Icll.LE.2 ) THEN
            IF ( Icll.EQ.1 ) CALL CEGRY(Chisq,itemp,Chilo,Idr,nowr,7,
     &                                  issp,0)
         ENDIF
         CALL BRANR(Chisq,NWR,Chilo)
         CALL MIXR(NWR,0,Chisq,Chilo)
         CALL CHMEM(NWR,Chisq,Chilo)
         NWR = NWR + NLIft
         Chisq = Chisq/NWR
         IF ( INTr.NE.0 ) THEN
            chx = Chisq
            Chisq = Chilo
            Chilo = chx
         ENDIF
      ENDIF
      END
