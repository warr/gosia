 
C----------------------------------------------------------------------
 
      SUBROUTINE CEGRY(Chisq,Itemp,Chilo,Idr,Nwyr,Icall,Issp,Iredv)
      IMPLICIT NONE
      REAL*8 ACCa , ACCur , AGEli , AKS , BETar , CC , ccc , ccd , 
     &       Chilo , Chisq , CNOr , cnr , cocos , CORf , d , decen , 
     &       DELta , DEV , DIPol , DIX
      REAL*8 dl , DQ , DSIgs , DYEx , effi , EG , EMH , EN , ENDec , 
     &       ENZ , EP , EPS , EROot , fi0 , fi1 , fic , FIEx , figl , 
     &       fm , g
      REAL*8 gth , ODL , part , partl , Q , QCEn , rik , rl , rx , ry , 
     &       rys , rz , sf , sgm , SGW , SPIn , SUBch1 , SUBch2 , sum3 , 
     &       SUMcl
      REAL*8 sumpr , TACOS , TAU , TETacm , tetrc , tfac , thc , TLBdg , 
     &       TREp , UPL , VACdp , VINf , wf , XA , XA1 , XNOr , YEXp , 
     &       YGN , YGP , YNRm
      REAL*8 ZPOl
      INTEGER*4 iabc , IAXs , IBYp , Icall , ICLust , id , idc , Idr , 
     &          IDRn , IEXp , ifdu , IFMo , ifxd , IGRd , ii , ILE , 
     &          ile2 , IMIn , inclus , INM
      INTEGER*4 INNr , ipd , IPRm , IRAwex , Iredv , ISO , Issp , 
     &          Itemp , ITMa , ITS , iva , iw , IWF , ixl , ixm , IY , 
     &          iyex , IZ , IZ1 , jj
      INTEGER*4 jj1 , jk , jpc , JSKip , k , k9 , kc , kj , kk , KSEq , 
     &          KVAr , l , l1 , LAStcl , LFL , LFL1 , LFL2 , lic , 
     &          licz , ll1
      INTEGER*4 LNOrm , LP1 , LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , 
     &          LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , lth , lu , luu , 
     &          na , NANg , NDIm
      INTEGER*4 NDSt , NEXpt , nf , nf1 , ni , ni1 , NICc , NLIft , 
     &          NMAx , NMAx1 , Nwyr , NYLde
      CHARACTER*4 wupl , war
      DIMENSION part(32,50,2) , lic(32) , lth(500) , cnr(32,50) , 
     &          partl(32,50,2)
      COMMON /CLUST / ICLust(50,200) , LAStcl(50,20) , SUMcl(20,500) , 
     &                IRAwex(50)
      COMMON /ODCH  / DEV(500)
      COMMON /COEX2 / NMAx , NDIm , NMAx1
      COMMON /TRA   / DELta(500,3) , ENDec(500) , ITMa(50,200) , 
     &                ENZ(200)
      COMMON /BREC  / BETar(50)
      COMMON /DIMX  / DIX(4) , ODL(200)
      COMMON /VAC   / VACdp(3,75) , QCEn , DQ , XNOr , AKS(6,75) , IBYp
      COMMON /CINIT / CNOr(32,75) , INNr
      COMMON /PRT   / IPRm(20)
      COMMON /LIFE  / NLIft
      COMMON /LEV   / TAU(75) , KSEq(500,4)
      COMMON /IGRAD / IGRd
      COMMON /CX    / NEXpt , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBdg(50) , VINf(50)
      COMMON /MINNI / IMIn , LNOrm(50)
      COMMON /LCZP  / EMH , INM , LFL1 , LFL2 , LFL
      COMMON /YTEOR / YGN(500) , YGP(500) , IFMo
      COMMON /SEL   / KVAr(500)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /CCC   / EG(50) , CC(50,5) , AGEli(50,200,2) , Q(3,200,8) , 
     &                NICc , NANg(200)
      COMMON /YEXPT / YEXp(32,1500) , IY(1500,32) , CORf(1500,32) , 
     &                DYEx(32,1500) , NYLde(50,32) , UPL(32,50) , 
     &                YNRm(32,50) , IDRn , ILE(32)
      COMMON /KIN   / EPS(50) , EROot(50) , FIEx(50,2) , IEXp , IAXs(50)
      COMMON /WARN  / SGW , SUBch1 , SUBch2 , IWF
      COMMON /COEX  / EN(75) , SPIn(75) , ACCur , DIPol , ZPOl , ACCa , 
     &                ISO
      COMMON /SKP   / JSKip(50)
      COMMON /TRB   / ITS
      COMMON /TCM   / TETacm(50) , TREp(50) , DSIgs(50)
      COMMON /CCCDS / NDSt(50)
      ifxd = 0
      tetrc = TREp(IEXp)
      IF ( Icall.EQ.4 .AND. IPRm(13).EQ.-2 ) THEN
         IPRm(13) = 0
         WRITE (22,99001)
99001    FORMAT (1X//20X,'NORMALIZATION CONSTANTS'//2X,'EXPERIMENT',5X,
     &           'DETECTORS(1-32)')
         DO jpc = 1 , NEXpt
            k = NDSt(jpc)
            WRITE (22,99012) jpc , (CNOr(l,jpc),l=1,k)
         ENDDO
         WRITE (22,99002)
99002    FORMAT (1X//20X,'RECOMMENDED RELATIVE GE(LI) EFFICIENCIES'//2X,
     &           'EXPERIMENT')
         DO jpc = 1 , NEXpt
            IF ( ABS(cnr(1,jpc)).LT.1.E-9 ) cnr(1,jpc) = 1.
            k = NDSt(jpc)
            WRITE (22,99012) jpc , (cnr(l,jpc)/cnr(1,jpc),l=1,k)
         ENDDO
      ENDIF
      DO jpc = 1 , LP6
         lic(jpc) = 0
      ENDDO
      IF ( Icall.NE.7 ) THEN
         IF ( Itemp.EQ.0 ) THEN
            Nwyr = 0
            IF ( IGRd.NE.1 ) THEN
               IF ( IEXp.EQ.1 ) sumpr = 0.
               IF ( IEXp.EQ.1 ) sum3 = 0.
               DO jj = 1 , LP6
                  DO jk = 1 , 2
                     partl(jj,IEXp,jk) = 0.
                     part(jj,IEXp,jk) = 0.
                  ENDDO
               ENDDO
            ENDIF
            CALL DECAY(Chisq,NLIft,Chilo)
            IF ( Icall.EQ.4 .AND. IPRm(14).EQ.-1 ) THEN
               IF ( IEXp.EQ.NEXpt ) IPRm(14) = 0
               WRITE (22,99003)
99003          FORMAT (1X//20X,'VACUUM DEPOLARIZATION COEFFICIENTS '//)
               WRITE (22,99004) IEXp
99004          FORMAT (5X,'EXPERIMENT',1X,1I2/5X,'LEVEL',10X,'G2',10X,
     &                 'G4',10X,'G6'/)
               DO iva = 2 , NMAx
                  WRITE (22,99005) iva , (VACdp(ii,iva),ii=1,3)
99005             FORMAT (7X,1I2,9X,3(1F6.4,6X))
               ENDDO
            ENDIF
            fi0 = FIEx(IEXp,1)
            fi1 = FIEx(IEXp,2)
            na = NANg(IEXp)
            DO k = 1 , LP2
               DO k9 = 1 , 20
                  SUMcl(k9,k) = 0.
               ENDDO
            ENDDO
            k9 = 0
            DO k = 1 , na
               gth = AGEli(IEXp,k,1)
               figl = AGEli(IEXp,k,2)
               ifxd = 0
               fm = (fi0+fi1)/2.
               IF ( Icall.EQ.4 ) ifxd = 1
               CALL ANGULA(YGN,Idr,ifxd,fi0,fi1,tetrc,gth,figl,k)
               IF ( IFMo.NE.0 ) THEN
                  id = ITMa(IEXp,k)
                  d = ODL(id)
                  rx = d*SIN(gth)*COS(figl-fm) - .25*SIN(tetrc)*COS(fm)
                  ry = d*SIN(gth)*SIN(figl-fm) - .25*SIN(tetrc)*SIN(fm)
                  rz = d*COS(gth) - .25*COS(tetrc)
                  rl = SQRT(rx*rx+ry*ry+rz*rz)
                  sf = d*d/rl/rl
                  thc = TACOS(rz/rl)
                  fic = ATAN2(ry,rx)
                  CALL ANGULA(YGP,Idr,ifxd,fi0,fi1,tetrc,thc,fic,k)
                  DO ixl = 1 , Idr
                     ixm = KSEq(ixl,3)
                     tfac = TAU(ixm)
                     YGN(ixl) = YGN(ixl) + .01199182*tfac*BETar(IEXp)
     &                          *(sf*YGP(ixl)-YGN(ixl))
                  ENDDO
               ENDIF
               IF ( IRAwex(IEXp).NE.0 ) THEN
                  ipd = ITMa(IEXp,k)
                  DO l = 1 , Idr
                     decen = ENDec(l)
                     cocos = SIN(tetrc)*SIN(gth)*COS(fm-figl)
     &                       + COS(tetrc)*COS(gth)
                     decen = decen*(1.+BETar(IEXp)*cocos)
                     CALL EFFIX(ipd,decen,effi)
                     YGN(l) = YGN(l)*effi
                  ENDDO
                  inclus = ICLust(IEXp,k)
                  IF ( inclus.NE.0 ) THEN
                     DO l = 1 , Idr
                        SUMcl(inclus,l) = SUMcl(inclus,l) + YGN(l)
                     ENDDO
                     IF ( k.NE.LAStcl(IEXp,inclus) ) GOTO 20
                     DO l = 1 , Idr
                        YGN(l) = SUMcl(inclus,l)
                     ENDDO
                  ENDIF
               ENDIF
               k9 = k9 + 1
               IF ( Icall.EQ.4 .AND. IPRm(8).EQ.-2 ) THEN
                  WRITE (22,99006) IEXp , k9
99006             FORMAT (1X//5X,
     &                 'CALCULATED AND EXPERIMENTAL YIELDS   EXPERIMENT'
     &                 ,1X,1I2,1X,'DETECTOR',1X,1I2//6X,'NI',5X,'NF',7X,
     &                 'II',8X,'IF',9X,'ENERGY(MEV)',6X,'YCAL',8X,
     &                 'YEXP',7X,'PC. DIFF.',2X,'(YE-YC)/SIGMA')
               ENDIF
               lu = ILE(k9)
               DO iabc = 1 , LP2
                  lth(iabc) = 0
               ENDDO
               DO l = 1 , Idr
                  ni = KSEq(l,3)
                  nf = KSEq(l,4)
                  IF ( l.EQ.IY(lu,k9) .OR. l.EQ.(IY(lu,k9)/1000) ) THEN
                     ifdu = 0
                     lic(k9) = lic(k9) + 1
                     licz = lic(k9)
                     Nwyr = Nwyr + 1
                     wf = CORf(lu,k9)
                     IF ( Icall.EQ.4 ) wf = 1.
                     IF ( Icall.EQ.1 .AND. Issp.EQ.1 ) wf = 1.
                     IF ( IY(lu,k9).GE.1000 ) THEN
                        ifdu = 1
                        l1 = IY(lu,k9)/1000
                        l1 = IY(lu,k9) - 1000*l1
                        YGN(l) = YGN(l) + YGN(l1)
                        lth(l1) = 1
                        IF ( Icall.EQ.4 .AND. IPRm(8).EQ.-2 ) THEN
                           war = '    '
                           sgm = (YEXp(k9,lu)-YGN(l)*CNOr(k9,IEXp))
     &                           /DYEx(k9,lu)
                           ni1 = KSEq(l1,3)
                           nf1 = KSEq(l1,4)
                           WRITE (22,99007) ni , ni1 , nf , nf1 , 
     &                            SPIn(ni) , SPIn(ni1) , SPIn(nf) , 
     &                            SPIn(nf1) , ENDec(l) , ENDec(l1) , 
     &                            YGN(l)*CNOr(k9,IEXp) , YEXp(k9,lu) , 
     &                            100.*(YEXp(k9,lu)-YGN(l)*CNOr(k9,IEXp)
     &                            )/YEXp(k9,lu) , sgm , war
99007                      FORMAT (4X,1I2,'+',1I2,'--',1I2,'+',1I2,3X,
     &                             1F4.1,'+',1F4.1,'--',1F4.1,'+',1F4.1,
     &                             3X,1F6.4,'+',1F6.4,2X,1E9.4,6X,1E9.4,
     &                             3X,1F6.1,5X,1F4.1,10X,1A4)
                           SUBch1 = SUBch1 + sgm*sgm
                        ENDIF
                     ENDIF
                     ry = YGN(l)*wf*CNOr(k9,IEXp) - YEXp(k9,lu)
                     IF ( ifdu.NE.1 ) THEN
                        IF ( Icall.EQ.4 .AND. IPRm(8).EQ.-2 ) THEN
                           war = '    '
                           sgm = (YEXp(k9,lu)-YGN(l)*CNOr(k9,IEXp))
     &                           /DYEx(k9,lu)
                           WRITE (22,99013) ni , nf , SPIn(ni) , 
     &                            SPIn(nf) , ENDec(l) , YGN(l)
     &                            *CNOr(k9,IEXp) , YEXp(k9,lu) , 
     &                            100.*(YEXp(k9,lu)-YGN(l)*CNOr(k9,IEXp)
     &                            )/YEXp(k9,lu) , sgm , war
                           SUBch1 = SUBch1 + sgm*sgm
                        ENDIF
                     ENDIF
                     rys = ry*ry
                     IF ( IGRd.EQ.1 ) Chisq = Chisq + rys/DYEx(k9,lu)
     &                    /DYEx(k9,lu)
                     IF ( k9.EQ.1 .AND. Iredv.EQ.1 ) DEV(licz) = ry
                     IF ( Iredv.NE.1 ) THEN
                        IF ( LFL.EQ.1 ) THEN
                           IF ( k9.EQ.1 ) THEN
                              luu = 6*licz - 5
                              jk = (luu-1)/LP10 + 1
                              kk = luu - LP10*(jk-1)
                              rik = DEV(licz) + YEXp(k9,lu)
                              sgm = -DEV(licz)/DYEx(k9,lu)
                              IF ( ITS.EQ.1 .AND. KVAr(INM).NE.0 )
     &                             WRITE (17,*) ni , nf , sgm , YGN(l)
     &                             *CNOr(k9,IEXp)/DYEx(k9,lu)
                              IF ( ITS.EQ.1 .AND. INM.EQ.1 )
     &                             WRITE (15,*) IEXp , rik/CNOr(1,IEXp)
     &                             , CNOr(1,IEXp) , DYEx(k9,lu) , 
     &                             YEXp(k9,lu)
                              CALL SIXEL(rik,ry,EMH,jk,kk,INM,licz)
                           ENDIF
                        ENDIF
                     ENDIF
                     IF ( IGRd.NE.1 ) THEN
                        IF ( JSKip(IEXp).NE.0 ) THEN
                           dl = DYEx(k9,lu)*DYEx(k9,lu)
                           part(k9,IEXp,1) = part(k9,IEXp,1) + YGN(l)
     &                        *YGN(l)*wf*wf/dl
                           part(k9,IEXp,2) = part(k9,IEXp,2) - 2.*YGN(l)
     &                        *wf*YEXp(k9,lu)/dl
                           sumpr = sumpr + YEXp(k9,lu)*YEXp(k9,lu)/dl
                           partl(k9,IEXp,1) = partl(k9,IEXp,1)
     &                        + YEXp(k9,lu)*YEXp(k9,lu)/dl
                           partl(k9,IEXp,2) = partl(k9,IEXp,2)
     &                        + LOG(wf*YGN(l)/YEXp(k9,lu))*YEXp(k9,lu)
     &                        *YEXp(k9,lu)/dl
                           sum3 = sum3 + YEXp(k9,lu)*YEXp(k9,lu)
     &                            *LOG(wf*YGN(l)/YEXp(k9,lu))**2/dl
                        ENDIF
                     ENDIF
                     lu = lu + 1
                  ELSE
                     IF ( JSKip(IEXp).EQ.0 ) YGN(IDRn) = 1.E+10
                     ry = YGN(l)/YGN(IDRn)
                     IF ( Icall.EQ.4 .AND. IPRm(8).EQ.-2 ) THEN
                        wupl = '    '
                        IF ( ry.GT.UPL(k9,IEXp) .AND. lth(l).EQ.0 )
     &                       wupl = 'UPL!'
                        IF ( IPRm(16).NE.0 .OR. wupl.NE.'    ' ) THEN
                           IF ( wupl.EQ.'    ' ) WRITE (22,99008) ni , 
     &                          nf , SPIn(ni) , SPIn(nf) , ENDec(l) , 
     &                          YGN(l)*CNOr(k9,IEXp) , wupl
99008                      FORMAT (6X,1I2,5X,1I2,7X,1F4.1,6X,1F4.1,9X,
     &                             1F6.4,6X,1E9.4,10X,1A4)
                           IF ( wupl.NE.'    ' ) THEN
                              sgm = (ry-UPL(k9,IEXp))/UPL(k9,IEXp)
                              WRITE (22,99013) ni , nf , SPIn(ni) , 
     &                               SPIn(nf) , ENDec(l) , YGN(l)
     &                               *CNOr(k9,IEXp) , UPL(k9,IEXp)
     &                               *CNOr(k9,IEXp)*YGN(IDRn) , 
     &                               100.*(1.-YGN(l)/UPL(k9,IEXp)
     &                               /YGN(IDRn)) , sgm , wupl
                              SUBch1 = SUBch1 + sgm*sgm
                           ENDIF
                        ENDIF
                     ENDIF
                     IF ( ry.GE.UPL(k9,IEXp) .AND. lth(l).NE.1 ) THEN
                        Chisq = Chisq + (ry-UPL(k9,IEXp))
     &                          *(ry-UPL(k9,IEXp))/UPL(k9,IEXp)
     &                          /UPL(k9,IEXp)
                        Chilo = Chilo + LOG(ry/UPL(k9,IEXp))**2
                        IF ( IWF.NE.0 ) THEN
                           WRITE (22,99009) IEXp , ni , nf , 
     &                            ry/UPL(k9,IEXp)
99009                      FORMAT (5X,'WARNINIG-EXP.',1I2,2X,'TRANS. ',
     &                             1I2,'--',1I2,5X,
     &                             'EXCEEDS UPPER LIMIT (RATIO=',1E14.6,
     &                             ')')
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
               IF ( IEXp.EQ.NEXpt ) IWF = 0
               IF ( Icall.EQ.4 .AND. IPRm(8).EQ.-2 ) THEN
                  WRITE (22,99010) SUBch1 - SUBch2
99010             FORMAT (1X/50X,'CHISQ SUBTOTAL = ',E14.6)
                  SUBch2 = SUBch1
               ENDIF
 20         ENDDO
            IF ( IGRd.EQ.1 ) RETURN
            IF ( IEXp.NE.NEXpt ) RETURN
            IF ( Icall.EQ.1 ) RETURN
         ELSE
            ifxd = 1
            IF ( Itemp.NE.2 ) ifxd = 0
            Nwyr = 1
            CALL DECAY(ccd,0,ccc)
            fi0 = FIEx(IEXp,1)
            fi1 = FIEx(IEXp,2)
            na = NANg(IEXp)
            DO k = 1 , LP2
               DO kj = 1 , 20
                  SUMcl(kj,k) = 0
               ENDDO
            ENDDO
            k9 = 0
            DO k = 1 , na
               gth = AGEli(IEXp,k,1)
               figl = AGEli(IEXp,k,2)
               fm = (fi0+fi1)/2.
               CALL ANGULA(YGN,Idr,ifxd,fi0,fi1,tetrc,gth,figl,k)
               IF ( IFMo.NE.0 ) THEN
                  id = ITMa(IEXp,k)
                  d = ODL(id)
                  rx = d*SIN(gth)*COS(figl-fm) - .25*SIN(tetrc)*COS(fm)
                  ry = d*SIN(gth)*SIN(figl-fm) - .25*SIN(tetrc)*SIN(fm)
                  rz = d*COS(gth) - .25*COS(tetrc)
                  rl = SQRT(rx*rx+ry*ry+rz*rz)
                  sf = d*d/rl/rl
                  thc = TACOS(rz/rl)
                  fic = ATAN2(ry,rx)
                  CALL ANGULA(YGP,Idr,ifxd,fi0,fi1,tetrc,thc,fic,k)
                  DO ixl = 1 , Idr
                     ixm = KSEq(ixl,3)
                     tfac = TAU(ixm)
                     IF ( tfac.GT.1.E+4 ) GOTO 25
                     YGN(ixl) = YGN(ixl) + .01199182*tfac*BETar(IEXp)
     &                          *(sf*YGP(ixl)-YGN(ixl))
                  ENDDO
 25               IFMo = 0
                  WRITE (22,99011)
99011             FORMAT (1X,/,2X,'DURING THE MINIMIZATION',1X,
     &    'IT WAS NECESSARY TO SWITCH OFF THE TIME-OF-FLIGHT CORRECTION'
     &    )
               ENDIF
               IF ( IRAwex(IEXp).NE.0 ) THEN
                  ipd = ITMa(IEXp,k)
                  DO l = 1 , Idr
                     decen = ENDec(l)
                     cocos = SIN(tetrc)*SIN(gth)*COS(fm-figl)
     &                       + COS(tetrc)*COS(gth)
                     decen = decen*(1.+BETar(IEXp)*cocos)
                     CALL EFFIX(ipd,decen,effi)
                     YGN(l) = YGN(l)*effi
                  ENDDO
                  inclus = ICLust(IEXp,k)
                  IF ( inclus.NE.0 ) THEN
                     DO l = 1 , Idr
                        SUMcl(inclus,l) = SUMcl(inclus,l) + YGN(l)
                     ENDDO
                     IF ( k.NE.LAStcl(IEXp,inclus) ) GOTO 40
                     DO l = 1 , Idr
                        YGN(l) = SUMcl(inclus,l)
                     ENDDO
                  ENDIF
               ENDIF
               k9 = k9 + 1
               iyex = NYLde(IEXp,k9) + ILE(k9) - 1
               ile2 = ILE(k9)
               DO l = ile2 , iyex
                  IF ( JSKip(IEXp).NE.0 ) THEN
                     idc = IY(l,k9)
                     IF ( idc.GE.1000 ) THEN
                        idc = idc/1000
                        ll1 = IY(l,k9) - idc*1000
                        YGN(idc) = YGN(idc) + YGN(ll1)
                     ENDIF
                     IF ( Itemp.EQ.1 ) THEN
                        CORf(l,k9) = CORf(l,k9)/(YGN(idc)+1.E-24)
                     ELSE
                        CORf(l,k9) = YGN(idc)
                        IF ( IMIn.LE.1 .AND. l.EQ.iyex ) CNOr(k9,IEXp)
     &                       = YEXp(k9,l)/YGN(idc)
                     ENDIF
                  ENDIF
               ENDDO
 40         ENDDO
            RETURN
         ENDIF
      ENDIF
      DO jj = 1 , NEXpt
         IF ( JSKip(jj).NE.0 ) THEN
            kc = NDSt(jj)
            DO jk = 1 , kc
               cnr(jk,jj) = -.5*part(jk,jj,2)/part(jk,jj,1)
               IF ( INNr.NE.0 ) CNOr(jk,jj) = cnr(jk,jj)
            ENDDO
            IF ( INNr.NE.1 ) THEN
               d = 0.
               g = 0.
               DO jj1 = jj , NEXpt
                  IF ( LNOrm(jj1).EQ.jj ) THEN
                     k = NDSt(jj1)
                     DO jk = 1 , k
                        d = d + YNRm(jk,jj1)*part(jk,jj1,1)*YNRm(jk,jj1)
                        g = g - .5*YNRm(jk,jj1)*part(jk,jj1,2)
                     ENDDO
                  ENDIF
               ENDDO
               IF ( LNOrm(jj).EQ.jj ) THEN
                  CNOr(1,jj) = g*YNRm(1,jj)/d
                  k = NDSt(jj)
                  IF ( k.NE.1 ) THEN
                     DO jk = 2 , k
                        CNOr(jk,jj) = YNRm(jk,jj)*CNOr(1,jj)/YNRm(1,jj)
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      IF ( INNr.NE.1 ) THEN
         DO jj = 1 , NEXpt
            IF ( LNOrm(jj).NE.jj ) THEN
               iw = LNOrm(jj)
               k = NDSt(jj)
               DO jk = 1 , k
                  CNOr(jk,jj) = CNOr(1,iw)*YNRm(jk,jj)/YNRm(1,iw)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      IF ( Icall.EQ.7 ) Chisq = 0.
      DO jj = 1 , NEXpt
         k = NDSt(jj)
         DO jk = 1 , k
            Chilo = Chilo + partl(jk,jj,1)*LOG(CNOr(jk,jj))
     &              **2 + partl(jk,jj,2)*2.*LOG(CNOr(jk,jj))
            Chisq = Chisq + CNOr(jk,jj)*CNOr(jk,jj)*part(jk,jj,1)
     &              + CNOr(jk,jj)*part(jk,jj,2)
         ENDDO
      ENDDO
      Chisq = Chisq + sumpr
      Chilo = Chilo + sum3
      RETURN
99012 FORMAT (1X,1I2,2X,32(1E8.2,1X))
99013 FORMAT (6X,1I2,5X,1I2,7X,1F4.1,6X,1F4.1,9X,1F6.4,6X,1E9.4,6X,
     &        1E9.4,3X,1F6.1,5X,1F4.1,10X,1A4)
      END
