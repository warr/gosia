 
C----------------------------------------------------------------------
 
      SUBROUTINE ADHOC(Oph,Idr,Nfd,Ntap,Iyr)
      IMPLICIT NONE
      REAL*8 ACCa , ACCur , AGEli , BRAt , CC , CORf , DELta , DIPol , 
     &       DIX , DMIx , DMIxe , DYEx , EAMx , EG , EN , ENDec , ENZ , 
     &       EP , ODL , Q
      REAL*8 SPIn , TAU , TIMel , TLBdg , UPL , VINf , wamx , wbra , 
     &       wdl , wlf , XA , XA1 , YEXp , YGN , YGP , YNRm , ZPOl
      INTEGER*4 IAMx , IAMy , iax , IBRc , Idr , IDRn , iexp1 , IFMo , 
     &          ILE , ilft , IMIx , iosr , ipri , IPRm , ISO , isrt1 , 
     &          ITMa , ITS , iuf , IVAr
      INTEGER*4 IY , Iyr , IZ , IZ1 , jic , jicc , juf , KSEq , lb , 
     &          li , licc , LIFct , llia , LMAxe , lxt , MAGexc , MEM , 
     &          MEMax , MEMx6 , n1
      INTEGER*4 n2 , NAMx , NANg , NBRa , ndas , NDL , NDSt , ndtp , 
     &          NEXpt , Nfd , NICc , nistr , NLIft , ns1 , ns2 , ns3 , 
     &          ns4 , Ntap , nvare , NYLde
      CHARACTER*4 Oph
      COMMON /CCCDS / NDSt(50)
      COMMON /DIMX  / DIX(4) , ODL(200)
      COMMON /TRA   / DELta(500,3) , ENDec(500) , ITMa(50,200) , 
     &                ENZ(200)
      COMMON /LIFE  / NLIft
      COMMON /MIXD  / DMIxe(20,2) , DMIx(20) , IMIx(20) , NDL
      COMMON /ME2D  / EAMx(100,2) , NAMx , IAMx(100) , IAMy(100,2)
      COMMON /LIFE1 / LIFct(50) , TIMel(2,50)
      COMMON /BRNCH / BRAt(50,2) , IBRc(2,50) , NBRa
      COMMON /YEXPT / YEXp(32,1500) , IY(1500,32) , CORf(1500,32) , 
     &                DYEx(32,1500) , NYLde(50,32) , UPL(32,50) , 
     &                YNRm(32,50) , IDRn , ILE(32)
      COMMON /YTEOR / YGN(500) , YGP(500) , IFMo
      COMMON /LEV   / TAU(75) , KSEq(500,4)
      COMMON /CCC   / EG(50) , CC(50,5) , AGEli(50,200,2) , Q(3,200,8) , 
     &                NICc , NANg(200)
      COMMON /COEX  / EN(75) , SPIn(75) , ACCur , DIPol , ZPOl , ACCa , 
     &                ISO
      COMMON /CX    / NEXpt , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBdg(50) , VINf(50)
      COMMON /CEXC  / MAGexc , MEMax , LMAxe , MEMx6 , IVAr(500)
      COMMON /PRT   / IPRm(20)
      COMMON /TRB   / ITS
      iosr = 0
      READ * , IFMo
      READ * , NICc , nistr
      READ * , (EG(jicc),jicc=1,NICc)
      Iyr = 1
      DO jic = 1 , nistr
         READ * , isrt1
         IF ( isrt1.GT.6 ) isrt1 = isrt1 - 3
         READ * , (CC(jicc,isrt1),jicc=1,NICc)
      ENDDO
      READ * , (NANg(jicc),jicc=1,NEXpt)
      REWIND 9
      READ (9,*) Nfd
      DO jicc = 1 , Nfd
         READ (9,*) ODL(jicc)
         READ (9,*) ENZ(jicc)
         DO isrt1 = 1 , 8
            READ (9,*) (Q(licc,jicc,isrt1),licc=1,3)
         ENDDO
      ENDDO
      DO jic = 1 , NEXpt
         juf = NANg(jic)
         IF ( juf.LT.0 ) THEN
            juf = ABS(juf)
            DO jicc = 1 , juf
               AGEli(jic,jicc,1) = AGEli(jic-1,jicc,1)
               AGEli(jic,jicc,2) = AGEli(jic-1,jicc,2)
               ITMa(jic,jicc) = ITMa(jic-1,jicc)
            ENDDO
            IF ( Oph.NE.'GOSI' ) NANg(jic) = ABS(NANg(jic))
         ELSE
            READ * , (ITMa(jic,jicc),jicc=1,juf)
            READ * , (AGEli(jic,jicc,1),jicc=1,juf)
            READ * , (AGEli(jic,jicc,2),jicc=1,juf)
         ENDIF
      ENDDO
      CALL SEQ(Idr)
      DO jic = 1 , NEXpt
         juf = NANg(jic)
         juf = ABS(juf)
         DO jicc = 1 , juf
            DO lxt = 1 , 2
               AGEli(jic,jicc,lxt) = AGEli(jic,jicc,lxt)*.0174532925
            ENDDO
         ENDDO
      ENDDO
      TAU(1) = 1.E+25
      READ * , ns1 , ns2
      DO li = 1 , Idr
         IF ( KSEq(li,3).EQ.ns1 .AND. KSEq(li,4).EQ.ns2 ) GOTO 100
      ENDDO
 100  IDRn = li
      IF ( Oph.NE.'GOSI' ) RETURN
      DO li = 1 , NEXpt
         juf = NANg(li)
         IF ( juf.LT.0 ) THEN
            juf = ABS(juf)
            NANg(li) = juf
            NDSt(li) = NDSt(li-1)
            DO jicc = 1 , juf
               UPL(jicc,li) = UPL(jicc,li-1)
               YNRm(jicc,li) = YNRm(jicc,li-1)
            ENDDO
         ELSE
            READ * , NDSt(li)
            ndas = NDSt(li)
            READ * , (UPL(jicc,li),jicc=1,ndas)
            READ * , (YNRm(jicc,li),jicc=1,ndas)
         ENDIF
      ENDDO
      READ * , Ntap
      IF ( Ntap.NE.0 ) THEN
         ipri = IPRm(2)
         CALL READY(Idr,Ntap,ipri)
         ndtp = 0
         DO iexp1 = 1 , NEXpt
            juf = NDSt(iexp1)
            DO iuf = 1 , juf
               ndtp = ndtp + NYLde(iexp1,iuf)
            ENDDO
         ENDDO
         nvare = 0
         DO iexp1 = 1 , MEMax
            IF ( IVAr(iexp1).EQ.1 .OR. IVAr(iexp1).EQ.2 )
     &           nvare = nvare + 1
         ENDDO
         WRITE (22,99001) ndtp , nvare
99001    FORMAT (1X//5X,1I4,1X,'EXPERIMENTAL YIELDS',10X,1I3,1X,
     &           'MATRIX ELEMENTS TO BE VARIED'///)
      ENDIF
      READ * , NBRa , wbra
      IF ( ITS.EQ.2 ) THEN
         REWIND 18
         WRITE (18,*) MEMax
      ENDIF
      IF ( NBRa.NE.0 ) THEN
         WRITE (22,99002)
99002    FORMAT (40X,'BRANCHING RATIOS',//5X,'NS1',5X,'NF1',5X,'NS2',5X,
     &           'NF2',5X,'RATIO(1:2)',9X,'ERROR')
         DO lb = 1 , NBRa
            READ * , ns1 , ns2 , ns3 , ns4 , BRAt(lb,1) , BRAt(lb,2)
            BRAt(lb,2) = BRAt(lb,2)/(SQRT(wbra)+1.E-10)
            WRITE (22,99003) ns1 , ns2 , ns3 , ns4 , BRAt(lb,1) , 
     &                       BRAt(lb,2)
99003       FORMAT (5X,1I2,6X,1I2,6X,1I2,6X,1I2,5X,1F10.5,5X,1F10.5)
            DO li = 1 , Idr
               IF ( KSEq(li,3).EQ.ns3 .AND. KSEq(li,4).EQ.ns4 ) THEN
                  IBRc(2,lb) = li
               ELSEIF ( KSEq(li,3).EQ.ns1 .AND. KSEq(li,4).EQ.ns2 ) THEN
                  IBRc(1,lb) = li
               ENDIF
            ENDDO
            IF ( ITS.EQ.2 ) THEN
               n1 = IBRc(1,lb)
               n2 = IBRc(2,lb)
               WRITE (18,*) KSEq(n1,1) , KSEq(n2,1)
               WRITE (18,*) KSEq(n1,1) , KSEq(n2,2)
               WRITE (18,*) KSEq(n1,1) , KSEq(n1,2)
               WRITE (18,*) KSEq(n2,1) , KSEq(n1,2)
               WRITE (18,*) KSEq(n2,1) , KSEq(n2,2)
               IF ( KSEq(n1,2).NE.0 .AND. KSEq(n2,2).NE.0 ) WRITE (18,*)
     &              KSEq(n1,2) , KSEq(n2,2)
            ENDIF
         ENDDO
         WRITE (22,99004) wbra
99004    FORMAT (5X,'BRANCHING RATIOS ARE TAKEN WITH WEIGHT',2X,1E14.6)
      ENDIF
      READ * , NLIft , wlf
      IF ( NLIft.NE.0 ) THEN
         WRITE (22,99005)
99005    FORMAT (1X///30X,'LIFETIMES(PSEC)'///5X,'LEVEL',9X,'LIFETIME',
     &           5X,'ERROR'/)
         DO ilft = 1 , NLIft
            READ * , LIFct(ilft) , TIMel(1,ilft) , TIMel(2,ilft)
            TIMel(2,ilft) = TIMel(2,ilft)/(SQRT(wlf)+1.E-10)
            WRITE (22,99006) LIFct(ilft) , TIMel(1,ilft) , TIMel(2,ilft)
99006       FORMAT (7X,1I2,6X,1F10.2,3X,1F10.2)
         ENDDO
         WRITE (22,99007) wlf
99007    FORMAT (1X/10X,'LIFETIMES ARE TAKEN WITH WEIGHT',2X,1E14.6)
      ENDIF
      READ * , NDL , wdl
      IF ( NDL.NE.0 ) THEN
         WRITE (22,99008)
99008    FORMAT (1X//20X,'EXPERIMENTAL E2/M1 MIXING RATIOS'///10X,
     &           'TRANSITION',12X,'DELTA',10X,'ERROR'/)
         DO li = 1 , NDL
            READ * , ns1 , ns2 , DMIxe(li,1) , DMIxe(li,2)
            DMIxe(li,2) = DMIxe(li,2)/(SQRT(wdl)+1.E-10)
            WRITE (22,99012) ns1 , ns2 , DMIxe(li,1) , DMIxe(li,2)
            DO lb = 1 , Idr
               IF ( KSEq(lb,3).EQ.ns1 .AND. KSEq(lb,4).EQ.ns2 ) THEN
                  IMIx(li) = lb
                  DMIx(li) = .8326*(EN(ns1)-EN(ns2))
                  IF ( ITS.EQ.2 ) WRITE (18,*) KSEq(lb,1) , KSEq(lb,2)
               ENDIF
            ENDDO
         ENDDO
         WRITE (22,99009) wdl
99009    FORMAT (/10X,'E2/M1 MIXING RATIOS ARE TAKEN WITH WEIGHT',2X,
     &           1E14.6)
      ENDIF
      IF ( ITS.EQ.2 ) WRITE (18,*) iosr , iosr
      READ * , NAMx , wamx
      IF ( NAMx.EQ.0 ) RETURN
      WRITE (22,99010)
99010 FORMAT (1X//30X,'EXPERIMENTAL MATRIX ELEMENT(S)'///10X,
     &        'TRANSITION',10X,'MAT.EL.',10X,'ERROR'/)
      DO iax = 1 , NAMx
         READ * , llia , ns1 , ns2 , EAMx(iax,1) , EAMx(iax,2)
         IAMy(iax,1) = ns1
         IAMy(iax,2) = ns2
         EAMx(iax,2) = EAMx(iax,2)/(SQRT(wamx)+1.E-10)
         WRITE (22,99012) ns1 , ns2 , EAMx(iax,1) , EAMx(iax,2)
         IAMx(iax) = MEM(ns1,ns2,llia)
      ENDDO
      WRITE (22,99011) wamx
99011 FORMAT (/10X,' MATRIX ELEMENT(S) ARE TAKEN WITH WEIGHT',2X,1E14.6)
99012 FORMAT (10X,1I2,'---',1I2,14X,1F9.4,8X,1F9.4)
      END
