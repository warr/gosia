 
C----------------------------------------------------------------------
 
      SUBROUTINE KLOPOT(K,Rlr)
      IMPLICIT NONE
      REAL*8 a , al , al1 , b , c , ch , CORf , d , dy , DYEx , e , 
     &       ELM , ELMl , ELMu , EP , g , g1 , g2 , rl , Rlr
      REAL*8 SA , sgm , TLBdg , u , umm , ump , UPL , ux , VINf , XA , 
     &       XA1 , YEXp , YNRm , ZETa
      INTEGER*4 i , IDRn , iex , iexh , iexp , ILE , indx , inh , ipf , 
     &          IVAr , IY , IZ , IZ1 , j , jm , jp , K , KVAr , l , lc
      INTEGER*4 ll , LMAxe , lngt , loc , LP1 , LP10 , LP11 , LP12 , 
     &          LP13 , LP14 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          lu , LZEta , MAGexc
      INTEGER*4 MEMax , MEMx6 , NEXpt , nf , ni , nm , np , NYLde
      COMMON /COMME / ELM(500) , ELMu(500) , ELMl(500) , SA(500)
      COMMON /YEXPT / YEXp(32,1500) , IY(1500,32) , CORf(1500,32) , 
     &                DYEx(32,1500) , NYLde(50,32) , UPL(32,50) , 
     &                YNRm(32,50) , IDRn , ILE(32)
      COMMON /CEXC  / MAGexc , MEMax , LMAxe , MEMx6 , IVAr(500)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /CCOUP / ZETa(50000) , LZEta(8)
      COMMON /CX    / NEXpt , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBdg(50) , VINf(50)
      COMMON /SEL   / KVAr(500)
      REWIND 14
      REWIND 18
      ipf = 1
      lngt = 0
      indx = 1
      REWIND 15
      REWIND 17
      DO i = 1 , MEMax
         ELM(i) = 0.
         ELMu(i) = 0.
         ELMl(i) = 0.
      ENDDO
      iexh = 1
 100  g = 0.
      d = 0.
 200  READ (15,*) iex , a , b , c , e
      IF ( iex.NE.iexh ) THEN
         EP(iexh) = g/d
         TLBdg(iexh) = g
         VINf(iexh) = d
         iexh = iex
         BACKSPACE 15
         IF ( iex.NE.0 ) GOTO 100
         REWIND 15
         iexp = 1
      ELSE
         g = g + e*a/c/c
         d = d + a*a/c/c
         lngt = lngt + 1
         GOTO 200
      ENDIF
 300  g1 = 0.
      g2 = 0.
      inh = indx
      iexh = iexp
 400  READ (18,*) lu , indx , iexp , al
      IF ( indx.NE.0 ) THEN
         READ (17,*) ni , nf , sgm , al1
         READ (15,*) iex , a , b , c , e
         IF ( iexp.NE.1 .AND. ipf.NE.1 ) THEN
 420        READ (15,*) iex , a , b , c , e
            IF ( iexp.NE.iex ) GOTO 420
            ipf = 1
         ENDIF
         IF ( indx.EQ.inh ) THEN
            dy = al*al1/b/c
            g1 = e*dy + g1
            g2 = -2.*dy*a + g2
            WRITE (14,*) indx , iexp , ni , nf , dy , a , e , c
            GOTO 400
         ENDIF
      ENDIF
      loc = (iexh-1)*LP2 + inh
      ipf = 0
      ZETa(loc) = (VINf(iexh)*g1+TLBdg(iexh)*g2)/VINf(iexh)/VINf(iexh)
      inh = indx
      REWIND 15
      BACKSPACE 17
      BACKSPACE 18
      IF ( indx.NE.0 ) GOTO 300
      WRITE (14,*) indx , iexp , ni , nf , dy , a , e , b
      REWIND 14
      REWIND 17
 500  READ (14,*) indx , iexp , ni , nf , dy , a , e , b
      IF ( indx.EQ.0 ) THEN
         WRITE (17,*) indx , iexp , sgm , ni , nf , u
         REWIND 17
         ll = 0
         ch = 0.
 550     READ (17,*) indx , iexp , sgm , ni , nf , u
         IF ( indx.EQ.0 ) THEN
            WRITE (22,99001)
99001       FORMAT (2X////40X,'TROUBLESHOOTING ROUT',
     &              'INE HAS BEEN ACTIVATED...'//5X,
     &              'LOCAL MINIMUM ANALYSIS FOLLOWS:'//)
            WRITE (22,99002) ch/ll
99002       FORMAT (2X//5X,'CHISQ FOR FIRST GE(LI)S ONLY ',
     &              'WITH INDEPENDENT NORMALIZATION=',1E12.4//5X,
     &              'NORM.CONSTANTS:'//)
            DO i = 1 , NEXpt
               WRITE (22,99003) i , EP(i)
99003          FORMAT (5X,'EXP.',1X,1I2,5X,'C=',1E14.6)
            ENDDO
            WRITE (22,99004)
99004       FORMAT (1X//5X,'M.E.',20X,'RL',20X,'STRENGTH',//)
            DO i = 1 , MEMax
               IF ( KVAr(i).NE.0 ) THEN
                  rl = LOG10(ELMu(i)/ABS(ELM(i)))
                  IF ( rl.GE.Rlr ) ELMl(i) = 1.
                  WRITE (22,99005) i , rl , ELMu(i)/lngt
99005             FORMAT (6X,1I3,18X,1F4.1,20X,1E7.2)
               ENDIF
            ENDDO
            WRITE (22,99006)
99006       FORMAT (2X////40X,'ANALYSIS OF SIGNIFICANT DEPENDENCES'//)
            DO i = 1 , MEMax
               IF ( KVAr(i).NE.0 ) THEN
                  lc = 0
                  IF ( ELMl(i).GE..5 ) THEN
                     REWIND 17
 552                 READ (17,*) indx , iexp , sgm , ni , nf , al
                     IF ( indx.EQ.0 ) THEN
                        np = 0
                        nm = 0
                        DO j = 1 , lc
                           u = CORf(j,1)*CORf(j,2)*2.
                           IF ( ABS(u)/ELMu(i).GE..05 ) THEN
                              IF ( u.LT.0. ) nm = nm + 1
                              IF ( u.GT.0. ) np = np + 1
                           ENDIF
                        ENDDO
                        WRITE (22,99007) i , np , nm
99007                   FORMAT (1X/5X,10('*'),5X,'M.E.',1X,1I3,5X,1I3,
     &                          1X,'POSITIVE COMPONENTS',20X,1I3,1X,
     &                          'NEGATIVE COMPONENTS'///30X,'POSITIVE',
     &                          52X,'NEGATIVE'//5X,'EXP',2X,
     &                          'TRANSITION',2X,'SIGMA',3X,'DERIVATIVE',
     &                          3X,'D(SIGMA**2)/D(ME)',4X,'I',1X,'EXP',
     &                          2X,'TRANSITION',2X,'SIGMA',3X,
     &                          'DERIVATIVE',3X,'D(SIGMA**2)/D(ME)')
                        DO l = 1 , K
                           ump = 0.
                           umm = 0.
                           DO j = 1 , lc
                              u = 2.*CORf(j,1)*CORf(j,2)
                              IF ( u.LT.0. ) THEN
                                 IF ( u.LE.umm ) THEN
                                    umm = u
                                    jm = j
                                 ENDIF
                              ELSEIF ( u.GE.ump ) THEN
                                 ump = u
                                 jp = j
                              ENDIF
                           ENDDO
                           WRITE (22,99008) IY(jp,1) , IY(jp,2) , 
     &                            IY(jp,3) , CORf(jp,1) , CORf(jp,2) , 
     &                            ump , IY(jm,1) , IY(jm,2) , IY(jm,3) , 
     &                            CORf(jm,1) , CORf(jm,2) , umm
99008                      FORMAT (6X,1I2,3X,1I2,'--',1I2,5X,1F4.1,4X,
     &                             1E9.2,7X,1E9.2,9X,'I',2X,1I2,3X,1I2,
     &                             '--',1I2,5X,1F4.1,4X,1E9.2,7X,1E9.2)
                           CORf(jp,1) = 0.
                           CORf(jm,1) = 0.
                        ENDDO
                     ELSE
                        IF ( indx.EQ.i ) THEN
                           lc = lc + 1
                           IY(lc,1) = iexp
                           IY(lc,2) = ni
                           IY(lc,3) = nf
                           CORf(lc,1) = sgm
                           CORf(lc,2) = al
                        ENDIF
                        GOTO 552
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
            RETURN
         ELSE
            ll = ll + 1
            ch = ch + sgm*sgm
            ux = 2.*sgm*u
            ELM(indx) = ELM(indx) + ux
            ELMu(indx) = ELMu(indx) + ABS(ux)
            GOTO 550
         ENDIF
      ELSE
         loc = (iexp-1)*LP2 + indx
         sgm = (e-a*EP(iexp))/b
         u = dy*EP(iexp)*b + a*ZETa(loc)/b
         WRITE (17,*) indx , iexp , sgm , ni , nf , u
         GOTO 500
      ENDIF
      END
