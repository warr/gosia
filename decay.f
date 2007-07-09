 
C----------------------------------------------------------------------
 
      SUBROUTINE DECAY(Chisq,Nlift,Chilo)
      IMPLICIT NONE
      REAL*8 AKS , bsum , Chilo , Chisq , DELla , DELta , df , DQ , 
     &       el1 , ELM , ELMl , ELMu , emt , emt1 , ENDec , ENZ , EPS , 
     &       EROot , FIEx , FP
      REAL*8 gk , GKP , QCEn , SA , TAU , TIMel , VACdp , vcd , XNOr , 
     &       ZETa
      INTEGER*4 i , IAXs , ibra , IBYp , idr , idrh , IEXp , ifn , il , 
     &          inx , inx1 , ITMa , iu , j , jlt , k , kl , KLEc , kq , 
     &          KSEq
      INTEGER*4 l , l1 , lc1 , lc2 , LIFct , LZEta , n1 , n2 , NDIm , 
     &          Nlift , NMAx , NMAx1
      COMMON /TRA   / DELta(500,3) , ENDec(500) , ITMa(50,200) , 
     &                ENZ(200)
      COMMON /LIFE1 / LIFct(50) , TIMel(2,50)
      COMMON /VAC   / VACdp(3,75) , QCEn , DQ , XNOr , AKS(6,75) , IBYp
      COMMON /CCOUP / ZETa(50000) , LZEta(8)
      COMMON /LEV   / TAU(75) , KSEq(500,4)
      COMMON /COEX2 / NMAx , NDIm , NMAx1
      COMMON /COMME / ELM(500) , ELMu(500) , ELMl(500) , SA(500)
      COMMON /KIN   / EPS(50) , EROot(50) , FIEx(50,2) , IEXp , IAXs(50)
      COMMON /CATLF / FP(4,500,3) , GKP(4,500,2) , KLEc(75)
      COMMON /LCDL  / DELla(500,3)
      DIMENSION gk(4)
      idr = 1
      DO il = 1 , NMAx1
         l = KSEq(idr,3)
         n1 = 28*(l-1)
         ibra = KLEc(l)
         bsum = 0.
         idrh = idr
         DO j = 1 , ibra
            inx = KSEq(idr,1)
            inx1 = KSEq(idr,2)
            el1 = 0.
            IF ( inx.NE.0 ) el1 = ELM(inx)
            emt = el1*el1
            DELla(idr,1) = emt
            IF ( inx1.NE.0 ) emt1 = ELM(inx1)*ELM(inx1)
            bsum = bsum + DELta(idr,1)*emt
            IF ( inx1.NE.0 ) THEN
               DELla(idr,3) = el1*ELM(inx1)
               DELla(idr,2) = emt1
               bsum = bsum + DELta(idr,2)*emt1
            ENDIF
            idr = idr + 1
         ENDDO
         idr = idrh
         TAU(l) = 1./bsum
         CALL GKVAC(l)
         DO j = 1 , ibra
            l1 = KSEq(idr,4)
            n2 = 28*(l1-1)
            inx1 = KSEq(idr,2)
            DO i = 1 , 4
               gk(i) = GKP(i,idr,1)*DELla(idr,1)
            ENDDO
            IF ( inx1.NE.0 ) THEN
               DO i = 1 , 4
                  gk(i) = gk(i) + GKP(i,idr,2)*DELla(idr,2)
               ENDDO
            ENDIF
            DO i = 1 , 4
               vcd = 1.
               IF ( i.NE.1 ) vcd = VACdp(i-1,l)
               gk(i) = gk(i)*TAU(l)
               ifn = 2*i - 1
               iu = (i-1)*7
               IF ( IAXs(IEXp).EQ.0 ) ifn = 1
               DO kq = 1 , ifn
                  lc1 = n1 + iu + kq
                  lc2 = n2 + iu + kq
                  ZETa(lc2) = ZETa(lc2) + gk(i)*vcd*ZETa(lc1)
               ENDDO
            ENDDO
            idr = idr + 1
         ENDDO
      ENDDO
      IBYp = 1
      IF ( Nlift.NE.0 .AND. IEXp.EQ.1 ) THEN
         DO jlt = 1 , Nlift
            kl = LIFct(jlt)
            df = (TAU(kl)-TIMel(1,jlt))/TIMel(2,jlt)
            Chilo = Chilo + (LOG(TAU(kl)/TIMel(1,jlt))*TIMel(1,jlt)
     &              /TIMel(2,jlt))**2
            Chisq = Chisq + df*df
         ENDDO
      ENDIF
      DO l = 2 , NMAx
         IF ( KLEc(l).NE.0 ) THEN
            n1 = 28*(l-1)
            DO j = 1 , 4
               vcd = 1.
               IF ( j.NE.1 ) vcd = VACdp(j-1,l)
               ifn = 2*j - 1
               iu = (j-1)*7
               DO k = 1 , ifn
                  lc1 = n1 + iu + k
                  ZETa(lc1) = ZETa(lc1)*vcd
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      END
