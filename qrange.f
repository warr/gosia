 
C----------------------------------------------------------------------
 
      SUBROUTINE QRANGE(Icnt,Nlm,Lloc,Ibm,Icm,Idm,Irl)
      IMPLICIT NONE
      INTEGER*4 Ibm , Icm , Icnt , Idm , IRA , Irl , is , k , ke , km , 
     &          l , LAMda , LAMmax , ld , LDNum , LEAd , Lloc , ls , 
     &          MAXla , MULti
      INTEGER*4 nlend , Nlm
      DIMENSION Lloc(8) , Irl(8)
      COMMON /RNG   / IRA(8) , MAXla
      COMMON /CLCOM / LAMda(8) , LEAd(2,500) , LDNum(8,75) , LAMmax , 
     &                MULti(8)
      IF ( Icnt.EQ.1 ) THEN
         Nlm = 0
         DO l = 1 , 8
            Lloc(l) = 0
            Irl(l) = 0
         ENDDO
         DO k = 1 , 6
            ke = 7 - k
            km = 13 - k
            IF ( km.LE.8 ) THEN
               IF ( MULti(km).NE.0 ) THEN
                  Nlm = Nlm + 1
                  Lloc(Nlm) = km
                  Irl(Nlm) = IRA(km)
               ENDIF
            ENDIF
            IF ( MULti(ke).NE.0 ) THEN
               Nlm = Nlm + 1
               Lloc(Nlm) = ke
               Irl(Nlm) = IRA(ke)
            ENDIF
         ENDDO
         nlend = INT((DBLE(Nlm)+1.1)/2.)
         DO k = 1 , nlend
            ke = Nlm - k + 1
            ls = Lloc(ke)
            is = Irl(ke)
            Lloc(ke) = Lloc(k)
            Irl(ke) = Irl(k)
            Lloc(k) = ls
            Irl(k) = is
         ENDDO
         l = 0
         DO k = 1 , 6
            IF ( MULti(k).NE.0 ) l = k
         ENDDO
         Icm = MIN(4,l)
         Ibm = 2*l
         Idm = l
         l = 0
         DO k = 7 , 8
            ke = k - 6
            IF ( MULti(k).NE.0 ) l = ke
         ENDDO
         Ibm = MAX(Ibm,2*l)
         Idm = MAX(Idm,l)
         IF ( Icm.EQ.1 .AND. l.GT.1 ) Icm = 2
         MAXla = Lloc(1)
         RETURN
      ELSE
         IF ( Irl(Nlm).GE.Icnt ) RETURN
         ld = Lloc(Nlm)
         Lloc(Nlm) = 0
         Nlm = Nlm - 1
         IF ( Nlm.EQ.0 ) RETURN
         IF ( ld.GT.6 ) RETURN
         l = Lloc(Nlm)
         IF ( l.GT.6 ) l = l - 6
         Icm = MIN(2,l)
         Ibm = 2*l
         Idm = l
      ENDIF
      END
