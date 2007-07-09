 
C----------------------------------------------------------------------
 
      SUBROUTINE EFFIX(Ipd,En,Effi)
      IMPLICIT NONE
      REAL*8 ABC , AKAvka , d , Effi , En , enl , pw , s , t , THIck , 
     &       w , xx , yy
      INTEGER*4 i , Ipd , j , l , ll , n
      DIMENSION xx(51) , yy(51)
      COMMON /EFCAL / ABC(8,10) , AKAvka(8,200) , THIck(200,7)
      Effi = 1.E-6
      En = En + 1.E-24
      enl = LOG(En)
      DO i = 1 , 10
         ll = 11 - i
         j = ll
         IF ( enl.GE.ABC(8,ll) ) GOTO 100
         j = -1
      ENDDO
 100  IF ( j.EQ.-1 ) Effi = 1.E-10
      IF ( j.EQ.-1 ) RETURN
      IF ( j.EQ.1 .OR. j.EQ.10 ) THEN
         s = 0.
         DO l = 1 , 7
            IF ( ABS(THIck(Ipd,l)).GE.1.E-9 ) THEN
               t = EXP(ABC(l,j))
               d = THIck(Ipd,l)
               s = s + t*d
            ENDIF
         ENDDO
      ELSE
         IF ( j.EQ.9 ) THEN
            xx(1) = ABC(8,8)
            xx(2) = ABC(8,9)
            xx(3) = ABC(8,10)
         ELSE
            xx(1) = ABC(8,j)
            xx(2) = ABC(8,j+1)
            xx(3) = ABC(8,j+2)
         ENDIF
         s = 0.
         DO l = 1 , 7
            IF ( ABS(THIck(Ipd,l)).GE.1.E-9 ) THEN
               IF ( j.EQ.9 ) THEN
                  yy(1) = ABC(l,8)
                  yy(2) = ABC(l,9)
                  yy(3) = ABC(l,10)
               ELSE
                  yy(1) = ABC(l,j)
                  yy(2) = ABC(l,j+1)
                  yy(3) = ABC(l,j+2)
               ENDIF
               CALL LAGRAN(xx,yy,3,0,enl,t,1,1)
               s = s + EXP(t)*THIck(Ipd,l)
            ENDIF
         ENDDO
      ENDIF
      Effi = EXP(-s)
c FITEFF or GREMLIN check
      IF ( AKAvka(5,Ipd).GT.0. .AND. AKAvka(5,Ipd).LT.10. ) THEN
c FITEFF eff. calib. by P.Olbratowski use
c PJN@2000
         w = LOG(En/AKAvka(5,Ipd))
         pw = AKAvka(2,Ipd)*w
         IF ( En.LT.AKAvka(5,Ipd) ) pw = pw + 
     &        w*w*(AKAvka(3,Ipd)+w*AKAvka(4,Ipd))
         Effi = Effi*EXP(pw)*AKAvka(1,Ipd)
         RETURN
      ELSEIF ( AKAvka(5,Ipd).GE.10. ) THEN
c     JAERI calibration - TC, Nov.2000
         w = LOG(En/.511)
         Effi = EXP(AKAvka(1,Ipd)+AKAvka(2,Ipd)
     &          *w-EXP(AKAvka(3,Ipd)+AKAvka(4,Ipd)*w))
         GOTO 99999
      ELSE
c GREMLIN
         w = LOG(20.*En)
         pw = AKAvka(1,Ipd) + AKAvka(2,Ipd)*w + AKAvka(3,Ipd)
     &        *w*w + AKAvka(4,Ipd)*w*w*w
         Effi = Effi*EXP(pw)
         IF ( ABS(AKAvka(5,Ipd)).GE.1.E-9 ) THEN
            n = INT(AKAvka(6,Ipd)+.1)
            pw = w**n
            w = AKAvka(5,Ipd)/pw
            Effi = Effi*EXP(w)
         ENDIF
      ENDIF
      IF ( ABS(AKAvka(8,Ipd)).LT.1.E-9 ) RETURN
      w = (AKAvka(7,Ipd)-1000.*En)/AKAvka(8,Ipd)
      pw = EXP(w)
      IF ( ABS(pw-1.).LT.1.E-6 ) WRITE (22,99001)
99001 FORMAT (5x,'***** CRASH - EFFIX *****')
      Effi = Effi/(1.-pw)
99999 END
