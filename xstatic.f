 
C----------------------------------------------------------------------
 
      SUBROUTINE XSTATIC(Iz,Ido,Iup,Beta)
      IMPLICIT NONE
      REAL*8 AKS , Beta , DQ , h , QCEn , VACdp , XNOr
      INTEGER*4 IBYp , Ido , Iup , Iz , lq
      COMMON /VAC   / VACdp(3,75) , QCEn , DQ , XNOr , AKS(6,75) , IBYp
      h = 1./(1.+(Iz**.45*.012008/Beta)**1.666667)
      QCEn = Iz*h**.6
      DQ = SQRT(QCEn*(1.-h))/2.
      Iup = INT(QCEn+3.*DQ+.5)
      Ido = INT(QCEn-3.*DQ-.5)
      IF ( Iup.GT.Iz ) Iup = Iz
      IF ( Ido.LT.1 ) Ido = 1
      XNOr = 0.
      DO lq = Ido , Iup
         XNOr = XNOr + EXP(-((QCEn-DBLE(lq))/DQ)**2/2.)
      ENDDO
      END
