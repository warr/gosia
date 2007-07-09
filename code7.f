 
C----------------------------------------------------------------------
 
      SUBROUTINE CODE7(Ir,Is,N,Mt,Inqa,Indx)
      IMPLICIT NONE
      INTEGER*4 IAPr , idm , idn , Indx , Inqa , IPAth , Ir , Is , 
     &          ISEx , ism , MAGa , Mt , N
      REAL*8 QAPr
      COMMON /PTH   / IPAth(75) , MAGa(75)
      COMMON /APRCAT/ QAPr(500,2,7) , IAPr(500,2) , ISEx(75)
      IAPr(Indx,1) = N
      IAPr(Indx,2) = Mt
      IF ( IPAth(N).EQ.0 .OR. IPAth(Mt).EQ.0 ) THEN
         Inqa = -1
         GOTO 99999
      ELSE
         idn = Ir - IPAth(N)
         idm = Is - IPAth(Mt)
         ism = idn + idm + 3
         IF ( ism.EQ.2 ) THEN
            Inqa = 2
            IF ( idn.GT.idm ) Inqa = 3
            RETURN
         ELSEIF ( ism.EQ.3 ) THEN
            Inqa = 4
            RETURN
         ELSEIF ( ism.EQ.4 ) THEN
            Inqa = 5
            IF ( idn.GT.idm ) Inqa = 6
            RETURN
         ELSEIF ( ism.NE.5 ) THEN
            Inqa = 1
            RETURN
         ENDIF
      ENDIF
      Inqa = 7
99999 END
