 
C----------------------------------------------------------------------
 
      SUBROUTINE CHMEM(Nw,Chi,Chilo)
      IMPLICIT NONE
      REAL*8 Chi , Chilo , di , EAMx , ELM , ELMl , ELMu , SA
      INTEGER*4 ia , IAMx , IAMy , ib , NAMx , Nw
      COMMON /ME2D  / EAMx(100,2) , NAMx , IAMx(100) , IAMy(100,2)
      COMMON /COMME / ELM(500) , ELMu(500) , ELMl(500) , SA(500)
      IF ( NAMx.EQ.0 ) RETURN
      Nw = Nw + NAMx
      DO ia = 1 , NAMx
         ib = IAMx(ia)
         IF ( IAMy(ia,1).NE.IAMy(ia,2) ) THEN
            di = (ELM(ib)-EAMx(ia,1))/EAMx(ia,2)
            Chilo = Chilo + 
     &              (LOG(ABS(ELM(ib)/EAMx(ia,1)))*ABS(EAMx(ia,1))
     &              /EAMx(ia,2))**2
            Chi = Chi + di*di
         ELSE
            di = (ELM(ib)-EAMx(ia,1))/EAMx(ia,2)
            Chilo = Chilo + 
     &              (LOG(ABS(ELM(ib)/EAMx(ia,1)))*ABS(EAMx(ia,1))
     &              /EAMx(ia,2))**2
            Chi = Chi + di*di
         ENDIF
      ENDDO
      END
