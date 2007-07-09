 
C----------------------------------------------------------------------
 
      SUBROUTINE MIXUP
      IMPLICIT NONE
      REAL*8 ELM , ELMl , ELMu , RNDM , SA , SE
      INTEGER*4 IVAr , k , k1 , LMAxe , MAGexc , MEMax , MEMx6
      COMMON /COMME / ELM(500) , ELMu(500) , ELMl(500) , SA(500)
      COMMON /CEXC  / MAGexc , MEMax , LMAxe , MEMx6 , IVAr(500)
      COMMON /XRA   / SE
      DO k = 1 , MEMax
         IF ( IVAr(k).NE.0 .AND. IVAr(k).LE.999 ) ELM(k) = ELMl(k)
     &        + RNDM(SE)*(ELMu(k)-ELMl(k))
      ENDDO
      DO k = 1 , MEMax
         IF ( IVAr(k).GE.999 ) THEN
            k1 = IVAr(k) - 1000
            IF ( ABS(ELMu(k1)).LT.1.E-9 ) THEN
               ELM(k) = 0.
            ELSE
               ELM(k) = ELM(k1)*SA(k)
            ENDIF
         ENDIF
      ENDDO
      END
