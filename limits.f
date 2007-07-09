 
C----------------------------------------------------------------------
 
      SUBROUTINE LIMITS
      IMPLICIT NONE
      REAL*8 ELM , ELMl , ELMu , SA
      INTEGER*4 IVAr , j , LMAxe , MAGexc , MEMax , MEMx6
      COMMON /CEXC  / MAGexc , MEMax , LMAxe , MEMx6 , IVAr(500)
      COMMON /COMME / ELM(500) , ELMu(500) , ELMl(500) , SA(500)
      DO j = 1 , MEMax
         IF ( IVAr(j).NE.0 ) THEN
            IF ( ELM(j).GT.ELMu(j) .OR. ELM(j).LT.ELMl(j) ) THEN
               IF ( ELM(j).GT.ELMu(j) ) THEN
                  ELM(j) = ELMu(j)
                  WRITE (22,99001) j , ELM(j)
               ELSE
                  ELM(j) = ELMl(j)
                  WRITE (22,99001) j , ELM(j)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
99001 FORMAT (2X,'Warning - matrix element ',1I3,' reset to ',1F10.6)
      END
