 
C----------------------------------------------------------------------
 
      INTEGER*4 FUNCTION LEADF(N1,N2,N3)
      IMPLICIT NONE
      INTEGER*4 k , LAMda , LAMmax , LDNum , LEAd , lsum , MULti , N1 , 
     &          n1m , N2 , N3 , n3m
      COMMON /CLCOM / LAMda(8) , LEAd(2,500) , LDNum(8,75) , LAMmax , 
     &                MULti(8)
      lsum = 0
      n3m = N3 - 1
      IF ( n3m.NE.0 ) THEN
         DO k = 1 , n3m
            lsum = lsum + MULti(k)
         ENDDO
      ENDIF
      n1m = N1 - 1
      IF ( n1m.NE.0 ) THEN
         DO k = 1 , n1m
            lsum = lsum + LDNum(N3,k)
         ENDDO
      ENDIF
      n1m = lsum + N2
      LEADF = LEAd(2,n1m)
      END
