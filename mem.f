 
C----------------------------------------------------------------------
 
      INTEGER*4 FUNCTION MEM(N1,N2,N3)
      IMPLICIT NONE
      INTEGER*4 k , LAMda , LAMmax , LDNum , LEAd , msum , MULti , N1 , 
     &          n1m , N2 , N3 , n3m
      COMMON /CLCOM / LAMda(8) , LEAd(2,500) , LDNum(8,75) , LAMmax , 
     &                MULti(8)
      msum = 0
      IF ( N3.NE.1 ) THEN
         n3m = N3 - 1
         DO k = 1 , n3m
            msum = msum + MULti(k)
         ENDDO
      ENDIF
      n1m = N1 - 1
      IF ( n1m.NE.0 ) THEN
         DO k = 1 , n1m
            msum = msum + LDNum(N3,k)
         ENDDO
      ENDIF
      n1m = msum + 1
      n3m = n1m + LDNum(N3,N1)
      DO k = n1m , n3m
         msum = msum + 1
         IF ( LEAd(2,k).EQ.N2 ) GOTO 100
      ENDDO
 100  MEM = msum
      END
