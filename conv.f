 
C----------------------------------------------------------------------
C FUNCTION CONV
C
C Called by: BRANR, PTICC, SEQ
C Calls:     LAGRAN
C
C Purpose: calculate the conversion coefficient at a particular energy by
C interpolating over the values provided by the user.
C
C Formal parameters:
C      Ega    - gamma energy
C      N      - multipolarity N=1,2,3 = E1,2,3 and N=4,5 = M1,2 (not as elsewhere!)
C
C Return value:
C      conversion coefficient interpolated to energy Ega

      REAL*8 FUNCTION CONV(Ega,N)
      IMPLICIT NONE

      INTEGER*4 isfirst, i, j, N, nenergies
      DATA isfirst/1/
      REAL*8 energies(1500), bricc(1500, 5), Ega
      SAVE energies, bricc, isfirst, nenergies

C     The first time, we need to read the data
      IF ( isfirst.eq.1 ) THEN
        isfirst = 0
        DO nenergies = 1, 1500
          READ(29,*,END=100) energies(nenergies),
     &      (bricc(nenergies,j),j=1,5)
        ENDDO
      ENDIF

C     Check multipolarity is valid
 100  IF ( N.LT.1.OR.N.GT.5 ) THEN
         CONV = 0.0
         RETURN
      ENDIF

C     Search for the energy in the list

      DO i = 1, nenergies
        IF (ABS(Ega - energies(i)) .LT. 1E-3) THEN
           CONV = bricc(i,N)
           return
        ENDIF
      ENDDO

C     We get here if the energy isn't in the list, so stop with an error
C     message
      WRITE (*,'(A,F7.3,A)')
     & 'Unable to find conversion coefficients for ',
     &  Ega, ' MeV'
      STOP 'Missing conversion coefficients'

      END
