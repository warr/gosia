 
C----------------------------------------------------------------------
C FUNCTION FAZA
C
C Called by: LAISUM
C
C Purpose: calculate phase.
 
      COMPLEX*16 FUNCTION FAZA(La,Mi,Rmu,Rsg)
      IMPLICIT NONE
      INTEGER*4 ieven , La , Mi
      REAL*8 Rmu , Rsg
      COMPLEX*16 ci
      DATA ci/(0.,1.)/ ! sqrt(-1)

      IF ( La.GT.6 ) THEN ! M1, M2 multipolarity
         FAZA = -ci
         IF ( Rmu.LT.0. ) FAZA = -FAZA
         IF ( La.EQ.7 ) RETURN
         IF ( Mi.EQ.2 ) RETURN
         FAZA = CMPLX(Rsg,0.)
         IF ( Rmu.LT.0. ) FAZA = -FAZA
         GOTO 99999
      ELSE                ! E1...6
         ieven = (-1)**Mi
         IF ( ieven.LE.0 ) THEN
            FAZA = -ci
            RETURN
         ENDIF
      ENDIF
      FAZA = CMPLX(Rsg,0.)
99999 END
