 
C----------------------------------------------------------------------
C SUBROUTINE RESET
C
C Called by: INTG
C
C Purpose: to advance by one step. This means f(n-3) is set to the old value
C          of f(n-2), f(n-2) is set to the old value of f(n-1) and f(n-1) is
C          set to the old value of f(n).
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      ISMAX  - number of substates used
C      NMAX   - number of levels
C      NSTART - index in CAT of first substate associated with a level
C
C Formal parameters:
C      Iso    - isotropic flag
 
      SUBROUTINE RESET(Iso)
      IMPLICIT NONE
      REAL*8 CAT
      INTEGER*4 ir , ISMAX , Iso , j , NSTART , NSTOP
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      INCLUDE 'az.inc'
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      INCLUDE 'coex2.inc'
      
      IF ( Iso.EQ.0 ) THEN
         DO j = 1 , NMAX ! Loop over levels
            ir = NSTART(j) - 1 ! Index of first substate of level - 1
 20         ir = ir + 1
            ARM(ir,1) = ARM(ir,2)
            ARM(ir,2) = ARM(ir,3)
            ARM(ir,3) = ARM(ir,4)
            IF ( CAT(ir,3).LT.-.1 ) GOTO 20 ! m quantum number of substate ir
         ENDDO
         GOTO 99999
      ENDIF
       
      DO j = 1 , ISMAX ! Loop over substates
         ARM(j,1) = ARM(j,2)
         ARM(j,2) = ARM(j,3)
         ARM(j,3) = ARM(j,4)
      ENDDO
99999 END
