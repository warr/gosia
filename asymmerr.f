C----------------------------------------------------------------------
C SUBROUTINE ASYMMERR
C
C Called by: BRANR, CHMEM, DECAY, MIXR
C
C Purpose: calculate the contribution to chisqr and chilo due taking
C          asymmetric errorbars into account
C
C Formal parameters:
C      X      - the value obtained by the fit
C      X0     - the best value from literature
C      Lo     - the lower error bar (i.e. limit is X0 - Lo)
C      Hi     - the higher error bar (i.e. limit is X0 + Hi)
C      Chisq  - chi squared
C      Chilo  - chi squared of logs
C
      SUBROUTINE ASYMERR(X, X0, Lo, Hi, Chisq, Chilo)
      IMPLICIT NONE
      REAL*8 X, X0, Lo, Hi, Chisq, Chilo
      REAL*8 u, sigma, A

C     Calculate average sigma and asymmetry
      sigma = (Lo + Hi) * 0.5D0
      A = (hi - lo) / (hi + lo)
      
C     Use Taylor expansion to fourth term (must be an even power). If the
C     asymmetry of the errorbars is zero, this reduces to the usual
C     definition of chisq
      u = (X - X0) / sigma
      Chisq = Chisq + u * u * (1.D0 - 2.D0 * A * u + 5 * A * A * u * u)

      u = X0 * LOG(ABS(X/X0)) / sigma
      Chilo = Chilo + u * u * (1.D0 - 2.D0 * A * u + 5 * A * A * u * u)
      END
