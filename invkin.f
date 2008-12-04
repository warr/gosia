 
C----------------------------------------------------------------------
C SUBROUTINE CMLAB
C
C Called by: GOSIA
C
C Purpose: calculate the angle of the scattered projectile in the lab frame
C          when the user gave the angle of the recoiling target nucleus in
C          the lab frame. There are two solutions to this problem, so Iflag
C          = 1 selects the larger angle one and Iflag = 2 the smaller one.
C          Note that the smaller angle (Iflag = 2) corresponds to very low
C          energies of the recoiling target nucleus, which probably either
C          don't get out of the target or don't get detected. So Iflag = 2
C          is probably not very useful!
C
C Formal parameters:
C      E_p     - Beam energy in MeV (readonly)
C      E_x     - energy of excited state to use for kinematic in MeV (readonly)
C      M_p     - mass of projectiel nuclei in AMU (readonly)
C      M_t     - mass of target nuclei in AMU (readonly)
C      Theta_t - theta of recoiling target nucleus in lab frame (readonly)
C      Theta_p - theta of scattered projectile in lab frame (writeonly)
C      Iflag   - flag to select one of two possible solutions (readonly)
C      Ikin    - kinematic flag (writeonly)
C
C  Note: Athough this code calculates the appropriate kinematic flag for use
C        by CMLAB, it should be noted that this flag isn't quite what it is
C        documented to be. Specifically, it is NOT the point corresponding to
C        90° in the centre of mass frame, but the point where the projectile
C        is maximally scattered. So for the moment, meshpoints should not be
C        given between these two points.

      SUBROUTINE INVKIN(E_p, E_x , M_p, M_t , Theta_t , Theta_p ,
     &                  Iflag , Ikin)

      REAL*8 E_p , M_p , M_t , Theta_t , Theta_p , E_x
      REAL*8 ared , epmin , t , x(2), y
      INTEGER*4 Iflag , Ikin

C     Reduced mass

      ared = 1 + M_p / M_t

C     Excitation energy of inelastically scattered particle when state at
C     energy E_x is excited

      epmin = E_p - E_x * ared

C     Tau

      taup = sqrt(E_p / epmin)
      tau = taup * M_p / M_t

C     Calculate the two solutions

      y = tan(theta_t/57.2957795)
      t = taup * taup * y * y * y * y -
     &      (1 + y * y) * (taup * taup * y * y - 1)
      t = sqrt(t)
      x(1) = (taup * y * y + t) / (1 + y * y)
      x(1) = atan2(sqrt(1 - x(1) * x(1)), x(1) + tau)
      x(2) = (taup * y * y - t) / (1 + y * y)
      x(2) = atan2(sqrt(1 - x(2) * x(2)), x(2) + tau)

C     Figure out which kinematic we need for this angle
      
      IF ( y*taup.GT.1. ) THEN
         Ikin = 1
      ELSE
         Ikin = 0
      ENDIF

C     Select the solution we want according to the flag. Note that the
C     solution with the lower angle corresponds to target recoils which
C     are probably undetectable.

      IF ( Iflag.EQ.1 ) THEN
         theta_p = MAX(x(1),x(2))*57.2957795
      ELSE
         theta_p = MIN(x(1),x(2))*57.2957795
      ENDIF

      END
