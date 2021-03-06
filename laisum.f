
C----------------------------------------------------------------------
C SUBROUTINE LAISUM
C
C Called by: AMPDER, STING
C Calls:     FAZA
C
C Purpose: evaluate the sum  over matrix elements.
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      ELM    - matrix elements
C      EXPO   - adiabatic exponential
C      ISG    - sign of omega
C      ISG1   - index of omega
C      ISHA   - is half-integer spin
C      ISO    - isotropic flag
C      ISSTAR - first substate for given matrix element index
C      ISSTO  - last substate for given matrix element index
C      KDIV   - index for division
C      LOCQ   - location of collision functions in ZETA array
C      LP7    - start of collision functions in ZETA (45100)
C      MSTORE - index of final level number and index of matrix element
C      NDIV   - number of divisions
C      NPT    - index in ADB array (this is omega / 0.03)
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C      ZETA   - various coefficients
C
C Formal parameters:
C      Ir     - index of substate
C      N      - index of level
C      Rsg    - sign of omega
C      Lam    - multipolarity
C      Ld     - number of levels connected to level N by this multipolarity Lam
C      Nz     - index into ZETA array for this multipolarity
C      I57    - switch which is either 5 or 7 so we access ARM(I,5) or ARM(I,7)
C
C   \sum_{lmn} \zeta^{lm}_{kn} . M^(1)_{kn} f_{lm}(\omega) a_n(\omega)
C where
C   f_{lm} = -i Q_{lm} exp(i \xi_{kn} (\epsilon \sinh(\omega) + \omega))
C and
C   M^(1)_{kn} = <k||E(M)\lambda||n>
C
C EXPO is exp(i * xi * sinh(w) + w) calculated in function EXPON.
C ARM are the excitation amplitudes of the substates.
C q is the Qe or Qm calculated by the functions QE and QM, respectively and
C stored in ZETA array in the function SNAKE.
C z is the coupling parameter zeta, calculated in the function LSLOOP.

      SUBROUTINE LAISUM(Ir,N,Rsg,Lam,Ld,Nz,I57)
      IMPLICIT NONE
      REAL*8 q , rmir , rmis , rmu , Rsg , z
      INTEGER*4 i2 , i3 , I57 , iii , indq , indx , Ir , irs , is ,
     &          is1 , is2 , ismin , isplus
      INTEGER*4 la , Lam , Ld , m , mrange , mua , N , Nz
      COMPLEX*16 FAZA , pamp , pamp1
      INCLUDE 'pspin.inc'
      INCLUDE 'az.inc'
      INCLUDE 'coex.inc'
      INCLUDE 'caux.inc'
      INCLUDE 'pint.inc'
      INCLUDE 'adbxi.inc'
      INCLUDE 'ccoup.inc'
      INCLUDE 'mgn.inc'
      INCLUDE 'clcom8.inc'
      INCLUDE 'comme.inc'
      INCLUDE 'allc.inc'
      INCLUDE 'cexc0.inc'

      rmir = CAT(Ir,3) ! m quantum number of substate Ir
      iii = 0
      IF ( Lam.GT.6 ) iii = 1
      la = Lam
      IF ( Lam.GT.6 ) Lam = Lam - 6
      DO i2 = 1 , Ld
         pamp = (0.,0.)
         m = MSTORE(1,i2) ! Index of final level
         indx = MSTORE(2,i2) ! Index of matrix element in ELM
         ismin = 0
         is1 = NSTART(m) ! Index of first substate for level m
         IF ( is1.NE.0 ) THEN
            isplus = INT(rmir-CAT(is1,3)) - Lam
            IF ( isplus.LT.0 ) THEN
               ismin = isplus
               isplus = 0
            ENDIF
            is2 = is1 + isplus - 1
            mrange = 2*Lam + 1 + ismin
            IF ( is2+mrange.GT.NSTOP(m) ) mrange = NSTOP(m) - is2
            IF ( mrange.GT.0 ) THEN
               DO i3 = 1 , mrange
                  is = is2 + i3
                  rmis = CAT(is,3) ! m quantum number of substate is
                  IF ( ISO.NE.0 .OR. rmir.LE..1 .OR. rmis.LE..1 ) THEN
                     rmu = rmis - rmir
                     mua = INT(ABS(rmu) + 1.1) ! delta-m + 1
C                    Only consider electromagnetic and delta-m .NE. 0 magnetic
C                    contribution
                     IF ( la.LE.6 .OR. mua.NE.1 ) THEN
                        indq = LOCQ(Lam,mua) + NPT ! Index to Q function
                        Nz = Nz + 1                ! Index to Zeta
                        z = ZETA(Nz)               ! Zeta
                        q = ZETA(indq+LP7)         ! Q-function
                        IF ( NDIV.NE.0 ) q = ZETA(indq+LP7) + DBLE(KDIV)
     &                       *(ZETA(indq+LP7+ISG1)-ZETA(indq+LP7))
     &                       /DBLE(NDIV)
                        pamp1 = FAZA(la,mua,rmu,Rsg)*q*z
                        IF ( ISO.NE.0 .OR. rmir.LE..1 ) THEN
                           pamp = pamp1*ARM(is,I57) + pamp
                           IF ( ISO.EQ.0 .AND. rmis.GT..1 ) GOTO 10
                        ENDIF
                        IF ( N.NE.m ) THEN ! Not same level
                           irs = (-1)**(INT(rmir+rmis)-ISHA+iii) ! ISHA = 1 if half-integer spins, iii=0 for E, 1 for M
                           ARM(is,6) = ARM(is,6) + irs*pamp1*ARM(Ir,I57)
                           ISSTAR(i2) = MIN(is,ISSTAR(i2))
                           ISSTO(i2) = MAX(is,ISSTO(i2))
                        ENDIF
                     ENDIF
                  ENDIF ! If ISO.NE. 0 or either substate is spin 1
 10               CONTINUE
               ENDDO ! Loop on mrange
               IF ( N.EQ.m ) THEN ! N and m are level numbers, so if it is the same level, EXPO = 1
                  ARM(Ir,4) = ARM(Ir,4) + pamp*ELM(indx)
               ELSE
                  ARM(Ir,4) = ARM(Ir,4) + pamp*ELM(indx)*EXPO(indx)
               ENDIF
            ENDIF
         ENDIF
      ENDDO ! Loop on levels
      Lam = la
      END
