
C----------------------------------------------------------------------
C SUBROUTINE INTG
C
C Called by: FTBM, GOSIA
C Calls:     AMPDER, DOUBLE, HALF, RESET
C
C Purpose: the main integration routine.
C
C Uses global variables:
C      ACC50  - accuracy required for integration
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      D2W    - step in omega (= 0.03)
C      IFAC   - spin/parity phase factor
C      IFLG   - flag to determine whether to calculate exponential (so we don't calculate twice)
C      INTERV - default accuracy check parameter (see OP,CONT:INT)
C      IPATH  - index of substate in level with same m as substate Irld
C      IRA    - limit of omega for integration for each multipolarity
C      ISG    - sign of omega
C      ISMAX  - number of substates used
C      ISO    - Isotropic flag
C      KDIV   - index for division
C      LAMR   - flag = 1 if we should calculate this multipolarity
C      MAXLA  - multipolarity to calculate
C      NDIV   - number of divisions
C      NMAX   - number of levels
C      NPT    - index in ADB array (this is omega / 0.03)
C      NSTART - index in CAT of first substate associated with a level
C      NSW    - step in omega
C
C Formal parameters:
C      Ien    - experiment number
C
C Note that if it finds that the step size for the integral is too small, it
C calls DOUBLE to increase it by a factor of two, or if it finds that the
C step size is too big, it decreases it by a factor of two by calling HALF.
C
C We use the the 4th order Adams-Moulton predictor-corrector method for
C solving an ordinary differential equation. We use an adaptive version, which
C can change the step size (increase or decrease) in order to get the desired
C accuracy.
C
C The predictor is given as:
C
C y(n+1)_p = y(n) + h/24 * {55*f(n) - 59*f(n-1) + 37*f(n-2) - 9*f(n-3)}
C
C and the corrector is:
C
C y(n+1)_c = y(n) + h/24 * {9*f_p(n+1) + 19*f(n) - 5*f(n-1) + f(n-2)}
C
C The error is |E(n+1)| ~ 19/270 * {y_p(n+1) - y_c(n+1)}
C
C In this function:
C                   D2W        = h
C                   ARM(ir, 1) = f(n-3)
C                   ARM(ir, 2) = f(n-2)
C                   ARM(ir, 3) = f(n-1)
C                   ARM(ir, 4) = f(n)
C                   ARM(ir, 5) = y(n) initially
C                   ARM(ir, 5) = y_c(n+1) finally
C                   ARM(ir, 6) is not used
C                   ARM(ir, 7) = y_p(n+1)
C
C The function RESET is called to advance n by one. i.e. f(n-3) is set to the
C old value of f(n-2), f(n-2) to the old value of f(n-1) and f(n-1) to the old
C value of f(n).


      SUBROUTINE INTG(Ien)
      IMPLICIT NONE
      REAL*8 f , rim , rl , srt
      INTEGER*4 i , i57 , Ien , ihold , intend , ir , ir1 , k , kast ,
     &          mir , n
      COMPLEX*16 hold
      INCLUDE 'coex.inc'
      INCLUDE 'az.inc'
      INCLUDE 'rng.inc'
      INCLUDE 'a50.inc'
      INCLUDE 'clcom0.inc'
      INCLUDE 'clcom8.inc'
      INCLUDE 'coex2.inc'
      INCLUDE 'caux.inc'
      INCLUDE 'fla.inc'
      INCLUDE 'cexc0.inc'
      INCLUDE 'pth.inc'
      INCLUDE 'cexc9.inc'

      real*8 adamtemp,TCABS                  ! Rachel modification
c      logical diderrcheck                    ! Rachel modification
      integer*4 MEM                          ! Rachel modification
      integer*4 indx                         ! Rachel modification
      INCLUDE 'adbxi.inc'                    ! Rachel modification
      INCLUDE 'prt.inc'                      ! Rachel modification

c-------ADDITIONAL OUTPUT TO BE USED BY RACHEL.PY. MAR. 16 2011------------
      ! Rachel modification: print the collision functions to unit 99
      IF (IPRM(9).LT.0) CALL SPITQ(Ien,ABS(IPRM(9)))
      IF (IPRM(9).eq.11) THEN                ! If option was to print exc. amp. of substates
c     Write a blank line to mark the start of the experiment.
        WRITE(99,14617)
14617   FORMAT("   ")
      END IF
c-------END OF ADDITIONAL RACHEL OUTPUT.-----------------------------------

      intend = INTERV(Ien) ! Default accuracy set by INT option of OP,CONT
      D2W = .03D0 ! We use steps of 0.03 in omega
      NSW = 1
      kast = 0
      NDIV = 0
      KDIV = 0
 100  IF ( (NPT+NSW).GT.IRA(MAXLA) .AND. ISG.GT.0 ) RETURN
      DO i = 1 , 8
         LAMR(i) = 0
         IF ( (NPT+NSW).LT.IRA(i) ) LAMR(i) = 1
      ENDDO
C     Predictor
      IF ( ISO.EQ.0 ) THEN
         DO n = 1 , NMAX ! For each level
            ir = NSTART(n) - 1 ! First substate - 1
 120        ir = ir + 1
            ARM(ir,7) = ARM(ir,5)
     &                  + D2W/24.D0*(55.D0*ARM(ir,4)-59.D0*ARM(ir,3)
     &                  +37.D0*ARM(ir,2)-9.D0*ARM(ir,1))
            mir = INT(CAT(ir,3)) ! m quantum number of substate ir
            ir1 = ir - 2*mir
            ARM(ir1,7) = IFAC(n)*ARM(ir,7)
            IF ( DBLE(mir).LT.-0.1D0 ) GOTO 120
         ENDDO
      ELSE
         DO ir = 1 , ISMAX
            ARM(ir,7) = ARM(ir,5)
     &                  + D2W/24.D0*(55.D0*ARM(ir,4)-59.D0*ARM(ir,3)
     &                  +37.D0*ARM(ir,2)-9.D0*ARM(ir,1))
         ENDDO
      ENDIF
      NPT = NPT + NSW*ISG ! NPT loops over omega values, ISG is -1 at first then +1
      IF ( NPT.GT.0 ) THEN
         IF ( NDIV.EQ.0 ) GOTO 200
         KDIV = KDIV + 1
         IF ( KDIV.LT.NDIV ) GOTO 200
         KDIV = 0
         NPT = NPT + ISG
         IF ( NPT.GT.0 ) GOTO 200
      ENDIF
      NPT = -NPT + 2 ! We decreased omega to zero, so now start increasing
      ISG = 1
 200  CALL RESET(ISO)
      IFLG = 1
      i57 = 7 ! Tell LAISUM to use ARM(I,7) for excitation amplitudes

C     Calculate derivatives of amplitudes
      CALL AMPDER(i57)

C     Corrector
      IF ( ISO.EQ.0 ) THEN
         DO n = 1 , NMAX ! For each level
            ir = NSTART(n) - 1 ! First substate - 1
 220        ir = ir + 1
            ARM(ir,5) = ARM(ir,5)
     &                  + D2W/24.D0*(9.D0*ARM(ir,4)+19.D0*ARM(ir,3)
     &                  -5.D0*ARM(ir,2)+ARM(ir,1))
            mir = INT(CAT(ir,3)) ! m quantum number of substate ir
            ir1 = ir - 2*mir
            ARM(ir1,5) = IFAC(n)*ARM(ir,5)
            IF ( DBLE(mir).LT.-0.1D0 ) GOTO 220
         ENDDO
      ELSE
         DO ir = 1 , ISMAX ! For each substate
            ARM(ir,5) = ARM(ir,5)
     &                  + D2W/24.D0*(9.D0*ARM(ir,4)+19.D0*ARM(ir,3)
     &                  -5.D0*ARM(ir,2)+ARM(ir,1))
         ENDDO
      ENDIF
      kast = kast + 1
      IFLG = 0
      i57 = 5 ! Tell LAISUM to use ARM(I,5) for excitation amplitudes

C     Calculate derivatives of amplitudes
      CALL AMPDER(i57)
      IF ( (LAMR(2)+LAMR(3)).NE.0 ) THEN
         IF ( kast.GE.intend ) THEN
            kast = 0
            f = 0.
            DO k = 1 , NMAX ! For each level
               ihold = IPATH(k)
               IF ( ihold.NE.0 ) THEN
                  hold = ARM(ihold,5) - ARM(ihold,7)
                  rl = DBLE(hold)
                  rim = DIMAG(hold)
                  srt = rl*rl + rim*rim
                  f = MAX(f,srt)
               ENDIF
            ENDDO

C           Decide if we have appropriate accuracy (strictly it should be
C           f = SQRT(f)*19./270. but the difference is not all that large).
C
            f = SQRT(f)/14.D0
            IF ( f.GT.ACCUR .OR. f.LT.ACC50 ) THEN
               IF ( f.LT.ACC50 ) THEN
                  CALL DOUBLE(ISO) ! Double step size
                  D2W = 2.D0*D2W
                  NSW = 2*NSW
                  intend = INT((DBLE(intend)+.01D0)/2.D0)
                  IF ( intend.EQ.0 ) intend = 1
                  IF ( NSW.LT.1 ) THEN
                     NDIV = INT((DBLE(NDIV)+.01D0)/2.D0)
                     IF ( NDIV.LT.2 ) THEN
                        NDIV = 0
                        NSW = 1
                     ENDIF
                  ENDIF
               ELSE
                  CALL HALF(ISO) ! Halve step size
                  D2W = D2W/2.D0
                  NSW = INT((DBLE(NSW)+.01D0)/2.D0)
                  intend = 2*intend
                  IF ( NSW.LT.1 ) THEN
                     NDIV = 2*NDIV
                     IF ( NDIV.EQ.0 ) NDIV = 2
                  ENDIF
               ENDIF
            ENDIF

         ENDIF ! if kast>=intend
      ENDIF

c-------ADDITIONAL OUTPUT TO BE USED BY RACHEL.PY. MAR. 16 2011------------
c     I am trying to output the probabilities and amplitudes at each step

      if(IPRM(9).eq.11) then                ! If option was to print exc. amp. of substates
        write(99,14619) NPT,D2W
14619   format(2X,I4,2x,f5.3,$)
c       if(diderrcheck) then
c         write(99,14621) sqrt(f)/14.              ! print out the error term
c14621     format(2x,E11.4,$)
c       else
c         write(99,14622)                    !  if no error (f) term calc'd
c14622     format(' none',$)                  ! error wasn't checked
c       end if
      end if
      if((IPRM(9).GT.0).and.(IPRM(9).le.6)) then ! if option was to print adiab exp
c       Note that it looks up the the terms by multipolarity
c       so I select it by lambda = IPRM(9)
        indx = MEM(1,2,IPRM(9))                  ! Index for matrix element from level N to level m with multipolarity La
        write(99,14699)DBLE(EXPO(indx)),DIMAG(EXPO(indx))
14699   format(2x,'2',2x,D11.4,2x,D11.4,$)
        write(99,14623)    ! close the line
      else if(IPRM(9).eq.11) then                 ! if option was to print excitation amplitudes of substates
        do ir = 1, ismax
          adamtemp = TCABS(ARM(ir,i57))**2    ! probability(step)
          write(99,14618)ir,DBLE(ARM(ir,i57)),DIMAG(ARM(ir,i57)),
     &                   adamtemp
        enddo
14618   format(2x,I4,2x,D11.4,2x,D11.4,2x,D11.4,2x,$)
        write(99,14623)    ! close the line
      endif


14623 format('')
c-------END OF ADDITIONAL RACHEL OUTPUT.-----------------------------------


      GOTO 100
      END
