
C----------------------------------------------------------------------
C SUBROUTINE SPITQ
C
C Called by: INTG
C Calls:
C
C Purpose: write out collision functions to a file
C
C Uses global variables:
C
C Formal parameters:
C      Ixpt   - experiment number
C      Mult   - multipolarity

      SUBROUTINE SPITQ(Ixpt,Mult)
      IMPLICIT NONE
      INCLUDE 'ccoup.inc'
      INCLUDE 'kin.inc'
      INCLUDE 'allc.inc'
      INCLUDE 'rng.inc'
      INCLUDE 'mgn.inc'
      INTEGER*4 mimmex
      INTEGER*4 Ixpt , Mult
      INTEGER*4 ibm , icm , icnt , idm , irl  , k , lloc ,
     &           nind , nlm
      DIMENSION lloc(8) , irl(8)
      REAL*8 collfunc

c     Passed in the experiment number (Ixpt) and the multipolarity (Mult).
c     I am doing only for electric for now.
c     This routine is to read the electric collision functions from the ZETA array and print
c     them to the standard output file (23)

c     Look at how LAIAMP uses the ZETA array.
c     Are these the collision functions?

      WRITE(99,72072) Ixpt , Mult
72072 FORMAT('The *stored* collision functions ',
     &       'for experiment ' , i2 , ' and mult E' , i1)
      WRITE(99,71972)
71972 FORMAT('   |--------mu=0----------|--------mu=2---------',
     &       '|--------mu=2...')
      WRITE(99,71672)
71672 FORMAT('   |icnt nind  Qe         |icnt nind  Qe        ',
     &       '|icnt nind  Qe...')
      icnt = 0
 100  icnt = icnt + 1


      mimmex = Mult + 1
c     WRITE(99,71772) mimmex

C     With QRANGE I think I am checking that icnt is still in range to index the QE values
c     for this multipolarity and experiment.  There must be only one epsilon
c     value per experiment.
      CALL QRANGE(icnt , nlm , lloc , ibm , icm , idm , irl)
c     WRITE(99,71772) mimmex
c71772 FORMAT(i4)
c     test to here

      IF ( nlm.EQ.0 ) RETURN

      DO k = 1 , mimmex
         nind = LOCQ(Mult,k) + icnt
         collfunc = ZETA(nind+LP7)  ! These are the collision functions
         WRITE(99,71572) icnt , nind , collfunc
71572    FORMAT(1X,I4,1X,I6,1X,D11.4,$)
      ENDDO
      WRITE(99,71872)
71872 FORMAT('')
      GOTO 100

      END
