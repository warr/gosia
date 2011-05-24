
C----------------------------------------------------------------------

      subroutine spitq(ixpt,mult)
      integer*4 ixpt,mult
      INCLUDE 'ccoup.inc'
      INCLUDE 'kin.inc'
      INCLUDE 'allc.inc'
      INCLUDE 'rng.inc'
      integer*4 mimmex
      INTEGER*4 ibm , icm , icnt , idm , irl  , k , lloc , 
     &           nind , nlm
      DIMENSION lloc(8) , irl(8)
      INCLUDE 'mgn.inc'
      real*8 collfunc,w0

c     Passed in the experiment number (ixpt) and the multipolarity (mult).
c     I am doing only for electric for now.
c     This routine is to read the electric collision functions from the ZETA array and print
c     them to the standard output file (23)  

C     Borrowed this line from subroutine below.
      w0 = IRA(MAXLA)*.03 + .03 ! Maximum omega to calculate for (steps of 0.03)

c     Look at how LAIAMP uses the ZETA array.
c     Are these the collision functions?

      write(99,72072)ixpt,mult
72072 format('The *stored* collision functions ',
     &       'for experiment ',i2,' and mult E',i1)
      write(99,71972)
71972 format('   |--------mu=0----------|--------mu=2---------',
     &       '|--------mu=2...')
      write(99,71672)
71672 format('   |icnt nind  Qe         |icnt nind  Qe        ',
     &       '|icnt nind  Qe...')
      icnt = 0
 100  icnt = icnt + 1


      mimmex = mult + 1
c     write(99,71772)mimmex

C     With QRANGE I think I am checking that icnt is still in range to index the QE values
c     for this multipolarity and experiment.  There must be only one epsilon 
c     value per experiment.
      CALL QRANGE(icnt,nlm,lloc,ibm,icm,idm,irl)
c     write(99,71772)mimmex
71772 format(i4)
c     test to here

      IF ( nlm.EQ.0 ) RETURN

      DO k = 1 , mimmex
         nind = LOCQ(mult,k) + icnt
         collfunc = ZETA(nind+LP7)  ! These are the collision functions
         write(99,71572)icnt,nind,collfunc
71572    format(1x,i4,1x,i6,1x,d11.4,$)
      ENDDO
      write(99,71872)
71872 format('')
      goto 100 

      end


C----------------------------------------------------------------------

