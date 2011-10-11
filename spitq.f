CDEBUG
CDEBUGC----------------------------------------------------------------------
CDEBUG
CDEBUG      SUBROUTINE SPITQ(Ixpt,Mult)
CDEBUG      INTEGER*4 ixpt,mult
CDEBUG      INCLUDE 'ccoup.inc'
CDEBUG      INCLUDE 'kin.inc'
CDEBUG      INCLUDE 'allc.inc'
CDEBUG      INCLUDE 'rng.inc'
CDEBUG      INTEGER*4 mimmex
CDEBUG      INTEGER*4 ibm , icm , icnt , idm , irl  , k , lloc , 
CDEBUG     &           nind , nlm
CDEBUG      DIMENSION lloc(8) , irl(8)
CDEBUG      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
CDEBUG     &          LP10 , LP11 , LP12
CDEBUG      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
CDEBUG     &                LP10 , LP11 , LP12
CDEBUG      REAL*8 collfunc,w0
CDEBUG
CDEBUGc     Passed in the experiment number (ixpt) and the multipolarity (mult).
CDEBUGc     I am doing only for electric for now.
CDEBUGc     This routine is to read the electric collision functions from the COLLIS array and print
CDEBUGc     them to the standard output file (23)  
CDEBUG
CDEBUGC     Borrowed this line from subroutine below.
CDEBUG      w0 = IRA(MAXLA)*DOMEGA + DOMEGA ! Maximum omega to calculate for (steps of DOMEGA)
CDEBUG
CDEBUGc     Look at how LAIAMP uses the COLLIS array.
CDEBUGc     Are these the collision functions?
CDEBUG
CDEBUG      write(22,72072)ixpt,mult
CDEBUG72072 format('The *stored* collision functions ',
CDEBUG     &       'for experiment ',i2,' and mult E',i1)
CDEBUG      write(22,71972)
CDEBUG71972 format('|--------mu=1----------|--------mu=2---------',
CDEBUG     &       '|--------mu=2...')
CDEBUG      write(22,71672)
CDEBUG71672 format('|icnt nind  Qe         |icnt nind  Qe        ',
CDEBUG     &       '|icnt nind  Qe...')
CDEBUG      icnt = 0
CDEBUG 100  icnt = icnt + 1
CDEBUG
CDEBUG
CDEBUG      mimmex = mult + 1
CDEBUGc     write(22,71772)mimmex
CDEBUG
CDEBUGC     With QRANGE I think I am checking that icnt is still in range to index the QE values
CDEBUGc     for this multipolarity and experiment.  There must be only one epsilon 
CDEBUGc     value per experiment.
CDEBUG      CALL QRANGE(icnt,nlm,lloc,ibm,icm,idm,irl)
CDEBUGc     write(22,71772)mimmex
CDEBUG71772 format(i4)
CDEBUGc     test to here
CDEBUG
CDEBUG      IF ( nlm.EQ.0 ) RETURN
CDEBUG
CDEBUG      DO k = 1 , mimmex
CDEBUG         nind = LOCQ(mult,k) + icnt
CDEBUG         collfunc = COLLIS(nind)  ! These are the collision functions
CDEBUG         write(22,71572)icnt,nind,collfunc
CDEBUG71572    format(1x,i3,1x,i5,1x,d11.4,$)
CDEBUG      ENDDO
CDEBUG      write(22,71872)
CDEBUG71872 format('')
CDEBUG      goto 100 
CDEBUG
CDEBUG      end
CDEBUG
CDEBUG
CDEBUG
