 
C----------------------------------------------------------------------
C SUBROUTINE SNAKE
C
C Called by: FTBM, GOSIA
C Calls:     QE, QM, QRANGE
C
C Purpose: evaluate and store the dimensionless collision functions Qe and Qm.
C
C Uses global variables:
C      CH     - table of cosh values
C      COLLIS - collision functions
CDEBUGC      DOMEGA - Initial step in omega (default = 0.03)
C      EPS    - epsilon
C      EROOT  - sqrt(epsilon^2 -1)
C      LOCQ   - location of collision function in COLLIS array
C      SH     - table of sinh values
C
C Formal parameters:
C      Nexp   - experiment number
C      Zpol   - dipole term (GDR excitation)
C
C The function QE is used to calculate Qe and QM to calculate Qm, but first
C we call QRANGE to determine the range over which we need to calculate them.
C
C The results are stored in the COLLIS array.
C
C LOCQ (in ALLC) is used as an index to these values.
C
C EROOT is set in CMLAB to \sqrt(\epsilon^2 - 1).
C
C Note that when we call QE and QM that lmda = 1...6 for E1...6 and 7,8 for
C M1, M2.
 
      SUBROUTINE SNAKE(Nexp,Zpol)
      IMPLICIT NONE
      REAL*8 b10 , b12 , b2 , b4 , b6 , b8 , c , c2 , c4 , c6 , 
     &       chi , cq , d , d2 , d3 , d4 , d5 , d6
      REAL*8 ert , pol , shi , Zpol
      INTEGER*4 ibm , icm , icnt , idm , irl , j , k , lloc , lmd , 
     &          lmda , mimx , Nexp , nind , nlm
      DIMENSION lloc(8) , cq(7) , irl(8)
      INCLUDE 'kin.inc'
      INCLUDE 'ccoup.inc'
      INCLUDE 'mgn.inc'
      INCLUDE 'allc.inc'
      INCLUDE 'hiper.inc'
CDEBUG      INCLUDE 'wvary.inc'
      
      icnt = 0
 100  icnt = icnt + 1

C     Calculate range over which we will want Qe and Qm
      CALL QRANGE(icnt,nlm,lloc,ibm,icm,idm,irl)
      IF ( nlm.EQ.0 ) RETURN

C     Calculate some parameters, which we will pass to QE or QM
      chi = CH(icnt) ! \cosh(\omega)
      shi = SH(icnt) ! \sinh(\omega)
      b2 = EPS(Nexp)*chi + 1.
      pol = 1. - Zpol/b2 ! E1 polarisation term
      b2 = b2*b2 ! b^2 = (\epsilon \cosh(\omega) + 1)^2
      IF ( ibm.NE.2 ) THEN
         b4 = b2*b2
         IF ( ibm.NE.4 ) THEN
            b6 = b4*b2
            IF ( ibm.NE.6 ) THEN
               b8 = b4*b4
               IF ( ibm.NE.8 ) THEN
                  b10 = b6*b4
                  IF ( ibm.NE.10 ) b12 = b6*b6
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF ( icm.NE.0 ) THEN
         c = chi + EPS(Nexp) ! c = \cosh(\omega) + \epsilon
         IF ( icm.NE.1 ) THEN
            c2 = c*c
            IF ( icm.NE.2 ) THEN
               c4 = c2*c2
               IF ( icm.NE.4 ) c6 = c2*c4
            ENDIF
         ENDIF
      ENDIF
      IF ( idm.NE.0 ) THEN
         d = EROOT(Nexp)*shi ! d = \sinh(\omega) * \sqrt(epsilon^2 - 1)
         IF ( idm.NE.1 ) THEN
            d2 = d*d
            IF ( idm.NE.2 ) THEN
               d3 = d*d2
               IF ( idm.NE.3 ) THEN
                  d4 = d2*d2
                  IF ( idm.NE.4 ) THEN
                     d5 = d3*d2
                     IF ( idm.NE.5 ) d6 = d3*d3
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      DO j = 1 , nlm
         lmda = lloc(j)
         IF ( lmda.GT.6 ) THEN
            lmd = lmda
            lmda = lmda - 6
            ert = EROOT(Nexp)
            CALL QM(c,d,b2,b4,ert,lmda,cq)
            mimx = lmda
            DO k = 1 , mimx
               nind = LOCQ(lmd,k) + icnt
               COLLIS(nind) = cq(k) ! These are the collision functions
            ENDDO
         ELSE
            CALL QE(c,d,b2,c2,d2,b4,b6,d3,b8,c4,d4,b10,d5,b12,d6,lmda,
     &              pol,cq)
            mimx = lmda + 1
            DO k = 1 , mimx
               nind = LOCQ(lmda,k) + icnt
               COLLIS(nind) = cq(k) ! These are the collision functions
            ENDDO
         ENDIF
      ENDDO
      IF ( icnt.EQ.365) THEN
        STOP 'Sorry, I can only do 365 steps in omega. You need more!'
      ENDIF
      GOTO 100
      END
CDEBUG
CDEBUGC----------------------------------------------------------------------
CDEBUG
CDEBUG      subroutine spitq(ixpt,mult)
CDEBUG      integer*4 ixpt,mult
CDEBUG      REAL*8 COLLIS , ZETA
CDEBUG      INTEGER*4 LZETA
CDEBUG      COMMON /CCOUP / ZETA(155600) , LZETA(8) , COLLIS(4900)
CDEBUG      REAL*8 EPS, EROOT, FIEX
CDEBUG      INTEGER*4 IEXP, IAXS
CDEBUG      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
CDEBUG      INTEGER*4 LOCQ
CDEBUG      COMMON /ALLC  / LOCQ(8,7)
CDEBUG      INTEGER*4 IRA , MAXLA
CDEBUG      COMMON /RNG   / IRA(8) , MAXLA
CDEBUG      integer*4 mimmex
CDEBUG      INTEGER*4 ibm , icm , icnt , idm , irl  , k , lloc , 
CDEBUG     &           nind , nlm
CDEBUG      DIMENSION lloc(8) , irl(8)
CDEBUG      INTEGER*4 LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP9 , 
CDEBUG     &          LP10 , LP11 , LP12 , LP14
CDEBUG      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP9 , 
CDEBUG     &                LP10 , LP11 , LP12 , LP14
CDEBUG      real*8 collfunc,w0
CDEBUG
CDEBUGc     Passed in the experiment number (ixpt) and the multipolarity (mult).
CDEBUGc     I am doing only for electric for now.
CDEBUGc     This routine is to read the electric collision functions from the ZETA array and print
CDEBUGc     them to the standard output file (23)  
CDEBUG
CDEBUGC     Borrowed this line from subroutine below.
CDEBUG      w0 = IRA(MAXLA)*DOMEGA + DOMEGA ! Maximum omega to calculate for (steps of DOMEGA)
CDEBUG
CDEBUGc     Look at how LAIAMP uses the ZETA array.
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
CDEBUG         collfunc = ZETA(nind+LP7)  ! These are the collision functions
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
