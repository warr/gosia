      INCLUDE 'header.txt'
C
C PROGRAM GOSIA
C
C Calls: ADHOC, ALLOC, ANGULA, ARCCOS, ARCTG, CMLAB, COORD, DECAY, DJMM,
C        EFFIX, ELMT, FAKP, FHIP, FTBM, INTG, INVKIN, KLOPOT, KONTUR, LAGRAN,
C        LOAD, MINI, MIXR, MIXUP, OPENF, PATH, PRELM, PTICC, QFIT, READY,
C        SETIN, SIMIN, SNAKE, SPLNER, STING, TACOS, TAPMA, TEMB, TENS, WSIXJ,
C        WTHREJ
C
C Uses global variables:
C      ABC    - absorption coefficients
C      ACCA   - accuracy
C      ACCUR  - accuracy required
C      AGELI  - angles of Ge detectors
C      AKAVKA - efficiency curve parameters
C      ARM    - excitation amplitudes of substates.
C      AVJI   - average J (N.B. here it is G(1))
C      B      - table of factorials
C      BEQ    - identifier for angle for rotations
C      BETAR  - recoil beta
C      CAT    - substates of levels (n_level, J, m)
C      CC     - conversion coefficients
C      CNOR   - normalization factors
C      CORF   - internal correction factors
C      DEVD   -
C      DEVU   -
C      DIPOL  - E1 polarization parameter
C      DIX    - Ge parameters (inner & outer radius, length, distance)
C      DLOCK  - limit derivative below which matrix element is fixed if LOCKS=1
C      DS     - integrated rutherford cross-section
C      DSE    - rutherford cross section at given energy integrated over angles
C      DSG    - differential gamma-ray yield at meshpoints
C      DSIGS  - dsigma for each experiment
C      DYEX   - error on experimental yield
C      EAMX   - known matrix elements and their errors
C      ELM    - matrix elements
C      ELMH   -
C      ELML   - lower limit on matrix elements
C      ELMU   - upper limit on matrix elements
C      EMMA   - Controls number of magnetic substates in full coulex calc.
C      EN     - energy of level
C      EP     - bombarding energy
C      ERR    - error flag
C      EXPO   - adiabatic exponential
C      FIEL   - K (N.B. here it is G(6))
C      FIEX   - phi range of particle detector
C      GAMMA  - Gamma (N.B. here it is G(2))
C      GFAC   - g (N.B. here it is G(5))
C      GRAD   - partial derivative of chi squared wrt. each matrix element
C      HLM    - matrix elements before minimisation
C      HLMLM  - old value of matrix element or chi squared
C      IAMX   - index of matrix element for known matrix element
C      IAMY   - level indices of pair of levels for which matrix element is known
C      IAX    - axial symmetry flag
C      IBYP   - flag to indicate whether we calculate <\alpha_k>
C      ICLUST - cluster number for each experiment and detector
C      ICS    - read internal correction factors flag (OP,CONT switch CRF,)
C      IDIVE  - number of subdivisions
C      IDRN   - index of normalising transition for yields
C      IEXP   - experiment number
C      IFAC   - spin/parity phase factor
C      IFBFL  - calculate derivatives with forward-backward method
C      IFMO   - include correction to angular distance for finite recoil distance.
C      ILE    - yield number for each detector
C      IMIN   -
C      INHB   - inhibit error flag (LERF) setting in POMNOZ
C      INNR   - independent normalisation switch (see OP,CONT INR,)
C      INTERV - default accuracy check parameter for Adams-Moulton (see OP,CONT:INT)
C      INTR   - flag to swap chisqr and log(chisqr)
C      IP     - table of prime numbers
C      IPRM   - various flags to control output
C      IPS1   - terminate after calculating and storing internal correction factors
C      IRAWEX - flag to indicate raw uncorrected yield
C      ISEX   -
C      ISKIN  - kinematic flag (0,1)
C      ISMAX  - number of substates used
C      ISO    - isotropic flag
C      ITMA   - identify detectors according to OP,GDET
C      ITS    - create tape 18 file (OP,CONT switch SEL,)
C      ITTE   - thick target experiment flag
C      IUNIT3 - unit for TAPE3
C      IVAR   - indicates a limit or correlation is set
C      IWF    - warning flag
C      IY     - index for yields
C      IZ     - Z of investigated nucleus
C      IZ1    - Z of non-investigated nucleus
C      JENTR  - flag set to 0 normally, 1 in OP,ERRO
C      JSKIP  - Experiments to skip during minimisation.
C      JZB    - unit to read from
C      KFERR  - error flag for minimization
C      KSEQ   - index of level
C      KVAR   -
C      LAMAX  - number of multipolarities to calculate
C      LAMBDA - list of multipolarities to calculate
C      LASTCL - index of last detector in cluster
C      LDNUM  - number of matrix elements with each multipolarity populating each level
C      LEAD   - pair of levels involved in each matrix element
C      LIFCT  - index for lifetimes
C      LMAX   - ground-state spin + 1
C      LMAXE  - maximum multipolarity needed for calculation
C      LNORM  - normalisation constant control
C      LNY    - use logs to calculate chi squared
C      LOCKF  - flag to fix matrix elements with most significant derivative
C      LOCKS  - lock flag. If LOCKS=1, fix at first stage of minimization
C      LP1    - maximum number of experiments (50)
C      LP10   - maximum number of substates (1200)
C      LP11   - LP8 - 1 (103)
C      LP12   - number of steps of omega (365)
C      LP13   - LP9 + 1 (47901)
C      LP14   - maximum space for collision functions (4900)
C      LP2    - maximum number of matrix elements (1500)
C      LP3    - maximum number of levels (100)
C      LP4    - maximum number of yields (1500)
C      LP6    - maximum number of gamma detectors (32)
C      LP7    - start of collision functions (45100)
C      LP8    - (104)
C      LP9    - length of ZETA - 2100 (47900)
C      MAGA   - number of magnetic substates in approximate calculation
C      MAGEXC - flag: 0 means no magnetic excitations, 1 means with mag. exc.
C      MEMAX  - number of matrix elements
C      MEMX6  - number of matrix elements with E1...6 multipolarity
C      MULTI  - number of matrix elements having given multipolarity
C      NAMX   - number of known matrix elements
C      NANG   - number of gamma-ray detectors for each experiment
C      NBRA   - number of branching ratios
C      NCM    - calculate kinematics assuming this state for final state (default = 2)
C      NDIM   - maximum number of levels
C      NDST   - number of data sets
C      NEXPT  - number of experiments
C      NLIFT  - number of lifetimes
C      NLOCK  - number of elemnts to fix if LOCKF=1
C      NMAX   - number of levels
C      NMAX1  - number of levels with decays
C      NYLDE  - number of yields
C      ODL    - results of OP,GDET calculation
C      PARX   - [for maps]
C      PARXM  - [for maps]
C      POWER  - x (N.B. here it is G(7))
C      QAPR   - approximate Coulomb amplitudes
C      SA     - ratio of matrix elements for correlated elements
C      SE     - seed for random number generator of OP,RAND
C      SGW    - number of standard deviations to generate warning (see control option WRN,X)
C      SPIN   - spin of level
C      SUBCH1 - partial chisqr
C      SUBCH2 - partial chisqr
C      SUMCL  - sum of yields for clusters
C      TAU    - lifetime in picoseconds
C      THICK  - thickness of each absorber type
C      TIMEC  - Tau_C (N.B. here it is G(4))
C      TIMEL  - lifetimes and their errors
C      TLBDG  - theta of particle detector in lab frame (in degrees)
C      TREP   - theta of recoiling nucleus (in radians)
C      UPL    - upper limits for all gamma detectors
C      VINF   - speed of projectile at infinty
C      XA     - A of investigated nucleus
C      XA1    - A of non-investigated nucleus
C      XI     - xi coupling coefficients
C      XIR    - [for maps]
C      XLAMB  - Lambda* (N.B. here it is G(3))
C      XV     - energy meshpoints (sometimes theta meshpoints) where we calculate exact Coulex
C      YEXP   - experimental yields
C      YGN    - gamma yield calculated without correction to angular distribution from finite recoil distance
C      YGP    - gamma yield calculated with correction to angular distribution from finite recoil distance
C      YNRM   - relative normalisation for gamma detectors
C      YV     - scattering angle meshpoints where we calculate exact Coulex
C      ZETA   - various coefficients
C      ZPOL   - dipole term
C      ZV     - energy meshpoints

      PROGRAM GOSIA
      IMPLICIT NONE
      REAL*8 acof , ap , ARCCOS , ARCTG , arg , ax , bcof , be2 , 
     &       be2a , be2b , be2c
      REAL*8 bk , bl , bm , bmx , bten , bu , ccc , 
     &       ccd , cf , chilo , chiok , chis0 , chisl , chisq , chiss , 
     &       cnst
      REAL*8 cocos , conu , d , decen , dedx , dsd , dsig , dst
      REAL*8 dsx , dsxm , effi , eh1 , elmi , ELMT , emhl1 , emn , emx , 
     &       enb
      REAL*8 eng , enh , esd , esp , ess , 
     &       fi0 , fi1 , fic , fiex1 , figl , fipo1 , fm , gth
      REAL*8 hen , het , p , pfi , 
     &       ph1 , ph2 , pi , po1 , po2 , polm , pop1 , pr , pv
      REAL*8 q1 , q2 , qc , qfac , qr , qui , r , r1 , r2 , r3 , r4 , 
     &       rem , remax , rl , rlr , rm , rx , ry
      REAL*8 rz , s , s11 , s12 , s21 , s22 , sbe , sf , sh , sh1 , 
     &       sh2 , SIMIN , slim
      REAL*8 summm , sz1 , sz2 , TACOS , tau1 , tau2 , test , 
     &       tetrc , tfac , thc , title , tmn , tmx , todfi
      REAL*8 tta , tth , tting , ttttt , txx , u , 
     &       val , waga , wph , wpi , WSIXJ , wth , wthh , 
     &       WTHREJ
      REAL*8 xep , xi1 , xi2 , xk1 , xk2 , xl1 , xlevb , 
     &       xlk , xm1 , xm2 , xm3 , xtest , xw , xx , xxi , 
     &       ycorr
      REAL*8 yy , yyd1 , yydd , yyy , zmir , zp , zz
      REAL*8 ttttx ! Only gosia1 and pawel
      INTEGER*4 i , i122 , iapx , ib , ibaf , icg , icll , ict , ictl , 
     &          id , ideff , idf
      INTEGER*4 idr , iecd , ient , ifbp , ifc , ifm , ifwd , 
     &          ig1 , ig2 , ih1 , ih2 , ihlm , ihuj , ii , ij
      INTEGER*4 ija0 , ijaja , ijan , ijk , ijx , ile1 , ilevls , 
     &          ilx , im , imode , in1 , in2 , inclus , ind , 
     &          ind1 , ind2 , indx
      INTEGER*4 inko , inm1 , inm2 , inn , inpo , intend , intvh , 
     &          inva , inx1 , iobl , iocc , iopri , iosr , ipd , iph
      INTEGER*4 ipine , ipinf , ipo1 , ipo2 , ipo3 , ipp , iprc , 
     &          ipri , irea , irep , irfix , irix , isip , iske , iskf
      INTEGER*4 isko , iskok , isoh , ispa , ispb , itno , 
     &          itp , iuy , iva , iva1 , ivarh , ivari , ivrh
      INTEGER*4 ixj , ixl , ixm , iyr , izcap , j , ja , 
     &          jan , jan1 , jb , jb1 , jb2 , jd , jde , jdy , je
      INTEGER*4 jex , jexp , jfi , jfre , jgd , jgl , jgl1 , jgr , jgs , 
     &          jj , jj1 , jjjj , jjlx , jjx , jk , jkloo , jktt , jl , 
     &          jmm , jmpin
      INTEGER*4 jp , jphd , jpin , jrls , js , jt , jtp , jyi , jyi1 , 
     &          jyi2 , jyv , jz , k , kb , kclust , kerf , kex
      INTEGER*4 kh , kh1 , kh2 , kk , kk1 , kk2 , kkk , kl , kloop , 
     &          kmat , kq , ktt , kuku , l , la , la1 , lam , lamd
      INTEGER*4 lamh , lb , lck1 , lck2 , levl , lex , lexp , 
     &          lfagg , lfini , lh1 , lh2 , liscl , lkj
      INTEGER*4 lkj1 , ll , lli , lll , lmax1 , lmaxh , locat , 
     &          loct , lp0 , lpin
      INTEGER*4 ltrn , ltrn1 , ltrn2 , lu , lx , lxd , magh , MEM
      INTEGER*4 memax1 , memh , memx4 , mend , mexl , 
     &          mfla , mlt , mm , mpin , ms , n , na , na1 , naa , 
     &          nallow
      INTEGER*4 naxfl , nb1 , nb2 , nbands , nch , ndima , ndum , 
     &          ne , nf , nfd , nfdd , 
     &          nfi , nflr , nft , nged
      INTEGER*4 ngpr , ni , nksi , nl , nmaxh , nmemx , nnl , 
     &          nogeli , npce , npce1 , npct , npct1 , 
     &          npt , nptl , nptx , ns1
      INTEGER*4 ns2 , ntap , ntt , numcl , nval , nz
      INTEGER*4 iskin_protect
      CHARACTER*4 oph , op1 , opcja , op2
      CHARACTER*1 prp
      DIMENSION ihlm(32) , esp(20) , dedx(20) , bten(1600) , ! bten dimension = 16 * maxlevels
     &          fiex1(100,100,2) , title(20) , pfi(101) , zmir(6,2,50) , 
     &          iecd(50) , wpi(100,2) , tau1(10) , eng(10) , 
     &          tau2(10,7) , xl1(7) , qui(8,10) , cf(8,2) , 
     &          ivarh(1500) , liscl(200) , dsxm(100,100,100) , 
     &          levl(50) , xlevb(50,2) , bm(8,20,20,3) , mlt(1500) , 
     &          ivari(1500) , jpin(50) , ideff(50) , iskin_protect(50)
      INCLUDE 'clust.inc'
      INCLUDE 'cccds.inc'
      INCLUDE 'inhi.inc'
      INCLUDE 'ident.inc'
      INCLUDE 'efcal.inc'
      INCLUDE 'tcm.inc'
      INCLUDE 'brec.inc'
      INCLUDE 'adbxi.inc'
      INCLUDE 'dimx.inc'
      INCLUDE 'tra.inc'
      INCLUDE 'cinit.inc'
      INCLUDE 'xra.inc'
      INCLUDE 'hhh.inc'
      INCLUDE 'vac.inc'
      INCLUDE 'me2d.inc'
      INCLUDE 'life1.inc'
      INCLUDE 'dftb.inc'
      INCLUDE 'erran.inc'
      INCLUDE 'mgn.inc'
      INCLUDE 'seck.inc'
      INCLUDE 'vlin.inc'
      INCLUDE 'dumm.inc'
      INCLUDE 'brnch.inc'
      INCLUDE 'yexpt.inc'
      INCLUDE 'yteor.inc'
      INCLUDE 'lev.inc'
      INCLUDE 'map.inc'
      INCLUDE 'ccc.inc'
      INCLUDE 'ggg.inc'
      INCLUDE 'az.inc'
      INCLUDE 'kin.inc'
      INCLUDE 'cxi.inc'
      INCLUDE 'clcom.inc'
      INCLUDE 'coex.inc'
      INCLUDE 'minni.inc'
      INCLUDE 'cx.inc'
      INCLUDE 'cexc.inc'
      INCLUDE 'prt.inc'
      INCLUDE 'ccoup.inc'
      INCLUDE 'cb.inc'
      INCLUDE 'clm.inc'
      INCLUDE 'clcom0.inc'
      INCLUDE 'clcom8.inc'
      INCLUDE 'clcom9.inc'
      INCLUDE 'comme.inc'
      INCLUDE 'coex2.inc'
      INCLUDE 'cexc9.inc'
      INCLUDE 'caux0.inc'
      INCLUDE 'pth.inc'
      INCLUDE 'aprcat.inc'
      INCLUDE 'warn.inc'
      INCLUDE 'thtar.inc'
      INCLUDE 'fit.inc'
      INCLUDE 'aprx.inc'
      INCLUDE 'skp.inc'
      INCLUDE 'trb.inc'
      INCLUDE 'sel.inc'
      INCLUDE 'ercal.inc'
      INCLUDE 'logy.inc'
      INCLUDE 'fakul.inc'
      INCLUDE 'life.inc'
      INCLUDE 'switch.inc'
      DATA (eng(k),k=1,10)/.05 , .06 , .08 , .1 , .15 , .2 , .3 , .5 , 
     &      1. , 1.5/
C     Absorption coefficients in units of 1/cm for Ge
      DATA (tau1(k),k=1,10)/17.656 , 10.726 , 5.076 , 2.931 , 1.3065 , 
     &      .8828 , .5959 , .4357 , .3041 , .2472/
C     Absorption coefficients in units of 1/cm for Al, C, Fe, Cu, Ag/Cd/Sn, Ta
C     and Pb at the energies 0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5, 1, 1.5
C     MeV
      DATA (tau2(k,1),k=1,10)/.9883 , .7473 , .5442 , .4592 , .3718 , 
     &      .3302 , .2814 , .2278 , .1657 , .1350/
      DATA (tau2(k,2),k=1,10)/1.014 , .7443 , .5195 , .4261 , .3362 , 
     &      .2967 , .2518 , .2038 , .1479 , .1204/
      DATA (tau2(k,3),k=1,10)/15.167 , 9.405 , 4.652 , 2.889 , 1.525 , 
     &      1.135 , .8643 , .6592 , .4703 , .3830/
      DATA (tau2(k,4),k=1,10)/23.184 , 14.182 , 6.777 , 4.059 , 1.970 , 
     &      1.384 , .9936 , .7473 , .5274 , .4297/
      DATA (tau2(k,5),k=1,10)/84.351 , 51.445 , 23.822 , 13.070 , 
     &      4.774 , 2.605 , 1.339 , .7925 , .5005 , .4032/
      DATA (tau2(k,6),k=1,10)/93.364 , 58.559 , 125.96 , 70.713 , 
     &      25.302 , 12.541 , 5.193 , 2.215 , 1.077 , .8176/
      DATA (tau2(k,7),k=1,10)/89.809 , 56.338 , 27.009 , 62.966 , 
     &      22.933 , 11.334 , 4.540 , 1.813 , .8020 , .5900/
      DATA q1/0./,q2/0./,iph/0/
      DATA cnst/0./,sh1/0./,irfix/0/,jfre/0/ ! Only gosia1 and pawel

C     Initialize prime numbers
      IP(1) = 2
      IP(2) = 3
      IP(3) = 5
      IP(4) = 7
      IP(5) = 11
      IP(6) = 13
      IP(7) = 17
      IP(8) = 19
      IP(9) = 23
      IP(10) = 29
      IP(11) = 31
      IP(12) = 37
      IP(13) = 41
      IP(14) = 43
      IP(15) = 47
      IP(16) = 53
      IP(17) = 59
      IP(18) = 61
      IP(19) = 67
      IP(20) = 71
      IP(21) = 73
      IP(22) = 79
      IP(23) = 83
      IP(24) = 89
      IP(25) = 97
      IP(26) = 101

C     Initialize pointers
      lp0 = 50000 ! Size of ZETA array
      LP1 = 50 ! Maximum number of experiments
      LP2 = 1500 ! Maximum number of matrix elements
      LP3 = 100 ! Maximum number of levels
      LP4 = 1500
      LP6 = 32 ! Maximum number of gamma detectors
      LP7 = lp0 - 4900 ! Start of collision coefficients in ZETA
      LP8 = LP3*28 + 1
      LP9 = lp0 - LP3*28
      LP10 = 1200 ! Maximum number of substates
      LP11 = LP8 - 1
      LP12 = 365 ! Maximum number of steps of omega (dimension of ADB, SH, CH)
      LP13 = LP9 + 1
      LP14 = 4900 ! Maximum number of collision coefficients

      JZB = 5

C     Initialize normalization to 1.
      DO i = 1 , LP3 ! LP3 = 100 (maximum number of levels)
         DO j = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
            CNOR(j,i) = 1.
         ENDDO
      ENDDO

      IUNIT3 = 3 ! Is 33 in gosia2
      IBYP = 0
      INHB = 0
      BEQ = -983872.
      ipinf = 0
      iyr = 0
      pi = 3.141592654
      INNR = 0
      itno = 0
      chisq = 0.
      chilo = 0.
      IWF = 1 ! Turn on warnings
      ifm = 0 ! Fast minimisation switch off by default
      IPS1 = 11
      ifwd = -1
      INTR = 0
      LNY = 0
      JENTR = 0 ! Flag to indicate we are not in OP,ERRO
      ICS = 0
      ISPL = 0 ! Flag to indicate we should use LAGRAN not SPLNER

      DO i = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         ideff(i) = 0
         jpin(i) = 0
         iecd(i) = 0
      ENDDO
      txx = 0.
      SGW = 3.
      SUBCH1 = 0.
      SUBCH2 = 0.
      ITS = 0 ! Create tape 18 flag
      iosr = 0
      LOCKS = 0
      DLOCK = 1.1
      kerf = 0
      IFBFL = 0
      NLOCK = 0
      LOCKF = 0
      DO i = 1 , LP4 ! LP4 = 1500
         DO j = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
            CORF(i,j) = 1.
         ENDDO
      ENDDO
      DO i = 1 , 20
         IPRM(i) = 1
      ENDDO
      DO i = 1 , 50
         DO j = 1 , 5
            CC(i,j) = 0.
         ENDDO
      ENDDO
      IPRM(4) = -2
      IPRM(5) = 11111
      IPRM(6) = 11111
      IPRM(7) = 0
      IPRM(16) = 0
      IPRM(17) = 0
      IPRM(18) = 0
      IPRM(19) = 0
      IPRM(20) = 0
      DO i = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         DO j = 1 , 5
            IF ( j.NE.5 ) THEN
               DO k = 1 , 10
                  DO kuku = 1 , 6
                     PARXM(i,j,k,kuku) = 0.
                  ENDDO
               ENDDO
            ENDIF
            DO k = 1 , 12
               PARX(i,k,j) = 0.
            ENDDO
         ENDDO
      ENDDO
      DO k = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         IDIVE(k,1) = 1
         IDIVE(k,2) = 1
         DO iuy = 1 , 6
            XIR(iuy,k) = 0.
         ENDDO
      ENDDO
      iobl = 0
      lfagg = 0
      izcap = 12800
      KFERR = 0
      NDIM = LP3 ! LP3 = 100 (maximum number of levels)
      ISO = 1
      B(1) = 1.
      DO i = 2 , 20
         B(i) = B(i-1)*(i-1)
      ENDDO
      LMAXE = 0
      CALL FAKP
      CALL FHIP
      NCM = 2 ! Default final state for kinematics calculation (OP,CONT NCM,)
      DO ijx = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         INTERV(ijx) = 1
      ENDDO
      la = 0
      ipo3 = 1
      indx = 0
      ACCUR = .00001
      icg = 1
      ient = 1
      jphd = 1 ! Print header flag
      DIPOL = 0.005
      MAGEXC = 0 ! Initially flag that we don't need magnetic excitations
      LAMMAX = 0
      DO lam = 1 , 8
         DO lexp = 1 , LP3 ! LP3 = 100 (maximum number of levels)
            LDNUM(lam,lexp) = 0
         ENDDO
         MULTI(lam) = 0
         LAMDA(lam) = 0
      ENDDO
      DO j = 1 , LP2 ! LP2 = 1500 (maximum number of matrix elements)
         EXPO(j) = (1.,0.)
         KVAR(j) = 1
         ELM(j) = 0.
      ENDDO
      DO j = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         JSKIP(j) = 1
         ISKIN(j) = 0
      ENDDO
      DO j = 1 , LP3 ! LP3 = 100 (maximum number of levels)
         ISEX(j) = 1111
      ENDDO
      ISEX(1) = 0
      ACCA = .00001
      oph = '    '
      nmemx = LP2 + 9 ! LP2 = 1500 (maximum number of matrix elements)
      IEXP = 1
      IMIN = 0
      i122 = 0
      DO j = 1 , LP2 ! LP2 = 1500 (maximum number of matrix elements)
         DO k = 1 , 2
            DO l = 1 , 7
               QAPR(j,k,l) = 0.
            ENDDO
         ENDDO
      ENDDO
      ERR = .FALSE.
      intend = 0 ! End of initialization

C.............................................................................
C     Start reading input file.
 100  READ (JZB,99001) op1 , op2
99001 FORMAT (1A3,1A4)
      
      IF ( op1.EQ.'OP, ' ) THEN
         IF ( op2.EQ.'GOSI' ) oph = op2
         IF ( op2.EQ.'GOSI' ) opcja = op2

C        Treat OP,FILE (attach files to fortran units)
         IF ( op2.EQ.'FILE' ) THEN
            CALL OPENF
            GOTO 100 ! End of OP,FILE - back to input loop
         ENDIF

C        Print header         
         IF ( jphd.EQ.1 ) WRITE (22,99002)
99002    FORMAT ('1'/1X,125('*')/1X,125('*')/1X,50('*'),25X,50('*')/1X,
     &           50('*'),10X,'GOSIA',10X,50('*')/1X,50('*'),25X,50('*')
     &           /1X,125('*')/1X,125('*')////)
         IF ( jphd.EQ.1 ) WRITE (22,99003)
99003    FORMAT (1X/20X,'ROCHESTER COULOMB EXCITATION DATA ANALYSIS ',
     &           'CODE BY T.CZOSNYKA,D.CLINE AND C.Y.WU'/50X,
     &           'LATEST REVISION- JUNE  2006'//////)
         jphd = 0 ! Set print header flag to zero, so we don't repeat header

C        Handle OP,GDET (germanium detectors)
         IF ( op2.EQ.'GDET' ) THEN
            nl = 7
            READ (JZB,*) nfdd ! number of physical detectors

            nfd = ABS(nfdd) ! Negative value means graded absorber
            IF ( nfdd.LE.0 ) THEN
               REWIND 8
               DO i = 1 , nl
                  WRITE (8,*) (tau2(l,i),l=1,10)
               ENDDO
               WRITE (8,*) (eng(l),l=1,10)
            ENDIF

C           Write file for gamma-ray energy dependence of Ge solid-angle
C           attenuation coefficients
            REWIND 9
            WRITE (9,*) nfd
            DO i = 1 , nfd ! For each detector
               READ (JZB,*) (DIX(k),k=1,4) ! radius of core, outer radius, length, distance
               READ (JZB,*) (xl1(k),k=1,nl) ! thicknesses of 7 kinds of absorber
               IF ( DIX(1).LE.0. ) DIX(1) = .01
               WRITE (9,*) DIX(4) ! length
               IF ( nfdd.LE.0 ) WRITE (8,*) (xl1(k),k=1,nl)
               ind = 1
               IF ( xl1(5).GT.0. ) ind = 3
               IF ( xl1(6).GT.0. ) ind = 4
               IF ( xl1(7).GT.0. ) ind = 5
               WRITE (9,*) eng(ind) ! First energy
               CALL QFIT(qui,tau1,tau2,eng,xl1,cf,nl,ind)
               WRITE (22,99004) i
99004          FORMAT (10X,'DETECTOR',1X,1I2)
               DO k = 1 , 8
                  WRITE (22,99005) k , cf(k,1) , cf(k,2)
99005             FORMAT (1X,//5X,'K=',1I1,2X,'C1=',1E14.6,2X,'C2=',
     &                    1E14.6/5X,'ENERGY(MEV)',5X,'FITTED QK',5X,
     &                    'CALC.QK',5X,'PC.DIFF.'/)
                  WRITE (9,*) cf(k,1) , cf(k,2) , qui(k,ind)
                  DO l = 1 , 10
                     arg = (eng(l)-eng(ind))**2
                     qc = (qui(k,ind)*cf(k,2)+cf(k,1)*arg)/(cf(k,2)+arg)
                     WRITE (22,99006) eng(l) , qc , qui(k,l) , 
     &                                100.*(qc-qui(k,l))/qui(k,l)
99006                FORMAT (8X,1F4.2,6X,1F9.4,5X,1F9.4,3X,1E10.2)
                  ENDDO
               ENDDO
            ENDDO
            GOTO 100 ! End of OP,GDET - back to input loop

C        Treat OP,RAND (randomise matrix elements)
         ELSEIF ( op2.EQ.'RAND' ) THEN
            READ (JZB,*) SE ! Seed for random number generator
            CALL MIXUP
            WRITE (22,99007)
99007       FORMAT (1X///5X,'MATRIX ELEMENTS RANDOMIZED...'///)
            CALL PRELM(2)
            GOTO 100 ! End of OP,RAND - back to input loop

C        Treat OP,TROU (troubleshooting)
         ELSEIF ( op2.EQ.'TROU' ) THEN
            ITS = 1 ! Create tape 18 flag
            READ (JZB,*) kmat , rlr
            GOTO 100 ! End of OP,TROU - back to input loop

C        Treat OP,REST (restart)
         ELSEIF ( op2.EQ.'REST' ) THEN
            irix = 12
            REWIND irix
            memax1 = MEMAX + 1
            DO lkj = 1 , MEMAX
               READ (irix,*) ELM(lkj)
            ENDDO
            DO lkj = 1 , memax1
               READ (JZB,*) lkj1 , xlk
               IF ( lkj1.EQ.0 ) GOTO 120
               ELM(lkj1) = xlk
            ENDDO
 120        WRITE (22,99008)
99008       FORMAT (1X///5X,'*****',2X,
     &              'RESTART-MATRIX ELEMENTS OVERWRITTEN',2X,'*****'///)
            DO kk = 1 , MEMAX
               la = mlt(kk)
               IF ( ivari(kk).GE.10000 ) THEN
                  kk1 = ivari(kk)/10000
                  kk2 = ivari(kk) - 10000*kk1
                  la1 = la
                  IF ( kk2.GE.100 ) THEN
                     la1 = kk2/100
                     kk2 = kk2 - 100*la1
                  ENDIF
                  inx1 = MEM(kk1,kk2,la1)
C      ELML(KK)=ELML(INX1)*ELM(KK)/ELM(INX1)
C      ELMU(KK)=ELMU(INX1)*ELM(KK)/ELM(INX1)
                  SA(kk) = ELM(kk)/ELM(inx1)
                  IVAR(kk) = 1000 + inx1
                  IF ( ELMU(kk).LE.ELML(kk) ) THEN
                     elmi = ELMU(kk)
                     ELMU(kk) = ELML(kk)
                     ELML(kk) = elmi
                  ENDIF
               ENDIF
            ENDDO
            CALL PRELM(2) ! Parameter is 4 in gosia2
            GOTO 100 ! End of OP,REST - back to input loop

C     Treat OP,SELE
         ELSEIF ( op2.EQ.'SELE' ) THEN
            CALL SELECT
            GOTO 2000 ! End of execution

C     Treat OP,BRIC
         ELSEIF ( op2.EQ.'BRIC' ) THEN
            CALL BRICC
            GOTO 100 ! End of OP,BRIC - back to input loop

C        Treat other options
         ELSE

C           Treat OP,RE,A (release A)
            IF ( op2.EQ.'RE,A' ) GOTO 900
           
C           Treat OP,RE,F (release F)
            IF ( op2.EQ.'RE,F' ) GOTO 900

C           Treat OP,ERRO (calculate errors)
            IF ( op2.EQ.'ERRO' ) THEN
               READ (JZB,*) idf , ms , mend , irep , ifc , remax
               rem = LOG(remax)
               LOCKS = 0
               LOCKF = 0
               JENTR = 1 ! Flag to indicate we are in OP,ERRO
               sh = 1.
               ifbp = 0
               inpo = 1
               inko = 1
               IF ( iosr.NE.0 .AND. idf.NE.0 ) THEN
                  inn = 0
                  ij = MULTI(1)
                  IF ( ij.NE.0 ) THEN
                     DO ij = 1 , NMAX
                        lxd = LDNUM(1,ij)
                        IF ( lxd.NE.0 ) THEN
                           DO ijk = 1 , lxd
                              inn = inn + 1
                           ENDDO
                        ENDIF
                     ENDDO
                     inpo = inn + 1
                  ENDIF
                  DO ij = 1 , NMAX
                     lxd = LDNUM(2,ij)
                     IF ( lxd.NE.0 ) THEN
                        DO ijk = 1 , lxd
                           inn = inn + 1
                        ENDDO
                     ENDIF
                  ENDDO
                  inko = inn
                  IF ( irep.NE.2 ) THEN
                     WRITE (IUNIT3,*) NMAX , MEMAX , inpo , inko
                     DO inn = 1 , NMAX
                        WRITE (IUNIT3,*) inn , SPIN(inn) , EN(inn)
                     ENDDO
                     DO inn = 1 , MEMAX
                        WRITE (IUNIT3,*) inn , LEAD(1,inn) , LEAD(2,inn)
                     ENDDO
                     DO inn = 1 , MEMAX
                        WRITE (IUNIT3,*) inn , ELM(inn)
                     ENDDO
                  ENDIF ! IF ( irep.NE.2 )
               ENDIF ! IF ( iosr.NE.0 .AND. idf.NE.0 )
               IF ( irep.NE.0 ) THEN
                  REWIND 15
                  READ (15,*) (DEVD(kh1),DEVU(kh1),kh1=1,MEMAX)
               ELSE
                  DO kh1 = 1 , MEMAX
                     DEVD(kh1) = ELML(kh1) - ELM(kh1)
                     DEVU(kh1) = ELMU(kh1) - ELM(kh1)
                  ENDDO
               ENDIF
               IF ( IMIN.EQ.0 ) CALL CMLAB(0,dsig,ttttt)
               IF ( ERR ) GOTO 2000 ! Error
               IF ( IMIN.NE.0 ) GOTO 400
               GOTO 1300 ! End of OP,ERRO

C           Treat OP,RE,C (release C)
            ELSEIF ( op2.EQ.'RE,C' ) THEN
               jfre = 1
               irfix = 0
               GOTO 1000 ! End of OP,RE,C

C           Treat OP,TITL (title)
            ELSEIF ( op2.EQ.'TITL' ) THEN
               READ (JZB,99009) (title(k),k=1,20)
99009          FORMAT (20A4)
               WRITE (22,99010) (title(k),k=1,20)
99010          FORMAT (10X,20A4/10X,100('-'))
               GOTO 100 ! End of OP,TITL - back to input loop

            ELSE

C              Treat OP,GOSI
               IF ( op2.EQ.'GOSI' ) GOTO 200

C              Treat OP,COUL
               IF ( op2.EQ.'COUL' ) GOTO 200

C              Treat OP,EXIT
               IF ( op2.EQ.'EXIT' ) THEN
                  GOTO 430 ! End of OP,EXIT

C              Treat OP,MINI
               ELSEIF ( op2.EQ.'MINI' ) THEN
                  READ (JZB,*) imode , nptl , chiok , conu , xtest , 
     &                 LOCKF , NLOCK , IFBFL , LOCKS , DLOCK
                  op2 = opcja
                  IMIN = IMIN + 1
                  IF ( IMIN.NE.1 ) GOTO 1400
                  GOTO 1200 ! End of OP,MINI

C              Treat OP,THEO
               ELSEIF ( op2.EQ.'THEO' ) THEN
                  irix = 12
                  REWIND (irix)
                  ibaf = 1
                  DO jb = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
                     DO lb = 1 , 2
                        xlevb(jb,lb) = 0
                     ENDDO
                  ENDDO
                  READ (JZB,*) nbands ! Number of bands
                  IF ( nbands.LE.0 ) ibaf = 0
                  nbands = ABS(nbands)
                  DO nl = 1 , 8
                     DO jb = 1 , nbands
                        DO jl = 1 , nbands
                           DO kl = 1 , 3
                              bm(nl,jb,jl,kl) = 0.
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
                  DO jb = 1 , nbands
                     READ (JZB,*) bk , ilevls ! K of band, number of levels in band
                     READ (JZB,*) (levl(ib),ib=1,ilevls) ! Level list for band
                     DO kb = 1 , ilevls
                        inva = levl(kb)
                        xlevb(inva,2) = bk
                        xlevb(inva,1) = DBLE(jb)
                     ENDDO
                  ENDDO
                  DO nl = 1 , 8
                     READ (JZB,*) nnl ! Multipolarity
 126                 IF ( nnl.LE.0 ) GOTO 130
                     READ (JZB,*) jb1 , jb2 ! band indices
                     IF ( jb1.NE.0 ) THEN
                        READ (JZB,*) (bm(nnl,jb1,jb2,j),j=1,3) ! intrinsic moments
                        DO j = 1 , 3
                           bm(nnl,jb2,jb1,j) = bm(nnl,jb1,jb2,j)
                        ENDDO
                        GOTO 126
                     ENDIF
                  ENDDO
 130              DO kb = 1 , MEMAX
                     IF ( ibaf.NE.0 ) THEN
                        ind1 = LEAD(1,kb)
                        ind2 = LEAD(2,kb)
                        xi1 = SPIN(ind1)
                        xi2 = SPIN(ind2)
                        lamd = mlt(kb)
                        nb1 = INT(xlevb(ind1,1)+.1)
                        nb2 = INT(xlevb(ind2,1)+.1)
                        xk1 = xlevb(ind1,2)
                        xk2 = xlevb(ind2,2)
                        xm1 = bm(lamd,nb1,nb2,1)
                        xm2 = bm(lamd,nb1,nb2,2)
                        xm3 = bm(lamd,nb1,nb2,3)
                        ELM(kb) = ELMT(xi1,xi2,lamd,nb1,nb2,xk1,xk2,xm1,
     &                            xm2,xm3)
                        IF ( ABS(ELM(kb)).LT.1E-6 ) ELM(kb) = 1.E-6
                        irix = 12
                        WRITE (irix,*) ELM(kb)
                     ENDIF
                  ENDDO
                  GOTO 100 ! End of OP,THEO - back to input loop

C              Treat OP,YIEL
               ELSEIF ( op2.EQ.'YIEL' ) THEN
                  CALL ADHOC(oph,idr,nfd,ntap,iyr)
                  GOTO 100 ! End of OP,YIEL - back to input loop

C              Treat OP,INTG
               ELSEIF ( op2.EQ.'INTG' ) THEN
                  REWIND 14
                  lfagg = 1
                  IF ( SPIN(1).LT..25 ) ISO = 0
                  DO lx = 1 , NEXPT ! For each experiment
                     lpin = 1
                     IF ( ipinf.NE.0 ) THEN
                        IF ( jpin(lx).NE.0 ) lpin = jpin(lx)
                     ENDIF
                     IEXP = lx
                     tth = TLBDG(lx)
                     enh = EP(lx)
                     DO mpin = 1 , lpin ! For each pin diode
                        IF ( iecd(lx).EQ.1 ) THEN ! Circular detector
                           READ (JZB,*) ne , ntt , emn , emx , wth , 
     &                          wph , wthh
                           mfla = 1
                           CALL COORD(wth,wph,wthh,ntt,0,pfi,wpi,tth,lx,
     &                                tmn,tmx)
                        ELSE
                           READ (JZB,*) ne , ntt , emn , emx , tmn , tmx
                           mfla = 0
                           IF ( ntt.LT.0 ) mfla = 1
                        ENDIF
                        ntt = ABS(ntt)
                        jan = NANG(lx)
                        jan1 = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) jan1 = jan
                        IF ( iecd(lx).EQ.1 ) THEN ! Circular detector
                           WRITE (14,*) ne , ntt , emn , emx , tmn , 
     &                                  tmx , jan1 , wth , wph , wthh
                        ELSE
                           WRITE (14,*) ne , ntt , emn , emx , tmn , 
     &                                  tmx , jan1 , tmx , tmx , tmx
                        ENDIF
                        READ (JZB,*) (XV(i),i=1,ne)
                        IF ( iecd(lx).NE.1 ) READ (JZB,*)
     &                       (YV(i),i=1,ntt)
                        IF ( tth.LT.0. ) ELMH(2*lx-1) = YV(1)
                        IF ( tth.LT.0. ) ELMH(2*lx) = YV(ntt)
                        DO kloop = 1 , ne ! For each energy meshpoint
                           enb = XV(kloop)
                           EP(lx) = enb
                           DO ktt = 1 , ntt
                              tta = SIGN(YV(ktt),tth)
                              IF ( IAXS(lx).NE.0 ) THEN ! If not axial symmetry
                                 IF ( iecd(lx).NE.1 ) THEN
                                    IF ( kloop.EQ.1 ) THEN
                                       READ (JZB,*) nfi ! Number of phi ranges
                                       READ (JZB,*) 
     &                                    (fiex1(ktt,jfi,1),fiex1(ktt,
     &                                    jfi,2),jfi=1,nfi)
                                       IF ( tth.LT.0. ) THEN
                                         DO jfi = 1 , nfi ! For each phi angle
                                         fiex1(ktt,jfi,1)
     &                                      = fiex1(ktt,jfi,1) + 180.
                                         fiex1(ktt,jfi,2)
     &                                      = fiex1(ktt,jfi,2) + 180.
                                         ENDDO
                                       ENDIF
                                    ENDIF
                                 ENDIF
                              ENDIF ! If not axial symmetry
                              TLBDG(lx) = tta
                              IF ( kloop.EQ.1 ) THEN
                                 IF ( iecd(lx).NE.0 ) THEN
                                    nfi = 1
                                    fiex1(ktt,1,1) = wpi(ktt,1) ! Lower phi limit
                                    fiex1(ktt,1,2) = wpi(ktt,2) ! Upper phi limit
                                 ENDIF
                              ENDIF
                              CALL CMLAB(lx,dsig,tetrc)
                              IF ( ERR ) GOTO 2000 ! Error
                              tting = TLBDG(lx)
                              IF ( ERR ) GOTO 1900 ! Troubleshoot
                              CALL LOAD(lx,1,1,0.D0,jj)
                              CALL ALLOC(ACCUR)
                              CALL SNAKE(lx,ZPOL)
                              CALL SETIN
                              DO j = 1 , LMAX ! For each spin up to ground-state spin + 1
                                 polm = DBLE(j-1) - SPIN(1)
                                 CALL LOAD(lx,2,1,polm,jj)
                                 CALL STING(jj)
                                 CALL PATH(jj)
                                 CALL INTG(IEXP)
                                 CALL TENB(j,bten,LMAX)
                              ENDDO
                              CALL TENS(bten)
                              CALL DECAY(ccd,0,ccc)
                              DO j = 1 , LP2 ! LP2 = 1500 (maximum number of matrix elements)
                                 DO ijan = 1 , 20
                                    SUMCL(ijan,j) = 0.
                                 ENDDO
                              ENDDO
                              ija0 = 0
                              DO ijan = 1 , jan ! For each detector angle
                                 IF ( IAXS(lx).EQ.0 ) nfi = 1
                                 DO jyi = 1 , idr
                                    GRAD(jyi) = 0.
                                 ENDDO
                                 todfi = 0.
                                 DO jfi = 1 , nfi ! For each phi angle
                                    fi0 = fiex1(ktt,jfi,1)/57.2957795
                                    fi1 = fiex1(ktt,jfi,2)/57.2957795
                                    gth = AGELI(IEXP,ijan,1)
                                    fm = (fi0+fi1)/2.
                                    figl = AGELI(IEXP,ijan,2)
                                    CALL ANGULA(YGN,idr,1,fi0,fi1,tetrc,
     &                                 gth,figl,ijan,op2)
                                    IF ( IFMO.NE.0 ) THEN ! If correction due to recoil
                                       id = ITMA(IEXP,ijan) ! Get detector identity
                                       d = ODL(id) ! Get result of OP,GDET calculation
                                       rx = d*SIN(gth)*COS(figl-fm)
     &                                    - .25*SIN(tetrc)*COS(fm)
                                       ry = d*SIN(gth)*SIN(figl-fm)
     &                                    - .25*SIN(tetrc)*SIN(fm)
                                       rz = d*COS(gth) - .25*COS(tetrc)
                                       rl = SQRT(rx*rx+ry*ry+rz*rz)
                                       sf = d*d/rl/rl
                                       thc = TACOS(rz/rl)
                                       fic = ATAN2(ry,rx)
                                       CALL ANGULA(YGP,idr,1,fi0,fi1,
     &                                    tetrc,thc,fic,ijan,op2)
                                       DO ixl = 1 , idr ! For each decay
                                         ixm = KSEQ(ixl,3)
                                         tfac = TAU(ixm)
                                         YGN(ixl) = YGN(ixl)
     &                                      + .01199182*tfac*BETAR(IEXP)
     &                                      *(sf*YGP(ixl)-YGN(ixl))
                                       ENDDO ! Loop on decays
                                    ENDIF ! If correction due to recoil
                                    IF ( IRAWEX(lx).NE.0 ) THEN
                                       ipd = ITMA(lx,ijan) ! Get identity of detector
                                       DO jyi = 1 , idr ! For each decay
                                         ni = KSEQ(jyi,3)
                                         nf = KSEQ(jyi,4)
                                         decen = EN(ni) - EN(nf)
                                         cocos = SIN(tetrc)*SIN(gth)
     &                                      *COS(fm-figl) + COS(tetrc)
     &                                      *COS(gth)
                                         decen = decen*(1.+BETAR(lx)
     &                                      *cocos)
                                         CALL EFFIX(ipd,decen,effi)
                                         YGN(jyi) = YGN(jyi)*effi
                                       ENDDO
                                       inclus = ICLUST(lx,ijan) ! Cluster number for detector ijan
                                       IF ( inclus.NE.0 ) THEN
                                         DO jyi = 1 , idr ! For each decay
                                         SUMCL(inclus,jyi)
     &                                      = SUMCL(inclus,jyi)
     &                                      + YGN(jyi)
                                         ENDDO
                                         IF ( ijan.NE.LASTCL(lx,inclus)
     &                                      ) GOTO 132 ! If it is not the last detector in the cluster
                                         DO jyi = 1 , idr ! For each decay
                                         YGN(jyi) = SUMCL(inclus,jyi)
                                         ENDDO
                                       ENDIF
                                    ENDIF
                                    IF ( jfi.EQ.1 ) ija0 = ija0 + 1
                                    DO jyi = 1 , idr ! For each decay
                                       GRAD(jyi) = GRAD(jyi) + YGN(jyi)
                                    ENDDO ! Loop on decays jyi
                                    todfi = todfi + ABS(fi1-fi0)
                                 ENDDO ! For each phi angle jfi
                                 IF ( IAXS(lx).EQ.0 ) todfi = 6.283185
                                 ax = 1.
                                 IF ( mfla.EQ.1 ) ax = 1./todfi
                                 dsx = dsig
                                 IF ( mfla.NE.1 ) dsx = dsig*todfi
                                 dsxm(mpin,kloop,ktt) = dsx
                                 WRITE (17,*) lx , mpin , kloop , ktt , 
     &                                  dsx
                                 WRITE (14,*) lx , enb , tting , ija0 , 
     &                                  dsx , 
     &                                  (GRAD(jyi)*dsig*ax,jyi=1,idr)
                                 IF ( IPRM(11).EQ.1 ) THEN
                                    WRITE (22,99048) lx , ija0 , enb , 
     &                                 tta
                                    IF ( tta.LT.0. ) WRITE (22,99017)
     &                                 tting
99017                               FORMAT (5X,
     &                             'RESPECTIVE TARGET SCATTERING ANGLE='
     &                             ,1F7.3,1X,'DEG'/)
                                    DO jyi = 1 , idr
                                       ni = KSEQ(jyi,3)
                                       nf = KSEQ(jyi,4)
                                       WRITE (22,99049) ni , nf , 
     &                                    SPIN(ni) , SPIN(nf) , 
     &                                    GRAD(jyi)*dsig*ax , GRAD(jyi)
     &                                    /GRAD(IDRN)
                                    ENDDO ! Loop on decays jyi
                                 ENDIF ! If printout of yields at meshpoints
 132                             CONTINUE
                              ENDDO ! Loop on detector angles ijan
                           ENDDO ! Loop on theta angles ktt
                        ENDDO ! Loop on energy meshpoints kloop
                     ENDDO ! Loop on pin diodes mpin
                      
                     EP(lx) = enh
                     TLBDG(lx) = tth
                  ENDDO ! Loop on experiments lx
                  REWIND 14
                  REWIND 15
                  iske = 0
                  DO na = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
                     ILE(na) = 1
                  ENDDO
                  ilx = 0
C                 We have now performed the full coulex calculation at each of the
C                 meshpoints, so now we start the integration
                  DO lx = 1 , NEXPT ! Loop over experiments
C                    Read tape 17
                     REWIND 17
                     DO ijaja = 1 , 300000
                        READ (17,*,END=134) jjlx , jmpin , jkloo , 
     &                        jktt , dsx
                        IF ( jjlx.EQ.lx ) dsxm(jmpin,jkloo,jktt) = dsx
                     ENDDO
 134                 na = NANG(lx)
                     IF ( lx.NE.1 ) THEN
                        DO na1 = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
                           ILE(na1) = ILE(na1) + NYLDE(lx-1,na1)
                        ENDDO
                     ENDIF
                     READ (JZB,*) nptx ! Number of meshpoints for stopping powers
                     IF ( nptx.NE.0 ) THEN
                        READ (JZB,*) (esp(i),i=1,nptx) ! Energy
                        READ (JZB,*) (dedx(i),i=1,nptx) ! Stopping power
                        npt = nptx
                     ENDIF
                     READ (JZB,*) npce , npct
                     mfla = 0
                     IF ( npct.LT.0 ) mfla = 1
                     IF ( iecd(lx).EQ.1 ) mfla = 1
                     npct = ABS(npct)
                     IF ( npct.GT.100 )
     &                  STOP 'ABS(NI2) is limited to 100!'
                     npce = npce + MOD(npce,2)
                     npct = npct + MOD(npct,2)
                     mpin = 1
                     IF ( ipinf.NE.0 ) THEN
                        IF ( jpin(lx).NE.0 ) mpin = jpin(lx)
                     ENDIF
                     dst = 0.
                     DO lpin = 1 , mpin ! Loop over pin diodes
                        ilx = ilx + 1
                        IF ( ilx.NE.1 )
     &                       CALL TAPMA(lx,iske,isko,iskf,nflr,idr,0,
     &                       nft,enb)
                        READ (14,*) ne , ntt , emn , emx , tmn , tmx , 
     &                              jan , wth , wph , wthh
                        iocc = (ne+ntt)*idr
                        IF ( iocc.GT.izcap ) GOTO 1800
                        hen = (emx-emn)/npce
                        npce1 = npce + 1
                        het = (tmx-tmn)/npct ! Step in theta in degrees
                        npct1 = npct + 1
                        IF ( iecd(lx).EQ.1 ) ! Circular detector
     &                       CALL COORD(wth,wph,wthh,npct1,1,pfi,wpi,
     &                       TLBDG(lx),lx,tmn,tmx)
                        IF ( iecd(lx).NE.1 ) THEN
                           IF ( mfla.EQ.1 ) READ (JZB,*)
     &                          (pfi(j),j=1,npct1)
                        ENDIF
                        het = het/57.2957795 ! Step in theta in radians
                        
C                       Interpolate stopping power for each of the energies
C                       that we need. esp is an array of energies and dedx is
C                       an array containing the stopping powers at those
C                       energies. Function is unweighted sqrt. The energies
C                       are not the energies we gave for the meshpoints, but
C                       the range over which we integrate the bombarding energy
C                       with the number of steps specified.
                        DO j = 1 , npce1
                           xx = (j-1)*hen + emn
                           IF ( ISPL.EQ.0 )
     &                        CALL LAGRAN(esp,dedx,npt,1,xx,yy,3,1)
                           IF ( ISPL.EQ.1 )
     &                        CALL SPLNER(esp,dedx,npt,xx,yy,3)
                           HLMLM(j) = 1./yy
                        ENDDO
                         
C                       Now we calculate for all the mesh points. 
                        naa = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) naa = NANG(lx)
                        iskf = naa - 1
                        DO ja = 1 , naa ! Loop over detector angles
                           icll = 3 ! Weighting mode
                           DO je = 1 , ne ! ne = number of energy mesh points
                              lu = ILE(ja)
                              isko = (je-1)*naa*ntt + ja - 1
                              CALL TAPMA(lx,iske,isko,iskf,ntt,idr,1,
     &                           nft,enb)
                              IF ( nft.EQ.1 ) GOTO 1900 ! Troubleshoot
                              DO jd = 1 , idr ! For each decay
                                 DO jtp = 1 , ntt ! ntt = number of theta meshpoints
                                    IF ( jd.EQ.1 .AND. ja.EQ.1 )
     &                                 DSG(jtp) = dsxm(lpin,je,jtp)
                                    jyv = (jtp-1)*idr + jd
                                    YV(jtp) = ZETA(jyv) ! Point yield
                                 ENDDO ! Loop on theta meshpoints jtp
                                 DO jt = 1 , npct1 ! number of equal divisions in theta for interpolation
                                    xx = (jt-1)*het + tmn/57.2957795
                                    IF ( ISPL.EQ.0 )
     &                                 CALL LAGRAN(XV,YV,ntt,jt,xx,yy,2,
     &                                 icll) ! interpolate point yield at theta = xx
                                    IF ( ISPL.EQ.1 )
     &                                 CALL SPLNER(XV,YV,ntt,xx,yy,2) ! interpolate point yield at theta = xx
                                    IF ( ISPL.EQ.0 )
     &                                 CALL LAGRAN(XV,DSG,ntt,jt,xx,zz,
     &                                 2,icll) ! interpolate gamma yield at theta = xx
                                    IF ( ISPL.EQ.1 )
     &                                 CALL SPLNER(XV,DSG,ntt,xx,zz,
     &                                 2) ! interpolate gamma yield at theta = xx
                                    IF ( mfla.EQ.1 ) yy = yy*pfi(jt)
     &                                 /57.2957795
                                    IF ( yy.LE.0. ) yy = 1.E-15
                                    IF ( mfla.EQ.1 ) zz = zz*pfi(jt)
     &                                 /57.2957795
                                    XI(jt) = yy*SIN(xx) ! yy = integral of point yields over phi
                                    IF ( jd.EQ.1 .AND. ja.EQ.1 ) HLM(jt)
     &                                 = zz*SIN(xx) ! zz = integral over phi of Rutherford cross section
                                 ENDDO ! Loop on equal theta divisions jt
                                 icll = 4
                                 locat = ntt*idr + (je-1)*idr + jd
C                                Integrate point yields over theta using Simpson's rule
                                 ZETA(locat) = SIMIN(npct1,het,XI)
C                                If it is first decay and angle, integrate Rutherford cross section over theta
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 ) DSE(je)
     &                                = SIMIN(npct1,het,HLM)
                                 ZV(je) = enb
                              ENDDO ! Loop on decays jd
                           ENDDO ! Loop on energy meshpoints je

C    Interpolation over energy:
C    The array ZV contains the energies of the meshpoints and the elements of the YV
C    array are set to the angle-integrated yield for each decay at the corresponding
C    energy, while DSE contains the Rutherford cross section for those energies. Since
C    the energies of the meshpoints are not necessarily equally spaced, we need to
C    interpolate to a set of equally spaced energies separated by "hen" starting from
C    "emn". To get the contribution from each energy, dE = 1 / (stopping power). Note
C    that we only evaluate the Rutherford cross section for the first decay and first
C    angle, since it is the same for all.

                           icll = 3
                           DO jd = 1 , idr ! For each decay
                              DO jtp = 1 , ne ! For each energy meshpoint
                                 jyv = (jtp-1)*idr + jd + ntt*idr
                                 YV(jtp) = ZETA(jyv)
                              ENDDO ! Loop on energy meshpoints jtp
                              DO jt = 1 , npce1 ! npce1 is number of equal energy steps
                                 xx = (jt-1)*hen + emn

C                                Interpolate the angle-integrated yield for this energy
                                 IF ( ISPL.EQ.0 )
     &                                CALL LAGRAN(ZV,YV,ne,jt,xx,yy,2,
     &                                icll)
                                 IF ( ISPL.EQ.1 )
     &                                CALL SPLNER(ZV,YV,ne,xx,yy,2)

C                                Interpolate Rutherford cross-section for this energy
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 .AND. ! Only for first decay and angle
     &                                ISPL.EQ.0 )
     &                                CALL LAGRAN(ZV,DSE,ne,jt,xx,zz,2,
     &                                icll) ! Interpolate for this energy
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 .AND.
     &                                ISPL.EQ.1 )
     &                                CALL SPLNER(ZV,DSE,ne,xx,zz,2) ! Interpolate for this energy
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 ) HLM(jt)
     &                             = zz*HLMLM(jt) ! HLMLM = 1 / stopping power
                                 XI(jt) = yy*HLMLM(jt)
                              ENDDO ! Loop on equal energy steps

C   So now after this loop, we have XI containing the angle-integrated yield times dE for 
C   a set of equally spaced energies, so we use Simpson's rule to integrate them and store
C   in GRAD(jd). The first time, we also have in HLM a set of Rutherford cross-sections for
C   equally spaced energies, which we integrate in the same way.
                              icll = 4
                              IF ( jd.EQ.1 .AND. ja.EQ.1 )
     &                             DS = SIMIN(npce1,hen,HLM) ! integrate
                              GRAD(jd) = SIMIN(npce1,hen,XI)
                           ENDDO ! Loop over decays jd

                           IF ( ja.EQ.1 ) dst = dst + DS
                           IF ( ja.EQ.1 ) WRITE (22,99018) DS , lx
99018                      FORMAT (1X/////5X,
     &                            'INTEGRATED RUTHERFORD CROSS SECTION='
     &                            ,1E9.4,2X,'FOR EXP.',1I2///)

                           WRITE (22,99019) lx , ja , emn , emx , tmn , 
     &                            tmx
99019                      FORMAT (1X,//50X,'INTEGRATED YIELDS'//5X,
     &                             'EXPERIMENT ',1I2,2X,'DETECTOR ',
     &                             1I2/5X,'ENERGY RANGE ',1F8.3,'---',
     &                             1F8.3,1X,'MEV',3X,
     &                             'SCATTERING ANGLE RANGE ',1F7.3,
     &                             '---',1F7.3,1X,'DEG'//5X,'NI',5X,
     &                             'NF',5X,'II',5X,'IF',5X,'YIELD',5X,
     &                             'NORMALIZED YIELD'/)
                           DO jd = 1 , idr
                              WRITE (15,*) GRAD(jd)
                           ENDDO
                           DO jd = 1 , idr
                              ni = KSEQ(jd,3)
                              nf = KSEQ(jd,4)
                              WRITE (22,99049) ni , nf , SPIN(ni) , 
     &                               SPIN(nf) , GRAD(jd) , GRAD(jd)
     &                               /GRAD(IDRN) ! IDRN is the normalising transition
                           ENDDO
                        ENDDO ! Loop over detector angles ja

                        IF ( iecd(lx).EQ.1 ) THEN ! Circular detector
                           IF ( jpin(lx).EQ.0 ) THEN
                              CALL COORD(wth,wph,wthh,1,2,pfi,wpi,
     &                           TLBDG(lx),lx,txx,txx)
                              WRITE (22,99020) FIEX(lx,1)*57.2957795 , 
     &                               FIEX(lx,2)*57.2957795 , lx
99020                         FORMAT (//5X,
     &                          'WARNING: THE PHI ANGLE WAS REPLACED BY'
     &                          ,1X,F8.3,1X,'TO',F8.3,3X,
     &                          'FOR EXPERIMENT',2X,I3)
                              IF ( TLBDG(lx).LT.0 ) THEN
                                 FIEX(lx,1) = FIEX(lx,1) + 3.14159265
                                 FIEX(lx,2) = FIEX(lx,2) + 3.14159265
                              ENDIF ! If theta_lab < 0
                           ENDIF ! If no pin diodes
                        ENDIF ! If circular detector
                        iske = iske + ne*ntt*naa
                     ENDDO ! Loop over pin diodes
                     IF ( mpin.GT.1 ) WRITE (22,99021) dst , lx
99021                FORMAT (1x//2x,
     &                      'Total integrated Rutherford cross section='
     &                      ,1E8.3,' for exp. ',1I2/)
                  ENDDO
                  REWIND 17 ! Added PJN (17Jul2009)
                  IF ( ipinf.NE.0 ) THEN
                     ngpr = 0
                     DO lx = 1 , NEXPT ! For each experiment
                        nged = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) nged = NANG(lx)
                        IF ( lx.NE.1 ) ngpr = ngpr + idr*jpin(lx-1)
     &                       *NDST(lx-1)
                        lpin = jpin(lx)
                        IF ( lpin.EQ.0 ) lpin = 1
                        DO jgd = 1 , nged ! For each angle or dataset
                           DO jd = 1 , idr
                              GRAD(jd) = 0.
                           ENDDO
                           DO mpin = 1 , lpin ! For each pin diode
                              REWIND 15
                              ndum = ngpr + (jgd-1)*idr + (mpin-1)
     &                          *nged*idr ! Was jgd instead of nged (PJN 17Jul2009)
                              IF ( ndum.NE.0 ) THEN
                                 DO jd = 1 , ndum
                                    READ (15,*) xx
                                 ENDDO
                              ENDIF
                              DO jd = 1 , idr ! For each decay
                                 READ (15,*) xx
                                 GRAD(jd) = GRAD(jd) + xx
                              ENDDO ! Loop on decays jd
                           ENDDO ! Loop on pin diodes mpin
                           WRITE (17,*) (GRAD(jd),jd=1,idr)
                        ENDDO ! Loop on angle or dataset jgd
                     ENDDO ! Loop on experiment lx
                     REWIND 15
                     REWIND 17
                     DO lx = 1 , NEXPT ! For each experiment
                        nged = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) nged = NANG(lx)
                        DO ija0 = 1 , nged ! For each angle or dataset
                           READ (17,*) (GRAD(jdy),jdy=1,idr)
                           DO jd = 1 , idr ! For each decay
                              WRITE (15,*) GRAD(jd)
                           ENDDO ! Loop on decays jd
                        ENDDO ! Loop on angle or dataset ija0
                     ENDDO ! Loop on experiments lx
                  ENDIF
                  GOTO 100 ! End of OP,INTG - back to input loop

C              Treat OP,INTI
               ELSEIF ( op2.EQ.'INTI' ) THEN
                  DO lx = 1 , NEXPT ! For each experiment store original ISKIN
                     iskin_protect(lx) = ISKIN(lx)
                  ENDDO
                  REWIND 14
                  lfagg = 1
                  IF ( SPIN(1).LT..25 ) ISO = 0
                  DO lx = 1 , NEXPT ! For each experiment
                     lpin = 1
                     IF ( ipinf.NE.0 ) THEN
                        IF ( jpin(lx).NE.0 ) lpin = jpin(lx)
                     ENDIF
                     IEXP = lx
                     tth = TLBDG(lx)
                     enh = EP(lx)
                     DO mpin = 1 , lpin ! For each pin diode
                        IF ( iecd(lx).EQ.1 ) THEN ! Circular detector
                           READ (JZB,*) ne , ntt , emn , emx , wth , 
     &                          wph , wthh
                           mfla = 1
                           CALL COORD(wth,wph,wthh,ntt,0,pfi,wpi,tth,lx,
     &                                tmn,tmx)
                        ELSE
                           READ (JZB,*) ne , ntt , emn , emx , tmn , tmx
                           mfla = 0
                           IF ( ntt.LT.0 ) mfla = 1
                        ENDIF
                        ntt = ABS(ntt)
                        jan = NANG(lx)
                        jan1 = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) jan1 = jan
                        IF ( iecd(lx).EQ.1 ) THEN ! Circular detector
                           WRITE (14,*) ne , ntt , emn , emx , tmn , 
     &                                  tmx , jan1 , wth , wph , wthh
                        ELSE
                           WRITE (14,*) ne , ntt , emn , emx , tmn , 
     &                                  tmx , jan1 , tmx , tmx , tmx
                        ENDIF
                        READ (JZB,*) (XV(i),i=1,ne)
                        IF ( iecd(lx).NE.1 ) READ (JZB,*)
     &                       (YV(i),i=1,ntt)
                        IF ( tth.LT.0. ) ELMH(2*lx-1) = YV(1)
                        IF ( tth.LT.0. ) ELMH(2*lx) = YV(ntt)
                        DO kloop = 1 , ne ! For each energy meshpoint
                           enb = XV(kloop)
                           EP(lx) = enb
                           DO ktt = 1 , ntt
                              tta = YV(ktt)
                              IF ( tth.LT.0 )
     &                           CALL INVKIN(EP(lx),EN(NCM),IZ1(lx),
     &                                       XA,XA1(lx),YV(ktt),tta,
     &                                       1,ISKIN(lx))
                              tta = SIGN(tta, tth)
                              IF ( IAXS(lx).NE.0 ) THEN ! If not axial symmetry
                                 IF ( iecd(lx).NE.1 ) THEN
                                    IF ( kloop.EQ.1 ) THEN
                                       READ (JZB,*) nfi ! Number of phi ranges
                                       READ (JZB,*) 
     &                                    (fiex1(ktt,jfi,1),fiex1(ktt,
     &                                    jfi,2),jfi=1,nfi)
                                       IF ( tth.LT.0. ) THEN
                                         DO jfi = 1 , nfi ! For each phi angle
                                         fiex1(ktt,jfi,1)
     &                                      = fiex1(ktt,jfi,1) + 180.
                                         fiex1(ktt,jfi,2)
     &                                      = fiex1(ktt,jfi,2) + 180.
                                         ENDDO
                                       ENDIF
                                    ENDIF
                                 ENDIF
                              ENDIF ! If not axial symmetry
                              TLBDG(lx) = tta
                              IF ( kloop.EQ.1 ) THEN
                                 IF ( iecd(lx).NE.0 ) THEN
                                    nfi = 1
                                    fiex1(ktt,1,1) = wpi(ktt,1) ! Lower phi limit
                                    fiex1(ktt,1,2) = wpi(ktt,2) ! Upper phi limit
                                 ENDIF
                              ENDIF
                              CALL CMLAB(lx,dsig,tetrc)
                              IF ( ERR ) GOTO 2000 ! Error
                              tting = TLBDG(lx)
                              IF ( ERR ) GOTO 1900 ! Troubleshoot
                              CALL LOAD(lx,1,1,0.D0,jj)
                              CALL ALLOC(ACCUR)
                              CALL SNAKE(lx,ZPOL)
                              CALL SETIN
                              DO j = 1 , LMAX ! For each spin up to ground-state spin + 1
                                 polm = DBLE(j-1) - SPIN(1)
                                 CALL LOAD(lx,2,1,polm,jj)
                                 CALL STING(jj)
                                 CALL PATH(jj)
                                 CALL INTG(IEXP)
                                 CALL TENB(j,bten,LMAX)
                              ENDDO
                              CALL TENS(bten)
                              CALL DECAY(ccd,0,ccc)
                              DO j = 1 , LP2 ! LP2 = 1500 (maximum number of matrix elements)
                                 DO ijan = 1 , 20
                                    SUMCL(ijan,j) = 0.
                                 ENDDO
                              ENDDO
                              ija0 = 0
                              DO ijan = 1 , jan ! For each detector angle
                                 IF ( IAXS(lx).EQ.0 ) nfi = 1
                                 DO jyi = 1 , idr
                                    GRAD(jyi) = 0.
                                 ENDDO
                                 todfi = 0.
                                 DO jfi = 1 , nfi ! For each phi angle
                                    fi0 = fiex1(ktt,jfi,1)/57.2957795
                                    fi1 = fiex1(ktt,jfi,2)/57.2957795
                                    gth = AGELI(IEXP,ijan,1)
                                    fm = (fi0+fi1)/2.
                                    figl = AGELI(IEXP,ijan,2)
                                    CALL ANGULA(YGN,idr,1,fi0,fi1,tetrc,
     &                                 gth,figl,ijan,op2)
                                    IF ( IFMO.NE.0 ) THEN ! If correction due to recoil
                                       id = ITMA(IEXP,ijan) ! Get detector identity
                                       d = ODL(id) ! Get result of OP,GDET calculation
                                       rx = d*SIN(gth)*COS(figl-fm)
     &                                    - .25*SIN(tetrc)*COS(fm)
                                       ry = d*SIN(gth)*SIN(figl-fm)
     &                                    - .25*SIN(tetrc)*SIN(fm)
                                       rz = d*COS(gth) - .25*COS(tetrc)
                                       rl = SQRT(rx*rx+ry*ry+rz*rz)
                                       sf = d*d/rl/rl
                                       thc = TACOS(rz/rl)
                                       fic = ATAN2(ry,rx)
                                       CALL ANGULA(YGP,idr,1,fi0,fi1,
     &                                    tetrc,thc,fic,ijan,op2)
                                       DO ixl = 1 , idr ! For each decay
                                         ixm = KSEQ(ixl,3)
                                         tfac = TAU(ixm)
                                         YGN(ixl) = YGN(ixl)
     &                                      + .01199182*tfac*BETAR(IEXP)
     &                                      *(sf*YGP(ixl)-YGN(ixl))
                                       ENDDO ! Loop on decays
                                    ENDIF ! If correction due to recoil
                                    IF ( IRAWEX(lx).NE.0 ) THEN
                                       ipd = ITMA(lx,ijan) ! Get identity of detector
                                       DO jyi = 1 , idr ! For each decay
                                         ni = KSEQ(jyi,3)
                                         nf = KSEQ(jyi,4)
                                         decen = EN(ni) - EN(nf)
                                         cocos = SIN(tetrc)*SIN(gth)
     &                                      *COS(fm-figl) + COS(tetrc)
     &                                      *COS(gth)
                                         decen = decen*(1.+BETAR(lx)
     &                                      *cocos)
                                         CALL EFFIX(ipd,decen,effi)
                                         YGN(jyi) = YGN(jyi)*effi
                                       ENDDO
                                       inclus = ICLUST(lx,ijan) ! Cluster number for detector ijan
                                       IF ( inclus.NE.0 ) THEN
                                         DO jyi = 1 , idr ! For each decay
                                         SUMCL(inclus,jyi)
     &                                      = SUMCL(inclus,jyi)
     &                                      + YGN(jyi)
                                         ENDDO
                                         IF ( ijan.NE.LASTCL(lx,inclus)
     &                                      ) GOTO 432 ! If it is not the last detector in the cluster
                                         DO jyi = 1 , idr ! For each decay
                                         YGN(jyi) = SUMCL(inclus,jyi)
                                         ENDDO
                                       ENDIF
                                    ENDIF
                                    IF ( jfi.EQ.1 ) ija0 = ija0 + 1
                                    DO jyi = 1 , idr ! For each decay
                                       GRAD(jyi) = GRAD(jyi) + YGN(jyi)
                                    ENDDO ! Loop on decays jyi
                                    todfi = todfi + ABS(fi1-fi0)
                                 ENDDO ! For each phi angle jfi
                                 IF ( IAXS(lx).EQ.0 ) todfi = 6.283185
                                 ax = 1.
                                 IF ( mfla.EQ.1 ) ax = 1./todfi
                                 dsx = dsig
                                 IF ( mfla.NE.1 ) dsx = dsig*todfi
                                 dsxm(mpin,kloop,ktt) = dsx
                                 WRITE (17,*) lx , mpin , kloop , ktt , 
     &                                  dsx
                                 WRITE (14,*) lx , enb , tting , ija0 , 
     &                                  dsx , 
     &                                  (GRAD(jyi)*dsig*ax,jyi=1,idr)
                                 IF ( IPRM(11).EQ.1 ) THEN
                                    WRITE (22,99048) lx , ija0 , enb , 
     &                                 tta
                                    IF ( tta.LT.0. ) WRITE (22,99017)
     &                                 tting
                                    DO jyi = 1 , idr
                                       ni = KSEQ(jyi,3)
                                       nf = KSEQ(jyi,4)
                                       WRITE (22,99049) ni , nf , 
     &                                    SPIN(ni) , SPIN(nf) , 
     &                                    GRAD(jyi)*dsig*ax , GRAD(jyi)
     &                                    /GRAD(IDRN)
                                    ENDDO ! Loop on decays jyi
                                 ENDIF ! If printout of yields at meshpoints
 432                             CONTINUE
                              ENDDO ! Loop on detector angles ijan
                           ENDDO ! Loop on theta angles ktt
                        ENDDO ! Loop on energy meshpoints kloop
                     ENDDO ! Loop on pin diodes mpin
                      
                     EP(lx) = enh
                     TLBDG(lx) = tth
                  ENDDO ! Loop on experiments lx
                  REWIND 14
                  REWIND 15
                  iske = 0
                  DO na = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
                     ILE(na) = 1
                  ENDDO
                  ilx = 0
C                 We have now performed the full coulex calculation at each of the
C                 meshpoints, so now we start the integration
                  DO lx = 1 , NEXPT ! Loop over experiments
C                    Read tape 17
                     REWIND 17
                     DO ijaja = 1 , 300000
                        READ (17,*,END=434) jjlx , jmpin , jkloo , 
     &                        jktt , dsx
                        IF ( jjlx.EQ.lx ) dsxm(jmpin,jkloo,jktt) = dsx
                     ENDDO
 434                 na = NANG(lx)
                     IF ( lx.NE.1 ) THEN
                        DO na1 = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
                           ILE(na1) = ILE(na1) + NYLDE(lx-1,na1)
                        ENDDO
                     ENDIF
                     READ (JZB,*) nptx ! Number of meshpoints for stopping powers
                     IF ( nptx.NE.0 ) THEN
                        READ (JZB,*) (esp(i),i=1,nptx) ! Energy
                        READ (JZB,*) (dedx(i),i=1,nptx) ! Stopping power
                        npt = nptx
                     ENDIF
                     READ (JZB,*) npce , npct
                     mfla = 0
                     IF ( npct.LT.0 ) mfla = 1
                     IF ( iecd(lx).EQ.1 ) mfla = 1
                     npct = ABS(npct)
                     IF ( npct.GT.100 )
     &                  STOP 'ABS(NI2) is limited to 100!'
                     npce = npce + MOD(npce,2)
                     npct = npct + MOD(npct,2)
                     mpin = 1
                     IF ( ipinf.NE.0 ) THEN
                        IF ( jpin(lx).NE.0 ) mpin = jpin(lx)
                     ENDIF
                     dst = 0.
                     DO lpin = 1 , mpin ! Loop over pin diodes
                        ilx = ilx + 1
                        IF ( ilx.NE.1 )
     &                       CALL TAPMA(lx,iske,isko,iskf,nflr,idr,0,
     &                       nft,enb)
                        READ (14,*) ne , ntt , emn , emx , tmn , tmx , 
     &                              jan , wth , wph , wthh
                        iocc = (ne+ntt)*idr
                        IF ( iocc.GT.izcap ) GOTO 1800
                        hen = (emx-emn)/npce
                        npce1 = npce + 1
                        het = (tmx-tmn)/npct ! Step in theta in degrees
                        npct1 = npct + 1
                        IF ( iecd(lx).EQ.1 ) ! Circular detector
     &                       CALL COORD(wth,wph,wthh,npct1,1,pfi,wpi,
     &                       TLBDG(lx),lx,tmn,tmx)
                        IF ( iecd(lx).NE.1 ) THEN
                           IF ( mfla.EQ.1 ) READ (JZB,*)
     &                          (pfi(j),j=1,npct1)
                        ENDIF
                        het = het/57.2957795 ! Step in theta in radians
                        
C                       Interpolate stopping power for each of the energies
C                       that we need. esp is an array of energies and dedx is
C                       an array containing the stopping powers at those
C                       energies. Function is unweighted sqrt. The energies
C                       are not the energies we gave for the meshpoints, but
C                       the range over which we integrate the bombarding energy
C                       with the number of steps specified.
                        DO j = 1 , npce1
                           xx = (j-1)*hen + emn
                           IF ( ISPL.EQ.0 )
     &                        CALL LAGRAN(esp,dedx,npt,1,xx,yy,3,1)
                           IF ( ISPL.EQ.1 )
     &                        CALL SPLNER(esp,dedx,npt,xx,yy,3)
                           HLMLM(j) = 1./yy
                        ENDDO
                         
C                       Now we calculate for all the mesh points. 
                        naa = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) naa = NANG(lx)
                        iskf = naa - 1
                        DO ja = 1 , naa ! Loop over detector angles
                           icll = 3 ! Weighting mode
                           DO je = 1 , ne ! ne = number of energy mesh points
                              lu = ILE(ja)
                              isko = (je-1)*naa*ntt + ja - 1
                              CALL TAPMA(lx,iske,isko,iskf,ntt,idr,1,
     &                           nft,enb)
                              IF ( nft.EQ.1 ) GOTO 1900 ! Troubleshoot
                              DO jd = 1 , idr ! For each decay
                                 DO jtp = 1 , ntt ! ntt = number of theta meshpoints
                                    IF ( jd.EQ.1 .AND. ja.EQ.1 )
     &                                 DSG(jtp) = dsxm(lpin,je,jtp)
                                    jyv = (jtp-1)*idr + jd
                                    YV(jtp) = ZETA(jyv) ! Point yield
                                 ENDDO ! Loop on theta meshpoints jtp
                                 DO jt = 1 , npct1 ! number of equal divisions in theta for interpolation
                                    xx = (jt-1)*het + tmn/57.2957795
                                    IF ( ISPL.EQ.0 )
     &                                 CALL LAGRAN(XV,YV,ntt,jt,xx,yy,2,
     &                                 icll) ! interpolate point yield at theta = xx
                                    IF ( ISPL.EQ.1 )
     &                                 CALL SPLNER(XV,YV,ntt,xx,yy,2) ! interpolate point yield at theta = xx
                                    IF ( ISPL.EQ.0 )
     &                                 CALL LAGRAN(XV,DSG,ntt,jt,xx,zz,
     &                                 2,icll) ! interpolate gamma yield at theta = xx
                                    IF ( ISPL.EQ.1 )
     &                                 CALL SPLNER(XV,DSG,ntt,xx,zz,
     &                                 2) ! interpolate gamma yield at theta = xx
                                    IF ( mfla.EQ.1 ) yy = yy*pfi(jt)
     &                                 /57.2957795
                                    IF ( yy.LE.0. ) yy = 1.E-15
                                    IF ( mfla.EQ.1 ) zz = zz*pfi(jt)
     &                                 /57.2957795
                                    XI(jt) = yy*SIN(xx) ! yy = integral of point yields over phi
                                    IF ( jd.EQ.1 .AND. ja.EQ.1 ) HLM(jt)
     &                                 = zz*SIN(xx) ! zz = integral over phi of Rutherford cross section
                                 ENDDO ! Loop on equal theta divisions jt
                                 icll = 4
                                 locat = ntt*idr + (je-1)*idr + jd
C                                Integrate point yields over theta using Simpson's rule
                                 ZETA(locat) = SIMIN(npct1,het,XI)
C                                If it is first decay and angle, integrate Rutherford cross section over theta
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 ) DSE(je)
     &                                = SIMIN(npct1,het,HLM)
                                 ZV(je) = enb
                              ENDDO ! Loop on decays jd
                           ENDDO ! Loop on energy meshpoints je

C    Interpolation over energy:
C    The array ZV contains the energies of the meshpoints and the elements of the YV
C    array are set to the angle-integrated yield for each decay at the corresponding
C    energy, while DSE contains the Rutherford cross section for those energies. Since
C    the energies of the meshpoints are not necessarily equally spaced, we need to
C    interpolate to a set of equally spaced energies separated by "hen" starting from
C    "emn". To get the contribution from each energy, dE = 1 / (stopping power). Note
C    that we only evaluate the Rutherford cross section for the first decay and first
C    angle, since it is the same for all.

                           icll = 3
                           DO jd = 1 , idr ! For each decay
                              DO jtp = 1 , ne ! For each energy meshpoint
                                 jyv = (jtp-1)*idr + jd + ntt*idr
                                 YV(jtp) = ZETA(jyv)
                              ENDDO ! Loop on energy meshpoints jtp
                              DO jt = 1 , npce1 ! npce1 is number of equal energy steps
                                 xx = (jt-1)*hen + emn

C                                Interpolate the angle-integrated yield for this energy
                                 IF ( ISPL.EQ.0 )
     &                                CALL LAGRAN(ZV,YV,ne,jt,xx,yy,2,
     &                                icll)
                                 IF ( ISPL.EQ.1 )
     &                                CALL SPLNER(ZV,YV,ne,xx,yy,2)

C                                Interpolate Rutherford cross-section for this energy
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 .AND. ! Only for first decay and angle
     &                                ISPL.EQ.0 )
     &                                CALL LAGRAN(ZV,DSE,ne,jt,xx,zz,2,
     &                                icll) ! Interpolate for this energy
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 .AND.
     &                                ISPL.EQ.1 )
     &                                CALL SPLNER(ZV,DSE,ne,xx,zz,2) ! Interpolate for this energy
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 ) HLM(jt)
     &                             = zz*HLMLM(jt) ! HLMLM = 1 / stopping power
                                 XI(jt) = yy*HLMLM(jt)
                              ENDDO ! Loop on equal energy steps

C   So now after this loop, we have XI containing the angle-integrated yield times dE for 
C   a set of equally spaced energies, so we use Simpson's rule to integrate them and store
C   in GRAD(jd). The first time, we also have in HLM a set of Rutherford cross-sections for
C   equally spaced energies, which we integrate in the same way.
                              icll = 4
                              IF ( jd.EQ.1 .AND. ja.EQ.1 )
     &                             DS = SIMIN(npce1,hen,HLM) ! integrate
                              GRAD(jd) = SIMIN(npce1,hen,XI)
                           ENDDO ! Loop over decays jd

                           IF ( ja.EQ.1 ) dst = dst + DS
                           IF ( ja.EQ.1 ) WRITE (22,99018) DS , lx

                           WRITE (22,99019) lx , ja , emn , emx , tmn , 
     &                            tmx
                           DO jd = 1 , idr
                              WRITE (15,*) GRAD(jd)
                           ENDDO
                           DO jd = 1 , idr
                              ni = KSEQ(jd,3)
                              nf = KSEQ(jd,4)
                              WRITE (22,99049) ni , nf , SPIN(ni) , 
     &                               SPIN(nf) , GRAD(jd) , GRAD(jd)
     &                               /GRAD(IDRN) ! IDRN is the normalising transition
                           ENDDO
                        ENDDO ! Loop over detector angles ja

                        IF ( iecd(lx).EQ.1 ) THEN ! Circular detector
                           IF ( jpin(lx).EQ.0 ) THEN
                              CALL COORD(wth,wph,wthh,1,2,pfi,wpi,
     &                           TLBDG(lx),lx,txx,txx)
                              WRITE (22,99020) FIEX(lx,1)*57.2957795 , 
     &                               FIEX(lx,2)*57.2957795 , lx
                              IF ( TLBDG(lx).LT.0 ) THEN
                                 FIEX(lx,1) = FIEX(lx,1) + 3.14159265
                                 FIEX(lx,2) = FIEX(lx,2) + 3.14159265
                              ENDIF ! If theta_lab < 0
                           ENDIF ! If no pin diodes
                        ENDIF ! If circular detector
                        iske = iske + ne*ntt*naa
                     ENDDO ! Loop over pin diodes
                     IF ( mpin.GT.1 ) WRITE (22,99021) dst , lx
                  ENDDO
                  REWIND 17 ! Added PJN (17Jul2009)
                  IF ( ipinf.NE.0 ) THEN
                     ngpr = 0
                     DO lx = 1 , NEXPT ! For each experiment
                        nged = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) nged = NANG(lx)
                        IF ( lx.NE.1 ) ngpr = ngpr + idr*jpin(lx-1)
     &                       *NDST(lx-1)
                        lpin = jpin(lx)
                        IF ( lpin.EQ.0 ) lpin = 1
                        DO jgd = 1 , nged ! For each angle or dataset
                           DO jd = 1 , idr
                              GRAD(jd) = 0.
                           ENDDO
                           DO mpin = 1 , lpin ! For each pin diode
                              REWIND 15
                              ndum = ngpr + (jgd-1)*idr + (mpin-1)
     &                          *nged*idr ! Was jgd instead of nged (PJN 17Jul2009)
                              IF ( ndum.NE.0 ) THEN
                                 DO jd = 1 , ndum
                                    READ (15,*) xx
                                 ENDDO
                              ENDIF
                              DO jd = 1 , idr ! For each decay
                                 READ (15,*) xx
                                 GRAD(jd) = GRAD(jd) + xx
                              ENDDO ! Loop on decays jd
                           ENDDO ! Loop on pin diodes mpin
                           WRITE (17,*) (GRAD(jd),jd=1,idr)
                        ENDDO ! Loop on angle or dataset jgd
                     ENDDO ! Loop on experiment lx
                     REWIND 15
                     REWIND 17
                     DO lx = 1 , NEXPT ! For each experiment
                        nged = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) nged = NANG(lx)
                        DO ija0 = 1 , nged ! For each angle or dataset
                           READ (17,*) (GRAD(jdy),jdy=1,idr)
                           DO jd = 1 , idr ! For each decay
                              WRITE (15,*) GRAD(jd)
                           ENDDO ! Loop on decays jd
                        ENDDO ! Loop on angle or dataset ija0
                     ENDDO ! Loop on experiments lx
                  ENDIF
                  DO lx = 1 , NEXPT ! For each experiment restore original ISKIN
                     ISKIN(lx) = iskin_protect(lx)
                  ENDDO
                  GOTO 100 ! End of OP,INTI - back to input loop

C              Treat OP,CORR
               ELSEIF ( op2.EQ.'CORR' ) THEN
                  CALL READY(idr,ntap,0)
                  REWIND 3
                  REWIND 15
                  REWIND 4
                  GOTO 1200 ! End of OP,CORR
               ELSE

C                 Treat OP,POIN
                  IF ( op2.EQ.'POIN' ) GOTO 1200

C                 Treat OP,MAP
                  IF ( op2.EQ.'MAP ' ) iobl = 1

C                 Treat OP,STAR
                  IF ( op2.EQ.'STAR' ) GOTO 1200

C                 Treat OP,SIXJ
                  IF ( op2.EQ.'SIXJ' ) THEN
                     DO k = 1 , 2
                        l = 4*k
                        DO j = 1 , 80
                           ixj = j - 1
                           DO ms = 1 , 5
                              mend = 2*(ms-3) + ixj
                              WRITE (14,*) WSIXJ(l,4,4,ixj,mend,ixj-4) , 
     &                               WSIXJ(l,4,4,ixj,mend,ixj-2) , 
     &                               WSIXJ(l,4,4,ixj,mend,ixj) , 
     &                               WSIXJ(l,4,4,ixj,mend,ixj+2) , 
     &                               WSIXJ(l,4,4,ixj,mend,ixj+4)
                           ENDDO
                        ENDDO
                     ENDDO
                     GOTO 2000 ! End of OP,SIXJ - normal end of execution

C                 Treat OP,RAW (raw uncorrected gamma yields)
                  ELSEIF ( op2.EQ.'RAW ' ) THEN
C                    Read absorber coefficients from unit 8
                     REWIND 8
                     DO l = 1 , 8
                        READ (8,*) (ABC(l,j),j=1,10) ! Absorption coefficients
                        DO j = 1 , 10
                           ABC(l,j) = LOG(ABC(l,j))
                        ENDDO
                     ENDDO
                     DO l = 1 , nfd
                        READ (8,*) (THICK(l,j),j=1,7) ! thickness of absorbers
                     ENDDO
                     DO l = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
                        DO j = 1 , 200
                           ICLUST(l,j) = 0
                        ENDDO
                        DO j = 1 , 20
                           LASTCL(l,j) = 0
                        ENDDO
                        IRAWEX(l) = 0
                     ENDDO

C                    Read input from standard input
                     DO l = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
                        READ (JZB,*) mexl ! experiment number
                        IF ( mexl.EQ.0 ) GOTO 100 ! Back to input loop
                        IRAWEX(mexl) = 1
                        n = NANG(mexl)
                        DO j = 1 , n
                           jj = ITMA(mexl,j) ! Get identity of detector
                           READ (JZB,*) (AKAVKA(k,jj),k=1,8) ! efficiency curve parameters
                           AKAVKA(9,jj) = ideff(mexl)
                        ENDDO
                        READ (JZB,*) kclust ! number of clusters
                        IF ( kclust.NE.0 ) THEN
                           DO j = 1 , kclust
                              READ (JZB,*) numcl ! Number of detectors for this cluster
                              READ (JZB,*) (liscl(k),k=1,numcl) ! Indices of logical detectors
                              LASTCL(l,j) = liscl(numcl) ! Index of last detector in cluster
                              DO k = 1 , numcl
                                 kk = liscl(k)
                                 ICLUST(l,kk) = j ! Set cluster number
                              ENDDO
                           ENDDO
                        ENDIF
                     ENDDO
                     GOTO 100 ! End of OP,RAW - back to input loop

C                 Treat OP,MAP
                  ELSEIF ( op2.EQ.'MAP ' ) THEN
                     GOTO 1200 ! End of OP,MAP 
                  ENDIF ! IF ( op2.EQ.'SIXJ' )
               ENDIF
            ENDIF
         ENDIF
      ENDIF ! End of if (op1.eq."OP, ") if statement

      WRITE (22,99022) op1 , op2
99022 FORMAT (5X,'UNRECOGNIZED OPTION',1X,1A3,1A4)
      GOTO 2000 ! Normal end of execution

C     Treat suboptions of OP,COUL and OP,GOSI
 200  READ (JZB,99023) op1 ! Read the suboption
99023 FORMAT (1A4)
      IF ( op1.EQ.'    ' ) GOTO 100 ! Back to input loop

C     Treat suboption LEVE (levels)
      IF ( op1.EQ.'LEVE' ) THEN
         NMAX = 0
         IF ( ABS(IPRM(1)).EQ.1 ) WRITE (22,99024)
99024    FORMAT (1X/40X,'LEVELS',//5X,'INDEX',5X,'PARITY',9X,'SPIN',11X,
     &           'ENERGY(MEV)')
         ndima = NDIM + 1
         DO k = 1 , ndima
            READ (JZB,*) ipo1 , ipo2 , po2 , po1 ! level number, parity, spin, energy
            IF ( ipo1.EQ.0 ) GOTO 200
            IF ( ipo1.EQ.1 .AND. ABS(po2).LT.1.E-6 ) ISO = 0
            NMAX = NMAX + 1
            SPIN(ipo1) = po2
            IF ( k.EQ.1 ) iph = ipo2
            iprc = ipo2 - iph
            IF ( iprc.NE.0 ) iprc = 1
            IFAC(ipo1) = (-1)**(iprc-INT(po2-SPIN(1)))
            EN(ipo1) = po1
            prp = '+'
            IF ( ipo2.EQ.-1 ) prp = '-'
            IF ( ABS(IPRM(1)).EQ.1 ) WRITE (22,99025) ipo1 , prp , 
     &           SPIN(ipo1) , EN(ipo1)
99025       FORMAT (5X,1I3,11X,1A1,10X,1F4.1,8X,1F10.4)
         ENDDO

C     Treat suboption ME (matrix elements)
      ELSEIF ( op1.EQ.'ME  ' ) THEN
         DO k = 1 , nmemx
            IF ( op2.EQ.'GOSI' ) THEN
               READ (JZB,*) ipo1 , ipo2 , po1 , bl , bu ! lamda, 0, 0, 0, 0 OR ind1, ind2, me, lo, hi
               iopri = 2
               icg = 2
            ELSE
               iopri = 1
               READ (JZB,*) ipo1 , ipo2 , po1 ! lambda, 0, 0 OR ind1, ind2, me
            ENDIF
            IF ( ipo1.NE.0 ) THEN
               IF ( ipo2.EQ.0 ) THEN
                  IF ( ipo1.LE.la ) GOTO 1600 ! Error - wrong sequence of multipolarities
                  LAMMAX = LAMMAX + 1
                  LAMDA(LAMMAX) = ipo1
                  ipo3 = 0
                  IF ( indx.EQ.0 ) GOTO 220
               ELSE
                  MULTI(la) = MULTI(la) + 1
                  indx = indx + 1
                  IF ( ipo1.GT.ABS(ipo2) ) GOTO 1500 ! Error - M.E. does not belong to the upper triangle
                  IF ( ipo1.NE.ipo3 ) THEN
                     IF ( ipo1.LT.ipo3 ) GOTO 1700 ! Error - repeated appearance of the state
                     ipo3 = ipo1
                  ENDIF
                  ELM(indx) = po1
                  mlt(indx) = la
                  LEAD(1,indx) = ipo1
                  LEAD(2,indx) = ABS(ipo2)
                  LDNUM(la,ipo1) = LDNUM(la,ipo1) + 1
                  IF ( op2.EQ.'GOSI' ) THEN
                     IF ( ipo2.LT.0 ) THEN ! If negative, bl and bu are indices
                                           ! to which we fix this element
                        IVAR(indx) = 10000*INT(bl) + INT(bu)
                     ELSE                 ! Otherwise they are limits
                        ELMU(indx) = bu
                        ELML(indx) = bl
                        IF ( ABS(bl-bu).LT.1.E-6 ) THEN
                           IVAR(indx) = 0 ! Fixed
                        ELSE
                           IVAR(indx) = 2
                           IF ( la.GT.4 ) IVAR(indx) = 1
                        ENDIF
                     ENDIF
                     isip = ISEX(ipo1) + 1
                     ISEX(ABS(ipo2)) = MIN(isip,ISEX(ABS(ipo2)))
                  ENDIF
                  GOTO 250
               ENDIF
            ENDIF
            DO kk = 1 , indx
               IF ( ABS(ELM(kk)).LE.1.E-6 ) ELM(kk) = 1.E-6
               IF ( IVAR(kk).GE.10000 ) THEN ! Correlated
                  kk1 = IVAR(kk)/10000
                  kk2 = IVAR(kk) - 10000*kk1
                  la1 = la
                  IF ( kk2.GE.100 ) THEN
                     la1 = kk2/100
                     kk2 = kk2 - 100*la1
                  ENDIF
                  inx1 = MEM(kk1,kk2,la1)
                  ELML(kk) = ELML(inx1)*ELM(kk)/ELM(inx1)
                  ELMU(kk) = ELMU(inx1)*ELM(kk)/ELM(inx1)
                  SA(kk) = ELM(kk)/ELM(inx1)
                  ivari(kk) = IVAR(kk)
                  IVAR(kk) = 1000 + inx1
                  IF ( ELMU(kk).LE.ELML(kk) ) THEN
                     elmi = ELMU(kk)
                     ELMU(kk) = ELML(kk)
                     ELML(kk) = elmi
                  ENDIF
               ENDIF
            ENDDO
            IF ( ipo1.EQ.0 ) GOTO 300
 220        la = ipo1
            IF ( la.GT.LMAXE .AND. la.LE.6 ) LMAXE = la
 250        CONTINUE
         ENDDO
 300     MEMAX = indx
         IF ( la.GT.6 ) MAGEXC = 1 ! Flag that we need magnetic excitations
         memx4 = MULTI(1) + MULTI(2) + MULTI(3) + MULTI(4)
         MEMX6 = memx4 + MULTI(5) + MULTI(6)
         IF ( ABS(IPRM(1)).EQ.1 ) CALL PRELM(iopri)
         DO kh = 1 , NMAX
            IF ( ISEX(kh).EQ.1111 ) ISEX(kh) = 1
         ENDDO
         DO kh = 1 , MEMAX
            ivarh(kh) = IVAR(kh)
         ENDDO

C     Treat suboption CONT (control)
      ELSEIF ( op1.EQ.'CONT' ) THEN
 350     READ (JZB,99026) op1 , fipo1
99026    FORMAT (1A4,1F7.1)
         ipo1 = INT(fipo1)
         IF ( op1.EQ.'ACP,' ) ACCA = 10.**(-fipo1)
         IF ( op1.EQ.'SEL,' ) ITS = 2
         IF ( op1.EQ.'SMR,' ) iosr = 1
         IF ( op1.EQ.'SPL,' ) ISPL = ipo1
         IF ( op1.EQ.'EFF,' ) THEN
            DO jjx = 1 , ipo1
               READ (JZB,*) ipo2 , ijx
               ideff(ipo2) = ijx
            ENDDO
         ENDIF
         IF ( op1.EQ.'FMI,' ) ifm = 1
         IF ( op1.EQ.'TEN,' ) itno = 1
         IF ( op1.EQ.'NCM,' ) NCM = ipo1
         IF ( op1.EQ.'WRN,' ) SGW = fipo1
         IF ( op1.EQ.'INT,' ) THEN
            DO jjx = 1 , ipo1
               READ (JZB,*) ipo2 , ijx
               INTERV(ipo2) = ijx
            ENDDO
         ELSE
            IF ( op1.EQ.'VAC,' ) THEN
               DO jjx = 1 , 7
                  READ (JZB,*) ijx , val
                  IF ( ijx.EQ.0 ) GOTO 350
                  G(ijx) = val
               ENDDO
            ELSE
               IF ( op1.EQ.'DIP,' ) DIPOL = 0.001*fipo1
               IF ( op1.EQ.'ACC,' ) ACCUR = 10.**(-fipo1)
               IF ( op1.EQ.'PRT,' ) THEN
                  DO jjx = 1 , 20
                     READ (JZB,*) inm1 , inm2
                     IF ( inm1.EQ.0 ) GOTO 350
                     IPRM(inm1) = inm2
                  ENDDO
                  GOTO 350
               ELSEIF ( op1.NE.'FIX,' ) THEN
                  IF ( op1.EQ.'SKP,' ) THEN
                     DO jjx = 1 , ipo1
                        READ (JZB,*) ijx
                        JSKIP(ijx) = 0
                     ENDDO
                     GOTO 350
                  ELSE
                     IF ( op1.EQ.'CRF,' ) ICS = 1
                     IF ( op1.EQ.'LCK,' ) THEN
 352                    READ (JZB,*) lck1 , lck2
                        IF ( lck1.EQ.0 ) GOTO 350
                        DO jjx = lck1 , lck2
                           ivarh(jjx) = 0
                           IVAR(jjx) = 0
                        ENDDO
                        GOTO 352
                     ELSE
                        IF ( op1.EQ.'INR,' ) INNR = 1
                        IF ( op1.EQ.'CRD,' ) THEN
                           DO jjx = 1 , ipo1
                              READ (JZB,*) ipo2
                              iecd(ipo2) = 1
                           ENDDO
                           GOTO 350
                        ELSE
                           IF ( op1.EQ.'CCF,' ) IPS1 = ipo1
                           IF ( op1.EQ.'PIN,' ) THEN
                              ipine = ipo1
                              ipinf = 1
                              DO ipp = 1 , ipine
                                 READ (JZB,*) ig1 , ig2
                                 jpin(ig1) = ig2
                              ENDDO
                              GOTO 350
                           ELSE
                              IF ( op1.NE.'END,' ) GOTO 350
                              GOTO 200
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            READ (JZB,*) nallow
            DO jjx = 1 , nallow
               READ (JZB,*) ijk
               IVAR(ijk) = -IVAR(ijk)
            ENDDO
            DO jjx = 1 , MEMAX
               IF ( IVAR(jjx).GE.0 ) THEN
                  IF ( IVAR(jjx).LE.999 ) IVAR(jjx) = 0
               ENDIF
            ENDDO
            DO jjx = 1 , MEMAX
               IF ( IVAR(jjx).LT.0 ) IVAR(jjx) = -IVAR(jjx)
               ivarh(jjx) = IVAR(jjx)
            ENDDO
         ENDIF
         GOTO 350 ! Back to beginning of CONT loop

C     Treat suboption EXPT
      ELSEIF ( op1.EQ.'EXPT' ) THEN
         READ (JZB,*) NEXPT , IZ , XA
         G(1) = 3.             ! AVJI
         G(2) = .02            ! GAMMA
         G(3) = .0345          ! XLAMB
         G(4) = 3.5            ! TIMEC
         G(5) = DBLE(IZ)/XA    ! GFAC
         G(6) = 6.E-06         ! FIEL
         G(7) = .6             ! POWER
         DO k = 1 , NEXPT ! Zn, An, E_p, THETA_lab, M_c, M_A, IAX, phi1, phi2, ikin, ln
            READ (JZB,*) IZ1(k) , XA1(k) , EP(k) , TLBDG(k) , EMMA(k) ,
     &           MAGA(k) , IAXS(k) , fi0 , fi1 , ISKIN(k) , LNORM(k)
            ITTE(k) = 0
            IF ( XA1(k).LT.0. ) ITTE(k) = 1
            XA1(k) = ABS(XA1(k))
            FIEX(k,1) = fi0/57.2957795 ! Convert to radians
            FIEX(k,2) = fi1/57.2957795
            IF ( TLBDG(k).LT.0. ) THEN
               FIEX(k,1) = FIEX(k,1) + 3.14159265
               FIEX(k,2) = FIEX(k,2) + 3.14159265
            ENDIF
         ENDDO

C     Else we don't recognize the suboption
      ELSE
         WRITE (22,99027) op1
99027    FORMAT (5X,'UNRECOGNIZED SUBOPTION',1X,1A4)
         GOTO 2000 ! Normal end of execution
      ENDIF
      GOTO 200 ! Get next suboption

C     Handle OP,ERRO      
 400  IF ( ICS.EQ.1 ) THEN
         REWIND 11
         DO kh1 = 1 , LP4
            READ (11) (CORF(kh1,kh2),kh2=1,LP6) ! LP6 = 32 (maximum number of gamma detectors)
         ENDDO
      ELSE
         CALL FTBM(0,chiss,idr,0,chilo,bten)
         REWIND 11
         DO kh1 = 1 , LP4
            WRITE (11) (CORF(kh1,kh2),kh2=1,LP6) ! LP6 = 32 (maximum number of gamma detectors)
         ENDDO
      ENDIF

      CALL FTBM(3,chiss,idr,1,chilo,bten)
      chis0 = chiss
      WRITE (22,99028) chis0
99028 FORMAT (1X///10X,'***** CENTRAL CHISQ=',1E12.4,1X,'*****'//)
      INHB = 1
      chisl = chiss
      DO kh = 1 , MEMAX
         HLM(kh) = ELM(kh)
      ENDDO
      IF ( idf.EQ.1 ) THEN
         IFBFL = 1
         IF ( irep.NE.2 ) GOTO 700
         IF ( iosr.EQ.0 ) GOTO 700
         REWIND IUNIT3
         READ (IUNIT3,*) ll , mm , kk , inn
         DO inn = 1 , ll
            READ (IUNIT3,*) mm , yyy , zz
         ENDDO
         DO inn = 1 , MEMAX
            READ (IUNIT3,*) mm , ll , kk
         ENDDO
         DO inn = 1 , MEMAX
            READ (IUNIT3,*) mm , yyy
         ENDDO
 450     READ (IUNIT3,*) mm , ll
         IF ( mm.EQ.0 ) THEN
            BACKSPACE IUNIT3
            GOTO 700
         ELSE
            READ (IUNIT3,*) kk , ll , yyy
            READ (IUNIT3,*) (SA(mm),mm=1,MEMAX)
            GOTO 450
         ENDIF
      ELSE
         naxfl = 0
         IF ( ms.EQ.0 ) mend = MEMAX
         IF ( ms.EQ.0 ) ms = 1
         DO kh = ms , mend ! Loop over matrix elements
            DO ij = 1 , 2
               pv = (ELMU(kh)-ELML(kh))/100.
               IF ( ij.NE.1 .OR. (ELM(kh)-ELML(kh)).GE.pv ) THEN
                  IF ( ij.NE.2 .OR. (ELMU(kh)-ELM(kh)).GE.pv ) THEN
                     DO kh1 = 1 , MEMAX
                        SA(kh1) = 0.
                     ENDDO
                     IF ( IVAR(kh).EQ.0 ) GOTO 500
                     SA(kh) = 1.*(-1)**ij
                     kh1 = kh
                     CALL KONTUR(idr,chis0,chisl,ifbp,-1,kh1,sh,bten,
     &                           rem)
                     ELM(kh) = HLM(kh)
                  ENDIF
               ENDIF
            ENDDO
            REWIND 15
            WRITE (15,*) (DEVD(ij),DEVU(ij),ij=1,MEMAX)
 500        CONTINUE
         ENDDO
      ENDIF
 600  IF ( ifbp.EQ.1 ) THEN
         REWIND 17
         DO lkj = 1 , MEMAX
            READ (17,*) ELM(lkj)
         ENDDO
         WRITE (22,99029)
99029    FORMAT (1X///20X,'*** BEST POINT FOUND (TAPE17) ***'///)
         CALL PRELM(3)
      ENDIF
      IF ( naxfl.EQ.0 ) WRITE (22,99051)
      IF ( naxfl.NE.0 ) WRITE (22,99050)
      WRITE (22,99030)
99030 FORMAT (40X,'ESTIMATED ERRORS'//5X,'INDEX',5X,'NI',5X,'NF',5X,
     &        'ME AND ERRORS'//)
      DO kh1 = 1 , MEMAX
         IF ( IVAR(kh1).NE.0 .AND. IVAR(kh1).LE.999 ) THEN
            WRITE (22,99031) kh1 , LEAD(1,kh1) , LEAD(2,kh1) , HLM(kh1)
     &                       , DEVD(kh1) , DEVU(kh1) , DEVD(kh1)
     &                       *100./ABS(HLM(kh1)) , DEVU(kh1)
     &                       *100./ABS(HLM(kh1))
99031       FORMAT (6X,1I3,5X,1I3,4X,1I3,5X,1F9.5,2X,'(',1F9.5,' ,',
     &              1F9.5,')','......',1F7.1,' ,',1F7.1,1X,'PC')
         ENDIF
      ENDDO
      IF ( naxfl.NE.0 ) WRITE (22,99050)
      IF ( naxfl.EQ.0 ) WRITE (22,99051)
      WRITE (22,99032)
99032 FORMAT (40X,'ESTIMATED ERRORS',//5X,'INDEX',5X,'NI',5X,'NF',5X,
     &        'B(E,ML)(OR QUADRUPOLE MOMENT)',' AND ERRORS'//)

      DO kh2 = 1 , MEMAX
         IF ( IVAR(kh2).NE.0 .AND. IVAR(kh2).LE.999 ) THEN
            ispa = LEAD(2,kh2)
            IF ( LEAD(1,kh2).NE.LEAD(2,kh2) ) THEN
               sbe = 2.*SPIN(ispa) + 1.
               be2 = HLM(kh2)*HLM(kh2)/sbe
               be2a = HLM(kh2) + DEVD(kh2)
               be2b = HLM(kh2) + DEVU(kh2)
               be2c = be2b
               IF ( ABS(be2a).GT.ABS(be2b) ) be2b = be2a
               IF ( ABS(be2a-be2c).LT.1.E-6 ) be2a = be2c
               IF ( be2a/HLM(kh2).LE.0. .OR. be2b/HLM(kh2).LE.0. )
     &              be2a = 0.
               be2a = be2a**2/sbe
               be2b = be2b**2/sbe
               WRITE (22,99052) kh2 , LEAD(2,kh2) , LEAD(1,kh2) , be2 , 
     &                          be2a - be2 , be2b - be2
            ELSE
               ispb = INT(SPIN(ispa))*2
               qfac = 3.170662*WTHREJ(ispb,4,ispb,-ispb,0,ispb)
               WRITE (22,99052) kh2 , LEAD(2,kh2) , LEAD(1,kh2) , 
     &                          HLM(kh2)*qfac , DEVD(kh2)*qfac , 
     &                          DEVU(kh2)*qfac
            ENDIF
         ENDIF
      ENDDO
      GOTO 2000 ! Normal end of execution

 700  irea = 0
      IF ( ms.LT.0 ) irea = 1
      IF ( ms.EQ.0 ) mend = MEMAX
      IF ( ms.EQ.0 ) ms = 1
 800  naxfl = 1
      IF ( irea.EQ.1 ) READ (JZB,*) ms , mend
      IF ( ms.NE.0 ) THEN
         DO kh = ms , mend ! For matrix elements
            IF ( ifc.NE.1 ) THEN
               REWIND 18
               DO kh1 = 1 , kh
                  READ (18,*) (KVAR(jyi),jyi=1,MEMAX)
               ENDDO
               DO kh1 = 1 , MEMAX ! For each matrix element
                  ivrh = IVAR(kh1)
                  IF ( KVAR(kh1).EQ.0 ) IVAR(kh1) = 0
                  KVAR(kh1) = ivrh
               ENDDO
            ENDIF
            DO ij = 1 , 2
               sh = DEVU(kh)
               IF ( ij.EQ.1 ) sh = DEVD(kh)
               IF ( ABS(sh).LT.1.E-6 ) sh = (-1)**ij*ABS(HLM(kh))/10.
               ELM(kh) = HLM(kh) + 1.5*sh
               mm = 0
               DO kh1 = 1 , MEMAX ! For each matrix element
                  IF ( ifc.EQ.1 ) KVAR(kh1) = IVAR(kh1)
                  mm = mm + IVAR(kh1)
               ENDDO
               IF ( mm.EQ.0 ) WRITE (22,99033) kh
99033          FORMAT (10X,'ME=',1I3,5X,'NO FREE MATRIX ELEMENTS')
               IF ( mm.NE.0 ) THEN
                  KFERR = 1
                  IF ( iosr.EQ.1 ) WRITE (IUNIT3,*) kh , kh ! For sigma program
                  IF ( iosr.EQ.1 ) WRITE (IUNIT3,*) kh , ij , ELM(kh)
                  LOCKS = 1
                  DLOCK = .05
                  CALL MINI(chiss,-1.D0,2,.0001D0,1000,idr,100000.D0,0,
     &                      iosr,kh,bten)
                  DO kh1 = 1 , MEMAX ! For each matrix element
                     SA(kh1) = (ELM(kh1)-HLM(kh1))/ABS(sh)
                  ENDDO
                  CALL KONTUR(idr,chis0,chisl,ifbp,inpo,kh,sh,bten,rem)
               ENDIF
               DO kh1 = 1 , MEMAX ! For each matrix element
                  IF ( ifc.EQ.1 ) IVAR(kh1) = KVAR(kh1)
                  ELM(kh1) = HLM(kh1)
               ENDDO
            ENDDO
            IF ( ifc.NE.1 ) THEN
               DO kh1 = 1 , MEMAX ! For each matrix element
                  IVAR(kh1) = KVAR(kh1)
               ENDDO
            ENDIF
            REWIND 15
            WRITE (15,*) (DEVD(kh1),DEVU(kh1),kh1=1,MEMAX)
         ENDDO ! Loop on matrix elements kh
         IF ( irea.EQ.1 ) GOTO 800
      ENDIF
      IF ( iosr.NE.0 ) THEN
         im = 0
         WRITE (IUNIT3,*) im , im
      ENDIF
      GOTO 600

 900  jfre = 0
      irfix = 0
      IF ( op2.EQ.'RE,F' ) irfix = 1
 1000 DO jrls = 1 , MEMAX ! For each matrix element
         IF ( IVAR(jrls).NE.0 .OR. irfix.NE.1 ) THEN
            IF ( IVAR(jrls).GT.999 ) THEN
               IF ( jfre.EQ.1 ) GOTO 1100
            ENDIF
            IVAR(jrls) = 2
            ELML(jrls) = -ABS(ELML(jrls))
            ELMU(jrls) = ABS(ELMU(jrls))
            IF ( jrls.GT.MEMX6 ) IVAR(jrls) = 1
         ENDIF
 1100    CONTINUE
      ENDDO ! For each matrix element jrls
      DO jrls = 1 , MEMAX
         ivarh(jrls) = IVAR(jrls)
      ENDDO
      GOTO 100 ! Back to input loop

 1200 CALL CMLAB(0,dsig,ttttt) ! Options MAP, STAR, POINT, MINI etc.
      IF ( ERR ) GOTO 2000 ! Error
      IF ( op2.EQ.'POIN' ) READ (JZB,*) ifwd , slim
      ient = 1
      icg = 1
      IF ( SPIN(1).LT.1.E-6 ) ISO = 0
      IF ( iobl.LT.1 ) THEN
         IF ( op2.NE.'GOSI' ) THEN
            iapx = 0
            DO ii = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
               ILE(ii) = 1
            ENDDO
            nch = 0
            DO jexp = 1 , NEXPT ! For each experiment
               IEXP = jexp
               ttttt = TREP(IEXP)
               dsig = DSIGS(IEXP)
               IF ( op2.NE.'STAR' ) THEN
                  jmm = IEXP
                  IF ( IEXP.NE.1 ) THEN
                     DO lli = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
                        ILE(lli) = ILE(lli) + NYLDE(IEXP-1,lli)
                     ENDDO
                  ENDIF
               ENDIF
               fi0 = FIEX(IEXP,1) ! Lower phi limit
               fi1 = FIEX(IEXP,2) ! Upper phi limit
               CALL LOAD(IEXP,1,icg,0.D0,jj)
               CALL ALLOC(ACCUR)
               CALL SNAKE(IEXP,ZPOL)
               CALL SETIN
               DO j = 1 , LMAX ! For each spin up to ground-state spin + 1
                  polm = DBLE(j-1) - SPIN(1)
                  CALL LOAD(IEXP,2,icg,polm,jj)
                  CALL STING(jj)
                  CALL PATH(jj)
                  CALL INTG(IEXP)
                  CALL TENB(j,bten,LMAX)
                  pr = 0.
                  IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                 WRITE (22,99034) (DBLE(j)-1.-SPIN(1)) , IEXP
99034             FORMAT (1X//40X,'EXCITATION AMPLITUDES'//10X,'M=',
     &                    1F5.1,5X,'EXPERIMENT',1X,1I2//5X,'LEVEL',2X,
     &                    'SPIN',2X,'M',5X,'REAL AMPLITUDE',2X,
     &                    'IMAGINARY AMPLITUDE'//)
                  DO k = 1 , ISMAX ! For substates
                     pr = pr + DBLE(ARM(k,5))**2 + DIMAG(ARM(k,5))**2
                     IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                    WRITE (22,99035) INT(CAT(k,1)) , CAT(k,2) , 
     &                    CAT(k,3) , DBLE(ARM(k,5)) , DIMAG(ARM(k,5))
99035                FORMAT (7X,1I2,3X,1F4.1,2X,1F5.1,2X,1E14.6,2X,
     &                       1E14.6)
                  ENDDO ! Loop on substates k
                  IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                 WRITE (22,99036) pr
99036             FORMAT (1X/5X,'SUM OF PROBABILITIES=',1E14.6)
               ENDDO ! Loop over spins j
               CALL TENS(bten)
               IF ( itno.NE.0 ) THEN ! write statistical tensors on tape 17
                  DO k = 2 , NMAX
                     WRITE (17,*) k
                     DO kk = 1 , 4
                        in1 = (k-1)*28 + 1 + (kk-1)*7
                        in2 = in1 + 2*kk - 2
                        WRITE (17,*) (ZETA(kkk),kkk=in1,in2)
                     ENDDO
                  ENDDO
               ENDIF
               summm = 0.
               DO jgl = 2 , NMAX
                  loct = (jgl-1)*28 + 1
                  summm = summm + ZETA(loct)
               ENDDO
               pop1 = 1. - summm
               jgl = 1
               IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 ) WRITE (22,99053)
     &              jgl , pop1
               DO jgl = 2 , NMAX
                  loct = (jgl-1)*28 + 1
                  IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                 WRITE (22,99053) jgl , ZETA(loct)
               ENDDO
               IF ( op2.NE.'STAR' ) THEN
                  CALL DECAY(ccd,0,ccc)
                  nogeli = NANG(IEXP) ! Number of detector angles for expt
                  jgl1 = 0
                  DO js = 1 , LP2 ! LP2 = 1500 (maximum number of matrix elements)
                     DO jgl = 1 , 20
                        SUMCL(jgl,js) = 0.
                     ENDDO
                  ENDDO
                  DO jgl = 1 , nogeli ! For each detector angle
                     IF ( IRAWEX(IEXP).NE.0 ) THEN
                        IF ( op2.EQ.'POIN' .AND. IPRM(20).EQ.1 )
     &                       WRITE (23,99037) IEXP , jgl , EP(IEXP) , 
     &                       TLBDG(IEXP)
99037                   FORMAT (1x//50x,'CALCULATED YIELDS'//5x,
     &                          'EXPERIMENT ',1I2,2x,'DETECTOR ',1I2/5x,
     &                          'ENERGY ',1F10.3,1x,'MEV',2x,'THETA ',
     &                          1F7.3,1x,'DEG'//5x,'NI',5x,'NF',5x,'II',
     &                          5x,'IF',5x,'E(MeV)',5x,'EFFICIENCY'/)
                     ENDIF
                     gth = AGELI(IEXP,jgl,1)
                     figl = AGELI(IEXP,jgl,2)
                     fm = (fi0+fi1)/2.
                     CALL ANGULA(YGN,idr,1,fi0,fi1,ttttt,gth,figl,jgl,
     &                 op2)
                     IF ( IFMO.NE.0 ) THEN
                        id = ITMA(IEXP,jgl) ! Get identity of detector
                        d = ODL(id) ! Get results of OP,GDET for that detector
                        rx = d*SIN(gth)*COS(figl-fm) - .25*SIN(ttttt)
     &                       *COS(fm)
                        ry = d*SIN(gth)*SIN(figl-fm) - .25*SIN(ttttt)
     &                       *SIN(fm)
                        rz = d*COS(gth) - .25*COS(ttttt)
                        rl = SQRT(rx*rx+ry*ry+rz*rz)
                        thc = TACOS(rz/rl)
                        sf = d*d/rl/rl
                        fic = ATAN2(ry,rx)
                        CALL ANGULA(YGP,idr,1,fi0,fi1,ttttt,thc,fic,jgl,
     &                    op2)
                        DO ixl = 1 , idr
                           ixm = KSEQ(ixl,3)
                           tfac = TAU(ixm)
                           YGN(ixl) = YGN(ixl)
     &                                + .01199182*tfac*BETAR(IEXP)
     &                                *(sf*YGP(ixl)-YGN(ixl))
                        ENDDO
                     ENDIF
                     IF ( IRAWEX(IEXP).NE.0 ) THEN
                        ipd = ITMA(IEXP,jgl) ! Get identity of detector
                        DO jyi = 1 , idr ! For each decay
                           ni = KSEQ(jyi,3)
                           nf = KSEQ(jyi,4)
                           decen = EN(ni) - EN(nf)
                           cocos = SIN(ttttt)*SIN(gth)*COS(fm-figl)
     &                             + COS(ttttt)*COS(gth)
                           decen = decen*(1.+BETAR(IEXP)*cocos)
                           CALL EFFIX(ipd,decen,effi)
                           IF ( op2.EQ.'POIN' .AND. IPRM(20).EQ.1 )
     &                          WRITE (23,99049) ni , nf , SPIN(ni) , 
     &                                 SPIN(nf) , decen , effi
                           YGN(jyi) = YGN(jyi)*effi
                        ENDDO
                        inclus = ICLUST(IEXP,jgl) ! Cluster number for detector jgl
                        IF ( inclus.NE.0 ) THEN
                           DO jyi = 1 , idr ! For each decay
                              SUMCL(inclus,jyi) = SUMCL(inclus,jyi)
     &                           + YGN(jyi)
                           ENDDO
                           IF ( jgl.NE.LASTCL(IEXP,inclus) ) GOTO 1205 ! If it is not the last detector in the cluster
                           DO jyi = 1 , idr ! For each decay
                              YGN(jyi) = SUMCL(inclus,jyi)
                           ENDDO
                        ENDIF
                     ENDIF
                     jgl1 = jgl1 + 1
                     lu = ILE(jgl1)
                     IF ( op2.EQ.'POIN' .OR. IPRM(11).EQ.1 )
     &                    WRITE (22,99048) IEXP , jgl1 , EP(IEXP) , 
     &                    TLBDG(IEXP)
                     jmm = 0
C---- this bit removed in gosia2 start
                     ttttx = TLBDG(IEXP)/57.2957795
                     YGN(IDRN) = YGN(IDRN)*dsig*SIN(ttttx)
                     DO jyi = 1 , idr
                        IF ( jyi.NE.IDRN ) YGN(jyi) = YGN(jyi)
     &                       *dsig*SIN(ttttx)
                     ENDDO
C---- this bit removed in gosia2 end
                     DO jyi = 1 , idr
                        ni = KSEQ(jyi,3)
                        nf = KSEQ(jyi,4)
                        IF ( op2.EQ.'POIN' .OR. IPRM(11).EQ.1 )
     &                       WRITE (22,99049) ni , nf , SPIN(ni) , 
     &                       SPIN(nf) , YGN(jyi) , YGN(jyi)/YGN(IDRN)
                        IF ( ifwd.EQ.1 ) THEN
                           IF ( (YGN(jyi)/YGN(IDRN)).GE.slim ) THEN
                              IF ( jgl1.EQ.1 ) sh1 = YGN(IDRN)
                              jmm = jmm + 1
                              CORF(jmm,1) = DBLE(ni)
                              CORF(jmm,2) = DBLE(nf)
                              CORF(jmm,3) = YGN(jyi)/sh1 ! Not divided by sh1 in gosia2
                              IF ( YGN(jyi).GE.YGN(IDRN) ) CORF(jmm,4)
     &                             = CORF(jmm,3)/20.
                              IF ( YGN(jyi).LT.YGN(IDRN) ) CORF(jmm,4)
     &                             = CORF(jmm,3)
     &                             *(.05+.2*(1.-YGN(jyi)/YGN(IDRN)))
                           ENDIF
                        ENDIF
                        IF ( op2.EQ.'CORR' ) THEN
                           READ (15,*) yydd
                           nch = nch + 1
                           jjjj = IY(lu,jgl1)/1000
                           jyi1 = IY(lu,jgl1) - jjjj*1000
                           IF ( IY(lu,jgl1).EQ.jyi .OR. jjjj.EQ.jyi .OR. 
     &                          jyi1.EQ.jyi ) THEN
                              IF ( IY(lu,jgl1).GE.1000 ) THEN
                                 jyi2 = jyi1 - jjjj
                                 IF ( jyi2.LE.0 ) GOTO 1202
                                 DO ihuj = 1 , jyi2
                                    READ (15,*) yyd1
                                 ENDDO
                                 yydd = yydd + yyd1
                                 YGN(jyi) = YGN(jyi) + YGN(jyi1)
                                 REWIND 15
                                 DO ihuj = 1 , nch
                                    READ (15,*) yyd1
                                 ENDDO
                              ENDIF
                              IF ( IEXP.EQ.1 .AND. lu.EQ.NYLDE(1,1)
     &                             .AND. jgl1.EQ.1 )
     &                             cnst = yydd/YGN(jyi)
                              CORF(lu,jgl1) = YEXP(jgl1,lu)
                              YEXP(jgl1,lu) = YEXP(jgl1,lu)
     &                           /yydd*YGN(jyi)
                              DYEX(jgl1,lu) = DYEX(jgl1,lu)
     &                           /yydd*YGN(jyi)
                              lu = lu + 1
                           ENDIF
                        ENDIF
 1202                   CONTINUE
                     ENDDO
                     IF ( ifwd.EQ.1 ) THEN
                        xw = 1.
                        WRITE (4,*) IEXP , jgl1 , ABS(IZ1(IEXP)) , 
     &                              ABS(XA1(IEXP)) , ABS(EP(IEXP)) , 
     &                              jmm , xw
                        DO jyi = 1 , jmm
                           WRITE (4,*) INT(CORF(jyi,1)) , 
     &                                 INT(CORF(jyi,2)) , CORF(jyi,3) , 
     &                                 CORF(jyi,4)
                        ENDDO
                     ENDIF
 1205                CONTINUE
                  ENDDO ! Loop on detector angles jgl
                  IF ( op2.EQ.'CORR' ) THEN
                     jgl1 = 0
                     DO jgl = 1 , nogeli ! For each detector
                        IF ( IRAWEX(jexp).NE.0 ) THEN
                           inclus = ICLUST(jexp,jgl) ! Cluster number for detector jgl
                           IF ( inclus.NE.0 ) THEN
                              IF ( jgl.NE.LASTCL(jexp,inclus) ) ! If detector is not the last in the cluster
     &                             GOTO 1206
                           ENDIF
                        ENDIF
                        jgl1 = jgl1 + 1
                        READ (3,*) ne , na , zp , ap , xep , nval , waga
                        WRITE (4,*) ne , na , zp , ap , EP(IEXP) , 
     &                              nval , waga
                        WRITE (22,99038) IEXP , jgl1
99038                   FORMAT (///10X,'EXPERIMENT',1X,I2,8X,'DETECTOR',
     &                          1X,I2,//9X,'NI',5X,'NF',5X,'YEXP',8X,
     &                          'YCOR',8X,'COR.F'/)
                        ile1 = ILE(jgl1)
                        DO itp = 1 , nval
                           READ (3,*) ns1 , ns2 , fiex1(1,1,1) , 
     &                                fiex1(1,1,2)
                           ltrn = IY(ile1+itp-1,jgl1)
                           IF ( ltrn.LT.1000 ) THEN
                              ns1 = KSEQ(ltrn,3)
                              ns2 = KSEQ(ltrn,4)
                           ELSE
                              ltrn1 = ltrn/1000
                              ns1 = KSEQ(ltrn1,3)*100
                              ns2 = KSEQ(ltrn1,4)*100
                              ltrn2 = ltrn - ltrn1*1000
                              ns1 = ns1 + KSEQ(ltrn2,3)
                              ns2 = ns2 + KSEQ(ltrn2,4)
                           ENDIF
                           ycorr = YEXP(jgl1,ile1+itp-1)*cnst ! Not multiplied by cnst in gosia2
                           WRITE (4,*) ns1 , ns2 , ycorr , 
     &                                 DYEX(jgl1,ile1+itp-1)*cnst ! Not multiplied by cnst in gosia2
                           WRITE (22,99039) ns1 , ns2 , 
     &                            CORF(ile1+itp-1,jgl1) , ycorr , 
     &                            ycorr/CORF(ile1+itp-1,jgl1)
99039                      FORMAT (5X,I4,5X,I4,3X,E8.3,4X,E8.3,4X,E8.3)
                        ENDDO ! Loop over itp
 1206                   CONTINUE
                     ENDDO ! Loop over jgl
                  ENDIF ! if ( op2.EQ. 'CORR')
               ENDIF
            ENDDO ! Loop over jexp
            IF ( op2.EQ.'STAR' ) oph = op2
            IF ( op2.NE.'STAR' ) THEN
               IF ( op2.EQ.'CORR' ) THEN
                  ntap = 4
                  CALL READY(idr,ntap,ipri)
                  REWIND ntap
               ENDIF
            ENDIF
            GOTO 100 ! Back to input loop
         ENDIF ! if (op2 .NE. 'GOSI') if statement
      ENDIF ! if ( iobl.LT.1 ) if statement

 1300 IF ( iobl.GE.1 ) THEN ! OP,ERRO
         ient = 1
         icg = 2
         nmaxh = NMAX
         lmax1 = LMAX
         sh1 = SPIN(1) ! Save ground-state spin
         sh2 = SPIN(2) ! Save spin of first excited state
         ih1 = IFAC(1)
         ih2 = IFAC(2)
         magh = MAGEXC
         lmaxh = LMAXE
         isoh = ISO
         ISO = 0
         eh1 = ELM(1)
         lh1 = LEAD(1,1)
         lh2 = LEAD(2,1)
         lamh = LAMMAX
         memh = MEMAX
         DO kh = 1 , 8 ! For each multipolarity
            ihlm(kh) = MULTI(kh)
            ihlm(kh+24) = LDNUM(kh,2)
            ihlm(kh+8) = LAMDA(kh)
            ihlm(kh+16) = LDNUM(kh,1)
         ENDDO
         DO jexp = 1 , NEXPT ! For each experiment
            IEXP = jexp
            intvh = INTERV(IEXP)
            DO jgs = 1 , MEMAX
               DO jgr = 1 , 7
                  QAPR(jgs,1,jgr) = 0.
               ENDDO
            ENDDO
            DO iuy = 1 , 6
               XIR(iuy,IEXP) = 0.
            ENDDO
            emhl1 = EMMA(IEXP)
            EMMA(IEXP) = DBLE(MAGA(IEXP))
            jde = 2
            IF ( MAGA(IEXP).EQ.0 ) jde = 1
            DO iuy = 1 , 6
               zmir(iuy,1,IEXP) = 0.
               zmir(iuy,2,IEXP) = 0.
            ENDDO
            CALL LOAD(IEXP,1,2,0.D0,jj)
            DO jgs = 1 , LMAX ! For each spin up to ground-state spin + 1
               polm = DBLE(jgs-1) - SPIN(1)
               CALL LOAD(IEXP,3,2,polm,jj)
               CALL PATH(jj)
               CALL LOAD(IEXP,2,2,polm,jj)
               ictl = 1
               DO kk = 1 , 6
                  ll = ihlm(kk)
                  IF ( ll.NE.0 ) THEN
                     lfini = ll + ictl - 1
                     ict = ictl
                     DO lll = ict , lfini
                        ictl = ictl + 1
                        IF ( jgs.EQ.1 ) XIR(kk,IEXP)
     &                       = MAX(XIR(kk,IEXP),ABS(XI(lll)))
                        r1 = ABS(QAPR(lll,1,1))
                        r2 = ABS(QAPR(lll,1,4))
                        r3 = ABS(QAPR(lll,1,7))
                        rm = MAX(r1,r2,r3)
                        bmx = MAX(ABS(ELMU(lll)),ABS(ELML(lll)))
                        zmir(kk,2,IEXP)
     &                     = MAX(zmir(kk,2,IEXP),rm*bmx/ABS(ELM(lll)),
     &                     rm)
                        r1 = ABS(QAPR(lll,1,2))
                        r2 = ABS(QAPR(lll,1,3))
                        r3 = ABS(QAPR(lll,1,5))
                        r4 = ABS(QAPR(lll,1,6))
                        rm = MAX(r1,r2,r3,r4)
                        zmir(kk,1,IEXP)
     &                     = MAX(zmir(kk,1,IEXP),rm*bmx/ABS(ELM(lll)),
     &                     rm)
                     ENDDO
                     IF ( zmir(kk,1,IEXP).LT..5 ) zmir(kk,1,IEXP) = .5
                     IF ( zmir(kk,2,IEXP).LT..5 ) zmir(kk,2,IEXP) = .5
                  ENDIF
               ENDDO
            ENDDO
            DO kk = 1 , 6
               XIR(kk,IEXP) = XIR(kk,IEXP)*1.01
               DO kh = 1 , 8
                  MULTI(kh) = 0
                  LAMDA(kh) = 0
                  LDNUM(kh,2) = 0
                  LDNUM(kh,1) = 0
               ENDDO
               NMAX = 2
               ELM(1) = 1.
               LEAD(1,1) = 1
               LEAD(2,1) = 2
               SPIN(1) = 0.
               IFAC(1) = 1
               LAMMAX = 1
               MEMAX = 1
               MAGEXC = 0
               kkk = 0
               icg = 1
               IF ( ihlm(kk).NE.0 ) THEN
                  MULTI(kk) = 1
                  LAMDA(1) = kk
                  SPIN(2) = DBLE(kk)
                  IFAC(2) = 1
                  LDNUM(kk,1) = 1
                  icg = 1
                  CALL LOAD(IEXP,1,icg,0.D0,jj)
                  CALL LOAD(IEXP,2,icg,0.D0,jj)
                  CALL PATH(1)
                  sz1 = MIN(zmir(kk,1,IEXP),10.D0)
                  sz2 = zmir(kk,2,IEXP)/50.
                  acof = 2.4009604E-3/zmir(kk,2,IEXP)
                  bcof = 8.163265E-4
                  DO jd = 1 , jde
                     nksi = 5
                     IF ( jd.EQ.2 ) nksi = 10
                     IF ( MAGA(IEXP).EQ.0 ) nksi = 10
                     DO jk = 1 , 3
                        ZETA(jk) = 0.
                     ENDDO
                     nz = 50
                     IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) nz = 1
                     DO jk = 1 , nksi
                        XI(1) = XIR(kk,IEXP)*(jk-1)/(nksi-1)
                        IF ( jk.EQ.1 ) XI(1) = .02
                        s11 = 0.
                        s21 = 0.
                        s12 = 0.
                        s22 = 0.
                        ph1 = 0.
                        ph2 = 0.
                        DO jz = 1 , nz
                           ZETA(jd) = sz2*jz
                           IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) ZETA(jd)
     &                          = sz1
                           IF ( ZETA(jd).LT..1 ) INTERV(IEXP) = 1000
                           IF ( ZETA(jd).GE..1 ) INTERV(IEXP) = intvh
                           CALL ALLOC(ACCUR)
                           CALL SNAKE(IEXP,ZPOL)
                           CALL SETIN
                           CALL STING(1)
                           IF ( kk.GT.2 ) THEN
                              ARM(1,5) = (.9999999,0.)
                              ARM(2,5) = (1.2E-6,0.)
                              ARM(1,6) = (.9999998,0.)
                              ARM(2,6) = (.9E-6,0.)
                              DO kh = 1 , 4
                                 ARM(1,kh) = (-1.E-6,0.)
                                 ARM(2,kh) = (1.E-6,0.)
                              ENDDO
                           ENDIF
                           CALL INTG(IEXP)
                           jp = 2
                           IF ( MAGA(IEXP).NE.0 .AND. jd.EQ.2 ) jp = 3
                           p = DBLE(ARM(1,5))
                           r = DIMAG(ARM(1,5))
                           qr = DBLE(ARM(jp,5))
                           s = DIMAG(ARM(jp,5))
                           test = p*p + r*r + qr*qr + s*s
                           p = p/SQRT(test)
                           s = ABS(r/s)
                           IF ( jk.EQ.1 ) THEN
                              IF ( MAGA(IEXP).EQ.0 ) THEN
                                 q1 = 0.
                                 GOTO 1302
                              ELSEIF ( jd.EQ.2 .OR. MAGA(IEXP).EQ.0 )
     &                                 THEN
                                 q1 = 0.
                                 GOTO 1302
                              ENDIF
                           ENDIF
                           q1 = ARCTG(s,ph1,pi)
                           ph1 = q1
 1302                      IF ( jk.EQ.1 ) THEN
                              IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) THEN
                                 q2 = 0.
                                 GOTO 1304
                              ENDIF
                           ENDIF
                           q2 = ARCCOS(p,ph2,pi)
                           ph2 = q2
 1304                      q1 = q1/ZETA(jd)/2.
                           q2 = q2/ZETA(jd)
                           IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) q2 = -q2
                           IF ( jd.NE.1 .OR. MAGA(IEXP).EQ.0 ) THEN
                              s11 = s11 + q1
                              s12 = s12 + q1*jz
                              s21 = s21 + q2
                              s22 = s22 + jz*q2
                           ENDIF
                        ENDDO
                        IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) THEN
                           PARX(IEXP,2*kk-1,jk) = q1
                           PARX(IEXP,2*kk,jk) = q2
                        ELSE
                           PARXM(IEXP,1,jk,kk) = acof*(2.*s12-51.*s11)
                           PARXM(IEXP,2,jk,kk) = bcof*(101.*s11-3.*s12)
                           PARXM(IEXP,3,jk,kk) = acof*(2.*s22-51.*s21)
                           PARXM(IEXP,4,jk,kk) = bcof*(101.*s21-3.*s22)
                        ENDIF
                     ENDDO ! Loop over jk
                  ENDDO ! Loop over jd
               ENDIF
            ENDDO ! Loop over kk
            EMMA(IEXP) = emhl1
            NMAX = nmaxh
            SPIN(1) = sh1 ! Restore ground-state spin
            SPIN(2) = sh2 ! Restore spin of first excited state
            IFAC(1) = ih1
            IFAC(2) = ih2
            MAGEXC = magh
            ISO = isoh
            ELM(1) = eh1
            LEAD(1,1) = lh1
            LEAD(2,1) = lh2
            LAMMAX = lamh
            MEMAX = memh
            DO kh = 1 , 8 ! For each multipolarity
               LDNUM(kh,2) = ihlm(kh+24)
               MULTI(kh) = ihlm(kh)
               LAMDA(kh) = ihlm(kh+8)
               LDNUM(kh,1) = ihlm(kh+16)
            ENDDO
            INTERV(IEXP) = intvh
         ENDDO ! Loop over experiments jexp

         irix = 7
         REWIND irix
         DO iuy = 1 , 6
            WRITE (irix,*) (XIR(iuy,jj),jj=1,NEXPT)
            WRITE (irix,*) (zmir(iuy,1,jj),zmir(iuy,2,jj),jj=1,NEXPT)
         ENDDO
         DO jj = 1 , NEXPT ! For each experiment
            DO jk = 1 , 4
               DO kuku = 1 , 6
                  WRITE (irix,*) (PARXM(jj,jk,jl,kuku),jl=1,10)
               ENDDO
            ENDDO
            DO jk = 1 , 12
               WRITE (irix,*) (PARX(jj,jk,jl),jl=1,5)
            ENDDO
         ENDDO
         DO jj = 1 , 2
            DO jj1 = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
               IDIVE(jj1,jj) = 1
            ENDDO
         ENDDO
      ELSE ! iobl .lt. 1
         irix = 7
         REWIND irix
         DO iuy = 1 , 6
            READ (irix,*) (XIR(iuy,jj),jj=1,NEXPT)
            READ (irix,*) (zmir(iuy,1,jj),zmir(iuy,2,jj),jj=1,NEXPT)
         ENDDO
         DO jj = 1 , NEXPT ! For each experiment
            DO jk = 1 , 4
               DO kuku = 1 , 6
                  READ (irix,*) (PARXM(jj,jk,jl,kuku),jl=1,10)
               ENDDO
            ENDDO
            DO jk = 1 , 12
               READ (irix,*) (PARX(jj,jk,jl),jl=1,5)
            ENDDO
         ENDDO
         DO jgs = 1 , MEMAX ! For each matrix element
            DO jgr = 1 , 7
               QAPR(jgs,1,jgr) = 0.
            ENDDO
         ENDDO
      ENDIF

C     Handle map
      IF ( IPRM(12).NE.0 ) THEN ! gosia2 has additional .OR. op2 .eq 'MAP '
         IPRM(12) = 0
         DO jex = 1 , NEXPT
            DO lex = 1 , 6
               IF ( MULTI(lex).NE.0 ) THEN
                  WRITE (22,99040) jex , XIR(lex,jex)
99040             FORMAT (1X//30X,'EXPERIMENT',1X,1I2,10X,'MAX.XI=',
     &                    1F6.4)
                  WRITE (22,99041) lex , zmir(lex,2,jex)
99041             FORMAT (1X/30X,'E',1I1,8X,'MI=0',5X,'MAX.ZETA=',
     &                    1F6.3//)
                  WRITE (22,99054)
                  DO kex = 1 , 10
                     xxi = XIR(lex,jex)*(kex-1)/9.
                     WRITE (22,99055) xxi , 
     &                                (PARXM(jex,ilx,kex,lex),ilx=1,4)
                  ENDDO
                  IF ( MAGA(jex).NE.0 ) THEN
                     WRITE (22,99042) lex , zmir(lex,1,jex)
99042                FORMAT (1X//30X,'E',1I1,8X,'MI=+/-1',5X,
     &                       'MAX.ZETA=',1F6.3//)
                     WRITE (22,99054)
                     DO kex = 1 , 5
                        xxi = XIR(lex,jex)*(kex-1)/4.
                        u = 0.
                        WRITE (22,99055) xxi , u , PARX(jex,2*lex-1,kex)
     &                         , u , PARX(jex,2*lex,kex)
                     ENDDO ! Loop on kex
                  ENDIF ! if maga(jex).ne.0
               ENDIF ! if multi(lex).ne.0
            ENDDO ! Loop on lex
         ENDDO ! Loop on jex
      ENDIF ! IPRM(12).ne.0
      IF ( op2.NE.'GOSI' .AND. op2.NE.'ERRO' ) GOTO 100 ! Back to input loop
      IF ( op2.EQ.'ERRO' ) GOTO 400

 1400 DO kh1 = 1 , MEMAX
         HLM(kh1) = ELM(kh1)
      ENDDO
      lfagg = 0
      DO kh1 = 1 , MEMAX
         IVAR(kh1) = ivarh(kh1)
      ENDDO
      CALL MINI(chisq,chiok,nptl,conu,imode,idr,xtest,0,0,0,bten)
      IF ( IPS1.EQ.0 ) GOTO 2000 ! Normal end of execution
      IMIN = IMIN + 1
      DO iva = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         JSKIP(iva) = 1
      ENDDO
      irix = 12
      REWIND irix
      DO lkj = 1 , MEMAX
         WRITE (irix,*) ELM(lkj)
      ENDDO
      IF ( ifm.EQ.1 ) CALL PRELM(3) ! ifm = fast minimisation switch
      IF ( ifm.NE.1 ) GOTO 100 ! Back to input loop
      GOTO 2000 ! Normal end of execution

C.............................................................................
 1500 WRITE (22,99043)
99043 FORMAT (5X,'ERROR-M.E. DOES NOT BELONG TO THE UPPER TRIANGLE')
      GOTO 1900 ! Troubleshoot

C.............................................................................
 1600 WRITE (22,99044)
99044 FORMAT (5X,'ERROR-WRONG SEQUENCE OF MULTIPOLARITIES')
      GOTO 1900 ! Troubleshoot

C.............................................................................
 1700 WRITE (22,99045)
99045 FORMAT (5X,'ERROR-REPEATED APPEARANCE OF THE STATE')
      GOTO 1900 ! Troubleshoot

C.............................................................................
 1800 WRITE (22,99046)
99046 FORMAT (1X///10X,'ERROR-INSUFFICIENT SPACE FOR E-THETA INTEGR ',
     &        'ATION')
      GOTO 1900 ! Troubleshoot

C.............................................................................
C     Troubleshooting
 1900 IF ( ITS.NE.0 ) THEN
         iva = 0
         WRITE (18,*) iva , iva , iva , chisq
         IF ( ITS.NE.2 ) THEN
            WRITE (15,*) iva , chisq , chisq , chisq , chisq
            CALL KLOPOT(kmat,rlr) ! Troubleshooting
         ENDIF
      ENDIF

C     End of execution
 2000 WRITE (22,99047)
99047 FORMAT (15X,'********* END OF EXECUTION **********')
      STOP

C     Handle OP,EXIT
 430  IF ( IPRM(18).NE.0 ) CALL PTICC(idr)
      IF ( oph.EQ.'GOSI' ) THEN
         IF ( lfagg.NE.1 ) THEN
            IF ( IMIN.NE.0 ) THEN
               IF ( IPRM(4).EQ.-1 ) IPRM(4) = 111111
               iskok = IPRM(7) + IPRM(8) + IPRM(13) + IPRM(14)
               IF ( iskok.NE.0 .OR. IPRM(4).NE.111111 ) THEN
                  IF ( iskok.NE.0 ) THEN
                     IF ( IPRM(7).EQ.1 ) IPRM(7) = -1
                     IF ( IPRM(8).EQ.1 ) IPRM(8) = -1
                     IF ( IPRM(3).EQ.1 .AND. NBRA.NE.0 ) IPRM(3) = -1
                     IF ( IPRM(13).EQ.1 ) IPRM(13) = -1
                     IF ( IPRM(14).EQ.1 ) IPRM(14) = -1
                  ENDIF
                  CALL MINI(chisq,chiok,+1,conu,2000,idr,xtest,2,0,0,
     &                      bten)
               ENDIF
            ENDIF
            CALL MIXR(iva,1,chisq,chilo)
            IF ( IPRM(15).NE.0 .AND. KFERR.NE.1 .AND. iyr.NE.0 ) THEN
               WRITE (22,99011)
99011          FORMAT (1X//20X,'CALCULATED LIFETIMES'//5X,'LEVEL',5X,
     &                 'LIFETIME(PSEC)',5X,'EXP',8X,'ERROR'/)
               DO iva = 2 , NMAX
                  DO iva1 = 1 , NLIFT
                     IF ( LIFCT(iva1).EQ.iva ) GOTO 122
                  ENDDO
                  WRITE (22,99012) iva , TAU(iva)
99012             FORMAT (6X,1I3,7X,1E10.4)
                  GOTO 124
 122              WRITE (22,99013) iva , TAU(iva) , TIMEL(1,iva1) , 
     &                   TIMEL(2,iva1)
99013             FORMAT (6X,1I3,7X,1E10.4,5X,1E10.4,4X,1E10.4)
 124              IF ( iva.EQ.NMAX ) THEN
                     IF ( NAMX.GE.1 ) THEN
                        WRITE (22,99014)
99014                   FORMAT (5x,//,
     &         'CALCULATED AND EXPERIMENTAL MATRIX ELEMENTS'
     &         ,//)
                        WRITE (22,99015)
99015                   FORMAT (5x,'NI ','NF ',' EXP. ME   ',
     &                     'CURRENT ME','   SIGMA')
                        DO kq = 1 , NAMX
                           ni = IAMY(kq,1)
                           nf = IAMY(kq,2)
                           ind = IAMX(kq)
                           ess = ELM(ind)
                           esd = EAMX(kq,1)
                           dsd = EAMX(kq,2)
                           WRITE (22,99016) ni , nf , esd , ess , 
     &                        (ess-esd)/dsd
99016                      FORMAT (4x,1I3,1x,1I3,1x,1F9.4,1x,1F9.4,1x,
     &                              1F9.4)
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
            IF ( IMIN.NE.0 ) CALL PRELM(3)
         ENDIF
      ENDIF
      GOTO 1900 ! End of OP,EXIT - troubleshoot

99048 FORMAT (1X//50X,'CALCULATED YIELDS'//5X,'EXPERIMENT ',1I2,2X,
     &        'DETECTOR ',1I2/5X,'ENERGY ',1F10.3,1X,'MEV',2X,'THETA ',
     &        1F7.3,1X,'DEG'//5X,'NI',5X,'NF',5X,'II',5X,'IF',5X,
     &        'YIELD',5X,'NORMALIZED YIELD'/)
99049 FORMAT (4X,1I3,4X,1I3,3X,1F4.1,3X,1F4.1,3X,1E11.5,3X,1E11.5)
99050 FORMAT (1X///44X,'OVERALL')
99051 FORMAT (1X///43X,'DIAGONAL')
99052 FORMAT (6X,1I3,5X,1I3,4X,1I3,5X,1F10.5,2X,'(',1F10.5,' ,',1F10.5,
     &        ')')
99053 FORMAT (2X,'LEVEL',1X,1I2,10X,'POPULATION',1X,1E14.6)
99054 FORMAT (5X,'XI',13X,'Q1',22X,'Q2'///13X,'SLOPE',2X,'INTERCEPT',7X,
     &        'SLOPE',5X,'INTERCEPT'//)
99055 FORMAT (2X,1F6.4,3X,1E8.2,2X,1E8.2,6X,1E8.2,2X,1E8.2)
      END
