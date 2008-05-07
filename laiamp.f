 
C----------------------------------------------------------------------
C SUBROUTINE LAIAMP
C
C Called by: STING
C Calls:     FAZA1, LEADF, STAMP, TCABS
C
C Purpose: calculate excitation amplitudes
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      ELM    - matrix elements
C      EPS    - epsilon
C      EROOT  - sqrt(epsilon^2 - 1)
C      IEXP   - number of experiment
C      ISG    - index for sigma
C      LAMDA  - list of multipolarities to calculate
C      LAMMAX - number of multipolarities to calculate
C      LAMR   - flag = 1 if we should calculate this multipolarity
C      LDNUM  - number of matrix elements with each multipolarity populating levels
C      LZETA  - index in ZETA to coupling coefficients for given multipolarity
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C      XI     - xi coupling coefficients
C      ZETA   - various coefficients
C
C Formal parameters:
C      Ir     - index of substate
C      W0     - omega limit
      
      SUBROUTINE LAIAMP(Ir,W0)
      IMPLICIT NONE
      REAL*8 CAT , D2W , ELM , ELML , ELMU , EPS , epsi , EROOT , errt , 
     &       FIEX , pm , ppp , rmir , rmis , rmu , SA , TCABS , W0 , 
     &       XI , xiv
      REAL*8 z , ZETA
      INTEGER*4 i1 , i2 , i3 , IAXS , IEXP , indx , Ir , is , is1 , 
     &          is2 , ISG , ISG1 , ISMAX , ismin , isplus , KDIV , la , 
     &          lam , LAMDA , LAMMAX
      INTEGER*4 LAMR , ld , LDNUM , LEAD , LEADF , LZETA , m , MEM , 
     &          mrange , mua , MULTI , NDIV , NPT , NSTART , NSTOP , 
     &          NSW , nz
      COMPLEX*16 ARM , STAMP , dis , uhuj
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /AZ    / ARM(600,7)
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /CXI   / XI(1500)

      ppp = 0.
      epsi = EPS(IEXP) ! epsilon
      errt = EROOT(IEXP) ! sqrt(epsilon^2 - 1)
      rmir = CAT(Ir,3) ! m quantum number of substate Ir
      
      DO i1 = 1 , LAMMAX ! Loop on multipolarity
         lam = LAMDA(i1) ! Get multipolarity
         nz = LZETA(lam) ! nz is an index into ZETA array for this multipolarity
         IF ( LAMR(lam).NE.0 ) THEN
            la = lam
            IF ( lam.GT.6 ) lam = lam - 6 ! la = 7,8 for M1,M2
            ld = LDNUM(la,1) ! Number of matrix elements with multipolarity la, connecting to ground state
            IF ( ld.NE.0 ) THEN
               DO i2 = 1 , ld ! Loop on matrix elements of that multipolarity connected to ground state
                  m = LEADF(1,i2,la) ! m is level index connected to ground state by element i2, mul. la
                  indx = MEM(1,m,la)
                  xiv = XI(indx) ! xi value
                  ismin = 0
                  is1 = NSTART(m) ! Index of first substate for level m
                  IF ( NSTART(m).NE.0 ) THEN
                     isplus = INT(rmir-CAT(is1,3)) - lam
                     IF ( isplus.LT.0 ) THEN
                        ismin = isplus
                        isplus = 0
                     ENDIF
                     is2 = is1 + isplus - 1
                     mrange = 2*lam + 1 + ismin
                     IF ( is2+mrange.GT.NSTOP(m) ) mrange = NSTOP(m)
     &                    - is2
                     IF ( mrange.GT.0 ) THEN ! If there are substates for level m
                        DO i3 = 1 , mrange
                           is = is2 + i3
                           nz = nz + 1
                           z = ZETA(nz) ! zeta coefficient
                           rmis = CAT(is,3) ! m quantum number of substate is
                           rmu = rmis - rmir
                           mua = ABS(rmu) + 1.1 ! delta-mu + 1

C                          Only consider electromagnetic and delta-mu = 0 magnetic
C                          contribution
                           IF ( lam.LE.6 .OR. mua.NE.1 ) THEN
C                             calculate complex phase (dis)
                              CALL FAZA1(la,mua,rmir,rmis,dis,rmu)
                              pm = ELM(indx)*z ! Matrix element * zeta
C                             estimate amplitude
                              uhuj = STAMP(epsi,errt,xiv,.03D0,W0,lam,
     &                               mua)
                              ARM(is,5) = dis*pm*uhuj
                              ppp = ppp + TCABS(ARM(is,5))
     &                              *TCABS(ARM(is,5))
                           ENDIF

                        ENDDO ! Loop over substates
                     ENDIF ! If there are substates
                  ENDIF ! If there are substates for level m
               ENDDO ! Loop on matrix elements connected to ground state with multipolarity la
            ENDIF ! If there are matrix elements of this multipolarity connecting to the ground state
         ENDIF
      ENDDO ! Loop on lambda
      ARM(Ir,5) = CMPLX(SQRT(1.-ppp),0.)
      END
