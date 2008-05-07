 
C----------------------------------------------------------------------
C SUBROUTINE PRELM
C
C Called by: GOSIA
C Calls:     KONTUR, MINI
C
C Purpose: print matrix elements
C
C Uses global variables:
C      ELM    - matrix elements
C      ELML   - lower limits on matrix elements
C      ELMU   - upper limits on matrix elements
C      HLM    - matrix elements before minimisation
C      IVAR   - indicates a limit or correlation is set
C      LDNUM  - number of matrix elements with each multipolarity populating levels
C      LEAD   - pair of levels involved in each matrix element
C      MULTI  - number of matrix elements having a given multipolarity
C      NMAX   - number of levels
C      SPIN   - spin of level
C
C Formal parameters:
C      Iop    - print flag (controls what is written to output).
 
      SUBROUTINE PRELM(Iop)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , b , DIPOL , ELM , ELML , ELMU , EN , HLM , 
     &       pv , SA , SPIN , ste , ZPOL
      INTEGER*4 inx , Iop , ISO , isp , IVAR , j , k , kk , l , LAMDA , 
     &          LAMMAX , LDNUM , LEAD , LMAXE , m , MAGEXC , MEMAX , 
     &          MEMX6 , MULTI , NDIM
      INTEGER*4 NMAX , NMAX1
      CHARACTER*3 wrn
      COMMON /HHH   / HLM(1500)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA , 
     &                ISO
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /COEX2 / NMAX , NDIM , NMAX1

      inx = 0
      WRITE (22,99001)
99001 FORMAT (2X/40X,'MATRIX ELEMENTS',//)

      DO j = 1 , 8
         m = MULTI(j)
         IF ( m.NE.0 ) THEN
            WRITE (22,99002) j
99002       FORMAT (5X,'MULTIPOLARITY=',1I1)
            IF ( Iop.EQ.1 ) WRITE (22,99003)
99003       FORMAT (4X,'INDEX',3X,'NF',5X,'NS',10X,'ME')
            IF ( Iop.EQ.2 ) WRITE (22,99004)
99004       FORMAT (4X,'INDEX',3X,'NF',5X,'NS',10X,'ME',15X,'LIMITS')
            IF ( Iop.EQ.3 ) WRITE (22,99005)
99005       FORMAT (4X,'INDEX',3X,'NF',5X,'NS',10X,'ME',10X,'PC CHANGE',
     &              5X,'RED. TRANS. PROB.')

            DO k = 1 , NMAX
               l = LDNUM(j,k)
               IF ( l.NE.0 ) THEN
                  DO kk = 1 , l
                     inx = inx + 1

                     IF ( Iop.EQ.2 ) THEN ! Iop is 2
                        IF ( IVAR(inx).EQ.0 ) THEN ! Fixed
                           WRITE (22,99006) inx , LEAD(1,inx) , 
     &                            LEAD(2,inx) , ELM(inx)
                        ELSEIF ( IVAR(inx).GT.1000 ) THEN ! Correlation
                           WRITE (22,99007) inx , LEAD(1,inx) , 
     &                            LEAD(2,inx) , ELM(inx) , 
     &                            (IVAR(inx)-1000)
                        ELSE ! Limit
                           WRITE (22,99009) inx , LEAD(1,inx) , 
     &                            LEAD(2,inx) , ELM(inx) , ELML(inx) , 
     &                            ELMU(inx)
                        ENDIF

                     ELSEIF ( Iop.EQ.3 ) THEN ! Iop is 3
                        isp = LEAD(2,inx)
                        pv = (ELMU(inx)-ELML(inx))/100.
                        wrn = '   '
                        IF ( (ELM(inx)-ELML(inx)).LT.pv ) wrn = '*?*'
                        IF ( (ELMU(inx)-ELM(inx)).LT.pv ) wrn = '*?*'
                        ste = HLM(inx)
                        b = ELM(inx)*ELM(inx)/(2.*SPIN(isp)+1.)
                        IF ( LEAD(1,inx).EQ.LEAD(2,inx) ) b = 9999999.
                        WRITE (22,99009) inx , LEAD(1,inx) , LEAD(2,inx)
     &                         , ELM(inx) , 100.*(ELM(inx)-ste)/ste , 
     &                         b , wrn
                     ELSE ! Iop is 1
                        WRITE (22,99008) inx , LEAD(1,inx) , LEAD(2,inx)
     &                         , ELM(inx)
                     ENDIF

                  ENDDO ! Loop on kk
               ENDIF ! If l .ne. 0
            ENDDO ! Loop on k
         ENDIF ! If m .ne. 0
      ENDDO ! Loop on j

99006 FORMAT (5X,1I3,5X,1I2,5X,1I2,5X,1F10.5,5X,'FIXED')
99007 FORMAT (5X,1I3,5X,1I2,5X,1I2,5X,1F10.5,5X,'COUPLED TO',1X,1I3)
99008 FORMAT (5X,1I3,5X,1I2,5X,1I2,5X,1F10.5)
99009 FORMAT (5X,1I3,5X,1I2,5X,1I2,3(5X,1F10.5),1A3)
      END
