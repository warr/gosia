 
C----------------------------------------------------------------------
C SUBROUTINE SIXEL
C
C Called by: CEGRY
C
C Purpose:
C
C Uses global variables:
C      ARM    - reduced matrix elements
C      DEV    -
C      IEXP   - experiment number
C      ITS    -
C      KVAR   -
C
C Formal parameters:
C      Rik    -
C      Rv     -
C      Em     -
C      Jk     -
C      Kk     -
C      Indx   -
C      Lu     -
 
      SUBROUTINE SIXEL(Rik,Rv,Em,Jk,Kk,Indx,Lu)
      IMPLICIT NONE
      REAL*8 a1 , al , al1 , c1 , c2 , DEV , Em , EPS , EROOT , FIEX , 
     &       Rik , rn , Rv , rx
      INTEGER*4 IAXS , IEXP , Indx , ITS , j , j1 , Jk , Kk , kk6 , 
     &          KVAR , l , l1 , Lu
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(600,7)
      COMMON /ODCH  / DEV(500)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /TRB   / ITS
      COMMON /SEL   / KVAR(500)

      kk6 = Kk + 5
      rn = DEV(Lu)
      al = (Rv-rn)*20./Rik
      IF ( ITS.EQ.1 .AND. KVAR(Indx).NE.0 ) WRITE (18,*) Lu , Indx , 
     &     IEXP , al/Em
      al1 = ABS(al)
      IF ( ITS.EQ.2 ) WRITE (18,*) Lu , Indx , IEXP , al1
      IF ( al1.LE.ABS(IMAG(ARM(kk6,Jk))) ) RETURN

      DO j = Kk , kk6
         a1 = ABS(IMAG(ARM(j,Jk)))
         IF ( al1.GT.a1 ) THEN
            j1 = j + 1
            DO l = j1 , kk6
               l1 = kk6 + j1 - l
               c1 = DBLE(ARM(l1-1,Jk))
               c2 = IMAG(ARM(l1-1,Jk))
               ARM(l1,Jk) = CMPLX(c1,c2)
            ENDDO
            rx = DBLE(Indx)
            ARM(j,Jk) = CMPLX(rx,al)
            GOTO 99999
         ENDIF
      ENDDO

99999 END
