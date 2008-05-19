 
C----------------------------------------------------------------------
C SUBROUTINE GKVAC
C
C Called from: DECAY
C Calls:       GKK
C
C Purpose: calculate the nuclear deorientation and store the results in the
C GKI array of common GVAC
C
C Uses global variables:
C      BETAR  - recoil beta
C      GKI    - G_k for a single level
C      IEXP   - experiment number
C      ITTE   - thick target experiment flag
C      SPIN   - spin of level
C      TAU    - lifetime in picoseconds
C      VACDP  - G_k for each level
C      XLAMB  - Lambda*       (this is G(3) in GOSIA)
C
C Formal parameters:
C      Il     - level index
 
      SUBROUTINE GKVAC(Il)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , AVJI , beta , DIPOL , EN , EP , FIEL ,
     &       GAMMA , GFAC , GKI , POWER , sp
      REAL*8 SPIN , SUM , TAU , time , TIMEC , TLBDG , VINF , 
     &       XA , XA1 , XLAMB , ZPOL
      INTEGER*4 i , Il , ISO , ITTE , IZ , IZ1 , 
     &          KSEQ , NEXPT
      COMMON /LEV   / TAU(75) , KSEQ(1500,4)
      INCLUDE 'brec.inc'
      COMMON /GGG   / AVJI , GAMMA , XLAMB , TIMEC , GFAC , FIEL , POWER
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /GVAC  / GKI(3) , SUM(3)
      INCLUDE 'kin.inc'
      INCLUDE 'vac.inc'
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA , 
     &                ISO
      COMMON /THTAR / ITTE(50)

      IF ( ABS(XLAMB).GE.1.E-9 ) THEN
         IF ( ITTE(IEXP).EQ.0 ) THEN
            sp = SPIN(Il) ! Spin of level
            beta = BETAR(IEXP)
            time = TAU(Il) ! lifetime of level
            CALL GKK(IZ,beta,sp,time,Il)
            VACDP(1,Il) = GKI(1)
            VACDP(2,Il) = GKI(2)
            VACDP(3,Il) = GKI(3)
            GOTO 99999
         ENDIF
      ENDIF
      DO i = 1 , 3
         VACDP(i,Il) = 1.
      ENDDO
99999 END
