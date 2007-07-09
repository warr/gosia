 
C----------------------------------------------------------------------
 
      SUBROUTINE GKVAC(Il)
      IMPLICIT NONE
      REAL*8 ACCa , ACCur , AKS , AVJi , beta , BETar , DIPol , DQ , 
     &       EN , EP , EPS , EROot , FIEl , FIEx , GAMma , GFAc , GKI , 
     &       POWer , QCEn , sp
      REAL*8 SPIn , SUM , TAU , time , TIMec , TLBdg , VACdp , VINf , 
     &       XA , XA1 , XLAmb , XNOr , ZPOl
      INTEGER*4 i , IAXs , IBYp , IEXp , Il , ISO , ITTe , IZ , IZ1 , 
     &          KSEq , NEXpt
      COMMON /LEV   / TAU(75) , KSEq(500,4)
      COMMON /BREC  / BETar(50)
      COMMON /GGG   / AVJi , GAMma , XLAmb , TIMec , GFAc , FIEl , POWer
      COMMON /CX    / NEXpt , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBdg(50) , VINf(50)
      COMMON /GVAC  / GKI(3) , SUM(3)
      COMMON /KIN   / EPS(50) , EROot(50) , FIEx(50,2) , IEXp , IAXs(50)
      COMMON /VAC   / VACdp(3,75) , QCEn , DQ , XNOr , AKS(6,75) , IBYp
      COMMON /COEX  / EN(75) , SPIn(75) , ACCur , DIPol , ZPOl , ACCa , 
     &                ISO
      COMMON /THTAR / ITTe(50)
      IF ( ABS(XLAmb).GE.1.E-9 ) THEN
         IF ( ITTe(IEXp).EQ.0 ) THEN
            sp = SPIn(Il)
            beta = BETar(IEXp)
            time = TAU(Il)
            CALL GKK(IZ,beta,sp,time,Il)
            VACdp(1,Il) = GKI(1)
            VACdp(2,Il) = GKI(2)
            VACdp(3,Il) = GKI(3)
            GOTO 99999
         ENDIF
      ENDIF
      DO i = 1 , 3
         VACdp(i,Il) = 1.
      ENDDO
99999 END
