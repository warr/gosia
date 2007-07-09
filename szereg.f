 
C----------------------------------------------------------------------
 
      SUBROUTINE SZEREG(Lst,Ls,L)
      IMPLICIT NONE
      REAL*8 CORf , DYEx , dyh , UPL , YEXp , yh , YNRm
      INTEGER*4 ia , ib , IDRn , ih , ILE , inx , IY , k , L , Ls , 
     &          lsp , Lst , lst1 , NYLde
      COMMON /YEXPT / YEXp(32,1500) , IY(1500,32) , CORf(1500,32) , 
     &                DYEx(32,1500) , NYLde(50,32) , UPL(32,50) , 
     &                YNRm(32,50) , IDRn , ILE(32)
      IF ( Lst.EQ.Ls ) RETURN
      lst1 = Lst
      lsp = Ls - 1
 100  ia = IY(lst1,L)
      IF ( ia.GT.1000 ) ia = ia/1000
      inx = lst1
      DO k = lst1 , lsp
         ib = IY(k+1,L)
         IF ( ib.GT.1000 ) ib = ib/1000
         ia = MIN(ia,ib)
         IF ( ia.EQ.ib ) inx = k + 1
      ENDDO
      IF ( inx.NE.lst1 ) THEN
         ih = IY(lst1,L)
         IY(lst1,L) = IY(inx,L)
         IY(inx,L) = ih
         yh = YEXp(L,lst1)
         dyh = DYEx(L,lst1)
         YEXp(L,lst1) = YEXp(L,inx)
         DYEx(L,lst1) = DYEx(L,inx)
         YEXp(L,inx) = yh
         DYEx(L,inx) = dyh
      ENDIF
      lst1 = lst1 + 1
      IF ( lst1.GT.lsp ) RETURN
      GOTO 100
      END
