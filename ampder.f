 
C----------------------------------------------------------------------
 
      SUBROUTINE AMPDER(I57)
      IMPLICIT NONE
      REAL*8 CAT , D2W , ELM , ELMl , ELMu , rsg , SA , ZETa
      INTEGER*4 i1 , I57 , ibg , iend , iflg , indx , ir , is2 , ISG , 
     &          ISG1 , ISMax , ISStar , ISSto , k , KDIv , lam , LAMda , 
     &          LAMmax , LAMr , lax
      INTEGER*4 ld , LDNum , LEAd , LZEta , m , mm , MSTore , MULti , 
     &          n , NDIm , NDIv , nhold , NMAx , NMAx1 , NPT , NSTart , 
     &          NSTop , NSW , nz
      COMPLEX*16 ARM , EXPo
      COMMON /AZ    / ARM(600,7)
      COMMON /COMME / ELM(500) , ELMu(500) , ELMl(500) , SA(500)
      COMMON /CLCOM / LAMda(8) , LEAd(2,500) , LDNum(8,75) , LAMmax , 
     &                MULti(8)
      COMMON /COEX2 / NMAx , NDIm , NMAx1
      COMMON /CAUX  / NPT , NDIv , KDIv , LAMr(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /PINT  / ISStar(76) , ISSto(75) , MSTore(2,75)
      COMMON /ADBXI / EXPo(500)
      COMMON /CCOUP / ZETa(50000) , LZEta(8)
      COMMON /CLCOM8/ CAT(600,3) , ISMax
      COMMON /CEXC0 / NSTart(76) , NSTop(75)
      DO k = 1 , ISMax
         ARM(k,6) = (0.,0.)
         ARM(k,4) = (0.,0.)
      ENDDO
      ISG1 = ISG
      IF ( NPT.EQ.1 ) ISG1 = ABS(ISG1)
      rsg = DBLE(ISG)
      DO i1 = 1 , LAMmax
         lam = LAMda(i1)
         lax = lam
         nz = LZEta(lam)
         IF ( LAMr(lam).NE.0 ) THEN
            iflg = 1
            nhold = 1
 20         CALL NEWLV(nhold,ld,lam)
            IF ( ld.EQ.0 ) THEN
 30            nhold = nhold + 1
               IF ( NSTart(nhold).EQ.0 ) GOTO 30
               GOTO 20
            ELSE
               ir = NSTart(nhold) - 1
 40            ir = ir + 1
               IF ( ir.LE.ISMax ) THEN
                  n = CAT(ir,1)
                  IF ( n.NE.nhold ) THEN
                     DO mm = 1 , ld
                        m = MSTore(1,mm)
                        IF ( m.NE.nhold ) THEN
                           indx = MSTore(2,mm)
                           ibg = ISStar(mm)
                           iend = ISSto(mm)
                           DO is2 = ibg , iend
                              ARM(is2,4) = ARM(is2,4) + ARM(is2,6)
     &                           *ELM(indx)/EXPo(indx)
                              ARM(is2,6) = (0.,0.)
                           ENDDO
                        ENDIF
                     ENDDO
 42                  CALL NEWLV(n,ld,lam)
                     IF ( ld.EQ.0 ) THEN
                        ir = ir + NSTop(n) - NSTart(n) + 1
                        n = n + 1
                        IF ( n.GT.NMAx ) GOTO 100
                        GOTO 42
                     ELSE
                        nhold = n
                     ENDIF
                  ENDIF
                  CALL LAISUM(ir,n,rsg,lax,ld,nz,I57)
                  GOTO 40
               ENDIF
            ENDIF
         ENDIF
 100  ENDDO
      END
