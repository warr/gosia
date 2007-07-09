 
C----------------------------------------------------------------------
 
      SUBROUTINE STING(Irld)
      IMPLICIT NONE
      REAL*8 CAT , D2W , ELM , ELMl , ELMu , rsg , SA , w0 , ZETa
      INTEGER*4 i , i57 , ibg , iend , IFLg , indx , IRA , Irld , is2 , 
     &          ISG , ISG1 , ISMax , ISStar , ISSto , j , j1 , jj , 
     &          KDIv , lam , LAMda
      INTEGER*4 LAMmax , LAMr , ld , LDNum , LEAd , LZEta , maxh , 
     &          MAXla , mm , MSTore , MULti , n , NDIv , NPT , NSW , nz
      COMPLEX*16 ARM , EXPo
      COMMON /CLCOM / LAMda(8) , LEAd(2,500) , LDNum(8,75) , LAMmax , 
     &                MULti(8)
      COMMON /AZ    / ARM(600,7)
      COMMON /COMME / ELM(500) , ELMu(500) , ELMl(500) , SA(500)
      COMMON /ADBXI / EXPo(500)
      COMMON /FLA   / IFLg
      COMMON /PINT  / ISStar(76) , ISSto(75) , MSTore(2,75)
      COMMON /CCOUP / ZETa(50000) , LZEta(8)
      COMMON /CAUX  / NPT , NDIv , KDIv , LAMr(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /CLCOM8/ CAT(600,3) , ISMax
      COMMON /RNG   / IRA(8) , MAXla
      maxh = MAXla
 100  ISG = -1
      n = 1
      rsg = -1.
      IFLg = 1
      w0 = IRA(MAXla)*.03 + .03
      DO j = 1 , ISMax
         DO jj = 1 , 6
            ARM(j,jj) = (0.,0.)
         ENDDO
      ENDDO
      ARM(Irld,5) = (1.,0.)
      DO j = 1 , 8
         LAMr(j) = 0
      ENDDO
      LAMr(MAXla) = 1
      NPT = IRA(MAXla) + 1
      IF ( MAXla.EQ.7 .AND. IRA(2).NE.0 ) THEN
         LAMr(2) = 1
         NPT = NPT - 1
         w0 = w0 - .03
      ENDIF
      NDIv = 0
      KDIv = 0
      DO j = 1 , 4
         NPT = NPT - 1
         DO j1 = 1 , LAMmax
            lam = LAMda(j1)
            IF ( LAMr(lam).NE.0 ) THEN
               CALL NEWLV(n,ld,lam)
               IF ( ld.NE.0 ) THEN
                  nz = LZEta(lam)
                  ld = LDNum(lam,1)
                  i57 = 5
                  CALL LAISUM(Irld,n,rsg,lam,ld,nz,i57)
                  DO mm = 1 , ld
                     indx = MSTore(2,mm)
                     ibg = ISStar(mm)
                     iend = ISSto(mm)
                     DO is2 = ibg , iend
                        ARM(is2,4) = ARM(is2,4) + ARM(is2,6)*ELM(indx)
     &                               /EXPo(indx)
                        ARM(is2,6) = (0.,0.)
                     ENDDO
                  ENDDO
               ELSEIF ( j1.EQ.MAXla ) THEN
                  IRA(MAXla) = -IRA(MAXla)
                  DO jj = 1 , LAMmax
                     lam = LAMda(jj)
                     IF ( IRA(lam).GT.0 ) GOTO 105
                  ENDDO
 105              MAXla = LAMda(jj)
                  GOTO 100
               ENDIF
            ENDIF
         ENDDO
         IF ( j.EQ.4 ) GOTO 200
         DO i = 1 , ISMax
            ARM(i,j) = ARM(i,4)
            ARM(i,4) = (0.,0.)
         ENDDO
      ENDDO
 200  CALL LAIAMP(Irld,w0)
      MAXla = maxh
      DO jj = 1 , 8
         IRA(jj) = ABS(IRA(jj))
      ENDDO
      END
