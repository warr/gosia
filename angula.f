 
C----------------------------------------------------------------------
 
      SUBROUTINE ANGULA(Ygn,Idr,Iful,Fi0,Fi1,Trec,Gth,Figl,Ngl)
      IMPLICIT NONE
      REAL*8 AGELI , alab , arg , at , attl , BETAR , bt , CC , DELLA , 
     &       DELTA , EG , ENDEC , ENZ , EPS , EROOT , f , Fi0 , fi01 , 
     &       Fi1 , fi11
      REAL*8 FIEX , Figl , FP , GKP , Gth , Q , qv , sm , TAU , Trec , 
     &       Ygn , ylmr , ZETA
      INTEGER*4 IAXS , Idr , IEXP , ifn , Iful , ig , il , inat , inx1 , 
     &          ipd , is , ITMA , ITTE , iu , ixs , j , ji , jj , jm , k
      INTEGER*4 KLEC , kq , KSEQ , l , lf , lf1 , LZETA , mind , NANG , 
     &          Ngl , NICC , nlv
      DIMENSION f(4) , ylmr(9,9) , at(28) , alab(9,9) , attl(9,9) , 
     &          Ygn(500)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /TRA   / DELTA(500,3) , ENDEC(500) , ITMA(50,200) , 
     &                ENZ(200)
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) , 
     &                NICC , NANG(200)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /LCDL  / DELLA(500,3)
      COMMON /CATLF / FP(4,500,3) , GKP(4,500,2) , KLEC(75)
      COMMON /BREC  / BETAR(50)
      COMMON /THTAR / ITTE(50)
      DO l = 1 , Idr
         nlv = KSEQ(l,3)
         il = (nlv-1)*28
         inx1 = KSEQ(l,2)
         DO j = 1 , 4
            f(j) = FP(j,l,1)*DELLA(l,1)
         ENDDO
         IF ( inx1.NE.0 ) THEN
            DO j = 1 , 4
               f(j) = f(j) + 2.*FP(j,l,3)*DELLA(l,3) + FP(j,l,2)
     &                *DELLA(l,2)
            ENDDO
         ENDIF
         DO j = 1 , 4
            f(j) = f(j)*TAU(nlv)
            iu = (j-1)*7
            ifn = 2*j - 1
            IF ( IAXS(IEXP).EQ.0 ) ifn = 1
            DO kq = 1 , ifn
               is = iu + kq
               ig = is + il
               at(is) = ZETA(ig)*f(j)
            ENDDO
         ENDDO
         IF ( Iful.EQ.1 ) THEN
            DO j = 1 , 9
               DO k = 1 , 9
                  alab(j,k) = 0.
                  attl(j,k) = 0.
               ENDDO
            ENDDO
            DO j = 1 , 4
               lf = 2*j - 1
               lf1 = lf
               IF ( IAXS(IEXP).EQ.0 ) lf1 = 1
               DO k = 1 , lf1
                  inat = (j-1)*7 + k
                  alab(lf,k) = at(inat)
               ENDDO
            ENDDO
            bt = BETAR(IEXP)
            IF ( ITTE(IEXP).NE.1 ) CALL RECOIL(alab,attl,bt,Trec)
            IF ( l.EQ.1 ) CALL YLM1(Gth,ylmr)
            ixs = IAXS(IEXP)
            fi01 = Fi0 - Figl
            fi11 = Fi1 - Figl
            CALL FIINT1(fi01,fi11,alab,ixs)
            Ygn(l) = alab(1,1)*.0795774715
            DO j = 2 , 9
               sm = ylmr(j,1)*alab(j,1)
               IF ( IAXS(IEXP).NE.0 ) THEN
                  DO k = 2 , j
                     sm = sm + 2.*ylmr(j,k)*alab(j,k)
                  ENDDO
               ENDIF
               ipd = ITMA(IEXP,Ngl)
               arg = (ENDEC(l)-ENZ(ipd))**2
               qv = (Q(3,ipd,j-1)*Q(2,ipd,j-1)+Q(1,ipd,j-1)*arg)
     &              /(Q(2,ipd,j-1)+arg)
               Ygn(l) = Ygn(l) + sm*qv
            ENDDO
         ELSE
            ixs = IAXS(IEXP)
            fi01 = Fi0 - Figl
            fi11 = Fi1 - Figl
            CALL FIINT(fi01,fi11,at,ixs)
            IF ( l.EQ.1 ) CALL YLM(Gth,ylmr)
            Ygn(l) = at(1)*.0795774715
            DO jj = 1 , 3
               ji = jj*7 + 1
               sm = ylmr(jj,1)*at(ji)
               IF ( IAXS(IEXP).NE.0 ) THEN
                  mind = 2*jj + 1
                  DO jm = 2 , mind
                     ji = ji + 1
                     sm = ylmr(jj,jm)*at(ji)*2. + sm
                  ENDDO
               ENDIF
               ipd = ITMA(IEXP,Ngl)
               arg = (ENDEC(l)-ENZ(ipd))**2
               qv = (Q(3,ipd,2*jj)*Q(2,ipd,2*jj)+Q(1,ipd,2*jj)*arg)
     &              /(Q(2,ipd,2*jj)+arg)
               Ygn(l) = Ygn(l) + sm*qv
            ENDDO
         ENDIF
      ENDDO
      END
