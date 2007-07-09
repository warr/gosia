 
C----------------------------------------------------------------------
 
      SUBROUTINE INTG(Ien)
      IMPLICIT NONE
      REAL*8 ACC50 , ACCa , ACCur , CAT , D2W , DIPol , EN , f , rim , 
     &       rl , SPIn , srt , ZPOl
      INTEGER*4 i , i57 , Ien , IFAc , IFLg , ihold , intend , INTerv , 
     &          IPAth , ir , ir1 , IRA , ISG , ISG1 , ISMax , ISO , k , 
     &          kast , KDIv , LAMr
      INTEGER*4 MAGa , MAXla , mir , n , NDIm , NDIv , NMAx , NMAx1 , 
     &          NPT , NSTart , NSTop , NSW
      COMPLEX*16 ARM , hold
      COMMON /COEX  / EN(75) , SPIn(75) , ACCur , DIPol , ZPOl , ACCa , 
     &                ISO
      COMMON /AZ    / ARM(600,7)
      COMMON /RNG   / IRA(8) , MAXla
      COMMON /A50   / ACC50
      COMMON /CLCOM0/ IFAc(75)
      COMMON /CLCOM8/ CAT(600,3) , ISMax
      COMMON /COEX2 / NMAx , NDIm , NMAx1
      COMMON /CAUX  / NPT , NDIv , KDIv , LAMr(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /FLA   / IFLg
      COMMON /CEXC0 / NSTart(76) , NSTop(75)
      COMMON /PTH   / IPAth(75) , MAGa(75)
      COMMON /CEXC9 / INTerv(50)
      intend = INTerv(Ien)
      D2W = .03
      NSW = 1
      kast = 0
      NDIv = 0
      KDIv = 0
 100  IF ( (NPT+NSW).GT.IRA(MAXla) .AND. ISG.GT.0 ) RETURN
      DO i = 1 , 8
         LAMr(i) = 0
         IF ( (NPT+NSW).LT.IRA(i) ) LAMr(i) = 1
      ENDDO
      IF ( ISO.EQ.0 ) THEN
         DO n = 1 , NMAx
            ir = NSTart(n) - 1
 120        ir = ir + 1
            ARM(ir,7) = ARM(ir,5)
     &                  + D2W/24.*(55.0*ARM(ir,4)-59.0*ARM(ir,3)
     &                  +37.0*ARM(ir,2)-9.0*ARM(ir,1))
            mir = CAT(ir,3)
            ir1 = ir - 2*mir
            ARM(ir1,7) = IFAc(n)*ARM(ir,7)
            IF ( DBLE(mir).LT.-0.1 ) GOTO 120
         ENDDO
      ELSE
         DO ir = 1 , ISMax
            ARM(ir,7) = ARM(ir,5)
     &                  + D2W/24.*(55.0*ARM(ir,4)-59.0*ARM(ir,3)
     &                  +37.0*ARM(ir,2)-9.0*ARM(ir,1))
         ENDDO
      ENDIF
      NPT = NPT + NSW*ISG
      IF ( NPT.GT.0 ) THEN
         IF ( NDIv.EQ.0 ) GOTO 200
         KDIv = KDIv + 1
         IF ( KDIv.LT.NDIv ) GOTO 200
         KDIv = 0
         NPT = NPT + ISG
         IF ( NPT.GT.0 ) GOTO 200
      ENDIF
      NPT = -NPT + 2
      ISG = 1
 200  CALL RESET(ISO)
      IFLg = 1
      i57 = 7
      CALL AMPDER(i57)
      IF ( ISO.EQ.0 ) THEN
         DO n = 1 , NMAx
            ir = NSTart(n) - 1
 220        ir = ir + 1
            ARM(ir,5) = ARM(ir,5)
     &                  + D2W/24.*(9.0*ARM(ir,4)+19.0*ARM(ir,3)
     &                  -5.0*ARM(ir,2)+ARM(ir,1))
            mir = CAT(ir,3)
            ir1 = ir - 2*mir
            ARM(ir1,5) = IFAc(n)*ARM(ir,5)
            IF ( DBLE(mir).LT.-0.1 ) GOTO 220
         ENDDO
      ELSE
         DO ir = 1 , ISMax
            ARM(ir,5) = ARM(ir,5)
     &                  + D2W/24.*(9.0*ARM(ir,4)+19.0*ARM(ir,3)
     &                  -5.0*ARM(ir,2)+ARM(ir,1))
         ENDDO
      ENDIF
      kast = kast + 1
      IFLg = 0
      i57 = 5
      CALL AMPDER(i57)
      IF ( (LAMr(2)+LAMr(3)).NE.0 ) THEN
         IF ( kast.GE.intend ) THEN
            kast = 0
            f = 0.
            DO k = 1 , NMAx
               ihold = IPAth(k)
               IF ( ihold.NE.0 ) THEN
                  hold = ARM(ihold,5) - ARM(ihold,7)
                  rl = DBLE(hold)
                  rim = IMAG(hold)
                  srt = rl*rl + rim*rim
                  f = MAX(f,srt)
               ENDIF
            ENDDO
            f = SQRT(f)/14.
            IF ( f.GT.ACCur .OR. f.LT.ACC50 ) THEN
               IF ( f.LT.ACC50 ) THEN
                  CALL DOUBLE(ISO)
                  D2W = 2.*D2W
                  NSW = 2*NSW
                  intend = (DBLE(intend)+.01)/2.
                  IF ( intend.EQ.0 ) intend = 1
                  IF ( NSW.LT.1 ) THEN
                     NDIv = (DBLE(NDIv)+.01)/2.
                     IF ( NDIv.LT.2 ) THEN
                        NDIv = 0
                        NSW = 1
                     ENDIF
                  ENDIF
               ELSE
                  CALL HALF(ISO)
                  D2W = D2W/2.
                  NSW = (DBLE(NSW)+.01)/2.
                  intend = 2*intend
                  IF ( NSW.LT.1 ) THEN
                     NDIv = 2*NDIv
                     IF ( NDIv.EQ.0 ) NDIv = 2
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      GOTO 100
      END
