 
C----------------------------------------------------------------------
 
      SUBROUTINE OPENF
      IMPLICIT NONE
      INTEGER*4 i , j , k
      CHARACTER name*60 , opt1*20 , opt2*20
 100  READ * , i , j , k
      IF ( i.EQ.0 ) RETURN
      IF ( j.EQ.1 ) opt1 = 'OLD'
      IF ( j.EQ.2 ) opt1 = 'NEW'
      IF ( j.EQ.3 ) opt1 = 'UNKNOWN'
      IF ( k.EQ.1 ) opt2 = 'FORMATTED'
      IF ( k.EQ.2 ) opt2 = 'UNFORMATTED'
      READ 99001 , name
99001 FORMAT (A)
      OPEN (i,IOSTAT=k,FILE=name,STATUS=opt1,FORM=opt2)
      IF ( k.EQ.0 ) WRITE (6,99002) 'OPENED ' , name
99002 FORMAT (1X,2A)
      WRITE (6,99003) ' IO-num = ' , i , opt1 , opt2
99003 FORMAT (1X,A,I4,2(1x,A))
      IF ( k.EQ.0 ) GOTO 100
      WRITE (6,99004) 'PROBLEMS OPENING ' , name , k
99004 FORMAT (A,A,I6)
      END
