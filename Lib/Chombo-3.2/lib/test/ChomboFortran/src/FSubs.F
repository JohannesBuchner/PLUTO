      subroutine TESTCALL_2( int1,status )
      integer int1,status
      if( int1 .NE. 23 )then
        status = 1
        write(6,*) 'error: TESTCALL_2: got int arg 1 = ',int1
     &            ,' instead of 23.'
      else
        status = 0
      endif
      int1 = 32
      return
      end
