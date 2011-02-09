!
!      ******************************************************************
!      *                                                                *
!      * File:          writeVarsToBuffer.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-20-2005                                      *
!      * Last modified: 03-23-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeVarsToBuffer(buf, var, nn)
!
!      ******************************************************************
!      *                                                                *
!      * writeVarsToBuffer copies the given number of bytes from the    *
!      * array var into the array buf.                                  *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: nn

       character, dimension(*), intent(in)  :: var
       character, dimension(*), intent(out) :: buf
!
!      Local variables.
!
       integer(kind=intType) :: i
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Perform the copy.

       do i=1,nn
         buf(i) = var(i)
       enddo

       end subroutine writeVarsToBuffer
