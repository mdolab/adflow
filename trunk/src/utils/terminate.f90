!
!      ******************************************************************
!      *                                                                *
!      * File:          terminate.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-10-2002                                      *
!      * Last modified: 10-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine terminate(routineName, errorMessage)
!
!      ******************************************************************
!      *                                                                *
!      * terminate writes an error message to standard output and       *
!      * terminates the execution of the program.                       *
!      *                                                                *
!      ******************************************************************
!
       use precision
       use communication
       use constants
       implicit none
!
!      Subroutine arguments
!
       character(len=*), intent(in) :: routineName
       character(len=*), intent(in) :: errorMessage
!
!      Local parameter
!
       integer, parameter :: maxCharLine = 55
!
!      Local variables
!
       integer :: ierr, len, i2
       logical :: firstTime

       character(len=len_trim(errorMessage)) :: message
       character(len=8) :: integerString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Copy the errorMessage into message. It is not possible to work
       ! with errorMessage directly, because it is modified in this
       ! routine. Sometimes a constant string is passed to this routine
       ! and some compilers simply fail then.

       message = errorMessage

       ! Print a nice error message. In case of a parallel executable
       ! also the processor id is printed.

       print "(a)", "#"
       print "(a)", "#=========================== !!! Error !!! &
                    &============================"

       if(.not. SU_MPI_isSequential) then
         write(integerString,"(i8)") myID
         integerString = adjustl(integerString)

         print "(2a)", "#* Terminate called by processor ", &
                       trim(integerString)
       endif

       ! Write the header of the error message.

       print "(2a)", "#* Run-time error in procedure ", &
                     trim(routineName)

       ! Loop to write the error message. If the message is too long it
       ! is split over several lines.

       firstTime = .true.
       do
         ! Determine the remaining error message to be written.
         ! If longer than the maximum number of characters allowed
         ! on a line, it is attempted to split the message.

         message = adjustl(message)
         len = len_trim(message)
         i2  = min(maxCharLine,len)

         if(i2 < len) i2 = index(message(:i2), " ", .true.) - 1
         if(i2 < 0)   i2 = index(message, " ") - 1
         if(i2 < 0)   i2 = len

         ! Write this part of the error message. If it is the first
         ! line of the message some additional stuff is printed.

         if( firstTime ) then
           print "(2a)", "#* Error message: ", &
                         trim(message(:i2))
           firstTime = .false.
         else
           print "(2a)", "#*                ", &
                         trim(message(:i2))
         endif

         ! Exit the loop if the entire message has been written.

         if(i2 == len) exit

         ! Adapt the string for the next part to be written.

         message = message(i2+1:)

       enddo

       ! Write the trailing message.

       print "(a)", "#*"
       print "(a)", "#* Now exiting"
       print "(a)", "#==========================================&
                    &============================"
       print "(a)", "#"

       ! Call abort and stop the program. This stop should be done in
       ! abort, but just to be sure.

       call mpi_abort(SUmb_comm_world, 1, ierr)
       stop

       end subroutine terminate
