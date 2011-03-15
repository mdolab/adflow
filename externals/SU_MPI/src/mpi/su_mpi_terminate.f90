!
!      ******************************************************************
!      *                                                                *
!      * File:          su_mpi_terminate.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-17-2005                                      *
!      * Last modified: 04-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine su_mpi_terminate(routineName, errorMessage)
!
!      ******************************************************************
!      *                                                                *
!      * The subroutine su_mpi_terminate writes an error message to     *
!      * standard output terminates the execution of the program.       *
!      *                                                                *
!      ******************************************************************
!
!      Subroutine arguments
!
       character(len=*), intent(in) :: routineName
       character(len=*), intent(in) :: errorMessage
!
!      Local parameter.
!
       integer, parameter :: max_char_line = 55
!
!      Local variables
!
       integer :: len, i2
       logical :: first_time

       character(len=len_trim(errorMessage)) :: message
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

       ! Print a nice error message.

       print "(a)", "#"
       print "(a)", "#=========================== !!! Error !!! &
                    &============================"

       ! Write the header of the error message.

       print "(2a)", "#* Run-time error in procedure ", &
                     trim(routineName)

       ! Loop to write the error message. If the message is too long it
       ! is split over several lines.

       first_time = .true.
       do
         ! Determine the remaining error message to be written.
         ! If longer than the maximum number of characters allowed
         ! on a line, it is attempted to split the message.

         message = adjustl(message)
         len = len_trim(message)
         i2  = min(max_char_line,len)

         if(i2 < len) i2 = index(message(:i2), " ", .true.) - 1
         if(i2 < 0)   i2 = index(message, " ") - 1
         if(i2 < 0)   i2 = len

         ! Write this part of the error message. If it is the first
         ! line of the message some additional stuff is printed.

         if( first_time ) then
           print "(2a)", "#* Error message: ", &
                         trim(message(:i2))
           first_time = .false.
         else
           print "(2a)", "#*                ", &
                         trim(message(:i2))
         endif

         ! Exit the loop if the entire message has been written.

         if(i2 == len) exit

         ! Adapt the string for the next part to be written.

         message = message(i2+1:)

       enddo

       ! write the trailing message.

       print "(a)", "#*"
       print "(a)", "#* Now exiting"
       print "(a)", "#==========================================&
                    &============================"
       print "(a)", "#"

       ! And just stop.

       stop

       end subroutine su_mpi_terminate
