       subroutine unsteadyHeader
!
!       unsteadyHeader writes a header to stdout when a new time step  
!       is started.                                                    
!
       use constants
       use monitor, only : nTimeStepsRestart, timeUnsteadyRestart, &
            timeUnsteady, timeStepUnsteady
       implicit none
!
!      Local variables
!
       character(len=7)  :: integerString
       character(len=12) :: realString

       ! Write the time step number to the integer string and the
       ! physical time to the real string.

       write(integerString,"(i7)") timeStepUnsteady + &
                                   nTimeStepsRestart
       write(realString,"(e12.5)") timeUnsteady + &
                                   timeUnsteadyRestart

       integerString = adjustl(integerString)
       realString    = adjustl(realString)

       ! Write the header to stdout.

       print "(a)", "#"
       print 100
       print 101
       print 102, trim(integerString), trim(realString)
       print 101
       print 100
       print "(a)", "#"

 100   format("#*************************************************&
              &*************************")
 101   format("#")
 102   format("# Unsteady time step ",a,", physical time ",a, " seconds")

       end subroutine unsteadyHeader
