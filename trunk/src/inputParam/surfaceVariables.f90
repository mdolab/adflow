!
!      ******************************************************************
!      *                                                                *
!      * File:          surfaceVariables.f90                            *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 03-26-2003                                      *
!      * Last modified: 07-13-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine surfaceVariables(variables)
!
!      ******************************************************************
!      *                                                                *
!      * surfaceVariables extracts from the given string the surface    *
!      * variables to be written to the solution file.                  *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       use extraOutput
       use allInputParam
       implicit none
!
!      Subroutine arguments.
!
       character(len=*), intent(inout) :: variables
!
!      Local variables.
!
       integer :: nVarSpecified, pos

       character(len=15) :: keyword
       character(len=maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Convert the string variables to lower case.

       call convertToLowerCase(variables)

       ! Initialize all the surface output variables to .false.

       surfWriteRho  = .false.
       surfWriteP    = .false.
       surfWriteTemp = .false.
       surfWriteVx   = .false.
       surfWriteVy   = .false.
       surfWriteVz   = .false.

       surfWriteCp       = .false.
       surfWritePtotloss = .false.
       surfWriteMach     = .false.

       surfWriteCf    = .false.
       surfWriteCh    = .false.
       surfWriteYplus = .false.
       surfWriteCfx   = .false.
       surfWriteCfy   = .false.
       surfWriteCfz   = .false.

       surfWriteBlank = .false.

       ! Initialize nVarSpecified to 0. This serves as a test
       ! later on.

       nVarSpecified = 0

       ! Loop to extract the info from the string variables.

       do
         ! Condition to exit the loop.

         if(len_trim(variables) == 0) exit

         ! Locate the first occurance of the _ in the string and
         ! determine the string keyword.

         pos = index(variables, "_")
         if(pos == 0) then
           keyword   = variables
           variables = ""
         else
           keyword   = variables(:pos-1)
           variables = variables(pos+1:)
         endif

         ! Check the keyword.

         select case (keyword)
           case ("")
             ! Multiple occurence of "_". Just ignore it.

           case ("rho")
             surfWriteRho = .true.
             nVarSpecified = nVarSpecified + 1

           case ("p")
             surfWriteP = .true.
             nVarSpecified = nVarSpecified + 1

           case ("temp")
             surfWriteTemp = .true.
             nVarSpecified = nVarSpecified + 1

           case ("vx")
             surfWriteVx = .true.
             nVarSpecified = nVarSpecified + 1

           case ("vy")
             surfWriteVy = .true.
             nVarSpecified = nVarSpecified + 1

           case ("vz")
             surfWriteVz = .true.
             nVarSpecified = nVarSpecified + 1

           case ("cp")
             surfWriteCp = .true.
             nVarSpecified = nVarSpecified + 1

           case ("ptloss")
             surfWritePtotloss = .true.
             nVarSpecified = nVarSpecified + 1

           case ("mach")
             surfWriteMach = .true.
             nVarSpecified = nVarSpecified + 1

           case ("cf")
             surfWriteCf = .true.
             nVarSpecified = nVarSpecified + 1

           case ("ch")
             surfWriteCh = .true.
             nVarSpecified = nVarSpecified + 1

           case ("yplus")
             surfWriteYplus = .true.
             nVarSpecified = nVarSpecified + 1

           case ("cfx")
             surfWriteCfx = .true.
             nVarSpecified = nVarSpecified + 1

           case ("cfy")
             surfWriteCfy = .true.
             nVarSpecified = nVarSpecified + 1

           case ("cfz")
             surfWriteCfz = .true.
             nVarSpecified = nVarSpecified + 1

           case ("blank")
             surfWriteBlank = .true.
             nVarSpecified = nVarSpecified + 1

           case default
             pos = len_trim(keyword)
             write(errorMessage,"(3a)") "Unknown surface output &
                                        &variable, ", trim(keyword), &
                                        ", specified"
             if(myID == 0) &
               call terminate("surfaceVariables", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)

         end select

       enddo

       ! Set surfaceOutSpecified to .true. if variables were specified.
       ! If not, later on the defaults will be set.

       if(nVarSpecified > 0) surfaceOutSpecified = .true.

       end subroutine surfaceVariables
