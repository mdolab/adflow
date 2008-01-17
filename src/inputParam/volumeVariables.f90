!
!      ******************************************************************
!      *                                                                *
!      * File:          volumeVariables.f90                             *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 03-26-2003                                      *
!      * Last modified: 07-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine volumeVariables(variables)
!
!      ******************************************************************
!      *                                                                *
!      * volumeVariables extracts from the given string the extra       *
!      * volume variables to be written to the solution file.           *
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

       ! Initialize all the volume output variables to .False.

       volWriteMx    = .false.
       volWriteMy    = .false.
       volWriteMz    = .false.
       volWriteRhoe  = .false.
       volWriteTemp  = .false.
       volWriteVort  = .false.
       volWriteVortx = .false.
       volWriteVorty = .false.
       volWriteVortz = .false.

       volWriteCp       = .false.
       volWriteMach     = .false.
       volWriteMachTurb = .false.
       volWritePtotloss = .false.

       volWriteEddyVis      = .false.
       volWriteRatioEddyVis = .false.
       volWriteDist         = .false.

       volWriteResRho  = .false.
       volWriteResMom  = .false.
       volWriteResRhoe = .false.
       volWriteResTurb = .false.

       volWriteBlank = .false.

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

           case ("mx")
             volWriteMx = .true.
             nVarSpecified = nVarSpecified + 1

           case ("my")
             volWriteMy = .true.
             nVarSpecified = nVarSpecified + 1

           case ("mz")
             volWriteMz = .true.
             nVarSpecified = nVarSpecified + 1

           case ("rhoe")
             volWriteRhoe = .true.
             nVarSpecified = nVarSpecified + 1

           case ("temp")
             volWriteTemp = .true.
             nVarSpecified = nVarSpecified + 1

           case ("vort")
             volWriteVort = .true.
             nVarSpecified = nVarSpecified + 1

           case ("vortx")
             volWriteVortx = .true.
             nVarSpecified = nVarSpecified + 1

           case ("vorty")
             volWriteVorty = .true.
             nVarSpecified = nVarSpecified + 1

           case ("vortz")
             volWriteVortz = .true.
             nVarSpecified = nVarSpecified + 1

           case ("cp")
             volWriteCp = .true.
             nVarSpecified = nVarSpecified + 1

           case ("mach")
             volWriteMach = .true.
             nVarSpecified = nVarSpecified + 1

           case ("macht")
             volWriteMachTurb = .true.
             nVarSpecified = nVarSpecified + 1

           case ("ptloss")
             volWritePtotloss = .true.
             nVarSpecified = nVarSpecified + 1

           case ("eddy")
             volWriteEddyVis = .true.
             nVarSpecified = nVarSpecified + 1

           case ("eddyratio")
             volWriteRatioEddyVis = .true.
             nVarSpecified = nVarSpecified + 1

           case ("dist")
             volWriteDist = .true.
             nVarSpecified = nVarSpecified + 1

           case ("resrho")
             volWriteResRho = .true.
             nVarSpecified = nVarSpecified + 1

           case ("resmom")
             volWriteResMom = .true.
             nVarSpecified = nVarSpecified + 1

           case ("resrhoe")
             volWriteResRhoe = .true.
             nVarSpecified = nVarSpecified + 1

           case ("resturb")
             volWriteResTurb = .true.
             nVarSpecified = nVarSpecified + 1

           case ("blank")
             volWriteBlank = .true.
             nVarSpecified = nVarSpecified + 1

           case default
             pos = len_trim(keyword)
             write(errorMessage,"(3a)" ) "Unknown extra volume output &
                                         &variable, ", trim(keyword), &
                                         ", specified"
             if(myID == 0) &
               call terminate("volumeVariables", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)

         end select

       enddo

       ! Set volumeOutSpecified to .true. if variables were specified.
       ! If not, later on the defaults will be set.

       if(nVarSpecified > 0) volumeOutSpecified = .true.

       end subroutine volumeVariables
