!
!      ******************************************************************
!      *                                                                *
!      * File:          isoVariables.f90                                *
!      * Author:        Gaetan Kenway                                   *
!      * Starting date: 07-21-2013                                      *
!      * Last modified: 07-21-2013                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine isoVariables(variables)
!
!      ******************************************************************
!      *                                                                *
!      * isoVariables extracts from the given string the extra          *
!      * iso surface variables to be written to the solution file.      *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       use extraOutput
       use allInputParam
       use utils, only : convertToLowerCase, terminate
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
 
       ! Initialize all the iso output variables to .False.
       isoWriteRho   = .false.
       isoWriteVx    = .false.
       isoWriteVy    = .false.
       isoWriteVz    = .false.
       isoWriteP     = .false.
       isoWriteTurb  = .false.

       isoWriteMx    = .false.
       isoWriteMy    = .false.
       isoWriteMz    = .false.
       isoWriteRhoe  = .false.
       isoWriteTemp  = .false.
       isoWriteVort  = .false.
       isoWriteVortx = .false.
       isoWriteVorty = .false.
       isoWriteVortz = .false.

       isoWriteCp       = .false.
       isoWriteMach     = .false.
       isoWriteMachTurb = .false.
       isoWritePtotloss = .false.

       isoWriteEddyVis      = .false.
       isoWriteRatioEddyVis = .false.
       isoWriteDist         = .false.

       isoWriteResRho  = .false.
       isoWriteResMom  = .false.
       isoWriteResRhoe = .false.
       isoWriteResTurb = .false.

       isoWriteShock = .false.
       isoWriteFilteredShock = .false.
       
       isoWriteBlank = .false.

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

           case("rho")
              isoWriteRho   = .true.
              nVarSpecified = nVarSpecified + 1
                 
           case("vx")
              isoWriteVx   = .true.
              nVarSpecified = nVarSpecified + 1

           case("vy")
              isoWriteVy   = .true.
              nVarSpecified = nVarSpecified + 1
              
           case("vz")
              isoWriteVz   = .true.
              nVarSpecified = nVarSpecified + 1
              
           case("P")
              isoWriteP   = .true.
              nVarSpecified = nVarSpecified + 1

           case("turb")
              isoWriteTurb = .true.
              nVarSpecified = nVarSpecified + 1

           case ("mx")
             isoWriteMx = .true.
             nVarSpecified = nVarSpecified + 1

           case ("my")
             isoWriteMy = .true.
             nVarSpecified = nVarSpecified + 1

           case ("mz")
             isoWriteMz = .true.
             nVarSpecified = nVarSpecified + 1

           case ("rvx")
             isoWriteRVx = .true.
             nVarSpecified = nVarSpecified + 1

           case ("rvy")
             isoWriteRVy = .true.
             nVarSpecified = nVarSpecified + 1

           case ("rvz")
             isoWriteRVz = .true.
             nVarSpecified = nVarSpecified + 1

           case ("rhoe")
             isoWriteRhoe = .true.
             nVarSpecified = nVarSpecified + 1

           case ("temp")
             isoWriteTemp = .true.
             nVarSpecified = nVarSpecified + 1

           case ("vort")
             isoWriteVort = .true.
             nVarSpecified = nVarSpecified + 1

           case ("vortx")
             isoWriteVortx = .true.
             nVarSpecified = nVarSpecified + 1

           case ("vorty")
             isoWriteVorty = .true.
             nVarSpecified = nVarSpecified + 1

           case ("vortz")
             isoWriteVortz = .true.
             nVarSpecified = nVarSpecified + 1

           case ("cp")
             isoWriteCp = .true.
             nVarSpecified = nVarSpecified + 1

           case ("mach")
             isoWriteMach = .true.
             nVarSpecified = nVarSpecified + 1
          
           case ("rmach")
             isoWriteRMach = .true.
             nVarSpecified = nVarSpecified + 1

           case ("macht")
             isoWriteMachTurb = .true.
             nVarSpecified = nVarSpecified + 1

           case ("ptloss")
             isoWritePtotloss = .true.
             nVarSpecified = nVarSpecified + 1

           case ("eddy")
             isoWriteEddyVis = .true.
             nVarSpecified = nVarSpecified + 1

           case ("eddyratio")
             isoWriteRatioEddyVis = .true.
             nVarSpecified = nVarSpecified + 1

           case ("dist")
             isoWriteDist = .true.
             nVarSpecified = nVarSpecified + 1

           case ("resrho")
             isoWriteResRho = .true.
             nVarSpecified = nVarSpecified + 1

           case ("resmom")
             isoWriteResMom = .true.
             nVarSpecified = nVarSpecified + 1

           case ("resrhoe")
             isoWriteResRhoe = .true.
             nVarSpecified = nVarSpecified + 1

           case ("resturb")
             isoWriteResTurb = .true.
             nVarSpecified = nVarSpecified + 1

           case ("blank")
             isoWriteBlank = .true.
             nVarSpecified = nVarSpecified + 1

          case("shock")
             isoWriteShock = .true.
             nVarSpecified = nVarSpecified + 1
             
          case("filteredshock")
             isoWriteFilteredShock = .true.
             nVarSpecified = nVarSpecified + 1

             
           case default
             pos = len_trim(keyword)
             write(errorMessage,"(3a)" ) "Unknown extra iso output &
                                         &variable, ", trim(keyword), &
                                         ", specified"
             if(myID == 0) &
               call terminate("isoVariables", errorMessage)
             call mpi_barrier(SUmb_comm_world, pos)

         end select

       enddo

       ! Set this to true regardless...it is possible no varibles were
       ! specified
       isoOutSpecified = .true.

     end subroutine isoVariables
