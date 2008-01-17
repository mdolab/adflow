!
!      ******************************************************************
!      *                                                                *
!      * File:          monitorVariables.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-26-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine monitorVariables(variables)
!
!      ******************************************************************
!      *                                                                *
!      * monitorVariables extracts from the given string the variables  *
!      * to be monitored during the convergence.                        *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use cgnsNames
       use monitor
       use allInputParam
       implicit none
!
!      Subroutine arguments.
!
       character(len=*), intent(inout) :: variables
!
!      Local parameter.
!
       integer(kind=intType), parameter :: nVarMax = 21
!
!      Local variables.
!
       integer :: pos, ierr

       character(len=15) :: keyword
       character(len=maxStringLen) :: errorMessage

       character(len=maxCGNSNameLen), dimension(nVarMax) :: tmpNames

       logical :: monDrho
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check if the monitoring names have already been allocated.
       ! This happens when multiple lines for the monitoring variables
       ! are specified in the parameter file. If this happens the last
       ! value is taken and thus release the memory of previously
       ! specified names.

       if( allocated(monNames) ) then
         deallocate(monNames, stat=ierr)
         if(ierr /= 0) call terminate("monitorVariables", &
                                      "Deallocation error for monNames")
       endif

       ! Initialize monDrho, monDturb and showCPU to .false.

       monDrho  = .false.
       monDturb = .false.
       showCPU  = .false.

       ! Initialize nMonSum, nMonMax and nMon to 0.

       nMonSum = 0
       nMonMax = 0
       nMon    = 0

       ! Convert the string variables to lower case.

       call convertToLowerCase(variables)

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

           case ("cpu")    ! only written to stdout.
             showCPU = .true.

           case ("resrho")
             monDrho = .true.
             nMon = nMon + 1; nMonSum = nMonSum + 1
             tmpNames(nMon) = cgnsL2resRho

           case ("resmom")
             nMon = nMon + 3; nMonSum = nMonSum + 3
             tmpNames(nMon-2) = cgnsL2resMomx
             tmpNames(nMon-1) = cgnsL2resMomy
             tmpNames(nMon)   = cgnsL2resMomz

           case ("resrhoe")
             nMon = nMon + 1; nMonSum = nMonSum + 1
             tmpNames(nMon) = cgnsL2resRhoe

           case ("resturb")  ! special case, because the turbulence model
                             ! Is not yet known. See checkMonitor.
             monDturb = .true.

           case ("cl")
             nMon = nMon + 1; nMonSum = nMonSum + 1
             tmpNames(nMon) = cgnsCl

           case ("clp")
             nMon = nMon + 1; nMonSum = nMonSum + 1
             tmpNames(nMon) = cgnsClp

           case ("clv")
             nMon = nMon + 1; nMonSum = nMonSum + 1
             tmpNames(nMon) = cgnsClv

           case ("cd")
             nMon = nMon + 1; nMonSum = nMonSum + 1
             tmpNames(nMon) = cgnsCd

           case ("cdp")
             nMon = nMon + 1; nMonSum = nMonSum + 1
             tmpNames(nMon) = cgnsCdp

           case ("cdv")
             nMon = nMon + 1; nMonSum = nMonSum + 1
             tmpNames(nMon) = cgnsCdv

           case ("cfx")
             nMon = nMon + 1; nMonSum = nMonSum + 1
             tmpNames(nMon) = cgnsCfx

           case ("cfy")
             nMon = nMon + 1; nMonSum = nMonSum + 1
             tmpNames(nMon) = cgnsCfy

           case ("cfz")
             nMon = nMon + 1; nMonSum = nMonSum + 1
             tmpNames(nMon) = cgnsCfz

           case ("cmx")
             nMon = nMon + 1; nMonSum = nMonSum + 1
             tmpNames(nMon) = cgnsCmx

           case ("cmy")
             nMon = nMon + 1; nMonSum = nMonSum + 1
             tmpNames(nMon) = cgnsCmy

           case ("cmz")
             nMon = nMon + 1; nMonSum = nMonSum + 1
             tmpNames(nMon) = cgnsCmz

           case ("hdiff")
             nMon = nMon + 1; nMonMax = nMonMax + 1
             tmpNames(nMon) = cgnsHdiffMax

           case ("mach")
             nMon = nMon + 1; nMonMax = nMonMax + 1
             tmpNames(nMon) = cgnsMachMax

           case ("yplus")
             nMon = nMon + 1; nMonMax = nMonMax + 1
             tmpNames(nMon) = cgnsYplusMax

           case ("eddyv")
             nMon = nMon + 1; nMonMax = nMonMax + 1
             tmpNames(nMon) = cgnsEddyMax

           case default
             write(errorMessage,"(3a)") "Unknown monitoring variable, ", &
                                        trim(keyword), ", specified"
             if(myID == 0) &
               call terminate("monitorVariables", errorMessage)
             call mpi_barrier(SUmb_comm_world, ierr)

         end select

       enddo

       ! If the density residual was not specified to be monitored,
       ! add it to tmpNames.

       if(.not. monDrho) then
         nMon = nMon + 1; nMonSum = nMonSum + 1
         tmpNames(nMon) = cgnsL2resRho
       endif

       ! Allocate the memory for monNames. If the turbulent residuals
       ! must be monitored allocate some extra place.

       pos = nMon
       if( monDturb ) pos = nMon + 4
       allocate(monNames(pos), stat=ierr)
       if(ierr /= 0)                         &
         call terminate("monitorVariables", &
                        "Memory allocation failure for monNames")

       ! Copy the monitoring names into monNames.

       do pos=1,nMon
         monNames(pos) = tmpNames(pos)
       enddo

       ! Set monitorSpecified to .true. to indicate that monitoring
       ! variables have been specified.

       monitorSpecified = .true.

       end subroutine monitorVariables
