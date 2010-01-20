!
!      ******************************************************************
!      *                                                                *
!      * File:          SUmb_fortran.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-14-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine SUmb_fortran(execName, lenExec, paramName, &
                               lenParam, nArgs, sizeC_int)
!
!      ******************************************************************
!      *                                                                *
!      * SUmb_fortran is the fortran main, which is called from the     *
!      * C main. The reason for splitting the main in a C and a Fortran *
!      * part is that in c there is a standard to specify command line  *
!      * arguments; in Fortran there is not.                            *
!      * SUmb_fortran is a high level interface to the libraries        *
!      * that perform the individual subtasks.                          *
!      *                                                                *
!      ******************************************************************
!
       use inputPhysics
       use inputTimeSpectral
       use inputADjoint
       use inputTSStabDeriv
       use communication
       implicit none
!
!      Subroutine arguments
!
       integer, intent(in) :: lenExec, lenParam, nArgs
       integer, intent(in) :: sizeC_int

       character(len=lenExec),  intent(in) :: execName
       character(len=lenParam), intent(in) :: paramName

       ! variable for timing
       real(kind=realType) :: timer(10)

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Test if the program was called correctly and store the name of
       ! the parameter file.

       !set reference time
!!$       if(myID==0)then
!!$          print *, 'Set reference time ....'
       call cpu_time(timer(1))
!!$       endif

       call initExec(execName, lenExec, paramName, lenParam, &
                     nArgs, sizeC_int)

       ! Read the parameter file.

       call readParamFile

       ! Partition the blocks and read the grid.

       call partitionAndReadGrid

       ! Perform the preprocessing task.

       call preprocessing

       ! Initialize of the flow variables.

       call initFlow

       ! Solve the equations.

       call solver

       if(myID==0)then
          call cpu_time(timer(2))
          print *, "       time to end of solver = ", timer(2)-timer(1) 
       endif

       ! Solve ADjoint
       
       if(solveADjoint) then
          call solverADjoint
          
          if(myID==0)then
             call cpu_time(timer(3))
             print *, "       time for adjoint = ", timer(3)-timer(2)
          endif
       end if
       !print *,'tsstability',TSStability
!stop
       if(TSStability)then
          call stabilityDerivativeDriver
          !call computeTSDerivatives
          
          if(myID==0)then
             call cpu_time(timer(4))
             print *, "       time for tsDerivatives = ", timer(4)-timer(3)
          endif
       end if

       ! First part to release the memory.
       !print *,'releasing memory'
       call releaseMemoryPart1
       !print *,'memory released'
       ! Check if for the time spectral method additional solution
       ! files must be written.
       
       if(equationMode == timeSpectral) then
        !  print *,'writing ts 1'
         if( writeUnsteadyRestartSpectral ) &
           call writeUnsteadyFromSpectral
         !print *,'writing ts 2'
         if(writeUnsteadyVolSpectral .or. &
            writeUnsteadySurfSpectral)    &
           call writeInterpolFromSpectral
         !print *,'written'
       endif

       ! Second part to release the memory.
       !print *,'releasing memory 2'
       call releaseMemoryPart2
       !print *,'memory released 2'
       ! Write the parameters used for this run to stdout.

       call writeInputParam

       end subroutine SUmb_fortran
