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
       implicit none
!
!      Subroutine arguments
!
       integer, intent(in) :: lenExec, lenParam, nArgs
       integer, intent(in) :: sizeC_int

       character(len=lenExec),  intent(in) :: execName
       character(len=lenParam), intent(in) :: paramName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Test if the program was called correctly and store the name of
       ! the parameter file.

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

       ! First part to release the memory.

       call releaseMemoryPart1

       ! Check if for the time spectral method additional solution
       ! files must be written.

       if(equationMode == timeSpectral) then

         if( writeUnsteadyRestartSpectral ) &
           call writeUnsteadyFromSpectral

         if(writeUnsteadyVolSpectral .or. &
            writeUnsteadySurfSpectral)    &
           call writeInterpolFromSpectral

       endif

       ! Second part to release the memory.

       call releaseMemoryPart2

       ! Write the parameters used for this run to stdout.

       call writeInputParam

       end subroutine SUmb_fortran
