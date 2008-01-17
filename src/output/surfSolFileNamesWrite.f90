!
!      ******************************************************************
!      *                                                                *
!      * File:          surfSolFileNamesWrite.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-11-2005                                      *
!      * Last modified: 01-24-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine surfSolFileNamesWrite
!
!      ******************************************************************
!      *                                                                *
!      * surfSolFileNamesWrite determines the names and number of       *
!      * surface solution files to be written.                          *
!      *                                                                *
!      ******************************************************************
!
       use inputIO
       use inputPhysics
       use inputTimeSpectral
       use monitor
       use outputMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn

       character(len=7) :: intString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! In contrast to the grids and volume solutions, possible states
       ! in the past don't need to be written for the surface. Therefore
       ! the memory allocation can be done independent of the test of
       ! the equation mode we are solving for.

       allocate(surfSolFileNames(nTimeIntervalsSpectral), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("surfSolFileNamesWrite", &
                        "Memory allocation failure for surfSolFileNames")

       ! Set the number of surface solution files to be written.

       if( writeSurface ) then
         nSurfSolToWrite = nTimeIntervalsSpectral
       else
         nSurfSolToWrite = 0
       endif
 
       ! Determine the name(s) of the solution file(s), depending on
       ! the situation.

       select case (equationMode)

         case (steady)

           ! Steady state computation. Possible previous files will
           ! be overwritten.

           surfSolFileNames(1) = surfaceSolFile

         !===============================================================

         case (unsteady)

           ! Unsteady computation. A suffix is added depending on the
           ! time step.

           write(intString,"(i7)") timeStepUnsteady + nTimeStepsRestart
           intString = adjustl(intString)

           surfSolFileNames(1) = trim(surfaceSolFile)//"&
                                 &Timestep"//trim(intString)

         !===============================================================

         case (timeSpectral)

           ! Time spectral computation. A suffix is added depending on
           ! the time instance.

           do nn=1,nTimeIntervalsSpectral
             write(intString,"(i7)") nn
             intString = adjustl(intString)

             surfSolFileNames(nn) = trim(surfaceSolFile)//"&
                                    &Spectral"//trim(intString)
           enddo

       end select

       end subroutine surfSolFileNamesWrite
