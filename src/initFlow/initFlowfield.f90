!
!      ******************************************************************
!      *                                                                *
!      * File:          initFlowfield.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-07-2003                                      *
!      * Last modified: 03-22-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initFlowfield
!
!      ******************************************************************
!      *                                                                *
!      * initFlowfield initializes the flow field to a uniform flow on  *
!      * the start level grid. Exception may be some turbulence         *
!      * variables, which are initialized a bit smarter.                *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use communication, only : myID
       use inputIteration, only : nCycles, nSGStartup
       use inputPhysics, only : equationMode
       use inputUnsteady, only: nTimeStepsFine
       use iteration, only : nOldSolAvail
       use monitor, only : nTimeStepsRestart, timeUnsteadyRestart
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize nTimeStepsRestart to 0 (no restart
       ! is performed) and allocate the memory for the arrays to store
       ! convergence history. This allocation is only to be done by
       ! processor 0. For an unsteady computation the entire convergence
       ! history is stored and the values after every time step is
       ! stored in a separate array.

       nTimeStepsRestart   = 0
       timeUnsteadyRestart = zero

       nn = nsgStartup + nCycles
       if(equationMode == unsteady) nn = nTimeStepsFine*nn
       if(myID == 0) call allocConvArrays(nn)
       if(equationMode == unsteady .and. myId == 0) &
         call allocTimeArrays(nTimeStepsFine)

       ! Set nOldSolAvail to 1, to indicate that an unsteady
       ! computation should be started with a lower order
       ! time integration scheme, because not enough states
       ! in the past are available.

       nOldSolAvail = 1

       ! Initialize the flow field to uniform flow.

       call setUniformFlow

       end subroutine initFlowfield
