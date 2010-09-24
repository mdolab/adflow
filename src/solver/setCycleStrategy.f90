!
!      ******************************************************************
!      *                                                                *
!      * File:          setCycleStrategy.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-12-2003                                      *
!      * Last modified: 03-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setCycleStrategy
!
!      ******************************************************************
!      *                                                                *
!      * setCycleStrategy sets the multigrid cycling strategy for the   *
!      * multigrid level groundLevel. It is cycle strategy for the      *
!      * fine grid cut off at the current grid level. If the grid level *
!      * is not in the range of the fine grid cycle strategy, cycling   *
!      * will be set to a single grid strategy.                         *
!      *                                                                *
!      ******************************************************************
!
       use inputIteration
       use iteration
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i
       integer(kind=intType) :: thisLevel, maxLevel
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize thisLevel and maxLevel to 1, i.e. the finest grid.

       thisLevel = 1
       maxLevel  = 1

       ! Determine the cycling strategy for groundLevel by looping over
       ! the fine grid cycling strategy and picking the correct entries.

       nStepsCycling = 0
       do i=1,nMGSteps
         thisLevel = thisLevel + cycleStrategy(i)
         maxLevel  = max(maxLevel, thisLevel)

         ! Store this entry in cycling if a) we are on a coarser grid
         ! than groundLevel or b) if we are on groundLevel and
         ! cycleStrategy(i) does not correspond to a restriction,
         ! i.e. 1.

         if((thisLevel == groundLevel .and. cycleStrategy(i) /= 1) .or. &
            thisLevel > groundLevel) then
           nStepsCycling = nStepsCycling + 1
           cycling(nstepsCycling) = cycleStrategy(i)
         endif

         ! Break the loop if a cycle on the current grid level has
         ! been completed.

         if(thisLevel == groundLevel .and. cycleStrategy(i) == -1) exit

       enddo

       ! Take care of the case that groundLevel >= maxLevel.
       ! In this case a single grid strategy is used.

       if(groundLevel >= maxLevel) then
         nStepsCycling = 1
         cycling(1)    = 0
       endif

       ! Check in debug mode if the multigrid strategy created is
       ! a valid one.

       if( debug ) then

         thisLevel = 0
         do i=1,nstepsCycling
           thisLevel = thisLevel + cycling(i)
         enddo

         if(thisLevel /= 0) &
           call terminate("setCyleStrategy", "Invalid strategy created")

       endif

       end subroutine setCycleStrategy
