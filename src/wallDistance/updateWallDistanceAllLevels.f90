!
!      ******************************************************************
!      *                                                                *
!      * File:          updateWallDistanceAllLevels.f90                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-14-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine updateWallDistanceAllLevels
!
!      ******************************************************************
!      *                                                                *
!      * updateWallDistanceAllLevels updates the wall distances for     *
!      * the cell centers on all grid levels. This routine is typically *
!      * called when grid parts have been moved, either due to a        *
!      * physical motion of some parts or due to deformation.           *
!      *                                                                *
!      ******************************************************************
!
       use block
       use inputPhysics
       use iteration
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: nLevels, nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if the rans equations are not solved.

       if(equations /= RANSEquations) return

       ! Loop over the grid levels and call wallDistance.

       nLevels = ubound(flowDoms,2)
       do nn=groundLevel,nLevels
         call wallDistance(nn, .false.)
       enddo

       end subroutine updateWallDistanceAllLevels
