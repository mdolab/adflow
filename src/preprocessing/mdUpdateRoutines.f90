!
!      ******************************************************************
!      *                                                                *
!      * File:          mdUpdateRoutines.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-14-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine updateCoordinatesAllLevels
!
!      ******************************************************************
!      *                                                                *
!      * updateCoordinatesAllLevels updates the coordinates of all      *
!      * grid levels, assuming that the owned coordinates of the fine   *
!      * grid are known.                                                *
!      *                                                                *
!      ******************************************************************
!
       use block
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
       ! Determine the halo coordinates of the fine level.

       call xhalo(groundLevel)

       ! Loop over the coarse grid levels; first the owned coordinates
       ! are determined, followed by the halo's.

       nLevels = ubound(flowDoms,2)
       do nn=(groundLevel+1),nLevels
         call coarseOwnedCoordinates(nn)
         call xhalo(nn)
       enddo

       end subroutine updateCoordinatesAllLevels

!      ==================================================================

       subroutine updateMetricsAllLevels
!
!      ******************************************************************
!      *                                                                *
!      * updateMetricsAllLevels recomputes the metrics on all grid      *
!      * levels. This routine is typically called when the coordinates  *
!      * have changed, but the connectivity remains the same, i.e. for  *
!      * moving or deforming mesh problems.                             *
!      *                                                                *
!      ******************************************************************
!
       use block
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
       ! Loop over the grid levels and call metric and checkSymmetry.

       nLevels = ubound(flowDoms,2)
       do nn=groundLevel,nLevels
         call metric(nn)
         call checkSymmetry(nn)
       enddo

       end subroutine updateMetricsAllLevels
