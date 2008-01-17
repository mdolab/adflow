!
!      ******************************************************************
!      *                                                                *
!      * File:          updateSlidingAllLevels.f90                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-14-2004                                      *
!      * Last modified: 06-12-2004                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine updateSlidingAllLevels
!
!      ******************************************************************
!      *                                                                *
!      * updateSlidingAllLevels updates the sliding mesh                *
!      * communication pattern on all grid levels. This routine is      *
!      * typically called when grid parts have been moved, either due   *
!      * to a physical motion of some parts or due to deformation.      *
!      *                                                                *
!      ******************************************************************
!
       use block
       use iteration
       use interfaceGroups
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
       ! Find out if there are sliding mesh interfaces in the grid.
       ! If not, return immediately.

       if(nInterfaceGroups == 0) return

       ! Loop over the grid levels and call slidingComm.

       nLevels = ubound(flowDoms,2)
       do nn=groundLevel,nLevels
         call slidingComm(nn, .false.)
       enddo

       end subroutine updateSlidingAllLevels
