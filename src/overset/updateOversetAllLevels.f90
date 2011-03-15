!
!      ******************************************************************
!      *                                                                *
!      * File:          updateOversetAllLevels.f90                      *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 04-08-2005                                      *
!      * Last modified: 10-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine updateOversetAllLevels
!
!      ******************************************************************
!      *                                                                *
!      * UpdateOversetAllLevels updates the overset mesh                *
!      * communication pattern on all grid levels. This routine is      *
!      * typically called when grid parts have been moved, either due   *
!      * to a physical motion of some parts or due to deformation.      *
!      *                                                                *
!      ******************************************************************
!
       use block
       use inputTimeSpectral
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: nLevels, nn, mm, ll
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Update the fine level's communication, which is done by
       ! re-cutting the holes and re-creating the boundaries.

       call oversetComm(1_intType, .false., .false.)

       ! Loop over the coarse levels and update the communication.

       nLevels = ubound(flowDoms,2)
       do nn=2,nLevels
         call oversetComm(nn, .false., .true.)
       enddo

       ! Before heading to the solver, set all the boundary iblanks
       ! for the levels just updated to 0.

       do ll=1,nTimeIntervalsSpectral
         do mm=1,nLevels
           do nn=1,nDom
             call setPointers(nn, mm, ll)
             call changeIblanks(.false., 0_intType)
           end do
         end do
       end do

       end subroutine updateOversetAllLevels
