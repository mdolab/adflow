!
!      ******************************************************************
!      *                                                                *
!      * File:          indexHalo.f90                                  *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 01-18-2008                                      *
!      * Last modified: 01-18-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine indexHalo1(level)
!
!      ******************************************************************
!      *                                                                *
!      * whalo1 exchanges all the 1st level internal halo's for the     *
!      * cell centered variables.                                       *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use commSliding
       use communication
       use flowVarRefState
       use inputPhysics
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, ll

       logical :: correctForK
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
 

       ! Exchange the 1 to 1 matching 1st level cell halo's.

       call indexHalo1to1(level, commPatternCell_1st,  &
                      internalCell_1st)

       ! Exchange the 1 to 1 matching 1st level Node halo's.

       call indexHaloNode1to1(level)

       ! Exchange the sliding mesh 1st level cell halo's.

!!$       !!! NOT Implemented Yet
!!$
!!$       mm = ubound(commSlidingCell_1st,1)
!!$       call whaloSliding(level, start, end, commPressure,       &
!!$                         commVarGamma, commLamVis, commEddyVis, &
!!$                         commSlidingCell_1st, intSlidingCell_1st, mm)
!!$
!!$       ! Exchange the mixing plane 1st level cell halo's.
!!$
!!$       call whaloMixing(level, start, end, commPressure, commVarGamma, &
!!$                        commLamVis, commEddyVis, 1_intType)
!!$
!!$       ! Exchange the overset boundary cell data.
!!$
!!$       call wOverset(level, start, end, commPressure,       &
!!$                     commVarGamma, commLamVis, commEddyVis, &
!!$                     commPatternOverset, internalOverset, mm)
!!$
     end subroutine indexHalo1

!      ==================================================================

     subroutine indexHalo2(level)
!
!      ******************************************************************
!      *                                                                *
!      * whalo2 exchanges all the 2nd level internal halo's for the     *
!      * cell centered variables.                                       *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use commSliding
       use communication
       use flowVarRefState
       use inputPhysics
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, ll
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

       logical :: correctForK
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
 

       ! Exchange the 1 to 1 matching 2nd level cell halo's.

       call indexHalo1to1(level, commPatternCell_2nd,  &
                      internalCell_2nd)

       ! Exchange the sliding mesh 2nd level cell halo's.
       
!!$       !!! NOT Implemented Yet
!!$
!!$       mm = ubound(commSlidingCell_2nd,1)
!!$       call whaloSliding(level, start, end, commPressure,       &
!!$                         commVarGamma, commLamVis, commEddyVis, &
!!$                         commSlidingCell_2nd, intSlidingCell_2nd, mm)
!!$
!!$       ! Exchange the mixing plane 2nd level cell halo's.
!!$
!!$       call whaloMixing(level, start, end, commPressure, commVarGamma, &
!!$                        commLamVis, commEddyVis, 2_intType)
!!$
!!$       ! Exchange the overset boundary cell data.
!!$
!!$       call wOverset(level, start, end, commPressure,       &
!!$                     commVarGamma, commLamVis, commEddyVis, &
!!$                     commPatternOverset, internalOverset, mm)


     end subroutine indexHalo2
