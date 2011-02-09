!
!      ******************************************************************
!      *                                                                *
!      * File:          tmpSliding.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-26-2003                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module tmpSliding
!
!      ******************************************************************
!      *                                                                *
!      * tmpSliding stores the derived datatype for the temporary       *
!      * storage of the distribution of the sliding mesh interfaces     *
!      * over the processors.                                           *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

       ! The definition of the derived datatype tmpSlidingType.

       type tmpSlidingType

         ! nProcs1:         Number of processors that participate to
         !                  part 1 of this sliding mesh interface.
         ! nProcs2:         Number of processors that participate to
         !                  part 2 of this sliding mesh interface.
         ! nFaceMax:        Maximum number of faces on a certain
         !                  processor for this sliding mesh interface.
         ! procs1(nProcs1): The processor ID's of part 1.
         ! procs2(nProcs2): The processor ID's of part 2.


         integer(kind=intType) :: nProcs1, nProcs2
         integer(kind=intType) :: nFaceMax
         integer(kind=intType), dimension(:), pointer :: procs1, procs2

       end type tmpSlidingType

       ! Variables to store the intersections between sliding mesh
       ! interfaces.

       ! nIntersecting(0:cgnsNsliding) : Number of intersections with
       !                                 lower interfaces, cumulative
       !                                 storage format.
       ! intersecting():                 The corresponding interface
       !                                 ID's.

       integer(kind=intType), dimension(:), allocatable :: nIntersecting
       integer(kind=intType), dimension(:), allocatable :: intersecting

       end module tmpSliding
