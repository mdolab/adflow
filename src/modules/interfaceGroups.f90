!
!      ******************************************************************
!      *                                                                *
!      * File:          interfaceGroups.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-26-2003                                      *
!      * Last modified: 03-22-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module interfaceGroups
!
!      ******************************************************************
!      *                                                                *
!      * Local module which stores the information of sliding mesh      *
!      * interfaces to which this processor contributes.                *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the data type to store the processors having *
!      * faces on the sliding mesh interface.                           *
!      *                                                                *
!      ******************************************************************
!
       type interfaceInfoType

         ! globalSlideId: The global number of the current sliding
         !                mesh interface.

         integer(kind=intType) :: globalSlideId

         ! procContributes: Whether or not the current processor
         !                  contributes to this interface.

         logical :: procContributes

         ! nSlices1: Number of periodic slices to complete a full
         !           rotation for side 1 of the interface.
         ! nSlices2: Idem for side 2.

         integer(kind=intType) :: nSlices1, nSlices2

         ! nProcs1: Number of processors which have faces on side 1 of
         !           the interface.
         ! nProcs2: Idem for side 2.

         integer(kind=intType) :: nProcs1, nProcs2

         ! procs1(nProcs1): The corresponding processor ids of side 1.
         ! procs2(nProcs2): Idem for side 2.

         integer(kind=intType), dimension(:), pointer :: procs1
         integer(kind=intType), dimension(:), pointer :: procs2

         ! Number of processors, my processor number and communicator
         ! for the entire sliding mesh interface.

         integer :: nProcSlide
         integer :: myIdSlide
         integer :: commSlide

         ! Minimum and maximum value of the radius and axial coordinate
         ! in this slide and a scale factor for the radius and axial
         ! coordinate, such that these coordinates can be made
         ! dimensionless, as the angle.

         real(kind=realType) :: rMin, rMax, axMin, axMax, scale

         ! The rotation center.

         real(kind=realType), dimension(3) :: rotCenter

         ! The unit vector of the rotation axis (axial direction) and
         ! the two unit vectors which define the radial plane.

         real(kind=realType), dimension(3) :: rotAxis
         real(kind=realType), dimension(3) :: radVec1, radVec2

       end type interfaceInfoType

       ! nInterfaceGroups: Number of interface groups, such that as many
       !                   interfaces are treated simultaneously, but a
       !                   processor can work on only 1 interface.

       integer(kind=intType) :: nInterfaceGroups

       ! myInterfaces(ninterfaceGroups): The interface info to which
       !                                 this processor possibly
       !                                 contributes.

       type(interfaceInfoType), dimension(:), allocatable :: &
                                                        myInterfaces

       end module interfaceGroups
