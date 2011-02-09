!
!      ******************************************************************
!      *                                                                *
!      * File:          thisSlide.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-14-2003                                      *
!      * Last modified: 02-10-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       module thisSlide
!
!      ******************************************************************
!      *                                                                *
!      * Local module which stores the information of the currently     *
!      * active sliding mesh interface.                                 *
!      *                                                                *
!      ******************************************************************
!
       use adtAPI
       use precision
       implicit none
       save

       ! Logicals to indicate whether or not this processor has subfaces
       ! on part 1 and part 2 of the sliding mesh interface.

       logical :: partOf1, partOf2

       ! nQuad1: # of quadrilateral faces on part 1 of the interface.
       ! nQuad2: Idem for part 2.
       ! nNode1: # of points on part 1 of the interface.
       ! nNode2: Idem for part 2.

       integer(kind=adtIntType) :: nQuad1, nQuad2
       integer(kind=adtIntType) :: nNode1, nNode2

       ! thetapMin1: The minimum positive polar angle that is present
       !             in part 1 of the interface.
       ! thetapMin2: Idem for part 2.
       ! thetapMax1: The maximum positive polar angle that is present
       !             in part 1 of the interface.
       ! thetapMax2: Idem for part 2.
       ! thetanMin1: The minimum negative polar angle that is present
       !             in part 1 of the interface.
       ! thetanMin2: Idem for part 2.
       ! thetanMax1: The maximum negative polar angle that is present
       !             in part 1 of the interface.
       ! thetanMax2: Idem for part 2.

       real(kind=realType) :: thetapMin1, thetapMin2
       real(kind=realType) :: thetapMax1, thetapMax2
       real(kind=realType) :: thetanMin1, thetanMin2
       real(kind=realType) :: thetanMax1, thetanMax2

       ! conn1(4,nQuad1): The 4 nodes of the quadrilateral of part 1.
       ! conn2(4,nQuad2): Idem for part 2.

       integer(kind=adtIntType), dimension(:,:), allocatable :: conn1
       integer(kind=adtIntType), dimension(:,:), allocatable :: conn2

       ! subface1(nQuad1): Subface ID of the quad of part 1 to which
       !                   it belongs.
       ! subface2(nQuad2): Idem for part 2.
       ! quadID1(nQuad1):  Corresponding ID in the subface of the
       !                   quad of part 1.
       ! quadID2(nQuad2):  Idem for part 2.

       integer(kind=intType), dimension(:), allocatable :: subface1
       integer(kind=intType), dimension(:), allocatable :: subface2
       integer(kind=intType), dimension(:), allocatable :: quadID1
       integer(kind=intType), dimension(:), allocatable :: quadID2

       ! coor1(3,nNode1): The coordinates of the points of part 1 of
       !                  the interface.
       ! coor2(3,nNode2): Idem for part 2.

       real(kind=adtRealType), dimension(:,:), allocatable :: coor1
       real(kind=adtRealType), dimension(:,:), allocatable :: coor2

       ! coorInt1(3,nNode1): The coordinates of the nodes of part 1
       !                     one layer into the block. These are
       !                     needed for the interpolation of the
       !                     coordinates of the halo nodes.
       ! coorInt2(3,nNode2): Idem for part 2.

       real(kind=adtRealType), dimension(:,:), allocatable :: coorInt1
       real(kind=adtRealType), dimension(:,:), allocatable :: coorInt2

       end module thisSlide
