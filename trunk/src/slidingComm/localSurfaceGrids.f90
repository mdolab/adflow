!
!      ******************************************************************
!      *                                                                *
!      * File:          localSurfaceGrids.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-21-2005                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine localSurfaceGrids(color, gridType)
!
!      ******************************************************************
!      *                                                                *
!      * localSurfaceGrids stores the local surface mesh of the         *
!      * active interface in a useable format for the ADT.              *
!      *                                                                *
!      ******************************************************************
!
       use interfaceGroups
       use localSubfacesMod
       use thisSlide
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: color, gridType
!
!      Local variables.
!
       integer :: ierr, comm

       real(kind=realType) :: thetapMin, thetanMin
       real(kind=realType) :: thetapMax, thetanMax
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the communicator of this sliding mesh a bit easier.

       comm = myInterfaces(color)%commSlide
!
!      ******************************************************************
!      *                                                                *
!      * Construct part 1 of the sliding mesh interface.                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the local number of nodes and quadrilaterals for
       ! part 1 of the sliding interface. Although this information is
       ! already known and stored in the module localSubfaces, it must
       ! be adapted, because in localSubfaces the halo cells are taken
       ! into account and it is possible that some cells and thus nodes
       ! must be duplicated, because of crossing the line theta = pi and
       ! radial singularities.

       call nLocalNodesAndQuads(nMySubfaces1, mySubfaces1, &
                                gridType, nQuad1, nNode1)

       ! Allocate the memory for the variables to store the information
       ! of part 1 of the sliding interface.

       allocate(conn1(4,nQuad1),    subface1(nQuad1), &
                quadID1(nQuad1),    coor1(3,nNode1),  &
                coorInt1(3,nNode1), stat=ierr)
       if(ierr /= 0)                           &
         call terminate("localSurfaceGrids", &
                        "Memory allocation error for part 1")

       ! Determine the data for part 1 of the interface.

       call localNodesAndQuads(nMySubfaces1, mySubfaces1, &
                               nQuad1,       nNode1,      &
                               gridType,     conn1,       &
                               subface1,     quadID1,     &
                               coor1,        coorInt1,    &
                               thetapMin,    thetanMin,   &
                               thetapMax,    thetanMax)

       ! Determine the global minima and maxima of the polar angles for
       ! part 1 of the interface.

       call mpi_allreduce(thetapMin, thetapMin1, 1, sumb_real, &
                          mpi_min,   comm, ierr)
       call mpi_allreduce(thetanMin, thetanMin1, 1, sumb_real, &
                          mpi_min,   comm, ierr)
       call mpi_allreduce(thetapMax, thetapMax1, 1, sumb_real, &
                          mpi_max,   comm, ierr)
       call mpi_allreduce(thetanMax, thetanMax1, 1, sumb_real, &
                          mpi_max,   comm, ierr)
!
!      ******************************************************************
!      *                                                                *
!      * Construct part 2 of the sliding mesh interface.                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the local number of nodes and quadrilaterals for
       ! part 2 of the sliding interface.

       call nLocalNodesAndQuads(nMySubfaces2, mySubfaces2, &
                                gridType, nQuad2, nNode2)

       ! Allocate the memory for the variables to store the information
       ! of part 2 of the sliding interface.

       allocate(conn2(4,nQuad2),    subface2(nQuad2), &
                quadID2(nQuad2),    coor2(3,nNode2),  &
                coorInt2(3,nNode2), stat=ierr)
       if(ierr /= 0)                         &
         call terminate("localSurfaceGrids", &
                        "Memory allocation error for part 2")

       ! Determine the data for part 2 of the interface.

       call localNodesAndQuads(nMySubfaces2, mySubfaces2, &
                               nQuad2,       nNode2,      &
                               gridType,     conn2,       &
                               subface2,     quadID2,     &
                               coor2,        coorInt2,    &
                               thetapMin,    thetanMin,   &
                               thetapMax,    thetanMax)

       ! Determine the global minima and maxima of the polar angles for
       ! part 2 of the interface.

       call mpi_allreduce(thetapMin, thetapMin2, 1, sumb_real, &
                          mpi_min, comm, ierr)
       call mpi_allreduce(thetanMin, thetanMin2, 1, sumb_real, &
                          mpi_min, comm, ierr)
       call mpi_allreduce(thetapMax, thetapMax2, 1, sumb_real, &
                          mpi_max, comm, ierr)
       call mpi_allreduce(thetanMax, thetanMax2, 1, sumb_real, &
                          mpi_max, comm, ierr)

       end subroutine localSurfaceGrids
