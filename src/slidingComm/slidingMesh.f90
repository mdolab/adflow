!
!      ******************************************************************
!      *                                                                *
!      * File:          slidingMesh.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-02-2003                                      *
!      * Last modified: 04-08-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine slidingMesh(level, sps, color)
!
!      ******************************************************************
!      *                                                                *
!      * slidingMesh updates the sliding mesh communication pattern     *
!      * for the given grid level and spectral solution when this       *
!      * processor has to search a sliding mesh interface for the       *
!      * given color.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use interfaceGroups
       use localSubfacesMod
       use thisSlide
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps, color
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check if this processor actually participates in a search.
       ! If not, return immediately.

       if(.not. myInterfaces(color)%procContributes) return

       ! Determine whether or not this processor has subfaces in part 1
       ! and part 2 of the sliding mesh interface.

       partOf1 = .false.
       do nn=1,myInterfaces(color)%nProcs1
         if(myInterfaces(color)%procs1(nn) == myID) partOf1 = .true.
       enddo

       partOf2 = .false.
       do nn=1,myInterfaces(color)%nProcs2
         if(myInterfaces(color)%procs2(nn) == myID) partOf2 = .true.
       enddo

       ! Determine and store the info of the local subfaces of this
       ! sliding mesh interface.

       call mySubfacesSlide(level, sps, color)
!
!      ******************************************************************
!      *                                                                *
!      * Construct the coordinates of the halo nodes.                   *
!      *                                                                *
!      ******************************************************************
!
       ! Store both primary surface grids in a form useable by the ADT
       ! routines. The last argument, 1, indicates that the primary
       ! grid must be created.

       call localSurfaceGrids(color, 1_intType)

       ! Interpolate the halo coordinates of part 1 of the interface.
       ! The second last argument, 1, indicates that coordinates must
       ! be interpolated.

       call interpolateSlide(nMySubfaces1, mySubfaces1, nNode2,     &
                             nQuad2,       thetapMin2,  thetapMax2, &
                             thetanMin2,   thetanMax2,  conn2,      &
                             coor2,        coorInt2,    level,      &
                             sps,          color,       1_intType,  &
                             myInterfaces(color)%nSlices2)

       ! And for part 2.

       call interpolateSlide(nMySubfaces2, mySubfaces2, nNode1,     &
                             nQuad1,       thetapMin1,  thetapMax1, &
                             thetanMin1,   thetanMax1,  conn1,      &
                             coor1,        coorInt1,    level,      &
                             sps,          color,       1_intType,  &
                             myInterfaces(color)%nSlices1)

       ! Release the memory of mySubfaces1, which is only used to
       ! interpolate the coordinates.

       do nn=1,nMySubfaces1
         deallocate(mySubfaces1(nn)%searchNode,          &
                    mySubfaces1(nn)%statusQuad,          &
                    mySubfaces1(nn)%indHaloN,            &
                    mySubfaces1(nn)%connQuad,            &
                    mySubfaces1(nn)%storeQuad,           &
                    mySubfaces1(nn)%coorN,               &
                    mySubfaces1(nn)%coorNInt, stat=ierr)
         if(ierr /= 0)                   &
           call terminate("slidingMesh", &
                          "Deallocation failure of the primary grid &
                          &arrays for mySubfaces1.")
       enddo

       ! Idem for mySubfaces2.

       do nn=1,nMySubfaces2
         deallocate(mySubfaces2(nn)%searchNode,          &
                    mySubfaces2(nn)%statusQuad,          &
                    mySubfaces2(nn)%indHaloN,            &
                    mySubfaces2(nn)%connQuad,            &
                    mySubfaces2(nn)%storeQuad,           &
                    mySubfaces2(nn)%coorN,               &
                    mySubfaces2(nn)%coorNInt, stat=ierr)
         if(ierr /= 0)                   &
           call terminate("slidingMesh", &
                          "Deallocation failure of the primary grid &
                          &arrays for mySubfaces2.")
       enddo

       ! Release the memory of the local surface grids sent to the
       ! ADT routines.

       deallocate(conn1, coor1, coorInt1, subface1, quadID1, &
                  conn2, coor2, coorInt2, subface2, quadID2, &
                  stat=ierr)
       if(ierr /= 0)                   &
         call terminate("slidingMesh", &
                        "Deallocation error for the local surface grids")
!
!      ******************************************************************
!      *                                                                *
!      * Find the donor cells and interpolation weights for the cell    *
!      * centered values.                                               *
!      *                                                                *
!      ******************************************************************
!
       ! Store both dual surface grids in a form useable by the ADT
       ! routines. The last argument, 2, indicates that the dual
       ! grid must be created.

       call localSurfaceGrids(color, 2_intType)

       ! Interpolate the donor cells and weights for the quadrilateral
       ! faces of part 1 of the interface. The second last argument, 2,
       ! indicates that faces must be interpolated.

       call interpolateSlide(nMySubfaces1, mySubfaces1, nNode2,     &
                             nQuad2,       thetapMin2,  thetapMax2, &
                             thetanMin2,   thetanMax2,  conn2,      &
                             coor2,        coorInt2,    level,      &
                             sps,          color,       2_intType,  &
                             myInterfaces(color)%nSlices2)

       ! Idem for part 2.

       call interpolateSlide(nMySubfaces2, mySubfaces2, nNode1,     &
                             nQuad1,       thetapMin1,  thetapMax1, &
                             thetanMin1,   thetanMax1,  conn1,      &
                             coor1,        coorInt1,    level,      &
                             sps,          color,       2_intType,  &
                             myInterfaces(color)%nSlices1)

       ! Release the memory of mySubfaces1, which was only used to
       ! interpolate the faces.

       do nn=1,nMySubfaces1
         deallocate(mySubfaces1(nn)%statusDual,          &
                    mySubfaces1(nn)%storeDual,           &
                    mySubfaces1(nn)%coorQuad, stat=ierr)
         if(ierr /= 0)                   &
           call terminate("slidingMesh", &
                          "Deallocation failure of the dual grid &
                          &arrays for mySubfaces1.")
       enddo

       ! Idem for mySubfaces2.

       do nn=1,nMySubfaces2
         deallocate(mySubfaces2(nn)%statusDual,          &
                    mySubfaces2(nn)%storeDual,           &
                    mySubfaces2(nn)%coorQuad, stat=ierr)
         if(ierr /= 0)                   &
           call terminate("slidingMesh", &
                          "Deallocation failure of the dual grid &
                          &arrays for mySubfaces2.")
       enddo

       ! Release the memory of the local surface grids sent to the
       ! ADT routines, which is not needed anymore.

       deallocate(conn1, coor1, coorInt1, conn2, coor2, coorInt2, &
                  stat=ierr)
       if(ierr /= 0)                   &
         call terminate("slidingMesh", &
                        "Deallocation error for the partial local &
                        &surface grids")

       ! Update the rotation matrix info stored in intSlidingCell1st
       ! and intSlidingCell2nd.

       call updateRotationInfo(level, sps, color)

       ! Update the sliding mesh communication pattern for the currently
       ! active sliding interface.

       call updateSlidingCommPattern(level, sps, color)

       end subroutine slidingMesh
