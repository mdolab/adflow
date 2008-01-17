!
!      ******************************************************************
!      *                                                                *
!      * File:          readOversetHoles.F90                            *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 12-24-2004                                      *
!      * Last modified: 07-11-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readOversetHoles(cgnsInd, cgnsBase, nZone)
!
!      ******************************************************************
!      *                                                                *
!      * ReadOversetHoles reads in the hole sets for this block.        *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use su_cgns
       implicit none
!
!      Subroutine arguments
!
       integer, intent(in) :: cgnsInd, cgnsBase, nZone

#ifdef USE_NO_CGNS

       call terminate("readOversetHoles", &
                      "Routine should not be called if no cgns support &
                      &is selected.")

#else
!
!      Local variables.
!
       integer :: i, nholes, ierr
       integer :: location, ptsetType, npntsets, npnts

       integer, dimension(:,:), allocatable :: myData

       character(len=maxCGNSNameLen) :: holeName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of hole sets.

       call cg_nholes_f(cgnsInd, cgnsBase, nZone, nholes, ierr)
       if(ierr /= all_ok)                   &
         call terminate("readOversetHoles", &
                        "Something wrong when calling cg_nholes_f")

       ! Copy the number of holes to the cgns type.

       cgnsDoms(nZone)%nholes = nholes

       ! Allocate memory for the holes.

       allocate(cgnsDoms(nZone)%hole(cgnsDoms(nZone)%nholes), &
                stat=ierr)
       if(ierr /= 0) &
         call terminate("readOversetHoles", &
                        "Memory allocation failure for hole")

       ! Loop over the holes.

       holes: do i=1,nholes

         ! Read the information of this hole.

         call cg_hole_info_f(cgnsInd, cgnsBase, nZone, i, holeName,  &
                             location, ptsetType, npntsets, npnts, ierr)
         if(ierr /= all_ok)                   &
           call terminate("readOversetHoles", &
                          "Something wrong when calling cg_hole_info_f")

         ! Check if this is a valid hole.

         testHole: if(location         == CellCenter   .and. &
                      ptsetType        == PointList) then

           ! Copy some information.

           cgnsDoms(nZone)%hole(i)%holeName = holeName
           cgnsDoms(nZone)%hole(i)%npnts     = npnts

           ! Allocate the memory for indices.

           allocate(myData(3,npnts), &
                    cgnsDoms(nZone)%hole(i)%indices(3,npnts), &
                    stat=ierr)
           if(ierr /= 0)                            &
             call terminate("readOversetHoles", &
                            "Memory allocation failure for hole data")

           ! Read the indices of the hole.

           call cg_hole_read_f(cgnsInd, cgnsBase, nZone, i, myData, &
                               ierr)
           if(ierr /= all_ok)                   &
             call terminate("readOversetHoles", &
                          "Something wrong when calling cg_hole_read_f")

           ! Store the indices of the current hole.

           cgnsDoms(nZone)%hole(i)%indices = myData

           ! Release the memory of the current hole.

           deallocate(myData, stat=ierr)
           if(ierr /= 0)                            &
             call terminate("readOversetHoles", &
                            "Deallocation error for hole data")

         else testHole

           ! Something wrong with overset cgns data

           call terminate("readOversetHoles", &
                      "Something wrong with overset holes")

         endif testHole

       enddo holes
#endif

       end subroutine readOversetHoles
