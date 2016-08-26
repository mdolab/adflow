       subroutine read1to1Conn(cgnsInd, cgnsBase, nZone)
!
!       read1to1Conn reads the 1 to 1 block to block, i.e.             
!       continuous grid lines across block boundaries, connectivities  
!       for the given zone/block.                                      
!
       use constants
       use su_cgns
       use cgnsGrid, only : cgnsDoms, cgnsNDom
       use communication, only : sumb_comm_world, myID
       use utils, only : terminate
       implicit none
!
!      Subroutine arguments
!
       integer, intent(in) :: cgnsInd, cgnsBase, nZone
!
!      Local variables
!
       integer :: cgnsN1to1
       integer :: i, ierr

       integer, dimension(3,2) :: zoneRange, donorRange
       integer, dimension(3)   :: transform

       ! Determine the number of 1 to 1 connectivities stored in the
       ! CGNS file for this zone.

       cgnsN1to1 = cgnsDoms(nZone)%n1to1 - cgnsDoms(nZone)%n1to1General

       ! Loop over the number of 1 to 1 connectivities stored in the
       ! grid file.

       do i=1,cgnsN1to1

         ! Read the 1 to 1 connectivity info from the cgns file.

         call cg_1to1_read_f(cgnsInd, cgnsBase, nZone, i,             &
                             cgnsDoms(nZone)%conn1to1(i)%connectName, &
                             cgnsDoms(nZone)%conn1to1(i)%donorName,   &
                             zoneRange, donorRange, transform, ierr)
         if(ierr /= CG_OK)                 &
           call terminate("read1to1Conn", &
                          "Something wrong when calling cg_1to1_read_f")

         ! Store the zone range and donor range in cgnsDoms.

         cgnsDoms(nZone)%conn1to1(i)%iBeg = zoneRange(1,1)
         cgnsDoms(nZone)%conn1to1(i)%jBeg = zoneRange(2,1)
         cgnsDoms(nZone)%conn1to1(i)%kBeg = zoneRange(3,1)

         cgnsDoms(nZone)%conn1to1(i)%iEnd = zoneRange(1,2)
         cgnsDoms(nZone)%conn1to1(i)%jEnd = zoneRange(2,2)
         cgnsDoms(nZone)%conn1to1(i)%kEnd = zoneRange(3,2)

         cgnsDoms(nZone)%conn1to1(i)%diBeg = donorRange(1,1)
         cgnsDoms(nZone)%conn1to1(i)%djBeg = donorRange(2,1)
         cgnsDoms(nZone)%conn1to1(i)%dkBeg = donorRange(3,1)

         cgnsDoms(nZone)%conn1to1(i)%diEnd = donorRange(1,2)
         cgnsDoms(nZone)%conn1to1(i)%djEnd = donorRange(2,2)
         cgnsDoms(nZone)%conn1to1(i)%dkEnd = donorRange(3,2)

         ! Check the transformation matrices between this zone and the
         ! donor zone and store it in l1, L2 and l3.

         call checkTransform(transform, nZone, i, .true.)

         cgnsDoms(nZone)%conn1to1(i)%l1 = transform(1)
         cgnsDoms(nZone)%conn1to1(i)%l2 = transform(2)
         cgnsDoms(nZone)%conn1to1(i)%l3 = transform(3)

         ! Subface is a normal boundary. Set periodic to .false. and
         ! initialize the periodic data to zero to avoid possible
         ! problems due to uninitialized data.

         cgnsDoms(nZone)%conn1to1(i)%periodic = .false.

         cgnsDoms(nZone)%conn1to1(i)%rotationCenter = zero
         cgnsDoms(nZone)%conn1to1(i)%rotationAngles = zero
         cgnsDoms(nZone)%conn1to1(i)%translation    = zero

       enddo



       end subroutine read1to1Conn
