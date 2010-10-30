!
!      ******************************************************************
!      *                                                                *
!      * File:          initSuggarInterface.f90                         *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 09-02-2005                                      *
!      * Last modified: 10-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initSuggarInterface
!
!      ******************************************************************
!      *                                                                *
!      * initSuggarInterface initializes the data stored in the modue   *
!      * suggarData, needed for interfacing with SUGGAR. A sorted list  *
!      * of zone names is stored such that grid names from a DCI output *
!      * file can be quickly translated to zone numbers and a Python    *
!      * dictionary can be formed. An array of parallel IO types for    *
!      * each zone is also allocated. In addition, data for the way all *
!      * zones are split up among the processes is also gathered and    *
!      * stored such that the collective communication does not have to *
!      * be repeated at every time step.                                *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsGrid
       use communication
       use suggarData
       implicit none
!
!      Local variables.
!
       integer :: ierr, nSend

       integer, dimension(nProc) :: recvCounts, displs

       integer(kind=intType) :: i, j, nn, mm, nDomTotal

       integer(kind=intType), dimension(nProc)    :: nDomPerProc
       integer(kind=intType), dimension(7,nDom)   :: bufferLocal
       integer(kind=intType), dimension(cgnsNDom) :: ii

       integer(kind=intType), dimension(:,:), allocatable :: buffer
!
!      Function definitions.
!
       integer(kind=intType) :: bsearchStrings
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Copy the number of zones to the suggarData module.

       nZones = cgnsNDom

       ! Allocate memory for the sorted zone names, their zone numbers
       ! and the split data for each zone.

       allocate(zoneNames(nZones), unsortedZone(nZones), &
                splitData(nZones), stat=ierr)
       if (ierr /= 0) &
         call terminate("initSuggarInterface", &
                        "Memory allocation failure for zoneNames, &
                        &unsortedZone, and splitData")

       ! Copy the zone name from the derived data type into zoneNames
       ! and sort the list.

       do i=1,cgnsNDom
         zoneNames(i) = cgnsDoms(i)%zoneName
       enddo

       call qsortStrings(zoneNames, cgnsNDom)

       ! Initialize unsortedZone to -1. This serves as a check during
       ! the search.

       unsortedZone = -1

       ! Find the original zone ids for the sorted zone names.

       do i=1,cgnsNDom
         j = bsearchStrings(cgnsDoms(i)%zoneName, zoneNames, cgnsNDom)

         ! Check if the zone number is not already taken. If this is the
         ! case, this means that the grid file contains two identical
         ! zone names.

         if(unsortedZone(j) /= -1) &
           call terminate("initSuggarInterface", &
                          "Error occurs only when two identical zone &
                          &names are present")

         ! And set the zone number.

         unsortedZone(j) = i
       enddo

       ! Gather the number of domains on each processor.

       call mpi_allgather(nDom, 1, sumb_integer, nDomPerProc, 1, &
                          sumb_integer, SUmb_comm_world, ierr)

       ! Fill the local buffer with the split data for the domains.

       do i=1,nDom
         bufferLocal(1,i) = flowDoms(i,1,1)%cgnsBlockID

         bufferLocal(2,i) = flowDoms(i,1,1)%iBegOr
         bufferLocal(3,i) = flowDoms(i,1,1)%jBegOr
         bufferLocal(4,i) = flowDoms(i,1,1)%kBegOr

         bufferLocal(5,i) = flowDoms(i,1,1)%iEndOr
         bufferLocal(6,i) = flowDoms(i,1,1)%jEndOr
         bufferLocal(7,i) = flowDoms(i,1,1)%kEndOr
       end do

       ! Allocate memory for the global buffer.

       nDomTotal = sum(nDomPerProc)
       allocate(buffer(7,nDomTotal), stat=ierr)
       if (ierr /= 0) &
         call terminate("initSuggarInterface", &
                        "Memory allocation failure for buffer")

       ! Determine the receive counts and displacements and gather the
       ! subblock data for all processors.

       displs(1) = 0
       recvCounts = 7*nDomPerProc(1)
       do i=2,nProc
         recvCounts(i) = 7*nDomPerProc(i)
         displs(i)     = displs(i-1) + recvCounts(i-1)
       end do

       nSend = 7*nDom
       call mpi_allgatherv(bufferLocal, nSend, sumb_integer,         &
                           buffer, recvCounts, displs, sumb_integer, &
                           SUmb_comm_world, ierr)

       ! Initialize the number of subblocks for each zone to 0, then
       ! loop over the total number of domains to counts them.

       do i=1,cgnsNDom
         splitData(i)%nSubBlocks = 0
       end do

       do i=1,nDomTotal
         j = buffer(1,i)
         splitData(j)%nSubBlocks = splitData(j)%nSubBlocks + 1
       end do

       ! Allocate the split data arrays for each zone.

       do i=1,cgnsNDom
         j = splitData(i)%nSubBlocks
         allocate(splitData(i)%proc(j),   splitData(i)%localBlock(j), &
                  splitData(i)%iBegOr(j), splitData(i)%iEndOr(j),     &
                  splitData(i)%jBegOr(j), splitData(i)%jEndOr(j),     &
                  splitData(i)%kBegOr(j), splitData(i)%kEndOr(j),     &
                  stat=ierr)
         if (ierr /= 0) &
           call terminate("initSuggarInterface", &
                          "Memory allocation failure for split arrays")
       end do

       ! Fill in the split arrays for each zone from the buffer data.

       ii = 0
       nn = 0

       do i=1,nProc
         do j=1,nDomPerProc(i)

           ! Update the global counter for the buffer, abbreviate the 
           ! zone number, and update the zone's subblock counter.

           nn     = nn + 1
           mm     = buffer(1,nn)
           ii(mm) = ii(mm) + 1

           ! Store the processor, local block number, and the nodal
           ! range for this subblock.

           splitData(mm)%proc(ii(mm))       = i
           splitData(mm)%localBlock(ii(mm)) = j

           splitData(mm)%iBegOr(ii(mm)) = buffer(2,nn)
           splitData(mm)%jBegOr(ii(mm)) = buffer(3,nn)
           splitData(mm)%kBegOr(ii(mm)) = buffer(4,nn)

           splitData(mm)%iEndOr(ii(mm)) = buffer(5,nn)
           splitData(mm)%jEndOr(ii(mm)) = buffer(6,nn)
           splitData(mm)%kEndOr(ii(mm)) = buffer(7,nn)

         end do
       end do

       ! Deallocate the global buffer.

       deallocate(buffer, stat=ierr)
       if (ierr /= 0) &
         call terminate("initSuggarInterface", &
                        "Deallocation failure for buffer")

       end subroutine initSuggarInterface
