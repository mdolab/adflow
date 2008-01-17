!
!      ******************************************************************
!      *                                                                *
!      * File:          writeOversetCGNSZone.F90                        *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 08-05-2005                                      *
!      * Last modified: 08-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeOversetCGNSZone(zone, cgnsZone)
!
!      ******************************************************************
!      *                                                                *
!      * writeOversetCGNSZone writes the overset connectivities and     *
!      * holes of the given zone to the cgns file(s).                   *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use inputPhysics
       use iteration
       use outputMod
       implicit none
!
!      Subroutine arguments.
!
       integer,               intent(in) :: cgnsZone
       integer(kind=intType), intent(in) :: zone
!
!      Local variables.
!
       integer(kind=intType) :: mm
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if there are no overset grids.

       if (.not. oversetPresent) return

       ! How and what to write depends on the equation mode. In any case,
       ! only processor 0 writes the data.

       select case (equationMode)
         case (steady, unsteady)

           ! In steady mode, there can only be 1 file to write to. For
           ! unsteady cases, the overset data at old time steps is not
           ! stored and so it is gathered and only written to the 1st
           ! file which is the current time step.

           call gatherOversetCGNSZone(zone, 1_intType)
           if (myID == 0) then
             call writeCGNSOversetConn(zone, 1_intType, cgnsZone)
             call writeCGNSNodalIblank(zone, 1_intType, cgnsZone)
             call releaseOversetCGNSZone(zone)
           end if

         case (timeSpectral)

           ! For time spectral cases where the overset connectivities
           ! are different for each mode, the data needs to be gathered
           ! for that mode before writing to its file. Otherwise, time
           ! can be saved by gathering once for the first mode and then
           ! writing to all files.

           if (changingOverset) then

             do mm = 1,nGridsToWrite
               call gatherOversetCGNSZone(zone, mm)
               if (myID == 0) then
                 call writeCGNSOversetConn(zone, mm, cgnsZone)
                 call writeCGNSNodalIblank(zone, mm, cgnsZone)
                 call releaseOversetCGNSZone(zone)
               end if
             end do

           else

             call gatherOversetCGNSZone(zone, 1_intType)
             if (myID == 0) then
               do mm = 1,nGridsToWrite
                 call writeCGNSOversetConn(zone, mm, cgnsZone)
                 call writeCGNSNodalIblank(zone, mm, cgnsZone)
               end do
               call releaseOversetCGNSZone(zone)
             end if

           end if

       end select

       end subroutine writeOversetCGNSZone

!      ==================================================================

       subroutine gatherOversetCGNSZone(nn, sps)
!
!      ******************************************************************
!      *                                                                *
!      * gatherOversetCGNSZone collects all of overset connectivities   *
!      * and holes for CGNS domain nn and spectral mode sps to process  *
!      * 0 and stores it in the cgnsGrid module.                        *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsGrid
       use communication
       use inputOverset
       use outputMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn, sps
!
!      Local variables.
!
       integer :: ierr, msize

       integer, dimension(nProc) :: recvCounts, displs
       integer, dimension(2)     :: nSend

       integer(kind=intType) :: mm, ii, jj, i, j, k
       integer(kind=intType) :: offset, nBlocks, nCellsGlobal, nInterp

       integer(kind=intType), dimension(2,nProc)  :: nRecv
       integer(kind=intType), dimension(cgnsNDom) :: bndryPerDonor
       integer(kind=intType), dimension(cgnsNDom) :: bndryPerDonorGlobal
       integer(kind=intType), dimension(cgnsNDom) :: donorIndex

       integer(kind=intType), dimension(:,:), allocatable :: buffer
       integer(kind=intType), dimension(:,:), allocatable :: bufferGlobal

       real(kind=realType), dimension(:,:), allocatable :: bufInt
       real(kind=realType), dimension(:,:), allocatable :: bufIntGlobal
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the number of interpolants easier.

       nInterp = nDonorWeights(oversetInterpType)

       ! Store the number of local blocks and the offset in
       ! blocksCGNSblock for this zone a bit easier.

       offset  = nBlocksCGNSblock(nn-1)
       nBlocks = nBlocksCGNSblock(nn) - offset

       ! Count the number of boundary cells that I have for this CGNS
       ! block to be sent. The first index stores the number of
       ! boundary cells, and the second is the number of holes.

       nSend = 0
       do mm = 1,nBlocks
         i = blocksCGNSblock(mm+offset)
         nSend(1) = nSend(1) + flowDoms(i,1,sps)%nCellsOverset
         nSend(2) = nSend(2) + flowDoms(i,1,sps)%nHoles
       end do

       ! Process 0 gathers the message sizes being received.

       call mpi_gather(nSend, 2, sumb_integer, nRecv, 2, &
                       sumb_integer, 0, SUmb_comm_world, ierr)

       ! Allocate a buffer to store the converted boundary indices,
       ! donor CGNS block, donor indices, and interpolants to be sent.

       allocate(buffer(7,nSend(1)), bufInt(nInterp,nSend(1)), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("gatherOversetCGNSZone", &
                        "Memory allocation failure for buffers")

       ! Loop over the local boundaries for the blocks belonging to
       ! this CGNS block.

       k = 0
       bndryPerDonor = 0

       do mm = 1,nBlocks
         i = blocksCGNSblock(mm+offset)
         do j = 1,flowDoms(i,1,sps)%nCellsOverset

           ! Find the CGNS block ID for the donor block using the map
           ! stored in outputMod, store it in the buffer, and update
           ! my counter for the number of cells going to each donor.

           jj = nDomPerProc(flowDoms(i,1,sps)%neighProcOver(j)) &
              + flowDoms(i,1,sps)%neighBlockOver(j)
           ii = IDsBegOrAllDoms(1,jj)

           k = k + 1
           buffer(1,k) = ii
           bndryPerDonor(ii) = bndryPerDonor(ii) + 1

           ! Convert the local boundary indices to the global CGNS
           ! block and store them in the buffer.

           buffer(2,k) = flowDoms(i,1,sps)%ibndry(1,j) &
                       + flowDoms(i,1,1)%ibegor - 2
           buffer(3,k) = flowDoms(i,1,sps)%ibndry(2,j) &
                       + flowDoms(i,1,1)%jbegor - 2
           buffer(4,k) = flowDoms(i,1,sps)%ibndry(3,j) &
                       + flowDoms(i,1,1)%kbegor - 2

           ! Convert the donor indices to CGNS and store them in the
           ! buffer as well. Also copy the interpolants.

           buffer(5:7,k) = flowDoms(i,1,sps)%idonor(:,j) &
                         + IDsBegOrAllDoms(2:4,jj) - 2
           bufInt(:,k)   = flowDoms(i,1,sps)%overint(:,j)

         end do
       end do

       ! Reduce the global number of boundary cells per CGNS donor.

       msize = cgnsNDom
       call mpi_reduce(bndryPerDonor, bndryPerDonorGlobal, msize, &
                       sumb_integer, mpi_sum, 0, SUmb_comm_world, ierr)

       ! Process 0 needs to do some extra work before gathering.

       rootProc1: if (myID == 0) then

         ! Determine the number of CGNS connectivites that are to be
         ! written for this zone, and the index used to help the
         ! distribution into each connectivity.

         j = 0
         do i = 1,cgnsNDom
           if (bndryPerDonorGlobal(i) > 0) then
             j = j + 1
             donorIndex(i) = j
           else
             donorIndex(i) = 0
           end if
         end do

         ! Allocate the connectivites array in the cgnsGrid.

         cgnsDoms(nn)%nOverset = j
         allocate(cgnsDoms(nn)%connOver(j), stat=ierr)
         if(ierr /= 0)                             &
           call terminate("gatherOversetCGNSZone", &
                          "Memory allocation failure for connOver")

         ! Fill in the donor zone name and number for each of the
         ! connectivities along with the number of points.

         j = 0
         do i = 1,cgnsNDom
           if (bndryPerDonorGlobal(i) > 0) then
             j = j + 1
             cgnsDoms(nn)%connOver(j)%npnts = bndryPerDonorGlobal(i)
             cgnsDoms(nn)%connOver(j)%donorBlock = i
             cgnsDoms(nn)%connOver(j)%donorName = cgnsDoms(i)%zoneName
             cgnsDoms(nn)%connOver(j)%connectName = &
                                   "Overset with "//cgnsDoms(i)%zoneName
           end if
         end do

         ! Allocate buffers to store the gathered boundary data.

         nCellsGlobal = sum(nRecv(1,:))
         allocate(bufferGlobal(      7,nCellsGlobal), &
                  bufIntGlobal(nInterp,nCellsGlobal), stat=ierr)
         if(ierr /= 0)                             &
           call terminate("gatherOversetCGNSZone", &
                          "Memory allocation failure for global buffers")

         ! Determine the receive counts and displacements vectors for the
         ! gathering call. Note they need to multiplied by the amount of
         ! info per cell.

         recvCounts = 7*nRecv(1,:)
         displs(1)  = 0
         do i = 2,nProc
           displs(i) = displs(i-1) + recvCounts(i-1)
         end do

       end if rootProc1

       ! Gather both the indical and interpolants data. The recvCounts
       ! and displacement arrays need to be redetermined.

       msize = 7*nSend(1)
       call mpi_gatherv(buffer, msize, sumb_integer, bufferGlobal, &
                        recvCounts, displs, sumb_integer, 0,       &
                        SUmb_comm_world, ierr)

       rootProc2: if (myID == 0) then
         recvCounts = nInterp*nRecv(1,:)
         displs(1)  = 0
         do i = 2,nProc
           displs(i) = displs(i-1) + recvCounts(i-1)
         end do
       end if rootProc2

       msize = nInterp*nSend(1)
       call mpi_gatherv(bufInt, msize, sumb_real, bufIntGlobal, &
                        recvCounts, displs, sumb_real, 0,       &
                        SUmb_comm_world, ierr)

       ! Deallocate the memory taken up by the sending buffers.

       deallocate(buffer, bufInt, stat=ierr)
       if(ierr /= 0)                             &
         call terminate("gatherOversetCGNSZone", &
                        "Deallocation failure for buffers")

       ! Process 0 again requires some extra work.

       rootProc3: if (myID == 0) then

         ! Loop over the connectivities again and allocate the arrays
         ! to store the indices and interpolants. Reset the number of
         ! points to 0 so it can serve as a counter during distribution.

         do j = 1,cgnsDoms(nn)%nOverset
           i = cgnsDoms(nn)%connOver(j)%npnts
           allocate(cgnsDoms(nn)%connOver(j)%ibndry(3,i),       &
                    cgnsDoms(nn)%connOver(j)%idonor(3,i),       &
                    cgnsDoms(nn)%connOver(j)%interp(nInterp,i), &
                    stat=ierr)
           if(ierr /= 0)                             &
             call terminate("gatherOversetCGNSZone", &
                            "Memory allocation failure for conn. data")
           cgnsDoms(nn)%connOver(j)%npnts = 0
         end do

         ! Loop over the global buffer received and distribute the cells
         ! into the correct connectivity.

         do i = 1,nCellsGlobal

           ! Store the index of the connectivity and update the number
           ! of points.

           j = donorIndex(bufferGlobal(1,i))
           k = cgnsDoms(nn)%connOver(j)%npnts + 1
           cgnsDoms(nn)%connOver(j)%npnts = k

           ! Copy the data from the buffers to the connectivity.

           cgnsDoms(nn)%connOver(j)%ibndry(:,k) = bufferGlobal(2:4,i)
           cgnsDoms(nn)%connOver(j)%idonor(:,k) = bufferGlobal(5:7,i)
           cgnsDoms(nn)%connOver(j)%interp(:,k) = bufIntGlobal( : ,i)

         end do

         ! Deallocate the global buffers.

         deallocate(bufferGlobal, bufIntGlobal, stat=ierr)
         if(ierr /= 0)                             &
           call terminate("gatherOversetCGNSZone", &
                          "Deallocation failure for global buffers")

         ! Compute the global number of holes and allocate a buffer.

         nCellsGlobal = sum(nRecv(2,:))
         allocate(bufferGlobal(3,nCellsGlobal), stat=ierr)
         if(ierr /= 0)                             &
           call terminate("gatherOversetCGNSZone", &
                          "Memory allocation failure for global buffer")

         ! Determine the receive counts and displacements vectors for the
         ! gathering call.

         recvCounts = 3*nRecv(2,:)
         displs(1)  = 0
         do i = 2,nProc
           displs(i) = displs(i-1) + recvCounts(i-1)
         end do

       end if rootProc3

       ! Allocate the buffer to store the indices of holes to send.

       allocate(buffer(3,nSend(2)), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("gatherOversetCGNSZone", &
                        "Memory allocation failure for hole buffer")

       ! Loop over my local blocks for this CGNS zone again and extract
       ! the hole indices into the buffer.

       jj = 0
       do mm = 1,nBlocks
         ii = blocksCGNSblock(mm+offset)

         ! Quickly skip the block if it has no holes.

         if (flowDoms(ii,1,sps)%nHoles == 0) cycle

         ! Loop over the owned cells and convert the indices of each
         ! hole and store in the buffer to send.

         do k = 2,flowDoms(ii,1,1)%kl
           do j = 2,flowDoms(ii,1,1)%jl
             do i = 2,flowDoms(ii,1,1)%il

               if (flowDoms(ii,1,sps)%iblank(i,j,k) == 0) then
                 jj = jj + 1
                 buffer(1,jj) = i + flowDoms(ii,1,1)%iBegor - 2
                 buffer(2,jj) = j + flowDoms(ii,1,1)%jBegor - 2
                 buffer(3,jj) = k + flowDoms(ii,1,1)%kBegor - 2
               end if

             end do
           end do
         end do
       end do

       ! Gather all the hole indices into the one hole set that process
       ! 0 created. Note a global buffer isn't needed here.

       msize = 3*nSend(2)
       call mpi_gatherv(buffer, msize, sumb_integer, bufferGlobal, &
                        recvCounts, displs, sumb_integer, 0,       &
                        SUmb_comm_world, ierr)

       ! Deallocate the memory taken up by the sending buffer.

       deallocate(buffer, stat=ierr)
       if(ierr /= 0)                             &
         call terminate("gatherOversetCGNSZone", &
                        "Deallocation failure for hole buffer")

       ! One last time process 0 does more work.

       rootProc4: if (myID == 0) then

         ! If the total number of holes is > 0, then allocate a hole
         ! set and fill in all the information. Otherwise, allocate
         ! with 0 size to be consistent.

         i = min(nCellsGlobal, 1_intType)
         allocate(cgnsDoms(nn)%hole(i), stat=ierr)
         if(ierr /= 0)                             &
           call terminate("gatherOversetCGNSZone", &
                          "Memory allocation failure for hole set")

         createdHoleSet: if (nCellsGlobal > 0) then

           cgnsDoms(nn)%nHoles = 1
           cgnsDoms(nn)%hole(1)%holeName = "Holes"
           cgnsDoms(nn)%hole(1)%npnts    = nCellsGlobal

           allocate(cgnsDoms(nn)%hole(1)%indices(3,nCellsGlobal), &
                    stat=ierr)
           if(ierr /= 0)                             &
             call terminate("gatherOversetCGNSZone", &
                            "Memory allocation failure for hole indices")

           cgnsDoms(nn)%hole(1)%indices = bufferGlobal

         end if createdHoleSet

         ! Deallocate the global buffer.

         deallocate(bufferGlobal, stat=ierr)
         if(ierr /= 0)                             &
           call terminate("gatherOversetCGNSZone", &
                          "Deallocation failure for global buffer")

       end if rootProc4

       end subroutine gatherOversetCGNSZone

!      ==================================================================

       subroutine writeCGNSOversetConn(nn, ind, cgnsZone)
!
!      ******************************************************************
!      *                                                                *
!      * writeCGNSOversetConn writes the overset connectivity and holes *
!      * for CGNS domain nn to the grid file gridNames(ind) with zone   *
!      * ID in the file given by cgnsZone. The connectivity and holes   *
!      * are taken from the module cgnsGrid.                            *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use su_cgns
       use outputMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn, ind
       integer,               intent(in) :: cgnsZone
!
!      Local variables.
!
       integer :: ierr, ii, jj, cgnsInd, cgnsBase, realTypeCGNS

       integer(kind=intType) :: mm

       integer, dimension(2)                :: dimVector
       integer, dimension(:,:), allocatable :: myData, donorData
!
!      Function definition.
!
       integer :: setCGNSRealType
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS

       call terminate("writeCGNSOversetConn", &
                      "Routine should not be called if no cgns support &
                      &is selected.")

#else
       ! Set the cgns real type and store the file and base IDs easier.

       realTypeCGNS = setCGNSRealType()
       cgnsInd      = fileIDs(ind)
       cgnsBase     = cgnsBases(ind)

       ! Loop over the overset connectivities.

       loopOverset: do mm=1,cgnsDoms(nn)%nOverset

         ! Allocate memory for the indices of this zone and the donor
         ! just so that default integers are given to cgns calls.
         ! Copy the indices to the temporary arrays.

         ii = cgnsDoms(nn)%connOver(mm)%npnts
         allocate(myData(3,ii), donorData(3,ii), stat=ierr)
         if(ierr /= 0)                            &
           call terminate("writeCGNSOversetConn", &
                          "Memory allocation failure for overset data")

         myData    = cgnsDoms(nn)%connOver(mm)%ibndry
         donorData = cgnsDoms(nn)%connOver(mm)%idonor

         ! Write the general connectivity.

         call cg_conn_write_f(cgnsInd, cgnsBase, cgnsZone,           &
                              cgnsDoms(nn)%connOver(mm)%connectName, &
                              CellCenter, Overset, PointList,        &
                              ii, myData,                            &
                              cgnsDoms(nn)%connOver(mm)%donorName,   &
                              Structured, CellListDonor, Integer,    &
                              ii, donorData, jj, ierr)
         if(ierr /= all_ok)                       &
           call terminate("writeCGNSOversetConn", &
                          "Something wrong when calling cg_conn_write_f")

         ! Release the memory of the temporary data arrays.

         deallocate(myData, donorData, stat=ierr)
         if(ierr /= 0)                            &
           call terminate("writeCGNSOversetConn", &
                          "Deallocation error for temp overset data")

         ! Goto this connectivity's node in the file and write the
         ! interpolants array.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t",      &
                        cgnsZone, "ZoneGridConnectivity_t", 1,  &
                        "GridConnectivity_t", jj, "end")
         if(ierr /= all_ok)                       &
           call terminate("writeCGNSOversetConn", &
                          "Something wrong when calling cg_goto_f")

         dimVector(1) = ubound(cgnsDoms(nn)%connOver(mm)%interp, 1)
         dimVector(2) = ii

         call cg_array_write_f("InterpolantsDonor", realTypeCGNS, &
                               2, dimVector,                      &
                               cgnsDoms(nn)%connOver(mm)%interp, ierr)
         if(ierr /= all_ok)                     &
           call terminate("writeCGNSOversetConn", &
                         "Something wrong when calling cg_array_write_f")

       end do loopOverset

       ! Loop over the holes.

       loopHoles: do mm=1,cgnsDoms(nn)%nHoles

         ! Allocate memory for the indices just so that default 
         ! integers are given to cgns calls. Then copy the indices to
         ! the temporary array.

         ii = cgnsDoms(nn)%hole(mm)%npnts
         allocate(myData(3,ii), stat=ierr)
         if(ierr /= 0)                            &
           call terminate("writeCGNSOversetConn", &
                          "Memory allocation failure for hole data")

         myData = cgnsDoms(nn)%hole(mm)%indices

         ! Write the hole and then release the temporary memory.

         call cg_hole_write_f(cgnsInd, cgnsBase, cgnsZone,     &
                              cgnsDoms(nn)%hole(mm)%holename,  &
                              CellCenter, PointList, 1,        &
                              ii, myData, jj, ierr)
         if(ierr /= all_ok)                       &
           call terminate("writeCGNSOversetConn", &
                          "Something wrong when calling cg_hole_write_f")

         deallocate(myData, stat=ierr)
         if(ierr /= 0)                            &
           call terminate("writeCGNSOversetConn", &
                          "Deallocation error for temp hole data")

       end do loopHoles

#endif

       end subroutine writeCGNSOversetConn

!      ==================================================================

       subroutine writeCGNSNodalIblank(nn, ind, cgnsZone)
!
!      ******************************************************************
!      *                                                                *
!      * writeCGNSNodalIblank writes a vertex-based iblank array for    *
!      * CGNS domain nn to the grid file gridNames(ind) with zone ID in *
!      * in the file given by cgnsZone. The array is written as an      *
!      * integer solution variable, in a separate solution node from    *
!      * the actual flow solution. This is done for visualization       *
!      * purposes.                                                      *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use cgnsNames
       use su_cgns
       use outputMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn, ind
       integer,               intent(in) :: cgnsZone
!
!      Local variables.
!
       integer :: ierr, cgnsInd, cgnsBase, cgnsSol, ii

       integer(kind=intType) :: mm, l, i, j, k

       integer, dimension(:,:,:), allocatable :: iblank
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS

       call terminate("writeCGNSNodalIblank", &
                      "Routine should not be called if no cgns support &
                      &is selected.")

#else
       ! Store the file and base IDs easier.

       cgnsInd  = fileIDs(ind)
       cgnsBase = cgnsBases(ind)

       ! Allocate memory to store the iblank array to be written and
       ! initialize it to 1 everywhere.

       allocate(iblank(cgnsDoms(nn)%il, &
                       cgnsDoms(nn)%jl, &
                       cgnsDoms(nn)%kl), stat=ierr)
       if(ierr /= 0)                            &
         call terminate("writeCGNSNodalIblank", &
                        "Memory allocation failure for iblank")

       iblank = 1

       ! Loop over the hole sets and their points.

       do mm = 1,cgnsDoms(nn)%nHoles
         do l = 1,cgnsDoms(nn)%hole(mm)%npnts

         ! Store the cell indices a bit easier.

         i = cgnsDoms(nn)%hole(mm)%indices(1,l)
         j = cgnsDoms(nn)%hole(mm)%indices(2,l)
         k = cgnsDoms(nn)%hole(mm)%indices(3,l)

         ! Set the iblank of the 8 nodes making up this hole cell
         ! to 0 (i.e. blank them).

         iblank(i:i+1,j:j+1,k:k+1) = 0

         end do
       end do

       ! Create the flow solution node in the CGNS file and write
       ! the iblank array as an integer solution variable.

       call cg_sol_write_f(cgnsInd, cgnsBase, cgnsZone, &
                           "Nodal Blanks", Vertex, cgnsSol, ierr)
       if(ierr /= all_ok)                       &
         call terminate("writeCGNSNodalIblank", &
                        "Something wrong when calling cg_sol_write_f")

       call cg_field_write_f(cgnsInd, cgnsBase, cgnsZone, cgnsSol, &
                             Integer, cgnsBlank, iblank, ii, ierr)
       if(ierr /= all_ok)                       &
         call terminate("writeCGNSNodalIblank", &
                        "Something wrong when calling cg_field_write_f")

       ! Deallocate the iblank array.

       deallocate(iblank, stat=ierr)
       if(ierr /= 0)                            &
         call terminate("writeCGNSNodalIblank", &
                        "Deallocation failure for iblank")

#endif

       end subroutine writeCGNSNodalIblank
