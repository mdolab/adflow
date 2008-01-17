!
!      ******************************************************************
!      *                                                                *
!      * File:          loadSuggarDCIFile.f90                           *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 08-28-2005                                      *
!      * Last modified: 10-11-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine loadSuggarDCIFile(dciFile, sps)
!
!      ******************************************************************
!      *                                                                *
!      * loadSuggarDCIFile reads the given SUGGAR ouput DCI file into   *
!      * the given spectral mode. Any previous overset boundary and     *
!      * hole data for the CGNS zones affected by the DCI file are      *
!      * erradicated and replaced with the data read from the file.     *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsGrid
       use communication
       use suggarData
       implicit none
!
!      Local parameters.
!
       integer, parameter :: readUnit = 24
!
!      Subroutine arguments.
!
       character(len=*), intent(in) :: dciFile

       integer(kind=intType), intent(in) :: sps
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, k, n, m, nGrids, nPairs, nIblank
       integer(kind=intType) :: nBlocks, igrid, izone, idom, ii
       integer(kind=intType) :: myIndex, donorGrid, nDonor, donorZone
 
       integer(kind=intType), dimension(nDom) :: blockNum
       integer(kind=intType), dimension(8)    :: donorIndices

       integer(kind=intType), dimension(:), allocatable :: zoneNumber
       integer(kind=intType), dimension(:,:,:), allocatable :: zoneBlank
       integer(kind=intType), dimension(:,:,:), pointer :: iblank

       real(kind=realType), dimension(8) :: donorInterp

       logical, dimension(cgnsNDom) :: readZone

       character(len=64)             :: string
       character(len=maxStringLen)   :: errorMessage
       character(len=maxCGNSNameLen) :: gridName
!
!      Function definitions.
!
       integer(kind=intType) :: bsearchStrings
!
!      Interfaces
!
       interface
         subroutine reallocateInteger(intArray, newSize, oldSize, &
                                       alwaysFreeMem)
           use precision
           implicit none

           integer(kind=intType), dimension(:), pointer :: intArray
           integer(kind=intType), intent(in) :: newSize, oldSize
           logical, intent(in) :: alwaysFreeMem
         end subroutine reallocateInteger

         !===============================================================

         subroutine reallocateInteger2(intArray,            &
                                        newSize1, newSize2, &
                                        oldSize1, oldSize2, &
                                        alwaysFreeMem)
           use precision
           implicit none

           integer(kind=intType), dimension(:,:), pointer :: intArray
           integer(kind=intType), intent(in) :: newSize1, newSize2, &
                                                   oldSize1, oldSize2
           logical, intent(in) :: alwaysFreeMem
         end subroutine reallocateInteger2

         !===============================================================

         subroutine reallocateReal2(realArray,           &
                                     newSize1, newSize2, &
                                     oldSize1, oldSize2, &
                                     alwaysFreeMem)
           use precision
           implicit none

           real(kind=realType), dimension(:,:), pointer :: realArray
           integer(kind=intType), intent(in) :: newSize1, newSize2, &
                                                   oldSize1, oldSize2
           logical, intent(in) :: alwaysFreeMem
         end subroutine reallocateReal2
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Write a message that the file is being loaded.

       if (myID == 0) then
         print "( a)", "#"
         print "(3a)", "# Loading SUGGAR DCI file: ", &
                       trim(adjustl(dciFile)), "..."
       end if

       ! Open the SUGGAR DCI file for reading, and set the format for
       ! reading keyword lines from the file.

       open (unit=readUnit, file=dciFile, action="read", status="old", &
             iostat=ierr)
       if (ierr /= 0)                        &
         call terminate("loadSuggarDCIFile", &
                        "Error opening SUGGAR DCI file")

 100   format (a64)

       ! Read the number of grids in the file.

       do
         read (readUnit,100) string
         if (string(1:13) == ":number_grids") exit
       end do
       read (string(14:),*) nGrids

       ! The DCI file contains a mapping at the beginning of the file to
       ! state which index refers to which grid name. Allocate memory
       ! for this mapping.

       allocate(zoneNumber(nGrids), stat=ierr)
       if (ierr /= 0)                        &
         call terminate("loadSuggarDCIFile", &
                        "Memory allocation failure for zoneNumber")

       ! Loop over the grids and read the mapping.

       do i=1,nGrids

         ! Read the grid name for the current index.

         do
           read (readUnit,100) string
           if (string(1:10) == ":grid_name") exit
         end do
         gridName = trim(adjustl(string(11:)))

         ! Find which CGNS zone has the same name as this grid and
         ! store its index. Print an error if the grid name does not
         ! match any of the zone names.

         j = bsearchStrings(gridName, zoneNames, cgnsNDom)

 101     format ("Grid name ",a," not found in internal zone names")
         if (j == 0) then
           write(errorMessage,101) gridName
           call terminate("loadSuggarDCIFile", errorMessage)
         end if

         zoneNumber(i) = unsortedZone(j)

       end do

       ! Determine the global CGNS zones for which I have subblocks so
       ! that zones can be skipped for faster reading.

       readZone = .false.
       do i=1,nDom
         readZone(flowDoms(i,1,1)%cgnsBlockID) = .true.
       end do

       ! Now that the zones for this process and those affected by the
       ! file are known, determine if I need to continue reading the file
       ! at all by seeing if any match up.

       do i=1,nGrids
         if(readZone(zoneNumber(i))) exit
       end do

       if (i > nGrids) then
         close (readUnit)
         return
       end if

       ! Loop over the grids in the file and either skip it or read the
       ! donor-receptor pairs and iblanks depending on if I have a part
       ! of that zone.

       gridLoop: do igrid = 1,nGrids

         izone = zoneNumber(igrid)

         checkSkipGrid: if (readZone(izone)) then

           ! Determine the number of local blocks and their numbers that
           ! belong to izone.

           nBlocks = 0
           do i=1,nDom
             if (flowDoms(i,1,1)%cgnsBlockID == izone) then
               nBlocks = nBlocks + 1
               blockNum(nBlocks) = i
             end if
           end do

           ! Find the keyword for the donor-receptor pairs for this grid
           ! and read the number of pairs.

           do
             read (readUnit,100) string
             if (string(1:32) == ":general_cell_centered_drt_pairs") exit
           end do
           read (string(33:),*) nPairs

           ! Clear the overset boundary data and iblank array for the
           ! local blocks that are part of izone.

           clearLoop: do i=1,nBlocks
             j = blockNum(i)

             ! Deallocate the arrays to store the overset boundary and
             ! its donor information.

             deallocate(flowDoms(j,1,sps)%ibndry,         &
                        flowDoms(j,1,sps)%neighProcOver,  &
                        flowDoms(j,1,sps)%neighBlockOver, &
                        flowDoms(j,1,sps)%idonor,         &
                        flowDoms(j,1,sps)%overint, stat=ierr)
             if (ierr /= 0) &
               call terminate("loadSuggarDCIFile", &
                              "Deallocation failure for block data")

             ! Reset the number of overset cells and holes to 0.

             flowDoms(j,1,sps)%nCellsOverset    = 0
             flowDoms(j,1,sps)%nCellsOversetAll = 0
             flowDoms(j,1,sps)%nOrphans         = 0
             flowDoms(j,1,sps)%nHoles           = 0

             ! Reset all iblanks for this zone to 1 (all field).

             flowDoms(j,1,sps)%iblank = 1

             ! Allocate enough space assuming each block gets all pairs. 
             ! The sizes will be fixed after reading the pairs.

             allocate(flowDoms(j,1,sps)%ibndry(3,nPairs),       &
                      flowDoms(j,1,sps)%neighProcOver(nPairs),  &
                      flowDoms(j,1,sps)%neighBlockOver(nPairs), &
                      flowDoms(j,1,sps)%idonor(3,nPairs),       &
                      flowDoms(j,1,sps)%overint(8,nPairs), stat=ierr)
             if (ierr /= 0) &
               call terminate("loadSuggarDCIFile", &
                              "Memory allocation failure for block data")
           end do clearLoop

           ! Loop over the donor-receptor pairs.

           drPairLoop: do n = 1,nPairs

             ! Read a line for the donor-receptor pair. The line
             ! contains an index for the receptor cell followed by a
             ! donor grid number and an arbitrary number of donor
             ! indices and interpolants which is assumed to be 8 here.

             read (readUnit,*) myIndex, donorGrid, nDonor, &
                               donorIndices, donorInterp

             ! Convert the receptor index to i,j,k coordinates and loop
             ! over the local blocks for izone to determine which block
             ! (if any) contains the cell.

             call indexToIJK(izone, myIndex, i, j, k)

             do m = 1,nBlocks
               idom = blockNum(m)
               if (i >= flowDoms(idom,1,1)%iBegOr .and. &
                   i <  flowDoms(idom,1,1)%iEndOr .and. &
                   j >= flowDoms(idom,1,1)%jBegOr .and. &
                   j <  flowDoms(idom,1,1)%jEndOr .and. &
                   k >= flowDoms(idom,1,1)%kBegOr .and. &
                   k <  flowDoms(idom,1,1)%kEndOr) exit
             end do

             ! If the above loop was exhausted then none of the local
             ! blocks for izeon contain the cell so go to the next one.

             if (m > nBlocks) cycle

             ! Convert the zone indices to that of the local block.

             i = i - flowDoms(idom,1,1)%iBegOr + 2
             j = j - flowDoms(idom,1,1)%jBegOr + 2
             k = k - flowDoms(idom,1,1)%kBegOr + 2

             ! Update the number of overset cells for block idom and
             ! copy the indices to its boundary.

             ii = flowDoms(idom,1,sps)%nCellsOverset + 1
             flowDoms(idom,1,sps)%nCellsOverset = ii

             flowDoms(idom,1,sps)%ibndry(1,ii) = i
             flowDoms(idom,1,sps)%ibndry(2,ii) = j
             flowDoms(idom,1,sps)%ibndry(3,ii) = k

             ! Find the donor zone from the grid number and convert the
             ! donor index to i,j,k coordinates. Note we only need to
             ! use the first donor index which is the base or origin
             ! for the stencil.

             donorZone = zoneNumber(donorGrid)
             call indexToIJK(donorZone, donorIndices(1), i, j, k)

             ! Loop over the number of subblocks for the donor zone and
             ! determine which one conatains the donor cell.

             do m = 1,splitData(donorZone)%nSubBlocks
               if (i >= splitData(donorZone)%iBegOr(m) .and. &
                   i <  splitData(donorZone)%iEndOr(m) .and. &
                   j >= splitData(donorZone)%jBegOr(m) .and. &
                   j <  splitData(donorZone)%jEndOr(m) .and. &
                   k >= splitData(donorZone)%kBegOr(m) .and. &
                   k <  splitData(donorZone)%kEndOr(m)) exit
             end do

             ! Convert the donor zone indices to that of the local block.

             i = i - splitData(donorZone)%iBegOr(m) + 2
             j = j - splitData(donorZone)%jBegOr(m) + 2
             k = k - splitData(donorZone)%kBegOr(m) + 2

             ! Fill in the donor processor and local block, and copy the
             ! donor indices and interpolants to the boundary info.

             flowDoms(idom,1,sps)%neighProcOver(ii) = & 
                                      splitData(donorZone)%proc(m)
             flowDoms(idom,1,sps)%neighBlockOver(ii) = & 
                                      splitData(donorZone)%localBlock(m)

             flowDoms(idom,1,sps)%idonor(1,ii) = i
             flowDoms(idom,1,sps)%idonor(2,ii) = j
             flowDoms(idom,1,sps)%idonor(3,ii) = k

             flowDoms(idom,1,sps)%overint(:,ii) = donorInterp

           end do drPairLoop

           ! Reallocate the boundary and donor info arrays for each of
           ! the blocks of this zone since all pairs may not belong to
           ! one block. Note that ibndry is reallocated below after any
           ! orphans have been accounted for.

           do i = 1,nBlocks
             j  = blockNum(i)
             ii = flowDoms(j,1,sps)%nCellsOverset

             call reallocateInteger(flowDoms(j,1,sps)%neighProcOver, &
                                    ii, nPairs, .true.)
             call reallocateInteger(flowDoms(j,1,sps)%neighBlockOver, &
                                    ii, nPairs, .true.)
             call reallocateInteger2(flowDoms(j,1,sps)%idonor, &
                                     3_intType, ii, 3_intType, &
                                     nPairs, .true.)
             call reallocateReal2(flowDoms(j,1,sps)%overint, &
                                  8_intType, ii, 8_intType,  &
                                  nPairs, .true.)
           end do

           ! Find the keyword that starts the iblank section for this
           ! grid and read the number of iblank values in the file.

           do
             read (readUnit,100) string
             if (string(1:33) == ":general_cell_centered_drt_iblank") exit
           end do
           read (string(34:),*) nIblank

           ! Allocate memory to store an iblank array for the entire
           ! zone. The memory penalty is allowed here to avoid the if
           ! statement and loop to determine which indices belong to any 
           ! of my local blocks, which would be much slower.

           allocate(zoneBlank(cgnsDoms(izone)%nx, &
                              cgnsDoms(izone)%ny, &
                              cgnsDoms(izone)%nz), stat=ierr)
           if (ierr /= 0)                        &
             call terminate("loadSuggarDCIFile", &
                            "Memory allocation failure for zoneBlank")

           ! Read the entire iblank array for this zone. Ignore
           ! differences in iblank values for efficiency - this is
           ! handled below.

           do k = 1,cgnsDoms(izone)%nz
             do j = 1,cgnsDoms(izone)%ny
               do i =1,cgnsDoms(izone)%nx
                 read (readUnit,*) zoneBlank(i,j,k)
               end do
             end do
           end do

           ! Loop over the number of local blocks that belong to izone
           ! and copy the section of the zone's iblank array to that of
           ! the local block.

           blockLoop: do m = 1,nBlocks
             idom = blockNum(m)

             ! Set a pointer for this blocks iblank array and store
             ! the block dimensions easier.

             iblank => flowDoms(idom,1,sps)%iblank

             i = flowDoms(idom,1,1)%il
             j = flowDoms(idom,1,1)%jl
             k = flowDoms(idom,1,1)%kl

             ! Copy the section of the zone iblank array belonging to
             ! this local block.

             iblank(2:i,2:j,2:k) = zoneBlank( &
                flowDoms(idom,1,1)%iBegOr:flowDoms(idom,1,1)%iEndOr-1, &
                flowDoms(idom,1,1)%jBegOr:flowDoms(idom,1,1)%jEndOr-1, &
                flowDoms(idom,1,1)%kBegOr:flowDoms(idom,1,1)%kEndOr-1)

             ! Loop over the iblank array for this block and count the
             ! number of holes and orphans. Also, reset the iblank values
             ! to avoid any compatibility issues.

             do k = 2,flowDoms(idom,1,1)%kl
               do j = 2,flowDoms(idom,1,1)%jl
                 do i = 2,flowDoms(idom,1,1)%il

                   if (iblank(i,j,k) == 0) then

                     flowDoms(idom,1,sps)%nHoles = &
                     flowDoms(idom,1,sps)%nHoles + 1

                   else if (iblank(i,j,k) == 101) then

                     flowDoms(idom,1,sps)%nOrphans = &
                     flowDoms(idom,1,sps)%nOrphans + 1

                   else

                     iblank(i,j,k) = 1

                   end if

                 end do
               end do
             end do

             ! Add the number of orphans to the overset cells.

             n = flowDoms(idom,1,sps)%nCellsOverset
             flowDoms(idom,1,sps)%nCellsOversetAll = n + &
                                           flowDoms(idom,1,sps)%nOrphans

             ! Reallocate the boundary array for this block since
             ! memory was originally allocated for all pairs of izone.

             ii = flowDoms(idom,1,sps)%nCellsOversetAll

             call reallocateInteger2(flowDoms(idom,1,sps)%ibndry, &
                                     3_intType, ii, 3_intType,    &
                                     nPairs, .true.)

             ! If there were any orphans for this block then loop back
             ! over the iblank array to extract the indices.

             if (flowDoms(idom,1,sps)%nOrphans == 0) cycle

             do k = 2,flowDoms(idom,1,1)%kl
               do j = 2,flowDoms(idom,1,1)%jl
                 do i = 2,flowDoms(idom,1,1)%il
                   if (iblank(i,j,k) == 101) then

                     n = n + 1
                     flowDoms(idom,1,sps)%ibndry(1,n) = i
                     flowDoms(idom,1,sps)%ibndry(2,n) = j
                     flowDoms(idom,1,sps)%ibndry(3,n) = k

                   end if
                 end do
               end do
             end do

           end do blockLoop

           ! Deallocate the memory for zoneBlank.

           deallocate(zoneBlank, stat=ierr)
           if (ierr /= 0)                        &
             call terminate("loadSuggarDCIFile", &
                            "Deallocation failure for zoneBlank")

         else checkSkipGrid

           ! This processor may skip this grid. If this is the last grid 
           ! in the file then exit the loop. Otherwise, read the number 
           ! of pairs and iblanks and skip over each section.

           if (igrid == nGrids) exit

           do
             read (readUnit,100) string
             if (string(1:32) == ":general_cell_centered_drt_pairs") exit
           end do
           read (string(33:),*) nPairs

           do i=1,nPairs
             read (readUnit,*)
           end do

           do
             read (readUnit,100) string
             if (string(1:33) == ":general_cell_centered_drt_iblank") exit
           end do
           read (string(34:),*) nIblank

           do i=1,nIblank
             read (readUnit,*)
           end do

         end if checkSkipGrid
       end do gridLoop

       ! Close the DCI file and deallocate the memory stored for the
       ! grid to zone mapping.

       close (readUnit)

       deallocate(zoneNumber, stat=ierr)
       if (ierr /= 0)                        &
         call terminate("loadSuggarDCIFile", &
                        "Deallocation failure for zoneNumber")

!      ==================================================================

       contains

         subroutine indexToIJK(nn, ind, i, j, k)
!
!        ****************************************************************
!        *                                                              *
!        * indexToIJK converts the index ind for zone nn into its       *
!        * corresponding (i,j,k) indices.                               *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=intType), intent(in)  :: nn, ind
         integer(kind=intType), intent(out) :: i, j, k
!
!        Local variables.
!
         integer(kind=intType) :: ii, jj
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ii = ind - 1
         jj = cgnsDoms(nn)%nx*cgnsDoms(nn)%ny
         k  = ii/jj + 1
         jj = ii - (k - 1)*jj
         j  = jj/cgnsDoms(nn)%nx + 1
         i  = jj - (j - 1)*cgnsDoms(nn)%nx + 1

         end subroutine indexToIJK

       end subroutine loadSuggarDCIFile
