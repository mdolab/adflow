!
!      ******************************************************************
!      *                                                                *
!      * File:          distributeOversetCells.F90                      *
!      * Author:        Steve repsher                                   *
!      * Starting date: 12-30-2004                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine distributeOversetCells(cgnsId, nsubPerCGNS, &
                                         splitInfo)
!
!      ******************************************************************
!      *                                                                *
!      * DistributeOversetCells distributes the overset connectivity    *
!      * for an original cgns block cgnsID into its sublocks. It also   *
!      * determines the amount of communication from its donors.        *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use partitionMod
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in)    :: cgnsId

       integer(kind=intType), dimension(0:cgnsNDom), intent(in) :: &
                                                           nsubPerCGNS
       type(splitCGNSType), dimension(cgnsNDom), intent(in) :: splitInfo
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, k, l, m, ii, jj
       integer(kind=intType) :: ival, jval, kval, donorId

       integer(kind=intType), allocatable :: ncells(:)
       integer(kind=intType), allocatable :: cgnsOver(:,:)
       integer(kind=intType), allocatable :: ipntOver(:,:)
       integer(kind=intType), allocatable :: neighOver(:,:)

       character(len=maxCGNSNameLen) :: zonename, overName
       character(len=2*maxStringLen) :: errorMessage

       logical :: found
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the array that will count the number of overset cells
       ! in each of the sublocks and initialize to 0.

       allocate(ncells(splitInfo(cgnsId)%nsubblocks), stat=ierr)
       if(ierr /= 0) &
         call terminate("distributeOversetCells", &
                        "Memory allocation failure for ncells")
       ncells = 0

       ! Allocate enough memory to store the cells in each sublock
       ! by assumming all are in each sublock.

       ii = splitInfo(cgnsId)%nsubblocks
       jj = cgnsDoms(cgnsId)%ncellsOverset
       allocate(cgnsOver(jj,ii), ipntOver(jj,ii), neighOver(jj,ii), &
                stat=ierr)
       if(ierr /= 0) &
         call terminate("distributeOversetCells", &
                        "Memory allocation failure for cell data")

       ! Count the number of possible donor blocks and allocate
       ! storage for the amount of communication to this block.
       ! The donor block numbers are initialized to itself (this
       ! helps the graph partiioning setup ignore it if the
       ! number of cells remains 0).

       jj = 0
       do j = 1, cgnsDoms(cgnsId)%noverset
         donorID = cgnsDoms(cgnsId)%connOver(j)%donorBlock
         jj       = jj + splitInfo(donorId)%nsubblocks
       end do
       do k = 1, splitInfo(cgnsId)%nsubblocks
         ii = nsubPerCGNS(cgnsID - 1) + k
         allocate(blocks(ii)%overComm(jj,2), stat=ierr)
         if(ierr /= 0) &
           call terminate("distributeOversetCells", &
                        "Memory allocation failure for overComm")
         blocks(ii)%overComm(:,1) = ii
         blocks(ii)%overComm(:,2) = 0
       end do

       ! Loop over the overset boundaries of the original block.

       m = 0
       noverLoop: do j = 1, cgnsDoms(cgnsId)%noverset

         ! Store the cgns donor block for this overset connectivity.

         donorID = cgnsDoms(cgnsId)%connOver(j)%donorBlock

         ! Loop over the cells to be distributed.

         cellLoop: do i = 1, cgnsDoms(cgnsId)%connOver(j)%npnts

           ! Store the indices easier.

           ival = cgnsDoms(cgnsId)%connOver(j)%ibndry(1,i)
           jval = cgnsDoms(cgnsId)%connOver(j)%ibndry(2,i)
           kval = cgnsDoms(cgnsId)%connOver(j)%ibndry(3,i)

           ! Loop over the sublocks to find which one this cell is in.

           found = .false.
           sublockLoop: do k = 1, splitInfo(cgnsId)%nsubblocks

             ! Find the compute block number for this sublock.

             ii = nsubPerCGNS(cgnsID - 1) + k

             ! Check if the cell is in this sublock.

             if(ival <  blocks(ii)%iEndor .and. &
                ival >= blocks(ii)%iBegor .and. &
                jval <  blocks(ii)%jEndor .and. &
                jval >= blocks(ii)%jBegor .and. &
                kval <  blocks(ii)%kEndor .and. &
                kval >= blocks(ii)%kBegor) then

               ! Update the number of overset cells for this sublock and
               ! store the info for this cell in the temporary set.

               ncells(k)              = ncells(k) + 1
               cgnsOver(ncells(k),k) = j
               ipntOver(ncells(k),k) = i

               ! Cell has been found so break out of the loop.

               found = .true.
               exit

             end if
           end do sublockLoop

           ! Cell has not been found in any sublock so something is
           ! wrong. Processor 0 prints an error message, while the
           ! others wait until they are killed.

           if(.not. found) then
             if(myID == 0) then
               zonename  = cgnsDoms(cgnsId)%zonename
               overName = cgnsDoms(cgnsId)%connOver(j)%connectName
               write(errorMessage,140) trim(zonename), trim(overName), i
 140           format("Zone",1X,A,", overset block connectivity",1X,A, &
                      ": Cell",1X,I5,1X,"could not be found.")
               call terminate("distributeOversetCells", errorMessage)
             endif

             call mpi_barrier(SUmb_comm_world, ierr)
           endif

           ! Store the indices easier and repeat the
           ! process to find the donor compute block.

           ival = cgnsDoms(cgnsId)%connOver(j)%idonor(1,i)
           jval = cgnsDoms(cgnsId)%connOver(j)%idonor(2,i)
           kval = cgnsDoms(cgnsId)%connOver(j)%idonor(3,i)

           ! Loop over the sublocks to find which one this cell is in.

           found = .false.
           donorSublockLoop: do l = 1, splitInfo(donorId)%nsubblocks

             ! Find the compute block number for this sublock.

             jj = nsubPerCGNS(donorID - 1) + l

             ! Check if the cell is in this sublock.

             if(ival <  blocks(jj)%iEndor .and. &
                ival >= blocks(jj)%iBegor .and. &
                jval <  blocks(jj)%jEndor .and. &
                jval >= blocks(jj)%jBegor .and. &
                kval <  blocks(jj)%kEndor .and. &
                kval >= blocks(jj)%kBegor) then

               ! Store the donor compute block id and increment the
               ! amount of communication.

               neighOver(ncells(k),k)     = jj
               blocks(ii)%overComm(m+l,1) = jj
               blocks(ii)%overComm(m+l,2) = blocks(ii)%overComm(m+l,2) + 1

               ! Cell has been found so break out of the loop.

               found = .true.
               exit

             end if
           end do donorSublockLoop

           ! Cell has not been found in any sublock so something is
           ! wrong. Processor 0 prints an error message, while the
           ! others wait until they are killed.

           if(.not. found) then
             if(myID == 0) then
               zonename  = cgnsDoms(cgnsId)%zonename
               overName = cgnsDoms(cgnsId)%connOver(j)%connectName
               write(errorMessage,150) trim(zonename), trim(overName), i
 150           format("Zone",1X,A,", overset block connectivity",1X,A, &
                      ": Donor",1X,I5,1X,"could not be found.")
               call terminate("distributeOversetCells", errorMessage)
             endif

             call mpi_barrier(SUmb_comm_world, ierr)
           endif

         end do cellLoop

         ! Update the index m for the next overset connectivity.

         m = m + splitInfo(donorId)%nsubblocks

       end do noverLoop

       ! Loop over the sublocks and copy the number of cells.  Allocate
       ! the correct amount of memory and copy the overset data.

       copyDataLoop: do k = 1, splitInfo(cgnsId)%nsubblocks

         ii = nsubPerCGNS(cgnsID - 1) + k

         allocate(blocks(ii)%cgnsOver(ncells(k)), &
                  blocks(ii)%ipntOver(ncells(k)),      &
                  blocks(ii)%neighOver(ncells(k)), stat=ierr)
         if(ierr /= 0) &
           call terminate("distributeOversetCells", &
                          "Memory allocation failure for block data")

         blocks(ii)%ncellsOverset = ncells(k)
         blocks(ii)%cgnsOver      = cgnsOver(1:ncells(k),k)
         blocks(ii)%ipntOver      = ipntOver(1:ncells(k),k)
         blocks(ii)%neighOver     = neighOver(1:ncells(k),k)

       end do copyDataLoop

       ! Release the memory for the temporary data.

       deallocate(ncells, cgnsOver, ipntOver, neighOver, &
                  stat=ierr)
       if(ierr /= 0) &
         call terminate("distributeOversetCells", &
                        "Memory release failure for cell data")

       end subroutine distributeOversetCells
