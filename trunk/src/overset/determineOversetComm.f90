!
!      ******************************************************************
!      *                                                                *
!      * File:          determineOversetComm.f90                        *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 04-22-2005                                      *
!      * Last modified: 10-16-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineOversetComm(level, sps, coarseLevel)
!
!      ******************************************************************
!      *                                                                *
!      * determineOversetComm determines the overset communication data *
!      * structures for the given level and spectral mode using the     *
!      * boundary list. If this is a coarse level or the input overset  *
!      * donors are being treated as guesses, then a donor search is    *
!      * also performed. At the end of this routine is a complete set   *
!      * of iblanks and boundary info for each block along with the     *
!      * communication data structures.                                 *
!      *                                                                *
!      ******************************************************************
!
       use block
       use boundaryList
       use cgnsGrid
       use communication
       use inputOverset
       use searchMod
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: level, sps
       logical,               intent(in) :: coarseLevel
!
!      Local variables.
!
       integer               :: ierr, size
       integer(kind=intType) :: b, i, j, k, m, n, l, fineLevel

       integer(kind=intType), dimension(2) :: nResFail, nResFailGlobal

       integer(kind=intType), dimension(nDom) :: cc

       integer(kind=intType), dimension(nProc) :: nOrphansPerProc
       integer,               dimension(nProc) :: recvcounts, displs

       integer(kind=intType), allocatable :: orphanGlobal(:,:)
       integer(kind=intType), allocatable :: orphan(:,:)

       type(commType)         :: commRestarts
       type(internalCommType) :: internalRestarts
!
!      Function definitions.
!
       integer(kind=intType) :: numAdjacentField
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Sort the boundary list in increasing order.

       call qsortHaloListType(oversetHalo, nHaloOver)

       ! Determine the communication structures from the sorted list.
       ! Note the size of interp is 3 if this is a coarse level or the
       ! fine inputs were guesses because the coordinates are stored.

       i = nInterp
       if (coarseLevel .or. oversetDonorsAreGuesses) i = 3

       call finalCommStructures(oversetHalo, nHaloOver,        &
                                commPatternOverset(level,sps), &
                                internalOverset(level,sps), i)

       ! Nullify the periodic stuff for the communication.

       commPatternOverset(level,sps)%nPeriodic = 0
       commRestarts%nPeriodic                  = 0
       internalOverset(level,sps)%nPeriodic    = 0
       internalRestarts%nPeriodic              = 0

       nullify(commPatternOverset(level,sps)%periodicData)
       nullify(commRestarts%periodicData)
       nullify(internalOverset(level,sps)%periodicData)
       nullify(internalRestarts%periodicData)

       ! The coarse boundary list built contained donor index estimates
       ! based on the next finer level; therefore, they need to be
       ! brought up to this level. Loop over the communcation lists
       ! and do this. Note there is a check if the donor indices are
       ! equal to 1 in which case the origin is not an owned cell and
       ! is not within the bound of the coarsening variables.

       fixDonorIndices: if (coarseLevel) then

         fineLevel = level - 1

         do n = 1,commPatternOverset(level,sps)%nProcSend
           do m = 1,commPatternOverset(level,sps)%nSend(n)

             b = commPatternOverset(level,sps)%sendList(n)%block(m)
             i = commPatternOverset(level,sps)%sendList(n)%indices(m,1)
             j = commPatternOverset(level,sps)%sendList(n)%indices(m,2)
             k = commPatternOverset(level,sps)%sendList(n)%indices(m,3)

             if (i > 1) &
             commPatternOverset(level,sps)%sendList(n)%indices(m,1) = &
                                flowDoms(b,fineLevel,sps)%mgICoarse(i,1)
             if (j > 1) &
             commPatternOverset(level,sps)%sendList(n)%indices(m,2) = &
                                flowDoms(b,fineLevel,sps)%mgJCoarse(j,1)
             if (k > 1) &
             commPatternOverset(level,sps)%sendList(n)%indices(m,3) = &
                                flowDoms(b,fineLevel,sps)%mgKCoarse(k,1)
           end do
         end do

         do n = 1,internalOverset(level,sps)%ncopy

           b = internalOverset(level,sps)%donorBlock(n)
           i = internalOverset(level,sps)%donorIndices(n,1)
           j = internalOverset(level,sps)%donorIndices(n,2)
           k = internalOverset(level,sps)%donorIndices(n,3)

           if (i > 1) &
           internalOverset(level,sps)%donorIndices(n,1) = &
                               flowDoms(b,fineLevel,sps)%mgICoarse(i,1)
           if (j > 1) &
           internalOverset(level,sps)%donorIndices(n,2) = &
                               flowDoms(b,fineLevel,sps)%mgJCoarse(j,1)
           if (k > 1) &
           internalOverset(level,sps)%donorIndices(n,3) = &
                               flowDoms(b,fineLevel,sps)%mgKCoarse(k,1)
         end do

       end if fixDonorIndices

       ! If this is a coarse level or the external inputs are being 
       ! ignored for the fine level, then perform the steps to find
       ! the donor stencils and interpolants.

       performDonorSearch: if (coarseLevel .or. &
                               oversetDonorsAreGuesses) then

         ! Tell the search to ignore the donor quality cutoff input such
         ! that as many "perfect" donors are found as possible.

         ignoreCutoffQuality = .true.

         ! Initializate the number of donor info sets to 0 and allocate
         ! the array with 0 size. This will be changed during a donor
         ! search cycle if needed.

         nDonorInfo = 0
         allocate(donorInfo(5,nDonorInfo), stat=ierr)
         if (ierr /= 0)                           &
           call terminate("determineOversetComm", &
                          "Memory allocation failure for donorInfo")

         ! Initialize the number of failures (orphans). This is the total
         ! number of halos sorted above minus those with donors that were
         ! added to the communication structures.

         i = commPatternOverset(level,sps)%nProcRecv
         nResFail(2) = nHaloOver - internalOverset(level,sps)%nCopy &
                     - commPatternOverset(level,sps)%nRecvCum(i)

         ! Execute a donor search cycle with the main list.

         call donorSearchCycle(level, sps, nResFail(2),       &
                               commPatternOverset(level,sps), &
                               internalOverset(level,sps))

         ! Print a header in debug mode.

 101     format ("# Searching for donors on level ",i2,", spectral ",i2)
 102     format ("#",3x,i5,4x,i8,4x,i7)
 103     format ("#",3x,i5,4x,i8,4x,i7,4x,a)

         if (debug .and. myID == 0) then
           print "(a)", "#"
           print  101 , level, sps
           print "(a)", "#   Cycle    Restarts    Orphans"
         end if
 
         ! Initialize the cycle counter and begin iterating.

         n = 1
         restartLoop: do

           ! Compute the local number of restarts. This is the difference
           ! between the number of halos left without donors and the
           ! number of failures.

           nResFail(1) = nHaloOver - nResFail(2)

           ! Find the global number of restarts for every process and
           ! print it along with the number of orphans in debug mode.

           call mpi_allreduce(nResFail, nResFailGlobal, 2, sumb_integer, &
                              mpi_sum, SUmb_comm_world, ierr)

           if (debug .and. myID == 0) then
             print 102, n, nResFailGlobal
           end if

           ! Exchange the iblanks on the 1-to-1 boundaries since there
           ! may be new field or "restarts" as a result of fringe changes
           ! due to poor donor quality. This also refreshes the iblanks
           ! on the halos so a scan for new fringe doesn't find additions
           ! from a previous cycle.

           call exchangeIblanks(level, sps, commPatternCell_2nd, &
                                internalCell_2nd)

           ! If the total sum of restarts is 0, then possibly reduce the
           ! number of orphans by using the donor quality input.

           checkApplyQuality: if (nResFailGlobal(1) == 0) then

             ! If the quality cutoff has alreay been applied then exit.

             if (.not. ignoreCutoffQuality) exit

             ! Use the quality for the next search cycle. Loop over the
             ! halo list and look for orphans that were declared as bad
             ! donors. This is indicated by the level of indirectness set
             ! to something > 0 becasue it is actually the donor
             ! processor. Also, require that the donor block be > 0. If
             ! the donor block and indices are 0 then the cell is a 
             ! permanent orphan. Just swap the value and update the local number
             ! of orphans.

             ignoreCutoffQuality = .false.

             do i = 1,nHaloOver
               if (oversetHalo(i)%levOfInd  >= 0 .and. &
                   oversetHalo(i)%donorBlock > 0) then
                 oversetHalo(i)%donorProc = oversetHalo(i)%levOfInd
                 oversetHalo(i)%levOfInd  = 0
                 nResFail(2)              = nResFail(2) - 1
               end if
             end do

             ! Recalculate the local and global numbers of restarts and
             ! orphans. Exit if still no restarts.

             nResFail(1) = nHaloOver - nResFail(2)

             call mpi_allreduce(nResFail, nResFailGlobal, 2, &
                                sumb_integer, mpi_sum,       &
                                SUmb_comm_world, ierr)

             if (nResFailGlobal(1) == 0) exit

             ! Print another line for this cycle in debug mode.

             if (debug .and. myID == 0) then
               print 103, n, nResFailGlobal, "Appyling input qulity"
             end if

           end if checkApplyQuality

           ! Sort the boundary list in increasing order.

           call qsortHaloListType(oversetHalo, nHaloOver)

           ! Determine the communication structures for the restarts.

           call finalCommStructures(oversetHalo, nHaloOver,         &
                                    commRestarts, internalRestarts, &
                                    3_intType)

           ! Reset the number of donor info sets stored to 0 since those
           ! from the previous cycle are not needed.

           nDonorInfo = 0

           ! Execute a donor search cycle with the restart lists.

           call donorSearchCycle(level, sps, nResFail(2), &
                                 commRestarts, internalRestarts)

           ! Merge the successes of the restart list into the main list.
           ! Note that all failures were removed so the entire list can
           ! be merged. Then release the memory of the restart list.

           call mergeComm(commPatternOverset(level,sps), &
                          internalOverset(level,sps),    &
                          commRestarts, internalRestarts)

           call releaseCommPattern(commRestarts)
           call releaseInternalComm(internalRestarts)

           ! Update the number of cycles.

           n = n + 1

         end do restartLoop

         ! Release the memory for the donor info sets.

         deallocate(donorInfo, stat=ierr)
         if (ierr /= 0)                           &
           call terminate("determineOversetComm", &
                          "Deallocation failure for donorInfo")

         ! Initialize the number of orphans for each local block to 0.
         ! Then loop over the list and count the orphans per block.

         do n=1,nDom
           flowDoms(n,level,sps)%nOrphans = 0
         end do

         do i=1,nHaloOver
           b = oversetHalo(i)%myBlock
           flowDoms(b,level,sps)%nOrphans = &
                                      flowDoms(b,level,sps)%nOrphans + 1
         end do

         ! Subtract the number of orphans from nCellsOverset so that a
         ! distinction can be made between boundary cells with donors
         ! and those which are orphans.

         do n=1,nDom
           flowDoms(n,level,sps)%nCellsOverset = & 
                                 flowDoms(n,level,sps)%nCellsOverset - &
                                 flowDoms(n,level,sps)%nOrphans
         end do

         ! Use the final communication pattern to form the arrays of the
         ! flow domains which are needed to process the next level and
         ! for output.

         call setBlockOversetData(level, sps)

         ! Now that the flow domain data has been set, change the
         ! iblanks on the boundary to their appropriate values.

         do n=1,nDom
           call setPointers(n, level, sps)
           call changeIblanks(.false., 10_intType)
         end do

       else performDonorSearch

         ! Since a donor search isn't performed, modify nHaloOver to be
         ! the number of cells in the list which still need to be dealt
         ! with (i.e. the total number of orphans).

         nHaloOver = 0
         do n=1,nDom
           nHaloOver = nHaloOver + flowDoms(n,level,sps)%nOrphans
         end do

       end if performDonorSearch

       ! Sort the orphans by the number of adjacent field cells (or
       ! fringe cells with donors), in descending order. This is done by
       ! putting that number in the donorProc variable of the list
       ! because this is the primary sort variable for a halo list.

       do i=1,nHaloOver
         oversetHalo(i)%donorProc = -numAdjacentField(level, sps, &
                                     oversetHalo(i)%myBlock,      &
                                     oversetHalo(i)%myI,          &
                                     oversetHalo(i)%myJ,          &
                                     oversetHalo(i)%myK)
       end do

       call qsortHaloListType(oversetHalo, nHaloOver)

       ! Now extract the orphan indices to the blocks in the sorted
       ! order. In this way, when orphans are averaged, the ones with
       ! the most adjacent field cells are done first. This is so that
       ! when orphans are averaged using other orphans, things are done
       ! as efficiently as possible.

       do n=1,nDom
         cc(n) = flowDoms(n,level,sps)%nCellsOverset
       end do

       do i=1,nHaloOver
         n = oversetHalo(i)%myBlock
         cc(n) = cc(n) + 1

         flowDoms(n,level,sps)%ibndry(1,cc(n)) = oversetHalo(i)%myI
         flowDoms(n,level,sps)%ibndry(2,cc(n)) = oversetHalo(i)%myJ
         flowDoms(n,level,sps)%ibndry(3,cc(n)) = oversetHalo(i)%myK
       end do

       ! Determine the number of orphans per processor, and reset the
       ! global number in case a search wasn't performed.

       call mpi_allgather(nHaloOver, 1, sumb_integer, nOrphansPerProc, &
                          1, sumb_integer, SUmb_comm_world, ierr)

       nResFailGlobal(2) = sum(nOrphansPerProc)

       ! If the global number of orphans is > 0 then gather the
       ! orphans from all processors and print a fairly detailed report.

       printOrphanReport: if (nResFailGlobal(2) > 0) then

         ! Allocate the memory to store all my orphan info.

         allocate(orphan(8,nHaloOver), stat=ierr)
         if (ierr /= 0)                           &
           call terminate("determineOversetComm", &
                          "Memory allocation failure for orphan")
 
         ! Extract the orphans from the boundary list, and convert the
         ! block and indices to the original cgns blocks. The indices
         ! are converted to the "closest fine level" indices so that
         ! the user can roughly locate the problem area on the original
         ! blocks. Note this step is omitted for the terminal donor
         ! info because that was taken care of by the donor search
         ! cycles.

         do n = 1,nHaloOver

           ! Store the cgns block id.

           b           = oversetHalo(n)%myBlock
           orphan(1,n) = flowDoms(b,1,1)%cgnsBlockId

           orphan(5,n) = oversetHalo(n)%donorBlock

           ! Store the indices.

           i = oversetHalo(n)%myI
           j = oversetHalo(n)%myJ
           k = oversetHalo(n)%myK

           orphan(6,n) = oversetHalo(n)%dI
           orphan(7,n) = oversetHalo(n)%dJ
           orphan(8,n) = oversetHalo(n)%dK
 
           ! Take the indices to the finest level. Note the subscript
           ! 1 is used in the restriction variable so that the final
           ! result will be the origin of the cells that the coarse
           ! orphan contained.

           if (i == flowDoms(b,level,1)%ib) then
             i = flowDoms(b,1,1)%ib
           else if (i == flowDoms(b,level,1)%ie) then
             i = flowDoms(b,1,1)%ie
           else if (i <= 1) then
           else
             do l = level,2,-1
               i = flowDoms(b,l,1)%mgIFine(i,1)
             end do
           end if

           if (j == flowDoms(b,level,1)%jb) then
             j = flowDoms(b,1,1)%jb
           else if (j == flowDoms(b,level,1)%je) then
             j = flowDoms(b,1,1)%je
           else if (j <= 1) then
           else
             do l = level,2,-1
               j = flowDoms(b,l,1)%mgJFine(j,1)
             end do
           end if

           if (k == flowDoms(b,level,1)%kb) then
             k = flowDoms(b,1,1)%kb
           else if (k == flowDoms(b,level,1)%ke) then
             k = flowDoms(b,1,1)%ke
           else if (k <= 1) then
           else
             do l = level,2,-1
               k = flowDoms(b,l,1)%mgKFine(k,1)
             end do
           end if

           ! Finally, convert the cell indices to the cgns index style.

           orphan(2,n) = i + flowDoms(b,1,1)%iBegor - 2
           orphan(3,n) = j + flowDoms(b,1,1)%jBegor - 2
           orphan(4,n) = k + flowDoms(b,1,1)%kBegor - 2

         end do

         ! Processor 0 needs to do some extra work.

         if (myID == 0) then

           ! Determine arrays recvcounts and displs for the call to 
           ! gatherv.

           recvcounts(1) = 8*nOrphansPerProc(1)
           displs(1)     = 0

           do n=2,nProc
             recvcounts(n) = 8*nOrphansPerProc(n)
             displs(n)     = displs(n-1) + recvcounts(n-1)
           enddo

           ! Allocate memory to store all of the orphans from every
           ! processor.

           allocate(orphanGlobal(8,nResFailGlobal(2)), stat=ierr)
           if (ierr /= 0)                           &
             call terminate("determineOversetComm", &
                         "Memory allocation failure for orphanGlobal")

         end if

         ! Gather the data to process 0.

         size = 8*nHaloOver
         call mpi_gatherv(orphan, size, sumb_integer, orphanGlobal, &
                          recvcounts, displs, sumb_integer, 0,      &
                          SUmb_comm_world, ierr)

         ! Processor 0 prints a header and the gathered info.

 104     format ("# Orphan Report for Level ",i2, &
                 ", Spectral Mode ",i2," (Total = ",i6,")")
 105     format ("# ",a18,3(1x,i4),1x,a18,3(1x,i4))
 106     format ("# ",a18,3(1x,i4),1x,"No donor info available")
 107     format ("# ",a18,1x,a14,1x,a18,1x,a14)

         if (myID == 0) then

           print "(a)", "#"
           print "(a)", "#*==================== !!! Warning !!! &
                        &=============================="
           print  104 , level, sps, nResFailGlobal(2)
           print "(a)", "#"

           print  107 , "Zone              ", &
                        "Cell Indices  ",     &
                        "Last Donor Zone   ", &
                        "Last Indices  "
           print  107 , "------------------", &
                        "--------------",     &
                        "------------------", &
                        "--------------"

           do i = 1,nResFailGlobal(2)
             if (orphanGlobal(5,i) == 0) then
               print 106,                                              &
                  trim(adjustl(cgnsDoms(orphanGlobal(1,i))%zonename)), &
                  orphanGlobal(2:4,i)
             else
               print 105,                                              &
                  trim(adjustl(cgnsDoms(orphanGlobal(1,i))%zonename)), &
                  orphanGlobal(2:4,i),                                 &
                  trim(adjustl(cgnsDoms(orphanGlobal(5,i))%zonename)), &
                  orphanGlobal(6:8,i)
             end if
           end do

           print "(a)", "#"
           print "(a)", "#*=====================================&
                        &=============================="
           print "(a)", "#"

         end if
       end if printOrphanReport

       ! Deallocate the memory for the interpolants of the halo list.

       do i = 1,ubound(oversetHalo,1)
         deallocate(oversetHalo(i)%interp, stat=ierr)
         if(ierr /= 0)                           &
          call terminate("determineOversetComm", &
                         "Deallocation failure for interp")
       enddo

       ! Deallocate the memory for the halo list.

       deallocate(oversetHalo, stat=ierr)
       if(ierr /= 0)                           &
        call terminate("determineOversetComm", &
                       "Deallocation error for oversetHalo")
 
       end subroutine determineOversetComm
