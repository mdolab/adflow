!
!      ******************************************************************
!      *                                                                *
!      * File:          checkOverset.f90                                *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 03-12-2005                                      *
!      * Last modified: 07-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine checkOverset(level, sps, checkOverlap, checkField, &
                               fringeSize)
!
!      ******************************************************************
!      *                                                                *
!      * CheckOverset checks the integrity of the overset connectivity  *
!      * and holes on the given level. There are 2 checks available to  *
!      * be performed, which are controlled by the logical arguments.   *
!      *                                                                *
!      * (1) overlap check:                                             *
!      * This check will essentially look at every block for cells that *
!      * are both donors and receivers and thus there is a breakdown in *
!      * the communication, in that region which would create non-      *
!      * physical results from the incomplete boundary.                 *
!      *                                                                *
!      * (2) field cell check:                                          *
!      * This check will identify situations where the residual of any  *
!      * field cell will be invalid because its stencil has a hole.     *
!      *                                                                *
!      * The results of all processors are gathered and if even 1 bad   *
!      * entity is found an error message is printed and the program    *
!      * will terminate.                                                *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsGrid
       use communication
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: level, sps, fringeSize
       logical,               intent(in) :: checkField, checkOverlap
!
!      Local variables.
!
       integer               :: ierr, size
       integer(kind=intType) :: i, j, k, l, n, m, ibval, icomO, icomF
       integer(kind=intType) :: ncom, nbad, nbadGlobal, b, i1, j1, k1
       integer(kind=intType) :: magic

       integer(kind=intType), dimension(nProc)  :: counts
       integer,               dimension(nProc)  :: recvcounts, displs
       integer(kind=intType), dimension(nDom)   :: badIndex
       integer(kind=intType), dimension(3)      :: total

       integer(kind=intType), allocatable :: badGlobal(:,:), bad(:,:)

       character(len=8) :: intString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine how many variables will be gathered for each domain
       ! and the indices to store them in the bad arrays.
 
       ncom  = 1
       icomO = 0
       icomF = 0
       if (checkOverlap) then
         ncom = ncom + 1
         icomO = ncom
       end if
       if (checkField) then
         ncom = ncom + 1
         icomF = ncom
       end if
 
       ! Allocate memory for the bad array counter for this processor.
 
       allocate(bad(ncom,nDom), stat=ierr)
       if(ierr /= 0)                    &
         call terminate("checkOverset", &
                        "Memory allocation failure for bad")
 
       ! Initialize the bad block and cell counters and index to 0.
 
       nbad     = 0
       bad      = 0
       badIndex = 0
 
!      ==================================================================
!      | (1) - Check the overlap of local domains                       |
!      ==================================================================
 
       overlapCheck: if (checkOverlap) then

        ! Loop over the overset send list for this processor - these are
        ! the donor cells.
 
        do n = 1,commPatternOverset(level,sps)%nProcSend
          do m = 1,commPatternOverset(level,sps)%nsend(n)

            ! Store the block and indices for the stencil corner point.
 
            b = commPatternOverset(level,sps)%sendList(n)%block(m)
            i = commPatternOverset(level,sps)%sendList(n)%indices(m,1)
            j = commPatternOverset(level,sps)%sendList(n)%indices(m,2)
            k = commPatternOverset(level,sps)%sendList(n)%indices(m,3)
 
            ! Look at the iblank values for the 8 cell stencil. To be valid,
            ! all 8 must be part of the field. There are situations where
            ! the donor stencil may include cells on the 1st halo level.
            ! One is if blocks were split, so it is assummed that these
            ! iblanks have already been exchanged. The other is when a
            ! cell is adjacent to a boco (not oversetOuterBound!),
            ! where the iblank is 2. The math below takes this into
            ! account
 
            ibval = maxval(abs(flowDoms(b,level,sps)%iblank(i:i+1, &
                                                            j:j+1, &
                                                            k:k+1) - 2))
 
            if (ibval > 1) then
              if (badIndex(b) == 0) then
                nbad = nbad + 1
                badIndex(b) = nbad
                bad(1,nbad) = flowDoms(b,level,1)%cgnsBlockId
              end if
              bad(icomO,badIndex(b)) = bad(icomO,badIndex(b)) + 1
            end if
 
          end do
        end do

        ! Loop over the internal overset list and repeat the process.
 
        do n = 1,internalOverset(level,sps)%ncopy
 
          b = internalOverset(level,sps)%donorBlock(n)
          i = internalOverset(level,sps)%donorIndices(n,1)
          j = internalOverset(level,sps)%donorIndices(n,2)
          k = internalOverset(level,sps)%donorIndices(n,3)
 
          ibval = maxval(abs(flowDoms(b,level,sps)%iblank(i:i+1, &
                                                         j:j+1, &
                                                         k:k+1) - 2))
 
          if (ibval > 1) then
            if (badIndex(b) == 0) then
              nbad = nbad + 1
              badIndex(b) = nbad
              bad(1,nbad) = flowDoms(b,level,1)%cgnsBlockId
            end if
            bad(icomO,badIndex(b)) = bad(icomO,badIndex(b)) + 1
          end if
 
        end do
 
       end if overlapCheck

!      ==================================================================
!      | (2) - Check stencils of field cells on local domains           |
!      ==================================================================
 
       fieldCheck: if (checkField) then

        ! Compute the magic number that should be the total from the 
        ! sum used below. Basically this is the number of cells which
        ! contribute to a cell's residual including itself.

        magic = 27 + 6*(fringeSize - 1)
 
        ! Loop over the domains.
 
        m = 1

        domains: do b = 1,nDom
 
          ! Initialize the counters for bad cells to 0.
 
          n = 0
 
          ! Loop over the iblank array and for each field cell do a
          ! special sum such that each stencil point should contribute
          ! a value of 1 unless it is a hole, which contributes 0. An
          ! oversetOuterBoundary may also contribute -1.
 
          do k = 2,flowDoms(b,level,1)%kl
            do j = 2,flowDoms(b,level,1)%jl
              do i = 2,flowDoms(b,level,1)%il
                if (flowDoms(b,level,sps)%iblank(i,j,k) == 1) then
 
                  ibval = 0
 
                  ! First layer of cells around this one (= 27
                  ! including itself)
 
                  do k1 = k-1,k+1
                   do j1 = j-1,j+1
                    do i1 = i-1,i+1
                     ibval = ibval &
                           + min(m, flowDoms(b,level,sps)%iblank(i1,j1,k1))
                    end do
                   end do
                  end do
 
                  ! Now the 2nd and higher level cells in the coordinate
                  ! directions (= 6*(fringeSize-1) more cells)

                  do l = 2,fringeSize 
                    ibval = ibval &
                          + min(m, flowDoms(b,level,sps)%iblank(i+l,j,k)) &
                          + min(m, flowDoms(b,level,sps)%iblank(i-l,j,k)) &
                          + min(m, flowDoms(b,level,sps)%iblank(i,j+l,k)) &
                          + min(m, flowDoms(b,level,sps)%iblank(i,j-l,k)) &
                          + min(m, flowDoms(b,level,sps)%iblank(i,j,k+l)) &
                          + min(m, flowDoms(b,level,sps)%iblank(i,j,k-l))
                  end do
 
                  ! The total should be 33, otherwise update the counter.
 
                  if (ibval /= magic) n = n + 1
 
                end if
              end do
            end do
          end do
 
          ! If the counter is non-zero, then update the counter nbad
          ! and store the number of bad cells.
 
          if (n > 0) then
            if (badIndex(b) == 0) then
              nbad = nbad + 1
              badIndex(b) = nbad
              bad(1,nbad) = flowDoms(b,level,1)%cgnsBlockId
            end if
            bad(icomF,badIndex(b)) = n
          end if
 
        end do domains
 
       end if fieldCheck
 
!      ==================================================================
!      | Gather results from all processors and write errors.           |
!      ==================================================================
 
       ! Determine the number of bad blocks per processor.

       call mpi_allgather(nbad, 1, sumb_integer, counts, 1, &
                          sumb_integer, SUmb_comm_world, ierr)

       ! Determine the global number of bad blocks and the arrays
       ! recvcounts and displs needed for the call to allgatherv.

       nbadGlobal    = counts(1)
       recvcounts(1) = ncom*counts(1)
       displs(1)     = 0

       do n=2,nProc
         nbadGlobal    = nbadGlobal + counts(n)
         recvcounts(n) = ncom*counts(n)
         displs(n)     = displs(n-1) + recvcounts(n-1)
       enddo

       ! Allocate the memory to store all the bad blocks.

       allocate(badGlobal(ncom,nbadGlobal), stat=ierr)
       if(ierr /= 0)                   &
         call terminate("checkOverset", &
                        "Memory allocation failure for badGlobal")

       ! Gather the data.

       size = ncom*nbad
       call mpi_allgatherv(bad, size, sumb_integer, badGlobal, &
                           recvcounts, displs, sumb_integer,   &
                           SUmb_comm_world, ierr)
 
       ! If bad blocks are present, print them to stdout and terminate.

       testBadPresent: if(nbadGlobal > 0) then

         ! The data is only written by processor 0.

         testRootProc: if(myID == 0) then

           ! Write a header.

           print "(a)", "#"
           print  101,  level, sps
           print "(a)", "# and/or connectivity information is invalid."
           print "(a)", "# Here is the list of bad blocks..."
           print "(a)", "#"

           ! Write the bad cgns blocks by adding up the contributions
           ! from each processor. Not really going for efficiency here.

           do i = 1,cgnsNDom

             ! Reset the counters and total up the bad cells for this
             ! cgns domain.
 
             total = 0
 
             do j = 1,nbadGlobal
               if (badGlobal(1,j) == i) then
                 do k = 2,ncom
                   total(k) = total(k) + badGlobal(k,j)
                 end do
               end if
             end do
 
             ! Write the error message if each counter is non-zero.
 
             if (icomF > 0 .and. total(icomF) > 0) then
               write(intString,"(i7)") total(icomF)
               intString = adjustl(intString)
               print 102, trim(cgnsDoms(i)%zonename), &
                          trim(intString)
             end if
 
             if (icomO > 0 .and. total(icomO) > 0) then
               write(intString,"(i7)") total(icomO)
               intString = adjustl(intString)
               print 103, trim(cgnsDoms(i)%zonename), &
                          trim(intString)
             end if
 
           enddo

 101       format("# Found blocks on level ",i1,", spectral ",i2, &
                  " for which the overset holes")
 102       format("# Zone ",a, ": ",a," field cell stencils contain &
                   &holes (bad fringe)")
 103       format("# Zone ",a, ": ",a," donor stencils contain holes or &
                   &fringe")

           ! Terminate.

           call terminate("checkOverset", &
                          "Bad holes or overset data present")

         endif testRootProc

         ! The other processors wait to get killed.

         call mpi_barrier(SUmb_comm_world, ierr)

       endif testBadPresent

       ! Deallocate the memory of bad and badGlobal.

       deallocate(bad, badGlobal, stat=ierr)
       if(ierr /= 0)                   &
         call terminate("checkOverset", &
                        "Deallocation failure for bad,badGlobal")

       end subroutine checkOverset
