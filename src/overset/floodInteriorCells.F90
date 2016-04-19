subroutine floodInteriorCells(level, sps)
  use communication
  use blockPointers
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: level, sps

  ! Working
  integer(kind=intType) :: nn, i, j, k, nSeed, iSeed, ierr
  integer(kind=intType), dimension(:, :), allocatable :: stack, floodSeeds
  integer(kind=intType) :: nChanged, nChangedLocal, stackPointer, loopIter
  logical :: tmpSave, isCompute, isWallDonor, isFloodSeed, isHole
  integer(kind=intType), dimension(:, :, :), pointer :: changed

  ! Allocate pointer space for the integer flag communication
  do nn=1, nDom
     call setPointers(nn, level, sps)
     ! Note that it has to start at 1 since this is normally a pointer
     ! so the exchange routine expects the ordering to start at 1.
     allocate(flowDoms(nn, level, sps)%intCommVars(1)%var(1:ib+1, 1:jb+1, 1:kb+1))
     flowDoms(nn, level, sps)%intCommVars(1)%var = 0
  end do

  ! Keep track of the total number of loops
  loopIter = 1

  parallelSyncLoop: do 

     ! Keep track of the total number of fringes we've modified
     nChangedLocal = 0

     do nn=1,nDom
        call setPointers(nn, level, sps)
        changed => flowDoms(nn, level, sps)%intCommVars(1)%var

        ! Allocate space for our queue (stack). It needs to be 6*nx*ny*nz + 1:
        ! 6 for each of the 6 coordinate directions plus our extra
        ! seed. It should never come close to this unless the entire
        ! block will be blanked.

        allocate(stack(3, nx*ny*nz*6 + 1))

        ! Also allocate space for our flood seeds. Make it big enough to
        ! include the first level halos. 
        allocate(floodSeeds(3, 6*ie*je*je))

        ! These are the seeds we have directly. We will only use these on the first iteration:
        nSeed = 0

        if (loopIter == 1) then 

           ! Make the -3 and -2 cells, those inside the body,
           ! "compute" cells. This allows the flooding algorithm to
           ! work the same on subsequent outer iterations as the first
           ! iteration. 
           do k=2, kl
              do j=2, jl
                 do i=2, il
                    if (iblank(i, j, k) == -3 .or. iblank(i, j, k) == -2) then 
                       call setIsCompute(fringes(i, j, k)%status, .True.)
                       fringes(i, j, k)%donorProc = -1
                    end if
                 end do
              end do
           end do

           do k=2, kl
              do j=2, jl
                 do i=2, il
                       if (isWallDonor(fringes(i, j, k)%status) .and. &
                            isCompute(fringes(i, j, k)%status)) then 
                       call addSeed(i,j ,k)
                    end if
                 end do
              end do
           end do
        else
           ! On the second and subsequent passes, check each 1st
           ! non-corner halos in the 6 faces to see if we received
           ! "changed" info from neighbour proc. This will allow us to
           ! continue the flooding on this processor/block. Note that
           ! even in a single processor case, the halo exchange in
           ! necessary to communicate between two local blocks

           ! iMin/iMax
           do k=2, kl
              do j=2, jl
                 if (changed(1, j, k) == 1) then
                    call addSeed(2, j, k)
                 end if
                 if (changed(ie, j, k) == 1) then
                    call addSeed(il, j, k)
                 end if
              end do
           end do

           ! jMin/jMax
           do k=2, kl
              do i=2, il
                 if (changed(i, 1, k) == 1) then 
                    call addSeed(i, 2, k)
                 end if
                 if (changed(i, je, k) == 1) then 
                    call addSeed(i, jl, k)
                 end if
              end do
           end do

           ! kMin:
           do j=2, jl
              do i=2, il
                 if (changed(i, j, 1) == 1) then
                    call addSeed(i, j, 2)
                 end if
                 if (changed(i, j, ke) == 1) then
                    call addSeed(i, j, kl)
                 end if
              end do
           end do
        end if

        ! Loop over our seeds we currently have
        do iSeed = 1, nSeed

           ! Put the particular seed in the first slot of the stack
           stack(:, 1) = floodSeeds(:, iSeed)

           ! Reset the stack pointer length back to just the one seed
           ! we have
           stackPointer = 1

           ! flag the seed points --- only on first pass
           if (loopIter == 1) then 
              i = stack(1, stackPointer)
              j = stack(2, stackPointer)
              k = stack(3, stackPointer)
              call setIsFloodSeed(fringes(i, j, k)%status, .True. )
           end if

           ! Start the flooding (stacked based, not recursive)
           do while (stackPointer > 0 )

              ! 'Pop' the current point off the stack
              i = stack(1, stackPointer)
              j = stack(2, stackPointer)
              k = stack(3, stackPointer)
              stackPointer = stackPointer - 1

              if (isCompute(fringes(i, j, k)%status) .and. fringes(i, j, k)%donorProc == -1) then 
                 ! Flag the cell (using changed) as being changed
                 changed(i, j, k) = 1

                 ! Keep track of the total number we've changed.
                 nChangedLocal = nChangedLocal + 1

                 ! Pure compute cell, convert to hole
                 tmpSave = isFloodSeed(fringes(i, j, k)%status)
                 call emptyFringe(fringes(i, j, k))
                 call setIsHole(fringes(i, j, k)%status, .True.)
                 call setIsFlooded(fringes(i, j, k)%status, .True.)
                 call setIsCompute(fringes(i, j, k)%status, .False.)
                 call setIsFloodSeed(fringes(i, j, k)%status, tmpSave)
                 
                 ! Now add the six nearest neighbours to the stack
                 ! provided they are in the owned cell range:

                 if (i-1 >= 2) then 
                    stackPointer = stackPointer + 1
                    stack(:, stackPointer) = (/i-1, j  , k  /)
                 end if

                 if (i+1 <= il) then 
                    stackPointer = stackPointer + 1
                    stack(:, stackPointer) = (/i+1, j  , k  /)
                 end if

                 if (j-1 >= 2) then 
                    stackPointer = stackPointer + 1
                    stack(:, stackPointer) = (/i  , j-1, k  /)
                 end if

                 if (j+1 <= jl) then 
                    stackPointer = stackPointer + 1
                    stack(:, stackPointer) = (/i  , j+1, k  /)
                 end if

                 if (k-1 >= 2) then 
                    stackPointer = stackPointer + 1
                    stack(:, stackPointer) = (/i  , j  , k-1/)
                 end if

                 if (k+1 <= kl) then 
                    stackPointer = stackPointer + 1
                    stack(:, stackPointer) = (/i  , j , k+1 /)
                 end if
              end if

           end do
        end do

        deallocate(stack, floodSeeds)
     end do

    ! Exchange "changed"
     call wHalo1to1IntGeneric(1, level, sps, commPatternCell_1st, internalCell_1st)

     ! Determine if cells got changd. If so do another loop.
     call mpi_allreduce(nChangedLocal, nChanged, 1, sumb_integer, MPI_SUM, &
          sumb_comm_world, ierr)
     call ECHK(ierr, __FILE__, __LINE__)

     if (myid == 0) then 
        print *, 'Flood Iteration:', loopIter, 'Blanked ', nChanged, 'Interior Cells.'
     end if

     if (nChanged == 0) then 
        exit parallelSyncLoop
     end if

     loopIter = loopIter + 1

  end do parallelSyncLoop

  ! deallocate the temporary int space
  do nn=1, nDom
     call setPointers(nn, level, sps)
     deallocate(flowDoms(nn, level, sps)%intCommVars(1)%var)
  end do

contains
  ! Simple routine to make code easier to read above
  subroutine addSeed(i, j, k)
    use precision
    implicit none
    integer(kind=intType), intent(in) :: i, j, k
    nSeed = nSeed + 1
    floodSeeds(:, nSeed) = (/i, j, k/)
  end subroutine addSeed
end subroutine floodInteriorCells
