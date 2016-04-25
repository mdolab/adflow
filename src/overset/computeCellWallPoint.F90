subroutine computeCellWallPoint(level, sps)

  !  This routine is vastly more complex that it really should
  !  be. Essentially what we want to do is to determine the "wall
  !  point" for every cell. Essentially for every wall surface, we wan
  !  to record the coordintes of the wall surface cell center along
  !  cell centers eminating from the surface. The reason why this gets
  !  complex is that block can get cut in the off-wall direction which
  !  breaks the propagation. If this propatation isn't continued, the
  !  overset hole cut will be dependent on the block distribution and
  !  hense the numbe rof processors. 

  use blockPointers
  use BCTypes
  use communication
  implicit none 

  ! Input Params
  integer(kind=intType), intent(in) :: level, sps

  ! Working paramters
  integer(kind=intType) :: i, j, k, nn, mm, ierr, nFlaggedLocal
  integer(kind=intType) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
  integer(kind=intType) :: loopIter, nSeedsSet, nChanged, curCell
  real(kind=realType) :: curSeed(3)
  logical :: isWallType

  ! Since we need to have the 'x' and 'nearWall' arrays allocated at
  ! once, these need to go into block. Allocate the two x arrays the
  ! near wall arrays
 
  nSeedsSet = 0
  do nn=1, nDom
     call setPointers(nn, level, sps)

     if (.not. associated(flowDoms(nn, level, sps)%xSeed)) then 
        allocate(flowDoms(nn, level, sps)%XSeed(0:ib, 0:jb, 0:kb, 3))
        allocate(flowDoms(nn, level, sps)%wallInd(0:ib, 0:jb, 0:kb))
        ! Manaully set the pointer for xSeed
        xSeed => flowDoms(nn, level, sps)%xSeed
        wallInd => flowDoms(nn, level, sps)%wallInd
     end if

     ! Initialize to large to indicate that nothing has been changed. 
     xSeed = large
     wallInd = -1

     ! Flag all nodes that are within nearWallDist as being nearWall
     do mm=1, nBocos
        if (isWallType(BCType(mm))) then 
    
           ! Set the bounds for the triple loop. We loop over all
           ! cells that are 
           select case (BCFaceID(mm))
           case (iMin, iMax)

              do k=BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                 do j=BCData(mm)%inBeg+1,BCData(mm)%inEnd
                    if (BCFaceID(mm) == iMin) then 
                       curSeed = fourth*(x(1, j-1, k-1, :) + x(1, j, k-1, :) + &
                            x(1, j, k, :) + x(1, j-1, k, :))
                       curCell = globalCell(2, j, k)
                    else
                       curSeed = fourth*(x(il, j-1, k-1, :) + x(il, j, k-1, :) + &
                            x(il, j, k, :) + x(il, j-1, k, :))
                       curCell = globalCell(il, j, k)
                    end if
                    
                    do i=2, il
                       call addSeed()
                    end do
                 end do
              end do

           case (jMin, jMax)

              do k=BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                 do i=BCData(mm)%inBeg+1,BCData(mm)%inEnd
                    if (BCFaceID(mm) == jMin) then 
                       curSeed = fourth*(x(i-1, 1, k-1, :) + x(i, 1, k-1, :) + &
                            x(i, 1, k, :) + x(i-1, 1, k, :))
                       curCell = globalCell(i, 2, k)
                    else
                       curSeed = fourth*(x(i-1, jl, k-1, :) + x(i, jl, k-1, :) + &
                            x(i, jl, k, :) + x(i-1, jl, k, :))
                       curCell = globalCell(i, jl, k)
                    end if
                    do j=2, jl
                       call addSeed()
                    end do
                 end do
              end do
           case (kMin, kMax)

              do j=BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                 do i=BCData(mm)%inBeg+1,BCData(mm)%inEnd
                    if (BCFaceID(mm) == kMin) then 
                       curSeed = fourth*(x(i-1, j-1, 1, :) + x(i, j-1, 1, :) + &
                            x(i, j, 1, :) + x(i-1, j, 1, :))
                       curCell = globalCell(i, j, 2)
                    else
                       curSeed = fourth*(x(i-1, j-1, kl, :) + x(i, j-1, kl, :) + &
                            x(i, j, kl, :) + x(i-1, j, kl, :))
                       curCell = globalCell(i, j, kl)
                    end if
                    do k=2, kl
                       call addSeed()
                    end do
                 end do
              end do
           end select
        end if
     end do
  end do

  ! Exchange the xSeeds 
  do nn=1, nDom
     flowDoms(nn, level, sps)%realCommVars(1)%var => &
          flowDoms(nn, level, sps)%xSeed(:, :, :, 1)
     flowDoms(nn, level, sps)%realCommVars(2)%var => & 
          flowDoms(nn, level, sps)%xSeed(:, :, :, 2) 
     flowDoms(nn, level, sps)%realCommVars(3)%var => &    
          flowDoms(nn, level, sps)%xSeed(:, :, :, 3) 
     flowDoms(nn, level, sps)%intCommVars(1)%var => & 
          flowDoms(nn, level, sps)%wallInd(:, :, :)
  end do
  
  ! Run the generic halo exchange. Only need to use the first level
  ! halos. 
  call wHalo1to1RealGeneric(3, level, sps, commPatternCell_1st, internalCell_1st)
  call wHalo1to1IntGeneric(1, level, sps, commPatternCell_1st, internalCell_1st)

   ! Iterative loop
  loopIter = 1
  parallelSyncLoop: do 
     nSeedsSet = 0
     do nn=1, nDom
        call setPointers(nn, level, sps)

        ! Loop over all halos. We won't be selective here, since the
        ! seeds cound show up anywhere! :-)

        ! Generic loop over all boundaries
        do mm=1, 6
    
           select case (BCFaceID(mm))
           case (iMin, iMax)

              do k=2, kl
                 do j=2, jl
                    if (BCFaceID(mm) == iMin) then 
                       ! Extract the seed from the halo
                       curSeed = xSeed(1, j, k, :)
                       curCell = wallInd(1, j, k)
                    else
                       curSeed = xSeed(ie, j, k, :)
                       curCell = wallInd(ie, j, k)
                    end if
                    if (curSeed(1) /= large) then 
                       do i=2, il
                          call addSeed()
                       end do
                    end if
                 end do
              end do

           case (jMin, jMax)

              do k=2, kl
                 do i=2, il
                    if (BCFaceID(mm) == jMin) then 
                       curSeed = xSeed(i, 1, k, :)
                       curCell = wallInd(i, 1, k)
                    else
                       curSeed = xSeed(i, je, k, :)
                       curCell = wallInd(i, je, k)
                    end if
                    if (curSeed(1) /= large) then 
                       do j=2, jl
                          call addSeed()
                       end do
                    end if
                 end do
              end do
           case (kMin, kMax)

              do j=2, jl
                 do i=2, il
                    if (BCFaceID(mm) == kMin) then 
                       curSeed = xSeed(i, j, 1, :)
                       curCell = wallInd(i, j, 1)
                    else
                       curSeed = xSeed(i, j, ke, :)
                       curCell = wallInd(i, j, ke)
                    end if
                    if (curSeed(1) /= large) then 
                       do k=2, kl
                          call addSeed()
                       end do
                    end if
                 end do
              end do
           end select
        end do
     end do

     ! Exchange the xSeeds to make sure everything is up to date.
     do nn=1, nDom
        flowDoms(nn, level, sps)%realCommVars(1)%var => &
             flowDoms(nn, level, sps)%xSeed(:, :, :, 1)
        flowDoms(nn, level, sps)%realCommVars(2)%var => & 
             flowDoms(nn, level, sps)%xSeed(:, :, :, 2) 
        flowDoms(nn, level, sps)%realCommVars(3)%var => &    
             flowDoms(nn, level, sps)%xSeed(:, :, :, 3) 
        flowDoms(nn, level, sps)%intCommVars(1)%var => & 
             flowDoms(nn, level, sps)%wallInd(:, :, :)
     end do
     
     ! Run the generic halo exchange. Only need to use the first level
     ! halos. 
     call wHalo1to1RealGeneric(3, level, sps, commPatternCell_1st, internalCell_1st)
     call wHalo1to1IntGeneric(1, level, sps, commPatternCell_1st, internalCell_1st)

     ! Determine if any cells made it to the other side of a face. If so, we have to keep going:
     call mpi_allreduce(nSeedsSet, nChanged, 1, sumb_integer, MPI_SUM, &
          sumb_comm_world, ierr)
     call ECHK(ierr, __FILE__, __LINE__)

     if (nChanged == 0) then 
        exit parallelSyncLoop
     end if

     if (myid == 0) then 
        print *, 'Flag Near Wall Iteration:', loopIter, 'nAtBoundary', nChanged
     end if

     loopIter = loopIter + 1
 
  end do parallelSyncLoop

contains 

  subroutine addSeed

    implicit none
    real(kind=realType) :: curDist, newDist, xp(3)

    if (xSeed(i, j, k, 1) == large) then 
       ! Not added yet. Just set. 
       xSeed(i, j, k, :) = curSeed
       wallInd(i, j, k) = curCell
       nSeedsSet = nSeedsSet + 1
    else
       ! We have multiple faces and the seed
       ! already exists. Take the one that is
       ! closest. Reconstruct the cell center point
       ! first. 
       xp = eighth*(&
            x(i-1, j-1, k-1, :) + &
            x(i  , j-1, k-1, :) + &
            x(i-1, j  , k-1, :) + &
            x(i  , j  , k-1, :) + &
            x(i-1, j-1, k  , :) + &
            x(i  , j-1, k  , :) + &
            x(i-1, j  , k  , :) + &
            x(i  , j  , k  , :))
       curDist = norm2(xp - xSeed(i, j, k, :))
       newDist = norm2(xp - curSeed)
       
       if (newDist < curDist) then 
          xSeed(i, j, k, :) = curSeed
          wallInd(i, j, k) = curCell
          nSeedsSet = nSeedsSet + 1
       end if
    end if
    
  end subroutine addSeed
      
end subroutine computeCellWallPoint
