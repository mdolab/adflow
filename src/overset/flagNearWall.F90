subroutine flagNearWallCells(level, sps)

  ! This routine is vastly more complex that it really should
  !  be. Essentially what we want to do is flag cells that are
  !  approximately within a distance "nearWallDistance" of a
  !  wall. Nominally, this just means looking at planes parallel to a
  !  wall boundary and determining the distance to the point on the
  !  surface. Very quick (no searching). However, if a block is split
  !  parallel to a wall BC which is not desirable, but will happen
  !  occasionally for block partitioning purposes, we can have "near
  !  wall" cells that are on a split block and have no idea they are
  !  near a wall. So what we have to do is similar to the flooding
  !  algorithm: Starting at the wall, you proceed outwards. If you get
  !  to the other side of the block and it is still near a wall, you
  !  communicate and then continue on the other side. For the fringe
  !  search, (which is built from the dual mesh), we want to know of
  !  the *cells* of the dual mesh are near a wall. That means we
  !  actually want check the distnaces of the dual of the
  !  dual. Basically that means we're back to the primal. 

  use constants
  use overset
  use inputOverset
  use blockPointers
  use adtAPI
  use BCTypes
  use cgnsGrid
  use communication
  implicit none 

  ! Input Params
  integer(kind=intType), intent(in) :: level, sps

  ! Working paramters
  integer(kind=intType) :: i, j, k, nn, mm, ierr, nFlaggedLocal
  integer(kind=intType) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
  integer(kind=intType) :: loopIter, nAtBoundaryLocal, nAtBoundary
  logical :: isWallType, tostop
  type(Xplane), dimension(:), allocatable :: planes
  real(Kind=realType), dimension(:, :, :, :), pointer :: xx, xSeed
  integer(kind=intType), dimension(:, :, :), pointer :: nearWall
  real(kind=realType) :: dist
  ! Since we need to have the 'x' and 'nearWall' arrays allocated at
  ! once, these need to go into block. Allocate the two x arrays the
  ! near wall arrays
  nAtBoundaryLocal = 0
  nFlaggedLocal = 0
  do nn=1, nDom
     call setPointers(nn, level, sps)
     allocate(flowDoms(nn, level, sps)%XSeed(0:ie, 0:je, 0:ke, 3), &
          flowDoms(nn, level, sps)%nearWall(1:il, 1:jl, 1:kl))

     ! Manaully set the three pointers
     xSeed => flowDoms(nn, level, sps)%xSeed
     nearWall => flowDoms(nn, level, sps)%nearWall

     ! Use large to denote no seed presnet
     xSeed = large

     ! Flag all nodes that are within nearWallDist as being nearWall
     do mm=1, nBocos
        if (isWallType(BCType(mm))) then 
    
           call setBoundaryPointers(mm, BCFaceID(mm), .False.)
           ! Loop over the generalized plane
           do j=jStart, jEnd
              do i=iStart, iEnd

                 ! Loop over the 'k' ie offwall direction
                 do k=kStart, kEnd

                    dist = norm2(planes(k)%x(i, j, :) - planes(1)%x(i, j, :))
                    if (dist < nearWallDist) then 
                       planes(k)%nearWall(i, j) = 1
                       nFlaggedLocal = nFlaggedLocal + 1

                       ! If we made it all the way to the other side
                       ! of the block, set the wall point into the
                       ! xSeed array. Note that we need the *SECOND
                       ! LAST NODE*. This is the value that is
                       ! actually transfered as it is the halo node
                       ! for the other block. 

                       if (k == kEnd-1) then
                          planes(k)%xSeed(i, j, :) = planes(1)%x(i, j, :)
                          nAtBoundaryLocal = nAtBoundaryLocal + 1
                       end if
                    end if
                 end do
              end do
           end do
           deallocate(planes)
        end if
     end do ! BocoLoop
  end do

  ! Determine if any cells made it to the other side of a face. If so, we have to keep going:
  call mpi_allreduce(nAtBoundaryLocal, nAtBoundary, 1, sumb_integer, MPI_SUM, &
       sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)
 
  ! Iterative loop
  loopIter = 1
  parallelSyncLoop: do while (nAtBoundary > 0)
     if (myid == 0) then 
        print *, 'Flag Near Wall Iteration:', loopIter, 'nAtBoundary', nAtBoundary
     end if

     ! Reset the counter
     nAtBoundaryLocal = 0

     ! Exchange the xSeeds after the initial flooding.
     call exchangeXSeeds(level, sps, commPatternNode_1st, internalNode_1st)

     do nn=1, nDom
        call setPointers(nn, level, sps)

        ! Manaully set the additional pointers
        xSeed => flowDoms(nn, level, sps)%xSeed
        nearWall => flowDoms(nn, level, sps)%nearWall

        ! Loop over all halos. We won't be selective here, since the
        ! seeds cound show up anywhere! :-)

        ! Generic loop over all boundaries
        do mm=1, 6
           call setBoundaryPointers(mm, mm, .True.)

           ! Loop over the size of the generalized plane
           do j=jStart, jEnd
              do i=iStart, iEnd

                 ! Determine if we need to do the generalized 'k'
                 ! direction at all. We only need to do if a valid
                 ! seed has shown up in the xSeed:
                 if (planes(0)%xSeed(i, j, 1) < large) then 
        
                    ! Loop over the 'k' ie offwall direction
                    do k=kStart, kEnd

                       dist = norm2(planes(k)%x(i, j, :) - planes(0)%xSeed(i, j, :))
                       if (dist < nearWallDist .and. planes(k)%nearWall(i, j) == 0) then 
                          planes(k)%nearWall(i, j) = 1
                                               
                          ! If we made it all the way to the other side
                          ! of the block, copy over the seed for the next exchange. 
                          if (k == kEnd-1) then 
                             planes(k)%xSeed(i, j, :) = planes(0)%xSeed(i, j, :)
                             nAtBoundaryLocal = nAtBoundaryLocal + 1
                          end if
                       end if
                    end do
                 end if
              end do
           end do
           deallocate(planes)
        end do 
     end do

     ! Determine if any cells made it to the other side of a face. If so, we have to keep going:
     call mpi_allreduce(nAtBoundaryLocal, nAtBoundary, 1, sumb_integer, MPI_SUM, &
          sumb_comm_world, ierr)
     call ECHK(ierr, __FILE__, __LINE__)

     loopIter = loopIter + 1
 
  end do parallelSyncLoop

  ! Deallocate X and XSeed since they are no longer needed. We
  ! have to hold onto nearWall a little longer since we need to use it
  ! in initializeOBlock. It will deallocate it. 
 do nn=1, nDom
    deallocate(flowDoms(nn, level, sps)%XSeed)
  end do

contains 
  subroutine setBoundaryPointers(mm, faceID, fullFaces)
    implicit none

    integer(kind=intType), intent(in) :: mm, faceID
    logical, intent(in) :: fullFaces

    if (.not. fullFaces) then 
       iStart=BCData(mm)%inBeg; iEnd=BCData(mm)%inEnd
       jStart=BCData(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
    else
       iStart=1; jStart=1
       select case (mm)
       case (iMin, iMax)
          iEnd=jl; jEnd=kl
       case (jMin, jMax) 
          iEnd=il; jEnd=kl
       case(kMin, kMax)
          iEnd=il; jEnd=jl
       end select
    end if

    select case (faceID)
    case (iMin)
       kStart=1; kend=il
       allocate(planes(0:il))
       planes(0)%xSeed => xSeed(0, 1:jl, 1:kl, :)
       do i=1, il
          planes(i)%x => x(i, 1:jl, 1:kl, :)
          planes(i)%xSeed => xSeed(i, 1:jl, 1:kl, :)
          planes(i)%nearWall => nearWall(i, :, :)
       end do
    case (iMax)
       kStart=1; kend=il
       allocate(planes(0:il))
       planes(0)%xSeed => xSeed(ie, 1:jl, 1:kl, :)
       do i=1, il
          planes(i)%x => x( il-i+1, 1:jl, 1:kl, :)
          planes(i)%xSeed => xSeed(il-i+1, 1:jl, 1:kl, :)
          planes(i)%nearWall => nearWall(il-i+1, :, :)
       end do
    case (jMin)
       kStart=1; kend=jl
       allocate(planes(0:jl))
       planes(0)%xSeed => xSeed(1:il, 0, 1:kl, :)
       do j=1, jl
          planes(j)%x => x(1:il, j, 1:kl, :)
          planes(j)%xSeed => xSeed(1:il, j, 1:kl, :)
          planes(j)%nearWall => nearWall(:, j, :)
       end do
    case (jMax)
       kStart=1; kend=jl
       allocate(planes(0:jl))
       planes(0)%xSeed => xSeed(1:il, je, 1:kl, :)
       do j=1, jl
          planes(j)%x => x(1:il, jl-j+1, 1:kl, :)
          planes(j)%xSeed => xSeed(1:il, jl-j+1, 1:kl, :)
          planes(j)%nearWall => nearWall(:, jl-j+1, :)
       end do
    case (kMin)
       kStart=1; kend=kl
       allocate(planes(0:kl))
       planes(0)%xSeed => xSeed(1:il, 1:jl, 0, :)
       do k=1, kl
          planes(k)%x => x(1:il, 1:jl, k, :)
          planes(k)%xSeed => xSeed(1:il, 1:jl, k, :)
          planes(k)%nearWall => nearWall(:, :, k)
       end do
    case (kMax)
       kStart=1; kend=kl
       allocate(planes(0:kl))
       planes(0)%xSeed => xSeed(1:il, 1:jl, ke, :)
       do k=1, kl
          planes(k)%x => x(1:il, 1:jl, kl-k+1, :)
          planes(k)%xSeed => xSeed(1:il, 1:jl, kl-k+1, :)
          planes(k)%nearWall => nearWall(:, :, kl-k+1)
       end do
    end select
  end subroutine setBoundaryPointers
end subroutine flagNearWallCells
      
 
