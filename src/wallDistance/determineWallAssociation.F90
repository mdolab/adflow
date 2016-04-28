subroutine determineWallAssociation(level, sps)

  ! This routine will determine the closest surface point for every
  ! field cell. Special treatment is required for overlapping surfaces. 

  use adtAPI
  use blockPointers
  use wallDistanceData
  use BCTypes
  use communication
  use inputphysics
  use inputTimeSpectral
  use overset, only : oversetPresent, oversetWall, nClusters, clusters, cumDomProc
  use inputOverset
  use adjointVars
  implicit none

  ! Input Variables
  integer(kind=intType), intent(in) :: level, sps

  ! Local Variables
  integer(kind=intType) :: i, j, k, l, ii, jj, kk, nn, mm, iNode, iCell, c
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ni, nj, nUnique, cellID, cellID2
  integer(kind=intType) :: ierr, iDim

  ! Data for local surface
  integer(kind=intType) :: nNodes, nCells
  logical :: gridHasOverset

  ! Overset Walls for storing the surface ADT's
  type(oversetWall), dimension(:), allocatable, target :: walls
  type(oversetWall), target :: fullWall
  integer(kind=intType), dimension(:),  allocatable :: link, indicesToGet

  ! Data for the ADT
  integer(kind=intType) :: intInfo(3), intInfo2(3)
  real(kind=realType) :: coor(4), uvw(5), uvw2(5)
  real(kind=realType), dimension(3, 2) :: dummy
  real(kind=realType), parameter :: tol=1e-12
  integer(kind=intType), dimension(:), pointer :: frontLeaves, frontLeavesNew, BBint
  type(adtBBoxTargetType), dimension(:), pointer :: BB
  real(kind=realType), dimension(3) :: xp

  ! The first thing we do is gather all the surface nodes to
  ! each processor such that every processor can make it's own copy of
  ! the complete surface mesh to use to search. Note that this
  ! procedure *DOES NOT SCALE IN MEMORY*...ie eventually the surface
  ! mesh will become too large to store on a single processor,
  ! although this will probably not happen until the sizes get up in
  ! the hundreds of millions of cells. 

  allocate(walls(nClusters))
  call buildClusterWalls(level, sps, .False., walls)

  if (oversetPresent) then 
     ! Finally build up a "full wall" that is made up of all the cluster
     ! walls. 

     nNodes = 0
     nCells = 0
     do i=1, nClusters
        nNodes = nNodes+ walls(i)%nNodes
        nCells = nCells + walls(i)%nCells
     end do

     allocate(fullWall%x(3, nNodes))
     allocate(fullWall%conn(4, nCells))
     allocate(fullWall%ind(nNodes))

     nNodes = 0
     nCells = 0
     ii = 0
     do i=1, nClusters

        ! Add in the nodes/elements from this cluster

        do j=1, walls(i)%nNodes
           nNodes = nNodes + 1
           fullWall%x(:, nNodes) = walls(i)%x(:, j)
           fullWall%ind(nNodes) = walls(i)%ind(j)
        end do

        do j=1, walls(i)%nCells
           nCells = nCells + 1
           fullWall%conn(:, nCells) = walls(i)%conn(:, j) + ii
        end do

        ! Increment the node offset
        ii = ii + walls(i)%nNodes
     end do

     ! Finish the setup of the full wall.
     fullWall%nCells = nCells
     fullWall%nNodes = nNodes
     call buildSerialQuad(nCells, nNodes, fullWall%x, fullWall%conn, fullWall%ADT)
  end if

  ! Allocate the (pointer) memory that may be resized as necessary for
  ! the singlePoint search routine. 
  allocate(stack(100), BB(20), BBint(20), frontLeaves(25), frontLeavesNew(25))

  ! We need to store the 4 global node indices defining the quad that
  ! each point has the closest point wrt. We also ned to store the uv
  ! values. This allows us to recompute the exact surface point, after
  ! the rquired nodes are fetched from (a possibly) remote proc. 

  do nn=1,nDom
     call setPointers(nn, level, sps)

     ! Check if elemID and uv are allocated yet.
     if (.not. associated(flowDoms(nn,level,sps)%surfNodeIndices)) then
        allocate(flowDoms(nn,level,sps)%surfNodeIndices(4, 2:il, 2:jl, 2:kl))
        allocate(flowDoms(nn,level,sps)%uv(2, 2:il, 2:jl, 2:kl))
     end if

     ! Set the cluster for this block
     c = clusters(cumDomProc(myid) + nn)

     do k=2, kl
        do j=2, jl
           do i=2, il

              ! Compute the coordinates of the cell center 
              coor(1) = eighth*(x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
                   +         x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
                   +         x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                   +         x(i-1,j,  k,  1) + x(i,j,  k,  1))

              coor(2) = eighth*(x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
                   +         x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
                   +         x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                   +         x(i-1,j,  k,  2) + x(i,j,  k,  2))

              coor(3) = eighth*(x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
                   +         x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
                   +         x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                   +         x(i-1,j,  k,  3) + x(i,j,  k,  3))

              if (.not. oversetPresent) then 
                 ! No overset present. Simply search our own wall,
                 ! walls(c), (the only one we have) up to the wall
                 ! cutoff.
                 coor(4) = wallDistCutoff**2
                 intInfo(3) = 0 ! Must be initialized since the search
                 ! may not find closer point.
                 call minDistancetreeSearchSinglePoint(walls(c)%ADT, coor, intInfo, &
                      uvw, dummy, 0, BB, frontLeaves, frontLeavesNew)

                 cellID = intInfo(3)
                 if (cellID > 0) then 
                    do kk=1,4
                       flowDoms(nn, level, sps)%surfNodeIndices(kk, i, j, k) = &
                            walls(c)%ind(walls(c)%conn(kk, cellID))
                    end do
                    flowDoms(nn, level, sps)%uv(:, i, j, k) = uvw(1:2)
                 else
                    ! Just set dummy values. These will never be used. 
                    flowDoms(nn, level, sps)%surfNodeIndices(:, i, j, k) = 0
                    flowDoms(nn, level, sps)%uv(:, i, j, k) = 0
                 end if

                 ! We are done with this point. 
                 cycle
              end if

              ! This is now the overset (possibly) overlapping surface
              ! mesh case. It is somewhat more complex since we use
              ! the same searches to flag cells that are inside the
              ! body.

              coor(4) = wallDistCutoff**2
              intInfo(3) = 0
              call minDistancetreeSearchSinglePoint(fullWall%ADT, coor, &
                   intInfo, uvw, dummy, 0, BB, frontLeaves, frontLeavesNew)
              cellID = intInfo(3)

              if (cellID > 0) then
                 ! We found the cell:

                 ! If the cell is outside of near-wall distance or our
                 ! cluster doesn't have any owned cells. Just accept it. 
                 if (uvw(4) > nearWallDist**2 .or. walls(c)%nCells == 0) then 

                    do kk=1,4
                       flowDoms(nn, level, sps)%surfNodeIndices(kk, i, j, k) = &
                            fullWall%ind(fullWall%conn(kk, cellID))
                    end do
                    flowDoms(nn, level, sps)%uv(:, i, j, k) = uvw(1:2)

                 else

                    ! This point is *closer* than the nearWallDist AND
                    ! it has a wall. Search on our own wall.

                    coor(4) = large
                    call minDistancetreeSearchSinglePoint(walls(c)%ADT, coor, &
                         intInfo2, uvw2, dummy, 0, BB, frontLeaves, frontLeavesNew)
                    cellID2 = intInfo2(3)

                    if (uvw2(4) < nearWallDist**2) then 
                       ! Both are close to the wall. Accept the one
                       ! from our own wall unconditionally.
                       do kk=1,4
                          flowDoms(nn, level, sps)%surfNodeIndices(kk, i, j, k) = &
                               walls(c)%ind(walls(c)%conn(kk, cellID2))
                       end do
                       flowDoms(nn, level, sps)%uv(:, i, j, k) = uvw2(1:2)
                    else
                       ! The full wall distance is better. Take that. 

                       do kk=1,4
                          flowDoms(nn, level, sps)%surfNodeIndices(kk, i, j, k) = &
                               fullWall%ind(fullWall%conn(kk, cellID))
                       end do
                       flowDoms(nn, level, sps)%uv(:, i, j, k) = uvw(1:2)

                    end if
                 end if
              else

                 ! What happend here is a cell is outside the
                 ! wallDistCutoff. We don't care about wall distance
                 ! info here so just set dummy info.

                 flowDoms(nn, level, sps)%surfNodeIndices(:, i, j, k) = 0
                 flowDoms(nn, level, sps)%uv(:, i, j, k) = 0

              end if
           end do
        end do
     end do
  end do

  ! Now determine all the node indices this processor needs to get. 
  mm = 0
  allocate(indicesToGet(nCellsLocal(level)*4), link(nCellsLocal(level)*4))
  do nn=1, nDom
     call setPointers(nn, level, sps)
     do k=2, kl
        do j=2, jl
           do i=2, il
              do kk=1,4
                 mm = mm + 1
                 indicesToGet(mm) = flowDoms(nn, level, sps)%surfNodeIndices(kk, i, j, k)
              end do
           end do
        end do
     end do
  end do

  ! This unique-ifies the indices. 
  call unique(indicesToGet, 4*nCellsLocal(level), nUnique, link)

  ! we need to update the stored indices to use the ordering of the nodes we will receive. 
  mm = 0
  do nn=1, nDom
     call setPointers(nn, level, sps)
     do k=2, kl
        do j=2, jl
           do i=2, il
              do kk=1,4
                 mm = mm + 1
                 flowDoms(nn, level, sps)%surfNodeIndices(kk, i, j, k) = link(mm)
              end do
           end do
        end do
     end do
  end do
  deallocate(link)

  ! Now create the index set for the nodes we need to get. We have to
  ! expand "indices to get" to include the DOF. Use link for this
  ! temporary array operation.

  allocate(link(nUnique*3))
  do i=1, nUnique
     link((i-1)*3+1) = indicesToGet(i)*3
     link((i-1)*3+2) = indicesToGet(i)*3+1
     link((i-1)*3+3) = indicesToGet(i)*3+2
  end do

  call ISCreateGeneral(sumb_comm_world, nUnique*3, link, PETSC_COPY_VALUES, IS1, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  deallocate(link)

  ! Create the volume vector the nodes will be scatter from. Note that
  ! this vector contains all the spectal instances. It is therefore
  ! only allocated on the first call with sps=1
  if (sps == 1) then 
     call VecCreateMPI(SUMB_COMM_WORLD, 3*nNodesLocal(level)*nTimeIntervalsSpectral, &
          PETSC_DETERMINE, xVolumeVec(level), ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end if

  ! This is the vector we will scatter the nodes into. 
  call VecCreateMPI(SUMB_COMM_WORLD, 3*nUnique, PETSC_DETERMINE, &
       xSurfVec(level, sps), ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecGetOwnershipRange(xSurfVec(level, sps), i, j, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call ISCreateStride(SUMB_COMM_WORLD, j-i, i, 1, IS2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Create the actual final scatter context.
  call  VecScatterCreate(xVolumeVec(level), IS1, xSurfVec(level, sps), IS2, &
       wallScatter(level, sps), ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call ISDestroy(IS1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call ISDestroy(IS2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Deallocate all the remaining temporary data
  deallocate(stack, BB, frontLeaves, frontLeavesNew, BBint)

  do i=1, nClusters
     deallocate(walls(i)%x, walls(i)%conn, walls(i)%ind)
     call destroySerialQuad(walls(i)%ADT)
  end do
  deallocate(walls)

  if (oversetPresent) then 
     deallocate(fullWall%x, fullWall%conn, fullWall%ind)
     call destroySerialQuad(fullWall%ADT)
  end if

end subroutine determineWallAssociation
