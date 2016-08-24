subroutine computeHolesInsideBody(level, sps)

  ! This routine will flag the iBlank values in oBlocks with 0 if the
  ! the cell center falls inside the body. This is a (semi) parallel
  ! implementation: The global surface mesh is communicated to all
  ! processors who then search accordingly. This scalable in terms of
  ! computation but not strictly memory. 

  use adtAPI
  use blockPointers
  use wallDistanceData
  use communication
  use inputphysics
  use inputTimeSpectral
  use overset
  use inputOverset
  use adjointVars, only : totalVolumeNodes => nNodesLocal, totalVolumeCells => nCellsLocal
  implicit none

  ! Input Variables
  integer(kind=intType), intent(in) :: level, sps

  ! Local Variables
  integer(kind=intType) :: i, j, k, l, ii, jj, kk, nn, mm, iNode, iCell, c
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ni, nj, nUnique, cellID, cellID2
  integer(kind=intType) :: ierr, iDim

  ! Data for local surface
  integer(kind=intType) :: nNodes, nCells
  integer(kind=intType) :: nNodesLocal, nCellsLocal
  integer(kind=intType), dimension(:, :), allocatable :: connLocal
  integer(kind=intType), dimension(:), allocatable :: clusterNodeLocal
  integer(kind=intType), dimension(:), allocatable :: clusterCellLocal
  real(kind=realType), dimension(:, :), allocatable :: nodesLocal
  real(kind=realType), dimension(:,:,:), pointer :: xx
  integer(kind=intType), dimension(:,:), pointer :: ind
  real(kind=realType) :: myDist, otherDist, timeA
  logical :: regularOrdering, nearAWall

  ! Overset Walls for storing the surface ADT's
  type(oversetWall), dimension(:), allocatable, target :: walls
  type(oversetWall), target :: fullWall

  ! Data for global surface
  integer(kind=intTYpe) :: nNodesGlobal, nCellsGlobal
  integer(kind=intType), dimension(:, :), allocatable, target :: connGlobal
  real(kind=realType), dimension(:, :), allocatable, target :: nodesGlobal
  integer(kind=intType), dimension(:), allocatable, target :: nodeIndicesGlobal

  integer(kind=intType), dimension(:), allocatable :: nodesPerCluster, cellsPerCluster, cnc, ccc
  integer(kind=intType), dimension(:), allocatable :: clusterNodeGlobal
  integer(kind=intType), dimension(:), allocatable :: clusterCellGlobal
  integer(kind=intType), dimension(:), allocatable :: localNodeNums
  integer(kind=intType), dimension(:), allocatable :: nodeIndicesLocal

  integer(kind=intType), dimension(:),    allocatable :: nCellProc, cumCellProc
  integer(kind=intType), dimension(:),    allocatable :: nNodeProc, cumNodeProc
  real(kind=realType),   dimension(:, :), allocatable :: uniqueNodes
  integer(kind=intType), dimension(:),    allocatable :: link, indicesToGet

  ! Pointers for easier readibility
  integer(kind=intType), dimension(:, :), pointer :: conn
  real(kind=realType), dimension(:, :), pointer :: nodes, norm
  integer(kind=intType), dimension(:), pointer :: tmpInd

  ! Data for the ADT
  integer(kind=intType) :: intInfo(3), intInfo2(3)
  real(kind=realType) :: coor(4), uvw(5), uvw2(5)
  real(kind=realType), dimension(3, 2) :: dummy
  real(kind=realType), parameter :: tol=1e-12
  integer(kind=intType), dimension(:), pointer :: frontLeaves, frontLeavesNew, BBint
  type(adtBBoxTargetType), dimension(:), pointer :: BB
  
  ! Misc
  real(kind=realType) :: dp, shp(4)
  real(kind=realType), dimension(3) ::xp, normal, v1

  ! The first thing we do is gather all the surface nodes to
  ! each processor such that every processor can make it's own copy of
  ! the complete surface mesh to use to search. Note that this
  ! procedure *DOES NOT SCALE IN MEMORY*...ie eventually the surface
  ! mesh will become too large to store on a single processor,
  ! although this will probably not happen until the sizes get up in
  ! the hundreds of millions of cells. 
  timea = mpi_wtime()

  nNodesLocal = 0
  nCellsLocal = 0

  do nn=1,nDom
     call setPointers(nn, level, sps)
     do mm=1, nBocos
        if(BCType(mm) == NSWallAdiabatic .or. &
           BCType(mm) == NSWallIsothermal .or. &
           BCType(mm) == EulerWall) then
           iBeg = bcData(mm)%inBeg
           iEnd = bcData(mm)%inEnd
           jBeg = bcData(mm)%jnBeg
           jEnd = bcData(mm)%jnEnd

           nNodesLocal = nNodesLocal + &
                (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
           nCellsLocal = nCellsLocal + & 
                (iEnd - iBeg)*(jEnd - jBeg)
        end if
     end do
  end do

  ! Now communicate these sizes with everyone
  allocate(nCellProc(nProc), cumCellProc(0:nProc), &
           nNodeProc(nProc), cumNodeProc(0:nProc))

  call mpi_allgather(nCellsLocal, 1, sumb_integer, nCellProc, 1, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call mpi_allgather(nNodesLocal, 1, sumb_integer, nNodeProc, 1, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now make cumulative versions of these
  cumCellProc(0) = 0
  cumNodeProc(0) = 0
  do i=1,nProc
     cumCellProc(i) = cumCellProc(i-1) + nCellProc(i)
     cumNodeProc(i) = cumNodeProc(i-1) + nNodeProc(i)
  end do

  ! And save the total number of nodes and cells for reference
  nCellsGlobal = cumCellProc(nProc)
  nNodesGlobal = cumNodeProc(nProc)

  ! Allocate the space for the local nodes and element connectivity
  allocate(nodesLocal(3, nNodesLocal), connLocal(4, nCellsLocal), &
       clusterCellLocal(nCellsLocal), clusterNodeLocal(NNodesLocal), &
       nodeIndicesLocal(nNodesLocal))

  iCell = 0
  iNode = 0
  ! Second loop over the local walls
  do nn=1,nDom
     call setPointers(nn, level, sps)
     c = clusters(cumDomProc(myid) + nn)

     do mm=1,nBocos
        if(  BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal .or. &
             BCType(mm) == EulerWall) then

           select case (BCFaceID(mm))
           case (iMin)
              xx   => x(1,:,:,:)
              ind  => globalNode(1, :, :)
                
           case (iMax)
              xx   => x(il,:,:,:)
              ind  => globalNode(il, :, :)
                 
           case (jMin)
              xx   => x(:,1,:,:)
              ind  => globalNode(:, 1, :)

           case (jMax)
              xx   => x(:,jl,:,:)
              ind  => globalNode(:, jl, :)
              
           case (kMin)
              xx   => x(:,:,1,:)
              ind  => globalNode(:, :, 1)

           case (kMax)
              xx   => x(:,:,kl,:)
              ind  => globalNode(:, :, kl)

           end select

           ! We want to ensure that all the normals of the faces are
           ! consistent. To ensure this, we enforce that all normals
           ! are "into" the domain. Therefore we must treat difference
           ! faces of a block differently. For example for an iLow
           ! face, when looping over j-k in the regular way, results
           ! in in a domain inward pointing normal for iLow but
           ! outward pointing normal for iHigh. The same is true for
           ! kMin and kMax. However, it is reverse for the J-faces:
           ! This is becuase the way the pointers are extracted i then
           ! k is the reverse of what "should" be for consistency. The
           ! other two, the pointers are cyclic consistent: i,j->k,
           ! j,k (wrap) ->i, but for the j-direction is is i,k->j when
           ! to be consistent with the others it should be
           ! k,i->j. Hope that made sense. 

           select case(BCFaceID(mm))
           case(iMin, jMax, kMin)
              regularOrdering = .True.
           case default
              regularOrdering = .False.
           end select

           ! Now this can be reversed *again* if we have a block that
           ! is left handed. 
           if (.not. rightHanded) then 
              regularOrdering = .not. (regularOrdering)
           end if

           ! Start and end bounds for NODES
           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

           ! ni, nj are the number of NODES
           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1

           ! Loop over the faces....this is the node sizes - 1
           if (regularOrdering) then 
              do j=1,nj-1
                 do i=1,ni-1
                    iCell = iCell + 1
                    connLocal(1, iCell) = cumNodeProc(myid) + iNode + (j-1)*ni + i
                    connLocal(2, iCell) = cumNodeProc(myid) + iNode + (j-1)*ni + i + 1
                    connLocal(3, iCell) = cumNodeProc(myid) + iNode + (j)*ni + i + 1 
                    connLocal(4, iCell) = cumNodeProc(myid) + iNode + (j)*ni + i
                    ! Set the cluster
                    clusterCellLocal(iCell) = c
                 end do
              end do
           else
              ! Do the reverse ordering
              do j=1,nj-1
                 do i=1,ni-1
                    iCell = iCell + 1
                    connLocal(1, iCell) = cumNodeProc(myid) + iNode + (j-1)*ni + i
                    connLocal(2, iCell) = cumNodeProc(myid) + iNode + (j  )*ni + i
                    connLocal(3, iCell) = cumNodeProc(myid) + iNode + (j)  *ni + i + 1 
                    connLocal(4, iCell) = cumNodeProc(myid) + iNode + (j-1)*ni + i + 1   

                    ! Set the cluster
                    clusterCellLocal(iCell) = c
                 end do
              end do
           end if
                   
           ! Loop over the nodes
           do j=jBeg,jEnd
              do i=iBeg,iEnd
                 iNode = iNode + 1
                 ! The plus one is for the pointer offset
                 nodesLocal(:, iNode) = xx(i+1, j+1, :)
                 clusterNodeLocal(iNode) = c
                 nodeIndicesLocal(iNode) = ind(i+1, j+1)
              end do
           end do
        end if
     end do
  end do


  ! Allocate space for the global reduced surface
  allocate(nodesGlobal(3, nNodesGlobal), connGlobal(4, nCellsGlobal), &
       clusterCellGlobal(nCellsGlobal), clusterNodeGlobal(nNodesGlobal), &
       nodeIndicesGlobal(nNodesGlobal))
         
  ! Communicate the nodes, connectivity and cluster information to everyone
  call mpi_allgatherv(nodesLocal, 3*nNodesLocal, sumb_real, & 
       nodesGlobal, nNodeProc*3, cumNodeProc*3, sumb_real, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call mpi_allgatherv(clusterNodeLocal, nNodesLocal, sumb_integer, & 
       clusterNodeGlobal, nNodeProc, cumNodeProc, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call mpi_allgatherv(nodeIndicesLocal, nNodesLocal, sumb_integer, & 
       nodeIndicesGlobal, nNodeProc, cumNodeProc, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call mpi_allgatherv(connLocal, 4*nCellsLocal, sumb_integer, &
       connGlobal, nCellProc*4, cumCellProc*4, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call mpi_allgatherv(clusterCellLocal, nCellsLocal, sumb_integer, &
       clusterCellGlobal, nCellProc, cumCellProc, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Free the local data we do not need anymore
  deallocate(nodesLocal, connLocal, clusterCellLocal, clusterNodeLocal, &
       nCellProc, cumCellProc, nNodeProc, cumNodeProc, nodeIndicesLocal)

  ! We will now build separate trees for each cluster. 
  allocate(nodesPerCluster(nClusters), cellsPerCluster(nClusters), &
       cnc(nClusters), ccc(nClusters))
  nodesPerCluster = 0
  cellsPerCluster = 0

  ! Count up the total number of elements and nodes on each cluster
  do i=1, nCellsGlobal
     cellsPerCluster(clusterCellGlobal(i)) = cellsPerCluster(clusterCellGlobal(i)) + 1
  end do
  
  do i=1, nNodesGlobal
     nodesPerCluster(clusterNodeGlobal(i)) = nodesPerCluster(clusterNodeGlobal(i)) + 1
  end do

  ! Create the list of the walls. We are reusing the overset wall derived type here.
  allocate(walls(nClusters))
  allocate(localNodeNums(nNodesGlobal))

  ! Allocate the memory for each of the cluster nodes
  do i=1, nClusters
     nNodes = nodesPerCluster(i)
     nCells = cellsPerCluster(i)
     walls(i)%nCells = nCells
     walls(i)%nNodes = nNodes

     allocate(walls(i)%x(3, nNodes), walls(i)%conn(4, nCells), &
          walls(i)%ind(nNodes))
  end do

  ! We now loop through the master list of nodes and elements and
  !  "push" them back to where they should go. We also keep track of
  !  the local node numbers so that the cluster surcells can update
  !  their own conn.
  localNodeNums = 0
  cnc = 0
  do i=1, nNodesGlobal
     c = clusterNodeGlobal(i) ! Cluter this node belongs to
     cnc(c) = cnc(c) + 1 ! "cluster node count:" the 'nth' node for this  cluster

     walls(c)%x(:, cnc(c))= nodesGlobal(:, i)
     walls(c)%ind(cnc(c)) = nodeIndicesGlobal(i)
     localNodeNums(i) = cnc(c)
     
  end do


  ccc = 0
  do i=1, nCellsGlobal
     c = clusterCellGlobal(i)
     ccc(c) = ccc(c) + 1 ! "Cluster cell count" the 'nth' cell for this cluster
     walls(c)%conn(:, ccc(c)) = connGlobal(:, i)
  end do


  do i=1, nClusters

     nCells = walls(i)%nCells
     nNodes = walls(i)%nNodes 

     ! Fistly we need to update the conn to use our local node ordering. 
     do j=1, nCells
        do k=1, 4
           walls(i)%conn(k, j) = localNodeNums(walls(i)%conn(k, j))
        end do
     end do

     ! Allocate temporary space for doing the point reduction. 
     allocate(uniqueNodes(3, nNodes), link(nNodes))

     call pointReduce(walls(i)%x, nNodes, tol, uniqueNodes, link, nUnique)


            

     ! Update the global indices. Use the returned link
     tmpInd => walls(i)%ind
     allocate(walls(i)%ind(nUnique))
     do j=1, walls(i)%nNodes
        walls(i)%ind(link(j)) = tmpInd(j)
     end do
     deallocate(tmpInd)

     ! Reset the number of nodes to be number of unique nodes
     nNodes = nUnique
     walls(i)%nNodes = nNodes
     
     ! Update the nodes with the unique ones. 
     do j=1, nUnique
        walls(i)%x(:, j) = uniqueNodes(:, j)
     end do

     ! Update conn using the link:
     do j=1, nCells
        do k=1, 4
           walls(i)%conn(k, j) = link(walls(i)%conn(k, j))
        end do
     end do

     ! Unique nodes and link are no longer needed
     deallocate(link, uniqueNodes)

     call buildSerialQuad(nCells, nNodes, walls(i)%x, walls(i)%conn, walls(i)%ADT)
     call buildUniqueNormal(walls(i))
  end do

  if (oversetPresent) then 
     ! Finally build up a "full wall" that is made up of all the cluster
     ! walls. Note that we can reuse the space previously allocated for
     ! the global data, namely, nodes and conn. These will be slightly
     ! larger than necessary becuase of the point reduce. 
     
     fullWall%x => nodesGlobal
     fullWall%conn => connGlobal
     fullWall%ind => nodeIndicesGlobal
     allocate(fullWall%norm(3, nNodesGlobal))

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
     call buildUniqueNormal(fullWall)
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
     conn => fullWall%conn
     nodes => fullWall%x
     norm => fullWall%norm
     
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
                 ! walls(c), up to the wall cutoff. 
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

              ! This is now the (possibly) overlapping surface mesh
              ! case. It is somewhat more complex since we use the
              ! same searches to flag cells that are inside the body. 
              
              coor(4) = wallDistCutoff**2
              intInfo(3) = 0
              call minDistancetreeSearchSinglePoint(fullWall%ADT, coor, &
                   intInfo, uvw, dummy, 0, BB, frontLeaves, frontLeavesNew)
              cellID = intInfo(3)

              if (cellID > 0) then
                 

                 if (uvw(4) > nearWallDist**2 .or. walls(c)%nCells == 0) then 

                    call checkInside()
 
                    ! We found a point within the wallDist cutoff OR the
                    ! cell is from a cluster with no walls, ie a
                    ! background cell. Accept it's wall distance, since
                    ! it guaranteed to be correct.
                    
                    do kk=1,4
                       flowDoms(nn, level, sps)%surfNodeIndices(kk, i, j, k) = &
                            fullWall%ind(fullWall%conn(kk, cellID))
                    end do
                    flowDoms(nn, level, sps)%uv(:, i, j, k) = uvw(1:2)
                    
                 else
                    
                    ! This point is *closer* than the nearWallDist AND
                    ! it has a wall. Search for our own wall. 

                    coor(4) = large
                    call minDistancetreeSearchSinglePoint(walls(c)%ADT, coor, &
                         intInfo2, uvw2, dummy, 0, BB, frontLeaves, frontLeavesNew)
                    cellID2 = intInfo2(3)

                    if (uvw2(4) < nearWallDist**2) then 
                       ! Both are close to the wall. Accept the one from our own wall. 
                       do kk=1,4
                          flowDoms(nn, level, sps)%surfNodeIndices(kk, i, j, k) = &
                               walls(c)%ind(walls(c)%conn(kk, cellID2))
                       end do
                       flowDoms(nn, level, sps)%uv(:, i, j, k) = uvw2(1:2)
                    else
                       ! We have already found a closer point from the
                       ! full wall.  This means we need to check if it
                       ! is inside.

                       call checkInside()

                       ! And save the wall-dist info we already had
                       ! computed from the full wall search
        
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

                 ! HOWEVER, It is possible that this cell is actually
                 ! inside the body. To quickly check, run the ray cast
                 ! algo.

                 call intersectionTreeSearchSinglePoint(fullWall%ADT, coor(1:3), &
                      intInfo(1), BBint, frontLeaves, frontLeavesNew)
              
                 ! If we never found *any* intersections, cannot
                 ! possibly be inside, and there is nothing else to do. 
                 if (intInfo(1) == 0) then 
                    cycle
                 else
                    ! We found a ray cast intersection. Looks like it
                    ! might actually be inside the surface after
                    ! all. Re-run the full distance search with no
                    ! dist cutoff.
                    coor(4) = large
                    call minDistancetreeSearchSinglePoint(fullWall%ADT, coor, &
                         intInfo, uvw, dummy, 0, BB, frontLeaves, frontLeavesNew)
                    cellID = intInfo(3)

                    ! Determine if it is inside:
                    call checkInside()
                    
                 end if
              end if
           end do
        end do
     end do
  end do

  ! Now determine all the node indices this processor needs to get. 
  mm = 0
  allocate(indicesToGet(totalVolumeCells(level)*4), link(totalVolumeCells(level)*4))
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
  call unique(indicesToGet, 4*totalVolumeCells(level), nUnique, link)

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
     call VecCreateMPI(SUMB_COMM_WORLD, 3*totalVolumeNodes(level)*nTimeIntervalsSpectral, &
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
     deallocate(walls(i)%x, walls(i)%norm, walls(i)%conn)
     call destroySerialQuad(walls(i)%ADT)
  end do
  deallocate(walls)

  deallocate(nodesGlobal, connGlobal, clusterCellGlobal, &
       clusterNodeGlobal, localNodeNums)

  if (oversetPresent) then 
     call destroySerialQuad(fullWall%ADT)
     deallocate(fullWall%norm)
  end if

  ! Finally communicate the updated iBlanks
  domainLoop:do nn=1, nDom
     flowDoms(nn, level, sps)%intCommVars(1)%var => &
          flowDoms(nn, level, sps)%iblank(:, :, :)
  end do domainLoop
  
  ! Run the generic integer exchange
  call wHalo1to1IntGeneric(1, level, sps, commPatternCell_2nd, internalCell_2nd)

contains

  subroutine checkInside()

    implicit none

    ! bi-linear shape functions (CCW ordering)
    shp(1) = (one-uvw(1))*(one-uvw(2))
    shp(2) = (    uvw(1))*(one-uvw(2))
    shp(3) = (    uvw(1))*(    uvw(2))
    shp(4) = (one-uvw(1))*(    uvw(2))
                    
    xp = zero
    normal = zero
    do jj=1, 4
       xp = xp + shp(jj)*nodes(:, conn(jj, cellID))
       normal = normal + shp(jj)*norm(:, conn(jj, cellID))
    end do
                    
    ! Compute the dot product of normal with cell center
    ! (stored in coor) with the point on the surface.
    v1 = coor(1:3) - xp
    dp = normal(1)*v1(1) + normal(2)*v1(2) + normal(3)*v1(3)
    
    if (dp < zero) then 
       ! We're inside so blank this cell. Set it to -3 as
       ! being a flood seed. 
       
       iBlank(i, j, k) = -3
    end if
  end subroutine checkInside

end subroutine computeHolesInsideBody


subroutine buildUniqueNormal(wall)

  use overset
  implicit none

  ! Input/Output parameters
  type(oversetWall), intent(inout) :: wall
  
  ! Working
  integer(kind=intType), dimension(:),    allocatable :: link, normCount  
  real(kind=realType),   dimension(:, :), pointer :: nodes, norm
  integer(kind=intType), dimension(:, :), pointer :: conn
  real(kind=realType), dimension(3) :: sss,  v1, v2
  integer(kind=intTYpe) :: i, j
  ! Compute the (averaged) uniqe nodal vectors:
  allocate(wall%norm(3, wall%nNodes), normCount(wall%nNodes))
  nodes => wall%x
  conn => wall%conn
  norm => wall%norm

  norm = zero
  normCount = 0

  do i=1, wall%nCells
        
     ! Compute cross product normal and normize
     v1 = nodes(:, conn(3, i)) -  nodes(:, conn(1, i))
     v2 =  nodes(:, conn(4, i)) -  nodes(:, conn(2, i))
     
     sss(1) = (v1(2)*v2(3) - v1(3)*v2(2))
     sss(2) = (v1(3)*v2(1) - v1(1)*v2(3))
     sss(3) = (v1(1)*v2(2) - v1(2)*v2(1))
     sss = sss / sqrt(sss(1)**2 + sss(2)**2 + sss(3)**2)
     
     ! Add to each of the four nodes and increment the number added
     do j=1, 4
        norm(:, conn(j, i)) = norm(:, conn(j, i)) + sss
        normCount(conn(j, i)) = normCount(conn(j, i)) + 1
     end do
  end do
  
  ! Now just divide by the norm count
  do i=1, wall%nNodes
     norm(:, i) = norm(:, i) / normCount(i)
  end do
  
  ! Node count is no longer needed
  deallocate(normCount)

end subroutine buildUniqueNormal
