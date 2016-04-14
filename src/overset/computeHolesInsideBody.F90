subroutine computeHolesInsideBody(level, sps)

  ! This routine will flag the iBlank values in oBlocks with 0 if the
  ! the cell center falls inside the body. This is a (semi) parallel
  ! implementation: The global surface mesh is communicated to all
  ! processors who then search accordingly. This scalable in terms of
  ! computation but not strictly memory. 

  use adtAPI
  use blockPointers
  use wallDistanceData
  use BCTypes
  use communication
  use inputTimeSpectral
  use overset
  use inputOverset
  implicit none

  ! Input Variables
  integer(kind=intType), intent(in) :: level, sps

  ! Local Variables
  integer(kind=intType) :: i, j, k, l, ii, jj, kk, nn, mm, iNode, iCell, c
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ni, nj, nUnique, cellID
  integer(kind=intType) :: ierr, iDim, nCluster

  ! Data for local surface
  integer(kind=intType) :: nNodes, nCells
  integer(kind=intType) :: nNodesLocal, nCellsLocal
  integer(kind=intType), dimension(:, :), allocatable :: connLocal
  integer(kind=intType), dimension(:), allocatable :: clusterNodeLocal
  integer(kind=intType), dimension(:), allocatable :: clusterCellLocal
  real(kind=realType), dimension(:, :), allocatable :: nodesLocal
  real(kind=realType), dimension(:,:,:), pointer :: xx
  real(kind=realType) :: myDist, otherDist, timeA
  logical :: regularOrdering, nearAWall
  integer(kind=intType), dimension(:), allocatable :: clusters

  ! Overset Walls for storing the surface ADT's
  type(oversetWall), dimension(:), allocatable, target :: walls
  type(oversetWall), target :: fullWall

  ! Data for global surface
  integer(kind=intTYpe) :: nNodesGlobal, nCellsGlobal
  integer(kind=intType), dimension(:, :), allocatable, target :: connGlobal
  real(kind=realType), dimension(:, :), allocatable, target :: nodesGlobal

  integer(kind=intType), dimension(:), allocatable :: nodesPerCluster, cellsPerCluster, cnc, ccc
  integer(kind=intType), dimension(:), allocatable :: clusterNodeGlobal
  integer(kind=intType), dimension(:), allocatable :: clusterCellGlobal
  integer(kind=intType), dimension(:), allocatable :: localNodeNums

  integer(kind=intType), dimension(:),    allocatable :: nCellProc, cumCellProc
  integer(kind=intType), dimension(:),    allocatable :: nNodeProc, cumNodeProc
  real(kind=realType),   dimension(:, :), allocatable :: uniqueNodes
  integer(kind=intType), dimension(:),    allocatable :: link

  ! Pointers for easier readibility
  integer(kind=intType), dimension(:, :), pointer :: conn
  real(kind=realType), dimension(:, :), pointer :: nodes, norm

  ! Data for the ADT
  integer(kind=intType) :: intInfo(3)
  real(kind=realType) :: coor(4), uvw(5) 
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
  allocate(clusters(nDomTotal))
  call determineClusters(clusters, nDomTotal, cumDomProc, nCluster)

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
       clusterCellLocal(nCellsLocal), clusterNodeLocal(NNodesLocal))

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
           case (iMax)
              xx   => x(il,:,:,:)
           case (jMin)
              xx   => x(:,1,:,:)
           case (jMax)
              xx   => x(:,jl,:,:)
           case (kMin)
              xx   => x(:,:,1,:)
           case (kMax)
              xx   => x(:,:,kl,:)
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
              end do
           end do
        end if
     end do
  end do

  ! Allocate space for the global reduced surface
  allocate(nodesGlobal(3, nNodesGlobal), connGlobal(4, nCellsGlobal), &
       clusterCellGlobal(nCellsGlobal), clusterNodeGlobal(nNodesGlobal))
  
  ! Communicate the nodes, connectivity and cluster information to everyone
  call mpi_allgatherv(nodesLocal, 3*nNodesLocal, sumb_real, & 
       nodesGlobal, nNodeProc*3, cumNodeProc*3, sumb_real, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call mpi_allgatherv(clusterNodeLocal, nNodesLocal, sumb_integer, & 
       clusterNodeGlobal, nNodeProc, cumNodeProc, sumb_integer, &
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
       nCellProc, cumCellProc, nNodeProc, cumNodeProc)

  ! We will now build separate trees for each cluster. 
  allocate(nodesPerCluster(nCluster), cellsPerCluster(nCluster), &
       cnc(nCluster), ccc(nCluster))
  nodesPerCluster = 0
  cellsPerCluster = 0

  ! Count up the total number of elements and nodes on each cluster
  do i=1, nCellsGlobal
     cellsPerCluster(clusterCellGlobal(i)) = cellsPerCluster(clusterCellGlobal(i)) + 1
  end do
  
  do i=1, nNodesGlobal
     nodesPerCluster(clusterNodeGlobal(i)) = nodesPerCluster(clusterNodeGlobal(i)) + 1
  end do

  ! Create the list of the walls. We are reusing the overset walls here. 
  allocate(walls(nCluster))
  allocate(localNodeNums(nNodesGlobal))

  ! Allocate the memory for each of the cluster nodes
  do i=1, nCluster
     nNodes = nodesPerCluster(i)
     nCells = cellsPerCluster(i)
     walls(i)%nCells = nCells
     walls(i)%nNodes = nNodes

     allocate(walls(i)%x(3, nNodes), walls(i)%conn(4, nCells))
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
     localNodeNums(i) = cnc(c)
  end do

  ccc = 0
  do i=1, nCellsGlobal
     c = clusterCellGlobal(i)
     ccc(c) = ccc(c) + 1 ! "Cluster cell count" the 'nth' cell for this cluster
     walls(c)%conn(:, ccc(c)) = connGlobal(:, i)
  end do


  do i=1, nCluster

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

  ! Finally build up a "full wall" that is made up of all the cluster
  ! walls. Note that we can reuse the space previously allocated for
  ! the global data, namely, nodes and conn. These will be slightly
  ! larger than necessary becuase of the point reduce. 

  fullWall%x => nodesGlobal
  fullWall%conn => connGlobal
  allocate(fullWall%norm(3, nNodesGlobal))

  nNodes = 0
  nCells = 0
  ii = 0
  do i=1, nCluster

     ! Add in the nodes/elements from this cluster

     do j=1, walls(i)%nNodes
        nNodes = nNodes + 1
        fullWall%x(:, nNodes) = walls(i)%x(:, j)
     end do

     do j=1, walls(i)%nCells
        nCells = nCells + 1
        fullWall%conn(:, nCells) = walls(i)%conn(:, j) + ii
     end do

     ! Increment the node offset
     ii = ii + walls(i)%nNodes
  end do

  ! Finish the setup of the full surcell
  fullWall%nCells = nCells
  fullWall%nNodes = nNodes
  call buildSerialQuad(nCells, nNodes, fullWall%x, fullWall%conn, fullWall%ADT)
  call buildUniqueNormal(fullWall)

  ! Allocate the (pointer) memory that may be resized as necessary for
  ! the singlePoint search routine. 
  allocate(stack(100), BB(20), BBint(20), frontLeaves(25), frontLeavesNew(25))

  ! Now determine the iblank status:
  kk = 0
  do nn=1, nDom

     ! Set the cluster for this block
     c = clusters(cumDomProc(myid) + nn)
     conn => fullWall%conn
     nodes => fullWall%x
     norm => fullWall%norm
     call setPointers(nn, level, sps)

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

              ! Compute our exact wall distance using our own surface
              ! ADT. This can be used for the actual wall distance for
              ! RANS.
              coor(4) = large
              uvw(4) = large
              if (walls(c)%nCells > 0) then 
                 kk = kk + 1
                 call minDistancetreeSearchSinglePoint(walls(c)%ADT, coor, intInfo, &
                      uvw, dummy, 0, BB, frontLeaves, frontLeavesNew)
                 cellID = intInfo(3)
              end if

              myDist = sqrt(uvw(4))
              if (myDist < nearWallDist) then 
                 cycle
              end if

              ! Do the very fast ray cast method-intersection search. 
              call intersectionTreeSearchSinglePoint(fullWall%ADT, coor(1:3), &
                   intInfo(1), BBint, frontLeaves, frontLeavesNew)
              
              ! If we never found *any* intersections again, cannot
              ! possibly be inside.
              if (intInfo(1) == 0) then 
                 cycle
               end if

              ! Otherwise do the full min distance search
              coor(4) = large
              call minDistancetreeSearchSinglePoint(fullWall%ADT, coor, &
                   intInfo, uvw, dummy, 0, BB, frontLeaves, frontLeavesNew)
              otherDist = sqrt(uvw(4))
              cellID = intInfo(3)

              if (otherDist < myDist .or. myDist > nearWallDist) then 
                 
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
              end if
           end do
        end do
     end do
  end do
  print *,'myid kk:', myid, kk
  ! Deallocate all the remaining temporary data
  deallocate(stack, BB, frontLeaves, frontLeavesNew, BBint)

  do i=1, nCluster
     deallocate(walls(i)%x, walls(i)%norm, walls(i)%conn)
     call destroySerialQuad(walls(i)%ADT)
  end do
  deallocate(walls)

  call destroySerialQuad(fullWall%ADT)
  deallocate(nodesGlobal, connGlobal, fullWall%norm, &
       clusterCellGlobal, clusterNodeGlobal, localNodeNums)

  ! Finally communicate the updated iBlanks
  ! Update the iblank info. 
  domainLoop:do nn=1, nDom
     flowDoms(nn, level, sps)%intCommVars(1)%var => &
          flowDoms(nn, level, sps)%iblank(:, :, :)
  end do domainLoop
  
  ! Run the generic integer exchange
  call wHalo1to1IntGeneric(1, level, sps, commPatternCell_2nd, internalCell_2nd)
  print *,' new time:', mpi_wtime()-timeA
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
