subroutine buildClusterWalls(level, sps, useDual, walls)

  ! This routine will will build a global reduced surface mesh and ADT
  ! for each cluster. It can build using either the primal mesh or the
  ! dual mesh depending on the  useDual option. 

  use adtBuild, only : buildSerialQuad
  use blockPointers
  use communication
  use inputphysics
  use inputTimeSpectral
  use overset
  use inputOverset
  use utils, only : setPointers, EChk, pointReduce
  use warping, only : getCGNSMeshIndices
  implicit none

  ! Input Variables
  integer(kind=intType), intent(in) :: level, sps
  logical :: useDual
  type(oversetWall), intent(inout), dimension(nClusters), target :: walls

  ! Local Variables
  integer(kind=intType) :: i, j, k, l, ii, jj, kk, nn, mm, iNode, iCell, c
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ni, nj, nUnique, cellID, cellID2
  integer(kind=intType) :: ierr, iDim, lj

  ! Data for local surface
  integer(kind=intType) :: nNodes, nCells
  integer(kind=intType) :: nNodesLocal, nCellsLocal
  integer(kind=intType), dimension(:, :), allocatable :: connLocal
  integer(kind=intType), dimension(:), allocatable :: clusterNodeLocal
  integer(kind=intType), dimension(:), allocatable :: clusterCellLocal
  real(kind=realType), dimension(:, :), allocatable :: nodesLocal
  real(kind=realType), dimension(:,:,:), pointer :: xx, xx1, xx2, xx3, xx4
  integer(kind=intType), dimension(:,:,:), pointer :: globalCGNSNode
  integer(kind=intType), dimension(:,:), pointer :: ind, indCGNS
  integer(kind=intType), dimension(:,:), pointer :: indCell
  logical :: regularOrdering

  ! Data for global surface
  integer(kind=intTYpe) :: nNodesGlobal, nCellsGlobal
  integer(kind=intType), dimension(:, :), allocatable, target :: connGlobal
  real(kind=realType), dimension(:, :), allocatable, target :: nodesGlobal
  integer(kind=intType), dimension(:), allocatable, target :: nodeIndicesGlobal
  integer(kind=intType), dimension(:), allocatable, target :: nodeIndicesCGNSGlobal
  integer(kind=intType), dimension(:), allocatable, target :: cellIndicesGlobal

  integer(kind=intType), dimension(:), allocatable :: nodesPerCluster, cellsPerCluster, cnc, ccc
  integer(kind=intType), dimension(:), allocatable :: clusterNodeGlobal
  integer(kind=intType), dimension(:), allocatable :: clusterCellGlobal
  integer(kind=intType), dimension(:), allocatable :: localNodeNums
  integer(kind=intType), dimension(:), allocatable :: nodeIndicesLocal
  integer(kind=intType), dimension(:), allocatable :: nodeIndicesCGNSLocal
  integer(kind=intType), dimension(:), allocatable :: cellIndicesLocal
  integer(kind=intType), dimension(:), allocatable :: cgnsIndices, curCGNSNode

  integer(kind=intType), dimension(:),    allocatable :: nCellProc, cumCellProc
  integer(kind=intType), dimension(:),    allocatable :: nNodeProc, cumNodeProc
  real(kind=realType),   dimension(:, :), allocatable :: uniqueNodes
  integer(kind=intType), dimension(:),    allocatable :: link
  real(kind=realType), parameter :: tol=1e-12

  ! Pointers for easier readibility
  integer(kind=intType), dimension(:, :), pointer :: conn
  integer(kind=intType), dimension(:), pointer :: tmpInd

  ! The first thing we do is gather all the surface nodes to
  ! each processor such that every processor can make it's own copy of
  ! the complete surface mesh to use to search. Note that this
  ! procedure *DOES NOT SCALE IN MEMORY*...ie eventually the surface
  ! mesh will become too large to store on a single processor,
  ! although this will probably not happen until the sizes get up in
  ! the hundreds of millions of cells. 

  nNodesLocal = 0
  nCellsLocal = 0

  ! Before we start generate a local node indices for the globalCGNS index
  ii = 0
  do nn=1, nDom
     call setPointers(nn, level, sps)
     allocate(flowDoms(nn, level, sps)%globalCGNSNode(1:il, 1:jl, 1:kl))
     flowDoms(nn, level, sps)%globalCGNSNode = 0
     ii = ii + il*jl*kl
  end do

  if (level == 1) then 
     allocate(cgnsIndices(3*ii))
     call getCGNSMeshIndices(size(cgnsIndices), cgnsIndices)
     ii = 0
     do nn=1, nDom
        call setPointers(nn, level, sps)
        do k=1, kl
           do j=1, jl
              do i=1, il
                 ii = ii + 3
                 flowDoms(nn, level, sps)%globalCGNSNode(i,j,k) = cgnsIndices(ii)/3
              end do
           end do
        end do
     end do
     deallocate(cgnsIndices)
  end if

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
           if (useDual) then 
              nNodesLocal = nNodesLocal + &
                   (iEnd - iBeg + 2)*(jEnd - jBeg + 2)
              nCellsLocal = nCellsLocal + & 
                   (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
           else
              nNodesLocal = nNodesLocal + &
                   (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
              nCellsLocal = nCellsLocal + & 
                   (iEnd - iBeg)*(jEnd - jBeg)
           end if
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
       nodeIndicesLocal(nNodesLocal), nodeIndicesCGNSLocal(nNodesLocal), &
       cellIndicesLocal(nCellsLocal))

  iCell = 0
  iNode = 0
  ! Second loop over the local walls
  do nn=1, nDom
     call setPointers(nn, level, sps)
     c = clusters(cumDomProc(myid) + nn)
     globalCGNSNode => flowDoms(nn, level, sps)%globalCGNSNode
     do mm=1,nBocos
        if(  BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal .or. &
             BCType(mm) == EulerWall) then

           jBeg = BCData(mm)%jnBeg-1 ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg-1 ; iEnd = BCData(mm)%inEnd

           if (useDual) then 
              ! For the dual we have to allocate the pointer, xx.
              select case (BCFaceID(mm))
              case (iMin)
                 xx1 => x(1,0:jl,0:kl,:)
                 xx2 => x(1,1:je,0:kl,:)
                 xx3 => x(1,0:jl,1:ke,:)
                 xx4 => x(1,1:je,1:ke,:)
                 ind => globalCell(2, 1:je, 1:ke)
              case (iMax)
                 xx1 => x(il,0:jl,0:kl,:)
                 xx2 => x(il,1:je,0:kl,:)
                 xx3 => x(il,0:jl,1:ke,:)
                 xx4 => x(il,1:je,1:ke,:)
                 ind => globalCell(il, 1:je, 1:ke)
                 
              case (jMin)
                 xx1  => x(0:il,1,0:kl,:)
                 xx2  => x(1:ie,1,0:kl,:)
                 xx3  => x(0:il,1,1:ke,:)
                 xx4  => x(1:ie,1,1:ke,:)
                 ind  => globalCell(1:ie, 2, 1:ke)
                
              case (jMax)
                 xx1 => x(0:il,jl,0:kl,:)
                 xx2 => x(1:ie,jl,0:kl,:)
                 xx3 => x(0:il,jl,1:ke,:)
                 xx4 => x(1:ie,jl,1:ke,:)
                 ind => globalCell(1:ie, jl, 1:ke)
              
              case (kMin)
                 xx1 => x(0:il,0:jl,1,:)
                 xx2 => x(1:ie,0:jl,1,:)
                 xx3 => x(0:il,1:je,1,:)
                 xx4 => x(1:ie,1:je,1,:)
                 ind => globalCell(1:ie, 1:je, 2)

              case (kMax)
                 xx1 => x(0:il,0:jl,kl,:)
                 xx2 => x(1:ie,0:jl,kl,:)
                 xx3 => x(0:il,1:je,kl,:)
                 xx4 => x(1:ie,1:je,kl,:)
                 ind => globalCell(1:ie, 1:je, kl)
              end select

           else
              select case (BCFaceID(mm))
              case (iMin)
                 xx   => x(1,:,:,:)
                 ind  => globalNode(1, :, :)
                 indCGNS => globalCGNSNode(1, :, :)
                 ! Pointer to owned global cell indices
                 indCell => globalCell(2, :, :)
                
              case (iMax)
                 xx   => x(il,:,:,:)
                 ind  => globalNode(il, :, :)
                 indCGNS => globalCGNSNode(il, :, :)
                 
                 ! Pointer to owned global cell indices
                 indCell => globalCell(il, :, :)

              case (jMin)
                 xx   => x(:,1,:,:)
                 ind  => globalNode(:, 1, :)
                 indCGNS => globalCGNSNode(:, 1, :)
                 ! Pointer to owned global cell indices
                 indCell => globalCell(:, 2, :)
                
              case (jMax)
                 xx   => x(:,jl,:,:)
                 ind  => globalNode(:, jl, :)
                 indCGNS => globalCGNSNode(:, jl, :)
                 ! Pointer to owned global cell indices
                 indCell => globalCell(:, jl, :)

              case (kMin)
                 xx   => x(:,:,1,:)
                 ind  => globalNode(:, :, 1)
                 indCGNS => globalCGNSNode(:, :, 1)
                 ! Pointer to owned global cell indices
                 indCell => globalCell(:, :, 2)
                
              case (kMax)
                 xx   => x(:,:,kl,:)
                 ind  => globalNode(:, :, kl)
                 indCGNS => globalCGNSNode(:, :, kl)
                 
                 ! Pointer to owned global cell indices
                 indCell => globalCell(:, :, kl)
                
              end select
              
              ! Just set hte 4 other pointers to xx so we can use the
              ! same quarter-summation code below:
              xx1 => xx
              xx2 => xx
              xx3 => xx
              xx4 => xx

           end if
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

           if (useDual) then 
              ! Start and end bounds for NODES
              jBeg = BCData(mm)%jnBeg-1; jEnd = BCData(mm)%jnEnd
              iBeg = BCData(mm)%inBeg-1; iEnd = BCData(mm)%inEnd
           else
              ! Start and end bounds for NODES
              jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
              iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
           end if

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
        
                    ! Save the global cell index
                    if (useDual) then
                       cellIndicesLocal(iCell) = 0
                    else
                       ! Valid only when using primary nodes
                       cellIndicesLocal(iCell) = indCell(iBeg+i+1, jBeg+j+1)
                    end if
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

                    ! Save the global cell index
                    if (useDual) then
                       cellIndicesLocal(iCell) = 0
                    else
                       ! Valid only when using primary nodes
                       cellIndicesLocal(iCell) = indCell(iBeg+i+1, jBeg+j+1)
                    end if
                 end do
              end do
           end if
                   
           ! Loop over the nodes
           do j=jBeg, jEnd
              do i=iBeg, iEnd
                 iNode = iNode + 1
                 ! The plus one is for the pointer offset
                 nodesLocal(:, iNode) = fourth*(&
                      xx1(i+1, j+1, :) + xx2(i+1, j+1, :) + &
                      xx3(i+1, j+1, :) + xx4(i+1, j+1, :))

                 clusterNodeLocal(iNode) = c
                 nodeIndicesLocal(iNode) = ind(i+1, j+1) ! +1 for pointer offset
                 nodeIndicesCGNSLocal(iNode) = indCGNS(i, j) ! No pointer offset
              end do
           end do
        end if
     end do
  end do

  ! Allocate space for the global reduced surface
  allocate(nodesGlobal(3, nNodesGlobal), connGlobal(4, nCellsGlobal), &
       clusterCellGlobal(nCellsGlobal), clusterNodeGlobal(nNodesGlobal), &
       nodeIndicesGlobal(nNodesGlobal), nodeIndicesCGNSGlobal(nNodesGlobal), &
       cellIndicesGlobal(nCellsGlobal))
         
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

  call mpi_allgatherv(nodeIndicesCGNSLocal, nNodesLocal, sumb_integer, & 
       nodeIndicesCGNSGlobal, nNodeProc, cumNodeProc, sumb_integer, &
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

  call mpi_allgatherv(cellIndicesLocal, nCellsLocal, sumb_integer, &
       cellIndicesGlobal, nCellProc, cumCellProc, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Free the local data we do not need anymore
  deallocate(nodesLocal, connLocal, clusterCellLocal, clusterNodeLocal, &
       nCellProc, cumCellProc, nNodeProc, cumNodeProc, nodeIndicesLocal, &
       nodeIndicesCGNSLocal, cellIndicesLocal)

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
  allocate(localNodeNums(nNodesGlobal))

  ! Allocate the memory for each of the cluster nodes
  do i=1, nClusters
     nNodes = nodesPerCluster(i)
     nCells = cellsPerCluster(i)
     walls(i)%nCells = nCells
     walls(i)%nNodes = nNodes

     allocate(walls(i)%x(3, nNodes), walls(i)%conn(4, nCells), &
          walls(i)%ind(nNodes))
     allocate(walls(i)%indCell(nCells))
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
     
     walls(c)%indCell(ccc(c)) = cellIndicesGlobal(i)
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
     allocate(curCGNSNode(nUnique))
     walls(i)%ind = -1
     curCGNSNode = -1
     do j=1, walls(i)%nNodes
        ! Insted of blinding setting the index, we use the the
        ! nodeIndicesCGNSGlobal to only set the the globalNode with
        ! the smallest CGNS index. This guarantees that the same node
        ! ID is always selected independent of the block
        ! partitioning/splitting. Note that this will work even for
        ! the coarse levels when nodeIndicesCGNSGLobal are all 0's. In
        ! that case the first time wall(i)%ind(link(j)) is touched,
        ! that index is taken. 
        lj = link(j)
        if (walls(i)%ind(lj) == -1 .or. & ! Not set yet
             nodeIndicesCGNSGlobal(j) < curCGNSNode(jl)) then  
           ! OR then potential gloabl CGNS node index is LOWER
           ! than the one I already have

           walls(i)%ind(lj) = tmpInd(j)
           curCGNSNode(lj) = nodeIndicesCGNSGlobal(j)
        end if
     end do
     deallocate(tmpInd, curCGNSNode)

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
  end do

  ! Clean up memeory
  deallocate(nodesGlobal, connGlobal, clusterCellGlobal, &
       clusterNodeGlobal, localNodeNums, nodeIndicesGlobal, &
       nodeIndicesCGNSGlobal)

  do nn=1, nDom
     deallocate(flowDoms(nn, level, sps)%globalCGNSNode)
  end do

end subroutine buildClusterWalls
