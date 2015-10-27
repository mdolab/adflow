subroutine computeHolesInsideBody

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
  implicit none

  ! Local Variables
  integer(kind=intType) :: i, j, k, l, ii, jj, nn, mm, level, sps, iNode, iFace
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ni, nj, nUnique, faceID
  integer(kind=intType) :: ierr, iDim

  ! Data for local surface
  integer(kind=intType) :: nNodeLocal, nFaceLocal
  integer(kind=intType), dimension(:, :), allocatable :: connLocal
  real(kind=realType), dimension(:, :), allocatable :: nodesLocal
  logical :: regularOrdering
  real(kind=realType), dimension(:,:,:), pointer :: xx

  ! Data for global surface
  integer(kind=intTYpe) :: nNode, nFace
  integer(kind=intType), dimension(:),    allocatable :: nFaceProc, cumFaceProc
  integer(kind=intType), dimension(:),    allocatable :: nNodeProc, cumNodeProc
  integer(kind=intType), dimension(:, :), allocatable :: conn
  real(kind=realType),   dimension(:, :), allocatable :: nodes
  real(kind=realType),   dimension(:, :), allocatable :: uniqueNodes, norm
  integer(kind=intType), dimension(:),    allocatable :: link, normCount
  real(kind=realType), dimension(3) :: xMin, xMax

  ! Data for the ADT
  character(len=10), parameter :: adtName = "holeCutADT"
  integer(kind=intType), dimension(3, 0) :: connTria
  real(kind=realType), dimension(0, 0) :: dummy
  logical :: useBBox 
  integer(kind=intType) :: nTria, intInfo(3), jjADT, nAlloc, nSearch
  real(kind=realType) :: coor(4), uvw(5) 
  real(kind=realType), parameter :: tol=1e-12
  integer(kind=intType), dimension(:), pointer :: frontLeaves, frontLeavesNew, BBint
  type(adtBBoxTargetType), dimension(:), pointer :: BB
  
  ! Misc
  real(kind=realType), dimension(3) :: sss, xp, normal, v1, v2
  real(kind=realType) :: dp, shp(4)

  interface
     subroutine pointReduce(pts, N, tol, uniquePts, link, nUnique)
       use precision
       implicit none
       real(kind=realType), dimension(:, :) :: pts
       integer(kind=intType), intent(in) :: N
       real(kind=realType), intent(in) :: tol
       real(kind=realType), dimension(:, :) :: uniquePts
       integer(kind=intType), dimension(:) :: link
       integer(kind=intType) :: nUnique
     end subroutine pointReduce
  end interface

  ! We explicitly use level=1, and sps=1 here. This should change at some point
  level = 1
  sps = 1

  ! The first thing we do is gather all the surface nodes to
  ! each processor such that every processor can make it's own copy of
  ! the complex surface mesh to use to search. Note that this
  ! procedure *DOES NOT SCALE IN MEMORY*...ie eventually the surface
  ! mesh will become too large to store on a single processor,
  ! although this will probably not happen until the sizes get up in
  ! the hundreds of millions of cells. 

  nNodeLocal = 0
  nFaceLocal = 0

  do nn=1,nDom
     call setPointers(nn, level, sps)
     do mm=1,nBocos
        if(BCType(mm) == NSWallAdiabatic .or. &
           BCType(mm) == NSWallIsothermal .or. &
           BCType(mm) == EulerWall) then
           iBeg = bcData(mm)%inBeg
           iEnd = bcData(mm)%inEnd
           jBeg = bcData(mm)%jnBeg
           jEnd = bcData(mm)%jnEnd

           nNodeLocal = nNodeLocal + &
                (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
           nFaceLocal = nFaceLocal + & 
                (iEnd - iBeg)*(jEnd - jBeg)
        end if
     end do
  end do

  ! Now communicate these sizes with everyone
  allocate(nFaceProc(nProc), cumFaceProc(0:nProc), &
           nNodeProc(nProc), cumNodeProc(0:nProc))

  call mpi_allgather(nFaceLocal, 1, sumb_integer, nFaceProc, 1, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call mpi_allgather(nNodeLocal, 1, sumb_integer, nNodeProc, 1, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now make cumulative versions of these
  cumFaceProc(0) = 0
  cumNodeProc(0) = 0
  do i=1,nProc
     cumFaceProc(i) = cumFaceProc(i-1) + nFaceProc(i)
     cumNodeProc(i) = cumNodeProc(i-1) + nNodeProc(i)
  end do

  ! And save the total number of nodes and faces for reference
  nFace = cumFaceProc(nProc)
  nNode = cumNodeProc(nProc)

  ! Allocate the space for the local nodes and element connectivity
  allocate(nodesLocal(3, nNodeLocal), connLocal(4, nFaceLocal))

  iFace = 0
  iNode = 0
  ! Second loop over the local walls
  do nn=1,nDom
     call setPointers(nn, level, sps)
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
                    iFace = iFace + 1
                    connLocal(1, iFace) = cumNodeProc(myid) + iNode + (j-1)*ni + i
                    connLocal(2, iFace) = cumNodeProc(myid) + iNode + (j-1)*ni + i + 1
                    connLocal(3, iFace) = cumNodeProc(myid) + iNode + (j)*ni + i + 1 
                    connLocal(4, iFace) = cumNodeProc(myid) + iNode + (j)*ni + i
                 end do
              end do
           else
              ! Do the reverse ordering
              do j=1,nj-1
                 do i=1,ni-1
                    iFace = iFace + 1
                    connLocal(1, iFace) = cumNodeProc(myid) + iNode + (j-1)*ni + i
                    connLocal(2, iFace) = cumNodeProc(myid) + iNode + (j  )*ni + i
                    connLocal(3, iFace) = cumNodeProc(myid) + iNode + (j)  *ni + i + 1 
                    connLocal(4, iFace) = cumNodeProc(myid) + iNode + (j-1)*ni + i + 1
                 end do
              end do
           end if
           ! Loop over the nodes
           do j=jBeg,jEnd
              do i=iBeg,iEnd
                 iNode = iNode + 1
                 ! The plus one is for the pointer offset
                 nodesLocal(:, iNode) = xx(i+1, j+1, :)

              end do
           end do
        end if
     end do
  end do

  ! Allocate space for the global reduced surface
  allocate(nodes(3, nNode), conn(4, nFace))

  ! Communicate the nodes and connectivity to everyone
  call mpi_allgatherv(nodesLocal, 3*nNodeLocal, sumb_real, & 
       nodes, nNodeProc*3, cumNodeProc*3, sumb_real, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call mpi_allgatherv(connLocal, 4*nFaceLocal, sumb_integer, &
       conn, nFaceProc*4, cumFaceProc*4, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Free the local data we do not need anymore
  deallocate(nodesLocal, connLocal)

  ! Before we can build the ADTTree we need to compute a unique set of
  ! surface nodes. 

  ! Maximum space for the unique coordinates and link array
  allocate(uniqueNodes(3, nNode))
  allocate(link(nNode))

  call pointReduce(nodes, nNode, tol, uniqueNodes, link, nUnique)

  ! Reset nNode to be nUnique
  nNode = nUnique

  ! Overwrite nodes with the uniqueNodes. Compute the max/min while
  ! we're accessing the memory
  xMin = large
  xMax = eps

  do i=1, nUnique
     nodes(:, i) = uniqueNodes(:, i)

     do iDim=1,3
        xMin(iDim) = min(xMin(iDim), uniqueNodes(iDim, i))
        xMax(iDim) = max(xMax(iDim), uniqueNodes(iDim, i))
     end do
  end do

  ! Update conn using the link:
  do i=1, nFace
     do j=1, 4
        conn(j, i) = link(conn(j, i))
     end do
  end do

  ! No longer need the unique data
  deallocate(uniqueNodes, link)

  ! Compute the (averaged) uniqe nodal vectors:
  allocate(norm(3, nNode), normCount(nNode))
  norm = zero
  normCount = 0

  do i=1, nFace

     ! Compute cross product normal and normize
     v1 = nodes(:, conn(3, i)) -  nodes(:, conn(1, i))
     v2 =  nodes(:, conn(4, i)) -  nodes(:, conn(2, i))

     sss(1) = (v1(2)*v2(3) - v1(3)*v2(2))
     sss(2) = (v1(3)*v2(1) - v1(1)*v2(3))
     sss(3) = (v1(1)*v2(2) - v1(2)*v2(1))
     sss = sss / sqrt(sss(1)**2 + sss(2)**2 + sss(3)**2)

     ! Add to each of the four nodes and increment the number added
     do j=1,4
        norm(:, conn(j, i)) = norm(:, conn(j, i)) + sss
        normCount(conn(j, i)) = normCount(conn(j, i)) + 1
     end do
  end do

  ! Now just divide by the norm count
  do i=1, nNode
     norm(:, i) = norm(:, i) / normCount(i)
  end do

  ! Node count is no longer needed
  deallocate(normCount)

  ! Dummy triangular data for the ADT
  nTria    = 0
  useBBox = .False.

  call adtBuildSurfaceADT(nTria, nFace, nNode, nodes, connTria,  conn,  &
       dummy, useBBox, MPI_comm_self, adtName)
  
  nAlloc = ubound(ADTs, 1)
  do jjAdt=1,nAlloc
     if(adtName == ADTs(jjAdt)%adtID) exit
  enddo
  
  ! Allocate the (pointer) memory that may be resized as necessary for
  ! the singlePoint search routine. 
  allocate(stack(100), BB(20), BBint(20), frontLeaves(25), frontLeavesNew(25))
   
  ! Now determine the iblank status:
  nSearch = 0
  do nn=1, nDom
  
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

              ! We can do a bounding box check here...if our query
              ! point is outside the entire surface bouding box,
              ! there's no way it can be inside. 
              if ( coor(1) > xMax(1) .or. &
                   coor(1) < xMin(1) .or. &
                   coor(2) > xMax(2) .or. &
                   coor(2) < xMin(2) .or. &
                   coor(3) > xMax(3) .or. &
                   coor(3) < xMin(3)) then 
                 cycle
              end if

              ! ! Do the very fast ray cast method-intersection search. 
              call intersectionTreeSearchSinglePoint(ADTs(jjAdt), coor(1:3), &
                   intInfo(1), BBint, frontLeaves, frontLeavesNew)
              
              ! If we never found *any* intersections again, cannot
              ! possibly be inside.
              if (intInfo(1) == 0) then 
                  cycle
               end if

              nSearch = nSearch + 1
              ! Otherwise do the full min distance search
              coor(4) = 1e30
              call minDistancetreeSearchSinglePoint(ADTs(jjAdt), coor, intInfo, &
                   uvw, dummy, 0, BB, frontLeaves, frontLeavesNew)
              
              faceID = intInfo(3)

              ! bi-linear shape functions (CCW ordering)
              shp(1) = (one-uvw(1))*(one-uvw(2))
              shp(2) = (    uvw(1))*(one-uvw(2))
              shp(3) = (    uvw(1))*(    uvw(2))
              shp(4) = (one-uvw(1))*(    uvw(2))

              xp = zero
              normal = zero
              do jj=1, 4
                 xp = xp + shp(jj)*nodes(:, conn(jj, faceID))
                 normal = normal + shp(jj)*norm(:, conn(jj, faceID))
              end do

              ! Compute the dot product of normal with cell center
              ! (stored in coor) with the point on the surface.
              v1 = coor(1:3) - xp
              dp = normal(1)*v1(1) + normal(2)*v1(2) + normal(3)*v1(3)

              if (dp < zero) then 
                 ! We're inside so blank this cell. 
                 iBlank(i, j, k) = 0
              end if
           end do
        end do
     end do
  end do

  ! Destroy ADT since we're done
  call adtDeallocateADTs(adtName)

  ! Deallocate all the remaining temporary data
  deallocate(nodes, conn, norm)
  deallocate(stack, BB, frontLeaves, frontLeavesNew, BBint)

  ! Finally communicate the updated iBlanks
  call exchangeIBlanks(level, sps, commPatternCell_2nd, internalCell_2nd)

end subroutine computeHolesInsideBody
