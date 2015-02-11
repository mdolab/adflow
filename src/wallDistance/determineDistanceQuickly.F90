#ifndef USE_TAPENADE
subroutine determineWallAssociation(level)

  ! This is part one of the two step modified wall distance
  ! computation. In this stage, we determine the u-v coordinates of
  ! the surface cell that yields the minimim distance to each cell. In
  ! addition to this, a petsc VecScatter is generated that pulls nodes
  ! from the full set of x-coordinates jus the nodes needed on a given
  ! processor to form the cells required by a particular
  ! processor. This is done on a per-level/sps basis. 
  use adtAPI
  use blockPointers
  use wallDistanceData
  use BCTypes
  use communication
  use inputTimeSpectral
  use adjointVars
  implicit none

  ! Subroutine Arguments
  integer(kind=intType), intent(in) :: level

  ! Local Variables
  integer(kind=intType) :: ierr
  integer(kind=intType) :: nNodeViscLocal, nFaceViscLocal
  integer(kind=intType) :: nNodeVisc, nFaceVisc
  integer(kind=intType) :: nCell
  integer(kind=intType) :: i,j,k,l,ii,nn,mm, sps, iNode, iFace, iProc
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ni, nj, nUnique, eID, nID, faceID
  integer(kind=intType), dimension(:), allocatable :: &
       nFaceViscProc, cumFaceViscProc, nNodeViscProc, cumNodeViscProc
  real(kind=realType), dimension(:,:), allocatable :: nodesViscLocal, nodesVisc
  integer(kind=intType), dimension(:,:), allocatable :: connViscLocal, connVisc
  integer(kind=intType), dimension(:), allocatable :: nodeIndicesLocal, nodeIndices
  integer(kind=intTYpe), dimension(:), allocatable :: indicesToGet
  integer(kind=intType) :: totalUniqueFaces
  real(kind=realType), dimension(:,:,:), pointer :: xx
  integer(kind=intType), dimension(:,:), pointer :: ind
  integer(kind=intType), dimension(:), allocatable :: nIndicesProc, cumIndicesProc
  ! Data for the ADT
  character(len=10), parameter :: viscAdt = "ViscousADT"
  integer(kind=intType), dimension(1,1) :: connTria
  real(kind=realType),   dimension(3,2) :: dummy
  integer(kind=intType) :: nTria
  integer(kind=intType), dimension(:), allocatable :: procID, elementID, elemInverse
  real(kind=realType), dimension(:),   allocatable :: dist2
  real(kind=realType), dimension(:,:), allocatable :: coor, uvw
  integer(kind=adtElementType), dimension(:), allocatable :: &
       elementType

  ! The first thing we do is gather all the viscous surface nodes to
  ! each processor such that every processor can make it's own copy of
  ! the complex surface mesh to use to search. Note that this
  ! procedure *DOES NOT SCALE IN MEMORY*...ie eventually the surface
  ! mesh will become too large to store on a single processor,
  ! although this will probably not happen until the sizes get up in
  ! the hundreds of millions of cells. 

  nNodeViscLocal = 0
  nFaceViscLocal = 0

  do nn=1,nDom
     call setPointers(nn, level, 1)
     do mm=1,nBocos
        if(BCType(mm) == NSWallAdiabatic .or. BCType(mm) == NSWallIsothermal) then
           iBeg = bcData(mm)%inBeg
           iEnd = bcData(mm)%inEnd
           jBeg = bcData(mm)%jnBeg
           jEnd = bcData(mm)%jnEnd

           nNodeViscLocal = nNodeViscLocal + &
                (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
           nFaceViscLocal = nFaceViscLocal + & 
                (iEnd - iBeg)*(jEnd - jBeg)
        end if
     end do
  end do

  ! Now communicate these sizes with everyone
  allocate(nFaceViscProc(nProc), cumFaceViscProc(0:nProc), &
       nNodeViscProc(nProc), cumNodeViscProc(0:nProc))

  call mpi_allgather(nFaceViscLocal, 1, sumb_integer, nFaceViscProc, 1, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call mpi_allgather(nNodeViscLocal, 1, sumb_integer, nNodeViscProc, 1, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now make cumulative versions of these
  cumFaceViscProc(0) = 0
  cumNodeViscProc(0) = 0
  do i=1,nProc
     cumFaceViscProc(i) = cumFaceViscProc(i-1) + nFaceViscProc(i)
     cumNodeViscProc(i) = cumNodeViscProc(i-1) + nNodeViscProc(i)
  end do

  ! And save the total number of nodes and faces for reference
  nFaceVisc = cumFaceViscProc(nProc)
  nNodeVisc = cumNodeViscProc(nProc)

  ! Allocate the space for the local nodes, elements and the node
  ! indices
  allocate(nodesViscLocal(3, nNodeViscLocal), &
       connViscLocal(4, nFaceViscLocal),&
       nodeIndicesLocal(nNodeViscLocal))

  ! Now allocate the sizes for the global data
  allocate(nodesVisc(3, nNodeVisc), connVisc(4, nFaceVisc), nodeIndices(nNodeVisc))

  ! Determine the number of cell centers for which the distance
  ! must be computed...on a per-sps basis
  nCell = 0
  do nn=1,nDom
     nCell = nCell + flowDoms(nn,level,1)%nx * flowDoms(nn,level,1)%ny*flowDoms(nn,level,1)%nz
  enddo

  ! Allocate the memory for the arrays needed by the ADT.
  allocate(coor(3, nCell), procID(nCell), elementType(nCell), &
       elementID(nCell), uvw(3, nCell), dist2(nCell), elemInverse(nCell), &
       stat=ierr)
  if(ierr /= 0)                         &
       call terminate("determineDistance", &
       "Memory allocation failure for the variables &
       &needed by the adt.")

  ! Master loop over the spectral instances here
  totalUniqueFaces = 0
  allocate(indicesForSPS(nTimeIntervalsSpectral))
  masterSpectralLoop: do sps=1,nTimeIntervalsSpectral
     iNode = 0
     iFace = 0
     ! Second loop over the walls
     do nn=1,nDom
        call setPointers(nn, level, sps)
        do mm=1,nBocos
           if(BCType(mm) == NSWallAdiabatic .or. BCType(mm) == NSWallIsothermal) then

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

              ! Start and end bounds for NODES
              jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
              iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

              ! ni, nj are the number of NODES
              ni = iEnd - iBeg + 1
              nj = jEnd - jBeg + 1

              ! Loop over the faces....this is the node sizes - 1
              do j=1,nj-1
                 do i=1,ni-1
                    iFace = iFace + 1
                    connViscLocal(1, iFace) = cumNodeViscProc(myid) + iNode + (j-1)*ni + i
                    connViscLocal(2, iFace) = cumNodeViscProc(myid) + iNode + (j-1)*ni + i + 1
                    connViscLocal(3, iFace) = cumNodeViscProc(myid) + iNode + (j)*ni + i + 1 
                    connViscLocal(4, iFace) = cumNodeViscProc(myid) + iNode + (j)*ni + i
                 end do
              end do

              ! Loop over the node
              do j=jBeg,jEnd
                 do i=iBeg,iEnd
                    iNode = iNode + 1
                    ! The plus one is for the pointer offset
                    nodesViscLocal(:, iNode) = xx(i+1, j+1, :)
                    nodeIndicesLocal(iNode) = ind(i+1, j+1)
                 end do
              end do
           end if
        end do
     end do

     ! Communicate the nodes, the indices and the elements
     call mpi_allgatherv(nodesViscLocal, 3*nNodeViscLocal, sumb_real, & 
          nodesVisc, nNodeViscProc*3, cumNodeViscProc*3, sumb_real, &
          sumb_comm_world, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call mpi_allgatherv(nodeIndicesLocal, nNodeViscLocal, sumb_integer, &
          nodeIndices, nNodeViscProc, cumNodeViscProc, sumb_integer, &
          sumb_comm_world, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call mpi_allgatherv(connViscLocal, 4*nFaceViscLocal, sumb_integer, &
          connVisc, nFaceViscProc*4, cumFaceViscProc*4, sumb_integer, &
          sumb_comm_world, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Now we have the a full copy of the (viscous) surface mesh on each
     ! processor. Now we can build the ADT on each processor. A few dummy
     ! arguments need to be defined since we don't have triangles.
     nTria    = 0
     connTria = 0
     dummy    = zero

     call adtBuildSurfaceADT(nTria, nFaceVisc, nNodeVisc,&
          nodesVisc, connTria,  connVisc, &
          dummy,    .false.,   MPI_comm_self, &
          viscAdt)

     ! Fill up the center of the owned cells we need to find the distance of
     mm = 0
     domains: do nn=1,nDom
        call setPointers(nn, level, sps)
        do k=2,kl
           do j=2,jl
              do i=2,il
                 mm = mm + 1
                 ! Compute the coordinates of the cell center 
                 coor(1, mm) = eighth*(x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
                      +         x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
                      +         x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                      +         x(i-1,j,  k,  1) + x(i,j,  k,  1))

                 coor(2, mm) = eighth*(x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
                      +         x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
                      +         x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                      +         x(i-1,j,  k,  2) + x(i,j,  k,  2))

                 coor(3, mm) = eighth*(x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
                      +         x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
                      +         x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                      +         x(i-1,j,  k,  3) + x(i,j,  k,  3))

                 ! Initialize the distance squared, because this is an
                 ! inout argument in the call to adtMinDistanceSearch.

                 dist2(mm) = d2Wall(i,j,k)
              enddo
           enddo
        enddo
     end do domains

     ! Perform the search. As no no interpolations are required,
     ! some dummies are passed.
     call adtMinDistanceSearch(nCell, coor, viscAdt, &
          procID, elementType, elementID, &
          uvw,    dist2, 0_intType, &
          dummy,  dummy)

     ! Destroy ADT since we're done
     call adtDeallocateADTs(viscAdt)

     ! What we are doing here is determinng the unique set of element ID
     ! that the cells on this processor *actually* used. This will be
     ! used to determine which nodes need to communicated to update the
     ! wallDistance quickly. 

     call unique(elementID, nCell, nUnique, elemInverse)

     ! We need to STORE elem Inverse and uv, and set d2wall. elemInverse
     ! now refers to the index of the unique faces.
     mm = 0
     do nn=1,nDom
        call setPointers(nn, level, sps)

        ! Check if elemID and uv are allocated yet.
        if (.not. associated(flowDoms(nn,level,sps)%elemID)) then
           allocate(flowDoms(nn,level,sps)%elemID(2:il,2:jl,2:kl),stat=ierr)
           allocate(flowDoms(nn,level,sps)%uv(2,2:il,2:jl,2:kl),stat=ierr)
        end if

        do k=2,kl
           do j=2,jl
              do i=2,il
                 mm = mm + 1
                 flowDoms(nn,level,sps)%elemID(i,j,k) =  elemInverse(mm) + totalUniqueFaces
                 flowDoms(nn,level,sps)%uv(1,i,j,k)   = uvw(1,mm)
                 flowDoms(nn,level,sps)%uv(2,i,j,k)   = uvw(2,mm)
                 flowDoms(nn,level,sps)%d2wall(i,j,k) = sqrt(dist2(mm))
              end do
           end do
        end do
     end do
     ! Increment totalUniqueFaces such that for the next spectral
     ! instance, the face indices are unique
     totalUniqueFaces = totalUniqueFaces + nUnique

     ! Now we have enough information to determine the set of global
     ! nodes we need to form the nUnique face indices that this processor
     ! needs. We will be slightly inefficient here by expanding the faces
     ! (elements) out such that the number of nodes is 4*3*nFace. This
     ! way, we can index off the nodes just my knowing the face
     ! number. This eliminates the need to redoing a local connectivity
     ! with a reduced set of nodes. 
     allocate(indicesForSPS(sps)%id(4*3*nUnique))
     indicesForSPS(sps)%n = nUnique
     ii = 0
     do k=1, nUnique
        ! Element ID in the global surface: NOTE the -ve sign. This is
        ! becuase the the element numbers come out as -ve from the ADT
        ! search of the point we are looking for isn't "in" the
        ! cell. Of course since we have quad cells, none of the points
        ! we will be looking for are actually "in" the
        ! cell. Therefore, all elementIDs are -ve and we account for
        ! that here.
        eID = -elementID(k)

        ! 4 Nodes on each quad
        do j=1, 4 
           nID = connVisc(j, eID)
           do i=1,3 ! 3 spatial dimensions on each node
              ii = ii + 1
              ! Note that ind and indicesToGet are both 0-based for PETSc
              indicesForSPS(sps)%id(ii) = 3*nodeIndices(nID) + i -1
           end do
        end do
     end do
  end do masterSpectralLoop

  ! We need to communicate the sizes each processor wants such that we know 
  allocate(nIndicesProc(nProc), cumIndicesProc(0:nProc))

  call mpi_allgather(totalUniqueFaces*12, 1, sumb_integer, nIndicesProc, 1, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now make cumulative version
  cumIndicesProc(0) = 0
  do i=1,nProc
     cumIndicesProc(i) = cumIndicesProc(i-1) + nIndicesProc(i)
  end do

  ! Now create the full set of indesToGet, deallocating the sps stuff as we go
  allocate(indicesToGet(totalUniqueFaces*12))

  ii = 0
  do sps=1,nTimeIntervalsSpectral
     do i=1,indicesForSPS(sps)%n*12
        ii = ii + 1
        indicesToGet(ii) =  indicesForSPS(sps)%id(i)
     end do
     deallocate(indicesForSPS(sps)%id)
  end do
  deallocate(indicesForSPS)

  ! Now create the index set for this
  call ISCreateGeneral(sumb_comm_world, 4*3*totalUniqueFaces, indicesToGet, &
       PETSC_COPY_VALUES, IS1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Create the volume vector the nodes will be scatter from 
  call VecCreateMPI(SUMB_COMM_WORLD, 3*nNodesLocal(level)*nTimeIntervalsSpectral, PETSC_DETERMINE, &
       xVolumeVec(level), ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Create the vector for where the nodes will be scattered into. We
  ! will make this a parallel vector, even though we will never
  ! actually use it in the parallel sense. Note that it has no data
  ! storage; we will be vec-placing an array into it before it is used. 

  call VecCreateMPI(SUMB_COMM_WORLD, 4*3*totalUniqueFaces, PETSC_DETERMINE, &
       xSurfVec(level), ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now create an index set for this array, which is just all of them
  call ISCreateStride(SUMB_COMM_WORLD, 4*3*totalUniqueFaces, cumIndicesProc(myid), 1, IS2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! We now have enough information to make a vecScatter that goes from
  ! the global x vector to a localSurface vector.
  call VecScatterCreate(xVolumeVec(level), IS1, xSurfVec(level), IS2, &
       wallScatter(level), ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! And dont' forget to destroy the index sets
  call ISDestroy(IS1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call ISDestroy(IS2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Deallocate all the temporary data
  deallocate(nFaceViscProc, cumFaceViscProc, nNodeViscProc, cumNodeViscProc)
  deallocate(NodesViscLocal, nodesVisc)
  deallocate(connViscLocal, connVisc)
  deallocate(nodeIndicesLocal, nodeIndices)
  deallocate(nIndicesProc, cumIndicesProc)
  deallocate(indicesToGet)
  deallocate(procID, elementID, elemInverse, dist2, coor, uvw, elementType)

end subroutine determineWallAssociation
#endif

subroutine updateWallDistancesQuickly(nn, level, sps)

  ! This is the actual update routine that uses xSurf. It is done on
  ! block-level-sps basis.  This is the used to update the wall
  ! distance. Most importantly, this routine is included in the
  ! reverse mode AD routines, but NOT the forward mode. Since it is
  ! done on a per-block basis, it is assumed that the required block
  ! pointers are already set. 

  use blockPointers
  use wallDistanceData
  implicit none

  ! Subroutine arguments
  integer(kind=intType) :: nn, level, sps

  ! Local Variables
  integer(kind=intType) :: i, j, k, ii, faceID
  real(kind=realType) :: xp(3), xc(3), u, v

#ifdef TAPENADE_FAST
  !$AD II-LOOP
  do ii=0,nx*ny*nz-1
     i = mod(ii, nx) + 2
     j = mod(ii/nx, ny) + 2
     k = ii/(nx*ny) + 2
#else
     do k=2,kl
        do j=2,jl
           do i=2,il
#endif
              ! Extract elemID and u-v position for the association of
              ! this cell:
              faceID = flowDoms(nn,level,sps)%elemID(i,j,k)
              u      = flowDoms(nn,level,sps)%uv(1,i,j,k)
              v      = flowDoms(nn,level,sps)%uv(2,i,j,k)

              ! Now we have the 4 corners, use bi-linear shape
              ! functions o to get target: (CCW ordering remember!)
              
              xp(:) = &
                   (one-u)*(one-v)*xSurf(12*(faceID-1) + 1:12*(faceID-1) + 3) + &
                   (    u)*(one-v)*xSurf(12*(faceID-1) + 4:12*(faceID-1) + 6) + &
                   (    u)*(    v)*xSurf(12*(faceID-1) + 7:12*(faceID-1) + 9) + & 
                   (one-u)*(    v)*xSurf(12*(faceID-1) +10:12*(faceID-1) +12)
              
              ! Get the cell center
              xc(1) = eighth*(x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
                   +         x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
                   +         x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                   +         x(i-1,j,  k,  1) + x(i,j,  k,  1))
              
              xc(2) = eighth*(x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
                   +         x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
                   +         x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                   +         x(i-1,j,  k,  2) + x(i,j,  k,  2))
              
              xc(3) = eighth*(x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
                   +         x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
                   +         x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                   +         x(i-1,j,  k,  3) + x(i,j,  k,  3))
              
              ! Now we have the two points...just take the norm of the
              ! distance between them
              
              d2wall(i,j,k) = sqrt(&
                   (xc(1)-xp(1))**2 + (xc(2)-xp(2))**2 + (xc(3)-xp(3))**2)
#ifdef TAPENADE_FAST
           end do
#else
        enddo
     enddo
  enddo
#endif
  
end subroutine updateWallDistancesQuickly
