subroutine computeHolesInsideBody

  ! This routine will flag the iBlank values in oBlocks with 0 if the
  ! the cell center falls inside the body. Currently only a serial implementation
  use adtAPI
  use blockPointers
  use wallDistanceData
  use BCTypes
  use communication
  use inputTimeSpectral
  use overset
  implicit none

  ! Local Variables
  integer(kind=intType) :: ierr
  integer(kind=intType) :: nNodeVisc, nFaceVisc
  integer(kind=intType) :: nCell
  integer(kind=intType) :: i,j,k,l,ii,nn,mm, ll,sps, level, iNode, iFace
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ni, nj, nUnique, faceID
  integer(kind=intType), dimension(:,:), allocatable ::connVisc

  real(kind=realType), dimension(:,:,:), pointer :: xx
  integer(kind=intType), dimension(:,:), pointer :: ind
  logical :: regularOrdering
  ! Data for the ADT
  character(len=10), parameter :: viscAdt = "ViscousADT"
  integer(kind=intType), dimension(1,1) :: connTria
  real(kind=realType),   dimension(3,2) :: dummy
  integer(kind=intType) :: nTria
  integer(kind=intType), dimension(:), allocatable :: procID, elementID, elemInverse
  real(kind=realType), dimension(:),   allocatable :: dist2
  real(kind=realType), dimension(:,:), allocatable :: coor, uvw
  integer(kind=adtElementType), dimension(:), allocatable :: elementType
  real(kind=realType), parameter :: tol=1e-8
  integer(kind=intType), dimension(:), allocatable :: link, normCount
  real(kind=realType), dimension(:, :), allocatable :: uniqueNodes
  real(kind=realType),dimension(3) :: sss, pt1, pt2, pt3, pt4, v1, v2, v3, v4, xp, xc, normal
  real(kind=realType) :: dp, u, v 
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

  ! Compute the surface mesh
  nNodeVisc = 0
  nFaceVisc = 0
  level = 1
  do nn=1,nDom
     call setPointers(nn, level, 1)
     do mm=1,nBocos
        if(BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal .or. &
             BCType(mm) == EulerWall) then
           iBeg = bcData(mm)%inBeg
           iEnd = bcData(mm)%inEnd
           jBeg = bcData(mm)%jnBeg
           jEnd = bcData(mm)%jnEnd

           nNodeVisc = nNodeVisc + &
                (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
           nFaceVisc = nFaceVisc + & 
                (iEnd - iBeg)*(jEnd - jBeg)
        end if
     end do
  end do

  ! Now allocate the sizes for the global data
  allocate(nodesVisc(3, nNodeVisc), connVisc(4, nFaceVisc))

  ! Determine the number of cell centers for which the distance
  ! must be computed...on a per-sps basis
  nCell = 0
  do nn=1,nDom
     nCell = nCell + flowDoms(nn,level,1)%ie * flowDoms(nn,level,1)%je * flowDoms(nn,level,1)%ke
  enddo

  ! Allocate the memory for the arrays needed by the ADT.
  allocate(coor(3, nCell), procID(nCell), elementType(nCell), &
       elementID(nCell), uvw(3, nCell), dist2(nCell), elemInverse(nCell), &
       stat=ierr)

  ! Hard code these for now
  sps = 1
  level = 1

  iNode = 0
  iFace = 0
  ! Second loop over the walls
  do nn=1,nDom
     call setPointers(nn, level, sps)
     do mm=1,nBocos
        if(BCType(mm) == NSWallAdiabatic .or. &
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
                    connVisc(1, iFace) =  iNode + (j-1)*ni + i
                    connVisc(2, iFace) =  iNode + (j-1)*ni + i + 1
                    connVisc(3, iFace) =  iNode + (j)*ni + i + 1 
                    connVisc(4, iFace) =  iNode + (j)*ni + i
                 end do
              end do
           else
              ! Do the reverse ordering
              do j=1,nj-1
                 do i=1,ni-1
                    iFace = iFace + 1
                    connVisc(1, iFace) = iNode + (j-1)*ni + i
                    connVisc(2, iFace) = iNode + (j  )*ni + i
                    connVisc(3, iFace) = iNode + (j)  *ni + i + 1 
                    connVisc(4, iFace) = iNode + (j-1)*ni + i + 1
                 end do
              end do
           end if
           ! Loop over the nodes
           do j=jBeg,jEnd
              do i=iBeg,iEnd
                 iNode = iNode + 1
                 ! The plus one is for the pointer offset
                 nodesVisc(:, iNode) = xx(i+1, j+1, :)

              end do
           end do
        end if
     end do
  end do

  ! Maximum space for the unique coordinates and link array
  allocate(uniqueNodes(3, nNodeVisc))
  allocate(link(nNodeVisc))

  call pointReduce(nodesVisc, nNodeVisc, tol, uniqueNodes, link, nUnique)

  ! Reset nNodeVisc to be nUnique
  nNodeVisc = nUnique

  ! Overwrite nodesVisc with the uniqueNodes
  do i=1,nUnique
     nodesVisc(:, i) = uniqueNodes(:, i)
  end do

  ! Update connVisc using the link:
  do i=1,nFaceVisc
     do j=1,4
        connVisc(j, i) = link(connVisc(j, i))
     end do
  end do

  ! No longer need the unique data
  deallocate(uniqueNodes, link)

  ! Compute the uniqe nodal vectors:
  allocate(normVisc(3, nNodeVisc), normCount(nNodeVisc))
  normVisc = zero
  normCount = 0

  do i=1, nFaceVisc
     ! Extract the four nodes for this face
     pt1 = nodesVisc(:, connVisc(1, i))
     pt2 = nodesVisc(:, connVisc(2, i))
     pt3 = nodesVisc(:, connVisc(3, i))
     pt4 = nodesVisc(:, connVisc(4, i))

     ! Compute cross product normal and normize
     v1 = pt3 - pt1
     v2 = pt4 - pt2
     sss(1) = (v1(2)*v2(3) - v1(3)*v2(2))
     sss(2) = (v1(3)*v2(1) - v1(1)*v2(3))
     sss(3) = (v1(1)*v2(2) - v1(2)*v2(1))
     sss = sss / sqrt(sss(1)**2 + sss(2)**2 + sss(3)**2)

     ! Add to each of the four nodes and increment the number added
     do j=1,4
        normVisc(:, connVisc(j, i)) = normVisc(:, connVisc(j, i)) + sss
        normCount(connVisc(j, i)) = normCount(connVisc(j, i)) + 1
     end do
  end do

  ! Now just divide by the nomr count
  do i=1,nNodeVisc
     normVisc(:, i) = normVisc(:, i) / normCount(i)
  end do

  ! Dummy triangular data for the ADT
  nTria    = 0
  connTria = 0
  dummy    = zero

  call adtBuildSurfaceADT(nTria, nFaceVisc, nNodeVisc, &
       nodesVisc, connTria,  connVisc,  dummy,    .false., &
       MPI_comm_self, viscAdt)

  ! Initialize the distance squared, because this is an
  ! inout argument in the call to adtMinDistanceSearch.
  dist2 = 1e30
    
  ! Fill up the center of the owned cells we need to find the distance of
  mm = 0
  domains: do nn=1,nDom
     ll = 0
     do k=1,oBlocks(nn)%ke
        do j=1,oBlocks(nn)%je
           do i=1,oBlocks(nn)%ie
              mm = mm + 1
              ll = ll + 1
              ! Compute the coordinates of the cell center 
              coor(:, mm) = oBlocks(nn)%xDual(:, ll)
              
           enddo
        enddo
     enddo
  end do domains

  ! Perform the search. As no interpolations are required,
  ! some dummies are passed.
  call adtMinDistanceSearch(nCell, coor, viscAdt, &
       procID, elementType, elementID, &
       uvw,    dist2, 0_intType, &
       dummy,  dummy)

  ! Destroy ADT since we're done
  call adtDeallocateADTs(viscAdt)

  ! Now determine the iblank status:
  mm = 0
  do nn=1, nDom
     call setPointers(nn, level, sps)
     ll = 0
     do k=1,ke
        do j=1,je
           do i=1,ie
              mm = mm + 1
              ll = ll + 1
              if (elementID(mm) < 0) then 
                 faceID = -elementID(mm) ! Note the negative 
              else
                 ! In the rare case that search node is actually exactly on
                 ! the surface, elementID will be positive. 
                 faceID = elementID(mm)
              end if

              u = uvw(1, mm)
              v = uvw(2, mm)

              ! Extract the four corners and four vectors
              pt1 = nodesVisc(:, connVisc(1, faceID))
              pt2 = nodesVisc(:, connVisc(2, faceID))
              pt3 = nodesVisc(:, connVisc(3, faceID))
              pt4 = nodesVisc(:, connVisc(4, faceID))

              v1 = normVisc(:, connVisc(1, faceID))
              v2 = normVisc(:, connVisc(2, faceID))
              v3 = normVisc(:, connVisc(3, faceID))
              v4 = normVisc(:, connVisc(4, faceID))

              ! Bilinearly interpoolate the location on the face as
              ! well as the normal. 
              xp(:) = &
                   (one-u)*(one-v)*pt1 + &
                   (    u)*(one-v)*pt2 + & 
                   (    u)*(    v)*pt3 + &
                   (one-u)*(    v)*pt4 

              normal = &
                   (one-u)*(one-v)*v1 + &
                   (    u)*(one-v)*v2 + & 
                   (    u)*(    v)*v3 + &
                   (one-u)*(    v)*v4 

              ! Get the cell center
              xc = oBlocks(nn)%xDual(:, ll)

              ! Compute the dot product of normal with xc - xp
              v1 = xc - xp
              dp = normal(1)*v1(1) + normal(2)*v1(2) + normal(3)*v1(3)

              if (dp < zero) then 
                 ! We're inside!
                 oBlocks(nn)%iBlank(i, j, k) = 0
              end if
           end do
        end do
     end do
  end do

  ! Deallocate all the temporary data
  deallocate(connVisc, procID, elementID, elemInverse, dist2, coor, uvw, elementType, normCount)

end subroutine computeHolesInsideBody
