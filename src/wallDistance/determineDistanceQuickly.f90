!
!      ******************************************************************
!      *                                                                *
!      * File:          determineDistance.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-03-2003                                      *
!      * Last modified: 02-10-2006                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine determineDistance2(level, sps)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * determineDistance determines the distance from the center      *
  !      * of the cell to the nearest viscous wall for owned cells.       *
  !      *                                                                *
  !      ******************************************************************
  !
  use adtAPI
  use blockPointers
  use communication
  use inputPhysics
  use section
  use viscSurface
  implicit none
  !
  !      Subroutine arguments
  !
  integer(kind=intType), intent(in) :: level, sps
  !
  !      Local parameter, which defines the name of the adt to be create
  !      and used here. Only needed because of the api of the adt library
  !
  character(len=10), parameter :: viscAdt = "ViscousADT"
  !
  !      Local variables.
  !
  integer :: ierr

  integer, dimension(:), allocatable :: procID

  integer(kind=intType) :: nCell, nTria

  integer(kind=intType), dimension(1,1) :: connTria
  real(kind=realType),   dimension(3,2) :: dummy

  integer(kind=intType), dimension(:), allocatable :: elementID
  integer(kind=intType), dimension(:), allocatable :: elemInverse

  real(kind=realType), dimension(:),   allocatable :: dist2
  real(kind=realType), dimension(:,:), allocatable :: coor, uvw

  integer(kind=adtElementType), dimension(:), allocatable :: elementType

  integer(kind=intType) :: nn, mm, ll, ii, jj, i, j, k, idim, faceID
  real(kind=realType), dimension(3) :: xc,xp

  integer(kind=intType), dimension(:), allocatable :: nFaceViscProc
  integer(kind=intType), dimension(:), allocatable :: cumFaceViscProc

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Build the adt of the surface grid. As the api requires to
  ! specify both the quadrilateral and the triangular connectivity,
  ! some dummy variables must be passed.

  nTria    = 0
  connTria = 0
  dummy    = zero

  call adtBuildSurfaceADT(nTria,    nquadVisc, nNodeVisc,       &
       coorVisc, connTria,  connVisc,        &
       dummy,    .false.,   SUmb_comm_world, &
       viscAdt)

  ! Determine the number of cell centers for which the distance
  ! must be computed. Also determine how many of them are periodic.

  nCell    = 0
  do nn=1,nDom
     ll    = flowDoms(nn,level,sps)%nx * flowDoms(nn,level,sps)%ny &
          * flowDoms(nn,level,sps)%nz
     nCell = nCell + ll
  enddo

  ! Allocate the memory for the arrays needed by the ADT.

  allocate(coor(3,nCell), procID(nCell), elementType(nCell), &
       elementID(nCell), uvw(3,nCell), dist2(nCell), elemInverse(nCell), &
       stat=ierr)
  if(ierr /= 0)                         &
       call terminate("determineDistance", &
       "Memory allocation failure for the variables &
       &needed by the adt.")

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Step 1: The search of the original coordinates                 *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Loop over the domains to store the coordinates of the cell
  ! centers of the owned cells. Apply the transformation such that
  ! the sections are aligned.

  mm = 0
  domains: do nn=1,nDom

     ! Set the pointers for this block, store the section id a bit
     ! easier in ll and loop over the cell centers.

     call setPointers(nn, level, sps)
     ll = sectionID

     do k=2,kl
        do j=2,jl
           do i=2,il

              ! Compute the coordinates of the cell center relative
              ! to the rotation center of this section.

              xc(1) = eighth*(x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
                   +         x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
                   +         x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                   +         x(i-1,j,  k,  1) + x(i,j,  k,  1)) &
                   - sections(ll)%rotCenter(1)

              xc(2) = eighth*(x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
                   +         x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
                   +         x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                   +         x(i-1,j,  k,  2) + x(i,j,  k,  2)) &
                   - sections(ll)%rotCenter(2)

              xc(3) = eighth*(x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
                   +         x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
                   +         x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                   +         x(i-1,j,  k,  3) + x(i,j,  k,  3)) &
                   - sections(ll)%rotCenter(3)

              ! Apply the periodic transformation for this section to
              ! align it with other sections and store this coordinate
              ! in the appropriate place in coor.

              mm = mm + 1
              coor(1,mm) = rotMatrixSections(ll,1,1)*xc(1) &
                   + rotMatrixSections(ll,1,2)*xc(2) &
                   + rotMatrixSections(ll,1,3)*xc(3) &
                   + sections(ll)%rotCenter(1)

              coor(2,mm) = rotMatrixSections(ll,2,1)*xc(1) &
                   + rotMatrixSections(ll,2,2)*xc(2) &
                   + rotMatrixSections(ll,2,3)*xc(3) &
                   + sections(ll)%rotCenter(2)

              coor(3,mm) = rotMatrixSections(ll,3,1)*xc(1) &
                   + rotMatrixSections(ll,3,2)*xc(2) &
                   + rotMatrixSections(ll,3,3)*xc(3) &
                   + sections(ll)%rotCenter(3)

              ! Initialize the distance squared, because this is an
              ! inout argument in the call to adtMinDistanceSearch.

              dist2(mm) = d2Wall(i,j,k)

           enddo
        enddo
     enddo
  enddo domains

  ! Perform the search. As no no interpolations are required,
  ! some dummies are passed.

  call adtMinDistanceSearch(nCell,  coor,        viscAdt,      &
       procID, elementType, elementID,    &
       uvw,    dist2,       0_intType, &
       dummy,  dummy)

  
  ! Communicate the sizes of the viscous surface patches on each proc
  allocate(nFaceViscProc(nProc))
  allocate(cumFaceViscProc(0:nProc))
  
  call mpi_allgather(nquadVisc,1,sumb_integer,nFaceViscProc,1,sumb_integer, &
       sumb_comm_world,ierr)
  
  cumFaceViscProc(0) = 0
  do i=1,nProc
     cumFaceViscProc(i) = cumFaceViscProc(i-1) + NFaceViscProc(i)
  end do

  ! Modify elementID by the processor that the it is on. This converts
  ! elementID to a global ordering
  do i=1,nCell
     elementID(i) = -elementID(i) + cumFaceViscProc(procID(i))
  end do

  ! Now generate a sorted, unique list of the element ID's this
  ! processor requires
 
  call unique(elementID, nCell, unique_face_info(level,sps)%n, elemInverse)
 
  ! Now that we know n_unique_face do the additional memory. This must
  ! be called every time since n_unique_face may change and thus
  ! unique_elem_id must be reallocated accordingly 
  call allocExtraWallDistanceMemory(level, sps)

  ! Copy the n_unique_face values on elementID into unique_elem_id as
  ! this will be stored
  do i=1,unique_face_info(level,sps)%n
     unique_face_info(level,sps)%id(i) = elementID(i)
  end do

  ! We need to STORE elem Inverse and uv, and set d2wall
  mm = 1
  do nn=1,nDom
     call setPointers(nn, level, sps)
     do k=2,kl
        do j=2,jl
           do i=2,il
              flowDoms(nn,level,sps)%elemID(i,j,k) = elemInverse(mm)
              flowDoms(nn,level,sps)%uv(1,i,j,k)   = uvw(1,mm)
              flowDoms(nn,level,sps)%uv(2,i,j,k)   = uvw(2,mm)
              flowDoms(nn,level,sps)%d2wall(i,j,k) = sqrt(dist2(mm)) + wallOffset
              mm = mm + 1
           end do
        end do
     end do
  end do

  ! Release the memory of the ADT and the arrays of the module
  ! viscSurface.

  call adtDeallocateADTs(viscAdt)

  deallocate(connVisc, coorVisc, rotMatrixSections, stat=ierr)
  if(ierr /= 0)                         &
       call terminate("determineDistance", &
       "Deallocation error for the arrays &
       &of viscSurface")

  ! Release the variables needed by the ADT.

  deallocate(coor, procID, elementType, elementID, uvw, dist2, &
       elemInverse, stat=ierr)
  if(ierr /= 0)                          &
       call terminate("determineDistance", &
       "Deallocation failure for the variables &
       &needed by the adt.")

end subroutine determineDistance2

subroutine determineDistance3(level, sps)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * determineDistance determines the distance from the center      *
  !      * of the cell to the nearest viscous wall for owned cells.       *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use communication
  use inputPhysics
  use section
  use viscSurface
  implicit none
  !
  !      Subroutine arguments
  !
  integer(kind=intType), intent(in) :: level, sps
  !
  !      Local parameter, which defines the name of the adt to be create
  !      and used here. Only needed because of the api of the adt library
  !
  !      Local variables.
  !
  integer :: ierr
  integer, dimension(mpi_status_size) :: status

  integer(kind=intType) :: nn, mm, ll, ii, jj, i, j, k, idim, index
  integer(kind=intType) :: iProc, faceID

  real(kind=realType), dimension(3) :: xc, xp
  real(kind=realType) :: u,v, newDist

  integer(kind=intType), dimension(:), allocatable :: nFaceViscProc
  integer(kind=intType), dimension(:), allocatable :: cumFaceViscProc

  integer(kind=intType), dimension(:), allocatable :: n_proc_send
  integer(kind=intType), dimension(:), allocatable :: n_proc_recv

  integer(kind=intType), dimension(:), allocatable :: cum_proc_recv
  integer(kind=intType), dimension(:), allocatable :: cum_proc_send

  integer(kind=intType), dimension(:), allocatable :: iRecvBuff
  real(kind=realType), dimension(:), allocatable :: rSendBuff
  real(kind=realType), dimension(:), allocatable :: rRecvBuff

  real(kind=realType), dimension(:,:,:), allocatable :: local_coor
  
  ! Communicate the sizes of the viscous surface patches on each proc
  allocate(nFaceViscProc(nProc))
  allocate(cumFaceViscProc(0:nProc))
  
  call mpi_allgather(nquadVisc,1,sumb_integer,nFaceViscProc,1,sumb_integer, &
       sumb_comm_world,ierr)
  
  cumFaceViscProc(0) = 0
  do i=1,nProc
     cumFaceViscProc(i) = cumFaceViscProc(i-1) + NFaceViscProc(i)
  end do

  ! Allocate arrays that store the size of the sends and receives and zero

  allocate(local_coor(3,4,unique_face_info(level,sps)%n))
  allocate(n_proc_send(nProc),n_proc_recv(nProc))
  allocate(cum_proc_recv(0:nProc),cum_proc_send(0:nProc))
  n_proc_send = 0_intType
  n_proc_recv = 0_intType

  ! Loop over the unique faces and determine the number of elements we
  ! actually need from each processor. 
  
  n_proc_recv(:) = 0_intType
  do i=1,unique_face_info(level,sps)%n

     ! Determine which processor this face this is on:
     procLoop: do iProc=1,nProc
        if (unique_face_info(level,sps)%id(i) > cumFaceViscProc(iProc-1) .and.&
             unique_face_info(level,sps)%id(i) <= cumFaceViscProc(iProc)) then
           exit procLoop 
        end if
     end do procLoop
     n_proc_recv(iProc) = n_proc_recv(iProc) + 1
  end do
 
  ! ------------------------------------------------------------
  !     Step 1: Communicate the sizes of data each processor
  !             will Send and Recieve
  ! ------------------------------------------------------------

  call  MPI_Alltoall(&
       n_proc_recv, 1, sumb_integer, &
       n_proc_send, 1, sumb_integer, &
       sumb_comm_world, ierr)

  ! Create a cumulative form of this as well
  
  cum_proc_recv(0) = 0
  cum_proc_send(0) = 0
  do iProc=1,nProc
     cum_proc_recv(iProc) = cum_proc_recv(iProc-1) + n_proc_recv(iProc)
     cum_proc_send(iProc) = cum_proc_send(iProc-1) + n_proc_send(iProc)
  end do

  ! ------------------------------------------------------------
  !     Step 2: Each processors sends the list of elements it 
  !             wants to every processor
  ! ------------------------------------------------------------

  allocate(iRecvBuff(cum_proc_send(nProc)))
  
  ! Note: we are sending unique_elem_id which has the global number of
  ! viscous faces. The offset is subtracted when we form coordinates
  ! in the next step

  call MPI_Alltoallv(&
       unique_face_info(level,sps)%id, n_proc_recv, cum_proc_recv, &
       sumb_integer, &
       iRecvBuff                     , n_proc_send, cum_proc_send, &
       sumb_integer, &
       sumb_comm_world, ierr)

  ! ------------------------------------------------------------
  !     Step 3: Each processors takes the faces in iRecvBuff and 
  !             assembles the coordinates it must send off
  ! ------------------------------------------------------------

  allocate(rSendBuff(cum_proc_send(nProc)*12))
  
  ii = 1
  do i=1,cum_proc_send(nProc)  ! size of the recvBuff
     index = iRecvBuff(i) - cumFaceViscProc(myid)
     do j=1,4
        do iDim=1,3
           rSendBuff(ii) = coorVisc(iDim,connVisc(j,index))
           ii = ii + 1
        end do
     end do
  end do

  ! ------------------------------------------------------------
  !     Step 4: Each processors sends off the coordiantes
  ! ------------------------------------------------------------
  
  allocate(rRecvBuff(cum_proc_recv(nProc)*12))

  call MPI_Alltoallv(&
       rSendBuff, n_proc_send*12, cum_proc_send*12, sumb_real, &
       rRecvBuff, n_proc_recv*12, cum_proc_recv*12, sumb_real, &
       sumb_comm_world,ierr)

  ! ------------------------------------------------------------
  !     Step 5: Copy recvBuff Each processors sends off the coordiantes
  ! ------------------------------------------------------------

  local_coor = zero
  ii = 1
  do i=1,cum_proc_recv(nProc)
     do j=1,4
        do iDim=1,3
           local_coor(iDim,j,i) = rRecvBuff(ii)
           ii = ii + 1
        end do
     end do
  end do

  ! ------------------------------------------------------------
  !     Step 6: Finally we can update the distances
  ! ------------------------------------------------------------
  
  do nn=1,nDom
     call setPointers(nn,level,sps)
     ll = sectionID
     do k=2,kl
        do j=2,jl
           do i=2,il

              ! Extract elemID and u-v position for the association of
              ! this cell:
              faceID = flowDoms(nn,level,sps)%elemID(i,j,k)
              u      = flowDoms(nn,level,sps)%uv(1,i,j,k)
              v      = flowDoms(nn,level,sps)%uv(2,i,j,k)

              ! Now we have the 4 corners, use bi-linear shape
              ! functions o to get target: (CCW ordering remember!)

              xp(:) = &
                   (one-u)*(one-v)*local_coor(:,1,faceID) + &
                   (    u)*(one-v)*local_coor(:,2,faceID) + &
                   (    u)*(    v)*local_coor(:,3,faceID) + &
                   (one-u)*(    v)*local_coor(:,4,faceID)
       
              ! Get the cell center
              xc(1) = eighth*(x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
                   +         x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
                   +         x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                   +         x(i-1,j,  k,  1) + x(i,j,  k,  1))&
              - sections(ll)%rotCenter(1)

              xc(2) = eighth*(x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
                   +         x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
                   +         x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                   +         x(i-1,j,  k,  2) + x(i,j,  k,  2))&
              - sections(ll)%rotCenter(2)

              xc(3) = eighth*(x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
                   +         x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
                   +         x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                   +         x(i-1,j,  k,  3) + x(i,j,  k,  3))&
              - sections(ll)%rotCenter(3)
         
              ! Now we have the two points...just take the norm of the
              ! distance between them and add possible wallOffset
              
              newDist =  sqrt(&
                   (xc(1)-xp(1))**2 + (xc(2)-xp(2))**2 + (xc(3)-xp(3))**2) + &
                   wallOffset
              
              d2wall(i,j,k) = newDist
           end do
        end do
     end do
  end do

  ! Deallocate all the data we used in this routine.
  deallocate(nFaceViscProc, cumFaceViscProc)
  deallocate(n_proc_send, n_proc_recv)
  deallocate(cum_proc_recv, cum_proc_send)
  deallocate(iRecvBuff, rSendBuff, rRecvBuff)
  deallocate(local_coor)

  ! We ALSO need to deallocate data that was allocated in
  ! visccousSurfaceMesh
  deallocate(connVisc, coorVisc, rotMatrixSections, stat=ierr)
  if(ierr /= 0)                         &
       call terminate("determineDistance", &
       "Deallocation error for the arrays &
       &of viscSurface")

end subroutine determineDistance3
