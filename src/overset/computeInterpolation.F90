!
!      ******************************************************************
!      *                                                                *
!      * computeOversetInterpolation is the top level routine that      *
!      * implements the implicit hole cutting method for determing      *
!      * overset grid connectivitiies. For now it only operates on the  *
!      * finest grid and first spectral instance. It also only works    *
!      * with a single processor                                        *
!      *                                                                *
!      ******************************************************************

subroutine computeOversetInterpolation

  use constants
  use communication
  use blockPointers
  use BCTypes
  use overset
  use adtAPI
  implicit none

  ! Local Variables
  integer(kind=intType) ::i, j, k,  level, sps, iil, jjl, iDom, jDom, iDim, nn, mm, tmp
  real(kind=realType), dimension(3, nDom) :: xMinLocal
  integer(kind=intType), dimension(3, nDom) :: localDim
  real(kind=realType), dimension(:, :), allocatable :: xMin, xMax
  integer(kind=intType), dimension(:), allocatable :: nDomProc, cumDomProc
  integer(kind=intType), dimension(:, :), allocatable :: dims
  integer(kind=intType), dimension(:, :, :), allocatable :: iBlankTmp
  logical ,dimension(:, :), allocatable :: iOverlap
  integer(kind=intType), dimension(4, 0) :: tetraConn
  integer(kind=intType), dimension(5, 0) :: pyraConn
  integer(kind=intType), dimension(6, 0) :: prismsConn
  integer(kind=intType), parameter :: nPyra=0
  integer(kind=intType), parameter :: nTetra=0
  integer(kind=intType), parameter :: nPrisms=0
  integer(kind=intType) :: nDualNodes, nHexa, nSearchCells
  integer(kind=intType) :: planeOffset
  real(kind=realType), dimension(3, 2) :: BBox
  logical :: useBBox
  character*40 :: tmpStr

  real(kind=realType) :: timeA, timeB

  ! Explictly set level and sps to 1. This will be removed in the future.
  level = 1
  sps =1 
  print *,'doing overset interpolation'
  if (nProc > 1) then 
     print *,'Only in serial for now...'
     stop
  end if

   ! Step 1: The first step is to identify the "clusters" that are
   ! present in the grid. This is currently not implemented so we will
   ! skip this. 


   ! Step 2: The next step is to do a very rough connecitivty check
   ! using the 3d bounding box of each domain. This information is
   ! shared gloablly with all processors. This information is then used
   ! to determine which blocks must be communicated between
   ! processors. This should be a conservative estimate; if the
   ! bounding boxes do not intersect, there is no possible way any of
   ! the cells in either domain can intersect so we can completely
   ! eliminate checking those. Hoewver, it possible that the bounding
   ! boxes intersect but no cells inersect. 

  allocate(xMin(3, nDom), xMax(3, nDom))
   do nn=1,nDom
     call setPointers(nn, level, sps)

     xMin(1, nn) = minval(x(:, :, :, 1))
     xMin(2, nn) = minval(x(:, :, :, 2))
     xMin(3, nn) = minval(x(:, :, :, 3))

     xMax(1, nn) = maxval(x(:, :, :, 1))
     xMax(2, nn) = maxval(x(:, :, :, 2))
     xMax(3, nn) = maxval(x(:, :, :, 3))

  end do
  
  ! Initially assume that all meshes are overlapped by setting
  ! iOverlap to .True.. When we can determine for sure a pair do not
  ! overlap, we can set to .False. Each processor essentially only
  ! stores the local rows of iOverlap that it owns. 

  allocate(iOverlap(nDom, nDom))
  iOverlap = .True.

  ! Loop over each of my owned blocks
  do iDom=1, nDom
     
     ! Now Loop over *all* of the other blocks
     do jDom=1, nDom

        ! Now do the box overlap check
        if ( &
             xMin(1, iDom) >= xMax(1, jDom) .or. &
             xMax(1, iDom) <= xMin(1, jDom) .or. &
             xMin(2, iDom) >= xMax(2, jDom) .or. &
             xMax(2, iDom) <= xMin(2, jDom) .or. &     
             xMin(3, iDom) >= xMax(3, jDom) .or. &
             xMax(3, iDom) <= xMin(3, jDom)) then 
           
           ! These bounding boxes do not intersect. 
           iOverlap(iDom, jDom) = .False.
        end if
     end do
  end do
  
  ! Now allocate a list of the 'oBlocks'. These are the overset blocks
  allocate(oBlocks(nn))

  ! Now allocate space for all the oBlocks. This is a special data
  ! structure used for doing overset matching. It has copies from the
  ! regular blocks so it is inefficient,including the volDonor/volRecv
  ! including the modifications necessary for the boundary conditions
  
  timeA = mpi_wtime()
  do nn=1, nDom
     call setPointers(nn, level, sps)

     ! Set some sizes
     oBlocks(nn)%il = il
     oBlocks(nn)%jl = jl
     oBlocks(nn)%kl = kl

     oBlocks(nn)%ie = ie
     oBlocks(nn)%je = je
     oBlocks(nn)%ke = ke

     oBlocks(nn)%nx = nx
     oBlocks(nn)%ny = ny
     oBlocks(nn)%nz = nz

     !nDualNodes = The number of cells of the primal mesh to be used
     !to form the dual mesh

     nDualNodes   = oBlocks(nn)%ie * oBLocks(nn)%je * oBlocks(nn)%ke
     nHexa        = oBlocks(nn)%il * oBlocks(nn)%jl * oBlocks(nn)%kl
     nSearchCells = oBlocks(nn)%nx * oBLocks(nn)%ny * oBlocks(nn)%nz
     obLocks(nn)%globalID  = nn
     allocate( &
          oBlocks(nn)%xDual(3, nDualNodes), &
          oBlocks(nn)%xSearch(3, nSearchCells), &
          oBlocks(nn)%qualDonor(1, nDualNodes), &
          oBlocks(nn)%qualRecv(1, nSearchCells), &
          oBlocks(nn)%iblank(1:ie, 1:je, 1:ke), &
          oBlocks(nn)%globalCell(1:ie, 1:je, 1:ke), &
          oBlocks(nn)%hexaConn(8, nHexa), &
          oBlocks(nn)%ind(3, nDualNodes), &
          ! Maximum possible number of donors is total number of search cells. 
          oBlocks(nn)%donorIndices(8, nSearchCells), &
          oBlocks(nn)%donorIndices2(9, 3, nSearchCells), &
          oBlocks(nn)%donorFrac(3, nSearchCells), &
          oBlocks(nn)%fringeIndices(3, nSearchCells))

     oBlocks(nn)%nDonor = 0
     oBlocks(nn)%nFringe = 0

     ! Fill up the xDual for the cell duals as well as the quality for the
     ! donor and receiver. 
     
     mm = 0
     do k=1 ,ke
        do j=1, je
           do i=1, ie
              mm = mm + 1
              do iDim=1, 3
                 oBlocks(nn)%xDual(iDim, mm) = eighth*(&
                      x(i-1, j-1, k-1, iDim) + &
                      x(i  , j-1, k-1, iDim) + &
                      x(i-1, j  , k-1, iDim) + &
                      x(i  , j  , k-1, iDim) + &
                      x(i-1, j-1, k  , iDim) + &
                      x(i  , j-1, k  , iDim) + &
                      x(i-1, j  , k  , iDim) + &
                      x(i  , j  , k  , iDim))
              end do

              ! Just copy out the volumes
              oBlocks(nn)%qualDonor(1, mm) = vol(i,j,k)
        
              ! And the global cell
              oBlocks(nn)%globalCell(i, j, k) = globalCell(i, j, k)
           end do
        end do
     end do

     ! And get the search cells. *Does not include halos*
     mm =0
     do k=2 ,kl
        do j=2, jl
           do i=2, il
              mm = mm + 1
              do iDim=1, 3
                 oBlocks(nn)%xSearch(iDim, mm) = eighth*(&
                      x(i-1, j-1, k-1, iDim) + &
                      x(i  , j-1, k-1, iDim) + &
                      x(i-1, j  , k-1, iDim) + &
                      x(i  , j  , k-1, iDim) + &
                      x(i-1, j-1, k  , iDim) + &
                      x(i  , j-1, k  , iDim) + &
                      x(i-1, j  , k  , iDim) + &
                      x(i  , j  , k  , iDim))
              end do
              oBlocks(nn)%qualRecv(1, mm) = vol(i,j,k)
           end do
        end do
     end do
     
     ! Initialize al the iblanks to 1..everything is comptue until we
     ! determine otherwise. forceRecv is false until we look at the BCs.
     oBlocks(nn)%iBlank = 1

     ! We will build an ADT Tree for this block...this will prevent us
     ! having to create/destroy it many times:
     
     mm = 0
     ! These are the 'elements' of the dual mesh.
     planeOffset = ie*je
     do k=2, ke
        do j=2, je
           do i=2, ie
              mm = mm + 1
              oBlocks(nn)%hexaConn(1, mm) = (k-2)*ie*je + (j-2)*ie + (i-2) + 1
              oBlocks(nn)%hexaConn(2, mm) = oBlocks(nn)%hexaConn(1, mm) + 1 
              oBlocks(nn)%hexaConn(3, mm) = oBlocks(nn)%hexaConn(2, mm) + ie 
              oBlocks(nn)%hexaConn(4, mm) = oBlocks(nn)%hexaConn(3, mm) - 1 
              

              oBlocks(nn)%hexaConn(5, mm) = oBlocks(nn)%hexaConn(1, mm) + planeOffset
              oBlocks(nn)%hexaConn(6, mm) = oBlocks(nn)%hexaConn(2, mm) + planeOffset
              oBlocks(nn)%hexaConn(7, mm) = oBlocks(nn)%hexaConn(3, mm) + planeOffset
              oBlocks(nn)%hexaConn(8, mm) = oBlocks(nn)%hexaConn(4, mm) + planeOffset
           end do
        end do
     end do

     BBox = zero ! BBox is not used so values here not meaningful
     useBBox = .False.

     ! Make a name for the ADT
     write(tmpStr, *) nn
     tmpStr = adjustl(tmpStr)
     oBlocks(nn)%adtName = 'domain.'//tmpStr

     call adtbuildVolumeADT(nTetra, nPyra, nPrisms, nHexa, nDualNodes, &
          oBlocks(nn)%xDual, tetraConn, pyraConn, prismsConn, &
          oBlocks(nn)%hexaConn, BBox,  useBBox, MPI_COMM_SELF, &
          oBlocks(nn)%adtName)
  end do

 ! Now we will run the closest distance routine here to find the
 ! points inside. Again, don't worry about paraellel stuff yet

  call computeHolesInsideBody

 ! We now know that some of the cells are *already* blanked with
 ! iBlank = 0. For these cells. These cells cannot be a donor *or* a
 ! receiver. Adjust the quality measure on these cells accordingly. 
 do nn=1,nDom
    mm = 0
    do k=1, oBlocks(nn)%ke
       do j=1, oBlocks(nn)%je
          do i=1, oBlocks(nn)%ie
             mm = mm + 1
             if (oBlocks(nn)%iBlank(i, j, k) == 0) then 
                oBlocks(nn)%qualDonor(1, mm) = 1e30 
             end if
          end do
       end do
    end do

    ! Go through the boudary conditions and flag the overset ones Note
    ! that we are setting three slabs: The last two levels of real
    ! cells as well as the first level of halos. 
    call setPointers(nn, level, sps)
    ! do mm=1,nBocos
    !    if(BCType(mm) == OversetOuterBound) then 
    !       select case (BCFaceID(mm))
    !       case (iMin)
    !          oBlocks(nn)%qualRecv(1:3, :, :) = 1e30
    !          oBlocks(nn)%qualDonor(1:3, :, :) = -1e30
    !       case (iMax)
    !          oBlocks(nn)%qualRecv(nx:ie, :, :) = 1e30
    !          oBlocks(nn)%qualDonor(nx:ie, :, :) = -1e30
    !       case (jMin)
    !          oBlocks(nn)%qualRecv(:, 1:3, :) = 1e30
    !          oBlocks(nn)%qualDonor(:, 1:3, :) = -1e30
    !       case (jMax)
    !          oBlocks(nn)%qualRecv(:, ny:je, :) = 1e30
    !          oBlocks(nn)%qualDonor(:, ny:je, :) = -1e30
    !       case (kMin)
    !          oBlocks(nn)%qualRecv(:, :, 1:3) = 1e30
    !          oBlocks(nn)%qualDonor(:, :, 1:3) = -1e30
    !       case (kMax)
    !          oBlocks(nn)%qualRecv(:, :, nz:ke) = 1e30
    !          oBlocks(nn)%qualDonor(:, :, nz:ke) = -1e30
    !       end select
    !    end if
    ! end do
 end do

 ! Master loop over the blocks...
 do iDom=1,nDom
    do jDom=1, nDom
       if (iOverlap(iDom, jDom) .and. iDom /= jDom) then 
          call pairSearch(oBlocks(iDom), oBlocks(jDom))
       end if
    end do
 end do

 ! Copy the iblank in oBlocks to the real iblank
 do nn=1,nDom
    call setPointers(nn, level, sps)
    do k=2, kl
       do j=2, jl
          do i=2, il
             iBlank(i, j, k) = oBlocks(nn)%iBlank(i, j, k)
          end do
       end do
    end do
 end do

 call exchangeIBlanks(level, sps, commPatternCell_2nd, internalCell_2nd)

 ! Compute the overset interpolation required. 
 call initializeOversetComm
 
end subroutine computeOversetInterpolation


! subroutine test


! use overset
!   use blockPointers
! implicit none

! integer(kind=intType) :: nn, i, j, k, mm
! real(kind=realType) :: val1, val2
! ! Put a scalar field  on domain 1 and interpolate to domain 2

! !flowDoms(2, 1, 1)%w(:, :, :, 1) = one

! do nn=1,nDom

!    call setPointers(nn, 1, 1)
!    mm = 0
!    do k=1, ke
!       do j=1, je
!          do i=1, ie
!             mm = mm + 1
!             w(i, j, k, iVx) = oBlocks(nn)%xDual(1, mm) + &
!                  2*oBlocks(nn)%xDual(2, mm)
!          end do
!       end do
!    end do
! end do

! call wOverset(1, 2, 2, .False., .False., .False., .False.) 

! ! Check things on the first domain

! do nn=1,1
!    call setPointers(nn, 1,1 )
!    mm = 0
!    do k=2, kl
!       do j=2, jl
!          do i=2, il
!             mm = mm + 1
!             val1 = w(i, j, k, ivx)
!             val2 = oBlocks(nn)%xsearch(1, mM) + 2*oBlocks(nn)%xsearch(2, mm)

!             if (abs(val1- val2) > 1e-8) then 
!                print *, 'error:', i, j, k, val1, val2
!             end if

!          end do
!       end do
!    end do
! end do
! end subroutine test


     ! open(unit=1,file='vis.dat',form='formatted',status='unknown')
     ! write(1,*) "Variables = X Y Z"

     ! write(1,*) "ZONE NODES=",ie*je*ke, " Elements=", il*jl*kl, "DATAPACKING=POINT, ZONETYPE=FEBRICK"
     ! do mm=1,ie*je*ke
     !    write(1, *) oBlocks(nn)%xdual(1, mm), &
     !         oBlocks(nn)%xDual(2, mm), & 
     !         oBlocks(nn)%xDual(3, mm)
     ! end do
     ! do mm=1,il*jl*kl
     !    write(1, *) oBlocks(nn)%hexaConn(1, mm), &
     !         oBlocks(nn)%hexaConn(2, mm), &
     !         oBlocks(nn)%hexaConn(3, mm), &
     !         oBlocks(nn)%hexaConn(4, mm), &
     !         oBlocks(nn)%hexaConn(5, mm), &
     !         oBlocks(nn)%hexaConn(6, mm), &
     !         oBlocks(nn)%hexaConn(7, mm), &
     !         oBlocks(nn)%hexaConn(8, mm)
     ! end do
     ! close(1)
