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
  use inputOverset
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
  integer(kind=intType) :: nDualNodes, nPrimalNodes, nHexa, nSearchCells
  integer(kind=intType) :: planeOffset
  real(kind=realType), dimension(3, 2) :: BBox
  logical :: useBBox
  character*40 :: tmpStr

  real(kind=realType) :: timeA, timeB

  ! Explictly set level and sps to 1. This will be removed in the future.
  level = 1
  sps =1 

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

     ! Depending on the order of the overset interpolation we need
     ! slightly different data: For linear interpolation we build the
     ! ADT tree directly from the dual mesh. For the quadratic
     ! interpolation we build the tree from the primal mesh and then
     ! use that as the starting guess for a secondary quadratic search.

     nDualNodes   = oBlocks(nn)%ie * oBLocks(nn)%je * oBlocks(nn)%ke
     nPrimalNodes = oBlocks(nn)%il * oBlocks(nn)%jl * oBlocks(nn)%kl
     nSearchCells = oBlocks(nn)%nx * oBLocks(nn)%ny * oBlocks(nn)%nz

     if (oversetInterpolation == linear) then 
        nHexa        = oBlocks(nn)%il * oBlocks(nn)%jl * oBlocks(nn)%kl
     else
        nHexa        = oBlocks(nn)%nx * oBLocks(nn)%ny * oBlocks(nn)%nz
     end if

     obLocks(nn)%globalID  = nn
     allocate( &
          oBlocks(nn)%xDual(3, nDualNodes), &
          oBlocks(nn)%xPrimal(3, nPrimalNodes), &
          oBlocks(nn)%xSearch(3, nSearchCells), &
          oBlocks(nn)%qualRecv(2:il, 2:jl, 2:kl), &
          oBlocks(nn)%iblank(1:ie, 1:je, 1:ke), &
          oBlocks(nn)%globalCell(1:ie, 1:je, 1:ke), &
          oBlocks(nn)%hexaConn(8, nHexa), &
          oBlocks(nn)%qualDonor(1:ie, 1:je, 1:ke), &
          ! Maximum possible number of donors is total number of search cells. 
          oBlocks(nn)%donorFrac(3, nSearchCells), &
          oBlocks(nn)%fringeIndices(3, nSearchCells))

     if (oversetInterpolation == linear) then 
        allocate(oBlocks(nn)%donorIndices(8, nSearchCells))
     else
        allocate(oBlocks(nn)%donorIndices(27, nSearchCells))
     end if

     oBlocks(nn)%nDonor = 0
     oBlocks(nn)%nFringe = 0

     ! Fill up the xDual for the dual cells
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

              ! And the global cell
              oBlocks(nn)%globalCell(i, j, k) = globalCell(i, j, k)
           end do
        end do
     end do

     ! Simply copy over the primal nodes
     mm = 0
     do k=1 ,kl
        do j=1, jl
           do i=1, il
              mm = mm + 1
              do iDim=1, 3
                 oBlocks(nn)%xPrimal(iDim, mm) = x(i, j, k, iDim)
              end do
           end do
        end do
     end do

     do k=1, ke
        do j=1, je
           do i=1, ie
              ! Just copy out the volumes
              oBlocks(nn)%qualDonor(i, j, k) = vol(i,j,k)
           end do
        end do
     end do

     ! And get the search cells. *Does not include halos*. This is the
     ! same in boht cases. 
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
              oBlocks(nn)%qualRecv(i, j, k) = vol(i,j,k)
           end do
        end do
     end do


     ! Initialize al the iblanks to 1..everything is comptue until we
     ! determine otherwise. forceRecv is false until we look at the BCs.
     oBlocks(nn)%iBlank = 1

     ! We will build an ADT Tree for this block...this will prevent us
     ! having to create/destroy it many times:

     if (oversetInterpolation == linear) then 
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
     else
        mm = 0
        ! These are the 'elements' of the primal mesh
        planeOffset = il*jl
        do k=2, kl
           do j=2, jl
              do i=2, il
                 mm = mm + 1
                 oBlocks(nn)%hexaConn(1, mm) = (k-2)*il*jl + (j-2)*il + (i-2) + 1
                 oBlocks(nn)%hexaConn(2, mm) = oBlocks(nn)%hexaConn(1, mm) + 1 
                 oBlocks(nn)%hexaConn(3, mm) = oBlocks(nn)%hexaConn(2, mm) + il
                 oBlocks(nn)%hexaConn(4, mm) = oBlocks(nn)%hexaConn(3, mm) - 1 

                 oBlocks(nn)%hexaConn(5, mm) = oBlocks(nn)%hexaConn(1, mm) + planeOffset
                 oBlocks(nn)%hexaConn(6, mm) = oBlocks(nn)%hexaConn(2, mm) + planeOffset
                 oBlocks(nn)%hexaConn(7, mm) = oBlocks(nn)%hexaConn(3, mm) + planeOffset
                 oBlocks(nn)%hexaConn(8, mm) = oBlocks(nn)%hexaConn(4, mm) + planeOffset
              end do
           end do
        end do
     end if

     BBox = zero ! BBox is not used so values here not meaningful
     useBBox = .False.

     ! Make a name for the ADT
     write(tmpStr, *) nn
     tmpStr = adjustl(tmpStr)
     oBlocks(nn)%adtName = 'domain.'//tmpStr
     if (oversetInterpolation == linear) then 
        call adtbuildVolumeADT(nTetra, nPyra, nPrisms, nHexa, nDualNodes, &
             oBlocks(nn)%xDual, tetraConn, pyraConn, prismsConn, &
             oBlocks(nn)%hexaConn, BBox,  useBBox, MPI_COMM_SELF, &
             oBlocks(nn)%adtName)
     else
        call adtbuildVolumeADT(nTetra, nPyra, nPrisms, nHexa, nPrimalNodes, &
             oBlocks(nn)%xPrimal, tetraConn, pyraConn, prismsConn, &
             oBlocks(nn)%hexaConn, BBox,  useBBox, MPI_COMM_SELF, &
             oBlocks(nn)%adtName)
     end if
  end do


  ! Now we will run the closest distance routine here to find the
  ! points inside. Again, don't worry about paraellel stuff yet

  call computeHolesInsideBody

  ! We now know that some of the cells are *already* blanked with
  ! iBlank = 0. For these cells. These cells cannot be a donor *or* a
  ! receiver. Adjust the quality measure on these cells accordingly. 
 

     do nn=1,nDom
        do k=1, oBlocks(nn)%ke
           do j=1, oBlocks(nn)%je
              do i=1, oBlocks(nn)%ie
                 if (oBlocks(nn)%iBlank(i, j, k) == 0) then 
                    oBlocks(nn)%qualDonor(i, j, k) = 1e30 
                 end if
              end do
           end do
        end do


        ! **************************************************************
        ! *WARNING*: i Think this only works for having at least 2
        ! cells in the direction off of a outer boundary. That should
        ! always be the case
        ! **************************************************************

        ! Go through the boudary conditions and flag the overset ones Note
        ! that we are setting three slabs: The last two levels of real
        ! cells as well as the first level of halos. 
        call setPointers(nn, level, sps)
        do mm=1,nBocos
           if(BCType(mm) == OversetOuterBound) then 
              select case (BCFaceID(mm))
              case (iMin)
                 oBlocks(nn)%qualRecv(2:3, :, :) = 1e30
                 oBlocks(nn)%qualDonor(1:3, :, :) = 1e30
              case (iMax)
                 oBlocks(nn)%qualRecv(nx:il, :, :) = 1e30
                 oBlocks(nn)%qualDonor(nx:ie, :, :) = 1e30
              case (jMin)
                 oBlocks(nn)%qualRecv(:, 2:3, :) = 1e30
                 oBlocks(nn)%qualDonor(:, 1:3, :) = 1e30
              case (jMax)
                 oBlocks(nn)%qualRecv(:, ny:jl, :) = 1e30
                 oBlocks(nn)%qualDonor(:, ny:je, :) = 1e30
              case (kMin)
                 oBlocks(nn)%qualRecv(:, :, 2:3) = 1e30
                 oBlocks(nn)%qualDonor(:, :, 1:3) = 1e30
              case (kMax)
                 oBlocks(nn)%qualRecv(:, :, nz:kl) = 1e30
                 oBlocks(nn)%qualDonor(:, :, nz:ke) = 1e30
              end select
           end if
        end do
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
