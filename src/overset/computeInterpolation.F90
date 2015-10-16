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
  use cgnsGrid
  use stencils
  implicit none

  ! Local Variables
  integer(kind=intType) ::i, j, k,  level, sps, ii, iDom, jDom, iDim, nn, mm, tmp
  integer(kind=intType) ::  jj, kk, ip2, im2, jp2, jm2, kp2, km2, ind(3)
  integer(kind=intType) :: donorBlockID, i_stencil
  real(kind=realType), dimension(3, nDom) :: xMinLocal
  real(kind=realType), dimension(:, :), allocatable :: xMin, xMax
  logical, dimension(:, :), allocatable :: iOverlap
  integer(kind=intType), dimension(3, nDom) :: localDim
  integer(kind=intTYpe) :: iProc, ierr, blockID
  ! Buffers variables
  integer(kind=intType) :: buffSize, nRealSend, nIntSend, tag
  integer(kind=intType) :: nRealRecv, nIntRecv
  real(kind=realType), dimension(:), allocatable :: realBuffer
  integer(kind=intType), dimension(:), allocatable :: intBuffer
  integer(kind=intTYpe), dimension(:, :, :), allocatable :: iblankTmp
  ! Data for ADT building
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
  integer status(MPI_STATUS_SIZE) 
  real(kind=realType) :: timeA, timeB, timeC, timeD, time1
  real(kind=realType) :: procTimes(0:nProc-1)
  logical :: computeCellFound

  ! Explictly set level and sps to 1. This will be removed in the future.
  level = 1
  sps = 1 

  call initialize_stencils()
  ! Step 1: The first step is to identify the "clusters" that are
  ! present in the grid. We will use the original CGNS information
  ! since this is slightly more readily available. Futhermore, there
  ! will generally be less blocks in the CGNS file that compute blocks
  ! due to splitting. Once we have labelled the CGNS blocks, we can
  ! simply label the compute blocks by looking at what CGNS block they
  ! came from.  Since all procs have this information we can simply
  ! run on all procs
  timeA = mpi_wtime()

  call determineClusters

  ! We will first perform the search for the cells that are inside the
  ! body.  This is done in parallel in the same manner as the wall
  ! distance calc. 
  call computeHolesInsideBody
  call exchangeIBlanks(level, sps, commPatternCell_2nd, internalCell_2nd)
  timeB = mpi_wtime()

  ! Step two: We need to gather all the block information to the root
  ! processor who will do all the hard work. The minimum amount of
  ! information we need to transfer is the primal mesh (including the
  ! first level of halos), the globalCell, iBlank and a cellStatus
  ! array. Everything else can be recomputed from this information. 

  ! Gather the dimensions of all blocks to everyone
  call mpi_allreduce(nDom, nDomTotal, 1, sumb_integer, MPI_SUM, &
       sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Store the sizes of the local blocks
  do nn=1,nDom
     call setPointers(nn, level, sps)

     ! Store the 'l' sizes for transferring
     localDim(1, nn) = il
     localDim(2, nn) = jl
     localDim(3, nn) = kl
  end do

  ! Allocate the space we need for the numbers and cumulative form
  allocate(nDomProc(0:nProc-1), cumDomProc(0:nProc), dims(3, nDomTotal))

  ! Receive the number of domains from each proc using an allgather.
  call mpi_allgather(nDom, 1, sumb_integer, nDomProc, 1, sumb_integer, &
       sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Compute the cumulative format:
  cumDomProc(0) = 0
  do iProc=1, nProc
     cumDomProc(iProc) = cumDomProc(iProc-1) + nDomProc(iProc-1)
  end do

  ! We will also allgather all of the block sizes which will make
  ! things a little easier since everyone will know the proper sizes
  ! for the sends
  call mpi_allgatherV(localDim, nDom*3, sumb_integer, dims, 3*nDomProc, &
       3*cumDomProc, sumb_integer, sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Determine the maximum size of the buffer we will need for any one
  ! block. It will be the 0:ie etc. So it is 2 more than the 'l' sizes
  ! currently in dims.

  buffSize = 0
  do nn=1, nDomTotal
     buffSize = max(buffSize, (dims(1, nn) + 2)*(dims(2, nn) + 2)*(dims(3, nn) + 2))
  end do

  ! While not strictly necessary, all procs are using the max sized
  ! buffer. Each local proc *could* use its own (potentially) smaller
  ! buffer

  allocate(realBuffer(4*buffSize), intBuffer(3*buffSize+2))

  ! Now we do the actual communication of the blocks to the root:

  if (myid > 0) then 

     do nn=1,nDom
        call setPointers(nn, level, sps)
        call packBlock

        ! tag is the global block index
        tag = cumDomProc(myID) + nn 

        call mpi_send(realBuffer, nRealSend, sumb_real, 0, tag, &
             sumb_comm_world, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        call mpi_send(intBuffer, nIntSend, mpi_integer, 0, tag, &
             sumb_comm_world, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

     end do

  else

     ! Allocate space for the 'oBlocks' the data structure used
     ! storing all the blocks
     allocate(oBlocks(nDomTotal))

     ! Loop over all processors except myself
     do iProc=1, nProc-1

        do blockID=cumDomProc(iProc)+1, cumDomProc(iProc+1)
           tag = blockID
           call MPI_Recv(realBuffer, 4*buffSize, sumb_real, iProc, tag, &
                sumb_comm_world, status, ierr)

           call MPI_Get_count(status, sumb_real, nRealRecv, ierr)

           call MPI_Recv(intBuffer, 3*buffSize+2, sumb_integer, iProc, tag, &
                sumb_comm_world, status, ierr)
           call MPI_Get_count(status, sumb_real, nIntRecv, ierr)

           call unpackBlock
        end do
     end do

     ! And we have to do our own blocks:
     do nn=1,nDom
        call setPointers(nn, level, sps)
        call packBlock
        blockID = nn
        call unpackBlock
     end do

     ! Now we have all the oBlocks setup here on the root
     ! processor. We still have a few additional things to do with the
     ! oBlocks

     do nn=1, nDomTotal

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

        allocate( &
             oBlocks(nn)%xDual(3, nDualNodes), &
             oBlocks(nn)%xPrimal(3, nPrimalNodes), &
             oBlocks(nn)%xSearch(3, 2:oBlocks(nn)%il, 2:oBlocks(nn)%jl, 2:oBlocks(nn)%kl), &
             oBlocks(nn)%hexaConn(8, nHexa), &

             ! Maximum possible number of donors is total number of search cells. 
             oBlocks(nn)%donorFrac(3, 2*nSearchCells), &
             oBlocks(nn)%fringeIndices(3, 2*nSearchCells))

        if (oversetInterpolation == linear) then 
           allocate(oBlocks(nn)%donorIndices(8, 2*nSearchCells))
        else
           allocate(oBlocks(nn)%donorIndices(27, 2*nSearchCells))
        end if

        oBlocks(nn)%nDonor = 0
        oBlocks(nn)%nFringe = 0
        oBlocks(nn)%globalBlockID = nn
        ! Fill up the xDual for the dual cells
        mm = 0
        do k=1, oBlocks(nn)%ke
           do j=1, oBlocks(nn)%je
              do i=1, oBlocks(nn)%ie
                 mm = mm + 1
                 oBlocks(nn)%xDual(:, mm) = eighth*(&
                      oBlocks(nn)%x(:, i-1, j-1, k-1) + &
                      oBlocks(nn)%x(:, i  , j-1, k-1) + &
                      oBlocks(nn)%x(:, i-1, j  , k-1) + &
                      oBlocks(nn)%x(:, i  , j  , k-1) + &
                      oBlocks(nn)%x(:, i-1, j-1, k  ) + &
                      oBlocks(nn)%x(:, i  , j-1, k  ) + &
                      oBlocks(nn)%x(:, i-1, j  , k  ) + &
                      oBlocks(nn)%x(:, i  , j  , k  ))
              end do
           end do
        end do

        ! Simply copy over the primal nodes. In this case, just the owned ones. 
        mm = 0
        do k=1, oBlocks(nn)%kl
           do j=1, oBlocks(nn)%jl
              do i=1, oBlocks(nn)%il
                 mm = mm + 1
                 oBlocks(nn)%xPrimal(:, mm) = oBlocks(nn)%x(:, i, j, k)
              end do
           end do
        end do

        ! And get the search cells. *Does not include halos*.
        do k=2, oBlocks(nn)%kl
           do j=2, oBlocks(nn)%jl
              do i=2, oBlocks(nn)%il
                 oBlocks(nn)%xSearch(:, i, j, k) = eighth*(&
                      oBlocks(nn)%x(:, i-1, j-1, k-1) + &
                      oBlocks(nn)%x(:, i  , j-1, k-1) + &
                      oBlocks(nn)%x(:, i-1, j  , k-1) + &
                      oBlocks(nn)%x(:, i  , j  , k-1) + &
                      oBlocks(nn)%x(:, i-1, j-1, k  ) + &
                      oBlocks(nn)%x(:, i  , j-1, k  ) + &
                      oBlocks(nn)%x(:, i-1, j  , k  ) + &
                      oBlocks(nn)%x(:, i  , j  , k  ))
              end do
           end do
        end do

        ! We will build an ADT Tree for this block...this will prevent us
        ! having to create/destroy it many times:

        if (oversetInterpolation == linear) then 
           mm = 0
           ! These are the 'elements' of the dual mesh.
           planeOffset = oBlocks(nn)%ie*oBlocks(nn)%je
           do k=2, oBlocks(nn)%ke
              do j=2, oBlocks(nn)%je
                 do i=2, oBlocks(nn)%ie
                    mm = mm + 1
                    oBlocks(nn)%hexaConn(1, mm) = (k-2)*planeOffset + (j-2)*oBlocks(nn)%ie + (i-2) + 1
                    oBlocks(nn)%hexaConn(2, mm) = oBlocks(nn)%hexaConn(1, mm) + 1 
                    oBlocks(nn)%hexaConn(3, mm) = oBlocks(nn)%hexaConn(2, mm) + oBlocks(nn)%ie
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
           planeOffset = oBlocks(nn)%il*oBlocks(nn)%jl
           do k=2, oBlocks(nn)%kl
              do j=2, oBlocks(nn)%jl
                 do i=2, oBlocks(nn)%il
                    mm = mm + 1
                    oBlocks(nn)%hexaConn(1, mm) = (k-2)*planeOffset + (j-2)*oBlocks(nn)%il + (i-2) + 1
                    oBlocks(nn)%hexaConn(2, mm) = oBlocks(nn)%hexaConn(1, mm) + 1 
                    oBlocks(nn)%hexaConn(3, mm) = oBlocks(nn)%hexaConn(2, mm) + oBlocks(nn)%il
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

        ! Before we do the searching we need to modify the the
        ! receiver/donor quality depending on the information we have
        ! in the cellStatus array. We also have to initialize the
        ! oBlocks's iBlank data.
        
        oBlocks(nn)%iBlank = 1

        do k=1, oBlocks(nn)%ke
           do j=1, oBlocks(nn)%je
              do i=1, oBlocks(nn)%ie

                 select case (oBlocks(nn)%cellStatus(i, j, k))

                 case (0) 
                    ! Blanekd Cell: Set the iblank value and ensure it
                    ! cannot be a donor
                    oBlocks(nn)%iBlank(i, j, k) = 0
                    oBlocks(nn)%invalidDonor(i, j, k) = .True. 

                    ! We have to be careful with holes: We need to
                    ! make sure that there is sufficient frige padding
                    ! around them. Therefore we loop over the stencil
                    ! for this cell and for all the cells in its
                    ! stencil that are not already holes, flag them as
                    ! forceRecv. Be careful of bounds

                    stencilLoop: do i_stencil=1,  N_visc_drdw
                       ii = visc_drdw_stencil(i_stencil, 1) + i
                       jj = visc_drdw_stencil(i_stencil, 2) + j
                       kk = visc_drdw_stencil(i_stencil, 3) + k
                          
                       ii = max(min(oBlocks(nn)%ie, ii), 1)
                       jj = max(min(oBlocks(nn)%je, jj), 1)
                       kk = max(min(oBlocks(nn)%ke, kk), 1)

                       if (oBlocks(nn)%iblank(ii, jj, kk) /= 0) then 
                          oblocks(nn)%forceRecv(ii, jj, kk) = .True. 
                       end if
                    end do stencilLoop
                 case (-1)
                    ! Cell as ben flagged as necessairly being a
                    ! receiver. We therefore force it to be a receiver
                    ! and prevent it from being a donor
                    oBlocks(nn)%invalidDonor(i, j, k) = .True. 
                    oBlocks(nn)%forceRecv(i, j, k) = .True. 
                    oBlocks(nn)%iBlank(i, j, k) = -1
                 end select

              end do
           end do
        end do
     end do ! All domain looop

     ! Step 2: The next step is to do a very rough connecitivty check
     ! using the 3d bounding box of each domain. This should be a
     ! conservative estimate; if the bounding boxes do not intersect,
     ! there is no possible way any of the cells in either domain can
     ! intersect so we can completely eliminate checking
     ! those. Hoewver, it possible that the bounding boxes intersect
     ! but no cells inersect.

     allocate(xMin(3, nDomTotal), xMax(3, nDomTotal))
     do nn=1,nDomTotal

        xMin(1, nn) = minval(oBlocks(nn)%x(1, :, :, :))
        xMin(2, nn) = minval(oBlocks(nn)%x(2, :, :, :))
        xMin(3, nn) = minval(oBlocks(nn)%x(3, :, :, :))

        xMax(1, nn) = maxval(oBlocks(nn)%x(1, :, :, :))
        xMax(2, nn) = maxval(oBlocks(nn)%x(2, :, :, :))
        xMax(3, nn) = maxval(oBlocks(nn)%x(3, :, :, :))

     end do

     ! Initially assume that all meshes are overlapped by setting
     ! iOverlap to .True.. When we can determine for sure a pair do
     ! not overlap, we can set to .False. Each processor essentially
     ! only stores the local rows of iOverlap that it owns.  Note that
     ! the cluster anlysis is guaranteed to set the diagonal to False
     ! so we don't have to treat it specially later on.

     allocate(iOverlap(nDomTotal, nDomTotal))
     iOverlap = .True.

     ! Loop over all blocks
     do iDom=1, nDomTotal

        ! Now Loop over *all* of the other blocks
        do jDom=1, nDomTotal 

           ! We can eliminate some of pairs using the cluser analysis:
           if (oBlocks(iDom)%cluster == oBlocks(jDom)%cluster) then 
              iOverlap(iDom, jDom) = .False.
           end if

           ! Only do the spatial check if we haven't elminated the
           ! connection through the cluster check
           if (iOverlap(iDom, jDom)) then 

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
           end if
        end do
     end do

     ! ii = 0
     ! do i=1,nDomTotal
     !    print *,ioverlap(i, :)
     !    do j=1,nDomTotal
     !       if (iOverlap(i,j)) then 
     !          ii = ii + 1
     !       end if
     !    end do
     ! end do
     ! print *,ii, ' nonzeros in connectivity matrix of size:', shape(ioverlap)

     ! Master loop over the blocks for doing the implicit hole cutting
     ! technique
     
     iDom = 0
     do iProc=0, nProc-1
        do i=1,nDomProc(iProc)
           iDom = iDom + 1
           oBlocks(iDom)%procID = iProc
        end do
     end do

     do iDom=1,nDomTotal
        timec = mpi_wtime()

        do jDom=1, nDomTotal
           if (iOverlap(iDom, jDom)) then 
              call pairSearch(oBlocks(iDom), oBlocks(jDom))
           end if
        end do
     end do

     ! Flag cells that are actually donors...this will be much more
     ! complex in parallel...here we can just read the index and set
     ! directly.
     do nn=1, nDomTotal
        do k=2, oBlocks(nn)%kl
           do j=2, oBlocks(nn)%jl
              do i=2, oBlocks(nn)%il
                 if (oBLocks(nn)%donors(i, j, k)%donorProcID >= 0) then 

                    ! copy for easier code reading
                    donorBlockID = oBLocks(nn)%donors(i, j, k)%donorBlockID
                    ind = oBlocks(nn)%donors(i, j, k)%ind

                    do kk=0,1
                       do jj=0,1
                          do ii=0,1

                             oBlocks(donorBlockID)%donorStatus(&
                                  ind(1)+ii, ind(2)+jj, ind(3)+kk) = .True.
                          end do
                       end do
                    end do
                 end if
              end do
           end do
        end do
     end do
  end if ! myid == 0

  ! Update the iBlank stuff
  call tmpUpdateIblanks

  if (myid == 0) then 
     ! This is the "onera" correction...if a cell is both a receiver
     ! and a donor, force it to be a compute cell and remove it's
     ! receiver status. Do this as long as it isn't a forced
     ! receiver...we can't touch the forced receivers. This operation
      ! can be done in parallel. 
     do nn=1, nDomTotal
        do k=2, oBlocks(nn)%kl
           do j=2, oBlocks(nn)%jl
              do i=2, oBlocks(nn)%il
                 if (oBlocks(nn)%iblank(i,j,k) == -1) then
                 if ( oBLocks(nn)%recvStatus(i, j, k) .and. &
                      oBlocks(nn)%donorStatus(i, j, k) .and. &
                      .not. oBlocks(nn)%forceRecv(i, j, k)  ) then

                    ! Force it to be compute
                    oBlocks(nn)%iBlank(i,j,k) = 1
                    oBlocks(nn)%recvStatus(i, j, k) = .False. 

                 end if
              end if
              end do
           end do
        end do
     end do
  end if

  ! Exchange iblanks again
  call tmpUpdateIblanks

  if (myid == 0) then 
  ! iBlank = 0 flooding to reduce fringes
     do nn=1, nDomTotal

        ! Allocate temp variable and initialize to actual values. This
        ! is necessary since we cannot update the iblank array as we
        ! go, since that will affect the result after it. 

         allocate(iBlankTmp(0:oBlocks(nn)%ie+1, 0:oBlocks(nn)%je+1, 0:oBlocks(nn)%ke+1))
        iBlankTmp = oBlocks(nn)%iBlank

        do k=2, oBlocks(nn)%kl
           do j=2, oBlocks(nn)%jl
              do i=2, oBlocks(nn)% il

                 ! See if we can make this fringe a actual hole
                 if (oBlocks(nn)%iBlank(i,j,k) == -1) then 

                    computeCellFound = .False.

                    stencilLoop2: do i_stencil=1, N_visc_drdw
                       ii = visc_drdw_stencil(i_stencil, 1) + i
                       jj = visc_drdw_stencil(i_stencil, 2) + j
                       kk = visc_drdw_stencil(i_stencil, 3) + k

                       if (oBlocks(nn)%iBlank(ii, jj, kk) == 1) then 
                          computeCellFound = .True.
                       end if
                    end do stencilLoop2

                    if (.not. computeCellFound) then 
                       ! Everything was interpoolated so hard blank to zero
                       iBlankTmp(i, j, k) = 0
                    end if
                 end if
              end do
           end do
        end do

        ! Finally copy back and deallocate
         oBlocks(nn)%iBlank = iBlankTmp
        deallocate(iBlanktmp)
     end do
  end if

  ! Exchange iblanks again --this time, what is left is the actual
  ! iblank array is what we want. 

  call tmpUpdateIblanks

  if (myid == 0) then 
     ! A final overset check.
     call checkOverset
  end if


  ! if (myid == 0) then 
  !    timeB = mpi_wtime()
  !    print *, 'total time:', timeB-timeA
    
  !    open(unit=1,file='vis.dat',form='formatted',status='unknown')
  !    write(1, *) "Variables = X Y Z iBlank"
  !    print *,'writing...'
  !    do nn=1,nDomTotal

  !       write(1,*) "ZONE I=", oBlocks(nn)%il, " J=",oBlocks(nn)%jl , "K=", oBlocks(nn)%kl
  !       write(1, *) "DATAPACKING=BLOCK"
  !       write(1, *) "VARLOCATION=([1,2,3]=NODAL, [4]=CELLCENTERED)"

  !       do iDim=1,3
  !          do k=1,oBlocks(nn)%kl
  !             do j=1,oBlocks(nn)%jl
  !                do i=1,oBlocks(nn)%il
  !                   write(1, *) oBlocks(nn)%x(idim,i,j,k)
  !                end do
  !             end do
  !          end do
  !       end do
        
  !       do k=2,oBlocks(nn)%kl
  !          do j=2,oBlocks(nn)%jl
  !             do i=2,oBlocks(nn)%il
  !                write(1, *) oBlocks(nn)%iblank(i,j,k)
  !             end do
  !          end do
  !       end do
  !    end do
  !    close(1)
  ! end if

  ! Compute the overset interpolation required. 
  call initializeOversetComm


contains

  subroutine packBlock
    implicit none
    ! Pack the required data into the realBuffer and intBuffer using
    ! the blockPointers which are assumed already set. Note that we
    ! take this opportunity to reorder 'x' to have the x-y-z values
    ! next to each other in memory.

    integer(kind=intType), allocatable, dimension(:, :, :) :: cellStatus
    logical dump
    allocate(cellStatus(1:ie, 1:je, 1:ke))
    cellStatus = 1
    nRealSend = 0

    do k=0,ke
       do j=0,je
          do i=0,ie
             do iDim=1,3
                nRealSend = nRealSend + 1
                realBuffer(nRealSend) = x(i, j, k, iDim)
             end do
          end do
       end do
    end do
    dump = .false.
    do k=1,ke
       do j=1,je
          do i=1,ie
             nRealSend = nRealSend + 1
             realBuffer(nRealSend) = vol(i, j, k)
          end do
       end do
    end do

    nIntSend = 0

    ! We need to flag the status of cells according to the
    ! oversetOuterBoundary condition. cellStatus is by default set
    ! to 1

    do mm=1,nBocos
       if(BCType(mm) == OversetOuterBound) then 
          select case (BCFaceID(mm))
          case (iMin)
             cellStatus(1:3, :, :) = -1
          case (iMax)
             cellStatus(nx:ie, :, :) = -1
          case (jMin)
             cellStatus(:, 1:3, :) = -1
          case (jMax)
             cellStatus(:, ny:je, :) = -1
          case (kMin)
             cellStatus(:, :, 1:3) = -1
          case (kMax)
             cellStatus(:, :, nz:ke) = -1
          end select
       end if
    end do

    ! Compbine the cellStatus information array with our current
    ! iBlank data we got from the holesInsideBody calc.
    do k=1,ke
       do j=1,je
          do i=1,ie
             if (iblank(i, j, k) == 0) then 
                cellStatus(i, j, k) = iBlank(i, j, k)
             end if
          end do
       end do
    end do

    ! Now pack the cellStatus and globalCell
    do k=1, ke
       do j=1, je
          do i=1, ie
             nIntSend = nIntSend + 1
             intBuffer(nIntSend) = cellStatus(i, j, k)

             nIntSend = nIntSend + 1
             intBuffer(nIntSend) = globalCell(i, j, k)
          end do
       end do
    end do

    ! Finally add the cluster ID and the local index of the block
    nIntSend = nIntSend + 1
    intBuffer(nIntSend) = cgnsDoms(nbkglobal)%cluster

    nIntSend = nIntSend + 1
    intBuffer(nIntSend) = nn
 
    deallocate(cellStatus)

  end subroutine packBlock


  subroutine unpackBlock
    implicit none
    integer(kind=intType) :: ii, ie, je, ke, il, jl, kl

    ! Unpack the data in the real and iteger buffers to the oBlocks 
    ! data structure
    il = dims(1, blockID)
    jl = dims(2, blockID)
    kl = dims(3, blockID)

    ie = il + 1
    je = jl + 1
    ke = kl + 1

    oBlocks(blockID)%il = il
    oBlocks(blockID)%jl = jl
    oBlocks(blockID)%kl = kl

    oBlocks(blockID)%ie = ie
    oBlocks(blockID)%je = je
    oBlocks(blockID)%ke = ke

    oBlocks(blockID)%ib = ie + 1
    oBlocks(blockID)%jb = je + 1
    oBlocks(blockID)%kb = ke + 1

    oBlocks(blockID)%nx = il - 1
    oBlocks(blockID)%ny = jl - 1 
    oBlocks(blockID)%nz = kl - 1


    nPrimalNodes = (ie+1) * (je +1)* (ke+1)

    allocate( &
         oBlocks(blockID)%x(3, 0:ie, 0:je, 0:ke), &
         oBlocks(blockID)%iblank(0:ie+1, 0:je+1, 0:ke+1), &
         oBlocks(blockID)%globalCell(1:ie, 1:je, 1:ke), &
         oBlocks(blockID)%qualRecv(1:ie, 1:je, 1:ke), &
         oBlocks(blockID)%qualDonor(1:ie, 1:je, 1:ke), &
         oBlocks(blockID)%cellStatus(1:ie, 1:je, 1:ke), &
         oBlocks(blockID)%recvStatus(1:ie, 1:je, 1:ke), &
         oBlocks(blockID)%donorStatus(1:ie, 1:je, 1:ke), &
         oBlocks(blockID)%forceRecv(1:ie, 1:je, 1:ke), &
         oBlocks(blockID)%invalidDonor(1:ie, 1:je, 1:ke), &
         oBlocks(blockID)%donors(2:il, 2:jl, 2:kl))

    oBlocks(blockID)%recvStatus = .False.
    oBlocks(blockID)%donorStatus = .False.
    oBlocks(blockID)%forceRecv = .False. 
    oBlocks(blockID)%invalidDonor = .False. 

    ! Initialize the data in the donor derived type
    do k=2,kl
       do j=2,jl
          do i=2,il
             oBlocks(blockID)%donors(i, j, k)%donorProcID = -1
             oBlocks(blockID)%donors(i, j, k)%donorBlockID = -1
             oBlocks(blockID)%donors(i, j, k)%frac = zero
          end do
       end do
    end do

    ! Copy out the primal nodes
    ii = 0
    do k=0, ke
       do j=0, je
          do i=0, ie
             do iDim=1,3
                ii = ii + 1
                oBlocks(blockID)%x(iDim, i, j, k) = realBuffer(ii)
             end do
          end do
       end do
    end do

    ! Copy out the volume. Also compute the min volume for this block
    ! as we go through the loop since this is basically free.
    oBlocks(blockID)%minVolume = large
    do k=1, ke
       do j=1, je
          do i=1, ie
             ii = ii + 1
             oBlocks(blockID)%qualDonor(i,j,k) = realBuffer(ii)
             oBlocks(blockID)%qualRecv(i,j,k) = realBuffer(ii)
             oBlocks(blockID)%minVolume = min(oBlocks(blockID)%minVolume, realBuffer(ii))
          end do
       end do
    end do

    ! Copy out the cellStatus and globalCell
    ii = 0
    do k=1, ke
       do j=1, je
          do i=1, ie

             ii = ii + 1
             oBlocks(blockID)%cellStatus(i, j, k) = intBuffer(ii)

             ii = ii + 1
             oBlocks(blockID)%globalCell(i, j, k) = intBuffer(ii)
          end do
       end do
    end do

    ! Finally pull out the clusterID and local index
    ii = ii + 1
    oBlocks(blockID)%cluster = intBuffer(ii)

    ii = ii + 1
    oBlocks(blockID)%localBlockID = intBuffer(ii)

  end subroutine unpackBlock

end subroutine computeOversetInterpolation

subroutine tmpUpdateIblanks

  ! Temporary routine that takes all the iblank values on the
  ! 'oBlocks' on the root processor, sends them back to the proper
  ! process, does an iBlank exchange and then sends them all back to
  ! the root processor again. This is *horrendously* inefficient and
  ! will be replaced when the fully parallel connectivity algorithm is
  ! complete. 
  use communication
  use overset
  use blockPointers
  implicit none
  integer(kind=intType) :: nn, i, j, k, ii, tag, dest, iDom, ierr
  integer(kind=intType), dimension(:), allocatable :: intBuffer
  integer status(MPI_STATUS_SIZE) 
  
  if (myid == 0) then 
     
     nDomTotal = size(oBlocks)

     do iDom=1, nDomTotal
        if (oBlocks(iDom)%procID /= 0) then 

           ! We need to pack and send
           allocate(intBuffer((1+oBlocks(iDom)%ib)*(1+oBlocks(iDom)%jb)*(1+oBlocks(iDom)%kb)))
           ii = 0
           do k=0, oBlocks(iDom)%kb 
              do j=0, oBlocks(iDom)%jb
                 do i=0, oBlocks(iDom)%ib
                    ii = ii + 1
                    intBuffer(ii) = oBlocks(iDom)%iBlank(i, j, k)
                 end do
              end do
           end do
           
           ! Send 
           tag = iDom
           dest = oBlocks(iDom)%procID
           call mpi_send(intBuffer, size(intBuffer), sumb_integer, dest, tag, &
                sumb_comm_world, ierr)
           call ECHK(ierr, __FILE__, __LINE__)

           deallocate(intBuffer)

        else
           ! Do the local copy:
           flowDoms(iDom, 1, 1)%iBlank = oBlocks(iDom)%iBlank
        end if
     end do
     
  else
     
     ! Recieve the iblank
     do nn=1, nDom
         call setPointers(nn, 1, 1)
         allocate(intBuffer((ib+1)*(jb+1)*(kb+1)))

         tag = cumDomProc(myid) + nn

         call MPI_Recv(intBuffer, size(intBuffer), sumb_integer, 0, tag, &
              sumb_comm_world, status, ierr)

         ii = 0
         do k=0, kb
            do j=0, jb
               do i=0, ib
                  ii = ii + 1
                  iBlank(i, j, k) = intBuffer(ii)
               end do
            end do
         end do
         deallocate(intBuffer)
      end do
   end if

   ! Now call the actual iblank communication
   call exchangeIBlanks(1, 1, commPatternCell_2nd, internalCell_2nd)

   ! And finally reverse the processes...gather everything back to the oBLocks

   if (myid == 0) then 
      do iDom=1, nDomTotal
         if (oBlocks(iDom)%procID /= 0) then 
            
           ! We need to pack and send
           allocate(intBuffer((1+oBlocks(iDom)%ib)*(1+oBlocks(iDom)%jb)*(1+oBlocks(iDom)%kb)))
           tag = iDom
           dest = oBlocks(iDom)%procID
           call mpi_recv(intBuffer, size(intBuffer), sumb_integer, dest, tag, &
                sumb_comm_world, status, ierr)
           call ECHK(ierr, __FILE__, __LINE__)

           ii = 0
           do k=0, oBlocks(iDom)%kb 
              do j=0, oBlocks(iDom)%jb
                 do i=0, oBlocks(iDom)%ib
                    ii = ii + 1
                    oBlocks(iDom)%iBlank(i, j, k) = intBuffer(ii)
                 end do
              end do
           end do
           
           deallocate(intBuffer)

        else
           ! Do the local copy:
           oBlocks(iDom)%iBlank = flowDoms(iDom, 1, 1)%iBlank
        end if
     end do
     
  else
     
     ! Pack and send the iblank
     do nn=1, nDom
         call setPointers(nn, 1, 1)
         allocate(intBuffer((ib+1)*(jb+1)*(kb+1)))

         ii = 0
         do k=0, kb
            do j=0, jb
               do i=0, ib
                  ii = ii + 1
                   intBuffer(ii) = iBlank(i, j, k)
               end do
            end do
         end do

         tag = cumDomProc(myid) + nn
         
         call MPI_Send(intBuffer, size(intBuffer), sumb_integer, 0, tag, &
              sumb_comm_world, ierr)
         deallocate(intBuffer)
      end do
   end if

   
 end subroutine tmpUpdateIblanks
 
 subroutine test


 use overset
   use blockPointers
 implicit none

 integer(kind=intType) :: nn, i, j, k, mm
 real(kind=realType) :: val1, val2, coor(3)

 do nn=1,nDom

    call setPointers(nn, 1, 1)
    do k=1, ke
       do j=1, je
          do i=1, ie

             ! Get cell center coordinate
             coor = eighth*(&
                  x(i-1, j-1, k-1, :) + &
                  x(i  , j-1, k-1, :) + &
                  x(i-1, j  , k-1, :) + &
                  x(i  , j  , k-1, :) + &
                  x(i-1, j-1, k  , :) + &
                  x(i  , j-1, k  , :) + &
                  x(i-1, j  , k  , :) + &
                  x(i  , j  , k  , :))
             w(i, j, k, iVx) = coor(1) + 2*coor(2)
          end do
       end do
    end do
 end do
 !call  whalo2(1, 1, 5, .True., .True. , .True. )

 call wOverset(1, 2, 2, .False., .False., .False., .False.) 




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
 end subroutine test
