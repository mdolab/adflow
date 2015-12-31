!      ******************************************************************
!      *                                                                *
!      * computeOversetInterpolation is the top level routine that      *
!      * implements the implicit hole cutting method for determing      *
!      * overset grid connectivitiies. It operates on a given multigrid *
!      * level and spectral instance                                    *
!      *                                                                *
!      ******************************************************************

subroutine oversetComm(level, firstTime, coarseLevel)
  use constants
  use communication
  use blockPointers
  use overset
  use stencils
  use BCTypes
  use inputTimeSpectral
  use ADTapi
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: level
  logical :: firstTime, coarseLevel
  integer(kind=intType) :: sps

  ! Local Variables
  integer(kind=intType) :: i, ii, iii, j, jj, jjj, k, kk, kkk,  iDom, jDom, iDim
  integer(kind=intType) :: nCopy, m, iSize, oldSize
  integer(kind=intType) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
  integer(kind=intType) :: nn, mm, n, ierr, iProc
  integer(kind=intTYpe) :: rowStart, rowEnd, nnRow, nUniqueProc
  integer(kind=intType) :: sendCount, recvCount, currentProc, nLocalFringe
  integer(kind=intType) :: iWork, nWork, nFringeProc
  real(kind=realType) :: startTime, endTime, timeA, timeB
  real(kind=realType), dimension(:, :), allocatable :: xMin, xMax
  real(kind=realType), dimension(:), allocatable :: minVol
  integer(kind=intType), dimension(:), allocatable :: clusters, nProcSend
  integer(kind=intType), dimension(:), allocatable :: fringesReceivedFromProc
  integer(kind=intType), dimension(:,:), allocatable :: work
  integer(kind=intType), dimension(:), allocatable :: procsForThisRow, inverse
  integer(kind=intType), dimension(:), allocatable :: fringeProc, cumFringeProc
  logical, dimension(:), allocatable :: oBlockReady, fringesReady
  type(fringeType), dimension(:), allocatable :: localFringes

  type(CSRMatrix), pointer :: overlap
  type(CSRMatrix) :: overlapTranspose

  ! MPI/Communication related
  integer status(MPI_STATUS_SIZE) 
  logical :: barrierActive, barrierDone, flag
  integer :: barrierRequest
  integer(kind=intType) :: MAGIC, source, tag

  ! -----------------------------------------------------------------
  ! Step 0: Initializeion of communication data structures. 
  ! -----------------------------------------------------------------
  timeA = mpi_wtime()

  ! Determine the magic number which is actually the same as
  ! nDomTotal. nDomTotal is guaranteed to be greater than or equal to
  ! nProc. This will be used to space out tags on communications to
  ! make sure they do not overlap. 
  call mpi_allreduce(nDom, MAGIC, 1, sumb_integer, MPI_SUM, &
       sumB_comm_world, ierr)

  ! -----------------------------------------------------------------
  ! Step 1: Initializaion: We need to register our fringe type with
  ! MPI so we can use it directly in communication of fringes. This
  ! will make our life a lot simplier. We also need to initialize
  ! stencils since this may be first time they are used. 
  ! -----------------------------------------------------------------

  call initialize_stencils()
  call registerOversetFringeType()

  ! -----------------------------------------------------------------
  ! Step 2: Communicate the block size info to everyone. Also generate
  ! cumDomProc which is the cumulative form. This will make our lives
  ! easier a little later on with indexing. (Routine below)
  ! -----------------------------------------------------------------

  allocate(nDomProc(0:nProc-1), cumDomProc(0:nProc))
  call blockSizeInfo

  ! Master SPS loop
  spectralLoop: do sps=1,nTimeIntervalsSpectral

     ! Set a pointer to make the code easier to read
     overlap => overlapMatrix(level, sps)

     ! -----------------------------------------------------------------
     ! Step 3: Next we determine the "clusters" that are present in the
     ! grid. We will use the original CGNS information since this is
     ! slightly more readily available. Futhermore, there will generally
     ! be fewer blocks in the CGNS file that compute blocks due to
     ! splitting. Once we have labelled the CGNS blocks, we can simply
     ! label the compute blocks by looking at what CGNS block they came
     ! from. Since all procs have this information we can simply run on
     ! all procs. 
     ! -----------------------------------------------------------------
     allocate(clusters(nDomTotal))
     call determineClusters(clusters, nDomTotal, cumDomProc)

     ! -----------------------------------------------------------------
     ! Step 4: Compute the 3D axis oriented bounding boxes for each block
     ! and communicate then to everyone. Also determine communicate the
     ! minimum volume for each block to everyone.  (Routine below)
     ! -----------------------------------------------------------------

     allocate(xMin(3, nDomTotal), xMax(3, nDomTotal), minVol(nDomTotal))
     call computeDomainBoundingBoxes
     ! -----------------------------------------------------------------
     ! Step 8: Build a global sparse matrix representation of the overlap
     ! matrix.  Every processor will have the same sparse matrix
     ! representation when it is finished. (Routine below)
     ! -----------------------------------------------------------------

     call deallocateCSRMatrix(overlap)
     call buildGlobalSparseOverlap(overlap)

     ! Done with the clusters
     deallocate(clusters)

     ! -----------------------------------------------------------------
     ! Step 8: This is going to put the number of searches (coordinates)
     ! in for the costs.  Eventually we need to use the previous time,
     ! but that is going to be tricky since the sparsity structure of the
     ! overlap matrix could change....:-(
     ! -----------------------------------------------------------------

     overlap%data = zero
     ! Loop over by blocks and the owned cells on my blocks. This will
     ! determine which of my coordinates need to be searched on a given
     ! block. 

     do nn=1, nDom
        call setPointers(nn, level, sps)
        iDom = cumDomProc(myid) + nn

        do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
           overlap%data(jj) = nx*ny*nz
        end do
     end do

     ! Determine the total costs for everyone. 
     call mpi_allreduce(MPI_IN_PLACE, overlap%data, overlap%nnz, sumb_real, MPI_SUM, &
          sumb_comm_world, ierr)

     ! -----------------------------------------------------------------
     ! Step 8: We are now ready to partiaion and loadbalance the work
     ! based on the costs stored in the overlap matrix. These costs may
     ! be the search estimates from initializeOverlapCosts OR they may be
     ! actual timings from a previous assembly.
     ! -----------------------------------------------------------------

     call oversetLoadBalance(overlap)
     
     ! Sending the fringes will be easier with the transpose of the
     ! overlap matrix
     call transposeOverlap(overlap, overlapTranspose)

     ! Create a local array with just the information about the work I need to do:
     nWork = 0
     do jj=1,overlap%nnz
        if (overlap%assignedProc(jj) == myid) then 
           nWork = nWork + 1
        end if
     end do
     allocate(work(4, nWork))

     nWork = 0
     do iDom=1, nDomTotal
        do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
           jDom = overlap%colInd(jj)
           if (overlap%assignedProc(jj) == myID) then 
              nWork = nWork + 1
              work(1, nWork) = iDom
              work(2, nWork) = jDom
              work(3, nWork) = jj
              work(4, nWork) = 0
           end if
        end do
     end do

     ! Send/recv Requests: As an overestimate if my rows and columns are
     ! dense, I would need (nDom*nDomTotal).  We have at most 3 sets of
     ! messages in flight at a time so we need 6*nDom*nDomTotal. This
     ! will always be an overestimate and it is guaranteed to be larger
     ! than nProc since nDomTotal>nProc
     if (allocated(sendRequests)) deallocate(sendRequests)
     if (allocated(recvRequests)) deallocate(recvRequests)
     i = nDom*nDomTotal*6
     allocate(sendRequests(i), recvRequests(i))

     ! -----------------------------------------------------------------
     ! Step 8: Allocation of temporary data structures: oBlocks and
     ! fringeList. The oBlocks contains the ADTrees while the fringe list
     ! contains the coordines we need to search for.  Note that we
     ! allocate both to nDomTotal...which is not technically scalable,
     ! but there since there are only a few scattered variables and no
     ! large arrays it should be fine. 
     ! -----------------------------------------------------------------

     allocate(oBlocks(nDomTotal))
     allocate(fringeList(nDomTotal))

     ! Build all the trees for my local procs. Note that we call
     ! setPointers and oBlock is build with whatever is pointed to with
     ! blockPointers. Also build all the fringes. When we are finished
     ! building each, immediately send them. 

     sendCount = 0
     allocate(procsForThisRow(nDomTotal), inverse(nDomTotal))
     allocate(oBlockReady(iDom), fringesReady(iDom))
     oBlockReady = .False.
     fringesReady = .False.

     do nn=1, nDom
        call setPointers(nn, level, sps)
        iDom = cumDomProc(myid) + nn

        ! -----------------------------------------------------------------
        ! Step 8a: Initialize and pack the oBlock corresponding to the
        ! block that I own
        ! -----------------------------------------------------------------

        call initializeOBlock(oBlocks(iDom), nn)
        call packOBlock(oBlocks(iDom))
        oBlockReady(iDom) = .True.

        ! -----------------------------------------------------------------
        ! Step 8b: Send this block to all other processors that neeed it:
        ! -----------------------------------------------------------------

        ! Number of non-zeros in this row
        nnRow = overlap%rowPtr(iDom+1) - overlap%rowPtr(iDom)

        ! The processors for this row. Copy info the first nnRow entries.
        procsForThisRow(1:nnRow) = overlap%assignedProc(overlap%rowPtr(iDom) : overlap%rowPtr(iDom+1)-1)

        ! Unique-ify them since we could end up sending the same block to
        ! the same processor more than once, which is unncecessary. We
        ! can do this in place and nUniqueProc are the unique number of
        ! processors to send this block to. inverse is unused.
        call unique(procsForThisRow, nnRow, nUniqueProc, inverse)

        do jj = 1, nUniqueProc
           if (procsForThisRow(jj) /= myid) then 

              tag = 1*MAGIC + nn
              sendCount = sendCount + 1
              call mpi_issend(oBlocks(iDom)%rBuffer, size(oBlocks(iDom)%rBuffer), &
                   sumb_real, procsForThisRow(jj), tag, SUmb_comm_world, &
                   sendRequests(sendCount), ierr)
              call ECHK(ierr, __FILE__, __LINE__)

              tag = 2*MAGIC + nn
              sendCount = sendCount + 1
              call mpi_issend(oBlocks(iDom)%iBuffer, size(oBlocks(iDom)%iBuffer), &
                   sumb_integer, procsForThisRow(jj), tag, SUmb_comm_world, &
                   sendRequests(sendCount), ierr)
              call ECHK(ierr, __FILE__, __LINE__)

           end if
        end do

        ! -----------------------------------------------------------------
        ! Step 8c: Initialize the fringes originating from this
        ! processor's block
        ! -----------------------------------------------------------------

        ! Now add our owned cells only to the fringe list.
        allocate(fringeList(iDom)%arr(nx*ny*nz))
        mm = 0
        do k=2, kl
           do j=2, jl
              do i=2, il
                 mm = mm + 1

                 ! Set empty (default) values for this fringe. 
                 call emptyFringe(fringeList(iDom)%arr(mm))

                 do iDim=1, 3
                    fringeList(iDom)%arr(mm)%x(iDim) = eighth*(&
                         x(i-1, j-1, k-1, iDim) + &
                         x(i  , j-1, k-1, iDim) + &
                         x(i-1, j  , k-1, iDim) + &
                         x(i  , j  , k-1, iDim) + &
                         x(i-1, j-1, k  , iDim) + &
                         x(i  , j-1, k  , iDim) + &
                         x(i-1, j  , k  , iDim) + &
                         x(i  , j  , k  , iDim))
                 end do

                 fringeList(iDom)%arr(mm)%myBlock = nn
                 fringeList(iDom)%arr(mm)%myI = i
                 fringeList(iDom)%arr(mM)%myJ = j
                 fringeList(iDom)%arr(mm)%myK = k
                 fringeList(iDom)%arr(mm)%quality = vol(i, j, k)
              end do
           end do
        end do

        ! Now loop over this block's boundary condiitons and we need
        ! to set a litle info: We need to flag the two levels next to
        ! an overset outer boundary as being forceRecv. This will be
        ! used when determining donors since we need to accept a donor
        ! even if it is worse quality than my own cell. We also need
        ! to flag a single layer of cells next a wall boundary
        ! condition as being "isWall". Knowing the fringes next to
        ! walls will be necessary for determine the overap wall
        ! distance correction as well as the flooding procedure. 

        do mm=1,nBocos
           ! Just record the ranges necessary and we'll add in a generic loop
           select case (BCFaceID(mm))
           case (iMin)
              iStart=2; iEnd=3;
              jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
              kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
           case (iMax)
              iStart=nx; iEnd=il;
              jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
              kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
           case (jMin)
              iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
              jStart=2; jEnd=3;
              kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
           case (jMax)
              iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
              jStart=ny; jEnd=jl;
              kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
           case (kMin)
              iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
              jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
              kStart=2; kEnd=3;
           case (kMax)
              iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
              jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
              kStart=nz; kEnd=kl;
           end select

           if (BCType(mm) == OversetOuterBound) then
              do k=kStart, kEnd
                 do j=jStart, jEnd
                    do i=iStart, iEnd
                       ! Compute the index
                       kk = (k-2)*nx*ny + (j-2)*nx + (i-2) + 1
                       fringeList(iDom)%arr(kk)%forceRecv = .True.
                       fringeList(iDom)%arr(kk)%quality = large
                    end do
                 end do
              end do
           end if

           ! For the wall check, we only need 1 layer next to
           ! wall. Modify the required bounds for this
           select case (BCFaceID(mm))
           case (iMin)
              iEnd=2
           case (iMax)
              iStart=il
           case (jMin)
              jEnd=2
           case (jMax)
              jStart=jl
           case (kMin)
              kEnd=2
           case (kMax)
              kStart=kl
           end select

           if (BCType(mm) == NSWallAdiabatic .or. &
               BCType(mm) == NSWallIsoThermal .or. &
               BCType(mm) == EulerWall) then 
              do k=kStart, kEnd
                 do j=jStart, jEnd
                    do i=iStart, iEnd
                       ! Compute the index
                       kk = (k-2)*nx*ny + (j-2)*nx + (i -2) + 1
                       fringeList(idom)%arr(kk)%isWall = .True.
                    end do
                 end do
              end do
           end if
        end do ! BocoLoop

        ! Flag my own set of fringes as being ready
        fringesReady(iDom) = .True.

        ! -----------------------------------------------------------------
        ! Step 8d: Send these fringes to the processors that need
        ! them. Note we use the transpose overlap matrix here. 
        ! -----------------------------------------------------------------

        ! Number of non-zeros in this row
        rowEnd = overlapTranspose%rowPtr(iDom+1)
        rowStart = overlapTranspose%rowPtr(iDom)
        nnRow =  rowEnd - rowStart

        ! The processors for this row. Copy info the first nnRow entries.
        procsForThisRow(1:nnRow) = overlapTranspose%assignedProc(rowStart:rowEnd-1)

        ! Unique-ify them since we could end up sending the same fringes
        ! to the same processor more than once, which is unncecessary. We
        ! can do this in place and nUniqueProc are the unique number of
        ! processors to send this block to. inverse is unused.
        call unique(procsForThisRow, nnRow, nUniqueProc, inverse)

        do jj = 1, nUniqueProc
           if (procsForThisRow(jj) /= myid) then 

              ! This intersection requires fringes from me
              tag = 3*MAGIC + nn
              sendCount = sendCount + 1

              call mpi_issend(fringeList(iDom)%arr, size(fringeList(iDom)%arr), &
                   oversetMPIFringe, procsForThisRow(jj), tag, &
                   SUmb_comm_world, sendRequests(sendCount), ierr)
              call ECHK(ierr, __FILE__, __LINE__)
           end if
        end do
     end do

     ! -----------------------------------------------------------------
     ! Step 8d: Start of the second part of the NBX algorithm
     ! -----------------------------------------------------------------

     ! Done with the transposed matrix
     call deallocateCSRMatrix(overlapTranspose)

     ! Keep track of where we receive fringes from 
     allocate(fringesReceivedFromProc(nDomTotal))
     fringesReceivedFromProc = -1
     overlap%data = zero

     ! Allocate space for the localWallFringes.. localWallFringes should not be
     ! allocated yet.
     allocate(localWallFringes(1000))
     nLocalWallFringe = 0
     barrierDone = .False.
     barrierActive = .False.

     do while (.not. barrierDone)
        flag = .False.
        call MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, SUmb_comm_world, flag, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! Check if a message is ready
        messageReceived: if (flag) then

           tag = status(MPI_TAG)
           source = status(MPI_SOURCE)

           ! What we do depends on the tag: Note that we don't start
           ! the tags at 1 becuase we will eventually call
           ! exchangeFringes() and exchangeIBlanks(), both of which
           ! have iSends with tags from 0 to nProc (nProc <=
           ! MAGIC). There is *very* small chance that one of those
           ! communications could interfere with our NBX exchanges. To
           ! make sure this never happens, every comm set has a new and
           ! unique MAGIC range and iProbes check the tag to ensure we
           ! are receiving the message from the correct communiation. 

           !   MAGIC + 1 <= tag <= 2*MAGIC : real receive for oBlock
           ! 2*MAGIC + 1 <= tag <= 3*MAGIC : integer receive for oBlock
           ! 3*MAGIC + 1 <= tag <= 4*MAGIC : fringe receive

           if (tag >= MAGIC+1 .and. tag <= 2*MAGIC) then 

              ! Get size
              call MPI_Get_count(status, sumb_real, n, ierr)
              call ECHK(ierr, __FILE__, __LINE__)

              ! Back out the domain from the tag.
              iDom = cumDomProc(source) + (tag - MAGIC)

              ! Allocate space for the real receiver
              allocate(oBlocks(iDom)%rBuffer(n))

              ! Now actually receive the real buffer. Blocking is fine
              call mpi_recv(oBlocks(iDom)%rBuffer, n, sumb_real, source, &
                   tag, SUmb_comm_world, status, ierr)
              call ECHK(ierr, __FILE__, __LINE__)

           else if (tag >= 2*MAGIC+1 .and. tag <= 3*MAGIC) then 

              ! Get size
              call MPI_Get_count(status, sumb_integer, n, ierr)
              call ECHK(ierr, __FILE__, __LINE__)

              ! Back out the domain from the tag.
              iDom = cumDomProc(source) + (tag - 2*MAGIC)

              ! Allocate space for integer receiver
              allocate(oBlocks(iDom)%iBuffer(n))

              ! Now actually receive the integer buffer. Blocking is fine
              call mpi_recv(oBlocks(iDom)%iBuffer, n, sumb_integer, source, &
                   tag, SUmb_comm_world, status, ierr)
              call ECHK(ierr, __FILE__, __LINE__)

           else if (tag >= 3*MAGIC+1 .and. tag <= 4*MAGIC) then 

              ! Get size
              call MPI_Get_count(status, oversetMPIFringe, n, ierr)
              call ECHK(ierr, __FILE__, __LINE__)

              ! Back out the domain from the tag.
              iDom = cumDomProc(source) + (tag - 3*MAGIC)

              ! Allocate space for the integer receiver
              allocate(fringeList(iDom)%arr(n))

              ! Now actually receive the real buffer. Blocking is fine
              call mpi_recv(fringeList(iDom)%arr, n, oversetMPIFringe, status(MPI_SOURCE), &
                   tag, SUmb_comm_world, status, ierr)
              call ECHK(ierr, __FILE__, __LINE__)

              ! Save the proc number of where we got this from
              fringesReceivedFromProc(iDom) = status(MPI_SOURCE)

              ! Flag it to say we now have these fringes
              fringesReady(iDom) = .True.

           end if

           if (tag >= MAGIC+1 .and. tag <= 3*MAGIC) then 
              ! Have received either reals or integers for oBlock. Check
              ! if we now have *both* so that we can unpack the block and
              ! delclare it ready for searching. 

              if (allocated(oBlocks(iDom)%rBuffer) .and. allocated(oBlocks(iDom)%iBuffer)) then 
                 call unpackOBlock(oBlocks(iDom))           
                 oBlockReady(iDom) = .True.
              end if
           end if

        end if messageReceived

        ! -----------------------------------------------------------------
        ! Step 8: We can now loop over the number of searches we have to
        ! do and do them provided we have the oBlock and fringes we
        ! need. Note that we are looping over only the list of
        ! intersections we have to do, that is extracted from the
        ! full overlap matrix
        ! -----------------------------------------------------------------

        do iWork=1, nWork
           
           iDom = work(1, iWork)
           jDom = work(2, iWork)
           jj   = work(3, iWork) ! Index from the overlap matrix

           ! Check if I have the oBlock and fringes i need to do this
           ! intersection and I haven't already done it.
           if (oBlockReady(iDom) .and. fringesReady(jDom) .and. work(4, iWork) == 0) then 

              startTime = mpi_wtime()
              call fringeSearch(oBlocks(iDom), fringeList(jDom)%arr, size(fringeList(jDom)%arr))
              endTime = mpi_wtime()
              overlap%data(jj) = endTime - startTime

              ! Flag this work as being done:
              work(4, iWork) = 1

           end if
        end do

        ! Complete NBX algorithm
        if (.not. barrierActive) then
           call MPI_testAll(sendCount, sendRequests, flag, MPI_STATUSES_IGNORE, ierr)
           call ECHK(ierr, __FILE__, __LINE__)

           if (flag) then 
              call MPI_Ibarrier(sumb_comm_world, barrierRequest, ierr)
              call ECHK(ierr, __FILE__, __LINE__)
              barrierActive = .True.
           end if
        else
           call MPI_test(barrierRequest, barrierDone, MPI_STATUS_IGNORE, ierr)
           call ECHK(ierr, __FILE__, __LINE__)
        end if
     end do

     ! -----------------------------------------------------------------
     ! Step 9: Well, all the searches are done, so now we will send the
     ! fringes back to where they came from. We were smart and recorded
     ! the proc from where we recieved the fringes in the first place, so
     ! it is easy to just reverse the comm. Again we do the dynamic
     ! sparse data exchange.
     ! -----------------------------------------------------------------

     sendCount = 0
     do iDom=1, nDomTotal
        ! Check if we have recived this one. This would be indicated by
        ! fringesReceivedFromProc being different from the default value
        ! of -1
        if (fringesReceivedFromProc(iDom) /= -1) then 
           tag =  4*MAGIC + iDom
           sendCount = sendCount + 1
           call mpi_issend(fringeList(iDom)%arr, size(fringeList(iDom)%arr), &
                oversetMPIFringe, fringesReceivedFromProc(iDom), tag, &
                SUmb_comm_world, sendRequests(sendCount), ierr)
           call ECHK(ierr, __FILE__, __LINE__)
        end  if
     end do 

     barrierDone = .False.
     barrierActive = .False.

     ! Allocate tmpFringes to be an arbitrary size...it will get larger
     ! if necessary
     allocate(tmpFringes(1000))

     do while (.not. barrierDone)

        call MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, SUmb_comm_world, flag, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        tag = status(MPI_TAG)

        ! Check if a message is ready
        if (flag .and. tag >= 4*MAGIC+1 .and. tag <= 5*MAGIC) then 

           call receiveFringes(status, n)

           ! Back out the domain from the tag.
           iDom = (tag - 4*MAGIC)

           ! This is essentially the same logical operation that is done
           ! in fringeSearch, except this is now for the fringes that
           ! were on other processors that we are now receiving. 

           do i=1, n
              if (tmpFringes(i)%quality < fringeList(iDom)%arr(i)%quality) then 
                 ! Accept the new fringe. It is now no longer a compute cell
                 fringeList(iDom)%arr(i) = tmpFringes(i)
                 fringeList(iDom)%arr(i)%isCompute = .False.
              end if
           end do
        end if

        if (.not. barrierActive) then
           call MPI_testAll(sendCount, sendRequests, flag, MPI_STATUSES_IGNORE, ierr)
           call ECHK(ierr, __FILE__, __LINE__)

           if (flag) then 
              call MPI_Ibarrier(sumb_comm_world, barrierRequest, ierr)
              call ECHK(ierr, __FILE__, __LINE__)
              barrierActive = .True.
           end if
        else
           call MPI_test(barrierRequest, barrierDone, MPI_STATUS_IGNORE, ierr)
           call ECHK(ierr, __FILE__, __LINE__)
        end if
     end do

     ! -----------------------------------------------------------------
     ! Step 9: We now have computed all the fringes that we can. Some of
     ! them may be 'irrigular' as described in "A highly automated
     ! parallel Chimera method for overset grids based on the implicit
     ! hole cutting technique". If a cell is both a fringe and a donor,
     ! we need to force it to be a compute cell and remove it's receiver
     ! status. Do this as long as it isn't a forced receiver...we can't
     ! touch the forced receivers. 
     !
     ! How we do this is as follows: On each processor we gather up the
     ! fringes on each block on each processor that are actual fringes,
     ! ie, they have donorProc /= -1. That is, there is a
     ! potential donor that is better than its own cell. Once we get a
     ! list of these, we sort them, by processor, then block, then
     ! index. Then we can send this fringe information back to the
     ! processors where the donors came from.
     ! -----------------------------------------------------------------

     !-----------------------------------------------------------------
     ! Step 15: At this point we now have to take our fringe list and
     ! expand it out to a block-wise 3D array WITH HALOS. At this point,
     ! we have only been dealing with the owned cells in a 1D array
     ! which. By putting it in a 3D arry, we can communicate the fringe
     ! information to the halos, very easily using the same comm
     ! structure as the 1:1 exchange.
     ! -----------------------------------------------------------------

     nLocalFringe = 0
     do nn=1, nDom
        iDom = cumDomProc(myid) + nn

        call setPointers(nn, level, sps) ! Just to get sizes...
        allocate(flowDoms(nn, level, sps)%fringes(0:ib, 0:jb, 0:kb))

        ! Call again for the pointer to the fringes to be set
        call setPointers(nn, level, sps)

        ! Loop over all cells, (but we really just want the halo) setting
        ! fringes on the halos to not be compute cells. 
        ! 
        do k=0, kb
           do j=0, jb
              do i=0, ib
                 call emptyFringe(fringes(i, j, k))
                 fringes(i, j, k)%myI = i
                 fringes(i, j, k)%myJ = j
                 fringes(i, j, k)%myK = k
                 fringes(i, j, k)%myBlock = nn
                 fringes(i, j, k)%isCompute = .False.
              end do
           end do
        end do

        ! Now copy in the ones we have from the list into the owned cells
        mm = 0
        do k=2, kl
           do j=2, jl
              do i=2, il
                 mm = mm + 1
                 fringes(i, j, k) = fringeList(iDom)%arr(mm)
                 if (fringes(i, j, k)%donorProc /= -1) then 
                    nLocalFringe = nLocalFringe + 1
                 end if
              end do
           end do
        end do
     end do

     ! We can now deallocate the fringeList since we've copied all the
     ! info out and we no longer need it. 

     do iDom=1, nDomTotal
        if (allocated(fringeList(iDom)%arr)) then 
           deallocate(fringeList(iDom)%arr)
        end if
     end do

     deallocate(fringeList)

     ! Allocate the a new local 1D fringe list that just has our local
     ! fringes (from all local blocks) that are *actually* fringes. Do
     ! not include the halos.
     allocate(localFringes(nLocalFringe))

     ! Fill up these fringes
     nLocalFringe = 0
     do nn=1, nDom
        call setPointers(nn, level, sps)
        do k=2, kl
           do j=2, jl
              do i=2, il
                 if (fringes(i, j, k)%donorProc /= -1) then 
                    nLocalFringe = nLocalFringe + 1
                    localFringes(nLocalFringe) = fringes(i, j, k)
                 end if
              end do
           end do
        end do
     end do

     ! Now just sort the fringes such that they are grouped by destination procesor.
     call qsortFringeType(localFringes, nLocalFringe)

     !-----------------------------------------------------------------
     ! Step 15: Now can scan through the fringeList which is guaranteed
     ! to be sorted by processor. We scan though the (sorted) list
     ! sequentally, detecting where the processor splits are. Then we
     ! fire that off to the processor that needs it, again using the
     ! dynamic sparse data exchange. For the section that is
     ! on-processor, we can do that donor flagging overlapped with the
     ! communication. 
     ! -----------------------------------------------------------------

     allocate(fringeProc(nProc), cumFringeProc(1:nProc+1))
     call computeFringeProcArray(localFringes, nLocalFringe, &
          fringeProc, cumFringeProc, nFringeProc)

     ! nFringeProc is the total number of donor processors from localFringes
     ! fringeProc(1:nProcFringe) are the donor processors for the fringes
     ! cumFringeProc(1:nFringeProc) are the cumulative offset from in the localFringe Array
     ! Note that we have over-estimated the size as nProc.

     sendCount = 0
     do j=1, nFringeProc

        iProc = fringeProc(j)
        iStart = cumFringeProc(j)
        iEnd = cumFringeProc(j+1)-1

        if (iProc == myid) then 
           do i=iStart, iEnd
              nn = localFringes(i)%donorBlock

              do kk=0, 1
                 do jj=0, 1
                    do ii=0, 1

                       iii = ii + localFringes(i)%dI
                       jjj = jj + localFringes(i)%dJ
                       kkk = kk + localFringes(i)%dK

                       flowDoms(nn, level, sps)%fringes(iii, jjj, kkk)%isDonor = .True.
                    end do
                 end do
              end do
           end do

        else

           ! I need to send these fringes to the donor proc
           sendCount = sendCount + 1

           ! Use a new set of tags starting at 5*MAGIC
           tag = 5*MAGIC + iProc + 1
           call mpi_issend(localFringes(iStart:iEnd), iEnd-iStart+1, &
                oversetMPIFringe, iProc, tag, &
                SUmb_comm_world, sendRequests(sendCount), ierr)
           call ECHK(ierr, __FILE__, __LINE__)
        end if
     end do

     barrierDone = .False.
     barrierActive = .False.
     do while (.not. barrierDone)

        call MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, SUmb_comm_world, flag, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        tag = status(MPI_TAG)
        source = status(MPI_SOURCE)

        ! Check if a message is ready
        if (flag .and. tag >= 5*MAGIC+1 .and. tag <= 6*MAGIC) then 

           call receiveFringes(status, n)

           ! Do the flagging of the donors on the fly here:
           do j=1, n
              nn = tmpFringes(j)%donorBlock
              do kk=0, 1
                 do jj=0, 1
                    do ii=0, 1
                       iii = ii + tmpFringes(j)%dI
                       jjj = jj + tmpFringes(j)%dJ
                       kkk = kk + tmpFringes(j)%dK

                       flowDoms(nn, level, sps)%fringes(iii, jjj, kkk)%isDonor = .True.
                    end do
                 end do
              end do
           end do
        end if

        if (.not. barrierActive) then
           call MPI_testAll(sendCount, sendRequests, flag, MPI_STATUSES_IGNORE, ierr)
           call ECHK(ierr, __FILE__, __LINE__)

           if (flag) then 
              call MPI_Ibarrier(sumb_comm_world, barrierRequest, ierr)
              call ECHK(ierr, __FILE__, __LINE__)
              barrierActive = .True.
           end if
        else
           call MPI_test(barrierRequest, barrierDone, MPI_STATUS_IGNORE, ierr)
           call ECHK(ierr, __FILE__, __LINE__)
        end if
     end do

     ! Done with the localFringe array for now
     deallocate(localFringes)

!==================================================================================

     ! -----------------------------------------------------------------
     ! Step 10: We also have to send the wall fringes to their
     ! respective donor procs. We use a very similar communication
     ! structure as we used for the regular fringes. The difference
     ! here is that we unlike the regular case we keep *all* fringes
     ! we receive and just put them in a big list:
     ! wallFringes. wallFringes is the list of fringes ON THE DONOR
     ! processor whos receiver is a wall point. This is precisely the
     ! information we need to start the flooding process later on. 
     ! -----------------------------------------------------------------

     call determineWallDonors(level, sps, MAGIC)

!=================================================================================
     ! -----------------------------------------------------------------
     ! Step 10: We can now locally perform the irregular cell correction
     ! by looping over the fringes on my proc and just checking if
     ! donorProc is not -1 and isDonor is True.  If so, we force it back
     ! to be compute, by cancelling the donor information. Update the
     ! fringes when we're done so everyone has up to date information.
     ! -----------------------------------------------------------------
     call exchangeDonorStatus(level, sps, commPatternCell_2nd, internalCell_2nd)
     
     call irregularCellCorrection(level, sps)

     call exchangeFringes(level, sps, commPatternCell_2nd, internalCell_2nd)

     ! Next we have to perfrom the interior cell flooding. We already
     ! have the information we need: we have isWallFringe defined in
     ! the fringes as well as knowing if a cell is a compute. We
     ! should probably only flood compute cells that are not also
     ! donors, since that would get a little complicated. 

     call floodInteriorCells(level, sps)

     call exchangeFringes(level, sps, commPatternCell_2nd, internalCell_2nd)

     !-----------------------------------------------------------------
     ! Step 15: Reduction of the number of fringes. What we do is look at
     ! all the fringes and see if all the cells in its stencil are also
     ! fringes or holes. If so we can flag that particular cell as a
     ! hole.
     ! -----------------------------------------------------------------
     
     call fringeReduction(level, sps)

     call exchangeFringes(level, sps, commPatternCell_2nd, internalCell_2nd)

     ! -----------------------------------------------------------------
     ! Step 17: We can now create the final required comm structures
     ! based on our interpolation. This is relatively straight forward:
     ! The fringes we have local are the cells the "receiver" cells. We
     ! sort those by donor processor (as we did before) and then send
     ! them to donor proc. This then form the sending information. The
     ! internal copy is formed from the part that is on-processor. 
     ! -----------------------------------------------------------------

     call finalOversetCommStructures(level, sps, MAGIC)
     
     ! VERY last thing is to update iBlank based on the status of our local fringes. 
     call setIblankArray(level, sps)

     ! -----------------------------------------------------------------
     ! Step 16: The algorithm is now complete. Run the checkOverset
     ! algorithm to verify that we actually have a valid interpolation
     ! -----------------------------------------------------------------
     call checkOverset(level, sps)

     ! Deallocate some data we no longer need
     deallocate(Xmin, Xmax, minVol, fringesReceivedFromProc, &
          tmpFringes, oBlockReady, fringesReady, work, &
          procsForThisRow, inverse, fringePRoc, cumFringeProc)

     ! oBlocks is a little tricker since we need to deallocate the
     ! individual parts first:
     do iDom=1, nDomTotal
        if (oblocks(iDom)%allocated) then 
           deallocate(&
                oBlocks(iDom)%hexaConn, &
                oBlocks(iDom)%globalCell, &
                oBLocks(iDOm)%nearWall, &
                oBlocks(iDom)%qualDonor, &
                oBlocks(iDom)%xADT, &
                oBlocks(iDom)%rBuffer, &
                oBlocks(iDom)%iBuffer)
           call destroySerialHex(oBlocks(iDom)%ADT)
        end if
     end do

     ! Finally deallocate the oBlocks array itself
     deallocate(oBlocks)

     ! Setup the buffer sizes
     call setBufferSizes(level, sps, .false., .false., .true.)

  end do spectralLoop

  ! Free the buffer and make a new one that includes necessary sizes
  ! for the overset comm
  deallocate(sendBuffer, recvBuffer)
  allocate(sendBuffer(sendBufferSize), recvBuffer(recvBufferSize))
  
  print *,' DONE! interpolation', myid, mpi_wtime()-timeA
  deallocate(nDomProc, cumDomProc)

contains

  ! Simple utility-type routines that make the main subroutine
  ! easier to read

  subroutine blockSizeInfo

    implicit none

    ! Gather the dimensions of all blocks to everyone
    call mpi_allreduce(nDom, nDomTotal, 1, sumb_integer, MPI_SUM, &
         sumb_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    ! Receive the number of domains from each proc using an allgather.
    call mpi_allgather(nDom, 1, sumb_integer, nDomProc, 1, sumb_integer, &
         sumb_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    ! Compute the cumulative format:
    cumDomProc(0) = 0
    do iProc=1, nProc
       cumDomProc(iProc) = cumDomProc(iProc-1) + nDomProc(iProc-1)
    end do

  end subroutine blockSizeInfo

  subroutine computeDomainBoundingBoxes

    implicit none

    ! Working Variables
    real(kind=realType), dimension(3, nDom) :: xMinLocal, xMaxLocal
    real(kind=realType), dimension(nDom) :: minVolLocal

    do nn=1,nDom
       call setPointers(nn, level, sps)

       xMinLocal(1, nn) = minval(x(:, :, :, 1))
       xMinLocal(2, nn) = minval(x(:, :, :, 2))
       xMinLocal(3, nn) = minval(x(:, :, :, 3))

       xMaxLocal(1, nn) = maxval(x(:, :, :, 1))
       xMaxLocal(2, nn) = maxval(x(:, :, :, 2))
       xMaxLocal(3, nn) = maxval(x(:, :, :, 3))

       minVolLocal(nn) = minval(vol(2:il, 2:jl, 2:kl))

    end do

    ! Now we can allgather the xMin, xMax and minVolume from each
    ! processor to everyone
    call mpi_allgatherV(xMinLocal, nDom*3, sumb_real, xMin, 3*nDomProc, &
         3*cumDomProc, sumb_real, sumb_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    call mpi_allgatherV(xMaxLocal, nDom*3, sumb_real, xMax, 3*nDomProc, &
         3*cumDomProc, sumb_real, sumb_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    call mpi_allgatherV(minVolLocal, nDom, sumb_real, minVol, nDomProc, &
         cumDomProc, sumb_real, sumb_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

  end subroutine computeDomainBoundingBoxes

  subroutine buildGLobalSparseOverlap(overlap)

    use overset
    implicit none

    ! Input/Output
    type(CSRMatrix), intent(inout) :: overlap

    ! Working
    integer(kind=intType), dimension(:), allocatable :: colIndLocal, rowPtrLocal
    logical , dimension(:, :), allocatable :: localOverlap
    integer(kind=intType) :: nnzLocal
    integer(kind=intType), dimension(:), allocatable :: nnzProc, cumNNzProc

    ! Allocate the space for my processor's compoent of iOverlap. The
    ! number of rows is the number of blocks I own (in my own numbering)
    ! and the coluns are the total number of blocks. It is initialized
    ! to True, and remains that way unless we can conclusively prove the
    ! blocks don't overlap. 
    allocate(localOverlap(nDom, nDomTotal))
    localOverlap = .True.

    ! Assume all rows are dense, and decrement when we eliminate an
    ! entry
    nnzLocal = nDom*nDomTotal

    ! Bounding box check: We will check our *own* blocks against the
    ! bounding boxes for for each domain.  This will fill up the
    ! localOverlapArray and determine nnzLocal -- the number of
    ! non-zero local entries in the global sparse matrix.

    do nn=1, nDom
       iDom = cumDomProc(myid) + nn

       ! Now Loop over *all* of the other blocks
       do jDom=1, nDomTotal 

          ! We can eliminate some of pairs using the cluser
          ! analysis. Note that the cluster analysis will take care
          ! of making sure that a block doesn't interset itself. ie,
          ! the diagonal of the full matrix must be .False.

          if (clusters(iDom) == clusters(jDom)) then 
             localOverlap(nn, jDom) = .False.
             nnzLocal = nnzLocal - 1
          end if

          ! Only do the spatial check if we haven't elminated the
          ! connection through the cluster check
          if (localOverlap(nn, jDom)) then 

             ! Now do the box overlap check
             if ( &
                  xMin(1, iDom) >= xMax(1, jDom) .or. &
                  xMax(1, iDom) <= xMin(1, jDom) .or. &
                  xMin(2, iDom) >= xMax(2, jDom) .or. &
                  xMax(2, iDom) <= xMin(2, jDom) .or. &     
                  xMin(3, iDom) >= xMax(3, jDom) .or. &
                  xMax(3, iDom) <= xMin(3, jDom)) then 

                ! These bounding boxes do not intersect. 
                localOverlap(nn, jDom) = .False.
                nnzLocal = nnzLocal -1 
             end if
          end if
       end do
    end do

    ! Now, create a sparse matrix representation of the local part

    allocate(colIndLocal(nnzLocal), rowPtrLocal(nDom + 1))
    rowPtrLocal(1) = 1
    i = 0
    do nn=1, nDom
       do jDom=1,nDomTotal
          if (localOverlap(nn, jDom)) then 
             i = i + 1
             colIndLocal(i) = jDom
          end if
       end do
       rowPtrLocal(nn+1)  = i + 1 
    end do

    ! Now we want to assemble the global (sparse!) connectivity
    ! matrix for all processors. This is going to require a little
    ! communication of sizes first followed by the actual data.

    ! Determine distribution of non-zero locations
    allocate(nnzProc(nProc), cumNNZProc(0:nProc))
    call mpi_allgather(nnzLocal, 1, sumb_integer, nnzProc, 1, sumb_integer, &
         sumb_comm_world, ierr)

    overlap%nnz = sum(nnzProc)
    overlap%nRow = nDomTotal
    overlap%nCol = nDomTotal
    overlap%nnzLocal = nnzLocal
    ! We can now know how big the global data will be
    allocate(&
         overlap%data(overlap%nnz), &
         overlap%colInd(overlap%nnz), &
         overlap%rowPtr(overlap%nRow + 1), &
         overlap%assignedProc(overlap%nnz))
    overlap%allocated = .True.
    cumNNZProc(0) = 0
    do iProc=1,nProc
       cumNNZProc(iProc) = cumNNZProc(iProc-1) + nnzProc(iProc)
    end do

    ! Correct for rowPtrLocal for the cumulative number of nnz up to
    ! this proc:
    rowPtrLocal = rowPtrLocal + cumNNZProc(myid)

    ! Gather the column indicies
    call mpi_allgatherV(colIndLocal, nnzLocal, sumb_integer, overlap%colInd, &
         nnzProc, cumNNZProc, sumb_integer, sumb_comm_world, ierr)

    ! Now we gather the rowPtr to everone
    overlap%rowPtr(1) = 1
    call mpi_allgatherV(rowPtrLocal(2:nDom), nDom, sumb_integer, &
         overlap%rowPtr(2:nDomTotal+1), &
         nDomProc, cumDomProc, sumb_integer, sumb_comm_world, ierr)

    deallocate(colIndLocal, rowPtrLocal, localOverlap)

  end subroutine buildGlobalSparseOverlap
end subroutine oversetComm

! ----------------------------------------------------------------
! Other subroutines
! ----------------------------------------------------------------

subroutine initializeOBlock(oBlock, nn)

  ! This routine allocates the data for the supplied oBlock using the
  !  data currently in blockPointers
  use constants
  use overset
  use blockPointers
  use adtAPI
  use BCTypes
  use cgnsGrid
  use communication
  implicit none 

  ! Input Params
  type(oversetBlock), intent(inout) :: oBlock
  integer(kind=intType) :: nn, kk

  ! Working paramters
  integer(kind=intType) :: i, j, k, mm, nADT, nHexa, planeOffset
  integer(kind=intType) :: iStart, iEnd, jStart, jEnd, kStart, kEnd

  ! Set all the sizes for this block.
  oBlock%il = il
  oBlock%jl = jl
  oBlock%kl = kl

  oBlock%ie = ie
  oBlock%je = je
  oBlock%ke = ke

  oBlock%ib = ib
  oBlock%jb = jb
  oBlock%kb = kb

  oBlock%nx = nx
  oBlock%ny = ny
  oBlock%nz = nz

  oBlock%proc = myID
  oBlock%block = nn

  allocate( &
       oBlock%qualDonor(1, ie*je*ke), &
       oBlock%globalCell(0:ib, 0:jb, 0:kb), &
       oBlock%nearWall(1:ie, 1:je, 1:ke))
  oBlock%nearWall = 0

  kk = 5
  do mm=1,nBocos
     select case (BCFaceID(mm))
     case (iMin)
        iStart=1; iEnd=2+kk;
        jStart=BCData(mm)%icBeg; jEnd=BCData(mm)%icEnd
        kStart=BCData(mm)%jcBeg; kEnd=BCData(mm)%jcEnd
     case (iMax)
        iStart=ie-kk; iEnd=ie;
        jStart=BCData(mm)%icBeg; jEnd=BCData(mm)%icEnd
        kStart=BCData(mm)%jcBeg; kEnd=BCData(mm)%jcEnd
     case (jMin)
        iStart=BCData(mm)%icBeg; iEnd=BCData(mm)%icEnd
        jStart=1; jEnd=1+kk;
        kStart=BCData(mm)%jcBeg; kEnd=BCData(mm)%jcEnd
     case (jMax)
        iStart=BCData(mm)%icBeg; iEnd=BCData(mm)%icEnd
        jStart=je-kk; jEnd=je;
        kStart=BCData(mm)%jcBeg; kEnd=BCData(mm)%jcEnd
     case (kMin)
        iStart=BCData(mm)%icBeg; iEnd=BCData(mm)%icEnd
        jStart=BCData(mm)%jcBeg; jEnd=BCData(mm)%jcEnd
        kStart=1; kEnd=1+kk;
     case (kMax)
        iStart=BCData(mm)%icBeg; iEnd=BCData(mm)%icEnd
        jStart=BCData(mm)%jcBeg; jEnd=BCData(mm)%jcEnd
        kStart=ke-kk; kEnd=ke;
     end select

     if (BCType(mm) == NSWallAdiabatic .or. &
          BCType(mm) == NSWallIsoThermal .or. &
          BCType(mm) == EulerWall) then 

        do k=kStart, kEnd
           do j=jStart, jEnd
              do i=iStart, iEnd
                 oBlock%nearWall(i, j, k) = 1
              end do
           end do
        end do
     end if
  end do ! BocoLoop

  ! Copy Volume to qualDonor and do minVol while we're at it
  oBlock%minVol = Large
  mm = 0
  do k=1,ke
     do j=1,je
        do i=1,ie
           mm = mm + 1
           oBlock%qualDonor(1, mm) = vol(i, j, k)
           oBlock%minVol = min(oBlock%minVol, vol(i, j, k))
        end do
     end do
  end do
  oBlock%globalCell = globalCell

  ! Now setup the data for the ADT
  nHexa = il * jl * kl
  nADT = ie * je * ke

  allocate(oBlock%xADT(3, nADT), oBlock%hexaConn(8, nHexa))

  ! Fill up the xADT using cell centers (dual mesh)
  mm = 0
  do k=1, ke
     do j=1, je
        do i=1, ie
           mm = mm + 1
           oBlock%xADT(:, mm) = eighth*(&
                x(i-1, j-1, k-1, :) + &
                x(i  , j-1, k-1, :) + &
                x(i-1, j  , k-1, :) + &
                x(i  , j  , k-1, :) + &
                x(i-1, j-1, k  , :) + &
                x(i  , j-1, k  , :) + &
                x(i-1, j  , k  , :) + &
                x(i  , j  , k  , :))
        end do
     end do
  end do

  mm = 0
  ! These are the 'elements' of the dual mesh.
  planeOffset = ie * je
  do k=2, ke
     do j=2, je
        do i=2, ie
           mm = mm + 1
           oBlock%hexaConn(1, mm) = (k-2)*planeOffset + (j-2)*ie + (i-2) + 1
           oBlock%hexaConn(2, mm) = oBlock%hexaConn(1, mm) + 1 
           oBlock%hexaConn(3, mm) = oBlock%hexaConn(2, mm) + ie
           oBlock%hexaConn(4, mm) = oBlock%hexaConn(3, mm) - 1 

           oBlock%hexaConn(5, mm) = oBlock%hexaConn(1, mm) + planeOffset
           oBlock%hexaConn(6, mm) = oBlock%hexaConn(2, mm) + planeOffset
           oBlock%hexaConn(7, mm) = oBlock%hexaConn(3, mm) + planeOffset
           oBlock%hexaConn(8, mm) = oBlock%hexaConn(4, mm) + planeOffset
        end do
     end do
  end do

  ! Call the custom build routine -- Serial only, only Hexa volumes,
  ! we supply our own ADT Type

  call buildSerialHex(nHexa, nADT, oBlock%xADT, oBlock%hexaConn, oBlock%ADT)

  ! Flag this block as being allocated
  oBlock%allocated = .True.

end subroutine initializeOBlock

  
subroutine test
  use overset
  use blockPointers
  implicit none

  integer(kind=intType) :: nn, i, j, k, mm
  real(kind=realType) :: val1, val2, coor(3)

  do nn=1,1

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
              w(i, j, k, iVx) = coor(1) + 2*coor(2) + 3*coor(3)
           end do
        end do
     end do
  end do

  call setPointers(2, 1, 1)
  w(:, :, :, iVx) = 0.0

  call whalo2(1, 1, 5, .False., .False. , .False. )

end subroutine test


subroutine writePartionedMesh(fileName)

  ! This is a debugging routine for writing out meshes *as they are
  ! partioned*. This can be useful for debugging overset issues.

  use communication
  use blockPointers

  implicit none

  character(len=*), intent(in) :: fileName
  integer(kind=intType) :: nDomTotal, iProc, nn, i, j, k, iDim, iDom, ierr, ii
  integer(kind=intType) :: bufSize, maxSize, ibufSize, imaxSize
  integer(kind=intType), dimension(3, nDom) :: localDim
  integer(kind=intType), dimension(:), allocatable :: nDomProc, cumDomProc
  integer(kind=intType), dimension(:, :), allocatable :: dims
  real(kind=realType), dimension(:), allocatable :: buffer
  integer(kind=intType), dimension(:), allocatable :: ibuffer !iblank
  character*40 :: tmpStr

  integer status(MPI_STATUS_SIZE) 
  ! Gather the dimensions of all blocks to everyone
  call mpi_allreduce(nDom, nDomTotal, 1, sumb_integer, MPI_SUM, &
       sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Store the sizes of the local blocks
  do nn=1,nDom
     call setPointers(nn, 1, 1)

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

  maxSize = 0
  imaxSize = 0
  do i=1,nDomTotal
     maxSize = max(maxSize, dims(1, i)*dims(2,i)*dims(3,i)*3)
     imaxSize = max(imaxSize, dims(1, i)*dims(2,i)*dims(3,i))
  end do

  allocate(buffer(maxSize))
  allocate(ibuffer(imaxSize)) !iblank

  if (myid == 0) then 
     print *,'writing mesh...'
     ! Root proc does all the writing. Just dump to ascii tecplot
     ! file---really slow.

     open(unit=1, file=fileName, form='formatted', status='unknown')
     !write(1, *) "Variables = X Y Z"
     write(1, *) "Variables = X Y Z IBLANK"

     ! Write my own blocks first
     do nn=1,nDom
        call setPointers(nn, 1, 1)
        !write(tmpStr, *) "Proc ", 0, " Local ID", nn
        write(tmpStr, "(a,I2.2,a,I2.2,a)"), """Proc ", 0, " Local ID ", nn , """"
        write(1,*) "ZONE I=", il, " J=",jl, "K=", kl, "T=", trim(tmpStr)
        write(1, *) "DATAPACKING=BLOCK"
        !write(1, *) "VARLOCATION=([1,2,3]=NODAL)"
        write(1, *) "VARLOCATION=([1,2,3,4]=NODAL)"

        do iDim=1, 3
           do k=1, kl
              do j=1, jl
                 do i=1, il
                    write(1, *) x(i, j, k, idim)
                 end do
              end do
           end do
        end do

        ! Iblanks are cell center values,  imposed on primal nodes, for plotting
        ! purpose only
        do k=1, kl
           do j=1, jl
              do i=1, il
                 write(1, *) iBlank(i+1, j+1, k+1) 
              end do
           end do
        end do
     end do

     ! Now loop over the remaining blocks...receiving each and writing:

     do iProc=1, nProc-1
        do nn=1, nDomProc(iProc)
           iDom = cumDomProc(iProc) + nn
           bufSize = dims(1, iDom)*dims(2, iDom)*dims(3,iDom)*3
           ibufSize = dims(1, iDom)*dims(2, iDom)*dims(3,iDom)

           call MPI_Recv(buffer, bufSize, sumb_real, iProc, iProc, &
                sumb_comm_world, status, ierr)

           call MPI_Recv(ibuffer, ibufSize, sumb_integer, iProc, 10*iProc, &
                sumb_comm_world, status, ierr)

           write(tmpStr, "(a,I2.2,a,I2.2,a)"), """Proc ", iProc, " Local ID ", nn ,""""
           write(1,*) "ZONE I=", dims(1, iDom), " J=", dims(2, iDom), "K=", dims(3, iDom), "T=", trim(tmpStr)
           write(1, *) "DATAPACKING=BLOCK"
           !write(1, *) "VARLOCATION=([1,2,3]=NODAL)"
           write(1, *) "VARLOCATION=([1,2,3,4]=NODAL)"

           ! Dump directly...already in the right order
           do i=1, bufSize
              write(1, *), buffer(i)
           end do
           do i=1, ibufSize
              write(1, *), ibuffer(i)
           end do
        end do
     end do
     close(1)
     print *,'Done writing mesh...'
  else 

     ! Pand and send my stuff:
     do nn=1, nDom
        call setPointers(nn, 1, 1)
        ii = 0
        do iDim=1, 3
           do k=1, kl
              do j=1, jl
                 do i=1, il
                    ii = ii + 1
                    buffer(ii) = x(i, j, k, iDim)
                 end do
              end do
           end do
        end do

        call mpi_send(buffer, ii, sumb_real, 0, myid, &
             sumb_comm_world, ierr)

        ! Iblanks are cell center values,  imposed on primal nodes, for plotting
        ! purpose only
        ii = 0
        do k=1, kl
           do j=1, jl
              do i=1, il
                 ii = ii + 1
                 ibuffer(ii) = iBlank(i+1, j+1, k+1)
              end do
           end do
        end do
        call mpi_send(ibuffer, ii, sumb_integer, 0, 10*myid, &
             sumb_comm_world, ierr)
     end do
  end if

  deallocate(buffer, ibuffer, nDomProc, cumDomProc, dims)

end subroutine writePartionedMesh
