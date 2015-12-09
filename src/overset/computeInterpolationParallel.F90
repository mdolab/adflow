!      * computeOversetInterpolation is the top level routine that      *
!      * implements the implicit hole cutting method for determing      *
!      * overset grid connectivitiies. For now it only operates on the  *
!      * finest grid and first spectral instance.                       *
!      *                                                                *
!      ******************************************************************

subroutine computeOversetInterpolationParallel

  use constants
  use communication
  use blockPointers
  use overset
  use stencils
  use BCTypes
  implicit none

  ! Local Variables
  integer(kind=intType) ::i, j, k,  level, sps, ii, iDom, jDom, iDim, nn, mm, tmp
  integer(kind=intType) ::  jj, kk, ip2, im2, jp2, jm2, kp2, km2, ind(3), n
  integer(kind=intType) ::  iii, jjj, kkk, nfrng
  integer(kind=intType) :: donorBlockID, i_stencil, rSize, iSize, jj1
  integer(kind=intType) :: startDom, endDom, nblkRecv, source, nblkprocRecv
  real(kind=realType), dimension(:, :), allocatable :: xMin, xMax
  real(kind=realType), dimension(:), allocatable :: minVol
  logical, dimension(:, :), allocatable :: localOverlap
  logical, dimension(:), allocatable :: receiveBlocks, ownReceiveBlocks
  logical, dimension(:, :), allocatable :: sendBlocks
  logical, dimension(:), allocatable :: touchOBlocks
  logical, dimension(:), allocatable :: usedBlocks
  integer(kind=intType), dimension(:), allocatable :: recvBlocks
  logical, dimension(:, :), allocatable :: invrecvBlocks
  integer(kind=intType), dimension(:), allocatable :: invsendBlocks 

  real(kind=realType), dimension(3, nDom) :: xMinLocal, xMaxLocal
  real(kind=realType), dimension(nDom) :: minVolLocal

  integer(kind=intType), dimension(3, nDom) :: localDim
  integer(kind=intTYpe) :: iProc, ierr, blockID
  integer(kind=intType), dimension(:), allocatable :: clusters, blkProc, blkLocalID
  real(kind=realType), dimension(:, :), pointer :: real2DPtr
  integer(kind=intType), dimension(:, :), pointer :: int2DPtr

  ! Buffers variables
  integer(kind=intType) :: nIntSend, tag
  real(kind=realType), dimension(:), allocatable :: realSendBuffer, realRecvBuffer
  integer(kind=intType), dimension(:), allocatable :: intSendBuffer, intRecvBuffer
  integer(kind=intTYpe), dimension(:, :, :), allocatable :: iblankTmp
  integer(kind=intType), dimension(1:nDom+1) :: realSendOffset, intSendOffset
  integer(kind=intType), dimension(:), allocatable :: realRecvOffset, intRecvOffset
  integer(kind=intType) :: index, sendCount, recvCount, totalRecv
  integer status(MPI_STATUS_SIZE) 
  logical :: computeCellFound, needBlk
  real(kind=realType) :: timeA, timeB, xpt(3), time1, time2

  integer(kind=intType) :: oldSize, newSize

  integer(kind=intType), dimension(:), allocatable :: nnzProc, cumNNzProc
  integer(kind=intType), dimension(:), allocatable :: colIndLocal, rowPtrLocal
  integer(Kind=intType) :: nnzLocal
  ! Final global CSR matrix of overlap
  type(CSRMatrix) :: overlap

  ! Explictly set level and sps to 1. This will be removed in the future.
  level = 1
  sps = 1 

  ! This may be the first routine that uses stencils so initialize
  call initialize_stencils()

  ! -----------------------------------------------------------------
  ! Step 1: Determine which cells are inside the computational
  ! domain. Essential we determine the cells that must have iblank=0.
  ! This is done in parallel using a hybrid of the ray-cast and min
  ! distance method.
  ! -----------------------------------------------------------------

  call computeHolesInsideBody
  call MPi_barrier(sumb_comm_world, ierr)
  timeA = mpi_wtime()
  ! -----------------------------------------------------------------
  ! Step 2: We will communcation all the block size information. We
  ! also generate cumDomProc which is cumulative form of the number
  ! of blocks on each processor. This will make our lives easier a
  ! little later on.
  ! -----------------------------------------------------------------

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

  ! -----------------------------------------------------------------
  ! Step 5: Next we determine the "clusters" that are present in the
  ! grid. We will use the original CGNS information since this is
  ! slightly more readily available. Futhermore, there will generally
  ! be fewer blocks in the CGNS file that compute blocks due to
  ! splitting. Once we have labelled the CGNS blocks, we can simply
  ! label the compute blocks by looking at what CGNS block they came
  ! from. Since all procs have this information we can simply run on
  ! all procs
  ! -----------------------------------------------------------------
  allocate(clusters(nDomTotal))
  call determineClusters(clusters, nDomTotal, cumDomProc)

  ! -----------------------------------------------------------------
  ! Step 3: Here we will compute the inverse mapping of of which
  ! processor a global block belongs to as well what the the local
  ! index is for that block on its owner processor. Again, this
  ! information will make the code cleaner, easier to understand later
  ! on. 
  ! -----------------------------------------------------------------
  allocate(blkProc(nDomTotal), blkLocalID(nDomTotal))
  ii = 0
  do iProc=0, nProc-1
     do i=1, nDomProc(iProc)
        ii = ii + 1
        blkProc(ii) = iProc
        blkLocalID(ii) = i
     end do
  end do

  ! -----------------------------------------------------------------
  ! Step 4: Bounding box check + minVolume communication: Communicate
  ! all the xMin, xMax bounds + min block volume for all blocks to
  ! everyone.
  ! -----------------------------------------------------------------

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

  ! Allocate the space we need for the numbers and cumulative
  allocate(xMin(3, nDomTotal), xMax(3, nDomTotal), minVol(nDomTotal))

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

  ! -----------------------------------------------------------------
  ! Step 6: Bounding box check: We will check check our *own* blocks
  ! against the bounding boxes for for each domain. 
  ! -----------------------------------------------------------------

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

  ! Loop over my all blocks
  do nn=1, nDom
     iDom = cumDomProc(myid) + nn

     ! Now Loop over *all* of the other blocks
     do jDom=1, nDomTotal 

        ! We can eliminate some of pairs using the cluser analysis:
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

  ! Done with the clusters
  deallocate(clusters)


  !a ! Debug output for localOverlap (still kinda square-ish) 
  !a if (myid ==  0) then 
  !a    print *, '------------------------------' 
  !a end if
  !a do iProc=0, nProc-1
  !a    if (iProc == myid) then 
  !a       do nn=1,nDom 
  !a          call setPointers(nn, 1,1) 
  !a          print *, nn+cumDomProc(myid), ' ', localOverlap(nn, :) 
  !a       end do
  !a       print *, '------------------------------'
  !a       call sleep(1) !am
  !a    end if
  !a    call MPI_barrier(sumb_comm_world, ierr)
  !a end do

  ! Need this barrier here
  call MPI_barrier(sumb_comm_world, ierr)

  ! -----------------------------------------------------------------
  ! Step 7: Build a sparse matrix representation of the local part of
  ! the overlap matrix. The number of rows is nDom, the number of
  ! columns is nDomTotal. We now know the number of non-zero entries
  ! as well which is great. 
  ! -----------------------------------------------------------------

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

  ! -----------------------------------------------------------------
  ! Step 7: Now we want to assemble the global (sparse!) connectivity
  ! matrix for all processors. This is going to require a little
  ! communication of sizes first followed by the actual data. 
  ! -----------------------------------------------------------------

  ! Determine distribution of non-zero locations
  allocate(nnzProc(nProc), cumNNZProc(0:nProc))
  call mpi_allgather(nnzLocal, 1, sumb_integer, nnzProc, 1, sumb_integer, &
       sumb_comm_world, ierr)

  overlap%nnz = sum(nnzProc)
  overlap%nRow = nDomTotal
  overlap%nCol = nDomTotal
  ! We can now know how big the global data will be
  allocate(&
       overlap%data(overlap%nnz), &
       overlap%colInd(overlap%nnz), &
       overlap%rowPtr(overlap%nRow + 1), &
       overlap%assignedProc(overlap%nnz))

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

  ! The overlap matrix only tells us which blocks (potentially)
  ! overlap with each other. That's not a whole lot of
  ! information. Here we try to generate a little more information by
  ! estimating how many searches we will have to do for each
  ! block. 

  overlap%data = zero
  ! Loop over by blocks and the owned cells on my blocks. This will
  ! determine which of my coordinates need to be searched on a given
  ! block. 

  startDom = cumDomProc(myid) + 1
  endDom   = cumDomProc(myid+1)

  do nn=1, nDom
     call setPointers(nn, level, sps)
     iDom = cumDomProc(myid) + nn

     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
        overlap%data(jj) = nx*ny*nz
     end do
  end do

  !a timeB = mpi_wtime()
  !a print *,'myid time 1:', myid, timeB-timeA

  ! Determine the total costs for everyone. 
  call mpi_allreduce(MPI_IN_PLACE, overlap%data, overlap%nnz, sumb_real, MPI_SUM, &
       sumB_comm_world, ierr)


  ! Now run the load balancing algorithm.
  call oversetLoadBalance(overlap)

  ! ! Overwrite so that there is no ADT movement
  ! do iProc=0,nProc-1
  !    overlap%assignedProc(cumNNZProc(iProc)+1:cumNNZProc(iProc+1)) = iProc
  ! end do

  ! -----------------------------------------------------------------
  ! Step 8: Allocation of oBlocks data structure
  ! -----------------------------------------------------------------

  ! Note that we allocate oBlocks to nDomTotal...which is not
  ! technically scalable, but there since there are only a few
  ! scattered integers in the type, it should be fine. We only
  ! allocate the array's for the blocks that are going to be sent to
  ! us.
  allocate(oBlocks(nDomTotal))

  ! ! Build all the trees for my local procs. Note that we call
  ! ! setPointers and oBlock is build with whatever is pointed to with
  ! ! blockPointers. 
  i = 0
  do nn=1,nDom
     call setPointers(nn, level, sps)
     i = i + il * jl * kl
     iDom = cumDomProc(myid) + nn
     call initializeOBlock(oBlocks(iDom))
     call packOBlock(oBlocks(iDom))
  end do
  
  !a if (myid == 0) then 
  !a    ! Now dump out who owns what:
  !a    do i=1,nDomTotal
  !a       write(*, "(a,I4, a)", advance='no'), 'Row:', i, "   "
  !a       do jj=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
  !a          write(*, "(a,I2, a, I6)", advance='no'), "(", overlap%colInd(jj), ")", int(overlap%data(jj))
  !a       end do
  !a       write(*, *) " "
  !a    end do
  !a    
  !a    print *, '--------------------------------------'
  !a    ! Now dump out who owns what:
  !a    do i=1,nDomTotal
  !a       write(*, "(a,I4, a)", advance='no'), 'Row:', i, "   "
  !a       do jj=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
  !a          write(*, "(a,I2, a, I6)", advance='no'), "(", overlap%colInd(jj), ")", int(overlap%assignedProc(jj))
  !a       end do
  !a       write(*, *) " "
  !a    end do
  !a end if

  ! -----------------------------------------------------------------
  ! On-block preprocessing: Flag cells near holes as
  ! forceRecv as well as cells that already have iblank=-1
  ! This info will be used to set 'qualRecv' in searchCoords intialization
  ! -----------------------------------------------------------------

  ! Initializes donor derived type of the compute block, i.e. 
  ! flowdoms(nn, level, sps)%donors, and forceRecv, recvStatus, qualRecv etc.
  do nn=1, nDom
     call allocateOverset(nn, level, sps)
  end do

  do nn=1,nDom
     call setPointers(nn, level, sps)

     ! We need to flag the status of cells according to the
     ! oversetOuterBoundary condition.
     ! iBlank refers to block pointed to in setPointers

     do mm=1,nBocos
        if(BCType(mm) == OversetOuterBound) then 
           select case (BCFaceID(mm))
           case (iMin)
              iBlank(1:3, :, :) = -1
           case (iMax)
              iBlank(nx:ie, :, :) = -1
           case (jMin)
              iBlank(:, 1:3, :) = -1
           case (jMax)
              iBlank(:, ny:je, :) = -1
           case (kMin)
              iBlank(:, :, 1:3) = -1
           case (kMax)
              iBlank(:, :, nz:ke) = -1
           end select
        end if
     end do

     do k=2, kl
        do j=2, jl
           do i=2, il

              select case (iBlank(i, j, k))

              case (0) 

                 ! We have to be careful with holes: We need to make
                 ! sure that there is sufficient frige padding around
                 ! them. Therefore for each real cell we loop over the
                 ! cells in its stencil. If a cell its stencil is a
                 ! hole, then the cell in question MUST be forced to
                 ! be a receiver. Note there is no issue with bounds,
                 ! since iBlank is double haloed and we're only
                 ! looping over the owned cells. 

                 stencil_Loop: do i_stencil=1,  N_visc_drdw
                    ii = visc_drdw_stencil(i_stencil, 1) + i
                    jj = visc_drdw_stencil(i_stencil, 2) + j
                    kk = visc_drdw_stencil(i_stencil, 3) + k

                    if (iblank(ii, jj, kk) /= 0) then 
                       forceRecv(i, j, k) = .True. 
                    end if
                 end do stencil_Loop
              case (-1)
                 ! Cell as ben flagged as necessairly being a
                 ! receiver so force it to be. 
                 forceRecv(i, j, k) = .True. 
              end select
           end do
        end do
     end do
  end do !nn=1, nDom



  ! -----------------------------------------------------------------
  ! Step 8: Allocation of searchCoords structure
  ! -----------------------------------------------------------------
  
  ! Like oBlocks, this is allocated to size nDomTotal and we allocate
  ! only own blocks. Any additional coordinates that will be
  ! communicated to us will be allocated only as necessary.
  
  allocate(searchCoords(nDomTotal))

  i = 0
  do nn=1,nDom
     call setPointers(nn, level, sps)
     iDom = nn + cumDomProc(myid)
     call initializeSearchCoords(searchCoords(iDom))
  end do

  ! Now for the *really* fun part...doing the communication of the
  ! blocks....there are blocks going *EVERYWHERE*!. We already have
  ! the blocks packed and we know the offsets to the send buffer for
  ! each of ourowned blocks. We can now use the assignedProc structure
  ! to determine who I need to send blocks to, and who I need to
  ! receive blocks from. 

  if (allocated(sendRequests)) deallocate(sendRequests)
  if (allocated(recvRequests)) deallocate(recvRequests)

  ! The send/recv request arrays are nnz. This is quite a large
  ! overestimate. Try to fix this later. 
  allocate(sendRequests(overlap%nnz), recvRequests(overlap%nnz))

  ! -----------------------------------------------------------------
  ! Step 8: Send the oBlocks as necessary
  ! -----------------------------------------------------------------

  ! Problem with the matrix approach here is that a single processor
  ! may need the same search coordintes *twice*. Since we don't want
  ! to send the same information twice, we compact back being with
  ! respect to processor. 
  
  allocate(sendBlocks(nDom, 0:nProc-1))
  sendBlocks = .False. 

  ! Loop over the rows I own...that is the ADTrees I have built
  do nn=1,nDom

     iDom = cumDomProc(myid) + nn 
     
     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
        jDom = overlap%colInd(jj)
        if (overlap%assignedProc(jj) /= myid) then

           ! This intersection is computed by a different processor
           ! than me.  I need to send my tree to him.
           sendBlocks(nn, overlap%assignedProc(jj)) = .True.

        end if
     end do
  end do

  ! This will store the processor of the blocks I need to receive. 
  allocate(recvBlocks(nDomTotal))
  recvBlocks = -1

  ! Loop over all the rows
  do iDom=1, nDomTotal
     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
        jDom = overlap%colInd(jj)

        if (overlap%assignedProc(jj) == myid) then ! I have to do this one.
        
           ! If I don't already own this block....I need to receive
           ! it. Use blkProc to know when it is coming from.
           if (iDom < startDom .or. iDom > endDom) then 
              !recvBlocks(iDom) = blkProc(jDom) !bug?
              recvBlocks(iDom) = blkProc(iDom)
           end if
        end if
     end do
  end do

  nBlkRecv = 0
  do iDom=1,nDomTotal
     if (recvBlocks(iDom)/=-1) then 
        nBlkRecv = nBlkRecv + 1
     end if
  end do
  !a timeB = mpi_wtime()
  !a print *,'myid before comm:', myid, timeB-timeA

  sendCount = 0
  do nn=1, nDom
     ! block in gloabl ordering. 
     iDom = cumDomProc(myid) + nn  
     
     do iProc = 0, nProc-1
        if (sendBlocks(nn, iProc)) then 
           tag = 2*(iDom-1) + 1
           sendCount = sendCount + 1

           call mpi_isend(oBlocks(iDom)%rBuffer, size(oBlocks(iDom)%rBuffer), &
                sumb_real, iProc, tag, SUmb_comm_world, &
                sendRequests(sendCount), ierr)
           call ECHK(ierr, __FILE__, __LINE__)

           tag = 2*(iDom-1) + 2
           sendCount = sendCount + 1
           call mpi_isend(oBlocks(iDom)%iBuffer, size(oBlocks(iDom)%iBuffer), &
                sumb_integer, iProc, tag, SUmb_comm_world, &
                sendRequests(sendCount), ierr)
           call ECHK(ierr, __FILE__, __LINE__)
        end if
     end do
  end do

 ! Now post all the receives...these block, but we use MPI_ANY_SOURCE
  ! and MPI_ANY_TAG

  !twice because of two kinds of (int and real) sends/recvs?
  do i=1, nBlkRecv*2
     
     ! Probe the message
     call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, sumb_comm_world, status, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
     tag = status(MPI_TAG)
     source = status(MPI_SOURCE)

     ! Real receive because the tag is ODD
     if (mod(tag, 2) == 1) then 

        ! Get size
        call MPI_Get_count(status, sumb_real, n, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! Back out the domain from the tag.
        iDom = (tag - 1)/2 + 1
        
        ! Allocate space for the integer receiver
        allocate(oBlocks(iDom)%rBuffer(n))

        ! Now actually receive the real buffer. Blocking is fine
        call mpi_recv(oBlocks(iDom)%rBuffer, n, sumb_real, source, &
             tag, SUmb_comm_world, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

     else ! Even tag so it is an integer receive
        
        ! Get size
        call MPI_Get_count(status, sumb_integer, n, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! Back out the domain from the tag.
        iDom = (tag-2)/2 + 1

        ! Allocate space for the integer receiver
        allocate(oBlocks(iDom)%iBuffer(n))

        ! Now actually receive the integer buffer. Blocking is fine
        call mpi_recv(oBlocks(iDom)%iBuffer, n, sumb_integer, source, &
             tag, SUmb_comm_world, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)
     end if
  end do

  ! Wait for all the send requests for finish. 
  do i=1, sendCount
     call mpi_waitany(sendCount, sendRequests, index, status, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
  enddo

  ! Now we can unpack the blocks that were received:
  do iDom=1, nDomTotal
     if (recvBlocks(iDom) /= -1) then
        call unpackOBlock(oBlocks(iDom))
     end if
  end do

  ! I think we actually need a barrier here.
  call MPi_barrier(sumb_comm_world, ierr)

  ! -----------------------------------------------------------------
  !              Communicate the required coordinates
  ! -----------------------------------------------------------------

  ! Reinitialzie sendBlocks
  sendBlocks = .False. 

  ! Loop over everything
  do iDom=1, nDomTotal
     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
        jDom = overlap%colInd(jj)

        if (overlap%assignedProc(jj) /= myid) then

           if (jDom >= startDom .and. jDom <= endDom) then 

              ! This intersection requires nodes from me
              nn = jDom - cumDomProc(myid)
              sendBlocks(nn, overlap%assignedProc(jj)) = .True.
              
           end if
        end if
     end do
  end do

  ! Reinitialzie recvBlocks
  recvBlocks = -1

  ! Loop over all the rows
  do iDom=1, nDomTotal
     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
        jDom = overlap%colInd(jj)

        if (overlap%assignedProc(jj) == myid) then ! I have to do this one.
        
           ! If I don't already own this block....I need to receive
           ! it. Use blkProc to know when it is coming from.
           if (jDom < startDom .or. jDom > endDom) then 
              recvBlocks(jDom) = blkProc(jDom)
           end if
        end if
     end do
  end do

  nBlkRecv = 0
  do iDom=1,nDomTotal
     if (recvBlocks(iDom)/=-1) then 
        nBlkRecv = nBlkRecv + 1
     end if
  end do

  sendCount = 0
  do nn=1, nDom
     ! block in gloabl ordering. 
     iDom = cumDomProc(myid) + nn  
     
     do iProc = 0, nProc-1
        if (sendBlocks(nn, iProc)) then 
           ! earlier only x-coord was sent
           ! -----------------------------------------
           !tag = iDom
           !sendCount = sendCount + 1

           !call mpi_isend(searchCoords(iDom)%x, searchCoords(iDom)%n*3, &
           !     sumb_real, iProc, tag, SUmb_comm_world, &
           !     sendRequests(sendCount), ierr)
           !call ECHK(ierr, __FILE__, __LINE__)
           ! -----------------------------------------

           !odd send is x-coord
           tag = 2*(iDom-1) + 1
           sendCount = sendCount + 1

           call mpi_isend(searchCoords(iDom)%x, searchCoords(iDom)%n*3, &
                sumb_real, iProc, tag, SUmb_comm_world, &
                sendRequests(sendCount), ierr)
           call ECHK(ierr, __FILE__, __LINE__)


           !even send is qualRecv
           tag = 2*(iDom-1) + 2
           sendCount = sendCount + 1

           call mpi_isend(searchCoords(iDom)%qualRecv, searchCoords(iDom)%n, &
                sumb_real, iProc, tag, SUmb_comm_world, &
                sendRequests(sendCount), ierr)
           call ECHK(ierr, __FILE__, __LINE__)
        end if
     end do
  end do

 ! Now post all the receives...these block, but we use MPI_ANY_SOURCE
  ! and MPI_ANY_TAG

  ! twice nBlkRecv because of two receives, x and qualRecv
  do i=1, nBlkRecv*2
     
     ! Probe the message
     call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, sumb_comm_world, status, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
     tag = status(MPI_TAG)
     iDom = tag
     source = status(MPI_SOURCE)

     ! receive search x-coord because the tag is ODD
     if (mod(tag, 2) == 1) then

        ! Get size
        call MPI_Get_count(status, sumb_real, n, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! Back out the domain from the tag.
        iDom = (tag - 1)/2 + 1

        ! Allocate space for the receive
        allocate(&
             searchCoords(iDom)%x(3, n/3), &
             !searchCoords(iDom)%dInd(3, n/3), &
             searchCoords(iDom)%dInd(4, n/3), & !store donorBlockId + donor indices
             searchCoords(iDom)%gInd(8, n/3), &
             searchCoords(iDom)%frac(3, n/3))
        searchCoords(iDom)%n = n/3
        searchCoords(iDom)%arrSize = n/3

        ! Now actually receive the real buffer. Blocking is fine
        call mpi_recv(searchCoords(iDom)%x, n, sumb_real, source, &
             tag, SUmb_comm_world, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

     else  ! Even tag is search qualRecv
  
        ! Get size
        call MPI_Get_count(status, sumb_real, n, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! Back out the domain from the tag.
        iDom = (tag-2)/2 + 1

        ! Allocate space for the receive qualRecv
        ! Note: n is different in this EVEN tag than the ODD tag n value
        allocate(searchCoords(iDom)%qualRecv(1,n))

        ! Now actually receive the real buffer. Blocking is fine
        call mpi_recv(searchCoords(iDom)%qualRecv, n, sumb_real, source, &
             tag, SUmb_comm_world, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

     end if !tag
  end do !nblkRecv

  ! Wait for all the send requests for finish. 
  do i=1, sendCount
     call mpi_waitany(sendCount, sendRequests, index, status, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
  enddo

  ! Update invalidDonor info for all OBlocks I am assigned (myid) using
  ! cellStatus arrays which were read in unpackOBlock and which I already know.
  ! cellStatus arrays for blocks that I already know/have were defined in
  ! initializeOBlock.

  ! ***
  ! Note that every proc assigned an OBlock does the same update,
  ! although does it in parallel. Need to do it before packing blocks?
  ! ***
 
  ! Avoid repeat update on the same OBlocks I am assigned
  allocate(touchOBlocks(nDomTotal))
  touchOblocks = .False.

  do iDom=1, nDomTotal
     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
        if (overlap%assignedProc(jj) == myid .and. .not.touchOBlocks(iDom)) then 
        
           touchOBlocks(iDom) = .TRUE.

           do k=1, oBlocks(iDom)%ke
              do j=1, oBlocks(iDom)%je
                 do i=1, oBlocks(iDom)%ie

                 select case (oBlocks(iDom)%cellStatus(i, j, k))

                 case (0) 
                    ! Blanked Cell: Ensure it  cannot be a donor

                    oBlocks(iDom)%invalidDonor(i, j, k) = 1 !integer 

                 case (-1)
                    ! Cell has been flagged as necessairly being a
                    ! receiver. We therefore prevent it from being a donor

                    oBlocks(iDom)%invalidDonor(i, j, k) = 1 !integer
                 end select

                 end do !i
              end do !j
           end do !k

        end if !assigned proc
     end do !jj
  end do !All domain loop

  ! Do not need touchOblocks any more
  deallocate(touchOBlocks)

  ! Finally we can do the searches:
  overlap%data = zero
  do iDom=1, nDomTotal
     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
        jDom = overlap%colInd(jj)
        if (overlap%assignedProc(jj) == myid) then 
           time1 = mpi_wtime()
           ! globalBlockId of donor: iDom
           call newSearch(oBlocks(iDom), searchCoords(jDom), iDom, jDom)
           time2 = mpi_wtime()
           overlap%data(jj) = time2-time1
        end if
     end do
  end do
  !a timeB = mpi_wtime()
  !a print *,'myid time 3:', myid, timeB-timeA

  ! This barrier is needed
  call mpi_barrier(sumb_comm_world, ierr)

  ! -----------------------------------------------------------------
  ! Update donorstatus now from oBlocks(iDom) that were updated from
  ! newSearch.
  ! Reverse communication of oBlock, of only oBlock()%donorStatus
  ! -----------------------------------------------------------------

  ! Find out the procs who used my oBlocks blocks that I owned

  ! Initialize invsendBlocks. I could have got assigned the oBlocks from
  ! all the nDomTotal blocks. This will store the original proc info of those blocks.
  allocate(invsendBlocks(nDomTotal))
  invsendBlocks = -1


  ! Loop over all the rows
  do iDom=1, nDomTotal
     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
        jDom = overlap%colInd(jj)

        if (overlap%assignedProc(jj) == myid) then ! I did this one.
        
           ! If I don't already own this block....I need to send it back
           ! to blkProc(iDom), from where it is came from. NOT JDOM
           if (iDom < startDom .or. iDom > endDom) then 
              invsendBlocks(iDom) = blkProc(iDom)
           end if
        end if
     end do
  end do

  ! Initialize invrecvBlocks. I own max nDom and max nDom only will be
  ! received from other assign procs. Further, I could have sent my oBlocks
  ! to max all assigned procs (0:nProc-1), from where I will receive them back.
  allocate(invrecvBlocks(nDom,0:nProc-1))
  invrecvBlocks = .False.

  nblkprocRecv = 0
  ! Loop over the rows I own...that is the ADTrees I have built
  do nn=1,nDom

     iDom = cumDomProc(myid) + nn 
     
     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
        jDom = overlap%colInd(jj)
        if (overlap%assignedProc(jj) /= myid) then

           ! This intersection was computed by a different processor
           ! than me.  I need to get back my tree from him.
           invrecvBlocks(nn, overlap%assignedProc(jj)) = .True.
  
        end if
     end do
  end do

  nblkprocRecv = 0
  do nn=1, nDom
     do iProc=0, nProc-1
        if (invrecvBlocks(nn, iProc))  nblkprocRecv = nblkprocRecv + 1
     end do
  end do

  ! Now do the send to the compute blocks owner of procs from all the copies of
  ! oBlocks from all the assigned procs
   
  sendCount = 0
  do iDom=1, nDomTotal
     
     if (invsendBlocks(iDom) /= -1) then

        ! I was assigned iDom but I didn't own it
        
        iProc = invsendBlocks(iDom) !send from myid to this iProc
        tag = iDom

        sendCount = sendCount + 1
        
        call mpi_isend(oBlocks(iDom)%donorStatus, size(oBlocks(iDom)%donorStatus), &
             MPI_LOGICAL, iProc, tag, SUmb_comm_world, &
             sendRequests(sendCount), ierr)
        call ECHK(ierr, __FILE__, __LINE__)
        !print*,'myid to iProc, size: ',myid, iProc, size(oBlocks(iDom)%donorStatus)
 
     end if
  end do !iDom=1, nDomTotal 

  ! Now receive oBlocks. I, the compute block proc could receive from 
  ! multiple procs (except for myid) for multiple blocks that I own. 
  ! Use MPI_ANY_SOURCE, MPI_ANY_TAG
  
  loop_invrecv1: do iii=1, nblkprocRecv

     ! Probe message to iProc
     call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, sumb_comm_world, status, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
     tag = status(MPI_TAG)
     source = status(MPI_SOURCE)

     do nn=1, nDom
        iDom = nn + cumDomProc(myid) 

        if (iDom == tag) then

           ! Get size
           call MPI_Get_count(status, MPI_LOGICAL, n, ierr)
           call ECHK(ierr, __FILE__, __LINE__)

           call setPointers(nn, level, sps)

           !print*,'myid from iProc, size: ',myid, source, n, size(flowDoms(nn, level, sps)%donorStatus)

           if (n /= size(donorStatus)) then
              print*, 'Fatal error: received size is not same as flowDoms()%donorStatus '
              stop
           end if
           ! Now actually receive the real buffer. Blocking is fine
           call mpi_recv(donorStatus, n, MPI_LOGICAL, source, &
                tag, SUmb_comm_world, status, ierr)
           call ECHK(ierr, __FILE__, __LINE__)

        end if

     end do

  end do loop_invrecv1

  ! Wait for all the send requests for finish. 
  do i=1, sendCount
     call mpi_waitany(sendCount, sendRequests, index, status, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
  enddo

!am  ! Determine the total costs for everyone. 
!am  call mpi_allreduce(MPI_IN_PLACE, overlap%data, overlap%nnz, sumb_real, MPI_SUM, &
!am       sumB_comm_world, ierr)
!am
!am  call oversetLoadBalance(overlap)
!am
!am  if (myid == 0) then 
!am     ! Now dump out who owns what:
!am     do i=1,nDomTotal
!am        write(*, "(a,I4, a)", advance='no'), 'Row:', i, "   "
!am        do jj=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
!am           write(*, "(a,I2, a, F6.3)", advance='no'), "(", overlap%colInd(jj), ")", overlap%data(jj)
!am        end do
!am        write(*, *) " "
!am     end do
!am     
!am     print *, '--------------------------------------'
!am     ! Now dump out who owns what:
!am     do i=1,nDomTotal
!am        write(*, "(a,I4, a)", advance='no'), 'Row:', i, "   "
!am        do jj=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
!am           write(*, "(a,I2, a, I6)", advance='no'), "(", overlap%colInd(jj), ")", int(overlap%assignedProc(jj))
!am        end do
!am        write(*, *) " "
!am     end do
!am
!am     do i=1,overlap%nnz
!am        print *,overlap%data(i)
!am     end  do
!am
!am  end if


  ! This barrier is very important
  ! ------------------------------
  call MPi_barrier(sumb_comm_world, ierr)
  ! ------------------------------

  ! ---------------------------------------------------------------------------
  ! Step 9: Send the donor info back to actual procs owning search coordinates
  ! ---------------------------------------------------------------------------

  ! Do inverse of searchCoords communication
  ! One global block can be used by many procs for overlap intersections
  ! Now, receive copies of searchCoords from all those assigned procs
  ! and update my own local copy with each of those assigned procs copies.
  ! Need updated info of dInd, frac, gInd, donorblockID (?) and qualRecv
  ! I already know x and fInd info from my own initializeSearchCoords

  ! Find out the procs who used my searchCoords blocks that I owned

  ! Initialize invsendBlocks. I could have got assigned the searchCoords from
  ! all the nDomTotal blocks. This will store the original proc info of those blocks.
  invsendBlocks = -1

  ! Loop over all the rows
  do iDom=1, nDomTotal
     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
        jDom = overlap%colInd(jj)

        if (overlap%assignedProc(jj) == myid) then ! I did this one.
        
           ! If I don't already own this block....I need to send it back
           ! to blkProc(jDom), from where it is came from.
           if (jDom < startDom .or. jDom > endDom) then 
              invsendBlocks(jDom) = blkProc(jDom)
           end if
        end if
     end do
  end do

  ! Initialize invrecvBlocks. I own max nDom and max nDom only will be
  ! received from other assign procs. Further, I could have sent my searchCoords
  ! to max all assigned procs (0:nProc-1), from where I will receive them back.
  invrecvBlocks = .False.

  ! Loop over everything
  do iDom=1, nDomTotal
     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
        jDom = overlap%colInd(jj)

        if (overlap%assignedProc(jj) /= myid) then

           if (jDom >= startDom .and. jDom <= endDom) then 

              ! This intersection required nodes from me
              ! Now I want them back
              nn = jDom - cumDomProc(myid)
              invrecvBlocks(nn, overlap%assignedProc(jj)) = .True.

           end if
        end if
     end do
  end do

  nblkprocRecv = 0
  do nn=1, nDom
     do iProc=0, nProc-1
        if (invrecvBlocks(nn, iProc))  nblkprocRecv = nblkprocRecv + 1
     end do
  end do

  ! Now do the sends to the compute blocks owner procs from all the copies of 
  ! searchCoords from all the assigned procs

  ! Send/recvs are done in two phases. 
  !
  ! 1st phase send/recvs just the qualRecvs and fracs. Min qualRecvs/fracs received 
  !  are saved, and so are the global block info of the corresponding copy of
  !  searcCoords from where qualRecvs were received. A 'special array' of current
  !  searchCoords size is created in my proc which stores the assigned 'proc'  
  !  of the searchCoords copy received.
  !
  ! 2nd phase send/recvs then goes through the recvs and picks out only the
  !  dInds/gInds from only those procIds of searchCoords received which
  !  are noted on the 'special array'. It picks out dInds/gInds corresponding
  !  to the current fInds.
   

  ! ---------------------------------------------------------------------
  ! 1st phase of sends of searchCoords copies from assigned procs to the
  ! original owners of the searchCoords (or to the compute blocks)
  ! ---------------------------------------------------------------------

  sendCount = 0
  do iDom=1, nDomTotal
     
     if (invsendBlocks(iDom) /= -1) then 
       
        ! I was assigned iDom but I didn't own it
        
        iProc = invsendBlocks(iDom) !send from myid to this iProc

        ! Define buffer 
        ! -------------------------------------------
        n = searchCoords(iDom)%n

        ! Allocate buffer
        rSize = size(searchCoords(iDom)%frac) + size(searchCoords(iDom)%qualRecv)
        allocate(searchCoords(iDom)%rBuffer(rSize))

        rSize = 0

        ! Real buffer
        ! Save frac 
        do j=1, n
           do i=1, 3
              rSize = rSize + 1
              searchCoords(iDom)%rBuffer(rSize) = searchCoords(iDom)%frac(i,j)
           end do
        end do

        ! Save qualRecv 
        do j=1, n
           rSize = rSize + 1
           searchCoords(iDom)%rBuffer(rSize) = searchCoords(iDom)%qualRecv(1,j)
        end do

        ! --- end define buffer ---------------------

        ! Real buffer
        tag = 1000*iProc + 2*(iDom-1) + 1   
        sendCount = sendCount + 1

        call mpi_isend(searchCoords(iDom)%rBuffer, size(searchCoords(iDom)%rBuffer), &
             sumb_real, iProc, tag, SUmb_comm_world, &
             sendRequests(sendCount), ierr)
        call ECHK(ierr, __FILE__, __LINE__)
       
     end if
  end do !nDomTotal

  ! 1st phase recvs: now only rBuffer (qualRecv/frac)
  !
  ! Now receive searchCoords. I, the compute block proc could receive from 
  ! multiple procs (except for myid) for multiple blocks that I own. 
  ! Use MPI_ANY_SOURCE, MPI_ANY_TAG

  ! Loop over my owned searchCoords and receive only those that were sent from 
  ! assigned procs

  ! Create and allocate 'special array' to save assigned 'proc' info of the
  ! incoming searchCoords. I own nDom blocks.
  allocate(specialArray(nDom))
  do nn=1, nDom
     iDom = nn + cumDomProc(myid)

     n = searchCoords(iDom)%n
     !allocate(specialArray(nn)%procId(3,n)) !size(fInd) !BUGS?
     allocate(specialArray(nn)%procId(1,n)) 

     ! Initialize to default value. Means it doesn't receive any searchCoords
     specialArray(nn)%procId = -1
  end do
  
  ! Before receiving begin by updating flowDoms from my own searchCoords for 
  ! potential donors. I already have access to the corresponding dInd/gInd too.
  ! -----------------------------------------------------------------

  ! See if my blocks that I owned was assigned to me for checking donor
  ! search for any blocks in the overlap matrix
  allocate(usedBlocks(nDom))
  usedBlocks = .False.
  ! Loop over everything
  do iDom=1, nDomTotal
     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
        jDom = overlap%colInd(jj)

        if ( (overlap%assignedProc(jj) == myid) .and. &
             (jDom >= startDom .and. jDom <= endDom) ) then

           ! This intersection required nodes from me
           ! and I owned this jDom block
           nn = jDom - cumDomProc(myid) !my local block order
           usedBlocks(nn) = .True.
        end if
     end do
  end do

  do nn=1, nDom
     ! cycle if my block was not used for my assigned overlap intersections
     if (.not.usedBlocks(nn)) cycle

     call setPointers(nn, level, sps) !points to flowDoms arrays
     iDom = nn + cumDomProc(myid) !my global block order

     do j=1, searchCoords(iDom)%n
        ! Current receiver cell's indices
        ii = searchCoords(iDom)%fInd(1,j)
        jj = searchCoords(iDom)%fInd(2,j)
        kk = searchCoords(iDom)%fInd(3,j)

        if (searchCoords(iDom)%qualRecv(1,j) < qualRecv(ii, jj, kk) .or. &
            (forceRecv(ii, jj, kk) .and. &
             searchCoords(iDom)%qualRecv(1, j) < large) ) then 

           ! if a the donor quality is better than existing one or
           ! if a donor has been found for the forced receiver already

           qualRecv(ii, jj, kk) = searchCoords(iDom)%qualRecv(1,j)
           donors(ii, jj, kk)%frac(1:3) = searchCoords(iDom)%frac(1:3, j)

           donors(ii, jj, kk)%donorBlockId = searchCoords(iDom)%dInd(1, j)
           donors(ii, jj, kk)%ind(1:3)  = searchCoords(iDom)%dInd(2:4, j)
           donors(ii, jj, kk)%gInd(1:8) = searchCoords(iDom)%gInd(1:8, j)

           iBlank(ii, jj, kk) = -1 
           recvStatus(ii, jj, kk) = .True.

        end if
     end do
     
  end do !nn=1, nDom
  deallocate(usedBlocks) ! don't need it anymore

  ! -----------------------------------------------------------------
  ! Caution: 
  ! searchCoords(iDom) (iDom: global block order) is now going to be overwritten
  ! by what I receive from other assigned procs
  ! -----------------------------------------------------------------
  

  !for real (odd), so 1 per each of nblkprocRecv
  loop_invrecvphase1: do iii=1,1*nblkprocRecv 

     ! Probe message to iProc
     call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, sumb_comm_world, status, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
     tag = status(MPI_TAG)
     source = status(MPI_SOURCE)

     if (mod(tag, 2) == 1) then !odd

        tmp = mod(tag, 1000)
        iDom = (tmp-1)/2 + 1 !my current global block order
        nn = iDom - cumDomProc(myid) !my local block order

        if (iDom < StartDom .or. iDom >endDom) then
           print*,'bad receive in myid, iDom, source ',myid, idom, source
        end if

        ! Point to flowDoms arrays
        call setPointers(nn, level, sps)

        ! Get size of rBuffer from iProc
        call MPI_Get_count(status, sumb_real, n, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! We have already allocated rBuffer in initializeSearchCoords
        if (n /= size(searchCoords(iDom)%rBuffer) ) then
           print*, 'Received wrong size rBuffer from what was allocated '
           print*, 'myid, iDom, received size, allocated size: ', &
                    myid, iDom, n, size(searchCoords(iDom)%rBuffer)
           stop
        end if

        searchCoords(iDom)%rBuffer = -one

        ! Now receive the real buffer. Blocking is fine
        call mpi_recv(searchCoords(iDom)%rBuffer, n, sumb_real, source, &
             tag, SUmb_comm_world, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)
   
        ! Copy real buffer. Overwrite existing searchCoords(iDom) arrays
        ! Save frac 
        rsize = 0
        do j=1, searchCoords(iDom)%n
           do i=1, 3
              rSize = rSize + 1
              searchCoords(iDom)%frac(i,j) =  searchCoords(iDom)%rBuffer(rSize)
           end do
        end do
        ! Save qualRecv 
        do j=1, searchCoords(iDom)%n
           rSize = rSize + 1
           searchCoords(iDom)%qualRecv(1,j) = searchCoords(iDom)%rBuffer(rSize)
        end do
     
        ! Compare qualRecv with flowDoms()%qualRecv and save globalBlockId 
        ! of potential searchCoords blocks
        do j=1, searchCoords(iDom)%n

           ! Current receiver cell's indices
           ii = searchCoords(iDom)%fInd(1,j)
           jj = searchCoords(iDom)%fInd(2,j)
           kk = searchCoords(iDom)%fInd(3,j)

           if (searchCoords(iDom)%qualRecv(1,j) < qualRecv(ii, jj, kk) .or. &
               (forceRecv(ii, jj, kk) .and. &
                searchCoords(iDom)%qualRecv(1, j) < large) ) then 

              ! if a the donor quality is better than existing one or
              ! if a donor has been found for the forced receiver already

              qualRecv(ii, jj, kk) = searchCoords(iDom)%qualRecv(1,j)
              donors(ii, jj, kk)%frac(1:3) = searchCoords(iDom)%frac(1:3, j)

              specialArray(nn)%procId(1,j) = source !procId of incoming searchCoords
           end if
        end do

     end if !odd tag

  end do loop_invrecvphase1

  ! Wait for all the send requests for finish. 
  do i=1, sendCount
     call mpi_waitany(sendCount, sendRequests, index, status, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
  enddo

  ! --------------------------------------------------------------
  ! End of 1st phase of send/recvs of real buffer (qualRecv/frac)
  ! --------------------------------------------------------------
  call MPi_barrier(sumb_comm_world, ierr)

  ! ------------------------------------------------------
  ! Begin 2nd phase of sends/recvs of integers (dInd/gInd)
  ! ------------------------------------------------------

  ! ---------------------------------------------------------------------
  ! 2nd phase of sends of searchCoords copies from assigned procs to the
  ! original owners of the searchCoords (or to the compute blocks)
  ! ---------------------------------------------------------------------

  !print*,'Begin 2nd phase,myid: ',myid
  sendCount = 0
  do iDom=1, nDomTotal
     
     if (invsendBlocks(iDom) /= -1) then 
       
        ! I was assigned iDom but I didn't own it
        
        iProc = invsendBlocks(iDom) !send from myid to this iProc

        ! Define buffer 
        ! -------------------------------------------
        n = searchCoords(iDom)%n

        ! Allocate buffer
        iSize = size(searchCoords(iDom)%dInd) + size(searchCoords(iDom)%gInd)
        allocate(searchCoords(iDom)%iBuffer(iSize))

        iSize = 0

        ! Integer buffer
        ! Save donor index
        do j=1, n
           !write(7000+myid,*)iDom,searchCoords(iDom)%dInd(1,j)
           do i=1, 4 !donorBlockId and donor indices
              iSize = iSize + 1
              searchCoords(iDom)%iBuffer(iSize) = searchCoords(iDom)%dInd(i,j)
           end do
        end do

        ! Save global donor index
        do j=1, n
           do i=1, 8
              iSize = iSize + 1
              searchCoords(iDom)%iBuffer(iSize) = searchCoords(iDom)%gInd(i,j)
           end do
        end do
        ! --- end define buffer ---------------------


        ! Integer buffer
        tag = 1000*iProc + 2*(iDom-1) + 2  
        sendCount = sendCount + 1

        call mpi_isend(searchCoords(iDom)%iBuffer, size(searchCoords(iDom)%iBuffer), &
             sumb_integer, iProc, tag, SUmb_comm_world, &
             sendRequests(sendCount), ierr)
        call ECHK(ierr, __FILE__, __LINE__)
        !print*,'myid to tag, iProc -> ',myid, tag, iProc
     end if
  end do !nDomTotal


  ! 2nd phase recvs: now only iBuffer (dInd/gInd)
  !
  ! Use 'specialArray' to sift through procs that I need to receive
  ! dInd/gInd from, corresponding to the identified minimum qualRecv in 1st phase.


  !for integer (even), so 1 per each of nblkprocRecv
  loop_invrecvphase2: do iii=1,1*nblkprocRecv 

     ! Probe message to iProc
     call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, sumb_comm_world, status, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
     tag = status(MPI_TAG)
     source = status(MPI_SOURCE)
     !print*,'myid from tag, source ', myid, tag, source

     if (mod(tag, 2) == 0) then !even

        tmp = mod(tag, 1000)
        iDom = (tmp-2)/2 + 1 !my current global block order
        nn = iDom - cumDomProc(myid) !my local block order

        ! Point to flowDoms arrays
        call setPointers(nn, level, sps)

        ! Get size of rBuffer from iProc
        call MPI_Get_count(status, sumb_integer, n, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! We have already allocated iBuffer in initializeSearchCoords
        if (n /= size(searchCoords(iDom)%iBuffer) ) then
           print*, 'Received wrong size iBuffer from what was allocated '
           print*, 'myid, iDom, received size, allocated size: ',&
                    myid, iDom, n, size(searchCoords(iDom)%iBuffer)
           stop
        end if

        searchCoords(iDom)%iBuffer = -1

        ! Now receive the integer buffer. Blocking is fine
        call mpi_recv(searchCoords(iDom)%iBuffer, n, sumb_integer, source, &
             tag, SUmb_comm_world, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        !print*,'myid from iProc iSize <-',myid, source, n

        ! Copy integer buffer. Overwrite existing searchCoords(iDom) arrays
        ! Save donor index
        iSize = 0
        do j=1, searchCoords(iDom)%n
           do i=1, 4 !donorBlockId and donor indices
              iSize = iSize + 1
              searchCoords(iDom)%dInd(i,j) = searchCoords(iDom)%iBuffer(iSize)
           end do
        end do
        ! Save global donor index
        do j=1, searchCoords(iDom)%n
           do i=1, 8
              iSize = iSize + 1
              searchCoords(iDom)%gInd(i,j) = searchCoords(iDom)%iBuffer(iSize)
           end do
        end do

        ! Use specialArray to see if this globalBlockId (iDom) is where
        ! the potential donor info is to be received. Skip if not.
        do j=1, searchCoords(iDom)%n

           ! Current receiver cell's indices
           ii = searchCoords(iDom)%fInd(1,j)
           jj = searchCoords(iDom)%fInd(2,j)
           kk = searchCoords(iDom)%fInd(3,j)

           if ( (specialArray(nn)%procId(1,j) == source) .and. &
                (iBlank(ii, jj, kk) /= 0) ) then
              ! this source proc had the right donor candidates

              !corr to flowDoms(nn, level, sps)%donors
              donors(ii, jj, kk)%donorBlockID  = searchCoords(iDom)%dInd(1, j)
              donors(ii, jj, kk)%ind(1:3)  = searchCoords(iDom)%dInd(2:4, j)
              donors(ii, jj, kk)%gInd(1:8) = searchCoords(iDom)%gInd(1:8, j)

              iBlank(ii, jj, kk) = -1 
              recvStatus(ii, jj, kk) = .True.

           end if !specialArray
        end do !j=1, searchCoords(iDom)%n
     
     end if !even tag

  end do loop_invrecvphase2

  ! Wait for all the send requests for finish. 
  do i=1, sendCount
     call mpi_waitany(sendCount, sendRequests, index, status, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
  enddo

  ! --------------------------------------------------------------
  ! End of 2nd phase of send/recvs of integer buffer (dInd/gInd)
  ! --------------------------------------------------------------
  !print*,'End 2nd phase,myid: ',myid

  ! This barrier is very important
  call MPi_barrier(sumb_comm_world, ierr)
  !---------------------------------------------------------------------
  call exchangeIBlanks(level, sps, commPatternCell_2nd, internalCell_2nd)
  !---------------------------------------------------------------------

  ! This is the "onera" correction...if a cell is both a receiver
  ! and a donor, force it to be a compute cell and remove it's
  ! receiver status. Do this as long as it isn't a forced
  ! receiver...we can't touch the forced receivers. 
  ! donorstatus was updated earlier from oBlock info after donor search
  
  do nn=1, nDom
     call setPointers(nn, level, sps)

     do k=2, kl
        do j=2, jl
           do i=2, il
             
              if(recvStatus(i, j, k) .and. donorStatus(i, j, k) .and. &
                 .not. forceRecv(i, j, k) ) then

                iBlank(i, j, k) = 1
                recvStatus(i, j, k) = .False.
              end if

           end do
        end do
     end do

  end do !nn=1, nDom


  ! Now call the actual iblank communication
  ! All procs have their iblank info and now the halos would be exchanged
  !---------------------------------------------------------------------
  call MPi_barrier(sumb_comm_world, ierr)
  call exchangeIBlanks(level, sps, commPatternCell_2nd, internalCell_2nd)
  !---------------------------------------------------------------------
  
  !! -----------------------------------------------------------------
  !! Step 15: Reduction of the number of fringes. What we do is look at
  !! all the fringes and see if al the cells in its stencil are also
  !! fringes or holes. If so we can flag that particular cell as a
  !! hole.
  !! -----------------------------------------------------------------

  !do nn=1, nDom
  !   call setPointers(nn, level, sps)

  !   ! Allocate temp variable and initialize to actual values. This
  !   ! is necessary since we cannot update the iblank array as we
  !   ! go, since that will affect the result after it. 

  !   allocate(iBlankTmp(0:ib, 0:jb, 0:kb))
  !   iBlankTmp = iBlank

  !   do k=2, kl
  !      do j=2, jl
  !         do i=2, il

  !            ! See if we can make this fringe a actual hole
  !            if (iBlank(i,j,k) == -1) then 

  !               computeCellFound = .False.

  !               stencilLoop2: do i_stencil=1, N_visc_drdw
  !                  ii = visc_drdw_stencil(i_stencil, 1) + i
  !                  jj = visc_drdw_stencil(i_stencil, 2) + j
  !                  kk = visc_drdw_stencil(i_stencil, 3) + k

  !                  if (iBlank(ii, jj, kk) == 1) then 
  !                     computeCellFound = .True.
  !                  end if
  !               end do stencilLoop2

  !               if (.not. computeCellFound) then 
  !                  ! Everything was interpoolated so hard blank to zero
  !                  iBlankTmp(i, j, k) = 0
  !               end if
  !            end if
  !         end do
  !      end do
  !   end do

  !   ! Finally copy back and deallocate
  !   iBlank = iBlankTmp
  !   deallocate(iBlanktmp)
  !end do
  !
  !! Exchange iblanks again
  !call exchangeIBlanks(level, sps, commPatternCell_2nd, internalCell_2nd)
  !
  !call MPi_barrier(sumb_comm_world, ierr)
  !call writePartionedMesh('test.dat')
  !call mpi_barrier(sumb_comm_world, ierr)

  ! -----------------------------------------------------------------
  ! Step 16: The algorithm is now complete. Run the checkOverset
  ! algorithm to verify that we actually have a valid interpolation
  ! -----------------------------------------------------------------
  call checkOverset

  ! -----------------------------------------------------------------
  ! Step 17: We can now create the required comm structures based on
  ! our interpolation.
  ! -----------------------------------------------------------------
  call initializeOversetComm
  timeB = mpi_wtime()
  print *,'myid, total time:', myid, timeB-timeA

  ! Deallocate local arrays
  ! ---------------------------------
  deallocate(blkProc, blkLocalID)
  deallocate(Xmin, Xmax, minVol)
  deallocate(nnzProc, cumNNZProc)

  deallocate(overlap%data, overlap%colInd, overlap%rowPtr, &
             overlap%assignedProc)

  deallocate(sendBlocks, recvBlocks)
  deallocate(invsendBlocks, invrecvBlocks)

  ! Deallocate other arrays
  call deallocateInterpolationArrays
  ! --- end deallocate local arrays ---

end subroutine computeOversetInterpolationParallel


! ----------------------------------------------------------------
! Other subroutines
! ----------------------------------------------------------------

subroutine initializeOBlock(oBlock)

  ! This routine allocates the data for the supplied oBlock using the
  !  data currently in blockPointers
  use constants
  use overset
  use blockPointers
  use adtAPI
  use BCTypes
  use cgnsGrid
  implicit none 

  ! Input Params
  type(oversetBlock), intent(inout) :: oBlock

  ! Working paramters
  integer(kind=intType) :: i, j, k, mm, nADT, nHexa, planeOffset

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

  allocate( &
       oBlock%qualDonor(1, ie*je*ke), &
       oBlock%globalCell(0:ib, 0:jb, 0:kb), &
       oBlock%iblank(0:ib, 0:jb, 0:kb), &
       oBlock%cellStatus(1:ie, 1:je, 1:ke), & 
       oBlock%donorStatus(1:ie, 1:je, 1:ke), & 
       oBlock%recvStatus(1:ie, 1:je, 1:ke), & 
       oBlock%forceRecv(1:ie, 1:je, 1:ke), & 
       oBlock%invalidDonor(1:ie, 1:je, 1:ke))

  ! Initialize invalidDonor to 0 (false)
  oBlock%invalidDonor = 0

  ! Initialize cellStatus to 1 (iblank=1 or compute point)
  oBlock%cellStatus = 1

  ! Initialize status values
  oBlock%donorStatus = .False.
  oBlock%recvStatus = .False.
  oBlock%forceRecv = .False.
  

  ! Define cellstatus:
  ! We need to flag the status of cells according to the
  ! oversetOuterBoundary condition. cellStatus is by default set
  ! to 1

  do mm=1,nBocos
     if(BCType(mm) == OversetOuterBound) then 
        select case (BCFaceID(mm))
        case (iMin)
           oBlock%cellStatus(1:3, :, :) = -1
        case (iMax)
           oBlock%cellStatus(nx:ie, :, :) = -1
        case (jMin)
           oBlock%cellStatus(:, 1:3, :) = -1
        case (jMax)
           oBlock%cellStatus(:, ny:je, :) = -1
        case (kMin)
           oBlock%cellStatus(:, :, 1:3) = -1
        case (kMax)
           oBlock%cellStatus(:, :, nz:ke) = -1
        end select
     end if
  end do

  ! Combine the cellStatus information array with our current
  ! iBlank data we got from the holesInsideBody calc.
  do k=1,ke
     do j=1,je
        do i=1,ie
           if (iblank(i, j, k) == 0) then 
              oBlock%cellStatus(i, j, k) = iBlank(i, j, k)
           end if
        end do
     end do
  end do


  ! Copy Volume to qualDonor and do minVol while we're at it
  oBlock%minVolume = Large
  mm = 0
  do k=1,ke
     do j=1,je
        do i=1,ie
           mm = mm + 1
           oBlock%qualDonor(1, mm) = vol(i, j, k)
           oBlock%minVolume = min(oBlock%minVolume, vol(i, j, k))
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

end subroutine initializeOBlock

subroutine packOBlock(oBlock)

  use overset
  use constants
  implicit none
  
  ! Pack up everything we need for this block into its own buffer
  ! inlucding the data required for hte ADTree

  ! Input/Output Parameters
  type(oversetBlock), intent(inout) :: oBlock

  ! Working paramters
  integer(kind=intType) :: rSize, iSize, i, j, k, nHexa, nADT

  iSize = 0
  rSize = 0

  ! Count up the integers we want to send:

  iSize = iSize + 12 ! All block indices
  
  iSize = iSize + size(oBlock%hexaConn)

  iSize = iSize + size(oBlock%invalidDonor)

  iSize = iSize + size(oBlock%cellStatus)

  iSize = iSize + size(oBlock%globalCell)

  iSize = iSize + oBlock%ADT%nLeaves*2 ! The two itegers for the
                                       ! childer in each leaf

  ! Count up the reals we ned to send:
  rSize = rSize + size(oBlock%qualDonor)
  
  rSize = rSize + size(oBlock%xADT)
  
  rSize = rSize + oBlock%ADT%nBBoxes*6 ! Cell bounding boxes

  ! Bounding boxes for leaves
  rSize = rSize + oBlock%ADT%nLeaves*12

  rSize = rSize + 1 ! Min block volume
  
  ! Allocate the buffers
  allocate(oBlock%rBuffer(rSize), oBlock%iBuffer(iSize))

  ! Reset the integer counter and add all the integers on this pass
  iSize = 0
  
  oBlock%iBuffer(1) = oBlock%ib
  oBlock%iBuffer(2) = oBlock%jb
  oBlock%iBuffer(3) = oBlock%kb
  oBlock%iBuffer(4) = oBlock%ie
  oBlock%iBuffer(5) = oBlock%je
  oBlock%iBuffer(6) = oBlock%ke
  oBlock%iBuffer(7) = oBlock%il
  oBlock%iBuffer(8) = oBlock%jl
  oBlock%iBuffer(9) = oBlock%kl
  oBlock%iBuffer(10)= oBlock%nx
  oBlock%iBuffer(11) = oBlock%ny
  oBlock%iBuffer(12) = oBlock%nz

  iSize = iSize + 12

  nHexa = oBlock%il * oBlock%jl * oBlock%kl
  nADT = oBlock%ie * oBlock%je * oBlock%ke

  do j=1, nHexa
     do i=1, 8
        iSize = iSize + 1
        oBlock%iBuffer(iSize) = oBlock%hexaConn(i, j)
     end do
  end do

  do k=1,oBlock% ke
     do j=1, oBlock%je
        do i=1, oBlock%ie
           iSize = iSize + 1
           oBlock%iBuffer(iSize) = oBlock%invalidDonor(i, j, k)
        end do
     end do
  end do

  ! Add cellstatus info into buffer
  do k=1,oBlock%ke
     do j=1, oBlock%je
        do i=1, oBlock%ie
           iSize = iSize + 1
           oBlock%iBuffer(iSize) = oBlock%cellStatus(i, j, k)
        end do
     end do
  end do

  do k=0, oBlock%kb
     do j=0, oBlock%jb
        do i=0, oBlock%ib
           iSize = iSize + 1
           oBlock%iBuffer(iSize) = oBlock%globalCell(i, j, k)
        end do
     end do
  end do
  
  do i=1, oBlock%ADT%nLeaves
     iSize = iSize + 1
     oBlock%iBuffer(iSize) = oBlock%ADT%ADTree(i)%children(1)
     iSize = iSize + 1
     oBlock%iBuffer(iSize) = oBlock%ADT%ADTree(i)%children(2)  
  end do

  ! Reset the real counter and add all the real values on this pass.
  rSize = 0
  
  do i=1, oBlock%ie * oBlock%je * oBlock%ke
     rSize = rSize + 1
     oBlock%rBuffer(rSize) = oBlock%qualDonor(1, i)
  end do

  do j=1, oBlock%ie * oBlock%je * oBlock%ke
     do i=1, 3
        rSize = rSize + 1
        oBlock%rBuffer(rSize) = oBlock%xADT(i, j)
     end do
  end do

  do i=1, oBlock%ADT%nBboxes
     oBlock%rBuffer(rSize+1:rSize+6) = oBlock%ADT%xBBox(:, i)
     rSize = rSize + 6
  end do
  
  do i=1, oBlock%ADT%nLeaves
     oBlock%rBuffer(rSize+1:rSize+6) = oBlock%ADT%ADTree(i)%xMin(:)
     rSize = rSize + 6

     oBlock%rBuffer(rSize+1:rSize+6) = oBlock%ADT%ADTree(i)%xMax(:)
     rSize = rSize + 6
  end do

  rSize = rSize + 1
  oBlock%rBuffer(rSize) = oBlock%minVolume

end subroutine packOBlock

subroutine unpackOBlock(oBlock)

  use adtData
  use overset
  implicit none

  ! unPack up everything we need for this block from its own buffer
  ! and reconstitute the data required for the ADTree. It is assumed
  ! the buffers are already allocated and the data is available. This
  ! does the exact *OPPOSITE* operation as the packBlock() routine

  ! Input/Output Parameters
  type(oversetBlock), intent(inout) :: oBlock

  ! Working paramters
  integer(kind=intType) :: rSize, iSize, i, j, k, nHexa, nADT

  ! Reset the integer counter and add all the integers on this pass
  iSize = 0
  
  oBlock%ib = oBlock%iBuffer(1) !ie+1?
  oBlock%jb = oBlock%iBuffer(2)
  oBlock%kb = oBlock%iBuffer(3)
  oBlock%ie = oBlock%iBuffer(4)
  oBlock%je = oBlock%iBuffer(5)
  oBlock%ke = oBlock%iBuffer(6)
  oBlock%il = oBlock%iBuffer(7)
  oBlock%jl = oBlock%iBuffer(8)
  oBlock%kl = oBlock%iBuffer(9)
  oBlock%nx = oBlock%iBuffer(10)
  oBlock%ny = oBlock%iBuffer(11)
  oBlock%nz = oBlock%iBuffer(12)

  iSize = iSize + 12

  nHexa = oBlock%il * oBlock%jl * oBlock%kl
  nADT = oBlock%ie * oBlock%je * oBlock%ke

  ! Allocate the remainder of the arrays in oBlock.
  allocate(oBlock%hexaConn(8, nHexa))
  allocate(oBlock%invalidDonor(1:oBlock%ie, 1:oBlock%je, 1:oBlock%ke))
  allocate(oBlock%cellStatus(1:oBlock%ie, 1:oBlock%je, 1:oBlock%ke))
  allocate(oBlock%iblank(0:oBlock%ib, 0:oBlock%jb, 0:oBlock%kb)) !check size, ib=ie+1?
  allocate(oBlock%donorStatus(1:oBlock%ie, 1:oBlock%je, 1:oBlock%ke))
  allocate(oBlock%recvStatus(1:oBlock%ie, 1:oBlock%je, 1:oBlock%ke))
  allocate(oBlock%forceRecv(1:oBlock%ie, 1:oBlock%je, 1:oBlock%ke))
  allocate(oBlock%globalCell(0:oBlock%ib, 0:oBlock%jb, 0:oBlock%kb))
  allocate(oBlock%qualDonor(1, oBlock%ie * oBlock%je * oBlock%ke))
  allocate(oBlock%xADT(3, nADT))

  ! Initialize arrays and status values
  oBlock%iblank = 1
  oBlock%donorStatus = .False.
  oBlock%recvStatus = .False.
  oBlock%forceRecv = .False.

  ! -------------------------------------------------------------------
  ! Once we know the sizes, allocate all the arrays in the
  ! ADTree. Since we are not going to call the *actual* build routine
  ! for the ADT, we need to set all the information ourselves. This
  ! essentially does the same thing as buildSerialHex.
  oBlock%ADT%adtType = adtVolumeADT
  oBlock%ADT%nNodes = nADT
  oBlock%ADT%nTetra = 0
  oBlock%ADT%nPyra = 0
  oBlock%ADT%nPrisms = 0
  oBlock%ADT%nTria = 0
  oBlock%ADT%nQuads = 0
  oBlock%ADT%coor => oBlock%xADT
  oBlock%ADT%hexaConn => oBlock%hexaConn
  nullify(oBlock%ADT%tetraConn, oBlock%ADT%pyraConn, oBlock%ADT%prismsConn)
  oBlock%ADT%nBBoxes = nHexa
  allocate(oBlock%ADT%xBBOX(6, nHexa))
  allocate(oBlock%ADT%elementType(nHexa))
  allocate(oBlock%ADT%elementID(nHexa))

  ! All hexas
  oBlock%ADT%elementType = adtHexahedron

  do i=1,nHexa
     oBlock%ADT%elementID(i) = i
  end do

  oBlock%ADT%nLeaves = oBlock%ADT%nBBoxes - 1
  if(oBlock%ADT%nBBoxes <= 1) oBlock%ADT%nLeaves = oBlock%ADT%nLeaves + 1
  allocate(oBlock%ADT%ADTree(oBlock%ADT%nLeaves))

  ! -------------------------------------------------------------------

  ! Now continue copying out the integer values
  do i=1, nHexa
     do j=1, 8
        iSize = iSize + 1
        oBlock%hexaConn(j, i) = oBlock%iBuffer(iSize)
     end do
  end do

  do k=1, oBlock%ke
     do j=1, oBlock%je
        do i=1, oBlock%ie
           iSize = iSize + 1
           oBlock%invalidDonor(i, j, k) = oBlock%iBuffer(iSize)
        end do
     end do
  end do

  ! Copy cellStatus integer now
  do k=1, oBlock%ke
     do j=1, oBlock%je
        do i=1, oBlock%ie
           iSize = iSize + 1
           oBlock%cellStatus(i, j, k) = oBlock%iBuffer(iSize)
        end do
     end do
  end do

  do k=0, oBlock%kb
     do j=0, oBlock%jb
        do i=0, oBlock%ib
           iSize = iSize + 1
           oBlock%globalCell(i, j, k) = oBlock%iBuffer(iSize)
        end do
     end do
  end do
   
  do i=1, oBlock%ADT%nLeaves
     iSize = iSize + 1
     oBlock%ADT%ADTree(i)%children(1) = oBlock%iBuffer(iSize)
     iSize = iSize + 1
     oBlock%ADT%ADTree(i)%children(2) = oBlock%iBuffer(iSize)
  end do
  
  ! Now copy out the real values
  rSize = 0
    
  do i=1, oBlock%ie * oBlock%je * oBlock%ke
     rSize = rSize + 1
     oBlock%qualDonor(1, i) =  oBlock%rBuffer(rSize)
  end do

  do j=1, oBlock%ie * oBlock%je * oBlock%ke
       do i=1, 3
        rSize = rSize + 1
        oBlock%xADT(i, j) = oBlock%rBuffer(rSize)
     end do
  end do

  do i=1, oBlock%ADT%nBboxes
     oBlock%ADT%xBBox(:, i) = oBlock%rBuffer(rSize+1:rSize+6)
     rSize = rSize + 6
  end do
  
  do i=1, oBlock%ADT%nLeaves
      oBlock%ADT%ADTree(i)%xMin(:) = oBlock%rBuffer(rSize+1:rSize+6)
     rSize = rSize + 6

     oBlock%ADT%ADTree(i)%xMax(:) = oBlock%rBuffer(rSize+1:rSize+6)
     rSize = rSize + 6
  end do

  rSize = rSize + 1
  oBlock%minVolume = oBlock%rBuffer(rSize)

end subroutine unpackOBlock

subroutine initializeSearchCoords(sBlock)

  ! This routine allocates and initializes the requested searchCoord
  ! structure. We add 'x' using the x currently set in blockPointers
  use blockPointers
  use overset
  implicit none

  ! Input Params
  type(oversetSearchCoords), intent(inout) :: sBlock

  ! Working parameters
  integer(kind=intType) :: i, j, k, mm, n, iDim, iSize, rSize

  !nx = il-1 or 2:il or 1:ie

  n = nx*ny*nz
  allocate(sBlock%x(3, n), sBlock%fInd(3, n), sBlock%dInd(4, n), &
       sBlock%gInd(8, n), sBlock%frac(3,n))
  allocate(sBlock%qualRecv(1,n))

  sBlock%dInd = -1
  sBlock%gInd = -1
  sBlock%frac = -one
  sBlock%n = n
  mm = 0
  do k=2, kl
     do j=2, jl
        do i=2, il
           mm = mm + 1
           do iDim=1, 3
              sBlock%x(iDim, mm) = eighth*(&
                   x(i-1, j-1, k-1, iDim) + &
                   x(i  , j-1, k-1, iDim) + &
                   x(i-1, j  , k-1, iDim) + &
                   x(i  , j  , k-1, iDim) + &
                   x(i-1, j-1, k  , iDim) + &
                   x(i  , j-1, k  , iDim) + &
                   x(i-1, j  , k  , iDim) + &
                   x(i  , j  , k  , iDim))
              sBlock%fInd(:, mm) = (/i, j, k/)
           end do
           sBlock%qualRecv(1, mm) = vol(i, j, k) !check
   
           ! if forceRecv before sBlock initialization make its qualRecv large
           if (forceRecv(i, j, k)) sBlock%qualRecv(1, mm) = large !1.0e20
        end do
     end do
  end do
  !sBlock%qualRecv = large !1.e20 ! assign large value to be eligible to receive

  !initialize buffer sizes 
  iSize = 0
  rSize = 0
  
  ! Add integer buffer size
  iSize = iSize + size(sBlock%dInd) ! dInd
  iSize = iSize + size(sBlock%gInd) ! gInd
 
  ! Add real buffer size
  rSize = rSize + size(sBlock%frac) ! frac
  rSize = rSize + size(sBlock%qualRecv) ! qualRecv

  allocate(sBlock%iBuffer(iSize), sBlock%rBuffer(rSize))
  
  
end subroutine initializeSearchCoords

subroutine newSearch(oBlock, sCoords, globalBlockId, iDomRecv)

  use constants
  use overset
  use inputOverset
  use adtLocalSearch
  implicit none

  type(oversetBlock), intent(inout) :: oBlock
  type(oversetSearchCoords), intent(inout) :: sCoords
  integer(kind=intType) :: globalBlockId, iDomRecv

  ! Working Varaibles
  integer(kind=intType) :: i, j, k, ii, jj, kk, iii, jjj, kkk, l, mm, mmm, tmp
  integer(kind=intType) :: ir, jr, kr
  integer(kind=intType) :: nCoor, nHexa, nInterpol, elemID, nalloc
  logical :: invalidDonors
  real(kind=realType) :: uvw(4)
  integer(Kind=intType) ::intInfo(3), jjADT, localSearches
  real(kind=realType) :: donorQual

  ! Variables we have to pass the ADT search routine
  integer(kind=intType), dimension(:), pointer :: BB
  integer(kind=intType), dimension(:), pointer :: frontLeaves
  integer(kind=intType), dimension(:), pointer :: frontLeavesNew

  nInterpol = 1 ! we get the ADT to compute the interpolated volume for us. 

  ! Allocate the (pointer) memory that may be resized as necessary for
  ! the singlePoint search routine. 
  allocate(BB(20), frontLeaves(25), frontLeavesNew(25))

  ! Search the cells one at a time:
  do i=1, sCoords%n

     ! *** To add later ***
     ! Short cut: Can we skip search, if the min cell volume of 
     ! donor (oBlock) is larger than current cell volume (sCoords), 
     ! But, this requires the knowledge of sCoords cell forcedReceiver info
     ! which is expensive to keep now. Ignored for now...

     call containmentTreeSearchSinglePoint(oBlock%ADT, &
          sCoords%x(:, i), intInfo, uvw, oBlock%qualDonor, &
          nInterpol, BB, frontLeaves, frontLeavesNew)

     elemFound: if (.not. intInfo(1) < 0) then 

        donorQual = uvw(4)

        ! Two conditions are check: IF the donorQual is
        ! better than my quality, we accept the donor. If we
        ! have a forced receiver that does not yet have donor
        ! assigned to it, we also accept it, since it may be
        ! only one we get. 

        ! Unwind the element indices for the donor
        ! Remember we have (il, jl, kl) elements in the dual
        ! mesh. 

        elemID = intInfo(3) 
        tmp = elemID - 1
        ii = mod(tmp, oBlock%il) + 1
        jj = mod(tmp/oBlock%il, oBlock%jl) + 1
        kk = tmp/(oBlock%il*oBlock%jl) + 1

        ! We need to check if any of the 8 dual nodes to make
        ! sure that they are valid donors

        invalidDonors = .False.
        do iii=ii, ii+1
           do jjj=jj, jj+1
              do kkk=kk, kk+1
                 if (oBlock%invalidDonor(iii, jjj, kkk) == 1) then 
                    invalidDonors = .True.
                 end if
              end do
           end do
        end do

        if ( (.not. invalidDonors) .and. &
             (donorQual < sCoords%qualRecv(1,i)) ) then 

           ! Save the necessary all donor information about
           ! the donor on the receiving processor (an on-proc
           ! block
           sCoords%dInd(1, i)   = globalBlockId
           sCoords%dInd(2:4, i) = (/ii, jj, kk/)
           sCoords%frac(:, i) = uvw(1:3)

           ! Save the global indices as well. 
           sCoords%gInd(1, i) = oBlock%globalCell(ii  , jj  , kk  )
           sCoords%gInd(2, i) = oBlock%globalCell(ii+1, jj  , kk  )
           sCoords%gInd(3, i) = oBlock%globalCell(ii  , jj+1, kk  )
           sCoords%gInd(4, i) = oBlock%globalCell(ii+1, jj+1, kk  )
           sCoords%gInd(5, i) = oBlock%globalCell(ii  , jj  , kk+1)
           sCoords%gInd(6, i) = oBlock%globalCell(ii+1, jj  , kk+1)
           sCoords%gInd(7, i) = oBlock%globalCell(ii  , jj+1, kk+1)
           sCoords%gInd(8, i) = oBlock%globalCell(ii+1, jj+1, kk+1)

           ! Set my receiver quality from this donor to compare
           ! against any potential new donor
           sCoords%qualRecv(1,i) = donorQual

           ! Also update the donor's donorstatus while at it
           ! Which is later used in 'onera' corrections
           do iii=ii, ii+1
              do jjj=jj, jj+1
                 do kkk=kk, kk+1
                    oBlock%donorStatus(iii, jjj, kkk) = .TRUE.
                 end do
              end do
           end do
           
        end if
     end if elemFound
  end do
end subroutine newSearch


subroutine allocateOverset(nn, level, sps)

  use blockPointers
  implicit none

  ! Allocate (if necessary) and initialize the on-block overset
  ! information for the requested block. 

  ! Input Parameters
  integer(kind=intType) :: nn, level, sps

  ! Working parameters
  integer(kind=intType) :: i, j, k, iDim

  call setPointers(nn, level, sps)

  if (.not. associated(flowDoms(nn, level, sps)%donors)) then 
     allocate(flowDoms(nn, level, sps)%donors(0:ib, 0:jb, 0:kb))
  end if

  if (.not. associated(flowDoms(nn, level, sps)%forceRecv)) then 
     allocate(flowDoms(nn, level, sps)%forceRecv(2:il, 2:jl, 2:kl))
  end if

  if (.not. associated(flowDoms(nn, level, sps)%recvStatus)) then 
     allocate(flowDoms(nn, level, sps)%recvStatus(2:il, 2:jl, 2:kl))
  end if

  if (.not. associated(flowDoms(nn, level, sps)%donorStatus)) then 
     allocate(flowDoms(nn, level, sps)%donorStatus(1:ie, 1:je, 1:ke))
  end if

  if (.not. associated(flowDoms(nn, level, sps)%xSearch)) then 
     allocate(flowDoms(nn, level, sps)%xSearch(3, 2:il, 2:jl, 2:kl))
  end if

  if (.not. associated(flowDoms(nn, level, sps)%qualRecv)) then 
     allocate(flowDoms(nn, level, sps)%qualRecv(2:il, 2:jl, 2:kl))
  end if

  ! While we're here, compute xSearch which is just the cell centers. 
  do k=2, kl
     do j=2, jl
        do i=2, il
           do iDim=1, 3
              flowDoms(nn, level, sps)%xSearch(iDim, i, j, k) = eighth*(&
                   x(i-1, j-1, k-1, iDim) + &
                   x(i  , j-1, k-1, iDim) + &
                   x(i-1, j  , k-1, iDim) + &
                   x(i  , j  , k-1, iDim) + &
                   x(i-1, j-1, k  , iDim) + &
                   x(i  , j-1, k  , iDim) + &
                   x(i-1, j  , k  , iDim) + &
                   x(i  , j  , k  , iDim))
           end do

           ! And set qualRecv which is just our own volume
           flowDoms(nn, level, sps)%qualRecv(i, j, k) = vol(i, j, k)
        end do
     end do
  end do

  ! Initialize all the data in the donor derived type
  do k=2, kl
     do j=2, jl
        do i=2, il
           flowDoms(nn, level, sps)%donors(i, j, k)%donorProcID = -1
           flowDoms(nn, level, sps)%donors(i, j, k)%donorBlockID = -1
           flowDoms(nn, level, sps)%donors(i, j, k)%frac = zero
           flowDoms(nn, level, sps)%donors(i, j, k)%ind = -1
           flowDoms(nn, level, sps)%donors(i, j, k)%gInd = -1
        end do
     end do
  end do

  ! These default to false. Will be modified as necessary before
  ! searching.
  flowDoms(nn, level, sps)%recvStatus = .False.
  flowDOms(nn, level, sps)%forceRecv = .False. 
  flowDOms(nn, level, sps)%donorStatus = .False. 

end subroutine allocateOverset


subroutine deallocateInterpolationArrays

  use constants
  use communication
  use blockPointers
  use overset
  implicit none

  ! Local Variables

  ! Global block related
  if (allocated(nDomProc)) deallocate(nDomProc)
  if (allocated(cumDomProc)) deallocate(cumDomProc)
  if (allocated(dims)) deallocate(dims)

  ! ADT tree blocks
  if (allocated(oBlocks)) deallocate(oBlocks)

  ! Search coordinates blocks
  if (allocated(searchCoords)) deallocate(searchCoords)

  ! Special array for receiving searchCoords
  if (allocated(specialArray)) deallocate(specialArray)

end subroutine deallocateInterpolationArrays


! ! ============= DEBUG CODE FOR CHECKING SPARSE MATRIX ASSEMBLY ==================
! ! Debug code to check that we did that correctly...
!  allocate(globalOverlap(nDomTotal, nDomTotal))
!  globalOverlap = .False.
!  do nn=1,nDom
!     iDom = cumDomProc(myid) + nn
!     do jDom=1, nDomtotal
!        globalOverlap(iDom, jDom) = localOverlap(nn, jDom)
!     end do
!  end do

! call mpi_allreduce(MPI_IN_PLACE, globalOverlap, nDomTotal**2, &
!      mpi_logical, MPI_LOR, sumb_comm_world, ierr)


! if (myid == 0) then 
!    do i=1, nDomTotal
!       do j=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
!          if (globalOverlap(i, overlap%colInd(j))) then 
!             ! We're good
!          else
!             print *,'Error in sparse matrix: (i,j) =', i, overlap%colInd(j)
!          end if
!       end do
!    end do
! end if
! ! ================================================================================


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


