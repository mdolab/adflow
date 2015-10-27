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
  implicit none

  ! Local Variables
  integer(kind=intType) ::i, j, k,  level, sps, ii, iDom, jDom, iDim, nn, mm, tmp
  integer(kind=intType) ::  jj, kk, ip2, im2, jp2, jm2, kp2, km2, ind(3), n
  integer(kind=intType) :: donorBlockID, i_stencil, rSize, iSize, jj1
  integer(kind=intType) :: startDom, endDom, nblkRecv, source
  real(kind=realType), dimension(:, :), allocatable :: xMin, xMax
  real(kind=realType), dimension(:), allocatable :: minVol
  logical, dimension(:, :), allocatable :: localOverlap
  logical, dimension(:), allocatable :: receiveBlocks, ownReceiveBlocks
  logical, dimension(:, :), allocatable :: sendBlocks
  integer(kind=intType), dimension(:), allocatable :: recvBlocks

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
  ! also generate cumDomsProc which is cumulative form of the number
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


  ! Debug output for localOverlap (still kinda square-ish) 
  if (myid ==  0) then 
     print *, '------------------------------' 
  end if
  do iProc=0, nProc-1
     if (iProc == myid) then 
        do nn=1,nDom 
           call setPointers(nn, 1,1) 
           print *, nn+cumDomProc(myid), ' ', localOverlap(nn, :) 
        end do
        print *, '------------------------------'
     end if
     call MPI_barrier(sumb_comm_world, ierr)
  end do

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

  timeB = mpi_wtime()
  print *,'myid time 1:', myid, timeB-timeA

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
  
  if (myid == 0) then 
     ! Now dump out who owns what:
     do i=1,nDomTotal
        write(*, "(a,I4, a)", advance='no'), 'Row:', i, "   "
        do jj=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
           write(*, "(a,I2, a, I6)", advance='no'), "(", overlap%colInd(jj), ")", int(overlap%data(jj))
        end do
        write(*, *) " "
     end do
     
     print *, '--------------------------------------'
     ! Now dump out who owns what:
     do i=1,nDomTotal
        write(*, "(a,I4, a)", advance='no'), 'Row:', i, "   "
        do jj=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
           write(*, "(a,I2, a, I6)", advance='no'), "(", overlap%colInd(jj), ")", int(overlap%assignedProc(jj))
        end do
        write(*, *) " "
     end do
  end if

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
              recvBlocks(iDom) = blkProc(jDom)
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
  timeB = mpi_wtime()
  print *,'myid before comm:', myid, timeB-timeA

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
           tag = iDom
           sendCount = sendCount + 1

           call mpi_isend(searchCoords(iDom)%x, searchCoords(iDom)%n*3, &
                sumb_real, iProc, tag, SUmb_comm_world, &
                sendRequests(sendCount), ierr)
           call ECHK(ierr, __FILE__, __LINE__)
        end if
     end do
  end do

 ! Now post all the receives...these block, but we use MPI_ANY_SOURCE
  ! and MPI_ANY_TAG
  do i=1, nBlkRecv
     
     ! Probe the message
     call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, sumb_comm_world, status, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
     tag = status(MPI_TAG)
     iDom = tag
     source = status(MPI_SOURCE)

     ! Get size
     call MPI_Get_count(status, sumb_real, n, ierr)
     call ECHK(ierr, __FILE__, __LINE__)

     ! Allocate space for the receive
     allocate(&
          searchCoords(iDom)%x(3, n/3), &
          searchCoords(iDom)%dInd(3, n/3), &
          searchCoords(iDom)%gInd(8, n/3), &
          searchCoords(iDom)%frac(3, n/3))
     searchCoords(iDom)%n = n/3
     searchCoords(iDom)%arrSize = n/3

     ! Now actually receive the real buffer. Blocking is fine
     call mpi_recv(searchCoords(iDom)%x, n, sumb_real, source, &
          tag, SUmb_comm_world, status, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
  end do

  ! Wait for all the send requests for finish. 
  do i=1, sendCount
     call mpi_waitany(sendCount, sendRequests, index, status, ierr)
     call ECHK(ierr, __FILE__, __LINE__)
  enddo
  
  ! Finally we can do the searches:
  overlap%data = zero
  do iDom=1, nDomTotal
     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
        jDom = overlap%colInd(jj)
        if (overlap%assignedProc(jj) == myid) then 
           time1 = mpi_wtime()
           call newSearch(oBlocks(iDom), searchCoords(jDom))
           time2 = mpi_wtime()
           overlap%data(jj) = time2-time1
        end if
     end do
  end do
  timeB = mpi_wtime()
  print *,'myid time 3:', myid, timeB-timeA
  call mpi_barrier(sumb_comm_world, ierr)


  ! Determine the total costs for everyone. 
  call mpi_allreduce(MPI_IN_PLACE, overlap%data, overlap%nnz, sumb_real, MPI_SUM, &
       sumB_comm_world, ierr)

  call oversetLoadBalance(overlap)

  if (myid == 0) then 
     ! Now dump out who owns what:
     do i=1,nDomTotal
        write(*, "(a,I4, a)", advance='no'), 'Row:', i, "   "
        do jj=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
           write(*, "(a,I2, a, F6.3)", advance='no'), "(", overlap%colInd(jj), ")", overlap%data(jj)
        end do
        write(*, *) " "
     end do
     
     print *, '--------------------------------------'
     ! Now dump out who owns what:
     do i=1,nDomTotal
        write(*, "(a,I4, a)", advance='no'), 'Row:', i, "   "
        do jj=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
           write(*, "(a,I2, a, I6)", advance='no'), "(", overlap%colInd(jj), ")", int(overlap%assignedProc(jj))
        end do
        write(*, *) " "
     end do

     do i=1,overlap%nnz
        print *,overlap%data(i)
     end  do

  end if



  call MPi_barrier(sumb_comm_world, ierr)
  stop

  

  


  
  ! ! In theory we can now *FINALLY* do the search in parallel....We
  ! ! have our own ADTrees along with the search coords from myself and
  ! ! the search coords from the other procs. Essentially we are
  ! ! searching for stuff in our proc's COLUMN

  ! do iDom=1, nDomTotal
  !    do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
  !       jDom = overlap%colInd(jj)

  !       if (jDom >= startDom .and. jDom <= endDom) then 
  !          if (overlap%data(jj) /= 0) then

  !             if (blkProc(iDom) == myid) then 
  !                ! We can use my own search Coords:

  !                call newSearch(oBlocks(jDom), mySearchCoords(jj))

  !             else
  !                ! Otherwise use the ones that we have got from
  !                ! another block
  !                call newSearch(oBlocks(jDom), otherSearchCoords(jj))

  !             end if
  !          end if
  !       end if
  !    end do
  ! end do



  ! ! Now essentially do the reverse of the communication...we send
  ! ! stuff from otherSearchCoords to mySearchCoords.

  ! sendCount = 0

  ! ! Now post all the receives...these block, but we use MPI_ANY_SOURCE
  ! ! and MPI_ANY_TAG
  ! do iDom=1, nDomTotal
  !    if (blkProc(iDom) /= myid) then 

  !       do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
  !          jDom = overlap%colInd(jj)

  !          if (jDom >= startDom .and. jDom <= endDom) then 
  !             if (overlap%data(jj) /= 0) then 
  !                n = otherSearchCoords(jj)%n
  !                call mpi_isend(otherSearchCoords(jj)%frac, 3*n, sumb_real, blkProc(iDom), &
  !                     jj, SUmb_comm_world, sendRequests(sendCount+1), ierr)
  !                call ECHK(ierr, __FILE__, __LINE__)
  !                sendCount = sendCount + 1

  !             end if
  !          end if
  !       end do
  !    end if
  ! end do

  ! do nn=1, nDom 
  !    ! block in gloabl ordering. 
  !    iDom = cumDomProc(myid) + nn  

  !    do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
  !       jDom = overlap%colInd(jj)

  !       if (overlap%data(jj) /= 0 .and. blkProc(jDom) /= myid) then 


  !          call mpi_recv(realRecvBuffer, rSize*3, sumb_real, MPI_ANY_SOURCE, &
  !               MPI_ANY_TAG, SUmb_comm_world, status, ierr)
  !          call ECHK(ierr, __FILE__, __LINE__)
           
  !          jj1 = status(MPI_TAG)
  !          n = mySearchCoords(jj1)%n
  !          do i=1, n
  !             mySearchCoords(jj1)%frac(:, i) = realRecvBuffer(3*i-2:3*i)
  !          end do
  !       end if
  !    end do
  ! end do

  ! ! Wait for all the send requests for finish. Perhaps could use an
  ! ! mpi_waitall here? Probably doesn't really matter.
  ! do i=1, sendCount
  !    call mpi_waitany(sendCount, sendRequests, index, status, ierr)
  !    call ECHK(ierr, __FILE__, __LINE__)
  ! enddo

  














  !  ! Flag cells that are actually donors...this will be much more
  ! ! complex in parallel...here we can just read the index and set
  ! ! directly.
  ! do nn=1, nDomTotal
  !    do k=2, oBlocks(nn)%kl
  !       do j=2, oBlocks(nn)%jl
  !          do i=2, oBlocks(nn)%il
  !             if (oBLocks(nn)%donors(i, j, k)%donorProcID >= 0) then 

  !                ! copy for easier code reading
  !                donorBlockID = oBLocks(nn)%donors(i, j, k)%donorBlockID
  !                ind = oBlocks(nn)%donors(i, j, k)%ind

  !                do kk=0,1
  !                   do jj=0,1
  !                      do ii=0,1

  !                         oBlocks(donorBlockID)%donorStatus(&
  !                              ind(1)+ii, ind(2)+jj, ind(3)+kk) = .True.
  !                      end do
  !                   end do
  !                end do
  !             end if
  !          end do
  !       end do
  !    end do
  ! end do

  ! Update the iBlanks
  call exchangeIBlanks(level, sps, commPatternCell_2nd, internalCell_2nd)

  ! ! -----------------------------------------------------------------
  ! ! Step 14: This is the "onera" correction...if a cell is both a
  ! ! receiver and a donor, force it to be a compute cell and remove
  ! ! it's receiver status. Do this as long as it isn't a forced
  ! ! receiver...we can't touch the forced receivers. 
  ! ! -----------------------------------------------------------------
  ! do nn=1, nDom
  !    call setPointers(nn, level, sps)
  !    do k=2, kl
  !       do j=2, jl
  !          do i=2, il
  !             if (iblank(i,j,k) == -1) then
  !                if ( recvStatus(i, j, k) .and. &
  !                     donorStatus(i, j, k) .and. &
  !                     .not. forceRecv(i, j, k)  ) then

  !                   ! Force it to be compute
  !                   iBlank(i,j,k) = 1
  !                   recvStatus(i, j, k) = .False. 

  !                end if
  !             end if
  !          end do
  !       end do
  !    end do
  ! end do

  ! Exchange iblanks again
  call exchangeIBlanks(level, sps, commPatternCell_2nd, internalCell_2nd)

  ! -----------------------------------------------------------------
  ! Step 15: Reduction of the number of fringes. What we do is look at
  ! all the fringes and see if al the cells in its stencil are also
  ! fringes or holes. If so we can flag that particular cell as a
  ! hole.
  ! -----------------------------------------------------------------

  do nn=1, nDom
     call setPointers(nn, level, sps)

     ! Allocate temp variable and initialize to actual values. This
     ! is necessary since we cannot update the iblank array as we
     ! go, since that will affect the result after it. 

     allocate(iBlankTmp(0:ib, 0:jb, 0:kb))
     iBlankTmp = iBlank

     do k=2, kl
        do j=2, jl
           do i=2, il

              ! See if we can make this fringe a actual hole
              if (iBlank(i,j,k) == -1) then 

                 computeCellFound = .False.

                 stencilLoop2: do i_stencil=1, N_visc_drdw
                    ii = visc_drdw_stencil(i_stencil, 1) + i
                    jj = visc_drdw_stencil(i_stencil, 2) + j
                    kk = visc_drdw_stencil(i_stencil, 3) + k

                    if (iBlank(ii, jj, kk) == 1) then 
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
     iBlank = iBlankTmp
     deallocate(iBlanktmp)
  end do

  ! Exchange iblanks again
  call exchangeIBlanks(level, sps, commPatternCell_2nd, internalCell_2nd)

  ! -----------------------------------------------------------------
  ! Step 16: The algorithm is now complete. Run the checkOverset
  ! algorithm to verify that we actually have a valid interpolation
  ! -----------------------------------------------------------------
  call checkOverset

  ! -----------------------------------------------------------------
  ! Step 17: We can now create the required comm structures based on
  ! our interpolation.
  ! -----------------------------------------------------------------
  !call initializeOversetComm

end subroutine computeOversetInterpolationParallel

subroutine initializeOBlock(oBlock)

  ! This routine allocates the data for the supplied oBlock using the
  !  data currently in blockPointers
  use constants
  use overset
  use blockPointers
  use adtAPI
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
       oBlock%invalidDonor(1:ie, 1:je, 1:ke))

  ! Initialize invalidDonor to 0 (false)
  oBlock%invalidDonor = 0
  
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
  
  oBlock%ib = oBlock%iBuffer(1)
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
  allocate(oBlock%globalCell(0:oBlock%ib, 0:oBlock%jb, 0:oBlock%kb))
  allocate(oBlock%qualDonor(1, oBlock%ie * oBlock%je * oBlock%ke))
  allocate(oBlock%xADT(3, nADT))

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
  integer(kind=intType) :: i, j, k, mm, n, iDim

  n = nx*ny*nz
  allocate(sBlock%x(3, n), sBlock%fInd(3, n), sBlock%dInd(3, n), &
       sBlock%gInd(8, n), sBlock%frac(3,n))

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
        end do
     end do
  end do
end subroutine initializeSearchCoords

subroutine newSearch(oBlock, sCoords)

  use constants
  use overset
  use inputOverset
  use adtLocalSearch
  implicit none

  type(oversetBlock), intent(inout) :: oBlock
  type(oversetSearchCoords), intent(inout) :: sCoords

  ! Working Varaibles
  integer(kind=intType) :: i, j, k, ii, jj, kk, iii, jjj, kkk, l, mm, mmm, tmp
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

        if (.not. invalidDonors) then 

           ! Save the necessary all donor information about
           ! the donor on the receiving processor (an on-proc
           ! block

           sCoords%dInd(:, i) = (/ii, jj, kk/)
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
        end if
     end if elemFound
  end do
end subroutine newSearch


! if (myid == 0) then 

!    ! do i=1, nDomTotal
!    !    print *,'blk seach:', i, blkSearch(i)
!    ! end do
!    ! print *,'Total search:', sum(blkSearch)
!    ! print *, 'data sum:', sum(data)

!    ! ! Processor breakdwon
!    ! do iProc=0, nProc-1
!    !    j =0
!    !    do nn=1, nDomProc(iProc)
!    !       j = j + blkSearch(cumDomProc(iProc) + nn)
!    !    end do
!    !    print *,'Proc Search:', iProc, j
!    ! end do

!    ! Dump the all costs array:
!    open(unit=1, file='allcosts.dat', form='formatted', status='unknown')
!    open(unit=2, file='allcosts.txt', form='formatted', status='unknown')
!    write(1, *) "Variables = I J COST"
!    write(1, *) "Zone I=", nDomTotal, " J=", nDomTotal
!    write(1, *) "DATAPACKING=POINT"
!    i = 0

!    allocate(allCosts(nDOmTotal, nDomTotal))

!    do iDom=1,nDomTotal ! Row Loop
!       do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
!          jDom = overlap%colInd(jj)
!          print *, 'data jj:', jj, data(jj)
!          allCosts(iDom, jDom) = data(jj)
!       end do
!    end do


!    do iDom=1,nDomTotal ! Row Loop
!       do jDom=1,nDomTotal
!          write(1, *), iDom, jDom, allCosts(iDom, jDom), '' 
!          write(2, "(F12.1)", advance='no'), allCosts(iDom, jDom)

!       end do
!       write(2, *), ' '
!    end do


!    close(1)
!    close(2)



! end if


! ! i = i + (ie+1)*(je+1)*(ke+1)*3 + ie*je*ke 
! ! j = j + (ib+1)*(jb+1)*(kb+1)*2           
! iDom = cumDomProc(myid) + nn
! call unpackBlockPart1(realSendBuffer(i+1), dims(:, iDom), iDom)
! call unPackBlockPart2(intSendBuffer(j+1), dims(:, iDom), iDom)
! call buildADTForBlock(iDom)



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




!  sendCount = 0
!   do nn=1, nDom 
!      ! block in gloabl ordering. 
!      iDom = cumDomProc(myid) + nn  

!      ! Size and offset in the buffers
!      rSize = realSendOffset(nn+1) - realSendOffset(nn)
!      i = realSendOffset(nn) + 1

!      iSize = intSendOffset(nn+1) - intSendOffset(nn)
!      j = intSendOffset(nn) + 1

!      do iProc=0, nProc-1
!         if (sendBlocks(iProc, nn)) then 

!            call mpi_isend(realSendBuffer(i), rSize, sumb_real, iProc, &
!                 iDom, SUmb_comm_world, sendRequests(sendCount+1), ierr)
!            call ECHK(ierr, __FILE__, __LINE__)

!            call mpi_isend(intSendBuffer(j), iSize, sumb_integer, iProc,  &
!                 iDom, SUmb_comm_world, sendRequests(sendCount+2), ierr)
!            call ECHK(ierr, __FILE__, __LINE__)

!            sendCount = sendCount + 2

!         end if
!      end do
!   end do

!   ! Now post all the (i)receives
!   recvCount = 0
!   do iDom=1, nDomTotal

!      if (receiveBlocks(iDom)) then 

!         ! Size and offset in the buffers
!         rSize = realRecvOffset(iDom+1) - realRecvOffset(iDOm)
!         i = realRecvOffset(iDom) + 1

!         iSize = intRecvOffset(iDom+1) - intRecvOffset(iDom)
!         j = intRecvOffset(iDom) + 1

!         call mpi_irecv(realRecvBuffer(i), rSize, sumb_real, blkProc(iDom), &
!              iDom, SUmb_comm_world, recvRequests(recvCount+1), ierr)
!         call ECHK(ierr, __FILE__, __LINE__)

!         call mpi_irecv(intRecvBuffer(j), iSize, sumb_integer, blkProc(iDom), &
!              iDom, SUmb_comm_world, recvRequests(totalRecv + recvCount+1), ierr)
!         call ECHK(ierr, __FILE__, __LINE__)


!         ! call mpi_recv(realRecvBuffer(i), rSize, sumb_real, blkProc(iDom), &
!         !      iDom, SUmb_comm_world, ierr)
!         ! call ECHK(ierr, __FILE__, __LINE__)

!         ! call mpi_recv(intRecvBuffer(j), iSize, sumb_integer, blkProc(iDom), &
!         !      iDom, SUmb_comm_world, ierr)
!         ! call ECHK(ierr, __FILE__, __LINE__)

!         recvCount = recvCount + 1

!      end if
!   end do

! ! Now complete the Real receives
!   do iDom=1, nDomTotal

!      if (receiveBlocks(iDom)) then 

!         call mpi_waitany(recvCount, recvRequests(1:totalRecv), index, status, ierr)
!         call ECHK(ierr, __FILE__, __LINE__)

!         ! Extract which block this is by looking at the status:
!         jDom = status(MPI_TAG)
!         i = realRecvOffset(jDom) + 1

!         ! Unpack the real part of the block
!         call unpackBlockPart1(realRecvBuffer(i), dims(:, jDom), jDom)
!      end if
!   end do

!   ! Now complete the Integer receives
!   do iDom=1, nDomTotal

!      if (receiveBlocks(iDom)) then 

!         call mpi_waitany(recvCount, recvRequests(totalRecv+1:2*totalRecv), &
!              index, status, ierr)
!         call ECHK(ierr, __FILE__, __LINE__)

!         ! Extract which block this is by looking at the status:
!         jDom = status(MPI_TAG)
!         j = intRecvOffset(jDom) + 1

!         ! Unpack the integer part of the block here
!         call unpackBlockPart2(intRecvBuffer(j), dims(:, jDom), jDom)

!      end if
!   end do

!   ! Wait for all the send requests for finish. Perhaps could use an
!   ! mpi_waitall here? Probably doesn't really matter.
!   do i=1, sendCount
!      call mpi_waitany(sendCount, sendRequests, index, status, ierr)
!      call ECHK(ierr, __FILE__, __LINE__)
!   enddo

!   i = 0
!   do nn=1,nDomTotal
!      if (allocated(oBlocks(nn)%x)) then 
!         i = i + 1
!      end if
!   end do

!   print *, 'myid allocated:', myid, i

!   call MPi_barrier(sumb_comm_world, ierr)
!   stop



  ! ! Debug output for localOverlap (still kinda square-ish) 
  ! if (myid ==  0) then 
  !    print *, '------------------------------' 
  ! end if
  ! do iProc=0, nProc-1
  !    if (iProc == myid) then 
  !       do nn=1,nDom 
  !          call setPointers(nn, 1,1) 
  !          print *, nn+cumDomProc(myid), ' ', localOverlap(nn, :) 
  !       end do
  !       print *, '------------------------------'
  !    end if
  !    call MPI_barrier(sumb_comm_world, ierr)
  ! end do


  ! ! -----------------------------------------------------------------
  ! ! Step 7: Block Packing: Here we simply pack all of my block
  ! ! information into a single real and integer buffer. At this point
  ! ! we don't know or care who they will be sent to, just do it
  ! ! all. This particular calc is scalable.
  ! ! -----------------------------------------------------------------

  ! ! Determine the size for our send buffers
  ! i = 0
  ! j = 0
  ! realSendOffset(1) = 0
  ! intSendOffset(1) = 0
  ! do nn=1, nDom
  !    call setPointers(nn, level, sps)
  !    i = i + (ie+1)*(je+1)*(ke+1)*3 + ie*je*ke ! real size: coordinates + volume
  !    j = j + (ib+1)*(jb+1)*(kb+1)*2            ! int size: iblank + globalCell
  !    realSendOffset(nn+1) = i
  !    intSendOffset(nn+1) = j
  ! end do

  ! allocate(realSendBuffer(i), intSendBuffer(j))

  ! i = 0
  ! j = 0
  ! do nn=1, nDom
  !    call setPointers(nn, level, sps)
  !    call packBlock(realSendBuffer(i+1), intSendBuffer(j+1))
  !    i = i + (ie+1)*(je+1)*(ke+1)*3 + ie*je*ke 
  !    j = j + (ib+1)*(jb+1)*(kb+1)*2       
  ! end do

  !call writePartionedMesh('test.dat')


  ! else

  !    ! Number of hexa are the real cells.
  !    nHexa = nx * ny * nz
  !    nADT = il * jl * kl

  !    mm = 0
  !    do k=1, kl
  !       do j=1, jl
  !          do i=1, il
  !             mm = mm + 1
  !             oBlock%xADT(:, mm) = x(i, j, k, :)
  !          end do
  !       end do
  !    end do

  !    mm = 0
  !    ! These are the 'elements' of the primal mesh
  !    planeOffset = il*jl
  !    do k=2, kl
  !       do j=2, jl
  !          do i=2, il
  !             mm = mm + 1
  !             oBlock%hexaConn(1, mm) = (k-2)*planeOffset + (j-2)*il + (i-2) + 1
  !             oBlock%hexaConn(2, mm) = oBlock%hexaConn(1, mm) + 1 
  !             oBlock%hexaConn(3, mm) = oBlock%hexaConn(2, mm) + il
  !             oBlock%hexaConn(4, mm) = oBlock%hexaConn(3, mm) - 1 

  !             oBlock%hexaConn(5, mm) = oBlock%hexaConn(1, mm) + planeOffset
  !             oBlock%hexaConn(6, mm) = oBlock%hexaConn(2, mm) + planeOffset
  !             oBlock%hexaConn(7, mm) = oBlock%hexaConn(3, mm) + planeOffset
  !             oBlock%hexaConn(8, mm) = oBlock%hexaConn(4, mm) + planeOffset
  !          end do
  !       end do
  !    end do
  ! end if



! subroutine packBlock(realBuffer, intBuffer)

!   ! Pack the requested block into realBuffer and intBuffer. It is
!   ! assumed that blockPointers are already set. 
!   use BCTypes
!   use blockPointers
!   implicit none

!   ! Input/Output Arguments
!   real(kind=realType), dimension(*), intent(out) :: realBuffer
!   integer(kind=intType), dimension(*), intent(out) :: intBuffer

!   ! Working
!   integer(kind=intType) :: i, j, k, mm, iDim, nRealSend, nIntSend

!   nRealSend = 0

!   ! Note that we take this opportunity to reorder 'x' to have the
!   ! x-y-z values next to each other in memory.
!   do k=0,ke
!      do j=0,je
!         do i=0,ie
!            do iDim=1,3
!               nRealSend = nRealSend + 1
!               realBuffer(nRealSend) = x(i, j, k, iDim)
!            end do
!         end do
!      end do
!   end do

!   ! Pack in the volumes
!   do k=1,ke
!      do j=1,je
!         do i=1,ie
!            nRealSend = nRealSend + 1
!            realBuffer(nRealSend) = vol(i, j, k)
!         end do
!      end do
!   end do

!   nIntSend = 0

!   ! We need to flag the status of cells according to the
!   ! oversetOuterBoundary condition.

!   do mm=1,nBocos
!      if(BCType(mm) == OversetOuterBound) then 
!         select case (BCFaceID(mm))
!         case (iMin)
!            iBlank(1:3, :, :) = -1
!         case (iMax)
!            iBlank(nx:ie, :, :) = -1
!         case (jMin)
!            iBlank(:, 1:3, :) = -1
!         case (jMax)
!            iBlank(:, ny:je, :) = -1
!         case (kMin)
!            iBlank(:, :, 1:3) = -1
!         case (kMax)
!            iBlank(:, :, nz:ke) = -1
!         end select
!      end if
!   end do

!   ! Now pack the cellStatus and globalCell
!   do k=0, kb
!      do j=0, jb
!         do i=0, ib
!            nIntSend = nIntSend + 1
!            intBuffer(nIntSend) = iBlank(i, j, k)

!            nIntSend = nIntSend + 1
!            intBuffer(nIntSend) = globalCell(i, j, k)
!         end do
!      end do
!   end do

! end subroutine packBlock


! subroutine unpackBlockPart1(realBuffer, sizes, blockID)

!   ! This routine allocates and unpacks the data from the
!   ! realBuffer. The remaining operations are in unpackBlockPart2 which
!   ! unpacks the integer buffer. 
!   use constants
!   use overset

!   implicit none 
!   ! Input Params
!   real(kind=realType), dimension(*), intent(in) :: realBuffer
!   integer(kind=intType), intent(in) :: sizes(3), blockID

!   ! Working
!   integer(kind=intType) :: i, j, k, idim
!   integer(kind=intType) :: ii, ie, je, ke, il, jl, kl, ib, jb, kb

!   il = sizes(1)
!   jl = sizes(2)
!   kl = sizes(3)

!   ie = il + 1
!   je = jl + 1
!   ke = kl + 1

!   ib = ie + 1
!   jb = je + 1
!   kb = ke + 1

!   ! Set all the sizes for this block.
!   oBlocks(blockID)%il = il
!   oBlocks(blockID)%jl = jl
!   oBlocks(blockID)%kl = kl

!   oBlocks(blockID)%ie = ie
!   oBlocks(blockID)%je = je
!   oBlocks(blockID)%ke = ke

!   oBlocks(blockID)%ib = ie + 1
!   oBlocks(blockID)%jb = je + 1
!   oBlocks(blockID)%kb = ke + 1

!   oBlocks(blockID)%nx = il - 1
!   oBlocks(blockID)%ny = jl - 1 
!   oBlocks(blockID)%nz = kl - 1

!   allocate( &
!        oBlocks(blockID)%x(3, 0:ie, 0:je, 0:ke), &
!        oBlocks(blockID)%qualDonor(1:ie, 1:je, 1:ke), &
!        oBlocks(blockID)%iblank(0:ib, 0:jb, 0:kb), &
!        oBlocks(blockID)%globalCell(0:ib, 0:jb, 0:kb), &
!        oBlocks(blockID)%invalidDonor(1:ie, 1:je, 1:ke))

!   ! Initialize invalidDonor to False
!   oBlocks(blockID)%invalidDonor = .False. 

!   ! Copy out the primal nodes
!   ii = 0
!   do k=0, ke
!      do j=0, je
!         do i=0, ie
!            do iDim=1,3
!               ii = ii + 1
!               oBlocks(blockID)%x(iDim, i, j, k) = realBuffer(ii)
!            end do
!         end do
!      end do
!   end do

!   ! Copy out the volume. Also compute the min volume for this block
!   ! as we go through the loop since this is basically free.
!   oBlocks(blockID)%minVolume = large
!   do k=1, ke
!      do j=1, je
!         do i=1, ie
!            ii = ii + 1
!            oBlocks(blockID)%qualDonor(i,j,k) = realBuffer(ii)
!            oBlocks(blockID)%minVolume = min(oBlocks(blockID)%minVolume, realBuffer(ii))
!         end do
!      end do
!   end do
! end subroutine unpackBlockPart1

! subroutine unpackBlockPart2(intBuffer, sizes, nn)
!   use overset

!   implicit none

!   ! Input Params
!   integer(kind=intType), dimension(*), intent(in) :: intBuffer
!   integer(kind=intType), intent(in) :: sizes(3), nn

!   ! Working
!   integer(kind=intType) :: i, j, k, idim, ib, jb, kb, ii

!   ib = sizes(1) + 2
!   jb = sizes(2) + 2
!   kb = sizes(3) + 2
!   ! Copy out the iBlank and global cell
!   ii = 0
!   do k=0, kb
!      do j=0, jb
!         do i=0, ib
!            oBlocks(nn)%iBlank(i, j, k) = intBuffer(ii+1)
!            oBlocks(nn)%globalCell(i, j, k) = intBuffer(ii+2)
!            ii = ii + 2
!         end do
!      end do
!   end do

!   ! For this newly received block, flag the holes and fringes as being
!   ! invalid donors. 
!   do k=1, oBlocks(nn)%ke
!      do j=1, oBlocks(nn)%je
!         do i=1, oBlocks(nn)%ie

!            select case (oBlocks(nn)%iBlank(i, j, k))

!            case (0) 
!               ! Blanekd Cell: Set the iblank value and ensure it
!               ! cannot be a donor
!               oBlocks(nn)%invalidDonor(i, j, k) = .True. 

!            case (-1)
!               ! Cell as ben flagged as necessairly being a
!               ! receiver. We therefore force it to be a receiver
!               ! and prevent it from being a donor
!               oBlocks(nn)%invalidDonor(i, j, k) = .True. 
!            end select
!         end do
!      end do
!   end do
! end subroutine unpackBlockPart2



subroutine writePartionedMesh(fileName)

  ! This is a debugging routine for writing out meshes *as they are
  ! partioned*. This can be useful for debugging overset issues.

  use communication
  use blockPointers

  implicit none

  character(len=*), intent(in) :: fileName
  integer(kind=intType) :: nDomTotal, iProc, nn, i, j, k, iDim, iDom, ierr, ii
  integer(kind=intType) :: bufSize, maxSize
  integer(kind=intType), dimension(3, nDom) :: localDim
  integer(kind=intType), dimension(:), allocatable :: nDomProc, cumDomProc
  integer(kind=intType), dimension(:, :), allocatable :: dims
  real(kind=realType), dimension(:), allocatable :: buffer
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
  do i=1,nDomTotal
     maxSize = max(maxSize, dims(1, i)*dims(2,i)*dims(3,i)*3)
  end do

  allocate(buffer(maxSize))

  if (myid == 0) then 
     print *,'writing mesh...'
     ! Root proc does all the writing. Just dump to ascii tecplot
     ! file---really slow.

     open(unit=1, file=fileName, form='formatted', status='unknown')
     write(1, *) "Variables = X Y Z"

     ! Write my own blocks first
     do nn=1,nDom
        call setPointers(nn, 1, 1)
        !write(tmpStr, *) "Proc ", 0, " Local ID", nn
        write(tmpStr, "(a,I2.2,a,I2.2,a)"), """Proc ", 0, " Local ID ", nn , """"
        write(1,*) "ZONE I=", il, " J=",jl, "K=", kl, "T=", trim(tmpStr)
        write(1, *) "DATAPACKING=BLOCK"
        write(1, *) "VARLOCATION=([1,2,3]=NODAL)"

        do iDim=1, 3
           do k=1, kl
              do j=1, jl
                 do i=1, il
                    write(1, *) x(i, j, k, idim)
                 end do
              end do
           end do
        end do
     end do

     ! Now loop over the remaining blocks...receiving each and writing:

     do iProc=1, nProc-1
        do nn=1, nDomProc(iProc)
           iDom = cumDomProc(iProc) + nn
           bufSize = dims(1, iDom)*dims(2, iDom)*dims(3,iDom)*3
           call MPI_Recv(buffer, bufSize, sumb_real, iProc, iProc, &
                sumb_comm_world, status, ierr)

           write(tmpStr, "(a,I2.2,a,I2.2,a)"), """Proc ", iProc, " Local ID ", nn ,""""
           write(1,*) "ZONE I=", dims(1, iDom), " J=", dims(2, iDom), "K=", dims(3, iDom), "T=", trim(tmpStr)
           write(1, *) "DATAPACKING=BLOCK"
           write(1, *) "VARLOCATION=([1,2,3]=NODAL)"

           ! Dump directly...already in the right order
           do i=1, bufSize
              write(1, *), buffer(i)
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
     end do
  end if

  deallocate(buffer, nDomProc, cumDomProc, dims)

end subroutine writePartionedMesh


 ! ! ALlocate to the range of the data poitner in the CSR matrix. In
 !  ! theory this should make indexing into the array easier. (hopefully)
 !  allocate(mySearchCoords(overlap%rowPtr(startDom): overlap%rowPtr(endDom+1)-1))

 !  ! Set all the sizes to 0:
 !  do i=overlap%rowPtr(startDom), overlap%rowPtr(endDom+1)-1
 !     mySearchCoords(i)%n = 0
 !  end do

 !  do nn=1, nDom
 !     call setPointers(nn, level, sps)
 !     iDom = cumDomProc(myid) + nn

 !     ! Allocate approximate space for the coordinates and local fringes
 !     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
 !        ! Allocate to 1/8 of the maximum number...we will double it if
 !        ! necessary. That means at most 3 doublings. 

 !        ii = (il-1)*(jl-1)*(kl-1) / 8
 !        allocate(mySearchCoords(jj)%x(3, ii), &
 !             mySearchCoords(jj)%fInd(3, ii))
 !        mySearchCoords(jj)%arrSize = ii
 !     end do

 !     do k=2, kl
 !        do j=2, jl
 !           do i=2, il

 !              ! Generate the search point.
 !              do iDim=1, 3
 !                 xpt(iDim) = eighth*(&
 !                      x(i-1, j-1, k-1, iDim) + &
 !                      x(i  , j-1, k-1, iDim) + &
 !                      x(i-1, j  , k-1, iDim) + &
 !                      x(i  , j  , k-1, iDim) + &
 !                      x(i-1, j-1, k  , iDim) + &
 !                      x(i  , j-1, k  , iDim) + &
 !                      x(i-1, j  , k  , iDim) + &
 !                      x(i  , j  , k  , iDim))
 !              end do

 !              ! Loop over the other overlapping blocks
 !              do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
 !                 jDom = overlap%colInd(jj)


 !                 ! Do the search point overlap the bounding box?
 !                 if ( xpt(1) >= xMin(1, jDom) .and. &
 !                      xpt(1) <= xMax(1, jDom) .and. &
 !                      xpt(2) >= xMin(2, jDom) .and. &
 !                      xpt(2) <= xMax(2, jDom) .and. &
 !                      xpt(3) >= xMin(3, jDom) .and. &
 !                      xpt(3) <= xMax(3, jDom)) then 

 !                    ! There is a chance that we can find a donor
 !                    ! since there is least one cell with a larger
 !                    ! volume than mine. Or if we're a forced
 !                    ! receiver, we got to do it anyway. 
 !                    if (vol(i, j, k) > minVol(jDom)  .or. forceRecv(i, j, k)) then

 !                       ! This will need to be searched
 !                       overlap%data(jj) = overlap%data(jj) + 1 

 !                       ! Since we did a non-trivial amount of work to
 !                       ! figure this out...we will save the information.                        
 !                       mySearchCoords(jj)%n = mySearchCoords(jj)%n + 1

 !                       if (mySearchCoords(jj)%n > mySearchCoords(jj)%arrSize) then 

 !                          ! Double the length: We we a little pointer
 !                          ! magic here to do only 1 copy instead of 2. 
 !                          real2DPtr => mySearchCoords(jj)%x
 !                          int2DPtr  => mySearchCoords(jj)%fInd
 !                          oldSize = mySearchCoords(jj)%arrSize
 !                          newSize = oldSize * 2

 !                          ! This allocated *new* memory.
 !                          allocate(mySearchCoords(jj)%x(3, newSize), & 
 !                               mySearchCoords(jj)%fInd(3, newSize))

 !                          ! Copy old values. 
 !                          mySearchCoords(jj)%x(:, 1:oldSize) = real2DPtr
 !                          mySearchCoords(jj)%fInd(:, 1:oldSize) = int2DPtr

 !                          ! The magic is we deallocate real2Dptr and
 !                          ! int2Dptr, which is really deallocating the
 !                          ! *original* memory. 
 !                          deallocate(real2DPtr, int2DPtr)

 !                          ! Set the new size:
 !                          mySearchCoords(jj)%arrSize = newSize

 !                       end if

 !                       ! Now it is safe to add
 !                       mySearchCoords(jj)%x(:, mySearchCoords(jj)%n) = xpt
 !                       mySearchCoords(jj)%fInd(:, mySearchCoords(jj)%n) = (/i, j, k/)

 !                    end if
 !                 end if
 !              end do
 !           end do
 !        end do
 !     end do

 !     ! We can now allocate the remainder of the space requied
 !     ! exactly. We won't worry about the additional unused space in x and fInd.
 !     do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
 !        n = mySearchCoords(jj)%n
 !        allocate(mySearchCoords(jj)%dInd(4, n), mySearchCoords(jj)%gInd(8, n))
 !        allocate(mySearchCoords(jj)%frac(3, n))

 !        ! Initialize so we can check later if a donor wasn't found
 !        mySearchCoords(jj)%dInd = -1
 !        mySearchCoords(jj)%gInd = -1
 !        mySearchCoords(jj)%frac = -one
        
 !     end do
 !  end do

! subroutine allocateOverset(nn, level, sps)

!   use blockPointers
!   implicit none

!   ! Allocate (if necessary) and initialize the on-block overset
!   ! information for the requested block. 

!   ! Input Parameters
!   integer(kind=intType) :: nn, level, sps

!   ! Working parameters
!   integer(kind=intType) :: i, j, k, iDim

!   call setPointers(nn, level, sps)

!   if (.not. associated(flowDoms(nn, level, sps)%donors)) then 
!      allocate(flowDoms(nn, level, sps)%donors(0:ib, 0:jb, 0:kb))
!   end if

!   if (.not. associated(flowDoms(nn, level, sps)%forceRecv)) then 
!      allocate(flowDoms(nn, level, sps)%forceRecv(2:il, 2:jl, 2:kl))
!   end if

!   if (.not. associated(flowDoms(nn, level, sps)%recvStatus)) then 
!      allocate(flowDoms(nn, level, sps)%recvStatus(2:il, 2:jl, 2:kl))
!   end if

!   ! if (.not. associated(flowDoms(nn, level, sps)%xSearch)) then 
!   !    allocate(flowDoms(nn, level, sps)%xSearch(3, 2:il, 2:jl, 2:kl))
!   ! end if

!   if (.not. associated(flowDoms(nn, level, sps)%qualRecv)) then 
!      allocate(flowDoms(nn, level, sps)%qualRecv(2:il, 2:jl, 2:kl))
!   end if

!   ! While we're here, compute xSearch which is just the cell centers. 
!   do k=2, kl
!      do j=2, jl
!         do i=2, il
!            ! do iDim=1, 3
!            !    flowDoms(nn, level, sps)%xSearch(iDim, i, j, k) = eighth*(&
!            !         x(i-1, j-1, k-1, iDim) + &
!            !         x(i  , j-1, k-1, iDim) + &
!            !         x(i-1, j  , k-1, iDim) + &
!            !         x(i  , j  , k-1, iDim) + &
!            !         x(i-1, j-1, k  , iDim) + &
!            !         x(i  , j-1, k  , iDim) + &
!            !         x(i-1, j  , k  , iDim) + &
!            !         x(i  , j  , k  , iDim))
!            ! end do

!            ! And set qualRecv which is just our own volume
!            flowDoms(nn, level, sps)%qualRecv(i, j, k) = vol(i, j, k)
!         end do
!      end do
!   end do

!   ! Initialize all the data in the donor derived type
!   do k=2, kl
!      do j=2, jl
!         do i=2, il
!            flowDoms(nn, level, sps)%donors(i, j, k)%donorProcID = -1
!            flowDoms(nn, level, sps)%donors(i, j, k)%donorBlockID = -1
!            flowDoms(nn, level, sps)%donors(i, j, k)%frac = zero
!            flowDoms(nn, level, sps)%donors(i, j, k)%ind = -1
!            flowDoms(nn, level, sps)%donors(i, j, k)%gInd = -1
!         end do
!      end do
!   end do

!   ! These default to false. Will be modified as necessary before
!   ! searching.
!   flowDoms(nn, level, sps)%recvStatus = .False.
!   flowDOms(nn, level, sps)%forceRecv = .False. 

! end subroutine allocateOverset


  ! ! -----------------------------------------------------------------
  ! ! Step 12: On-block preprocessing: Flag cells near holes as
  ! ! forceRecv as well as cells that already have iblank=-1
  ! ! -----------------------------------------------------------------

  ! ! We have to allocate some data for our own blocks we will be
  ! ! searching. Do that while we're waiting for the data to arrive:

  ! do nn=1, nDom
  !    call allocateOverset(nn, level, sps)
  ! end do

  ! do nn=1,nDom
  !    call setPointers(nn, level, sps)
  !    do k=2, kl
  !       do j=2, jl
  !          do i=2, il

  !             select case (iBlank(i, j, k))

  !             case (0) 

  !                ! We have to be careful with holes: We need to make
  !                ! sure that there is sufficient frige padding around
  !                ! them. Therefore for each real cell we loop over the
  !                ! cells in its stencil. If a cell its stencil is a
  !                ! hole, then the cell in question MUST be forced to
  !                ! be a receiver. Note there is no issue with bounds,
  !                ! since iBlank is double haloed and we're only
  !                ! looping over the owned cells. 

  !                stencilLoop: do i_stencil=1,  N_visc_drdw
  !                   ii = visc_drdw_stencil(i_stencil, 1) + i
  !                   jj = visc_drdw_stencil(i_stencil, 2) + j
  !                   kk = visc_drdw_stencil(i_stencil, 3) + k

  !                   if (iblank(ii, jj, kk) /= 0) then 
  !                      forceRecv(i, j, k) = .True. 
  !                   end if
  !                end do stencilLoop
  !             case (-1)
  !                ! Cell as ben flagged as necessairly being a
  !                ! receiver so force it to be. 
  !                forceRecv(i, j, k) = .True. 
  !             end select
  !          end do
  !       end do
  !    end do
  ! end do




  ! sendCount = 0
  ! do nn=1, nDom 
  !    ! block in gloabl ordering. 
  !    iDom = cumDomProc(myid) + nn  

  !    do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
  !       jDom = overlap%colInd(jj)

  !       if (overlap%data(jj) /= 0 .and. blkProc(jDom) /= myid) then 
  !          n = mySearchCoords(jj)%n
  !          call mpi_isend(mySearchCoords(jj)%x, 3*n, sumb_real, blkProc(jDom), &
  !               jj, SUmb_comm_world, sendRequests(sendCount+1), ierr)
  !          call ECHK(ierr, __FILE__, __LINE__)
  !          sendCount = sendCount + 1
  !       end if
  !    end do
  ! end do

  ! ! Now post all the receives...these block, but we use MPI_ANY_SOURCE
  ! ! and MPI_ANY_TAG
  ! do iDom=1, nDomTotal
  !    if (blkProc(iDom) /= myid) then 

  !       do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
  !          jDom = overlap%colInd(jj)

  !          if (jDom >= startDom .and. jDom <= endDom) then 
  !             if (overlap%data(jj) /= 0) then 

  !                call mpi_recv(realRecvBuffer, rSize*3, sumb_real, MPI_ANY_SOURCE, &
  !                     MPI_ANY_TAG, SUmb_comm_world, status, ierr)
  !                call ECHK(ierr, __FILE__, __LINE__)

  !                ! I don't care where it came from...the TAG, which is
  !                ! the jj from the sending process tells me exactly
  !                ! where to put it. 

  !                jj1 = status(MPI_TAG)
  !                call MPI_Get_count(status, sumb_real, n, ierr)
  !                call ECHK(ierr, __FILE__, __LINE__)

  !                allocate(&
  !                     otherSearchCoords(jj1)%x(3, n/3), &
  !                     otherSearchCoords(jj1)%dInd(3, n/3), &
  !                     otherSearchCoords(jj1)%gInd(8, n/3), &
  !                     otherSearchCoords(jj1)%frac(3, n/3))
  !                otherSearchCoords(jj1)%n = n/3
  !                otherSearchCoords(jj1)%arrSize = n/3
                 
  !                ! Initialize so we know when something wasn't found
  !                otherSearchCoords(jj1)%dInd = -1
  !                otherSearchCoords(jj1)%gInd = -1
  !                otherSearchCoords(jj1)%frac = -one

  !                ! Copy the data in so we can reuse the buffer:
  !                do i=1, n/3
  !                   otherSearchCoords(jj1)%x(:, i) = realRecvBuffer(3*i-2:3*i)
  !                end do
  !             end if
  !          end if
  !       end do
  !    end if
  ! end do

  ! ! Wait for all the send requests for finish. Perhaps could use an
  ! ! mpi_waitall here? Probably doesn't really matter.
  ! do i=1, sendCount
  !    call mpi_waitany(sendCount, sendRequests, index, status, ierr)
  !    call ECHK(ierr, __FILE__, __LINE__)
  ! enddo

  ! ! In theory we can now *FINALLY* do the search in parallel....We
  ! ! have our own ADTrees along with the search coords from myself and
  ! ! the search coords from the other procs. Essentially we are
  ! ! searching for stuff in our proc's COLUMN

  ! do iDom=1, nDomTotal
  !    do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
  !       jDom = overlap%colInd(jj)

  !       if (jDom >= startDom .and. jDom <= endDom) then 
  !          if (overlap%data(jj) /= 0) then

  !             if (blkProc(iDom) == myid) then 
  !                ! We can use my own search Coords:

  !                call newSearch(oBlocks(jDom), mySearchCoords(jj))

  !             else
  !                ! Otherwise use the ones that we have got from
  !                ! another block
  !                call newSearch(oBlocks(jDom), otherSearchCoords(jj))

  !             end if
  !          end if
  !       end if
  !    end do
  ! end do



  ! ! Now essentially do the reverse of the communication...we send
  ! ! stuff from otherSearchCoords to mySearchCoords.

  ! sendCount = 0

  ! ! Now post all the receives...these block, but we use MPI_ANY_SOURCE
  ! ! and MPI_ANY_TAG
  ! do iDom=1, nDomTotal
  !    if (blkProc(iDom) /= myid) then 

  !       do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
  !          jDom = overlap%colInd(jj)

  !          if (jDom >= startDom .and. jDom <= endDom) then 
  !             if (overlap%data(jj) /= 0) then 
  !                n = otherSearchCoords(jj)%n
  !                call mpi_isend(otherSearchCoords(jj)%frac, 3*n, sumb_real, blkProc(iDom), &
  !                     jj, SUmb_comm_world, sendRequests(sendCount+1), ierr)
  !                call ECHK(ierr, __FILE__, __LINE__)
  !                sendCount = sendCount + 1

  !             end if
  !          end if
  !       end do
  !    end if
  ! end do

  ! do nn=1, nDom 
  !    ! block in gloabl ordering. 
  !    iDom = cumDomProc(myid) + nn  

  !    do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
  !       jDom = overlap%colInd(jj)

  !       if (overlap%data(jj) /= 0 .and. blkProc(jDom) /= myid) then 


  !          call mpi_recv(realRecvBuffer, rSize*3, sumb_real, MPI_ANY_SOURCE, &
  !               MPI_ANY_TAG, SUmb_comm_world, status, ierr)
  !          call ECHK(ierr, __FILE__, __LINE__)
           
  !          jj1 = status(MPI_TAG)
  !          n = mySearchCoords(jj1)%n
  !          do i=1, n
  !             mySearchCoords(jj1)%frac(:, i) = realRecvBuffer(3*i-2:3*i)
  !          end do
  !       end if
  !    end do
  ! end do

  ! ! Wait for all the send requests for finish. Perhaps could use an
  ! ! mpi_waitall here? Probably doesn't really matter.
  ! do i=1, sendCount
  !    call mpi_waitany(sendCount, sendRequests, index, status, ierr)
  !    call ECHK(ierr, __FILE__, __LINE__)
  ! enddo

  
  ! timeB = mpi_wtime()
  ! print *,'myid time 2:', myid, timeB-timeA
  ! call MPi_barrier(sumb_comm_world, ierr)
  ! stop
