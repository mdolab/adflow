
subroutine emptyFringe(fringe)

  use block
  implicit none
  ! Input/Output
  type(fringeType), intent(inout) :: fringe

  ! Initialize data in empty fringe
  fringe%quality = large
  fringe%donorProc = -1
  fringe%donorBlock = -1
  fringe%dI = -1
  fringe%dJ = -1
  fringe%dK = -1
  fringe%myBlock = -1
  fringe%myI = -1
  fringe%myJ = -1
  fringe%myK = -1
  fringe%donorFrac = -one
  fringe%gInd = -1
  fringe%status = 0
  call setIsCompute(fringe%status, .True.)
end subroutine emptyFringe

subroutine printOverlapMatrix(overlap)

  ! This is a debugging routine to print out the overlap matrix.

  use communication
  use overset
  implicit none

  ! Input/output
  type(CSRMatrix), intent(in) :: overlap

  ! Working
  integer(kind=intType) :: i, jj

  if (myid == 0) then 
     ! Now dump out who owns what:
     do i=1, overlap%nrow
        write(*, "(a,I4, a)", advance='no'), 'Row:', i, "   "
        do jj=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
           write(*, "(a,I2, a, e10.5)", advance='no'), "(", overlap%colInd(jj), ")", overlap%data(jj)
        end do
        write(*, *) " "
     end do

     print *, '--------------------------------------'
     ! Now dump out who owns what:
     do i=1, overlap%nRow
        write(*, "(a,I4, a)", advance='no'), 'Row:', i, "   "
        do jj=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
           write(*, "(a,I2, a, I8)", advance='no'), "(", overlap%colInd(jj), ")", int(overlap%assignedProc(jj))
        end do
        write(*, *) " "
     end do
  end if
end subroutine printOverlapMatrix

subroutine getCumulativeForm(sizeArray, n, cumArray)
  use precision
  implicit none


  ! Input/Output
  integer(kind=intType), dimension(n), intent(in) :: sizeArray
  integer(kind=intType), intent(in) :: n
  integer(kind=intType), dimension(0:n), intent(out) :: cumArray

  ! Working
  integer(kind=intType) :: i

  cumArray(0) = 0
  do i=1, n
     cumArray(i) = cumArray(i-1) + sizeArray(i)
  end do

end subroutine getCumulativeForm

subroutine transposeOverlap(A, B)

  ! Create the matrix Create the matrix transpose. 
  ! Inspired by: https://people.sc.fsu.edu/~jburkardt/f_src/sparsekit/sparsekit.f90
  use overset
  implicit none

  ! Input/Output
  type(CSRMatrix), intent(in) :: A
  type(CSRMatrix), intent(inout) :: B

  ! Working
  integer(kind=intType) :: col, colp1, i, k, next

  ! A CSR matrix is the same as a CSC matrix of the transpose. So
  ! essentially the algorithm is convert A as a CSC matrix to B
  ! (a CSR matrix)

  ! Allocate space for everything in B
  B%nnz = A%nnz
  B%nRow = A%nCol
  B%nCol = A%nRow
  allocate(B%data(B%nnz), B%colInd(B%nnz), &
       B%assignedProc(B%nnz), B%rowPtr(B%nRow + 1))
  B%allocated = .True.
  !  Compute lengths of rows of B (ie the columns of A)

  B%rowPtr = 0

  do i = 1, A%nRow
     do k = A%rowPTr(i), A%rowPtr(i+1)-1
        colp1 = A%colInd(k) +1
        B%rowPtr(colp1) = B%rowPtr(colp1) + 1

     end do
  end do
  !
  !  Compute pointers from lengths.
  !
  B%rowPtr(1) = 1
  do i = 1, A%nRow
     B%rowPtr(i+1) = B%rowPtr(i) + B%rowPtr(i+1)
  end do
  !
  !  Do the actual copying.
  !
  do i = 1, A%nRow
     do k = A%rowPtr(i), A%rowPtr(i+1)-1
        col = A%colInd(k)
        next = B%rowPtr(col)

        B%data(next) = A%data(k)
        B%assignedProc(next) = A%assignedProc(k) 

        B%colInd(next) = i 
        B%rowPtr(col) = next + 1
     end do
  end do

  !  Reshift B%rowPtr
  do i = A%nRow, 1, -1
     B%rowPtr(i+1) = B%rowPtr(i)
  end do
  B%rowPtr(1) = 1 

end subroutine transposeOverlap

subroutine deallocateCSRMatrix(mat1)

  use overset
  implicit none

  type(CSRMatrix), intent(inout) :: mat1

  if (mat1%allocated) then 
     deallocate(&
          mat1%data, &
          mat1%colInd, &
          mat1%rowPtr, &
          mat1%assignedProc)
  end if

end subroutine deallocateCSRMatrix

subroutine computeFringeProcArray(fringes, n, fringeProc, cumFringeProc, nFringeProc)
  ! Compute the breakpoints "cumFringeProc" for a list of sorted n
  ! fringes "fringes". nFringeProc is the total number of unique
  ! processors. fringeProc is the processor number for each section.
  use block
  use communication
  implicit none

  ! Input/Output
  type(fringeType), intent(in), dimension(n) :: fringes
  integer(kind=intType), intent(in) :: n
  integer(kind=intType), intent(out) :: nFringeProc
  integer(kind=intType), intent(out) :: fringeProc(nProc), cumFringeProc(1:nProc+1)

  ! Working
  integer(kind=intType) :: i, currentProc

  fringeProc = -1
  nFringeProc = 0
  cumFringeProc(1) = 1
  currentProc = -1

  do i=1, n
     if (currentProc /= fringes(i)%donorProc) then 
        nFringeProc = nFringeProc + 1
        cumFringeProc(nFringeProc) = i
        fringeProc(nFringeProc) = fringes(i)%donorProc
        currentProc = fringes(i)%donorProc
     end if
  end do

  ! Finally, the nFringeProc+1 entry is always n+1
  cumFringeProc(nFringeProc+1) = n + 1

end subroutine computeFringeProcArray

subroutine fracToWeights(frac, weights)
  use constants
  implicit none
  real(kind=realType), intent(in), dimension(3) :: frac
  real(kind=realType), intent(out), dimension(8) :: weights

  weights(1) = (one-frac(1))*(one-frac(2))*(one-frac(3))
  weights(2) = (    frac(1))*(one-frac(2))*(one-frac(3))
  weights(3) = (one-frac(1))*(    frac(2))*(one-frac(3))
  weights(4) = (    frac(1))*(    frac(2))*(one-frac(3))
  weights(5) = (one-frac(1))*(one-frac(2))*(    frac(3))
  weights(6) = (    frac(1))*(one-frac(2))*(    frac(3))
  weights(7) = (one-frac(1))*(    frac(2))*(    frac(3))
  weights(8) = (    frac(1))*(    frac(2))*(    frac(3))
end subroutine fracToWeights

subroutine getCommPattern(oMat,  sendList, size1, nSend, recvList, size2, nRecv)

  use overset
  use communication 
  implicit none

  ! Input/output
  type(CSRMatrix), intent(in) :: oMat
  integer(kind=intType), intent(in) :: size1, size2
  integer(kind=intType), intent(out) :: sendList(2, size1), recvList(2, size2)
  integer(kind=intType), intent(out) :: nSend, nRecv

  ! Working:
  integer(kind=intType) :: nn, iDom, nnRow, i, jj, ii, nUniqueProc, iProc
  integer(kind=intType), dimension(:), allocatable :: procsForThisRow, inverse, blkProc
  logical :: added

  ! Generic routine to determine what I need to send/recv based on the
  ! data provided in the overlap matrix. The '2' in the send and
  ! receive lists will record the processor and the global 'idom'
  ! index, which is suffient to use for the subsequent communication
  ! structure.

  nSend = 0
  nRecv = 0

  ! These variables are used to compact the sending of
  ! blocks/fringes. The logic is as follows: A a different
  ! processor may need a block/fringe for more than 1 search. This
  ! is manifested by having two or more entries is the rows (or
  ! columns) I own. It would be inefficient to send the same data
  ! to the same processor more than once, so we "uniquify" the
  ! processors before we send. There is also another salient
  ! reason: If we were to send the same data twice, and the other
  ! processor started using the data, we could get a race condition
  ! as it was modified the received fringes (during a search) while
  ! the same fringes were being overwritten by the receive operation. 

  allocate(procsForThisRow(nDomTotal), inverse(nDomTotal), blkProc(nDomTotal))

  ii = 0
  do iProc=0, nProc-1
     do i=1, nDomProc(iProc)
        ii = ii + 1
        blkProc(ii) = iProc
     end do
  end do

  ! Loop over the owned rows of the normal matrix
  do nn=1, nDom
     iDom = cumDomProc(myid) + nn
     nnRow = oMat%rowPtr(iDom+1) - oMat%rowPtr(iDom)
     procsForThisRow(1:nnRow) = oMat%assignedProc(oMat%rowPtr(iDom) : oMat%rowPtr(iDom+1)-1)
     call unique(procsForThisRow, nnRow, nUniqueProc, inverse)

     do jj = 1, nUniqueProc
        if (procsForThisRow(jj) /= myid) then 
           ! This intersection requires a row quantity from me
           nSend = nSend + 1
           sendList(1, nSend) = procsForThisRow(jj) 
           sendList(2, nSend) = iDom 
        end if
     end do
  end do

  ! Now we loop back through the whole matrix looking at what I have
  ! to do. If there is a row I don't own, I will have to receive it:
 
  do iDom=1, oMat%nRow
     added = .False.
     rowLoop: do jj=oMat%rowPtr(iDom), oMat%rowPtr(iDom+1)-1

        ! I have to do this intersection
        if (oMat%assignedProc(jj) == myID) then 

           ! But I don't know the row entry
           if (.not. (iDom > cumDomProc(myid) .and. iDom <= cumDomProc(myid+1))) then 

              ! Need to back out what proc the iDom correponds to:
              nRecv = nRecv + 1
              recvList(1, nRecv) = blkProc(iDom)
              recvList(2, nRecv) = iDom
              added = .True.
           end if
        end if

        ! Just move on to the next row since we only need to receive it once. 
        if (added) then 
           exit rowLoop
        end if
     end do rowLoop
  end do

  deallocate(procsForThisRow, inverse, blkProc)
end subroutine getCommPattern

subroutine getOWallCommPattern(oMat, oMatT, sendList, size1, nSend, &
     recvList, size2, nRecv, rBufSize)
  
  ! This subroutine get the the comm pattern to send the oWall types. 
  
  use overset
  use communication 
  implicit none

  ! Input/output
  type(CSRMatrix), intent(in) :: oMat, oMatT
  integer(kind=intType), intent(in) :: size1, size2
  integer(kind=intType), intent(out) :: sendList(2, size1), recvList(2, size2)
  integer(kind=intType), intent(out) :: nSend, nRecv
  integer(kind=intType), intent(in) :: rBufSize(nDomTotal)

  ! Working:
  integer(kind=intType) :: nn, iDom, jDom, nnRow, nnRowT, i, jj, ii, nUniqueProc, iProc
  integer(kind=intType), dimension(:), allocatable :: blkProc, toRecv
  integer(kind=intType), dimension(:), allocatable :: procsForThisRow, inverse


  nSend = 0
  nRecv = 0

  allocate(procsForThisRow(2*nDomTotal), inverse(2*nDomTotal), blkProc(nDomTotal))

  ii = 0
  do iProc=0, nProc-1
     do i=1, nDomProc(iProc)
        ii = ii + 1
        blkProc(ii) = iProc
     end do
  end do

  ! Loop over the owned rows of the regular matrix and the rows
  ! transposed matrix (ie columns of the regular matrix)
  do nn=1, nDom

     iDom = cumDomProc(myid) + nn

     ! Only deal with this block if the rbuffer size for the oWall is
     ! greater than zero. If it is zero, it is empty we don't need to
     ! deal with it.
     if (rBufSize(iDom) > 0) then 

        nnRow = oMat%rowPtr(iDom+1) - oMat%rowPtr(iDom)
        procsForThisRow(1:nnRow) = oMat%assignedProc(oMat%rowPtr(iDom) : oMat%rowPtr(iDom+1)-1)
        
        nnRowT = oMatT%rowPtr(iDom+1) - oMatT%rowPtr(iDom)
        procsForThisRow(nnRow+1:nnRow+nnRowT) = oMatT%assignedProc(oMatT%rowPtr(iDom) : oMatT%rowPtr(iDom+1)-1)

        call unique(procsForThisRow, nnRow+nnRowT, nUniqueProc, inverse)
        
        do jj = 1, nUniqueProc
           if (procsForThisRow(jj) /= myid) then 
              ! This intersection requires a row quantity from me
              nSend = nSend + 1
              sendList(1, nSend) = procsForThisRow(jj) 
              sendList(2, nSend) = iDom 
           end if
        end do
     end if
  end do

  ! Now we loop back through the whole matrix looking at what I have
  ! to do. If there is a row or column I don't own, I will have to receive it:
  allocate(toRecv(nDomTotal))
  toRecv = 0
  do iDom=1, nDomTotal
     do jj=oMat%rowPtr(iDom), oMat%rowPtr(iDom+1)-1
        jDom = oMat%colInd(jj)
        if (oMat%assignedProc(jj) == myID) then 
           ! I don't have the row entry:
           if (.not. (iDom > cumDomProc(myid) .and. iDom <= cumDomProc(myid+1))) then 
              toRecv(iDom) = 1
           end if
           ! Don't have the column entry:
           if (.not. (jDom > cumDomProc(myid) .and. jDom <= cumDomProc(myid+1))) then 
              toRecv(jDom) = 1
           end if
        end if
     end do
  end do
     
  ! Now loop back through and set my recvList. Only add if the
  ! rBufferSize is larger than zero.
  do iDom=1, nDomTotal
     if (toRecv(iDom) == 1 .and. rBufSize(iDom) > 0) then 
        nRecv = nRecv + 1
        recvList(1, nRecv) = blkProc(iDom)
        recvList(2, nRecv) = iDom
     end if
  end do

  deallocate(procsForThisRow, inverse, blkProc, toRecv)
end subroutine getOWallCommPattern

subroutine sendOBlock(oBlock, iDom, iProc, tagOffset, sendCount)

  use communication
  use overset
  implicit none

  ! Input/Output
  type(oversetBlock), intent(inout) :: oBlock
  integer(kind=intType), intent(in) :: iProc, iDom, tagOffset
  integer(kind=intType), intent(inout) :: sendCount
  
  ! Working
  integer(kind=intType) :: tag, ierr

  tag = tagOffset + iDom
  sendCount = sendCount + 1
  call mpi_isend(oBlock%rBuffer, size(oBlock%rbuffer), sumb_real, &
       iProc, tag, SUmb_comm_world, sendRequests(sendCount), ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  
  sendCount = sendCount + 1
  call mpi_isend(oBlock%iBuffer, size(oBlock%iBuffer), sumb_integer, &
       iProc, tag, SUmb_comm_world, sendRequests(sendCount), ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  
end subroutine sendOBlock

subroutine sendOFringe(oFringe, iDom, iProc, tagOffset, sendCount)

  use communication
  use overset
  implicit none

  ! Input/Output
  type(oversetFringe), intent(inout) :: oFringe
  integer(kind=intType), intent(in) :: iProc, iDom, tagOffset
  integer(kind=intType), intent(inout) :: sendCount

  ! Working
  integer(kind=intType) :: tag, ierr

  tag = iDom + tagOffset
  sendCount = sendCount + 1
  call mpi_isend(oFringe%rBuffer, size(oFringe%rbuffer), sumb_real, &
       iProc, tag, SUmb_comm_world, sendRequests(sendCount), ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  
  sendCount = sendCount + 1
  call mpi_isend(oFringe%iBuffer, size(oFringe%iBuffer), sumb_integer, &
       iProc, tag, SUmb_comm_world, sendRequests(sendCount), ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  
end subroutine sendOFringe

subroutine sendOWall(oWall, iDom, iProc, tagOffset, sendCount)

  use communication
  use overset
  implicit none

  ! Input/Output
  type(oversetWall), intent(inout) :: oWall
  integer(kind=intType), intent(in) :: iProc, iDom, tagOffset
  integer(kind=intType), intent(inout) :: sendCount

  ! Working
  integer(kind=intType) :: tag, ierr

  tag = iDom + tagOffset
  sendCount = sendCount + 1
  call mpi_isend(oWall%rBuffer, size(oWall%rbuffer), sumb_real, &
       iProc, tag, SUmb_comm_world, sendRequests(sendCount), ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  
  sendCount = sendCount + 1
  call mpi_isend(oWall%iBuffer, size(oWall%iBuffer), sumb_integer, &
       iProc, tag, SUmb_comm_world, sendRequests(sendCount), ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  
end subroutine sendOWall

subroutine recvOBlock(oBlock, iDom, iProc, tagOffset, iSize, rSize, &
     recvCount, recvInfo)

  use communication
  use overset
  implicit none

  ! Input/Output
  type(oversetBlock), intent(inout) :: oBlock
  integer(kind=intType), intent(in) :: iDom, iProc, tagOffset, rSize, iSize
  integer(kind=intType), intent(inout) :: recvCount
  integer(kind=intType), intent(inout) :: recvInfo(2, recvCount+2)

  ! Working
  integer(kind=intType) :: tag, ierr

  tag = tagOffset + iDom
  allocate(oBLock%rBuffer(rSize), oBlock%iBuffer(iSize))
 
  recvCount = recvCount + 1
  call mpi_irecv(oBlock%rBuffer, rSize, sumb_real, &
       iProc, tag, SUmb_comm_world, recvRequests(recvCount), ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  recvInfo(:, recvCount) = (/iDom, 1/) 

  recvCount = recvCount + 1
  call mpi_irecv(oBlock%iBuffer, iSize, sumb_integer, &
       iProc, tag, SUmb_comm_world, recvRequests(recvCount), ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  recvInfo(:, recvCount) = (/iDom, 2/) 
  
end subroutine recvOBlock

subroutine recvOFringe(oFringe, iDom, iProc, tagOffset, iSize, rSize, &
     recvCount, recvInfo)

  use communication
  use overset
  implicit none

  ! Input/Output
  type(oversetFringe), intent(inout) :: oFringe
  integer(kind=intType), intent(in) :: iDom, iProc, tagOffset, rSize, iSize
  integer(kind=intType), intent(inout) :: recvCount
  integer(kind=intType), intent(inout) :: recvInfo(2, recvCount+2)

  ! Working
  integer(kind=intType) :: tag, ierr

  tag = tagOffset + iDom
  allocate(oFringe%rBuffer(rSize), oFringe%iBuffer(iSize))
 
  recvCount = recvCount + 1
  call mpi_irecv(oFringe%rBuffer, rSize, sumb_real, &
       iProc, tag, SUmb_comm_world, recvRequests(recvCount), ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  recvInfo(:, recvCount) = (/iDom, 3/) 

  recvCount = recvCount + 1
  call mpi_irecv(oFringe%iBuffer, iSize, sumb_integer, &
       iProc, tag, SUmb_comm_world, recvRequests(recvCount), ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  recvInfo(:, recvCount) = (/iDom, 4/) 
  
end subroutine recvOFringe

subroutine recvOWall(oWall, iDom, iProc, tagOffset, iSize, rSize, &
     recvCount, recvInfo)

  use communication
  use overset
  implicit none

  ! Input/Output
  type(oversetWall), intent(inout) :: oWall
  integer(kind=intType), intent(in) :: iDom, iProc, tagOffset, rSize, iSize
  integer(kind=intType), intent(inout) :: recvCount
  integer(kind=intType), intent(inout) :: recvInfo(2, recvCount+2)

  ! Working
  integer(kind=intType) :: tag, ierr

  tag = tagOffset + iDom
  allocate(oWall%rBuffer(rSize), oWall%iBuffer(iSize))
 
  recvCount = recvCount + 1
  call mpi_irecv(oWall%rBuffer, rSize, sumb_real, &
       iProc, tag, SUmb_comm_world, recvRequests(recvCount), ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  recvInfo(:, recvCount) = (/iDom, 5/) 

  recvCount = recvCount + 1

  call mpi_irecv(oWall%iBuffer, iSize, sumb_integer, &
       iProc, tag, SUmb_comm_world, recvRequests(recvCount), ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  recvInfo(:, recvCount) = (/iDom, 6/) 
  
end subroutine recvOWall

subroutine deallocateOBlocks(oBlocks, n)

  ! This subroutine deallocates all data stores in a list of oBlocks
  use adtAPI
  use overset
  implicit none

  ! Input Params
  type(oversetBlock), dimension(n), intent(inout) :: oBLocks
  integer(kind=intType) :: n

  ! Working Parameters
  integer(kind=intType) :: i

  do i=1, n

     ! oBlock:
     if (oblocks(i)%allocated) then 
        deallocate(&
             oBlocks(i)%hexaConn, &
             oBlocks(i)%globalCell, &
             oBLocks(i)%nearWall, &
             oBLocks(i)%invalidDonor, &
             oBlocks(i)%qualDonor, &
             oBlocks(i)%xADT)
        if (allocated(oblocks(i)%rbuffer)) then 
           deallocate(oBlocks(i)%rBuffer, &
                oBlocks(i)%iBuffer)
        end if
        call destroySerialHex(oBlocks(i)%ADT)
     end if
  end do
end subroutine deallocateOBlocks

subroutine deallocateOFringes(oFringes, n)

  ! This subroutine deallocates all data stores in a list of oFringes

  use adtAPI
  use overset
  implicit none

  ! Input Params
  type(oversetFringe), dimension(n), intent(inout) :: oFringes
  integer(kind=intType) :: n

  ! Working Parameters
  integer(kind=intType) :: i

  do i=1, n
     if (oFringes(i)%allocated) then 
        deallocate(&
             oFringes(i)%x, &
             oFringes(i)%quality, &
             oFringes(i)%myBlock, &
             oFringes(i)%myIndex, &
             oFringes(i)%donorProc, &
             oFringes(i)%donorBlock, &
             oFringes(i)%dI, &
             oFringes(i)%dJ, &
             oFringes(i)%dK, &
             oFringes(i)%donorFrac, &
             oFringes(i)%gInd, &
             oFringes(i)%isWall, &
             oFringes(i)%xSeed, &
             oFringes(i)%wallInd)
        if (allocated(oFringes(i)%rbuffer)) then 
           deallocate(oFringes(i)%rBuffer,&
                oFringes(i)%iBuffer)
        end if
     end if
     oFringes(i)%allocated = .False. 
  end do
end subroutine deallocateOFringes

subroutine deallocateOWalls(oWalls, n)

  ! This subroutine deallocates all data stores in a list of oWalls

  use adtAPI
  use overset
  implicit none

  ! Input Params
  type(oversetWall), dimension(n), intent(inout) :: oWalls
  integer(kind=intType) :: n

  ! Working Parameters
  integer(kind=intType) :: i

  do i=1, n
     if (oWalls(i)%allocated) then 
        deallocate(&
             oWalls(i)%x, &
             oWalls(i)%conn, &
             oWalls(i)%iblank, &
             oWalls(i)%cellPtr)
        call destroySerialQuad(oWalls(i)%ADT)
        if (oWalls(i)%nNodes > 0) then 
           call kdtree2_destroy(oWalls(i)%tree)
        end if
     end if
     oWalls(i)%allocated = .False.
  end do
end subroutine deallocateOWalls

subroutine wallsOnBlock(wallsPresent) 

  use blockPointers
  use bcTypes
  use cgnsGrid
  implicit none

  logical, intent(out) :: wallsPresent
  integer(kind=intType) :: mm
  wallsPresent = .False.
  ! Check THE ORIGINAL CGNS blocks for BCs, because the block may have
  ! been split. 
  do mm=1, cgnsDoms(nbkGlobal)%nBocos
     if (cgnsDoms(nbkGlobal)%bocoInfo(mm)%BCType == NSWallAdiabatic .or. &
          cgnsDoms(nbkGlobal)%bocoInfo(mm)%BCType == NSWallIsothermal .or. &
          cgnsDoms(nbkGlobal)%bocoInfo(mm)%BCType == EulerWall) then
        wallsPresent = .True.
     end if
  end do
end subroutine wallsOnBlock

subroutine flagForcedReceivers(tmp)

  use blockPointers
  use BCTypes
  implicit none

  ! This is generic routine for filling up a 3D array of 1st level halos
  ! cells (1:ie, 1:je, 1:ke) indicating cells that are forced
  ! receivers. BlockPointers must have already been set.

  integer(kind=intType), intent(out), dimension(1:ie, 1:je, 1:ke) :: tmp
  integer(kind=intType) :: i, j, k, mm, iStart, iEnd, jStart, jEnd, kStart, kEnd

  tmp = 0
  do mm=1,nBocos
     ! Just record the ranges necessary and we'll add in a generic
     ! loop. Why is it the first three? Well, the first level of halos
     ! off of an overset outer bound is completely
     ! meaningless. Essentially we ignore those. So the outer two
     ! layers of cells are indices 2 and 3. Therefore the first 3 on
     ! either side need to be flagged as invalid.

     select case (BCFaceID(mm))
     case (iMin)
        iStart=1; iEnd=3;
        jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
        kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
     case (iMax)
        iStart=nx; iEnd=ie;
        jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
        kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
     case (jMin)
        iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jStart=1; jEnd=3;
        kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
     case (jMax)
        iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jStart=ny; jEnd=je;
        kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
     case (kMin)
        iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
        kStart=1; kEnd=3;
     case (kMax)
        iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
        kStart=nz; kEnd=ke;
     end select

     if (BCType(mm) == OversetOuterBound) then
        do k=kStart, kEnd
           do j=jStart, jEnd
              do i=iStart, iEnd
                 tmp(i, j, k) = 1
              end do
           end do
        end do
     end if
  end do
end subroutine flagForcedReceivers

function isWallType(bType)
  
  use BCTypes
  implicit none
  integer(kind=intType) :: bType
  logical :: isWallType

  isWallType = .False.
  if (bType == NSWallAdiabatic .or. &
       bType == NSWallIsoThermal .or. &
       bType == EulerWall) then 
     isWallType = .True.
  end if

end function isWallType

! Utility function for unpacking/accessing the status variable

function isDonor(i)
  use precision
  implicit none
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
  integer(kind=intType), intent(in) :: i
  call getStatus(i, isDonor, isDonor, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
end function isDonor

function isHole(i)
  use precision
  implicit none
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
  integer(kind=intType), intent(in) :: i
  call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
end function isHole

function isCompute(i)
  use precision
  implicit none
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
  integer(kind=intType), intent(in) :: i
  call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
end function isCompute

function isFloodSeed(i)
  use precision
  implicit none
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
  integer(kind=intType), intent(in) :: i
  call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
end function isFloodSeed

function isFlooded(i)
  use precision
  implicit none
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
  integer(kind=intType), intent(in) :: i
  call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
end function isFlooded

function isWall(i)
  use precision
  implicit none
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
  integer(kind=intType), intent(in) :: i
  call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
end function isWall

function isWallDonor(i)
  use precision
  implicit none
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
  integer(kind=intType), intent(in) :: i
  call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
end function isWallDonor

subroutine setIsDonor(i, flag)
  use precision
  implicit none
  integer(kind=intType), intent(inout) :: i
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor, flag
  call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  call setStatus(i, flag   , isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
end subroutine setIsDonor

subroutine setIsHole(i, flag)
  use precision
  implicit none
  integer(kind=intType), intent(inout) :: i
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor, flag
  call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  call setStatus(i, isDonor, flag  , isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
end subroutine setIsHole

subroutine setIsCompute(i, flag)
  use precision
  implicit none
  integer(kind=intType), intent(inout) :: i
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor, flag
  call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  call setStatus(i, isDonor, isHole, flag     , isFloodSeed, isFlooded, isWall, isWallDonor)
end subroutine setIsCompute

subroutine setIsFloodSeed(i, flag)
  use precision
  implicit none
  integer(kind=intType), intent(inout) :: i
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor, flag
  call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  call setStatus(i, isDonor, isHole, isCompute, flag       , isFlooded, isWall, isWallDonor)
end subroutine setIsFloodSeed

subroutine setIsFlooded(i, flag)
  use precision
  implicit none
  integer(kind=intType), intent(inout) :: i
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor, flag
  call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  call setStatus(i, isDonor, isHole, isCompute, isFloodSeed, flag     , isWall, isWallDonor)
end subroutine setIsFlooded

subroutine setIsWall(i, flag)
  use precision
  implicit none
  integer(kind=intType), intent(inout) :: i
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor, flag
  call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  call setStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, flag  , isWallDonor)
end subroutine setIsWall

subroutine setIsWallDonor(i, flag)
  use precision
  implicit none
  integer(kind=intType), intent(inout) :: i
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor, flag
  call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  call setStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, flag)
end subroutine setIsWallDonor

subroutine setStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)

  use precision
  implicit none
  integer(kind=intType), intent(out) :: i
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
  i = 0

  if (isDonor  )   i = i + 1
  if (isHole   )   i = i + 2
  if (isCompute)   i = i + 4
  if (isFloodSeed) i = i + 8
  if (isFlooded  ) i = i + 16
  if (isWall     ) i = i + 32
  if (isWallDonor) i = i + 64

end subroutine setStatus

subroutine getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  
  use precision
  implicit none
  logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
  integer(kind=intType) :: i, j
  j = i

  isDonor = .False.
  isHole = .False.
  isCompute = .False.
  isFloodSeed = .False.
  isFlooded = .False.
  isWall = .False.
  isWallDonor = .False.

  if (j/64 > 0) then 
     isWallDonor = .True. 
     j = j - 64
  end if

  if (j/32 > 0) then 
     isWall = .True. 
     j = j - 32
  end if

  if (j/16 > 0) then 
     isFlooded = .True. 
     j = j - 16
  end if

  if (j/8 > 0) then 
     isFloodSeed = .True. 
     j = j - 8
  end if

  if (j/4 > 0) then 
     isCompute = .True. 
     j = j - 4
  end if

  if (j/2 > 0) then 
     isHole = .True. 
     j = j - 2
  end if

  if (j/1 > 0) then 
     isDonor = .True. 
     j = j - 1
  end if
end subroutine getStatus
