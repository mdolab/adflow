module oversetCommUtilities

contains

  subroutine getCommPattern(oMat,  sendList,  nSend, recvList, nRecv)

    use constants
    use oversetData, only : cumDomProc, nDomProc, nDomTotal, CSRMatrix
    use blockPointers, only : nDom
    use communication , only : nProc, myid
    use sorting, only : unique
    implicit none

    ! Input/output
    type(CSRMatrix), intent(in) :: oMat
    integer(kind=intType), intent(out) :: sendList(:, :) , recvList(:,:)
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

  subroutine getOSurfCommPattern(oMat, oMatT, sendList, nSend, &
       recvList, nRecv, rBufSize)

    ! This subroutine get the the comm pattern to send the oWall types.
    use constants
    use blockPointers, only : nDom
    use oversetData, only : nDomTotal, CSRMatrix, cumDomProc, nDomProc
    use communication, only : myid, nProc
    use sorting, only : unique
    implicit none

    ! Input/output
    type(CSRMatrix), intent(in) :: oMat, oMatT
    integer(kind=intType), intent(out) :: sendList(:, :), recvList(:, :)
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
  end subroutine getOSurfCommPattern

  subroutine sendOBlock(oBlock, iDom, iProc, tagOffset, sendCount)

    use constants
    use communication, only : adflow_comm_world, sendRequests
    use oversetData, only : oversetBlock
    use utils, only : EChk
    implicit none

    ! Input/Output
    type(oversetBlock), intent(inout) :: oBlock
    integer(kind=intType), intent(in) :: iProc, iDom, tagOffset
    integer(kind=intType), intent(inout) :: sendCount

    ! Working
    integer(kind=intType) :: tag, ierr

    tag = tagOffset + iDom
    sendCount = sendCount + 1
    call mpi_isend(oBlock%rBuffer, size(oBlock%rbuffer), adflow_real, &
         iProc, tag, ADflow_comm_world, sendRequests(sendCount), ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    sendCount = sendCount + 1
    call mpi_isend(oBlock%iBuffer, size(oBlock%iBuffer), adflow_integer, &
         iProc, tag, ADflow_comm_world, sendRequests(sendCount), ierr)
    call ECHK(ierr, __FILE__, __LINE__)

  end subroutine sendOBlock

  subroutine sendOFringe(oFringe, iDom, iProc, tagOffset, sendCount)

    use constants
    use communication, only : adflow_comm_world, sendRequests
    use oversetData, only : oversetFringe
    use utils, only : EChk
    implicit none

    ! Input/Output
    type(oversetFringe), intent(inout) :: oFringe
    integer(kind=intType), intent(in) :: iProc, iDom, tagOffset
    integer(kind=intType), intent(inout) :: sendCount

    ! Working
    integer(kind=intType) :: tag, ierr

    tag = iDom + tagOffset
    sendCount = sendCount + 1
    call mpi_isend(oFringe%rBuffer, size(oFringe%rbuffer), adflow_real, &
         iProc, tag, ADflow_comm_world, sendRequests(sendCount), ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    sendCount = sendCount + 1
    call mpi_isend(oFringe%iBuffer, size(oFringe%iBuffer), adflow_integer, &
         iProc, tag, ADflow_comm_world, sendRequests(sendCount), ierr)
    call ECHK(ierr, __FILE__, __LINE__)

  end subroutine sendOFringe

  subroutine sendOSurf(oWall, iDom, iProc, tagOffset, sendCount)

    use constants
    use communication, only : sendRequests, adflow_comm_world
    use oversetData, only : oversetWall
    use utils, only : EChk
    implicit none

    ! Input/Output
    type(oversetWall), intent(inout) :: oWall
    integer(kind=intType), intent(in) :: iProc, iDom, tagOffset
    integer(kind=intType), intent(inout) :: sendCount

    ! Working
    integer(kind=intType) :: tag, ierr

    tag = iDom + tagOffset
    sendCount = sendCount + 1
    call mpi_isend(oWall%rBuffer, size(oWall%rbuffer), adflow_real, &
         iProc, tag, ADflow_comm_world, sendRequests(sendCount), ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    sendCount = sendCount + 1
    call mpi_isend(oWall%iBuffer, size(oWall%iBuffer), adflow_integer, &
         iProc, tag, ADflow_comm_world, sendRequests(sendCount), ierr)
    call ECHK(ierr, __FILE__, __LINE__)

  end subroutine sendOSurf

  subroutine recvOBlock(oBlock, iDom, iProc, tagOffset, iSize, rSize, &
       recvCount, recvInfo)

    use constants
    use communication, only : adflow_comm_world, recvRequests
    use oversetData, only : oversetBlock
    use utils, only : EChk
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
    call mpi_irecv(oBlock%rBuffer, rSize, adflow_real, &
         iProc, tag, ADflow_comm_world, recvRequests(recvCount), ierr)
    call ECHK(ierr, __FILE__, __LINE__)
    recvInfo(:, recvCount) = (/iDom, 1/)

    recvCount = recvCount + 1
    call mpi_irecv(oBlock%iBuffer, iSize, adflow_integer, &
         iProc, tag, ADflow_comm_world, recvRequests(recvCount), ierr)
    call ECHK(ierr, __FILE__, __LINE__)
    recvInfo(:, recvCount) = (/iDom, 2/)

  end subroutine recvOBlock

  subroutine recvOFringe(oFringe, iDom, iProc, tagOffset, iSize, rSize, &
       recvCount, recvInfo)

    use constants
    use communication, only : adflow_comm_world, recvRequests
    use oversetData, only : oversetFringe
    use utils, only : EChk
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
    call mpi_irecv(oFringe%rBuffer, rSize, adflow_real, &
         iProc, tag, ADflow_comm_world, recvRequests(recvCount), ierr)
    call ECHK(ierr, __FILE__, __LINE__)
    recvInfo(:, recvCount) = (/iDom, 3/)

    recvCount = recvCount + 1
    call mpi_irecv(oFringe%iBuffer, iSize, adflow_integer, &
         iProc, tag, ADflow_comm_world, recvRequests(recvCount), ierr)
    call ECHK(ierr, __FILE__, __LINE__)
    recvInfo(:, recvCount) = (/iDom, 4/)

  end subroutine recvOFringe

  subroutine recvOSurf(oWall, iDom, iProc, tagOffset, iSize, rSize, &
       recvCount, recvInfo)

    use constants
    use communication, only : adflow_comm_world, recvRequests
    use oversetData, only : oversetWall
    use utils, only : EChk
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
    call mpi_irecv(oWall%rBuffer, rSize, adflow_real, &
         iProc, tag, ADflow_comm_world, recvRequests(recvCount), ierr)
    call ECHK(ierr, __FILE__, __LINE__)
    recvInfo(:, recvCount) = (/iDom, 5/)

    recvCount = recvCount + 1

    call mpi_irecv(oWall%iBuffer, iSize, adflow_integer, &
         iProc, tag, ADflow_comm_world, recvRequests(recvCount), ierr)
    call ECHK(ierr, __FILE__, __LINE__)
    recvInfo(:, recvCount) = (/iDom, 6/)

  end subroutine recvOSurf

  subroutine getFringeReturnSizes(oFringeSendList, oFringeRecvList, &
       nOFringeSend, nOfringeRecv, oFringes, &
       fringeRecvSizes, cumFringeRecv)

    ! For this data exchange we use the exact *reverse* of fringe
    ! communication pattern. This communiation simply determines the
    ! number of fringes that must be returned to the owning process.

    use constants
    use communication , only : sendRequests, recvRequests, adflow_comm_world
    use utils, only : EChk
    use oversetData, onlY : oversetFringe
    implicit none

    ! Input/output
    type(oversetFringe), dimension(:) :: oFringes
    integer(kind=intType), dimension(:, :) :: oFringeSendList, oFringeRecvList
    integer(kind=intType), dimension(:), allocatable :: cumFringeRecv, fringeRecvSizes
    integer(kind=intType) :: nOFringeSend, nOfringeRecv
    ! Working
    integer(kind=intType) :: sendCount, recvCount
    integer(kind=intType) :: iDom, iProc, jj, ierr, index, i
    integer mpiStatus(MPI_STATUS_SIZE)

    ! Post all the fringe iSends
    sendCount = 0
    do jj=1, nOFringeRecv

       iProc = oFringeRecvList(1, jj)
       iDom = oFringeRecvList(2, jj)
       sendCount = sendCount + 1
       call mpi_isend(oFringes(iDom)%fringeReturnSize, 1, adflow_integer, &
            iproc, iDom, adflow_comm_world, sendRequests(sendCount), ierr)
       call ECHK(ierr, __FILE__, __LINE__)
    end do

    allocate(fringeRecvSizes(nOfringeSend))

    ! Non-blocking receives
    recvCount = 0
    do jj=1, nOFringeSend

       iProc = oFringeSendList(1, jj)
       iDom = oFringeSendList(2, jj)
       recvCount = recvCount + 1

       call mpi_irecv(fringeRecvSizes(jj), 1, adflow_integer, &
            iProc, iDom, adflow_comm_world, recvRequests(recvCount), ierr)
       call ECHK(ierr, __FILE__, __LINE__)
    end do

    ! Last thing to do wait for all the sends and receives to finish
    do i=1,sendCount
       call mpi_waitany(sendCount, sendRequests, index, mpiStatus, ierr)
       call ECHK(ierr, __FILE__, __LINE__)
    end do

    do i=1,recvCount
       call mpi_waitany(recvCount, recvRequests, index, mpiStatus, ierr)
       call ECHK(ierr, __FILE__, __LINE__)
    end do

    ! Compute the cumulative form of the fringeRecvSizes

    allocate(cumFringeRecv(1:nOFringeSend+1))
    cumFringeRecv(1) = 1
    do jj=1, nOFringeSend ! These are the fringes we *sent*
       ! originally, now are going to receive them
       ! back
       cumFringeRecv(jj+1) = cumFringeRecv(jj) + fringeRecvSizes(jj)
    end do

  end subroutine getFringeReturnSizes



  !
  !       oversetLoadBalance determine the deistributation of donor and
  !       receiver blocks that will result in approximate even load
  !       balancing. The sparse matrix structrue of the overla is
  !       provided. This computation runs on all processors.

  subroutine oversetLoadBalance(overlap)

    use constants
    use communication, only : nProc
    use oversetData, only : CSRMatrix
    implicit none

    ! Input/Output
    type(CSRMatrix), intent(inout) :: overlap

    ! Working paramters
    integer(kind=intType) :: curRow, jj, jj1, iProc, iRow
    real(kind=realType) :: evenCost, potentialSum, targetCost
    real(Kind=realType) :: totalSearch, totalBuild

    real(kind=realType), dimension(0:nProc-1) :: procCosts
    real(kind=realType), dimension(0:nProc) :: cumProcCosts
    real(kind=realType), dimension(overlap%nRow) :: buildCost
    real(kind=realType), parameter :: tol=0.1_realType
    !  real(kind=realType), parameter :: K=10_realType
    logical, dimension(overlap%nnz) :: blockTaken
    logical :: increment

    ! Pointers to make code a litte easier to read
    integer(kind=intType), pointer, dimension(:) :: rowPtr, assignedProc
    real(kind=realType), pointer, dimension(:) :: data

    ! Set the couple  of pointers
    rowPtr => overlap%rowPtr
    assignedProc => overlap%assignedProc
    data => overlap%data

    ! Determine the total search cost:
    totalSearch = sum(overlap%data)

    ! Target amount of work for each processor
    evenCost = totalSearch / nProc

    ! Initialize the taken processor to False
    blockTaken = .False.

    ! Initialzie assignedProc to -1 since there could be entries we can
    ! ignore.
    assignedProc(:) = -1
    procCosts = zero
    cumProcCosts(0) = zero

    ! Initialize the starting point
    jj = 1
    iProc = 0

    ! Find the first row with non-zeros
    curRow = 1
    do while(rowPtr(curRow+1)-rowPtr(curRow) == 0)
       curRow = curRow + 1
    end do

    masterLoop: do while (curRow <= overlap%nRow .and. iProc <= nProc)

       ! Normally we increment
       increment = .True.

       ! This is our current target cost.
       targetCost = evenCost*(iProc + 1)

       ! It is still possible that data(jj) is zero. That's ok...we'll
       ! explictly ignore them.
       if (data(jj) /= zero .and. .not. (blockTaken(jj))) then

          if (procCosts(iProc) == 0 .or. iProc == nProc-1) then
             ! Must be added
             procCosts(iProc) = procCosts(iProc) + data(jj)
             blockTaken(jj) = .True.
             assignedProc(jj) = iProc

          else

             ! There is already something in there. See what the
             !  potential sum will be:
             potentialSum = cumProcCosts(iProc) + procCosts(iProc) + data(jj)

             if (potentialSum < targetCost - tol*evenCost) then
                !  We are not close to our limit yet so just add it normally
                procCosts(iProc) = procCosts(iProc) + data(jj)
                blockTaken(jj) = .True.
                assignedProc(jj) = iProc

             else if (potentialSum >= targetCost - tol*evenCost  .and. &
                  potentialSum <= targetCost + tol*evenCost) then

                ! This one looks perfect. Call it a day...add it and
                !  move on to the next proc

                procCosts(iProc) = procCosts(iProc) + data(jj)
                blockTaken(jj) = .True.
                assignedProc(jj) = iProc

                ! Processor can be incremented
                cumProcCosts(iProc+1) = cumProcCosts(iProc) + procCosts(iProc)
                iProc = iProc + 1
             else
                ! This means potentialSum > targetCost + tol*evenCost

                ! This is somewhat bad news...this may be *horrendly*
                ! load balanced. The algorithm dictates we *MUST*
                ! finish this proc no matter what before we go back to
                ! the outer loop. Essentially we know jj is bad,
                ! instead scan over the rest of the row and see if we
                ! can add something else that is decent.
                increment = .False.

                restOfRow: do jj1=jj+1, rowPtr(curRow+1)-1

                   potentialSum = cumProcCosts(iProc) + procCosts(iProc) + data(jj1)

                   if (data(jj1) /= zero .and. .not. (blockTaken(jj1))) then

                      if (potentialSum < targetCost - tol*evenCost) then
                         !Huh...that one fit in without going
                         ! over....add it and kep going in the loop

                         procCosts(iProc) = procCosts(iProc) + data(jj1)
                         blockTaken(jj1) = .True.
                         assignedProc(jj1) = iProc

                      else if (potentialSum >= targetCost - tol*evenCost  .and. &
                           potentialSum <= targetCost + tol*evenCost) then

                         ! This one fit in perfectly.
                         procCosts(iProc) = procCosts(iProc) + data(jj1)
                         blockTaken(jj1) = .True.
                         assignedProc(jj1) = iProc

                         ! No need to keep going
                         exit  restOfRow

                      end if
                   end if
                end do restOfRow

                ! Well, the loop finished, we may or may not have
                ! added something. If so great...if not, oh well. We
                ! just keep going to the next proc. That's the greedy
                ! algorithm for you.

                ! Processor can be incremented
                cumProcCosts(iProc+1) = cumProcCosts(iProc) + procCosts(iProc)
                iProc = iProc + 1
             end if
          end if
       end if

       ! Move 1 in jj, until we reach the end and wrap around.
       if (increment) then
          jj = jj + 1

          ! Switch to the next row:
          if (jj == rowPtr(curRow+1)) then

             ! This is really tricky...we know we're at the end of the
             ! row, but we have to SKIP OVER THE EMPTY rows, or else the
             ! algorithm will crap out. Keep incrementing the curRow
             ! until we get a row with something in it. Make sure we
             ! don't go out the end, so check again nRow

             findNextNonZeroRow: do while(jj == rowPtr(curRow+1))
                curRow = curRow + 1
                if (curRow > overlap%nRow) then
                   exit findNextNonZeroRow
                end if
             end do findNextNonZeroRow
          end if
       end if
    end do masterLoop

  end subroutine oversetLoadBalance

  subroutine exchangeFringes(level, sps, commPattern, internal)
    !
    !       ExchangeFringes exchanges the donorInformation of the fringes:
    !       donorProc, donorBlock, dIndex and donorFrac. It does this
    !       the 1:1 halos for the given level and spectral instance. Since
    !       we have real values and integer values we will do all the ints
    !       first and then the reals.
    !
    use constants
    use block, only : fringeType
    use blockPointers, only : flowDoms
    use communication, only : commType, internalCommType, recvBuffer, sendBuffer, myid, &
         adflow_comm_world, sendRequests, recvRequests
    use oversetUtilities, only : addToFringeList, windIndex
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps

    type(commType),          dimension(*), intent(in) :: commPattern
    type(internalCommType), dimension(*), intent(in) :: internal
    !
    !      Local variables.
    !
    integer :: size, procId, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: i, j, ii, jj, nVar, iFringe, jFringe
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2
    integer(kind=intType) :: il, jl, kl, myIndex
    type(fringeType) :: fringe
    integer(kind=intType), dimension(:), allocatable :: sendBufInt
    integer(kind=intType), dimension(:), allocatable :: recvBufInt

    ! Allocate the memory for the sending and receiving buffers.
    nVar = 3
    ii = commPattern(level)%nProcSend
    ii = commPattern(level)%nsendCum(ii)
    jj = commPattern(level)%nProcRecv
    jj = commPattern(level)%nrecvCum(jj)

    allocate(sendBufInt(ii*nVar), recvBufInt(jj*nVar), stat=ierr)

    ! Send the variables. The data is first copied into
    ! the send buffer after which the buffer is sent asap.

    ii = 1
    intSends: do i=1,commPattern(level)%nProcSend

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level)%sendProc(i)
       size   = nVar*commPattern(level)%nsend(i)

       ! Copy the data in the correct part of the send buffer.

       jj = ii
       do j=1,commPattern(level)%nsend(i)

          ! Store the block id and the indices of the donor
          ! a bit easier.

          d1 = commPattern(level)%sendList(i)%block(j)
          i1 = commPattern(level)%sendList(i)%indices(j,1)
          j1 = commPattern(level)%sendList(i)%indices(j,2)
          k1 = commPattern(level)%sendList(i)%indices(j,3)

          ! Copy integer values to buffer
          iFringe = flowDoms(d1, level, sps)%fringePtr(1, i1, j1, k1)
          if (iFringe > 0) then
             sendBufInt(jj  ) = flowDoms(d1,level,sps)%fringes(iFringe)%donorProc
             sendBufInt(jj+1) = flowDoms(d1,level,sps)%fringes(iFringe)%donorBlock
             sendBufInt(jj+2) = flowDoms(d1,level,sps)%fringes(iFringe)%dIndex
          else
             sendBufInt(jj  ) = -1
             sendBufInt(jj+1) = 0
             sendBufInt(jj+2) = 0
          end if

          jj = jj + nVar

       enddo

       ! Send the data.

       call mpi_isend(sendBufInt(ii), size, adflow_integer, procId, &
            procId, ADflow_comm_world, sendRequests(i),   &
            ierr)

       ! Set ii to jj for the next processor.

       ii = jj

    enddo intSends

    ! Post the nonblocking receives.

    ii = 1
    intReceives: do i=1,commPattern(level)%nProcRecv

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level)%recvProc(i)
       size    = nVar*commPattern(level)%nrecv(i)

       ! Post the receive.

       call mpi_irecv(recvBufInt(ii), size, adflow_integer, procId, &
            myId, ADflow_comm_world, recvRequests(i), ierr)

       ! And update ii.

       ii = ii + size

    enddo intReceives

    ! Copy the local data.

    intLocalCopy: do i=1,internal(level)%ncopy

       ! Store the block and the indices of the donor a bit easier.

       d1 = internal(level)%donorBlock(i)
       i1 = internal(level)%donorIndices(i,1)
       j1 = internal(level)%donorIndices(i,2)
       k1 = internal(level)%donorIndices(i,3)

       ! Idem for the halo's.

       d2 = internal(level)%haloBlock(i)
       i2 = internal(level)%haloIndices(i,1)
       j2 = internal(level)%haloIndices(i,2)
       k2 = internal(level)%haloIndices(i,3)


       iFringe = flowDoms(d1, level, sps)%fringePtr(1, i1, j1, k1)
       if (iFringe > 0) then
          ! The sender has an actual fringe. Nowe check if the
          ! receiver has somewhere already to put it:

          jFringe = flowDoms(d2, level, sps)%fringePtr(1, i2, j2, k2)

          ! Setup the new fringe:
          fringe%myBlock = d2

          il = flowDoms(d2, level, sps)%il
          jl = flowDoms(d2, level, sps)%jl
          kl = flowDoms(d2, level, sps)%kl
          fringe%myIndex = windIndex(i2, j2, k2, il, jl, kl)

          fringe%donorProc = flowDoms(d1, level, sps)%fringes(iFringe)%donorProc
          fringe%donorBlock= flowDoms(d1, level, sps)%fringes(iFringe)%donorBlock
          fringe%dIndex    = flowDoms(d1, level, sps)%fringes(iFringe)%dIndex
          fringe%donorFrac = flowDoms(d1, level, sps)%fringes(iFringe)%donorFrac

          if (jFringe > 0) then
             ! Just copy the fringe into the slot. No need to update
             ! the pointer since it is already correct.
             flowDoms(d2, level, sps)%fringes(jFringe) = fringe
          else

             ! There is no slot available yet. Tack the fringe onto
             ! the end of the d2 fringe list and set the pointers
             ! accordingly.

             call addToFringeList(flowDoms(d2, level, sps)%fringes, &
                  flowDoms(d2, level, sps)%nDonors,  fringe)

             ! Note that all three values (pointer, start and end) are
             ! all the same.
             flowDoms(d2, level, sps)%fringePtr(:, i2, j2, k2) = &
                  flowDoms(d2, level, sps)%nDonors
          end if
       else
          ! The donor isn't a receiver so make sure the halo isn't
          ! either. Just set the fringePtr to 0

          flowDoms(d2, level, sps)%fringePtr(1, i2, j2, k2) = 0
       end if

    enddo intLocalCopy

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern(level)%nProcRecv
    intCompleteRecvs: do i=1,commPattern(level)%nProcRecv

       ! Complete any of the requests.

       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)

       ! Copy the data just arrived in the halo's.

       ii = index
       jj = nVar*commPattern(level)%nrecvCum(ii-1)
       do j=1,commPattern(level)%nrecv(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern(level)%recvList(ii)%block(j)
          i2 = commPattern(level)%recvList(ii)%indices(j,1)
          j2 = commPattern(level)%recvList(ii)%indices(j,2)
          k2 = commPattern(level)%recvList(ii)%indices(j,3)

          fringe%myBlock = d2
          ! Recompute my Index:
          il = flowDoms(d2, level, sps)%il
          jl = flowDoms(d2, level, sps)%jl
          kl = flowDoms(d2, level, sps)%kl
          fringe%myIndex = windIndex(i2, j2, k2, il, jl, kl)

          fringe%donorProc  = recvBufInt(jj+1)
          fringe%donorBlock = recvBufInt(jj+2)
          fringe%dIndex     = recvBufInt(jj+3)

          iFringe = flowDoms(d2, level, sps)%fringePtr(1, i2, j2, k2)
          if (iFringe > 0) then
             ! We have somehwere to to put the data already:
             flowDoms(d2, level, sps)%fringes(iFringe) = fringe
          else
             ! We don't somehwhere to put the fringe to add to the list:
             call addToFringeList(flowDoms(d2, level, sps)%fringes, &
                  flowDoms(d2, level, sps)%nDonors, fringe)

             ! Note that all three values (pointer, start and end) are
             ! all the same.
             flowDoms(d2, level, sps)%fringePtr(:, i2, j2, k2) = &
                  flowDoms(d2, level, sps)%nDonors
          end if
          jj = jj + nVar
       enddo

    enddo intCompleteRecvs

    ! Complete the nonblocking sends.

    size = commPattern(level)%nProcSend
    do i=1,commPattern(level)%nProcSend
       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
    enddo

    ! Done with the integer memory.

    deallocate(sendBufInt, recvBufInt)

    ! Now do the real exchange. We can use the regular real buffers here
    ! since they are large enough


    ! ================================================================================


    ! Allocate the memory for the sending and receiving buffers.
    nVar = 3

    ! Send the variables. The data is first copied into
    ! the send buffer after which the buffer is sent asap.

    ii = 1
    sends: do i=1,commPattern(level)%nProcSend

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level)%sendProc(i)
       size   = nVar*commPattern(level)%nsend(i)

       ! Copy the data in the correct part of the send buffer.

       jj = ii
       do j=1,commPattern(level)%nsend(i)

          ! Store the block id and the indices of the donor
          ! a bit easier.

          d1 = commPattern(level)%sendList(i)%block(j)
          i1 = commPattern(level)%sendList(i)%indices(j,1)
          j1 = commPattern(level)%sendList(i)%indices(j,2)
          k1 = commPattern(level)%sendList(i)%indices(j,3)

          ! Copy real values to buffer
          iFringe = flowDoms(d1, level, sps)%fringePtr(1, i1, j1, k1)
          if (iFringe > 0) then
             sendBuffer(jj:jj+2) = flowDoms(d1,level,sps)%fringes(iFringe)%donorFrac
          else
             sendBuffer(jj:jj+2) = zero
          end if
          jj = jj + nVar

       enddo

       ! Send the data.

       call mpi_isend(sendBuffer(ii), size, adflow_real, procId, &
            procId, ADflow_comm_world, sendRequests(i),   &
            ierr)

       ! Set ii to jj for the next processor.

       ii = jj

    enddo sends

    ! Post the nonblocking receives.

    ii = 1
    receives: do i=1,commPattern(level)%nProcRecv

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level)%recvProc(i)
       size    = nVar*commPattern(level)%nrecv(i)

       ! Post the receive.

       call mpi_irecv(recvBuffer(ii), size, adflow_real, procId, &
            myId, ADflow_comm_world, recvRequests(i), ierr)

       ! And update ii.

       ii = ii + size

    enddo receives

    ! ***********************************************************
    ! No local copy since we copied the fringes directly and the
    ! donorFrac info is already there.
    ! ***********************************************************

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern(level)%nProcRecv
    completeRecvs: do i=1,commPattern(level)%nProcRecv

       ! Complete any of the requests.

       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)

       ! Copy the data just arrived in the halo's.

       ii = index
       jj = nVar*commPattern(level)%nrecvCum(ii-1)
       do j=1,commPattern(level)%nrecv(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern(level)%recvList(ii)%block(j)
          i2 = commPattern(level)%recvList(ii)%indices(j,1)
          j2 = commPattern(level)%recvList(ii)%indices(j,2)
          k2 = commPattern(level)%recvList(ii)%indices(j,3)

          ! Now, there should already be a spot available for the
          ! donorFrac since it was created if necessary in the integer exchange.
          iFringe = flowDoms(d2, level, sps)%fringePtr(1, i2, j2, k2)
          flowDoms(d2, level, sps)%fringes(iFringe)%donorFrac = recvBuffer(jj+1:jj+3)

          jj = jj + nVar
       enddo

    enddo completeRecvs

    ! Complete the nonblocking sends.

    size = commPattern(level)%nProcSend
    do i=1,commPattern(level)%nProcSend
       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
    enddo

  end subroutine exchangeFringes

  subroutine exchangeStatus(level, sps, commPattern, internal)
    !
    !       ExchangeIsCompute exchanges the isCompute flag for the 1 to 1
    !       connectivity for the given level and sps instance.
    !
    use constants
    use blockPointers, only : nDom, ib, jb, kb, flowDoms
    use communication, only : commType, internalCommType
    use utils, only : setPointers
    use haloExchange, only : whalo1to1intgeneric

    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps

    type(commType),          dimension(*), intent(in) :: commPattern
    type(internalCommType), dimension(*), intent(in) :: internal
    integer(kind=intType) :: nn

    domainLoop:do nn=1, nDom
       flowDoms(nn, level, sps)%intCommVars(1)%var => flowDoms(nn, level, sps)%status(:, :, :)
    end do domainLoop

    ! Run the generic integer exchange
    call wHalo1to1IntGeneric(1, level, sps, commPattern, internal)

  end subroutine exchangeStatus

  subroutine exchangeStatusTranspose(level, sps, commPattern, internal)

    ! exchangeStatusTranspose performs the *TRANSPOSE* of the normal
    ! halo exchange. That means it takes information *in the halo cells*
    ! and accumulate it into the *owned cells*. In this particular case,
    ! we are transmitting the isDonor and isWallDonor information from
    ! the halos to the owned cells. The "accumulate" operation will be
    ! an MPI_LOR.  Note that this actually hast he same comm structure
    ! as 'whalo1to1_b'.

    use constants
    use blockPointers, only : flowDoms
    use communication, only : commType, internalCommType, recvBuffer, sendBuffer, myid, &
         adflow_comm_world, sendRequests, recvRequests
    use oversetUtilities, only : getStatus, setIsDonor, setIsWallDonor
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps
    type(commType), dimension(*), intent(in)         :: commPattern
    type(internalCommType), dimension(*), intent(in) :: internal
    !
    !      Local variables.
    !
    integer :: size, procID, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: mm
    integer(kind=intType) :: i, j, k, ii, jj
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2
    integer(kind=intType), dimension(:), allocatable :: sendBuf, recvBuf
    logical :: CisDonor, CisHole, CisCompute, CisFloodSeed, CisFlooded, CisWallDonor, CisReceiver
    logical :: DisDonor, DisHole, DisCompute, DisFloodSeed, DisFlooded, DisWallDonor, DisReceiver
    integer(kind=intType) :: cellStatus, donorStatus



    ii = commPattern(level)%nProcSend
    ii = commPattern(level)%nsendCum(ii)
    jj = commPattern(level)%nProcRecv
    jj = commPattern(level)%nrecvCum(jj)

    ! We are exchanging 1 piece of information
    allocate(sendBuf(ii), recvBuf(jj), stat=ierr)

    ! Gather up the seeds into the *recv* buffer. Note we loop
    ! over nProcRECV here! After the buffer is assembled it is
    ! sent off.

    jj = 1
    ii = 1
    recvs: do i=1,commPattern(level)%nProcRecv

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level)%recvProc(i)
       size    = commPattern(level)%nrecv(i)

       ! Copy the data into the buffer

       do j=1,commPattern(level)%nrecv(i)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern(level)%recvList(i)%block(j)
          i2 = commPattern(level)%recvList(i)%indices(j,1)
          j2 = commPattern(level)%recvList(i)%indices(j,2)
          k2 = commPattern(level)%recvList(i)%indices(j,3)

          recvBuf(jj) = flowDoms(d2, level, sps)%status(i2, j2, k2)
          jj = jj + 1

       enddo

       ! Send the data.
       call mpi_isend(recvBuf(ii), size, adflow_integer, procID,  &
            procID, ADflow_comm_world, sendRequests(i), &
            ierr)

       ! Set ii to jj for the next processor.

       ii = jj

    enddo recvs

    ! Post the nonblocking receives.

    ii = 1
    sends: do i=1,commPattern(level)%nProcSend

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level)%sendProc(i)
       size    = commPattern(level)%nsend(i)

       ! Post the receive.

       call mpi_irecv(sendBuf(ii), size, adflow_integer, procID, &
            myID, ADflow_comm_world, recvRequests(i), ierr)

       ! And update ii.

       ii = ii + size

    enddo sends

    ! Copy the local data.

    localCopy: do i=1,internal(level)%ncopy

       ! Store the block and the indices of the donor a bit easier.

       d1 = internal(level)%donorBlock(i)
       i1 = internal(level)%donorIndices(i,1)
       j1 = internal(level)%donorIndices(i,2)
       k1 = internal(level)%donorIndices(i,3)

       ! Idem for the halo's.

       d2 = internal(level)%haloBlock(i)
       i2 = internal(level)%haloIndices(i,1)
       j2 = internal(level)%haloIndices(i,2)
       k2 = internal(level)%haloIndices(i,3)

       ! OR operation. Note we modify the '1' values ie. the 'donors'
       ! which are now receivers because of the transpose operation.
       cellStatus = flowDoms(d1, level, sps)%status(i1, j1, k1)
       call getStatus(cellStatus, CisDonor, CisHole, CisCompute, &
            CisFloodSeed, CisFlooded, CisWallDonor, CisReceiver)

       donorStatus = flowDoms(d2, level, sps)%status(i2, j2, k2)
       call getStatus(donorStatus, DisDonor, DisHole, DisCompute, &
            DisFloodSeed, DisFlooded, DisWallDonor, DisReceiver)

       call setIsDonor(flowDoms(d1, level, sps)%status(i1, j1, k1), &
            CIsDonor .or. DisDonor)

       call setIsWallDonor(flowDoms(d1, level, sps)%status(i1, j1, k1), &
            CIsWallDonor .or. DisWallDonor)

    enddo localCopy

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern(level)%nProcSend
    completeSends: do i=1,commPattern(level)%nProcSend

       ! Complete any of the requests.

       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)

       ! Copy the data just arrived in the halo's.

       ii = index

       jj = commPattern(level)%nsendCum(ii-1)

       do j=1,commPattern(level)%nsend(ii)

          ! Store the block and the indices of the halo a bit easier.

          d1 = commPattern(level)%sendList(ii)%block(j)
          i1 = commPattern(level)%sendList(ii)%indices(j,1)
          j1 = commPattern(level)%sendList(ii)%indices(j,2)
          k1 = commPattern(level)%sendList(ii)%indices(j,3)


          cellStatus = flowDoms(d1, level, sps)%status(i1, j1, k1)
          call getStatus(cellStatus, CisDonor, CisHole, CisCompute, &
               CisFloodSeed, CisFlooded, CisWallDonor, CisReceiver)
          jj = jj + 1
          donorStatus = sendBuf(jj)
          call getStatus(donorStatus, DisDonor, DisHole, DisCompute, &
               DisFloodSeed, DisFlooded, DisWallDonor, DisReceiver)

          call setIsDonor(flowDoms(d1, level, sps)%status(i1, j1, k1), &
               CIsDonor .or. DisDonor)

          call setIsWallDonor(flowDoms(d1, level, sps)%status(i1, j1, k1), &
               CIsWallDonor .or. DisWallDonor)
       enddo
    enddo completeSends

    ! Complete the nonblocking sends.

    size = commPattern(level)%nProcRecv
    do i=1,commPattern(level)%nProcRecv
       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
    enddo

    deallocate(recvBuf, sendBuf)

  end subroutine exchangeStatusTranspose

  subroutine setupFringeGlobalInd(level, sps)

    use constants
    use blockPointers
    use communication
    use utils, only : EChk
    implicit none

    ! This subroutine is used to record the global index of each of
    ! the donors for overset fringes. It has the same comm structure
    ! as  wOverset and flagInvalidDonors.

    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps

    !
    !      Local variables.
    !
    integer :: size, procId, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: nVar
    integer(kind=intType) :: i, j, k, ii, jj, iii, jjj, kkk, iFringe
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2, ind
    integer(kind=intType), dimension(:), allocatable :: sendBufInt
    integer(kind=intType), dimension(:), allocatable :: recvBufInt
    logical :: invalid
    type(commType), pointer :: commPattern
    type(internalCommType), pointer :: internal

    commPattern => commPatternOverset(level, sps)
    internal => internalOverset(level, sps)

    ii = commPattern%nProcSend
    ii = commPattern%nsendCum(ii)
    jj = commPattern%nProcRecv
    jj = commPattern%nrecvCum(jj)
    nVar = 8
    allocate(sendBufInt(ii*nVar), recvBufInt(jj*nVar), stat=ierr)

    ! Send the variables. The data is first copied into
    ! the send buffer after which the buffer is sent asap.

    ii = 1
    sends: do i=1,commPattern%nProcSend

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern%sendProc(i)
       size    = nVar*commPattern%nsend(i)

       ! Copy the data in the correct part of the send buffer.

       jj = ii
       do j=1,commPattern%nsend(i)

          ! Store the block id and the indices of the donor
          ! a bit easier.

          d1 = commPattern%sendList(i)%block(j)
          i1 = commPattern%sendList(i)%indices(j,1)
          j1 = commPattern%sendList(i)%indices(j,2)
          k1 = commPattern%sendList(i)%indices(j,3)

          ! Loop over the 8 donors:
          do kkk=k1, k1+1
             do jjj=j1, j1+1
                do iii=i1, i1+1
                   sendBufInt(jj) = flowDoms(d1, level, sps)%globalCell(iii,jjj,kkk)
                   jj =jj + 1
                end do
             end do
          end do
       enddo

       ! Send the data.

       call mpi_isend(sendBufInt(ii), size, adflow_integer, procId,  &
            procId, ADflow_comm_world, sendRequests(i), &
            ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Set ii to jj for the next processor.

       ii = jj

    enddo sends

    ! Post the nonblocking receives.

    ii = 1
    receives: do i=1,commPattern%nProcRecv

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern%recvProc(i)
       size    = nVar*commPattern%nrecv(i)

       ! Post the receive.

       call mpi_irecv(recvBufInt(ii), size, adflow_integer, procId, &
            myId, ADflow_comm_world, recvRequests(i), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! And update ii.

       ii = ii + size

    enddo receives

    ! Do the local interpolation.
    localInterp: do i=1,internal%ncopy

       ! Store the block and the indices of the donor a bit easier.

       d1 = internal%donorBlock(i)
       i1 = internal%donorIndices(i, 1)
       j1 = internal%donorIndices(i, 2)
       k1 = internal%donorIndices(i, 3)

       ! Idem for the halo's.

       d2 = internal%haloBlock(i)
       i2 = internal%haloIndices(i, 1)
       j2 = internal%haloIndices(i, 2)
       k2 = internal%haloIndices(i, 3)

       ! Loop over the 8 donors:
       ind = 0
       do kkk=k1, k1+1
          do jjj=j1, j1+1
             do iii=i1, i1+1
                ind = ind + 1
                flowDoms(d2, level, sps)%gInd(ind, i2, j2, k2) = &
                     flowDoms(d1, level, sps)%globalCell(iii,jjj,kkk)
             end do
          end do
       end do
    enddo localInterp

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern%nProcRecv
    completeRecvs: do i=1,commPattern%nProcRecv

       ! Complete any of the requests.

       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Copy the data just arrived in the halo's.

       ii = index
       jj = nVar*commPattern%nrecvCum(ii-1)
       do j=1,commPattern%nrecv(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern%recvList(ii)%block(j)
          i2 = commPattern%recvList(ii)%indices(j,1)
          j2 = commPattern%recvList(ii)%indices(j,2)
          k2 = commPattern%recvList(ii)%indices(j,3)

          ! Just set the 8 values
          do ind=1,8
             flowDoms(d2, level, sps)%gInd(ind, i2, j2, k2) = &
                  recvBufInt(jj+ind)
          end do

          jj = jj + 8
       enddo
    end do completeRecvs

    ! Complete the nonblocking sends.

    size = commPattern%nProcSend
    do i=1,commPattern%nProcSend
       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    enddo
    deallocate(sendBufInt, recvBufInt)

  end subroutine setupFringeGlobalInd

  subroutine exchangeSurfaceDelta(zipperFamList, level, sps, commPattern, internal)
    !
    !       ExchangeSurfaceDelta exchanges surface delta to fill up halo
    !       surface cells from adjacent blocks.
    !
    use constants
    use blockPointers, onlY : nDom, flowDoms, nBocos, BCType, BCFaceID, BCData, &
         ib, il, jb, jl, kb, kl
    use communication, only : commType, internalCOmmType, commPatternCell_1st, &
         internalCell_1st
    use utils, only : setPointers
    use haloExchange, only : whalo1to1RealGeneric
    use sorting, only : famInList
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in), dimension(:) :: zipperFamList
    integer(kind=intType), intent(in) :: level, sps

    type(commType),          dimension(*), intent(in) :: commPattern
    type(internalCommType), dimension(*), intent(in) :: internal

    ! Local
    integer(kind=intType) :: i, j, k, ii, nn, mm
    real(kind=realType), dimension(:), allocatable :: pSave
    real(kind=realType), dimension(:, :), pointer :: deltaPtr

    ! Just cheat by exchangint pressure. saving Pressure, dumping deltaPtr into the pressure,
    ! exchanging that and then restoring the pressure

    do nn=1, nDom
       call setPointers(nn, level, sps)

       ! Allocate pointer space for the integer flag communication
       allocate(flowDoms(nn, level, sps)%realCommVars(1)%var(1:ib+1, 1:jb+1, 1:kb+1))

       ! Push the surface iblank back to the generic volume variable rVar1
       bocoLoop: do mm=1, nBocos
          famInclude: if (famInList(BCData(mm)%famID, zipperFamList)) then

             select case (BCFaceID(mm))
             case (iMin)
                deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(2+1, :, :)
             case (iMax)
                deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(il+1, :, :)
             case (jMin)
                deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, 2+1,  :)
             case (jMax)
                deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, jl+1, :)
             case (kMin)
                deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, :, 2+1 )
             case (kMax)
                deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, :, kl+1)
             end select

             ! NO HALOS!
             do j=BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd

                   ! Remember to account for the pointer offset since
                   ! the iblank starts at zero
                   deltaPtr(i+1, j+1) = BCData(mm)%delta(i, j)
                end do
             end do
          end if famInclude
       end do bocoLoop
    end do

    ! Exchange the variane
    call whalo1to1RealGeneric(1, level, sps, commPatternCell_1st, internalCell_1st)

    ! Copy back out
    ii = 0
    do nn=1, nDom
       call setPointers(nn, level, sps)

       ! Extract the surface iblank from the volume.
       bocoLoop2: do mm=1, nBocos
          famInclude2: if (famInList(BCData(mm)%famID, zipperFamList)) then

             select case (BCFaceID(mm))
             case (iMin)
                deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(2+1, :, :)
             case (iMax)
                deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(il+1, :, :)
             case (jMin)
                deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, 2+1,  :)
             case (jMax)
                deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, jl+1, :)
             case (kMin)
                deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, :, 2+1 )
             case (kMax)
                deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, :, jl+1)
             end select

             ! INCLUDE THE HALOS!
             do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                   ! Remember to account for the pointer offset since
                   ! the iblank starts at zero
                   BCData(mm)%delta(i,j) = deltaPtr(i+1, j+1)
                end do
             end do
          end if famInclude2
       end do bocoLoop2

       ! Now deallocate this pointer
       deallocate(flowDoms(nn, level, sps)%realCommVars(1)%var)
    end do
  end subroutine exchangeSurfaceDelta

  subroutine exchangeSurfaceIblanks(zipperFamList, level, sps, commPattern, internal)
    !
    !       ExchangeIblank exchanges the 1 to 1 internal halo's for the
    !       given level and sps instance.
    !
    use constants
    use blockPointers, only : nDom, ib, jb, kb, iBlank, BCData, &
         il, jl, kl, BCFaceID, nBocos, flowDoms, BCType
    use communication, only : commType, internalCommType
    use utils, only : setPointers
    use haloExchange, only : whalo1to1intgeneric
    use sorting, only : famInList
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps
    integer(kind=intType), intent(in), dimension(:) :: zipperFamList

    type(commType),          dimension(*), intent(in) :: commPattern
    type(internalCommType), dimension(*), intent(in) :: internal

    ! Local
    integer(kind=intType) :: i, j, k, ii, nn, mm
    integer(kind=intType), dimension(:), allocatable :: iBlankSave
    integer(kind=intType), dimension(:, :), pointer :: ibp

    ! Just cheat by saving iBlank iblank array, resusing itand
    ii = 0
    do nn=1, nDom
       call setPointers(nn, level, sps)
       ii = ii + (ib+1)*(jb+1)*(kb+1)
    end do

    allocate(iBlankSave(ii))
    ii = 0
    do nn=1, nDom
       call setPointers(nn, level, sps)
       do k=0,kb
          do j=0,jb
             do i=0,ib
                ii =ii + 1
                iBlankSave(ii) = iblank(i,j,k)
                ! The following algorithm uses the volume iblank array to communicate surface blocking
                ! across block boundaries. However, for h-topology meshes, for example at a sharp
                ! trailing edge this breaks, because the surface on the top an bottom of the shape are
                ! not topologically adjacent in the block structure. If we leave the existing volume
                ! blanking in place the correct values are communicated, if we zero it as done originally
                ! the connection is broken. Therefore the following line is commented out.
                !iblank(i,j,k) = 0 !commented out to fix issue with h-topology blocks on the zipper
                
             end do
          end do
       end do

       ! Push the surface iblank back to the volume:
       bocoLoop: do mm=1, nBocos
          famInclude: if (famInList(BCData(mm)%famID, zipperFamList)) then

             select case (BCFaceID(mm))
             case (iMin)
                ibp => iblank(2, :, :)
             case (iMax)
                ibp => iblank(il, :, :)
             case (jMin)
                ibp => iblank(:, 2, :)
             case (jMax)
                ibp => iblank(:, jl, :)
             case (kMin)
                ibp => iblank(:, :, 2)
             case (kMax)
                ibp => iblank(:, :, kl)
             end select

             ! NO HALOS!
             do j=BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd

                   ! Remember to account for the pointer offset since
                   ! the iblank starts at zero
                   ibp(i+1, j+1) = BCData(mm)%iBlank(i,j)
                end do
             end do
          end if famInclude
       end do bocoLoop
    end do

    ! Exchange iblanks
    domainLoop:do nn=1, nDom
       flowDoms(nn, level, sps)%intCommVars(1)%var => &
            flowDoms(nn, level, sps)%iblank(:, :, :)
    end do domainLoop

    ! Run the generic integer exchange
    call wHalo1to1IntGeneric(1, level, sps, commPattern, internal)

    ii = 0
    do nn=1, nDom
       call setPointers(nn, level, sps)

       ! Extract the surface iblank from the volume.
       bocoLoop2: do mm=1, nBocos
          famInclude2: if (famInList(BCData(mm)%famID, zipperFamList)) then

             select case (BCFaceID(mm))
             case (iMin)
                ibp => iblank(2, :, :)
             case (iMax)
                ibp => iblank(il, :, :)
             case (jMin)
                ibp => iblank(:, 2, :)
             case (jMax)
                ibp => iblank(:, jl, :)
             case (kMin)
                ibp => iblank(:, :, 2)
             case (kMax)
                ibp => iblank(:, :, kl)
             end select

             ! INCLUDE THE HALOS!
             do j=BCData(mm)%jnBeg, BCData(mm)%jnEnd+1
                do i=BCData(mm)%inBeg, BCData(mm)%inEnd+1
                   ! Remember to account for the pointer offset since
                   ! the iblank starts at zero
                   BCData(mm)%iBlank(i,j) = ibp(i+1, j+1)
                end do
             end do
          end if famInclude2
       end do bocoLoop2

       ! Restore the saved array
       do k=0,kb
          do j=0,jb
             do i=0,ib
                ii =ii + 1
                iBlank(i,j,k) = iBlankSave(ii)
             end do
          end do
       end do
    end do
    deallocate(iblankSave)
  end subroutine exchangeSurfaceIblanks

  subroutine emptyOversetComm(level, sps)

    !  Short cut function to make empty overset comm structure for
    !  problems that do not use overset meshes.

    use constants
    use communication, only : commPatternOverset, internalOverset
    implicit none

    ! Function
    integer(kind=intType), intent(in) :: level, sps

    ! Working
    integer(Kind=intType) :: nn, mm, ierr

    commPatternOverset(level, sps)%nProcRecv = 0
    allocate(commPatternOverset(level, sps)%recvProc(0))
    allocate(commPatternOverset(level, sps)%nRecv(0))
    allocate(commPatternOverset(level, sps)%recvList(0))

    commPatternOverset(level, sps)%nProcSend = 0
    allocate(commPatternOverset(level, sps)%sendProc(0))
    allocate(commPatternOverset(level, sps)%nSend(0))
    allocate(commPatternOverset(level, sps)%sendList(0))

    internalOverset(level, sps)%nCopy = 0
    allocate(internalOverset(level, sps)%donorBlock(0))
    allocate(internalOverset(level, sps)%donorIndices(0, 3))
    allocate(internalOverset(level, sps)%donorInterp(0, 8))
    allocate(internalOverset(level, sps)%haloBlock(0))
    allocate(internalOverset(level, sps)%haloIndices(0, 3))

  end subroutine emptyOversetComm

  subroutine updateOversetConnectivity(level, sps)

    ! This subroutine updates the overset connectivity for a perturbed
    ! mesh. It does *not* completely redo the connectivity. Rather, a
    ! newton search on the existing donors are performed using the
    ! updated coordinates. This type of update is only applicable if the
    ! entire volume mesh is warped as one a-la pyWarpUstruct. This
    ! actually ends up being a fairly small correction most of the time,
    ! however, desipite looks to the contrary is actually quite fast to
    ! run.

    use constants
    use communication
    use blockPointers, only : nDom, il, jl, kl, xSeed, flowDoms, x, ib, jb, kb, &
         ie, je, ke, fringes,scratch
    use haloExchange, only : whalo1to1RealGeneric
    use oversetUtilities, only : newtonUpdate, fracToWeights
    use utils, only : setPointers

    implicit none

    ! Input
    integer(kind=intType), intent(in) :: level, sps
    type(commType), pointer :: commPattern
    type(internalCommType), pointer :: internal

    ! Working
    integer(kind=intType) :: nn, ii,jj, ierr,  i, j, k, d1, i1, j1, k1, d2, i2, j2, k2
    integer(kind=intType) :: size, procID, index, iii,jjj
    integer, dimension(mpi_status_size) :: mpiStatus
    real(kind=realType) :: frac(3), frac0(3), xCen(3)
    integer(kind=intType), dimension(8), parameter :: indices=(/1,2,4,3,5,6,8,7/)

    ! Pointers to the overset comms to make it easier to read
    commPattern => commPatternOverset(level, sps)
    internal => internalOverset(level, sps)

    ! Step 1: Since we need to update donors for all cells including the
    ! donors for double halos, we must know the new cell center
    ! locations for the all of these receivers. Unfortunately, the
    ! double halos don't have coordinates, so we must first perform a
    ! (forward) block-to-block halo exchange to populate the xSeed
    ! values for all cells, including double halos.

    do nn=1, nDom
       call setPointers(nn, level, sps)

       if (.not. associated(flowDoms(nn, level, sps)%xSeed)) then
          allocate(flowDoms(nn, level, sps)%XSeed(0:ib, 0:jb, 0:kb, 3))
       end if
       xSeed => flowDoms(nn, level, sps)%xSeed
       xSeed = zero
       do k=2,kl
          do j=2, jl
             do i=2, il
                xSeed(i, j, k, :) = eighth*(&
                     x(i-1, j-1, k-1, :) + &
                     x(i  , j-1, k-1, :) + &
                     x(i-1, j  , k-1, :) + &
                     x(i  , j  , k-1, :) + &
                     x(i-1, j-1, k  , :) + &
                     x(i  , j-1, k  , :) + &
                     x(i-1, j  , k  , :) + &
                     x(i  , j  , k  , :)) !+ fringes(i,j,k)%offset
             end do
          end do
       end do
    end do

    ! Exchange the xSeeds.
    do nn=1, nDom
       flowDoms(nn, level, sps)%realCommVars(1)%var => flowDoms(nn, level, sps)%xSeed(:, :, :, 1)
       flowDoms(nn, level, sps)%realCommVars(2)%var => flowDoms(nn, level, sps)%xSeed(:, :, :, 2)
       flowDoms(nn, level, sps)%realCommVars(3)%var => flowDoms(nn, level, sps)%xSeed(:, :, :, 3)
    end do

    ! Run the (foward) generic halo exchange.
    call wHalo1to1RealGeneric(3, level, sps, commPatternCell_2nd, internalCell_2nd)

    ! Step 2: Next we need to communicate the xSeeds to their donor
    ! procs. This means running the overset exchange in REVERSE (ie from
    ! receiver to donor).  Most of this code will look like
    ! wOverset_b. We will runt he newtonUpdate code (below) on the fly
    ! as we receive the data, which should hide some of the comm time.

    ! Gather up the seeds into the *recv* buffer. Note we loop over
    ! nProcRECV here! After the buffer is assembled it is send off.

    jj = 1
    ii = 1
    recvs: do i=1,commPattern%nProcRecv

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern%recvProc(i)
       size    = 3*commPattern%nrecv(i)

       ! Copy the data into the buffer

       do j=1,commPattern%nrecv(i)

          ! Store the block and the indices to make code a bit easier to read

          d2 = commPattern%recvList(i)%block(j)
          i2 = commPattern%recvList(i)%indices(j,1)
          j2 = commPattern%recvList(i)%indices(j,2)
          k2 = commPattern%recvList(i)%indices(j,3)

          ! Copy the xSeed
          recvBuffer(jj)   = flowDoms(d2,level,sps)%xSeed(i2,j2,k2,1)
          recvBuffer(jj+1) = flowDoms(d2,level,sps)%xSeed(i2,j2,k2,2)
          recvBuffer(jj+2) = flowDoms(d2,level,sps)%xSeed(i2,j2,k2,3)
          jj = jj + 3
       end do

       ! Send the data.
       call mpi_isend(recvBuffer(ii), size, adflow_real, procID,  &
            procID, ADflow_comm_world, recvRequests(i), &
            ierr)

       ! Set ii to jj for the next processor.

       ii = jj

    end do recvs

    ! Post the nonblocking receives.

    ii = 1
    sends: do i=1,commPattern%nProcSend

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern%sendProc(i)
       size    = 3*commPattern%nsend(i)

       ! Post the receive.

       call mpi_irecv(sendBuffer(ii), size, adflow_real, procId, &
            myId, ADflow_comm_world, sendRequests(i), ierr)

       ! And update ii.

       ii = ii + size

    enddo sends

    ! Do the local interpolation.

    localInterp: do i=1,internal%ncopy

       ! Store the block and the indices of the donor a bit easier.

       d1 = internal%donorBlock(i)
       i1 = internal%donorIndices(i, 1)
       j1 = internal%donorIndices(i, 2)
       k1 = internal%donorIndices(i, 3)

       ! Idem for the halo's.

       d2 = internal%haloBlock(i)
       i2 = internal%haloIndices(i, 1)
       j2 = internal%haloIndices(i, 2)
       k2 = internal%haloIndices(i, 3)

       ! xCen is the '2'. This was the receiver, but since this is
       ! reverse, it's now the "input"
       xCen = flowDoms(d2, level, sps)%xSeed(i2, j2, k2, :)

       ! Store in the comm structure. We need it for the derivatives.
       internal%xCen(i, :) = xCen

       ! Do newton update
       frac0 = (/half, half, half/)
       call newtonUpdate(xCen, &
            flowDoms(d1, level, sps)%x(i1-1:i1+1, j1-1:j1+1, k1-1:k1+1, :), frac0, frac)

       ! Set the new weights
       call fracToWeights(frac, internal%donorInterp(i, :))
    enddo localInterp

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern%nProcSend
    completeSends: do i=1,commPattern%nProcSend

       ! Complete any of the requests.

       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)


       ii = index

       jj = 3*commPattern%nsendCum(ii-1)
       do j=1,commPattern%nsend(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern%sendList(ii)%block(j)
          i2 = commPattern%sendList(ii)%indices(j,1)
          j2 = commPattern%sendList(ii)%indices(j,2)
          k2 = commPattern%sendList(ii)%indices(j,3)

          xCen = sendBuffer(jj+1:jj+3)
          jj = jj + 3

          ! Store in the comm structure. We need it for derivatives.
          commPattern%sendList(ii)%xCen(j, :) = xCen

          ! Compute new fraction
          frac0 = (/half, half, half/)
          call newtonUpdate(xCen, &
               flowDoms(d2, level, sps)%x(i2-1:i2+1, j2-1:j2+1, k2-1:k2+1, :), frac0, frac)

          ! Set the new weights
          call fracToWeights(frac, commPattern%sendList(ii)%interp(j, :))
       enddo

    enddo completeSends

    ! Complete the nonblocking sends.

    size = commPattern%nProcRecv
    do i=1,commPattern%nProcRecv
       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)
    enddo

  end subroutine updateOversetConnectivity
#ifndef USE_COMPLEX
  subroutine updateOversetConnectivity_d(level, sps)

    ! Forward mode linearization of updateOversetConnectivity

    use constants
    use communication
    use blockPointers, only : nDom, il, jl, kl, xSeed, flowDoms, x, ib, jb, kb, &
         ie, je, ke, fringes, scratch, flowDomsd, xd
    use haloExchange, only : whalo1to1RealGeneric
    use oversetUtilities_d, only : newtonUpdate_d, fracToWeights_d
    use utils, only : setPointers_d

    implicit none

    ! Input
    integer(kind=intType), intent(in) :: level, sps
    type(commType), pointer :: commPattern
    type(internalCommType), pointer :: internal

    ! Working
    integer(kind=intType) :: nn, ii,jj, ierr,  i, j, k, d1, i1, j1, k1, d2, i2, j2, k2
    integer(kind=intType) :: size, procID, index, iii,jjj
    integer, dimension(mpi_status_size) :: mpiStatus
    real(kind=realType) :: frac(3), fracd(3), frac0(3), xCen(3), xCend(3), weight(8)
    integer(kind=intType), dimension(8), parameter :: indices=(/1,2,4,3,5,6,8,7/)

    ! Pointers to the overset comms to make it easier to read
    commPattern => commPatternOverset(level, sps)
    internal => internalOverset(level, sps)

    ! Step 1: Since we need to update donors for all cells including the
    ! donors for double halos, we must know the new cell center
    ! locations for the all of these receivers. Unfortunately, the
    ! double halos don't have coordinates, so we must first perform a
    ! (forward) block-to-block halo exchange to populate the xSeed
    ! values for all cells, including double halos.

    do nn=1, nDom
       call setPointers_d(nn, level, sps)

       if (.not. associated(flowDoms(nn, level, sps)%xSeed)) then
          allocate(flowDoms(nn, level, sps)%XSeed(0:ib, 0:jb, 0:kb, 3))
       end if
       xSeed => flowDoms(nn, level, sps)%xSeed
       xSeed = zero
       do k=2,kl
          do j=2, jl
             do i=2, il
                xSeed(i, j, k, :) = eighth*(&
                     x(i-1, j-1, k-1, :) + &
                     x(i  , j-1, k-1, :) + &
                     x(i-1, j  , k-1, :) + &
                     x(i  , j  , k-1, :) + &
                     x(i-1, j-1, k  , :) + &
                     x(i  , j-1, k  , :) + &
                     x(i-1, j  , k  , :) + &
                     x(i  , j  , k  , :)) !+ fringes(i,j,k)%offset

                ! Offset is not active so the xSeed_d just has the x
                ! part. Just dump the values into scratch so we don't
                ! acllocate any additional memory.
                scratch(i, j, k, 1:3) = eighth*(&
                     xd(i-1, j-1, k-1, :) + &
                     xd(i  , j-1, k-1, :) + &
                     xd(i-1, j  , k-1, :) + &
                     xd(i  , j  , k-1, :) + &
                     xd(i-1, j-1, k  , :) + &
                     xd(i  , j-1, k  , :) + &
                     xd(i-1, j  , k  , :) + &
                     xd(i  , j  , k  , :))
             end do
          end do
       end do
    end do

    ! Exchange the xSeeds.
    do nn=1, nDom
       flowDoms(nn, level, sps)%realCommVars(1)%var => flowDoms(nn, level, sps)%xSeed(:, :, :, 1)
       flowDoms(nn, level, sps)%realCommVars(2)%var => flowDoms(nn, level, sps)%xSeed(:, :, :, 2)
       flowDoms(nn, level, sps)%realCommVars(3)%var => flowDoms(nn, level, sps)%xSeed(:, :, :, 3)
       flowDoms(nn, level, sps)%realCommVars(4)%var => flowDoms(nn, level, sps)%scratch(:, :, :, 1)
       flowDoms(nn, level, sps)%realCommVars(5)%var => flowDoms(nn, level, sps)%scratch(:, :, :, 2)
       flowDoms(nn, level, sps)%realCommVars(6)%var => flowDoms(nn, level, sps)%scratch(:, :, :, 3)
    end do

    ! Run the (foward) generic halo exchange.
    call wHalo1to1RealGeneric(6, level, sps, commPatternCell_2nd, internalCell_2nd)

    ! Step 2: Next we need to communicate the xSeeds to their donor
    ! procs. This means running the overset exchange in REVERSE (ie from
    ! receiver to donor).  Most of this code will look like
    ! wOverset_b. We will runt he newtonUpdate code (below) on the fly
    ! as we receive the data, which should hide some of the comm time.

    ! Gather up the seeds into the *recv* buffer. Note we loop over
    ! nProcRECV here! After the buffer is assembled it is send off.

    jj = 1
    ii = 1
    recvs: do i=1,commPattern%nProcRecv

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern%recvProc(i)
       size    = 6*commPattern%nrecv(i)

       ! Copy the data into the buffer

       do j=1,commPattern%nrecv(i)

          ! Store the block and the indices to make code a bit easier to read

          d2 = commPattern%recvList(i)%block(j)
          i2 = commPattern%recvList(i)%indices(j,1)
          j2 = commPattern%recvList(i)%indices(j,2)
          k2 = commPattern%recvList(i)%indices(j,3)

          ! Copy the xSeed and it's derivative
          recvBuffer(jj)   = flowDoms(d2,level,sps)%xSeed(i2,j2,k2,1)
          recvBuffer(jj+1) = flowDoms(d2,level,sps)%xSeed(i2,j2,k2,2)
          recvBuffer(jj+2) = flowDoms(d2,level,sps)%xSeed(i2,j2,k2,3)
          recvBuffer(jj+3) = flowDoms(d2,level,sps)%scratch(i2,j2,k2,1)
          recvBuffer(jj+4) = flowDoms(d2,level,sps)%scratch(i2,j2,k2,2)
          recvBuffer(jj+5) = flowDoms(d2,level,sps)%scratch(i2,j2,k2,3)

          jj = jj + 6
       end do

       ! Send the data.
       call mpi_isend(recvBuffer(ii), size, adflow_real, procID,  &
            procID, ADflow_comm_world, recvRequests(i), &
            ierr)

       ! Set ii to jj for the next processor.

       ii = jj

    end do recvs

    ! Post the nonblocking receives.

    ii = 1
    sends: do i=1,commPattern%nProcSend

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern%sendProc(i)
       size    = 6*commPattern%nsend(i)

       ! Post the receive.

       call mpi_irecv(sendBuffer(ii), size, adflow_real, procId, &
            myId, ADflow_comm_world, sendRequests(i), ierr)

       ! And update ii.

       ii = ii + size

    enddo sends

    ! Do the local interpolation.

    localInterp: do i=1,internal%ncopy

       ! Store the block and the indices of the donor a bit easier.

       d1 = internal%donorBlock(i)
       i1 = internal%donorIndices(i, 1)
       j1 = internal%donorIndices(i, 2)
       k1 = internal%donorIndices(i, 3)

       ! Idem for the halo's.

       d2 = internal%haloBlock(i)
       i2 = internal%haloIndices(i, 1)
       j2 = internal%haloIndices(i, 2)
       k2 = internal%haloIndices(i, 3)

       xCen = flowDoms(d2, level, sps)%xSeed(i2, j2, k2, :)
       xCend = flowDoms(d2, level, sps)%scratch(i2, j2, k2, 1:3)
       frac0 = (/half, half, half/)
       call newtonUpdate_d(xCen, xCend, &
            flowDoms(d1, level, sps)%x(i1-1:i1+1, j1-1:j1+1, k1-1:k1+1, :), &
            flowDomsd(d1, level, sps)%x(i1-1:i1+1, j1-1:j1+1, k1-1:k1+1, :), &
            frac0, frac, fracd)

       ! Set the new weights
       call fracToWeights_d(frac, fracd, weight, internal%donorInterpd(i, :))

    enddo localInterp

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern%nProcSend
    completeSends: do i=1,commPattern%nProcSend

       ! Complete any of the requests.

       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)


       ii = index

       jj = 6*commPattern%nsendCum(ii-1)
       do j=1,commPattern%nsend(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern%sendList(ii)%block(j)
          i2 = commPattern%sendList(ii)%indices(j,1)
          j2 = commPattern%sendList(ii)%indices(j,2)
          k2 = commPattern%sendList(ii)%indices(j,3)

          xCen = sendBuffer(jj+1:jj+3)
          xCend = sendBuffer(jj+4:jj+6)
          jj = jj + 6

          ! Compute new fraction
          frac0 = (/half, half, half/)
          call newtonUpdate_d(xCen, xCend, &
               flowDoms(d2, level, sps)%x(i2-1:i2+1, j2-1:j2+1, k2-1:k2+1, :), &
               flowDomsd(d2, level, sps)%x(i2-1:i2+1, j2-1:j2+1, k2-1:k2+1, :), &
               frac0, frac, fracd)

          ! Set the new weights
          call fracToWeights_d(frac, fracd, weight, &
               commPattern%sendList(ii)%interpd(j, :))
       enddo

    enddo completeSends

    ! Complete the nonblocking sends.

    size = commPattern%nProcRecv
    do i=1,commPattern%nProcRecv
       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)
    enddo

  end subroutine updateOversetConnectivity_d

  subroutine updateOversetConnectivity_b(level, sps)

    ! Reverse  mode linearization of updateOversetConnectivity

    use constants
    use communication
    use blockPointers, only : nDom, il, jl, kl, xSeed, flowDoms, x, ib, jb, kb, &
         ie, je, ke, fringes, scratch, flowDomsd, xd
    use haloExchange, only : whalo1to1RealGeneric_b
    use oversetUtilities_b, only : newtonUpdate_b, fracToWeights_b, newtonUpdate, fracToWeights
    use utils, only : setPointers_b

    implicit none

    ! Input
    integer(kind=intType), intent(in) :: level, sps
    type(commType), pointer :: commPattern
    type(internalCommType), pointer :: internal

    ! Working
    integer(kind=intType) :: nn, ii,jj, kk, ierr,  i, j, k, d1, i1, j1, k1, d2, i2, j2, k2
    integer(kind=intType) :: size, procID, index, iii,jjj
    integer, dimension(mpi_status_size) :: mpiStatus
    real(kind=realType) :: frac(3), fracd(3), frac0(3), xCen(3), xCend(3), weight(8), add(3)
    integer(kind=intType), dimension(8), parameter :: indices=(/1,2,4,3,5,6,8,7/)

    ! Pointers to the overset comms to make it easier to read
    commPattern => commPatternOverset(level, sps)
    internal => internalOverset(level, sps)

    ! Zero out xSeedd (in scratch)
    do nn=1, nDom
       flowDoms(nn, 1, sps)%scratch(:, :, :, 1:3) = zero
    end do

    ! In reverse the fist thing we must do is the compute the
    ! sensitivites of xCen and send it to the receiving
    ! processor. This comm pattern is the same as wOverset forward.

    ii = 1
    sends: do i=1,commPattern%nProcSend

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern%sendProc(i)
       size    = 3*commPattern%nsend(i)

       ! Copy the data in the correct part of the send buffer.

       jj = ii
       do j=1,commPattern%nsend(i)

          ! Store the block id and the indices of the donor
          ! a bit easier.
          d1 = commPattern%sendList(i)%block(j)
          i1 = commPattern%sendList(i)%indices(j,1)
          j1 = commPattern%sendList(i)%indices(j,2)
          k1 = commPattern%sendList(i)%indices(j,3)

          ! -------- Recompute forward pass -------------
          xCen = commPattern%sendList(i)%xCen(j, :)

          ! Do newton update
          frac0 = (/half, half, half/)
          call newtonUpdate(xCen, &
               flowDoms(d1, level, sps)%x(i1-1:i1+1, j1-1:j1+1, k1-1:k1+1, :), frac0, frac)

          ! Set the new weights
          call fracToWeights(frac, weight)

          ! ------------- Reverse pass -----------

          ! Transfer the weights back to the frac
          call fracToWeights_b(frac, fracd, weight, commPattern%sendList(i)%interpd(j, :))

          ! Run the reverse newton update.  Note that we are
          ! accumulating into the local xd here in the newton_b call.
          frac0 = (/half, half, half/)
          call newtonUpdate_b(xCen, xCend, &
               flowDoms(d1, level, sps)%x(i1-1:i1+1, j1-1:j1+1, k1-1:k1+1, :), &
               flowDomsd(d1, level, sps)%x(i1-1:i1+1, j1-1:j1+1, k1-1:k1+1, :), &
               frac0, frac, fracd)

          ! -------------------------------------

          ! We want to send xCend to the receiver
          sendBuffer(jj:jj+2) = xCend
          jj = jj + 3

       enddo
       ! Send the data.

       call mpi_isend(sendBuffer(ii), size, adflow_real, procId,  &
            procId, ADflow_comm_world, sendRequests(i), &
            ierr)

       ! Set ii to jj for the next processor.

       ii = jj

    enddo sends

    ! Post the nonblocking receives.

    ii = 1
    receives: do i=1,commPattern%nProcRecv

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern%recvProc(i)
       size    = 3*commPattern%nrecv(i)

       ! Post the receive.

       call mpi_irecv(recvBuffer(ii), size, adflow_real, procId, &
            myId, ADflow_comm_world, recvRequests(i), ierr)

       ! And update ii.

       ii = ii + size

    enddo receives

    ! Do the local stuff while we're waiting

    localInterp: do i=1,internal%ncopy

       ! Store the block and the indices of the donor a bit easier.
       d1 = internal%donorBlock(i)
       i1 = internal%donorIndices(i, 1)
       j1 = internal%donorIndices(i, 2)
       k1 = internal%donorIndices(i, 3)

       d2 = internal%haloBlock(i)
       i2 = internal%haloIndices(i, 1)
       j2 = internal%haloIndices(i, 2)
       k2 = internal%haloIndices(i, 3)

       ! -------- Recompute forward pass -------------
       xCen = internal%XCen(i, :)

       ! Do newton update
       frac0 = (/half, half, half/)
       call newtonUpdate(xCen, &
            flowDoms(d1, level, sps)%x(i1-1:i1+1, j1-1:j1+1, k1-1:k1+1, :), frac0, frac)

       ! Set the new weights
       call fracToWeights(frac, weight)

      ! ------------- Reverse pass -----------

       ! Transfer the weights back to the frac
       call fracToWeights_b(frac, fracd, weight, internal%Donorinterpd(i, :))

       ! Run the reverse newton update.  Note that we are
       ! accumulating into the local xd here in the newton_b call.
       frac0 = (/half, half, half/)
       call newtonUpdate_b(xCen, xCend, &
            flowDoms(d1, level, sps)%x(i1-1:i1+1, j1-1:j1+1, k1-1:k1+1, :), &
            flowDomsd(d1, level, sps)%x(i1-1:i1+1, j1-1:j1+1, k1-1:k1+1, :), &
            frac0, frac, fracd)

       ! Accumulate in to the xSeedd (which we're using scratch for)
       flowDoms(d2, level, sps)%scratch(i2, j2, k2, 1:3) =  &
            flowDoms(d2, level, sps)%scratch(i2, j2, k2, 1:3) + xCend

    end do localInterp

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern%nProcRecv
    completeRecvs: do i=1,commPattern%nProcRecv

       ! Complete any of the requests.

       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)

       ! Copy the data just arrived in the halo's.

       ii = index
       jj = 3*commPattern%nrecvCum(ii-1)
       do j=1,commPattern%nrecv(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern%recvList(ii)%block(j)
          i2 = commPattern%recvList(ii)%indices(j,1)
          j2 = commPattern%recvList(ii)%indices(j,2)
          k2 = commPattern%recvList(ii)%indices(j,3)

          flowDoms(d2, level, sps)%scratch(i2, j2, k2, 1:3) = &
               flowDoms(d2, level, sps)%scratch(i2, j2, k2, 1:3) + recvBuffer(jj+1:jj+3)

          jj =jj + 3
       enddo
    end do completeRecvs

    ! Complete the nonblocking sends.

    size = commPattern%nProcSend
    do i=1,commPattern%nProcSend
       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
    enddo

    ! Now we have accumulated back as far as xSeedd (stored in
    ! scratch). We can now do the whalo1to1_b

    ! Exchange the xSeeds in reverse.
    do nn=1, nDom
       flowDoms(nn, level, sps)%realCommVars(1)%var => flowDoms(nn, level, sps)%scratch(:, :, :, 1)
       flowDoms(nn, level, sps)%realCommVars(2)%var => flowDoms(nn, level, sps)%scratch(:, :, :, 2)
       flowDoms(nn, level, sps)%realCommVars(3)%var => flowDoms(nn, level, sps)%scratch(:, :, :, 3)
    end do

    call wHalo1to1RealGeneric_b(3, level, sps, commPatternCell_2nd, internalCell_2nd)

    ! Finaly we can push back to the local x
    do nn=1, nDom
       call setPointers_b(nn, level, sps)

       do k=2,kl
          do j=2, jl
             do i=2, il
                ! Add is accumulate seed for xSeed (stored in scratch )
                add = eighth * scratch(i, j, k, 1:3)
                do kk=k-1,k
                   do jj=j-1,j
                      do ii=i-1,i
                         xd(ii, jj, kk, :) = xd(ii, jj, kk, :) + add
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine updateOversetConnectivity_b
#endif
end module oversetCommUtilities
