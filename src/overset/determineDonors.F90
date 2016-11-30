!       This subroutine takes a list of all local fringes, sorts it,   
!       then communicates the required fringes back to their donor     
!       procs, and flags all required cells as donors. This routine    
!       does double duty: If the flag useWall is set, it flags the     
!       donor cells as wallDonors, otherwise it flags the cells as     
!       regular donors. We do this since the communication structure   
!       for these two operations are identical                         

subroutine determineDonors(level, sps, fringeList, nFringe, useWall)

  use constants
  use block, only : fringeType, flowDoms
  use communication, only : adflow_comm_world, myid, nProc, recvRequests, sendRequests
  use utils, only : Echk, setPointers
  use overset, only : clusters, nDomTotal, nClusters
  use oversetUtilities, only : qsortFringeType, setIsDonor, setIsWallDonor, &
       computeFringeProcArray, unwindIndex
  implicit none

  ! Input Params
  integer(kind=intType), intent(in) :: level, sps, nFringe
  type(fringeType), intent(inout), dimension(nFringe) :: fringeList
  logical, intent(in) :: useWall

  ! Working
  integer(kind=intType), dimension(:), allocatable :: fringeProc, cumFringeProc
  integer(kind=intType), dimension(:), allocatable :: tmpInt
  integer(kind=intType), dimension(:), allocatable :: recvSizes
  integer(kind=intType), dimension(:), allocatable :: intSendBuf, intRecvBuf
  integer(kind=intType) :: i, j, k, ii, jj, kk, iii, jjj,  kkk, nn, index
  integer(kind=intType) :: il, jl, kl, dIndex
  integer(kind=intType) :: iStart, iEnd, iProc, iSize, nFringeProc
  integer(kind=intType) :: sendCount, recvCount, ierr, totalRecvSize
  integer mpiStatus(MPI_STATUS_SIZE) 

  ! First sort the fringes such that they are grouped by destination
  ! procesor.
  call qsortFringeType(fringeList, nFringe, sortByDonor)

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
  call computeFringeProcArray(fringeList, nFringe, &
       fringeProc, cumFringeProc, nFringeProc)
  
  ! nFringeProc is the total number of donor processors from
  ! fringeList fringeProc(1:nProcFringe) are the donor processors
  ! for the fringes cumFringeProc(1:nFringeProc) are the cumulative
  ! offset from in the localFringe Array. Note that we have
  ! over-estimated their size as nProc.

  allocate(tmpInt(0:nProc-1), recvSizes(0:nProc-1))
  tmpInt = 0
  do j=1, nFringeProc
     iProc = fringeProc(j)
     if (iProc /= myid) then 
        tmpInt(iProc) = (cumFringeProc(j+1) - cumFringeProc(j))*2
     end if
  end do

  ! Sum how much data we must receive from each processor. 
  call mpi_alltoall(tmpInt, 1, adflow_integer, recvSizes, 1, adflow_integer, &
       adflow_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! We will need to send 2 integers to the donor processor:
  ! donorBlock, dIndex
  
  ! Allocate space for the sending and receiving buffers
  totalRecvSize = sum(recvSizes)
  allocate(intSendBuf(2*nFringe), intRecvBuf(totalRecvSize))
  
     ! Pack the full buffer with donorBlock, dI, dJ, and dK
  do j=1, nFringe
     intSendBuf(2*j-1) = fringeList(j)%donorBlock
     intSendbuf(2*j  ) = fringeList(j)%dIndex
  end do

  ! Send the donors back to their own processors.
  sendCount = 0
  do j=1, nFringeProc
     
     iProc = fringeProc(j)
     iStart = (cumFringeProc(j)-1)*2 + 1
     iSize = (cumFringeProc(j+1) - cumFringeProc(j))*2
     
     if (iProc /= myid) then 
        sendCount = sendCount + 1
        call mpi_isend(intSendBuf(iStart), iSize, adflow_integer, iProc, myid, &
             adflow_comm_world, sendRequests(sendCount), ierr)
        call ECHK(ierr, __FILE__, __LINE__)
     end if
  end do
  
  ! Non-blocking receives
  recvCount = 0
  ii = 1
  do iProc=0, nProc-1
     
     if (recvSizes(iProc) > 0) then
        recvCount = recvCount + 1
        call mpi_irecv(intRecvBuf(ii), recvSizes(iProc), adflow_integer, &
             iProc, iProc, adflow_comm_world, recvRequests(recvCount), ierr) 
        call ECHK(ierr, __FILE__, __LINE__) 
        
        ii = ii + recvSizes(iProc)
     end if
  end do

  ! Local Work to do while we wait for data to send/recv
  do j=1, nFringeProc

     iProc = fringeProc(j)
     iStart = cumFringeProc(j)
     iEnd = cumFringeProc(j+1)-1
     
     if (iProc == myid) then 
        do i=iStart, iEnd
           nn = fringeList(i)%donorBlock
           il = flowDoms(nn, level, sps)%il
           jl = flowDoms(nn, level, sps)%jl
           kl = flowDoms(nn, level, sps)%kl
           dIndex = fringeList(i)%dIndex
           call unwindIndex(dIndex, il, jl, kl, iii, jjj, kkk)

           ! For the wall donors, we just flag the 1 cell that was
           ! identified in the fringe search based on the octant. 
           if (useWall) then 
              call setIsWallDonor(flowDoms(nn, level, sps)%status(iii, jjj, kkk), .True. )
           else
              do kk=0, 1
                 do jj=0, 1
                    do ii=0, 1
                       call setIsDonor(&
                            flowDoms(nn, level, sps)%status(iii+ii, jjj+jj, kkk+kk), .True. )
                    end do
                 end do
              end do
           end if
        end do
     end if
  end do
  
  ! Complete all the sends/receives. We could do overlapping here
  ! like the frist comm for the fringes/blocks. 
  do i=1, recvCount
     call mpi_waitany(recvCount, recvRequests, index, mpiStatus, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  enddo

  do i=1, sendCount
     call mpi_waitany(sendCount, sendRequests, index, mpiStatus, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  enddo
  
  ! Loop over the full receive buffer that should now be full. 
  do j=1, totalRecvSize/2

     nn = intRecvBuf(2*j-1)
     il = flowDoms(nn, level, sps)%il
     jl = flowDoms(nn, level, sps)%jl
     kl = flowDoms(nn, level, sps)%kl

     dIndex = intRecvBuf(2*j)
     call unwindIndex(dIndex, il, jl, kl, iii, jjj, kkk)

     if (useWall) Then 
        call setIsWallDonor(flowDoms(nn, level, sps)%status(iii, jjj, kkk), .True. )
     else
        do kk=0, 1
           do jj=0, 1
              do ii=0, 1
                 call setIsDonor(flowDoms(nn, level, sps)%status(iii+ii, jjj+jj, kkk+kk), .True. )
              end do
           end do
        end do
     end if
  end do
  ! Finished with the buffers and allocatable arrays
  deallocate(intSendBuf, intRecvBuf, fringeProc, cumFringeProc, tmpInt, recvSizes)

end subroutine determineDonors
