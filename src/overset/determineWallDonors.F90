subroutine determineWallDonors(level, sps, MAGIC)

  ! This subroutine communicates the "localWallFringes" array and
  ! determine which of the fringes are wall donors. 
  use communication
  use overset
  implicit none


  ! Input params
  integer(kind=intType), intent(in) :: level, sps, MAGIC

  ! Local Variables
  integer(kind=intType) :: i, j, k, ii, jj, kk, iii, jjj, kkk, nn, tag, ierr, n
  integer(kind=intType) :: nCopy,  iSize, iStart, iEnd, iProc, sendCount, nWallFringeProc
  integer(kind=intType), dimension(:), allocatable :: fringeProc, cumFringeProc
  integer(kind=intType), dimension(:), allocatable :: wallFringeProc, cumWallFringeProc

  ! MPI Stuff
  logical :: barrierActive, barrierDone, flag
  integer :: barrierRequest
  integer status(MPI_STATUS_SIZE) 

  ! Sort the wall fringes and determine the processor splits. 
  call qsortFringeType(localWallFringes, nLocalWallFringe)
  allocate(wallFringeProc(nProc), cumWallFringeProc(1:nProc+1))
  call computeFringeProcArray(localWallFringes, nLocalWallFringe, &
       wallFringeProc, cumWallFringeProc, nWallFringeProc)

  sendCount = 0
  do j=1, nWallFringeProc

     iProc = wallFringeProc(j)
     iStart = cumWallFringeProc(j)
     iEnd = cumWallFringeProc(j+1)-1
     if (iProc == myid) then 
        do i=iStart, iEnd
           nn = localWallFringes(i)%donorBlock

           do kk=0, 1
              do jj=0, 1
                 do ii=0, 1

                    iii = ii + localWallFringes(i)%dI
                    jjj = jj + localWallFringes(i)%dJ
                    kkk = kk + localWallFringes(i)%dK

                    flowDoms(nn, level, sps)%fringes(iii, jjj, kkk)%isWallDonor = .True.
                 end do
              end do
           end do
        end do
     else

        ! I need to send these wall fringes to the donor proc
        sendCount = sendCount + 1

        ! Note the tag here: The previous comm used tags
        ! 1:nDomTotal. These start at nDomTotal + 1 such that there is
        ! no overlap.

        tag = 6*MAGIC + iProc + 1

        call mpi_issend(localWallFringes(iStart:iEnd), iEnd-iStart+1, &
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

     ! Check if a message is ready
     tag = status(MPI_TAG)
     if (flag .and. tag >= 6*MAGIC+1 .and. tag <= 7*MAGIC) then 

        call receiveFringes(status, n)

        ! Do the flagging of wall donors on the fly
        do j=1, n
           nn = tmpFringes(j)%donorBlock

           do ii=0, 1
              do jj=0, 1
                 do kk=0, 1
                    iii = ii + tmpFringes(j)%dI
                    jjj = jj + tmpFringes(j)%dJ
                    kkk = kk + tmpFringes(j)%dK

                    flowDoms(nn, level, sps)%fringes(iii, jjj, kkk)%isWallDonor = .True.
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
  deallocate(wallFringeProc, cumWallFringeProc)

end subroutine determineWallDonors
