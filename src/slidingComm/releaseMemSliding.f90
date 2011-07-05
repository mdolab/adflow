!
!      ******************************************************************
!      *                                                                *
!      * File:          releaseMemSliding.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-26-2003                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine releaseMemSliding(level)
!
!      ******************************************************************
!      *                                                                *
!      * releaseMemSliding releases the memory of the sliding mesh      *
!      * communication pattern for the given grid level. This must be   *
!      * done, because every new time step the communication pattern    *
!      * changes and must be recomputed.                                *
!      *                                                                *
!      ******************************************************************
!
       use commSliding
       use inputTimeSpectral
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer(kind=intType) :: mm
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of spectral solutions.

       do mm=1,nTimeIntervalsSpectral

         ! Release the memory for the external communication pattern for
         ! the 1st and 2nd level sliding mesh halo's.

         call releaseCommSlidingCell(commSlidingCell_1st(level,mm))
         call releaseCommSlidingCell(commSlidingCell_2nd(level,mm))

         ! Release the memory for the internal communication pattern for
         ! the 1st and 2nd level sliding mesh halo's.

         call releaseIntSlidingCell(intSlidingCell_1st(level,mm))
         call releaseIntSlidingCell(intSlidingCell_2nd(level,mm))

       enddo

       end subroutine releaseMemSliding

!      ==================================================================

       subroutine releaseCommSlidingCell(commSlidingCell)

!      ******************************************************************
!      *                                                                *
!      * releaseCommSlidingCell releases the memory of                  *
!      * commSlidingCell, the external communication pattern sliding    *
!      * mesh halo's.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use commSliding
       implicit none
!
!      Subroutine arguments.
!
       type(slidingCommType), intent(inout) :: commSlidingCell
!
!      Local variables.
!
       integer :: ierr
       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       do nn=1,commSlidingCell%nProcSend

         deallocate(commSlidingCell%sendList(nn)%block,   &
                    commSlidingCell%sendList(nn)%indices, &
                    stat=ierr)
         if(ierr /= 0) call terminate("releaseCommSlidingCell", &
                                      "Deallocation error for sendList")
       enddo

       do nn=1,commSlidingCell%nProcRecv

         deallocate(commSlidingCell%recvList(nn)%indRecv, &
                    commSlidingCell%recvList(nn)%block,   &
                    commSlidingCell%recvList(nn)%indices, &
                    commSlidingCell%recvList(nn)%weight,  &
                    stat=ierr)
         if(ierr /= 0) call terminate("releaseCommSlidingCell", &
                                      "Deallocation error for recvList")

       enddo

       deallocate(commSlidingCell%sendProc, &
                  commSlidingCell%recvProc, &
                  commSlidingCell%nsend,    &
                  commSlidingCell%nrecv,    &
                  commSlidingCell%nsendCum, &
                  commSlidingCell%nrecvCum, &
                  commSlidingCell%sendList, &
                  commSlidingCell%recvList, &
                  stat=ierr)
       if(ierr /= 0)                              &
         call terminate("releaseCommSlidingCell", &
                        "Deallocation error for commSlidingCell")

       end subroutine releaseCommSlidingCell

!      ==================================================================

       subroutine releaseIntSlidingCell(intSlidingCell)
!
!      ******************************************************************
!      *                                                                *
!      * releaseIntSlidingCell releases the memory of                   *
!      * intSlidingCell, the internal communication pattern for         *
!      * sliding mesh halo's.                                           *
!      *                                                                *
!      ******************************************************************
!
       use commSliding
       implicit none
!
!      Subroutine arguments.
!
       type(internalSlidingCommType), intent(inout) :: intSlidingCell
!
!      Local variables.
!
       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       deallocate(intSlidingCell%slidingHaloList%block,   &
                  intSlidingCell%slidingHaloList%indices, &
                  intSlidingCell%donorList%block,         &
                  intSlidingCell%donorList%indices,       &
                  intSlidingCell%haloList%block,          &
                  intSlidingCell%haloList%indices,        &
                  intSlidingCell%rotIndex,                &
                  intSlidingCell%weight,                  &
                  stat=ierr)
       if(ierr /= 0)                             &
         call terminate("releaseIntSlidingCell", &
                        "Deallocation error for intSlidingCell")

       end subroutine releaseIntSlidingCell
