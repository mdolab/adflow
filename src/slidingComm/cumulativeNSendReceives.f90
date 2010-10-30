!
!      ******************************************************************
!      *                                                                *
!      * File:          cumulativeNSendReceives.f90                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-17-2003                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine cumulativeNSendReceives(commSlidingCell)
!
!      ******************************************************************
!      *                                                                *
!      * cumulativeNSendReceives allocates the memory for and           *
!      * determines the values of the cumulative storage arrays of the  *
!      * number of sends and receives for the given external sliding    *
!      * mesh communication pattern.                                    *
!      *                                                                *
!      ******************************************************************
!
       use commSliding
       implicit none
!
!      Subroutine argument
!
       type(slidingCommType), intent(inout) :: commSlidingCell
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Abbreviate the number of sending and receiving processors a
       ! bit easier and allocate the memory for the cumulative storage
       ! arrays.

       nn = commSlidingCell%nProcSend
       mm = commSlidingCell%nProcRecv

       allocate(commSlidingCell%nSendCum(0:nn), &
                commSlidingCell%nRecvCum(0:mm), stat=ierr)
       if(ierr /= 0)                                 &
         call terminate("cumulativeNSendReceives", &
                        "Memory allocation failure for nSendCum &
                        &and nRecvCum.")

       ! Determine the values of these arrays.

       commSlidingCell%nSendCum(0) = 0
       commSlidingCell%nRecvCum(0) = 0

       do nn=1,commSlidingCell%nProcSend
         commSlidingCell%nSendCum(nn) =          &
                commSlidingCell%nSendCum(nn-1) + &
                commSlidingCell%nSend(nn)
       enddo

       do nn=1,commSlidingCell%nProcRecv
         commSlidingCell%nRecvCum(nn) =          &
                commSlidingCell%nRecvCum(nn-1) + &
                commSlidingCell%nRecv(nn)
       enddo

       end subroutine cumulativeNSendReceives
