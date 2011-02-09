!
!      ******************************************************************
!      *                                                                *
!      * File:          initMemSliding.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-26-2003                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initMemSliding(level)
!
!      ******************************************************************
!      *                                                                *
!      * initMemSliding initializes the data structures for the         *
!      * communication pattern of the sliding mesh interfaces for the   *
!      * given grid level.                                              *
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
       integer(kind=intType) :: sps
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of spectral solutions.

       do sps=1,nTimeIntervalsSpectral

         ! Initialize the data structure for the external communication
         ! for both the first and second level halo's.

         call initCommSlidingCell(commSlidingCell_1st(level,sps))
         call initCommSlidingCell(commSlidingCell_2nd(level,sps))

         ! Initialize the data structure for the internal communication
         ! for both the first and second level halo's.

         call initIntSlidingCell(intSlidingCell_1st(level,sps))
         call initIntSlidingCell(intSlidingCell_2nd(level,sps))

       enddo

       end subroutine initMemSliding

!      ==================================================================

       subroutine initCommSlidingCell(commSlidingCell)

!      ******************************************************************
!      *                                                                *
!      * initCommSlidingCell initializes the memory of                  *
!      * commSlidingCell, the external communication pattern for        *
!      * sliding mesh halo's.                                           *
!      *                                                                *
!      ******************************************************************
!
       use commSliding
       implicit none
!
!      Subroutine arguments.
!
       type(slidingCommType), intent(out) :: commSlidingCell
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
       ! Initialization of the integer data.

       commSlidingCell%nProcSend = 0
       commSlidingCell%nProcRecv = 0

       ! Initial allocation of the arrays of commSlidingCell.

       allocate(commSlidingCell%sendProc(0), &
                commSlidingCell%recvProc(0), &
                commSlidingCell%nSend(0),    &
                commSlidingCell%nRecv(0),    &
                commSlidingCell%sendList(0), &
                commSlidingCell%recvList(0), stat=ierr)
       if(ierr /= 0)                           &
         call terminate("initCommSlidingCell", &
                        "Memory allocation failure for commSlidingCell")

       end subroutine initCommSlidingCell

!      ==================================================================

       subroutine initIntSlidingCell(intSlidingCell)
!      ******************************************************************
!      *                                                                *
!      * initIntSlidingCell initializes the memory of                   *
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
       type(internalSlidingCommType), intent(out) :: intSlidingCell
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
       ! Initialization of the integer data.

       intSlidingCell%nSlidingHalos = 0
       intSlidingCell%nCopy          = 0

       ! Initial allocation of the arrays of intSlidingCell

       allocate(intSlidingCell%slidingHaloList%block(0),     &
                intSlidingCell%slidingHaloList%indices(0,0), &
                intSlidingCell%donorList%block(0),           &
                intSlidingCell%donorList%indices(0,0),       &
                intSlidingCell%haloList%block(0),            &
                intSlidingCell%haloList%indices(0,0),        &
                intSlidingCell%rotIndex(0),                  &
                intSlidingCell%weight(0), stat=ierr)
       if(ierr /= 0)                          &
         call terminate("initIntSlidingCell", &
                        "Memory allocation failure for intSlidingCell")

       end subroutine initIntSlidingCell
