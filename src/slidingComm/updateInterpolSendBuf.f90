!
!      ******************************************************************
!      *                                                                *
!      * File:          updateInterpolSendBuf.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-17-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine updateInterpolSendBuf(intSend, realSend, donorInfo, &
                                        level,   sps,      recvProc)
!
!      ******************************************************************
!      *                                                                *
!      * updateInterpolSendBuf creates the integer and real send        *
!      * buffers, which contain the interpolation information needed by *
!      * recvProc, from the info present in donorInfo; recvProc is      *
!      * the processor ID of the communicator SUmb_comm_world. This     *
!      * routine furthermore updates the sending part of the 1st and    *
!      * 2nd level sliding mesh communication pattern for the given     *
!      * grid level. The memory of donorInfo is released afterwards.    *
!      *                                                                *
!      ******************************************************************
!
       use commSliding
       use updateComm
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: level, sps, recvProc

       integer(kind=intType), dimension(*), intent(out) :: intSend
       real(kind=realType),   dimension(*), intent(out) :: realSend

       type(updateCommType), intent(inout) :: donorInfo
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, nn, mm, nDom1st, nDom2nd

       type(sortedDonorCommType) :: tmp
       type(sortedDonorCommType), dimension(donorInfo%nCopy2nd) :: &
                                                            sortedDonors
!
!      Function definition.
!
       integer(kind=intType) :: bsearchSortedDonorCommType
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Copy the halo info in sortedDonors and sort this array in
       ! increasing order.

       do i=1,donorInfo%nCopy1st
         sortedDonors(i)%haloLevel = 1
         sortedDonors(i)%block      = donorInfo%block(i)
         sortedDonors(i)%indices(1) = donorInfo%indices(i,1)
         sortedDonors(i)%indices(2) = donorInfo%indices(i,2)
         sortedDonors(i)%indices(3) = donorInfo%indices(i,3)
       enddo

       do i=(donorInfo%nCopy1st+1),donorInfo%nCopy2nd
         sortedDonors(i)%haloLevel = 2
         sortedDonors(i)%block     = donorInfo%block(i)
         sortedDonors(i)%indices(1) = donorInfo%indices(i,1)
         sortedDonors(i)%indices(2) = donorInfo%indices(i,2)
         sortedDonors(i)%indices(3) = donorInfo%indices(i,3)
       enddo

       call qsortSortedDonorCommType(sortedDonors, donorInfo%nCopy2nd)

       ! Get rid of the multiple entries and determine the number of
       ! different donors for the 1st and 2nd level.

       nDom1st = 0
       if(sortedDonors(1)%haloLevel == 1) nDom1st = 1
       nDom2nd = 1
       do i=2,donorInfo%nCopy2nd
         if(sortedDonors(nDom2nd) < sortedDonors(i)) then
           nDom2nd = nDom2nd + 1
           sortedDonors(nDom2nd) = sortedDonors(i)
           if(sortedDonors(i)%haloLevel == 1) nDom1st = nDom1st + 1
         endif
       enddo

       ! Determine the contents of the integer and real send buffer.
       ! First store the number of 1st and 2nd level interpolation
       ! entities and donor cells in the first 4 entries of intSend.

       intSend(1) = donorInfo%nCopy1st
       intSend(2) = donorInfo%nCopy2nd
       intSend(3) = nDom1st
       intSend(4) = nDom2nd

       ! Loop over the number of entities stored in donorInfo and
       ! determine the corresponding info for intSend and realSend.

       nn = 5
       do i=1,donorInfo%nCopy2nd

         ! Create the variable of type sortedDonorCommType needed
         ! for the binary search.

         tmp%haloLevel = 1
         if(i > donorInfo%nCopy1st) tmp%haloLevel = 2

         tmp%block      = donorInfo%block(i)
         tmp%indices(1) = donorInfo%indices(i,1)
         tmp%indices(2) = donorInfo%indices(i,2)
         tmp%indices(3) = donorInfo%indices(i,3)

         ! Search for tmp in sortedDonors. Entry must be found, which
         ! is checked in debug mode.

         mm = bsearchSortedDonorCommType(tmp, sortedDonors, nDom2nd)

         if( debug ) then
           if(mm == 0)                               &
             call terminate("updateInterpolSendBuf", &
                            "Entry not found in sortedDonors.")
         endif

         ! Store the information in intSend and realSend. In intSend
         ! the original index in the send buffer is stored as well as
         ! the index in the receive buffer of the final communication
         ! pattern, which is mm. The interpolation weight is simply
         ! copied in realSend.

         intSend(nn)   = donorInfo%indBuf(i)
         intSend(nn+1) = mm
         nn = nn + 2

         realSend(i) = donorInfo%weight(i)

       enddo

       ! Release the memory of the member variables of donorInfo.

       deallocate(donorInfo%indBuf,  donorInfo%block, &
                  donorInfo%indices, donorInfo%weight, stat=ierr)
       if(ierr /= 0)                             &
         call terminate("updateInterpolSendBuf", &
                        "Deallocation error for the member variables &
                        &of donorInfo.")

       ! Update the sending part of commSlidingCell1st and
       ! commSlidingCell2nd for this level.

       call updateCommSlidingCellSend(nDom1st, &
                                      commSlidingCell_1st(level,sps))
       call updateCommSlidingCellSend(nDom2nd, &
                                      commSlidingCell_2nd(level,sps))

       !=================================================================

       contains

         !===============================================================

         subroutine updateCommSlidingCellSend(nDom, commSlidingCell)
!
!        ****************************************************************
!        *                                                              *
!        * updateCommSlidingCellSend updates the sending part of        *
!        * the external communication pattern for sliding mesh          *
!        * interfaces. It stores the information of sortedDonors in     *
!        * the correct place of commSlidingCell.                        *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=intType), intent(in)    :: nDom
         type(slidingCommType), intent(inout) :: commSlidingCell
!
!        Local variables.
!
         type(slidingCommListType), dimension(:), pointer :: &
                                                           tmpSendList
!
!        Interfaces
!
         interface
           subroutine reallocateInteger(intArray, newSize, oldSize, &
                                        alwaysFreeMem)
             use precision
             implicit none

             integer(kind=intType), dimension(:), pointer :: intArray
             integer(kind=intType), intent(in) :: newSize, oldSize
             logical, intent(in) :: alwaysFreeMem
           end subroutine reallocateInteger

           !=============================================================

           subroutine reallocateInteger2(intArray,            &
                                          newSize1, newSize2, &
                                          oldSize1, oldSize2, &
                                          alwaysFreeMem)
             use precision
             implicit none
 
             integer(kind=intType), dimension(:,:), pointer :: intArray
             integer(kind=intType), intent(in) :: newSize1, newSize2, &
                                                  oldSize1, oldSize2
             logical, intent(in) :: alwaysFreeMem
           end subroutine reallocateInteger2
         end interface
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Return immediately if nDom == 0. This will probably never
         ! occur, but it does not hurt to test it.

         if(nDom == 0) return

         ! Find out if recvProc is already stored in the processors to
         ! which data must be sent. As sendProc is not necessarily
         ! sorted a linear search algorithm is used. Normally this is
         ! okay, because the number of processors to communicate with is
         ! limited.

         do nn=1,commSlidingCell%nProcSend
           if(recvProc == commSlidingCell%sendProc(nn)) exit
         enddo

         ! If recvProc is not in the list the memory of sendProc,
         ! nsend and sendList is reallocated and the new entries are
         ! initialized accordingly. NsendCum will be constructed later.

         testRealloc: if(nn > commSlidingCell%nProcSend) then

           ! First reallocate the memory of the integer arrays.

           call reallocateInteger(commSlidingCell%sendProc, nn, &
                                  commSlidingCell%nProcSend, .true.)
           call reallocateInteger(commSlidingCell%nsend, nn, &
                                  commSlidingCell%nProcSend, .true.)

           commSlidingCell%sendProc(nn) = recvProc
           commSlidingCell%nsend(nn)    = 0

           ! The reallocation of sendList. This is slightly more
           ! complicated, because it is a derived datatype.

           ! First set the pointer for tmpSendList.

           tmpSendList => commSlidingCell%sendList

           ! Allocate the memory for sendList with one more entry.

           allocate(commSlidingCell%sendList(nn), stat=ierr)
           if(ierr /= 0)                                 &
             call terminate("updateCommSlidingCellSend", &
                            "Memory allocation failure for sendList.")

           ! Set the pointers back for the previously stored data
           ! of sendList.

           do mm=1,commSlidingCell%nProcSend

             commSlidingCell%sendList(mm)%block   => &
                                             tmpSendList(mm)%block
             commSlidingCell%sendList(mm)%indices => &
                                             tmpSendList(mm)%indices
           enddo

           ! Nullify the pointers of the new entry.

           nullify(commSlidingCell%sendList(nn)%block, &
                   commSlidingCell%sendList(nn)%indices)

           ! Release the memory of tmpSendList.

           deallocate(tmpSendList, stat=ierr)
           if(ierr /= 0)                                 &
             call terminate("updateCommSlidingCellSend", &
                            "Deallocation error for tmpSendList.")

           ! Set the new value of nProcSend.

           commSlidingCell%nProcSend = nn

         endif testRealloc

         ! Determine the new value of nsend(nn) and reallocate the
         ! memory of block and indices of sendList. Store the old
         ! value in mm.

         mm = commSlidingCell%nsend(nn)
         commSlidingCell%nsend(nn) = mm + nDom

         call reallocateInteger(commSlidingCell%sendList(nn)%block, &
                                commSlidingCell%nsend(nn), mm, .false.)

         call reallocateInteger2(commSlidingCell%sendList(nn)%indices, &
                                 commSlidingCell%nsend(nn), 3_intType, &
                                 mm, 3_intType, .false.)

         ! Loop over nDom and store the donor indices in sendList.

         do i=1,nDom
           mm = mm + 1

           commSlidingCell%sendList(nn)%block(mm) = &
                                             sortedDonors(i)%block

           commSlidingCell%sendList(nn)%indices(mm,1) = &
                                             sortedDonors(i)%indices(1)
           commSlidingCell%sendList(nn)%indices(mm,2) = &
                                             sortedDonors(i)%indices(2)
           commSlidingCell%sendList(nn)%indices(mm,3) = &
                                             sortedDonors(i)%indices(3)
         enddo

         end subroutine updateCommSlidingCellSend

       end subroutine updateInterpolSendBuf
