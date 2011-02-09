!
!      ******************************************************************
!      *                                                                *
!      * File:          updateExternalCommSlide.f90                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-17-2003                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine updateExternalCommSlide(haloInd, intRecv, realRecv, &
                                          level,   sps,     sendProc)
!
!      ******************************************************************
!      *                                                                *
!      * updateExternalCommSlide updates the receiving part of the      *
!      * external sliding mesh communication pattern for both the 1st   *
!      * and 2nd level halo's for the given multigrid level and         *
!      * spectral solution using the information provided in haloInd,   *
!      * intRecv and realRecv.                                          *
!      *                                                                *
!      ******************************************************************
!
       use commSliding
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps, sendProc
       integer(kind=intType), dimension(*), intent(in) :: haloInd
       integer(kind=intType), dimension(*), intent(in) :: intRecv
       real(kind=realType),   dimension(*), intent(in) :: realRecv
!
!      Local variables.
!
       integer(kind=intType) :: nCopy1st, nCopy2nd
       integer(kind=intType) :: nHalo1st, nHalo2nd
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the number of entities to copy and the number of halo's in
       ! the receive buffer for both the 1st and 2nd halo cells a bit
       ! easier.

       nCopy1st = intRecv(1)
       nCopy2nd = intRecv(2)
       nHalo1st = intRecv(3)
       nHalo2nd = intRecv(4)

       ! Update the receiving part of commSlidingCell1st and
       ! commSlidingCell2nd for this level.

       call updateCommSlidingCellRecv(nCopy1st, nHalo1st, &
                                      commSlidingCell_1st(level,sps))
       call updateCommSlidingCellRecv(nCopy2nd, nHalo2nd, &
                                      commSlidingCell_2nd(level,sps))

       !=================================================================

       contains

         !===============================================================

         subroutine updateCommSlidingCellRecv(nCopy, nHalo, &
                                              commSlidingCell)
!
!        ****************************************************************
!        *                                                              *
!        * UpdateCommSlidingCellRecv updates the receiving part of      *
!        * the external communication pattern for sliding mesh          *
!        * interfaces.                                                  *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=intType), intent(in) :: nCopy, nHalo
         type(slidingCommType), intent(inout) :: commSlidingCell
!
!        Local variables.
!
         integer :: ierr

         integer(kind=intType) :: i, jj, nn, mm, ll

         type(interpolHaloListType), pointer, dimension(:) :: &
                                                          tmpRecvList
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

           subroutine reallocateInteger2(intArray,           &
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

           !=============================================================

           subroutine reallocateReal(realArray, newSize, oldSize, &
                                     alwaysFreeMem)
             use precision
             implicit none

             real(kind=realType), dimension(:), pointer :: realArray
             integer(kind=intType), intent(in) :: newSize, oldSize
             logical, intent(in) :: alwaysFreeMem
           end subroutine reallocateReal
         end interface
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Return immediately if nCopy == 0. This will probably never
         ! occur, but it does not hurt to test it.

         if(nCopy == 0) return

         ! Find out if sendProc is already stored in the processors from
         ! which data must be received. As recvProc is not necessarily
         ! sorted a linear search algorithm is used. Normally this is
         ! okay, because the number of processors to communicate with is
         ! limited.

         do nn=1,commSlidingCell%nProcRecv
           if(sendProc == commSlidingCell%recvProc(nn)) exit
         enddo

         ! If sendProc is not in the list the memory of recvProc,
         ! nRecv and recvList is reallocated and the new entries are
         ! initialized accordingly. NRecvCum will be constructed later.

         testRealloc: if(nn > commSlidingCell%nProcRecv) then

           ! First reallocate the memory for the integer arrays.

           call reallocateInteger(commSlidingCell%recvProc, nn, &
                                  commSlidingCell%nProcRecv, .true.)
           call reallocateInteger(commSlidingCell%nRecv, nn, &
                                  commSlidingCell%nProcRecv, .true.)

           commSlidingCell%recvProc(nn) = sendProc
           commSlidingCell%nRecv(nn)    = 0

           ! The reallocation of recvList. This is slightly more
           ! complicated, because it is a derived datatype.

           ! First set the pointer for tmpRecvList.

           tmpRecvList => commSlidingCell%recvList

           ! Allocate the memory for recvList with one more entry.

           allocate(commSlidingCell%recvList(nn), stat=ierr)
           if(ierr /= 0)                                 &
             call terminate("updateCommSlidingCellRecv", &
                            "Memory allocation failure for recvList.")

           ! Set the values of nCopy and the pointers back to the
           ! previously stored data of recvList.

           do mm=1,commSlidingCell%nProcRecv

             commSlidingCell%recvList(mm)%nCopy = &
                                             tmpRecvList(mm)%nCopy
             commSlidingCell%recvList(mm)%indRecv => &
                                             tmpRecvList(mm)%indRecv
             commSlidingCell%recvList(mm)%block    => &
                                             tmpRecvList(mm)%block
             commSlidingCell%recvList(mm)%indices  => &
                                             tmpRecvList(mm)%indices
             commSlidingCell%recvList(mm)%weight   => &
                                             tmpRecvList(mm)%weight
           enddo

           ! Set nCopy to 0 and nullify the pointers of the new entry.

           commSlidingCell%recvList(nn)%nCopy = 0

           nullify(commSlidingCell%recvList(nn)%indRecv, &
                   commSlidingCell%recvList(nn)%block,   &
                   commSlidingCell%recvList(nn)%indices, &
                   commSlidingCell%recvList(nn)%weight)

           ! Release the memory of tmpRecvList.

           deallocate(tmpRecvList, stat=ierr)
           if(ierr /= 0)                                 &
             call terminate("updateCommSlidingCellRecv", &
                            "Deallocation error for tmpRecvList.")

           ! Set the new value of nProcRecv.

           commSlidingCell%nProcRecv = nn

         endif testRealloc

         ! Determine the new value of nCopy and reallocate the memory
         ! of indRecv, block, indices and weight of recvList. Store
         ! the old value of nCopy in mm.

         mm = commSlidingCell%recvList(nn)%nCopy
         ll = mm + nCopy
         commSlidingCell%recvList(nn)%nCopy = ll

         call reallocateInteger(commSlidingCell%recvList(nn)%indRecv, &
                                ll, mm, .false.)
         call reallocateInteger(commSlidingCell%recvList(nn)%block, &
                                ll, mm, .false.)
         call reallocateInteger2(commSlidingCell%recvList(nn)%indices, &
                                 ll, 3_intType, mm, 3_intType, .false.)
         call reallocateReal(commSlidingCell%recvList(nn)%weight, &
                             ll, mm, .false.)

         ! Loop over nCopy and store the data of intRecv and realRecv
         ! in recvList. Note that the actual data of intRecv starts at
         ! position 5, indicated by the counter ll.

         ll = 5
         do i=1,nCopy

           mm = mm + 1

           ! Determine the offset in haloInd and store the block and
           ! indices of the halo.

           jj = 4*(intRecv(ll)-1)

           commSlidingCell%recvList(nn)%block(mm)     = haloInd(jj+1)
           commSlidingCell%recvList(nn)%indices(mm,1) = haloInd(jj+2)
           commSlidingCell%recvList(nn)%indices(mm,2) = haloInd(jj+3)
           commSlidingCell%recvList(nn)%indices(mm,3) = haloInd(jj+4)

           ! Store the index in the actual receive buffer and its
           ! corresponding interpolation weight. Note that the offset due
           ! to contributions of previous interfaces must be taken into
           ! account in indRecv.

           commSlidingCell%recvList(nn)%indRecv(mm) = intRecv(ll+1) &
                                           + commSlidingCell%nRecv(nn)
           commSlidingCell%recvList(nn)%weight(mm)   = realRecv(i)

           ! Update ll by 2.

           ll = ll + 2
         enddo

         ! Determine the new value of nRecv from this processor.

         commSlidingCell%nRecv(nn) = commSlidingCell%nRecv(nn) &
                                   + nHalo

         end subroutine updateCommSlidingCellRecv

       end subroutine updateExternalCommSlide

