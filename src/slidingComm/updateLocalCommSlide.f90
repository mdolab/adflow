!
!      ******************************************************************
!      *                                                                *
!      * File:          updateLocalCommSlide.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-14-2003                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine updateLocalCommSlide(haloInd, donorInfo, level, sps)
!
!      ******************************************************************
!      *                                                                *
!      * updateLocalCommSlide updates intSlidingCell_1st and            *
!      * intSlidingCell_2nd of the given level and spectral solution    *
!      * using the information stored in haloInd and donorInfo.         *
!      * Afterwards the memory of the member variables of donorInfo     *
!      * is deallocated.                                                *
!      *                                                                *
!      ******************************************************************
!
       use commSliding
       use updateComm
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps
       integer(kind=intType), dimension(*), intent(in) :: haloInd
       type(updateCommType), intent(inout) :: donorInfo
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
       ! Update the internal communication pattern for the 1st and 2nd
       ! level halo cells.

       call updateIntSlidingCell(donorInfo%nCopy1st, &
                                 intSlidingCell_1st(level,sps))
       call updateIntSlidingCell(donorInfo%nCopy2nd, &
                                 intSlidingCell_2nd(level,sps))

       ! Release the memory of the member variables of donorInfo.

       deallocate(donorInfo%indBuf,  donorInfo%block, &
                  donorInfo%indices, donorInfo%weight, stat=ierr)
       if(ierr /= 0)                            &
         call terminate("updateLocalCommSlide", &
                        "Deallocation error for the member variables &
                        &of donorInfo.")

       !=================================================================

       contains

         !===============================================================

         subroutine updateIntSlidingCell(nn, intSlidingCell)
!
!        ****************************************************************
!        *                                                              *
!        * updateIntSlidingCell updates the internal communication      *
!        * pattern stored in intSlidingCell using the first nn          *
!        * elements of donorInfo.                                       *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=intType), intent(in) :: nn
         type(internalSlidingCommType), intent(inout) :: intSlidingCell
!
!        Local variables.
!
         integer(kind=intType) :: nOld, ii, jj, i
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
         ! Store the old value of nCopy and determine the new value.

         nOld = intSlidingCell%nCopy
         ii   = nOld + nn
         intSlidingCell%nCopy = ii

         ! Reallocate the memory for the indices and the weight.

         call reallocateInteger(intSlidingCell%donorList%block, &
                                ii, nOld, .true.)
         call reallocateInteger(intSlidingCell%haloList%block, &
                                ii, nOld, .true.)

         call reallocateInteger2(intSlidingCell%donorList%indices, &
                                 ii, 3_intType, nOld, 3_intType, .true.)
         call reallocateInteger2(intSlidingCell%haloList%indices, &
                                 ii, 3_intType, nOld, 3_intType, .true.)

         call reallocateReal(intSlidingCell%weight, ii, nOld, .true.)

         ! Copy the info from donorInfo and haloInd into
         ! intSlidingCell.

         do i=1,nn

           ! Store the halo indices in haloList.

           ii = i + nOld
           jj = 4*(donorInfo%indBuf(i)-1)

           intSlidingCell%haloList%block(ii)     = haloInd(jj+1)
           intSlidingCell%haloList%indices(ii,1) = haloInd(jj+2)
           intSlidingCell%haloList%indices(ii,2) = haloInd(jj+3)
           intSlidingCell%haloList%indices(ii,3) = haloInd(jj+4)

           ! Store the donor info in donorList.

           intSlidingCell%donorList%block(ii) = donorInfo%block(i)

           intSlidingCell%donorList%indices(ii,1) = &
                                               donorInfo%indices(i,1)
           intSlidingCell%donorList%indices(ii,2) = &
                                               donorInfo%indices(i,2)
           intSlidingCell%donorList%indices(ii,3) = &
                                               donorInfo%indices(i,3)

           ! And finally the weights.

           intSlidingCell%weight(ii) = donorInfo%weight(i)

         enddo

         end subroutine updateIntSlidingCell

       end subroutine updateLocalCommSlide
