!
!      ******************************************************************
!      *                                                                *
!      * File:          sortRangesSplitInfo.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-04-2005                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module sortSubRange
!
!      ******************************************************************
!      *                                                                *
!      * This local module contains the derived datatype as well as the *
!      * functions needed to sort the subranges.                        *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

       public
       private :: lessEqualSortSubRangeType
       private :: lessSortSubRangeType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived datatype sortSubRangeType.       *
!      *                                                                *
!      ******************************************************************
!
       type sortSubRangeType

         ! iMin: minimum i-index in the subrange.
         ! jMin: minimum j-index in the subrange.
         ! kMin: minimum k-index in the subrange.
         ! iMax: maximum i-index in the subrange.
         ! jMax: maximum j-index in the subrange.
         ! kMax: maximum k-index in the subrange.

         integer(kind=intType) :: iMin, jMin, kMin
         integer(kind=intType) :: iMax, jMax, kMax

       end type sortSubRangeType

       ! Interfaces for the definitions of the operators <=, < and /=.
       ! These are needed for the sorting of this derived data type.
       ! Note that the = operator does not need to be defined, because
       ! sortSubRangeType only contains primitive types.

       interface operator(<=)
         module procedure lessEqualSortSubRangeType
       end interface

       interface operator(<)
         module procedure lessSortSubRangeType
       end interface

       !=================================================================

       contains

         !===============================================================

         logical function lessEqualSortSubRangeType(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * lessEqualSortSubRangeType defines the operator <= for the    *
!        * derived datatype sortSubRangeType. The comparison is first   *
!        * based on kMin, followed by jMin and finally iMin.            *
!        * The comparison is therefore not based on the max values.     *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(sortSubRangeType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin executation.                                           *
!        *                                                              *
!        ****************************************************************
!
         ! Compare the kMin index and return .true. or .false. if they
         ! differ.

         if(g1%kMin < g2%kMin) then
           lessEqualSortSubRangeType = .true.
           return
         else if(g1%kMin > g2%kMin) then
           lessEqualSortSubRangeType = .false.
           return
         endif

         ! kMin indices are equal. Compare the jMin's.

         if(g1%jMin < g2%jMin) then
           lessEqualSortSubRangeType = .true.
           return
         else if(g1%jMin > g2%jMin) then
           lessEqualSortSubRangeType = .false.
           return
         endif

         ! Also the jMin's are equal. Compare iMin's.

         if(g1%iMin < g2%iMin) then
           lessEqualSortSubRangeType = .true.
           return
         else if(g1%iMin > g2%iMin) then
           lessEqualSortSubRangeType = .false.
           return
         endif

         ! g1 equals g2. Return .true.

         lessEqualSortSubRangeType = .true.

         end function lessEqualSortSubRangeType

         !===============================================================

         logical function lessSortSubRangeType(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * lessSortSubRangeType defines the operator < for the derived  *
!        * datatype sortSubRangeType. The comparison is first based on  *
!        * kMin, followed by jMin and finally iMin.                     *
!        * The comparison is therefore not based on the max values.     *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(sortSubRangeType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin executation.                                           *
!        *                                                              *
!        ****************************************************************
!
         ! Compare the kMin index and return .true. or .false. if they
         ! differ.

         if(g1%kMin < g2%kMin) then
           lessSortSubRangeType = .true.
           return
         else if(g1%kMin > g2%kMin) then
           lessSortSubRangeType = .false.
           return
         endif

         ! kMin indices are equal. Compare the jMin's.

         if(g1%jMin < g2%jMin) then
           lessSortSubRangeType = .true.
           return
         else if(g1%jMin > g2%jMin) then
           lessSortSubRangeType = .false.
           return
         endif

         ! Also the jMin's are equal. Compare iMin's.

         if(g1%iMin < g2%iMin) then
           lessSortSubRangeType = .true.
           return
         else if(g1%iMin > g2%iMin) then
           lessSortSubRangeType = .false.
           return
         endif

         ! g1 equals g2. Return .false.

         lessSortSubRangeType = .false.

         end function lessSortSubRangeType

       end module sortSubRange

!      ==================================================================

       subroutine sortRangesSplitInfo(splitInfo)
!
!      ******************************************************************
!      *                                                                *
!      * sortRangesSplitInfo sort the ranges of the given subblocks in  *
!      * increasing order such that a unique ordering is obtained,      *
!      * independent of the history of the splitting.                   *
!      *                                                                *
!      ******************************************************************
!
       use partitionMod
       use sortSubRange
       implicit none
!
!      Subroutine arguments.
!
       type(splitCGNSType), intent(inout) :: splitInfo
!
!      Local variables.
!
       integer(kind=intType) :: i, nSubBlocks

       type(sortSubRangeType), dimension(splitInfo%nSubBlocks) :: subRanges
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Copy the subface range from splitInfo into subRanges.

       nSubBlocks = splitInfo%nSubBlocks

       do i=1,nSubBlocks
         subRanges(i)%iMin = splitInfo%ranges(i,1,1)
         subRanges(i)%jMin = splitInfo%ranges(i,2,1)
         subRanges(i)%kMin = splitInfo%ranges(i,3,1)

         subRanges(i)%iMax = splitInfo%ranges(i,1,2)
         subRanges(i)%jMax = splitInfo%ranges(i,2,2)
         subRanges(i)%kMax = splitInfo%ranges(i,3,2)
       enddo

       ! Sort subRanges in increasing order.

       call qsortSortSubRangeType(subRanges, nSubBlocks)

       ! Copy the data back into splitInfo.

       do i=1,nSubBlocks
         splitInfo%ranges(i,1,1) = subRanges(i)%iMin
         splitInfo%ranges(i,2,1) = subRanges(i)%jMin
         splitInfo%ranges(i,3,1) = subRanges(i)%kMin

         splitInfo%ranges(i,1,2) = subRanges(i)%iMax
         splitInfo%ranges(i,2,2) = subRanges(i)%jMax
         splitInfo%ranges(i,3,2) = subRanges(i)%kMax
       enddo

       end subroutine sortRangesSplitInfo

!      ==================================================================

       subroutine qsortSortSubRangeType(arr, nn)
!
!      ******************************************************************
!      *                                                                *
!      * qsortSortSubRangeType sorts the given number of halo's in      *
!      * increasing order based on the <= operator for this derived     *
!      * data type.                                                     *
!      *                                                                *
!      ******************************************************************
!
       use sortSubRange
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn

       type(sortSubRangeType), dimension(*), intent(inout) :: arr
!
!      Local variables.
!
       integer(kind=intType), parameter :: m = 7

       integer(kind=intType) :: nStack
       integer(kind=intType) :: i, j, k, r, l, jStack, ii

       integer :: ierr

       type(sortSubRangeType) :: a, tmp

       integer(kind=intType), allocatable, dimension(:) :: stack
       integer(kind=intType), allocatable, dimension(:) :: tmpStack
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for stack.

       nStack = 100
       allocate(stack(nStack), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("qsortSortSubRangeType", &
                        "Memory allocation failure for stack")

       ! Initialize the variables that control the sorting.

       jStack = 0
       l      = 1
       r      = nn

       ! Start of the algorithm

       do

         ! Check for the size of the subarray.

         if((r-l) < m) then

           ! Perform insertion sort

           do j=l+1,r
             a = arr(j)
             do i=(j-1),l,-1
               if(arr(i) <= a) exit
               arr(i+1) = arr(i)
             enddo
             arr(i+1) = a
           enddo

           ! In case there are no more elements on the stack, exit from
           ! the outermost do-loop. Algorithm has finished.

           if(jStack == 0) exit

           ! Pop stack and begin a new round of partitioning.

           r = stack(jStack)
           l = stack(jStack-1)
           jStack = jStack - 2

         else

           ! Subarray is larger than the threshold for a linear sort.
           ! Choose median of left, center and right elements as
           ! partitioning element a.
           ! Also rearrange so that (l) <= (l+1) <= (r).

           k = (l+r)/2
           tmp      = arr(k)             ! Swap the elements
           arr(k)   = arr(l+1)           ! k and l+1.
           arr(l+1) = tmp

           if(arr(r) < arr(l)) then
             tmp    = arr(l)             ! Swap the elements
             arr(l) = arr(r)             ! r and l.
             arr(r) = tmp
           endif

           if(arr(r) < arr(l+1)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(r)           ! r and l+1.
             arr(r)   = tmp
           endif

           if(arr(l+1) < arr(l)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(l)           ! l and l+1.
             arr(l)   = tmp
           endif

           ! Initialize the pointers for partitioning.

           i = l+1
           j = r
           a = arr(l+1)

           ! The innermost loop

           do

             ! Scan up to find element >= a.
             do
               i = i+1
               if(a <= arr(i)) exit
             enddo

             ! Scan down to find element <= a.
             do
               j = j-1
               if(arr(j) <= a) exit
             enddo

             ! Exit the loop in case the pointers i and j crossed.

             if(j < i) exit

             ! Swap the element i and j.

             tmp    = arr(i)
             arr(i) = arr(j)
             arr(j) = tmp
           enddo

           ! Swap the entries j and l+1. Remember that a equals
           ! arr(l+1).

           arr(l+1) = arr(j)
           arr(j)   = a

           ! Push pointers to larger subarray on stack,
           ! process smaller subarray immediately.

           jStack = jStack + 2
           if(jStack > nStack) then

             ! Storage of the stack is too small. Reallocate.

             allocate(tmpStack(nStack), stat=ierr)
             if(ierr /= 0)                             &
               call terminate("qsortSortSubRangeType", &
                              "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                             &
               call terminate("qsortSortSubRangeType", &
                              "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                             &
               call terminate("qsortSortSubRangeType", &
                              "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                             &
               call terminate("qsortSortSubRangeType", &
                              "Deallocation error for tmpStack")
           endif

           if((r-i+1) >= (j-l)) then
             stack(jStack)   = r
             r               = j-1
             stack(jStack-1) = j
           else
             stack(jStack)   = j-1
             stack(jStack-1) = l
             l               = j
           endif

         endif
       enddo

       ! Release the memory of stack.

       deallocate(stack, stat=ierr)
       if(ierr /= 0)                             &
         call terminate("qsortSortSubRangeType", &
                        "Deallocation error for stack")

       ! Check in debug mode whether the array is really sorted.

       if( debug ) then
         do i=1,(nn-1)
           if(arr(i+1) < arr(i))                     &
             call terminate("qsortSortSubRangeType", &
                            "Array is not sorted correctly")
         enddo
       endif

       end subroutine qsortSortSubRangeType
