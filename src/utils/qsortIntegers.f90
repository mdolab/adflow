!
!      ******************************************************************
!      *                                                                *
!      * File:          qsortIntegers.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-19-2002                                      *
!      * Last modified: 03-23-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine qsortIntegers(arr, nn)
!
!      ******************************************************************
!      *                                                                *
!      * qsortIntegers sorts the given number of integers in            *
!      * increasing order.                                              *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), dimension(*), intent(inout) :: arr
       integer(kind=intType), intent(in)                  :: nn
!
!      Local variables
!
       integer(kind=intType), parameter :: m = 7

       integer(kind=intType) :: nStack
       integer(kind=intType) :: i, j, k, r, l, jStack, ii

       integer :: ierr

       integer(kind=intType) :: a, tmp

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
       if(ierr /= 0)                     &
         call terminate("qsortIntegers", &
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
           ! Choose median of left, center and right elements as partitioning
           ! element a. Also rearrange so that (l) <= (l+1) <= (r).

           k = (l+r)/2
           tmp      = arr(k)      ! Wwap the elements
           arr(k)   = arr(l+1)    ! k and l+1.
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
             if(ierr /= 0)                     &
               call terminate("qsortIntegers", &
                              "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                     &
               call terminate("qsortIntegers", &
                              "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                     &
               call terminate("qsortIntegers", &
                              "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                     &
               call terminate("qsortIntegers", &
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
       if(ierr /= 0)                      &
         call terminate("qsortIntegers", &
                        "Deallocation error for stack")

       ! Check in debug mode whether the array is really sorted.

       if( debug ) then
         do i=1,(nn-1)
           if(arr(i+1) < arr(i))              &
             call terminate("qsortIntegers", &
                            "Array is not sorted correctly")
         enddo
       endif

       end subroutine qsortIntegers
