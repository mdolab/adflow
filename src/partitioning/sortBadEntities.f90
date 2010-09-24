!
!      ******************************************************************
!      *                                                                *
!      * File:          sortBadEntities.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-30-2004                                      *
!      * Last modified: 04-20-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module fourIntPlusRealDataType
!
!      ******************************************************************
!      *                                                                *
!      * This local module contains the derived data type, which        *
!      * consists of four integers and one real.                        *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

       public
       private :: lessEqualFourIntPlusRealType
       private :: lessFourIntPlusRealType
       private :: notEqualFourIntPlusRealType

       ! Definition of the derived data type.

       type fourIntPlusRealType
         integer(kind=intType) :: n1, n2, n3, n4
         real(kind=realType)   :: dist
       end type fourIntPlusRealType

       ! Interfaces for the definitions of the operators <=, < and /=.
       ! These are needed for the sorting of this derived data type.
       ! Note that the = operator does not need to be defined, because
       ! fourIntPlusRealType only contains primitive types.

       interface operator(<=)
         module procedure lessEqualFourIntPlusRealType
       end interface

       interface operator(<)
         module procedure lessFourIntPlusRealType
       end interface

       interface operator(/=)
         module procedure notEqualFourIntPlusRealType
       end interface

       !=================================================================

       contains

         !===============================================================
!
!        ****************************************************************
!        *                                                              *
!        * Functions to define the operators <, <= and /=.              *
!        * Note that the comparison is only based on the integers.      *
!        * The real contains additional info, the maximum deviation,    *
!        * which is normally different even if the subfaces are         *
!        * identical.                                                   *
!        *                                                              *
!        ****************************************************************
!
         logical function lessEqualFourIntPlusRealType(g1, g2)
         implicit none
         type(fourIntPlusRealType), intent(in) :: g1, g2

         ! Compare the first element.

         if(g1%n1 < g2%n1) then
           lessEqualFourIntPlusRealType = .true.
           return
         else if(g1%n1 > g2%n1) then
           lessEqualFourIntPlusRealType = .false.
           return
         endif

         ! Compare the second element.

         if(g1%n2 < g2%n2) then
           lessEqualFourIntPlusRealType = .true.
           return
         else if(g1%n2 > g2%n2) then
           lessEqualFourIntPlusRealType = .false.
           return
         endif

         ! Compare the third element.

         if(g1%n3 < g2%n3) then
           lessEqualFourIntPlusRealType = .true.
           return
         else if(g1%n3 > g2%n3) then
           lessEqualFourIntPlusRealType = .false.
           return
         endif

         ! Compare the fourth element.

         if(g1%n4 < g2%n4) then
           lessEqualFourIntPlusRealType = .true.
           return
         else if(g1%n4 > g2%n4) then
           lessEqualFourIntPlusRealType = .false.
           return
         endif

         ! g1 equals g2. Return .true.

         lessEqualFourIntPlusRealType = .true.

         end function lessEqualFourIntPlusRealType

         !===============================================================

         logical function lessFourIntPlusRealType(g1, g2)
         implicit none
         type(fourIntPlusRealType), intent(in) :: g1, g2

         ! Compare the first element.

         if(g1%n1 < g2%n1) then
           lessFourIntPlusRealType = .true.
           return
         else if(g1%n1 > g2%n1) then
           lessFourIntPlusRealType = .false.
           return
         endif

         ! Compare the second element.

         if(g1%n2 < g2%n2) then
           lessFourIntPlusRealType = .true.
           return
         else if(g1%n2 > g2%n2) then
           lessFourIntPlusRealType = .false.
           return
         endif

         ! Compare the third element.

         if(g1%n3 < g2%n3) then
           lessFourIntPlusRealType = .true.
           return
         else if(g1%n3 > g2%n3) then
           lessFourIntPlusRealType = .false.
           return
         endif

         ! Compare the fourth element.

         if(g1%n4 < g2%n4) then
           lessFourIntPlusRealType = .true.
           return
         else if(g1%n4 > g2%n4) then
           lessFourIntPlusRealType = .false.
           return
         endif

         ! g1 equals g2. Return .false.

         lessFourIntPlusRealType = .false.

         end function lessFourIntPlusRealType

         !===============================================================

         logical function notEqualFourIntPlusRealType(g1, g2)
         implicit none
         type(fourIntPlusRealType), intent(in) :: g1, g2

         notEqualFourIntPlusRealType = .true.
         if(g1%n1 == g2%n1 .and. g1%n2 == g2%n2 .and. &
            g1%n3 == g2%n3 .and. g1%n4 == g2%n4) &
            notEqualFourIntPlusRealType = .false.

         end function notEqualFourIntPlusRealType

       end module fourIntPlusRealDataType

!      ==================================================================

       subroutine sortBadEntities(nEntities, entities, dist, sortDist)
!
!      ******************************************************************
!      *                                                                *
!      * sortBadEntities sorts the given number of entities in          *
!      * increasing order and gets rid of the multiple entries.         *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use fourIntPlusRealDataType
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(inout) :: nEntities
       integer(kind=intType), dimension(4,*), intent(inout) :: entities

       real(kind=realType),   dimension(*),   intent(inout) :: dist

       logical, intent(in) :: sortDist
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm

       type(fourIntPlusRealType), dimension(nEntities) :: tmp
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if there are no entities to be sorted.

       if(nEntities == 0) return

       ! Copy the info into tmp. If the distances must be sorted as
       ! well, copy the info. Otherwise simply put zero.

       do nn=1,nEntities
         tmp(nn)%n1   = entities(1,nn)
         tmp(nn)%n2   = entities(2,nn)
         tmp(nn)%n3   = entities(3,nn)
         tmp(nn)%n4   = entities(4,nn)
         if( sortDist ) then
           tmp(nn)%dist = dist(nn)
         else
           tmp(nn)%dist = zero
         endif
       enddo

       ! Sort tmp in increasing order.

       call qsortFourIntPlusRealType(tmp, nEntities)

       ! Get rid of the multiple entries. Note that the exceptional
       ! case of zero entities does not to be considered, because in
       ! that case this part of the subroutine is not executed.
       ! If multiple entries are present the distance is taken as
       ! the maximum of the two.

       mm = 1
       do nn=2,nEntities
         if(tmp(nn) /= tmp(mm)) then
           mm = mm + 1
           tmp(mm) = tmp(nn)
         else
           tmp(mm)%dist = max(tmp(mm)%dist, tmp(nn)%dist)
         endif
       enddo

       ! Copy the data back info entities and dist. The latter
       ! only if the distances should be sorted as well.

       nEntities = mm

       do nn=1,nEntities
         entities(1,nn) = tmp(nn)%n1
         entities(2,nn) = tmp(nn)%n2
         entities(3,nn) = tmp(nn)%n3
         entities(4,nn) = tmp(nn)%n4
         if( sortDist ) dist(nn) = tmp(nn)%dist
       enddo

       end subroutine sortBadEntities

!      ==================================================================

       subroutine qsortFourIntPlusRealType(arr, nn)
!
!      ******************************************************************
!      *                                                                *
!      * qsortFourIntPlusRealType sorts the given number of halo's in   *
!      * increasing order based on the <= operator for this derived     *
!      * data type.                                                     *
!      *                                                                *
!      ******************************************************************
!
       use fourIntPlusRealDataType
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn

       type(fourIntPlusRealType), dimension(*), intent(inout) :: arr
!
!      Local variables.
!
       integer(kind=intType), parameter :: m = 7

       integer(kind=intType) :: nStack
       integer(kind=intType) :: i, j, k, r, l, jStack, ii

       integer :: ierr

       type(fourIntPlusRealType) :: a, tmp

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
       if(ierr /= 0)                                &
         call terminate("qsortFourIntPlusRealType", &
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
             if(ierr /= 0)                                &
               call terminate("qsortFourIntPlusRealType", &
                              "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                                &
               call terminate("qsortFourIntPlusRealType", &
                              "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                                &
               call terminate("qsortFourIntPlusRealType", &
                              "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                                &
               call terminate("qsortFourIntPlusRealType", &
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
       if(ierr /= 0)                                &
         call terminate("qsortFourIntPlusRealType", &
                        "Deallocation error for stack")

       ! Check in debug mode whether the array is really sorted.

       if( debug ) then
         do i=1,(nn-1)
           if(arr(i+1) < arr(i))                        &
             call terminate("qsortFourIntPlusRealType", &
                            "Array is not sorted correctly")
         enddo
       endif

       end subroutine qsortFourIntPlusRealType
