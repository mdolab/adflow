module sorting

  use utils, only : terminate
contains

  function bsearchIntegers(key, base)
    !
    !       bsearchIntegers returns the index in base where key is stored. 
    !       A binary search algorithm is used here, so it is assumed that  
    !       base is sorted in increasing order. In case key appears more   
    !       than once in base, the result is arbitrary. If key is not      
    !       found, a zero is returned.                                     
    !
    use precision
    implicit none
    !
    !      Function type
    !
    integer(kind=intType) :: bsearchIntegers
    !
    !      Function arguments.
    !
    integer(kind=intType), intent(in)               :: key
    integer(kind=intType), dimension(:), intent(in) :: base
    integer(kind=intType)                           :: nn
    !
    !      Local variables.
    !
    integer(kind=intType) :: ii, pos, start
    logical               :: entryFound

    ! Initialize some values.

    start       = 1
    ii          = size(base)
    entryFound = .false.

    ! Binary search to find key.

    do
       ! Condition for breaking the loop

       if(ii == 0) exit

       ! Determine the position in the array to compare.

       pos = start + ii/2

       ! In case this is the entry, break the search loop.

       if(base(pos) == key) then
          entryFound = .true.
          exit
       endif

       ! In case the search key is larger than the current position,
       ! only parts to the right must be searched. Remember that base
       ! is sorted in increasing order. Nothing needs to be done if the
       ! key is smaller than the current element.

       if(key > base(pos)) then
          start = pos +1
          ii    = ii -1
       endif

       ! Modify ii for the next branch to search.

       ii = ii/2
    enddo

    ! Set bsearchIntegers. This depends whether the key was found.

    if( entryFound ) then
       bsearchIntegers = pos
    else
       bsearchIntegers = 0
    endif

  end function bsearchIntegers

  ! ----------------------------------------------------------------------
  !                                                                      |
  !                    No Tapenade Routine below this line               |
  !                                                                      |
  ! ----------------------------------------------------------------------

#ifndef  USE_TAPENADE


  subroutine qsortIntegers(arr, nn)
    !
    !       qsortIntegers sorts the given number of integers in            
    !       increasing order.                                              
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

  subroutine qsortReals(arr, nn)
    !
    !       qsortReals sorts the given number of reals in increasing       
    !       order.                                                         
    !
    use constants
    implicit none
    !
    !      Subroutine arguments
    !
    real(kind=realType), dimension(*), intent(inout) :: arr
    integer(kind=intType), intent(in)                :: nn
    !
    !      Local variables
    !
    integer(kind=intType), parameter :: m = 7

    integer(kind=intType) :: nStack
    integer(kind=intType) :: i, j, k, r, l, jStack, ii

    integer :: ierr

    real(kind=realType) :: a, tmp

    integer(kind=intType), allocatable, dimension(:) :: stack
    integer(kind=intType), allocatable, dimension(:) :: tmpStack

    ! Allocate the memory for stack.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                  &
         call terminate("qsortReals", &
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
          tmp      = arr(k)      ! Swap the elements
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
             if(ierr /= 0)                  &
                  call terminate("qsortReals", &
                  "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                  &
                  call terminate("qsortReals", &
                  "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                  &
                  call terminate("qsortReals", &
                  "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                  &
                  call terminate("qsortReals", &
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
    if(ierr /= 0)                  &
         call terminate("qsortReals", &
         "Deallocation error for stack")

    ! Check in debug mode whether the array is really sorted.

    if( debug ) then
       do i=1,(nn-1)
          if(arr(i+1) < arr(i))          &
               call terminate("qsortReals", &
               "Array is not sorted correctly")
       enddo
    endif

  end subroutine qsortReals

  subroutine qsortStrings(arr, nn)
    !
    !       qsortStrings sorts the given number of strings in increasing   
    !       order.                                                         
    !
    use constants
    implicit none
    !
    !      Subroutine arguments
    !
    character(len=*), dimension(*), intent(inout) :: arr
    integer(kind=intType), intent(in)             :: nn
    !
    !      Local variables
    !
    integer(kind=intType), parameter :: m = 7

    integer(kind=intType) :: nStack
    integer(kind=intType) :: i, j, k, r, l, jStack, ii

    integer :: ierr

    character(len=maxStringLen) :: a, tmp

    integer(kind=intType), allocatable, dimension(:) :: stack
    integer(kind=intType), allocatable, dimension(:) :: tmpStack

    ! Allocate the memory for stack.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                    &
         call terminate("qsortStrings", &
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
          tmp      = arr(k)      ! Swap the elements
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
             if(ierr /= 0)                    &
                  call terminate("qsortStrings", &
                  "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                    &
                  call terminate("qsortStrings", &
                  "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                    &
                  call terminate("qsortStrings", &
                  "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                    &
                  call terminate("qsortStrings", &
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
    if(ierr /= 0)                    &
         call terminate("qsortStrings", &
         "Deallocation error for stack")

    ! Check in debug mode whether the array is really sorted.

    if( debug ) then
       do i=1,(nn-1)
          if(arr(i+1) < arr(i))            &
               call terminate("qsortStrings", &
               "Array is not sorted correctly")
       enddo
    endif

  end subroutine qsortStrings


  function bsearchReals(key, base)
    !
    !       bsearchReals returns the index in base where key is stored.    
    !       A binary search algorithm is used here, so it is assumed that  
    !       base is sorted in increasing order. In case key appears more   
    !       than once in base, the result is arbitrary. If key is not      
    !       found, a zero is returned.                                     
    !
    use precision
    implicit none
    !
    !      Function type
    !
    integer(kind=intType) :: bsearchReals
    !
    !      Function arguments.
    !
    real(kind=realType), intent(in)               :: key
    real(kind=realType), dimension(:), intent(in) :: base
    integer(kind=intType)                         :: nn
    !
    !      Local variables.
    !
    integer(kind=intType) :: ii, pos, start
    logical               :: entryFound

    ! Initialize some values.

    start       = 1
    ii          = size(base)
    entryFound = .false.

    ! Binary search to find key.

    do
       ! Condition for breaking the loop

       if(ii == 0) exit

       ! Determine the position in the array to compare.

       pos = start + ii/2

       ! In case this is the entry, break the search loop.

       if(base(pos) == key) then
          entryFound = .true.
          exit
       endif

       ! In case the search key is larger than the current position,
       ! only parts to the right must be searched. Remember that base
       ! is sorted in increasing order. Nothing needs to be done if the
       ! key is smaller than the current element.

       if(key > base(pos)) then
          start = pos +1
          ii    = ii -1
       endif

       ! Modify ii for the next branch to search.

       ii = ii/2
    enddo

    ! Set bsearchReals. This depends whether the key was found.

    if( entryFound ) then
       bsearchReals = pos
    else
       bsearchReals = 0
    endif

  end function bsearchReals

  function bsearchStrings(key, base)
    !
    !       bsearchStrings returns the index in base where key is stored.  
    !       A binary search algorithm is used here, so it is assumed that  
    !       base is sorted in increasing order. In case key appears more   
    !       than once in base, the result is arbitrary. If key is not      
    !       found, a zero is returned.                                     
    !
    use precision
    implicit none
    !
    !      Function type
    !
    integer(kind=intType) :: bsearchStrings
    !
    !      Function arguments.
    !
    character(len=*), intent(in)               :: key
    character(len=*), dimension(:), intent(in) :: base
    integer(kind=intType)                      :: nn
    !
    !      Local variables.
    !
    integer(kind=intType) :: ii, pos, start
    logical               :: entryFound

    ! Initialize some values.

    start       = 1
    ii          = size(base)
    entryFound = .false.

    ! Binary search to find key.

    do
       ! Condition for breaking the loop

       if(ii == 0) exit

       ! Determine the position in the array to compare.

       pos = start + ii/2

       ! In case this is the entry, break the search loop.

       if(base(pos) == key) then
          entryFound = .true.
          exit
       endif

       ! In case the search key is larger than the current position,
       ! only parts to the right must be searched. Remember that base
       ! is sorted in increasing order. Nothing needs to be done if the
       ! key is smaller than the current element.

       if(key > base(pos)) then
          start = pos +1
          ii    = ii -1
       endif

       ! Modify ii for the next branch to search.

       ii = ii/2
    enddo

    ! Set bsearchStrings. This depends whether the key was found.

    if( entryFound ) then
       bsearchStrings = pos
    else
       bsearchStrings = 0
    endif

  end function bsearchStrings
  ! -----------------------------------------------------------------
  !
  ! This file contains two functions that are used to find the unique
  ! values from a list of integers as well as the inverse mapping. These
  ! routines have been generously borrowed from 
  !
  ! Michel Olagnon
  ! http://www.fortran-2000.com/rank/
  ! 
  ! Slight modifications of I_unirnk and I_uniinv for use with ADflow by
  ! Gaetan Kenway. The 'unique' subroutine was added that combines the
  ! functionality of both routines in a single call.

  subroutine unique(arr, nn, n_unique, inverse)
    use precision
    implicit none
    ! Input: 
    ! arr: Array of integers to be unique-sorted. Overwirtten on output 
    !      with unique values
    ! nn : Number of input values
    ! n_unique: Number of unique output values in arr. Only the first n_unique
    !           values of arr on output are meaningful
    ! inverse: size(nn): For each origianl entry in arr, this gives the index
    !          into to unique, sorted array. 

    ! Input Arguments
    integer(kind=intType), intent(in) :: nn
    integer(kind=intType), intent(inout), dimension(nn) :: arr

    ! Output Arguments
    integer(kind=intType), intent(out) :: n_unique
    integer(kind=intType), intent(out), dimension(nn) :: inverse

    ! Local Arguments
    integer(kind=intType), dimension(:), allocatable :: temp_arr, irngt
    integer(kind=intType) :: i

    allocate(temp_arr(nn), irngt(nn))
    ! Copy arr to temp array:
    temp_arr(:) = arr(:)

    call i_uniinv(arr, nn, inverse)
    call i_unirnk(arr, nn, irngt, n_unique)

    ! Since unirank is an arg sort, fill in sorted values into arr_unique

    do i=1,n_unique
       arr(i) = temp_arr(irngt(i))
    end do

    ! Fill remaining values of array with zeros since these have no
    ! meaning
    do i=n_unique+1,nn
       arr(i) = 0_intType
    end do
    deallocate(temp_arr, irngt)
  end subroutine unique


  Subroutine I_unirnk (XVALT, NVAL, IRNGT, NUNI)
    use precision
    implicit none

    ! __________________________________________________________
    !   UNIRNK = Merge-sort ranking of an array, with removal of
    !   duplicate entries.
    !   The routine is similar to pure merge-sort ranking, but on
    !   the last pass, it discards indices that correspond to
    !   duplicate entries.
    !   For performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.
    ! __________________________________________________________
    ! __________________________________________________________
    integer(kind=intType), intent(in) :: NVAL
    integer(kind=intType), Dimension (NVAL), Intent (In) :: XVALT

    integer(kind=intType), Dimension (NVAL), Intent (Out) :: IRNGT
    integer(kind=intType), Intent (Out) :: NUNI
    ! __________________________________________________________
    integer(kind=intType), Dimension (:), allocatable :: JWRKT
    integer(kind=intType) :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2
    integer(kind=intType) :: IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
    integer(kind=intType) :: XTST, XVALA, XVALB


    NUNI = NVAL
    !
    Select Case (NVAL)
    Case (:0)
       Return
    Case (1)
       IRNGT (1) = 1
       Return
    Case Default
       Continue
    End Select
    allocate(JWRKT(NVAL))
    !
    !  Fill-in the index array, creating ordered couples
    !
    Do IIND = 2, NVAL, 2
       If (XVALT(IIND-1) < XVALT(IIND)) Then
          IRNGT (IIND-1) = IIND - 1
          IRNGT (IIND) = IIND
       Else
          IRNGT (IIND-1) = IIND
          IRNGT (IIND) = IIND - 1
       End If
    End Do
    If (Modulo(NVAL, 2) /= 0) Then
       IRNGT (NVAL) = NVAL
    End If
    !
    !  We will now have ordered subsets A - B - A - B - ...
    !  and merge A and B couples into     C   -   C   - ...
    !
    LMTNA = 2
    LMTNC = 4
    !
    !  First iteration. The length of the ordered subsets goes from 2 to 4
    !
    Do
       If (NVAL <= 4) Exit
       !
       !   Loop on merges of A and B into C
       !
       Do IWRKD = 0, NVAL - 1, 4
          If ((IWRKD+4) > NVAL) Then
             If ((IWRKD+2) >= NVAL) Exit
             !
             !   1 2 3
             !
             If (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) Exit
             !
             !   1 3 2
             !
             If (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) Then
                IRNG2 = IRNGT (IWRKD+2)
                IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                IRNGT (IWRKD+3) = IRNG2
                !
                !   3 1 2
                !
             Else
                IRNG1 = IRNGT (IWRKD+1)
                IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                IRNGT (IWRKD+2) = IRNG1
             End If
             Exit
          End If
          !
          !   1 2 3 4
          !
          If (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) Cycle
          !
          !   1 3 x x
          !
          If (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) Then
             IRNG2 = IRNGT (IWRKD+2)
             IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
             If (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) Then
                !   1 3 2 4
                IRNGT (IWRKD+3) = IRNG2
             Else
                !   1 3 4 2
                IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                IRNGT (IWRKD+4) = IRNG2
             End If
             !
             !   3 x x x
             !
          Else
             IRNG1 = IRNGT (IWRKD+1)
             IRNG2 = IRNGT (IWRKD+2)
             IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
             If (XVALT(IRNG1) <= XVALT(IRNGT(IWRKD+4))) Then
                IRNGT (IWRKD+2) = IRNG1
                If (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) Then
                   !   3 1 2 4
                   IRNGT (IWRKD+3) = IRNG2
                Else
                   !   3 1 4 2
                   IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                   IRNGT (IWRKD+4) = IRNG2
                End If
             Else
                !   3 4 1 2
                IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                IRNGT (IWRKD+3) = IRNG1
                IRNGT (IWRKD+4) = IRNG2
             End If
          End If
       End Do
       !
       !  The Cs become As and Bs
       !
       LMTNA = 4
       Exit
    End Do
    !
    !  Iteration loop. Each time, the length of the ordered subsets
    !  is doubled.
    !
    Do
       If (2*LMTNA >= NVAL) Exit
       IWRKF = 0
       LMTNC = 2 * LMTNC
       !
       !   Loop on merges of A and B into C
       !
       Do
          IWRK = IWRKF
          IWRKD = IWRKF + 1
          JINDA = IWRKF + LMTNA
          IWRKF = IWRKF + LMTNC
          If (IWRKF >= NVAL) Then
             If (JINDA >= NVAL) Exit
             IWRKF = NVAL
          End If
          IINDA = 1
          IINDB = JINDA + 1
          !
          !  One steps in the C subset, that we create in the final rank array
          !
          !  Make a copy of the rank array for the iteration
          !
          JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
          XVALA = XVALT (JWRKT(IINDA))
          XVALB = XVALT (IRNGT(IINDB))
          !
          Do
             IWRK = IWRK + 1
             !
             !  We still have unprocessed values in both A and B
             !
             If (XVALA > XVALB) Then
                IRNGT (IWRK) = IRNGT (IINDB)
                IINDB = IINDB + 1
                If (IINDB > IWRKF) Then
                   !  Only A still with unprocessed values
                   IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                   Exit
                End If
                XVALB = XVALT (IRNGT(IINDB))
             Else
                IRNGT (IWRK) = JWRKT (IINDA)
                IINDA = IINDA + 1
                If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                XVALA = XVALT (JWRKT(IINDA))
             End If
             !
          End Do
       End Do
       !
       !  The Cs become As and Bs
       !
       LMTNA = 2 * LMTNA
    End Do
    !
    !   Last merge of A and B into C, with removal of duplicates.
    !
    IINDA = 1
    IINDB = LMTNA + 1
    NUNI = 0
    !
    !  One steps in the C subset, that we create in the final rank array
    !
    JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
    If (IINDB <= NVAL) Then
       XTST = I_NEARLESS (Min(XVALT(JWRKT(1)), XVALT(IRNGT(IINDB))))
    Else
       XTST = I_NEARLESS (XVALT(JWRKT(1)))
    Endif
    Do IWRK = 1, NVAL
       !
       !  We still have unprocessed values in both A and B
       !
       If (IINDA <= LMTNA) Then
          If (IINDB <= NVAL) Then
             If (XVALT(JWRKT(IINDA)) > XVALT(IRNGT(IINDB))) Then
                IRNG = IRNGT (IINDB)
                IINDB = IINDB + 1
             Else
                IRNG = JWRKT (IINDA)
                IINDA = IINDA + 1
             End If
          Else
             !
             !  Only A still with unprocessed values
             !
             IRNG = JWRKT (IINDA)
             IINDA = IINDA + 1
          End If
       Else
          !
          !  Only B still with unprocessed values
          !
          IRNG = IRNGT (IWRK)
       End If
       If (XVALT(IRNG) > XTST) Then
          XTST = XVALT (IRNG)
          NUNI = NUNI + 1
          IRNGT (NUNI) = IRNG
       End If
       !
    End Do
    deallocate(JWRKT)
    Return
    !
  End Subroutine I_unirnk

  Subroutine I_uniinv (XDONT, NVAL, IGOEST)
    use precision
    implicit none
    ! __________________________________________________________
    !   UNIINV = Merge-sort inverse ranking of an array, with removal of
    !   duplicate entries.
    !   The routine is similar to pure merge-sort ranking, but on
    !   the last pass, it sets indices in IGOEST to the rank
    !   of the value in the ordered set with duplicates removed.
    !   For performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.
    ! __________________________________________________________
    ! __________________________________________________________
    Integer(kind=IntType), intent(in) :: NVAL
    Integer(kind=IntType), Dimension (NVAL), Intent (In)  :: XDONT
    Integer(kind=IntType), Dimension (NVAL), Intent (Out) :: IGOEST
    ! __________________________________________________________
    Integer(kind=IntType) :: XTST, XDONA, XDONB
    !
    ! __________________________________________________________
    Integer(kind=IntType), allocatable, Dimension (:) :: JWRKT, IRNGT
    Integer(kind=IntType) :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2, NUNI
    Integer(kind=IntType) :: IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB

    !
    Select Case (NVAL)
    Case (:0)
       Return
    Case (1)
       IGOEST (1) = 1
       Return
    Case Default
       Continue
    End Select
    allocate(JWRKT(NVAL), IRNGT(NVAL))
    !
    !  Fill-in the index array, creating ordered couples
    !
    Do IIND = 2, NVAL, 2
       If (XDONT(IIND-1) < XDONT(IIND)) Then
          IRNGT (IIND-1) = IIND - 1
          IRNGT (IIND) = IIND
       Else
          IRNGT (IIND-1) = IIND
          IRNGT (IIND) = IIND - 1
       End If
    End Do
    If (Modulo (NVAL, 2) /= 0) Then
       IRNGT (NVAL) = NVAL
    End If
    !
    !  We will now have ordered subsets A - B - A - B - ...
    !  and merge A and B couples into     C   -   C   - ...
    !
    LMTNA = 2
    LMTNC = 4
    !
    !  First iteration. The length of the ordered subsets goes from 2 to 4
    !
    Do
       If (NVAL <= 4) Exit
       !
       !   Loop on merges of A and B into C
       !
       Do IWRKD = 0, NVAL - 1, 4
          If ((IWRKD+4) > NVAL) Then
             If ((IWRKD+2) >= NVAL) Exit
             !
             !   1 2 3
             !
             If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
             !
             !   1 3 2
             !
             If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                IRNG2 = IRNGT (IWRKD+2)
                IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                IRNGT (IWRKD+3) = IRNG2
                !
                !   3 1 2
                !
             Else
                IRNG1 = IRNGT (IWRKD+1)
                IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                IRNGT (IWRKD+2) = IRNG1
             End If
             Exit
          End If
          !
          !   1 2 3 4
          !
          If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
          !
          !   1 3 x x
          !
          If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
             IRNG2 = IRNGT (IWRKD+2)
             IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
             If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
                !   1 3 2 4
                IRNGT (IWRKD+3) = IRNG2
             Else
                !   1 3 4 2
                IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                IRNGT (IWRKD+4) = IRNG2
             End If
             !
             !   3 x x x
             !
          Else
             IRNG1 = IRNGT (IWRKD+1)
             IRNG2 = IRNGT (IWRKD+2)
             IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
             If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                IRNGT (IWRKD+2) = IRNG1
                If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
                   !   3 1 2 4
                   IRNGT (IWRKD+3) = IRNG2
                Else
                   !   3 1 4 2
                   IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                   IRNGT (IWRKD+4) = IRNG2
                End If
             Else
                !   3 4 1 2
                IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                IRNGT (IWRKD+3) = IRNG1
                IRNGT (IWRKD+4) = IRNG2
             End If
          End If
       End Do
       !
       !  The Cs become As and Bs
       !
       LMTNA = 4
       Exit
    End Do
    !
    !  Iteration loop. Each time, the length of the ordered subsets
    !  is doubled.
    !
    Do
       If (2*LMTNA >= NVAL) Exit
       IWRKF = 0
       LMTNC = 2 * LMTNC
       !
       !   Loop on merges of A and B into C
       !
       Do
          IWRK = IWRKF
          IWRKD = IWRKF + 1
          JINDA = IWRKF + LMTNA
          IWRKF = IWRKF + LMTNC
          If (IWRKF >= NVAL) Then
             If (JINDA >= NVAL) Exit
             IWRKF = NVAL
          End If
          IINDA = 1
          IINDB = JINDA + 1
          !
          !  One steps in the C subset, that we create in the final rank array
          !
          !  Make a copy of the rank array for the iteration
          !
          JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
          XDONA = XDONT (JWRKT(IINDA))
          XDONB = XDONT (IRNGT(IINDB))
          !
          Do
             IWRK = IWRK + 1
             !
             !  We still have unprocessed values in both A and B
             !
             If (XDONA > XDONB) Then
                IRNGT (IWRK) = IRNGT (IINDB)
                IINDB = IINDB + 1
                If (IINDB > IWRKF) Then
                   !  Only A still with unprocessed values
                   IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                   Exit
                End If
                XDONB = XDONT (IRNGT(IINDB))
             Else
                IRNGT (IWRK) = JWRKT (IINDA)
                IINDA = IINDA + 1
                If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                XDONA = XDONT (JWRKT(IINDA))
             End If
             !
          End Do
       End Do
       !
       !  The Cs become As and Bs
       !
       LMTNA = 2 * LMTNA
    End Do
    !
    !   Last merge of A and B into C, with removal of duplicates.
    !
    IINDA = 1
    IINDB = LMTNA + 1
    NUNI = 0
    !
    !  One steps in the C subset, that we create in the final rank array
    !
    JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
    If (IINDB <= NVAL) Then
       XTST = I_NEARLESS (Min(XDONT(JWRKT(1)), XDONT(IRNGT(IINDB))))
    Else
       XTST = I_NEARLESS (XDONT(JWRKT(1)))
    Endif
    Do IWRK = 1, NVAL
       !
       !  We still have unprocessed values in both A and B
       !
       If (IINDA <= LMTNA) Then
          If (IINDB <= NVAL) Then
             If (XDONT(JWRKT(IINDA)) > XDONT(IRNGT(IINDB))) Then
                IRNG = IRNGT (IINDB)
                IINDB = IINDB + 1
             Else
                IRNG = JWRKT (IINDA)
                IINDA = IINDA + 1
             End If
          Else
             !
             !  Only A still with unprocessed values
             !
             IRNG = JWRKT (IINDA)
             IINDA = IINDA + 1
          End If
       Else
          !
          !  Only B still with unprocessed values
          !
          IRNG = IRNGT (IWRK)
       End If
       If (XDONT(IRNG) > XTST) Then
          XTST = XDONT (IRNG)
          NUNI = NUNI + 1
       End If
       IGOEST (IRNG) = NUNI
       !
    End Do
    deallocate(JWRKT, IRNGT)
    Return
    !
  End Subroutine I_uniinv

  Function I_nearless (XVAL) result (I_nl)
    use precision
    implicit none
    !  Nearest value less than given value
    ! __________________________________________________________
    Integer(kind=intType), Intent (In) :: XVAL
    Integer(kind=intType) :: I_nl
    ! __________________________________________________________
    I_nl = XVAL - 1
    return
    !
  End Function I_nearless
#endif
end module sorting
