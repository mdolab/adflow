!
!     ******************************************************************
!     *                                                                *
!     * File:          adtUtils.F90                                    *
!     * Author:        Edwin van der Weide                             *
!     * Starting date: 02-09-2006                                      *
!     * Last modified: 05-18-2006                                      *
!     *                                                                *
!     ******************************************************************
!
      module adtUtils
!
!     ******************************************************************
!     *                                                                *
!     * Module, which contains small subroutines which perform useful  *
!     * tasks.                                                         *
!     *                                                                *
!     ******************************************************************
!
      use adtData
      implicit none
!
!     ******************************************************************
!     *                                                                *
!     *           Variables stored in this module.                     *
!     *                                                                *
!     ******************************************************************
!
      ! nStack:   Number of elements allocated in the stack array;
      !           needed for a more efficient implementation of the
      !           local qsort routines.
      ! stack(:): The corresponding array to store the stack.
      !           This is a pointer, such that the reallocation
      !           is easier.

      integer(kind=intType) :: nStack
      integer(kind=intType), dimension(:), pointer :: stack

      !=================================================================

      contains

        !===============================================================

        subroutine adtTerminate(jj, routineName, errorMessage)
!
!       ****************************************************************
!       *                                                              *
!       * This routine writes the given error message to standard      *
!       * output and terminates the executation of the program.        *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * routineName: Name of the routine where the error occured.    *
!       * jj:          Currently active ADT.                           *
!       *                                                              *
!       * Subroutine intent(inout) arguments.                          *
!       * -----------------------------------                          *
!       * errorMessage: On input it contains the error message to be   *
!       *               written. It is modified in this routine, such  *
!       *               that it fits on one line. On output its        *
!       *               contents is undefined, which does not matter   *
!       *               a whole lot.                                   *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments
!
        integer(kind=intType), intent(in) :: jj

        character(len=*), intent(in) :: routineName
        character(len=*), intent(in) :: errorMessage
!
!       Local parameter
!
        integer, parameter :: maxCharLine = 55
!
!       Local variables
!
        integer :: ierr, len, i2
        logical :: firstTime

        character(len=len_trim(errorMessage)) :: message
        character(len=8) :: integerString
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Copy the errorMessage into message. It is not possible to work
        ! with errorMessage directly, because it is modified in this
        ! routine. Sometimes a constant string is passed to this routine
        ! and some compilers simply fail then.

        message = errorMessage

        ! Print a nice error message. In case of a parallel executable
        ! also the processor ID is printed.

        print "(a)", "#"
        print "(a)", "#=========================== !!! Error !!! &
                      &============================"

#ifndef SEQUENTIAL_MODE
        write(integerString,"(i8)") ADTs(jj)%myID
        integerString = adjustl(integerString)

        print "(2a)", "#* adtTerminate called by processor ", &
                      trim(integerString)
#endif

        print "(2a)", "#* Run-time error in procedure ", &
                      trim(routineName)

        ! Loop to write the error message. If the message is too long it
        ! is split over several lines.

        firstTime = .true.
        do
          ! Determine the remaining error message to be written.
          ! If longer than the maximum number of characters allowed
          ! on a line, it is attempted to split the message.

           message = adjustl(message)
           len = len_trim(message)
           i2  = min(maxCharLine,len)

           if(i2 < len) i2 = index(message(:i2), " ", .true.) - 1
           if(i2 < 0)   i2 = index(message, " ") - 1
           if(i2 < 0)   i2 = len

           ! Write this part of the error message. If it is the first
           ! line of the message some additional stuff is printed.

          if( firstTime ) then
            print "(2a)", "#* Error message: ", trim(message(:i2))
            firstTime = .false.
          else
            print "(2a)", "#*                ", trim(message(:i2))
          endif

          ! Exit the loop if the entire message has been written.

          if(i2 == len) exit

          ! Adapt the string for the next part to be written.

          message = message(i2+1:)

        enddo

        ! Write the trailing message.

        print "(a)", "#*"
        print "(a)", "#* Now exiting"
        print "(a)", "#==========================================&
                     &============================"
        print "(a)", "#"

        ! Call abort and stop the program. This stop should be done in
        ! abort, but just to be sure.

        call mpi_abort(ADTs(jj)%comm, 1, ierr)
        stop

        end subroutine adtTerminate

        !***************************************************************
        !***************************************************************

        subroutine allocateADTs
!
!       ****************************************************************
!       *                                                              *
!       * This routine allocates the memory for the first ADT and is   *
!       * only called when no other ADT's are present.                 *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Local variables.
!
        integer :: ierr
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Allocate the memory for 1 ADT. Note that adtTerminate is not
        ! called when the memory allocation fails. The reason is that
        ! the processor ID for the current tree is used in this routine
        ! and that value has not been set yet.

        allocate(ADTs(1), stat=ierr)
        if(ierr /= 0) stop "Allocation failure for ADTs"

        ! Nullify the pointers of ADTs(1).

        nullify(ADTs(1)%coor)

        nullify(ADTs(1)%triaConn)
        nullify(ADTs(1)%quadsConn)
        nullify(ADTs(1)%tetraConn)
        nullify(ADTs(1)%pyraConn)
        nullify(ADTs(1)%prismsConn)
        nullify(ADTs(1)%hexaConn)

        nullify(ADTs(1)%rootLeavesProcs)
        nullify(ADTs(1)%rootBBoxes)

        nullify(ADTs(1)%elementType)
        nullify(ADTs(1)%elementID)
        nullify(ADTs(1)%xBBox)

        nullify(ADTs(1)%ADTree)

        end subroutine allocateADTs

        !***************************************************************
        !***************************************************************

        subroutine deallocateADTs(adtID)
!
!       ****************************************************************
!       *                                                              *
!       * This routine deallocates the memory for the given entry in   *
!       * the array ADTs and it tries to reallocate ADTs itself        *
!       * accordingly.                                                 *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * adtID: The entry in ADTs to be deallocated.                  *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        character(len=*), intent(in) :: adtID
!
!       Local variables.
!
        integer :: ierr

        integer(kind=intType) :: jj, nn, nAlloc, nAllocNew

        type(adtType), dimension(:), allocatable :: tmpADTs
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Determine the index in the array ADTs, which stores the given
        ! ID. As the number of trees stored is limited, a linear search
        ! algorithm is okay.

        if( allocated(ADTs) ) then
          nAlloc = ubound(ADTs, 1)
        else
          nAlloc = 0
        endif

        do jj=1,nAlloc
          if(adtID == ADTs(jj)%adtID) exit
        enddo

        ! Return immediately if the ID is not present.

        if(jj > nAlloc) return

        ! Deallocate the data for this ADT entry. Note that the memory
        ! for the nodal coordinates and the connectivity should not be
        ! deallocated, because these pointers are just set to the given
        ! input. The deallocation only takes place if the tree is active.

        if( ADTs(jj)%isActive ) then

          deallocate(ADTs(jj)%rootLeavesProcs, ADTs(jj)%rootBBoxes, &
                     ADTs(jj)%elementType,     ADTs(jj)%elementID,  &
                     ADTs(jj)%xBBox,           ADTs(jj)%ADTree,     &
                     stat=ierr)
          if(ierr /= 0)                             &
            call adtTerminate(jj, "deallocateADTs", &
                              "Deallocation failure for the ADT data")
        endif

        ! Make sure the ADT is inactive and nullify the pointers.

        ADTs(jj)%isActive = .false.

        nullify(ADTs(jj)%coor)

        nullify(ADTs(jj)%triaConn)
        nullify(ADTs(jj)%quadsConn)
        nullify(ADTs(jj)%tetraConn)
        nullify(ADTs(jj)%pyraConn)
        nullify(ADTs(jj)%prismsConn)
        nullify(ADTs(jj)%hexaConn)

        nullify(ADTs(jj)%rootLeavesProcs)
        nullify(ADTs(jj)%rootBBoxes)

        nullify(ADTs(jj)%elementType)
        nullify(ADTs(jj)%elementID)
        nullify(ADTs(jj)%xBBox)

        nullify(ADTs(jj)%ADTree)

        ! Determine the highest entry in ADTs which is still valid.

        do nn=nAlloc,1,-1
          if( ADTs(nn)%isActive ) exit
        enddo

        ! Determine the situation we are having here.

        if(nn == 0) then

          ! No active ADT's anymore. Deallocte the entire array.
          ! Note that adtTerminate cannot be called when something
          ! goes wrong.

          deallocate(ADTs, stat=ierr)
          if(ierr /= 0) stop "Deallocation failure for ADTs"

        else if(nn < nAlloc) then

          ! There are still some active ADT's, but the highest ones
          ! are inactive. Therefore ADTs is reallocated. First allocate
          ! the memory for tmpADTs to be able to retrieve the currently
          ! stored data later on.

          nAllocNew = nn
          allocate(tmpADTs(nAllocNew), stat=ierr)
          if(ierr /= 0) &
            call adtTerminate(jj, "adtDeallocateADTs", &
                              "Memory allocation failure for tmpADTs")

          ! Copy the data from ADTs to tmpADTs.

          do nn=1,nAllocNew
            tmpADTs(nn) = ADTs(nn)
          enddo

          ! Deallocate and allocate the memory for ADTs. Note that
          ! adtTerminate is not called when the memory allocation fails.
          ! The reason is that the processor ID for the current tree is
          ! used in this routine and that value may not be available
          ! anymore.

          deallocate(ADTs, stat=ierr)
          if(ierr /= 0) stop "Deallocation failure for ADTs"

          allocate(ADTs(nAllocNew), stat=ierr)
          if(ierr /= 0) stop "Allocation failure for ADTs"

          ! Copy the data back into ADTs and release the memory of
          ! tmpADTs afterwards.

          do nn=1,nAllocNew
            ADTs(nn) = tmpADTs(nn)
          enddo

          deallocate(tmpADTs, stat=ierr)
          if(ierr /= 0) stop "Deallocation failure for tmpADTs"

        endif

        end subroutine deallocateADTs

        !***************************************************************
        !***************************************************************

        subroutine qsortBBoxes(arr, nn, jj, dir)
!
!       ****************************************************************
!       *                                                              *
!       * This routine sorts the integer array arr, such that the      *
!       * coordinate of the corresponding bounding box in the          *
!       * direction dir is in increasing order. Note that the array to *
!       * store the stack is stored in this module. The reason is that *
!       * this routine is called quite often and in this way constant  *
!       * allocation, deallocation and reallocation of stack is        *
!       * avoided.                                                     *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nn:  Size of the array to be sorted.                         *
!       * jj:  The local index of the ADT from which the coordinate of *
!       *      the bounding box must be taken.                         *
!       * dir: Index of the coordinate, which must be sorted.          *
!       *                                                              *
!       * Subroutine intent(inout) arguments.                          *
!       * -----------------------------------                          *
!       * arr(:): On input it contains the bounding box ID's which     *
!       *         must be sorted. On output these ID's are sorted,     *
!       *         such that the given coordinate is in increasing      *
!       *         order.                                               *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=intType), intent(in)  :: nn, jj, dir

        integer(kind=intType), dimension(:), intent(inout) :: arr
!
!       Local parameters.
!
        integer(kind=intType), parameter :: m = 7
!
!       Local variables.
!
        integer(kind=intType) :: i, j, k, r, l, jStack
        integer(kind=intType) :: a, tmp

        real(kind=realType) :: ra
        real(kind=realType), dimension(:,:), pointer :: xBBox
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Set the pointer for the coordinates of the bounding boxes.

        xBBox => ADTs(jj)%xBBox

        ! Initialize the variables that control the sorting.

        jStack = 0
        l      = 1
        r      = nn

        ! Start of the algorithm.

        sortLoop: do

          ! Check for the size of the subarray.

          testInsertion: if((r-l) < m) then

            ! Perform the insertion sort.

            do j=(l+1),r
              a  = arr(j)
              ra = xBBox(dir,a)
              do i=(j-1),l,-1
                if(xBBox(dir,arr(i)) <= ra) exit
                arr(i+1) = arr(i)
              enddo
              arr(i+1) = a
            enddo

            ! In case there are no more elements on the stack, exit from
            ! the outermost do-loop. Algorithm has finished.

            if(jStack == 0) exit sortLoop

            ! Pop stack and begin a new round of partitioning.

            r = stack(jStack)
            l = stack(jStack-1)
            jStack = jStack - 2

          else testInsertion

            ! Subarray is larger than the threshold for a linear sort.
            ! Choose median of left, center and right elements as
            ! partitioning element a. Also rearrange so that 
            ! (l) <= (l+1) <= (r).

            k = (l+r)/2
            tmp      = arr(k)      ! Swap the elements
            arr(k)   = arr(l+1)    ! k and l+1.
            arr(l+1) = tmp

            if(xBBox(dir,arr(r)) < xBBox(dir,arr(l))) then
              tmp    = arr(l)             ! Swap the elements
              arr(l) = arr(r)             ! r and l.
              arr(r) = tmp
            endif

            if(xBBox(dir,arr(r)) < xBBox(dir,arr(l+1))) then
              tmp      = arr(l+1)         ! Swap the elements
              arr(l+1) = arr(r)           ! r and l+1.
              arr(r)   = tmp
            endif

            if(xBBox(dir,arr(l+1)) < xBBox(dir,arr(l))) then
              tmp      = arr(l+1)         ! Swap the elements
              arr(l+1) = arr(l)           ! l and l+1.
              arr(l)   = tmp
            endif

            ! Initialize the pointers for partitioning.

            i  = l+1
            j  = r
            a  = arr(l+1)
            ra = xBBox(dir,a)

            ! The innermost loop.

            innerLoop: do

              ! Scan up to find element >= a.
              do
                i = i+1
                if(ra <= xBBox(dir,arr(i))) exit
              enddo

              ! Scan down to find element <= a.
              do
                j = j-1
                if(xBBox(dir,arr(j)) <= ra) exit
              enddo

              ! Exit the loop in case the pointers i and j crossed.

              if(j < i) exit innerLoop

              ! Swap the element i and j.

              tmp    = arr(i)
              arr(i) = arr(j)
              arr(j) = tmp

            enddo innerLoop

            ! Swap the entries j and l+1. Remember that a equals
            ! arr(l+1).

            arr(l+1) = arr(j)
            arr(j)   = a

            ! Push pointers to larger subarray on stack; process smaller
            ! subarray immediately. Check if enough memory is available.
            ! If not reallocate it.

            jStack = jStack + 2

            if(jStack > nStack) call reallocPlus(stack, nStack, 100, jj)

            if((r-i+1) >= (j-l)) then
              stack(jStack)   = r
              r               = j-1
              stack(jStack-1) = j
            else
              stack(jStack)   = j-1
              stack(jStack-1) = l
              l               = j
            endif

          endif testInsertion
        enddo sortLoop

        ! Check in debug mode if the sort has been done correctly.

        if( debug ) then
          do i=1,(nn-1)
            if(xBBox(dir,arr(i+1)) < xBBox(dir,arr(i))) then
              call adtTerminate(jj, "qsortBBoxes", &
                                "Array is not sorted correctly")
            endif
          enddo
        endif

        end subroutine qsortBBoxes

        !***************************************************************
        !***************************************************************

        subroutine qsortBBoxTargets(arr, nn, jj)
!
!       ****************************************************************
!       *                                                              *
!       * This routine sorts the given number of bounding box targets  *
!       * in increasing order, based on the generalized < operator.    *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments
!
        integer(kind=intType), intent(in) :: nn, jj

        type(adtBBoxTargetType), dimension(:), pointer :: arr
!
!       Local variables
!
        integer(kind=intType), parameter :: m = 7

        integer(kind=intType) :: i, j, k, r, l, jStack

        type(adtBBoxTargetType) :: a, tmp
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution                                              *
!       *                                                              *
!       ****************************************************************
!
        ! Initialize the variables that control the sorting.

        jStack = 0
        l      = 1
        r      = nn

        ! Start of the algorithm

        sortLoop: do

          ! Check for the size of the subarray.

          testInsertion: if((r-l) < m) then

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

            if(jStack == 0) exit sortLoop

            ! Pop stack and begin a new round of partitioning.

            r = stack(jStack)
            l = stack(jStack-1)
            jStack = jStack - 2

          else testInsertion

            ! Subarray is larger than the threshold for a linear sort.
            ! Choose median of left, center and right elements as
            ! partitioning element a. Also rearrange so that
            ! (l) <= (l+1) <= (r).

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

            innerLoop: do

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

              if(j < i) exit innerLoop

              ! Swap the element i and j.

              tmp    = arr(i)
              arr(i) = arr(j)
              arr(j) = tmp

            enddo innerLoop

            ! Swap the entries j and l+1. Remember that a equals
            ! arr(l+1).

            arr(l+1) = arr(j)
            arr(j)   = a

            ! Push pointers to larger subarray on stack; process smaller
            ! subarray immediately. Check if enough memory is available.
            ! If not reallocate it.

            jStack = jStack + 2

            if(jStack > nStack) call reallocPlus(stack, nStack, 100, jj)

            if((r-i+1) >= (j-l)) then
              stack(jStack)   = r
              r               = j-1
              stack(jStack-1) = j
            else
              stack(jStack)   = j-1
              stack(jStack-1) = l
              l               = j
            endif

          endif testInsertion
        enddo sortLoop

        ! Check in debug mode whether the array is really sorted.

        if( debug ) then
          do i=1,(nn-1)
            if(arr(i+1) < arr(i))                       &
              call adtTerminate(jj, "qsortBBoxTargets", &
                                "Array is not sorted correctly")
          enddo
        endif

        end subroutine qsortBBoxTargets

        !***************************************************************
        !***************************************************************

        subroutine reallocateADTs(adtID, jj)
!
!       ****************************************************************
!       *                                                              *
!       * This routine reallocates the memory for the ADTs array, such *
!       * that it is possible to store a new ADT. First it is tried to *
!       * find an empty spot in the currently allocated array. If this *
!       * is not present a true reallocation takes place.              *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * adtID: The ID of the ADT.                                    *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * jj: The index in the array ADTs, where the new entry will be *
!       *     stored.                                                  *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        character(len=*), intent(in)          :: adtID
        integer(kind=intType), intent(out) :: jj
!
!       Local variables.
!
        integer :: ierr

        integer(kind=intType) :: nn, nAlloc

        type(adtType), dimension(:), allocatable :: tmpADTs
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Determine the current size of ADTs and look for an empty spot
        ! in the currently allocated array. Also check if an ADT with
        ! the given ID is not already active. A linear search algorithm
        ! is used, because the number of ADT's stored is limited.

        nAlloc = ubound(ADTs, 1)
        jj = nAlloc + 1
        do nn=1,nAlloc
          if( ADTs(nn)%isActive ) then
            if(adtID == ADTs(nn)%adtID) exit
          else if(jj > nAlloc) then
            jj = nn
          endif
        enddo

        ! If the given ID corresponds to an already active tree,
        ! terminate. To avoid a messy output only processor 0 prints
        ! an error message while the other ones wait to get killed.

        if(nn <= nAlloc) then
          if(ADTs(nn)%myID == 0)                    &
            call adtTerminate(nn, "reallocateADTs", &
                              "Given ID corresponds to an already &
                              &active ADT")
          call mpi_barrier(ADTs(nn)%comm, ierr)
        endif

        ! Check if a reallocate must be done.

        checkReallocate: if(jj > nAlloc) then

          ! No empty spot present in ADTs. A true reallocation must be
          ! performed. First allocate the memory for tmpADTs to be able
          ! to retrieve the currently stored data later on. Note that
          ! adtTerminate is not called when the memory allocation fails.
          ! The reason is that the processor ID for the current tree is
          ! used in this routine and that value has not been set yet.

          allocate(tmpADTs(nAlloc), stat=ierr)
          if(ierr /= 0) stop "Allocation failure for tmpADTs"

          ! Copy the data from ADTs to tmpADTs.

          do nn=1,nAlloc
            tmpADTs(nn) = ADTs(nn)
          enddo

          ! Release the memory of ADTs and allocate it again with
          ! increased size.

          deallocate(ADTs, stat=ierr)
          if(ierr /= 0) stop "Deallocation failure for ADTs"

          allocate(ADTs(jj), stat=ierr)
          if(ierr /= 0) stop "Allocation failure for ADTs"

          ! Copy the data back from tmpADTs.

          do nn=1,nAlloc
            ADTs(nn) = tmpADTs(nn)
          enddo

          ! Release the memory of tmpADTs.

          deallocate(tmpADTs, stat=ierr)
          if(ierr /= 0) stop "Deallocation failure for tmpADTs"

          ! Nullify the pointers of the new entry.

          nullify(ADTs(jj)%coor)

          nullify(ADTs(jj)%triaConn)
          nullify(ADTs(jj)%quadsConn)
          nullify(ADTs(jj)%tetraConn)
          nullify(ADTs(jj)%pyraConn)
          nullify(ADTs(jj)%prismsConn)
          nullify(ADTs(jj)%hexaConn)

          nullify(ADTs(jj)%rootLeavesProcs)
          nullify(ADTs(jj)%rootBBoxes)

          nullify(ADTs(jj)%elementType)
          nullify(ADTs(jj)%elementID)
          nullify(ADTs(jj)%xBBox)

          nullify(ADTs(jj)%ADTree)

        endif checkReallocate

        end subroutine reallocateADTs

        !***************************************************************
        !***************************************************************

        subroutine reallocBBoxTargetTypePlus(arr, nSize, nInc, jj)
!
!       ****************************************************************
!       *                                                              *
!       * This routine reallocates the memory of the given             *
!       * adtBBoxTargetType pointer array.                             *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * jj:   Entry in the array ADTs, whose ADT is being searched.  *
!       * nInc: Increment of the size of the array.                    *
!       *                                                              *
!       * Subroutine intent(inout) arguments.                          *
!       * -----------------------------------                          *
!       * nSize: On input it contains the size of the given array.     *
!       *        On output this value is incremented by nInc.          *
!       *                                                              *
!       * Subroutine pointer arguments.                                *
!       * -----------------------------                                *
!       * arr: Array to be reallocated.                                *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer,                  intent(in)    :: nInc
        integer(kind=intType), intent(in)    :: jj
        integer(kind=intType), intent(inout) :: nSize

        type(adtBBoxTargetType), dimension(:), pointer :: arr
!
!       Local variables.
!
        integer :: ierr
        integer(kind=intType) :: i, nOld

        type(adtBBoxTargetType), dimension(:), pointer :: tmp
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Store the input value of nSize and set the pointer tmp to the
        ! original array.

        nOld = nSize
        tmp => arr

        ! Allocate the new memory for the array.

        nSize = nSize + nInc
        allocate(arr(nSize), stat=ierr)
        if(ierr /= 0) &
          call adtTerminate(jj, "reallocBBoxTargetTypePlus", &
                            "Memory allocation failure for arr.")

        ! Copy the data from the original array into arr.

        nOld = min(nOld,nSize)
          do i=1,nOld
          arr(i) = tmp(i)
        enddo

        ! Release the memory of tmp, which points to the original
        ! memory of the given array.

        deallocate(tmp, stat=ierr)
        if(ierr /= 0) &
          call adtTerminate(jj, "reallocBBoxTargetTypePlus", &
                            "Deallocation failure for tmp.")

        end subroutine reallocBBoxTargetTypePlus

        !***************************************************************
        !***************************************************************

        subroutine reallocPlus(arr, nSize, nInc, jj)
!
!       ****************************************************************
!       *                                                              *
!       * This internal routine reallocates the memory of the given    *
!       * pointer array.                                               *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * jj:   Entry in the array ADTs, whose ADT is being searched.  *
!       * nInc: Increment of the size of the array.                    *
!       *                                                              *
!       * Subroutine intent(inout) arguments.                          *
!       * -----------------------------------                          *
!       * nSize: On input it contains the size of the given array.     *
!       *        On output this value is incremented by nInc.          *
!       *                                                              *
!       * Subroutine pointer arguments.                                *
!       * -----------------------------                                *
!       * arr: Array to be reallocated.                                *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer,                  intent(in)    :: nInc
        integer(kind=intType), intent(in)    :: jj
        integer(kind=intType), intent(inout) :: nSize

        integer(kind=intType), dimension(:), pointer :: arr
!
!       Local variables.
!
        integer :: ierr
        integer(kind=intType) :: i, nOld

        integer(kind=intType), dimension(:), pointer :: tmp
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Store the input value of nSize and set the pointer tmp to the
        ! original array.

        nOld = nSize
        tmp => arr

        ! Allocate the new memory for the array.

        nSize = nSize + nInc
        allocate(arr(nSize), stat=ierr)
        if(ierr /= 0)                          &
          call adtTerminate(jj, "reallocPlus", &
                            "Memory allocation failure for arr.")

        ! Copy the data from the original array into arr.

        nOld = min(nOld,nSize)
        do i=1,nOld
          arr(i) = tmp(i)
        enddo

        ! Release the memory of tmp, which points to the original
        ! memory of the given array.

        deallocate(tmp, stat=ierr)
        if(ierr /= 0)                          &
          call adtTerminate(jj, "reallocPlus", &
                            "Deallocation failure for tmp.")

        end subroutine reallocPlus

      end module adtUtils
