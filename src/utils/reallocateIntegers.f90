!
!      ******************************************************************
!      *                                                                *
!      * File:          reallocateIntegers.F90                          *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 12-23-2002                                      *
!      * Last modified: 06-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine reallocateInteger(intArray, newSize, oldSize, &
                                    alwaysFreeMem)
!
!      ******************************************************************
!      *                                                                *
!      * reallocateInteger reallocates the given integer array to the   *
!      * given new size. The old values of the array are copied. Note   *
!      * that newSize can be both smaller and larger than oldSize.      *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), dimension(:), pointer :: intArray
       integer(kind=intType), intent(in) :: newSize, oldSize
       logical, intent(in) :: alwaysFreeMem
!
!      Local variables.
!
       integer(kind=intType), dimension(:), pointer :: tmp

       integer(kind=intType) :: i, nn, ll

       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the minimum of newSize and oldSize.

       nn = min(newSize, oldSize)

       ! Set the pointer for tmp to intArray.

       tmp => intArray

       ! Allocate the memory for intArray in case newSize is larger
       ! than 0 or if alwaysFreeMem is .true. And copy the old data
       ! into it. Preserve the lower bound.

       if(newSize > 0 .or. alwaysFreeMem) then

         ll = 1
         if (associated(intArray)) ll = lbound(intArray,1)

         allocate(intArray(ll:newSize+ll-1), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("reallocateInteger", &
                          "Memory allocation failure for intArray")
         do i=ll,ll+nn-1
           intArray(i) = tmp(i)
         enddo
       endif

       ! Release the memory for tmp in case oldSize is larger than 0 or
       ! if alwaysFreeMem is .true.

       if(oldSize > 0 .or. alwaysFreeMem) then
         deallocate(tmp, stat=ierr)
         if(ierr /= 0)                         &
           call terminate("reallocateInteger", &
                          "Deallocation error for tmp")
       endif

       end subroutine reallocateInteger

       !================================================================

       subroutine reallocateMpiOffsetKindInteger(intArray, newSize, &
                                                 oldSize, alwaysFreeMem)
!
!      ******************************************************************
!      *                                                                *
!      * reallocateMpiOffsetKindInteger reallocates the given           *
!      * mpi_offset_kind integer array to the given new size. The old   *
!      * values of the array are copied. Note that newSize can be both  *
!      * smaller and larger than oldSize.                               *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=mpi_offset_kind), dimension(:), pointer :: intArray
       integer(kind=intType), intent(in) :: newSize, oldSize
       logical, intent(in) :: alwaysFreeMem
!
!      Local variables.
!
       integer(kind=mpi_offset_kind), dimension(:), pointer :: tmp

       integer(kind=intType) :: i, nn, ll

       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the minimum of newSize and oldSize.

       nn = min(newSize, oldSize)

       ! Set the pointer for tmp to intArray.

       tmp => intArray

       ! Allocate the memory for intArray in case newSize is larger
       ! than 0 or if alwaysFreeMem is .true. And copy the old data
       ! into it. Preserve the lower bound.

       if(newSize > 0 .or. alwaysFreeMem) then

         ll = 1
         if (associated(intArray)) ll = lbound(intArray,1)

         allocate(intArray(ll:newSize+ll-1), stat=ierr)
         if(ierr /= 0)                                      &
           call terminate("reallocateMpiOffsetKindInteger", &
                          "Memory allocation failure for intArray")
         do i=ll,ll+nn-1
           intArray(i) = tmp(i)
         enddo
       endif

       ! Release the memory for tmp in case oldSize is larger than 0 or
       ! if alwaysFreeMem is .true.

       if(oldSize > 0 .or. alwaysFreeMem) then
         deallocate(tmp, stat=ierr)
         if(ierr /= 0)                                      &
           call terminate("reallocateMpiOffsetKindInteger", &
                          "Deallocation error for tmp")
       endif

       end subroutine reallocateMpiOffsetKindInteger

       !================================================================

       subroutine reallocateInteger2(intArray, newSize1, newSize2, &
                                     oldSize1, oldSize2,           &
                                     alwaysFreeMem)
!
!      ******************************************************************
!      *                                                                *
!      * reallocateInteger2 reallocates the given 2D integer array to   *
!      * the given new sizes. The old values of the array are copied.   *
!      * Note that the newSizes can be both smaller and larger than     *
!      * the oldSizes.                                                  *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), dimension(:,:), pointer :: intArray
       integer(kind=intType), intent(in) :: newSize1, newSize2, &
                                            oldSize1, oldSize2
       logical, intent(in) :: alwaysFreeMem
!
!      Local variables.
!
       integer(kind=intType), dimension(:,:), pointer :: tmp

       integer(kind=intType) :: newSize, oldSize
       integer(kind=intType) :: nn1, nn2, nn

       integer(kind=intType) :: i, j

       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the total new and old size.

       newSize = newSize1*newSize2
       oldSize = oldSize1*oldSize2

       ! Determine for each of the 2 components the minimum of the new
       ! and the old size. Multiply these values to obtain the total
       ! amount of data that must be copied.

       nn1 = min(newSize1, oldSize1)
       nn2 = min(newSize2, oldSize2)

       nn = nn1*nn2

       ! Set the pointer for tmp.

       tmp => intArray

       ! Allocate the memory for intArray in case newSize is larger
       ! than 0 or if alwaysFreeMem is .true. and copy the old data
       ! into it.

       if(newSize > 0 .or. alwaysFreeMem) then
         allocate(intArray(newSize1,newSize2), stat=ierr)
         if(ierr /= 0)                          &
           call terminate("reallocateInteger2", &
                          "Memory allocation failure for intArray")
         do j=1,nn2
           do i=1,nn1
             intArray(i,j) = tmp(i,j)
           enddo
         enddo
       endif

       ! Release the memory of tmp in case oldSize is larger than 0
       ! or if alwaysFreeMem is .true..

       if(oldSize > 0 .or. alwaysFreeMem) then
         deallocate(tmp, stat=ierr)
         if(ierr /= 0)                          &
           call terminate("reallocateInteger2", &
                          "Deallocation error for tmp")
       endif

       end subroutine reallocateInteger2
