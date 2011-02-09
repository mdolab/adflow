!
!      ******************************************************************
!      *                                                                *
!      * File:          reallocateReals.F90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-18-2003                                      *
!      * Last modified: 05-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine reallocateReal(realArray, newSize, oldSize, &
                                  alwaysFreeMem)
!
!      ******************************************************************
!      *                                                                *
!      * ReallocateReal reallocates the given real array to the given   *
!      * new size. The old values of the array are copied. Note that    *
!      * newSize can be both smaller and larger than oldSize.           *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(:), pointer :: realArray
       integer(kind=intType), intent(in) :: newSize, oldSize
       logical, intent(in) :: alwaysFreeMem
!
!      Local variables.
!
       real(kind=realType), dimension(:), pointer :: tmp

       integer(kind=intType) :: i, nn

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

       ! Set the pointer for tmp to realArray.

       tmp => realArray

       ! Allocate the memory for realArray in case newSize is larger
       ! than 0 or if alwaysFreeMem is .True. And copy the old data
       ! into it.

       if(newSize > 0 .or. alwaysFreeMem) then
         allocate(realArray(newSize), stat=ierr)
         if(ierr /= 0)                       &
           call terminate("reallocateReal", &
                          "Memory allocation failure for realArray")
         do i=1,nn
           realArray(i) = tmp(i)
         enddo
       endif

       ! Release the memory for tmp in case oldSize is larger than 0 or
       ! if alwaysFreeMem is .True.

       if(oldSize > 0 .or. alwaysFreeMem) then
         deallocate(tmp, stat=ierr)
         if(ierr /= 0)                       &
           call terminate("reallocateReal", &
                          "Deallocation error for tmp")
       endif

       end subroutine reallocateReal

       !================================================================

       subroutine reallocateReal2(realArray, newSize1, newSize2, &
                                   oldSize1, oldSize2,            &
                                   alwaysFreeMem)
!
!      ******************************************************************
!      *                                                                *
!      * ReallocateReal2 reallocates the given 2d integer array to      *
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
       real(kind=realType), dimension(:,:), pointer :: realArray
       integer(kind=intType), intent(in) :: newSize1, newSize2, &
                                             oldSize1, oldSize2
       logical, intent(in) :: alwaysFreeMem
!
!      Local variables.
!
       real(kind=realType), dimension(:,:), pointer :: tmp

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

       tmp => realArray

       ! Allocate the memory for realArray in case newSize is larger
       ! than 0 or if alwaysFreeMem is .True. And copy the old data
       ! into it.

       if(newSize > 0 .or. alwaysFreeMem) then
         allocate(realArray(newSize1,newSize2), stat=ierr)
         if(ierr /= 0)                           &
           call terminate("reallocateReal2", &
                          "Memory allocation failure for realArray")
         do j=1,nn2
           do i=1,nn1
             realArray(i,j) = tmp(i,j)
           enddo
         enddo
       endif

       ! Release the memory of tmp in case oldSize is larger than 0
       ! or if alwaysFreeMem is .True..

       if(oldSize > 0 .or. alwaysFreeMem) then
         deallocate(tmp, stat=ierr)
         if(ierr /= 0)                           &
           call terminate("reallocateReal2", &
                          "Deallocation error for tmp")
       endif

       end subroutine reallocateReal2
