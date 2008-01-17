!
!      ******************************************************************
!      *                                                                *
!      * File:          checkSizeBoundaryList.f90                       *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 08-25-2005                                      *
!      * Last modified: 08-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine checkSizeBoundaryList(nStored, nAdd)
!
!      ******************************************************************
!      *                                                                *
!      * checkSizeBoundaryList checks to make sure that the size of the *
!      * of the oversetHalo boundary list is sufficient enough to store *
!      * all required data. The current portion of the array being used *
!      * is given by nStored (may be less than the actual size) and the *
!      * amount to be added is given by nAdd. A reallocation process is *
!      * performed if the required size is > the current size.          *
!      *                                                                *
!      ******************************************************************
!
       use boundaryLIst
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nStored, nAdd
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, currentSize, requiredSize, newSize

       type(haloListType), dimension(:), allocatable :: tmp
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the current and required size of the list and return if
       ! no reallocation is necessary.

       currentSize  = ubound(oversetHalo,1)
       requiredSize = nStored + nAdd

       if (requiredSize <= currentSize) return

       ! A reallocation of the list is needed. First allocate a temporary
       ! list for the current data and copy all current data to it.

       allocate(tmp(nStored), stat=ierr)
       if (ierr /= 0)                            &
         call terminate("checkSizeBoundaryList", &
                        "Memory allocation failure for tmp")

       do i=1,nStored
         allocate(tmp(i)%interp(3), stat=ierr)
         if (ierr /= 0)                            &
           call terminate("checkSizeBoundaryList", &
                          "Memory allocation failure for interp")

         tmp(i) = oversetHalo(i)
       end do

       ! Deallocate the current boundary list.

       do i=1,currentSize
         deallocate(oversetHalo(i)%interp, stat=ierr)
         if (ierr /= 0)                            &
           call terminate("checkSizeBoundaryList", &
                          "Deallocation failure for interp")
       end do

       deallocate(oversetHalo, stat=ierr)
       if (ierr /= 0)                            &
         call terminate("checkSizeBoundaryList", &
                        "Deallocation failure for oversetHalo")

       ! Set the new size of the boundary list and reallocate. Note the
       ! new size is at least that required but could be more to avoid
       ! lots of reallocation.

       newSize = nStored + max(nAdd, 100_intType)

       allocate(oversetHalo(newSize), stat=ierr)
       if (ierr /= 0)                            &
         call terminate("checkSizeBoundaryList", &
                        "Memory allocation failure for oversetHalo")

       do i=1,newSize
         allocate(oversetHalo(i)%interp(3), stat=ierr)
         if (ierr /= 0)                            &
           call terminate("checkSizeBoundaryList", &
                          "Memory allocation failure for interp")
       end do

       ! Copy the data back from the temporary array.

       do i=1,nStored
         oversetHalo(i) = tmp(i)
       end do

       ! Deallocate the temporary list.

       do i=1,nStored
         deallocate(tmp(i)%interp, stat=ierr)
         if (ierr /= 0)                            &
           call terminate("checkSizeBoundaryList", &
                          "Deallocation failure for interp")
       end do

       deallocate(tmp, stat=ierr)
       if (ierr /= 0)                            &
         call terminate("checkSizeBoundaryList", &
                        "Deallocation failure for tmp")

       end subroutine checkSizeBoundaryList
