       subroutine allocTimeArrays(nTimeTot)
!
!      ******************************************************************
!      *                                                                *
!      * allocTimeArrays allocates the memory for the arrays to store   *
!      * the time history of the unsteady computation. The number of    *
!      * time steps specified is enought to store the total number of   *
!      * time steps of the current computation plus possible earlier    *
!      * computations.                                                  *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use monitor, only : timeArray, timeDataArray, nMon
       use utils, only : terminate
       implicit none
!
!      Subroutine argument.
!
       integer(kind=intType) :: nTimeTot
!
!      Local variables.
!
       integer :: ierr

       ! Allocate the memory for both the time array as well as the
       ! data array.

       if (allocated(timeArray)) then 
          deallocate(timeArray) 
       end if
       if (allocated(timeDataArray)) then
          deallocate(timeDataArray) 
       end if
          
       allocate(timeArray(nTimeTot), &
                timeDataArray(nTimeTot,nMon), stat=ierr)
       if(ierr /= 0)                       &
         call terminate("allocTimeArrays", &
                        "Memory allocation failure for timeArray &
                        &and timeDataArray")
       
       end subroutine allocTimeArrays
