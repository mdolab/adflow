!
!      ******************************************************************
!      *                                                                *
!      * File:          allocTimeArrays.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-20-2004                                      *
!      * Last modified: 03-22-2005                                      *
!      *                                                                *
!      ******************************************************************
!
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
       use monitor
       implicit none
!
!      Subroutine argument.
!
       integer(kind=intType) :: nTimeTot
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
       ! Allocate the memory for both the time array as well as the
       ! data array.

       allocate(timeArray(nTimeTot), &
                timeDataArray(nTimeTot,nMon), stat=ierr)
       if(ierr /= 0)                       &
         call terminate("allocTimeArrays", &
                        "Memory allocation failure for timeArray &
                        &and timeDataArray")

       end subroutine allocTimeArrays
