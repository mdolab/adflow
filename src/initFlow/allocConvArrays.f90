!
!      ******************************************************************
!      *                                                                *
!      * File:          allocConvArrays.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-23-2003                                      *
!      * Last modified: 07-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine allocConvArrays(nIterTot)
!
!      ******************************************************************
!      *                                                                *
!      * allocConvArrays allocates the memory for the convergence       *
!      * arrays. The number of iterations allocated, nIterTot, is       *
!      * enough to store the maximum number of iterations specified     *
!      * plus possible earlier iterations read from the restart file.   *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use inputIO
       use inputTimeSpectral
       use monitor
       implicit none
!
!      Subroutine argument.
!
       integer(kind=intType) :: nIterTot
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
       ! Return immediately if the convergence history (of the inner
       ! iterations) does not need to be stored. This logical can
       ! only be .false. for an unsteady computation.

       if(.not. storeConvInnerIter) return 

       ! Allocate the memory for convArray and initialize them,
       ! just to be sure.

       allocate(convArray(0:nIterTot,nTimeIntervalsSpectral,nMon), &
                stat=ierr)
       if(ierr /= 0)                         &
         call terminate("allocConvArrays", &
                        "Memory allocation failure for convArray")

       convArray = zero

       end subroutine allocConvArrays
