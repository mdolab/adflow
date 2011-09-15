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
       real(kind=realType),dimension(nTimeIntervalsSpectral)::convTemp
       logical :: storingPrev=.false.
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
       if( allocated(convArray)) then
          convTemp(:) = convArray(0,:,1)
          storingPrev = .True.
       end if
       
       if( allocated(convArray)) call deallocConvArrays

       ! Allocate the memory for convArray and initialize them,
       ! just to be sure.

       allocate(convArray(0:nIterTot,nTimeIntervalsSpectral,nMon), &
                stat=ierr)
       if(ierr /= 0)                         &
         call terminate("allocConvArrays", &
                        "Memory allocation failure for convArray")

       convArray = zero
       if (storingPrev)then
          convArray(0,:,1)=convTemp(:)
       end if

       end subroutine allocConvArrays

       subroutine deallocConvArrays
!
!      ******************************************************************
!      *                                                                *
!      * deallocConvArrays deallocates the memory for the convergence   *
!      * arrays.                                                        *
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

       deallocate(convArray,stat=ierr)
       
       if(ierr /= 0)                         &
         call terminate("deallocConvArrays", &
                        "Memory deallocation failure for convArray")


     end subroutine deallocConvArrays
