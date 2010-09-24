!
!      ******************************************************************
!      *                                                                *
!      * File:          initFlowDoms.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-24-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initFlowDoms
!
!      ******************************************************************
!      *                                                                *
!      * initFlowDoms allocates the memory for flowDoms and initializes *
!      * its pointers to null pointers, such that they do not have      *
!      * random targets.                                                *
!      *                                                                *
!      ******************************************************************
!
       use block
       use inputIteration
       use inputTimeSpectral
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, k, nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for flowDoms. Set nn to the maximum of the
       ! number of mg levels needed in the cycle and mg start level.
       ! This is namely the amount of grid levels the solver needs.

       nn = max(nMGLevels, mgStartlevel)
       allocate(flowDoms(nDom, nn, nTimeIntervalsSpectral), stat=ierr)
       if(ierr /= 0)                    &
         call terminate("initFlowDoms", &
                        "Memory allocation failure for flowDoms")

       ! Loop over all the blocks and initialize its pointers to the
       ! null-pointer.

       do k=1,nTimeIntervalsSpectral
         do j=1,nn
           do i=1,nDom
             call nullifyFlowDomPointers(i,j,k)
           enddo
         enddo
       enddo

       end subroutine initFlowDoms
