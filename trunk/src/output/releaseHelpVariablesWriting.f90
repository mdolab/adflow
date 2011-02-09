!
!      ******************************************************************
!      *                                                                *
!      * File:          releaseHelpVariablesWriting.f90                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-10-2005                                      *
!      * Last modified: 10-13-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine releaseHelpVariablesWriting
!
!      ******************************************************************
!      *                                                                *
!      * releaseHelpVariablesWriting releases the memory of the         *
!      * variables, which were needed to write the CGNS files.          *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use monitor
       use outputMod
       implicit none
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
       ! Release the memory of the allocatable arrays in outputMod.

       deallocate(nBlocksCGNSblock, blocksCGNSblock, stat=ierr)
       if(ierr /= 0)                                   &
         call terminate("releaseHelpVariablesWriting", &
                        "Deallocation failure for nBlocksCGNSblock, &
                        &etc.")

       if (writeGrid .and. oversetPresent) then
         deallocate(nDomPerProc, IDsBegOrAllDoms, stat=ierr)
         if(ierr /= 0)                                   &
           call terminate("releaseHelpVariablesWriting", &
                          "Deallocation failure for nDomPerProc, etc.")
       end if

       end subroutine releaseHelpVariablesWriting
