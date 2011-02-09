!
!      ******************************************************************
!      *                                                                *
!      * File:          storeFreestreamSubface.f90                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-13-2005                                      *
!      * Last modified: 04-13-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine storeFreestreamSubface(boco)
!
!      ******************************************************************
!      *                                                                *
!      * storeFreestreamSubface stores the currently active subface in  *
!      * the array freestreamSubfaces, such that the primitive flow     *
!      * field variables can be set to the free stream values later on  *
!      * in setSupersonicInletFreeStream.                               *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCDataMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: boco
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn

       integer(kind=intType), dimension(:,:), allocatable :: tmp
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the situation we are dealing with here.

       testAllocated: if( allocated(freestreamSubfaces) ) then

         ! freestreamSubfaces has already been allocated and thus
         ! contains information. It must be reallocated and the current
         ! subface should be added.

         ! Allocate the memory for tmp and copy the data from
         ! freestreamSubfaces.

         allocate(tmp(nFreestreamSubfaces,3), stat=ierr)
         if(ierr /= 0)                              &
           call terminate("storeFreestreamSubface", &
                          "Memory allocation failure for tmp")

         tmp = freestreamSubfaces

         ! Release freestreamSubfaces and allocate it again
         ! with an increased dimension.

         deallocate(freestreamSubfaces, stat=ierr)
         if(ierr /= 0)                              &
           call terminate("storeFreestreamSubface", &
                          "Deallocation failure for freestreamSubfaces")

         nFreestreamSubfaces = nFreestreamSubfaces + 1

         allocate(freestreamSubfaces(nFreestreamSubfaces,3), &
                  stat=ierr)
         if(ierr /= 0)                              &
           call terminate("storeFreestreamSubface", &
                          "Memory allocation failure for &
                          &freestreamSubfaces")

         ! Copy the data back from tmp into freestreamSubfaces
         ! and release tmp again.

         do nn=1,(nFreestreamSubfaces-1)
           freestreamSubfaces(nn,1) = tmp(nn,1)
           freestreamSubfaces(nn,2) = tmp(nn,2)
           freestreamSubfaces(nn,3) = tmp(nn,3)
         enddo

         deallocate(tmp, stat=ierr)
         if(ierr /= 0)                              &
           call terminate("storeFreestreamSubface", &
                          "Deallocation failure for tmp")

       else testAllocated

         ! freestreamSubfaces has not been allocated yet. This is the
         ! first subface to store in this array. Allocate the array and
         ! set nFreestreamSubfaces to 1.

         nFreestreamSubfaces = 1
         allocate(freestreamSubfaces(nFreestreamSubfaces,3), &
                  stat=ierr)
         if(ierr /= 0)                              &
           call terminate("storeFreestreamSubface", &
                          "Memory allocation failure for &
                          &freestreamSubfaces")

       endif testAllocated

       ! Store the current subface in freestreamSubfaces.

       nn = nFreestreamSubfaces
       freestreamSubfaces(nn,1) = nbkLocal
       freestreamSubfaces(nn,2) = boco
       freestreamSubfaces(nn,3) = spectralSol

       end subroutine storeFreestreamSubface
