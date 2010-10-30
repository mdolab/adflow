!
!      ******************************************************************
!      *                                                                *
!      * File:          mdFamilyNames.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-05-2004                                      *
!      * Last modified: 03-24-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mdGetFamilyNames
!
!      ******************************************************************
!      *                                                                *
!      * mdGetFamilyNames allocates the memory for mdFamilyNames        *
!      * and stores the family names in it.                             *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use mdData
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! If the memory of mdFamilyNames is already allocated, release it.

       call mdDeleteFamilyNames

       ! Set the number of families to all the families present in the
       ! grid, even if they do not have a solid wall boundary condition.
       ! Allocate the memory for mdFamilyNames.

       mdNFamilies = cgnsNfamilies
       allocate(mdFamilyNames(mdNFamilies), stat=ierr)
       if(ierr /= 0)                        &
         call terminate("mdGetFamilyNames", &
                        "Memory allocation failure for mdFamilyNames")

       ! And copy the family names in mdFamilyNames.

       do nn=1,mdNFamilies
         mdFamilyNames(nn) = cgnsFamilies(nn)%familyName
       enddo

       end subroutine mdGetFamilyNames

!      ==================================================================

       subroutine mdDeleteFamilyNames
!
!      ******************************************************************
!      *                                                                *
!      * mdDeleteFamilyNames releases the memory of mdNFamilies.        *
!      *                                                                *
!      ******************************************************************
!
       use mdData
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
       ! Release the memory of mdNFamilies; only if it is allocated.

       if( allocated(mdFamilyNames) ) then
         deallocate(mdFamilyNames, stat=ierr)
         if(ierr /= 0)                           &
           call terminate("mdDeleteFamilyNames", &
                          "Deallocation failure for mdFamilyNames")
       endif

       ! Set mdNFamilies to 0.

       mdNFamilies = 0

       end subroutine mdDeleteFamilyNames
