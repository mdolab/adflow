
!
!      ******************************************************************
!      *                                                                *
!      * File:          storeReferenceMesh.f90                          *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 12-11-2008                                      *
!      * Last modified: 12-11-2008                                      *
!      *                                                                *
!      ******************************************************************
!

subroutine storeReferenceMesh

use blockPointers
implicit none

!
! Local Variables
!

integer(kind=intType)::i,j,k,nn

!
!Begin Execution
!

do nn=1,nDom
   call setPointers(nn,1,1)
   do i = 0,ie
      do j=0,je
         do k=0,ke
            xInit(i,j,k,:) = x(i,j,k,:)
         enddo
      enddo
   enddo
enddo


end subroutine storeReferenceMesh
