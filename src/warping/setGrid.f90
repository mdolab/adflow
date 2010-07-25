!
! ***********************************
! *  File: setGrid.f90
! *  Author: Gaetan Kenway
! *  Started: 05-29-2010
! *  Modified: 05-29-2010
! ***********************************

subroutine setGrid(cgnsdof,ndofcgns)

! The purpose of this routine is to take the flattened grid of points
! from the cgns grid and perform a VecScatter Operation to set the
! grid in the sumb block format

use precision 
use blockPointers
use communication
use warpingPetsc
implicit none

integer(kind=intType),intent(in) :: ndofcgns
real(kind=realType) ,intent(in) :: cgnsdof(ndofcgns)

! Local Variables

integer(kind=intType) :: nn,ii,i,j,k,ierr
integer(kind=intType) :: cgnsindices(ndofcgns)
integer(kind=intType) :: lowInd,HighInd 

! Set the cgnsdof into the cgnsGridVec 
call VecGetOwnershipRange(cgnsgridvec,lowInd,highInd,ierr)

ii = 0
do i=lowInd,highInd-1
   ii = ii + 1
   cgnsindices(ii) = i
end do

call VecSetValues(cgnsGridVec,ndofcgns,cgnsindices,cgnsdof,INSERT_VALUES,ierr)

call VecAssemblyBegin(cgnsGridVec,ierr)
call VecAssemblyEnd(cgnsGridVec,ierr)

call VecScatterBegin(cgnsTOsumbGrid,cgnsGridVec,sumbGridVec,INSERT_VALUES,SCATTER_FORWARD,ierr)
call VecScatterEnd  (cgnsTOsumbGrid,cgnsGridVec,sumbGridVec,INSERT_VALUES,SCATTER_FORWARD,ierr)

! Now simply loop over the domains and set the values

call VecGetOwnershipRange(sumbGridVec,lowInd,HighInd,ierr)

ii = 0
do nn=1,nDom
   call setPointers(nn,1_intType,1_intType)

   do k=1,kl
      do j=1,jl
         do i=1,il
            call VecGetValues(sumbGridVec,1,lowInd+ii  , X(i,j,k,1),ierr)
            call VecGetValues(sumbGridVec,1,lowInd+ii+1, X(i,j,k,2),ierr)
            call VecGetValues(sumbGridVec,1,lowInd+ii+2, X(i,j,k,3),ierr)
            ii = ii + 3
         end do
      end do
   end do
end do

call xhalo(1_intType)
end subroutine setGrid


