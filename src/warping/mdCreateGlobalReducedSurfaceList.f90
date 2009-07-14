!
! ***********************************
! *  File: generateReducedSurfaceList.f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 07-13-2009
! *  Modified: 07-13-2009
! ***********************************

subroutine mdCreateGlobalReducedSurfaceList

use blockPointers
use mdData           !mdNSurfNodesCompact,mdGlobalSurfCoords
use mdDataLocal    
use communication !myid

implicit none

!Subroutine Arguments



!Local Variables
integer(kind=intType) :: startInd, endInd,famID
integer(kind=intType)::level=1,sps=1,nn,mm,i,j,k,ierr


!Begin Execution
if(.not. allocated(mdGlobalSurfxx))then
   allocate(mdGlobalSurfxx(3,mdNSurfNodesCompact))
endif
!initialize coords to -inf
mdGlobalSurfxx = -large

!loop over reduced list pulling the values from the global list
do i = 1,mdNSurfNodesCompact
   do j = 1,size(mdSurfGlobalIndLocal(5,:))
      if (int(mdSurfGlobalIndLocal(5,j))==i-1)then
         do k = 1,3
            call setPointers(mdSurfGlobalIndLocal(4,j),level,sps)
            mdGlobalSurfxx(k,i) = x(mdSurfGlobalIndLocal(1,j),mdSurfGlobalIndLocal(2,j),mdSurfGlobalIndLocal(3,j),k)
         end do
         
      endif
   enddo
   call mpi_barrier(sumb_comm_world, ierr)
   do k = 1,3
      !take the maximum value found on any processor
      call mpi_allreduce(mdGlobalSurfxx(k,i),mdGlobalSurfxx(k,i),1,sumb_real,MPI_MAX,sumb_comm_world,ierr)
   enddo
   call mpi_barrier(sumb_comm_world, ierr)
   
enddo
if (myID==0) then
   print *,'global reduced surface coordinates: x',mdGlobalSurfxx(1,:),'y',mdGlobalSurfxx(2,:),'z',mdGlobalSurfxx(3,:)
endif



end subroutine mdCreateGlobalReducedSurfaceList
