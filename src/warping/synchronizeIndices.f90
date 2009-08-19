!
! ***********************************
! *  File: synchronizeIndices.f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 05-20-2009
! *  Modified: 05-20-2009
! ***********************************

subroutine synchronizeIndices

use blockPointers
use mdData
use mdDataLocal
use communication !myid

implicit none

!Subroutine Arguments



!Local Variables
integer(kind=intType) :: startInd, endInd,famID
integer(kind=intType)::level=1,sps=1,nn,mm,i,j,k


!Begin Execution

famID  = 0
print *,'creating global list'
call mdCreateSurfIndList(famID,startInd,endInd)

print *,'creating local list'
call mdCreateSurfIndListLocal(famID,startInd,endInd)
!print *,'local indices',mdNSurfNodesLocal,'global',mdNSurfNodes!(nn)

!temporarily reset the coordinates to accomodate the surface synch
!set x to -5
do nn=1,nDom
   call setPointers(nn,level,sps)
   x(:,:,:,:) = -5
enddo

!set global index in corresponding x location.
!loop over domains
do nn = 1,nDom
   call setPointers(nn,1,sps)
   !loop over new coordinates array
   do mm = 1,mdNSurfNodesLocal(1)
      !Check to see that coordinate is in this block. if so, update
      if(mdSurfIndLocal(4,mm)==nn)then
         x(mdSurfIndLocal(1,mm),mdSurfIndLocal(2,mm),mdSurfIndLocal(3,mm),1) = mdSurfIndLocal(5,mm)
         !print *,'localblock',flowdoms(nn,level,sps)%cgnsblockid,myid,mdSurfIndLocal(1,mm),mdSurfIndLocal(2,mm),mdSurfIndLocal(3,mm),mdSurfIndLocal(5,mm)
      endif
   end do
end do
!stop
!run syncronize faces
!what about duplicate nodes at split surface boundaries? Use sychronization to set common nodes to the lower of the two index values.
call synchronizeSurfaceIndices(level,sps)
print *,'indices synchronized'

call storeGlobalSurfaceIndices
print *,'indices stored'
!run through each block, counting the number of global indices
!create a varible only large enough to hold local "globalindex" values. 
!store global index with corresponding local block and i,j,k info. 


!to differentiate loop through local list, storing values in column corresponding to global index


!reset the mesh coordinates to initial values to prepare for warp
do nn=1,nDom
   call setPointers(nn,level,sps)
   DO I=1,il!IMAX
      DO J=1,jl!JMAX
         DO K=1,kl!KMAX
            X(I,J,K,1) = Xinit(I,J,K,1)
            X(I,J,K,2) = Xinit(I,J,K,2)
            X(I,J,K,3) = Xinit(I,J,K,3)
         END DO
      END DO
   END DO
enddo
call xhalo(level)

!Now generate the list for interfacing with python
call mdCreateGlobalReducedSurfaceList
print *,'surface list generated'
!stop

end subroutine synchronizeIndices
