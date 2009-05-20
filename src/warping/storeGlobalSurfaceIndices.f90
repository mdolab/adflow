!==========================================================
!
!       File: storeGlobalSurfaceIndices.f90
!       Author: C.A.(Sandy) Mader
!       Start Date: 05-20-2009
!       Modified Date: 05-20-2009
!
!==========================================================

!     ******************************************************************
!     *   storeGlobalSurfacaeIndices stores the global surface indices *
!     *   for the parallel derivative computation after sychronization.*
!     *                                                                *
!     *   05/20/09  C.A.Mader   Initial Implementation                 *
!     *                                                                *
!     *   C.A.(Sandy) Mader, UTIAS,Toronto,Canada                      *
!     ******************************************************************



subroutine storeGlobalSurfaceIndices
  !""" Store the synchronized surface indices
  !      """
  use blockPointers
  use communication, only: myID,sumb_comm_world,nProc
  use BCTypes
  use mdDataLocal
  implicit none


!Subroutine Variables



!Local Variables
integer(kind=intType)::level=1,sps=1,nn,i,j,k,ierr,counter


!begin Execution

!Allocate storage for counter tracking the number of global surface 
!nodes on each processor
  if(.not. allocated(mdNSurfNodesGlobal) ) then
     allocate(mdNSurfNodesGlobal(nProc), stat=ierr)
     if(ierr /= 0)                          &
          call terminate("storeGlobalSurfaceIndices", &
          "Memory allocation failure for mdNSurfNodesGlobal")
  endif


!Count the number of surface indices on this processor
counter = 0
do nn=1,nDom
   call setPointers(nn,level,sps)
   DO I=1,il!IMAX
      DO J=1,jl!JMAX
         DO K=1,kl!KMAX
            if (X(I,J,K,1) /= -5)then
               counter = counter+1
            endif
         end DO
      end DO
   end DO
end DO

!Collect this info on all processors
call mpi_allgather(counter, 1, sumb_integer, &
     mdNSurfNodesGlobal,1, sumb_integer, &
     SUmb_comm_world, ierr)

!Allocate storage for counter tracking the number of global surface 
!nodes on each processor
  if(.not. allocated(mdSurfGlobalIndLocal) ) then
     allocate(mdSurfGlobalIndLocal(5,counter), stat=ierr)
     if(ierr /= 0)                          &
          call terminate("storeGlobalSurfaceIndices", &
          "Memory allocation failure for mdSurfGlobalIndLocal")
  endif
counter = 0
do nn=1,nDom
   call setPointers(nn,level,sps)
   DO I=1,il!IMAX
      DO J=1,jl!JMAX
         DO K=1,kl!KMAX
            if (X(I,J,K,1) /= -5)then
               mdSurfGlobalIndLocal(1,counter) = i
               mdSurfGlobalIndLocal(2,counter) = j
               mdSurfGlobalIndLocal(3,counter) = k
               mdSurfGlobalIndLocal(4,counter) = nn
               mdSurfGlobalIndLocal(5,counter) = X(i,j,k,1)
               counter = counter+1
            endif
         end DO
      end DO
   end DO
end DO

end subroutine storeGlobalSurfaceIndices
