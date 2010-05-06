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
  use mdData        !mdNSurfNodesCompact
  implicit none


!Subroutine Variables



!Local Variables
integer(kind=intType)::level=1,sps=1,nn,i,j,k,ierr,counter,compresscounter
logical::indexFound,test
!testing vars
real(kind=realtype),dimension(:,:),allocatable::testindices
!begin Execution

!Allocate storage for counter tracking the number of global surface 
!nodes on each processor
  if(.not. allocated(mdNGlobalSurfNodesLocal) ) then
     allocate(mdNGlobalSurfNodeslocal(nProc), stat=ierr)
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
               !print *,'counter',counter
            endif
         end DO
      end DO
   end DO
end DO

!Collect this info on all processors
call mpi_allgather(counter, 1, sumb_integer, &
     mdNGlobalSurfNodesLocal,1, sumb_integer, &
     SUmb_comm_world, ierr)
!print *,'N Global Surface nodes',myID,counter,mdNGlobalSurfNodesLocal
!Allocate storage for counter tracking the number of global surface 
!nodes on each processor
  if(.not. allocated(mdSurfGlobalIndLocal) ) then
     allocate(mdSurfGlobalIndLocal(5,counter), stat=ierr)
     if(ierr /= 0)                          &
          call terminate("storeGlobalSurfaceIndices", &
          "Memory allocation failure for mdSurfGlobalIndLocal")
  endif
  allocate(testindices(5,counter), stat=ierr)

!put synchronized indices in mdSurfGlobalIndLocal
counter = 0
do nn=1,nDom
   call setPointers(nn,level,sps)
   DO I=1,il!IMAX
      DO J=1,jl!JMAX
         DO K=1,kl!KMAX
            if (X(I,J,K,1) /= -5)then
               !print *,'indices',nn,ndom,i,il,j,jl,k,kl,counter
               mdSurfGlobalIndLocal(1,counter+1) = i
               mdSurfGlobalIndLocal(2,counter+1) = j
               mdSurfGlobalIndLocal(3,counter+1) = k
               mdSurfGlobalIndLocal(4,counter+1) = nn
               mdSurfGlobalIndLocal(5,counter+1) = X(i,j,k,1)-1
!!$               if( nn == 11)then
!!$                  print *,'block 11',X(i,j,k,1),i,j,k,nn
!!$               end if
               !print *,'global index', mdSurfGlobalIndLocal(5,counter+1),counter,counter+1
               counter = counter+1
            endif
         end DO
      end DO
   end DO
end DO

!!$!check current index numbering
!!$do i = 1,counter!+1
!!$   do j = 1,counter!+1
!!$      print *,'i,j',i,j
!!$      if (int(mdSurfGlobalIndLocal(5,j))==i-1)then
!!$         print *,'i',i-1,j
!!$      endif
!!$   end do
!!$end do
!!$!stop

!compress to be continuous numbering
compresscounter = 0
do i = 1,mdNSurfNodes(nProc,1)!maximum number of possible surface nodes??counter+1
   !if (myid==0) print *,'global surface node...',i,mdNSurfNodes(nProc,1)
   !call sleep(0.2)
   indexfound = .False.
   do j = 1,counter!+1
      if (int(mdSurfGlobalIndLocal(5,j))==i-1)then
         mdSurfGlobalIndLocal(5,j)=compresscounter
!!$         if( i>25 .and. i<51)then
!!$            print *,'compressing',i,compresscounter
!!$         endif
         !testindices(5,j)=compresscounter
         !print *,'indices',int(mdSurfGlobalIndLocal(5,j)),compresscounter,myid,j
         !call sleep(0.2)
         indexfound = .True.
      endif
   enddo
   !print *,'indexfound before',myID, indexFound
   call mpi_barrier(sumb_comm_world, ierr)
   !if the current I index is found on any processor, 
   !identify index as having been found.
   !call mpi_allreduce(indexFound,indexFound,1,MPI_LOGICAL,MPI_LOR,sumb_comm_world,ierr)
   call mpi_allreduce(indexFound,test,1,MPI_LOGICAL,MPI_LOR,sumb_comm_world,ierr)
   !print *,'indexfound after',myID, indexFound
   call mpi_barrier(sumb_comm_world, ierr)
   !if the index was found (and hence replaced), increment the index counter
   !if (indexfound)then
   if(test) then
      compresscounter = compresscounter+1
   endif
enddo
!print *,'compresscounter',compresscounter
!reprint index
mdNSurfNodesCompact= compresscounter!-1
!print *,'storing indices',mdnsurfnodescompact
!!$do i = 1,counter!+1
!!$   do j = 1,counter!+1
!!$      if (int(mdSurfGlobalIndLocal(5,j))==i-1)then
!!$         !if (int(testindices(5,j))==i-1)then
!!$         !print *,'i2',myID,i-1,int(testindices(5,j)),j,int(mdSurfGlobalIndLocal(5,j))
!!$         print *,'i2',myID,i-1,j,int(mdSurfGlobalIndLocal(5,j))
!!$      endif
!!$      !call sleep(0.2)
!!$   end do
!!$end do
call mpi_barrier(sumb_comm_world, ierr)     

end subroutine storeGlobalSurfaceIndices
