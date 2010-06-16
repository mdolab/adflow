!
! ***********************************
! *  File: synchronizeIndices.f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 05-20-2009
! *  Modified: 05-20-2009
! ***********************************

subroutine synchronizeIndices

  use BCTypes
  use blockPointers
  use mddatalocal
  use mddata
  use cgnsGrid
  use communication !myid

  implicit none

  !Local Variables
  integer(kind=intType)::level=1,sps=1,nn,mm,i,j,k,ierr
  integer(kind=intType) :: famBCType,ifam,mstart,mstop
  logical :: checkFamily
  integer(kind=intType) :: currBC,famID

  !temporarily reset the coordinates to accomodate the surface synch -- Reset at end
  do nn=1,nDom
     call setPointers(nn,level,sps)
     x(:,:,:,1) = -5
  enddo

  if (mdnfamilieslocal == 1) then
     checkFamily = .False.
  else
     checkFamily = .True.
  end if
!   print *,'mdnfamilieslocal',mdnfamilieslocal
!   print *,'mdnsurfnodes:',mdNsurfNodesLocal
!   print *,'checkFamily:',checkFamily
!   print *,'mdNSurfNodesProc(myID)',mdNSurfNodesProc(myID)


  if (checkFamily) then
     do nn = 1,nDom
        call setPointers(nn,1,sps)
        
        do mm=1,mdNSurfNodesLocal(mdnfamilieslocal)
           !Check to see that coordinate is in this block. if so, update
           if(mdSurfIndLocal(4,mm)==nn) then
              famID = mdSurfIndLocal(6,mm) ! Check this famid
              currBC = cgnsFamilies(famID)%BCType
              if(currBC == EulerWall.or.currBC== NSWallAdiabatic .or.currBC==NSWallIsothermal) then
                 x(mdSurfIndLocal(1,mm),mdSurfIndLocal(2,mm),mdSurfIndLocal(3,mm),1) = mdSurfIndLocal(5,mm)
              endif
           end if
        end do
     end do

  else ! No family BCs -> All surface nodes used for warping
     do nn = 1,nDom
        call setPointers(nn,1,sps)
        do mm=1,mdSumNsurfNodesLocal
           !Check to see that coordinate is in this block. if so, update
           if(mdSurfIndLocal(4,mm)==nn) then
              x(mdSurfIndLocal(1,mm),mdSurfIndLocal(2,mm),mdSurfIndLocal(3,mm),1) = mdSurfIndLocal(5,mm)
           endif
        end do
     end do
  end if


  call mpi_barrier(sumb_comm_world, ierr)

!run syncronize faces
!what about duplicate nodes at split surface boundaries? Use sychronization to set common nodes to the lower of the two index values.
call synchronizeSurfaceIndices(level,sps)
call mdCreateNsurfNodes ! Needed for storeglobalsurfaceindices
call storeGlobalSurfaceIndices

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

end subroutine synchronizeIndices
