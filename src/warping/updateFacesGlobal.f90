!
! ***********************************
! *  File: updateFacesGlobal.f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 07-13-2009
! *  Modified: 07-13-2009
! ***********************************

subroutine updateFacesGlobal(ncoords,xyz_new)

  use blockPointers
!  use mdData       !mdNSurfNodesCompact
  use mdDataLocal   !mdSurfGlobalIndLocal
  implicit none

  !Subroutine Arguments
  integer(kind=intType)::ncoords
  real(kind=realType), dimension(3,ncoords)::xyz_new
  !integer(kind=intType),dimension(5,ncoords)::indices_new

  !Local variables

  integer(kind=intType)::level=1,sps=1,i,j,k

    
  ! update working block coordinates with new coords
  !loop over the incomping coordinates
  do i =1 ,ncoords
     !Loop over the surface nodes on this processor
     do j = 1,size(mdSurfGlobalIndLocal(5,:))
        !update if they match the current global surface node
        !print *,'surfaceindex',int(mdSurfGlobalIndLocal(5,j)),i-1,mdSurfGlobalIndLocal(4,j)
        if (int(mdSurfGlobalIndLocal(5,j))==i-1)then
           do k = 1,3
              call setPointers(mdSurfGlobalIndLocal(4,j),level,sps)
              !print *,'xyzold',x(mdSurfGlobalIndLocal(1,j),mdSurfGlobalIndLocal(2,j),mdSurfGlobalIndLocal(3,j),k)
              x(mdSurfGlobalIndLocal(1,j),mdSurfGlobalIndLocal(2,j),mdSurfGlobalIndLocal(3,j),k) = xyz_new(k,i)
              !print *,'xyznew',xyz_new(k,i)
           end do
           
        endif
     enddo
  end do
  
end subroutine updateFacesGlobal
