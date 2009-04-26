!
! ***********************************
! *  File: updateFaces.f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 12-11-2008
! *  Modified: 12-11-2008
! ***********************************

subroutine updateFaces(ncoords,xyz_new,indices_new)

  use blockPointers
  implicit none

  !Subroutine Arguments
  integer(kind=intType)::ncoords
  real(kind=realType), dimension(3,ncoords)::xyz_new
  integer(kind=intType),dimension(4,ncoords)::indices_new

  !Local variables

  integer(kind=intType)::nn,mm,sps=1

    
  ! update working block coordinates with new coords

  !loop over domains
  do nn = 1,nDom
     call setPointers(nn,1,sps)
     !loop over new coordinates array
     do mm = 1,ncoords
        !Check to see that coordinate is in this block. if so, update
        if(indices_new(4,mm)==nn)then
           !print *,'indices',indices_new(1,mm),indices_new(2,mm),indices_new(3,mm),indices_new(4,mm)
           x(indices_new(1,mm),indices_new(2,mm),indices_new(3,mm),:) = xyz_new(:,mm)
        endif
     end do
  end do
  !stop
end subroutine updateFaces
