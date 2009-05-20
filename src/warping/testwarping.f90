!
! ***********************************
! *  File: testwarping..f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 12-12-2008
! *  Modified: 12-12-2008
! ***********************************

subroutine testwarping

use blockPointers
use communication
use monitor
use inputIO
implicit none

! Local Arguments

integer(kind=intType)::ncoords,nn
real(kind=realType), dimension(3,1)::xyzface
integer(kind=intType),dimension(4,1)::indices
real(kind=realType)::dx
!real(kind=realType), dimension(3,ncoords)::xyzface
!integer(kind=intType),dimension(4,ncoords)::indices

!
! Begin execution
!

nn = 1

call setPointers(nn,1,1)

dx = 0.1!1.e-15

ncoords = 1
if (myid==0)then
   xyzface(1,1) = x(1,1,1,1)+ dx
   xyzface(2,1) = x(1,1,1,2)+dx
   xyzface(3,1) = x(1,1,1,3)+dx
else
   xyzface(1,1) = x(1,1,1,1)!+ dx
   xyzface(2,1) = x(1,1,1,2)!+dx
   xyzface(3,1) = x(1,1,1,3)!+dx
endif
!4th index is block
indices(:,1) = (/1,1,1,1/)
newgridfile = 'testwarpbefore.cgns'
writeGrid = .true.
writeVolume  = .true.
if(writeGrid .or. writeVolume .or. writeSurface) &
     call writeSol
print *,'warping parameters',ncoords,xyzface,indices
call integratedWarp(ncoords,xyzface,indices)

newgridfile = 'testwarpafter.cgns'
writeGrid = .true.
writeVolume  = .true.
if(writeGrid .or. writeVolume .or. writeSurface) &
     call writeSol
print *,'warping parameters2',ncoords,xyzface,indices
call integratedWarpDeriv(ncoords,xyzface,indices)
!stop
print *,'warping parameters3',ncoords,xyzface,indices
call integratedWarpDerivFD(ncoords,xyzface,indices)
print *,'testing global index creation'
call synchronizeIndices
end subroutine testwarping
