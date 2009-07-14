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
use mdDataLocal
use mdData
implicit none

! Local Arguments

integer(kind=intType)::ncoords,nn,ierr,level = 1,i
integer(kind=intType)::sps = 1,famid = 0,startind,endInd
!real(kind=realType), dimension(3,1)::xyzface
!integer(kind=intType),dimension(4,1)::indices
real(kind=realType), dimension(:,:),allocatable::xyzface
integer(kind=intType),dimension(:,:),allocatable::indices
real(kind=realType)::dx,dx1
!real(kind=realType), dimension(3,ncoords)::xyzface
!integer(kind=intType),dimension(4,ncoords)::indices

!derivative variables
real(kind=realType), dimension(:,:),allocatable::surfacePoints,surfacePointsRef
integer(kind=intType)::j,k,nnn,index
!
! Begin execution
!



print *,'getting local surface'
call mdCreateSurfCoorListLocal(sps,famid,startInd,endInd)
print *,'getting local indices'
call mdCreateSurfIndListLocal(famID,startInd,endInd)

print *,'count total number of nodes'
call mdCreateNsurfNodes

!stop
print *,'testing global index creation'
call synchronizeIndices
print *,'synchronizing done'
!do i = 1,mdNSurfNodesLocal(1)
!   print *,'index',mdSurfxxLocal(1,i),mdSurfIndLocal(5,i)
!end do
! Assertion testing, memory allocation and global number indexing.
print *,'setting up global volume numbering'
call preprocessingADjoint(level)
!initialize Petsc
call initializePETSc
print *,'petscInitialized'
call warpingInitializePETSc
print *,'warpingPETSc initialized'
call createWarpingPETScMat
print *,'Warping PETSC Mat Created'
!print *,'mdNSurfNodesLocal',mdNSurfNodesLocal(1),'whole',mdNSurfNodesLocal
allocate(xyzface(3,mdNSurfNodesLocal(1)),indices(5,mdNSurfNodesLocal(1)))


!allocate temporary storage for the Surface points
allocate(surfacePoints(3,mdNSurfNodesCompact),surfacePointsRef(3,mdNSurfNodesCompact), stat=ierr)
if(ierr /= 0)                         &
     call terminate("verifyWarpDerivFD", &
     "Memory allocation failure for surfacePoints")

call mdCreateGlobalReducedSurfaceList

surfacePoints=mdGlobalSurfxx


!do nn = 1,ndom
!   call setPointers(nn,1,1)

dx = 0!0.5!0.1!1.e-15
dx1 = 0!1.e-4!0!1.e-10
ncoords = mdNSurfNodesLocal(1)

xyzface(:,:) = mdSurfxxLocal(:,:)+dx1
!xyzface(1,:) = mdSurfxxLocal(1,:)!+dx1
indices(:,:) = mdSurfIndLocal(:,:)
   xyzface(1,1) = xyzface(1,1)+dx
   xyzface(2,1) = xyzface(2,1)+dx
   xyzface(3,1) = xyzface(3,1)+dx
!!$if (myid==0)then
!!$   xyzface(1,1) = x(1,1,1,1)+ dx
!!$   xyzface(2,1) = x(1,1,1,2)+dx
!!$   xyzface(3,1) = x(1,1,1,3)+dx
!!$else
!!$   xyzface(1,1) = x(1,1,1,1)!+ dx
!!$   xyzface(2,1) = x(1,1,1,2)!+dx
!!$   xyzface(3,1) = x(1,1,1,3)!+dx
!!$endif
!4th index is block
!indices(:,1) = (/1,1,1,1/)
newgridfile = 'testwarpbefore.cgns'
writeGrid = .true.
writeVolume  = .true.
if(writeGrid .or. writeVolume .or. writeSurface) &
     call writeSol
!new warping method
do nnn=1,nDom
   call setpointersadj(nnn,level,sps)
   !print *,'xplus2',xplus(1,1,1,1)
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

!Store the current Face
surfacePointsRef = Surfacepoints
index = 20
dx =-0.0
surfacePoints(1,index) = surfacePointsRef(1,index)+dx
surfacePoints(2,index) = surfacePointsRef(2,index)+dx
surfacePoints(3,index) = surfacePointsRef(3,index)+dx
call updateFacesGlobal(mdNSurfNodesCompact,surfacePoints)
!now warp the domain
call warpMesh

!print *,'warping parameters',ncoords!,xyzface,indices
!call integratedWarp(ncoords,xyzface,indices)

newgridfile = 'testwarpafter.cgns'
writeGrid = .true.
writeVolume  = .true.
if(writeGrid .or. writeVolume .or. writeSurface) &
     call writeSol
!stop
call verifyWarpDerivFD
!stop
print *,'warping parameters3',ncoords!,xyzface,indices
call integratedWarpDerivFD(ncoords,xyzface,indices)
print *,'FD Derivatives done',myID
stop
call mpi_barrier(sumb_comm_world, ierr)

print *,'warping parameters2',ncoords!,xyzface,indices
call integratedWarpDeriv(ncoords,xyzface,indices)
print *,'AD Derivatives done',myID
!stop
call integratedWarpDerivParallel
call mpi_barrier(sumb_comm_world, ierr)
print *,'parallel derivatives done'
stop
end subroutine testwarping
