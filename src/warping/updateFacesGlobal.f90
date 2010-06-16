
!
! ***********************************
! *  File: updateFacesGlobal.f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 07-13-2009
! *  Modified: 07-13-2009
! ***********************************

subroutine updateFacesGlobal(ncoords,xyz_new,reinitialize)

  use blockPointers
  use mdData       !mdNSurfNodesCompact
  use mdDataLocal   !mdSurfGlobalIndLocal
  use inputTimeSpectral !nTimeIntervalsSpectral
  use section       !nSection, section%
  use monitor       !timeUnsteadyRestart
  implicit none

  !Subroutine Arguments
  logical  :: reinitialize
  integer(kind=intType)::ncoords
  real(kind=realType), dimension(3,ncoords)::xyz_new
  !integer(kind=intType),dimension(5,ncoords)::indices_new

  !Local variables

  integer(kind=intType)::level=1,sps,i,j,k,nn
  real(kind=realType), dimension(nSections) :: dt, t
  real(kind=realType) :: displX, displY, displZ
  real(kind=realType) :: tNew, tOld
  real(kind=realType) :: xp, yp, zp

  real(kind=realType), dimension(3)   :: rotationPoint
  real(kind=realType), dimension(3,3) :: rotationMatrix

 !reset block values
  if (reinitialize) then
     do sps = 1,nTimeIntervalsSpectral
        do nn = 1,ndom
           call setPointers(nn,1,sps)
           x = xInit
        enddo
     end do
     call xhalo(level)
  end if
  
!
!      ******************************************************************
!      *                                                                *
!      * Step 1. Perform a rigid body motion of the coordinates of the  *
!      *         1st time instance to the other instances.              *
!      *                                                                *
!      ******************************************************************
!
  
  ! Determine the delta t for every section. Remember it is possible
  ! that every section has a different periodic time.
  
  do nn=1,nSections
     dt(nn) = sections(nn)%timePeriod &
                / real(nTimeIntervalsSpectral,realType)
  enddo

  timeUnsteady = zero
  
  do sps = 1,nTimeIntervalsSpectral
     ! Compute the corresponding times for this spectral solution
     ! and call updateCoorFineMesh to determine the coordinates.
     
     do nn=1,nSections
        t(nn) = (sps-1)*dt(nn)
     enddo
     
     ! Compute the displacements due to the rigid motion of the mesh.
     
     displX = zero
     displY = zero
     displZ = zero
     
     ! Determine the time values of the old and new time level.
     ! It is assumed that the rigid body rotation of the mesh is only
     ! used when only 1 section is present.
     
     tNew = timeUnsteady + timeUnsteadyRestart
     tOld = tNew - t(1)
     
     ! Compute the rotation matrix of the rigid body rotation as
     ! well as the rotation point; the latter may vary in time due
     ! to rigid body translation.
     
     call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)
     
     !print *,'ncoords',ncoords,shape(xyz_new),xyz_new
     ! update working block coordinates with new coords
     !loop over the incoming coordinates
     !print *,'ncoords',ncoords,mdNSurfNodesCompact
     !stop
     do i =1 ,ncoords
        !Loop over the surface nodes on this processor
        do j = 1,size(mdSurfGlobalIndLocal(5,:))
           !update if they match the current global surface node
           !print *,'surfaceindex',int(mdSurfGlobalIndLocal(5,j)),i-1,mdSurfGlobalIndLocal(4,j)
           if (int(mdSurfGlobalIndLocal(5,j))==i-1)then
              call setPointers(mdSurfGlobalIndLocal(4,j),level,sps)
              !if( mdSurfGlobalIndLocal(4,j)==11)then
              !   print *,'xyzold',x(mdSurfGlobalIndLocal(1,j),mdSurfGlobalIndLocal(2,j),mdSurfGlobalIndLocal(3,j),k),mdSurfGlobalIndLocal(1,j),mdSurfGlobalIndLocal(2,j),mdSurfGlobalIndLocal(3,j)
              !end if
              
              ! Determine the vector relative to the rotation point.
              
              xp = xyz_new(1,i) - rotationPoint(1)
              yp = xyz_new(2,i) - rotationPoint(2)
              zp = xyz_new(3,i) - rotationPoint(3)

              ! Apply the transformation matrix to the vector (xp,yp,zp)
              ! and set the new coordinates.
              
              x(mdSurfGlobalIndLocal(1,j),mdSurfGlobalIndLocal(2,j),mdSurfGlobalIndLocal(3,j),1)&
                   = rotationMatrix(1,1)*xp &
                   + rotationMatrix(1,2)*yp &
                   + rotationMatrix(1,3)*zp + rotationPoint(1)+ displX
              x(mdSurfGlobalIndLocal(1,j),mdSurfGlobalIndLocal(2,j),mdSurfGlobalIndLocal(3,j),2)&
                   = rotationMatrix(2,1)*xp &
                   + rotationMatrix(2,2)*yp &
                   + rotationMatrix(2,3)*zp + rotationPoint(2)+ displY
              x(mdSurfGlobalIndLocal(1,j),mdSurfGlobalIndLocal(2,j),mdSurfGlobalIndLocal(3,j),3)&
                   = rotationMatrix(3,1)*xp &
                   + rotationMatrix(3,2)*yp &
                   + rotationMatrix(3,3)*zp + rotationPoint(3)+ displZ

              mdGlobalSurfxx(1,i,sps) = rotationMatrix(1,1)*xp &
                   + rotationMatrix(1,2)*yp &
                   + rotationMatrix(1,3)*zp + rotationPoint(1)+ displX
              mdGlobalSurfxx(2,i,sps) = rotationMatrix(2,1)*xp &
                   + rotationMatrix(2,2)*yp &
                   + rotationMatrix(2,3)*zp + rotationPoint(2)+ displY
              mdGlobalSurfxx(3,i,sps) = rotationMatrix(3,1)*xp &
                   + rotationMatrix(3,2)*yp &
                   + rotationMatrix(3,3)*zp + rotationPoint(3)+ displZ
              !if(isnan(sum(mdGlobalSurfxx(:,i,sps))))then
              !     print *,'globalsurfacenan',mdGlobalSurfxx(:,i,sps),i
              !endif
              !print *,'xsurface',x(mdSurfGlobalIndLocal(1,j),mdSurfGlobalIndLocal(2,j),mdSurfGlobalIndLocal(3,j),k), xyz_new(k,i),sps
              !x(mdSurfGlobalIndLocal(1,j),mdSurfGlobalIndLocal(2,j),mdSurfGlobalIndLocal(3,j),k) = xyz_new(k,i)
              
              !if( mdSurfGlobalIndLocal(4,j)==11)then
              !   print *,'xyznew',xyz_new(k,i)
              !endif
              
           endif
        enddo
     end do
  end do
end subroutine updateFacesGlobal
