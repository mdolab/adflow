!
! ***********************************
! *  File: setGrid.f90
! *  Author: Gaetan Kenway
! *  Started: 05-29-2010
! *  Modified: 05-29-2010
! ***********************************

subroutine setGrid(grid,ndof)

  ! The purpose of this routine is to set the grid dof as returned by
  ! the external warping. This function takes the "Base" grid at the
  ! first time instance and does rotation/translation operations to
  ! get the grid at subsequent time instances

  use precision 
  use blockPointers
  use communication
  use section
  use inputTimeSpectral
  use monitor 

  implicit none

  integer(kind=intType),intent(in) :: ndof
  real(kind=realType) ,intent(in) :: grid(ndof)

  ! Local Variables

  integer(kind=intType) :: nn,i,j,k,counter,idim,sps
  real(kind=realType) :: t(nSections),dt(nSections)
  real(kind=realType) :: displ(3)
  real(kind=realType) :: tOld,tNew

  real(kind=realType), dimension(3)   :: rotationPoint,r
  real(kind=realType), dimension(3,3) :: rotationMatrix

  do nn=1,nSections
     dt(nn) = sections(nn)%timePeriod &
          / real(nTimeIntervalsSpectral,realType)
  enddo
  
  timeUnsteady = zero
  
  ! This is very straight forward...loop over all domains and set all elements

  do sps = 1,nTimeIntervalsSpectral
     do nn=1,nSections
        t(nn) = (sps-1)*dt(nn)
     enddo
     
     ! Compute the displacements due to the rigid motion of the mesh.
     
     displ(:) = zero

     tNew = timeUnsteady + timeUnsteadyRestart
     tOld = tNew - t(1)

     call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)
     counter = 0
     do nn=1,nDom
        call setPointers(nn,1_intType,sps)
        do k=1,kl
           do j=1,jl
              do i=1,il
                 ! r is distance from grid point to rotationPoint
                 r = grid(3*counter+1:3*counter+3) - rotationPoint
                    
                 X(i,j,k,:) = rotationPoint + matmul(rotationMatrix,r) + displ
                 counter = counter + 1

              end do
           end do
        end do
     end do
     call xhalo(1_intType)
  end do

end subroutine setGrid

subroutine getGrid(grid,ndof)

  ! Opposite of setGrid. This is ONLY a debugging function. NOT used
  ! in regular usage. Really only useful for direct mesh manipulation
  ! on single block and a single processor. s

  use precision 
  use blockPointers
  use inputTimeSpectral
  implicit none
  integer(kind=intType),intent(in) :: ndof
  real(kind=realType) ,intent(out) :: grid(ndof)

  ! Local Variables
  integer(kind=intType) :: nn,i,j,k,l,counter,sps
  
  ! This is very straight forward...loop over all domains and copy out
  counter = 1
  do sps = 1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointers(nn,1_intType,sps)
        do k=1,kl
           do j=1,jl
              do i=1,il
                 do l=1,3
                    grid(counter) = X(i,j,k,l)
                    counter = counter + 1
                 end do
              end do
           end do
        end do
     end do
  end do

end subroutine getGrid

subroutine setGridSimple(grid,ndof)

  ! The purpose of this routine is to set the grid dof as returned by
  ! the external warping. This function takes the "Base" grid at the
  ! first time instance and does rotation/translation operations to
  ! get the grid at subsequent time instances

  use precision 
  use blockPointers
  use communication
  use section
  use inputTimeSpectral
  use monitor 

  implicit none

  integer(kind=intType),intent(in) :: ndof
  real(kind=realType) ,intent(in) :: grid(ndof)

  ! Local Variables

  integer(kind=intType) :: nn,i,j,k,l,counter,sps

  do sps = 1,nTimeIntervalsSpectral

     counter = 1
     do nn=1,nDom
        call setPointers(nn,1_intType,sps)
        do k=1,kl
           do j=1,jl
              do i=1,il
                 do l=1,3
                    X(i,j,k,l) = grid(counter)
                    counter = counter + 1
                 end do
              end do
           end do
        end do
     end do
     call xhalo(1_intType)
  end do

end subroutine setGridSimple
