!
!       File:          adjustInflowAngle.f90                           
!       Author:        C.A.(Sandy) Mader                               
!       Starting date: 07-13-2011                                      
!       Last modified: 07-13-2011                                      
!
subroutine adjustInflowAngle(alpha, beta, liftIndex)

  use constants
  use inputPhysics

  implicit none

  !Subroutine Vars
  real(kind=realType), intent(in) :: alpha, beta
  integer(kind=intType), intent(in) :: liftIndex

  !Local Vars
  real(kind=realType), dimension(3) :: refDirection

  ! Velocity direction given by the rotation of a unit vector
  ! initially aligned along the positive x-direction (1,0,0)
  ! 1) rotate alpha radians cw about y or z-axis
  ! 2) rotate beta radians ccw about z or y-axis

  refDirection(:) = zero
  refDirection(1) = one
  call getDirVector(refDirection, alpha, beta, velDirFreestream,&
       liftIndex)

  ! Drag direction given by the rotation of a unit vector
  ! initially aligned along the positive x-direction (1,0,0)
  ! 1) rotate alpha radians cw about y or z-axis
  ! 2) rotate beta radians ccw about z or y-axis

  refDirection(:) = zero
  refDirection(1) = one
  call getDirVector(refDirection, alpha, beta, dragDirection, &
       liftIndex)

  ! Lift direction given by the rotation of a unit vector
  ! initially aligned along the positive z-direction (0,0,1)
  ! 1) rotate alpha radians cw about y or z-axis
  ! 2) rotate beta radians ccw about z or y-axis

  refDirection(:) = zero
  refDirection(liftIndex) = one

  call getDirVector(refDirection, alpha, beta,liftDirection, liftIndex)

end subroutine adjustInflowAngle
