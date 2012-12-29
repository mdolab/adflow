!
!      ******************************************************************
!      *                                                                *
!      * File:          adjustInflowAngle.f90                           *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 07-13-2011                                      *
!      * Last modified: 07-13-2011                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine adjustInflowAngle(alpha,beta,liftIndex)

  use constants
  use allInputParam !inputPhysics !veldirFreestream, liftDiretion, dragDirection
  implicit none

  !Subroutine Vars
  real(kind=realType) :: alpha, beta

  integer(kind=intType)::liftIndex

  !Local Vars

  real(kind=realType), dimension(3) :: refDirection
  !Begin Execution

  ! Velocity direction given by the rotation of a unit vector
  ! initially aligned along the positive x-direction (1,0,0)
  ! 1) rotate alpha radians cw about y or z-axis
  ! 2) rotate beta radians ccw about z or y-axis

  refDirection(:) = zero
  refDirection(1) = one
  call getDirVectorTS(refDirection, alpha, beta, velDirFreestream,&
       liftIndex)


  ! Drag direction given by the rotation of a unit vector
  ! initially aligned along the positive x-direction (1,0,0)
  ! 1) rotate alpha radians cw about y or z-axis
  ! 2) rotate beta radians ccw about z or y-axis


  refDirection(:) = zero
  refDirection(1) = one
  call getDirVectorTS(refDirection, alpha, beta, dragDirection, &
       liftIndex)


  ! Lift direction given by the rotation of a unit vector
  ! initially aligned along the positive z-direction (0,0,1)
  ! 1) rotate alpha radians cw about y or z-axis
  ! 2) rotate beta radians ccw about z or y-axis


  refDirection(:) = zero
  refDirection(liftIndex) = one

  call getDirVectorTS(refDirection, alpha, beta,liftDirection, &
       liftIndex)

  liftDirSpecified=.True.

end subroutine adjustInflowAngle
