!
!      ******************************************************************
!      *                                                                *
!      * File:          adjustInflowAngleAdj.f90                        *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 05-14-2008                                      *
!      * Last modified: 05-14-2008                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine adjustInflowAngleAdj(alphaAdj,betaAdj,velDirFreestreamAdj,&
     liftDirectionAdj,dragDirectionAdj,liftIndex)

  use constants
  implicit none

  !Subroutine Vars
  real(kind=realType) :: alphaAdj, betaAdj

  real(kind=realType), dimension(3) :: velDirFreestreamAdj
  real(kind=realType), dimension(3) :: liftDirectionAdj
  real(kind=realType), dimension(3) :: dragDirectionAdj

  integer(kind=intType)::liftIndex

  !Local Vars

  real(kind=realType) ::temp1, temp2, temp3
  real(kind=realType), dimension(3) :: refDirection
  !Begin Execution

  ! Velocity direction given by the rotation of a unit vector
  ! initially aligned along the positive x-direction (1,0,0)
  ! 1) rotate alpha radians cw about y or z-axis
  ! 2) rotate beta radians ccw about z or y-axis

  refDirection(:) = zero
  refDirection(1) = one
  call getDirVector(refDirection, alphaAdj, betaAdj, velDirFreestreamAdj,&
       liftIndex)


  ! Drag direction given by the rotation of a unit vector
  ! initially aligned along the positive x-direction (1,0,0)
  ! 1) rotate alpha radians cw about y or z-axis
  ! 2) rotate beta radians ccw about z or y-axis


  refDirection(:) = zero
  refDirection(1) = one
  call getDirVector(refDirection, alphaAdj, betaAdj, dragDirectionAdj(1), &
       liftIndex)


  ! Lift direction given by the rotation of a unit vector
  ! initially aligned along the positive z-direction (0,0,1)
  ! 1) rotate alpha radians cw about y or z-axis
  ! 2) rotate beta radians ccw about z or y-axis


  refDirection(:) = zero
  refDirection(liftIndex) = one

  call getDirVector(refDirection, alphaAdj, betaAdj,liftDirectionAdj, &
       liftIndex)


end subroutine adjustInflowAngleAdj
