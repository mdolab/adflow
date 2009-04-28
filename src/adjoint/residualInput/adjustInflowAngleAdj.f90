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
   real(kind=realType), dimension(3) :: velDirFreestreamAdj
   real(kind=realType), dimension(3) :: liftDirectionAdj
   real(kind=realType), dimension(3) :: dragDirectionAdj

   real(kind=realType) :: alphaAdj, betaAdj
   integer(kind=intType)::liftIndex

!Local Vars
   
   real(kind=realType) ::temp1, temp2, temp3
   real(kind=realType), dimension(3) :: refDirection
!Begin Execution

      ! Velocity direction given by the rotation of a unit vector
      ! initially aligned along the positive x-direction (1,0,0)
      ! 1) rotate alpha radians cw about y or z-axis
      ! 2) rotate beta radians ccw about z or y-axis

      !temp1 = velDirFreestreamAdj(1)
      !temp2 = velDirFreestreamAdj(2)
      !temp3 = velDirFreestreamAdj(3)

!      call getDirVector(one, zero, zero, alphaAdj, betaAdj, &
!                        velDirFreestreamAdj(1), &
!                        velDirFreestreamAdj(2), &
!                        velDirFreestreamAdj(3))
      refDirection(:) = zero
      refDirection(1) = one
      call getDirVector(refDirection, alphaAdj, betaAdj, velDirFreestreamAdj,&
           liftIndex)
      
!!$      velDirFreestreamAdj(1) = temp1
!!$      velDirFreestreamAdj(2) = temp2
!!$      velDirFreestreamAdj(3) = temp3

      
      ! Drag direction given by the rotation of a unit vector
      ! initially aligned along the positive x-direction (1,0,0)
      ! 1) rotate alpha radians cw about y or z-axis
      ! 2) rotate beta radians ccw about z or y-axis

!      call getDirVector(one, zero, zero, alphaAdj, betaAdj,   &
!                        dragDirectionAdj(1), dragDirectionAdj(2), &
!                        dragDirectionAdj(3))

      !temp1 = dragDirectionAdj(1)
      !temp2 = dragDirectionAdj(2)
      !temp3 = dragDirectionAdj(3)
      refDirection(:) = zero
      refDirection(1) = one
      call getDirVector(refDirection, alphaAdj, betaAdj, dragDirectionAdj(1), &
           liftIndex)

!!$      dragDirectionAdj(1)= temp1
!!$      dragDirectionAdj(2)= temp2
!!$      dragDirectionAdj(3)= temp3

      ! Lift direction given by the rotation of a unit vector
      ! initially aligned along the positive z-direction (0,0,1)
      ! 1) rotate alpha radians cw about y or z-axis
      ! 2) rotate beta radians ccw about z or y-axis

!      call getDirVector(zero,one, zero, alphaAdj, betaAdj,   &
!                        liftDirectionAdj(1), liftDirectionAdj(2), &
!                        liftDirectionAdj(3))

      !temp1 = liftDirectionAdj(1)
      !temp2 = liftDirectionAdj(2)
      !temp3 = liftDirectionAdj(3)

     
         refDirection(:) = zero
         refDirection(liftIndex) = one
         
         call getDirVector(refDirection, alphaAdj, betaAdj,liftDirectionAdj, &
              liftIndex)


    end subroutine adjustInflowAngleAdj
