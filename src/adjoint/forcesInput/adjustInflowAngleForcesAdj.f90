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
subroutine adjustInflowAngleForcesAdj(alphaAdj,betaAdj,velDirFreestreamAdj,&
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
      ! 1) rotate alpha radians cw about z-axis
      ! 2) rotate beta radians ccw about y-axis

      refDirection(:) = zero
      refDirection(1) = one
      call getDirVectorForces(refDirection, alphaAdj, betaAdj,&
           velDirFreestreamAdj,liftIndex)

      refDirection(:) = zero
      refDirection(1) = one
      call getDirVectorForces(refDirection, alphaAdj, betaAdj, &
           dragDirectionAdj(1),liftIndex)

      refDirection(:) = zero
      refDirection(liftIndex) = one
      
      call getDirVectorForces(refDirection, alphaAdj, betaAdj, &
           liftDirectionAdj,liftIndex)


    end subroutine adjustInflowAngleForcesAdj
