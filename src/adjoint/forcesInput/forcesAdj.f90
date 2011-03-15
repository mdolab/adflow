!
!      ******************************************************************
!      *                                                                *
!      * File:          forcesAdj.f90                                   *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader           *
!      * Starting date: 08-17-2008                                      *
!      * Last modified: 08-17-2008                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine forcesAdj(pAdj,pts,normAdj,refPoint,force,moment,fact,&
     iBeg,iEnd,jBeg,jEnd,iNode,jNode)
  
  use flowVarRefState
  implicit none

  ! Subroutine Arguments
  real(kind=realType), intent(in) :: pAdj(3,2,2)
  real(kind=realType), intent(in) :: pts(3,3,3)
  real(kind=realType), intent(in) :: normAdj(3,2,2)
  real(kind=realType), intent(in) :: refPoint(3)
  real(kind=realType), intent(out):: force(3)
  real(kind=realType), intent(out):: moment(3)
  real(kind=realType), intent(in) :: fact
  integer(kind=intType),intent(in):: iBeg,iEnd,jBeg,jEnd,iNode,jNode

  ! Local Variables
  integer(kind=intTYpe) :: i,j
  real(kind=realType) :: pp,scaleDim,xc(3),r(3),addForce(3)
  real(kind=realType) :: tauXx, tauYy, tauZz
  real(kind=realType) :: tauXy, tauXz, tauYz

  scaleDim = pRef
  force(:) = 0.0
  moment(:) = 0.0
  ! Force is the contribution of each of 4 cells

  do j=1,2
     do i=1,2

        if (.not.(iNode+i-2 < iBeg .or. iNode+i-1 > iEnd .or. &
                  jNode+j-2 < jBeg .or. jNode+j-1 > jEnd)) then

           ! Calculate the Pressure
           pp = half*(pAdj(1,i,j) + pAdj(2,i,j)) - Pinf
           pp = fourth*fact*scaleDim*pp
           
           ! Incremental Force to Add
           addForce = pp*normAdj(:,i,j)

           ! Add increment to total Force for this node
           force(:) = force(:) + addForce(:)

           ! Cell Center, xc
           xc(:) = fourth*(pts(:,i,j)   + pts(:,i+1,j) + &
                           pts(:,i,j+1) + pts(:,i+1,j+1))

           ! Vector from center to refPoint
           r(:) = xc(:)-refPoint(:)

           ! Moment is F x r ( F cross r)
           moment(1) = moment(1) + r(2)*addForce(3) - r(3)*addForce(2)
           moment(2) = moment(2) + r(3)*addForce(1) - r(1)*addForce(3)
           moment(3) = moment(3) + r(1)*addForce(2) - r(2)*addForce(1)

  !          ! Viscous Force: Below is the code one would use to compute
!            ! the viscous forces. This code has NOT been verified, and
!            ! will not run. tau MUST be passed to this function, which
!            ! will have been computed by a routine similar to
!            ! viscousFlux.  The outer driving routines, will remain the
!            ! same (since this is just an extra w-dependance). 
!            if  (viscousSubface) then
!               tauXx = tau(i,j,1)
!               tauYy = tau(i,j,2)
!               tauZz = tau(i,j,3)
!               tauXy = tau(i,j,4)
!               tauXz = tau(i,j,5)
!               tauYz = tau(i,j,6)

!               ! Compute the viscous force on the face. A minus sign
!               ! is now present, due to the definition of this force.
              
!               addForce(1)= -fact*(tauXx*normAdj(1,i,j) + tauXy*normAdj(2,i,j) &
!                    +        tauXz*normAdj(3,i,j))
!               addForce(2) = -fact*(tauXy*normAdj(i,j,1) + tauYy*normAdj(2,i,j) &
!                    +        tauYz*normAdj(3,i,j))
!               addForce(3) = -fact*(tauXz*normAdj(1,i,j) + tauYz*normAdj(2,i,j) &
!                    +        tauZz*normAdj(3,i,j))

!               ! Increment the Force
!               force(:) = force(:) + addForce(:)
              
!               ! Increment the Moment
!               moment(1) = moment(1) + r(2)*addForce(3) - r(3)*addForce(2)
!               moment(2) = moment(2) + r(3)*addForce(1) - r(1)*addForce(3)
!               moment(3) = moment(3) + r(1)*addForce(2) - r(2)*addForce(1)
!            end if

         end if
      end do
   end do
 end subroutine forcesAdj
