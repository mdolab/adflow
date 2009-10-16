!
!      ******************************************************************
!      *                                                                *
!      * File:          secondDerivativeRigidRotAngle.f90               *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 10-14-2009                                      *
!      * Last modified: 10-14-2009                                      *
!      *                                                                *
!      ******************************************************************
!
       function secondDerivativeRigidRotAngle(degreePolRot,   &
                                           coefPolRot,     &
                                           degreeFourRot,  &
                                           omegaFourRot,   &
                                           cosCoefFourRot, &
                                           sinCoefFourRot, t)
!
!      ******************************************************************
!      *                                                                *
!      * 2ndderivativeRigidRotAngle computes the 2nd time derivative of *
!      * the rigid body rotation angle at the given time for the given  *
!      * arguments. The angle is described by a combination of a        *
!      * polynomial and fourier series.                                 *
!      *                                                                *
!      ******************************************************************
!
       use flowVarRefState
       use inputPhysics
       implicit none
!
!      Function type
!
       real(kind=realType) :: secondDerivativeRigidRotAngle
!
!      Function arguments.
!
       integer(kind=intType), intent(in) :: degreePolRot
       integer(kind=intType), intent(in) :: degreeFourRot

       real(kind=realType), intent(in) :: omegaFourRot, t

       real(kind=realType), dimension(0:*), intent(in) :: coefPolRot
       real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourRot
       real(kind=realType), dimension(*),   intent(in) :: sinCoefFourRot
!
!      Local variables.
!
       integer(kind=intType) :: nn

       real(kind=realType) :: dPhi, val
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if this is a steady computation.

       if(equationMode == steady) then
         secondDerivativeRigidRotAngle = zero
         return
       endif

       ! Compute the polynomial contribution.

       dPhi = zero
       do nn=2,degreePolRot
         dPhi = dPhi + (nn-1)*nn*coefPolRot(nn)*(t**(nn-2))
       enddo

       ! Compute the fourier contribution.

       do nn=1,degreeFourRot
         val = nn*omegaFourRot
         dPhi = dPhi - val**2*cosCoefFourRot(nn)*sin(val*t)
         dPhi = dPhi - val**2*sinCoefFourRot(nn)*cos(val*t)
       enddo

       ! Set derivativeRigidRotAngle to dPhi. Multiply by timeRef
       ! to obtain the correct non-dimensional value.

       secondDerivativeRigidRotAngle = timeRef**2*dPhi

     end function secondDerivativeRigidRotAngle
