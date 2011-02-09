!
!      ******************************************************************
!      *                                                                *
!      * File:          derivativeRigidRotAngle.f90                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-01-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       function derivativeRigidRotAngle(degreePolRot,   &
                                        coefPolRot,     &
                                        degreeFourRot,  &
                                        omegaFourRot,   &
                                        cosCoefFourRot, &
                                        sinCoefFourRot, t)
!
!      ******************************************************************
!      *                                                                *
!      * derivativeRigidRotAngle computes the time derivative of the    *
!      * rigid body rotation angle at the given time for the given      *
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
       real(kind=realType) :: derivativeRigidRotAngle
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
         derivativeRigidRotAngle = zero
         return
       endif

       ! Compute the polynomial contribution.

       dPhi = zero
       do nn=1,degreePolRot
         dPhi = dPhi + nn*coefPolRot(nn)*(t**(nn-1))
       enddo

       ! Compute the fourier contribution.

       do nn=1,degreeFourRot
         val = nn*omegaFourRot
         dPhi = dPhi - val*cosCoefFourRot(nn)*sin(val*t)
         dPhi = dPhi + val*sinCoefFourRot(nn)*cos(val*t)
       enddo

       ! Set derivativeRigidRotAngle to dPhi. Multiply by timeRef
       ! to obtain the correct non-dimensional value.

       derivativeRigidRotAngle = timeRef*dPhi

       end function derivativeRigidRotAngle
