!
!      ******************************************************************
!      *                                                                *
!      * File:          rigidRotAngle.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-01-2004                                      *
!      * Last modified: 03-23-2004                                      *
!      *                                                                *
!      ******************************************************************
!
       function rigidRotAngle(degreePolRot,   coefPolRot,       &
                              degreeFourRot,  omegaFourRot,     &
                              cosCoefFourRot, sinCoefFourRot, t)
!
!      ******************************************************************
!      *                                                                *
!      * rigidRotAngle computes the rigid body rotation angle at the    *
!      * given time for the given arguments. The angle is described by  *
!      * a combination of a polynomial and fourier series.              *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use inputPhysics
       implicit none
!
!      Function type
!
       real(kind=realType) :: rigidRotAngle
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

       real(kind=realType) :: phi, val
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if this is a steady computation.

       if(equationMode == steady) then
         rigidRotAngle = zero
         return
       endif

       ! Compute the polynomial contribution. If no polynomial was
       ! specified, the value of index 0 is set to zero automatically.

       phi = coefPolRot(0)
       do nn=1,degreePolRot
         phi = phi + coefPolRot(nn)*(t**nn)
       enddo

       ! Compute the fourier contribution. Again the cosine coefficient
       ! of index 0 is defaulted to zero if not specified.

       phi = phi + cosCoefFourRot(0)
       do nn=1,degreeFourRot
         val = nn*omegaFourRot*t
         phi = phi + cosCoefFourRot(nn)*cos(val) &
                   + sinCoefFourRot(nn)*sin(val)
       enddo

       ! Set rigidRotAngle to phi.

       rigidRotAngle = phi

       end function rigidRotAngle
