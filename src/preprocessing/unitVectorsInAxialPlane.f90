!
!      ******************************************************************
!      *                                                                *
!      * File:          unitVectorsInAxialPlane.f90                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-10-2005                                      *
!      * Last modified: 12-20-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine unitVectorsInAxialPlane(axis, vecR1, vecR2)
!
!      ******************************************************************
!      *                                                                *
!      * unitVectorsInAxialPlane computes from the given unit vector    *
!      * axis the two unit vectors which describe the plane normal to   *
!      * axis. There is of course an ambiguity in this choice, but this *
!      * is not a problem as long as the choice is consistent           *
!      * throughout the code.                                           *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(3), intent(in)  :: axis
       real(kind=realType), dimension(3), intent(out) :: vecR1, vecR2
!
!      Local variables.
!
       real(kind=realType) :: dot
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! The vectors which span the axial plane must be normal to axis.
       ! For the first vector try first the y-axis. If not good enough
       ! use the z-axis.

       if(abs(axis(2)) < 0.707107_realType) then
         vecR1(1) = zero
         vecR1(2) = one
         vecR1(3) = zero
       else
         vecR1(1) = zero
         vecR1(2) = zero
         vecR1(3) = one
       endif

       ! Make sure that vecR1 is normal to axis. Create a unit
       ! vector again.

       dot = vecR1(1)*axis(1) + vecR1(2)*axis(2) + vecR1(3)*axis(3)
       vecR1(1) = vecR1(1) - dot*axis(1)
       vecR1(2) = vecR1(2) - dot*axis(2)
       vecR1(3) = vecR1(3) - dot*axis(3)

       dot = one/sqrt(vecR1(1)**2 + vecR1(2)**2 + vecR1(3)**2)
       vecR1(1) = vecR1(1)*dot
       vecR1(2) = vecR1(2)*dot
       vecR1(3) = vecR1(3)*dot

       ! Create the second vector which spans the axial plane. This must
       ! be normal to both axis and vecR1, i.e. the cross-product.

       vecR2(1) = axis(2)*vecR1(3) - axis(3)*vecR1(2)
       vecR2(2) = axis(3)*vecR1(1) - axis(1)*vecR1(3)
       vecR2(3) = axis(1)*vecR1(2) - axis(2)*vecR1(1)

       end subroutine unitVectorsInAxialPlane
