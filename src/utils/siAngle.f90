!
!      ******************************************************************
!      *                                                                *
!      * File:          siAngle.f90                                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-11-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine siAngle(angle, mult, trans)
!
!      ******************************************************************
!      *                                                                *
!      * siAngle computes the conversion from the given angle unit to   *
!      * the SI-unit radian. The conversion will look like:             *
!      * angle in radians = mult*(angle in NCU) + trans.                *
!      * NCU means non-christian units, i.e. everything that is not SI. *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use su_cgns
       implicit none
!
!      Subroutine arguments.
!
       integer, intent(in)              :: angle
       real(kind=realType), intent(out) :: mult, trans
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the situation we are having here.

       if(angle == Radian) then

         ! Angle is already given in radIans. No need for a conversion.

         mult  = one
         trans = zero

       else if(angle == Degree) then

         ! Angle is given in degrees. A multiplication must be performed.

         mult  = pi/180.0_realType
         trans = zero

       else

         call terminate("siAngle", &
                        "No idea how to convert this to SI units")

       endif

       end subroutine siAngle
