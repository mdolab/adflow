!
!      ******************************************************************
!      *                                                                *
!      * File:          siDensity.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-11-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine siDensity(mass, len, mult, trans)
!
!      ******************************************************************
!      *                                                                *
!      * siDensity computes the conversion from the given density       *
!      * unit, which can be constructed from mass and length, to the    *
!      * SI-unit kg/m^3. The conversion will look like:                 *
!      * density in kg/m^3 = mult*(density in NCU) + trans.             *
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
       integer, intent(in)              :: mass, len
       real(kind=realType), intent(out) :: mult, trans
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the situation we are having here.

       if(mass == Kilogram .and. len == Meter) then

         ! Density is given in kg/m^3, i.e. no need for a conversion.

         mult  = one
         trans = zero

       else

         call terminate("siDensity", &
                        "No idea how to convert this to SI units")

       endif

       end subroutine siDensity
