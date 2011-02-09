!
!      ******************************************************************
!      *                                                                *
!      * File:          siPressure.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-10-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine siPressure(mass, len, time, mult, trans)
!
!      ******************************************************************
!      *                                                                *
!      * siPressure computes the conversion from the given pressure     *
!      * unit, which can be constructed from mass, length and time, to  *
!      * the SI-unit Pa. The conversion will look like:                 *
!      * pressure in Pa = mult*(pressure in NCU) + trans.               *
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
       integer, intent(in)              :: mass, len, time
       real(kind=realType), intent(out) :: mult, trans
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the situation we are having here.

       if(mass == Kilogram .and. len == Meter .and. time == Second) then

         ! Pressure is given in Pa, i.e. no need for a conversion.

         mult  = one
         trans = zero

       else

         call terminate("siPressure", &
                        "No idea how to convert this to SI units")

       endif

       end subroutine siPressure
