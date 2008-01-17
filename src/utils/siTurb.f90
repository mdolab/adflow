!
!      ******************************************************************
!      *                                                                *
!      * File:          siTurb.f90                                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-14-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine siTurb(mass, len, time, temp, turbName, mult, trans)
!
!      ******************************************************************
!      *                                                                *
!      * siTurb computes the conversion from the given turbulence       *
!      * unit, which can be constructed from mass, len, time and temp,  *
!      * to the SI-unit for the given variable. The conversion will     *
!      * look like: var in SI = mult*(var in NCU) + trans.              *
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
       integer, intent(in)              :: mass, len, time, temp
       character(len=*), intent(in)     :: turbName
       real(kind=realType), intent(out) :: mult, trans
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the situation we are having here.

       if(mass == Kilogram .and. len  == Meter .and. &
          time == Second   .and. temp == Kelvin) then

         ! Everthing is already in SI units. No conversion needed.

         mult  = one
         trans = zero

       else

         call terminate("siTurb", &
                        "No idea how to convert this to SI units")

       endif

       end subroutine siTurb
