!
!      ******************************************************************
!      *                                                                *
!      * File:          siTemperature.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-10-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine siTemperature(temp, mult, trans)
!
!      ******************************************************************
!      *                                                                *
!      * siTemperature computes the conversion from the given           *
!      * temperature unit to the SI-unit kelvin. The conversion will    *
!      * look like:                                                     *
!      * temperature in K = mult*(temperature in NCU) + trans.          *
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
       integer, intent(in)              :: temp
       real(kind=realType), intent(out) :: mult, trans
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the situation we are having here.

       select case (temp)
 
         case (Kelvin)

           ! Temperature is already given in Kelvin. No need to convert.

           mult  = one
           trans = zero

         case (Celcius)      ! is it Celcius or Celsius?

           ! Temperature is in Celsius. Only an offset must be applied.

           mult  = one
           trans = 273.16_realType

         case (Rankine)

           ! Temperature is in Rankine. Only a multiplication needs to
           ! be performed.

           mult  = 5.0_realType/9.0_realType
           trans = zero

         case (Fahrenheit)

           ! Temperature is in Fahrenheit. Both a multiplication and an
           ! offset must be applied.

           mult  = 5.0_realType/9.0_realType
           trans = 255.382

         case default

           ! Unknown temperature unit.

           call terminate("siTemperature", &
                          "No idea how to convert this to SI units")

       end select

       end subroutine siTemperature
