!
!      ******************************************************************
!      *                                                                *
!      * File:          siVelocity.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-27-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine siVelocity(length, time, mult, trans)
!
!      ******************************************************************
!      *                                                                *
!      * siVelocity computes the conversion from the given velocity     *
!      * unit, which can be constructed from length and time, to the    *
!      * SI-unit m/s. The conversion will look like:                    *
!      * velocity in m/s = mult*(velocity in ncu) + trans.              *
!      * Ncu means non-christian units, i.e. everything that is not SI. *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use su_cgns
       implicit none
!
!      Subroutine arguments.
!
       integer, intent(in)              :: length, time
       real(kind=realType), intent(out) :: mult, trans
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the situation we are having here.
       ! First the length.

       select case (length)

         case (Meter)
           mult = one; trans = zero

         case (CenTimeter)
           mult = 0.01_realType; trans = zero

         case (Millimeter)
           mult = 0.001_realType; trans = zero

         case (Foot)
           mult = 0.3048_realType; trans = zero

         case (Inch)
           mult = 0.0254_realType; trans = zero

         case default
           call terminate("siVelocity", &
                          "No idea how to convert this length to &
                           &SI units")

       end select

       ! And the time.

       select case (time)

         case (Second)
           mult = mult

         case default
           call terminate("siVelocity", &
                          "No idea how to convert this time to &
                           &SI units")

       end select

       end subroutine siVelocity
