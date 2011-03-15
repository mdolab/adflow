!
!      ******************************************************************
!      *                                                                *
!      * File:          siLen.f90                                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-05-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine siLen(len, mult, trans)
!
!      ******************************************************************
!      *                                                                *
!      * siLen computes the conversion from the given length unit to    *
!      * the SI-unit meter. The conversion will look like:              *
!      * length in meter = mult*(length in NCU) + trans.                *
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
       integer, intent(in)              :: len
       real(kind=realType), intent(out) :: mult, trans
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the situation we are having here.

       select case (len)

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
           call terminate("siLen", &
                          "No idea how to convert this to SI units")

       end select

       end subroutine siLen
