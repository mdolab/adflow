module commonFormats
!
!       This module contains common format specifiers for printing output.
!       The convention for adding a format variable in this module or elsewhere is as follows:
!
!       commonFormats module variable: formats spanning multiple modules
!       other module variable: formats spanning multiple subroutines within a module
!       subroutine variable: formats that are repeated within a subroutine
!       string: one-off cases
!
    use constants, only: maxStringLen
    implicit none
    save

    ! The * in these patterns means that it is applied as often as there are input arguments matching this pattern.
    ! For example, (*(A, 1X)) means that this formatting pattern is applied N times if you have N strings that
    ! you want to be separated by a space.
    ! Similarly, (*(A, ES12.5)) expects N (string, float) pairs as an input.

    ! Strings
    character(len=maxStringLen) :: strings = '(*(A))'

    ! Strings followed by one space
    character(len=maxStringLen) :: stringSpace = '(*(A, 1X))'

    ! Strings followed by a number in scientific notation with 5 decimal places
    character(len=maxStringLen) :: stringSci5 = '(*(A, ES12.5))'

    ! Strings followed by a one-digit integer
    character(len=maxStringLen) :: stringInt1 = '(*(A, I1))'

    ! Numbers in scientific notation with 12 decimal places
    character(len=maxStringLen) :: sci12 = '(*(ES20.12))'

    ! Integers written with 5 characters
    character(len=maxStringLen) :: int5 = '(*(I5))'

end module commonFormats
