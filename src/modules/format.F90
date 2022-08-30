module format
!
!       This module contains common format specifiers for printing output.
!
       use constants, only : maxStringLen
       implicit none
       save

       ! An arbitrary number of strings
       character(len=maxStringLen) :: strings = '(*(A))'

       ! Strings followed by one space
       character(len=maxStringLen) :: stringSpace = '(*(A, 1X))'

       ! Strings followed by a number in scientific notation with 5 decimal places
       character(len=maxStringLen) :: stringSci5 = '(*(A, ES12.5))'

       ! Strings followed by a one-digit integer
       character(len=maxStringLen) :: stringInt1 = '(*(A, I1))'

end module format
