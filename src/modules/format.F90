module format
!
!       This module contains format specifiers for printing output.
!
       use constants, only : maxStringLen
       implicit none
       save

       ! Strings followed by one space
       character(len=maxStringLen) :: stringSpace = '(*(A, 1X))'

       ! Strings followed by a number in scientific notation with 5 decimal places
       character(len=maxStringLen) :: stringSci5 = '(*(A, ES12.5))'

end module format
