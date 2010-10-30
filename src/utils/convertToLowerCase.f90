!
!      ******************************************************************
!      *                                                                *
!      * File:          convertToLowerCase.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-12-2002                                      *
!      * Last modified: 03-23-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine convertToLowerCase(string)
!
!      ******************************************************************
!      *                                                                *
!      * convertToLowerCase converts the given string to lower case.    *
!      *                                                                *
!      ******************************************************************
!
       implicit none
!
!      Subroutine arguments
!
       character (len=*), intent(inout) :: string
!
!      Local variables
!
       integer, parameter :: upperToLower = iachar("a") - iachar("A")

       integer :: i, lenString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the length of the given string and convert the upper
       ! case characters to lower case.

       lenString = len_trim(string)
       do i=1,lenString
         if("A" <= string(i:i) .and. string(i:i) <= "Z")    &
         string(i:i) = achar(iachar(string(i:i)) + upperToLower)
       enddo

       end subroutine convertToLowerCase
