!
!      ******************************************************************
!      *                                                                *
!      * File:          replaceTabsAndReturns.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-28-2005                                      *
!      * Last modified: 04-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine replaceTabsAndReturns(string)
!
!      ******************************************************************
!      *                                                                *
!      * replaceTabsAndReturns replaces the tab and return characters   *
!      * in the given string by spaces.                                 *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
!
!      Subroutine arguments.
!
       character (len=*), intent(inout) :: string
!
!      Local variables.
!
       integer :: pos
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Replace the tabs.

       do
         pos = index(string, tabChar)
         if(pos == 0) exit
         string(pos:pos) = " "
       enddo

       ! Replace the returns.

       do
         pos = index(string, retChar)
         if(pos == 0) exit
         string(pos:pos) = " "
       enddo

       end subroutine replaceTabsAndReturns
