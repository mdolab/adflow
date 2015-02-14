!
!      ******************************************************************
!      *                                                                *
!      * File:          terminate.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-10-2002                                      *
!      * Last modified: 10-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine terminate(routineName, errorMessage)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * terminate writes an error message to standard output and       *
  !      * terminates the execution of the program.                       *
  !      *                                                                *
  !      ******************************************************************
  !
  use precision
  use communication
  use constants
  implicit none
  !
  !      Subroutine arguments
  !
  character(len=*), intent(in) :: routineName
  character(len=*), intent(in) :: errorMessage

end subroutine terminate
