!
!      ******************************************************************
!      *                                                                *
!      * File:          resetgamma.f90                                  *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 10-15-2014                                      *
!      * Last modified: 10-21-2014                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine resetgamma(nn, gamma1, gamma2)

  use BCTypes
  use blockPointers
  use flowVarRefState
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  real(kind=realType), dimension(:,:),   pointer :: gamma1, gamma2   

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Determine the face id on which the subface is located and set
  ! the pointers accordinly.

  !nullify(gamma1, gamma2)

end subroutine resetgamma





