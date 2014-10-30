
!
!      ******************************************************************
!      *                                                                *
!      * File:          resetxxssrhodd2Wall.f90                         *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 10-28-2014                                      *
!      * Last modified: 10-28-2014                                      *
!      *                                                                *
!      ******************************************************************

subroutine resetxxssrhodd2Wall(nn, xx, ss, rho1, rho2, dd2Wall)

  use BCTypes
  use blockPointers
  use flowVarRefState
  use inputPhysics
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  real(kind=realType), dimension(:,:),   pointer :: rho2, rho1
  real(kind=realType), dimension(:,:),   pointer :: dd2Wall
  real(kind=realType), dimension(:,:,:), pointer :: ss, xx
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !

  ! nullify(xx, ss, rho1, rho2, dd2Wall)

end subroutine resetxxssrhodd2Wall


