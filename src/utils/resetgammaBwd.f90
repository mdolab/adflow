!
!      ******************************************************************
!      *                                                                *
!      * File:          resetgammaBwd.f90                               *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 10-21-2014                                      *
!      * Last modified: 10-21-2014                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine resetgammaBwd(nn, gamma1, gamma2)

  use BCTypes
  use blockPointers
  use flowVarRefState
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  real(kind=realType), dimension(imaxDim,jmaxDim) :: gamma1, gamma2  

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Determine the face id on which the subface is located and set
  ! the pointers accordinly.

  select case (BCFaceID(nn))
  case (iMin)
      gamma(1, 1:je,1:ke) = gamma1(1:je,1:ke)
      gamma(2, 1:je,1:ke) = gamma2(1:je,1:ke)
  case (iMax)
      gamma(ie,1:je,1:ke) = gamma1(1:je,1:ke)
      gamma(il,1:je,1:ke) = gamma2(1:je,1:ke)
  case (jMin)
      gamma(1:ie,1, 1:ke) = gamma1(1:ie,1:ke)
      gamma(1:ie,2, 1:ke) = gamma2(1:ie,1:ke)
  case (jMax)
      gamma(1:ie,je,1:ke) = gamma1(1:ie,1:ke)
      gamma(1:ie,jl,1:ke) = gamma2(1:ie,1:ke)
  case (kMin)
      gamma(1:ie,1:je,1 ) = gamma1(1:ie,1:je)
      gamma(1:ie,1:je,2 ) = gamma2(1:ie,1:je)
  case (kMax)
      gamma(1:ie,1:je,ke) = gamma1(1:ie,1:je)
      gamma(1:ie,1:je,kl) = gamma2(1:ie,1:je)
  end select

end subroutine resetgammaBwd

