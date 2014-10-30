!
!      ******************************************************************
!      *                                                                *
!      * File:          setgammaBwd.f90                                 *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 10-15-2014                                      *
!      * Last modified: 10-21-2014                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine setgammaBwd(nn, gamma1, gamma2)

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
     gamma1(1:je,1:ke) = gamma(1, 1:je,1:ke)
     gamma2(1:je,1:ke) = gamma(2, 1:je,1:ke)
  case (iMax)
     gamma1(1:je,1:ke) = gamma(ie,1:je,1:ke)
     gamma2(1:je,1:ke) = gamma(il,1:je,1:ke)
  case (jMin)
     gamma1(1:ie,1:ke) = gamma(1:ie,1, 1:ke)
     gamma2(1:ie,1:ke) = gamma(1:ie,2, 1:ke)
  case (jMax)
     gamma1(1:ie,1:ke) = gamma(1:ie,je,1:ke)
     gamma2(1:ie,1:ke) = gamma(1:ie,jl,1:ke)
  case (kMin)
     gamma1(1:ie,1:je) = gamma(1:ie,1:je,1 ) 
     gamma2(1:ie,1:je) = gamma(1:ie,1:je,2 )
  case (kMax)
     gamma1(1:ie,1:je) = gamma(1:ie,1:je,ke) 
     gamma2(1:ie,1:je) = gamma(1:ie,1:je,kl)
  end select

end subroutine setgammaBwd

