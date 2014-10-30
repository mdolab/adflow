!
!      ******************************************************************
!      *                                                                *
!      * File:          setgamma.f90                                    *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 10-15-2014                                      *
!      * Last modified: 10-21-2014                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine setgamma(nn, gamma1, gamma2)

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

  select case (BCFaceID(nn))
  case (iMin)
     gamma1 => gamma(1, 1:,1:); gamma2 => gamma(2, 1:,1:)
  case (iMax)
     gamma1 => gamma(ie,1:,1:); gamma2 => gamma(il,1:,1:)
  case (jMin)
     gamma1 => gamma(1:,1, 1:); gamma2 => gamma(1:,2, 1:)
  case (jMax)
     gamma1 => gamma(1:,je,1:); gamma2 => gamma(1:,jl,1:)
  case (kMin)
     gamma1 => gamma(1:,1:,1 ); gamma2 => gamma(1:,1:,2 )
  case (kMax)
     gamma1 => gamma(1:,1:,ke); gamma2 => gamma(1:,1:,kl)
  end select

end subroutine setgamma





