
!
!      ******************************************************************
!      *                                                                *
!      * File:          resetssMetricBwd.f90                            *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 11-03-2014                                      *
!      * Last modified: 11-03-2014                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine resetssMetricBwd(nn, ss)

  use BCTypes
  use blockPointers
  use flowVarRefState
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  real(kind=realType), dimension(imaxDim,jmaxDim,3) :: ss

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
      si(1,1:je,1:ke,:) = ss(1:je,1:ke,:)
     
     !=======================================================

  case (iMax)
      si(il,1:je,1:ke,:) = ss(1:je,1:ke,:)

     !=======================================================

  case (jMin)
      sj(1:ie,1,1:ke,:) = ss(1:ie,1:ke,:)

     !=======================================================

  case (jMax)
      sj(1:ie,jl,1:ke,:) = ss(1:ie,1:ke,:)

     !=======================================================

  case (kMin)
      sk(1:ie,1:je,1,:) = ss(1:ie,1:je,:)

     !=======================================================

  case (kMax)
      sk(1:ie,1:je,kl,:) = ss(1:ie,1:je,:)

  end select

end subroutine resetssMetricBwd


