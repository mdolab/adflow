
!
!      ******************************************************************
!      *                                                                *
!      * File:          setssMetircBwd.f90                              *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 11-03-2014                                      *
!      * Last modified: 11-03-2014                                      *
!      *                                                                *
!      ******************************************************************

subroutine setssMetricBwd(nn, ss)
  
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
     ss(1:je,1:ke,:) = si(1,1:je,1:ke,:)
    
     !=======================================================

  case (iMax)
     ss(1:je,1:ke,:) = si(il,1:je,1:ke,:)
    
     !=======================================================

  case (jMin)
     ss(0:ie,1:ke,:) = sj(0:ie,1,1:ke,:)
    
     !=======================================================

  case (jMax)
     ss(0:ie,1:ke,:) = sj(0:ie,jl,1:ke,:)
     
     !=======================================================

  case (kMin)
     ss(0:ie,1:je,:) = sk(0:ie,1:je,1,:)
     
     !=======================================================

  case (kMax)
     ss(0:ie,1:je,:) = sk(0:ie,1:je,kl,:)
     
  end select
end subroutine setssMetricBwd


