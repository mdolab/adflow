
!
!      ******************************************************************
!      *                                                                *
!      * File:          setallxBwd.f90                                  *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 11-03-2014                                      *
!      * Last modified: 11-03-2014                                      *
!      *                                                                *
!      ******************************************************************

subroutine setallxBwd(nn, x0, x1, x2)                 

  use BCTypes
  use blockPointers
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  real(kind=realType), dimension(0:imaxDim,0:jmaxDim,3) :: x0, x1, x2
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Determine the face id on which the subface is located and set
  ! the pointers accordinly.
  x0 = zero
  x1 = zero
  x2 = zero

  select case (BCFaceID(nn))

  case (iMin)
     x0(0:je,0:ke,:) = x(0,0:je,0:ke,:)
     x1(0:je,0:ke,:) = x(1,0:je,0:ke,:)
     x2(0:je,0:ke,:) = x(2,0:je,0:ke,:)

     !===========================================================

  case (iMax)
     x0(0:je,0:ke,:) = x(ie,0:je,0:ke,:)
     x1(0:je,0:ke,:) = x(il,0:je,0:ke,:)
     x2(0:je,0:ke,:) = x(nx,0:je,0:ke,:)

     !===========================================================

  case (jMin)
     x0(0:ie,0:ke,:) = x(0:ie,0,0:ke,:)
     x1(0:ie,0:ke,:) = x(0:ie,1,0:ke,:)
     x2(0:ie,0:ke,:) = x(0:ie,2,0:ke,:)
 
     !===========================================================

  case (jMax)
     x0(0:ie,0:ke,:) = x(0:ie,je,0:ke,:)
     x1(0:ie,0:ke,:) = x(0:ie,jl,0:ke,:)
     x2(0:ie,0:ke,:) = x(0:ie,ny,0:ke,:)
 
     !===========================================================

  case (kMin)
     x0(0:ie,0:je,:) = x(0:ie,0:je,0,:)
     x1(0:ie,0:je,:) = x(0:ie,0:je,1,:)
     x2(0:ie,0:je,:) = x(0:ie,0:je,2,:)

     !===========================================================

  case (kMax)
     x0(0:ie,0:je,:) = x(0:ie,0:je,ke,:)
     x1(0:ie,0:je,:) = x(0:ie,0:je,kl,:)
     x2(0:ie,0:je,:) = x(0:ie,0:je,nz,:)

  end select
end subroutine setallxBwd


