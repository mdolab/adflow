
!
!      ******************************************************************
!      *                                                                *
!      * File:          resetallxBwd.f90                                *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 11-03-2014                                      *
!      * Last modified: 11-03-2014                                      *
!      *                                                                *
!      ******************************************************************

subroutine resetallxBwd(nn, x0, x1, x2)                 

  use BCTypes
  use blockPointers
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  real(kind=realType), dimension(imaxDim+1,jmaxDim+1,3) :: x0, x1, x2
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
     x(0,0:je,0:ke,:) = x0(1:je+1,1:ke+1,:)
     x(1,0:je,0:ke,:) = x1(1:je+1,1:ke+1,:)
     x(2,0:je,0:ke,:) = x2(1:je+1,1:ke+1,:)

     !===========================================================

  case (iMax)
     x(ie,0:je,0:ke,:) = x0(1:je+1,1:ke+1,:)
     x(il,0:je,0:ke,:) = x1(1:je+1,1:ke+1,:)
     x(nx,0:je,0:ke,:) = x2(1:je+1,1:ke+1,:)

     !===========================================================

  case (jMin)
     x(0:ie,0,0:ke,:) = x0(1:ie+1,1:ke+1,:)
     x(1:ie,0,0:ke,:) = x1(1:ie+1,1:ke+1,:)
     x(0:ie,2,0:ke,:) = x2(1:ie+1,1:ke+1,:)
 
     !===========================================================

  case (jMax)
     x(0:ie,je,0:ke,:) = x0(1:ie+1,1:ke+1,:)
     x(0:ie,jl,0:ke,:) = x1(1:ie+1,1:ke+1,:)
     x(0:ie,ny,0:ke,:) = x2(1:ie+1,1:ke+1,:)
 
     !===========================================================

  case (kMin)
     x(0:ie,0:je,0,:) = x0(1:ie+1,1:je+1,:)
     x(0:ie,0:je,1,:) = x1(1:ie+1,1:je+1,:)
     x(0:ie,0:je,2,:) = x2(1:ie+1,1:je+1,:)

     !===========================================================

  case (kMax)
      x(0:ie,0:je,ke,:) = x0(1:ie+1,1:je+1,:)
      x(0:ie,0:je,kl,:) = x1(1:ie+1,1:je+1,:)
      x(0:ie,0:je,nz,:) = x2(1:ie+1,1:je+1,:)

  end select
end subroutine resetallxBwd


