
!
!      ******************************************************************
!      *                                                                *
!      * File:          setallx.f90                                     *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 11-03-2014                                      *
!      * Last modified: 11-03-2014                                      *
!      *                                                                *
!      ******************************************************************

subroutine setallx(nn, x0, x1, x2)

  use BCTypes
  use blockPointers
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  real(kind=realType), dimension(:,:,:), pointer :: x0, x1, x2
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
     x0   => x(0,:,:,:)
     x1   => x(1,:,:,:)
     x2   => x(2,:,:,:)

     !===========================================================

  case (iMax)
     x0   => x(ie,:,:,:)
     x1   => x(il,:,:,:)
     x2   => x(nx,:,:,:)

     !===========================================================

  case (jMin)
     x0   => x(:,0,:,:)
     x1   => x(:,1,:,:)
     x2   => x(:,2,:,:)

     !===========================================================

  case (jMax)
     x0   => x(:,je,:,:)
     x1   => x(:,jl,:,:)
     x2   => x(:,ny,:,:)

     !===========================================================

  case (kMin)
     x0   => x(:,:,0,:)
     x1   => x(:,:,1,:)
     x2   => x(:,:,2,:)

     !===========================================================

  case (kMax)
     x0   => x(:,:,ke,:)
     x1   => x(:,:,kl,:)
     x2   => x(:,:,nz,:)

  end select

end subroutine setallx


