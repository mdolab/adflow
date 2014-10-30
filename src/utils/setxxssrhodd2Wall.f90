
!
!      ******************************************************************
!      *                                                                *
!      * File:          setxxssrhodd2Wall.f90                           *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 10-28-2014                                      *
!      * Last modified: 10-28-2014                                      *
!      *                                                                *
!      ******************************************************************

subroutine setxxssrhodd2Wall(nn, xx, ss, rho1, rho2, dd2Wall)

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
  ! Determine the face id on which the subface is located and set
  ! the pointers accordinly.
  select case (BCFaceID(nn))

  case (iMin)
     rho2 => w(2,1:,1:,irho); rho1 => w(1,1:,1:,irho)
     ss   => si(1,:,:,:);     xx   => x(1,:,:,:)

     if(equations == RANSEquations) dd2Wall => d2Wall(2,:,:)

     !===========================================================

  case (iMax)
     rho2 => w(il,1:,1:,irho); rho1 => w(ie,1:,1:,irho)
     ss   => si(il,:,:,:);     xx   => x(il,:,:,:)

     if(equations == RANSEquations) dd2Wall => d2Wall(il,:,:)

     !===========================================================

  case (jMin)
     rho2 => w(1:,2,1:,irho); rho1 => w(1:,1,1:,irho)
     ss   => sj(:,1,:,:);     xx   => x(:,1,:,:)

     if(equations == RANSEquations) dd2Wall => d2Wall(:,2,:)

     !===========================================================

  case (jMax)
     rho2 => w(1:,jl,1:,irho); rho1 => w(1:,je,1:,irho)
     ss   => sj(:,jl,:,:);     xx   => x(:,jl,:,:)

     if(equations == RANSEquations) dd2Wall => d2Wall(:,jl,:)

     !===========================================================

  case (kMin)
     rho2 => w(1:,1:,2,irho); rho1 => w(1:,1:,1,irho)
     ss   => sk(:,:,1,:);     xx   => x(:,:,1,:)

     if(equations == RANSEquations) dd2Wall => d2Wall(:,:,2)

     !===========================================================

  case (kMax)
     rho2 => w(1:,1:,kl,irho); rho1 => w(1:,1:,ke,irho)
     ss   => sk(:,:,kl,:);     xx   => x(:,:,kl,:)

     if(equations == RANSEquations) dd2Wall => d2Wall(:,:,kl)

  end select

end subroutine setxxssrhodd2Wall


