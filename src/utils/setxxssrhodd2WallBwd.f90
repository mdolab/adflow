
!
!      ******************************************************************
!      *                                                                *
!      * File:          setxxssrhodd2WallBwd.f90                        *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 10-28-2014                                      *
!      * Last modified: 10-28-2014                                      *
!      *                                                                *
!      ******************************************************************

subroutine setxxssrhodd2WallBwd(nn, xx, ss, rho1, rho2, dd2Wall)

  use BCTypes
  use blockPointers
  use flowVarRefState
  use inputPhysics
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  real(kind=realType), dimension(imaxDim,jmaxDim) :: rho2, rho1
  real(kind=realType), dimension(imaxDim,jmaxDim) :: dd2Wall
  real(kind=realType), dimension(imaxDim,jmaxDim,3) :: ss, xx
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
     rho2(1:je,1:ke) = w(2,1:je,1:ke,irho); rho1(1:je,1:ke) = w(1,1:je,1:ke,irho)
     ss(1:je,1:ke,:)  = si(1,1:je,1:ke,:);     xx(1:je,1:ke,:) = x(1,1:je,1:ke,:)

     if(equations == RANSEquations) dd2Wall(1:je,1:ke) = d2Wall(2,1:je,1:ke)

     !===========================================================

  case (iMax)
     rho2(1:je,1:ke) = w(il,1:je,1:ke,irho); rho1(1:je,1:ke) = w(ie,1:je,1:ke,irho)
     ss(1:je,1:ke,:) = si(il,1:je,1:ke,:);     xx(1:je,1:ke,:) = x(il,1:je,1:ke,:)

     if(equations == RANSEquations) dd2Wall(1:je,1:ke) = d2Wall(il,1:je,1:ke)

     !===========================================================

  case (jMin)
     rho2(1:ie,1:ke) = w(1:ie,2,1:ke,irho); rho1(1:ie,1:ke) = w(1:ie,1,1:ke,irho)
     ss(1:ie,1:ke,:) = sj(1:ie,1,1:ke,:);     xx(1:ie,1:ke,:) = x(1:ie,1,1:ke,:)

     if(equations == RANSEquations) dd2Wall(1:ie,1:ke) = d2Wall(1:ie,2,1:ke)

     !===========================================================

  case (jMax)
     rho2(1:ie,1:ke) = w(1:ie,jl,1:ke,irho); rho1(1:ie,1:ke) = w(1:ie,je,1:ke,irho)
     ss(1:ie,1:ke,:) = sj(1:ie,jl,1:ke,:);     xx(1:ie,1:ke,:) = x(1:ie,jl,1:ke,:)

     if(equations == RANSEquations) dd2Wall(1:ie,1:ke) = d2Wall(1:ie,jl,1:ke)

     !===========================================================

  case (kMin)
     rho2(1:ie,1:je) = w(1:ie,1:je,2,irho); rho1(1:ie,1:je) = w(1:ie,1:je,1,irho)
     ss(1:ie,1:je,:) = sk(1:ie,1:je,1,:);     xx(1:ie,1:je,:) = x(1:ie,1:je,1,:)

     if(equations == RANSEquations) dd2Wall(1:ie,1:je) = d2Wall(1:ie,1:je,2)

     !===========================================================

  case (kMax)
     rho2(1:ie,1:je) = w(1:ie,1:je,kl,irho); rho1(1:ie,1:je) = w(1:ie,1:je,ke,irho)
     ss(1:ie,1:je,:) = sk(1:ie,1:je,kl,:);     xx(1:ie,1:je,:) = x(1:ie,1:je,kl,:)

     if(equations == RANSEquations) dd2Wall(1:ie,1:je) = d2Wall(1:ie,1:je,kl)

  end select

end subroutine setxxssrhodd2WallBwd


