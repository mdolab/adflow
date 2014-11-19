
!
!      ******************************************************************
!      *                                                                *
!      * File:          resetxxssrhodd2WallBwd.f90                      *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 10-28-2014                                      *
!      * Last modified: 10-28-2014                                      *
!      *                                                                *
!      ******************************************************************

subroutine resetxxssrhodd2WallBwd(nn, xx, ss, rho1, rho2, dd2Wall)

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
  real(kind=realType), dimension(imaxDim-2,jmaxDim-2) :: dd2Wall
  real(kind=realType), dimension(imaxDim,jmaxDim,3) :: ss
  real(kind=realType), dimension(imaxDim+1,jmaxDim+1,3) :: xx
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
      w(2,1:je,1:ke,irho) = rho2(1:je,1:ke)
      w(1,1:je,1:ke,irho) = rho1(1:je,1:ke)
      si(1,1:je,1:ke,:) = ss(1:je,1:ke,:)
      x(1,0:je,0:ke,:) = xx(1:je+1,1:ke+1,:)

     if(equations == RANSEquations)  d2Wall(2,2:jl,2:kl) = dd2Wall(1:jl-1,1:kl-1)

     !===========================================================

  case (iMax)
      w(il,1:je,1:ke,irho) = rho2(1:je,1:ke)
      w(ie,1:je,1:ke,irho) = rho1(1:je,1:ke)
      si(il,1:je,1:ke,:)  = ss(1:je,1:ke,:)
      x(il,0:je,0:ke,:) = xx(1:je+1,1:ke+1,:)

     if(equations == RANSEquations)  d2Wall(il,2:jl,2:kl) = dd2Wall(1:jl-1,1:kl-1)

     !===========================================================

  case (jMin)
      w(1:ie,2,1:ke,irho) = rho2(1:ie,1:ke)
      w(1:ie,1,1:ke,irho) = rho1(1:ie,1:ke)
      sj(1:ie,1,1:ke,:) = ss(1:ie,1:ke,:)
      x(0:ie,1,0:ke,:) = xx(1:ie+1,1:ke+1,:)
     
     if(equations == RANSEquations)  d2Wall(2:il,2,2:kl) = dd2Wall(1:il-1,1:kl-1)

     !===========================================================

  case (jMax)
      w(1:ie,jl,1:ke,irho) = rho2(1:ie,1:ke)
      w(1:ie,je,1:ke,irho) = rho1(1:ie,1:ke)
      sj(1:ie,jl,1:ke,:) = ss(1:ie,1:ke,:)
      x(0:ie,jl,0:ke,:) = xx(1:ie+1,1:ke+1,:)

     if(equations == RANSEquations)  d2Wall(2:il,jl,2:kl) = dd2Wall(1:il-1,1:kl-1)

     !===========================================================

  case (kMin)
      w(1:ie,1:je,2,irho) = rho2(1:ie,1:je)
      w(1:ie,1:je,1,irho) = rho1(1:ie,1:je)
      sk(1:ie,1:je,1,:) = ss(1:ie,1:je,:)
      x(0:ie,0:je,1,:) = xx(1:ie+1,1:je+1,:)

     if(equations == RANSEquations)  d2Wall(2:il,2:jl,2) = dd2Wall(1:il-1,1:jl-1)

     !===========================================================

  case (kMax)
      w(1:ie,1:je,kl,irho) = rho2(1:ie,1:je)
      w(1:ie,1:je,ke,irho) = rho1(1:ie,1:je)
      sk(1:ie,1:je,kl,:) = ss(1:ie,1:je,:)
      x(0:ie,0:je,kl,:) = xx(1:ie+1,1:je+1,:)

     if(equations == RANSEquations)  d2Wall(2:il,2:jl,kl) = dd2Wall(1:il-1,1:jl-1)

  end select

end subroutine resetxxssrhodd2WallBwd


