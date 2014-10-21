
!
!      ******************************************************************
!      *                                                                *
!      * File:          setssBwd.f90                                    *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 10-21-2014                                      *
!      * Last modified: 10-21-2014                                      *
!      *                                                                *
!      ******************************************************************

subroutine setssBwd(nn, ssi, ssj, ssk, ss)
  
  use BCTypes
  use blockPointers
  use flowVarRefState
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  real(kind=realType), dimension(imaxDim,jmaxDim,3) :: ssi, ssj, ssk
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
     ssi(1:je,1:ke,:) = si(1,1:je,1:ke,:)
     ssj(0:je,1:ke,:) = sj(2,0:je,1:ke,:)
     ssk(1:je,0:ke,:) = sk(2,1:je,0:ke,:)

     if( addGridVelocities ) ss(1:je,1:ke,:) = s(2,1:je,1:ke,:)

     !=======================================================

  case (iMax)
     ssi(1:je,1:ke,:) = si(il,1:je,1:ke,:)
     ssj(0:je,1:ke,:) = sj(il,0:je,1:ke,:)
     ssk(1:je,0:ke,:) = sk(il,1:je,0:ke,:)

     if( addGridVelocities ) ss(1:je,1:ke,:) = s(il,1:je,1:ke,:)

     !=======================================================

  case (jMin)
     ssi(0:ie,1:ke,:) = sj(0:ie,1,1:ke,:)
     ssj(1:ie,1:ke,:) = si(1:ie,2,1:ke,:)
     ssk(1:ie,0:ke,:) = sk(1:ie,2,0:ke,:)

     if( addGridVelocities ) ss(1:ie,1:ke,:) = s(1:ie,2,1:ke,:)

     !=======================================================

  case (jMax)
     ssi(0:ie,1:ke,:) = sj(0:ie,jl,1:ke,:)
     ssj(1:ie,1:ke,:) = si(1:ie,jl,1:ke,:)
     ssk(1:ie,0:ke,:) = sk(1:ie,jl,0:ke,:)

     if( addGridVelocities ) ss(1:ie,1:ke,:) = s(1:ie,jl,1:ke,:)

     !=======================================================

  case (kMin)
     ssi(0:ie,1:je,:) = sk(0:ie,1:je,1,:)
     ssj(1:ie,0:je,:) = si(1:ie,0:je,2,:)
     ssk(1:ie,1:je,:) = sj(1:ie,1:je,2,:)

     if( addGridVelocities ) ss(1:ie,1:je,:) = s(1:ie,1:je,2,:)

     !=======================================================

  case (kMax)
     ssi(0:ie,1:je,:) = sk(0:ie,1:je,kl,:)
     ssj(1:ie,0:je,:) = si(1:ie,0:je,kl,:)
     ssk(1:ie,1:je,:) = sj(1:ie,1:je,kl,:)

     if( addGridVelocities ) ss(1:ie,1:je,:) = s(1:ie,1:je,kl,:)

  end select
end subroutine setssBwd


