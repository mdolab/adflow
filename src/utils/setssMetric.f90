
!
!      ******************************************************************
!      *                                                                *
!      * File:          setssMetric.f90                                 *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 11-03-2014                                      *
!      * Last modified: 11-03-2014                                      *
!      *                                                                *
!      ******************************************************************

subroutine setssMetric(nn, ss)

  use BCTypes
  use blockPointers
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  real(kind=realType), dimension(:,:,:), pointer :: ss

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
     ss => si(1,:,:,:)
    
     !=======================================================

  case (iMax)
     ss => si(il,:,:,:)
    
     !=======================================================

  case (jMin)
     ss => sj(:,1,:,:)

     !=======================================================

  case (jMax)
     ss => sj(:,jl,:,:)
     
     !=======================================================

  case (kMin)
     ss => sk(:,:,1,:)
    
     !=======================================================

  case (kMax)
     ss => sk(:,:,kl,:)
    
  end select
end subroutine setssMetric


