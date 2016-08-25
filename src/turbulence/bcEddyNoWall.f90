!
!      ******************************************************************
!      *                                                                *
!      * File:          bcEddyNoWall.f90                                *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 06-11-2003                                      *
!      * Last modified: 04-11-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine bcEddyNoWall(nn)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * bcEddyNoWall sets the eddy viscosity in the halo cells of      *
  !      * subface nn of the block given in blockPointers. The boundary   *
  !      * condition on the subface can be anything but a viscous wall.   *
  !      * A homogeneous neumann condition is applied, which means that   *
  !      * the eddy viscosity is simply copied from the interior cell.    *
  !      *                                                                *
  !      ******************************************************************
  !
  use constants
  use blockPointers
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Determine the face id on which the subface and copy

  select case (BCFaceid(nn))
  case (iMin)
     do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
        do i=BCData(nn)%icBeg, BCData(nn)%icEnd
           rev(1,i,j) = rev(2,i,j)
        enddo
     enddo

  case (iMax)
     do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
        do i=BCData(nn)%icBeg, BCData(nn)%icEnd
           rev(ie,i,j) = rev(il,i,j)
        enddo
     enddo

  case (jMin)
     do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
        do i=BCData(nn)%icBeg, BCData(nn)%icEnd
           rev(i,1,j) = rev(i,2,j)
        enddo
     enddo

  case (jMax)
     do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
        do i=BCData(nn)%icBeg, BCData(nn)%icEnd
           rev(i,je,j) = rev(i,jl,j)
        enddo
     enddo

  case (kMin)
     do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
        do i=BCData(nn)%icBeg, BCData(nn)%icEnd
           rev(i,j,1) = rev(i,j,2)
        enddo
     enddo

  case (kMax)
     do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
        do i=BCData(nn)%icBeg, BCData(nn)%icEnd
           rev(i,j,ke) = rev(i,j,kl)
        enddo
     enddo
  end select

end subroutine bcEddyNoWall
