!
!      ******************************************************************
!      *                                                                *
!      * File:          bcTurbOutflow.f90                               *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 06-11-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine bcTurbOutflow(nn)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * bcTurbOutflow applies the implicit treatment of the outflow    *
  !      * boundary conditions to subface nn. As the outflow boundary     *
  !      * condition is independent of the turbulence model, either       *
  !      * extrapolation or zero Neumann, this routine is valid for all   *
  !      * models. It is assumed that the pointers in blockPointers are   *
  !      * already set to the correct block on the correct grid level.    *
  !      *                                                                *
  !      ******************************************************************
  !
  use constants
  use blockPointers
  use flowVarRefState
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, l
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Loop over the faces of the subfaces and set the values of bmt
  ! for an implicit treatment. For an outflow the turbulent variable
  ! variable is either extrapolated or zero Neumann. As constant
  ! extrapolation is used this leads to an identical treatment, i.e.
  ! the halo value is identical to the value of the internal cell.

  do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
     do i=BCData(nn)%icBeg, BCData(nn)%icEnd
        do l=nt1,nt2
           select case (BCFaceID(nn))
           case (iMin)
              bmti1(i,j,l,l) = -one
           case (iMax)
              bmti2(i,j,l,l) = -one
           case (jMin)
              bmtj1(i,j,l,l) = -one
           case (jMax)
              bmtj2(i,j,l,l) = -one
           case (kMin)
              bmtk1(i,j,l,l) = -one
           case (kMax)
              bmtk2(i,j,l,l) = -one
           end select
        enddo
     enddo
  enddo

end subroutine bcTurbOutflow
