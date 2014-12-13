!
!      ******************************************************************
!      *                                                                *
!      * File:          bcTurbInflow.f90                                *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 06-11-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine bcTurbInflow(nn)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * bcTurbInflow applies the implicit treatment of the inflow      *
  !      * boundary conditions to subface nn. As the inflow boundary      *
  !      * condition is independent of the turbulence model, this routine *
  !      * is valid for all models. It is assumed that the pointers in    *
  !      * blockPointers are already set to the correct block on the      *
  !      * correct grid level.                                            *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use BCTypes
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

  ! Loop over the faces of the subfaces and set the values of
  ! bvt and bmt such that the inflow state is linearly extrapolated
  ! with a fixed state at the face.

  do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
     do i=BCData(nn)%icBeg, BCData(nn)%icEnd

        ! Loop over the number of turbulent variables.

        do l=nt1,nt2
           select case (BCFaceID(nn))
           case (iMin)
              bvti1(i,j,l)   = two*BCData(nn)%turbInlet(i,j,l)
              bmti1(i,j,l,l) = one
           case (iMax)
              bvti2(i,j,l)   = two*BCData(nn)%turbInlet(i,j,l)
              bmti2(i,j,l,l) = one
           case (jMin)
              bvtj1(i,j,l)   = two*BCData(nn)%turbInlet(i,j,l)
              bmtj1(i,j,l,l) = one
           case (jMax)
              bvtj2(i,j,l)   = two*BCData(nn)%turbInlet(i,j,l)
              bmtj2(i,j,l,l) = one
           case (kMin)
              bvtk1(i,j,l)   = two*BCData(nn)%turbInlet(i,j,l)
              bmtk1(i,j,l,l) = one
           case (kMax)
              bvtk2(i,j,l)   = two*BCData(nn)%turbInlet(i,j,l)
              bmtk2(i,j,l,l) = one
           end select
        end do
     enddo
  enddo
end subroutine bcTurbInflow
