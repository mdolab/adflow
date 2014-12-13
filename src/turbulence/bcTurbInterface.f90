!
!      ******************************************************************
!      *                                                                *
!      * File:          bcTurbInterface.f90                             *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 01-09-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine bcTurbInterface(nn)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * bcTurbInterface applies the halo treatment for interface halo  *
  !      * cells, sliding mesh interface and domain interface. As these   *
  !      * are not really boundary conditions, the variable bvt is simply *
  !      * set to keep the current value.                                 *
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
  !
  ! Loop over the faces of the subfaces and set the values of
  ! bvt to keep the current value.

  do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
     do i=BCData(nn)%icBeg, BCData(nn)%icEnd
        do l=nt1,nt2
           select case (BCFaceID(nn))
           case (iMin)
              bvti1(i,j,l) = w(1,i,j,l)
           case (iMax)
              bvti2(i,j,l) = w(ie,i,j,l)
           case (jMin)
              bvtj1(i,j,l) = w(i,1,j,l)
           case (jMax)
              bvtj2(i,j,l) = w(i,je,j,l)
           case (kMin)
              bvtk1(i,j,l) = w(i,j,1,l)
           case (kMax)
              bvtk2(i,j,l) = w(i,j,ke,l)
           end select
        enddo
     enddo
  enddo

  ! Note that the original code had an error in the pointers...they
  ! were pointing to {il,jl,kl} and not {ie, je, ke}.

end subroutine bcTurbInterface
