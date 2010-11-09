!
!      ******************************************************************
!      *                                                                *
!      * File:          bcEulerWallAdj.f90                              *
!      * Author:        Edwin van der Weide                             *
!      *                Seongim Choi,C.A.(Sandy) Mader                  *
!      * Starting date: 03-21-2006                                      *
!      * Last modified: 06-09-2008                                      *
!      *                                                                *
!      ******************************************************************
!

subroutine bcEulerWallForceCouplingAdj(wAdj,pAdj)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * bcEulerWallAdj applies inviscid wall bcoundary condition       *
  !      * to the small set of cells wAdj (2,2,2) cube.                   *
  !      * This function is based on BCEulerWarll.f90 in src/solver       *
  !      * It only works for constant and linear pressure extrapolation   *
  !      *                                                                *
  !      ******************************************************************
  !


  use flowVarRefState
  use inputDiscretization

  implicit none

  ! Subroutine Arguments

  real(kind=realType), intent(in) :: wAdj(2,2,2,nw)
  real(kind=realType), intent(inout) :: pAdj(3,2,2)

  ! Local variables.

  integer(kind=intType) :: i, j, k

  BCTreatment: select case (wallBCTreatment)

  case (constantPressure)

     do j=1,2
        do i=1,2
           pAdj(1,i,j) = zero
        enddo
     enddo

  case (linExtrapolPressure)
     do j=1,2
        do i=1,2
           pAdj(1,i,j) = pAdj(3,i,j) - pAdj(2,i,j)
        enddo
     enddo
     
  case (quadExtrapolPressure)
     call terminate("bcEulerWallForceCouplingAdj", "Adjoint Quadratic extrapolation not implemented")
  case (normalMomentum)
     call terminate("bcEulerWallForceCouplingAdj", "Normal Momentum not implemented")
  end select BCTreatment

  ! Determine the state in the halo cell. Again loop over
  ! the cell range for this subface.

  do j=1,2
     do i=1,2
        pAdj(1,i,j) = max(zero, pAdj(2,i,j)-pAdj(1,i,j))
     enddo
  enddo

end subroutine bcEulerWallForceCouplingAdj
