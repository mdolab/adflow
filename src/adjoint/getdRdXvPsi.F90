!
!     ******************************************************************
!     *                                                                *
!     * File:          getdRdXvPsi.F90                                 *
!     * Author:        Gaetan Kenway                                   *
!     * Starting date: 10-08-2010                                      *
!     * Last modified: 10-08-2010                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine getdRdXvPsi
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Multiply the current adjoint vector by dRdXv to get a vector   *
  !     * of length Xv. This is collective communication part            *
  !     *                                                                *
  !     ******************************************************************
  !
  use communication
  use ADjointPETSc, only: dRdX,psi
  use warpingPETSC 
  use inputADjoint
  implicit none

  integer(kind=intType) :: ierr


  call MatMultTranspose(dRdX,psi,sumbGridVec,ierr)

  call VecScatterBegin(cgnsTOsumbGrid,sumbGridVec,cgnsGridVec,INSERT_VALUES,SCATTER_REVERSE,ierr)
  call VecScatterEnd  (cgnsTOsumbGrid,sumbGridVec,cgnsGridVec,INSERT_VALUES,SCATTER_REVERSE,ierr)

  ! Now we have dRdXv*psi in cgnsGridVec and we can call getCGNSData to retrieve it

end subroutine getdRdXvPsi


