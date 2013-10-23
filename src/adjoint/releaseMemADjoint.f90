!
!     ******************************************************************
!     *                                                                *
!     * File:          releaseMemADjoint.f90                           *
!     * Author:        Andre C. Marta                                  *
!     * Starting date: 07-20-2006                                      *
!     * Last modified: 02-07-2007                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine releaseMemADjoint()
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Release all memory allocated for the adjoint solver by         *
  !     * "allocMemADjoint".                                             *
  !     *                                                                *
  !     * This routine was based on:                                     *
  !     * /utils/releaseMemory.f90                                       *
  !     *                                                                *
  !     ******************************************************************
  !
  use BCTypes
  use ADjointPETSc
  use blockPointers
  use inputTimeSpectral
  use costFUnctions
  implicit none

  deallocate(costFuncMat, costFuncMatd)

  ! Destroy the empty vectors:
  call vecDestroy(fVec1,PETScIerr)
  call vecDestroy(fVec2,PETScIerr)

  call vecDestroy(w_like1,PETScIerr)
  call vecDestroy(w_like2,PETScIerr)

  call vecDestroy(psi_like1,PETScIerr)
  call vecDestroy(psi_like2,PETScIerr)

end subroutine releaseMemADjoint
