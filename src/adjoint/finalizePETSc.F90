!
!     ******************************************************************
!     *                                                                *
!     * File:          finalizePETSc.F90                               *
!     * Author:        Andre C. Marta                                  *
!     * Starting date: 08-24-2005                                      *
!     * Last modified: 03-01-2007                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine finalizePETSc
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Finalize PETSc by calling the appropriate routine              *
  !     * PetscFinalize provided in the PETSc library. This              *
  !     * automatically calls MPI_Finalize().                            *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc, only : PETScIerr
  implicit none
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !
#ifndef USE_NO_PETSC

  call PetscFinalize(PETScIerr)

#endif
end subroutine finalizePETSc
