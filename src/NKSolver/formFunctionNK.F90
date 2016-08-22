subroutine FormFunction_mf(ctx, wVec, rVec, ierr)
#ifndef USE_NO_PETSC
  ! This is basically a copy of FormFunction, however it has a
  ! different calling sequence from PETSc. It performs the identical
  ! function. This is used for linear solve application for the
  ! aerostructural system pre-conditioner

  use constants
  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

  ! PETSc Variables
  PetscFortranAddr ctx(*)
  Vec     wVec, rVec
  integer(kind=intType) :: ierr
 
  ! This is just a shell routine that runs the more broadly useful
  ! computeResidualNK subroutine
 
  call setW(wVec)
  call computeResidualNK()
  call setRVec(rVec)
  ! We don't check an error here, so just pass back zero
  ierr = 0
#endif
end subroutine FormFunction_mf
