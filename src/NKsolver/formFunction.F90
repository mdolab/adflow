subroutine FormFunction_mf(ctx,wVec,rVec,ierr)
#ifndef USE_NO_PETSC
  ! This is basically a copy of FormFunction, however it has a
  ! different calling sequence from PETSc. It performs the identical
  ! function. This is used for linear solve application for the
  ! aerostructural system pre-conditioner

  use precision
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

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
