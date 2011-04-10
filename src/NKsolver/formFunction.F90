subroutine FormFunction(snes,wVec,rVec,ctx,ierr)
  ! ---------------------------------------------------------------------
  !
  !  FormFunction - Evaluates nonlinear function, f(x).
  !
  !  Input Parameters:
  !  snes  - the SNES context
  !  wVec  - input vector
  !
  !  Currnet ctx( ) doesn't have anythign in it
  !
  !  Output Parameter:
  !  f     - vector with newly computed function
  use precision
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! PETSc Variables
  SNES    snes
  Vec     wVec, rVec
  PetscFortranAddr ctx(*)
  integer(kind=intType) :: ierr

  ! This is just a shell routine that runs the more broadly useful
  ! computeResidualNK subroutine
  call setW(wVec)
  call computeResidualNK()
  call setRVec(rVec)

  ! We don't check an error here, so just pass back zero
  ierr = 0

end subroutine FormFunction


subroutine FormFunction2(ctx,wVec,rVec,ierr)
  ! This is basically a copy of FormFunction, however it has a
  ! different calling sequence from PETSc. It performs the identical
  ! function. This is used for linear solve application for the
  ! aerostructural system pre-conditioner

  use precision
  use flowVarRefState
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

end subroutine FormFunction2
