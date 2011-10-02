subroutine FormFunction_snes(snes,wVec,rVec,ctx,ierr)
#ifndef USE_NO_PETSC
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
  ! computeResidualNK subroutin

  call setW(wVec)
    call computeResidualNK()
  call setRVec(rVec)

  ! We don't check an error here, so just pass back zero
  ierr = 0
#endif
end subroutine FormFunction_snes

subroutine FormFunction_mf(ctx,wVec,rVec,ierr)
#ifndef USE_NO_PETSC
  ! This is basically a copy of FormFunction, however it has a
  ! different calling sequence from PETSc. It performs the identical
  ! function. This is used for linear solve application for the
  ! aerostructural system pre-conditioner

  use communication
  use precision
  use flowVarRefState
  use inputtimespectral
  use blockPointers
  use nksolvervars, only : times
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! PETSc Variables
  PetscFortranAddr ctx(*)
  Vec     wVec, rVec
  integer(kind=intType) :: ierr,sps,nn
  real(kind=realType) :: ovr_CFL
  ! This is just a shell routine that runs the more broadly useful
  ! computeResidualNK subroutine
 
  ! Also try doing the built-in vec scats

  call setW(wVec)
  call computeResidualNK()
  call setRVec(rVec)
  ! We don't check an error here, so just pass back zero
  ierr = 0
#endif
end subroutine FormFunction_mf

subroutine FormFunction_ts(pts,t,wVec,rVec,ctx,ierr)
#ifndef USE_NO_PETSC
  ! This is basically a copy of FormFunction, however it has a
  ! different calling sequence from PETSc. It performs the identical
  ! function. This is used for Pseudo time stepping. 

  use communication
  use precision
  use flowVarRefState
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
#include "include/finclude/petscts.h"

  ! PETSc Variables
  TS pts
  real(kind=realType) :: t
  PetscFortranAddr ctx(*)
  Vec     wVec, rVec
  integer(kind=intType) :: ierr
  
  ! This is just a shell routine that runs the more broadly useful
  ! computeResidualNK subroutine
  if (myid == 0) then
     print *,'Form Func 3'
  end if
  call setW(wVec)
  call computeResidualNK()
  ! We don't check an error here, so just pass back zero
  ierr = 0
#endif
end subroutine FormFunction_ts
