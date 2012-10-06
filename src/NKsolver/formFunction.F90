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

subroutine FormFunction_snes(snes, wVec, rVec, ctx, ierr)
#ifndef USE_NO_PETSC
  ! This is basically a copy of FormFunction, however it has a
  ! different calling sequence from PETSc. It performs the identical
  ! function. This is used for linear solve application for the
  ! aerostructural system pre-conditioner

  use precision
  use communication
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! PETSc Variables
  SNES snes
  PetscFortranAddr ctx(*)
  Vec     wVec, rVec
  integer(kind=intType) :: ierr
  real(kind=realType) :: val
  ! This is just a shell routine that runs the more broadly useful
  ! computeResidualNK subroutine
  call setW2(wVec)
  call computeResidualNK()
  call convergenceInfo
  call setRVec2(rVec)
  call vecnorm(rVec, NORM_2, val, ierr)


  ! We don't check an error here, so just pass back zero
  ierr = 0
#endif
end subroutine FormFunction_snes

subroutine psnes_func(psnes, x, ierr)
#ifndef USE_NO_PETSC
  use communication
  use precision
  use flowvarrefstate

  implicit none
#define PETSC_AVOID_MPIF_H
#include "finclude/petsc.h"

  ! Input/Output
  SNES psnes
  Vec x, rvec
  PetscErrorCode ierr

  integer(kind=intType) :: i
  real(kind=realType) :: val
  ! Basically we have 'x' which is sumb's w. We just solve the sucker using RK

  call setW2(x)
  
  do i=1,3
     call executeMGCycle
   !  call convergenceInfo
     end do

  ! Set our current soultion back into x
  call setWVec2(x)

  call SNESGetFunction(psnes, rVec, PETSC_NULL_FUNCTION, PETSC_NULL_OBJECT, ierr)

  call setRvec2(rVec)
  call vecnorm(rVec, NORM_2, val, ierr)
  call SNESSetFunctionNorm(psnes, val, ierr)

  ierr = 0
#endif
end subroutine Psnes_func
