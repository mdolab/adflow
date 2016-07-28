!
!     ******************************************************************
!     *                                                                *
!     * File:          ADjointPETSc.F90                                *
!     * Author:        Andre C. Marta                                  *
!     * Starting date: 12-15-2005                                      *
!     * Last modified: 06-14-2007                                      *
!     *                                                                *
!     ******************************************************************
!
module ADjointPETSc

  !     ******************************************************************
  !     *                                                                *
  !     * This module contains the objects used by PETSc for the         *
  !     * solution of the discrete adjoint equations.                    *
  !     *                                                                *
  !     ******************************************************************
  !
  use constants
  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

  Mat     dRdWT, dRdWPreT

  ! These are empty vectors
  Vec     w_like1, w_like2, psi_like1, psi_like2, psi_like3, x_like
  PetscErrorCode PETScIerr
  PetscFortranAddr   matfreectx(1)

  !adjointKSP   Linear solver (Krylov subspace method) context
  KSP     adjointKSP


  ! Initial, start and final adjoint residuals
  real(kind=realType) :: adjResInit
  real(kind=realType) :: adjResStart
  real(kind=realType) :: adjResFinal
  logical :: adjointPETScVarsAllocated

end module ADjointPETSc
