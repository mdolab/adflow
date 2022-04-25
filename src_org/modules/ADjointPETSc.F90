module ADjointPETSc

  !      This module contains the objects used by PETSc for the
  !      solution of the discrete adjoint equations.
  !
  use constants
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none

  Mat     dRdWT, dRdWPreT

  ! These are empty vectors
  Vec     w_like1, w_like2, psi_like1, psi_like2, psi_like3, x_like
  PetscErrorCode PETScIerr
  PetscFortranAddr   matfreectx(1)

  !adjointKSP   Linear solver (Krylov subspace method) context
  KSP     adjointKSP


  ! Initial, start and final adjoint residuals
  real(kind=alwaysRealType) :: adjResInit
  real(kind=alwaysRealType) :: adjResStart
  real(kind=alwaysRealType) :: adjResFinal
  logical :: adjointPETScVarsAllocated

end module ADjointPETSc
