module ANKsolverVars

  use constants
  implicit none


#ifndef USE_NO_PETSC

#define PETSC_AVOID_MPIF_H
#include "include/petscversion.h"

#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

  Mat  dRdwPre
  Vec wVec, rVec, deltaW
  KSP  ANK_KSP
  PetscFortranAddr   ctx(1)
#endif

  ! Options for ANK Solver
  logical :: useANKSolver
  integer(kind=intType) :: ANK_jacobianLag
  integer(kind=intType) :: ANK_subSpace
  integer(kind=intType) :: ANK_asmOverlap
  integer(kind=intType) :: ANK_iluFill
  integer(kind=intType) :: ANK_innerPreConIts
  real(kind=realType)   :: ANK_rtol
  real(kind=realType)   :: ANK_switchTol
  real(kind=realType)   :: ANK_divTol = 10

  ! Misc variables
  real(kind=realType) :: ANK_CFL
  logical :: ANK_solverSetup
  integer(kind=intTYpe) :: ANK_iter

end module ANKsolverVars
