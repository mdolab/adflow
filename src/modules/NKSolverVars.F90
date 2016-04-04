module NKsolverVars
 
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
  ! PETSc Variables

  ! PETSc Matrices:
  ! dRdw: This is the actual matrix-free matrix computed with FD
  ! dRdwPre: The preconditoner matrix for NK method. This matrix is stored.
  ! dRdwPseudo: Shell matrix used with the pseudo-transient 
  !             continuation method. 

  Mat  dRdw, dRdwPre, dRdwPseudo

  ! PETSc Vectors:
  ! wVec: PETsc version of SUmb 'w'
  ! rVec: PETSc version of SUmb 'dw', but divided by volume
  ! deltaW: Update to the wVec from linear solution
  ! diagV: Diagonal lumping term

  Vec wVec, rVec, deltaW, work, g, wBase

  ! NK_KSP: The ksp object for solving the newton udpate
  KSP  NK_KSP

  PetscFortranAddr   ctx(1)
#endif

  ! Options for NK Slver
  logical :: useNKSolver
  integer(kind=intType) :: NK_jacobianLag
  integer(kind=intType) :: NK_subspace
  integer(kind=intType) :: NK_asmOverlap
  integer(kind=intType) :: NK_iluFill
  integer(kind=intType) :: NK_innerPreConIts
  integer(kind=intType) :: NK_outerPreConIts
  integer(kind=intType) :: NK_LS
  logical :: NK_useEW 
  logical :: NK_ADPC
  logical :: NK_viscPC
  real(kind=realType) :: NK_CFL0
  real(kind=realType) :: NK_switchTol
  real(kind=realType) :: NK_rtolInit
  real(kind=realType) :: NK_divTol = 10

  ! Misc variables
  logical :: NK_solverSetup
  integer(kind=intType) :: NK_iter

  ! Eisenstat-Walker Parameters
  integer(kind=intType) :: ew_version
  real(kind=realType) :: ew_rtol_0
  real(kind=realType) :: ew_rtol_max
  real(kind=realType) :: ew_gamma
  real(kind=realType) :: ew_alpha
  real(kind=realType) :: ew_alpha2
  real(kind=realType) :: ew_threshold
  real(kind=realType) :: rtolLast, oldNorm
  
  ! Misc Parameters
  real(kind=realType) :: totalR0, totalRStart, totalRFinal, totalR
  real(kind=realType) :: rhoRes0, rhoResStart, rhoResFinal, rhoRes
  logical :: freeStreamResSet
  real(kind=realType) :: NK_CFL

  ! Variables for non-monotone line search
  real(kind=realType), dimension(:), allocatable :: func_evals
  integer(kind=intType) :: Mmax
  integer(kind=intType) :: iter_k
  integer(kind=intType) :: iter_m

  ! Line search parameters
  integer(kind=intType), parameter :: noLineSearch = 0_intType, &
                                      cubicLineSearch = 1_intType, &
                                      nonMonotoneLineSearch = 2_intType

  ! Parameter for external preconditioner
  integer(kind=intType) :: applyPCSubSpaceSize

end module NKsolverVars
