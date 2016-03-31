!
!     ******************************************************************
!     *                                                                *
!     * File:          NKsolverVars.F90                                *
!     * Author:        Gaetan K.W. Kenway                              *
!     * Starting date: 12-01-2010                                      *
!     * Last modified: 12-01-2010                                      *
!     *                                                                *
!     ******************************************************************
!
module NKsolverVars
  !
  !     ******************************************************************
  !     *                                                                *
  !     * This module contains the variables needed by the Newton-Krylov *
  !     * solver.                                                        *
  !     *                                                                *
  !     ******************************************************************
  !
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

  ! PETSc SNES: 
  ! snes: Non-linear solver solution context for NK problem. Currently All Unsued.

  SNES  snes, psnes, outer_snes

  ! PETSc Matrices:
  ! dRdw: This is a matrix-free matrix for the state residal matrix. 
  !       Matrix-vector products are computed with finite differences
  ! dRdwPre: The preconditoner matrix for NK method. This matrix is stored
  ! dRdwPseudo: Shell matrix used with the pseudo-transient 
  !             continuation method. Currently unused.

  Mat  dRdw, dRdwPre, dRdwPseudo

  ! PETSc Vectors:
  ! wVec: PETsc version of SUmb 'w'
  ! rVec: PETSc versin of SUmb 'dw', but divided by volume
  ! deltaW: Update to the wVec from linear solution
  ! diagV: Diagonal lumping term

  Vec wVec, rVec, rvec2, deltaW, diagV, work, g, scaleVec, wBase, rBase, diag, rhs
  Vec w_like1, w_like2
  ! PETSc KSP/PC 
  ! NK_KSP: The ksp object for solving the newton udpate

  KSP  newtonKrylovKSP

  ! PETSc Misc:
  SNESConvergedReason reason
  PetscFortranAddr   ctx(3), ctx2(3)
#endif

  ! Non-linear Solver Options
  integer(kind=intType) :: jacobian_lag
  integer(kind=intType) :: NKSolveCount
  logical :: NKuseEW
  logical :: useNKSolver
  logical :: NKSolverSetup
  logical :: NKPCSetup
  logical :: NKADPC
  logical :: NKViscPC
  logical :: RKreset

  real(kind=realType) :: resSum(8)

  ! Non-linear Solver Tolerances
  real(kind=realType) :: snes_atol 
  real(kind=realType) :: snes_rtol
  real(kind=realType) :: snes_stol
  integer(kind=intType) :: snes_max_its
  integer(kind=intType) :: snes_max_funcs

  ! Eisenstat-Walker Paramtersr
  integer(kind=intType) :: ew_version
  real(kind=realType) :: ew_rtol_0
  real(kind=realType) :: ew_rtol_max
  real(kind=realType) :: ew_gamma
  real(kind=realType) :: ew_alpha
  real(kind=realType) :: ew_alpha2
  real(kind=realType) :: ew_threshold
  
  ! Krylov-Solver Options
  character(maxStringLen) :: ksp_solver_type
  integer(kind=intType) :: ksp_subspace
  integer(kind=intType) :: asm_overlap
  character(maxStringLen):: global_pc_type
  character(maxStringLen):: global_pc_side
  character(maxStringLen) :: local_pc_ordering
  integer(kind=intType)   :: local_pc_ilu_level
  integer(kind=intType) :: innerPreConIts
  integer(kind=intType) :: outerPreConIts

  ! Krylov-Solver Tolerances
  real(kind=realType) :: ksp_atol
  real(kind=realType) :: ksp_rtol
  real(kind=realType) :: ksp_rtol_init
  integer(kind=intTYpe) :: ksp_max_it
  real(kind=realType) :: ksp_div_tol = 10

  ! Parameter from switching from RK to NK
  real(kind=realType) :: NK_switch_tol

  ! Misc Parameters
  real(kind=realType) :: totalR0, totalRStart, totalRFinal
  real(kind=realType) :: rhoRes0, rhoResStart, rhoResFinal
  real(kind=realType) :: rhoResL1Start
  logical :: freeStreamResSet
  real(kind=realType) :: CFL0
  integer(kind=intType) :: iterTot0
  integer(kind=intType) :: applyPCSubSpaceSize

  real(kind=realType), dimension(:), allocatable :: func_evals
  integer(kind=intType) :: Mmax
  integer(kind=intType) :: iter_k
  integer(kind=intType) :: iter_m

  ! Line search parameters
  integer(kind=intType), parameter :: noLineSearch = 0_intType, &
                                      cubicLineSearch = 1_intType, &
                                      nonMonotoneLineSearch = 2_intType
  integer(kind=intType) :: NKLS

end module NKsolverVars
