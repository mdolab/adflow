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
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! PETSc Variables
  SNES               snes 
  PetscFortranAddr   ctx(3)
  Mat                dRdw,dRdwPre 
  Vec wVec,rVec
  KSP                global_ksp,local_ksp
  PC                 global_pc ,local_pc
  SNESConvergedReason reason


  ! Non-linear Solver Options
  integer(kind=intType) :: jacobian_lag
  logical :: useEW
  logical :: useNKSolver
  logical :: NKSolverSetup
  logical :: NKSolvedOnce
  logical :: NKPCSetup
  logical :: NKFiniteDifferencePC
  logical :: RKreset
  integer(kind=intType) :: nRKreset
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
  character(maxStringLen) :: local_pc_type
  character(maxStringLen) :: local_pc_ordering
  integer(kind=intType)   :: local_pc_ilu_level

  ! Krylov-Solver Tolerances
  real(kind=realType) :: ksp_atol
  real(kind=realType) :: ksp_rtol
  integer(kind=intTYpe) :: ksp_max_it
  real(kind=realType) :: ksp_div_tol = 10

  ! Parameter from switching from RK to NK
  real(kind=realType) :: NK_switch_tol

  ! Define some of the named constants for PETSc 

  ! KSP Types
  character, parameter :: ksp_gmres = "gmres"
  character, parameter :: ksP_fgmres = "fgmres"
  character, parameter :: ksp_bicgstab = "bcgs"

  ! Global PC Types
  character, parameter :: pc_blockjacobi = "bjacobi"
  character, parameter :: pc_jacobi = "jacobi"
  character, parameter :: pc_asm    = "asm"

  ! Local PC Types
  character, parameter :: pc_ilu = "ilu"
  character, parameter :: pc_lu  = "lu"
  
  ! Local Orderings
  character, parameter :: ord_natural = "natural"
  character, parameter :: ord_rcm     = "rcm"
  character, parameter :: ord_nd      = "nd"
  character, parameter :: ord_owd     = "owd"

  ! PC Side
  character, parameter :: side_left = "left"
  character, parameter :: side_right = "right"
  
  ! Misc Parameters
  real(kind=realType) :: totalR0,totalRStart,totalRFinal
  real(kind=realType) :: rhoRes0,rhoResStart,rhoResFinal
  real(kind=realType) :: CFL0
  integer(kind=intType) :: iterTot0

end module NKsolverVars
