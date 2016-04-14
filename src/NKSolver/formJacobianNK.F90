subroutine FormJacobianNK
#ifndef USE_NO_PETSC

  use NKSolverVars, only : work, dRdw, dRdWPre, NK_ADPC, NK_asmOverlap,&
       NK_innerPreConIts, NK_outerPreconIts, MAT_FINAL_ASSEMBLY, NK_viscPC, NK_iluFill, &
       NK_subspace, KSP_GMRES_CGS_REFINE_NEVER, NK_KSP, NK_CFL
  use inputADjoint, only : viscPC
  use inputiteration
  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate

  implicit none
#define PETSC_AVOID_MPIF_H


#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#endif

  ! Local Variables
  character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
  integer(kind=intType) :: ierr
  logical :: useAD, usePC, useTranspose, useObjective, tmp
  integer(kind=intType) :: i, j, k, l, ii, nn, sps
  real(kind=realType) :: dt, ovv
  real(kind=realType), pointer :: diag(:)
 interface
     subroutine setupStateResidualMatrix(matrix, useAD, usePC, useTranspose, &
          useObjective, frozenTurb, level, matrixTurb)
       use precision
       implicit none
#define PETSC_AVOID_MPIF_H
#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

       Mat :: matrix
       Mat, optional :: matrixTurb
       ! Input Variables
       logical, intent(in) :: useAD, usePC, useTranspose, useObjective, frozenTurb
       integer(kind=intType), intent(in) :: level
     end subroutine setupStateResidualMatrix
  end interface

  ! Dummy assembly begin/end calls for the matrix-free Matrx
  call MatAssemblyBegin(dRdw, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call MatAssemblyEnd(dRdw, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Assemble the approximate PC (fine leve, level 1)
  useAD = NK_ADPC
  usePC = .True.
  useTranspose = .False.
  useObjective = .False.
  tmp = viscPC ! Save what is in viscPC and set to the NKvarible
  viscPC = NK_viscPC

  call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, &
       useObjective, .False., 1_intType)
  ! Reset saved value
  viscPC = tmp

  call VecGetArrayF90(work, diag, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  diag(:) = one/NK_CFL

  call VecRestoreArrayF90(work, diag, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MatDiagonalSet(dRdwPre, work, ADD_VALUES, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Setup KSP Options
  preConSide = 'right'
  localPCType = 'ilu'
  kspObjectType = 'gmres'
  globalPCType = 'asm'
  localOrdering = 'rcm'

  ! Setup the KSP using the same code as used for the adjoint
  call setupStandardKSP(NK_KSP, kspObjectType, NK_subSpace, &
       preConSide, globalPCType, NK_asmOverlap, NK_outerPreConIts, localPCType, &
       localOrdering, NK_iluFill, NK_innerPreConIts)

  ! Don't do iterative refinement for the NKSolver.
  call KSPGMRESSetCGSRefinementType(NK_KSP, &
       KSP_GMRES_CGS_REFINE_NEVER, ierr)
  call EChk(ierr, __FILE__, __LINE__)

#endif

end subroutine FormJacobianNK
