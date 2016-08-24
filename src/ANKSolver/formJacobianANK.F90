subroutine FormJacobianANK

  use constants
  use ANKSolverVars, only : dRdWPre, ANK_asmOverlap, &
       ANK_innerPreConIts,  MAT_FINAL_ASSEMBLY, ANK_iluFill, &
       ANK_subSpace, KSP_GMRES_CGS_REFINE_NEVER, ANK_KSP, ANK_CFL, deltaW, &
       ANK_KSPTurb, deltawTurb, dRdwPreTurb, ANK_useTurbDADI
  use flowVarRefState, only : nw, nwf
  use inputADjoint, only : viscPC
  use utils, only : EChk
  implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  ! Local Variables
  character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
  integer(kind=intType) ::ierr
  logical :: useAD, usePC, useTranspose, useObjective, tmp
  real(kind=realType) ::  dt
  integer(kind=intType) :: i, j, k, l, ii, nn, sps, outerPreConIts
  real(kind=realType), pointer :: diag(:)
  external :: myKSPMonitor
  interface
     subroutine setupStateResidualMatrix(matrix, useAD, usePC, useTranspose, &
          useObjective, frozenTurb, level, matrixTurb)
       use precision
       implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
       Mat :: matrix
       Mat, optional :: matrixTurb
       ! Input Variables
       logical, intent(in) :: useAD, usePC, useTranspose, useObjective, frozenTurb
       integer(kind=intType), intent(in) :: level
     end subroutine setupStateResidualMatrix
  end interface

  ! Assemble the approximate PC (fine leve, level 1)
  useAD = .False.
  usePC = .True.
  useTranspose = .False.
  useObjective = .False.
  tmp = viscPC ! Save what is in viscPC and set to the NKvarible
  viscPC = .False.
  if (ANK_useTurbDADI) then 
     call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, &
          useObjective, .True., 1_intType)
  else
     ! The turbulence jacobian will only be accuate with AD.
     useAD = .True.
     call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, &
          useObjective, .False., 1_intType, dRdwPreTurb)
  end if
  ! Reset saved value
  viscPC = tmp

  ! ----------- Setup Flow KSP ----------
  call VecGetArrayF90(deltaW, diag, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  diag(:) = one/ANK_CFL

  call VecRestoreArrayF90(deltaW, diag, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MatDiagonalSet(dRdwPre, deltaW, ADD_VALUES, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Setup KSP Options
  preConSide = 'right'
  localPCType = 'ilu'
  kspObjectType = 'gmres'
  globalPCType = 'asm'
  localOrdering = 'rcm'
  outerPreConIts = 1 
  ! Setup the KSP using the same code as used for the adjoint
  call setupStandardKSP(ANK_KSP, kspObjectType,  ANK_subSpace, &
       preConSide, globalPCType, ANK_asmOverlap, outerPreConIts, localPCType, &
       localOrdering, ANK_iluFill, ANK_innerPreConIts)

  ! Don't do iterative refinement for the NKSolver.
  call KSPGMRESSetCGSRefinementType(ANK_KSP, &
       KSP_GMRES_CGS_REFINE_NEVER, ierr)
  call EChk(ierr, __FILE__, __LINE__)

   ! ----------- Setup Turb KSP ----------
  if (.not. ANK_useTurbDADI .and. nw > nwf) then 
     call VecGetArrayF90(deltaWTurb, diag, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     diag(:) = one/ANK_CFL
     
     call VecRestoreArrayF90(deltaWTurb, diag, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call MatDiagonalSet(dRdwPreTurb, deltaWTurb, ADD_VALUES, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Setup KSP Options
     preConSide = 'right'
     localPCType = 'ilu'
     kspObjectType = 'gmres'
     globalPCType = 'asm'
     localOrdering = 'rcm'
     outerPreConIts = 1 
     ! Setup the KSP using the same code as used for the adjoint
     call setupStandardKSP(ANK_KSPTurb, kspObjectType,  ANK_SubSpace, &
          preConSide, globalPCType, ANK_asmOverlap, outerPreConIts, localPCType, &
          localOrdering, ANK_iluFill*2, ANK_innerPreConIts)

     ! Don't do iterative refinement for the NKSolver.
     call KSPGMRESSetCGSRefinementType(ANK_KSPTurb, &
          KSP_GMRES_CGS_REFINE_NEVER, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

end subroutine FormJacobianANK
