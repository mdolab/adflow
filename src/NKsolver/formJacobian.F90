subroutine FormJacobian()
#ifndef USE_NO_PETSC

  use NKSolverVars

  ! Local Variables
  character(len=maxStringLen) :: preConSide, localPCType
  integer(kind=intType) ::ierr
  logical :: useAD, usePC, useTranspose

  ! Dummy assembly begin/end calls for the matrix-free Matrx
  call MatAssemblyBegin(dRdw, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call MatAssemblyEnd(dRdw, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Assemble the approximate PC (fine leve, level 1)
  useAD = .False.
  usePC = .True.
  useTranspose = .False.
  call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, 1_intType)

  ! Setup KSP Options
  preConSide = 'right'
  localPCType = 'ilu'
 
  ! There is an issue currently with this: When outerPreConIts and
  ! innerPreConIts both are 1, it produces the same sequence of
  ! iterates as the orignal, (as is expected) but it takes 50%
  ! longer. Not sure what is causing this. Leave the old on in for
  ! now.
!   call setupStandardKSP(newtonKrylovKSP, ksp_solver_type, ksp_subspace, preConSide, &
!        global_pc_type, asm_overlap, outerPreConIts, localPCType, &
!        local_pc_ordering, local_pc_ilu_level, innerPreConIts)

  call setupNKKSP()
#endif
end subroutine FormJacobian


subroutine setupNKKSP()

#ifndef USE_NO_PETSC

  use communication
  use iteration
  use NKSolverVars, only : newtonKrylovKSP, ksp_solver_type, &
       asm_overlap, local_pc_ilu_level, local_pc_ordering, &
       ksp_subspace, global_pc_type

  implicit none

#define PETSC_AVOID_MPIF_H
#include "finclude/petsc.h"


  ! Local Variables
  KSP local_ksp
  PC  local_PC, global_PC
  integer(kind=intType) ::ierr, nlocal, first

  ! ----------------------------------------------
  ! Setup the required options for the KSP object
  ! ----------------------------------------------

  ! Set Solver Type
  call KSPSetType(newtonKrylovKSP,ksp_solver_type,ierr);
  call EChk(ierr,__FILE__,__LINE__)

  ! Set Subspace Size
  call KSPGMRESSetRestart(newtonKrylovKSP, ksp_subspace,ierr); 
  call EChk(ierr,__FILE__,__LINE__)

  ! Set PC Side as RIGHT only
  call KSPSetPCSide(newtonKrylovKSP,PC_RIGHT,ierr);
  call EChk(ierr,__FILE__,__LINE__)

  ! Get the PC Handle to make modifications:
  call KSPGetPC(newtonKrylovKSP,global_pc,ierr);             
  call EChk(ierr,__FILE__,__LINE__)


  if (trim(global_pc_type) == 'asm') then
     call PCSetType(global_pc,global_pc_type,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call PCASMSetOverlap(global_pc,asm_overlap,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call PCSetup(global_pc,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call PCASMGetSubKSP(global_pc, nlocal, first, local_ksp, ierr )
     call EChk(ierr,__FILE__,__LINE__)  
  end if

  if (trim(global_pc_type) == 'bjacobi') then
     call PCSetType(global_pc,global_pc_type,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call PCSetup(global_pc,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call PCBJacobiGetSubKSP(global_pc,nlocal,first,local_ksp,ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end if

  if (trim(global_pc_type) == 'mg') then
!      call PCSetType(global_pc,'shell',ierr)
!      call EChk(ierr,__FILE__,__LINE__)

!      ! Set function for doing preconditioner application
!      !call PCShellSetApply(global_pc,myShellPCApply,ierr)
!      !call PCShellSetApplyTranspose(global_pc,myShellPCApply,ierr)
!      call PCShellSetApply(global_pc,myShellPCApply,ierr)
!      call EChk(ierr,__FILE__,__LINE__)
  end if

  ! Setup local KSP objects for asm and bjacobi conditioners
  if (trim(global_pc_type) == 'asm' .or. trim(global_pc_type) == 'bjacobi') then

     ! Setup the required options for the Local PC
     call KSPGetPC(local_ksp, local_pc, ierr )
     call EChk(ierr,__FILE__,__LINE__)

     call PCSetType(local_pc, 'ilu', ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call PCFactorSetLevels(local_pc, local_pc_ilu_level, ierr)
     call EChk(ierr,__FILE__,__LINE__)  
     call PCFactorSetMatOrderingtype(local_pc, local_pc_ordering, ierr )
     call EChk(ierr,__FILE__,__LINE__) 
     call KSPSetType(local_ksp, 'preonly', ierr)
     call EChk(ierr,__FILE__,__LINE__)  
  end if

#endif
   end subroutine SetupNKKSP
