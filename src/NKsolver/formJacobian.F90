subroutine FormJacobian(snes,wVec,dRdw,dRdwPre,flag,ctx,ierr)
  use communication
  use precision 
  use NKSolverVars,only : ksp_solver_type,ksp_subspace,global_pc_type,&
       asm_overlap,local_pc_ilu_level,local_pc_ordering,NKfinitedifferencepc
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  SNES           snes
  Mat            dRdw,dRdwPre 
  KSP            ksp,subksp
  PC             pc,subpc
  Vec            wVec,rVec
  PetscFortranAddr ctx(3)
  MatStructure   flag 
  PetscInt nlocal,first,Nsub,length
  integer(kind=intType) ::ierr

  ! Local Variables
  logical secondHalo
  logical :: useAD,usePC,useTranspose
 
  ! Dummy assembly begin/end calls for the matrix-free Matrx
  call MatAssemblyBegin(dRdw,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call MatAssemblyEnd(dRdw,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Assemble the approximate PC

  if (NKFiniteDifferencePC) then
     useAD = .False.
     usePC = .True.
     useTranspose = .False.
     call setupStateResidualMatrix(dRdwPre,useAD,usePC,useTranspose)
  else
     call setupNK_PC(dRdwPre)
  end if


  flag = SAME_NONZERO_PATTERN
  ! Setup the required options for the KSP solver
  call SNESGetKSP(snes,ksp,ierr);                   call EChk(ierr,__FILE__,__LINE__)
  call KSPSetType(ksp,ksp_solver_type,ierr);       call EChk(ierr,__FILE__,__LINE__)
  call KSPGMRESSetRestart(ksp, ksp_subspace,ierr);  call EChk(ierr,__FILE__,__LINE__)
  call KSPSetPreconditionerSide(ksp,PC_RIGHT,ierr); call EChk(ierr,__FILE__,__LINE__)

  ! Setup the required options for the Global PC
  call KSPGetPC(ksp,pc,ierr);                 call EChk(ierr,__FILE__,__LINE__)
!   call PCSetType(pc,'hypre',ierr);            call EChk(ierr,__FILE__,__LINE__)
!   call PCHYPRESetType(pc,'euclid',ierr);   call EChk(ierr,__FILE__,__LINE__)
!   call PCFactorSetMatOrderingtype(pc, local_pc_ordering, ierr ); call EChk(ierr,__FILE__,__LINE__) 


  call PCSetType(pc,global_pc_type,ierr);     call EChk(ierr,__FILE__,__LINE__)

  if (trim(global_pc_type) == 'asm') then
     call PCASMSetOverlap(pc,asm_overlap,ierr);  call EChk(ierr,__FILE__,__LINE__)
     call PCSetup(pc,ierr);                      call EChk(ierr,__FILE__,__LINE__)
     call PCASMGetSubKSP( pc, nlocal,  first, subksp, ierr );          call EChk(ierr,__FILE__,__LINE__)  
  end if

  if (trim(global_pc_type) == 'bjacobi') then
     call PCSetup(pc,ierr);                      call EChk(ierr,__FILE__,__LINE__)
     call PCBJacobiGetSubKSP(pc,nlocal,first,subksp,ierr);   call EChk(ierr,__FILE__,__LINE__)
  end if


  ! Setup the required options for the Local PC
  call KSPGetPC(subksp, subpc, ierr );                              call EChk(ierr,__FILE__,__LINE__)
  call PCSetType(subpc, 'ilu', ierr);                       call EChk(ierr,__FILE__,__LINE__)
  call PCFactorSetLevels(subpc, local_pc_ilu_level, ierr);          call EChk(ierr,__FILE__,__LINE__)  
  call PCFactorSetMatOrderingtype(subpc, local_pc_ordering, ierr ); call EChk(ierr,__FILE__,__LINE__) 
  call KSPSetType(subksp, KSPPREONLY, ierr);    
  call EChk(ierr,__FILE__,__LINE__)  

  ierr = 0
  
end subroutine FormJacobian
