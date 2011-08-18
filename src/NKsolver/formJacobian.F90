subroutine FormJacobian()
  use communication
  use precision 
  use iteration
  use NKSolverVars, only: dRdw,dRdwPre,NKFiniteDifferencePC,ksp,&
       global_pc_type,asm_overlap,local_pc_ordering,local_pc_ilu_level
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  KSP subksp
  PC  subpc,pc
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

  useAD = .False.
  usePC = .True.
  useTranspose = .False.
  call setupStateResidualMatrix(dRdwPre,useAD,usePC,useTranspose)

  ! Setup the required options for the KSP object
  call NKSeutp_KSP(ksp)
  
end subroutine FormJacobian

! subroutine NKSetup_KSP(ksp)

!   ! Common routine for setting up the KSP solver for NK solver. This
!   ! can be called from multiple places.

!   use NKSolverVars, only: dRdw,dRdwPre,NKFiniteDifferencePC,ksp,&
!        global_pc_type,asm_overlap,local_pc_ordering,local_pc_ilu_level


!   implicit none
! #define PETSC_AVOID_MPIF_H
! #include "include/finclude/petsc.h"

!   ! Input variables
!   KSP ksp
!   integer(kind=intType) :: ierr

!   ! Set Solver Type
!   call KSPSetType(ksp,ksp_solver_type,ierr);
!   call EChk(ierr,__FILE__,__LINE__)

!   ! Set Subspace Size
!   call KSPGMRESSetRestart(ksp, ksp_subspace,ierr); 
!   call EChk(ierr,__FILE__,__LINE__)

!   ! Set PC Side as RIGHT only
!   call KSPSetPreconditionerSide(ksp,PC_RIGHT,ierr);
!   call EChk(ierr,__FILE__,__LINE__)

!   ! Get the PC Handle to make modifications:
!   call KSPGetPC(ksp,pc,ierr);             
!   call EChk(ierr,__FILE__,__LINE__)

!   call PCSetType(pc,global_pc_type,ierr)
!   call EChk(ierr,__FILE__,__LINE__)

!   if (trim(global_pc_type) == 'asm') then
!      call PCASMSetOverlap(pc,asm_overlap,ierr)
!      call EChk(ierr,__FILE__,__LINE__)
!      call PCSetup(pc,ierr)
!      call EChk(ierr,__FILE__,__LINE__)
!      call PCASMGetSubKSP( pc, nlocal,  first, subksp, ierr )
!      call EChk(ierr,__FILE__,__LINE__)  
!   end if

!   if (trim(global_pc_type) == 'bjacobi') then
!      call PCSetup(pc,ierr)
!      call EChk(ierr,__FILE__,__LINE__)
!      call PCBJacobiGetSubKSP(pc,nlocal,first,subksp,ierr)
!      call EChk(ierr,__FILE__,__LINE__)
!   end if

!   ! Setup the required options for the Local PC
!   call KSPGetPC(subksp, subpc, ierr )
!   call EChk(ierr,__FILE__,__LINE__)

!   call PCSetType(subpc, 'ilu', ierr)
!   call EChk(ierr,__FILE__,__LINE__)
!   call PCFactorSetLevels(subpc, local_pc_ilu_level, ierr)
!   call EChk(ierr,__FILE__,__LINE__)  
!   call PCFactorSetMatOrderingtype(subpc, local_pc_ordering, ierr )
!   call EChk(ierr,__FILE__,__LINE__) 
!   call KSPSetType(subksp, KSPPREONLY, ierr)
!   call EChk(ierr,__FILE__,__LINE__)  

! end subroutine NKSetup_KSP
