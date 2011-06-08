subroutine FormJacobian_custom()
  use communication
  use precision 
  use NKSolverVars, only: dRdw,dRdwPre,NKFiniteDifferencePC,ksp,&
       global_pc_type,asm_overlap,local_pc_ordering,local_pc_ilu_level,&
       diagV

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

  if (NKFiniteDifferencePC) then
     useAD = .False.
     usePC = .True.
     useTranspose = .False.
     call setupStateResidualMatrix(dRdwPre,useAD,usePC,useTranspose)

  else
     call setupNK_PC(dRdwPre)
  end if

  ! Add diagV to the diagonal
  !call MatDiagonalSet(dRdwPre,diagV,ADD_VALUES,ierr)
  !call EChk(ierr,__FILE__,__LINE__)
  
!   ! Redo Assembly
!   call MatAssemblyBegin(dRdwPre,MAT_FINAL_ASSEMBLY,ierr)
!   call EChk(ierr,__FILE__,__LINE__)
!   call MatAssemblyEnd  (dRdwPre,MAT_FINAL_ASSEMBLY,ierr)
!   call EChk(ierr,__FILE__,__LINE__)


  ! Setup the required options for the Global PC
  call KSPGetPC(ksp,pc,ierr);             
  call EChk(ierr,__FILE__,__LINE__)

  call PCSetType(pc,global_pc_type,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (trim(global_pc_type) == 'asm') then
     call PCASMSetOverlap(pc,asm_overlap,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call PCSetup(pc,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call PCASMGetSubKSP( pc, nlocal,  first, subksp, ierr )
     call EChk(ierr,__FILE__,__LINE__)  
  end if

  if (trim(global_pc_type) == 'bjacobi') then
     call PCSetup(pc,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call PCBJacobiGetSubKSP(pc,nlocal,first,subksp,ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end if

  ! Setup the required options for the Local PC
  call KSPGetPC(subksp, subpc, ierr )
  call EChk(ierr,__FILE__,__LINE__)

  call PCSetType(subpc, 'ilu', ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call PCFactorSetLevels(subpc, local_pc_ilu_level, ierr)
  call EChk(ierr,__FILE__,__LINE__)  
  call PCFactorSetMatOrderingtype(subpc, local_pc_ordering, ierr )
  call EChk(ierr,__FILE__,__LINE__) 
  call KSPSetType(subksp, KSPPREONLY, ierr)
  call EChk(ierr,__FILE__,__LINE__)  

  ierr = 0
  
end subroutine FormJacobian_custom

