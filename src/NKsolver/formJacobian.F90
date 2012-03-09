subroutine FormJacobian()
#ifndef USE_NO_PETSC

  use NKSolverVars

  ! Local Variables

  integer(kind=intType) ::ierr
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

  ! Setup KSP Options
  call setupNKKSP()
#endif
end subroutine FormJacobian
  
subroutine setupNKKSP()

#ifndef USE_NO_PETSC

  use communication
  use iteration
  use NKSolverVars

  implicit none

  ! Local Variables
  integer(kind=intType) ::ierr, nlocal, first
 
  ! ----------------------------------------------
  ! Setup the required options for the KSP object
  ! ----------------------------------------------

  ! Set Solver Type
  call KSPSetType(global_ksp,ksp_solver_type,ierr);
  call EChk(ierr,__FILE__,__LINE__)

  ! Set Subspace Size
  call KSPGMRESSetRestart(global_ksp, ksp_subspace,ierr); 
  call EChk(ierr,__FILE__,__LINE__)

  ! Set PC Side as RIGHT only
  call KSPSetPCSide(global_ksp,PC_RIGHT,ierr);
  call EChk(ierr,__FILE__,__LINE__)

  ! Get the PC Handle to make modifications:
  call KSPGetPC(global_ksp,global_pc,ierr);             
  call EChk(ierr,__FILE__,__LINE__)

  call PCSetType(global_pc,global_pc_type,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (trim(global_pc_type) == 'asm') then
     call PCASMSetOverlap(global_pc,asm_overlap,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call PCSetup(global_pc,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call PCASMGetSubKSP(global_pc, nlocal, first, local_ksp, ierr )
     call EChk(ierr,__FILE__,__LINE__)  
  end if

  if (trim(global_pc_type) == 'bjacobi') then
     call PCSetup(global_pc,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call PCBJacobiGetSubKSP(global_pc,nlocal,first,local_ksp,ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end if

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
#endif
end subroutine SetupNKKSP
