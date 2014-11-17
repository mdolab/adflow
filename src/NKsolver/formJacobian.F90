subroutine FormJacobian()
#ifndef USE_NO_PETSC

  use NKSolverVars
  use inputADjoint, only : viscPC
  implicit none
  ! Local Variables
  character(len=maxStringLen) :: preConSide, localPCType
  integer(kind=intType) ::ierr
  logical :: useAD, usePC, useTranspose, useObjective, tmp

  ! Dummy assembly begin/end calls for the matrix-free Matrx
  call MatAssemblyBegin(dRdw, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call MatAssemblyEnd(dRdw, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Assemble the approximate PC (fine leve, level 1)
  useAD = NKADPC
  usePC = .True.
  useTranspose = .False.
  useObjective = .False.
  tmp = viscPC ! Save what is in viscPC and set to the NKvarible
  viscPC = NKViscPC
  call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, &
       useObjective, 1_intType)
  ! Reset saved value
  viscPC = tmp

  ! Setup KSP Options
  preConSide = 'right'
  localPCType = 'ilu'

  ! Setup the KSP using the same code as used for the adjoint
  call setupStandardKSP(newtonKrylovKSP, ksp_solver_type, ksp_subspace, &
       preConSide, global_pc_type, asm_overlap, outerPreConIts, localPCType, &
       local_pc_ordering, local_pc_ilu_level, innerPreConIts)

  ! Don't do iterative refinement for the NKSolver.
  call KSPGMRESSetCGSRefinementType(newtonKrylovKSP, &
       KSP_GMRES_CGS_REFINE_NEVER, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
#endif

end subroutine FormJacobian
