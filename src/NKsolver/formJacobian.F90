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
