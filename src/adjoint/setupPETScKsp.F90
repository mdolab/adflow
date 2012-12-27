!
!     ******************************************************************
!     *                                                                *
!     * File:          setupPETScKsp.F90                               *
!     * Author:        Gaetan Kenway der                               *
!     * Starting date: 26-12-2012                                      *
!     * Last modified: 26-12-2012                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine setupPETScKsp

  use ADjointPETSc
  use ADjointVars
  use inputADjoint
  use communication
  implicit none

  !     Local variables.
  character(len=10)  :: pcType
  logical :: useAD, usePC, useTranspose 
  integer(kind=intType) :: ierr
  external MyKSPMonitor

#ifndef USE_NO_PETSC

  if (ApproxPC)then
     !setup the approximate PC Matrix
     useAD = .False.
     useTranspose = .True.
     usePC = .True.
     call setupStateResidualMatrix(drdwpret, useAD, usePC, useTranspose)

     !now set up KSP Context
     call KSPSetOperators(ksp, dRdWT, dRdWPreT, &
          DIFFERENT_NONZERO_PATTERN, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  else
     ! Use the exact jacobian.  Here the matrix that defines the
     ! linear system also serves as the preconditioning matrix.
     call KSPSetOperators(ksp, dRdWT, dRdWT, &
          SAME_NONZERO_PATTERN, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  ! First, KSPSetFromOptions MUST be called
  call KSPSetFromOptions(ksp, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the type of solver to use:
  call KSPSetType(ksp, ADJointSolverType, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! If we're using GMRES set the possible gmres restart
  call KSPGMRESSetRestart(ksp, adjRestart, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! If you're using GMRES, set refinement type
  call KSPGMRESSetCGSRefinementType(ksp, KSP_GMRES_CGS_REFINE_IFNEEDED, ierr)   
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the preconditioner side from option:
  if (trim(PCSide) == 'right') then
     call KSPSetPCSide(ksp, PC_RIGHT, ierr)
  else
     call KSPSetPCSide(ksp, PC_LEFT, ierr)
  end if
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the main tolerance for th main adjoint KSP solver
  call KSPSetTolerances(ksp, adjRelTol, adjAbsTol, adjDivTol, &
       adjMaxIter, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Extract preconditioning context for main KSP solver: (master_PC)
  call KSPGetPC(ksp, master_PC, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call KSPGetPC(ksp, pc, ierr)
  call EChk(ierr, __FILE__, __LINE__)

   ! Set the type of master_PC to ksp. This lets us do multiple
   ! iterations of preconditioner application
  call PCSetType(master_PC, 'ksp', ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Get the ksp context from master_PC which is the actual preconditioner:
  call PCKSPGetKSP(master_PC, master_PC_KSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! master_PC_KSP type will always be of type richardson. If the
  ! number  of iterations is set to 1, this ksp object is transparent. 
  call KSPSetType(master_PC_KSP, 'richardson', ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Important to set the norm-type to None for efficiency.
  call kspsetnormtype(master_PC_KSP, KSP_NORM_NONE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Do one iteration of the outer ksp preconditioners. Note the
  ! tolerances are unsued since we have set KSP_NORM_NON
  call KSPSetTolerances(master_PC_KSP, PETSC_DEFAULT_DOUBLE_PRECISION, &
       PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, &
       outerPreConIts, ierr)

  ! Get the 'preconditioner for master_PC_KSP, called 'pc'. This
  ! preconditioner is potentially run multiple times. 
  call KSPgetPC(master_PC_KSP, pc, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the type of 'pc'. This will almost always be additive schwarty
  call PCSetType(pc, 'asm', ierr)!PreCondType, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the overlap required
  call PCASMSetOverlap(pc, overlap, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  !Setup the main ksp context before extracting the subdomains
  call KSPSetUp(ksp, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Extract the ksp objects for each subdomain
  call PCASMGetSubKSP(pc, nlocal, first, subksp, ierr )
  call EChk(ierr, __FILE__, __LINE__)

  ! This 'subksp' object will ALSO be of type richardson so we can do
  ! multiple iterations on the sub-domains
  call KSPSetType(subksp, 'richardson', ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Again, norm_type is NONE since we don't want to check error
  call kspsetnormtype(subksp, KSP_NORM_NONE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the number of iterations to do on local blocks. Tolerances are ignored. 
  call KSPSetTolerances(subksp, PETSC_DEFAULT_DOUBLE_PRECISION, &
       PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, &
       innerPreConIts, ierr)

  call EChk(ierr, __FILE__, __LINE__)

  ! Extract the preconditioner for subksp object.
  call KSPGetPC(subksp, subpc, ierr )
  call EChk(ierr, __FILE__, __LINE__)

  ! The subpc type will almost always be ILU
  call PCSetType(subpc, LocalPCType, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Setup the matrix ordering for the subpc object:
  call PCFactorSetMatOrderingtype(subpc, matrixOrdering, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the ILU parameters
  call PCFactorSetLevels(subpc, fillLevel , ierr)
  call EChk(ierr, __FILE__, __LINE__) 

  ! Setup monitor if necessary:
  if (setMonitor) then
     call KSPMonitorSet(ksp, MyKSPMonitor, PETSC_NULL_OBJECT, &
          PETSC_NULL_FUNCTION, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  endif

  ! Send some feedback to screen
  if( myid==0 ) then
     write(*, 20) adjRelTol, adjAbsTol, adjDivTol
     write(*, 21) adjMaxIter
     call PCGetType(pc, pcType, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     write(*, 30) pcType
  endif

  ! Output formats.
10 format("# ... KSP properties:", /, &
       "#", 7x, "type        :", 1x, a)
20 format("#", 7x, "tolerances  :", 1x, "rel =", 1x, es7.1, /, &
       "#", 21x, "abs =", 1x, es7.1, /, &
       "#", 21x, "div =", 1x, es7.1)
21 format("#", 7x, "max.iter.   :", 1x, i4)
30 format("#", 7x, "precond.type:", 1x, a)
40 format(a, 1x, i3, a, 1x, i6, a, 1x, i6, 1x, a, 1x, i6) ! ownership
#endif

end subroutine setupPETScKsp
!
!     ******************************************************************
!
subroutine MyKSPMonitor(myKsp, n, rnorm, dummy, ierr)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * This is a user-defined routine for monitoring the KSP          *
  !     * iterative solvers. Instead of outputing the L2-norm at every   *
  !     * iteration (default PETSc monitor), it only does it every       *
  !     * 'adjMonStep' iterations.                                       *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc
  use inputADjoint
  use communication
  implicit none
  !
  !     Subroutine arguments.
  !
  ! myKsp - Iterative context
  ! n     - Iteration number
  ! rnorm - 2-norm (preconditioned) residual value
  ! dummy - Optional user-defined monitor context (unused here)
  ! ierr  - Return error code

  real(kind=realType), pointer, dimension(:, :) :: myKsp
  integer(kind=intType) :: n, dummy, ierr
  real(kind=realType)   :: rnorm
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !
#ifndef USE_NO_PETSC

  ! Write the residual norm to stdout every adjMonStep iterations.

  if( mod(n, adjMonStep)==0 ) then
     if( myid==0 ) write(*, 10) n, rnorm
  end if

  ierr = 0

  ! Output format.

10 format(i4, 1x, 'KSP Residual norm', 1x, e16.10)

#endif

end subroutine MyKSPMonitor

