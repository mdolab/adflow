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

  use ADjointPETSc, only: drdwpret, drdwt, adjointKSP
  use ADjointPETSc, only : coarsedRdWPreT, restrictionOperator
  use ADjointPETSc, only : prolongationOperator
  use ADjointVars
  use inputADjoint
  use communication
  use blockPointers
  implicit none
#ifndef USE_NO_PETSC
#define PETSC_AVOID_MPIF_H
#include "finclude/petsc.h"

  !     Local variables.
  logical :: useAD, usePC, useTranspose 
  integer(kind=intType) :: ierr, nLevels, i, l
  integer(kind=intType) :: nlocal, first
  integer(kind=intType), allocatable, dimension(:) :: comms

  PC master_PC, coarsePC, finePC, levelPC, subpc
  KSP coarseKSPSolver, fineKSPSolver, levelKSP, subksp
  external MyKSPMonitor

  if (ApproxPC)then
     !setup the approximate PC Matrix
     useAD = .False.
     useTranspose = .True.
     usePC = .True.
     call setupStateResidualMatrix(drdwpret, useAD, usePC, useTranspose, &
          1_intType)

     !now set up KSP Context
     call KSPSetOperators(adjointKSP, dRdWT, dRdWPreT, &
          DIFFERENT_NONZERO_PATTERN, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  else
     ! Use the exact jacobian.  Here the matrix that defines the
     ! linear system also serves as the preconditioning matrix.
     call KSPSetOperators(adjointKSP, dRdWT, dRdWT, &
          SAME_NONZERO_PATTERN, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  ! Assemble the coarse grid operators for the mg preconditioner:
!   if (PreCondType == 'mg') then
!      nLevels = ubound(flowDoms, 2)
!      useAD = .False.
!      useTranspose = .True.
!      usePC = .True.
     
!      ! Setup the restriction matrices and the residual matrices on
!      ! the lower grid levels
!      do i=2, nLevels
!         if (myid == 0) then
!            print *, 'Assembling coarse-grid interpolation and residual matrix for Level:',i
!         end if
!         call setupStateResidualMatrix(coarsedRdWpreT(i), useAD, usePC, useTranspose, i)
!         call setupRestrictionMatrix(restrictionOperator(i), i)
!         call setupProlongationMatrix(prolongationOperator(i), i)
!      end do

!      ! Neet to convert the coarse grid to aij format inorder to use superlu_dist
!      print *,'Converting Coarse...'
!      call matConvert(coarsedRdWPreT(nLevels), "aij", MAT_REUSE_MATRIX, &
!           coarsedRdWPreT(nLevels), ierr)
!      call EChk(ierr, __FILE__, __LINE__)
!      print *,'Done converting.'
  
!      ! The mg preconditioner gets exact ksp operators regardless of
!      ! approxPC variable
!      call KSPSetOperators(adjointKSP, dRdWT, dRdWT, &
!           SAME_NONZERO_PATTERN, ierr)
!      call EChk(ierr, __FILE__, __LINE__)
!   end if
  
  ! Set the main tolerance for th main adjoint KSP solver
  call KSPSetTolerances(adjointKSP, adjRelTol, adjAbsTol, adjDivTol, &
       adjMaxIter, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  if (PreCondType == 'asm') then
     ! Run the super-dee-duper function to setup the ksp object:
     call setupStandardKSP(adjointKSP, ADjointSolverType, adjRestart, adjointpcside, &
          PreCondType, overlap, outerPreConIts, localPCType, &
          matrixOrdering, FillLevel, innerPreConIts)
  else if (PreCondType == 'mg') then

     print *, 'Converting to aij'
     call MatConvert(dRdwT, 'aij', MAT_REUSE_MATRIX, drdwt, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     print *,'Done converting to aij'

     call KSPSetFromOptions(adjointKSP, ierr)
     call KSPGMRESSetRestart(adjointKSP, adjRestart, ierr)
     call KSPGMRESSetCGSRefinementType(adjointKSP, &
       KSP_GMRES_CGS_REFINE_IFNEEDED, ierr)   
     call KSPSetPCSide(adjointKSP, PC_RIGHT, ierr)
     call KSPGetPC(adjointKSP, master_PC, ierr)
     call PCSetType(master_PC, PCML, ierr)
     !call pcMGsetNumberSmoothDown(master_PC, 10, ierr)
     !call EChk(ierr, __FILE__, __LINE__)
  end if

  ! Setup monitor if necessary:
  if (setMonitor) then
     call KSPMonitorSet(adjointKSP, MyKSPMonitor, PETSC_NULL_OBJECT, &
          PETSC_NULL_FUNCTION, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  endif

  ! Send some feedback to screen
  if( myid==0 ) then
     write(*, 20) adjRelTol, adjAbsTol, adjDivTol
     write(*, 21) adjMaxIter
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

subroutine setupStandardKSP(kspObject, kspObjectType, gmresRestart, preConSide, &
     globalPCType, ASMOverlap, globalPreConIts, localPCType, &
     localMatrixOrdering, localFillLevel, localPreConIts)

  ! This function sets up the supplied kspObject in the followin
  ! specific fashion. The reason this setup is in
  ! its own function is that it is used in the following places:
  ! 1. Setting up the preconditioner to use for the NKsolver
  ! 2. Setting up the preconditioner to use for the adjoint solver
  ! 3. Setting up the smoothers on the coarse multigrid levels. 
  ! 
  ! The hierarchy of the setup is:
  !  kspObject --> Supplied KSP object
  !  |
  !  --> master_PC --> Preconditioner type set to KSP
  !      |
  !      --> master_PC_KSP --> KSP type set to Richardson with 'globalPreConIts'
  !          |
  !           --> globalPC --> PC type set to 'globalPCType'
  !               |            Usually Additive Schwartz and overlap is set
  !               |            with 'ASMOverlap'. Use 0 to get BlockJacobi
  !               |
  !               --> subKSP --> KSP type set to Richardon with 'LocalPreConIts'
  !                   |
  !                   --> subPC -->  PC type set to 'loclaPCType'.
  !                                  Usually ILU. 'localFillLevel' is 
  !                                  set and 'localMatrixOrder' is used.

  use constants
  
  implicit none

#ifndef USE_NO_PETSC
#define PETSC_AVOID_MPIF_H
#include "finclude/petsc.h"

  ! Input Params
  KSP kspObject
  character(len=maxStringLen), intent(in) :: kspObjectType, preConSide
  character(len=maxStringLen), intent(in) :: globalPCType, localPCType
  character(len=maxStringLen), intent(in) :: localMatrixOrdering
  integer(kind=intType), intent(in) :: ASMOverlap, localFillLevel, gmresRestart
  integer(kind=intType), intent(in) :: globalPreConIts, localPreConIts

  ! Working Variables
  PC  master_PC, globalPC, subpc
  KSP master_PC_KSP, subksp
  character(len=10)  :: pcType
  integer(kind=intType) :: nlocal, first, ierr

  ! First, KSPSetFromOptions MUST be called
  call KSPSetFromOptions(kspObject, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the type of solver to use:
  call KSPSetType(kspObject, kspObjectType, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! If we're using GMRES set the possible gmres restart
  call KSPGMRESSetRestart(kspObject, gmresRestart, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! If you're using GMRES, set refinement type
  call KSPGMRESSetCGSRefinementType(kspObject, &
       KSP_GMRES_CGS_REFINE_IFNEEDED, ierr)   
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the preconditioner side from option:
  if (trim(preConSide) == 'right') then
     call KSPSetPCSide(kspObject, PC_RIGHT, ierr)
  else
     call KSPSetPCSide(kspObject, PC_LEFT, ierr)
  end if
  call EChk(ierr, __FILE__, __LINE__)

  ! Since there is an extraneous matMult required when using the
  ! richardson precondtiter with only 1 iteration, only use it we need
  ! to do more than 1 iteration.
  if (globalPreConIts > 1) then
     ! Extract preconditioning context for main KSP solver: (master_PC)
     call KSPGetPC(kspObject, master_PC, ierr)
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
          globalPreConIts, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     ! Get the 'preconditioner for master_PC_KSP, called 'globalPC'. This
     ! preconditioner is potentially run multiple times. 
     call KSPgetPC(master_PC_KSP, globalPC, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  else
     ! Just pull out the pc-object if we are not using kspRichardson
     call KSPGetPC(kspObject, globalPC, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  ! Set the type of 'globalPC'. This will almost always be additive schwartz
  call PCSetType(globalPC, 'asm', ierr)!globalPCType, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the overlap required
  call PCASMSetOverlap(globalPC, ASMOverlap, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  !Setup the main ksp context before extracting the subdomains
  call KSPSetUp(kspObject, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Extract the ksp objects for each subdomain
  call PCASMGetSubKSP(globalPC, nlocal, first, subksp, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Since there is an extraneous matMult required when using the
  ! richardson precondtiter with only 1 iteration, only use it we need
  ! to do more than 1 iteration.
  if (localPreConIts > 1) then
     ! This 'subksp' object will ALSO be of type richardson so we can do
     ! multiple iterations on the sub-domains
     call KSPSetType(subksp, 'richardson', ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Set the number of iterations to do on local blocks. Tolerances are ignored. 
     call KSPSetTolerances(subksp, PETSC_DEFAULT_DOUBLE_PRECISION, &
          PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, &
          localPreConIts, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Again, norm_type is NONE since we don't want to check error
     call kspsetnormtype(subksp, KSP_NORM_NONE, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  else
     call KSPSetType(subksp, 'preonly', ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  ! Extract the preconditioner for subksp object.
  call KSPGetPC(subksp, subpc, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! The subpc type will almost always be ILU
  call PCSetType(subpc, localPCType, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Setup the matrix ordering for the subpc object:
  call PCFactorSetMatOrderingtype(subpc, localMatrixOrdering, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the ILU parameters
  call PCFactorSetLevels(subpc, localFillLevel , ierr)
  call EChk(ierr, __FILE__, __LINE__) 

#endif
end subroutine setupStandardKSP


! subroutine mgsolver()

!   use ADjointPETSc, only: drdwpret, drdwt, adjointKSP
!   use ADjointPETSc, only : coarsedRdWPreT, restrictionOperator
!   use ADjointPETSc, only : prolongationOperator, dJdw, adjointRes, psi
!   use ADjointVars
!   use inputADjoint
!   use communication
!   use blockPointers
!   implicit none
! #ifndef USE_NO_PETSC
! #define PETSC_AVOID_MPIF_H
! #include "finclude/petsc.h"

!   !     Local variables.
!   logical :: useAD, usePC, useTranspose 
!   integer(kind=intType) :: ierr, nLevels, i, l, m, n
!   real(kind=realType) :: value
!   Vec coarseR, coarsePsi, coarseF
!   PC coarsePC
!   KSP coarseKSP


!   useAD = .False.
!   useTranspose = .True.
!   usePC = .True.
  
!   print *,'Setting coarse stuff...'
!   !call setupStateResidualMatrix(coarsedRdWpreT(2), useAD, usePC, useTranspose, 2)
!   call setupStateResidualMatrix(dRdwPreT, useAD, usePC, useTranspose, 1)
!   call setupRestrictionMatrix(restrictionOperator(2), 2)
!   call setupProlongationMatrix(prolongationOperator(2), 2)


!   call matConvert(dRdwPreT, "aij", MAT_REUSE_MATRIX, drdwpret, ierr)
!   call EChk(ierr, __FILE__, __LINE__)
!   call matConvert(prolongationOperator(2), "aij", MAT_REUSE_MATRIX, prolongationOperator(2), ierr)
!   call EChk(ierr, __FILE__, __LINE__)

!   call MatPtAP(dRdwPreT, prolongationOperator(2), MAT_INITIAL_MATRIX, 5.0, coarsedRdwPreT(2), ierr)
!   !call MatRARt(dRdwPreT, restrictionOperator(2), MAT_INITIAL_MATRIX, 5.0, coarsedRdwPreT(2), ierr)
!   call EChk(ierr, __FILE__, __LINE__)
!   print *, 'Done setting coarse stuff.'

!   call MatGetVecs(coarsedRdWPreT(2), coarseR, coarsePsi, ierr)
!   call EChk(ierr, __FILE__, __LINE__)
  
!   call VecDuplicate(coarseR, coarseF, ierr)
!   call EChk(ierr,  __FILE__, __LINE__)


!   ! Need to setup the corase direct solve
!   call KSPCreate(SUMB_PETSC_COMM_WORLD, coarseKSP, ierr)
!   call EChk(ierr, __FILE__, __LINE__)

! !   ! Convert coarse grid to aij
! !   call matConvert(coarsedRdWPreT(2), "aij", MAT_REUSE_MATRIX, &
! !        coarsedRdWPreT(2), ierr)
! !   call EChk(ierr, __FILE__, __LINE__)
!   print *,'Done ksp create'
!   ! Set coarse grid operators
!   call KSPSetOperators(coarseKSP, coarsedRdwPreT(2), &
!        coarsedRdwPreT(2), SAME_NONZERO_PATTERN, ierr)
!   call EChk(ierr, __FILE__, __LINE__)
!   print *,'done set ops'

!   call KSPSetFromOptions(coarseKSP, ierr)
!   call EChk(ierr, __FILE__, __LINE__)

!   call KSPSetType(coarseKSP, 'preonly', ierr)
!   call EChk(ierr, __FILE__, __LINE__, ierr)

! !   call KSPSetPCSide(coarseKSP, PC_RIGHT, ierr)
! !   call EChk(ierr, __FILE__, __LINE__)

!   ! Extract the PC so we can setup superLU_dist
!   call KSPGetPC(coarseKSP, coarsePC, ierr)
!   call EChk(ierr, __FILE__, __LINE__, ierr)

!   ! Set the PC type to 'lu'
!   call PCSetType(coarsePC, "lu", ierr)
!   call EChk(ierr, __FILE__, __LINE__)

!   ! Set 'lu' solver to SuperLU_dist
!   call PCFactorSetMatSolverPackage(coarsePC, "superlu_dist", ierr)
!   call EChk(ierr, __FILE__, __LINE__)
!   print *,'before setup'
!   call KSPSetup(coarseKSP, ierr)
!   call EChk(ierr, __FILE__, __LINE__)
!   print *,'done setup'
!   ! Make sure the ksp is right
!   call KSPView(coarseKSP, PETSC_VIEWER_STDOUT_WORLD, ierr)
!   call EChk(ierr, __FILE__, __LINE__)

!   ! Set zero initial guess
!   call VecSet(psi, zero, ierr)
!   call EChk(ierr,  __FILE__, __LINE__)

!   ! Compute the initial residual:
!   call MatMult(dRdwt, psi, adjointRes, ierr)
!   call EChk(ierr,  __FILE__, __LINE__)

!   call VecAYPX(adjointRes, -one, dJdw, ierr)
!   call EChk(ierr,  __FILE__, __LINE__)

!   call VecScale(adjointRes, -one, ierr)
!   call EChk(ierr,  __FILE__, __LINE__)

!   ! Start iteration loop:
!   do i=1,20
     
!      ! Compute the norm of this residual:
!      call VecNorm(adjointRes, NORM_2, value, ierr)
!      call EChk(ierr,  __FILE__, __LINE__)

!      if (myid == 0) then
!         print *,'Start Iter, value:',i, value
!      end if

!      ! Restrict the reisudal to coarse level:
!      call VecScale(adjointRes,1/8.0, ierr)

!      call MatMultTranspose(prolongationOperator(2), adjointRes, coarseR, ierr)
!      call EChk(ierr,  __FILE__, __LINE__)

!      ! Now solve on the coarse grid:
!      call kspSolve(coarseKSP, coarseR, coarsePsi, ierr)
!      call EChk(ierr,  __FILE__, __LINE__)

! !      ! Check coarse grid residual:
! !      call matmult(coarsedrdwpret(2),coarsePsi, coarseF, ierr)
! !      call vecAYPX(coarseF, -one, coarseR, ierr)

! !      call VecNorm(coarseF, NORM_2, value, ierr)
! !      call EChk(ierr,  __FILE__, __LINE__)

! !      if (myid == 0) then
! !         print *,'  coarse-sol res:',value
! !      end if



!      ! Prolong the solution back to fine level
!       call MatMult(prolongationOperator(2), coarsePsi, adjointRes, ierr)
!       call EChk(ierr,  __FILE__, __LINE__)

! !      call MatMultTranspose(restrictionOperator(2), coarsePsi, adjointRes, ierr)
! !      call EChk(ierr,  __FILE__, __LINE__)

!      ! Advance the psi solution with the correction stored in adjointRes
!       call VecAXPY(psi, one, adjointRes, ierr)
!       call EChk(ierr,  __FILE__, __LINE__)

!      ! Compute the new residual:
!      call MatMult(dRdwt, psi, adjointRes, ierr)
!      call EChk(ierr,  __FILE__, __LINE__)
     
!      call VecAYPX(adjointRes, -one, dJdw, ierr)
!      call EChk(ierr,  __FILE__, __LINE__)

!      call VecScale(adjointRes, -one, ierr)
!      call EChk(ierr,  __FILE__, __LINE__)

! !      ! Check norm at end as well
! !      call VecNorm(adjointRes, NORM_2, value, ierr)
! !      call EChk(ierr,  __FILE__, __LINE__)

! !      if (myid == 0) then
! !         print *,'End Iter, value:',i, value
! !      end if
     
!   end do


!   ! Destroy the coarse vectors. 
!   call VecDestroy(coarseR, ierr)
!   call VecDestroy(coarseF, ierr)
!   call VecDestroy(coarsePsi, ierr)


! #endif
! end subroutine mgsolver



! ------------------------- BACKUP ------------

!      ! We have more work to do since we have to setup coarse grid ops and the like.
!      ! First, KSPSetFromOptions MUST be called
!      print *,'Setting up multigrid preconditioner'
!      call KSPSetFromOptions(adjointKSP, ierr)
!      call EChk(ierr, __FILE__, __LINE__)

!      ! Set the type of solver to use:
!      call KSPSetType(adjointKSP, ADjointSolverType, ierr)
!      call EChk(ierr, __FILE__, __LINE__)

!      ! If we're using GMRES set the possible gmres restart
!      call KSPGMRESSetRestart(adjointKSP, adjRestart, ierr)
!      call EChk(ierr, __FILE__, __LINE__)

!      ! If you're using GMRES, set refinement type
!      call KSPGMRESSetCGSRefinementType(adjointKSP, &
!           KSP_GMRES_CGS_REFINE_IFNEEDED, ierr)   
!      call EChk(ierr, __FILE__, __LINE__)
     
!      ! Set the preconditioner side from option:
!      call KSPSetPCSide(adjointKSP, PC_RIGHT, ierr)
!      call EChk(ierr, __FILE__, __LINE__)
     
!      ! Extract preconditioning context for main KSP solver: (master_PC)
!      call KSPGetPC(adjointKSP, master_PC, ierr)
!      call EChk(ierr, __FILE__, __LINE__)

!      ! Explictly set to mg
!      call PCSetType(master_PC, 'mg', ierr)
!      call EChk(ierr, __FILE__, __LINE__, ierr)

!      call PCMGSetLevels(master_PC, nLevels, PETSC_NULL_OBJECT, ierr)
!      call EChk(ierr, __FILE__, __LINE__, ierr)

!      ! Set to multiplicative (default)
!      call PCMGSetType(master_PC, PC_MG_MULTIPLICATIVE, ierr)
!      call EChk(ierr, __FILE__, __LINE__, ierr)

!      ! For now, hard-code a V-cycle...probably best anyway
!      call PCMGSetCycleType(master_PC, PC_MG_CYCLE_V, ierr)
!      call EChk(ierr, __FILE__, __LINE__, ierr)

!      ! Set the number of application of the MG we use
!      !  call PCMGMultiplicativeSetCycles(master_PC, niter, ierr)
!      !  call EChk(ierr, __FILE__, __LINE__, ierr)

!      ! Now loop over the levels to supply the restriction
!      ! operators. NOTE: Petsc labels the mg level OPPOSITE of sumb. In
!      ! petsc, grid level 0 is coarsest, up to nLevels-1 is finest.
!      do l=1,nLevels-1
!          call PCMGSetRestriction(master_PC, l, &
!               restrictionOperator(nLevels-l+1), ierr)
!          call EChk(ierr, __FILE__, __LINE__, ierr)

! !        call PCMGSetInterpolation(master_PC, l, &
! !             prolongationOperator(nLevels-l+1), ierr)
! !        call EChk(ierr, __FILE__, __LINE__, ierr)
!      end do
     
!      ! Now extract the coarse KSP operator
!      call PCMGGetCoarseSolve(master_PC, coarseKSPSolver, ierr)
!      call EChk(ierr, __FILE__, __LINE__, ierr)

!      ! Must call KSPSetup
! !      call KSPSetup(coarseKSPSolver, ierr)
! !      call EChk(ierr, __FILE__, __LINE__, ierr)

!      ! This will be a direct solve with superLU_dist so set preonly to ksp type
!      call KSPSetType(coarseKSPSolver, 'preonly', ierr)
!      call EChk(ierr, __FILE__, __LINE__, ierr)

!      call KSPSetPCSide(coarseKSPSolver, PC_RIGHT, ierr)
!      call EChk(ierr, __FILE__, __LINE__)

!      ! Set the operator for the coarse grid:
!      call KSPSetOperators(coarseKSPSolver, coarsedRdWPreT(nLevels), &
!           coarsedRdWPreT(nlevels), SAME_NONZERO_PATTERN, ierr)
!      call EChk(ierr, __FILE__, __LINE__, ierr)

!      ! Extract the PC so we can setup superLU_dist
!      call KSPGetPC(coarseKSPSolver, coarsePC, ierr)
!      call EChk(ierr, __FILE__, __LINE__, ierr)

!      ! Set the PC type to 'lu'
!      call PCSetType(coarsePC, "lu", ierr)
!      call EChk(ierr, __FILE__, __LINE__)

! !      ! Set 'lu' solver to SuperLU_dist
!      call PCFactorSetMatSolverPackage(coarsePC, "superlu_dist", ierr)
!      call EChk(ierr, __FILE__, __LINE__)
     
!      ! Use nested dissection ordering for direct solve as this reduces
!      ! fill for the direct solve.
!      call PCFactorSetMatOrderingType(coarsePC, "nd", ierr)
!      call EChk(ierr, __FILE__, __LINE__)
     
!      ! Now extract the finest grid smoother and set that up explictly
!      ! --- Basically the fineGrid ksp does absolutely nothing: no pre
!      ! or post smoothing. 
!      call PCMGGetSmoother(master_PC, nLevels-1, fineKSPSolver, ierr)
!      call EChk(ierr, __FILE__, __LINE__)

!      call KSPSetType(fineKSPSolver, 'richardson', ierr)
!      call EChk(ierr, __FILE__, __LINE__)

!      ! Must call KSP setup
! !      call KSPSetup(fineKSPSolver, ierr)
! !      call EChk(ierr, __FILE__, __LINE__, ierr)

!      ! Set the norm type to none on the fine grid and also skip the
!      ! convergence test. 
! !     call KSPSetNormType(fineKSPSolver, KSP_NORM_NONE, ierr)
! !     call EChk(ierr, __FILE__, __LINE__)

!      call KSPGetPC(fineKSPSolver, finePC, ierr)
!      call EChk(ierr, __FILE__, __LINE__)

!      call PCSetType(finePC, PCNONE, ierr)
!      call EChk(ierr, __FILE__, __LINE__)

!      ! Now we have to go though the remainder of the interior levels
!      ! (nlevels-2) of them extract the smoother objects and do
!      ! something with them.

!      do l=2,nLevels-1

!         call PCMGGetSmoother(master_PC, l, levelKSP, ierr)
!         call EChk(ierr, __FILE__, __LINE__)
        
!         ! Must call KSPSetup for this level
! !         call KSPSetup(levelKSP, ierr)
! !         call EChk(ierr, __FILE__, __LINE__)

!         ! Set it to richardson
!         call KSPSetType(levelKSP, 'richardson', ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! Set the operator for this level
!         call KSPSetOperators(levelKSP, coarsedRdWPreT(nLevels-1+1), &
!           coarsedRdWPreT(nlevels-1+1), SAME_NONZERO_PATTERN, ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! Set the norm type to none on the fine grid and also skip the
!         ! convergence test. 
!         call KSPSetNormType(levelKSP, KSP_NORM_NONE, ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! Do one iteration of the outer ksp preconditioners. Note the
!         ! tolerances are unsued since we have set KSP_NORM_NON
!         call KSPSetTolerances(levelKSP, PETSC_DEFAULT_DOUBLE_PRECISION, &
!              PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, &
!              1_intType, ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! Extract PC for this level
!         call KSPGetPC(levelKSP, levelPC, ierr)
!         call EChk(ierr, __FILE__, __LINE__)
        
!         ! Set it to asm
!         call PCSetType(levelPC, 'asm', ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! Overlap
!         call PCASMSetOverlap(levelPC, 0, ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! Setup so we can get the subobjects
!         call KSPSetUp(levelKSP, ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! Get the subksp
!         call PCASMGetSubKSP(levelPC, nlocal, first, subksp, ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         call KSPSetType(subksp, 'richardson', ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! Set the number of iterations to do on local blocks. Tolerances are ignored. 
!         call KSPSetTolerances(subksp, PETSC_DEFAULT_DOUBLE_PRECISION, &
!              PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, &
!              1_intType, ierr)
!         call EChk(ierr, __FILE__, __LINE__)
        
!         ! Again, norm_type is NONE since we don't want to check error
!         call kspsetnormtype(subksp, KSP_NORM_NONE, ierr)
!         call EChk(ierr, __FILE__, __LINE__)
        
!         ! Extract the preconditioner for subksp object.
!         call KSPGetPC(subksp, subpc, ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! The subpc type will almost always be ILU
!         call PCSetType(subpc, 'ilu', ierr)
!         call EChk(ierr, __FILE__, __LINE__)

!         ! Setup the matrix ordering for the subpc object:
!         call PCFactorSetMatOrderingtype(subpc, 'rcm', ierr)
!         call EChk(ierr, __FILE__, __LINE__)
        
!         ! Set the ILU parameters
!         call PCFactorSetLevels(subpc, 1_intType, ierr)
!         call EChk(ierr, __FILE__, __LINE__) 
!      end do

!      call PCView(master_PC, PETSC_VIEWER_STDOUT_WORLD, ierr)
!      call EChk(ierr, __FILE__, __LINE__) 
