!
!     ******************************************************************
!     *                                                                *
!     * File:          setupPETScKsp.F90                               *
!     * Author:        Gaetan Kenway                                   *
!     * Starting date: 26-12-2012                                      *
!     * Last modified: 26-12-2012                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine setupPETScKsp

  use ADjointPETSc, only: drdwpret, drdwt, adjointKSP
  use ADjointVars
  use inputADjoint
  use communication
  use blockPointers
  implicit none

#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

  !     Local variables.
  logical :: useAD, usePC, useTranspose, useObjective
  integer(kind=intType) :: ierr

  PC master_PC, coarsePC, finePC, levelPC, subpc
  KSP coarseKSPSolver, fineKSPSolver, levelKSP, subksp
  Mat tmp1
  external MyKSPMonitor

  if (ApproxPC)then
     !setup the approximate PC Matrix
     useAD = ADPC
     useTranspose = .True.
     usePC = .True.
     useObjective = .False.
     call setupStateResidualMatrix(drdwpret, useAD, usePC, useTranspose, &
          useObjective, frozenTurbulence, 1_intType)
     call KSPSetOperators(adjointKSP, dRdwT, dRdWPreT, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  else
     ! Use the exact jacobian.  Here the matrix that defines the
     ! linear system also serves as the preconditioning matrix. This
     ! is only valid if useMatrixFree is flase. 
     if (useMatrixfreedRdw) then 
        call terminate("setupPETScKSP", "useMatrixFreedRdW option cannot be true when the approxPC option is False")
     end if
     call KSPSetOperators(adjointKSP, dRdWt, dRdWT, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

  if (PreCondType == 'asm') then
     ! Run the super-dee-duper function to setup the ksp object:
     call setupStandardKSP(adjointKSP, ADjointSolverType, adjRestart, adjointpcside, &
          PreCondType, overlap, outerPreConIts, localPCType, &
          matrixOrdering, FillLevel, innerPreConIts)
  else if (PreCondType == 'mg') then
     print *,'Only ASM precondtype is usable'
     stop
  end if

  ! Setup monitor if necessary:
  if (setMonitor) then
     call KSPMonitorSet(adjointKSP, MyKSPMonitor, PETSC_NULL_OBJECT, &
          PETSC_NULL_FUNCTION, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  endif

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

  if(mod(n, adjMonStep) ==0 ) then
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
  !
  ! Note that if globalPreConIts=1 then maser_PC_KSP is NOT created and master_PC=globalPC
  ! and if localPreConIts=1 then subKSP is set to preOnly. 
  use constants
  
  implicit none

#ifndef USE_NO_PETSC
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

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

  if (trim(kspObjectType) == 'richardson') then
     call KSPSetPCSide(kspObject, PC_LEFT, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if

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
     call KSPSetTolerances(master_PC_KSP, PETSC_DEFAULT_REAL, &
          PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, &
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

     call KSPSetTolerances(subksp, PETSC_DEFAULT_REAL, &
          PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, &
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
