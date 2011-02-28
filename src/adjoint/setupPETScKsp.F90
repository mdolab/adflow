!
!     ******************************************************************
!     *                                                                *
!     * File:          setupPETScKsp.F90                               *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 11-05-2010                                      *
!     * Last modified: 11-05-2010                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine setupPETScKsp
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create the Krylov subspace linear solver context,              *
  !     * the preconditioner context and set their various options.      *
  !     * This defines the linear solver to be used to solve the adjoint *
  !     * system of equations.                                           *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc
  use ADjointVars
  use blockPointers
  use flowVarRefState  !nw
  use inputADjoint
  use communication
  implicit none

  !
  !     User-defined Fortran routine.
  !
  external MyKSPMonitor

  !
  !     Local variables.
  !
  real(kind=realType)   :: rTol, aTol, dTol
  integer(kind=intType) :: mIts
  character(len=10)     :: kspType, pcType

#ifndef USE_NO_PETSC

  ! PETSc macros are lost and have to be redefined.
  ! They were extracted from: <PETSc root dir>/include/petscksp.h
  !                                                   /petscpc.h
  !Solvers
#define KSPGMRES      "gmres"   
#define KSPBCGS       "bcgs"
#define KSPCG         "cg"
#define KSPFGMRES     "fgmres"

  !Global Preconditioners
#define PCJACOBI      "jacobi"
#define PCBJACOBI     "bjacobi"
#define PCASM         "asm"

  !Local Preconditioners
#define PCILU         "ilu"
#define PCICC         "icc"
#define PCLU          "lu"
#define PCCHOLESKY    "cholesky"

  !Matrix Reorderings
#define MATORDERING_NATURAL       "natural"
#define MATORDERING_RCM       "rcm"
#define MATORDERING_ND        "nd"
#define MATORDERING_1WD       "1wd"
#define MATORDERING_QMD       "qmd"

  !Other miscellaneous definitions
#define PETSC_NULL           0
#define PETSC_DEFAULT        -2
#define KSPPREONLY    "preonly"

  if (ApproxPC)then

     !setup the approximate PC Matrix
     !call setupADjointPCMatrix(level)
     call setupADjointPCMatrixTranspose()

     !now set up KSP Context
     !call KSPSetOperators(ksp,dRdW,dRdWPre, &
     !                  DIFFERENT_NONZERO_PATTERN,PETScIerr)
     call KSPSetOperators(ksp,dRdWT,dRdWPreT, &
          DIFFERENT_NONZERO_PATTERN,PETScIerr)
     call EChk(PETScIerr,__file__,__line__)
  else

     ! Use the exact jacobian.
     ! Here the matrix that defines the linear system
     ! also serves as the preconditioning matrix.

     call KSPSetOperators(ksp,dRdWT,dRdWT, &
          DIFFERENT_NONZERO_PATTERN,PETScIerr)
     call EChk(PETScIerr,__file__,__line__)
  end if

  !     ******************************************************************
  !     *                                                                *
  !     * Set the various options for the KSP.                           *
  !     *                                                                *
  !     ******************************************************************

  call KSPSetFromOptions(ksp, PETScIerr)
  call EChk(PETScIerr,__file__,__line__)


  !     *****************************************************************
  !     *                                                               *
  !     * Now setup the specific options for the selected solver        *
  !     *                                                               *
  !     *****************************************************************
  !

  select case(ADjointSolverType)

  case(PETSCGMRES)

     call KSPSetType(ksp, KSPGMRES, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     call KSPGMRESSetRestart(ksp, adjRestart, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     call KSPGMRESSetCGSRefinementType(ksp, & 
          KSP_GMRES_CGS_REFINE_IFNEEDED, PETScIerr)   
     call EChk(PETScIerr,__file__,__line__)

     select case(PCSide)
     case(Right)
        call KSPSetPreconditionerSide(ksp, PC_RIGHT, PETScIerr)
        call EChk(PETScIerr,__file__,__line__)
     case(Left)
        call KSPSetPreconditionerSide(ksp, PC_LEFT, PETScIerr)
        call EChk(PETScIerr,__file__,__line__)
     end select

  case(PETSCBICGStab)

     call KSPSetType(ksp, KSPBCGS, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     select case(PCSide)
     case(Right)
        call KSPSetPreconditionerSide(ksp, PC_RIGHT, PETScIerr)
        call EChk(PETScIerr,__file__,__line__)
     case(Left)
        call KSPSetPreconditionerSide(ksp, PC_LEFT, PETScIerr)
        call EChk(PETScIerr,__file__,__line__)
     end select

  case(PETSCCG)

     call KSPSetType(ksp, KSPCG, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     select case(PCSide)
     case(Right)
        call KSPSetPreconditionerSide(ksp, PC_RIGHT, PETScIerr)
        call EChk(PETScIerr,__file__,__line__)
     case(Left)
        call KSPSetPreconditionerSide(ksp, PC_LEFT, PETScIerr)
        call EChk(PETScIerr,__file__,__line__)
     end select

  case(PETSCFGMRES)

     call KSPSetType(ksp, KSPFGMRES, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     call KSPGMRESSetRestart(ksp, adjRestart, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     select case(PCSide)
     case(Right)
        call KSPSetPreconditionerSide(ksp, PC_RIGHT, PETScIerr)
        call EChk(PETScIerr,__file__,__line__)
     case(Left)
        call KSPSetPreconditionerSide(ksp, PC_LEFT, PETScIerr)
        call EChk(PETScIerr,__file__,__line__)
     end select
  end select

  ! Non-solver specific options

  call KSPSetTolerances(ksp, adjRelTol, adjAbsTol, adjDivTol, &
       adjMaxIter, PETScIerr)
  call EChk(PETScIerr,__file__,__line__)

  if( myid==0 ) then
     call KSPGetTolerances(ksp, rTol, aTol, dTol, mIts, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)
     write(*,20) rTol, aTol, dTol
     write(*,21) mIts
  endif

  ! Setup Preconditioning 
  call KSPGetPC(ksp, pc, PETScIerr)
  call EChk(PETScIerr,__file__,__line__)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Select the preconditioning method.                             *
  !     *                                                                *
  !     ****************************************************************** 
  !

  select case(PreCondType)

  case(Jacobi)

     call PCSetType( pc, PCJACOBI, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)
     ! set tolerances on subksp (is this really needed???)
     call KSPSetTolerances(ksp, 1.e-8, adjAbsTol, adjDivTol, &
          adjMaxIter, PETScIerr)

     !set Scaling Factor Type
     select case(ScaleType)
     case(Normal)
        continue
     case(RowMax)
        call PCJacobiSetUseRowMax(pc, PETScIerr)
        call EChk(PETScIerr,__file__,__line__)
     case(RowSum)
        call PCJacobiSetUseRowSum(pc, PETScIerr)
        call EChk(PETScIerr,__file__,__line__)
     case(RowAbs)
        call PCJacobiSetUseAbs(pc, PETScIerr)
        call EChk(PETScIerr,__file__,__line__)
     end select

     select case(MatrixOrdering)
     case(Natural)
        call PCFactorSetMatOrderingtype( pc, MATORDERING_NATURAL, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(ReverseCuthillMckee)
        call PCFactorSetMatOrderingtype( pc, MATORDERING_RCM, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(NestedDissection)
        call PCFactorSetMatOrderingtype( pc, MATORDERING_ND, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(OnewayDissection )
        call PCFactorSetMatOrderingtype( pc, MATORDERING_1WD, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case( QuotientMinimumDegree)
        call PCFactorSetMatOrderingtype( pc, MATORDERING_QMD, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     end select

     !Set the iteration Monitor
     if (setMonitor) then
        call KSPMonitorSet(ksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
             PETSC_NULL_FUNCTION, PETScIerr)
     endif

  case(BlockJacobi)

     call PCSetType( pc, PCBJACOBI, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)
     call PCSetUp(pc,PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     Nsub = 1!_intType
     length = nCellsLocal*nw
     call PCBJacobiSetLocalBlocks(pc,Nsub,length,PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     !Setup KSP context before calling local subdomain ksp contexts
     call KSPSetUp(ksp, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     !Setup local subcontexts
     call PCBJacobiGetSubKSP(pc,nlocal,first,subksp,PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     !Setup Local ILU precondtioner
     call KSPGetPC( subksp, subpc, PETScIerr )
     call EChk(PETScIerr,__file__,__line__)

     select case(LocalPCType)
     case(ILU)
        call PCSetType( subpc, PCILU, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(ICC)
        call PCSetType( subpc, PCICC, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(LU)
        call PCSetType( subpc, PCLU, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(Cholesky)
        call PCSetType( subpc, PCCHOLESKY, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     end select

     !set matrix ordering

     select case(MatrixOrdering)
     case(Natural)
        call PCFactorSetMatOrderingtype( subpc, MATORDERING_NATURAL, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(ReverseCuthillMckee)
        call PCFactorSetMatOrderingtype( subpc, MATORDERING_RCM, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(NestedDissection)
        call PCFactorSetMatOrderingtype( subpc, MATORDERING_ND, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(OnewayDissection )
        call PCFactorSetMatOrderingtype( subpc, MATORDERING_1WD, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case( QuotientMinimumDegree)
        call PCFactorSetMatOrderingtype( subpc, MATORDERING_QMD, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     end select

     !Set ILU parameters
     call PCFactorSetLevels( subpc, fillLevel , PETScIerr)!
     call EChk(PETScIerr,__file__,__line__) 

     !Set local contexts to preconditioner's only
     call KSPSetType(subksp, KSPPREONLY, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     ! set tolerances on subksp (is this really needed???)
     call KSPSetTolerances(subksp, 1.e-8, adjAbsTol, adjDivTol, &
          adjMaxIter, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     if(setMonitor)then
        !Set the convergence monitors
        call KSPMonitorSet(subksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
             PETSC_NULL_FUNCTION, PETScIerr)
        call EChk(PETScIerr,__file__,__line__)
        call KSPMonitorSet(ksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
             PETSC_NULL_FUNCTION, PETScIerr)
        call EChk(PETScIerr,__file__,__line__)
     endif

  case(AdditiveSchwartz)

     call PCSetType( pc, PCASM, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     !setup a basic overlaping scheme, more detailed scheme to be used later
     call PCASMSetOverlap(pc,overlap,PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     !Setup KSP context before calling local subdomain ksp contexts
     call KSPSetUp(ksp, PETScIerr);
     call EChk(PETScIerr,__file__,__line__)

     !get the subdomain contexts
     call PCASMGetSubKSP( pc, nlocal,  first, subksp, PETScIerr )
     call EChk(PETScIerr,__file__,__line__)
        
     !Setup the subdomain preconditioner
     call KSPGetPC( subksp, subpc, PETScIerr )
     call EChk(PETScIerr,__file__,__line__)

     !Set subdomain preconditioner type

     select case(LocalPCType)
     case(ILU)
        call PCSetType( subpc, PCILU, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(ICC)
        call PCSetType( subpc, PCICC, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(LU)
        call PCSetType( subpc, PCLU, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(Cholesky)
        call PCSetType( subpc, PCCHOLESKY, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     end select

     !set matrix ordering

     select case(MatrixOrdering)
     case(Natural)
        call PCFactorSetMatOrderingtype( subpc, MATORDERING_NATURAL, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(ReverseCuthillMckee)
        call PCFactorSetMatOrderingtype( subpc, MATORDERING_RCM, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(NestedDissection)
        call PCFactorSetMatOrderingtype( subpc, MATORDERING_ND, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case(OnewayDissection )
        call PCFactorSetMatOrderingtype( subpc, MATORDERING_1WD, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     case( QuotientMinimumDegree)
        call PCFactorSetMatOrderingtype( subpc, MATORDERING_QMD, PETScIerr )
        call EChk(PETScIerr,__file__,__line__)
     end select

     !Set ILU parameters
     call PCFactorSetLevels( subpc, fillLevel , PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     !Set local contexts to preconditioner's only
     call KSPSetType(subksp, KSPPREONLY, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     ! set tolerances on subksp (is this really needed???)
     call KSPSetTolerances(subksp, 1.e-8, adjAbsTol, adjDivTol, &
          adjMaxIter, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)

     if(setMonitor)then
        call KSPMonitorSet(ksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
             PETSC_NULL_FUNCTION, PETScIerr)
        call EChk(PETScIerr,__file__,__line__)
     endif
  end select

  if( myid==0 ) then
     call PCGetType(pc, pcType, PETScIerr)
     call EChk(PETScIerr,__file__,__line__)
     write(*,30) pcType
  endif

  ! Output formats.

10 format("# ... KSP properties:",/, &
       "#",7x,"type        :",1x,a)
20 format("#",7x,"tolerances  :",1x,"rel =",1x,es7.1,/, &
       "#",21x,"abs =",1x,es7.1,/, &
       "#",21x,"div =",1x,es7.1)
21 format("#",7x,"max.iter.   :",1x,i4)
30 format("#",7x,"precond.type:",1x,a)
40 format(a,1x,i3,a,1x,i6,a,1x,i6,1x,a,1x,i6) ! ownership
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

  real(kind=realType), pointer, dimension(:,:) :: myKsp
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

  if( mod(n,adjMonStep)==0 ) then
     if( myid==0 ) write(*,10) n, rnorm
  end if

  ierr = 0

  ! Output format.

10 format(i4,1x,'KSP Residual norm',1x,e16.10)

#endif

end subroutine MyKSPMonitor

