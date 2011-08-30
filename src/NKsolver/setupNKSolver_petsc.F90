subroutine setupNKsolver_petsc

  ! Setup the PETSc objects for the Newton-Krylov
  ! solver. destroyNKsolver can be used to destroy the objects created
  ! in this function

  use communication
  use constants
  use inputTimeSpectral
  use flowVarRefState
  use ADjointVars , only: nCellsLocal
  use NKSolverVars, only: snes,dRdw,dRdwPre,ctx,jacobian_lag,NKsolvedOnce, &
       snes_stol,snes_max_its,snes_max_funcs,nksolversetup,wVec,rVec,itertot0

  use InputIO ! L2conv,l2convrel
  use inputIteration
  use monitor
  use killSignals
  use iteration
  use stencils
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Working Variables
  integer(kind=intType) :: ierr,nDimw
  integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  real(kind=realType) :: rhoRes,rhoRes1,totalRRes
  integer(kind=intType) :: version
  integer(kind=intType) :: n_stencil,nCellTotal
  integer(kind=intType), dimension(:,:), allocatable :: stencil

  external FormFunction,FormJacobian,snes_monitor,LSCheck

  if (.not. NKsolverSetup) then

     !  Create nonlinear solver context
     call SNESCreate(SUMB_PETSC_COMM_WORLD,snes,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     !  Create residual and state vectors
     nDimW = nw * nCellsLocal * nTimeIntervalsSpectral
     call VecCreateMPI(SUMB_PETSC_COMM_WORLD,nDimw,PETSC_DETERMINE,wVec,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call VecSetBlockSize(wVec,nw,ierr);  call EChk(ierr,__FILE__,__LINE__)
     call VecDuplicate(wVec, rVec, ierr);  call EChk(ierr,__FILE__,__LINE__)
     
     !  Set Non-linear Function
     call SNESSetFunction(snes,rVec,FormFunction,ctx,ierr)
     call EChk(ierr,__FILE__,__LINE__)

   
     !  Create Jacobian and Approximate Jacobian Matrices
     call MatCreateSNESMF(snes,dRdw,ierr);  call EChk(ierr,__FILE__,__LINE__)
     
     ! Need to get correct Pre-allocation for dRdwPre; we can re-use
     ! adjoint routines for this
     allocate( nnzDiagonal(nDimW/nw),nnzOffDiag(nDimW/nw))

     call initialize_stencils
     if (not(viscous)) then
        n_stencil = N_euler_pc
        allocate(stencil(n_stencil,3))
        stencil = euler_pc_stencil
     else
        n_stencil = N_visc_pc
        allocate(stencil(n_stencil,3))
        stencil = visc_pc_stencil
     end if

     call statePreAllocation(nnzDiagonal,nnzOffDiag,nDimW/nw,stencil,n_stencil)
  
     call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,&
          nDimW, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dRdWPre, ierr); call EChk(ierr,__FILE__,__LINE__)
     
     deallocate(nnzDiagonal,nnzOffDiag,stencil)
     
#ifdef USE_PETSC_3
     call MatSetOption(dRdWPre, MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)
#else
     call MatSetOption(dRdWPre, MAT_COLUMN_ORIENTED, ierr)
     call EChk(ierr,__FILE__,__LINE__)
#endif
     
     !  Set Jacobian Function 
     call SNESSetJacobian(snes,dRdw,dRdwPre,FormJacobian,ctx,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     ! Set SNES Options
     ! Store the number of iterations completed by the RK solver
     iterTot0 = iterTot
     call SNESMonitorSet(snes,snes_monitor,ctx,PETSC_NULL_FUNCTION,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Use Eisenstat-Walker convergence criteria for KSP solver. Recommended
     call SNESKSPSetUseEW(snes,.True.,ierr)  
     call EChk(ierr,__FILE__,__LINE__)

     call SNESSetFromOptions(snes,ierr); call EChk(ierr,__FILE__,__LINE__)

     call SNESLineSearchSet(snes,SNESLineSearchNo,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set the Checking Function to use at the start of line search to
     ! make sure we dont have nans
!      print *,'Setting preCheck'
!      call SNESLineSearchSetPreCheck(snes,LSCheck,lsctx,ierr)
!      call EChk(ierr,__FILE__,__LINE__)

     ! See the monitor function for more information as to why this is -2
     call SNESSetLagJacobian(snes, -2_intType, ierr); call EChk(ierr,__FILE__,__LINE__)

     ! Since we're limiting the gmres to no restarts...there's a good
     ! chance that we're get lots of solve failues which is OK. Set
     ! this to the ncycles....basically large enough that it never happens
     call SNESSetMaxLinearSolveFailures(snes, ncycles,ierr); call EChk(ierr,__FILE__,__LINE__)
     
     ! We are going to have to compute what the tolerances should be
     ! since we are going to be using the same convergence criteria as
     ! SUmb originally uses, that is L2Conv and L2ConvRel. This however,
     ! gets a little trickier, since the NKsolver will always be called
     ! after the RK solver has been run at least once to get a good
     ! starting point. 

     NKSolverSetup = .True.
     NKSolvedOnce = .False.
  end if
end subroutine setupNKsolver_petsc


! subroutine setupNKsolverPseudo

!   ! Setup the PETSc objects for the Newton-Krylov
!   ! solver. destroyNKsolver can be used to destroy the objects created
!   ! in this function

!   use communication
!   use constants
!   use inputTimeSpectral
!   use flowVarRefState
!   use ADjointVars , only: nCellsLocal
!   use NKSolverVars, only: snes,dRdw,dRdwPre,ctx,jacobian_lag,NKsolvedOnce, &
!        snes_stol,snes_max_its,snes_max_funcs,nksolversetup,wVec,rVec,itertot0, &
!        lsctx,pts
!   use InputIO ! L2conv,l2convrel
!   use inputIteration
!   use monitor
!   use killSignals
!   use iteration
!   implicit none
! #define PETSC_AVOID_MPIF_H
! #include "include/finclude/petsc.h"
! #include "include/finclude/petscts.h"
!   ! Working Variables
!   integer(kind=intType) :: ierr,nDimw
!   integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
!   real(kind=realType) :: rhoRes,rhoRes1,totalRRes

!   external FormFunction2,FormFunction3,FormJacobian2,ts_monitor,LSCheck

!   if (.not. NKsolverSetup) then

!      ! Create Time-stepping context
!      call TSCreate(SUMB_PETSC_COMM_WORLD,pts,ierr)
!      call EChk(ierr,__FILE__,__LINE__)


!      ! Set Problem Type
!      call TSSetProblemType(pts,TS_NONLINEAR,ierr)
!      call EChk(ierr,__FILE__,__LINE__)
     
!      !  Create residual and state vectors
!      nDimW = nw * nCellsLocal * nTimeIntervalsSpectral
!      call VecCreateMPI(SUMB_PETSC_COMM_WORLD,nDimw,PETSC_DETERMINE,wVec,ierr)
!      call EChk(ierr,__FILE__,__LINE__)
     
!      call VecSetBlockSize(wVec,nw,ierr)
!      call EChk(ierr,__FILE__,__LINE__)

!      call VecDuplicate(wVec, rVec, ierr)
!      call EChk(ierr,__FILE__,__LINE__)
     
!      ! Tell TS where to compute solutions
!      call TSSetSolution(pts,wVec,ierr)
!      call EChk(ierr,__FILE__,__LINE__)

!      !  Set Non-linear Function
!      call TSSetRHSFunction(pts,formFunction3,ctx,ierr)
!      call EChk(ierr,__FILE__,__LINE__)

!      !  Create Jacobian and Approximate Jacobian Matrices
!     ! Setup Matrix-Free dRdw matrix
!      call MatCreateMFFD(sumb_comm_world,nDimW,nDimW,&
!           PETSC_DETERMINE,PETSC_DETERMINE,dRdw,ierr)
!      call EChk(ierr,__FILE__,__LINE__)

!      call MatSetOption(dRdW   , MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
!      call EChk(ierr,__FILE__,__LINE__)
     
!      ! Need to get correct Pre-allocation for dRdwPre; we can re-use
!      ! adjoint routines for this
!      allocate( nnzDiagonal(nCellsLocal*nTimeIntervalsSpectral),&
!           nnzOffDiag(nCellsLocal*nTimeIntervalsSpectral) )

!      call drdwPCPreAllocation(nnzDiagonal,nnzOffDiag,&
!           nCellsLocal*nTimeIntervalsSpectral)
 
!      call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,&
!           nDimW, nDimW,                     &
!           PETSC_DETERMINE, PETSC_DETERMINE, &
!           0, nnzDiagonal,         &
!           0, nnzOffDiag,            &
!           dRdWPre, ierr); call EChk(ierr,__FILE__,__LINE__)
     
!      deallocate(nnzDiagonal,nnzOffDiag)
     
!      call MatSetOption(dRdWPre, MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
!      call EChk(ierr,__FILE__,__LINE__)
     
!      !  Set Jacobian Function 
!      call TSSetRHSJacobian(pts,dRdw,dRdwPre,FormJacobian2,ctx,ierr)
!      call EChk(ierr,__FILE__,__LINE__)

!      ! Set the type to Pseudo Transient
!      call TSSetType(pts,TSPSEUDO,ierr)
!      call EChk(ierr,__FILE__,__LINE__)
     
!      ! Set initial time step:
!      call TSSetInitialTimeStep(pts,0.0,0.1,ierr)
!      call EChk(ierr,__FILE__,__LINE__)

!      ! Set Duration
!      call TSSetDuration(pts,1000,1.0e12,ierr)
!      call EChk(ierr,__FILE__,__LINE__)

!      ! Set options from database
!      call TSSetFromOptions(pts,ierr)
!      call EChk(ierr,__FILE__,__LINE__)

!      ! Setup TS object
!      call TSSetUp(pts,ierr)
!      call EChk(ierr,__FILE__,__LINE__)

!      ! Set TS Options
!      ! Store the number of iterations completed by the RK solver
!      iterTot0 = iterTot
!      call TSMonitorSet(pts,ts_monitor,ctx,PETSC_NULL_FUNCTION,ierr)
!      call EChk(ierr,__FILE__,__LINE__)

!      ! Set the pseudo time step function -- Using a default one here
!      !call TSPseudoSetTimeStep(pts,TSPseudoDefaultTimeStep,ctx,ierr)
!      !call EChk(ierr,__FILE__,__LINE__)

!      ! Pull out the SNES Context to set a few of those options:
!      !call TSGetSNES(pts,snes,ierr)
!      !call EChk(ierr,__FILE__,__LINE__)

!      ! Use Eisenstat-Walker convergence criteria for KSP solver. Recommended
!      !call SNESKSPSetUseEW(snes,.True.,ierr)  
!      !call EChk(ierr,__FILE__,__LINE__)

!      !call SNESSetFromOptions(snes,ierr)
!      !call EChk(ierr,__FILE__,__LINE__)

!      ! See the monitor function for more information as to why this is -2
!      !call SNESSetLagJacobian(snes, -2_intType, ierr); call EChk(ierr,__FILE__,__LINE__)

!      ! Since we're limiting the gmres to no restarts...there's a good
!      ! chance that we're get lots of solve failues which is OK. Set
!      ! this to the ncycles....basically large enough that it never happens
!      !call SNESSetMaxLinearSolveFailures(snes, ncycles,ierr); call EChk(ierr,__FILE__,__LINE__)
     
!      ! We are going to have to compute what the tolerances should be
!      ! since we are going to be using the same convergence criteria as
!      ! SUmb originally uses, that is L2Conv and L2ConvRel. This however,
!      ! gets a little trickier, since the NKsolver will always be called
!      ! after the RK solver has been run at least once to get a good
!      ! starting point. 

!      NKSolverSetup = .True.
!      NKSolvedOnce = .False.
!   end if
!   print *,'Done Pseudo Setup'
! end subroutine setupNKsolverPseudo
