subroutine setupNKsolver

  ! Setup the PETSc objects for the Newton-Krylov
  ! solver. destroyNKsolver can be used to destroy the objects created
  ! in this function

  use communication
  use inputTimeSpectral
  use flowVarRefState
  use iteration
  use inputPhysics
  use stencils
  use ADjointVars , only: nCellsLocal
  use NKSolverVars, only: dRdw,dRdwPre,dRdwPseudo, ctx, wVec,rVec,deltaW,&
       NKsolvedOnce,nksolversetup,ksp_subspace,ksp_solver_type,global_ksp

  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Working Variables
  integer(kind=intType) :: ierr,nDimw,totalCells
  integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  integer(kind=intType) :: n_stencil
  integer(kind=intType), dimension(:,:), allocatable :: stencil

  external FormFunction_mf

  if (not(NKSolverSetup)) then
     
     !  Create residual and state vectors
     nDimW = nw * nCellsLocal * nTimeIntervalsSpectral
     call VecCreateMPI(SUMB_PETSC_COMM_WORLD,nDimw,PETSC_DETERMINE,wVec,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call VecSetBlockSize(wVec,nw,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Use the wVec Template to create deltaW and rVec
     call VecDuplicate(wVec, rVec, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDuplicate(wVec, deltaW, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Setup Pre-Conditioning Matrix
     totalCells = nCellsLocal*nTimeIntervalsSpectral
     allocate( nnzDiagonal(totalCells),nnzOffDiag(totalCells))

     call initialize_stencils
     if (not(viscous)) then
        n_stencil = N_euler_drdw
        allocate(stencil(n_stencil,3))
        stencil = euler_drdw_stencil
     else
        n_stencil = N_visc_pc
        allocate(stencil(n_stencil,3))
        stencil = visc_pc_stencil
     end if

     call statePreAllocation(nnzDiagonal,nnzOffDiag,nDimW/nw,stencil,n_stencil)
  
     call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
          nDimW, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dRdWPre, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     deallocate(nnzDiagonal,nnzOffDiag,stencil)

     ! Setup Matrix-Free dRdw matrix and its function
     call MatCreateMFFD(sumb_comm_world,nDimW,nDimW,&
          PETSC_DETERMINE,PETSC_DETERMINE,dRdw,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call MatMFFDSetFunction(dRdw,FormFunction_mf,ctx,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set the mat_row_oriented option to false so that dense
     ! subblocks can be passed in in fortran column-oriented format
     call MatSetOption(dRdWPre, MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call MatSetOption(dRdW   , MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     !  Create the linear solver context
     call KSPCreate(SUMB_PETSC_COMM_WORLD,global_ksp,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set operators for the solver
     call KSPSetOperators(global_ksp,dRdw,dRdWPre, DIFFERENT_NONZERO_PATTERN,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call SNESSetFromOptions(snes,ierr); call EChk(ierr,__FILE__,__LINE__)

     !call SNESLineSearchSet(snes,SNESLineSearchNo,ierr)
     !call EChk(ierr,__FILE__,__LINE__)

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
     if(equations == RANSEquations)  then
        turbCoupled = .True.
        turbSegregated = .False.
     end if
  end if
end subroutine setupNKsolver


subroutine MyMult(matrix,X,F,ierr)

  !   Input Parameters:
  !.  X - input vector
  !
  !   Output Parameter:
  !.  F - function vector
  !
  use precision
  use NKSolverVars, only: dRdw,diagV
  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
  
  ! Input/Output Vars
  Mat matrix
  Vec X,F
  integer(kind=intType) :: ierr

  ! Do a matmult followed by an addition

  call MatMult(dRdw,X,F,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! VecAXPY : Computes y = alpha x + y. 
  ! VecAXPY(Vec y,PetscScalar alpha,Vec x)

  !call VecAXPY(F,1.0,diagV,ierr)
  !call EChk(ierr,__FILE__,__LINE__)

end subroutine MyMult

