subroutine setupNKsolver
#ifndef USE_NO_PETSC
  ! Setup the PETSc objects for the Newton-Krylov
  ! solver. destroyNKsolver can be used to destroy the objects created
  ! in this function
  use blockPointers
  use communication
  use inputTimeSpectral
  use flowVarRefState
  use iteration
  use inputPhysics
  use stencils
  use ADjointVars , only: nCellsLocal
  use NKSolverVars, only: dRdw,dRdwPre,dRdwPseudo, ctx, wVec,rVec,deltaW,&
       NKsolvecount,nksolversetup,ksp_subspace,ksp_solver_type,global_ksp

  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Working Variables
  integer(kind=intType) :: ierr,nDimw,totalCells
  integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  integer(kind=intType) :: n_stencil
  integer(kind=intType), dimension(:,:), allocatable :: stencil

  integer(kind=intType) :: i,j,k,nn,i2,j2,k2,d2,l,sps
  external FormFunction_mf, mykspmonitor

  if (not(NKSolverSetup)) then
     nDimW = nw * nCellsLocal * nTimeIntervalsSpectral

     call VecCreateMPI(SUMB_PETSC_COMM_WORLD,nDimw,PETSC_DETERMINE,wVec,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call VecSetBlockSize(wVec,nw,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     !  Create residual and state vectors
     call VecCreateMPI(SUMB_PETSC_COMM_WORLD,nDimw,PETSC_DETERMINE,rVec,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call VecSetBlockSize(rVec,nw,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Use the rVec Template to create deltaW 
     call VecDuplicate(rVec, deltaW, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Create Pre-Conditioning Matrix
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

     ! Uncomment if you want to look at the KSP residual during the
     ! solve...this is a hack for debuggin only...it uses the adjoint
     ! monitor function 

     ! call KSPMonitorSet(global_ksp,MyKSPMonitor, PETSC_NULL_OBJECT, &
     !      PETSC_NULL_FUNCTION, ierr) call EChk(ierr,__FILE__,__LINE__)

     NKSolverSetup = .True.
     NKSolveCount = 0
     if(equations == RANSEquations)  then
        turbCoupled = .True.
        turbSegregated = .False.
     end if
  end if
#endif
end subroutine setupNKsolver
