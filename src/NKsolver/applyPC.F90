subroutine applyPC(in_vec,out_vec,N)
#ifndef USE_NO_PETSC
  ! Apply the NK PC to the in_vec. This subroutine is ONLY used as a
  ! preconditioner for a global Aero-Structural Newton-Krylov Method
  use inputTimeSpectral
  use flowVarRefState
  use ADjointVars , only: nCellsLocal
  use NKSolverVars, only: dRdw,dRdwPre,NKPCSetup,&
       global_ksp,local_ksp,global_pc,local_pc, &
       ksp_rtol,ksp_atol,ksp_div_tol,ksp_max_it,ctx,&
       ksp_subspace,global_pc_type,asm_overlap,local_pc_ilu_level,&
       local_pc_ordering,nkfinitedifferencepc
  use communication 
  use inputIteration
  use stencils
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Input/Output
  integer(kind=intType) :: N
  real(kind=realType), dimension(N), intent(in) :: in_vec(N)
  real(kind=realTYpe), dimension(N), intent(out):: out_vec(N)
  
  ! PETSc 
  Vec VecA, VecB,wVec,rVec
  KSPConvergedReason ksp_reason

  external KSPSkipConverged
  external FormFunction_mf

  ! Working Variables
  integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  integer(kind=intTYpe) :: ierr,ndimw,nlocal,first,i,ilow,ihigh,size
  real(kind=realType) :: value
  integer(kind=intType) :: blksize
  logical :: useAD,usePC,useTranspose
  integer(kind=intType) :: n_stencil,totalCells
  integer(kind=intType), dimension(:,:), allocatable :: stencil

  nDimW = nw * nCellsLocal * nTimeIntervalsSpectral
  
  ! Put a petsc wrapper around the input and output vectors
  call VecCreateMPIWithArray(sumb_comm_world,N,PETSC_DETERMINE,in_vec,VecA,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecSetBlockSize(vecA,nw,ierr);
  call EChk(ierr,__FILE__,__LINE__)

  call VecCreateMPIWithArray(sumb_comm_world,N,PETSC_DETERMINE,out_vec,VecB,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecSetBlockSize(vecB,nw,ierr);
  call EChk(ierr,__FILE__,__LINE__)

  if (not(NKPCSetup)) then

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

     ! Setup Matrix-Free dRdw matrix
     call MatCreateMFFD(sumb_comm_world,nDimW,nDimW,&
          PETSC_DETERMINE,PETSC_DETERMINE,dRdw,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call MatMFFDSetFunction(dRdw,FormFunction_mf,ctx,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call MatAssemblyBegin(dRdw,MAT_FINAL_ASSEMBLY,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call MatAssemblyEnd(dRdw,MAT_FINAL_ASSEMBLY,ierr)
     call EChk(ierr,__FILE__,__LINE__)
             
     call MatSetOption(dRdWPre, MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call MatSetOption(dRdW   , MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     useAD = .False.
     usePC = .True.
     useTranspose = .False.
     call setupStateResidualMatrix(dRdwPre,useAD,usePC,useTranspose)
     
     call KSPCreate(SUMB_COMM_WORLD, global_ksp, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call KSPSetOperators(global_ksp,dRdW,dRdWPre,DIFFERENT_NONZERO_PATTERN,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call KSPSetFromOptions(global_ksp,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call KSPGMRESSetRestart(global_ksp,ksp_subspace, ierr)
     call EChk(ierr ,__FILE__,__LINE__)

     call KSPSetPCSide(global_ksp, PC_RIGHT, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call NKSetup_KSP(global_ksp)

     NKPCSetup = .True.
  end if

  ! Set the base vec
  call VecDuplicate(vecA,wVec,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  call setwVec(wVec)
  
  call MatMFFDSetBase(dRdW,wVec,PETSC_NULL_OBJECT,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Actually do the Linear Krylov Solve
  call KSPSolve(global_ksp,vecA,vecB,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Destroy the three petsc vectors
   call VecDestroy(wVec,ierr);  call EChk(ierr,__FILE__,__LINE__)
   call VecDestroy(VecA,ierr);  call EChk(ierr,__FILE__,__LINE__)
   call VecDestroy(VecB,ierr);  call EChk(ierr,__FILE__,__LINE__)

#endif
end subroutine applyPC
