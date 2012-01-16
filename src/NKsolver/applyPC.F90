subroutine applyPC(in_vec,out_vec,N)
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
  external FormFunction2

  ! Working Variables
  integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  integer(kind=intTYpe) :: ierr,ndimw,nlocal,first,i,ilow,ihigh,size
  real(kind=realType) :: value
  integer(kind=intType) :: blksize
  logical :: useAD,usePC,useTranspose
  
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

     allocate( nnzDiagonal(nCellsLocal*nTimeIntervalsSpectral),&
          nnzOffDiag(nCellsLocal*nTimeIntervalsSpectral) )

     call drdwPCPreAllocation(nnzDiagonal,nnzOffDiag,nCellsLocal*nTimeIntervalsSpectral)
     call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
          nDimW, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dRdWPre, ierr); call EChk(ierr,__FILE__,__LINE__)
     
     deallocate(nnzDiagonal,nnzOffDiag)

     ! Setup Matrix-Free dRdw matrix
     call MatCreateMFFD(sumb_comm_world,nDimW,nDimW,&
          PETSC_DETERMINE,PETSC_DETERMINE,dRdw,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call MatMFFDSetFunction(dRdw,FormFunction2,ctx,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call MatAssemblyBegin(dRdw,MAT_FINAL_ASSEMBLY,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call MatAssemblyEnd(dRdw,MAT_FINAL_ASSEMBLY,ierr)
     call EChk(ierr,__FILE__,__LINE__)
             
     call MatSetOption(dRdWPre, MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call MatSetOption(dRdW   , MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     if (NKFiniteDifferencePC) then
        useAD = .False.
        usePC = .True.
        useTranspose = .False.
        call setupStateResidualMatrix(dRdwPre,useAD,usePC,useTranspose)
     else
        call setupNK_PC(dRdwPre)
     end if
     
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

     ! Setup the required options for the Global PC
     call KSPSetup(global_ksp,ierr);      call EChk(ierr,__FILE__,__LINE__)
     call KSPGetPC(global_ksp,global_pc,ierr);                 call EChk(ierr,__FILE__,__LINE__)
     call PCSetType(global_pc,global_pc_type,ierr);     call EChk(ierr,__FILE__,__LINE__)
 
     call KSPSetTolerances(global_ksp,L2ConvRel,1e-10,1e5,ksp_subspace,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     if (trim(global_pc_type) == 'asm') then
        call PCASMSetOverlap(global_pc,asm_overlap,ierr);  call EChk(ierr,__FILE__,__LINE__)
        call PCSetup(global_pc,ierr);                      call EChk(ierr,__FILE__,__LINE__)
        call PCASMGetSubKSP(global_pc, nlocal,  first, local_ksp, ierr );          call EChk(ierr,__FILE__,__LINE__)  
     end if

     if (trim(global_pc_type) == 'bjacobi') then
        call PCSetup(global_pc,ierr);                      call EChk(ierr,__FILE__,__LINE__)
        call PCBJacobiGetSubKSP(global_pc,nlocal,first,local_ksp,ierr);   call EChk(ierr,__FILE__,__LINE__)
     end if

     ! Setup the required options for the Local PC
     call KSPGetPC(local_ksp, local_pc, ierr );                              call EChk(ierr,__FILE__,__LINE__)
     call PCSetType(local_pc, 'ilu', ierr);                       call EChk(ierr,__FILE__,__LINE__)
     call PCFactorSetLevels(local_pc, local_pc_ilu_level, ierr);          call EChk(ierr,__FILE__,__LINE__)  
     call PCFactorSetMatOrderingtype(local_pc, local_pc_ordering, ierr ); call EChk(ierr,__FILE__,__LINE__) 
     call KSPSetType(local_ksp, KSPPREONLY, ierr);         call EChk(ierr,__FILE__,__LINE__)  

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

end subroutine applyPC
