subroutine FormJacobian_petsc(snes,wVec,dRdw,dRdwPre,flag,ctx,ierr)
  use communication
  use precision 
  use NKSolverVars,only : ksp_solver_type,ksp_subspace,global_pc_type,&
       asm_overlap,local_pc_ilu_level,local_pc_ordering,NKfinitedifferencepc!,rvec
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  SNES           snes
  Mat            dRdw,dRdwPre 
  KSP            ksp
  Vec            wVec,rVec
  PetscFortranAddr ctx(3)
  MatStructure   flag 
  PetscInt nlocal,first,Nsub,length
  integer(kind=intType) ::ierr

  ! Local Variables
  logical secondHalo
  logical :: useAD,usePC,useTranspose
 
  ! Dummy assembly begin/end calls for the matrix-free Matrx
  call MatAssemblyBegin(dRdw,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call MatAssemblyEnd(dRdw,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Assemble the approximate PC
  useAD = .False.
  usePC = .True.
  useTranspose = .False.
  call setupStateResidualMatrix(dRdwPre,useAD,usePC,useTranspose)

  flag = SAME_NONZERO_PATTERN
  ! Setup the required options for the KSP solver
  call SNESGetKSP(snes,ksp,ierr);  
  call EChk(ierr,__FILE__,__LINE__)

  call NKSetup_KSP(ksp)

end subroutine FormJacobian_petsc


















subroutine FormJacobian2(pts,t,wVec,dRdw,dRdwPre,flag,ctx,ierr)
  ! This is a wrapper for the formJacobian function. The formJacobian
  ! for the time stepping routine is slightly different. Basically it
  ! calls with the TimeStepping Context instead of SNES. We are just
  ! going to pull ou the SNES context and call the formJacobian
  ! function which will do all the required calculations/option
  ! setting.
  use NKSolverVars, only: snes_stol,snes_max_its,snes_max_funcs,snes_rtol,&
       snes_atol

  use precision 
  use inputIteration
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
#include "include/finclude/petscts.h"
  TS             pts
  SNES           snes
  Mat            dRdw,dRdwPre 
  real(kind=realType) :: t
  Vec            wVec
  PetscFortranAddr ctx(3)
  MatStructure   flag 
  integer(kind=intType) ::ierr
  external formfunction2

  call TSGetSNES(pts,snes,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MatMFFDSetFunction(dRdw,FormFunction2,ctx,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  call MatSetOption(dRdW   , MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MatMFFDSetBase(dRdW,wVec,PETSC_NULL_OBJECT,ierr)
  call EChk(ierr,__FILE__,__LINE__)


  !Use Eisenstat-Walker convergence criteria for KSP solver. Recommended
  call SNESKSPSetUseEW(snes,.True.,ierr)  
  call EChk(ierr,__FILE__,__LINE__)

  ! See the monitor function for more information as to why this is -2
  call SNESSetLagJacobian(snes, -2_intType, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call SNESSetTolerances(snes,snes_atol,snes_rtol,snes_stol,1,&
       100,ierr); call EChk(ierr,__FILE__,__LINE__)

  ! Since we're limiting the gmres to no restarts...there's a good
  ! chance that we're get lots of solve failues which is OK. Set
  ! this to the ncycles....basically large enough that it never happens
  call SNESSetMaxLinearSolveFailures(snes, ncycles,ierr)
  call EChk(ierr,__FILE__,__LINE__)
     
  

  call FormJacobian(snes,wVec,dRdw,dRdwPre,flag,ctx,ierr)
  call EChk(ierr,__FILE__,__LINE__)
end subroutine FormJacobian2









! subroutine setupNK_KSP(ksp)

! ! This is a shell-type routine to setup the KSP object for the newton
! ! Krylov solver. Since there are a number of formJacobian functions,
! ! this eliminates repetition. 

!   KSP ksp

!   use NKSolverVars,only : ksp_solver_type,ksp_subspace,global_pc_type,&
!        asm_overlap,local_pc_ilu_level,local_pc_ordering,NKfinitedifferencepc
!   implicit none
! #define PETSC_AVOID_MPIF_H
! #include "include/finclude/petsc.h"

!   call KSPSetType(ksp,ksp_solver_type,ierr);
!   call EChk(ierr,__FILE__,__LINE__)
!   call KSPGMRESSetRestart(ksp, ksp_subspace,ierr); 
!   call EChk(ierr,__FILE__,__LINE__)
!   call KSPSetPreconditionerSide(ksp,PC_RIGHT,ierr);
!   call EChk(ierr,__FILE__,__LINE__)

!   ! Setup the required options for the Global PC
!   call KSPGetPC(ksp,pc,ierr);             
!   call EChk(ierr,__FILE__,__LINE__)

!   ! Set the global PC Type:
!   call PCSetType(pc,global_pc_type,ierr)
!   call EChk(ierr,__FILE__,__LINE__)

!   if (trim(global_pc_type) == 'asm') then
!      ! Set overlap
!      call PCASMSetOverlap(pc,asm_overlap,ierr)
!      call EChk(ierr,__FILE__,__LINE__)

!      ! Setup the global PC so we can extract subKSP objects
!      call PCSetup(pc,ierr)
!      call EChk(ierr,__FILE__,__LINE__)
     
!      call PCASMGetSubKSP( pc, nlocal,  first, subksp, ierr )
!      call EChk(ierr,__FILE__,__LINE__)  

!   else  if (trim(global_pc_type) == 'bjacobi') then
!      ! Setup the global PC so we can extract subKSP objects
!      call PCSetup(pc,ierr)
!      call EChk(ierr,__FILE__,__LINE__)
     
!      call PCBJacobiGetSubKSP(pc,nlocal,first,subksp,ierr)
!      call EChk(ierr,__FILE__,__LINE__)
!   end if

!   ! Additional setup for subpc when block jacobi or asm are used:
!   if (trim(global_pc_type) == 'bjacobi' .or. &
!       trim(global_pc_type) == 'asm') then

!      ! Setup the required options for the Local PC
!      call KSPGetPC(subksp, subpc, ierr )
!      call EChk(ierr,__FILE__,__LINE__)

!      call PCSetType(subpc, 'ilu', ierr)
!      call EChk(ierr,__FILE__,__LINE__)
!      call PCFactorSetLevels(subpc, local_pc_ilu_level, ierr)
!      call EChk(ierr,__FILE__,__LINE__)  

!      call PCFactorSetMatOrderingtype(subpc, local_pc_ordering, ierr )
!      call EChk(ierr,__FILE__,__LINE__) 
!      call KSPSetType(subksp, KSPPREONLY, ierr)
!      call EChk(ierr,__FILE__,__LINE__)  

!   end if

!   ierr = 0

! end subroutine setupNK_KSP
