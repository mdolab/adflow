!
!      ******************************************************************
!      *                                                                *
!      * File:          NKsolver.F90                                    *
!      * Author:        Gaetan Kenway                                   *
!      * Starting date: 11-27-2010                                      *
!      * Last modified: 11-27-2010                                      *
!      *                                                                *
!      ******************************************************************

!
!      ******************************************************************
!      *                                                                *
!      * File:          NKsolver.F90                                    *
!      * Author:        Gaetan Kenway                                   *
!      * Starting date: 11-27-2010                                      *
!      * Last modified: 11-27-2010                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine setupNKsolver

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
  implicit none
#define PETSC_AVOID_MPIF_H
#if PETSC_VERSION_MINOR>=1
#include "include/finclude/petsc.h"
#else
#include "include/finclude/petscall.h"
#endif

  ! Working Variables
  integer(kind=intType) :: ierr,nDimw
  integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  real(kind=realType) :: rhoRes,rhoRes1,totalRRes

  external FormFunction,FormJacobian,snes_monitor

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
     allocate( nnzDiagonal(nCellsLocal*nTimeIntervalsSpectral),&
          nnzOffDiag(nCellsLocal*nTimeIntervalsSpectral) )
     
     call drdwPCPreAllocation(nnzDiagonal,nnzOffDiag,nCellsLocal)
     call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
          nDimW, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dRdWPre, ierr); call EChk(ierr,__FILE__,__LINE__)
     
     deallocate(nnzDiagonal,nnzOffDiag)
     
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
     !call SNESSetLagJacobian(snes, jacobian_lag, ierr); call EChk(ierr,__FILE__,__LINE__)
     call SNESSetLagJacobian(snes, -2_intType, ierr); call EChk(ierr,__FILE__,__LINE__)
     call SNESSetMaxLinearSolveFailures(snes, 25,ierr); call EChk(ierr,__FILE__,__LINE__)
     
     ! We are going to have to compute what the tolerances should be
     ! since we are going to be using the same convergence criteria as
     ! SUmb originally uses, that is L2Conv and L2ConvRel. This however,
     ! gets a little trickier, since the NKsolver will always be called
     ! after the RK solver has been run at least once to get a good
     ! starting point. 

     snes_stol = 1e-10
     snes_max_its = 1250_intType
     snes_max_funcs = snes_max_its * 2
     
     NKSolverSetup = .True.
     NKSolvedOnce = .False.
  end if
end subroutine setupNKsolver

subroutine destroyNKsolver

  use NKsolverVars
  implicit none
  integer(kind=intType) :: ierr
  ! We will destroy the PETSc variables created in setupNKsolver
  call SNESDestroy(snes,ierr) ! Also destroys the underlying ksp and
                              ! pc contexts
  call EChk(ierr,__FILE__,__LINE__)
  call MatDestroy(dRdw,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call MatDestroy(dRdwPre,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecDestroy(wVec,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecDestroy(rVec,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  NKSolverSetup = .False.

end subroutine DestroyNKsolver

subroutine NKsolver
  use communication
  use constants
  use inputTimeSpectral
  use flowVarRefState
  use ADjointVars , only: nCellsLocal
  use NKSolverVars, only: snes,dRdw,dRdwPre,ctx,jacobian_lag,&
       snes_stol,snes_max_funcs,nksolversetup,rhoRes0,&
       snes_rtol,snes_atol,totalRes0,itertot0,wVec,rVec,reason,NKSolvedOnce
  use InputIO ! L2conv,l2convrel
  use inputIteration
  use monitor
  use killSignals
  use iteration
  implicit none
#define PETSC_AVOID_MPIF_H
#if PETSC_VERSION_MINOR>=1
#include "include/finclude/petsc.h"
#else
#include "include/finclude/petscall.h"
#endif

  integer(kind=intTYpe) :: sns_max_its,ierr,snes_max_its,temp
  real(kind=realType) :: rhoRes,totalRRes,rhoRes1


 ! We are going to have to compute what the tolerances should be
  ! since we are going to be using the same convergence criteria as
  ! SUmb originally uses, that is L2Conv and L2ConvRel. This however,
  ! gets a little trickier, since the NKsolver will always be called
  ! after the RK solver has been run at least once to get a good
  ! starting point. 

  snes_stol = 1e-10
  snes_max_its = 1250_intType
  snes_max_funcs = snes_max_its * 2

  ! Determine the current level of convergence of the solution

  call getCurrentResidual(rhoRes,totalRRes)

  if (myid == 0) then

     rhoRes1 = convArray(1,1,1) ! Second Density residual for sps 1 (startResidual)

     
     if (rhoRes1 == 0) then ! This will happen on subsequent solves
        rhoRes1 = rhoRes0
     end if

     ! We need to compute two convergences: One coorsponding to L2ConvRel
     ! and one for L2Conv
     
     snes_rtol = (rhoRes1 * L2ConvRel)/ rhoRes  ! Target / Current
  
     ! Absolute Tol is the original totalR * L2conv

     snes_atol = totalRes0 * L2Conv
  end if

  ! BroadCast
  call mpi_bcast(snes_rtol, 1, sumb_real, 0, SUmb_comm_world, ierr)
  call mpi_bcast(snes_atol, 1, sumb_real, 0, SUmb_comm_world, ierr)

  call SNESSetTolerances(snes,snes_atol,snes_rtol,snes_stol,snes_max_its,&
       snes_max_funcs,ierr); call EChk(ierr,__FILE__,__LINE__)

  ! Note: the krylov linear solver options are set in FormJacobian

  ! Form the initial guess from the current w-vector
  call setwVec(wVec)

  ! Solve IT!
  call SNESSolve(snes,PETSC_NULL_OBJECT,wVec,ierr) ! PETSC_NULL_OBJECT
                                                   ! MAY GIVE MEMORY
                                                   ! LEAK!!!!!!
  
  NKSolvedOnce = .True.
  call EChk(ierr,__FILE__,__LINE__)
  iterTot = iterTot0

  call SNESGetConvergedReason(snes,reason,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  if (myid == 0) then
     if (printIterations) then
        print *,'Reason:',reason
     end if
  end if
  if (reason == SNES_CONVERGED_FNORM_ABS .or. &
      reason == SNES_CONVERGED_FNORM_RELATIVE .or. &
      reason == SNES_CONVERGED_PNORM_RELATIVE) then
     routineFailed = .False. 
  else
     routineFailed = .True. 
  end if


end subroutine NKsolver

subroutine FormFunction(snes,wVec,rVec,ctx,ierr)
  ! ---------------------------------------------------------------------
  !
  !  FormFunction - Evaluates nonlinear function, f(x).
  !
  !  Input Parameters:
  !  snes  - the SNES context
  !  wVec  - input vector
  !
  !  Currnet ctx( ) doesn't have anythign in it
  !
  !  Output Parameter:
  !  f     - vector with newly computed function
  use precision
  implicit none
#define PETSC_AVOID_MPIF_H
#if PETSC_VERSION_MINOR>=1
#include "include/finclude/petsc.h"
#else
#include "include/finclude/petscall.h"
#endif

  ! PETSc Variables
  SNES    snes
  Vec     wVec, rVec
  PetscFortranAddr ctx(*)
  integer(kind=intType) :: ierr

  ! This is just a shell routine that runs the more broadly useful
  ! computeResidualNK subroutine
  call setW(wVec)
  call computeResidualNK()
  call setRVec(rVec)

end subroutine FormFunction

subroutine computeResidualNK()
  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate
  use iteration
  use inputPhysics 
  implicit none
  ! Local Variables
  integer(kind=intType) :: ierr,i,j,k,l,sps,nn
  logical secondHalo ,correctForK
  real(kind=realType) :: gm1,v2,val

  ! Next we Run the equilivent of the residual routine
  secondHalo = .True. 

  currentLevel = 1_intType
  groundLevel = 1_intTYpe
  ! Next we need to compute the pressures
  gm1 = gammaConstant - one
  spectralLoop: do sps=1,nTimeIntervalsSpectral
     domainsState: do nn=1,nDom
        ! Set the pointers to this block.
        call setPointers(nn, currentLevel, sps)

        do k=2,kl 
           do j=2,jl
              do i=2,il

                 v2 = w(i,j,k,ivx)**2 + w(i,j,k,ivy)**2 &
                      + w(i,j,k,ivz)**2

                 p(i,j,k) = gm1*(w(i,j,k,irhoE) &
                      - half*w(i,j,k,irho)*v2)
                 p(i,j,k) = max(p(i,j,k), 1.e-4_realType*pInfCorr)
              enddo
           enddo
        enddo
     end do domainsState
  end do spectralLoop

  !   Apply BCs
  call applyAllBC(secondHalo)

  ! Exchange solution -- always the fine level
   call whalo2(1_intType, 1_intType, nMGVar, .true., &
       .true., .true.)

   ! Why does this need to be set?
  rkStage = 0
 
  call timestep(.false.)
  call initres(1_intType, nwf)
  call residual 
  
end subroutine computeResidualNK

subroutine FormJacobian(snes,wVec,dRdw,dRdwPre,flag,ctx,ierr)
  use communication
  use precision 
  use NKSolverVars,only : ksp_solver_type,ksp_subspace,global_pc_type,&
       asm_overlap,local_pc_ilu_level,local_pc_ordering
  implicit none
#define PETSC_AVOID_MPIF_H
#if PETSC_VERSION_MINOR>=1
#include "include/finclude/petsc.h"
#else
#include "include/finclude/petscall.h"
#endif
  SNES           snes
  Mat            dRdw,dRdwPre 
  KSP            ksp,subksp
  PC             pc,subpc
  Vec            wVec,rVec
  PetscFortranAddr ctx(3)
  MatStructure   flag 
  PetscInt nlocal,first,Nsub,length
  integer(kind=intType) ::ierr

  ! Local Variables
  logical secondHalo

  ! Dummy assembly begin/end calls for the matrix-free Matrx
  call MatAssemblyBegin(dRdw,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call MatAssemblyEnd(dRdw,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Assemble the approximate PC
  call setupNK_KSP_PC(dRdwPre)
  flag = SAME_NONZERO_PATTERN
  ! Setup the required options for the KSP solver
  call SNESGetKSP(snes,ksp,ierr);                   call EChk(ierr,__FILE__,__LINE__)
  call KSPSetType(ksp,ksp_solver_type,ierr);       call EChk(ierr,__FILE__,__LINE__)
  call KSPGMRESSetRestart(ksp, ksp_subspace,ierr);  call EChk(ierr,__FILE__,__LINE__)
  call KSPSetPreconditionerSide(ksp,PC_RIGHT,ierr); call EChk(ierr,__FILE__,__LINE__)

  ! Setup the required options for the Global PC
  call KSPGetPC(ksp,pc,ierr);                 call EChk(ierr,__FILE__,__LINE__)
  call PCSetType(pc,global_pc_type,ierr);     call EChk(ierr,__FILE__,__LINE__)

  if (trim(global_pc_type) == 'asm') then
     call PCASMSetOverlap(pc,asm_overlap,ierr);  call EChk(ierr,__FILE__,__LINE__)
     call PCSetup(pc,ierr);                      call EChk(ierr,__FILE__,__LINE__)
     call PCASMGetSubKSP( pc, nlocal,  first, subksp, ierr );          call EChk(ierr,__FILE__,__LINE__)  
  end if

  if (trim(global_pc_type) == 'bjacobi') then
     call PCSetup(pc,ierr);                      call EChk(ierr,__FILE__,__LINE__)
     call PCBJacobiGetSubKSP(pc,nlocal,first,subksp,ierr);   call EChk(ierr,__FILE__,__LINE__)
  end if


  ! Setup the required options for the Local PC
  call KSPGetPC(subksp, subpc, ierr );                              call EChk(ierr,__FILE__,__LINE__)
  call PCSetType(subpc, 'ilu', ierr);                       call EChk(ierr,__FILE__,__LINE__)
  call PCFactorSetLevels(subpc, local_pc_ilu_level, ierr);          call EChk(ierr,__FILE__,__LINE__)  
  call PCFactorSetMatOrderingtype(subpc, local_pc_ordering, ierr ); call EChk(ierr,__FILE__,__LINE__) 
  call KSPSetType(subksp, KSPPREONLY, ierr);    
  call EChk(ierr,__FILE__,__LINE__)  

  ierr = 0
  
end subroutine FormJacobian

subroutine setWVec(wVec)

  ! Set the petsc vector wVec from the current SUmb soltuion in w

  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate 
  implicit none
#define PETSC_AVOID_MPIF_H
#if PETSC_VERSION_MINOR>=1
#include "include/finclude/petsc.h"
#else
#include "include/finclude/petscall.h"
#endif
  
  Vec     wVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l
  real(kind=realType) :: states(nw)

  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointersAdj(nn,1_intType,sps)
        ! Copy off w to wVec
        do k=2,kl
           do j=2,jl
              do i=2,il
                 do l=1,nw
                    states(l) = w(i,j,k,l)
                 end do
                 
                 call VecSetValuesBlocked(wVec,1,globalCell(i,j,k),states,&
                      INSERT_VALUES,ierr)
                 call EChk(ierr,__FILE__,__LINE__)
              end do
           end do
        end do
     end do
  end do
  call VecAssemblyBegin(wVec,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecAssemblyEnd(wVec,ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine setWVec

subroutine setRVec(rVec)
  ! Set the current residual in dw into the PETSc Vector

  use communication
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration
  implicit none
#define PETSC_AVOID_MPIF_H
#if PETSC_VERSION_MINOR>=1
#include "include/finclude/petsc.h"
#else
#include "include/finclude/petscall.h"
#endif

  Vec     rVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k
  real(kind=realType) :: ovv
  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointersAdj(nn,1_intType,sps)
        ! Copy off dw/vol to rVec
        do k=2,kl
           do j=2,jl
              do i=2,il
                 ovv = 1/vol(i,j,k)
                 call VecSetValuesBlocked(rVec,1,globalCell(i,j,k),&
                      dw(i,j,k,:)*ovv, INSERT_VALUES,ierr)
                 call EChk(ierr,__FILE__,__LINE__)
              end do
           end do
        end do
     end do
  end do

  call VecAssemblybegin(rVec,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecAssemblyEnd(rVec,ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine setRVec

subroutine setW(wVec)

  ! Set the SUmb state vector, w, from the petsc vec wVec

  use blockPointers
  use inputTimeSpectral
  use flowVarRefState
  implicit none

!#define PETSC_AVOID_MPIF_H
!#include "include/finclude/petsc.h"

  Vec     wVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l

  ! Note this is not ideal memory access but the values are stored by
  ! block in PETSc (grouped in nw) but are stored separately in SUmb
  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointersAdj(nn,1_intType,sps)
        ! Copy off w to wVec
        do k=2,kl
           do j=2,jl
              do i=2,il
                 do l=1,nw
                    call VecGetValues(wVec,1,globalCell(i,j,k)*nw+l-1,&
                         w(i,j,k,l),ierr)
                    call EChk(ierr,__FILE__,__LINE__)
                 
                 end do
              end do
           end do
        end do
     end do
  end do
end subroutine setW

subroutine setupNK_KSP_PC(dRdwPre)

  !     ******************************************************************
  !     *                                                                *
  !     * Compute the dRdWPre matrix for the NK Solver                   *
  !     ******************************************************************
  !
  use ADjointVars
  use blockPointers       ! i/j/kl/b/e, i/j/k/Min/MaxBoundaryStencil
  use communication       ! procHalo(currentLevel)%nProcSend
  use inputDiscretization ! spaceDiscr
  USE inputTimeSpectral   ! nTimeIntervalsSpectral
  use iteration           ! overset, currentLevel
  use flowVarRefState     ! nw
  use inputTimeSpectral   ! spaceDiscr
  use inputADjoint        !sigma
  implicit none
#define PETSC_AVOID_MPIF_H
#if PETSC_VERSION_MINOR>=1
#include "include/finclude/petsc.h"
#else
#include "include/finclude/petscall.h"
#endif

  Mat dRdwPre

 !
  !     Local variables.

  integer(kind=intType) :: iCell, jCell, kCell, nn, level, m
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw, nTimeIntervalsSpectral) :: wAdj,wadjb
  real(kind=realType), dimension(-3:2,-3:2,-3:2,3,&
       nTimeIntervalsSpectral) :: siAdj, sjAdj, skAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,&
       nTimeIntervalsSpectral) ::sFaceIAdj,sFaceJAdj,sFaceKAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,3,&
       nTimeIntervalsSpectral) :: sAdj
  real(kind=realType),dimension(nTimeIntervalsSpectral) :: volAdj  
  real(kind=realType),dimension(3) :: rotRateAdj
  real(kind=realType), dimension(nw,nTimeIntervalsSpectral)  :: dwAdj,dwadjb
  real(kind=realType), dimension(2) :: time
  real(kind=realType)               ::setupTime

  logical :: correctForK, secondHalo, exchangeTurb

  ! dR/dw stencil

  real(kind=realType), dimension(nw,nw,nTimeIntervalsSpectral) :: Aad, Bad, Cad, Dad, Ead, Fad, Gad
  real(kind=realType) ::eye(nw,nw)
  ! idxmgb - global block row index
  ! idxngb - global block column index

  integer(kind=intType) :: idxmgb, idxngb,ierr, sps, sps2,ilow,ihigh

  !Reference values of the dissipation coeff for the preconditioner
  real(kind=realType) :: vis2_ref, vis4_ref

  ! Set the grid level of the current MG cycle, the value of the
  ! discretization and the logical correctForK.
  level = 1_intType
  currentLevel = level
  time(1) = mpi_wtime()
  rkStage = 0
  currentLevel = groundLevel

  !     ******************************************************************
  !     *                                                                *
  !     * Compute the ADjoint matrix dR/dW using Tapenade's reverse mode *
  !     * of Automatic Differentiation.  NOTE: This is the reason I have *
  !     * been writing the word "ADjoint" with A and D capitalized. A    *
  !     * simple play with letter so that:                               *
  !     *                                                                *
  !     * ADjoint = Automatically Differentiated adjoint                 *
  !     *                                                                *
  !     ******************************************************************
  !
  ! Send some feedback to screen.

  if (myid == 0) then
     print * ,"Assembling NK KSP PC matrix..."
  end if
 
  !store the current values of vis2,vis4 and reset vis2 for preconditioner
  !method based on (Hicken and Zingg,2008) AIAA journal,vol46,no.11
  vis2_ref = vis2
  vis4_ref = vis4
  lumpedDiss=.True.
  
  call MatZeroEntries(dRdwPre,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  domainLoopAD: do nn=1,nDom

     ! Loop over the number of time instances for this block.
     spectralLoop: do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,level,sps)
        ! Loop over location of output (R) cell of residual
        do kCell = 2, kl
           do jCell = 2, jl
              do iCell = 2, il
                 ! Copy the state w to the wAdj array in the stencil

                 call copyNKPCStencil(iCell, jCell, kCell, nn, level, sps, wAdj, &
                      siAdj, sjAdj, skAdj, sAdj, sfaceIAdj, sfaceJAdj, sfaceKAdj, rotRateAdj,&
                      voladj)
                 Aad(:,:,:)  = zero
                 Bad(:,:,:)  = zero
                 Cad(:,:,:)  = zero
                 Dad(:,:,:)  = zero
                 Ead(:,:,:)  = zero
                 Fad(:,:,:)  = zero
                 Gad(:,:,:)  = zero

                 mLoop: do m = 1, nw      ! Loop over output cell residuals (R)
                    ! Initialize the seed for the reverse mode
                    dwAdjb(:,:) = 0.
                    dwAdjb(m,sps) = 1.
                    dwAdj(:,:)  = 0.
                    wAdjb(:,:,:,:,:)  = 0.  !dR(m)/dw

                    ! Call the reverse mode of residual computation.
                    !
                    !                          dR(iCell,jCell,kCell,l)
                    ! wAdjb(ii,jj,kk,n) = --------------------------------
                    !                     dW(iCell+ii,jCell+jj,kCell+kk,n)

                    ! Call reverse mode of residual computation

                    call COMPUTERNKPC_B(wadj, wadjb, dwadj, dwadjb, siadj, sjadj, &
                         &  skadj, sadj, voladj, sfaceiadj, sfacejadj, sfacekadj, rotrateadj, &
                         &  icell, jcell, kcell, nn, level, sps)


                !     call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
!                          &  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
!                          &  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
!                          &  icell, jcell, kcell, nn, level, sps, correctfork, secondhalo, prefadj&
!                          &  , rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj&
!                          &  , rotrateadjb, rotcenteradj, rotcenteradjb, pointrefadj, pointrefadjb&
!                          &  , rotpointadj, rotpointadjb, murefadj, timerefadj, pinfcorradj, &
!                          &  liftindex)


                    ! Store the block Jacobians (by rows).

                    Aad(m,:,:)  = wAdjB( 0, 0, 0,:,:)
                    Bad(m,:,:)  = wAdjB(-1, 0, 0,:,:)
                    Cad(m,:,:)  = wAdjB( 1, 0, 0,:,:)
                    Dad(m,:,:)  = wAdjB( 0,-1, 0,:,:)
                    Ead(m,:,:)  = wAdjB( 0, 1, 0,:,:)
                    Fad(m,:,:)  = wAdjB( 0, 0,-1,:,:)
                    Gad(m,:,:)  = wAdjB( 0, 0, 1,:,:)
                 enddo mLoop
                 ! Global matrix block row mgb function of node indices.

                 idxmgb = globalCell(iCell,jCell,kCell)

                 ! >>> center block A < W(i,j,k)
                 do sps2 = 1,nTimeIntervalsSpectral
                    idxngb = flowDoms(nn,level,sps2)%globalCell(iCell,jCell,kCell)
                    call MatSetValuesBlocked(dRdWPre, 1, idxmgb, 1, idxngb,Aad(:,:,sps2), &
                         ADD_VALUES,ierr)
                    call EChk(ierr,__FILE__,__LINE__)
                 enddo

                 ! >>> west block B < W(i-1,j,k)
                 if( (iCell-1) >= 0 ) then
                    idxngb = globalCell(iCell-1,jCell,kCell)
                    if (idxngb >=0 .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdWPre, 1, idxmgb, 1, idxngb, Bad(:,:,sps), &
                            ADD_VALUES,ierr)
                       call EChk(ierr,__FILE__,__LINE__)
                    endif
                 endif

                 ! >>> east block C < W(i+1,j,k)
                 if( (iCell+1) <= ib ) then
                    idxngb = globalCell(iCell+1,jCell,kCell)
                    if (idxngb<nCellsGlobal*nTimeIntervalsSpectral .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdWPre, 1, idxmgb, 1, idxngb, Cad(:,:,sps), &
                            ADD_VALUES,ierr)
                       call EChk(ierr,__FILE__,__LINE__)
                    endif
                 end if

                 ! >>> south block D < W(i,j-1,k)
                 if( (jCell-1) >= 0 ) then
                    idxngb = globalCell(iCell,jCell-1,kCell)
                    if (idxngb>=0 .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdWPre, 1, idxmgb, 1, idxngb, Dad(:,:,sps), &
                            ADD_VALUES,ierr)
                       call EChk(ierr,__FILE__,__LINE__)
                    endif
                 endif

                 ! >>> north block E < W(i,j+1,k)
                 if( (jCell+1) <= jb ) then
                    idxngb = globalCell(iCell,jCell+1,kCell)
                    if (idxngb<nCellsGlobal*nTimeIntervalsSpectral .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdWPre, 1, idxmgb, 1, idxngb, Ead(:,:,sps), &
                            ADD_VALUES,ierr)
                       call EChk(ierr,__FILE__,__LINE__)
                    endif
                 end if

                 ! >>> back block F < W(i,j,k-1)
                 if( (kCell-1) >= 0 ) then
                    idxngb = globalCell(iCell,jCell,kCell-1)
                    if (idxngb>=0 .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdWPre, 1, idxmgb, 1, idxngb, Fad(:,:,sps), &
                            ADD_VALUES,ierr)
                       call EChk(ierr,__FILE__,__LINE__)
                    endif
                 endif

                 ! >>> front block G < W(i,j,k+1)
                 if( (kCell+1) <= kb ) then
                    idxngb = globalCell(iCell,jCell,kCell+1)
                    if (idxngb<nCellsGlobal*nTimeIntervalsSpectral .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdWPre, 1, idxmgb, 1, idxngb, Gad(:,:,sps), &
                            ADD_VALUES,ierr)
                       call EChk(ierr,__FILE__,__LINE__)
                    endif
                 end if
              enddo
           enddo
        enddo
     enddo spectralLoop
  enddo domainLoopad

  !Return dissipation Parameters to normal
  vis2 = vis2_ref
  vis4 = vis4_ref

  call MatAssemblyBegin(dRdWPre,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call MatAssemblyEnd  (dRdWPre,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#ifdef USE_PETSC_3
  call MatSetOption(dRdWPre,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#else
  call MatSetOption(dRdWPre,MAT_NO_NEW_NONZERO_LOCATIONS,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif

  time(2) = mpi_wtime()
  call mpi_reduce(time(2)-time(1),setupTime,1,sumb_real,mpi_max,0,&
       SUmb_comm_world, ierr)

  if (myid == 0) then
     print *,'Done PC Assembly'
     print *,'Time:',setupTime
  end if
end subroutine setupNK_KSP_PC

subroutine getCurrentResidual(rhoRes,totalRRes)
  use communication
  use blockPointers
  use flowVarRefState
  use inputTimeSpectral
  use iteration
  use inputIteration
  implicit none
  ! Compute the current resdiual of w
  real(kind=realType), intent(out) :: rhoRes,totalRRes
  real(kind=realType) :: ovv,r_sum,rho_sum
  integer(kind=intType) :: sps,nn,i,j,k,l,ierr
  currentLevel = 1
  groundLevel = 1
  rkStage = 0

  call timestep(.false.)
  call initres(1_intType, nwf)
  call residual 
 
  r_sum = 0.0
  rho_sum = 0.0
  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointers(nn,1_intType,sps)
        ! Copy off dw/vol to rVec
        do k=2,kl
           do j=2,jl
              do i=2,il
                 ovv = 1/vol(i,j,k)
                 do l=1,nw
                    r_sum = r_sum + (dw(i,j,k,l)*ovv)**2
                 end do
                 rho_sum = rho_sum + (dw(i,j,k,irho)*ovv)**2
              end do
           end do
        end do
     end do
  end do

  call mpi_allreduce(r_sum,totalRRes,1,sumb_real,mpi_sum,&
       SUmb_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call mpi_allreduce(rho_sum,rhoRes,1,sumb_real,mpi_sum,&
       SUmb_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ! curRes now has the inverse-volume weighted sum squared of the
  ! residuals, finally take the squareRoot

 totalRRes = sqrt(totalRRes)
 rhoRes = sqrt(rhoRes/nCellGlobal(currentLevel))

end subroutine getCurrentResidual

subroutine getFreeStreamResidual(rhoRes,totalRRes)
  use communication
  use ADjointVars
  use precision
  use inputTimeSpectral
  use flowVarRefState
  use inputIteration
  use blockpointers
  implicit none

  real(kind=realType), intent(out) :: rhoRes,totalRRes
  real(kind=realType),dimension(:), allocatable :: wtemp
  integer(kind=intType) :: nDimW,ierr,tempStartLevel,counter
  integer(kind=intType) :: nn,sps,i,j,k,l
  ! Get the residual cooresponding to the free-stream on the fine grid-level

  ! We need to copy the current wvector temporirly since it may be a
  ! restart and actually useful
  nDimW = nw * nCellsLocal * nTimeIntervalsSpectral
  allocate(wtemp(nDimW))

  ! Copy w to wTemp
  counter = 0
  spectralLoop: do sps=1,nTimeIntervalsSpectral
     domains: do nn=1,nDom
        call setPointers(nn,1,sps)
        do l=1,nw
           do k=2,kl
              do j=2,jl
                 do i=2,il
                    counter = counter + 1
                    wtemp(counter) = w(i,j,k,l)
                 enddo
              enddo
           enddo
        enddo
     end do domains
  end do spectralLoop
  
  tempStartLevel = mgStartLevel
  mgStartLevel = 1

  call setUniformFlow
  call getCurrentResidual(rhoRes,totalRRes)
  counter = 0
  spectralLoop2: do sps=1,nTimeIntervalsSpectral
     domains2: do nn=1,nDom
        call setPointers(nn,1,sps)
        do l=1,nw
           do k=2,kl
              do j=2,jl
                 do i=2,il
                    counter = counter + 1
                    w(i,j,k,l) = wtemp(counter) 
                 enddo
              enddo
           enddo
        enddo
     end do domains2
  end do spectralLoop2

  mgStartLevel = tempStartLevel
  deallocate(wtemp)
end subroutine getFreeStreamResidual

subroutine snes_monitor(snes,its,norm,ctx,ierr)
  use communication
  use precision 
  use iteration
  use inputIteration
  use NKsolverVars, only: ksp_rtol,ksp_atol,ksp_div_tol,ksp_max_it,&
       snes_atol,itertot0,jacobian_lag
  implicit none
#define PETSC_AVOID_MPIF_H
#if PETSC_VERSION_MINOR>=1
#include "include/finclude/petsc.h"
#else
#include "include/finclude/petscall.h"
#endif

  SNES snes
  KSP  ksp
  PetscInt its, ierr
  PetscReal norm
  PetscFortranAddr ctx(*) ! This is probably going to be empty

  integer(kind=intType) :: ksp_its,temp
  real(kind=realType) :: CFLNew,rhoRes,totalRRes
  ! We want to get the number of iterations of the last KSP and
  ! increment iterTot by this. This will give the user an indication
  ! of how long each KSP is taking to solver

  call SNESGetKSP(snes,ksp,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call KSPGetTolerances(ksp,ksp_rtol,ksp_atol,ksp_div_tol,ksp_max_it,ierr)
  ksp_atol = snes_atol
  ksp_max_it = 100
  call KSPSetTolerances(ksp,ksp_rtol,ksp_atol,ksp_div_tol,ksp_max_it,ierr)

  if (its == 1) then
     ! Reset the value of the jacobianLag to what we actually want. It
     ! had been set to -1 or -2 depending on if we wanted to recompute
     ! the preconditioner on the first entry or not. 
     call SNESSetLagJacobian(snes, jacobian_lag, ierr); call EChk(ierr,__FILE__,__LINE__)
  end if

  if (its > 0 .or. iterTot0 == 0) then
     call SNESGetLinearSolveIterations(snes,ksp_its,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ksp_its = max(ksp_its,1)
     iterTot = iterTot0 + ksp_its

     call convergenceInfo

     ierr = 0
  end if
end subroutine snes_monitor

