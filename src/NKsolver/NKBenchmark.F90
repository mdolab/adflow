! This file contains routine to test the parallel scaling of various
! parts of the Newton-Krylov algorithim implementation. It is strictly
! concerned ONLY with parallel communication efficiency, not
! algorithmic efficiency. As such, no "problems" are actually
! solved. Only the "core" of the solver is run multiple times to produce
! timing results. 

! The main component of the NK algorithim is the "Function
! Evaluation". This is essentially a single fine grid residual
! evaluation along with a double-halo exchange of data. Additionally,
! we also test the application of the global preconditioner as that
! involves communication if the Additive Schwartz Preconditioner. 

! The subroutine takes various options to test different components:

! EP, logical: Run just the residual in an embarrassing parallel
!      fashion. All other logicals are meaningless
! res, logical: Run the residual evaluation routine
! applyPC, logical: Run the applyPC routine of the preconditioner. 
! niter: Number of iterations to run


subroutine NKBenchmark(NKRes,niter)
#ifndef USE_NO_PETSC
  use communication
  use flowVarRefState
  use inputtimespectral
  use blockPointers
  use nksolvervars, only :wVec, rVec, ctx, drdw, deltaw
  use inputPhysics
  use iteration
  use inputiteration
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Input
  logical, intent(in) :: NKRes
  integer(kind=intType), intent(in) :: niter

  ! Working
  integer(kind=intType) :: iter, ierr
  real(kind=realType) :: timeA,timeB,timeX,timeY
  real(kind=realType) :: resTime, totalTime
  real(kind=realType),dimension(:), allocatable :: all_times
  ! Setup NK solver in case it isn't already done:
  call setupNKsolver
  
  ! Set the current w in to the wvec
  call setwVec(wVec)
  
  ! Reset all times

  resTime = zero
 
  if (.not. NKres) then
     allocate(cycling(nMGSteps), stat=ierr)
     call setCycleStrategy
  else
     call setupNKSolver()
     call FormJacobian()
     call MatMFFDSetBase(dRdW,wVec,PETSC_NULL_OBJECT,ierr)
  end if
  timeX = mpi_wtime()
  do iter=1,niter
   
     if (NKres) then

        ! Run the 'core' computeResidualNK
        timeA = mpi_wtime()

        call formFunction_mf(ctx,wVec,rVec,ierr)
        !call PCApply(global_pc,rVec,deltaW,ierr)
        
        timeB = mpi_wtime()

        resTime = resTime + timeB - timeA

     else ! Use the MG Res

        ! Run the core of the MG solver...executeMGcycle
        timeA = mpi_wtime()
        call executeMGCycle()
        timeB = mpi_wtime()
        
        resTime = resTime + timeB - timeA

     end if

  end do
  timeY = mpi_wtime()
  allocate(all_times(nproc))

  ! Collect all the Restimes
  call MPI_Gather (resTime,1,sumb_real,all_times,1,sumb_real,0,&
       SUMB_COMM_WORLD,ierr) 
  if (myid == 0) then
     !print *,'Min Res Time:',minval(all_times)/niter,minloc(all_times)-1
     !print *,'Max Res Time:',maxval(all_times)/niter,maxloc(all_times)-1
     !print *,'Average Res Time:',sum(all_times)/nproc/niter
  end if
  deallocate(all_times)

  ! Determine total Time:
  call MPI_reduce(timeY-timeX,totalTime,1,sumb_real,MPI_MAX,0,&
       SUMB_COMM_WORLD,ierr)

  if (myid == 0) then
     print *,'Total Time:',totalTime/niter
   print *, ' '
  end if

  if (.not. NKres) then
     deallocate(cycling)
  end if
#endif
end subroutine NKBenchmark



