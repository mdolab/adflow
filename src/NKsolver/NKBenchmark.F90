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


subroutine NKBenchmark(EP,res,applyPC,niter)

  use communication
  use precision
  use flowVarRefState
  use inputtimespectral
  use blockPointers
  use nksolvervars, only : petscComm,wVec
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Input
  logical, intent(in) :: EP,res,applyPC
  integer(kind=intType), intent(in) :: niter

  ! Working
  integer(kind=intType) :: iter,i,ierr
  real(kind=realType) :: timeA,timeB,timeX,timeY
  real(kind=realType) :: haloTime,resTime,totalTime
  real(kind=realType),dimension(:), allocatable :: all_times
  ! Setup NK solver in case it isn't already done:
  call setupNKsolver
  
  ! Set the current w in to the wvec
  call setwVec(wVec)
  
  ! Reset all times
  haloTime = 0.0
  resTime = 0.0
  timeX = mpi_wtime()
  do iter=1,niter
   
     if (res) then

        timeA = mpi_wtime()
        if (.not. EP) then
           call MPI_barrier(sumb_comm_world,ierr)

           if (petscComm) then
              call setW_ghost(wVec)
           else
              call setW(wVec)
           end if

        end if
        timeB = mpi_wtime()
        haloTime = haloTime + timeB - timeA

        ! Run the 'core' computeResidualNK
     
        call MPI_barrier(sumb_comm_world,ierr)

        timeA = mpi_wtime()
        call computeResidualNK()
        timeB = mpi_wtime()

        resTime = resTime + timeB - timeA

     end if

     if (applyPC .and. .not. EP) then
        print *,'Applying PC...'
     end if
  end do
  timeY = mpi_wtime()
  allocate(all_times(nproc))

  ! Collect all the Halotimes
  call MPI_Gather (haloTime,1,sumb_real,all_times,1,sumb_real,0,&
       SUMB_COMM_WORLD,ierr) 
  if (myid == 0) then
     do i=1,nproc
        print *,i-1,'Halo Time: ',all_times(i)/niter
     end do
     print *,' '
     print *,'Average Halo Time:',sum(all_times)/nproc/niter
     print *,' '

  end if
  ! Collect all the Restimes
  call MPI_Gather (resTime,1,sumb_real,all_times,1,sumb_real,0,&
       SUMB_COMM_WORLD,ierr) 
  if (myid == 0) then
     do i=1,nproc
        print *,i-1,'Res Time: ',all_times(i)/niter
     end do
     print *,' '
     print *,'Average Res Time:',sum(all_times)/nproc/niter
     print *,' '
  end if
  deallocate(all_times)

  ! Determine total Time:
  call MPI_reduce(timeY-timeX,totalTime,1,sumb_real,MPI_MAX,0,&
       SUMB_COMM_WORLD,ierr)

  if (myid == 0) then
     print *, ' '
     print *,'Total Time:',totalTime/niter
     print *, ' '
  end if

end subroutine NKBenchmark



