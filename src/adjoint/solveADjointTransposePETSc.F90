!
!      File:          solveAdjoint.F90                                
!      Author:        Gaetan K.W. Kenway                              
!      Starting date: 03-12-2014                                      
!      Last modified: 03-12-2014                                      
!
subroutine solveAdjoint(RHS, psi, checkSolution, nState)
  !
  !      Solve the linear discrete ADjoint system of equations          
  !          [dR/dW]T . psi = {RHS}                                     
  !      using preconditioned GMRES provided by PETSc. The values in psi
  !      are significant as they are used as the inital guess.          
  !
  use constants
  use ADjointPETSc, only : dRdwT, psi_like1, psi_like2, adjointKSP, &
       adjResInit, adjResStart, adjResFinal
 
  use killsignals, only : adjointFailed
  use inputADjoint, only : adjAbsTol, adjDivTol, adjMaxIter, adjRelTol, &
       adjRelTolRel, printTiming
  use adjointVars, only: derivVarsAllocated
  use communication, only : myid, sumb_comm_world
  use blockPointers, only : nDom
  use inputTimeSpectral, only : nTimeIntervalsSpectral
  use utils, only : EChk
  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

  ! Input Parameters
  real(kind=realType), dimension(nState) :: RHS, psi
  integer(kind=intType) :: nState
  logical :: checkSolution 
  !
  !     Local variables.
  real(kind=realType)   :: norm
  real(kind=realType), dimension(2) :: time
  real(kind=realType)               :: timeAdjLocal, timeAdj
  real(kind=realType) :: l2abs, l2rel
  integer(kind=intType) :: ierr, nn, sps
  integer(kind=intType) :: adjConvIts
  KSPConvergedReason adjointConvergedReason
  Vec adjointRes, RHSVec

  ! Send some feedback to screen.

  if(myid ==0 .and. printTiming)  &
       write(*,10) "Solving ADjoint Transpose with PETSc..."

  call cpu_time(time(1))

  ! Make sure the derivative memory is allocated and zeroed. 
  if (.not. derivVarsAllocated) then 
     call allocDerivativeValues(1_intType)
  end if

  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call zeroADSeeds(nn, 1_intType, sps)
     end do
  end do

  ! Dump psi into psi_like1 and RHS into psi_like2
  call VecPlaceArray(psi_like1, psi, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecPlaceArray(psi_like2, RHS, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDuplicate(psi_like1, adjointRes, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  if (checkSolution) then 
     call VecDuplicate(psi_like1, RHSVec, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call vecCopy(psi_like2, RHSVec, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end if

  ! Get the RHS norm....this is the 'init' norm:
  call VecNorm(psi_like2, NORM_2, adjResInit, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Get Current Residual -- we always solve for the delta
  call MatMult(dRdWT, psi_like1, adjointRes, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! AdjointRes = AdjointRes - adjointRHS
  call VecAXPY(adjointRes, -one, psi_like2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Norm of adjoint Residual
  call VecNorm(adjointRes, NORM_2, adjResStart,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! The way we use tolerances are as follows: The residual must
  ! statify:
  ! res < adjRelTol * adjResInit OR 
  ! res < adjRelTolRel * adjResStart OR
  ! res < adjAbsTol

  ! L2Abs is used to stipulate an exit criteria for adjreltolrel
  L2abs = adjResStart * adjreltolrel

  ! If L2Abs is less that what we actually want as the absolute
  ! tolerance, clip it
  if (L2Abs < adjAbsTol) then
     L2abs = adjabstol
  end if

  ! L2Rel is a little tricky since if the start residual is *larger*
  ! than the init residual, it won't converge enough. While this seems
  ! strange this is *always* the case for restarted RANS-based
  ! adjoints.
  L2Rel = (adjReltol * adjResInit) / adjResStart

  ! We need to clip L2Rel such that it can never be greater than one. 
  L2Rel = min(L2Rel, 0.9)

  ! Set the tolerances
  call KSPSetTolerances(adjointKSP, L2Rel, L2Abs, adjDivTol, &
       adjMaxIter, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Solve the update (psi_like2)
  call KSPSolve(adjointKSP, adjointRes, psi_like2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now compute the update to psi_like1 (psi)
  call VecAXPY(psi_like1, -one, psi_like2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (checkSolution) then 

     ! Get new time and compute the elapsed time.
     call cpu_time(time(2))
     timeAdjLocal = time(2)-time(1)

     ! Determine the maximum time using MPI reduce
     ! with operation mpi_max.

     call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
          mpi_max, 0, SUMB_COMM_WORLD, ierr)
   
     call MatMult(dRdWT, psi_like1, adjointRes, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecAXPY(adjointRes, -one, RHSVec, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecNorm(adjointRes, NORM_2, norm,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     adjResFinal = norm

     call KSPGetIterationNumber(adjointKSP,adjConvIts,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     ! Use the root processor to display the output summary, such as
     ! the norm of error and the number of iterations
     
     if( myid ==0 .and. printTiming) then
        write(*,20) "Solving ADjoint Transpose with PETSc time (s) =", timeAdj
        write(*,30) "Norm of error =",norm,"Iterations =",adjConvIts
        write(*,*) "------------------------------------------------"
        if( adjConvIts.lt.0 ) then
           write(*,40) "PETSc solver diverged after", -adjConvIts, &
                "iterations..."
        else
           write(*,40) "PETSc solver converged after", adjConvIts, &
                "iterations."
        endif
          write(*,*) "------------------------------------------------"
       endif

       call VecDestroy(RHSVec, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end if
    
    ! Destroy the temporary vector and reset the arrays
    call VecDestroy(adjointRes, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecResetArray(psi_like1, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecResetArray(psi_like2, ierr)
    call EChk(ierr,__FILE__,__LINE__)
    
    ! Get the petsc converged reason and set the fail flag
    call KSPGetConvergedReason(adjointKSP, adjointConvergedReason,ierr)
    call EChk(ierr,__FILE__,__LINE__)
    
    if (adjointConvergedReason ==  KSP_CONVERGED_RTOL .or. &
         adjointConvergedReason ==  KSP_CONVERGED_ATOL .or. &
         adjointConvergedReason ==  KSP_CONVERGED_HAPPY_BREAKDOWN) then
       adjointFailed = .False.
    else
       adjointFailed = .True.
    end if
    
  ! Output formats.

10 format(a)
20 format(a,1x,f8.2)
30 format(1x,a,1x,e10.4,4x,a,1x,i4)
40 format(1x,a,1x,i5,1x,a)

end subroutine solveAdjoint

