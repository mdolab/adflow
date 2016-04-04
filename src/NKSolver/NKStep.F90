subroutine NKStep(firstCall)
#ifndef USE_NO_PETSC

  use communication
  use constants
  use inputTimeSpectral
  use flowVarRefState
  use NKSolverVars, only: dRdw, dRdwPre, NK_jacobianLag, totalR0, wVec, rVec, &
       deltaW, NK_KSP, NK_subspace, NK_divTol, NK_LS, NK_useEw, NK_iter, &
       nolinesearch, cubiclinesearch, nonmonotonelinesearch, &
       work, g, NK_rtolInit, NK_CFL, NK_CFL0, oldNorm, rtolLast
  use InputIO 
  use inputIteration
  use inputPhysics
  use monitor
  use killSignals
  use iteration
  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

  ! Input Variables
  logical, intent(in) :: firstCall

  ! Working Variables
  integer(kind=intType) :: iter, ierr, kspIterations
  integer(kind=intType) :: maxNonLinearIts, nfevals, maxIt
  real(kind=realType) :: norm, rtol, atol
  real(kind=realType) :: fnorm, ynorm, gnorm
  logical :: flag

  if (firstCall) then 
     call setupNKSolver()

     ! Copy the sumb 'w' into the petsc wVec
     call setwVec(wVec)
  
     ! Evaluate the residual before we start and put the residual in
     ! 'g', which is what would be the case after a linesearch.
     call computeResidualNK()
     call setRVec(rVec)
  else
     NK_iter = NK_iter + 1
  end if
 
  ! Compute the norm of rVec for use in EW Criteria
  call VecNorm(rVec, NORM_2, norm, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Determine if if we need to form the Preconditioner
  if (mod(NK_iter, NK_jacobianLag) == 0) then
     NK_CFL = NK_CFL0 * totalR0 / norm !* nk_switch_tol
     call FormJacobianNK()
  else
     call MatAssemblyBegin(dRdw, MAT_FINAL_ASSEMBLY, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     call MatAssemblyEnd(dRdw, MAT_FINAL_ASSEMBLY, ierr)
     call EChk(ierr, __FILE__, __LINE__)
  end if
 
  ! set the BaseVector of the matrix-free matrix:
  call MatMFFDSetBase(dRdW, wVec, PETSC_NULL_OBJECT, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  if (NK_iter == 1 .or. .not. NK_useEW) then 
     rtol = NK_rtolInit
  else
    call getEWTol(norm, oldNorm, rtolLast, rtol)
 end if

  ! Save the old rtol and norm for the next iteration
  oldNorm = norm
  rtolLast = rtol

  ! Set all tolerances for linear solve:
  atol = totalR0*L2Conv
  maxIt = NK_subspace

  call KSPSetTolerances(NK_KSP, real(rtol), &
       real(atol), real(NK_divTol), maxIt, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ! Actually do the Linear Krylov Solve
  call KSPSolve(NK_KSP, rVec, deltaW, ierr)

  ! DON'T just check the error. We want to catch error code 72
  ! which is a floating point error. This is ok, we just reset and
  ! keep going
  if (ierr == 72) then
     ! The convergence check will get the nan
  else
     call EChk(ierr, __FILE__, __LINE__)
  end if

  nfevals = 0
  if (NK_LS == noLineSearch) then
     call LSNone(wVec, rVec, g, deltaW, work, nfevals, flag)
  else if(NK_LS == cubicLineSearch) then
     call LSCubic(wVec, rVec, g, deltaW, work, fnorm, ynorm, gnorm, &
          nfevals, flag)
  else if (NK_LS == nonMonotoneLineSearch) then

     print *, 'Not implemented yet'
     stop
     call LSNM(wVec, rVec, g, deltaW, work, fnorm, ynorm, gnorm, &
          nfevals, flag)
  end if
  
  if (.not. flag) then
     routineFailed = .True.
  end if

  ! Copy the work vector to wVec. This is our new state vector
  call VecCopy(work, wVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Use the result from the line sesarch for the residual
  call vecCopy(g, rVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)


  ! Update the approximate iteration counter. The +nFevals is for the
  ! iterations taken during the linesearch

  call KSPGetIterationNumber(NK_KSP, kspIterations, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  approxTotalIts = approxTotalIts + nfEvals + kspIterations

#endif
end subroutine NKStep
