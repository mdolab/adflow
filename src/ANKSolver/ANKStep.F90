subroutine ANKStep(firstCall)

  use constants
  use flowVarRefState
  use NKSolverVars, only : totalR0
  use ANKSolverVars, only: ANK_jacobianLag,  wVec, rVec, deltaW, &
       ANK_KSP, ANK_rTol, ANK_subSpace, ANK_divTol, &
       ANK_Iter, ANK_CFL, NORM_2
  use inputIteration
  use inputPhysics
  use monitor
  use iteration
  implicit none

  ! Input Variables
  logical, intent(in) :: firstCall

  ! Working Variables
  integer(kind=intType) :: ierr, maxIt, kspIterations
  real(kind=realType) :: norm, atol

  if (firstCall) then 
     call setupANKSolver()

     ! Copy the sumb 'w' into the petsc wVec
     call setwVecANK(wVec)
  
     ! Evaluate the residual before we start and put the residual in
     ! 'g', which is what would be the case after a linesearch.
     call computeResidualANK()
     call setRVecANK(rVec)
  else
     ANK_iter = ANK_iter + 1
  end if

  ! Determine if if we need to form the Preconditioner
  if (mod(ANK_iter, ANK_jacobianLag) == 0) then

     ! Compute the norm of rVec to update the CFL
     call VecNorm(rVec, NORM_2, norm, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ANK_CFL = min(CFL * totalR0 / norm, 100000.0)

     call FormJacobianANK()
  end if

  ! Set all tolerances for linear solve:
  atol = totalR0*L2Conv
  call KSPSetTolerances(ANK_KSP, real(ANK_rtol), &
       real(atol), real(ANK_divTol), ANK_subSpace, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ! Actually do the Linear Krylov Solve
  call KSPSolve(ANK_KSP, rVec, deltaW, ierr)

  ! DON'T just check the error. We want to catch error code 72
  ! which is a floating point error. This is ok, we just reset and
  ! keep going
  if (ierr == 72) then
     ! The convergence check will get the nan
  else
     call EChk(ierr, __FILE__, __LINE__)
  end if

  ! No line search...just take the new solution
  call VecAXPY(wVec, -1.0_realType, deltaW, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Set the updated state variables
  call setWANK(wVec)

  ! Compute new function for next iteration.
  call computeResidualANK()
  call setRVecANK(rVec)

  if (equations==RANSEquations) then 
     call computeUtau
     call turbSolveSegregated
  end if

  ! Update the approximate iteration counter. The +1 is for the
  ! residual evaluation. 

  call KSPGetIterationNumber(ANK_KSP, kspIterations, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  approxTotalIts = approxTotalIts + 1 + kspIterations

end subroutine ANKStep
