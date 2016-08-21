subroutine ANKStep(firstCall)

  use constants
  use flowVarRefState, only : nw
  use NKSolverVars, only : totalR0
  use ANKSolverVars, only: ANK_jacobianLag,  wVec, rVec, deltaW, &
       ANK_KSP, ANK_rTol, ANK_subSpace, ANK_divTol, &
       ANK_Iter, ANK_CFL, NORM_2, wVecTurb, rVecTurb, deltaWTurb, ANK_KSPTurb, &
       ANK_useTurbDADI, ANK_CFL0
  use inputPhysics, only : equations
  use flowVarRefState, only :  nw, nwf
  use inputIteration, only : L2conv
  use iteration, only : approxTotalIts
  implicit none

  ! Input Variables
  logical, intent(in) :: firstCall

  ! Working Variables
  integer(kind=intType) :: ierr, maxIt, kspIterations, j
  real(kind=realType) :: norm, atol, val

  if (firstCall) then 
     call setupANKSolver()
     call destroyANKSolver()
     call setupANKSolver()

     ! Copy the sumb 'w' into the petsc wVec
     call setwVecANK(wVec)

     ! Evaluate the residual before we start and put the residual in
     ! 'g', which is what would be the case after a linesearch.
     if (ANK_useTurbDADI) then 
        call computeResidualANK() ! Only flow residual
     else
        ! Compute rull residual
        call computeResidualNK()
        if (nw > nwf) then 
           call setwVecANKTurb(wVecTurb)
        end if
     end if

     call setRVecANK(rVec)
  else
     ANK_iter = ANK_iter + 1
  end if

  ! Determine if if we need to form the Preconditioner
  if (mod(ANK_iter, ANK_jacobianLag) == 0) then

     ! Compute the norm of rVec to update the CFL
     call VecNorm(rVec, NORM_2, norm, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ANK_CFL = min(ANK_CFL0 * (totalR0 / norm)**1.5, 100000.0)

     call FormJacobianANK()
  end if

  ! ============== Flow Update =============
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

  ! ============== Turb Update =============
  if (equations==RANSEquations) then 
     if (ANK_useTurbDADI) then 
        call computeUtau
        call turbSolveSegregated
     else
        atol = 1e-30
        call KSPSetTolerances(ANK_KSPTurb, real(ANK_rtol), &
             real(atol), real(ANK_divTol), ANK_subSpace, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        
        ! Actually do the Linear Krylov Solve
        call KSPSolve(ANK_KSPTurb, rVecTurb, deltaWTurb, ierr)
     
        ! DON'T just check the error. We want to catch error code 72
        ! which is a floating point error. This is ok, we just reset and
        ! keep going
        if (ierr == 72) then
           ! The convergence check will get the nan
        else
           call EChk(ierr, __FILE__, __LINE__)
        end if
     
        ! No line search...just take the new solution
        call VecAXPY(wVecTurb, -.250_realType, deltaWTurb, ierr)
        call EChk(ierr, __FILE__, __LINE__)
     
        ! Set the updated turbulent state variables
        call setWANKTurb(wVecTurb)
     end if
  end if

  ! ==============================================

  ! Compute new function for next iteration.
  if (ANK_useTurbDADI) then 
     ! Only Flow Residual
     call computeResidualANK()
  else
     ! Ful Residual
     call computeResidualNK()
     call setRVecANKTurb(rVecTurb)
  end if
  call setRVecANK(rVec)

  ! Update the approximate iteration counter. The +1 is for the
  ! residual evaluation. 

  call KSPGetIterationNumber(ANK_KSP, kspIterations, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  approxTotalIts = approxTotalIts + 1 + kspIterations

end subroutine ANKStep

