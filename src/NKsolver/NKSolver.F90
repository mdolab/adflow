!
!      ******************************************************************
!      *                                                                *
!      * File:          NKsolver.F90                                    *
!      * Author:        Gaetan Kenway                                   *
!      * Starting date: 11-27-2010                                      *
!      * Last modified: 11-27-2010                                      *
!      *                                                                *
!      ******************************************************************

subroutine NKsolver
#ifndef USE_NO_PETSC
  use communication
  use constants
  use inputTimeSpectral
  use flowVarRefState
  use ADjointVars , only: nCellsLocal
  use NKSolverVars, only: dRdw, dRdwPre, jacobian_lag, &
       totalR0, totalRStart, wVec, rVec, deltaW, global_ksp, &
       ksp_rtol, ksp_atol, func_evals, ksp_max_it, ksp_subspace, ksp_div_tol, &
       nksolvecount, Mmax, iter_k, iter_m, NKLS, &
       nolinesearch, cubiclinesearch, nonmonotonelinesearch, rhores0, &
       NK_switch_tol, rhoresstart, work, g

  use InputIO ! L2conv,l2convrel
  use inputIteration
  use monitor
  use killSignals
  use iteration
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Working Variables
  integer(kind=intType) :: iter,ierr,ksp_iterations
  integer(kind=intType) :: maxNonLinearIts,nfevals
  real(kind=realType) :: norm,old_norm,rtol_last
  real(kind=realType) :: fnorm,ynorm,gnorm
  logical :: flag
  ! maxNonLinearIts is (far) larger that necessary. The "iteration"
  ! limit is really set from the maxmimum number of funcEvals
  maxNonLinearIts = ncycles-iterTot
  ksp_iterations = 0
  norm = 0.0
  old_norm=0.0
  rtol_last =0.0
  nfevals = 0

  Mmax = 10
  iter_k = 1
  iter_m = 0

  allocate(func_evals(1000))
  func_evals = 0.0
  ! Set the inital wVec
  call setwVec(wVec)

   ! Evaluate the residual before we start and copy the value into g
  call setW(wVec)
  call computeResidualNK()
  call setRVec(rVec)

  call vecCopy(rVec,g,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  iterTot = iterTot + 1 ! Add this function evaluation

  ! Master Non-Linear Loop:
  NonLinearLoop: do iter= 1,maxNonLinearIts

     ! Increment the function evals from the Krylov Iterations and the
     ! line search iterations
     if (iter .ne. 1) then
        iterTot = iterTot + ksp_iterations + nfevals 
        
        ! If iterTot is anywhere near ncycles...we're cooked
        if (iterTot >= ncycles) then
           iterTot = ncycles
           ! We need to call convergence Info since this has the
           ! "approximate" convergence check
           call convergenceInfo
           exit NonLinearLoop
        else
           call convergenceInfo
        end if
     end if

     ! Use the result from the last line search
     call vecCopy(g,rVec,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Determine if if we need to form the Preconditioner: For a
     ! normal aerodynamic solution, form the preconditioner every
     ! 'jacobian_lag' iterations.  However, for an aerostructural
     ! analysis (with Gauss Sidel) the NK solver is called
     ! repedily. In this case we want the 'jacobian_lag' to refer to
     ! the number of time the solver was called. So on first
     ! iteration, we check the variable NKsolveCount. If THIS is a
     ! multiple of jacobian lag, we then reform the preconditioner. 


     if (mod(iter-1,jacobian_lag) == 0) then
        if (iter == 1) then ! Special case check:
           if (mod(NKsolveCount,jacobian_lag) == 0) then
              call FormJacobian()
           else
              call MatAssemblyBegin(dRdw,MAT_FINAL_ASSEMBLY,ierr)
              call EChk(ierr,__FILE__,__LINE__)
              call MatAssemblyEnd(dRdw,MAT_FINAL_ASSEMBLY,ierr)
              call EChk(ierr,__FILE__,__LINE__)
           end if
        else
           call FormJacobian()
        end if
     else
        ! Else just call assmebly begin/end on dRdW
        call MatAssemblyBegin(dRdw,MAT_FINAL_ASSEMBLY,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call MatAssemblyEnd(dRdw,MAT_FINAL_ASSEMBLY,ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if
     
     ! Set the BaseVector of the matrix-free matrix:
     call MatMFFDSetBase(dRdW,wVec,PETSC_NULL_OBJECT,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Compute the norm of rVec for use in EW Criteria
     old_norm = norm
     rtol_last = ksp_rtol
     call VecNorm(rVec,NORM_2,norm,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Check to see if we're converged: We need to check if we've meet
     ! L2Conv or L2ConvRel
     if (norm / totalR0 < L2Conv) then
        routineFailed = .False.
        exit NonLinearLoop
     end if

     if (norm / totalRStart < L2ConvRel) then
        routineFailed = .False.
        exit NonLinearLoop
     end if

     ! Check to see if we've done too many function Evals:
     if (iterTot >= ncycles) then
        iterTot = ncycles 
        exit NonLinearLoop
     end if

     ! Get the EW Forcing tolerance ksp_rtol
     call getEWTol(iter,norm,old_norm,rtol_last,ksp_rtol)
  
     ! Set all tolerances for linear solve:
     ksp_atol = totalR0*L2Conv
     ksp_max_it = min(ksp_subspace,ncycles-iterTot)
     ksp_max_it = max(ksp_max_it,1) ! At least one iteration!
    
     call KSPSetTolerances(global_ksp,ksp_rtol,ksp_atol,ksp_div_tol,&
          ksp_max_it,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Actually do the Linear Krylov Solve
     call KSPSolve(global_ksp,rVec,deltaW,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Linesearching:
     iter_k = iter
     iter_m = min(iter_m+1,Mmax)

     if (iter == 1) then
        call LSCubic(wVec,rVec,g,deltaW,work,fnorm,ynorm,gnorm,nfevals,flag)
     else
        if (NKLS == noLineSearch) then
           call LSNone(wVec,rVec,g,deltaW,work,nfevals,flag)
        else if(NKLS == cubicLineSearch) then
           call LSCubic(wVec,rVec,g,deltaW,work,fnorm,ynorm,gnorm,nfevals,flag)
        else if (NKLS == nonMonotoneLineSearch) then
           call LSNM(wVec,rVec,g,deltaW,work,fnorm,ynorm,gnorm,nfevals,flag)
        end if
     end if

     if (.not. flag) then
        routineFailed = .True.
        exit NonLinearLoop
     end if

     ! Copy the work vector to wVec
     call VecCopy(work,wVec,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Get the number of iterations to use with Convergence Info
     call KSPGetIterationNumber(global_ksp,ksp_iterations,ierr)
     call EChk(ierr,__FILE__,__LINE__)

  end do NonLinearLoop

  ! Not really anything else to do...

  NKSolveCount = NKSolveCount + 1

  deallocate(func_evals)

#endif
end subroutine NKsolver

subroutine LSCubic(x,f,g,y,w,fnorm,ynorm,gnorm,nfevals,flag)
#ifndef USE_NO_PETSC
  use constants
  use communication
  use NKSolverVars, only: dRdw
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Input/Output
  Vec x,f,g,y,w
  !x 	- current iterate
  !f 	- residual evaluated at x
  !y 	- search direction
  !w 	- work vector -> On output, new iterate
  !g    - residual evaluated at new iterate y

  real(kind=realType) :: fnorm,gnorm,ynorm
  real(kind=realType) :: alpha
  logical :: flag
  integer(kind=intType) :: nfevals
  !   Note that for line search purposes we work with with the related
  !   minimization problem:
  !      min  z(x):  R^n -> R,
  !   where z(x) = .5 * fnorm*fnorm, and fnorm = || f ||_2.
  !         

  real(kind=realType) :: initslope,lambdaprev,gnormprev,a,b,d,t1,t2,rellength
  real(kind=realType) :: minlambda,lambda,lambdatemp

  integer(kind=intType) :: ierr

  ! Set some defaults:
  alpha		= 1.e-2
  minlambda     = 1.e-5
  nfevals = 0
  flag = .True. 

  ! Compute the two norms we need:
  call VecNorm(y,NORM_2,ynorm,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecNorm(f,NORM_2,fnorm,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecMaxPointwiseDivide(y,x,rellength,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  minlambda = minlambda/rellength ! Fix this
  call MatMult(dRdw,y,w,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  nfevals = nfevals + 1

  call VecDot(f,w,initslope,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (initslope > 0.0)  then
     initslope = -initslope
  end if

  if (initslope == 0.0) then
     initslope = -1.0
  end if

  call VecWAXPY(w,-1.0,y,x,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Compute Function:
  call setW(w)
  call computeResidualNK()

  nfevals = nfevals + 1

  call VecNorm(g,NORM_2,gnorm,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Sufficient reduction 
  if (half*gnorm*gnorm <= half*fnorm*fnorm + alpha*initslope) then
     goto 100
  end if

  ! Fit points with quadratic 
  lambda     = 1.0
  lambdatemp = -initslope/(gnorm*gnorm - fnorm*fnorm - 2.0*initslope)
  lambdaprev = lambda
  gnormprev  = gnorm

  if (lambdatemp > half*lambda) then
     lambdatemp = half*lambda
  end if

  if (lambdatemp <= .1*lambda) then
     lambda = .1*lambda
  else                 
     lambda = lambdatemp
  end if

  call VecWAXPY(w,-lambda,y,x,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Compute new function again:
  call setW(w)
  call computeResidualNK()

  call setRVec(g)

  nfevals = nfevals + 1

  call VecNorm(g,NORM_2,gnorm,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Sufficient reduction 
  if (half*gnorm*gnorm <= half*fnorm*fnorm + lambda*alpha*initslope) then
     goto 100
  end if

  ! Fit points with cubic 
  cubic_loop: do while (.True.) 
     
     if (lambda <= minlambda) then 
        flag = .False.
        exit cubic_loop
     end if
     t1 = half*(gnorm*gnorm - fnorm*fnorm) - lambda*initslope
     t2 = half*(gnormprev*gnormprev  - fnorm*fnorm) - lambdaprev*initslope

     a  = (t1/(lambda*lambda) - t2/(lambdaprev*lambdaprev))/(lambda-lambdaprev)
     b  = (-lambdaprev*t1/(lambda*lambda) + lambda*t2/(lambdaprev*lambdaprev))/(lambda-lambdaprev)
     d  = b*b - 3*a*initslope
     if (d < 0.0) then
        d = 0.0
     end if

     if (a == 0.0) then
        lambdatemp = -initslope/(2.0*b)
     else
        lambdatemp = (-b + sqrt(d))/(3.0*a)
     end if

     lambdaprev = lambda
     gnormprev  = gnorm

     if (lambdatemp > half*lambda)  then
        lambdatemp = half*lambda
     end if
     if (lambdatemp <= .1*lambda) then
        lambda = .1*lambda
     else           
        lambda = lambdatemp
     end if

     call VecWAXPY(w,-lambda,y,x,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Compute new function again:
     call setW(w)
     call computeResidualNK()

     call setRVec(g)
     nfevals = nfevals + 1

     call VecNorm(g,NORM_2,gnorm,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Is reduction enough?
     if (.5*gnorm*gnorm <= .5*fnorm*fnorm + lambda*alpha*initslope) then
        exit cubic_loop
     end if
  end do cubic_loop

100 continue

  ! Optional user-defined check for line search step validity */
#endif
end subroutine LSCubic

subroutine LSNone(x,f,g,y,w,nfevals,flag)
#ifndef USE_NO_PETSC
  use precision 
  use communication
  use NKSolverVars, only: dRdw
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Input/Output
  Vec x,f,g,y,w
  !x 	- current iterate
  !f 	- residual evaluated at x
  !y 	- search direction
  !w 	- work vector -> On output, new iterate
  !g    - residual evaluated at new iterate y

  integer(kind=intType) :: nfevals
  integer(kind=intType) :: ierr
  logical :: flag
  flag = .True. 
  ! We just accept the step and compute the new residual at the new iterate
  nfevals = 0
  call VecWAXPY(w,-1.0,y,x,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Compute new function:
  call setW(w)
  call computeResidualNK()
  
  call setRVec(g)
  nfevals = nfevals + 1
#endif
end subroutine LSNone

subroutine LSNM(x,f,g,y,w,fnorm,ynorm,gnorm,nfevals,flag)
#ifndef USE_NO_PETSC
  use precision 
  use communication
  use NKSolverVars, only: dRdw, func_evals, iter_k,iter_m
  use constants
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Input/Output
  Vec x,f,g,y,w
  !x 	- current iterate
  !f 	- residual evaluated at x
  !y 	- search direction
  !w 	- work vector -> On output, new iterate
  !g    - residual evaluated at new iterate y

  real(kind=realType) :: fnorm,gnorm,ynorm
  logical :: flag
  integer(kind=intType) :: nfevals
  !   Note that for line search purposes we work with with the related
  !   minimization problem:
  !      min  z(x):  R^n -> R,
  !   where z(x) = .5 * fnorm*fnorm, and fnorm = || f ||_2.
  !         

  real(kind=realType) :: initslope,c1,c2,gamma,sigma,alpha,max_val
  integer(kind=intType) :: ierr,iter,j

  ! Set some defaults:
  c1 = 1e-5
  c2 = 1e5
  gamma = 1e-3
  sigma = 0.5

  nfevals = 0
  flag = .True. 

  ! Compute the two norms we need:
  call VecNorm(y,NORM_2,ynorm,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecNorm(f,NORM_2,fnorm,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  func_evals(iter_k) = half*fnorm*fnorm

  call MatMult(dRdw,y,w,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  nfevals = nfevals + 1

  call VecDot(f,w,initslope,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (initslope > 0.0)  then
     initslope = -initslope
  end if

  if (initslope == 0.0) then
     initslope = -1.0
  end if

  alpha = 1.0 ! Initial step length:
  backtrack: do iter=1,10

     ! Compute new x value:
     call VecWAXPY(w,-alpha,y,x,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Compute Function @ new x (w is the work vector
     call setW(w)
     call computeResidualNK()
     call setRVec(g)
     nfevals = nfevals + 1
     
     ! Compute the norm at the new trial location
     call VecNorm(g,NORM_2,gnorm,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     max_val = func_evals(iter_k) + alpha*gamma*initSlope
     ! Loop over the previous, m function values and find the max:
     do j=iter_k-1,iter_k-iter_m,-1
        max_val = max(max_val,func_evals(j) + alpha*gamma*initSlope)
     end do
        
     ! Sufficient reduction 
     if (half*gnorm*gnorm <= max_val) then
        exit backtrack
     else
        alpha = alpha * sigma
     end if
  end do backtrack

#endif
end subroutine LSNM

subroutine getEWTol(iter,norm,old_norm,rtol_last,rtol)

  use constants
  implicit none

  ! There are the default EW Parameters from PETSc. They seem to work well
  !version:           2
  !rtol_0:  0.300000000000000     
  !rtol_max:  0.900000000000000     
  !gamma:   1.00000000000000     
  !alpha:   1.61803398874989     
  !alpha2:   1.61803398874989     
  !threshold:  0.100000000000000     

  integer(kind=intType) :: iter
  real(kind=realType), intent(in) :: norm,old_norm,rtol_last
  real(kind=realType), intent(out) :: rtol
  real(kind=realType) :: rtol_0,rtol_max,gamma,alpha,alpha2,threshold,stol

  rtol_0    = 0.30_realType
  rtol_max  = 0.9_realType
  gamma     = 1.0_realType
  alpha     = (one+sqrt(five))/two
  alpha2    = (one+sqrt(five))/two
  threshold = 0.10_realType

  if (iter == 1) then
     rtol = rtol_0
  else
     ! We use version 2:
     rtol = gamma*(norm/old_norm)**alpha
     stol = gamma*rtol_last**alpha

     if (stol > threshold) then
        rtol = max(rtol,stol)
     end if

     ! Safeguard: avoid rtol greater than one
     rtol = min(rtol,rtol_max)
  end if

end subroutine getEWTol
