
subroutine LSCubic(x, f, g, y, w, fnorm, ynorm, gnorm, nfevals, flag)
#ifndef USE_NO_PETSC
  use constants
  use communication
  use NKSolverVars, only: dRdw
  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif


  ! Input/Output
  Vec x, f, g, y, w
  !x 	- current iterate
  !f 	- residual evaluated at x
  !y 	- search direction
  !w 	- work vector -> On output, new iterate
  !g    - residual evaluated at new iterate y

  real(kind=realType) :: fnorm, gnorm, ynorm
  real(kind=realType) :: alpha
  logical :: flag
  integer(kind=intType) :: nfevals
  !   Note that for line search purposes we work with with the related
  !   minimization problem:
  !      min  z(x):  R^n -> R, 
  !   where z(x) = .5 * fnorm*fnorm, and fnorm = || f ||_2.
  !         

  real(kind=realType) :: initslope, lambdaprev, gnormprev, a, b, d, t1, t2
  real(kind=realType) :: minlambda, lambda, lambdatemp, rellength
#ifdef USE_COMPLEX
  complex(kind=realType) :: cinitslope
#endif
  integer(kind=intType) :: ierr, iter

  ! Set some defaults:
  alpha		= 1.e-2_realType
  minlambda     = 1.e-7_realType
  nfevals = 0
  flag = .True. 

  ! Compute the two norms we need:
  call VecNorm(y, NORM_2, ynorm, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecNorm(f, NORM_2, fnorm, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecMaxPointwiseDivide(y, x, rellength, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  minlambda = minlambda/rellength ! Fix this
  call MatMult(dRdw, y, w, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  nfevals = nfevals + 1

#ifdef USE_COMPLEX
  call VecDot(f, w, cinitslope, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  initslope = real(cinitslope)
#else
  call VecDot(f, w, initslope, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif

  if (initslope > 0.0_realType)  then
     initslope = -initslope
  end if

  if (initslope == 0.0_realType) then
     initslope = -1.0_realType
  end if
#ifdef USE_COMPLEX
  call VecWAXPY(w, cmplx(-1.0, 0.0), y, x, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#else
  call VecWAXPY(w, -one, y, x, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif

  ! Compute Function:
  call setW(w)
  call computeResidualNK()
  call setRVec(g)  

  nfevals = nfevals + 1

  call VecNorm(g, NORM_2, gnorm, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  if (isnan(gnorm)) then 
     ! Special testing for nans
     lambda = 0.1
     backtrack: do iter=1, 10
        ! Compute new x value:
#ifdef USE_COMPLEX
        call VecWAXPY(w, cmplx(-lambda, 0.0), y, x, ierr)
        call EChk(ierr, __FILE__, __LINE__)
#else
        call VecWAXPY(w, -lambda, y, x, ierr)
        call EChk(ierr, __FILE__, __LINE__)
#endif
        call VecNorm(w, NORM_2, gnorm, ierr)

        ! Compute Function @ new x (w is the work vector
        call setW(w)
        call computeResidualNK()
        call setRVec(g)
        nfevals = nfevals + 1
        
        ! Compute the norm at the new trial location
        call VecNorm(g, NORM_2, gnorm, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        
        if (isnan(gnorm)) then 
           ! Just apply the step limit and keep going (back to the loop start)
           lambda = lambda * .1
        else
           ! We don't care what the value is...its screwed anyway
           exit backtrack
        end if
     end do backtrack
     return
  end if

  ! Sufficient reduction 
  if (0.5_realType*gnorm*gnorm <= 0.5_realType*fnorm*fnorm + alpha*initslope) then
     goto 100
  end if

  ! Fit points with quadratic 
  lambda     = 1.0_realType
  lambdatemp = -initslope/(gnorm*gnorm - fnorm*fnorm - 2.0_realType*initslope)
  lambdaprev = lambda
  gnormprev  = gnorm
  if (lambdatemp > 0.5_realType*lambda) then
     lambdatemp = 0.5_realType*lambda
  end if
  
  if (lambdatemp <= .1_realType*lambda) then
     lambda = .1_realType*lambda
  else                 
     lambda = lambdatemp
  end if

#ifdef USE_COMPLEX
  call VecWAXPY(w, -cmplx(lambda, 0.0), y, x, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#else
  call VecWAXPY(w, -lambda, y, x, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif

  ! Compute new function again:
  call setW(w)
  call computeResidualNK()
  call setRVec(g)

  nfevals = nfevals + 1

  call VecNorm(g, NORM_2, gnorm, ierr)
  if (ierr == PETSC_ERR_FP) then
     flag = .False.
     return
  end if
  call EChk(ierr, __FILE__, __LINE__)

  ! Sufficient reduction 
    if (0.5_realType*gnorm*gnorm <= 0.5_realType*fnorm*fnorm + lambda*alpha*initslope) then
     goto 100
  end if

  ! Fit points with cubic 
  cubic_loop: do while (.True.) 

     if (lambda <= minlambda) then 
        flag = .False.
        exit cubic_loop
     end if
     t1 = 0.5_realType*(gnorm*gnorm - fnorm*fnorm) - lambda*initslope
     t2 = 0.5_realType*(gnormprev*gnormprev  - fnorm*fnorm) - lambdaprev*initslope

     a  = (t1/(lambda*lambda) - t2/(lambdaprev*lambdaprev))/(lambda-lambdaprev)
     b  = (-lambdaprev*t1/(lambda*lambda) + lambda*t2/(lambdaprev*lambdaprev))/(lambda-lambdaprev)
     d  = b*b - three*a*initslope
     if (d < 0.0_realType) then
        d = 0.0_realType
     end if

     if (a == 0.0_realType) then
        lambdatemp = -initslope/(2.0_realType*b)
     else
        lambdatemp = (-b + sqrt(d))/(3.0_realType*a)
     end if

     lambdaprev = lambda
     gnormprev  = gnorm

     if (lambdatemp > 0.5_realType*lambda)  then
        lambdatemp = 0.5_realType*lambda
     end if
     if (lambdatemp <= .1_realType*lambda) then
        lambda = .1_realType*lambda
     else           
        lambda = lambdatemp
     end if

     if (isnan(lambda)) then
        flag = .False.
        exit cubic_loop
     end if

#ifdef USE_COMPLEX
     call VecWAXPY(w, cmplx(-lambda, 0.0), y, x, ierr)
     call EChk(ierr, __FILE__, __LINE__)
#else
     call VecWAXPY(w, -lambda, y, x, ierr)
     call EChk(ierr, __FILE__, __LINE__)
#endif
     ! Compute new function again:
     call setW(w)
     call computeResidualNK()
     call setRVec(g)
     nfevals = nfevals + 1

     call VecNorm(g, NORM_2, gnorm, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Is reduction enough?
     if (0.5_realType*gnorm*gnorm <= 0.5_realType*fnorm*fnorm + lambda*alpha*initslope) then
        exit cubic_loop
     end if
  end do cubic_loop

100 continue
#endif
end subroutine LSCubic

subroutine LSNone(x, f, g, y, w, nfevals, flag)
#ifndef USE_NO_PETSC
  use constants
  implicit none

#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif


  ! Input/Output
  Vec x, f, g, y, w
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
  call VecWAXPY(w, -1.0_realType, y, x, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Compute new function:
  call setW(w)
  call computeResidualNK()
  call setRVec(g)
  
  nfevals = nfevals + 1

#endif
end subroutine LSNone

subroutine LSNM(x, f, g, y, w, fnorm, ynorm, gnorm, nfevals, flag)
#ifndef USE_NO_PETSC
  use precision 
  use communication
  use NKSolverVars, only: dRdw, func_evals, iter_k, iter_m
  use constants
  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif


  ! Input/Output
  Vec x, f, g, y, w
  !x 	- current iterate
  !f 	- residual evaluated at x
  !y 	- search direction
  !w 	- work vector -> On output, new iterate
  !g    - residual evaluated at new iterate y

  real(kind=realType) :: fnorm, gnorm, ynorm
  real(kind=realType) :: alpha
  logical :: flag
  integer(kind=intType) :: nfevals
  !   Note that for line search purposes we work with with the related
  !   minimization problem:
  !      min  z(x):  R^n -> R, 
  !   where z(x) = .5 * fnorm*fnorm, and fnorm = || f ||_2.
  !         
#ifdef USE_COMPLEX
  complex(kind=realType) :: cinitslope
#endif
  real(kind=realType) :: initslope, gamma, sigma,  max_val
  integer(kind=intType) :: ierr, iter, j

  ! Set some defaults:
  gamma = 1e-3_realType
  sigma = 0.5_realType

  nfevals = 0
  flag = .True. 

  ! Compute the two norms we need:
  call VecNorm(y, NORM_2, ynorm, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecNorm(f, NORM_2, fnorm, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  func_evals(iter_k) = 0.5_realType*fnorm*fnorm

  call MatMult(dRdw, y, w, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  nfevals = nfevals + 1

#ifdef USE_COMPLEX
  call VecDot(f, w, cinitslope, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  initslope = real(cinitslope)
#else
  call VecDot(f, w, initslope, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif

  if (initslope > 0.0_realType)  then
     initslope = -initslope
  end if

  if (initslope == 0.0_realType) then
     initslope = -1.0_realType
  end if

  alpha = 1.0 ! Initial step length:
  backtrack: do iter=1, 10

     ! Compute new x value:
#ifdef USE_COMPLEX
     call VecWAXPY(w, cmplx(-alpha, 0.0), y, x, ierr)
     call EChk(ierr, __FILE__, __LINE__)
#else
     call VecWAXPY(w, -alpha, y, x, ierr)
     call EChk(ierr, __FILE__, __LINE__)
#endif

     ! Compute Function @ new x (w is the work vector
     call setW(w)
     call computeResidualNK()
     call setRVec(g)
     nfevals = nfevals + 1
     
     ! Compute the norm at the new trial location
     call VecNorm(g, NORM_2, gnorm, ierr)
     if (ierr == PETSC_ERR_FP) then ! Error code 72 floating point error
        ! Just apply the step limit and keep going (back to the loop start)
        alpha = alpha * sigma
     else
        call EChk(ierr, __FILE__, __LINE__)

        max_val = func_evals(iter_k) + alpha*gamma*initSlope

        ! Loop over the previous, m function values and find the max:
        do j=iter_k-1, iter_k-iter_m+1, -1
           max_val = max(max_val, func_evals(j) + alpha*gamma*initSlope)
        end do
        
        ! Sufficient reduction 
        if (0.5_realType*gnorm*gnorm <= max_val) then
           exit backtrack
        else
           alpha = alpha * sigma
        end if
     end if
  end do backtrack

#endif
end subroutine LSNM
