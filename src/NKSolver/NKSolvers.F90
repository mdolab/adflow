module NKSolver

  use constants
  implicit none

  ! MPI comes from constants, so we need to avoid MPIF_H in PETSc
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"

  ! PETSc Matrices:
  ! dRdw: This is the actual matrix-free matrix computed with FD
  ! dRdwPre: The preconditoner matrix for NK method. This matrix is stored.
  ! dRdwPseudo: Shell matrix used with the pseudo-transient 
  !             continuation method. 

  Mat  dRdw, dRdwPre, dRdwPseudo

  ! PETSc Vectors:
  ! wVec: PETsc version of SUmb 'w'
  ! rVec: PETSc version of SUmb 'dw', but divided by volume
  ! deltaW: Update to the wVec from linear solution
  ! diagV: Diagonal lumping term

  Vec wVec, rVec, deltaW, work, g, wBase

  ! NK_KSP: The ksp object for solving the newton udpate
  KSP  NK_KSP

  PetscFortranAddr   ctx(1)

  ! Options for NK Slver
  logical :: useNKSolver
  integer(kind=intType) :: NK_jacobianLag
  integer(kind=intType) :: NK_subspace
  integer(kind=intType) :: NK_asmOverlap
  integer(kind=intType) :: NK_iluFill
  integer(kind=intType) :: NK_innerPreConIts
  integer(kind=intType) :: NK_outerPreConIts
  integer(kind=intType) :: NK_LS
  logical :: NK_useEW 
  logical :: NK_ADPC
  logical :: NK_viscPC
  real(kind=realType) :: NK_CFL0
  real(kind=realType) :: NK_switchTol
  real(kind=realType) :: NK_rtolInit
  real(kind=realType) :: NK_divTol = 10

  ! Misc variables
  logical :: NK_solverSetup=.False.
  integer(kind=intType) :: NK_iter
  logical :: rkReset
  ! Eisenstat-Walker Parameters
  integer(kind=intType) :: ew_version
  real(kind=realType) :: ew_rtol_0
  real(kind=realType) :: ew_rtol_max
  real(kind=realType) :: ew_gamma
  real(kind=realType) :: ew_alpha
  real(kind=realType) :: ew_alpha2
  real(kind=realType) :: ew_threshold
  real(kind=realType) :: rtolLast, oldNorm

  ! Misc Parameters
  real(kind=realType) :: totalR0, totalRStart, totalRFinal, totalR
  real(kind=realType) :: rhoRes0, rhoResStart, rhoResFinal, rhoRes
  logical :: freeStreamResSet=.False.
  real(kind=realType) :: NK_CFL

  ! Variables for non-monotone line search
  real(kind=realType), dimension(:), allocatable :: NKLSFuncEvals
  integer(kind=intType) :: Mmax=5
  integer(kind=intType) :: iter_k
  integer(kind=intType) :: iter_m

  ! Parameter for external preconditioner
  integer(kind=intType) :: applyPCSubSpaceSize

contains

  subroutine setupNKsolver

    ! Setup the PETSc objects for the Newton-Krylov
    ! solver. destroyNKsolver can be used to destroy the objects created
    ! in this function

    use constants
    use stencils, only : visc_pc_stencil, euler_pc_stencil, N_visc_pc, N_euler_pc
    use communication, only : sumb_comm_world
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowVarRefState, only : nw, viscous
    use InputAdjoint, only: viscPC
    use ADjointVars , only: nCellsLocal
    use utils, only : EChk
    implicit none

    ! Working Variables
    integer(kind=intType) :: ierr, nDimw
    integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
    integer(kind=intType) :: n_stencil
    integer(kind=intType), dimension(:, :), pointer :: stencil
    integer(kind=intType) :: level

    ! Make sure we don't have memory for the approximate and exact
    ! Newton solvers kicking around at the same time.
    !call destroyANKSolver()

    if (.not. NK_solverSetup) then
       nDimW = nw * nCellsLocal(1_intTYpe) * nTimeIntervalsSpectral

       call VecCreate(SUMB_COMM_WORLD, wVec, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecSetSizes(wVec, nDimW, PETSC_DECIDE, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecSetBlockSize(wVec, nw, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecSetType(wVec, VECMPI, ierr) 
       call EChk(ierr, __FILE__, __LINE__)

       !  Create duplicates for residual and delta
       call VecDuplicate(wVec, rVec, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecDuplicate(wVec, deltaW, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Create the two additional work vectors for the line search:
       call VecDuplicate(wVec, g, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecDuplicate(wVec, work, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Create Pre-Conditioning Matrix
       allocate(nnzDiagonal(nCellsLocal(1_intType)*nTimeIntervalsSpectral), &
            nnzOffDiag(nCellsLocal(1_intType)*nTimeIntervalsSpectral) )

       if (viscous .and. NK_viscPC) then
          stencil => visc_pc_stencil
          n_stencil = N_visc_pc
       else
          stencil => euler_pc_stencil
          n_stencil = N_euler_pc
       end if

       level = 1
       call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nw, stencil, n_stencil, &
            level, .False.)
       call myMatCreate(dRdwPre, nw, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
            __FILE__, __LINE__)

       call matSetOption(dRdwPre, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
       call EChk(ierr, __FILE__, __LINE__)
       deallocate(nnzDiagonal, nnzOffDiag)

       ! Setup Matrix-Free dRdw matrix and its function
       call MatCreateMFFD(sumb_comm_world, nDimW, nDimW, &
            PETSC_DETERMINE, PETSC_DETERMINE, dRdw, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call MatMFFDSetFunction(dRdw, FormFunction_mf, ctx, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Setup a matrix free matrix for drdw
       call MatCreateShell(SUMB_COMM_WORLD, nDimW, nDimW, PETSC_DETERMINE, &
            PETSC_DETERMINE, ctx, dRdwPseudo, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set the shell operation for doing matrix vector multiplies
       call MatShellSetOperation(dRdwPseudo, MATOP_MULT, NKMatMult, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set the mat_row_oriented option to false so that dense
       ! subblocks can be passed in in fortran column-oriented format
       call MatSetOption(dRdWPre, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call MatSetOption(dRdW, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       !  Create the linear solver context
       call KSPCreate(SUMB_COMM_WORLD, NK_KSP, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set operators for the solver
       call KSPSetOperators(NK_KSP, dRdw, dRdwPre, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       NK_solverSetup = .True.
       NK_iter = 0
    end if

  end subroutine setupNKsolver

  subroutine NKMatMult(A, vecX,  vecY, ierr)

    ! PETSc user-defied call back function for computing the product of
    ! dRdw with a vector. Here we just call the much more broadly
    ! useful routine computeMatrixFreeProductFwd()

    use constants
    use utils, only : EChk
    implicit none

    ! PETSc Arguments
    Mat   A
    Vec   vecX, vecY
    integer(kind=intType) ::ierr, i, j, k, l, nn, sps, ii
    real(kind=realType) :: dt
    real(kind=realType), pointer :: yPtr(:), xPtr(:)

    ! Frist run the underlying matrix-free mult
    call matMult(dRdw, vecX, vecY, ierr)

    call VecGetArrayF90(vecY, yPtr, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecGetArrayReadF90(vecX, xPtr, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    yPtr = yPtr + one/NK_CFL*xPtr

    call VecRestorearrayF90(vecY, yPtr, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecRestorearrayReadF90(vecX, xPtr, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine NKMatMult

  subroutine getFreeStreamResidual(rhoRes, totalRRes)

    use constants
    use blockPointers, only : nDom, ib, jb, kb, w
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowVarRefState, only : nw, winf
    use utils, only : setPointers
    implicit none

    real(kind=realType), intent(out) :: rhoRes, totalRRes
    real(kind=realType), dimension(:), allocatable :: tmp
    integer(kind=intType) :: nDimW, nDimP, counter
    integer(kind=intType) :: nn, sps, i, j, k, l, n

    ! Get the residual cooresponding to the free-stream on the fine
    ! grid-level --- This saves the current values in W, P, rlv and rev
    ! and restores them when finished. 

    call getInfoSize(n)
    allocate(tmp(n))
    call getInfo(tmp, n)

    ! Set the w-variables to the ones of the uniform flow field.
    spectralLoop4b: do sps=1, nTimeIntervalsSpectral
       domains4b: do nn=1, nDom
          call setPointers(nn, 1, sps)
          do l=1,nw
             do k=0, kb
                do j=0, jb
                   do i=0, ib
                      w(i,j,k,l) = winf(l)
                   enddo
                enddo
             enddo
          end do
       end do domains4b
    end do spectralLoop4b

    ! Evaluate the residual now
    call computeResidualNK()
    call getCurrentResidual(rhoRes, totalRRes)

    ! Put everything back
    call setInfo(tmp, n)

    deallocate(tmp)

  end subroutine getFreeStreamResidual

  subroutine getCurrentResidual(rhoRes,totalRRes)

    use constants
    use communication, only : sumb_comm_world 
    use block, only : nCellGlobal
    use blockPointers, only : nDom
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use iteration, only : currentLevel
    use monitor, only: monLoc, monGlob, nMonSum
    use utils, only : setPointers
    implicit none

    ! Compute rhoRes and totalR. The actual residual must have already
    ! been evaluated

    real(kind=realType), intent(out) :: rhoRes,totalRRes
    integer(kind=intType) :: sps,nn,ierr

    monLoc = zero
    do sps=1, nTimeIntervalsSpectral
       do nn=1, nDom
          call setPointers(nn, currentLevel, sps)
          call sumResiduals(1, 1) ! Sum 1st state res into first mon location
          call sumAllResiduals(2) ! Sum into second mon location
       end do
    end do

    ! This is the same calc as in convergence info, just for rehoRes and
    ! totalR only. 
    call mpi_allreduce(monLoc, monGlob, nMonSum, sumb_real, &
         mpi_sum, SUmb_comm_world, ierr)

    rhoRes = sqrt(monGlob(1)/nCellGlobal(currentLevel))
    totalRRes = sqrt(monGlob(2))

  end subroutine getCurrentResidual

  subroutine FormJacobianNK

    use constants
    use inputADjoint, only : viscPC
    use utils, only : EChk
    implicit none

    ! Local Variables
    character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
    integer(kind=intType) :: ierr
    logical :: useAD, usePC, useTranspose, useObjective, tmp
    integer(kind=intType) :: i, j, k, l, ii, nn, sps
    real(kind=realType) :: dt
    real(kind=realType), pointer :: diag(:)
    interface
       subroutine setupStateResidualMatrix(matrix, useAD, usePC, useTranspose, &
            useObjective, frozenTurb, level, matrixTurb)
         use precision
         implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
         Mat :: matrix
         Mat, optional :: matrixTurb
         ! Input Variables
         logical, intent(in) :: useAD, usePC, useTranspose, useObjective, frozenTurb
         integer(kind=intType), intent(in) :: level
       end subroutine setupStateResidualMatrix
    end interface

    ! Dummy assembly begin/end calls for the matrix-free Matrx
    call MatAssemblyBegin(dRdw, MAT_FINAL_ASSEMBLY, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call MatAssemblyEnd(dRdw, MAT_FINAL_ASSEMBLY, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Assemble the approximate PC (fine leve, level 1)
    useAD = NK_ADPC
    usePC = .True.
    useTranspose = .False.
    useObjective = .False.
    tmp = viscPC ! Save what is in viscPC and set to the NKvarible
    viscPC = NK_viscPC

    call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, &
         useObjective, .False., 1_intType)
    ! Reset saved value
    viscPC = tmp

    call VecGetArrayF90(work, diag, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    diag(:) = one/NK_CFL

    call VecRestoreArrayF90(work, diag, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call MatDiagonalSet(dRdwPre, work, ADD_VALUES, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Setup KSP Options
    preConSide = 'right'
    localPCType = 'ilu'
    kspObjectType = 'gmres'
    globalPCType = 'asm'
    localOrdering = 'rcm'

    ! Setup the KSP using the same code as used for the adjoint
    call setupStandardKSP(NK_KSP, kspObjectType, NK_subSpace, &
         preConSide, globalPCType, NK_asmOverlap, NK_outerPreConIts, localPCType, &
         localOrdering, NK_iluFill, NK_innerPreConIts)

    ! Don't do iterative refinement for the NKSolver.
    call KSPGMRESSetCGSRefinementType(NK_KSP, &
         KSP_GMRES_CGS_REFINE_NEVER, ierr)
    call EChk(ierr, __FILE__, __LINE__)

  end subroutine FormJacobianNK

  subroutine FormFunction_mf(ctx, wVec, rVec, ierr)

    ! This is basically a copy of FormFunction, however it has a
    ! different calling sequence from PETSc. It performs the identical
    ! function. This is used for linear solve application for the
    ! aerostructural system pre-conditioner

    use constants
    implicit none

    ! PETSc Variables
    PetscFortranAddr ctx(*)
    Vec     wVec, rVec
    integer(kind=intType) :: ierr

    ! This is just a shell routine that runs the more broadly useful
    ! computeResidualNK subroutine

    call setW(wVec)
    call computeResidualNK()
    call setRVec(rVec)
    ! We don't check an error here, so just pass back zero
    ierr = 0

  end subroutine FormFunction_mf

  subroutine destroyNKsolver

    ! Destroy all the PETSc objects for the Newton-Krylov
    ! solver. 

    use constants
    use utils, only: EChk
    implicit none
    integer(kind=intType) :: ierr

    if (NK_solverSetup) then 

       call MatDestroy(dRdw, ierr) 
       call EChk(ierr, __FILE__, __LINE__)

       call MatDestroy(dRdwPre, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call MatDestroy(dRdwPseudo, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecDestroy(wVec, ierr)  
       call EChk(ierr, __FILE__, __LINE__)

       call VecDestroy(rVec, ierr) 
       call EChk(ierr, __FILE__, __LINE__)

       call VecDestroy(deltaW, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecDestroy(g, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecDestroy(work, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call KSPDestroy(NK_KSP, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       NK_solverSetup = .False.
    end if

  end subroutine destroyNKsolver

  subroutine NKStep(firstCall)

    use constants
    use flowVarRefState, only : nw
    use inputPhysics, only : equations
    use flowVarRefState, only :  nw, nwf
    use inputIteration, only : L2conv
    use iteration, only : approxTotalIts
    use utils, only : EChk
    use killSignals, only : routineFailed
    implicit none

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
       iter_k = 1
       iter_m = 0
    else
       NK_iter = NK_iter + 1

       ! Increment counter for the nonmonotne line serach
       iter_k = iter_k + 1
    end if

    ! Compute the norm of rVec for use in EW Criteria
    call VecNorm(rVec, NORM_2, norm, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Determine if if we need to form the Preconditioner
    if (mod(NK_iter, NK_jacobianLag) == 0) then
       NK_CFL = NK_CFL0 * (totalR0 / norm)**1.5
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
       iter_m = min(iter_m+1, mMax)
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


  end subroutine NKStep

  subroutine LSCubic(x, f, g, y, w, fnorm, ynorm, gnorm, nfevals, flag)

    use constants
    use utils, only : EChk
    implicit none

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

  end subroutine LSCubic

  subroutine LSNone(x, f, g, y, w, nfevals, flag)

    use constants
    use utils, only : EChk
    implicit none

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

  end subroutine LSNone

  subroutine LSNM(x, f, g, y, w, fnorm, ynorm, gnorm, nfevals, flag)

    use constants
    use utils, only : EChk
    implicit none

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

    NKLSFuncEvals(iter_k) = 0.5_realType*fnorm*fnorm

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

          max_val = NKLSFuncEvals(iter_k) + alpha*gamma*initSlope

          ! Loop over the previous, m function values and find the max:
          do j=iter_k-1, iter_k-iter_m+1, -1
             max_val = max(max_val, NKLSFuncEvals(j) + alpha*gamma*initSlope)
          end do

          ! Sufficient reduction 
          if (0.5_realType*gnorm*gnorm <= max_val) then
             exit backtrack
          else
             alpha = alpha * sigma
          end if
       end if
    end do backtrack
  end subroutine LSNM

  subroutine computeResidualNK()

    ! This is the residual evaluation driver for the NK solver. The
    ! actual function that is used for the matrix free matrix-vector
    ! products is FormFunction_mf (see formFunction.F90). This the
    ! routine that actually computes the residual. This works with
    ! Euler, Laminar and RANS equation modes. 

    ! This function uses the w that is currently stored in the flowDoms
    ! datastructure and leaves the resulting residual dw, in the same
    ! structure. setW() and setRVec() is used in formFunction to
    ! set/extract these values for communication with PETSc. 

    use constants
    use blockPointers, only : nDom, ib, jb, kb, p, w, gamma
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nw, pInfCorr, nwf, kPresent, nt1, nt2
    use iteration, only : groundLevel, currentLevel, rkStage
    use inputPhysics, only : equations, gammaConstant
    use utils, only : setPointers
    use haloExchange, only : whalo2
    use turbUtils, only : computeEddyViscosity
    use turbAPI, only : turbResidual
    use turbBCRoutines, only : applyAllTurbBCThisBLock, bcturbTreatment
    implicit none

    ! Local Variables
    integer(kind=intType) :: i, j, k, sps,nn
    logical secondHalo, correctForK
    real(kind=realType) :: gm1, factK, v2

    gm1 = gammaConstant - one
    rkStage = 0

    secondHalo = .false.
    correctForK = .false.
    if(currentLevel <= groundLevel) then
       secondHalo = .true.
       if (kPresent) then 
          correctForK = .True.
       end if
    end if

    ! Recompute pressure on ALL cells 
    spectralLoop: do sps=1, nTimeIntervalsSpectral
       domainsState: do nn=1, nDom
          ! Set the pointers to this block.
          call setPointers(nn, currentLevel, sps)
          factK = zero
          do k=0, kb
             do j=0, jb
                do i=0, ib

                   gm1  = gamma(i, j, k) - one
                   factK = five*third - gamma(i, j ,k)
                   v2 = w(i,j,k,ivx)**2 + w(i,j,k,ivy)**2 &
                        + w(i,j,k,ivz)**2

                   p(i,j,k) = gm1*(w(i,j,k,irhoE) &
                        - half*w(i,j,k,irho)*v2) 

                   if( correctForK ) then 
                      p(i, j ,K) = p(i,j, k) + factK*w(i, j, k, irho) &
                           * w(i, j, k, itu1)
                   end if

                   ! Clip to make sure it is positive.
                   p(i,j,k) = max(p(i,j,k), 1.e-4_realType*pInfCorr)
                end do
             end do
          end do

          ! Compute Viscosities
          call computeLamViscosity  
          call computeEddyViscosity 
       end do domainsState
    end do spectralLoop

    ! Apply BCs
    call applyAllBC(secondHalo)

    if (equations == RANSequations) then 
       do nn=1,nDom
          do sps=1,nTimeIntervalsSpectral
             call setPointers(nn, currentLevel, sps)
             call bcTurbTreatment
             call applyAllTurbBCThisBLock(.True.)
          end do
       end do
    end if

    ! Exchange halos
    call whalo2(currentLevel, 1_intType, nw, .true., &
         .true., .true.)

    ! Compute time step (spectral radius is actually what we need)
    call timestep(.false.)

    ! Possible Turblent Equations
    if( equations == RANSEquations ) then
       ! Compute the skin-friction velocity (wall functions only)
       call computeUtau
       call initres(nt1, nt2) ! Initialize only the Turblent Variables
       call turbResidual
    endif

    ! Initialize Flow residuals
    call initres(1_intType, nwf)

    ! Actual Residual Calc
    call residual 

  end subroutine computeResidualNK

  subroutine applyPC(in_vec, out_vec, ndof)

    ! Apply the NK PC to the in_vec. This subroutine is ONLY used as a
    ! preconditioner for a global Aero-Structural Newton-Krylov Method

    use constants
    use utils, only : EChk

    implicit none

    ! Input/Output
    integer(kind=intType) :: ndof
    real(kind=realType), dimension(ndof), intent(in)    :: in_vec
    real(kind=realTYpe), dimension(ndof), intent(inout) :: out_vec

    ! Working Variables
    integer(kind=intType) :: ierr

    ! Setup the NKsolver if not already done so
    if (.not. NK_solverSetup) then
       call setupNKSolver
    end if

    ! We possibly need to re-form the jacobian
    if (mod(NK_iter, NK_jacobianLag) == 0) then 
       call FormJacobianNK()
    end if

    ! Place the two arrays into two vectos. We reuse 'work' and 'g'. 
    call VecPlaceArray(work, in_vec, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecPlaceArray(g, out_vec, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set the base vec
    call setwVec(wVec)

    call MatMFFDSetBase(dRdW, wVec, PETSC_NULL_OBJECT, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! This needs to be a bit better...
    call KSPSetTolerances(NK_KSP, 1e-8, 1e-16, 10.0, &
         applyPCSubSpaceSize, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Actually do the Linear Krylov Solve
    call KSPSolve(NK_KSP, work, g, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Reset the array pointers:
    call VecResetArray(work, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecResetArray(g, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    NK_iter = NK_iter + 1

  end subroutine applyPC

  subroutine applyAdjointPC(in_vec, out_vec, ndof)

    ! Apply the Adjoint PC to the in_vec. This subroutine is ONLY used as a
    ! preconditioner for a global Aero-Structural Krylov Method

    use constants
    use ADjointPETSc, only : adjointKSP, KSP_NORM_NONE, PETSC_DEFAULT_REAL, &
         psi_like1, psi_like2
    use inputAdjoint, only : applyAdjointPCSubSpaceSize
    use utils, only : EChk
    implicit none

    ! Input/Output
    integer(kind=intType) :: ndof
    real(kind=realType), dimension(ndof), intent(in)    :: in_vec
    real(kind=realTYpe), dimension(ndof), intent(inout) :: out_vec

    ! Working Variables
    integer(kind=intType) :: ierr

    ! Hijack adjoint and adjointRes with in_vec and out_vec
    call VecPlaceArray(psi_like1, in_vec, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecPlaceArray(psi_like2, out_vec, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set KSP_NORM Type to none. Implictly turns off convergence
    ! check. Since we just want to run a fixed number of iterations this
    ! is fine. The should be set regardless of the KSPType.

    call KSPSetNormType(adjointKSP, KSP_NORM_NONE, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! This needs to be a bit better...
    call KSPSetTolerances(adjointKSP, PETSC_DEFAULT_REAL, &
         PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, &
         applyAdjointPCSubSpaceSize, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Actually do the Linear Krylov Solve
    call KSPSolve(adjointKSP, psi_like1, psi_like2, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Reset the array pointers:
    call VecResetArray(psi_like1, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecResetArray(psi_like2, ierr)
    call EChk(ierr, __FILE__, __LINE__)

  end subroutine applyAdjointPC

  subroutine setWVec(wVec)

    ! Set the current residual in dw into the PETSc Vector

    use constants
    use blockPointers, only : nDom, il, jl, kl, w
    use inputtimespectral, only : ntimeIntervalsSpectral
    use flowvarrefstate, only : nw
    use utils, only : setPointers, EChk
    implicit none

    Vec   wVec
    integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
    real(kind=realType),pointer :: wvec_pointer(:)

    call VecGetArrayF90(wVec,wvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)
    ii = 1
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn,1_intType,sps)
          ! Copy off w to wVec
          do k=2,kl
             do j=2,jl
                do i=2,il
                   do l=1,nw
                      wvec_pointer(ii) = w(i,j,k,l)
                      ii = ii + 1
                   end do
                end do
             end do
          end do
       end do
    end do

    call VecRestoreArrayF90(wVec,wvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine setWVec

  subroutine setRVec(rVec)

    ! Set the current residual in dw into the PETSc Vector
    use constants
    use blockPointers, only : nDom, volRef, il, jl, kl, dw
    use inputtimespectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nw, nwf, nt1, nt2
    use inputIteration, only : turbResScale
    use utils, only : setPointers, EChk
    implicit none

    Vec    rVec
    integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
    real(kind=realType),pointer :: rvec_pointer(:)
    real(Kind=realType) :: ovv

    call VecGetArrayF90(rVec,rvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)
    ii = 1

    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn,1_intType,sps)
          ! Copy off dw/vol to rVec
          do k=2, kl
             do j=2, jl
                do i=2, il
                   ovv = 1/volRef(i, j, k)
                   do l=1,nwf
                      rvec_pointer(ii) = dw(i, j, k, l)*ovv
                      ii = ii + 1
                   end do
                   do l=nt1,nt2
                      rvec_pointer(ii) = dw(i, j, k, l)*ovv*turbResScale(l-nt1+1)
                      ii = ii + 1
                   end do
                end do
             end do
          end do
       end do
    end do

    call VecRestoreArrayF90(rVec,rvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine setRVec

  subroutine setW(wVec)

    use constants
    use blockPointers, only : nDom, il, jl, kl, w
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowVarRefState, only : nw
    use utils, only : setPointers, EChk

    implicit none

    Vec  wVec
    integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
    real(kind=realType),pointer :: wvec_pointer(:)


    call VecGetArrayReadF90(wVec,wvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ii = 1
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn,1_intType,sps)

          do k=2,kl
             do j=2,jl
                do i=2,il
                   do l=1,nw
                      w(i,j,k,l) = wvec_pointer(ii) 
                      ii = ii + 1
                   end do
                end do
             end do
          end do
       end do
    end do

    call VecRestoreArrayReadF90(wVec,wvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine setW

  subroutine getStates(states,ndimw)
    ! Return the state vector, w to Python

    use constants
    use blockPointers, only : il, jl, kl, nDom, w
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nw
    use utils, only : setPointers

    implicit none

    integer(kind=intType),intent(in):: ndimw
    real(kind=realType),dimension(ndimw),intent(out) :: states(ndimw)

    ! Local Variables
    integer(kind=intType) :: nn,i,j,k,l,counter,sps

    counter = 0 
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn,1,sps)
          do k=2,kl
             do j=2,jl
                do i=2,il
                   do l=1,nw
                      counter = counter + 1
                      states(counter) = w(i,j,k,l)
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine getStates

  subroutine getRes(res,ndimw)

    ! Compute the residual and return result to Python
    use constants
    use blockPointers, only : il, jl, kl, nDom, dw, volRef
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nw
    use utils, only : setPointers

    implicit none

    integer(kind=intType),intent(in):: ndimw
    real(kind=realType),dimension(ndimw),intent(inout) :: res(ndimw)

    ! Local Variables
    integer(kind=intType) :: nn,i,j,k,l,counter,sps
    real(kind=realType) :: ovv

    call computeResidualNK()
    counter = 0 
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn,1,sps)
          do k=2,kl
             do j=2,jl
                do i=2,il
                   ovv = one/volRef(i,j,k)
                   do l=1,nw
                      counter = counter + 1
                      res(counter) = dw(i,j,k,l)*ovv
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine getRes

  subroutine setStates(states,ndimw)

    ! Take in externallly generated states and set them in SUmb
    use constants
    use blockPointers, only : il, jl, kl, nDom, w
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nw
    use utils, only : setPointers

    implicit none

    integer(kind=intType),intent(in):: ndimw
    real(kind=realType),dimension(ndimw),intent(in) :: states(ndimw)

    ! Local Variables
    integer(kind=intType) :: nn,i,j,k,l,counter,sps

    counter = 0 
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn,1,sps)
          do k=2,kl
             do j=2,jl
                do i=2,il
                   do l=1,nw
                      counter = counter + 1
                      w(i,j,k,l) = states(counter)
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine setStates

  subroutine getInfoSize(iSize)
    use constants
    use blockPointers, only : ib, jb, kb, nDom
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nw, viscous, eddymodel
    use utils, only : setPointers

    implicit none
    integer(kind=intType), intent(out) :: iSize
    integer(kind=intType) :: nn, sps, nc
    ! Determine the size of a flat array needed to store w, P, ( and
    ! rlv, rev if necessary) with full double halos. 
    iSize = 0
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn,1_intType,sps)
          nc = (kb+1)*(jb+1)*(ib+1)
          iSize = iSize + nc*(nw + 1) ! plus 1 for the P
          if (viscous) then
             iSize = iSize + nc
          end if
          if (eddyModel) then
             iSize = iSize + nc
          end if
       end do
    end do
  end subroutine getInfoSize

  subroutine setInfo(info, iSize)

    use constants
    use blockPointers, only : w, p, ib, jb, kb, rlv, rev, nDom
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nw, viscous, eddymodel
    use utils, only : setPointers
    implicit none

    real(kind=realType), intent(in), dimension(iSize) :: info
    integer(kind=intType), intent(in) :: iSize
    integer(kind=intType) :: nn, counter, i, j, k, l, sps
    ! Determine the size of a flat array needed to store w, P, ( and
    ! rlv, rev if necessary) with full double halos. 
    counter = 0
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn,1,sps)
          do k=0,kb
             do j=0,jb
                do i=0,ib
                   do l=1,nw
                      counter = counter + 1
                      w(i,j,k,l) = info(counter)
                   end do

                   counter = counter + 1
                   P(i,j,k) = info(counter)

                   if (viscous) then
                      counter = counter + 1
                      rlv(i,j,k) = info(counter)
                   end if

                   if (eddyModel) then
                      counter = counter + 1
                      rev(i,j,k) = info(counter)
                   end if
                end do
             end do
          end do
       end do
    end do
  end subroutine setInfo

  subroutine getInfo(info, iSize)

    use constants
    use blockPointers, only : w, p, ib, jb, kb, rlv, rev, nDom
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nw, viscous, eddymodel
    use utils, only : setPointers

    implicit none

    real(kind=realType), intent(out), dimension(iSize) :: info
    integer(kind=intType), intent(in) :: iSize
    integer(kind=intType) ::  nn, counter, i, j, k, l, sps
    ! Determine the size of a flat array needed to store w, P, ( and
    ! rlv, rev if necessary) with full double halos. 
    counter = 0
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn,1,sps)
          do k=0,kb
             do j=0,jb
                do i=0,ib
                   do l=1,nw
                      counter = counter + 1
                      info(counter) = w(i,j,k,l)
                   end do

                   counter = counter + 1
                   info(counter) = P(i,j,k) 

                   if (viscous) then
                      counter = counter + 1
                      info(counter) = rlv(i,j,k) 
                   end if

                   if (eddyModel) then
                      counter = counter + 1
                      info(counter) = rev(i,j,k)
                   end if
                end do
             end do
          end do
       end do
    end do
  end subroutine getInfo

  subroutine getEWTol(norm, old_norm, rtol_last, rtol)

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

    real(kind=realType), intent(in) :: norm, old_norm, rtol_last
    real(kind=realType), intent(out) :: rtol
    real(kind=realType) :: rtol_max, gamma, alpha, alpha2, threshold, stol

    rtol_max  = 0.5_realType
    gamma     = 1.0_realType
    alpha     = (1.0_realType+sqrt(five))/2.0_realType
    alpha2    = (1.0_realType+sqrt(five))/2.0_realType
    threshold = 0.10_realType
    ! We use version 2:
    rtol = gamma*(norm/old_norm)**alpha
    stol = gamma*rtol_last**alpha

    if (stol > threshold) then
       rtol = max(rtol, stol)
    end if

    ! Safeguard: avoid rtol greater than one
    rtol = min(rtol, rtol_max)

  end subroutine getEWTol
end module NKSolver


module ANKSolver

  use constants
  implicit none

  ! MPI comes from constants, so we need to avoid MPIF_H in PETSc
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"

  Mat  dRdwPre
  Vec wVec, rVec, deltaW
  KSP  ANK_KSP

  Mat  dRdwPreTurb
  Vec wVecTurb, rVecTurb, deltaWTurb
  KSP  ANK_KSPTurb

  ! Options for ANK Solver
  logical :: useANKSolver
  integer(kind=intType) :: ANK_jacobianLag
  integer(kind=intType) :: ANK_subSpace
  integer(kind=intType) :: ANK_asmOverlap
  integer(kind=intType) :: ANK_iluFill
  integer(kind=intType) :: ANK_innerPreConIts
  real(kind=realType)   :: ANK_rtol
  real(kind=realType)   :: ANK_switchTol
  real(kind=realType)   :: ANK_divTol = 10
  logical :: ANK_useTurbDADI

  ! Misc variables
  real(kind=realType) :: ANK_CFL, ANK_CFL0
  logical :: ANK_solverSetup=.False.
  logical :: ANK_turbSetup=.False.
  integer(kind=intTYpe) :: ANK_iter

contains

  subroutine setupANKsolver

    ! Setup the PETSc objects for the Newton-Krylov
    ! solver. destroyNKsolver can be used to destroy the objects created
    ! in this function

    use constants
    use stencils, only : euler_PC_stencil, N_euler_PC
    use communication, only : sumb_comm_world
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowVarRefState, only : nw, viscous, nwf
    use ADjointVars , only: nCellsLocal
    use NKSolver, only : destroyNKSolver
    use utils, only : EChk
    implicit none

    ! Working Variables
    integer(kind=intType) :: ierr, nDimw
    integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
    integer(kind=intType) :: n_stencil
    integer(kind=intType), dimension(:, :), pointer :: stencil
    integer(kind=intType) :: level

    ! Make sure we don't have memory for the approximate and exact
    ! Newton solvers kicking around at the same time.
    call destroyNKSolver()

    if (.not. ANK_solverSetup) then
       nDimW = nwf * nCellsLocal(1_intTYpe) * nTimeIntervalsSpectral

       call VecCreate(SUMB_COMM_WORLD, wVec, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecSetSizes(wVec, nDimW, PETSC_DECIDE, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecSetBlockSize(wVec, nwf, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecSetType(wVec, VECMPI, ierr) 
       call EChk(ierr, __FILE__, __LINE__)

       !  Create duplicates for residual and delta
       call VecDuplicate(wVec, rVec, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecDuplicate(wVec, deltaW, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Create Pre-Conditioning Matrix
       allocate(nnzDiagonal(nCellsLocal(1_intType)*nTimeIntervalsSpectral), &
            nnzOffDiag(nCellsLocal(1_intType)*nTimeIntervalsSpectral) )

       stencil => euler_pc_stencil
       n_stencil = N_euler_pc

       level = 1
       call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nwf, stencil, n_stencil, &
            level)
       call myMatCreate(dRdwPre, nwf, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
            __FILE__, __LINE__)

       call matSetOption(dRdwPre, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
       call EChk(ierr, __FILE__, __LINE__)
       deallocate(nnzDiagonal, nnzOffDiag)

       ! Set the mat_row_oriented option to false so that dense
       ! subblocks can be passed in in fortran column-oriented format
       call MatSetOption(dRdWPre, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       !  Create the linear solver context
       call KSPCreate(SUMB_COMM_WORLD, ANK_KSP, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set operators for the solver

       call KSPSetOperators(ANK_KSP, dRdwPre, dRdwPre, ierr)
       call EChk(ierr, __FILE__, __LINE__)


       ! =================== Turbulence Setup =====================
       if (.not. ANK_useTurbDADI .and. nw > nwf) then 
          nDimW = nCellsLocal(1_intTYpe) * nTimeIntervalsSpectral

          call VecCreate(SUMB_COMM_WORLD, wVecTurb, ierr)
          call EChk(ierr, __FILE__, __LINE__)

          call VecSetSizes(wVecTurb, nDimW, PETSC_DECIDE, ierr)
          call EChk(ierr, __FILE__, __LINE__)

          call VecSetBlockSize(wVecTurb, 1, ierr)
          call EChk(ierr, __FILE__, __LINE__)

          call VecSetType(wVecTurb, VECMPI, ierr) 
          call EChk(ierr, __FILE__, __LINE__)

          !  Create duplicates for residual and delta
          call VecDuplicate(wVecTurb, rVecTurb, ierr)
          call EChk(ierr, __FILE__, __LINE__)

          call VecDuplicate(wVecTurb, deltaWTurb, ierr)
          call EChk(ierr, __FILE__, __LINE__)

          ! Create Pre-Conditioning Matrix
          allocate(nnzDiagonal(nCellsLocal(1_intType)*nTimeIntervalsSpectral), &
               nnzOffDiag(nCellsLocal(1_intType)*nTimeIntervalsSpectral) )

          stencil => euler_pc_stencil
          n_stencil = N_euler_pc

          level = 1
          call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW, stencil, n_stencil, &
               level)
          call myMatCreate(dRdwPreTurb, 1, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
               __FILE__, __LINE__)

          call matSetOption(dRdwPreTurb, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
          call EChk(ierr, __FILE__, __LINE__)
          deallocate(nnzDiagonal, nnzOffDiag)

          ! Set the mat_row_oriented option to false so that dense
          ! subblocks can be passed in in fortran column-oriented format
          call MatSetOption(dRdWPreTurb, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
          call EChk(ierr, __FILE__, __LINE__)

          !  Create the linear solver context
          call KSPCreate(SUMB_COMM_WORLD, ANK_KSPTurb, ierr)
          call EChk(ierr, __FILE__, __LINE__)

          ! Set operators for the solver
          call KSPSetOperators(ANK_KSPTurb, dRdwPreTurb, dRdwPreTurb, ierr)
          call EChk(ierr, __FILE__, __LINE__)

          ANK_turbSetup = .True.
       else
          ANK_turbSetup = .False.
       end if

       ANK_solverSetup = .True.
       ANK_iter = 0
    end if

  end subroutine setupANKsolver


  subroutine FormJacobianANK

    use constants
    use flowVarRefState, only : nw, nwf
    use inputADjoint, only : viscPC
    use utils, only : EChk
    implicit none

    ! Local Variables
    character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
    integer(kind=intType) ::ierr
    logical :: useAD, usePC, useTranspose, useObjective, tmp
    real(kind=realType) ::  dt
    integer(kind=intType) :: i, j, k, l, ii, nn, sps, outerPreConIts
    real(kind=realType), pointer :: diag(:)
    external :: myKSPMonitor
    interface
       subroutine setupStateResidualMatrix(matrix, useAD, usePC, useTranspose, &
            useObjective, frozenTurb, level, matrixTurb)
         use precision
         implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
         Mat :: matrix
         Mat, optional :: matrixTurb
         ! Input Variables
         logical, intent(in) :: useAD, usePC, useTranspose, useObjective, frozenTurb
         integer(kind=intType), intent(in) :: level
       end subroutine setupStateResidualMatrix
    end interface

    ! Assemble the approximate PC (fine leve, level 1)
    useAD = .False.
    usePC = .True.
    useTranspose = .False.
    useObjective = .False.
    tmp = viscPC ! Save what is in viscPC and set to the NKvarible
    viscPC = .False.
    if (ANK_useTurbDADI) then 
       call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, &
            useObjective, .True., 1_intType)
    else
       ! The turbulence jacobian will only be accuate with AD.
       useAD = .True.
       call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, &
            useObjective, .False., 1_intType, dRdwPreTurb)
    end if
    ! Reset saved value
    viscPC = tmp

    ! ----------- Setup Flow KSP ----------
    call VecGetArrayF90(deltaW, diag, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    diag(:) = one/ANK_CFL

    call VecRestoreArrayF90(deltaW, diag, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call MatDiagonalSet(dRdwPre, deltaW, ADD_VALUES, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Setup KSP Options
    preConSide = 'right'
    localPCType = 'ilu'
    kspObjectType = 'gmres'
    globalPCType = 'asm'
    localOrdering = 'rcm'
    outerPreConIts = 1 
    ! Setup the KSP using the same code as used for the adjoint
    call setupStandardKSP(ANK_KSP, kspObjectType,  ANK_subSpace, &
         preConSide, globalPCType, ANK_asmOverlap, outerPreConIts, localPCType, &
         localOrdering, ANK_iluFill, ANK_innerPreConIts)

    ! Don't do iterative refinement for the NKSolver.
    call KSPGMRESSetCGSRefinementType(ANK_KSP, &
         KSP_GMRES_CGS_REFINE_NEVER, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! ----------- Setup Turb KSP ----------
    if (.not. ANK_useTurbDADI .and. nw > nwf) then 
       call VecGetArrayF90(deltaWTurb, diag, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       diag(:) = one/ANK_CFL

       call VecRestoreArrayF90(deltaWTurb, diag, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call MatDiagonalSet(dRdwPreTurb, deltaWTurb, ADD_VALUES, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Setup KSP Options
       preConSide = 'right'
       localPCType = 'ilu'
       kspObjectType = 'gmres'
       globalPCType = 'asm'
       localOrdering = 'rcm'
       outerPreConIts = 1 
       ! Setup the KSP using the same code as used for the adjoint
       call setupStandardKSP(ANK_KSPTurb, kspObjectType,  ANK_SubSpace, &
            preConSide, globalPCType, ANK_asmOverlap, outerPreConIts, localPCType, &
            localOrdering, ANK_iluFill*2, ANK_innerPreConIts)

       ! Don't do iterative refinement for the NKSolver.
       call KSPGMRESSetCGSRefinementType(ANK_KSPTurb, &
            KSP_GMRES_CGS_REFINE_NEVER, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end if

  end subroutine FormJacobianANK

  subroutine destroyANKsolver

    ! Destroy all the PETSc objects for the Newton-Krylov
    ! solver. 

    use constants
    use utils, only : EChk
    implicit none
    integer(kind=intType) :: ierr

    if (ANK_SolverSetup) then

       call MatDestroy(dRdwPre, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecDestroy(wVec, ierr)  
       call EChk(ierr, __FILE__, __LINE__)

       call VecDestroy(rVec, ierr) 
       call EChk(ierr, __FILE__, __LINE__)

       call VecDestroy(deltaW, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call KSPDestroy(ANK_KSP, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       if (ANK_turbSetup) then 
          call MatDestroy(dRdwPreTurb, ierr)
          call EChk(ierr, __FILE__, __LINE__)

          call VecDestroy(wVecTurb, ierr)  
          call EChk(ierr, __FILE__, __LINE__)

          call VecDestroy(rVecTurb, ierr) 
          call EChk(ierr, __FILE__, __LINE__)

          call VecDestroy(deltaWTurb, ierr)
          call EChk(ierr, __FILE__, __LINE__)

          call KSPDestroy(ANK_KSPTurb, ierr)
          call EChk(ierr, __FILE__, __LINE__)
       end if
       ANK_SolverSetup = .False.
       ANK_TurbSetup = .False.
    end if
  end subroutine destroyANKsolver

  subroutine computeResidualANK()

    ! This is the residual evaluation driver for the ANK solver. It
    ! computes the residual for the mean flow but does not compute the
    ! turbulent residuals. 
    use constants
    use blockPointers
    use inputTimeSpectral
    use flowvarrefstate
    use iteration
    use inputPhysics 
    use utils, only : setPointers
    use haloExchange, only : whalo2
    use turbUtils, only : computeEddyViscosity
    implicit none

    ! Local Variables
    integer(kind=intType) :: i, j, k, sps,nn
    logical secondHalo, correctForK
    real(kind=realType) :: gm1, factK, v2

    gm1 = gammaConstant - one
    rkStage = 0

    secondHalo = .false.
    correctForK = .false.
    if(currentLevel <= groundLevel) then
       secondHalo = .true.
       if (kPresent) then 
          correctForK = .True.
       end if
    end if

    ! Recompute pressure on ALL cells 
    spectralLoop: do sps=1, nTimeIntervalsSpectral
       domainsState: do nn=1, nDom
          ! Set the pointers to this block.
          call setPointers(nn, currentLevel, sps)
          factK = zero
          do k=0, kb
             do j=0, jb
                do i=0, ib

                   gm1  = gamma(i, j, k) - one
                   factK = five*third - gamma(i, j ,k)
                   v2 = w(i,j,k,ivx)**2 + w(i,j,k,ivy)**2 &
                        + w(i,j,k,ivz)**2

                   p(i,j,k) = gm1*(w(i,j,k,irhoE) &
                        - half*w(i,j,k,irho)*v2) 

                   if( correctForK ) then 
                      p(i, j ,K) = p(i,j, k) + factK*w(i, j, k, irho) &
                           * w(i, j, k, itu1)
                   end if

                   ! Clip to make sure it is positive.
                   p(i,j,k) = max(p(i,j,k), 1.e-4_realType*pInfCorr)
                end do
             end do
          end do

          ! Compute Viscosities
          call computeLamViscosity  
          call computeEddyViscosity 
       end do domainsState
    end do spectralLoop

    ! Apply BCs
    call applyAllBC(secondHalo)

    ! Exchange halos
    call whalo2(currentLevel, 1_intType, nwf, .true., &
         .true., .true.)

    ! Compute time step (spectral radius is actually what we need)
    call timestep(.false.)

    ! Initialize Flow residuals
    call initres(1_intType, nwf)

    ! Actual Residual Calc
    call residual 

  end subroutine computeResidualANK


  ! subroutine computeResidualANKTurb()

  !   ! This is the residual evaluation driver for the ANK solver. It
  !   ! computes the residual for the mean flow but does not compute the
  !   ! turbulent residuals. 

  !   use blockPointers
  !   use inputTimeSpectral
  !   use flowvarrefstate
  !   use iteration
  !   use inputPhysics 
  !   implicit none

  !   ! Local Variables
  !   integer(kind=intType) :: nn, sps


  !   call whalo2(currentLevel, nt1, nt2, .False., .False., .False.)
  !   spectralLoop: do sps=1, nTimeIntervalsSpectral
  !      domainsState: do nn=1, nDom
  !         ! Set the pointers to this block.
  !         call setPointers(nn, currentLevel, sps)
  !         call computeEddyViscosity 
  !      end do domainsState
  !   end do spectralLoop
  !   if (equations == RANSEquations) then 

  !      call computeUTau
  !      call initRes(nt1, nt2)
  !      call turbResidual
  !   end if

  ! end subroutine computeResidualANKTurb

  ! This files contains several utilitiy functions that are used with
  ! the ANK solver.

  subroutine setWVecANK(wVec)
    ! Set the current FLOW variables in the PETSc Vector

    use constants
    use blockPointers, only : nDom, il, jl, kl, w
    use inputtimespectral, only : ntimeIntervalsSpectral
    use flowvarrefstate, only : nwf
    use utils, only : setPointers, EChk
    implicit none

    Vec   wVec
    integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
    real(kind=realType),pointer :: wvec_pointer(:)

    call VecGetArrayF90(wVec,wvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)
    ii = 0
    do nn=1, nDom
       do sps=1, nTimeIntervalsSpectral
          call setPointers(nn, 1_intType, sps)
          ! Copy off w to wVec
          do k=2, kl
             do j=2, jl
                do i=2, il
                   do l=1, nwf
                      ii = ii + 1
                      wvec_pointer(ii) = w(i, j, k, l)
                   end do
                end do
             end do
          end do
       end do
    end do

    call VecRestoreArrayF90(wVec, wvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine setWVecANK

  subroutine setRVecANK(rVec)

    ! Set the current FLOW residual in dw into the PETSc Vector
    use constants
    use blockPointers, only : nDom, volRef, il, jl, kl, dw 
    use inputtimespectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nwf
    use utils, only : setPointers, EChk
    implicit none
    Vec    rVec
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
    real(kind=realType),pointer :: rvec_pointer(:)
    real(Kind=realType) :: ovv
    call VecGetArrayF90(rVec,rvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)
    ii = 0
    do nn=1, nDom
       do sps=1, nTimeIntervalsSpectral
          call setPointers(nn,1_intType,sps)
          ! Copy off dw/vol to rVec
          do k=2, kl
             do j=2, jl
                do i=2, il
                   ovv = one/volRef(i,j,k)
                   do l=1, nwf
                      ii = ii + 1        
                      rvec_pointer(ii) = dw(i, j, k, l)*ovv
                   end do
                end do
             end do
          end do
       end do
    end do

    call VecRestoreArrayF90(rVec, rvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine setRVecANK

  subroutine setWANK(wVec)

    use constants
    use blockPointers, only : nDom, vol, il, jl, kl, w
    use inputtimespectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nwf
    use utils, only : setPointers, EChk
    implicit none

    Vec  wVec
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
    real(kind=realType), pointer :: wvec_pointer(:)
    call VecGetArrayReadF90(wVec, wvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ii = 0
    do nn=1, nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn, 1_intType, sps)

          do k=2, kl
             do j=2, jl
                do i=2, il
                   do l=1, nwf
                      ii = ii + 1
                      w(i, j, k, l) = wvec_pointer(ii) 
                   end do
                end do
             end do
          end do
       end do
    end do
    call VecRestoreArrayReadF90(wVec, wvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine setWANK


  subroutine setWVecANKTurb(wVec)
    ! Set the current Turbulence variables in the PETSc Vector
    use constants
    use blockPointers, only : nDom, vol, il, jl, kl, w
    use inputtimespectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nt1, nt2
    use utils, only : setPointers, EChk
    implicit none

    Vec   wVec
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
    real(kind=realType),pointer :: wvec_pointer(:)

    call VecGetArrayF90(wVec, wvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)
    ii = 0
    do nn=1, nDom
       do sps=1, nTimeIntervalsSpectral
          call setPointers(nn, 1_intType, sps)
          ! Copy off w to wVec
          do k=2, kl
             do j=2, jl
                do i=2, il
                   do l=nt1, nt2
                      ii = ii + 1
                      wvec_pointer(ii) = w(i, j, k, l)
                   end do
                end do
             end do
          end do
       end do
    end do

    call VecRestoreArrayF90(wVec, wvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine setWVecANKTurb

  subroutine setRVecANKTurb(rVec)

    ! Set the current FLOW residual in dw into the PETSc Vector
    use constants
    use blockPointers, only : nDom, vol, il, jl, kl, dw, volRef
    use inputtimespectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nt1, nt2
    use inputIteration, only : turbResScale
    use utils, only : setPointers, EChk
    implicit none

    Vec    rVec
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
    real(kind=realType),pointer :: rvec_pointer(:)

    call VecGetArrayF90(rVec,rvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)
    ii = 0
    do nn=1, nDom
       do sps=1, nTimeIntervalsSpectral
          call setPointers(nn, 1_intType, sps)
          ! Copy off dw/vol to rVec
          do k=2, kl
             do j=2, jl
                do i=2, il
                   do l=nt1, nt2
                      ii = ii + 1
                      rvec_pointer(ii) = dw(i,j,k,l)*turbResScale(l-nt1+1)/volRef(i, j, k)
                   end do
                end do
             end do
          end do
       end do
    end do

    call VecRestoreArrayF90(rVec, rvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine setRVecANKTurb

  subroutine setWANKTurb(wVec)
    use constants
    use blockPointers, only : nDom, vol, il, jl, kl, w
    use inputtimespectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nt1, nt2
    use utils, only : setPointers, EChk
    implicit none

    Vec  wVec
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
    real(kind=realType), pointer :: wvec_pointer(:)

    call VecGetArrayReadF90(wVec, wvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ii = 0
    do nn=1, nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn, 1_intType, sps)

          do k=2, kl
             do j=2, jl
                do i=2, il
                   do l=nt1, nt2
                      ii = ii + 1
                      w(i, j, k, l) = wvec_pointer(ii) 
                   end do
                end do
             end do
          end do
       end do
    end do

    call VecRestoreArrayReadF90(wVec, wvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine setWANKTurb

  subroutine ANKStep(firstCall)

    use constants
    use flowVarRefState, only : nw
    use NKSolver, only : totalR0, computeResidualNK
    use inputPhysics, only : equations
    use flowVarRefState, only :  nw, nwf
    use inputIteration, only : L2conv
    use iteration, only : approxTotalIts
    use utils, only : EChk
    use turbAPI, only : turbSolveSegregated
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
end module ANKSolver
