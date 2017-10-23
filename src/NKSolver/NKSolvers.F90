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
  ! wVec: PETsc version of ADflow 'w'
  ! rVec: PETSc version of ADflow 'dw', but divided by volume
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
  real(kind=realType) :: NK_fixedStep

  ! Misc variables
  logical :: NK_solverSetup=.False.
  integer(kind=intType) :: NK_iter

  ! Eisenstat-Walker Parameters
  integer(kind=intType) :: ew_version
  real(kind=realType) :: ew_rtol_0
  real(kind=realType) :: ew_rtol_max
  real(kind=realType) :: ew_gamma
  real(kind=realType) :: ew_alpha
  real(kind=realType) :: ew_alpha2
  real(kind=realType) :: ew_threshold
  real(kind=alwaysRealType) :: rtolLast, oldNorm

  ! Misc Parameters
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
    use communication, only : adflow_comm_world
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowVarRefState, only : nw, viscous
    use InputAdjoint, only: viscPC
    use ADjointVars , only: nCellsLocal
    use utils, only : EChk
    use adjointUtils, only : myMatCreate, statePreAllocation
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

       call VecCreate(ADFLOW_COMM_WORLD, wVec, ierr)
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
       call MatCreateMFFD(adflow_comm_world, nDimW, nDimW, &
            PETSC_DETERMINE, PETSC_DETERMINE, dRdw, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call MatMFFDSetFunction(dRdw, FormFunction_mf, ctx, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Setup a matrix free matrix for drdw
       call MatCreateShell(ADFLOW_COMM_WORLD, nDimW, nDimW, PETSC_DETERMINE, &
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
       call KSPCreate(ADFLOW_COMM_WORLD, NK_KSP, ierr)
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
    use communication, only : adflow_comm_world
    use block, only : nCellGlobal
    use blockPointers, only : nDom
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use iteration, only : currentLevel
    use monitor, only: monLoc, monGlob, nMonSum
    use utils, only : setPointers, sumResiduals, sumAllResiduals
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
    call mpi_allreduce(monLoc, monGlob, nMonSum, adflow_real, &
         mpi_sum, ADflow_comm_world, ierr)

    rhoRes = sqrt(monGlob(1)/nCellGlobal(currentLevel))
    totalRRes = sqrt(monGlob(2))

  end subroutine getCurrentResidual

  subroutine FormJacobianNK

    use constants
    use inputADjoint, only : viscPC
    use utils, only : EChk
    use adjointUtils, only :setupStateResidualMatrix, setupStandardKSP
    implicit none

    ! Local Variables
    character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
    integer(kind=intType) :: ierr
    logical :: useAD, usePC, useTranspose, useObjective, tmp
    integer(kind=intType) :: i, j, k, l, ii, nn, sps
    real(kind=realType) :: dt
    real(kind=realType), pointer :: diag(:)


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
    use iteration, only : approxTotalIts, totalR0, stepMonitor, linResMonitor
    use utils, only : EChk
    use killSignals, only : routineFailed
    implicit none

    ! Input Variables
    logical, intent(in) :: firstCall

    ! Working Variables
    integer(kind=intType) :: iter, ierr, kspIterations
    integer(kind=intType) :: maxNonLinearIts, nfevals, maxIt
    real(kind=alwaysRealType) :: norm, rtol, atol
    real(kind=alwaysrealType) :: fnorm, ynorm, gnorm
    logical :: flag
    real(kind=alwaysRealType) :: resHist(NK_subspace+1)

    if (firstCall) then
       call setupNKSolver()

       ! Copy the adflow 'w' into the petsc wVec
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

    if (NK_iter == 0 .or. .not. NK_useEW) then
       rtol = NK_rtolInit
    else
       call getEWTol(norm, oldNorm, rtolLast, rtol)
    end if

    ! Save the old rtol and norm for the next iteration
    oldNorm = norm
    rtolLast = rtol

    ! Set all tolerances for linear solver.

    ! The 0.01 requires some explaination: The linear residual is
    ! roughly the same magnitude as the non-linear one. However, it
    ! very rare situations, it can happen that the non-linear residual
    ! is *just* above the convergence criteria, while the linear
    ! residual is *just* below. What happens is that the linear sover
    ! kicks up and doesnt' do anything and then the non-linear
    ! convergnce check can't do anything either. By multiplying by
    ! 0.5, we make sure that the linear solver actually has to do
    ! *something* and not just kick out immediately.
    atol = totalR0*L2Conv*0.01_realType
    maxIt = NK_subspace

    call KSPSetTolerances(NK_KSP, real(rtol), &
         real(atol), real(NK_divTol), maxIt, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call KSPSetResidualHistory(NK_KSP, resHist, maxIt+1, PETSC_TRUE, ierr)
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
       call LSNone(wVec, rVec, g, deltaW, work, nfevals, flag, stepMonitor)
    else if(NK_LS == cubicLineSearch) then
       call LSCubic(wVec, rVec, g, deltaW, work, fnorm, ynorm, gnorm, &
            nfevals, flag, stepMonitor)
    else if (NK_LS == nonMonotoneLineSearch) then
       iter_m = min(iter_m+1, mMax)
       call LSNM(wVec, rVec, g, deltaW, work, fnorm, ynorm, gnorm, &
            nfevals, flag, stepMonitor)
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

    linResMonitor = resHist(kspIterations+1)/resHist(1)

    approxTotalIts = approxTotalIts + nfEvals + kspIterations

  end subroutine NKStep

  subroutine LSCubic(x, f, g, y, w, fnorm, ynorm, gnorm, nfevals, flag, lambda)

    use constants
    use utils, only : EChk
    use communication, only : myid
    use initializeFlow, only : setUniformFlow
    use iteration, only : totalR0
    implicit none

    ! Input/Output
    Vec x, f, g, y, w
    !x 	- current iterate
    !f 	- residual evaluated at x
    !y 	- search direction
    !w 	- work vector -> On output, new iterate
    !g    - residual evaluated at new iterate y

    real(kind=alwaysrealType) :: fnorm, gnorm, ynorm
    real(kind=realType) :: alpha
    logical :: flag
    integer(kind=intType) :: nfevals
    !   Note that for line search purposes we work with with the related
    !   minimization problem:
    !      min  z(x):  R^n -> R,
    !   where z(x) = .5 * fnorm*fnorm, and fnorm = || f ||_2.
    !

    real(kind=realType) :: initslope, lambdaprev, gnormprev, a, b, d, t1, t2
    real(kind=alwaysRealType) :: minlambda, lambda, lambdatemp
    real(kind=alwaysRealType) :: rellength
    integer(kind=intType) :: ierr, iter
    real(kind=alwaysRealType) :: turbRes1, turbRes2, flowRes1, flowRes2, totalRes1, totalRes2
    logical :: hadANan
    ! Call to get the split norms
    call setRVec(g, flowRes1, turbRes1, totalRes1)

    ! Set some defaults:
    alpha		= 1.e-2_realType
    minlambda     = .01
    nfevals = 0
    flag = .True.
    lambda     = 1.0_realType
    ! Compute the two norms we need:
    call VecNorm(y, NORM_2, ynorm, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecNorm(f, NORM_2, fnorm, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call MatMult(dRdw, y, w, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    nfevals = nfevals + 1

    call VecDot(f, w, initslope, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (initslope > zero)  then
       initslope = -initslope
    end if

    if (initslope == 0.0_realType) then
       initslope = -1.0_realType
    end if
#ifdef USE_COMPLEX
    call VecWAXPY(w, cmplx(-lambda, 0.0), y, x, ierr)
    call EChk(ierr, __FILE__, __LINE__)
#else
    call VecWAXPY(w, -lambda, y, x, ierr)
    call EChk(ierr, __FILE__, __LINE__)
#endif

    ! Compute Function:
    call setW(w)
    call computeResidualNK()
    call setRVec(g, flowRes2, turbRes2, gnorm)

    nfevals = nfevals + 1

    ! Before we get to the actual line search we do two additional
    ! checks:

    ! 1. If the full step has a nan, we backtrack until we get a valid
    ! step. We then lower the NK switch tol such that the solver is
    ! forced back up to ANK or DADI/RK to keep going a bit further.
    !
    ! 2. If the turbulence residual goes up by large-ish factor (2.0),
    ! we pre-limit the step. The reason for this is that a unit step
    ! might lower the total residual, but the turb res could go up an
    ! order of magnitude or more.

    hadANan = .False.
    if (isnan(gnorm) .or. turbRes2 > 2.0*turbRes1) then
       ! Special testing for nans

       if (isnan(gnorm)) then
          hadANan = .True.
          call setUniformFlow()
          lambda = 0.5
       else
          ! Large turb jump
          lambda = lambda * (turbRes1 / turbRes2)
          lambda = max(lambda, 0.1)
       end if

       backtrack: do iter=1, 10
          ! Compute new x value:
#ifdef USE_COMPLEX
           call VecWAXPY(w, cmplx(-lambda, 0.0), y, x, ierr)
           call EChk(ierr, __FILE__, __LINE__)
#else
          call VecWAXPY(w, -lambda, y, x, ierr)
          call EChk(ierr, __FILE__, __LINE__)
#endif

          ! Compute Function
          call setW(w)
          call computeResidualNK()
          call setRVec(g, flowRes2, turbRes2, gnorm)

          nfevals = nfevals + 1

          if (isnan(gnorm)) then
             ! Just reset the flow, adjust the step back and keep
             ! going
             call setUniformFlow()
             lambda = lambda * .5
          else

             ! Sufficient reduction! Whoo! This is great we're done!
             if (0.5_realType*gnorm*gnorm <= 0.5_realType*fnorm*fnorm + alpha*initslope) then
                exit
             end if

             ! If we're less than min lambda, just take it. This could
             ! let the residual go up slightly. That's ok.
             if (lambda < minlambda)  then
                exit
             end if

             ! Otherwise, cut back the lambda
             lambda = lambda * 0.5

          end if
       end do backtrack

       if (hadANan) then
          ! Adjust the NK switch tolerance such that the ANK or DADI
          ! goes a little further.
          nk_switchtol = 0.8*(gnorm/totalR0)
       end if

       ! All finished with this "pre" line search.
       return
    end if

    ! Sufficient reduction from the basic step. This is the return for
    ! a unit step. This is what we want.
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

  subroutine LSNone(x, f, g, y, w, nfevals, flag, step)

    use constants
    use utils, only : EChk
    use communication
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
    real(kind=alwaysRealType) :: step
    real(kind=realType) :: tmp
    flag = .True.
    ! We just accept the step and compute the new residual at the new iterate
    nfevals = 0
    step = nk_fixedStep
    tmp = -step
    call VecWAXPY(w, tmp, y, x, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Compute new function:
    call setW(w)
    call computeResidualNK()
    call setRVec(g)

    nfevals = nfevals + 1
  end subroutine LSNone

  subroutine LSNM(x, f, g, y, w, fnorm, ynorm, gnorm, nfevals, flag, step)

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

    real(kind=alwaysRealType) :: fnorm, gnorm, ynorm
    real(kind=realType) :: alpha
    real(kind=alwaysRealType) :: step
    logical :: flag
    integer(kind=intType) :: nfevals
    !   Note that for line search purposes we work with with the related
    !   minimization problem:
    !      min  z(x):  R^n -> R,
    !   where z(x) = .5 * fnorm*fnorm, and fnorm = || f ||_2.
    !
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

    call VecDot(f, w, initslope, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (initslope > 0.0_realType)  then
       initslope = -initslope
    end if

    if (initslope == 0.0_realType) then
       initslope = -1.0_realType
    end if

    alpha = 1.0 ! Initial step length:
    backtrack: do iter=1, 10

       ! Compute new x value:
       call VecWAXPY(w, -alpha, y, x, ierr)
       call EChk(ierr, __FILE__, __LINE__)

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
    step = alpha
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
    use flowUtils, only : computeLamViscosity
    use BCRoutines, only : applyAllBC, applyAllBC_block
    use solverUtils, only : timeStep, computeUtau
    use residuals, only :residual, initRes, sourceTerms
    use oversetData, only : oversetPresent
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
          call computeLamViscosity(.False.)
          call computeEddyViscosity (.False.)
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

    ! Need to re-apply the BCs. The reason is that BC halos behind
    ! interpolated cells need to be recomputed with their new
    ! interpolated values from actual compute cells. Only needed for
    ! overset.
    if (oversetPresent) then
       do sps=1,nTimeIntervalsSpectral
          do nn=1,nDom
             call setPointers(nn, 1, sps)
             if (equations == RANSequations) then
                call BCTurbTreatment
                call applyAllTurbBCthisblock(.True.)
             end if
             call applyAllBC_block(.True.)
          end do
       end do
    end if

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
    call sourceTerms

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

  subroutine setRVec(rVec, flowRes, turbRes, totalRes)

    ! Set the current residual in dw into the PETSc Vector
    use constants
    use blockPointers, only : nDom, volRef, il, jl, kl, dw
    use inputtimespectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nw, nwf, nt1, nt2
    use inputIteration, only : turbResScale
    use utils, only : setPointers, EChk
    use communication, only : adflow_comm_world
    implicit none

    Vec    rVec
    integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
    real(kind=realType),pointer :: rvec_pointer(:)
    real(Kind=realType) :: ovv
    real(kind=alwaysRealType), intent(out), optional :: flowRes, turbRes, totalRes
    real(kind=realType) :: tmp, tmp2(2), flowResLocal, turbResLocal

    flowResLocal = zero
    turbResLocal = zero

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
                      tmp =  dw(i, j, k, l)*ovv
                      rvec_pointer(ii) = tmp
                      ii = ii + 1
                      flowResLocal = flowResLocal + tmp**2
                   end do
                   do l=nt1,nt2
                      tmp = dw(i, j, k, l)*ovv*turbResScale(l-nt1+1)
                      rvec_pointer(ii) = tmp
                      ii = ii + 1
                      turbResLocal = turbResLocal + tmp**2
                   end do
                end do
             end do
          end do
       end do
    end do

    call VecRestoreArrayF90(rVec,rvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    if (present(flowRes) .and. present(turbRes) .and. present(totalRes)) then
       call mpi_allreduce((/flowResLocal, turbResLocal/), tmp2, 2, adflow_real, &
         mpi_sum, ADflow_comm_world, ierr)
       flowRes = sqrt(tmp2(1))
       totalRes = sqrt(tmp2(1) + tmp2(2))
       if (tmp2(2) > zero) then
          turbRes = sqrt(tmp2(2))
       else
          turbRes = zero
       end if
    end if

  end subroutine setRVec

  subroutine setW(wVec)

    use constants
    use blockPointers, only : nDom, il, jl, kl, w
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowVarRefState, only : nwf, nt1, nt2, winf
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
                   do l=1,nwf
                      w(i,j,k,l) = wvec_pointer(ii)
                      ii = ii + 1
                   end do
                   ! Clip the turb to prevent negative turb SA
                   ! values. This is similar to the pressure
                   ! clip. Need to check this for other Turb models.
                   do l=nt1, nt2
                      w(i, j, k, l) = max(1e-6*winf(l), wvec_pointer(ii))
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

    ! Take in externallly generated states and set them in ADflow
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

    real(kind=alwaysrealType), intent(in) :: norm, old_norm, rtol_last
    real(kind=alwaysrealType), intent(out) :: rtol
    real(kind=alwaysrealType) :: rtol_max, gamma, alpha, alpha2, threshold, stol

    rtol_max  = 0.8_realType
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
#include "petsc/finclude/petscvec.h90"


  Mat  dRdw, dRdwPre
  Vec wVec, rVec, deltaW
  KSP  ANK_KSP

  ! Additional stuff required for KSP for turbulence only
  Mat dRdWTurb, dRdwPreTurb
  Vec wVecTurb, rVecTurb, deltaWTurb
  KSP ANK_KSP_Turb

  PetscFortranAddr   ctx(1)

  ! Options for ANK Solver
  logical :: useANKSolver
  integer(kind=intType) :: ANK_jacobianLag
  integer(kind=intType) :: ANK_subSpace
  integer(kind=intType) :: ANK_maxIter
  integer(kind=intType) :: ANK_asmOverlap
  integer(kind=intType) :: ANK_iluFill
  integer(kind=intType) :: ANK_innerPreConIts
  real(kind=realType)   :: ANK_rtol
  real(kind=realType)   :: ANK_switchTol
  real(kind=realType)   :: ANK_divTol = 10
  logical :: ANK_useTurbDADI
  logical :: ANK_coupled=.False.
  real(kind=realType) :: ANK_saRelax
  real(kind=realType) :: ANK_turbSwitchTol
  integer(kind=intType) :: ANK_nSubIterTurb
  real(kind=realType) :: ANK_turbcflscale

  ! Misc variables
  real(kind=realType) :: ANK_CFL, ANK_CFL0, ANK_CFLLimit, ANK_StepFactor, lambda
  real(kind=realType) :: ANK_stepInit, ANK_stepCutback, ANK_stepMin, ANK_stepExponent, ANK_CFLExponent
  real(kind=realType) :: ANK_secondOrdSwitchTol, ANK_coupledSwitchTol
  logical :: ANK_solverSetup=.False.
  logical :: ANK_turbSetup=.False.
  integer(kind=intTYpe) :: ANK_iter, ANK_Turb_iter
  integer(kind=intType) :: nState, nStateTurb
  real(kind=alwaysRealType) :: totalR_old ! for recording the previous residual
  real(kind=alwaysRealType) :: rtolLast ! for recording the previous relativel tolerance for Eisenstat-Walker

contains

  subroutine setupANKsolver

    ! Setup the PETSc objects for the Newton-Krylov
    ! solver. destroyNKsolver can be used to destroy the objects created
    ! in this function

    use constants
    use stencils, only : euler_PC_stencil, N_euler_PC
    use communication, only : adflow_comm_world
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowVarRefState, only : nw, viscous, nwf
    use ADjointVars , only: nCellsLocal
    use NKSolver, only : destroyNKSolver
    use utils, only : EChk
    use adjointUtils, only : myMatCreate, statePreAllocation
    implicit none

    ! Working Variables
    integer(kind=intType) :: ierr, nDimw, nDimWTurb
    integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
    integer(kind=intType) :: n_stencil
    integer(kind=intType), dimension(:, :), pointer :: stencil
    integer(kind=intType) :: level

    ! Make sure we don't have memory for the approximate and exact
    ! Newton solvers kicking around at the same time.
    call destroyNKSolver()

    if (.not. ANK_solverSetup) then

        if (ANK_coupled) then ! NK solver for flow variables, DADI for turbulence
            nState = nw ! coupled ank uses all variables
        else
            nState = nwf ! Frozen Turbulence
            nStateTurb = nw-nwf ! Number of states in turbulence model
        endif

       nDimW = nState * nCellsLocal(1_intTYpe) * nTimeIntervalsSpectral

       call VecCreate(ADFLOW_COMM_WORLD, wVec, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecSetSizes(wVec, nDimW, PETSC_DECIDE, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecSetBlockSize(wVec, nState, ierr)
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
       call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nState, stencil, n_stencil, &
            level, .False.)
       call myMatCreate(dRdwPre, nState, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
            __FILE__, __LINE__)

       call matSetOption(dRdwPre, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
       call EChk(ierr, __FILE__, __LINE__)
       deallocate(nnzDiagonal, nnzOffDiag)

       ! Set the mat_row_oriented option to false so that dense
       ! subblocks can be passed in in fortran column-oriented format
       call MatSetOption(dRdWPre, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Setup Matrix-Free dRdw matrix and its function
       call MatCreateMFFD(ADFLOW_COMM_WORLD, nDimW, nDimW, &
            PETSC_DETERMINE, PETSC_DETERMINE, dRdw, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call MatMFFDSetFunction(dRdw, FormFunction_mf, ctx, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call MatSetOption(dRdW, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       !  Create the linear solver context
       call KSPCreate(ADFLOW_COMM_WORLD, ANK_KSP, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set operators for the solver
       call KSPSetOperators(ANK_KSP, dRdw, dRdwPre, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Setup vectors and KSP for turbulence
       if ((.not. ANK_coupled) .and. (.not. ANK_useTurbDADI)) then
           nDimWTurb = nStateTurb * nCellsLocal(1_intTYpe) * nTimeIntervalsSpectral

           call VecCreate(ADFLOW_COMM_WORLD, wVecTurb, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call VecSetSizes(wVecTurb, nDimWTurb, PETSC_DECIDE, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call VecSetBlockSize(wVecTurb, nStateTurb, ierr)
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
           call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimWTurb/nStateTurb, stencil, n_stencil, &
                level, .False.)
           call myMatCreate(dRdwPreTurb, nStateTurb, nDimWTurb, nDimWTurb, nnzDiagonal, nnzOffDiag, &
                __FILE__, __LINE__)

           call matSetOption(dRdwPreTurb, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
           call EChk(ierr, __FILE__, __LINE__)
           deallocate(nnzDiagonal, nnzOffDiag)

           ! Set the mat_row_oriented option to false so that dense
           ! subblocks can be passed in in fortran column-oriented format
           call MatSetOption(dRdWPreTurb, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! Setup Matrix-Free dRdw matrix and its function
           call MatCreateMFFD(ADFLOW_COMM_WORLD, nDimWTurb, nDimWTurb, &
                PETSC_DETERMINE, PETSC_DETERMINE, dRdwTurb, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call MatMFFDSetFunction(dRdwTurb, FormFunction_mf_Turb, ctx, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call MatSetOption(dRdWTurb, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           !  Create the linear solver context
           call KSPCreate(ADFLOW_COMM_WORLD, ANK_KSP_Turb, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! Set operators for the solver
           call KSPSetOperators(ANK_KSP_Turb, dRdwTurb, dRdwPreTurb, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ANK_turbSetup = .True.
       end if

       ANK_solverSetup = .True.
       ANK_iter = 0
       ANK_Turb_iter = -1
    end if

  end subroutine setupANKsolver

  subroutine FormJacobianANK

    use constants
    use flowVarRefState, only : nw, nwf, nt1, nt2
    use blockPointers, only : nDom, volRef, il, jl, kl, dw
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputIteration, only : turbResScale
    use inputADjoint, only : viscPC
    use iteration, only : totalR0
    use utils, only : EChk, setPointers
    use adjointUtils, only :setupStateResidualMatrix, setupStandardKSP
    implicit none

    ! Local Variables
    character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
    integer(kind=intType) ::ierr
    logical :: useAD, usePC, useTranspose, useObjective, tmp
    real(kind=realType) ::  dt
    integer(kind=intType) :: i, j, k, l, ii, nn, sps, outerPreConIts, subspace
    real(kind=realType), pointer :: diag(:)

    ! Assemble the approximate PC (fine leve, level 1)
    useAD = .False.
    usePC = .True.
    useTranspose = .False.
    useObjective = .False.
    tmp = viscPC ! Save what is in viscPC and set to the NKvarible
    viscPC = .False.

    if (ANK_coupled) then
       ! The turbulence jacobian will only be accuate with AD.
       useAD = .True.
       ! get the full jacobian with turbulent variables included
       call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, &
            useObjective, .False., 1_intType)
    else
       call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, &
            useObjective, .True., 1_intType)
    end if
    ! Reset saved value
    viscPC = tmp

    ! ----------- Setup Flow KSP ----------

    !!! The routine of setting CFL can be done just by petsc functions, may consider changing here !!!
    ! the turbulent cfl can be scaled separately by VecStrideScale, did not help

    call VecGetArrayF90(deltaW, diag, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Calculate the contribution from the time-stepping term
    dt = one/ANK_CFL

    if (.not. ANK_coupled) then
        ! For the segragated solver, no need for scaling, each variable gets the same value
        diag(:) = dt
    else
        ! For the coupled solver, CFL number for the turbulent variable needs scaling
        ii = 1
        do nn=1, nDom
           do sps=1, nTimeIntervalsSpectral
              call setPointers(nn,1_intType,sps)
              do k=2, kl
                 do j=2, jl
                    do i=2, il
                        do l = 1, nwf
                          diag(ii) = dt
                          ii = ii + 1
                        end do
                        do l = nt1, nt2
                          diag(ii) = turbResScale(l-nt1+1)*dt
                          ii = ii + 1
                        end do
                    end do
                 end do
              end do
           end do
        end do
    end if

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
    if (ank_subspace < 0) then
       subspace = ANK_maxIter
    else
       subspace = ANK_subspace
    end if
    call setupStandardKSP(ANK_KSP, kspObjectType, subSpace, &
         preConSide, globalPCType, ANK_asmOverlap, outerPreConIts, localPCType, &
         localOrdering, ANK_iluFill, ANK_innerPreConIts)

    ! Don't do iterative refinement for the NKSolver.
    call KSPGMRESSetCGSRefinementType(ANK_KSP, &
         KSP_GMRES_CGS_REFINE_NEVER, ierr)
    call EChk(ierr, __FILE__, __LINE__)

  end subroutine FormJacobianANK

  subroutine FormJacobianANKTurb

    use constants
    use flowVarRefState, only : nw, nwf, nt1, nt2
    use blockPointers, only : nDom, volRef, il, jl, kl, dw
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputIteration, only : turbResScale
    use inputADjoint, only : viscPC
    use iteration, only : totalR0
    use utils, only : EChk, setPointers
    use adjointUtils, only :setupStateResidualMatrix, setupStandardKSP
    implicit none

    ! Local Variables
    character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
    integer(kind=intType) ::ierr
    logical :: useAD, usePC, useTranspose, useObjective, tmp, secondOrdSave
    real(kind=realType) ::  dt
    integer(kind=intType) :: i, j, k, l, ii, nn, sps, outerPreConIts, subspace
    real(kind=realType), pointer :: diag(:)

    ! Assemble the approximate PC (fine leve, level 1)
    useAD = .True. ! Turbulence Jacobian will only be accurate with AD
    usePC = .True.
    useTranspose = .False.
    useObjective = .False.
    tmp = viscPC ! Save what is in viscPC and set to the NKvarible
    viscPC = .False.

    ! Only form the preconditioner matrix for turbulence variable
    call setupStateResidualMatrix(dRdwPreTurb, useAD, usePC, useTranspose, &
         useObjective, .False., 1_intType, .True.)
    ! Reset saved values
    viscPC = tmp

    ! ----------- Setup Turbulence KSP ----------

    !!! The routine of setting CFL can be done just by petsc functions, may consider changing here !!!
    ! the turbulent cfl can be scaled separately by VecStrideScale, did not help

    call VecGetArrayF90(deltaWTurb, diag, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Calculate the contribution from the time-stepping term
    dt = turbResScale(1)/(ANK_CFL*ANK_turbCFLScale)
    diag(:) = dt

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
    if (ank_subspace < 0) then
       subspace = ANK_maxIter
    else
       subspace = ANK_subspace
    end if
    call setupStandardKSP(ANK_KSP_Turb, kspObjectType, subSpace, &
         preConSide, globalPCType, ANK_asmOverlap, outerPreConIts, localPCType, &
         localOrdering, ANK_iluFill+2, ANK_innerPreConIts)

    ! Don't do iterative refinement for the NKSolver.
    call KSPGMRESSetCGSRefinementType(ANK_KSP_Turb, &
         KSP_GMRES_CGS_REFINE_NEVER, ierr)
    call EChk(ierr, __FILE__, __LINE__)

  end subroutine FormJacobianANKTurb

  subroutine FormFunction_mf(ctx, wVec, rVec, ierr)

    ! This is the function used for the matrix-free matrix-vector products
    ! for the GMRES solver used in ANK

    use constants
    use blockPointers, only : nDom, volRef, il, jl, kl, dw
    use inputtimespectral, only : nTimeIntervalsSpectral
    use inputIteration, only : turbResScale
    use flowvarrefstate, only : nwf, nt1, nt2
    use NKSolver, only : computeResidualNK, setRvec
    use utils, only : setPointers, EChk
    implicit none

    ! PETSc Variables
    PetscFortranAddr ctx(*)
    Vec     wVec, rVec
    real(kind=realType) :: dt
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
    real(kind=realType),pointer :: rvec_pointer(:)
    real(kind=realType),pointer :: wvec_pointer(:)

    ! get the input vector
    call setWANK(wVec)

    if (ANK_coupled) then ! calculate all variables
      call computeResidualNK()
      call setRVec(rVec)
    else ! Flow variables only
      call computeResidualANK()
      call setRVecANK(rVec)
    end if

    ! Calculate the contribution from the time stepping term
    dt = one/ANK_CFL

    ! Add the contribution from the diagonal time stepping term
    if (.not. ANK_coupled) then
      ! For the segragated solver each variable gets the same dt value
      call VecAXPY(rVec, dt, wVec, ierr)
      call EChk(ierr,__FILE__,__LINE__)
    else
      ! For the coupled solver, time stepping term for turbulence needs to be scaled
      call VecGetArrayF90(rVec,rvec_pointer,ierr)
      call EChk(ierr,__FILE__,__LINE__)

      call VecGetArrayReadF90(wVec,wvec_pointer,ierr)
      call EChk(ierr,__FILE__,__LINE__)

      ii = 1
      do nn=1, nDom
         do sps=1, nTimeIntervalsSpectral
            call setPointers(nn,1_intType,sps)
            ! read the density residuals and set local CFL
            do k=2, kl
               do j=2, jl
                  do i=2, il
                      do l = 1, nwf
                        rvec_pointer(ii) = rvec_pointer(ii) + wvec_pointer(ii)*dt
                        ii = ii + 1
                      end do
                      do l = nt1, nt2
                        rvec_pointer(ii) = rvec_pointer(ii) + wvec_pointer(ii)*turbResScale(l-nt1+1)*dt
                        ii = ii + 1
                      end do
                  end do
               end do
            end do
         end do
      end do

      call VecRestoreArrayF90(rVec, rvec_pointer, ierr)
      call EChk(ierr,__FILE__,__LINE__)

      call VecRestoreArrayReadF90(wVec, wvec_pointer, ierr)
      call EChk(ierr,__FILE__,__LINE__)
    end if
    ! We don't check an error here, so just pass back zero
    ierr = 0

  end subroutine FormFunction_mf

  subroutine FormFunction_mf_Turb(ctx, wVecTurb, rVecTurb, ierr)
      ! This is the function used for the matrix-free matrix-vector products
      ! for the GMRES solver used in ANK

      use constants
      use blockPointers, only : nDom, volRef, il, jl, kl, dw
      use inputtimespectral, only : nTimeIntervalsSpectral
      use inputIteration, only : turbResScale
      use flowvarrefstate, only : nwf, nt1, nt2
      use NKSolver, only : computeResidualNK, setRvec
      use utils, only : setPointers, EChk
      implicit none

      ! PETSc Variables
      PetscFortranAddr ctx(*)
      Vec     wVecTurb, rVecTurb
      real(kind=realType) :: dt
      integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
      real(kind=realType),pointer :: rvecTurb_pointer(:)
      real(kind=realType),pointer :: wvecTurb_pointer(:)

      ! get the input vector
      call setWANKTurb(wVecTurb)

      ! Compute turbulent residual only
      call computeResidualANKTurb()

      ! set the R vec for turbulence in petsc
      call setRVecANKTurb(rVecTurb)

      ! Calculate the contribution from the time stepping term
      dt = turbResScale(1)/(ANK_CFL*ANK_turbcflscale)

      ! For the segragated solver each variable gets the same dt value
      call VecAXPY(rVecTurb, dt, wVecTurb, ierr)
      call EChk(ierr,__FILE__,__LINE__)

      ! We don't check an error here, so just pass back zero
      ierr = 0

  end subroutine FormFunction_mf_Turb

  subroutine destroyANKsolver

    ! Destroy all the PETSc objects for the Newton-Krylov
    ! solver.

    use constants
    use utils, only : EChk
    implicit none
    integer(kind=intType) :: ierr

    if (ANK_SolverSetup) then

       call MatDestroy(dRdw, ierr)
       call EChk(ierr, __FILE__, __LINE__)

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

           call MatDestroy(dRdwTurb, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call MatDestroy(dRdwPreTurb, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call VecDestroy(wVecTurb, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call VecDestroy(rVecTurb, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call VecDestroy(deltaWTurb, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call KSPDestroy(ANK_KSP_Turb, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ANK_turbSetup = .False.
       end if

       ANK_SolverSetup = .False.
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
    use flowUtils, only : computeLamViscosity
    use BCRoutines, only : applyAllBC
    use solverUtils, only : timeStep, computeUtau
    use residuals, only :residual, initRes, sourceTerms
    use oversetData, only : oversetPresent

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
          call computeLamViscosity(.False.)
          call computeEddyViscosity(.False.)
       end do domainsState
    end do spectralLoop

    ! Apply BCs
    call applyAllBC(secondHalo)

    ! Exchange halos
    call whalo2(currentLevel, 1_intType, nwf, .true., &
         .true., .true.)

    ! Need to re-apply the BCs. The reason is that BC halos behind
    ! interpolated cells need to be recomputed with their new
    ! interpolated values from actual compute cells. Only needed for
    ! overset.
    if (oversetPresent) then
       call applyAllBC(secondHalo)
    end if

    ! Compute time step (spectral radius is actually what we need)
    call timestep(.false.)

    ! Initialize Flow residuals
    call initres(1_intType, nwf)
    call sourceTerms

    ! Actual Residual Calc
    call residual

  end subroutine computeResidualANK

  subroutine computeResidualANKTurb

    ! This is the residual evaluation driver for the KSP for Turbulence in ANK

    ! This routine might be cleaned up a bit, depending on what values we are using during turbulence residual calc

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
    use flowUtils, only : computeLamViscosity
    use BCRoutines, only : applyAllBC, applyAllBC_block
    use solverUtils, only : timeStep, computeUtau
    use residuals, only :residual, initRes, sourceTerms
    use oversetData, only : oversetPresent
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
          ! Compute Viscosities
          call computeLamViscosity(.False.)
          call computeEddyViscosity (.False.)
       end do domainsState
    end do spectralLoop

    ! Apply BCs
    do nn=1,nDom
      do sps=1,nTimeIntervalsSpectral
         call setPointers(nn, currentLevel, sps)
         call bcTurbTreatment
         call applyAllTurbBCThisBLock(.True.)
      end do
    end do

    ! Exchange halos
    call whalo2(currentLevel, nt1, nt2, .false., &
         .false., .True.)

    ! Need to re-apply the BCs. The reason is that BC halos behind
    ! interpolated cells need to be recomputed with their new
    ! interpolated values from actual compute cells. Only needed for
    ! overset.
    if (oversetPresent) then
       do sps=1,nTimeIntervalsSpectral
          do nn=1,nDom
             call setPointers(nn, 1, sps)
             call BCTurbTreatment
             call applyAllTurbBCthisblock(.True.)
          end do
       end do
    end if

    ! Compute time step (spectral radius is actually what we need)
    call timestep(.false.)

    ! Compute the skin-friction velocity (wall functions only)
    call computeUtau
    call initres(nt1, nt2) ! Initialize only the Turblent Variables
    call turbResidual

  end subroutine computeResidualANKTurb

  subroutine setWVecANK(wVec)
    ! Set the current FLOW variables in the PETSc Vector

    use constants
    use blockPointers, only : nDom, il, jl, kl, w
    use inputtimespectral, only : ntimeIntervalsSpectral
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
                   do l=1, nState
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

  subroutine setWVecANKTurb(wVecTurb)
    ! Set the current Turbulence variables in the PETSc Vector

    use constants
    use blockPointers, only : nDom, il, jl, kl, w
    use inputtimespectral, only : ntimeIntervalsSpectral
    use flowvarrefstate, only : nt1, nt2
    use utils, only : setPointers, EChk
    implicit none

    Vec   wVecTurb
    integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
    real(kind=realType),pointer :: wvecTurb_pointer(:)

    call VecGetArrayF90(wVecTurb,wvecTurb_pointer,ierr)
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
                      wvecTurb_pointer(ii) = w(i, j, k, l)
                   end do
                end do
             end do
          end do
       end do
    end do

    call VecRestoreArrayF90(wVecTurb, wvecTurb_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine setWVecANKTurb

  subroutine setRVecANK(rVec)

    ! Set the current FLOW residual in dw into the PETSc Vector
    use constants
    use blockPointers, only : nDom, volRef, il, jl, kl, dw
    use inputtimespectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nwf, nt1, nt2
    use inputIteration, only : turbResScale
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

  subroutine setRVecANKTurb(rVecTurb)

    ! Set the current Turb residual in dw into the PETSc Vector
    use constants
    use blockPointers, only : nDom, volRef, il, jl, kl, dw
    use inputtimespectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nt1, nt2
    use inputIteration, only : turbResScale
    use utils, only : setPointers, EChk
    implicit none
    Vec    rVecTurb
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
    real(kind=realType),pointer :: rvecTurb_pointer(:)
    real(Kind=realType) :: ovv
    call VecGetArrayF90(rVecTurb,rvecTurb_pointer,ierr)
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
                   do l=nt1, nt2
                      ii = ii + 1
                      rvecTurb_pointer(ii) = dw(i, j, k, l)*ovv*turbResScale(l-nt1+1)
                   end do
                end do
             end do
          end do
       end do
    end do

     call VecRestoreArrayF90(rVecTurb, rvecTurb_pointer, ierr)
     call EChk(ierr,__FILE__,__LINE__)

 end subroutine setRVecANKTurb

  subroutine setWANK(wVec)
    ! Get the updated solution from the PETSc Vector

    use constants
    use blockPointers, only : nDom, vol, il, jl, kl, w
    use inputtimespectral, only : nTimeIntervalsSpectral
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
                   do l=1, nState
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

  subroutine setWANKTurb(wVecTurb)
    ! Get the updated solution from the PETSc Vector

    use constants
    use blockPointers, only : nDom, vol, il, jl, kl, w
    use inputtimespectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nt1, nt2
    use utils, only : setPointers, EChk
    implicit none

    Vec  wVecTurb
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
    real(kind=realType), pointer :: wvecTurb_pointer(:)
    call VecGetArrayReadF90(wVecTurb, wvecTurb_pointer, ierr)
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
                      w(i, j, k, l) = wvecTurb_pointer(ii)
                   end do
                end do
             end do
          end do
       end do
    end do
    call VecRestoreArrayReadF90(wVecTurb, wvecTurb_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine setWANKTurb

  ! for debugging/experimenting
  subroutine scaleDeltaWTurb(deltaWTurb, wVecTurb)
    ! Get the updated solution from the PETSc Vector

    use constants
    use blockPointers, only : nDom, vol, il, jl, kl, w, volRef
    use inputtimespectral, only : nTimeIntervalsSpectral
    use flowvarrefstate, only : nt1, nt2
    use utils, only : setPointers, EChk
    implicit none

    Vec  deltaWTurb, wVecTurb
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii, counter
    real(kind=realType), pointer :: dvecTurb_pointer(:), wVecTurb_pointer(:)
    call VecGetArrayReadF90(deltaWTurb, dvecTurb_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecGetArrayReadF90(wVecTurb, wVecTurb_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ii = 0
    counter = 0
    do nn=1, nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn, 1_intType, sps)
          do k=2, kl
             do j=2, jl
                do i=2, il
                   do l=nt1, nt2
                      ii = ii + 1
                      if (abs(max(dvecturb_pointer(ii), 0.0_realType)) > wVecTurb_pointer(ii)) then
                        dvecturb_pointer(ii) = 0.9_realType*wVecTurb_pointer(ii)
                        counter = counter + 1
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do
    call VecRestoreArrayReadF90(deltaWTurb, dvecTurb_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecRestoreArrayReadF90(wVecTurb, wVecTurb_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    if (counter > 0)&
        write(*,*) "Number of turb updates clipped: ",counter

  end subroutine scaleDeltaWTurb

  subroutine ANKStep(firstCall)

    use constants
    use blockPointers, only : nDom, flowDoms, shockSensor, ib, jb, kb, p, w, gamma
    use communication, only : myid
    use inputPhysics, only : equations
    use inputIteration, only : L2conv, nsubiterturb, turbResScale
    use inputDiscretization, only : lumpedDiss, sa_relax
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use iteration, only : approxTotalIts, totalR0, totalR, stepMonitor, linResMonitor, currentlevel
    use utils, only : EChk, setpointers
    use turbAPI, only : turbSolveSegregated
    use turbMod, only : secondOrd
    use solverUtils, only : computeUTau
    use adjointUtils, only : referenceShockSensor
    use NKSolver, only : setRVec, computeResidualNK, getEwTol
    use initializeFlow, only : setUniformFlow
    use haloexchange, only : whalo2
    use BCRoutines, only : applyAllBC, applyAllBC_block
    use flowvarrefstate, only : nwf, nw, kPresent, pInfCorr
    use oversetData, only : oversetPresent
    use flowUtils, only : computeLamViscosity
    use turbUtils, only : computeEddyViscosity
    use communication

    implicit none

    ! Input Variables
    logical, intent(in) :: firstCall

    ! Working Variables
    integer(kind=intType) :: ierr, maxIt, kspIterations, kspIterations_Turb, nn, sps, reason, nHist, iter, n_turb
    integer(kind=intType) :: i,j,k
    real(kind=realType) :: atol, val, lambdaBT, v2, factK, gm1
    real(kind=alwaysRealType) :: rtol, totalR_dummy, linearRes, norm
    real(kind=alwaysRealType) :: resHist(ank_maxIter+1)
    real(kind=alwaysRealType) :: resHist_turb(ank_maxIter+1)
    logical :: secondOrdSave, correctForK

    ! Enter this check if this is the first ANK step OR we are switching to the coupled ANK solver
    if (firstCall .or. (totalR < ANK_coupledSwitchTol * totalR0 .and. (.not. ANK_coupled) ) &
        .or. (totalR < ANK_turbSwitchTol * totalR0 .and. ANK_useTurbDADI) ) then

       ! If using segragated ANK and below the coupled switch tol, set ANK_useTurbDADI
       ! to .False. to create the PETSc objets required for the coupled ANK solver
       if (totalR < ANK_coupledSwitchTol * totalR0) then
         ANK_coupled = .True.
         call destroyANKSolver()
       else if (totalR < ANK_turbSwitchTol * totalR0) then
         ANK_useTurbDADI = .False.
         ANK_coupled = .False.
         call destroyanksolver()
       else
         ANK_coupled = .False.
       end if

       call setupANKSolver()

       ! Copy the adflow 'w' into the petsc wVec
       call setwVecANK(wVec)

       ! Evaluate the residual before we start
       if (ANK_coupled) then ! Compute and set the full residual
          call computeResidualNK()
          call setRVec(rVec)
       else
          call computeResidualANK() ! Compute only flow residuals
          call setRVecANK(rVec) ! Set flow residuals

          ! If not DADI, Turb gets updates via KSP
          if (.not. ANK_useTurbDADI) then ! Segragated ANK with KSP for turbulence
              call computeResidualANKTurb() ! Compute turbulent residuals only
              call setRVecANKTurb(rVecTurb) ! Set turbulence residual vector in PETSc
              call setWVecANKTurb(wVecTurb) ! Set turbulence state vector in PETSc
          end if
       end if

       totalR_old = totalR ! Record the old residual for the first iteration
       rtolLast = ANK_rtol ! Set the previous relative convergence tolerance for the first iteration

       if (firstCall) then
         ! Start with the selected fraction of the ANK_StepFactor
         lambda = ANK_stepInit*ANK_StepFactor
       end if
    else
       ANK_iter = ANK_iter + 1
    end if

    ! ANK CFL calculation, use relative convergence w.r.t. totalR0 for better performance with restarts
    ANK_CFL = min(ANK_CFL0 * (totalR0 / totalR)**ANK_CFLExponent, ANK_CFLLimit)

    ! Determine if if we need to form the Preconditioner
    if (mod(ANK_iter, ANK_jacobianLag) == 0) then
       call FormJacobianANK()
    end if

    ! Dummy matrix assembly for the matrix-free matrix
    call MatAssemblyBegin(dRdw, MAT_FINAL_ASSEMBLY, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call MatAssemblyEnd(dRdw, MAT_FINAL_ASSEMBLY, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set the BaseVector of the matrix-free matrix:
    call MatMFFDSetBase(dRdw, wVec, PETSC_NULL_OBJECT, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Continuation for step factor:
    ! If total residual have increased in the previous iteration,
    ! reduce the step wrt the cutback factor, until ANK_StepMin factor of the step factor is reached
    if (totalR > totalR_old) then
      lambda = max(lambda*ANK_stepCutback, ANK_StepFactor*ANK_stepMin)
    ! If total residual have decreased, slowly ramp the step up
    else
      lambda = min(lambda*(totalR_old/totalR)**ANK_stepExponent, ANK_StepFactor)
    end if

    ! ============== Flow Update =============

    ! For the approximate solver, we need the approximate flux routines
    ! We set the variables required for approximate fluxes here and they will be used
    ! for the matrix-free matrix-vector product routines when the KSP solver calls it
    ! Very important to set the variables back to their original values after each
    ! KSP solve because we want actual flux functions when calculating residuals
    if (totalR > ANK_secondOrdSwitchTol*totalR0) then
      ! Setting lumped dissipation to true gives approximate fluxes
      lumpedDiss =.True.

      ! Save if second order turbulence is used, we will only use 1st order during ANK (only matters for the coupled solver)
      secondOrdSave = secondOrd
      secondOrd =.False.

      ! Calculate the shock sensor here because the approximate routines do not
      call referenceShockSensor()

      ! Determine the relative convergence for the KSP solver
      rtol = ANK_rtol ! Just use the input relative tolerance for approximate fluxes

   else
      ! If the second order fluxes are used, Eisenstat-Walker algorithm to determine relateive
      ! convergence tolerance helps with performance.
      totalR_dummy = totalR
      call getEWTol(totalR_dummy, totalR_old, rtolLast, rtol)

   end if

    ! Record the total residual and relative convergence for next iteration
    totalR_old = totalR
    rtolLast = rtol

    ! Set all tolerances for linear solve:
    atol = totalR0*L2Conv

    ! Set the iteration limit to maxIt, determined by which fluxes are used.
    ! This is because ANK step require 0.1 convergence for stability during initial stages.
    ! Due to an outdated preconditioner, the KSP solve might take more iterations.
    ! If this happens, the preconditioner is re-computed and because of this,
    ! ANK iterations usually don't take more than 2 times number of ANK_subSpace size iterations
    call KSPSetTolerances(ANK_KSP, rtol, &
         real(atol), real(ANK_divTol), ank_maxIter, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call KSPSetResidualHistory(ANK_KSP, resHist, ank_maxIter+1, PETSC_TRUE, ierr)
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

    ! Return previously changed variables back to normal, VERY IMPORTANT
    if (totalR > ANK_secondOrdSwitchTol*totalR0) then
      ! Set lumpedDiss back to False to go back to using actual flux routines
      lumpedDiss =.False.

      ! Replace the second order turbulence option
      secondOrd = secondOrdSave

      ! Deallocate the memory used for the shock sensor
      do nn=1, nDom
          do sps=1, nTimeIntervalsSpectral
              deallocate(flowDoms(nn,1,sps)%shockSensor)
          end do
      end do
    end if

    ! No line search...just take the new solution, possibly (fixed)
    ! limited
    call VecAXPY(wVec, -lambda, deltaW, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    stepMonitor = lambda
    ! Set the updated state variables
    call setWANK(wVec)

    ! ============== Turb Update =============
    if ((.not. ANK_coupled) .and. equations == RANSEquations) then ! Turb either gets dadi or KSP
        ! Update the intermediate variables required for turb. solve
        if (kPresent) then
           correctForK = .True.
        end if
        ! because the flow variables has changed.
        do sps=1,nTimeIntervalsSpectral
           do nn=1,nDom
             call setPointers(nn, 1, sps)
             ! Calculate pressure
             factK = zero
             do k=0, kb
                do j=0, jb
                   do i=0, ib

                      gm1  = gamma(i, j, k) - one
                      v2 = w(i,j,k,ivx)**2 + w(i,j,k,ivy)**2 &
                           + w(i,j,k,ivz)**2

                      p(i,j,k) = gm1*(w(i,j,k,irhoE) &
                           - half*w(i,j,k,irho)*v2)

                      if( correctForK ) then
                         factK = five*third - gamma(i, j ,k)
                         p(i, j ,K) = p(i,j, k) + factK*w(i, j, k, irho) &
                              * w(i, j, k, itu1)
                      end if

                      ! Clip to make sure it is positive.
                      p(i,j,k) = max(p(i,j,k), 1.e-4_realType*pInfCorr)
                   end do
                end do
             end do
             ! Compute Viscosities
             call computeLamViscosity(.False.)
             call computeEddyViscosity (.False.)
           end do
        end do

        ! Apply BCs
        call applyAllBC(.true.)

        ! Exchange halos
        call whalo2(currentLevel, 1_intType, nwf, .True., &
             .false., .true.)

        ! Need to re-apply the BCs. The reason is that BC halos behind
        ! interpolated cells need to be recomputed with their new
        ! interpolated values from actual compute cells. Only needed for
        ! overset.
        if (oversetPresent) then
           do sps=1,nTimeIntervalsSpectral
              do nn=1,nDom
                 call setPointers(nn, 1, sps)
                 call applyAllBC_block(.True.)
              end do
           end do
        end if

        if (ANK_useTurbDADI) then ! Do DDADI update
            ! parameter to control the approximations in the ddadi jacobian for turbulence
            !sa_relax = min(ANK_saRelax*totalR0/totalR, one)
            call computeUtau
            call turbSolveSegregated
        else ! Do ksp update
          do n_turb = 1, ANK_nsubiterturb ! Repeat for desired turbulent sub-iterations

            ! Re-calculate the turbulent residuals
            call computeResidualANKTurb()

            ! Set the residual vector in Petsc
            call setRVecANKTurb(rVecTurb)

            ! Increment the iteration counter
            ANK_Turb_iter = ANK_Turb_iter + 1

            ! Check if PC needs a refresh
            if (mod(ANK_Turb_iter, ANK_jacobianLag) == 0) then
               call FormJacobianANKTurb()
            end if

            ! Dummy matrix assembly for the mat-free matrix
            call MatAssemblyBegin(dRdwTurb, MAT_FINAL_ASSEMBLY, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            call MatAssemblyEnd(dRdwTurb, MAT_FINAL_ASSEMBLY, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Set base vector of the MFFD matrix
            call MatMFFDSetBase(dRdwTurb, wVecTurb, PETSC_NULL_OBJECT, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Set different stuff required for MFFD operations
            if (totalR > ANK_secondOrdSwitchTol*totalR0) then
              ! Setting lumped dissipation to true gives approximate fluxes
              lumpedDiss =.True. ! Doesn't do anything for turbulence, might remove

              ! Save if second order turbulence is used, we will only use 1st order until 2nd ord switchtol is reached
              secondOrdSave = secondOrd
              secondOrd =.False.

              ! Determine the relative convergence for the KSP solver
              rtol = ANK_rtol ! Use the same relative tolerance for the turbulence
            else
              ! If the second order fluxes are used, Eisenstat-Walker algorithm to determine relateive
              ! convergence tolerance helps with performance.
              !totalR_dummy = totalR
              !call getEWTol(totalR_dummy, totalR_old, rtolLast, rtol)

              rtol = ANK_rtol ! Use the same value for now
            end if

            ! Set R tol for the KSP. Turb residuals should already be scaled, so they should be same order wrt flow residuals
            atol = totalR0*L2Conv

            ! Set KSP tolerances, get residual history

            ! Set the iteration limit to maxIt, determined by which fluxes are used.
            ! This is because ANK step require 0.1 convergence for stability during initial stages.
            ! Due to an outdated preconditioner, the KSP solve might take more iterations.
            ! If this happens, the preconditioner is re-computed and because of this,
            ! ANK iterations usually don't take more than 2 times number of ANK_subSpace size iterations
            call KSPSetTolerances(ANK_KSP_Turb, rtol, &
                 real(atol), real(ANK_divTol), ank_maxIter, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            call KSPSetResidualHistory(ANK_KSP_turb, resHist_turb, ank_maxIter+1, PETSC_TRUE, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Actually do the KSP Solve
            call KSPSolve(ANK_KSP_turb, rVecTurb, deltaWTurb, ierr)

            ! Check for NaN
            if (ierr == 72) then
               ! The convergence check will get the nan
            else
               call EChk(ierr, __FILE__, __LINE__)
            end if

            ! Return previously changed variables back to normal, VERY IMPORTANT
            if (totalR > ANK_secondOrdSwitchTol*totalR0) then
              ! Set lumpedDiss back to False to go back to using actual flux routines
              lumpedDiss =.False. ! again, shouldn't do anything for turbulence

              ! Replace the second order turbulence option
              secondOrd = secondOrdSave
            end if

            ! Scale the update to prevent negative turbulence variables
            call scaleDeltaWTurb(deltaWTurb, wVecTurb)

            ! take the update
            call VecAXPY(wVecTurb, -lambda, deltaWTurb, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Set the updated state variables
            call setWANKTurb(wVecTurb)

            ! Get the number of iterations from the KSP solver
            call KSPGetIterationNumber(ANK_KSP_Turb, kspIterations_Turb, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Get reason for convergence, for debugging purposes only
            call KSPGetConvergedReason(ANK_KSP_Turb, reason, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Print some info about the turbulence ksp
            if (myid == 0) then
              Write(*,*) "LIN RES, ITER, INITRES, REASON", resHist_turb(kspIterations_turb+1)/resHist_turb(1),kspIterations_turb, &
              reshist_turb(1), reason
            end if
            if (kspIterations_turb > .5 * ank_maxIter .and. totalR > ANK_secondOrdSwitchTol*totalR0) then
               ! We should reform the PC since it took longer than we want.
               ANK_Turb_iter = -1
            end if
          end do
        end if
    end if

    ! Calculate the residual with the new values and set the R vec in PETSc
    ! We calculate the full residual using NK routine because the dw values
    ! for turbulent variable has the update, not the residual and this
    ! gives the wrong turbulent residual norm. This causes issues when
    ! switching to NK or restarts after ANK. To fix, we calculate turbulent
    ! residuals after each ANK step, which the turbSolveSegregated routine
    ! does not do on its own
    call computeResidualNK()
    if (ANK_coupled) then
       call setRVec(rVec)
    else
       call setRVecANK(rVec)
    end if

    ! Check if the norm of the rVec is bad:
    call VecNorm(rVec, NORM_2, norm, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (isnan(norm)) then

       ! Do backtracking linesearch:
       call setUniformFlow()

       ! Restore the starting (old) w value
       call VecAXPY(wVec, lambda, deltaW, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set the initial new lambda
       lambdaBT = 0.5 * lambda

       backtrack: do iter=1, 10

          ! Apply the new step
          call VecAXPY(wVec, -lambdaBT, deltaW, ierr)
          call EChk(ierr, __FILE__, __LINE__)

          ! Set and recompute
          call setWANK(wVec)
          call computeResidualNK()
          if (ANK_coupled) then
             call setRVec(rVec)
          else
             call setRVecANK(rVec)
          end if
          call VecNorm(rVec, NORM_2, norm, ierr)
          call EChk(ierr, __FILE__, __LINE__)

         if (isnan(norm)) then

            ! Restore back to the original wVec
            call VecAXPY(wVec, lambdaBT, deltaW, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Haven't backed off enough yet....keep going
            call setUniformFlow()
            lambdaBT = lambdaBT * .5
         else
            ! We don't have an nan anymore...break out
            exit
         end if
      end do backtrack
      stepMonitor = lambdaBT
   end if

    ! Get the number of iterations from the KSP solver
    call KSPGetIterationNumber(ANK_KSP, kspIterations, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call KSPGetConvergedReason(ANK_KSP, reason, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    linResMonitor = resHist(kspIterations+1)/resHist(1)

    if (kspIterations > .5 * ank_maxIter .and. totalR > ANK_secondOrdSwitchTol*totalR0) then
       ! We should reform the PC since it took longer than we want.
       ANK_iter = -1
    end if

    ! Update the approximate iteration counter. The +1 is for the
    ! residual evaluations.
    approxTotalIts = approxTotalIts + 1 + kspIterations

  end subroutine ANKStep
end module ANKSolver
