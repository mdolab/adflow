module NKSolver

  use constants

  ! MPI comes from constants, so we need to avoid MPIF_H in PETSc
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
#else
  implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif

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

  Vec wVec, rVec, deltaW, work, g, baseRes

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
    use inputIteration, only : useLinResMonitor
    use flowVarRefState, only : nw, viscous
    use InputAdjoint, only: viscPC, precondtype
    use ADjointVars , only: nCellsLocal
    use utils, only : EChk
    use adjointUtils, only : myMatCreate, statePreAllocation
    use agmg, only : setupAGMG
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

       call VecDuplicate(wVec, baseRes, ierr)
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

       if (preCondType == 'mg') then
          call setupAGMG(drdwpre, nDimW/nw, nw)
       end if

       !  Create the linear solver context
       call KSPCreate(ADFLOW_COMM_WORLD, NK_KSP, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set operators for the solver
       call KSPSetOperators(NK_KSP, dRdw, dRdwPre, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       if (useLinResMonitor) then
#if PETSC_VERSION_GE(3,8,0)
          ! This could be wrong. There is no petsc_null_context???
          call KSPMonitorSet(NK_KSP, linearResidualMonitor, PETSC_NULL_FUNCTION, &
               PETSC_NULL_FUNCTION, ierr)
#else
          call KSPMonitorSet(NK_KSP, linearResidualMonitor, PETSC_NULL_OBJECT, &
               PETSC_NULL_FUNCTION, ierr)
#endif
          call EChk(ierr, __FILE__, __LINE__)
       end if

       NK_solverSetup = .True.
       NK_iter = 0
    end if

  end subroutine setupNKsolver

  subroutine linearResidualMonitor(myKSP, n, rnorm, dummy, ierr)
    use communication, only : myid
    implicit none
    !
    !     Subroutine arguments.
    !
    ! myKsp - Iterative context
    ! n     - Iteration number
    ! rnorm - 2-norm (preconditioned) residual value
    ! dummy - Optional user-defined monitor context (unused here)
    ! ierr  - Return error code

    KSP myKSP
    integer(kind=intType) :: n, dummy, ierr
    real(kind=alwaysRealType)   :: rnorm

    ! Write the residual norm to stdout every adjMonStep iterations.
    if (myid == 0) then
       print *, n, rnorm
    end if
    ierr = 0
  end subroutine LinearResidualMonitor

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
    use inputADjoint, only : viscPC, precondType
    use utils, only : EChk
    use adjointUtils, only :setupStateResidualMatrix, setupStandardKSP, setupStandardMultigrid
    implicit none

    ! Local Variables
    character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
    integer(kind=intType) :: ierr
    logical :: useAD, usePC, useTranspose, useObjective, tmp
    integer(kind=intType) :: i, j, k, l, ii, nn, sps
    logical :: useCoarseMats

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

    if (preCondType == 'mg') then
       useCoarseMats = .True.
    else
       useCoarseMats = .False.
    end if

    call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, &
         useObjective, .False., 1_intType, useCoarseMats=useCoarseMats)
    ! Reset saved value
    viscPC = tmp

    ! Setup KSP Options
    preConSide = 'right'
    localPCType = 'ilu'
    kspObjectType = 'gmres'
    globalPCType = 'asm'
    localOrdering = 'rcm'

    ! Setup the KSP using the same code as used for the adjoint
    if (PreCondType == 'asm') then
       call setupStandardKSP(NK_KSP, kspObjectType, NK_subSpace, &
            preConSide, globalPCType, NK_asmOverlap, NK_outerPreConIts, localPCType, &
            localOrdering, NK_iluFill, NK_innerPreConIts)
    else
       call setupStandardMultigrid(NK_KSP, kspObjectType, NK_subSpace, &
            preConSide, NK_asmOverlap, NK_outerPreConIts, &
            localOrdering, NK_iluFill)
    end if


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
    use agmg, only : destroyAGMG
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

       call VecDestroy(baseRes, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecDestroy(g, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call VecDestroy(work, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call KSPDestroy(NK_KSP, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call destroyAGMG()

       NK_solverSetup = .False.
    end if

  end subroutine destroyNKsolver

  subroutine NKStep(firstCall)

    use constants
    use flowVarRefState, only : nw
    use inputPhysics, only : equations
    use flowVarRefState, only :  nw, nwf
    use inputIteration, only : L2conv
    use iteration, only : approxTotalIts, totalR0, stepMonitor, LinResMonitor, iterType
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
       iterType = "     *NK"
       call FormJacobianNK()
    else
       call MatAssemblyBegin(dRdw, MAT_FINAL_ASSEMBLY, ierr)
       call EChk(ierr, __FILE__, __LINE__)
       call MatAssemblyEnd(dRdw, MAT_FINAL_ASSEMBLY, ierr)
       call EChk(ierr, __FILE__, __LINE__)
       iterType = "      NK"
    end if

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

    ! set the BaseVector of the matrix-free matrix
    call formFunction_mf(ctx, wVec, baseRes, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call MatMFFDSetBase(dRdW, wVec, baseRes, ierr)
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
    use utils, only : EChk, myisnan
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
    alpha = 1.e-2_realType
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

    use constants
    use blockette, only : blocketteRes
    implicit none

    logical :: updateDt

    ! We want to update the time step
    updateDt = .true.

    ! Shell function to maintain backward compatibility with code using computeResidualNK
    call blocketteRes(useUpdateDT = updateDt)

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

    ! Set the base vec
    call setwVec(wVec)
    call formFunction_mf(ctx, wVec, baseRes, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call MatMFFDSetBase(dRdW, wVec, baseRes, ierr)
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
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none
#else
  implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#endif

  Mat  dRdw, dRdwPre
  Vec wVec, rVec, deltaW, baseRes
  KSP  ANK_KSP

  ! Turb KSP related PETSc objects
  Mat  dRdwTurb, dRdwPreTurb
  Vec wVecTurb, rVecTurb, deltaWTurb, baseResTurb
  KSP  ANK_KSPTurb

  PetscFortranAddr   ctx(1)

  ! Options for ANK Solver
  logical :: useANKSolver
  integer(kind=intType) :: ANK_jacobianLag
  integer(kind=intType) :: ANK_subSpace
  integer(kind=intType) :: ANK_maxIter
  integer(kind=intType) :: ANK_asmOverlap
  integer(kind=intType) :: ANK_iluFill
  integer(kind=intType) :: ANK_innerPreConIts
  integer(kind=intType) :: ANK_outerPreConIts
  real(kind=realType)   :: ANK_rtol
  real(kind=realType)   :: ANK_linResMax
  real(kind=realType)   :: ANK_switchTol
  real(kind=realType)   :: ANK_divTol = 10
  logical :: ANK_useTurbDADI
  real(kind=realType) :: ANK_turbcflscale
  logical :: ANK_useFullVisc
  logical :: ANK_ADPC
  logical :: ANK_turbDebug
  logical :: ANK_useMatrixFree
  integer(kind=intType) :: ANK_nsubIterTurb

  ! Misc variables
  real(kind=realType) :: ANK_CFL, ANK_CFL0, ANK_CFLLimit, ANK_CFLFactor, ANK_CFLCutback
  real(kind=realType) :: ANK_CFLMin0, ANK_CFLMin, ANK_CFLMinBase, ANK_CFLExponent
  real(kind=realType) :: ANK_stepMin, ANK_StepFactor, ANK_constCFLStep
  real(kind=realType) :: ANK_secondOrdSwitchTol, ANK_coupledSwitchTol
  real(kind=realType) :: ANK_physLSTol, ANK_unstdyLSTol
  real(kind=realType) :: ANK_pcUpdateTol
  real(kind=realType) :: lambda
  logical :: ANK_solverSetup=.False.
  integer(kind=intTYpe) :: ANK_iter
  integer(kind=intType) :: nState
  real(kind=alwaysRealType) :: totalR_old, totalR_pcUpdate ! for recording the previous residual
  real(kind=alwaysRealType) :: rtolLast, linResOld ! for recording the previous relativel tolerance for Eisenstat-Walker
  logical :: ANK_useDissApprox

  ! Turb KSP related modifications
  logical :: ANK_coupled=.False.
  logical :: ANK_turbSetup=.False.
  integer(kind=intType) :: ANK_iterTurb, nStateTurb
  real(kind=realType) :: lambdaTurb
  real(kind=alwaysRealType) :: linResOldTurb, ANK_physLSTolTurb

contains

  subroutine setupANKsolver

    ! Setup the PETSc objects for the Newton-Krylov
    ! solver. destroyNKsolver can be used to destroy the objects created
    ! in this function

    use constants
    use stencils, only : euler_PC_stencil, N_euler_PC
    use communication, only : adflow_comm_world, myid
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputIteration, only : useLinResMonitor
    use inputPhysics, only : equations
    use flowVarRefState, only : nw, viscous, nwf, nt1, nt2
    use ADjointVars , only: nCellsLocal
    use NKSolver, only : destroyNKSolver, linearResidualMonitor
    use utils, only : EChk
    use adjointUtils, only : myMatCreate, statePreAllocation
    use inputadjoint, only : precondtype
    use agmg, only : setupAGMG
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

       ! Determine if we are in coupled mode
       if (ANK_coupled) then
          nState = nw
       else
          nState = nwf
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

       call VecDuplicate(wVec, baseRes, ierr)
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

       if (preCondType == 'mg') then
          call setupAGMG(drdwpre, nDimW/nState, nState)
       end if

       !  Create the linear solver context
       call KSPCreate(ADFLOW_COMM_WORLD, ANK_KSP, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set operators for the solver
       if (ANK_useMatrixFree) then
           ! Matrix free drdw
           call KSPSetOperators(ANK_KSP, dRdw, dRdwPre, ierr)
       else
           ! Matrix based drdw = drdwpre
           call KSPSetOperators(ANK_KSP, dRdwPre, dRdwPre, ierr)
       end if
       call EChk(ierr, __FILE__, __LINE__)

       if (useLinResMonitor) then

#if PETSC_VERSION_GE(3,8,0)
          ! This is probably wrong. NO petsc_null_context
          call KSPMonitorSet(ANK_KSP, LinearResidualMonitor, PETSC_NULL_FUNCTION, &
               PETSC_NULL_FUNCTION, ierr)
#else
          call KSPMonitorSet(ANK_KSP, LinearResidualMonitor, PETSC_NULL_OBJECT, &
               PETSC_NULL_FUNCTION, ierr)

#endif
          call EChk(ierr, __FILE__, __LINE__)
       end if

       ANK_solverSetup = .True.
       ANK_iter = 0
       ANK_useDissApprox = .False.

       ! Check if we need to set up the Turb KSP
       if ((.not. ANK_coupled) .and. (.not. ANK_useTurbDADI) .and. equations==RANSEquations) then
           nStateTurb = nt2-nt1+1

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

           call VecDuplicate(wVecTurb, baseResTurb, ierr)
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
           call KSPCreate(ADFLOW_COMM_WORLD, ANK_KSPTurb, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ! Set operators for the solver
           if (ANK_useMatrixFree) then
              ! Matrix free
              call KSPSetOperators(ANK_KSPTurb, dRdwTurb, dRdwPreTurb, ierr)
           else
              ! Matrix based
              call KSPSetOperators(ANK_KSPTurb, dRdwPreTurb, dRdwPreTurb, ierr)
           end if
           call EChk(ierr, __FILE__, __LINE__)

           ANK_turbSetup = .True.
           ANK_iterTurb = 0
       end if
    end if

  end subroutine setupANKsolver

  subroutine FormJacobianANK

    use constants
    use flowVarRefState, only : nw, nwf, nt1, nt2
    use blockPointers, only : nDom, volRef, il, jl, kl, w, dw, dtl, globalCell, iblank
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputIteration, only : turbResScale
    use inputADjoint, only : viscPC
    use inputDiscretization, only : approxSA
    use iteration, only : totalR0, totalR
    use utils, only : EChk, setPointers
    use adjointUtils, only :setupStateResidualMatrix, setupStandardKSP, setupStandardMultigrid
    use communication
    use agmg, only : setupShellPC, destroyShellPC, applyShellPC, agmgLevels, coarseIndices, A
    use inputadjoint, only : precondtype
    implicit none

    ! Local Variables
    character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
    integer(kind=intType) ::ierr
    logical :: useAD, usePC, useTranspose, useObjective, tmp, frozenTurb
    real(kind=realType) :: dtinv, rho
    integer(kind=intType) :: i, j, k, l, ii, irow, nn, sps, outerPreConIts, subspace, lvl
    integer(kind=intType), dimension(2:10) :: coarseRows
    real(kind=realType), dimension(:,:), allocatable :: blk
    logical :: useCoarseMats
    PC shellPC

    if (preCondType == 'mg') then
       useCoarseMats = .True.
    else
       useCoarseMats = .False.
    end if

    ! Assemble the approximate PC (fine leve, level 1)
    useAD = ANK_ADPC
    frozenTurb = (.not. ANK_coupled)
    usePC = .True.
    useTranspose = .False.
    useObjective = .False.
    tmp = viscPC ! Save what is in viscPC and set to the NKvarible
    viscPC = .False.

    if (totalR > ANK_secondOrdSwitchTol*totalR0) &
       approxSA = .True.

    ! Create the preconditoner matrix
    call setupStateResidualMatrix(dRdwPre, useAD, usePC, useTranspose, &
         useObjective, frozenTurb, 1_intType, useCoarseMats=useCoarseMats)

    ! Reset saved value
    viscPC = tmp
    approxSA = .False.

    ! Add the contribution from the time step term

    ! Generic block to use while setting values
    allocate(blk(nState, nState))

    ! Zero the block once, since the previous entries will be overwritten
    ! for each cell, and zero entries will remain zero.
    blk = zero

    if (.not. ANK_coupled) then
       ! For the segragated solver, only calculate the time step for flow variables
       do nn=1, nDom
          do sps=1, nTimeIntervalsSpectral
             call setPointers(nn,1_intType,sps)
             do k=2, kl
                do j=2, jl
                   do i=2, il
                      ! Calculate one over time step for this cell. Multiply
                      ! the dtl by cell volume to get the actual time step
                      ! required for a CFL of one, then multiply with the
                      ! actual cfl number in the solver
                      dtinv = one/(ANK_CFL * dtl(i,j,k) * volRef(i,j,k))

                      ! We need to convert the momentum residuals to velocity
                      ! residuals to get the desired effect from time steps.
                      ! To do this, save a "pseudo" jacobian for this cell,
                      ! that has dU/du, where U is the vector of conservative
                      ! variables, and u are the primitive variables. For this
                      ! jacobian, only the velocity entries are modified,
                      ! since ADflow saves density, velocities and total
                      ! energy in the state vector w(:,:,:,:).

                      ! Density and energy updates are unchanged.
                      blk(iRho, iRho) = dtinv
                      blk(iRhoE, iRhoE) = dtinv

                      ! save the density
                      rho = w(i,j,k,iRho)

                      ! x-velocity
                      blk(ivx, iRho) = w(i,j,k,ivx)*dtinv
                      blk(ivx, ivx)  = rho*dtinv

                      ! y-velocity
                      blk(ivy, iRho) = w(i,j,k,ivy)*dtinv
                      blk(ivy, ivy)  = rho*dtinv

                      ! z-velocity
                      blk(ivz, iRho) = w(i,j,k,ivz)*dtinv
                      blk(ivz, ivz)  = rho*dtinv

                      ! get the global cell index
                      irow = globalCell(i, j, k)

                      if (useCoarseMats) then
                         do lvl=1, agmgLevels-1
                            coarseRows(lvl+1) = coarseIndices(nn, lvl)%arr(i, j, k)
                         end do
                      end if

                      ! Add the contribution to the matrix in PETSc
                      call setBlock()
                   end do
                end do
             end do
          end do
       end do
    else
       ! For the coupled solver, CFL number for the turbulent variable needs scaling
       ! because the residuals are scaled, and additional scaling of the time step
       ! for the turbulence variable might be required.
       ii = 1
       do nn=1, nDom
          do sps=1, nTimeIntervalsSpectral
             call setPointers(nn,1_intType,sps)
             do k=2, kl
                do j=2, jl
                   do i=2, il
                      ! See the comment for the same calculation above
                      dtinv = one/(ANK_CFL * dtl(i,j,k) * volRef(i,j,k))

                      ! We need to convert the momentum residuals to velocity
                      ! residuals to get the desired effect from time steps.
                      ! To do this, save a "pseudo" jacobian for this cell,
                      ! that has dU/du, where U is the vector of conservative
                      ! variables, and u are the primitive variables. For this
                      ! jacobian, only the velocity entries are modified,
                      ! since ADflow saves density, velocities and total
                      ! energy in the state vector w(:,:,:,:).

                      ! Density update is unchanged.
                      blk(iRho, iRho) = dtinv

                      ! save the density
                      rho = w(i,j,k,iRho)

                      ! x-velocity
                      blk(ivx, iRho) = w(i,j,k,ivx)*dtinv
                      blk(ivx, ivx)  = rho*dtinv

                      ! y-velocity
                      blk(ivy, iRho) = w(i,j,k,ivy)*dtinv
                      blk(ivy, ivy)  = rho*dtinv

                      ! z-velocity
                      blk(ivz, iRho) = w(i,j,k,ivz)*dtinv
                      blk(ivz, ivz)  = rho*dtinv

                      ! Energy update is unchanged
                      blk(iRhoE, iRhoE) = dtinv

                      ! For the turbulence variable, additionally scale the cfl.
                      ! turbresscale is required because the turbulent residuals
                      ! are scaled with it. Furthermore, the turbulence variable
                      ! can get a different CFL number. Scale it by turbCFLScale
                      blk(nt1, nt1) = dtinv*turbResScale(1)/ANK_turbCFLScale

                      ! get the global cell index
                      irow = globalCell(i, j, k)

                      ! Add the contribution to the matrix in PETSc
                      call setBlock()
                   end do
                end do
             end do
          end do
       end do
    end if

    ! PETSc Matrix Assembly begin
    call MatAssemblyBegin(dRdwPre, MAT_FINAL_ASSEMBLY, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Setup KSP Options
    preConSide = 'right'
    localPCType = 'ilu'
    kspObjectType = 'gmres'
    globalPCType = 'asm'
    localOrdering = 'rcm'
    outerPreConIts = ank_outerPreconIts

    ! Setup the KSP using the same code as used for the adjoint
    if (ank_subspace < 0) then
       subspace = ANK_maxIter
    else
       subspace = ANK_subspace
    end if

    ! de-allocate the generic block
    deallocate(blk)

    ! Complete the matrix assembly.
    call MatAssemblyEnd  (dRdwPre, MAT_FINAL_ASSEMBLY, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (useCoarseMats) then
       do lvl=2, agmgLevels
          call MatAssemblyBegin(A(lvl), MAT_FINAL_ASSEMBLY, ierr)
          call EChk(ierr, __FILE__, __LINE__)
          call MatAssemblyEnd(A(lvl), MAT_FINAL_ASSEMBLY, ierr)
          call EChk(ierr, __FILE__, __LINE__)

       end do
    end if

    if (PreCondType == 'asm') then
       ! Run the super-dee-duper function to setup the ksp object:

       call setupStandardKSP(ANK_KSP, kspObjectType, subSpace, &
            preConSide, globalPCType, ANK_asmOverlap, outerPreConIts, localPCType, &
            localOrdering, ANK_iluFill, ANK_innerPreConIts)
    else if (PreCondType == 'mg') then

       ! Setup the MG preconditioner!
       call setupStandardMultigrid(ANK_KSP, kspObjectType, subSpace, &
            preConSide, ANK_asmOverlap, outerPreConIts, &
            localOrdering, ANK_iluFill)
    end if

    ! Don't do iterative refinement for the NKSolver.
    call KSPGMRESSetCGSRefinementType(ANK_KSP, &
         KSP_GMRES_CGS_REFINE_NEVER, ierr)
    call EChk(ierr, __FILE__, __LINE__)

  contains
    subroutine setBlock()
      ! This subroutine is used to set the diagonal time stepping terms
      ! for the Jacobians in ANK. It is only used to set diagonal blocks

      implicit none

      call MatSetValuesBlocked(dRdwPre, 1, irow, 1, irow, blk, &
           ADD_VALUES, ierr)
      call EChk(ierr, __FILE__, __LINE__)

      ! Extension for setting coarse grids:
      if (useCoarseMats) then
         do lvl=2, agmgLevels
            call MatSetValuesBlocked(A(lvl), 1, coarseRows(lvl), 1, coarseRows(lvl), &
                 blk, ADD_VALUES, ierr)
         end do
      end if


    end subroutine setBlock
  end subroutine FormJacobianANK

  subroutine FormJacobianANKTurb

    use constants
    use flowVarRefState, only : nw, nwf, nt1, nt2
    use blockPointers, only : nDom, volRef, il, jl, kl, w, dw, dtl, globalCell
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputIteration, only : turbResScale
    use inputADjoint, only : viscPC
    use inputDiscretization, only : approxSA
    use iteration, only : totalR0, totalR
    use utils, only : EChk, setPointers
    use adjointUtils, only :setupStateResidualMatrix, setupStandardKSP
    use communication
    implicit none

    ! Local Variables
    character(len=maxStringLen) :: preConSide, localPCType, kspObjectType, globalPCType, localOrdering
    integer(kind=intType) ::ierr
    logical :: useAD, usePC, useTranspose, useObjective, tmp, frozenTurb
    real(kind=realType) :: dtinv, rho
    integer(kind=intType) :: i, j, k, l, l1, ii, irow, nn, sps, outerPreConIts, subspace
    real(kind=realType), dimension(:,:), allocatable :: blk

    ! Assemble the approximate PC (fine leve, level 1)
    useAD = ANK_ADPC
    frozenTurb = .False.
    usePC = .True.
    useTranspose = .False.
    useObjective = .False.
    tmp = viscPC ! Save what is in viscPC and set to the NKvarible
    viscPC = .False.

    if (totalR > ANK_secondOrdSwitchTol*totalR0) &
       approxSA = .True.

    ! Create the preconditoner matrix
    call setupStateResidualMatrix(dRdwPreTurb, useAD, usePC, useTranspose, &
         useObjective, frozenTurb, 1_intType, .True.)

    ! Reset saved value
    viscPC = tmp
    approxSA = .False.

    ! Add the contribution from the time step term

    ! Generic block to use while setting values
    allocate(blk(nStateTurb, nStateTurb))

    ! Zero the block once, since the previous entries will be overwritten
    ! for each cell, and zero entries will remain zero.
    blk = zero

    ! For the coupled solver, CFL number for the turbulent variable needs scaling
    ! because the residuals are scaled, and additional scaling of the time step
    ! for the turbulence variable might be required.
    ii = 1
    do nn=1, nDom
        do sps=1, nTimeIntervalsSpectral
            call setPointers(nn,1_intType,sps)
            do k=2, kl
                do j=2, jl
                    do i=2, il

                      ! See the comment for the same calculation above
                      dtinv = one/(ANK_CFL * dtl(i,j,k) * volRef(i,j,k))

                      do l=nt1, nt2

                        ! l1 is just l that starts with 1 on the turb variables
                        l1 = l-nt1+1

                        ! For the turbulence variable, additionally scale the cfl.
                        ! turbresscale is required because the turbulent residuals
                        ! are scaled with it. Furthermore, the turbulence variable
                        ! can get a different CFL number. Scale it by turbCFLScale
                        blk(l1, l1) = dtinv*turbResScale(l1)/ANK_turbCFLScale
                      end do

                      ! get the global cell index
                      irow = globalCell(i, j, k)

                      ! Add the contribution to the matrix in PETSc
                      call setBlock()
                    end do
                end do
            end do
        end do
    end do

    ! PETSc Matrix Assembly begin
    call MatAssemblyBegin(dRdwPreTurb, MAT_FINAL_ASSEMBLY, ierr)
    call EChk(ierr, __FILE__, __LINE__)

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

    ! de-allocate the generic block
    deallocate(blk)

    ! Complete the matrix assembly.
    call MatAssemblyEnd  (dRdwPreTurb, MAT_FINAL_ASSEMBLY, ierr)
    call EChk(ierr, __FILE__, __LINE__)


    call setupStandardKSP(ANK_KSPTurb, kspObjectType, subSpace, &
         preConSide, globalPCType, ANK_asmOverlap, outerPreConIts, localPCType, &
         localOrdering, ANK_iluFill, ANK_innerPreConIts)

    ! Don't do iterative refinement for the NKSolver.
    call KSPGMRESSetCGSRefinementType(ANK_KSPTurb, &
         KSP_GMRES_CGS_REFINE_NEVER, ierr)
    call EChk(ierr, __FILE__, __LINE__)

  contains
    subroutine setBlock()
      ! This subroutine is used to set the diagonal time stepping terms
      ! for the Jacobians in ANK. It is only used to set diagonal blocks

      implicit none

      call MatSetValuesBlocked(dRdwPreTurb, 1, irow, 1, irow, blk, &
           ADD_VALUES, ierr)
      call EChk(ierr, __FILE__, __LINE__)

    end subroutine setBlock
  end subroutine FormJacobianANKTurb

  subroutine FormFunction_mf(ctx, inVec, rVec, ierr)

    ! This is the function used for the matrix-free matrix-vector products
    ! for the GMRES solver used in ANK

    use constants
    use blockPointers, only : nDom, volRef, il, jl, kl, dw, dtl
    use inputtimespectral, only : nTimeIntervalsSpectral
    use inputIteration, only : turbResScale
    use flowvarrefstate, only : nwf, nt1, nt2
    use NKSolver, only : setRvec
    use utils, only : setPointers, EChk
    use blockette, only : blocketteRes
    implicit none

    ! PETSc Variables
    PetscFortranAddr ctx(*)
    Vec     inVec, rVec
    real(kind=realType) :: dtinv, rho
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii, iiRho
    real(kind=realType),pointer :: rvec_pointer(:)
    real(kind=realType),pointer :: invec_pointer(:)
    real(kind=realType),pointer :: wvec_pointer(:)
    logical :: useViscApprox

    ! get the input vector
    call setWANK(inVec,1,nState)

    ! determine if we want the approximate viscous fluxes
    useViscApprox = (.not. ANK_useFullVisc) .and. ANK_useDissApprox

    ! Determine if we want the turb residuals
    call blocketteRes(useDissApprox=ANK_useDissApprox, useViscApprox=useViscApprox, &
         useTurbRes=ANK_coupled, useStoreWall=.False.)

    ! Copy the residuals to rVec in petsc
    if (ANK_coupled) then
       call setRVec(rVec)
    else
       call setRVecANK(rVec)
    end if

    ! Add the contribution from the time stepping term

    call VecGetArrayF90(rVec,rvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! inVec contains the perturbed state vector
    call VecGetArrayReadF90(inVec,invec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Also read the wVec to access the un-perturbed state vector.
    call VecGetArrayReadF90(wVec,wvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    if (.not. ANK_coupled) then
       ! Only flow variables
       ii = 1
       do nn=1, nDom
          do sps=1, nTimeIntervalsSpectral
             call setPointers(nn,1_intType,sps)
             ! read the density residuals and set local CFL
             do k=2, kl
                do j=2, jl
                   do i=2, il
                      dtinv = one/(ANK_CFL * dtl(i,j,k) * volRef(i,j,k))

                      ! Update the first entry in this block, corresponds to
                      ! density. Also save this density value.
                      rvec_pointer(ii) = rvec_pointer(ii) + invec_pointer(ii)*dtinv
                      rho = wvec_pointer(ii)
                      iirho = ii
                      ii = ii + 1

                      ! updates 2nd-4th are velocities. They need to get converted
                      ! to momentum residuals.

                      rvec_pointer(ii) = rvec_pointer(ii) + dtinv*( &
                           wvec_pointer(ii)*invec_pointer(iiRho)+&
                           rho * invec_pointer(ii))
                      ii = ii + 1

                      rvec_pointer(ii) = rvec_pointer(ii) + dtinv*( &
                           wvec_pointer(ii)*invec_pointer(iiRho)+&
                           rho * invec_pointer(ii))
                      ii = ii + 1

                      rvec_pointer(ii) = rvec_pointer(ii) + dtinv*( &
                           wvec_pointer(ii)*invec_pointer(iiRho)+&
                           rho * invec_pointer(ii))
                      ii = ii + 1

                      ! Finally energy gets the same update
                      rvec_pointer(ii) = rvec_pointer(ii) + invec_pointer(ii)*dtinv
                      ii = ii + 1
                   end do
                end do
             end do
          end do
       end do
    else
       ! Include time step for turbulence
       ii = 1
       do nn=1, nDom
          do sps=1, nTimeIntervalsSpectral
             call setPointers(nn,1_intType,sps)
             ! read the density residuals and set local CFL
             do k=2, kl
                do j=2, jl
                   do i=2, il
                      dtinv = one/(ANK_CFL * dtl(i,j,k) * volRef(i,j,k))

                      ! Update the first entry in this block, corresponds to
                      ! density. Also save this density value.
                      rvec_pointer(ii) = rvec_pointer(ii) + invec_pointer(ii)*dtinv
                      rho = wvec_pointer(ii)
                      iirho = ii
                      ii = ii + 1

                      ! updates 2nd-4th are velocities. They need to get converted
                      ! to momentum residuals.

                      rvec_pointer(ii) = rvec_pointer(ii) + dtinv*( &
                           wvec_pointer(ii)*invec_pointer(iiRho)+&
                           rho * invec_pointer(ii))
                      ii = ii + 1

                      rvec_pointer(ii) = rvec_pointer(ii) + dtinv*( &
                           wvec_pointer(ii)*invec_pointer(iiRho)+&
                           rho * invec_pointer(ii))
                      ii = ii + 1

                      rvec_pointer(ii) = rvec_pointer(ii) + dtinv*( &
                           wvec_pointer(ii)*invec_pointer(iiRho)+&
                           rho * invec_pointer(ii))
                      ii = ii + 1

                      ! energy gets the same update
                      rvec_pointer(ii) = rvec_pointer(ii) + invec_pointer(ii)*dtinv
                      ii = ii + 1

                      do l=nt1, nt2
                         ! turbulence variable needs additional scaling, and it may
                         ! get a different CFL number
                         rvec_pointer(ii) = rvec_pointer(ii) + invec_pointer(ii)* &
                              dtinv*turbResScale(l-nt1+1)/ANK_turbCFLScale
                         ii = ii + 1
                      end do
                   end do
                end do
             end do
          end do
       end do
    end if

    call VecRestoreArrayF90(rVec, rvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecRestoreArrayReadF90(wVec, wvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecRestoreArrayReadF90(inVec, invec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! We don't check an error here, so just pass back zero
    ierr = 0

  end subroutine FormFunction_mf

  subroutine FormFunction_mf_turb(ctx, inVec, rVec, ierr)

    ! This is the function used for the matrix-free matrix-vector products
    ! for the GMRES solver used in ANK

    use constants
    use blockPointers, only : nDom, volRef, il, jl, kl, dw, dtl
    use inputtimespectral, only : nTimeIntervalsSpectral
    use inputIteration, only : turbResScale
    use flowvarrefstate, only : nwf, nt1, nt2
    use NKSolver, only : setRvec
    use utils, only : setPointers, EChk
    use blockette, only : blocketteRes
    implicit none

    ! PETSc Variables
    PetscFortranAddr ctx(*)
    Vec     inVec, rVec
    real(kind=realType) :: dtinv, rho
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii, iiRho
    real(kind=realType),pointer :: rvec_pointer(:)
    real(kind=realType),pointer :: invec_pointer(:)

    ! get the input vector
    call setWANK(inVec,nt1,nt2)

    call blocketteRes(useFlowRes=.False., useStoreWall=.False.)
    call setRVecANKTurb(rVec)

    ! Add the contribution from the time stepping term

    call VecGetArrayF90(rVec,rvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! inVec contains the perturbed state vector
    call VecGetArrayReadF90(inVec,invec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Include time step for turbulence
    ii = 1
    do nn=1, nDom
        do sps=1, nTimeIntervalsSpectral
            call setPointers(nn,1_intType,sps)
            ! read the density residuals and set local CFL
            do k=2, kl
                do j=2, jl
                    do i=2, il
                     ! needs to be modified
                        dtinv = one/(ANK_CFL * dtl(i,j,k) * volRef(i,j,k))

                        do l=nt1, nt2
                           ! turbulence variable needs additional scaling, and it may
                           ! get a different CFL number
                           rvec_pointer(ii) = rvec_pointer(ii) + invec_pointer(ii)* &
                           dtinv*turbResScale(l-nt1+1)/ANK_turbCFLScale
                           ii = ii + 1
                        end do
                    end do
                end do
            end do
        end do
    end do


    call VecRestoreArrayF90(rVec, rvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecRestoreArrayReadF90(inVec, invec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! We don't check an error here, so just pass back zero
    ierr = 0

  end subroutine FormFunction_mf_turb

  subroutine computeUnsteadyResANK(omega)

    ! This routine calculates the unsteady residual in a given iteration.
    ! It needs the following variables/vectors:
    !
    !   omega:      This is the step size taken in the last update to the state
    !   deltaW:     Vector that contains the full update given from the
    !               Newton/Euler iteration.
    !   w(:,:,:,:): Should contain the updated state with the given step size
    !               lambdaLS and given update deltaW
    !   ANK_CFL:    The CFL number used for this non-linear iteration
    !   dtl:        Array containing time step values giving a CFL number of 1
    !               on each cell.
    !
    ! The routine calculates the unsteady residual and leaves the result in
    ! rVec, which was previously used to keep the steady residual only. This
    ! is done because the norm of this vector can easily be calculated with
    ! PETSc, however, after the line search, the rVec vector needs to be
    ! updated to contain only the steady state residuals. This can be done with
    ! setRVecANK/setRVec, with a dw(:,:,:,:) that is also up to date.

    use constants
    use blockPointers, only : nDom, volRef, il, jl, kl, w, dw, dtl
    use inputtimespectral, only : nTimeIntervalsSpectral
    use inputIteration, only : turbResScale
    use flowvarrefstate, only : nwf, nt1, nt2
    use NKSolver, only : setRvec
    use utils, only : setPointers, EChk
    use blockette, only : blocketteRes
    implicit none

    real(kind=realType), intent(in) :: omega

    real(kind=realType) :: dtinv, rho, uu, vv, ww
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii, iiRho
    real(kind=realType),pointer :: rvec_pointer(:)
    real(kind=realType),pointer :: dvec_pointer(:)

    ! Calculate the steady residuals
    call blocketteRes(useTurbRes=ANK_coupled)
    call setRVecANK(rVec)

    ! Add the contribution from the time stepping term

    call VecGetArrayF90(rVec,rvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! deltaW contains the full update to the state
    call VecGetArrayReadF90(deltaW,dvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    if (.not. ANK_coupled) then
       ! Only flow variables
       ii = 1
       do nn=1, nDom
          do sps=1, nTimeIntervalsSpectral
             call setPointers(nn,1_intType,sps)
             ! read the density residuals and set local CFL
             do k=2, kl
                do j=2, jl
                   do i=2, il
                      dtinv = one/(ANK_CFL * dtl(i,j,k) * volRef(i,j,k))

                      ! Update the first entry in this block, corresponds to
                      ! density. Also save this density value.
                      rvec_pointer(ii) = rvec_pointer(ii) - omega*dtinv*dvec_pointer(ii)

                      ! calculate the density in the previous non-linear iteration
                      rho = w(i,j,k,iRho) + omega*dvec_pointer(ii)
                      iiRho = ii
                      ii = ii + 1

                      ! updates 2nd-4th are velocities. They need to get converted
                      ! to momentum residuals.

                      ! Calculate the u velocity in the previous non-linear iter.
                      uu = w(i,j,k,ivx) + omega*dvec_pointer(ii)

                      rvec_pointer(ii) = rvec_pointer(ii) - omega*dtinv*( &
                           uu*dvec_pointer(iiRho)+&
                           rho * dvec_pointer(ii))
                      ii = ii + 1

                      vv = w(i,j,k,ivx) + omega*dvec_pointer(ii)

                      rvec_pointer(ii) = rvec_pointer(ii) - omega*dtinv*( &
                           vv*dvec_pointer(iiRho)+&
                           rho * dvec_pointer(ii))
                      ii = ii + 1

                      ww = w(i,j,k,ivx) + omega*dvec_pointer(ii)

                      rvec_pointer(ii) = rvec_pointer(ii) - omega*dtinv*( &
                           ww*dvec_pointer(iiRho)+&
                           rho * dvec_pointer(ii))
                      ii = ii + 1

                      ! Finally energy gets the same update
                      rvec_pointer(ii) = rvec_pointer(ii) - omega*dtinv*dvec_pointer(ii)
                      ii = ii + 1
                   end do
                end do
             end do
          end do
       end do
    else
       ! Include time step for turbulence
       ii = 1
       do nn=1, nDom
          do sps=1, nTimeIntervalsSpectral
             call setPointers(nn,1_intType,sps)
             ! read the density residuals and set local CFL
             do k=2, kl
                do j=2, jl
                   do i=2, il
                      dtinv = one/(ANK_CFL * dtl(i,j,k) * volRef(i,j,k))

                      ! Update the first entry in this block, corresponds to
                      ! density. Also save this density value.
                      rvec_pointer(ii) = rvec_pointer(ii) - omega*dtinv*dvec_pointer(ii)

                      ! calculate the density in the previous non-linear iteration
                      rho = w(i,j,k,iRho) + omega*dvec_pointer(ii)
                      iiRho = ii
                      ii = ii + 1

                      ! updates 2nd-4th are velocities. They need to get converted
                      ! to momentum residuals.

                      ! Calculate the u velocity in the previous non-linear iter.
                      uu = w(i,j,k,ivx) + omega*dvec_pointer(ii)

                      rvec_pointer(ii) = rvec_pointer(ii) - omega*dtinv*( &
                           uu*dvec_pointer(iiRho)+&
                           rho * dvec_pointer(ii))
                      ii = ii + 1

                      vv = w(i,j,k,ivx) + omega*dvec_pointer(ii)

                      rvec_pointer(ii) = rvec_pointer(ii) - omega*dtinv*( &
                           vv*dvec_pointer(iiRho)+&
                           rho * dvec_pointer(ii))
                      ii = ii + 1

                      ww = w(i,j,k,ivx) + omega*dvec_pointer(ii)

                      rvec_pointer(ii) = rvec_pointer(ii) - omega*dtinv*( &
                           ww*dvec_pointer(iiRho)+&
                           rho * dvec_pointer(ii))
                      ii = ii + 1

                      ! Finally energy gets the same update
                      rvec_pointer(ii) = rvec_pointer(ii) - omega*dtinv*dvec_pointer(ii)
                      ii = ii + 1

                      do l=nt1, nt2
                         ! turbulence variable needs additional scaling, and it may
                         ! get a different CFL number
                         rvec_pointer(ii) = rvec_pointer(ii) - omega*dvec_pointer(ii)* &
                              dtinv*turbResScale(l-nt1+1)/ANK_turbCFLScale
                         ii = ii + 1
                      end do
                   end do
                end do
             end do
          end do
       end do
    end if

    call VecRestoreArrayF90(rVec, rvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecRestoreArrayReadF90(deltaW, dvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! We don't check an error here, so just pass back zero
    ierr = 0

  end subroutine computeUnsteadyResANK

  subroutine computeUnsteadyResANKTurb(omega)

    ! This routine calculates the unsteady residual in a given iteration.
    ! It needs the following variables/vectors:
    !
    !   omega:      This is the step size taken in the last update to the state
    !   deltaWTurb: Vector that contains the full update given from the
    !               Newton/Euler iteration.
    !   w(:,:,:,:): Should contain the updated state with the given step size
    !               lambdaLS and given update deltaW
    !   ANK_CFL:    The CFL number used for this non-linear iteration
    !   dtl:        Array containing time step values giving a CFL number of 1
    !               on each cell.
    !
    ! The routine calculates the unsteady residual and leaves the result in
    ! rVecTurb, which was previously used to keep the steady residual only. This
    ! is done because the norm of this vector can easily be calculated with
    ! PETSc, however, after the line search, the rVec vector needs to be
    ! updated to contain only the steady state residuals. This can be done with
    ! setRVecANK/setRVec, with a dw(:,:,:,:) that is also up to date.

    use constants
    use blockPointers, only : nDom, volRef, il, jl, kl, w, dw, dtl
    use inputtimespectral, only : nTimeIntervalsSpectral
    use inputIteration, only : turbResScale
    use flowvarrefstate, only : nwf, nt1, nt2
    use NKSolver, only : setRvec
    use utils, only : setPointers, EChk
    use blockette, only : blocketteRes
    implicit none

    real(kind=realType), intent(in) :: omega

    real(kind=realType) :: dtinv, rho, uu, vv, ww
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii, iiRho
    real(kind=realType),pointer :: rvec_pointer(:)
    real(kind=realType),pointer :: dvec_pointer(:)

    ! Calculate the steady residuals
    call blocketteRes(useFlowRes=.False.)
    call setRVecANKTurb(rVecTurb)

    ! Add the contribution from the time stepping term

    call VecGetArrayF90(rVecTurb,rvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! deltaW contains the full update to the state
    call VecGetArrayReadF90(deltaWTurb,dvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Include time step for turbulence
    ii = 1
    do nn=1, nDom
        do sps=1, nTimeIntervalsSpectral
            call setPointers(nn,1_intType,sps)
            ! read the density residuals and set local CFL
            do k=2, kl
                do j=2, jl
                    do i=2, il
                        dtinv = one/(ANK_CFL * dtl(i,j,k) * volRef(i,j,k))

                        do l=nt1, nt2
                           ! turbulence variable needs additional scaling, and it may
                           ! get a different CFL number
                           rvec_pointer(ii) = rvec_pointer(ii) - omega*dvec_pointer(ii)* &
                                dtinv*turbResScale(l-nt1+1)/ANK_turbCFLScale
                           ii = ii + 1
                        end do
                    end do
                end do
            end do
        end do
    end do

    call VecRestoreArrayF90(rVecTurb, rvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecRestoreArrayReadF90(deltaWTurb, dvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! We don't check an error here, so just pass back zero
    ierr = 0

  end subroutine computeUnsteadyResANKTurb

  subroutine destroyANKsolver

    ! Destroy all the PETSc objects for the Newton-Krylov
    ! solver.

    use constants
    use utils, only : EChk
    use agmg, only : destroyAGMG
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

       call VecDestroy(baseRes, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call KSPDestroy(ANK_KSP, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call destroyAGMG()

       ANK_SolverSetup = .False.

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

           call VecDestroy(baseResTurb, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           call KSPDestroy(ANK_KSPTurb, ierr)
           call EChk(ierr, __FILE__, __LINE__)

           ANK_turbSetup = .False.
       end if
    end if
  end subroutine destroyANKsolver

  subroutine setWVecANK(wVec,lStart,lEnd)
    ! Set the current FLOW variables in the PETSc Vector

    use constants
    use blockPointers, only : nDom, il, jl, kl, w
    use inputtimespectral, only : ntimeIntervalsSpectral
    use utils, only : setPointers, EChk
    implicit none

    Vec   wVec
    integer(kind=intType),intent(in) :: lStart, lEnd
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
                   do l=lStart, lEnd
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
    real(kind=realType),pointer :: rvec_pointer(:)
    real(Kind=realType) :: ovv
    call VecGetArrayF90(rVecTurb,rvec_pointer,ierr)
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
                      rvec_pointer(ii) = dw(i, j, k, l)*ovv*turbResScale(1)
                   end do
                end do
             end do
          end do
       end do
    end do

    call VecRestoreArrayF90(rVecTurb, rvec_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine setRVecANKTurb

  subroutine setWANK(wVec,lStart,lEnd)
    ! Get the updated solution from the PETSc Vector

    use constants
    use blockPointers, only : nDom, vol, il, jl, kl, w
    use inputtimespectral, only : nTimeIntervalsSpectral
    use utils, only : setPointers, EChk
    implicit none

    Vec  wVec
    integer(kind=intType),intent(in) :: lStart, lEnd
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
                   do l=lStart, lEnd
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

  subroutine physicalityCheckANK(lambdaP)

    use constants
    use blockPointers, only : ndom, il, jl, kl
    use flowVarRefState, only : nw, nwf, nt1, nt2
    use inputtimespectral, only : nTimeIntervalsSpectral
    use utils, only : setPointers, EChk, myisnan
    use communication, only : ADflow_comm_world
    implicit none

    ! input variable
    real(kind=realType) , intent(inout) :: lambdaP

    ! local variables
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
    real(kind=realType), pointer :: wvec_pointer(:)
    real(kind=realType), pointer :: dvec_pointer(:)
    real(kind=alwaysRealType) :: lambdaL ! L is for local
    real(kind=realType) :: ratio


    ! Determine the maximum step size that would yield
    ! a maximum change of 10% in density, total energy,
    ! and turbulence variable after a KSP solve.

    ! Initialize the local step size as ANK_stepFactor
    ! because the initial step is likely to be equal to this.
    lambdaL = real(lambdaP)

    ! First we need to read both the update and the state
    ! from PETSc because the w in ADFlow currently contains
    ! the state that is perturbed during the matrix-free
    ! operations.

    ! wVec contains the state vector
    call VecGetArrayF90(wVec,wvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! deltaW contains the full update
    call VecGetArrayF90(deltaW,dvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    if(.not. ANK_coupled) then
       ii = 1
       do nn=1, nDom
          do sps=1, nTimeIntervalsSpectral
             call setPointers(nn,1_intType,sps)
             do k=2, kl
                do j=2, jl
                   do i=2, il
                      ! multiply the ratios by 10 to check if the change in a
                      ! variable is greater than 10% of the variable itself.

                      ! check density
                      ratio = abs(wvec_pointer(ii)/(dvec_pointer(ii)+eps))*ANK_physLSTol
                      lambdaL = min(lambdaL, ratio)

                      ! increment by 4 because we want to skip momentum variables
                      ii = ii + 4

                      ! check energy
                      ratio = abs(wvec_pointer(ii)/(dvec_pointer(ii)+eps))*ANK_physLSTol
                      lambdaL = min(lambdaL, ratio)
                      ii = ii + 1
                   end do
                end do
             end do
          end do
       end do
    else
       ii = 1
       do nn=1, nDom
          do sps=1, nTimeIntervalsSpectral
             call setPointers(nn,1_intType,sps)
             do k=2, kl
                do j=2, jl
                   do i=2, il
                      ! multiply the ratios by 10 to check if the change in a
                      ! variable is greater than 10% of the variable itself.

                      ! check density
                      ratio = abs(wvec_pointer(ii)/(dvec_pointer(ii)+eps))*ANK_physLSTol
                      lambdaL = min(lambdaL, ratio)

                      ! increment by 4 because we want to skip momentum variables
                      ii = ii + 4

                      ! check energy
                      ratio = abs(wvec_pointer(ii)/(dvec_pointer(ii)+eps))*ANK_physLSTol
                      lambdaL = min(lambdaL, ratio)
                      ii = ii + 1

                        ! needs to be modified
                      ! if coupled ank is used, nstate = nw and this loop is executed
                      ! if no turbulence variables, this loop will be automatically skipped
                      ! check turbulence variable
                      ratio = (wvec_pointer(ii)/(dvec_pointer(ii)+eps))*ANK_physLSTolTurb
                      ! if the ratio is less than min step, the update is either
                      ! in the positive direction, therefore we do not clip it,
                      ! or the update is very limiting, so we just clip the
                      ! individual update for this cell.
                      if (ratio .lt. ANK_stepFactor*ANK_stepMin) then
                        ! The update was very limiting, so just clip this
                        ! individual update and dont change the overall
                        ! step size. To select the new update, instead of
                        ! clipping to zero, we clip to 1 percent of the original.
                        if (ratio .gt. zero) &
                          dvec_pointer(ii) = wvec_pointer(ii)*ANK_physLSTolTurb

                        ! Either case, set the ratio to one. Positive updates
                        ! do not limit the step, negative updates below minimum
                        ! step were already clipped.
                        ratio = one
                      end if
                      lambdaL = min(lambdaL, ratio)
                      ii = ii + 1

                      ! TODO: Do we need physicality checks for the additional turbulence model variables?
                     !  do l=nt1+1, nt2
                     !     ii = ii + 1
                     !  end do
                      ! do this instead of the above loop for now...
                      ! Will need to modify this if we want physicality check
                      ! for the new turb model variables.
                      ii = ii + (nt2-nt1)
                   end do
                end do
             end do
          end do
       end do
    end if

    ! Restore the pointers to PETSc vectors

    call VecRestoreArrayF90(wVec,wvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecRestoreArrayF90(deltaW,dvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Make sure that we did not get any NaN's in the process
    if (isnan(lambdaL)) lambdaL = zero

    ! Finally, communicate the step size across processes and return
    call mpi_allreduce(lambdaL, lambdaP, 1_intType, adflow_real, &
         mpi_min, ADflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine physicalityCheckANK

  subroutine physicalityCheckANKTurb(lambdaP)

    use constants
    use blockPointers, only : ndom, il, jl, kl
    use flowVarRefState, only : nw, nwf, nt1,nt2
    use inputtimespectral, only : nTimeIntervalsSpectral
    use utils, only : setPointers, EChk, myisnan
    use communication, only : ADflow_comm_world
    implicit none

    ! input variable
    real(kind=realType) , intent(inout) :: lambdaP

    ! local variables
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
    real(kind=realType), pointer :: wvec_pointer(:)
    real(kind=realType), pointer :: dvec_pointer(:)
    real(kind=alwaysRealType) :: lambdaL ! L is for local
    real(kind=realType) :: ratio


    ! Determine the maximum step size that would yield
    ! a maximum change of 10% in density, total energy,
    ! and turbulence variable after a KSP solve.

    ! Initialize the local step size as ANK_stepFactor
    ! because the initial step is likely to be equal to this.
    lambdaL = real(lambdaP)

    ! First we need to read both the update and the state
    ! from PETSc because the w in ADFlow currently contains
    ! the state that is perturbed during the matrix-free
    ! operations.

    ! wVec contains the state vector
    call VecGetArrayF90(wVecTurb,wvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! deltaW contains the full update
    call VecGetArrayF90(deltaWTurb,dvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ii = 1
    do nn=1, nDom
        do sps=1, nTimeIntervalsSpectral
            call setPointers(nn,1_intType,sps)
            do k=2, kl
                do j=2, jl
                    do i=2, il
                        ! multiply the ratios by 10 to check if the change in a
                        ! variable is greater than 10% of the variable itself.

                        ! needs to be modified
                        ! if coupled ank is used, nstate = nw and this loop is executed
                        ! if no turbulence variables, this loop will be automatically skipped
                        ! check turbulence variable
                        ratio = (wvec_pointer(ii)/(dvec_pointer(ii)+eps))* ANK_physLSTolTurb
                        ! if the ratio is less than min step, the update is either
                        ! in the positive direction, therefore we do not clip it,
                        ! or the update is very limiting, so we just clip the
                        ! individual update for this cell.
                        if (ratio .lt. ANK_stepFactor*ANK_stepMin) then
                          ! The update was very limiting, so just clip this
                          ! individual update and dont change the overall
                          ! step size. To select the new update, instead of
                          ! clipping to zero, we clip to 1 percent of the original.
                          if (ratio .gt. zero) &
                            dvec_pointer(ii) = wvec_pointer(ii)*ANK_physLSTolTurb

                          ! Either case, set the ratio to one. Positive updates
                          ! do not limit the step, negative updates below minimum
                          ! step were already clipped.
                          ratio = one
                        end if
                        lambdaL = min(lambdaL, ratio)
                        ii = ii + 1

                        ! TODO: Do we need physicality checks for the additional turbulence model variables?
                        !  do l=nt1+1, nt2
                        !     ii = ii + 1
                        !  end do
                        ! do this instead of the above loop for now...
                        ! Will need to modify this if we want physicality check
                        ! for the new turb model variables.
                        ii = ii + (nt2-nt1)
                    end do
                end do
            end do
        end do
    end do

    ! Restore the pointers to PETSc vectors

    call VecRestoreArrayF90(wVecTurb,wvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecRestoreArrayF90(deltaWTurb,dvec_pointer,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Make sure that we did not get any NaN's in the process
    if (isnan(lambdaL)) lambdaL = zero

    ! Finally, communicate the step size across processes and return
    call mpi_allreduce(lambdaL, lambdaP, 1_intType, adflow_real, &
         mpi_min, ADflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine physicalityCheckANKTurb

  subroutine ANKTurbSolveKSP

    ! This routine solves the turbulence model equation using
    ! a similar approach to the main ank solver.

    use constants
    use blockPointers, only : nDom, flowDoms
    use inputIteration, only : L2conv
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputDiscretization, only : approxSA
    use iteration, only : approxTotalIts, totalR0, totalR, currentLevel
    use utils, only : EChk, setPointers, myisnan
    use turbMod, only : secondOrd
    use solverUtils, only : computeUTau
    use NKSolver, only : getEwTol
    use BCRoutines, only : applyAllBC, applyAllBC_block
    use haloExchange, only : whalo2
    use oversetData, only : oversetPresent
    use flowVarRefState, only : nw, nwf, nt1,nt2 , kPresent, pInfCorr
    use communication
    use blockette, only : blocketteRes
    implicit none

    ! Working Variables
    integer(kind=intType) :: ierr, maxIt, kspIterations, nn, sps, reason, nHist, iter, feval
    integer(kind=intType) :: i,j,k,n
    real(kind=realType) :: atol, val, v2, factK, gm1
    real(kind=alwaysRealType) :: rtol, totalR_dummy, linearRes, norm
    real(kind=alwaysRealType) :: resHist(ANK_maxIter+1)
    real(kind=alwaysRealType) :: unsteadyNorm, unsteadyNorm_old
    real(kind=alwaysRealType) :: linResMonitorTurb, totalRTurb
    logical :: secondOrdSave, correctForK, LSFailed

    ! Calculate the residuals and set rVecTurb before the first iteration
    call blocketteRes(useFlowRes=.False.,useStoreWall=.False.)
    call setRVecANKTurb(rVecTurb)

    do n = 1,ANK_nsubIterTurb

        ! Compute the norm of rVecTurb, which is identical to the
        ! norm of the unsteady residual vector.
        call VecNorm(rVecTurb, NORM_2, totalRTurb, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Determine if we need to form the Preconditioner
        if (mod(ANK_iterTurb, ANK_jacobianLag) == 0) then

            ! Actually form the preconditioner and factorize it.
            if (myid .eq. 0 .and. ANK_turbDebug) &
            write(*,*) "Re-doing turb PC"
            call FormJacobianANKTurb()
            ANK_iterTurb = 0
        end if

        ! Increment the iteration counter
        ANK_iterTurb = ANK_iterTurb + 1

        ! Start with trying to take the full step set by the user.
        lambdaTurb = ANK_StepFactor

        ! Dummy matrix assembly for the matrix-free matrix
        call MatAssemblyBegin(dRdwTurb, MAT_FINAL_ASSEMBLY, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call MatAssemblyEnd(dRdwTurb, MAT_FINAL_ASSEMBLY, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        if (totalR > ANK_secondOrdSwitchTol*totalR0) then
            ! Save if second order turbulence is used, we will only use 1st order during ANK (only matters for the coupled solver)
            approxSA = .True.
            secondOrdSave = secondOrd
            secondOrd =.False.

            ! Determine the relative convergence for the KSP solver
            rtol = ANK_rtol ! Just use the input relative tolerance for approximate fluxes
        else
            ! If the second order fluxes are used, Eisenstat-Walker algorithm to determine relateive
            ! convergence tolerance helps with performance.
            totalR_dummy = totalR
            call getEWTol(totalR_dummy, totalR_old, rtolLast, rtol)

            ! Use the ANK rtol if E-W algorithm is not picking anything lower
            rtol = min(ANK_rtol, rtol)
        end if

        ! Record the total residual and relative convergence for next iteration
        totalR_old = totalR
        rtolLast = rtol

        ! Set all tolerances for linear solve:
        atol = totalR0*L2Conv*0.01_realType

        ! Set the iteration limit to maxIt, determined by which fluxes are used.
        ! This is because ANK step require 0.1 convergence for stability during initial stages.
        ! Due to an outdated preconditioner, the KSP solve might take more iterations.
        ! If this happens, the preconditioner is re-computed and because of this,
        ! ANK iterations usually don't take more than 2 times number of ANK_subSpace size iterations
        call KSPSetTolerances(ANK_KSPTurb, rtol, &
        real(atol), real(ANK_divTol), ank_maxIter, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call KSPSetResidualHistory(ANK_KSPTurb, resHist, ank_maxIter+1, PETSC_TRUE, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Set the BaseVector of the matrix-free matrix:
        call formFunction_mf_turb(ctx, wVecTurb, baseResTurb, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call MatMFFDSetBase(dRdWTurb, wVecTurb, baseResTurb, ierr)
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

        ! Get the number of iterations from the KSP solver
        call KSPGetIterationNumber(ANK_KSPTurb, kspIterations, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call KSPGetConvergedReason(ANK_KSPTurb, reason, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Return previously changed variables back to normal, VERY IMPORTANT
        if (totalR > ANK_secondOrdSwitchTol*totalR0) then
            ! Replace the second order turbulence option
            secondOrd = secondOrdSave
            approxSA = .False.
        end if

        ! Compute the maximum step that will limit the change
        ! in SA variable to some user defined fraction.
        call physicalityCheckANKTurb(lambdaTurb)
        !if (myid .eq. 0) write(*,*)"physicality check lambda: ",lambdaTurb
        !lambdaTurb = max(ANK_stepMin, lambdaTurb)

        ! Take the uodate after the physicality check.
        call VecAXPY(wVecTurb, -lambdaTurb, deltaWTurb, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! Set the updated state variables
        call setWANK(wVecTurb,nt1,nt2)

        ! Compute the unsteady residuals. The actual residuals
        ! also get calculated in the process, and are stored in
        ! dw. Make sure to call setRVec/setRVecANK after this
        ! routine because rVec contains the unsteady residuals,
        ! and we need the steady residuals for the next iteration.
        call computeUnsteadyResANKTurb(lambdaTurb)

        ! Count the number of of residual evaluations outside the KSP solve
        feval = 1_intType

        ! Check if the norm of the rVec is bad:
        call VecNorm(rVecTurb, NORM_2, unsteadyNorm, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        ! initialize this outside the ls
        LSFailed = .False.

        if ((unsteadyNorm > totalRTurb*ANK_unstdyLSTol .or. isnan(unsteadyNorm))) then
            ! The unsteady residual is too high or we have a NAN. Do a
            ! backtracking line search until we get a residual that is lower.

            LSFailed = .True.

            ! Restore the starting (old) w value by adding lamda*deltaW
            call VecAXPY(wVecTurb, lambdaTurb, deltaWTurb, ierr)
            call EChk(ierr, __FILE__, __LINE__)

            ! Set the initial new lambda. This is working off the
            ! potentially already physically limited step.
            lambdaTurb = 0.7_realType * lambdaTurb

            backtrack: do iter=1, 12

                ! Apply the new step
                call VecAXPY(wVecTurb, -lambdaTurb, deltaWTurb, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                ! Set and recompute
                call setWANK(wVecTurb,nt1,nt2)

                ! Compute the unsteady residuals with the current step
                call computeUnsteadyResANKTurb(lambdaTurb)
                feval = feval + 1

                call VecNorm(rVecTurb, NORM_2, unsteadyNorm, ierr)
                call EChk(ierr, __FILE__, __LINE__)

                if (unsteadyNorm > totalRTurb*ANK_unstdyLSTol .or. isnan(unsteadyNorm)) then

                    ! Restore back to the original wVec
                    call VecAXPY(wVecTurb, lambdaTurb, deltaWTurb, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    ! Haven't backed off enough yet....keep going
                    lambdaTurb = lambdaTurb * 0.7_realType
                else
                    ! We have succefssfully reduced the norm
                    LSFailed = .False.
                    exit
                end if
            end do backtrack

            if (LSFailed .or. isnan(unsteadyNorm)) then
                ! the line search wasn't much help.

                if (ANK_CFL > ANK_CFLMin) then
                    ! the cfl number is not already at the lower limit.  We
                    ! can cut the CFL back and try again. Set lambda to zero
                    ! to indicate we never took a step.
                    lambdaTurb = zero
                else
                    ! cfl is as low as it goes, try taking the step
                    ! anyway. We can't do  anything else
                    call VecAXPY(wVecTurb, -lambdaTurb, deltaWTurb, ierr)
                    call EChk(ierr, __FILE__, __LINE__)
                end if

                ! Set the state vec and compute the new residual
                call setWANK(wVecTurb,nt1,nt2)
                call blocketteRes(useFlowRes=.False., &
                                  useStoreWall=.False.)
                feval = feval + 1
            end if
        end if

        call setRvecANKTurb(rVecTurb)

        linResMonitorTurb = resHist(kspIterations+1)/resHist(1)

        if ((linResMonitorTurb .ge. ANK_rtol .and. &
            totalR > ANK_secondOrdSwitchTol*totalR0 .and.&
            linResOldTurb .le. ANK_rtol) &
            !.or. LSFailed) then
!            .or. lambdaTurb .le. ANK_stepMin) then
            .or. (lambdaTurb .eq. zero)) then

            ! We should reform the PC since it took longer than we want,
            ! or we need to adjust the CFL because the last update was bad,
            ! or convergence since the last PC update was good enough and we
            ! would benefit from re-calculating the PC.
            ANK_iterTurb = 0
        end if

        ! update the linear residual for next iteration
        linResOldTurb = linResMonitorTurb

        ! Update step monitor
        ! stepMonitor = lambda

        ! Update the approximate iteration counter. The +1 is for the
        ! residual evaluations.
        ! approxTotalIts = approxTotalIts + feval + kspIterations

                    ! Print some info about the turbulence ksp
        if (myid == 0 .and. ANK_turbDebug) then
          Write(*,*) "LIN RES, ITER, INITRES, REASON, STEP", linResMonitorTurb, kspIterations, &
          reshist(1), reason, lambdaTurb
        end if

    end do

  end subroutine ANKTurbSolveKSP

  subroutine ANKStep(firstCall)

    use constants
    use blockPointers, only : nDom, flowDoms, shockSensor, ib, jb, kb, p, w, gamma
    use inputPhysics, only : equations
    use inputIteration, only : L2conv
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputDiscretization, only : lumpedDiss, approxSA
    use iteration, only : approxTotalIts, totalR0, totalR, stepMonitor, linResMonitor, currentLevel, iterType
    use utils, only : EChk, setPointers, myisnan
    use turbAPI, only : turbSolveSegregated
    use turbMod, only : secondOrd
    use solverUtils, only : computeUTau
    use adjointUtils, only : referenceShockSensor
    use NKSolver, only : setRVec, getEwTol
    use initializeFlow, only : setUniformFlow
    use BCRoutines, only : applyAllBC, applyAllBC_block
    use haloExchange, only : whalo2
    use oversetData, only : oversetPresent
    use flowVarRefState, only : nw, nwf, nt1,nt2 , kPresent, pInfCorr
    use flowUtils, only : computeLamViscosity
    use turbUtils, only : computeEddyViscosity
    use communication
    use blockette, only : blocketteRes
    implicit none

    ! Input Variables
    logical, intent(in) :: firstCall

    ! Working Variables
    integer(kind=intType) :: ierr, maxIt, kspIterations, nn, sps, reason, nHist, iter, feval
    integer(kind=intType) :: i,j,k
    real(kind=realType) :: atol, val, v2, factK, gm1
    real(kind=alwaysRealType) :: rtol, totalR_dummy, linearRes, norm
    real(kind=alwaysRealType) :: resHist(ANK_maxIter+1)
    real(kind=alwaysRealType) :: unsteadyNorm, unsteadyNorm_old
    logical :: secondOrdSave, correctForK, LSFailed

    ! Enter this check if this is the first ANK step OR we are switching to the coupled ANK solver
    if (firstCall .or. &
       ((totalR .le. ANK_coupledSwitchTol * totalR0) .and. (.not. ANK_coupled) &
        .and. (equations .eq. RANSEquations))) then

       ! If this is a first call, we need to change the coupled switch
       ! to the correct value.
       if (firstCall) then

         ! Check if we are above or below the coupled switch tolerance
         if (totalR .le. ANK_coupledSwitchTol * totalR0 .and. equations .eq. RANSEquations) then
           ANK_coupled = .True.
         else
           ANK_coupled = .False.
         end if

       ! This is not a first call, and the only option left is that,
       ! we may be switching from uncoupled to coupled
       else
         ANK_coupled = .True.
       end if

       ! If we are in here, destroy the solver regardless,
       ! and set up with the correct coupling mode.
       call destroyANKSolver()
       call setupANKSolver()

       ! Copy the adflow 'w' into the petsc wVec
       call setwVecANK(wVec,1,nstate)

       ! Evaluate the residual before we start
       call blocketteRes(useUpdateDt=.True.)
       if (ANK_coupled) then
          call setRvec(rVec)
       else
          call setRVecANK(rVec)
       end if

       ! Check if we are using the turb KSP
       if ((.not. ANK_coupled) .and. (.not. ANK_useTurbDADI) .and. equations == RANSEquations) then
          call setwVecANK(wVecTurb,nt1,nt2)
          call setRVecANKTurb(rVecTurb)
       end if

       if (firstCall) then
          ! Start with the selected fraction of the ANK_StepFactor
          lambda = ANK_StepFactor
          lambdaTurb = ANK_stepFactor

          ! Initialize some variables
          totalR_old = totalR ! Record the old residual for the first iteration
          rtolLast = ANK_rtol ! Set the previous relative convergence tolerance for the first iteration
          ANK_CFL = ANK_CFL0 ! only set the initial cfl for the first iteration
          ANK_CFLMinBase = ANK_CFLMin0
          totalR_pcUpdate = totalR ! only update the residual at last PC calculation for the first iteration
          linResOld = zero
          linResOldTurb=zero
          ANK_iter = 0
       end if
    else
       ANK_iter = ANK_iter + 1
    end if

    ! Compute the norm of rVec, which is identical to the
    ! norm of the unsteady residual vector.
    call VecNorm(rVec, NORM_2, unsteadyNorm_old, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Determine if if we need to form the Preconditioner
    if (mod(ANK_iter, ANK_jacobianLag) == 0 .or. totalR/totalR_pcUpdate < ANK_pcUpdateTol) then

       ! First of all, update the minimum cfl wrt the overall convergence
       ANK_CFLMin = min(ANK_CFLLimit, ANK_CFLMinBase*(totalR0/totalR)**ANK_CFLExponent)

       ! Update the CFL number depending on the outcome of the last iteration
       if (lambda < ANK_stepMin * ANK_stepFactor) then

          ! The step was too small, cut back the cfl
          ANK_CFL = max(ANK_CFL*ANK_CFLCutback, ANK_CFLMin)

       else if (totalR < totalR_pcUpdate .and. lambda .ge. ANK_constCFLStep * ANK_stepFactor) then

          ! total residuals have decreased since the last cfl
          ! change, or the step was large enough, we can ramp
          ! the cfl up
          ANK_CFL = max(min(ANK_CFL * ANK_CFLFactor**&
               ((totalR_pcUpdate-totalR)/totalR_pcUpdate), ANK_CFLLimit), ANK_CFLMin)

       else

          ! The step was not small, but it was not large enough
          ! to ramp up the cfl, so we keep it constant. Just
          ! make sure that the cfl does not go below the
          ! minimum value.
          ANK_CFL = max(ANK_CFL, ANK_CFLMin)

       end if

       ! Record the total residuals when the PC is calculated.
       totalR_pcUpdate = totalR

       ! Actually form the preconditioner and factorize it.

       call FormJacobianANK()
       if (totalR .le. ANK_secondOrdSwitchTol*totalR0) then
           if (ANK_coupled) then
               iterType = "  *CSANK"
           else
               iterType = "   *SANK"
           end if
       else
           if (ANK_coupled) then
               iterType = "   *CANK"
           else
               iterType = "    *ANK"
           end if
       end if
       ANK_iter = 0

       ! Also update the turb PC bec. the CFL has changed
       ANK_iterTurb = 0
    else
       if (totalR .le. ANK_secondOrdSwitchTol*totalR0) then
           if (ANK_coupled) then
               iterType = "   CSANK"
           else
               iterType = "    SANK"
           end if
       else
           if (ANK_coupled) then
               iterType = "    CANK"
           else
               iterType = "     ANK"
           end if
       end if
    end if

    ! Start with trying to take the full step set by the user.
    lambda = ANK_StepFactor

    ! Dummy matrix assembly for the matrix-free matrix
    call MatAssemblyBegin(dRdw, MAT_FINAL_ASSEMBLY, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call MatAssemblyEnd(dRdw, MAT_FINAL_ASSEMBLY, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! ============== Flow Update =============

    ! For the approximate solver, we need the approximate flux routines
    ! We set the variables required for approximate fluxes here and they will be used
    ! for the matrix-free matrix-vector product routines when the KSP solver calls it
    ! Very important to set the variables back to their original values after each
    ! KSP solve because we want actual flux functions when calculating residuals
    if (totalR > ANK_secondOrdSwitchTol*totalR0) then
       ! Setting lumped dissipation to true gives approximate fluxes
       ANK_useDissApprox =.True.
       lumpedDiss = .True.
       approxSA = .True.

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

       ! Use the ANK rtol if E-W algorithm is not picking anything lower
       rtol = min(ANK_rtol, rtol)
    end if

    ! Record the total residual and relative convergence for next iteration
    totalR_old = totalR
    rtolLast = rtol

    ! Set all tolerances for linear solve:
    atol = totalR0*L2Conv*0.01_realType

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

    ! Set the BaseVector of the matrix-free matrix:
    call formFunction_mf(ctx, wVec, baseRes, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call MatMFFDSetBase(dRdW, wVec, baseRes, ierr)
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

    ! Get the number of iterations from the KSP solver
    call KSPGetIterationNumber(ANK_KSP, kspIterations, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call KSPGetConvergedReason(ANK_KSP, reason, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Return previously changed variables back to normal, VERY IMPORTANT
    if (totalR > ANK_secondOrdSwitchTol*totalR0) then
       ! Set ANK_useDissApprox back to False to go back to using actual flux routines
       ANK_useDissApprox =.False.
       lumpedDiss = .False.
       approxSA = .False.

       ! Replace the second order turbulence option
       secondOrd = secondOrdSave

    end if

    ! Compute the maximum step that will limit the change in pressure
    ! and energy to some user defined fraction.
    call physicalityCheckANK(lambda)
    if (ANK_CFL .gt. ANK_CFLMin .and. lambda .lt. ANK_stepMin) &
        lambda = zero

    ! Take the uodate after the physicality check.
    call VecAXPY(wVec, -lambda, deltaW, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set the updated state variables
    call setWANK(wVec,1,nState)

    ! Compute the unsteady residuals. The actual residuals
    ! also get calculated in the process, and are stored in
    ! dw. Make sure to call setRVec/setRVecANK after this
    ! routine because rVec contains the unsteady residuals,
    ! and we need the steady residuals for the next iteration.
    call computeUnsteadyResANK(lambda)

    ! Count the number of of residual evaluations outside the KSP solve
    feval = 1_intType

    ! Check if the norm of the rVec is bad:
    call VecNorm(rVec, NORM_2, unsteadyNorm, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! initialize this outside the ls
    LSFailed = .False.

    if ((unsteadyNorm > unsteadyNorm_old*ANK_unstdyLSTol .or. isnan(unsteadyNorm))) then
       ! The unsteady residual is too high or we have a NAN. Do a
       ! backtracking line search until we get a residual that is lower.

       LSFailed = .True.

       ! Restore the starting (old) w value by adding lamda*deltaW
       call VecAXPY(wVec, lambda, deltaW, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set the initial new lambda. This is working off the
       ! potentially already physically limited step.
       lambda = 0.7_realType * lambda

       backtrack: do iter=1, 12

          ! Apply the new step
          call VecAXPY(wVec, -lambda, deltaW, ierr)
          call EChk(ierr, __FILE__, __LINE__)

          ! Set and recompute
          call setWANK(wVec,1,nState)

          ! Compute the unsteady residuals with the current step
          call computeUnsteadyResANK(lambda)
          feval = feval + 1

          call VecNorm(rVec, NORM_2, unsteadyNorm, ierr)
          call EChk(ierr, __FILE__, __LINE__)

          if (unsteadyNorm > unsteadyNorm_old*ANK_unstdyLSTol .or. isnan(unsteadyNorm)) then

             ! Restore back to the original wVec
             call VecAXPY(wVec, lambda, deltaW, ierr)
             call EChk(ierr, __FILE__, __LINE__)

             ! Haven't backed off enough yet....keep going
             lambda = lambda * 0.7_realType
          else
             ! We have succefssfully reduced the norm
             LSFailed = .False.
             exit
          end if
       end do backtrack

       if (LSFailed .or. isnan(unsteadyNorm)) then
          ! the line search wasn't much help.

          if (ANK_CFL > ANK_CFLMin) then
             ! the cfl number is not already at the lower limit.  We
             ! can cut the CFL back and try again. Set lambda to zero
             ! to indicate we never took a step.
             lambda = zero
          else
             ! cfl is as low as it goes, try taking the step
             ! anyway. We can't do  anything else
             call VecAXPY(wVec, -lambda, deltaW, ierr)
             call EChk(ierr, __FILE__, __LINE__)
          end if

          ! Set the state vec and compute the new residual
          call setWANK(wVec,1,nState)
          if (.not. ANK_coupled) then
             call blocketteRes(useTurbRes=.False., useStoreWall=.False.)
          else
             call blocketteRes()
          end if
          feval = feval + 1
       else
       end if
    end if

    ! ============== Turb Update =============
    if ((.not. ANK_coupled) .and. equations==RANSEquations .and. lambda > zero) then

        if (ANK_useTurbDADI) then
            ! actually do the turbulence update
            call computeUtau
            call turbSolveSegregated
        else
            call ANKTurbSolveKSP
        end if
    end if

    ! We need to now compute the residual for the next iteration.  We
    ! also need the to update the update the time step and the
    ! viscWall pointer stuff

    call blocketteRes(useUpdateDt=.True.)

    feval = feval + 1
    if (ANK_coupled) then
       call setRvec(rVec)
    else
       call setRVecANK(rVec)
    end if

    linResMonitor = resHist(kspIterations+1)/resHist(1)

    if ((linResMonitor .ge. ANK_rtol .and. &
         totalR > ANK_secondOrdSwitchTol*totalR0 .and.&
         linResOld .le. ANK_rtol) &
         !.or. LSFailed) then
         !.or. lambda .le. ANK_stepMin) then
         .or. (lambda .eq. zero)) then
       ! We should reform the PC since it took longer than we want,
       ! or we need to adjust the CFL because the last update was bad,
       ! or convergence since the last PC update was good enough and we
       ! would benefit from re-calculating the PC.
       ANK_iter = -1
    end if

    ! update the linear residual for next iteration
    linResOld = linResMonitor

    ! Update step monitor
    stepMonitor = lambda

    ! Check if the linear solutions are failing.
    ! If the lin res is above .5 or so, the solver
    ! might stall, so we might be better off just
    ! reducing the CFL and keep going. We Modify
    ! the CFLMin by altering CFLMinBase.
    if (linResMonitor .gt. ANK_linResMax) then
      ! This will adjust MinBase such that we can halve the cfl
      ! based on the current CFL.
      ANK_CFLMinBase = ANK_CFLCutback*ANK_CFL*((totalR/totalR0)**ANK_CFLExponent)
      ! flags to refresh the Jacobian and cut back the CFL
      ANK_iter = -1
      lambda = zero
    end if

    ! Update the approximate iteration counter. The +1 is for the
    ! residual evaluations.
    approxTotalIts = approxTotalIts + feval + kspIterations

  end subroutine ANKStep
end module ANKSolver
