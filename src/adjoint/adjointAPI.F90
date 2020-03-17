module adjointAPI

contains
#ifndef USE_COMPLEX
  subroutine computeMatrixFreeProductFwd(xvdot, extradot, wdot, bcDataValuesdot, useSpatial, &
       useState, famLists, bcDataNames, bcDataValues, bcDataFamLists, bcVarsEmpty, dwdot, funcsDot, fDot, &
       costSize, fSize, nTime)

    ! This is the main matrix-free forward mode computation
    use constants
    use adjointvars
    use blockPointers, only : nDom
    use communication, only : adflow_comm_world
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputPhysics, only :pointRefd, alphad, betad, equations, machCoefd, &
         machd, machGridd, rgasdimd
    use iteration, only : currentLevel, groundLevel
    use flowVarRefState, only : pInfDimd, rhoInfDimd, TinfDimd
    use adjointUtils, only : allocDerivativeValues, zeroADSeeds
    use masterRoutines, only : master_d
    implicit none

    ! Input Variables
    real(kind=realType), dimension(:), intent(in) :: xvdot
    real(kind=realType), dimension(:), intent(in) :: extradot
    real(kind=realType), dimension(:), intent(in) :: wdot
    logical, intent(in) :: useSpatial, useState
    integer(kind=intType), dimension(:, :) :: famLists
    integer(kind=intType) :: costSize, fSize, nTime

    character, dimension(:, :), intent(in) :: bcDataNames
    real(kind=realType), dimension(:), intent(in) :: bcDataValues, bcDataValuesDot
    integer(kind=intType), dimension(:, :) :: bcDataFamLists
    logical, intent(in) :: BCVarsEmpty

    ! Ouput Variables
    real(kind=realType), dimension(size(wdot)), intent(out) :: dwDot
    real(kind=realType), dimension(costSize, size(famLists,1)), intent(out) :: funcsDot
    real(kind=realType), dimension(3, fSize, nTime), intent(out) :: fDot

    ! Working Variables
    integer(kind=intType) :: nn,sps, level
    real(kind=realType), dimension(costSize, size(famLists,1)) :: funcs

    ! Need to trick the residual evalution to use coupled (mean flow and
    ! turbulent) together.
    level = 1
    currentLevel = level
    groundLevel = level

    ! Allocate the memory we need for derivatives if not done so
    ! already. Note this isn't deallocated until the adflow is
    ! destroyed.
    if (.not. derivVarsAllocated) then
       call allocDerivativeValues(level)
    end if

    ! Zero all AD seesd.
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call zeroADSeeds(nn,level, sps)
       end do
    end do

    ! Set the extra seeds now do the extra ones. Note that we are assuming the
    ! machNumber used for the coefficients follows the Mach number,
    ! not the grid mach number.
    alphad = extraDot(iAlpha)
    betad = extraDot(iBeta)
    machd = extraDot(iMach)
    machCoefd = extraDot(iMach)
    machGridd = extraDot(iMachGrid)
    PinfDimd = extraDot(iPressure)
    rhoinfDimd = extraDot(iDensity)
    tinfdimd = extraDot(iTemperature)
    pointrefd(1) = extraDot(iPointRefX)
    pointrefd(2) = extraDot(iPointRefY)
    pointrefd(3) = extraDot(iPointRefZ)
    rgasdimd = zero

    ! Run the super-dee-duper master forward rotuine
   if (bcVarsEmpty) then
      call master_d(wDot, xVDot, fDot, dwDot, famLists, funcs, funcsDot)
   else
      call master_d(wDot, xVDot, fDot, dwDot, &
           famLists, funcs, funcsDot, bcDataNames, bcDataValues, bcDataValuesdot, bcDataFamLists)
    end if

  end subroutine computeMatrixFreeProductFwd

  subroutine computeMatrixFreeProductBwd(dwbar, funcsBar, fbar, useSpatial, useState, xvbar, &
       extrabar, wbar, spatialSize, extraSize, stateSize, famLists, &
       bcDataNames, bcDataValues, bcDataValuesbar, bcDataFamLists, BCVarsEmpty)
    use constants
    use communication, only : adflow_comm_world
    use blockPointers, only : nDom, dwd, il, jl, kl
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputPhysics, only : equations
    use iteration, only : currentLevel, groundLevel
    use flowVarRefState, only : nw, nwf
    use inputAdjoint, only : frozenTurbulence
    use ADjointPETSc, only : x_like, psi_like3
    use adjointvars, only : derivVarsAllocated
    use utils, only : setPointers_d, EChk
    use adjointUtils, only : allocDerivativeValues, zeroADSeeds
    use masterRoutines, only : master_b
    implicit none

    ! Input Variables
    integer(kind=intType), intent(in) :: stateSize, extraSize, spatialSize
    real(kind=realType), dimension(:), intent(in) :: dwbar
    real(kind=realType), dimension(:, :), intent(in) :: funcsBar
    real(kind=realType), dimension(:, :, :) :: fBar
    logical, intent(in) :: useSpatial, useState
    integer(kind=intType), intent(in) :: famLists(:, :)
    character, dimension(:, :), intent(in) :: bcDataNames
    real(kind=realType), dimension(:), intent(in) :: bcDataValues
    integer(kind=intType), dimension(:, :) :: bcDataFamLists
    logical, intent(in) :: BCVarsEmpty

    ! Ouput Variables
    real(kind=realType), dimension(stateSize), intent(out) :: wbar
    real(kind=realType), dimension(extraSize), intent(out) :: extrabar
    real(kind=realType), dimension(spatialSize), intent(out) :: xvbar
    real(kind=realType), dimension(size(bcDataValues)), intent(out) :: bcDataValuesbar

    ! Working variables
    integer(kind=intType) :: nn, sps, i, j, k, l, ii, level, nState, mm
    logical :: resetToRans
    real(kind=realType), dimension(size(funcsBar,1), size(funcsBar, 2)) :: funcs

    ! Setup number of state variable based on turbulence assumption
    if ( frozenTurbulence ) then
       nState = nwf
    else
       nState = nw
    endif

    ! Need to trick the residual evalution to use coupled (mean flow and
    ! turbulent) together.
    level = 1
    currentLevel = level
    groundLevel = level

    ! Determine if we want to use frozenTurbulent Adjoint
    resetToRANS = .False.
    if (frozenTurbulence .and. equations == RANSEquations) then
       equations = NSEquations
       resetToRANS = .True.
    end if

    ! Allocate the memory for reverse
    if (.not. derivVarsAllocated) then
       call allocDerivativeValues(level)
    end if
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
         call zeroADSeeds(nn,level, sps)
       end do
    end do

    if (bcVarsEmpty) then
       call master_b(wbar, xvbar, extraBar, fBar, dwbar, nState, famLists, &
            funcs, funcsBar)
    else
       call master_b(wbar, xvbar, extraBar, fBar, dwbar, nState, famLists, &
            funcs, funcsBar, bcDataNames, bcDataValues, bcDataValuesbar, bcDataFamLists)
    end if

    ! Reset the correct equation parameters if we were useing the frozen
    ! Turbulent
    if (resetToRANS) then
       equations = RANSEquations
    end if

  end subroutine computeMatrixFreeProductBwd

  subroutine computeMatrixFreeProductBwdFast(dwbar, wbar, stateSize)
    ! This is the "Fast" ie. State variable only version of the reverse
    ! mode computation. It is intended to compute dRdw^T product
    ! ONLY. The main purpose is for fast matrix-vector products for the
    ! actual adjoint solve.
    use constants
    use inputPhysics, only : equations
    use inputAdjoint, only : frozenTurbulence
    use flowVarRefState, only : nw, nwf
    use iteration, only : currentLevel, groundLevel
    use masterRoutines, only : master_state_b, master_b
    use blockpointers
    use inputtimespectral
    use adjointutils
    implicit none

    ! Input Variables
    integer(kind=intType), intent(in) :: stateSize
    real(kind=realType), dimension(stateSize), intent(in) :: dwbar

    ! Ouput Variables
    real(kind=realType), dimension(stateSize), intent(out) :: wbar


    real(kind=realType), dimension(:), allocatable :: extrabar
    real(kind=realType), dimension(:), allocatable :: xvbar
    real(kind=realType), dimension(:, :, :), allocatable :: fBar


    ! Working variables
    integer(kind=intType) :: nState, level, nn, sps
    logical :: resetToRans

    ! Setup number of state variable based on turbulence assumption
    if ( frozenTurbulence ) then
       nState = nwf
    else
       nState = nw
    endif

    ! Assembling matrix on coarser levels is not entirely implemented yet.
    level = 1
    currentLevel = level
    groundLevel = level

    ! Determine if we want to use frozenTurbulent Adjoint
    resetToRANS = .False.
    if (frozenTurbulence .and. equations == RANSEquations) then
       equations = NSEquations
       resetToRANS = .True.
    end if

    ! Note: The calling routine is responsible for ensuring that the
    ! derivative values are allocated AND ZEROED! This routine makes use
    ! of the fact that only wbar needs to be zeroed since all other
    ! required seeds are zeroed in the individual fast routines. This is
    ! slightly unsafe, but it necessary for speed.
     do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call zeroADSeeds(nn,level, sps)
       end do
    end do
    ! allocate(xvbar(1000000), extraBar(100), fBar(3, 1466, 1))
    ! extraBar = zero
    ! xvbar = zero
    ! fbar = zero
    call master_state_b(wBar, dwBar, nState)
    !call master_b(wbar, xvbar, extraBar, fBar, dwBar, nstate)

    ! Reset the correct equation parameters if we are using the frozen
    ! Turbulent
    if (resetToRANS) then
       equations = RANSEquations
    end if
  end subroutine computeMatrixFreeProductBwdFast


#endif
  ! if def for complex

  subroutine solveAdjointForRHS(inVec, outVec, nDOF, relativeTolerance)

    use ADJointPETSc
    use inputADjoint
    use adjointvars
    use killsignals
    use constants
    use blockPointers
    use inputTimeSpectral
    use utils, only : EChk
    use adjointUtils, only : allocDerivativeValues, zeroADSeeds
    implicit none

    ! Input Variables
    real(kind=realType), dimension(ndof), intent(in) :: inVec
    real(kind=realType), dimension(ndof), intent(out) :: outVec
    real(kind=realType), intent(in) :: relativeTolerance
    integer(kind=intType), intent(in) :: nDOF
    integer(kind=intTYpe) :: adjointConvergedReason
    ! Working variables
    integer(kind=intType) :: ierr, nn, sps
    real(kind=realType) :: val

#ifndef USE_COMPLEX

    ! Place the arrays
    call VecPlaceArray(psi_like1, inVec, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecPlaceArray(psi_like2, outVec, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Zero out initial solution
    call VecSet(psi_like2, zero, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Reset normType to default in case the globalPrecon function has been called
    ! by the user.-1 is KSP_NORM_DEFAULT, which isn't in the fortran
    ! header for some reason
    call KSPSetNormType(adjointKSP, -1, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set desired realtive tolerance
    call KSPSetTolerances(adjointKSP, relativeTolerance, adjAbsTol, adjDivTol, &
         adjMaxIter, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Make sure the derivative memory is allocated and zeroed.
    if (.not. derivVarsAllocated) then
       call allocDerivativeValues(1_intType)
    end if

    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call zeroADSeeds(nn, 1_intType, sps)
       end do
    end do

    ! Solve (remember this is actually a transpose solve)
    call KSPSolve(adjointKSP, psi_like1, psi_like2, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call KSPGetConvergedReason(adjointKSP, adjointConvergedReason, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (adjointConvergedReason ==  KSP_CONVERGED_RTOL .or. &
         adjointConvergedReason ==  KSP_CONVERGED_ATOL .or. &
         adjointConvergedReason ==  KSP_CONVERGED_HAPPY_BREAKDOWN) then
       adjointFailed = .False.
    else
       adjointFailed = .True.
    end if

    ! Rest arrays
    call VecResetArray(psi_like1,  ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecResetArray(psi_like2,  ierr)
    call EChk(ierr,__FILE__,__LINE__)

#endif

  end subroutine solveAdjointForRHS

  subroutine solveDirectForRHS(inVec, outVec, nDOF, relativeTolerance)


    use ADJointPETSc
    use inputADjoint
    use adjointVars
    use constants
    use killsignals
    use blockPointers
    use inputTimeSpectral
    use utils, only : EChk
    use adjointUtils, only : allocDerivativeValues, zeroADSeeds
    implicit none

    ! Input Variables
    real(kind=realType), dimension(ndof), intent(in) :: inVec
    real(kind=realType), dimension(ndof), intent(out) :: outVec
    real(kind=realType), intent(in) :: relativeTolerance
    integer(kind=intType), intent(in) :: nDOF
    integer(kind=intTYpe) :: adjointConvergedReason
    ! Working variables
    integer(kind=intType) :: ierr, nn, sps

#ifndef USE_COMPLEX

    ! Place the arrays
    call VecPlaceArray(psi_like1, inVec, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecPlaceArray(psi_like2, outVec, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Zero out initial solution
    call VecSet(psi_like2, zero, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Reset normType to default in case the globalPrecon function has been called
    ! by the user.-1 is KSP_NORM_DEFAULT, which isn't in the fortran
    ! header for some reason
    call KSPSetNormType(adjointKSP, -1, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set desired realtive tolerance
    call KSPSetTolerances(adjointKSP, relativeTolerance, adjAbsTol, adjDivTol, &
         adjMaxIter, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Make sure the derivative memory is allocated and zeroed.
    if (.not. derivVarsAllocated) then
       call allocDerivativeValues(1_intType)
    end if

    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call zeroADSeeds(nn, 1_intType, sps)
       end do
    end do

    ! Solve (this is the transpose solve of a transpose matrix, so it's direct)
    call KSPSolveTranspose(adjointKSP, psi_like1, psi_like2, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call KSPGetConvergedReason(adjointKSP, adjointConvergedReason, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (adjointConvergedReason ==  KSP_CONVERGED_RTOL .or. &
         adjointConvergedReason ==  KSP_CONVERGED_ATOL .or. &
         adjointConvergedReason ==  KSP_CONVERGED_HAPPY_BREAKDOWN) then
       adjointFailed = .False.
    else
       adjointFailed = .True.
    end if

    ! Rest arrays
    call VecResetArray(psi_like1, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecResetArray(psi_like2, ierr)
    call EChk(ierr,__FILE__,__LINE__)

#endif

  end subroutine solveDirectForRHS

  subroutine saveADjointMatrix(fileName)

    use constants
    use ADjointPETSc, only: drdwt
    use communication, only : adflow_comm_world
    use utils, only : EChk
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input params
    character*(*), intent(in) :: fileName

    ! Working parameters
    PetscViewer binViewer
    integer(kind=intType) :: ierr

    call PetscViewerBinaryOpen(adflow_comm_world, fileName, FILE_MODE_WRITE, binViewer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call MatView(dRdwT, binViewer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call PetscViewerDestroy(binViewer,ierr)
    call EChk(ierr, __FILE__, __LINE__)

  end subroutine saveADjointMatrix

  subroutine saveAdjointPC(fileName)

    use constants
    use ADjointPETSc, only: drdwpret
    use communication, only : adflow_comm_world
    use utils, only : EChk
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input params
    character*(*), intent(in) :: fileName

    ! Working parameters
    PetscViewer binViewer
    integer(kind=intType) :: ierr

    call PetscViewerBinaryOpen(adflow_comm_world, fileName, FILE_MODE_WRITE, binViewer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call MatView(dRdwPreT, binViewer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call PetscViewerDestroy(binViewer,ierr)
    call EChk(ierr, __FILE__, __LINE__)

  end subroutine saveAdjointPC

  subroutine saveAdjointRHS(RHS, fileName, nstate)

    use constants
    use ADjointPETSc, only: psi_like1
    use communication, only : adflow_comm_world
    use utils, only : EChk
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input params
    character*(*), intent(in) :: fileName
    real(kind=realType), dimension(nState) :: RHS
    integer(kind=intType) :: nstate

    ! Working parameters
    PetscViewer binViewer
    integer(kind=intType) :: ierr

    ! Dump RHS into psi_like1
    call VecPlaceArray(psi_like1, RHS, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call PetscViewerBinaryOpen(adflow_comm_world, fileName, FILE_MODE_WRITE, binViewer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecView(psi_like1, binViewer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call PetscViewerDestroy(binViewer,ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecResetArray(psi_like1, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine saveAdjointRHS


  subroutine spectralPrecscribedMotion(input, nin, dXv, nout)

    use constants
    use blockPointers, only : il, jl, kl, nDom
    use section, only : sections, nSections
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use monitor , only : timeUnsteadyRestart, timeUnsteady
    use utils, only : setPointers, rotMatrixRigidBody
    implicit none
    ! Input/Output Variables
    integer(kind=intType), intent(in) :: nin, nout
    real(kind=realType), intent(out)  :: dXv(nout)
    real(kind=realType), intent(in)   :: input(nin)

    ! Local Variables
    integer(kind=intType) :: ierr, sps, i, nn, mm, counter0, counter1
    integer(kind=intType) :: nodes_on_block, cum_nodes_on_block
    real(kind=realType), dimension(3)   :: rotationPoint, r
    real(kind=realType), dimension(3, 3) :: rotationMatrix
    real(kind=realType) :: t(nSections), dt(nSections)
    real(kind=realType) :: tOld, tNew, pt(3)
    real(kind=realType), pointer :: xvec_pointer(:)
    real(kind=realType) :: time(3)

    !       For the TimeSpectral case, we need to include    *
    !      the operation that rotates the base grid to each time instance
    !      This is basically the reverse of the operation that is done in
    !      setGrid.f90
    !      The operation in setGrid.f90 is the following
    !      X_sps = M(X - rotPoint) + rotPoint
    !      where
    !      X_sps is the set of coordinates at each time instance
    !      M is the rotation matrix calculated by rotMatrixRigidBody
    !      rotPoint is the point about which the motion takes place
    !      It is easy to see dX_sps/dX = M
    !      What we are actually computing is the following:
    !                 T          T
    !        /dX_sps \ /   dR   \
    !        |-------| |------- |  psi
    !        \  dX   / \ dX_sps /

    ! Zero dXv for time spectral case since we add to array.
    dXv = zero

    do nn=1, nSections
       dt(nn) = sections(nn)%timePeriod &
            / real(nTimeIntervalsSpectral, realType)
    enddo

    timeUnsteady = zero
    counter0 = 0
    cum_nodes_on_block = 0
    ! The nDom loop followed by the sps loop is required to follow
    ! the globalNode ordering such that we can use the pointer from
    ! vecGetArrayF90

    do nn=1, nDom
       do sps = 1, nTimeIntervalsSpectral

          call setPointers(nn, 1, sps)
          nodes_on_block = il*jl*kl

          do mm=1, nSections
             t(mm) = (sps-1)*dt(mm)
          enddo

          ! Compute the displacements due to the rigid motion of the mesh.

          tNew = timeUnsteady + timeUnsteadyRestart
          tOld = tNew - t(1)

          call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

          ! Take rotation Matrix Transpose
          rotationMatrix = transpose(rotationMatrix)

          counter1 = cum_nodes_on_block

          ! Loop over the localally owned nodes:
          do i=1, nodes_on_block
             pt = (/input(3*counter0+1), &
                  input(3*counter0+2), &
                  input(3*counter0+3)/)

             dXv(3*counter1+1:3*counter1+3) = &
                  dXv(3*counter1+1:3*counter1+3) + &
                  matmul(rotationMatrix, pt)

             counter0 = counter0 + 1
             counter1 = counter1 + 1
          end do

       end do
       ! Increment the cumulative number of nodes by the nodes on the
       ! block we just did
       cum_nodes_on_block = cum_nodes_on_block + nodes_on_block
    end do

  end subroutine spectralPrecscribedMotion

  subroutine setupAllResidualMatricesfwd

    use constants
    use ADjointPETSc, only : dRdwT
    use communication, only : adflow_comm_world, myid
    use inputADjoint, only : frozenTurbulence, useMatrixFreedRdw
    use adjointUtils, only : setupStateResidualMatrix
    use utils, only : EChk
    implicit none

    logical :: useAD, useTranspose, usePC, useObjective
    real(kind=realType) :: timeAdjLocal, timeAdj, time(2)
    integer(kind=intType) :: ierr

    ! If we are assembling matrices...we ned to assemble the
    ! 'transpose', with 'AD', we want the exact matrix not the 'PC',
    ! and will compute objective RHS
    useAD = .True.
    usePC = .False.
    useTranspose = .True.
    useObjective = .True.

    if (.not. useMatrixFreedRdw) then
       if( myid ==0 ) then
          write(*, 10) "Assembling State Residual Matrix in Forward mode..."
       end if
       time(1) = mpi_wtime()
       call setupStateResidualMatrix(drdwT, useAD, usePC, useTranspose, &
            useObjective, frozenTurbulence, 1_intType)
       time(2) = mpi_wtime()
       timeAdjLocal = time(2)-time(1)

       call mpi_reduce(timeAdjLocal, timeAdj, 1, adflow_real, &
            mpi_max, 0, ADFLOW_COMM_WORLD, ierr)
       call EChk(ierr,  __FILE__, __LINE__)

       if(myid ==0)  then
          write(*, 20) "Assembling State Residaul Matrices Fwd time (s) = ", timeAdj
       end if
    end if

    ! Output formats.
10  format(a)
20  format(a, 1x, f8.2)

  end subroutine setupAllResidualMatricesfwd

  subroutine solveAdjoint(RHS, psi, checkSolution, nState)
    !
    !      Solve the linear discrete ADjoint system of equations
    !          [dR/dW]T . psi = {RHS}
    !      using preconditioned GMRES provided by PETSc. The values in psi
    !      are significant as they are used as the inital guess.
    !

    use constants, only : realType, intType, alwaysRealType, one, adflow_real, &
         mpi_max, mpi_sum, mpi_double_precision, mpi_integer, mpi_double_complex
    use ADjointPETSc, only : dRdwT, psi_like1, psi_like2, adjointKSP, &
         adjResInit, adjResStart, adjResFinal

    use killsignals, only : adjointFailed
    use inputADjoint, only : adjAbsTol, adjDivTol, adjMaxIter, adjRelTol, &
         adjRelTolRel, printTiming
    use adjointVars, only: derivVarsAllocated
    use communication, only : myid, adflow_comm_world
    use blockPointers, only : nDom
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use adjointUtils, only : allocDerivativeValues, zeroADSeeds
    use utils, only : EChk
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input Parameters
    real(kind=realType), dimension(nState) :: RHS, psi
    integer(kind=intType) :: nState
    logical :: checkSolution
    !
    !     Local variables.
    real(kind=alwaysRealType)   :: norm
    real(kind=alwaysRealType), dimension(2) :: time
    real(kind=alwaysRealType)               :: timeAdjLocal, timeAdj
    real(kind=alwaysRealType) :: l2abs, l2rel
    integer(kind=intType) :: ierr, nn, sps
    integer(kind=intType) :: adjConvIts
    KSPConvergedReason adjointConvergedReason
    Vec adjointRes, RHSVec

#ifndef USE_COMPLEX
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

    ! Reset normType to default in case the globalPrecon function has been called
    ! by the user. -1 is KSP_NORM_DEFAULT, which isn't in the fortran
    ! header for some reason
    call KSPSetNormType(adjointKSP, -1, ierr)
    call EChk(ierr, __FILE__, __LINE__)

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

       ! call mpi_reduce(timeAdjLocal, timeAdj, 1, adflow_real, &
       !      mpi_max, 0, ADFLOW_COMM_WORLD, ierr)

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

10  format(a)
20  format(a,1x,f8.2)
30  format(1x,a,1x,e10.4,4x,a,1x,i4)
40  format(1x,a,1x,i5,1x,a)

#endif

  end subroutine solveAdjoint

  subroutine setupPETScKsp

    use ADjointPETSc, only: drdwpret, drdwt, adjointKSP
    use inputADjoint
    use utils, only : ECHk, terminate
    use adjointUtils, only : mykspmonitor
    use adjointUtils, only : setupStateResidualMatrix, setupStandardKSP, setupStandardMultigrid
    use communication
    use agmg, only : setupShellPC, destroyShellPC, applyShellPC
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    !     Local variables.
    logical :: useAD, usePC, useTranspose, useObjective, useCoarseMats
    integer(kind=intType) :: ierr
    real(kind=realType) :: timeA
    PC shellPC
    if (ApproxPC)then
       !setup the approximate PC Matrix
       useAD = ADPC
       useTranspose = .True.
       usePC = .True.
       useObjective = .False.
       useCoarseMats = .False.
       if (preCondType == 'mg') then
          useCoarseMats = .True.
       end if

       call setupStateResidualMatrix(drdwpret, useAD, usePC, useTranspose, &
            useObjective, frozenTurbulence, 1_intType, useCoarseMats=useCoarseMats)

       call KSPSetOperators(adjointKSP, dRdwT, dRdWPreT, ierr)
       call EChk(ierr, __FILE__, __LINE__)

    else
       ! Use the exact jacobian.  Here the matrix that defines the
       ! linear system also serves as the preconditioning matrix. This
       ! is only valid if useMatrixFree is flase.
       if (useMatrixfreedRdw) then
          call terminate("setupPETScKSP", "useMatrixFreedRdW option cannot be true when the approxPC option is False")
       end if
       call KSPSetOperators(adjointKSP, dRdWt, dRdWT, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end if

    if (PreCondType == 'asm') then
       ! Run the super-dee-duper function to setup the ksp object:

       call setupStandardKSP(adjointKSP, ADjointSolverType, adjRestart, adjointpcside, &
            PreCondType, overlap, outerPreConIts, localPCType, &
            matrixOrdering, FillLevel, innerPreConIts)
    else if (PreCondType == 'mg') then

       call setupStandardMultigrid(adjointKSP, ADjointSolverType, adjRestart, &
            adjointPCSide, overlap, outerPreconIts, matrixOrdering,  fillLevel)
    end if

    ! Setup monitor if necessary:
    if (setMonitor) then
       ! PETSC_NULL_CONTEXT doesn't exit...
       call KSPMonitorSet(adjointKSP, MyKSPMonitor, PETSC_NULL_FUNCTION, &
            PETSC_NULL_FUNCTION, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    endif

  end subroutine setupPETScKsp

  subroutine saveCellCenters(fileName)

    use blockPointers
    use iteration
    use inputTimeSpectral
    use adjointVars, only: nCellsLocal
    use communication
    use utils, only : setPointers, EChk
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input params
    character*(*), intent(in) :: fileName

    ! Working parameters
    PetscViewer binViewer
    integer(kind=intType) :: nn,sps,i,j,k,n
    integer(kind=intType) :: ierr, iRow, level
    real(kind=realType),dimension(3)::cellCenter

    Vec cellCenters

    call VecCreateMPI(ADFLOW_COMM_WORLD, nCellsLocal(1)*3, &
         PETSC_DETERMINE, cellCenters, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call VecSetBlockSize(cellCenters, 3, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! compute and store the cell centers
    level = 1
    domainLoop: do nn=1, nDom
       spectralLoop: do sps=1, nTimeIntervalsSpectral
          call setPointers(nn, level, sps)
          do k=2,kl
             do j=2,jl
                do i=2,il
                   iRow = flowDoms(nn, level, sps)%globalCell(i, j, k)
                   ! The location of the cell center is determined
                   ! by averaging the cell coordinates.
                   do n=1,3
                      cellCenter(n) = (x(i-1,j-1,k-1,n) + x(i,j-1,k-1,n)  &
                           +  x(i-1,j,  k-1,n) + x(i,j,  k-1,n)  &
                           +  x(i-1,j-1,k,  n) + x(i,j-1,k,  n)  &
                           +  x(i-1,j,  k,  n) + x(i,j,  k,  n))/8
                   end do
                   call VecSetValues(cellCenters, 1, [iRow], [cellCenter], INSERT_VALUES, ierr)
                   call EChk(ierr, __FILE__, __LINE__)
                end do
             end do
          end do
       end do spectralLoop
    end do domainLoop

    ! PETSc Matrix Assembly begin
    call VecAssemblyBegin(cellCenters, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecAssemblyEnd  (cellCenters, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call PetscViewerBinaryOpen(adflow_comm_world, fileName, FILE_MODE_WRITE, binViewer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecView(cellCenters, binViewer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call PetscViewerDestroy(binViewer,ierr)
    call EChk(ierr, __FILE__, __LINE__)

  end subroutine saveCellCenters

 subroutine dRdwTMatMult(A, vecX,  vecY, ierr)

    ! PETSc user-defied call back function for computing the product of
    ! dRdwT with a vector. Here we just call the much more broadly
    ! useful routine computeMatrixFreeProductBwdFast()

    use constants
    use communication
    use blockPointers
    use iteration
    use flowVarRefState
    use inputAdjoint
    use ADjointVars
    use inputTimeSpectral
    use utils, only : EChk
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none


    ! PETSc Arguments
    Mat   A
    Vec   vecX, vecY
    integer(kind=intType) ::ierr

    real(kind=realType), pointer :: dwb_pointer(:)
    real(kind=realType), pointer :: wb_pointer(:)

#ifndef USE_COMPLEX
    call VecGetArrayReadF90(vecX, dwb_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecGetArrayF90(VecY, wb_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call computeMatrixFreeProductBwdFast(dwb_pointer, wb_pointer, size(wb_pointer))

    call VecRestoreArrayF90(vecX, dwb_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ierr = 0
#endif

  end subroutine dRdwTMatMult

  subroutine dRdwMatMult(A, vecX,  vecY, ierr)

    ! PETSc user-defied call back function for computing the product of
    ! dRdw with a vector. Here we just call the much more broadly
    ! useful routine computeMatrixFreeProductFwd()

    use constants
    use communication
    use blockPointers
    use iteration
    use flowVarRefState
    use inputAdjoint
    use ADjointVars
    use inputTimeSpectral
    use surfaceFamilies, only : fullFamList, BCFamGroups
    use utils, only : EChk
    use surfaceUtils, only : getSurfaceSize
    use adjointUtils, only : allocDerivativeValues, zeroADSeeds
#ifndef USE_COMPLEX
    use masterRoutines, only : master_d
#endif
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! PETSc Arguments
    Mat   A
    Vec   vecX, vecY
    integer(kind=intType) ::ierr

    real(kind=realType), pointer :: wd_pointer(:)
    real(kind=realType), pointer :: dwd_pointer(:)
    integer(kind=intType) :: stateSize, costSize, fSize, fSIzeCell, spatialSize
    real(kind=realType), dimension(:), allocatable :: Xvdot
    real(kind=realType), dimension(:, :, :), allocatable :: fDot
    integer(kind=intType) :: nn, sps
    integer(kind=intType), dimension(:), pointer :: walLFamList

#ifndef USE_COMPLEX

    call VecGetArrayReadF90(vecX, wd_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecGetArrayF90(VecY, dwd_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    if (.not. derivVarsAllocated) then
       call allocDerivativeValues(1)
    end if

    ! Zero all AD seesd.
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call zeroADSeeds(nn, 1, sps)
       end do
    end do

    wallFamList => BCFamGroups(iBCGroupWalls)%famList
    call getSurfaceSize(fSize, fSizeCell, wallFamList, size(wallFamList), .True.)
    spatialSize =  3 * nNodesLocal(1_intType)*nTimeIntervalsSpectral

    allocate(xvdot(spatialSize))
    allocate(fdot(3, fSize, nTimeIntervalsSpectral))

    xvdot = zero
    fdot = zero

    call master_d(wd_pointer, xvDot, fDot, dwd_pointer)

    deallocate(xvDot, Fdot)

    call VecRestoreArrayReadF90(vecX, wd_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecRestoreArrayF90(VecY, dwd_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ierr = 0
#endif
  end subroutine dRdwMatMult

  subroutine createPETScVars
    !
    !      Create the matrices/vectors that are required for the adjoint
    !
    use constants
    use ADjointPETSc, only: dRdwT, dRdwPreT, &
         adjointKSP, matfreectx, x_like, psi_like1, adjointPETScVarsAllocated
    use ADjointVars
    use communication, only : adflow_comm_world, myid
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowVarRefState, only : nwf, nw, viscous
    use inputADjoint, only : approxPC, frozenTurbulence, useMatrixFreedRdw, viscPC
    use stencils, only : N_visc_drdw, n_euler_drdw, visc_drdw_stencil,  euler_drdw_stencil, &
         visc_drdw_stencil, visc_pc_stencil, N_visc_PC, N_euler_PC, euler_PC_stencil
    use utils, only : EChk, setPointers
    use adjointUtils, only : myMatCreate, destroyPETScVars, statePreAllocation
    use agmg, only : setupAGMG
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    !     Local variables.
    integer(kind=intType)  :: nDimW, nDimX
    integer(kind=intType) :: i, n_stencil, nState
    integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
    integer(kind=intType), dimension(:), allocatable :: nnzDiagonal2, nnzOffDiag2
    integer(kind=intType), dimension(:, :), pointer :: stencil
    integer(kind=intType) :: level, ierr, nlevels
    integer(kind=intType) :: rows(4), iCol, nn, sps, ii
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, iDim, iStride, j, mm
    integer(kind=intType) :: npts, ncells, nTS

    ! Destroy variables if they already exist
    call destroyPETScVars()

    ! DETERMINE ALL SIZES HERE!
    if ( frozenTurbulence ) then
       nState = nwf
    else
       nState = nw
    endif

    nDimW = nState * nCellsLocal(1_intType)*nTimeIntervalsSpectral
    nDimX = 3 * nNodesLocal(1_intType)*nTimeIntervalsSpectral


    if (.not. useMatrixFreedRdw) then
       ! Setup matrix-based dRdwT
       allocate(nnzDiagonal(nCellsLocal(1_intType)*nTimeIntervalsSpectral), &
            nnzOffDiag(nCellsLocal(1_intType)*nTimeIntervalsSpectral) )

       if (viscous) then
          n_stencil = N_visc_drdw
          stencil => visc_drdw_stencil
       else
          n_stencil = N_euler_drdw
          stencil => euler_drdw_stencil
       end if

       level = 1

       call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nState, stencil, n_stencil, &
            level, .true.)
       call myMatCreate(dRdwT, nState, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
            __FILE__, __LINE__)

       call matSetOption(dRdwT, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       deallocate(nnzDiagonal, nnzOffDiag)
    else
       ! Setup matrix-free dRdwT
       call MatCreateShell(ADFLOW_COMM_WORLD, nDimW, nDimW, PETSC_DETERMINE, &
            PETSC_DETERMINE, matfreectx, dRdwT, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set the shell operation for doing matrix vector multiplies
       call MatShellSetOperation(dRdwT, MATOP_MULT, dRdwTMatMult, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set the shell operation for doing TRNASPOSE matrix vector
       ! multiplies
       call MatShellSetOperation(dRdwT, MATOP_MULT_TRANSPOSE, dRdwMatMult, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call MatSetup(dRdwT, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end if

    ! Create the approxPC if required
    if (ApproxPC) then
       ! ------------------- Determine Preallocation for dRdwPre -------------
       allocate(nnzDiagonal(nCellsLocal(1_intType)*nTimeIntervalsSpectral), &
            nnzOffDiag(nCellsLocal(1_intType)*nTimeIntervalsSpectral) )

       if (viscous .and. viscPC) then
          stencil => visc_pc_stencil
          n_stencil = N_visc_pc
       else
          stencil => euler_pc_stencil
          n_stencil = N_euler_pc
       end if

       level = 1
       call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nState, stencil, n_stencil, &
            level, .true.)
       call myMatCreate(dRdwPreT, nState, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
            __FILE__, __LINE__)

       call matSetOption(dRdwPreT, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       deallocate(nnzDiagonal, nnzOffDiag)
    end if


    call setupAGMG(drdwpret, nDimW/nState, nState)

    ! Create the KSP Object
    call KSPCreate(ADFLOW_COMM_WORLD, adjointKSP, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    adjointPETScVarsAllocated = .True.
  end subroutine createPETScVars

end module adjointAPI
