subroutine createPETScVars

#ifndef USE_NO_PETSC
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create the matrices/vectors that are required for the adjoint  *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc, only: dRdwT, dRdwPreT, dJdw, psi, adjointRHS, adjointRes, &
       FMw, dFcdw, dFcdx, dFndFc, dFdx, dFdw, dRdx, xVec, dJdx, FMx, dRda, &
       adjointKSP, dFMdExtra, dRda_data, overArea, fCell, fNode, doAdx, nFM, &
       dRdwTShell, matfreectx
  use ADjointVars   
  use BCTypes
  use communication  
  use inputTimeSpectral 
  use flowVarRefState 
  use inputADjoint    
  use stencils
  use blockPointers
  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
#include "include/petscversion.h"

  !     Local variables.
  integer(kind=intType)  :: nDimW, nDimX, nDimPt, nDimCell
  integer(kind=intType) :: i, n_stencil, nState
  integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  integer(kind=intType), dimension(:), allocatable :: nnzDiagonal2, nnzOffDiag2
  integer(kind=intType), dimension(:, :), pointer :: stencil
  integer(kind=intType) :: level, ierr, nlevels
  integer(kind=intType) :: rows(4), iCol, nn, sps, ii
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, iDim, iStride, j, mm
  integer(kind=intType) :: npts, ncells, nTS
  external dRdwTMatMult
  ! DETERMINE ALL SIZES HERE!
  if ( frozenTurbulence ) then
     nState = nwf
  else
     nState = nw
  endif

  nDimW = nState * nCellsLocal(1_intType)*nTimeIntervalsSpectral
  nDimX = 3 * nNodesLocal(1_intType)*nTimeIntervalsSpectral

  call getForceSize(npts, ncells)
  nDimPt = npts * 3 * nTimeIntervalsSpectral
  nDimCell = nCells * 3 * nTimeIntervalsSpectral

  ! ------------------- Determine Preallocation for dRdw --------------

  allocate(nnzDiagonal(nCellsLocal(1_intType)*nTimeIntervalsSpectral), &
       nnzOffDiag(nCellsLocal(1_intType)*nTimeIntervalsSpectral) )

  call initialize_stencils
  if (.not. viscous) then
     n_stencil = N_euler_drdw
     stencil => euler_drdw_stencil 
  else
     n_stencil = N_visc_drdw
     stencil => visc_drdw_stencil
  end if

  level = 1

  call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nState, stencil, n_stencil, &
       level)
  call myMatCreate(dRdwT, nState, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
       __FILE__, __LINE__)

  call matSetOption(dRdwT, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  deallocate(nnzDiagonal, nnzOffDiag)

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
          level)
     call myMatCreate(dRdwPreT, nState, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
          __FILE__, __LINE__)

     call matSetOption(dRdwPreT, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     deallocate(nnzDiagonal, nnzOffDiag)
     ! --------------------------------------------------------------------
  end if ! Approx PC

  ! Vectors:
  ! Get dJdw and psi from one MatGetVecs Call
  call MatGetVecs(dRdwT, dJdW, psi, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! adjointRes is the same size as dJdw, psi
  call VecDuplicate(dJdW, adjointRes, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDuplicate(dJdW, adjointRHS, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! ------------------- Setup matrix free matix  --------------

  call VecGetSize(psi, mm,ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call MatCreateShell(SUMB_COMM_WORLD, nDimW, nDimW, mm,mm, &
       matfreectx, dRdwTShell, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! call matCreate(SUMB_COMM_WORLD, dRdwTShell, ierr)
  ! call EChk(ierr, __FILE__, __LINE__)
  
  ! call matSetSizes(dRdwTShell, nDimW, nDimW, PETSC_DETERMINE, PETSC_DETERMINE, ierr)
  ! call EChk(ierr, __FILE__, __LINE__)

  ! call matSetType(dRdwTShell, MATSHELL, ierr)
  ! call EChk(ierr, __FILE__, __LINE__)

  call MatShellSetOperation(dRdwTShell, MATOP_MULT, dRdwTMatMult, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call MatSetup(dRdwTShell, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Create the nFM * nTimeIntervalsSpectral vectors for d{F,M}/dw plus
  ! any additional functions 
  allocate(FMw(nFM, nTimeIntervalsSpectral))
  do sps=1,nTimeIntervalsSpectral
     do i=1,nFM
        call VecDuplicate(dJdw, FMw(i, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end do

  ! Create dFcdw, dFcdx, dFcdx2, dFcdFn

  ! dFcdw: Derivative of the face centered forces with respect
  ! states. For Euler each face depends on the two cells above it so
  ! have 2*nw non-zeros per row. For viscous, stecel is 3x3x2 where
  ! two is in the off-wall direction. The larger stencil is due to the
  ! evaluation of viscous fluxes

  ! dFcdx: Derivative of the face centered forces with respect to the
  ! nodes lying on the face. Euler has a stencil of 3x3 (face
  ! only). For viscous forces, the stencil is also 3x3.

  ! dFndFc: Derivative of the nodal forces with respect to the face
  ! centered forces. This matrix will consist entirely of rows with 4
  ! values each of which is exactly 1/4 

  ! dFcdw
  allocate( nnzDiagonal(nDimCell), nnzOffDiag(nDimCell))
  if (.not. viscous) then
     nnzDiagonal = 2 * nState

     ! OffDiag is an overestimate...at most one cell is from neighbour
     ! proc
     nnzOffDiag  = 1 * nState 
  else
     nnzDiagonal = 3*3*2*nState

     ! In general should not have much more than 6 cells on off
     ! proc. If there is a malloc or two that isn't the end of the world
     nnzOffDiag = 5*2*nState
  end if

  call myMatCreate(dFcdw, 1, nDimCell, nDimw, nnzDiagonal, nnzOffDiag, &
       __FILE__, __LINE__)

  ! dFcdx
  if (.not. viscous) then
     nnzDiagonal = 4 * 3
     nnzOffDiag  = 1
  else
     nnzDiagonal = 4*4*2*3
     nnzOffDiag = 7*2*3
  end if

  call myMatCreate(dFcdx, 1, nDimCell, nDimX, nnzDiagonal, nnzOffDiag, &
       __FILE__, __LINE__)
  deallocate(nnzDiagonal, nnzOffDiag)

  ! doAdx -> Derviative of 1/area wrt the spatial nodes.
  allocate( nnzDiagonal(nDimPt), nnzOffDiag(nDimPt))
  nnzDiagonal = 3*3*3
  nnzOffDiag = 0

  call myMatCreate(doAdx, 1, nDimPt, nDimx, nnzDiagonal, nnzOffDiag, &
       __FILE__, __LINE__)

  deallocate(nnzDiagonal, nnzOffDiag)

  ! Finally we need dFndFc
  allocate( nnzDiagonal(nDimPt), nnzOffDiag(nDimPt))

  nnzDiagonal = 4
  nnzOffDiag  = 4 ! This should be enough...might get a couple of mallocs
  call myMatCreate(dFndFc, 1, nDimPt, nDimCell, nnzDiagonal, nnzOffDiag, &
       __FILE__, __LINE__)

  deallocate(nnzDiagonal, nnzOffDiag)

  ! Get a right hand and left hand vec. We need both:
  call MatGetVecs(dFndFc, fCell, fNode, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecDuplicate(fNode, overArea, ierr)
  call EChk(ierr, __FILE__, __LINE__)


  ! We will also take this opportunity to assemble dFndFc. 
  spectral: do sps=1,nTimeIntervalsSpectral
     domains: do nn=1,nDom
        call setPointers(nn,1_intType,sps)

        ! Loop over the number of boundary subfaces of this block.
        bocos: do mm=1,nBocos
           if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) then

              jBeg = BCData(mm)%jnBeg + 1; jEnd = BCData(mm)%jnEnd
              iBeg = BCData(mm)%inBeg + 1; iEnd = BCData(mm)%inEnd

              do j=jBeg, jEnd ! Face Loop
                 do i=iBeg, iEnd ! Face Loop
                    do iDim = 0,2
                       iCol = bcData(mm)%FMCellIndex(i,j)*3 + iDim 
                       rows(1) = bcData(mm)%FMNodeIndex(i-1, j-1)*3 + iDim 
                       rows(2) = bcData(mm)%FMNodeIndex(i  , j-1)*3 + iDim 
                       rows(3) = bcData(mm)%FMNodeIndex(i-1, j  )*3 + iDim 
                       rows(4) = bcData(mm)%FMNodeIndex(i  , j  )*3 + iDim 

                       do ii=1,4
                          call MatSetValues(dFndFc, 1, rows(ii), 1, iCol, &
                               fourth, INSERT_VALUES, ierr) 
                          call EChk(ierr, __FILE__, __LINE__)
                       end do

                    end do
                 end do
              end do
           end if
        end do bocos
     end do domains
  end do spectral

  call MatAssemblyBegin(dFndFc, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call MatAssemblyEnd  (dFndFc, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr,  __FILE__, __LINE__)

  ! For now, leave dFdw, and dFdx in.
  allocate( nnzDiagonal(nDimPt), nnzOffDiag(nDimPt) )
  nnzDiagonal = 8*nState
  nnzOffDiag  = 8*nState! Make the off diagonal the same

  call myMatCreate(dFdw, 1, nDimPt, nDimW, nnzDiagonal, nnzOffDiag, &
       __FILE__, __LINE__)

  ! Create the matrix dFdx
  nnzDiagonal = 27
  nnzOffDiag = 27

  call myMatCreate(dFdx, 1, nDimPt, nDimPt, nnzDiagonal, nnzOffDiag, &
       __FILE__, __LINE__)


  deallocate(nnzDiagonal, nnzOffDiag)

  !     ******************************************************************
  !     *                                                                *
  !     * Create matrix dRdx that is used to compute the total cost /    *
  !     * constraint function sensitivity with respect to the spatial    *
  !     * design variables 'x' as dIdx = dJdx - psi^T dRdx.              *
  !     *                                                                *
  !     * Matrix dRdx has size [nDimW, nDimX] and is generally            *
  !     * sparse for the coordinate design variables.                    *
  !     *                                                                *
  !     * The local dimensions are specified so that the spatial         *
  !     * coordinates x (a) are placed in the local processor. This has  *
  !     * to be consistent with the vectors dIdx and dJdx.               *
  !     *                                                                *
  !     ******************************************************************

  allocate( nnzDiagonal(nDimX), nnzOffDiag(nDimX) )
  ! Create the matrix dRdx.
  level = 1_intType
  call drdxPreAllocation(nnzDiagonal, nnzOffDiag, nDimX, level)
  ! Sanity check on diagonal portion: For 2D cases nnzDiagon may be
  ! too large because of the x-block corners:
  do i=1, nDimx
     nnzDiagonal(i) = min(nnzDiagonal(i), nDimw)
  end do

  ! Note we are creating the TRANPOSE of dRdx. It is size dDimX by nDimW
  call myMatCreate(dRdx, 1, nDimX, nDimW, nnzDiagonal, nnzOffDiag, &
       __FILE__, __LINE__)
  deallocate( nnzDiagonal, nnzOffDiag )

  call MatSetOption(dRdx, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call getForceSize(npts, ncells)

  ! xVec
  call VecCreate(SUMB_COMM_WORLD, xVec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecSetSizes(xVec, nDimX, PETSC_DECIDE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecSetBlockSize(xVec, 3, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecSetType(xVec, "mpi", ierr) 
  call EChk(ierr, __FILE__, __LINE__)

  ! Vectors
  call VecCreate(SUMB_COMM_WORLD, dJdx, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecSetSizes(dJdx, nDimX, PETSC_DECIDE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecSetBlockSize(dJdx, 3, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecSetType(dJdx, "mpi", ierr) 
  call EChk(ierr, __FILE__, __LINE__)

  allocate(FMx(nFM, nTimeIntervalsSpectral))
  do sps=1,nTimeIntervalsSpectral
     do i=1,nFM
        call VecDuplicate(dJdx, FMx(i, sps), ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end do
  end do

  ! Create the KSP Object
  call KSPCreate(SUMB_COMM_WORLD, adjointKSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Once again, PETSC is royally screwed up. You CANNOT use PETSC_NULL
  ! arguments. They simply do NOT work in Fortran. The PETSc
  ! documentation lies to you. We have to allocate our own data. 
  if (allocated(dRda_data)) then
     deallocate(dRda_data)
  end if
  allocate(dRda_data(nDimw, nDesignExtra))

  if (allocated(dFMdExtra)) then
     deallocate(dFMdExtra)
  end if
  allocate(dFMdExtra(nFM, nDesignExtra, nTimeIntervalsSpectral))

#if PETSC_VERSION_MINOR < 3 
     call MatCreateMPIDense(SUMB_COMM_WORLD, nDimW, PETSC_DECIDE, &
          PETSC_DETERMINE, nDesignExtra, dRda_data, dRda, ierr)
#else
     call MatCreateDense(SUMB_COMM_WORLD, nDimW, PETSC_DECIDE, &
          PETSC_DETERMINE, nDesignExtra, dRda_data, dRda, ierr)
#endif
  call EChk(ierr, __FILE__, __LINE__)
#endif
end subroutine createPETScVars

subroutine myMatCreate(matrix, blockSize, m, n, nnzDiagonal, nnzOffDiag, &
     file, line)
  ! Function to create petsc matrix to make stuff a little cleaner in
  ! the code above. Also, PETSc always thinks is a good idea to
  ! RANDOMLY change syntax between versions so this way there is only
  ! one place to make a change based on petsc version. 

  use communication
  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  Mat matrix
  integer(kind=intType), intent(in) :: blockSize, m, n
  integer(kind=intType), intent(in), dimension(*) :: nnzDiagonal, nnzOffDiag
  character*(*) :: file
  integer(kind=intType) :: ierr, line
  if (blockSize > 1) then
#if PETSC_VERSION_MINOR <  3
     call MatCreateMPIBAIJ(SUMB_COMM_WORLD, blockSize, &
          m, n, PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal, 0, nnzOffDiag, matrix, ierr)
#else
     call MatCreateBAIJ(SUMB_COMM_WORLD, blockSize, &
          m, n, PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal, 0, nnzOffDiag, matrix, ierr)
#endif
  else
     
#if PETSC_VERSION_MINOR <  3
     call MatCreateMPIAIJ(SUMB_COMM_WORLD, &
          m, n, PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal, 0, nnzOffDiag, matrix, ierr)
#else
     call MatCreateAIJ(SUMB_COMM_WORLD,&
          m, n, PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal, 0, nnzOffDiag, matrix, ierr)
#endif
     call EChk(ierr, file, line)
  end if
  
  ! Warning: The array values is logically two-dimensional, 
  ! containing the values that are to be inserted. By default the
  ! values are given in row major order, which is the opposite of
  ! the Fortran convention, meaning that the value to be put in row
  ! idxm[i] and column idxn[j] is located in values[i*n+j]. To allow
  ! the insertion of values in column major order, one can call the
  ! command MatSetOption(Mat A, MAT COLUMN ORIENTED);

  call MatSetOption(matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)


end subroutine myMatCreate

subroutine dRdwTMatMult(A, vecX,  vecY, ierr)
  use constants
  use communication
  use blockPointers
  use blockPointers_b
  use inputDiscretization 
  use inputTimeSpectral 
  use inputPhysics
  use iteration         
  use flowVarRefState     
  use inputAdjoint       
  use stencils
  use diffSizes
  use ADjointPETSc, only : dRdwT
  implicit none
#define PETSC_AVOID_MPIF_H
#include "finclude/petsc.h"
  
  Mat   A
  Vec   vecX, vecY
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii, sps2
  real(kind=realType) :: alpha, beta, force(3), moment(3), sepSensor, cavitation, cavitationb
  real(kind=realType) :: alphab, betab, forceb(3), momentb(3), sepSensorb
  real(kind=realType),pointer :: dwb_pointer(:)
  integer(kind=intType) :: nState, level, irow, liftIndex
  logical :: resetToRans

#ifndef USE_COMPLEX

  !call VecGetArrayF90(vecX, dwb_pointer, ierr)
  !call EChk(ierr,__FILE__,__LINE__)

  call VecSet(vecY, zero, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Setup number of state variable based on turbulence assumption
  if ( frozenTurbulence ) then
     nState = nwf
  else
     nState = nw
  end if
  ! This routine will not use the extra variables to block_res or the
  ! extra outputs, so we must zero them here
  call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)
  
  ! Need to trick the residual evalution to use coupled (mean flow and
  ! turbulent) together.
  
  level = 1
  currentLevel = level
  groundLevel = level

  ! If we are computing the jacobian for the RANS equations, we need
  ! to make block_res think that we are evauluating the residual in a
  ! fully coupled sense.  This is reset after this routine is
  ! finished.
  if (equations == RANSEquations) then
     nMGVar = nw
     nt1MG = nt1
     nt2MG = nt2

     turbSegregated = .False.
     turbCoupled = .True.
  end if

  ! Determine if we want to use frozenTurbulent Adjoint
  resetToRANS = .False. 
  if (frozenTurbulence .and. equations == RANSEquations) then
     equations = NSEquations 
     resetToRANS = .True.
  end if

  ! Zero the seeds
  alphab = zero
  betab = zero
  forceb = zero
  momentb = zero
  sepSensorb = zero
  cavitationb = zero

  do nn=1,nDom

     ! Just to get sizes
     call setPointers(nn,1_intType,1)
     call setDiffSizes
     
     ! Allocate the memory we need for this block to do the forward
     ! mode derivatives and copy reference values
     call alloc_derivative_values_bwd(nn, level)
        
     do sps=1,nTimeIntervalsSpectral
        ! Set pointers and derivative pointers
        call setPointers_b(nn, level, sps)

        ii = 0
        do k=2, kl
           do j=2,jl
              do i=2,il
                 do l = 1, nstate
                    ii = ii + 1
                    ! flowdomsb(nn,1,sps)%dw(i, j, k, l) = dwb_pointer(ii)
                    call VecGetValues(vecX, 1, (/ii-1/), &
                          flowdomsb(nn,1,sps)%dw(i, j, k, l), ierr)
                    call EChk(ierr,__FILE__,__LINE__)
                 end do
              end do
           end do
        end do

     call BLOCK_RES_B(nn, 1, .False., alpha, alphab, beta, betab, &
          & liftindex, force, forceb, moment, momentb, sepsensor, sepsensorb, &
          & cavitation, cavitationb)

        do sps2=1,nTimeIntervalsSpectral
           do k=2, kl
              do j=2,jl
                 do i=2,il
                    do l = 1, nstate
           
                       irow = flowDoms(nn, 1, sps2)%globalCell(i,j,k)*nstate + l -1
                       call VecSetValues(vecY, 1, (/irow/), &
                            flowdomsb(1,1,1)%w(i, j, k, l), ADD_VALUES, ierr)
                       call EChk(ierr,__FILE__,__LINE__)
                    end do
                 end do
              end do
           end do
        end do
     end do

     call dealloc_derivative_values_bwd(nn, 1)
  end do

  ! Reset the correct equation parameters if we were useing the frozen
  ! Turbulent 
  if (resetToRANS) then
     equations = RANSEquations
  end if

  ! Reset the paraters to use segrated turbulence solve. 
  if (equations == RANSEquations) then
     nMGVar = nwf
     nt1MG = nwf + 1
     nt2MG = nwf

     turbSegregated = .True.
     turbCoupled = .False.
     restrictEddyVis = .false.
     if( eddyModel ) restrictEddyVis = .true.
  end if
  
  call VecAssemblyBegin(vecY, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  call VecAssemblyEnd(vecY, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Assume no error
  ierr = 0

  !call VecRestoreArrayF90(vecX, dwb_pointer, ierr)
  !call EChk(ierr,__FILE__,__LINE__)

#endif
end subroutine dRdwTMatMult
