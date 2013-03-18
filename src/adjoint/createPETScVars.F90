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
       adjointKSP, dFMdExtra, dRda_data
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
  integer(kind=intType) :: i, n_stencil
  integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  integer(kind=intType), dimension(:), allocatable :: nnzDiagonal2, nnzOffDiag2
  integer(kind=intType), dimension(:, :), pointer :: stencil
  integer(kind=intType) :: level, ierr, nlevels
  integer(kind=intType) :: rows(4), iCol, iCellCount, iNodeCount, nn, sps, ii
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, iDim, iStride, j, mm
  integer(kind=intType) :: npts, ncells, nTS

  ! DETERMINE ALL SIZES HERE!
  if ( frozenTurbulence ) then
     nDimW = nwf * nCellsLocal(1_intType)*nTimeIntervalsSpectral
  else
     nDimW = nw * nCellsLocal(1_intType)*nTimeIntervalsSpectral
  endif
  nDimX = 3 * nNodesLocal(1_intType)*nTimeIntervalsSpectral

  call getForceSize(npts, ncells, nTS)
  nDimPt = npts * 3 * nTS
  nDimCell = nCells * 3 * nTS

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
  if ( frozenTurbulence ) then
     call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nwf, stencil, n_stencil, &
          level)
       call myMatCreate(dRdwT, nwf, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
            __FILE__, __LINE__)
  else
     call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nw, stencil, n_stencil, &
          level)
     call myMatCreate(dRdwT, nw, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
          __FILE__, __LINE__)
  endif
     
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
     if ( frozenTurbulence ) then
        call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nwf, stencil, n_stencil, &
             level)
        call myMatCreate(dRdwPreT, nwf, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
             __FILE__, __LINE__)
    else
       call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nw, stencil, n_stencil, &
            level)
       call myMatCreate(dRdwPreT, nw, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
            __FILE__, __LINE__)
    endif

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

  ! Create the 6 vectors for d{F,M}/dw
  do i=1,6
     call VecDuplicate(dJdw, FMw(i), ierr)
     call EChk(ierr, __FILE__, __LINE__)
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
     nnzDiagonal = 2 * nw

     ! OffDiag is an overestimate...at most one cell is from neighbour
     ! proc
     nnzOffDiag  = 1 * nw 
  else
     nnzDiagonal = 3*3*2*nw

     ! In general should not have much more than 6 cells on off
     ! proc. If there is a malloc or two that isn't the end of the world
     nnzOffDiag = 6*nw
  end if

  call myMatCreate(dFcdw, 1, nDimCell, nDimw, nnzDiagonal, nnzOffDiag, &
       __FILE__, __LINE__)

  ! dFcdx
  if (.not. viscous) then
     nnzDiagonal = 4 * 3
     nnzOffDiag  = 1
  else
     nnzDiagonal = 9 * 3
     nnzOffDiag = 6*nw
  end if

  call myMatCreate(dFcdx, 1, nDimCell, nDimX, nnzDiagonal, nnzOffDiag, &
       __FILE__, __LINE__)

  deallocate(nnzDiagonal, nnzOffDiag)
  
  ! Finally we need dFndFc
  allocate( nnzDiagonal(nDimPt), nnzOffDiag(nDimPt))

  nnzDiagonal = 4
  nnzOffDiag  = 1 ! This should be enough...might get a couple of mallocs
  call myMatCreate(dFndFc, 1, nDimPt, nDimCell, nnzDiagonal, nnzOffDiag, &
       __FILE__, __LINE__)
  deallocate(nnzDiagonal, nnzOffDiag)

  ! We will also take this opportunity to assemble dFndFc. 
  iCellCount = 0
  iNodeCount = 0
 
  spectral: do sps=1,nTimeIntervalsSpectral
     domains: do nn=1,nDom
        call setPointers(nn,1_intType,sps)
        
        ! Loop over the number of boundary subfaces of this block.
        bocos: do mm=1,nBocos
           if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) then
              
              jBeg = BCData(mm)%jnBeg + 1; jEnd = BCData(mm)%jnEnd
              iBeg = BCData(mm)%inBeg + 1; iEnd = BCData(mm)%inEnd
       
              iStride = iEnd-iBeg+2
              do j=jBeg, jEnd ! Face Loop
                 do i=iBeg, iEnd ! Face Loop
                    do iDim = 0,2
                       
                       iCol = iCellCount*3 + iDim
                       rows(1) = 3*iNodeCount + 3*(j-2)*iStride + 3*(i-2) + iDim
                       rows(2) = 3*iNodeCount + 3*(j-1)*iStride + 3*(i-2) + iDim
                       rows(3) = 3*iNodeCount + 3*(j-2)*iStride + 3*(i-1) + iDim
                       rows(4) = 3*iNodeCount + 3*(j-1)*iStride + 3*(i-1) + iDim
                       do ii=1,4
                          call MatSetValues(dFndFc, 1, rows(ii), 1, iCol, &
                               fourth, INSERT_VALUES, ierr) 
                          call EChk(ierr, __FILE__, __LINE__)
                       end do

                    end do
                    iCellCount = iCellCount + 1
                 end do
              end do
             iNodeCount = iNodeCount + (iEnd-iBeg+2)*(jEnd-jBeg+2)
          end if
       end do bocos
    end do domains
 end do spectral

 call MatAssemblyBegin(dFndFc, MAT_FINAL_ASSEMBLY, ierr)
 call EChk(ierr, __FILE__, __LINE__)
 call MatAssemblyEnd  (dFndFc, MAT_FINAL_ASSEMBLY, ierr)
 call EChk(ierr,  __FILE__, __LINE__)

 ! For the tractions we also need dAdx. This is done in a similar fashion.



 ! For now, leave dFdw, and dFdx in.
 allocate( nnzDiagonal(nDimPt), nnzOffDiag(nDimPt) )
 nnzDiagonal = 8*nw
 nnzOffDiag  = 8*nw! Make the off diagonal the same

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
 
  ! Note we are creating the TRANPOSE of dRdx. It is size dDimX by nDimW
  call myMatCreate(dRdx, 1, nDimX, nDimW, nnzDiagonal, nnzOffDiag, &
       __FILE__, __LINE__)
  deallocate( nnzDiagonal, nnzOffDiag )
  
  call MatSetOption(dRdx, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call getForceSize(npts, ncells, nTS)

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

  ! Create the vectors for the FMx
  do i=1,6
     call VecDuplicate(dJdx, FMx(i), ierr)
     call EChk(ierr, __FILE__, __LINE__)
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
  allocate(dFMdExtra(6, nDesignExtra))

  if (PETSC_VERSION_MINOR < 3 ) then
     call MatCreateMPIDense(SUMB_COMM_WORLD, nDimW, PETSC_DECIDE, &
          PETSC_DETERMINE, nDesignExtra, dRda_data, dRda, ierr)
  else
     call MatCreateDense(SUMB_COMM_WORLD, nDimW, PETSC_DECIDE, &
          PETSC_DETERMINE, nDesignExtra, dRda_data, dRda, ierr)
  end if
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
  character*(*) :: file, line
  integer(kind=intType) :: ierr
  if (blockSize > 1) then
     if (PETSC_VERSION_MINOR <  3) then
        call MatCreateMPIBAIJ(SUMB_COMM_WORLD, blockSize, &
             m, n, PETSC_DETERMINE, PETSC_DETERMINE, &
             0, nnzDiagonal, 0, nnzOffDiag, matrix, ierr)
     else
        call MatCreateBAIJ(SUMB_COMM_WORLD, blockSize, &
             m, n, PETSC_DETERMINE, PETSC_DETERMINE, &
             0, nnzDiagonal, 0, nnzOffDiag, matrix, ierr)
     end if
  else
     if (PETSC_VERSION_MINOR <  3) then
        call MatCreateMPIAIJ(SUMB_COMM_WORLD, &
             m, n, PETSC_DETERMINE, PETSC_DETERMINE, &
             0, nnzDiagonal, 0, nnzOffDiag, matrix, ierr)
     else
        call MatCreateAIJ(SUMB_COMM_WORLD,&
             m, n, PETSC_DETERMINE, PETSC_DETERMINE, &
             0, nnzDiagonal, 0, nnzOffDiag, matrix, ierr)
     end if
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

end subroutine myMatCreate
