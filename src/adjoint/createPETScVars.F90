subroutine createPETScVars

#ifndef USE_NO_PETSC
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create the matrices/vectors that are required for the adjoint  *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc, only: dRdwT, dRdwPreT, &
       adjointKSP, matfreectx, x_like, psi_like1, adjointPETScVarsAllocated
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

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif


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
  external dRdwTMatMult, dRdwMatMult

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

  call getForceSize(npts, ncells)
  nDimPt = npts * 3 * nTimeIntervalsSpectral
  nDimCell = nCells * 3 * nTimeIntervalsSpectral

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
          level, .True.)
     call myMatCreate(dRdwT, nState, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
       __FILE__, __LINE__)
     
     call matSetOption(dRdwT, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     deallocate(nnzDiagonal, nnzOffDiag)
  else
       ! Setup matrix-free dRdwT
     call MatCreateShell(SUMB_COMM_WORLD, nDimW, nDimW, PETSC_DETERMINE, &
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
          level, .True.)
     call myMatCreate(dRdwPreT, nState, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
          __FILE__, __LINE__)

     call matSetOption(dRdwPreT, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     deallocate(nnzDiagonal, nnzOffDiag)
  end if 

  ! Create the KSP Object
  call KSPCreate(SUMB_COMM_WORLD, adjointKSP, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  adjointPETScVarsAllocated = .True.
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
#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

  Mat matrix
  integer(kind=intType), intent(in) :: blockSize, m, n
  integer(kind=intType), intent(in), dimension(*) :: nnzDiagonal, nnzOffDiag
  character*(*) :: file
  integer(kind=intType) :: ierr, line
  if (blockSize > 1) then
     call MatCreateBAIJ(SUMB_COMM_WORLD, blockSize, &
          m, n, PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal, 0, nnzOffDiag, matrix, ierr)
  else     
     call MatCreateAIJ(SUMB_COMM_WORLD,&
          m, n, PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal, 0, nnzOffDiag, matrix, ierr)
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
