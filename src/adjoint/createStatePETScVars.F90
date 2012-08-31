subroutine createStatePETScVars

#ifndef USE_NO_PETSC
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create the matrices/vectors that are required for the adjoint  *
  !     * solution: Matrices: dRdwT (dRdwPreT)                           *
  !     *           Vectors: dRdw,psi,adjointRes,adjointRHS,             *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc, only: dRdwT, drdwPreT, dJdw, psi, adjointRes, adjointRHS,&
       PETScIerr, PETScBlockMatrix
  use ADjointVars     ! nCellsLocal,nNodesLocal, nDesignExtra
  use communication   ! myID, nProc
  use inputTimeSpectral !nTimeIntervalsSpectral
  use flowVarRefState ! 
  use inputADjoint    !ApproxPC
  use stencils
  
  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  !
  !     Local variables.
  !
  integer(kind=intType)  :: nDimW,nDimX
  integer(kind=intType) :: i,n_stencil
  integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  integer(kind=intType), dimension(:), allocatable :: nnzDiagonal2, nnzOffDiag2
  integer(kind=intType), dimension(:,:), allocatable :: stencil

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !

  ! Define matrix dRdW local size, taking into account the total
  ! number of Cells owned by the processor and the number of 
  ! equations.

  nDimW = nw * nCellsLocal*nTimeIntervalsSpectral
  nDimX = 3 * nNodesLocal*nTimeIntervalsSpectral
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create matrix dRdW that define the adjoint linear system of    *
  !     * equations, dRdW^T psi = dJdW. Matrix dRdW has size [nDimW,nDimW]*
  !     * but is very sparse because of the computational stencil R=R(W).*
  !     *                                                                *
  !     ******************************************************************
  !

  ! ------------------- Determine Preallocation for dRdw --------------

  allocate( nnzDiagonal(nCellsLocal*nTimeIntervalsSpectral),&
       nnzOffDiag(nCellsLocal*nTimeIntervalsSpectral) )

  call initialize_stencils
  if (.not. viscous) then
     n_stencil = N_euler_drdw
     allocate(stencil(n_stencil,3))
     stencil = euler_drdw_stencil
  else
     n_stencil = N_visc_pc
     allocate(stencil(n_stencil,3))
     stencil = visc_pc_stencil
  end if

  call statePreAllocation(nnzDiagonal,nnzOffDiag,nDimW/nw,stencil,n_stencil)

  ! --------------------------------------------------------------------

  if( nw <= 7 ) then

     PETScBlockMatrix = .true.

     ! Create a Block AIJ Matrix with block size nw, (number of states)
     if (PETSC_VERSION_MINOR == 2) then
        call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
             nDimW, nDimW,                     &
             PETSC_DETERMINE, PETSC_DETERMINE, &
             0, nnzDiagonal,         &
             0, nnzOffDiag,            &
             dRdWT, PETScIerr)
     else
        call MatCreateBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
             nDimW, nDimW,                     &
             PETSC_DETERMINE, PETSC_DETERMINE, &
             0, nnzDiagonal,         &
             0, nnzOffDiag,            &
             dRdWT, PETScIerr)
     end if
     call EChk(PETScIerr,__FILE__,__LINE__)
  else

     PETScBlockMatrix = .false.

     allocate(nnzDiagonal2(nDimw),nnzOffDiag2(nDimw))
     ! The drdw prealloc function is done per block, which we can use
     ! to compute the correct preallocation for the non-block matrix:

     do i=1,nCellsLocal*nTimeIntervalsSpectral
        nnzDiagonal2((i-1)*nw+1:(i-1)*nw+nw) = nnzDiagonal(i)
        nnzOffDiag((i-1)*nw+1:(i-1)*nw+nw) = nnzOffDiag(i)
     end do
     if (PETSC_VERSION_MINOR == 2) then
        call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
             nDimW, nDimW,                     &
             PETSC_DETERMINE, PETSC_DETERMINE, &
             8, nnzDiagonal2,         &
             8, nnzOffDiag2,            &
             dRdWT, PETScIerr)
     else
        call MatCreateAIJ(SUMB_PETSC_COMM_WORLD,                 &
             nDimW, nDimW,                     &
             PETSC_DETERMINE, PETSC_DETERMINE, &
             8, nnzDiagonal2,         &
             8, nnzOffDiag2,            &
             dRdWT, PETScIerr)
     end if
     call EChk(PETScIerr,__FILE__,__LINE__)

     deallocate(nnzDiagonal2,nnzOffDiag2)
  endif

  deallocate(nnzDiagonal, nnzOffDiag, stencil)

  ! Set the matrix dRdW options.

  ! Warning: The array values is logically two-dimensional, 
  ! containing the values that are to be inserted. By default the
  ! values are given in row major order, which is the opposite of
  ! the Fortran convention, meaning that the value to be put in row
  ! idxm[i] and column idxn[j] is located in values[i*n+j]. To allow
  ! the insertion of values in column major order, one can call the
  ! command MatSetOption(Mat A,MAT COLUMN ORIENTED);

  call MatSetOption(dRdWt, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)


  !****************
  !create dRdWPre
  !***************

  if (ApproxPC) then

     !     ******************************************************************
     !     *                                                                *
     !     * Create matrix dRdWPre for the ADjoint approximate preconitioner*
     !     * This matrix is sparse with a narrower bandwidth than drdw      *
     !     *                                                                *
     !     ******************************************************************
     !


     ! ------------------- Determine Preallocation for dRdwPre -------------

     allocate( nnzDiagonal(nCellsLocal*nTimeIntervalsSpectral),&
          nnzOffDiag(nCellsLocal*nTimeIntervalsSpectral) )

     n_stencil = N_euler_PC
     allocate(stencil(n_stencil,3))
     stencil = euler_PC_stencil

     call statePreAllocation(nnzDiagonal,nnzOffDiag,nDimW/nw,stencil,n_stencil)

     ! --------------------------------------------------------------------

     if( nw <= 7 ) then

        PETScBlockMatrix = .true.
        if (PETSC_VERSION_MINOR == 2) then
           call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
                nDimW, nDimW,                     &
                PETSC_DETERMINE, PETSC_DETERMINE, &
                0, nnzDiagonal,         &
                0, nnzOffDiag,            &
                dRdWPreT, PETScIerr)
        else
           call MatCreateBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
                nDimW, nDimW,                     &
                PETSC_DETERMINE, PETSC_DETERMINE, &
                0, nnzDiagonal,         &
                0, nnzOffDiag,            &
                dRdWPreT, PETScIerr)
        end if
        call EChk(PETScIerr,__FILE__,__LINE__)
     else

        PETScBlockMatrix = .false.
        allocate(nnzDiagonal2(nDimw),nnzOffDiag2(nDimw))
        ! The drdwPC prealloc function is done per block, which we can use
        ! to compute the correct preallocation for the non-block matrix:

        do i=1,nCellsLocal*nTimeIntervalsSpectral
           nnzDiagonal2((i-1)*nw+1:(i-1)*nw+nw) = nnzDiagonal(i)
           nnzOffDiag((i-1)*nw+1:(i-1)*nw+nw) = nnzOffDiag(i)
        end do
        if (PETSC_VERSION_MINOR == 2) then
           call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
                nDimW, nDimW,                     &
                PETSC_DETERMINE, PETSC_DETERMINE, &
                0, nnzDiagonal2,         &
                0,nnzOffDiag2,            &
                dRdWPret, PETScIerr)
        else
           call MatCreateAIJ(SUMB_PETSC_COMM_WORLD,                 &
                nDimW, nDimW,                     &
                PETSC_DETERMINE, PETSC_DETERMINE, &
                0, nnzDiagonal2,         &
                0,nnzOffDiag2,            &
                dRdWPret, PETScIerr)
        end if
        call EChk(PETScIerr,__FILE__,__LINE__)

        deallocate(nnzDiagonal2,nnzOffDiag2)
     endif

     deallocate( nnzDiagonal, nnzOffDiag, stencil )

     ! Set the matrix dRdWPre options.

     call MatSetOption(dRdWPret, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)

  end if ! Approx PC

  ! Vectors:

  ! Get dJdw and psi from one MatGetVecs Call
  call MatGetVecs(dRdwT,dJdW,psi,PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  ! adjointRes is the same size as dJdw,psi
  call VecDuplicate(dJdW, adjointRes, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call VecDuplicate(dJdW, adjointRHS, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)


#endif

end subroutine createStatePETScVars
