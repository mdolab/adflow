
!     ******************************************************************
!     *                                                                *
!     * File:          createPETScMat.F90                              *
!     * Author:        Andre C. Marta, C.A.(Sandy) Mader               *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 01-14-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine createPETScMat
!
!     ******************************************************************
!     *                                                                *
!     * Create the global matrix dRdW of the adjoint system of         *
!     * equations as a PETSc matrix object. Also, create the auxiliar  *
!     * matrices dRda and dRdx used to compute the gradients.          *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use ADjointVars     ! nCellsLocal,nNodesLocal, nDesignExtra
      use communication   ! myID, nProc
      use inputTimeSpectral !nTimeIntervalsSpectral
      use flowVarRefState ! 
      use inputADjoint    !ApproxPC
      implicit none
!
!     Local variables.
!
      ! nDim   - local dimension (processor owned)
      ! iLow   - first component owned by the local process
      ! iHigh  - one more than last component owned by the local process

      integer       :: nn, iLow, iHigh
      integer       :: nDimW, nDimX,nDimS
      integer       :: matBlockSize, matRows, matCols
      character(15) :: matTypeStr

      integer   :: nzDiagonalW, nzDiagonalX, nzOffDiag,nzDiagonalWPC
      integer, dimension(:), allocatable :: nnzDiagonal, nnzOffDiag

      character(len=2*maxStringLen) :: errorMessage
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
#ifndef USE_NO_PETSC

      ! PETSc macros are lost and have to be redefined.
      ! They were extracted from: <PETSc root dir>/include/petscmat.h
#define MATMPIDENSE        "mpidense"
#define MATMPIAIJ          "mpiaij"
      ! Define matrix dRdW local size, taking into account the total
      ! number of Cells owned by the processor and the number of 
      ! equations.

      nDimW = nw * nCellsLocal*nTimeIntervalsSpectral

      ! Define matrix dRdx local size (number of columns) for the
      ! spatial derivatives.

      nDimX = 3 * nNodesLocal*nTimeIntervalsSpectral

      ! Define matrix dSdx local size (number of Rows) for the
      ! Coupling derivatives.

      nDimS = 3 * nSurfNodesLocal

      ! Number of non-zero blocks per residual row in dRdW
      ! >>> This depends on the stencil being used R=R(W)
      !     1st order stencil ->  7 cells
      !     2nd order stencil -> 13 cells

      ! Stencil of W
      ! 1  - center cell
      ! 6  - 1st level cells along directions i,j,k
      ! 6  - 2nd level cells along directions i,j,k

      nzDiagonalW = 13+1*(nTimeIntervalsSpectral-1)
      !should be addition not multiplication but this will require filtering 
      !of some of the TS results for non-zero terms. 
      !i.e.:13 +nTimeIntervalsSpectral! 1 + 6 + 6  check!!!

      ! Stencil of W pc
      ! 1  - center cell
      ! 6  - 1st level cells along directions i,j,k

      nzDiagonalWPC = 7 *nTimeIntervalsSpectral
      !should be addition not multiplication but this will require filtering 
      !of some of the TS results for non-zero terms. 
      !i.e.:7 +nTimeIntervalsSpectral! 1 + 6  check!!!

      ! Stencil of X
      ! 1  - center node
      ! 6  - 1st level nodes along directions i,j,k
      ! 6  - 2nd level nodes along directions i,j,k
      ! 12 - 1st level nodes along diagonals (i,j),(i,k),(j,k) 

      nzDiagonalX = 33+8*(nTimeIntervalsSpectral-1)!25+4*nTimeIntervalsSpectral! 1 + 6 + 6 + 12 Check

      ! Average number of off processor contributions per Cell
      ! (average number of donor cells that come from other processor)

      nzOffDiag  = 12+5*(nTimeIntervalsSpectral-1)
!
!     ******************************************************************
!     *                                                                *
!     * Create matrix dRdW that define the adjoint linear system of    *
!     * equations, dRdW^T psi = dJdW. Matrix dRdW has size [nDim,nDim] *
!     * but is very sparse because of the computational stencil R=R(W).*
!     *                                                                *
!     ******************************************************************
!
      !call PetscOptionsSetValue('-mat_view_info',PETSC_NULL_CHARACTER,PETScIerr)

      ! Create the matrix dRdW.

      ! >>> option #1 : sparse parallel matrix in block AIJ format
      !                 Only works for block up to size 7...

      if( nw <= 7 ) then
!      if( nw == 7 ) then

         PETScBlockMatrix = .true.

        ! MatCreateMPIBAIJ - Creates a sparse parallel matrix in block
        !   AIJ format (block compressed row). For good matrix assembly
        !   performance the user should preallocate the matrix storage
        !   by setting the parameters d_nz (or d_nnz) and o_nz (or
        !   o_nnz). By setting these parameters accurately, performance
        !   can be increased by more than a factor of 50.
        !
        ! Synopsis
        !
        ! #include "petscmat.h"  
        ! call MatCreateMPIBAIJ(MPI_Comm comm,PetscInt bs,             &
        !                 PetscInt m,PetscInt n,PetscInt M,PetscInt N, &
        !                 PetscInt d_nz,const PetscInt d_nnz[],        &
        !                 PetscInt o_nz,const PetscInt o_nnz[],        &
        !                 Mat *A,PetscErrorCode ierr)
        !
        ! Collective on MPI_Comm
        !
        ! Input Parameters
        !   comm  - MPI communicator
        !   bs    - size of block
        !   m     - number of local rows (or PETSC_DECIDE to have 
        !           calculated if M is given). This value should be the
        !           same as the local size used in creating the y
        !           vector for the matrix-vector product y = Ax.
        !   n     - number of local columns (or PETSC_DECIDE to have
        !           calculated if N is given). This value should be the
        !           same as the local size used in creating the x
        !           vector for the matrix-vector product y = Ax.
        !   M     - number of global rows (or PETSC_DETERMINE to have
        !           calculated if m is given).
        !   N     - number of global columns (or PETSC_DETERMINE to
        !           have calculated if n is given).
        !   d_nz  - number of nonzero blocks per block row in diagonal
        !           portion of local submatrix (same for all local rows)
        !   d_nnz - array containing the number of nonzero blocks in the
        !           various block rows of the in diagonal portion of the
        !           local (possibly different for each block row) or
        !           PETSC_NULL. You must leave room for the diagonal
        !           entry even if it is zero.
        !   o_nz  - number of nonzero blocks per block row in the
        !           off-diagonal portion of local submatrix (same for
        !           all local rows).
        !   o_nnz - array containing the number of nonzero blocks in the
        !           various block rows of the off-diagonal portion of
        !           the local submatrix (possibly different for each
        !           block row) or PETSC_NULL.
        !
        ! Output Parameter
        !   A     - the matrix 
        !
        ! See .../petsc/docs/manualpages/Mat/MatCreateMPIBAIJ.html
  
        allocate( nnzDiagonal(nCellsLocal*nTimeIntervalsSpectral),&
             nnzOffDiag(nCellsLocal*nTimeIntervalsSpectral) )

        nnzDiagonal = nzDiagonalW
        nnzOffDiag  = nzOffDiag
        !print *,'nnzDiagonal',nnzDiagonal,'ofdiag',nnzOffDiag 
!!$        call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
!!$                              nDimW, nDimW,                     &
!!$                              PETSC_DETERMINE, PETSC_DETERMINE, &
!!$                              nzDiagonalW, nnzDiagonal,         &
!!$                              nzOffDiag, nnzOffDiag,            &
!!$                              dRdW, PETScIerr)
        call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
                              nDimW, nDimW,                     &
                              PETSC_DETERMINE, PETSC_DETERMINE, &
                              nzDiagonalW, nnzDiagonal,         &
                              nzOffDiag, nnzOffDiag,            &
                              dRdWT, PETScIerr)

      ! >>> option #2 : sparse parallel matrix in AIJ format
      !                 General case...

      else

         PETScBlockMatrix = .false.

        ! MatCreateMPIAIJ - Creates a sparse parallel matrix in AIJ
        !   format (the default parallel PETSc format). For good matrix
        !   assembly performance the user should preallocate the matrix
        !   storage by setting the parameters d_nz (or d_nnz) and o_nz
        !   (or o_nnz). By setting these parameters accurately,
        !   performance can be increased by more than a factor of 50.
        !
        ! Synopsis
        !
        ! #include "petscmat.h" 
        ! call MatCreateMPIAIJ(MPI_Comm comm,                           &
        !                  PetscInt m,PetscInt n,PetscInt M,PetscInt N, &
        !                  PetscInt d_nz,const PetscInt d_nnz[],        &
        !                  PetscInt o_nz,const PetscInt o_nnz[],        &
        !                  Mat *A, PetscErrorCode ierr)
        !
        ! Collective on MPI_Comm
        !
        ! Input Parameters
        !   comm  - MPI communicator
        !   m     - number of local rows (or PETSC_DECIDE to have
        !           calculated if M is given) This value should be the
        !           same as the local size used in creating the y vector
        !           for the matrix-vector product y = Ax.
        !   n     - This value should be the same as the local size used
        !           in creating the x vector for the matrix-vector
        !           product y = Ax. (or PETSC_DECIDE to have calculated
        !           if N is given) For square matrices n is almost
        !           always m.
        !   M     - number of global rows (or PETSC_DETERMINE to have
        !           calculated if m is given)
        !   N     - number of global columns (or PETSC_DETERMINE to have
        !           calculated if n is given)
        !   d_nz  - number of nonzeros per row in DIAGONAL portion of
        !           local submatrix (same value is used for all local
        !           rows)
        !   d_nnz - array containing the number of nonzeros in the
        !           various rows of the DIAGONAL portion of the local
        !           submatrix (possibly different for each row) or
        !           PETSC_NULL, if d_nz is used to specify the nonzero
        !           structure. The size of this array is equal to the
        !           number of local rows, i.e 'm'. You must leave room
        !           for the diagonal entry even if it is zero.
        !   o_nz  - number of nonzeros per row in the OFF-DIAGONAL
        !           portion of local submatrix (same value is used for
        !           all local rows).
        !   o_nnz - array containing the number of nonzeros in the
        !           various rows of the OFF-DIAGONAL portion of the
        !           local submatrix (possibly different for each row) or
        !           PETSC_NULL, if o_nz is used to specify the nonzero
        !           structure. The size of this array is equal to the
        !           number of local rows, i.e 'm'.
        !
        ! Output Parameter
        !   A     - the matrix
        !
        ! Notes
        ! The parallel matrix is partitioned such that the first m0 rows
        !  belong to process 0, the next m1 rows belong to process 1,
        !  the next m2 rows belong to process 2 etc.. where m0,m1,m2...
        !  are the input parameter 'm'.
        !
        ! The DIAGONAL portion of the local submatrix of a processor can
        !  be defined as the submatrix which is obtained by extraction
        !  the part corresponding to the rows r1-r2 and columns r1-r2 of
        !  the global matrix, where r1 is the first row that belongs to
        !  the processor, and r2 is the last row belonging to the this
        !  processor. This is a square mxm matrix. The remaining portion
        !  of the local submatrix (mxN) constitute the OFF-DIAGONAL
        !  portion.
        !
        ! If o_nnz, d_nnz are specified, then o_nz, and d_nz are ignored.
        !
        ! When calling this routine with a single process communicator,
        !  a matrix of type SEQAIJ is returned.
        !
        ! See .../petsc/docs/manualpages/Mat/MatCreateMPIAIJ.html

        nzDiagonalW = nzDiagonalW * nw
        nzOffDiag   = nzOffDiag   * nw

        allocate( nnzDiagonal(nDimW), nnzOffDiag(nDimW) )

        nnzDiagonal = nzDiagonalW
        nnzOffDiag  = nzOffDiag

!!$        call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
!!$                             nDimW, nDimW,                     &
!!$                             PETSC_DETERMINE, PETSC_DETERMINE, &
!!$                             nzDiagonalW, nnzDiagonal,         &
!!$                             nzOffDiag, nnzOffDiag,            &
!!$                             dRdW, PETScIerr)

        call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
                             nDimW, nDimW,                     &
                             PETSC_DETERMINE, PETSC_DETERMINE, &
                             nzDiagonalW, nnzDiagonal,         &
                             nzOffDiag, nnzOffDiag,            &
                             dRdWT, PETScIerr)

      endif

      deallocate( nnzDiagonal, nnzOffDiag )

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
                     "Could not create matrix dRdW of local size", nDimW
        call terminate("createPETScMat", errorMessage)
      endif

      ! Set the matrix dRdW options.

      ! Warning: The array values is logically two-dimensional, 
      ! containing the values that are to be inserted. By default the
      ! values are given in row major order, which is the opposite of
      ! the Fortran convention, meaning that the value to be put in row
      ! idxm[i] and column idxn[j] is located in values[i*n+j]. To allow
      ! the insertion of values in column major order, one can call the
      ! command MatSetOption(Mat A,MAT COLUMN ORIENTED);

      ! MatSetOption - Sets a parameter option for a matrix.
      !   Some options may be specific to certain storage formats.
      !   Some options determine how values will be inserted (or added).
      !   Sorted, row-oriented input will generally assemble the fastest.
      !   The default is row-oriented, nonsorted input.
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! call MatSetOption(Mat mat,MatOption op,PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat    - the matrix
      !   option - the option, one of those listed below (and possibly
      !     others), e.g., MAT_ROWS_SORTED, MAT_NEW_NONZERO_LOCATION_ERR
      !
      ! Options For Use with MatSetValues()
      ! Insert a logically dense subblock, which can be
      !   MAT_ROW_ORIENTED     - row-oriented (default)
      !   MAT_COLUMN_ORIENTED  - column-oriented
      !   MAT_ROWS_SORTED      - sorted by row
      !   MAT_ROWS_UNSORTED    - not sorted by row (default)
      !   MAT_COLUMNS_SORTED   - sorted by column
      !   MAT_COLUMNS_UNSORTED - not sorted by column (default)
      !
      ! Note these options reflect the data you pass in with
      !   MatSetValues(); it has nothing to do with how the data
      !   is stored internally in the matrix data structure.
      !
      ! When (re)assembling a matrix, we can restrict the input for
      !   efficiency/debugging purposes. These options include
      !     MAT_NO_NEW_NONZERO_LOCATIONS  - additional insertions will
      !       not be allowed if they generate a new nonzero
      !     MAT_YES_NEW_NONZERO_LOCATIONS - additional insertions will
      !       be allowed
      !     MAT_NO_NEW_DIAGONALS          - additional insertions will
      !       not be allowed if they generate a nonzero in a new
      !       diagonal (for block diagonal format only)
      !     MAT_YES_NEW_DIAGONALS         - new diagonals will be
      !       allowed (for block diagonal format only)
      !     MAT_IGNORE_OFF_PROC_ENTRIES   - drops off-processor entries
      !     MAT_NEW_NONZERO_LOCATION_ERR  - generates an error for new
      !       matrix entry
      !     MAT_USE_HASH_TABLE            - uses a hash table to speed
      !       up matrix assembly
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetOption.html
      ! or PETSc users manual, pp.51-52


#ifdef USE_PETSC_3
    !  call MatSetOption(dRdW, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)
      call MatSetOption(dRdWt, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)
      !call MatSetOption(dRdWt,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE,PETScIErr)
      !call MatSetOption(dRdWt,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE,PETScIErr)
      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dRdW")
#else
  !    call MatSetOption(dRdW, MAT_COLUMN_ORIENTED, PETScIerr)
      call MatSetOption(dRdWt, MAT_COLUMN_ORIENTED, PETScIerr)
      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dRdW")
#endif
!****************
!create dRdWPre
!***************

if (ApproxPC) then

!
!     ******************************************************************
!     *                                                                *
!     * Create matrix dRdWPre for the ADjoint approximate preconitioner*
!     * This matrix is sparse with a narrower bandwidth than drdw      *
!     *                                                                *
!     ******************************************************************
!
      ! Create the matrix dRdWPre.

      ! >>> option #1 : sparse parallel matrix in block AIJ format
      !                 Only works for block up to size 7...

      if( nw <= 7 ) then
!      if( nw == 7 ) then

         PETScBlockMatrix = .true.

        ! MatCreateMPIBAIJ - Creates a sparse parallel matrix in block
        !   AIJ format (block compressed row). For good matrix assembly
        !   performance the user should preallocate the matrix storage
        !   by setting the parameters d_nz (or d_nnz) and o_nz (or
        !   o_nnz). By setting these parameters accurately, performance
        !   can be increased by more than a factor of 50.
        !
        ! Synopsis
        !
        ! #include "petscmat.h"  
        ! call MatCreateMPIBAIJ(MPI_Comm comm,PetscInt bs,             &
        !                 PetscInt m,PetscInt n,PetscInt M,PetscInt N, &
        !                 PetscInt d_nz,const PetscInt d_nnz[],        &
        !                 PetscInt o_nz,const PetscInt o_nnz[],        &
        !                 Mat *A,PetscErrorCode ierr)
        !
        ! Collective on MPI_Comm
        !
        ! Input Parameters
        !   comm  - MPI communicator
        !   bs    - size of block
        !   m     - number of local rows (or PETSC_DECIDE to have 
        !           calculated if M is given). This value should be the
        !           same as the local size used in creating the y
        !           vector for the matrix-vector product y = Ax.
        !   n     - number of local columns (or PETSC_DECIDE to have
        !           calculated if N is given). This value should be the
        !           same as the local size used in creating the x
        !           vector for the matrix-vector product y = Ax.
        !   M     - number of global rows (or PETSC_DETERMINE to have
        !           calculated if m is given).
        !   N     - number of global columns (or PETSC_DETERMINE to
        !           have calculated if n is given).
        !   d_nz  - number of nonzero blocks per block row in diagonal
        !           portion of local submatrix (same for all local rows)
        !   d_nnz - array containing the number of nonzero blocks in the
        !           various block rows of the in diagonal portion of the
        !           local (possibly different for each block row) or
        !           PETSC_NULL. You must leave room for the diagonal
        !           entry even if it is zero.
        !   o_nz  - number of nonzero blocks per block row in the
        !           off-diagonal portion of local submatrix (same for
        !           all local rows).
        !   o_nnz - array containing the number of nonzero blocks in the
        !           various block rows of the off-diagonal portion of
        !           the local submatrix (possibly different for each
        !           block row) or PETSC_NULL.
        !
        ! Output Parameter
        !   A     - the matrix 
        !
        ! See .../petsc/docs/manualpages/Mat/MatCreateMPIBAIJ.html
  
        allocate( nnzDiagonal(nCellsLocal*nTimeIntervalsSpectral),&
             nnzOffDiag(nCellsLocal*nTimeIntervalsSpectral) )

        nnzDiagonal = nzDiagonalWPC
        nnzOffDiag  = nzOffDiag
        !print *,'nnzDiagonal',nnzDiagonal,'ofdiag',nnzOffDiag 
!!$        call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
!!$                              nDimW, nDimW,                     &
!!$                              PETSC_DETERMINE, PETSC_DETERMINE, &
!!$                              nzDiagonalW, nnzDiagonal,         &
!!$                              nzOffDiag, nnzOffDiag,            &
!!$                              dRdWPre, PETScIerr)
        call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
                              nDimW, nDimW,                     &
                              PETSC_DETERMINE, PETSC_DETERMINE, &
                              nzDiagonalW, nnzDiagonal,         &
                              nzOffDiag, nnzOffDiag,            &
                              dRdWPreT, PETScIerr)

      ! >>> option #2 : sparse parallel matrix in AIJ format
      !                 General case...

      else

         PETScBlockMatrix = .false.

        ! MatCreateMPIAIJ - Creates a sparse parallel matrix in AIJ
        !   format (the default parallel PETSc format). For good matrix
        !   assembly performance the user should preallocate the matrix
        !   storage by setting the parameters d_nz (or d_nnz) and o_nz
        !   (or o_nnz). By setting these parameters accurately,
        !   performance can be increased by more than a factor of 50.
        !
        ! Synopsis
        !
        ! #include "petscmat.h" 
        ! call MatCreateMPIAIJ(MPI_Comm comm,                           &
        !                  PetscInt m,PetscInt n,PetscInt M,PetscInt N, &
        !                  PetscInt d_nz,const PetscInt d_nnz[],        &
        !                  PetscInt o_nz,const PetscInt o_nnz[],        &
        !                  Mat *A, PetscErrorCode ierr)
        !
        ! Collective on MPI_Comm
        !
        ! Input Parameters
        !   comm  - MPI communicator
        !   m     - number of local rows (or PETSC_DECIDE to have
        !           calculated if M is given) This value should be the
        !           same as the local size used in creating the y vector
        !           for the matrix-vector product y = Ax.
        !   n     - This value should be the same as the local size used
        !           in creating the x vector for the matrix-vector
        !           product y = Ax. (or PETSC_DECIDE to have calculated
        !           if N is given) For square matrices n is almost
        !           always m.
        !   M     - number of global rows (or PETSC_DETERMINE to have
        !           calculated if m is given)
        !   N     - number of global columns (or PETSC_DETERMINE to have
        !           calculated if n is given)
        !   d_nz  - number of nonzeros per row in DIAGONAL portion of
        !           local submatrix (same value is used for all local
        !           rows)
        !   d_nnz - array containing the number of nonzeros in the
        !           various rows of the DIAGONAL portion of the local
        !           submatrix (possibly different for each row) or
        !           PETSC_NULL, if d_nz is used to specify the nonzero
        !           structure. The size of this array is equal to the
        !           number of local rows, i.e 'm'. You must leave room
        !           for the diagonal entry even if it is zero.
        !   o_nz  - number of nonzeros per row in the OFF-DIAGONAL
        !           portion of local submatrix (same value is used for
        !           all local rows).
        !   o_nnz - array containing the number of nonzeros in the
        !           various rows of the OFF-DIAGONAL portion of the
        !           local submatrix (possibly different for each row) or
        !           PETSC_NULL, if o_nz is used to specify the nonzero
        !           structure. The size of this array is equal to the
        !           number of local rows, i.e 'm'.
        !
        ! Output Parameter
        !   A     - the matrix
        !
        ! Notes
        ! The parallel matrix is partitioned such that the first m0 rows
        !  belong to process 0, the next m1 rows belong to process 1,
        !  the next m2 rows belong to process 2 etc.. where m0,m1,m2...
        !  are the input parameter 'm'.
        !
        ! The DIAGONAL portion of the local submatrix of a processor can
        !  be defined as the submatrix which is obtained by extraction
        !  the part corresponding to the rows r1-r2 and columns r1-r2 of
        !  the global matrix, where r1 is the first row that belongs to
        !  the processor, and r2 is the last row belonging to the this
        !  processor. This is a square mxm matrix. The remaining portion
        !  of the local submatrix (mxN) constitute the OFF-DIAGONAL
        !  portion.
        !
        ! If o_nnz, d_nnz are specified, then o_nz, and d_nz are ignored.
        !
        ! When calling this routine with a single process communicator,
        !  a matrix of type SEQAIJ is returned.
        !
        ! See .../petsc/docs/manualpages/Mat/MatCreateMPIAIJ.html

        nzDiagonalW = nzDiagonalW * nw
        nzOffDiag   = nzOffDiag   * nw

        allocate( nnzDiagonal(nDimW), nnzOffDiag(nDimW) )

        nnzDiagonal = nzDiagonalWPC
        nnzOffDiag  = nzOffDiag

!!$        call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
!!$                             nDimW, nDimW,                     &
!!$                             PETSC_DETERMINE, PETSC_DETERMINE, &
!!$                             nzDiagonalW, nnzDiagonal,         &
!!$                             nzOffDiag, nnzOffDiag,            &
!!$                             dRdWPre, PETScIerr)
        call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
                             nDimW, nDimW,                     &
                             PETSC_DETERMINE, PETSC_DETERMINE, &
                             nzDiagonalW, nnzDiagonal,         &
                             nzOffDiag, nnzOffDiag,            &
                             dRdWPret, PETScIerr)

      endif

      deallocate( nnzDiagonal, nnzOffDiag )

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
                     "Could not create matrix dRdWPre of local size", nDimW
        call terminate("createPETScMat", errorMessage)
      endif

      ! Set the matrix dRdWPre options.

      ! Warning: The array values is logically two-dimensional, 
      ! containing the values that are to be inserted. By default the
      ! values are given in row major order, which is the opposite of
      ! the Fortran convention, meaning that the value to be put in row
      ! idxm[i] and column idxn[j] is located in values[i*n+j]. To allow
      ! the insertion of values in column major order, one can call the
      ! command MatSetOption(Mat A,MAT COLUMN ORIENTED);

      ! MatSetOption - Sets a parameter option for a matrix.
      !   Some options may be specific to certain storage formats.
      !   Some options determine how values will be inserted (or added).
      !   Sorted, row-oriented input will generally assemble the fastest.
      !   The default is row-oriented, nonsorted input.
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! call MatSetOption(Mat mat,MatOption op,PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat    - the matrix
      !   option - the option, one of those listed below (and possibly
      !     others), e.g., MAT_ROWS_SORTED, MAT_NEW_NONZERO_LOCATION_ERR
      !
      ! Options For Use with MatSetValues()
      ! Insert a logically dense subblock, which can be
      !   MAT_ROW_ORIENTED     - row-oriented (default)
      !   MAT_COLUMN_ORIENTED  - column-oriented
      !   MAT_ROWS_SORTED      - sorted by row
      !   MAT_ROWS_UNSORTED    - not sorted by row (default)
      !   MAT_COLUMNS_SORTED   - sorted by column
      !   MAT_COLUMNS_UNSORTED - not sorted by column (default)
      !
      ! Note these options reflect the data you pass in with
      !   MatSetValues(); it has nothing to do with how the data
      !   is stored internally in the matrix data structure.
      !
      ! When (re)assembling a matrix, we can restrict the input for
      !   efficiency/debugging purposes. These options include
      !     MAT_NO_NEW_NONZERO_LOCATIONS  - additional insertions will
      !       not be allowed if they generate a new nonzero
      !     MAT_YES_NEW_NONZERO_LOCATIONS - additional insertions will
      !       be allowed
      !     MAT_NO_NEW_DIAGONALS          - additional insertions will
      !       not be allowed if they generate a nonzero in a new
      !       diagonal (for block diagonal format only)
      !     MAT_YES_NEW_DIAGONALS         - new diagonals will be
      !       allowed (for block diagonal format only)
      !     MAT_IGNORE_OFF_PROC_ENTRIES   - drops off-processor entries
      !     MAT_NEW_NONZERO_LOCATION_ERR  - generates an error for new
      !       matrix entry
      !     MAT_USE_HASH_TABLE            - uses a hash table to speed
      !       up matrix assembly
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetOption.html
      ! or PETSc users manual, pp.51-52
#ifdef USE_PETSC_3
   !   call MatSetOption(dRdWPre, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)
      call MatSetOption(dRdWPret, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)
      !call MatSetOption(dRdWPret,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE,PETScIErr)
     ! call MatSetOption(dRdWPret,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE,PETScIErr)
      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dRdW")
#else
    !  call MatSetOption(dRdWPre, MAT_COLUMN_ORIENTED, PETScIerr)
      call MatSetOption(dRdWPret, MAT_COLUMN_ORIENTED, PETScIerr)
      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dRdWPre")
#endif
   end if

if(Debug) then

      !******************************************
      !Create dRdWFD for debugging
      !*****************************************

       ! MatCreateMPIAIJ - Creates a sparse parallel matrix in AIJ
        !   format (the default parallel PETSc format). For good matrix
        !   assembly performance the user should preallocate the matrix
        !   storage by setting the parameters d_nz (or d_nnz) and o_nz
        !   (or o_nnz). By setting these parameters accurately,
        !   performance can be increased by more than a factor of 50.
        !
        ! Synopsis
        !
        ! #include "petscmat.h" 
        ! call MatCreateMPIAIJ(MPI_Comm comm,                           &
        !                  PetscInt m,PetscInt n,PetscInt M,PetscInt N, &
        !                  PetscInt d_nz,const PetscInt d_nnz[],        &
        !                  PetscInt o_nz,const PetscInt o_nnz[],        &
        !                  Mat *A, PetscErrorCode ierr)
        !
        ! Collective on MPI_Comm
        !
        ! Input Parameters
        !   comm  - MPI communicator
        !   m     - number of local rows (or PETSC_DECIDE to have
        !           calculated if M is given) This value should be the
        !           same as the local size used in creating the y vector
        !           for the matrix-vector product y = Ax.
        !   n     - This value should be the same as the local size used
        !           in creating the x vector for the matrix-vector
        !           product y = Ax. (or PETSC_DECIDE to have calculated
        !           if N is given) For square matrices n is almost
        !           always m.
        !   M     - number of global rows (or PETSC_DETERMINE to have
        !           calculated if m is given)
        !   N     - number of global columns (or PETSC_DETERMINE to have
        !           calculated if n is given)
        !   d_nz  - number of nonzeros per row in DIAGONAL portion of
        !           local submatrix (same value is used for all local
        !           rows)
        !   d_nnz - array containing the number of nonzeros in the
        !           various rows of the DIAGONAL portion of the local
        !           submatrix (possibly different for each row) or
        !           PETSC_NULL, if d_nz is used to specify the nonzero
        !           structure. The size of this array is equal to the
        !           number of local rows, i.e 'm'. You must leave room
        !           for the diagonal entry even if it is zero.
        !   o_nz  - number of nonzeros per row in the OFF-DIAGONAL
        !           portion of local submatrix (same value is used for
        !           all local rows).
        !   o_nnz - array containing the number of nonzeros in the
        !           various rows of the OFF-DIAGONAL portion of the
        !           local submatrix (possibly different for each row) or
        !           PETSC_NULL, if o_nz is used to specify the nonzero
        !           structure. The size of this array is equal to the
        !           number of local rows, i.e 'm'.
        !
        ! Output Parameter
        !   A     - the matrix
        !
        ! Notes
        ! The parallel matrix is partitioned such that the first m0 rows
        !  belong to process 0, the next m1 rows belong to process 1,
        !  the next m2 rows belong to process 2 etc.. where m0,m1,m2...
        !  are the input parameter 'm'.
        !
        ! The DIAGONAL portion of the local submatrix of a processor can
        !  be defined as the submatrix which is obtained by extraction
        !  the part corresponding to the rows r1-r2 and columns r1-r2 of
        !  the global matrix, where r1 is the first row that belongs to
        !  the processor, and r2 is the last row belonging to the this
        !  processor. This is a square mxm matrix. The remaining portion
        !  of the local submatrix (mxN) constitute the OFF-DIAGONAL
        !  portion.
        !
        ! If o_nnz, d_nnz are specified, then o_nz, and d_nz are ignored.
        !
        ! When calling this routine with a single process communicator,
        !  a matrix of type SEQAIJ is returned.
        !
        ! See .../petsc/docs/manualpages/Mat/MatCreateMPIAIJ.html

        nzDiagonalW = nzDiagonalW * nw
        nzOffDiag   = nzOffDiag   * nw

        allocate( nnzDiagonal(nDimW), nnzOffDiag(nDimW) )

        nnzDiagonal = nzDiagonalW
        nnzOffDiag  = nzOffDiag

        call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
                             nDimW, nDimW,                     &
                             PETSC_DETERMINE, PETSC_DETERMINE, &
                             nzDiagonalW, nnzDiagonal,         &
                             nzOffDiag, nnzOffDiag,            &
                             dRdWFD, PETScIerr)

      deallocate( nnzDiagonal, nnzOffDiag )

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
                     "Could not create matrix dRdWFD of local size", nDimW
        call terminate("createPETScMat", errorMessage)
      endif

      ! Set the matrix dRdWFD options.

      ! Warning: The array values is logically two-dimensional, 
      ! containing the values that are to be inserted. By default the
      ! values are given in row major order, which is the opposite of
      ! the Fortran convention, meaning that the value to be put in row
      ! idxm[i] and column idxn[j] is located in values[i*n+j]. To allow
      ! the insertion of values in column major order, one can call the
      ! command MatSetOption(Mat A,MAT COLUMN ORIENTED);

      ! MatSetOption - Sets a parameter option for a matrix.
      !   Some options may be specific to certain storage formats.
      !   Some options determine how values will be inserted (or added).
      !   Sorted, row-oriented input will generally assemble the fastest.
      !   The default is row-oriented, nonsorted input.
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! call MatSetOption(Mat mat,MatOption op,PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat    - the matrix
      !   option - the option, one of those listed below (and possibly
      !     others), e.g., MAT_ROWS_SORTED, MAT_NEW_NONZERO_LOCATION_ERR
      !
      ! Options For Use with MatSetValues()
      ! Insert a logically dense subblock, which can be
      !   MAT_ROW_ORIENTED     - row-oriented (default)
      !   MAT_COLUMN_ORIENTED  - column-oriented
      !   MAT_ROWS_SORTED      - sorted by row
      !   MAT_ROWS_UNSORTED    - not sorted by row (default)
      !   MAT_COLUMNS_SORTED   - sorted by column
      !   MAT_COLUMNS_UNSORTED - not sorted by column (default)
      !
      ! Note these options reflect the data you pass in with
      !   MatSetValues(); it has nothing to do with how the data
      !   is stored internally in the matrix data structure.
      !
      ! When (re)assembling a matrix, we can restrict the input for
      !   efficiency/debugging purposes. These options include
      !     MAT_NO_NEW_NONZERO_LOCATIONS  - additional insertions will
      !       not be allowed if they generate a new nonzero
      !     MAT_YES_NEW_NONZERO_LOCATIONS - additional insertions will
      !       be allowed
      !     MAT_NO_NEW_DIAGONALS          - additional insertions will
      !       not be allowed if they generate a nonzero in a new
      !       diagonal (for block diagonal format only)
      !     MAT_YES_NEW_DIAGONALS         - new diagonals will be
      !       allowed (for block diagonal format only)
      !     MAT_IGNORE_OFF_PROC_ENTRIES   - drops off-processor entries
      !     MAT_NEW_NONZERO_LOCATION_ERR  - generates an error for new
      !       matrix entry
      !     MAT_USE_HASH_TABLE            - uses a hash table to speed
      !       up matrix assembly
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetOption.html
      ! or PETSc users manual, pp.51-52
#ifdef USE_PETSC_3
      call MatSetOption(dRdWFD, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dRdW")
#else
      call MatSetOption(dRdWFD, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dRdWFD")
#endif
      !*****************************************
      ! end of create dRdWFD
      !*****************************************
   end if

      ! Extract info from the global matrix (only processor 0 does it).

      if( PETScRank==0 .and. debug ) then
    !  if( PETScRank==0 ) then

        ! Get the matrix block size.

        ! MatGetBlockSize - Returns the matrix block size; useful
        !   especially for the block row and block diagonal formats.
        !
        ! Synopsis
        !
        ! #include "petscmat.h" 
        ! call MatGetBlockSize(Mat mat,PetscInt *bs,PetscErrorCode ierr)
        !
        ! Not Collective
        !
        ! Input Parameter
        !   mat - the matrix
        !
        ! Output Parameter
        !   bs  - block size
        !
        ! see .../petsc/docs/manualpages/Mat/MatGetBlockSize.html

        call MatGetBlockSize(dRdWt, matBlockSize, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetBlockSize dRdW")

        write(*,10) "# MATRIX: dRdWt block size  =", matBlockSize

        ! Get the global number of rows and columns.

        ! MatGetSize - Returns the numbers of rows and
        !              columns in a matrix.
        ! Synopsis
        !
        ! #include "petscmat.h" 
        ! call MatGetSize(Mat mat,PetscInt *m,PetscInt* n, &
        !                 PetscErrorCode ierr)
        !
        ! Not Collective
        !
        ! Input Parameter
        !   mat - the matrix
        !
        ! Output Parameters
        !   m   - the number of global rows
        !   n   - the number of global columns
        !
        ! see .../petsc/docs/manualpages/Mat/MatGetSize.html

        call MatGetSize(dRdWt, matRows, matCols, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatGetSize dRdWt")

        write(*,20) "# MATRIX: dRdWt global size =", &
                    matRows, "x", matCols

        ! Gets the matrix type as a string from the matrix object.

        ! MatGetType - Gets the matrix type as a string
        !              from the matrix object.
        ! Synopsis
        !
        ! #include "petscmat.h" 
        ! call MatGetType(Mat mat,MatType *type,PetscErrorCode ierr)
        !
        ! Not Collective
        !
        ! Input Parameter
        !   mat  - the matrix
        !
        ! Output Parameter
        !   name - name of matrix type
        !
        ! see .../petsc/docs/manualpages/Mat/MatGetType.html

        call MatGetType(dRdWt, matTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatGetType dRdWt")

        write(*,30) "# MATRIX: dRdWt type        =", matTypeStr

      endif

      ! Query about the ownership range.

      ! MatGetOwnershipRange - Returns the range of matrix rows owned
      !   by this processor, assuming that the matrix is laid out with
      !   the first n1 rows on the first processor, the next n2 rows on
      !   the second, etc. For certain parallel layouts this range may
      !   not be well defined.
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! call MatGetOwnershipRange(Mat mat,PetscInt *m,PetscInt* n, &
      !                           PetscErrorCode ierr)
      ! Not Collective
      !
      ! Input Parameters
      !   mat - the matrix
      !
      ! Output Parameters
      !   m - the global index of the first local row
      !   n - one more than the global index of the last local row
      !
      ! see .../petsc/docs/manualpages/Mat/MatGetOwnershipRange.html
      ! or PETSc users manual, pp.56

      if( debug ) then
        call MatGetOwnershipRange(dRdWt, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetOwnershipRange dRdWt")

        write(*,40) "# MATRIX: dRdWt Proc", PETScRank, "; #rows =", &
                    nDimW, "; ownership =", iLow, "to", iHigh-1

      endif

!
!     ******************************************************************
!     *                                                                *
!     * Create matrix dRda that is used to compute the total cost /    *
!     * constraint function sensitivity with respect to the extra      *
!     * design variables 'a' as dIda = dJda - psi^T dRda.              *
!     *                                                                *
!     * Matrix dRda has size [nDim,nDesignExtra] and can be either     *
!     * dense or sparse depending on the choosen design vaiables.      *
!     *                                                                *
!     * The global dimension is specified so that the extra design     *
!     * variables are scattered among the processors.                  *
!     * This has to be consistent with the vectors dIda and dJda.      *
!     *                                                                *
!     ******************************************************************
!
      ! Create the matrix. Depending on either this is a sequential or 
      ! parallel run,  PETSc automatically generates the apropriate
      ! matrix type over all processes in SUMB_PETSC_COMM_WORLD

      ! MatCreate - Creates a matrix where the type is determined from
      !   either a call to MatSetType() or from the options database
      !   with a call to MatSetFromOptions(). The default matrix type is
      !   AIJ, using the routines MatCreateSeqAIJ() or MatCreateMPIAIJ()
      !   if you do not set a type in the options database. If you never
      !   call MatSetType() or MatSetFromOptions() it will generate an
      !   error when you try to use the matrix.
      !
      ! Synopsis
      !
      ! #include "petscmat.h"  
      ! call MatCreate(MPI_Comm comm,Mat *A,PetscErrorCode ierr)
      !
      ! Collective on MPI_Comm
      !
      ! Input Parameter
      !   comm - MPI communicator
      ! Output Parameter
      !   A    - the matrix 
      !
      ! see .../petsc/docs/manualpages/Mat/MatCreate.html
      ! or PETSc users manual, pp.51-52

      call MatCreate(SUMB_PETSC_COMM_WORLD, dRda, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatCreate dRda")

      ! R(# rows): Set the local size (PETSc determine the global size)
      ! a(# columns): Set the global size (PETSc decide the local size)
      !
      ! This is done for sigma sensitivity (volume). Any additional
      ! design variable can be appended to the root processor later on.
      ! It has to be consistent with vector dJda local size...

      ! MatSetSizes - Sets the local and global sizes, and checks to
      !               determine compatibility
      ! Synopsis
      !
      ! #include "petscmat.h"  
      ! call MatSetSizes(Mat A, PetscInt m, PetscInt n, &
      !                  PetscInt M, PetscInt N, PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   A - the matrix
      !   m - number of local rows (or PETSC_DECIDE)
      !   n - number of local columns (or PETSC_DECIDE)
      !   M - number of global rows (or PETSC_DETERMINE)
      !   N - number of global columns (or PETSC_DETERMINE)
      !
      ! Notes
      ! m (n) and M (N) cannot be both PETSC_DECIDE If one processor
      !   calls this with M (N) of PETSC_DECIDE then all processors
      !   must, otherwise the program will hang.
      !
      ! If PETSC_DECIDE is not used for the arguments 'm' and 'n', then
      !   the user must ensure that they are chosen to be compatible
      !   with the vectors. To do this, one first considers the
      !   matrix-vector product 'y = A x'. The 'm' that is used in the
      !   above routine must match the local size used in the vector
      !   creation routine VecCreateMPI() for 'y'. Likewise, the 'n'
      !   used must match that used as the local size in VecCreateMPI()
      !   for 'x'. 
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetSizes.html
      ! or PETSc users manual, pp.51-52

      call MatSetSizes(dRda, nDimW, PETSC_DECIDE, &
                       PETSC_DETERMINE, nDesignExtra, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetSizes dRda")

      ! MatSetType - Builds matrix object for a particular matrix type
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! call MatSetType(Mat mat, MatType matype, PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat    - the matrix object
      !   matype - matrix type
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetType.html

      call MatSetType(dRda,MATMPIDENSE,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetFromOptions dRda")

      ! Set column major order for the matrix dRda.
#ifdef USE_PETSC_3
      call MatSetOption(dRda, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dRda")
#else
      call MatSetOption(dRda, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dRda")
#endif
      ! Extract info from the global matrix (only processor 0 does it).

      if( PETScRank==0 .and. debug ) then

        ! Get the global number of rows and columns.

        call MatGetSize(dRda, matRows, matCols, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatGetSize dRda")

        write(*,20) "# MATRIX: dRda global size =", &
                    matRows, " x ", matCols

        ! Gets the matrix type as a string from the matrix object.

        call MatGetType(dRda, matTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatGetType dRda")

        write(*,30) "# MATRIX: dRda type        =", matTypeStr

      endif

      ! Query about the ownership range.

      if( debug ) then
        call MatGetOwnershipRange(dRda, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetOwnershipRange dRda")

        write(*,40) "# MATRIX: dRda Proc", PETScRank, "; #rows =", &
                    nDimW, "; ownership =", iLow, "to", iHigh-1
      endif
!
!     ******************************************************************
!     *                                                                *
!     * Create matrix dRdx that is used to compute the total cost /    *
!     * constraint function sensitivity with respect to the spatial    *
!     * design variables 'x' as dIdx = dJdx - psi^T dRdx.              *
!     *                                                                *
!     * Matrix dRdx has size [nDimW,nDimX] and is generally            *
!     * sparse for the coordinate design variables.                    *
!     *                                                                *
!     * The local dimensions are specified so that the spatial         *
!     * coordinates x (a) are placed in the local processor. This has  *
!     * to be consistent with the vectors dIdx and dJdx.               *
!     *                                                                *
!     ******************************************************************
!
      allocate( nnzDiagonal(nDimW), nnzOffDiag(nDimW) )

      nnzDiagonal = nzDiagonalX * 3!*nTimeIntervalsSpectral
      nnzOffDiag  = nzOffDiag   * 3*nTimeIntervalsSpectral

      ! Create the matrix dRdx.

      call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
                           nDimW, nDimX,                     &
                           PETSC_DETERMINE, PETSC_DETERMINE, &
                           nzDiagonalX, nnzDiagonal,     &
                           nzOffDiag, nnzOffDiag,            &
                           dRdx, PETScIerr)

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
                     "Could not create matrix dRdx of local size", nDimW
        call terminate("createPETScMat", errorMessage)
      endif

      deallocate( nnzDiagonal, nnzOffDiag )

      ! Create the matrix. Depending on either this is a sequential or
      ! parallel run,  PETSc automatically generates the apropriate
      ! matrix type over all processes in SUMB_PETSC_COMM_WORLD
      ! MatCreateMPIAIJ
      ! Creates a sparse parallel matrix in AIJ format (the default 
      ! parallel PETSc format). For good matrix assembly performance the 
      ! user should preallocate the matrix storage by setting the parameters 
      ! d_nz (or d_nnz) and o_nz (or o_nnz). By setting these parameters
      ! accurately, performance can be increased by more than a factor of 50.
      !      
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! PetscErrorCode PETSCMAT_DLLEXPORT MatCreateMPIAIJ(MPI_Comm comm,
      ! PetscInt m,PetscInt n,PetscInt M,PetscInt N,PetscInt d_nz,const 
      ! PetscInt d_nnz[],PetscInt o_nz,const PetscInt o_nnz[],Mat *A)
      
      ! Collective on MPI_Comm
      
      ! Input Parameters
      !	comm 	- MPI communicator
      !	m 	- number of local rows (or PETSC_DECIDE to have calculated 
      !   if M is given) This value should be the same as the local size 
      !   used in creating the y vector for the matrix-vector product y = Ax.
      !	n 	- This value should be the same as the local size used in
      ! creating the x vector for the matrix-vector product y = Ax.
      !	 (or PETSC_DECIDE to have calculated if N is given) For square matrices
      !  n is almost always m.
      !	M 	- number of global rows (or PETSC_DETERMINE to have calculated 
      ! if m is given)
      !	N 	- number of global columns (or PETSC_DETERMINE to have 
      ! calculated if n is given)
      !	d_nz 	- number of nonzeros per row in DIAGONAL portion of local 
      ! submatrix (same value is used for all local rows)
      !	d_nnz 	- array containing the number of nonzeros in the various rows 
      ! of the DIAGONAL portion of the local submatrix (possibly different for
      ! each row) or PETSC_NULL, if d_nz is used to specify the nonzero 
      ! structure. The size of this array is equal to the number of local 
      ! rows, i.e 'm'. You must leave room for the diagonal entry even if it is
      ! zero.
      !	o_nz 	- number of nonzeros per row in the OFF-DIAGONAL portion of 
      ! local submatrix (same value is used for all local rows).
      !	o_nnz 	- array containing the number of nonzeros in the various rows
      ! of the OFF-DIAGONAL portion of the local submatrix (possibly different
      ! for each row) or PETSC_NULL, if o_nz is used to specify the nonzero 
      ! structure. The size of this array is equal to the number of local rows,
      ! i.e 'm'.
 
      ! see .../petsc/docs/manualpages/Mat/MatCreateMPIAIJ.html
      ! or PETSc users manual, pp.??
      
!###      allocate( nnzDiagonal(nDim), nnzOffDiag(nDim) )
      
      !the size 375 is chosen to cover the possibility that all nodes in 
      !the stencil affect the residual. Also allows for consistant block
      ! placement
!###      nnzDiagonal(:) = nw**3*3!375 for euler ! ### for 2nd order stencil ###
!###      nnzOffDiag(:)  = nw**3*3!375
      
!###      call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD, nDim, nDimx, PETSC_DETERMINE,&
!###           PETSC_DETERMINE, nw**3*3, nnzDiagonal,nw**3*3, nnzOffDiag, dRdx,&
!###	   PETScIerr)

!###      if( PETScIerr/=0 ) &
!###        call terminate("createPETScMat", "Error in MatCreate dRdx")
    
!###      deallocate( nnzDiagonal, nnzOffDiag)

      ! R(# rows): Set the local size (PETSc determine the global size)
      ! x(# columns): Set the local size (PETSc determine the global size)
      !
      ! This is done for x sensitivity (volume). Any additional
      ! design variable can be appended to the root processor later on.
      ! It has to be consistent with vector dJdx local size...

!!$      ! MatCreate - Creates a matrix where the type is determined from
!!$      !   either a call to MatSetType() or from the options database
!!$      !   with a call to MatSetFromOptions(). The default matrix type
!!$      !   is AIJ, using the routines MatCreateSeqAIJ() or
!!$      !   MatCreateMPIAIJ() if you do not set a type in the options
!!$      !   database. If you never call MatSetType() or
!!$      !   MatSetFromOptions() it will generate an error when you try
!!$      !   to use the matrix.
!!$      !
!!$      ! Synopsis
!!$      !
!!$      ! #include "petscmat.h"
!!$      ! call MatCreate(MPI_Comm comm,Mat *A,PetscErrorCode ierr)
!!$      !
!!$      ! Collective on MPI_Comm
!!$      !
!!$      ! Input Parameter
!!$      !   comm - MPI communicator
!!$      ! Output Parameter
!!$      !   A    - the matrix
!!$      !
!!$      ! see .../petsc/docs/manualpages/Mat/MatCreate.html
!!$      ! or PETSc users manual, pp.51-52
!!$
!!$      call MatCreate(SUMB_PETSC_COMM_WORLD, dRdx, PETScIerr)
!!$
!!$        if( PETScIerr/=0 ) &
!!$          call terminate("createPETScMat","Error in MatCreate dRdx")
!!$
!!$      ! R(#rows): Set the local size (PETSc determine the global size)
!!$      ! a(#cols): Set the local size (PETSc determine the global size)
!!$      !
!!$      ! This is done for coordinates x (3*volume). Any additional
!!$      ! design variable can later be appended to the root processor.
!!$      ! It has to be consistent with vector dJda local size...
!!$
!!$      ! MatSetSizes - Sets the local and global sizes, and checks to
!!$      !               determine compatibility
!!$      ! Synopsis
!!$      !
!!$      ! #include "petscmat.h"
!!$      ! call MatSetSizes(Mat A, PetscInt m, PetscInt n, &
!!$      !                  PetscInt M, PetscInt N, PetscErrorCode ierr)
!!$      !
!!$      ! Collective on Mat
!!$      !
!!$      ! Input Parameters
!!$      !   A - the matrix
!!$      !   m - number of local rows (or PETSC_DECIDE)
!!$      !   n - number of local columns (or PETSC_DECIDE)
!!$      !   M - number of global rows (or PETSC_DETERMINE)
!!$      !   N - number of global columns (or PETSC_DETERMINE)
!!$      !
!!$      ! Notes
!!$      ! m (n) and M (N) cannot be both PETSC_DECIDE If one processor
!!$      !   calls this with M (N) of PETSC_DECIDE then all processors
!!$      !   must, otherwise the program will hang.
!!$      !
!!$      ! If PETSC_DECIDE is not used for the arguments 'm' and 'n', then
!!$      !   the user must ensure that they are chosen to be compatible
!!$      !   with the vectors. To do this, one first considers the
!!$      !   matrix-vector product 'y = A x'. The 'm' that is used in the
!!$      !   above routine must match the local size used in the vector
!!$      !   creation routine VecCreateMPI() for 'y'. Likewise, the 'n'
!!$      !   used must match that used as the local size in VecCreateMPI()
!!$      !   for 'x'.
!!$      !
!!$      ! see .../petsc/docs/manualpages/Mat/MatSetSizes.html
!!$      ! or PETSc users manual, pp.51-52
!!$
!!$      call MatSetSizes(dRdx, nDim, nDimX, &
!!$                       PETSC_DETERMINE, PETSC_DETERMINE, PETScIerr)
!!$
!!$      if( PETScIerr/=0 ) &
!!$        call terminate("createPETScMat", "Error in MatSetSizes dRdx")

      ! MatSetFromOptions - Creates a matrix where the type is
      !   determined from the options database. Generates a parallel
      !   MPI matrix if the communicator has more than one processor.
      !   The default matrix type is AIJ, using the routines
      !   MatCreateSeqAIJ() and MatCreateMPIAIJ() if you do not select
      !   a type in the options database.
      !
      ! Synopsis
      !
      ! #include "petscmat.h"
      ! call MatSetFromOptions(Mat B, PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameter
      !   A - the matrix
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetFromOptions.html

      call MatSetFromOptions(dRdx, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetFromOptions dRdx")

      ! Set column major order for the matrix dRdx.
#ifdef USE_PETSC_3
      call MatSetOption(dRdx, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)
      call MatSetOption(dRdx,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE,PETScIErr)
      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dRdx")
#else
      call MatSetOption(dRdx, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetOption dRdx")
#endif
      ! Extract info from the global matrix (only processor 0).

      if( PETScRank==0 .and. debug ) then

        ! Get the global number of rows and columns.

        call MatGetSize(dRdx, matRows, matCols, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetSize dRdx")

        write(*,20) "# MATRIX: dRdx global size =", &
                    matRows, " x ", matCols

        ! Gets the matrix type as a string from the matrix object.

        call MatGetType(dRdx, matTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetType dRdx")

        write(*,30) "# MATRIX: dRdx type        =", matTypeStr

      endif

      ! Query about the ownership range.

      if( debug ) then
        call MatGetOwnershipRange(dRdx, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetOwnershipRange dRdx")

        write(*,40) "# MATRIX: dRdx Proc", PETScRank,   &
                    "; #rows =", nDimW,                 &
                    "; ownership =", iLow, "to", iHigh-1
      endif


if(debug)then
!***********************
! Create dRdxFD
!***********************

!
!     ******************************************************************
!     *                                                                *
!     * Create matrix dRdx that is used to compute the total cost /    *
!     * constraint function sensitivity with respect to the spatial    *
!     * design variables 'x' as dIdx = dJdx - psi^T dRdx.              *
!     *                                                                *
!     * Matrix dRdx has size [nDimW,nDimX] and is generally            *
!     * sparse for the coordinate design variables.                    *
!     *                                                                *
!     * The local dimensions are specified so that the spatial         *
!     * coordinates x (a) are placed in the local processor. This has  *
!     * to be consistent with the vectors dIdx and dJdx.               *
!     *                                                                *
!     ******************************************************************
!
      allocate( nnzDiagonal(nDimW), nnzOffDiag(nDimW) )

      nnzDiagonal = nzDiagonalX * 3*nTimeIntervalsSpectral
      nnzOffDiag  = nzOffDiag   * 3*nTimeIntervalsSpectral

      ! Create the matrix dRdx.

      call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
                           nDimW, nDimX,                     &
                           PETSC_DETERMINE, PETSC_DETERMINE, &
                           nzDiagonalX, nnzDiagonal,         &
                           nzOffDiag, nnzOffDiag,            &
                           dRdxFD, PETScIerr)

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
                     "Could not create matrix dRdxFD of local size", nDimW
        call terminate("createPETScMat", errorMessage)
      endif

      deallocate( nnzDiagonal, nnzOffDiag )

      ! MatSetFromOptions - Creates a matrix where the type is
      !   determined from the options database. Generates a parallel
      !   MPI matrix if the communicator has more than one processor.
      !   The default matrix type is AIJ, using the routines
      !   MatCreateSeqAIJ() and MatCreateMPIAIJ() if you do not select
      !   a type in the options database.
      !
      ! Synopsis
      !
      ! #include "petscmat.h"
      ! call MatSetFromOptions(Mat B, PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameter
      !   A - the matrix
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetFromOptions.html

      call MatSetFromOptions(dRdxFD, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetFromOptions dRdxFD")

      ! Set column major order for the matrix dRdxFD.
#ifdef USE_PETSC_3
      call MatSetOption(dRdxFD, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dRdxfd")
#else
      call MatSetOption(dRdxFD, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetOption dRdxFD")
#endif
      ! Extract info from the global matrix (only processor 0).

      if( PETScRank==0 .and. debug ) then

        ! Get the global number of rows and columns.

        call MatGetSize(dRdxFD, matRows, matCols, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetSize dRdxFD")

        write(*,20) "# MATRIX: dRdxFD global size =", &
                    matRows, " x ", matCols

        ! Gets the matrix type as a string from the matrix object.

        call MatGetType(dRdxFD, matTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetType dRdxFD")

        write(*,30) "# MATRIX: dRdx type        =", matTypeStr

      endif

      ! Query about the ownership range.

      if( debug ) then
        call MatGetOwnershipRange(dRdxFD, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetOwnershipRange dRdxFD")

        write(*,40) "# MATRIX: dRdxFD Proc", PETScRank,   &
                    "; #rows =", nDimW,                 &
                    "; ownership =", iLow, "to", iHigh-1
      endif

!********************
! end create dRdxFD
!********************
   endif

if(.false.)then
!
!     ******************************************************************
!     *                                                                *
!     * Create matrix dSdx that is used to compute the coupling        *
!     * with the structures                                            *
!     *                                                                *
!     *                                                                *
!     * Matrix dRdx has size [nDimS,nDimX] and is generally            *
!     * sparse for the coordinate design variables.                    *
!     *                                                                *
!     * The local dimensions are specified so that the spatial         *
!     * coordinates x (a) are placed in the local processor. This has  *
!     * to be consistent with the vectors dIdx and dJdx.               *
!     *                                                                *
!     ******************************************************************
!
      allocate( nnzDiagonal(nDimS), nnzOffDiag(nDimS) )

      nnzDiagonal = nzDiagonalX * 3
      nnzOffDiag  = nzOffDiag   * 3

      ! Create the matrix dSdx.

      call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
                           nDimS, nDimX,                     &
                           PETSC_DETERMINE, PETSC_DETERMINE, &
                           nzDiagonalX, nnzDiagonal,         &
                           nzOffDiag, nnzOffDiag,            &
                           dSdx, PETScIerr)

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
                     "Could not create matrix dSdx of local size", nDimX
        call terminate("createPETScMat", errorMessage)
      endif

      deallocate( nnzDiagonal, nnzOffDiag )

      ! S(# rows): Set the local size (PETSc determine the global size)
      ! x(# columns): Set the local size (PETSc determine the global size)
      !
      ! This is done for x sensitivity (volume). Any additional
      ! design variable can be appended to the root processor later on.
      ! It has to be consistent with vector dJdx local size...

      ! MatSetFromOptions - Creates a matrix where the type is
      !   determined from the options database. Generates a parallel
      !   MPI matrix if the communicator has more than one processor.
      !   The default matrix type is AIJ, using the routines
      !   MatCreateSeqAIJ() and MatCreateMPIAIJ() if you do not select
      !   a type in the options database.
      !
      ! Synopsis
      !
      ! #include "petscmat.h"
      ! call MatSetFromOptions(Mat B, PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameter
      !   A - the matrix
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetFromOptions.html

      call MatSetFromOptions(dSdx, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetFromOptions dSdx")

      ! Set column major order for the matrix dSdx.
#ifdef USE_PETSC_3
      call MatSetOption(dSdx, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dsdx")
#else
      call MatSetOption(dSdx, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetOption dSdx")
#endif

   endif
      ! Extract info from the global matrix (only processor 0).

      if( PETScRank==0 .and. debug ) then

        ! Get the global number of rows and columns.

        call MatGetSize(dSdx, matRows, matCols, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetSize dSdx")

        write(*,20) "# MATRIX: dSdx global size =", &
                    matRows, " x ", matCols

        ! Gets the matrix type as a string from the matrix object.

        call MatGetType(dSdx, matTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetType dSdx")

        write(*,30) "# MATRIX: dSdx type        =", matTypeStr

      endif

      ! Query about the ownership range.

      if( debug ) then
        call MatGetOwnershipRange(dSdx, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetOwnershipRange dSdx")

        write(*,40) "# MATRIX: dSdx Proc", PETScRank,   &
                    "; #rows =", nDimW,                 &
                    "; ownership =", iLow, "to", iHigh-1
      endif

if( debug) then
!**************
!
!dsdxfd for debug
!
!**************

!
!     ******************************************************************
!     *                                                                *
!     * Create matrix dSdxfd that is used to compute the coupling        *
!     * with the structures                                            *
!     *                                                                *
!     *                                                                *
!     * Matrix dRdx has size [nDimS,nDimX] and is generally            *
!     * sparse for the coordinate design variables.                    *
!     *                                                                *
!     * The local dimensions are specified so that the spatial         *
!     * coordinates x (a) are placed in the local processor. This has  *
!     * to be consistent with the vectors dIdx and dJdx.               *
!     *                                                                *
!     ******************************************************************
!
      allocate( nnzDiagonal(nDimS), nnzOffDiag(nDimS) )

      nnzDiagonal = nzDiagonalX * 3
      nnzOffDiag  = nzOffDiag   * 3

      ! Create the matrix dSdxfd.

      call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
                           nDimS, nDimX,                     &
                           PETSC_DETERMINE, PETSC_DETERMINE, &
                           nzDiagonalX, nnzDiagonal,         &
                           nzOffDiag, nnzOffDiag,            &
                           dSdxfd2, PETScIerr)

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
                     "Could not create matrix dSdxfd of local size", nDimX
        call terminate("createPETScMat", errorMessage)
      endif

      deallocate( nnzDiagonal, nnzOffDiag )

      ! S(# rows): Set the local size (PETSc determine the global size)
      ! x(# columns): Set the local size (PETSc determine the global size)
      !
      ! This is done for x sensitivity (volume). Any additional
      ! design variable can be appended to the root processor later on.
      ! It has to be consistent with vector dJdx local size...

      ! MatSetFromOptions - Creates a matrix where the type is
      !   determined from the options database. Generates a parallel
      !   MPI matrix if the communicator has more than one processor.
      !   The default matrix type is AIJ, using the routines
      !   MatCreateSeqAIJ() and MatCreateMPIAIJ() if you do not select
      !   a type in the options database.
      !
      ! Synopsis
      !
      ! #include "petscmat.h"
      ! call MatSetFromOptions(Mat B, PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameter
      !   A - the matrix
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetFromOptions.html

      call MatSetFromOptions(dSdxfd2, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetFromOptions dSdxfd")

      ! Set column major order for the matrix dSdx.
#ifdef USE_PETSC_3
      call MatSetOption(dSdxfd2, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dsdxfd2")
#else
      call MatSetOption(dSdxfd2, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetOption dSdxfd")
#endif
      ! Extract info from the global matrix (only processor 0).

      if( PETScRank==0 .and. debug ) then

        ! Get the global number of rows and columns.

        call MatGetSize(dSdxfd2, matRows, matCols, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetSize dSdxfd")

        write(*,20) "# MATRIX: dSdxfd global size =", &
                    matRows, " x ", matCols

        ! Gets the matrix type as a string from the matrix object.

        call MatGetType(dSdxfd2, matTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetType dSdxfd")

        write(*,30) "# MATRIX: dSdxfd type        =", matTypeStr

      endif

      ! Query about the ownership range.

      if( debug ) then
        call MatGetOwnershipRange(dSdxfd2, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetOwnershipRange dSdxfD")

        write(*,40) "# MATRIX: dSdxFD Proc", PETScRank,   &
                    "; #rows =", nDimS,                 &
                    "; ownership =", iLow, "to", iHigh-1
      endif

!*************
!
! end dsdxfd
!
!************
   endif

if(.false.)then
!
!     ******************************************************************
!     *                                                                *
!     * Create matrix dSdx that is used to compute the coupling        *
!     * with the structures                                            *
!     *                                                                *
!     *                                                                *
!     * Matrix dRdx has size [nDimS,nDimX] and is generally            *
!     * sparse for the coordinate design variables.                    *
!     *                                                                *
!     * The local dimensions are specified so that the spatial         *
!     * coordinates x (a) are placed in the local processor. This has  *
!     * to be consistent with the vectors dIdx and dJdx.               *
!     *                                                                *
!     ******************************************************************
!
      allocate( nnzDiagonal(nDimS), nnzOffDiag(nDimS) )

      nnzDiagonal = nzDiagonalX * 3
      nnzOffDiag  = nzOffDiag   * 3

      ! Create the matrix dSdx.

      call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
                           nDimS, nDimW,                     &
                           PETSC_DETERMINE, PETSC_DETERMINE, &
                           nzDiagonalX, nnzDiagonal,         &
                           nzOffDiag, nnzOffDiag,            &
                           dSdw, PETScIerr)

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
                     "Could not create matrix dSdw of local size", nDimW
        call terminate("createPETScMat", errorMessage)
      endif

      deallocate( nnzDiagonal, nnzOffDiag )

      ! S(# rows): Set the local size (PETSc determine the global size)
      ! x(# columns): Set the local size (PETSc determine the global size)
      !
      ! This is done for x sensitivity (volume). Any additional
      ! design variable can be appended to the root processor later on.
      ! It has to be consistent with vector dJdx local size...

      ! MatSetFromOptions - Creates a matrix where the type is
      !   determined from the options database. Generates a parallel
      !   MPI matrix if the communicator has more than one processor.
      !   The default matrix type is AIJ, using the routines
      !   MatCreateSeqAIJ() and MatCreateMPIAIJ() if you do not select
      !   a type in the options database.
      !
      ! Synopsis
      !
      ! #include "petscmat.h"
      ! call MatSetFromOptions(Mat B, PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameter
      !   A - the matrix
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetFromOptions.html

      call MatSetFromOptions(dSdw, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetFromOptions dSdw")

      ! Set column major order for the matrix dSdw.
#ifdef USE_PETSC_3
      call MatSetOption(dSdw, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dsdW")
#else
      call MatSetOption(dSdw, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetOption dSdw")
#endif

   endif
      ! Extract info from the global matrix (only processor 0).

      if( PETScRank==0 .and. debug ) then

        ! Get the global number of rows and columns.

        call MatGetSize(dSdw, matRows, matCols, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetSize dSdw")

        write(*,20) "# MATRIX: dSdw global size =", &
                    matRows, " x ", matCols

        ! Gets the matrix type as a string from the matrix object.

        call MatGetType(dSdw, matTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetType dSdw")

        write(*,30) "# MATRIX: dSdw type        =", matTypeStr

      endif

      ! Query about the ownership range.

      if( debug ) then
        call MatGetOwnershipRange(dSdw, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetOwnershipRange dSdw")

        write(*,40) "# MATRIX: dSdw Proc", PETScRank,   &
                    "; #rows =", nDimW,                 &
                    "; ownership =", iLow, "to", iHigh-1
      endif
!
!     ******************************************************************
!     dCdw

!
!     ******************************************************************
!     *                                                                *
!     * Create matrix dCdW that is used to compute the RHS of the      *
!     * ADjoint for the time spectral case                             *
!     *                                                                *
!     * Matrix dCdw has size [nTimeIntervalsSpectral,nDimW] and is     *
!     * generally sparse.                                              *
!     *                                                                *
!     *                                                                *
!     ******************************************************************
!
      ! Create the matrix. Depending on either this is a sequential or 
      ! parallel run,  PETSc automatically generates the apropriate
      ! matrix type over all processes in SUMB_PETSC_COMM_WORLD

      ! MatCreate - Creates a matrix where the type is determined from
      !   either a call to MatSetType() or from the options database
      !   with a call to MatSetFromOptions(). The default matrix type is
      !   AIJ, using the routines MatCreateSeqAIJ() or MatCreateMPIAIJ()
      !   if you do not set a type in the options database. If you never
      !   call MatSetType() or MatSetFromOptions() it will generate an
      !   error when you try to use the matrix.
      !
      ! Synopsis
      !
      ! #include "petscmat.h"  
      ! call MatCreate(MPI_Comm comm,Mat *A,PetscErrorCode ierr)
      !
      ! Collective on MPI_Comm
      !
      ! Input Parameter
      !   comm - MPI communicator
      ! Output Parameter
      !   A    - the matrix 
      !
      ! see .../petsc/docs/manualpages/Mat/MatCreate.html
      ! or PETSc users manual, pp.51-52

      call MatCreate(SUMB_PETSC_COMM_WORLD, dCdw, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatCreate dCdw")

      ! C(# rows): Set the global size (PETSc determine the global size)
      ! w(# columns): Set the local size (PETSc decide the local size)
      !

      ! MatSetSizes - Sets the local and global sizes, and checks to
      !               determine compatibility
      ! Synopsis
      !
      ! #include "petscmat.h"  
      ! call MatSetSizes(Mat A, PetscInt m, PetscInt n, &
      !                  PetscInt M, PetscInt N, PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   A - the matrix
      !   m - number of local rows (or PETSC_DECIDE)
      !   n - number of local columns (or PETSC_DECIDE)
      !   M - number of global rows (or PETSC_DETERMINE)
      !   N - number of global columns (or PETSC_DETERMINE)
      !
      ! Notes
      ! m (n) and M (N) cannot be both PETSC_DECIDE If one processor
      !   calls this with M (N) of PETSC_DECIDE then all processors
      !   must, otherwise the program will hang.
      !
      ! If PETSC_DECIDE is not used for the arguments 'm' and 'n', then
      !   the user must ensure that they are chosen to be compatible
      !   with the vectors. To do this, one first considers the
      !   matrix-vector product 'y = A x'. The 'm' that is used in the
      !   above routine must match the local size used in the vector
      !   creation routine VecCreateMPI() for 'y'. Likewise, the 'n'
      !   used must match that used as the local size in VecCreateMPI()
      !   for 'x'. 
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetSizes.html
      ! or PETSc users manual, pp.51-52

      call MatSetSizes(dCdw, PETSC_DECIDE,nDimW, &
                       nTimeIntervalsSpectral,PETSC_DETERMINE, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetSizes dCdw")

      ! MatSetType - Builds matrix object for a particular matrix type
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! call MatSetType(Mat mat, MatType matype, PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat    - the matrix object
      !   matype - matrix type
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetType.html

      call MatSetType(dCdw,MATMPIAIJ,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetFromOptions dCdw")

      ! Set column major order for the matrix dRda.
#ifdef USE_PETSC_3
      call MatSetOption(dCdw, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dcdw")
#else
      call MatSetOption(dCdw, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dCdw")
#endif
      ! Extract info from the global matrix (only processor 0 does it).

      if( PETScRank==0 .and. debug ) then

        ! Get the global number of rows and columns.

        call MatGetSize(dCdw, matRows, matCols, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatGetSize dCdw")

        write(*,20) "# MATRIX: dCdw global size =", &
                    matRows, " x ", matCols

        ! Gets the matrix type as a string from the matrix object.

        call MatGetType(dCdw, matTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatGetType dCdw")

        write(*,30) "# MATRIX: dCdw type        =", matTypeStr

      endif

      ! Query about the ownership range.

      if( debug ) then
        call MatGetOwnershipRange(dCdw, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetOwnershipRange dCdw")

        write(*,40) "# MATRIX: dCdw Proc", PETScRank, "; #rows =", &
                    nTimeIntervalsSpectral, "; ownership =", iLow, "to", iHigh-1
      endif

!     ******************************************************************
!     dCdx

!
!     ******************************************************************
!     *                                                                *
!     * Create matrix dCdx that is used to compute the Gradient        *
!     * partial derivative for the ADjoint.                            *
!     *                                                                *
!     * Matrix dCdx has size [nTimeIntervalsSpectral,nDimX] and is     *
!     * generally sparse.                                              *
!     *                                                                *
!     *                                                                *
!     ******************************************************************
!
      ! Create the matrix. Depending on either this is a sequential or 
      ! parallel run,  PETSc automatically generates the apropriate
      ! matrix type over all processes in SUMB_PETSC_COMM_WORLD

      ! MatCreate - Creates a matrix where the type is determined from
      !   either a call to MatSetType() or from the options database
      !   with a call to MatSetFromOptions(). The default matrix type is
      !   AIJ, using the routines MatCreateSeqAIJ() or MatCreateMPIAIJ()
      !   if you do not set a type in the options database. If you never
      !   call MatSetType() or MatSetFromOptions() it will generate an
      !   error when you try to use the matrix.
      !
      ! Synopsis
      !
      ! #include "petscmat.h"  
      ! call MatCreate(MPI_Comm comm,Mat *A,PetscErrorCode ierr)
      !
      ! Collective on MPI_Comm
      !
      ! Input Parameter
      !   comm - MPI communicator
      ! Output Parameter
      !   A    - the matrix 
      !
      ! see .../petsc/docs/manualpages/Mat/MatCreate.html
      ! or PETSc users manual, pp.51-52

      call MatCreate(SUMB_PETSC_COMM_WORLD, dCdx, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatCreate dCdx")

      ! C(# rows): Set the global size (PETSc determine the local size)
      ! x(# columns): Set the local size (PETSc decide the global size)
      !
      ! This is done for sigma sensitivity (volume). Any additional
      ! design variable can be appended to the root processor later on.
      ! It has to be consistent with vector dJda local size...

      ! MatSetSizes - Sets the local and global sizes, and checks to
      !               determine compatibility
      ! Synopsis
      !
      ! #include "petscmat.h"  
      ! call MatSetSizes(Mat A, PetscInt m, PetscInt n, &
      !                  PetscInt M, PetscInt N, PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   A - the matrix
      !   m - number of local rows (or PETSC_DECIDE)
      !   n - number of local columns (or PETSC_DECIDE)
      !   M - number of global rows (or PETSC_DETERMINE)
      !   N - number of global columns (or PETSC_DETERMINE)
      !
      ! Notes
      ! m (n) and M (N) cannot be both PETSC_DECIDE If one processor
      !   calls this with M (N) of PETSC_DECIDE then all processors
      !   must, otherwise the program will hang.
      !
      ! If PETSC_DECIDE is not used for the arguments 'm' and 'n', then
      !   the user must ensure that they are chosen to be compatible
      !   with the vectors. To do this, one first considers the
      !   matrix-vector product 'y = A x'. The 'm' that is used in the
      !   above routine must match the local size used in the vector
      !   creation routine VecCreateMPI() for 'y'. Likewise, the 'n'
      !   used must match that used as the local size in VecCreateMPI()
      !   for 'x'. 
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetSizes.html
      ! or PETSc users manual, pp.51-52

      call MatSetSizes(dCdx, PETSC_DECIDE,nDimX, &
                       nTimeIntervalsSpectral,PETSC_DETERMINE, PETScIerr)
      !call MatSetSizes(dCdx,PETSC_DETERMINE, nTimeIntervalsSpectral , &
      !                 nDimX, PETSC_DECIDE, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetSizes dCdx")


      ! MatSetType - Builds matrix object for a particular matrix type
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! call MatSetType(Mat mat, MatType matype, PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat    - the matrix object
      !   matype - matrix type
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetType.html

      call MatSetType(dCdx,MATMPIAIJ,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetFromOptions dCdw")


      ! Set column major order for the matrix dRda.
#ifdef USE_PETSC_3
      call MatSetOption(dCdx, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dcdx")
#else
      call MatSetOption(dCdx, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dCdx")
#endif
      ! Extract info from the global matrix (only processor 0 does it).

      if( PETScRank==0 .and. debug ) then

        ! Get the global number of rows and columns.

        call MatGetSize(dCdx, matRows, matCols, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatGetSize dCdx")

        write(*,20) "# MATRIX: dCdx global size =", &
                    matRows, " x ", matCols

        ! Gets the matrix type as a string from the matrix object.

        call MatGetType(dCdx, matTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatGetType dCdx")

        write(*,30) "# MATRIX: dCdx type        =", matTypeStr

      endif

      ! Query about the ownership range.

      if( debug ) then
        call MatGetOwnershipRange(dCdx, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetOwnershipRange dCdx")

        write(*,40) "# MATRIX: dCdx Proc", PETScRank, "; #rows =", &
                    nTimeIntervalsSpectral, "; ownership =", iLow, "to", iHigh-1
      endif

!     ******************************************************************
!     dCda

!
!     ******************************************************************
!     *                                                                *
!     * Create matrix dCda that is used to compute the partial         *
!     * derivative of the Extra variabled for the the time spectral    *
!     * case                                                           *
!     *                                                                *
!     * Matrix dRda has size [nTimeIntervals,nDesignExtra] and is      *
!     * generally dense.                                               *
!     *                                                                *

!     *                                                                *
!     ******************************************************************
!
      ! Create the matrix. Depending on either this is a sequential or 
      ! parallel run,  PETSc automatically generates the apropriate
      ! matrix type over all processes in SUMB_PETSC_COMM_WORLD

      ! MatCreate - Creates a matrix where the type is determined from
      !   either a call to MatSetType() or from the options database
      !   with a call to MatSetFromOptions(). The default matrix type is
      !   AIJ, using the routines MatCreateSeqAIJ() or MatCreateMPIAIJ()
      !   if you do not set a type in the options database. If you never
      !   call MatSetType() or MatSetFromOptions() it will generate an
      !   error when you try to use the matrix.
      !
      ! Synopsis
      !
      ! #include "petscmat.h"  
      ! call MatCreate(MPI_Comm comm,Mat *A,PetscErrorCode ierr)
      !
      ! Collective on MPI_Comm
      !
      ! Input Parameter
      !   comm - MPI communicator
      ! Output Parameter
      !   A    - the matrix 
      !
      ! see .../petsc/docs/manualpages/Mat/MatCreate.html
      ! or PETSc users manual, pp.51-52

      call MatCreate(SUMB_PETSC_COMM_WORLD, dCda, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatCreate dCda")

      ! C(# rows): Set the global size (PETSc determine the local size)
      ! a(# columns): Set the global size (PETSc decide the local size)
  

      ! MatSetSizes - Sets the local and global sizes, and checks to
      !               determine compatibility
      ! Synopsis
      !
      ! #include "petscmat.h"  
      ! call MatSetSizes(Mat A, PetscInt m, PetscInt n, &
      !                  PetscInt M, PetscInt N, PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   A - the matrix
      !   m - number of local rows (or PETSC_DECIDE)
      !   n - number of local columns (or PETSC_DECIDE)
      !   M - number of global rows (or PETSC_DETERMINE)
      !   N - number of global columns (or PETSC_DETERMINE)
      !
      ! Notes
      ! m (n) and M (N) cannot be both PETSC_DECIDE If one processor
      !   calls this with M (N) of PETSC_DECIDE then all processors
      !   must, otherwise the program will hang.
      !
      ! If PETSC_DECIDE is not used for the arguments 'm' and 'n', then
      !   the user must ensure that they are chosen to be compatible
      !   with the vectors. To do this, one first considers the
      !   matrix-vector product 'y = A x'. The 'm' that is used in the
      !   above routine must match the local size used in the vector
      !   creation routine VecCreateMPI() for 'y'. Likewise, the 'n'
      !   used must match that used as the local size in VecCreateMPI()
      !   for 'x'. 
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetSizes.html
      ! or PETSc users manual, pp.51-52

      call MatSetSizes(dCda, PETSC_DECIDE,PETSC_DECIDE, &
                       nTimeIntervalsSpectral,nDesignExtra, PETScIerr)
     ! call MatSetSizes(dCda,PETSC_DETERMINE, nTimeIntervalsSpectral, &
     !                  PETSC_DETERMINE, nDesignExtra, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetSizes dCda")

      ! MatSetType - Builds matrix object for a particular matrix type
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! call MatSetType(Mat mat, MatType matype, PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat    - the matrix object
      !   matype - matrix type
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetType.html

      call MatSetType(dCda,MATMPIDENSE,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetType dCda")

      ! Set column major order for the matrix dRda.
#ifdef USE_PETSC_3
      call MatSetOption(dCda, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dcda")
#else
      call MatSetOption(dCda, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dCda")
#endif
      ! Extract info from the global matrix (only processor 0 does it).

      if( PETScRank==0 .and. debug ) then

        ! Get the global number of rows and columns.

        call MatGetSize(dCda, matRows, matCols, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatGetSize dCda")

        write(*,20) "# MATRIX: dCda global size =", &
                    matRows, " x ", matCols

        ! Gets the matrix type as a string from the matrix object.

        call MatGetType(dCda, matTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatGetType dCda")

        write(*,30) "# MATRIX: dCda type        =", matTypeStr

      endif

      ! Query about the ownership range.

      if( debug ) then
        call MatGetOwnershipRange(dCda, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetOwnershipRange dCda")

        write(*,40) "# MATRIX: dCda Proc", PETScRank, "; #rows =", &
             nTimeIntervalsSpectral, "; ownership =", iLow, "to", iHigh-1
     endif
       

     ! Synchronize the processors.

     call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

!     ******************************************************************
      ! Output formats.

   10 format(a,1x,i2)                            ! block size
   20 format(a,1x,i6,1x,a,1x,i6)                 ! global size
   30 format(a,1x,a)                             ! type
   40 format(a,1x,i3,a,1x,i6,a,1x,i6,1x,a,1x,i6) ! ownership
   99 format(a,1x,i6)                            ! error

#endif

      end subroutine createPETScMat
