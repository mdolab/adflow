!
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
      use flowVarRefState ! 
      implicit none
!
!     Local variables.
!
      ! nDim   - local dimension (processor owned)
      ! iLow   - first component owned by the local process
      ! iHigh  - one more than last component owned by the local process

      integer       :: nn, iLow, iHigh
      integer       :: nDimW, nDimX
      integer       :: matBlockSize, matRows, matCols
      character(15) :: matTypeStr

      integer   :: nzDiagonalW, nzDiagonalX, nzOffDiag
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

      ! Define matrix dRdW local size, taking into account the total
      ! number of Cells owned by the processor and the number of 
      ! equations.

      nDimW = nw * nCellsLocal     

      ! Define matrix dRdx local size (number of columns) for the
      ! spatial derivatives.

      nDimX = 3 * nNodesLocal

      ! Number of non-zero blocks per residual row in dRdW
      ! >>> This depends on the stencil being used R=R(W)
      !     1st order stencil ->  7 cells
      !     2nd order stencil -> 13 cells

      ! Stencil of W
      ! 1  - center cell
      ! 6  - 1st level cells along directions i,j,k
      ! 6  - 2nd level cells along directions i,j,k

      nzDiagonalW = 13 ! 1 + 6 + 6  check!!!

      ! Stencil of X
      ! 1  - center node
      ! 6  - 1st level nodes along directions i,j,k
      ! 6  - 2nd level nodes along directions i,j,k
      ! 12 - 1st level nodes along diagonals (i,j),(i,k),(j,k) 

      nzDiagonalX = 25 ! 1 + 6 + 6 + 12 Check

      ! Average number of off processor contributions per Cell
      ! (average number of donor cells that come from other processor)

      nzOffDiag  = 3
!
!     ******************************************************************
!     *                                                                *
!     * Create matrix dRdW that define the adjoint linear system of    *
!     * equations, dRdW^T psi = dJdW. Matrix dRdW has size [nDim,nDim] *
!     * but is very sparse because of the computational stencil R=R(W).*
!     *                                                                *
!     ******************************************************************
!
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
  
        allocate( nnzDiagonal(nCellsLocal), nnzOffDiag(nCellsLocal) )

        nnzDiagonal = nzDiagonalW
        nnzOffDiag  = nzOffDiag

        call MatCreateMPIBAIJ(PETSC_COMM_WORLD, nw,             &
                              nDimW, nDimW,                     &
                              PETSC_DETERMINE, PETSC_DETERMINE, &
                              nzDiagonalW, nnzDiagonal,         &
                              nzOffDiag, nnzOffDiag,            &
                              dRdW, PETScIerr)

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

        call MatCreateMPIAIJ(PETSC_COMM_WORLD,                 &
                             nDimW, nDimW,                     &
                             PETSC_DETERMINE, PETSC_DETERMINE, &
                             nzDiagonalW, nnzDiagonal,         &
                             nzOffDiag, nnzOffDiag,            &
                             dRdW, PETScIerr)

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

      call MatSetOption(dRdW, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dRdW")

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

        call MatCreateMPIAIJ(PETSC_COMM_WORLD,                 &
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

      call MatSetOption(dRdWFD, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dRdWFD")

      !*****************************************
      ! end of create dRdWFD
      !*****************************************


      ! Extract info from the global matrix (only processor 0 does it).

      if( PETScRank==0 .and. debug ) then

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

        call MatGetBlockSize(dRdW, matBlockSize, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetBlockSize dRdW")

        write(*,10) "# MATRIX: dRdW block size  =", matBlockSize

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

        call MatGetSize(dRdW, matRows, matCols, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatGetSize dRdW")

        write(*,20) "# MATRIX: dRdW global size =", &
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

        call MatGetType(dRdW, matTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatGetType dRdW")

        write(*,30) "# MATRIX: dRdW type        =", matTypeStr

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
        call MatGetOwnershipRange(dRdW, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetOwnershipRange dRdW")

        write(*,40) "# MATRIX: dRdW Proc", PETScRank, "; #rows =", &
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
      ! matrix type over all processes in PETSC_COMM_WORLD

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

      call MatCreate(PETSC_COMM_WORLD, dRda, PETScIerr)

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

      call MatSetOption(dRda, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dRda")

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

      nnzDiagonal = nzDiagonalX * 3
      nnzOffDiag  = nzOffDiag   * 3

      ! Create the matrix dRdx.

      call MatCreateMPIAIJ(PETSC_COMM_WORLD,                 &
                           nDimW, nDimX,                     &
                           PETSC_DETERMINE, PETSC_DETERMINE, &
                           nzDiagonalX, nnzDiagonal,         &
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
      ! matrix type over all processes in PETSC_COMM_WORLD
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
      
!###      call MatCreateMPIAIJ(PETSC_COMM_WORLD, nDim, nDimx, PETSC_DETERMINE,&
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
!!$      call MatCreate(PETSC_COMM_WORLD, dRdx, PETScIerr)
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

      call MatSetOption(dRdx, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetOption dRdx")

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

      nnzDiagonal = nzDiagonalX * 3
      nnzOffDiag  = nzOffDiag   * 3

      ! Create the matrix dRdx.

      call MatCreateMPIAIJ(PETSC_COMM_WORLD,                 &
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

      call MatSetOption(dRdxFD, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", &
                       "Error in MatSetOption dRdxFD")

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

!
!     ******************************************************************
!
      ! Output formats.

   10 format(a,1x,i2)                            ! block size
   20 format(a,1x,i6,1x,a,1x,i6)                 ! global size
   30 format(a,1x,a)                             ! type
   40 format(a,1x,i3,a,1x,i6,a,1x,i6,1x,a,1x,i6) ! ownership
   99 format(a,1x,i6)                            ! error

#endif

      end subroutine createPETScMat
