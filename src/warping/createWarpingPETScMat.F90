
!     ******************************************************************
!     *                                                                *
!     * File:          createWarpingPETScMat.F90                       *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 05-21-2009                                      *
!     * Last modified: 05-21-2009                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine createWarpingPETScMat
!
!     ******************************************************************
!     *                                                                *
!     * Create the global matrix dRdW of the adjoint system of         *
!     * equations as a PETSc matrix object. Also, create the auxiliar  *
!     * matrices dRda and dRdx used to compute the gradients.          *
!     *                                                                *
!     ******************************************************************
!
      !use ADjointPETSc
      use warpingPETSc
      use ADjointVars     ! nCellsLocal,nNodesLocal, nDesignExtra
      use communication   ! myID, nProc
      use flowVarRefState ! 
      use mdData          !mdNSurfNodes,mdNSurfNodesCompact
      use mdDataLocal     !mdSurfGlobalIndLocal
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

      integer   ::  nzDiagonalXs, nzOffDiag,nzDiagonal
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
   

      ! Define matrix dXvdXs local size (number of columns) for the mesh
      ! volume coordinates.
      nDimW = nw*nCellsLocal

      nDimX = 3 * nNodesLocal
      ! Define matrix dXvdXs global size (number of Rows) for the
      ! surface coordinates.
      nDimS = 3 * mdNSurfNodesCompact

      ! Number of non-zero blocks per residual row in dRdW
      ! >>> This depends on the stencil being used R=R(W)
      !     1st order stencil ->  7 cells
      !     2nd order stencil -> 13 cells

      ! Stencil of Xs
      ! 1  - corner node
      
	
      nzDiagonalXs = 1
      nzDiagonal = mdNSurfNodesCompact/(nproc*2)!13 ! 1 + 6 + 6  check!!!
      print *,'nzdiagonal', nzDiagonal,mdNSurfNodesCompact

      ! Average number of off processor contributions per Cell
      ! (average number of donor cells that come from other processor)

      nzOffDiag  = 1!3
!
!     ******************************************************************
!     *                                                                *
!     * Create matrix dXvdXs to store teh derivative of the volume     *
!     * coordinates wrt the surface coodinates.                        *
!     * Matrix dXvdXs has size [nDimXv,nDimS].                          *
!     *                                                                *
!     ******************************************************************
!
      ! Create the matrix dXvdXs.

      ! sparse parallel matrix in AIJ format
      !                 General case...

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

        nzDiagonalXs = nzDiagonalXs * 3
        nzOffDiag   = nzOffDiag   * 3

        allocate( nnzDiagonal(nDimX), nnzOffDiag(nDimX) )

        nnzDiagonal = nzDiagonalXs
        nnzOffDiag  = nzOffDiag

	!print *,'petscnull',PETSC_NULL

        !call MatCreateMPIAIJ(PETSC_COMM_WORLD,                 &
        !                     nDimX, nDimS,                     &
        !                     PETSC_DETERMINE, PETSC_DETERMINE, &
        !                     nzDiagonalW, nnzDiagonal,         &
        !                     nzOffDiag, nnzOffDiag,            &
        !                     dXvdXs, PETScIerr)
        !call MatCreateMPIAIJ(PETSC_COMM_WORLD,                 &
        !                     nDimX,PETSC_DECIDE,        &
	!	             PETSC_DETERMINE, nDimS,                     &
        !                     0,PETSC_NULL,         &
        !                     0, PETSC_NULL,            &
        !                     dXvdXs, PETScIerr)
        call MatCreateMPIAIJ(PETSC_COMM_WORLD,                 &
                             nDimX,PETSC_DECIDE,               &
		             PETSC_DETERMINE, nDimS,           &
                             nzDiagonalXs, nnzDiagonal,         &
                             nzOffDiag, nnzOffDiag,            &
                             dXvdXs, PETScIerr)

      deallocate( nnzDiagonal, nnzOffDiag )

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
                     "Could not create matrix dXvdXs of local size", nDimX
        call terminate("createWarpingPETScMat", errorMessage)
      endif

      ! Set the matrix dXvdXs options.

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

      call MatSetOption(dXvdXs, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dXvdXs")

      if(debug)then
      !******************************************
      !Create dXvdXsFD for debugging
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

        nzDiagonalXs = nzDiagonalXs * 3
        nzOffDiag   = nzOffDiag   * 3

        allocate( nnzDiagonal(nDimX), nnzOffDiag(nDimX) )

        nnzDiagonal = nzDiagonalXs
        nnzOffDiag  = nzOffDiag
        !print *,'petscnull 2',PETSC_NULL	
        !call MatCreateMPIAIJ(PETSC_COMM_WORLD,                 &
        !                     nDimX, nDimS,                     &
        !                     PETSC_DETERMINE, PETSC_DETERMINE, &
        !                     nzDiagonalW, nnzDiagonal,         &
        !                     nzOffDiag, nnzOffDiag,            &
        !                     dXvdXsFD, PETScIerr)
	call MatCreateMPIAIJ(PETSC_COMM_WORLD,                 &
                             nDimX,PETSC_DECIDE,        &
		             PETSC_DETERMINE, nDimS,                     &
                             nzDiagonalXs, nnzDiagonal,         &
                             nzOffDiag, nnzOffDiag,            &
                             dXvdXsFD, PETScIerr)

      deallocate( nnzDiagonal, nnzOffDiag )

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
                     "Could not create matrix dXvdXsFD of local size", nDimX
        call terminate("createPETScMat", errorMessage)
      endif

      ! Set the matrix dXvdXsFD options.

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

      call MatSetOption(dXvdXsFD, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createwarpingPETScMat", "Error in MatSetOption dXvdXsFD")

      !*****************************************
      ! end of create dXvdXsFD
      !*****************************************
      !******************************************
      !Create dXvdXsPara for debugging
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

        nzDiagonalXs = nzDiagonalXs * 3
        nzOffDiag   = nzOffDiag   * 3

        allocate( nnzDiagonal(nDimX), nnzOffDiag(nDimX) )

        nnzDiagonal = nzDiagonalXs
        nnzOffDiag  = nzOffDiag
        !print *,'petscnull 2',PETSC_NULL	
   	call MatCreateMPIAIJ(PETSC_COMM_WORLD,                 &
                             nDimX,PETSC_DECIDE,        &
		             PETSC_DETERMINE, nDimS,                     &
                             nzDiagonalXs, nnzDiagonal,         &
                             nzOffDiag, nnzOffDiag,            &
                             dXvdXsPara, PETScIerr)

      deallocate( nnzDiagonal, nnzOffDiag )

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
                     "Could not create matrix dXvdXsFD of local size", nDimX
        call terminate("createPETScMat", errorMessage)
      endif

      ! Set the matrix dXvdXsPara options.

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

      call MatSetOption(dXvdXsPara, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createwarpingPETScMat", "Error in MatSetOption dXvdXsFD")

      !*****************************************
      ! end of create dXvdXsPara
      !*****************************************
   endif

!
!     ******************************************************************
!     *                                                                *
!     * Create matrix dRdXs to store the derivative of the residuals   *
!     * coordinates wrt the surface coodinates.                        *
!     * Matrix dRdXs has size [nDimw,nDimS].                         *
!     *                                                                *
!     ******************************************************************
!
      ! Create the matrix dRdXs.

      ! sparse parallel matrix in AIJ format
      !                 General case...

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

        nzDiagonal = nzDiagonal * 3!nw
        nzOffDiag   = nzOffDiag   * 3
        
        allocate( nnzDiagonal(nDimW), nnzOffDiag(nDimW) )

        nnzDiagonal = nzDiagonal
        nnzOffDiag  = nzOffDiag

	!print *,'petscnull',PETSC_NULL

        !call MatCreateMPIAIJ(PETSC_COMM_WORLD,                 &
        !                     nDimX, nDimS,                     &
        !                     PETSC_DETERMINE, PETSC_DETERMINE, &
        !                     nzDiagonalW, nnzDiagonal,         &
        !                     nzOffDiag, nnzOffDiag,            &
        !                     dXvdXs, PETScIerr)
        !call MatCreateMPIAIJ(PETSC_COMM_WORLD,                 &
        !                     nDimX,PETSC_DECIDE,        &
	!	             PETSC_DETERMINE, nDimS,                     &
        !                     0,PETSC_NULL,         &
        !                     0, PETSC_NULL,            &
        !                     dXvdXs, PETScIerr)
	print *,'creating dxvdxs',nDimW,PETSC_DECIDE,               &
		             PETSC_DETERMINE, nDimS
        call MatCreateMPIAIJ(PETSC_COMM_WORLD,                 &
                             nDimW,PETSC_DECIDE,               &
		             PETSC_DETERMINE, nDimS,           &
                             nzDiagonal, nnzDiagonal,         &
                             nzOffDiag, nnzOffDiag,            &
                             dRdXs, PETScIerr)


      deallocate( nnzDiagonal, nnzOffDiag )

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
                     "Could not create matrix dRdXs of local size", nDimW
        call terminate("createWarpingPETScMat", errorMessage)
      endif

      ! Set the matrix dRdXs options.

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

      call MatSetOption(dRdXs, MAT_COLUMN_ORIENTED, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScMat", "Error in MatSetOption dRdXs")


   
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

        call MatGetBlockSize(dXvdXs, matBlockSize, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetBlockSize dXvdXs")

        write(*,10) "# MATRIX: dXvdXs block size  =", matBlockSize

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

        call MatGetSize(dXvdXs, matRows, matCols, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatGetSize dXvdXs")

        write(*,20) "# MATRIX: dXvdXs global size =", &
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

        call MatGetType(dXvdXs, matTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatGetType dXvdXs")

        write(*,30) "# MATRIX: dXvdXs type        =", matTypeStr

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
        call MatGetOwnershipRange(dXvdXs, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScMat", &
                         "Error in MatGetOwnershipRange dXvdXs")

        write(*,40) "# MATRIX: dXvdXs Proc", PETScRank, "; #rows =", &
                    nDimX, "; ownership =", iLow, "to", iHigh-1

      endif

!********************************************************
      ! Output formats.

   10 format(a,1x,i2)                            ! block size
   20 format(a,1x,i6,1x,a,1x,i6)                 ! global size
   30 format(a,1x,a)                             ! type
   40 format(a,1x,i3,a,1x,i6,a,1x,i6,1x,a,1x,i6) ! ownership
   99 format(a,1x,i6)                            ! error

#endif

      end subroutine createWarpingPETScMat
