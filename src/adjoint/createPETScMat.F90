
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

  call getForceSize(nDimS)
  nDimS = nDimS * 3

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

  nzDiagonalWPC = 7+ 1*(nTimeIntervalsSpectral-1)
  !should be addition not multiplication but this will require filtering 
  !of some of the TS results for non-zero terms. 
  !i.e.:7 +nTimeIntervalsSpectral! 1 + 6  check!!!

  ! Stencil of X
  ! 1  - center node
  ! 6  - 1st level nodes along directions i,j,k
  ! 6  - 2nd level nodes along directions i,j,k
  ! 12 - 1st level nodes along diagonals (i,j),(i,k),(j,k) 

  nzDiagonalX = 32+8*(nTimeIntervalsSpectral-1)!25+4*nTimeIntervalsSpectral! 1 + 6 + 6 + 12 Check

  ! Average number of off processor contributions per Cell
  ! (average number of donor cells that come from other processor)

  !      nzOffDiag  = 12+1*(nTimeIntervalsSpectral-1)

  nzOffDiag  = 3+5*(nTimeIntervalsSpectral-1)

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create matrix dRdW that define the adjoint linear system of    *
  !     * equations, dRdW^T psi = dJdW. Matrix dRdW has size [nDim,nDim] *
  !     * but is very sparse because of the computational stencil R=R(W).*
  !     *                                                                *
  !     ******************************************************************
  !
  if( nw <= 7 ) then

     PETScBlockMatrix = .true.

     allocate( nnzDiagonal(nCellsLocal*nTimeIntervalsSpectral),&
          nnzOffDiag(nCellsLocal*nTimeIntervalsSpectral) )

     nnzDiagonal = nzDiagonalW
     nnzOffDiag  = nzOffDiag

     call drdwPreAllocation(nnzDiagonal,nnzOffDiag,nCellsLocal)

     call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
          nDimW, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dRdWT, PETScIerr)


  else

     PETScBlockMatrix = .false.
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



#ifdef USE_PETSC_3
  call MatSetOption(dRdWt, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)
  if( PETScIerr/=0 ) &
       call terminate("createPETScMat", "Error in MatSetOption dRdW")
#else
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

     if( nw <= 7 ) then

        PETScBlockMatrix = .true.
        allocate( nnzDiagonal(nCellsLocal*nTimeIntervalsSpectral),&
             nnzOffDiag(nCellsLocal*nTimeIntervalsSpectral) )

        nnzDiagonal = nzDiagonalWPC
        nnzOffDiag  = nzOffDiag

        call drdwPCPreAllocation(nnzDiagonal,nnzOffDiag,nCellsLocal)

        call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
             nDimW, nDimW,                     &
             PETSC_DETERMINE, PETSC_DETERMINE, &
             0, nnzDiagonal,         &
             0, nnzOffDiag,            &
             dRdWPreT, PETScIerr)

     else

        PETScBlockMatrix = .false.
        nzDiagonalW = nzDiagonalW * nw
        nzOffDiag   = nzOffDiag   * nw

        allocate( nnzDiagonal(nDimW), nnzOffDiag(nDimW) )

        nnzDiagonal = nzDiagonalWPC
        nnzOffDiag  = nzOffDiag
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

#ifdef USE_PETSC_3
     call MatSetOption(dRdWPret, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)
     if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatSetOption dRdW")
#else
     call MatSetOption(dRdWPret, MAT_COLUMN_ORIENTED, PETScIerr)
     if( PETScIerr/=0 ) &
          call terminate("createPETScMat", "Error in MatSetOption dRdWPre")
#endif
  end if ! Approx PC

  !  if(Debug) then

  !      !******************************************
  !      !Create dRdWFD for debugging
  !      !*****************************************
  !      nzDiagonalW = nzDiagonalW * nw
  !      nzOffDiag   = nzOffDiag   * nw

  !      allocate( nnzDiagonal(nDimW), nnzOffDiag(nDimW) )

  !      nnzDiagonal = nzDiagonalW
  !      nnzOffDiag  = nzOffDiag

  !      call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
  !           nDimW, nDimW,                     &
  !           PETSC_DETERMINE, PETSC_DETERMINE, &
  !           nzDiagonalW, nnzDiagonal,         &
  !           nzOffDiag, nnzOffDiag,            &
  !           dRdWFD, PETScIerr)

  !      deallocate( nnzDiagonal, nnzOffDiag )

  !      if( PETScIerr/=0 ) then
  !         write(errorMessage,99) &
  !              "Could not create matrix dRdWFD of local size", nDimW
  !         call terminate("createPETScMat", errorMessage)
  !      endif

  ! #ifdef USE_PETSC_3
  !      call MatSetOption(dRdWFD, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

  !      if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", "Error in MatSetOption dRdW")
  ! #else
  !      call MatSetOption(dRdWFD, MAT_COLUMN_ORIENTED, PETScIerr)

  !      if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", "Error in MatSetOption dRdWFD")
  ! #endif
  !      !*****************************************
  !      ! end of create dRdWFD
  !      !*****************************************
  !   end if

  !   ! Extract info from the global matrix (only processor 0 does it).

  !   if( PETScRank==0 .and. debug ) then

  !      call MatGetBlockSize(dRdWt, matBlockSize, PETScIerr)

  !      if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", &
  !           "Error in MatGetBlockSize dRdW")

  !      write(*,10) "# MATRIX: dRdWt block size  =", matBlockSize


  !      call MatGetSize(dRdWt, matRows, matCols, PETScIerr)

  !      if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", "Error in MatGetSize dRdWt")

  !      write(*,20) "# MATRIX: dRdWt global size =", &
  !           matRows, "x", matCols

  !      call MatGetType(dRdWt, matTypeStr, PETScIerr)

  !      if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", "Error in MatGetType dRdWt")

  !      write(*,30) "# MATRIX: dRdWt type        =", matTypeStr

  !   endif

  !   if( debug ) then
  !      call MatGetOwnershipRange(dRdWt, iLow, iHigh, PETScIerr)

  !      if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", &
  !           "Error in MatGetOwnershipRange dRdWt")

  !      write(*,40) "# MATRIX: dRdWt Proc", PETScRank, "; #rows =", &
  !           nDimW, "; ownership =", iLow, "to", iHigh-1

  !   endif


  call MatCreate(SUMB_PETSC_COMM_WORLD, dRda, PETScIerr)

  if( PETScIerr/=0 ) &
       call terminate("createPETScMat", "Error in MatCreate dRda")

  call MatSetSizes(dRda, nDimW, PETSC_DECIDE, &
       PETSC_DETERMINE, nDesignExtra, PETScIerr)

  if( PETScIerr/=0 ) &
       call terminate("createPETScMat", "Error in MatSetSizes dRda")


  call MatSetType(dRda,MATMPIDENSE,PETScIerr)

  if( PETScIerr/=0 ) &
       call terminate("createPETScMat", &
       "Error in MatSetFromOptions dRda")

  ! Set column major order for the matrix dRda.
#ifdef USE_PETSC_3
  call MatSetOption(dRda, MAT_ROW_ORIENTED,PETSC_TRUE, PETScIerr)

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

  !   if( debug ) then
  !      call MatGetOwnershipRange(dRda, iLow, iHigh, PETScIerr)

  !      if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", &
  !           "Error in MatGetOwnershipRange dRda")

  !      write(*,40) "# MATRIX: dRda Proc", PETScRank, "; #rows =", &
  !           nDimW, "; ownership =", iLow, "to", iHigh-1
  !   endif
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
  nnzOffDiag  = 12!nzOffDiag   * 3*nTimeIntervalsSpectral

  ! Create the matrix dRdx.

  call drdxPreAllocation(nnzDiagonal,nnzOffDiag,nDimW)
  
  call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
       nDimW, nDimX,                     &
       PETSC_DETERMINE, PETSC_DETERMINE, &
       0, nnzDiagonal,     &
       0, nnzOffDiag,            &
       dRdx, PETScIerr)

  if( PETScIerr/=0 ) then
     write(errorMessage,99) &
          "Could not create matrix dRdx of local size", nDimW
     call terminate("createPETScMat", errorMessage)
  endif

  deallocate( nnzDiagonal, nnzOffDiag )

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

  !  if( PETScRank==0 .and. debug ) then

  !      ! Get the global number of rows and columns.

  !      call MatGetSize(dRdx, matRows, matCols, PETScIerr)

  !      if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", &
  !           "Error in MatGetSize dRdx")

  !      write(*,20) "# MATRIX: dRdx global size =", &
  !           matRows, " x ", matCols

  !      ! Gets the matrix type as a string from the matrix object.

  !      call MatGetType(dRdx, matTypeStr, PETScIerr)

  !      if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", &
  !           "Error in MatGetType dRdx")

  !      write(*,30) "# MATRIX: dRdx type        =", matTypeStr

  !   endif

  ! Query about the ownership range.

  !   if( debug ) then
  !      call MatGetOwnershipRange(dRdx, iLow, iHigh, PETScIerr)

  !      if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", &
  !           "Error in MatGetOwnershipRange dRdx")

  !      write(*,40) "# MATRIX: dRdx Proc", PETScRank,   &
  !           "; #rows =", nDimW,                 &
  !           "; ownership =", iLow, "to", iHigh-1
  !   endif


  !   if(debug)then
  !      !***********************
  !      ! Create dRdxFD
  !      !***********************

  !      !
  !      !     ******************************************************************
  !      !     *                                                                *
  !      !     * Create matrix dRdx that is used to compute the total cost /    *
  !      !     * constraint function sensitivity with respect to the spatial    *
  !      !     * design variables 'x' as dIdx = dJdx - psi^T dRdx.              *
  !      !     *                                                                *
  !      !     * Matrix dRdx has size [nDimW,nDimX] and is generally            *
  !      !     * sparse for the coordinate design variables.                    *
  !      !     *                                                                *
  !      !     * The local dimensions are specified so that the spatial         *
  !      !     * coordinates x (a) are placed in the local processor. This has  *
  !      !     * to be consistent with the vectors dIdx and dJdx.               *
  !      !     *                                                                *
  !      !     ******************************************************************
  !      !
  !      allocate( nnzDiagonal(nDimW), nnzOffDiag(nDimW) )

  !      nnzDiagonal = nzDiagonalX * 3*nTimeIntervalsSpectral
  !      nnzOffDiag  = nzOffDiag   * 3*nTimeIntervalsSpectral

  !      ! Create the matrix dRdx.

  !      call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
  !           nDimW, nDimX,                     &
  !           PETSC_DETERMINE, PETSC_DETERMINE, &
  !           nzDiagonalX, nnzDiagonal,         &
  !           nzOffDiag, nnzOffDiag,            &
  !           dRdxFD, PETScIerr)

  !      if( PETScIerr/=0 ) then
  !         write(errorMessage,99) &
  !              "Could not create matrix dRdxFD of local size", nDimW
  !         call terminate("createPETScMat", errorMessage)
  !      endif

  !      deallocate( nnzDiagonal, nnzOffDiag )


  !      call MatSetFromOptions(dRdxFD, PETScIerr)

  !      if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", &
  !           "Error in MatSetFromOptions dRdxFD")

  !      ! Set column major order for the matrix dRdxFD.
  ! #ifdef USE_PETSC_3
  !      call MatSetOption(dRdxFD, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

  !      if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", "Error in MatSetOption dRdxfd")
  ! #else
  !      call MatSetOption(dRdxFD, MAT_COLUMN_ORIENTED, PETScIerr)

  !      if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", &
  !           "Error in MatSetOption dRdxFD")
  ! #endif
  !      ! Extract info from the global matrix (only processor 0).

  !      if( PETScRank==0 .and. debug ) then

  !         ! Get the global number of rows and columns.

  !         call MatGetSize(dRdxFD, matRows, matCols, PETScIerr)

  !         if( PETScIerr/=0 ) &
  !              call terminate("createPETScMat", &
  !              "Error in MatGetSize dRdxFD")

  !         write(*,20) "# MATRIX: dRdxFD global size =", &
  !              matRows, " x ", matCols

  !         ! Gets the matrix type as a string from the matrix object.

  !         call MatGetType(dRdxFD, matTypeStr, PETScIerr)

  !         if( PETScIerr/=0 ) &
  !              call terminate("createPETScMat", &
  !              "Error in MatGetType dRdxFD")

  !         write(*,30) "# MATRIX: dRdx type        =", matTypeStr

  !      endif

  !      ! Query about the ownership range.

  !      if( debug ) then
  !         call MatGetOwnershipRange(dRdxFD, iLow, iHigh, PETScIerr)

  !         if( PETScIerr/=0 ) &
  !              call terminate("createPETScMat", &
  !              "Error in MatGetOwnershipRange dRdxFD")

  !         write(*,40) "# MATRIX: dRdxFD Proc", PETScRank,   &
  !              "; #rows =", nDimW,                 &
  !              "; ownership =", iLow, "to", iHigh-1
  !      endif

  !      !********************
  !      ! end create dRdxFD
  !      !********************
  !   endif


  ! !
  ! !     ******************************************************************
  ! !     *                                                                *
  ! !     * Create matrix dSdx that is used to compute the coupling        *
  ! !     * with the structures                                            *
  ! !     *                                                                *
  ! !     *                                                                *
  ! !     * Matrix dRdx has size [nDimS,nDimX] and is generally            *
  ! !     * sparse for the coordinate design variables.                    *
  ! !     *                                                                *
  ! !     * The local dimensions are specified so that the spatial         *
  ! !     * coordinates x (a) are placed in the local processor. This has  *
  ! !     * to be consistent with the vectors dIdx and dJdx.               *
  ! !     *                                                                *
  ! !     ******************************************************************
  ! !
  !       allocate( nnzDiagonal(nDimS), nnzOffDiag(nDimS) )

  !       nnzDiagonal = nzDiagonalX * 3
  !       nnzOffDiag  = nzOffDiag   * 3

  !       ! Create the matrix dSdx.

  !       call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
  !                            nDimS, nDimX,                     &
  !                            PETSC_DETERMINE, PETSC_DETERMINE, &
  !                            nzDiagonalX, nnzDiagonal,         &
  !                            nzOffDiag, nnzOffDiag,            &
  !                            dSdx, PETScIerr)

  !       if( PETScIerr/=0 ) then
  !         write(errorMessage,99) &
  !                      "Could not create matrix dSdx of local size", nDimX
  !         call terminate("createPETScMat", errorMessage)
  !       endif

  !       deallocate( nnzDiagonal, nnzOffDiag )


  !       call MatSetFromOptions(dSdx, PETScIerr)

  !       if( PETScIerr/=0 ) &
  !         call terminate("createPETScMat", &
  !                        "Error in MatSetFromOptions dSdx")

  !       ! Set column major order for the matrix dSdx.
  ! #ifdef USE_PETSC_3
  !       call MatSetOption(dSdx, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

  !       if( PETScIerr/=0 ) &
  !         call terminate("createPETScMat", "Error in MatSetOption dsdx")
  ! #else
  !       call MatSetOption(dSdx, MAT_COLUMN_ORIENTED, PETScIerr)

  !       if( PETScIerr/=0 ) &
  !         call terminate("createPETScMat", &
  !                        "Error in MatSetOption dSdx")
  ! #endif


  !       ! Extract info from the global matrix (only processor 0).

  !       if( PETScRank==0 .and. debug ) then

  !         ! Get the global number of rows and columns.

  !         call MatGetSize(dSdx, matRows, matCols, PETScIerr)

  !         if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", &
  !                          "Error in MatGetSize dSdx")

  !         write(*,20) "# MATRIX: dSdx global size =", &
  !                     matRows, " x ", matCols

  !         ! Gets the matrix type as a string from the matrix object.

  !         call MatGetType(dSdx, matTypeStr, PETScIerr)

  !         if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", &
  !                          "Error in MatGetType dSdx")

  !         write(*,30) "# MATRIX: dSdx type        =", matTypeStr

  !       endif

  !       ! Query about the ownership range.

  !       if( debug ) then
  !         call MatGetOwnershipRange(dSdx, iLow, iHigh, PETScIerr)

  !         if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", &
  !                          "Error in MatGetOwnershipRange dSdx")

  !         write(*,40) "# MATRIX: dSdx Proc", PETScRank,   &
  !                     "; #rows =", nDimW,                 &
  !                     "; ownership =", iLow, "to", iHigh-1
  !       endif

  ! if( debug) then
  ! !**************
  ! !
  ! !dsdxfd for debug
  ! !
  ! !**************

  ! !
  ! !     ******************************************************************
  ! !     *                                                                *
  ! !     * Create matrix dSdxfd that is used to compute the coupling        *
  ! !     * with the structures                                            *
  ! !     *                                                                *
  ! !     *                                                                *
  ! !     * Matrix dRdx has size [nDimS,nDimX] and is generally            *
  ! !     * sparse for the coordinate design variables.                    *
  ! !     *                                                                *
  ! !     * The local dimensions are specified so that the spatial         *
  ! !     * coordinates x (a) are placed in the local processor. This has  *
  ! !     * to be consistent with the vectors dIdx and dJdx.               *
  ! !     *                                                                *
  ! !     ******************************************************************
  ! !
  !       allocate( nnzDiagonal(nDimS), nnzOffDiag(nDimS) )

  !       nnzDiagonal = nzDiagonalX * 3
  !       nnzOffDiag  = nzOffDiag   * 3

  !       ! Create the matrix dSdxfd.

  !       call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
  !                            nDimS, nDimX,                     &
  !                            PETSC_DETERMINE, PETSC_DETERMINE, &
  !                            nzDiagonalX, nnzDiagonal,         &
  !                            nzOffDiag, nnzOffDiag,            &
  !                            dSdxfd2, PETScIerr)

  !       if( PETScIerr/=0 ) then
  !         write(errorMessage,99) &
  !                      "Could not create matrix dSdxfd of local size", nDimX
  !         call terminate("createPETScMat", errorMessage)
  !       endif

  !       deallocate( nnzDiagonal, nnzOffDiag )

  !       call MatSetFromOptions(dSdxfd2, PETScIerr)

  !       if( PETScIerr/=0 ) &
  !         call terminate("createPETScMat", &
  !                        "Error in MatSetFromOptions dSdxfd")

  !       ! Set column major order for the matrix dSdx.
  ! #ifdef USE_PETSC_3
  !       call MatSetOption(dSdxfd2, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

  !       if( PETScIerr/=0 ) &
  !         call terminate("createPETScMat", "Error in MatSetOption dsdxfd2")
  ! #else
  !       call MatSetOption(dSdxfd2, MAT_COLUMN_ORIENTED, PETScIerr)

  !       if( PETScIerr/=0 ) &
  !         call terminate("createPETScMat", &
  !                        "Error in MatSetOption dSdxfd")
  ! #endif
  !       ! Extract info from the global matrix (only processor 0).

  !       if( PETScRank==0 .and. debug ) then

  !         ! Get the global number of rows and columns.

  !         call MatGetSize(dSdxfd2, matRows, matCols, PETScIerr)

  !         if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", &
  !                          "Error in MatGetSize dSdxfd")

  !         write(*,20) "# MATRIX: dSdxfd global size =", &
  !                     matRows, " x ", matCols

  !         ! Gets the matrix type as a string from the matrix object.

  !         call MatGetType(dSdxfd2, matTypeStr, PETScIerr)

  !         if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", &
  !                          "Error in MatGetType dSdxfd")

  !         write(*,30) "# MATRIX: dSdxfd type        =", matTypeStr

  !       endif

  !       ! Query about the ownership range.

  !       if( debug ) then
  !         call MatGetOwnershipRange(dSdxfd2, iLow, iHigh, PETScIerr)

  !         if( PETScIerr/=0 ) &
  !           call terminate("createPETScMat", &
  !                          "Error in MatGetOwnershipRange dSdxfD")

  !         write(*,40) "# MATRIX: dSdxFD Proc", PETScRank,   &
  !                     "; #rows =", nDimS,                 &
  !                     "; ownership =", iLow, "to", iHigh-1
  !       endif

  ! !*************
  ! !
  ! ! end dsdxfd
  ! !
  ! !************
  !    endif

  ! Create dSdx and dSdw

  allocate( nnzDiagonal(nDimS), nnzOffDiag(nDimS) )

  ! Each nodal force is contribed by (nominally) 4 quadrilateral cells
  ! surrounding it. The pressure on each of these cells depend on the
  ! average of the pressure of the cell above and the halo below. The
  ! 1st halo is computed by linear pressure extrapolation from the
  ! first two cells ABOVE the surface. The results in each node being
  ! affected by 8 cells. All of these cells on on-processors so there
  ! should be zero offdiag entries. For dSdx, the coordinates of each
  ! of the 4 quadrilateral affects the force, so this results in 9
  ! points spatial points affecting the force. Each coordiante has 3
  ! dimension which results in 3x3x3=27 nonzeros per row
  nnzDiagonal = 8*nw
  nnzOffDiag  = 8*nw! Make the off diagonal the same, since we
  ! don't know where the non-zeros will end up

  ! Create the matrix dSdw

  call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,&
       nDimS, nDimW,                     &
       PETSC_DETERMINE, PETSC_DETERMINE, &
       nzDiagonalX, nnzDiagonal,         &
       nzOffDiag, nnzOffDiag,            &
       dSdw, PETScIerr)

  ! Create the matrix dSdx
  nnzDiagonal = 27
  nnzOffDiag = 27
  call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,&
       nDimS, nDimS,                     &
       PETSC_DETERMINE, PETSC_DETERMINE, &
       nzDiagonalX, nnzDiagonal,         &
       nzOffDiag, nnzOffDiag,            &
       dSdx, PETScIerr)

  deallocate( nnzDiagonal, nnzOffDiag )

  ! Set column major order for the matrix dSdw.
#ifdef USE_PETSC_3
  call MatSetOption(dSdw, MAT_ROW_ORIENTED,PETSC_TRUE, PETScIerr)
  call MatSetOption(dSdx, MAT_ROW_ORIENTED,PETSC_TRUE, PETScIerr)
#else
  call MatSetOption(dSdw, MAT_ROW_ORIENTED, PETScIerr)
  call MatSetOption(dSdx, MAT_ROW_ORIENTED, PETScIerr)
#endif

  
! !       ******************************************************************
! !       dCdw

  
! !       ******************************************************************
! !       *                                                                *
! !       * Create matrix dCdW that is used to compute the RHS of the      *
! !       * ADjoint for the time spectral case                             *
! !       *                                                                *
! !       * Matrix dCdw has size [nTimeIntervalsSpectral,nDimW] and is     *
! !       * generally sparse.                                              *
! !       *                                                                *
! !       *                                                                *
! !       ******************************************************************
  

!     call MatCreate(SUMB_PETSC_COMM_WORLD, dCdw, PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", "Error in MatCreate dCdw")


!         call MatSetSizes(dCdw, PETSC_DECIDE,nDimW, &
!                          nTimeIntervalsSpectral,PETSC_DETERMINE, PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", "Error in MatSetSizes dCdw")

!         call MatSetType(dCdw,MATMPIAIJ,PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", &
!                          "Error in MatSetFromOptions dCdw")

!         ! Set column major order for the matrix dRda.
! #ifdef USE_PETSC_3
!         call MatSetOption(dCdw, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", "Error in MatSetOption dcdw")
! #else
!         call MatSetOption(dCdw, MAT_COLUMN_ORIENTED, PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", "Error in MatSetOption dCdw")
! #endif
!         ! Extract info from the global matrix (only processor 0 does it).

!         if( PETScRank==0 .and. debug ) then

!           ! Get the global number of rows and columns.

!           call MatGetSize(dCdw, matRows, matCols, PETScIerr)

!           if( PETScIerr/=0 ) &
!             call terminate("createPETScMat", "Error in MatGetSize dCdw")

!           write(*,20) "# MATRIX: dCdw global size =", &
!                       matRows, " x ", matCols

!           ! Gets the matrix type as a string from the matrix object.

!           call MatGetType(dCdw, matTypeStr, PETScIerr)

!           if( PETScIerr/=0 ) &
!             call terminate("createPETScMat", "Error in MatGetType dCdw")

!           write(*,30) "# MATRIX: dCdw type        =", matTypeStr

!         endif

!         ! Query about the ownership range.

!         if( debug ) then
!           call MatGetOwnershipRange(dCdw, iLow, iHigh, PETScIerr)

!           if( PETScIerr/=0 ) &
!             call terminate("createPETScMat", &
!                            "Error in MatGetOwnershipRange dCdw")

!           write(*,40) "# MATRIX: dCdw Proc", PETScRank, "; #rows =", &
!                       nTimeIntervalsSpectral, "; ownership =", iLow, "to", iHigh-1
!         endif

!   !     ******************************************************************
!   !     dCdx

!         call MatCreate(SUMB_PETSC_COMM_WORLD, dCdx, PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", "Error in MatCreate dCdx")


!         call MatSetSizes(dCdx, PETSC_DECIDE,nDimX, &
!                          nTimeIntervalsSpectral,PETSC_DETERMINE, PETScIerr)
!         !call MatSetSizes(dCdx,PETSC_DETERMINE, nTimeIntervalsSpectral , &
!         !                 nDimX, PETSC_DECIDE, PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", "Error in MatSetSizes dCdx")


!         call MatSetType(dCdx,MATMPIAIJ,PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", &
!                          "Error in MatSetFromOptions dCdw")


!         ! Set column major order for the matrix dRda.
! #ifdef USE_PETSC_3
!         call MatSetOption(dCdx, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", "Error in MatSetOption dcdx")
! #else
!         call MatSetOption(dCdx, MAT_COLUMN_ORIENTED, PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", "Error in MatSetOption dCdx")
! #endif
!         ! Extract info from the global matrix (only processor 0 does it).

!         if( PETScRank==0 .and. debug ) then

!           ! Get the global number of rows and columns.

!           call MatGetSize(dCdx, matRows, matCols, PETScIerr)

!           if( PETScIerr/=0 ) &
!             call terminate("createPETScMat", "Error in MatGetSize dCdx")

!           write(*,20) "# MATRIX: dCdx global size =", &
!                       matRows, " x ", matCols

!           ! Gets the matrix type as a string from the matrix object.

!           call MatGetType(dCdx, matTypeStr, PETScIerr)

!           if( PETScIerr/=0 ) &
!             call terminate("createPETScMat", "Error in MatGetType dCdx")

!           write(*,30) "# MATRIX: dCdx type        =", matTypeStr

!         endif

!         ! Query about the ownership range.

!         if( debug ) then
!           call MatGetOwnershipRange(dCdx, iLow, iHigh, PETScIerr)

!           if( PETScIerr/=0 ) &
!             call terminate("createPETScMat", &
!                            "Error in MatGetOwnershipRange dCdx")

!           write(*,40) "# MATRIX: dCdx Proc", PETScRank, "; #rows =", &
!                       nTimeIntervalsSpectral, "; ownership =", iLow, "to", iHigh-1
!         endif

!   !     ******************************************************************
!   !     dCda

!   !
!   !     ******************************************************************
!   !     *                                                                *
!   !     * Create matrix dCda that is used to compute the partial         *
!   !     * derivative of the Extra variabled for the the time spectral    *
!   !     * case                                                           *
!   !     *                                                                *
!   !     * Matrix dRda has size [nTimeIntervals,nDesignExtra] and is      *
!   !     * generally dense.                                               *
!   !     *                                                                *

!   !     *                                                                *
!   !     ******************************************************************
!   !

!         call MatCreate(SUMB_PETSC_COMM_WORLD, dCda, PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", "Error in MatCreate dCda")


!         call MatSetSizes(dCda, PETSC_DECIDE,PETSC_DECIDE, &
!                          nTimeIntervalsSpectral,nDesignExtra, PETScIerr)
!        ! call MatSetSizes(dCda,PETSC_DETERMINE, nTimeIntervalsSpectral, &
!        !                  PETSC_DETERMINE, nDesignExtra, PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", "Error in MatSetSizes dCda")

!         call MatSetType(dCda,MATMPIDENSE,PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", &
!                          "Error in MatSetType dCda")

!         ! Set column major order for the matrix dRda.
! #ifdef USE_PETSC_3
!         call MatSetOption(dCda, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", "Error in MatSetOption dcda")
! #else
!         call MatSetOption(dCda, MAT_COLUMN_ORIENTED, PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScMat", "Error in MatSetOption dCda")
! #endif
!         ! Extract info from the global matrix (only processor 0 does it).

!         if( PETScRank==0 .and. debug ) then

!           ! Get the global number of rows and columns.

!           call MatGetSize(dCda, matRows, matCols, PETScIerr)

!           if( PETScIerr/=0 ) &
!             call terminate("createPETScMat", "Error in MatGetSize dCda")

!           write(*,20) "# MATRIX: dCda global size =", &
!                       matRows, " x ", matCols

!           ! Gets the matrix type as a string from the matrix object.

!           call MatGetType(dCda, matTypeStr, PETScIerr)

!           if( PETScIerr/=0 ) &
!             call terminate("createPETScMat", "Error in MatGetType dCda")

!           write(*,30) "# MATRIX: dCda type        =", matTypeStr

!         endif

!         ! Query about the ownership range.

!         if( debug ) then
!           call MatGetOwnershipRange(dCda, iLow, iHigh, PETScIerr)

!           if( PETScIerr/=0 ) &
!             call terminate("createPETScMat", &
!                            "Error in MatGetOwnershipRange dCda")

!           write(*,40) "# MATRIX: dCda Proc", PETScRank, "; #rows =", &
!                nTimeIntervalsSpectral, "; ownership =", iLow, "to", iHigh-1
!        endif


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


subroutine drdwPreAllocation(nnzDiag,nnzOffDiag,wSize)

  ! Get a good estimate of the number of non zero rows for the
  ! on-diagonal and off-diagonal portions of the matrix
  use blockPointers
  use ADjointPETSc
  use ADjointVars    
  use communication   
  use inputTimeSpectral 
  use flowVarRefState 
  use inputADjoint    

  implicit none

  ! Subroutine Arguments
  integer(kind=intType),intent(in)  :: wSize
  integer(kind=intType),intent(out) :: nnzDiag(wSize),nnzOffDiag(Wsize)


  ! Local Variables

  integer(kind=intType) :: nn,i,j,k,counter,sps

  counter = 1
  nnzDiag(:) = 13+(nTimeIntervalsSpectral-1) ! ALWAYS have the center cell ON-PROCESSOR
  nnzOffDiag(:) = 0_intType 
  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointersAdj(nn,1_intType,sps)
        ! Loop over each Cell
        do k=2,kl
           do j=2,jl
              do i=2,il 
                 if (i-2 < 2) nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (i-1 < 2) nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (i+1 > il)nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (i+2 > il)nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (j-2 < 2) nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (j-1 < 2) nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (j+1 > jl)nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (j+2 > jl)nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (k-2 < 2) nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (k-1 < 2) nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (k+1 > kl)nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (k+2 > kl)nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 counter = counter + 1
              end do ! I loop
           end do ! J loop
        end do ! K loop
     end do ! Domain Loop
  end do ! sps Loop
end subroutine drdwPreAllocation
subroutine drdwPCPreAllocation(nnzDiag,nnzOffDiag,wSize)

  ! Get a good estimate of the number of non zero rows for the
  ! on-diagonal and off-diagonal portions of the matrix
  use blockPointers
  use ADjointPETSc
  use ADjointVars    
  use communication   
  use inputTimeSpectral 
  use flowVarRefState 
  use inputADjoint    

  implicit none

  ! Subroutine Arguments
  integer(kind=intType),intent(in)  :: wSize
  integer(kind=intType),intent(out) :: nnzDiag(wSize),nnzOffDiag(Wsize)


  ! Local Variables

  integer(kind=intType) :: nn,i,j,k,counter,sps

  counter = 1
  nnzDiag(:) = 7+(nTimeIntervalsSpectral-1) ! ALWAYS have the center cell ON-PROCESSOR
  nnzOffDiag(:) = 0_intType 
  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointersAdj(nn,1_intType,sps)
        ! Loop over each Cell
        do k=2,kl
           do j=2,jl
              do i=2,il 
                 if (i-1 < 2) nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (i+1 > il)nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (j-1 < 2) nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (j+1 > jl)nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (k-1 < 2) nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 if (k+1 > kl)nnzOffDiag(counter) = nnzOffDiag(counter) + 1
                 counter = counter + 1
              end do ! I loop
           end do ! J loop
        end do ! K loop
     end do ! Domain Loop
  end do ! sps loop
end subroutine drdwPCPreAllocation

subroutine drdxPreAllocation(nnzDiag,nnzOffDiag,wSize)

  ! Get a good estimate of the number of non zero rows for the
  ! on-diagonal and off-diagonal portions of the matrix
  use blockPointers
  use ADjointPETSc
  use ADjointVars    
  use communication   
  use inputTimeSpectral 
  use flowVarRefState 
  use inputADjoint    

  implicit none

  ! Subroutine Arguments
  integer(kind=intType),intent(in)  :: wSize
  integer(kind=intType),intent(out) :: nnzDiag(wSize),nnzOffDiag(Wsize)

  ! Local Variables

  integer(kind=intType) :: nn,i,j,k,l,counter,sps

  counter = 1
  nnzDiag(:) = 32*3 + 8*3*(nTimeIntervalsSpectral-1) ! ALWAYS have the center cell ON-PROCESSOR
  nnzOffDiag(:) = 0_intType 
  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointersAdj(nn,1_intType,sps)

        ! Loop over each Cell
        do k=2,kl
           do j=2,jl
              do i=2,il 
                 do l=1,nw
                    if (i-1 < 2) nnzOffDiag(counter) = nnzOffDiag(counter) + 12
                    if (i+1 > il)nnzOffDiag(counter) = nnzOffDiag(counter) + 12
                    if (j-1 < 2) nnzOffDiag(counter) = nnzOffDiag(counter) + 12
                    if (j+1 > jl)nnzOffDiag(counter) = nnzOffDiag(counter) + 12
                    if (k-1 < 2) nnzOffDiag(counter) = nnzOffDiag(counter) + 12
                    if (k+1 > kl)nnzOffDiag(counter) = nnzOffDiag(counter) + 12
                    counter = counter + 1
                 end do ! l loop
              end do ! I loop
           end do ! J loop
        end do ! K loop
     end do ! Domain Loop
  end do ! sps loop
end subroutine drdxPreAllocation
