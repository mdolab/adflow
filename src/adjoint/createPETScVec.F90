!
!     ******************************************************************
!     *                                                                *
!     * File:          createPETScVec.F90                              *
!     * Author:        C.A.(Sandy) Mader, Andre C. Marta               *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 02-02-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine createPETScVec
!
!     ******************************************************************
!     *                                                                *
!     * Create the global RHS vector of the adjoint system of          *
!     * equations, the solution vector and the residual vector         *
!     * as PETSc vector objects, which are used to compute the adjoint *
!     * solution. Also, create the vectors of partial and total        *
!     * function gradients - dJda,dJdx - and - dIda,dIdx - respectively*
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use ADjointVars     ! nCellsLocal, nDesignExtra
      use communication   ! myID, nProc 
      use inputTimeSpectral !nTimeIntervalsSpectral
      use flowVarRefState ! 
      use mdData          ! mdNSurfNodesCompact
      implicit none
!
!     Local variables.
!
      ! nDim  - local dimension (processor owned)
      ! iLow  - first component owned by the local process
      ! iHigh - one more than last component owned by the local process

      integer       :: nn, iLow, iHigh
      integer       :: nDimW, nDimX,nDimS

      integer       :: vecBlockSize, vecRows
      character(15) :: vecTypeStr

      character(len=2*maxStringLen) :: errorMessage
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
#ifndef USE_NO_PETSC

      ! Define the vector dJdW local size, taking into account the total
      ! number  of nodes owned by the processor and the number of
      ! equations.

      nDimW = nw * nCellsLocal*nTimeIntervalsSpectral

      ! Define matrix dJdx local size for the spatial derivatives.

      nDimX = 3 * nNodesLocal*nTimeIntervalsSpectral

      ! Define vec phic local size (number of Rows) for the
      ! Coupling derivatives.

      call getForceSize(nDimS)
      nDimS = nDimS * 3
!
!     ******************************************************************
!     *                                                                *
!     * Create the right-hand-side vector dJdW that define the linear  *
!     * system, dRdW^T psi = dJdW, from scratch.                       *
!     *                                                                *
!     * Vector dJdW has size [nDim] but can be very sparse depending   *
!     * the cost/constraint function definition J=J(W).                *
!     *                                                                *
!     ******************************************************************
!
      call VecCreate(SUMB_PETSC_COMM_WORLD, dJdW, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", "Error in VecCreate dJdW")
 
      call VecSetSizes(dJdW, nDimW, PETSC_DECIDE, PETScIerr)

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
              "Error in VecSetSizes dJdW for local size", nDimW
        call terminate("createPETScVec", errorMessage)
      endif
      call VecSetFromOptions(dJdW, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", &
                       "Error in VecSetFromOptions dJdW")
     
      call VecSetBlockSize(dJdW, nw, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", &
                       "Error in VecSetBlockSize dJdW")

      ! Extract info from the vector.

      if( PETScRank==0 .and. debug ) then

         call VecGetBlockSize(dJdW, vecBlockSize, PETScIerr)

        if( PETScIerr/=0 ) &
             call terminate("createPETScVec", &
             "Error in VecGetBlockSize dJdW")
        
        write(*,10) "# VECTOR: dJdW block size  =", vecBlockSize

        
        call VecGetSize(dJdW, vecRows, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScVec", "Error in VecGetSize dJdW")

        write(*,20) "# VECTOR: dJdW global size =", vecRows

        call VecGetType(dJdW, vecTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScVec", "Error in VecGetType dJdW")

        write(*,30) "# VECTOR: dJdW type        =", vecTypeStr

      endif
!
!     ******************************************************************
!     *                                                                *
!     * Query about the ownership range.                               *
!     *                                                                *
!     ******************************************************************
!
      if( debug ) then
        call VecGetOwnershipRange(dJdW, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScVec", &
                         "Error in VecGetOwnershipRange dJdW")

        write(*,40) "# VECTOR: dJdW Proc", PETScRank, "; #rows =", &
                    nDimW, "; ownership =", iLow, "to", iHigh-1
      endif
!
!     ******************************************************************
!     *                                                                *
!     * Create the solution vector psi and the residual vector pvr     *
!     * by duplicating the existing vector.                            *
!     *                                                                *
!     ******************************************************************
!
      call VecDuplicate(dJdW, psi, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", &
                       "Error in VecDuplicate dJdW->psi")
      !now zero the psi vector
      
      call VecSet(psi,PETScZero,PETScIerr)
         
         if( PETScIerr/=0 ) &
              call terminate("creatPETScVec", "Error in VecSet:psi")

      !continue with the rest of the vec create operations

      call VecDuplicate(dJdW, pvr, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", &
                       "Error in VecDuplicate dJdW->pvr")
  
      call VecDuplicate(dJdW, dJcdW, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", &
                       "Error in VecDuplicate dJdW->dJcdW")

!
!     ******************************************************************
!     *                                                                *
!     * Create the vectors of partial and total function gradient      *
!     * with respect to the extra design variables dJda and dIda,      *
!     * respectively. Vector dJda has size [nDesignExtra].             *
!     *                                                                *
!     * The global dimension is specified so that the extra design     *
!     * variables are scattered among the processors.                  *
!     * This has to be consistent with the matrix dRda.                *
!     *                                                                *
!     ******************************************************************
!
      call VecCreate(SUMB_PETSC_COMM_WORLD, dJda, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", "Error in VecCreate dJda")

      ! Set the global size and let PETSc decide its local size

      call VecSetSizes(dJda, PETSC_DECIDE, &
                       nDesignExtra, PETScIerr)

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
              "Error in VecSetSizes dJda for global size", nDesignExtra
        call terminate("createPETScVec", errorMessage)
      endif

      ! Set the vector from options.

      call VecSetFromOptions(dJda, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", &
                       "Error in VecSetFromOptions dJda")

      ! Query about the ownership range.

      if( debug ) then
        call VecGetOwnershipRange(dJda, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScVec", &
                         "Error in VecGetOwnershipRange dJda")

        write(*,40) "# VECTOR: dJda Proc", PETScRank, "; #rows =", &
                    nDesignExtra, "; ownership =", iLow, "to", iHigh-1
      endif

      ! Extract info from the vector.

      if( PETScRank==0 .and. debug ) then

        ! Get the vector global size.

        call VecGetSize(dJda, vecRows, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScVec", "Error in VecGetSize dJda")

        write(*,20) "# VECTOR: dJda global size =", vecRows

        ! Gets the vector type as a string from the vector object.

        call VecGetType(dJda, vecTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScVec", "Error in VecGetType dJda")

        write(*,30) "# VECTOR: dJda type        =", vecTypeStr

      endif

      ! Duplicate the vectors as needed.

      call VecDuplicate(dJda, dIda, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", &
                       "Error in VecDuplicate dJda->dIda")
!
!     ******************************************************************
!     *                                                                *
!     * Create the vectors of partial and total function gradient with *
!     * respect to the spatial coordinates dJdx and dIdx, respectively.*
!     *                                                                *
!     * Vector dJdx has size [3*nNodesGlobal] but can be very sparse   *
!     * depending the cost/constraint function definition J=J(x).      *
!     *                                                                *
!     * The local dimensions are specified so that the coordinates x   *
!     * are placed in the local processor.                             *
!     * This has to be consistent with the matrix dRdx.                *
!     *                                                                *
!     ******************************************************************
!
      ! Create the vector. Depending on either this is a sequential or 
      ! parallel run,  PETSc automatically generates the apropriate
      ! vector type over all processes in SUMB_PETSC_COMM_WORLD.

      call VecCreate(SUMB_PETSC_COMM_WORLD, dJdx, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", "Error in VecCreate dJdx")
      
      ! Set the local size and let PETSc determine its global size

      call VecSetSizes(dJdx,nDimX,PETSC_DETERMINE,PETScIerr)

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
              "Error in VecSetSizes dJdx for local size", nDimX
        call terminate("createPETScVec", errorMessage)
      endif

      ! Set the vector from options.

      call VecSetFromOptions(dJdx, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", &
                       "Error in VecSetFromOptions dJdx")

      ! Set the vector block size, for future calls to
      ! VecSetValuesBlocked() and VecSetValuesBlockedLocal().

      call VecSetBlockSize(dJdx, 3, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", &
                       "Error in VecSetBlockSize dJdx")

      ! Query about the ownership range.

      if( debug ) then
        call VecGetOwnershipRange(dJdx, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScVec", &
                         "Error in VecGetOwnershipRange dJdx")

        write(*,40) "# VECTOR: dJdx Proc", PETScRank, "; #rows =", &
                    nDimX, "; ownership =", iLow, "to", iHigh-1
      endif

      ! Extract info from the vector.

      if( PETScRank==0 .and. debug ) then

        ! Get the vector block size.

        call VecGetBlockSize(dJdx, vecBlockSize, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScVec", &
                         "Error in VecGetBlockSize dJdx")

        write(*,10) "# VECTOR: dJdx block size  =", vecBlockSize

        ! Get the vector global size.

        call VecGetSize(dJdx, vecRows, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScVec", "Error in VecGetSize dJdx")

        write(*,20) "# VECTOR: dJdx global size =", vecRows

        ! Gets the vector type as a string from the vector object.

        call VecGetType(dJdx, vecTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScVec", "Error in VecGetType dJdx")

        write(*,30) "# VECTOR: dJdx type        =", vecTypeStr

      endif

      ! Duplicate the vectors as needed.

      call VecDuplicate(dJdx, dIdx, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", &
                       "Error in VecDuplicate dJdx->dIdx")

      call VecDuplicate(dJdx, dJcdx, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", &
                       "Error in VecDuplicate dJdx->dJcdx")

     !
!     ******************************************************************
!     *                                                                *
!     * Create the phic vector for  aero-structural coupling           *
!     *                                                                *
!     ******************************************************************
!
!       call VecCreate(SUMB_PETSC_COMM_WORLD, phic, PETScIerr)
!       call VecSetSizes(phic, nDimS, PETSC_DECIDE, PETScIerr)

!
!     ******************************************************************
!     *                                                                *
!     * Query about the ownership range.                               *
!     *                                                                *
!     ******************************************************************
!
!       if( debug ) then
!         call VecGetOwnershipRange(phic, iLow, iHigh, PETScIerr)

!         if( PETScIerr/=0 ) &
!           call terminate("createPETScVec", &
!                          "Error in VecGetOwnershipRange phic")

!         write(*,40) "# VECTOR: phic Proc", PETScRank, "; #rows =", &
!                     nDimS, "; ownership =", iLow, "to", iHigh-1
!       endif
!
!     ******************************************************************
!     dJdc

!
!     ******************************************************************
!     *                                                                *
!     * Create the right-hand-side vector dJdC that agregates the      *
!     * Time spectral derivatives of the various other components.     *
!     *                                                                *
!     * Vector dJdC has size [nTimeIntervalsSpectral] and is dense.    *
!     *                                                                *
!     ******************************************************************
!
      call VecCreate(SUMB_PETSC_COMM_WORLD, dJdC, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", "Error in VecCreate dJdC")
      
      call VecSetSizes(dJdC, PETSC_DETERMINE, nTimeIntervalsSpectral,&
           PETScIerr)

      if( PETScIerr/=0 ) then
         write(errorMessage,99) &
              "Error in VecSetSizes dJdC for local size", nTimeIntervalsSpectral
        call terminate("createPETScVec", errorMessage)
      endif

      call VecSetFromOptions(dJdC, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", &
                      "Error in VecSetFromOptions dJdC")
      ! Extract info from the vector.

      if( PETScRank==0 .and. debug ) then

        write(*,10) "# VECTOR: dJdC block size  =", vecBlockSize

        call VecGetSize(dJdC, vecRows, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScVec", "Error in VecGetSize dJdC")

        write(*,20) "# VECTOR: dJdc global size =", vecRows

        ! Gets the vector type as a string from the vector object.

        call VecGetType(dJdC, vecTypeStr, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScVec", "Error in VecGetType dJdC")

        write(*,30) "# VECTOR: dJdC type        =", vecTypeStr

      endif
!
!     ******************************************************************
!     *                                                                *
!     * Query about the ownership range.                               *
!     *                                                                *
!     ******************************************************************
!

      if( debug ) then
        call VecGetOwnershipRange(dJdC, iLow, iHigh, PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("createPETScVec", &
                         "Error in VecGetOwnershipRange dJdC")

        write(*,40) "# VECTOR: dJdC Proc", PETScRank, "; #rows =", &
             nTimeIntervalsSpectral, "; ownership =", iLow, "to", iHigh-1
      endif

! ************************************************
! dIdxsDV
! ************************************************

      call VecCreate(SUMB_PETSC_COMM_WORLD, dIdxsDV, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("computeADjointGradientSurface", "Error in VecCreate dIdxs")
      
      ! Set the local size and let PETSc determine its global size
      call VecSetSizes(dIdxsDV,PETSC_DECIDE,3*mdNSurfNodesCompact,PETScIerr)
!

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
              "Error in VecSetSizes dIdxs for global size", mdNSurfNodesCompact
        call terminate("createPETScVec", errorMessage)
      endif

      ! Set the vector from options.

      call VecSetFromOptions(dIdxsDV, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", &
                       "Error in VecSetFromOptions dIdxsDV")

!     ******************************************************************
!
      ! Output formats.

   10 format(a,1x,i2)                            ! block size
   20 format(a,1x,i6)                            ! global size
   30 format(a,1x,a)                             ! type
   40 format(a,1x,i3,a,1x,i6,a,1x,i6,1x,a,1x,i6) ! ownership
   99 format(a,1x,i6)                            ! error

#endif

      end subroutine createPETScVec
