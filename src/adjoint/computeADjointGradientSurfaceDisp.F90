!
!     ******************************************************************
!     *                                                                *
!     * File:          computeADjointGradientSurfaceDisp.F90           *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 07-15-2009                                      *
!     * Last modified: 01-25-2010                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine computeADjointGradientSurfaceDisp(costFunction)
!
!     ******************************************************************
!     *                                                                *
!     * Computes the total gradient of the cost function J with        *
!     * respect to the surface coordinates using the adjoint sol. *
!     * In addition, it also copies the gradient from the PETSc object *
!     * to the storage array functionGradSurface in the root processor.*
!     *                                                                *
!     * Compute the total cost/constraint function sensitivity as      *
!     * dIdx = (dJdxv - psi^T dRdxv)*dxvdxs, where                     *
!     *                                                                *
!     *  dIdx [3*nNodesGlobal*nTimeIntervalsSpectral] is the total     *
!     *            gradient of the cost                                *
!     *            function with respect to the design variables,      *
!     *                                                                *
!     *  dJdx [3*nNodesglobal*nTimeIntervalsSpectral] is the partial   *
!     *            gradient of the cost                                *
!     *            function with respect to the design variables,      *
!     *                                                                *
!     *  psi [nDim] is the adjoint solution,                           *
!     *                                                                *
!     *  dRdx [nDim,3*nNodesGlobal*nTimeIntervalsSpectral] is the      *
!     *             partial gradient of the                            *
!     *             residual with respect to the design variables,     *
!     *                                                                *
!     *  dxvdxs [3*nNodesGlobal*nTimeIntervalsSpectral,                *
!     *          3*ndNSurfNodesCompact*nTimeIntervalsSpectral] is the  *
!     *             derivative of the volume mesh wrt the surface      *
!     *             coordinates(meshwarping)                           *
!     *                                                                *
!     *  with                                                          *
!     *    nNodesGlobal = number of grid nodes                         *
!     *    nDim    = global number of grid nodes X number of equations *
!     *                                                                *
!     ******************************************************************
!
#ifndef USE_NO_PETSC
      use ADjointPETSc
      use ADjointVars
      !use warpingPETSc
      use WarpingPETSc, only: dXvdXsDisp,drdxsDisp,dIdxsDisp,didxs2,djdxs2
      use mdData, only: mdNSurfNodesCompact
      use blockpointers !globalnode
      use inputTimeSpectral !nTimeIntervalsSpectral
      implicit none

!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: costFunction
!
!     Local variables.
!
      real(kind=realType), dimension(2) :: time
      real(kind=realType)               :: timeAdjLocal, timeAdj

!      integer(kind=intType) :: idx!tmp for fd

      character(len=2*maxStringLen) :: errorMessage

      integer(kind=intType) :: idxmg, iLow, iHigh, n 		
      integer(kind=intType) :: nDesignLocal, nDisplsLocal 		
      integer(kind=intType), dimension(:),allocatable :: nDesignGlobal
      integer(kind=intType), dimension(:),allocatable :: nDisplsGlobal 
      real(kind=realType),dimension(:),allocatable :: functionGradLocal
      real(kind=realType),dimension(:,:),allocatable :: functionGradSurfaceDisp2
      logical :: designVarPresent 		
      integer(kind=intType) :: idx
      integer(kind=intType) :: idxlocal,i
      integer(kind=intType),dimension(PETSCSize) :: idxglobal
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!

!
      ! Send some feedback to screen.

      if( PETScRank==0 ) &
        write(*,10) "Computing total sensitivity wrt Surface..."

      !Allocate storage for completed array
      if(.not. allocated(functionGradSurfaceDisp))&
           allocate(functionGradSurfaceDisp(nCostFunction,3*mdNSurfNodesCompact*nTimeIntervalsSpectral))
      if(.not. allocated(functionGradSurfaceDisp2))&
           allocate(functionGradSurfaceDisp2(nCostFunction,3*mdNSurfNodesCompact*nTimeIntervalsSpectral))
      !create the PETSc Vector aswell
      ! Create the vector. Depending on either this is a sequential or 
      ! parallel run,  PETSc automatically generates the apropriate
      ! vector type over all processes in PETSC_COMM_WORLD.

      call VecCreate(PETSC_COMM_WORLD, dIdxsDisp, PETScIerr)
      call VecCreate(PETSC_COMM_WORLD, dIdxs2, PETScIerr)
      call VecCreate(PETSC_COMM_WORLD, dJdxs2, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("computeADjointGradientSurface", "Error in VecCreate dIdxs")
      
      ! Set the local size and let PETSc determine its global size

      call VecSetSizes(dIdxsDisp,PETSC_DECIDE,3*mdNSurfNodesCompact*nTimeIntervalsSpectral,PETScIerr)
      call VecSetSizes(dIdxs2,PETSC_DECIDE,3*mdNSurfNodesCompact*nTimeIntervalsSpectral,PETScIerr)
      call VecSetSizes(dJdxs2,PETSC_DECIDE,3*mdNSurfNodesCompact*nTimeIntervalsSpectral,PETScIerr)

      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
              "Error in VecSetSizes dIdxs for global size", mdNSurfNodesCompact
        call terminate("createPETScVec", errorMessage)
      endif

      ! Set the vector from options.

      call VecSetFromOptions(dIdxsDisp, PETScIerr)
      call VecSetFromOptions(dIdxs2, PETScIerr)
      call VecSetFromOptions(dJdxs2, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("createPETScVec", &
                       "Error in VecSetFromOptions dIdxs")

      ! Get the initial time.

      call cpu_time(time(1))

      ! Compute the second contribution term psi^T dRdx, which requires
      ! a matrix-vector multiplication.

      call MatMultTranspose(dRdx,psi,dIdx,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("computeADjointGradientSurface", &
                       "Error in MatMultTranspose X")

      ! Add the first dJ/dx term to get the total gradient and take into
      ! account the negative sign in front of the second contribution.

      call VecAYPX(dIdx,PETScNegOne,dJdx,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("computeADjointGradientSurface", &
                       "Error in VecAYPX X")

      !now multiply by the volume surface derivative
      call MatMultTranspose(dXvdXsDisp,dIdx,dIdxsDisp,PETScIerr)

!      !now multiply by the volume surface derivative
!      call MatMultTranspose(dXvdXs,dJdx,dIdxs,PETScIerr)

      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Computing total sensitivity wrt Surface time (s) =", &
                    timeAdj

      !Alternate solution      
      !multiply and store in drdxs
      call MatMatMult(dRdx,dXvdXsDisp,MAT_INITIAL_MATRIX,PETSC_DEFAULT_DOUBLE_PRECISION,dRdXsDisp, PETScIerr) 
      
      call MatMultTranspose(dXvdXsDisp,dJdx,dJdxs2,PETScIerr)
      call MatMultTranspose(dRdxsDisp,psi,dIdxs2,PETScIerr)
      call VecAYPX(dIdxs2,PETScNegOne,dJdxs2,PETScIerr)
      ! View the solution vector dIdx.
 
     if(  debug ) then

        if( PETScRank==0 ) then
          write(*,*) "# ============================ "
          write(*,*) "#  dIdxs total gradient vector  "
          write(*,*) "# ============================ "
        endif

        call VecView(dIdxsDisp,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
!        call VecView(dIdx,PETSC_VIEWER_DRAW_WORLD,PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("computeADjointGradientSurface", &
                         "Error in VecView")
	pause

      endif
!
!     ******************************************************************
!     *                                                                *
!     * Transfer solution from PETSc context.                          *
!     *                                                                *
!     ******************************************************************
!

      ! Query about the ownership range.
      ! iHigh is one more than the last element stored locally.

      call VecGetOwnershipRange(dIdxsDisp, iLow, iHigh, PETScIerr)
	
!      print *,'irange',ilow,ihigh,petscrank
!      print *,'irange',iLow,iHigh

      if( PETScIerr/=0 ) &
        call terminate("computeADjointGradientSurface", &
                       "Error in VecGetOwnershipRange dIdxs")
!      RETURN!!!!!!!!	


      ! Determine the number of variables stores in the local processor.
      ! Allocate memory to store the local function gradient values.
      ! Note: it might result in a zero-size array but it's ok!

      nDesignLocal = dim(iHigh,iLow)
	
      allocate( functionGradLocal(nDesignLocal) )

      ! VecGetValues - Gets values from certain locations of a vector.
      !           Currently can only get values on the same processor.

      n = 0
      designVarPresent = .false.

      do idxmg=iLow, iHigh-1

        n = n + 1
        designVarPresent = .true.

        call VecGetValues(dIdxsDisp, 1, idxmg, &
                          functionGradLocal(n), PETScIerr)
	
	!functionGradLocal(n) = idxmg !to debug interface with meshwarping derivatives	

        if( PETScIerr/=0 ) then
           write(errorMessage,99) &
                "Error in VecGetValues for global node", idxmg
           call terminate("computeADjointGradientSurfaceDisp", errorMessage)
        endif

     enddo

      ! Gather the number of design variables stored per processor
      ! in the root processor.

      allocate( nDesignGlobal(PETScSize) )

      !call mpi_gather(nDesignLocal, 1, sumb_integer, &
      !                nDesignGlobal, 1, sumb_integer, &
      !                0, PETSC_COMM_WORLD, PETScIerr)
      call mpi_allgather(nDesignLocal, 1, sumb_integer, &
                      nDesignGlobal, 1, sumb_integer, &
                       PETSC_COMM_WORLD, PETScIerr)

      ! Gather the displacement of the number of design variables
      ! per processor in the root processor.

      allocate( nDisplsGlobal(PETScSize) )

      nDisplsLocal = iLow

      !call mpi_gather(nDisplsLocal, 1, sumb_integer, &
      !                nDisplsGlobal, 1, sumb_integer, &
       !               0, PETSC_COMM_WORLD, PETScIerr)
      call mpi_allgather(nDisplsLocal, 1, sumb_integer, &
                      nDisplsGlobal, 1, sumb_integer, &
                       PETSC_COMM_WORLD, PETScIerr)

      ! Gather the total gradients in the root processor.
      ! Note: if the local processor does not hold any design variable
      !   then nDesignLocal = 0 and nothing is actually sent to the
      !   processor.

!      print *,'costfunction',costFunction,nDisplsGlobal,'shapes',shape(functionGradSpatial),'s2',shape(functionGradLocal)
!	call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)	
!	stop

!      call mpi_gatherv(functionGradLocal, nDesignLocal, sumb_real, &
!                       functionGradSpatial(costFunction,:), nDesignGlobal,&
!                       nDisplsGlobal, sumb_real, &
!                       0, PETSC_COMM_WORLD, PETScIerr)

!	if (PETScRank==0) then
!  		print *,'gathering solution',functionGradLocal, nDesignLocal, sumb_real, &
!                       functionGradSpatial(costFunction,:), nDesignGlobal,&
!                       nDisplsGlobal, sumb_real, &
!                        PETSC_COMM_WORLD, PETScIerr
!	endif


       call mpi_allgatherv(functionGradLocal, nDesignLocal, sumb_real, &
                       functionGradSurfaceDisp(costFunction,:), nDesignGlobal,&
                       nDisplsGlobal, sumb_real, &
                        PETSC_COMM_WORLD, PETScIerr)

      ! Release memory to store the local function gradient values.

      if (allocated(functionGradLocal)) deallocate(functionGradLocal)
      
      if (allocated(nDesignGlobal)) deallocate(nDesignGlobal)
      if (allocated(nDisplsGlobal)) deallocate(nDisplsGlobal)
!!******************
! Get alternate solution
!************
!
!     ******************************************************************
!     *                                                                *
!     * Transfer solution from PETSc context.                          *
!     *                                                                *
!     ******************************************************************
!

      ! Query about the ownership range.
      ! iHigh is one more than the last element stored locally.

      call VecGetOwnershipRange(dIdxs2, iLow, iHigh, PETScIerr)
	
!      print *,'irange',ilow,ihigh,petscrank
!      print *,'irange',iLow,iHigh

      if( PETScIerr/=0 ) &
        call terminate("computeADjointGradientSurface", &
                       "Error in VecGetOwnershipRange dIdxsDisp")
!      RETURN!!!!!!!!	


      ! Determine the number of variables stores in the local processor.
      ! Allocate memory to store the local function gradient values.
      ! Note: it might result in a zero-size array but it's ok!

      nDesignLocal = dim(iHigh,iLow)
	
      allocate( functionGradLocal(nDesignLocal) )

      ! VecGetValues - Gets values from certain locations of a vector.
      !           Currently can only get values on the same processor.

      n = 0
      designVarPresent = .false.

      do idxmg=iLow, iHigh-1

        n = n + 1
        designVarPresent = .true.

        call VecGetValues(dIdxs2, 1, idxmg, &
                          functionGradLocal(n), PETScIerr)
	
	!functionGradLocal(n) = idxmg !to debug interface with meshwarping derivatives	

        if( PETScIerr/=0 ) then
           write(errorMessage,99) &
                "Error in VecGetValues for global node", idxmg
           call terminate("computeADjointGradientSurfaceDisp", errorMessage)
        endif

     enddo

      ! Gather the number of design variables stored per processor
      ! in the root processor.

      allocate( nDesignGlobal(PETScSize) )

      !call mpi_gather(nDesignLocal, 1, sumb_integer, &
      !                nDesignGlobal, 1, sumb_integer, &
      !                0, PETSC_COMM_WORLD, PETScIerr)
      call mpi_allgather(nDesignLocal, 1, sumb_integer, &
                      nDesignGlobal, 1, sumb_integer, &
                       PETSC_COMM_WORLD, PETScIerr)

      ! Gather the displacement of the number of design variables
      ! per processor in the root processor.

      allocate( nDisplsGlobal(PETScSize) )

      nDisplsLocal = iLow

      !call mpi_gather(nDisplsLocal, 1, sumb_integer, &
      !                nDisplsGlobal, 1, sumb_integer, &
       !               0, PETSC_COMM_WORLD, PETScIerr)
      call mpi_allgather(nDisplsLocal, 1, sumb_integer, &
                      nDisplsGlobal, 1, sumb_integer, &
                       PETSC_COMM_WORLD, PETScIerr)

      ! Gather the total gradients in the root processor.
      ! Note: if the local processor does not hold any design variable
      !   then nDesignLocal = 0 and nothing is actually sent to the
      !   processor.


       call mpi_allgatherv(functionGradLocal, nDesignLocal, sumb_real, &
                       functionGradSurfaceDisp2(costFunction,:), nDesignGlobal,&
                       nDisplsGlobal, sumb_real, &
                        PETSC_COMM_WORLD, PETScIerr)
!*********
!Print Both solutions
!**************


!!$      if(PetscRank == 0)then
!!$        
!!$	 print *,'printing result'
!!$         do i = 1,3*mdNSurfNodesCompact
!!$            print *,'costfunction derivative',functionGradSurface(costFunction,i),functionGradSurface2(costfunction,i),i
!!$         enddo
!!$      endif
      

      ! Release memory to store the local function gradient values.

      if (allocated(functionGradLocal)) deallocate(functionGradLocal)
      
      if (allocated(nDesignGlobal)) deallocate(nDesignGlobal)
      if (allocated(nDisplsGlobal)) deallocate(nDisplsGlobal)
      if (allocated(functionGradSurfaceDisp2)) deallocate(functionGradSurfaceDisp2)
      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)

      ! Output format.

   10 format(a)
   20 format(a,1x,f8.2)
   99 format(a,1x,i6)
!
#endif

    end subroutine computeADjointGradientSurfaceDisp
