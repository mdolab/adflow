!
!     ******************************************************************
!     *                                                                *
!     * File:          computeAeroCoupling.F90                         *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 08-14-2008                                      *
!     * Last modified: 08-14-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine computeAeroCoupling(costFunction)
!
!     ******************************************************************
!     *                                                                *
!     * Computes the coupling derivative of the objective function J   *
!     * with respect to the volume mesh coordinates.                   *
!     * This result is passed back to python where it is combined with *
!     * the other half of the coupling derivative.                     *
!     * Compute the coupling cost/constraint function sensitivity as   *
!     * dIdx = psi^T dRdx, where                                       *
!     *                                                                *
!     *  dIdx [3*nNodesGlobal] is the total gradient of the cost       *
!     *            function with respect to the design variables,      *
!     *                                                                *
!     *                                                                *
!     *  psi [nDim] is the adjoint solution,                           *
!     *                                                                *
!     *  dRdx [nDim,3*nNodesGlobal] is the partial gradient of the     *
!     *             residual with respect to the design variables,     *
!     *                                                                *
!     *  with                                                          *
!     *    nNodesGlobal = number of grid nodes                         *
!     *    nDim    = global number of grid nodes X number of equations *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use ADjointVars

      use blockpointers !globalnode

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
      integer(kind=intType), dimension(:),allocatable :: & 		
           nDesignGlobal, nDisplsGlobal 		
      real(kind=realType),dimension(:),allocatable :: functionGradLocal
		
      logical :: designVarPresent 		
#ifndef USE_NO_PETSC
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


      ! Send some feedback to screen.

      if( PETScRank==0 ) &
        write(*,10) "Computing Aero Coupling sensitivity..."

      ! Get the initial time.

      call cpu_time(time(1))

      ! Compute the second contribution term psi^T dRdx, which requires
      ! a matrix-vector multiplication.

      call MatMultTranspose(dRdx,psi,dIdx,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("computeAeroCoupling", &
                       "Error in MatMultTranspose X")
     ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Computing Aero Coupling sensitivity time (s) =", &
                    timeAdj

      ! View the solution vector dIdx.
 
     if(  debug ) then

        if( PETScRank==0 ) then
          write(*,*) "# ============================ "
          write(*,*) "#  dIdx Coupling gradient vector  "
          write(*,*) "# ============================ "
        endif

        call VecView(dIdx,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
!        call VecView(dIdx,PETSC_VIEWER_DRAW_WORLD,PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("computeAeroCoupling", &
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

      call VecGetOwnershipRange(dIdx, iLow, iHigh, PETScIerr)
	
      if( PETScIerr/=0 ) &
        call terminate("computeAeroCoupling", &
                       "Error in VecGetOwnershipRange dIdx")

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

        call VecGetValues(dIdx, 1, idxmg, &
                          functionGradLocal(n), PETScIerr)

        if( PETScIerr/=0 ) then
          write(errorMessage,99) &
                "Error in VecGetValues for global node", idxmg
          call terminate("computeAeroCoupling", errorMessage)
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

!      call mpi_gatherv(functionGradLocal, nDesignLocal, sumb_real, &
!                       functionGradSpatial(costFunction,:), nDesignGlobal,&
!                       nDisplsGlobal, sumb_real, &
!                       0, PETSC_COMM_WORLD, PETScIerr)


       call mpi_allgatherv(functionGradLocal, nDesignLocal, sumb_real, &
                       functionGradCoupling(costFunction,:), nDesignGlobal,&
                       nDisplsGlobal, sumb_real, &
                        PETSC_COMM_WORLD, PETScIerr)

      ! Release memory to store the local function gradient values.

      if (allocated(functionGradLocal)) deallocate(functionGradLocal)

      if (allocated(nDesignGlobal)) deallocate(nDesignGlobal)
      if (allocated(nDisplsGlobal)) deallocate(nDisplsGlobal)

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)

      ! Output format.

   10 format(a)
   20 format(a,1x,f8.2)
   99 format(a,1x,i6)

#endif

    end subroutine computeAeroCoupling
