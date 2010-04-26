!
!     ******************************************************************
!     *                                                                *
!     * File:          computeADjointGradientExtra.F90                 *
!     * Author:        Andre C. Marta                                  *
!     * Starting date: 07-27-2006                                      *
!     * Last modified: 03-06-2007                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine computeADjointGradientExtra(costFunction)
!
!     ******************************************************************
!     *                                                                *
!     * Computes the total gradient of the cost function J with        *
!     * respect to the extra design variables using the adjoint sol.   *
!     * In addition, it also copies the gradient from the PETSc object *
!     * to the storage array functionGrad in all processors.           *
!     *                                                                *
!     * Compute the total cost/constraint function sensitivity as      *
!     * dIda = dJda - psi^T dRda, where                                *
!     *                                                                *
!     *  dIda [nDesignExtra] is the total gradient of the cost         *
!     *              function with respect to the design variables,    *
!     *                                                                *
!     *  dJda [nDesignExtra] is the partial gradient of the cost       *
!     *              function with respect to the design variables,    *
!     *                                                                *
!     *  psi  [nDim] is the adjoint solution,                          *
!     *                                                                *
!     *  dRda [nDim,nDesignExtra] is the partial gradient of the       *
!     *              residual with respect to the design variables,    *
!     *                                                                *
!     *  with                                                          *
!     *    nDesignExtra = number of extra design variables             *
!     *    nDim         = global number of grid nodes nNodesGlobal     *
!     *                   times the number of equations nw             *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use ADjointVars
      use flowvarrefstate !Timeref
      use inputPhysics
      use constants
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: costFunction
!
!     Local variables.
!
      integer(kind=intType) :: idxmg, iLow, iHigh, n
      integer(kind=intType) :: nDesignLocal, nDisplsLocal
      integer(kind=intType), dimension(:),allocatable :: &
                                           nDesignGlobal, nDisplsGlobal
      real(kind=realType),dimension(:),allocatable :: functionGradLocal

      logical :: designVarPresent

      real(kind=realType), dimension(2) :: time
      real(kind=realType)               :: timeAdjLocal, timeAdj

      character(len=2*maxStringLen) :: errorMessage
      
      !Values for derivative normalization
      real(kind=realType)::a !speed of sound
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
#ifndef USE_NO_PETSC

      ! Send some feedback to screen.

      if( PETScRank==0 ) &
        write(*,10) "Computing total sensitivity wrt Extra..."

      ! Get the initial time.

      call cpu_time(time(1))

      ! Compute the second contribution term psi^T dRda, which requires
      ! a matrix-vector multiplication.

      call MatMultTranspose(dRda,psi,dIda,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("computeADjointGradientExtra", &
                       "Error in MatMultTranspose")

      ! Add the first dJ/da term to get the total gradient and take into
      ! account the negative sign in front of the second contribution.

      call VecAYPX(dIda,PETScNegOne,dJda,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("computeADjointGradientExtra", &
                       "Error in VecAYPX")

      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Computing total sensitivity wrt Extra time (s) =", &
                    timeAdj

      ! View the solution vector dIda.
 
      if(  debug ) then

        if( PETScRank==0 ) then
          write(*,*) "# ============================ "
          write(*,*) "#  dJda partial gradient vect. "
          write(*,*) "# ============================ "
        endif

        call VecView(dJda,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("computeADjointGradientExtra", &
                         "Error in VecView")

        if( PETScRank==0 ) then
          write(*,*) "# ============================ "
          write(*,*) "#  dIda total gradient vector  "
          write(*,*) "# ============================ "
        endif

        call VecView(dIda,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)

        if( PETScIerr/=0 ) &
          call terminate("computeADjointGradientExtra", &
                         "Error in VecView")

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

      call VecGetOwnershipRange(dIda, iLow, iHigh, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("computeADjointGradientExtra", &
                       "Error in VecGetOwnershipRange dJda")

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

        call VecGetValues(dIda, 1, idxmg, &
                          functionGradLocal(n), PETScIerr)

        if( PETScIerr/=0 ) then
          write(errorMessage,99) &
                "Error in VecGetValues for global node", idxmg
          call terminate("computeADjointGradientExtra", errorMessage)
        endif

      enddo

      ! Gather the number of design variables stored per processor
      ! in the root processor.

      allocate( nDesignGlobal(PETScSize) )

      call mpi_gather(nDesignLocal, 1, sumb_integer, &
                      nDesignGlobal, 1, sumb_integer, &
                      0, PETSC_COMM_WORLD, PETScIerr)

      ! Gather the displacement of the number of design variables
      ! per processor in the root processor.

      allocate( nDisplsGlobal(PETScSize) )

      nDisplsLocal = iLow

      call mpi_gather(nDisplsLocal, 1, sumb_integer, &
                      nDisplsGlobal, 1, sumb_integer, &
                      0, PETSC_COMM_WORLD, PETScIerr)

      ! Gather the total gradients in the root processor.
      ! Note: if the local processor does not hold any design variable
      !   then nDesignLocal = 0 and nothing is actually sent to the
      !   processor.

      call mpi_gatherv(functionGradLocal, nDesignLocal, sumb_real,  &
                       functionGrad(costFunction,:), nDesignGlobal, &
                       nDisplsGlobal, sumb_real,                    &
                       0, PETSC_COMM_WORLD, PETScIerr)

      ! Broadcast the total gradients from the root processor to
      ! all processors.

      call mpi_bcast(functionGrad(costFunction,1:nDesignExtra), nDesignExtra, &
                     sumb_real, 0, PETSC_COMM_WORLD, PETScIerr)

      ! Release memory to store the local function gradient values.

      if (allocated(functionGradLocal)) deallocate(functionGradLocal)

      if (allocated(nDesignGlobal))     deallocate(nDesignGlobal)
      if (allocated(nDisplsGlobal))     deallocate(nDisplsGlobal)

      !Divide Timeref out of the solution for the rotational derivatives
      functionGrad(costFunction,nDesignRotX:nDesignRotZ) = functionGrad(costFunction,nDesignRotX:nDesignRotZ)*timeref
!      print *,'timeref',timeref
      !Convert radians to degrees to match design variables
      functionGrad(costFunction,nDesignAOA:nDesignSSA) = functionGrad(costFunction,nDesignAOA:nDesignSSA)*(pi/180.0_realType)

      if( PETScRank==0 ) then	
	print *,'Other Derivatives'
	print *,'Alpha',functionGrad(costFunction,nDesignAoA)
	print *,'Beta',functionGrad(costFunction,nDesignSSA)
	print *,'Mach',functionGrad(costFunction,nDesignMach)
	print *,'MachGrid',functionGrad(costFunction,nDesignMachGrid)	
	print *,'Corrected rotational derivatives...'
	print *,'Rotx:',functionGrad(costFunction,nDesignRotX)
	print *,'RotY:',functionGrad(costFunction,nDesignRotY)
	print *,'RotZ:',functionGrad(costFunction,nDesignRotZ)
      endif

      if( PETScRank==0 ) then	
	print *,'Stability Derivatives per Rad'
	print *,'Alpha',functionGrad(costFunction,nDesignAoA)*(180.0_realType/pi)
	print *,'Beta',functionGrad(costFunction,nDesignSSA)*(180.0_realType/pi)
	print *,'Mach',functionGrad(costFunction,nDesignMach)
	print *,'MachGrid',functionGrad(costFunction,nDesignMachGrid)	
	print *,'Corrected rotational derivatives...'
        a  = sqrt(gammaInf*pInfDim/rhoInfDim)
        print *,'Speed of Sound',a,MachGrid*a,lengthRef
	print *,'Rotx:',functionGrad(costFunction,nDesignRotX)
	print *,'RotY:',functionGrad(costFunction,nDesignRotY)
	print *,'RotZ:',functionGrad(costFunction,nDesignRotZ)*(machGrid*a)/lengthRef
      endif

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)

      ! Output formats.

   10 format(a)
   20 format(a,1x,f8.2)
   99 format(a,1x,i6)

#endif

      end subroutine computeADjointGradientExtra
